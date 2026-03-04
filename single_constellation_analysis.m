%% Single Constellation Detailed Delta-V Analysis
% Performs detailed delta-V analysis for a single constellation using live
% TLEs
% 
%
% Author: Yaman Saran
% Team: OSCAR@VT
% 
%
% RAAN change model: Edelbaum approximation
%   deltaV = 2 * V * sin(i) * sin(deltaRAAN/2)
%   t      = deltaV * m / F
%
% J2 precession during spiral altitude change (integrated form):
%   deltaRAAN = (3/16) * J2 * RE^2 * mu * cos(i) * (m/F) * (af^-4 - a0^-4)
%
% TODO
% Daytime burntime to save battery life

clear; clc; close all;

diary('single_constellation_output.txt');
diary on;

fprintf('=======================================================================\n');
fprintf('SINGLE CONSTELLATION DETAILED DELTA-V ANALYSIS\n');
fprintf('Run date: %s\n', datestr(now));
fprintf('=======================================================================\n\n');

%% Parameters
CONSTELLATION_FILE = 'Selected.csv';
RAAN_TOLERANCE = 0.1;
MU_EARTH = 3.986004418e14;  % m^3/s^2
R_EARTH = 6378.137e3;       % m
TLE_CACHE_HOURS = 3;
J2 = 0.00108263;

% Iom Thruster Specifications
THRUST = 0.320;          % Newton
ISP = 1900;              % seconds
SPACECRAFT_MASS = 1662;  % kg
G0 = 9.80665;            % m/s^2

MIN_SPIRAL_REVOLUTIONS = 5;
MAX_SPIRAL_REVOLUTIONS = 500;
CONTINUOUS_BURN_EFFICIENCY = 0.95;
DV_MARGIN_FACTOR = 1.0;

% =========================================================================
% THRUST DUTY CYCLE
%
% Fraction of the orbit the thruster is allowed to fire (0 to 1).
% Models eclipse periods or other downtime.  0.5 = thruster fires
% roughly half the time (daytime only in a dawn-dusk orbit).
% Set to 1.0 for continuous burn (no duty-cycle penalty).
%
% Affects wall-clock duration of all maneuvers and the RAAN
% accumulated during J2 spirals (precession keeps going while
% the thruster is off).
% =========================================================================
THRUST_DUTY_CYCLE = 0.50;  % <-- EDIT THIS (0-1, 1 = always on)

% =========================================================================
% J2 PRECESSION ALTITUDE CHANGE OPTIONS (km) - EDIT THIS VECTOR
% Positive = raise orbit, Negative = lower orbit
% =========================================================================
J2_ALTITUDE_CHANGES_KM = [-300:10:-20,20:10:300];

% =========================================================================
% J2 SELECTED ALTITUDES PER PLANE CHANGE (km)
%
% Specify exactly ONE altitude change per plane change in visit order.
% Length must equal (number of orbital planes - 1).
%
% Example for 4 planes (3 plane changes):
%   J2_SELECTED_ALTITUDES_KM_PER_PC = [150, -75, 200];
%
% Positive = raise orbit, Negative = lower orbit
% Set to [] to skip the selected-altitude analysis entirely.
% =========================================================================
J2_SELECTED_ALTITUDES_KM_PER_PC = [225,200];  % <-- EDIT THIS

% =========================================================================
% INSERTION COST (m/s)
%
% Hard-coded delta-V for orbit insertion to reach the first satellite
% in the visit chain. Propellant is calculated automatically from the
% rocket equation using SPACECRAFT_MASS and ISP.
% Set to 0 to ignore.
% =========================================================================
INSERTION_DV = 300;       % m/s  <-- EDIT THIS
INSERTION_TIME_DAYS = 90;      % days  <-- EDIT THIS (time spent on insertion)
% Propellant for insertion is now computed from the rocket equation:
INSERTION_PROP_KG = SPACECRAFT_MASS * (1 - 1/exp(INSERTION_DV / (ISP * G0)));

% =========================================================================
% OPERATION COST PER SATELLITE
%
% Delta-V and time spent performing operations (docking, servicing,
% inspection, etc.) at each client satellite.
%
% Two parallel arrays indexed by VISIT ORDER (1 = first satellite
% visited, 2 = second, etc.).  Length must equal the total number of
% satellites.  Use 0 for satellites that need no active operation.
%
% Set both to [] to disable.
%
% Example for 7 satellites:
%   OPERATION_DV_PER_SAT        = [5, 5, 5, 10, 5, 5, 5];  % m/s
%   OPERATION_TIME_DAYS_PER_SAT = [1, 1, 1, 2,  1, 1, 1];  % days
% =========================================================================
OPERATION_DV_PER_SAT        = [110,110,110,110];  % <-- EDIT THIS (length = nSatellites)
OPERATION_TIME_DAYS_PER_SAT = [5,5,5,5];  % <-- EDIT THIS (length = nSatellites)

% =========================================================================
% CARGO DROP SCHEDULE
%
% Nx2 matrix: each row is [NORAD_ID, cargo_mass_kg].
% At each listed satellite the chaser drops that cargo mass (kg) after
% rendezvousing, reducing the chaser's wet mass for subsequent maneuvers.
% Set to [] to disable.
%
% Example:
%   CARGO_DROP_SCHEDULE = [25544, 50; 28654, 30];  % drop 50 kg at ISS, 30 kg at NORAD 28654
% =========================================================================
CARGO_DROP_SCHEDULE = [];  % <-- EDIT THIS

LOG_TIME_SCALE = true;
if LOG_TIME_SCALE, scaleStr = 'log'; else, scaleStr = 'linear'; end

fprintf('ION PROPULSION SYSTEM:\n');
fprintf('  Thrust:            %.0f mN\n', THRUST * 1000);
fprintf('  Specific Impulse:  %.0f s\n', ISP);
fprintf('  Exhaust Velocity:  %.0f m/s\n', ISP * G0);
fprintf('  Spacecraft Mass:   %.0f kg\n', SPACECRAFT_MASS);
fprintf('  Acceleration:      %.4f mm/s^2 (%.2e m/s^2)\n', THRUST/SPACECRAFT_MASS*1000, THRUST/SPACECRAFT_MASS);
fprintf('  Delta-V Margin:    %.1f%% (NO SAFETY MARGIN)\n', (DV_MARGIN_FACTOR-1)*100);
fprintf('  Thrust Duty Cycle: %.0f%% (thruster ON fraction)\n', THRUST_DUTY_CYCLE*100);
if INSERTION_DV > 0
    fprintf('  Insertion dV:      %.0f m/s  ->  Prop = %.2f kg (rocket eq.)\n', INSERTION_DV, INSERTION_PROP_KG);
end
fprintf('\n');

%% Load Constellation Data
fprintf('Loading constellation: %s\n', CONSTELLATION_FILE);
fprintf('%s\n', repmat('=', 1, 80));
if ~isfile(CONSTELLATION_FILE), error('Constellation file %s not found!', CONSTELLATION_FILE); end
try
    satellites = loadJCATFormat(CONSTELLATION_FILE, TLE_CACHE_HOURS);
catch ME
    error('Error loading constellation: %s', ME.message);
end
nSatellites = length(satellites);
if nSatellites == 0, error('No valid satellites found!'); end
fprintf('Successfully loaded %d satellites\n\n', nSatellites);

% Validate OPERATION arrays
hasOperations = ~isempty(OPERATION_DV_PER_SAT) || ~isempty(OPERATION_TIME_DAYS_PER_SAT);
if hasOperations
    if isempty(OPERATION_DV_PER_SAT)
        OPERATION_DV_PER_SAT = zeros(1, nSatellites);
    end
    if isempty(OPERATION_TIME_DAYS_PER_SAT)
        OPERATION_TIME_DAYS_PER_SAT = zeros(1, nSatellites);
    end
    if length(OPERATION_DV_PER_SAT) ~= nSatellites
        error('OPERATION_DV_PER_SAT has %d elements but there are %d satellites.', length(OPERATION_DV_PER_SAT), nSatellites);
    end
    if length(OPERATION_TIME_DAYS_PER_SAT) ~= nSatellites
        error('OPERATION_TIME_DAYS_PER_SAT has %d elements but there are %d satellites.', length(OPERATION_TIME_DAYS_PER_SAT), nSatellites);
    end
    fprintf('OPERATION COSTS PER SATELLITE (visit order):\n');
    fprintf('  dV:   ['); fprintf('%.0f ', OPERATION_DV_PER_SAT); fprintf('] m/s\n');
    fprintf('  Time: ['); fprintf('%.1f ', OPERATION_TIME_DAYS_PER_SAT); fprintf('] days\n\n');
end

%% Display Satellite Information
fprintf('SATELLITE CATALOG\n');
fprintf('%s\n', repmat('=', 1, 90));
fprintf('%-4s %-35s %10s %10s %8s %8s %8s %8s\n', 'ID','Name','NORAD','Alt(km)','Inc(deg)','RAAN(deg)','AoP(deg)','MA(deg)');
fprintf('%s\n', repmat('-', 1, 90));
for i = 1:nSatellites
    sat = satellites(i); alt = (sat.a / 1e3) - R_EARTH/1e3;
    fprintf('%-4d %-35s %10d %10.1f %8.2f %8.2f %8.2f %8.2f\n', i, sat.name, sat.norad, alt, sat.i, sat.raan, sat.aop, sat.ma);
end
fprintf('\n');

%% Constellation Statistics
avgAlt = mean([satellites.a]) / 1e3 - R_EARTH/1e3;
minAlt = min([satellites.a]) / 1e3 - R_EARTH/1e3;
maxAlt = max([satellites.a]) / 1e3 - R_EARTH/1e3;
avgInc = mean([satellites.i]); stdInc = std([satellites.i]);
fprintf('CONSTELLATION STATISTICS\n');
fprintf('%s\n', repmat('=', 1, 80));
fprintf('Number of Satellites:     %d\n', nSatellites);
fprintf('Average Altitude:         %.1f km\n', avgAlt);
fprintf('Altitude Range:           %.1f - %.1f km\n', minAlt, maxAlt);
fprintf('Average Inclination:      %.2f deg\n', avgInc);
fprintf('Inclination Std Dev:      %.2f deg\n\n', stdInc);

%% Group into Orbital Planes
fprintf('ORBITAL PLANE ANALYSIS\n');
fprintf('%s\n', repmat('=', 1, 80));
planes = groupIntoPlanes(satellites, RAAN_TOLERANCE);
nPlanes = length(planes);
fprintf('Identified %d orbital planes (RAAN tolerance: %.1f deg)\n\n', nPlanes, RAAN_TOLERANCE);
for p = 1:nPlanes
    fprintf('PLANE %d: RAAN = %.2f deg, Inc = %.2f deg (%d satellites)\n', p, planes(p).raan, planes(p).inc, length(planes(p).satellites));
    fprintf('%s\n', repmat('-', 1, 80));
    fprintf('  %-4s %-35s %10s %10s %8s\n', 'ID', 'Name', 'Alt(km)', 'Ecc', 'MA(deg)');
    fprintf('  %s\n', repmat('-', 1, 75));
    for j = 1:length(planes(p).satellites)
        sat = planes(p).satellites(j); alt = (sat.a / 1e3) - R_EARTH/1e3;
        satID = find(arrayfun(@(s) strcmp(s.name, sat.name), satellites), 1);
        fprintf('  %-4d %-35s %10.1f %10.6f %8.2f\n', satID, sat.name, alt, sat.e, sat.ma);
    end
    fprintf('\n');
end

%% Optimize Plane Visit Order
fprintf('OPTIMIZING PLANE VISIT ORDER\n');
fprintf('%s\n', repmat('=', 1, 80));
[planeOrder, totalPlaneChangeCost] = solvePlaneTSP(planes);
fprintf('Optimized plane visit sequence: ');
for p = 1:nPlanes, fprintf('%d', planeOrder(p)); if p < nPlanes, fprintf(' -> '); end; end
fprintf('\nRAANs in visit order:            ');
for p = 1:nPlanes, fprintf('%.1f deg', planes(planeOrder(p)).raan); if p < nPlanes, fprintf(' -> '); end; end
fprintf('\nPlane transition cost estimate:  %.2f (cost units)\n\n', totalPlaneChangeCost);

%% Time Allocation (Continuous Thrust)
fprintf('OPTIMIZED TIME ALLOCATION (CONTINUOUS THRUST SPIRALS)\n');
fprintf('%s\n', repmat('=', 1, 80));
nPlaneChanges = nPlanes - 1;
phasingTimePerPlane = zeros(nPlanes, 1); planeChangeTime = zeros(nPlaneChanges, 1);
propellantPerPlane = zeros(nPlanes, 1); propellantPerPlaneChange = zeros(nPlaneChanges, 1);
currentMass = SPACECRAFT_MASS; v_e = ISP * G0;

for p = 1:nPlanes
    planeIdx = planeOrder(p); plane = planes(planeIdx);
    nSatsInPlane = length(plane.satellites); nPhasingManeuvers = max(nSatsInPlane - 1, 0);
    if nPhasingManeuvers > 0
        avgDeltaVperManeuver = 40; massRatio = exp(avgDeltaVperManeuver / v_e);
        propellantPerManeuver = currentMass * (1 - 1/massRatio);
        burnTimePerManeuver = currentMass * v_e / (THRUST * THRUST_DUTY_CYCLE) * (1 - 1/massRatio);
        phasingTimePerPlane(p) = burnTimePerManeuver * 2.5 * nPhasingManeuvers;
        propellantPerPlane(p) = propellantPerManeuver * nPhasingManeuvers;
        currentMass = currentMass - propellantPerPlane(p);
    end
end

for i = 1:nPlaneChanges
    currentPlane = planes(planeOrder(i)); nextPlane = planes(planeOrder(i+1));
    avgA = (mean([currentPlane.satellites.a]) + mean([nextPlane.satellites.a])) / 2;
    avgV = sqrt(MU_EARTH / avgA);
    i1_rad = deg2rad(mean([currentPlane.satellites.i])); i2_rad = deg2rad(mean([nextPlane.satellites.i]));
    i_avg_rad = (i1_rad + i2_rad) / 2;
    dRAAN = nextPlane.raan - currentPlane.raan;
    if dRAAN > 180, dRAAN = dRAAN - 360; end
    if dRAAN < -180, dRAAN = dRAAN + 360; end
    raanChange_rad = deg2rad(abs(dRAAN)); incChange_rad = abs(i2_rad - i1_rad);
    dv_raan = 2 * avgV * sin(i_avg_rad) * sin(raanChange_rad / 2);
    dv_inc = avgV * incChange_rad;
    deltaV_planeChange = sqrt(dv_raan^2 + dv_inc^2);
    planeChangeTime(i) = deltaV_planeChange * currentMass / (THRUST * THRUST_DUTY_CYCLE);
    massRatio = exp(deltaV_planeChange / v_e);
    propellantPerPlaneChange(i) = currentMass * (1 - 1/massRatio);
    currentMass = currentMass - propellantPerPlaneChange(i);
end

totalPhasingTime = sum(phasingTimePerPlane); totalPlaneChangeTime = sum(planeChangeTime);
totalMissionTime = totalPhasingTime + totalPlaneChangeTime;
totalPropellant = sum(propellantPerPlane) + sum(propellantPerPlaneChange);
totalPhasingManeuvers = nSatellites - nPlanes;
avgTimePerPhasing = 0;
if totalPhasingManeuvers > 0, avgTimePerPhasing = totalPhasingTime / totalPhasingManeuvers; end

fprintf('Continuous Thrust Strategy (Edelbaum RAAN Model):\n');
fprintf('  - Spiral transfers for all maneuvers\n');
fprintf('  - RAAN change: Edelbaum dV = 2V*sin(i)*sin(dRAAN/2)\n');
fprintf('  - RAAN time:   t = dV*m/F\n');
fprintf('  - Delta-V margin: %.1f%%\n', (DV_MARGIN_FACTOR - 1) * 100);
fprintf('  - Thrust duty cycle: %.0f%%\n\n', THRUST_DUTY_CYCLE * 100);
fprintf('ESTIMATED MISSION TIMELINE:\n');
fprintf('Total phasing time:           %.2f days (%.1f%%)\n', totalPhasingTime/(24*3600), 100*totalPhasingTime/totalMissionTime);
fprintf('Total plane change time:      %.2f days (%.1f%%)\n', totalPlaneChangeTime/(24*3600), 100*totalPlaneChangeTime/totalMissionTime);
fprintf('TOTAL MISSION TIME:           %.2f days (%.2f months, %.2f years)\n', totalMissionTime/(24*3600), totalMissionTime/(24*3600*30), totalMissionTime/(24*3600*365));
fprintf('\nESTIMATED PROPELLANT BUDGET:\n');
fprintf('Total propellant consumed:    %.2f kg (%.1f%% of spacecraft mass)\n', totalPropellant, 100*totalPropellant/SPACECRAFT_MASS);
fprintf('Final spacecraft mass:        %.2f kg\n', currentMass);
fprintf('Propellant for phasing:       %.2f kg\n', sum(propellantPerPlane));
fprintf('Propellant for plane changes: %.2f kg\n', sum(propellantPerPlaneChange));
fprintf('Average time per phasing:     %.2f days\n', avgTimePerPhasing/(24*3600));
if nPlaneChanges > 0, fprintf('Average time per plane change: %.2f days\n', mean(planeChangeTime)/(24*3600)); end
fprintf('\n');

%% Detailed Maneuver Sequence
fprintf('DETAILED MANEUVER SEQUENCE (EDELBAUM CONTINUOUS THRUST)\n');
fprintf('%s\n', repmat('=', 1, 110));
totalPhasingDV = 0; totalPlaneChangeDV = 0; maneuverLog = [];
satVisitOrder = []; cumulativeDV = []; cumulativeTime = []; cumulativePropellant = [];
currentTime = 0; currentMass = SPACECRAFT_MASS; maneuverID = 0;

for planeIdx = 1:nPlanes
    currentPlaneIdx = planeOrder(planeIdx); plane = planes(currentPlaneIdx);
    nSatsInPlane = length(plane.satellites);
    fprintf('\n=== PLANE %d (RAAN = %.2f deg, Inc = %.2f deg) ===\n', currentPlaneIdx, plane.raan, plane.inc);
    fprintf('%s\n', repmat('-', 1, 110));
    [phasingDV, visitOrder, phasingDetails] = calculateContinuousThrustPhasingDV(plane.satellites, currentMass, THRUST, ISP, G0, MU_EARTH, MIN_SPIRAL_REVOLUTIONS, MAX_SPIRAL_REVOLUTIONS, THRUST_DUTY_CYCLE);
    fprintf('Visiting %d satellites in this plane...\n', nSatsInPlane);
    fprintf('Allocated time: %.2f days\n\n', phasingTimePerPlane(planeIdx)/(24*3600));
    fprintf('%-5s %-35s %-35s %10s %10s %10s %10s %10s\n', 'Man#','From Satellite','To Satellite','Phase(deg)','dV(m/s)','Prop(kg)','Time(days)','Revs');
    fprintf('%s\n', repmat('-', 1, 140));
    for i = 1:length(phasingDetails)
        maneuverID = maneuverID + 1; detail = phasingDetails(i);
        currentSat = plane.satellites(visitOrder(i)); nextSat = plane.satellites(visitOrder(i+1));
        currentSatID = find(arrayfun(@(s) strcmp(s.name, currentSat.name), satellites), 1);
        nextSatID = find(arrayfun(@(s) strcmp(s.name, nextSat.name), satellites), 1);
        satVisitOrder = [satVisitOrder, currentSatID];
        if i == length(phasingDetails), satVisitOrder = [satVisitOrder, nextSatID]; end
        detail.dv = detail.dv * DV_MARGIN_FACTOR; detail.propellant = detail.propellant * DV_MARGIN_FACTOR;
        fprintf('%-5d %-35s %-35s %10.2f %10.2f %10.3f %10.2f %10.1f\n', maneuverID, truncateString(sprintf('[%d] %s', currentSatID, currentSat.name), 35), truncateString(sprintf('[%d] %s', nextSatID, nextSat.name), 35), detail.phaseAngle, detail.dv, detail.propellant, detail.burnTime/(24*3600), detail.revolutions);
        maneuverLog(maneuverID).type = 'Phasing'; maneuverLog(maneuverID).fromSatID = currentSatID;
        maneuverLog(maneuverID).toSatID = nextSatID; maneuverLog(maneuverID).fromSatName = currentSat.name;
        maneuverLog(maneuverID).toSatName = nextSat.name; maneuverLog(maneuverID).dv = detail.dv;
        maneuverLog(maneuverID).time = detail.burnTime; maneuverLog(maneuverID).propellant = detail.propellant;
        maneuverLog(maneuverID).revolutions = detail.revolutions; maneuverLog(maneuverID).planeID = currentPlaneIdx;
        maneuverLog(maneuverID).phaseAngle = detail.phaseAngle;
        maneuverLog(maneuverID).massBefore = currentMass; maneuverLog(maneuverID).massAfter = currentMass - detail.propellant;
        totalPhasingDV = totalPhasingDV + detail.dv; currentTime = currentTime + detail.burnTime;
        currentMass = currentMass - detail.propellant;
        cumulativeDV = [cumulativeDV, totalPhasingDV + totalPlaneChangeDV];
        cumulativeTime = [cumulativeTime, currentTime];
        cumulativePropellant = [cumulativePropellant, SPACECRAFT_MASS - currentMass];
    end
    fprintf('\nPlane %d Summary: %d sats, %.2f m/s dV, %.2f kg prop, %.2f days\n', currentPlaneIdx, nSatsInPlane, phasingDV, propellantPerPlane(planeIdx), phasingTimePerPlane(planeIdx)/(24*3600));

    if planeIdx < nPlanes
        nextPlaneIdx = planeOrder(planeIdx + 1); nextPlane = planes(nextPlaneIdx);
        currentSat = plane.satellites(end); nextSat = nextPlane.satellites(1);
        currentSatID = find(arrayfun(@(s) strcmp(s.name, currentSat.name), satellites), 1);
        nextSatID = find(arrayfun(@(s) strcmp(s.name, nextSat.name), satellites), 1);
        [planeChangeDV, pcDetails] = calculateEdelbaumPlaneChange(currentSat, nextSat, currentMass, THRUST, ISP, G0, MU_EARTH, J2, R_EARTH, J2_ALTITUDE_CHANGES_KM, THRUST_DUTY_CYCLE);
        planeChangeDV = planeChangeDV * DV_MARGIN_FACTOR;
        pcDetails.propellant = pcDetails.propellant * DV_MARGIN_FACTOR;
        pcDetails.burnTime = pcDetails.burnTime * DV_MARGIN_FACTOR;
        maneuverID = maneuverID + 1;
        fprintf('\n%s\n', repmat('=', 1, 140));
        fprintf('PLANE CHANGE #%d: Plane %d -> Plane %d (EDELBAUM CONTINUOUS THRUST)\n', planeIdx, currentPlaneIdx, nextPlaneIdx);
        fprintf('%s\n', repmat('-', 1, 140));
        fprintf('From: [%d] %s (Plane %d)\n', currentSatID, currentSat.name, currentPlaneIdx);
        fprintf('To:   [%d] %s (Plane %d)\n', nextSatID, nextSat.name, nextPlaneIdx);
        fprintf('RAAN change:    %.2f deg -> %.2f deg (d = %.2f deg)\n', plane.raan, nextPlane.raan, pcDetails.raanChange);
        fprintf('Inc change:     %.2f deg -> %.2f deg (d = %.2f deg)\n', plane.inc, nextPlane.inc, pcDetails.incChange);
        fprintf('Combined angle: %.2f deg\n', pcDetails.combinedAngle);
        fprintf('Altitude diff:  %.2f km\n\n', pcDetails.altitudeChange);
        fprintf('EDELBAUM DELTA-V BREAKDOWN:\n');
        fprintf('  RAAN rotation:      %.2f m/s  [dV = 2V*sin(i)*sin(dRAAN/2)]\n', pcDetails.dv_raan);
        fprintf('  Inclination change: %.2f m/s  [dV = V*di]\n', pcDetails.dv_inc);
        fprintf('  Combined (vector):  %.2f m/s\n', sqrt(pcDetails.dv_raan^2 + pcDetails.dv_inc^2));
        fprintf('  Total:              %.2f m/s\n\n', planeChangeDV);
        fprintf('EDELBAUM TIME: t = dV * m / F\n');
        fprintf('  Burn time:          %.2f days (%.2f hours)\n', pcDetails.burnTime/(24*3600), pcDetails.burnTime/3600);
        fprintf('  Propellant:         %.2f kg\n', pcDetails.propellant);
        fprintf('  Spiral revolutions: %.1f\n', pcDetails.revolutions);
        fprintf('  Mass before:        %.2f kg\n', currentMass);
        fprintf('  Mass after:         %.2f kg\n', currentMass - pcDetails.propellant);
        fprintf('  Strategy:           %s\n', pcDetails.strategy);

        fprintf('\nJ2 DIFFERENTIAL PRECESSION ALTERNATIVES:\n');
        fprintf('(Slow continuous-thrust spiral to change altitude, wait, spiral back)\n');
        fprintf('Nominal precession rate: %.4f deg/day at %.1f km altitude\n', pcDetails.j2Analysis.nominalPrecessionRate, pcDetails.j2Analysis.nominalAltitude);
        if pcDetails.j2Analysis.raanChangeRequired > 0, raanDirStr = 'eastward'; else, raanDirStr = 'westward'; end
        fprintf('RAAN change required:    %+.2f degrees (%s)\n', pcDetails.j2Analysis.raanChangeRequired, raanDirStr);
        fprintf('RAAN spiral integral:    dRAAN = (3/16)*J2*RE^2*mu*cos(i)*(m/F)*(af^-4 - a0^-4)\n\n');
        fprintf('%-12s %-12s %-14s %-14s %-18s %-14s %-14s %-12s %-24s\n', 'Alt d(km)','New Alt(km)','dV Total(m/s)','Spiral(days)','RAAN Spiral(deg)','Wait(days)','Total(days)','Prop(kg)','Direction');
        fprintf('%s\n', repmat('-', 1, 170));
        for opt = 1:length(pcDetails.j2Analysis.options)
            j2opt = pcDetails.j2Analysis.options(opt);
            if j2opt.totalTime < Inf
                fprintf('%-+12.0f %-12.1f %-14.1f %-14.1f %-+18.4f %-14.1f %-14.1f %-12.2f %-24s\n', j2opt.deltaAlt, j2opt.newAlt, j2opt.dv_total, j2opt.spiralTime, j2opt.spiralRAAN, j2opt.waitTime, j2opt.totalTime, j2opt.propellant, j2opt.direction);
            else
                fprintf('%-+12.0f %-12.1f %-14s %-14s %-18s %-14s %-14s %-12s %-24s\n', j2opt.deltaAlt, j2opt.newAlt, 'Infeasible', '-', '-', '-', '-', '-', j2opt.direction);
            end
        end
        fprintf('\nCOMPARISON: Active Edelbaum vs Best J2 Strategy\n');
        feasibleIdx = find([pcDetails.j2Analysis.options.totalTime] < Inf);
        if ~isempty(feasibleIdx)
            feasibleTimes = [pcDetails.j2Analysis.options(feasibleIdx).totalTime];
            [~, minIdx] = min(feasibleTimes); bestJ2 = pcDetails.j2Analysis.options(feasibleIdx(minIdx));
            fprintf('  Active Edelbaum: %.0f m/s, %.1f days, %.1f kg propellant\n', planeChangeDV, pcDetails.burnTime/(24*3600), pcDetails.propellant);
            fprintf('  J2 (%+.0f km):   %.0f m/s, %.1f days, %.1f kg propellant [%s]\n', bestJ2.deltaAlt, bestJ2.dv_total, bestJ2.totalTime, bestJ2.propellant, bestJ2.direction);
            fprintf('  SAVINGS:         %.1f%% dV, %.1f%% propellant\n', 100*(1-bestJ2.dv_total/planeChangeDV), 100*(1-bestJ2.propellant/pcDetails.propellant));
        else
            fprintf('  No feasible J2 strategy found.\n');
        end
        fprintf('%s\n', repmat('=', 1, 140));

        maneuverLog(maneuverID).type = 'Plane Change'; maneuverLog(maneuverID).fromSatID = currentSatID;
        maneuverLog(maneuverID).toSatID = nextSatID; maneuverLog(maneuverID).fromSatName = currentSat.name;
        maneuverLog(maneuverID).toSatName = nextSat.name; maneuverLog(maneuverID).dv = planeChangeDV;
        maneuverLog(maneuverID).time = pcDetails.burnTime; maneuverLog(maneuverID).propellant = pcDetails.propellant;
        maneuverLog(maneuverID).revolutions = pcDetails.revolutions;
        maneuverLog(maneuverID).fromPlaneID = currentPlaneIdx; maneuverLog(maneuverID).toPlaneID = nextPlaneIdx;
        maneuverLog(maneuverID).raanChange = pcDetails.raanChange; maneuverLog(maneuverID).incChange = pcDetails.incChange;
        maneuverLog(maneuverID).combinedAngle = pcDetails.combinedAngle; maneuverLog(maneuverID).strategy = pcDetails.strategy;
        maneuverLog(maneuverID).massBefore = currentMass; maneuverLog(maneuverID).massAfter = currentMass - pcDetails.propellant;
        totalPlaneChangeDV = totalPlaneChangeDV + planeChangeDV; currentTime = currentTime + pcDetails.burnTime;
        currentMass = currentMass - pcDetails.propellant;
        cumulativeDV = [cumulativeDV, totalPhasingDV + totalPlaneChangeDV];
        cumulativeTime = [cumulativeTime, currentTime];
        cumulativePropellant = [cumulativePropellant, SPACECRAFT_MASS - currentMass];
    end
end

totalDV = totalPhasingDV + totalPlaneChangeDV; totalPropellantUsed = SPACECRAFT_MASS - currentMass;
fprintf('\n%s\n', repmat('=', 1, 110));
fprintf('MISSION SUMMARY (EDELBAUM CONTINUOUS THRUST WITH ION PROPULSION)\n');
fprintf('%s\n', repmat('=', 1, 110));
fprintf('Total satellites visited:     %d\n', nSatellites);
fprintf('Total phasing maneuvers:      %d\n', nSatellites - nPlanes);
fprintf('Total plane change maneuvers: %d\n', nPlaneChanges);
fprintf('Total maneuvers:              %d\n\n', length(maneuverLog));
fprintf('DELTA-V BUDGET:\n');
fprintf('Total phasing delta-V:        %.2f m/s (%.3f km/s)\n', totalPhasingDV, totalPhasingDV/1000);
fprintf('Total plane change delta-V:   %.2f m/s (%.3f km/s)\n', totalPlaneChangeDV, totalPlaneChangeDV/1000);
fprintf('TOTAL MISSION DELTA-V:        %.2f m/s (%.3f km/s)\n', totalDV, totalDV/1000);
fprintf('Average delta-V per satellite: %.2f m/s\n\n', totalDV/nSatellites);
fprintf('PROPELLANT BUDGET:\n');
fprintf('Initial spacecraft mass:      %.2f kg\n', SPACECRAFT_MASS);
fprintf('Total propellant consumed:    %.2f kg (%.1f%% of initial mass)\n', totalPropellantUsed, 100*totalPropellantUsed/SPACECRAFT_MASS);
fprintf('Final spacecraft mass:        %.2f kg\n', currentMass);
phasingPropTotal = sum([maneuverLog(strcmp({maneuverLog.type}, 'Phasing')).propellant]);
pcPropTotal = sum([maneuverLog(strcmp({maneuverLog.type}, 'Plane Change')).propellant]);
fprintf('Propellant for phasing:       %.2f kg (%.1f%%)\n', phasingPropTotal, 100*phasingPropTotal/totalPropellantUsed);
fprintf('Propellant for plane changes: %.2f kg (%.1f%%)\n\n', pcPropTotal, 100*pcPropTotal/totalPropellantUsed);
fprintf('MISSION DURATION:\n');
fprintf('Total mission time:           %.2f days (%.2f months, %.2f years)\n', currentTime/(24*3600), currentTime/(24*3600*30), currentTime/(24*3600*365));
fprintf('Total spiral revolutions:     %.1f\n\n', sum([maneuverLog.revolutions]));
fprintf('ION THRUSTER PERFORMANCE:\n');
fprintf('Thrust:                       %.0f mN\n', THRUST*1000);
fprintf('Specific impulse:             %.0f s\n', ISP);
fprintf('Exhaust velocity:             %.0f m/s\n', ISP*G0);
fprintf('Average acceleration:         %.4f mm/s^2\n', THRUST/((SPACECRAFT_MASS+currentMass)/2)*1000);
fprintf('Thrust duty cycle:            %.0f%%\n', THRUST_DUTY_CYCLE*100);
fprintf('%s\n', repmat('=', 1, 110));

%% Export
maneuverTable = struct2table(maneuverLog);
writetable(maneuverTable, 'maneuver_log.csv');
fprintf('\nManeuver log exported to: maneuver_log.csv\n');

%% Visualizations (refactored into function calls)
plotMissionTimeline(cumulativeTime, cumulativeDV, cumulativePropellant, maneuverLog, totalPhasingDV, totalPlaneChangeDV, nPlanes, CONSTELLATION_FILE);
plotConstellation3D(planes, planeOrder, nPlanes, R_EARTH, MU_EARTH);
plotVisitSequence(satVisitOrder, satellites);

%% J2 Trade Study Plots (per plane change)
planeChangeManeuvers = find(strcmp({maneuverLog.type}, 'Plane Change'));
for pcIdx = 1:length(planeChangeManeuvers)
    mIdx = planeChangeManeuvers(pcIdx);
    fromPlaneIdx = maneuverLog(mIdx).fromPlaneID; toPlaneIdx = maneuverLog(mIdx).toPlaneID;
    fromSat = planes(fromPlaneIdx).satellites(end); toSat = planes(toPlaneIdx).satellites(1);
    [~, pcDet] = calculateEdelbaumPlaneChange(fromSat, toSat, maneuverLog(mIdx).massBefore, THRUST, ISP, G0, MU_EARTH, J2, R_EARTH, J2_ALTITUDE_CHANGES_KM, THRUST_DUTY_CYCLE);
    plotJ2TradeStudySingle(pcDet, maneuverLog(mIdx), fromPlaneIdx, toPlaneIdx, LOG_TIME_SCALE);
end

%% Aggregated J2 Trade Study
if ~isempty(planeChangeManeuvers)
    plotJ2TradeStudyAggregate(planeChangeManeuvers, maneuverLog, planes, THRUST, ISP, G0, MU_EARTH, J2, R_EARTH, J2_ALTITUDE_CHANGES_KM, LOG_TIME_SCALE, THRUST_DUTY_CYCLE);
end

%% =====================================================================
%  SELECTED J2 ALTITUDE MISSION ANALYSIS
%  Evaluates the full mission using user-chosen J2 altitudes per plane
%  change.  Tracks per-satellite costs: RAAN, inclination, phasing,
%  insertion, and operations.
%  Phasing after RAAN changes uses mean-anomaly propagation so the
%  chaser's position at the end of the plane change is accounted for.
% =====================================================================

if ~isempty(J2_SELECTED_ALTITUDES_KM_PER_PC)
    if length(J2_SELECTED_ALTITUDES_KM_PER_PC) ~= nPlaneChanges
        warning('J2_SELECTED_ALTITUDES_KM_PER_PC has %d elements but there are %d plane changes. Skipping.', length(J2_SELECTED_ALTITUDES_KM_PER_PC), nPlaneChanges);
    else
        fprintf('\n%s\n', repmat('=', 1, 130));
        fprintf('SELECTED J2 ALTITUDE MISSION ANALYSIS  (per-satellite breakdown)\n');
        fprintf('Using user-specified J2 altitude for each plane change\n');
        fprintf('%s\n', repmat('=', 1, 130));
        fprintf('Selected altitudes: [');
        for k = 1:nPlaneChanges
            fprintf('%+.0f', J2_SELECTED_ALTITUDES_KM_PER_PC(k));
            if k < nPlaneChanges, fprintf(', '); end
        end
        fprintf('] km\n\n');

        selMass = SPACECRAFT_MASS; selTime = 0; selTotalDV = 0;
        v_e_sel = ISP * G0;

        % Per-satellite tracking
        perSat = struct('name',{},'satID',{},'planeID',{},...
                        'dv_raan',{},'dv_inc',{},'dv_phasing',{},...
                        'dv_insertion',{},'dv_operation',{},...
                        'time_raan',{},'time_inc',{},'time_phasing',{},...
                        'time_insertion',{},'time_operation',{},...
                        'prop_raan',{},'prop_inc',{},'prop_phasing',{},...
                        'prop_insertion',{},'prop_operation',{},...
                        'cargo_dropped',{});
        satCounter = 0;
        selCumDV = []; selCumTime = []; selCumProp = [];
        satArrivalDV = [];    % cumulative dV at each satellite arrival (after phasing, before ops)
        satArrivalTime = [];  % mission time (s) at each satellite arrival
        satArrivalMass = [];  % chaser mass (kg) at each satellite arrival
        massTimeline = [0, SPACECRAFT_MASS];  % [time_s, mass_kg] for mass-vs-time plot
        pcDetailLog = struct('pcNumber',{},'fromPlane',{},'toPlane',{},...
            'chosenAlt',{},'dv',{},'time_days',{},'propellant',{},...
            'spiralTime',{},'waitTime',{},'spiralRAAN',{},...
            'direction',{},'strategy',{});
        pcCounter = 0;
        chaserMA = NaN;  % will be set once we start visiting

        for planeVisit = 1:nPlanes
            pIdx = planeOrder(planeVisit); plane = planes(pIdx);
            nSatsInPlane = length(plane.satellites);
            a_plane = mean([plane.satellites.a]);

            % Sort satellites in this plane by MA (ascending)
            [~, sortOrder] = sort([plane.satellites.ma]);

            % Determine chaser MA at entry to this plane
            if planeVisit == 1
                % First plane: start at first satellite -> no phasing cost
                chaserMA = plane.satellites(sortOrder(1)).ma;
            end

            for si = 1:nSatsInPlane
                sIdx = sortOrder(si);
                sat = plane.satellites(sIdx);
                satCounter = satCounter + 1;
                satGlobalID = find(arrayfun(@(s) strcmp(s.name, sat.name), satellites), 1);
                perSat(satCounter).name   = sat.name;
                perSat(satCounter).satID  = satGlobalID;
                perSat(satCounter).planeID = pIdx;

                % ---------------------------------------------------------
                % Insertion cost: first satellite only
                % Propellant computed from rocket equation
                % ---------------------------------------------------------
                if satCounter == 1 && (INSERTION_DV > 0 || INSERTION_TIME_DAYS > 0)
                    ins_prop = selMass * (1 - 1/exp(INSERTION_DV / v_e_sel));
                    perSat(satCounter).dv_insertion   = INSERTION_DV;
                    perSat(satCounter).prop_insertion  = ins_prop;
                    perSat(satCounter).time_insertion  = INSERTION_TIME_DAYS;
                    selMass = selMass - ins_prop;
                    t_ins = INSERTION_TIME_DAYS * 86400;
                    selTime = selTime + t_ins;
                    selTotalDV = selTotalDV + INSERTION_DV;
                    if INSERTION_DV > 0 || INSERTION_TIME_DAYS > 0
                        selCumDV = [selCumDV, selTotalDV];
                        selCumTime = [selCumTime, selTime];
                        selCumProp = [selCumProp, SPACECRAFT_MASS - selMass];
                        massTimeline = [massTimeline; selTime, selMass];
                    end
                else
                    perSat(satCounter).dv_insertion   = 0;
                    perSat(satCounter).prop_insertion  = 0;
                    perSat(satCounter).time_insertion  = 0;
                end

                % ---------------------------------------------------------
                %    RAAN cost: attribute to first satellite after a
                %    plane change
                % ---------------------------------------------------------
                if si == 1 && planeVisit > 1
                    pc = pcDetailLog(pcCounter);
                    perSat(satCounter).dv_raan   = pc.dv;
                    perSat(satCounter).time_raan  = pc.time_days;
                    perSat(satCounter).prop_raan  = pc.propellant;
                else
                    perSat(satCounter).dv_raan   = 0;
                    perSat(satCounter).time_raan  = 0;
                    perSat(satCounter).prop_raan  = 0;
                end

                % Inclination correction (plane avg -> this sat)
                [dv_inc, tInc, pInc] = calculateSmallIncChange(...
                    plane.inc, sat.i, sat.a, MU_EARTH, selMass, THRUST, ISP, G0, THRUST_DUTY_CYCLE);
                dv_inc = dv_inc * DV_MARGIN_FACTOR;
                pInc = pInc * DV_MARGIN_FACTOR;
                perSat(satCounter).dv_inc   = dv_inc;
                perSat(satCounter).time_inc  = tInc / 86400;
                perSat(satCounter).prop_inc  = pInc;
                selMass = selMass - pInc;
                selTime = selTime + tInc;
                selTotalDV = selTotalDV + dv_inc;
                selCumDV = [selCumDV, selTotalDV];
                selCumTime = [selCumTime, selTime];
                selCumProp = [selCumProp, SPACECRAFT_MASS - selMass];
                massTimeline = [massTimeline; selTime, selMass];

                % Phasing to this satellite
                if si == 1 && planeVisit == 1
                    dv_ph = 0; tPh = 0; pPh = 0;
                else
                    targetMA = sat.ma;
                    phaseAngle = mod(targetMA - chaserMA, 360);
                    if phaseAngle > 180, phaseAngle = 360 - phaseAngle; end
                    [dv_ph, tPh, pPh] = calculateSinglePhasingDV(...
                        phaseAngle, sat.a, MU_EARTH, selMass, THRUST, ISP, G0, MIN_SPIRAL_REVOLUTIONS, MAX_SPIRAL_REVOLUTIONS, THRUST_DUTY_CYCLE);
                    dv_ph = dv_ph * DV_MARGIN_FACTOR;
                    pPh = pPh * DV_MARGIN_FACTOR;
                end
                perSat(satCounter).dv_phasing   = dv_ph;
                perSat(satCounter).time_phasing  = tPh / 86400;
                perSat(satCounter).prop_phasing  = pPh;
                selMass = selMass - pPh;
                selTime = selTime + tPh;
                selTotalDV = selTotalDV + dv_ph;
                selCumDV = [selCumDV, selTotalDV];
                selCumTime = [selCumTime, selTime];
                selCumProp = [selCumProp, SPACECRAFT_MASS - selMass];
                massTimeline = [massTimeline; selTime, selMass];

                % Update chaser MA: after phasing, chaser is at the
                % target satellite's MA (they have rendezvoused).
                chaserMA = sat.ma;

                % Record arrival state (after all transfer maneuvers,
                % before operational costs at this satellite)
                satArrivalDV   = [satArrivalDV,   selTotalDV];
                satArrivalTime = [satArrivalTime,  selTime];
                satArrivalMass = [satArrivalMass,  selMass];

                % ---------------------------------------------------------
                % Operation cost at this satellite
                % ---------------------------------------------------------
                if hasOperations && satCounter <= length(OPERATION_DV_PER_SAT)
                    opDV   = OPERATION_DV_PER_SAT(satCounter);
                    opTime = OPERATION_TIME_DAYS_PER_SAT(satCounter);
                else
                    opDV   = 0;
                    opTime = 0;
                end
                if opDV > 0
                    opProp = selMass * (1 - 1/exp(opDV / v_e_sel));
                else
                    opProp = 0;
                end
                perSat(satCounter).dv_operation   = opDV;
                perSat(satCounter).time_operation  = opTime;
                perSat(satCounter).prop_operation  = opProp;
                if opDV > 0 || opTime > 0
                    selMass = selMass - opProp;
                    selTime = selTime + opTime * 86400;
                    selTotalDV = selTotalDV + opDV;
                    selCumDV  = [selCumDV, selTotalDV];
                    selCumTime = [selCumTime, selTime];
                    selCumProp = [selCumProp, SPACECRAFT_MASS - selMass];
                    massTimeline = [massTimeline; selTime, selMass];
                end

                % Cargo drop at this satellite (if scheduled)
                cargoDrop = getCargoDropMass(sat.norad, CARGO_DROP_SCHEDULE);
                perSat(satCounter).cargo_dropped = cargoDrop;
                if cargoDrop > 0
                    selMass = selMass - cargoDrop;
                    massTimeline = [massTimeline; selTime, selMass];
                end
            end  % satellites in plane

            % J2 RAAN change to next plane
            if planeVisit < nPlanes
                pcCounter = pcCounter + 1;
                nextPIdx = planeOrder(planeVisit + 1);
                nextPlane = planes(nextPIdx);
                fromSat = plane.satellites(sortOrder(end));
                toSat   = nextPlane.satellites(1);
                chosenAlt = J2_SELECTED_ALTITUDES_KM_PER_PC(pcCounter);
                [edDV_ref, pcDet] = calculateEdelbaumPlaneChange(...
                    fromSat, toSat, selMass, THRUST, ISP, G0, MU_EARTH, J2, R_EARTH, chosenAlt, THRUST_DUTY_CYCLE);
                j2opt = pcDet.j2Analysis.options(1);

                if j2opt.totalTime >= Inf
                    fprintf('  WARNING: J2 at %+.0f km INFEASIBLE for PC#%d (Plane %d->%d). Falling back to Edelbaum.\n',...
                        chosenAlt, pcCounter, pIdx, nextPIdx);
                    pcDV = edDV_ref * DV_MARGIN_FACTOR;
                    pcTime_s = pcDet.burnTime * DV_MARGIN_FACTOR;
                    pcProp = pcDet.propellant * DV_MARGIN_FACTOR;
                    pcDetailLog(pcCounter).strategy = 'Edelbaum fallback';
                    pcDetailLog(pcCounter).spiralTime = 0;
                    pcDetailLog(pcCounter).waitTime = 0;
                    pcDetailLog(pcCounter).spiralRAAN = 0;
                    pcDetailLog(pcCounter).direction = 'N/A';
                    chaserMA = propagateMeanAnomaly(chaserMA, a_plane, MU_EARTH, pcTime_s);
                else
                    pcDV = j2opt.dv_total;
                    pcTime_s = j2opt.totalTime * 86400;
                    pcProp = j2opt.propellant;
                    pcDetailLog(pcCounter).spiralTime = j2opt.spiralTime;
                    pcDetailLog(pcCounter).waitTime   = j2opt.waitTime;
                    pcDetailLog(pcCounter).spiralRAAN  = j2opt.spiralRAAN;
                    pcDetailLog(pcCounter).direction   = j2opt.direction;
                    pcDetailLog(pcCounter).strategy    = 'J2 precession';
                    a_new = a_plane + chosenAlt * 1e3;
                    t_spiral_s = j2opt.spiralTime * 86400;
                    t_wait_s   = j2opt.waitTime * 86400;
                    chaserMA = propagateMAThroughJ2Maneuver(...
                        chaserMA, a_plane, a_new, t_spiral_s, t_wait_s, MU_EARTH);
                end

                pcDetailLog(pcCounter).pcNumber   = pcCounter;
                pcDetailLog(pcCounter).fromPlane   = pIdx;
                pcDetailLog(pcCounter).toPlane     = nextPIdx;
                pcDetailLog(pcCounter).chosenAlt   = chosenAlt;
                pcDetailLog(pcCounter).dv          = pcDV;
                pcDetailLog(pcCounter).time_days   = pcTime_s / 86400;
                pcDetailLog(pcCounter).propellant  = pcProp;
                selMass = selMass - pcProp;
                selTime = selTime + pcTime_s;
                selTotalDV = selTotalDV + pcDV;
                selCumDV  = [selCumDV, selTotalDV];
                selCumTime = [selCumTime, selTime];
                selCumProp = [selCumProp, SPACECRAFT_MASS - selMass];
                massTimeline = [massTimeline; selTime, selMass];
            end
        end  % plane visit loop

        % Console output — widened table with Operation columns
        fprintf('\nPER-SATELLITE COST BREAKDOWN:\n');
        fprintf('%s\n', repmat('=', 1, 210));
        fprintf('%-4s %-30s %-6s %10s %10s %10s %10s %10s %12s %10s %10s %10s %10s %10s %10s\n',...
            '#','Satellite','Plane','Insert dV','RAAN dV','Inc dV','Phase dV','Oper dV','Total dV','Insert(d)','RAAN(d)','Inc(d)','Phase(d)','Oper(d)','Cargo(kg)');
        fprintf('%s\n', repmat('-', 1, 210));
        for k = 1:satCounter
            ps = perSat(k);
            totDV_k = ps.dv_insertion + ps.dv_raan + ps.dv_inc + ps.dv_phasing + ps.dv_operation;
            fprintf('%-4d %-30s %-6d %10.1f %10.1f %10.1f %10.1f %10.1f %12.1f %10.2f %10.2f %10.2f %10.2f %10.2f %10.1f\n',...
                ps.satID, truncateString(ps.name,30), ps.planeID,...
                ps.dv_insertion, ps.dv_raan, ps.dv_inc, ps.dv_phasing, ps.dv_operation, totDV_k,...
                ps.time_insertion, ps.time_raan, ps.time_inc, ps.time_phasing, ps.time_operation, ps.cargo_dropped);
        end
        fprintf('%s\n', repmat('-', 1, 210));
        totInsDV = sum([perSat.dv_insertion]); totRaanDV = sum([perSat.dv_raan]);
        totIncDV = sum([perSat.dv_inc]); totPhDV = sum([perSat.dv_phasing]);
        totOpDV  = sum([perSat.dv_operation]);
        totCargo = sum([perSat.cargo_dropped]);
        fprintf('%-4s %-30s %-6s %10.1f %10.1f %10.1f %10.1f %10.1f %12.1f %10.2f %10.2f %10.2f %10.2f %10.2f %10.1f\n',...
            '','TOTALS','', totInsDV, totRaanDV, totIncDV, totPhDV, totOpDV, selTotalDV,...
            sum([perSat.time_insertion]), sum([perSat.time_raan]), sum([perSat.time_inc]),...
            sum([perSat.time_phasing]), sum([perSat.time_operation]), totCargo);
        fprintf('%s\n', repmat('=', 1, 210));

        fprintf('\nDETAILED PLANE CHANGE BREAKDOWN (J2 Selected Altitudes):\n');
        fprintf('%s\n', repmat('-', 1, 130));
        fprintf('%-4s %-10s %-12s %-12s %-12s %-12s %-12s %-14s %-20s\n',...
            'PC#','Planes','Alt(km)','dV(m/s)','Time(days)','Spiral(d)','Wait(d)','RAAN Sprl(deg)','Direction');
        fprintf('%s\n', repmat('-', 1, 130));
        for k = 1:length(pcDetailLog)
            pd = pcDetailLog(k);
            fprintf('%-4d %-10s %-+12.0f %-12.1f %-12.1f %-12.1f %-12.1f %-+14.2f %-20s\n',...
                pd.pcNumber, sprintf('%d->%d',pd.fromPlane,pd.toPlane),...
                pd.chosenAlt, pd.dv, pd.time_days, pd.spiralTime,...
                pd.waitTime, pd.spiralRAAN, pd.direction);
        end
        fprintf('%s\n', repmat('-', 1, 130));

        fprintf('\nMISSION SUMMARY:\n');
        if totInsDV > 0
        fprintf('  Insertion dV:            %8.1f m/s  (%5.1f%%)   Prop: %8.2f kg\n', totInsDV, 100*totInsDV/selTotalDV, sum([perSat.prop_insertion]));
        end
        fprintf('  RAAN change total dV:    %8.1f m/s  (%5.1f%%)   Prop: %8.2f kg\n', totRaanDV, 100*totRaanDV/selTotalDV, sum([perSat.prop_raan]));
        fprintf('  Inclination total dV:    %8.1f m/s  (%5.1f%%)   Prop: %8.2f kg\n', totIncDV, 100*totIncDV/selTotalDV, sum([perSat.prop_inc]));
        fprintf('  Phasing total dV:        %8.1f m/s  (%5.1f%%)   Prop: %8.2f kg\n', totPhDV, 100*totPhDV/selTotalDV, sum([perSat.prop_phasing]));
        if totOpDV > 0 || (hasOperations && any(OPERATION_TIME_DAYS_PER_SAT > 0))
        fprintf('  Operation total dV:      %8.1f m/s  (%5.1f%%)   Prop: %8.2f kg\n', totOpDV, 100*totOpDV/selTotalDV, sum([perSat.prop_operation]));
        end
        fprintf('  TOTAL MISSION dV:        %8.1f m/s\n', selTotalDV);
        if totCargo > 0
        fprintf('  Total cargo dropped:     %8.1f kg  (%d deliveries)\n', totCargo, nnz([perSat.cargo_dropped]));
        end
        fprintf('  MISSION DURATION:        %8.1f days  (%.2f months, %.2f years)\n',...
            selTime/86400, selTime/86400/30, selTime/86400/365);
        fprintf('  FINAL MASS:              %8.2f kg  (%.1f%% propellant used)\n',...
            selMass, 100*(SPACECRAFT_MASS-selMass)/SPACECRAFT_MASS);
        fprintf('%s\n', repmat('=', 1, 130));

        plotSelectedJ2Mission(perSat, selCumDV, selCumTime, SPACECRAFT_MASS, J2_SELECTED_ALTITUDES_KM_PER_PC, nPlaneChanges, hasOperations);
        plotCumulativeDVSatArrivals(selCumDV, selCumTime, satArrivalDV, satArrivalTime, perSat, J2_SELECTED_ALTITUDES_KM_PER_PC, nPlaneChanges);
        plotMassTimeline(massTimeline, perSat, SPACECRAFT_MASS, satArrivalTime, satArrivalMass);
    end
end

fprintf('\n=======================================================================\n');
fprintf('Analysis complete.\n');
fprintf('=======================================================================\n');
diary off;

%% ==================== PLOTTING FUNCTIONS ====================

function plotMissionTimeline(cumulativeTime, cumulativeDV, ~, maneuverLog, totalPhasingDV, totalPlaneChangeDV, nPlanes, constellationFile)
    figure('Position', [150, 50, 1200, 800]);
    phasingMask = strcmp({maneuverLog.type}, 'Phasing');
    pcMask = strcmp({maneuverLog.type}, 'Plane Change');
    subplot(2,2,1);
    plot(cumulativeTime/(24*3600), cumulativeDV/1000, 'b-', 'LineWidth', 2); hold on;
    if any(pcMask), scatter(cumulativeTime(pcMask)/(24*3600), cumulativeDV(pcMask)/1000, 100, 'r', 'filled', 'MarkerEdgeColor', 'k'); end
    hold off; xlabel('Mission Time (days)'); ylabel('Cumulative Delta-V (km/s)');
    title('Delta-V Budget Over Time'); legend('Cumulative dV', 'Plane Changes', 'Location', 'northwest'); grid on;
    subplot(2,2,2);
    histogram([maneuverLog.dv], 20, 'FaceColor', [0.3, 0.6, 0.9]);
    xlabel('Delta-V (m/s)'); ylabel('Number of Maneuvers'); title('Distribution of Maneuver Delta-V'); grid on;
    subplot(2,2,3);
    dvByType = [totalPhasingDV, totalPlaneChangeDV] / 1000;
    bar(dvByType, 'FaceColor', 'flat', 'CData', [0.2,0.6,0.8; 0.8,0.4,0.2]);
    set(gca, 'XTickLabel', {'Phasing', 'Plane Changes'}); ylabel('Delta-V (km/s)'); title('Delta-V Breakdown'); grid on;
    subplot(2,2,4);
    dvPerPlane = zeros(nPlanes, 1);
    for ii = 1:length(maneuverLog)
        if strcmp(maneuverLog(ii).type, 'Phasing'), dvPerPlane(maneuverLog(ii).planeID) = dvPerPlane(maneuverLog(ii).planeID) + maneuverLog(ii).dv; end
    end
    bar(dvPerPlane/1000, 'FaceColor', [0.4, 0.7, 0.4]); xlabel('Plane Number'); ylabel('Phasing dV (km/s)');
    title('Phasing Delta-V by Orbital Plane'); grid on;
    sgtitle(sprintf('Constellation Analysis: %s (Edelbaum Model)', constellationFile), 'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, 'single_constellation_timeline.png');
end

function plotConstellation3D(planes, planeOrder, nPlanes, R_EARTH, MU_EARTH)
    figure('Position', [150, 50, 1200, 800]);
    [x, y, z] = sphere(50);
    surf(R_EARTH/1e3*x, R_EARTH/1e3*y, R_EARTH/1e3*z, 'FaceColor', [0.3,0.5,0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on; axis equal; xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
    title('3D Orbital Constellation View'); grid on; colors = lines(nPlanes);
    for p = 1:nPlanes
        plane = planes(planeOrder(p));
        for s = 1:length(plane.satellites)
            sat = plane.satellites(s);
            [r_vec, ~] = orbital_elements_to_cartesian(sat.a, sat.e, sat.i, sat.raan, sat.aop, sat.ma, MU_EARTH);
            scatter3(r_vec(1)/1e3, r_vec(2)/1e3, r_vec(3)/1e3, 100, colors(p,:), 'filled');
        end
    end
    view(45, 30); lighting gouraud; light('Position', [1, 1, 1]);
    saveas(gcf, 'single_constellation_3d.png');
end

function plotVisitSequence(satVisitOrder, satellites)
    figure('Position', [50, 100, 1400, 600]);
    subplot(1,2,1); scatter(1:length(satVisitOrder), satVisitOrder, 100, 'filled');
    hold on; plot(1:length(satVisitOrder), satVisitOrder, 'k--', 'LineWidth', 0.5); hold off;
    xlabel('Visit Order'); ylabel('Satellite ID'); title('Satellite Visit Sequence'); grid on;
    subplot(1,2,2); polarscatter(deg2rad([satellites.raan]), [satellites.ma], 100, 'filled');
    title('Satellite Distribution (RAAN vs Mean Anomaly)'); grid on;
    sgtitle('Visit Sequence Analysis', 'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, 'single_constellation_sequence.png');
end

function plotJ2TradeStudySingle(pcDet, mLog, fromPlaneIdx, toPlaneIdx, logTimeScale)
    if logTimeScale, scaleStr = 'log'; else, scaleStr = 'linear'; end
    j2opts = pcDet.j2Analysis.options; feasible = [j2opts.totalTime] < Inf;
    if ~any(feasible), return; end
    altVals = [j2opts(feasible).deltaAlt]; dvVals = [j2opts(feasible).dv_total];
    timeVals = [j2opts(feasible).totalTime]; propVals = [j2opts(feasible).propellant];
    waitVals = [j2opts(feasible).waitTime]; spiralVals = [j2opts(feasible).spiralTime];
    raanVals = [j2opts(feasible).spiralRAAN]; diffRates = [j2opts(feasible).differentialRate];
    negMask = altVals < 0; posMask = altVals > 0;
    activeDV = mLog.dv; activeTime = mLog.time / 86400; activeProp = mLog.propellant;
    figure('Position', [100, 50, 1400, 750]);
    sgtitle(sprintf('J2 Precession Trade Study: Plane %d -> %d  (dRAAN = %+.1f deg)', fromPlaneIdx, toPlaneIdx, pcDet.j2Analysis.raanChangeRequired), 'FontSize', 14, 'FontWeight', 'bold');
    subplot(2,3,1);
    if any(negMask), bar(altVals(negMask), dvVals(negMask), 0.4, 'FaceColor', [0.2,0.5,0.8]); hold on; end
    if any(posMask), bar(altVals(posMask), dvVals(posMask), 0.4, 'FaceColor', [0.85,0.4,0.2]); hold on; end
    yline(activeDV, 'r--', 'LineWidth', 2, 'Label', sprintf('Active: %.0f m/s', activeDV)); hold off;
    xlabel('Altitude Change (km)'); ylabel('Total Delta-V (m/s)'); title('Delta-V: J2 vs Active');
    legend('Lower orbit','Raise orbit','Location','best'); grid on;
    subplot(2,3,2);
    if any(negMask), bar(altVals(negMask), timeVals(negMask), 0.4, 'FaceColor', [0.2,0.5,0.8]); hold on; end
    if any(posMask), bar(altVals(posMask), timeVals(posMask), 0.4, 'FaceColor', [0.85,0.4,0.2]); hold on; end
    yline(activeTime, 'r--', 'LineWidth', 2, 'Label', sprintf('Active: %.0f d', activeTime)); hold off;
    if logTimeScale, set(gca, 'YScale', 'log'); end
    xlabel('Altitude Change (km)'); ylabel('Total Time (days)'); title(['Total Time (' scaleStr ')']); grid on;
    subplot(2,3,3);
    if any(negMask), bar(altVals(negMask), propVals(negMask), 0.4, 'FaceColor', [0.2,0.5,0.8]); hold on; end
    if any(posMask), bar(altVals(posMask), propVals(posMask), 0.4, 'FaceColor', [0.85,0.4,0.2]); hold on; end
    yline(activeProp, 'r--', 'LineWidth', 2, 'Label', sprintf('Active: %.1f kg', activeProp)); hold off;
    xlabel('Altitude Change (km)'); ylabel('Propellant (kg)'); title('Propellant'); grid on;
    subplot(2,3,4);
    bh = bar(altVals, [spiralVals', waitVals'], 'stacked');
    bh(1).FaceColor = [0.4,0.7,0.4]; bh(2).FaceColor = [0.95,0.75,0.3];
    if logTimeScale, set(gca, 'YScale', 'log'); end
    xlabel('Altitude Change (km)'); ylabel('Time (days)'); title(['Time Breakdown (' scaleStr ')']);
    legend('Spiral','Wait','Location','best'); grid on;
    subplot(2,3,5);
    yyaxis left; plot(altVals, raanVals, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    ylabel('RAAN from Spiral (deg)');
    yline(pcDet.j2Analysis.raanChangeRequired, 'b--', 'LineWidth', 1, 'Label', sprintf('Required: %+.1f', pcDet.j2Analysis.raanChangeRequired));
    yline(0, 'b:', 'LineWidth', 0.5);
    yyaxis right; plot(altVals, diffRates, 'rs-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    ylabel('Diff. Rate (deg/day)'); yline(0, 'r:', 'LineWidth', 0.5);
    xlabel('Altitude Change (km)'); title('Spiral RAAN & Diff. Rate'); grid on;
    subplot(2,3,6);
    if any(negMask), scatter(timeVals(negMask), dvVals(negMask), 80, [0.2,0.5,0.8], 'filled'); hold on;
        for k = find(negMask), text(timeVals(k), dvVals(k)+max(dvVals)*0.02, sprintf('%+.0f',altVals(k)), 'FontSize',7,'HorizontalAlignment','center','Color',[0.2,0.5,0.8]); end
    end
    if any(posMask), scatter(timeVals(posMask), dvVals(posMask), 80, [0.85,0.4,0.2], 'filled'); hold on;
        for k = find(posMask), text(timeVals(k), dvVals(k)+max(dvVals)*0.02, sprintf('%+.0f',altVals(k)), 'FontSize',7,'HorizontalAlignment','center','Color',[0.85,0.4,0.2]); end
    end
    scatter(activeTime, activeDV, 150, 'rp', 'filled', 'MarkerEdgeColor', 'k'); hold off;
    if logTimeScale, set(gca, 'YScale', 'log'); end
    xlabel('Total Time (days)'); ylabel('Delta-V (m/s)'); title(['dV vs Time (' scaleStr ')']);
    legend('Lower','Raise','Active','Location','best'); grid on;
    saveas(gcf, sprintf('j2_trade_plane%d_to_%d.png', fromPlaneIdx, toPlaneIdx));
    fprintf('J2 trade study plot saved: j2_trade_plane%d_to_%d.png\n', fromPlaneIdx, toPlaneIdx);
end

function plotJ2TradeStudyAggregate(planeChangeManeuvers, maneuverLog, planes, THRUST, ISP, G0, MU_EARTH, J2, R_EARTH, J2_ALT, logTimeScale, dutyCycle)
    if logTimeScale, scaleStr = 'log'; else, scaleStr = 'linear'; end
    nAlt = length(J2_ALT); nPC = length(planeChangeManeuvers);
    aggDV=zeros(1,nAlt); aggTime=zeros(1,nAlt); aggProp=zeros(1,nAlt);
    aggSpiralTime=zeros(1,nAlt); aggWaitTime=zeros(1,nAlt); aggSpiralRAAN=zeros(1,nAlt);
    aggDiffRate=zeros(1,nAlt); aggFeasible=true(1,nAlt);
    totalActiveDV=0; totalActiveTime=0; totalActiveProp=0; totalRAAN_req=0;
    for pcIdx = 1:nPC
        mIdx = planeChangeManeuvers(pcIdx);
        totalActiveDV = totalActiveDV + maneuverLog(mIdx).dv;
        totalActiveTime = totalActiveTime + maneuverLog(mIdx).time/86400;
        totalActiveProp = totalActiveProp + maneuverLog(mIdx).propellant;
        fPI = maneuverLog(mIdx).fromPlaneID; tPI = maneuverLog(mIdx).toPlaneID;
        [~, pcDet] = calculateEdelbaumPlaneChange(planes(fPI).satellites(end), planes(tPI).satellites(1), maneuverLog(mIdx).massBefore, THRUST, ISP, G0, MU_EARTH, J2, R_EARTH, J2_ALT, dutyCycle);
        totalRAAN_req = totalRAAN_req + pcDet.j2Analysis.raanChangeRequired;
        j2o = pcDet.j2Analysis.options;
        for k = 1:nAlt
            if j2o(k).totalTime >= Inf, aggFeasible(k)=false;
            else
                aggDV(k)=aggDV(k)+j2o(k).dv_total; aggTime(k)=aggTime(k)+j2o(k).totalTime;
                aggProp(k)=aggProp(k)+j2o(k).propellant; aggSpiralTime(k)=aggSpiralTime(k)+j2o(k).spiralTime;
                aggWaitTime(k)=aggWaitTime(k)+j2o(k).waitTime; aggSpiralRAAN(k)=aggSpiralRAAN(k)+j2o(k).spiralRAAN;
                aggDiffRate(k)=aggDiffRate(k)+j2o(k).differentialRate;
            end
        end
    end
    aggDiffRate(aggFeasible) = aggDiffRate(aggFeasible)/nPC;
    aggDV(~aggFeasible)=NaN; aggTime(~aggFeasible)=NaN; aggProp(~aggFeasible)=NaN;
    aggSpiralTime(~aggFeasible)=NaN; aggWaitTime(~aggFeasible)=NaN;
    aggSpiralRAAN(~aggFeasible)=NaN; aggDiffRate(~aggFeasible)=NaN;
    negM = J2_ALT<0 & aggFeasible; posM = J2_ALT>0 & aggFeasible;
    figure('Position', [80,40,1400,750]);
    sgtitle(sprintf('Aggregated J2 Trade — %d Plane Changes (Total dRAAN=%+.1f deg)', nPC, totalRAAN_req), 'FontSize', 14, 'FontWeight', 'bold');
    subplot(2,3,1);
    if any(negM), bar(J2_ALT(negM),aggDV(negM),0.4,'FaceColor',[0.2,0.5,0.8]); hold on; end
    if any(posM), bar(J2_ALT(posM),aggDV(posM),0.4,'FaceColor',[0.85,0.4,0.2]); hold on; end
    yline(totalActiveDV,'r--','LineWidth',2,'Label',sprintf('Active: %.0f',totalActiveDV)); hold off;
    xlabel('Alt Change (km)'); ylabel('dV (m/s)'); title('Delta-V'); legend('Lower','Raise','Location','best'); grid on;
    subplot(2,3,2);
    if any(negM), bar(J2_ALT(negM),aggTime(negM),0.4,'FaceColor',[0.2,0.5,0.8]); hold on; end
    if any(posM), bar(J2_ALT(posM),aggTime(posM),0.4,'FaceColor',[0.85,0.4,0.2]); hold on; end
    yline(totalActiveTime,'r--','LineWidth',2,'Label',sprintf('Active: %.0f d',totalActiveTime)); hold off;
    if logTimeScale, set(gca,'YScale','log'); end
    xlabel('Alt Change (km)'); ylabel('Time (days)'); title(['Time (' scaleStr ')']); grid on;
    subplot(2,3,3);
    if any(negM), bar(J2_ALT(negM),aggProp(negM),0.4,'FaceColor',[0.2,0.5,0.8]); hold on; end
    if any(posM), bar(J2_ALT(posM),aggProp(posM),0.4,'FaceColor',[0.85,0.4,0.2]); hold on; end
    yline(totalActiveProp,'r--','LineWidth',2,'Label',sprintf('Active: %.1f kg',totalActiveProp)); hold off;
    xlabel('Alt Change (km)'); ylabel('Propellant (kg)'); title('Propellant'); grid on;
    subplot(2,3,4);
    fI = find(aggFeasible); bh=bar(J2_ALT(fI),[aggSpiralTime(fI)',aggWaitTime(fI)'],'stacked');
    bh(1).FaceColor=[0.4,0.7,0.4]; bh(2).FaceColor=[0.95,0.75,0.3];
    if logTimeScale, set(gca,'YScale','log'); end
    xlabel('Alt Change (km)'); ylabel('Time (days)'); title(['Breakdown (' scaleStr ')']); legend('Spiral','Wait','Location','best'); grid on;
    subplot(2,3,5);
    yyaxis left; plot(J2_ALT(aggFeasible),aggSpiralRAAN(aggFeasible),'bo-','LineWidth',1.5,'MarkerFaceColor','b');
    ylabel('Sum Spiral RAAN (deg)');
    yline(totalRAAN_req,'b--','LineWidth',1,'Label',sprintf('Required: %+.1f',totalRAAN_req)); yline(0,'b:','LineWidth',0.5);
    yyaxis right; plot(J2_ALT(aggFeasible),aggDiffRate(aggFeasible),'rs-','LineWidth',1.5,'MarkerFaceColor','r');
    ylabel('Avg Diff. Rate (deg/day)'); yline(0,'r:','LineWidth',0.5);
    xlabel('Alt Change (km)'); title('RAAN & Rate'); grid on;
    subplot(2,3,6);
    if any(negM), scatter(aggTime(negM),aggDV(negM),80,[0.2,0.5,0.8],'filled'); hold on;
        for k=find(negM), text(aggTime(k),aggDV(k)+max(aggDV(aggFeasible))*0.02,sprintf('%+.0f',J2_ALT(k)),'FontSize',7,'HorizontalAlignment','center','Color',[0.2,0.5,0.8]); end
    end
    if any(posM), scatter(aggTime(posM),aggDV(posM),80,[0.85,0.4,0.2],'filled'); hold on;
        for k=find(posM), text(aggTime(k),aggDV(k)+max(aggDV(aggFeasible))*0.02,sprintf('%+.0f',J2_ALT(k)),'FontSize',7,'HorizontalAlignment','center','Color',[0.85,0.4,0.2]); end
    end
    scatter(totalActiveTime,totalActiveDV,150,'rp','filled','MarkerEdgeColor','k'); hold off;
    if logTimeScale, set(gca,'YScale','log'); end
    xlabel('Time (days)'); ylabel('dV (m/s)'); title(['Pareto (' scaleStr ')']); legend('Lower','Raise','Active','Location','best'); grid on;
    saveas(gcf, 'j2_trade_aggregate.png'); fprintf('Aggregated J2 trade study plot saved: j2_trade_aggregate.png\n');
    fI = find(aggFeasible);
    if ~isempty(fI)
        [~,bDVi]=min(aggDV(fI)); [~,bTi]=min(aggTime(fI));
        fprintf('\nAGGREGATED BEST J2 OPTIONS:\n');
        fprintf('  Active total: %.0f m/s, %.1f days, %.1f kg\n', totalActiveDV, totalActiveTime, totalActiveProp);
        fprintf('  Best dV  (%+.0f km): %.0f m/s, %.1f d, %.1f kg [%.1f%% dV savings]\n', J2_ALT(fI(bDVi)), aggDV(fI(bDVi)), aggTime(fI(bDVi)), aggProp(fI(bDVi)), 100*(1-aggDV(fI(bDVi))/totalActiveDV));
        fprintf('  Best time (%+.0f km): %.0f m/s, %.1f d, %.1f kg [%.1f%% dV savings]\n', J2_ALT(fI(bTi)), aggDV(fI(bTi)), aggTime(fI(bTi)), aggProp(fI(bTi)), 100*(1-aggDV(fI(bTi))/totalActiveDV));
    end
end

function plotSelectedJ2Mission(perSat, selCumDV, selCumTime, SPACECRAFT_MASS, selectedAlts, nPlaneChanges, hasOperations)
    nS = length(perSat);
    dvIns  = [perSat.dv_insertion];
    dvRaan = [perSat.dv_raan]; dvInc = [perSat.dv_inc]; dvPh = [perSat.dv_phasing];
    dvOp   = [perSat.dv_operation];
    tRaan  = [perSat.time_raan]; tInc = [perSat.time_inc]; tPh = [perSat.time_phasing];
    tIns   = [perSat.time_insertion]; tOp = [perSat.time_operation];
    pIns   = [perSat.prop_insertion];
    pRaan  = [perSat.prop_raan]; pInc = [perSat.prop_inc]; pPh = [perSat.prop_phasing];
    pOp    = [perSat.prop_operation];
    dvTotal = dvIns + dvRaan + dvInc + dvPh + dvOp;

    % X-axis labels
    xLabels = cell(1, nS);
    for k = 1:nS
        rawName = perSat(k).name;
        if length(rawName) > 18, rawName = [rawName(1:15), '...']; end
        xLabels{k} = rawName;
    end

    altStr = '[';
    for k = 1:nPlaneChanges
        altStr = [altStr, sprintf('%+.0f', selectedAlts(k))];
        if k < nPlaneChanges, altStr = [altStr, ', ']; end
    end
    altStr = [altStr, '] km'];

    % Determine which optional categories are present
    hasInsertion  = any(dvIns > 0) || any(pIns > 0) || any(tIns > 0);
    hasOperation  = hasOperations && (any(dvOp > 0) || any(pOp > 0) || any(tOp > 0));

    % Build stacked data and colors dynamically
    stackDV   = []; stackTime = []; stackProp = [];
    barColors = {}; legLabels = {};

    if hasInsertion
        stackDV   = [stackDV,   dvIns'];
        stackTime = [stackTime, tIns'];
        stackProp = [stackProp, pIns'];
        barColors{end+1} = [0.15, 0.70, 0.15];
        legLabels{end+1} = 'Insertion';
    end
    % Always include core categories
    stackDV   = [stackDV,   dvRaan', dvInc', dvPh'];
    stackTime = [stackTime, tRaan',  tInc',  tPh'];
    stackProp = [stackProp, pRaan',  pInc',  pPh'];
    barColors = [barColors, {[0.85,0.35,0.15], [0.55,0.20,0.70], [0.20,0.60,0.85]}];
    legLabels = [legLabels, {'RAAN', 'Inclination', 'Phasing'}];

    if hasOperation
        stackDV   = [stackDV,   dvOp'];
        stackTime = [stackTime, tOp'];
        stackProp = [stackProp, pOp'];
        barColors{end+1} = [0.95, 0.75, 0.15];
        legLabels{end+1} = 'Operation';
    end
    nCat = length(barColors);

    figure('Position', [40, 20, 1200, 800]);
    sgtitle(sprintf('Selected J2 Mission — Per-Satellite Cost Breakdown | Altitudes: %s', altStr),...
        'FontSize', 14, 'FontWeight', 'bold');

    % Panel 1: Cumulative dV over time
    subplot(2,3,1);
    plot(selCumTime/86400, selCumDV/1000, 'b-', 'LineWidth', 2); hold on;
    raanSats = find(dvRaan > 0);
    for rs = raanSats
        cumDV_at = sum(dvTotal(1:rs));
        [~, bI] = min(abs(selCumDV - cumDV_at));
        scatter(selCumTime(bI)/86400, selCumDV(bI)/1000, 120, 'r', 'filled', 'MarkerEdgeColor', 'k');
    end
    hold off; xlabel('Time (days)'); ylabel('Cumulative dV (km/s)');
    title('Cumulative Delta-V'); grid on;
    if ~isempty(raanSats), legend('J2 Mission','Plane Changes','Location','northwest'); end

    % Panel 2: Per-satellite dV stacked bar
    subplot(2,3,2);
    bh = bar(1:nS, stackDV, 'stacked');
    for c = 1:nCat, bh(c).FaceColor = barColors{c}; end
    set(gca, 'XTick', 1:nS, 'XTickLabel', xLabels, 'XTickLabelRotation', 60, 'FontSize', 7);
    ylabel('Delta-V (m/s)'); title('Per-Satellite Delta-V');
    legend(legLabels{:}, 'Location', 'best'); grid on;

    % Panel 3: Per-satellite time stacked bar
    subplot(2,3,3);
    bh = bar(1:nS, stackTime, 'stacked');
    for c = 1:nCat, bh(c).FaceColor = barColors{c}; end
    set(gca, 'XTick', 1:nS, 'XTickLabel', xLabels, 'XTickLabelRotation', 60, 'FontSize', 7);
    ylabel('Time (days)'); title('Per-Satellite Time');
    legend(legLabels{:}, 'Location', 'best'); grid on;

    % Panel 4: Per-satellite propellant stacked bar
    subplot(2,3,4);
    bh = bar(1:nS, stackProp, 'stacked');
    for c = 1:nCat, bh(c).FaceColor = barColors{c}; end
    set(gca, 'XTick', 1:nS, 'XTickLabel', xLabels, 'XTickLabelRotation', 60, 'FontSize', 7);
    ylabel('Propellant (kg)'); title('Per-Satellite Propellant');
    legend(legLabels{:}, 'Location', 'best'); grid on;

    % Panel 5: Pie chart — total dV by category
    subplot(2,3,5);
    totI = sum(dvIns); totR = sum(dvRaan); totInc = sum(dvInc); totP = sum(dvPh); totO = sum(dvOp);
    totAll = totI + totR + totInc + totP + totO;

    pieVals  = []; pieLabels = {}; pieColors = [];
    if hasInsertion
        pieVals(end+1)  = totI;
        pieLabels{end+1} = sprintf('Insert: %.0f m/s (%.0f%%)', totI, 100*totI/totAll);
        pieColors = [pieColors; 0.15,0.70,0.15];
    end
    pieVals  = [pieVals, totR, totInc, totP];
    pieLabels = [pieLabels, ...
        {sprintf('RAAN: %.0f m/s (%.0f%%)', totR, 100*totR/totAll)}, ...
        {sprintf('Inc: %.0f m/s (%.0f%%)', totInc, 100*totInc/totAll)}, ...
        {sprintf('Phasing: %.0f m/s (%.0f%%)', totP, 100*totP/totAll)}];
    pieColors = [pieColors; 0.85,0.35,0.15; 0.55,0.20,0.70; 0.20,0.60,0.85];
    if hasOperation
        pieVals(end+1)  = totO;
        pieLabels{end+1} = sprintf('Oper: %.0f m/s (%.0f%%)', totO, 100*totO/totAll);
        pieColors = [pieColors; 0.95,0.75,0.15];
    end
    % Remove zero-valued slices
    nonzero = pieVals > 0;
    pie(pieVals(nonzero), pieLabels(nonzero));
    title('Delta-V Breakdown by Type');
    colormap(gca, pieColors(nonzero,:));

    % Panel 6: Per-satellite total dV bar (colored by plane)
    subplot(2,3,6);
    planeIDs = [perSat.planeID]; uniquePlanes = unique(planeIDs, 'stable');
    nUP = length(uniquePlanes); planeColors = lines(nUP);
    bh = bar(1:nS, dvTotal, 'FaceColor', 'flat');
    cdata = zeros(nS, 3);
    for k = 1:nS
        pci = find(uniquePlanes == planeIDs(k), 1);
        cdata(k,:) = planeColors(pci,:);
    end
    bh.CData = cdata;
    set(gca, 'XTick', 1:nS, 'XTickLabel', xLabels, 'XTickLabelRotation', 60, 'FontSize', 7);
    ylabel('Total dV (m/s)'); title('Total dV per Satellite (colored by plane)');
    hold on;
    for pp = 1:nUP
        bar(NaN, NaN, 'FaceColor', planeColors(pp,:), 'DisplayName', sprintf('Plane %d', uniquePlanes(pp)));
    end
    hold off; legend('Location','best'); grid on;

    saveas(gcf, 'selected_j2_mission_analysis.png');
    fprintf('Selected J2 mission analysis plot saved: selected_j2_mission_analysis.png\n');
end

%% ===================== HELPER FUNCTIONS =====================

function str = truncateString(inputStr, maxLen)
    if length(inputStr) > maxLen, str = [inputStr(1:maxLen-3), '...']; else, str = inputStr; end
end

function satellites = loadJCATFormat(filename, tleCacheHours)
    if nargin < 2, tleCacheHours = 4; end
    opts = detectImportOptions(filename); opts.VariableNamingRule = 'preserve';
    for i = 1:length(opts.VariableTypes)
        if strcmp(opts.VariableTypes{i}, 'datetime'), opts.VariableTypes{i} = 'char'; end
    end
    data = readtable(filename, opts); colNames = data.Properties.VariableNames;
    nameCol = findColumn(colNames, {'Name','name','TF'});
    jcatCol = findColumn(colNames, {'#JCAT','JCAT','NORAD','CatalogNumber','TOp'});
    periCol = findColumn(colNames, {'DispPeri','UNPerigee'});
    apoCol = findColumn(colNames, {'DispApo','UNApogee'});
    incCol = findColumn(colNames, {'DispInc','UNInc'});
    raanCol = findColumn(colNames, {'RAAN','raan','DispRAAN','UNRaan'});
    nSats = height(data);
    satellites = struct('name',cell(1,nSats),'a',cell(1,nSats),'e',cell(1,nSats),'i',cell(1,nSats),'raan',cell(1,nSats),'aop',cell(1,nSats),'ma',cell(1,nSats),'norad',cell(1,nSats));
    noradIDs = []; validRows = [];
    for k = 1:nSats
        peri = getNumericValue(data, k, periCol); apo = getNumericValue(data, k, apoCol); inc = getNumericValue(data, k, incCol);
        if isnan(peri) || isnan(apo) || isnan(inc) || peri <= 0 || apo <= 0, continue; end
        if ~isempty(jcatCol)
            jcatVal = data{k, jcatCol}; if iscell(jcatVal), jcatStr = jcatVal{1}; else, jcatStr = char(jcatVal); end
            noradID = extractNoradID(jcatStr);
        else, noradID = NaN; end
        if ~isnan(noradID), noradIDs(end+1) = noradID; validRows(end+1) = k; end
    end
    tleData = fetchTLEsFromCelesTrak(noradIDs, tleCacheHours);
    validCount = 0;
    for idx = 1:length(validRows)
        k = validRows(idx); noradID = noradIDs(idx); validCount = validCount + 1;
        if ~isempty(nameCol)
            nameVal = data{k, nameCol};
            if iscell(nameVal), satellites(validCount).name = nameVal{1}; else, satellites(validCount).name = char(nameVal); end
        else, satellites(validCount).name = sprintf('SAT-%d', noradID); end
        satellites(validCount).name = strrep(strrep(satellites(validCount).name, '[', ''), ']', '');
        satellites(validCount).norad = noradID;
        if isfield(tleData, sprintf('n%d', noradID))
            tle = tleData.(sprintf('n%d', noradID));
            satellites(validCount).a = tle.a; satellites(validCount).e = tle.e;
            satellites(validCount).i = tle.i; satellites(validCount).raan = tle.raan;
            satellites(validCount).aop = tle.aop; satellites(validCount).ma = tle.ma;
        else
            peri = getNumericValue(data, k, periCol); apo = getNumericValue(data, k, apoCol); inc = getNumericValue(data, k, incCol);
            R_EARTH_KM = 6378.137; r_peri = (peri + R_EARTH_KM) * 1e3; r_apo = (apo + R_EARTH_KM) * 1e3;
            satellites(validCount).a = (r_peri + r_apo) / 2; satellites(validCount).e = (r_apo - r_peri) / (r_apo + r_peri);
            satellites(validCount).i = inc;
            if ~isempty(raanCol), satellites(validCount).raan = getNumericValue(data, k, raanCol); else, satellites(validCount).raan = 0; end
            satellites(validCount).aop = 0; satellites(validCount).ma = mod(validCount * 45, 360);
        end
    end
    satellites = satellites(1:validCount);
end

function noradID = extractNoradID(jcatStr)
    jcatStr = strtrim(char(jcatStr));
    if ~isempty(jcatStr) && upper(jcatStr(1)) == 'S', jcatStr = jcatStr(2:end); end
    numStr = regexp(jcatStr, '\d+', 'match', 'once');
    if ~isempty(numStr), noradID = str2double(numStr); else, noradID = NaN; end
end

function tleData = fetchTLEsFromCelesTrak(noradIDs, tleCacheHours)
    if nargin < 2, tleCacheHours = 4; end
    cacheFile = 'tle_cache.txt'; tleData = struct(); nSats = length(noradIDs);
    cacheValid = false; cachedTLEs = struct(); originalCacheTimestamp = '';
    if isfile(cacheFile)
        fid = fopen(cacheFile, 'r');
        if fid ~= -1
            firstLine = fgetl(fid); fclose(fid);
            if startsWith(firstLine, 'CACHE_DATE:')
                cacheDateStr = strtrim(firstLine(12:end)); originalCacheTimestamp = cacheDateStr;
                try
                    cacheDateTime = datenum(cacheDateStr, 'yyyy-mm-dd HH:MM:SS');
                    hoursSinceCache = (now - cacheDateTime) * 24;
                    if hoursSinceCache < tleCacheHours
                        cacheValid = true; cachedTLEs = loadTLECache(cacheFile);
                        fprintf('Using cached TLE data (%.1f hours old)\n', hoursSinceCache);
                    else, fprintf('Cache expired (%.1f hours old), fetching new TLE data...\n', hoursSinceCache); end
                catch, fprintf('Invalid cache date, fetching new TLE data...\n'); end
            end
        end
    end
    if cacheValid
        missingIDs = [];
        for i = 1:nSats, if ~isfield(cachedTLEs, sprintf('n%d', noradIDs(i))), missingIDs(end+1) = noradIDs(i); end; end
        tleData = cachedTLEs;
    else, missingIDs = noradIDs; end
    if ~isempty(missingIDs)
        options = weboptions('Timeout', 30);
        for i = 1:length(missingIDs)
            noradID = missingIDs(i);
            url = sprintf('https://celestrak.org/NORAD/elements/gp.php?CATNR=%d&FORMAT=TLE', noradID);
            fprintf('Fetching TLE for NORAD ID %d\n', noradID);
            try
                tleText = webread(url, options); lines = strsplit(tleText, '\n');
                if length(lines) >= 3
                    line0 = strtrim(lines{1}); line1 = strtrim(lines{2}); line2 = strtrim(lines{3});
                    if length(line1) >= 69 && line1(1) == '1' && length(line2) >= 69 && line2(1) == '2'
                        inc = str2double(line2(9:16)); raan = str2double(line2(18:25));
                        ecc = str2double(['0.' line2(27:33)]); aop = str2double(line2(35:42));
                        ma = str2double(line2(44:51)); n = str2double(line2(53:63));
                        n_rad_s = n * 2 * pi / 86400; a = (3.986004418e14 / n_rad_s^2)^(1/3);
                        tleData.(sprintf('n%d', noradID)) = struct('a',a,'e',ecc,'i',inc,'raan',raan,'aop',aop,'ma',ma,'name',line0,'line1',line1,'line2',line2);
                    end
                end
            catch, end
            if i < length(missingIDs), pause(0.2); end
        end
        if cacheValid && ~isempty(originalCacheTimestamp), saveTLECacheWithTimestamp(cacheFile, tleData, originalCacheTimestamp);
        else, saveTLECacheWithTimestamp(cacheFile, tleData, datestr(now, 'yyyy-mm-dd HH:MM:SS')); end
    end
end

function cachedTLEs = loadTLECache(cacheFile)
    cachedTLEs = struct(); fid = fopen(cacheFile, 'r'); if fid == -1, return; end
    fgetl(fid);
    while ~feof(fid)
        line0 = fgetl(fid);
        if ~ischar(line0) || isempty(strtrim(line0)) || startsWith(line0, '---'), continue; end
        line1 = fgetl(fid); line2 = fgetl(fid);
        if ~ischar(line1) || ~ischar(line2), break; end
        line0 = strtrim(line0); line1 = strtrim(line1); line2 = strtrim(line2);
        if length(line1) >= 69 && line1(1) == '1' && length(line2) >= 69 && line2(1) == '2'
            noradID = str2double(line1(3:7));
            inc = str2double(line2(9:16)); raan = str2double(line2(18:25));
            ecc = str2double(['0.' line2(27:33)]); aop = str2double(line2(35:42));
            ma = str2double(line2(44:51)); n = str2double(line2(53:63));
            n_rad_s = n * 2 * pi / 86400; a = (3.986004418e14 / n_rad_s^2)^(1/3);
            cachedTLEs.(sprintf('n%d', noradID)) = struct('a',a,'e',ecc,'i',inc,'raan',raan,'aop',aop,'ma',ma,'name',line0,'line1',line1,'line2',line2);
        end
    end
    fclose(fid);
end

function saveTLECacheWithTimestamp(cacheFile, tleData, timestamp)
    fid = fopen(cacheFile, 'w'); if fid == -1, return; end
    fprintf(fid, 'CACHE_DATE: %s\n', timestamp);
    fprintf(fid, '-----------------------------------------------------------\n');
    fields = fieldnames(tleData);
    for i = 1:length(fields)
        tle = tleData.(fields{i});
        if isfield(tle, 'line1') && isfield(tle, 'line2'), fprintf(fid, '%s\n%s\n%s\n', tle.name, tle.line1, tle.line2); end
    end
    fclose(fid);
end

function col = findColumn(colNames, possibleNames)
    col = [];
    for i = 1:length(possibleNames), idx = find(strcmpi(colNames, possibleNames{i}), 1); if ~isempty(idx), col = idx; return; end; end
end

function val = getNumericValue(data, row, col)
    if isempty(col), val = NaN; return; end
    rawVal = data{row, col}; if iscell(rawVal), rawVal = rawVal{1}; end
    if isnumeric(rawVal), val = rawVal;
    elseif ischar(rawVal) || isstring(rawVal)
        rawVal = strtrim(char(rawVal)); rawVal = strrep(strrep(rawVal, ',', ''), ' ', '');
        if isempty(rawVal) || strcmp(rawVal, '-') || strcmp(rawVal, '*'), val = NaN; else, val = str2double(rawVal); end
    else, val = NaN; end
end

function planes = groupIntoPlanes(satellites, raanTolerance)
% Clusters satellites into orbital planes based on RAAN proximity.
% Walks through sats sorted by RAAN; if a sat's RAAN is within tolerance
% of an existing plane (and inclination is close), it gets lumped in.
% Otherwise it starts a new plane. Plane RAAN/inc are the averages of
% all members once grouping is done.
    nSats = length(satellites); raans = [satellites.raan]; [~, sortIdx] = sort(raans);
    planes = struct('raan', {}, 'inc', {}, 'satellites', {}); assigned = false(1, nSats);
    for i = 1:nSats
        idx = sortIdx(i); if assigned(idx), continue; end
        planeIdx = length(planes) + 1;
        planes(planeIdx).raan = satellites(idx).raan; planes(planeIdx).inc = satellites(idx).i;
        planes(planeIdx).satellites = satellites(idx); assigned(idx) = true;
        for j = i+1:nSats
            jdx = sortIdx(j); if assigned(jdx), continue; end
            raanDiff = abs(satellites(jdx).raan - planes(planeIdx).raan); raanDiff = min(raanDiff, 360 - raanDiff);
            incDiff = abs(satellites(jdx).i - planes(planeIdx).inc);
            if raanDiff <= raanTolerance && incDiff <= 2.0, planes(planeIdx).satellites(end+1) = satellites(jdx); assigned(jdx) = true; end
        end
        % update the plane's RAAN and inc to the group average
        planes(planeIdx).raan = mean([planes(planeIdx).satellites.raan]);
        planes(planeIdx).inc = mean([planes(planeIdx).satellites.i]);
    end
end

function [bestOrder, bestCost] = solvePlaneTSP(planes)
% Finds the cheapest order to visit all orbital planes (like a traveling
% salesman problem but for RAAN/inclination jumps). For small numbers of
% planes it uses an exact dynamic programming solver (Held-Karp), for
% bigger constellations it falls back to greedy + 2-opt local search.
% Held-Karp reduces complexity from O(n!) to O(n^2 2^n)
    nPlanes = length(planes);
    if nPlanes <= 1, bestOrder = 1:nPlanes; bestCost = 0; return; end
    if nPlanes == 2, bestOrder = [1, 2]; bestCost = calculatePlaneTransitionCost(planes(1), planes(2)); return; end
    costMatrix = buildCostMatrix(planes);
    if nPlanes <= 15, [bestOrder, bestCost] = heldKarpTSP(costMatrix);
    else, [bestOrder, bestCost] = greedyTwoOptTSP(costMatrix); end
end

function costMatrix = buildCostMatrix(planes)
% Builds an NxN cost matrix where each entry is the "expense" of jumping
% from one plane to another (based on RAAN gap + weighted inclination diff).
    nPlanes = length(planes); costMatrix = zeros(nPlanes);
    for i = 1:nPlanes, for j = i+1:nPlanes
        cost = calculatePlaneTransitionCost(planes(i), planes(j)); costMatrix(i,j) = cost; costMatrix(j,i) = cost;
    end; end
end

function cost = calculatePlaneTransitionCost(plane1, plane2)
% Simple heuristic cost: RAAN difference + 3x inclination difference.
% The 3x weight on inclination reflects how much more expensive inc
% changes are compared to RAAN changes (inc burns are brutal).
    mu = 3.986004418e14;
    a_avg = (mean([plane1.satellites.a]) + mean([plane2.satellites.a])) / 2;
    v_avg = sqrt(mu / a_avg);
    i1_rad = deg2rad(plane1.inc);
    i2_rad = deg2rad(plane2.inc);
    i_avg_rad = (i1_rad + i2_rad) / 2;
    dRAAN = plane2.raan - plane1.raan;
    if dRAAN > 180, dRAAN = dRAAN - 360; end
    if dRAAN < -180, dRAAN = dRAAN + 360; end
    dRAAN_rad = deg2rad(abs(dRAAN));
    dv_raan = 2 * v_avg * sin(i_avg_rad) * sin(dRAAN_rad / 2);
    dv_inc = v_avg * abs(i2_rad - i1_rad);
    cost = sqrt(dv_raan^2 + dv_inc^2);
end

function [bestOrder, bestCost] = heldKarpTSP(costMatrix)
% Exact TSP solver using Held-Karp dynamic programming. Uses bitmasks
% to track which planes have been visited. Guaranteed optimal but
% exponential in memory/time, so only use for <= ~15 planes.
    n = size(costMatrix, 1); numSubsets = 2^n;

    % dp(mask, last) = cheapest cost to visit the set of planes encoded
    % in 'mask', ending at plane 'last'
    dp = inf(numSubsets, n); parent = zeros(numSubsets, n);

    % base case: starting at each plane individually costs 0
    for i = 1:n, dp(bitshift(1, i-1) + 1, i) = 0; end

    % fill out DP table and try extending each partial tour by one more plane
    for mask = 1:numSubsets-1
        for last = 1:n
            if bitand(mask, bitshift(1, last-1)) == 0 || isinf(dp(mask + 1, last)), continue; end
            for next = 1:n
                if bitand(mask, bitshift(1, next-1)) ~= 0, continue; end
                newMask = bitor(mask, bitshift(1, next-1));
                newCost = dp(mask + 1, last) + costMatrix(last, next);
                if newCost < dp(newMask + 1, next), dp(newMask + 1, next) = newCost; parent(newMask + 1, next) = last; end
            end
        end
    end

    % backtrack from the full bitmask to reconstruct the best visit order
    fullMask = numSubsets - 1; bestCost = inf; bestEnd = 1;
    for i = 1:n, if dp(fullMask + 1, i) < bestCost, bestCost = dp(fullMask + 1, i); bestEnd = i; end; end
    bestOrder = zeros(1, n); mask = fullMask; current = bestEnd;
    for pos = n:-1:1, bestOrder(pos) = current; prev = parent(mask + 1, current); mask = bitxor(mask, bitshift(1, current-1)); current = prev; end
end

function [bestOrder, bestCost] = greedyTwoOptTSP(costMatrix)
% Greedy nearest-neighbor construction followed by 2-opt edge swaps.
% Not guaranteed optimal but runs fast for large plane counts.
% Tries every possible starting node and keeps the best result.
    n = size(costMatrix, 1); bestOrder = []; bestCost = inf;
    for startNode = 1:n
        order = greedyConstruct(costMatrix, startNode); cost = calculateTourCost(costMatrix, order);
        [order, cost] = twoOptImprove(costMatrix, order, cost);
        if cost < bestCost, bestCost = cost; bestOrder = order; end
    end
end

function order = greedyConstruct(costMatrix, startNode)
% Nearest-neighbor heuristic: always go to the cheapest unvisited plane.
    n = size(costMatrix, 1); visited = false(1, n); order = zeros(1, n);
    current = startNode; order(1) = current; visited(current) = true;
    for i = 2:n
        minCost = inf; nextNode = -1;
        for j = 1:n, if ~visited(j) && costMatrix(current, j) < minCost, minCost = costMatrix(current, j); nextNode = j; end; end
        order(i) = nextNode; visited(nextNode) = true; current = nextNode;
    end
end

function cost = calculateTourCost(costMatrix, order)
    cost = 0; for i = 1:length(order)-1, cost = cost + costMatrix(order(i), order(i+1)); end
end

function [bestOrder, bestCost] = twoOptImprove(costMatrix, order, currentCost)
% Classic 2-opt: try reversing every sub-segment of the tour.
% If flipping a segment reduces cost, keep it. Repeat until no
% improvement is found (local optimum).
    n = length(order); improved = true; bestOrder = order; bestCost = currentCost;
    while improved, improved = false;
        for i = 1:n-2, for j = i+2:n
            oldCost = costMatrix(order(i), order(i+1));
            if j < n, oldCost = oldCost + costMatrix(order(j), order(j+1)); end
            newCost = costMatrix(order(i), order(j));
            if j < n, newCost = newCost + costMatrix(order(i+1), order(j+1)); end
            if newCost < oldCost - 1e-10
                bestOrder = [order(1:i), fliplr(order(i+1:j)), order(j+1:end)];
                order = bestOrder; bestCost = bestCost - oldCost + newCost; improved = true;
            end
        end; end
    end
end

function [totalDV, visitOrder, details] = calculateContinuousThrustPhasingDV(satellites, currentMass, thrust, isp, g0, mu, minRevs, maxRevs, dutyCycle)
% Figures out the dV to visit every satellite in one orbital plane using
% low-thrust phasing maneuvers. Sorts them by mean anomaly and hops from
% one to the next. Each hop uses a slightly different orbit period to
% drift and catch up to the next sat.
%
% dutyCycle scales wall-clock time: thruster only fires this fraction
% of the time, so everything takes longer in real life.
    nSats = length(satellites);
    if nSats <= 1, totalDV = 0; visitOrder = 1;
        details = struct('dv',{},'from',{},'to',{},'phaseAngle',{},'revolutions',{},'burnTime',{},'propellant',{}); return; end

    % sort by mean anomaly so we visit them in order around the orbit
    mas = [satellites.ma]; [~, visitOrder] = sort(mas);
    totalDV = 0; v_e = isp * g0;
    details = struct('dv',{},'from',{},'to',{},'phaseAngle',{},'revolutions',{},'burnTime',{},'propellant',{});
    mass = currentMass;
    for i = 1:nSats-1
        currentSat = satellites(visitOrder(i)); nextSat = satellites(visitOrder(i+1));

        % how far apart are these two sats in the orbit (in degrees)
        phaseAngle = nextSat.ma - currentSat.ma;
        if phaseAngle < 0, phaseAngle = phaseAngle + 360; end
        if phaseAngle > 180, phaseAngle = 360 - phaseAngle; end

        a_avg = (currentSat.a + nextSat.a) / 2;
        v_circ = sqrt(mu / a_avg); T_orbit = 2 * pi * sqrt(a_avg^3 / mu);

        if phaseAngle < 5
            % really close together, just give it a small nudge
            dv = 10; nRevs = minRevs;
        else
            % pick how many revolutions to spend drifting based on phase angle.
            % bigger angle = fewer revs needed (you're already close-ish in the
            % "wrong" direction, so a faster transfer works)
            nRevs = max(minRevs, min(maxRevs, round(minRevs + (maxRevs - minRevs) * (1 - phaseAngle/180))));
            phaseRad = deg2rad(phaseAngle);

            % the servicer flies a slightly different orbit period so it
            % completes a different number of orbits in the same elapsed time,
            % ending up at the target's position
            if phaseAngle <= 90, servicerOrbits = nRevs + phaseRad/(2*pi);
            else, servicerOrbits = nRevs - (2*pi - phaseRad)/(2*pi); end
            t_transfer = nRevs * T_orbit; T_servicer = t_transfer / servicerOrbits;

            % back out what SMA gives us that period (Kepler's 3rd law)
            a_servicer = (mu * (T_servicer/(2*pi))^2)^(1/3);

            % safety check: make sure the phasing orbit doesn't dip into
            % the atmosphere. if it does, double the revs to use a gentler orbit
            if (2*a_servicer - a_avg) < 6478e3
                nRevs = nRevs * 2; servicerOrbits = nRevs + phaseRad/(2*pi);
                T_servicer = (nRevs * T_orbit) / servicerOrbits;
                a_servicer = (mu * (T_servicer/(2*pi))^2)^(1/3);
            end

            % dV is twice the Hohmann-ish burn (raise then lower back),
            % with a 0.9 fudge factor since continuous thrust is a bit
            % more efficient than impulsive
            v_servicer = sqrt(mu * (2/a_avg - 1/a_servicer));
            dv_per_burn = abs(v_circ - v_servicer); dv = 2 * dv_per_burn * 0.9;

            % if the two sats are at noticeably different altitudes,
            % RSS in a Hohmann-style altitude correction
            alt_diff = abs(currentSat.a - nextSat.a);
            if alt_diff > 100, dv = sqrt(dv^2 + (abs(sqrt(mu/nextSat.a) - sqrt(mu/currentSat.a)))^2); end
        end
        dv = max(dv, 10);  % floor at 10 m/s for numerical sanity

        % rocket equation for propellant (doesn't depend on duty cycle
        % same total impulse regardless of when you fire)
        massRatio = exp(dv / v_e); propellant = mass * (1 - 1/massRatio);

        % wall-clock burn time: thrust-on time divided by duty cycle
        % (if you can only fire 50% of the time, it takes 2x as long)
        burnTime = mass * v_e / (thrust * dutyCycle) * (1 - 1/massRatio);
        totalDV = totalDV + dv;
        details(i).dv = dv; details(i).from = currentSat.name; details(i).to = nextSat.name;
        details(i).phaseAngle = phaseAngle; details(i).revolutions = nRevs;
        details(i).burnTime = burnTime * 2.5; details(i).propellant = propellant;
        mass = mass - propellant;
    end
end

function [dv, details] = calculateEdelbaumPlaneChange(sat1, sat2, currentMass, thrust,...
        isp, g0, mu, J2, R_EARTH, altChangeOptionsKm, dutyCycle)
% Computes the Edelbaum continuous-thrust dV to change orbital planes
% (RAAN and/or inclination) between two satellites. Also kicks off the
% J2 precession analysis as an alternative strategy.
%
% The key idea: Edelbaum says you can rotate the orbital plane by
% thrusting continuously along the velocity direction with a yaw bias.
% The dV formula is essentially 2*V*sin(i)*sin(dRAAN/2) for RAAN and
% V*di for inclination, then RSS them together.
%
% dutyCycle scales wall-clock burn time (thruster off during eclipse etc.)
    v_e = isp * g0; i1 = sat1.i; i2 = sat2.i;

    % figure out the signed RAAN change, wrapping to [-180, 180]
    raanChange_signed = sat2.raan - sat1.raan;
    if raanChange_signed > 180, raanChange_signed = raanChange_signed - 360; end
    if raanChange_signed < -180, raanChange_signed = raanChange_signed + 360; end
    raanChange_abs = abs(raanChange_signed); raanChange_rad = deg2rad(raanChange_abs);
    incChange_abs = abs(i2 - i1); incChange_rad = deg2rad(incChange_abs);

    % compute the "combined angle" the actual 3D rotation between
    % the two orbital planes using the spherical trig formula
    i1_rad = deg2rad(i1); i2_rad = deg2rad(i2);
    dRAAN_rad = deg2rad(raanChange_signed);
    cosTheta = cos(i1_rad)*cos(i2_rad) + sin(i1_rad)*sin(i2_rad)*cos(dRAAN_rad);
    combinedAngle = rad2deg(acos(max(-1, min(1, cosTheta))));

    a_avg = (sat1.a + sat2.a) / 2; v_avg = sqrt(mu / a_avg);
    altDiff = abs(sat1.a - sat2.a); i_avg_rad = deg2rad((i1 + i2) / 2);

    % Edelbaum RAAN dV: thrust along velocity with a yaw steering program
    if raanChange_abs > 0.1, dv_raan = 2 * v_avg * sin(i_avg_rad) * sin(raanChange_rad / 2);
        else, dv_raan = 0; end

    % Inclination dV: small-angle approximation (just V * delta_i)
    if incChange_abs > 0.01, dv_inc = v_avg * incChange_rad; else, dv_inc = 0; end

    % RSS the two components since they're roughly orthogonal
    dv_plane = sqrt(dv_raan^2 + dv_inc^2);

    % label the strategy for the output log
    if raanChange_abs < 5 && incChange_abs < 5, strategy = sprintf(...
            'Small (dRAAN=%.2f, di=%.2f) - Edelbaum', raanChange_abs, incChange_abs);
    elseif incChange_abs < 1.0, strategy = sprintf('RAAN dominant (%.2f deg) - Edelbaum', raanChange_abs);
    else, strategy = sprintf('Combined RAAN (%.2f) + inc (%.2f) - Edelbaum', raanChange_abs, incChange_abs); end

    % if the sats are at different altitudes, tack on a Hohmann-ish cost
    dv_altitude = 0;
    if altDiff > 1e3, dv_altitude = abs(sqrt(mu/sat2.a) - sqrt(mu/sat1.a)); end
    if altDiff < 1e3, dv = dv_plane;
    else, dv = sqrt(dv_plane^2 + dv_altitude^2 * 0.8); strategy = [strategy, ' + alt change']; end

    % wall-clock burn time: divide by duty cycle since we can't fire
    % during eclipse. dV and propellant don't change same total impulse.
    burnTime = dv * currentMass / (thrust * dutyCycle);
    massRatio = exp(dv / v_e); propellant = currentMass * (1 - 1/massRatio);
    T_orbit = 2 * pi * sqrt(a_avg^3 / mu); revolutions = burnTime / T_orbit;

    % now run the J2 precession trade study as an alternative
    j2Analysis = analyzeJ2Precession_Integrated(sat1, sat2, raanChange_signed, (i1+i2)/2, a_avg, currentMass, thrust, isp, g0, mu, J2, R_EARTH, altChangeOptionsKm, dutyCycle);
    details.raanChange = raanChange_abs; details.incChange = incChange_abs;
    details.combinedAngle = combinedAngle; details.dv_raan = dv_raan;
    details.dv_inc = dv_inc; details.altitudeChange = altDiff / 1e3;
    details.strategy = strategy; details.propellant = propellant;
    details.burnTime = burnTime; details.revolutions = revolutions;
    details.j2Analysis = j2Analysis;
end

function j2Analysis = analyzeJ2Precession_Integrated(sat1, sat2, raanChange_signed,...
        i_avg, a_nominal, currentMass, thrust, isp, g0, mu, J2, R_EARTH, altChangeOptionsKm, dutyCycle)
% The "free lunch" alternative to active RAAN changes: instead of
% burning to rotate the orbit, we spiral up/down to a different altitude
% where J2 precesses the RAAN at a different rate, wait for the RAAN
% to drift to where we want it, then spiral back.
%
% The tradeoff: way less propellant (only need dV for the altitude
% change, not the plane rotation), but potentially way more time
% (waiting for differential precession to accumulate).
%
% Uses the integrated form of the RAAN accumulated during the spiral
% itself (not just during the coast), since the ion thruster spirals
% are slow enough that significant precession happens en route.
%
% dutyCycle affects spiral wall-clock time AND the RAAN accumulated
% during spirals: if the thruster is only on half the time, the spiral
% takes 2x as long, and precession keeps ticking the whole time (it
% doesn't care if the thruster is on or off). Mathematically the RAAN
% integral picks up a 1/dutyCycle factor because da/dt is slower but
% dOmega/dt stays the same.
    v_e = isp * g0; i_rad = deg2rad(i_avg);

    % J2 RAAN precession rate as a function of SMA (the a^(-7/2) dependence)
    raanDotFunc = @(a) -1.5 * sqrt(mu/a^3) * J2 * (R_EARTH/a)^2 * cos(i_rad);
    raanDot_nominal = raanDotFunc(a_nominal);
    raanDot_nominal_degPerDay = rad2deg(raanDot_nominal) * 86400;

    altChanges_km = sort(altChangeOptionsKm); nOptions = length(altChanges_km);
    j2Analysis.options = struct('deltaAlt',{},'newAlt',{},'dv_total',{},'spiralTime',{},...
        'spiralRAAN',{},'waitTime',{},'totalTime',{},'propellant',{},'differentialRate',{},'direction',{},'description',{});

    for idx = 1:nOptions
        deltaAlt_m = altChanges_km(idx) * 1e3; a_new = a_nominal + deltaAlt_m;

        % don't let the orbit dip below ~150 km that's atmospheric drag city
        if a_new < (R_EARTH + 150e3)
            j2Analysis.options(idx) = struct('deltaAlt',altChanges_km(idx),'newAlt',...
                (a_new-R_EARTH)/1e3,'dv_total',Inf,'spiralTime',Inf,'spiralRAAN',0,'waitTime',Inf,'totalTime',Inf,'propellant',Inf,'differentialRate',0,'direction','Infeasible','description','Infeasible');
            continue;
        end

        % dV to change altitude: just the vis-viva velocity difference
        % (one way), doubled for the round trip (go up, come back down)
        v_nominal = sqrt(mu / a_nominal); v_new = sqrt(mu / a_new);
        dv_oneWay = abs(v_new - v_nominal); dv_total = 2 * dv_oneWay;

        % rocket equation for total propellant (independent of duty cycle)
        massRatio_total = exp(dv_total / v_e); propellant = currentMass * (1 - 1/massRatio_total);

        % wall-clock spiral times: divide by duty cycle because the thruster
        % is only on for a fraction of each orbit, so the altitude change
        % takes proportionally longer in real elapsed time
        t_spiral_up = dv_oneWay * currentMass / (thrust * dutyCycle);
        massRatio_oneWay = exp(dv_oneWay / v_e); mass_after_up = currentMass / massRatio_oneWay;
        t_spiral_down = dv_oneWay * mass_after_up / (thrust * dutyCycle);
        t_spiral_total = t_spiral_up + t_spiral_down;

        % here's the cool part: RAAN accumulated DURING the spiral itself.
        % This is the integrated form: dRAAN = (3/16)*J2*RE^2*mu*cos(i)*(m/F)*(af^-4 - a0^-4)
        % It comes from integrating the precession rate over the slow spiral
        % where da/dt = 2a^2*F/(mu*m) and accounting for the changing SMA.
        %
        % With duty cycle: the effective da/dt slows by eta, so
        % the spacecraft spends 1/eta times longer at each intermediate
        % altitude. Precession doesn't care about the thruster, it just
        % keeps ticking so the integral picks up a 1/dutyCycle factor:
        %   coeff * (m / (F * dutyCycle)) * (af^-4 - a0^-4)
        coeff = (3.0/16.0) * J2 * R_EARTH^2 * mu * cos(i_rad);
        if a_new > a_nominal, sign_leg1 = +1; sign_leg2 = -1; else, sign_leg1 = -1; sign_leg2 = +1; end

        % leg 1: spiral out (or in), leg 2: spiral back
        % note the (m / (thrust * dutyCycle)) that's the key duty cycle change
        dRAAN_leg1_rad = sign_leg1 * coeff * (currentMass / (thrust * dutyCycle)) * (a_new^(-4) - a_nominal^(-4));
        dRAAN_leg2_rad = sign_leg2 * coeff * (mass_after_up / (thrust * dutyCycle)) * (a_nominal^(-4) - a_new^(-4));
        dRAAN_spiral_total_rad = dRAAN_leg1_rad + dRAAN_leg2_rad;

        % subtract out what would have precessed at the nominal altitude
        % anyway during the same wall-clock time (we want the NET extra precession)
        dRAAN_nominal_during_spiral_rad = raanDot_nominal * t_spiral_total;
        spiralRAAN_net_deg = rad2deg(dRAAN_spiral_total_rad - dRAAN_nominal_during_spiral_rad);

        % differential precession rate at the new altitude vs nominal
        raanDot_new_degPerDay = rad2deg(raanDotFunc(a_new)) * 86400;
        signedDiffRate = raanDot_new_degPerDay - raanDot_nominal_degPerDay;

        % how much RAAN do we still need after the spiral gets us partway there?
        remainingRAAN_signed = mod(raanChange_signed - spiralRAAN_net_deg + 180, 360) - 180;

        % figure out the coast/wait time at the new altitude
        % (thruster is off during wait, but precession keeps going no
        %  duty cycle correction needed here since we're just coasting)
        if abs(remainingRAAN_signed) < 0.1
            % the spiral alone did the job no waiting needed
            waitTime = 0; dirStr = 'Correct (spiral sufficient)';
        elseif abs(signedDiffRate) < 1e-6
            % diff rate is basically zero we'd wait forever
            waitTime = Inf; dirStr = 'No differential rate';
        elseif sign(signedDiffRate) == sign(remainingRAAN_signed)
            % the diff rate is pushing us the right way, just divide
            waitTime = abs(remainingRAAN_signed) / abs(signedDiffRate) * 86400; dirStr = 'Correct';
        else
            % uh oh, diff rate is going the WRONG direction. We have to
            % go the long way around (360 - remaining) degrees
            longWay = 360 - abs(remainingRAAN_signed);
            waitTime = longWay / abs(signedDiffRate) * 86400; dirStr = sprintf('Wrong (%.0f deg)', longWay);
        end
        totalTime_s = t_spiral_total + waitTime; totalTime_days = totalTime_s / 86400;
        if totalTime_s < Inf, descStr = sprintf('%+.0f km: %.0f m/s, %.0f d (%s)', altChanges_km(idx), dv_total, totalTime_days, dirStr);
        else, descStr = 'Infeasible'; end
        j2Analysis.options(idx) = struct('deltaAlt',altChanges_km(idx),'newAlt',(a_new-R_EARTH)/1e3,'dv_total',dv_total,'spiralTime',t_spiral_total/86400,'spiralRAAN',spiralRAAN_net_deg,'waitTime',waitTime/86400,'totalTime',totalTime_days,'propellant',propellant,'differentialRate',signedDiffRate,'direction',dirStr,'description',descStr);
    end
    j2Analysis.nominalPrecessionRate = raanDot_nominal_degPerDay;
    j2Analysis.raanChangeRequired = raanChange_signed;
    j2Analysis.nominalAltitude = (a_nominal - R_EARTH) / 1e3;
end

function [r_vec, v_vec] = orbital_elements_to_cartesian(a, e, i, raan, aop, ma, mu)
% Classic Keplerian elements -> ECI position & velocity.
% Steps: mean anomaly -> eccentric anomaly (Newton iteration) ->
% true anomaly -> perifocal frame -> rotate by AoP, inc, RAAN into ECI.
    i = deg2rad(i); raan = deg2rad(raan); aop = deg2rad(aop); ma = deg2rad(ma);

    % solve Kepler's equation M = E - e*sin(E) iteratively
    E = ma; for iter = 1:10, E = ma + e * sin(E); end

    % eccentric anomaly -> true anomaly
    nu = 2 * atan2(sqrt(1 + e) * sin(E/2), sqrt(1 - e) * cos(E/2));

    % position in the perifocal (orbit-plane) frame
    r = a * (1 - e^2) / (1 + e * cos(nu));
    r_orb = r * [cos(nu); sin(nu); 0]; h = sqrt(mu * a * (1 - e^2));
    v_orb = (mu / h) * [-sin(nu); e + cos(nu); 0];

    % 3-1-3 rotation sequence (RAAN, inc, AoP) to get from perifocal to ECI
    R3_raan = [cos(raan),-sin(raan),0; sin(raan),cos(raan),0; 0,0,1];
    R1_i = [1,0,0; 0,cos(i),-sin(i); 0,sin(i),cos(i)];
    R3_aop = [cos(aop),-sin(aop),0; sin(aop),cos(aop),0; 0,0,1];
    Q = R3_raan * R1_i * R3_aop; r_vec = Q * r_orb; v_vec = Q * v_orb;
end

function ma_new = propagateMeanAnomaly(ma0, a, mu, dt)
% Just advances the mean anomaly forward in time. Mean anomaly increases
% linearly at the mean motion rate n = sqrt(mu/a^3). Simple but important
% for tracking where the chaser ends up after long maneuvers.
    n = sqrt(mu / a^3);           % mean motion (rad/s)
    ma_new = mod(ma0 + rad2deg(n * dt), 360);
end

function chaserMA = propagateMAThroughJ2Maneuver(ma0, a_nom, a_new, t_spiral_s, t_wait_s, mu)
% Propagates the chaser's mean anomaly through the full J2 altitude
% maneuver: spiral up (at average SMA), coast at the new altitude,
% then spiral back down. Each phase has a different mean motion since
% the SMA is different, so the chaser drifts at different rates.
    t_half = t_spiral_s / 2;
    a_avg = (a_nom + a_new) / 2;

    % different mean motions at each altitude
    n_avg  = sqrt(mu / a_avg^3);
    n_new  = sqrt(mu / a_new^3);

    % leg 1: spiral up (half the total spiral time, at average altitude)
    ma_after_spiral1 = ma0 + rad2deg(n_avg * t_half);
    % coast at new altitude
    ma_after_wait    = ma_after_spiral1 + rad2deg(n_new * t_wait_s);
    % leg 2: spiral back down
    ma_after_spiral2 = ma_after_wait + rad2deg(n_avg * t_half);
    chaserMA = mod(ma_after_spiral2, 360);
end

function [dv, burnTime, propellant] = calculateSmallIncChange(inc_plane, inc_sat, a, mu, mass, thrust, isp, g0, dutyCycle)
% Inclination correction for when a satellite's actual inclination
% doesn't perfectly match its plane average.
%
% Uses the exact single-impulse plane change formula:
%   dV = 2 * V * sin(di/2)
% instead of the small-angle approximation (V * di). These agree for
% tiny angles but diverge noticeably above ~5 deg. The exact form is
% what you'd actually fly for a dedicated inclination burn at the
% ascending/descending node.
    di_rad = deg2rad(abs(inc_sat - inc_plane));
    if di_rad < 1e-9
        dv = 0; burnTime = 0; propellant = 0; return;
    end
    v = sqrt(mu / a);
    % exact plane change: burn at the node, split velocity triangle
    dv = 2 * v * sin(di_rad / 2);
    v_e = isp * g0;
    massRatio = exp(dv / v_e);
    propellant = mass * (1 - 1/massRatio);
    % wall-clock time accounts for duty cycle
    burnTime = dv * mass / (thrust * dutyCycle);
end

function [dv, burnTime, propellant] = calculateSinglePhasingDV(phaseAngle, a, mu, mass, thrust, isp, g0, minRevs, maxRevs, dutyCycle)
% Same idea as calculateContinuousThrustPhasingDV but for a single hop.
% Given a phase angle to close, figures out the dV by computing a drift
% orbit with a slightly different period. The servicer flies N revolutions
% while the target flies N +/- a fraction, and they meet up.
    if phaseAngle < 0.5
        dv = 0; burnTime = 0; propellant = 0; return;
    end
    v_e = isp * g0;
    v_circ = sqrt(mu / a);
    T_orbit = 2 * pi * sqrt(a^3 / mu);

    if phaseAngle < 5
        % tiny angle just give it a minimum dV kick
        dv = 10; nRevs = minRevs;
    else
        % pick number of drift revs (more revs = gentler transfer)
        nRevs = max(minRevs, min(maxRevs, round(minRevs + (maxRevs - minRevs) * (1 - phaseAngle/180))));
        phaseRad = deg2rad(phaseAngle);

        % servicer completes a slightly different number of orbits
        % so it ends up at the target's position after nRevs target orbits
        if phaseAngle <= 90
            servicerOrbits = nRevs + phaseRad / (2*pi);
        else
            servicerOrbits = nRevs - (2*pi - phaseRad) / (2*pi);
        end
        t_transfer = nRevs * T_orbit;
        T_servicer = t_transfer / servicerOrbits;

        % Kepler's 3rd law gives us the SMA for that period
        a_servicer = (mu * (T_servicer / (2*pi))^2)^(1/3);

        % make sure we don't clip the atmosphere
        if (2 * a_servicer - a) < 6478e3
            nRevs = nRevs * 2;
            servicerOrbits = nRevs + phaseRad / (2*pi);
            T_servicer = (nRevs * T_orbit) / servicerOrbits;
            a_servicer = (mu * (T_servicer / (2*pi))^2)^(1/3);
        end

        % dV is twice the burn (enter and exit the phasing orbit)
        v_servicer = sqrt(mu * (2/a - 1/a_servicer));
        dv_per_burn = abs(v_circ - v_servicer);
        dv = 2 * dv_per_burn * 0.9;  % 0.9 fudge for continuous thrust efficiency
    end
    dv = max(dv, 5);  % floor at 5 m/s

    % rocket equation (propellant is duty-cycle independent)
    massRatio = exp(dv / v_e);
    propellant = mass * (1 - 1/massRatio);
    % wall-clock burn time: divide by duty cycle, then 2.5x for coast arcs
    burnTime = mass * v_e / (thrust * dutyCycle) * (1 - 1/massRatio);
    burnTime = burnTime * 2.5;  % account for coast arcs
end

function cargoMass = getCargoDropMass(satNoradID, cargoSchedule)
% Simple lookup: checks if we're supposed to drop cargo at this particular
% satellite. Returns the mass to drop, or 0 if nothing scheduled here.
    cargoMass = 0;
    if isempty(cargoSchedule), return; end
    idx = find(cargoSchedule(:,1) == satNoradID, 1);
    if ~isempty(idx)
        cargoMass = cargoSchedule(idx, 2);
    end
end

function plotCumulativeDVSatArrivals(selCumDV, selCumTime, satArrivalDV, satArrivalTime, perSat, selectedAlts, nPlaneChanges)
% Cumulative delta-V vs mission time with a red circle at every satellite
% arrival.  Arrival is defined as the moment after all transfer maneuvers
% (RAAN, inclination, phasing) have completed but BEFORE operational costs
% are incurred at the client.
    nS = length(perSat);

    % Build altitude string for the title
    altStr = '[';
    for k = 1:nPlaneChanges
        altStr = [altStr, sprintf('%+.0f', selectedAlts(k))];
        if k < nPlaneChanges, altStr = [altStr, ', ']; end
    end
    altStr = [altStr, '] km'];

    figure('Position', [80, 80, 1100, 550]);
    plot(selCumTime/86400, selCumDV/1000, 'b-', 'LineWidth', 2); hold on;

    % Red circles at every satellite arrival
    for k = 1:nS
        scatter(satArrivalTime(k)/86400, satArrivalDV(k)/1000, 100, 'r', 'filled', 'MarkerEdgeColor', 'k');
    end

    % Label each arrival with truncated satellite name
    for k = 1:nS
        rawName = perSat(k).name;
        if length(rawName) > 15, rawName = [rawName(1:12), '...']; end
        text(satArrivalTime(k)/86400, satArrivalDV(k)/1000 + max(selCumDV/1000)*0.02, ...
            rawName, 'FontSize', 12, 'Rotation', 75, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
            'Color', [0.5, 0, 0]);
    end

    hold off;
    xlabel('Mission Time (days)'); ylabel('Cumulative Delta-V (km/s)');
    title(sprintf('Cumulative Delta-V with Satellite Arrivals | Altitudes: %s', altStr));
    legend('Cumulative dV', 'Satellite Arrival', 'Location', 'northwest');
    grid on;

    saveas(gcf, 'cumulative_dv_sat_arrivals.png');
    fprintf('Cumulative dV with satellite arrivals plot saved: cumulative_dv_sat_arrivals.png\n');
end

function plotMassTimeline(massTimeline, perSat, SPACECRAFT_MASS, satArrivalTime, satArrivalMass)
% Plot chaser mass vs mission time.
%   massTimeline    - Nx2 array [time_s, mass_kg]
%   perSat          - per-satellite struct array (for cargo drop annotations)
%   SPACECRAFT_MASS - initial mass (kg)
%   satArrivalTime  - 1xN vector, mission time (s) at each satellite arrival
%   satArrivalMass  - 1xN vector, chaser mass (kg) at each satellite arrival
    figure('Position', [100, 100, 1000, 500]);
    timeDays = massTimeline(:,1) / 86400;
    massKg   = massTimeline(:,2);
    plot(timeDays, massKg, 'b-', 'LineWidth', 2); hold on;

    % =====================================================================
    % SATELLITE ARRIVAL MARKERS
    % Red circles at the moment each satellite is reached (after phasing,
    % before operational costs).  Comment out or delete this block to
    % remove the circles from the figure.
    % =====================================================================
    hasArrivalMarkers = false;
    if nargin >= 5 && ~isempty(satArrivalTime) && ~isempty(satArrivalMass)
        scatter(satArrivalTime/86400, satArrivalMass, 80, 'r', 'filled', 'MarkerEdgeColor', 'k');
        hasArrivalMarkers = true;
    end
    % =====================================================================
    % END SATELLITE ARRIVAL MARKERS
    % =====================================================================

    % Cargo drop annotations (red triangles)
    cargoSats = find([perSat.cargo_dropped] > 0);
    if ~isempty(cargoSats)
        for cs = cargoSats
            ps = perSat(cs);
            for r = 2:size(massTimeline,1)
                if massTimeline(r,1) == massTimeline(r-1,1) && ...
                   abs((massTimeline(r-1,2) - massTimeline(r,2)) - ps.cargo_dropped) < 0.5
                    scatter(timeDays(r), massKg(r), 120, 'rv', 'filled', 'MarkerEdgeColor', 'k');
                    text(timeDays(r), massKg(r) - SPACECRAFT_MASS*0.01, ...
                        sprintf('  %s\n  -%.0f kg', truncateString(ps.name, 15), ps.cargo_dropped), ...
                        'FontSize', 7, 'Color', [0.7,0,0], 'VerticalAlignment', 'top');
                    break;
                end
            end
        end
    end

    hold off;
    xlabel('Mission Time (days)'); ylabel('Chaser Mass (kg)');
    title('Spacecraft Mass Over Mission Duration');
    grid on;
    ylim([min(massKg)*0.95, SPACECRAFT_MASS*1.02]);

    % Build legend dynamically based on what's plotted
    legEntries = {'Mass'};
    if hasArrivalMarkers, legEntries{end+1} = 'Satellite Arrival'; end
    if ~isempty(cargoSats), legEntries{end+1} = 'Cargo Drop'; end
    if length(legEntries) > 1
        legend(legEntries{:}, 'Location', 'best');
    end
    saveas(gcf, 'mass_timeline.png');
    fprintf('Mass timeline plot saved: mass_timeline.png\n');
end