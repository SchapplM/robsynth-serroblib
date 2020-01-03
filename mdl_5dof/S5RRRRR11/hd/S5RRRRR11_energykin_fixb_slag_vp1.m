% Calculate kinetic energy for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR11_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR11_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:46
% EndTime: 2019-12-31 22:38:48
% DurationCPUTime: 1.78s
% Computational Cost: add. (2034->306), mult. (4641->486), div. (0->0), fcn. (5662->12), ass. (0->139)
t423 = sin(qJ(2));
t424 = sin(qJ(1));
t426 = cos(qJ(2));
t427 = cos(qJ(1));
t449 = cos(pkin(5));
t436 = t427 * t449;
t400 = t423 * t424 - t426 * t436;
t401 = t423 * t436 + t424 * t426;
t420 = sin(pkin(5));
t444 = t427 * t420;
t364 = Icges(3,5) * t401 - Icges(3,6) * t400 - Icges(3,3) * t444;
t437 = t424 * t449;
t402 = t427 * t423 + t426 * t437;
t403 = -t423 * t437 + t427 * t426;
t445 = t424 * t420;
t365 = Icges(3,5) * t403 - Icges(3,6) * t402 + Icges(3,3) * t445;
t454 = (t364 * t427 - t365 * t424) * t420;
t452 = cos(qJ(3));
t425 = cos(qJ(4));
t451 = pkin(4) * t425;
t421 = sin(qJ(4));
t448 = t400 * t421;
t447 = t402 * t421;
t446 = t420 * t426;
t376 = pkin(2) * t401 + pkin(8) * t400;
t377 = pkin(2) * t403 + pkin(8) * t402;
t441 = qJD(2) * t420;
t412 = t424 * t441;
t438 = t427 * t441;
t443 = t376 * t412 + t377 * t438;
t389 = qJD(3) * t402 + t412;
t442 = qJD(1) * (pkin(1) * t424 - pkin(7) * t444);
t413 = qJD(2) * t449 + qJD(1);
t440 = t421 * t446;
t422 = sin(qJ(3));
t439 = t420 * t452;
t387 = t403 * t422 - t424 * t439;
t349 = qJD(4) * t387 + t389;
t390 = qJD(3) * t400 - t438;
t385 = t401 * t422 + t427 * t439;
t386 = t401 * t452 - t422 * t444;
t347 = pkin(3) * t386 + pkin(9) * t385;
t388 = t403 * t452 + t422 * t445;
t348 = pkin(3) * t388 + pkin(9) * t387;
t434 = t389 * t347 - t348 * t390 + t443;
t350 = qJD(4) * t385 + t390;
t405 = -qJD(3) * t446 + t413;
t404 = (pkin(2) * t423 - pkin(8) * t426) * t420;
t406 = qJD(1) * (pkin(1) * t427 + pkin(7) * t445);
t433 = t413 * t377 - t404 * t412 + t406;
t398 = t420 * t422 * t423 - t449 * t452;
t378 = qJD(4) * t398 + t405;
t432 = -t376 * t413 - t404 * t438 - t442;
t399 = t449 * t422 + t423 * t439;
t375 = pkin(3) * t399 + pkin(9) * t398;
t431 = t405 * t348 - t375 * t389 + t433;
t430 = -t347 * t405 + t390 * t375 + t432;
t419 = qJ(4) + qJ(5);
t418 = cos(t419);
t417 = sin(t419);
t409 = rSges(2,1) * t427 - rSges(2,2) * t424;
t408 = rSges(2,1) * t424 + rSges(2,2) * t427;
t394 = t449 * rSges(3,3) + (rSges(3,1) * t423 + rSges(3,2) * t426) * t420;
t393 = Icges(3,5) * t449 + (Icges(3,1) * t423 + Icges(3,4) * t426) * t420;
t392 = Icges(3,6) * t449 + (Icges(3,4) * t423 + Icges(3,2) * t426) * t420;
t391 = Icges(3,3) * t449 + (Icges(3,5) * t423 + Icges(3,6) * t426) * t420;
t384 = t399 * t425 - t440;
t383 = -t399 * t421 - t425 * t446;
t380 = t399 * t418 - t417 * t446;
t379 = -t399 * t417 - t418 * t446;
t372 = rSges(3,1) * t403 - rSges(3,2) * t402 + rSges(3,3) * t445;
t371 = rSges(3,1) * t401 - rSges(3,2) * t400 - rSges(3,3) * t444;
t369 = Icges(3,1) * t403 - Icges(3,4) * t402 + Icges(3,5) * t445;
t368 = Icges(3,1) * t401 - Icges(3,4) * t400 - Icges(3,5) * t444;
t367 = Icges(3,4) * t403 - Icges(3,2) * t402 + Icges(3,6) * t445;
t366 = Icges(3,4) * t401 - Icges(3,2) * t400 - Icges(3,6) * t444;
t363 = rSges(4,1) * t399 - rSges(4,2) * t398 - rSges(4,3) * t446;
t362 = Icges(4,1) * t399 - Icges(4,4) * t398 - Icges(4,5) * t446;
t361 = Icges(4,4) * t399 - Icges(4,2) * t398 - Icges(4,6) * t446;
t360 = Icges(4,5) * t399 - Icges(4,6) * t398 - Icges(4,3) * t446;
t359 = qJD(5) * t398 + t378;
t358 = t388 * t425 + t447;
t357 = -t388 * t421 + t402 * t425;
t356 = t386 * t425 + t448;
t355 = -t386 * t421 + t400 * t425;
t354 = t388 * t418 + t402 * t417;
t353 = -t388 * t417 + t402 * t418;
t352 = t386 * t418 + t400 * t417;
t351 = -t386 * t417 + t400 * t418;
t344 = rSges(4,1) * t388 - rSges(4,2) * t387 + rSges(4,3) * t402;
t343 = rSges(4,1) * t386 - rSges(4,2) * t385 + rSges(4,3) * t400;
t342 = Icges(4,1) * t388 - Icges(4,4) * t387 + Icges(4,5) * t402;
t341 = Icges(4,1) * t386 - Icges(4,4) * t385 + Icges(4,5) * t400;
t340 = Icges(4,4) * t388 - Icges(4,2) * t387 + Icges(4,6) * t402;
t339 = Icges(4,4) * t386 - Icges(4,2) * t385 + Icges(4,6) * t400;
t338 = Icges(4,5) * t388 - Icges(4,6) * t387 + Icges(4,3) * t402;
t337 = Icges(4,5) * t386 - Icges(4,6) * t385 + Icges(4,3) * t400;
t336 = rSges(5,1) * t384 + rSges(5,2) * t383 + rSges(5,3) * t398;
t335 = Icges(5,1) * t384 + Icges(5,4) * t383 + Icges(5,5) * t398;
t334 = Icges(5,4) * t384 + Icges(5,2) * t383 + Icges(5,6) * t398;
t333 = Icges(5,5) * t384 + Icges(5,6) * t383 + Icges(5,3) * t398;
t332 = -pkin(4) * t440 + pkin(10) * t398 + t451 * t399;
t331 = rSges(6,1) * t380 + rSges(6,2) * t379 + rSges(6,3) * t398;
t330 = qJD(5) * t385 + t350;
t329 = qJD(5) * t387 + t349;
t327 = Icges(6,1) * t380 + Icges(6,4) * t379 + Icges(6,5) * t398;
t326 = Icges(6,4) * t380 + Icges(6,2) * t379 + Icges(6,6) * t398;
t325 = Icges(6,5) * t380 + Icges(6,6) * t379 + Icges(6,3) * t398;
t324 = t372 * t413 - t394 * t412 + t406;
t323 = -t371 * t413 - t394 * t438 - t442;
t322 = (t371 * t424 + t372 * t427) * t441;
t321 = rSges(5,1) * t358 + rSges(5,2) * t357 + rSges(5,3) * t387;
t320 = rSges(5,1) * t356 + rSges(5,2) * t355 + rSges(5,3) * t385;
t319 = Icges(5,1) * t358 + Icges(5,4) * t357 + Icges(5,5) * t387;
t318 = Icges(5,1) * t356 + Icges(5,4) * t355 + Icges(5,5) * t385;
t317 = Icges(5,4) * t358 + Icges(5,2) * t357 + Icges(5,6) * t387;
t316 = Icges(5,4) * t356 + Icges(5,2) * t355 + Icges(5,6) * t385;
t315 = Icges(5,5) * t358 + Icges(5,6) * t357 + Icges(5,3) * t387;
t314 = Icges(5,5) * t356 + Icges(5,6) * t355 + Icges(5,3) * t385;
t313 = rSges(6,1) * t354 + rSges(6,2) * t353 + rSges(6,3) * t387;
t312 = rSges(6,1) * t352 + rSges(6,2) * t351 + rSges(6,3) * t385;
t311 = Icges(6,1) * t354 + Icges(6,4) * t353 + Icges(6,5) * t387;
t310 = Icges(6,1) * t352 + Icges(6,4) * t351 + Icges(6,5) * t385;
t309 = Icges(6,4) * t354 + Icges(6,2) * t353 + Icges(6,6) * t387;
t308 = Icges(6,4) * t352 + Icges(6,2) * t351 + Icges(6,6) * t385;
t307 = Icges(6,5) * t354 + Icges(6,6) * t353 + Icges(6,3) * t387;
t306 = Icges(6,5) * t352 + Icges(6,6) * t351 + Icges(6,3) * t385;
t305 = pkin(4) * t447 + pkin(10) * t387 + t451 * t388;
t304 = pkin(4) * t448 + pkin(10) * t385 + t451 * t386;
t303 = t344 * t405 - t363 * t389 + t433;
t302 = -t343 * t405 + t363 * t390 + t432;
t301 = t343 * t389 - t344 * t390 + t443;
t300 = t321 * t378 - t336 * t349 + t431;
t299 = -t320 * t378 + t336 * t350 + t430;
t298 = t320 * t349 - t321 * t350 + t434;
t297 = t305 * t378 + t313 * t359 - t329 * t331 - t332 * t349 + t431;
t296 = -t304 * t378 - t312 * t359 + t330 * t331 + t332 * t350 + t430;
t295 = t304 * t349 - t305 * t350 + t312 * t329 - t313 * t330 + t434;
t1 = m(3) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + ((t391 * t445 - t392 * t402 + t393 * t403) * t413 + (-(-t366 * t402 + t368 * t403) * t427 + (-t367 * t402 + t369 * t403 - t454) * t424) * t441) * t412 / 0.2e1 - ((-t391 * t444 - t392 * t400 + t393 * t401) * t413 + ((-t367 * t400 + t369 * t401) * t424 + (t366 * t400 - t368 * t401 + t454) * t427) * t441) * t438 / 0.2e1 + t413 * ((t449 * t365 + (t367 * t426 + t369 * t423) * t420) * t412 - (t449 * t364 + (t366 * t426 + t368 * t423) * t420) * t438 + (t449 * t391 + (t392 * t426 + t393 * t423) * t420) * t413) / 0.2e1 + m(4) * (t301 ^ 2 + t302 ^ 2 + t303 ^ 2) / 0.2e1 + t389 * ((t338 * t402 - t340 * t387 + t342 * t388) * t389 + (t337 * t402 - t339 * t387 + t341 * t388) * t390 + (t360 * t402 - t361 * t387 + t362 * t388) * t405) / 0.2e1 + t390 * ((t338 * t400 - t340 * t385 + t342 * t386) * t389 + (t337 * t400 - t339 * t385 + t341 * t386) * t390 + (t360 * t400 - t361 * t385 + t362 * t386) * t405) / 0.2e1 + t405 * ((-t338 * t446 - t340 * t398 + t342 * t399) * t389 + (-t337 * t446 - t339 * t398 + t341 * t399) * t390 + (-t360 * t446 - t361 * t398 + t362 * t399) * t405) / 0.2e1 + m(5) * (t298 ^ 2 + t299 ^ 2 + t300 ^ 2) / 0.2e1 + t349 * ((t315 * t387 + t317 * t357 + t319 * t358) * t349 + (t314 * t387 + t316 * t357 + t318 * t358) * t350 + (t333 * t387 + t334 * t357 + t335 * t358) * t378) / 0.2e1 + t350 * ((t315 * t385 + t317 * t355 + t319 * t356) * t349 + (t314 * t385 + t316 * t355 + t318 * t356) * t350 + (t333 * t385 + t334 * t355 + t335 * t356) * t378) / 0.2e1 + t378 * ((t315 * t398 + t317 * t383 + t319 * t384) * t349 + (t314 * t398 + t316 * t383 + t318 * t384) * t350 + (t333 * t398 + t334 * t383 + t335 * t384) * t378) / 0.2e1 + m(6) * (t295 ^ 2 + t296 ^ 2 + t297 ^ 2) / 0.2e1 + t329 * ((t307 * t387 + t309 * t353 + t311 * t354) * t329 + (t306 * t387 + t308 * t353 + t310 * t354) * t330 + (t325 * t387 + t326 * t353 + t327 * t354) * t359) / 0.2e1 + t330 * ((t307 * t385 + t309 * t351 + t311 * t352) * t329 + (t306 * t385 + t308 * t351 + t310 * t352) * t330 + (t325 * t385 + t326 * t351 + t327 * t352) * t359) / 0.2e1 + t359 * ((t307 * t398 + t309 * t379 + t311 * t380) * t329 + (t306 * t398 + t308 * t379 + t310 * t380) * t330 + (t325 * t398 + t326 * t379 + t327 * t380) * t359) / 0.2e1 + (m(2) * (t408 ^ 2 + t409 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
