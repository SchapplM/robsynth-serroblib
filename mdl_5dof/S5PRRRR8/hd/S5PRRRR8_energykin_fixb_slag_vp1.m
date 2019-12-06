% Calculate kinetic energy for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:46
% EndTime: 2019-12-05 17:14:48
% DurationCPUTime: 2.02s
% Computational Cost: add. (2068->299), mult. (3965->479), div. (0->0), fcn. (4746->12), ass. (0->137)
t456 = qJD(2) ^ 2;
t429 = cos(qJ(3));
t455 = pkin(3) * t429;
t421 = sin(pkin(10));
t422 = sin(pkin(5));
t453 = t421 * t422;
t427 = sin(qJ(2));
t452 = t422 * t427;
t451 = t422 * t429;
t430 = cos(qJ(2));
t450 = t422 * t430;
t423 = cos(pkin(10));
t449 = t423 * t422;
t424 = cos(pkin(5));
t426 = sin(qJ(3));
t448 = t424 * t426;
t447 = t424 * t427;
t446 = t424 * t430;
t445 = qJ(3) + qJ(4);
t408 = t421 * t446 + t423 * t427;
t444 = qJD(2) * t422;
t417 = t421 * t444;
t393 = qJD(3) * t408 + t417;
t419 = qJD(2) * t424;
t442 = t426 * t453;
t441 = t426 * t449;
t361 = qJD(4) * t408 + t393;
t440 = t423 * t444;
t406 = t421 * t427 - t423 * t446;
t407 = t421 * t430 + t423 * t447;
t381 = t407 * pkin(2) + pkin(7) * t406;
t409 = -t421 * t447 + t423 * t430;
t382 = t409 * pkin(2) + t408 * pkin(7);
t439 = t381 * t417 + t382 * t440 + qJD(1);
t438 = cos(t445);
t437 = t422 * t438;
t394 = qJD(3) * t406 - t440;
t412 = (pkin(2) * t427 - pkin(7) * t430) * t422;
t436 = t382 * t419 - t412 * t417;
t362 = qJD(4) * t406 + t394;
t395 = t419 + (-qJD(3) - qJD(4)) * t450;
t336 = -pkin(3) * t441 + pkin(8) * t406 + t407 * t455;
t337 = pkin(3) * t442 + pkin(8) * t408 + t409 * t455;
t435 = t393 * t336 - t337 * t394 + t439;
t434 = (-t381 * t424 - t412 * t449) * qJD(2);
t378 = pkin(3) * t448 + (-pkin(8) * t430 + t427 * t455) * t422;
t413 = -qJD(3) * t450 + t419;
t433 = t413 * t337 - t378 * t393 + t436;
t432 = -t336 * t413 + t394 * t378 + t434;
t428 = cos(qJ(5));
t425 = sin(qJ(5));
t420 = sin(t445);
t411 = t427 * t451 + t448;
t410 = t424 * t429 - t426 * t452;
t401 = t424 * t420 + t427 * t437;
t400 = t420 * t452 - t424 * t438;
t399 = rSges(3,3) * t424 + (rSges(3,1) * t427 + rSges(3,2) * t430) * t422;
t398 = Icges(3,5) * t424 + (Icges(3,1) * t427 + Icges(3,4) * t430) * t422;
t397 = Icges(3,6) * t424 + (Icges(3,4) * t427 + Icges(3,2) * t430) * t422;
t396 = Icges(3,3) * t424 + (Icges(3,5) * t427 + Icges(3,6) * t430) * t422;
t392 = t409 * t429 + t442;
t391 = -t409 * t426 + t421 * t451;
t390 = t407 * t429 - t441;
t389 = -t407 * t426 - t429 * t449;
t388 = t401 * t428 - t425 * t450;
t387 = -t401 * t425 - t428 * t450;
t386 = t409 * t438 + t420 * t453;
t385 = t409 * t420 - t421 * t437;
t384 = t407 * t438 - t420 * t449;
t383 = t407 * t420 + t423 * t437;
t379 = qJD(5) * t400 + t395;
t377 = pkin(4) * t401 + pkin(9) * t400;
t376 = rSges(4,1) * t411 + rSges(4,2) * t410 - rSges(4,3) * t450;
t375 = Icges(4,1) * t411 + Icges(4,4) * t410 - Icges(4,5) * t450;
t374 = Icges(4,4) * t411 + Icges(4,2) * t410 - Icges(4,6) * t450;
t373 = Icges(4,5) * t411 + Icges(4,6) * t410 - Icges(4,3) * t450;
t370 = rSges(3,1) * t409 - rSges(3,2) * t408 + rSges(3,3) * t453;
t369 = rSges(3,1) * t407 - rSges(3,2) * t406 - rSges(3,3) * t449;
t368 = Icges(3,1) * t409 - Icges(3,4) * t408 + Icges(3,5) * t453;
t367 = Icges(3,1) * t407 - Icges(3,4) * t406 - Icges(3,5) * t449;
t366 = Icges(3,4) * t409 - Icges(3,2) * t408 + Icges(3,6) * t453;
t365 = Icges(3,4) * t407 - Icges(3,2) * t406 - Icges(3,6) * t449;
t364 = Icges(3,5) * t409 - Icges(3,6) * t408 + Icges(3,3) * t453;
t363 = Icges(3,5) * t407 - Icges(3,6) * t406 - Icges(3,3) * t449;
t360 = rSges(5,1) * t401 - rSges(5,2) * t400 - rSges(5,3) * t450;
t359 = Icges(5,1) * t401 - Icges(5,4) * t400 - Icges(5,5) * t450;
t358 = Icges(5,4) * t401 - Icges(5,2) * t400 - Icges(5,6) * t450;
t357 = Icges(5,5) * t401 - Icges(5,6) * t400 - Icges(5,3) * t450;
t356 = t386 * t428 + t408 * t425;
t355 = -t386 * t425 + t408 * t428;
t354 = t384 * t428 + t406 * t425;
t353 = -t384 * t425 + t406 * t428;
t352 = pkin(4) * t386 + pkin(9) * t385;
t351 = pkin(4) * t384 + pkin(9) * t383;
t349 = (-t369 * t424 - t399 * t449) * qJD(2);
t348 = (t370 * t424 - t399 * t453) * qJD(2);
t347 = qJD(5) * t383 + t362;
t346 = qJD(5) * t385 + t361;
t345 = rSges(4,1) * t392 + rSges(4,2) * t391 + rSges(4,3) * t408;
t344 = rSges(4,1) * t390 + rSges(4,2) * t389 + rSges(4,3) * t406;
t343 = Icges(4,1) * t392 + Icges(4,4) * t391 + Icges(4,5) * t408;
t342 = Icges(4,1) * t390 + Icges(4,4) * t389 + Icges(4,5) * t406;
t341 = Icges(4,4) * t392 + Icges(4,2) * t391 + Icges(4,6) * t408;
t340 = Icges(4,4) * t390 + Icges(4,2) * t389 + Icges(4,6) * t406;
t339 = Icges(4,5) * t392 + Icges(4,6) * t391 + Icges(4,3) * t408;
t338 = Icges(4,5) * t390 + Icges(4,6) * t389 + Icges(4,3) * t406;
t335 = rSges(5,1) * t386 - rSges(5,2) * t385 + rSges(5,3) * t408;
t334 = rSges(5,1) * t384 - rSges(5,2) * t383 + rSges(5,3) * t406;
t333 = Icges(5,1) * t386 - Icges(5,4) * t385 + Icges(5,5) * t408;
t332 = Icges(5,1) * t384 - Icges(5,4) * t383 + Icges(5,5) * t406;
t331 = Icges(5,4) * t386 - Icges(5,2) * t385 + Icges(5,6) * t408;
t330 = Icges(5,4) * t384 - Icges(5,2) * t383 + Icges(5,6) * t406;
t329 = Icges(5,5) * t386 - Icges(5,6) * t385 + Icges(5,3) * t408;
t328 = Icges(5,5) * t384 - Icges(5,6) * t383 + Icges(5,3) * t406;
t327 = rSges(6,1) * t388 + rSges(6,2) * t387 + rSges(6,3) * t400;
t326 = Icges(6,1) * t388 + Icges(6,4) * t387 + Icges(6,5) * t400;
t325 = Icges(6,4) * t388 + Icges(6,2) * t387 + Icges(6,6) * t400;
t324 = Icges(6,5) * t388 + Icges(6,6) * t387 + Icges(6,3) * t400;
t322 = qJD(1) + (t369 * t421 + t370 * t423) * t444;
t320 = rSges(6,1) * t356 + rSges(6,2) * t355 + rSges(6,3) * t385;
t319 = rSges(6,1) * t354 + rSges(6,2) * t353 + rSges(6,3) * t383;
t318 = Icges(6,1) * t356 + Icges(6,4) * t355 + Icges(6,5) * t385;
t317 = Icges(6,1) * t354 + Icges(6,4) * t353 + Icges(6,5) * t383;
t316 = Icges(6,4) * t356 + Icges(6,2) * t355 + Icges(6,6) * t385;
t315 = Icges(6,4) * t354 + Icges(6,2) * t353 + Icges(6,6) * t383;
t314 = Icges(6,5) * t356 + Icges(6,6) * t355 + Icges(6,3) * t385;
t313 = Icges(6,5) * t354 + Icges(6,6) * t353 + Icges(6,3) * t383;
t312 = -t344 * t413 + t376 * t394 + t434;
t311 = t345 * t413 - t376 * t393 + t436;
t310 = t344 * t393 - t345 * t394 + t439;
t309 = -t334 * t395 + t360 * t362 + t432;
t308 = t335 * t395 - t360 * t361 + t433;
t307 = t334 * t361 - t335 * t362 + t435;
t306 = -t319 * t379 + t327 * t347 - t351 * t395 + t362 * t377 + t432;
t305 = t320 * t379 - t327 * t346 + t352 * t395 - t361 * t377 + t433;
t304 = t319 * t346 - t320 * t347 + t351 * t361 - t352 * t362 + t435;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t322 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 - t456 * ((-t364 * t449 - t366 * t406 + t368 * t407) * t453 - (-t363 * t449 - t365 * t406 + t367 * t407) * t449 + (-t396 * t449 - t397 * t406 + t398 * t407) * t424) * t449 / 0.2e1 + m(4) * (t310 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + t393 * ((t408 * t339 + t391 * t341 + t392 * t343) * t393 + (t338 * t408 + t340 * t391 + t342 * t392) * t394 + (t373 * t408 + t374 * t391 + t375 * t392) * t413) / 0.2e1 + t394 * ((t339 * t406 + t341 * t389 + t343 * t390) * t393 + (t406 * t338 + t389 * t340 + t390 * t342) * t394 + (t373 * t406 + t374 * t389 + t375 * t390) * t413) / 0.2e1 + t413 * ((-t339 * t450 + t341 * t410 + t343 * t411) * t393 + (-t338 * t450 + t340 * t410 + t342 * t411) * t394 + (-t373 * t450 + t374 * t410 + t375 * t411) * t413) / 0.2e1 + m(5) * (t307 ^ 2 + t308 ^ 2 + t309 ^ 2) / 0.2e1 + t361 * ((t408 * t329 - t385 * t331 + t386 * t333) * t361 + (t328 * t408 - t330 * t385 + t332 * t386) * t362 + (t357 * t408 - t358 * t385 + t359 * t386) * t395) / 0.2e1 + t362 * ((t329 * t406 - t331 * t383 + t333 * t384) * t361 + (t406 * t328 - t383 * t330 + t384 * t332) * t362 + (t357 * t406 - t358 * t383 + t359 * t384) * t395) / 0.2e1 + t395 * ((-t329 * t450 - t331 * t400 + t333 * t401) * t361 + (-t328 * t450 - t330 * t400 + t332 * t401) * t362 + (-t357 * t450 - t400 * t358 + t401 * t359) * t395) / 0.2e1 + m(6) * (t304 ^ 2 + t305 ^ 2 + t306 ^ 2) / 0.2e1 + t346 * ((t385 * t314 + t355 * t316 + t356 * t318) * t346 + (t313 * t385 + t315 * t355 + t317 * t356) * t347 + (t324 * t385 + t325 * t355 + t326 * t356) * t379) / 0.2e1 + t347 * ((t314 * t383 + t316 * t353 + t318 * t354) * t346 + (t383 * t313 + t353 * t315 + t354 * t317) * t347 + (t324 * t383 + t325 * t353 + t326 * t354) * t379) / 0.2e1 + t379 * ((t314 * t400 + t316 * t387 + t318 * t388) * t346 + (t313 * t400 + t315 * t387 + t317 * t388) * t347 + (t400 * t324 + t387 * t325 + t388 * t326) * t379) / 0.2e1 + (((t364 * t453 - t366 * t408 + t368 * t409) * t453 - (t363 * t453 - t365 * t408 + t367 * t409) * t449 + (t396 * t453 - t397 * t408 + t398 * t409) * t424) * t453 + t424 * (t424 ^ 2 * t396 + (((t366 * t430 + t368 * t427) * t421 - (t365 * t430 + t367 * t427) * t423) * t422 + (-t363 * t423 + t364 * t421 + t397 * t430 + t398 * t427) * t424) * t422)) * t456 / 0.2e1;
T = t1;
