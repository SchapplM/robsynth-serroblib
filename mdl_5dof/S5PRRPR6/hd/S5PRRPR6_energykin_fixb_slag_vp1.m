% Calculate kinetic energy for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:04
% EndTime: 2019-12-05 16:30:06
% DurationCPUTime: 2.23s
% Computational Cost: add. (1886->293), mult. (4381->447), div. (0->0), fcn. (5358->12), ass. (0->131)
t453 = Icges(4,2) + Icges(5,3);
t415 = sin(pkin(9));
t418 = cos(pkin(9));
t423 = cos(qJ(2));
t419 = cos(pkin(5));
t422 = sin(qJ(2));
t437 = t419 * t422;
t398 = t415 * t423 + t418 * t437;
t421 = sin(qJ(3));
t416 = sin(pkin(5));
t438 = t418 * t416;
t445 = cos(qJ(3));
t385 = t398 * t445 - t421 * t438;
t436 = t419 * t423;
t397 = t415 * t422 - t418 * t436;
t414 = sin(pkin(10));
t417 = cos(pkin(10));
t355 = -t385 * t414 + t397 * t417;
t443 = t397 * t414;
t356 = t385 * t417 + t443;
t431 = t416 * t445;
t384 = t398 * t421 + t418 * t431;
t452 = -Icges(4,4) * t385 + Icges(5,5) * t356 - Icges(4,6) * t397 + Icges(5,6) * t355 + t453 * t384;
t400 = -t415 * t437 + t418 * t423;
t440 = t416 * t421;
t387 = t400 * t445 + t415 * t440;
t399 = t415 * t436 + t418 * t422;
t357 = -t387 * t414 + t399 * t417;
t442 = t399 * t414;
t358 = t387 * t417 + t442;
t386 = t400 * t421 - t415 * t431;
t451 = -Icges(4,4) * t387 + Icges(5,5) * t358 - Icges(4,6) * t399 + Icges(5,6) * t357 + t453 * t386;
t402 = t419 * t421 + t422 * t431;
t439 = t416 * t423;
t382 = -t402 * t414 - t417 * t439;
t432 = t414 * t439;
t383 = t402 * t417 - t432;
t401 = -t419 * t445 + t422 * t440;
t450 = -Icges(4,4) * t402 + Icges(5,5) * t383 + Icges(4,6) * t439 + Icges(5,6) * t382 + t453 * t401;
t449 = qJD(2) ^ 2;
t444 = pkin(4) * t417;
t441 = t415 * t416;
t434 = qJD(2) * t416;
t407 = t415 * t434;
t388 = qJD(3) * t399 + t407;
t410 = qJD(2) * t419;
t430 = t418 * t434;
t374 = pkin(2) * t398 + pkin(7) * t397;
t375 = pkin(2) * t400 + pkin(7) * t399;
t429 = t374 * t407 + t375 * t430 + qJD(1);
t389 = qJD(3) * t397 - t430;
t404 = -qJD(3) * t439 + t410;
t347 = pkin(3) * t385 + qJ(4) * t384;
t428 = qJD(4) * t401 + t388 * t347 + t429;
t403 = (pkin(2) * t422 - pkin(7) * t423) * t416;
t427 = t375 * t410 - t403 * t407;
t348 = pkin(3) * t387 + qJ(4) * t386;
t426 = qJD(4) * t384 + t404 * t348 + t427;
t425 = (-t374 * t419 - t403 * t438) * qJD(2);
t376 = t402 * pkin(3) + t401 * qJ(4);
t424 = qJD(4) * t386 + t389 * t376 + t425;
t413 = pkin(10) + qJ(5);
t412 = cos(t413);
t411 = sin(t413);
t393 = t419 * rSges(3,3) + (rSges(3,1) * t422 + rSges(3,2) * t423) * t416;
t392 = Icges(3,5) * t419 + (Icges(3,1) * t422 + Icges(3,4) * t423) * t416;
t391 = Icges(3,6) * t419 + (Icges(3,4) * t422 + Icges(3,2) * t423) * t416;
t390 = Icges(3,3) * t419 + (Icges(3,5) * t422 + Icges(3,6) * t423) * t416;
t381 = qJD(5) * t401 + t404;
t378 = t402 * t412 - t411 * t439;
t377 = -t402 * t411 - t412 * t439;
t372 = t402 * rSges(4,1) - t401 * rSges(4,2) - rSges(4,3) * t439;
t371 = Icges(4,1) * t402 - Icges(4,4) * t401 - Icges(4,5) * t439;
t369 = Icges(4,5) * t402 - Icges(4,6) * t401 - Icges(4,3) * t439;
t366 = rSges(3,1) * t400 - rSges(3,2) * t399 + rSges(3,3) * t441;
t365 = rSges(3,1) * t398 - rSges(3,2) * t397 - rSges(3,3) * t438;
t364 = Icges(3,1) * t400 - Icges(3,4) * t399 + Icges(3,5) * t441;
t363 = Icges(3,1) * t398 - Icges(3,4) * t397 - Icges(3,5) * t438;
t362 = Icges(3,4) * t400 - Icges(3,2) * t399 + Icges(3,6) * t441;
t361 = Icges(3,4) * t398 - Icges(3,2) * t397 - Icges(3,6) * t438;
t360 = Icges(3,5) * t400 - Icges(3,6) * t399 + Icges(3,3) * t441;
t359 = Icges(3,5) * t398 - Icges(3,6) * t397 - Icges(3,3) * t438;
t354 = t387 * t412 + t399 * t411;
t353 = -t387 * t411 + t399 * t412;
t352 = t385 * t412 + t397 * t411;
t351 = -t385 * t411 + t397 * t412;
t350 = qJD(5) * t384 + t389;
t349 = qJD(5) * t386 + t388;
t344 = (-t365 * t419 - t393 * t438) * qJD(2);
t343 = (t366 * t419 - t393 * t441) * qJD(2);
t342 = rSges(5,1) * t383 + rSges(5,2) * t382 + rSges(5,3) * t401;
t341 = rSges(4,1) * t387 - rSges(4,2) * t386 + rSges(4,3) * t399;
t340 = rSges(4,1) * t385 - rSges(4,2) * t384 + rSges(4,3) * t397;
t339 = Icges(5,1) * t383 + Icges(5,4) * t382 + Icges(5,5) * t401;
t338 = Icges(5,4) * t383 + Icges(5,2) * t382 + Icges(5,6) * t401;
t336 = Icges(4,1) * t387 - Icges(4,4) * t386 + Icges(4,5) * t399;
t335 = Icges(4,1) * t385 - Icges(4,4) * t384 + Icges(4,5) * t397;
t332 = Icges(4,5) * t387 - Icges(4,6) * t386 + Icges(4,3) * t399;
t331 = Icges(4,5) * t385 - Icges(4,6) * t384 + Icges(4,3) * t397;
t330 = -pkin(4) * t432 + pkin(8) * t401 + t402 * t444;
t329 = rSges(6,1) * t378 + rSges(6,2) * t377 + rSges(6,3) * t401;
t328 = Icges(6,1) * t378 + Icges(6,4) * t377 + Icges(6,5) * t401;
t327 = Icges(6,4) * t378 + Icges(6,2) * t377 + Icges(6,6) * t401;
t326 = Icges(6,5) * t378 + Icges(6,6) * t377 + Icges(6,3) * t401;
t324 = qJD(1) + (t365 * t415 + t366 * t418) * t434;
t323 = rSges(5,1) * t358 + rSges(5,2) * t357 + rSges(5,3) * t386;
t322 = rSges(5,1) * t356 + rSges(5,2) * t355 + rSges(5,3) * t384;
t321 = Icges(5,1) * t358 + Icges(5,4) * t357 + Icges(5,5) * t386;
t320 = Icges(5,1) * t356 + Icges(5,4) * t355 + Icges(5,5) * t384;
t319 = Icges(5,4) * t358 + Icges(5,2) * t357 + Icges(5,6) * t386;
t318 = Icges(5,4) * t356 + Icges(5,2) * t355 + Icges(5,6) * t384;
t315 = rSges(6,1) * t354 + rSges(6,2) * t353 + rSges(6,3) * t386;
t314 = rSges(6,1) * t352 + rSges(6,2) * t351 + rSges(6,3) * t384;
t313 = Icges(6,1) * t354 + Icges(6,4) * t353 + Icges(6,5) * t386;
t312 = Icges(6,1) * t352 + Icges(6,4) * t351 + Icges(6,5) * t384;
t311 = Icges(6,4) * t354 + Icges(6,2) * t353 + Icges(6,6) * t386;
t310 = Icges(6,4) * t352 + Icges(6,2) * t351 + Icges(6,6) * t384;
t309 = Icges(6,5) * t354 + Icges(6,6) * t353 + Icges(6,3) * t386;
t308 = Icges(6,5) * t352 + Icges(6,6) * t351 + Icges(6,3) * t384;
t307 = pkin(4) * t442 + pkin(8) * t386 + t387 * t444;
t306 = pkin(4) * t443 + pkin(8) * t384 + t385 * t444;
t305 = -t340 * t404 + t372 * t389 + t425;
t304 = t341 * t404 - t372 * t388 + t427;
t303 = t340 * t388 - t341 * t389 + t429;
t302 = t342 * t389 + (-t322 - t347) * t404 + t424;
t301 = t323 * t404 + (-t342 - t376) * t388 + t426;
t300 = t322 * t388 + (-t323 - t348) * t389 + t428;
t299 = -t314 * t381 + t329 * t350 + t330 * t389 + (-t306 - t347) * t404 + t424;
t298 = t307 * t404 + t315 * t381 - t329 * t349 + (-t330 - t376) * t388 + t426;
t297 = t306 * t388 + t314 * t349 - t315 * t350 + (-t307 - t348) * t389 + t428;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t324 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 - t449 * ((-t360 * t438 - t362 * t397 + t364 * t398) * t441 - (-t359 * t438 - t361 * t397 + t363 * t398) * t438 + (-t390 * t438 - t391 * t397 + t392 * t398) * t419) * t438 / 0.2e1 + m(4) * (t303 ^ 2 + t304 ^ 2 + t305 ^ 2) / 0.2e1 + m(5) * (t300 ^ 2 + t301 ^ 2 + t302 ^ 2) / 0.2e1 + m(6) * (t297 ^ 2 + t298 ^ 2 + t299 ^ 2) / 0.2e1 + t349 * ((t386 * t309 + t353 * t311 + t354 * t313) * t349 + (t308 * t386 + t310 * t353 + t312 * t354) * t350 + (t326 * t386 + t327 * t353 + t328 * t354) * t381) / 0.2e1 + t350 * ((t309 * t384 + t311 * t351 + t313 * t352) * t349 + (t384 * t308 + t351 * t310 + t352 * t312) * t350 + (t326 * t384 + t327 * t351 + t328 * t352) * t381) / 0.2e1 + t381 * ((t309 * t401 + t311 * t377 + t313 * t378) * t349 + (t308 * t401 + t310 * t377 + t312 * t378) * t350 + (t401 * t326 + t377 * t327 + t378 * t328) * t381) / 0.2e1 + ((t338 * t357 + t339 * t358 + t369 * t399 + t371 * t387 + t450 * t386) * t404 + (t318 * t357 + t320 * t358 + t331 * t399 + t335 * t387 + t452 * t386) * t389 + (t357 * t319 + t358 * t321 + t399 * t332 + t387 * t336 + t451 * t386) * t388) * t388 / 0.2e1 + ((t338 * t355 + t339 * t356 + t369 * t397 + t371 * t385 + t450 * t384) * t404 + (t355 * t318 + t356 * t320 + t397 * t331 + t385 * t335 + t452 * t384) * t389 + (t319 * t355 + t321 * t356 + t332 * t397 + t336 * t385 + t451 * t384) * t388) * t389 / 0.2e1 + ((t382 * t338 + t383 * t339 - t369 * t439 + t402 * t371 + t450 * t401) * t404 + (t318 * t382 + t320 * t383 - t331 * t439 + t402 * t335 + t452 * t401) * t389 + (t319 * t382 + t321 * t383 - t332 * t439 + t402 * t336 + t451 * t401) * t388) * t404 / 0.2e1 + (((t360 * t441 - t362 * t399 + t364 * t400) * t441 - (t359 * t441 - t361 * t399 + t363 * t400) * t438 + (t390 * t441 - t391 * t399 + t392 * t400) * t419) * t441 + t419 * (t419 ^ 2 * t390 + (((t362 * t423 + t364 * t422) * t415 - (t361 * t423 + t363 * t422) * t418) * t416 + (-t359 * t418 + t360 * t415 + t391 * t423 + t392 * t422) * t419) * t416)) * t449 / 0.2e1;
T = t1;
