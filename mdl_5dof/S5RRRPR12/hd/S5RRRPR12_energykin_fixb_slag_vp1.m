% Calculate kinetic energy for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR12_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR12_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:37:02
% EndTime: 2019-12-31 21:37:04
% DurationCPUTime: 2.23s
% Computational Cost: add. (1938->301), mult. (4425->456), div. (0->0), fcn. (5386->12), ass. (0->134)
t462 = Icges(4,2) + Icges(5,3);
t426 = sin(qJ(2));
t427 = sin(qJ(1));
t428 = cos(qJ(2));
t429 = cos(qJ(1));
t451 = cos(pkin(5));
t437 = t429 * t451;
t401 = t426 * t427 - t428 * t437;
t402 = t426 * t437 + t427 * t428;
t422 = sin(pkin(5));
t446 = t422 * t429;
t365 = Icges(3,5) * t402 - Icges(3,6) * t401 - Icges(3,3) * t446;
t438 = t427 * t451;
t403 = t429 * t426 + t428 * t438;
t404 = -t426 * t438 + t429 * t428;
t448 = t422 * t427;
t366 = Icges(3,5) * t404 - Icges(3,6) * t403 + Icges(3,3) * t448;
t461 = (t365 * t429 - t366 * t427) * t422;
t425 = sin(qJ(3));
t453 = cos(qJ(3));
t387 = t402 * t453 - t425 * t446;
t421 = sin(pkin(10));
t423 = cos(pkin(10));
t357 = -t387 * t421 + t401 * t423;
t450 = t401 * t421;
t358 = t387 * t423 + t450;
t440 = t422 * t453;
t386 = t402 * t425 + t429 * t440;
t460 = -Icges(4,4) * t387 + Icges(5,5) * t358 - Icges(4,6) * t401 + Icges(5,6) * t357 + t462 * t386;
t389 = t404 * t453 + t425 * t448;
t359 = -t389 * t421 + t403 * t423;
t449 = t403 * t421;
t360 = t389 * t423 + t449;
t388 = t404 * t425 - t427 * t440;
t459 = -Icges(4,4) * t389 + Icges(5,5) * t360 - Icges(4,6) * t403 + Icges(5,6) * t359 + t462 * t388;
t400 = t425 * t451 + t426 * t440;
t447 = t422 * t428;
t384 = -t400 * t421 - t423 * t447;
t441 = t421 * t447;
t385 = t400 * t423 - t441;
t399 = t422 * t425 * t426 - t451 * t453;
t458 = -Icges(4,4) * t400 + Icges(5,5) * t385 + Icges(4,6) * t447 + Icges(5,6) * t384 + t462 * t399;
t452 = pkin(4) * t423;
t377 = pkin(2) * t402 + pkin(8) * t401;
t378 = pkin(2) * t404 + pkin(8) * t403;
t442 = qJD(2) * t422;
t413 = t427 * t442;
t439 = t429 * t442;
t444 = t377 * t413 + t378 * t439;
t390 = qJD(3) * t403 + t413;
t443 = qJD(1) * (pkin(1) * t427 - pkin(7) * t446);
t414 = qJD(2) * t451 + qJD(1);
t349 = pkin(3) * t387 + qJ(4) * t386;
t436 = qJD(4) * t399 + t390 * t349 + t444;
t391 = qJD(3) * t401 - t439;
t406 = -qJD(3) * t447 + t414;
t405 = (pkin(2) * t426 - pkin(8) * t428) * t422;
t407 = qJD(1) * (pkin(1) * t429 + pkin(7) * t448);
t434 = t414 * t378 - t405 * t413 + t407;
t350 = pkin(3) * t389 + qJ(4) * t388;
t433 = qJD(4) * t386 + t406 * t350 + t434;
t432 = -t377 * t414 - t405 * t439 - t443;
t376 = pkin(3) * t400 + qJ(4) * t399;
t431 = qJD(4) * t388 + t391 * t376 + t432;
t420 = pkin(10) + qJ(5);
t419 = cos(t420);
t418 = sin(t420);
t410 = rSges(2,1) * t429 - rSges(2,2) * t427;
t409 = rSges(2,1) * t427 + rSges(2,2) * t429;
t395 = t451 * rSges(3,3) + (rSges(3,1) * t426 + rSges(3,2) * t428) * t422;
t394 = Icges(3,5) * t451 + (Icges(3,1) * t426 + Icges(3,4) * t428) * t422;
t393 = Icges(3,6) * t451 + (Icges(3,4) * t426 + Icges(3,2) * t428) * t422;
t392 = Icges(3,3) * t451 + (Icges(3,5) * t426 + Icges(3,6) * t428) * t422;
t381 = qJD(5) * t399 + t406;
t380 = t400 * t419 - t418 * t447;
t379 = -t400 * t418 - t419 * t447;
t373 = rSges(3,1) * t404 - rSges(3,2) * t403 + rSges(3,3) * t448;
t372 = rSges(3,1) * t402 - rSges(3,2) * t401 - rSges(3,3) * t446;
t370 = Icges(3,1) * t404 - Icges(3,4) * t403 + Icges(3,5) * t448;
t369 = Icges(3,1) * t402 - Icges(3,4) * t401 - Icges(3,5) * t446;
t368 = Icges(3,4) * t404 - Icges(3,2) * t403 + Icges(3,6) * t448;
t367 = Icges(3,4) * t402 - Icges(3,2) * t401 - Icges(3,6) * t446;
t364 = rSges(4,1) * t400 - rSges(4,2) * t399 - rSges(4,3) * t447;
t363 = Icges(4,1) * t400 - Icges(4,4) * t399 - Icges(4,5) * t447;
t361 = Icges(4,5) * t400 - Icges(4,6) * t399 - Icges(4,3) * t447;
t356 = t389 * t419 + t403 * t418;
t355 = -t389 * t418 + t403 * t419;
t354 = t387 * t419 + t401 * t418;
t353 = -t387 * t418 + t401 * t419;
t352 = qJD(5) * t386 + t391;
t351 = qJD(5) * t388 + t390;
t346 = rSges(4,1) * t389 - rSges(4,2) * t388 + rSges(4,3) * t403;
t345 = rSges(4,1) * t387 - rSges(4,2) * t386 + rSges(4,3) * t401;
t344 = Icges(4,1) * t389 - Icges(4,4) * t388 + Icges(4,5) * t403;
t343 = Icges(4,1) * t387 - Icges(4,4) * t386 + Icges(4,5) * t401;
t340 = Icges(4,5) * t389 - Icges(4,6) * t388 + Icges(4,3) * t403;
t339 = Icges(4,5) * t387 - Icges(4,6) * t386 + Icges(4,3) * t401;
t338 = rSges(5,1) * t385 + rSges(5,2) * t384 + rSges(5,3) * t399;
t337 = Icges(5,1) * t385 + Icges(5,4) * t384 + Icges(5,5) * t399;
t336 = Icges(5,4) * t385 + Icges(5,2) * t384 + Icges(5,6) * t399;
t333 = -pkin(4) * t441 + pkin(9) * t399 + t400 * t452;
t332 = rSges(6,1) * t380 + rSges(6,2) * t379 + rSges(6,3) * t399;
t331 = Icges(6,1) * t380 + Icges(6,4) * t379 + Icges(6,5) * t399;
t330 = Icges(6,4) * t380 + Icges(6,2) * t379 + Icges(6,6) * t399;
t329 = Icges(6,5) * t380 + Icges(6,6) * t379 + Icges(6,3) * t399;
t328 = t373 * t414 - t395 * t413 + t407;
t327 = -t372 * t414 - t395 * t439 - t443;
t326 = (t372 * t427 + t373 * t429) * t442;
t325 = rSges(5,1) * t360 + rSges(5,2) * t359 + rSges(5,3) * t388;
t324 = rSges(5,1) * t358 + rSges(5,2) * t357 + rSges(5,3) * t386;
t323 = Icges(5,1) * t360 + Icges(5,4) * t359 + Icges(5,5) * t388;
t322 = Icges(5,1) * t358 + Icges(5,4) * t357 + Icges(5,5) * t386;
t321 = Icges(5,4) * t360 + Icges(5,2) * t359 + Icges(5,6) * t388;
t320 = Icges(5,4) * t358 + Icges(5,2) * t357 + Icges(5,6) * t386;
t317 = rSges(6,1) * t356 + rSges(6,2) * t355 + rSges(6,3) * t388;
t316 = rSges(6,1) * t354 + rSges(6,2) * t353 + rSges(6,3) * t386;
t315 = Icges(6,1) * t356 + Icges(6,4) * t355 + Icges(6,5) * t388;
t314 = Icges(6,1) * t354 + Icges(6,4) * t353 + Icges(6,5) * t386;
t313 = Icges(6,4) * t356 + Icges(6,2) * t355 + Icges(6,6) * t388;
t312 = Icges(6,4) * t354 + Icges(6,2) * t353 + Icges(6,6) * t386;
t311 = Icges(6,5) * t356 + Icges(6,6) * t355 + Icges(6,3) * t388;
t310 = Icges(6,5) * t354 + Icges(6,6) * t353 + Icges(6,3) * t386;
t309 = pkin(4) * t449 + pkin(9) * t388 + t389 * t452;
t308 = pkin(4) * t450 + pkin(9) * t386 + t387 * t452;
t307 = t346 * t406 - t364 * t390 + t434;
t306 = -t345 * t406 + t364 * t391 + t432;
t305 = t345 * t390 - t346 * t391 + t444;
t304 = t325 * t406 + (-t338 - t376) * t390 + t433;
t303 = t338 * t391 + (-t324 - t349) * t406 + t431;
t302 = t324 * t390 + (-t325 - t350) * t391 + t436;
t301 = t309 * t406 + t317 * t381 - t332 * t351 + (-t333 - t376) * t390 + t433;
t300 = -t316 * t381 + t332 * t352 + t333 * t391 + (-t308 - t349) * t406 + t431;
t299 = t308 * t390 + t316 * t351 - t317 * t352 + (-t309 - t350) * t391 + t436;
t1 = m(3) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + ((t392 * t448 - t393 * t403 + t394 * t404) * t414 + (-(-t367 * t403 + t369 * t404) * t429 + (-t403 * t368 + t404 * t370 - t461) * t427) * t442) * t413 / 0.2e1 - ((-t392 * t446 - t393 * t401 + t394 * t402) * t414 + ((-t368 * t401 + t370 * t402) * t427 + (t401 * t367 - t402 * t369 + t461) * t429) * t442) * t439 / 0.2e1 + t414 * ((t451 * t366 + (t368 * t428 + t370 * t426) * t422) * t413 - (t451 * t365 + (t367 * t428 + t369 * t426) * t422) * t439 + (t451 * t392 + (t393 * t428 + t394 * t426) * t422) * t414) / 0.2e1 + m(4) * (t305 ^ 2 + t306 ^ 2 + t307 ^ 2) / 0.2e1 + m(5) * (t302 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + m(6) * (t299 ^ 2 + t300 ^ 2 + t301 ^ 2) / 0.2e1 + t351 * ((t388 * t311 + t355 * t313 + t356 * t315) * t351 + (t310 * t388 + t312 * t355 + t314 * t356) * t352 + (t329 * t388 + t330 * t355 + t331 * t356) * t381) / 0.2e1 + t352 * ((t311 * t386 + t313 * t353 + t315 * t354) * t351 + (t386 * t310 + t353 * t312 + t354 * t314) * t352 + (t329 * t386 + t330 * t353 + t331 * t354) * t381) / 0.2e1 + t381 * ((t311 * t399 + t313 * t379 + t315 * t380) * t351 + (t310 * t399 + t312 * t379 + t314 * t380) * t352 + (t399 * t329 + t379 * t330 + t380 * t331) * t381) / 0.2e1 + ((t336 * t359 + t337 * t360 + t361 * t403 + t363 * t389 + t388 * t458) * t406 + (t320 * t359 + t322 * t360 + t339 * t403 + t343 * t389 + t388 * t460) * t391 + (t359 * t321 + t360 * t323 + t403 * t340 + t389 * t344 + t388 * t459) * t390) * t390 / 0.2e1 + ((t336 * t357 + t337 * t358 + t361 * t401 + t363 * t387 + t386 * t458) * t406 + (t357 * t320 + t358 * t322 + t401 * t339 + t387 * t343 + t386 * t460) * t391 + (t321 * t357 + t323 * t358 + t340 * t401 + t344 * t387 + t386 * t459) * t390) * t391 / 0.2e1 + ((t384 * t336 + t385 * t337 - t361 * t447 + t400 * t363 + t399 * t458) * t406 + (t320 * t384 + t322 * t385 - t339 * t447 + t343 * t400 + t399 * t460) * t391 + (t321 * t384 + t323 * t385 - t340 * t447 + t344 * t400 + t399 * t459) * t390) * t406 / 0.2e1 + (m(2) * (t409 ^ 2 + t410 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
