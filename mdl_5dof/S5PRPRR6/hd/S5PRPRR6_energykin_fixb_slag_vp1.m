% Calculate kinetic energy for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:25
% EndTime: 2019-12-05 15:56:27
% DurationCPUTime: 2.06s
% Computational Cost: add. (1963->299), mult. (3755->474), div. (0->0), fcn. (4494->12), ass. (0->135)
t415 = sin(pkin(10));
t418 = cos(pkin(10));
t420 = cos(pkin(5));
t417 = sin(pkin(5));
t423 = sin(qJ(2));
t449 = t417 * t423;
t400 = -t415 * t449 + t418 * t420;
t452 = t415 * t420;
t401 = t418 * t449 + t452;
t425 = cos(qJ(2));
t448 = t417 * t425;
t358 = Icges(4,5) * t401 + Icges(4,6) * t400 - Icges(4,3) * t448;
t391 = Icges(3,6) * t420 + (Icges(3,4) * t423 + Icges(3,2) * t425) * t417;
t455 = t358 - t391;
t454 = qJD(2) ^ 2;
t453 = pkin(3) * t418;
t416 = sin(pkin(9));
t451 = t416 * t417;
t419 = cos(pkin(9));
t450 = t417 * t419;
t447 = t420 * t423;
t446 = t420 * t425;
t404 = t416 * t446 + t419 * t423;
t405 = -t416 * t447 + t419 * t425;
t376 = pkin(2) * t405 + qJ(3) * t404;
t402 = t416 * t423 - t419 * t446;
t413 = qJD(2) * t420;
t444 = qJD(3) * t402 + t376 * t413;
t443 = qJD(2) * t417;
t411 = t416 * t443;
t388 = qJD(4) * t404 + t411;
t442 = qJD(3) * t425;
t440 = pkin(10) + qJ(4);
t439 = t415 * t451;
t438 = t415 * t450;
t437 = t419 * t443;
t403 = t416 * t425 + t419 * t447;
t375 = pkin(2) * t403 + qJ(3) * t402;
t436 = t375 * t411 + t376 * t437 + qJD(1);
t434 = cos(t440);
t406 = (pkin(2) * t423 - qJ(3) * t425) * t417;
t433 = (-t401 * rSges(4,1) - t400 * rSges(4,2) + rSges(4,3) * t448 - t406) * t417;
t432 = (-pkin(3) * t452 - (-pkin(7) * t425 + t423 * t453) * t417 - t406) * t417;
t389 = qJD(4) * t402 - t437;
t407 = -qJD(4) * t448 + t413;
t429 = t417 * t434;
t334 = -pkin(3) * t438 + pkin(7) * t402 + t403 * t453;
t335 = pkin(3) * t439 + pkin(7) * t404 + t405 * t453;
t428 = t334 * t411 + t335 * t437 - t417 * t442 + t436;
t427 = qJD(2) * t416 * t432 + t335 * t413 + t444;
t399 = qJD(3) * t404;
t426 = t399 + ((-t334 - t375) * t420 + t419 * t432) * qJD(2);
t424 = cos(qJ(5));
t422 = sin(qJ(5));
t414 = sin(t440);
t395 = t420 * rSges(3,3) + (rSges(3,1) * t423 + rSges(3,2) * t425) * t417;
t394 = t420 * t414 + t423 * t429;
t393 = t414 * t449 - t420 * t434;
t392 = Icges(3,5) * t420 + (Icges(3,1) * t423 + Icges(3,4) * t425) * t417;
t390 = Icges(3,3) * t420 + (Icges(3,5) * t423 + Icges(3,6) * t425) * t417;
t387 = t405 * t418 + t439;
t386 = -t405 * t415 + t418 * t451;
t385 = t403 * t418 - t438;
t384 = -t403 * t415 - t418 * t450;
t383 = t394 * t424 - t422 * t448;
t382 = -t394 * t422 - t424 * t448;
t381 = t405 * t434 + t414 * t451;
t380 = t405 * t414 - t416 * t429;
t379 = t403 * t434 - t414 * t450;
t378 = t403 * t414 + t419 * t429;
t377 = qJD(5) * t393 + t407;
t373 = rSges(3,1) * t405 - rSges(3,2) * t404 + rSges(3,3) * t451;
t372 = rSges(3,1) * t403 - rSges(3,2) * t402 - rSges(3,3) * t450;
t368 = pkin(4) * t394 + pkin(8) * t393;
t366 = Icges(3,1) * t405 - Icges(3,4) * t404 + Icges(3,5) * t451;
t365 = Icges(3,1) * t403 - Icges(3,4) * t402 - Icges(3,5) * t450;
t364 = Icges(3,4) * t405 - Icges(3,2) * t404 + Icges(3,6) * t451;
t363 = Icges(3,4) * t403 - Icges(3,2) * t402 - Icges(3,6) * t450;
t362 = Icges(3,5) * t405 - Icges(3,6) * t404 + Icges(3,3) * t451;
t361 = Icges(3,5) * t403 - Icges(3,6) * t402 - Icges(3,3) * t450;
t360 = Icges(4,1) * t401 + Icges(4,4) * t400 - Icges(4,5) * t448;
t359 = Icges(4,4) * t401 + Icges(4,2) * t400 - Icges(4,6) * t448;
t357 = t394 * rSges(5,1) - t393 * rSges(5,2) - rSges(5,3) * t448;
t356 = Icges(5,1) * t394 - Icges(5,4) * t393 - Icges(5,5) * t448;
t355 = Icges(5,4) * t394 - Icges(5,2) * t393 - Icges(5,6) * t448;
t354 = Icges(5,5) * t394 - Icges(5,6) * t393 - Icges(5,3) * t448;
t353 = t381 * t424 + t404 * t422;
t352 = -t381 * t422 + t404 * t424;
t351 = t379 * t424 + t402 * t422;
t350 = -t379 * t422 + t402 * t424;
t349 = qJD(5) * t378 + t389;
t348 = qJD(5) * t380 + t388;
t347 = pkin(4) * t381 + pkin(8) * t380;
t346 = pkin(4) * t379 + pkin(8) * t378;
t345 = (-t372 * t420 - t395 * t450) * qJD(2);
t344 = (t373 * t420 - t395 * t451) * qJD(2);
t343 = rSges(4,1) * t387 + rSges(4,2) * t386 + rSges(4,3) * t404;
t342 = rSges(4,1) * t385 + rSges(4,2) * t384 + rSges(4,3) * t402;
t341 = Icges(4,1) * t387 + Icges(4,4) * t386 + Icges(4,5) * t404;
t340 = Icges(4,1) * t385 + Icges(4,4) * t384 + Icges(4,5) * t402;
t339 = Icges(4,4) * t387 + Icges(4,2) * t386 + Icges(4,6) * t404;
t338 = Icges(4,4) * t385 + Icges(4,2) * t384 + Icges(4,6) * t402;
t337 = Icges(4,5) * t387 + Icges(4,6) * t386 + Icges(4,3) * t404;
t336 = Icges(4,5) * t385 + Icges(4,6) * t384 + Icges(4,3) * t402;
t333 = rSges(5,1) * t381 - rSges(5,2) * t380 + rSges(5,3) * t404;
t332 = rSges(5,1) * t379 - rSges(5,2) * t378 + rSges(5,3) * t402;
t331 = Icges(5,1) * t381 - Icges(5,4) * t380 + Icges(5,5) * t404;
t330 = Icges(5,1) * t379 - Icges(5,4) * t378 + Icges(5,5) * t402;
t329 = Icges(5,4) * t381 - Icges(5,2) * t380 + Icges(5,6) * t404;
t328 = Icges(5,4) * t379 - Icges(5,2) * t378 + Icges(5,6) * t402;
t327 = Icges(5,5) * t381 - Icges(5,6) * t380 + Icges(5,3) * t404;
t326 = Icges(5,5) * t379 - Icges(5,6) * t378 + Icges(5,3) * t402;
t325 = rSges(6,1) * t383 + rSges(6,2) * t382 + rSges(6,3) * t393;
t324 = Icges(6,1) * t383 + Icges(6,4) * t382 + Icges(6,5) * t393;
t323 = Icges(6,4) * t383 + Icges(6,2) * t382 + Icges(6,6) * t393;
t322 = Icges(6,5) * t383 + Icges(6,6) * t382 + Icges(6,3) * t393;
t318 = qJD(1) + (t372 * t416 + t373 * t419) * t443;
t317 = rSges(6,1) * t353 + rSges(6,2) * t352 + rSges(6,3) * t380;
t316 = rSges(6,1) * t351 + rSges(6,2) * t350 + rSges(6,3) * t378;
t315 = Icges(6,1) * t353 + Icges(6,4) * t352 + Icges(6,5) * t380;
t314 = Icges(6,1) * t351 + Icges(6,4) * t350 + Icges(6,5) * t378;
t313 = Icges(6,4) * t353 + Icges(6,2) * t352 + Icges(6,6) * t380;
t312 = Icges(6,4) * t351 + Icges(6,2) * t350 + Icges(6,6) * t378;
t311 = Icges(6,5) * t353 + Icges(6,6) * t352 + Icges(6,3) * t380;
t310 = Icges(6,5) * t351 + Icges(6,6) * t350 + Icges(6,3) * t378;
t309 = t399 + ((-t342 - t375) * t420 + t419 * t433) * qJD(2);
t308 = (t343 * t420 + t416 * t433) * qJD(2) + t444;
t307 = (-t442 + (t342 * t416 + t343 * t419) * qJD(2)) * t417 + t436;
t306 = -t332 * t407 + t357 * t389 + t426;
t305 = t333 * t407 - t357 * t388 + t427;
t304 = t388 * t332 - t389 * t333 + t428;
t303 = -t316 * t377 + t325 * t349 - t346 * t407 + t368 * t389 + t426;
t302 = t317 * t377 - t325 * t348 + t347 * t407 - t368 * t388 + t427;
t301 = t348 * t316 - t349 * t317 + t388 * t346 - t389 * t347 + t428;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t318 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(4) * (t307 ^ 2 + t308 ^ 2 + t309 ^ 2) / 0.2e1 + m(5) * (t304 ^ 2 + t305 ^ 2 + t306 ^ 2) / 0.2e1 + t388 * ((t404 * t327 - t380 * t329 + t381 * t331) * t388 + (t326 * t404 - t328 * t380 + t330 * t381) * t389 + (t354 * t404 - t355 * t380 + t356 * t381) * t407) / 0.2e1 + t389 * ((t327 * t402 - t329 * t378 + t331 * t379) * t388 + (t402 * t326 - t378 * t328 + t379 * t330) * t389 + (t354 * t402 - t355 * t378 + t356 * t379) * t407) / 0.2e1 + t407 * ((-t327 * t448 - t393 * t329 + t394 * t331) * t388 + (-t326 * t448 - t393 * t328 + t394 * t330) * t389 + (-t354 * t448 - t393 * t355 + t394 * t356) * t407) / 0.2e1 + m(6) * (t301 ^ 2 + t302 ^ 2 + t303 ^ 2) / 0.2e1 + t348 * ((t380 * t311 + t352 * t313 + t353 * t315) * t348 + (t310 * t380 + t312 * t352 + t314 * t353) * t349 + (t322 * t380 + t323 * t352 + t324 * t353) * t377) / 0.2e1 + t349 * ((t311 * t378 + t313 * t350 + t315 * t351) * t348 + (t378 * t310 + t350 * t312 + t351 * t314) * t349 + (t322 * t378 + t323 * t350 + t324 * t351) * t377) / 0.2e1 + t377 * ((t311 * t393 + t313 * t382 + t315 * t383) * t348 + (t310 * t393 + t312 * t382 + t314 * t383) * t349 + (t393 * t322 + t382 * t323 + t383 * t324) * t377) / 0.2e1 - ((-t362 * t450 - t364 * t402 + t366 * t403) * t451 - (-t361 * t450 - t363 * t402 + t365 * t403) * t450 + ((t337 * t402 + t339 * t384 + t341 * t385) * t416 - (t336 * t402 + t338 * t384 + t340 * t385) * t419) * t417 + (t359 * t384 + t360 * t385 - t390 * t450 + t392 * t403 + t455 * t402) * t420) * t454 * t450 / 0.2e1 + ((((t364 * t425 + t366 * t423) * t416 - (t363 * t425 + t365 * t423) * t419) * t417 ^ 2 + (-t337 * t448 + t400 * t339 + t401 * t341) * t451 - (-t336 * t448 + t400 * t338 + t401 * t340) * t450 + ((-t361 * t419 + t362 * t416 + t391 * t425 + t392 * t423) * t417 - t358 * t448 + t400 * t359 + t401 * t360 + t420 * t390) * t420) * t420 + ((t362 * t451 - t364 * t404 + t366 * t405) * t451 - (t361 * t451 - t363 * t404 + t365 * t405) * t450 + ((t337 * t404 + t339 * t386 + t341 * t387) * t416 - (t336 * t404 + t338 * t386 + t340 * t387) * t419) * t417 + (t359 * t386 + t360 * t387 + t390 * t451 + t392 * t405 + t404 * t455) * t420) * t451) * t454 / 0.2e1;
T = t1;
