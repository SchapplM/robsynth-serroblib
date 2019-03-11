% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:06:43
% EndTime: 2019-03-08 21:06:50
% DurationCPUTime: 6.13s
% Computational Cost: add. (3100->483), mult. (7299->635), div. (0->0), fcn. (5749->14), ass. (0->221)
t426 = sin(pkin(11));
t433 = sin(qJ(3));
t511 = qJD(2) * t433;
t429 = cos(pkin(11));
t436 = cos(qJ(3));
t528 = t429 * t436;
t381 = -qJD(2) * t528 + t426 * t511;
t435 = cos(qJ(6));
t432 = sin(qJ(6));
t509 = qJD(3) * t432;
t356 = -t435 * t381 + t509;
t390 = t426 * t436 + t429 * t433;
t384 = t390 * qJD(2);
t569 = qJD(6) + t384;
t570 = t356 * t569;
t358 = qJD(3) * t435 + t381 * t432;
t477 = t569 * t358;
t475 = MDP(23) * t569;
t549 = qJ(4) + pkin(8);
t430 = cos(pkin(6));
t428 = sin(pkin(6));
t434 = sin(qJ(2));
t532 = t428 * t434;
t387 = t430 * t436 - t433 * t532;
t568 = MDP(12) + MDP(14);
t567 = t432 * t569;
t483 = qJD(3) * t549;
t374 = qJD(4) * t436 - t433 * t483;
t375 = -qJD(4) * t433 - t436 * t483;
t437 = cos(qJ(2));
t530 = t428 * t437;
t490 = qJD(1) * t530;
t518 = t374 * t426 - t429 * t375 - t390 * t490;
t536 = t426 * t433;
t563 = t528 - t536;
t517 = t374 * t429 + t375 * t426 - t563 * t490;
t547 = cos(pkin(10));
t480 = t547 * t437;
t427 = sin(pkin(10));
t534 = t427 * t434;
t377 = -t430 * t480 + t534;
t481 = t547 * t434;
t533 = t427 * t437;
t379 = t430 * t533 + t481;
t472 = g(1) * t379 + g(2) * t377;
t455 = g(3) * t530 - t472;
t439 = qJD(2) ^ 2;
t497 = qJDD(2) * t437;
t565 = -t434 * t439 + t497;
t500 = qJD(2) * qJD(3);
t487 = t433 * t500;
t401 = t426 * t487;
t486 = t436 * t500;
t338 = qJDD(2) * t390 + t429 * t486 - t401;
t564 = -qJ(5) * t338 - qJD(5) * t384;
t422 = qJ(3) + pkin(11);
t419 = sin(t422);
t420 = cos(t422);
t562 = pkin(4) * t420 + qJ(5) * t419;
t512 = qJD(1) * t434;
t491 = t428 * t512;
t476 = qJD(2) * t549 + t491;
t513 = qJD(1) * t430;
t351 = t433 * t513 + t436 * t476;
t343 = t426 * t351;
t350 = -t433 * t476 + t436 * t513;
t316 = t350 * t429 - t343;
t561 = -qJD(5) + t316;
t560 = pkin(3) * t487 + qJDD(4);
t347 = qJD(3) * pkin(3) + t350;
t529 = t429 * t351;
t312 = t426 * t347 + t529;
t307 = -qJD(3) * qJ(5) - t312;
t551 = pkin(5) * t381;
t297 = -t307 - t551;
t314 = t350 * t426 + t529;
t333 = qJDD(6) + t338;
t554 = pkin(3) * t429;
t415 = -pkin(4) - t554;
t411 = -pkin(9) + t415;
t559 = t411 * t333 + (t297 - t314 + t551) * t569;
t438 = qJD(3) ^ 2;
t501 = qJD(1) * qJD(2);
t488 = t434 * t501;
t402 = t428 * t488;
t484 = qJDD(1) * t530;
t466 = t402 - t484;
t558 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t438 + t428 * (-g(3) * t437 + t488) - t466 + t472;
t376 = t384 ^ 2;
t556 = pkin(4) + pkin(9);
t555 = pkin(3) * t426;
t383 = t390 * qJD(3);
t498 = qJDD(2) * t436;
t404 = t429 * t498;
t499 = qJDD(2) * t433;
t337 = qJD(2) * t383 + t426 * t499 - t404;
t553 = pkin(4) * t337;
t550 = pkin(5) * t384;
t548 = qJD(2) * pkin(2);
t543 = qJDD(3) * pkin(4);
t505 = qJD(6) * t435;
t492 = t435 * qJDD(3) + t432 * t337 + t381 * t505;
t506 = qJD(6) * t432;
t309 = -qJD(3) * t506 + t492;
t542 = t309 * t435;
t416 = pkin(3) * t436 + pkin(2);
t461 = -qJ(5) * t390 - t416;
t322 = -t556 * t563 + t461;
t541 = t322 * t333;
t540 = t333 * t432;
t539 = t356 * t381;
t538 = t358 * t381;
t537 = t563 * t432;
t535 = t427 * t428;
t531 = t428 * t436;
t526 = t549 * t434;
t525 = t432 * t437;
t330 = t435 * t333;
t523 = t435 * t437;
t522 = t436 * t439;
t521 = qJDD(1) - g(3);
t496 = t430 * qJDD(1);
t405 = t436 * t496;
t362 = qJDD(2) * pkin(8) + (qJDD(1) * t434 + t437 * t501) * t428;
t448 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t513 + t362;
t462 = t476 * qJD(3);
t303 = qJDD(3) * pkin(3) - t433 * t448 - t436 * t462 + t405;
t304 = (-t462 + t496) * t433 + t448 * t436;
t293 = t426 * t303 + t429 * t304;
t508 = qJD(3) * t433;
t386 = qJD(3) * t528 - t426 * t508;
t520 = pkin(5) * t386 + t518;
t519 = -pkin(5) * t383 + t517;
t378 = t430 * t481 + t533;
t516 = -t377 * t416 + t378 * t549;
t380 = -t430 * t534 + t480;
t515 = -t379 * t416 + t380 * t549;
t424 = t433 ^ 2;
t514 = -t436 ^ 2 + t424;
t510 = qJD(2) * t434;
t504 = qJD(6) * t437;
t503 = t550 - t561;
t495 = g(3) * t532;
t418 = pkin(3) * t508;
t417 = pkin(3) * t511;
t493 = qJDD(3) * qJ(5) + t293;
t489 = qJD(2) * t530;
t485 = t437 * t500;
t482 = t428 * t547;
t479 = qJ(5) * t381 + t417;
t292 = t303 * t429 - t426 * t304;
t311 = t347 * t429 - t343;
t397 = t549 * t433;
t398 = t549 * t436;
t354 = t429 * t397 + t398 * t426;
t478 = t435 * t569;
t474 = qJDD(3) * t432 - t435 * t337;
t473 = (-t380 * t433 + t427 * t531) * pkin(3);
t471 = g(1) * t380 + g(2) * t378;
t457 = -qJ(5) * t386 - qJD(5) * t390 + t418;
t470 = -t383 * t556 - t457 + t491;
t320 = pkin(4) * t383 + t457;
t469 = -t320 + t491;
t468 = qJD(5) - t311;
t465 = qJDD(5) - t292;
t296 = -qJD(3) * t556 + t468 + t550;
t371 = -qJD(2) * t416 + qJD(4) - t490;
t445 = -qJ(5) * t384 + t371;
t308 = t381 * t556 + t445;
t289 = t296 * t435 - t308 * t432;
t290 = t296 * t432 + t308 * t435;
t355 = -t397 * t426 + t398 * t429;
t463 = t387 * pkin(3);
t288 = -qJD(3) * qJD(5) - t493;
t460 = -g(1) * t427 + g(2) * t547;
t388 = t430 * t433 + t434 * t531;
t458 = t383 * t432 - t505 * t563;
t340 = t378 * t420 - t419 * t482;
t342 = t380 * t420 + t419 * t535;
t366 = t419 * t430 + t420 * t532;
t456 = -g(1) * t342 - g(2) * t340 - g(3) * t366;
t395 = -t490 - t548;
t454 = -qJD(2) * t395 - t362 + t471;
t453 = (-t378 * t433 - t436 * t482) * pkin(3);
t452 = t455 * t419;
t451 = t484 - t455;
t447 = -pkin(8) * qJDD(3) + (t395 + t490 - t548) * qJD(3);
t287 = -pkin(5) * t337 - t288;
t326 = pkin(5) * t390 + t354;
t446 = -t287 * t563 + t297 * t383 - t326 * t333 + t471;
t336 = -qJDD(2) * t416 + t466 + t560;
t444 = t402 - t451 + t560;
t443 = t287 + (-qJD(6) * t411 + t384 * t556 + t479) * t569 + t456;
t321 = pkin(4) * t381 + t445;
t339 = t378 * t419 + t420 * t482;
t341 = t380 * t419 - t420 * t535;
t365 = t419 * t532 - t420 * t430;
t442 = -g(1) * t341 - g(2) * t339 - g(3) * t365 + t321 * t384 + t465;
t441 = t336 + t564;
t440 = -t337 * t355 + t338 * t354 - t381 * t517 + t384 * t518 - t471 - t495;
t413 = qJ(5) + t555;
t403 = t428 * t525;
t393 = t416 * t530;
t372 = qJD(3) * t381;
t349 = -qJD(3) * t388 - t433 * t489;
t348 = qJD(3) * t387 + t436 * t489;
t334 = -pkin(4) * t563 + t461;
t329 = t387 * t426 + t388 * t429;
t328 = -t429 * t387 + t388 * t426;
t327 = pkin(5) * t563 + t355;
t323 = pkin(4) * t384 + t479;
t315 = t348 * t429 + t349 * t426;
t313 = t348 * t426 - t429 * t349;
t310 = qJD(6) * t358 + t474;
t306 = -qJD(3) * pkin(4) + t468;
t295 = t441 + t553;
t294 = t337 * t556 + t441;
t291 = t465 - t543;
t286 = pkin(5) * t338 - qJDD(3) * t556 + t465;
t285 = t435 * t286;
t1 = [t521 * MDP(1) + (qJD(3) * t349 + qJDD(3) * t387) * MDP(10) + (-qJD(3) * t348 - qJDD(3) * t388) * MDP(11) + (-t292 * t328 + t293 * t329 - t311 * t313 + t312 * t315 - g(3)) * MDP(13) + (qJD(3) * t313 + qJDD(3) * t328) * MDP(15) + (qJD(3) * t315 + qJDD(3) * t329) * MDP(16) + (-t288 * t329 + t291 * t328 + t306 * t313 - t307 * t315 - g(3)) * MDP(17) + ((t313 * t435 - t328 * t506) * t569 + (t328 * t435 + t403) * t333 + t315 * t356 + t329 * t310) * MDP(23) + (-(t313 * t432 + t328 * t505) * t569 - t328 * t540 + t315 * t358 + t329 * t309) * MDP(24) + (t565 * MDP(3) + (-qJDD(2) * t434 - t437 * t439) * MDP(4) + (-t433 * t485 - t434 * t522 + t436 * t497) * MDP(10) + (-t433 * t565 - t436 * t485) * MDP(11) + (-t336 * t437 + t371 * t510) * MDP(13) + (t337 * t437 - t381 * t510) * MDP(15) + (t338 * t437 - t384 * t510) * MDP(16) + (-t295 * t437 + t321 * t510) * MDP(17) + (-t432 * t510 + t435 * t504) * t475 + (-(t432 * t504 + t435 * t510) * t569 + t333 * t523) * MDP(24)) * t428 + t568 * (t313 * t384 - t315 * t381 + t328 * t338 - t329 * t337); qJDD(2) * MDP(2) + t451 * MDP(3) + (-t521 * t532 + t471) * MDP(4) + (qJDD(2) * t424 + 0.2e1 * t433 * t486) * MDP(5) + 0.2e1 * (t433 * t498 - t500 * t514) * MDP(6) + (qJDD(3) * t433 + t436 * t438) * MDP(7) + (qJDD(3) * t436 - t433 * t438) * MDP(8) + (t447 * t433 + t436 * t558) * MDP(10) + (-t433 * t558 + t447 * t436) * MDP(11) + (-t292 * t390 + t293 * t563 - t311 * t386 - t312 * t383 + t440) * MDP(12) + (t293 * t355 - t292 * t354 - t336 * t416 - g(1) * t515 - g(2) * t516 - g(3) * (t428 * t526 + t393) + (-t491 + t418) * t371 + t517 * t312 - t518 * t311) * MDP(13) + (-t288 * t563 + t291 * t390 + t306 * t386 + t307 * t383 + t440) * MDP(14) + (qJD(3) * t518 + qJDD(3) * t354 + t295 * t563 - t321 * t383 - t334 * t337 + t381 * t469 + t420 * t455) * MDP(15) + (qJD(3) * t517 + qJDD(3) * t355 - t295 * t390 - t321 * t386 - t334 * t338 + t384 * t469 - t452) * MDP(16) + (t295 * t334 + t321 * t320 - t288 * t355 + t291 * t354 - g(1) * (-t379 * t562 + t515) - g(2) * (-t377 * t562 + t516) - g(3) * t393 - t517 * t307 + t518 * t306 + (-t321 * t512 - g(3) * (t437 * t562 + t526)) * t428) * MDP(17) + (-t309 * t537 + t358 * t458) * MDP(18) + ((-t356 * t432 + t358 * t435) * t383 - (t542 - t310 * t432 + (-t356 * t435 - t358 * t432) * qJD(6)) * t563) * MDP(19) + (t309 * t390 - t333 * t537 + t358 * t386 + t458 * t569) * MDP(20) + (-t563 * t330 - t310 * t390 - t356 * t386 + (t383 * t435 + t506 * t563) * t569) * MDP(21) + (t333 * t390 + t386 * t569) * MDP(22) + (t285 * t390 + t289 * t386 + t327 * t310 + (-t294 * t390 + t419 * t472 - t541) * t432 - t446 * t435 - g(3) * (t419 * t525 + t434 * t435) * t428 + (t470 * t432 + t520 * t435) * t569 + t519 * t356 + ((-t322 * t435 - t326 * t432) * t569 - t290 * t390 - t297 * t537) * qJD(6)) * MDP(23) + (-t290 * t386 + t327 * t309 + t519 * t358 + (-t541 - (qJD(6) * t296 + t294) * t390 - t297 * qJD(6) * t563 + (-qJD(6) * t326 + t470) * t569 - t452) * t435 + (-(-qJD(6) * t308 + t286) * t390 + t495 + (qJD(6) * t322 - t520) * t569 + t446) * t432) * MDP(24); -t433 * MDP(5) * t522 + t514 * MDP(6) * t439 + MDP(7) * t499 + MDP(8) * t498 + qJDD(3) * MDP(9) + (-g(3) * t387 + t433 * t454 + t460 * t531 + t405) * MDP(10) + (g(3) * t388 + (-t428 * t460 - t496) * t433 + t454 * t436) * MDP(11) + ((t312 - t314) * t384 + (-t311 + t316) * t381 + (-t337 * t426 - t338 * t429) * pkin(3)) * MDP(12) + (-g(1) * t473 - g(2) * t453 - g(3) * t463 + t292 * t554 + t293 * t555 + t311 * t314 - t312 * t316 - t371 * t417) * MDP(13) + (-t337 * t413 + t338 * t415 + (-t307 - t314) * t384 + (t306 + t561) * t381) * MDP(14) + (-qJD(3) * t314 + t323 * t381 + (-pkin(4) + t415) * qJDD(3) + t442) * MDP(15) + (qJDD(3) * t413 - t321 * t381 + t323 * t384 + (0.2e1 * qJD(5) - t316) * qJD(3) + t456 + t493) * MDP(16) + (-t288 * t413 + t291 * t415 - t321 * t323 - t306 * t314 - g(1) * (-pkin(4) * t341 + qJ(5) * t342 + t473) - g(2) * (-t339 * pkin(4) + t340 * qJ(5) + t453) - g(3) * (-pkin(4) * t365 + qJ(5) * t366 + t463) + t561 * t307) * MDP(17) + (-t432 * t477 + t542) * MDP(18) + ((-t310 - t477) * t435 + (-t309 + t570) * t432) * MDP(19) + (-t567 * t569 + t330 + t538) * MDP(20) + (-t478 * t569 - t539 - t540) * MDP(21) + t569 * t381 * MDP(22) + (t289 * t381 + t413 * t310 + t503 * t356 + t443 * t432 + t435 * t559) * MDP(23) + (-t290 * t381 + t413 * t309 + t503 * t358 - t432 * t559 + t443 * t435) * MDP(24); (t311 * t384 + t312 * t381 + t444) * MDP(13) + t404 * MDP(15) + (t372 + t401) * MDP(16) + (-t306 * t384 - t307 * t381 + t444 + t553 + t564) * MDP(17) + (t539 - t540) * MDP(23) + (-t330 + t538) * MDP(24) + (MDP(24) * t567 - t435 * t475) * t569 + (-t384 * MDP(15) + (-MDP(15) * t390 - MDP(16) * t528) * qJD(2)) * qJD(3) + (-MDP(15) * t536 - MDP(16) * t390 + (-MDP(13) - MDP(17)) * t416) * qJDD(2) + t568 * (-t381 ^ 2 - t376); (t372 + t338) * MDP(14) + (-t381 * t384 + qJDD(3)) * MDP(15) + (-t376 - t438) * MDP(16) + (qJD(3) * t307 + t442 - t543) * MDP(17) + (-qJD(3) * t356 + t330) * MDP(23) + (-qJD(3) * t358 - t540) * MDP(24) + (-MDP(24) * t478 - t432 * t475) * t569; t358 * t356 * MDP(18) + (-t356 ^ 2 + t358 ^ 2) * MDP(19) + (t492 + t570) * MDP(20) + (-t474 + t477) * MDP(21) + t333 * MDP(22) + (-t432 * t294 + t285 + t290 * t569 - t297 * t358 - g(1) * (t341 * t435 - t379 * t432) - g(2) * (t339 * t435 - t377 * t432) - g(3) * (t365 * t435 + t403)) * MDP(23) + (-t435 * t294 - t432 * t286 + t289 * t569 + t297 * t356 - g(1) * (-t341 * t432 - t379 * t435) - g(2) * (-t339 * t432 - t377 * t435) - g(3) * (-t365 * t432 + t428 * t523)) * MDP(24) + (-MDP(20) * t509 - MDP(21) * t358 - MDP(23) * t290 - MDP(24) * t289) * qJD(6);];
tau  = t1;
