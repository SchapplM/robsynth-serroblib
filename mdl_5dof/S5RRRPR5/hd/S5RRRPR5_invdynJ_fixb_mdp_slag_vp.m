% Calculate vector of inverse dynamics joint torques for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:26
% EndTime: 2019-12-31 21:14:33
% DurationCPUTime: 4.46s
% Computational Cost: add. (4445->387), mult. (10730->514), div. (0->0), fcn. (7951->14), ass. (0->193)
t451 = qJ(2) + qJ(3);
t443 = sin(t451);
t444 = cos(t451);
t457 = sin(qJ(1));
t461 = cos(qJ(1));
t492 = g(1) * t461 + g(2) * t457;
t567 = -g(3) * t444 + t443 * t492;
t448 = qJD(2) + qJD(3);
t454 = sin(qJ(5));
t458 = cos(qJ(5));
t459 = cos(qJ(3));
t460 = cos(qJ(2));
t523 = qJD(1) * t460;
t510 = t459 * t523;
t455 = sin(qJ(3));
t456 = sin(qJ(2));
t524 = qJD(1) * t456;
t511 = t455 * t524;
t383 = -t510 + t511;
t385 = -t455 * t523 - t459 * t524;
t452 = sin(pkin(9));
t453 = cos(pkin(9));
t487 = -t383 * t452 - t453 * t385;
t346 = -t458 * t448 + t454 * t487;
t504 = -t453 * t383 + t385 * t452;
t558 = qJD(5) - t504;
t566 = t346 * t558;
t348 = t448 * t454 + t458 * t487;
t565 = t348 * t558;
t515 = qJDD(1) * t460;
t517 = qJD(1) * qJD(2);
t508 = t460 * t517;
t516 = qJDD(1) * t456;
t562 = t508 + t516;
t343 = qJD(3) * t510 - t448 * t511 + t455 * t515 + t562 * t459;
t446 = qJDD(2) + qJDD(3);
t555 = pkin(6) + pkin(7);
t363 = qJDD(2) * pkin(2) - t555 * t562;
t509 = t456 * t517;
t365 = t555 * (-t509 + t515);
t411 = t555 * t460;
t404 = qJD(1) * t411;
t390 = t459 * t404;
t410 = t555 * t456;
t402 = qJD(1) * t410;
t546 = qJD(2) * pkin(2);
t392 = -t402 + t546;
t486 = -t392 * t455 - t390;
t469 = qJD(3) * t486 + t459 * t363 - t455 * t365;
t302 = pkin(3) * t446 - qJ(4) * t343 + qJD(4) * t385 + t469;
t397 = t455 * t460 + t456 * t459;
t362 = t448 * t397;
t490 = t455 * t516 - t459 * t515;
t344 = qJD(1) * t362 + t490;
t522 = qJD(3) * t455;
t557 = (qJD(3) * t392 + t365) * t459 + t455 * t363 - t404 * t522;
t306 = -qJ(4) * t344 - qJD(4) * t383 + t557;
t289 = t302 * t453 - t306 * t452;
t287 = -pkin(4) * t446 - t289;
t442 = pkin(9) + t451;
t431 = cos(t442);
t563 = g(3) * t431 + t287;
t500 = t558 * t458;
t318 = -t343 * t452 - t453 * t344;
t315 = qJDD(5) - t318;
t535 = t454 * t315;
t561 = -t558 * t500 - t535;
t445 = t460 * pkin(2);
t548 = pkin(1) + t445;
t554 = pkin(3) * t385;
t326 = pkin(4) * t487 - pkin(8) * t504 - t554;
t432 = pkin(3) * t452 + pkin(8);
t560 = t558 * (qJD(5) * t432 + t326);
t379 = t385 * qJ(4);
t386 = t455 * t404;
t503 = t459 * t392 - t386;
t341 = t379 + t503;
t527 = -t455 * t410 + t459 * t411;
t396 = t455 * t456 - t459 * t460;
t512 = qJD(2) * t555;
t403 = t456 * t512;
t405 = t460 * t512;
t521 = qJD(3) * t459;
t475 = -t459 * t403 - t455 * t405 - t410 * t521 - t411 * t522;
t322 = -qJ(4) * t362 - qJD(4) * t396 + t475;
t361 = t448 * t396;
t468 = -qJD(3) * t527 + t403 * t455 - t459 * t405;
t465 = qJ(4) * t361 - qJD(4) * t397 + t468;
t299 = t453 * t322 + t452 * t465;
t336 = pkin(3) * t448 + t341;
t545 = qJ(4) * t383;
t342 = -t486 - t545;
t540 = t342 * t452;
t311 = t336 * t453 - t540;
t309 = -pkin(4) * t448 - t311;
t358 = t453 * t396 + t397 * t452;
t359 = -t396 * t452 + t397 * t453;
t493 = pkin(3) * t396 - t548;
t327 = pkin(4) * t358 - pkin(8) * t359 + t493;
t352 = -qJ(4) * t396 + t527;
t501 = -t459 * t410 - t411 * t455;
t479 = -qJ(4) * t397 + t501;
t329 = t453 * t352 + t452 * t479;
t334 = -t361 * t453 - t362 * t452;
t290 = t452 * t302 + t453 * t306;
t288 = pkin(8) * t446 + t290;
t409 = t548 * qJD(1);
t364 = pkin(3) * t383 + qJD(4) - t409;
t321 = -pkin(4) * t504 - pkin(8) * t487 + t364;
t498 = qJD(5) * t321 + t288;
t556 = t287 * t359 + t309 * t334 - t329 * t315 - (qJD(5) * t327 + t299) * t558 - t358 * t498;
t430 = sin(t442);
t551 = g(3) * t430;
t547 = pkin(2) * qJD(3);
t319 = t343 * t453 - t344 * t452;
t519 = qJD(5) * t458;
t514 = t458 * t319 + t454 * t446 + t448 * t519;
t520 = qJD(5) * t454;
t304 = -t487 * t520 + t514;
t544 = t304 * t454;
t543 = t309 * t504;
t542 = t309 * t359;
t541 = t327 * t315;
t539 = t346 * t487;
t538 = t348 * t487;
t537 = t452 * t455;
t337 = t453 * t342;
t536 = t453 * t455;
t534 = t454 * t457;
t533 = t454 * t461;
t532 = t457 * t458;
t313 = t458 * t315;
t531 = t458 * t461;
t312 = t452 * t336 + t337;
t528 = -t459 * t402 - t386;
t349 = t379 + t528;
t502 = t402 * t455 - t390;
t480 = t502 + t545;
t530 = -t349 * t452 + t453 * t480 + (t452 * t459 + t536) * t547;
t529 = -t453 * t349 - t452 * t480 + (t453 * t459 - t537) * t547;
t438 = pkin(2) * t459 + pkin(3);
t378 = pkin(2) * t536 + t452 * t438;
t526 = pkin(3) * t444 + t445;
t449 = t456 ^ 2;
t525 = -t460 ^ 2 + t449;
t441 = t456 * t546;
t507 = pkin(3) * t362 + t441;
t380 = pkin(2) * t509 - qJDD(1) * t548;
t470 = pkin(3) * t344 + qJDD(4) + t380;
t295 = -pkin(4) * t318 - pkin(8) * t319 + t470;
t310 = pkin(8) * t448 + t312;
t497 = qJD(5) * t310 - t295;
t373 = pkin(8) + t378;
t440 = pkin(2) * t524;
t495 = qJD(5) * t373 + t326 + t440;
t491 = g(1) * t457 - g(2) * t461;
t489 = -t315 * t373 - t543;
t488 = t311 * t504 + t312 * t487;
t297 = t310 * t458 + t321 * t454;
t485 = t297 * t487 + t309 * t519 + t563 * t454;
t296 = -t310 * t454 + t321 * t458;
t484 = -t296 * t487 + t309 * t520 + (g(1) * t531 + g(2) * t532) * t430;
t483 = t313 + (t454 * t504 - t520) * t558;
t377 = -pkin(2) * t537 + t438 * t453;
t482 = t492 * t430;
t481 = -0.2e1 * pkin(1) * t517 - pkin(6) * qJDD(2);
t478 = t334 * t458 - t359 * t520;
t317 = t341 * t453 - t540;
t473 = -t315 * t432 + t317 * t558 - t543;
t462 = qJD(2) ^ 2;
t472 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t462 + t491;
t463 = qJD(1) ^ 2;
t471 = pkin(1) * t463 - pkin(6) * qJDD(1) + t492;
t467 = g(3) * t443 - t409 * t383 + t492 * t444 - t557;
t427 = t458 * t446;
t305 = qJD(5) * t348 + t319 * t454 - t427;
t466 = -t385 * t383 * MDP(11) - t558 * t487 * MDP(24) + ((t304 - t566) * t458 + (-t305 - t565) * t454) * MDP(21) + (t483 + t539) * MDP(23) + (-t538 - t561) * MDP(22) + (t348 * t500 + t544) * MDP(20) + (t383 * t448 + t343) * MDP(13) + (-t490 + (-qJD(1) * t397 - t385) * t448) * MDP(14) + (-t383 ^ 2 + t385 ^ 2) * MDP(12) + t446 * MDP(15);
t464 = -t409 * t385 + t469 + t567;
t447 = -qJ(4) - t555;
t433 = -pkin(3) * t453 - pkin(4);
t401 = pkin(1) + t526;
t372 = -pkin(4) - t377;
t371 = t431 * t531 + t534;
t370 = -t431 * t533 + t532;
t369 = -t431 * t532 + t533;
t368 = t431 * t534 + t531;
t333 = -t361 * t452 + t453 * t362;
t328 = t352 * t452 - t453 * t479;
t316 = t341 * t452 + t337;
t301 = pkin(4) * t333 - pkin(8) * t334 + t507;
t298 = t322 * t452 - t453 * t465;
t294 = t458 * t295;
t1 = [0.2e1 * (t456 * t515 - t517 * t525) * MDP(5) + (qJDD(2) * t456 + t460 * t462) * MDP(6) + (qJDD(2) * t460 - t456 * t462) * MDP(7) + (t456 * t481 + t460 * t472) * MDP(9) + (-t456 * t472 + t460 * t481) * MDP(10) + (t343 * t397 + t361 * t385) * MDP(11) + (-t343 * t396 - t344 * t397 + t361 * t383 + t362 * t385) * MDP(12) + (-t361 * t448 + t397 * t446) * MDP(13) + (-t362 * t448 - t396 * t446) * MDP(14) + (-t344 * t548 - t409 * t362 + t380 * t396 + t383 * t441 + t444 * t491 + t446 * t501 + t448 * t468) * MDP(16) + (-t343 * t548 + t409 * t361 + t380 * t397 - t385 * t441 - t443 * t491 - t446 * t527 - t448 * t475) * MDP(17) + (-t289 * t359 - t290 * t358 + t298 * t487 + t299 * t504 - t311 * t334 - t312 * t333 + t318 * t329 + t319 * t328 - t492) * MDP(18) + (t290 * t329 + t312 * t299 - t289 * t328 - t311 * t298 + t470 * t493 + t364 * t507 - g(1) * (-t401 * t457 - t447 * t461) - g(2) * (t401 * t461 - t447 * t457)) * MDP(19) + (t304 * t359 * t458 + t348 * t478) * MDP(20) + ((-t346 * t458 - t348 * t454) * t334 + (-t544 - t305 * t458 + (t346 * t454 - t348 * t458) * qJD(5)) * t359) * MDP(21) + (t304 * t358 + t313 * t359 + t333 * t348 + t478 * t558) * MDP(22) + (-t359 * t535 - t305 * t358 - t333 * t346 + (-t334 * t454 - t359 * t519) * t558) * MDP(23) + (t315 * t358 + t333 * t558) * MDP(24) + (-g(1) * t369 - g(2) * t371 + t294 * t358 + t296 * t333 + t298 * t346 + t328 * t305 + (t301 * t558 + t541 + (-t310 * t358 - t329 * t558 + t542) * qJD(5)) * t458 + t556 * t454) * MDP(25) + (-g(1) * t368 - g(2) * t370 - t297 * t333 + t298 * t348 + t328 * t304 + (-(-qJD(5) * t329 + t301) * t558 - t541 + t497 * t358 - qJD(5) * t542) * t454 + t556 * t458) * MDP(26) + qJDD(1) * MDP(1) + (qJDD(1) * t449 + 0.2e1 * t456 * t508) * MDP(4) + t491 * MDP(2) + t492 * MDP(3); (t318 * t378 - t319 * t377 + t487 * t530 + t504 * t529 + t488) * MDP(18) + (-t502 * t448 + (-t383 * t524 + t446 * t459 - t448 * t522) * pkin(2) + t464) * MDP(16) + qJDD(2) * MDP(8) + (-g(3) * t460 + t456 * t471) * MDP(9) + (g(3) * t456 + t460 * t471) * MDP(10) + (t290 * t378 + t289 * t377 - t364 * (t440 - t554) - g(3) * t526 - t492 * (-pkin(2) * t456 - pkin(3) * t443) + t529 * t312 - t530 * t311) * MDP(19) + (t528 * t448 + (t385 * t524 - t455 * t446 - t448 * t521) * pkin(2) + t467) * MDP(17) + (t372 * t305 - t563 * t458 + t489 * t454 + t530 * t346 + (-t454 * t529 - t458 * t495) * t558 + t484) * MDP(25) + (t372 * t304 + t489 * t458 - t454 * t482 + t530 * t348 + (t454 * t495 - t458 * t529) * t558 + t485) * MDP(26) + MDP(6) * t516 + t466 + MDP(7) * t515 + (-t456 * t460 * MDP(4) + MDP(5) * t525) * t463; (t433 * t305 - t316 * t346 + t473 * t454 + (-t563 - t560) * t458 + t484) * MDP(25) + (-t448 * t486 + t464) * MDP(16) + (t448 * t503 + t467) * MDP(17) + (-t316 * t487 - t317 * t504 + (t318 * t452 - t319 * t453) * pkin(3) + t488) * MDP(18) + (t311 * t316 - t312 * t317 + (t289 * t453 + t290 * t452 + t364 * t385 + t567) * pkin(3)) * MDP(19) + (t433 * t304 - t316 * t348 + t473 * t458 + (-t482 + t560) * t454 + t485) * MDP(26) + t466; (-t487 ^ 2 - t504 ^ 2) * MDP(18) + (t311 * t487 - t312 * t504 + t470 - t491) * MDP(19) + (t483 - t539) * MDP(25) + (-t538 + t561) * MDP(26); t348 * t346 * MDP(20) + (-t346 ^ 2 + t348 ^ 2) * MDP(21) + (t514 + t566) * MDP(22) + (t427 + t565) * MDP(23) + t315 * MDP(24) + (-g(1) * t370 + g(2) * t368 + t297 * t558 - t309 * t348 + t294) * MDP(25) + (g(1) * t371 - g(2) * t369 + t296 * t558 + t309 * t346) * MDP(26) + ((-t288 + t551) * MDP(26) + (-MDP(23) * t487 - MDP(25) * t310 - MDP(26) * t321) * qJD(5)) * t458 + (-qJD(5) * t487 * MDP(22) + (-qJD(5) * t448 - t319) * MDP(23) + (-t498 + t551) * MDP(25) + t497 * MDP(26)) * t454;];
tau = t1;
