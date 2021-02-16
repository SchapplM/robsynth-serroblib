% Calculate vector of inverse dynamics joint torques for
% S5RRRPR1
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:49:12
% EndTime: 2021-01-15 22:49:29
% DurationCPUTime: 6.39s
% Computational Cost: add. (4804->400), mult. (11642->514), div. (0->0), fcn. (8650->16), ass. (0->189)
t496 = sin(qJ(3));
t500 = cos(qJ(3));
t501 = cos(qJ(2));
t553 = qJD(1) * t501;
t497 = sin(qJ(2));
t554 = qJD(1) * t497;
t419 = -t496 * t553 - t500 * t554;
t413 = t419 * qJ(4);
t583 = pkin(6) + pkin(7);
t450 = t583 * t501;
t436 = qJD(1) * t450;
t420 = t496 * t436;
t449 = t583 * t497;
t434 = qJD(1) * t449;
t574 = qJD(2) * pkin(2);
t426 = -t434 + t574;
t537 = t500 * t426 - t420;
t372 = t413 + t537;
t489 = qJD(2) + qJD(3);
t362 = pkin(3) * t489 + t372;
t493 = sin(pkin(9));
t424 = t500 * t436;
t523 = -t426 * t496 - t424;
t542 = t500 * t553;
t543 = t496 * t554;
t417 = -t542 + t543;
t573 = qJ(4) * t417;
t373 = -t523 - t573;
t494 = cos(pkin(9));
t567 = t494 * t373;
t326 = t493 * t362 + t567;
t390 = -t494 * t417 + t419 * t493;
t580 = pkin(8) * t390;
t318 = t326 + t580;
t499 = cos(qJ(5));
t564 = t390 * t499;
t495 = sin(qJ(5));
t524 = -t417 * t493 - t494 * t419;
t569 = t524 * t495;
t347 = t564 - t569;
t486 = t501 * pkin(2);
t576 = pkin(1) + t486;
t448 = t576 * qJD(1);
t401 = pkin(3) * t417 + qJD(4) - t448;
t357 = -pkin(4) * t390 + t401;
t492 = qJ(2) + qJ(3);
t483 = pkin(9) + t492;
t474 = qJ(5) + t483;
t463 = sin(t474);
t464 = cos(t474);
t498 = sin(qJ(1));
t502 = cos(qJ(1));
t532 = g(1) * t502 + g(2) * t498;
t550 = qJD(5) * t495;
t592 = g(3) * t463 + t318 * t550 - t357 * t347 + t464 * t532;
t547 = qJDD(1) * t501;
t549 = qJD(1) * qJD(2);
t540 = t501 * t549;
t548 = qJDD(1) * t497;
t586 = t540 + t548;
t374 = qJD(3) * t542 - t489 * t543 + t496 * t547 + t500 * t586;
t487 = qJDD(2) + qJDD(3);
t400 = qJDD(2) * pkin(2) - t583 * t586;
t541 = t497 * t549;
t402 = t583 * (-t541 + t547);
t509 = t523 * qJD(3) + t500 * t400 - t496 * t402;
t313 = pkin(3) * t487 - qJ(4) * t374 + qJD(4) * t419 + t509;
t429 = t496 * t501 + t497 * t500;
t399 = t489 * t429;
t530 = t496 * t548 - t500 * t547;
t375 = t399 * qJD(1) + t530;
t552 = qJD(3) * t496;
t584 = (qJD(3) * t426 + t402) * t500 + t496 * t400 - t436 * t552;
t315 = -qJ(4) * t375 - qJD(4) * t417 + t584;
t299 = t494 * t313 - t315 * t493;
t333 = t374 * t494 - t375 * t493;
t295 = pkin(4) * t487 - pkin(8) * t333 + t299;
t300 = t493 * t313 + t494 * t315;
t332 = t374 * t493 + t494 * t375;
t296 = -pkin(8) * t332 + t300;
t591 = -t495 * t295 - t499 * t296 + t592;
t482 = qJD(5) + t489;
t571 = t347 * t482;
t344 = t390 * t495 + t499 * t524;
t481 = qJDD(5) + t487;
t590 = t481 * MDP(26) + (t344 ^ 2 - t347 ^ 2) * MDP(23) - t347 * MDP(22) * t344;
t572 = t344 * t482;
t511 = -g(3) * t464 + t499 * t295 - t495 * t296 - t357 * t344 + t463 * t532;
t484 = sin(t492);
t485 = cos(t492);
t588 = -g(3) * t485 + t532 * t484;
t385 = t524 * pkin(8);
t536 = t434 * t496 - t424;
t377 = t536 + t573;
t559 = -t500 * t434 - t420;
t378 = t413 + t559;
t566 = t494 * t496;
t575 = pkin(2) * qJD(3);
t561 = -t494 * t377 + t378 * t493 + (-t493 * t500 - t566) * t575;
t568 = t493 * t496;
t560 = -t493 * t377 - t494 * t378 + (t494 * t500 - t568) * t575;
t558 = -t496 * t449 + t500 * t450;
t539 = t499 * t332 + t333 * t495;
t303 = qJD(5) * t344 + t539;
t582 = pkin(3) * t419;
t581 = pkin(3) * t493;
t363 = t493 * t373;
t325 = t494 * t362 - t363;
t316 = pkin(4) * t489 + t325 - t385;
t565 = t499 * t316;
t563 = t580 + t561;
t562 = -t385 - t560;
t428 = t496 * t497 - t500 * t501;
t544 = qJD(2) * t583;
t435 = t497 * t544;
t437 = t501 * t544;
t551 = qJD(3) * t500;
t517 = -t500 * t435 - t496 * t437 - t449 * t551 - t450 * t552;
t338 = -qJ(4) * t399 - qJD(4) * t428 + t517;
t398 = t489 * t428;
t508 = -t558 * qJD(3) + t435 * t496 - t500 * t437;
t339 = qJ(4) * t398 - qJD(4) * t429 + t508;
t309 = t494 * t338 + t493 * t339;
t331 = t494 * t372 - t363;
t535 = -t500 * t449 - t450 * t496;
t386 = -qJ(4) * t429 + t535;
t387 = -qJ(4) * t428 + t558;
t343 = t493 * t386 + t494 * t387;
t557 = pkin(3) * t485 + t486;
t490 = t497 ^ 2;
t556 = -t501 ^ 2 + t490;
t480 = t497 * t574;
t546 = qJD(5) * t564 - t495 * t332 + t499 * t333;
t392 = pkin(3) * t399 + t480;
t308 = -t338 * t493 + t494 * t339;
t330 = -t372 * t493 - t567;
t342 = t494 * t386 - t387 * t493;
t404 = pkin(3) * t428 - t576;
t477 = pkin(2) * t500 + pkin(3);
t411 = -pkin(2) * t568 + t494 * t477;
t361 = pkin(4) * t524 - t582;
t531 = g(1) * t498 - g(2) * t502;
t529 = -t495 * t316 - t499 * t318;
t396 = -t428 * t493 + t429 * t494;
t323 = -pkin(8) * t396 + t342;
t395 = t494 * t428 + t429 * t493;
t324 = -pkin(8) * t395 + t343;
t528 = t323 * t499 - t324 * t495;
t527 = t323 * t495 + t324 * t499;
t526 = t325 * t390 + t326 * t524;
t349 = t499 * t395 + t396 * t495;
t350 = -t395 * t495 + t396 * t499;
t471 = pkin(3) * t494 + pkin(4);
t522 = t471 * t495 + t499 * t581;
t521 = t471 * t499 - t495 * t581;
t520 = -0.2e1 * pkin(1) * t549 - pkin(6) * qJDD(2);
t414 = pkin(2) * t541 - qJDD(1) * t576;
t302 = -t524 * t550 + t546;
t469 = sin(t483);
t470 = cos(t483);
t516 = g(3) * t469 - t390 * t401 + t470 * t532 - t300;
t503 = qJD(2) ^ 2;
t513 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t503 + t531;
t504 = qJD(1) ^ 2;
t512 = pkin(1) * t504 - pkin(6) * qJDD(1) + t532;
t354 = pkin(3) * t375 + qJDD(4) + t414;
t510 = -g(3) * t470 - t401 * t524 + t469 * t532 + t299;
t507 = g(3) * t484 - t448 * t417 + t485 * t532 - t584;
t506 = -t419 * t417 * MDP(11) + (t302 - t571) * MDP(24) + (-t303 + t572) * MDP(25) + (t417 * t489 + t374) * MDP(13) + (-t530 + (-qJD(1) * t429 - t419) * t489) * MDP(14) + (-t417 ^ 2 + t419 ^ 2) * MDP(12) + t487 * MDP(15) + t590;
t505 = -t448 * t419 + t509 + t588;
t488 = -qJ(4) - t583;
t479 = pkin(2) * t554;
t433 = pkin(1) + t557;
t412 = pkin(2) * t566 + t477 * t493;
t407 = pkin(4) + t411;
t403 = t479 - t582;
t367 = pkin(4) * t395 + t404;
t358 = t361 + t479;
t356 = -t398 * t494 - t399 * t493;
t355 = -t398 * t493 + t494 * t399;
t337 = pkin(4) * t355 + t392;
t320 = -t385 + t331;
t319 = t330 - t580;
t310 = pkin(4) * t332 + t354;
t307 = t350 * qJD(5) + t499 * t355 + t356 * t495;
t306 = -t349 * qJD(5) - t355 * t495 + t356 * t499;
t305 = -pkin(8) * t355 + t309;
t304 = -pkin(8) * t356 + t308;
t1 = [(-t302 * t349 - t303 * t350 + t306 * t347 - t307 * t344) * MDP(23) + (-t337 * t347 + t367 * t303 + t310 * t349 + t357 * t307 + (-t527 * qJD(5) + t304 * t499 - t305 * t495) * t482 + t528 * t481 + t531 * t464) * MDP(27) + (-t375 * t576 - t448 * t399 + t414 * t428 + t417 * t480 + t531 * t485 + t535 * t487 + t508 * t489) * MDP(16) + (qJDD(1) * t490 + 0.2e1 * t497 * t540) * MDP(4) + (-t299 * t396 - t300 * t395 - t308 * t524 + t309 * t390 - t325 * t356 - t326 * t355 - t332 * t343 - t333 * t342 - t532) * MDP(20) + (t308 * t489 + t332 * t404 + t342 * t487 + t354 * t395 + t355 * t401 - t390 * t392 + t531 * t470) * MDP(18) + (-t374 * t576 + t448 * t398 + t414 * t429 - t419 * t480 - t531 * t484 - t558 * t487 - t517 * t489) * MDP(17) + (-t398 * t489 + t429 * t487) * MDP(13) + (-t374 * t428 - t375 * t429 + t398 * t417 + t399 * t419) * MDP(12) + (t374 * t429 + t398 * t419) * MDP(11) + (t337 * t344 + t367 * t302 + t310 * t350 + t357 * t306 - (t528 * qJD(5) + t304 * t495 + t305 * t499) * t482 - t527 * t481 - t531 * t463) * MDP(28) + (t302 * t350 + t306 * t344) * MDP(22) + (-t309 * t489 + t333 * t404 - t343 * t487 + t354 * t396 + t356 * t401 + t392 * t524 - t531 * t469) * MDP(19) + 0.2e1 * (t497 * t547 - t556 * t549) * MDP(5) + t532 * MDP(3) + t531 * MDP(2) + qJDD(1) * MDP(1) + (t300 * t343 + t326 * t309 + t299 * t342 + t325 * t308 + t354 * t404 + t401 * t392 - g(1) * (-t433 * t498 - t488 * t502) - g(2) * (t433 * t502 - t488 * t498)) * MDP(21) + (qJDD(2) * t497 + t501 * t503) * MDP(6) + (qJDD(2) * t501 - t497 * t503) * MDP(7) + (-t399 * t489 - t428 * t487) * MDP(14) + (t306 * t482 + t350 * t481) * MDP(24) + (-t307 * t482 - t349 * t481) * MDP(25) + (t520 * t497 + t513 * t501) * MDP(9) + (-t513 * t497 + t520 * t501) * MDP(10); (t300 * t412 + t299 * t411 - t401 * t403 - g(3) * t557 - t532 * (-pkin(2) * t497 - pkin(3) * t484) + t560 * t326 + t561 * t325) * MDP(21) + (t559 * t489 + (t419 * t554 - t496 * t487 - t489 * t551) * pkin(2) + t507) * MDP(17) + (t390 * t403 + t411 * t487 + t561 * t489 + t510) * MDP(18) + (-t403 * t524 - t412 * t487 - t560 * t489 + t516) * MDP(19) + (-t332 * t412 - t333 * t411 + t390 * t560 - t524 * t561 + t526) * MDP(20) + (-g(3) * t501 + t512 * t497) * MDP(9) + (g(3) * t497 + t512 * t501) * MDP(10) + (-t536 * t489 + (-t417 * t554 + t487 * t500 - t489 * t552) * pkin(2) + t505) * MDP(16) + ((t407 * t499 - t412 * t495) * t481 + t358 * t347 + (t562 * t495 + t563 * t499) * t482 + ((-t407 * t495 - t412 * t499) * t482 + t529) * qJD(5) + t511) * MDP(27) + t506 + qJDD(2) * MDP(8) + MDP(6) * t548 + MDP(7) * t547 + (-t358 * t344 + (-t407 * t481 - t295 + (qJD(5) * t412 - t563) * t482) * t495 + (-qJD(5) * t316 - t412 * t481 - t296 + (-qJD(5) * t407 + t562) * t482) * t499 + t592) * MDP(28) + (-t497 * t501 * MDP(4) + t556 * MDP(5)) * t504; (-t330 * t489 + (-t390 * t419 + t487 * t494) * pkin(3) + t510) * MDP(18) + (t331 * t489 + (t419 * t524 - t487 * t493) * pkin(3) + t516) * MDP(19) + (-t325 * t330 - t326 * t331 + (t299 * t494 + t300 * t493 + t401 * t419 + t588) * pkin(3)) * MDP(21) + (t330 * t524 - t331 * t390 + (-t332 * t493 - t333 * t494) * pkin(3) + t526) * MDP(20) + (-t523 * t489 + t505) * MDP(16) + (t521 * t481 + t361 * t347 - (t319 * t499 - t320 * t495) * t482 + (-t522 * t482 + t529) * qJD(5) + t511) * MDP(27) + (-t522 * t481 - t361 * t344 + (t319 * t495 + t320 * t499) * t482 + (-t521 * t482 - t565) * qJD(5) + t591) * MDP(28) + t506 + (t537 * t489 + t507) * MDP(17); (t489 * t524 + t332) * MDP(18) + (t390 * t489 + t333) * MDP(19) + (-t390 ^ 2 - t524 ^ 2) * MDP(20) + (t325 * t524 - t326 * t390 + t354 - t531) * MDP(21) + (t303 + t572) * MDP(27) + (t302 + t571) * MDP(28); (t546 - t571) * MDP(24) + (-t539 + t572) * MDP(25) + (-t529 * t482 + t511) * MDP(27) + ((-t318 * t495 + t565) * t482 + t591) * MDP(28) + (-MDP(24) * t569 - MDP(25) * t344 + t529 * MDP(27) - MDP(28) * t565) * qJD(5) + t590;];
tau = t1;
