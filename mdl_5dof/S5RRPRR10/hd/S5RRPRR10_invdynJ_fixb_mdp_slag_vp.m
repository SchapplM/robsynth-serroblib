% Calculate vector of inverse dynamics joint torques for
% S5RRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:53
% EndTime: 2019-12-31 20:27:04
% DurationCPUTime: 6.55s
% Computational Cost: add. (5260->513), mult. (15123->721), div. (0->0), fcn. (12484->12), ass. (0->241)
t493 = sin(pkin(5));
t502 = cos(qJ(2));
t488 = pkin(2) * t502 + pkin(1);
t537 = t488 * qJDD(1);
t498 = sin(qJ(2));
t584 = qJD(2) * t498;
t563 = t493 * t584;
t542 = qJD(1) * t563;
t572 = pkin(2) * t542 + qJDD(3);
t639 = t493 * t537 - t572;
t492 = sin(pkin(10));
t585 = qJD(1) * t498;
t564 = t493 * t585;
t494 = cos(pkin(10));
t594 = t502 * t494;
t567 = t493 * t594;
t446 = qJD(1) * t567 - t492 * t564;
t442 = qJD(4) - t446;
t527 = t492 * t502 + t494 * t498;
t519 = qJD(1) * t527;
t449 = t493 * t519;
t501 = cos(qJ(4));
t495 = cos(pkin(5));
t586 = qJD(1) * t495;
t547 = qJD(2) + t586;
t470 = t501 * t547;
t497 = sin(qJ(4));
t421 = t449 * t497 - t470;
t420 = qJD(5) + t421;
t499 = sin(qJ(1));
t596 = t499 * t502;
t503 = cos(qJ(1));
t597 = t498 * t503;
t458 = -t495 * t596 - t597;
t593 = t502 * t503;
t598 = t498 * t499;
t525 = t495 * t593 - t598;
t638 = -g(1) * t458 - g(2) * t525;
t453 = t527 * t495;
t461 = t492 * t498 - t594;
t529 = t453 * t503 - t461 * t499;
t602 = t493 * t503;
t400 = -t497 * t602 + t501 * t529;
t522 = t461 * t495;
t413 = -t499 * t527 - t503 * t522;
t496 = sin(qJ(5));
t500 = cos(qJ(5));
t637 = t400 * t496 + t413 * t500;
t636 = t400 * t500 - t413 * t496;
t508 = -t501 * t449 - t497 * t547;
t611 = t508 * t496;
t373 = -t500 * t442 - t611;
t635 = t373 * t442;
t626 = pkin(1) * t495;
t483 = t502 * t626;
t623 = pkin(7) + qJ(3);
t560 = t623 * t498;
t543 = t493 * t560;
t514 = t495 * pkin(2) - t543;
t434 = t483 + t514;
t600 = t495 * t498;
t482 = pkin(1) * t600;
t603 = t493 * t502;
t588 = pkin(7) * t603 + t482;
t443 = qJ(3) * t603 + t588;
t390 = t492 * t434 + t494 * t443;
t381 = pkin(8) * t495 + t390;
t605 = t493 * t498;
t451 = t492 * t605 - t567;
t452 = t527 * t493;
t541 = t488 * t493;
t396 = pkin(3) * t451 - pkin(8) * t452 - t541;
t633 = t501 * t381 + t497 * t396;
t489 = t493 ^ 2;
t627 = pkin(1) * t489;
t632 = t495 * t588 + t498 * t627;
t573 = t502 * qJDD(1);
t556 = t493 * t573;
t466 = t494 * t556;
t518 = t527 * qJD(2);
t574 = t498 * qJDD(1);
t558 = t492 * t574;
t409 = qJDD(4) - t466 + (qJD(1) * t518 + t558) * t493;
t478 = qJD(1) * t483;
t424 = qJD(2) * pkin(2) + qJD(1) * t514 + t478;
t436 = (t603 * t623 + t482) * qJD(1);
t601 = t494 * t436;
t369 = t492 * t424 + t601;
t362 = pkin(8) * t547 + t369;
t455 = -qJD(1) * t541 + qJD(3);
t379 = -pkin(3) * t446 - pkin(8) * t449 + t455;
t341 = t362 * t501 + t379 * t497;
t569 = pkin(1) * t573;
t477 = t495 * t569;
t575 = qJDD(1) * t495;
t480 = qJDD(2) + t575;
t570 = qJD(2) * t626;
t544 = qJD(1) * t570;
t555 = qJD(2) * t623;
t583 = qJD(3) * t498;
t366 = -t498 * t544 + pkin(2) * t480 + t477 + (-qJDD(1) * t560 + (-t502 * t555 - t583) * qJD(1)) * t493;
t513 = qJD(3) * t502 - t498 * t555;
t565 = pkin(7) * t556 + qJDD(1) * t482 + t502 * t544;
t376 = (qJ(3) * t573 + qJD(1) * t513) * t493 + t565;
t339 = t492 * t366 + t494 * t376;
t336 = pkin(8) * t480 + t339;
t410 = t466 + (-qJD(2) * t519 - t558) * t493;
t576 = qJD(1) * qJD(2);
t559 = t502 * t576;
t411 = -t492 * t542 + (qJDD(1) * t527 + t494 * t559) * t493;
t351 = -pkin(3) * t410 - pkin(8) * t411 - t639;
t552 = t336 * t497 - t501 * t351;
t629 = qJD(4) * t341 + t552;
t322 = -pkin(4) * t409 + t629;
t528 = -t453 * t499 - t461 * t503;
t604 = t493 * t499;
t402 = -t497 * t528 + t501 * t604;
t431 = t452 * t497 - t495 * t501;
t526 = t497 * t529 + t501 * t602;
t516 = g(1) * t402 - g(2) * t526 - g(3) * t431;
t631 = t420 * (-pkin(4) * t508 + t420 * pkin(9)) + t322 + t516;
t356 = -qJD(4) * t508 + t497 * t411 - t501 * t480;
t354 = qJDD(5) + t356;
t435 = -qJD(1) * t543 + t478;
t382 = t435 * t492 + t601;
t487 = -pkin(2) * t494 - pkin(3);
t460 = -pkin(4) * t501 - pkin(9) * t497 + t487;
t630 = t420 * (t382 - t442 * (pkin(4) * t497 - pkin(9) * t501)) - t460 * t354;
t628 = 0.2e1 * t489;
t504 = qJD(1) ^ 2;
t624 = g(3) * t502;
t580 = qJD(4) * t497;
t355 = qJD(4) * t470 + t501 * t411 - t449 * t580 + t497 * t480;
t577 = qJD(5) * t500;
t566 = t500 * t355 + t496 * t409 + t442 * t577;
t578 = qJD(5) * t496;
t329 = t508 * t578 + t566;
t622 = t329 * t496;
t621 = t373 * t420;
t375 = t442 * t496 - t500 * t508;
t620 = t375 * t420;
t619 = t375 * t446;
t486 = pkin(2) * t492 + pkin(8);
t616 = t420 * t486;
t548 = t420 * t500;
t615 = t421 * t442;
t614 = t421 * t449;
t613 = t508 * t442;
t612 = t508 * t449;
t610 = t442 * t497;
t609 = t446 * t501;
t607 = t480 * MDP(8);
t606 = t489 * t504;
t428 = t492 * t436;
t599 = t496 * t354;
t595 = t500 * t354;
t383 = t435 * t494 - t428;
t392 = pkin(2) * t564 + pkin(3) * t449 - pkin(8) * t446;
t590 = t501 * t383 + t497 * t392;
t589 = t501 * t409 + t446 * t610;
t490 = t498 ^ 2;
t587 = -t502 ^ 2 + t490;
t582 = qJD(4) * t420;
t581 = qJD(4) * t486;
t579 = qJD(4) * t501;
t568 = t502 * t606;
t562 = t493 * t495 * t504;
t561 = t623 * t493;
t557 = t493 * t574;
t521 = t501 * t336 + t497 * t351 - t362 * t580 + t379 * t579;
t321 = pkin(9) * t409 + t521;
t338 = t366 * t494 - t492 * t376;
t335 = -pkin(3) * t480 - t338;
t324 = pkin(4) * t356 - pkin(9) * t355 + t335;
t553 = -t496 * t321 + t500 * t324;
t551 = t355 * t496 - t500 * t409;
t479 = t502 * t570;
t425 = t493 * t513 + t479;
t426 = -t493 * t583 + (-t502 * t561 - t482) * qJD(2);
t363 = t425 * t492 - t494 * t426;
t368 = t494 * t424 - t428;
t389 = t434 * t494 - t492 * t443;
t549 = t442 * t501;
t546 = qJD(2) + 0.2e1 * t586;
t545 = t480 + t575;
t539 = g(1) * t503 + g(2) * t499;
t538 = g(1) * t499 - g(2) * t503;
t536 = t500 * t321 + t496 * t324;
t333 = pkin(9) * t442 + t341;
t361 = -pkin(3) * t547 - t368;
t337 = t421 * pkin(4) + pkin(9) * t508 + t361;
t326 = t333 * t500 + t337 * t496;
t535 = t333 * t496 - t337 * t500;
t347 = pkin(9) * t451 + t633;
t380 = -pkin(3) * t495 - t389;
t432 = t452 * t501 + t495 * t497;
t350 = pkin(4) * t431 - pkin(9) * t432 + t380;
t534 = t347 * t500 + t350 * t496;
t533 = -t347 * t496 + t350 * t500;
t340 = -t362 * t497 + t379 * t501;
t364 = t425 * t494 + t426 * t492;
t447 = t493 * t518;
t448 = t461 * t493 * qJD(2);
t393 = pkin(2) * t563 + pkin(3) * t447 + pkin(8) * t448;
t532 = -t364 * t497 + t393 * t501;
t530 = -t381 * t497 + t396 * t501;
t395 = t432 * t500 + t451 * t496;
t394 = t432 * t496 - t451 * t500;
t524 = -t420 * t577 - t599;
t520 = t501 * t364 - t381 * t580 + t497 * t393 + t396 * t579;
t517 = t361 * t442 - t486 * t409;
t416 = t499 * t522 - t503 * t527;
t515 = g(1) * t416 + g(2) * t413 - g(3) * t451;
t511 = -t335 - t515;
t398 = t449 * t496 + t500 * t609;
t509 = -t497 * t578 + t500 * t579 - t398;
t332 = -pkin(4) * t442 - t340;
t507 = -pkin(9) * t354 + (t332 + t340) * t420;
t506 = qJD(5) * t616 + t515;
t505 = -g(1) * t528 - g(2) * t529 - g(3) * t452 + (pkin(9) * t449 - qJD(5) * t460 + t590) * t420;
t459 = -t495 * t598 + t593;
t457 = -t495 * t597 - t596;
t454 = pkin(2) * t600 - t561;
t403 = t497 * t604 + t501 * t528;
t397 = -t500 * t449 + t496 * t609;
t388 = -qJD(4) * t431 - t448 * t501;
t387 = qJD(4) * t432 - t448 * t497;
t367 = t375 * t580;
t358 = t403 * t500 - t416 * t496;
t357 = -t403 * t496 - t416 * t500;
t346 = -pkin(4) * t451 - t530;
t345 = -qJD(5) * t394 + t388 * t500 + t447 * t496;
t344 = qJD(5) * t395 + t388 * t496 - t447 * t500;
t342 = -pkin(4) * t449 + t383 * t497 - t392 * t501;
t331 = pkin(4) * t387 - pkin(9) * t388 + t363;
t330 = qJD(5) * t375 + t551;
t328 = -pkin(4) * t447 + qJD(4) * t633 - t532;
t327 = pkin(9) * t447 + t520;
t320 = -t326 * qJD(5) + t553;
t319 = -t535 * qJD(5) + t536;
t1 = [(t339 * t390 + t369 * t364 + t338 * t389 - t368 * t363 - g(1) * (-t454 * t503 - t488 * t499) - g(2) * (-t454 * t499 + t488 * t503)) * MDP(12) + (t355 * t451 + t388 * t442 + t409 * t432 - t447 * t508) * MDP(15) + (-t355 * t431 - t356 * t432 + t387 * t508 - t388 * t421) * MDP(14) + (t355 * t432 - t388 * t508) * MDP(13) + (-(-pkin(7) * t563 + t479) * t547 - t588 * t480 - (-pkin(7) * t542 + t565) * t495 + g(1) * t525 - g(2) * t458 + 0.2e1 * (-t559 - t574) * t627) * MDP(10) + (t569 * t628 + (-pkin(7) * t605 + t483) * t480 + (-pkin(7) * t557 + t477) * t495 - g(1) * t457 - g(2) * t459 - t588 * qJD(2) ^ 2 - 0.2e1 * t632 * t576) * MDP(9) + (t532 * t442 + t530 * t409 - t552 * t451 + t340 * t447 + t363 * t421 + t380 * t356 + t335 * t431 + t361 * t387 + g(1) * t400 - g(2) * t403 + (-t341 * t451 - t442 * t633) * qJD(4)) * MDP(18) + (-g(1) * t526 - g(2) * t402 + t335 * t432 - t341 * t447 + t380 * t355 + t361 * t388 - t363 * t508 - t409 * t633 - t442 * t520 - t451 * t521) * MDP(19) + (-t338 * t452 - t339 * t451 + t363 * t449 + t364 * t446 + t368 * t448 - t369 * t447 - t389 * t411 + t390 * t410) * MDP(11) + ((-qJD(5) * t534 - t327 * t496 + t331 * t500) * t420 + t533 * t354 + t320 * t431 - t535 * t387 + t328 * t373 + t346 * t330 + t322 * t394 + t332 * t344 + g(1) * t636 - g(2) * t358) * MDP(25) + (-(qJD(5) * t533 + t327 * t500 + t331 * t496) * t420 - t534 * t354 - t319 * t431 - t326 * t387 + t328 * t375 + t346 * t329 + t322 * t395 + t332 * t345 - g(1) * t637 - g(2) * t357) * MDP(26) + (-t356 * t451 - t387 * t442 - t409 * t431 - t421 * t447) * MDP(16) + (t409 * t451 + t442 * t447) * MDP(17) + (t329 * t431 + t345 * t420 + t354 * t395 + t375 * t387) * MDP(22) + (t354 * t431 + t387 * t420) * MDP(24) + (-t330 * t431 - t344 * t420 - t354 * t394 - t373 * t387) * MDP(23) + (-t329 * t394 - t330 * t395 - t344 * t375 - t345 * t373) * MDP(21) + (t329 * t395 + t345 * t375) * MDP(20) + qJDD(1) * MDP(1) + (qJDD(1) * t490 + 0.2e1 * t498 * t559) * t489 * MDP(4) + ((pkin(2) * t455 * t584 + t639 * t488) * MDP(12) - t539 * MDP(11) + (qJD(2) * t502 * t546 + t498 * t545) * MDP(6) + (t502 * t545 - t546 * t584) * MDP(7)) * t493 + (t498 * t573 - t576 * t587) * MDP(5) * t628 + t538 * MDP(2) + t539 * MDP(3) + t495 * t607; -t498 * MDP(4) * t568 + t587 * MDP(5) * t606 + (-t502 * t562 + t557) * MDP(6) + (t498 * t562 + t556) * MDP(7) + t607 + (t477 + (-pkin(7) * t574 - t624) * t493 + t632 * t504 + t638) * MDP(9) + (pkin(1) * t568 + (-pkin(7) * t564 + t478) * t586 + g(1) * t459 - g(2) * t457 + g(3) * t605 + t478 * qJD(2) - t565) * MDP(10) + ((t369 - t382) * t449 + (t368 - t383) * t446 + (t410 * t492 - t411 * t494) * pkin(2)) * MDP(11) + (t368 * t382 - t369 * t383 + (t339 * t492 + t338 * t494 + (-t455 * t585 - t624) * t493 + t638) * pkin(2)) * MDP(12) + (t355 * t497 - t508 * t549) * MDP(13) + ((t355 - t615) * t501 + (-t356 + t613) * t497) * MDP(14) + (t497 * t409 + t442 * t549 + t612) * MDP(15) + (-t442 * t580 + t589 + t614) * MDP(16) - t442 * t449 * MDP(17) + (-t340 * t449 + t487 * t356 - t382 * t421 + (t383 * t442 + t517) * t497 + ((-t392 - t581) * t442 + t511) * t501) * MDP(18) + (t487 * t355 + t590 * t442 + t341 * t449 + t382 * t508 + t517 * t501 + (t442 * t581 - t511) * t497) * MDP(19) + (t329 * t497 * t500 + t375 * t509) * MDP(20) + (t373 * t398 + t375 * t397 + (-t373 * t500 - t375 * t496) * t579 + (-t622 - t330 * t500 + (t373 * t496 - t375 * t500) * qJD(5)) * t497) * MDP(21) + (-t329 * t501 + t367 + (t595 - t619) * t497 + t509 * t420) * MDP(22) + (t330 * t501 + (-t496 * t579 + t397) * t420 + (t524 - t635) * t497) * MDP(23) + (-t354 * t501 + t420 * t610) * MDP(24) + (-t332 * t397 - t342 * t373 - t630 * t500 + t505 * t496 + (-t486 * t599 - t320 + (t332 * t496 + t373 * t486) * qJD(4) - t506 * t500) * t501 + (t332 * t577 + t322 * t496 + t535 * t446 + t486 * t330 + (t496 * t616 - t535) * qJD(4)) * t497) * MDP(25) + (-t332 * t398 - t342 * t375 + t630 * t496 + t505 * t500 + (-t486 * t595 + t319 + (t332 * t500 + t375 * t486) * qJD(4) + t506 * t496) * t501 + (-t332 * t578 + t322 * t500 + t326 * t446 + t486 * t329 + (t486 * t548 - t326) * qJD(4)) * t497) * MDP(26); (-t446 ^ 2 - t449 ^ 2) * MDP(11) + (t589 - t614) * MDP(18) + MDP(19) * t612 + t397 * t420 * MDP(25) + (t398 * t420 + t367) * MDP(26) + ((-t496 * t582 - t330) * MDP(25) + (-t500 * t582 - t329) * MDP(26) - t442 ^ 2 * MDP(19)) * t501 + (-qJD(4) * t442 * MDP(18) - t409 * MDP(19) + (t524 + t635) * MDP(25) + (t420 * t578 - t595 - t619) * MDP(26)) * t497 + (-g(3) * t495 + t368 * t449 - t369 * t446 + t572 + (-t537 - t538) * t493) * MDP(12); -t421 ^ 2 * MDP(14) + (t355 + t615) * MDP(15) + (-t356 - t613) * MDP(16) + t409 * MDP(17) + (t341 * t442 - t516 - t629) * MDP(18) + (g(1) * t403 + g(2) * t400 + g(3) * t432 + t340 * t442 + t361 * t421 - t521) * MDP(19) + (t375 * t548 + t622) * MDP(20) + ((t329 - t621) * t500 + (-t330 - t620) * t496) * MDP(21) + (t420 * t548 + t599) * MDP(22) + (-t420 ^ 2 * t496 + t595) * MDP(23) + (-pkin(4) * t330 - t341 * t373 + t507 * t496 - t500 * t631) * MDP(25) + (-pkin(4) * t329 - t341 * t375 + t496 * t631 + t507 * t500) * MDP(26) - (MDP(13) * t421 - MDP(14) * t508 - t361 * MDP(18) - MDP(22) * t375 + t373 * MDP(23) - t420 * MDP(24) + MDP(25) * t535 + t326 * MDP(26)) * t508; t375 * t373 * MDP(20) + (-t373 ^ 2 + t375 ^ 2) * MDP(21) + (t566 + t621) * MDP(22) + (-t551 + t620) * MDP(23) + t354 * MDP(24) + (-g(1) * t357 + g(2) * t637 + g(3) * t394 + t326 * t420 - t332 * t375 + t553) * MDP(25) + (g(1) * t358 + g(2) * t636 + g(3) * t395 + t332 * t373 - t535 * t420 - t536) * MDP(26) + (MDP(22) * t611 - MDP(23) * t375 - MDP(25) * t326 + MDP(26) * t535) * qJD(5);];
tau = t1;
