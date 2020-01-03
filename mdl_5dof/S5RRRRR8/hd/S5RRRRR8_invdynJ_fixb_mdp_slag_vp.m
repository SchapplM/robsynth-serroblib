% Calculate vector of inverse dynamics joint torques for
% S5RRRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:26:04
% EndTime: 2019-12-31 22:26:16
% DurationCPUTime: 7.98s
% Computational Cost: add. (5448->490), mult. (12019->643), div. (0->0), fcn. (9010->14), ass. (0->226)
t548 = qJD(2) + qJD(3);
t555 = sin(qJ(3));
t560 = cos(qJ(2));
t671 = cos(qJ(3));
t616 = qJD(1) * t671;
t556 = sin(qJ(2));
t636 = qJD(1) * t556;
t688 = -t555 * t636 + t560 * t616;
t692 = t688 * t548;
t552 = qJ(2) + qJ(3);
t545 = cos(t552);
t668 = g(3) * t545;
t543 = sin(t552);
t557 = sin(qJ(1));
t561 = cos(qJ(1));
t598 = g(1) * t561 + g(2) * t557;
t690 = t598 * t543;
t691 = -t690 + t668;
t674 = qJD(4) + qJD(5);
t689 = -t688 + t674;
t673 = pkin(6) + pkin(7);
t515 = t673 * t560;
t500 = qJD(1) * t515;
t487 = t555 * t500;
t513 = t673 * t556;
t498 = qJD(1) * t513;
t446 = -t671 * t498 - t487;
t615 = qJD(3) * t671;
t681 = pkin(2) * t615 - t446;
t646 = t555 * t560;
t486 = -qJD(1) * t646 - t556 * t616;
t439 = -pkin(3) * t486 - pkin(8) * t688;
t428 = pkin(2) * t636 + t439;
t554 = sin(qJ(4));
t559 = cos(qJ(4));
t687 = -t559 * t428 - t554 * t681;
t458 = -t486 * t554 - t559 * t548;
t558 = cos(qJ(5));
t553 = sin(qJ(5));
t586 = t486 * t559 - t548 * t554;
t656 = t586 * t553;
t403 = t558 * t458 - t656;
t480 = qJD(4) - t688;
t473 = qJD(5) + t480;
t686 = t403 * t473;
t587 = t458 * t553 + t558 * t586;
t685 = t473 * t587;
t547 = qJDD(2) + qJDD(3);
t629 = qJD(1) * qJD(2);
t613 = t560 * t629;
t628 = qJDD(1) * t556;
t455 = qJDD(2) * pkin(2) + t673 * (-t613 - t628);
t614 = t556 * t629;
t627 = qJDD(1) * t560;
t457 = t673 * (-t614 + t627);
t665 = qJD(2) * pkin(2);
t489 = -t498 + t665;
t635 = qJD(3) * t555;
t602 = -t671 * t455 + t555 * t457 + t489 * t635 + t500 * t615;
t381 = -pkin(3) * t547 + t602;
t684 = t381 + t668;
t650 = t553 * t559;
t496 = t554 * t558 + t650;
t683 = t689 * t496;
t494 = t553 * t554 - t558 * t559;
t642 = t689 * t494;
t541 = -pkin(2) * t560 - pkin(1);
t511 = t541 * qJD(1);
t426 = -pkin(3) * t688 + pkin(8) * t486 + t511;
t488 = t671 * t500;
t443 = t555 * t489 + t488;
t430 = pkin(8) * t548 + t443;
t393 = t426 * t554 + t430 * t559;
t374 = -pkin(9) * t458 + t393;
t632 = qJD(5) * t553;
t371 = t374 * t632;
t442 = t671 * t489 - t487;
t429 = -t548 * pkin(3) - t442;
t400 = t458 * pkin(4) + t429;
t551 = qJ(4) + qJ(5);
t542 = sin(t551);
t544 = cos(t551);
t652 = t545 * t557;
t466 = t542 * t561 - t544 * t652;
t651 = t545 * t561;
t468 = t542 * t557 + t544 * t651;
t534 = g(3) * t543;
t679 = g(1) * t468 - g(2) * t466 + t400 * t403 + t544 * t534 + t371;
t465 = t542 * t652 + t544 * t561;
t467 = -t542 * t651 + t544 * t557;
t610 = qJDD(1) * t671;
t414 = t555 * t627 + t556 * t610 + t692;
t497 = t671 * t556 + t646;
t451 = t548 * t497;
t590 = t555 * t628 - t560 * t610;
t415 = t451 * qJD(1) + t590;
t481 = pkin(2) * t614 + t541 * qJDD(1);
t372 = pkin(3) * t415 - pkin(8) * t414 + t481;
t370 = t559 * t372;
t568 = t555 * t455 + t671 * t457 + t489 * t615 - t500 * t635;
t380 = t547 * pkin(8) + t568;
t633 = qJD(4) * t559;
t634 = qJD(4) * t554;
t390 = t559 * t414 + t486 * t634 + t554 * t547 + t548 * t633;
t413 = qJDD(4) + t415;
t348 = pkin(4) * t413 - pkin(9) * t390 - t393 * qJD(4) - t380 * t554 + t370;
t391 = -t586 * qJD(4) + t414 * t554 - t559 * t547;
t574 = t554 * t372 + t559 * t380 + t426 * t633 - t430 * t634;
t350 = -pkin(9) * t391 + t574;
t609 = t558 * t348 - t553 * t350;
t678 = -g(1) * t467 + g(2) * t465 + t400 * t587 + t542 * t534 + t609;
t411 = qJDD(5) + t413;
t677 = t411 * MDP(29) + (-t403 ^ 2 + t587 ^ 2) * MDP(26) - t403 * t587 * MDP(25);
t434 = t496 * t497;
t445 = -t555 * t498 + t488;
t600 = pkin(2) * t635 - t445;
t655 = t688 * t554;
t676 = (t634 - t655) * pkin(4);
t675 = t554 * t428 - t681 * t559;
t608 = t390 * t553 + t558 * t391;
t359 = -t587 * qJD(5) + t608;
t672 = -pkin(8) - pkin(9);
t667 = t559 * pkin(4);
t546 = t559 * pkin(9);
t538 = pkin(2) * t555 + pkin(8);
t666 = -pkin(9) - t538;
t392 = t559 * t426 - t430 * t554;
t373 = pkin(9) * t586 + t392;
t368 = pkin(4) * t480 + t373;
t664 = t368 * t558;
t663 = t374 * t558;
t662 = t390 * t554;
t661 = t429 * t688;
t580 = -t555 * t556 + t671 * t560;
t450 = t548 * t580;
t660 = t450 * t554;
t659 = t450 * t559;
t658 = t458 * t480;
t657 = t586 * t480;
t654 = t497 * t554;
t653 = t497 * t559;
t649 = t554 * t413;
t648 = t554 * t557;
t647 = t554 * t561;
t645 = t557 * t559;
t462 = -t555 * t513 + t671 * t515;
t456 = t559 * t462;
t644 = t559 * t561;
t641 = t554 * t439 + t559 * t442;
t441 = -pkin(3) * t580 - pkin(8) * t497 + t541;
t639 = t554 * t441 + t456;
t638 = t676 + t600;
t549 = t556 ^ 2;
t637 = -t560 ^ 2 + t549;
t631 = qJD(5) * t558;
t626 = pkin(9) * t655;
t625 = t556 * t665;
t622 = qJD(4) * pkin(8) * t480;
t621 = t558 * t390 - t553 * t391 - t458 * t631;
t620 = qJD(2) * t673;
t619 = qJD(4) * t672;
t617 = t497 * t634;
t421 = t429 * t633;
t611 = qJD(4) * t666;
t606 = t480 * t559;
t605 = -qJD(4) * t426 - t380;
t604 = qJD(5) * t368 + t350;
t539 = -t671 * pkin(2) - pkin(3);
t601 = -t486 * pkin(4) - t546 * t688;
t599 = -t443 + t676;
t597 = g(1) * t557 - g(2) * t561;
t596 = -t430 * t633 + t370;
t491 = t538 * t559 + t546;
t595 = qJD(5) * t491 - t559 * t611 + t601 - t687;
t433 = t559 * t439;
t514 = pkin(8) * t559 + t546;
t594 = qJD(5) * t514 - t442 * t554 - t559 * t619 + t433 + t601;
t490 = t666 * t554;
t593 = -qJD(5) * t490 - t554 * t611 - t626 + t675;
t512 = t672 * t554;
t592 = -qJD(5) * t512 - t554 * t619 - t626 + t641;
t591 = -pkin(8) * t413 - t661;
t354 = t368 * t553 + t663;
t588 = -t413 * t538 - t661;
t585 = -t393 * t486 + t554 * t684 + t421;
t584 = t392 * t486 + t429 * t634 + t559 * t690;
t582 = -0.2e1 * pkin(1) * t629 - pkin(6) * qJDD(2);
t581 = -t671 * t513 - t555 * t515;
t579 = t497 * t633 + t660;
t578 = -t617 + t659;
t399 = pkin(3) * t451 - pkin(8) * t450 + t625;
t499 = t556 * t620;
t501 = t560 * t620;
t407 = t581 * qJD(3) - t671 * t499 - t555 * t501;
t573 = t554 * t399 + t559 * t407 + t441 * t633 - t462 * t634;
t358 = t586 * t632 + t621;
t562 = qJD(2) ^ 2;
t571 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t562 + t597;
t563 = qJD(1) ^ 2;
t570 = pkin(1) * t563 - pkin(6) * qJDD(1) + t598;
t569 = t511 * t486 - t602 - t691;
t353 = -t374 * t553 + t664;
t362 = pkin(4) * t391 + t381;
t567 = t353 * t486 + t362 * t494 + t400 * t683 - t544 * t691;
t408 = t462 * qJD(3) - t555 * t499 + t671 * t501;
t566 = -t354 * t486 + t362 * t496 - t642 * t400 + t691 * t542;
t565 = g(1) * t651 + g(2) * t652 - t511 * t688 + t534 - t568;
t564 = (-t358 * t494 - t359 * t496 + t642 * t403 + t587 * t683) * MDP(26) + (t358 * t496 + t587 * t642) * MDP(25) + ((t390 - t658) * t559 + (-t391 + t657) * t554) * MDP(19) + (t411 * t496 - t642 * t473 - t486 * t587) * MDP(27) + (-t403 * t486 - t411 * t494 - t473 * t683) * MDP(28) + (-t586 * t606 + t662) * MDP(18) + (-t480 ^ 2 * t554 + t413 * t559 - t458 * t486) * MDP(21) + (t480 * t606 - t486 * t586 + t649) * MDP(20) + (t414 - t692) * MDP(13) + (-t590 + (-qJD(1) * t497 - t486) * t548) * MDP(14) + (t486 ^ 2 - t688 ^ 2) * MDP(12) + t547 * MDP(15) + (MDP(11) * t688 + t480 * MDP(22) + t473 * MDP(29)) * t486;
t540 = -pkin(3) - t667;
t510 = t539 - t667;
t478 = t545 * t644 + t648;
t477 = -t545 * t647 + t645;
t476 = -t545 * t645 + t647;
t475 = t545 * t648 + t644;
t437 = t559 * t441;
t435 = t494 * t497;
t427 = pkin(4) * t654 - t581;
t396 = t559 * t399;
t394 = -pkin(9) * t654 + t639;
t383 = -pkin(4) * t580 - pkin(9) * t653 - t462 * t554 + t437;
t377 = t579 * pkin(4) + t408;
t367 = t450 * t650 - t553 * t617 - t632 * t654 + (t674 * t653 + t660) * t558;
t366 = -t674 * t434 - t494 * t450;
t355 = -t579 * pkin(9) + t573;
t352 = -pkin(9) * t659 + pkin(4) * t451 - t407 * t554 + t396 + (-t456 + (pkin(9) * t497 - t441) * t554) * qJD(4);
t1 = [(-t408 * t548 + t415 * t541 + t451 * t511 - t481 * t580 + t597 * t545 + t547 * t581 - t625 * t688) * MDP(16) + (t414 * t580 - t415 * t497 + t450 * t688 + t451 * t486) * MDP(12) + (qJDD(1) * t549 + 0.2e1 * t556 * t613) * MDP(4) + (-t407 * t548 + t414 * t541 + t450 * t511 - t462 * t547 + t481 * t497 - t486 * t625 - t597 * t543) * MDP(17) + 0.2e1 * (t556 * t627 - t637 * t629) * MDP(5) + ((t352 * t558 - t355 * t553) * t473 + (t383 * t558 - t394 * t553) * t411 - t609 * t580 + t353 * t451 + t377 * t403 + t427 * t359 + t362 * t434 + t400 * t367 - g(1) * t466 - g(2) * t468 + ((-t383 * t553 - t394 * t558) * t473 + t354 * t580) * qJD(5)) * MDP(30) + (t391 * t580 - t451 * t458 - t579 * t480 - t497 * t649) * MDP(21) + (-t358 * t434 + t359 * t435 - t366 * t403 + t367 * t587) * MDP(26) + (-t358 * t435 - t366 * t587) * MDP(25) + ((-t458 * t559 + t554 * t586) * t450 + (-t662 - t391 * t559 + (t458 * t554 + t559 * t586) * qJD(4)) * t497) * MDP(19) + (t390 * t653 - t578 * t586) * MDP(18) + (-g(1) * t465 - g(2) * t467 - t354 * t451 + t427 * t358 - t362 * t435 + t400 * t366 - t371 * t580 - t377 * t587 + (-(-qJD(5) * t394 + t352) * t473 - t383 * t411 + t348 * t580) * t553 + (-(qJD(5) * t383 + t355) * t473 - t394 * t411 + t604 * t580) * t558) * MDP(31) + (-t358 * t580 + t366 * t473 - t411 * t435 - t451 * t587) * MDP(27) + (-t390 * t580 + t413 * t653 - t451 * t586 + t578 * t480) * MDP(20) + (-g(1) * t475 - g(2) * t477 + t381 * t653 - t390 * t581 - t393 * t451 - t408 * t586 - t639 * t413 + t578 * t429 - t573 * t480 + t574 * t580) * MDP(24) + (t408 * t458 - t581 * t391 + t497 * t421 + (-t462 * t633 + t396) * t480 + t437 * t413 - t596 * t580 + t392 * t451 - g(1) * t476 - g(2) * t478 + (t381 * t497 + t429 * t450 + (-qJD(4) * t441 - t407) * t480 - t462 * t413 - t605 * t580) * t554) * MDP(23) + (-t451 * t548 + t547 * t580) * MDP(14) + (-t413 * t580 + t451 * t480) * MDP(22) + (t359 * t580 - t367 * t473 - t403 * t451 - t411 * t434) * MDP(28) + (-t411 * t580 + t451 * t473) * MDP(29) + qJDD(1) * MDP(1) + (t582 * t556 + t571 * t560) * MDP(9) + (-t571 * t556 + t582 * t560) * MDP(10) + t597 * MDP(2) + t598 * MDP(3) + (qJDD(2) * t556 + t560 * t562) * MDP(6) + (qJDD(2) * t560 - t556 * t562) * MDP(7) + (t450 * t548 + t497 * t547) * MDP(13) + (t414 * t497 - t450 * t486) * MDP(11); (-g(3) * t560 + t556 * t570) * MDP(9) + (g(3) * t556 + t560 * t570) * MDP(10) + t564 + (t446 * t548 + (t486 * t636 - t547 * t555 - t548 * t615) * pkin(2) + t565) * MDP(17) + (t539 * t390 + t588 * t559 - t554 * t690 - t600 * t586 + (t538 * t634 + t675) * t480 + t585) * MDP(24) + (t445 * t548 + (t671 * t547 - t548 * t635 + t636 * t688) * pkin(2) + t569) * MDP(16) + (t539 * t391 - t684 * t559 + t588 * t554 + t600 * t458 + (-t538 * t633 + t687) * t480 + t584) * MDP(23) + qJDD(2) * MDP(8) + (-(t490 * t553 + t491 * t558) * t411 + t510 * t358 + (t553 * t595 + t558 * t593) * t473 - t638 * t587 + t566) * MDP(31) + ((t490 * t558 - t491 * t553) * t411 + t510 * t359 + (t553 * t593 - t558 * t595) * t473 + t638 * t403 + t567) * MDP(30) + MDP(6) * t628 + MDP(7) * t627 + (-t556 * t560 * MDP(4) + t637 * MDP(5)) * t563; t564 + ((t512 * t558 - t514 * t553) * t411 + t540 * t359 + (t553 * t592 - t558 * t594) * t473 + t599 * t403 + t567) * MDP(30) + (-(t512 * t553 + t514 * t558) * t411 + t540 * t358 + (t553 * t594 + t558 * t592) * t473 - t599 * t587 + t566) * MDP(31) + (-pkin(3) * t390 + t443 * t586 + t641 * t480 + t591 * t559 + (-t690 + t622) * t554 + t585) * MDP(24) + (-pkin(3) * t391 - t433 * t480 - t443 * t458 + (t442 * t480 + t591) * t554 + (-t684 - t622) * t559 + t584) * MDP(23) + (t442 * t548 + t565) * MDP(17) + (t443 * t548 + t569) * MDP(16); -t586 * t458 * MDP(18) + (-t458 ^ 2 + t586 ^ 2) * MDP(19) + (t390 + t658) * MDP(20) + (-t391 - t657) * MDP(21) + t413 * MDP(22) + (-g(1) * t477 + g(2) * t475 + t393 * t480 + t429 * t586 + (t605 + t534) * t554 + t596) * MDP(23) + (g(1) * t478 - g(2) * t476 + t392 * t480 + t429 * t458 + t559 * t534 - t574) * MDP(24) + (t358 + t686) * MDP(27) + (-t359 - t685) * MDP(28) + (-(-t373 * t553 - t663) * t473 - t354 * qJD(5) + (t403 * t586 + t558 * t411 - t473 * t632) * pkin(4) + t678) * MDP(30) + ((-t374 * t473 - t348) * t553 + (t373 * t473 - t604) * t558 + (-t553 * t411 - t473 * t631 - t586 * t587) * pkin(4) + t679) * MDP(31) + t677; (t621 + t686) * MDP(27) + (-t608 - t685) * MDP(28) + (t354 * t473 + t678) * MDP(30) + (-t553 * t348 - t558 * t350 + t353 * t473 + t679) * MDP(31) + (MDP(27) * t656 + MDP(28) * t587 - MDP(30) * t354 - MDP(31) * t664) * qJD(5) + t677;];
tau = t1;
