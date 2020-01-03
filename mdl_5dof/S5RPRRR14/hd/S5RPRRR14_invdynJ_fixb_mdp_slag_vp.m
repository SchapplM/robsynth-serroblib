% Calculate vector of inverse dynamics joint torques for
% S5RPRRR14
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR14_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR14_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:20:02
% EndTime: 2019-12-31 19:20:16
% DurationCPUTime: 9.88s
% Computational Cost: add. (8222->548), mult. (26366->774), div. (0->0), fcn. (23167->14), ass. (0->240)
t545 = sin(pkin(5));
t544 = sin(pkin(6));
t670 = cos(pkin(5));
t618 = t670 * t544;
t673 = cos(qJ(3));
t582 = t673 * t618;
t547 = cos(pkin(6));
t546 = cos(pkin(11));
t629 = t673 * t546;
t605 = t547 * t629;
t691 = t545 * t605 + t582;
t543 = sin(pkin(11));
t551 = sin(qJ(1));
t554 = cos(qJ(1));
t620 = t554 * t670;
t503 = t551 * t543 - t546 * t620;
t504 = t543 * t620 + t551 * t546;
t550 = sin(qJ(3));
t655 = t547 * t550;
t656 = t545 * t554;
t659 = t544 * t550;
t432 = -t503 * t655 + t504 * t673 - t656 * t659;
t478 = -t503 * t544 + t547 * t656;
t549 = sin(qJ(4));
t553 = cos(qJ(4));
t406 = t432 * t553 - t478 * t549;
t632 = t544 * t673;
t607 = t545 * t632;
t631 = t547 * t673;
t431 = t503 * t631 + t504 * t550 + t554 * t607;
t548 = sin(qJ(5));
t552 = cos(qJ(5));
t690 = t406 * t548 - t431 * t552;
t689 = t406 * t552 + t431 * t548;
t565 = t545 * (-t543 * t655 + t629);
t491 = qJD(1) * t565;
t625 = qJD(3) * t673;
t600 = t544 * t625;
t688 = t491 - t600;
t598 = t550 * t618;
t630 = t673 * t543;
t476 = t545 * (t546 * t655 + t630) + t598;
t468 = t476 * qJD(3);
t637 = qJDD(1) * t545;
t623 = t543 * t637;
t421 = qJD(1) * t468 - qJDD(1) * t691 + t550 * t623;
t648 = qJD(1) * t545;
t628 = t543 * t648;
t684 = t691 * qJD(1) - t550 * t628;
t459 = qJD(4) - t684;
t469 = t476 * qJD(1);
t615 = qJD(1) * t670;
t627 = t546 * t648;
t639 = -t544 * t627 + qJD(3);
t570 = -t547 * t615 - t639;
t494 = t553 * t570;
t425 = t469 * t549 + t494;
t424 = qJD(5) + t425;
t683 = t432 * t549 + t478 * t553;
t541 = t545 ^ 2;
t682 = t541 * (t543 ^ 2 + t546 ^ 2);
t633 = pkin(1) * t670;
t537 = t546 * t633;
t661 = t543 * t545;
t559 = t670 * pkin(2) + (-pkin(8) * t547 - qJ(2)) * t661;
t477 = t537 + t559;
t495 = (-pkin(8) * t543 * t544 - pkin(2) * t546 - pkin(1)) * t545;
t428 = -t477 * t544 + t547 * t495;
t660 = t543 * t550;
t475 = t545 * t660 - t691;
t391 = pkin(3) * t475 - pkin(9) * t476 + t428;
t617 = t670 * t547;
t658 = t545 * t546;
t502 = t544 * t658 - t617;
t507 = qJ(2) * t658 + t543 * t633;
t563 = (t547 * t658 + t618) * pkin(8);
t472 = t563 + t507;
t458 = t673 * t472;
t634 = t477 * t655 + t495 * t659 + t458;
t395 = -pkin(9) * t502 + t634;
t681 = t549 * t391 + t553 * t395;
t508 = -t547 * t553 + t549 * t659;
t604 = t544 * t628;
t680 = qJD(4) * t508 + t549 * t604 + t688 * t553;
t557 = -t550 * t472 + t477 * t631 + t495 * t632;
t509 = t547 * t549 + t553 * t659;
t679 = qJD(4) * t509 - t688 * t549 + t553 * t604;
t564 = t545 * (t546 * t550 + t547 * t630);
t490 = qJD(1) * t564;
t647 = qJD(3) * t550;
t626 = t544 * t647;
t678 = t490 - t626;
t419 = qJDD(4) + t421;
t601 = pkin(1) * t615;
t532 = t546 * t601;
t460 = qJD(1) * t559 + t532;
t485 = qJD(1) * t495 + qJD(2);
t423 = -t460 * t544 + t547 * t485;
t382 = -pkin(3) * t684 - pkin(9) * t469 + t423;
t446 = t460 * t655;
t501 = qJ(2) * t627 + t543 * t601;
t453 = qJD(1) * t563 + t501;
t397 = t673 * t453 + t485 * t659 + t446;
t384 = -pkin(9) * t570 + t397;
t359 = t382 * t549 + t384 * t553;
t610 = qJDD(1) * t670;
t636 = qJDD(1) * t546;
t622 = t545 * t636;
t498 = t544 * t622 - t547 * t610 - qJDD(3);
t599 = pkin(1) * t610;
t638 = qJD(1) * qJD(2);
t624 = t545 * t638;
t487 = qJ(2) * t622 + t543 * t599 + t546 * t624;
t441 = qJDD(1) * t563 + t487;
t530 = t546 * t599;
t442 = qJDD(1) * t559 - t543 * t624 + t530;
t481 = qJDD(1) * t495 + qJDD(2);
t606 = t460 * t631;
t566 = -qJD(3) * t606 - t673 * t441 - t442 * t655 + t453 * t647 - t481 * t659 - t485 * t600;
t363 = -pkin(9) * t498 - t566;
t412 = -t442 * t544 + t547 * t481;
t420 = t684 * qJD(3) + qJDD(1) * t598 + t622 * t655 + t673 * t623;
t371 = pkin(3) * t421 - pkin(9) * t420 + t412;
t613 = t363 * t549 - t553 * t371;
t674 = t359 * qJD(4) + t613;
t346 = -pkin(4) * t419 + t674;
t427 = t553 * t469 - t549 * t570;
t621 = t551 * t670;
t505 = -t543 * t621 + t554 * t546;
t569 = t554 * t543 + t546 * t621;
t562 = t569 * t547;
t657 = t545 * t551;
t436 = t505 * t673 + (t544 * t657 - t562) * t550;
t480 = t544 * t569 + t547 * t657;
t408 = -t436 * t549 + t480 * t553;
t429 = t476 * t549 + t502 * t553;
t577 = g(1) * t408 - g(2) * t683 - g(3) * t429;
t676 = t424 * (pkin(4) * t427 + t424 * pkin(10)) + t346 + t577;
t611 = t549 * t420 + t553 * t498;
t378 = qJD(4) * t427 + t611;
t376 = qJDD(5) + t378;
t524 = -pkin(4) * t553 - pkin(10) * t549 - pkin(3);
t675 = t424 * (t397 - t459 * (pkin(4) * t549 - pkin(10) * t553)) - t524 * t376;
t672 = pkin(1) * t541;
t671 = pkin(9) * qJD(4);
t645 = qJD(4) * t549;
t377 = -qJD(4) * t494 + t553 * t420 - t469 * t645 - t549 * t498;
t641 = qJD(5) * t552;
t635 = t552 * t377 + t548 * t419 + t459 * t641;
t642 = qJD(5) * t548;
t353 = -t427 * t642 + t635;
t669 = t353 * t548;
t664 = t427 * t548;
t398 = -t552 * t459 + t664;
t668 = t398 * t424;
t400 = t427 * t552 + t459 * t548;
t667 = t400 * t424;
t666 = t425 * t459;
t665 = t427 * t459;
t663 = t684 * t553;
t654 = t548 * t376;
t653 = t552 * t376;
t396 = -t550 * t453 + t485 * t632 + t606;
t422 = pkin(3) * t469 - pkin(9) * t684;
t650 = t553 * t396 + t549 * t422;
t646 = qJD(4) * t548;
t644 = qJD(4) * t552;
t643 = qJD(4) * t553;
t640 = t498 * MDP(12);
t555 = qJD(1) ^ 2;
t619 = t555 * t670;
t579 = t553 * t363 + t549 * t371 + t382 * t643 - t384 * t645;
t345 = pkin(10) * t419 + t579;
t581 = -qJD(3) * t446 - t550 * t441 + t442 * t631 - t453 * t625 + t481 * t632 - t485 * t626;
t364 = pkin(3) * t498 - t581;
t348 = pkin(4) * t378 - pkin(10) * t377 + t364;
t614 = -t548 * t345 + t552 * t348;
t612 = t377 * t548 - t552 * t419;
t609 = t459 * t553;
t608 = t552 * t424;
t602 = qJD(2) * t544 * t661;
t596 = g(1) * t554 + g(2) * t551;
t595 = g(1) * t551 - g(2) * t554;
t594 = qJD(2) * t615;
t414 = t469 * t548 + t552 * t663;
t592 = t552 * t643 - t414;
t590 = t552 * t345 + t548 * t348;
t356 = pkin(10) * t459 + t359;
t383 = pkin(3) * t570 - t396;
t365 = t425 * pkin(4) - t427 * pkin(10) + t383;
t350 = t356 * t552 + t365 * t548;
t589 = t356 * t548 - t365 * t552;
t362 = pkin(10) * t475 + t681;
t394 = t502 * pkin(3) - t557;
t430 = t476 * t553 - t502 * t549;
t370 = t429 * pkin(4) - t430 * pkin(10) + t394;
t588 = t362 * t552 + t370 * t548;
t587 = -t362 * t548 + t370 * t552;
t358 = t382 * t553 - t384 * t549;
t386 = qJD(2) * t565 + qJD(3) * t557;
t467 = (t582 + (t605 - t660) * t545) * qJD(3);
t411 = pkin(3) * t468 - pkin(9) * t467 + t602;
t586 = -t386 * t549 + t411 * t553;
t585 = t391 * t553 - t395 * t549;
t404 = t430 * t552 + t475 * t548;
t403 = t430 * t548 - t475 * t552;
t583 = (-qJ(2) * t628 + t532) * t543 - t501 * t546;
t580 = -pkin(9) * t419 + t383 * t459;
t578 = t553 * t386 + t391 * t643 - t395 * t645 + t549 * t411;
t435 = t505 * t550 - t551 * t607 + t562 * t673;
t576 = g(1) * t435 + g(2) * t431 + g(3) * t475;
t575 = -g(1) * t436 - g(2) * t432 - g(3) * t476;
t573 = -t548 * t509 - t552 * t632;
t571 = -t552 * t509 + t548 * t632;
t567 = -t364 + t576;
t355 = -pkin(4) * t459 - t358;
t561 = -pkin(10) * t376 + (t355 + t358) * t424;
t558 = pkin(9) * qJD(5) * t424 - t576;
t556 = (pkin(10) * t469 - qJD(5) * t524 + t650) * t424 + t575;
t387 = qJD(2) * t564 + (t458 + (t477 * t547 + t495 * t544) * t550) * qJD(3);
t533 = -pkin(1) * t637 + qJDD(2);
t506 = -qJ(2) * t661 + t537;
t486 = t530 + (-qJ(2) * qJDD(1) - t638) * t661;
t413 = -t552 * t469 + t548 * t663;
t409 = t436 * t553 + t480 * t549;
t402 = -qJD(4) * t429 + t467 * t553;
t401 = qJD(4) * t430 + t467 * t549;
t380 = t409 * t552 + t435 * t548;
t379 = -t409 * t548 + t435 * t552;
t373 = -qJD(5) * t403 + t402 * t552 + t468 * t548;
t372 = qJD(5) * t404 + t402 * t548 - t468 * t552;
t366 = -pkin(4) * t469 + t396 * t549 - t422 * t553;
t361 = -pkin(4) * t475 - t585;
t357 = t401 * pkin(4) - t402 * pkin(10) + t387;
t354 = qJD(5) * t400 + t612;
t352 = -pkin(4) * t468 + qJD(4) * t681 - t586;
t351 = pkin(10) * t468 + t578;
t344 = -t350 * qJD(5) + t614;
t343 = -t589 * qJD(5) + t590;
t1 = [(t586 * t459 + t585 * t419 - t613 * t475 + t358 * t468 + t387 * t425 + t394 * t378 + t364 * t429 + t383 * t401 + g(1) * t406 - g(2) * t409 + (-t359 * t475 - t459 * t681) * qJD(4)) * MDP(20) + (-t420 * t475 - t421 * t476 + t467 * t684 - t468 * t469) * MDP(9) + (g(1) * t432 - g(2) * t436 + t387 * t570 + t412 * t475 + t428 * t421 + t423 * t468 - t498 * t557 - t502 * t581 - t602 * t684) * MDP(13) + ((-qJD(5) * t588 - t351 * t548 + t357 * t552) * t424 + t587 * t376 + t344 * t429 - t589 * t401 + t352 * t398 + t361 * t354 + t346 * t403 + t355 * t372 + g(1) * t689 - g(2) * t380) * MDP(27) + (-(qJD(5) * t587 + t351 * t552 + t357 * t548) * t424 - t588 * t376 - t343 * t429 - t350 * t401 + t352 * t400 + t361 * t353 + t346 * t404 + t355 * t373 - g(1) * t690 - g(2) * t379) * MDP(28) + (-t487 * t670 - g(1) * t503 + g(2) * t569 + (t533 * t543 - t546 * t594) * t545 + (-t507 * t670 - t543 * t672) * qJDD(1)) * MDP(5) + (t486 * t670 + g(1) * t504 - g(2) * t505 + (-t533 * t546 - t543 * t594) * t545 + (t506 * t670 + t546 * t672) * qJDD(1)) * MDP(4) + (-g(1) * t683 - g(2) * t408 - t359 * t468 + t364 * t430 + t394 * t377 + t383 * t402 + t387 * t427 - t681 * t419 - t578 * t459 - t579 * t475) * MDP(21) + t595 * MDP(2) + t596 * MDP(3) + (t486 * t506 + t487 * t507 + t595 * pkin(1) + (-t533 * pkin(1) - qJ(2) * t596 - qJD(2) * t583) * t545) * MDP(7) + t502 * t640 + (t420 * t476 + t467 * t469) * MDP(8) + (-t378 * t475 - t401 * t459 - t419 * t429 - t425 * t468) * MDP(18) + (t377 * t475 + t402 * t459 + t419 * t430 + t427 * t468) * MDP(17) + (t419 * t475 + t459 * t468) * MDP(19) + (-t354 * t429 - t372 * t424 - t376 * t403 - t398 * t401) * MDP(25) + (t353 * t429 + t373 * t424 + t376 * t404 + t400 * t401) * MDP(24) + (t376 * t429 + t401 * t424) * MDP(26) + (-t377 * t429 - t378 * t430 - t401 * t427 - t402 * t425) * MDP(16) + (t377 * t430 + t402 * t427) * MDP(15) + (-t353 * t403 - t354 * t404 - t372 * t400 - t373 * t398) * MDP(23) + (t353 * t404 + t373 * t400) * MDP(22) + (t638 * t682 + (-t486 * t543 + t487 * t546 + (-t506 * t543 + t507 * t546) * qJDD(1) - t596) * t545) * MDP(6) + qJDD(1) * MDP(1) + (-g(1) * t431 + g(2) * t435 + t386 * t570 + t412 * t476 + t428 * t420 + t423 * t467 + t469 * t602 + t498 * t634 - t502 * t566) * MDP(14) + (-t420 * t502 - t467 * t570 - t476 * t498) * MDP(10) + (t421 * t502 + t468 * t570 + t475 * t498) * MDP(11); -t555 * MDP(6) * t682 + (-g(3) * t670 + qJDD(2)) * MDP(7) + (t547 * t421 - t490 * t570 + (-t498 * t673 + t570 * t647 + t628 * t684) * t544) * MDP(13) + (t547 * t420 - t491 * t570 + (-t469 * t628 + t550 * t498 + t570 * t625) * t544) * MDP(14) + (-t378 * t632 - t508 * t419 - t425 * t678 - t459 * t679) * MDP(20) + (-t377 * t632 - t509 * t419 - t427 * t678 + t680 * t459) * MDP(21) + (t508 * t354 + t573 * t376 + (qJD(5) * t571 + t548 * t680 - t552 * t678) * t424 + t679 * t398) * MDP(27) + (t508 * t353 + t571 * t376 + (-qJD(5) * t573 + t548 * t678 + t552 * t680) * t424 + t679 * t400) * MDP(28) + ((t543 * t619 - t636) * MDP(4) + (qJDD(1) * t543 + t546 * t619) * MDP(5) + (-pkin(1) * qJDD(1) + qJD(1) * t583 - t595) * MDP(7)) * t545; -t684 ^ 2 * MDP(9) + (t570 * t684 + t420) * MDP(10) - t421 * MDP(11) - t640 + (-t397 * t570 + t576 + t581) * MDP(13) + (-t396 * t570 - t423 * t684 + t566 - t575) * MDP(14) + (t377 * t549 + t427 * t609) * MDP(15) + ((t377 - t666) * t553 + (-t378 - t665) * t549) * MDP(16) + (t549 * t419 + t459 * t609) * MDP(17) + (-t459 ^ 2 * t549 + t553 * t419) * MDP(18) + (-pkin(3) * t378 - t397 * t425 + (t396 * t459 + t580) * t549 + ((-t422 - t671) * t459 + t567) * t553) * MDP(20) + (-pkin(3) * t377 + t650 * t459 - t397 * t427 + t580 * t553 + (t459 * t671 - t567) * t549) * MDP(21) + (t353 * t549 * t552 + (-t549 * t642 + t592) * t400) * MDP(22) + (t398 * t414 + t400 * t413 + (-t398 * t552 - t400 * t548) * t643 + (-t669 - t354 * t552 + (t398 * t548 - t400 * t552) * qJD(5)) * t549) * MDP(23) + (-t353 * t553 + t592 * t424 + (t400 * t459 - t424 * t642 + t653) * t549) * MDP(24) + (t354 * t553 + (-t548 * t643 + t413) * t424 + (-t398 * t459 - t424 * t641 - t654) * t549) * MDP(25) + (t424 * t459 * t549 - t376 * t553) * MDP(26) + (-t355 * t413 - t366 * t398 - t675 * t552 + t556 * t548 + (t355 * t646 - t344 + (qJD(4) * t398 - t654) * pkin(9) - t558 * t552) * t553 + (t355 * t641 + t346 * t548 - t459 * t589 + (t424 * t646 + t354) * pkin(9)) * t549) * MDP(27) + (-t355 * t414 - t366 * t400 + t675 * t548 + t556 * t552 + (t355 * t644 + t343 + (qJD(4) * t400 - t653) * pkin(9) + t558 * t548) * t553 + (-t355 * t642 + t346 * t552 - t459 * t350 + (t424 * t644 + t353) * pkin(9)) * t549) * MDP(28) + (-t684 * MDP(8) + (qJD(1) * t617 + t639) * MDP(11) - t423 * MDP(13) - t427 * MDP(17) + t425 * MDP(18) - t459 * MDP(19) - t358 * MDP(20) + t359 * MDP(21) + MDP(9) * t469) * t469; -t425 ^ 2 * MDP(16) + (t377 + t666) * MDP(17) + (-t611 + t665) * MDP(18) + t419 * MDP(19) + (t359 * t459 - t577 - t674) * MDP(20) + (g(1) * t409 + g(2) * t406 + g(3) * t430 + t358 * t459 + t383 * t425 - t579) * MDP(21) + (t400 * t608 + t669) * MDP(22) + ((t353 - t668) * t552 + (-t354 - t667) * t548) * MDP(23) + (t424 * t608 + t654) * MDP(24) + (-t424 ^ 2 * t548 + t653) * MDP(25) + (-pkin(4) * t354 - t359 * t398 + t561 * t548 - t552 * t676) * MDP(27) + (-pkin(4) * t353 - t359 * t400 + t548 * t676 + t561 * t552) * MDP(28) + (MDP(15) * t425 + MDP(16) * t427 - MDP(18) * qJD(4) - MDP(20) * t383 - MDP(24) * t400 + MDP(25) * t398 - MDP(26) * t424 + MDP(27) * t589 + MDP(28) * t350) * t427; t400 * t398 * MDP(22) + (-t398 ^ 2 + t400 ^ 2) * MDP(23) + (t635 + t668) * MDP(24) + (-t612 + t667) * MDP(25) + t376 * MDP(26) + (-g(1) * t379 + g(2) * t690 + g(3) * t403 + t350 * t424 - t355 * t400 + t614) * MDP(27) + (g(1) * t380 + g(2) * t689 + g(3) * t404 + t355 * t398 - t589 * t424 - t590) * MDP(28) + (-MDP(24) * t664 - MDP(25) * t400 - MDP(27) * t350 + MDP(28) * t589) * qJD(5);];
tau = t1;
