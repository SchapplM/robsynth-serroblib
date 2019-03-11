% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:08:52
% EndTime: 2019-03-09 18:09:05
% DurationCPUTime: 7.25s
% Computational Cost: add. (10168->506), mult. (26010->671), div. (0->0), fcn. (20298->10), ass. (0->250)
t598 = sin(qJ(3));
t599 = sin(qJ(2));
t602 = cos(qJ(3));
t603 = cos(qJ(2));
t562 = t598 * t599 - t602 * t603;
t546 = t562 * qJD(1);
t674 = qJD(1) * t603;
t675 = qJD(1) * t599;
t548 = -t598 * t674 - t602 * t675;
t595 = sin(pkin(11));
t711 = cos(pkin(11));
t649 = -t711 * t546 + t548 * t595;
t720 = qJD(5) + qJD(6);
t742 = t649 - t720;
t619 = -t595 * t546 - t548 * t711;
t717 = pkin(3) * t548;
t461 = pkin(4) * t619 - pkin(9) * t649 - t717;
t589 = pkin(2) * t675;
t458 = t461 + t589;
t597 = sin(qJ(5));
t601 = cos(qJ(5));
t718 = pkin(7) + pkin(8);
t573 = t718 * t603;
t567 = qJD(1) * t573;
t553 = t602 * t567;
t572 = t718 * t599;
t565 = qJD(1) * t572;
t645 = t565 * t598 - t553;
t710 = qJ(4) * t546;
t502 = t645 + t710;
t541 = t548 * qJ(4);
t549 = t598 * t567;
t678 = -t602 * t565 - t549;
t503 = t541 + t678;
t690 = t595 * t598;
t713 = pkin(2) * qJD(3);
t679 = -t595 * t502 - t503 * t711 + (t602 * t711 - t690) * t713;
t741 = -t601 * t458 - t597 * t679;
t671 = qJD(5) * t597;
t695 = t649 * t597;
t740 = t671 - t695;
t592 = qJD(2) + qJD(3);
t499 = -t601 * t592 + t597 * t619;
t600 = cos(qJ(6));
t501 = t592 * t597 + t601 * t619;
t596 = sin(qJ(6));
t696 = t501 * t596;
t450 = t600 * t499 + t696;
t721 = qJD(5) - t649;
t506 = qJD(6) + t721;
t739 = t450 * t506;
t738 = t499 * t721;
t628 = t499 * t596 - t600 * t501;
t737 = t506 * t628;
t689 = t596 * t601;
t563 = t597 * t600 + t689;
t735 = t742 * t563;
t561 = t596 * t597 - t600 * t601;
t734 = t742 * t561;
t666 = qJD(1) * qJD(2);
t660 = t603 * t666;
t672 = qJD(3) * t602;
t518 = -t592 * t598 * t675 + t602 * t660 + t672 * t674;
t564 = t598 * t603 + t599 * t602;
t722 = qJD(1) * t564;
t608 = t592 * t722;
t476 = t518 * t595 + t608 * t711;
t644 = t721 * t601;
t733 = -t476 * t597 - t721 * t644;
t712 = qJD(2) * pkin(2);
t555 = -t565 + t712;
t648 = t602 * t555 - t549;
t494 = t541 + t648;
t486 = pkin(3) * t592 + t494;
t627 = -t555 * t598 - t553;
t495 = -t627 - t710;
t655 = t711 * t495;
t437 = t595 * t486 + t655;
t435 = pkin(9) * t592 + t437;
t588 = -pkin(2) * t603 - pkin(1);
t571 = t588 * qJD(1);
t526 = pkin(3) * t546 + qJD(4) + t571;
t445 = -pkin(4) * t649 - pkin(9) * t619 + t526;
t408 = t435 * t601 + t445 * t597;
t398 = -pkin(10) * t499 + t408;
t669 = qJD(6) * t596;
t396 = t398 * t669;
t487 = t595 * t495;
t436 = t486 * t711 - t487;
t434 = -t592 * pkin(4) - t436;
t421 = t499 * pkin(5) + t434;
t732 = t421 * t450 + t396;
t477 = t518 * t711 - t595 * t608;
t670 = qJD(5) * t601;
t427 = t601 * t477 + t592 * t670 - t619 * t671;
t662 = qJD(2) * t718;
t641 = qJD(1) * t662;
t557 = t603 * t641;
t673 = qJD(3) * t598;
t646 = -t598 * t557 - t567 * t673;
t556 = t599 * t641;
t725 = t602 * (qJD(3) * t555 - t556);
t430 = -qJ(4) * t608 - t546 * qJD(4) + t646 + t725;
t647 = t598 * t556 - t602 * t557;
t612 = qJD(3) * t627 + t647;
t609 = -qJ(4) * t518 + qJD(4) * t548 + t612;
t403 = t430 * t711 + t595 * t609;
t607 = pkin(3) * t608 + qJD(2) * t589;
t415 = t476 * pkin(4) - t477 * pkin(9) + t607;
t414 = t601 * t415;
t611 = -qJD(5) * t408 - t403 * t597 + t414;
t373 = pkin(5) * t476 - pkin(10) * t427 + t611;
t701 = t477 * t597;
t428 = qJD(5) * t501 + t701;
t617 = t601 * t403 + t597 * t415 - t435 * t671 + t445 * t670;
t374 = -pkin(10) * t428 + t617;
t652 = t600 * t373 - t596 * t374;
t731 = t421 * t628 + t652;
t730 = t476 * MDP(31) + (-t450 ^ 2 + t628 ^ 2) * MDP(28) - t450 * MDP(27) * t628;
t729 = -0.2e1 * t666;
t728 = MDP(5) * (t599 ^ 2 - t603 ^ 2);
t521 = -t595 * t562 + t564 * t711;
t478 = t563 * t521;
t726 = t599 * MDP(4);
t654 = t711 * t598;
t680 = t711 * t502 - t503 * t595 + (t595 * t602 + t654) * t713;
t724 = t740 * pkin(5);
t723 = t597 * t458 - t601 * t679;
t719 = t476 * t563 + t506 * t734;
t651 = t427 * t596 + t600 * t428;
t387 = -qJD(6) * t628 + t651;
t716 = t601 * pkin(5);
t591 = t601 * pkin(10);
t587 = pkin(2) * t602 + pkin(3);
t540 = pkin(2) * t654 + t595 * t587;
t534 = pkin(9) + t540;
t715 = -pkin(10) - t534;
t583 = pkin(3) * t595 + pkin(9);
t714 = -pkin(10) - t583;
t407 = -t435 * t597 + t601 * t445;
t397 = -pkin(10) * t501 + t407;
t389 = pkin(5) * t721 + t397;
t709 = t389 * t600;
t708 = t398 * t600;
t707 = t427 * t597;
t706 = t434 * t649;
t705 = t450 * t619;
t704 = t628 * t619;
t524 = t592 * t562;
t615 = t564 * qJD(3);
t525 = qJD(2) * t564 + t615;
t482 = -t524 * t711 - t595 * t525;
t700 = t482 * t597;
t699 = t482 * t601;
t698 = t499 * t619;
t697 = t501 * t619;
t694 = t521 * t597;
t693 = t521 * t601;
t691 = t571 * t548;
t604 = qJD(2) ^ 2;
t688 = t599 * t604;
t511 = -qJ(4) * t564 - t572 * t602 - t573 * t598;
t626 = t572 * t598 - t573 * t602;
t512 = -qJ(4) * t562 - t626;
t471 = t595 * t511 + t512 * t711;
t466 = t601 * t471;
t472 = t601 * t476;
t687 = t603 * t604;
t605 = qJD(1) ^ 2;
t686 = t603 * t605;
t443 = t494 * t711 - t487;
t685 = t601 * t443 + t597 * t461;
t520 = t562 * t711 + t564 * t595;
t625 = pkin(3) * t562 + t588;
t469 = pkin(4) * t520 - pkin(9) * t521 + t625;
t683 = t597 * t469 + t466;
t681 = t724 + t680;
t668 = qJD(6) * t600;
t665 = pkin(10) * t695;
t590 = t599 * t712;
t663 = t600 * t427 - t596 * t428 - t499 * t668;
t661 = t521 * t671;
t432 = t434 * t670;
t659 = -pkin(2) * t592 - t555;
t658 = pkin(3) * t525 + t590;
t657 = qJD(5) * t715;
t656 = qJD(5) * t714;
t653 = pkin(1) * t729;
t402 = t430 * t595 - t711 * t609;
t650 = t443 * t597 - t601 * t461;
t566 = t599 * t662;
t568 = t603 * t662;
t618 = -t602 * t566 - t598 * t568 - t572 * t672 - t573 * t673;
t454 = -qJ(4) * t525 - qJD(4) * t562 + t618;
t610 = qJD(3) * t626 + t566 * t598 - t602 * t568;
t455 = qJ(4) * t524 - qJD(4) * t564 + t610;
t411 = t454 * t595 - t711 * t455;
t442 = t494 * t595 + t655;
t470 = -t711 * t511 + t512 * t595;
t642 = qJD(6) * t389 + t374;
t640 = pkin(5) * t619 - t591 * t649;
t584 = -pkin(3) * t711 - pkin(4);
t639 = -t442 + t724;
t638 = t402 * t597 + t408 * t619 + t432;
t637 = -t561 * t476 + t506 * t735;
t559 = t583 * t601 + t591;
t636 = qJD(6) * t559 - t601 * t656 + t640 - t650;
t528 = t534 * t601 + t591;
t635 = qJD(6) * t528 - t601 * t657 + t640 - t741;
t558 = t714 * t597;
t634 = -qJD(6) * t558 - t597 * t656 - t665 + t685;
t527 = t715 * t597;
t633 = -qJD(6) * t527 - t597 * t657 - t665 + t723;
t539 = -pkin(2) * t690 + t587 * t711;
t379 = t389 * t596 + t708;
t632 = t402 * t521 - t471 * t476;
t631 = -t476 * t534 - t706;
t630 = -t476 * t583 - t706;
t629 = t436 * t649 + t437 * t619;
t533 = -pkin(4) - t539;
t624 = -t721 * t740 + t472;
t623 = -t402 * t601 - t407 * t619 + t434 * t671;
t622 = t571 * t546 - t646;
t621 = t521 * t670 + t700;
t620 = -t661 + t699;
t412 = t454 * t711 + t595 * t455;
t481 = -t524 * t595 + t525 * t711;
t420 = pkin(4) * t481 - pkin(9) * t482 + t658;
t616 = t601 * t412 + t597 * t420 + t469 * t670 - t471 * t671;
t386 = -t501 * t669 + t663;
t378 = -t398 * t596 + t709;
t385 = pkin(5) * t428 + t402;
t614 = -t378 * t619 + t385 * t561 - t421 * t735;
t613 = t379 * t619 + t385 * t563 + t421 * t734;
t606 = -t548 * t546 * MDP(11) + (-t386 * t561 - t387 * t563 - t450 * t734 - t628 * t735) * MDP(28) + (t386 * t563 - t628 * t734) * MDP(27) + ((t427 - t738) * t601 + (-t501 * t721 - t428) * t597) * MDP(21) + (t704 + t719) * MDP(29) + (t637 + t705) * MDP(30) + (t624 + t698) * MDP(23) + (-t697 - t733) * MDP(22) + (t501 * t644 + t707) * MDP(20) + (t546 * t592 + t518) * MDP(13) + (-t548 * t592 - t608) * MDP(14) + (-t546 ^ 2 + t548 ^ 2) * MDP(12) + (-MDP(24) * t721 - MDP(31) * t506) * t619;
t569 = t584 - t716;
t530 = t533 - t716;
t479 = t561 * t521;
t464 = t601 * t469;
t439 = t476 * t520;
t433 = pkin(5) * t694 + t470;
t419 = t601 * t420;
t409 = -pkin(10) * t694 + t683;
t406 = pkin(5) * t520 - pkin(10) * t693 - t471 * t597 + t464;
t391 = t482 * t689 - t596 * t661 - t669 * t694 + (t693 * t720 + t700) * t600;
t390 = -t478 * t720 - t561 * t482;
t388 = pkin(5) * t621 + t411;
t377 = -pkin(10) * t621 + t616;
t375 = -pkin(10) * t699 + pkin(5) * t481 - t412 * t597 + t419 + (-t466 + (pkin(10) * t521 - t469) * t597) * qJD(5);
t1 = [(t546 * t590 + t571 * t525 + (t588 * t615 + (t599 * pkin(2) * t562 + t564 * t588) * qJD(2)) * qJD(1)) * MDP(16) + (-t403 * t520 + t411 * t619 + t412 * t649 - t436 * t482 - t437 * t481 + t470 * t477 + t632) * MDP(18) + ((-t499 * t601 - t501 * t597) * t482 + (-t707 - t428 * t601 + (t499 * t597 - t501 * t601) * qJD(5)) * t521) * MDP(21) + (t427 * t693 + t501 * t620) * MDP(20) - MDP(7) * t688 + (pkin(7) * t688 + t603 * t653) * MDP(10) + (-pkin(7) * t687 + t599 * t653) * MDP(9) + (t402 * t470 + t403 * t471 - t436 * t411 + t437 * t412 + t526 * t658 + t607 * t625) * MDP(19) + ((t375 * t600 - t377 * t596) * t506 + (t406 * t600 - t409 * t596) * t476 + t652 * t520 + t378 * t481 + t388 * t450 + t433 * t387 + t385 * t478 + t421 * t391 + ((-t406 * t596 - t409 * t600) * t506 - t379 * t520) * qJD(6)) * MDP(32) + (-t379 * t481 - t385 * t479 + t433 * t386 - t388 * t628 + t421 * t390 + t396 * t520 + (-(-qJD(6) * t409 + t375) * t506 - t406 * t476 - t373 * t520) * t596 + (-(qJD(6) * t406 + t377) * t506 - t409 * t476 - t642 * t520) * t600) * MDP(33) + (t386 * t520 + t390 * t506 - t476 * t479 - t481 * t628) * MDP(29) + (-t386 * t478 + t387 * t479 - t390 * t450 + t391 * t628) * MDP(28) + (-t386 * t479 - t390 * t628) * MDP(27) + (t402 * t693 - t408 * t481 + t411 * t501 + t470 * t427 + t434 * t620 - t476 * t683 - t520 * t617 - t616 * t721) * MDP(26) + (-t428 * t520 - t476 * t694 - t481 * t499 - t621 * t721) * MDP(23) + (t427 * t520 + t472 * t521 + t481 * t501 + t620 * t721) * MDP(22) + 0.2e1 * t660 * t726 + ((-t471 * t670 + t419) * t721 + t464 * t476 + (-t435 * t670 + t414) * t520 + t407 * t481 + t411 * t499 + t470 * t428 + t521 * t432 + ((-qJD(5) * t469 - t412) * t721 + (-qJD(5) * t445 - t403) * t520 + t434 * t482 + t632) * t597) * MDP(25) + (t481 * t721 + t439) * MDP(24) + (-t524 * MDP(13) - t525 * MDP(14) + MDP(16) * t610 - MDP(17) * t618) * t592 + (t588 * t518 - t571 * t524 + (-t548 + t722) * t590) * MDP(17) + (t518 * t564 + t524 * t548) * MDP(11) + (-t518 * t562 + t524 * t546 + t548 * t525 - t564 * t608) * MDP(12) + (-t387 * t520 - t391 * t506 - t450 * t481 - t476 * t478) * MDP(30) + (t481 * t506 + t439) * MDP(31) + t728 * t729 + MDP(6) * t687; (t548 * t589 + t678 * t592 + (qJD(3) * t659 + t556) * t602 + t622) * MDP(17) + (-t546 * t589 + t691 - t645 * t592 + (t598 * t659 - t553) * qJD(3) + t647) * MDP(16) + (-(t527 * t596 + t528 * t600) * t476 + t530 * t386 + (t596 * t635 + t600 * t633) * t506 - t681 * t628 + t613) * MDP(33) + t605 * t728 + t606 + (t533 * t427 + t631 * t601 + t680 * t501 + (t534 * t671 + t723) * t721 + t638) * MDP(26) + (-t476 * t540 - t477 * t539 + t619 * t680 + t649 * t679 + t629) * MDP(18) + (t533 * t428 + t631 * t597 + t680 * t499 + (-t534 * t670 + t741) * t721 + t623) * MDP(25) + ((t527 * t600 - t528 * t596) * t476 + t530 * t387 + (t596 * t633 - t600 * t635) * t506 + t681 * t450 + t614) * MDP(32) + (t403 * t540 - t402 * t539 - t526 * (t589 - t717) + t679 * t437 - t680 * t436) * MDP(19) - t686 * t726 + (MDP(9) * t599 * t605 + MDP(10) * t686) * pkin(1); (t592 * t648 + t622 - t725) * MDP(17) + t606 + ((t558 * t600 - t559 * t596) * t476 + t569 * t387 + (t596 * t634 - t600 * t636) * t506 + t639 * t450 + t614) * MDP(32) + (t436 * t442 - t437 * t443 + (-t402 * t711 + t403 * t595 + t526 * t548) * pkin(3)) * MDP(19) + (-t442 * t619 - t443 * t649 + (-t476 * t595 - t477 * t711) * pkin(3) + t629) * MDP(18) + (-t592 * t627 + t612 + t691) * MDP(16) + (-(t558 * t596 + t559 * t600) * t476 + t569 * t386 + (t596 * t636 + t600 * t634) * t506 - t639 * t628 + t613) * MDP(33) + (t584 * t428 - t442 * t499 + t630 * t597 + (-t583 * t670 + t650) * t721 + t623) * MDP(25) + (t584 * t427 - t442 * t501 + t630 * t601 + (t583 * t671 + t685) * t721 + t638) * MDP(26); (-t619 ^ 2 - t649 ^ 2) * MDP(18) + (t436 * t619 - t437 * t649 + t607) * MDP(19) + (t624 - t698) * MDP(25) + (-t697 + t733) * MDP(26) + (t637 - t705) * MDP(32) + (t704 - t719) * MDP(33); t501 * t499 * MDP(20) + (-t499 ^ 2 + t501 ^ 2) * MDP(21) + (t427 + t738) * MDP(22) + (-t701 + (-qJD(5) + t721) * t501) * MDP(23) + t476 * MDP(24) + (t408 * t721 - t434 * t501 + t611) * MDP(25) + (t407 * t721 + t434 * t499 - t617) * MDP(26) + (t386 + t739) * MDP(29) + (-t387 - t737) * MDP(30) + (-(-t397 * t596 - t708) * t506 - t379 * qJD(6) + (-t450 * t501 + t476 * t600 - t506 * t669) * pkin(5) + t731) * MDP(32) + ((-t398 * t506 - t373) * t596 + (t397 * t506 - t642) * t600 + (-t476 * t596 + t501 * t628 - t506 * t668) * pkin(5) + t732) * MDP(33) + t730; (t663 + t739) * MDP(29) + (-t651 - t737) * MDP(30) + (t379 * t506 + t731) * MDP(32) + (-t596 * t373 - t600 * t374 + t378 * t506 + t732) * MDP(33) + (-MDP(29) * t696 + MDP(30) * t628 - MDP(32) * t379 - MDP(33) * t709) * qJD(6) + t730;];
tauc  = t1;
