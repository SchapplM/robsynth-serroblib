% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:17:00
% EndTime: 2019-03-09 06:17:12
% DurationCPUTime: 8.72s
% Computational Cost: add. (8409->533), mult. (19873->669), div. (0->0), fcn. (15449->14), ass. (0->241)
t601 = cos(pkin(10));
t608 = cos(qJ(3));
t678 = qJD(1) * t608;
t576 = t601 * t678;
t600 = sin(pkin(10));
t605 = sin(qJ(3));
t679 = qJD(1) * t605;
t661 = t600 * t679;
t537 = t576 - t661;
t549 = t600 * t608 + t601 * t605;
t538 = t549 * qJD(1);
t494 = pkin(3) * t538 - pkin(8) * t537;
t607 = cos(qJ(4));
t480 = t607 * t494;
t604 = sin(qJ(4));
t737 = pkin(8) + pkin(9);
t664 = qJD(4) * t737;
t726 = pkin(7) + qJ(2);
t566 = t726 * t600;
t550 = qJD(1) * t566;
t567 = t726 * t601;
t551 = qJD(1) * t567;
t752 = -t550 * t608 - t605 * t551;
t768 = pkin(4) * t538 - t604 * t752 + t480 + (-pkin(9) * t537 + t664) * t607;
t685 = t604 * t494 + t607 * t752;
t714 = t537 * t604;
t767 = -pkin(9) * t714 + t604 * t664 + t685;
t672 = qJD(1) * qJD(2);
t740 = qJDD(1) * t726 + t672;
t521 = t740 * t600;
t522 = t740 * t601;
t676 = qJD(3) * t608;
t677 = qJD(3) * t605;
t635 = -t608 * t521 - t605 * t522 + t550 * t677 - t551 * t676;
t724 = qJDD(3) * pkin(3);
t446 = -t635 - t724;
t528 = qJD(4) - t537;
t598 = pkin(10) + qJ(3);
t588 = sin(t598);
t589 = cos(t598);
t606 = sin(qJ(1));
t609 = cos(qJ(1));
t646 = g(1) * t609 + g(2) * t606;
t619 = -g(3) * t589 + t588 * t646;
t766 = qJD(4) * pkin(8) * t528 + t446 - t619;
t675 = qJD(4) * t604;
t765 = t675 - t714;
t669 = qJDD(1) * t608;
t670 = qJDD(1) * t605;
t665 = qJD(3) * t576 + t600 * t669 + t601 * t670;
t622 = qJD(3) * t661 - t665;
t762 = qJD(3) * qJD(4) - t622;
t452 = t604 * qJDD(3) - t538 * t675 + t607 * t762;
t508 = qJD(3) * t607 - t604 * t538;
t509 = qJD(3) * t604 + t538 * t607;
t603 = sin(qJ(5));
t674 = qJD(4) * t607;
t666 = t538 * t674 + t604 * t762;
t633 = t607 * qJDD(3) - t666;
t736 = cos(qJ(5));
t657 = t736 * qJD(5);
t673 = qJD(5) * t603;
t409 = -t736 * t452 - t508 * t657 + t509 * t673 - t603 * t633;
t631 = t603 * t508 + t509 * t736;
t410 = qJD(5) * t631 + t603 * t452 - t736 * t633;
t456 = -t736 * t508 + t509 * t603;
t454 = t456 ^ 2;
t540 = t549 * qJD(3);
t573 = t601 * t669;
t641 = -t600 * t670 + t573;
t497 = qJD(1) * t540 - t641;
t491 = qJDD(4) + t497;
t488 = qJDD(5) + t491;
t524 = qJD(5) + t528;
t738 = t631 ^ 2;
t764 = t488 * MDP(26) + (t524 * t631 - t410) * MDP(25) + t456 * MDP(22) * t631 + (t456 * t524 - t409) * MDP(24) + (-t454 + t738) * MDP(23);
t763 = t456 * qJ(6);
t663 = t736 * t604;
t553 = t603 * t607 + t663;
t746 = qJD(4) + qJD(5);
t502 = t746 * t553;
t687 = t553 * t537 - t502;
t699 = t603 * t604;
t630 = t736 * t607 - t699;
t748 = t736 * qJD(4) + t657;
t686 = t630 * t537 - t748 * t607 + t746 * t699;
t548 = t600 * t605 - t608 * t601;
t539 = t548 * qJD(3);
t659 = t549 * t674;
t761 = -t539 * t604 + t659;
t751 = g(1) * t606 - g(2) * t609;
t760 = qJDD(2) - t751;
t492 = -qJD(3) * pkin(3) - t752;
t453 = -pkin(4) * t508 + t492;
t599 = qJ(4) + qJ(5);
t592 = cos(t599);
t704 = t592 * t606;
t591 = sin(t599);
t705 = t591 * t609;
t513 = -t589 * t704 + t705;
t703 = t592 * t609;
t706 = t591 * t606;
t515 = t589 * t703 + t706;
t580 = pkin(2) * t601 + pkin(1);
t561 = -qJD(1) * t580 + qJD(2);
t467 = -pkin(3) * t537 - pkin(8) * t538 + t561;
t499 = -t605 * t550 + t608 * t551;
t493 = qJD(3) * pkin(8) + t499;
t439 = t467 * t604 + t493 * t607;
t560 = -qJDD(1) * t580 + qJDD(2);
t448 = t497 * pkin(3) + pkin(8) * t622 + t560;
t444 = t607 * t448;
t640 = -t521 * t605 + t522 * t608;
t445 = qJDD(3) * pkin(8) + qJD(3) * t752 + t640;
t390 = pkin(4) * t491 - pkin(9) * t452 - qJD(4) * t439 - t445 * t604 + t444;
t626 = t607 * t445 + t604 * t448 + t467 * t674 - t493 * t675;
t393 = pkin(9) * t633 + t626;
t438 = t607 * t467 - t493 * t604;
t426 = -pkin(9) * t509 + t438;
t416 = pkin(4) * t528 + t426;
t427 = pkin(9) * t508 + t439;
t648 = -t603 * t390 - t736 * t393 - t416 * t657 + t427 * t673;
t728 = g(3) * t592;
t759 = g(1) * t515 - g(2) * t513 + t453 * t456 + t588 * t728 + t648;
t757 = qJ(6) * t631;
t419 = pkin(5) * t456 + qJD(6) + t453;
t756 = t419 * t631;
t755 = t588 * t751;
t496 = pkin(3) * t548 - pkin(8) * t549 - t580;
t485 = t607 * t496;
t507 = -t566 * t605 + t567 * t608;
t710 = t549 * t607;
t434 = pkin(4) * t548 - pkin(9) * t710 - t507 * t604 + t485;
t500 = t607 * t507;
t684 = t604 * t496 + t500;
t711 = t549 * t604;
t441 = -pkin(9) * t711 + t684;
t688 = t603 * t434 + t736 * t441;
t753 = pkin(4) * t765 - t499;
t506 = t566 * t608 + t605 * t567;
t568 = t737 * t604;
t569 = t737 * t607;
t683 = -t603 * t568 + t736 * t569;
t750 = -t683 * qJD(5) + t767 * t603 - t736 * t768;
t749 = t568 * t657 + t569 * t673 + t603 * t768 + t767 * t736;
t747 = qJ(2) * qJDD(1);
t512 = t589 * t706 + t703;
t514 = -t589 * t705 + t704;
t730 = g(3) * t588;
t745 = -g(1) * t514 + g(2) * t512 + t591 * t730;
t423 = t736 * t427;
t402 = t603 * t416 + t423;
t615 = -qJD(5) * t402 + t736 * t390 - t603 * t393;
t744 = -t453 * t631 + t615 + t745;
t743 = t488 * t553 - t524 * t686;
t742 = t646 * t589 + t730;
t741 = t409 * t630 - t631 * t687;
t593 = t607 * pkin(4);
t727 = pkin(3) + t593;
t725 = qJDD(1) * pkin(1);
t722 = t452 * t604;
t721 = t456 * t538;
t720 = t631 * t538;
t718 = t508 * t528;
t717 = t508 * t538;
t716 = t509 * t528;
t715 = t509 * t538;
t712 = t539 * t607;
t562 = pkin(4) * t604 + pkin(5) * t591;
t708 = t562 * t589;
t701 = t601 * MDP(4);
t421 = t603 * t427;
t698 = t604 * t491;
t697 = t604 * t606;
t696 = t604 * t609;
t695 = t606 * t607;
t477 = t607 * t491;
t694 = t607 * t609;
t401 = t736 * t416 - t421;
t395 = t401 - t757;
t394 = pkin(5) * t524 + t395;
t693 = -t395 + t394;
t692 = t736 * t426 - t421;
t691 = qJ(6) * t687 + qJD(6) * t630 - t749;
t690 = -pkin(5) * t538 + qJ(6) * t686 - t553 * qJD(6) + t750;
t682 = t562 + t726;
t563 = pkin(5) * t592 + t593;
t681 = t600 ^ 2 + t601 ^ 2;
t660 = t549 * t675;
t655 = -t426 * t603 - t423;
t653 = t736 * t434 - t441 * t603;
t651 = -t736 * t568 - t569 * t603;
t650 = t528 * t607;
t649 = -qJD(4) * t467 - t445;
t469 = qJD(2) * t549 - t566 * t677 + t567 * t676;
t647 = 0.2e1 * t681;
t644 = -t553 * t410 + t456 * t686;
t643 = t630 * t488 + t524 * t687;
t470 = pkin(4) * t711 + t506;
t642 = -t493 * t674 + t444;
t557 = pkin(3) + t563;
t595 = -qJ(6) - t737;
t638 = t557 * t589 - t588 * t595;
t636 = t725 - t760;
t634 = -t528 * t765 + t477;
t447 = pkin(4) * t761 + t469;
t632 = t580 + t638;
t628 = -t660 - t712;
t627 = -pkin(8) * t491 + t492 * t528;
t468 = -t548 * qJD(2) - t506 * qJD(3);
t495 = pkin(3) * t540 + pkin(8) * t539;
t625 = t607 * t468 + t604 * t495 + t496 * t674 - t507 * t675;
t481 = t607 * t495;
t408 = pkin(9) * t712 + pkin(4) * t540 - t468 * t604 + t481 + (-t500 + (pkin(9) * t549 - t496) * t604) * qJD(4);
t412 = -pkin(9) * t761 + t625;
t624 = t603 * t408 + t736 * t412 + t434 * t657 - t441 * t673;
t617 = t647 * t672 - t646;
t614 = -t688 * qJD(5) + t736 * t408 - t603 * t412;
t413 = -pkin(4) * t633 + t446;
t391 = t410 * pkin(5) + qJDD(6) + t413;
t585 = pkin(4) * t736 + pkin(5);
t532 = t589 * t694 + t697;
t531 = -t589 * t696 + t695;
t530 = -t589 * t695 + t696;
t529 = t589 * t697 + t694;
t483 = t630 * t549;
t482 = t553 * t549;
t475 = qJ(6) * t630 + t683;
t474 = -qJ(6) * t553 + t651;
t418 = -t539 * t663 - t603 * t660 - t673 * t711 + (-t539 * t603 + t748 * t549) * t607;
t417 = t502 * t549 + t630 * t539;
t404 = -qJ(6) * t482 + t688;
t403 = pkin(5) * t548 - qJ(6) * t483 + t653;
t398 = t692 - t757;
t397 = t655 + t763;
t396 = t402 - t763;
t387 = -qJ(6) * t418 - qJD(6) * t482 + t624;
t386 = t540 * pkin(5) + t417 * qJ(6) - t483 * qJD(6) + t614;
t385 = -qJ(6) * t410 - qJD(6) * t456 - t648;
t384 = t488 * pkin(5) + t409 * qJ(6) - qJD(6) * t631 + t615;
t1 = [(-qJD(3) * t469 - qJDD(3) * t506 - t497 * t580 + t540 * t561 + t548 * t560 + t589 * t751) * MDP(13) + t751 * MDP(2) + (-t468 * qJD(3) - t507 * qJDD(3) - t561 * t539 + t560 * t549 + t580 * t622 - t755) * MDP(14) + (-t384 * t483 - t385 * t482 - t386 * t631 - t387 * t456 + t394 * t417 - t396 * t418 + t403 * t409 - t404 * t410 + t755) * MDP(29) + t646 * MDP(3) + (t452 * t548 + t477 * t549 + t509 * t540 + t528 * t628) * MDP(17) + ((-t507 * t674 + t481) * t528 + t485 * t491 + t642 * t548 + t438 * t540 - t469 * t508 - t506 * t633 + t492 * t659 - g(1) * t530 - g(2) * t532 + ((-qJD(4) * t496 - t468) * t528 - t507 * t491 + t649 * t548 + t446 * t549 - t492 * t539) * t604) * MDP(20) + (-qJD(3) * t539 + qJDD(3) * t549) * MDP(10) + (-(t508 * t607 - t509 * t604) * t539 + (t607 * t633 - t722 + (-t604 * t508 - t509 * t607) * qJD(4)) * t549) * MDP(16) + (-g(1) * t512 - g(2) * t514 - t402 * t540 - t470 * t409 + t413 * t483 - t453 * t417 + t447 * t631 - t488 * t688 - t524 * t624 + t548 * t648) * MDP(28) + (-t409 * t548 - t417 * t524 + t483 * t488 + t540 * t631) * MDP(24) + (t409 * t482 - t410 * t483 + t417 * t456 - t418 * t631) * MDP(23) + (-t409 * t483 - t417 * t631) * MDP(22) + (-t538 * t539 - t549 * t622) * MDP(8) + (-t549 * t497 - t539 * t537 - t538 * t540 + t548 * t622) * MDP(9) + (-MDP(5) * t600 + t701) * (t636 + t725) + (pkin(1) * t636 + (t681 * t747 + t617) * qJ(2)) * MDP(7) + (t647 * t747 + t617) * MDP(6) + (t508 * t540 - t528 * t761 + t548 * t633 - t549 * t698) * MDP(18) + qJDD(1) * MDP(1) + (t385 * t404 + t396 * t387 + t384 * t403 + t394 * t386 + t391 * (pkin(5) * t482 + t470) + t419 * (pkin(5) * t418 + t447) + (-g(1) * t682 - g(2) * t632) * t609 + (g(1) * t632 - g(2) * t682) * t606) * MDP(30) + (t452 * t710 + t509 * t628) * MDP(15) + (-g(1) * t529 - g(2) * t531 - t439 * t540 + t446 * t710 + t506 * t452 + t469 * t509 - t491 * t684 + t492 * t628 - t528 * t625 - t548 * t626) * MDP(21) + (-g(1) * t513 - g(2) * t515 + t401 * t540 + t470 * t410 + t413 * t482 + t453 * t418 + t447 * t456 + t488 * t653 + t524 * t614 + t548 * t615) * MDP(27) + (t488 * t548 + t524 * t540) * MDP(26) + (-t410 * t548 - t418 * t524 - t456 * t540 - t482 * t488) * MDP(25) + (-qJD(3) * t540 - qJDD(3) * t548) * MDP(11) + (t491 * t548 + t528 * t540) * MDP(19); t760 * MDP(7) - t573 * MDP(13) + t665 * MDP(14) + (t634 + t717) * MDP(20) + (-t528 ^ 2 * t607 - t698 - t715) * MDP(21) + (t643 - t721) * MDP(27) + (-t720 - t743) * MDP(28) + (t644 + t741) * MDP(29) + (t384 * t630 + t385 * t553 + t394 * t687 - t396 * t686 - t419 * t538 - t751) * MDP(30) + (-t701 - pkin(1) * MDP(7) + (MDP(13) * t605 + MDP(5)) * t600) * qJDD(1) + ((t600 * t678 + t601 * t679 + t538) * MDP(13) + (t537 - t661) * MDP(14)) * qJD(3) + (-MDP(7) * qJ(2) - MDP(6)) * qJD(1) ^ 2 * t681; -t537 ^ 2 * MDP(9) + ((-t537 - t661) * qJD(3) + t665) * MDP(10) + t641 * MDP(11) + qJDD(3) * MDP(12) + (qJD(3) * t499 + t619 + t635) * MDP(13) + (-t537 * t561 - t640 + t742) * MDP(14) + (t509 * t650 + t722) * MDP(15) + ((t452 + t718) * t607 + (t633 - t716) * t604) * MDP(16) + (t528 * t650 + t698 - t715) * MDP(17) + (t634 - t717) * MDP(18) + (-pkin(3) * t666 - t480 * t528 + t499 * t508 + (t528 * t752 + t627) * t604 + (t724 - t766) * t607) * MDP(20) + (-pkin(3) * t452 - t499 * t509 + t685 * t528 + t604 * t766 + t627 * t607) * MDP(21) + (-t409 * t553 - t631 * t686) * MDP(22) + (t644 - t741) * MDP(23) + (-t720 + t743) * MDP(24) + (t643 + t721) * MDP(25) + (-t727 * t410 - t413 * t630 + t488 * t651 - t589 * t728 + (g(1) * t703 + g(2) * t704) * t588 + t750 * t524 + t753 * t456 - t687 * t453) * MDP(27) + (t409 * t727 + t413 * t553 - t686 * t453 - t683 * t488 + t749 * t524 - t591 * t619 + t753 * t631) * MDP(28) + (-t384 * t553 + t385 * t630 + t394 * t686 + t396 * t687 + t409 * t474 - t410 * t475 - t456 * t691 - t631 * t690 - t742) * MDP(29) + (t385 * t475 + t384 * t474 + t391 * (-pkin(5) * t630 - t727) - g(3) * t638 + (-pkin(5) * t687 + t753) * t419 + t691 * t396 + t690 * t394 + t646 * (t557 * t588 + t589 * t595)) * MDP(30) + (-t561 * MDP(13) - t528 * MDP(19) - t438 * MDP(20) + t439 * MDP(21) - t524 * MDP(26) - t401 * MDP(27) + t402 * MDP(28) - t537 * MDP(8) + MDP(9) * t538) * t538; -t509 * t508 * MDP(15) + (-t508 ^ 2 + t509 ^ 2) * MDP(16) + (t452 - t718) * MDP(17) + (t633 + t716) * MDP(18) + t491 * MDP(19) + (-g(1) * t531 + g(2) * t529 + t439 * t528 - t492 * t509 + (t649 + t730) * t604 + t642) * MDP(20) + (g(1) * t532 - g(2) * t530 + t438 * t528 - t492 * t508 + t607 * t730 - t626) * MDP(21) + (-t655 * t524 + (-t509 * t456 + t488 * t736 - t524 * t673) * pkin(4) + t744) * MDP(27) + (t692 * t524 + (-t603 * t488 - t509 * t631 - t524 * t657) * pkin(4) + t759) * MDP(28) + (-t394 * t456 + t396 * t631 + t397 * t631 + t398 * t456 + t585 * t409 + (-t410 * t603 + (-t456 * t736 + t603 * t631) * qJD(5)) * pkin(4)) * MDP(29) + (t384 * t585 - t396 * t398 - t394 * t397 - pkin(5) * t756 - g(1) * (t563 * t606 - t609 * t708) - g(2) * (-t563 * t609 - t606 * t708) + t562 * t730 + (t385 * t603 - t419 * t509 + (-t394 * t603 + t396 * t736) * qJD(5)) * pkin(4)) * MDP(30) + t764; (t402 * t524 + t744) * MDP(27) + (t401 * t524 + t759) * MDP(28) + (pkin(5) * t409 - t456 * t693) * MDP(29) + (t693 * t396 + (t384 + t745 - t756) * pkin(5)) * MDP(30) + t764; (-t454 - t738) * MDP(29) + (t394 * t631 + t396 * t456 + t391 - t619) * MDP(30);];
tau  = t1;
