% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:09:44
% EndTime: 2019-03-10 01:09:59
% DurationCPUTime: 8.60s
% Computational Cost: add. (10417->521), mult. (24547->673), div. (0->0), fcn. (18059->8), ass. (0->236)
t659 = cos(qJ(2));
t772 = -pkin(8) - pkin(7);
t632 = t772 * t659;
t622 = qJD(1) * t632;
t655 = sin(qJ(3));
t598 = t655 * t622;
t656 = sin(qJ(2));
t630 = t772 * t656;
t620 = qJD(1) * t630;
t770 = cos(qJ(3));
t552 = t620 * t770 + t598;
t707 = t770 * qJD(3);
t797 = -pkin(2) * t707 + t552;
t708 = qJD(1) * t770;
t726 = qJD(1) * t656;
t592 = t655 * t726 - t659 * t708;
t750 = t655 * t659;
t594 = -qJD(1) * t750 - t656 * t708;
t544 = -pkin(3) * t594 + pkin(9) * t592;
t526 = pkin(2) * t726 + t544;
t654 = sin(qJ(4));
t658 = cos(qJ(4));
t798 = -t658 * t526 + t797 * t654;
t796 = t654 * t526 + t658 * t797;
t756 = t592 * t658;
t690 = -t594 * pkin(4) + pkin(10) * t756;
t642 = pkin(2) * t655 + pkin(9);
t767 = -pkin(10) - t642;
t705 = qJD(4) * t767;
t795 = -t658 * t705 + t690 - t798;
t766 = qJD(2) * pkin(2);
t602 = t620 + t766;
t548 = t602 * t770 + t598;
t697 = t658 * t544 - t548 * t654;
t771 = -pkin(9) - pkin(10);
t712 = qJD(4) * t771;
t794 = -t658 * t712 + t690 + t697;
t757 = t592 * t654;
t717 = pkin(10) * t757;
t791 = -t654 * t705 + t717 + t796;
t734 = t654 * t544 + t658 * t548;
t790 = -t654 * t712 + t717 + t734;
t657 = cos(qJ(5));
t653 = sin(qJ(5));
t752 = t653 * t658;
t614 = t654 * t657 + t752;
t718 = qJD(4) + qJD(5);
t555 = t718 * t614;
t789 = t614 * t592 + t555;
t753 = t653 * t654;
t612 = -t657 * t658 + t753;
t721 = qJD(5) * t657;
t723 = qJD(4) * t658;
t735 = t612 * t592 - t657 * t723 - t658 * t721 + t718 * t753;
t650 = qJD(2) + qJD(3);
t564 = t654 * t594 + t650 * t658;
t565 = -t594 * t658 + t650 * t654;
t758 = t565 * t653;
t500 = -t657 * t564 + t758;
t498 = t500 ^ 2;
t615 = t656 * t770 + t750;
t557 = t650 * t615;
t543 = t557 * qJD(1);
t683 = t564 * t653 + t657 * t565;
t773 = t683 ^ 2;
t793 = t543 * MDP(29) + (-t498 + t773) * MDP(26);
t792 = qJ(6) * t500;
t679 = -t655 * t656 + t659 * t770;
t668 = t679 * qJD(3);
t556 = qJD(2) * t679 + t668;
t762 = t556 * t654;
t788 = t615 * t723 + t762;
t719 = qJD(1) * qJD(2);
t787 = -0.2e1 * t719;
t786 = MDP(5) * (t656 ^ 2 - t659 ^ 2);
t785 = qJ(6) * t683;
t529 = -t650 * pkin(3) - t548;
t496 = -t564 * pkin(4) + t529;
t452 = t500 * pkin(5) + qJD(6) + t496;
t784 = t452 * t683;
t783 = t496 * t683;
t588 = qJD(4) + t592;
t578 = qJD(5) + t588;
t782 = t578 * t683;
t646 = -pkin(2) * t659 - pkin(1);
t547 = -pkin(3) * t679 - pkin(9) * t615 + t646;
t540 = t658 * t547;
t567 = t655 * t630 - t632 * t770;
t754 = t615 * t658;
t474 = -pkin(4) * t679 - pkin(10) * t754 - t567 * t654 + t540;
t560 = t658 * t567;
t731 = t654 * t547 + t560;
t755 = t615 * t654;
t483 = -pkin(10) * t755 + t731;
t737 = t653 * t474 + t657 * t483;
t781 = -qJ(6) * t789 - t612 * qJD(6);
t601 = t770 * t622;
t551 = t655 * t620 - t601;
t725 = qJD(3) * t655;
t689 = pkin(2) * t725 - t551;
t780 = t795 * t657;
t724 = qJD(4) * t654;
t779 = (t724 + t757) * pkin(4);
t607 = t767 * t654;
t649 = t658 * pkin(10);
t608 = t642 * t658 + t649;
t730 = t653 * t607 + t657 * t608;
t778 = t794 * t657;
t629 = t771 * t654;
t631 = pkin(9) * t658 + t649;
t729 = t653 * t629 + t657 * t631;
t722 = qJD(5) * t653;
t777 = -t629 * t721 + t631 * t722 + t653 * t794 + t657 * t790;
t776 = -t607 * t721 + t608 * t722 + t653 * t795 + t791 * t657;
t775 = t594 * pkin(5) + qJ(6) * t735 - qJD(6) * t614;
t774 = t770 * t630 + t655 * t632;
t663 = t556 * qJD(1);
t495 = t594 * t724 + t650 * t723 + t658 * t663;
t714 = t594 * t723 - t650 * t724 - t654 * t663;
t698 = t495 * t653 - t657 * t714;
t438 = qJD(5) * t683 + t698;
t769 = t612 * pkin(5);
t768 = t658 * pkin(4);
t765 = qJ(6) * t614;
t764 = t495 * t654;
t763 = t543 * t658;
t761 = t556 * t658;
t760 = t564 * t588;
t759 = t565 * t588;
t628 = t646 * qJD(1);
t524 = pkin(3) * t592 + pkin(9) * t594 + t628;
t549 = t655 * t602 - t601;
t530 = pkin(9) * t650 + t549;
t479 = t524 * t654 + t530 * t658;
t460 = pkin(10) * t564 + t479;
t454 = t653 * t460;
t751 = t654 * t543;
t660 = qJD(2) ^ 2;
t749 = t656 * t660;
t456 = t657 * t460;
t747 = t659 * t660;
t661 = qJD(1) ^ 2;
t746 = t659 * t661;
t478 = t658 * t524 - t530 * t654;
t459 = -pkin(10) * t565 + t478;
t451 = pkin(4) * t588 + t459;
t429 = t657 * t451 - t454;
t417 = t429 - t785;
t414 = pkin(5) * t578 + t417;
t745 = t414 - t417;
t744 = -t776 + t781;
t743 = -qJD(5) * t730 + t653 * t791 + t775 - t780;
t742 = t657 * t459 - t454;
t741 = -t777 + t781;
t740 = -qJD(5) * t729 + t653 * t790 + t775 - t778;
t728 = t779 + t689;
t720 = t683 * MDP(25);
t716 = t656 * t766;
t715 = t657 * t495 + t564 * t721 + t653 * t714;
t645 = -pkin(3) - t768;
t713 = qJD(2) * t772;
t710 = t615 * t724;
t519 = t529 * t723;
t706 = t656 * t719;
t704 = pkin(1) * t787;
t484 = t543 * pkin(3) + (-pkin(9) * t668 + (t656 * pkin(2) - pkin(9) * t679) * qJD(2)) * qJD(1);
t482 = t658 * t484;
t691 = qJD(1) * t713;
t603 = t656 * t691;
t604 = t659 * t691;
t487 = t602 * t707 + t603 * t770 + t655 * t604 + t622 * t725;
t666 = -qJD(4) * t479 - t487 * t654 + t482;
t411 = pkin(4) * t543 - pkin(10) * t495 + t666;
t675 = t654 * t484 + t658 * t487 + t524 * t723 - t530 * t724;
t416 = pkin(10) * t714 + t675;
t702 = t657 * t411 - t653 * t416;
t494 = pkin(3) * t557 - pkin(9) * t556 + t716;
t490 = t658 * t494;
t621 = t656 * t713;
t623 = t659 * t713;
t505 = qJD(3) * t774 + t770 * t621 + t655 * t623;
t424 = -pkin(10) * t761 + pkin(4) * t557 - t505 * t654 + t490 + (-t560 + (pkin(10) * t615 - t547) * t654) * qJD(4);
t674 = t654 * t494 + t658 * t505 + t547 * t723 - t567 * t724;
t432 = -pkin(10) * t788 + t674;
t701 = t657 * t424 - t432 * t653;
t700 = -t459 * t653 - t456;
t699 = t657 * t474 - t483 * t653;
t696 = t657 * t607 - t608 * t653;
t695 = t657 * t629 - t631 * t653;
t694 = t588 * t658;
t692 = t653 * t411 + t657 * t416 + t451 * t721 - t460 * t722;
t488 = t602 * t725 + t655 * t603 - t770 * t604 - t622 * t707;
t644 = -pkin(2) * t770 - pkin(3);
t688 = -t549 + t779;
t687 = -t479 * t594 + t488 * t654 + t519;
t430 = t451 * t653 + t456;
t684 = t529 * t592 - t543 * t642;
t525 = pkin(4) * t755 - t774;
t681 = pkin(5) * t789 + t779;
t680 = t478 * t594 - t488 * t658 + t529 * t724;
t677 = -t710 + t761;
t627 = t644 - t768;
t676 = t594 * t628 - t488;
t506 = t655 * t621 - t623 * t770 + t630 * t725 - t632 * t707;
t673 = t653 * t424 + t657 * t432 + t474 * t721 - t483 * t722;
t437 = t565 * t722 - t715;
t444 = -pkin(4) * t714 + t488;
t670 = t429 * t594 + t444 * t612 + t496 * t789;
t669 = -t430 * t594 + t444 * t614 - t496 * t735;
t465 = pkin(4) * t788 + t506;
t665 = -qJD(5) * t430 + t702;
t401 = pkin(5) * t543 + qJ(6) * t437 - qJD(6) * t683 + t665;
t403 = -qJ(6) * t438 - qJD(6) * t500 + t692;
t418 = t430 - t792;
t667 = -t401 * t614 - t403 * t612 + t414 * t735 - t418 * t789;
t413 = t438 * pkin(5) + t444;
t664 = t628 * t592 - t487;
t662 = (t437 * t612 - t438 * t614 + t500 * t735 - t683 * t789) * MDP(26) + (-t437 * t614 - t683 * t735) * MDP(25) + ((t495 + t760) * t658 + (t714 - t759) * t654) * MDP(19) + (t543 * t614 - t578 * t735 + t594 * t683) * MDP(27) + (-t500 * t594 - t543 * t612 - t578 * t789) * MDP(28) + (t565 * t694 + t764) * MDP(18) + (-t588 ^ 2 * t654 + t564 * t594 + t763) * MDP(21) + (t565 * t594 + t588 * t694 + t751) * MDP(20) + t663 * MDP(13) + (-t592 ^ 2 + t594 ^ 2) * MDP(12) + (-MDP(11) * t592 + MDP(22) * t588 + MDP(29) * t578) * t594 + (t592 * MDP(13) + (-qJD(1) * t615 - t594) * MDP(14)) * t650;
t643 = pkin(4) * t657 + pkin(5);
t606 = t612 * qJ(6);
t537 = t612 * t615;
t536 = t614 * t615;
t528 = -t606 + t729;
t527 = t695 - t765;
t514 = t543 * t679;
t513 = -t606 + t730;
t512 = t696 - t765;
t442 = t556 * t752 - t653 * t710 - t722 * t755 + (t718 * t754 + t762) * t657;
t441 = t555 * t615 + t556 * t612;
t436 = -qJ(6) * t536 + t737;
t433 = -pkin(5) * t679 + qJ(6) * t537 + t699;
t420 = t742 - t785;
t419 = t700 + t792;
t405 = -qJ(6) * t442 - qJD(6) * t536 + t673;
t404 = pkin(5) * t557 + qJ(6) * t441 - qJD(5) * t737 + qJD(6) * t537 + t701;
t1 = [((-t567 * t723 + t490) * t588 + t540 * t543 - (-t530 * t723 + t482) * t679 + t478 * t557 - t506 * t564 + t774 * t714 + t615 * t519 + ((-qJD(4) * t547 - t505) * t588 - t567 * t543 - (-qJD(4) * t524 - t487) * t679 + t488 * t615 + t529 * t556) * t654) * MDP(23) + (-t479 * t557 + t488 * t754 - t495 * t774 + t506 * t565 + t529 * t677 - t543 * t731 - t588 * t674 + t675 * t679) * MDP(24) + (t701 * t578 + t699 * t543 - t702 * t679 + t429 * t557 + t465 * t500 + t525 * t438 + t444 * t536 + t496 * t442 + (t430 * t679 - t578 * t737) * qJD(5)) * MDP(30) + (MDP(13) * t556 - MDP(14) * t557 - MDP(16) * t506 - MDP(17) * t505) * t650 + (t437 * t679 - t441 * t578 - t537 * t543 + t557 * t683) * MDP(27) + (t543 * t646 + t557 * t628 + (-qJD(1) * t679 + t592) * t716) * MDP(16) + (-t615 * t543 - t556 * t592 + t594 * t557 + t663 * t679) * MDP(12) + (-t430 * t557 - t525 * t437 - t496 * t441 - t444 * t537 + t465 * t683 - t543 * t737 - t578 * t673 + t679 * t692) * MDP(31) + (t438 * t679 - t442 * t578 - t500 * t557 - t536 * t543) * MDP(28) + (-t495 * t679 + t543 * t754 + t557 * t565 + t588 * t677) * MDP(20) + (t557 * t588 - t514) * MDP(22) + (t557 * t578 - t514) * MDP(29) + (t403 * t436 + t418 * t405 + t401 * t433 + t414 * t404 + t413 * (t536 * pkin(5) + t525) + t452 * (t442 * pkin(5) + t465)) * MDP(33) + (t564 * t557 - t588 * t788 - t615 * t751 - t679 * t714) * MDP(21) + (t401 * t537 - t403 * t536 - t404 * t683 - t405 * t500 + t414 * t441 - t418 * t442 + t433 * t437 - t436 * t438) * MDP(32) + (t437 * t536 + t438 * t537 + t441 * t500 - t442 * t683) * MDP(26) + MDP(6) * t747 + (pkin(2) * t615 * t706 + t628 * t556 - t594 * t716 + t646 * t663) * MDP(17) + 0.2e1 * t659 * MDP(4) * t706 + (-t594 * t556 + t615 * t663) * MDP(11) + (t437 * t537 - t441 * t683) * MDP(25) + ((t564 * t658 - t565 * t654) * t556 + (t658 * t714 - t764 + (-t654 * t564 - t565 * t658) * qJD(4)) * t615) * MDP(19) + t786 * t787 + (-pkin(7) * t747 + t656 * t704) * MDP(9) - MDP(7) * t749 + (pkin(7) * t749 + t659 * t704) * MDP(10) + (t495 * t754 + t565 * t677) * MDP(18); t662 + (t552 * t650 + (t594 * t726 - t650 * t707) * pkin(2) + t664) * MDP(17) + (-t644 * t714 + t684 * t654 - t689 * t564 + (-t642 * t723 + t798) * t588 + t680) * MDP(23) + (t403 * t513 + t401 * t512 + t413 * (t627 + t769) + (t681 + t689) * t452 + t744 * t418 + t743 * t414) * MDP(33) + (t437 * t512 - t438 * t513 - t500 * t744 - t683 * t743 + t667) * MDP(32) + (-t627 * t437 - t730 * t543 + t578 * t776 + t728 * t683 + t669) * MDP(31) + (t644 * t495 + t684 * t658 + t689 * t565 + (t642 * t724 + t796) * t588 + t687) * MDP(24) + (t696 * t543 + t627 * t438 + (-t608 * t721 + (-qJD(5) * t607 + t791) * t653 - t780) * t578 + t728 * t500 + t670) * MDP(30) + (t551 * t650 + (-t592 * t726 - t650 * t725) * pkin(2) + t676) * MDP(16) + t661 * t786 - t656 * MDP(4) * t746 + (MDP(9) * t656 * t661 + MDP(10) * t746) * pkin(1); t662 + (t695 * t543 + t645 * t438 + (-t631 * t721 + (-qJD(5) * t629 + t790) * t653 - t778) * t578 + t688 * t500 + t670) * MDP(30) + (pkin(3) * t714 - t697 * t588 + t549 * t564 + t529 * t757 + (-t588 * t723 - t751) * pkin(9) + t680) * MDP(23) + (t403 * t528 + t401 * t527 + t413 * (t645 + t769) + (t681 - t549) * t452 + t741 * t418 + t740 * t414) * MDP(33) + (t437 * t527 - t438 * t528 - t500 * t741 - t683 * t740 + t667) * MDP(32) + (-t645 * t437 - t729 * t543 + t578 * t777 + t688 * t683 + t669) * MDP(31) + (t548 * t650 + t664) * MDP(17) + (t549 * t650 + t676) * MDP(16) + (-pkin(3) * t495 + t734 * t588 - t549 * t565 + t529 * t756 + (t588 * t724 - t763) * pkin(9) + t687) * MDP(24); -t565 * t564 * MDP(18) + (-t564 ^ 2 + t565 ^ 2) * MDP(19) + (t495 - t760) * MDP(20) + (t714 + t759) * MDP(21) + t543 * MDP(22) + (t479 * t588 - t529 * t565 + t666) * MDP(23) + (t478 * t588 - t529 * t564 - t675) * MDP(24) + t500 * t720 + (t500 * t578 - t437) * MDP(27) + (-t438 + t782) * MDP(28) + (-t700 * t578 - t783 + (-t500 * t565 + t543 * t657 - t578 * t722) * pkin(4) + t665) * MDP(30) + (t742 * t578 + t496 * t500 + (-t543 * t653 - t565 * t683 - t578 * t721) * pkin(4) - t692) * MDP(31) + (-t414 * t500 + t418 * t683 + t419 * t683 + t420 * t500 + t437 * t643 + (-t438 * t653 + (-t500 * t657 + t653 * t683) * qJD(5)) * pkin(4)) * MDP(32) + (-pkin(5) * t784 + t401 * t643 - t414 * t419 - t418 * t420 + (t403 * t653 - t452 * t565 + (-t414 * t653 + t418 * t657) * qJD(5)) * pkin(4)) * MDP(33) + t793; t715 * MDP(27) + (-t698 + t782) * MDP(28) + (t430 * t578 + t702 - t783) * MDP(30) + (t429 * t578 - t692) * MDP(31) + t745 * MDP(33) * t418 + (t437 * MDP(32) + (t401 - t784) * MDP(33)) * pkin(5) + (t578 * MDP(27) + t496 * MDP(31) - MDP(32) * t745 + t720) * t500 + (-MDP(27) * t758 - t683 * MDP(28) - MDP(30) * t430) * qJD(5) + t793; (-t498 - t773) * MDP(32) + (t414 * t683 + t418 * t500 + t413) * MDP(33);];
tauc  = t1;
