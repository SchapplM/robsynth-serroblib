% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRRP4
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
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:15:58
% EndTime: 2019-03-10 01:16:13
% DurationCPUTime: 9.81s
% Computational Cost: add. (13247->570), mult. (31122->715), div. (0->0), fcn. (22740->8), ass. (0->248)
t656 = cos(qJ(2));
t777 = -pkin(8) - pkin(7);
t629 = t777 * t656;
t619 = qJD(1) * t629;
t653 = sin(qJ(3));
t603 = t653 * t619;
t654 = sin(qJ(2));
t627 = t777 * t654;
t617 = qJD(1) * t627;
t775 = cos(qJ(3));
t562 = t617 * t775 + t603;
t707 = t775 * qJD(3);
t813 = -pkin(2) * t707 + t562;
t708 = qJD(1) * t775;
t730 = qJD(1) * t654;
t598 = t653 * t730 - t656 * t708;
t754 = t653 * t656;
t600 = -qJD(1) * t754 - t654 * t708;
t554 = -pkin(3) * t600 + pkin(9) * t598;
t534 = pkin(2) * t730 + t554;
t652 = sin(qJ(4));
t655 = cos(qJ(4));
t814 = -t655 * t534 + t813 * t652;
t812 = t652 * t534 + t655 * t813;
t648 = qJD(2) + qJD(3);
t571 = t652 * t600 + t648 * t655;
t572 = -t600 * t655 + t648 * t652;
t651 = sin(qJ(5));
t774 = cos(qJ(5));
t509 = -t774 * t571 + t572 * t651;
t683 = t651 * t571 + t572 * t774;
t811 = t509 * t683;
t762 = t598 * t655;
t696 = -t600 * pkin(4) + pkin(10) * t762;
t640 = pkin(2) * t653 + pkin(9);
t771 = -pkin(10) - t640;
t703 = qJD(4) * t771;
t810 = -t655 * t703 + t696 - t814;
t770 = qJD(2) * pkin(2);
t605 = t617 + t770;
t559 = t605 * t775 + t603;
t701 = t655 * t554 - t559 * t652;
t776 = -pkin(10) - pkin(9);
t717 = qJD(4) * t776;
t809 = -t655 * t717 + t696 + t701;
t763 = t598 * t652;
t722 = pkin(10) * t763;
t808 = -t652 * t703 + t722 + t812;
t735 = t652 * t554 + t655 * t559;
t807 = -t652 * t717 + t722 + t735;
t714 = t774 * t652;
t614 = t651 * t655 + t714;
t782 = qJD(4) + qJD(5);
t565 = t782 * t614;
t802 = t614 * t598 + t565;
t713 = t774 * t655;
t758 = t651 * t652;
t680 = t713 - t758;
t705 = t774 * qJD(5);
t784 = t774 * qJD(4) + t705;
t801 = t680 * t598 + t655 * t784 - t758 * t782;
t686 = -t653 * t654 + t656 * t775;
t668 = t686 * qJD(3);
t566 = qJD(2) * t686 + t668;
t663 = t566 * qJD(1);
t662 = t655 * t663;
t660 = qJD(4) * t571 + t662;
t727 = qJD(4) * t655;
t728 = qJD(4) * t652;
t719 = t600 * t727 - t648 * t728 - t652 * t663;
t726 = qJD(5) * t651;
t446 = -t571 * t705 + t572 * t726 - t651 * t719 - t774 * t660;
t594 = qJD(4) + t598;
t587 = qJD(5) + t594;
t429 = t509 * t587 - t446;
t447 = t571 * t726 + t572 * t705 + t651 * t660 - t774 * t719;
t778 = t683 ^ 2;
t615 = t654 * t775 + t754;
t567 = t648 * t615;
t805 = t567 * qJD(1);
t806 = t805 * MDP(29) + (t587 * t683 - t447) * MDP(28) + MDP(25) * t811 + (-t509 ^ 2 + t778) * MDP(26) + t429 * MDP(27);
t535 = -t648 * pkin(3) - t559;
t506 = -t571 * pkin(4) + t535;
t450 = t509 * pkin(5) - qJ(6) * t683 + t506;
t804 = t450 * t509;
t803 = t506 * t509;
t787 = (t728 + t763) * pkin(4);
t755 = t652 * t566;
t800 = t615 * t727 + t755;
t465 = pkin(5) * t683 + qJ(6) * t509;
t797 = t594 ^ 2;
t723 = qJD(1) * qJD(2);
t796 = -0.2e1 * t723;
t794 = MDP(5) * (t654 ^ 2 - t656 ^ 2);
t608 = t771 * t652;
t647 = t655 * pkin(10);
t760 = t640 * t655;
t609 = t647 + t760;
t682 = t608 * t774 - t651 * t609;
t793 = -qJD(5) * t682 + t810 * t651 + t808 * t774;
t553 = t651 * t608 + t609 * t774;
t792 = -qJD(5) * t553 + t808 * t651 - t810 * t774;
t644 = -pkin(2) * t656 - pkin(1);
t558 = -pkin(3) * t686 - pkin(9) * t615 + t644;
t547 = t655 * t558;
t576 = t653 * t627 - t629 * t775;
t750 = t655 * t615;
t481 = -pkin(4) * t686 - pkin(10) * t750 - t576 * t652 + t547;
t570 = t655 * t576;
t733 = t652 * t558 + t570;
t761 = t615 * t652;
t494 = -pkin(10) * t761 + t733;
t791 = t651 * t481 + t774 * t494;
t626 = t776 * t652;
t628 = pkin(9) * t655 + t647;
t681 = t626 * t774 - t651 * t628;
t790 = -qJD(5) * t681 + t809 * t651 + t807 * t774;
t575 = t651 * t626 + t628 * t774;
t789 = -qJD(5) * t575 + t807 * t651 - t809 * t774;
t604 = t775 * t619;
t561 = t653 * t617 - t604;
t729 = qJD(3) * t653;
t788 = -pkin(2) * t729 + t561;
t786 = -t572 * t728 + t660 * t655;
t785 = -t802 * pkin(5) + t801 * qJ(6) + qJD(6) * t614 - t787;
t783 = t775 * t627 + t653 * t629;
t549 = t805 * pkin(5);
t625 = t644 * qJD(1);
t532 = pkin(3) * t598 + pkin(9) * t600 + t625;
t495 = t805 * pkin(3) + (-pkin(9) * t668 + (t654 * pkin(2) - pkin(9) * t686) * qJD(2)) * qJD(1);
t492 = t655 * t495;
t560 = t653 * t605 - t604;
t536 = pkin(9) * t648 + t560;
t693 = -t536 * t727 + t492;
t718 = qJD(2) * t777;
t697 = qJD(1) * t718;
t606 = t654 * t697;
t607 = t656 * t697;
t498 = t605 * t707 + t606 * t775 + t653 * t607 + t619 * t729;
t757 = t652 * t498;
t420 = pkin(4) * t805 - pkin(10) * t660 - t532 * t728 + t693 - t757;
t676 = t652 * t495 + t655 * t498 + t532 * t727 - t536 * t728;
t423 = pkin(10) * t719 + t676;
t488 = t655 * t532 - t536 * t652;
t470 = -pkin(10) * t572 + t488;
t464 = pkin(4) * t594 + t470;
t489 = t652 * t532 + t655 * t536;
t471 = pkin(10) * t571 + t489;
t698 = -t774 * t420 + t651 * t423 + t464 * t726 + t471 * t705;
t409 = -t549 + t698;
t669 = t450 * t683 + t409;
t781 = -t506 * t683 - t698;
t720 = t654 * t770;
t505 = pkin(3) * t567 - pkin(9) * t566 + t720;
t501 = t655 * t505;
t618 = t654 * t718;
t620 = t656 * t718;
t516 = qJD(3) * t783 + t775 * t618 + t653 * t620;
t752 = t655 * t566;
t428 = -pkin(10) * t752 + pkin(4) * t567 - t516 * t652 + t501 + (-t570 + (pkin(10) * t615 - t558) * t652) * qJD(4);
t675 = t652 * t505 + t655 * t516 + t558 * t727 - t576 * t728;
t435 = -pkin(10) * t800 + t675;
t780 = -qJD(5) * t791 + t428 * t774 - t651 * t435;
t773 = t600 * pkin(5);
t772 = t655 * pkin(4);
t715 = t774 * t471;
t433 = t651 * t464 + t715;
t769 = t433 * t587;
t768 = t805 * t655;
t767 = t682 * t805;
t766 = t553 * t805;
t765 = t681 * t805;
t764 = t575 * t805;
t759 = t651 * t471;
t756 = t652 * t805;
t657 = qJD(2) ^ 2;
t753 = t654 * t657;
t751 = t655 * t571;
t749 = t656 * t657;
t658 = qJD(1) ^ 2;
t748 = t656 * t658;
t437 = t470 * t774 - t759;
t747 = -pkin(4) * t705 - qJD(6) + t437;
t592 = t600 * qJ(6);
t746 = -t592 + t793;
t745 = t773 + t792;
t744 = -t592 + t790;
t743 = t773 + t789;
t742 = t560 + t785;
t741 = t785 + t788;
t732 = t787 - t788;
t432 = t464 * t774 - t759;
t724 = qJD(6) - t432;
t721 = t775 * pkin(2);
t643 = -pkin(3) - t772;
t712 = t594 * t728;
t711 = t615 * t728;
t526 = t535 * t728;
t527 = t535 * t727;
t704 = t654 * t723;
t702 = pkin(1) * t796;
t699 = t651 * t420 + t774 * t423 + t464 * t705 - t471 * t726;
t499 = t605 * t729 + t653 * t606 - t775 * t607 - t619 * t707;
t642 = -t721 - pkin(3);
t695 = -t560 + t787;
t436 = t651 * t470 + t715;
t694 = pkin(4) * t726 - t436;
t533 = pkin(4) * t761 - t783;
t545 = t805 * qJ(6);
t580 = t587 * qJD(6);
t407 = t545 + t580 + t699;
t690 = t488 * t600 - t499 * t655 + t526;
t687 = t481 * t774 - t651 * t494;
t678 = t432 * t587 - t699;
t677 = t600 * t625 - t499;
t517 = t653 * t618 - t620 * t775 + t627 * t729 - t629 * t707;
t674 = t651 * t428 + t774 * t435 + t481 * t705 - t494 * t726;
t458 = -pkin(4) * t719 + t499;
t415 = t447 * pkin(5) + t446 * qJ(6) - qJD(6) * t683 + t458;
t430 = -t587 * pkin(5) + t724;
t673 = -t415 * t680 - t430 * t600 + t450 * t802;
t431 = t587 * qJ(6) + t433;
t672 = -t415 * t614 + t431 * t600 - t450 * t801;
t671 = t432 * t600 - t458 * t680 + t506 * t802;
t670 = -t433 * t600 + t458 * t614 + t506 * t801;
t555 = -pkin(5) * t680 - t614 * qJ(6) + t643;
t667 = -t489 * t600 + t499 * t652 + t535 * t762 + t527;
t475 = pkin(4) * t800 + t517;
t666 = t407 * t680 + t409 * t614 + t430 * t801 - t431 * t802;
t665 = t625 * t598 - t498;
t661 = (-t446 * t680 - t447 * t614 - t509 * t801 - t683 * t802) * MDP(26) + (-t446 * t614 + t683 * t801) * MDP(25) + (t571 * t727 - t572 * t763 + t598 * t751 + t652 * t719 + t786) * MDP(19) + (t587 * t801 + t600 * t683 + t614 * t805) * MDP(27) + (-t509 * t600 - t587 * t802 + t680 * t805) * MDP(28) + (t572 * t762 + (t652 * t571 + t572 * t655) * qJD(4) + t652 * t662) * MDP(18) + (t571 * t600 - t797 * t652 + t768) * MDP(21) + (t572 * t600 + t655 * t797 + t756) * MDP(20) + (t598 * t648 + t663) * MDP(13) + (-t600 * t648 - t805) * MDP(14) + (-t598 ^ 2 + t600 ^ 2) * MDP(12) + (-t598 * MDP(11) + MDP(22) * t594 + MDP(29) * t587) * t600;
t641 = -pkin(4) * t774 - pkin(5);
t637 = pkin(4) * t651 + qJ(6);
t624 = t642 - t772;
t543 = t680 * t615;
t542 = t614 * t615;
t539 = -t721 + t555;
t523 = t805 * t686;
t468 = t542 * pkin(5) - t543 * qJ(6) + t533;
t455 = pkin(4) * t572 + t465;
t454 = t566 * t714 - t651 * t711 - t726 * t761 + (t566 * t651 + t615 * t784) * t655;
t453 = t565 * t615 - t566 * t713 + t651 * t755;
t449 = pkin(5) * t686 - t687;
t448 = -qJ(6) * t686 + t791;
t416 = t454 * pkin(5) + t453 * qJ(6) - t543 * qJD(6) + t475;
t412 = -t567 * pkin(5) - t780;
t411 = qJ(6) * t567 - qJD(6) * t686 + t674;
t1 = [(t407 * t448 + t409 * t449 + t411 * t431 + t412 * t430 + t415 * t468 + t416 * t450) * MDP(35) + (-t489 * t567 + t499 * t750 + t517 * t572 - t526 * t615 + t535 * t752 - t594 * t675 - t660 * t783 + t676 * t686 - t733 * t805) * MDP(24) + ((-t576 * t727 + t501) * t594 + t547 * t805 - t693 * t686 + t488 * t567 - t517 * t571 + t783 * t719 + t615 * t527 + ((-qJD(4) * t558 - t516) * t594 - t576 * t805 - (-qJD(4) * t532 - t498) * t686 + t499 * t615 + t535 * t566) * t652) * MDP(23) + (t432 * t567 + t533 * t447 + t506 * t454 + t458 * t542 + t475 * t509 + t587 * t780 + t686 * t698 + t687 * t805) * MDP(30) + (t805 * t750 + t572 * t567 - t660 * t686 + (-t711 + t752) * t594) * MDP(20) + (t446 * t686 - t453 * t587 + t543 * t805 + t567 * t683) * MDP(27) + (-t407 * t686 + t411 * t587 - t415 * t543 - t416 * t683 + t431 * t567 + t446 * t468 + t448 * t805 + t450 * t453) * MDP(34) + (-t566 * t598 + t600 * t567 - t615 * t805 + t663 * t686) * MDP(12) + (t447 * t686 - t454 * t587 - t509 * t567 - t542 * t805) * MDP(28) + (t409 * t686 - t412 * t587 + t415 * t542 + t416 * t509 - t430 * t567 + t447 * t468 - t449 * t805 + t450 * t454) * MDP(32) + (t805 * t644 + t567 * t625 + (-qJD(1) * t686 + t598) * t720) * MDP(16) + (-t433 * t567 - t533 * t446 - t506 * t453 + t458 * t543 + t475 * t683 - t587 * t674 + t686 * t699 - t791 * t805) * MDP(31) + (t566 * t751 - t571 * t711 - t572 * t800 - t660 * t761 + t719 * t750) * MDP(19) + (t571 * t567 - t594 * t800 - t615 * t756 - t686 * t719) * MDP(21) + (t572 * t752 + t615 * t786) * MDP(18) + (MDP(13) * t566 - MDP(14) * t567 - MDP(16) * t517 - MDP(17) * t516) * t648 + (pkin(2) * t615 * t704 + t625 * t566 - t600 * t720 + t644 * t663) * MDP(17) + (-t446 * t543 - t453 * t683) * MDP(25) + (t446 * t542 - t447 * t543 + t453 * t509 - t454 * t683) * MDP(26) + (-t407 * t542 + t409 * t543 - t411 * t509 + t412 * t683 - t430 * t453 - t431 * t454 - t446 * t449 - t447 * t448) * MDP(33) + (-t600 * t566 + t615 * t663) * MDP(11) + 0.2e1 * t656 * MDP(4) * t704 - MDP(7) * t753 + (pkin(7) * t753 + t656 * t702) * MDP(10) + (-pkin(7) * t749 + t654 * t702) * MDP(9) + MDP(6) * t749 + t794 * t796 + (t567 * t594 - t523) * MDP(22) + (t567 * t587 - t523) * MDP(29); -t654 * MDP(4) * t748 + t661 + (-t788 * t572 + t812 * t594 + t640 * t712 + t642 * t660 - t805 * t760 + t667) * MDP(24) + (t561 * t648 + (-t598 * t730 - t648 * t729) * pkin(2) + t677) * MDP(16) + (t407 * t553 - t409 * t682 + t415 * t539 - t430 * t745 - t431 * t746 - t450 * t741) * MDP(35) + (t446 * t539 - t587 * t746 + t683 * t741 + t672 + t766) * MDP(34) + (t446 * t682 - t447 * t553 + t509 * t746 - t683 * t745 + t666) * MDP(33) + (-t624 * t446 + t587 * t793 + t732 * t683 + t670 - t766) * MDP(31) + (-t642 * t719 + (t535 * t598 - t640 * t805) * t652 + t788 * t571 + (-t640 * t727 + t814) * t594 + t690) * MDP(23) + t658 * t794 + (t562 * t648 + (t600 * t730 - t648 * t707) * pkin(2) + t665) * MDP(17) + (t624 * t447 + t732 * t509 + t587 * t792 + t671 + t767) * MDP(30) + (t447 * t539 - t509 * t741 + t587 * t745 + t673 + t767) * MDP(32) + (MDP(9) * t654 * t658 + MDP(10) * t748) * pkin(1); t661 + (t407 * t575 - t409 * t681 + t415 * t555 - t430 * t743 - t431 * t744 - t450 * t742) * MDP(35) + (-pkin(3) * t660 - t560 * t572 + t735 * t594 + t667 + (t712 - t768) * pkin(9)) * MDP(24) + (t643 * t447 + t695 * t509 + t587 * t789 + t671 + t765) * MDP(30) + (t559 * t648 + t665) * MDP(17) + (t560 * t648 + t677) * MDP(16) + (t446 * t555 - t587 * t744 + t683 * t742 + t672 + t764) * MDP(34) + (t446 * t681 - t447 * t575 + t509 * t744 - t683 * t743 + t666) * MDP(33) + (-t643 * t446 + t587 * t790 + t695 * t683 + t670 - t764) * MDP(31) + (pkin(3) * t719 - t701 * t594 + t560 * t571 + t535 * t763 + (-t594 * t727 - t756) * pkin(9) + t690) * MDP(23) + (t447 * t555 - t509 * t742 + t587 * t743 + t673 + t765) * MDP(32); -t572 * t571 * MDP(18) + (-t571 ^ 2 + t572 ^ 2) * MDP(19) + (-t571 * t594 + t660) * MDP(20) + (t572 * t594 + t719) * MDP(21) + t805 * MDP(22) + (-t535 * t572 + t492 - t757 + (-qJD(4) + t594) * t489) * MDP(23) + (t488 * t594 - t535 * t571 - t676) * MDP(24) + (t436 * t587 + (-t509 * t572 - t587 * t726 + t774 * t805) * pkin(4) + t781) * MDP(30) + (t437 * t587 + t803 + (-t572 * t683 - t587 * t705 - t651 * t805) * pkin(4) - t699) * MDP(31) + (-t455 * t509 - t587 * t694 - t641 * t805 - t669) * MDP(32) + (-t446 * t641 - t447 * t637 + (t431 + t694) * t683 + (t430 + t747) * t509) * MDP(33) + (t455 * t683 - t587 * t747 + t637 * t805 + t407 - t804) * MDP(34) + (t407 * t637 + t409 * t641 + t430 * t694 - t431 * t747 - t450 * t455) * MDP(35) + t806; (t769 + t781) * MDP(30) + (t678 + t803) * MDP(31) + (-t465 * t509 + t549 - t669 + t769) * MDP(32) + (pkin(5) * t446 - qJ(6) * t447 + (t431 - t433) * t683 + (t430 - t724) * t509) * MDP(33) + (t465 * t683 + 0.2e1 * t545 + 0.2e1 * t580 - t678 - t804) * MDP(34) + (-pkin(5) * t409 + qJ(6) * t407 - t430 * t433 + t431 * t724 - t450 * t465) * MDP(35) + t806; t429 * MDP(33) + (-t587 ^ 2 - t778) * MDP(34) + (-t431 * t587 + t669) * MDP(35) + (-t805 + t811) * MDP(32);];
tauc  = t1;
