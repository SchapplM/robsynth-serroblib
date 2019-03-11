% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR11_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR11_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:23
% EndTime: 2019-03-09 11:15:37
% DurationCPUTime: 11.64s
% Computational Cost: add. (6372->610), mult. (13420->777), div. (0->0), fcn. (9011->12), ass. (0->291)
t679 = cos(qJ(6));
t680 = cos(qJ(4));
t676 = sin(qJ(4));
t771 = qJD(2) * t676;
t681 = cos(qJ(2));
t774 = qJD(1) * t681;
t602 = t680 * t774 + t771;
t737 = t676 * t774;
t769 = qJD(2) * t680;
t604 = -t737 + t769;
t672 = sin(pkin(10));
t673 = cos(pkin(10));
t713 = -t602 * t673 - t604 * t672;
t536 = t602 * t672 - t604 * t673;
t675 = sin(qJ(6));
t809 = t536 * t675;
t481 = t679 * t713 + t809;
t677 = sin(qJ(2));
t775 = qJD(1) * t677;
t635 = qJD(4) + t775;
t627 = qJD(6) + t635;
t813 = t481 * t627;
t756 = qJD(1) * qJD(2);
t736 = t677 * t756;
t754 = qJDD(1) * t681;
t846 = t736 - t754;
t529 = -qJD(4) * t602 + t680 * qJDD(2) + t676 * t846;
t735 = t681 * t756;
t755 = qJDD(1) * t677;
t699 = t735 + t755;
t599 = qJDD(4) + t699;
t683 = -pkin(2) - pkin(8);
t658 = t677 * qJ(3);
t733 = -pkin(1) - t658;
t697 = t681 * t683 + t733;
t563 = t697 * qJD(1);
t650 = pkin(7) * t775;
t843 = qJD(3) + t650;
t758 = pkin(3) * t775 + t843;
t568 = qJD(2) * t683 + t758;
t507 = t563 * t680 + t568 * t676;
t634 = pkin(2) * t736;
t817 = qJ(3) * t681;
t716 = pkin(8) * t677 - t817;
t767 = qJD(3) * t677;
t691 = qJD(2) * t716 - t767;
t505 = qJD(1) * t691 + qJDD(1) * t697 + t634;
t633 = pkin(7) * t735;
t647 = pkin(7) * t755;
t734 = qJDD(3) + t633 + t647;
t533 = t699 * pkin(3) + qJDD(2) * t683 + t734;
t728 = -t505 * t676 + t680 * t533;
t689 = -qJD(4) * t507 + t728;
t443 = pkin(4) * t599 - qJ(5) * t529 - qJD(5) * t604 + t689;
t530 = -qJD(4) * t737 + qJDD(2) * t676 + (qJD(2) * qJD(4) - t846) * t680;
t764 = qJD(4) * t680;
t748 = -t680 * t505 - t676 * t533 - t568 * t764;
t765 = qJD(4) * t676;
t445 = -qJ(5) * t530 - qJD(5) * t602 - t563 * t765 - t748;
t431 = t673 * t443 - t445 * t672;
t475 = t529 * t673 - t530 * t672;
t429 = pkin(5) * t599 - pkin(9) * t475 + t431;
t432 = t672 * t443 + t673 * t445;
t474 = -t529 * t672 - t530 * t673;
t430 = pkin(9) * t474 + t432;
t506 = -t563 * t676 + t680 * t568;
t495 = -qJ(5) * t604 + t506;
t487 = pkin(4) * t635 + t495;
t496 = -qJ(5) * t602 + t507;
t804 = t673 * t496;
t453 = t672 * t487 + t804;
t841 = pkin(9) * t713;
t448 = t453 + t841;
t760 = qJD(6) * t675;
t447 = t448 * t760;
t651 = pkin(7) * t774;
t610 = pkin(3) * t774 + t651;
t669 = qJD(2) * qJ(3);
t583 = t669 + t610;
t541 = pkin(4) * t602 + qJD(5) + t583;
t488 = -pkin(5) * t713 + t541;
t657 = qJ(4) + pkin(10) + qJ(6);
t643 = sin(t657);
t644 = cos(t657);
t678 = sin(qJ(1));
t682 = cos(qJ(1));
t797 = t677 * t682;
t552 = t643 * t797 + t644 * t678;
t799 = t677 * t678;
t554 = -t643 * t799 + t644 * t682;
t666 = g(3) * t681;
t851 = g(1) * t552 - g(2) * t554 - t675 * t429 - t679 * t430 - t488 * t481 - t643 * t666 + t447;
t593 = qJDD(6) + t599;
t842 = -t679 * t536 + t675 * t713;
t850 = t593 * MDP(28) + (-t481 ^ 2 + t842 ^ 2) * MDP(25) - t481 * MDP(24) * t842;
t814 = t842 * t627;
t655 = pkin(2) * t775;
t574 = qJD(1) * t716 + t655;
t591 = t680 * t610;
t788 = qJ(5) - t683;
t801 = t676 * t677;
t848 = t574 * t676 - t591 - (pkin(4) * t681 - qJ(5) * t801) * qJD(1) - qJD(5) * t680 + t765 * t788;
t613 = t788 * t680;
t743 = t680 * t775;
t780 = t680 * t574 + t676 * t610;
t847 = qJ(5) * t743 + qJD(4) * t613 + qJD(5) * t676 + t780;
t711 = t672 * t680 + t673 * t676;
t781 = t635 * t711;
t551 = -t643 * t678 + t644 * t797;
t553 = t643 * t682 + t644 * t799;
t730 = t679 * t429 - t675 * t430;
t845 = -g(1) * t551 - g(2) * t553 - t488 * t842 + t644 * t666 + t730;
t829 = pkin(3) + pkin(7);
t844 = pkin(9) * t536;
t738 = t673 * t764;
t805 = t672 * t676;
t782 = -t672 * t765 + t673 * t743 - t775 * t805 + t738;
t785 = t672 * t847 + t673 * t848;
t784 = t672 * t848 - t673 * t847;
t838 = t781 * t679;
t663 = t681 * pkin(2);
t778 = t663 + t658;
t615 = -pkin(1) - t778;
t597 = -pkin(8) * t681 + t615;
t622 = t829 * t677;
t779 = t680 * t597 + t676 * t622;
t710 = -t673 * t680 + t805;
t535 = -t675 * t711 - t679 * t710;
t837 = pkin(4) * t764 + t843;
t646 = pkin(4) * t680 + pkin(3);
t763 = qJD(4) * t681;
t739 = t676 * t763;
t770 = qJD(2) * t677;
t836 = -pkin(4) * t739 + (-pkin(7) - t646) * t770;
t835 = t677 * t769 + t739;
t580 = t680 * t599;
t834 = -t635 * t765 + t580;
t721 = g(1) * t682 + g(2) * t678;
t833 = pkin(4) * t530 + qJDD(5);
t792 = t680 * t682;
t585 = -t676 * t678 + t677 * t792;
t795 = t678 * t680;
t587 = t676 * t682 + t677 * t795;
t793 = t680 * t681;
t832 = -g(1) * t585 - g(2) * t587 + g(3) * t793;
t831 = t583 * t635 + t683 * t599;
t729 = -t679 * t474 + t475 * t675;
t438 = qJD(6) * t842 + t729;
t489 = t672 * t496;
t452 = t673 * t487 - t489;
t745 = -g(1) * t797 - g(2) * t799 + t666;
t830 = -t431 * t710 + t432 * t711 - t452 * t781 + t453 * t782 + t745;
t827 = pkin(4) * t672;
t825 = g(1) * t678;
t821 = g(2) * t682;
t820 = g(3) * t677;
t659 = t676 * pkin(4);
t818 = pkin(7) * qJDD(2);
t816 = qJDD(2) * pkin(2);
t446 = pkin(5) * t635 + t452 + t844;
t815 = t446 * t679;
t812 = t529 * t680;
t534 = -t675 * t710 + t679 * t711;
t811 = t534 * t593;
t810 = t535 * t593;
t807 = t602 * t635;
t806 = t604 * t635;
t674 = -qJ(5) - pkin(8);
t803 = t674 * t681;
t802 = t676 * t599;
t800 = t676 * t681;
t798 = t677 * t680;
t685 = qJD(1) ^ 2;
t796 = t677 * t685;
t794 = t678 * t681;
t791 = t681 * t682;
t789 = qJ(3) + t659;
t654 = pkin(2) * t770;
t549 = t654 + t691;
t768 = qJD(2) * t681;
t611 = t829 * t768;
t592 = t680 * t611;
t731 = qJ(5) * t681 - t597;
t761 = qJD(5) * t681;
t464 = pkin(4) * t768 + t592 + t731 * t764 + (-qJ(5) * t770 - qJD(4) * t622 - t549 + t761) * t676;
t698 = t680 * t549 - t597 * t765 + t676 * t611 + t622 * t764;
t468 = qJ(5) * t835 - t680 * t761 + t698;
t440 = t672 * t464 + t673 * t468;
t787 = -qJD(6) * t534 - t675 * t782 - t838;
t786 = qJD(6) * t535 - t675 * t781 + t679 * t782;
t458 = t673 * t495 - t489;
t606 = t680 * t622;
t519 = pkin(4) * t677 + t676 * t731 + t606;
t524 = -qJ(5) * t793 + t779;
t470 = t672 * t519 + t673 * t524;
t720 = qJD(1) * t646;
t783 = pkin(5) * t782 + t677 * t720 + t837;
t612 = t788 * t676;
t543 = -t673 * t612 - t672 * t613;
t623 = t829 * t681;
t670 = t677 ^ 2;
t671 = t681 ^ 2;
t777 = t670 - t671;
t773 = qJD(2) * t602;
t772 = qJD(2) * t604;
t766 = qJD(4) * t563;
t762 = qJD(4) * t683;
t759 = qJD(6) * t679;
t753 = pkin(4) * t800;
t751 = t676 * t797;
t750 = t681 * t796;
t749 = t675 * t474 + t679 * t475 + t713 * t759;
t648 = pkin(7) * t754;
t667 = qJDD(2) * qJ(3);
t668 = qJD(2) * qJD(3);
t747 = t648 + t667 + t668;
t746 = pkin(4) * t793 + t623;
t744 = t829 * qJD(2);
t742 = t676 * t770;
t732 = -qJD(2) * pkin(2) + qJD(3);
t439 = t673 * t464 - t468 * t672;
t457 = -t495 * t672 - t804;
t469 = t673 * t519 - t524 * t672;
t542 = t612 * t672 - t673 * t613;
t726 = t682 * pkin(1) + pkin(2) * t791 + t678 * pkin(7) + qJ(3) * t797;
t725 = -t647 - t745;
t724 = pkin(3) * t754 + t747;
t723 = qJD(1) * t744;
t684 = qJD(2) ^ 2;
t722 = pkin(7) * t684 + t821;
t515 = -pkin(9) * t711 + t543;
t719 = pkin(5) * t774 - pkin(9) * t781 + qJD(6) * t515 - t785;
t514 = pkin(9) * t710 + t542;
t718 = pkin(9) * t782 - qJD(6) * t514 - t784;
t717 = qJD(6) * t710 - t782;
t434 = t446 * t675 + t448 * t679;
t571 = t710 * t681;
t572 = t711 * t681;
t714 = t679 * t571 + t572 * t675;
t517 = t571 * t675 - t572 * t679;
t614 = t650 + t732;
t621 = -t651 - t669;
t712 = t614 * t681 + t621 * t677;
t709 = t635 * t676;
t708 = t733 - t663;
t645 = pkin(4) * t673 + pkin(5);
t706 = t645 * t675 + t679 * t827;
t705 = t645 * t679 - t675 * t827;
t704 = -0.2e1 * pkin(1) * t756 - t818;
t584 = t708 * qJD(1);
t703 = t584 * t775 + qJDD(3) - t725;
t702 = -t635 * t764 - t802;
t700 = -qJ(3) * t768 - t767;
t437 = t536 * t760 + t749;
t696 = 0.2e1 * qJDD(1) * pkin(1) - t722;
t694 = t818 + (-qJD(1) * t615 - t584) * qJD(2);
t693 = -t681 * t721 - t820;
t540 = -t677 * t723 + t724;
t690 = t540 + t693;
t531 = qJD(1) * t700 + qJDD(1) * t708 + t634;
t576 = t654 + t700;
t688 = qJD(1) * t576 + qJDD(1) * t615 + t531 + t722;
t564 = pkin(7) * t736 - t747;
t573 = t734 - t816;
t687 = qJD(2) * t712 - t564 * t681 + t573 * t677;
t485 = t540 + t833;
t664 = t682 * pkin(7);
t639 = g(1) * t794;
t632 = qJ(3) * t791;
t630 = qJ(3) * t794;
t609 = t677 * t744;
t607 = -qJ(3) * t774 + t655;
t588 = -t676 * t799 + t792;
t586 = t751 + t795;
t562 = pkin(5) * t711 + t789;
t532 = -pkin(5) * t571 + t746;
t521 = -t672 * t835 - t673 * t742 + t681 * t738;
t520 = qJD(4) * t572 - t710 * t770;
t501 = pkin(4) * t604 - pkin(5) * t536;
t486 = -pkin(5) * t520 + t836;
t461 = pkin(9) * t571 + t470;
t460 = pkin(5) * t677 + pkin(9) * t572 + t469;
t455 = qJD(6) * t517 - t679 * t520 - t521 * t675;
t454 = qJD(6) * t714 + t520 * t675 - t521 * t679;
t451 = -pkin(5) * t474 + t485;
t450 = t458 + t844;
t449 = t457 - t841;
t436 = pkin(9) * t520 + t440;
t435 = pkin(5) * t768 + pkin(9) * t521 + t439;
t433 = -t448 * t675 + t815;
t1 = [(-t438 * t677 - t455 * t627 + t481 * t768 + t593 * t714) * MDP(27) + ((t435 * t679 - t436 * t675) * t627 + (t460 * t679 - t461 * t675) * t593 + t730 * t677 + t433 * t768 - t486 * t481 + t532 * t438 - t451 * t714 + t488 * t455 - g(1) * t554 - g(2) * t552 + ((-t460 * t675 - t461 * t679) * t627 - t434 * t677) * qJD(6)) * MDP(29) + (t437 * t714 - t438 * t517 + t454 * t481 - t455 * t842) * MDP(25) + (t677 * t704 + t681 * t696 + t639) * MDP(9) + (-g(2) * t791 + t431 * t572 + t432 * t571 + t439 * t536 + t440 * t713 + t452 * t521 + t453 * t520 - t469 * t475 + t470 * t474 + t639) * MDP(22) + (-t698 * t635 - t779 * t599 - t609 * t604 + t623 * t529 + g(1) * t587 - g(2) * t585 + ((qJD(2) * t583 + t766) * t676 + t748) * t677 + (-qJD(2) * t507 - t540 * t676 - t583 * t764) * t681) * MDP(21) + ((-t549 * t676 + t592) * t635 + (-t597 * t676 + t606) * t599 + t728 * t677 - t609 * t602 + t623 * t530 + t540 * t793 - g(1) * t588 - g(2) * t586 + (t506 * t681 - t583 * t798) * qJD(2) + (-t507 * t677 - t583 * t800 - t635 * t779) * qJD(4)) * MDP(20) + (qJDD(1) * t670 + 0.2e1 * t677 * t735) * MDP(4) + ((t635 * t769 - t530) * t677 + (-t773 - t834) * t681) * MDP(18) + (t432 * t470 + t453 * t440 + t431 * t469 + t452 * t439 + t485 * t746 - g(1) * (t646 * t682 + t664) - g(2) * (pkin(4) * t751 - t674 * t791 + t726) + (-g(1) * (-pkin(4) * t801 + t708 + t803) - g(2) * t646) * t678 + t836 * t541) * MDP(23) + (t599 * t677 + t635 * t768) * MDP(19) + (t593 * t677 + t627 * t768) * MDP(28) + ((t670 + t671) * qJDD(1) * pkin(7) + t687 - t721) * MDP(11) + (-t529 * t800 + (-t680 * t763 + t742) * t604) * MDP(15) + qJDD(1) * MDP(1) + (t677 * t694 + t681 * t688 - t639) * MDP(12) + (pkin(7) * t687 - g(1) * t664 - g(2) * t726 + t531 * t615 + t584 * t576 - t708 * t825) * MDP(14) + (t704 * t681 + (-t696 - t825) * t677) * MDP(10) + (t694 * t681 + (-t688 + t825) * t677) * MDP(13) + (-t821 + t825) * MDP(2) + ((t635 * t771 + t529) * t677 + (t702 + t772) * t681) * MDP(17) + 0.2e1 * (t677 * t754 - t756 * t777) * MDP(5) + ((-t602 * t676 + t604 * t680) * t770 + (-t812 + t530 * t676 + (t602 * t680 + t604 * t676) * qJD(4)) * t681) * MDP(16) + t721 * MDP(3) + (qJDD(2) * t677 + t681 * t684) * MDP(6) + (qJDD(2) * t681 - t677 * t684) * MDP(7) + (-t434 * t768 + g(1) * t553 - g(2) * t551 + t532 * t437 + t447 * t677 + t451 * t517 + t488 * t454 + t486 * t842 + (-(-qJD(6) * t461 + t435) * t627 - t460 * t593 - t429 * t677) * t675 + (-(qJD(6) * t460 + t436) * t627 - t461 * t593 - (qJD(6) * t446 + t430) * t677) * t679) * MDP(30) + (t437 * t677 + t454 * t627 + t517 * t593 + t768 * t842) * MDP(26) + (t437 * t517 + t454 * t842) * MDP(24); (pkin(1) * t796 + t725) * MDP(9) + ((t514 * t679 - t515 * t675) * t593 + t562 * t438 + t451 * t534 + (t675 * t718 - t679 * t719) * t627 + t786 * t488 - t783 * t481 + t693 * t643) * MDP(29) + (-t607 * MDP(12) - t635 * MDP(19) - t506 * MDP(20) + t507 * MDP(21) - MDP(26) * t842 - MDP(27) * t481 - t627 * MDP(28) - t433 * MDP(29) + t434 * MDP(30)) * t774 + (-t437 * t534 - t438 * t535 + t481 * t787 - t786 * t842) * MDP(25) + (qJ(3) * t530 - t591 * t635 + t758 * t602 + t831 * t680 + ((t574 - t762) * t635 + t690) * t676) * MDP(20) + (qJ(3) * t529 + t780 * t635 + t758 * t604 - t831 * t676 + (-t635 * t762 + t690) * t680) * MDP(21) + ((t602 * t681 - t635 * t798) * qJD(1) + t702) * MDP(18) + (-t564 * qJ(3) - t621 * qJD(3) - t573 * pkin(2) - t584 * t607 - g(1) * (-pkin(2) * t797 + t632) - g(2) * (-pkin(2) * t799 + t630) - g(3) * t778 - t712 * qJD(1) * pkin(7)) * MDP(14) + (t474 * t543 - t475 * t542 + t536 * t785 + t713 * t784 - t830) * MDP(22) + (-t627 * t786 - t811) * MDP(27) + ((-t604 * t681 - t635 * t801) * qJD(1) + t834) * MDP(17) + (t648 + 0.2e1 * t667 + 0.2e1 * t668 + (qJD(1) * t607 - g(3)) * t677 + (qJD(1) * t584 - t721) * t681) * MDP(13) + t777 * MDP(5) * t685 + qJDD(2) * MDP(8) + (t820 - t648 + (pkin(1) * t685 + t721) * t681) * MDP(10) + ((-pkin(2) * t677 + t817) * qJDD(1) + ((-t621 - t669) * t677 + (-t614 + t732) * t681) * qJD(1)) * MDP(11) + ((-t530 - t806) * t680 + (-t529 + t807) * t676) * MDP(16) - MDP(4) * t750 + (t627 * t787 + t810) * MDP(26) + (t703 - 0.2e1 * t816) * MDP(12) + (t432 * t543 + t431 * t542 + t485 * t789 - g(1) * (t682 * t753 + t632) - g(2) * (t678 * t753 + t630) - g(3) * (t778 - t803) + t837 * t541 + t784 * t453 + t785 * t452 + (-g(3) * t659 + t541 * t720 + t721 * (pkin(2) - t674)) * t677) * MDP(23) + (-t604 * t709 + t812) * MDP(15) + (-(t514 * t675 + t515 * t679) * t593 + t562 * t437 + t451 * t535 + (t675 * t719 + t679 * t718) * t627 + t787 * t488 + t783 * t842 + t693 * t644) * MDP(30) + (t437 * t535 + t787 * t842) * MDP(24) + MDP(7) * t754 + MDP(6) * t755; MDP(11) * t755 + (qJDD(2) + t750) * MDP(12) + (-t670 * t685 - t684) * MDP(13) + (qJD(2) * t621 + t633 + t703 - t816) * MDP(14) + (-t635 * t709 + t580 - t773) * MDP(20) + (-t635 ^ 2 * t680 - t772 - t802) * MDP(21) + (t474 * t711 + t475 * t710 - t536 * t781 + t713 * t782) * MDP(22) + (-t541 * qJD(2) + t830) * MDP(23) + (t810 + qJD(2) * t481 + (t675 * t717 - t711 * t759 - t838) * t627) * MDP(29) + (-t811 - qJD(2) * t842 + (t717 * t679 + (qJD(6) * t711 + t781) * t675) * t627) * MDP(30); t604 * t602 * MDP(15) + (-t602 ^ 2 + t604 ^ 2) * MDP(16) + (t529 + t807) * MDP(17) + (-t530 + t806) * MDP(18) + t599 * MDP(19) + (t507 * t635 - t583 * t604 + t689 + t832) * MDP(20) + (g(1) * t586 - g(2) * t588 + t506 * t635 + t583 * t602 + (t766 - t666) * t676 + t748) * MDP(21) + ((t474 * t672 - t475 * t673) * pkin(4) + (t452 - t458) * t713 + (-t453 - t457) * t536) * MDP(22) + (-t452 * t457 - t453 * t458 + (t431 * t673 + t432 * t672 - t541 * t604 + t832) * pkin(4)) * MDP(23) + (t437 - t813) * MDP(26) + (-t438 + t814) * MDP(27) + (t705 * t593 - (t449 * t679 - t450 * t675) * t627 + t501 * t481 + (-t627 * t706 - t434) * qJD(6) + t845) * MDP(29) + (-t706 * t593 + (t449 * t675 + t450 * t679) * t627 - t501 * t842 + (-t627 * t705 - t815) * qJD(6) + t851) * MDP(30) + t850; (-t536 ^ 2 - t713 ^ 2) * MDP(22) + (t438 + t814) * MDP(29) + (t437 + t813) * MDP(30) + (-g(1) * t791 - g(2) * t794 - t452 * t536 - t453 * t713 + t724 + (-g(3) - t723) * t677 + t833) * MDP(23); (t749 - t813) * MDP(26) + (-t729 + t814) * MDP(27) + (t434 * t627 + t845) * MDP(29) + (t433 * t627 + t851) * MDP(30) + (MDP(26) * t809 - MDP(27) * t842 - MDP(29) * t434 - MDP(30) * t815) * qJD(6) + t850;];
tau  = t1;
