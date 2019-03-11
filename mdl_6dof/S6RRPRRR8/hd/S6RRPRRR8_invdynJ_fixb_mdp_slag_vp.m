% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR8
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:05:57
% EndTime: 2019-03-09 14:06:22
% DurationCPUTime: 18.72s
% Computational Cost: add. (13594->699), mult. (31351->909), div. (0->0), fcn. (24783->18), ass. (0->279)
t741 = sin(pkin(11));
t746 = sin(qJ(2));
t833 = qJD(1) * t746;
t812 = t741 * t833;
t742 = cos(pkin(11));
t820 = t742 * qJD(2);
t685 = t812 - t820;
t811 = t742 * t833;
t831 = qJD(2) * t741;
t687 = t811 + t831;
t745 = sin(qJ(4));
t750 = cos(qJ(4));
t614 = t685 * t745 - t687 * t750;
t615 = t750 * t685 + t687 * t745;
t744 = sin(qJ(5));
t749 = cos(qJ(5));
t558 = t614 * t744 - t749 * t615;
t748 = cos(qJ(6));
t743 = sin(qJ(6));
t878 = t614 * t749 + t615 * t744;
t863 = t878 * t743;
t499 = t558 * t748 + t863;
t751 = cos(qJ(2));
t735 = t751 * qJDD(1);
t818 = qJD(1) * qJD(2);
t773 = t746 * t818 - t735;
t692 = qJDD(4) + t773;
t679 = qJDD(5) + t692;
t673 = qJDD(6) + t679;
t785 = -t558 * t743 + t748 * t878;
t902 = t673 * MDP(33) + (-t499 ^ 2 + t785 ^ 2) * MDP(30) + t499 * MDP(29) * t785;
t787 = pkin(2) * t746 - qJ(3) * t751;
t696 = t787 * qJD(1);
t676 = t741 * t696;
t857 = t742 * t746;
t858 = t741 * t751;
t778 = -pkin(7) * t857 - pkin(8) * t858;
t621 = qJD(1) * t778 + t676;
t901 = qJD(3) * t742 - t621;
t731 = t742 * qJDD(2);
t808 = t751 * t818;
t817 = qJDD(1) * t746;
t774 = t808 + t817;
t642 = t741 * t774 - t731;
t816 = qJDD(2) * t741;
t643 = t742 * t774 + t816;
t825 = qJD(4) * t750;
t827 = qJD(4) * t745;
t543 = -t745 * t642 + t750 * t643 - t685 * t825 - t687 * t827;
t544 = -qJD(4) * t614 + t750 * t642 + t643 * t745;
t823 = qJD(5) * t749;
t824 = qJD(5) * t744;
t482 = t749 * t543 - t744 * t544 + t614 * t824 - t615 * t823;
t756 = qJD(5) * t878 - t543 * t744 - t749 * t544;
t821 = qJD(6) * t748;
t813 = t748 * t482 + t558 * t821 + t743 * t756;
t822 = qJD(6) * t743;
t457 = t822 * t878 + t813;
t804 = t743 * t482 - t748 * t756;
t757 = qJD(6) * t785 - t804;
t832 = qJD(1) * t751;
t720 = -qJD(4) + t832;
t713 = -qJD(5) + t720;
t707 = -qJD(6) + t713;
t885 = t707 * t785;
t893 = t713 * t878;
t895 = t558 * t713;
t899 = t499 * t707;
t900 = t679 * MDP(26) + (-t558 ^ 2 + t878 ^ 2) * MDP(23) + (t482 + t895) * MDP(24) + (t756 + t893) * MDP(25) + t558 * MDP(22) * t878 + (t757 + t885) * MDP(32) + (t457 + t899) * MDP(31) + t902;
t860 = t741 * t745;
t693 = -t750 * t742 + t860;
t775 = t693 * t751;
t839 = qJD(1) * t775 - t693 * qJD(4);
t694 = t741 * t750 + t742 * t745;
t776 = t694 * t751;
t838 = -qJD(1) * t776 + t694 * qJD(4);
t780 = pkin(2) * t751 + qJ(3) * t746 + pkin(1);
t678 = t780 * qJD(1);
t729 = pkin(7) * t832;
t708 = qJD(2) * qJ(3) + t729;
t623 = -t742 * t678 - t708 * t741;
t815 = pkin(3) * t832;
t576 = -pkin(8) * t687 + t623 - t815;
t624 = -t741 * t678 + t742 * t708;
t580 = -pkin(8) * t685 + t624;
t521 = t750 * t576 - t580 * t745;
t508 = pkin(9) * t614 + t521;
t503 = -pkin(4) * t720 + t508;
t522 = t576 * t745 + t580 * t750;
t509 = -pkin(9) * t615 + t522;
t507 = t749 * t509;
t476 = t503 * t744 + t507;
t890 = pkin(10) * t558;
t468 = t476 + t890;
t466 = t468 * t822;
t700 = -qJD(2) * pkin(2) + pkin(7) * t833 + qJD(3);
t633 = pkin(3) * t685 + t700;
t571 = pkin(4) * t615 + t633;
t511 = -pkin(5) * t558 + t571;
t738 = pkin(11) + qJ(4);
t736 = qJ(5) + t738;
t725 = qJ(6) + t736;
t715 = sin(t725);
t716 = cos(t725);
t752 = cos(qJ(1));
t747 = sin(qJ(1));
t854 = t747 * t751;
t637 = t715 * t752 - t716 * t854;
t852 = t751 * t752;
t639 = t715 * t747 + t716 * t852;
t867 = g(3) * t746;
t882 = g(1) * t639 - g(2) * t637 - t511 * t499 + t716 * t867 + t466;
t634 = pkin(7) * t812 + t742 * t696;
t856 = t742 * t751;
t781 = pkin(3) * t746 - pkin(8) * t856;
t602 = qJD(1) * t781 + t634;
t865 = pkin(8) + qJ(3);
t705 = t865 * t741;
t706 = t865 * t742;
t782 = qJD(3) * t741 + qJD(4) * t706;
t894 = -t705 * t825 + t901 * t750 + (-t602 - t782) * t745;
t636 = t715 * t854 + t716 * t752;
t638 = -t715 * t852 + t716 * t747;
t669 = qJD(2) * t787 - qJD(3) * t746;
t613 = qJD(1) * t669 - qJDD(1) * t780;
t656 = -pkin(7) * t773 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t567 = t742 * t613 - t656 * t741;
t529 = pkin(3) * t773 - pkin(8) * t643 + t567;
t568 = t741 * t613 + t742 * t656;
t546 = -pkin(8) * t642 + t568;
t801 = t750 * t529 - t546 * t745;
t759 = -qJD(4) * t522 + t801;
t461 = pkin(4) * t692 - pkin(9) * t543 + t759;
t772 = t745 * t529 + t750 * t546 + t576 * t825 - t580 * t827;
t463 = -pkin(9) * t544 + t772;
t805 = t749 * t461 - t744 * t463;
t758 = -qJD(5) * t476 + t805;
t451 = pkin(5) * t679 - pkin(10) * t482 + t758;
t792 = -t744 * t461 - t749 * t463 - t503 * t823 + t509 * t824;
t452 = pkin(10) * t756 - t792;
t806 = t748 * t451 - t743 * t452;
t881 = -g(1) * t638 + g(2) * t636 + t511 * t785 + t715 * t867 + t806;
t594 = t750 * t602;
t837 = -t745 * t705 + t750 * t706;
t891 = pkin(4) * t833 + pkin(9) * t839 + t694 * qJD(3) + qJD(4) * t837 - t621 * t745 + t594;
t883 = -pkin(9) * t838 + t894;
t727 = pkin(7) * t817;
t667 = -qJDD(2) * pkin(2) + pkin(7) * t808 + qJDD(3) + t727;
t791 = g(1) * t752 + g(2) * t747;
t866 = g(3) * t751;
t767 = t746 * t791 - t866;
t874 = t667 - t767;
t889 = pkin(10) * t878;
t887 = t614 * t720;
t886 = t615 * t720;
t619 = t749 * t693 + t694 * t744;
t845 = -qJD(5) * t619 - t744 * t838 + t749 * t839;
t620 = -t693 * t744 + t694 * t749;
t844 = qJD(5) * t620 + t744 * t839 + t749 * t838;
t722 = sin(t736);
t723 = cos(t736);
t645 = t722 * t752 - t723 * t854;
t647 = t722 * t747 + t723 * t852;
t880 = g(1) * t647 - g(2) * t645 - t558 * t571 + t723 * t867 + t792;
t644 = t722 * t854 + t723 * t752;
t646 = -t722 * t852 + t723 * t747;
t879 = -g(1) * t646 + g(2) * t644 + t571 * t878 + t722 * t867 + t758;
t875 = t891 * t749;
t665 = t693 * t746;
t684 = t742 * t780;
t622 = -pkin(8) * t857 - t684 + (-pkin(7) * t741 - pkin(3)) * t751;
t649 = pkin(7) * t856 - t741 * t780;
t859 = t741 * t746;
t629 = -pkin(8) * t859 + t649;
t796 = t750 * t622 - t629 * t745;
t536 = -pkin(4) * t751 + pkin(9) * t665 + t796;
t664 = t694 * t746;
t841 = t745 * t622 + t750 * t629;
t542 = -pkin(9) * t664 + t841;
t846 = t744 * t536 + t749 * t542;
t794 = -t750 * t705 - t706 * t745;
t595 = -pkin(9) * t694 + t794;
t596 = -pkin(9) * t693 + t837;
t843 = t744 * t595 + t749 * t596;
t680 = t741 * t815 + t729;
t807 = pkin(4) * t838 - t680;
t873 = t595 * t823 - t596 * t824 - t744 * t891 + t749 * t883;
t790 = g(1) * t747 - g(2) * t752;
t871 = pkin(7) * t685;
t870 = pkin(7) * t687;
t505 = t744 * t509;
t475 = t749 * t503 - t505;
t467 = t475 + t889;
t465 = -pkin(5) * t713 + t467;
t864 = t465 * t748;
t862 = t673 * t744;
t726 = pkin(4) * t749 + pkin(5);
t861 = t726 * t673;
t855 = t744 * t748;
t853 = t748 * t468;
t562 = t748 * t619 + t620 * t743;
t851 = -qJD(6) * t562 - t743 * t844 + t748 * t845;
t563 = -t619 * t743 + t620 * t748;
t850 = qJD(6) * t563 + t743 * t845 + t748 * t844;
t849 = t749 * t508 - t505;
t847 = pkin(5) * t844 + t807;
t830 = qJD(2) * t746;
t814 = pkin(7) * t830;
t627 = t742 * t669 + t741 * t814;
t829 = qJD(2) * t751;
t681 = (pkin(3) * t741 + pkin(7)) * t829;
t697 = pkin(3) * t859 + t746 * pkin(7);
t739 = t746 ^ 2;
t836 = -t751 ^ 2 + t739;
t826 = qJD(4) * t746;
t819 = qJD(3) - t700;
t724 = -pkin(3) * t742 - pkin(2);
t810 = qJ(3) * t735;
t597 = -qJD(2) * t775 - t694 * t826;
t591 = qJD(2) * t781 + t627;
t654 = t741 * t669;
t604 = qJD(2) * t778 + t654;
t798 = t750 * t591 - t604 * t745;
t490 = pkin(4) * t830 - pkin(9) * t597 - qJD(4) * t841 + t798;
t598 = qJD(2) * t776 + t825 * t857 - t826 * t860;
t771 = t745 * t591 + t750 * t604 + t622 * t825 - t629 * t827;
t492 = -pkin(9) * t598 + t771;
t803 = t749 * t490 - t492 * t744;
t802 = -t508 * t744 - t507;
t800 = t749 * t536 - t542 * t744;
t797 = t749 * t595 - t596 * t744;
t793 = qJD(6) * t465 + t452;
t572 = pkin(4) * t598 + t681;
t626 = pkin(4) * t664 + t697;
t515 = -pkin(10) * t620 + t797;
t789 = -pkin(10) * t844 + qJD(6) * t515 + t873;
t516 = -pkin(10) * t619 + t843;
t788 = pkin(5) * t833 + pkin(10) * t845 + qJD(5) * t843 + qJD(6) * t516 + t744 * t883 + t875;
t454 = t743 * t465 + t853;
t589 = t749 * t664 - t665 * t744;
t590 = -t664 * t744 - t665 * t749;
t533 = t748 * t589 + t590 * t743;
t534 = -t589 * t743 + t590 * t748;
t653 = pkin(4) * t693 + t724;
t779 = -0.2e1 * pkin(1) * t818 - pkin(7) * qJDD(2);
t770 = t744 * t490 + t749 * t492 + t536 * t823 - t542 * t824;
t754 = qJD(1) ^ 2;
t768 = pkin(1) * t754 + t791;
t592 = pkin(3) * t642 + t667;
t765 = -t751 * t791 - t867;
t753 = qJD(2) ^ 2;
t763 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t753 + t790;
t514 = pkin(4) * t544 + t592;
t734 = cos(t738);
t733 = sin(t738);
t660 = t733 * t747 + t734 * t852;
t659 = -t733 * t852 + t734 * t747;
t658 = t733 * t752 - t734 * t854;
t657 = t733 * t854 + t734 * t752;
t648 = -pkin(7) * t858 - t684;
t635 = -pkin(7) * t811 + t676;
t628 = -t742 * t814 + t654;
t575 = pkin(5) * t619 + t653;
t552 = pkin(5) * t589 + t626;
t517 = -pkin(4) * t614 - pkin(5) * t878;
t513 = qJD(5) * t590 + t597 * t744 + t749 * t598;
t512 = -qJD(5) * t589 + t597 * t749 - t598 * t744;
t493 = pkin(5) * t513 + t572;
t485 = -pkin(10) * t589 + t846;
t484 = -pkin(5) * t751 - pkin(10) * t590 + t800;
t472 = t849 + t889;
t471 = t802 - t890;
t470 = qJD(6) * t534 + t512 * t743 + t748 * t513;
t469 = -qJD(6) * t533 + t512 * t748 - t513 * t743;
t464 = -pkin(5) * t756 + t514;
t456 = -pkin(10) * t513 + t770;
t455 = pkin(5) * t830 - pkin(10) * t512 - qJD(5) * t846 + t803;
t453 = -t468 * t743 + t864;
t1 = [(t482 * t590 - t512 * t878) * MDP(22) + (-g(1) * t644 - g(2) * t646 - t476 * t830 + t626 * t482 + t571 * t512 + t514 * t590 - t572 * t878 - t679 * t846 + t713 * t770 - t751 * t792) * MDP(28) + (-t482 * t751 - t512 * t713 + t590 * t679 - t830 * t878) * MDP(24) + (-t798 * t720 + t796 * t692 - t801 * t751 + t521 * t830 + t681 * t615 + t697 * t544 + t592 * t664 + t633 * t598 - g(1) * t658 - g(2) * t660 + (t522 * t751 + t720 * t841) * qJD(4)) * MDP(20) + (-t482 * t589 + t512 * t558 + t513 * t878 + t590 * t756) * MDP(23) + (-t803 * t713 + t800 * t679 - t805 * t751 + t475 * t830 - t572 * t558 - t626 * t756 + t514 * t589 + t571 * t513 - g(1) * t645 - g(2) * t647 + (t476 * t751 + t713 * t846) * qJD(5)) * MDP(27) + (t513 * t713 + t558 * t830 - t589 * t679 - t751 * t756) * MDP(25) + (t567 * t648 + t568 * t649 + t623 * t627 + t624 * t628 + (t667 * t746 + t700 * t829 - t791) * pkin(7) + t790 * t780) * MDP(14) + qJDD(1) * MDP(1) + (-t454 * t830 - g(1) * t636 - g(2) * t638 + t552 * t457 + t464 * t534 - t466 * t751 + t511 * t469 - t493 * t785 + ((-qJD(6) * t485 + t455) * t707 - t484 * t673 + t451 * t751) * t743 + ((qJD(6) * t484 + t456) * t707 - t485 * t673 + t793 * t751) * t748) * MDP(35) + (-t457 * t751 - t469 * t707 + t534 * t673 - t785 * t830) * MDP(31) + (t457 * t534 - t469 * t785) * MDP(29) + (-t457 * t533 + t469 * t499 + t470 * t785 + t534 * t757) * MDP(30) + (t470 * t707 + t499 * t830 - t533 * t673 - t751 * t757) * MDP(32) + (-(t455 * t748 - t456 * t743) * t707 + (t484 * t748 - t485 * t743) * t673 - t806 * t751 + t453 * t830 - t493 * t499 - t552 * t757 + t464 * t533 + t511 * t470 - g(1) * t637 - g(2) * t639 + (-(-t484 * t743 - t485 * t748) * t707 + t454 * t751) * qJD(6)) * MDP(34) + (qJDD(2) * t751 - t746 * t753) * MDP(7) + (qJDD(2) * t746 + t751 * t753) * MDP(6) + (-t791 * t742 + (pkin(7) * t643 + t667 * t742 + (-qJD(1) * t649 - t624) * qJD(2)) * t746 + (t628 * qJD(1) + t649 * qJDD(1) + t568 - t790 * t741 + (t700 * t742 + t870) * qJD(2)) * t751) * MDP(12) + (-t791 * t741 + (pkin(7) * t642 + t667 * t741 + (qJD(1) * t648 + t623) * qJD(2)) * t746 + (-t627 * qJD(1) - t648 * qJDD(1) - t567 + t790 * t742 + (t700 * t741 + t871) * qJD(2)) * t751) * MDP(11) + 0.2e1 * (t735 * t746 - t818 * t836) * MDP(5) + (-t692 * t751 - t720 * t830) * MDP(19) + (-t679 * t751 - t713 * t830) * MDP(26) + (-t673 * t751 - t707 * t830) * MDP(33) + (t544 * t751 + t598 * t720 - t615 * t830 - t664 * t692) * MDP(18) + (-t627 * t687 - t628 * t685 - t642 * t649 - t643 * t648 + (-t623 * t742 - t624 * t741) * t829 + (-t567 * t742 - t568 * t741 + t790) * t746) * MDP(13) + (qJDD(1) * t739 + 0.2e1 * t746 * t808) * MDP(4) + (-t543 * t664 + t544 * t665 - t597 * t615 + t598 * t614) * MDP(16) + (-t543 * t665 - t597 * t614) * MDP(15) + (-g(1) * t657 - g(2) * t659 - t522 * t830 + t697 * t543 - t592 * t665 + t633 * t597 - t614 * t681 - t692 * t841 + t720 * t771 + t751 * t772) * MDP(21) + (-t543 * t751 - t597 * t720 - t614 * t830 - t665 * t692) * MDP(17) + (t746 * t779 + t751 * t763) * MDP(9) + (-t746 * t763 + t751 * t779) * MDP(10) + t790 * MDP(2) + t791 * MDP(3); (t653 * t482 + t514 * t620 + t845 * t571 - t843 * t679 + t713 * t873 - t722 * t767 - t807 * t878) * MDP(28) + (t482 * t620 - t845 * t878) * MDP(22) + (-(t515 * t743 + t516 * t748) * t673 + t575 * t457 + t464 * t563 + (-t743 * t788 + t748 * t789) * t707 + t851 * t511 - t847 * t785 - t767 * t715) * MDP(35) + (t620 * t679 - t713 * t845) * MDP(24) + (t563 * t673 - t707 * t851) * MDP(31) + (t794 * t692 + t724 * t544 + t592 * t693 - t680 * t615 + (t594 + t782 * t750 + (-qJD(4) * t705 + t901) * t745) * t720 + t838 * t633 + t767 * t734) * MDP(20) + (-t562 * t673 + t707 * t850) * MDP(32) + (-t619 * t679 + t713 * t844) * MDP(25) + (t742 * t810 - pkin(2) * t643 + t874 * t741 + ((-qJ(3) * t820 + t624) * t746 + (t742 * t819 - t635 - t870) * t751) * qJD(1)) * MDP(12) + (t724 * t543 + t592 * t694 + t680 * t614 + t839 * t633 - t837 * t692 + t720 * t894 - t767 * t733) * MDP(21) + (-t482 * t619 + t558 * t845 + t620 * t756 + t844 * t878) * MDP(23) + (t797 * t679 - t653 * t756 + t514 * t619 + (t596 * t823 + (qJD(5) * t595 + t883) * t744 + t875) * t713 + t844 * t571 - t807 * t558 + t767 * t723) * MDP(27) + (-MDP(4) * t746 * t751 + MDP(5) * t836) * t754 + qJDD(2) * MDP(8) + (t457 * t563 - t785 * t851) * MDP(29) + ((t515 * t748 - t516 * t743) * t673 - t575 * t757 + t464 * t562 + (t743 * t789 + t748 * t788) * t707 + t850 * t511 - t847 * t499 + t767 * t716) * MDP(34) + (-t457 * t562 + t499 * t851 + t563 * t757 + t785 * t850) * MDP(30) + (MDP(17) * t614 + t615 * MDP(18) + t720 * MDP(19) - t521 * MDP(20) + t522 * MDP(21) + MDP(24) * t878 - MDP(25) * t558 + t713 * MDP(26) - t475 * MDP(27) + t476 * MDP(28) + MDP(31) * t785 - MDP(32) * t499 + t707 * MDP(33) - t453 * MDP(34) + t454 * MDP(35)) * t833 + (t741 * t810 - pkin(2) * t642 - t874 * t742 + ((-qJ(3) * t831 - t623) * t746 + (t741 * t819 + t634 - t871) * t751) * qJD(1)) * MDP(11) + (-t700 * t729 - t623 * t634 - t624 * t635 + (-t623 * t741 + t624 * t742) * qJD(3) - t874 * pkin(2) + (-t567 * t741 + t568 * t742 + t765) * qJ(3)) * MDP(14) + MDP(6) * t817 + MDP(7) * t735 + (t746 * t768 - t727 - t866) * MDP(9) + (t867 + (-pkin(7) * qJDD(1) + t768) * t751) * MDP(10) + (t692 * t694 - t720 * t839) * MDP(17) + (t634 * t687 + t635 * t685 + (-qJ(3) * t642 - qJD(3) * t685 + t623 * t832 + t568) * t742 + (qJ(3) * t643 + qJD(3) * t687 + t624 * t832 - t567) * t741 + t765) * MDP(13) + (-t692 * t693 + t720 * t838) * MDP(18) + (t543 * t694 - t614 * t839) * MDP(15) + (-t543 * t693 - t544 * t694 + t614 * t838 - t615 * t839) * MDP(16); (t741 * t817 - t731 + (-t687 + t831) * t832) * MDP(11) + (t742 * t817 + t816 + (t685 + t820) * t832) * MDP(12) + (-t685 ^ 2 - t687 ^ 2) * MDP(13) + (t623 * t687 + t624 * t685 + t874) * MDP(14) + (t544 + t887) * MDP(20) + (t543 + t886) * MDP(21) + (-t756 + t893) * MDP(27) + (t482 - t895) * MDP(28) + (-t757 + t885) * MDP(34) + (t457 - t899) * MDP(35); -t614 * t615 * MDP(15) + (t517 * t785 + (-t861 - t451 + (-t471 + (-qJD(5) - qJD(6)) * t744 * pkin(4)) * t707) * t743 + (-pkin(4) * t862 + (pkin(4) * t823 + qJD(6) * t726 - t472) * t707 - t793) * t748 + t882) * MDP(35) + t692 * MDP(19) + (g(1) * t660 - g(2) * t658 - t521 * t720 + t615 * t633 + t734 * t867 - t772) * MDP(21) + (t748 * t861 + (t471 * t748 - t472 * t743) * t707 + t517 * t499 + (-t743 * t862 - (-t743 * t749 - t855) * t707 * qJD(5)) * pkin(4) + (-(-pkin(4) * t855 - t726 * t743) * t707 - t454) * qJD(6) + t881) * MDP(34) + (t614 ^ 2 - t615 ^ 2) * MDP(16) + (-t544 + t887) * MDP(18) + (-g(1) * t659 + g(2) * t657 - t522 * t720 + t614 * t633 + t733 * t867 + t759) * MDP(20) + (t802 * t713 + (-t558 * t614 + t679 * t749 + t713 * t824) * pkin(4) + t879) * MDP(27) + (-t849 * t713 + (-t614 * t878 - t679 * t744 + t713 * t823) * pkin(4) + t880) * MDP(28) + (t543 - t886) * MDP(17) + t900; (-t476 * t713 + t879) * MDP(27) + (-t475 * t713 + t880) * MDP(28) + ((-t467 * t743 - t853) * t707 - t454 * qJD(6) + (-t499 * t878 + t673 * t748 + t707 * t822) * pkin(5) + t881) * MDP(34) + ((t468 * t707 - t451) * t743 + (-t467 * t707 - t793) * t748 + (-t673 * t743 + t707 * t821 - t785 * t878) * pkin(5) + t882) * MDP(35) + t900; (t813 + t899) * MDP(31) + (-t804 + t885) * MDP(32) + (-t454 * t707 + t881) * MDP(34) + (-t743 * t451 - t748 * t452 - t453 * t707 + t882) * MDP(35) + (MDP(31) * t863 + MDP(32) * t785 - MDP(34) * t454 - MDP(35) * t864) * qJD(6) + t902;];
tau  = t1;
