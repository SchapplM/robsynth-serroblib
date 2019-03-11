% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:09:41
% EndTime: 2019-03-09 22:09:55
% DurationCPUTime: 9.67s
% Computational Cost: add. (11901->537), mult. (28422->713), div. (0->0), fcn. (21487->10), ass. (0->257)
t741 = cos(qJ(2));
t853 = -pkin(8) - pkin(7);
t713 = t853 * t741;
t707 = qJD(1) * t713;
t737 = sin(qJ(3));
t687 = t737 * t707;
t738 = sin(qJ(2));
t712 = t853 * t738;
t705 = qJD(1) * t712;
t852 = cos(qJ(3));
t638 = t705 * t852 + t687;
t788 = qJD(3) * t852;
t880 = -pkin(2) * t788 + t638;
t789 = qJD(1) * t852;
t806 = qJD(1) * t738;
t868 = -t737 * t806 + t741 * t789;
t675 = qJD(4) - t868;
t739 = cos(qJ(6));
t826 = t737 * t741;
t684 = -qJD(1) * t826 - t738 * t789;
t730 = qJD(2) + qJD(3);
t736 = sin(qJ(4));
t740 = cos(qJ(4));
t647 = -t684 * t736 - t740 * t730;
t733 = sin(pkin(11));
t734 = cos(pkin(11));
t765 = t684 * t740 - t730 * t736;
t766 = -t647 * t734 + t733 * t765;
t823 = t739 * t766;
t578 = t647 * t733 + t734 * t765;
t735 = sin(qJ(6));
t840 = t578 * t735;
t519 = t823 + t840;
t662 = qJD(6) + t675;
t842 = t519 * t662;
t624 = t868 * t730;
t802 = qJD(4) * t740;
t803 = qJD(4) * t736;
t572 = t740 * t624 + t684 * t803 + t730 * t802;
t703 = t738 * t852 + t826;
t642 = t730 * t703;
t625 = t642 * qJD(1);
t725 = -pkin(2) * t741 - pkin(1);
t711 = t725 * qJD(1);
t602 = -pkin(3) * t868 + pkin(9) * t684 + t711;
t688 = t852 * t707;
t844 = qJD(2) * pkin(2);
t689 = t705 + t844;
t631 = t737 * t689 - t688;
t610 = pkin(9) * t730 + t631;
t549 = t602 * t736 + t610 * t740;
t800 = qJD(1) * qJD(2);
t787 = t738 * t800;
t559 = pkin(2) * t787 + pkin(3) * t625 - pkin(9) * t624;
t554 = t740 * t559;
t793 = qJD(2) * t853;
t777 = qJD(1) * t793;
t691 = t738 * t777;
t692 = t741 * t777;
t805 = qJD(3) * t737;
t564 = t689 * t788 + t691 * t852 + t737 * t692 + t707 * t805;
t746 = -qJD(4) * t549 - t564 * t736 + t554;
t463 = pkin(4) * t625 - qJ(5) * t572 + qJD(5) * t765 + t746;
t573 = -qJD(4) * t765 + t624 * t736;
t751 = t736 * t559 + t740 * t564 + t602 * t802 - t610 * t803;
t465 = -qJ(5) * t573 - qJD(5) * t647 + t751;
t449 = t734 * t463 - t465 * t733;
t508 = t572 * t734 - t573 * t733;
t444 = pkin(5) * t625 - pkin(10) * t508 + t449;
t450 = t733 * t463 + t734 * t465;
t507 = -t572 * t733 - t573 * t734;
t445 = pkin(10) * t507 + t450;
t548 = t740 * t602 - t610 * t736;
t527 = qJ(5) * t765 + t548;
t514 = pkin(4) * t675 + t527;
t528 = -qJ(5) * t647 + t549;
t829 = t734 * t528;
t489 = t733 * t514 + t829;
t862 = pkin(10) * t766;
t468 = t489 + t862;
t801 = qJD(6) * t735;
t467 = t468 * t801;
t630 = t689 * t852 + t687;
t609 = -t730 * pkin(3) - t630;
t569 = t647 * pkin(4) + qJD(5) + t609;
t515 = -pkin(5) * t766 + t569;
t879 = -t735 * t444 - t739 * t445 - t515 * t519 + t467;
t864 = -t739 * t578 + t735 * t766;
t878 = t625 * MDP(31) + (-t519 ^ 2 + t864 ^ 2) * MDP(28) - t519 * MDP(27) * t864;
t698 = t733 * t740 + t734 * t736;
t870 = t675 * t698;
t763 = t733 * t736 - t734 * t740;
t869 = t675 * t763;
t626 = -pkin(3) * t684 - pkin(9) * t868;
t608 = pkin(2) * t806 + t626;
t877 = t736 * t608 + t740 * t880;
t833 = t868 * t736;
t876 = -qJ(5) * t833 - t740 * qJD(5);
t843 = t864 * t662;
t601 = t740 * t608;
t729 = t740 * qJ(5);
t775 = -t684 * pkin(4) - t729 * t868;
t723 = pkin(2) * t737 + pkin(9);
t820 = -qJ(5) - t723;
t780 = qJD(4) * t820;
t874 = t740 * t780 - t601 - t775 + (-qJD(5) + t880) * t736;
t781 = t740 * t626 - t630 * t736;
t845 = -qJ(5) - pkin(9);
t785 = qJD(4) * t845;
t873 = -t736 * qJD(5) + t740 * t785 - t775 - t781;
t872 = -t736 * t780 + t876 + t877;
t813 = t736 * t626 + t740 * t630;
t871 = -t736 * t785 + t813 + t876;
t783 = t739 * t444 - t735 * t445;
t867 = -t515 * t864 + t783;
t866 = pkin(10) * t578;
t633 = t698 * t739 - t735 * t763;
t819 = -qJD(6) * t633 + t869 * t735 - t739 * t870;
t764 = -t698 * t735 - t739 * t763;
t818 = -qJD(6) * t764 + t735 * t870 + t869 * t739;
t637 = t737 * t705 - t688;
t776 = pkin(2) * t805 - t637;
t755 = -t737 * t738 + t741 * t852;
t641 = t730 * t755;
t790 = t703 * t802;
t865 = t641 * t736 + t790;
t863 = -0.2e1 * t800;
t861 = MDP(5) * (t738 ^ 2 - t741 ^ 2);
t817 = t733 * t872 + t734 * t874;
t816 = t733 * t874 - t734 * t872;
t859 = t870 * pkin(10);
t815 = t733 * t871 + t734 * t873;
t814 = t733 * t873 - t734 * t871;
t726 = pkin(4) * t803;
t858 = pkin(5) * t870 + t726;
t659 = pkin(4) * t833;
t857 = -t659 + t776;
t856 = -t684 * pkin(5) - pkin(10) * t869;
t855 = t852 * t712 + t737 * t713;
t854 = qJD(1) * t703;
t782 = -t739 * t507 + t508 * t735;
t460 = qJD(6) * t864 + t782;
t851 = pkin(4) * t733;
t849 = pkin(10) * t698;
t847 = t763 * pkin(5);
t846 = t740 * pkin(4);
t841 = t572 * t736;
t839 = t609 * t868;
t837 = t625 * t740;
t835 = t647 * t675;
t834 = t765 * t675;
t832 = t703 * t736;
t831 = t703 * t740;
t521 = t733 * t528;
t827 = t736 * t625;
t742 = qJD(2) ^ 2;
t825 = t738 * t742;
t488 = t734 * t514 - t521;
t466 = pkin(5) * t675 + t488 + t866;
t824 = t739 * t466;
t651 = t737 * t712 - t713 * t852;
t644 = t740 * t651;
t822 = t741 * t742;
t743 = qJD(1) ^ 2;
t821 = t741 * t743;
t799 = t738 * t844;
t568 = pkin(3) * t642 - pkin(9) * t641 + t799;
t567 = t740 * t568;
t706 = t738 * t793;
t708 = t741 * t793;
t584 = qJD(3) * t855 + t852 * t706 + t737 * t708;
t629 = -pkin(3) * t755 - pkin(9) * t703 + t725;
t762 = -qJ(5) * t641 - qJD(5) * t703;
t475 = pkin(4) * t642 - t584 * t736 + t567 + t762 * t740 + (-t644 + (qJ(5) * t703 - t629) * t736) * qJD(4);
t796 = t736 * t568 + t740 * t584 + t629 * t802;
t485 = -qJ(5) * t790 + (-qJD(4) * t651 + t762) * t736 + t796;
t456 = t733 * t475 + t734 * t485;
t492 = t734 * t527 - t521;
t619 = t740 * t629;
t542 = -pkin(4) * t755 - t651 * t736 - t703 * t729 + t619;
t810 = t736 * t629 + t644;
t557 = -qJ(5) * t832 + t810;
t502 = t733 * t542 + t734 * t557;
t811 = t857 + t858;
t795 = t659 + t631;
t809 = -t795 + t858;
t693 = t820 * t736;
t694 = t723 * t740 + t729;
t622 = t733 * t693 + t734 * t694;
t709 = t845 * t736;
t710 = pkin(9) * t740 + t729;
t646 = t733 * t709 + t734 * t710;
t804 = qJD(4) * t703;
t797 = qJD(6) * t823 + t735 * t507 + t739 * t508;
t794 = -pkin(3) - t846;
t599 = t609 * t802;
t784 = pkin(1) * t863;
t455 = t734 * t475 - t485 * t733;
t491 = -t527 * t733 - t829;
t501 = t734 * t542 - t557 * t733;
t621 = t734 * t693 - t694 * t733;
t645 = t734 * t709 - t710 * t733;
t779 = t675 * t740;
t565 = t689 * t805 + t737 * t691 - t852 * t692 - t707 * t788;
t724 = -pkin(2) * t852 - pkin(3);
t774 = -t549 * t684 + t565 * t736 + t599;
t690 = t763 * pkin(10);
t589 = -t690 + t622;
t773 = qJD(6) * t589 - t817 + t856;
t606 = -t690 + t646;
t772 = qJD(6) * t606 - t815 + t856;
t588 = t621 - t849;
t771 = -qJD(6) * t588 - t816 + t859;
t605 = t645 - t849;
t770 = -qJD(6) * t605 - t814 + t859;
t454 = t735 * t466 + t739 * t468;
t768 = -t625 * t723 - t839;
t613 = t698 * t703;
t614 = t763 * t703;
t767 = -t739 * t613 + t614 * t735;
t556 = -t613 * t735 - t614 * t739;
t760 = pkin(4) * t832 - t855;
t721 = pkin(4) * t734 + pkin(5);
t759 = t721 * t735 + t739 * t851;
t758 = t721 * t739 - t735 * t851;
t756 = t548 * t684 - t565 * t740 + t609 * t803;
t754 = t641 * t740 - t703 * t803;
t509 = pkin(4) * t573 + t565;
t753 = t724 - t846;
t752 = t684 * t711 - t565;
t585 = t737 * t706 - t708 * t852 + t712 * t805 - t713 * t788;
t459 = t578 * t801 + t797;
t453 = -t468 * t735 + t824;
t479 = -pkin(5) * t507 + t509;
t750 = t453 * t684 - t479 * t764 - t515 * t819;
t749 = -t454 * t684 + t479 * t633 - t515 * t818;
t748 = pkin(4) * t865 + t585;
t747 = -t449 * t698 - t450 * t763 + t488 * t869 - t489 * t870;
t745 = -t711 * t868 - t564;
t744 = (t459 * t764 - t460 * t633 - t519 * t818 + t819 * t864) * MDP(28) + (t459 * t633 - t818 * t864) * MDP(27) + (t625 * t633 - t662 * t818 + t684 * t864) * MDP(29) + (t519 * t684 + t625 * t764 + t662 * t819) * MDP(30) + ((t572 - t835) * t740 + (-t573 + t834) * t736) * MDP(19) + (-t765 * t779 + t841) * MDP(18) + (-t675 ^ 2 * t736 - t647 * t684 + t837) * MDP(21) + (t675 * t779 - t684 * t765 + t827) * MDP(20) + t624 * MDP(13) + (t684 ^ 2 - t868 ^ 2) * MDP(12) + (MDP(11) * t868 + MDP(22) * t675 + MDP(31) * t662) * t684 + (-t868 * MDP(13) + (-t684 - t854) * MDP(14)) * t730;
t655 = t794 + t847;
t652 = t753 + t847;
t591 = t625 * t755;
t558 = t613 * pkin(5) + t760;
t547 = -pkin(4) * t765 - pkin(5) * t578;
t534 = t641 * t763 + t698 * t804;
t533 = -t641 * t698 + t763 * t804;
t494 = -t533 * pkin(5) + t748;
t493 = -pkin(10) * t613 + t502;
t490 = -pkin(5) * t755 + pkin(10) * t614 + t501;
t477 = t492 + t866;
t476 = t491 - t862;
t472 = qJD(6) * t556 - t739 * t533 - t534 * t735;
t471 = qJD(6) * t767 + t533 * t735 - t534 * t739;
t452 = pkin(10) * t533 + t456;
t451 = pkin(5) * t642 + pkin(10) * t534 + t455;
t1 = [-MDP(7) * t825 + (pkin(7) * t825 + t741 * t784) * MDP(10) + (-pkin(7) * t822 + t738 * t784) * MDP(9) + (t449 * t501 + t450 * t502 + t488 * t455 + t489 * t456 + t509 * t760 + t569 * t748) * MDP(26) + (-t572 * t755 + t625 * t831 - t642 * t765 + t675 * t754) * MDP(20) + (t642 * t675 - t591) * MDP(22) + (t642 * t662 - t591) * MDP(31) + (t573 * t755 - t642 * t647 - t675 * t865 - t703 * t827) * MDP(21) + (-t454 * t642 + t558 * t459 - t467 * t755 + t515 * t471 + t479 * t556 + t494 * t864 + (-(-qJD(6) * t493 + t451) * t662 - t490 * t625 + t444 * t755) * t735 + (-(qJD(6) * t490 + t452) * t662 - t493 * t625 + (qJD(6) * t466 + t445) * t755) * t739) * MDP(33) + (-t459 * t755 + t471 * t662 + t556 * t625 + t642 * t864) * MDP(29) + (t459 * t556 + t471 * t864) * MDP(27) + (-(-t651 * t803 + t796) * t675 - t810 * t625 + t751 * t755 - t549 * t642 - t585 * t765 - t855 * t572 + t565 * t831 + t754 * t609) * MDP(24) + ((-t651 * t802 + t567) * t675 + t619 * t625 - (-t610 * t802 + t554) * t755 + t548 * t642 + t585 * t647 - t855 * t573 + t703 * t599 + ((-qJD(4) * t629 - t584) * t675 - t651 * t625 - (-qJD(4) * t602 - t564) * t755 + t565 * t703 + t609 * t641) * t736) * MDP(23) + 0.2e1 * t741 * MDP(4) * t787 + MDP(6) * t822 + t861 * t863 + (t624 * t725 + t641 * t711 + (-t684 + t854) * t799) * MDP(17) + (t449 * t614 - t450 * t613 + t455 * t578 + t456 * t766 + t488 * t534 + t489 * t533 - t501 * t508 + t502 * t507) * MDP(25) + ((-t647 * t740 + t736 * t765) * t641 + (-t841 - t573 * t740 + (t647 * t736 + t740 * t765) * qJD(4)) * t703) * MDP(19) + (t572 * t831 - t754 * t765) * MDP(18) + (t624 * t703 - t641 * t684) * MDP(11) + (t460 * t755 - t472 * t662 + t519 * t642 + t625 * t767) * MDP(30) + ((t451 * t739 - t452 * t735) * t662 + (t490 * t739 - t493 * t735) * t625 - t783 * t755 + t453 * t642 - t494 * t519 + t558 * t460 - t479 * t767 + t515 * t472 + ((-t490 * t735 - t493 * t739) * t662 + t454 * t755) * qJD(6)) * MDP(32) + (t459 * t767 - t460 * t556 + t471 * t519 - t472 * t864) * MDP(28) + (t624 * t755 - t625 * t703 + t641 * t868 + t642 * t684) * MDP(12) + (t625 * t725 + t642 * t711 + (-qJD(1) * t755 - t868) * t799) * MDP(16) + (MDP(13) * t641 - MDP(14) * t642 - MDP(16) * t585 - MDP(17) * t584) * t730; t744 + (t507 * t622 - t508 * t621 + t578 * t817 + t766 * t816 + t747) * MDP(25) + (t450 * t622 + t449 * t621 + t509 * t753 + (t726 + t857) * t569 + t816 * t489 + t817 * t488) * MDP(26) + (t637 * t730 + (-t730 * t805 + t806 * t868) * pkin(2) + t752) * MDP(16) + ((t588 * t739 - t589 * t735) * t625 + t652 * t460 + (t735 * t771 - t739 * t773) * t662 - t811 * t519 + t750) * MDP(32) + (t724 * t572 + t768 * t740 - t776 * t765 + (t723 * t803 + t877) * t675 + t774) * MDP(24) + (t638 * t730 + (t684 * t806 - t730 * t788) * pkin(2) + t745) * MDP(17) + (-(t588 * t735 + t589 * t739) * t625 + t652 * t459 + (t735 * t773 + t739 * t771) * t662 + t811 * t864 + t749) * MDP(33) + (t724 * t573 + t768 * t736 + t776 * t647 + (-t723 * t802 + t736 * t880 - t601) * t675 + t756) * MDP(23) - t738 * MDP(4) * t821 + t743 * t861 + (MDP(9) * t738 * t743 + MDP(10) * t821) * pkin(1); t744 + (t507 * t646 - t508 * t645 + t578 * t815 + t766 * t814 + t747) * MDP(25) + (t450 * t646 + t449 * t645 + t509 * t794 + (-t795 + t726) * t569 + t814 * t489 + t815 * t488) * MDP(26) + (-pkin(3) * t573 - t781 * t675 - t631 * t647 - t609 * t833 + (-t675 * t802 - t827) * pkin(9) + t756) * MDP(23) + (-(t605 * t735 + t606 * t739) * t625 + t655 * t459 + (t735 * t772 + t739 * t770) * t662 + t809 * t864 + t749) * MDP(33) + (-pkin(3) * t572 + t813 * t675 + t631 * t765 - t740 * t839 + (t675 * t803 - t837) * pkin(9) + t774) * MDP(24) + (t630 * t730 + t745) * MDP(17) + (t631 * t730 + t752) * MDP(16) + ((t605 * t739 - t606 * t735) * t625 + t655 * t460 + (t735 * t770 - t739 * t772) * t662 - t809 * t519 + t750) * MDP(32); -t765 * t647 * MDP(18) + (-t647 ^ 2 + t765 ^ 2) * MDP(19) + (t572 + t835) * MDP(20) + (-t573 - t834) * MDP(21) + t625 * MDP(22) + (t549 * t675 + t609 * t765 + t746) * MDP(23) + (t548 * t675 + t609 * t647 - t751) * MDP(24) + ((t507 * t733 - t508 * t734) * pkin(4) + (t488 - t492) * t766 + (-t489 - t491) * t578) * MDP(25) + (-t488 * t491 - t489 * t492 + (t449 * t734 + t450 * t733 + t569 * t765) * pkin(4)) * MDP(26) + (t459 - t842) * MDP(29) + (-t460 + t843) * MDP(30) + (t758 * t625 - (t476 * t739 - t477 * t735) * t662 + t547 * t519 + (-t662 * t759 - t454) * qJD(6) + t867) * MDP(32) + (-t759 * t625 + (t476 * t735 + t477 * t739) * t662 - t547 * t864 + (-t662 * t758 - t824) * qJD(6) + t879) * MDP(33) + t878; (-t578 ^ 2 - t766 ^ 2) * MDP(25) + (-t488 * t578 - t489 * t766 + t509) * MDP(26) + (t460 + t843) * MDP(32) + (t459 + t842) * MDP(33); (t797 - t842) * MDP(29) + (-t782 + t843) * MDP(30) + (t454 * t662 + t867) * MDP(32) + (t453 * t662 + t879) * MDP(33) + (MDP(29) * t840 - MDP(30) * t864 - MDP(32) * t454 - MDP(33) * t824) * qJD(6) + t878;];
tauc  = t1;
