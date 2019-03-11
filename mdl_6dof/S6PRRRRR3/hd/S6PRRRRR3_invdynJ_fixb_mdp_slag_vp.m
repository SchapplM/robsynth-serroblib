% Calculate vector of inverse dynamics joint torques for
% S6PRRRRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:51:56
% EndTime: 2019-03-09 00:52:19
% DurationCPUTime: 18.12s
% Computational Cost: add. (7601->673), mult. (17214->928), div. (0->0), fcn. (13695->18), ass. (0->273)
t690 = sin(qJ(4));
t695 = cos(qJ(4));
t779 = t695 * qJD(3);
t691 = sin(qJ(3));
t793 = qJD(2) * t691;
t631 = t690 * t793 - t779;
t790 = qJD(3) * t690;
t633 = t695 * t793 + t790;
t689 = sin(qJ(5));
t694 = cos(qJ(5));
t558 = t694 * t631 + t633 * t689;
t693 = cos(qJ(6));
t688 = sin(qJ(6));
t727 = t631 * t689 - t694 * t633;
t839 = t727 * t688;
t495 = -t693 * t558 + t839;
t696 = cos(qJ(3));
t676 = t696 * qJDD(2);
t777 = qJD(2) * qJD(3);
t853 = -t691 * t777 + t676;
t628 = qJDD(4) - t853;
t623 = qJDD(5) + t628;
t605 = qJDD(6) + t623;
t729 = t558 * t688 + t693 * t727;
t882 = t605 * MDP(30) + (-t495 ^ 2 + t729 ^ 2) * MDP(27) + t495 * MDP(26) * t729;
t785 = qJD(4) * t691;
t865 = -qJD(2) * t785 + qJDD(3);
t757 = t696 * t777;
t775 = qJDD(2) * t691;
t867 = t757 + t775;
t546 = qJD(4) * t779 + t690 * t865 + t695 * t867;
t792 = qJD(2) * t696;
t547 = t690 * (qJD(3) * (qJD(4) + t792) + t775) - t865 * t695;
t782 = qJD(5) * t694;
t783 = qJD(5) * t689;
t470 = t694 * t546 - t689 * t547 - t631 * t782 - t633 * t783;
t703 = qJD(5) * t727 - t546 * t689 - t694 * t547;
t780 = qJD(6) * t693;
t773 = t693 * t470 - t558 * t780 + t688 * t703;
t781 = qJD(6) * t688;
t443 = t727 * t781 + t773;
t663 = -qJD(4) + t792;
t656 = -qJD(5) + t663;
t748 = t470 * t688 - t693 * t703;
t704 = qJD(6) * t729 - t748;
t649 = -qJD(6) + t656;
t870 = t649 * t729;
t871 = t495 * t649;
t881 = t623 * MDP(23) + (-t558 ^ 2 + t727 ^ 2) * MDP(20) + (-t558 * t656 + t470) * MDP(21) + (t656 * t727 + t703) * MDP(22) - t558 * MDP(19) * t727 + (t704 + t870) * MDP(29) + (t443 + t871) * MDP(28) + t882;
t692 = sin(qJ(2));
t686 = sin(pkin(6));
t797 = qJD(1) * t686;
t644 = qJD(2) * pkin(8) + t692 * t797;
t687 = cos(pkin(6));
t796 = qJD(1) * t687;
t661 = t691 * t796;
t591 = t696 * t644 + t661;
t584 = qJD(3) * pkin(9) + t591;
t646 = -pkin(3) * t696 - pkin(9) * t691 - pkin(2);
t697 = cos(qJ(2));
t795 = qJD(1) * t697;
t770 = t686 * t795;
t593 = qJD(2) * t646 - t770;
t523 = -t584 * t690 + t695 * t593;
t499 = -pkin(10) * t633 + t523;
t482 = -pkin(4) * t663 + t499;
t524 = t584 * t695 + t593 * t690;
t500 = -pkin(10) * t631 + t524;
t489 = t694 * t500;
t464 = t482 * t689 + t489;
t873 = pkin(11) * t558;
t458 = t464 - t873;
t456 = t458 * t781;
t855 = -t691 * t644 + t696 * t796;
t583 = -qJD(3) * pkin(3) - t855;
t538 = pkin(4) * t631 + t583;
t484 = pkin(5) * t558 + t538;
t843 = cos(pkin(12));
t752 = t843 * t692;
t685 = sin(pkin(12));
t830 = t685 * t697;
t607 = t687 * t752 + t830;
t753 = t686 * t843;
t573 = t607 * t696 - t691 * t753;
t751 = t843 * t697;
t831 = t685 * t692;
t609 = -t687 * t831 + t751;
t575 = t685 * t686 * t691 + t609 * t696;
t606 = -t687 * t751 + t831;
t608 = t687 * t830 + t752;
t828 = t686 * t696;
t614 = t687 * t691 + t692 * t828;
t684 = qJ(4) + qJ(5);
t680 = qJ(6) + t684;
t669 = sin(t680);
t670 = cos(t680);
t827 = t686 * t697;
t864 = -t484 * t495 - g(1) * (-t575 * t670 - t608 * t669) - g(2) * (-t573 * t670 - t606 * t669) - g(3) * (-t614 * t670 + t669 * t827) + t456;
t738 = pkin(3) * t691 - pkin(9) * t696;
t642 = t738 * qJD(3);
t789 = qJD(3) * t691;
t818 = t696 * t697;
t845 = pkin(8) * t690;
t878 = (-t690 * t818 + t692 * t695) * t797 - t695 * t642 - t789 * t845;
t784 = qJD(4) * t695;
t877 = -(t690 * t692 + t695 * t818) * t797 + t690 * t642 + t646 * t784;
t778 = qJD(1) * qJD(2);
t596 = qJDD(2) * pkin(8) + (qJDD(1) * t692 + t697 * t778) * t686;
t776 = qJDD(1) * t687;
t755 = t691 * t776;
t508 = qJDD(3) * pkin(9) + qJD(3) * t855 + t596 * t696 + t755;
t759 = t692 * t778;
t733 = -qJDD(1) * t827 + t686 * t759;
t535 = qJD(2) * t642 + qJDD(2) * t646 + t733;
t531 = t695 * t535;
t705 = -qJD(4) * t524 - t690 * t508 + t531;
t449 = pkin(4) * t628 - pkin(10) * t546 + t705;
t772 = t695 * t508 + t690 * t535 + t593 * t784;
t786 = qJD(4) * t690;
t718 = t584 * t786 - t772;
t454 = -pkin(10) * t547 - t718;
t749 = t694 * t449 - t689 * t454;
t706 = -qJD(5) * t464 + t749;
t439 = pkin(5) * t623 - pkin(11) * t470 + t706;
t740 = -t689 * t449 - t694 * t454 - t482 * t782 + t500 * t783;
t440 = pkin(11) * t703 - t740;
t750 = t693 * t439 - t688 * t440;
t863 = t484 * t729 - g(1) * (-t575 * t669 + t608 * t670) - g(2) * (-t573 * t669 + t606 * t670) - g(3) * (-t614 * t669 - t670 * t827) + t750;
t819 = t695 * t696;
t665 = pkin(8) * t819;
t726 = pkin(4) * t691 - pkin(10) * t819;
t876 = t726 * qJD(3) + (-t665 + (pkin(10) * t691 - t646) * t690) * qJD(4) - t878;
t788 = qJD(3) * t696;
t763 = t690 * t788;
t866 = t691 * t784 + t763;
t875 = -t866 * pkin(10) + (-t691 * t779 - t696 * t786) * pkin(8) + t877;
t639 = t738 * qJD(2);
t621 = t695 * t639;
t847 = pkin(9) + pkin(10);
t771 = qJD(4) * t847;
t874 = qJD(2) * t726 - t690 * t855 + t695 * t771 + t621;
t765 = t690 * t792;
t806 = t690 * t639 + t695 * t855;
t868 = pkin(10) * t765 - t690 * t771 - t806;
t872 = pkin(11) * t727;
t634 = t689 * t690 - t694 * t695;
t719 = t634 * t696;
t851 = qJD(4) + qJD(5);
t808 = qJD(2) * t719 - t851 * t634;
t635 = t689 * t695 + t690 * t694;
t807 = (-t792 + t851) * t635;
t678 = sin(t684);
t679 = cos(t684);
t862 = t538 * t558 - g(1) * (-t575 * t679 - t608 * t678) - g(2) * (-t573 * t679 - t606 * t678) - g(3) * (-t614 * t679 + t678 * t827) + t740;
t861 = t538 * t727 - g(1) * (-t575 * t678 + t608 * t679) - g(2) * (-t573 * t678 + t606 * t679) - g(3) * (-t614 * t678 - t679 * t827) + t706;
t601 = t635 * t691;
t858 = t874 * t694;
t857 = t875 * t689 - t694 * t876;
t630 = t695 * t646;
t822 = t691 * t695;
t566 = -pkin(10) * t822 + t630 + (-pkin(4) - t845) * t696;
t801 = t690 * t646 + t665;
t823 = t690 * t691;
t580 = -pkin(10) * t823 + t801;
t856 = t566 * t782 - t580 * t783 + t689 * t876 + t875 * t694;
t739 = -t591 + (-t765 + t786) * pkin(4);
t809 = t689 * t566 + t694 * t580;
t651 = t847 * t690;
t652 = t847 * t695;
t802 = -t689 * t651 + t694 * t652;
t854 = -t651 * t782 - t652 * t783 - t689 * t874 + t694 * t868;
t852 = -t690 * t785 + t696 * t779;
t698 = qJD(3) ^ 2;
t737 = g(1) * t608 + g(2) * t606;
t850 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t698 + t686 * (-g(3) * t697 + t759) - t733 + t737;
t848 = -g(3) * t827 + qJD(4) * (pkin(8) * t663 + t584) + t737;
t846 = pkin(4) * t689;
t844 = qJD(2) * pkin(2);
t487 = t689 * t500;
t463 = t694 * t482 - t487;
t457 = t463 + t872;
t455 = -pkin(5) * t656 + t457;
t841 = t455 * t693;
t840 = t546 * t690;
t838 = t631 * t663;
t837 = t633 * t663;
t836 = t633 * t695;
t835 = t669 * t696;
t834 = t670 * t696;
t833 = t678 * t696;
t832 = t679 * t696;
t829 = t686 * t692;
t826 = t688 * t439;
t825 = t688 * t605;
t824 = t689 * t693;
t821 = t693 * t458;
t820 = t693 * t605;
t817 = qJDD(1) - g(3);
t510 = -qJD(3) * t719 - t601 * t851;
t816 = -pkin(5) * t789 + pkin(11) * t510 + qJD(5) * t809 + t857;
t511 = -t783 * t823 + (t822 * t851 + t763) * t694 + t852 * t689;
t815 = -pkin(11) * t511 + t856;
t564 = t693 * t634 + t635 * t688;
t814 = -qJD(6) * t564 - t688 * t807 + t693 * t808;
t565 = -t634 * t688 + t635 * t693;
t813 = qJD(6) * t565 + t688 * t808 + t693 * t807;
t812 = t694 * t499 - t487;
t810 = pkin(5) * t807 + t739;
t643 = pkin(4) * t823 + t691 * pkin(8);
t682 = t691 ^ 2;
t800 = -t696 ^ 2 + t682;
t794 = qJD(2) * t686;
t791 = qJD(3) * t631;
t787 = qJD(4) * t663;
t592 = pkin(4) * t866 + pkin(8) * t788;
t672 = -pkin(4) * t695 - pkin(3);
t768 = t691 * t795;
t767 = t692 * t794;
t766 = t697 * t794;
t764 = t663 * t779;
t747 = -t499 * t689 - t489;
t744 = t694 * t566 - t580 * t689;
t743 = -t694 * t651 - t652 * t689;
t742 = qJD(6) * t455 + t440;
t736 = g(1) * t609 + g(2) * t607;
t543 = -pkin(11) * t634 + t802;
t735 = pkin(5) * t793 + pkin(11) * t808 + qJD(5) * t802 + qJD(6) * t543 + t689 * t868 + t858;
t542 = -pkin(11) * t635 + t743;
t734 = -pkin(11) * t807 + qJD(6) * t542 + t854;
t442 = t688 * t455 + t821;
t602 = t634 * t691;
t477 = -pkin(5) * t696 + pkin(11) * t602 + t744;
t478 = -pkin(11) * t601 + t809;
t732 = t477 * t688 + t478 * t693;
t578 = -t614 * t690 - t695 * t827;
t723 = -t614 * t695 + t690 * t827;
t503 = t578 * t694 + t689 * t723;
t504 = t578 * t689 - t694 * t723;
t731 = t503 * t693 - t504 * t688;
t730 = t503 * t688 + t504 * t693;
t536 = t693 * t601 - t602 * t688;
t537 = -t601 * t688 - t602 * t693;
t613 = -t687 * t696 + t691 * t829;
t721 = t628 * t690 - t663 * t784;
t720 = t628 * t695 + t663 * t786;
t715 = g(1) * (-t609 * t691 + t685 * t828) + g(2) * (-t607 * t691 - t696 * t753) - g(3) * t613;
t714 = qJD(3) * t661 + t691 * t596 + t644 * t788 - t696 * t776;
t712 = -g(3) * t829 - t736;
t710 = -pkin(9) * t628 - t583 * t663;
t509 = -qJDD(3) * pkin(3) + t714;
t702 = pkin(9) * t787 - t509 - t715;
t474 = pkin(4) * t547 + t509;
t645 = -t770 - t844;
t701 = -pkin(8) * qJDD(3) + (t645 + t770 - t844) * qJD(3);
t699 = qJD(2) ^ 2;
t671 = pkin(4) * t694 + pkin(5);
t599 = pkin(5) * t634 + t672;
t577 = qJD(3) * t614 + t691 * t766;
t576 = -qJD(3) * t613 + t696 * t766;
t568 = pkin(5) * t601 + t643;
t532 = pkin(4) * t633 - pkin(5) * t727;
t486 = qJD(4) * t578 + t576 * t695 + t690 * t767;
t485 = qJD(4) * t723 - t576 * t690 + t695 * t767;
t479 = pkin(5) * t511 + t592;
t462 = qJD(6) * t537 + t510 * t688 + t693 * t511;
t461 = -qJD(6) * t536 + t510 * t693 - t511 * t688;
t460 = t812 + t872;
t459 = t747 + t873;
t452 = -pkin(5) * t703 + t474;
t451 = -qJD(5) * t504 + t485 * t694 - t486 * t689;
t450 = qJD(5) * t503 + t485 * t689 + t486 * t694;
t441 = -t458 * t688 + t841;
t1 = [t817 * MDP(1) + (-qJD(3) * t577 - qJDD(3) * t613) * MDP(10) + (-qJD(3) * t576 - qJDD(3) * t614) * MDP(11) + (-t485 * t663 + t547 * t613 + t577 * t631 + t578 * t628) * MDP(17) + (t486 * t663 + t546 * t613 + t577 * t633 + t628 * t723) * MDP(18) + (-t451 * t656 + t503 * t623 + t558 * t577 - t613 * t703) * MDP(24) + (t450 * t656 + t470 * t613 - t504 * t623 - t577 * t727) * MDP(25) + (-(-qJD(6) * t730 - t450 * t688 + t451 * t693) * t649 + t731 * t605 - t577 * t495 - t613 * t704) * MDP(31) + ((qJD(6) * t731 + t450 * t693 + t451 * t688) * t649 - t730 * t605 - t577 * t729 + t613 * t443) * MDP(32) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t696 + MDP(11) * t691 - MDP(3)) * t699) * t692 + (t853 * MDP(10) - MDP(11) * t867 + qJDD(2) * MDP(3) - t699 * MDP(4)) * t697) * t686; 0.2e1 * (t676 * t691 - t777 * t800) * MDP(6) + ((t663 * t790 + t547) * t696 + (-t721 - t791) * t691) * MDP(15) + (-t628 * t696 - t663 * t789) * MDP(16) + (-t623 * t696 - t656 * t789) * MDP(23) + (-t605 * t696 - t649 * t789) * MDP(30) + (t462 * t649 + t495 * t789 - t536 * t605 - t696 * t704) * MDP(29) + ((t477 * t693 - t478 * t688) * t605 - t750 * t696 + t441 * t789 - t479 * t495 - t568 * t704 + t452 * t536 + t484 * t462 - g(1) * (-t608 * t834 + t609 * t669) - g(2) * (-t606 * t834 + t607 * t669) + (t688 * t815 + t693 * t816) * t649 + (t442 * t696 + t649 * t732) * qJD(6) + (t495 * t768 - g(3) * (t669 * t692 + t670 * t818)) * t686) * MDP(31) + (-t443 * t536 + t461 * t495 + t462 * t729 + t537 * t704) * MDP(27) + ((-t631 * t695 - t633 * t690) * t788 + (-t840 - t547 * t695 + (t631 * t690 - t836) * qJD(4)) * t691) * MDP(13) + ((-t546 - t764) * t696 + (qJD(3) * t633 + t720) * t691) * MDP(14) + (qJDD(2) * t682 + 0.2e1 * t691 * t757) * MDP(5) + (-t817 * t829 + t736) * MDP(4) + (t546 * t822 + t633 * t852) * MDP(12) + (t817 * t827 + t737) * MDP(3) + (t701 * t691 + t696 * t850) * MDP(10) + (-t691 * t850 + t701 * t696) * MDP(11) + qJDD(2) * MDP(2) + (t511 * t656 - t558 * t789 - t601 * t623 - t696 * t703) * MDP(22) + (-t470 * t601 - t510 * t558 + t511 * t727 - t602 * t703) * MDP(20) + (t744 * t623 - t749 * t696 + t463 * t789 + t592 * t558 - t643 * t703 + t474 * t601 + t538 * t511 - g(1) * (-t608 * t832 + t609 * t678) - g(2) * (-t606 * t832 + t607 * t678) + t857 * t656 + (t464 * t696 + t656 * t809) * qJD(5) + (-t558 * t768 - g(3) * (t678 * t692 + t679 * t818)) * t686) * MDP(24) + (qJDD(3) * t691 + t696 * t698) * MDP(7) + (qJDD(3) * t696 - t691 * t698) * MDP(8) + (-t801 * t628 + t877 * t663 + t712 * t695 + ((pkin(8) * t633 + t583 * t695) * qJD(3) - t848 * t690 + t772) * t696 + (-t633 * t770 - t583 * t786 - t524 * qJD(3) + t509 * t695 + (t546 - t764) * pkin(8)) * t691) * MDP(18) + (t630 * t628 + t878 * t663 + (t646 * t787 + t712) * t690 + (pkin(8) * t791 - t531 + (-pkin(8) * t628 + qJD(3) * t583 + qJD(4) * t593 + t508) * t690 + t848 * t695) * t696 + (pkin(8) * t547 + qJD(3) * t523 + t509 * t690 + t583 * t784 - t631 * t770) * t691) * MDP(17) + (-t443 * t696 - t461 * t649 + t537 * t605 - t729 * t789) * MDP(28) + (t443 * t537 - t461 * t729) * MDP(26) + (-t732 * t605 + (t742 * t693 - t456 + t826) * t696 - t442 * t789 - t479 * t729 + t568 * t443 + t452 * t537 + t484 * t461 - g(1) * (t608 * t835 + t609 * t670) - g(2) * (t606 * t835 + t607 * t670) + ((qJD(6) * t477 + t815) * t693 + (-qJD(6) * t478 - t816) * t688) * t649 + (t729 * t768 - g(3) * (-t669 * t818 + t670 * t692)) * t686) * MDP(32) + (-t470 * t696 - t510 * t656 - t602 * t623 - t727 * t789) * MDP(21) + (-t470 * t602 - t510 * t727) * MDP(19) + (-t809 * t623 - t740 * t696 - t464 * t789 - t592 * t727 + t643 * t470 - t474 * t602 + t538 * t510 - g(1) * (t608 * t833 + t609 * t679) - g(2) * (t606 * t833 + t607 * t679) + t856 * t656 + (t727 * t768 - g(3) * (-t678 * t818 + t679 * t692)) * t686) * MDP(25); (qJD(3) * t591 - t714 - t715) * MDP(10) + ((t542 * t693 - t543 * t688) * t605 - t599 * t704 + t452 * t564 + (t688 * t734 + t693 * t735) * t649 - t810 * t495 + t813 * t484 - t715 * t670) * MDP(31) + (-t443 * t564 + t495 * t814 + t565 * t704 + t729 * t813) * MDP(27) + (-MDP(10) * t645 + t663 * MDP(16) - t523 * MDP(17) + MDP(18) * t524 + MDP(21) * t727 + t558 * MDP(22) + t656 * MDP(23) - t463 * MDP(24) + t464 * MDP(25) + MDP(28) * t729 - MDP(29) * t495 + t649 * MDP(30) - t441 * MDP(31) + t442 * MDP(32)) * t793 + (-pkin(3) * t546 - t591 * t633 - t663 * t806 - t690 * t702 + t695 * t710) * MDP(18) + (-t663 * t836 + t840) * MDP(12) + (-MDP(5) * t691 * t696 + MDP(6) * t800) * t699 + (-t623 * t634 + t656 * t807) * MDP(22) + (t623 * t635 - t656 * t808) * MDP(21) + (-t564 * t605 + t649 * t813) * MDP(29) + ((t546 + t838) * t695 + (-t547 + t837) * t690) * MDP(13) + (t565 * t605 - t649 * t814) * MDP(28) + MDP(8) * t676 + MDP(7) * t775 + (-t755 + g(1) * t575 + g(2) * t573 + g(3) * t614 + (-qJD(2) * t645 - t596) * t696) * MDP(11) + qJDD(3) * MDP(9) + (t743 * t623 - t672 * t703 + t474 * t634 + (t652 * t782 + (-qJD(5) * t651 + t868) * t689 + t858) * t656 + t739 * t558 + t807 * t538 - t715 * t679) * MDP(24) + (-t470 * t634 - t558 * t808 + t635 * t703 + t727 * t807) * MDP(20) + ((-t663 * t690 * t696 + t631 * t691) * qJD(2) + t720) * MDP(15) + ((-t633 * t691 + t663 * t819) * qJD(2) + t721) * MDP(14) + (-(t542 * t688 + t543 * t693) * t605 + t599 * t443 + t452 * t565 + (-t688 * t735 + t693 * t734) * t649 - t810 * t729 + t814 * t484 + t715 * t669) * MDP(32) + (t443 * t565 - t729 * t814) * MDP(26) + (t470 * t635 - t727 * t808) * MDP(19) + (t672 * t470 + t474 * t635 + t808 * t538 - t802 * t623 + t656 * t854 + t715 * t678 - t727 * t739) * MDP(25) + (-pkin(3) * t547 - t591 * t631 + t621 * t663 + (-t663 * t855 + t710) * t690 + t702 * t695) * MDP(17); (t546 - t838) * MDP(14) + (t671 * t820 + (t459 * t693 - t460 * t688) * t649 + t532 * t495 + (-t689 * t825 - (-t688 * t694 - t824) * t649 * qJD(5)) * pkin(4) + (-(-pkin(4) * t824 - t671 * t688) * t649 - t442) * qJD(6) + t863) * MDP(31) + (-t547 - t837) * MDP(15) + (-t524 * t663 - t583 * t633 - g(1) * (-t575 * t690 + t608 * t695) - g(2) * (-t573 * t690 + t606 * t695) - g(3) * t578 + t705) * MDP(17) + t628 * MDP(16) + (t532 * t729 + (-t671 * t605 - t439 + (-t459 + (-qJD(5) - qJD(6)) * t846) * t649) * t688 + (-t605 * t846 + (pkin(4) * t782 + qJD(6) * t671 - t460) * t649 - t742) * t693 + t864) * MDP(32) + (-t631 ^ 2 + t633 ^ 2) * MDP(13) + t633 * t631 * MDP(12) + (t747 * t656 + (-t558 * t633 + t623 * t694 + t656 * t783) * pkin(4) + t861) * MDP(24) + (-t812 * t656 + (-t623 * t689 + t633 * t727 + t656 * t782) * pkin(4) + t862) * MDP(25) + (-t523 * t663 + t583 * t631 - g(1) * (-t575 * t695 - t608 * t690) - g(2) * (-t573 * t695 - t606 * t690) - g(3) * t723 + t718) * MDP(18) + t881; (-t464 * t656 + t861) * MDP(24) + (-t463 * t656 + t862) * MDP(25) + ((-t457 * t688 - t821) * t649 - t442 * qJD(6) + (-t495 * t727 + t649 * t781 + t820) * pkin(5) + t863) * MDP(31) + ((t458 * t649 - t439) * t688 + (-t457 * t649 - t742) * t693 + (t649 * t780 - t727 * t729 - t825) * pkin(5) + t864) * MDP(32) + t881; (t773 + t871) * MDP(28) + (-t748 + t870) * MDP(29) + (-t442 * t649 + t863) * MDP(31) + (-t693 * t440 - t441 * t649 - t826 + t864) * MDP(32) + (MDP(28) * t839 + MDP(29) * t729 - MDP(31) * t442 - MDP(32) * t841) * qJD(6) + t882;];
tau  = t1;
