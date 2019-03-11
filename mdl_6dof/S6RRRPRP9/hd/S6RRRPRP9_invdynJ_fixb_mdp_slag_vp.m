% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP9_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:27:38
% EndTime: 2019-03-09 17:27:58
% DurationCPUTime: 15.35s
% Computational Cost: add. (7510->772), mult. (15719->909), div. (0->0), fcn. (10606->8), ass. (0->292)
t723 = sin(qJ(3));
t727 = cos(qJ(3));
t828 = t727 * qJD(2);
t724 = sin(qJ(2));
t844 = qJD(1) * t724;
t633 = t723 * t844 - t828;
t813 = t727 * t844;
t841 = qJD(2) * t723;
t635 = t813 + t841;
t721 = qJD(2) * pkin(2);
t655 = pkin(7) * t844 - t721;
t542 = pkin(3) * t633 - qJ(4) * t635 + t655;
t524 = -pkin(4) * t633 - t542;
t722 = sin(qJ(5));
t726 = cos(qJ(5));
t561 = t633 * t722 + t635 * t726;
t773 = -t633 * t726 + t635 * t722;
t489 = pkin(5) * t773 - qJ(6) * t561 + t524;
t729 = cos(qJ(1));
t869 = t727 * t729;
t725 = sin(qJ(1));
t728 = cos(qJ(2));
t871 = t725 * t728;
t613 = t723 * t871 + t869;
t867 = t729 * t723;
t870 = t727 * t728;
t614 = t725 * t870 - t867;
t547 = t613 * t722 + t614 * t726;
t872 = t725 * t727;
t615 = t728 * t867 - t872;
t868 = t728 * t729;
t616 = t723 * t725 + t727 * t868;
t553 = t615 * t722 + t616 * t726;
t823 = qJDD(1) * t724;
t668 = t727 * t823;
t825 = qJD(1) * qJD(2);
t802 = t728 * t825;
t836 = qJD(3) * t724;
t807 = t723 * t836;
t824 = qJD(2) * qJD(3);
t554 = qJD(1) * t807 - t723 * qJDD(2) - t668 + (-t802 - t824) * t727;
t703 = t728 * qJDD(1);
t907 = -t724 * t825 + t703;
t631 = qJDD(3) - t907;
t789 = pkin(2) * t724 - pkin(8) * t728;
t646 = t789 * qJD(2);
t709 = t724 * pkin(8);
t714 = t728 * pkin(2);
t816 = -pkin(1) - t714;
t769 = t816 - t709;
t574 = qJD(1) * t646 + qJDD(1) * t769;
t610 = pkin(7) * t907 + qJDD(2) * pkin(8);
t621 = t769 * qJD(1);
t843 = qJD(1) * t728;
t700 = pkin(7) * t843;
t656 = qJD(2) * pkin(8) + t700;
t835 = qJD(3) * t727;
t837 = qJD(3) * t723;
t793 = -t574 * t727 + t610 * t723 + t621 * t837 + t656 * t835;
t770 = qJDD(4) + t793;
t897 = pkin(3) + pkin(4);
t480 = pkin(9) * t554 - t631 * t897 + t770;
t765 = t823 + t824;
t745 = t765 + t802;
t702 = t727 * qJDD(2);
t801 = qJD(1) * t835;
t762 = t724 * t801 - t702;
t819 = t574 * t723 + t610 * t727 + t621 * t835;
t619 = t631 * qJ(4);
t678 = -qJD(3) + t843;
t659 = t678 * qJD(4);
t909 = t619 - t659;
t483 = t762 * pkin(9) + (pkin(9) * t745 - qJD(3) * t656) * t723 + t819 + t909;
t566 = t621 * t727 - t656 * t723;
t827 = qJD(4) - t566;
t929 = -pkin(9) * t635 + t827;
t512 = t678 * t897 + t929;
t567 = t621 * t723 + t656 * t727;
t528 = pkin(9) * t633 + t567;
t661 = t678 * qJ(4);
t522 = t528 - t661;
t832 = qJD(5) * t726;
t833 = qJD(5) * t722;
t795 = -t480 * t722 - t483 * t726 - t512 * t832 + t522 * t833;
t637 = t722 * t723 + t726 * t727;
t917 = t637 * t724;
t741 = g(1) * t553 + g(2) * t547 + g(3) * t917 + t795;
t944 = -t489 * t773 - t741;
t943 = t524 * t773 + t741;
t942 = t727 * t843 - t835;
t814 = t723 * t843;
t941 = t814 - t837;
t796 = qJD(3) + t843;
t768 = t796 * qJD(2);
t735 = (-t768 - t823) * t723 - t762;
t490 = t554 * t726 - t633 * t832 + t635 * t833 + t722 * t735;
t666 = qJD(5) + t678;
t932 = -t666 * t773 + t490;
t940 = t932 * MDP(30);
t834 = qJD(5) * t561;
t491 = -t554 * t722 + t726 * t735 + t834;
t620 = -qJDD(5) + t631;
t898 = t561 ^ 2;
t933 = t561 * t773;
t939 = (-t773 ^ 2 + t898) * MDP(23) + (t561 * t666 - t491) * MDP(25) + MDP(22) * t933 - t932 * MDP(24) - t620 * MDP(26);
t896 = pkin(8) - pkin(9);
t658 = t896 * t727;
t815 = -pkin(7) * t723 - pkin(3);
t748 = -pkin(9) * t870 + (-pkin(4) + t815) * t724;
t643 = t789 * qJD(1);
t878 = t643 * t727;
t938 = qJD(1) * t748 - qJD(3) * t658 - t878;
t617 = t723 * t643;
t855 = qJ(4) * t844 + t617;
t874 = t724 * t727;
t875 = t723 * t728;
t937 = (-pkin(7) * t874 + pkin(9) * t875) * qJD(1) + t855 + t896 * t837;
t934 = t524 * t561;
t569 = qJD(5) * t637 - t722 * t837 - t726 * t835;
t603 = t637 * t728;
t861 = qJD(1) * t603 + t569;
t860 = -t722 * t942 + t723 * t832 + t726 * t941 - t727 * t833;
t931 = qJ(4) * t942 - t723 * qJD(4) - t700;
t508 = pkin(5) * t561 + qJ(6) * t773;
t876 = t723 * t726;
t600 = t722 * t874 - t724 * t876;
t774 = -t615 * t726 + t616 * t722;
t794 = -t480 * t726 + t483 * t722 + t512 * t833 + t522 * t832;
t910 = -t613 * t726 + t614 * t722;
t739 = g(1) * t774 + g(2) * t910 + g(3) * t600 - t794;
t928 = -pkin(5) * t774 + qJ(6) * t553;
t926 = t489 * t561 + qJDD(6);
t925 = -pkin(5) * t910 + qJ(6) * t547;
t924 = t666 ^ 2;
t488 = t512 * t722 + t522 * t726;
t486 = qJ(6) * t666 + t488;
t920 = t486 * t666;
t771 = t722 * t727 - t876;
t756 = t724 * t771;
t784 = g(1) * t729 + g(2) * t725;
t915 = t724 * t784;
t914 = qJ(4) * t833 + t528 * t722 - t726 * t929 + t832 * t897;
t657 = t896 * t723;
t772 = t657 * t726 - t658 * t722;
t913 = -qJD(5) * t772 + t722 * t938 + t726 * t937;
t577 = t657 * t722 + t658 * t726;
t912 = -qJD(5) * t577 + t722 * t937 - t726 * t938;
t847 = -t709 - t714;
t648 = -pkin(1) + t847;
t680 = pkin(7) * t875;
t713 = t728 * pkin(3);
t555 = pkin(4) * t728 + t680 + t713 + (-pkin(9) * t724 - t648) * t727;
t682 = pkin(7) * t870;
t850 = t648 * t723 + t682;
t578 = -qJ(4) * t728 + t850;
t877 = t723 * t724;
t565 = pkin(9) * t877 + t578;
t911 = t555 * t722 + t565 * t726;
t820 = t897 * t723;
t859 = -qJD(3) * t820 + t814 * t897 - t931;
t848 = qJ(4) * t726 - t722 * t897;
t822 = -MDP(27) - MDP(29);
t821 = MDP(28) - MDP(31);
t707 = t723 * qJ(4);
t905 = pkin(3) * t727 + pkin(2) + t707;
t891 = pkin(8) * t631;
t903 = t542 * t678 + t891;
t885 = g(3) * t728;
t890 = pkin(8) * t678;
t902 = qJD(3) * t890 - t885 + t915;
t782 = -qJD(3) * t682 + t646 * t727 - t648 * t837;
t506 = pkin(9) * t807 + qJD(2) * t748 - t782;
t840 = qJD(2) * t724;
t854 = t646 * t723 + t648 * t835;
t818 = qJ(4) * t840 + t854;
t507 = (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t874 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t723) * t728 + t818;
t900 = -qJD(5) * t911 + t506 * t726 - t507 * t722;
t899 = -0.2e1 * pkin(1);
t732 = qJD(1) ^ 2;
t895 = pkin(1) * t732;
t894 = pkin(3) * t631;
t893 = pkin(5) * t620;
t892 = pkin(7) * t633;
t889 = g(1) * t725;
t886 = g(2) * t729;
t884 = qJ(6) * t620;
t883 = t488 * t666;
t535 = -t661 + t567;
t882 = t535 * t678;
t880 = t633 * t678;
t879 = t635 * t678;
t873 = t724 * t729;
t866 = pkin(5) * t860 + qJ(6) * t861 + qJD(6) * t771 + t859;
t865 = qJ(6) * t844 - t913;
t864 = -pkin(5) * t844 - t912;
t858 = -qJD(6) - t914;
t857 = qJD(5) * t848 + t528 * t726 + t722 * t929;
t856 = -pkin(3) * t941 + t931;
t573 = pkin(3) * t635 + qJ(4) * t633;
t853 = (g(1) * t869 + g(2) * t872) * t724;
t809 = t728 * t828;
t851 = qJ(4) * t809 + qJD(4) * t874;
t718 = t724 ^ 2;
t846 = -t728 ^ 2 + t718;
t842 = qJD(2) * t635;
t839 = qJD(2) * t728;
t838 = qJD(3) * t633;
t830 = t633 * qJD(2);
t829 = t655 * qJD(3);
t487 = t512 * t726 - t522 * t722;
t826 = qJD(6) - t487;
t697 = pkin(7) * t823;
t611 = -qJDD(2) * pkin(2) + pkin(7) * t802 + t697;
t817 = g(1) * t868 + g(2) * t871 + g(3) * t724;
t811 = t678 * t828;
t810 = t723 * t839;
t808 = t678 * t837;
t806 = t724 * t835;
t800 = -pkin(3) * t613 + qJ(4) * t614;
t799 = -pkin(3) * t615 + qJ(4) * t616;
t797 = t648 * t727 - t680;
t627 = pkin(4) * t727 + t905;
t792 = pkin(3) * t870 + qJ(4) * t875 - t847;
t791 = -pkin(7) - t820;
t529 = -pkin(4) * t635 - t573;
t790 = t815 * t724;
t788 = -g(1) * t910 + g(2) * t774;
t787 = g(1) * t547 - g(2) * t553;
t786 = -g(1) * t613 + g(2) * t615;
t785 = g(1) * t614 - g(2) * t616;
t783 = -pkin(3) * t614 + pkin(7) * t729 - qJ(4) * t613;
t781 = -pkin(5) * t600 + qJ(6) * t917;
t780 = t829 - t891;
t779 = -qJ(4) * t722 - t726 * t897;
t776 = t555 * t726 - t565 * t722;
t766 = -g(1) * t615 - g(2) * t613 - g(3) * t877;
t494 = -pkin(3) * t735 + qJ(4) * t554 - qJD(4) * t635 + t611;
t761 = -t727 * t897 - t707;
t760 = -pkin(7) * qJDD(2) + t825 * t899;
t759 = t631 * t727 + t808;
t754 = -t656 * t837 + t819;
t753 = t506 * t722 + t507 * t726 + t555 * t832 - t565 * t833;
t731 = qJD(2) ^ 2;
t750 = pkin(7) * t731 + qJDD(1) * t899 - t889;
t674 = qJ(4) * t874;
t575 = t724 * t791 + t674;
t749 = pkin(1) * t729 + pkin(2) * t868 + pkin(3) * t616 + pkin(7) * t725 + pkin(8) * t873 + qJ(4) * t615;
t582 = t725 * t917;
t584 = t729 * t917;
t747 = g(1) * t584 + g(2) * t582 - g(3) * t603 - t620 * t772;
t581 = t725 * t756;
t583 = t729 * t756;
t602 = t722 * t870 - t726 * t875;
t746 = g(1) * t583 + g(2) * t581 - g(3) * t602 - t577 * t620;
t744 = -t766 - t793;
t738 = -t567 * t678 + t744;
t737 = -t739 + t926;
t736 = -t666 * t857 - t739;
t733 = g(1) * t616 + g(2) * t614 + g(3) * t874 - t566 * t678 - t754;
t513 = t761 * t836 + t791 * t839 + t851;
t484 = pkin(4) * t735 - t494;
t689 = g(2) * t873;
t685 = pkin(8) * t868;
t681 = pkin(8) * t871;
t642 = pkin(5) - t779;
t641 = -qJ(6) + t848;
t593 = -t674 + (pkin(3) * t723 + pkin(7)) * t724;
t579 = t713 - t797;
t572 = qJD(1) * t790 - t878;
t571 = -pkin(7) * t813 + t855;
t533 = pkin(3) * t678 + t827;
t531 = pkin(5) * t637 + qJ(6) * t771 + t627;
t525 = qJ(4) * t807 + pkin(7) * t839 + (t806 + t810) * pkin(3) - t851;
t523 = qJD(2) * t790 - t782;
t521 = -t554 - t880;
t520 = -qJD(4) * t728 + (-t724 * t828 - t728 * t837) * pkin(7) + t818;
t515 = qJD(2) * t603 + (qJD(3) - qJD(5)) * t756;
t514 = t569 * t724 + t722 * t809 - t726 * t810;
t511 = t575 - t781;
t502 = -pkin(5) * t728 - t776;
t501 = qJ(6) * t728 + t911;
t495 = t770 - t894;
t493 = t754 + t909;
t492 = -t508 + t529;
t485 = -pkin(5) * t666 + t826;
t476 = pkin(5) * t514 - qJ(6) * t515 - qJD(6) * t917 + t513;
t475 = pkin(5) * t840 - t900;
t474 = -qJ(6) * t840 + qJD(6) * t728 + t753;
t473 = pkin(5) * t491 + qJ(6) * t490 - qJD(6) * t561 + t484;
t472 = qJDD(6) + t794 + t893;
t471 = qJD(6) * t666 - t795 - t884;
t1 = [(t724 * t750 + t728 * t760 + t689) * MDP(10) + ((-t728 * t830 + (t702 + (-t635 - t813) * qJD(3)) * t724) * t727 + (-t635 * t839 + (-t727 * t768 + t554 - t668 + t838) * t724) * t723) * MDP(12) + (-t472 * t728 + t473 * t600 - t475 * t666 + t476 * t773 + t485 * t840 + t489 * t514 + t491 * t511 + t502 * t620 + t787) * MDP(29) + (-t491 * t728 - t514 * t666 + t600 * t620 + t773 * t840) * MDP(25) + t784 * MDP(3) + (qJDD(1) * t718 + 0.2e1 * t724 * t802) * MDP(4) + (-t554 * t874 + (-t807 + t809) * t635) * MDP(11) + (t760 * t724 + (-t750 - t886) * t728) * MDP(9) + (-g(1) * t783 - g(2) * t749 + t493 * t578 + t494 * t593 + t495 * t579 + t535 * t520 + t533 * t523 + t542 * t525 - t769 * t889) * MDP(21) + (-t520 * t633 + t523 * t635 - t579 * t554 + t578 * t702 - t689 + (t533 * t870 + (-t535 * t728 - t578 * t796) * t723) * qJD(2) + (t889 + t495 * t727 + (-qJDD(1) * t578 - t493) * t723 + (-t533 * t723 + (-qJD(1) * t578 - t535) * t727) * qJD(3)) * t724) * MDP(19) + (-t886 + t889) * MDP(2) + (-t782 * t678 + t797 * t631 + ((t655 * t723 + t892) * qJD(2) + t793) * t728 + (t727 * t829 + t566 * qJD(2) + t611 * t723 + (-t702 + (qJDD(1) * t723 + t801) * t724 + (-t678 + t796) * t841) * pkin(7)) * t724 + t785) * MDP(16) + (qJDD(2) * t724 + t728 * t731) * MDP(6) + (qJDD(2) * t728 - t724 * t731) * MDP(7) + (t484 * t600 - t487 * t840 + t575 * t491 + t513 * t773 + t524 * t514 - t776 * t620 + t666 * t900 - t794 * t728 + t787) * MDP(27) + ((t554 - t811) * t728 + (t759 + t842) * t724) * MDP(13) + (-t702 * t728 + (-t830 + (t678 + t843) * t835) * t724 + ((-t631 + t703) * t724 + (t678 + t796) * t839) * t723) * MDP(14) + qJDD(1) * MDP(1) + 0.2e1 * (t703 * t724 - t825 * t846) * MDP(5) + (-t631 * t728 - t678 * t840) * MDP(15) + (-t620 * t728 - t666 * t840) * MDP(26) + (t495 * t728 + t523 * t678 + t525 * t633 - t579 * t631 - t593 * t702 + (-t533 * qJD(2) + (qJD(1) * t593 + t542) * t835) * t724 + ((qJDD(1) * t593 + t494) * t724 + (t542 * t728 + t593 * t796) * qJD(2)) * t723 + t785) * MDP(18) + (-t520 * t678 - t525 * t635 + t554 * t593 + t578 * t631 + (-t542 * t828 - t493) * t728 + (qJD(2) * t535 - t494 * t727 + t542 * t837) * t724 - t786) * MDP(20) + (t471 * t501 + t486 * t474 + t473 * t511 + t489 * t476 + t472 * t502 + t485 * t475 - g(1) * (-pkin(4) * t614 - pkin(5) * t547 - qJ(6) * t910 + t783) - g(2) * (pkin(4) * t616 + pkin(5) * t553 - pkin(9) * t873 + qJ(6) * t774 + t749) - (-t724 * t896 + t816) * t889) * MDP(32) + (-t490 * t917 + t515 * t561) * MDP(22) + (t490 * t600 - t491 * t917 - t514 * t561 - t515 * t773) * MDP(23) + (-t471 * t600 + t472 * t917 - t474 * t773 + t475 * t561 + t485 * t515 - t486 * t514 - t490 * t502 - t491 * t501 - t724 * t889 + t689) * MDP(30) + (-t490 * t728 + t515 * t666 - t561 * t840 - t620 * t917) * MDP(24) + (t471 * t728 - t473 * t917 + t474 * t666 - t476 * t561 - t486 * t840 - t489 * t515 + t490 * t511 - t501 * t620 - t788) * MDP(31) + (t484 * t917 + t488 * t840 - t575 * t490 + t513 * t561 + t524 * t515 + t620 * t911 - t666 * t753 + t728 * t795 + t788) * MDP(28) + (t854 * t678 - t850 * t631 + (t655 * t828 + (-t808 + t842) * pkin(7) + t754) * t728 + (-t723 * t829 - t567 * qJD(2) + t611 * t727 + (-t554 - t811) * pkin(7)) * t724 + t786) * MDP(17); (-t678 * t835 + t723 * t631 + (-t635 * t724 + t678 * t870) * qJD(1)) * MDP(13) + (t473 * t637 - t485 * t844 + t489 * t860 + t491 * t531 - t666 * t864 + t773 * t866 + t747) * MDP(29) + (t620 * t637 - t666 * t860 - t773 * t844) * MDP(25) + ((t633 * t724 - t678 * t875) * qJD(1) + t759) * MDP(14) + (-t554 * t723 - t727 * t879) * MDP(11) + ((-t554 + t880) * t727 + (-t635 * qJD(3) + t702 - t765 * t723 + (-t806 + (t635 - t841) * t728) * qJD(1)) * t723) * MDP(12) + (pkin(2) * t702 + (-t566 * t724 - t728 * t892) * qJD(1) + (-t885 + t643 * t678 - t611 + (-pkin(2) * t844 + t890) * qJD(3)) * t727 + (-pkin(2) * t765 + (pkin(7) * t678 * t724 + (-t655 - t721) * t728) * qJD(1) + t780) * t723 + t853) * MDP(16) + (-t885 - t697 + (t784 + t895) * t724) * MDP(9) + ((-pkin(7) * qJDD(1) + t895) * t728 + t817) * MDP(10) + (-t724 * t728 * MDP(4) + MDP(5) * t846) * t732 + (t571 * t633 - t572 * t635 + (t493 - t678 * t533 + (t702 + (t635 - t813) * qJD(3)) * pkin(8)) * t727 + (t495 + t882 + (-t727 * t745 - t554 + t838) * pkin(8)) * t723 - t817) * MDP(19) + t678 * MDP(15) * t844 + t666 * MDP(26) * t844 + (t484 * t637 + t487 * t844 + t627 * t491 + t860 * t524 + t666 * t912 + t773 * t859 + t747) * MDP(27) + qJDD(2) * MDP(8) + (-t535 * t571 - t533 * t572 - g(1) * t685 - g(2) * t681 - g(3) * t792 + t856 * t542 + (t493 * t727 + t495 * t723 + (t533 * t727 - t535 * t723) * qJD(3)) * pkin(8) + (-t494 + t915) * t905) * MDP(21) + (t533 * t844 - t572 * t678 + t905 * t702 + t856 * t633 + (-t885 - t494 + (-t844 * t905 + t890) * qJD(3)) * t727 + (-t745 * t905 - t903) * t723 + t853) * MDP(18) + (-t535 * t844 - t554 * t905 + t571 * t678 - t856 * t635 + t903 * t727 + (-t494 + t902) * t723) * MDP(20) + (pkin(2) * t554 - t617 * t678 + t780 * t727 + (-t655 * t870 + t567 * t724 + (-t635 * t728 + t678 * t874) * pkin(7)) * qJD(1) + (t611 - t902) * t723) * MDP(17) + (t490 * t637 + t491 * t771 - t561 * t860 + t773 * t861) * MDP(23) + (t473 * t771 + t486 * t844 + t489 * t861 + t490 * t531 - t561 * t866 + t666 * t865 + t746) * MDP(31) + (t561 * t844 + t620 * t771 - t666 * t861) * MDP(24) + (t490 * t771 - t561 * t861) * MDP(22) + (-t484 * t771 - t488 * t844 - t627 * t490 - t861 * t524 + t859 * t561 + t666 * t913 - t746) * MDP(28) + (-t471 * t637 - t472 * t771 - t485 * t861 - t486 * t860 + t490 * t772 - t491 * t577 + t561 * t864 - t773 * t865 + t817) * MDP(30) + (t471 * t577 + t473 * t531 - t472 * t772 - g(1) * (-pkin(5) * t584 - pkin(9) * t868 - qJ(6) * t583 + t685) - g(2) * (-pkin(5) * t582 - pkin(9) * t871 - qJ(6) * t581 + t681) - g(3) * (pkin(4) * t870 + pkin(5) * t603 + qJ(6) * t602 + t792) + t866 * t489 + t865 * t486 + t864 * t485 + (g(3) * pkin(9) + t784 * (pkin(2) - t761)) * t724) * MDP(32) + MDP(7) * t703 + MDP(6) * t823; t521 * MDP(13) + (t471 * t641 + t472 * t642 - t489 * t492 - g(1) * (-pkin(4) * t615 + t799 - t928) - g(2) * (-pkin(4) * t613 + t800 - t925) - g(3) * (-t724 * t820 + t674 - t781) + t858 * t486 + t857 * t485) * MDP(32) + (-t529 * t561 + t848 * t620 + t666 * t914 - t943) * MDP(28) + (t492 * t561 + (qJ(6) - t641) * t620 + (-qJD(6) + t858) * t666 - t944) * MDP(31) + (-t529 * t773 - t620 * t779 + t736 + t934) * MDP(27) + (-t492 * t773 + (pkin(5) + t642) * t620 + t736 + t926) * MDP(29) + (-t490 * t642 - t491 * t641 + (-t485 - t858) * t773 + (-t486 + t857) * t561) * MDP(30) + (-t635 * t655 + t738) * MDP(16) + (pkin(3) * t554 + (t535 - t567) * t635 + (t533 - t827) * t633 + t735 * qJ(4)) * MDP(19) + (t493 * qJ(4) - t495 * pkin(3) - t542 * t573 - t533 * t567 - g(1) * t799 - g(2) * t800 - g(3) * (-pkin(3) * t877 + t674) + t827 * t535) * MDP(21) + (t735 - t879) * MDP(14) + (-t542 * t633 + t573 * t635 + 0.2e1 * t619 - 0.2e1 * t659 - t733) * MDP(20) + (t633 * t655 + t733) * MDP(17) + (-t542 * t635 - t573 * t633 - qJDD(4) + t738 + 0.2e1 * t894) * MDP(18) + t631 * MDP(15) + (-t633 ^ 2 + t635 ^ 2) * MDP(12) + t635 * t633 * MDP(11) - t939; -t631 * MDP(18) + t521 * MDP(19) - t678 ^ 2 * MDP(20) + (qJDD(4) - t744 + t882 - t894) * MDP(21) + t766 * MDP(32) + (MDP(18) * t633 - MDP(20) * t635 + MDP(21) * t542 - MDP(32) * t489 - t561 * t821 + t773 * t822) * t635 + (t940 + (-t472 + t920) * MDP(32) + t822 * t620 - t821 * t924) * t726 + ((t561 * t678 - t491 + t834) * MDP(30) + (t485 * t666 + t471) * MDP(32) + t821 * t620 + t822 * t924) * t722; (t739 + t883 - t934) * MDP(27) + (t487 * t666 + t943) * MDP(28) + (-t508 * t773 - t737 + t883 - 0.2e1 * t893) * MDP(29) + (pkin(5) * t490 - qJ(6) * t491 + (t486 - t488) * t561 + (t485 - t826) * t773) * MDP(30) + (-0.2e1 * t884 + t508 * t561 + (0.2e1 * qJD(6) - t487) * t666 + t944) * MDP(31) + (-t472 * pkin(5) - g(1) * t928 - g(2) * t925 - g(3) * t781 + t471 * qJ(6) - t485 * t488 + t826 * t486 - t489 * t508) * MDP(32) + t939; (t620 + t933) * MDP(29) - t940 + (-t898 - t924) * MDP(31) + (t737 + t893 - t920) * MDP(32);];
tau  = t1;
