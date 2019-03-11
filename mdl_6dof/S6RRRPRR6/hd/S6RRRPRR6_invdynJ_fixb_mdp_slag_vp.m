% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:29:58
% EndTime: 2019-03-09 18:30:26
% DurationCPUTime: 19.11s
% Computational Cost: add. (14031->683), mult. (32181->896), div. (0->0), fcn. (24690->16), ass. (0->289)
t792 = sin(qJ(3));
t793 = sin(qJ(2));
t883 = qJD(1) * t793;
t860 = t792 * t883;
t797 = cos(qJ(3));
t868 = t797 * qJD(2);
t737 = t860 - t868;
t879 = qJD(2) * t792;
t739 = t797 * t883 + t879;
t787 = sin(pkin(11));
t788 = cos(pkin(11));
t664 = t737 * t787 - t739 * t788;
t791 = sin(qJ(5));
t796 = cos(qJ(5));
t826 = -t737 * t788 - t739 * t787;
t608 = t664 * t791 + t796 * t826;
t790 = sin(qJ(6));
t795 = cos(qJ(6));
t948 = -t796 * t664 + t791 * t826;
t537 = t608 * t790 + t795 * t948;
t921 = t948 * t790;
t541 = t608 * t795 - t921;
t798 = cos(qJ(2));
t781 = t798 * qJDD(1);
t867 = qJD(1) * qJD(2);
t938 = -t793 * t867 + t781;
t729 = qJDD(3) - t938;
t726 = qJDD(5) + t729;
t715 = qJDD(6) + t726;
t974 = t715 * MDP(31) + (t537 ^ 2 - t541 ^ 2) * MDP(28) - t541 * MDP(27) * t537;
t882 = qJD(1) * t798;
t765 = -qJD(3) + t882;
t853 = t798 * t867;
t866 = qJDD(1) * t793;
t875 = qJD(3) * t793;
t951 = -qJD(1) * t875 + qJDD(2);
t657 = qJD(3) * t868 + (t853 + t866) * t797 + t951 * t792;
t658 = t792 * (qJD(2) * (qJD(3) + t882) + t866) - t951 * t797;
t587 = -t657 * t787 - t658 * t788;
t588 = t657 * t788 - t658 * t787;
t871 = qJD(5) * t796;
t872 = qJD(5) * t791;
t519 = t791 * t587 + t796 * t588 + t664 * t872 + t826 * t871;
t520 = qJD(5) * t948 - t796 * t587 + t588 * t791;
t869 = qJD(6) * t795;
t863 = t795 * t519 - t790 * t520 + t608 * t869;
t870 = qJD(6) * t790;
t492 = -t870 * t948 + t863;
t848 = t519 * t790 + t795 * t520;
t803 = -qJD(6) * t537 - t848;
t755 = -qJD(5) + t765;
t922 = t608 * t755;
t923 = t948 * t755;
t748 = -qJD(6) + t755;
t971 = t541 * t748;
t972 = t537 * t748;
t973 = t726 * MDP(24) + (-t520 - t923) * MDP(23) + (-t608 ^ 2 + t948 ^ 2) * MDP(21) - t608 * MDP(20) * t948 + (t519 + t922) * MDP(22) + (t492 + t971) * MDP(29) + (t803 - t972) * MDP(30) + t974;
t909 = t797 * t798;
t933 = pkin(3) * t793;
t822 = -qJ(4) * t909 + t933;
t789 = -qJ(4) - pkin(8);
t851 = qJD(3) * t789;
t837 = pkin(2) * t793 - pkin(8) * t798;
t740 = t837 * qJD(1);
t890 = pkin(7) * t860 + t797 * t740;
t970 = -qJD(1) * t822 - qJD(4) * t792 + t797 * t851 - t890;
t722 = t792 * t740;
t873 = qJD(4) * t797;
t913 = t793 * t797;
t914 = t792 * t798;
t969 = t722 + (-pkin(7) * t913 - qJ(4) * t914) * qJD(1) - t792 * t851 - t873;
t752 = -qJD(2) * pkin(2) + pkin(7) * t883;
t680 = pkin(3) * t737 + qJD(4) + t752;
t617 = -pkin(4) * t826 + t680;
t557 = -pkin(5) * t608 + t617;
t782 = qJ(3) + pkin(11) + qJ(5);
t772 = qJ(6) + t782;
t762 = sin(t772);
t763 = cos(t772);
t799 = cos(qJ(1));
t794 = sin(qJ(1));
t912 = t794 * t798;
t681 = t762 * t912 + t763 * t799;
t908 = t798 * t799;
t683 = -t762 * t908 + t763 * t794;
t747 = -pkin(2) * t798 - pkin(8) * t793 - pkin(1);
t727 = t747 * qJD(1);
t776 = pkin(7) * t882;
t753 = qJD(2) * pkin(8) + t776;
t674 = t797 * t727 - t753 * t792;
t635 = -qJ(4) * t739 + t674;
t624 = -pkin(3) * t765 + t635;
t675 = t727 * t792 + t753 * t797;
t636 = -qJ(4) * t737 + t675;
t629 = t787 * t636;
t565 = t788 * t624 - t629;
t956 = pkin(9) * t664;
t545 = -pkin(4) * t765 + t565 + t956;
t916 = t788 * t636;
t566 = t787 * t624 + t916;
t946 = pkin(9) * t826;
t548 = t566 + t946;
t511 = t545 * t791 + t548 * t796;
t741 = t837 * qJD(2);
t676 = qJD(1) * t741 + qJDD(1) * t747;
t670 = t797 * t676;
t712 = pkin(7) * t938 + qJDD(2) * pkin(8);
t544 = pkin(3) * t729 - qJ(4) * t657 - qJD(3) * t675 - qJD(4) * t739 - t712 * t792 + t670;
t874 = qJD(3) * t797;
t876 = qJD(3) * t792;
t814 = t792 * t676 + t797 * t712 + t727 * t874 - t753 * t876;
t550 = -qJ(4) * t658 - qJD(4) * t737 + t814;
t508 = t788 * t544 - t550 * t787;
t496 = pkin(4) * t729 - pkin(9) * t588 + t508;
t509 = t787 * t544 + t788 * t550;
t500 = pkin(9) * t587 + t509;
t849 = t796 * t496 - t791 * t500;
t804 = -qJD(5) * t511 + t849;
t486 = pkin(5) * t726 - pkin(10) * t519 + t804;
t839 = -t791 * t496 - t796 * t500 - t545 * t871 + t548 * t872;
t487 = -pkin(10) * t520 - t839;
t850 = t795 * t486 - t790 * t487;
t927 = g(3) * t793;
t949 = -g(1) * t683 + g(2) * t681 - t557 * t537 + t762 * t927 + t850;
t967 = pkin(10) * t608;
t505 = t511 + t967;
t501 = t505 * t870;
t682 = t762 * t799 - t763 * t912;
t684 = t762 * t794 + t763 * t908;
t950 = g(1) * t684 - g(2) * t682 - t557 * t541 + t763 * t927 + t501;
t731 = t787 * t797 + t788 * t792;
t815 = t731 * t798;
t959 = qJD(1) * t815 - t731 * qJD(3);
t823 = t787 * t792 - t788 * t797;
t966 = t765 * t823;
t896 = t787 * t969 + t788 * t970;
t895 = t787 * t970 - t788 * t969;
t768 = sin(t782);
t769 = cos(t782);
t688 = t768 * t799 - t769 * t912;
t690 = t768 * t794 + t769 * t908;
t965 = g(1) * t690 - g(2) * t688 - t617 * t608 + t769 * t927 + t839;
t962 = pkin(10) * t948;
t960 = pkin(4) * t883 + pkin(9) * t966 - t896;
t953 = pkin(9) * t959 + t895;
t774 = pkin(7) * t866;
t713 = -qJDD(2) * pkin(2) + pkin(7) * t853 + t774;
t836 = g(1) * t799 + g(2) * t794;
t926 = g(3) * t798;
t810 = t793 * t836 - t926;
t958 = qJD(3) * pkin(8) * t765 - t713 + t810;
t687 = t768 * t912 + t769 * t799;
t689 = -t768 * t908 + t769 * t794;
t957 = -g(1) * t689 + g(2) * t687 - t617 * t948 + t768 * t927 + t804;
t827 = -t731 * t791 - t796 * t823;
t899 = qJD(5) * t827 + t791 * t959 + t796 * t966;
t663 = t731 * t796 - t791 * t823;
t898 = qJD(5) * t663 + t791 * t966 - t796 * t959;
t934 = pkin(3) * t792;
t952 = pkin(3) * t876 - t882 * t934 - t776;
t943 = t960 * t796;
t733 = t797 * t747;
t932 = pkin(7) * t792;
t671 = -qJ(4) * t913 + t733 + (-pkin(3) - t932) * t798;
t767 = pkin(7) * t909;
t889 = t792 * t747 + t767;
t915 = t792 * t793;
t677 = -qJ(4) * t915 + t889;
t611 = t788 * t671 - t677 * t787;
t704 = t823 * t793;
t580 = -pkin(4) * t798 + pkin(9) * t704 + t611;
t612 = t787 * t671 + t788 * t677;
t703 = t731 * t793;
t583 = -pkin(9) * t703 + t612;
t900 = t791 * t580 + t796 * t583;
t749 = t789 * t792;
t750 = t789 * t797;
t678 = t788 * t749 + t750 * t787;
t649 = -pkin(9) * t731 + t678;
t679 = t787 * t749 - t788 * t750;
t650 = -pkin(9) * t823 + t679;
t897 = t791 * t649 + t796 * t650;
t894 = -pkin(4) * t959 + t952;
t573 = -t635 * t787 - t916;
t558 = t573 - t946;
t574 = t788 * t635 - t629;
t559 = t574 + t956;
t770 = pkin(3) * t788 + pkin(4);
t935 = pkin(3) * t787;
t838 = t796 * t770 - t791 * t935;
t942 = t838 * qJD(5) - t791 * t558 - t796 * t559;
t708 = t770 * t791 + t796 * t935;
t941 = t708 * qJD(5) + t796 * t558 - t559 * t791;
t939 = t649 * t871 - t650 * t872 - t791 * t960 + t796 * t953;
t718 = t792 * t912 + t797 * t799;
t720 = -t792 * t908 + t794 * t797;
t937 = -g(1) * t720 + g(2) * t718;
t920 = t657 * t792;
t919 = t737 * t765;
t918 = t739 * t765;
t917 = t739 * t797;
t510 = t796 * t545 - t548 * t791;
t504 = t510 - t962;
t498 = -pkin(5) * t755 + t504;
t911 = t795 * t498;
t910 = t795 * t505;
t907 = t941 - t967;
t906 = t942 + t962;
t602 = t663 * t790 - t795 * t827;
t905 = -qJD(6) * t602 - t790 * t898 + t795 * t899;
t603 = t663 * t795 + t790 * t827;
t904 = qJD(6) * t603 + t790 * t899 + t795 * t898;
t901 = pkin(5) * t898 + t894;
t878 = qJD(2) * t793;
t891 = t797 * t741 + t878 * t932;
t599 = -t793 * t873 + t822 * qJD(2) + (-t767 + (qJ(4) * t793 - t747) * t792) * qJD(3) + t891;
t892 = t792 * t741 + t747 * t874;
t610 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t913 + (-qJD(4) * t793 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t798) * t792 + t892;
t535 = t787 * t599 + t788 * t610;
t887 = pkin(3) * t915 + t793 * pkin(7);
t785 = t793 ^ 2;
t886 = -t798 ^ 2 + t785;
t881 = qJD(2) * t737;
t880 = qJD(2) * t739;
t877 = qJD(2) * t798;
t858 = t792 * t877;
t862 = pkin(3) * t858 + pkin(7) * t877 + t874 * t933;
t773 = pkin(3) * t797 + pkin(2);
t861 = pkin(7) + t934;
t859 = t765 * t868;
t857 = t798 * t868;
t856 = t765 * t876;
t855 = t765 * t874;
t534 = t788 * t599 - t610 * t787;
t644 = t731 * t875 + t787 * t858 - t788 * t857;
t527 = pkin(4) * t878 + pkin(9) * t644 + t534;
t643 = -qJD(2) * t815 + t823 * t875;
t530 = pkin(9) * t643 + t535;
t847 = t796 * t527 - t530 * t791;
t845 = t796 * t580 - t583 * t791;
t843 = t796 * t649 - t650 * t791;
t842 = -qJD(3) * t727 - t712;
t841 = qJD(6) * t498 + t487;
t672 = pkin(4) * t703 + t887;
t626 = pkin(3) * t739 - pkin(4) * t664;
t835 = g(1) * t794 - g(2) * t799;
t834 = t753 * t874 - t670;
t560 = -pkin(10) * t663 + t843;
t833 = -pkin(10) * t898 + qJD(6) * t560 + t939;
t561 = pkin(10) * t827 + t897;
t832 = pkin(5) * t883 + pkin(10) * t899 + qJD(5) * t897 + qJD(6) * t561 + t791 * t953 + t943;
t831 = -pkin(8) * t729 + qJD(3) * t752;
t489 = t790 * t498 + t910;
t642 = -t703 * t791 - t704 * t796;
t828 = -t796 * t703 + t704 * t791;
t576 = t642 * t790 - t795 * t828;
t577 = t642 * t795 + t790 * t828;
t825 = t773 * t798 - t789 * t793;
t694 = pkin(4) * t823 - t773;
t613 = -pkin(4) * t643 + t862;
t820 = pkin(1) + t825;
t819 = -0.2e1 * pkin(1) * t867 - pkin(7) * qJDD(2);
t818 = t729 * t792 - t855;
t817 = t729 * t797 + t856;
t813 = t791 * t527 + t796 * t530 + t580 * t871 - t583 * t872;
t801 = qJD(1) ^ 2;
t811 = pkin(1) * t801 + t836;
t622 = pkin(3) * t658 + qJDD(4) + t713;
t800 = qJD(2) ^ 2;
t807 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t800 + t835;
t551 = -pkin(4) * t587 + t622;
t721 = t792 * t794 + t797 * t908;
t719 = t792 * t799 - t794 * t909;
t706 = pkin(5) + t838;
t623 = -pkin(5) * t827 + t694;
t595 = -pkin(5) * t828 + t672;
t562 = pkin(5) * t948 + t626;
t556 = qJD(5) * t642 - t796 * t643 - t644 * t791;
t555 = qJD(5) * t828 + t643 * t791 - t644 * t796;
t528 = pkin(5) * t556 + t613;
t524 = pkin(10) * t828 + t900;
t521 = -pkin(5) * t798 - pkin(10) * t642 + t845;
t503 = qJD(6) * t577 + t555 * t790 + t795 * t556;
t502 = -qJD(6) * t576 + t555 * t795 - t556 * t790;
t499 = pkin(5) * t520 + t551;
t491 = -pkin(10) * t556 + t813;
t490 = pkin(5) * t878 - pkin(10) * t555 - qJD(5) * t900 + t847;
t488 = -t505 * t790 + t911;
t1 = [(t520 * t798 + t556 * t755 + t608 * t878 + t726 * t828) * MDP(23) + (t519 * t828 - t520 * t642 + t555 * t608 - t556 * t948) * MDP(21) + (t503 * t748 + t541 * t878 - t576 * t715 - t798 * t803) * MDP(30) + (-(t490 * t795 - t491 * t790) * t748 + (t521 * t795 - t524 * t790) * t715 - t850 * t798 + t488 * t878 - t528 * t541 - t595 * t803 + t499 * t576 + t557 * t503 - g(1) * t682 - g(2) * t684 + (-(-t521 * t790 - t524 * t795) * t748 + t489 * t798) * qJD(6)) * MDP(32) + (-t492 * t576 + t502 * t541 - t503 * t537 + t577 * t803) * MDP(28) + (-t492 * t798 - t502 * t748 + t537 * t878 + t577 * t715) * MDP(29) + (-t489 * t878 - g(1) * t681 - g(2) * t683 + t595 * t492 + t499 * t577 - t501 * t798 + t557 * t502 + t528 * t537 + ((-qJD(6) * t524 + t490) * t748 - t521 * t715 + t486 * t798) * t790 + ((qJD(6) * t521 + t491) * t748 - t524 * t715 + t841 * t798) * t795) * MDP(33) + (t492 * t577 + t502 * t537) * MDP(27) + (t508 * t704 - t509 * t703 + t534 * t664 + t535 * t826 + t565 * t644 + t566 * t643 + t587 * t612 - t588 * t611 + t793 * t835) * MDP(18) + ((-t737 * t797 - t739 * t792) * t877 + (-t920 - t658 * t797 + (t737 * t792 - t917) * qJD(3)) * t793) * MDP(12) + (t657 * t913 + (-t792 * t875 + t857) * t739) * MDP(11) + 0.2e1 * (t781 * t793 - t867 * t886) * MDP(5) + (t509 * t612 + t566 * t535 + t508 * t611 + t565 * t534 + t622 * t887 + t680 * t862 + (-g(1) * t861 - g(2) * t820) * t799 + (g(1) * t820 - g(2) * t861) * t794) * MDP(19) + (-(-t747 * t876 + t891) * t765 + t733 * t729 - g(1) * t719 - g(2) * t721 + ((t855 + t881) * pkin(7) + (-pkin(7) * t729 + qJD(2) * t752 - t842) * t792 + t834) * t798 + (pkin(7) * t658 + qJD(2) * t674 + t713 * t792 + t752 * t874) * t793) * MDP(16) + (t892 * t765 - t889 * t729 - g(1) * t718 - g(2) * t720 + (t752 * t868 + (-t856 + t880) * pkin(7) + t814) * t798 + (-t752 * t876 - t675 * qJD(2) + t713 * t797 + (t657 - t859) * pkin(7)) * t793) * MDP(17) + (-t715 * t798 - t748 * t878) * MDP(31) + ((-t657 - t859) * t798 + (t817 + t880) * t793) * MDP(13) + ((t765 * t879 + t658) * t798 + (-t818 - t881) * t793) * MDP(14) + (-t729 * t798 - t765 * t878) * MDP(15) + (-t726 * t798 - t755 * t878) * MDP(24) + (qJDD(1) * t785 + 0.2e1 * t793 * t853) * MDP(4) + (qJDD(2) * t793 + t798 * t800) * MDP(6) + (qJDD(2) * t798 - t793 * t800) * MDP(7) + t835 * MDP(2) + t836 * MDP(3) + (t793 * t819 + t798 * t807) * MDP(9) + (-t793 * t807 + t798 * t819) * MDP(10) + (-g(1) * t687 - g(2) * t689 - t511 * t878 + t672 * t519 + t551 * t642 + t617 * t555 + t613 * t948 - t726 * t900 + t755 * t813 - t798 * t839) * MDP(26) + (-t519 * t798 - t555 * t755 + t642 * t726 + t878 * t948) * MDP(22) + (t519 * t642 + t555 * t948) * MDP(20) + qJDD(1) * MDP(1) + (-t847 * t755 + t845 * t726 - t849 * t798 + t510 * t878 - t613 * t608 + t672 * t520 - t551 * t828 + t617 * t556 - g(1) * t688 - g(2) * t690 + (t511 * t798 + t755 * t900) * qJD(5)) * MDP(25); (t843 * t726 + t694 * t520 - t551 * t827 + (t650 * t871 + (qJD(5) * t649 + t953) * t791 + t943) * t755 + t898 * t617 - t894 * t608 + t810 * t769) * MDP(25) + (t519 * t827 - t520 * t663 + t608 * t899 - t898 * t948) * MDP(21) + ((t560 * t795 - t561 * t790) * t715 - t623 * t803 + t499 * t602 + (t790 * t833 + t795 * t832) * t748 + t904 * t557 - t901 * t541 + t810 * t763) * MDP(32) + (t765 * MDP(15) - MDP(22) * t948 - MDP(23) * t608 + t755 * MDP(24) - t510 * MDP(25) + t511 * MDP(26) - MDP(29) * t537 - MDP(30) * t541 + t748 * MDP(31) - t488 * MDP(32) + t489 * MDP(33)) * t883 + (-t492 * t602 - t537 * t904 + t541 * t905 + t603 * t803) * MDP(28) + (t492 * t603 + t537 * t905) * MDP(27) + (-(t560 * t790 + t561 * t795) * t715 + t623 * t492 + t499 * t603 + (-t790 * t832 + t795 * t833) * t748 + t905 * t557 + t901 * t537 - t810 * t762) * MDP(33) + (t603 * t715 - t748 * t905) * MDP(29) + (t663 * t726 - t755 * t899) * MDP(22) + (-MDP(4) * t793 * t798 + MDP(5) * t886) * t801 + MDP(7) * t781 + MDP(6) * t866 + (-t602 * t715 + t748 * t904) * MDP(30) + (-t508 * t731 - t509 * t823 - t565 * t966 + t959 * t566 + t587 * t679 - t588 * t678 + t896 * t664 - t836 * t798 + t895 * t826 - t927) * MDP(18) + (t509 * t679 + t508 * t678 - t622 * t773 - g(3) * t825 + t952 * t680 + t895 * t566 + t896 * t565 + t836 * (t773 * t793 + t789 * t798)) * MDP(19) + (t726 * t827 + t755 * t898) * MDP(23) + (t793 * t811 - t774 - t926) * MDP(9) + (t927 + (-pkin(7) * qJDD(1) + t811) * t798) * MDP(10) + (-t765 * t917 + t920) * MDP(11) + ((t657 + t919) * t797 + (-t658 + t918) * t792) * MDP(12) + ((t737 * t793 - t765 * t914) * qJD(1) + t817) * MDP(14) + ((-t739 * t793 + t765 * t909) * qJD(1) + t818) * MDP(13) + (t519 * t663 + t899 * t948) * MDP(20) + (t694 * t519 + t551 * t663 + t899 * t617 - t897 * t726 + t755 * t939 - t768 * t810 + t894 * t948) * MDP(26) + qJDD(2) * MDP(8) + (-pkin(2) * t658 + t890 * t765 + t831 * t792 + (-t674 * t793 + (-pkin(7) * t737 - t752 * t792) * t798) * qJD(1) + t958 * t797) * MDP(16) + (-pkin(2) * t657 - t722 * t765 + t831 * t797 + (-t752 * t909 + t675 * t793 + (-t739 * t798 + t765 * t913) * pkin(7)) * qJD(1) - t958 * t792) * MDP(17); (t608 * t626 + t838 * t726 + t755 * t941 + t957) * MDP(25) + ((t587 * t787 - t588 * t788) * pkin(3) + (t565 - t574) * t826 + (-t566 - t573) * t664) * MDP(18) + (-t675 * t765 - t739 * t752 + (t842 + t927) * t792 - t834 + t937) * MDP(16) + (-t565 * t573 - t566 * t574 + (g(3) * t915 + t508 * t788 + t509 * t787 - t680 * t739 + t937) * pkin(3)) * MDP(19) + (-t562 * t537 + (-t706 * t715 - t486 + (-qJD(6) * t708 - t907) * t748) * t790 + (-t708 * t715 + (qJD(6) * t706 + t906) * t748 - t841) * t795 + t950) * MDP(33) + (-t658 - t918) * MDP(14) + t739 * t737 * MDP(11) + (t657 - t919) * MDP(13) + (g(1) * t721 - g(2) * t719 + g(3) * t913 - t674 * t765 + t737 * t752 - t814) * MDP(17) + ((t706 * t795 - t708 * t790) * t715 + t562 * t541 + (t906 * t790 + t907 * t795) * t748 + (-(-t706 * t790 - t708 * t795) * t748 - t489) * qJD(6) + t949) * MDP(32) + (-t737 ^ 2 + t739 ^ 2) * MDP(12) + t729 * MDP(15) + (-t626 * t948 - t708 * t726 + t755 * t942 + t965) * MDP(26) + t973; (-t664 ^ 2 - t826 ^ 2) * MDP(18) + (-t565 * t664 - t566 * t826 + t622 - t810) * MDP(19) + (t520 - t923) * MDP(25) + (t519 - t922) * MDP(26) + (-t803 - t972) * MDP(32) + (t492 - t971) * MDP(33); (-t511 * t755 + t957) * MDP(25) + (-t510 * t755 + t965) * MDP(26) + ((-t504 * t790 - t910) * t748 - t489 * qJD(6) + (t541 * t948 + t715 * t795 + t748 * t870) * pkin(5) + t949) * MDP(32) + ((t505 * t748 - t486) * t790 + (-t504 * t748 - t841) * t795 + (-t537 * t948 - t715 * t790 + t748 * t869) * pkin(5) + t950) * MDP(33) + t973; (t863 + t971) * MDP(29) + (-t848 - t972) * MDP(30) + (-t489 * t748 + t949) * MDP(32) + (-t790 * t486 - t795 * t487 - t488 * t748 + t950) * MDP(33) + (-MDP(29) * t921 - t537 * MDP(30) - MDP(32) * t489 - MDP(33) * t911) * qJD(6) + t974;];
tau  = t1;
