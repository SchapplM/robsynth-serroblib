% Calculate vector of inverse dynamics joint torques for
% S6RRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:36:16
% EndTime: 2019-03-10 03:36:40
% DurationCPUTime: 16.67s
% Computational Cost: add. (15748->682), mult. (36385->870), div. (0->0), fcn. (28775->18), ass. (0->322)
t1015 = qJD(5) + qJD(6);
t847 = cos(qJ(3));
t841 = sin(qJ(3));
t842 = sin(qJ(2));
t954 = qJD(1) * t842;
t927 = t841 * t954;
t848 = cos(qJ(2));
t953 = qJD(1) * t848;
t745 = t847 * t953 - t927;
t746 = -t841 * t953 - t847 * t954;
t840 = sin(qJ(4));
t846 = cos(qJ(4));
t704 = t846 * t745 + t746 * t840;
t1050 = -t1015 + t704;
t838 = sin(qJ(6));
t839 = sin(qJ(5));
t844 = cos(qJ(6));
t845 = cos(qJ(5));
t759 = t838 * t839 - t844 * t845;
t1048 = t1050 * t759;
t981 = t838 * t845;
t761 = t839 * t844 + t981;
t1049 = t1050 * t761;
t1033 = t704 * t845;
t947 = qJD(5) * t845;
t1047 = t1033 - t947;
t837 = qJ(2) + qJ(3);
t830 = qJ(4) + t837;
t814 = sin(t830);
t843 = sin(qJ(1));
t849 = cos(qJ(1));
t899 = g(1) * t849 + g(2) * t843;
t1046 = t899 * t814;
t1034 = t704 * t839;
t938 = pkin(11) * t1034;
t1007 = pkin(9) * t745;
t1013 = pkin(7) + pkin(8);
t789 = t1013 * t848;
t770 = qJD(1) * t789;
t751 = t847 * t770;
t787 = t1013 * t842;
t768 = qJD(1) * t787;
t999 = qJD(2) * pkin(2);
t753 = -t768 + t999;
t884 = -t753 * t841 - t751;
t683 = -t884 + t1007;
t675 = t840 * t683;
t741 = t746 * pkin(9);
t747 = t841 * t770;
t915 = t847 * t753 - t747;
t682 = t741 + t915;
t623 = t682 * t846 - t675;
t949 = qJD(4) * t846;
t1028 = -pkin(3) * t949 + t623;
t914 = t768 * t841 - t751;
t689 = t914 - t1007;
t962 = -t847 * t768 - t747;
t690 = t741 + t962;
t818 = pkin(2) * t847 + pkin(3);
t950 = qJD(4) * t840;
t978 = t840 * t841;
t964 = t689 * t840 + t690 * t846 - t818 * t949 - (-t841 * t950 + (t846 * t847 - t978) * qJD(3)) * pkin(2);
t1011 = pkin(3) * t746;
t885 = t745 * t840 - t846 * t746;
t657 = pkin(4) * t885 - pkin(10) * t704;
t641 = -t1011 + t657;
t821 = pkin(2) * t954;
t636 = t641 + t821;
t1045 = -t845 * t636 + t839 * t964;
t948 = qJD(5) * t839;
t1020 = (-t1034 + t948) * pkin(5);
t902 = pkin(5) * t885 - pkin(11) * t1033;
t1044 = t1028 * t839 - t845 * t641;
t943 = -qJD(5) + t704;
t1029 = t943 * t839;
t942 = qJD(1) * qJD(2);
t925 = t848 * t942;
t940 = qJDD(1) * t842;
t878 = -t925 - t940;
t941 = qJD(1) * qJD(3);
t1026 = -t848 * t941 + t878;
t833 = qJD(2) + qJD(3);
t939 = qJDD(1) * t848;
t684 = -t1026 * t847 - t833 * t927 + t841 * t939;
t760 = t841 * t842 - t847 * t848;
t762 = t841 * t848 + t842 * t847;
t716 = t833 * t762;
t855 = qJD(1) * t716;
t853 = -t760 * qJDD(1) - t855;
t607 = qJD(4) * t885 + t684 * t840 - t846 * t853;
t605 = qJDD(5) + t607;
t600 = t845 * t605;
t1043 = -t1029 * t943 + t600;
t836 = qJ(5) + qJ(6);
t828 = cos(t836);
t1002 = g(3) * t828;
t832 = qJDD(2) + qJDD(3);
t824 = qJDD(4) + t832;
t717 = qJDD(2) * pkin(2) + t1013 * t878;
t926 = t842 * t942;
t877 = -t926 + t939;
t718 = t1013 * t877;
t864 = qJD(3) * t884 + t847 * t717 - t841 * t718;
t594 = pkin(3) * t832 - pkin(9) * t684 + t864;
t952 = qJD(3) * t841;
t743 = t770 * t952;
t909 = -qJD(3) * t753 - t718;
t604 = -t743 + (pkin(9) * t1026 + t717) * t841 + ((-t842 * t941 + t877) * pkin(9) - t909) * t847;
t669 = pkin(3) * t833 + t682;
t903 = -t846 * t594 + t840 * t604 + t669 * t950 + t683 * t949;
t548 = -pkin(4) * t824 + t903;
t825 = qJD(4) + t833;
t687 = t825 * t839 + t845 * t885;
t606 = t846 * t684 + t745 * t949 + t746 * t950 + t840 * t853;
t917 = t606 * t839 - t845 * t824;
t585 = t687 * qJD(5) + t917;
t535 = pkin(5) * t585 + t548;
t618 = t669 * t846 - t675;
t615 = -pkin(4) * t825 - t618;
t685 = -t845 * t825 + t839 * t885;
t590 = pkin(5) * t685 + t615;
t815 = cos(t830);
t1042 = -t1002 * t815 + t1046 * t828 - t1049 * t590 + t535 * t759;
t602 = qJDD(6) + t605;
t694 = qJD(6) - t943;
t1041 = t1049 * t694 - t759 * t602;
t826 = sin(t836);
t1003 = g(3) * t826;
t1040 = t815 * t1003 - t826 * t1046 + t1048 * t590 + t535 * t761;
t1039 = t1048 * t694 + t761 * t602;
t599 = t839 * t605;
t1038 = t1047 * t943 + t599;
t584 = t845 * t606 + t839 * t824 + t825 * t947 - t885 * t948;
t945 = qJD(6) * t844;
t933 = t844 * t584 - t838 * t585 - t685 * t945;
t946 = qJD(6) * t838;
t541 = -t687 * t946 + t933;
t887 = t685 * t838 - t844 * t687;
t919 = t584 * t838 + t844 * t585;
t542 = -qJD(6) * t887 + t919;
t582 = t584 * t839;
t993 = t687 * t838;
t626 = t844 * t685 + t993;
t1037 = t824 * MDP(22) + t541 * t761 * MDP(32) + (-t1047 * t687 + t582) * MDP(25) + (t1029 * t687 + t1047 * t685 + t584 * t845 - t839 * t585) * MDP(26) + (-t1048 * t626 - t541 * t759 - t761 * t542) * MDP(33) + (-t1048 * MDP(32) - t1049 * MDP(33)) * t887;
t996 = t615 * t704;
t1036 = t626 * t694;
t1035 = t694 * t887;
t1004 = g(3) * t815;
t923 = -t548 - t1004;
t676 = t846 * t683;
t619 = t840 * t669 + t676;
t616 = pkin(10) * t825 + t619;
t820 = -t848 * pkin(2) - pkin(1);
t785 = t820 * qJD(1);
t719 = -pkin(3) * t745 + t785;
t633 = -pkin(4) * t704 - pkin(10) * t885 + t719;
t576 = t616 * t845 + t633 * t839;
t563 = -pkin(11) * t685 + t576;
t561 = t563 * t946;
t983 = t828 * t843;
t984 = t826 * t849;
t726 = -t815 * t983 + t984;
t982 = t828 * t849;
t985 = t826 * t843;
t728 = t815 * t982 + t985;
t1025 = g(1) * t728 - g(2) * t726 + t814 * t1002 + t590 * t626 + t561;
t725 = t815 * t985 + t982;
t727 = -t815 * t984 + t983;
t1014 = (qJD(4) * t669 + t604) * t846 + t840 * t594 - t683 * t950;
t547 = pkin(10) * t824 + t1014;
t731 = pkin(3) * t760 + t820;
t807 = pkin(2) * t926;
t664 = pkin(3) * t855 + qJDD(1) * t731 + t807;
t556 = t607 * pkin(4) - t606 * pkin(10) + t664;
t553 = t845 * t556;
t527 = pkin(5) * t605 - pkin(11) * t584 - qJD(5) * t576 - t547 * t839 + t553;
t875 = t845 * t547 + t839 * t556 - t616 * t948 + t633 * t947;
t528 = -pkin(11) * t585 + t875;
t920 = t844 * t527 - t838 * t528;
t1024 = -g(1) * t727 + g(2) * t725 + t814 * t1003 + t590 * t887 + t920;
t1023 = t602 * MDP(36) + (-t626 ^ 2 + t887 ^ 2) * MDP(33) - t626 * MDP(32) * t887;
t712 = -t760 * t840 + t762 * t846;
t658 = t761 * t712;
t622 = t682 * t840 + t676;
t901 = pkin(3) * t950 - t622;
t977 = t841 * t846;
t963 = t689 * t846 - t690 * t840 + t818 * t950 + (t841 * t949 + (t840 * t847 + t977) * qJD(3)) * pkin(2);
t961 = -t841 * t787 + t847 * t789;
t1019 = t1028 * t845 + t839 * t641;
t1018 = MDP(29) * t943 - MDP(36) * t694;
t1017 = t839 * t636 + t964 * t845;
t1016 = qJDD(1) * t820;
t1012 = -pkin(10) - pkin(11);
t1009 = pkin(3) * t846;
t1008 = pkin(5) * t845;
t806 = g(3) * t814;
t959 = pkin(2) * t977 + t840 * t818;
t740 = pkin(10) + t959;
t1001 = -pkin(11) - t740;
t816 = pkin(3) * t840 + pkin(10);
t1000 = -pkin(11) - t816;
t575 = -t616 * t839 + t845 * t633;
t562 = -pkin(11) * t687 + t575;
t558 = -pkin(5) * t943 + t562;
t998 = t558 * t844;
t997 = t563 * t844;
t711 = t846 * t760 + t762 * t840;
t715 = t833 * t760;
t637 = -qJD(4) * t711 - t715 * t846 - t716 * t840;
t995 = t637 * t839;
t994 = t637 * t845;
t987 = t712 * t839;
t986 = t712 * t845;
t980 = t839 * t843;
t979 = t839 * t849;
t976 = t843 * t845;
t913 = -t847 * t787 - t789 * t841;
t697 = -pkin(9) * t762 + t913;
t698 = -pkin(9) * t760 + t961;
t654 = t697 * t840 + t698 * t846;
t647 = t845 * t654;
t975 = t845 * t849;
t969 = t845 * t618 + t839 * t657;
t652 = pkin(4) * t711 - pkin(10) * t712 + t731;
t966 = t839 * t652 + t647;
t965 = t1020 + t963;
t960 = t1020 + t901;
t834 = t842 ^ 2;
t958 = -t848 ^ 2 + t834;
t951 = qJD(3) * t847;
t944 = t885 * MDP(18);
t823 = t842 * t999;
t935 = qJD(5) * pkin(10) * t943;
t613 = t615 * t947;
t934 = -t839 * t923 + t613;
t931 = t1046 * t845 + t615 * t948;
t819 = -pkin(4) - t1008;
t930 = qJD(2) * t1013;
t929 = qJD(5) * t1012;
t928 = t712 * t948;
t706 = pkin(3) * t716 + t823;
t922 = qJD(5) * t1001;
t921 = qJD(5) * t1000;
t912 = -pkin(2) * t978 + t818 * t846;
t907 = -qJD(5) * t633 - t547;
t906 = qJD(6) * t558 + t528;
t739 = -pkin(4) - t912;
t900 = t1020 - t619;
t898 = g(1) * t843 - g(2) * t849;
t897 = -t616 * t947 + t553;
t831 = t845 * pkin(11);
t757 = t816 * t845 + t831;
t896 = qJD(6) * t757 - t845 * t921 - t1044 + t902;
t721 = t740 * t845 + t831;
t895 = qJD(6) * t721 - t845 * t922 - t1045 + t902;
t656 = t845 * t657;
t788 = pkin(10) * t845 + t831;
t894 = qJD(6) * t788 - t618 * t839 - t845 * t929 + t656 + t902;
t756 = t1000 * t839;
t893 = -qJD(6) * t756 - t839 * t921 + t1019 - t938;
t720 = t1001 * t839;
t892 = -qJD(6) * t720 - t839 * t922 + t1017 - t938;
t786 = t1012 * t839;
t891 = -qJD(6) * t786 - t839 * t929 - t938 + t969;
t890 = -pkin(10) * t605 - t996;
t537 = t558 * t838 + t997;
t889 = -t605 * t740 - t996;
t888 = -t605 * t816 - t996;
t886 = t697 * t846 - t698 * t840;
t882 = -0.2e1 * pkin(1) * t942 - pkin(7) * qJDD(2);
t881 = t712 * t947 + t995;
t880 = -t928 + t994;
t769 = t842 * t930;
t771 = t848 * t930;
t876 = -t847 * t769 - t841 * t771 - t787 * t951 - t789 * t952;
t645 = -pkin(9) * t716 + t876;
t863 = -qJD(3) * t961 + t769 * t841 - t847 * t771;
t646 = pkin(9) * t715 + t863;
t566 = qJD(4) * t886 + t645 * t846 + t646 * t840;
t638 = qJD(4) * t712 - t715 * t840 + t846 * t716;
t573 = pkin(4) * t638 - pkin(10) * t637 + t706;
t874 = t845 * t566 + t839 * t573 + t652 * t947 - t654 * t948;
t871 = -t903 - t1004 + t1046;
t850 = qJD(2) ^ 2;
t869 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t850 + t898;
t851 = qJD(1) ^ 2;
t868 = pkin(1) * t851 - pkin(7) * qJDD(1) + t899;
t866 = -t719 * t885 + t871;
t862 = -t575 * t885 + t923 * t845 + t931;
t861 = t815 * t899 - t1014 + t806;
t536 = -t563 * t838 + t998;
t860 = -t536 * t885 + t1042;
t567 = qJD(4) * t654 + t645 * t840 - t646 * t846;
t859 = -t1046 * t839 + t576 * t885 + t934;
t858 = -t719 * t704 + t861;
t827 = sin(t837);
t829 = cos(t837);
t857 = g(3) * t827 - t841 * t717 - t785 * t745 + t829 * t899 + t847 * t909 + t743;
t856 = t537 * t885 + t1040;
t854 = -g(3) * t829 + t785 * t746 + t827 * t899 + t864;
t852 = t1037 + (-t704 ^ 2 + t885 ^ 2) * MDP(19) + (t626 * t885 + t1041) * MDP(35) + (t885 * t887 + t1039) * MDP(34) + (t825 * t885 - t607) * MDP(21) + (t685 * t885 + t1043) * MDP(28) + (-t687 * t885 + t1038) * MDP(27) + t1018 * t885 + (-t704 * t825 + t606) * MDP(20) + (-t745 * t833 + t684) * MDP(13) + (-t746 * t833 + t853) * MDP(14) + t746 * t745 * MDP(11) - t704 * t944 + (-t745 ^ 2 + t746 ^ 2) * MDP(12) + t832 * MDP(15);
t817 = -pkin(4) - t1009;
t778 = t819 - t1009;
t742 = t807 + t1016;
t737 = t815 * t975 + t980;
t736 = -t815 * t979 + t976;
t735 = -t815 * t976 + t979;
t734 = t815 * t980 + t975;
t729 = t739 - t1008;
t724 = t821 - t1011;
t659 = t759 * t712;
t644 = t845 * t652;
t614 = pkin(5) * t987 - t886;
t577 = -pkin(11) * t987 + t966;
t572 = t845 * t573;
t570 = pkin(5) * t711 - pkin(11) * t986 - t654 * t839 + t644;
t551 = t637 * t981 - t838 * t928 - t946 * t987 + (t1015 * t986 + t995) * t844;
t550 = -t1015 * t658 - t759 * t637;
t549 = pkin(5) * t881 + t567;
t531 = -pkin(11) * t881 + t874;
t529 = -pkin(11) * t994 + pkin(5) * t638 - t566 * t839 + t572 + (-t647 + (pkin(11) * t712 - t652) * t839) * qJD(5);
t1 = [(-g(1) * t725 - g(2) * t727 - t535 * t659 - t537 * t638 + t614 * t541 - t549 * t887 + t590 * t550 + t561 * t711 + (-(-qJD(6) * t577 + t529) * t694 - t570 * t602 - t527 * t711) * t838 + (-(qJD(6) * t570 + t531) * t694 - t577 * t602 - t906 * t711) * t844) * MDP(38) + (-t541 * t659 - t550 * t887) * MDP(32) + (-t541 * t658 + t542 * t659 - t550 * t626 + t551 * t887) * MDP(33) + (t541 * t711 + t550 * t694 - t602 * t659 - t638 * t887) * MDP(34) + (-(-t654 * t947 + t572) * t943 + t644 * t605 + t897 * t711 + t575 * t638 + t567 * t685 - t886 * t585 + t712 * t613 - g(1) * t735 - g(2) * t737 + (-(-qJD(5) * t652 - t566) * t943 - t654 * t605 + t907 * t711 + t548 * t712 + t615 * t637) * t839) * MDP(30) + (-g(1) * t734 - g(2) * t736 + t548 * t986 + t567 * t687 - t576 * t638 - t584 * t886 - t605 * t966 + t615 * t880 - t711 * t875 + t874 * t943) * MDP(31) + (t584 * t711 + t600 * t712 + t638 * t687 - t880 * t943) * MDP(27) + (-t585 * t711 - t599 * t712 - t638 * t685 + t881 * t943) * MDP(28) + (t605 * t711 - t638 * t943) * MDP(29) + (-t567 * t825 + t607 * t731 + t638 * t719 + t664 * t711 - t704 * t706 + t815 * t898 + t824 * t886) * MDP(23) + (-t606 * t711 - t607 * t712 + t637 * t704 - t638 * t885) * MDP(19) + (-t684 * t760 - t715 * t745 + t746 * t716 + t762 * t853) * MDP(12) + (t820 * t684 - t785 * t715 + t742 * t762 - t746 * t823 - t827 * t898 - t832 * t961 - t833 * t876) * MDP(17) + (t684 * t762 + t715 * t746) * MDP(11) + (-t715 * t833 + t762 * t832) * MDP(13) + (-t745 * t823 + t829 * t898 + t832 * t913 + t833 * t863 + (t742 + t1016) * t760 + 0.2e1 * t785 * t716) * MDP(16) + t898 * MDP(2) + t899 * MDP(3) + (t842 * t882 + t848 * t869) * MDP(9) + (-t842 * t869 + t848 * t882) * MDP(10) + ((-t685 * t845 - t687 * t839) * t637 + (-t582 - t585 * t845 + (t685 * t839 - t687 * t845) * qJD(5)) * t712) * MDP(26) + (qJDD(1) * t834 + 0.2e1 * t842 * t925) * MDP(4) + ((t529 * t844 - t531 * t838) * t694 + (t570 * t844 - t577 * t838) * t602 + t920 * t711 + t536 * t638 + t549 * t626 + t614 * t542 + t535 * t658 + t590 * t551 - g(1) * t726 - g(2) * t728 + ((-t570 * t838 - t577 * t844) * t694 - t537 * t711) * qJD(6)) * MDP(37) + 0.2e1 * (t842 * t939 - t942 * t958) * MDP(5) + (t584 * t986 + t687 * t880) * MDP(25) + (-t566 * t825 + t606 * t731 + t637 * t719 - t654 * t824 + t664 * t712 + t706 * t885 - t814 * t898) * MDP(24) + (t606 * t712 + t637 * t885) * MDP(18) + qJDD(1) * MDP(1) + (qJDD(2) * t842 + t848 * t850) * MDP(6) + (qJDD(2) * t848 - t842 * t850) * MDP(7) + (t602 * t711 + t638 * t694) * MDP(36) + (-t542 * t711 - t551 * t694 - t602 * t658 - t626 * t638) * MDP(35) + (-t638 * t825 - t711 * t824) * MDP(21) + (t637 * t825 + t712 * t824) * MDP(20) + (-t716 * t833 - t760 * t832) * MDP(14); (g(3) * t842 + t848 * t868) * MDP(10) + (-g(3) * t848 + t842 * t868) * MDP(9) + (t739 * t585 + t889 * t839 + t963 * t685 - (-t740 * t947 + t1045) * t943 + t862) * MDP(30) + (t739 * t584 + t889 * t845 + t963 * t687 - (t740 * t948 + t1017) * t943 + t859) * MDP(31) + (t962 * t833 + (t746 * t954 - t832 * t841 - t833 * t951) * pkin(2) + t857) * MDP(17) + (-t914 * t833 + (t745 * t954 + t832 * t847 - t833 * t952) * pkin(2) + t854) * MDP(16) + MDP(7) * t939 + (-(t720 * t838 + t721 * t844) * t602 + t729 * t541 + (t838 * t895 + t844 * t892) * t694 - t965 * t887 + t856) * MDP(38) + (t704 * t724 + t824 * t912 - t825 * t963 + t866) * MDP(23) + (-t724 * t885 - t824 * t959 + t825 * t964 + t858) * MDP(24) + ((t720 * t844 - t721 * t838) * t602 + t729 * t542 + (t838 * t892 - t844 * t895) * t694 + t965 * t626 + t860) * MDP(37) + MDP(6) * t940 + qJDD(2) * MDP(8) + t852 + (-MDP(4) * t842 * t848 + MDP(5) * t958) * t851; (t817 * t585 + t888 * t839 + t901 * t685 - (-t816 * t947 + t1044) * t943 + t862) * MDP(30) + (t817 * t584 + t888 * t845 + t901 * t687 - (t816 * t948 + t1019) * t943 + t859) * MDP(31) + (t622 * t825 + (-t704 * t746 + t824 * t846 - t825 * t950) * pkin(3) + t866) * MDP(23) + (t623 * t825 + (t746 * t885 - t824 * t840 - t825 * t949) * pkin(3) + t858) * MDP(24) + ((t756 * t844 - t757 * t838) * t602 + t778 * t542 + (t838 * t893 - t844 * t896) * t694 + t960 * t626 + t860) * MDP(37) + (-(t756 * t838 + t757 * t844) * t602 + t778 * t541 + (t838 * t896 + t844 * t893) * t694 - t960 * t887 + t856) * MDP(38) + (-t833 * t884 + t854) * MDP(16) + (t833 * t915 + t857) * MDP(17) + t852; t1037 + (MDP(19) * t885 + MDP(21) * t825 - MDP(23) * t719 - MDP(27) * t687 + MDP(28) * t685 - MDP(30) * t575 + MDP(31) * t576 + MDP(34) * t887 + MDP(35) * t626 - MDP(37) * t536 + MDP(38) * t537 + t1018) * t885 + (-pkin(4) * t585 - t619 * t685 + t656 * t943 + (-t618 * t943 + t890) * t839 + (t923 + t935) * t845 + t931) * MDP(30) + (-pkin(4) * t584 - t969 * t943 - t619 * t687 + t890 * t845 + (-t1046 - t935) * t839 + t934) * MDP(31) - (MDP(19) * t704 + t825 * MDP(20) + t719 * MDP(24) + t944) * t704 + t1041 * MDP(35) + ((t786 * t844 - t788 * t838) * t602 + t819 * t542 + (t838 * t891 - t844 * t894) * t694 + t900 * t626 + t1042) * MDP(37) + t1039 * MDP(34) + (-(t786 * t838 + t788 * t844) * t602 + t819 * t541 + (t838 * t894 + t844 * t891) * t694 - t900 * t887 + t1040) * MDP(38) - t607 * MDP(21) + t1043 * MDP(28) + t1038 * MDP(27) + t606 * MDP(20) + (t618 * t825 + t861) * MDP(24) + (t619 * t825 + t871) * MDP(23); t687 * t685 * MDP(25) + (-t685 ^ 2 + t687 ^ 2) * MDP(26) + (-t685 * t943 + t584) * MDP(27) + (-t917 + (-qJD(5) - t943) * t687) * MDP(28) + t605 * MDP(29) + (-g(1) * t736 + g(2) * t734 - t576 * t943 - t615 * t687 + (t907 + t806) * t839 + t897) * MDP(30) + (g(1) * t737 - g(2) * t735 - t575 * t943 + t615 * t685 + t806 * t845 - t875) * MDP(31) + (t541 + t1036) * MDP(34) + (-t542 - t1035) * MDP(35) + (-(-t562 * t838 - t997) * t694 - t537 * qJD(6) + (t602 * t844 - t626 * t687 - t694 * t946) * pkin(5) + t1024) * MDP(37) + ((-t563 * t694 - t527) * t838 + (t562 * t694 - t906) * t844 + (-t602 * t838 + t687 * t887 - t694 * t945) * pkin(5) + t1025) * MDP(38) + t1023; (t933 + t1036) * MDP(34) + (-t919 - t1035) * MDP(35) + (t537 * t694 + t1024) * MDP(37) + (-t838 * t527 - t844 * t528 + t536 * t694 + t1025) * MDP(38) + (-MDP(34) * t993 + MDP(35) * t887 - MDP(37) * t537 - MDP(38) * t998) * qJD(6) + t1023;];
tau  = t1;
