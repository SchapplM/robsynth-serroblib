% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:09:03
% EndTime: 2019-03-09 23:09:33
% DurationCPUTime: 22.30s
% Computational Cost: add. (14420->798), mult. (34799->1034), div. (0->0), fcn. (28026->14), ass. (0->335)
t796 = cos(qJ(2));
t790 = sin(pkin(6));
t998 = sin(qJ(3));
t934 = t790 * t998;
t891 = qJD(1) * t934;
t860 = t796 * t891;
t924 = qJD(3) * t998;
t1036 = t860 - t924;
t1000 = cos(qJ(3));
t797 = -pkin(10) - pkin(9);
t756 = t797 * t1000;
t793 = sin(qJ(2));
t992 = cos(pkin(6));
t916 = t992 * qJD(1);
t894 = pkin(1) * t916;
t955 = qJD(1) * t790;
t929 = t793 * t955;
t721 = -pkin(8) * t929 + t796 * t894;
t879 = pkin(2) * t793 - pkin(9) * t796;
t722 = t879 * t955;
t868 = t1000 * t722 - t721 * t998;
t932 = t796 * t1000;
t1038 = (pkin(3) * t793 - pkin(10) * t932) * t955 + t868 - t756 * qJD(3);
t755 = t797 * t998;
t961 = t1000 * t721 + t998 * t722;
t1037 = pkin(10) * t860 + t755 * qJD(3) - t961;
t1010 = qJD(3) + qJD(4);
t792 = sin(qJ(4));
t935 = t790 * t1000;
t892 = qJD(1) * t935;
t857 = t796 * t892;
t999 = cos(qJ(4));
t895 = t999 * t1000;
t964 = -t1010 * t895 + t857 * t999 + (qJD(4) * t998 - t1036) * t792;
t855 = t916 + qJD(2);
t704 = t793 * t892 + t998 * t855;
t1035 = qJD(3) * t704;
t829 = t1000 * t855;
t703 = -t793 * t891 + t829;
t843 = t792 * t703 + t704 * t999;
t1003 = t843 ^ 2;
t1011 = qJD(6) + t843;
t795 = cos(qJ(6));
t1029 = t1011 * t795;
t791 = sin(qJ(6));
t1030 = t1011 * t791;
t646 = -t999 * t703 + t704 * t792;
t928 = t796 * t955;
t754 = -qJD(3) + t928;
t744 = -qJD(4) + t754;
t983 = t744 * t791;
t618 = -t795 * t646 - t983;
t1031 = t1011 * t618;
t889 = qJD(2) * t932;
t823 = -t793 * t924 + t889;
t933 = t793 * t1000;
t1026 = qJD(1) * t823 + qJDD(1) * t933;
t904 = t992 * qJDD(1);
t850 = t904 + qJDD(2);
t806 = qJD(3) * t829 + t850 * t998;
t802 = t1026 * t790 + t806;
t931 = t796 * t998;
t898 = t790 * t931;
t859 = qJD(2) * t898;
t899 = t793 * t934;
t813 = qJD(1) * t859 + qJDD(1) * t899 - t1000 * t850 + t1035;
t574 = qJD(4) * t843 + t792 * t802 + t999 * t813;
t945 = qJDD(1) * t796;
t770 = t790 * t945;
t946 = qJD(1) * qJD(2);
t921 = t793 * t946;
t887 = t790 * t921;
t719 = qJDD(3) - t770 + t887;
t712 = qJDD(4) + t719;
t951 = qJD(6) * t795;
t938 = t791 * t574 + t646 * t951 + t795 * t712;
t952 = qJD(6) * t791;
t544 = t744 * t952 + t938;
t543 = t544 * t795;
t620 = t646 * t791 - t744 * t795;
t911 = -t795 * t574 + t712 * t791;
t545 = qJD(6) * t620 + t911;
t923 = qJD(4) * t999;
t953 = qJD(4) * t792;
t573 = -t703 * t923 + t704 * t953 + t792 * t813 - t999 * t802;
t570 = -qJDD(6) + t573;
t567 = t795 * t570;
t989 = t570 * t791;
t1034 = t1003 * MDP(19) + t712 * MDP(22) + (-t744 * t843 - t574) * MDP(21) + (-t1011 * t1029 + t989) * MDP(32) + (-t1011 * t1030 - t567) * MDP(31) + (-t1029 * t620 + (-t544 + t1031) * t791 - t795 * t545) * MDP(30) + (-t1030 * t620 + t543) * MDP(29);
t996 = pkin(5) * t646;
t991 = qJ(5) * t646;
t1033 = t1037 * t999 - t1038 * t792 + t755 * t923 + t756 * t953;
t738 = t1000 * t792 + t998 * t999;
t963 = (t1010 - t928) * t738;
t765 = pkin(8) * t928;
t724 = t793 * t894 + t765;
t849 = -pkin(3) * t1036 - t724;
t687 = pkin(9) * t855 + t724;
t852 = -pkin(2) * t796 - pkin(9) * t793 - pkin(1);
t718 = t852 * t790;
t697 = qJD(1) * t718;
t632 = t1000 * t697 - t687 * t998;
t612 = -t704 * pkin(10) + t632;
t604 = -t754 * pkin(3) + t612;
t833 = -t1000 * t687 - t697 * t998;
t613 = t703 * pkin(10) - t833;
t611 = t999 * t613;
t565 = t792 * t604 + t611;
t562 = qJ(5) * t744 - t565;
t548 = -t562 - t996;
t1032 = t1011 * t548;
t974 = t792 * t613;
t564 = -t999 * t604 + t974;
t950 = -qJD(5) - t564;
t794 = sin(qJ(1));
t1001 = cos(qJ(1));
t882 = t992 * t1001;
t729 = t793 * t882 + t794 * t796;
t789 = qJ(3) + qJ(4);
t784 = sin(t789);
t785 = cos(t789);
t936 = t790 * t1001;
t669 = t729 * t784 + t785 * t936;
t728 = t793 * t794 - t796 * t882;
t1028 = t669 * t791 + t728 * t795;
t1027 = t669 * t795 - t728 * t791;
t1002 = pkin(4) + pkin(11);
t1021 = pkin(5) * t843;
t949 = t1021 - t950;
t546 = t1002 * t744 + t949;
t686 = -pkin(2) * t855 - t721;
t651 = -t703 * pkin(3) + t686;
t805 = -qJ(5) * t843 + t651;
t559 = t1002 * t646 + t805;
t536 = t546 * t795 - t559 * t791;
t537 = t546 * t791 + t559 * t795;
t581 = t646 * pkin(4) + t805;
t1023 = MDP(18) * t843 - MDP(19) * t646 + t651 * MDP(24) - t581 * MDP(27) + t620 * MDP(31) - t618 * MDP(32) + MDP(33) * t1011 + t536 * MDP(34) - t537 * MDP(35);
t1022 = pkin(4) * t843;
t1020 = t548 * t843;
t576 = t792 * t612 + t611;
t883 = pkin(3) * t953 - t576;
t577 = t612 * t999 - t974;
t971 = -pkin(3) * t923 - qJD(5) + t577;
t967 = qJ(5) * t929 - t1033;
t692 = t792 * t755 - t756 * t999;
t1019 = qJD(4) * t692 + t1037 * t792 + t1038 * t999;
t702 = t712 * qJ(5);
t1018 = qJD(5) * t744 - t702;
t917 = t796 * t992;
t978 = t790 * t793;
t847 = pkin(1) * t917 - pkin(8) * t978;
t1017 = qJ(5) * t964 - qJD(5) * t738 + t849;
t1016 = (qJDD(2) + 0.2e1 * t904) * t790;
t925 = qJD(3) * t1000;
t1015 = -t857 + t925;
t1014 = -t1000 * t719 - t754 * t924;
t730 = t1001 * t793 + t794 * t917;
t1013 = g(1) * t730 + g(2) * t728;
t1012 = t1002 * t843;
t893 = pkin(1) * qJD(2) * t992;
t851 = qJD(1) * t893;
t886 = pkin(1) * t904;
t943 = t793 * qJDD(1);
t919 = t790 * t943;
t900 = pkin(8) * t919 + qJD(2) * t765 + t793 * t851 - t796 * t886;
t657 = -pkin(2) * t850 + t900;
t1009 = t657 - t1013;
t918 = t793 * t992;
t731 = t1001 * t796 - t794 * t918;
t977 = t790 * t794;
t673 = t731 * t784 - t785 * t977;
t714 = t784 * t978 - t785 * t992;
t937 = pkin(8) * t770 + t793 * t886 + t796 * t851;
t826 = -pkin(8) * t887 + t937;
t656 = pkin(9) * t850 + t826;
t844 = t879 * qJD(2);
t658 = (qJD(1) * t844 + qJDD(1) * t852) * t790;
t872 = t1000 * t658 - t656 * t998;
t551 = t719 * pkin(3) - pkin(10) * t802 - t687 * t925 - t697 * t924 + t872;
t832 = t1000 * t656 + t998 * t658 - t687 * t924 + t697 * t925;
t558 = -pkin(10) * t813 + t832;
t901 = -t999 * t551 + t792 * t558 + t604 * t953 + t613 * t923;
t820 = g(1) * t673 + g(2) * t669 + g(3) * t714 - t901;
t810 = t581 * t843 + qJDD(5) - t820;
t781 = -pkin(3) * t999 - pkin(4);
t778 = -pkin(11) + t781;
t1008 = t1011 * (t883 + t996) - t778 * t570;
t1007 = -t651 * t843 + t820;
t976 = t790 * t796;
t830 = g(3) * t976 - t1013;
t1006 = -t692 * t712 + t830 * t784;
t691 = -t755 * t999 - t792 * t756;
t1005 = -t691 * t712 - t830 * t785;
t959 = pkin(1) * t918 + pkin(8) * t976;
t717 = pkin(9) * t992 + t959;
t881 = t992 * t1000;
t822 = t899 - t881;
t809 = -qJD(3) * t822 + t790 * t889;
t723 = t790 * t844;
t725 = t847 * qJD(2);
t871 = t1000 * t723 - t725 * t998;
t954 = qJD(2) * t793;
t927 = t790 * t954;
t587 = pkin(3) * t927 - pkin(10) * t809 - t717 * t925 - t718 * t924 + t871;
t727 = t790 * t933 + t992 * t998;
t808 = qJD(3) * t727 + t859;
t831 = t1000 * t725 - t717 * t924 + t718 * t925 + t998 * t723;
t591 = -pkin(10) * t808 + t831;
t869 = t1000 * t718 - t717 * t998;
t617 = -pkin(3) * t976 - t727 * pkin(10) + t869;
t962 = t1000 * t717 + t998 * t718;
t625 = -pkin(10) * t822 + t962;
t838 = -t792 * t587 - t999 * t591 - t617 * t923 + t625 * t953;
t538 = -t790 * (qJ(5) * t954 - qJD(5) * t796) + t838;
t799 = qJD(1) ^ 2;
t997 = pkin(4) * t712;
t993 = g(3) * t790;
t902 = -t792 * t551 - t999 * t558 - t604 * t923 + t613 * t953;
t531 = t902 + t1018;
t529 = -pkin(5) * t574 - t531;
t528 = t529 * t795;
t990 = t565 * t744;
t737 = t792 * t998 - t895;
t984 = t737 * t791;
t981 = t784 * t791;
t980 = t784 * t795;
t786 = t790 ^ 2;
t979 = t786 * t799;
t975 = t791 * t796;
t973 = t795 * t796;
t972 = t529 * t791 + t548 * t951;
t970 = -pkin(5) * t963 - t967;
t969 = pkin(4) * t963 + t1017;
t968 = t792 * t617 + t999 * t625;
t966 = pkin(4) * t929 + t1019;
t926 = qJD(2) * t976;
t726 = pkin(8) * t926 + t793 * t893;
t960 = -t971 + t1021;
t787 = t793 ^ 2;
t958 = -t796 ^ 2 + t787;
t941 = t796 * t979;
t940 = t790 * t975;
t939 = t790 * t973;
t922 = 0.2e1 * pkin(1) * t786;
t920 = t796 * t946;
t670 = t729 * t785 - t784 * t936;
t915 = -t669 * pkin(4) + t670 * qJ(5);
t674 = t731 * t785 + t784 * t977;
t914 = -t673 * pkin(4) + t674 * qJ(5);
t715 = t784 * t992 + t785 * t978;
t913 = -t714 * pkin(4) + t715 * qJ(5);
t853 = qJDD(5) + t901;
t526 = -pkin(5) * t573 - t1002 * t712 + t853;
t598 = pkin(3) * t813 + t657;
t803 = t573 * qJ(5) - qJD(5) * t843 + t598;
t530 = t1002 * t574 + t803;
t912 = t795 * t526 - t791 * t530;
t910 = t791 * t929 + t795 * t963;
t909 = -t791 * t963 + t795 * t929;
t903 = t1002 * t978;
t782 = pkin(3) * t1000 + pkin(2);
t896 = t1001 * t998;
t880 = t790 * t799 * t992;
t878 = g(1) * t669 - g(2) * t673;
t877 = -g(1) * t670 + g(2) * t674;
t876 = -g(1) * t728 + g(2) * t730;
t875 = g(1) * t731 + g(2) * t729;
t874 = t617 * t999 - t792 * t625;
t870 = t1000 * t729 - t790 * t896;
t840 = -t738 * qJ(5) - t782;
t652 = t1002 * t737 + t840;
t867 = pkin(5) * t964 - qJD(1) * t903 + qJD(6) * t652 - t1019;
t661 = t738 * pkin(5) + t691;
t866 = -qJD(6) * t661 - t1002 * t963 - t1017;
t865 = pkin(3) * t704 + t991;
t864 = t791 * t526 + t795 * t530;
t863 = (t565 - t996) * t1011 - t1002 * t570;
t583 = pkin(4) * t976 - t874;
t818 = t792 * t822;
t664 = t727 * t999 - t818;
t560 = t664 * pkin(5) + pkin(11) * t976 + t583;
t816 = t999 * t822;
t663 = t727 * t792 + t816;
t716 = -pkin(2) * t992 - t847;
t819 = t822 * pkin(3);
t665 = t819 + t716;
t804 = -t664 * qJ(5) + t665;
t579 = t1002 * t663 + t804;
t862 = t560 * t795 - t579 * t791;
t861 = t560 * t791 + t579 * t795;
t854 = 0.2e1 * t916 + qJD(2);
t582 = qJ(5) * t976 - t968;
t641 = t663 * t795 + t940;
t839 = -t587 * t999 + t792 * t591 + t617 * t953 + t625 * t923;
t837 = -g(1) * t674 - g(2) * t670 - g(3) * t715;
t836 = t850 * MDP(8);
t835 = t737 * t952 - t910;
t834 = t737 * t951 - t909;
t682 = -t731 * t998 + t794 * t935;
t827 = -t719 * t998 + t754 * t925;
t825 = -t646 * t744 - t573;
t824 = t1001 * t935 + t729 * t998;
t821 = -t837 + t902;
t815 = (-qJD(6) * t778 + t1012 + t865) * t1011 + t837;
t814 = (qJD(6) * t1002 + t1012 + t991) * t1011 + t837;
t812 = -t821 - t1018;
t811 = qJD(3) * t833 + t872;
t653 = pkin(3) * t808 + t726;
t596 = qJD(4) * t816 + t727 * t953 + t792 * t808 - t809 * t999;
t597 = -qJD(4) * t818 + t727 * t923 + t792 * t809 + t808 * t999;
t541 = t597 * pkin(4) + t596 * qJ(5) - t664 * qJD(5) + t653;
t779 = pkin(3) * t792 + qJ(5);
t683 = t1000 * t731 + t794 * t934;
t666 = t737 * pkin(4) + t840;
t662 = -t737 * pkin(5) + t692;
t642 = t663 * t791 - t939;
t637 = t673 * t791 + t730 * t795;
t636 = t673 * t795 - t730 * t791;
t599 = t991 + t1022;
t592 = t663 * pkin(4) + t804;
t588 = t865 + t1022;
t572 = qJD(6) * t641 + t597 * t791 + t795 * t927;
t571 = -t795 * t597 - qJD(6) * t939 + (qJD(6) * t663 + t927) * t791;
t563 = -pkin(5) * t663 - t582;
t561 = pkin(4) * t744 - t950;
t540 = -pkin(4) * t927 + t839;
t539 = t597 * pkin(11) + t541;
t535 = t574 * pkin(4) + t803;
t534 = -pkin(5) * t597 - t538;
t533 = -t596 * pkin(5) - qJD(2) * t903 + t839;
t532 = t853 - t997;
t524 = -t537 * qJD(6) + t912;
t523 = qJD(6) * t536 + t864;
t1 = [(-t719 * t796 - t754 * t954) * t790 * MDP(15) + (-t712 * t796 - t744 * t954) * t790 * MDP(22) + (t703 * t809 - t704 * t808 - t727 * t813 - t802 * t822) * MDP(12) + (-g(1) * t824 - g(2) * t682 + t657 * t727 + t686 * t809 + t726 * t704 + t716 * t802 - t719 * t962 + t754 * t831 + t832 * t976 + t833 * t927) * MDP(17) + (-t1011 * t596 - t570 * t664) * MDP(33) + (-t1011 * t571 - t545 * t664 - t570 * t641 + t596 * t618) * MDP(32) + (t1011 * t572 + t544 * t664 - t570 * t642 - t596 * t620) * MDP(31) + (t535 * t592 + t581 * t541 + t531 * t582 + t562 * t538 + t532 * t583 + t561 * t540 - g(1) * (-t794 * pkin(1) - t670 * pkin(4) - t669 * qJ(5) + t728 * t797 - t729 * t782 + (pkin(3) * t896 + pkin(8) * t1001) * t790) - g(2) * (t1001 * pkin(1) + t674 * pkin(4) + t673 * qJ(5) - t730 * t797 + t731 * t782 + (pkin(3) * t998 + pkin(8)) * t977)) * MDP(28) + (t531 * t663 + t532 * t664 + t538 * t646 + t540 * t843 - t561 * t596 + t562 * t597 - t573 * t583 + t574 * t582 - t876) * MDP(25) + (-t573 * t664 - t596 * t843) * MDP(18) + (t573 * t663 - t574 * t664 + t596 * t646 - t597 * t843) * MDP(19) + (-t535 * t664 + t538 * t744 - t541 * t843 + t573 * t592 + t581 * t596 - t582 * t712 + (t531 * t796 - t562 * t954) * t790 + t878) * MDP(27) + (t596 * t744 + t664 * t712 + (t573 * t796 + t843 * t954) * t790) * MDP(20) + (-t838 * t744 - t968 * t712 + t653 * t843 - t665 * t573 + t598 * t664 - t651 * t596 + (-t565 * t954 - t796 * t902) * t790 - t878) * MDP(24) + (-(qJD(6) * t862 + t533 * t791 + t539 * t795) * t1011 + t861 * t570 - t523 * t664 + t537 * t596 + t534 * t620 + t563 * t544 + t529 * t642 + t548 * t572 + g(1) * t1027 - g(2) * t636) * MDP(35) + ((-qJD(6) * t861 + t533 * t795 - t539 * t791) * t1011 - t862 * t570 + t524 * t664 - t536 * t596 + t534 * t618 + t563 * t545 - t529 * t641 + t548 * t571 + g(1) * t1028 - g(2) * t637) * MDP(34) + (t1016 * t793 + t854 * t926) * MDP(6) + (t1016 * t796 - t854 * t927) * MDP(7) + (-(-qJD(3) * t962 + t871) * t754 + t869 * t719 - t811 * t976 + t632 * t927 - t726 * t703 + t716 * t813 + t657 * t822 + t686 * t808 + g(1) * t870 - g(2) * t683) * MDP(16) + qJDD(1) * MDP(1) + (t544 * t641 - t545 * t642 - t571 * t620 - t572 * t618) * MDP(30) + (t544 * t642 + t572 * t620) * MDP(29) + (g(1) * t794 - g(2) * t1001) * MDP(2) + (g(1) * t1001 + g(2) * t794) * MDP(3) + (0.2e1 * (t796 * t943 - t946 * t958) * MDP(5) + (qJDD(1) * t787 + 0.2e1 * t793 * t920) * MDP(4)) * t786 + t992 * t836 + (-t726 * t855 + t847 * t850 - t900 * t992 + g(1) * t729 - g(2) * t731 + (-t921 + t945) * t922) * MDP(9) + (-t725 * t855 - t959 * t850 - t826 * t992 + (-t920 - t943) * t922 + t876) * MDP(10) + (t839 * t744 + t874 * t712 + t653 * t646 + t665 * t574 + t598 * t663 + t651 * t597 + (-t564 * t954 + t796 * t901) * t790 - t877) * MDP(23) + (t597 * t744 - t663 * t712 + (t574 * t796 - t646 * t954) * t790) * MDP(21) + (-t535 * t663 - t540 * t744 - t541 * t646 - t574 * t592 - t581 * t597 + t583 * t712 + (-t532 * t796 + t561 * t954) * t790 + t877) * MDP(26) + (t806 * t727 + t881 * t1035 + (t1026 * t727 + t704 * t823) * t790) * MDP(11) + (t704 * t927 + t727 * t719 - t754 * t809 - t802 * t976) * MDP(13) + (t703 * t927 - t719 * t822 + t754 * t808 + t813 * t976) * MDP(14); (t793 * t880 + t770) * MDP(7) + (-pkin(2) * t813 + t686 * t924 + t868 * t754 + t724 * t703 + t827 * pkin(9) + (-g(3) * t932 + (-t632 * t793 - t686 * t931) * qJD(1)) * t790 - t1009 * t1000) * MDP(16) + (-t535 * t737 - t574 * t666 - t581 * t963 - t646 * t969 - t744 * t966 - t1005) * MDP(26) + (t1019 * t744 - t782 * t574 + t598 * t737 + t849 * t646 + t963 * t651 + t1005) * MDP(23) + (-pkin(2) * t802 + pkin(9) * t1014 + g(3) * t898 + t1009 * t998 + t1015 * t686 - t724 * t704 - t754 * t961) * MDP(17) - t793 * MDP(4) * t941 + t958 * MDP(5) * t979 + (-t535 * t738 + t573 * t666 + t581 * t964 + t744 * t967 - t843 * t969 - t1006) * MDP(27) + (t754 * MDP(15) - MDP(17) * t833 - MDP(20) * t843 + t646 * MDP(21) + t744 * MDP(22) + t564 * MDP(23) + t565 * MDP(24) - t561 * MDP(26) + t562 * MDP(27)) * t929 + (-(-t652 * t791 + t661 * t795) * t570 + t524 * t738 + t662 * t545 - t737 * t528 - g(1) * (-t730 * t981 + t731 * t795) - g(2) * (-t728 * t981 + t729 * t795) - (t784 * t975 + t793 * t795) * t993 + (t791 * t866 - t795 * t867) * t1011 + t970 * t618 - t964 * t536 + t835 * t548) * MDP(34) + ((t652 * t795 + t661 * t791) * t570 - t523 * t738 + t662 * t544 + t529 * t984 - g(1) * (-t730 * t980 - t731 * t791) - g(2) * (-t728 * t980 - t729 * t791) - (t784 * t973 - t791 * t793) * t993 + (t791 * t867 + t795 * t866) * t1011 + t970 * t620 + t964 * t537 + t834 * t548) * MDP(35) + (-t1011 * t964 - t570 * t738) * MDP(33) + (-t1011 * t835 - t545 * t738 - t567 * t737 + t618 * t964) * MDP(32) + (t1011 * t834 + t544 * t738 - t570 * t984 - t620 * t964) * MDP(31) + t836 + (t573 * t737 - t574 * t738 + t646 * t964 - t843 * t963) * MDP(19) + (-t573 * t738 - t843 * t964) * MDP(18) + (-g(3) * t978 + t531 * t737 + t532 * t738 - t561 * t964 + t562 * t963 - t573 * t691 - t574 * t692 + t646 * t967 + t843 * t966 - t875) * MDP(25) + (t1000 * t802 + t1015 * t703 + t1036 * t704 - t813 * t998) * MDP(12) + (t1033 * t744 + t782 * t573 + t598 * t738 - t964 * t651 + t849 * t843 + t1006) * MDP(24) + (-t531 * t692 + t532 * t691 + t535 * t666 + t966 * t561 + t967 * t562 + t969 * t581 + (t793 * t993 + t875) * t797 + (-t796 * t993 + t1013) * (pkin(4) * t785 + qJ(5) * t784 + t782)) * MDP(28) + ((-t703 * t793 - t754 * t931) * t955 - t1014) * MDP(14) + (t1015 * t704 + t802 * t998) * MDP(11) + (t712 * t738 + t744 * t964) * MDP(20) + (-t712 * t737 + t744 * t963) * MDP(21) + (-t796 * t880 + t919) * MDP(6) + ((-t704 * t793 + t754 * t932) * t955 - t827) * MDP(13) + (pkin(1) * t941 + t721 * t855 + (pkin(8) * t946 + g(3)) * t978 + t875 - t937) * MDP(10) + (pkin(1) * t793 * t979 + t724 * t855 - t830 - t900) * MDP(9) + (t544 * t984 + t620 * t834) * MDP(29) + (t910 * t620 + t909 * t618 + (t543 - t545 * t791 + (-t618 * t795 - t620 * t791) * qJD(6)) * t737) * MDP(30); (-t704 * t754 - t813) * MDP(14) - t573 * MDP(20) + (t588 * t646 - t883 * t744 + (-pkin(4) + t781) * t712 + t810) * MDP(26) + (-g(1) * t682 + g(2) * t824 + g(3) * t822 - t686 * t704 + t754 * t833 + t811) * MDP(16) + (g(1) * t683 + g(2) * t870 + g(3) * t727 - t632 * t754 - t686 * t703 - t832) * MDP(17) + (-MDP(20) * t744 + MDP(25) * t561 + t1023) * t646 + (t703 * t754 + t802) * MDP(13) + (-t576 * t744 + (-t646 * t704 + t712 * t999 + t744 * t953) * pkin(3) + t1007) * MDP(23) + (-t573 * t781 - t574 * t779 + t646 * t971 + (-t562 + t883) * t843) * MDP(25) + (t779 * t544 + t528 + t960 * t620 + t815 * t795 + (-t1008 - t1032) * t791) * MDP(35) + (-t531 * t779 + t532 * t781 - t581 * t588 - g(1) * (pkin(3) * t682 + t914) - g(2) * (-pkin(3) * t824 + t915) - g(3) * (-t819 + t913) + t971 * t562 + t883 * t561) * MDP(28) + t1034 + (t588 * t843 + t712 * t779 + t744 * t971 + t812) * MDP(27) + t719 * MDP(15) + (-t703 ^ 2 + t704 ^ 2) * MDP(12) - t704 * t703 * MDP(11) + (t779 * t545 + t960 * t618 + (t1008 + t1020) * t795 + t815 * t791 + t972) * MDP(34) + (-t577 * t744 + (-t704 * t843 - t712 * t792 + t744 * t923) * pkin(3) + t821) * MDP(24); t825 * MDP(20) + (t1007 - t990) * MDP(23) + (t564 * t744 + t821) * MDP(24) + (pkin(4) * t573 - qJ(5) * t574 + (-t562 - t565) * t843) * MDP(25) + (t810 + t990 - 0.2e1 * t997) * MDP(26) + (t599 * t843 + t744 * t950 + t702 + t812) * MDP(27) + (-t532 * pkin(4) - g(1) * t914 - g(2) * t915 - g(3) * t913 - t531 * qJ(5) - t561 * t565 + t562 * t950 - t581 * t599) * MDP(28) + (qJ(5) * t545 + t949 * t618 + (-t863 + t1020) * t795 + t814 * t791 + t972) * MDP(34) + (qJ(5) * t544 + t528 + t949 * t620 + (t863 - t1032) * t791 + t814 * t795) * MDP(35) + ((t561 + t950) * MDP(25) + t599 * MDP(26) + t1023) * t646 + t1034; t825 * MDP(25) + (-t646 * t843 + t712) * MDP(26) + (-t744 ^ 2 - t1003) * MDP(27) + (-t562 * t744 + t810 - t997) * MDP(28) + (t618 * t744 - t567) * MDP(34) + (t620 * t744 + t989) * MDP(35) + (-MDP(34) * t1030 - MDP(35) * t1029) * t1011; t620 * t618 * MDP(29) + (-t618 ^ 2 + t620 ^ 2) * MDP(30) + (t938 + t1031) * MDP(31) + (t1011 * t620 - t911) * MDP(32) - t570 * MDP(33) + (t537 * t1011 - t548 * t620 - g(1) * t636 - g(2) * t1027 - g(3) * (t714 * t795 + t940) + t912) * MDP(34) + (t536 * t1011 + t548 * t618 + g(1) * t637 + g(2) * t1028 - g(3) * (-t714 * t791 + t939) - t864) * MDP(35) + (MDP(31) * t983 - MDP(32) * t620 - MDP(34) * t537 - MDP(35) * t536) * qJD(6);];
tau  = t1;
