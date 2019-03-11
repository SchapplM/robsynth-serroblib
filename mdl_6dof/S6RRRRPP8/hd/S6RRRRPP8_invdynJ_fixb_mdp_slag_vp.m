% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:39:00
% EndTime: 2019-03-09 21:39:26
% DurationCPUTime: 20.19s
% Computational Cost: add. (13948->894), mult. (33861->1097), div. (0->0), fcn. (26704->10), ass. (0->336)
t856 = cos(qJ(2));
t850 = sin(pkin(6));
t998 = qJD(1) * t850;
t832 = t856 * t998;
t853 = sin(qJ(2));
t1053 = cos(pkin(6));
t954 = t1053 * qJD(1);
t933 = pkin(1) * t954;
t773 = pkin(8) * t832 + t853 * t933;
t852 = sin(qJ(3));
t855 = cos(qJ(3));
t925 = pkin(3) * t852 - pkin(10) * t855;
t677 = t832 * t925 + t773;
t798 = t925 * qJD(3);
t1087 = t677 - t798;
t914 = t954 + qJD(2);
t965 = t853 * t998;
t737 = t852 * t965 - t855 * t914;
t728 = qJD(4) + t737;
t892 = qJD(3) * t914;
t942 = t1053 * qJDD(1);
t910 = t942 + qJDD(2);
t996 = qJD(2) * t856;
t963 = t852 * t996;
t982 = qJDD(1) * t853;
t993 = qJD(3) * t855;
t637 = (qJD(1) * (t853 * t993 + t963) + t852 * t982) * t850 + t852 * t892 - t855 * t910;
t631 = qJDD(4) + t637;
t935 = t852 * t832;
t995 = qJD(3) * t852;
t1086 = t935 - t995;
t739 = t852 * t914 + t855 * t965;
t851 = sin(qJ(4));
t854 = cos(qJ(4));
t917 = t832 - qJD(3);
t681 = t854 * t739 - t851 * t917;
t1037 = t681 * t728;
t891 = t852 * t910;
t962 = t855 * t996;
t981 = qJDD(1) * t855;
t861 = t855 * t892 + (t853 * t981 + (-t853 * t995 + t962) * qJD(1)) * t850 + t891;
t983 = qJD(1) * qJD(2);
t958 = t853 * t983;
t932 = t850 * t958;
t980 = qJDD(1) * t856;
t831 = t850 * t980;
t978 = qJDD(3) - t831;
t889 = t932 + t978;
t992 = qJD(4) * t681;
t594 = t851 * t861 - t854 * t889 + t992;
t1085 = -t594 - t1037;
t1084 = -pkin(3) * t855 - pkin(2);
t1050 = qJ(5) * t851;
t1083 = pkin(3) + t1050;
t1082 = t631 * qJ(5) + t728 * qJD(5);
t770 = -pkin(8) * t965 + t856 * t933;
t926 = pkin(2) * t853 - pkin(9) * t856;
t771 = t926 * t998;
t1008 = t855 * t770 + t852 * t771;
t651 = pkin(10) * t965 + t1008;
t807 = -pkin(10) * t852 + t1084;
t1023 = t854 * t855;
t842 = pkin(9) * t1023;
t991 = qJD(4) * t851;
t1081 = -qJD(4) * t842 + t851 * t651 - t807 * t991;
t990 = qJD(4) * t854;
t1080 = t1087 * t851 + t854 * t651 - t807 * t990;
t1060 = sin(qJ(1));
t1061 = cos(qJ(1));
t929 = t1053 * t1061;
t788 = t1060 * t856 + t853 * t929;
t967 = t850 * t1061;
t703 = t788 * t855 - t852 * t967;
t787 = t1060 * t853 - t856 * t929;
t652 = t703 * t851 - t787 * t854;
t653 = t703 * t854 + t787 * t851;
t945 = t854 * t917;
t679 = t739 * t851 + t945;
t1079 = t594 * qJ(6) + t679 * qJD(6);
t1078 = 0.2e1 * t1082;
t713 = pkin(9) * t914 + t773;
t912 = -pkin(2) * t856 - pkin(9) * t853 - pkin(1);
t725 = t912 * t998;
t636 = t855 * t713 + t852 * t725;
t1077 = qJD(5) * t851 + t636;
t1026 = t850 * t856;
t1027 = t850 * t853;
t785 = t1027 * t852 - t1053 * t855;
t699 = -qJD(3) * t785 + t850 * t962;
t1076 = -qJD(4) * t1026 + t699;
t968 = pkin(1) * t1053;
t907 = -pkin(8) * t1027 + t856 * t968;
t1074 = -qJ(5) * t1086 - qJD(5) * t855 - t1080;
t1073 = (qJDD(2) + 0.2e1 * t942) * t850;
t934 = t855 * t832;
t1072 = t934 - t993;
t712 = -pkin(2) * t914 - t770;
t617 = t737 * pkin(3) - t739 * pkin(10) + t712;
t621 = -pkin(10) * t917 + t636;
t579 = t854 * t617 - t851 * t621;
t985 = qJD(5) - t579;
t1071 = t677 * t854 - t1081;
t1039 = t679 * t728;
t593 = qJD(4) * t945 + t739 * t991 - t851 * t889 - t854 * t861;
t1070 = -t593 - t1039;
t628 = t631 * pkin(4);
t1069 = t628 - qJDD(5);
t1048 = qJ(6) * t593;
t1062 = pkin(4) + pkin(5);
t635 = -t852 * t713 + t855 * t725;
t620 = pkin(3) * t917 - t635;
t883 = -t681 * qJ(5) + t620;
t561 = -t1062 * t679 + qJD(6) - t883;
t928 = t1053 * t1060;
t790 = t1061 * t856 - t853 * t928;
t966 = t850 * t1060;
t707 = t790 * t855 + t852 * t966;
t789 = t1061 * t853 + t856 * t928;
t656 = t707 * t851 - t789 * t854;
t1022 = t854 * t856;
t786 = t1027 * t855 + t1053 * t852;
t700 = t1022 * t850 + t786 * t851;
t911 = qJD(2) * t933;
t930 = pkin(1) * t942;
t971 = pkin(8) * t831 + t853 * t930 + t856 * t911;
t882 = -pkin(8) * t932 + t971;
t664 = pkin(9) * t910 + t882;
t904 = t926 * qJD(2);
t667 = (qJD(1) * t904 + qJDD(1) * t912) * t850;
t900 = -t855 * t664 - t852 * t667 + t713 * t995 - t725 * t993;
t569 = pkin(10) * t889 - t900;
t957 = t856 * t983;
t931 = t850 * t957;
t956 = t850 * t982;
t937 = t853 * t911 - t856 * t930 + (t931 + t956) * pkin(8);
t665 = -pkin(2) * t910 + t937;
t575 = t637 * pkin(3) - pkin(10) * t861 + t665;
t940 = t851 * t569 - t854 * t575 + t617 * t991 + t621 * t990;
t874 = g(1) * t656 + g(2) * t652 + g(3) * t700 - t940;
t871 = t874 + t1069;
t1068 = (qJD(6) + t561) * t681 - t1048 + t871;
t1056 = pkin(10) * t631;
t578 = t679 * pkin(4) + t883;
t1066 = t578 * t728 - t1056;
t1049 = qJ(5) * t854;
t1065 = t1062 * t851 - t1049;
t1064 = -t1062 * t854 - t1050;
t1063 = t679 ^ 2;
t678 = t681 ^ 2;
t858 = qJD(1) ^ 2;
t1058 = pkin(5) * t631;
t1057 = pkin(5) * t700;
t1055 = pkin(10) - qJ(6);
t1054 = pkin(10) * qJD(4);
t1052 = qJ(5) * t594;
t1051 = qJ(5) * t679;
t1047 = qJ(6) * t737;
t1046 = qJ(6) * t852;
t580 = t851 * t617 + t854 * t621;
t563 = qJ(6) * t679 + t580;
t719 = t728 * qJ(5);
t560 = t563 + t719;
t1045 = t560 * t728;
t572 = t719 + t580;
t1044 = t572 * t728;
t1043 = t580 * t728;
t1042 = t593 * t851;
t1041 = t631 * t851;
t1040 = t631 * t854;
t1038 = t681 * t679;
t948 = -t788 * t852 - t855 * t967;
t1036 = t948 * t854;
t706 = t790 * t852 - t855 * t966;
t1035 = t706 * t854;
t1034 = t728 * t851;
t776 = t785 * qJ(5);
t1033 = t785 * t854;
t1031 = t787 * t852;
t1029 = t789 * t852;
t847 = t850 ^ 2;
t1028 = t847 * t858;
t1025 = t851 * t855;
t1024 = t852 * t854;
t741 = (t1022 * t855 + t851 * t853) * t850;
t724 = qJD(1) * t741;
t969 = -pkin(9) * t851 - pkin(4);
t987 = qJD(6) * t854;
t1021 = qJ(6) * t724 + t1062 * t935 + (-qJ(6) * t993 - t798) * t854 + (qJ(6) * t991 - t987 + (-pkin(5) + t969) * qJD(3)) * t852 + t1071;
t723 = t851 * t934 - t854 * t965;
t1020 = -qJ(6) * t723 + (-pkin(9) * qJD(3) + qJ(6) * qJD(4)) * t1024 + (qJD(6) * t852 + (-pkin(9) * qJD(4) + qJ(6) * qJD(3)) * t855) * t851 + t1074;
t1009 = -t852 * t770 + t855 * t771;
t903 = pkin(3) * t965 + t1009;
t879 = qJ(5) * t724 + t903;
t897 = -pkin(9) - t1065;
t988 = qJD(5) * t854;
t1019 = t1062 * t723 + (qJD(4) * t1064 + t988) * t852 + t897 * t993 - t879;
t666 = pkin(3) * t739 + pkin(10) * t737;
t1018 = t854 * t635 + t851 * t666;
t994 = qJD(3) * t854;
t1017 = (-t852 * t994 - t855 * t991) * pkin(9) + t1074;
t760 = -pkin(2) * t1053 - t907;
t777 = t785 * pkin(3);
t642 = -t786 * pkin(10) + t760 + t777;
t1001 = pkin(8) * t1026 + t853 * t968;
t761 = pkin(9) * t1053 + t1001;
t1002 = -pkin(2) * t1026 - pkin(9) * t1027;
t762 = -pkin(1) * t850 + t1002;
t1011 = t855 * t761 + t852 * t762;
t644 = -pkin(10) * t1026 + t1011;
t1016 = t851 * t642 + t854 * t644;
t1014 = pkin(4) * t935 - t798 * t854 + t969 * t995 + t1071;
t918 = pkin(4) * t851 - t1049;
t909 = pkin(9) + t918;
t919 = pkin(4) * t854 + t1050;
t1013 = -pkin(4) * t723 + (qJD(4) * t919 - t988) * t852 + t909 * t993 + t879;
t1012 = -t852 * t761 + t855 * t762;
t1010 = -t1065 * t728 + t1077;
t584 = t739 * qJ(5) + t1018;
t1007 = t1047 * t851 - t1055 * t991 - t584 - t987;
t1006 = t728 * t918 - t1077;
t624 = t851 * t635;
t815 = t1055 * t854;
t1005 = qJD(4) * t815 - qJD(6) * t851 - t624 - (-t666 + t1047) * t854 + t1062 * t739;
t1000 = t851 * t807 + t842;
t848 = t853 ^ 2;
t999 = -t856 ^ 2 + t848;
t997 = qJD(2) * t853;
t562 = qJ(6) * t681 + t579;
t986 = qJD(5) - t562;
t836 = pkin(3) * t1026;
t977 = t856 * t1028;
t976 = t852 * t1026;
t975 = t851 * t1026;
t587 = t776 + t1016;
t974 = pkin(4) * t1036 + t1083 * t948;
t973 = -pkin(4) * t1035 - t1083 * t706;
t643 = t836 - t1012;
t972 = -pkin(4) * t1033 - t851 * t776 - t777;
t970 = -g(1) * t1029 - g(2) * t1031 + g(3) * t976;
t964 = t850 * t997;
t961 = t728 * t991;
t959 = 0.2e1 * pkin(1) * t847;
t953 = -t652 * pkin(4) + qJ(5) * t653;
t657 = t707 * t854 + t789 * t851;
t952 = -t656 * pkin(4) + qJ(5) * t657;
t701 = t786 * t854 - t975;
t951 = -t700 * pkin(4) + qJ(5) * t701;
t950 = t642 * t854 - t851 * t644;
t841 = pkin(9) * t1025;
t947 = t807 * t854 - t841;
t944 = t856 * t917;
t943 = t728 * t854;
t941 = qJD(3) * t917;
t939 = t852 * t664 - t855 * t667 + t713 * t993 + t725 * t995;
t772 = t850 * t904;
t774 = t907 * qJD(2);
t938 = -t761 * t993 - t762 * t995 + t855 * t772 - t852 * t774;
t927 = t850 * t858 * t1053;
t924 = -g(1) * t652 + g(2) * t656;
t923 = g(1) * t653 - g(2) * t657;
t922 = g(1) * t948 + g(2) * t706;
t716 = -qJ(5) * t855 + t1000;
t921 = t851 * t993 - t723;
t920 = t854 * t993 - t724;
t570 = -t889 * pkin(3) + t939;
t571 = -pkin(4) * t728 + t985;
t916 = t571 * t854 - t572 * t851;
t913 = 0.2e1 * t954 + qJD(2);
t899 = -t761 * t995 + t762 * t993 + t852 * t772 + t855 * t774;
t603 = pkin(10) * t964 + t899;
t698 = qJD(3) * t786 + t850 * t963;
t775 = t1001 * qJD(2);
t609 = t698 * pkin(3) - t699 * pkin(10) + t775;
t908 = -t851 * t603 + t609 * t854 - t642 * t991 - t644 * t990;
t549 = t940 - t1069;
t902 = t620 * t728 - t1056;
t596 = -t951 + t643;
t901 = t957 + t982;
t550 = t854 * t569 + t851 * t575 + t617 * t990 - t621 * t991;
t898 = t854 * t603 + t851 * t609 + t642 * t990 - t644 * t991;
t671 = -t1025 * t787 - t788 * t854;
t673 = -t1025 * t789 - t790 * t854;
t740 = -t1027 * t854 + t855 * t975;
t896 = g(1) * t673 + g(2) * t671 + g(3) * t740;
t672 = -t1023 * t787 + t788 * t851;
t674 = -t1023 * t789 + t790 * t851;
t895 = -g(1) * t674 - g(2) * t672 - g(3) * t741;
t894 = g(1) * t706 - g(2) * t948 + g(3) * t785;
t893 = g(1) * t707 + g(2) * t703 + g(3) * t786;
t890 = t910 * MDP(8);
t887 = t741 * pkin(4) + pkin(10) * t976 + qJ(5) * t740 + t855 * t836 - t1002;
t886 = g(1) * t789 + g(2) * t787 - g(3) * t1026;
t552 = t594 * pkin(4) + t593 * qJ(5) - t681 * qJD(5) + t570;
t547 = -pkin(5) * t594 + qJDD(6) - t552;
t885 = t547 + t894;
t884 = pkin(3) * t964 + t938;
t548 = t550 + t1082;
t556 = t698 * qJ(5) + t785 * qJD(5) + t898;
t881 = t672 * pkin(4) + pkin(9) * t788 - pkin(10) * t1031 + qJ(5) * t671 + t1084 * t787;
t880 = t674 * pkin(4) + pkin(9) * t790 - pkin(10) * t1029 + qJ(5) * t673 + t1084 * t789;
t878 = t593 - t1039;
t876 = t1061 * pkin(1) + t790 * pkin(2) + t707 * pkin(3) + t657 * pkin(4) + pkin(8) * t966 + pkin(9) * t789 + qJ(5) * t656;
t875 = -t1054 * t728 + t894;
t873 = -t552 + t875;
t872 = -pkin(1) * t1060 - t788 * pkin(2) - pkin(3) * t703 - pkin(4) * t653 + pkin(8) * t967 - t787 * pkin(9) - qJ(5) * t652;
t870 = g(1) * t657 + g(2) * t653 + g(3) * t701 - t550;
t615 = t1076 * t854 - t786 * t991 + t851 * t964;
t869 = qJ(5) * t615 + qJD(5) * t701 + t884;
t866 = t578 * t681 - t871;
t865 = t579 * t728 + t870;
t845 = t855 * pkin(4);
t814 = t1055 * t851;
t802 = -pkin(3) - t919;
t794 = pkin(3) - t1064;
t763 = t909 * t852;
t717 = t845 - t947;
t714 = t897 * t852;
t688 = t1046 * t851 + t716;
t683 = pkin(5) * t855 + t841 + t845 + (-t807 - t1046) * t854;
t614 = t1076 * t851 + t786 * t990 - t854 * t964;
t608 = pkin(4) * t681 + t1051;
t595 = -t1062 * t681 - t1051;
t588 = -pkin(4) * t785 - t950;
t586 = -pkin(4) * t739 - t666 * t854 + t624;
t582 = -t596 - t1057;
t576 = qJ(6) * t700 + t587;
t568 = -qJ(6) * t701 - t1062 * t785 - t950;
t559 = -t1062 * t728 + t986;
t558 = pkin(4) * t614 - t869;
t557 = -pkin(4) * t698 - t908;
t555 = -t1062 * t614 + t869;
t554 = qJ(6) * t614 + qJD(6) * t700 + t556;
t553 = -qJ(6) * t615 - qJD(6) * t701 - t1062 * t698 - t908;
t546 = t548 + t1079;
t545 = -qJD(6) * t681 + t1048 - t1058 + t549;
t1 = [(t548 * t587 + t572 * t556 + t552 * t596 + t578 * t558 + t549 * t588 + t571 * t557 - g(1) * (pkin(10) * t948 + t872) - g(2) * (pkin(10) * t706 + t876)) * MDP(28) + (-t1016 * t631 - t550 * t785 + t570 * t701 - t580 * t698 - t643 * t593 + t620 * t615 - t681 * t884 - t728 * t898 + t924) * MDP(24) + (t570 * t700 + t579 * t698 + t643 * t594 + t620 * t614 + t631 * t950 - t679 * t884 + t728 * t908 - t785 * t940 + t923) * MDP(23) + (-t938 * t917 + t1012 * t978 + t775 * t737 + t760 * t637 + t665 * t785 + t712 * t698 + g(1) * t703 - g(2) * t707 + (t939 * t856 + (qJD(1) * t1012 + t635) * t997) * t850) * MDP(16) + (-t548 * t700 + t549 * t701 - t556 * t679 + t557 * t681 + t571 * t615 - t572 * t614 - t587 * t594 - t588 * t593 - t922) * MDP(26) + (-t545 * t701 + t546 * t700 - t553 * t681 + t554 * t679 - t559 * t615 + t560 * t614 + t568 * t593 + t576 * t594 + t922) * MDP(31) + (-t549 * t785 + t552 * t700 - t557 * t728 + t558 * t679 - t571 * t698 + t578 * t614 - t588 * t631 + t594 * t596 + t923) * MDP(25) + (-t545 * t785 - t547 * t700 - t553 * t728 - t555 * t679 - t559 * t698 - t561 * t614 - t568 * t631 - t582 * t594 + t923) * MDP(29) + (t548 * t785 - t552 * t701 + t556 * t728 - t558 * t681 + t572 * t698 - t578 * t615 + t587 * t631 + t593 * t596 - t924) * MDP(27) + (t546 * t785 + t547 * t701 + t554 * t728 + t555 * t681 + t560 * t698 + t561 * t615 + t576 * t631 - t582 * t593 - t924) * MDP(30) + (-t775 * t914 + t907 * t910 - t937 * t1053 + g(1) * t788 - g(2) * t790 + (-t958 + t980) * t959) * MDP(9) + (-g(1) * t787 + g(2) * t789 - t1001 * t910 - t1053 * t882 - t774 * t914 - t901 * t959) * MDP(10) + t1053 * t890 + (t546 * t576 + t560 * t554 + t545 * t568 + t559 * t553 + t547 * t582 + t561 * t555 - g(1) * (-pkin(5) * t653 + t1055 * t948 + t872) - g(2) * (pkin(5) * t657 + t1055 * t706 + t876)) * MDP(32) + (-t699 * t917 + t786 * t978 + ((-t891 + (-t892 - t931) * t855) * t856 + (-(-qJD(1) * t995 + t981) * t1026 + (qJD(1) * t786 + t739) * qJD(2)) * t853) * t850) * MDP(13) + (-t1011 * t889 - t1026 * t900 - t636 * t964 + t665 * t786 + t712 * t699 + t775 * t739 + t760 * t861 + t899 * t917 + t922) * MDP(17) + (g(1) * t1060 - g(2) * t1061) * MDP(2) + (g(1) * t1061 + g(2) * t1060) * MDP(3) + (-t978 * t856 + (-t832 - t917) * t997) * t850 * MDP(15) + (0.2e1 * (t853 * t980 - t983 * t999) * MDP(5) + (qJDD(1) * t848 + 0.2e1 * t853 * t957) * MDP(4)) * t847 + (t698 * t917 - t785 * t978 + (t637 * t856 + (-qJD(1) * t785 - t737) * t997) * t850) * MDP(14) + (t1073 * t856 - t913 * t964) * MDP(7) + (t850 * t913 * t996 + t1073 * t853) * MDP(6) + (t739 * t699 + t786 * t861) * MDP(11) + (-t786 * t637 - t739 * t698 - t699 * t737 - t785 * t861) * MDP(12) + (-t594 * t785 - t614 * t728 - t631 * t700 - t679 * t698) * MDP(21) + (-t593 * t785 + t615 * t728 + t631 * t701 + t681 * t698) * MDP(20) + (t631 * t785 + t698 * t728) * MDP(22) + (t593 * t700 - t594 * t701 - t614 * t681 - t615 * t679) * MDP(19) + (-t593 * t701 + t615 * t681) * MDP(18) + qJDD(1) * MDP(1); (-pkin(2) * t637 + t1009 * t917 - t635 * t965 - t773 * t737 + (-pkin(9) * t889 - t712 * t917) * t852 + (pkin(9) * t941 - t665 + t886) * t855) * MDP(16) + (-t855 * t941 + t852 * t978 + (t855 * t944 + (qJD(2) * t852 - t739) * t853) * t998) * MDP(13) + (t852 * t941 + t855 * t978 + (-t852 * t944 + (qJD(2) * t855 + t737) * t853) * t998) * MDP(14) + (t546 * t688 + t545 * t683 + t547 * t714 - g(1) * (pkin(5) * t674 + qJ(6) * t1029 + t880) - g(2) * (pkin(5) * t672 + qJ(6) * t1031 + t881) - g(3) * (pkin(5) * t741 - qJ(6) * t976 + t887) + t1019 * t561 + t1020 * t560 + t1021 * t559) * MDP(32) + (-t728 * t852 * t917 - t631 * t855) * MDP(22) + (t853 * t927 + t831) * MDP(7) + (t679 * t724 + t681 * t723 + (-t679 * t854 - t681 * t851) * t993 + (t1042 - t594 * t854 + (t679 * t851 - t681 * t854) * qJD(4)) * t852) * MDP(19) + (-t546 * t855 - t593 * t714 + t631 * t688 + t1020 * t728 + t1019 * t681 + t920 * t561 + (t547 * t854 - t560 * t917 - t561 * t991) * t852 - t896) * MDP(30) + (t559 * t724 - t560 * t723 + t593 * t683 + t594 * t688 - t1021 * t681 + t1020 * t679 + (-t559 * t854 + t560 * t851) * t993 + (-t545 * t854 + t546 * t851 + (t559 * t851 + t560 * t854) * qJD(4)) * t852 + t970) * MDP(31) + (t545 * t855 - t594 * t714 - t631 * t683 - t1021 * t728 - t1019 * t679 - t921 * t561 + (-t547 * t851 + t559 * t917 - t561 * t990) * t852 + t895) * MDP(29) - t853 * MDP(4) * t977 + (pkin(1) * t977 + t770 * t914 + g(1) * t790 + g(2) * t788 + (pkin(8) * t983 + g(3)) * t1027 - t971) * MDP(10) + (pkin(1) * t1028 * t853 + t773 * t914 + t886 - t937) * MDP(9) + (-t593 * t1024 + (-t852 * t991 + t920) * t681) * MDP(18) + t999 * MDP(5) * t1028 + (-t856 * t927 + t956) * MDP(6) + ((-qJD(3) * t965 + t910) * t852 ^ 2 + ((t850 * t901 + t892) * t852 - t917 * t739) * t855) * MDP(11) + (t1072 * t737 + t1086 * t739 - t852 * t637 + t855 * t861) * MDP(12) + t917 * MDP(15) * t965 + (t593 * t855 + t920 * t728 + (-t681 * t917 + t1040 - t961) * t852) * MDP(20) + (t594 * t855 - t921 * t728 + (t679 * t917 - t728 * t990 - t1041) * t852) * MDP(21) + t890 + (t549 * t855 + t594 * t763 - t631 * t717 - t1014 * t728 + t1013 * t679 + t921 * t578 + (t552 * t851 + t571 * t917 + t578 * t990) * t852 + t895) * MDP(25) + (-t548 * t855 + t593 * t763 + t631 * t716 + t1017 * t728 - t1013 * t681 - t920 * t578 + (-t552 * t854 - t572 * t917 + t578 * t991) * t852 - t896) * MDP(27) + (-g(1) * t880 - g(2) * t881 - g(3) * t887 + t1013 * t578 + t1014 * t571 + t1017 * t572 + t548 * t716 + t549 * t717 + t552 * t763) * MDP(28) + (-t571 * t724 + t572 * t723 - t593 * t717 - t594 * t716 + t1014 * t681 - t1017 * t679 + t916 * t993 + (-t548 * t851 + t549 * t854 + (-t571 * t851 - t572 * t854) * qJD(4)) * t852 - t970) * MDP(26) + (t947 * t631 + t903 * t679 - t620 * t723 + (-t1087 * t854 + t1081) * t728 + (t940 + (pkin(9) * t679 + t620 * t851) * qJD(3)) * t855 + (t620 * t990 + t570 * t851 - t917 * t579 + (qJD(3) * t1034 + t594) * pkin(9)) * t852 + t895) * MDP(23) + (-pkin(2) * t861 - t1008 * t917 + t636 * t965 + t665 * t852 - t773 * t739 + t970 - t1072 * t712 + (-t855 * t889 - t917 * t995) * pkin(9)) * MDP(17) + (-t1000 * t631 + t903 * t681 - t620 * t724 + t1080 * t728 + (t620 * t994 + t550 + (qJD(3) * t681 + t961) * pkin(9)) * t855 + (-t620 * t991 + t570 * t854 + t917 * t580 + (t728 * t994 - t593) * pkin(9)) * t852 + t896) * MDP(24); -t737 ^ 2 * MDP(12) + (-t737 * t917 + t861) * MDP(13) - t637 * MDP(14) + t889 * MDP(15) + (-t636 * t917 + t894 - t939) * MDP(16) + (-t635 * t917 + t712 * t737 + t893 + t900) * MDP(17) + (t681 * t943 - t1042) * MDP(18) + (t1070 * t854 + t1085 * t851) * MDP(19) + (t728 * t943 + t1041) * MDP(20) + (-t1034 * t728 + t1040) * MDP(21) + (-pkin(3) * t594 + t624 * t728 - t636 * t679 + t902 * t851 + (-t570 + (-t666 - t1054) * t728 + t894) * t854) * MDP(23) + (pkin(3) * t593 + t1018 * t728 - t636 * t681 + t902 * t854 + (t570 - t875) * t851) * MDP(24) + (t1006 * t679 + t1066 * t851 + t586 * t728 + t594 * t802 + t873 * t854) * MDP(25) + (t584 * t679 - t586 * t681 + (t548 + t728 * t571 + (-t594 + t992) * pkin(10)) * t854 + (t549 - t1044 + (qJD(4) * t679 - t593) * pkin(10)) * t851 - t893) * MDP(26) + (-t1006 * t681 - t1066 * t854 - t584 * t728 + t593 * t802 + t873 * t851) * MDP(27) + (t552 * t802 - t572 * t584 - t571 * t586 - g(1) * t973 - g(2) * t974 - g(3) * t972 + t1006 * t578 + (qJD(4) * t916 + t548 * t854 + t549 * t851 - t893) * pkin(10)) * MDP(28) + (-t1005 * t728 - t1010 * t679 - t1034 * t561 - t594 * t794 - t631 * t814 + t854 * t885) * MDP(29) + (t1007 * t728 + t1010 * t681 + t561 * t943 - t593 * t794 + t631 * t815 + t851 * t885) * MDP(30) + (t593 * t814 + t594 * t815 - t1005 * t681 + t1007 * t679 + (-t559 * t728 - t546) * t854 + (-t545 + t1045) * t851 + t893) * MDP(31) + (t546 * t815 + t545 * t814 + t547 * t794 - g(1) * (-pkin(5) * t1035 + t1055 * t707 + t973) - g(2) * (pkin(5) * t1036 + t1055 * t703 + t974) - g(3) * (-pkin(5) * t1033 + t1055 * t786 + t972) + t1010 * t561 + t1007 * t560 + t1005 * t559) * MDP(32) + (MDP(11) * t737 + MDP(12) * t739 - MDP(14) * t917 - MDP(16) * t712 - MDP(20) * t681 + MDP(21) * t679 - MDP(22) * t728 - MDP(23) * t579 + MDP(24) * t580 + MDP(25) * t571 - MDP(27) * t572 + MDP(29) * t559 - MDP(30) * t560) * t739; MDP(18) * t1038 + (t678 - t1063) * MDP(19) - t878 * MDP(20) + (-t594 + t1037) * MDP(21) + t631 * MDP(22) + (-t620 * t681 + t1043 + t874) * MDP(23) + (t620 * t679 + t865) * MDP(24) + (-t608 * t679 + t1043 + t628 - t866) * MDP(25) + (pkin(4) * t593 - t1052 + (t572 - t580) * t681 + (t571 - t985) * t679) * MDP(26) + (-t578 * t679 + t608 * t681 + t1078 - t865) * MDP(27) + (-t549 * pkin(4) - g(1) * t952 - g(2) * t953 - g(3) * t951 + t548 * qJ(5) - t571 * t580 + t572 * t985 - t578 * t608) * MDP(28) + ((pkin(5) + t1062) * t631 + t563 * t728 + t595 * t679 + t1068) * MDP(29) + (t561 * t679 - t562 * t728 - t595 * t681 + t1078 + t1079 - t870) * MDP(30) + (t1052 - t593 * t1062 + (-t560 + t563) * t681 + (-t559 + t986) * t679) * MDP(31) + (t546 * qJ(5) - t545 * t1062 - t559 * t563 - t561 * t595 - g(1) * (-pkin(5) * t656 + t952) - g(2) * (-pkin(5) * t652 + t953) - g(3) * (t951 - t1057) + t986 * t560) * MDP(32); (t866 - t1044) * MDP(28) + (-t1045 - t1058 - t1068) * MDP(32) + (MDP(25) + MDP(29)) * (t1038 - t631) + (MDP(27) + MDP(30)) * (-t728 ^ 2 - t678) + (-MDP(26) + MDP(31)) * t878; t1085 * MDP(29) + t1070 * MDP(30) + (-t678 - t1063) * MDP(31) + (t559 * t681 - t560 * t679 + t885) * MDP(32);];
tau  = t1;
