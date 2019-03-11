% Calculate vector of inverse dynamics joint torques for
% S6RRRRRR3
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
%   see S6RRRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:43:09
% EndTime: 2019-03-10 03:43:35
% DurationCPUTime: 18.54s
% Computational Cost: add. (16424->718), mult. (35686->917), div. (0->0), fcn. (27894->18), ass. (0->321)
t902 = sin(qJ(6));
t905 = sin(qJ(3));
t911 = cos(qJ(2));
t1033 = t905 * t911;
t906 = sin(qJ(2));
t1064 = cos(qJ(3));
t987 = qJD(1) * t1064;
t811 = -qJD(1) * t1033 - t906 * t987;
t897 = qJD(2) + qJD(3);
t904 = sin(qJ(4));
t910 = cos(qJ(4));
t768 = -t811 * t904 - t910 * t897;
t903 = sin(qJ(5));
t909 = cos(qJ(5));
t946 = t811 * t910 - t897 * t904;
t947 = t768 * t903 + t909 * t946;
t1050 = t947 * t902;
t694 = t909 * t768 - t903 * t946;
t908 = cos(qJ(6));
t627 = -t908 * t694 + t1050;
t833 = t1064 * t906 + t1033;
t758 = t897 * t833;
t980 = qJDD(1) * t1064;
t999 = qJDD(1) * t906;
t952 = t905 * t999 - t911 * t980;
t708 = qJD(1) * t758 + t952;
t706 = qJDD(4) + t708;
t704 = qJDD(5) + t706;
t702 = qJDD(6) + t704;
t949 = t694 * t902 + t908 * t947;
t1109 = t702 * MDP(36) + (-t627 ^ 2 + t949 ^ 2) * MDP(33) + t627 * MDP(32) * t949;
t1066 = pkin(7) + pkin(8);
t858 = t1066 * t911;
t840 = qJD(1) * t858;
t815 = t905 * t840;
t856 = t1066 * t906;
t838 = qJD(1) * t856;
t753 = -t1064 * t838 - t815;
t986 = qJD(3) * t1064;
t1105 = -pkin(2) * t986 + t753;
t1067 = qJD(4) + qJD(5);
t1008 = qJD(1) * t906;
t1095 = -t905 * t1008 + t911 * t987;
t1108 = t1067 - t1095;
t741 = -pkin(3) * t811 - pkin(9) * t1095;
t727 = pkin(2) * t1008 + t741;
t1107 = t1105 * t904 - t910 * t727;
t803 = qJD(4) - t1095;
t791 = qJD(5) + t803;
t787 = qJD(6) + t791;
t1091 = t787 * t949;
t1092 = t627 * t787;
t1002 = qJD(6) * t902;
t1001 = qJD(6) * t908;
t1003 = qJD(5) * t909;
t1004 = qJD(5) * t903;
t1005 = qJD(4) * t910;
t1006 = qJD(4) * t904;
t1104 = t1095 * t897;
t998 = qJDD(1) * t911;
t707 = t905 * t998 + t906 * t980 + t1104;
t896 = qJDD(2) + qJDD(3);
t672 = t897 * t1005 + t1006 * t811 + t910 * t707 + t904 * t896;
t673 = -qJD(4) * t946 + t707 * t904 - t910 * t896;
t603 = -t768 * t1003 + t1004 * t946 + t909 * t672 - t903 * t673;
t923 = qJD(5) * t947 - t672 * t903 - t909 * t673;
t994 = -t694 * t1001 + t908 * t603 + t902 * t923;
t564 = t1002 * t947 + t994;
t976 = t603 * t902 - t908 * t923;
t924 = qJD(6) * t949 - t976;
t1106 = t704 * MDP(29) + (-t694 ^ 2 + t947 ^ 2) * MDP(26) + (t694 * t791 + t603) * MDP(27) + (-t791 * t947 + t923) * MDP(28) - t694 * MDP(25) * t947 + (t924 - t1091) * MDP(35) + (t564 - t1092) * MDP(34) + t1109;
t901 = qJ(2) + qJ(3);
t893 = cos(t901);
t1058 = g(3) * t893;
t891 = sin(t901);
t907 = sin(qJ(1));
t912 = cos(qJ(1));
t960 = g(1) * t912 + g(2) * t907;
t1101 = t960 * t891;
t1084 = t1058 - t1101;
t1094 = pkin(11) * t694;
t887 = -pkin(2) * t911 - pkin(1);
t854 = t887 * qJD(1);
t724 = -pkin(3) * t1095 + pkin(9) * t811 + t854;
t818 = t1064 * t840;
t1053 = qJD(2) * pkin(2);
t819 = -t838 + t1053;
t746 = t905 * t819 + t818;
t731 = pkin(9) * t897 + t746;
t675 = t910 * t724 - t731 * t904;
t640 = pkin(10) * t946 + t675;
t620 = pkin(4) * t803 + t640;
t676 = t724 * t904 + t731 * t910;
t641 = -pkin(10) * t768 + t676;
t637 = t909 * t641;
t594 = t620 * t903 + t637;
t583 = t594 - t1094;
t581 = t583 * t1002;
t745 = t1064 * t819 - t815;
t730 = -t897 * pkin(3) - t745;
t687 = t768 * pkin(4) + t730;
t621 = t694 * pkin(5) + t687;
t1041 = t893 * t907;
t900 = qJ(4) + qJ(5);
t894 = qJ(6) + t900;
t881 = sin(t894);
t882 = cos(t894);
t777 = -t1041 * t882 + t881 * t912;
t1040 = t893 * t912;
t779 = t1040 * t882 + t881 * t907;
t877 = g(3) * t891;
t1082 = g(1) * t779 - g(2) * t777 - t621 * t627 + t882 * t877 + t581;
t1037 = t903 * t910;
t832 = t904 * t909 + t1037;
t1086 = t1108 * t832;
t830 = t903 * t904 - t909 * t910;
t1019 = t1108 * t830;
t1100 = t1105 * t910 + t904 * t727;
t776 = t1041 * t881 + t882 * t912;
t778 = -t1040 * t881 + t882 * t907;
t1000 = qJD(1) * qJD(2);
t985 = t906 * t1000;
t805 = pkin(2) * t985 + qJDD(1) * t887;
t634 = pkin(3) * t708 - pkin(9) * t707 + t805;
t632 = t910 * t634;
t1007 = qJD(3) * t905;
t984 = t911 * t1000;
t762 = qJDD(2) * pkin(2) + t1066 * (-t984 - t999);
t767 = t1066 * (-t985 + t998);
t922 = -t840 * t1007 + t1064 * t767 + t905 * t762 + t819 * t986;
t649 = t896 * pkin(9) + t922;
t572 = pkin(4) * t706 - pkin(10) * t672 - qJD(4) * t676 - t649 * t904 + t632;
t935 = t724 * t1005 - t1006 * t731 + t904 * t634 + t910 * t649;
t576 = -pkin(10) * t673 + t935;
t978 = t909 * t572 - t903 * t576;
t925 = -qJD(5) * t594 + t978;
t556 = pkin(5) * t704 - pkin(11) * t603 + t925;
t965 = -t620 * t1003 + t641 * t1004 - t903 * t572 - t909 * t576;
t557 = pkin(11) * t923 - t965;
t979 = t908 * t556 - t902 * t557;
t1081 = -g(1) * t778 + g(2) * t776 + t621 * t949 + t881 * t877 + t979;
t895 = t910 * pkin(10);
t963 = -t811 * pkin(4) - t1095 * t895;
t883 = pkin(2) * t905 + pkin(9);
t1054 = -pkin(10) - t883;
t981 = qJD(4) * t1054;
t1099 = -t910 * t981 - t1107 + t963;
t735 = t910 * t741;
t1065 = -pkin(9) - pkin(10);
t992 = qJD(4) * t1065;
t1098 = -t745 * t904 - t910 * t992 + t735 + t963;
t1044 = t1095 * t904;
t997 = pkin(10) * t1044;
t1097 = -t904 * t981 + t1100 - t997;
t1017 = t904 * t741 + t910 * t745;
t1096 = -t904 * t992 + t1017 - t997;
t1093 = pkin(11) * t947;
t749 = t908 * t830 + t832 * t902;
t1027 = -qJD(6) * t749 - t1019 * t908 - t1086 * t902;
t750 = -t830 * t902 + t832 * t908;
t1089 = qJD(6) * t750 - t1019 * t902 + t1086 * t908;
t752 = -t905 * t838 + t818;
t962 = pkin(2) * t1007 - t752;
t964 = t819 * t1007 - t1064 * t762 + t905 * t767 + t840 * t986;
t650 = -pkin(3) * t896 + t964;
t982 = -t650 - t1058;
t941 = t1064 * t911 - t905 * t906;
t757 = t897 * t941;
t1048 = t757 * t904;
t940 = -t833 * t1005 - t1048;
t890 = sin(t900);
t892 = cos(t900);
t781 = -t1041 * t892 + t890 * t912;
t783 = t1040 * t892 + t890 * t907;
t1080 = g(1) * t783 - g(2) * t781 + t687 * t694 + t892 * t877 + t965;
t780 = t1041 * t890 + t892 * t912;
t782 = -t1040 * t890 + t892 * t907;
t1079 = -g(1) * t782 + g(2) * t780 + t687 * t947 + t890 * t877 + t925;
t736 = t832 * t833;
t1042 = t833 * t910;
t744 = -pkin(3) * t941 - pkin(9) * t833 + t887;
t739 = t910 * t744;
t772 = t1064 * t858 - t905 * t856;
t660 = -pkin(4) * t941 - pkin(10) * t1042 - t772 * t904 + t739;
t763 = t910 * t772;
t1015 = t904 * t744 + t763;
t1043 = t833 * t904;
t678 = -pkin(10) * t1043 + t1015;
t1023 = t903 * t660 + t909 * t678;
t1076 = t1086 * pkin(11);
t1075 = t1099 * t909;
t823 = t1054 * t904;
t824 = t883 * t910 + t895;
t1014 = t903 * t823 + t909 * t824;
t1074 = t1098 * t909;
t855 = t1065 * t904;
t857 = pkin(9) * t910 + t895;
t1013 = t903 * t855 + t909 * t857;
t788 = pkin(4) * t1044;
t1073 = -t788 + t962;
t1072 = -t1064 * t856 - t905 * t858;
t888 = pkin(4) * t1006;
t1071 = pkin(5) * t1086 + t888;
t1070 = -t811 * pkin(5) - pkin(11) * t1019;
t1069 = -t855 * t1003 + t1004 * t857 + t1096 * t909 + t1098 * t903;
t1068 = -t823 * t1003 + t1004 * t824 + t1097 * t909 + t1099 * t903;
t1063 = pkin(4) * t903;
t1061 = pkin(11) * t832;
t1056 = t830 * pkin(5);
t1055 = t910 * pkin(4);
t635 = t903 * t641;
t593 = t909 * t620 - t635;
t582 = t593 + t1093;
t579 = pkin(5) * t791 + t582;
t1052 = t579 * t908;
t1051 = t672 * t904;
t1049 = t730 * t1095;
t1047 = t757 * t910;
t1046 = t768 * t803;
t1045 = t946 * t803;
t1039 = t902 * t702;
t1038 = t903 * t908;
t1036 = t904 * t706;
t1035 = t904 * t907;
t1034 = t904 * t912;
t1032 = t907 * t910;
t1031 = t908 * t583;
t1030 = t908 * t702;
t1029 = t910 * t912;
t1028 = t909 * t640 - t635;
t1021 = t1071 + t1073;
t705 = t788 + t746;
t1018 = -t705 + t1071;
t1012 = t888 + t1073;
t898 = t906 ^ 2;
t1011 = -t911 ^ 2 + t898;
t996 = t906 * t1053;
t995 = qJD(4) * pkin(9) * t803;
t886 = -pkin(3) - t1055;
t993 = qJD(2) * t1066;
t989 = t833 * t1006;
t718 = t730 * t1005;
t686 = pkin(3) * t758 - pkin(9) * t757 + t996;
t683 = t910 * t686;
t839 = t906 * t993;
t841 = t911 * t993;
t698 = qJD(3) * t1072 - t1064 * t839 - t905 * t841;
t588 = -pkin(10) * t1047 + pkin(4) * t758 - t698 * t904 + t683 + (-t763 + (pkin(10) * t833 - t744) * t904) * qJD(4);
t934 = t744 * t1005 - t1006 * t772 + t904 * t686 + t910 * t698;
t597 = pkin(10) * t940 + t934;
t977 = t909 * t588 - t597 * t903;
t975 = -t640 * t903 - t637;
t974 = t909 * t660 - t678 * t903;
t971 = t909 * t823 - t824 * t903;
t970 = t909 * t855 - t857 * t903;
t969 = t803 * t910;
t968 = -qJD(4) * t724 - t649;
t967 = qJD(6) * t579 + t557;
t885 = -pkin(2) * t1064 - pkin(3);
t961 = -t705 + t888;
t959 = g(1) * t907 - g(2) * t912;
t958 = -t1005 * t731 + t632;
t822 = t830 * pkin(11);
t712 = -t822 + t1014;
t957 = qJD(5) * t1014 + qJD(6) * t712 - t1097 * t903 + t1070 + t1075;
t729 = -t822 + t1013;
t956 = qJD(5) * t1013 + qJD(6) * t729 - t1096 * t903 + t1070 + t1074;
t711 = t971 - t1061;
t955 = -qJD(6) * t711 + t1068 + t1076;
t728 = t970 - t1061;
t954 = -qJD(6) * t728 + t1069 + t1076;
t953 = -pkin(9) * t706 - t1049;
t563 = t902 * t579 + t1031;
t948 = -t706 * t883 - t1049;
t737 = t830 * t833;
t679 = t908 * t736 - t737 * t902;
t680 = -t736 * t902 - t737 * t908;
t945 = -t676 * t811 - t904 * t982 + t718;
t944 = t730 * t1006 + t1101 * t910 + t675 * t811;
t725 = pkin(4) * t1043 - t1072;
t942 = -0.2e1 * pkin(1) * t1000 - pkin(7) * qJDD(2);
t939 = -t989 + t1047;
t853 = t885 - t1055;
t699 = -t856 * t1007 + t1064 * t841 - t905 * t839 + t858 * t986;
t933 = t660 * t1003 - t1004 * t678 + t903 * t588 + t909 * t597;
t611 = pkin(4) * t673 + t650;
t646 = -pkin(4) * t940 + t699;
t913 = qJD(2) ^ 2;
t929 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t913 + t959;
t914 = qJD(1) ^ 2;
t928 = pkin(1) * t914 - pkin(7) * qJDD(1) + t960;
t926 = t854 * t811 - t1084 - t964;
t562 = -t583 * t902 + t1052;
t574 = -pkin(5) * t923 + t611;
t921 = -t1084 * t882 + t1089 * t621 + t562 * t811 + t574 * t749;
t920 = -t1084 * t892 + t1086 * t687 + t593 * t811 + t611 * t830;
t919 = t1027 * t621 + t1084 * t881 - t563 * t811 + t574 * t750;
t918 = -t1019 * t687 + t1084 * t890 - t594 * t811 + t611 * t832;
t916 = g(1) * t1040 + g(2) * t1041 - t1095 * t854 + t877 - t922;
t915 = (t1027 * t627 + t1089 * t949 - t564 * t749 + t750 * t924) * MDP(33) + (-t1027 * t949 + t564 * t750) * MDP(32) + (t1019 * t694 + t1086 * t947 - t603 * t830 + t832 * t923) * MDP(26) + (t1027 * t787 + t702 * t750 - t811 * t949) * MDP(34) + (-t1089 * t787 + t627 * t811 - t702 * t749) * MDP(35) + (t1019 * t947 + t603 * t832) * MDP(25) + ((t672 - t1046) * t910 + (-t673 + t1045) * t904) * MDP(19) + (-t1019 * t791 + t704 * t832 - t811 * t947) * MDP(27) + (-t1086 * t791 - t694 * t811 - t704 * t830) * MDP(28) + (-t946 * t969 + t1051) * MDP(18) + (-t803 ^ 2 * t904 + t706 * t910 - t768 * t811) * MDP(21) + (t803 * t969 - t811 * t946 + t1036) * MDP(20) + (t707 - t1104) * MDP(13) + (-t952 + (-qJD(1) * t833 - t811) * t897) * MDP(14) + (-t1095 ^ 2 + t811 ^ 2) * MDP(12) + t896 * MDP(15) + (MDP(11) * t1095 + t803 * MDP(22) + t791 * MDP(29) + t787 * MDP(36)) * t811;
t884 = pkin(4) * t909 + pkin(5);
t800 = t1029 * t893 + t1035;
t799 = -t1034 * t893 + t1032;
t798 = -t1032 * t893 + t1034;
t797 = t1035 * t893 + t1029;
t784 = t886 + t1056;
t775 = t853 + t1056;
t681 = t736 * pkin(5) + t725;
t674 = -pkin(4) * t946 - pkin(5) * t947;
t617 = t757 * t1037 - t903 * t989 - t1004 * t1043 + (t1042 * t1067 + t1048) * t909;
t616 = -t1067 * t736 - t830 * t757;
t608 = -pkin(11) * t736 + t1023;
t605 = -pkin(5) * t941 + pkin(11) * t737 + t974;
t592 = t617 * pkin(5) + t646;
t585 = t1028 + t1093;
t584 = t975 + t1094;
t578 = qJD(6) * t680 + t616 * t902 + t908 * t617;
t577 = -qJD(6) * t679 + t616 * t908 - t617 * t902;
t561 = -pkin(11) * t617 + t933;
t559 = pkin(5) * t758 - pkin(11) * t616 - qJD(5) * t1023 + t977;
t1 = [(t906 * t942 + t911 * t929) * MDP(9) + (-t906 * t929 + t911 * t942) * MDP(10) + (qJDD(1) * t898 + 0.2e1 * t906 * t984) * MDP(4) + (t564 * t680 - t577 * t949) * MDP(32) + (-t564 * t679 + t577 * t627 + t578 * t949 + t680 * t924) * MDP(33) + (-t578 * t787 + t627 * t758 - t679 * t702 - t924 * t941) * MDP(35) + ((t559 * t908 - t561 * t902) * t787 + (t605 * t908 - t608 * t902) * t702 - t979 * t941 + t562 * t758 - t592 * t627 - t681 * t924 + t574 * t679 + t621 * t578 - g(1) * t777 - g(2) * t779 + ((-t605 * t902 - t608 * t908) * t787 + t563 * t941) * qJD(6)) * MDP(37) + (t1095 * t757 + t707 * t941 - t708 * t833 + t758 * t811) * MDP(12) + (t1072 * t896 - t1095 * t996 - t699 * t897 + t708 * t887 + t758 * t854 - t805 * t941 + t893 * t959) * MDP(16) + 0.2e1 * (-t1000 * t1011 + t906 * t998) * MDP(5) + (-t603 * t737 - t616 * t947) * MDP(25) + ((-t768 * t910 + t904 * t946) * t757 + (-t1051 - t673 * t910 + (t768 * t904 + t910 * t946) * qJD(4)) * t833) * MDP(19) + (t1042 * t672 - t939 * t946) * MDP(18) + (-t1036 * t833 + t673 * t941 - t758 * t768 + t803 * t940) * MDP(21) + (-t758 * t897 + t896 * t941) * MDP(14) + (-t706 * t941 + t758 * t803) * MDP(22) + (-t704 * t941 + t758 * t791) * MDP(29) + (-t702 * t941 + t758 * t787) * MDP(36) + (-g(1) * t776 - g(2) * t778 - t563 * t758 + t681 * t564 + t574 * t680 + t621 * t577 - t581 * t941 - t592 * t949 + (-(-qJD(6) * t608 + t559) * t787 - t605 * t702 + t556 * t941) * t902 + (-(qJD(6) * t605 + t561) * t787 - t608 * t702 + t967 * t941) * t908) * MDP(38) + (-t564 * t941 + t577 * t787 + t680 * t702 - t758 * t949) * MDP(34) + (-g(1) * t780 - g(2) * t782 - t1023 * t704 - t594 * t758 + t725 * t603 - t611 * t737 + t687 * t616 - t646 * t947 - t791 * t933 - t941 * t965) * MDP(31) + (-t603 * t941 + t616 * t791 - t704 * t737 - t758 * t947) * MDP(27) + (t1042 * t706 - t672 * t941 - t758 * t946 + t803 * t939) * MDP(20) + (-t698 * t897 + t707 * t887 + t757 * t854 - t772 * t896 + t805 * t833 - t811 * t996 - t891 * t959) * MDP(17) + t959 * MDP(2) + t960 * MDP(3) + (-t603 * t736 - t616 * t694 + t617 * t947 - t737 * t923) * MDP(26) + (-t617 * t791 - t694 * t758 - t704 * t736 - t923 * t941) * MDP(28) + ((-t772 * t1005 + t683) * t803 + t739 * t706 - t958 * t941 + t675 * t758 + t699 * t768 - t1072 * t673 + t833 * t718 - g(1) * t798 - g(2) * t800 + ((-qJD(4) * t744 - t698) * t803 - t772 * t706 - t968 * t941 + t650 * t833 + t730 * t757) * t904) * MDP(23) + (-g(1) * t797 - g(2) * t799 - t1015 * t706 + t1042 * t650 - t1072 * t672 - t676 * t758 - t699 * t946 + t730 * t939 - t803 * t934 + t935 * t941) * MDP(24) + (t977 * t791 + t974 * t704 - t978 * t941 + t593 * t758 + t646 * t694 - t725 * t923 + t611 * t736 + t687 * t617 - g(1) * t781 - g(2) * t783 + (-t1023 * t791 + t594 * t941) * qJD(5)) * MDP(30) + (qJDD(2) * t906 + t911 * t913) * MDP(6) + (qJDD(2) * t911 - t906 * t913) * MDP(7) + (t757 * t897 + t833 * t896) * MDP(13) + (t707 * t833 - t757 * t811) * MDP(11) + qJDD(1) * MDP(1); (-g(3) * t911 + t906 * t928) * MDP(9) + (g(3) * t906 + t911 * t928) * MDP(10) + (t885 * t673 + t982 * t910 + t948 * t904 + t962 * t768 + (-t883 * t1005 + t1107) * t803 + t944) * MDP(23) + ((t711 * t908 - t712 * t902) * t702 - t775 * t924 + (t902 * t955 - t908 * t957) * t787 - t1021 * t627 + t921) * MDP(37) + (-(t711 * t902 + t712 * t908) * t702 + t775 * t564 + (t902 * t957 + t908 * t955) * t787 - t1021 * t949 + t919) * MDP(38) + (-t1012 * t947 - t1014 * t704 + t1068 * t791 + t853 * t603 + t918) * MDP(31) + t915 + (t885 * t672 + t948 * t910 - t904 * t1101 - t962 * t946 + (t1006 * t883 + t1100) * t803 + t945) * MDP(24) + MDP(6) * t999 + (t752 * t897 + (-t1007 * t897 + t1008 * t1095 + t1064 * t896) * pkin(2) + t926) * MDP(16) + MDP(7) * t998 + (t753 * t897 + (t1008 * t811 - t896 * t905 - t897 * t986) * pkin(2) + t916) * MDP(17) + qJDD(2) * MDP(8) + (t971 * t704 - t853 * t923 + (-t824 * t1003 + (-qJD(5) * t823 + t1097) * t903 - t1075) * t791 + t1012 * t694 + t920) * MDP(30) + (-MDP(4) * t906 * t911 + MDP(5) * t1011) * t914; ((t728 * t908 - t729 * t902) * t702 - t784 * t924 + (t902 * t954 - t908 * t956) * t787 - t1018 * t627 + t921) * MDP(37) + (t970 * t704 - t886 * t923 + (-t857 * t1003 + (-qJD(5) * t855 + t1096) * t903 - t1074) * t791 + t961 * t694 + t920) * MDP(30) + t915 + (-(t728 * t902 + t729 * t908) * t702 + t784 * t564 + (t902 * t956 + t908 * t954) * t787 - t1018 * t949 + t919) * MDP(38) + (-pkin(3) * t672 + t1017 * t803 + t746 * t946 + t953 * t910 + (-t1101 + t995) * t904 + t945) * MDP(24) + (-t1013 * t704 + t1069 * t791 + t886 * t603 - t947 * t961 + t918) * MDP(31) + (-pkin(3) * t673 - t735 * t803 - t746 * t768 + (t745 * t803 + t953) * t904 + (t982 - t995) * t910 + t944) * MDP(23) + (t745 * t897 + t916) * MDP(17) + (t746 * t897 + t926) * MDP(16); (t672 + t1046) * MDP(20) + (g(1) * t800 - g(2) * t798 + t675 * t803 + t730 * t768 + t877 * t910 - t935) * MDP(24) + (t884 * t1030 - (t584 * t908 - t585 * t902) * t787 + t674 * t627 + (-t903 * t1039 + (-t902 * t909 - t1038) * t787 * qJD(5)) * pkin(4) + ((-pkin(4) * t1038 - t884 * t902) * t787 - t563) * qJD(6) + t1081) * MDP(37) + (t674 * t949 + (-t884 * t702 - t556 + (t584 - (-qJD(5) - qJD(6)) * t1063) * t787) * t902 + (-t702 * t1063 + (-pkin(4) * t1003 - qJD(6) * t884 + t585) * t787 - t967) * t908 + t1082) * MDP(38) - t946 * t768 * MDP(18) + (-t673 - t1045) * MDP(21) + (-t768 ^ 2 + t946 ^ 2) * MDP(19) + (-g(1) * t799 + g(2) * t797 + t676 * t803 + t730 * t946 + (t968 + t877) * t904 + t958) * MDP(23) + (-t975 * t791 + (-t1004 * t791 + t694 * t946 + t704 * t909) * pkin(4) + t1079) * MDP(30) + (t1028 * t791 + (-t1003 * t791 - t704 * t903 - t946 * t947) * pkin(4) + t1080) * MDP(31) + t706 * MDP(22) + t1106; (t594 * t791 + t1079) * MDP(30) + (t593 * t791 + t1080) * MDP(31) + (-(-t582 * t902 - t1031) * t787 - t563 * qJD(6) + (-t1002 * t787 - t627 * t947 + t1030) * pkin(5) + t1081) * MDP(37) + ((-t583 * t787 - t556) * t902 + (t582 * t787 - t967) * t908 + (-t1001 * t787 - t947 * t949 - t1039) * pkin(5) + t1082) * MDP(38) + t1106; (t994 - t1092) * MDP(34) + (-t976 - t1091) * MDP(35) + (t563 * t787 + t1081) * MDP(37) + (-t902 * t556 - t908 * t557 + t562 * t787 + t1082) * MDP(38) + (MDP(34) * t1050 + MDP(35) * t949 - MDP(37) * t563 - MDP(38) * t1052) * qJD(6) + t1109;];
tau  = t1;
