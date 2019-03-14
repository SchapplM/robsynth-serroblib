% Calculate vector of inverse dynamics joint torques for
% S6PRRRRR6
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:25:54
% EndTime: 2019-03-09 01:26:37
% DurationCPUTime: 31.59s
% Computational Cost: add. (18728->910), mult. (52929->1316), div. (0->0), fcn. (47235->18), ass. (0->388)
t904 = cos(qJ(3));
t1037 = qJD(3) * t904;
t894 = cos(pkin(7));
t1007 = t894 * t1037;
t892 = sin(pkin(6));
t905 = cos(qJ(2));
t1063 = t904 * t905;
t899 = sin(qJ(3));
t900 = sin(qJ(2));
t1071 = t899 * t900;
t940 = -t1071 * t894 + t1063;
t823 = t940 * t892;
t1128 = -pkin(2) * t1007 + qJD(1) * t823;
t891 = sin(pkin(7));
t893 = cos(pkin(8));
t969 = t891 * (-pkin(11) * t893 - pkin(10));
t950 = t899 * t969;
t1127 = qJD(3) * t950 - t1128;
t903 = cos(qJ(4));
t1035 = qJD(4) * t903;
t1003 = t893 * t1035;
t890 = sin(pkin(8));
t1005 = t890 * t1035;
t898 = sin(qJ(4));
t1036 = qJD(4) * t898;
t1079 = t893 * t898;
t1087 = t890 * t898;
t1082 = t891 * t904;
t1076 = t894 * t899;
t883 = pkin(2) * t1076;
t1044 = pkin(10) * t1082 + t883;
t1077 = t893 * t904;
t783 = (t1077 * t891 + t890 * t894) * pkin(11) + t1044;
t1075 = t894 * t904;
t884 = pkin(2) * t1075;
t792 = pkin(3) * t894 + t884 + t950;
t794 = (t904 * t969 - t883) * qJD(3);
t1106 = pkin(11) * t890;
t949 = -pkin(3) * t904 - t1106 * t899;
t821 = (-pkin(2) + t949) * t891;
t929 = (pkin(3) * t899 - t1106 * t904) * t891;
t828 = qJD(3) * t929;
t1068 = t900 * t904;
t1069 = t899 * t905;
t942 = -t1068 * t894 - t1069;
t822 = t942 * t892;
t812 = qJD(1) * t822;
t1081 = t892 * t900;
t1013 = qJD(1) * t1081;
t978 = t891 * t1013;
t920 = t812 * t893 + t890 * t978;
t1126 = t792 * t1003 + t821 * t1005 - t1036 * t783 + t794 * t1079 + t828 * t1087 + t1127 * t903 - t898 * t920;
t1052 = (t828 - t978) * t893 + (-t794 + t812) * t890;
t1040 = qJD(2) * t891;
t1065 = t903 * t904;
t1073 = t898 * t899;
t944 = -t1073 * t893 + t1065;
t813 = t944 * t1040;
t1125 = -t813 + t1005;
t1078 = t893 * t903;
t879 = pkin(11) * t1087;
t1114 = pkin(3) * t1078 - t879;
t1084 = t891 * t893;
t1022 = pkin(11) * t1084;
t1041 = qJD(1) * t905;
t1104 = qJD(2) * pkin(2);
t862 = t1041 * t892 + t1104;
t895 = cos(pkin(6));
t1042 = qJD(1) * t895;
t1014 = t891 * t1042;
t867 = t904 * t1014;
t1046 = t862 * t1075 + t867;
t848 = pkin(10) * t1040 + t1013;
t735 = (-qJD(2) * t1022 - t848) * t899 + t1046;
t1010 = qJD(2) * t1077;
t1064 = t904 * t848;
t736 = -t862 * t1076 - t1064 + (-pkin(11) * t1010 - t1042 * t899) * t891;
t827 = qJD(2) * t929;
t1124 = t1114 * qJD(4) - t736 * t1079 - t827 * t1087 - t903 * t735;
t1038 = qJD(3) * t899;
t1009 = t891 * t1038;
t974 = t890 * t1009;
t1123 = -pkin(12) * t974 - t1126;
t1112 = t893 * t1065 - t1073;
t733 = t894 * t1005 + (t944 * qJD(3) + qJD(4) * t1112) * t891;
t1006 = t890 * t1036;
t1070 = t899 * t903;
t1072 = t898 * t904;
t945 = t1072 * t893 + t1070;
t946 = t1070 * t893 + t1072;
t907 = qJD(3) * t946 + qJD(4) * t945;
t734 = t1006 * t894 + t891 * t907;
t1122 = pkin(4) * t734 - pkin(12) * t733 + t1052;
t1012 = t899 * t1040;
t977 = t890 * t1012;
t1121 = pkin(12) * t977 - t1124;
t688 = -t736 * t890 + t893 * t827;
t811 = t946 * t1040;
t1119 = pkin(4) * t811 - pkin(12) * t813 + t688 - (pkin(4) * t898 - pkin(12) * t903) * t890 * qJD(4);
t1004 = t893 * t1036;
t1118 = t792 * t1004 + t821 * t1006 + t783 * t1035 + t1127 * t898 + t903 * t920;
t897 = sin(qJ(5));
t902 = cos(qJ(5));
t844 = t1087 * t897 - t902 * t893;
t1051 = -qJD(5) * t844 + t1125 * t902 - t897 * t977;
t1086 = t890 * t902;
t845 = t1086 * t898 + t893 * t897;
t1050 = qJD(5) * t845 + t1125 * t897 + t902 * t977;
t1085 = t890 * t903;
t1045 = pkin(3) * t1079 + pkin(11) * t1085;
t1117 = t1045 * qJD(4) - t898 * t735;
t973 = t891 * t1010;
t856 = t903 * t973;
t1039 = qJD(2) * t894;
t983 = qJD(3) + t1039;
t951 = t983 * t890;
t975 = t898 * t1012;
t767 = -t903 * t951 - t856 + t975;
t765 = qJD(5) + t767;
t924 = t945 * t891;
t769 = qJD(2) * t924 + t898 * t951;
t1021 = t890 * t1082;
t921 = qJD(2) * t1021 - t893 * t983 - qJD(4);
t807 = t902 * t921;
t717 = t769 * t897 + t807;
t713 = qJD(6) + t717;
t967 = -t811 + t1006;
t1025 = qJDD(2) * t899;
t1026 = qJDD(2) * t894;
t878 = qJDD(3) + t1026;
t938 = qJD(4) * t951;
t1024 = qJDD(2) * t904;
t998 = t891 * t1024;
t971 = t893 * t998;
t679 = -t878 * t1085 + (qJD(2) * t907 + t898 * t1025) * t891 + t898 * t938 - t903 * t971;
t1059 = -t794 * t1078 + (-pkin(4) * t1009 - t828 * t903) * t890 + t1118;
t1031 = qJD(5) * t902;
t1033 = qJD(5) * t897;
t729 = -t792 * t890 + t893 * t821;
t786 = -t1085 * t894 - t1112 * t891;
t789 = t1087 * t894 + t924;
t654 = pkin(4) * t786 - pkin(12) * t789 + t729;
t1015 = t792 * t1079 + t821 * t1087 + t903 * t783;
t838 = -t893 * t894 + t1021;
t667 = -pkin(12) * t838 + t1015;
t1116 = -t654 * t1031 + t1033 * t667 - t1122 * t897 + t1123 * t902;
t1115 = t897 * t654 + t902 * t667;
t831 = pkin(12) * t893 + t1045;
t832 = (-pkin(4) * t903 - pkin(12) * t898 - pkin(3)) * t890;
t1048 = t902 * t831 + t897 * t832;
t1047 = t736 * t1078 - (-pkin(4) * t1012 - t827 * t903) * t890 + t1117;
t1113 = -t832 * t1031 + t1033 * t831 + t1119 * t897 + t1121 * t902;
t939 = t862 * t894 + t1014;
t716 = t1064 + t939 * t899 + (t951 + t973) * pkin(11);
t722 = pkin(3) * t983 + t735;
t876 = t894 * t1042;
t752 = t876 + (qJD(2) * t949 - t862) * t891;
t626 = -t898 * t716 + t903 * (t722 * t893 + t752 * t890);
t776 = qJD(2) * t974 + t878 * t893 - t890 * t998 + qJDD(4);
t1027 = qJDD(1) * t895;
t1000 = t891 * t1027;
t816 = pkin(10) * qJDD(2) * t891 + (qJD(2) * t1041 + qJDD(1) * t900) * t892;
t1080 = t892 * t905;
t877 = qJDD(1) * t1080;
t1011 = qJD(2) * t1081;
t972 = qJD(1) * t1011;
t826 = qJDD(2) * pkin(2) + t877 - t972;
t961 = qJD(3) * t867 + t899 * t1000 + t862 * t1007 + t826 * t1076 + t904 * t816;
t918 = -t1038 * t848 + t961;
t1028 = qJD(2) * qJD(3);
t935 = t1028 * t899 - t1024;
t645 = (-t1084 * t935 + t878 * t890) * pkin(11) + t918;
t1049 = t904 * t1000 + t826 * t1075;
t646 = pkin(3) * t878 + (-qJDD(2) * t1022 - t816) * t899 + t736 * qJD(3) + t1049;
t875 = t894 * t1027;
t1001 = t904 * t1028;
t936 = t1001 + t1025;
t697 = t875 + (pkin(3) * t935 - t1106 * t936 - t826) * t891;
t917 = -t722 * t1003 - t752 * t1005 + t1036 * t716 - t646 * t1079 - t697 * t1087 - t903 * t645;
t575 = pkin(12) * t776 - t917;
t610 = -t646 * t890 + t893 * t697;
t999 = t891 * t1025;
t678 = qJD(4) * t856 + t878 * t1087 + t898 * t971 + (-qJD(3) * t893 - qJD(4)) * t975 + (t1001 * t891 + t938 + t999) * t903;
t588 = pkin(4) * t679 - pkin(12) * t678 + t610;
t627 = t722 * t1079 + t752 * t1087 + t903 * t716;
t614 = -pkin(12) * t921 + t627;
t668 = -t722 * t890 + t893 * t752;
t625 = pkin(4) * t767 - pkin(12) * t769 + t668;
t590 = t614 * t902 + t625 * t897;
t567 = -qJD(5) * t590 - t575 * t897 + t902 * t588;
t677 = qJDD(5) + t679;
t565 = -pkin(5) * t677 - t567;
t719 = t902 * t769 - t897 * t921;
t1083 = t891 * t899;
t941 = t1069 * t894 + t1068;
t788 = t1083 * t895 + t892 * t941;
t943 = t1063 * t894 - t1071;
t787 = t1082 * t895 + t892 * t943;
t839 = -t1080 * t891 + t894 * t895;
t954 = t787 * t893 + t839 * t890;
t683 = t788 * t903 + t898 * t954;
t738 = -t787 * t890 + t839 * t893;
t635 = t683 * t897 - t738 * t902;
t1102 = sin(pkin(14));
t992 = t1102 * t900;
t1103 = cos(pkin(14));
t993 = t1103 * t905;
t840 = t895 * t993 - t992;
t991 = t1102 * t905;
t994 = t1103 * t900;
t841 = t895 * t994 + t991;
t996 = t892 * t1103;
t970 = t891 * t996;
t739 = -t841 * t899 + (t840 * t894 - t970) * t904;
t740 = t1076 * t840 + t841 * t904 - t899 * t970;
t927 = -t840 * t891 - t894 * t996;
t915 = t927 * t890;
t638 = t740 * t903 + (t739 * t893 + t915) * t898;
t843 = -t895 * t992 + t993;
t842 = -t895 * t991 - t994;
t995 = t892 * t1102;
t925 = t842 * t894 + t891 * t995;
t741 = -t843 * t899 + t904 * t925;
t742 = t843 * t904 + t899 * t925;
t926 = -t842 * t891 + t894 * t995;
t914 = t926 * t890;
t640 = t742 * t903 + (t741 * t893 + t914) * t898;
t684 = -t739 * t890 + t893 * t927;
t685 = -t741 * t890 + t893 * t926;
t932 = g(1) * (-t640 * t897 + t685 * t902) + g(2) * (-t638 * t897 + t684 * t902) - g(3) * t635;
t1111 = t713 * (pkin(5) * t719 + pkin(13) * t713) + t565 + t932;
t1110 = -t898 * t783 + t903 * (t792 * t893 + t821 * t890);
t986 = t897 * t678 - t902 * t776;
t617 = qJD(5) * t719 + t986;
t615 = qJDD(6) + t617;
t871 = -pkin(5) * t902 - pkin(13) * t897 - pkin(4);
t1108 = (t627 - t765 * (pkin(5) * t897 - pkin(13) * t902)) * t713 - t871 * t615;
t1107 = -qJD(5) * t1115 + t1122 * t902 + t1123 * t897;
t906 = qJD(2) ^ 2;
t1105 = pkin(12) * qJD(5);
t901 = cos(qJ(6));
t1029 = qJD(6) * t901;
t616 = -qJD(5) * t807 - t1033 * t769 + t902 * t678 + t897 * t776;
t896 = sin(qJ(6));
t1017 = t765 * t1029 + t901 * t616 + t896 * t677;
t1030 = qJD(6) * t896;
t586 = -t1030 * t719 + t1017;
t1101 = t586 * t896;
t1096 = t719 * t896;
t669 = -t901 * t765 + t1096;
t1100 = t669 * t713;
t671 = t719 * t901 + t765 * t896;
t1099 = t671 * t713;
t1098 = t717 * t765;
t1097 = t719 * t765;
t1095 = t767 * t902;
t1094 = t788 * t898;
t1093 = t839 * t891;
t1091 = t878 * MDP(9);
t1089 = t890 * t891;
t1088 = t890 * t897;
t1074 = t896 * t615;
t1067 = t900 * t906;
t1066 = t901 * t615;
t1062 = qJDD(1) - g(3);
t1061 = -pkin(5) * t734 - t1107;
t695 = pkin(4) * t769 + pkin(12) * t767;
t1058 = t902 * t626 + t897 * t695;
t1055 = -t967 * pkin(5) + qJD(5) * t1048 + t1119 * t902 - t1121 * t897;
t797 = t1085 * t901 + t845 * t896;
t1054 = -qJD(6) * t797 + t1051 * t901 + t896 * t967;
t1020 = t896 * t1085;
t1053 = -qJD(6) * t1020 + t1029 * t845 + t1051 * t896 - t901 * t967;
t888 = t899 ^ 2;
t1043 = -t904 ^ 2 + t888;
t1034 = qJD(5) * t896;
t1032 = qJD(5) * t901;
t1019 = t891 * t1081;
t1008 = t891 * t1037;
t1002 = t891 * t894 * t906;
t566 = t625 * t1031 - t1033 * t614 + t902 * t575 + t897 * t588;
t564 = pkin(13) * t677 + t566;
t578 = -t722 * t1004 - t752 * t1006 - t716 * t1035 + t646 * t1078 + t697 * t1085 - t898 * t645;
t576 = -pkin(4) * t776 - t578;
t569 = pkin(5) * t617 - pkin(13) * t616 + t576;
t989 = -t896 * t564 + t901 * t569;
t987 = t616 * t896 - t901 * t677;
t985 = t765 * t902;
t984 = t713 * t901;
t982 = qJD(3) + 0.2e1 * t1039;
t981 = t878 + t1026;
t887 = t891 ^ 2;
t980 = t887 * t892 * t1067;
t976 = t891 * t1011;
t690 = -t1095 * t901 + t769 * t896;
t966 = t1031 * t901 - t690;
t666 = pkin(4) * t838 - t1110;
t743 = t789 * t897 + t902 * t838;
t744 = t789 * t902 - t838 * t897;
t611 = pkin(5) * t743 - pkin(13) * t744 + t666;
t965 = -pkin(13) * t734 - qJD(6) * t611 + t1116;
t599 = pkin(13) * t786 + t1115;
t650 = -qJD(5) * t743 + t733 * t902 + t897 * t974;
t651 = qJD(5) * t744 + t733 * t897 - t902 * t974;
t964 = -pkin(5) * t651 + pkin(13) * t650 + qJD(6) * t599 - t1059;
t830 = t879 + (-pkin(3) * t903 - pkin(4)) * t893;
t745 = pkin(5) * t844 - pkin(13) * t845 + t830;
t963 = -pkin(13) * t967 - qJD(6) * t745 + t1113;
t747 = -pkin(13) * t1085 + t1048;
t962 = -pkin(5) * t1050 + pkin(13) * t1051 + qJD(6) * t747 - t1047;
t960 = t901 * t564 + t896 * t569;
t583 = pkin(13) * t765 + t590;
t613 = pkin(4) * t921 - t626;
t591 = t717 * pkin(5) - t719 * pkin(13) + t613;
t571 = t583 * t901 + t591 * t896;
t959 = t583 * t896 - t591 * t901;
t589 = -t614 * t897 + t625 * t902;
t636 = t683 * t902 + t738 * t897;
t682 = -t1078 * t787 - t1085 * t839 + t1094;
t602 = t636 * t901 + t682 * t896;
t601 = -t636 * t896 + t682 * t901;
t957 = t654 * t902 - t667 * t897;
t687 = t744 * t901 + t786 * t896;
t686 = t744 * t896 - t901 * t786;
t952 = -t831 * t897 + t832 * t902;
t754 = -t1075 * t841 - t840 * t899;
t948 = t1089 * t841 + t754 * t893;
t756 = -t1075 * t843 - t842 * t899;
t947 = t1089 * t843 + t756 * t893;
t937 = -pkin(12) * t677 + t613 * t765;
t637 = -t1078 * t739 + t740 * t898 - t903 * t915;
t639 = -t1078 * t741 + t742 * t898 - t903 * t914;
t931 = g(1) * t639 + g(2) * t637 + g(3) * t682;
t930 = -g(1) * t640 - g(2) * t638 - g(3) * t683;
t928 = t1019 * t890 + t822 * t893;
t922 = -t576 + t931;
t731 = -t895 * t1009 + (qJD(2) * t942 - qJD(3) * t941) * t892;
t919 = t731 * t893 + t890 * t976;
t582 = -pkin(5) * t765 - t589;
t913 = -pkin(13) * t615 + (t582 + t589) * t713;
t912 = qJD(4) * t921;
t910 = pkin(12) * qJD(6) * t713 - t931;
t909 = t890 * t921 * t1083;
t908 = (pkin(13) * t769 - qJD(6) * t871 + t1058) * t713 + t930;
t817 = -t862 * t891 + t876;
t798 = t845 * t901 - t1020;
t782 = -t826 * t891 + t875;
t766 = t1019 * t893 - t822 * t890;
t757 = -t1076 * t843 + t842 * t904;
t755 = -t1076 * t841 + t840 * t904;
t746 = pkin(5) * t1085 - t952;
t732 = t895 * t1008 + (qJD(2) * t940 + qJD(3) * t943) * t892;
t721 = t823 * t903 + t898 * t928;
t720 = t823 * t898 - t903 * t928;
t715 = t1084 * t843 - t756 * t890;
t714 = t1084 * t841 - t754 * t890;
t712 = -t1079 * t788 + t787 * t903;
t711 = t1078 * t788 + t787 * t898;
t698 = -t731 * t890 + t893 * t976;
t689 = -t1095 * t896 - t901 * t769;
t673 = t1088 * t788 + t712 * t902;
t672 = t721 * t902 + t766 * t897;
t665 = -t1079 * t742 + t741 * t903;
t664 = t1078 * t742 + t741 * t898;
t663 = -t1079 * t740 + t739 * t903;
t662 = t1078 * t740 + t739 * t898;
t661 = t757 * t903 + t898 * t947;
t660 = t757 * t898 - t903 * t947;
t659 = t755 * t903 + t898 * t948;
t658 = t755 * t898 - t903 * t948;
t629 = t1088 * t742 + t665 * t902;
t628 = t1088 * t740 + t663 * t902;
t624 = t661 * t902 + t715 * t897;
t623 = t659 * t902 + t714 * t897;
t609 = t732 * t903 + t919 * t898 + (t903 * t954 - t1094) * qJD(4);
t608 = qJD(4) * t683 + t732 * t898 - t903 * t919;
t606 = t640 * t902 + t685 * t897;
t604 = t638 * t902 + t684 * t897;
t598 = -pkin(5) * t786 - t957;
t597 = qJD(6) * t687 + t650 * t896 - t901 * t734;
t596 = -qJD(6) * t686 + t650 * t901 + t734 * t896;
t594 = -pkin(5) * t769 + t626 * t897 - t695 * t902;
t587 = qJD(6) * t671 + t987;
t580 = -qJD(5) * t635 + t609 * t902 + t698 * t897;
t579 = qJD(5) * t636 + t609 * t897 - t698 * t902;
t563 = -t571 * qJD(6) + t989;
t562 = -qJD(6) * t959 + t960;
t1 = [t1062 * MDP(1) + (t1093 * t935 + t731 * t983 + t787 * t878 - t904 * t980) * MDP(10) + (t1093 * t936 - t732 * t983 - t788 * t878 + t899 * t980) * MDP(11) + (t608 * t921 + t738 * t679 - t682 * t776 + t698 * t767) * MDP(17) + (t609 * t921 + t738 * t678 - t683 * t776 + t698 * t769) * MDP(18) + (-t579 * t765 + t608 * t717 + t617 * t682 - t635 * t677) * MDP(24) + (-t580 * t765 + t608 * t719 + t616 * t682 - t636 * t677) * MDP(25) + ((-qJD(6) * t602 - t580 * t896 + t608 * t901) * t713 + t601 * t615 + t579 * t669 + t635 * t587) * MDP(31) + (-(qJD(6) * t601 + t580 * t901 + t608 * t896) * t713 - t602 * t615 + t579 * t671 + t635 * t586) * MDP(32) + ((qJDD(2) * t905 - t1067) * MDP(3) + (-qJDD(2) * t900 - t905 * t906) * MDP(4)) * t892; (-g(1) * t624 - g(2) * t623 - g(3) * t672 + t1059 * t717 + t1107 * t765 + t567 * t786 + t576 * t743 + t589 * t734 + t613 * t651 + t666 * t617 + t957 * t677) * MDP(24) + ((-t599 * t896 + t611 * t901) * t615 + t563 * t743 - t959 * t651 + t598 * t587 + t565 * t686 + t582 * t597 - g(1) * (t624 * t901 + t660 * t896) - g(2) * (t623 * t901 + t658 * t896) - g(3) * (t672 * t901 + t720 * t896) + (t896 * t965 - t901 * t964) * t713 + t1061 * t669) * MDP(31) + (-t1115 * t677 - t566 * t786 - t590 * t734 + t666 * t616 + t576 * t744 + t613 * t650 - g(1) * (-t661 * t897 + t715 * t902) - g(2) * (-t659 * t897 + t714 * t902) - g(3) * (-t721 * t897 + t766 * t902) + t1116 * t765 + t1059 * t719) * MDP(25) + (t884 * t878 + (-t848 * t1037 + t1049) * t894 - t812 * t983 - g(1) * t757 - g(2) * t755 - g(3) * t823 + (-pkin(10) * qJD(3) * t983 - t782) * t1082 + ((-t816 + ((-t862 - t1104) * t894 - pkin(2) * qJD(3)) * qJD(3)) * t894 + (-pkin(10) * t878 + (t817 - t876) * qJD(3)) * t891) * t899) * MDP(10) + ((qJDD(2) * t888 + 0.2e1 * t1001 * t899) * MDP(5) + (-pkin(2) * t936 - t899 * t972) * MDP(11) + 0.2e1 * (t1024 * t899 - t1028 * t1043) * MDP(6) + (-pkin(2) * t935 + t904 * t972) * MDP(10)) * t887 + (g(1) * t843 + g(2) * t841 - t1062 * t1081) * MDP(4) + (-g(1) * t842 - g(2) * t840 - g(3) * t1080 + t877) * MDP(3) + (t729 * t679 + t610 * t786 + t668 * t734 + t1110 * t776 - t578 * t838 + t626 * t974 - g(1) * t661 - g(2) * t659 - g(3) * t721 + t1052 * t767 + (-(t794 * t893 + t828 * t890) * t903 + t1118) * t921) * MDP(17) + (t1037 * t982 + t899 * t981) * t891 * MDP(7) + (-t1038 * t982 + t904 * t981) * t891 * MDP(8) + (-(t599 * t901 + t611 * t896) * t615 - t562 * t743 - t571 * t651 + t598 * t586 + t565 * t687 + t582 * t596 - g(1) * (-t624 * t896 + t660 * t901) - g(2) * (-t623 * t896 + t658 * t901) - g(3) * (-t672 * t896 + t720 * t901) + (t896 * t964 + t901 * t965) * t713 + t1061 * t671) * MDP(32) + (-t678 * t786 - t679 * t789 - t733 * t767 - t734 * t769) * MDP(13) + (t678 * t789 + t733 * t769) * MDP(12) + (t616 * t786 + t650 * t765 + t677 * t744 + t719 * t734) * MDP(21) + (-t617 * t786 - t651 * t765 - t677 * t743 - t717 * t734) * MDP(22) + (t677 * t786 + t734 * t765) * MDP(23) + (-t616 * t743 - t617 * t744 - t650 * t717 - t651 * t719) * MDP(20) + (t616 * t744 + t650 * t719) * MDP(19) + (-t587 * t743 - t597 * t713 - t615 * t686 - t651 * t669) * MDP(29) + (t586 * t743 + t596 * t713 + t615 * t687 + t651 * t671) * MDP(28) + (t615 * t743 + t651 * t713) * MDP(30) + (-t586 * t686 - t587 * t687 - t596 * t669 - t597 * t671) * MDP(27) + (t586 * t687 + t596 * t671) * MDP(26) + t894 * t1091 + (t838 * t679 + t734 * t921 - t767 * t974 - t776 * t786) * MDP(15) + (-t838 * t678 - t733 * t921 + t769 * t974 + t776 * t789) * MDP(14) + (-qJD(3) * t909 - t776 * t838) * MDP(16) + qJDD(2) * MDP(2) + (t782 * t1083 + t817 * t1008 - t1044 * t878 - t918 * t894 - g(1) * t756 - g(2) * t754 - g(3) * t822 + (pkin(10) * t1009 + t1128) * t983) * MDP(11) + (g(1) * t660 + g(2) * t658 + g(3) * t720 - t1015 * t776 + t1052 * t769 + t1126 * t921 + t610 * t789 - t627 * t974 + t668 * t733 + t729 * t678 - t917 * t838) * MDP(18); (t952 * t677 - t567 * t1085 + t830 * t617 + t576 * t844 - g(1) * t629 - g(2) * t628 - g(3) * t673 + ((-qJD(5) * t831 - t1119) * t902 + (-qJD(5) * t832 + t1121) * t897) * t765 + t1047 * t717 + t1050 * t613 + t967 * t589) * MDP(24) + (-MDP(5) * t899 * t904 + MDP(6) * t1043) * t887 * t906 + ((t745 * t901 - t747 * t896) * t615 + t563 * t844 + t746 * t587 + t565 * t797 - g(1) * (t629 * t901 + t664 * t896) - g(2) * (t628 * t901 + t662 * t896) - g(3) * (t673 * t901 + t711 * t896) + (t896 * t963 - t901 * t962) * t713 + t1055 * t669 + t1053 * t582 - t1050 * t959) * MDP(31) + (-t1048 * t677 + t566 * t1085 + t830 * t616 + t576 * t845 - g(1) * (t1086 * t742 - t665 * t897) - g(2) * (t1086 * t740 - t663 * t897) - g(3) * (t1086 * t788 - t712 * t897) + t1113 * t765 + t1047 * t719 + t1051 * t613 - t967 * t590) * MDP(25) + (-t1002 * t904 + t999) * MDP(7) + (t1002 * t899 + t998) * MDP(8) + (-t1085 * t677 + t765 * t967) * MDP(23) + (t1051 * t765 - t1085 * t616 + t677 * t845 + t719 * t967) * MDP(21) + (-t1050 * t765 + t1085 * t617 - t677 * t844 - t717 * t967) * MDP(22) + (g(1) * t742 + g(2) * t740 + g(3) * t788 + (-t817 * t1082 + (-t899 * t848 + t1046) * t894) * qJD(2) + t1046 * qJD(3) - t961) * MDP(11) + (-t890 * pkin(3) * t679 - t610 * t1085 + t1114 * t776 + t578 * t893 - t688 * t767 - t626 * t977 - g(1) * t665 - g(2) * t663 - g(3) * t712 + t967 * t668 + ((t736 * t893 + t827 * t890) * t903 + t1117) * t921) * MDP(17) + t1091 + (-g(1) * t741 - g(2) * t739 - g(3) * t787 - t899 * t816 + (t894 * t1064 + (-t817 * t891 + t894 * t939) * t899) * qJD(2) + t1049) * MDP(10) + (t1054 * t671 + t586 * t798) * MDP(26) + (-t1053 * t671 - t1054 * t669 - t586 * t797 - t587 * t798) * MDP(27) + (-(t745 * t896 + t747 * t901) * t615 - t562 * t844 + t746 * t586 + t565 * t798 - g(1) * (-t629 * t896 + t664 * t901) - g(2) * (-t628 * t896 + t662 * t901) - g(3) * (-t673 * t896 + t711 * t901) + (t896 * t962 + t901 * t963) * t713 + t1055 * t671 + t1054 * t582 - t1050 * t571) * MDP(32) + (t1051 * t719 + t616 * t845) * MDP(19) + (-t1050 * t719 - t1051 * t717 - t616 * t844 - t617 * t845) * MDP(20) + (-t1050 * t669 - t1053 * t713 - t587 * t844 - t615 * t797) * MDP(29) + (t1050 * t671 + t1054 * t713 + t586 * t844 + t615 * t798) * MDP(28) + (t1050 * t713 + t615 * t844) * MDP(30) + (qJD(2) * t909 + t776 * t893) * MDP(16) + (-t893 * t679 - t921 * t811 + (t1012 * t767 + t776 * t903 + t898 * t912) * t890) * MDP(15) + (t893 * t678 + t921 * t813 + (-t1012 * t769 + t776 * t898 - t903 * t912) * t890) * MDP(14) + (t767 * t813 + t769 * t811 + (t678 * t903 - t679 * t898 + (-t767 * t903 - t769 * t898) * qJD(4)) * t890) * MDP(13) + (-t1045 * t776 + t917 * t893 - t688 * t769 - t668 * t813 + g(1) * t664 + g(2) * t662 + g(3) * t711 + (-pkin(3) * t678 + t1012 * t627 + t1035 * t668 + t610 * t898) * t890 + t1124 * t921) * MDP(18) + (t678 * t1087 + t1125 * t769) * MDP(12); -t767 ^ 2 * MDP(13) + (-t767 * t921 + t678) * MDP(14) - t679 * MDP(15) + t776 * MDP(16) + (-t627 * t921 + t578 + t931) * MDP(17) + (-t626 * t921 + t668 * t767 + t917 - t930) * MDP(18) + (t616 * t897 + t719 * t985) * MDP(19) + ((t616 - t1098) * t902 + (-t617 - t1097) * t897) * MDP(20) + (t677 * t897 + t765 * t985) * MDP(21) + (-t765 ^ 2 * t897 + t677 * t902) * MDP(22) + (-pkin(4) * t617 - t627 * t717 + (t626 * t765 + t937) * t897 + ((-t695 - t1105) * t765 + t922) * t902) * MDP(24) + (-pkin(4) * t616 + t1058 * t765 - t627 * t719 + t937 * t902 + (t1105 * t765 - t922) * t897) * MDP(25) + (t586 * t897 * t901 + (-t1030 * t897 + t966) * t671) * MDP(26) + (t669 * t690 + t671 * t689 + (-t669 * t901 - t671 * t896) * t1031 + (-t1101 - t587 * t901 + (t669 * t896 - t671 * t901) * qJD(6)) * t897) * MDP(27) + (-t586 * t902 + t966 * t713 + (-t1030 * t713 + t671 * t765 + t1066) * t897) * MDP(28) + (t587 * t902 + (-t1031 * t896 + t689) * t713 + (-t1029 * t713 - t669 * t765 - t1074) * t897) * MDP(29) + (t713 * t765 * t897 - t615 * t902) * MDP(30) + (-t582 * t689 - t594 * t669 - t1108 * t901 + t908 * t896 + (t582 * t1034 - t563 + (qJD(5) * t669 - t1074) * pkin(12) - t910 * t901) * t902 + (t582 * t1029 + t565 * t896 - t765 * t959 + (t1034 * t713 + t587) * pkin(12)) * t897) * MDP(31) + (-t582 * t690 - t594 * t671 + t1108 * t896 + t908 * t901 + (t582 * t1032 + t562 + (qJD(5) * t671 - t1066) * pkin(12) + t910 * t896) * t902 + (-t582 * t1030 + t565 * t901 - t765 * t571 + (t1032 * t713 + t586) * pkin(12)) * t897) * MDP(32) + (MDP(12) * t767 + MDP(13) * t769 - MDP(15) * t921 - t668 * MDP(17) - t719 * MDP(21) + t717 * MDP(22) - t765 * MDP(23) - t589 * MDP(24) + t590 * MDP(25)) * t769; -t717 ^ 2 * MDP(20) + (t616 + t1098) * MDP(21) + (t1097 - t986) * MDP(22) + t677 * MDP(23) + (t590 * t765 + t567 - t932) * MDP(24) + (g(1) * t606 + g(2) * t604 + g(3) * t636 + t589 * t765 + t613 * t717 - t566) * MDP(25) + (t671 * t984 + t1101) * MDP(26) + ((t586 - t1100) * t901 + (-t587 - t1099) * t896) * MDP(27) + (t713 * t984 + t1074) * MDP(28) + (-t713 ^ 2 * t896 + t1066) * MDP(29) + (-pkin(5) * t587 - t1111 * t901 - t590 * t669 + t913 * t896) * MDP(31) + (-pkin(5) * t586 + t1111 * t896 - t590 * t671 + t913 * t901) * MDP(32) + (MDP(19) * t717 + MDP(20) * t719 - MDP(22) * qJD(5) - MDP(24) * t613 - MDP(28) * t671 + MDP(29) * t669 - MDP(30) * t713 + MDP(31) * t959 + MDP(32) * t571) * t719; t671 * t669 * MDP(26) + (-t669 ^ 2 + t671 ^ 2) * MDP(27) + (t1017 + t1100) * MDP(28) + (-t987 + t1099) * MDP(29) + t615 * MDP(30) + (t571 * t713 - t582 * t671 - g(1) * (-t606 * t896 + t639 * t901) - g(2) * (-t604 * t896 + t637 * t901) - g(3) * t601 + t989) * MDP(31) + (-t959 * t713 + t582 * t669 - g(1) * (-t606 * t901 - t639 * t896) - g(2) * (-t604 * t901 - t637 * t896) + g(3) * t602 - t960) * MDP(32) + (-MDP(28) * t1096 - MDP(29) * t671 - MDP(31) * t571 + MDP(32) * t959) * qJD(6);];
tau  = t1;