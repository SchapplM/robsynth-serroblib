% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 13:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRPRPR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_invdynB_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:04:28
% EndTime: 2019-05-06 13:05:12
% DurationCPUTime: 43.50s
% Computational Cost: add. (148674->831), mult. (351981->1220), div. (0->0), fcn. (262294->10), ass. (0->586)
t1044 = sin(qJ(2));
t1048 = cos(qJ(2));
t1040 = sin(pkin(10));
t1041 = cos(pkin(10));
t1047 = cos(qJ(4));
t1043 = sin(qJ(4));
t1036 = qJDD(2) + qJDD(4);
t1109 = qJD(1) * t1048;
t1110 = qJD(1) * t1044;
t995 = -t1040 * t1110 + t1041 * t1109;
t996 = (t1040 * t1048 + t1041 * t1044) * qJD(1);
t956 = t1043 * t996 - t1047 * t995;
t958 = t1043 * t995 + t1047 * t996;
t1142 = t958 * t956;
t1158 = t1142 + t1036;
t1175 = t1043 * t1158;
t1037 = qJD(2) + qJD(4);
t1035 = t1037 ^ 2;
t1148 = t958 ^ 2;
t928 = t1148 + t1035;
t830 = t1047 * t928 + t1175;
t1171 = t1047 * t1158;
t833 = t1043 * t928 - t1171;
t750 = t1040 * t833 - t1041 * t830;
t752 = t1040 * t830 + t1041 * t833;
t684 = t1044 * t752 + t1048 * t750;
t1219 = pkin(7) * t684;
t685 = t1044 * t750 - t1048 * t752;
t1218 = pkin(7) * t685;
t1045 = sin(qJ(1));
t1217 = t1045 * t685;
t881 = -t1142 + t1036;
t1191 = t1043 * t881;
t934 = -t1148 + t1035;
t838 = t1047 * t934 + t1191;
t1190 = t1047 * t881;
t843 = t1043 * t934 - t1190;
t757 = t1040 * t843 - t1041 * t838;
t761 = t1040 * t838 + t1041 * t843;
t690 = t1044 * t757 - t1048 * t761;
t1216 = t1045 * t690;
t1149 = t956 ^ 2;
t933 = t1149 - t1035;
t840 = t1043 * t933 + t1171;
t844 = t1047 * t933 - t1175;
t758 = t1040 * t844 + t1041 * t840;
t763 = t1040 * t840 - t1041 * t844;
t693 = t1044 * t758 + t1048 * t763;
t1215 = t1045 * t693;
t1049 = cos(qJ(1));
t1214 = t1049 * t685;
t1213 = t1049 * t690;
t1212 = t1049 * t693;
t1211 = pkin(1) * t684 + pkin(2) * t750 - pkin(3) * t830;
t879 = -t1035 - t1149;
t802 = t1043 * t879 + t1190;
t805 = -t1047 * t879 + t1191;
t718 = t1040 * t805 - t1041 * t802;
t720 = t1040 * t802 + t1041 * t805;
t647 = t1044 * t720 + t1048 * t718;
t1210 = pkin(7) * t647;
t648 = t1044 * t718 - t1048 * t720;
t1209 = pkin(7) * t648;
t1208 = t1045 * t648;
t1207 = t1049 * t648;
t1206 = pkin(1) * t647 + pkin(2) * t718 - pkin(3) * t802;
t1205 = qJ(3) * t750;
t1204 = qJ(3) * t752;
t1202 = t1044 * t761 + t1048 * t757;
t1201 = t1044 * t763 - t1048 * t758;
t1200 = qJ(3) * t718;
t1199 = qJ(3) * t720;
t1196 = pkin(8) * t830;
t1195 = pkin(8) * t833;
t1193 = pkin(8) * t802;
t1192 = pkin(8) * t805;
t1182 = 2 * qJD(5);
t1074 = qJD(2) * pkin(3) - pkin(8) * t996;
t1083 = qJD(2) * t1109;
t1095 = qJDD(1) * t1044;
t1005 = t1083 + t1095;
t1031 = t1048 * qJDD(1);
t1084 = qJD(2) * t1110;
t1006 = t1031 - t1084;
t1078 = t1005 * t1040 - t1041 * t1006;
t1039 = t1048 ^ 2;
t1051 = qJD(1) ^ 2;
t1065 = qJD(2) * pkin(2) - qJ(3) * t1110;
t1015 = t1045 * g(1) - t1049 * g(2);
t1066 = qJDD(1) * pkin(1) + t1015;
t916 = t1006 * pkin(2) - t1065 * t1110 - qJDD(3) + t1066 + (qJ(3) * t1039 + pkin(7)) * t1051;
t993 = t995 ^ 2;
t834 = -pkin(3) * t1078 + t993 * pkin(8) - t996 * t1074 + t916;
t940 = t1037 * t958;
t1189 = pkin(4) * t940 - t1182 * t958 - t834;
t1061 = (-t1043 * t956 - t1047 * t958) * t1037;
t1062 = (t1043 * t958 - t1047 * t956) * t1037;
t1154 = -t1040 * t1061 + t1041 * t1062;
t1155 = t1040 * t1062 + t1041 * t1061;
t1164 = -t1044 * t1155 + t1048 * t1154;
t1188 = -t1049 * t1036 + t1045 * t1164;
t1086 = t1049 * t1142;
t1101 = t1037 * t1047;
t968 = t1005 * t1041 + t1006 * t1040;
t1079 = t1043 * t968 + t1047 * t1078;
t864 = qJD(4) * t958 + t1079;
t1063 = t1043 * t864 + t1101 * t956;
t1102 = t1037 * t1043;
t1072 = -t1047 * t864 + t1102 * t956;
t1152 = t1040 * t1063 + t1041 * t1072;
t1153 = -t1040 * t1072 + t1041 * t1063;
t1165 = -t1044 * t1152 + t1048 * t1153;
t1187 = t1045 * t1165 + t1086;
t1058 = -t1043 * t1078 + t1047 * t968;
t865 = -t956 * qJD(4) + t1058;
t1070 = t1043 * t865 + t1101 * t958;
t1071 = t1047 * t865 - t1102 * t958;
t1150 = -t1040 * t1070 + t1041 * t1071;
t1151 = t1040 * t1071 + t1041 * t1070;
t1168 = -t1044 * t1151 + t1048 * t1150;
t1186 = t1045 * t1168 - t1086;
t1089 = t1045 * t1142;
t1185 = t1049 * t1165 - t1089;
t1184 = t1049 * t1168 + t1089;
t1183 = t1036 * t1045 + t1049 * t1164;
t1156 = -t1149 - t1148;
t1181 = pkin(1) * t1156;
t1180 = pkin(2) * t1156;
t1179 = pkin(3) * t1156;
t967 = t995 * t996;
t1159 = qJDD(2) + t967;
t1178 = t1040 * t1159;
t1177 = t1041 * t1159;
t1042 = sin(qJ(6));
t862 = qJDD(6) + t865;
t1046 = cos(qJ(6));
t917 = t1037 * t1042 - t1046 * t956;
t919 = t1037 * t1046 + t1042 * t956;
t871 = t919 * t917;
t1162 = -t871 + t862;
t1176 = t1042 * t1162;
t1174 = t1045 * t1156;
t1157 = -t1148 + t1149;
t1173 = t1045 * t1157;
t1172 = t1046 * t1162;
t1170 = t1049 * t1156;
t1169 = t1049 * t1157;
t1167 = t1044 * t1150 + t1048 * t1151;
t1166 = t1044 * t1153 + t1048 * t1152;
t1163 = t1044 * t1154 + t1048 * t1155;
t939 = t956 * t1037;
t826 = -t939 - t865;
t825 = -t939 + t865;
t1161 = qJ(5) * t825;
t914 = t917 ^ 2;
t915 = t919 ^ 2;
t951 = qJD(6) + t958;
t949 = t951 ^ 2;
t994 = t996 ^ 2;
t1147 = 2 * qJD(3);
t1146 = pkin(4) + pkin(9);
t1145 = pkin(4) * t1043;
t1144 = pkin(4) * t1047;
t1143 = t917 * t951;
t1099 = t1044 * t1051;
t1111 = qJD(1) * qJD(2);
t1016 = g(1) * t1049 + g(2) * t1045;
t999 = -pkin(1) * t1051 + qJDD(1) * pkin(7) - t1016;
t1122 = t1044 * t999;
t908 = qJDD(2) * pkin(2) - t1005 * qJ(3) - t1122 + (pkin(2) * t1099 + qJ(3) * t1111 - g(3)) * t1048;
t1033 = t1039 * t1051;
t982 = -t1044 * g(3) + t1048 * t999;
t909 = -pkin(2) * t1033 + t1006 * qJ(3) - qJD(2) * t1065 + t982;
t836 = t1040 * t909 - t1041 * t908 + t1147 * t996;
t1141 = qJD(2) * t995;
t925 = -t968 + t1141;
t779 = pkin(3) * t1159 + pkin(8) * t925 - t836;
t837 = t1040 * t908 + t1041 * t909 + t995 * t1147;
t790 = -t993 * pkin(3) - pkin(8) * t1078 - qJD(2) * t1074 + t837;
t702 = t1043 * t790 - t1047 * t779;
t703 = t1043 * t779 + t1047 * t790;
t1140 = qJD(2) * t996;
t1138 = t1036 * qJ(5);
t628 = t1043 * t703 - t1047 * t702;
t1137 = t1040 * t628;
t1136 = t1040 * t916;
t963 = qJDD(2) - t967;
t1135 = t1040 * t963;
t1134 = t1041 * t628;
t1133 = t1041 * t916;
t1132 = t1041 * t963;
t886 = pkin(4) * t956 - qJ(5) * t958;
t1067 = -t1035 * pkin(4) - t956 * t886 + t703;
t930 = pkin(5) * t958 - pkin(9) * t1037;
t636 = t1138 - t864 * pkin(5) - t1149 * pkin(9) + (t1182 + t930) * t1037 + t1067;
t1131 = t1042 * t636;
t787 = t871 + t862;
t1130 = t1042 * t787;
t1129 = t1042 * t951;
t1128 = t1043 * t834;
t754 = t1040 * t837 - t1041 * t836;
t1124 = t1044 * t754;
t998 = t1051 * pkin(7) + t1066;
t1123 = t1044 * t998;
t1121 = t1046 * t636;
t1120 = t1046 * t787;
t1119 = t1046 * t951;
t819 = t864 + t940;
t1118 = t1047 * t819;
t1117 = t1047 * t834;
t1113 = t1048 * t754;
t1112 = t1048 * t998;
t1108 = qJD(2) * t1040;
t1107 = qJD(2) * t1041;
t1022 = t1048 * t1099;
t1013 = qJDD(2) + t1022;
t1106 = t1013 * t1044;
t1014 = qJDD(2) - t1022;
t1105 = t1014 * t1044;
t1104 = t1014 * t1048;
t1038 = t1044 ^ 2;
t1100 = t1038 * t1051;
t1098 = qJD(4) - t1037;
t1097 = qJD(4) + t1037;
t1096 = t1038 + t1039;
t1094 = qJDD(1) * t1045;
t1093 = qJDD(1) * t1049;
t1092 = qJDD(2) * t1049;
t1091 = -t915 - t949;
t1090 = t1043 * t871;
t1088 = t1045 * t967;
t1087 = t1047 * t871;
t1085 = t1049 * t967;
t801 = -t917 * qJD(6) + t1046 * t1036 + t1042 * t864;
t1082 = qJ(5) * t1043 + pkin(3);
t755 = t1040 * t836 + t1041 * t837;
t629 = t1043 * t702 + t1047 * t703;
t981 = t1048 * g(3) + t1122;
t912 = t1044 * t981 + t1048 * t982;
t1077 = t1042 * t1036 - t1046 * t864;
t974 = -t1015 * t1045 - t1049 * t1016;
t1076 = t1045 * t1022;
t1075 = t1049 * t1022;
t1010 = -t1045 * t1051 + t1093;
t1073 = -pkin(6) * t1010 - g(3) * t1045;
t1064 = -t1036 * pkin(4) - t1035 * qJ(5) + t958 * t886 + qJDD(5) + t702;
t635 = -pkin(5) * t826 - pkin(9) * t881 + t1064;
t1052 = -t1161 + t1189;
t645 = -pkin(5) * t1149 + t1146 * t864 - t958 * t930 + t1052;
t584 = t1042 * t645 - t1046 * t635;
t585 = t1042 * t635 + t1046 * t645;
t546 = t1042 * t585 - t1046 * t584;
t547 = t1042 * t584 + t1046 * t585;
t910 = t1044 * t982 - t1048 * t981;
t973 = t1015 * t1049 - t1016 * t1045;
t1068 = t801 - t1143;
t923 = -t1078 + t1140;
t1060 = (-qJD(6) + t951) * t919 - t1077;
t1059 = -0.2e1 * qJD(5) * t1037 - t1067;
t682 = -t1059 + t1138;
t1056 = -t1097 * t956 + t1058;
t1053 = -t864 * pkin(4) - t1189;
t1050 = qJD(2) ^ 2;
t1030 = t1045 * qJDD(2);
t1021 = -t1033 - t1050;
t1020 = t1033 - t1050;
t1019 = -t1050 - t1100;
t1018 = t1050 - t1100;
t1012 = t1033 - t1100;
t1011 = t1033 + t1100;
t1009 = t1049 * t1051 + t1094;
t1008 = t1096 * qJDD(1);
t1007 = t1031 - 0.2e1 * t1084;
t1004 = 0.2e1 * t1083 + t1095;
t1002 = t1048 * t1013;
t1001 = t1096 * t1111;
t992 = -pkin(6) * t1009 + g(3) * t1049;
t987 = -t994 - t1050;
t986 = -t994 + t1050;
t985 = t993 - t1050;
t984 = t1005 * t1048 - t1038 * t1111;
t983 = -t1006 * t1044 - t1039 * t1111;
t980 = -t1019 * t1044 - t1104;
t979 = -t1018 * t1044 + t1002;
t978 = t1021 * t1048 - t1106;
t977 = t1020 * t1048 - t1105;
t976 = t1019 * t1048 - t1105;
t975 = t1021 * t1044 + t1002;
t971 = t1008 * t1049 - t1011 * t1045;
t970 = t1008 * t1045 + t1011 * t1049;
t969 = -t1004 * t1044 + t1007 * t1048;
t966 = -t994 + t993;
t961 = -t1050 - t993;
t946 = (t1040 * t996 + t1041 * t995) * qJD(2);
t945 = (t1040 * t995 - t1041 * t996) * qJD(2);
t944 = t1004 * t1045 + t1049 * t980;
t943 = -t1007 * t1045 + t1049 * t978;
t942 = -t1004 * t1049 + t1045 * t980;
t941 = t1007 * t1049 + t1045 * t978;
t932 = -pkin(7) * t976 - t1112;
t931 = -pkin(7) * t975 - t1123;
t927 = -pkin(1) * t976 + t982;
t926 = -pkin(1) * t975 + t981;
t924 = t968 + t1141;
t921 = t1078 + t1140;
t920 = -t993 - t994;
t907 = t1041 * t968 - t1108 * t996;
t906 = t1040 * t968 + t1107 * t996;
t905 = t1040 * t1078 - t1107 * t995;
t904 = -t1041 * t1078 - t1108 * t995;
t900 = -t1040 * t987 - t1132;
t899 = -t1040 * t986 + t1177;
t898 = t1041 * t985 - t1135;
t897 = t1041 * t987 - t1135;
t896 = t1041 * t986 + t1178;
t895 = t1040 * t985 + t1132;
t894 = -t1045 * t998 + t1049 * t912;
t893 = t1045 * t912 + t1049 * t998;
t892 = t1041 * t961 - t1178;
t891 = t1040 * t961 + t1177;
t878 = -t915 + t949;
t877 = t914 - t949;
t876 = -t1044 * t945 + t1048 * t946;
t870 = t915 - t914;
t869 = -t1040 * t925 + t1041 * t923;
t868 = -t1040 * t924 - t1041 * t921;
t867 = t1040 * t923 + t1041 * t925;
t866 = -t1040 * t921 + t1041 * t924;
t856 = -qJ(3) * t897 - t1133;
t853 = -t1044 * t906 + t1048 * t907;
t852 = -t1044 * t904 + t1048 * t905;
t850 = -t1044 * t897 + t1048 * t900;
t849 = -t1044 * t896 + t1048 * t899;
t848 = -t1044 * t895 + t1048 * t898;
t847 = t1044 * t900 + t1048 * t897;
t846 = -qJ(3) * t891 - t1136;
t828 = -t949 - t914;
t827 = -t914 - t915;
t821 = -t1098 * t956 + t1058;
t820 = -t864 + t940;
t818 = t1098 * t958 + t1079;
t817 = t1097 * t958 + t1079;
t808 = -t1044 * t891 + t1048 * t892;
t807 = t1044 * t892 + t1048 * t891;
t806 = -pkin(2) * t924 + qJ(3) * t900 - t1136;
t800 = -qJD(6) * t919 - t1077;
t799 = -pkin(2) * t921 + qJ(3) * t892 + t1133;
t798 = t1045 * t924 + t1049 * t850;
t797 = t1045 * t850 - t1049 * t924;
t796 = (-t1042 * t919 + t1046 * t917) * t951;
t795 = (t1042 * t917 + t1046 * t919) * t951;
t784 = t1045 * t921 + t1049 * t808;
t783 = -t1044 * t867 + t1048 * t869;
t782 = -t1044 * t866 + t1048 * t868;
t781 = t1045 * t808 - t1049 * t921;
t780 = t1044 * t869 + t1048 * t867;
t775 = t801 + t1143;
t771 = (qJD(6) + t951) * t919 + t1077;
t769 = -t1046 * t801 + t1129 * t919;
t768 = -t1042 * t801 - t1119 * t919;
t767 = t1042 * t800 - t1119 * t917;
t766 = -t1046 * t800 - t1129 * t917;
t765 = t1045 * t920 + t1049 * t783;
t764 = t1045 * t783 - t1049 * t920;
t753 = -t1117 + t1196;
t748 = -t1043 * t795 + t1047 * t862;
t747 = t1043 * t862 + t1047 * t795;
t746 = -t1046 * t877 + t1130;
t745 = t1042 * t878 - t1172;
t744 = -t1042 * t877 - t1120;
t743 = -t1046 * t878 - t1176;
t742 = -t1043 * t826 + t1047 * t820;
t741 = -t1043 * t1056 - t1118;
t740 = -t1043 * t825 - t1047 * t817;
t739 = t1043 * t821 - t1047 * t818;
t738 = t1043 * t820 + t1047 * t826;
t737 = -t1043 * t819 + t1047 * t1056;
t736 = -t1043 * t817 + t1047 * t825;
t735 = -t1043 * t818 - t1047 * t821;
t734 = -t1128 - t1193;
t725 = -pkin(1) * t780 - pkin(2) * t867;
t724 = -pkin(1) * t847 - pkin(2) * t897 + t837;
t723 = pkin(2) * t916 + qJ(3) * t755;
t722 = -t1042 * t1091 - t1120;
t721 = t1046 * t1091 - t1130;
t716 = t1046 * t828 - t1176;
t715 = t1042 * t828 + t1172;
t714 = -t1043 * t768 + t1087;
t713 = -t1043 * t766 - t1087;
t712 = t1047 * t768 + t1090;
t711 = t1047 * t766 - t1090;
t710 = -pkin(1) * t807 - pkin(2) * t891 + t836;
t709 = -qJ(3) * t867 - t754;
t706 = -pkin(7) * t847 - t1044 * t806 + t1048 * t856;
t705 = -pkin(3) * t825 - t1128 + t1195;
t704 = -pkin(2) * t920 + qJ(3) * t869 + t755;
t700 = -pkin(3) * t817 + t1117 - t1192;
t699 = -pkin(7) * t807 - t1044 * t799 + t1048 * t846;
t698 = t1053 + t1161;
t697 = t1042 * t775 + t1046 * t1060;
t696 = t1042 * t1068 + t1046 * t771;
t695 = t1042 * t1060 - t1046 * t775;
t694 = t1042 * t771 - t1046 * t1068;
t689 = t1048 * t755 - t1124;
t688 = t1044 * t755 + t1113;
t681 = -t1043 * t743 + t1047 * t775;
t680 = -t1043 * t744 + t1047 * t1060;
t679 = t1043 * t775 + t1047 * t743;
t678 = t1043 * t1060 + t1047 * t744;
t677 = -t1045 * t916 + t1049 * t689;
t676 = t1045 * t689 + t1049 * t916;
t675 = -t1040 * t747 + t1041 * t748;
t674 = t1040 * t748 + t1041 * t747;
t673 = t1043 * t721 + t1047 * t1068;
t672 = t1043 * t1068 - t1047 * t721;
t671 = t1043 * t715 + t1047 * t771;
t670 = t1043 * t771 - t1047 * t715;
t669 = -t1040 * t738 + t1041 * t742;
t668 = -t1040 * t737 + t1041 * t741;
t667 = -t1040 * t736 + t1041 * t740;
t666 = -t1040 * t735 + t1041 * t739;
t665 = t1040 * t742 + t1041 * t738;
t664 = t1040 * t741 + t1041 * t737;
t663 = t1040 * t740 + t1041 * t736;
t662 = t1040 * t739 + t1041 * t735;
t661 = -t1043 * t694 + t1047 * t870;
t660 = t1043 * t870 + t1047 * t694;
t659 = t1052 + (t864 + t819) * pkin(4);
t658 = t1053 + (t825 + t1056) * qJ(5);
t653 = -qJ(5) * t1156 + t1064;
t652 = t1043 * t695 + t1047 * t827;
t651 = t1043 * t827 - t1047 * t695;
t650 = -pkin(4) * t1156 + t682;
t644 = t1045 * t825 - t1214;
t643 = -t1045 * t1056 + t1214;
t642 = -t1049 * t825 - t1217;
t641 = t1049 * t1056 + t1217;
t640 = -t1040 * t712 + t1041 * t714;
t639 = -t1040 * t711 + t1041 * t713;
t638 = t1040 * t714 + t1041 * t712;
t637 = t1040 * t713 + t1041 * t711;
t634 = -pkin(1) * t688 - pkin(2) * t754;
t633 = -t1045 * t819 - t1207;
t632 = t1045 * t817 + t1207;
t631 = t1049 * t819 - t1208;
t630 = -t1049 * t817 + t1208;
t627 = t1047 * t658 - t1056 * t1145 - t1196;
t626 = qJ(5) * t1118 - t1043 * t659 + t1193;
t625 = -t1040 * t705 + t1041 * t753 - t1205;
t624 = pkin(3) * t834 + pkin(8) * t629;
t623 = -pkin(7) * t780 - t1044 * t704 + t1048 * t709;
t622 = pkin(5) * t695 - qJ(5) * t697;
t621 = -t1195 + t1043 * t658 + (pkin(3) + t1144) * t1056;
t620 = t1047 * t659 + t1082 * t819 + t1192;
t619 = -t1040 * t700 + t1041 * t734 + t1200;
t618 = t1043 * t1064 + t1047 * t682;
t617 = t1043 * t682 - t1047 * t1064;
t616 = -pkin(2) * t825 + t1040 * t753 + t1041 * t705 + t1204;
t615 = -t1040 * t679 + t1041 * t681;
t614 = -t1040 * t678 + t1041 * t680;
t613 = t1040 * t681 + t1041 * t679;
t612 = t1040 * t680 + t1041 * t678;
t611 = -pkin(7) * t688 - qJ(3) * t1113 - t1044 * t723;
t610 = -pkin(8) * t738 - t628;
t609 = -t1044 * t674 + t1048 * t675;
t608 = -t1040 * t672 + t1041 * t673;
t607 = t1040 * t673 + t1041 * t672;
t606 = -t1040 * t670 + t1041 * t671;
t605 = t1040 * t671 + t1041 * t670;
t604 = -t1044 * t665 + t1048 * t669;
t603 = -t1044 * t664 + t1048 * t668;
t602 = -t1044 * t663 + t1048 * t667;
t601 = -t1044 * t662 + t1048 * t666;
t600 = t1044 * t669 + t1048 * t665;
t599 = t1044 * t666 + t1048 * t662;
t598 = -t1040 * t660 + t1041 * t661;
t597 = t1040 * t661 + t1041 * t660;
t596 = -pkin(2) * t817 + t1040 * t734 + t1041 * t700 - t1199;
t595 = pkin(8) * t742 - t1179 + t629;
t594 = -t1040 * t651 + t1041 * t652;
t593 = t1040 * t652 + t1041 * t651;
t592 = t1049 * t604 + t1174;
t591 = t1049 * t601 + t1174;
t590 = t1045 * t604 - t1170;
t589 = t1045 * t601 - t1170;
t588 = -t1211 + t703;
t587 = -t1044 * t638 + t1048 * t640;
t586 = -t1044 * t637 + t1048 * t639;
t582 = t1206 + t702;
t581 = -pkin(4) * t928 + (-t1036 - t1158) * qJ(5) + t1059 + t1211;
t580 = -pkin(8) * t735 - t1043 * t650 + t1047 * t653;
t579 = pkin(5) * t1068 - t1146 * t722 - t1131;
t578 = t1041 * t629 - t1137;
t577 = t1040 * t629 + t1134;
t576 = pkin(5) * t771 - t1146 * t716 + t1121;
t575 = pkin(8) * t739 + t1043 * t653 + t1047 * t650 - t1179;
t574 = pkin(4) * t881 + qJ(5) * t879 - t1064 - t1206;
t573 = -pkin(8) * t617 + (qJ(5) * t1047 - t1145) * t698;
t572 = -pkin(1) * t600 - pkin(2) * t665 - pkin(3) * t738;
t571 = pkin(5) * t721 - qJ(5) * t722 - t585;
t570 = pkin(5) * t715 - qJ(5) * t716 - t584;
t569 = -t1040 * t617 + t1041 * t618;
t568 = t1040 * t618 + t1041 * t617;
t567 = -t1040 * t621 + t1041 * t627 + t1205;
t566 = -t1044 * t613 + t1048 * t615;
t565 = -t1044 * t612 + t1048 * t614;
t564 = -t1040 * t620 + t1041 * t626 - t1200;
t563 = pkin(2) * t1056 + t1040 * t627 + t1041 * t621 - t1204;
t562 = -t1044 * t607 + t1048 * t608;
t561 = t1044 * t608 + t1048 * t607;
t560 = -pkin(1) * t599 - pkin(2) * t662 - pkin(3) * t735 + pkin(4) * t821 + qJ(5) * t818;
t559 = -t1044 * t605 + t1048 * t606;
t558 = t1044 * t606 + t1048 * t605;
t557 = -t1044 * t597 + t1048 * t598;
t556 = pkin(8) * t618 + (t1082 + t1144) * t698;
t555 = pkin(2) * t819 + t1040 * t626 + t1041 * t620 + t1199;
t554 = -t1044 * t616 + t1048 * t625 - t1219;
t553 = -t1044 * t593 + t1048 * t594;
t552 = t1044 * t594 + t1048 * t593;
t551 = t1045 * t722 + t1049 * t562;
t550 = t1045 * t562 - t1049 * t722;
t549 = t1045 * t716 + t1049 * t559;
t548 = t1045 * t559 - t1049 * t716;
t545 = -qJ(3) * t665 - t1040 * t595 + t1041 * t610;
t544 = -t1044 * t596 + t1048 * t619 + t1210;
t543 = qJ(3) * t669 + t1040 * t610 + t1041 * t595 - t1180;
t542 = t1045 * t697 + t1049 * t553;
t541 = t1045 * t553 - t1049 * t697;
t540 = -t1044 * t577 + t1048 * t578;
t539 = t1044 * t578 + t1048 * t577;
t538 = -pkin(8) * t1134 - qJ(3) * t577 - t1040 * t624;
t537 = -t1045 * t834 + t1049 * t540;
t536 = t1045 * t540 + t1049 * t834;
t535 = pkin(2) * t834 - pkin(8) * t1137 + qJ(3) * t578 + t1041 * t624;
t534 = t1043 * t546 + t1047 * t636;
t533 = t1043 * t636 - t1047 * t546;
t532 = -qJ(3) * t662 - t1040 * t575 + t1041 * t580;
t531 = qJ(3) * t666 + t1040 * t580 + t1041 * t575 - t1180;
t530 = pkin(5) * t827 - t1146 * t697 - t547;
t529 = -t1044 * t568 + t1048 * t569;
t528 = t1044 * t569 + t1048 * t568;
t527 = -pkin(8) * t672 - t1043 * t579 + t1047 * t571;
t526 = -pkin(8) * t670 - t1043 * t576 + t1047 * t570;
t525 = -pkin(3) * t722 + pkin(8) * t673 + t1043 * t571 + t1047 * t579;
t524 = -pkin(3) * t716 + pkin(8) * t671 + t1043 * t570 + t1047 * t576;
t523 = -t1045 * t698 + t1049 * t529;
t522 = t1045 * t529 + t1049 * t698;
t521 = -t1044 * t563 + t1048 * t567 + t1219;
t520 = -t1044 * t555 + t1048 * t564 - t1210;
t519 = -pkin(8) * t651 - t1043 * t530 + t1047 * t622;
t518 = -pkin(1) * t561 - pkin(2) * t607 - pkin(3) * t672 - qJ(5) * t1068 + t1146 * t721 - t1121;
t517 = -pkin(1) * t539 - pkin(2) * t577 - pkin(3) * t628;
t516 = pkin(5) * t546 - qJ(5) * t547;
t515 = -pkin(1) * t558 - pkin(2) * t605 - pkin(3) * t670 - qJ(5) * t771 + t1146 * t715 - t1131;
t514 = -pkin(3) * t697 + pkin(8) * t652 + t1043 * t622 + t1047 * t530;
t513 = pkin(5) * t636 - t1146 * t547;
t512 = -qJ(3) * t568 - t1040 * t556 + t1041 * t573;
t511 = pkin(2) * t698 + qJ(3) * t569 + t1040 * t573 + t1041 * t556;
t510 = -pkin(7) * t600 - t1044 * t543 + t1048 * t545;
t509 = -t1040 * t533 + t1041 * t534;
t508 = t1040 * t534 + t1041 * t533;
t507 = -pkin(1) * t552 - pkin(2) * t593 - pkin(3) * t651 - qJ(5) * t827 + t1146 * t695 + t546;
t506 = -pkin(1) * t528 - pkin(2) * t568 - pkin(3) * t617 + pkin(4) * t1064 - qJ(5) * t682;
t505 = -pkin(7) * t599 - t1044 * t531 + t1048 * t532;
t504 = -qJ(3) * t607 - t1040 * t525 + t1041 * t527;
t503 = -qJ(3) * t605 - t1040 * t524 + t1041 * t526;
t502 = -pkin(7) * t539 - t1044 * t535 + t1048 * t538;
t501 = -pkin(2) * t722 + qJ(3) * t608 + t1040 * t527 + t1041 * t525;
t500 = -pkin(2) * t716 + qJ(3) * t606 + t1040 * t526 + t1041 * t524;
t499 = -qJ(3) * t593 - t1040 * t514 + t1041 * t519;
t498 = -pkin(2) * t697 + qJ(3) * t594 + t1040 * t519 + t1041 * t514;
t497 = -t1044 * t508 + t1048 * t509;
t496 = t1044 * t509 + t1048 * t508;
t495 = -pkin(8) * t533 - t1043 * t513 + t1047 * t516;
t494 = t1045 * t547 + t1049 * t497;
t493 = t1045 * t497 - t1049 * t547;
t492 = -pkin(7) * t528 - t1044 * t511 + t1048 * t512;
t491 = -pkin(3) * t547 + pkin(8) * t534 + t1043 * t516 + t1047 * t513;
t490 = -pkin(7) * t561 - t1044 * t501 + t1048 * t504;
t489 = -pkin(7) * t558 - t1044 * t500 + t1048 * t503;
t488 = -pkin(7) * t552 - t1044 * t498 + t1048 * t499;
t487 = -pkin(1) * t496 - pkin(2) * t508 - pkin(3) * t533 - qJ(5) * t636 + t1146 * t546;
t486 = -qJ(3) * t508 - t1040 * t491 + t1041 * t495;
t485 = -pkin(2) * t547 + qJ(3) * t509 + t1040 * t495 + t1041 * t491;
t484 = -pkin(7) * t496 - t1044 * t485 + t1048 * t486;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1009, -t1010, 0, t974, 0, 0, 0, 0, 0, 0, t943, t944, t971, t894, 0, 0, 0, 0, 0, 0, t784, t798, t765, t677, 0, 0, 0, 0, 0, 0, t632, t644, t592, t537, 0, 0, 0, 0, 0, 0, t591, t633, t643, t523, 0, 0, 0, 0, 0, 0, t549, t551, t542, t494; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1010, -t1009, 0, t973, 0, 0, 0, 0, 0, 0, t941, t942, t970, t893, 0, 0, 0, 0, 0, 0, t781, t797, t764, t676, 0, 0, 0, 0, 0, 0, t630, t642, t590, t536, 0, 0, 0, 0, 0, 0, t589, t631, t641, t522, 0, 0, 0, 0, 0, 0, t548, t550, t541, t493; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t975, t976, 0, t910, 0, 0, 0, 0, 0, 0, t807, t847, t780, t688, 0, 0, 0, 0, 0, 0, -t647, t684, t600, t539, 0, 0, 0, 0, 0, 0, t599, t647, -t684, t528, 0, 0, 0, 0, 0, 0, t558, t561, t552, t496; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1010, 0, -t1009, 0, t1073, -t992, -t973, -pkin(6) * t973, t1049 * t984 - t1076, -t1012 * t1045 + t1049 * t969, t1044 * t1094 + t1049 * t979, t1049 * t983 + t1076, t1031 * t1045 + t1049 * t977, t1001 * t1049 + t1030, -pkin(6) * t941 - t1045 * t926 + t1049 * t931, -pkin(6) * t942 - t1045 * t927 + t1049 * t932, -pkin(6) * t970 - t1049 * t910, -pkin(6) * t893 + (pkin(1) * t1045 - pkin(7) * t1049) * t910, t1049 * t853 - t1088, -t1045 * t966 + t1049 * t782, -t1045 * t925 + t1049 * t849, t1049 * t852 + t1088, t1045 * t923 + t1049 * t848, t1049 * t876 + t1030, -pkin(6) * t781 - t1045 * t710 + t1049 * t699, -pkin(6) * t797 - t1045 * t724 + t1049 * t706, -pkin(6) * t764 - t1045 * t725 + t1049 * t623, -pkin(6) * t676 - t1045 * t634 + t1049 * t611, t1184, t1049 * t602 - t1173, -t1045 * t826 + t1213, t1185, -t1045 * t818 - t1212, t1183, -pkin(6) * t630 - t1045 * t582 + t1049 * t544, -pkin(6) * t642 - t1045 * t588 + t1049 * t554, -pkin(6) * t590 - t1045 * t572 + t1049 * t510, -pkin(6) * t536 - t1045 * t517 + t1049 * t502, t1183, -t1045 * t821 - t1213, -t1045 * t820 + t1212, t1184, t1049 * t603 - t1173, t1185, -pkin(6) * t589 - t1045 * t560 + t1049 * t505, -pkin(6) * t631 - t1045 * t574 + t1049 * t520, -pkin(6) * t641 - t1045 * t581 + t1049 * t521, -pkin(6) * t522 - t1045 * t506 + t1049 * t492, -t1045 * t769 + t1049 * t587, -t1045 * t696 + t1049 * t557, -t1045 * t745 + t1049 * t566, -t1045 * t767 + t1049 * t586, -t1045 * t746 + t1049 * t565, -t1045 * t796 + t1049 * t609, -pkin(6) * t548 - t1045 * t515 + t1049 * t489, -pkin(6) * t550 - t1045 * t518 + t1049 * t490, -pkin(6) * t541 - t1045 * t507 + t1049 * t488, -pkin(6) * t493 - t1045 * t487 + t1049 * t484; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1009, 0, t1010, 0, t992, t1073, t974, pkin(6) * t974, t1045 * t984 + t1075, t1012 * t1049 + t1045 * t969, -t1044 * t1093 + t1045 * t979, t1045 * t983 - t1075, -t1031 * t1049 + t1045 * t977, t1001 * t1045 - t1092, pkin(6) * t943 + t1045 * t931 + t1049 * t926, pkin(6) * t944 + t1045 * t932 + t1049 * t927, pkin(6) * t971 - t1045 * t910, pkin(6) * t894 + (-pkin(1) * t1049 - pkin(7) * t1045) * t910, t1045 * t853 + t1085, t1045 * t782 + t1049 * t966, t1045 * t849 + t1049 * t925, t1045 * t852 - t1085, t1045 * t848 - t1049 * t923, t1045 * t876 - t1092, pkin(6) * t784 + t1045 * t699 + t1049 * t710, pkin(6) * t798 + t1045 * t706 + t1049 * t724, pkin(6) * t765 + t1045 * t623 + t1049 * t725, pkin(6) * t677 + t1045 * t611 + t1049 * t634, t1186, t1045 * t602 + t1169, t1049 * t826 + t1216, t1187, t1049 * t818 - t1215, t1188, pkin(6) * t632 + t1045 * t544 + t1049 * t582, pkin(6) * t644 + t1045 * t554 + t1049 * t588, pkin(6) * t592 + t1045 * t510 + t1049 * t572, pkin(6) * t537 + t1045 * t502 + t1049 * t517, t1188, t1049 * t821 - t1216, t1049 * t820 + t1215, t1186, t1045 * t603 + t1169, t1187, pkin(6) * t591 + t1045 * t505 + t1049 * t560, pkin(6) * t633 + t1045 * t520 + t1049 * t574, pkin(6) * t643 + t1045 * t521 + t1049 * t581, pkin(6) * t523 + t1045 * t492 + t1049 * t506, t1045 * t587 + t1049 * t769, t1045 * t557 + t1049 * t696, t1045 * t566 + t1049 * t745, t1045 * t586 + t1049 * t767, t1045 * t565 + t1049 * t746, t1045 * t609 + t1049 * t796, pkin(6) * t549 + t1045 * t489 + t1049 * t515, pkin(6) * t551 + t1045 * t490 + t1049 * t518, pkin(6) * t542 + t1045 * t488 + t1049 * t507, pkin(6) * t494 + t1045 * t484 + t1049 * t487; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1015, t1016, 0, 0, (t1005 + t1083) * t1044, t1004 * t1048 + t1007 * t1044, t1018 * t1048 + t1106, (t1006 - t1084) * t1048, t1020 * t1044 + t1104, 0, pkin(1) * t1007 + pkin(7) * t978 + t1112, -pkin(1) * t1004 + pkin(7) * t980 - t1123, pkin(1) * t1011 + pkin(7) * t1008 + t912, pkin(1) * t998 + pkin(7) * t912, t1044 * t907 + t1048 * t906, t1044 * t868 + t1048 * t866, t1044 * t899 + t1048 * t896, t1044 * t905 + t1048 * t904, t1044 * t898 + t1048 * t895, t1044 * t946 + t1048 * t945, -pkin(1) * t921 + pkin(7) * t808 + t1044 * t846 + t1048 * t799, -pkin(1) * t924 + pkin(7) * t850 + t1044 * t856 + t1048 * t806, -pkin(1) * t920 + pkin(7) * t783 + t1044 * t709 + t1048 * t704, pkin(1) * t916 + pkin(7) * t689 - qJ(3) * t1124 + t1048 * t723, t1167, t1044 * t667 + t1048 * t663, -t1202, t1166, -t1201, t1163, -pkin(1) * t817 + t1044 * t619 + t1048 * t596 + t1209, -pkin(1) * t825 + t1044 * t625 + t1048 * t616 - t1218, pkin(7) * t604 + t1044 * t545 + t1048 * t543 - t1181, pkin(1) * t834 + pkin(7) * t540 + t1044 * t538 + t1048 * t535, t1163, t1202, t1201, t1167, t1044 * t668 + t1048 * t664, t1166, pkin(7) * t601 + t1044 * t532 + t1048 * t531 - t1181, pkin(1) * t819 + t1044 * t564 + t1048 * t555 - t1209, pkin(1) * t1056 + t1044 * t567 + t1048 * t563 + t1218, pkin(1) * t698 + pkin(7) * t529 + t1044 * t512 + t1048 * t511, t1044 * t640 + t1048 * t638, t1044 * t598 + t1048 * t597, t1044 * t615 + t1048 * t613, t1044 * t639 + t1048 * t637, t1044 * t614 + t1048 * t612, t1044 * t675 + t1048 * t674, -pkin(1) * t716 + pkin(7) * t559 + t1044 * t503 + t1048 * t500, -pkin(1) * t722 + pkin(7) * t562 + t1044 * t504 + t1048 * t501, -pkin(1) * t697 + pkin(7) * t553 + t1044 * t499 + t1048 * t498, -pkin(1) * t547 + pkin(7) * t497 + t1044 * t486 + t1048 * t485;];
tauB_reg  = t1;
