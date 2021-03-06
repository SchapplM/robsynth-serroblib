% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 03:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RRPRRR14V3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14V3_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_invdynm_fixb_reg2_snew_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 03:54:48
% EndTime: 2019-05-07 03:55:07
% DurationCPUTime: 12.29s
% Computational Cost: add. (31842->581), mult. (62272->728), div. (0->0), fcn. (48066->10), ass. (0->407)
t1117 = sin(qJ(6));
t1119 = sin(qJ(4));
t1124 = cos(qJ(4));
t1120 = sin(qJ(2));
t1191 = qJD(1) * t1120;
t1080 = qJD(2) * t1119 + t1124 * t1191;
t1125 = cos(qJ(2));
t1190 = qJD(1) * t1125;
t1101 = qJD(4) - t1190;
t1118 = sin(qJ(5));
t1123 = cos(qJ(5));
t1049 = t1080 * t1123 + t1101 * t1118;
t1078 = -t1124 * qJD(2) + t1119 * t1191;
t1072 = qJD(5) + t1078;
t1122 = cos(qJ(6));
t1003 = t1049 * t1117 - t1122 * t1072;
t1005 = t1049 * t1122 + t1072 * t1117;
t952 = t1005 * t1003;
t1106 = qJD(2) * t1190;
t1168 = t1120 * qJDD(1);
t1084 = t1106 + t1168;
t1031 = -qJD(4) * t1078 + qJDD(2) * t1119 + t1084 * t1124;
t1105 = qJD(2) * t1191;
t1167 = t1125 * qJDD(1);
t1085 = -t1105 + t1167;
t1076 = -qJDD(4) + t1085;
t1148 = -t1118 * t1031 - t1123 * t1076;
t1134 = qJD(5) * t1049 - t1148;
t962 = qJDD(6) + t1134;
t1221 = -t952 + t962;
t1230 = t1117 * t1221;
t1050 = t1080 * t1078;
t1219 = -t1050 - t1076;
t1229 = t1119 * t1219;
t1228 = t1122 * t1221;
t1227 = t1124 * t1219;
t1121 = sin(qJ(1));
t1126 = cos(qJ(1));
t1093 = g(1) * t1121 - t1126 * g(2);
t1216 = 2 * qJD(3);
t1226 = t1084 + t1106;
t1218 = qJ(3) * t1226 + t1191 * t1216 + t1093;
t1129 = t1125 * t1218;
t1083 = 0.2e1 * t1106 + t1168;
t1086 = -0.2e1 * t1105 + t1167;
t1182 = t1086 * t1125;
t1039 = t1083 * t1120 - t1182;
t1116 = t1125 ^ 2;
t1127 = qJD(1) ^ 2;
t1110 = t1116 * t1127;
t1115 = t1120 ^ 2;
t1177 = t1115 * t1127;
t1089 = -t1110 + t1177;
t1225 = t1039 * t1121 + t1089 * t1126;
t1224 = t1039 * t1126 - t1089 * t1121;
t1100 = t1125 * t1127 * t1120;
t1092 = qJDD(2) - t1100;
t1217 = qJD(2) ^ 2;
t1097 = t1110 - t1217;
t1057 = t1092 * t1120 - t1097 * t1125;
t1155 = t1126 * t1167;
t1223 = t1057 * t1121 + t1155;
t1156 = t1121 * t1167;
t1222 = t1057 * t1126 - t1156;
t1147 = t1124 * qJDD(2) - t1119 * t1084;
t1135 = qJD(4) * t1080 - t1147;
t1023 = qJDD(5) + t1135;
t1047 = t1080 * t1118 - t1123 * t1101;
t998 = t1049 * t1047;
t957 = t1023 - t998;
t1220 = t1118 * t957;
t1094 = g(1) * t1126 + g(2) * t1121;
t1069 = -t1120 * g(3) - t1125 * t1094;
t1132 = -qJ(3) * t1100 + qJD(2) * t1216 + t1069;
t1033 = qJDD(2) * qJ(3) + t1132;
t978 = t1119 * t1033 + t1124 * t1218;
t1194 = t1124 * t978;
t979 = t1124 * t1033 - t1119 * t1218;
t1150 = -t1119 * t979 + t1194;
t1001 = t1003 ^ 2;
t1002 = t1005 ^ 2;
t1042 = qJD(6) + t1047;
t1040 = t1042 ^ 2;
t1045 = t1047 ^ 2;
t1046 = t1049 ^ 2;
t1071 = t1072 ^ 2;
t1074 = t1078 ^ 2;
t1098 = t1101 ^ 2;
t1215 = t1121 * g(3);
t1214 = t1126 * g(3);
t1068 = t1125 * g(3) - t1120 * t1094;
t1043 = qJDD(3) + t1068 + (-t1177 - t1217) * qJ(3);
t947 = t1043 * t1118 + t1123 * t979;
t896 = t1117 * t947 - t1122 * t978;
t897 = t1117 * t978 + t1122 * t947;
t1143 = -t1117 * t897 + t1122 * t896;
t946 = -t1123 * t1043 + t1118 * t979;
t1203 = t1118 * t946;
t839 = t1117 * t896 + t1122 * t897;
t826 = t1123 * t839 + t1203;
t793 = -t1119 * t1143 + t1124 * t826;
t1213 = qJ(3) * t793;
t1142 = t1123 * t947 + t1203;
t1200 = t1119 * t978;
t860 = t1124 * t1142 + t1200;
t1212 = qJ(3) * t860;
t913 = t1124 * t979 + t1200;
t1211 = qJ(3) * t913;
t1197 = t1122 * t946;
t1196 = t1123 * t946;
t857 = t1117 * t1196 - t1118 * t896;
t1210 = -t1119 * t1197 + t1124 * t857;
t1201 = t1119 * t946;
t858 = -t1118 * t897 + t1122 * t1196;
t1209 = t1117 * t1201 + t1124 * t858;
t1208 = qJ(3) * t1033;
t1207 = qJ(3) * t1120;
t1206 = qJ(3) * t1125;
t905 = t952 + t962;
t1205 = t1117 * t905;
t1204 = t1117 * t946;
t969 = t1118 * t978;
t889 = -t1118 * t947 + t1196;
t1202 = t1119 * t889;
t1198 = t1122 * t905;
t958 = t1023 + t998;
t1195 = t1123 * t958;
t971 = t1123 * t978;
t945 = t1124 * t946;
t1193 = t1125 * t978;
t1192 = qJD(1) * qJD(2);
t1017 = -t1050 + t1076;
t1189 = t1017 * t1119;
t1188 = t1017 * t1124;
t1187 = t1042 * t1117;
t1186 = t1042 * t1122;
t1185 = t1072 * t1118;
t1184 = t1072 * t1123;
t1183 = t1086 * t1120;
t1179 = t1101 * t1119;
t1178 = t1101 * t1124;
t1176 = t1118 * t1120;
t1175 = t1118 * t1125;
t1174 = t1119 * t1043;
t1173 = t1119 * t1123;
t1081 = t1120 * t1093;
t1172 = t1121 * t1125;
t1171 = t1123 * t1124;
t1035 = t1124 * t1043;
t1000 = qJ(3) * t1083 + t1218;
t1170 = t1125 * t1000;
t1082 = t1125 * t1093;
t1169 = t1125 * t1126;
t792 = t1119 * t826 + t1124 * t1143;
t1166 = t792 * t1207;
t859 = t1119 * t1142 - t1194;
t1165 = t859 * t1207;
t1164 = t1150 * t1207;
t1163 = qJ(3) * t1167;
t1162 = t1118 * t952;
t1161 = t1123 * t952;
t1160 = t1119 * t998;
t1159 = t1124 * t998;
t1158 = t1120 * t1050;
t1157 = t1125 * t1050;
t1154 = t1119 * t839 + t1143 * t1171;
t1153 = t1119 * t857 + t1122 * t945;
t1152 = -t1117 * t945 + t1119 * t858;
t1151 = -t1119 * t947 + t978 * t1171;
t1140 = -t1123 * t1031 + t1118 * t1076;
t964 = -qJD(5) * t1047 - t1140;
t1149 = t1122 * t1023 - t1117 * t964;
t1146 = t1121 * t1100;
t1145 = t1126 * t1100;
t1144 = t1118 * t1200 + t945;
t1141 = -t1117 * t1023 - t1122 * t964;
t1036 = t1083 * t1125 + t1183;
t1053 = t1092 * t1125 + t1097 * t1120;
t1138 = t1118 * t1194 - t1201;
t1137 = -t1124 * t839 + t1143 * t1173;
t1136 = t1124 * t947 + t1173 * t978;
t907 = -qJD(6) * t1005 + t1149;
t908 = -qJD(6) * t1003 - t1141;
t1130 = qJ(3) * t1218;
t1128 = t1120 * t1130;
t1096 = t1177 - t1217;
t1091 = qJDD(2) + t1100;
t1088 = qJDD(1) * t1126 - t1121 * t1127;
t1087 = qJDD(1) * t1121 + t1126 * t1127;
t1077 = (t1115 + t1116) * t1192;
t1075 = t1080 ^ 2;
t1066 = t1101 * t1080;
t1065 = t1101 * t1078;
t1064 = -t1075 + t1098;
t1063 = t1074 - t1098;
t1062 = qJDD(2) * t1121 + t1077 * t1126;
t1061 = t1084 * t1125 - t1115 * t1192;
t1060 = -qJDD(2) * t1126 + t1077 * t1121;
t1059 = -t1085 * t1120 - t1116 * t1192;
t1058 = t1091 * t1125 + t1096 * t1120;
t1055 = t1091 * t1120 - t1096 * t1125;
t1052 = t1226 * t1120;
t1051 = (t1085 - t1105) * t1125;
t1044 = t1075 - t1074;
t1041 = -t1075 - t1098;
t1032 = -t1098 - t1074;
t1029 = t1061 * t1126 - t1146;
t1028 = t1059 * t1126 + t1146;
t1027 = t1061 * t1121 + t1145;
t1026 = t1059 * t1121 - t1145;
t1025 = t1058 * t1126 + t1121 * t1168;
t1024 = t1058 * t1121 - t1126 * t1168;
t1022 = -qJ(3) * (-t1110 - t1217) + t1043;
t1016 = -qJ(3) * (-t1110 - t1177) + t1043;
t1014 = t1072 * t1049;
t1013 = t1072 * t1047;
t1012 = t1068 * t1125 - t1069 * t1120;
t1011 = t1068 * t1120 + t1069 * t1125;
t1010 = -t1046 + t1071;
t1009 = t1045 - t1071;
t1008 = (qJDD(2) + t1092) * qJ(3) + t1132;
t1007 = (-t1078 * t1124 + t1080 * t1119) * t1101;
t1006 = (-t1078 * t1119 - t1080 * t1124) * t1101;
t999 = t1120 * t1000;
t997 = t1031 + t1065;
t996 = t1031 - t1065;
t995 = (-qJD(4) + t1101) * t1080 + t1147;
t994 = -t1066 - t1135;
t993 = -t1066 + t1135;
t992 = t1046 - t1045;
t991 = qJ(3) * t1182 - t1120 * t1218;
t990 = qJ(3) * t1183 + t1129;
t989 = t1031 * t1124 - t1080 * t1179;
t988 = t1031 * t1119 + t1080 * t1178;
t987 = t1078 * t1178 + t1119 * t1135;
t986 = -t1078 * t1179 + t1124 * t1135;
t985 = t1007 * t1125 - t1076 * t1120;
t984 = t1007 * t1120 + t1076 * t1125;
t983 = t1063 * t1124 + t1189;
t982 = -t1064 * t1119 + t1227;
t981 = t1063 * t1119 - t1188;
t980 = t1064 * t1124 + t1229;
t977 = t1042 * t1005;
t976 = t1042 * t1003;
t975 = -t1002 + t1040;
t974 = t1001 - t1040;
t968 = t1016 * t1125 - t1033 * t1120;
t967 = t1016 * t1120 + t1033 * t1125;
t965 = t1045 + t1046;
t960 = (-t1047 * t1123 + t1049 * t1118) * t1072;
t959 = (-t1047 * t1118 - t1049 * t1123) * t1072;
t956 = t1125 * t989 + t1158;
t955 = t1125 * t987 - t1158;
t954 = t1120 * t989 - t1157;
t953 = t1120 * t987 + t1157;
t951 = t1002 - t1001;
t950 = -t1002 - t1040;
t949 = t1035 - qJ(3) * (t1041 * t1124 + t1189);
t948 = t1174 + qJ(3) * (-t1041 * t1119 + t1188);
t944 = -t1035 + qJ(3) * (t1032 * t1124 - t1229);
t943 = t1174 - qJ(3) * (t1032 * t1119 + t1227);
t938 = -t1119 * t996 + t1124 * t994;
t937 = t1119 * t994 + t1124 * t996;
t936 = (qJD(5) + t1072) * t1047 + t1140;
t935 = t1013 + t964;
t934 = -t1013 + t964;
t933 = -t1014 - t1134;
t932 = -t1014 + t1134;
t931 = -t1040 - t1001;
t930 = -t1120 * t993 + t1125 * t983;
t929 = t1120 * t997 + t1125 * t982;
t928 = t1120 * t983 + t1125 * t993;
t927 = t1120 * t982 - t1125 * t997;
t926 = -t1049 * t1185 + t1123 * t964;
t925 = t1049 * t1184 + t1118 * t964;
t924 = t1047 * t1184 + t1118 * t1134;
t923 = -t1047 * t1185 + t1123 * t1134;
t922 = t1023 * t1119 + t1124 * t960;
t921 = -t1023 * t1124 + t1119 * t960;
t920 = t1009 * t1123 - t1118 * t958;
t919 = -t1010 * t1118 + t1123 * t957;
t918 = t1009 * t1118 + t1195;
t917 = t1010 * t1123 + t1220;
t916 = t1044 * t1120 + t1125 * t938;
t915 = -t1044 * t1125 + t1120 * t938;
t914 = -t1195 - t1118 * (-t1046 - t1071);
t911 = (-t1003 * t1122 + t1005 * t1117) * t1042;
t910 = (-t1003 * t1117 - t1005 * t1122) * t1042;
t909 = t1123 * (-t1071 - t1045) - t1220;
t903 = -t1120 * t979 + t1125 * t949;
t902 = t1120 * t949 + t1125 * t979;
t901 = t1124 * t926 + t1160;
t900 = t1124 * t924 - t1160;
t899 = t1119 * t926 - t1159;
t898 = t1119 * t924 + t1159;
t894 = -t1120 * t978 + t1125 * t943;
t893 = t1120 * t943 + t1193;
t892 = t1120 * t959 + t1125 * t922;
t891 = t1120 * t922 - t1125 * t959;
t887 = t1124 * t889;
t886 = t1118 * t962 + t1123 * t911;
t885 = t1118 * t911 - t1123 * t962;
t884 = t908 + t976;
t883 = t908 - t976;
t882 = (-qJD(6) + t1042) * t1005 + t1149;
t881 = t907 - t977;
t880 = -t907 - t977;
t879 = t1123 * ((-qJD(5) + t1072) * t1049 + t1148) + t1118 * t935;
t878 = -t1118 * t934 + t1123 * t933;
t877 = t1118 * t933 + t1123 * t934;
t876 = -t1005 * t1187 + t1122 * t908;
t875 = t1005 * t1186 + t1117 * t908;
t874 = t1003 * t1186 - t1117 * t907;
t873 = -t1003 * t1187 - t1122 * t907;
t872 = t1122 * t974 - t1205;
t871 = -t1117 * t975 + t1228;
t870 = t1117 * t974 + t1198;
t869 = t1122 * t975 + t1230;
t868 = -t1119 * t932 + t1124 * t920;
t867 = t1119 * t935 + t1124 * t919;
t866 = t1119 * t920 + t1124 * t932;
t865 = t1119 * t919 - t1124 * t935;
t864 = -qJ(3) * (t1119 * t995 - t1124 * t997) + t1150;
t863 = qJ(3) * (t1119 * t997 + t1124 * t995) + t913;
t862 = t1120 * t864;
t861 = -t1122 * t950 + t1205;
t856 = t1118 * t1197 + t1123 * t897;
t855 = -t1117 * t1203 - t1123 * t896;
t854 = -t1117 * t931 - t1228;
t851 = t1119 * t992 + t1124 * t878;
t850 = t1119 * t878 - t1124 * t992;
t849 = t1120 * t925 + t1125 * t901;
t848 = -t1120 * t923 + t1125 * t900;
t847 = t1120 * t901 - t1125 * t925;
t846 = t1120 * t900 + t1125 * t923;
t845 = t1123 * t876 + t1162;
t844 = t1123 * t874 - t1162;
t843 = t1118 * t876 - t1161;
t842 = t1118 * t874 + t1161;
t838 = t1118 * t1143;
t836 = t1119 * t910 + t1124 * t886;
t835 = t1119 * t886 - t1124 * t910;
t834 = t1120 * t918 + t1125 * t868;
t833 = t1120 * t917 + t1125 * t867;
t832 = t1120 * t868 - t1125 * t918;
t831 = t1120 * t867 - t1125 * t917;
t830 = -qJ(3) * (t1119 * t914 + t1124 * t936) + t1151;
t829 = qJ(3) * (-t1119 * t936 + t1124 * t914) + t1136;
t828 = -qJ(3) * (t1119 * t909 + t1124 * t933) + t1138;
t827 = -qJ(3) * (-t1119 * t933 + t1124 * t909) - t1144;
t825 = -t1117 * t883 + t1122 * t881;
t824 = -t1117 * t882 + t1122 * t884;
t823 = t1117 * t881 + t1122 * t883;
t822 = -t1118 * t880 + t1123 * t872;
t821 = t1118 * t884 + t1123 * t871;
t820 = t1118 * t872 + t1123 * t880;
t819 = t1118 * t871 - t1123 * t884;
t818 = t1123 * (-t1117 * t950 - t1198) - t1118 * ((qJD(6) + t1042) * t1003 + t1141);
t817 = t1125 * t830 + t1176 * t978;
t816 = t1120 * t830 - t1175 * t978;
t815 = t1123 * (t1122 * t931 - t1230) - t1118 * t881;
t814 = t887 - qJ(3) * (t1119 * t879 + t1124 * t965);
t813 = t1202 + qJ(3) * (-t1119 * t965 + t1124 * t879);
t812 = -t1120 * t971 + t1125 * t828;
t811 = t1120 * t828 + t1123 * t1193;
t810 = t1120 * t877 + t1125 * t851;
t809 = t1120 * t851 - t1125 * t877;
t808 = t1118 * t951 + t1123 * t825;
t807 = t1118 * t825 - t1123 * t951;
t806 = t1119 * t875 + t1124 * t845;
t805 = -t1119 * t873 + t1124 * t844;
t804 = t1119 * t845 - t1124 * t875;
t803 = t1119 * t844 + t1124 * t873;
t802 = t1123 * (t1117 * t884 + t1122 * t882) - t1118 * (t1001 + t1002);
t801 = t1120 * t885 + t1125 * t836;
t800 = t1120 * t836 - t1125 * t885;
t799 = t1119 * t870 + t1124 * t822;
t798 = t1119 * t869 + t1124 * t821;
t797 = t1119 * t822 - t1124 * t870;
t796 = t1119 * t821 - t1124 * t869;
t795 = t1120 * t1142 + t1125 * t814;
t794 = t1120 * t814 - t1125 * t1142;
t791 = t1120 * t843 + t1125 * t806;
t790 = t1120 * t842 + t1125 * t805;
t789 = t1120 * t806 - t1125 * t843;
t788 = t1120 * t805 - t1125 * t842;
t787 = t1119 * t823 + t1124 * t808;
t786 = t1119 * t808 - t1124 * t823;
t785 = -qJ(3) * (t1119 * t818 + t1124 * t861) + t1209;
t784 = qJ(3) * (-t1119 * t861 + t1124 * t818) + t1152;
t783 = -qJ(3) * (t1119 * t815 + t1124 * t854) + t1210;
t782 = qJ(3) * (-t1119 * t854 + t1124 * t815) + t1153;
t781 = t1120 * t820 + t1125 * t799;
t780 = t1120 * t819 + t1125 * t798;
t779 = t1120 * t799 - t1125 * t820;
t778 = t1120 * t798 - t1125 * t819;
t777 = t1120 * t856 + t1125 * t785;
t776 = t1120 * t785 - t1125 * t856;
t775 = -t1120 * t855 + t1125 * t783;
t774 = t1120 * t783 + t1125 * t855;
t773 = t1120 * t807 + t1125 * t787;
t772 = t1120 * t787 - t1125 * t807;
t771 = -qJ(3) * (t1119 * t802 + t1124 * t824) + t1154;
t770 = qJ(3) * (-t1119 * t824 + t1124 * t802) + t1137;
t769 = t1125 * t771 + t1143 * t1176;
t768 = t1120 * t771 - t1143 * t1175;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1088, 0, -t1087, 0, -t1215, -t1214, -t1093 * t1126 + t1094 * t1121, 0, t1029, -t1224, t1025, t1028, -t1222, t1062, -t1068 * t1121 - t1081 * t1126, -t1069 * t1121 - t1093 * t1169, t1126 * t1012, 0, t1029, t1025, t1224, t1062, t1222, t1028, -t1022 * t1121 + t1126 * t991, qJ(3) * t1156 + t1126 * t968, t1000 * t1169 + t1008 * t1121, (t1121 * t1033 + t1126 * t1129) * qJ(3), t1121 * t988 + t1126 * t956, t1121 * t937 + t1126 * t916, t1121 * t980 + t1126 * t929, -t1121 * t986 + t1126 * t955, t1121 * t981 + t1126 * t930, t1006 * t1121 + t1126 * t985, t1121 * t944 + t1126 * t894, t1121 * t948 + t1126 * t903, t1121 * t863 + t1169 * t864, (t1121 * t913 + t1150 * t1169) * qJ(3), t1121 * t899 + t1126 * t849, t1121 * t850 + t1126 * t810, t1121 * t865 + t1126 * t833, t1121 * t898 + t1126 * t848, t1121 * t866 + t1126 * t834, t1121 * t921 + t1126 * t892, -t1121 * t827 + t1126 * t812, t1121 * t829 + t1126 * t817, t1121 * t813 + t1126 * t795, (t1121 * t860 - t1169 * t859) * qJ(3), t1121 * t804 + t1126 * t791, t1121 * t786 + t1126 * t773, t1121 * t796 + t1126 * t780, t1121 * t803 + t1126 * t790, t1121 * t797 + t1126 * t781, t1121 * t835 + t1126 * t801, t1121 * t782 + t1126 * t775, t1121 * t784 + t1126 * t777, t1121 * t770 + t1126 * t769, (t1121 * t793 - t1169 * t792) * qJ(3); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1087, 0, t1088, 0, t1214, -t1215, -t1093 * t1121 - t1094 * t1126, 0, t1027, -t1225, t1024, t1026, -t1223, t1060, t1068 * t1126 - t1081 * t1121, t1069 * t1126 - t1082 * t1121, t1121 * t1012, 0, t1027, t1024, t1225, t1060, t1223, t1026, t1022 * t1126 + t1121 * t991, -qJ(3) * t1155 + t1121 * t968, -t1008 * t1126 + t1121 * t1170, (-t1126 * t1033 + t1121 * t1129) * qJ(3), t1121 * t956 - t1126 * t988, t1121 * t916 - t1126 * t937, t1121 * t929 - t1126 * t980, t1121 * t955 + t1126 * t986, t1121 * t930 - t1126 * t981, -t1006 * t1126 + t1121 * t985, t1121 * t894 - t1126 * t944, t1121 * t903 - t1126 * t948, -t1126 * t863 + t1172 * t864, (-t1126 * t913 + t1150 * t1172) * qJ(3), t1121 * t849 - t1126 * t899, t1121 * t810 - t1126 * t850, t1121 * t833 - t1126 * t865, t1121 * t848 - t1126 * t898, t1121 * t834 - t1126 * t866, t1121 * t892 - t1126 * t921, t1121 * t812 + t1126 * t827, t1121 * t817 - t1126 * t829, t1121 * t795 - t1126 * t813, (-t1126 * t860 - t1172 * t859) * qJ(3), t1121 * t791 - t1126 * t804, t1121 * t773 - t1126 * t786, t1121 * t780 - t1126 * t796, t1121 * t790 - t1126 * t803, t1121 * t781 - t1126 * t797, t1121 * t801 - t1126 * t835, t1121 * t775 - t1126 * t782, t1121 * t777 - t1126 * t784, t1121 * t769 - t1126 * t770, (-t1126 * t793 - t1172 * t792) * qJ(3); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1093, t1094, 0, 0, t1052, t1036, t1055, t1051, t1053, 0, t1082, -t1081, t1011, 0, t1052, t1055, -t1036, 0, -t1053, t1051, t990, t967, t999, t1128, t954, t915, t927, t953, t928, t984, t893, t902, t862, t1164, t847, t809, t831, t846, t832, t891, t811, t816, t794, -t1165, t789, t772, t778, t788, t779, t800, t774, t776, t768, -t1166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1127, 0, 0, -g(3), -t1093, 0, t1061, -t1039, t1058, t1059, -t1057, t1077, -t1081, -t1082, t1012, 0, t1061, t1058, t1039, t1077, t1057, t1059, t991, t968, t1170, qJ(3) * t1129, t956, t916, t929, t955, t930, t985, t894, t903, t1125 * t864, t1150 * t1206, t849, t810, t833, t848, t834, t892, t812, t817, t795, -t859 * t1206, t791, t773, t780, t790, t781, t801, t775, t777, t769, -t792 * t1206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1127, 0, qJDD(1), 0, g(3), 0, -t1094, 0, t1100, -t1089, -t1168, -t1100, -t1167, -qJDD(2), t1068, t1069, 0, 0, t1100, -t1168, t1089, -qJDD(2), t1167, -t1100, t1022, -t1163, -t1008, -t1208, -t988, -t937, -t980, t986, -t981, -t1006, -t944, -t948, -t863, -t1211, -t899, -t850, -t865, -t898, -t866, -t921, t827, -t829, -t813, -t1212, -t804, -t786, -t796, -t803, -t797, -t835, -t782, -t784, -t770, -t1213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1093, t1094, 0, 0, t1052, t1036, t1055, t1051, t1053, 0, t1082, -t1081, t1011, 0, t1052, t1055, -t1036, 0, -t1053, t1051, t990, t967, t999, t1128, t954, t915, t927, t953, t928, t984, t893, t902, t862, t1164, t847, t809, t831, t846, t832, t891, t811, t816, t794, -t1165, t789, t772, t778, t788, t779, t800, t774, t776, t768, -t1166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1084, t1086, t1091, -t1106, t1097, t1106, 0, -t1093, t1068, 0, t1084, t1091, -t1086, t1106, -t1097, -t1106, qJ(3) * t1086, t1016, t1000, t1130, t989, t938, t982, t987, t983, t1007, t943, t949, t864, qJ(3) * t1150, t901, t851, t867, t900, t868, t922, t828, t830, t814, -qJ(3) * t859, t806, t787, t798, t805, t799, t836, t783, t785, t771, -qJ(3) * t792; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1105, t1083, -t1096, t1085, t1092, -t1105, t1093, 0, t1069, 0, t1105, -t1096, -t1083, -t1105, -t1092, t1085, t1218, t1033, 0, 0, -t1050, -t1044, -t997, t1050, t993, t1076, t978, t979, 0, 0, -t925, -t877, -t917, t923, -t918, -t959, t971, -t969, -t1142, 0, -t843, -t807, -t819, -t842, -t820, -t885, t855, -t856, -t838, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1100, t1089, t1168, t1100, t1167, qJDD(2), -t1068, -t1069, 0, 0, -t1100, t1168, -t1089, qJDD(2), -t1167, t1100, -t1022, t1163, t1008, t1208, t988, t937, t980, -t986, t981, t1006, t944, t948, t863, t1211, t899, t850, t865, t898, t866, t921, -t827, t829, t813, t1212, t804, t786, t796, t803, t797, t835, t782, t784, t770, t1213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1084, t1091, -t1086, t1106, -t1097, -t1106, 0, t1043, t1218, 0, t989, t938, t982, t987, t983, t1007, t1174, t1035, t1150, 0, t901, t851, t867, t900, t868, t922, t1138, t1151, t887, 0, t806, t787, t798, t805, t799, t836, t1210, t1209, t1154, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1100, t1168, -t1089, qJDD(2), -t1167, t1100, -t1043, 0, t1033, 0, t988, t937, t980, -t986, t981, t1006, -t1035, t1174, t913, 0, t899, t850, t865, t898, t866, t921, t1144, t1136, t1202, 0, t804, t786, t796, t803, t797, t835, t1153, t1152, t1137, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1105, t1096, t1083, t1105, t1092, -t1085, -t1218, -t1033, 0, 0, t1050, t1044, t997, -t1050, -t993, -t1076, -t978, -t979, 0, 0, t925, t877, t917, -t923, t918, t959, -t971, t969, t1142, 0, t843, t807, t819, t842, t820, t885, -t855, t856, t838, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1031, t994, t1219, t1065, t1063, -t1065, 0, t1043, t978, 0, t926, t878, t919, t924, t920, t960, t969, t971, t889, 0, t845, t808, t821, t844, t822, t886, t857, t858, t1123 * t1143, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1066, t996, t1064, -t1135, -t1017, -t1066, -t1043, 0, t979, 0, -t998, -t992, -t935, t998, t932, -t1023, t946, t947, 0, 0, -t875, -t823, -t869, t873, -t870, -t910, t1197, -t1204, -t839, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1050, t1044, t997, -t1050, -t993, -t1076, -t978, -t979, 0, 0, t925, t877, t917, -t923, t918, t959, -t971, t969, t1142, 0, t843, t807, t819, t842, t820, t885, -t855, t856, t838, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t964, t933, t957, t1013, t1009, -t1013, 0, t978, t946, 0, t876, t825, t871, t874, t872, t911, t1204, t1197, t1143, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1014, t934, t1010, -t1134, t958, -t1014, -t978, 0, t947, 0, -t952, -t951, -t884, t952, t880, -t962, t896, t897, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t998, t992, t935, -t998, -t932, t1023, -t946, -t947, 0, 0, t875, t823, t869, -t873, t870, t910, -t1197, t1204, t839, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t908, t881, t1221, t976, t974, -t976, 0, t946, t896, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t977, t883, t975, t907, t905, -t977, -t946, 0, t897, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t952, t951, t884, -t952, -t880, t962, -t896, -t897, 0, 0;];
m_new_reg  = t1;
