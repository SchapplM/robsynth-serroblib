% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPPPR2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:00:09
% EndTime: 2022-01-23 09:00:27
% DurationCPUTime: 19.56s
% Computational Cost: add. (65432->721), mult. (188231->992), div. (0->0), fcn. (132261->10), ass. (0->490)
t1147 = sin(pkin(9));
t1149 = sin(pkin(7));
t1139 = t1149 * qJDD(1);
t1148 = sin(pkin(8));
t1219 = t1148 * t1139;
t1150 = cos(pkin(9));
t1152 = cos(pkin(7));
t1151 = cos(pkin(8));
t1243 = t1149 * t1151;
t1189 = t1147 * t1243 + t1150 * t1152;
t1097 = t1189 * qJD(1);
t1227 = qJD(1) * t1243;
t1099 = -qJD(1) * t1147 * t1152 + t1150 * t1227;
t1251 = t1099 * t1097;
t1284 = t1219 - t1251;
t1286 = t1147 * t1284;
t1285 = t1150 * t1284;
t1153 = sin(qJ(5));
t1155 = cos(qJ(5));
t1245 = t1148 * t1149;
t1228 = qJD(1) * t1245;
t1063 = t1099 * t1153 - t1155 * t1228;
t1065 = t1099 * t1155 + t1153 * t1228;
t1007 = t1065 * t1063;
t1279 = t1189 * qJDD(1);
t1087 = qJDD(5) + t1279;
t1281 = -t1007 + t1087;
t1283 = t1153 * t1281;
t1282 = t1155 * t1281;
t1074 = t1099 * t1228;
t1036 = -t1279 + t1074;
t1157 = qJD(1) ^ 2;
t1237 = t1152 * t1157;
t1269 = pkin(2) * t1152;
t1199 = -qJ(3) * t1149 - t1269;
t1273 = 2 * qJD(2);
t1172 = (qJD(1) * t1199 + t1273) * qJD(1);
t1154 = sin(qJ(1));
t1156 = cos(qJ(1));
t1128 = t1156 * g(1) + t1154 * g(2);
t1174 = -t1157 * pkin(1) + qJDD(1) * qJ(2) - t1128;
t1280 = t1172 + t1174;
t1267 = pkin(3) * t1151;
t1198 = -qJ(4) * t1148 - t1267;
t1268 = pkin(3) * t1148;
t1197 = -qJ(4) * t1151 + t1268;
t1264 = t1152 * g(3);
t1213 = qJDD(3) + t1264;
t1161 = ((t1152 * t1198 - pkin(1)) * t1157 + (qJ(2) + t1197) * qJDD(1) + t1172 - t1128) * t1149 + t1213;
t1254 = qJD(1) * t1149;
t1100 = t1197 * t1254;
t1141 = t1152 * qJDD(1);
t1159 = t1152 ^ 2;
t1142 = t1159 * t1157;
t1265 = t1149 * g(3);
t1032 = t1152 * t1280 - t1265;
t1127 = t1154 * g(1) - g(2) * t1156;
t1184 = -t1157 * qJ(2) + qJDD(2) - t1127;
t1191 = -pkin(1) + t1199;
t1162 = qJDD(1) * t1191 + t1184;
t979 = -0.2e1 * qJD(3) * t1228 + t1032 * t1151 + t1148 * t1162;
t951 = -pkin(3) * t1142 - qJ(4) * t1141 - t1100 * t1228 + t979;
t1212 = t1147 * t951 - t1150 * t1161;
t1271 = 2 * qJD(4);
t898 = t1099 * t1271 + t1212;
t899 = -t1097 * t1271 + t1147 * t1161 + t1150 * t951;
t843 = t1147 * t898 + t1150 * t899;
t1210 = t1148 * t1032 - t1151 * t1162;
t1173 = pkin(3) * t1141 - qJ(4) * t1142 + qJDD(4) + t1210;
t1272 = 0.2e1 * qJD(3);
t1211 = (t1272 + t1100) * t1151;
t948 = t1211 * t1254 + t1173;
t826 = t1148 * t948 + t1151 * t843;
t842 = t1147 * t899 - t1150 * t898;
t1278 = -(pkin(2) - t1198) * t842 + qJ(3) * t826;
t1229 = pkin(4) * t1150 + pkin(3);
t1045 = pkin(4) * t1097 - pkin(6) * t1099;
t1145 = t1149 ^ 2;
t1244 = t1148 * t1157;
t872 = -t1097 * t1045 + (-pkin(4) * t1145 * t1244 + pkin(6) * t1139) * t1148 + t899;
t1217 = t1151 * t1139;
t1092 = -t1141 * t1147 + t1150 * t1217;
t903 = t1279 * pkin(4) - t1092 * pkin(6) + (t1211 + (pkin(4) * t1099 + pkin(6) * t1097) * t1148) * t1254 + t1173;
t833 = t1153 * t872 - t1155 * t903;
t834 = t1153 * t903 + t1155 * t872;
t800 = t1153 * t833 + t1155 * t834;
t1250 = t1145 * t1157;
t1274 = t1148 ^ 2;
t1134 = t1274 * t1250;
t871 = -pkin(4) * t1219 - pkin(6) * t1134 + (t1271 + t1045) * t1099 + t1212;
t792 = t1147 * t871 + t1150 * t800;
t799 = t1153 * t834 - t1155 * t833;
t1277 = -(pkin(6) * t1147 + t1229) * t799 + qJ(4) * t792;
t1088 = qJD(5) + t1097;
t1209 = -t1153 * t1092 + t1155 * t1219;
t958 = (-qJD(5) + t1088) * t1065 + t1209;
t1179 = -t1155 * t1092 - t1153 * t1219;
t1002 = -qJD(5) * t1063 - t1179;
t1021 = t1088 * t1063;
t960 = t1002 + t1021;
t892 = t1153 * t958 - t1155 * t960;
t793 = -pkin(6) * t892 - t799;
t894 = t1153 * t960 + t1155 * t958;
t1060 = t1063 ^ 2;
t1061 = t1065 ^ 2;
t966 = t1060 + t1061;
t859 = -t1147 * t966 + t1150 * t894;
t1276 = qJ(4) * t859 + t1147 * t793 - t1229 * t892;
t1206 = t1097 * t1228;
t1038 = t1206 + t1092;
t968 = t1036 * t1147 - t1038 * t1150;
t829 = -qJ(4) * t968 - t842;
t1089 = t1097 ^ 2;
t1090 = t1099 ^ 2;
t1018 = t1089 + t1090;
t970 = t1036 * t1150 + t1038 * t1147;
t924 = -t1018 * t1148 + t1151 * t970;
t1275 = qJ(3) * t924 + t1148 * t829 - (pkin(2) + t1267) * t968;
t1086 = t1088 ^ 2;
t1270 = pkin(2) * t1149;
t1266 = pkin(4) * t1147;
t1261 = qJDD(1) * pkin(1);
t1260 = t1147 * t948;
t1259 = t1150 * t948;
t1258 = t1153 * t871;
t985 = t1007 + t1087;
t1257 = t1153 * t985;
t1256 = t1155 * t871;
t1255 = t1155 * t985;
t1253 = t1088 * t1153;
t1252 = t1088 * t1155;
t1030 = t1251 + t1219;
t1249 = t1147 * t1030;
t1028 = t1149 * t1280 + t1213;
t1248 = t1148 * t1028;
t1223 = t1151 * t1244;
t1116 = t1145 * t1223;
t1102 = t1141 - t1116;
t1247 = t1148 * t1102;
t1103 = -t1141 - t1116;
t1246 = t1148 * t1103;
t1242 = t1149 * t1152;
t1241 = t1150 * t1030;
t1240 = t1151 * t1028;
t1239 = t1151 * t1102;
t1238 = t1151 * t1103;
t1101 = -t1184 + t1261;
t1236 = t1154 * t1101;
t1235 = t1156 * t1101;
t1234 = qJDD(1) * t1148;
t1233 = qJDD(1) * t1151;
t1232 = t1154 * qJDD(1);
t1231 = t1156 * qJDD(1);
t1146 = t1151 ^ 2;
t1226 = t1146 * t1250;
t1225 = t1147 * t1007;
t1224 = t1148 * t1251;
t1222 = t1150 * t1007;
t1221 = t1151 * t1251;
t1220 = t1151 * t1237;
t1135 = t1149 * t1237;
t1218 = t1149 * t1141;
t1216 = t1152 * t1232;
t1215 = t1152 * t1231;
t1214 = t1101 + t1261;
t978 = t1227 * t1272 + t1210;
t916 = t1148 * t978 + t1151 * t979;
t1163 = qJD(1) * t1273 + t1174;
t1068 = t1149 * t1163 + t1264;
t1069 = t1152 * t1163 - t1265;
t1006 = t1068 * t1149 + t1069 * t1152;
t1208 = -t1127 * t1154 - t1128 * t1156;
t1144 = t1149 * t1145;
t1207 = t1144 * t1223;
t1205 = t1148 * t1220;
t1204 = t1151 * t1135;
t1203 = -pkin(4) * t871 + pkin(6) * t800;
t1202 = -pkin(3) * t948 + qJ(4) * t843;
t1124 = -t1154 * t1157 + t1231;
t1201 = -pkin(5) * t1124 - g(3) * t1154;
t1200 = -pkin(2) * t1028 + qJ(3) * t916;
t915 = t1148 * t979 - t1151 * t978;
t1193 = t1145 * t1205;
t1005 = t1068 * t1152 - t1069 * t1149;
t1192 = t1127 * t1156 - t1128 * t1154;
t1123 = t1156 * t1157 + t1232;
t1001 = -qJD(5) * t1065 + t1209;
t1112 = (t1145 + t1159) * t1237;
t1188 = -t1112 * t1154 + t1215;
t1187 = t1112 * t1156 + t1216;
t994 = -t1061 - t1086;
t922 = -t1153 * t994 - t1255;
t961 = (qJD(5) + t1088) * t1063 + t1179;
t1186 = pkin(4) * t961 + pkin(6) * t922 + t1258;
t983 = -t1086 - t1060;
t918 = t1155 * t983 - t1283;
t1022 = t1088 * t1065;
t957 = t1001 - t1022;
t1185 = pkin(4) * t957 + pkin(6) * t918 - t1256;
t1035 = -t1074 - t1279;
t1033 = -t1134 - t1089;
t965 = t1033 * t1150 - t1286;
t1183 = pkin(3) * t1035 + qJ(4) * t965 - t1259;
t1037 = -t1206 + t1092;
t1072 = -t1090 - t1134;
t992 = -t1072 * t1147 - t1241;
t1182 = -pkin(3) * t1037 + qJ(4) * t992 + t1260;
t1108 = -t1134 - t1142;
t1055 = t1108 * t1151 - t1246;
t1094 = (t1220 - t1234) * t1149;
t1181 = pkin(2) * t1094 + qJ(3) * t1055 - t1240;
t1111 = -t1142 - t1226;
t1056 = -t1111 * t1148 + t1239;
t1118 = t1148 * t1135;
t1095 = t1118 + t1217;
t1180 = -pkin(2) * t1095 + qJ(3) * t1056 + t1248;
t1178 = pkin(4) * t966 + pkin(6) * t894 + t800;
t1176 = pkin(3) * t1018 + qJ(4) * t970 + t843;
t1093 = (t1220 + t1234) * t1149;
t1096 = -t1118 + t1217;
t1042 = -t1093 * t1151 + t1096 * t1148;
t1104 = t1134 + t1226;
t1175 = pkin(2) * t1104 + qJ(3) * t1042 + t916;
t791 = t1147 * t800 - t1150 * t871;
t771 = -qJ(4) * t791 + (-pkin(6) * t1150 + t1266) * t799;
t775 = -pkin(3) * t791 - t1203;
t778 = t1148 * t799 + t1151 * t792;
t1171 = -pkin(2) * t791 + qJ(3) * t778 + t1148 * t771 + t1151 * t775;
t858 = t1147 * t894 + t1150 * t966;
t782 = -qJ(4) * t858 + t1150 * t793 + t1266 * t892;
t784 = -pkin(3) * t858 - t1178;
t828 = t1148 * t892 + t1151 * t859;
t1170 = -pkin(2) * t858 + qJ(3) * t828 + t1148 * t782 + t1151 * t784;
t917 = t1153 * t983 + t1282;
t817 = -pkin(4) * t917 + t833;
t839 = -pkin(6) * t917 + t1258;
t873 = t1147 * t918 + t1150 * t957;
t787 = -qJ(4) * t873 - t1147 * t817 + t1150 * t839;
t805 = -pkin(3) * t873 - t1185;
t874 = -t1147 * t957 + t1150 * t918;
t841 = t1148 * t917 + t1151 * t874;
t1169 = -pkin(2) * t873 + qJ(3) * t841 + t1148 * t787 + t1151 * t805;
t921 = t1155 * t994 - t1257;
t818 = -pkin(4) * t921 + t834;
t844 = -pkin(6) * t921 + t1256;
t875 = t1147 * t922 + t1150 * t961;
t788 = -qJ(4) * t875 - t1147 * t818 + t1150 * t844;
t806 = -pkin(3) * t875 - t1186;
t876 = -t1147 * t961 + t1150 * t922;
t847 = t1148 * t921 + t1151 * t876;
t1168 = -pkin(2) * t875 + qJ(3) * t847 + t1148 * t788 + t1151 * t806;
t964 = t1033 * t1147 + t1285;
t860 = -pkin(3) * t964 + t898;
t901 = -qJ(4) * t964 + t1260;
t930 = -t1035 * t1148 + t1151 * t965;
t1167 = -pkin(2) * t964 + qJ(3) * t930 + t1148 * t901 + t1151 * t860;
t989 = t1072 * t1150 - t1249;
t867 = -pkin(3) * t989 + t899;
t910 = -qJ(4) * t989 + t1259;
t945 = t1037 * t1148 + t1151 * t992;
t1166 = -pkin(2) * t989 + qJ(3) * t945 + t1148 * t910 + t1151 * t867;
t1165 = -pkin(3) * t917 + qJ(4) * t874 + t1147 * t839 + t1150 * t817;
t1164 = -pkin(3) * t921 + qJ(4) * t876 + t1147 * t844 + t1150 * t818;
t1140 = t1159 * qJDD(1);
t1138 = t1145 * qJDD(1);
t1133 = t1274 * t1139;
t1130 = 0.2e1 * t1218;
t1120 = -t1142 + t1250;
t1119 = t1142 + t1250;
t1115 = t1140 - t1138;
t1114 = t1140 + t1138;
t1110 = t1142 - t1226;
t1109 = (t1149 * t1159 + t1144) * t1157;
t1107 = t1134 - t1142;
t1106 = -pkin(5) * t1123 + g(3) * t1156;
t1105 = -t1134 + t1226;
t1085 = t1124 * t1242;
t1084 = t1123 * t1242;
t1081 = (t1146 + t1274) * t1135;
t1080 = (qJDD(1) * t1146 + t1205) * t1149;
t1079 = (-t1146 * t1237 + t1148 * t1233) * t1149;
t1078 = (t1148 * t1237 + t1233) * t1245;
t1077 = -t1148 * t1204 + t1133;
t1076 = t1109 * t1156 + t1149 * t1232;
t1075 = t1109 * t1154 - t1149 * t1231;
t1071 = -t1090 + t1134;
t1070 = t1089 - t1134;
t1059 = t1080 * t1152 + t1207;
t1058 = t1077 * t1152 - t1207;
t1057 = -t1110 * t1148 + t1238;
t1054 = t1107 * t1151 + t1247;
t1053 = t1110 * t1151 + t1246;
t1052 = t1111 * t1151 + t1247;
t1051 = t1108 * t1148 + t1238;
t1050 = t1107 * t1148 - t1239;
t1048 = t1080 * t1149 - t1193;
t1047 = t1077 * t1149 + t1193;
t1046 = t1090 - t1089;
t1044 = -qJ(2) * t1112 + t1152 * t1214;
t1043 = qJ(2) * t1109 - t1149 * t1214;
t1041 = t1094 * t1151 - t1095 * t1148;
t1040 = -t1093 * t1148 - t1096 * t1151;
t1039 = t1094 * t1148 + t1095 * t1151;
t1026 = -t1074 * t1147 + t1092 * t1150;
t1025 = t1074 * t1150 + t1092 * t1147;
t1024 = t1147 * t1279 + t1150 * t1206;
t1023 = -t1147 * t1206 + t1150 * t1279;
t1020 = -t1061 + t1086;
t1019 = t1060 - t1086;
t1017 = (-t1097 * t1150 + t1099 * t1147) * t1228;
t1016 = (-t1097 * t1147 - t1099 * t1150) * t1228;
t1015 = t1057 * t1152 + t1096 * t1149;
t1014 = t1056 * t1152 + t1095 * t1149;
t1013 = t1055 * t1152 - t1094 * t1149;
t1012 = t1054 * t1152 - t1093 * t1149;
t1011 = t1057 * t1149 - t1096 * t1152;
t1010 = t1056 * t1149 - t1095 * t1152;
t1009 = t1055 * t1149 + t1094 * t1152;
t1008 = t1054 * t1149 + t1093 * t1152;
t1003 = t1061 - t1060;
t1000 = t1042 * t1152 - t1104 * t1149;
t999 = t1041 * t1152 + t1105 * t1149;
t998 = t1042 * t1149 + t1104 * t1152;
t997 = t1041 * t1149 - t1105 * t1152;
t996 = t1017 * t1151 + t1133;
t995 = (t1017 - t1217) * t1148;
t993 = t1070 * t1150 - t1249;
t991 = -t1071 * t1147 + t1285;
t990 = t1070 * t1147 + t1241;
t988 = t1071 * t1150 + t1286;
t987 = pkin(1) * t1101 + qJ(2) * t1006;
t982 = pkin(1) * t1119 + qJ(2) * t1114 + t1006;
t981 = -qJ(3) * t1052 + t1240;
t980 = -qJ(3) * t1051 + t1248;
t974 = t1026 * t1151 + t1224;
t973 = t1024 * t1151 - t1224;
t972 = t1026 * t1148 - t1221;
t971 = t1024 * t1148 + t1221;
t969 = t1035 * t1150 - t1037 * t1147;
t967 = t1035 * t1147 + t1037 * t1150;
t963 = (-t1063 * t1155 + t1065 * t1153) * t1088;
t962 = (-t1063 * t1153 - t1065 * t1155) * t1088;
t959 = t1002 - t1021;
t956 = -t1001 - t1022;
t955 = t1002 * t1155 - t1065 * t1253;
t954 = t1002 * t1153 + t1065 * t1252;
t953 = -t1001 * t1153 + t1063 * t1252;
t952 = -t1001 * t1155 - t1063 * t1253;
t950 = -pkin(2) * t1052 + t979;
t949 = -pkin(2) * t1051 + t978;
t946 = t1036 * t1148 + t1151 * t993;
t944 = t1038 * t1148 + t1151 * t991;
t943 = -t1036 * t1151 + t1148 * t993;
t942 = -t1037 * t1151 + t1148 * t992;
t941 = -t1038 * t1151 + t1148 * t991;
t940 = t1016 * t1149 + t1152 * t996;
t939 = -t1016 * t1152 + t1149 * t996;
t938 = t1087 * t1147 + t1150 * t963;
t937 = -t1087 * t1150 + t1147 * t963;
t936 = t1046 * t1148 + t1151 * t969;
t935 = -t1046 * t1151 + t1148 * t969;
t934 = t1019 * t1155 - t1257;
t933 = -t1020 * t1153 + t1282;
t932 = t1019 * t1153 + t1255;
t931 = t1020 * t1155 + t1283;
t929 = t1035 * t1151 + t1148 * t965;
t928 = t1025 * t1149 + t1152 * t974;
t927 = -t1023 * t1149 + t1152 * t973;
t926 = -t1025 * t1152 + t1149 * t974;
t925 = t1023 * t1152 + t1149 * t973;
t923 = t1018 * t1151 + t1148 * t970;
t920 = -pkin(1) * t1009 - t1181;
t919 = -pkin(1) * t1010 - t1180;
t914 = t1150 * t955 + t1225;
t913 = t1150 * t953 - t1225;
t912 = t1147 * t955 - t1222;
t911 = t1147 * t953 + t1222;
t909 = t1149 * t990 + t1152 * t946;
t908 = t1149 * t989 + t1152 * t945;
t907 = t1149 * t988 + t1152 * t944;
t906 = t1149 * t946 - t1152 * t990;
t905 = t1149 * t945 - t1152 * t989;
t904 = t1149 * t944 - t1152 * t988;
t900 = -qJ(3) * t1040 - t915;
t897 = t1028 * t1149 + t1152 * t916;
t896 = -t1028 * t1152 + t1149 * t916;
t893 = -t1153 * t959 + t1155 * t957;
t891 = t1153 * t957 + t1155 * t959;
t890 = t1149 * t967 + t1152 * t936;
t889 = t1149 * t936 - t1152 * t967;
t888 = t1148 * t962 + t1151 * t938;
t887 = t1148 * t938 - t1151 * t962;
t886 = t1149 * t964 + t1152 * t930;
t885 = t1149 * t930 - t1152 * t964;
t884 = t1149 * t968 + t1152 * t924;
t883 = t1149 * t924 - t1152 * t968;
t882 = -t1147 * t956 + t1150 * t934;
t881 = t1147 * t960 + t1150 * t933;
t880 = t1147 * t934 + t1150 * t956;
t879 = t1147 * t933 - t1150 * t960;
t878 = -qJ(2) * t1010 - t1149 * t950 + t1152 * t981;
t877 = -qJ(2) * t1009 - t1149 * t949 + t1152 * t980;
t870 = t1003 * t1147 + t1150 * t893;
t869 = -t1003 * t1150 + t1147 * t893;
t868 = -pkin(1) * t998 - t1175;
t866 = -pkin(1) * t1052 + qJ(2) * t1014 + t1149 * t981 + t1152 * t950;
t865 = -pkin(1) * t1051 + qJ(2) * t1013 + t1149 * t980 + t1152 * t949;
t864 = t1148 * t954 + t1151 * t914;
t863 = -t1148 * t952 + t1151 * t913;
t862 = t1148 * t914 - t1151 * t954;
t861 = t1148 * t913 + t1151 * t952;
t857 = -qJ(2) * t998 + t1040 * t1270 + t1152 * t900;
t856 = -pkin(2) * t942 - t1182;
t855 = t1149 * t937 + t1152 * t888;
t854 = t1149 * t888 - t1152 * t937;
t853 = qJ(2) * t1000 + t1149 * t900 + (-pkin(1) - t1269) * t1040;
t852 = t1148 * t932 + t1151 * t882;
t851 = t1148 * t931 + t1151 * t881;
t850 = t1148 * t882 - t1151 * t932;
t849 = t1148 * t881 - t1151 * t931;
t848 = -pkin(2) * t929 - t1183;
t846 = t1148 * t876 - t1151 * t921;
t845 = -pkin(1) * t896 - t1200;
t840 = t1148 * t874 - t1151 * t917;
t838 = t1149 * t912 + t1152 * t864;
t837 = t1149 * t911 + t1152 * t863;
t836 = t1149 * t864 - t1152 * t912;
t835 = t1149 * t863 - t1152 * t911;
t831 = t1148 * t891 + t1151 * t870;
t830 = t1148 * t870 - t1151 * t891;
t827 = t1148 * t859 - t1151 * t892;
t825 = t1148 * t843 - t1151 * t948;
t824 = -qJ(2) * t896 + (-qJ(3) * t1152 + t1270) * t915;
t823 = -qJ(3) * t942 - t1148 * t867 + t1151 * t910;
t822 = t1149 * t880 + t1152 * t852;
t821 = t1149 * t879 + t1152 * t851;
t820 = t1149 * t852 - t1152 * t880;
t819 = t1149 * t851 - t1152 * t879;
t816 = -qJ(3) * t929 - t1148 * t860 + t1151 * t901;
t815 = t1149 * t875 + t1152 * t847;
t814 = t1149 * t847 - t1152 * t875;
t813 = t1149 * t873 + t1152 * t841;
t812 = t1149 * t841 - t1152 * t873;
t811 = -pkin(2) * t923 - t1176;
t810 = qJ(2) * t897 + t1191 * t915;
t809 = -qJ(3) * t923 + t1151 * t829 + t1268 * t968;
t808 = t1149 * t869 + t1152 * t831;
t807 = t1149 * t831 - t1152 * t869;
t804 = t1149 * t858 + t1152 * t828;
t803 = t1149 * t828 - t1152 * t858;
t802 = -pkin(1) * t905 - t1166;
t801 = -pkin(1) * t885 - t1167;
t798 = t1149 * t842 + t1152 * t826;
t797 = t1149 * t826 - t1152 * t842;
t796 = -pkin(1) * t883 - t1275;
t795 = -pkin(2) * t825 - t1202;
t794 = -qJ(2) * t905 - t1149 * t856 + t1152 * t823;
t790 = -qJ(2) * t885 - t1149 * t848 + t1152 * t816;
t789 = -pkin(1) * t942 + qJ(2) * t908 + t1149 * t823 + t1152 * t856;
t786 = -pkin(1) * t929 + qJ(2) * t886 + t1149 * t816 + t1152 * t848;
t785 = -qJ(3) * t825 + t1197 * t842;
t783 = -qJ(2) * t883 - t1149 * t811 + t1152 * t809;
t781 = -pkin(1) * t923 + qJ(2) * t884 + t1149 * t809 + t1152 * t811;
t780 = -pkin(2) * t846 - t1164;
t779 = -pkin(2) * t840 - t1165;
t777 = t1148 * t792 - t1151 * t799;
t776 = -pkin(2) * t827 - t1276;
t774 = -qJ(3) * t846 - t1148 * t806 + t1151 * t788;
t773 = -qJ(3) * t840 - t1148 * t805 + t1151 * t787;
t772 = -pkin(1) * t797 - t1278;
t770 = -pkin(1) * t814 - t1168;
t769 = t1149 * t791 + t1152 * t778;
t768 = t1149 * t778 - t1152 * t791;
t767 = -qJ(2) * t797 - t1149 * t795 + t1152 * t785;
t766 = -pkin(1) * t812 - t1169;
t765 = -qJ(3) * t827 - t1148 * t784 + t1151 * t782;
t764 = -pkin(1) * t825 + qJ(2) * t798 + t1149 * t785 + t1152 * t795;
t763 = -pkin(1) * t803 - t1170;
t762 = -qJ(2) * t814 - t1149 * t780 + t1152 * t774;
t761 = -qJ(2) * t812 - t1149 * t779 + t1152 * t773;
t760 = -pkin(1) * t846 + qJ(2) * t815 + t1149 * t774 + t1152 * t780;
t759 = -pkin(1) * t840 + qJ(2) * t813 + t1149 * t773 + t1152 * t779;
t758 = -pkin(2) * t777 - t1277;
t757 = -qJ(2) * t803 - t1149 * t776 + t1152 * t765;
t756 = -pkin(1) * t827 + qJ(2) * t804 + t1149 * t765 + t1152 * t776;
t755 = -qJ(3) * t777 - t1148 * t775 + t1151 * t771;
t754 = -pkin(1) * t768 - t1171;
t753 = -qJ(2) * t768 - t1149 * t758 + t1152 * t755;
t752 = -pkin(1) * t777 + qJ(2) * t769 + t1149 * t755 + t1152 * t758;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1124, 0, -t1123, 0, t1201, -t1106, -t1192, -pkin(5) * t1192, t1085, t1115 * t1156 + t1120 * t1154, t1076, -t1085, t1187, 0, -pkin(5) * t1188 - t1154 * t1068 - t1149 * t1235, -pkin(5) * t1075 - t1069 * t1154 - t1152 * t1235, t1156 * t1005 - pkin(5) * (t1114 * t1154 + t1119 * t1156), -pkin(5) * (t1006 * t1154 + t1235) - (pkin(1) * t1154 - qJ(2) * t1156) * t1005, t1059 * t1156 + t1079 * t1154, t1039 * t1154 + t1156 * t999, t1015 * t1156 + t1053 * t1154, t1058 * t1156 - t1078 * t1154, t1012 * t1156 + t1050 * t1154, t1081 * t1154 - t1149 * t1215, t1156 * t877 - t1154 * t920 - pkin(5) * (t1013 * t1154 - t1051 * t1156), t1156 * t878 - t1154 * t919 - pkin(5) * (t1014 * t1154 - t1052 * t1156), t1156 * t857 - t1154 * t868 - pkin(5) * (t1000 * t1154 - t1040 * t1156), t1156 * t824 - t1154 * t845 - pkin(5) * (t1154 * t897 - t1156 * t915), t1154 * t972 + t1156 * t928, t1154 * t935 + t1156 * t890, t1154 * t941 + t1156 * t907, t1154 * t971 + t1156 * t927, t1154 * t943 + t1156 * t909, t1154 * t995 + t1156 * t940, t1156 * t790 - t1154 * t801 - pkin(5) * (t1154 * t886 - t1156 * t929), t1156 * t794 - t1154 * t802 - pkin(5) * (t1154 * t908 - t1156 * t942), t1156 * t783 - t1154 * t796 - pkin(5) * (t1154 * t884 - t1156 * t923), t1156 * t767 - t1154 * t772 - pkin(5) * (t1154 * t798 - t1156 * t825), t1154 * t862 + t1156 * t838, t1154 * t830 + t1156 * t808, t1154 * t849 + t1156 * t821, t1154 * t861 + t1156 * t837, t1154 * t850 + t1156 * t822, t1154 * t887 + t1156 * t855, t1156 * t761 - t1154 * t766 - pkin(5) * (t1154 * t813 - t1156 * t840), t1156 * t762 - t1154 * t770 - pkin(5) * (t1154 * t815 - t1156 * t846), t1156 * t757 - t1154 * t763 - pkin(5) * (t1154 * t804 - t1156 * t827), t1156 * t753 - t1154 * t754 - pkin(5) * (t1154 * t769 - t1156 * t777); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1123, 0, t1124, 0, t1106, t1201, t1208, pkin(5) * t1208, t1084, t1115 * t1154 - t1120 * t1156, t1075, -t1084, -t1188, 0, -pkin(5) * t1187 + t1156 * t1068 - t1149 * t1236, pkin(5) * t1076 + t1069 * t1156 - t1152 * t1236, t1154 * t1005 + pkin(5) * (t1114 * t1156 - t1119 * t1154), pkin(5) * (t1006 * t1156 - t1236) - (-pkin(1) * t1156 - qJ(2) * t1154) * t1005, t1059 * t1154 - t1079 * t1156, -t1039 * t1156 + t1154 * t999, t1015 * t1154 - t1053 * t1156, t1058 * t1154 + t1078 * t1156, t1012 * t1154 - t1050 * t1156, -t1081 * t1156 - t1149 * t1216, t1154 * t877 + t1156 * t920 + pkin(5) * (t1013 * t1156 + t1051 * t1154), t1154 * t878 + t1156 * t919 + pkin(5) * (t1014 * t1156 + t1052 * t1154), t1154 * t857 + t1156 * t868 + pkin(5) * (t1000 * t1156 + t1040 * t1154), t1154 * t824 + t1156 * t845 + pkin(5) * (t1154 * t915 + t1156 * t897), t1154 * t928 - t1156 * t972, t1154 * t890 - t1156 * t935, t1154 * t907 - t1156 * t941, t1154 * t927 - t1156 * t971, t1154 * t909 - t1156 * t943, t1154 * t940 - t1156 * t995, t1154 * t790 + t1156 * t801 + pkin(5) * (t1154 * t929 + t1156 * t886), t1154 * t794 + t1156 * t802 + pkin(5) * (t1154 * t942 + t1156 * t908), t1154 * t783 + t1156 * t796 + pkin(5) * (t1154 * t923 + t1156 * t884), t1154 * t767 + t1156 * t772 + pkin(5) * (t1154 * t825 + t1156 * t798), t1154 * t838 - t1156 * t862, t1154 * t808 - t1156 * t830, t1154 * t821 - t1156 * t849, t1154 * t837 - t1156 * t861, t1154 * t822 - t1156 * t850, t1154 * t855 - t1156 * t887, t1154 * t761 + t1156 * t766 + pkin(5) * (t1154 * t840 + t1156 * t813), t1154 * t762 + t1156 * t770 + pkin(5) * (t1154 * t846 + t1156 * t815), t1154 * t757 + t1156 * t763 + pkin(5) * (t1154 * t827 + t1156 * t804), t1154 * t753 + t1156 * t754 + pkin(5) * (t1154 * t777 + t1156 * t769); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1127, t1128, 0, 0, t1138, t1130, 0, t1140, 0, 0, t1044, t1043, t982, t987, t1048, t997, t1011, t1047, t1008, t1140, t865, t866, t853, t810, t926, t889, t904, t925, t906, t939, t786, t789, t781, t764, t836, t807, t819, t835, t820, t854, t759, t760, t756, t752; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1157, 0, 0, -g(3), -t1127, 0, t1218, t1115, t1109, -t1218, t1112, 0, -t1149 * t1101, -t1152 * t1101, t1005, qJ(2) * t1005, t1059, t999, t1015, t1058, t1012, -t1218, t877, t878, t857, t824, t928, t890, t907, t927, t909, t940, t790, t794, t783, t767, t838, t808, t821, t837, t822, t855, t761, t762, t757, t753; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1157, 0, qJDD(1), 0, g(3), 0, -t1128, 0, t1135, -t1120, -t1139, -t1135, -t1141, 0, t1068, t1069, 0, pkin(1) * t1005, -t1079, -t1039, -t1053, t1078, -t1050, -t1081, t920, t919, t868, t845, -t972, -t935, -t941, -t971, -t943, -t995, t801, t802, t796, t772, -t862, -t830, -t849, -t861, -t850, -t887, t766, t770, t763, t754; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1127, t1128, 0, 0, t1138, t1130, 0, t1140, 0, 0, t1044, t1043, t982, t987, t1048, t997, t1011, t1047, t1008, t1140, t865, t866, t853, t810, t926, t889, t904, t925, t906, t939, t786, t789, t781, t764, t836, t807, t819, t835, t820, t854, t759, t760, t756, t752; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1139, t1141, t1135, 0, t1142, 0, 0, -t1101, t1068, 0, t1080, t1041, t1057, t1077, t1054, 0, t980, t981, t900, -qJ(3) * t915, t974, t936, t944, t973, t946, t996, t816, t823, t809, t785, t864, t831, t851, t863, t852, t888, t773, t774, t765, t755; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1139, -t1250, t1141, -t1135, 0, t1101, 0, t1069, 0, -t1116, -t1105, -t1096, t1116, t1093, t1141, t949, t950, -pkin(2) * t1040, -pkin(2) * t915, -t1025, -t967, -t988, t1023, -t990, -t1016, t848, t856, t811, t795, -t912, -t869, -t879, -t911, -t880, -t937, t779, t780, t776, t758; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1135, t1120, t1139, t1135, t1141, 0, -t1068, -t1069, 0, 0, t1079, t1039, t1053, -t1078, t1050, t1081, t1181, t1180, t1175, t1200, t972, t935, t941, t971, t943, t995, t1167, t1166, t1275, t1278, t862, t830, t849, t861, t850, t887, t1169, t1168, t1170, t1171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1217, t1094, t1103, -t1118, t1107, t1118, 0, t1028, t978, 0, t1026, t969, t991, t1024, t993, t1017, t901, t910, t829, -qJ(4) * t842, t914, t870, t881, t913, t882, t938, t787, t788, t782, t771; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1204, t1095, t1110, -t1219, -t1102, t1204, -t1028, 0, t979, 0, -t1251, -t1046, -t1038, t1251, -t1036, -t1219, t860, t867, -pkin(3) * t968, -pkin(3) * t842, -t954, -t891, -t931, t952, -t932, -t962, t805, t806, t784, t775; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1116, t1105, t1096, -t1116, -t1093, -t1141, -t978, -t979, 0, 0, t1025, t967, t988, -t1023, t990, t1016, t1183, t1182, t1176, t1202, t912, t869, t879, t911, t880, t937, t1165, t1164, t1276, t1277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1092, t1035, t1284, t1206, t1070, -t1206, 0, t948, t898, 0, t955, t893, t933, t953, t934, t963, t839, t844, t793, -pkin(6) * t799; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1074, t1037, t1071, -t1279, t1030, -t1074, -t948, 0, t899, 0, -t1007, -t1003, -t960, t1007, t956, -t1087, t817, t818, -pkin(4) * t892, -pkin(4) * t799; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1251, t1046, t1038, -t1251, t1036, t1219, -t898, -t899, 0, 0, t954, t891, t931, -t952, t932, t962, t1185, t1186, t1178, t1203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1002, t957, t1281, t1021, t1019, -t1021, 0, t871, t833, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1022, t959, t1020, t1001, t985, -t1022, -t871, 0, t834, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1007, t1003, t960, -t1007, -t956, t1087, -t833, -t834, 0, 0;];
m_new_reg = t1;
