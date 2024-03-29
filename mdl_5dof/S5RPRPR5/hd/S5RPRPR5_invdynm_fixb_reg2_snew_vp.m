% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPRPR5_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:26:24
% EndTime: 2022-01-23 09:26:42
% DurationCPUTime: 18.64s
% Computational Cost: add. (115607->713), mult. (302492->993), div. (0->0), fcn. (213627->10), ass. (0->490)
t1196 = cos(pkin(8));
t1182 = t1196 * qJD(1) - qJD(3);
t1308 = t1182 ^ 2;
t1193 = sin(pkin(9));
t1195 = cos(pkin(9));
t1198 = sin(qJ(3));
t1201 = cos(qJ(3));
t1194 = sin(pkin(8));
t1291 = qJD(1) * t1194;
t1137 = (t1193 * t1201 + t1195 * t1198) * t1291;
t1309 = t1137 ^ 2;
t1077 = -t1308 - t1309;
t1290 = qJD(1) * t1201;
t1251 = t1194 * t1290;
t1276 = t1194 * t1198;
t1252 = qJD(1) * t1276;
t1138 = -t1193 * t1252 + t1195 * t1251;
t1094 = t1138 * t1137;
t1258 = t1196 * qJDD(1);
t1181 = -qJDD(3) + t1258;
t1311 = -t1094 - t1181;
t1316 = t1195 * t1311;
t1014 = t1193 * t1077 + t1316;
t1199 = sin(qJ(1));
t1202 = cos(qJ(1));
t1172 = t1202 * g(1) + t1199 * g(2);
t1203 = qJD(1) ^ 2;
t1152 = -t1203 * pkin(1) + qJDD(1) * qJ(2) - t1172;
t1302 = pkin(2) * t1196;
t1232 = -pkin(6) * t1194 - t1302;
t1306 = 2 * qJD(2);
t1262 = t1232 * qJD(1) + t1306;
t1227 = t1262 * qJD(1) + t1152;
t1300 = t1194 * g(3);
t1088 = t1196 * t1227 - t1300;
t1171 = t1199 * g(1) - t1202 * g(2);
t1220 = -t1203 * qJ(2) + qJDD(2) - t1171;
t1229 = -pkin(1) + t1232;
t1125 = qJDD(1) * t1229 + t1220;
t1112 = t1201 * t1125;
t1186 = t1194 * qJDD(1);
t1178 = t1201 * t1186;
t1144 = -qJD(3) * t1252 + t1178;
t1191 = t1194 ^ 2;
t1253 = t1182 * t1291;
t1264 = t1201 * t1203;
t1006 = -t1181 * pkin(3) - t1144 * qJ(4) + t1112 + (-pkin(3) * t1191 * t1264 + qJ(4) * t1253 - t1088) * t1198;
t1041 = t1201 * t1088 + t1198 * t1125;
t1141 = -t1182 * pkin(3) - qJ(4) * t1251;
t1281 = t1191 * t1203;
t1307 = t1198 ^ 2;
t1180 = t1307 * t1281;
t1259 = qJDD(1) * t1198;
t1224 = qJD(3) * t1290 + t1259;
t1217 = t1224 * t1194;
t1009 = -pkin(3) * t1180 - qJ(4) * t1217 + t1182 * t1141 + t1041;
t1230 = -0.2e1 * qJD(4) * t1138 + t1195 * t1006 - t1193 * t1009;
t1318 = pkin(3) * t1014 + t1230;
t1317 = t1193 * t1311;
t1197 = sin(qJ(5));
t1200 = cos(qJ(5));
t1084 = t1200 * t1137 + t1197 * t1138;
t1086 = -t1197 * t1137 + t1200 * t1138;
t1022 = t1086 * t1084;
t1170 = -qJDD(5) + t1181;
t1313 = -t1022 - t1170;
t1315 = t1197 * t1313;
t1314 = t1200 * t1313;
t1154 = t1182 * t1252;
t1116 = t1144 + t1154;
t1105 = t1201 * t1116;
t1097 = t1195 * t1144 - t1193 * t1217;
t1238 = t1193 * t1144 + t1195 * t1217;
t1002 = -t1084 * qJD(5) + t1200 * t1097 - t1197 * t1238;
t1175 = -qJD(5) + t1182;
t1071 = t1084 * t1175;
t1312 = t1071 + t1002;
t1124 = t1137 * t1182;
t1057 = -t1124 + t1097;
t1310 = t1124 + t1097;
t1239 = t1197 * t1097 + t1200 * t1238;
t970 = (qJD(5) + t1175) * t1086 + t1239;
t1081 = t1084 ^ 2;
t1082 = t1086 ^ 2;
t1134 = t1138 ^ 2;
t1169 = t1175 ^ 2;
t1289 = qJD(4) * t1137;
t1129 = -0.2e1 * t1289;
t1260 = t1193 * t1006 + t1195 * t1009;
t937 = t1129 + t1260;
t886 = t1193 * t937 + t1195 * t1230;
t1305 = pkin(3) * t886;
t1284 = t1182 * t1138;
t1053 = t1238 + t1284;
t998 = -t1053 * t1193 - t1195 * t1057;
t1304 = pkin(3) * t998;
t1303 = pkin(2) * t1194;
t1299 = t1196 * g(3);
t1298 = qJDD(1) * pkin(1);
t901 = pkin(4) * t1311 - pkin(7) * t1057 + t1230;
t1110 = -t1182 * pkin(4) - t1138 * pkin(7);
t913 = -t1309 * pkin(4) - pkin(7) * t1238 + t1182 * t1110 + t937;
t865 = t1197 * t913 - t1200 * t901;
t866 = t1197 * t901 + t1200 * t913;
t830 = t1197 * t866 - t1200 * t865;
t1297 = t1193 * t830;
t1296 = t1195 * t830;
t1021 = t1299 + qJDD(4) + pkin(3) * t1217 - qJ(4) * t1180 + (t1152 + (t1141 * t1201 + t1262) * qJD(1)) * t1194;
t956 = pkin(4) * t1238 - t1309 * pkin(7) + t1138 * t1110 + t1021;
t1295 = t1197 * t956;
t1294 = t1198 * t886;
t1293 = t1200 * t956;
t1292 = t1201 * t886;
t1287 = t1175 * t1086;
t1286 = t1175 * t1197;
t1285 = t1175 * t1200;
t1283 = t1182 * t1193;
t1282 = t1182 * t1195;
t1280 = t1193 * t1021;
t1078 = -t1094 + t1181;
t1279 = t1193 * t1078;
t1278 = t1194 * t1181;
t1277 = t1194 * t1196;
t1275 = t1195 * t1021;
t1274 = t1195 * t1078;
t1017 = -t1022 + t1170;
t1273 = t1197 * t1017;
t1087 = t1194 * t1227 + t1299;
t1272 = t1198 * t1087;
t1245 = t1198 * t1264;
t1168 = t1191 * t1245;
t1142 = -t1168 + t1181;
t1271 = t1198 * t1142;
t1143 = -t1168 - t1181;
t1270 = t1198 * t1143;
t1145 = -t1220 + t1298;
t1269 = t1199 * t1145;
t1268 = t1200 * t1017;
t1267 = t1201 * t1087;
t1266 = t1201 * t1142;
t1265 = t1201 * t1143;
t1263 = t1202 * t1145;
t1257 = t1199 * qJDD(1);
t1256 = t1202 * qJDD(1);
t831 = t1197 * t865 + t1200 * t866;
t807 = t1193 * t831 + t1296;
t829 = pkin(4) * t830;
t1255 = pkin(3) * t807 + t829;
t973 = -t1071 + t1002;
t910 = -t1197 * t970 - t1200 * t973;
t912 = t1197 * t973 - t1200 * t970;
t869 = t1193 * t912 + t1195 * t910;
t907 = pkin(4) * t910;
t1254 = pkin(3) * t869 + t907;
t1192 = t1201 ^ 2;
t1250 = t1192 * t1281;
t1249 = t1194 * t1022;
t1248 = t1194 * t1094;
t1247 = t1196 * t1022;
t1246 = t1196 * t1094;
t1244 = t1194 * t1258;
t1243 = t1145 + t1298;
t887 = -t1193 * t1230 + t1195 * t937;
t1242 = qJD(1) * (qJD(3) - t1182);
t1241 = qJD(1) * t1306 + t1152;
t1040 = t1198 * t1088 - t1112;
t989 = t1198 * t1040 + t1201 * t1041;
t1120 = t1241 * t1194 + t1299;
t1121 = t1241 * t1196 - t1300;
t1065 = t1194 * t1120 + t1196 * t1121;
t1237 = -t1199 * t1171 - t1202 * t1172;
t1190 = t1194 * t1191;
t1236 = t1190 * t1245;
t1235 = -pkin(2) * t1087 + pkin(6) * t989;
t1106 = -t1134 - t1308;
t1024 = t1195 * t1106 + t1279;
t1234 = pkin(3) * t1024 - t1260;
t1016 = -t1169 - t1081;
t954 = t1197 * t1016 + t1314;
t1233 = pkin(4) * t954 - t865;
t1167 = -t1199 * t1203 + t1256;
t1231 = -pkin(5) * t1167 - t1199 * g(3);
t1228 = t1196 * t1168;
t988 = -t1201 * t1040 + t1198 * t1041;
t1064 = t1196 * t1120 - t1194 * t1121;
t1226 = t1202 * t1171 - t1199 * t1172;
t1166 = t1202 * t1203 + t1257;
t1062 = -t1082 - t1169;
t978 = t1200 * t1062 + t1273;
t1225 = pkin(4) * t978 - t866;
t1205 = t1196 ^ 2;
t1157 = (t1191 + t1205) * t1196 * t1203;
t1223 = -t1199 * t1157 + t1196 * t1256;
t1222 = t1202 * t1157 + t1196 * t1257;
t955 = t1200 * t1016 - t1315;
t895 = t1193 * t955 + t1195 * t954;
t1221 = pkin(3) * t895 + t1233;
t1133 = -t1308 - t1250;
t1090 = -t1198 * t1133 + t1266;
t1117 = t1242 * t1276 - t1178;
t1219 = pkin(2) * t1117 + pkin(6) * t1090 + t1272;
t1148 = -t1180 - t1308;
t1103 = t1201 * t1148 - t1270;
t1155 = t1182 * t1251;
t1115 = t1155 - t1217;
t1218 = pkin(2) * t1115 + pkin(6) * t1103 - t1267;
t981 = -t1197 * t1062 + t1268;
t916 = t1193 * t981 + t1195 * t978;
t1216 = pkin(3) * t916 + t1225;
t1113 = -t1144 + t1154;
t1114 = t1155 + t1217;
t1060 = -t1198 * t1113 - t1201 * t1114;
t1150 = t1180 + t1250;
t1215 = pkin(2) * t1150 + pkin(6) * t1060 + t989;
t808 = t1195 * t831 - t1297;
t819 = -pkin(4) * t956 + pkin(7) * t831;
t791 = -pkin(3) * t956 - pkin(7) * t1297 + qJ(4) * t808 + t1195 * t819;
t794 = -pkin(7) * t1296 - qJ(4) * t807 - t1193 * t819;
t798 = -t1198 * t807 + t1201 * t808;
t1214 = -pkin(2) * t956 + pkin(6) * t798 + t1198 * t794 + t1201 * t791;
t994 = -t1081 - t1082;
t815 = -pkin(4) * t994 + pkin(7) * t912 + t831;
t817 = -pkin(7) * t910 - t830;
t871 = -t1193 * t910 + t1195 * t912;
t801 = -pkin(3) * t994 + qJ(4) * t871 + t1193 * t817 + t1195 * t815;
t802 = -qJ(4) * t869 - t1193 * t815 + t1195 * t817;
t835 = -t1198 * t869 + t1201 * t871;
t1213 = -pkin(2) * t994 + pkin(6) * t835 + t1198 * t802 + t1201 * t801;
t969 = (qJD(5) - t1175) * t1086 + t1239;
t883 = -pkin(4) * t969 + pkin(7) * t955 - t1293;
t896 = -t1193 * t954 + t1195 * t955;
t898 = -pkin(7) * t954 + t1295;
t824 = -pkin(3) * t969 + qJ(4) * t896 + t1193 * t898 + t1195 * t883;
t837 = -qJ(4) * t895 - t1193 * t883 + t1195 * t898;
t859 = -t1198 * t895 + t1201 * t896;
t1212 = -pkin(2) * t969 + pkin(6) * t859 + t1198 * t837 + t1201 * t824;
t885 = -pkin(4) * t1312 + pkin(7) * t981 + t1295;
t908 = -pkin(7) * t978 + t1293;
t917 = -t1193 * t978 + t1195 * t981;
t836 = -pkin(3) * t1312 + qJ(4) * t917 + t1193 * t908 + t1195 * t885;
t839 = -qJ(4) * t916 - t1193 * t885 + t1195 * t908;
t874 = -t1198 * t916 + t1201 * t917;
t1211 = -pkin(2) * t1312 + pkin(6) * t874 + t1198 * t839 + t1201 * t836;
t1052 = t1238 - t1284;
t1015 = t1195 * t1077 - t1317;
t938 = -pkin(3) * t1052 + qJ(4) * t1015 - t1275;
t953 = -t1198 * t1014 + t1201 * t1015;
t958 = -qJ(4) * t1014 + t1280;
t1210 = -pkin(2) * t1052 + pkin(6) * t953 + t1198 * t958 + t1201 * t938;
t1025 = -t1193 * t1106 + t1274;
t941 = -pkin(3) * t1310 + qJ(4) * t1025 + t1280;
t968 = -qJ(4) * t1024 + t1275;
t976 = -t1198 * t1024 + t1201 * t1025;
t1209 = -pkin(2) * t1310 + pkin(6) * t976 + t1198 * t968 + t1201 * t941;
t1066 = -t1134 - t1309;
t1000 = -t1053 * t1195 + t1193 * t1057;
t867 = -pkin(3) * t1066 + qJ(4) * t1000 + t887;
t872 = -qJ(4) * t998 - t886;
t932 = t1201 * t1000 - t1198 * t998;
t1208 = -pkin(2) * t1066 + pkin(6) * t932 + t1198 * t872 + t1201 * t867;
t843 = t1201 * t887 - t1294;
t881 = -pkin(3) * t1021 + qJ(4) * t887;
t1207 = -pkin(2) * t1021 + pkin(6) * t843 - qJ(4) * t1294 + t1201 * t881;
t1188 = t1205 * t1203;
t1187 = t1205 * qJDD(1);
t1185 = t1191 * qJDD(1);
t1177 = t1203 * t1277;
t1173 = 0.2e1 * t1244;
t1164 = -t1188 + t1281;
t1163 = t1188 + t1281;
t1162 = t1196 * t1181;
t1161 = t1187 - t1185;
t1160 = t1187 + t1185;
t1156 = (t1194 * t1205 + t1190) * t1203;
t1151 = -t1180 + t1250;
t1149 = t1308 - t1250;
t1147 = t1180 - t1308;
t1146 = -pkin(5) * t1166 + t1202 * g(3);
t1136 = t1167 * t1277;
t1135 = t1166 * t1277;
t1127 = t1202 * t1156 + t1194 * t1257;
t1126 = t1199 * t1156 - t1194 * t1256;
t1122 = (-t1192 - t1307) * t1253;
t1119 = -t1134 + t1308;
t1118 = -t1308 + t1309;
t1108 = t1198 * t1144 - t1192 * t1253;
t1107 = (-t1307 * t1182 * qJD(1) - t1201 * t1224) * t1194;
t1104 = (t1201 * t1242 + t1259) * t1276;
t1102 = t1201 * t1147 + t1271;
t1101 = -t1198 * t1149 + t1265;
t1100 = t1198 * t1148 + t1265;
t1099 = t1198 * t1147 - t1266;
t1098 = t1201 * t1149 + t1270;
t1096 = -qJ(2) * t1157 + t1243 * t1196;
t1095 = qJ(2) * t1156 - t1243 * t1194;
t1092 = t1134 - t1309;
t1089 = t1201 * t1133 + t1271;
t1076 = t1196 * t1105 + t1236;
t1075 = t1196 * t1104 - t1236;
t1074 = t1194 * t1105 - t1228;
t1073 = t1194 * t1104 + t1228;
t1070 = -t1082 + t1169;
t1069 = t1081 - t1169;
t1068 = (t1137 * t1195 - t1138 * t1193) * t1182;
t1067 = (t1137 * t1193 + t1138 * t1195) * t1182;
t1061 = t1201 * t1115 - t1198 * t1116;
t1059 = t1198 * t1115 + t1105;
t1058 = t1201 * t1113 - t1198 * t1114;
t1051 = t1195 * t1097 + t1138 * t1283;
t1050 = t1193 * t1097 - t1138 * t1282;
t1049 = -t1137 * t1282 + t1193 * t1238;
t1048 = -t1137 * t1283 - t1195 * t1238;
t1047 = t1196 * t1103 - t1194 * t1115;
t1046 = t1196 * t1102 - t1194 * t1114;
t1045 = t1196 * t1101 - t1194 * t1113;
t1044 = t1194 * t1103 + t1196 * t1115;
t1043 = t1194 * t1102 + t1196 * t1114;
t1042 = t1194 * t1101 + t1196 * t1113;
t1038 = t1196 * t1090 - t1194 * t1117;
t1037 = t1194 * t1090 + t1196 * t1117;
t1036 = pkin(1) * t1145 + qJ(2) * t1065;
t1035 = t1195 * t1118 + t1279;
t1034 = -t1193 * t1119 + t1316;
t1033 = t1193 * t1118 - t1274;
t1032 = t1195 * t1119 + t1317;
t1031 = pkin(1) * t1163 + qJ(2) * t1160 + t1065;
t1030 = -pkin(6) * t1100 + t1272;
t1029 = t1196 * t1061 + t1194 * t1151;
t1028 = t1196 * t1060 - t1194 * t1150;
t1027 = t1194 * t1061 - t1196 * t1151;
t1026 = t1194 * t1060 + t1196 * t1150;
t1023 = -pkin(6) * t1089 + t1267;
t1020 = t1082 - t1081;
t1013 = (t1084 * t1200 - t1086 * t1197) * t1175;
t1012 = (t1084 * t1197 + t1086 * t1200) * t1175;
t1011 = -pkin(2) * t1100 + t1040;
t1010 = -pkin(2) * t1089 + t1041;
t1008 = -t1198 * t1067 + t1201 * t1068;
t1007 = t1201 * t1067 + t1198 * t1068;
t1001 = -t1086 * qJD(5) - t1239;
t999 = -t1195 * t1052 - t1193 * t1310;
t997 = -t1193 * t1052 + t1195 * t1310;
t996 = t1196 * t1008 - t1278;
t995 = t1194 * t1008 + t1162;
t993 = -t1198 * t1050 + t1201 * t1051;
t992 = -t1198 * t1048 + t1201 * t1049;
t991 = t1201 * t1050 + t1198 * t1051;
t990 = t1201 * t1048 + t1198 * t1049;
t987 = t1200 * t1069 + t1273;
t986 = -t1197 * t1070 + t1314;
t985 = t1197 * t1069 - t1268;
t984 = t1200 * t1070 + t1315;
t983 = -t1198 * t1033 + t1201 * t1035;
t982 = -t1198 * t1032 + t1201 * t1034;
t980 = t1201 * t1033 + t1198 * t1035;
t979 = t1201 * t1032 + t1198 * t1034;
t975 = t1201 * t1024 + t1198 * t1025;
t967 = -pkin(1) * t1044 - t1218;
t966 = t1200 * t1002 + t1086 * t1286;
t965 = t1197 * t1002 - t1086 * t1285;
t964 = -t1197 * t1001 - t1084 * t1285;
t963 = t1200 * t1001 - t1084 * t1286;
t962 = t1196 * t993 + t1248;
t961 = t1196 * t992 - t1248;
t960 = t1194 * t993 - t1246;
t959 = t1194 * t992 + t1246;
t957 = -pkin(1) * t1037 - t1219;
t952 = t1194 * t1087 + t1196 * t989;
t951 = t1201 * t1014 + t1198 * t1015;
t950 = -t1196 * t1087 + t1194 * t989;
t948 = -pkin(6) * t1058 - t988;
t947 = -t1193 * t1012 + t1195 * t1013;
t946 = t1195 * t1012 + t1193 * t1013;
t945 = -t1194 * t1053 + t1196 * t983;
t944 = t1194 * t1057 + t1196 * t982;
t943 = t1196 * t1053 + t1194 * t983;
t942 = -t1196 * t1057 + t1194 * t982;
t940 = t1194 * t1310 + t1196 * t976;
t939 = t1194 * t976 - t1196 * t1310;
t935 = t1194 * t1052 + t1196 * t953;
t934 = -t1196 * t1052 + t1194 * t953;
t933 = -qJ(2) * t1044 - t1194 * t1011 + t1196 * t1030;
t931 = -t1198 * t997 + t1201 * t999;
t930 = t1198 * t1000 + t1201 * t998;
t929 = t1198 * t999 + t1201 * t997;
t927 = -qJ(2) * t1037 - t1194 * t1010 + t1196 * t1023;
t926 = -pkin(1) * t1100 + qJ(2) * t1047 + t1196 * t1011 + t1194 * t1030;
t925 = -pkin(1) * t1026 - t1215;
t924 = t1194 * t1092 + t1196 * t931;
t923 = -t1196 * t1092 + t1194 * t931;
t922 = -t1193 * t985 + t1195 * t987;
t921 = -t1193 * t984 + t1195 * t986;
t920 = t1193 * t987 + t1195 * t985;
t919 = t1193 * t986 + t1195 * t984;
t918 = -pkin(1) * t1089 + qJ(2) * t1038 + t1196 * t1010 + t1194 * t1023;
t915 = t1194 * t1066 + t1196 * t932;
t914 = -t1196 * t1066 + t1194 * t932;
t911 = -t1197 * t1312 - t1200 * t969;
t909 = -t1197 * t969 + t1200 * t1312;
t906 = -qJ(2) * t1026 + t1058 * t1303 + t1196 * t948;
t905 = -t1193 * t965 + t1195 * t966;
t904 = -t1193 * t963 + t1195 * t964;
t903 = t1193 * t966 + t1195 * t965;
t902 = t1193 * t964 + t1195 * t963;
t899 = -pkin(1) * t950 - t1235;
t897 = -pkin(2) * t930 - t1304;
t894 = qJ(2) * t1028 + t1194 * t948 + (-pkin(1) - t1302) * t1058;
t893 = -t1198 * t946 + t1201 * t947;
t892 = t1198 * t947 + t1201 * t946;
t891 = -t1194 * t1170 + t1196 * t893;
t890 = t1196 * t1170 + t1194 * t893;
t889 = -qJ(2) * t950 + (-pkin(6) * t1196 + t1303) * t988;
t888 = -pkin(2) * t975 + t1129 - t1234;
t884 = -pkin(2) * t951 - t1318;
t882 = -pkin(6) * t975 - t1198 * t941 + t1201 * t968;
t880 = -t1198 * t920 + t1201 * t922;
t879 = -t1198 * t919 + t1201 * t921;
t878 = t1198 * t922 + t1201 * t920;
t877 = t1198 * t921 + t1201 * t919;
t876 = -pkin(6) * t951 - t1198 * t938 + t1201 * t958;
t875 = qJ(2) * t952 + t1229 * t988;
t873 = t1198 * t917 + t1201 * t916;
t870 = -t1193 * t909 + t1195 * t911;
t868 = t1193 * t911 + t1195 * t909;
t863 = -t1198 * t903 + t1201 * t905;
t862 = -t1198 * t902 + t1201 * t904;
t861 = t1198 * t905 + t1201 * t903;
t860 = t1198 * t904 + t1201 * t902;
t858 = t1198 * t896 + t1201 * t895;
t857 = -t1194 * t970 + t1196 * t880;
t856 = t1194 * t973 + t1196 * t879;
t855 = t1194 * t880 + t1196 * t970;
t854 = t1194 * t879 - t1196 * t973;
t853 = t1196 * t863 + t1249;
t852 = t1196 * t862 - t1249;
t851 = t1194 * t863 - t1247;
t850 = t1194 * t862 + t1247;
t849 = t1194 * t1312 + t1196 * t874;
t848 = t1194 * t874 - t1196 * t1312;
t847 = -pkin(1) * t939 - t1209;
t846 = -pkin(1) * t934 - t1210;
t845 = t1194 * t969 + t1196 * t859;
t844 = t1194 * t859 - t1196 * t969;
t842 = t1198 * t887 + t1292;
t841 = t1194 * t1021 + t1196 * t843;
t840 = -t1196 * t1021 + t1194 * t843;
t838 = -qJ(2) * t939 - t1194 * t888 + t1196 * t882;
t834 = -t1198 * t868 + t1201 * t870;
t833 = t1198 * t871 + t1201 * t869;
t832 = t1198 * t870 + t1201 * t868;
t828 = -qJ(2) * t934 - t1194 * t884 + t1196 * t876;
t827 = -pkin(1) * t975 + qJ(2) * t940 + t1194 * t882 + t1196 * t888;
t826 = t1194 * t1020 + t1196 * t834;
t825 = -t1196 * t1020 + t1194 * t834;
t823 = t1194 * t994 + t1196 * t835;
t822 = t1194 * t835 - t1196 * t994;
t821 = -pkin(2) * t842 - t1305;
t820 = -pkin(6) * t930 - t1198 * t867 + t1201 * t872;
t818 = -pkin(1) * t951 + qJ(2) * t935 + t1194 * t876 + t1196 * t884;
t816 = -pkin(2) * t873 - t1216;
t814 = -pkin(1) * t914 - t1208;
t813 = -pkin(2) * t858 - t1221;
t812 = -pkin(6) * t842 - qJ(4) * t1292 - t1198 * t881;
t811 = -qJ(2) * t914 - t1194 * t897 + t1196 * t820;
t810 = -pkin(2) * t833 - t1254;
t809 = -pkin(1) * t930 + qJ(2) * t915 + t1194 * t820 + t1196 * t897;
t806 = -pkin(6) * t873 - t1198 * t836 + t1201 * t839;
t805 = -pkin(1) * t840 - t1207;
t804 = -pkin(6) * t858 - t1198 * t824 + t1201 * t837;
t803 = -pkin(1) * t848 - t1211;
t800 = -pkin(1) * t844 - t1212;
t799 = -qJ(2) * t840 - t1194 * t821 + t1196 * t812;
t797 = t1198 * t808 + t1201 * t807;
t796 = -qJ(2) * t848 - t1194 * t816 + t1196 * t806;
t795 = -pkin(1) * t842 + qJ(2) * t841 + t1194 * t812 + t1196 * t821;
t793 = t1194 * t956 + t1196 * t798;
t792 = t1194 * t798 - t1196 * t956;
t790 = -pkin(1) * t873 + qJ(2) * t849 + t1194 * t806 + t1196 * t816;
t789 = -qJ(2) * t844 - t1194 * t813 + t1196 * t804;
t788 = -pkin(1) * t858 + qJ(2) * t845 + t1194 * t804 + t1196 * t813;
t787 = -pkin(2) * t797 - t1255;
t786 = -pkin(6) * t833 - t1198 * t801 + t1201 * t802;
t785 = -pkin(1) * t822 - t1213;
t784 = -qJ(2) * t822 - t1194 * t810 + t1196 * t786;
t783 = -pkin(1) * t833 + qJ(2) * t823 + t1194 * t786 + t1196 * t810;
t782 = -pkin(6) * t797 - t1198 * t791 + t1201 * t794;
t781 = -pkin(1) * t792 - t1214;
t780 = -qJ(2) * t792 - t1194 * t787 + t1196 * t782;
t779 = -pkin(1) * t797 + qJ(2) * t793 + t1194 * t782 + t1196 * t787;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1167, 0, -t1166, 0, t1231, -t1146, -t1226, -pkin(5) * t1226, t1136, t1202 * t1161 + t1199 * t1164, t1127, -t1136, t1222, 0, -pkin(5) * t1223 - t1199 * t1120 - t1194 * t1263, -pkin(5) * t1126 - t1199 * t1121 - t1196 * t1263, t1202 * t1064 - pkin(5) * (t1199 * t1160 + t1202 * t1163), -pkin(5) * (t1199 * t1065 + t1263) - (t1199 * pkin(1) - t1202 * qJ(2)) * t1064, t1202 * t1076 + t1199 * t1108, t1202 * t1029 + t1199 * t1059, t1202 * t1045 + t1199 * t1098, t1202 * t1075 + t1199 * t1107, t1202 * t1046 + t1199 * t1099, -t1199 * t1122 - t1202 * t1278, t1202 * t933 - t1199 * t967 - pkin(5) * (t1199 * t1047 - t1202 * t1100), t1202 * t927 - t1199 * t957 - pkin(5) * (t1199 * t1038 - t1202 * t1089), t1202 * t906 - t1199 * t925 - pkin(5) * (t1199 * t1028 - t1202 * t1058), t1202 * t889 - t1199 * t899 - pkin(5) * (t1199 * t952 - t1202 * t988), t1199 * t991 + t1202 * t962, t1199 * t929 + t1202 * t924, t1199 * t979 + t1202 * t944, t1199 * t990 + t1202 * t961, t1199 * t980 + t1202 * t945, t1199 * t1007 + t1202 * t996, t1202 * t828 - t1199 * t846 - pkin(5) * (t1199 * t935 - t1202 * t951), t1202 * t838 - t1199 * t847 - pkin(5) * (t1199 * t940 - t1202 * t975), t1202 * t811 - t1199 * t814 - pkin(5) * (t1199 * t915 - t1202 * t930), t1202 * t799 - t1199 * t805 - pkin(5) * (t1199 * t841 - t1202 * t842), t1199 * t861 + t1202 * t853, t1199 * t832 + t1202 * t826, t1199 * t877 + t1202 * t856, t1199 * t860 + t1202 * t852, t1199 * t878 + t1202 * t857, t1199 * t892 + t1202 * t891, t1202 * t789 - t1199 * t800 - pkin(5) * (t1199 * t845 - t1202 * t858), t1202 * t796 - t1199 * t803 - pkin(5) * (t1199 * t849 - t1202 * t873), t1202 * t784 - t1199 * t785 - pkin(5) * (t1199 * t823 - t1202 * t833), t1202 * t780 - t1199 * t781 - pkin(5) * (t1199 * t793 - t1202 * t797); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1166, 0, t1167, 0, t1146, t1231, t1237, pkin(5) * t1237, t1135, t1199 * t1161 - t1202 * t1164, t1126, -t1135, -t1223, 0, -pkin(5) * t1222 + t1202 * t1120 - t1194 * t1269, pkin(5) * t1127 + t1202 * t1121 - t1196 * t1269, t1199 * t1064 + pkin(5) * (t1202 * t1160 - t1199 * t1163), pkin(5) * (t1202 * t1065 - t1269) - (-t1202 * pkin(1) - t1199 * qJ(2)) * t1064, t1199 * t1076 - t1202 * t1108, t1199 * t1029 - t1202 * t1059, t1199 * t1045 - t1202 * t1098, t1199 * t1075 - t1202 * t1107, t1199 * t1046 - t1202 * t1099, t1202 * t1122 - t1199 * t1278, t1199 * t933 + t1202 * t967 + pkin(5) * (t1202 * t1047 + t1199 * t1100), t1199 * t927 + t1202 * t957 + pkin(5) * (t1202 * t1038 + t1199 * t1089), t1199 * t906 + t1202 * t925 + pkin(5) * (t1202 * t1028 + t1199 * t1058), t1199 * t889 + t1202 * t899 + pkin(5) * (t1199 * t988 + t1202 * t952), t1199 * t962 - t1202 * t991, t1199 * t924 - t1202 * t929, t1199 * t944 - t1202 * t979, t1199 * t961 - t1202 * t990, t1199 * t945 - t1202 * t980, -t1202 * t1007 + t1199 * t996, t1199 * t828 + t1202 * t846 + pkin(5) * (t1199 * t951 + t1202 * t935), t1199 * t838 + t1202 * t847 + pkin(5) * (t1199 * t975 + t1202 * t940), t1199 * t811 + t1202 * t814 + pkin(5) * (t1199 * t930 + t1202 * t915), t1199 * t799 + t1202 * t805 + pkin(5) * (t1199 * t842 + t1202 * t841), t1199 * t853 - t1202 * t861, t1199 * t826 - t1202 * t832, t1199 * t856 - t1202 * t877, t1199 * t852 - t1202 * t860, t1199 * t857 - t1202 * t878, t1199 * t891 - t1202 * t892, t1199 * t789 + t1202 * t800 + pkin(5) * (t1199 * t858 + t1202 * t845), t1199 * t796 + t1202 * t803 + pkin(5) * (t1199 * t873 + t1202 * t849), t1199 * t784 + t1202 * t785 + pkin(5) * (t1199 * t833 + t1202 * t823), t1199 * t780 + t1202 * t781 + pkin(5) * (t1199 * t797 + t1202 * t793); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1171, t1172, 0, 0, t1185, t1173, 0, t1187, 0, 0, t1096, t1095, t1031, t1036, t1074, t1027, t1042, t1073, t1043, t1162, t926, t918, t894, t875, t960, t923, t942, t959, t943, t995, t818, t827, t809, t795, t851, t825, t854, t850, t855, t890, t788, t790, t783, t779; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1203, 0, 0, -g(3), -t1171, 0, t1244, t1161, t1156, -t1244, t1157, 0, -t1194 * t1145, -t1196 * t1145, t1064, qJ(2) * t1064, t1076, t1029, t1045, t1075, t1046, -t1278, t933, t927, t906, t889, t962, t924, t944, t961, t945, t996, t828, t838, t811, t799, t853, t826, t856, t852, t857, t891, t789, t796, t784, t780; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1203, 0, qJDD(1), 0, g(3), 0, -t1172, 0, t1177, -t1164, -t1186, -t1177, -t1258, 0, t1120, t1121, 0, pkin(1) * t1064, -t1108, -t1059, -t1098, -t1107, -t1099, t1122, t967, t957, t925, t899, -t991, -t929, -t979, -t990, -t980, -t1007, t846, t847, t814, t805, -t861, -t832, -t877, -t860, -t878, -t892, t800, t803, t785, t781; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1171, t1172, 0, 0, t1185, t1173, 0, t1187, 0, 0, t1096, t1095, t1031, t1036, t1074, t1027, t1042, t1073, t1043, t1162, t926, t918, t894, t875, t960, t923, t942, t959, t943, t995, t818, t827, t809, t795, t851, t825, t854, t850, t855, t890, t788, t790, t783, t779; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1186, t1258, t1177, 0, t1188, 0, 0, -t1145, t1120, 0, t1105, t1061, t1101, t1104, t1102, 0, t1030, t1023, t948, -pkin(6) * t988, t993, t931, t982, t992, t983, t1008, t876, t882, t820, t812, t863, t834, t879, t862, t880, t893, t804, t806, t786, t782; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1186, -t1281, t1258, -t1177, 0, t1145, 0, t1121, 0, -t1168, -t1151, t1113, t1168, t1114, t1181, t1011, t1010, -pkin(2) * t1058, -pkin(2) * t988, -t1094, -t1092, -t1057, t1094, t1053, t1181, t884, t888, t897, t821, -t1022, -t1020, -t973, t1022, t970, t1170, t813, t816, t810, t787; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1177, t1164, t1186, t1177, t1258, 0, -t1120, -t1121, 0, 0, t1108, t1059, t1098, t1107, t1099, -t1122, t1218, t1219, t1215, t1235, t991, t929, t979, t990, t980, t1007, t1210, t1209, t1208, t1207, t861, t832, t877, t860, t878, t892, t1212, t1211, t1213, t1214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1144, t1115, t1143, -t1154, t1147, t1154, 0, t1087, t1040, 0, t1051, t999, t1034, t1049, t1035, t1068, t958, t968, t872, -qJ(4) * t886, t905, t870, t921, t904, t922, t947, t837, t839, t802, t794; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1155, t1116, t1149, -t1217, -t1142, t1155, -t1087, 0, t1041, 0, t1050, t997, t1032, t1048, t1033, t1067, t938, t941, t867, t881, t903, t868, t919, t902, t920, t946, t824, t836, t801, t791; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1168, t1151, -t1113, -t1168, -t1114, -t1181, -t1040, -t1041, 0, 0, t1094, t1092, t1057, -t1094, -t1053, -t1181, t1318, t1234 + 0.2e1 * t1289, t1304, t1305, t1022, t1020, t973, -t1022, -t970, -t1170, t1221, t1216, t1254, t1255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1097, -t1052, t1311, -t1124, t1118, t1124, 0, t1021, -t1230, 0, t966, t911, t986, t964, t987, t1013, t898, t908, t817, -pkin(7) * t830; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1284, t1310, t1119, -t1238, -t1078, t1284, -t1021, 0, t937, 0, t965, t909, t984, t963, t985, t1012, t883, t885, t815, t819; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1094, t1092, t1057, -t1094, -t1053, -t1181, t1230, -t937, 0, 0, t1022, t1020, t973, -t1022, -t970, -t1170, t1233, t1225, t907, t829; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1002, -t969, t1313, -t1071, t1069, t1071, 0, t956, t865, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1287, t1312, t1070, t1001, -t1017, t1287, -t956, 0, t866, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1022, t1020, t973, -t1022, -t970, -t1170, -t865, -t866, 0, 0;];
m_new_reg = t1;
