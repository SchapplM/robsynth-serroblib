% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPRPP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPRPP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:32
% EndTime: 2019-12-31 18:09:48
% DurationCPUTime: 17.12s
% Computational Cost: add. (34723->563), mult. (76084->683), div. (0->0), fcn. (47816->8), ass. (0->390)
t1209 = sin(qJ(1));
t1211 = cos(qJ(1));
t1203 = sin(pkin(7));
t1205 = cos(pkin(7));
t1210 = cos(qJ(3));
t1288 = qJD(1) * t1210;
t1188 = qJD(3) * t1288;
t1208 = sin(qJ(3));
t1257 = t1208 * qJDD(1);
t1164 = t1188 + t1257;
t1202 = sin(pkin(8));
t1204 = cos(pkin(8));
t1289 = qJD(1) * t1208;
t1251 = qJD(3) * t1289;
t1256 = t1210 * qJDD(1);
t1220 = -t1251 + t1256;
t1101 = t1204 * t1164 + t1202 * t1220;
t1148 = t1202 * t1289 - t1204 * t1288;
t1287 = qJD(3) * t1148;
t1316 = -t1287 + t1101;
t1150 = (t1202 * t1210 + t1204 * t1208) * qJD(1);
t1147 = t1150 ^ 2;
t1212 = qJD(3) ^ 2;
t1084 = t1212 + t1147;
t1282 = t1150 * t1148;
t1312 = qJDD(3) + t1282;
t1268 = t1204 * t1312;
t1021 = t1202 * t1084 - t1268;
t1275 = t1202 * t1312;
t996 = t1204 * t1084 + t1275;
t958 = t1210 * t1021 + t1208 * t996;
t913 = t1203 * t958 - t1205 * t1316;
t919 = t1203 * t1316 + t1205 * t958;
t1411 = pkin(5) * (t1209 * t913 - t1211 * t919);
t1410 = pkin(5) * (t1209 * t919 + t1211 * t913);
t1409 = pkin(1) * t913;
t1408 = qJ(2) * t913;
t952 = t1208 * t1021 - t1210 * t996;
t1407 = pkin(1) * t952 - qJ(2) * t919;
t1303 = t1148 ^ 2;
t1318 = t1147 - t1303;
t1100 = t1202 * t1164 - t1204 * t1220;
t1143 = qJD(3) * t1150;
t1042 = t1100 + t1143;
t971 = -t1202 * t1042 + t1204 * t1316;
t1269 = t1204 * t1042;
t1277 = t1202 * t1316;
t973 = t1269 + t1277;
t895 = t1208 * t971 + t1210 * t973;
t873 = t1203 * t895 + t1205 * t1318;
t875 = -t1203 * t1318 + t1205 * t895;
t1406 = t1209 * t875 + t1211 * t873;
t1405 = t1209 * t873 - t1211 * t875;
t1043 = t1100 - t1143;
t1131 = t1303 - t1212;
t1010 = t1202 * t1131 + t1268;
t1016 = t1204 * t1131 - t1275;
t954 = t1208 * t1010 - t1210 * t1016;
t909 = -t1205 * t1043 + t1203 * t954;
t915 = t1203 * t1043 + t1205 * t954;
t1404 = t1209 * t915 + t1211 * t909;
t1403 = t1209 * t909 - t1211 * t915;
t1399 = pkin(2) * t952;
t1398 = pkin(6) * t952;
t1397 = pkin(2) * t1316 - pkin(6) * t958;
t1040 = -t1303 - t1147;
t1315 = t1287 + t1101;
t1347 = -t1204 * t1043 + t1202 * t1315;
t1348 = -t1202 * t1043 - t1204 * t1315;
t1363 = -t1208 * t1348 + t1210 * t1347;
t1378 = -t1205 * t1040 + t1203 * t1363;
t1396 = pkin(1) * t1378;
t1395 = qJ(2) * t1378;
t1362 = t1208 * t1347 + t1210 * t1348;
t1376 = t1203 * t1040 + t1205 * t1363;
t1394 = -pkin(1) * t1362 + qJ(2) * t1376;
t1132 = t1147 - t1212;
t1313 = qJDD(3) - t1282;
t1274 = t1202 * t1313;
t1351 = -t1204 * t1132 + t1274;
t1067 = t1204 * t1313;
t1352 = t1202 * t1132 + t1067;
t1360 = -t1208 * t1351 + t1210 * t1352;
t1377 = t1203 * t1315 + t1205 * t1360;
t1379 = t1203 * t1360 - t1205 * t1315;
t1393 = -t1209 * t1379 + t1211 * t1377;
t1392 = t1209 * t1377 + t1211 * t1379;
t1391 = pkin(5) * (-t1209 * t1378 + t1211 * t1376);
t1390 = pkin(5) * (t1209 * t1376 + t1211 * t1378);
t1389 = pkin(3) * t996;
t1314 = -t1303 - t1212;
t1327 = t1204 * t1314 - t1274;
t1332 = t1202 * t1314 + t1067;
t1346 = -t1208 * t1332 + t1210 * t1327;
t1365 = -t1205 * t1042 + t1203 * t1346;
t1387 = pkin(1) * t1365;
t1386 = pkin(2) * t1362;
t1385 = pkin(6) * t1362;
t1384 = qJ(4) * t996;
t1383 = qJ(2) * t1365;
t1382 = qJ(4) * t1021;
t1381 = -pkin(2) * t1040 + pkin(6) * t1363;
t1345 = t1208 * t1327 + t1210 * t1332;
t1364 = t1203 * t1042 + t1205 * t1346;
t1380 = -pkin(1) * t1345 + qJ(2) * t1364;
t891 = t1208 * t973 - t1210 * t971;
t947 = t1210 * t1010 + t1208 * t1016;
t1375 = pkin(5) * (-t1209 * t1365 + t1211 * t1364);
t1374 = pkin(5) * (t1209 * t1364 + t1211 * t1365);
t1300 = pkin(3) * t1348;
t1372 = pkin(6) * t1345;
t1371 = qJ(4) * t1348;
t1357 = pkin(3) * t1332;
t1368 = -pkin(2) * t1345 - t1357;
t1367 = -pkin(2) * t1042 + pkin(6) * t1346;
t1366 = -pkin(3) * t1040 + qJ(4) * t1347;
t1361 = t1208 * t1352 + t1210 * t1351;
t1356 = qJ(4) * t1327;
t1355 = qJ(4) * t1332;
t1178 = t1211 * g(1) + t1209 * g(2);
t1304 = qJD(1) ^ 2;
t1161 = -t1304 * pkin(1) - t1178;
t1177 = t1209 * g(1) - t1211 * g(2);
t1223 = qJDD(1) * pkin(1) + t1177;
t1098 = t1203 * t1161 - t1205 * t1223;
t1099 = t1205 * t1161 + t1203 * t1223;
t1245 = t1203 * t1098 + t1205 * t1099;
t1006 = t1205 * t1098 - t1203 * t1099;
t1260 = t1211 * t1006;
t1350 = -t1209 * t1245 + t1260;
t1264 = t1209 * t1006;
t1349 = t1211 * t1245 + t1264;
t1254 = t1203 * t1282;
t1253 = t1204 * t1287;
t1221 = t1202 * t1100 + t1253;
t1283 = t1148 * t1202;
t1228 = qJD(3) * t1283 - t1204 * t1100;
t1308 = -t1208 * t1228 + t1210 * t1221;
t1326 = t1205 * t1308 - t1254;
t1252 = t1205 * t1282;
t1329 = t1203 * t1308 + t1252;
t1344 = -t1209 * t1329 + t1211 * t1326;
t1343 = t1209 * t1326 + t1211 * t1329;
t1258 = t1205 * qJDD(3);
t1281 = t1150 * t1204;
t1218 = (-t1281 - t1283) * qJD(3);
t1126 = t1202 * t1143;
t1227 = t1126 - t1253;
t1309 = -t1208 * t1218 + t1210 * t1227;
t1328 = t1203 * t1309 - t1258;
t1190 = t1203 * qJDD(3);
t1330 = t1205 * t1309 + t1190;
t1342 = t1209 * t1330 + t1211 * t1328;
t1341 = -t1209 * t1328 + t1211 * t1330;
t1339 = qJ(5) * t1316;
t1167 = t1203 * qJDD(1) + t1205 * t1304;
t1199 = g(3) - qJDD(2);
t1130 = qJ(2) * t1167 - t1205 * t1199;
t1168 = t1205 * qJDD(1) - t1203 * t1304;
t1225 = -qJ(2) * t1168 - t1203 * t1199;
t1317 = t1211 * t1167 + t1209 * t1168;
t1331 = pkin(5) * t1317 + t1211 * t1130 - t1209 * t1225;
t1104 = -t1209 * t1167 + t1211 * t1168;
t1325 = -pkin(5) * t1104 + t1209 * t1130 + t1211 * t1225;
t1069 = -t1304 * pkin(2) + qJDD(1) * pkin(6) + t1099;
t1050 = t1208 * t1069 + t1210 * t1199;
t1052 = t1210 * t1069 - t1208 * t1199;
t980 = t1208 * t1050 + t1210 * t1052;
t1320 = -pkin(4) * t1313 - qJ(5) * t1314;
t1286 = qJD(4) * t1148;
t1138 = -0.2e1 * t1286;
t1284 = qJD(5) * qJD(3);
t1319 = t1138 + 0.2e1 * t1284;
t1310 = t1208 * t1227 + t1210 * t1218;
t1307 = t1208 * t1221 + t1210 * t1228;
t1026 = qJD(3) * t1281 + t1202 * t1101;
t1027 = t1204 * t1101 - t1126;
t968 = -t1208 * t1026 + t1210 * t1027;
t1229 = t1203 * t968 - t1252;
t1230 = t1205 * t968 + t1254;
t1306 = t1209 * t1230 + t1211 * t1229;
t1305 = -t1209 * t1229 + t1211 * t1230;
t1302 = t1210 ^ 2;
t1184 = t1210 * t1304 * t1208;
t1175 = qJDD(3) + t1184;
t992 = (-t1164 + t1188) * qJ(4) + t1175 * pkin(3) - t1050;
t1174 = qJD(3) * pkin(3) - qJ(4) * t1289;
t1192 = t1302 * t1304;
t994 = -pkin(3) * t1192 + qJ(4) * t1220 - qJD(3) * t1174 + t1052;
t1247 = t1202 * t994 - t1204 * t992;
t1285 = qJD(4) * t1150;
t920 = t1247 + 0.2e1 * t1285;
t1294 = t1202 * t992 + t1204 * t994;
t921 = t1138 + t1294;
t860 = t1202 * t921 - t1204 * t920;
t1301 = pkin(3) * t860;
t1297 = pkin(4) * t1204;
t1081 = t1148 * pkin(4) - t1150 * qJ(5);
t1226 = -t1212 * pkin(4) + qJDD(3) * qJ(5) - t1148 * t1081 + t1294;
t878 = t1226 + t1319;
t1219 = -qJDD(3) * pkin(4) - t1212 * qJ(5) + qJDD(5) + t1247;
t880 = (0.2e1 * qJD(4) + t1081) * t1150 + t1219;
t1295 = -pkin(4) * t880 + qJ(5) * t878;
t1293 = t1208 * t860;
t1292 = t1210 * t860;
t1068 = -qJDD(1) * pkin(2) - t1304 * pkin(6) + t1098;
t1291 = -pkin(2) * t1068 + pkin(6) * t980;
t1290 = qJD(1) * qJD(3);
t1198 = t1208 ^ 2;
t1280 = t1198 * t1304;
t1003 = -t1220 * pkin(3) - qJ(4) * t1192 + t1174 * t1289 + qJDD(4) + t1068;
t1279 = t1202 * t1003;
t1270 = t1204 * t1003;
t1063 = t1208 * t1068;
t1266 = t1208 * t1175;
t1176 = qJDD(3) - t1184;
t1265 = t1208 * t1176;
t1064 = t1210 * t1068;
t1165 = -0.2e1 * t1251 + t1256;
t1110 = t1210 * t1165;
t1261 = t1210 * t1176;
t1259 = -pkin(4) * t1315 - qJ(5) * t1043;
t1255 = t1198 + t1302;
t1181 = -t1212 - t1280;
t1119 = -t1208 * t1181 - t1261;
t1163 = 0.2e1 * t1188 + t1257;
t1250 = -pkin(2) * t1163 + pkin(6) * t1119 + t1063;
t1183 = -t1192 - t1212;
t1117 = t1210 * t1183 - t1266;
t1249 = pkin(2) * t1165 + pkin(6) * t1117 - t1064;
t1248 = -qJ(5) * t1202 - pkin(3);
t861 = t1202 * t920 + t1204 * t921;
t843 = t1202 * t880 + t1204 * t878;
t1214 = t1100 * pkin(4) + t1003 - t1339;
t905 = (pkin(4) * qJD(3) - 0.2e1 * qJD(5)) * t1150 + t1214;
t818 = qJ(4) * t843 + (t1248 - t1297) * t905;
t842 = t1202 * t878 - t1204 * t880;
t822 = -t1208 * t842 + t1210 * t843;
t826 = -qJ(4) * t842 + (pkin(4) * t1202 - qJ(5) * t1204) * t905;
t1246 = -pkin(2) * t905 + pkin(6) * t822 + t1208 * t826 + t1210 * t818;
t1244 = -t1209 * t1177 - t1211 * t1178;
t866 = -pkin(4) * t1040 + t878;
t867 = -qJ(5) * t1040 + t880;
t834 = t1202 * t867 + t1204 * t866 + t1366;
t837 = -t1202 * t866 + t1204 * t867 - t1371;
t1243 = t1208 * t837 + t1210 * t834 + t1381;
t840 = t1366 + t861;
t846 = -t860 - t1371;
t1242 = t1208 * t846 + t1210 * t840 + t1381;
t1213 = 0.2e1 * qJD(5) * t1150 - t1214;
t881 = -pkin(4) * t1143 + t1213 + t1339;
t852 = -t1382 + t1202 * t881 + (pkin(3) + t1297) * t1316;
t858 = -pkin(4) * t1277 + t1204 * t881 - t1384;
t1241 = t1208 * t858 + t1210 * t852 + t1397;
t882 = (-t1042 - t1143) * pkin(4) + t1213;
t856 = t1042 * t1248 + t1204 * t882 + t1356;
t864 = -qJ(5) * t1269 - t1202 * t882 - t1355;
t1240 = t1208 * t864 + t1210 * t856 + t1367;
t900 = -pkin(3) * t1042 - t1270 + t1356;
t938 = t1279 - t1355;
t1239 = t1208 * t938 + t1210 * t900 + t1367;
t906 = -pkin(3) * t1316 + t1279 + t1382;
t942 = t1270 + t1384;
t1238 = t1208 * t942 + t1210 * t906 - t1397;
t1237 = pkin(3) * t842 + t1295;
t1236 = -t1294 - t1389;
t1235 = t1203 * t1184;
t1234 = t1205 * t1184;
t1169 = t1255 * qJDD(1);
t1172 = t1192 + t1280;
t1233 = pkin(2) * t1172 + pkin(6) * t1169 + t980;
t1232 = t1259 + t1300;
t1171 = t1211 * qJDD(1) - t1209 * t1304;
t1231 = -pkin(5) * t1171 - t1209 * g(3);
t979 = t1210 * t1050 - t1208 * t1052;
t1224 = t1211 * t1177 - t1209 * t1178;
t833 = t1210 * t861 - t1293;
t850 = -pkin(3) * t1003 + qJ(4) * t861;
t1222 = -pkin(2) * t1003 + pkin(6) * t833 - qJ(4) * t1293 + t1210 * t850;
t1217 = pkin(4) * t1084 + qJ(5) * t1312 + t1226;
t1216 = t1217 + t1389;
t1140 = -0.2e1 * t1285;
t1215 = -t1150 * t1081 + t1140 - t1219 - t1320;
t1182 = t1192 - t1212;
t1180 = t1212 - t1280;
t1173 = -t1192 + t1280;
t1170 = t1209 * qJDD(1) + t1211 * t1304;
t1159 = t1210 * t1175;
t1158 = t1255 * t1290;
t1144 = -pkin(5) * t1170 + t1211 * g(3);
t1139 = 0.2e1 * t1286;
t1124 = t1210 * t1164 - t1198 * t1290;
t1123 = -t1208 * t1220 - t1302 * t1290;
t1122 = t1205 * t1158 + t1190;
t1121 = t1203 * t1158 - t1258;
t1118 = -t1208 * t1180 + t1159;
t1116 = t1210 * t1182 - t1265;
t1115 = t1210 * t1181 - t1265;
t1114 = t1210 * t1180 + t1266;
t1113 = t1208 * t1183 + t1159;
t1112 = t1208 * t1182 + t1261;
t1111 = (t1164 + t1188) * t1208;
t1107 = t1205 * t1169 - t1203 * t1172;
t1106 = t1203 * t1169 + t1205 * t1172;
t1103 = -t1208 * t1163 + t1110;
t1102 = t1210 * t1163 + t1208 * t1165;
t1080 = t1205 * t1124 - t1235;
t1079 = t1205 * t1123 + t1235;
t1078 = t1203 * t1124 + t1234;
t1077 = t1203 * t1123 - t1234;
t1075 = t1205 * t1118 + t1203 * t1257;
t1074 = t1205 * t1116 + t1203 * t1256;
t1073 = t1203 * t1118 - t1205 * t1257;
t1072 = t1203 * t1116 - t1205 * t1256;
t1058 = t1205 * t1119 + t1203 * t1163;
t1057 = t1205 * t1117 - t1203 * t1165;
t1056 = t1203 * t1119 - t1205 * t1163;
t1055 = t1203 * t1117 + t1205 * t1165;
t1054 = -pkin(1) * t1167 - t1099;
t1053 = pkin(1) * t1168 - t1098;
t1051 = t1205 * t1103 + t1203 * t1173;
t1049 = t1203 * t1103 - t1205 * t1173;
t1009 = -pkin(6) * t1115 + t1064;
t1008 = -pkin(6) * t1113 + t1063;
t1002 = pkin(1) * t1006;
t1001 = -pkin(2) * t1115 + t1052;
t1000 = -pkin(2) * t1113 + t1050;
t993 = pkin(1) * t1199 + qJ(2) * t1245;
t965 = t1210 * t1026 + t1208 * t1027;
t960 = pkin(1) * t1055 + t1249;
t959 = pkin(1) * t1056 + t1250;
t944 = -qJ(2) * t1106 + t1205 * t979;
t943 = qJ(2) * t1107 + t1203 * t979;
t940 = t1203 * t1068 + t1205 * t980;
t939 = -t1205 * t1068 + t1203 * t980;
t937 = pkin(1) * t1106 + t1233;
t903 = -qJ(2) * t1056 - t1203 * t1001 + t1205 * t1009;
t902 = -qJ(2) * t1055 - t1203 * t1000 + t1205 * t1008;
t884 = -pkin(1) * t1115 + qJ(2) * t1058 + t1205 * t1001 + t1203 * t1009;
t883 = -pkin(1) * t1113 + qJ(2) * t1057 + t1205 * t1000 + t1203 * t1008;
t868 = pkin(1) * t939 + t1291;
t865 = -t1300 - t1386;
t863 = -qJ(2) * t939 - (pkin(2) * t1203 - pkin(6) * t1205) * t979;
t859 = t1138 - t1236 - t1399;
t854 = t1368 + t920;
t853 = -t1232 - t1386;
t848 = -t1208 * t906 + t1210 * t942 - t1398;
t847 = qJ(2) * t940 - (-pkin(2) * t1205 - pkin(6) * t1203 - pkin(1)) * t979;
t844 = t1320 + t1368 + t880;
t841 = -t1208 * t900 + t1210 * t938 - t1372;
t838 = t1139 - t1216 - 0.2e1 * t1284 + t1399;
t836 = t1238 + t1409;
t832 = t1208 * t861 + t1292;
t829 = t1239 + t1387;
t828 = t1203 * t1003 + t1205 * t833;
t827 = -t1205 * t1003 + t1203 * t833;
t824 = -t1208 * t856 + t1210 * t864 - t1372;
t823 = -t1208 * t852 + t1210 * t858 + t1398;
t821 = t1208 * t843 + t1210 * t842;
t819 = -t1203 * t859 + t1205 * t848 - t1408;
t816 = -pkin(2) * t832 - t1301;
t815 = t1240 + t1387;
t814 = t1203 * t848 + t1205 * t859 - t1407;
t813 = -t1203 * t854 + t1205 * t841 - t1383;
t812 = t1203 * t905 + t1205 * t822;
t811 = t1203 * t822 - t1205 * t905;
t810 = -t1208 * t840 + t1210 * t846 - t1385;
t809 = t1241 - t1409;
t808 = t1203 * t841 + t1205 * t854 + t1380;
t807 = t1242 + t1396;
t806 = -pkin(6) * t832 - qJ(4) * t1292 - t1208 * t850;
t805 = -t1203 * t844 + t1205 * t824 - t1383;
t804 = -t1208 * t834 + t1210 * t837 - t1385;
t803 = t1203 * t824 + t1205 * t844 + t1380;
t802 = -t1203 * t838 + t1205 * t823 + t1408;
t801 = -t1203 * t865 + t1205 * t810 - t1395;
t800 = t1203 * t823 + t1205 * t838 + t1407;
t799 = -pkin(2) * t821 - t1237;
t798 = t1243 + t1396;
t797 = t1203 * t810 + t1205 * t865 + t1394;
t796 = -t1203 * t853 + t1205 * t804 - t1395;
t795 = pkin(1) * t827 + t1222;
t794 = t1203 * t804 + t1205 * t853 + t1394;
t793 = -pkin(6) * t821 - t1208 * t818 + t1210 * t826;
t792 = -qJ(2) * t827 - t1203 * t816 + t1205 * t806;
t791 = -pkin(1) * t832 + qJ(2) * t828 + t1203 * t806 + t1205 * t816;
t790 = pkin(1) * t811 + t1246;
t789 = -qJ(2) * t811 - t1203 * t799 + t1205 * t793;
t788 = -pkin(1) * t821 + qJ(2) * t812 + t1203 * t793 + t1205 * t799;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1171, 0, -t1170, 0, t1231, -t1144, -t1224, -pkin(5) * t1224, 0, 0, t1104, 0, -t1317, 0, t1325, t1331, t1350, pkin(5) * t1350 + qJ(2) * t1260 - t1209 * t993, -t1209 * t1078 + t1211 * t1080, -t1209 * t1049 + t1211 * t1051, -t1209 * t1073 + t1211 * t1075, -t1209 * t1077 + t1211 * t1079, -t1209 * t1072 + t1211 * t1074, -t1209 * t1121 + t1211 * t1122, t1211 * t902 - t1209 * t883 - pkin(5) * (t1211 * t1055 + t1209 * t1057), t1211 * t903 - t1209 * t884 - pkin(5) * (t1211 * t1056 + t1209 * t1058), t1211 * t944 - t1209 * t943 - pkin(5) * (t1211 * t1106 + t1209 * t1107), t1211 * t863 - t1209 * t847 - pkin(5) * (t1209 * t940 + t1211 * t939), t1305, t1405, t1393, t1344, t1403, t1341, -t1209 * t808 + t1211 * t813 - t1374, -t1209 * t814 + t1211 * t819 - t1410, -t1209 * t797 + t1211 * t801 - t1390, t1211 * t792 - t1209 * t791 - pkin(5) * (t1209 * t828 + t1211 * t827), t1305, t1393, -t1405, t1341, -t1403, t1344, -t1209 * t803 + t1211 * t805 - t1374, -t1209 * t794 + t1211 * t796 - t1390, -t1209 * t800 + t1211 * t802 + t1410, t1211 * t789 - t1209 * t788 - pkin(5) * (t1209 * t812 + t1211 * t811); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1170, 0, t1171, 0, t1144, t1231, t1244, pkin(5) * t1244, 0, 0, t1317, 0, t1104, 0, -t1331, t1325, t1349, pkin(5) * t1349 + qJ(2) * t1264 + t1211 * t993, t1211 * t1078 + t1209 * t1080, t1211 * t1049 + t1209 * t1051, t1211 * t1073 + t1209 * t1075, t1211 * t1077 + t1209 * t1079, t1211 * t1072 + t1209 * t1074, t1211 * t1121 + t1209 * t1122, t1209 * t902 + t1211 * t883 + pkin(5) * (-t1209 * t1055 + t1211 * t1057), t1209 * t903 + t1211 * t884 + pkin(5) * (-t1209 * t1056 + t1211 * t1058), t1209 * t944 + t1211 * t943 + pkin(5) * (-t1209 * t1106 + t1211 * t1107), t1209 * t863 + t1211 * t847 + pkin(5) * (-t1209 * t939 + t1211 * t940), t1306, -t1406, t1392, t1343, -t1404, t1342, t1209 * t813 + t1211 * t808 + t1375, t1209 * t819 + t1211 * t814 - t1411, t1209 * t801 + t1211 * t797 + t1391, t1209 * t792 + t1211 * t791 + pkin(5) * (-t1209 * t827 + t1211 * t828), t1306, t1392, t1406, t1342, t1404, t1343, t1209 * t805 + t1211 * t803 + t1375, t1209 * t796 + t1211 * t794 + t1391, t1209 * t802 + t1211 * t800 + t1411, t1209 * t789 + t1211 * t788 + pkin(5) * (-t1209 * t811 + t1211 * t812); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1177, t1178, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1053, t1054, 0, -t1002, t1111, t1102, t1114, t1110, t1112, 0, t960, t959, t937, t868, t965, -t891, t1361, t1307, t947, t1310, t829, t836, t807, t795, t965, t1361, t891, t1310, -t947, t1307, t815, t798, t809, t790; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1304, 0, 0, -g(3), -t1177, 0, 0, 0, t1168, 0, -t1167, 0, t1225, t1130, t1006, qJ(2) * t1006, t1080, t1051, t1075, t1079, t1074, t1122, t902, t903, t944, t863, t1230, -t875, t1377, t1326, -t915, t1330, t813, t819, t801, t792, t1230, t1377, t875, t1330, t915, t1326, t805, t796, t802, t789; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1304, 0, qJDD(1), 0, g(3), 0, -t1178, 0, 0, 0, t1167, 0, t1168, 0, -t1130, t1225, t1245, t993, t1078, t1049, t1073, t1077, t1072, t1121, t883, t884, t943, t847, t1229, -t873, t1379, t1329, -t909, t1328, t808, t814, t797, t791, t1229, t1379, t873, t1328, t909, t1329, t803, t794, t800, t788; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1177, t1178, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1053, t1054, 0, -t1002, t1111, t1102, t1114, t1110, t1112, 0, t960, t959, t937, t868, t965, -t891, t1361, t1307, t947, t1310, t829, t836, t807, t795, t965, t1361, t891, t1310, -t947, t1307, t815, t798, t809, t790; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1304, 0, 0, -t1199, t1098, 0, t1124, t1103, t1118, t1123, t1116, t1158, t1008, t1009, t979, pkin(6) * t979, t968, -t895, t1360, t1308, -t954, t1309, t841, t848, t810, t806, t968, t1360, t895, t1309, t954, t1308, t824, t804, t823, t793; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1304, 0, qJDD(1), 0, t1199, 0, t1099, 0, t1184, -t1173, -t1257, -t1184, -t1256, -qJDD(3), t1000, t1001, 0, pkin(2) * t979, -t1282, -t1318, -t1315, t1282, t1043, -qJDD(3), t854, t859, t865, t816, -t1282, -t1315, t1318, -qJDD(3), -t1043, t1282, t844, t853, t838, t799; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1098, -t1099, 0, 0, t1111, t1102, t1114, t1110, t1112, 0, t1249, t1250, t1233, t1291, t965, -t891, t1361, t1307, t947, t1310, t1239, t1238, t1242, t1222, t965, t1361, t891, t1310, -t947, t1307, t1240, t1243, t1241, t1246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1164, t1165, t1175, -t1188, t1182, t1188, 0, t1068, t1050, 0, t1027, -t973, t1352, t1221, t1016, t1227, t938, t942, t846, -qJ(4) * t860, t1027, t1352, t973, t1227, -t1016, t1221, t864, t837, t858, t826; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1251, t1163, t1180, t1220, t1176, -t1251, -t1068, 0, t1052, 0, t1026, t971, t1351, t1228, t1010, t1218, t900, t906, t840, t850, t1026, t1351, -t971, t1218, -t1010, t1228, t856, t834, t852, t818; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1184, t1173, t1257, t1184, t1256, qJDD(3), -t1050, -t1052, 0, 0, t1282, t1318, t1315, -t1282, -t1043, qJDD(3), t1140 - t1247 + t1357, t1139 + t1236, t1300, t1301, t1282, t1315, -t1318, qJDD(3), t1043, -t1282, t1215 + t1357, t1232, t1216 + t1319, t1237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1101, -t1042, t1313, t1287, t1131, -t1287, 0, t1003, t920, 0, t1101, t1313, t1042, -t1287, -t1131, t1287, -qJ(5) * t1042, t867, t881, -qJ(5) * t905; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1143, t1316, -t1132, -t1100, t1312, -t1143, -t1003, 0, t921, 0, t1143, -t1132, -t1316, -t1143, -t1312, -t1100, t882, t866, pkin(4) * t1316, -pkin(4) * t905; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1282, t1318, t1315, -t1282, -t1043, qJDD(3), -t920, -t921, 0, 0, t1282, t1315, -t1318, qJDD(3), t1043, -t1282, t1215, t1259, t1217 + t1319, t1295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1101, t1313, t1042, -t1287, -t1131, t1287, 0, t880, -t905, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1282, t1315, -t1318, qJDD(3), t1043, -t1282, -t880, 0, t878, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1143, t1132, t1316, t1143, t1312, t1100, t905, -t878, 0, 0;];
m_new_reg = t1;
