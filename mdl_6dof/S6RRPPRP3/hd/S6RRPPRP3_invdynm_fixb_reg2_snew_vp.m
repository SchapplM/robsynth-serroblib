% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RRPPRP3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP3_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:16:58
% EndTime: 2019-05-06 09:17:13
% DurationCPUTime: 15.52s
% Computational Cost: add. (38783->742), mult. (82085->699), div. (0->0), fcn. (43480->6), ass. (0->461)
t1189 = sin(qJ(5));
t1192 = cos(qJ(5));
t1193 = cos(qJ(2));
t1325 = qJD(1) * t1193;
t1134 = -t1192 * qJD(2) + t1189 * t1325;
t1135 = qJD(2) * t1189 + t1192 * t1325;
t1078 = t1134 * t1135;
t1324 = qJD(2) * t1193;
t1172 = qJD(1) * t1324;
t1190 = sin(qJ(2));
t1283 = t1190 * qJDD(1);
t1142 = t1172 + t1283;
t1126 = qJDD(5) + t1142;
t1412 = t1078 - t1126;
t1420 = pkin(5) * t1412;
t1184 = t1190 ^ 2;
t1195 = qJD(1) ^ 2;
t1297 = t1184 * t1195;
t1385 = qJD(2) ^ 2;
t1159 = t1297 + t1385;
t1399 = t1193 * t1195;
t1272 = t1190 * t1399;
t1154 = -qJDD(2) + t1272;
t1292 = t1193 * t1154;
t1093 = t1159 * t1190 + t1292;
t1141 = 0.2e1 * t1172 + t1283;
t1191 = sin(qJ(1));
t1194 = cos(qJ(1));
t1036 = pkin(6) * (t1093 * t1194 + t1141 * t1191);
t1365 = pkin(6) * (t1093 * t1191 - t1141 * t1194);
t1153 = qJDD(2) + t1272;
t1131 = t1193 * t1153;
t1384 = t1193 ^ 2;
t1291 = t1384 * t1195;
t1161 = t1291 + t1385;
t1091 = -t1161 * t1190 + t1131;
t1369 = pkin(1) * t1091;
t1304 = t1153 * t1190;
t1095 = t1161 * t1193 + t1304;
t1294 = t1190 * qJD(1);
t1171 = qJD(2) * t1294;
t1282 = t1193 * qJDD(1);
t1144 = -0.2e1 * t1171 + t1282;
t1037 = pkin(6) * (t1095 * t1194 + t1144 * t1191);
t1364 = pkin(6) * (t1095 * t1191 - t1144 * t1194);
t1361 = pkin(7) * t1091;
t1143 = -t1171 + t1282;
t1064 = qJD(5) * t1134 - t1189 * qJDD(2) - t1143 * t1192;
t1165 = qJD(5) + t1294;
t1110 = t1165 * t1134;
t1155 = t1191 * g(1) - t1194 * g(2);
t1118 = qJDD(1) * pkin(1) + t1195 * pkin(7) + t1155;
t1213 = -pkin(2) * t1171 + t1118;
t1392 = t1143 * pkin(3) - qJ(4) * t1291 + qJDD(4);
t1205 = t1213 + t1392;
t1152 = -qJD(2) * pkin(3) - qJ(4) * t1294;
t1290 = (2 * qJD(3)) + t1152;
t1266 = t1290 * t1190;
t1371 = pkin(4) + qJ(3);
t924 = (pkin(2) + pkin(8)) * t1143 + t1371 * t1142 + (t1266 + (-t1190 * pkin(8) + t1193 * t1371) * qJD(2)) * qJD(1) + t1205;
t1359 = t1193 * g(3);
t1230 = -qJDD(2) * pkin(2) - t1385 * qJ(3) + qJDD(3) + t1359;
t1156 = g(1) * t1194 + g(2) * t1191;
t1119 = -pkin(1) * t1195 + qJDD(1) * pkin(7) - t1156;
t1293 = t1190 * t1119;
t1414 = pkin(3) * t1153;
t1206 = -t1142 * qJ(4) + t1230 + t1293 - t1414;
t1249 = pkin(4) * t1190 + pkin(8) * t1193;
t1275 = qJ(4) * t1324;
t1245 = -pkin(2) * t1193 - qJ(3) * t1190;
t1139 = t1245 * qJD(1);
t1289 = -(2 * qJD(4)) + t1139;
t956 = -t1385 * pkin(4) - qJDD(2) * pkin(8) + (t1275 + (-qJD(1) * t1249 + t1289) * t1190) * qJD(1) + t1206;
t878 = t1189 * t956 - t1192 * t924;
t1211 = -qJ(6) * t1110 - 0.2e1 * qJD(6) * t1135 + t1420 + t878;
t864 = qJ(6) * t1064 + t1211;
t1207 = -t864 - t1420;
t1124 = t1134 ^ 2;
t1163 = t1165 ^ 2;
t1065 = -t1163 - t1124;
t1318 = t1412 * t1192;
t978 = t1065 * t1189 - t1318;
t973 = pkin(4) * t978;
t1419 = -t973 - t1207;
t1125 = t1135 ^ 2;
t1280 = t1125 + t1163;
t1233 = -t1192 * qJDD(2) + t1189 * t1143;
t1063 = qJD(5) * t1135 + t1233;
t1099 = pkin(5) * t1165 + qJ(6) * t1135;
t879 = t1189 * t924 + t1192 * t956;
t866 = -t1124 * pkin(5) + t1063 * qJ(6) + 0.2e1 * qJD(6) * t1134 - t1165 * t1099 + t879;
t1216 = -pkin(5) * t1280 - t866;
t1049 = t1078 + t1126;
t1296 = t1189 * t1049;
t990 = -t1192 * t1280 - t1296;
t983 = pkin(4) * t990;
t1418 = -t983 - t1216;
t1288 = pkin(1) * t1141 - pkin(7) * t1093;
t1303 = t1154 * t1190;
t1085 = -t1159 * t1193 + t1303;
t1415 = pkin(1) * t1085;
t1362 = pkin(7) * t1085;
t1295 = t1189 * t1412;
t1263 = qJD(1) * t1139 + t1119;
t1029 = t1190 * t1263 + t1230;
t1342 = qJ(3) * t1161;
t1413 = pkin(2) * t1153 - t1029 - t1342;
t1237 = t1263 * t1193;
t1181 = t1190 * g(3);
t1284 = t1385 * pkin(2) + t1181;
t1320 = qJDD(2) * qJ(3);
t1411 = t1237 - t1284 + t1320;
t1287 = pkin(1) * t1144 - pkin(7) * t1095;
t1308 = t1144 * t1193;
t1314 = t1141 * t1190;
t1071 = -t1308 + t1314;
t1151 = (t1184 - t1384) * t1195;
t1410 = t1071 * t1191 + t1151 * t1194;
t1409 = t1071 * t1194 - t1151 * t1191;
t1162 = t1291 - t1385;
t1096 = t1162 * t1193 + t1303;
t1056 = t1096 * t1191 - t1194 * t1282;
t1058 = t1096 * t1194 + t1191 * t1282;
t1279 = t1184 + t1384;
t1147 = t1279 * qJDD(1);
t1150 = t1279 * t1195;
t1073 = pkin(6) * (t1147 * t1194 - t1150 * t1191);
t1387 = -t1063 * pkin(5) - t1124 * qJ(6) - t1135 * t1099 + qJDD(6);
t1229 = pkin(3) * t1291 + t1143 * qJ(4) + t1284;
t1210 = qJD(2) * t1290 - t1229;
t1236 = qJD(1) * t1289 + t1119;
t1201 = t1193 * t1236 + t1210;
t1386 = -t1385 * pkin(8) + qJDD(2) * t1371 - t1249 * t1399;
t955 = t1201 + t1386;
t898 = t955 + t1387;
t886 = qJ(6) * t1280 + t898;
t1398 = t1064 + t1110;
t961 = -pkin(5) * t1398 - qJ(6) * t1049;
t1406 = t1189 * t886 + t1192 * t961;
t1111 = t1165 * t1135;
t1017 = -t1063 - t1111;
t1405 = -pkin(3) * t978 + qJ(4) * t1017;
t1404 = -pkin(3) * t990 + qJ(4) * t1398;
t979 = t1065 * t1192 + t1295;
t1403 = -pkin(2) * t979 + qJ(3) * t1017;
t1047 = -t1124 - t1125;
t1019 = (qJD(5) - t1165) * t1135 + t1233;
t1022 = t1064 - t1110;
t945 = t1019 * t1189 - t1022 * t1192;
t1402 = -pkin(3) * t945 + qJ(4) * t1047;
t947 = t1019 * t1192 + t1189 * t1022;
t1401 = -pkin(2) * t947 + qJ(3) * t1047;
t849 = t1189 * t878 + t1192 * t879;
t1323 = qJD(3) * qJD(2);
t1180 = -0.2e1 * t1323;
t1209 = -qJD(2) * t1152 + 0.2e1 * qJD(4) * t1325 + t1180 + t1229;
t1202 = -t1237 + t1209;
t869 = -pkin(5) * t1017 + qJ(6) * t1065 + t1202 - t1386 - t1387;
t1400 = qJ(6) * t1295 + t1192 * t869;
t848 = t1189 * t879 - t1192 * t878;
t1286 = pkin(1) * t1150 + pkin(7) * t1147;
t1396 = pkin(2) * t1159 - qJ(3) * t1154;
t1090 = -t1162 * t1190 + t1292;
t1395 = pkin(3) * t1144 + qJ(4) * t1161;
t1393 = -t1396 + t1415;
t1244 = t1142 + t1172;
t1281 = 0.2e1 * t1294;
t1265 = qJD(3) * t1281;
t1199 = qJ(3) * t1244 + t1213 + t1265;
t989 = (t1143 + t1144) * pkin(2) + t1199;
t1336 = t1189 * t864;
t834 = t1192 * t866 + t1336;
t845 = -pkin(5) * t898 + qJ(6) * t866;
t1264 = -pkin(4) * t898 + pkin(8) * t834 + qJ(6) * t1336 + t1192 * t845;
t1383 = pkin(2) + pkin(3);
t1391 = qJ(3) * t898 - t1383 * t834 - t1264;
t1345 = pkin(4) * t955 - pkin(8) * t849;
t1390 = qJ(3) * t955 - t1383 * t849 + t1345;
t969 = t1201 + t1320;
t970 = (t1190 * t1289 + t1275) * qJD(1) + t1206;
t1389 = qJ(3) * t969 - t1383 * t970;
t1388 = -t1153 * t1383 + t1342;
t1382 = -pkin(3) - pkin(8);
t905 = t1047 * t1190 - t1193 * t947;
t1381 = pkin(1) * t905;
t1380 = pkin(2) * t945;
t1379 = pkin(2) * t978;
t1378 = pkin(2) * t990;
t906 = t1047 * t1193 + t1190 * t947;
t1377 = pkin(6) * (t1191 * t906 + t1194 * t945);
t917 = t1017 * t1193 + t1190 * t979;
t1376 = pkin(6) * (t1191 * t917 + t1194 * t978);
t1319 = t1049 * t1192;
t991 = t1189 * t1280 - t1319;
t923 = t1190 * t991 + t1193 * t1398;
t1375 = pkin(6) * (t1191 * t923 + t1194 * t990);
t1374 = pkin(7) * t905;
t916 = t1017 * t1190 - t1193 * t979;
t1373 = pkin(7) * t916;
t922 = t1190 * t1398 - t1193 * t991;
t1372 = pkin(7) * t922;
t1366 = pkin(3) * t1159;
t1363 = pkin(6) * (t1147 * t1191 + t1150 * t1194);
t1360 = t1143 * pkin(2);
t1357 = qJ(3) * t978;
t1356 = qJ(3) * t990;
t1355 = qJ(4) * t955;
t1354 = qJ(4) * t969;
t1353 = qJ(4) * t970;
t1352 = qJ(4) * t979;
t1351 = qJ(4) * t991;
t1334 = t1192 * t864;
t833 = t1189 * t866 - t1334;
t863 = pkin(5) * t864;
t1350 = -pkin(4) * t833 + t863;
t1349 = pkin(1) * t945 + pkin(7) * t906;
t1348 = pkin(1) * t978 + pkin(7) * t917;
t1347 = pkin(1) * t990 + pkin(7) * t923;
t1332 = t1192 * t955;
t981 = pkin(8) * t990;
t1346 = t1332 - t981;
t1344 = qJ(3) * t1144;
t1343 = qJ(3) * t1150;
t1341 = qJ(4) * t1153;
t1340 = qJ(4) * t1154;
t1339 = qJ(4) * t1159;
t950 = t1189 * t955;
t1331 = -pkin(4) * t1047 + pkin(8) * t947;
t1330 = -pkin(4) * t1017 + pkin(8) * t979;
t1329 = -pkin(4) * t1398 + pkin(8) * t991;
t1328 = pkin(2) - t1382;
t1327 = qJ(4) * qJDD(1);
t1326 = qJD(1) * qJD(2);
t1317 = t1118 * t1190;
t1316 = t1118 * t1193;
t1312 = t1141 * t1193;
t1310 = t1144 * t1190;
t1299 = t1165 * t1189;
t1298 = t1165 * t1192;
t1278 = -t983 + t879;
t1277 = t950 + t1329;
t1276 = -t1332 + t1330;
t1274 = t1190 * t1078;
t1273 = t1193 * t1078;
t847 = pkin(4) * t848;
t1269 = -qJ(4) * t849 + t847;
t940 = pkin(4) * t945;
t1268 = -qJ(4) * t947 + t940;
t971 = pkin(8) * t978;
t1267 = -t971 + t950;
t1100 = t1293 + t1359;
t1101 = t1193 * t1119 - t1181;
t1025 = t1100 * t1190 + t1193 * t1101;
t1262 = -t1155 * t1191 - t1194 * t1156;
t853 = -pkin(5) * t1047 + qJ(6) * t1019 + t866;
t857 = (t1022 + t1064) * qJ(6) + t1211;
t1261 = t1189 * t857 + t1192 * t853 + t1331;
t1260 = t1329 + t1406;
t1259 = t1346 + t1404;
t1258 = t1330 + t1400;
t1257 = t1191 * t1272;
t1256 = t1194 * t1272;
t1254 = -qJ(4) * t834 - t1350;
t938 = pkin(8) * t945;
t1253 = -t1189 * t853 + t1192 * t857 - t938;
t1252 = -t938 - t848;
t1251 = -t1189 * t961 + t1192 * t886 - t981;
t1250 = t878 - t973;
t1149 = qJDD(1) * t1194 - t1191 * t1195;
t1248 = -pkin(6) * t1149 - g(3) * t1191;
t942 = pkin(3) * t947;
t1247 = -t942 - t1261;
t1026 = 0.2e1 * t1323 + t1411;
t1246 = -pkin(2) * t1029 + qJ(3) * t1026;
t975 = pkin(3) * t979;
t1242 = -t975 - t1330 + t1403;
t985 = pkin(3) * t991;
t1241 = -pkin(2) * t991 + qJ(3) * t1398 - t1329 - t985;
t1240 = -t975 - t1258;
t1238 = -t1278 - t1351;
t1024 = t1100 * t1193 - t1101 * t1190;
t1068 = t1310 + t1312;
t1234 = t1155 * t1194 - t1156 * t1191;
t1232 = t1267 + t1405;
t1231 = t1331 + t849;
t865 = qJ(3) * t945 + t1268;
t1228 = qJ(6) * t1334 - t1189 * t845;
t1226 = t1251 + t1404;
t1225 = t1253 + t1402;
t1224 = t1252 + t1402;
t1223 = -t1231 - t942;
t1221 = t1247 + t1401;
t1220 = -t1250 - t1352;
t1219 = qJ(6) * t1318 - t1189 * t869 - t971;
t1218 = -pkin(1) * t916 - t1242;
t1217 = -pkin(1) * t922 - t1241;
t1215 = -qJ(4) * t898 - t1228;
t1214 = t1223 + t1401;
t1212 = t1219 + t1405;
t1208 = -t1351 - t1418;
t1203 = -t1352 - t1419;
t1200 = t1360 + (t1141 + t1244) * qJ(3) + t1213;
t1198 = qJD(4) * t1281 + (-t1190 * t1139 - t1275) * qJD(1) - t1206;
t1197 = qJ(4) * t1283 + t1198;
t1196 = t969 + t1366;
t962 = t1360 + t1142 * qJ(3) + (qJ(3) * t1324 + t1266) * qJD(1) + t1205;
t1174 = qJ(3) * t1282;
t1160 = t1297 - t1385;
t1148 = qJDD(1) * t1191 + t1194 * t1195;
t1137 = pkin(2) * t1283 - t1174;
t1130 = t1279 * t1326;
t1113 = -pkin(6) * t1148 + g(3) * t1194;
t1112 = t1283 * t1383 - t1174;
t1107 = -t1125 + t1163;
t1106 = t1124 - t1163;
t1105 = qJDD(2) * t1191 + t1130 * t1194;
t1104 = t1142 * t1193 - t1184 * t1326;
t1103 = -qJDD(2) * t1194 + t1130 * t1191;
t1102 = -t1143 * t1190 - t1326 * t1384;
t1092 = t1160 * t1190 + t1131;
t1087 = t1244 * t1190;
t1084 = -t1160 * t1193 + t1304;
t1083 = (t1143 - t1171) * t1193;
t1081 = -t1341 - t1344;
t1075 = t1125 - t1124;
t1062 = t1104 * t1194 - t1257;
t1061 = t1102 * t1194 + t1257;
t1060 = t1104 * t1191 + t1256;
t1059 = t1102 * t1191 - t1256;
t1057 = t1092 * t1194 + t1191 * t1283;
t1055 = t1092 * t1191 - t1194 * t1283;
t1051 = t1141 * t1383 + t1340;
t1035 = -t1316 - t1362;
t1034 = -t1317 - t1361;
t1033 = (t1134 * t1192 - t1135 * t1189) * t1165;
t1032 = (t1134 * t1189 + t1135 * t1192) * t1165;
t1028 = t1100 - t1369;
t1027 = t1101 - t1415;
t1018 = -t1063 + t1111;
t1014 = pkin(5) * t1022;
t1009 = t1287 + t1316;
t1008 = -t1288 - t1317;
t1005 = t1029 + t1343;
t1004 = pkin(2) * t1150 + t1026;
t1003 = t1064 * t1192 + t1135 * t1299;
t1002 = t1064 * t1189 - t1135 * t1298;
t1001 = t1063 * t1189 + t1134 * t1298;
t1000 = -t1192 * t1063 + t1134 * t1299;
t999 = t1199 + t1360;
t998 = t1033 * t1190 + t1126 * t1193;
t997 = -t1033 * t1193 + t1126 * t1190;
t996 = t1106 * t1192 - t1296;
t995 = -t1189 * t1107 - t1318;
t994 = t1106 * t1189 + t1319;
t993 = t1107 * t1192 - t1295;
t992 = pkin(1) * t1118 + pkin(7) * t1025;
t988 = t1200 + t1265;
t980 = t1025 + t1286;
t968 = -t1369 - t1413;
t967 = t1003 * t1190 + t1273;
t966 = -t1001 * t1190 - t1273;
t965 = -t1003 * t1193 + t1274;
t964 = t1001 * t1193 - t1274;
t963 = t1180 + t1393 - t1411;
t959 = t1197 - t1343;
t958 = t1026 * t1193 + t1029 * t1190;
t957 = t1026 * t1190 - t1029 * t1193;
t949 = -t1320 - t1383 * t1150 + (-t1263 + t1327) * t1193 + t1209;
t948 = -t1017 * t1192 - t1189 * t1398;
t946 = -t1017 * t1189 + t1192 * t1398;
t937 = qJD(1) * t1266 + t1200 + t1339 + t1392;
t936 = -t1032 * t1191 + t1194 * t998;
t935 = t1032 * t1194 + t1191 * t998;
t934 = -pkin(2) * t1314 + t1193 * t988 + t1362;
t933 = qJ(3) * t1308 - t1190 * t989 - t1361;
t932 = -t1004 * t1190 + t1005 * t1193;
t931 = t1022 * t1193 + t1190 * t995;
t930 = -t1018 * t1193 + t1190 * t996;
t929 = t1022 * t1190 - t1193 * t995;
t928 = -t1018 * t1190 - t1193 * t996;
t927 = -t1152 * t1294 - t1392 - t1395 - t989;
t926 = pkin(2) * t1312 + t1190 * t988 + t1288;
t925 = qJ(3) * t1310 + t1193 * t989 + t1287;
t918 = t1202 - t1320 - t1366 + t1393;
t915 = t1198 - t1388 + t1369;
t913 = t1004 * t1193 + t1005 * t1190 + t1286;
t912 = t1075 * t1193 + t1190 * t948;
t911 = t1075 * t1190 - t1193 * t948;
t910 = -t1002 * t1191 + t1194 * t967;
t909 = t1000 * t1191 + t1194 * t966;
t908 = t1002 * t1194 + t1191 * t967;
t907 = -t1000 * t1194 + t1191 * t966;
t903 = t1190 * t970 + t1193 * t969;
t902 = t1190 * t969 - t1193 * t970;
t901 = qJ(3) * t962 - t1353;
t900 = t1081 * t1193 - t1190 * t927 + t1361;
t899 = -t1051 * t1190 + t1193 * t937 + t1362;
t896 = t1081 * t1190 + t1193 * t927 - t1287;
t895 = t1051 * t1193 + t1190 * t937 + t1288;
t894 = -t1191 * t993 + t1194 * t931;
t893 = -t1191 * t994 + t1194 * t930;
t892 = t1191 * t931 + t1194 * t993;
t891 = t1191 * t930 + t1194 * t994;
t890 = -pkin(1) * t957 - t1246;
t888 = pkin(6) * (-t1191 * t990 + t1194 * t923);
t887 = -t1190 * t949 + t1193 * t959;
t882 = pkin(6) * (-t1191 * t978 + t1194 * t917);
t881 = -pkin(7) * t957 + (-pkin(2) * t1190 + qJ(3) * t1193) * t999;
t880 = t1190 * t959 + t1193 * t949 - t1286;
t874 = -t1191 * t946 + t1194 * t912;
t873 = t1191 * t912 + t1194 * t946;
t872 = t1383 * t962 - t1354;
t870 = pkin(6) * (-t1191 * t945 + t1194 * t906);
t867 = pkin(7) * t958 + (pkin(1) - t1245) * t999;
t860 = -t1259 + t1378;
t859 = -t1014 + t865;
t858 = -t1232 + t1379;
t854 = -pkin(1) * t902 - t1389;
t851 = t1238 + t1356;
t850 = t1220 + t1357;
t843 = t1217 + t950;
t842 = t1208 + t1356;
t841 = t1218 - t1332;
t840 = -t1226 + t1378;
t839 = t1190 * t849 + t1193 * t955;
t838 = t1190 * t955 - t1193 * t849;
t837 = -t1212 + t1379;
t836 = t1203 + t1357;
t835 = -pkin(7) * t902 - t1190 * t872 + t1193 * t901;
t830 = pkin(1) * t962 + pkin(7) * t903 + t1190 * t901 + t1193 * t872;
t829 = t1217 + t1406;
t828 = t1218 + t1400;
t827 = -t1224 + t1380;
t826 = t1190 * t834 + t1193 * t898;
t825 = t1190 * t898 - t1193 * t834;
t824 = -t1190 * t860 + t1193 * t851 - t1372;
t823 = -t1190 * t858 + t1193 * t850 - t1373;
t822 = t1190 * t851 + t1193 * t860 + t1347;
t821 = -t1214 - t1381;
t820 = t1190 * t850 + t1193 * t858 + t1348;
t819 = -t1225 + t1380;
t818 = -t1190 * t840 + t1193 * t842 - t1372;
t817 = -t1190 * t827 + t1193 * t865 - t1374;
t816 = -t1190 * t837 + t1193 * t836 - t1373;
t815 = t1190 * t842 + t1193 * t840 + t1347;
t814 = -t1221 - t1381;
t813 = t1190 * t836 + t1193 * t837 + t1348;
t812 = qJ(3) * t848 + t1269;
t811 = t1190 * t865 + t1193 * t827 + t1349;
t810 = t1328 * t848 - t1355;
t809 = -t1190 * t819 + t1193 * t859 - t1374;
t808 = t1190 * t859 + t1193 * t819 + t1349;
t807 = -pkin(1) * t838 - t1390;
t806 = qJ(3) * t833 + t1254;
t805 = t1328 * t833 + t1215;
t804 = -pkin(7) * t838 - t1190 * t810 + t1193 * t812;
t803 = pkin(1) * t848 + pkin(7) * t839 + t1190 * t812 + t1193 * t810;
t802 = -pkin(1) * t825 - t1391;
t801 = -pkin(7) * t825 - t1190 * t805 + t1193 * t806;
t800 = pkin(1) * t833 + pkin(7) * t826 + t1190 * t806 + t1193 * t805;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1149, 0, -t1148, 0, t1248, -t1113, -t1234, -pkin(6) * t1234, t1062, -t1409, t1057, t1061, t1058, t1105, -t1028 * t1191 + t1034 * t1194 + t1364, -t1191 * t1027 + t1194 * t1035 - t1365, t1024 * t1194 - t1363, -pkin(6) * (t1025 * t1191 + t1118 * t1194) - (pkin(1) * t1191 - pkin(7) * t1194) * t1024, t1062, t1057, t1409, t1105, -t1058, t1061, -t1191 * t968 + t1194 * t933 + t1364, -t1137 * t1191 + t1194 * t932 - t1363, -t1191 * t963 + t1194 * t934 + t1365, t1194 * t881 - t1191 * t890 - pkin(6) * (t1191 * t958 + t1194 * t999), t1061, -t1409, t1058, t1062, t1057, t1105, -t1191 * t918 + t1194 * t899 + t1365, -t1191 * t915 + t1194 * t900 - t1364, t1191 * t1112 + t1194 * t887 + t1363, t1194 * t835 - t1191 * t854 - pkin(6) * (t1191 * t903 + t1194 * t962), t910, t874, t894, t909, t893, t936, -t1191 * t841 + t1194 * t823 - t1376, -t1191 * t843 + t1194 * t824 - t1375, -t1191 * t821 + t1194 * t817 - t1377, t1194 * t804 - t1191 * t807 - pkin(6) * (t1191 * t839 + t1194 * t848), t910, t874, t894, t909, t893, t936, -t1191 * t828 + t1194 * t816 - t1376, -t1191 * t829 + t1194 * t818 - t1375, -t1191 * t814 + t1194 * t809 - t1377, t1194 * t801 - t1191 * t802 - pkin(6) * (t1191 * t826 + t1194 * t833); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1148, 0, t1149, 0, t1113, t1248, t1262, pkin(6) * t1262, t1060, -t1410, t1055, t1059, t1056, t1103, t1028 * t1194 + t1034 * t1191 - t1037, t1194 * t1027 + t1191 * t1035 + t1036, t1024 * t1191 + t1073, pkin(6) * (t1025 * t1194 - t1118 * t1191) - (-pkin(1) * t1194 - pkin(7) * t1191) * t1024, t1060, t1055, t1410, t1103, -t1056, t1059, t1191 * t933 + t1194 * t968 - t1037, t1137 * t1194 + t1191 * t932 + t1073, t1191 * t934 + t1194 * t963 - t1036, t1191 * t881 + t1194 * t890 + pkin(6) * (-t1191 * t999 + t1194 * t958), t1059, -t1410, t1056, t1060, t1055, t1103, t1191 * t899 + t1194 * t918 - t1036, t1191 * t900 + t1194 * t915 + t1037, -t1194 * t1112 + t1191 * t887 - t1073, t1191 * t835 + t1194 * t854 + pkin(6) * (-t1191 * t962 + t1194 * t903), t908, t873, t892, t907, t891, t935, t1191 * t823 + t1194 * t841 + t882, t1191 * t824 + t1194 * t843 + t888, t1191 * t817 + t1194 * t821 + t870, t1191 * t804 + t1194 * t807 + pkin(6) * (-t1191 * t848 + t1194 * t839), t908, t873, t892, t907, t891, t935, t1191 * t816 + t1194 * t828 + t882, t1191 * t818 + t1194 * t829 + t888, t1191 * t809 + t1194 * t814 + t870, t1191 * t801 + t1194 * t802 + pkin(6) * (-t1191 * t833 + t1194 * t826); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1155, t1156, 0, 0, t1087, t1068, t1084, t1083, -t1090, 0, t1009, t1008, t980, t992, t1087, t1084, -t1068, 0, t1090, t1083, t925, t913, t926, t867, t1083, t1068, -t1090, t1087, t1084, 0, t895, t896, t880, t830, t965, t911, t929, t964, t928, t997, t820, t822, t811, t803, t965, t911, t929, t964, t928, t997, t813, t815, t808, t800; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1195, 0, 0, -g(3), -t1155, 0, t1104, -t1071, t1092, t1102, t1096, t1130, t1034, t1035, t1024, pkin(7) * t1024, t1104, t1092, t1071, t1130, -t1096, t1102, t933, t932, t934, t881, t1102, -t1071, t1096, t1104, t1092, t1130, t899, t900, t887, t835, t967, t912, t931, t966, t930, t998, t823, t824, t817, t804, t967, t912, t931, t966, t930, t998, t816, t818, t809, t801; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1195, 0, qJDD(1), 0, g(3), 0, -t1156, 0, t1272, -t1151, -t1283, -t1272, -t1282, -qJDD(2), t1028, t1027, 0, pkin(1) * t1024, t1272, -t1283, t1151, -qJDD(2), t1282, -t1272, t968, t1137, t963, t890, -t1272, -t1151, -t1282, t1272, -t1283, -qJDD(2), t918, t915, -t1112, t854, t1002, t946, t993, -t1000, t994, t1032, t841, t843, t821, t807, t1002, t946, t993, -t1000, t994, t1032, t828, t829, t814, t802; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1155, t1156, 0, 0, t1087, t1068, t1084, t1083, -t1090, 0, t1009, t1008, t980, t992, t1087, t1084, -t1068, 0, t1090, t1083, t925, t913, t926, t867, t1083, t1068, -t1090, t1087, t1084, 0, t895, t896, t880, t830, t965, t911, t929, t964, t928, t997, t820, t822, t811, t803, t965, t911, t929, t964, t928, t997, t813, t815, t808, t800; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1142, t1144, t1153, -t1172, t1162, t1172, 0, -t1118, t1100, 0, t1142, t1153, -t1144, t1172, -t1162, -t1172, t1344, t1005, t988, qJ(3) * t999, -t1172, t1144, t1162, t1142, t1153, t1172, t937, t1081, t959, t901, t1078, t1075, t1022, -t1078, -t1018, t1126, t850, t851, t865, t812, t1078, t1075, t1022, -t1078, -t1018, t1126, t836, t842, t859, t806; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1171, t1141, -t1160, t1143, -t1154, -t1171, t1118, 0, t1101, 0, t1171, -t1160, -t1141, -t1171, t1154, t1143, t989, t1004, pkin(2) * t1141, pkin(2) * t999, t1143, t1141, -t1154, t1171, -t1160, -t1171, t1051, t927, t949, t872, -t1003, -t948, -t995, t1001, -t996, -t1033, t858, t860, t827, t810, -t1003, -t948, -t995, t1001, -t996, -t1033, t837, t840, t819, t805; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1272, t1151, t1283, t1272, t1282, qJDD(2), -t1100, -t1101, 0, 0, -t1272, t1283, -t1151, qJDD(2), -t1282, t1272, t1413, -t1137, t1026 + t1396, t1246, t1272, t1151, t1282, -t1272, t1283, qJDD(2), t1196 + t1396, t1388 + t970, t1112, t1389, -t1002, -t946, -t993, t1000, -t994, -t1032, t1332 + t1242, -t950 + t1241, t1214, t1390, -t1002, -t946, -t993, t1000, -t994, -t1032, t1240 + t1403, t1241 - t1406, t1221, t1391; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1142, t1153, -t1144, t1172, -t1162, -t1172, 0, t1029, t999, 0, -t1172, t1144, t1162, t1142, t1153, t1172, t962 + t1339, -t1341, t1197, -t1353, t1078, t1075, t1022, -t1078, -t1018, t1126, t1220, t1238, t1268, t1269, t1078, t1075, t1022, -t1078, -t1018, t1126, t1203, t1208, -t1014 + t1268, t1254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1272, t1283, -t1151, qJDD(2), -t1282, t1272, -t1029, 0, t1026, 0, t1272, t1151, t1282, -t1272, t1283, qJDD(2), t1196, t970 - t1414, pkin(3) * t1283, -pkin(3) * t970, -t1002, -t946, -t993, t1000, -t994, -t1032, -t975 - t1276, -t985 - t1277, t1223, -pkin(3) * t849 + t1345, -t1002, -t946, -t993, t1000, -t994, -t1032, t1240, -t985 - t1260, t1247, -pkin(3) * t834 - t1264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1171, t1160, t1141, t1171, -t1154, -t1143, -t999, -t1026, 0, 0, -t1143, -t1141, t1154, -t1171, t1160, t1171, -pkin(3) * t1141 - t1340, t962 + t1395, pkin(3) * t1150 + t1320 + (t1236 - t1327) * t1193 + t1210, -pkin(3) * t962 + t1354, t1003, t948, t995, -t1001, t996, t1033, t1232, t1259, t1224, t1382 * t848 + t1355, t1003, t948, t995, -t1001, t996, t1033, t1212, t1226, t1225, t1382 * t833 - t1215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1143, -t1141, t1154, -t1171, t1160, t1171, 0, t962, t969, 0, t1003, t948, t995, -t1001, t996, t1033, t1267, t1346, t1252, -pkin(8) * t848, t1003, t948, t995, -t1001, t996, t1033, t1219, t1251, t1253, -pkin(8) * t833 + t1228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1172, -t1144, -t1162, -t1142, -t1153, -t1172, -t962, 0, t970, 0, -t1078, -t1075, -t1022, t1078, t1018, -t1126, t1250, t1278, -t940, -t847, -t1078, -t1075, -t1022, t1078, t1018, -t1126, t1419, t1418, t1014 - t940, t1350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1272, -t1151, -t1282, t1272, -t1283, -qJDD(2), -t969, -t970, 0, 0, t1002, t946, t993, -t1000, t994, t1032, t1276, t1277, t1231, -t1345, t1002, t946, t993, -t1000, t994, t1032, t1258, t1260, t1261, t1264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1064, -t1017, -t1412, -t1110, t1106, t1110, 0, t955, t878, 0, t1064, -t1017, -t1412, -t1110, t1106, t1110, qJ(6) * t1412, t886, t857, qJ(6) * t864; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1111, t1398, t1107, t1063, t1049, t1111, -t955, 0, t879, 0, -t1111, t1398, t1107, t1063, t1049, t1111, t869, t961, t853, t845; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1078, t1075, t1022, -t1078, -t1018, t1126, -t878, -t879, 0, 0, t1078, t1075, t1022, -t1078, -t1018, t1126, t1207, t1216, -t1014, -t863; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1064, -t1017, -t1412, -t1110, t1106, t1110, 0, t898, t864, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1111, t1398, t1107, t1063, t1049, t1111, -t898, 0, t866, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1078, t1075, t1022, -t1078, -t1018, t1126, -t864, -t866, 0, 0;];
m_new_reg  = t1;
