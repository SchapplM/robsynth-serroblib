% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RPPRRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RPPRRR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:51:13
% EndTime: 2019-05-05 15:51:16
% DurationCPUTime: 2.98s
% Computational Cost: add. (10713->248), mult. (21132->262), div. (0->0), fcn. (13401->8), ass. (0->152)
t1330 = sin(qJ(5));
t1331 = sin(qJ(4));
t1334 = cos(qJ(5));
t1335 = cos(qJ(4));
t1292 = (t1330 * t1335 + t1331 * t1334) * qJD(1);
t1291 = qJD(6) + t1292;
t1375 = qJD(6) + t1291;
t1338 = qJD(1) ^ 2;
t1332 = sin(qJ(1));
t1336 = cos(qJ(1));
t1309 = -t1336 * g(1) - t1332 * g(2);
t1342 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t1309;
t1367 = pkin(1) + qJ(3);
t1287 = t1338 * t1367 - qJDD(3) + t1342;
t1365 = qJD(1) * t1335;
t1366 = qJD(1) * t1331;
t1294 = -t1330 * t1366 + t1334 * t1365;
t1324 = qJD(4) + qJD(5);
t1329 = sin(qJ(6));
t1333 = cos(qJ(6));
t1278 = t1294 * t1329 - t1333 * t1324;
t1374 = t1278 ^ 2;
t1280 = t1294 * t1333 + t1324 * t1329;
t1373 = t1280 ^ 2;
t1372 = t1291 ^ 2;
t1371 = t1292 ^ 2;
t1370 = t1294 ^ 2;
t1369 = t1324 ^ 2;
t1368 = 2 * qJD(3);
t1364 = t1278 * t1280;
t1363 = t1292 * t1294;
t1326 = t1331 ^ 2;
t1362 = t1326 * t1338;
t1277 = -qJDD(1) * pkin(7) - t1287;
t1361 = t1335 * t1277;
t1360 = t1335 * t1338;
t1359 = qJD(5) - t1324;
t1358 = qJD(6) - t1291;
t1273 = -g(3) * t1335 + t1331 * t1277;
t1351 = qJD(4) * t1365;
t1356 = t1331 * qJDD(1);
t1298 = -t1351 - t1356;
t1307 = qJD(4) * pkin(4) - pkin(8) * t1365;
t1263 = -pkin(4) * t1362 + pkin(8) * t1298 - qJD(4) * t1307 + t1273;
t1353 = qJD(4) * t1366;
t1355 = t1335 * qJDD(1);
t1299 = -t1353 + t1355;
t1339 = qJDD(4) * pkin(4) - t1299 * pkin(8) + t1361 + (-pkin(8) * qJD(1) * qJD(4) - pkin(4) * t1360 + g(3)) * t1331;
t1241 = t1334 * t1263 + t1330 * t1339;
t1327 = t1335 ^ 2;
t1357 = t1326 + t1327;
t1354 = qJDD(4) + qJDD(5);
t1352 = t1331 * t1360;
t1308 = t1332 * g(1) - t1336 * g(2);
t1240 = -t1330 * t1263 + t1334 * t1339;
t1347 = -t1330 * t1298 - t1334 * t1299;
t1266 = -qJD(5) * t1292 - t1347;
t1350 = t1324 * t1292 - t1266;
t1349 = -t1329 * t1266 + t1333 * t1354;
t1348 = -t1334 * t1298 + t1330 * t1299;
t1345 = t1338 * qJ(2) - qJDD(2) + t1308;
t1344 = -t1333 * t1266 - t1329 * t1354;
t1343 = -qJD(5) * t1294 - qJDD(6) - t1348;
t1257 = (qJD(5) + t1324) * t1294 + t1348;
t1340 = qJDD(1) * t1367 + t1345;
t1286 = qJD(1) * t1368 + t1340;
t1328 = t1338 * pkin(7);
t1265 = -t1298 * pkin(4) - pkin(8) * t1362 - t1328 + (t1307 * t1335 + t1368) * qJD(1) + t1340;
t1337 = qJD(4) ^ 2;
t1311 = -t1327 * t1338 - t1337;
t1310 = -t1337 - t1362;
t1306 = -qJDD(4) - t1352;
t1305 = qJDD(4) - t1352;
t1304 = t1357 * t1338;
t1303 = qJDD(1) * t1332 + t1336 * t1338;
t1302 = qJDD(1) * t1336 - t1332 * t1338;
t1301 = t1357 * qJDD(1);
t1300 = -0.2e1 * t1353 + t1355;
t1297 = 0.2e1 * t1351 + t1356;
t1290 = qJDD(1) * pkin(1) + t1345;
t1289 = pkin(1) * t1338 + t1342;
t1285 = -t1369 - t1370;
t1284 = t1306 * t1335 - t1311 * t1331;
t1283 = -t1305 * t1331 + t1310 * t1335;
t1282 = t1306 * t1331 + t1311 * t1335;
t1281 = t1305 * t1335 + t1310 * t1331;
t1276 = -t1328 + t1286;
t1274 = pkin(5) * t1292 - pkin(9) * t1294;
t1272 = t1331 * g(3) + t1361;
t1271 = -t1354 - t1363;
t1270 = t1354 - t1363;
t1269 = -t1369 - t1371;
t1267 = -t1370 - t1371;
t1264 = -t1372 - t1373;
t1262 = t1271 * t1334 - t1285 * t1330;
t1261 = t1271 * t1330 + t1285 * t1334;
t1260 = t1292 * t1359 + t1347;
t1258 = -t1294 * t1359 - t1348;
t1253 = -t1372 - t1374;
t1252 = -t1272 * t1331 + t1273 * t1335;
t1251 = t1272 * t1335 + t1273 * t1331;
t1250 = t1269 * t1334 - t1270 * t1330;
t1249 = t1269 * t1330 + t1270 * t1334;
t1248 = -t1373 - t1374;
t1247 = t1343 - t1364;
t1246 = -t1343 - t1364;
t1245 = t1278 * t1358 + t1344;
t1244 = -t1278 * t1375 - t1344;
t1243 = -t1280 * t1358 + t1349;
t1242 = t1280 * t1375 - t1349;
t1239 = -t1261 * t1331 + t1262 * t1335;
t1238 = t1261 * t1335 + t1262 * t1331;
t1237 = t1258 * t1334 - t1260 * t1330;
t1236 = t1258 * t1330 + t1260 * t1334;
t1235 = -t1249 * t1331 + t1250 * t1335;
t1234 = t1249 * t1335 + t1250 * t1331;
t1233 = t1247 * t1333 - t1264 * t1329;
t1232 = t1247 * t1329 + t1264 * t1333;
t1231 = -t1246 * t1329 + t1253 * t1333;
t1230 = t1246 * t1333 + t1253 * t1329;
t1229 = pkin(5) * t1257 + pkin(9) * t1350 + t1265;
t1228 = -pkin(5) * t1369 + pkin(9) * t1354 - t1292 * t1274 + t1241;
t1227 = -pkin(5) * t1354 - pkin(9) * t1369 + t1294 * t1274 - t1240;
t1226 = t1243 * t1333 - t1245 * t1329;
t1225 = t1243 * t1329 + t1245 * t1333;
t1224 = -t1240 * t1330 + t1241 * t1334;
t1223 = t1240 * t1334 + t1241 * t1330;
t1222 = -t1236 * t1331 + t1237 * t1335;
t1221 = t1236 * t1335 + t1237 * t1331;
t1220 = t1233 * t1334 + t1244 * t1330;
t1219 = t1233 * t1330 - t1244 * t1334;
t1218 = t1231 * t1334 + t1242 * t1330;
t1217 = t1231 * t1330 - t1242 * t1334;
t1216 = t1226 * t1334 + t1248 * t1330;
t1215 = t1226 * t1330 - t1248 * t1334;
t1214 = t1228 * t1333 + t1229 * t1329;
t1213 = -t1228 * t1329 + t1229 * t1333;
t1212 = -t1223 * t1331 + t1224 * t1335;
t1211 = t1223 * t1335 + t1224 * t1331;
t1210 = -t1219 * t1331 + t1220 * t1335;
t1209 = t1219 * t1335 + t1220 * t1331;
t1208 = -t1217 * t1331 + t1218 * t1335;
t1207 = t1217 * t1335 + t1218 * t1331;
t1206 = -t1215 * t1331 + t1216 * t1335;
t1205 = t1215 * t1335 + t1216 * t1331;
t1204 = -t1213 * t1329 + t1214 * t1333;
t1203 = t1213 * t1333 + t1214 * t1329;
t1202 = t1204 * t1334 + t1227 * t1330;
t1201 = t1204 * t1330 - t1227 * t1334;
t1200 = -t1201 * t1331 + t1202 * t1335;
t1199 = t1201 * t1335 + t1202 * t1331;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1303, -t1302, 0, -t1308 * t1332 + t1309 * t1336, 0, 0, 0, 0, 0, 0, 0, t1303, t1302, -t1289 * t1336 - t1290 * t1332, 0, 0, 0, 0, 0, 0, 0, t1302, -t1303, -t1286 * t1332 - t1287 * t1336, 0, 0, 0, 0, 0, 0, t1281 * t1336 - t1297 * t1332, t1282 * t1336 - t1300 * t1332, -t1301 * t1336 + t1304 * t1332, t1251 * t1336 - t1276 * t1332, 0, 0, 0, 0, 0, 0, t1234 * t1336 - t1257 * t1332, t1238 * t1336 + t1332 * t1350, t1221 * t1336 - t1267 * t1332, t1211 * t1336 - t1265 * t1332, 0, 0, 0, 0, 0, 0, t1207 * t1336 - t1230 * t1332, t1209 * t1336 - t1232 * t1332, t1205 * t1336 - t1225 * t1332, t1199 * t1336 - t1203 * t1332; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1302, -t1303, 0, t1308 * t1336 + t1309 * t1332, 0, 0, 0, 0, 0, 0, 0, -t1302, t1303, -t1289 * t1332 + t1290 * t1336, 0, 0, 0, 0, 0, 0, 0, t1303, t1302, t1286 * t1336 - t1287 * t1332, 0, 0, 0, 0, 0, 0, t1281 * t1332 + t1297 * t1336, t1282 * t1332 + t1300 * t1336, -t1301 * t1332 - t1304 * t1336, t1251 * t1332 + t1276 * t1336, 0, 0, 0, 0, 0, 0, t1234 * t1332 + t1257 * t1336, t1238 * t1332 - t1336 * t1350, t1221 * t1332 + t1267 * t1336, t1211 * t1332 + t1265 * t1336, 0, 0, 0, 0, 0, 0, t1207 * t1332 + t1230 * t1336, t1209 * t1332 + t1232 * t1336, t1205 * t1332 + t1225 * t1336, t1199 * t1332 + t1203 * t1336; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1283, t1284, 0, t1252, 0, 0, 0, 0, 0, 0, t1235, t1239, t1222, t1212, 0, 0, 0, 0, 0, 0, t1208, t1210, t1206, t1200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1338, -qJDD(1), 0, t1309, 0, 0, 0, 0, 0, 0, 0, t1338, qJDD(1), -t1289, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1338, -t1287, 0, 0, 0, 0, 0, 0, t1281, t1282, -t1301, t1251, 0, 0, 0, 0, 0, 0, t1234, t1238, t1221, t1211, 0, 0, 0, 0, 0, 0, t1207, t1209, t1205, t1199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1338, 0, t1308, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t1338, t1290, 0, 0, 0, 0, 0, 0, 0, t1338, qJDD(1), t1286, 0, 0, 0, 0, 0, 0, t1297, t1300, -t1304, t1276, 0, 0, 0, 0, 0, 0, t1257, -t1350, t1267, t1265, 0, 0, 0, 0, 0, 0, t1230, t1232, t1225, t1203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1283, t1284, 0, t1252, 0, 0, 0, 0, 0, 0, t1235, t1239, t1222, t1212, 0, 0, 0, 0, 0, 0, t1208, t1210, t1206, t1200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1283, t1284, 0, t1252, 0, 0, 0, 0, 0, 0, t1235, t1239, t1222, t1212, 0, 0, 0, 0, 0, 0, t1208, t1210, t1206, t1200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1338, -qJDD(1), t1289, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t1338, t1287, 0, 0, 0, 0, 0, 0, -t1281, -t1282, t1301, -t1251, 0, 0, 0, 0, 0, 0, -t1234, -t1238, -t1221, -t1211, 0, 0, 0, 0, 0, 0, -t1207, -t1209, -t1205, -t1199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1338, -t1290, 0, 0, 0, 0, 0, 0, 0, -t1338, -qJDD(1), -t1286, 0, 0, 0, 0, 0, 0, -t1297, -t1300, t1304, -t1276, 0, 0, 0, 0, 0, 0, -t1257, t1350, -t1267, -t1265, 0, 0, 0, 0, 0, 0, -t1230, -t1232, -t1225, -t1203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1283, t1284, 0, t1252, 0, 0, 0, 0, 0, 0, t1235, t1239, t1222, t1212, 0, 0, 0, 0, 0, 0, t1208, t1210, t1206, t1200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1338, -qJDD(1), -t1286, 0, 0, 0, 0, 0, 0, -t1297, -t1300, t1304, -t1276, 0, 0, 0, 0, 0, 0, -t1257, t1350, -t1267, -t1265, 0, 0, 0, 0, 0, 0, -t1230, -t1232, -t1225, -t1203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1338, -t1287, 0, 0, 0, 0, 0, 0, t1281, t1282, -t1301, t1251, 0, 0, 0, 0, 0, 0, t1234, t1238, t1221, t1211, 0, 0, 0, 0, 0, 0, t1207, t1209, t1205, t1199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1310, t1306, -t1356, t1273, 0, 0, 0, 0, 0, 0, t1250, t1262, t1237, t1224, 0, 0, 0, 0, 0, 0, t1218, t1220, t1216, t1202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1305, t1311, -t1355, t1272, 0, 0, 0, 0, 0, 0, t1249, t1261, t1236, t1223, 0, 0, 0, 0, 0, 0, t1217, t1219, t1215, t1201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1297, t1300, -t1304, t1276, 0, 0, 0, 0, 0, 0, t1257, -t1350, t1267, t1265, 0, 0, 0, 0, 0, 0, t1230, t1232, t1225, t1203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1269, t1271, t1258, t1241, 0, 0, 0, 0, 0, 0, t1231, t1233, t1226, t1204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1270, t1285, t1260, t1240, 0, 0, 0, 0, 0, 0, -t1242, -t1244, -t1248, -t1227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1257, -t1350, t1267, t1265, 0, 0, 0, 0, 0, 0, t1230, t1232, t1225, t1203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1253, t1247, t1243, t1214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1246, t1264, t1245, t1213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1242, t1244, t1248, t1227;];
f_new_reg  = t1;
