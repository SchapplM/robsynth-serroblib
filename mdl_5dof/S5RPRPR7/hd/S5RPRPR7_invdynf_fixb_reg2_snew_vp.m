% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RPRPR7
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
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPRPR7_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:20:24
% EndTime: 2019-12-31 18:20:26
% DurationCPUTime: 2.47s
% Computational Cost: add. (10840->219), mult. (23625->314), div. (0->0), fcn. (15841->10), ass. (0->178)
t1419 = sin(pkin(9));
t1421 = cos(pkin(9));
t1429 = cos(qJ(3));
t1455 = qJD(1) * t1429;
t1426 = sin(qJ(3));
t1456 = qJD(1) * t1426;
t1378 = t1419 * t1456 - t1421 * t1455;
t1377 = qJD(5) + t1378;
t1463 = qJD(5) + t1377;
t1380 = (t1419 * t1429 + t1421 * t1426) * qJD(1);
t1425 = sin(qJ(5));
t1428 = cos(qJ(5));
t1367 = -t1428 * qJD(3) + t1425 * t1380;
t1462 = t1367 ^ 2;
t1369 = t1425 * qJD(3) + t1428 * t1380;
t1461 = t1369 ^ 2;
t1460 = t1377 ^ 2;
t1459 = t1378 ^ 2;
t1458 = t1380 ^ 2;
t1457 = -2 * qJD(4);
t1454 = qJD(3) * t1378;
t1453 = qJD(3) * t1380;
t1452 = t1369 * t1367;
t1451 = t1380 * t1378;
t1415 = t1429 ^ 2;
t1432 = qJD(1) ^ 2;
t1450 = t1415 * t1432;
t1427 = sin(qJ(1));
t1430 = cos(qJ(1));
t1402 = -t1430 * g(1) - t1427 * g(2);
t1387 = -t1432 * pkin(1) + t1402;
t1420 = sin(pkin(8));
t1422 = cos(pkin(8));
t1401 = t1427 * g(1) - t1430 * g(2);
t1434 = qJDD(1) * pkin(1) + t1401;
t1363 = t1422 * t1387 + t1420 * t1434;
t1357 = -t1432 * pkin(2) + qJDD(1) * pkin(6) + t1363;
t1449 = t1426 * t1357;
t1448 = t1426 * t1432;
t1447 = qJD(5) - t1377;
t1416 = -g(3) + qJDD(2);
t1349 = t1429 * t1357 + t1426 * t1416;
t1414 = t1426 ^ 2;
t1446 = t1414 + t1415;
t1445 = t1426 * qJDD(1);
t1444 = qJD(3) * t1456;
t1443 = qJD(3) * t1455;
t1412 = t1429 * qJDD(1);
t1390 = t1412 - t1444;
t1398 = qJD(3) * pkin(3) - qJ(4) * t1456;
t1332 = -pkin(3) * t1450 + t1390 * qJ(4) - qJD(3) * t1398 + t1349;
t1389 = t1443 + t1445;
t1433 = qJDD(3) * pkin(3) - t1389 * qJ(4) - t1449 + (qJ(4) * qJD(1) * qJD(3) + pkin(3) * t1448 + t1416) * t1429;
t1307 = t1421 * t1332 + t1378 * t1457 + t1419 * t1433;
t1364 = t1421 * t1389 + t1419 * t1390;
t1442 = -t1364 + t1454;
t1441 = t1419 * t1332 - t1421 * t1433;
t1440 = t1428 * qJDD(3) - t1425 * t1364;
t1362 = -t1420 * t1387 + t1422 * t1434;
t1439 = t1419 * t1389 - t1421 * t1390;
t1392 = -t1420 * qJDD(1) - t1422 * t1432;
t1393 = t1422 * qJDD(1) - t1420 * t1432;
t1438 = t1430 * t1392 - t1427 * t1393;
t1437 = -qJDD(5) - t1439;
t1436 = t1427 * t1392 + t1430 * t1393;
t1435 = -t1425 * qJDD(3) - t1428 * t1364;
t1356 = -qJDD(1) * pkin(2) - t1432 * pkin(6) - t1362;
t1344 = t1439 + t1453;
t1336 = -t1390 * pkin(3) - qJ(4) * t1450 + t1398 * t1456 + qJDD(4) + t1356;
t1431 = qJD(3) ^ 2;
t1406 = t1429 * t1448;
t1405 = -t1431 - t1450;
t1404 = -t1414 * t1432 - t1431;
t1400 = -qJDD(3) + t1406;
t1399 = qJDD(3) + t1406;
t1397 = t1446 * t1432;
t1396 = -t1427 * qJDD(1) - t1430 * t1432;
t1395 = t1430 * qJDD(1) - t1427 * t1432;
t1394 = t1446 * qJDD(1);
t1391 = t1412 - 0.2e1 * t1444;
t1388 = 0.2e1 * t1443 + t1445;
t1374 = -t1431 - t1458;
t1373 = t1429 * t1400 - t1426 * t1404;
t1372 = -t1426 * t1399 + t1429 * t1405;
t1371 = t1426 * t1400 + t1429 * t1404;
t1370 = t1429 * t1399 + t1426 * t1405;
t1366 = t1422 * t1394 - t1420 * t1397;
t1365 = t1420 * t1394 + t1422 * t1397;
t1361 = -qJDD(3) - t1451;
t1360 = qJDD(3) - t1451;
t1359 = -t1431 - t1459;
t1358 = t1378 * pkin(4) - t1380 * pkin(7);
t1353 = t1422 * t1373 + t1420 * t1388;
t1352 = t1422 * t1372 - t1420 * t1391;
t1351 = t1420 * t1373 - t1422 * t1388;
t1350 = t1420 * t1372 + t1422 * t1391;
t1348 = t1429 * t1416 - t1449;
t1347 = -t1364 - t1454;
t1345 = -t1439 + t1453;
t1343 = -t1458 - t1459;
t1342 = -t1460 - t1461;
t1341 = t1421 * t1361 - t1419 * t1374;
t1340 = t1419 * t1361 + t1421 * t1374;
t1339 = -t1460 - t1462;
t1338 = -t1420 * t1362 + t1422 * t1363;
t1337 = t1422 * t1362 + t1420 * t1363;
t1335 = -t1461 - t1462;
t1334 = t1421 * t1359 - t1419 * t1360;
t1333 = t1419 * t1359 + t1421 * t1360;
t1331 = t1437 - t1452;
t1330 = -t1437 - t1452;
t1326 = -t1426 * t1348 + t1429 * t1349;
t1325 = t1429 * t1348 + t1426 * t1349;
t1324 = t1421 * t1345 - t1419 * t1347;
t1323 = t1419 * t1345 + t1421 * t1347;
t1322 = t1447 * t1367 + t1435;
t1321 = -t1463 * t1367 - t1435;
t1320 = -t1447 * t1369 + t1440;
t1319 = t1463 * t1369 - t1440;
t1318 = -t1426 * t1340 + t1429 * t1341;
t1317 = t1429 * t1340 + t1426 * t1341;
t1316 = t1422 * t1326 + t1420 * t1356;
t1315 = t1420 * t1326 - t1422 * t1356;
t1314 = t1428 * t1331 - t1425 * t1342;
t1313 = t1425 * t1331 + t1428 * t1342;
t1312 = -t1425 * t1330 + t1428 * t1339;
t1311 = t1428 * t1330 + t1425 * t1339;
t1310 = -t1426 * t1333 + t1429 * t1334;
t1309 = t1429 * t1333 + t1426 * t1334;
t1308 = t1344 * pkin(4) + t1442 * pkin(7) + t1336;
t1306 = t1380 * t1457 - t1441;
t1305 = t1422 * t1318 - t1420 * t1442;
t1304 = t1420 * t1318 + t1422 * t1442;
t1303 = -t1426 * t1323 + t1429 * t1324;
t1302 = t1429 * t1323 + t1426 * t1324;
t1301 = t1422 * t1310 + t1420 * t1344;
t1300 = t1420 * t1310 - t1422 * t1344;
t1299 = t1428 * t1320 - t1425 * t1322;
t1298 = t1425 * t1320 + t1428 * t1322;
t1297 = -t1431 * pkin(4) + qJDD(3) * pkin(7) - t1378 * t1358 + t1307;
t1296 = -qJDD(3) * pkin(4) - t1431 * pkin(7) + ((2 * qJD(4)) + t1358) * t1380 + t1441;
t1295 = t1421 * t1314 + t1419 * t1321;
t1294 = t1419 * t1314 - t1421 * t1321;
t1293 = t1421 * t1312 + t1419 * t1319;
t1292 = t1419 * t1312 - t1421 * t1319;
t1291 = t1422 * t1303 + t1420 * t1343;
t1290 = t1420 * t1303 - t1422 * t1343;
t1289 = t1421 * t1299 + t1419 * t1335;
t1288 = t1419 * t1299 - t1421 * t1335;
t1287 = -t1419 * t1306 + t1421 * t1307;
t1286 = t1421 * t1306 + t1419 * t1307;
t1285 = t1428 * t1297 + t1425 * t1308;
t1284 = -t1425 * t1297 + t1428 * t1308;
t1283 = -t1426 * t1294 + t1429 * t1295;
t1282 = t1429 * t1294 + t1426 * t1295;
t1281 = -t1426 * t1292 + t1429 * t1293;
t1280 = t1429 * t1292 + t1426 * t1293;
t1279 = -t1426 * t1288 + t1429 * t1289;
t1278 = t1429 * t1288 + t1426 * t1289;
t1277 = t1422 * t1283 + t1420 * t1313;
t1276 = t1420 * t1283 - t1422 * t1313;
t1275 = t1422 * t1281 + t1420 * t1311;
t1274 = t1420 * t1281 - t1422 * t1311;
t1273 = -t1426 * t1286 + t1429 * t1287;
t1272 = t1429 * t1286 + t1426 * t1287;
t1271 = t1422 * t1273 + t1420 * t1336;
t1270 = t1420 * t1273 - t1422 * t1336;
t1269 = t1422 * t1279 + t1420 * t1298;
t1268 = t1420 * t1279 - t1422 * t1298;
t1267 = -t1425 * t1284 + t1428 * t1285;
t1266 = t1428 * t1284 + t1425 * t1285;
t1265 = t1421 * t1267 + t1419 * t1296;
t1264 = t1419 * t1267 - t1421 * t1296;
t1263 = -t1426 * t1264 + t1429 * t1265;
t1262 = t1429 * t1264 + t1426 * t1265;
t1261 = t1422 * t1263 + t1420 * t1266;
t1260 = t1420 * t1263 - t1422 * t1266;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1396, -t1395, 0, -t1427 * t1401 + t1430 * t1402, 0, 0, 0, 0, 0, 0, t1438, -t1436, 0, -t1427 * t1337 + t1430 * t1338, 0, 0, 0, 0, 0, 0, -t1427 * t1350 + t1430 * t1352, -t1427 * t1351 + t1430 * t1353, -t1427 * t1365 + t1430 * t1366, -t1427 * t1315 + t1430 * t1316, 0, 0, 0, 0, 0, 0, -t1427 * t1300 + t1430 * t1301, -t1427 * t1304 + t1430 * t1305, -t1427 * t1290 + t1430 * t1291, -t1427 * t1270 + t1430 * t1271, 0, 0, 0, 0, 0, 0, -t1427 * t1274 + t1430 * t1275, -t1427 * t1276 + t1430 * t1277, -t1427 * t1268 + t1430 * t1269, -t1427 * t1260 + t1430 * t1261; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1395, t1396, 0, t1430 * t1401 + t1427 * t1402, 0, 0, 0, 0, 0, 0, t1436, t1438, 0, t1430 * t1337 + t1427 * t1338, 0, 0, 0, 0, 0, 0, t1430 * t1350 + t1427 * t1352, t1430 * t1351 + t1427 * t1353, t1430 * t1365 + t1427 * t1366, t1430 * t1315 + t1427 * t1316, 0, 0, 0, 0, 0, 0, t1430 * t1300 + t1427 * t1301, t1430 * t1304 + t1427 * t1305, t1430 * t1290 + t1427 * t1291, t1430 * t1270 + t1427 * t1271, 0, 0, 0, 0, 0, 0, t1430 * t1274 + t1427 * t1275, t1430 * t1276 + t1427 * t1277, t1430 * t1268 + t1427 * t1269, t1430 * t1260 + t1427 * t1261; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1416, 0, 0, 0, 0, 0, 0, t1370, t1371, 0, t1325, 0, 0, 0, 0, 0, 0, t1309, t1317, t1302, t1272, 0, 0, 0, 0, 0, 0, t1280, t1282, t1278, t1262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1432, -qJDD(1), 0, t1402, 0, 0, 0, 0, 0, 0, t1392, -t1393, 0, t1338, 0, 0, 0, 0, 0, 0, t1352, t1353, t1366, t1316, 0, 0, 0, 0, 0, 0, t1301, t1305, t1291, t1271, 0, 0, 0, 0, 0, 0, t1275, t1277, t1269, t1261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1432, 0, t1401, 0, 0, 0, 0, 0, 0, t1393, t1392, 0, t1337, 0, 0, 0, 0, 0, 0, t1350, t1351, t1365, t1315, 0, 0, 0, 0, 0, 0, t1300, t1304, t1290, t1270, 0, 0, 0, 0, 0, 0, t1274, t1276, t1268, t1260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1416, 0, 0, 0, 0, 0, 0, t1370, t1371, 0, t1325, 0, 0, 0, 0, 0, 0, t1309, t1317, t1302, t1272, 0, 0, 0, 0, 0, 0, t1280, t1282, t1278, t1262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1432, -qJDD(1), 0, t1363, 0, 0, 0, 0, 0, 0, t1372, t1373, t1394, t1326, 0, 0, 0, 0, 0, 0, t1310, t1318, t1303, t1273, 0, 0, 0, 0, 0, 0, t1281, t1283, t1279, t1263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1432, 0, t1362, 0, 0, 0, 0, 0, 0, t1391, -t1388, t1397, -t1356, 0, 0, 0, 0, 0, 0, -t1344, t1442, -t1343, -t1336, 0, 0, 0, 0, 0, 0, -t1311, -t1313, -t1298, -t1266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1416, 0, 0, 0, 0, 0, 0, t1370, t1371, 0, t1325, 0, 0, 0, 0, 0, 0, t1309, t1317, t1302, t1272, 0, 0, 0, 0, 0, 0, t1280, t1282, t1278, t1262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1405, t1400, t1412, t1349, 0, 0, 0, 0, 0, 0, t1334, t1341, t1324, t1287, 0, 0, 0, 0, 0, 0, t1293, t1295, t1289, t1265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1399, t1404, -t1445, t1348, 0, 0, 0, 0, 0, 0, t1333, t1340, t1323, t1286, 0, 0, 0, 0, 0, 0, t1292, t1294, t1288, t1264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1391, t1388, -t1397, t1356, 0, 0, 0, 0, 0, 0, t1344, -t1442, t1343, t1336, 0, 0, 0, 0, 0, 0, t1311, t1313, t1298, t1266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1359, t1361, t1345, t1307, 0, 0, 0, 0, 0, 0, t1312, t1314, t1299, t1267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1360, t1374, t1347, t1306, 0, 0, 0, 0, 0, 0, -t1319, -t1321, -t1335, -t1296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1344, -t1442, t1343, t1336, 0, 0, 0, 0, 0, 0, t1311, t1313, t1298, t1266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1339, t1331, t1320, t1285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1330, t1342, t1322, t1284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1319, t1321, t1335, t1296;];
f_new_reg = t1;