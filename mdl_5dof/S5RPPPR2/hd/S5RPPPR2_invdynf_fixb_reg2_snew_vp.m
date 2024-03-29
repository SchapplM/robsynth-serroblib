% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
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
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPPPR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:00:27
% EndTime: 2022-01-23 09:00:30
% DurationCPUTime: 3.11s
% Computational Cost: add. (10423->260), mult. (30201->358), div. (0->0), fcn. (21367->10), ass. (0->204)
t1452 = sin(pkin(9));
t1455 = cos(pkin(9));
t1457 = cos(pkin(7));
t1454 = sin(pkin(7));
t1456 = cos(pkin(8));
t1498 = t1454 * t1456;
t1472 = t1452 * t1498 + t1455 * t1457;
t1513 = t1472 * qJDD(1);
t1515 = -qJDD(5) - t1513;
t1413 = t1472 * qJD(1);
t1462 = qJD(1) ^ 2;
t1449 = t1454 ^ 2;
t1451 = t1457 ^ 2;
t1495 = t1449 + t1451;
t1430 = t1495 * t1462;
t1407 = qJD(5) + t1413;
t1514 = qJD(5) + t1407;
t1415 = (-t1452 * t1457 + t1455 * t1498) * qJD(1);
t1458 = sin(qJ(5));
t1460 = cos(qJ(5));
t1453 = sin(pkin(8));
t1504 = qJD(1) * t1454;
t1488 = t1453 * t1504;
t1397 = t1458 * t1415 - t1460 * t1488;
t1512 = t1397 ^ 2;
t1399 = t1460 * t1415 + t1458 * t1488;
t1511 = t1399 ^ 2;
t1510 = t1407 ^ 2;
t1509 = t1413 ^ 2;
t1508 = t1415 ^ 2;
t1507 = 2 * qJD(2);
t1506 = 2 * qJD(4);
t1505 = t1457 * g(3);
t1503 = qJD(1) * t1457;
t1502 = t1399 * t1397;
t1501 = t1415 * t1413;
t1500 = t1449 * t1462;
t1499 = t1451 * t1462;
t1497 = t1457 * t1462;
t1496 = qJD(5) - t1407;
t1494 = qJDD(1) * t1456;
t1493 = t1454 * qJDD(1);
t1446 = t1457 * qJDD(1);
t1459 = sin(qJ(1));
t1492 = t1459 * qJDD(1);
t1461 = cos(qJ(1));
t1491 = t1461 * qJDD(1);
t1490 = -0.2e1 * t1504;
t1489 = qJD(1) * t1413 * t1453;
t1448 = t1453 ^ 2;
t1487 = t1448 * t1500;
t1486 = t1453 * t1497;
t1442 = t1454 * t1497;
t1438 = -t1461 * g(1) - t1459 * g(2);
t1422 = -t1462 * pkin(1) + qJDD(1) * qJ(2) + t1438;
t1403 = -t1454 * g(3) + t1457 * t1422 + t1503 * t1507;
t1475 = -pkin(2) * t1457 - qJ(3) * t1454;
t1427 = t1475 * qJD(1);
t1384 = t1427 * t1503 + t1403;
t1437 = t1459 * g(1) - t1461 * g(2);
t1470 = -t1462 * qJ(2) + qJDD(2) - t1437;
t1466 = (-pkin(1) + t1475) * qJDD(1) + t1470;
t1481 = qJD(3) * t1490;
t1363 = t1456 * t1384 + (t1466 + t1481) * t1453;
t1417 = (pkin(3) * t1453 - qJ(4) * t1456) * t1504;
t1350 = -pkin(3) * t1499 - qJ(4) * t1446 - t1417 * t1488 + t1363;
t1471 = t1486 + t1494;
t1474 = t1422 + (t1507 + t1427) * qJD(1);
t1483 = qJDD(3) + t1505;
t1465 = (-t1471 * qJ(4) + (qJDD(1) * t1453 - t1456 * t1497) * pkin(3) + t1474) * t1454 + t1483;
t1332 = t1455 * t1350 - t1413 * t1506 + t1452 * t1465;
t1485 = t1455 * t1494;
t1484 = t1453 * t1493;
t1482 = (0.2e1 * qJD(3) + t1417) * t1456;
t1480 = t1452 * t1350 - t1455 * t1465;
t1479 = t1453 * t1384 - t1456 * t1466;
t1441 = t1452 * t1446;
t1408 = t1454 * t1485 - t1441;
t1478 = -t1458 * t1408 + t1460 * t1484;
t1477 = t1453 * t1456 * t1500;
t1476 = t1415 * t1488;
t1473 = -qJ(4) * t1499 + qJDD(4) + t1479;
t1469 = -t1460 * t1408 - t1458 * t1484;
t1450 = t1456 ^ 2;
t1434 = -t1461 * t1462 - t1492;
t1433 = -t1459 * t1462 + t1491;
t1429 = t1456 * t1442;
t1428 = t1495 * qJDD(1);
t1426 = t1457 * t1430;
t1425 = (-t1449 * t1450 - t1451) * t1462;
t1424 = t1454 * t1430;
t1423 = (-t1448 * t1449 - t1451) * t1462;
t1421 = (t1448 + t1450) * t1500;
t1420 = -t1446 - t1477;
t1419 = t1446 - t1477;
t1418 = qJDD(1) * pkin(1) - t1470;
t1412 = (t1486 - t1494) * t1454;
t1411 = t1471 * t1454;
t1410 = -t1429 - t1484;
t1409 = -t1429 + t1484;
t1404 = -t1487 - t1508;
t1402 = qJD(2) * t1490 - t1454 * t1422 - t1505;
t1396 = t1456 * t1419 - t1453 * t1425;
t1395 = -t1453 * t1420 + t1456 * t1423;
t1394 = t1453 * t1419 + t1456 * t1425;
t1393 = t1456 * t1420 + t1453 * t1423;
t1392 = t1413 * pkin(4) - t1415 * pkin(6);
t1391 = t1456 * t1410 - t1453 * t1412;
t1390 = t1453 * t1410 + t1456 * t1412;
t1389 = t1441 + (-t1485 - t1489) * t1454;
t1388 = -t1441 + (t1485 - t1489) * t1454;
t1387 = -t1513 + t1476;
t1386 = t1513 + t1476;
t1385 = -t1487 - t1509;
t1383 = -t1484 - t1501;
t1382 = t1484 - t1501;
t1381 = t1474 * t1454 + t1483;
t1378 = -t1508 - t1509;
t1377 = t1457 * t1396 + t1454 * t1411;
t1376 = t1457 * t1395 + t1454 * t1409;
t1375 = t1454 * t1396 - t1457 * t1411;
t1374 = t1454 * t1395 - t1457 * t1409;
t1373 = -t1454 * t1402 + t1457 * t1403;
t1372 = t1457 * t1402 + t1454 * t1403;
t1371 = t1457 * t1391 - t1454 * t1421;
t1370 = t1454 * t1391 + t1457 * t1421;
t1369 = -t1510 - t1511;
t1368 = t1455 * t1383 - t1452 * t1404;
t1367 = t1452 * t1383 + t1455 * t1404;
t1366 = -t1502 + t1515;
t1365 = -t1502 - t1515;
t1364 = -t1510 - t1512;
t1362 = t1456 * t1481 - t1479;
t1359 = t1455 * t1387 - t1452 * t1389;
t1358 = t1452 * t1387 + t1455 * t1389;
t1357 = -t1511 - t1512;
t1356 = -t1452 * t1382 + t1455 * t1385;
t1355 = t1455 * t1382 + t1452 * t1385;
t1354 = t1496 * t1397 + t1469;
t1353 = -t1514 * t1397 - t1469;
t1352 = -t1496 * t1399 + t1478;
t1351 = t1514 * t1399 - t1478;
t1349 = pkin(3) * t1446 + t1482 * t1504 + t1473;
t1347 = t1456 * t1368 + t1453 * t1388;
t1346 = t1453 * t1368 - t1456 * t1388;
t1345 = t1456 * t1356 + t1453 * t1386;
t1344 = t1453 * t1356 - t1456 * t1386;
t1343 = t1456 * t1359 + t1453 * t1378;
t1342 = t1453 * t1359 - t1456 * t1378;
t1341 = t1460 * t1366 - t1458 * t1369;
t1340 = t1458 * t1366 + t1460 * t1369;
t1339 = t1460 * t1364 - t1458 * t1365;
t1338 = t1458 * t1364 + t1460 * t1365;
t1337 = -t1453 * t1362 + t1456 * t1363;
t1336 = t1456 * t1362 + t1453 * t1363;
t1335 = t1457 * t1347 + t1454 * t1367;
t1334 = t1454 * t1347 - t1457 * t1367;
t1333 = -t1408 * pkin(6) + (t1457 * pkin(3) + t1472 * pkin(4)) * qJDD(1) + (t1482 + (pkin(4) * t1415 + pkin(6) * t1413) * t1453) * t1504 + t1473;
t1331 = -0.2e1 * qJD(4) * t1415 - t1480;
t1330 = t1457 * t1337 + t1454 * t1381;
t1329 = t1454 * t1337 - t1457 * t1381;
t1328 = t1460 * t1352 - t1458 * t1354;
t1327 = t1458 * t1352 + t1460 * t1354;
t1326 = t1457 * t1345 + t1454 * t1355;
t1325 = t1454 * t1345 - t1457 * t1355;
t1324 = t1457 * t1343 + t1454 * t1358;
t1323 = t1454 * t1343 - t1457 * t1358;
t1322 = t1455 * t1341 + t1452 * t1353;
t1321 = t1452 * t1341 - t1455 * t1353;
t1320 = t1455 * t1339 + t1452 * t1351;
t1319 = t1452 * t1339 - t1455 * t1351;
t1318 = -pkin(4) * t1487 + pkin(6) * t1484 - t1413 * t1392 + t1332;
t1317 = -pkin(4) * t1484 - pkin(6) * t1487 + (t1506 + t1392) * t1415 + t1480;
t1316 = t1455 * t1328 + t1452 * t1357;
t1315 = t1452 * t1328 - t1455 * t1357;
t1314 = t1456 * t1322 + t1453 * t1340;
t1313 = t1453 * t1322 - t1456 * t1340;
t1312 = -t1452 * t1331 + t1455 * t1332;
t1311 = t1455 * t1331 + t1452 * t1332;
t1310 = t1456 * t1320 + t1453 * t1338;
t1309 = t1453 * t1320 - t1456 * t1338;
t1308 = t1460 * t1318 + t1458 * t1333;
t1307 = -t1458 * t1318 + t1460 * t1333;
t1306 = t1456 * t1316 + t1453 * t1327;
t1305 = t1453 * t1316 - t1456 * t1327;
t1304 = t1456 * t1312 + t1453 * t1349;
t1303 = t1453 * t1312 - t1456 * t1349;
t1302 = t1457 * t1314 + t1454 * t1321;
t1301 = t1454 * t1314 - t1457 * t1321;
t1300 = t1457 * t1310 + t1454 * t1319;
t1299 = t1454 * t1310 - t1457 * t1319;
t1298 = t1457 * t1306 + t1454 * t1315;
t1297 = t1454 * t1306 - t1457 * t1315;
t1296 = -t1458 * t1307 + t1460 * t1308;
t1295 = t1460 * t1307 + t1458 * t1308;
t1294 = t1457 * t1304 + t1454 * t1311;
t1293 = t1454 * t1304 - t1457 * t1311;
t1292 = t1455 * t1296 + t1452 * t1317;
t1291 = t1452 * t1296 - t1455 * t1317;
t1290 = t1456 * t1292 + t1453 * t1295;
t1289 = t1453 * t1292 - t1456 * t1295;
t1288 = t1457 * t1290 + t1454 * t1291;
t1287 = t1454 * t1290 - t1457 * t1291;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1434, -t1433, 0, -t1459 * t1437 + t1461 * t1438, 0, 0, 0, 0, 0, 0, -t1461 * t1426 - t1457 * t1492, t1461 * t1424 + t1454 * t1492, t1461 * t1428 - t1459 * t1430, t1461 * t1373 - t1459 * t1418, 0, 0, 0, 0, 0, 0, t1461 * t1376 + t1459 * t1393, t1461 * t1377 + t1459 * t1394, t1461 * t1371 + t1459 * t1390, t1461 * t1330 + t1459 * t1336, 0, 0, 0, 0, 0, 0, t1461 * t1326 + t1459 * t1344, t1461 * t1335 + t1459 * t1346, t1461 * t1324 + t1459 * t1342, t1461 * t1294 + t1459 * t1303, 0, 0, 0, 0, 0, 0, t1461 * t1300 + t1459 * t1309, t1461 * t1302 + t1459 * t1313, t1461 * t1298 + t1459 * t1305, t1461 * t1288 + t1459 * t1289; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1433, t1434, 0, t1461 * t1437 + t1459 * t1438, 0, 0, 0, 0, 0, 0, -t1459 * t1426 + t1457 * t1491, t1459 * t1424 - t1454 * t1491, t1459 * t1428 + t1461 * t1430, t1459 * t1373 + t1461 * t1418, 0, 0, 0, 0, 0, 0, t1459 * t1376 - t1461 * t1393, t1459 * t1377 - t1461 * t1394, t1459 * t1371 - t1461 * t1390, t1459 * t1330 - t1461 * t1336, 0, 0, 0, 0, 0, 0, t1459 * t1326 - t1461 * t1344, t1459 * t1335 - t1461 * t1346, t1459 * t1324 - t1461 * t1342, t1459 * t1294 - t1461 * t1303, 0, 0, 0, 0, 0, 0, t1459 * t1300 - t1461 * t1309, t1459 * t1302 - t1461 * t1313, t1459 * t1298 - t1461 * t1305, t1459 * t1288 - t1461 * t1289; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1372, 0, 0, 0, 0, 0, 0, t1374, t1375, t1370, t1329, 0, 0, 0, 0, 0, 0, t1325, t1334, t1323, t1293, 0, 0, 0, 0, 0, 0, t1299, t1301, t1297, t1287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1462, -qJDD(1), 0, t1438, 0, 0, 0, 0, 0, 0, -t1426, t1424, t1428, t1373, 0, 0, 0, 0, 0, 0, t1376, t1377, t1371, t1330, 0, 0, 0, 0, 0, 0, t1326, t1335, t1324, t1294, 0, 0, 0, 0, 0, 0, t1300, t1302, t1298, t1288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1462, 0, t1437, 0, 0, 0, 0, 0, 0, t1446, -t1493, t1430, t1418, 0, 0, 0, 0, 0, 0, -t1393, -t1394, -t1390, -t1336, 0, 0, 0, 0, 0, 0, -t1344, -t1346, -t1342, -t1303, 0, 0, 0, 0, 0, 0, -t1309, -t1313, -t1305, -t1289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1372, 0, 0, 0, 0, 0, 0, t1374, t1375, t1370, t1329, 0, 0, 0, 0, 0, 0, t1325, t1334, t1323, t1293, 0, 0, 0, 0, 0, 0, t1299, t1301, t1297, t1287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1499, t1442, t1446, t1403, 0, 0, 0, 0, 0, 0, t1395, t1396, t1391, t1337, 0, 0, 0, 0, 0, 0, t1345, t1347, t1343, t1304, 0, 0, 0, 0, 0, 0, t1310, t1314, t1306, t1290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1442, -t1500, -t1493, t1402, 0, 0, 0, 0, 0, 0, -t1409, -t1411, t1421, -t1381, 0, 0, 0, 0, 0, 0, -t1355, -t1367, -t1358, -t1311, 0, 0, 0, 0, 0, 0, -t1319, -t1321, -t1315, -t1291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1446, t1493, -t1430, -t1418, 0, 0, 0, 0, 0, 0, t1393, t1394, t1390, t1336, 0, 0, 0, 0, 0, 0, t1344, t1346, t1342, t1303, 0, 0, 0, 0, 0, 0, t1309, t1313, t1305, t1289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1423, t1419, t1410, t1363, 0, 0, 0, 0, 0, 0, t1356, t1368, t1359, t1312, 0, 0, 0, 0, 0, 0, t1320, t1322, t1316, t1292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1420, t1425, t1412, t1362, 0, 0, 0, 0, 0, 0, -t1386, -t1388, -t1378, -t1349, 0, 0, 0, 0, 0, 0, -t1338, -t1340, -t1327, -t1295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1409, t1411, -t1421, t1381, 0, 0, 0, 0, 0, 0, t1355, t1367, t1358, t1311, 0, 0, 0, 0, 0, 0, t1319, t1321, t1315, t1291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1385, t1383, t1387, t1332, 0, 0, 0, 0, 0, 0, t1339, t1341, t1328, t1296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1382, t1404, t1389, t1331, 0, 0, 0, 0, 0, 0, -t1351, -t1353, -t1357, -t1317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1386, t1388, t1378, t1349, 0, 0, 0, 0, 0, 0, t1338, t1340, t1327, t1295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1364, t1366, t1352, t1308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1365, t1369, t1354, t1307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1351, t1353, t1357, t1317;];
f_new_reg = t1;
