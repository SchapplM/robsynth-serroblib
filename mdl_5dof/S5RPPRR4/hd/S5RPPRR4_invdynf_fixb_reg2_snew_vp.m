% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPPRR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:17:37
% EndTime: 2022-01-23 09:17:40
% DurationCPUTime: 3.61s
% Computational Cost: add. (17260->245), mult. (47949->356), div. (0->0), fcn. (33988->10), ass. (0->202)
t1497 = sin(qJ(1));
t1500 = cos(qJ(1));
t1476 = t1497 * g(1) - t1500 * g(2);
t1501 = qJD(1) ^ 2;
t1507 = -t1501 * qJ(2) + qJDD(2) - t1476;
t1492 = sin(pkin(8));
t1494 = cos(pkin(8));
t1512 = -pkin(2) * t1494 - qJ(3) * t1492;
t1540 = qJD(1) * t1492;
t1522 = -0.2e1 * t1540;
t1551 = (-pkin(1) + t1512) * qJDD(1) + t1507 + qJD(3) * t1522;
t1488 = t1492 ^ 2;
t1490 = t1494 ^ 2;
t1528 = t1488 + t1490;
t1469 = t1528 * t1501;
t1534 = t1494 * qJD(1);
t1482 = -qJD(4) + t1534;
t1479 = -qJD(5) + t1482;
t1550 = qJD(5) - t1479;
t1491 = sin(pkin(9));
t1493 = cos(pkin(9));
t1496 = sin(qJ(4));
t1499 = cos(qJ(4));
t1511 = t1491 * t1499 + t1493 * t1496;
t1508 = t1511 * t1492;
t1448 = qJD(1) * t1508;
t1535 = t1492 * t1493;
t1519 = t1499 * t1535;
t1450 = -t1496 * t1491 * t1540 + qJD(1) * t1519;
t1495 = sin(qJ(5));
t1498 = cos(qJ(5));
t1424 = t1498 * t1448 + t1495 * t1450;
t1549 = t1424 ^ 2;
t1426 = -t1495 * t1448 + t1498 * t1450;
t1548 = t1426 ^ 2;
t1547 = t1448 ^ 2;
t1546 = t1450 ^ 2;
t1545 = t1479 ^ 2;
t1544 = t1482 ^ 2;
t1543 = t1491 ^ 2;
t1542 = 2 * qJD(2);
t1541 = t1494 * g(3);
t1539 = t1426 * t1424;
t1538 = t1448 * t1482;
t1537 = t1450 * t1448;
t1536 = t1488 * t1501;
t1533 = t1494 * t1501;
t1466 = t1512 * qJD(1);
t1532 = t1542 + t1466;
t1531 = qJD(4) + t1482;
t1530 = qJD(5) + t1479;
t1477 = -t1500 * g(1) - t1497 * g(2);
t1460 = -t1501 * pkin(1) + qJDD(1) * qJ(2) + t1477;
t1443 = -t1492 * g(3) + t1494 * t1460 + t1534 * t1542;
t1429 = t1466 * t1534 + t1443;
t1510 = -pkin(3) * t1494 - pkin(6) * t1535;
t1529 = t1551 * t1493;
t1393 = t1510 * qJDD(1) + (-t1429 + (-pkin(3) * t1488 * t1493 + pkin(6) * t1492 * t1494) * t1501) * t1491 + t1529;
t1403 = t1493 * t1429 + t1551 * t1491;
t1461 = t1510 * qJD(1);
t1521 = t1491 * t1536;
t1526 = t1492 * qJDD(1);
t1394 = t1461 * t1534 + (-pkin(3) * t1521 - pkin(6) * t1526) * t1491 + t1403;
t1369 = t1496 * t1393 + t1499 * t1394;
t1527 = qJDD(1) * t1493;
t1525 = t1497 * qJDD(1);
t1524 = t1500 * qJDD(1);
t1486 = t1494 * qJDD(1);
t1523 = t1486 - qJDD(4);
t1520 = t1491 * t1533;
t1481 = t1492 * t1533;
t1518 = t1491 * t1526;
t1517 = -qJDD(5) + t1523;
t1516 = qJDD(3) + t1541;
t1368 = t1499 * t1393 - t1496 * t1394;
t1509 = -qJDD(1) * t1519 + t1496 * t1518;
t1433 = -t1448 * qJD(4) - t1509;
t1506 = qJDD(1) * t1508;
t1505 = t1450 * qJD(4) + t1506;
t1514 = -t1495 * t1433 - t1498 * t1505;
t1513 = t1493 * t1521;
t1423 = -t1523 - t1537;
t1504 = -t1498 * t1433 + t1495 * t1505;
t1404 = -t1543 * pkin(6) * t1536 + (pkin(3) * qJDD(1) * t1491 + t1460 + (t1461 * t1493 + t1532) * qJD(1)) * t1492 + t1516;
t1489 = t1493 ^ 2;
t1473 = -t1500 * t1501 - t1525;
t1472 = -t1497 * t1501 + t1524;
t1468 = t1493 * t1481;
t1467 = t1528 * qJDD(1);
t1465 = t1494 * t1469;
t1464 = (-t1488 * t1489 - t1490) * t1501;
t1463 = t1492 * t1469;
t1462 = (-t1488 * t1543 - t1490) * t1501;
t1459 = (t1489 + t1543) * t1536;
t1458 = -t1486 - t1513;
t1457 = t1486 - t1513;
t1456 = qJDD(1) * pkin(1) - t1507;
t1454 = (t1520 - t1527) * t1492;
t1453 = (t1520 + t1527) * t1492;
t1452 = -t1468 - t1518;
t1451 = -t1468 + t1518;
t1442 = qJD(2) * t1522 - t1492 * t1460 - t1541;
t1441 = -t1482 * pkin(4) - t1450 * pkin(7);
t1438 = -t1544 - t1546;
t1437 = t1493 * t1457 - t1491 * t1464;
t1436 = -t1491 * t1458 + t1493 * t1462;
t1435 = t1491 * t1457 + t1493 * t1464;
t1434 = t1493 * t1458 + t1491 * t1462;
t1432 = t1493 * t1452 - t1491 * t1454;
t1431 = t1491 * t1452 + t1493 * t1454;
t1427 = (t1532 * qJD(1) + t1460) * t1492 + t1516;
t1422 = t1523 - t1537;
t1420 = -t1544 - t1547;
t1418 = t1494 * t1437 + t1492 * t1453;
t1417 = t1494 * t1436 + t1492 * t1451;
t1416 = t1492 * t1437 - t1494 * t1453;
t1415 = t1492 * t1436 - t1494 * t1451;
t1414 = -t1546 - t1547;
t1413 = -t1492 * t1442 + t1494 * t1443;
t1412 = t1494 * t1442 + t1492 * t1443;
t1411 = -t1545 - t1548;
t1410 = t1494 * t1432 - t1492 * t1459;
t1409 = t1492 * t1432 + t1494 * t1459;
t1408 = t1531 * t1448 + t1509;
t1407 = t1433 + t1538;
t1406 = -t1531 * t1450 - t1506;
t1405 = (qJD(4) - t1482) * t1450 + t1511 * t1526;
t1402 = -t1491 * t1429 + t1529;
t1401 = t1499 * t1422 - t1496 * t1438;
t1400 = t1496 * t1422 + t1499 * t1438;
t1399 = -t1517 - t1539;
t1398 = t1517 - t1539;
t1397 = t1499 * t1420 - t1496 * t1423;
t1396 = t1496 * t1420 + t1499 * t1423;
t1395 = -t1545 - t1549;
t1389 = -t1548 - t1549;
t1388 = t1499 * t1406 - t1496 * t1408;
t1387 = t1496 * t1406 + t1499 * t1408;
t1386 = t1498 * t1398 - t1495 * t1411;
t1385 = t1495 * t1398 + t1498 * t1411;
t1384 = -t1491 * t1402 + t1493 * t1403;
t1383 = t1493 * t1402 + t1491 * t1403;
t1382 = -t1491 * t1400 + t1493 * t1401;
t1381 = t1493 * t1400 + t1491 * t1401;
t1380 = t1530 * t1424 + t1504;
t1379 = -t1550 * t1424 - t1504;
t1378 = -t1530 * t1426 + t1514;
t1377 = t1550 * t1426 - t1514;
t1376 = t1505 * pkin(4) - t1547 * pkin(7) + t1450 * t1441 + t1404;
t1375 = t1498 * t1395 - t1495 * t1399;
t1374 = t1495 * t1395 + t1498 * t1399;
t1373 = -t1491 * t1396 + t1493 * t1397;
t1372 = t1493 * t1396 + t1491 * t1397;
t1371 = t1494 * t1384 + t1492 * t1427;
t1370 = t1492 * t1384 - t1494 * t1427;
t1367 = t1494 * t1382 + t1492 * t1407;
t1366 = t1492 * t1382 - t1494 * t1407;
t1365 = t1494 * t1373 + t1492 * t1405;
t1364 = t1492 * t1373 - t1494 * t1405;
t1363 = -t1491 * t1387 + t1493 * t1388;
t1362 = t1493 * t1387 + t1491 * t1388;
t1361 = -t1496 * t1385 + t1499 * t1386;
t1360 = t1499 * t1385 + t1496 * t1386;
t1359 = -t1547 * pkin(4) - t1505 * pkin(7) + t1482 * t1441 + t1369;
t1358 = t1494 * t1363 + t1492 * t1414;
t1357 = t1492 * t1363 - t1494 * t1414;
t1356 = (-t1433 + t1538) * pkin(7) + t1423 * pkin(4) + t1368;
t1355 = t1498 * t1378 - t1495 * t1380;
t1354 = t1495 * t1378 + t1498 * t1380;
t1353 = -t1496 * t1374 + t1499 * t1375;
t1352 = t1499 * t1374 + t1496 * t1375;
t1351 = -t1496 * t1368 + t1499 * t1369;
t1350 = t1499 * t1368 + t1496 * t1369;
t1349 = -t1491 * t1360 + t1493 * t1361;
t1348 = t1493 * t1360 + t1491 * t1361;
t1347 = t1495 * t1356 + t1498 * t1359;
t1346 = t1498 * t1356 - t1495 * t1359;
t1345 = -t1496 * t1354 + t1499 * t1355;
t1344 = t1499 * t1354 + t1496 * t1355;
t1343 = -t1491 * t1352 + t1493 * t1353;
t1342 = t1493 * t1352 + t1491 * t1353;
t1341 = t1494 * t1349 + t1492 * t1379;
t1340 = t1492 * t1349 - t1494 * t1379;
t1339 = -t1491 * t1350 + t1493 * t1351;
t1338 = t1493 * t1350 + t1491 * t1351;
t1337 = t1494 * t1343 + t1492 * t1377;
t1336 = t1492 * t1343 - t1494 * t1377;
t1335 = t1494 * t1339 + t1492 * t1404;
t1334 = t1492 * t1339 - t1494 * t1404;
t1333 = -t1495 * t1346 + t1498 * t1347;
t1332 = t1498 * t1346 + t1495 * t1347;
t1331 = -t1491 * t1344 + t1493 * t1345;
t1330 = t1493 * t1344 + t1491 * t1345;
t1329 = t1494 * t1331 + t1492 * t1389;
t1328 = t1492 * t1331 - t1494 * t1389;
t1327 = -t1496 * t1332 + t1499 * t1333;
t1326 = t1499 * t1332 + t1496 * t1333;
t1325 = -t1491 * t1326 + t1493 * t1327;
t1324 = t1493 * t1326 + t1491 * t1327;
t1323 = t1494 * t1325 + t1492 * t1376;
t1322 = t1492 * t1325 - t1494 * t1376;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1473, -t1472, 0, -t1497 * t1476 + t1500 * t1477, 0, 0, 0, 0, 0, 0, -t1500 * t1465 - t1494 * t1525, t1500 * t1463 + t1492 * t1525, t1500 * t1467 - t1497 * t1469, t1500 * t1413 - t1497 * t1456, 0, 0, 0, 0, 0, 0, t1500 * t1417 + t1497 * t1434, t1500 * t1418 + t1497 * t1435, t1500 * t1410 + t1497 * t1431, t1500 * t1371 + t1497 * t1383, 0, 0, 0, 0, 0, 0, t1500 * t1365 + t1497 * t1372, t1500 * t1367 + t1497 * t1381, t1500 * t1358 + t1497 * t1362, t1335 * t1500 + t1338 * t1497, 0, 0, 0, 0, 0, 0, t1337 * t1500 + t1342 * t1497, t1341 * t1500 + t1348 * t1497, t1329 * t1500 + t1330 * t1497, t1323 * t1500 + t1324 * t1497; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1472, t1473, 0, t1476 * t1500 + t1477 * t1497, 0, 0, 0, 0, 0, 0, -t1465 * t1497 + t1494 * t1524, t1463 * t1497 - t1492 * t1524, t1467 * t1497 + t1469 * t1500, t1413 * t1497 + t1456 * t1500, 0, 0, 0, 0, 0, 0, t1417 * t1497 - t1434 * t1500, t1418 * t1497 - t1435 * t1500, t1410 * t1497 - t1431 * t1500, t1371 * t1497 - t1383 * t1500, 0, 0, 0, 0, 0, 0, t1365 * t1497 - t1372 * t1500, t1367 * t1497 - t1381 * t1500, t1358 * t1497 - t1362 * t1500, t1335 * t1497 - t1338 * t1500, 0, 0, 0, 0, 0, 0, t1337 * t1497 - t1342 * t1500, t1341 * t1497 - t1348 * t1500, t1329 * t1497 - t1330 * t1500, t1323 * t1497 - t1324 * t1500; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1412, 0, 0, 0, 0, 0, 0, t1415, t1416, t1409, t1370, 0, 0, 0, 0, 0, 0, t1364, t1366, t1357, t1334, 0, 0, 0, 0, 0, 0, t1336, t1340, t1328, t1322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1501, -qJDD(1), 0, t1477, 0, 0, 0, 0, 0, 0, -t1465, t1463, t1467, t1413, 0, 0, 0, 0, 0, 0, t1417, t1418, t1410, t1371, 0, 0, 0, 0, 0, 0, t1365, t1367, t1358, t1335, 0, 0, 0, 0, 0, 0, t1337, t1341, t1329, t1323; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1501, 0, t1476, 0, 0, 0, 0, 0, 0, t1486, -t1526, t1469, t1456, 0, 0, 0, 0, 0, 0, -t1434, -t1435, -t1431, -t1383, 0, 0, 0, 0, 0, 0, -t1372, -t1381, -t1362, -t1338, 0, 0, 0, 0, 0, 0, -t1342, -t1348, -t1330, -t1324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1412, 0, 0, 0, 0, 0, 0, t1415, t1416, t1409, t1370, 0, 0, 0, 0, 0, 0, t1364, t1366, t1357, t1334, 0, 0, 0, 0, 0, 0, t1336, t1340, t1328, t1322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1490 * t1501, t1481, t1486, t1443, 0, 0, 0, 0, 0, 0, t1436, t1437, t1432, t1384, 0, 0, 0, 0, 0, 0, t1373, t1382, t1363, t1339, 0, 0, 0, 0, 0, 0, t1343, t1349, t1331, t1325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1481, -t1536, -t1526, t1442, 0, 0, 0, 0, 0, 0, -t1451, -t1453, t1459, -t1427, 0, 0, 0, 0, 0, 0, -t1405, -t1407, -t1414, -t1404, 0, 0, 0, 0, 0, 0, -t1377, -t1379, -t1389, -t1376; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1486, t1526, -t1469, -t1456, 0, 0, 0, 0, 0, 0, t1434, t1435, t1431, t1383, 0, 0, 0, 0, 0, 0, t1372, t1381, t1362, t1338, 0, 0, 0, 0, 0, 0, t1342, t1348, t1330, t1324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1462, t1457, t1452, t1403, 0, 0, 0, 0, 0, 0, t1397, t1401, t1388, t1351, 0, 0, 0, 0, 0, 0, t1353, t1361, t1345, t1327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1458, t1464, t1454, t1402, 0, 0, 0, 0, 0, 0, t1396, t1400, t1387, t1350, 0, 0, 0, 0, 0, 0, t1352, t1360, t1344, t1326; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1451, t1453, -t1459, t1427, 0, 0, 0, 0, 0, 0, t1405, t1407, t1414, t1404, 0, 0, 0, 0, 0, 0, t1377, t1379, t1389, t1376; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1420, t1422, t1406, t1369, 0, 0, 0, 0, 0, 0, t1375, t1386, t1355, t1333; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1423, t1438, t1408, t1368, 0, 0, 0, 0, 0, 0, t1374, t1385, t1354, t1332; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1405, t1407, t1414, t1404, 0, 0, 0, 0, 0, 0, t1377, t1379, t1389, t1376; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1395, t1398, t1378, t1347; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1399, t1411, t1380, t1346; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1377, t1379, t1389, t1376;];
f_new_reg = t1;
