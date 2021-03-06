% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RPRPP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPRPP3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:13:26
% EndTime: 2019-12-31 18:13:29
% DurationCPUTime: 2.75s
% Computational Cost: add. (3658->220), mult. (9213->216), div. (0->0), fcn. (6216->6), ass. (0->130)
t1481 = cos(qJ(3));
t1476 = sin(pkin(7));
t1479 = sin(qJ(3));
t1514 = t1476 * t1479;
t1477 = cos(pkin(7));
t1519 = qJD(1) * t1477;
t1448 = qJD(1) * t1514 - t1481 * t1519;
t1512 = t1477 * t1479;
t1494 = t1476 * t1481 + t1512;
t1450 = t1494 * qJD(1);
t1516 = t1450 * t1448;
t1417 = qJDD(3) + t1516;
t1445 = t1450 ^ 2;
t1483 = qJD(3) ^ 2;
t1499 = t1483 + t1445;
t1393 = -t1481 * t1417 + t1479 * t1499;
t1524 = t1479 * t1417 + t1481 * t1499;
t1371 = t1477 * t1393 + t1476 * t1524;
t1480 = sin(qJ(1));
t1544 = t1480 * t1371;
t1482 = cos(qJ(1));
t1543 = t1482 * t1371;
t1444 = t1448 ^ 2;
t1414 = -t1483 - t1444;
t1526 = t1516 - qJDD(3);
t1535 = t1481 * t1414 + t1479 * t1526;
t1536 = t1479 * t1414 - t1481 * t1526;
t1358 = t1476 * t1536 - t1477 * t1535;
t1542 = t1480 * t1358;
t1541 = t1482 * t1358;
t1442 = t1450 * qJD(3);
t1467 = t1477 * qJDD(1);
t1462 = t1481 * t1467;
t1503 = t1476 * qJDD(1);
t1446 = t1479 * t1503 - t1462;
t1421 = t1446 + 0.2e1 * t1442;
t1540 = -t1482 * t1421 - t1542;
t1539 = t1480 * t1421 - t1541;
t1368 = t1476 * t1393 - t1477 * t1524;
t1538 = t1476 * t1535 + t1477 * t1536;
t1441 = t1448 * qJD(3);
t1447 = t1494 * qJDD(1);
t1490 = t1447 - t1441;
t1398 = -t1441 + t1490;
t1537 = t1482 * t1398 - t1544;
t1534 = -t1480 * t1398 - t1543;
t1527 = -t1444 - t1445;
t1532 = t1480 * t1527;
t1530 = t1482 * t1527;
t1484 = qJD(1) ^ 2;
t1471 = t1476 ^ 2;
t1472 = t1477 ^ 2;
t1504 = t1471 + t1472;
t1456 = t1504 * t1484;
t1412 = t1448 * pkin(3) - t1450 * qJ(4);
t1525 = -t1483 * pkin(3) + qJDD(3) * qJ(4) - t1448 * t1412;
t1522 = 2 * qJD(4);
t1521 = t1477 * g(3);
t1520 = qJD(3) * qJ(5);
t1515 = t1472 * t1484;
t1511 = t1477 * t1484;
t1508 = t1479 * t1447;
t1502 = t1480 * qJDD(1);
t1501 = t1482 * qJDD(1);
t1500 = -0.2e1 * t1441;
t1498 = t1477 * pkin(2) + pkin(1);
t1459 = t1480 * g(1) - t1482 * g(2);
t1460 = -t1482 * g(1) - t1480 * g(2);
t1452 = -t1484 * pkin(1) + qJDD(1) * qJ(2) + t1460;
t1497 = -0.2e1 * qJD(1) * qJD(2) - t1452;
t1428 = -t1476 * g(3) + 0.2e1 * qJD(2) * t1519 + t1477 * t1452;
t1403 = -pkin(2) * t1515 + pkin(6) * t1467 + t1428;
t1492 = pkin(2) * t1511 + t1497;
t1489 = (-pkin(6) * qJDD(1) + t1492) * t1476;
t1372 = -t1479 * t1403 + t1481 * (t1489 - t1521);
t1495 = -qJDD(2) + t1459;
t1493 = -g(3) * t1512 + t1481 * t1403;
t1491 = t1446 + t1442;
t1365 = -qJDD(3) * pkin(3) - t1483 * qJ(4) + t1450 * t1412 + qJDD(4) - t1372;
t1488 = (t1504 * pkin(6) + qJ(2)) * t1484 + t1495;
t1373 = t1479 * t1489 + t1493;
t1487 = (t1494 * qJ(4) + t1498) * qJDD(1) + t1450 * t1522 + t1488 + (-t1514 * qJDD(1) - 0.2e1 * t1442 + t1462) * pkin(3);
t1461 = t1476 * t1511;
t1458 = -t1482 * t1484 - t1502;
t1457 = -t1480 * t1484 + t1501;
t1455 = t1504 * qJDD(1);
t1454 = t1477 * t1456;
t1453 = t1476 * t1456;
t1443 = qJDD(1) * pkin(1) + t1484 * qJ(2) + t1495;
t1431 = t1481 * t1447;
t1430 = t1481 * t1446;
t1429 = t1479 * t1446;
t1427 = t1497 * t1476 - t1521;
t1422 = t1447 + t1500;
t1416 = t1498 * qJDD(1) + t1488;
t1405 = -t1430 + t1508;
t1404 = -t1429 - t1431;
t1397 = t1441 + t1490;
t1396 = t1442 - t1491;
t1395 = t1442 + t1491;
t1387 = -t1476 * t1427 + t1477 * t1428;
t1386 = t1477 * t1427 + t1476 * t1428;
t1385 = t1479 * t1397 - t1430;
t1384 = t1481 * t1396 + t1508;
t1383 = -t1481 * t1397 - t1429;
t1382 = t1479 * t1396 - t1431;
t1375 = -t1476 * t1404 + t1477 * t1405;
t1374 = t1477 * t1404 + t1476 * t1405;
t1364 = qJD(3) * t1522 + t1373 + t1525;
t1363 = -t1476 * t1383 + t1477 * t1385;
t1362 = -t1476 * t1382 + t1477 * t1384;
t1361 = t1477 * t1383 + t1476 * t1385;
t1360 = t1477 * t1382 + t1476 * t1384;
t1353 = qJ(4) * t1500 + t1487;
t1352 = -t1479 * t1372 + t1481 * t1373;
t1351 = t1481 * t1372 + t1479 * t1373;
t1350 = t1462 * pkin(4) - t1444 * qJ(5) + qJDD(5) + ((-pkin(4) - pkin(6)) * qJDD(1) + t1492) * t1514 + (t1522 - t1520) * qJD(3) + t1493 + t1525;
t1349 = t1490 * pkin(4) + (pkin(4) * t1448 - (2 * qJD(5))) * qJD(3) + t1365 + t1526 * qJ(5);
t1348 = t1450 * (t1450 * pkin(4) - t1520) + t1444 * pkin(4) - t1491 * qJ(5) + 0.2e1 * (-qJD(3) * qJ(4) - qJD(5)) * t1448 + t1487;
t1347 = t1481 * t1364 + t1479 * t1365;
t1346 = t1479 * t1364 - t1481 * t1365;
t1345 = -t1476 * t1351 + t1477 * t1352;
t1344 = t1477 * t1351 + t1476 * t1352;
t1343 = t1479 * t1349 + t1481 * t1350;
t1342 = -t1481 * t1349 + t1479 * t1350;
t1341 = -t1476 * t1346 + t1477 * t1347;
t1340 = t1477 * t1346 + t1476 * t1347;
t1339 = -t1476 * t1342 + t1477 * t1343;
t1338 = t1477 * t1342 + t1476 * t1343;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1458, -t1457, 0, -t1480 * t1459 + t1482 * t1460, 0, 0, 0, 0, 0, 0, -t1482 * t1454 - t1477 * t1502, t1482 * t1453 + t1476 * t1502, t1482 * t1455 - t1480 * t1456, t1482 * t1387 - t1480 * t1443, 0, 0, 0, 0, 0, 0, t1539, t1480 * t1422 + t1543, t1482 * t1362 + t1532, t1482 * t1345 - t1480 * t1416, 0, 0, 0, 0, 0, 0, t1482 * t1363 + t1532, -t1480 * t1395 + t1541, t1534, t1482 * t1341 - t1480 * t1353, 0, 0, 0, 0, 0, 0, t1482 * t1375 + t1532, t1534, t1539, t1482 * t1339 - t1480 * t1348; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1457, t1458, 0, t1482 * t1459 + t1480 * t1460, 0, 0, 0, 0, 0, 0, -t1480 * t1454 + t1477 * t1501, t1480 * t1453 - t1476 * t1501, t1480 * t1455 + t1482 * t1456, t1480 * t1387 + t1482 * t1443, 0, 0, 0, 0, 0, 0, t1540, -t1482 * t1422 + t1544, t1480 * t1362 - t1530, t1480 * t1345 + t1482 * t1416, 0, 0, 0, 0, 0, 0, t1480 * t1363 - t1530, t1482 * t1395 + t1542, t1537, t1480 * t1341 + t1482 * t1353, 0, 0, 0, 0, 0, 0, t1480 * t1375 - t1530, t1537, t1540, t1480 * t1339 + t1482 * t1348; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1386, 0, 0, 0, 0, 0, 0, t1538, t1368, t1360, t1344, 0, 0, 0, 0, 0, 0, t1361, -t1538, -t1368, t1340, 0, 0, 0, 0, 0, 0, t1374, -t1368, t1538, t1338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1484, -qJDD(1), 0, t1460, 0, 0, 0, 0, 0, 0, -t1454, t1453, t1455, t1387, 0, 0, 0, 0, 0, 0, -t1358, t1371, t1362, t1345, 0, 0, 0, 0, 0, 0, t1363, t1358, -t1371, t1341, 0, 0, 0, 0, 0, 0, t1375, -t1371, -t1358, t1339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1484, 0, t1459, 0, 0, 0, 0, 0, 0, t1467, -t1503, t1456, t1443, 0, 0, 0, 0, 0, 0, -t1421, -t1422, -t1527, t1416, 0, 0, 0, 0, 0, 0, -t1527, t1395, t1398, t1353, 0, 0, 0, 0, 0, 0, -t1527, t1398, -t1421, t1348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1386, 0, 0, 0, 0, 0, 0, t1538, t1368, t1360, t1344, 0, 0, 0, 0, 0, 0, t1361, -t1538, -t1368, t1340, 0, 0, 0, 0, 0, 0, t1374, -t1368, t1538, t1338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1515, t1461, t1467, t1428, 0, 0, 0, 0, 0, 0, t1535, t1393, t1384, t1352, 0, 0, 0, 0, 0, 0, t1385, -t1535, -t1393, t1347, 0, 0, 0, 0, 0, 0, t1405, -t1393, t1535, t1343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1461, -t1471 * t1484, -t1503, t1427, 0, 0, 0, 0, 0, 0, t1536, -t1524, t1382, t1351, 0, 0, 0, 0, 0, 0, t1383, -t1536, t1524, t1346, 0, 0, 0, 0, 0, 0, t1404, t1524, t1536, t1342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1467, t1503, -t1456, -t1443, 0, 0, 0, 0, 0, 0, t1421, t1422, t1527, -t1416, 0, 0, 0, 0, 0, 0, t1527, -t1395, -t1398, -t1353, 0, 0, 0, 0, 0, 0, t1527, -t1398, t1421, -t1348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1414, -t1417, t1396, t1373, 0, 0, 0, 0, 0, 0, -t1446, -t1414, t1417, t1364, 0, 0, 0, 0, 0, 0, -t1446, t1417, t1414, t1350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1526, -t1499, -t1447, t1372, 0, 0, 0, 0, 0, 0, -t1397, t1526, t1499, -t1365, 0, 0, 0, 0, 0, 0, -t1447, t1499, -t1526, -t1349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1421, t1422, t1527, -t1416, 0, 0, 0, 0, 0, 0, t1527, -t1395, -t1398, -t1353, 0, 0, 0, 0, 0, 0, t1527, -t1398, t1421, -t1348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1527, -t1395, -t1398, -t1353, 0, 0, 0, 0, 0, 0, t1527, -t1398, t1421, -t1348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1446, t1414, -t1417, -t1364, 0, 0, 0, 0, 0, 0, t1446, -t1417, -t1414, -t1350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1397, -t1526, -t1499, t1365, 0, 0, 0, 0, 0, 0, t1447, -t1499, t1526, t1349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1527, -t1398, t1421, -t1348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1447, -t1499, t1526, t1349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1446, t1417, t1414, t1350;];
f_new_reg = t1;
