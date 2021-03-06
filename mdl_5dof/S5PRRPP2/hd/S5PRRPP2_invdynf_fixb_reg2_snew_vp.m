% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5PRRPP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5PRRPP2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynf_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:52
% EndTime: 2019-12-05 16:10:55
% DurationCPUTime: 3.06s
% Computational Cost: add. (5162->184), mult. (11581->224), div. (0->0), fcn. (7833->8), ass. (0->128)
t1513 = sin(pkin(8));
t1515 = cos(pkin(8));
t1520 = cos(qJ(3));
t1542 = qJD(2) * t1520;
t1518 = sin(qJ(3));
t1543 = qJD(2) * t1518;
t1479 = t1513 * t1543 - t1515 * t1542;
t1481 = (t1513 * t1520 + t1515 * t1518) * qJD(2);
t1540 = t1481 * t1479;
t1453 = qJDD(3) + t1540;
t1478 = t1481 ^ 2;
t1522 = qJD(3) ^ 2;
t1551 = -t1478 - t1522;
t1421 = t1513 * t1453 - t1515 * t1551;
t1423 = t1515 * t1453 + t1513 * t1551;
t1409 = t1518 * t1421 - t1520 * t1423;
t1519 = sin(qJ(2));
t1521 = cos(qJ(2));
t1477 = t1479 * qJD(3);
t1529 = qJD(3) * t1542;
t1533 = t1518 * qJDD(2);
t1488 = t1529 + t1533;
t1530 = qJD(3) * t1543;
t1532 = t1520 * qJDD(2);
t1525 = -t1530 + t1532;
t1524 = t1515 * t1488 + t1513 * t1525;
t1550 = -t1524 + t1477;
t1398 = t1521 * t1409 - t1519 * t1550;
t1403 = t1520 * t1421 + t1518 * t1423;
t1514 = sin(pkin(7));
t1516 = cos(pkin(7));
t1576 = t1514 * t1398 + t1516 * t1403;
t1575 = t1516 * t1398 - t1514 * t1403;
t1396 = t1519 * t1409 + t1521 * t1550;
t1454 = qJDD(3) - t1540;
t1458 = t1479 ^ 2;
t1552 = -t1458 - t1522;
t1555 = -t1513 * t1454 + t1515 * t1552;
t1556 = t1515 * t1454 + t1513 * t1552;
t1557 = t1518 * t1555 + t1520 * t1556;
t1527 = -t1513 * t1488 + t1515 * t1525;
t1541 = t1481 * qJD(3);
t1435 = -t1527 + t1541;
t1558 = -t1518 * t1556 + t1520 * t1555;
t1565 = t1519 * t1435 + t1521 * t1558;
t1572 = t1514 * t1565 - t1516 * t1557;
t1571 = t1514 * t1557 + t1516 * t1565;
t1568 = -t1521 * t1435 + t1519 * t1558;
t1440 = t1477 + t1524;
t1526 = t1527 + t1541;
t1548 = t1513 * t1440 + t1515 * t1526;
t1549 = -t1515 * t1440 + t1513 * t1526;
t1553 = t1518 * t1548 + t1520 * t1549;
t1434 = t1478 + t1458;
t1554 = -t1518 * t1549 + t1520 * t1548;
t1559 = -t1519 * t1434 + t1521 * t1554;
t1567 = t1514 * t1559 - t1516 * t1553;
t1566 = t1514 * t1553 + t1516 * t1559;
t1560 = t1521 * t1434 + t1519 * t1554;
t1547 = qJD(2) ^ 2;
t1546 = t1520 ^ 2;
t1545 = -2 * qJD(4);
t1544 = -g(3) + qJDD(1);
t1494 = t1514 * g(1) - t1516 * g(2);
t1537 = t1514 * t1494;
t1495 = -t1516 * g(1) - t1514 * g(2);
t1469 = t1521 * t1495 + t1519 * t1544;
t1463 = -t1547 * pkin(2) + qJDD(2) * pkin(6) + t1469;
t1536 = t1518 * t1463;
t1534 = t1546 * t1547;
t1446 = t1520 * t1463 - t1518 * t1494;
t1510 = t1518 ^ 2;
t1531 = t1510 + t1546;
t1496 = qJD(3) * pkin(3) - qJ(4) * t1543;
t1425 = -pkin(3) * t1534 + t1525 * qJ(4) - qJD(3) * t1496 + t1446;
t1523 = qJDD(3) * pkin(3) - t1488 * qJ(4) - t1536 + (-t1494 + (pkin(3) * t1543 + qJ(4) * qJD(3)) * qJD(2)) * t1520;
t1401 = t1515 * t1425 + t1479 * t1545 + t1513 * t1523;
t1528 = t1513 * t1425 - t1515 * t1523;
t1468 = -t1519 * t1495 + t1521 * t1544;
t1462 = -qJDD(2) * pkin(2) - t1547 * pkin(6) - t1468;
t1430 = -t1525 * pkin(3) - qJ(4) * t1534 + t1496 * t1543 + qJDD(4) + t1462;
t1502 = t1518 * t1547 * t1520;
t1501 = -t1522 - t1534;
t1500 = -t1510 * t1547 - t1522;
t1498 = -qJDD(3) + t1502;
t1497 = qJDD(3) + t1502;
t1493 = t1531 * t1547;
t1492 = t1521 * qJDD(2) - t1519 * t1547;
t1491 = -t1519 * qJDD(2) - t1521 * t1547;
t1490 = t1531 * qJDD(2);
t1489 = -0.2e1 * t1530 + t1532;
t1487 = 0.2e1 * t1529 + t1533;
t1483 = t1516 * t1494;
t1467 = t1520 * t1498 - t1518 * t1500;
t1466 = -t1518 * t1497 + t1520 * t1501;
t1465 = t1518 * t1498 + t1520 * t1500;
t1464 = t1520 * t1497 + t1518 * t1501;
t1461 = t1521 * t1490 - t1519 * t1493;
t1460 = t1519 * t1490 + t1521 * t1493;
t1450 = t1479 * pkin(4) - t1481 * qJ(5);
t1448 = t1521 * t1467 + t1519 * t1487;
t1447 = t1521 * t1466 - t1519 * t1489;
t1445 = t1519 * t1467 - t1521 * t1487;
t1444 = t1519 * t1466 + t1521 * t1489;
t1443 = -t1520 * t1494 - t1536;
t1442 = -t1519 * t1468 + t1521 * t1469;
t1441 = t1521 * t1468 + t1519 * t1469;
t1417 = -t1518 * t1443 + t1520 * t1446;
t1416 = t1520 * t1443 + t1518 * t1446;
t1411 = t1521 * t1417 + t1519 * t1462;
t1410 = t1519 * t1417 - t1521 * t1462;
t1400 = t1481 * t1545 - t1528;
t1399 = -t1527 * pkin(4) + (pkin(4) * qJD(3) - (2 * qJD(5))) * t1481 + t1430 + t1550 * qJ(5);
t1386 = qJDD(5) - t1522 * qJ(5) - qJDD(3) * pkin(4) + ((2 * qJD(4)) + t1450) * t1481 + t1528;
t1385 = -t1522 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t1479 * t1450 + t1401;
t1380 = -t1513 * t1400 + t1515 * t1401;
t1379 = t1515 * t1400 + t1513 * t1401;
t1378 = t1515 * t1385 + t1513 * t1386;
t1377 = t1513 * t1385 - t1515 * t1386;
t1376 = -t1518 * t1379 + t1520 * t1380;
t1375 = t1520 * t1379 + t1518 * t1380;
t1374 = t1521 * t1376 + t1519 * t1430;
t1373 = t1519 * t1376 - t1521 * t1430;
t1372 = -t1518 * t1377 + t1520 * t1378;
t1371 = t1520 * t1377 + t1518 * t1378;
t1370 = t1521 * t1372 + t1519 * t1399;
t1369 = t1519 * t1372 - t1521 * t1399;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1516 * t1495 - t1537, 0, 0, 0, 0, 0, 0, t1516 * t1491, -t1516 * t1492, 0, t1516 * t1442 - t1537, 0, 0, 0, 0, 0, 0, t1516 * t1447 + t1514 * t1464, t1516 * t1448 + t1514 * t1465, t1516 * t1461, t1516 * t1411 + t1514 * t1416, 0, 0, 0, 0, 0, 0, t1571, t1575, t1566, t1516 * t1374 + t1514 * t1375, 0, 0, 0, 0, 0, 0, t1571, t1566, -t1575, t1516 * t1370 + t1514 * t1371; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1514 * t1495 + t1483, 0, 0, 0, 0, 0, 0, t1514 * t1491, -t1514 * t1492, 0, t1514 * t1442 + t1483, 0, 0, 0, 0, 0, 0, t1514 * t1447 - t1516 * t1464, t1514 * t1448 - t1516 * t1465, t1514 * t1461, t1514 * t1411 - t1516 * t1416, 0, 0, 0, 0, 0, 0, t1572, t1576, t1567, t1514 * t1374 - t1516 * t1375, 0, 0, 0, 0, 0, 0, t1572, t1567, -t1576, t1514 * t1370 - t1516 * t1371; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1544, 0, 0, 0, 0, 0, 0, t1492, t1491, 0, t1441, 0, 0, 0, 0, 0, 0, t1444, t1445, t1460, t1410, 0, 0, 0, 0, 0, 0, t1568, t1396, t1560, t1373, 0, 0, 0, 0, 0, 0, t1568, t1560, -t1396, t1369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1495, 0, 0, 0, 0, 0, 0, t1491, -t1492, 0, t1442, 0, 0, 0, 0, 0, 0, t1447, t1448, t1461, t1411, 0, 0, 0, 0, 0, 0, t1565, t1398, t1559, t1374, 0, 0, 0, 0, 0, 0, t1565, t1559, -t1398, t1370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1494, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1494, 0, 0, 0, 0, 0, 0, -t1464, -t1465, 0, -t1416, 0, 0, 0, 0, 0, 0, -t1557, t1403, -t1553, -t1375, 0, 0, 0, 0, 0, 0, -t1557, -t1553, -t1403, -t1371; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1544, 0, 0, 0, 0, 0, 0, t1492, t1491, 0, t1441, 0, 0, 0, 0, 0, 0, t1444, t1445, t1460, t1410, 0, 0, 0, 0, 0, 0, t1568, t1396, t1560, t1373, 0, 0, 0, 0, 0, 0, t1568, t1560, -t1396, t1369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1547, -qJDD(2), 0, t1469, 0, 0, 0, 0, 0, 0, t1466, t1467, t1490, t1417, 0, 0, 0, 0, 0, 0, t1558, t1409, t1554, t1376, 0, 0, 0, 0, 0, 0, t1558, t1554, -t1409, t1372; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1547, 0, t1468, 0, 0, 0, 0, 0, 0, t1489, -t1487, t1493, -t1462, 0, 0, 0, 0, 0, 0, -t1435, t1550, t1434, -t1430, 0, 0, 0, 0, 0, 0, -t1435, t1434, -t1550, -t1399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1494, 0, 0, 0, 0, 0, 0, t1464, t1465, 0, t1416, 0, 0, 0, 0, 0, 0, t1557, -t1403, t1553, t1375, 0, 0, 0, 0, 0, 0, t1557, t1553, t1403, t1371; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1501, t1498, t1532, t1446, 0, 0, 0, 0, 0, 0, t1555, -t1423, t1548, t1380, 0, 0, 0, 0, 0, 0, t1555, t1548, t1423, t1378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1497, t1500, -t1533, t1443, 0, 0, 0, 0, 0, 0, t1556, -t1421, t1549, t1379, 0, 0, 0, 0, 0, 0, t1556, t1549, t1421, t1377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1489, t1487, -t1493, t1462, 0, 0, 0, 0, 0, 0, t1435, -t1550, -t1434, t1430, 0, 0, 0, 0, 0, 0, t1435, -t1434, t1550, t1399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1552, -t1453, t1526, t1401, 0, 0, 0, 0, 0, 0, t1552, t1526, t1453, t1385; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1454, t1551, -t1440, t1400, 0, 0, 0, 0, 0, 0, t1454, -t1440, -t1551, -t1386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1435, -t1550, -t1434, t1430, 0, 0, 0, 0, 0, 0, t1435, -t1434, t1550, t1399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1552, t1526, t1453, t1385; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1435, -t1434, t1550, t1399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1454, t1440, t1551, t1386;];
f_new_reg = t1;
