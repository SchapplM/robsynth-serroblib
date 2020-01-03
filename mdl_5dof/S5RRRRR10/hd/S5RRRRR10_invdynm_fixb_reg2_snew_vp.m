% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRRRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRRRR10_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR10_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_invdynm_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:37:06
% EndTime: 2019-12-31 22:37:42
% DurationCPUTime: 38.37s
% Computational Cost: add. (353373->877), mult. (762129->1328), div. (0->0), fcn. (613089->12), ass. (0->626)
t1565 = sin(qJ(3));
t1570 = cos(qJ(3));
t1562 = cos(pkin(5));
t1718 = qJD(1) * t1562;
t1672 = qJD(2) + t1718;
t1561 = sin(pkin(5));
t1566 = sin(qJ(2));
t1695 = t1561 * t1566;
t1682 = qJD(1) * t1695;
t1523 = t1565 * t1672 + t1570 * t1682;
t1610 = t1565 * t1682 - t1570 * t1672;
t1594 = t1523 * t1610;
t1571 = cos(qJ(2));
t1685 = qJDD(1) * t1561;
t1530 = -qJD(2) * t1682 + t1571 * t1685;
t1618 = -qJDD(3) + t1530;
t1733 = -t1594 - t1618;
t1740 = t1565 * t1733;
t1739 = t1570 * t1733;
t1563 = sin(qJ(5));
t1564 = sin(qJ(4));
t1569 = cos(qJ(4));
t1495 = t1569 * t1523 - t1564 * t1610;
t1684 = t1566 * qJDD(1);
t1717 = qJD(1) * t1571;
t1529 = (qJD(2) * t1717 + t1684) * t1561;
t1667 = qJDD(1) * t1562 + qJDD(2);
t1481 = -qJD(3) * t1610 + t1570 * t1529 + t1565 * t1667;
t1588 = -t1565 * t1529 + t1570 * t1667;
t1580 = -t1523 * qJD(3) + t1588;
t1670 = t1564 * t1481 - t1569 * t1580;
t1404 = -t1495 * qJD(4) - t1670;
t1403 = qJDD(5) - t1404;
t1568 = cos(qJ(5));
t1694 = t1561 * t1571;
t1675 = qJD(1) * t1694;
t1548 = -qJD(3) + t1675;
t1620 = -qJD(4) + t1548;
t1468 = t1495 * t1563 + t1568 * t1620;
t1470 = t1568 * t1495 - t1563 * t1620;
t1412 = t1470 * t1468;
t1731 = t1403 - t1412;
t1738 = t1563 * t1731;
t1493 = t1523 * t1564 + t1569 * t1610;
t1434 = t1495 * t1493;
t1524 = -qJDD(4) + t1618;
t1730 = -t1434 - t1524;
t1737 = t1564 * t1730;
t1736 = t1568 * t1731;
t1735 = t1569 * t1730;
t1405 = -t1493 * qJD(4) + t1569 * t1481 + t1564 * t1580;
t1477 = t1493 * t1620;
t1734 = t1477 + t1405;
t1609 = t1610 ^ 2;
t1617 = t1620 ^ 2;
t1666 = t1672 ^ 2;
t1558 = t1561 ^ 2;
t1573 = qJD(1) ^ 2;
t1696 = t1558 * t1573;
t1658 = t1672 * qJD(1);
t1732 = t1558 * (-t1562 * t1573 + t1658);
t1567 = sin(qJ(1));
t1572 = cos(qJ(1));
t1550 = t1567 * g(1) - t1572 * g(2);
t1722 = pkin(7) * t1561;
t1586 = qJDD(1) * pkin(1) + t1573 * t1722 + t1550;
t1614 = t1566 * t1658;
t1615 = t1571 * t1658;
t1720 = t1562 * g(3);
t1721 = t1529 * pkin(8);
t1576 = -t1530 * pkin(2) - t1721 - t1720 + (pkin(2) * t1614 - pkin(8) * t1615 - t1586) * t1561;
t1551 = g(1) * t1572 + g(2) * t1567;
t1525 = -pkin(1) * t1573 + pkin(7) * t1685 - t1551;
t1726 = pkin(2) * t1571;
t1662 = -pkin(8) * t1566 - t1726;
t1719 = qJD(1) * t1561;
t1528 = t1662 * t1719;
t1583 = t1562 * t1586;
t1579 = -g(3) * t1695 + t1566 * t1583;
t1578 = pkin(8) * t1667 + t1579;
t1577 = -t1666 * pkin(2) + (t1528 * t1719 + t1525) * t1571 + t1578;
t1388 = t1565 * t1576 + t1570 * t1577;
t1504 = -pkin(3) * t1548 - pkin(9) * t1523;
t1328 = -pkin(3) * t1609 + pkin(9) * t1580 + t1548 * t1504 + t1388;
t1387 = t1565 * t1577 - t1570 * t1576;
t1512 = t1610 * t1548;
t1729 = t1512 - t1481;
t1575 = pkin(3) * t1733 + pkin(9) * t1729 - t1387;
t1255 = t1569 * t1328 + t1564 * t1575;
t1428 = pkin(4) * t1493 - pkin(10) * t1495;
t1237 = -pkin(4) * t1617 - t1524 * pkin(10) - t1493 * t1428 + t1255;
t1669 = t1566 * t1525 - t1571 * t1583;
t1441 = -t1667 * pkin(2) - t1666 * pkin(8) + (t1566 * qJD(1) * t1528 + t1571 * g(3)) * t1561 + t1669;
t1367 = -t1580 * pkin(3) - t1609 * pkin(9) + t1523 * t1504 + t1441;
t1608 = t1620 * t1495;
t1252 = -t1734 * pkin(10) + (-t1404 - t1608) * pkin(4) + t1367;
t1182 = t1237 * t1563 - t1568 * t1252;
t1183 = t1237 * t1568 + t1252 * t1563;
t1134 = t1563 * t1182 + t1568 * t1183;
t1452 = t1512 + t1481;
t1537 = t1561 * t1614;
t1511 = t1530 - t1537;
t1693 = t1565 * t1523;
t1465 = (t1570 * t1610 - t1693) * t1548;
t1728 = t1566 * t1465 + t1571 * t1618;
t1466 = t1468 ^ 2;
t1467 = t1470 ^ 2;
t1490 = qJD(5) + t1493;
t1489 = t1490 ^ 2;
t1491 = t1493 ^ 2;
t1492 = t1495 ^ 2;
t1522 = t1523 ^ 2;
t1543 = t1548 ^ 2;
t1727 = pkin(2) * t1566;
t1254 = t1328 * t1564 - t1569 * t1575;
t1194 = -t1254 * t1569 + t1255 * t1564;
t1725 = pkin(3) * t1194;
t1369 = t1495 * t1548 + t1670;
t1372 = -t1477 + t1405;
t1286 = -t1369 * t1564 - t1372 * t1569;
t1724 = pkin(3) * t1286;
t1723 = pkin(4) * t1564;
t1716 = t1194 * t1565;
t1715 = t1194 * t1570;
t1338 = t1403 + t1412;
t1714 = t1338 * t1563;
t1713 = t1338 * t1568;
t1712 = t1367 * t1564;
t1711 = t1367 * t1569;
t1419 = -t1434 + t1524;
t1710 = t1419 * t1564;
t1709 = t1419 * t1569;
t1708 = t1441 * t1565;
t1707 = t1441 * t1570;
t1474 = -t1594 + t1618;
t1706 = t1474 * t1565;
t1705 = t1474 * t1570;
t1704 = t1490 * t1563;
t1703 = t1490 * t1568;
t1515 = t1561 * t1586 + t1720;
t1702 = t1515 * t1566;
t1701 = t1515 * t1571;
t1547 = t1571 * t1566 * t1696;
t1526 = -t1547 + t1667;
t1700 = t1526 * t1566;
t1699 = t1526 * t1571;
t1527 = t1547 + t1667;
t1698 = t1527 * t1566;
t1697 = t1527 * t1571;
t1236 = t1524 * pkin(4) - pkin(10) * t1617 + t1428 * t1495 + t1254;
t1233 = t1563 * t1236;
t1444 = -t1512 * t1570 - t1565 * t1580;
t1692 = t1566 * t1444;
t1446 = t1481 * t1570 + t1548 * t1693;
t1691 = t1566 * t1446;
t1234 = t1568 * t1236;
t1689 = t1570 * t1523;
t1688 = t1571 * t1525;
t1687 = -pkin(4) * t1236 + pkin(10) * t1134;
t1559 = t1566 ^ 2;
t1560 = t1571 ^ 2;
t1686 = t1559 + t1560;
t1683 = -pkin(4) * t1569 - pkin(3);
t1681 = t1564 * t1412;
t1680 = t1569 * t1412;
t1679 = t1566 * t1434;
t1678 = t1571 * t1434;
t1677 = t1559 * t1696;
t1676 = t1560 * t1696;
t1402 = -t1467 - t1489;
t1279 = -t1402 * t1563 - t1713;
t1636 = -t1405 * t1568 + t1524 * t1563;
t1319 = (qJD(5) + t1490) * t1468 + t1636;
t1674 = pkin(4) * t1319 + pkin(10) * t1279 + t1233;
t1395 = -t1489 - t1466;
t1274 = t1395 * t1568 - t1738;
t1671 = -t1405 * t1563 - t1568 * t1524;
t1347 = -qJD(5) * t1470 + t1671;
t1423 = t1490 * t1470;
t1315 = t1347 - t1423;
t1673 = pkin(4) * t1315 + pkin(10) * t1274 - t1234;
t1195 = t1254 * t1564 + t1569 * t1255;
t1303 = t1387 * t1565 + t1570 * t1388;
t1668 = -t1550 * t1567 - t1572 * t1551;
t1316 = (-qJD(5) + t1490) * t1470 + t1671;
t1348 = -qJD(5) * t1468 - t1636;
t1422 = t1490 * t1468;
t1318 = t1348 + t1422;
t1248 = t1316 * t1568 + t1318 * t1563;
t1382 = t1466 + t1467;
t1665 = pkin(4) * t1382 + pkin(10) * t1248 + t1134;
t1125 = t1134 * t1564 - t1236 * t1569;
t1664 = pkin(3) * t1125 + t1687;
t1461 = -t1492 - t1617;
t1385 = t1461 * t1569 + t1710;
t1663 = pkin(3) * t1385 - t1255;
t1545 = qJDD(1) * t1572 - t1567 * t1573;
t1661 = -pkin(6) * t1545 - g(3) * t1567;
t1519 = -t1677 - t1666;
t1496 = -t1519 * t1566 - t1699;
t1660 = pkin(7) * t1496 - t1702;
t1534 = -t1666 - t1676;
t1501 = t1534 * t1571 - t1698;
t1659 = pkin(7) * t1501 + t1701;
t1126 = t1134 * t1569 + t1236 * t1564;
t1094 = -t1125 * t1565 + t1126 * t1570;
t1133 = -t1182 * t1568 + t1183 * t1563;
t1657 = t1094 * t1566 - t1133 * t1571;
t1140 = t1195 * t1570 - t1716;
t1656 = t1140 * t1566 - t1367 * t1571;
t1219 = t1248 * t1564 + t1382 * t1569;
t1220 = t1248 * t1569 - t1382 * t1564;
t1161 = -t1219 * t1565 + t1220 * t1570;
t1246 = t1316 * t1563 - t1318 * t1568;
t1655 = t1161 * t1566 - t1246 * t1571;
t1317 = t1348 - t1422;
t1247 = t1315 * t1568 - t1317 * t1563;
t1411 = t1467 - t1466;
t1225 = t1247 * t1564 - t1411 * t1569;
t1226 = t1247 * t1569 + t1411 * t1564;
t1166 = -t1225 * t1565 + t1226 * t1570;
t1245 = t1315 * t1563 + t1317 * t1568;
t1654 = t1166 * t1566 - t1245 * t1571;
t1229 = t1274 * t1564 + t1315 * t1569;
t1230 = t1274 * t1569 - t1315 * t1564;
t1169 = -t1229 * t1565 + t1230 * t1570;
t1273 = t1395 * t1563 + t1736;
t1653 = t1169 * t1566 - t1273 * t1571;
t1231 = t1279 * t1564 + t1319 * t1569;
t1232 = t1279 * t1569 - t1319 * t1564;
t1171 = -t1231 * t1565 + t1232 * t1570;
t1278 = t1402 * t1568 - t1714;
t1652 = t1171 * t1566 - t1278 * t1571;
t1418 = -t1467 + t1489;
t1294 = -t1418 * t1563 + t1736;
t1240 = t1294 * t1564 - t1318 * t1569;
t1242 = t1294 * t1569 + t1318 * t1564;
t1177 = -t1240 * t1565 + t1242 * t1570;
t1292 = t1418 * t1568 + t1738;
t1651 = t1177 * t1566 - t1292 * t1571;
t1417 = t1466 - t1489;
t1295 = t1417 * t1568 - t1714;
t1314 = -t1347 - t1423;
t1241 = t1295 * t1564 + t1314 * t1569;
t1243 = t1295 * t1569 - t1314 * t1564;
t1178 = -t1241 * t1565 + t1243 * t1570;
t1293 = t1417 * t1563 + t1713;
t1650 = t1178 * t1566 - t1293 * t1571;
t1305 = -t1347 * t1563 + t1468 * t1703;
t1266 = t1305 * t1564 + t1680;
t1268 = t1305 * t1569 - t1681;
t1210 = -t1266 * t1565 + t1268 * t1570;
t1304 = -t1568 * t1347 - t1468 * t1704;
t1649 = t1210 * t1566 + t1304 * t1571;
t1307 = t1348 * t1568 - t1470 * t1704;
t1267 = t1307 * t1564 - t1680;
t1269 = t1307 * t1569 + t1681;
t1211 = -t1267 * t1565 + t1269 * t1570;
t1306 = t1348 * t1563 + t1470 * t1703;
t1648 = t1211 * t1566 - t1306 * t1571;
t1368 = (0.2e1 * qJD(4) - t1548) * t1495 + t1670;
t1285 = -t1368 * t1564 + t1569 * t1734;
t1287 = -t1368 * t1569 - t1564 * t1734;
t1223 = -t1285 * t1565 + t1287 * t1570;
t1433 = t1492 - t1491;
t1647 = t1223 * t1566 - t1433 * t1571;
t1288 = -t1369 * t1569 + t1372 * t1564;
t1224 = -t1286 * t1565 + t1288 * t1570;
t1406 = -t1491 - t1492;
t1646 = t1224 * t1566 - t1406 * t1571;
t1357 = (-t1468 * t1568 + t1470 * t1563) * t1490;
t1296 = t1357 * t1564 - t1403 * t1569;
t1297 = t1357 * t1569 + t1403 * t1564;
t1239 = -t1296 * t1565 + t1297 * t1570;
t1356 = (-t1468 * t1563 - t1470 * t1568) * t1490;
t1645 = t1239 * t1566 - t1356 * t1571;
t1425 = -t1617 - t1491;
t1358 = t1425 * t1564 + t1735;
t1359 = t1425 * t1569 - t1737;
t1277 = -t1358 * t1565 + t1359 * t1570;
t1644 = t1277 * t1566 - t1368 * t1571;
t1386 = -t1461 * t1564 + t1709;
t1301 = -t1385 * t1565 + t1386 * t1570;
t1643 = t1301 * t1566 - t1571 * t1734;
t1642 = t1303 * t1566 - t1441 * t1571;
t1473 = -t1492 + t1617;
t1391 = t1473 * t1569 + t1737;
t1393 = -t1473 * t1564 + t1735;
t1312 = -t1391 * t1565 + t1393 * t1570;
t1641 = t1312 * t1566 - t1372 * t1571;
t1472 = t1491 - t1617;
t1392 = t1472 * t1564 - t1709;
t1394 = t1472 * t1569 + t1710;
t1313 = -t1392 * t1565 + t1394 * t1570;
t1640 = t1313 * t1566 + t1369 * t1571;
t1590 = t1569 * t1608;
t1593 = t1564 * t1477;
t1413 = t1593 + t1590;
t1591 = t1569 * t1477;
t1592 = t1564 * t1608;
t1414 = t1591 - t1592;
t1350 = -t1413 * t1565 + t1414 * t1570;
t1639 = t1350 * t1566 + t1524 * t1571;
t1302 = -t1387 * t1570 + t1388 * t1565;
t1513 = t1548 * t1523;
t1450 = t1513 + t1580;
t1398 = t1450 * t1570 - t1452 * t1565;
t1497 = t1522 - t1609;
t1638 = t1398 * t1566 - t1497 * t1571;
t1451 = (-qJD(3) - t1548) * t1523 + t1588;
t1399 = t1451 * t1570 - t1565 * t1729;
t1471 = t1609 + t1522;
t1637 = t1399 * t1566 + t1471 * t1571;
t1486 = -t1543 - t1609;
t1416 = t1486 * t1570 - t1740;
t1635 = t1416 * t1566 + t1450 * t1571;
t1498 = -t1522 - t1543;
t1427 = -t1498 * t1565 + t1705;
t1634 = t1427 * t1566 - t1452 * t1571;
t1506 = -t1522 + t1543;
t1431 = -t1506 * t1565 + t1739;
t1633 = t1431 * t1566 + t1571 * t1729;
t1505 = t1609 - t1543;
t1432 = t1505 * t1570 + t1706;
t1449 = t1513 - t1580;
t1632 = t1432 * t1566 + t1449 * t1571;
t1487 = g(3) * t1694 + t1669;
t1488 = t1579 + t1688;
t1631 = -t1487 * t1571 + t1488 * t1566;
t1424 = t1487 * t1566 + t1488 * t1571;
t1538 = t1561 * t1615;
t1508 = t1538 + t1529;
t1630 = t1508 * t1571 + t1511 * t1566;
t1509 = -t1538 + t1529;
t1510 = t1530 + t1537;
t1629 = -t1509 * t1571 + t1510 * t1566;
t1628 = t1519 * t1571 - t1700;
t1533 = -t1666 + t1676;
t1627 = t1533 * t1566 + t1699;
t1532 = t1666 - t1677;
t1626 = t1532 * t1571 + t1698;
t1625 = t1534 * t1566 + t1697;
t1624 = t1550 * t1572 - t1551 * t1567;
t1623 = pkin(3) * t1229 + t1673;
t1622 = pkin(3) * t1231 + t1674;
t1621 = t1561 * t1667;
t1619 = pkin(3) * t1358 - t1254;
t1616 = t1561 * t1658;
t1360 = t1569 * t1404 - t1593;
t1361 = -t1564 * t1404 - t1591;
t1282 = -t1360 * t1565 + t1361 * t1570;
t1613 = t1282 * t1566 + t1678;
t1362 = t1564 * t1405 - t1590;
t1363 = t1569 * t1405 + t1592;
t1283 = -t1362 * t1565 + t1363 * t1570;
t1612 = t1283 * t1566 - t1678;
t1611 = pkin(3) * t1219 + t1665;
t1085 = pkin(9) * t1126 + (-pkin(10) * t1564 + t1683) * t1133;
t1091 = -pkin(9) * t1125 + (-pkin(10) * t1569 + t1723) * t1133;
t1093 = t1125 * t1570 + t1126 * t1565;
t1071 = -pkin(8) * t1093 - t1085 * t1565 + t1091 * t1570;
t1080 = -pkin(2) * t1093 - t1664;
t1088 = t1094 * t1571 + t1133 * t1566;
t1607 = pkin(7) * t1088 + t1071 * t1566 + t1080 * t1571;
t1127 = -pkin(10) * t1246 - t1133;
t1098 = pkin(9) * t1220 + t1127 * t1564 + t1246 * t1683;
t1105 = -pkin(9) * t1219 + t1127 * t1569 + t1246 * t1723;
t1160 = t1219 * t1570 + t1220 * t1565;
t1084 = -pkin(8) * t1160 - t1098 * t1565 + t1105 * t1570;
t1095 = -pkin(2) * t1160 - t1611;
t1136 = t1161 * t1571 + t1246 * t1566;
t1606 = pkin(7) * t1136 + t1084 * t1566 + t1095 * t1571;
t1153 = -pkin(4) * t1273 + t1182;
t1193 = -pkin(10) * t1273 + t1233;
t1109 = -pkin(3) * t1273 + pkin(9) * t1230 + t1153 * t1569 + t1193 * t1564;
t1118 = -pkin(9) * t1229 - t1153 * t1564 + t1193 * t1569;
t1168 = t1229 * t1570 + t1230 * t1565;
t1089 = -pkin(8) * t1168 - t1109 * t1565 + t1118 * t1570;
t1116 = -pkin(2) * t1168 - t1623;
t1147 = t1169 * t1571 + t1273 * t1566;
t1605 = pkin(7) * t1147 + t1089 * t1566 + t1116 * t1571;
t1154 = -pkin(4) * t1278 + t1183;
t1196 = -pkin(10) * t1278 + t1234;
t1110 = -pkin(3) * t1278 + pkin(9) * t1232 + t1154 * t1569 + t1196 * t1564;
t1119 = -pkin(9) * t1231 - t1154 * t1564 + t1196 * t1569;
t1170 = t1231 * t1570 + t1232 * t1565;
t1090 = -pkin(8) * t1170 - t1110 * t1565 + t1119 * t1570;
t1117 = -pkin(2) * t1170 - t1622;
t1148 = t1171 * t1571 + t1278 * t1566;
t1604 = pkin(7) * t1148 + t1090 * t1566 + t1117 * t1571;
t1139 = t1195 * t1565 + t1715;
t1179 = -pkin(3) * t1367 + pkin(9) * t1195;
t1104 = -pkin(8) * t1139 - pkin(9) * t1715 - t1179 * t1565;
t1120 = -pkin(2) * t1139 - t1725;
t1135 = t1140 * t1571 + t1367 * t1566;
t1603 = pkin(7) * t1135 + t1104 * t1566 + t1120 * t1571;
t1159 = -pkin(3) * t1406 + pkin(9) * t1288 + t1195;
t1167 = -pkin(9) * t1286 - t1194;
t1222 = t1286 * t1570 + t1288 * t1565;
t1111 = -pkin(8) * t1222 - t1159 * t1565 + t1167 * t1570;
t1189 = -pkin(2) * t1222 - t1724;
t1203 = t1224 * t1571 + t1406 * t1566;
t1602 = pkin(7) * t1203 + t1111 * t1566 + t1189 * t1571;
t1250 = -pkin(3) * t1368 + pkin(9) * t1359 - t1711;
t1276 = t1358 * t1570 + t1359 * t1565;
t1284 = -pkin(9) * t1358 + t1712;
t1173 = -pkin(8) * t1276 - t1250 * t1565 + t1284 * t1570;
t1190 = -pkin(2) * t1276 - t1619;
t1249 = t1277 * t1571 + t1368 * t1566;
t1601 = pkin(7) * t1249 + t1173 * t1566 + t1190 * t1571;
t1256 = -pkin(3) * t1734 + pkin(9) * t1386 + t1712;
t1291 = -pkin(9) * t1385 + t1711;
t1300 = t1385 * t1570 + t1386 * t1565;
t1184 = -pkin(8) * t1300 - t1256 * t1565 + t1291 * t1570;
t1197 = -pkin(2) * t1300 - t1663;
t1258 = t1301 * t1571 + t1566 * t1734;
t1600 = pkin(7) * t1258 + t1184 * t1566 + t1197 * t1571;
t1415 = t1486 * t1565 + t1739;
t1327 = t1565 * (t1528 * t1675 + t1578 + t1688) - t1570 * (-pkin(8) * t1538 - t1515 - t1721) + (t1511 * t1570 - t1565 * t1666 - t1415) * pkin(2);
t1374 = -pkin(8) * t1415 + t1708;
t1378 = t1416 * t1571 - t1450 * t1566;
t1599 = pkin(7) * t1378 + t1327 * t1571 + t1374 * t1566;
t1426 = t1498 * t1570 + t1706;
t1331 = -pkin(2) * t1426 + t1388;
t1379 = -pkin(8) * t1426 + t1707;
t1381 = t1427 * t1571 + t1452 * t1566;
t1598 = pkin(7) * t1381 + t1331 * t1571 + t1379 * t1566;
t1462 = t1509 * t1566 + t1510 * t1571;
t1597 = pkin(7) * t1462 + t1424;
t1397 = t1451 * t1565 + t1570 * t1729;
t1263 = -pkin(8) * t1397 - t1302;
t1351 = t1399 * t1571 - t1471 * t1566;
t1595 = pkin(7) * t1351 + t1263 * t1566 - t1397 * t1726;
t1270 = t1303 * t1571 + t1441 * t1566;
t1587 = pkin(7) * t1270 + t1302 * t1662;
t1585 = t1566 * t1594;
t1584 = t1571 * t1594;
t1582 = t1562 * t1584;
t1581 = t1561 * t1696 + t1562 * t1616;
t1544 = qJDD(1) * t1567 + t1572 * t1573;
t1542 = t1562 * t1667;
t1536 = t1686 * t1696;
t1535 = (t1559 - t1560) * t1696;
t1531 = -pkin(6) * t1544 + g(3) * t1572;
t1514 = t1672 * t1686 * t1719;
t1507 = (t1684 + (0.2e1 * qJD(2) + t1718) * t1717) * t1561;
t1503 = t1571 * t1529 - t1559 * t1616;
t1502 = -t1566 * t1530 - t1560 * t1616;
t1500 = t1533 * t1571 - t1700;
t1499 = -t1532 * t1566 + t1697;
t1485 = (t1562 * t1529 + t1571 * t1581) * t1566;
t1484 = (t1561 * t1529 + t1571 * t1732) * t1566;
t1483 = (t1561 * t1530 - t1566 * t1732) * t1571;
t1482 = (t1562 * t1530 - t1566 * t1581) * t1571;
t1464 = (t1565 * t1610 + t1689) * t1548;
t1463 = -t1508 * t1566 + t1511 * t1571;
t1460 = t1511 * t1561 + t1562 * t1625;
t1459 = -t1510 * t1561 + t1562 * t1627;
t1458 = -t1509 * t1561 + t1562 * t1626;
t1457 = -t1511 * t1562 + t1561 * t1625;
t1456 = t1510 * t1562 + t1561 * t1627;
t1455 = t1509 * t1562 + t1561 * t1626;
t1448 = -t1507 * t1561 + t1562 * t1628;
t1447 = t1507 * t1562 + t1561 * t1628;
t1445 = t1481 * t1565 - t1548 * t1689;
t1443 = -t1512 * t1565 + t1570 * t1580;
t1442 = t1571 * t1465 - t1566 * t1618;
t1440 = -t1535 * t1561 + t1562 * t1630;
t1439 = t1536 * t1561 + t1562 * t1629;
t1438 = t1535 * t1562 + t1561 * t1630;
t1437 = -t1536 * t1562 + t1561 * t1629;
t1430 = t1505 * t1565 - t1705;
t1429 = t1506 * t1570 + t1740;
t1410 = t1571 * t1446 + t1585;
t1409 = t1571 * t1444 - t1585;
t1408 = t1515 * t1561 + t1562 * t1631;
t1407 = -t1515 * t1562 + t1561 * t1631;
t1401 = -t1561 * t1464 + t1562 * t1728;
t1400 = t1562 * t1464 + t1561 * t1728;
t1396 = t1450 * t1565 + t1452 * t1570;
t1390 = t1432 * t1571 - t1449 * t1566;
t1389 = t1431 * t1571 - t1566 * t1729;
t1383 = -t1702 + (-t1457 * t1561 - t1460 * t1562) * pkin(7);
t1377 = -t1701 + (-t1447 * t1561 - t1448 * t1562) * pkin(7);
t1376 = -pkin(1) * t1457 + t1487 * t1561 + t1562 * t1659;
t1375 = pkin(1) * t1460 - t1487 * t1562 + t1561 * t1659;
t1366 = t1398 * t1571 + t1497 * t1566;
t1365 = -pkin(1) * t1447 + t1488 * t1561 + t1562 * t1660;
t1364 = pkin(1) * t1448 - t1488 * t1562 + t1561 * t1660;
t1355 = -t1561 * t1445 + t1562 * t1691 - t1582;
t1354 = -t1561 * t1443 + t1562 * t1692 + t1582;
t1353 = t1562 * t1445 + (-t1584 + t1691) * t1561;
t1352 = t1562 * t1443 + (t1584 + t1692) * t1561;
t1349 = t1413 * t1570 + t1414 * t1565;
t1346 = pkin(1) * t1408 + t1424 * t1722;
t1345 = pkin(7) * t1424 * t1562 - pkin(1) * t1407;
t1343 = t1350 * t1571 - t1524 * t1566;
t1342 = -pkin(1) * t1437 + t1562 * t1597;
t1341 = pkin(1) * t1439 + t1561 * t1597;
t1340 = -pkin(2) * t1452 + pkin(8) * t1427 + t1708;
t1336 = -t1430 * t1561 + t1562 * t1632;
t1335 = -t1429 * t1561 + t1562 * t1633;
t1334 = t1430 * t1562 + t1561 * t1632;
t1333 = t1429 * t1562 + t1561 * t1633;
t1332 = (-t1407 * t1561 - t1408 * t1562) * pkin(7);
t1330 = pkin(2) * t1450 + pkin(8) * t1416 - t1707;
t1329 = (-t1437 * t1561 - t1439 * t1562) * pkin(7) - t1631;
t1325 = -t1426 * t1561 + t1562 * t1634;
t1324 = t1426 * t1562 + t1561 * t1634;
t1323 = -t1415 * t1561 + t1562 * t1635;
t1322 = t1415 * t1562 + t1561 * t1635;
t1311 = t1392 * t1570 + t1394 * t1565;
t1310 = t1391 * t1570 + t1393 * t1565;
t1299 = -t1396 * t1561 + t1562 * t1638;
t1298 = t1396 * t1562 + t1561 * t1638;
t1290 = -t1397 * t1561 + t1562 * t1637;
t1289 = t1397 * t1562 + t1561 * t1637;
t1281 = t1362 * t1570 + t1363 * t1565;
t1280 = t1360 * t1570 + t1361 * t1565;
t1272 = -pkin(2) * t1441 + pkin(8) * t1303;
t1265 = -t1349 * t1561 + t1562 * t1639;
t1264 = t1349 * t1562 + t1561 * t1639;
t1262 = t1283 * t1571 + t1679;
t1261 = t1282 * t1571 - t1679;
t1260 = t1313 * t1571 - t1369 * t1566;
t1259 = t1312 * t1571 + t1372 * t1566;
t1257 = pkin(2) * t1471 + pkin(8) * t1399 + t1303;
t1238 = t1296 * t1570 + t1297 * t1565;
t1228 = -t1302 * t1561 + t1562 * t1642;
t1227 = t1302 * t1562 + t1561 * t1642;
t1221 = t1285 * t1570 + t1287 * t1565;
t1218 = -t1311 * t1561 + t1562 * t1640;
t1217 = -t1310 * t1561 + t1562 * t1641;
t1216 = t1311 * t1562 + t1561 * t1640;
t1215 = t1310 * t1562 + t1561 * t1641;
t1214 = t1223 * t1571 + t1433 * t1566;
t1213 = -t1300 * t1561 + t1562 * t1643;
t1212 = t1300 * t1562 + t1561 * t1643;
t1209 = t1267 * t1570 + t1269 * t1565;
t1208 = t1266 * t1570 + t1268 * t1565;
t1207 = -t1281 * t1561 + t1562 * t1612;
t1206 = -t1280 * t1561 + t1562 * t1613;
t1205 = t1281 * t1562 + t1561 * t1612;
t1204 = t1280 * t1562 + t1561 * t1613;
t1202 = t1239 * t1571 + t1356 * t1566;
t1201 = -t1331 * t1566 + t1379 * t1571 + (-t1324 * t1561 - t1325 * t1562) * pkin(7);
t1200 = -t1276 * t1561 + t1562 * t1644;
t1199 = t1276 * t1562 + t1561 * t1644;
t1198 = -t1327 * t1566 + t1374 * t1571 + (-t1322 * t1561 - t1323 * t1562) * pkin(7);
t1192 = -pkin(1) * t1324 - t1340 * t1561 + t1562 * t1598;
t1191 = pkin(1) * t1325 + t1340 * t1562 + t1561 * t1598;
t1188 = -pkin(1) * t1322 - t1330 * t1561 + t1562 * t1599;
t1187 = pkin(1) * t1323 + t1330 * t1562 + t1561 * t1599;
t1186 = t1211 * t1571 + t1306 * t1566;
t1185 = t1210 * t1571 - t1304 * t1566;
t1176 = t1241 * t1570 + t1243 * t1565;
t1175 = t1240 * t1570 + t1242 * t1565;
t1174 = t1397 * t1727 + t1263 * t1571 + (-t1289 * t1561 - t1290 * t1562) * pkin(7);
t1172 = -pkin(2) * t1734 + pkin(8) * t1301 + t1256 * t1570 + t1291 * t1565;
t1165 = t1225 * t1570 + t1226 * t1565;
t1164 = -pkin(2) * t1368 + pkin(8) * t1277 + t1250 * t1570 + t1284 * t1565;
t1163 = -t1238 * t1561 + t1562 * t1645;
t1162 = t1238 * t1562 + t1561 * t1645;
t1158 = -t1221 * t1561 + t1562 * t1647;
t1157 = t1221 * t1562 + t1561 * t1647;
t1156 = -t1222 * t1561 + t1562 * t1646;
t1155 = t1222 * t1562 + t1561 * t1646;
t1152 = t1178 * t1571 + t1293 * t1566;
t1151 = t1177 * t1571 + t1292 * t1566;
t1150 = -pkin(1) * t1289 - t1257 * t1561 + t1562 * t1595;
t1149 = pkin(1) * t1290 + t1257 * t1562 + t1561 * t1595;
t1146 = -t1209 * t1561 + t1562 * t1648;
t1145 = -t1208 * t1561 + t1562 * t1649;
t1144 = t1209 * t1562 + t1561 * t1648;
t1143 = t1208 * t1562 + t1561 * t1649;
t1142 = t1166 * t1571 + t1245 * t1566;
t1141 = (-pkin(8) * t1571 + t1727) * t1302 + (-t1227 * t1561 - t1228 * t1562) * pkin(7);
t1138 = -pkin(1) * t1227 - t1272 * t1561 + t1562 * t1587;
t1137 = pkin(1) * t1228 + t1272 * t1562 + t1561 * t1587;
t1131 = -t1176 * t1561 + t1562 * t1650;
t1130 = -t1175 * t1561 + t1562 * t1651;
t1129 = t1176 * t1562 + t1561 * t1650;
t1128 = t1175 * t1562 + t1561 * t1651;
t1124 = -t1170 * t1561 + t1562 * t1652;
t1123 = t1170 * t1562 + t1561 * t1652;
t1122 = -t1168 * t1561 + t1562 * t1653;
t1121 = t1168 * t1562 + t1561 * t1653;
t1115 = -t1165 * t1561 + t1562 * t1654;
t1114 = t1165 * t1562 + t1561 * t1654;
t1113 = -t1160 * t1561 + t1562 * t1655;
t1112 = t1160 * t1562 + t1561 * t1655;
t1108 = -pkin(2) * t1406 + pkin(8) * t1224 + t1159 * t1570 + t1167 * t1565;
t1107 = t1184 * t1571 - t1197 * t1566 + (-t1212 * t1561 - t1213 * t1562) * pkin(7);
t1106 = t1173 * t1571 - t1190 * t1566 + (-t1199 * t1561 - t1200 * t1562) * pkin(7);
t1103 = -t1139 * t1561 + t1562 * t1656;
t1102 = t1139 * t1562 + t1561 * t1656;
t1101 = -pkin(2) * t1367 + pkin(8) * t1140 - pkin(9) * t1716 + t1179 * t1570;
t1100 = -pkin(1) * t1212 - t1172 * t1561 + t1562 * t1600;
t1099 = pkin(1) * t1213 + t1172 * t1562 + t1561 * t1600;
t1097 = -pkin(1) * t1199 - t1164 * t1561 + t1562 * t1601;
t1096 = pkin(1) * t1200 + t1164 * t1562 + t1561 * t1601;
t1092 = t1111 * t1571 - t1189 * t1566 + (-t1155 * t1561 - t1156 * t1562) * pkin(7);
t1087 = -pkin(2) * t1278 + pkin(8) * t1171 + t1110 * t1570 + t1119 * t1565;
t1086 = -pkin(2) * t1273 + pkin(8) * t1169 + t1109 * t1570 + t1118 * t1565;
t1083 = -pkin(1) * t1155 - t1108 * t1561 + t1562 * t1602;
t1082 = pkin(1) * t1156 + t1108 * t1562 + t1561 * t1602;
t1081 = -pkin(2) * t1246 + pkin(8) * t1161 + t1098 * t1570 + t1105 * t1565;
t1079 = -t1093 * t1561 + t1562 * t1657;
t1078 = t1093 * t1562 + t1561 * t1657;
t1077 = t1104 * t1571 - t1120 * t1566 + (-t1102 * t1561 - t1103 * t1562) * pkin(7);
t1076 = t1090 * t1571 - t1117 * t1566 + (-t1123 * t1561 - t1124 * t1562) * pkin(7);
t1075 = t1089 * t1571 - t1116 * t1566 + (-t1121 * t1561 - t1122 * t1562) * pkin(7);
t1074 = -pkin(1) * t1102 - t1101 * t1561 + t1562 * t1603;
t1073 = pkin(1) * t1103 + t1101 * t1562 + t1561 * t1603;
t1072 = t1084 * t1571 - t1095 * t1566 + (-t1112 * t1561 - t1113 * t1562) * pkin(7);
t1070 = -pkin(1) * t1123 - t1087 * t1561 + t1562 * t1604;
t1069 = pkin(1) * t1124 + t1087 * t1562 + t1561 * t1604;
t1068 = -pkin(1) * t1121 - t1086 * t1561 + t1562 * t1605;
t1067 = pkin(1) * t1122 + t1086 * t1562 + t1561 * t1605;
t1066 = -pkin(2) * t1133 + pkin(8) * t1094 + t1085 * t1570 + t1091 * t1565;
t1065 = -pkin(1) * t1112 - t1081 * t1561 + t1562 * t1606;
t1064 = pkin(1) * t1113 + t1081 * t1562 + t1561 * t1606;
t1063 = t1071 * t1571 - t1080 * t1566 + (-t1078 * t1561 - t1079 * t1562) * pkin(7);
t1062 = -pkin(1) * t1078 - t1066 * t1561 + t1562 * t1607;
t1061 = pkin(1) * t1079 + t1066 * t1562 + t1561 * t1607;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1545, 0, -t1544, 0, t1661, -t1531, -t1624, -pkin(6) * t1624, -t1485 * t1567 + t1503 * t1572, -t1440 * t1567 + t1463 * t1572, -t1458 * t1567 + t1499 * t1572, -t1482 * t1567 + t1502 * t1572, -t1459 * t1567 + t1500 * t1572, t1572 * t1514 + t1567 * t1621, t1572 * t1383 - t1567 * t1376 - pkin(6) * (t1460 * t1572 + t1501 * t1567), t1572 * t1377 - t1567 * t1365 - pkin(6) * (t1448 * t1572 + t1496 * t1567), t1572 * t1329 - t1567 * t1342 - pkin(6) * (t1439 * t1572 + t1462 * t1567), t1572 * t1332 - t1567 * t1345 - pkin(6) * (t1408 * t1572 + t1424 * t1567), -t1355 * t1567 + t1410 * t1572, -t1299 * t1567 + t1366 * t1572, -t1335 * t1567 + t1389 * t1572, -t1354 * t1567 + t1409 * t1572, -t1336 * t1567 + t1390 * t1572, -t1401 * t1567 + t1442 * t1572, t1572 * t1198 - t1567 * t1188 - pkin(6) * (t1323 * t1572 + t1378 * t1567), t1572 * t1201 - t1567 * t1192 - pkin(6) * (t1325 * t1572 + t1381 * t1567), t1572 * t1174 - t1567 * t1150 - pkin(6) * (t1290 * t1572 + t1351 * t1567), t1572 * t1141 - t1567 * t1138 - pkin(6) * (t1228 * t1572 + t1270 * t1567), -t1207 * t1567 + t1262 * t1572, -t1158 * t1567 + t1214 * t1572, -t1217 * t1567 + t1259 * t1572, -t1206 * t1567 + t1261 * t1572, -t1218 * t1567 + t1260 * t1572, -t1265 * t1567 + t1343 * t1572, t1572 * t1106 - t1567 * t1097 - pkin(6) * (t1200 * t1572 + t1249 * t1567), t1572 * t1107 - t1567 * t1100 - pkin(6) * (t1213 * t1572 + t1258 * t1567), t1572 * t1092 - t1567 * t1083 - pkin(6) * (t1156 * t1572 + t1203 * t1567), t1572 * t1077 - t1567 * t1074 - pkin(6) * (t1103 * t1572 + t1135 * t1567), -t1146 * t1567 + t1186 * t1572, -t1115 * t1567 + t1142 * t1572, -t1130 * t1567 + t1151 * t1572, -t1145 * t1567 + t1185 * t1572, -t1131 * t1567 + t1152 * t1572, -t1163 * t1567 + t1202 * t1572, t1572 * t1075 - t1567 * t1068 - pkin(6) * (t1122 * t1572 + t1147 * t1567), t1572 * t1076 - t1567 * t1070 - pkin(6) * (t1124 * t1572 + t1148 * t1567), t1572 * t1072 - t1567 * t1065 - pkin(6) * (t1113 * t1572 + t1136 * t1567), t1572 * t1063 - t1567 * t1062 - pkin(6) * (t1079 * t1572 + t1088 * t1567); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1544, 0, t1545, 0, t1531, t1661, t1668, pkin(6) * t1668, t1485 * t1572 + t1503 * t1567, t1440 * t1572 + t1463 * t1567, t1458 * t1572 + t1499 * t1567, t1482 * t1572 + t1502 * t1567, t1459 * t1572 + t1500 * t1567, t1567 * t1514 - t1572 * t1621, t1567 * t1383 + t1572 * t1376 + pkin(6) * (-t1460 * t1567 + t1501 * t1572), t1567 * t1377 + t1572 * t1365 + pkin(6) * (-t1448 * t1567 + t1496 * t1572), t1567 * t1329 + t1572 * t1342 + pkin(6) * (-t1439 * t1567 + t1462 * t1572), t1567 * t1332 + t1572 * t1345 + pkin(6) * (-t1408 * t1567 + t1424 * t1572), t1355 * t1572 + t1410 * t1567, t1299 * t1572 + t1366 * t1567, t1335 * t1572 + t1389 * t1567, t1354 * t1572 + t1409 * t1567, t1336 * t1572 + t1390 * t1567, t1401 * t1572 + t1442 * t1567, t1567 * t1198 + t1572 * t1188 + pkin(6) * (-t1323 * t1567 + t1378 * t1572), t1567 * t1201 + t1572 * t1192 + pkin(6) * (-t1325 * t1567 + t1381 * t1572), t1567 * t1174 + t1572 * t1150 + pkin(6) * (-t1290 * t1567 + t1351 * t1572), t1567 * t1141 + t1572 * t1138 + pkin(6) * (-t1228 * t1567 + t1270 * t1572), t1207 * t1572 + t1262 * t1567, t1158 * t1572 + t1214 * t1567, t1217 * t1572 + t1259 * t1567, t1206 * t1572 + t1261 * t1567, t1218 * t1572 + t1260 * t1567, t1265 * t1572 + t1343 * t1567, t1567 * t1106 + t1572 * t1097 + pkin(6) * (-t1200 * t1567 + t1249 * t1572), t1567 * t1107 + t1572 * t1100 + pkin(6) * (-t1213 * t1567 + t1258 * t1572), t1567 * t1092 + t1572 * t1083 + pkin(6) * (-t1156 * t1567 + t1203 * t1572), t1567 * t1077 + t1572 * t1074 + pkin(6) * (-t1103 * t1567 + t1135 * t1572), t1146 * t1572 + t1186 * t1567, t1115 * t1572 + t1142 * t1567, t1130 * t1572 + t1151 * t1567, t1145 * t1572 + t1185 * t1567, t1131 * t1572 + t1152 * t1567, t1163 * t1572 + t1202 * t1567, t1567 * t1075 + t1572 * t1068 + pkin(6) * (-t1122 * t1567 + t1147 * t1572), t1567 * t1076 + t1572 * t1070 + pkin(6) * (-t1124 * t1567 + t1148 * t1572), t1567 * t1072 + t1572 * t1065 + pkin(6) * (-t1113 * t1567 + t1136 * t1572), t1567 * t1063 + t1572 * t1062 + pkin(6) * (-t1079 * t1567 + t1088 * t1572); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1550, t1551, 0, 0, t1484, t1438, t1455, t1483, t1456, t1542, t1375, t1364, t1341, t1346, t1353, t1298, t1333, t1352, t1334, t1400, t1187, t1191, t1149, t1137, t1205, t1157, t1215, t1204, t1216, t1264, t1096, t1099, t1082, t1073, t1144, t1114, t1128, t1143, t1129, t1162, t1067, t1069, t1064, t1061; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1573, 0, 0, -g(3), -t1550, 0, t1503, t1463, t1499, t1502, t1500, t1514, t1383, t1377, t1329, t1332, t1410, t1366, t1389, t1409, t1390, t1442, t1198, t1201, t1174, t1141, t1262, t1214, t1259, t1261, t1260, t1343, t1106, t1107, t1092, t1077, t1186, t1142, t1151, t1185, t1152, t1202, t1075, t1076, t1072, t1063; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1573, 0, qJDD(1), 0, g(3), 0, -t1551, 0, t1485, t1440, t1458, t1482, t1459, -t1621, t1376, t1365, t1342, t1345, t1355, t1299, t1335, t1354, t1336, t1401, t1188, t1192, t1150, t1138, t1207, t1158, t1217, t1206, t1218, t1265, t1097, t1100, t1083, t1074, t1146, t1115, t1130, t1145, t1131, t1163, t1068, t1070, t1065, t1062; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1550, t1551, 0, 0, t1484, t1438, t1455, t1483, t1456, t1542, t1375, t1364, t1341, t1346, t1353, t1298, t1333, t1352, t1334, t1400, t1187, t1191, t1149, t1137, t1205, t1157, t1215, t1204, t1216, t1264, t1096, t1099, t1082, t1073, t1144, t1114, t1128, t1143, t1129, t1162, t1067, t1069, t1064, t1061; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1529, t1511, t1527, -t1538, t1533, t1538, 0, -t1515, t1487, 0, t1446, t1398, t1431, t1444, t1432, t1465, t1374, t1379, t1263, -pkin(8) * t1302, t1283, t1223, t1312, t1282, t1313, t1350, t1173, t1184, t1111, t1104, t1211, t1166, t1177, t1210, t1178, t1239, t1089, t1090, t1084, t1071; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1537, t1508, t1532, t1530, t1526, -t1537, t1515, 0, t1488, 0, -t1594, -t1497, t1729, t1594, t1449, t1618, t1327, t1331, -pkin(2) * t1397, -pkin(2) * t1302, -t1434, -t1433, -t1372, t1434, t1369, t1524, t1190, t1197, t1189, t1120, -t1306, -t1245, -t1292, t1304, -t1293, -t1356, t1116, t1117, t1095, t1080; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1547, t1535, t1509, t1547, t1510, t1667, -t1487, -t1488, 0, 0, t1445, t1396, t1429, t1443, t1430, t1464, t1330, t1340, t1257, t1272, t1281, t1221, t1310, t1280, t1311, t1349, t1164, t1172, t1108, t1101, t1209, t1165, t1175, t1208, t1176, t1238, t1086, t1087, t1081, t1066; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1481, t1450, t1733, -t1512, t1505, t1512, 0, t1441, t1387, 0, t1363, t1287, t1393, t1361, t1394, t1414, t1284, t1291, t1167, -pkin(9) * t1194, t1269, t1226, t1242, t1268, t1243, t1297, t1118, t1119, t1105, t1091; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1513, t1452, t1506, t1580, -t1474, t1513, -t1441, 0, t1388, 0, t1362, t1285, t1391, t1360, t1392, t1413, t1250, t1256, t1159, t1179, t1267, t1225, t1240, t1266, t1241, t1296, t1109, t1110, t1098, t1085; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1594, t1497, -t1729, -t1594, -t1449, -t1618, -t1387, -t1388, 0, 0, t1434, t1433, t1372, -t1434, -t1369, -t1524, t1619, t1663, t1724, t1725, t1306, t1245, t1292, -t1304, t1293, t1356, t1623, t1622, t1611, t1664; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1405, -t1368, t1730, -t1477, t1472, t1477, 0, t1367, t1254, 0, t1307, t1247, t1294, t1305, t1295, t1357, t1193, t1196, t1127, -pkin(10) * t1133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1608, t1734, t1473, t1404, -t1419, t1608, -t1367, 0, t1255, 0, -t1412, -t1411, -t1318, t1412, t1314, -t1403, t1153, t1154, -pkin(4) * t1246, -pkin(4) * t1133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1434, t1433, t1372, -t1434, -t1369, -t1524, -t1254, -t1255, 0, 0, t1306, t1245, t1292, -t1304, t1293, t1356, t1673, t1674, t1665, t1687; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1348, t1315, t1731, t1422, t1417, -t1422, 0, t1236, t1182, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1423, t1317, t1418, t1347, t1338, -t1423, -t1236, 0, t1183, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1412, t1411, t1318, -t1412, -t1314, t1403, -t1182, -t1183, 0, 0;];
m_new_reg = t1;
