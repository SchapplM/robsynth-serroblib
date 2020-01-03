% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRPPR10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRPPR10_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR10_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_invdynf_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:45:28
% EndTime: 2019-12-31 19:45:31
% DurationCPUTime: 3.30s
% Computational Cost: add. (8435->235), mult. (18941->271), div. (0->0), fcn. (12372->8), ass. (0->158)
t1637 = sin(qJ(2));
t1671 = qJD(1) * t1637;
t1624 = qJD(2) * t1671;
t1640 = cos(qJ(2));
t1657 = t1640 * qJDD(1);
t1608 = -t1624 + t1657;
t1634 = sin(pkin(8));
t1635 = cos(pkin(8));
t1601 = -t1635 * qJD(2) + t1634 * t1671;
t1603 = t1634 * qJD(2) + t1635 * t1671;
t1669 = t1603 * t1601;
t1567 = t1608 - t1669;
t1599 = t1603 ^ 2;
t1632 = t1640 ^ 2;
t1642 = qJD(1) ^ 2;
t1628 = t1632 * t1642;
t1680 = -t1599 - t1628;
t1546 = t1635 * t1567 - t1634 * t1680;
t1661 = t1640 * qJD(1);
t1655 = qJD(2) * t1661;
t1658 = t1637 * qJDD(1);
t1607 = t1655 + t1658;
t1651 = t1634 * qJDD(2) + t1635 * t1607;
t1656 = t1601 * t1661;
t1645 = t1651 + t1656;
t1530 = t1640 * t1546 + t1637 * t1645;
t1544 = t1634 * t1567 + t1635 * t1680;
t1638 = sin(qJ(1));
t1641 = cos(qJ(1));
t1704 = t1638 * t1530 - t1641 * t1544;
t1703 = t1641 * t1530 + t1638 * t1544;
t1528 = t1637 * t1546 - t1640 * t1645;
t1646 = -t1651 + t1656;
t1591 = t1603 * t1661;
t1653 = -t1635 * qJDD(2) + t1634 * t1607;
t1681 = -t1591 - t1653;
t1687 = t1634 * t1681 + t1635 * t1646;
t1675 = t1601 ^ 2;
t1559 = t1599 + t1675;
t1685 = -t1634 * t1646 + t1635 * t1681;
t1689 = -t1637 * t1559 + t1640 * t1685;
t1700 = t1638 * t1689 - t1641 * t1687;
t1568 = t1608 + t1669;
t1679 = -t1675 - t1628;
t1688 = -t1635 * t1568 + t1634 * t1679;
t1561 = -t1591 + t1653;
t1686 = t1634 * t1568 + t1635 * t1679;
t1691 = t1637 * t1561 + t1640 * t1686;
t1699 = t1638 * t1691 - t1641 * t1688;
t1698 = t1638 * t1687 + t1641 * t1689;
t1697 = t1638 * t1688 + t1641 * t1691;
t1692 = -t1640 * t1561 + t1637 * t1686;
t1690 = t1640 * t1559 + t1637 * t1685;
t1621 = qJD(5) + t1661;
t1682 = qJD(5) + t1621;
t1678 = qJD(2) ^ 2;
t1636 = sin(qJ(5));
t1639 = cos(qJ(5));
t1571 = -t1639 * t1601 + t1636 * t1603;
t1677 = t1571 ^ 2;
t1573 = t1636 * t1601 + t1639 * t1603;
t1676 = t1573 ^ 2;
t1674 = t1621 ^ 2;
t1673 = -0.2e1 * t1603;
t1672 = t1640 * g(3);
t1670 = t1573 * t1571;
t1631 = t1637 ^ 2;
t1668 = t1631 * t1642;
t1617 = -t1641 * g(1) - t1638 * g(2);
t1597 = -t1642 * pkin(1) + qJDD(1) * pkin(6) + t1617;
t1662 = t1637 * t1597;
t1576 = t1601 * pkin(3) - t1603 * qJ(4);
t1660 = (2 * qJD(3)) + t1576;
t1659 = qJD(5) - t1621;
t1584 = -t1637 * g(3) + t1640 * t1597;
t1605 = (-pkin(2) * t1640 - qJ(3) * t1637) * qJD(1);
t1555 = -t1678 * pkin(2) + qJDD(2) * qJ(3) + t1605 * t1661 + t1584;
t1616 = t1638 * g(1) - t1641 * g(2);
t1596 = qJDD(1) * pkin(1) + t1642 * pkin(6) + t1616;
t1643 = (-t1607 - t1655) * qJ(3) + (-t1608 + t1624) * pkin(2) - t1596;
t1527 = -0.2e1 * qJD(3) * t1601 + t1635 * t1555 + t1634 * t1643;
t1654 = t1634 * t1555 - t1635 * t1643;
t1652 = -qJDD(5) - t1608;
t1650 = pkin(4) * t1661 - t1603 * pkin(7);
t1649 = -qJDD(2) * pkin(2) - t1678 * qJ(3) + qJDD(3) + t1672;
t1648 = t1608 * pkin(3) - qJ(4) * t1628 + qJDD(4) + t1654;
t1647 = -t1636 * t1651 + t1639 * t1653;
t1507 = -pkin(3) * t1628 - t1608 * qJ(4) - 0.2e1 * qJD(4) * t1661 - t1601 * t1576 + t1527;
t1644 = -t1636 * t1653 - t1639 * t1651;
t1511 = qJD(4) * t1673 + t1662 - t1651 * qJ(4) + (-t1601 * t1640 * qJ(4) + t1637 * t1605) * qJD(1) + t1649 + t1561 * pkin(3);
t1620 = t1640 * t1642 * t1637;
t1619 = -t1628 - t1678;
t1618 = -t1668 - t1678;
t1615 = -qJDD(2) + t1620;
t1614 = qJDD(2) + t1620;
t1613 = t1628 + t1668;
t1612 = -t1638 * qJDD(1) - t1641 * t1642;
t1611 = t1641 * qJDD(1) - t1638 * t1642;
t1610 = (t1631 + t1632) * qJDD(1);
t1609 = -0.2e1 * t1624 + t1657;
t1606 = 0.2e1 * t1655 + t1658;
t1583 = -t1662 - t1672;
t1582 = t1640 * t1615 - t1637 * t1618;
t1581 = -t1637 * t1614 + t1640 * t1619;
t1580 = t1637 * t1615 + t1640 * t1618;
t1579 = t1640 * t1614 + t1637 * t1619;
t1554 = (qJD(1) * t1605 + t1597) * t1637 + t1649;
t1553 = -t1674 - t1676;
t1551 = -t1637 * t1583 + t1640 * t1584;
t1550 = t1640 * t1583 + t1637 * t1584;
t1543 = -t1674 - t1677;
t1538 = t1652 - t1670;
t1537 = -t1652 - t1670;
t1532 = -t1676 - t1677;
t1526 = qJD(3) * t1673 - t1654;
t1521 = t1659 * t1571 + t1644;
t1520 = -t1682 * t1571 - t1644;
t1519 = -t1659 * t1573 + t1647;
t1518 = t1682 * t1573 - t1647;
t1513 = t1639 * t1538 - t1636 * t1553;
t1512 = t1636 * t1538 + t1639 * t1553;
t1510 = -t1636 * t1537 + t1639 * t1543;
t1509 = t1639 * t1537 + t1636 * t1543;
t1508 = t1660 * t1603 + t1648;
t1506 = t1653 * pkin(4) + pkin(7) * t1675 - t1603 * t1650 + t1511;
t1505 = -pkin(4) * t1675 + t1653 * pkin(7) - t1650 * t1661 + t1507;
t1504 = -t1634 * t1526 + t1635 * t1527;
t1503 = t1635 * t1526 + t1634 * t1527;
t1502 = t1639 * t1519 - t1636 * t1521;
t1501 = t1636 * t1519 + t1639 * t1521;
t1500 = t1608 * pkin(4) + t1646 * pkin(7) + (t1601 * pkin(4) + t1660) * t1603 + t1648;
t1499 = t1634 * t1512 + t1635 * t1513;
t1498 = -t1635 * t1512 + t1634 * t1513;
t1497 = t1640 * t1504 + t1637 * t1554;
t1496 = t1637 * t1504 - t1640 * t1554;
t1495 = t1634 * t1509 + t1635 * t1510;
t1494 = -t1635 * t1509 + t1634 * t1510;
t1493 = t1635 * t1507 + t1634 * t1508;
t1492 = t1634 * t1507 - t1635 * t1508;
t1491 = t1640 * t1499 - t1637 * t1520;
t1490 = t1637 * t1499 + t1640 * t1520;
t1489 = t1640 * t1495 - t1637 * t1518;
t1488 = t1637 * t1495 + t1640 * t1518;
t1487 = t1634 * t1501 + t1635 * t1502;
t1486 = -t1635 * t1501 + t1634 * t1502;
t1485 = t1636 * t1500 + t1639 * t1505;
t1484 = t1639 * t1500 - t1636 * t1505;
t1483 = t1640 * t1493 + t1637 * t1511;
t1482 = t1637 * t1493 - t1640 * t1511;
t1481 = t1640 * t1487 - t1637 * t1532;
t1480 = t1637 * t1487 + t1640 * t1532;
t1479 = -t1636 * t1484 + t1639 * t1485;
t1478 = t1639 * t1484 + t1636 * t1485;
t1477 = t1634 * t1478 + t1635 * t1479;
t1476 = -t1635 * t1478 + t1634 * t1479;
t1475 = t1640 * t1477 + t1637 * t1506;
t1474 = t1637 * t1477 - t1640 * t1506;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1612, -t1611, 0, -t1638 * t1616 + t1641 * t1617, 0, 0, 0, 0, 0, 0, t1641 * t1581 - t1638 * t1609, t1641 * t1582 + t1638 * t1606, t1641 * t1610 - t1638 * t1613, t1641 * t1551 - t1638 * t1596, 0, 0, 0, 0, 0, 0, t1697, t1703, t1698, t1641 * t1497 + t1638 * t1503, 0, 0, 0, 0, 0, 0, t1697, t1698, -t1703, t1641 * t1483 + t1638 * t1492, 0, 0, 0, 0, 0, 0, t1641 * t1489 + t1638 * t1494, t1641 * t1491 + t1638 * t1498, t1641 * t1481 + t1638 * t1486, t1641 * t1475 + t1638 * t1476; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1611, t1612, 0, t1641 * t1616 + t1638 * t1617, 0, 0, 0, 0, 0, 0, t1638 * t1581 + t1641 * t1609, t1638 * t1582 - t1641 * t1606, t1638 * t1610 + t1641 * t1613, t1638 * t1551 + t1641 * t1596, 0, 0, 0, 0, 0, 0, t1699, t1704, t1700, t1638 * t1497 - t1641 * t1503, 0, 0, 0, 0, 0, 0, t1699, t1700, -t1704, t1638 * t1483 - t1641 * t1492, 0, 0, 0, 0, 0, 0, t1638 * t1489 - t1641 * t1494, t1638 * t1491 - t1641 * t1498, t1638 * t1481 - t1641 * t1486, t1638 * t1475 - t1641 * t1476; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1579, t1580, 0, t1550, 0, 0, 0, 0, 0, 0, t1692, t1528, t1690, t1496, 0, 0, 0, 0, 0, 0, t1692, t1690, -t1528, t1482, 0, 0, 0, 0, 0, 0, t1488, t1490, t1480, t1474; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1642, -qJDD(1), 0, t1617, 0, 0, 0, 0, 0, 0, t1581, t1582, t1610, t1551, 0, 0, 0, 0, 0, 0, t1691, t1530, t1689, t1497, 0, 0, 0, 0, 0, 0, t1691, t1689, -t1530, t1483, 0, 0, 0, 0, 0, 0, t1489, t1491, t1481, t1475; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1642, 0, t1616, 0, 0, 0, 0, 0, 0, t1609, -t1606, t1613, t1596, 0, 0, 0, 0, 0, 0, -t1688, -t1544, -t1687, -t1503, 0, 0, 0, 0, 0, 0, -t1688, -t1687, t1544, -t1492, 0, 0, 0, 0, 0, 0, -t1494, -t1498, -t1486, -t1476; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1579, t1580, 0, t1550, 0, 0, 0, 0, 0, 0, t1692, t1528, t1690, t1496, 0, 0, 0, 0, 0, 0, t1692, t1690, -t1528, t1482, 0, 0, 0, 0, 0, 0, t1488, t1490, t1480, t1474; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1619, t1615, t1657, t1584, 0, 0, 0, 0, 0, 0, t1686, t1546, t1685, t1504, 0, 0, 0, 0, 0, 0, t1686, t1685, -t1546, t1493, 0, 0, 0, 0, 0, 0, t1495, t1499, t1487, t1477; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1614, t1618, -t1658, t1583, 0, 0, 0, 0, 0, 0, -t1561, -t1645, t1559, -t1554, 0, 0, 0, 0, 0, 0, -t1561, t1559, t1645, -t1511, 0, 0, 0, 0, 0, 0, t1518, t1520, t1532, -t1506; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1609, t1606, -t1613, -t1596, 0, 0, 0, 0, 0, 0, t1688, t1544, t1687, t1503, 0, 0, 0, 0, 0, 0, t1688, t1687, -t1544, t1492, 0, 0, 0, 0, 0, 0, t1494, t1498, t1486, t1476; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1679, t1567, t1681, t1527, 0, 0, 0, 0, 0, 0, t1679, t1681, -t1567, t1507, 0, 0, 0, 0, 0, 0, t1510, t1513, t1502, t1479; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1568, t1680, t1646, t1526, 0, 0, 0, 0, 0, 0, -t1568, t1646, -t1680, -t1508, 0, 0, 0, 0, 0, 0, -t1509, -t1512, -t1501, -t1478; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1561, t1645, -t1559, t1554, 0, 0, 0, 0, 0, 0, t1561, -t1559, -t1645, t1511, 0, 0, 0, 0, 0, 0, -t1518, -t1520, -t1532, t1506; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1679, t1681, -t1567, t1507, 0, 0, 0, 0, 0, 0, t1510, t1513, t1502, t1479; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1561, -t1559, -t1645, t1511, 0, 0, 0, 0, 0, 0, -t1518, -t1520, -t1532, t1506; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1568, -t1646, t1680, t1508, 0, 0, 0, 0, 0, 0, t1509, t1512, t1501, t1478; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1543, t1538, t1519, t1485; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1537, t1553, t1521, t1484; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1518, t1520, t1532, -t1506;];
f_new_reg = t1;