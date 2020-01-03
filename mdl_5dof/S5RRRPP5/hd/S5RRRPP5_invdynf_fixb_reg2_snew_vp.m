% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRRPP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRRPP5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP5_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:59:04
% EndTime: 2019-12-31 20:59:07
% DurationCPUTime: 2.90s
% Computational Cost: add. (5247->218), mult. (11281->208), div. (0->0), fcn. (7430->6), ass. (0->121)
t1652 = sin(qJ(3));
t1653 = sin(qJ(2));
t1655 = cos(qJ(3));
t1656 = cos(qJ(2));
t1619 = (t1652 * t1656 + t1653 * t1655) * qJD(1);
t1616 = t1619 ^ 2;
t1692 = qJD(1) * t1656;
t1693 = qJD(1) * t1653;
t1617 = t1652 * t1693 - t1655 * t1692;
t1695 = t1617 ^ 2;
t1581 = t1616 + t1695;
t1654 = sin(qJ(1));
t1657 = cos(qJ(1));
t1648 = qJD(2) + qJD(3);
t1674 = qJD(2) * t1692;
t1678 = t1653 * qJDD(1);
t1624 = t1674 + t1678;
t1646 = t1656 * qJDD(1);
t1675 = qJD(2) * t1693;
t1625 = t1646 - t1675;
t1673 = t1652 * t1624 - t1655 * t1625;
t1556 = (qJD(3) - t1648) * t1619 + t1673;
t1671 = t1655 * t1624 + t1652 * t1625;
t1666 = -t1617 * qJD(3) + t1671;
t1691 = t1648 * t1617;
t1663 = t1666 + t1691;
t1531 = t1652 * t1556 + t1655 * t1663;
t1534 = t1655 * t1556 - t1652 * t1663;
t1704 = t1653 * t1531 - t1656 * t1534;
t1709 = t1657 * t1581 + t1654 * t1704;
t1708 = -t1654 * t1581 + t1657 * t1704;
t1635 = t1648 ^ 2;
t1600 = t1616 + t1635;
t1594 = t1619 * t1617;
t1677 = qJDD(2) + qJDD(3);
t1699 = t1594 + t1677;
t1710 = t1655 * t1600 + t1652 * t1699;
t1711 = -t1652 * t1600 + t1655 * t1699;
t1541 = t1653 * t1710 - t1656 * t1711;
t1664 = t1666 - t1691;
t1724 = t1654 * t1541 - t1657 * t1664;
t1723 = t1657 * t1541 + t1654 * t1664;
t1698 = -t1695 - t1635;
t1700 = -t1594 + t1677;
t1712 = -t1652 * t1700 + t1655 * t1698;
t1713 = t1652 * t1698 + t1655 * t1700;
t1717 = -t1653 * t1713 + t1656 * t1712;
t1722 = t1654 * t1717;
t1721 = t1657 * t1717;
t1578 = t1619 * qJD(3) + t1673;
t1705 = -t1648 * t1619 - t1578;
t1720 = t1657 * t1705 + t1722;
t1719 = -t1654 * t1705 + t1721;
t1514 = t1656 * t1531 + t1653 * t1534;
t1718 = t1653 * t1711 + t1656 * t1710;
t1716 = t1653 * t1712 + t1656 * t1713;
t1694 = 2 * qJD(4);
t1651 = t1656 ^ 2;
t1659 = qJD(1) ^ 2;
t1689 = t1651 * t1659;
t1685 = t1653 * t1659;
t1680 = qJD(3) + t1648;
t1634 = -t1657 * g(1) - t1654 * g(2);
t1665 = -t1659 * pkin(1) + qJDD(1) * pkin(6) + t1634;
t1602 = -t1653 * g(3) + t1656 * t1665;
t1670 = qJD(2) * pkin(2) - pkin(7) * t1693;
t1574 = -pkin(2) * t1689 + t1625 * pkin(7) - qJD(2) * t1670 + t1602;
t1662 = t1653 * t1665;
t1661 = -t1662 - t1624 * pkin(7) + qJDD(2) * pkin(2) + (qJD(2) * pkin(7) * qJD(1) + pkin(2) * t1685 - g(3)) * t1656;
t1544 = t1655 * t1574 + t1652 * t1661;
t1650 = t1653 ^ 2;
t1679 = t1650 + t1651;
t1633 = t1654 * g(1) - t1657 * g(2);
t1543 = -t1652 * t1574 + t1655 * t1661;
t1592 = t1617 * pkin(3) - t1619 * qJ(4);
t1672 = t1677 * qJ(4) - t1617 * t1592 + t1648 * t1694 + t1544;
t1669 = qJDD(1) * pkin(1) + t1633;
t1668 = -t1677 * pkin(3) - t1635 * qJ(4) + qJDD(4) - t1543;
t1582 = t1625 * pkin(2) - t1670 * t1693 + (t1651 * pkin(7) + pkin(6)) * t1659 + t1669;
t1660 = (-t1680 * t1617 + t1671) * qJ(4) + t1582 + t1705 * pkin(3);
t1658 = qJD(2) ^ 2;
t1639 = t1656 * t1685;
t1637 = -t1658 - t1689;
t1636 = -t1650 * t1659 - t1658;
t1632 = -qJDD(2) + t1639;
t1631 = qJDD(2) + t1639;
t1630 = t1679 * t1659;
t1629 = -t1654 * qJDD(1) - t1657 * t1659;
t1628 = t1657 * qJDD(1) - t1654 * t1659;
t1627 = t1679 * qJDD(1);
t1626 = t1646 - 0.2e1 * t1675;
t1623 = 0.2e1 * t1674 + t1678;
t1620 = t1659 * pkin(6) + t1669;
t1603 = -t1648 * pkin(4) - t1619 * qJ(5);
t1601 = -t1656 * g(3) - t1662;
t1598 = t1656 * t1632 - t1653 * t1636;
t1597 = -t1653 * t1631 + t1656 * t1637;
t1596 = t1653 * t1632 + t1656 * t1636;
t1595 = t1656 * t1631 + t1653 * t1637;
t1573 = -t1653 * t1601 + t1656 * t1602;
t1572 = t1656 * t1601 + t1653 * t1602;
t1555 = t1680 * t1619 + t1673;
t1524 = t1619 * t1592 + t1668;
t1523 = -t1635 * pkin(3) + t1672;
t1522 = t1619 * t1694 + t1660;
t1521 = -t1652 * t1543 + t1655 * t1544;
t1520 = t1655 * t1543 + t1652 * t1544;
t1513 = -t1695 * pkin(4) + t1578 * qJ(5) + 0.2e1 * qJD(5) * t1617 + (-pkin(3) * t1648 + t1603) * t1648 + t1672;
t1512 = -t1677 * pkin(4) + (pkin(4) * t1617 - 0.2e1 * qJD(5) + t1592) * t1619 + t1668 - t1663 * qJ(5);
t1511 = qJDD(5) - t1578 * pkin(4) - t1695 * qJ(5) + (t1694 + t1603) * t1619 + t1660;
t1510 = t1655 * t1523 + t1652 * t1524;
t1509 = t1652 * t1523 - t1655 * t1524;
t1508 = -t1653 * t1520 + t1656 * t1521;
t1507 = t1656 * t1520 + t1653 * t1521;
t1506 = t1652 * t1512 + t1655 * t1513;
t1505 = -t1655 * t1512 + t1652 * t1513;
t1504 = -t1653 * t1509 + t1656 * t1510;
t1503 = t1656 * t1509 + t1653 * t1510;
t1502 = -t1653 * t1505 + t1656 * t1506;
t1501 = t1656 * t1505 + t1653 * t1506;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1629, -t1628, 0, -t1654 * t1633 + t1657 * t1634, 0, 0, 0, 0, 0, 0, t1657 * t1597 - t1654 * t1626, t1657 * t1598 + t1654 * t1623, t1657 * t1627 - t1654 * t1630, t1657 * t1573 - t1654 * t1620, 0, 0, 0, 0, 0, 0, t1654 * t1555 + t1721, t1723, t1708, t1657 * t1508 - t1654 * t1582, 0, 0, 0, 0, 0, 0, t1719, t1708, -t1723, t1657 * t1504 - t1654 * t1522, 0, 0, 0, 0, 0, 0, t1719, -t1723, -t1708, t1657 * t1502 - t1654 * t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1628, t1629, 0, t1657 * t1633 + t1654 * t1634, 0, 0, 0, 0, 0, 0, t1654 * t1597 + t1657 * t1626, t1654 * t1598 - t1657 * t1623, t1654 * t1627 + t1657 * t1630, t1654 * t1573 + t1657 * t1620, 0, 0, 0, 0, 0, 0, -t1657 * t1555 + t1722, t1724, t1709, t1654 * t1508 + t1657 * t1582, 0, 0, 0, 0, 0, 0, t1720, t1709, -t1724, t1654 * t1504 + t1657 * t1522, 0, 0, 0, 0, 0, 0, t1720, -t1724, -t1709, t1654 * t1502 + t1657 * t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1595, t1596, 0, t1572, 0, 0, 0, 0, 0, 0, t1716, -t1718, -t1514, t1507, 0, 0, 0, 0, 0, 0, t1716, -t1514, t1718, t1503, 0, 0, 0, 0, 0, 0, t1716, t1718, t1514, t1501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1659, -qJDD(1), 0, t1634, 0, 0, 0, 0, 0, 0, t1597, t1598, t1627, t1573, 0, 0, 0, 0, 0, 0, t1717, t1541, t1704, t1508, 0, 0, 0, 0, 0, 0, t1717, t1704, -t1541, t1504, 0, 0, 0, 0, 0, 0, t1717, -t1541, -t1704, t1502; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1659, 0, t1633, 0, 0, 0, 0, 0, 0, t1626, -t1623, t1630, t1620, 0, 0, 0, 0, 0, 0, -t1555, -t1664, t1581, t1582, 0, 0, 0, 0, 0, 0, t1705, t1581, t1664, t1522, 0, 0, 0, 0, 0, 0, t1705, t1664, -t1581, t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1595, t1596, 0, t1572, 0, 0, 0, 0, 0, 0, t1716, -t1718, -t1514, t1507, 0, 0, 0, 0, 0, 0, t1716, -t1514, t1718, t1503, 0, 0, 0, 0, 0, 0, t1716, t1718, t1514, t1501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1637, t1632, t1646, t1602, 0, 0, 0, 0, 0, 0, t1712, -t1711, -t1534, t1521, 0, 0, 0, 0, 0, 0, t1712, -t1534, t1711, t1510, 0, 0, 0, 0, 0, 0, t1712, t1711, t1534, t1506; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1631, t1636, -t1678, t1601, 0, 0, 0, 0, 0, 0, t1713, -t1710, -t1531, t1520, 0, 0, 0, 0, 0, 0, t1713, -t1531, t1710, t1509, 0, 0, 0, 0, 0, 0, t1713, t1710, t1531, t1505; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1626, t1623, -t1630, -t1620, 0, 0, 0, 0, 0, 0, t1555, t1664, -t1581, -t1582, 0, 0, 0, 0, 0, 0, -t1705, -t1581, -t1664, -t1522, 0, 0, 0, 0, 0, 0, -t1705, -t1664, t1581, -t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1698, -t1699, -t1556, t1544, 0, 0, 0, 0, 0, 0, t1698, -t1556, t1699, t1523, 0, 0, 0, 0, 0, 0, t1698, t1699, t1556, t1513; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1700, -t1600, -t1663, t1543, 0, 0, 0, 0, 0, 0, t1700, -t1663, t1600, -t1524, 0, 0, 0, 0, 0, 0, t1700, t1600, t1663, -t1512; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1555, t1664, -t1581, -t1582, 0, 0, 0, 0, 0, 0, -t1705, -t1581, -t1664, -t1522, 0, 0, 0, 0, 0, 0, -t1705, -t1664, t1581, -t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1698, -t1556, t1699, t1523, 0, 0, 0, 0, 0, 0, t1698, t1699, t1556, t1513; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1705, -t1581, -t1664, -t1522, 0, 0, 0, 0, 0, 0, -t1705, -t1664, t1581, -t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1700, t1663, -t1600, t1524, 0, 0, 0, 0, 0, 0, -t1700, -t1600, -t1663, t1512; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1698, t1699, t1556, t1513; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1700, -t1600, -t1663, t1512; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1705, t1664, -t1581, t1511;];
f_new_reg = t1;
