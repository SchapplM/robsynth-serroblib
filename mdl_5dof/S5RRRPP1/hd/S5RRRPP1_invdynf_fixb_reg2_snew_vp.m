% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRRPP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynf_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:50:36
% EndTime: 2019-12-31 20:50:39
% DurationCPUTime: 3.08s
% Computational Cost: add. (9196->191), mult. (12649->236), div. (0->0), fcn. (8077->8), ass. (0->135)
t1634 = sin(pkin(8));
t1635 = cos(pkin(8));
t1630 = qJD(1) + qJD(2);
t1639 = cos(qJ(3));
t1662 = t1630 * t1639;
t1636 = sin(qJ(3));
t1663 = t1630 * t1636;
t1595 = t1634 * t1663 - t1635 * t1662;
t1597 = (t1634 * t1639 + t1635 * t1636) * t1630;
t1665 = t1597 * t1595;
t1570 = qJDD(3) + t1665;
t1594 = t1597 ^ 2;
t1642 = qJD(3) ^ 2;
t1672 = -t1594 - t1642;
t1535 = t1634 * t1570 - t1635 * t1672;
t1537 = t1635 * t1570 + t1634 * t1672;
t1525 = t1636 * t1535 - t1639 * t1537;
t1637 = sin(qJ(2));
t1640 = cos(qJ(2));
t1593 = qJD(3) * t1595;
t1652 = qJD(3) * t1662;
t1629 = qJDD(1) + qJDD(2);
t1658 = t1636 * t1629;
t1601 = t1652 + t1658;
t1653 = qJD(3) * t1663;
t1656 = t1639 * t1629;
t1648 = -t1653 + t1656;
t1645 = t1635 * t1601 + t1634 * t1648;
t1671 = -t1645 + t1593;
t1513 = t1637 * t1525 + t1640 * t1671;
t1515 = t1640 * t1525 - t1637 * t1671;
t1638 = sin(qJ(1));
t1641 = cos(qJ(1));
t1695 = t1641 * t1513 + t1638 * t1515;
t1694 = t1638 * t1513 - t1641 * t1515;
t1650 = -t1634 * t1601 + t1635 * t1648;
t1666 = qJD(3) * t1597;
t1553 = -t1650 + t1666;
t1571 = qJDD(3) - t1665;
t1577 = t1595 ^ 2;
t1673 = -t1577 - t1642;
t1678 = -t1634 * t1571 + t1635 * t1673;
t1681 = t1635 * t1571 + t1634 * t1673;
t1683 = -t1636 * t1681 + t1639 * t1678;
t1690 = t1637 * t1553 + t1640 * t1683;
t1691 = -t1640 * t1553 + t1637 * t1683;
t1693 = -t1638 * t1691 + t1641 * t1690;
t1692 = t1638 * t1690 + t1641 * t1691;
t1519 = t1639 * t1535 + t1636 * t1537;
t1552 = t1594 + t1577;
t1558 = t1593 + t1645;
t1647 = t1650 + t1666;
t1669 = t1634 * t1558 + t1635 * t1647;
t1670 = -t1635 * t1558 + t1634 * t1647;
t1677 = -t1636 * t1670 + t1639 * t1669;
t1684 = -t1637 * t1552 + t1640 * t1677;
t1685 = t1640 * t1552 + t1637 * t1677;
t1689 = -t1638 * t1685 + t1641 * t1684;
t1688 = t1638 * t1684 + t1641 * t1685;
t1682 = t1636 * t1678 + t1639 * t1681;
t1628 = t1630 ^ 2;
t1609 = t1637 * t1628 - t1640 * t1629;
t1649 = -t1640 * t1628 - t1637 * t1629;
t1680 = t1638 * t1609 + t1641 * t1649;
t1679 = t1641 * t1609 - t1638 * t1649;
t1676 = t1636 * t1669 + t1639 * t1670;
t1668 = t1639 ^ 2;
t1667 = -2 * qJD(4);
t1664 = t1628 * t1636;
t1623 = -t1641 * g(1) - t1638 * g(2);
t1643 = qJD(1) ^ 2;
t1612 = -t1643 * pkin(1) + t1623;
t1622 = t1638 * g(1) - t1641 * g(2);
t1646 = qJDD(1) * pkin(1) + t1622;
t1581 = t1640 * t1612 + t1637 * t1646;
t1574 = -t1628 * pkin(2) + t1629 * pkin(7) + t1581;
t1659 = t1636 * t1574;
t1655 = t1668 * t1628;
t1632 = t1636 ^ 2;
t1654 = t1632 + t1668;
t1564 = -t1636 * g(3) + t1639 * t1574;
t1613 = qJD(3) * pkin(3) - qJ(4) * t1663;
t1541 = -pkin(3) * t1655 + t1648 * qJ(4) - qJD(3) * t1613 + t1564;
t1644 = qJDD(3) * pkin(3) - t1601 * qJ(4) - t1659 + (qJD(3) * t1630 * qJ(4) + pkin(3) * t1664 - g(3)) * t1639;
t1517 = t1635 * t1541 + t1595 * t1667 + t1634 * t1644;
t1651 = t1634 * t1541 - t1635 * t1644;
t1580 = -t1637 * t1612 + t1640 * t1646;
t1573 = -t1629 * pkin(2) - t1628 * pkin(7) - t1580;
t1542 = -t1648 * pkin(3) - qJ(4) * t1655 + t1613 * t1663 + qJDD(4) + t1573;
t1621 = t1639 * t1664;
t1620 = -t1642 - t1655;
t1619 = -t1632 * t1628 - t1642;
t1617 = -t1638 * qJDD(1) - t1641 * t1643;
t1616 = t1641 * qJDD(1) - t1638 * t1643;
t1615 = -qJDD(3) + t1621;
t1614 = qJDD(3) + t1621;
t1611 = t1654 * t1628;
t1606 = t1654 * t1629;
t1602 = -0.2e1 * t1653 + t1656;
t1600 = 0.2e1 * t1652 + t1658;
t1585 = t1639 * t1615 - t1636 * t1619;
t1584 = -t1636 * t1614 + t1639 * t1620;
t1583 = t1636 * t1615 + t1639 * t1619;
t1582 = t1639 * t1614 + t1636 * t1620;
t1579 = t1640 * t1606 - t1637 * t1611;
t1578 = t1637 * t1606 + t1640 * t1611;
t1567 = t1595 * pkin(4) - t1597 * qJ(5);
t1563 = -t1639 * g(3) - t1659;
t1562 = t1640 * t1585 + t1637 * t1600;
t1561 = t1640 * t1584 - t1637 * t1602;
t1560 = t1637 * t1585 - t1640 * t1600;
t1559 = t1637 * t1584 + t1640 * t1602;
t1548 = -t1637 * t1580 + t1640 * t1581;
t1547 = t1640 * t1580 + t1637 * t1581;
t1533 = -t1636 * t1563 + t1639 * t1564;
t1532 = t1639 * t1563 + t1636 * t1564;
t1527 = t1640 * t1533 + t1637 * t1573;
t1526 = t1637 * t1533 - t1640 * t1573;
t1516 = t1597 * t1667 - t1651;
t1511 = -t1650 * pkin(4) + (pkin(4) * qJD(3) - (2 * qJD(5))) * t1597 + t1542 + t1671 * qJ(5);
t1502 = qJDD(5) - t1642 * qJ(5) - qJDD(3) * pkin(4) + ((2 * qJD(4)) + t1567) * t1597 + t1651;
t1501 = -t1642 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t1595 * t1567 + t1517;
t1496 = -t1634 * t1516 + t1635 * t1517;
t1495 = t1635 * t1516 + t1634 * t1517;
t1494 = t1635 * t1501 + t1634 * t1502;
t1493 = t1634 * t1501 - t1635 * t1502;
t1492 = -t1636 * t1495 + t1639 * t1496;
t1491 = t1639 * t1495 + t1636 * t1496;
t1490 = t1640 * t1492 + t1637 * t1542;
t1489 = t1637 * t1492 - t1640 * t1542;
t1488 = -t1636 * t1493 + t1639 * t1494;
t1487 = t1639 * t1493 + t1636 * t1494;
t1486 = t1640 * t1488 + t1637 * t1511;
t1485 = t1637 * t1488 - t1640 * t1511;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1617, -t1616, 0, -t1638 * t1622 + t1641 * t1623, 0, 0, 0, 0, 0, 0, t1680, t1679, 0, -t1638 * t1547 + t1641 * t1548, 0, 0, 0, 0, 0, 0, -t1638 * t1559 + t1641 * t1561, -t1638 * t1560 + t1641 * t1562, -t1638 * t1578 + t1641 * t1579, -t1638 * t1526 + t1641 * t1527, 0, 0, 0, 0, 0, 0, t1693, -t1694, t1689, -t1638 * t1489 + t1641 * t1490, 0, 0, 0, 0, 0, 0, t1693, t1689, t1694, -t1638 * t1485 + t1641 * t1486; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1616, t1617, 0, t1641 * t1622 + t1638 * t1623, 0, 0, 0, 0, 0, 0, -t1679, t1680, 0, t1641 * t1547 + t1638 * t1548, 0, 0, 0, 0, 0, 0, t1641 * t1559 + t1638 * t1561, t1641 * t1560 + t1638 * t1562, t1641 * t1578 + t1638 * t1579, t1641 * t1526 + t1638 * t1527, 0, 0, 0, 0, 0, 0, t1692, t1695, t1688, t1641 * t1489 + t1638 * t1490, 0, 0, 0, 0, 0, 0, t1692, t1688, -t1695, t1641 * t1485 + t1638 * t1486; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1582, t1583, 0, t1532, 0, 0, 0, 0, 0, 0, t1682, -t1519, t1676, t1491, 0, 0, 0, 0, 0, 0, t1682, t1676, t1519, t1487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1643, -qJDD(1), 0, t1623, 0, 0, 0, 0, 0, 0, t1649, t1609, 0, t1548, 0, 0, 0, 0, 0, 0, t1561, t1562, t1579, t1527, 0, 0, 0, 0, 0, 0, t1690, t1515, t1684, t1490, 0, 0, 0, 0, 0, 0, t1690, t1684, -t1515, t1486; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1643, 0, t1622, 0, 0, 0, 0, 0, 0, -t1609, t1649, 0, t1547, 0, 0, 0, 0, 0, 0, t1559, t1560, t1578, t1526, 0, 0, 0, 0, 0, 0, t1691, t1513, t1685, t1489, 0, 0, 0, 0, 0, 0, t1691, t1685, -t1513, t1485; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1582, t1583, 0, t1532, 0, 0, 0, 0, 0, 0, t1682, -t1519, t1676, t1491, 0, 0, 0, 0, 0, 0, t1682, t1676, t1519, t1487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1628, -t1629, 0, t1581, 0, 0, 0, 0, 0, 0, t1584, t1585, t1606, t1533, 0, 0, 0, 0, 0, 0, t1683, t1525, t1677, t1492, 0, 0, 0, 0, 0, 0, t1683, t1677, -t1525, t1488; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1629, -t1628, 0, t1580, 0, 0, 0, 0, 0, 0, t1602, -t1600, t1611, -t1573, 0, 0, 0, 0, 0, 0, -t1553, t1671, t1552, -t1542, 0, 0, 0, 0, 0, 0, -t1553, t1552, -t1671, -t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1582, t1583, 0, t1532, 0, 0, 0, 0, 0, 0, t1682, -t1519, t1676, t1491, 0, 0, 0, 0, 0, 0, t1682, t1676, t1519, t1487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1620, t1615, t1656, t1564, 0, 0, 0, 0, 0, 0, t1678, -t1537, t1669, t1496, 0, 0, 0, 0, 0, 0, t1678, t1669, t1537, t1494; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1614, t1619, -t1658, t1563, 0, 0, 0, 0, 0, 0, t1681, -t1535, t1670, t1495, 0, 0, 0, 0, 0, 0, t1681, t1670, t1535, t1493; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1602, t1600, -t1611, t1573, 0, 0, 0, 0, 0, 0, t1553, -t1671, -t1552, t1542, 0, 0, 0, 0, 0, 0, t1553, -t1552, t1671, t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1673, -t1570, t1647, t1517, 0, 0, 0, 0, 0, 0, t1673, t1647, t1570, t1501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1571, t1672, -t1558, t1516, 0, 0, 0, 0, 0, 0, t1571, -t1558, -t1672, -t1502; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1553, -t1671, -t1552, t1542, 0, 0, 0, 0, 0, 0, t1553, -t1552, t1671, t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1673, t1647, t1570, t1501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1553, -t1552, t1671, t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1571, t1558, t1672, t1502;];
f_new_reg = t1;
