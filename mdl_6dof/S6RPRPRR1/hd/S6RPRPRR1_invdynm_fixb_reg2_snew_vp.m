% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RPRPRR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 18:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RPRPRR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_invdynm_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:15:42
% EndTime: 2019-05-05 18:16:22
% DurationCPUTime: 42.70s
% Computational Cost: add. (377201->941), mult. (839945->1355), div. (0->0), fcn. (603971->12), ass. (0->657)
t1697 = cos(qJ(1));
t1693 = sin(qJ(1));
t1660 = g(1) * t1697 + g(2) * t1693;
t1808 = qJD(1) ^ 2;
t1643 = -pkin(1) * t1808 - t1660;
t1686 = sin(pkin(10));
t1688 = cos(pkin(10));
t1659 = g(1) * t1693 - t1697 * g(2);
t1708 = qJDD(1) * pkin(1) + t1659;
t1591 = t1686 * t1643 - t1688 * t1708;
t1592 = t1688 * t1643 + t1686 * t1708;
t1742 = t1591 * t1686 + t1688 * t1592;
t1514 = t1591 * t1688 - t1592 * t1686;
t1784 = t1514 * t1693;
t1828 = t1697 * t1742 + t1784;
t1783 = t1514 * t1697;
t1827 = -t1693 * t1742 + t1783;
t1685 = sin(pkin(11));
t1687 = cos(pkin(11));
t1696 = cos(qJ(3));
t1799 = qJD(1) * t1696;
t1692 = sin(qJ(3));
t1800 = qJD(1) * t1692;
t1631 = -t1685 * t1800 + t1687 * t1799;
t1632 = (t1685 * t1696 + t1687 * t1692) * qJD(1);
t1590 = t1631 * t1632;
t1810 = qJDD(3) + t1590;
t1826 = t1685 * t1810;
t1825 = t1687 * t1810;
t1690 = sin(qJ(6));
t1691 = sin(qJ(5));
t1695 = cos(qJ(5));
t1571 = t1631 * t1691 + t1632 * t1695;
t1672 = qJD(3) * t1799;
t1764 = t1692 * qJDD(1);
t1646 = t1672 + t1764;
t1753 = qJD(3) * t1800;
t1763 = t1696 * qJDD(1);
t1706 = -t1753 + t1763;
t1593 = t1687 * t1646 + t1685 * t1706;
t1740 = t1646 * t1685 - t1687 * t1706;
t1741 = t1691 * t1593 + t1695 * t1740;
t1471 = -t1571 * qJD(5) - t1741;
t1470 = qJDD(6) - t1471;
t1694 = cos(qJ(6));
t1769 = qJD(3) + qJD(5);
t1535 = t1571 * t1690 - t1694 * t1769;
t1537 = t1694 * t1571 + t1690 * t1769;
t1484 = t1537 * t1535;
t1813 = t1470 - t1484;
t1824 = t1690 * t1813;
t1569 = -t1695 * t1631 + t1632 * t1691;
t1504 = t1571 * t1569;
t1681 = qJDD(3) + qJDD(5);
t1812 = -t1504 + t1681;
t1823 = t1691 * t1812;
t1822 = t1694 * t1813;
t1821 = t1695 * t1812;
t1472 = -t1569 * qJD(5) + t1695 * t1593 - t1691 * t1740;
t1558 = t1769 * t1569;
t1820 = -t1558 + t1472;
t1649 = qJDD(1) * t1686 + t1688 * t1808;
t1683 = g(3) - qJDD(2);
t1620 = qJ(2) * t1649 - t1683 * t1688;
t1650 = qJDD(1) * t1688 - t1686 * t1808;
t1714 = -qJ(2) * t1650 - t1683 * t1686;
t1811 = t1697 * t1649 + t1650 * t1693;
t1819 = pkin(6) * t1811 + t1697 * t1620 - t1693 * t1714;
t1596 = -t1649 * t1693 + t1697 * t1650;
t1818 = -pkin(6) * t1596 + t1693 * t1620 + t1697 * t1714;
t1762 = t1769 ^ 2;
t1797 = qJD(4) * t1631;
t1625 = 0.2e1 * t1797;
t1573 = -pkin(2) * t1808 + qJDD(1) * pkin(7) + t1592;
t1545 = t1692 * t1573 + t1696 * t1683;
t1666 = t1696 * t1808 * t1692;
t1657 = qJDD(3) + t1666;
t1500 = (-t1646 + t1672) * qJ(4) + t1657 * pkin(3) - t1545;
t1547 = t1696 * t1573 - t1692 * t1683;
t1656 = qJD(3) * pkin(3) - qJ(4) * t1800;
t1806 = t1696 ^ 2;
t1679 = t1806 * t1808;
t1503 = -pkin(3) * t1679 + qJ(4) * t1706 - qJD(3) * t1656 + t1547;
t1767 = t1685 * t1500 + t1687 * t1503;
t1408 = t1625 + t1767;
t1616 = qJD(3) * pkin(4) - pkin(8) * t1632;
t1807 = t1631 ^ 2;
t1377 = -pkin(4) * t1807 - pkin(8) * t1740 - qJD(3) * t1616 + t1408;
t1407 = 0.2e1 * qJD(4) * t1632 - t1687 * t1500 + t1685 * t1503;
t1627 = qJD(3) * t1631;
t1543 = -t1627 + t1593;
t1699 = pkin(4) * t1810 - pkin(8) * t1543 - t1407;
t1291 = t1695 * t1377 + t1691 * t1699;
t1499 = pkin(5) * t1569 - pkin(9) * t1571;
t1258 = -pkin(5) * t1762 + t1681 * pkin(9) - t1569 * t1499 + t1291;
t1572 = -qJDD(1) * pkin(2) - pkin(7) * t1808 + t1591;
t1511 = -t1706 * pkin(3) - qJ(4) * t1679 + t1656 * t1800 + qJDD(4) + t1572;
t1418 = t1740 * pkin(4) - t1807 * pkin(8) + t1632 * t1616 + t1511;
t1744 = t1769 * t1571;
t1305 = t1418 + (-t1471 + t1744) * pkin(5) - t1820 * pkin(9);
t1214 = t1258 * t1690 - t1694 * t1305;
t1215 = t1258 * t1694 + t1305 * t1690;
t1154 = t1690 * t1214 + t1694 * t1215;
t1480 = t1692 * t1545 + t1696 * t1547;
t1533 = t1535 ^ 2;
t1534 = t1537 ^ 2;
t1565 = qJD(6) + t1569;
t1564 = t1565 ^ 2;
t1566 = t1569 ^ 2;
t1567 = t1571 ^ 2;
t1630 = t1632 ^ 2;
t1316 = -t1407 * t1687 + t1408 * t1685;
t1805 = pkin(3) * t1316;
t1798 = qJD(3) * t1632;
t1541 = -t1740 + t1798;
t1474 = t1541 * t1685 - t1543 * t1687;
t1804 = pkin(3) * t1474;
t1803 = pkin(5) * t1691;
t1801 = qJD(1) * qJD(3);
t1290 = t1377 * t1691 - t1695 * t1699;
t1218 = -t1290 * t1695 + t1291 * t1691;
t1796 = t1218 * t1685;
t1795 = t1218 * t1687;
t1794 = t1316 * t1692;
t1793 = t1316 * t1696;
t1388 = t1470 + t1484;
t1792 = t1388 * t1690;
t1791 = t1388 * t1694;
t1790 = t1418 * t1691;
t1789 = t1418 * t1695;
t1496 = t1504 + t1681;
t1788 = t1496 * t1691;
t1787 = t1496 * t1695;
t1786 = t1511 * t1685;
t1785 = t1511 * t1687;
t1782 = t1565 * t1690;
t1781 = t1565 * t1694;
t1585 = qJDD(3) - t1590;
t1780 = t1585 * t1685;
t1779 = t1585 * t1687;
t1778 = t1631 * t1685;
t1777 = t1631 * t1687;
t1776 = t1632 * t1685;
t1775 = t1632 * t1687;
t1647 = -0.2e1 * t1753 + t1763;
t1602 = t1647 * t1696;
t1772 = t1657 * t1692;
t1658 = qJDD(3) - t1666;
t1771 = t1658 * t1692;
t1770 = t1658 * t1696;
t1257 = -t1681 * pkin(5) - pkin(9) * t1762 + t1499 * t1571 + t1290;
t1254 = t1690 * t1257;
t1561 = t1692 * t1572;
t1255 = t1694 * t1257;
t1562 = t1696 * t1572;
t1768 = -pkin(5) * t1257 + pkin(9) * t1154;
t1766 = -pkin(2) * t1572 + pkin(7) * t1480;
t1765 = qJDD(3) * t1688;
t1682 = t1692 ^ 2;
t1761 = t1682 + t1806;
t1760 = -pkin(5) * t1695 - pkin(4);
t1759 = t1691 * t1484;
t1758 = t1695 * t1484;
t1757 = t1686 * t1504;
t1756 = t1688 * t1504;
t1755 = t1686 * t1590;
t1754 = t1688 * t1590;
t1455 = -t1534 - t1564;
t1332 = -t1455 * t1690 - t1791;
t1712 = -t1694 * t1472 - t1690 * t1681;
t1373 = (qJD(6) + t1565) * t1535 + t1712;
t1752 = pkin(5) * t1373 + pkin(9) * t1332 + t1254;
t1443 = -t1564 - t1533;
t1326 = t1443 * t1694 - t1824;
t1492 = t1565 * t1537;
t1743 = -t1690 * t1472 + t1694 * t1681;
t1707 = qJD(6) * t1537 - t1743;
t1369 = -t1492 - t1707;
t1751 = pkin(5) * t1369 + pkin(9) * t1326 - t1255;
t1123 = t1154 * t1691 - t1257 * t1695;
t1750 = pkin(4) * t1123 + t1768;
t1548 = -t1567 - t1762;
t1448 = t1548 * t1695 - t1788;
t1749 = pkin(4) * t1448 - t1291;
t1678 = t1682 * t1808;
t1698 = qJD(3) ^ 2;
t1663 = -t1678 - t1698;
t1611 = -t1663 * t1692 - t1770;
t1645 = 0.2e1 * t1672 + t1764;
t1748 = -pkin(2) * t1645 + pkin(7) * t1611 + t1561;
t1665 = -t1679 - t1698;
t1609 = t1665 * t1696 - t1772;
t1747 = pkin(2) * t1647 + pkin(7) * t1609 - t1562;
t1219 = t1290 * t1691 + t1695 * t1291;
t1156 = t1219 * t1685 + t1795;
t1217 = pkin(4) * t1218;
t1746 = pkin(3) * t1156 + t1217;
t1434 = qJD(3) * t1571 - t1741;
t1437 = t1558 + t1472;
t1341 = t1434 * t1691 - t1437 * t1695;
t1343 = t1434 * t1695 + t1437 * t1691;
t1266 = t1341 * t1687 + t1343 * t1685;
t1339 = pkin(4) * t1341;
t1745 = pkin(3) * t1266 + t1339;
t1317 = t1407 * t1685 + t1687 * t1408;
t1739 = -t1659 * t1693 - t1697 * t1660;
t1738 = t1686 * t1666;
t1737 = t1688 * t1666;
t1124 = t1154 * t1695 + t1257 * t1691;
t1153 = -t1214 * t1694 + t1215 * t1690;
t1073 = pkin(8) * t1124 + (-pkin(9) * t1691 + t1760) * t1153;
t1085 = -pkin(8) * t1123 + (-pkin(9) * t1695 + t1803) * t1153;
t1088 = -t1123 * t1685 + t1124 * t1687;
t1049 = -pkin(3) * t1153 + qJ(4) * t1088 + t1073 * t1687 + t1085 * t1685;
t1087 = t1123 * t1687 + t1124 * t1685;
t1054 = -qJ(4) * t1087 - t1073 * t1685 + t1085 * t1687;
t1067 = -t1087 * t1692 + t1088 * t1696;
t1736 = -pkin(2) * t1153 + pkin(7) * t1067 + t1696 * t1049 + t1692 * t1054;
t1370 = (-qJD(6) + t1565) * t1537 + t1743;
t1411 = -qJD(6) * t1535 - t1712;
t1491 = t1565 * t1535;
t1372 = t1411 + t1491;
t1297 = t1370 * t1690 - t1372 * t1694;
t1130 = -pkin(9) * t1297 - t1153;
t1299 = t1370 * t1694 + t1372 * t1690;
t1439 = t1533 + t1534;
t1253 = t1299 * t1695 - t1439 * t1691;
t1096 = pkin(8) * t1253 + t1691 * t1130 + t1297 * t1760;
t1252 = t1299 * t1691 + t1439 * t1695;
t1106 = -pkin(8) * t1252 + t1130 * t1695 + t1297 * t1803;
t1181 = -t1252 * t1685 + t1253 * t1687;
t1069 = -pkin(3) * t1297 + qJ(4) * t1181 + t1096 * t1687 + t1106 * t1685;
t1180 = t1252 * t1687 + t1253 * t1685;
t1071 = -qJ(4) * t1180 - t1096 * t1685 + t1106 * t1687;
t1133 = -t1180 * t1692 + t1181 * t1696;
t1735 = -pkin(2) * t1297 + pkin(7) * t1133 + t1696 * t1069 + t1692 * t1071;
t1325 = t1443 * t1690 + t1822;
t1173 = -pkin(5) * t1325 + t1214;
t1223 = -pkin(9) * t1325 + t1254;
t1271 = t1326 * t1695 - t1369 * t1691;
t1113 = -pkin(4) * t1325 + pkin(8) * t1271 + t1173 * t1695 + t1223 * t1691;
t1270 = t1326 * t1691 + t1369 * t1695;
t1118 = -pkin(8) * t1270 - t1173 * t1691 + t1223 * t1695;
t1195 = -t1270 * t1685 + t1271 * t1687;
t1076 = -pkin(3) * t1325 + qJ(4) * t1195 + t1113 * t1687 + t1118 * t1685;
t1194 = t1270 * t1687 + t1271 * t1685;
t1080 = -qJ(4) * t1194 - t1113 * t1685 + t1118 * t1687;
t1142 = -t1194 * t1692 + t1195 * t1696;
t1734 = -pkin(2) * t1325 + pkin(7) * t1142 + t1696 * t1076 + t1692 * t1080;
t1331 = t1455 * t1694 - t1792;
t1174 = -pkin(5) * t1331 + t1215;
t1224 = -pkin(9) * t1331 + t1255;
t1274 = t1332 * t1695 - t1373 * t1691;
t1114 = -pkin(4) * t1331 + pkin(8) * t1274 + t1174 * t1695 + t1224 * t1691;
t1273 = t1332 * t1691 + t1373 * t1695;
t1119 = -pkin(8) * t1273 - t1174 * t1691 + t1224 * t1695;
t1197 = -t1273 * t1685 + t1274 * t1687;
t1078 = -pkin(3) * t1331 + qJ(4) * t1197 + t1114 * t1687 + t1119 * t1685;
t1196 = t1273 * t1687 + t1274 * t1685;
t1082 = -qJ(4) * t1196 - t1114 * t1685 + t1119 * t1687;
t1146 = -t1196 * t1692 + t1197 * t1696;
t1733 = -pkin(2) * t1331 + pkin(7) * t1146 + t1696 * t1078 + t1692 * t1082;
t1157 = t1219 * t1687 - t1796;
t1201 = -pkin(4) * t1418 + pkin(8) * t1219;
t1095 = -pkin(3) * t1418 - pkin(8) * t1796 + qJ(4) * t1157 + t1201 * t1687;
t1100 = -pkin(8) * t1795 - qJ(4) * t1156 - t1201 * t1685;
t1104 = -t1156 * t1692 + t1157 * t1696;
t1732 = -pkin(2) * t1418 + pkin(7) * t1104 + t1696 * t1095 + t1692 * t1100;
t1469 = -t1566 - t1567;
t1172 = -pkin(4) * t1469 + pkin(8) * t1343 + t1219;
t1179 = -pkin(8) * t1341 - t1218;
t1268 = -t1341 * t1685 + t1343 * t1687;
t1109 = -pkin(3) * t1469 + qJ(4) * t1268 + t1172 * t1687 + t1179 * t1685;
t1112 = -qJ(4) * t1266 - t1172 * t1685 + t1179 * t1687;
t1193 = -t1266 * t1692 + t1268 * t1696;
t1731 = -pkin(2) * t1469 + pkin(7) * t1193 + t1696 * t1109 + t1692 * t1112;
t1494 = -t1762 - t1566;
t1413 = t1494 * t1695 - t1823;
t1432 = (0.2e1 * qJD(5) + qJD(3)) * t1571 + t1741;
t1300 = -pkin(4) * t1432 + pkin(8) * t1413 - t1789;
t1412 = t1494 * t1691 + t1821;
t1328 = -t1412 * t1685 + t1413 * t1687;
t1333 = -pkin(8) * t1412 + t1790;
t1186 = -pkin(3) * t1432 + qJ(4) * t1328 + t1300 * t1687 + t1333 * t1685;
t1327 = t1412 * t1687 + t1413 * t1685;
t1205 = -qJ(4) * t1327 - t1300 * t1685 + t1333 * t1687;
t1249 = -t1327 * t1692 + t1328 * t1696;
t1730 = -pkin(2) * t1432 + pkin(7) * t1249 + t1696 * t1186 + t1692 * t1205;
t1729 = pkin(5) * t1439 + pkin(9) * t1299 + t1154;
t1449 = -t1548 * t1691 - t1787;
t1307 = -pkin(4) * t1820 + pkin(8) * t1449 + t1790;
t1350 = -pkin(8) * t1448 + t1789;
t1353 = -t1448 * t1685 + t1449 * t1687;
t1203 = -pkin(3) * t1820 + qJ(4) * t1353 + t1307 * t1687 + t1350 * t1685;
t1352 = t1448 * t1687 + t1449 * t1685;
t1221 = -qJ(4) * t1352 - t1307 * t1685 + t1350 * t1687;
t1283 = -t1352 * t1692 + t1353 * t1696;
t1728 = -pkin(2) * t1820 + pkin(7) * t1283 + t1696 * t1203 + t1692 * t1221;
t1727 = pkin(4) * t1273 + t1752;
t1726 = pkin(4) * t1270 + t1751;
t1476 = t1541 * t1687 + t1543 * t1685;
t1538 = -t1630 - t1807;
t1292 = -pkin(3) * t1538 + qJ(4) * t1476 + t1317;
t1302 = -qJ(4) * t1474 - t1316;
t1386 = -t1474 * t1692 + t1476 * t1696;
t1725 = -pkin(2) * t1538 + pkin(7) * t1386 + t1696 * t1292 + t1692 * t1302;
t1583 = -t1698 - t1807;
t1506 = t1583 * t1687 - t1826;
t1539 = t1740 + t1798;
t1391 = -pkin(3) * t1539 + qJ(4) * t1506 - t1785;
t1505 = t1583 * t1685 + t1825;
t1416 = -t1505 * t1692 + t1506 * t1696;
t1429 = -qJ(4) * t1505 + t1786;
t1724 = -pkin(2) * t1539 + pkin(7) * t1416 + t1696 * t1391 + t1692 * t1429;
t1623 = -t1630 - t1698;
t1523 = -t1623 * t1685 - t1779;
t1542 = t1627 + t1593;
t1397 = -pkin(3) * t1542 + qJ(4) * t1523 + t1786;
t1520 = t1623 * t1687 - t1780;
t1444 = -qJ(4) * t1520 + t1785;
t1461 = -t1520 * t1692 + t1523 * t1696;
t1723 = -pkin(2) * t1542 + pkin(7) * t1461 + t1696 * t1397 + t1692 * t1444;
t1651 = t1761 * qJDD(1);
t1654 = t1678 + t1679;
t1722 = pkin(2) * t1654 + pkin(7) * t1651 + t1480;
t1721 = pkin(3) * t1520 - t1767;
t1653 = qJDD(1) * t1697 - t1693 * t1808;
t1720 = -pkin(6) * t1653 - g(3) * t1693;
t1719 = t1691 * t1558;
t1718 = t1691 * t1744;
t1717 = t1695 * t1558;
t1716 = t1695 * t1744;
t1715 = pkin(4) * t1412 - t1290;
t1713 = pkin(4) * t1252 + t1729;
t1479 = t1545 * t1696 - t1547 * t1692;
t1711 = t1659 * t1697 - t1660 * t1693;
t1710 = pkin(3) * t1087 + t1750;
t1709 = pkin(3) * t1352 + t1749;
t1705 = pkin(3) * t1194 + t1726;
t1704 = pkin(3) * t1196 + t1727;
t1232 = t1317 * t1696 - t1794;
t1309 = -pkin(3) * t1511 + qJ(4) * t1317;
t1703 = -pkin(2) * t1511 + pkin(7) * t1232 - qJ(4) * t1794 + t1696 * t1309;
t1702 = pkin(3) * t1327 + t1715;
t1701 = pkin(3) * t1180 + t1713;
t1700 = pkin(3) * t1505 - t1407;
t1675 = t1686 * qJDD(3);
t1664 = t1679 - t1698;
t1662 = -t1678 + t1698;
t1655 = t1679 - t1678;
t1652 = qJDD(1) * t1693 + t1697 * t1808;
t1641 = t1696 * t1657;
t1640 = t1761 * t1801;
t1629 = -pkin(6) * t1652 + g(3) * t1697;
t1622 = -t1630 + t1698;
t1621 = -t1698 + t1807;
t1615 = t1646 * t1696 - t1682 * t1801;
t1614 = -t1692 * t1706 - t1801 * t1806;
t1613 = t1640 * t1688 + t1675;
t1612 = t1640 * t1686 - t1765;
t1610 = -t1662 * t1692 + t1641;
t1608 = t1664 * t1696 - t1771;
t1607 = t1663 * t1696 - t1771;
t1606 = t1662 * t1696 + t1772;
t1605 = t1665 * t1692 + t1641;
t1604 = t1664 * t1692 + t1770;
t1603 = (t1646 + t1672) * t1692;
t1599 = t1651 * t1688 - t1654 * t1686;
t1598 = t1651 * t1686 + t1654 * t1688;
t1595 = -t1645 * t1692 + t1602;
t1594 = t1645 * t1696 + t1647 * t1692;
t1588 = t1630 - t1807;
t1581 = t1615 * t1688 - t1738;
t1580 = t1614 * t1688 + t1738;
t1579 = t1615 * t1686 + t1737;
t1578 = t1614 * t1686 - t1737;
t1577 = t1610 * t1688 + t1686 * t1764;
t1576 = t1608 * t1688 + t1686 * t1763;
t1575 = t1610 * t1686 - t1688 * t1764;
t1574 = t1608 * t1686 - t1688 * t1763;
t1560 = (t1776 + t1777) * qJD(3);
t1559 = (-t1775 + t1778) * qJD(3);
t1556 = t1611 * t1688 + t1645 * t1686;
t1555 = t1609 * t1688 - t1647 * t1686;
t1554 = t1611 * t1686 - t1645 * t1688;
t1553 = t1609 * t1686 + t1647 * t1688;
t1552 = -pkin(1) * t1649 - t1592;
t1551 = pkin(1) * t1650 - t1591;
t1550 = -t1567 + t1762;
t1549 = t1566 - t1762;
t1546 = t1595 * t1688 - t1655 * t1686;
t1544 = t1595 * t1686 + t1655 * t1688;
t1527 = -qJD(3) * t1776 + t1593 * t1687;
t1526 = qJD(3) * t1775 + t1593 * t1685;
t1525 = -qJD(3) * t1777 + t1685 * t1740;
t1524 = -qJD(3) * t1778 - t1687 * t1740;
t1522 = -t1622 * t1685 + t1825;
t1521 = t1621 * t1687 - t1780;
t1519 = t1622 * t1687 + t1826;
t1518 = t1621 * t1685 + t1779;
t1517 = -pkin(7) * t1607 + t1562;
t1516 = -pkin(7) * t1605 + t1561;
t1510 = pkin(1) * t1514;
t1509 = -pkin(2) * t1607 + t1547;
t1508 = -pkin(2) * t1605 + t1545;
t1502 = pkin(1) * t1683 + qJ(2) * t1742;
t1501 = t1567 - t1566;
t1490 = -t1534 + t1564;
t1489 = t1533 - t1564;
t1488 = -t1559 * t1692 + t1560 * t1696;
t1487 = t1559 * t1696 + t1560 * t1692;
t1486 = -t1717 + t1718;
t1485 = -t1719 - t1716;
t1483 = t1534 - t1533;
t1482 = t1488 * t1688 + t1675;
t1481 = t1488 * t1686 - t1765;
t1475 = -t1539 * t1687 - t1542 * t1685;
t1473 = -t1539 * t1685 + t1542 * t1687;
t1468 = -t1526 * t1692 + t1527 * t1696;
t1467 = -t1524 * t1692 + t1525 * t1696;
t1466 = t1526 * t1696 + t1527 * t1692;
t1465 = t1524 * t1696 + t1525 * t1692;
t1463 = pkin(1) * t1553 + t1747;
t1462 = pkin(1) * t1554 + t1748;
t1460 = -t1519 * t1692 + t1522 * t1696;
t1459 = -t1518 * t1692 + t1521 * t1696;
t1458 = t1520 * t1696 + t1523 * t1692;
t1457 = t1519 * t1696 + t1522 * t1692;
t1456 = t1518 * t1696 + t1521 * t1692;
t1453 = t1549 * t1695 - t1788;
t1452 = -t1550 * t1691 + t1821;
t1451 = t1549 * t1691 + t1787;
t1450 = t1550 * t1695 + t1823;
t1447 = -qJ(2) * t1598 + t1479 * t1688;
t1446 = qJ(2) * t1599 + t1479 * t1686;
t1441 = t1480 * t1688 + t1572 * t1686;
t1440 = t1480 * t1686 - t1572 * t1688;
t1428 = pkin(1) * t1598 + t1722;
t1427 = t1468 * t1688 - t1755;
t1426 = t1467 * t1688 + t1755;
t1425 = t1468 * t1686 + t1754;
t1424 = t1467 * t1686 - t1754;
t1422 = t1695 * t1472 - t1718;
t1421 = t1691 * t1472 + t1716;
t1420 = -t1691 * t1471 + t1717;
t1419 = t1695 * t1471 + t1719;
t1415 = t1505 * t1696 + t1506 * t1692;
t1406 = (-t1535 * t1694 + t1537 * t1690) * t1565;
t1405 = (-t1535 * t1690 - t1537 * t1694) * t1565;
t1404 = t1461 * t1688 + t1542 * t1686;
t1403 = t1460 * t1688 + t1543 * t1686;
t1402 = t1459 * t1688 + t1541 * t1686;
t1401 = t1461 * t1686 - t1542 * t1688;
t1400 = t1460 * t1686 - t1543 * t1688;
t1399 = t1459 * t1686 - t1541 * t1688;
t1396 = -t1485 * t1685 + t1486 * t1687;
t1395 = t1485 * t1687 + t1486 * t1685;
t1393 = -qJ(2) * t1554 - t1509 * t1686 + t1517 * t1688;
t1392 = -qJ(2) * t1553 - t1508 * t1686 + t1516 * t1688;
t1385 = -t1473 * t1692 + t1475 * t1696;
t1384 = t1474 * t1696 + t1476 * t1692;
t1383 = t1473 * t1696 + t1475 * t1692;
t1382 = t1416 * t1688 + t1539 * t1686;
t1381 = t1416 * t1686 - t1539 * t1688;
t1379 = -pkin(1) * t1607 + qJ(2) * t1556 + t1509 * t1688 + t1517 * t1686;
t1378 = -pkin(1) * t1605 + qJ(2) * t1555 + t1508 * t1688 + t1516 * t1686;
t1376 = t1385 * t1688 + t1588 * t1686;
t1375 = t1385 * t1686 - t1588 * t1688;
t1371 = t1411 - t1491;
t1368 = -t1492 + t1707;
t1365 = t1411 * t1694 - t1537 * t1782;
t1364 = t1411 * t1690 + t1537 * t1781;
t1363 = t1535 * t1781 + t1690 * t1707;
t1362 = -t1535 * t1782 + t1694 * t1707;
t1359 = -t1451 * t1685 + t1453 * t1687;
t1358 = -t1450 * t1685 + t1452 * t1687;
t1357 = t1451 * t1687 + t1453 * t1685;
t1356 = t1450 * t1687 + t1452 * t1685;
t1355 = t1386 * t1688 + t1538 * t1686;
t1354 = t1386 * t1686 - t1538 * t1688;
t1351 = pkin(1) * t1440 + t1766;
t1349 = t1406 * t1695 + t1470 * t1691;
t1348 = t1406 * t1691 - t1470 * t1695;
t1347 = t1489 * t1694 - t1792;
t1346 = -t1490 * t1690 + t1822;
t1345 = t1489 * t1690 + t1791;
t1344 = t1490 * t1694 + t1824;
t1342 = -t1432 * t1695 - t1691 * t1820;
t1340 = -t1432 * t1691 + t1695 * t1820;
t1338 = -t1421 * t1685 + t1422 * t1687;
t1337 = -t1419 * t1685 + t1420 * t1687;
t1336 = t1421 * t1687 + t1422 * t1685;
t1335 = t1419 * t1687 + t1420 * t1685;
t1334 = -pkin(2) * t1384 - t1804;
t1322 = -qJ(2) * t1440 - (pkin(2) * t1686 - pkin(7) * t1688) * t1479;
t1321 = t1365 * t1695 + t1759;
t1320 = t1363 * t1695 - t1759;
t1319 = t1365 * t1691 - t1758;
t1318 = t1363 * t1691 + t1758;
t1315 = -pkin(2) * t1458 + t1625 - t1721;
t1314 = -t1395 * t1692 + t1396 * t1696;
t1313 = t1395 * t1696 + t1396 * t1692;
t1312 = t1314 * t1688 + t1681 * t1686;
t1311 = t1314 * t1686 - t1681 * t1688;
t1310 = -pkin(2) * t1415 - t1700;
t1306 = -pkin(7) * t1458 - t1397 * t1692 + t1444 * t1696;
t1304 = qJ(2) * t1441 - (-pkin(2) * t1688 - pkin(7) * t1686 - pkin(1)) * t1479;
t1298 = t1369 * t1694 - t1371 * t1690;
t1296 = t1369 * t1690 + t1371 * t1694;
t1295 = -pkin(7) * t1415 - t1391 * t1692 + t1429 * t1696;
t1287 = -t1357 * t1692 + t1359 * t1696;
t1286 = -t1356 * t1692 + t1358 * t1696;
t1285 = t1357 * t1696 + t1359 * t1692;
t1284 = t1356 * t1696 + t1358 * t1692;
t1282 = t1352 * t1696 + t1353 * t1692;
t1280 = t1347 * t1695 - t1368 * t1691;
t1279 = t1346 * t1695 + t1372 * t1691;
t1278 = t1347 * t1691 + t1368 * t1695;
t1277 = t1346 * t1691 - t1372 * t1695;
t1276 = -t1348 * t1685 + t1349 * t1687;
t1275 = t1348 * t1687 + t1349 * t1685;
t1267 = -t1340 * t1685 + t1342 * t1687;
t1265 = t1340 * t1687 + t1342 * t1685;
t1264 = t1298 * t1695 + t1483 * t1691;
t1263 = t1298 * t1691 - t1483 * t1695;
t1262 = -t1336 * t1692 + t1338 * t1696;
t1261 = -t1335 * t1692 + t1337 * t1696;
t1260 = t1336 * t1696 + t1338 * t1692;
t1259 = t1335 * t1696 + t1337 * t1692;
t1250 = pkin(1) * t1401 + t1723;
t1248 = t1327 * t1696 + t1328 * t1692;
t1246 = t1287 * t1688 + t1434 * t1686;
t1245 = t1286 * t1688 + t1437 * t1686;
t1244 = t1287 * t1686 - t1434 * t1688;
t1243 = t1286 * t1686 - t1437 * t1688;
t1242 = t1283 * t1688 + t1686 * t1820;
t1241 = t1283 * t1686 - t1688 * t1820;
t1240 = t1262 * t1688 + t1757;
t1239 = t1261 * t1688 - t1757;
t1238 = t1262 * t1686 - t1756;
t1237 = t1261 * t1686 + t1756;
t1236 = -t1319 * t1685 + t1321 * t1687;
t1235 = -t1318 * t1685 + t1320 * t1687;
t1234 = t1319 * t1687 + t1321 * t1685;
t1233 = t1318 * t1687 + t1320 * t1685;
t1231 = t1317 * t1692 + t1793;
t1229 = pkin(1) * t1381 + t1724;
t1228 = t1232 * t1688 + t1511 * t1686;
t1227 = t1232 * t1686 - t1511 * t1688;
t1226 = t1249 * t1688 + t1432 * t1686;
t1225 = t1249 * t1686 - t1432 * t1688;
t1222 = -qJ(2) * t1401 + t1306 * t1688 - t1315 * t1686;
t1216 = -pkin(2) * t1231 - t1805;
t1211 = -t1278 * t1685 + t1280 * t1687;
t1210 = -t1277 * t1685 + t1279 * t1687;
t1209 = t1278 * t1687 + t1280 * t1685;
t1208 = t1277 * t1687 + t1279 * t1685;
t1207 = -pkin(1) * t1458 + qJ(2) * t1404 + t1306 * t1686 + t1315 * t1688;
t1206 = -qJ(2) * t1381 + t1295 * t1688 - t1310 * t1686;
t1200 = -t1275 * t1692 + t1276 * t1696;
t1199 = t1275 * t1696 + t1276 * t1692;
t1198 = -pkin(7) * t1384 - t1292 * t1692 + t1302 * t1696;
t1192 = -t1265 * t1692 + t1267 * t1696;
t1191 = t1266 * t1696 + t1268 * t1692;
t1190 = t1265 * t1696 + t1267 * t1692;
t1188 = -t1263 * t1685 + t1264 * t1687;
t1187 = t1263 * t1687 + t1264 * t1685;
t1184 = -pkin(1) * t1415 + qJ(2) * t1382 + t1295 * t1686 + t1310 * t1688;
t1183 = t1192 * t1688 + t1501 * t1686;
t1182 = t1192 * t1686 - t1501 * t1688;
t1178 = t1193 * t1688 + t1469 * t1686;
t1177 = t1193 * t1686 - t1469 * t1688;
t1176 = t1200 * t1688 + t1405 * t1686;
t1175 = t1200 * t1686 - t1405 * t1688;
t1171 = -t1234 * t1692 + t1236 * t1696;
t1170 = -t1233 * t1692 + t1235 * t1696;
t1169 = t1234 * t1696 + t1236 * t1692;
t1168 = t1233 * t1696 + t1235 * t1692;
t1167 = -pkin(2) * t1282 - t1709;
t1166 = pkin(1) * t1354 + t1725;
t1165 = -pkin(7) * t1231 - qJ(4) * t1793 - t1309 * t1692;
t1164 = -pkin(2) * t1248 - t1702;
t1163 = t1171 * t1688 + t1364 * t1686;
t1162 = t1170 * t1688 - t1362 * t1686;
t1161 = t1171 * t1686 - t1364 * t1688;
t1160 = t1170 * t1686 + t1362 * t1688;
t1159 = -qJ(2) * t1354 + t1198 * t1688 - t1334 * t1686;
t1158 = -pkin(1) * t1384 + qJ(2) * t1355 + t1198 * t1686 + t1334 * t1688;
t1155 = -pkin(2) * t1191 - t1745;
t1150 = -t1209 * t1692 + t1211 * t1696;
t1149 = -t1208 * t1692 + t1210 * t1696;
t1148 = t1209 * t1696 + t1211 * t1692;
t1147 = t1208 * t1696 + t1210 * t1692;
t1145 = t1196 * t1696 + t1197 * t1692;
t1143 = pkin(1) * t1227 + t1703;
t1141 = t1194 * t1696 + t1195 * t1692;
t1139 = -t1187 * t1692 + t1188 * t1696;
t1138 = t1187 * t1696 + t1188 * t1692;
t1137 = t1150 * t1688 + t1345 * t1686;
t1136 = t1149 * t1688 + t1344 * t1686;
t1135 = t1150 * t1686 - t1345 * t1688;
t1134 = t1149 * t1686 - t1344 * t1688;
t1132 = t1180 * t1696 + t1181 * t1692;
t1129 = -pkin(7) * t1282 - t1203 * t1692 + t1221 * t1696;
t1128 = t1146 * t1688 + t1331 * t1686;
t1127 = t1146 * t1686 - t1331 * t1688;
t1126 = t1142 * t1688 + t1325 * t1686;
t1125 = t1142 * t1686 - t1325 * t1688;
t1121 = t1139 * t1688 + t1296 * t1686;
t1120 = t1139 * t1686 - t1296 * t1688;
t1117 = -pkin(7) * t1248 - t1186 * t1692 + t1205 * t1696;
t1116 = t1133 * t1688 + t1297 * t1686;
t1115 = t1133 * t1686 - t1297 * t1688;
t1110 = -qJ(2) * t1227 + t1165 * t1688 - t1216 * t1686;
t1107 = pkin(1) * t1241 + t1728;
t1105 = pkin(1) * t1225 + t1730;
t1103 = t1156 * t1696 + t1157 * t1692;
t1101 = -pkin(1) * t1231 + qJ(2) * t1228 + t1165 * t1686 + t1216 * t1688;
t1098 = t1104 * t1688 + t1418 * t1686;
t1097 = t1104 * t1686 - t1418 * t1688;
t1093 = -qJ(2) * t1241 + t1129 * t1688 - t1167 * t1686;
t1092 = -pkin(2) * t1145 - t1704;
t1091 = -pkin(1) * t1282 + qJ(2) * t1242 + t1129 * t1686 + t1167 * t1688;
t1090 = -pkin(2) * t1141 - t1705;
t1089 = -qJ(2) * t1225 + t1117 * t1688 - t1164 * t1686;
t1086 = -pkin(1) * t1248 + qJ(2) * t1226 + t1117 * t1686 + t1164 * t1688;
t1084 = -pkin(2) * t1103 - t1746;
t1083 = -pkin(2) * t1132 - t1701;
t1074 = -pkin(7) * t1191 - t1109 * t1692 + t1112 * t1696;
t1072 = pkin(1) * t1177 + t1731;
t1066 = t1087 * t1696 + t1088 * t1692;
t1064 = -qJ(2) * t1177 + t1074 * t1688 - t1155 * t1686;
t1063 = -pkin(1) * t1191 + qJ(2) * t1178 + t1074 * t1686 + t1155 * t1688;
t1062 = -pkin(7) * t1103 - t1095 * t1692 + t1100 * t1696;
t1061 = t1067 * t1688 + t1153 * t1686;
t1060 = t1067 * t1686 - t1153 * t1688;
t1059 = pkin(1) * t1097 + t1732;
t1058 = -pkin(7) * t1145 - t1078 * t1692 + t1082 * t1696;
t1057 = -pkin(7) * t1141 - t1076 * t1692 + t1080 * t1696;
t1056 = pkin(1) * t1127 + t1733;
t1055 = pkin(1) * t1125 + t1734;
t1052 = -pkin(7) * t1132 - t1069 * t1692 + t1071 * t1696;
t1051 = -pkin(2) * t1066 - t1710;
t1050 = -qJ(2) * t1097 + t1062 * t1688 - t1084 * t1686;
t1047 = -qJ(2) * t1127 + t1058 * t1688 - t1092 * t1686;
t1046 = pkin(1) * t1115 + t1735;
t1045 = -qJ(2) * t1125 + t1057 * t1688 - t1090 * t1686;
t1044 = -pkin(1) * t1145 + qJ(2) * t1128 + t1058 * t1686 + t1092 * t1688;
t1043 = -pkin(1) * t1103 + qJ(2) * t1098 + t1062 * t1686 + t1084 * t1688;
t1042 = -pkin(1) * t1141 + qJ(2) * t1126 + t1057 * t1686 + t1090 * t1688;
t1041 = -qJ(2) * t1115 + t1052 * t1688 - t1083 * t1686;
t1040 = -pkin(1) * t1132 + qJ(2) * t1116 + t1052 * t1686 + t1083 * t1688;
t1039 = -pkin(7) * t1066 - t1049 * t1692 + t1054 * t1696;
t1038 = pkin(1) * t1060 + t1736;
t1037 = -qJ(2) * t1060 + t1039 * t1688 - t1051 * t1686;
t1036 = -pkin(1) * t1066 + qJ(2) * t1061 + t1039 * t1686 + t1051 * t1688;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1653, 0, -t1652, 0, t1720, -t1629, -t1711, -pkin(6) * t1711, 0, 0, t1596, 0, -t1811, 0, t1818, t1819, t1827, pkin(6) * t1827 + qJ(2) * t1783 - t1693 * t1502, -t1579 * t1693 + t1581 * t1697, -t1544 * t1693 + t1546 * t1697, -t1575 * t1693 + t1577 * t1697, -t1578 * t1693 + t1580 * t1697, -t1574 * t1693 + t1576 * t1697, -t1612 * t1693 + t1613 * t1697, t1697 * t1392 - t1693 * t1378 - pkin(6) * (t1553 * t1697 + t1555 * t1693), t1697 * t1393 - t1693 * t1379 - pkin(6) * (t1554 * t1697 + t1556 * t1693), t1697 * t1447 - t1693 * t1446 - pkin(6) * (t1598 * t1697 + t1599 * t1693), t1697 * t1322 - t1693 * t1304 - pkin(6) * (t1440 * t1697 + t1441 * t1693), -t1425 * t1693 + t1427 * t1697, -t1375 * t1693 + t1376 * t1697, -t1400 * t1693 + t1403 * t1697, -t1424 * t1693 + t1426 * t1697, -t1399 * t1693 + t1402 * t1697, -t1481 * t1693 + t1482 * t1697, t1697 * t1206 - t1693 * t1184 - pkin(6) * (t1381 * t1697 + t1382 * t1693), t1697 * t1222 - t1693 * t1207 - pkin(6) * (t1401 * t1697 + t1404 * t1693), t1697 * t1159 - t1693 * t1158 - pkin(6) * (t1354 * t1697 + t1355 * t1693), t1697 * t1110 - t1693 * t1101 - pkin(6) * (t1227 * t1697 + t1228 * t1693), -t1238 * t1693 + t1240 * t1697, -t1182 * t1693 + t1183 * t1697, -t1243 * t1693 + t1245 * t1697, -t1237 * t1693 + t1239 * t1697, -t1244 * t1693 + t1246 * t1697, -t1311 * t1693 + t1312 * t1697, t1697 * t1089 - t1693 * t1086 - pkin(6) * (t1225 * t1697 + t1226 * t1693), t1697 * t1093 - t1693 * t1091 - pkin(6) * (t1241 * t1697 + t1242 * t1693), t1697 * t1064 - t1693 * t1063 - pkin(6) * (t1177 * t1697 + t1178 * t1693), t1697 * t1050 - t1693 * t1043 - pkin(6) * (t1097 * t1697 + t1098 * t1693), -t1161 * t1693 + t1163 * t1697, -t1120 * t1693 + t1121 * t1697, -t1134 * t1693 + t1136 * t1697, -t1160 * t1693 + t1162 * t1697, -t1135 * t1693 + t1137 * t1697, -t1175 * t1693 + t1176 * t1697, t1697 * t1045 - t1693 * t1042 - pkin(6) * (t1125 * t1697 + t1126 * t1693), t1697 * t1047 - t1693 * t1044 - pkin(6) * (t1127 * t1697 + t1128 * t1693), t1697 * t1041 - t1693 * t1040 - pkin(6) * (t1115 * t1697 + t1116 * t1693), t1697 * t1037 - t1693 * t1036 - pkin(6) * (t1060 * t1697 + t1061 * t1693); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1652, 0, t1653, 0, t1629, t1720, t1739, pkin(6) * t1739, 0, 0, t1811, 0, t1596, 0, -t1819, t1818, t1828, pkin(6) * t1828 + qJ(2) * t1784 + t1697 * t1502, t1579 * t1697 + t1581 * t1693, t1544 * t1697 + t1546 * t1693, t1575 * t1697 + t1577 * t1693, t1578 * t1697 + t1580 * t1693, t1574 * t1697 + t1576 * t1693, t1612 * t1697 + t1613 * t1693, t1693 * t1392 + t1697 * t1378 + pkin(6) * (-t1553 * t1693 + t1555 * t1697), t1693 * t1393 + t1697 * t1379 + pkin(6) * (-t1554 * t1693 + t1556 * t1697), t1693 * t1447 + t1697 * t1446 + pkin(6) * (-t1598 * t1693 + t1599 * t1697), t1693 * t1322 + t1697 * t1304 + pkin(6) * (-t1440 * t1693 + t1441 * t1697), t1425 * t1697 + t1427 * t1693, t1375 * t1697 + t1376 * t1693, t1400 * t1697 + t1403 * t1693, t1424 * t1697 + t1426 * t1693, t1399 * t1697 + t1402 * t1693, t1481 * t1697 + t1482 * t1693, t1693 * t1206 + t1697 * t1184 + pkin(6) * (-t1381 * t1693 + t1382 * t1697), t1693 * t1222 + t1697 * t1207 + pkin(6) * (-t1401 * t1693 + t1404 * t1697), t1693 * t1159 + t1697 * t1158 + pkin(6) * (-t1354 * t1693 + t1355 * t1697), t1693 * t1110 + t1697 * t1101 + pkin(6) * (-t1227 * t1693 + t1228 * t1697), t1238 * t1697 + t1240 * t1693, t1182 * t1697 + t1183 * t1693, t1243 * t1697 + t1245 * t1693, t1237 * t1697 + t1239 * t1693, t1244 * t1697 + t1246 * t1693, t1311 * t1697 + t1312 * t1693, t1693 * t1089 + t1697 * t1086 + pkin(6) * (-t1225 * t1693 + t1226 * t1697), t1693 * t1093 + t1697 * t1091 + pkin(6) * (-t1241 * t1693 + t1242 * t1697), t1693 * t1064 + t1697 * t1063 + pkin(6) * (-t1177 * t1693 + t1178 * t1697), t1693 * t1050 + t1697 * t1043 + pkin(6) * (-t1097 * t1693 + t1098 * t1697), t1161 * t1697 + t1163 * t1693, t1120 * t1697 + t1121 * t1693, t1134 * t1697 + t1136 * t1693, t1160 * t1697 + t1162 * t1693, t1135 * t1697 + t1137 * t1693, t1175 * t1697 + t1176 * t1693, t1693 * t1045 + t1697 * t1042 + pkin(6) * (-t1125 * t1693 + t1126 * t1697), t1693 * t1047 + t1697 * t1044 + pkin(6) * (-t1127 * t1693 + t1128 * t1697), t1693 * t1041 + t1697 * t1040 + pkin(6) * (-t1115 * t1693 + t1116 * t1697), t1693 * t1037 + t1697 * t1036 + pkin(6) * (-t1060 * t1693 + t1061 * t1697); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1659, t1660, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1551, t1552, 0, -t1510, t1603, t1594, t1606, t1602, t1604, 0, t1463, t1462, t1428, t1351, t1466, t1383, t1457, t1465, t1456, t1487, t1229, t1250, t1166, t1143, t1260, t1190, t1284, t1259, t1285, t1313, t1105, t1107, t1072, t1059, t1169, t1138, t1147, t1168, t1148, t1199, t1055, t1056, t1046, t1038; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1808, 0, 0, -g(3), -t1659, 0, 0, 0, t1650, 0, -t1649, 0, t1714, t1620, t1514, qJ(2) * t1514, t1581, t1546, t1577, t1580, t1576, t1613, t1392, t1393, t1447, t1322, t1427, t1376, t1403, t1426, t1402, t1482, t1206, t1222, t1159, t1110, t1240, t1183, t1245, t1239, t1246, t1312, t1089, t1093, t1064, t1050, t1163, t1121, t1136, t1162, t1137, t1176, t1045, t1047, t1041, t1037; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1808, 0, qJDD(1), 0, g(3), 0, -t1660, 0, 0, 0, t1649, 0, t1650, 0, -t1620, t1714, t1742, t1502, t1579, t1544, t1575, t1578, t1574, t1612, t1378, t1379, t1446, t1304, t1425, t1375, t1400, t1424, t1399, t1481, t1184, t1207, t1158, t1101, t1238, t1182, t1243, t1237, t1244, t1311, t1086, t1091, t1063, t1043, t1161, t1120, t1134, t1160, t1135, t1175, t1042, t1044, t1040, t1036; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1659, t1660, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1551, t1552, 0, -t1510, t1603, t1594, t1606, t1602, t1604, 0, t1463, t1462, t1428, t1351, t1466, t1383, t1457, t1465, t1456, t1487, t1229, t1250, t1166, t1143, t1260, t1190, t1284, t1259, t1285, t1313, t1105, t1107, t1072, t1059, t1169, t1138, t1147, t1168, t1148, t1199, t1055, t1056, t1046, t1038; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1808, 0, 0, -t1683, t1591, 0, t1615, t1595, t1610, t1614, t1608, t1640, t1516, t1517, t1479, pkin(7) * t1479, t1468, t1385, t1460, t1467, t1459, t1488, t1295, t1306, t1198, t1165, t1262, t1192, t1286, t1261, t1287, t1314, t1117, t1129, t1074, t1062, t1171, t1139, t1149, t1170, t1150, t1200, t1057, t1058, t1052, t1039; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1808, 0, qJDD(1), 0, t1683, 0, t1592, 0, t1666, t1655, -t1764, -t1666, -t1763, -qJDD(3), t1508, t1509, 0, pkin(2) * t1479, t1590, -t1588, -t1543, -t1590, -t1541, -qJDD(3), t1310, t1315, t1334, t1216, -t1504, -t1501, -t1437, t1504, -t1434, -t1681, t1164, t1167, t1155, t1084, -t1364, -t1296, -t1344, t1362, -t1345, -t1405, t1090, t1092, t1083, t1051; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1591, -t1592, 0, 0, t1603, t1594, t1606, t1602, t1604, 0, t1747, t1748, t1722, t1766, t1466, t1383, t1457, t1465, t1456, t1487, t1724, t1723, t1725, t1703, t1260, t1190, t1284, t1259, t1285, t1313, t1730, t1728, t1731, t1732, t1169, t1138, t1147, t1168, t1148, t1199, t1734, t1733, t1735, t1736; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1646, t1647, t1657, -t1672, t1664, t1672, 0, t1572, t1545, 0, t1527, t1475, t1522, t1525, t1521, t1560, t1429, t1444, t1302, -qJ(4) * t1316, t1338, t1267, t1358, t1337, t1359, t1396, t1205, t1221, t1112, t1100, t1236, t1188, t1210, t1235, t1211, t1276, t1080, t1082, t1071, t1054; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1753, t1645, t1662, t1706, t1658, -t1753, -t1572, 0, t1547, 0, t1526, t1473, t1519, t1524, t1518, t1559, t1391, t1397, t1292, t1309, t1336, t1265, t1356, t1335, t1357, t1395, t1186, t1203, t1109, t1095, t1234, t1187, t1208, t1233, t1209, t1275, t1076, t1078, t1069, t1049; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1666, -t1655, t1764, t1666, t1763, qJDD(3), -t1545, -t1547, 0, 0, -t1590, t1588, t1543, t1590, t1541, qJDD(3), t1700, t1721 - 0.2e1 * t1797, t1804, t1805, t1504, t1501, t1437, -t1504, t1434, t1681, t1702, t1709, t1745, t1746, t1364, t1296, t1344, -t1362, t1345, t1405, t1705, t1704, t1701, t1710; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1593, -t1539, t1810, -t1627, t1621, t1627, 0, t1511, t1407, 0, t1422, t1342, t1452, t1420, t1453, t1486, t1333, t1350, t1179, -pkin(8) * t1218, t1321, t1264, t1279, t1320, t1280, t1349, t1118, t1119, t1106, t1085; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1798, t1542, t1622, -t1740, t1585, -t1798, -t1511, 0, t1408, 0, t1421, t1340, t1450, t1419, t1451, t1485, t1300, t1307, t1172, t1201, t1319, t1263, t1277, t1318, t1278, t1348, t1113, t1114, t1096, t1073; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1590, t1588, t1543, t1590, t1541, qJDD(3), -t1407, -t1408, 0, 0, t1504, t1501, t1437, -t1504, t1434, t1681, t1715, t1749, t1339, t1217, t1364, t1296, t1344, -t1362, t1345, t1405, t1726, t1727, t1713, t1750; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1472, -t1432, t1812, t1558, t1549, -t1558, 0, t1418, t1290, 0, t1365, t1298, t1346, t1363, t1347, t1406, t1223, t1224, t1130, -pkin(9) * t1153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1744, t1820, t1550, t1471, t1496, -t1744, -t1418, 0, t1291, 0, -t1484, -t1483, -t1372, t1484, t1368, -t1470, t1173, t1174, -pkin(5) * t1297, -pkin(5) * t1153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1504, t1501, t1437, -t1504, t1434, t1681, -t1290, -t1291, 0, 0, t1364, t1296, t1344, -t1362, t1345, t1405, t1751, t1752, t1729, t1768; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1411, t1369, t1813, t1491, t1489, -t1491, 0, t1257, t1214, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1492, t1371, t1490, -t1707, t1388, -t1492, -t1257, 0, t1215, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1484, t1483, t1372, -t1484, -t1368, t1470, -t1214, -t1215, 0, 0;];
m_new_reg  = t1;
