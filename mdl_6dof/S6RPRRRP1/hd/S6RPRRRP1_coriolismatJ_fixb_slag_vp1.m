% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:34
% EndTime: 2019-03-09 05:56:46
% DurationCPUTime: 63.88s
% Computational Cost: add. (239809->1254), mult. (205754->1771), div. (0->0), fcn. (224903->10), ass. (0->815)
t999 = qJ(1) + pkin(10);
t995 = cos(t999);
t1510 = -0.2e1 * t995;
t1001 = sin(qJ(5));
t1004 = cos(qJ(5));
t1000 = qJ(3) + qJ(4);
t996 = sin(t1000);
t927 = (-rSges(6,1) * t1001 - rSges(6,2) * t1004) * t996;
t1193 = t927 * t1510;
t994 = sin(t999);
t1360 = t927 * t994;
t1502 = -pkin(8) - pkin(7);
t1005 = cos(qJ(3));
t1423 = pkin(3) * t1005;
t990 = pkin(2) + t1423;
t1229 = -t995 * t1502 - t994 * t990;
t1426 = sin(qJ(1)) * pkin(1);
t1107 = t1229 - t1426;
t997 = cos(t1000);
t1442 = pkin(4) * t997;
t1477 = rSges(6,3) + pkin(9);
t1560 = t1477 * t996;
t1200 = t997 * t1001;
t899 = t1004 * t995 + t1200 * t994;
t1201 = t995 * t1001;
t1202 = t994 * t1004;
t900 = t1202 * t997 - t1201;
t1565 = -t900 * rSges(6,1) + t899 * rSges(6,2);
t579 = (-t1442 - t1560) * t994 + t1107 + t1565;
t901 = t1200 * t995 - t1202;
t1211 = t1004 * t997;
t1217 = t1001 * t994;
t902 = t1211 * t995 + t1217;
t1098 = t902 * rSges(6,1) - t901 * rSges(6,2);
t1425 = cos(qJ(1)) * pkin(1);
t978 = t994 * t1502;
t1130 = -t978 + t1425;
t1137 = t990 + t1442;
t580 = (t1137 + t1560) * t995 + t1098 + t1130;
t1287 = t579 * t1193 - 0.2e1 * t580 * t1360;
t1479 = rSges(7,1) + pkin(5);
t1617 = rSges(7,3) + qJ(6);
t1233 = (-t1001 * t1479 + t1004 * t1617) * t996;
t721 = t1233 * t995;
t1519 = -0.2e1 * t721;
t1422 = t901 * rSges(7,3);
t1478 = rSges(7,2) + pkin(9);
t1559 = t1478 * t996;
t873 = t901 * qJ(6);
t528 = t1422 + t873 + t1479 * t902 + (t1137 + t1559) * t995 + t1130;
t1539 = 0.2e1 * t528;
t1019 = (-t1442 - t1559) * t994 + t1107;
t527 = -t1479 * t900 - t1617 * t899 + t1019;
t720 = t1233 * t994;
t1295 = t527 * t1519 - t720 * t1539;
t1503 = m(7) / 0.4e1;
t1505 = m(6) / 0.4e1;
t743 = -rSges(6,1) * t901 - rSges(6,2) * t902;
t1517 = -0.2e1 * t743;
t738 = -rSges(6,1) * t899 - rSges(6,2) * t900;
t1518 = -0.2e1 * t738;
t1092 = pkin(5) * t1004 + qJ(6) * t1001;
t1095 = rSges(7,1) * t1004 + rSges(7,3) * t1001;
t1649 = t1092 + t1095;
t1239 = -rSges(7,2) * t997 + t1649 * t996;
t945 = pkin(4) * t996 - pkin(9) * t997;
t1174 = t945 + t1239;
t640 = t1174 * t995;
t1528 = -0.2e1 * t640;
t638 = t1174 * t994;
t1529 = -0.2e1 * t638;
t603 = t1479 * t899 - t1617 * t900;
t741 = -t901 * pkin(5) + qJ(6) * t902;
t742 = -t901 * rSges(7,1) + t902 * rSges(7,3);
t604 = t741 + t742;
t1414 = rSges(6,1) * t1004;
t1096 = -rSges(6,2) * t1001 + t1414;
t1427 = rSges(6,3) * t997;
t870 = t1096 * t996 - t1427;
t1238 = t870 + t945;
t716 = t1238 * t994;
t718 = t1238 * t995;
t1567 = (t1528 * t603 + t1529 * t604 + t1295) * t1503 + (t1517 * t716 - t1518 * t718 + t1287) * t1505;
t1534 = 0.2e1 * t580;
t872 = rSges(6,3) * t996 + t1096 * t997;
t1081 = t870 * t997 + t872 * t996;
t1343 = t994 * t996;
t701 = rSges(6,3) * t1343 - t1565;
t1216 = t1001 * t996;
t1186 = rSges(6,2) * t1216;
t1168 = t996 * t1202;
t935 = rSges(6,1) * t1168;
t773 = -t935 + (t1186 + t1427) * t994;
t513 = t1081 * t994 - t701 * t996 + t773 * t997;
t1323 = t995 * t996;
t705 = rSges(6,3) * t1323 + t1098;
t1212 = t1004 * t996;
t1169 = t996 * t1201;
t1322 = t995 * t997;
t1230 = rSges(6,2) * t1169 + rSges(6,3) * t1322;
t775 = -rSges(6,1) * t1212 * t995 + t1230;
t514 = -t1081 * t995 + t705 * t996 - t775 * t997;
t1299 = t514 * t1534 + 0.2e1 * t513 * t579;
t1170 = t994 * t1216;
t1232 = rSges(7,1) * t1168 + rSges(7,3) * t1170;
t1342 = t994 * t997;
t772 = rSges(7,2) * t1342 - t1232;
t842 = -pkin(5) * t1168 - qJ(6) * t1170;
t1263 = t772 + t842;
t1237 = rSges(7,2) * t996 + t1649 * t997;
t1552 = t1237 * t996 + t1239 * t997;
t1234 = -t900 * rSges(7,1) - t899 * rSges(7,3);
t700 = rSges(7,2) * t1343 - t1234;
t740 = -t900 * pkin(5) - qJ(6) * t899;
t1637 = t700 - t740;
t379 = t1263 * t997 + t1552 * t994 - t1637 * t996;
t954 = rSges(7,2) * t1322;
t774 = -t1095 * t1323 + t954;
t843 = t1092 * t1323;
t1262 = -t774 + t843;
t704 = t902 * rSges(7,1) + rSges(7,2) * t1323 + t1422;
t744 = t902 * pkin(5) + t873;
t1267 = t704 + t744;
t380 = t1262 * t997 + t1267 * t996 - t1552 * t995;
t1308 = t380 * t1539 + 0.2e1 * t379 * t527;
t1314 = t997 * t705;
t591 = t1323 * t870 + t1314;
t1532 = -0.2e1 * t591;
t1315 = t997 * t701;
t589 = t870 * t1343 + t1315;
t1533 = 0.2e1 * t589;
t1561 = t1267 * t997;
t512 = t1239 * t1323 + t1561;
t1540 = -0.2e1 * t512;
t1126 = t1637 * t997;
t1659 = t1239 * t1343 + t1126;
t1541 = 0.2e1 * t1659;
t1027 = (-t1001 * t1617 - t1004 * t1479 - pkin(4)) * t996;
t956 = pkin(9) * t1322;
t1228 = t954 + t956;
t605 = t1027 * t995 + t1228;
t955 = pkin(4) * t1343;
t1097 = t955 - t842 + t1232;
t1171 = t1478 * t997;
t606 = -t1171 * t994 + t1097;
t1042 = -t1477 * t997 - t1186;
t1231 = t935 + t955;
t676 = t1042 * t994 + t1231;
t1100 = (-pkin(4) - t1414) * t996;
t1172 = t956 + t1230;
t677 = t1100 * t995 + t1172;
t1573 = (t1532 * t677 + t1533 * t676 + t1299) * t1505 + (t1540 * t605 + t1541 * t606 + t1308) * t1503;
t1678 = t1573 - t1567;
t1002 = sin(qJ(3));
t1424 = pkin(3) * t1002;
t1129 = t945 + t1424;
t1063 = t1129 + t1239;
t613 = t1063 * t995;
t1530 = -0.2e1 * t613;
t611 = t1063 * t994;
t1531 = -0.2e1 * t611;
t1106 = t1129 + t870;
t696 = t1106 * t994;
t698 = t1106 * t995;
t1568 = (t1530 * t603 + t1531 * t604 + t1295) * t1503 + (t1517 * t696 - t1518 * t698 + t1287) * t1505;
t593 = (t1027 - t1424) * t995 + t1228;
t594 = (-t1171 + t1424) * t994 + t1097;
t648 = (t1042 + t1424) * t994 + t1231;
t649 = (t1100 - t1424) * t995 + t1172;
t1574 = (t1532 * t649 + t1533 * t648 + t1299) * t1505 + (t1540 * t593 + t1541 * t594 + t1308) * t1503;
t1677 = t1574 - t1568;
t1481 = -t997 / 0.2e1;
t876 = Icges(7,5) * t902;
t680 = Icges(7,6) * t1323 + Icges(7,3) * t901 + t876;
t1408 = Icges(7,5) * t901;
t692 = Icges(7,1) * t902 + Icges(7,4) * t1323 + t1408;
t1069 = t1001 * t680 + t1004 * t692;
t686 = Icges(7,4) * t902 + Icges(7,2) * t1323 + Icges(7,6) * t901;
t1375 = t686 * t997;
t498 = t1069 * t996 - t1375;
t1411 = Icges(6,4) * t902;
t689 = -Icges(6,2) * t901 + Icges(6,6) * t1323 + t1411;
t879 = Icges(6,4) * t901;
t695 = Icges(6,1) * t902 + Icges(6,5) * t1323 - t879;
t1067 = -t1001 * t689 + t1004 * t695;
t683 = Icges(6,5) * t902 - Icges(6,6) * t901 + Icges(6,3) * t1323;
t1377 = t683 * t997;
t501 = t1067 * t996 - t1377;
t1596 = t498 + t501;
t1404 = Icges(7,6) * t997;
t975 = Icges(7,5) * t1212;
t857 = Icges(7,3) * t1216 - t1404 + t975;
t1220 = Icges(7,5) * t1001;
t1077 = Icges(7,1) * t1004 + t1220;
t865 = -Icges(7,4) * t997 + t1077 * t996;
t1066 = t1001 * t857 + t1004 * t865;
t1075 = Icges(7,4) * t1004 + Icges(7,6) * t1001;
t862 = Icges(7,2) * t996 + t1075 * t997;
t1045 = t1066 - t862;
t861 = -Icges(7,2) * t997 + t1075 * t996;
t1368 = t861 * t997;
t1021 = -t1045 * t996 + t1368;
t1071 = Icges(7,5) * t1004 + Icges(7,3) * t1001;
t858 = Icges(7,6) * t996 + t1071 * t997;
t866 = Icges(7,4) * t996 + t1077 * t997;
t435 = t1021 * t995 + t858 * t901 + t866 * t902;
t1221 = Icges(6,4) * t1004;
t1076 = -Icges(6,2) * t1001 + t1221;
t863 = -Icges(6,6) * t997 + t1076 * t996;
t1222 = Icges(6,4) * t1001;
t1078 = Icges(6,1) * t1004 - t1222;
t867 = -Icges(6,5) * t997 + t1078 * t996;
t1065 = -t1001 * t863 + t1004 * t867;
t1072 = Icges(6,5) * t1004 - Icges(6,6) * t1001;
t860 = Icges(6,3) * t996 + t1072 * t997;
t1044 = -t1065 + t860;
t859 = -Icges(6,3) * t997 + t1072 * t996;
t1369 = t859 * t997;
t1022 = t1044 * t996 + t1369;
t864 = Icges(6,6) * t996 + t1076 * t997;
t868 = Icges(6,5) * t996 + t1078 * t997;
t436 = t1022 * t995 - t864 * t901 + t868 * t902;
t1598 = t435 + t436;
t1623 = t996 / 0.2e1;
t875 = Icges(7,5) * t900;
t678 = Icges(7,6) * t1343 + Icges(7,3) * t899 + t875;
t874 = Icges(7,5) * t899;
t691 = -Icges(7,1) * t900 - Icges(7,4) * t1343 - t874;
t1070 = t1001 * t678 - t1004 * t691;
t684 = Icges(7,4) * t900 + Icges(7,2) * t1343 + Icges(7,6) * t899;
t1376 = t684 * t997;
t496 = t1070 * t996 - t1376;
t878 = Icges(6,4) * t900;
t688 = Icges(6,2) * t899 - Icges(6,6) * t1343 - t878;
t877 = Icges(6,4) * t899;
t693 = Icges(6,1) * t900 + Icges(6,5) * t1343 - t877;
t1068 = t1001 * t688 + t1004 * t693;
t681 = Icges(6,5) * t900 - Icges(6,6) * t899 + Icges(6,3) * t1343;
t1378 = t681 * t997;
t499 = t1068 * t996 - t1378;
t1663 = -t496 - t499;
t1048 = t861 * t995 + t1069;
t1029 = -t1071 * t996 + t1404;
t761 = t1029 * t995;
t769 = t865 * t995;
t375 = t1048 * t997 + (t1001 * t761 - t1004 * t769 + t686) * t996;
t1046 = -t859 * t995 - t1067;
t767 = t863 * t995;
t771 = t867 * t995;
t377 = -t1046 * t997 + (t1001 * t767 - t1004 * t771 + t683) * t996;
t503 = t1045 * t997 + (t1001 * t858 + t1004 * t866 + t861) * t996;
t504 = -t1044 * t997 + (-t1001 * t864 + t1004 * t868 + t859) * t996;
t1320 = t996 * t861;
t554 = t1320 * t994 + t857 * t899 + t865 * t900;
t1321 = t996 * t859;
t555 = t1321 * t994 - t863 * t899 + t867 * t900;
t558 = t1320 * t995 + t857 * t901 + t865 * t902;
t559 = t1321 * t995 - t863 * t901 + t867 * t902;
t595 = t1066 * t996 - t1368;
t596 = t1065 * t996 - t1369;
t1676 = (t595 + t596) * t1623 + (t503 + t504) * t1481 + (t554 + t555 - t1663) * t1342 / 0.4e1 + (t558 + t559 + t1596) * t1322 / 0.4e1 + t1323 * (t375 + t377 + t1598) / 0.4e1;
t726 = -Icges(6,5) * t899 - Icges(6,6) * t900;
t728 = -Icges(7,4) * t899 + Icges(7,6) * t900;
t1675 = t726 + t728;
t727 = -Icges(6,5) * t901 - Icges(6,6) * t902;
t729 = -Icges(7,4) * t901 + Icges(7,6) * t902;
t1674 = t727 + t729;
t1271 = Icges(6,2) * t902 - t695 + t879;
t1273 = Icges(7,3) * t902 - t1408 - t692;
t1673 = t1271 + t1273;
t1272 = Icges(6,2) * t900 - t693 + t877;
t1274 = Icges(7,3) * t900 + t691 - t874;
t1672 = t1272 + t1274;
t1275 = -Icges(6,1) * t901 - t1411 - t689;
t1277 = -Icges(7,1) * t901 + t680 + t876;
t1671 = t1275 + t1277;
t1276 = -Icges(6,1) * t899 + t688 - t878;
t1278 = -Icges(7,1) * t899 + t678 + t875;
t1670 = t1276 + t1278;
t1643 = -m(6) * qJD(1) / 0.4e1;
t444 = t1343 * t681 + t688 * t899 + t693 * t900;
t445 = t683 * t1343 - t899 * t689 + t900 * t695;
t1090 = t444 * t994 + t445 * t995;
t80 = -t1090 * t996 + t555 * t997;
t442 = t1343 * t684 + t678 * t899 - t691 * t900;
t443 = t686 * t1343 + t899 * t680 + t900 * t692;
t1091 = t442 * t994 + t443 * t995;
t79 = -t1091 * t996 + t554 * t997;
t446 = t684 * t1323 + t901 * t678 - t691 * t902;
t447 = t686 * t1323 + t901 * t680 + t902 * t692;
t1089 = t446 * t994 + t447 * t995;
t1657 = t1089 * t996 - t558 * t997;
t448 = t681 * t1323 + t688 * t901 + t902 * t693;
t449 = t683 * t1323 - t901 * t689 + t902 * t695;
t1088 = t448 * t994 + t449 * t995;
t1658 = t1088 * t996 - t559 * t997;
t1638 = t1657 + t1658;
t1668 = t1343 * t1675 + t1670 * t900 + t1672 * t899;
t1667 = t1343 * t1674 + t1671 * t900 + t1673 * t899;
t1666 = t1323 * t1675 + t1670 * t902 + t1672 * t901;
t1665 = t1323 * t1674 + t1671 * t902 + t1673 * t901;
t919 = (Icges(7,3) * t1004 - t1220) * t996;
t920 = (-Icges(6,5) * t1001 - Icges(6,6) * t1004) * t996;
t921 = (-Icges(7,4) * t1001 + Icges(7,6) * t1004) * t996;
t922 = (-Icges(6,2) * t1004 - t1222) * t996;
t923 = -Icges(7,1) * t1216 + t975;
t924 = (-Icges(6,1) * t1001 - t1221) * t996;
t1656 = -(-(t867 / 0.2e1 + t922 / 0.2e1 + t865 / 0.2e1 - t919 / 0.2e1) * t1001 + (t924 / 0.2e1 - t863 / 0.2e1 + t923 / 0.2e1 + t857 / 0.2e1) * t1004) * t996 + (t920 / 0.2e1 + t921 / 0.2e1) * t997;
t1164 = t1343 / 0.4e1;
t1166 = -t1343 / 0.4e1;
t1645 = (t447 + t449) * t994 + (-t448 - t446) * t995;
t1655 = (t1164 + t1166) * t1645;
t1241 = -t865 + t919;
t1243 = t857 + t923;
t480 = t1241 * t899 + t1243 * t900 + t1343 * t921;
t1240 = -t867 - t922;
t1242 = -t863 + t924;
t481 = t1240 * t899 + t1242 * t900 + t1343 * t920;
t1654 = -t481 - t480;
t482 = t1241 * t901 + t1243 * t902 + t1323 * t921;
t483 = t1240 * t901 + t1242 * t902 + t1323 * t920;
t1653 = t482 + t483;
t1196 = qJD(3) + qJD(4);
t1371 = t775 * t996;
t1084 = -t1314 - t1371;
t1372 = t773 * t996;
t1085 = t1315 + t1372;
t736 = -pkin(5) * t899 + qJ(6) * t900;
t737 = -rSges(7,1) * t899 + rSges(7,3) * t900;
t1266 = t736 + t737;
t1109 = 0.2e1 * t1266;
t1125 = t604 * t995;
t1197 = 0.2e1 * m(7);
t1133 = t1197 / 0.4e1;
t1631 = t1262 * t996 - t1561;
t1633 = t1263 * t996 + t1126;
t1644 = m(6) / 0.2e1;
t575 = t738 * t994 + t743 * t995;
t138 = t575 * t1644 + (t1109 * t994 + 0.2e1 * t1125) * t1503 + (t1085 * t1644 + t1133 * t1633) * t995 + (t1084 * t1644 + t1133 * t1631) * t994;
t1250 = t994 * (pkin(9) * t1342 - t955) + t995 * (-pkin(4) * t1323 + t956);
t1251 = 0.2e1 * t1250;
t749 = t994 * t772;
t751 = t995 * t774;
t809 = t994 * t842;
t810 = t995 * t843;
t1056 = -0.2e1 * t810 + 0.2e1 * t749 + 0.2e1 * t751 + 0.2e1 * t809 + t1251;
t750 = t994 * t773;
t752 = t995 * t775;
t1113 = 0.2e1 * t750 + 0.2e1 * t752 + t1251;
t943 = rSges(5,1) * t996 + rSges(5,2) * t997;
t894 = t943 * t994;
t895 = t943 * t995;
t707 = -t994 * t894 - t995 * t895;
t1252 = 0.2e1 * t707;
t1504 = m(7) / 0.2e1;
t1508 = m(5) / 0.2e1;
t321 = t1056 * t1504 + t1113 * t1644 + t1252 * t1508;
t1139 = t1200 / 0.2e1;
t1324 = t995 * t901;
t1344 = t994 * t899;
t708 = t1324 + t1344;
t661 = m(7) * t1139 + t1133 * t708;
t1648 = t321 * qJD(4) + t138 * qJD(5) + t661 * qJD(6);
t137 = ((-t741 / 0.2e1 - t742 / 0.2e1 + (t700 / 0.2e1 - t740 / 0.2e1) * t997 + (t772 / 0.2e1 + t842 / 0.2e1) * t996) * t995 + (-t736 / 0.2e1 - t737 / 0.2e1 + (-t704 / 0.2e1 - t744 / 0.2e1) * t997 + (-t774 / 0.2e1 + t843 / 0.2e1) * t996) * t994) * m(7) + ((t1315 / 0.2e1 + t1372 / 0.2e1 - t743 / 0.2e1) * t995 + (-t1314 / 0.2e1 - t1371 / 0.2e1 - t738 / 0.2e1) * t994) * m(6);
t660 = (t1139 - t1344 / 0.2e1 - t1324 / 0.2e1) * m(7);
t1646 = -t137 * qJD(5) - t660 * qJD(6);
t1433 = m(7) * qJD(1);
t1641 = t1667 * t994 - t1668 * t995;
t1640 = t1665 * t994 - t1666 * t995;
t989 = Icges(5,4) * t997;
t939 = -Icges(5,2) * t996 + t989;
t940 = Icges(5,1) * t996 + t989;
t1636 = t940 + t939;
t991 = t994 ^ 2;
t992 = t995 ^ 2;
t1225 = t991 + t992;
t1108 = t1225 * t1424;
t1507 = m(5) / 0.4e1;
t1486 = t994 / 0.2e1;
t1485 = -t995 / 0.2e1;
t1626 = -t995 / 0.4e1;
t1484 = t995 / 0.2e1;
t1624 = -t996 / 0.2e1;
t1622 = t997 / 0.2e1;
t1049 = t861 * t994 + t1070;
t1024 = -t1049 * t996 + t1376;
t760 = t1029 * t994;
t768 = t865 * t994;
t338 = t1024 * t994 + t760 * t899 - t768 * t900;
t1023 = -t1048 * t996 + t1375;
t339 = t1023 * t994 + t761 * t899 - t769 * t900;
t433 = t1021 * t994 + t858 * t899 + t866 * t900;
t65 = (t1091 - t433) * t997 + (t338 * t994 + t339 * t995 + t554) * t996;
t1047 = -t859 * t994 - t1068;
t1026 = t1047 * t996 + t1378;
t766 = t863 * t994;
t770 = t867 * t994;
t340 = t1026 * t994 + t766 * t899 - t770 * t900;
t1025 = t1046 * t996 + t1377;
t341 = t1025 * t994 + t767 * t899 - t771 * t900;
t434 = t1022 * t994 - t864 * t899 + t868 * t900;
t66 = (t1090 - t434) * t997 + (t340 * t994 + t341 * t995 + t555) * t996;
t1619 = t66 + t65;
t342 = t1024 * t995 + t760 * t901 - t768 * t902;
t343 = t1023 * t995 + t761 * t901 - t769 * t902;
t67 = (t1089 - t435) * t997 + (t342 * t994 + t343 * t995 + t558) * t996;
t344 = t1026 * t995 + t766 * t901 - t770 * t902;
t345 = t1025 * t995 + t767 * t901 - t771 * t902;
t68 = (t1088 - t436) * t997 + (t344 * t994 + t345 * t995 + t559) * t996;
t1618 = t68 + t67;
t1429 = m(7) * qJD(6);
t1605 = t1654 * t997 + (t1667 * t995 + t1668 * t994) * t996;
t1604 = -t1653 * t997 + (t1665 * t995 + t1666 * t994) * t996;
t1603 = (-t338 - t340) * t995 + (t339 + t341) * t994;
t1602 = (-t342 - t344) * t995 + (t343 + t345) * t994;
t1599 = t433 + t434;
t1595 = t868 + t866;
t1428 = rSges(5,1) * t997;
t944 = -rSges(5,2) * t996 + t1428;
t1358 = t944 * t994;
t1131 = -rSges(5,2) * t1323 + t994 * rSges(5,3);
t710 = (t990 + t1428) * t995 + t1130 + t1131;
t642 = -0.2e1 * t710 * t1358;
t1192 = t944 * t1510;
t830 = rSges(5,1) * t1342 - rSges(5,2) * t1343 - t995 * rSges(5,3);
t709 = t1107 - t830;
t643 = t709 * t1192;
t1279 = t642 + t643;
t946 = pkin(9) * t996 + t1442;
t1236 = -t872 - t946;
t717 = t1236 * t994;
t523 = t717 * t1534;
t719 = t1236 * t995;
t1520 = 0.2e1 * t719;
t524 = t579 * t1520;
t1289 = t523 + t524;
t1173 = -t946 - t1237;
t639 = t1173 * t994;
t453 = t639 * t1539;
t641 = t1173 * t995;
t1527 = 0.2e1 * t641;
t454 = t527 * t1527;
t1298 = t453 + t454;
t1521 = -0.2e1 * t718;
t1522 = -0.2e1 * t716;
t1041 = t943 + t1424;
t1563 = t1041 * t995;
t1564 = t1041 * t994;
t1182 = (t1528 * t594 + t1529 * t593 + t1298) * t1503 + (t1521 * t648 + t1522 * t649 + t1289) * t1505 + (0.2e1 * (t1563 * t994 - t1564 * t995) * t943 + t1279) * t1507;
t1524 = -0.2e1 * t698;
t1525 = -0.2e1 * t696;
t1183 = (t1530 * t606 + t1531 * t605 + t1298) * t1503 + (t1524 * t676 + t1525 * t677 + t1289) * t1505 + (-0.2e1 * t1563 * t894 + 0.2e1 * t1564 * t895 + t1279) * t1507;
t1585 = t1182 - t1183;
t417 = (-t1267 * t994 + t1637 * t995) * t996;
t1544 = 0.2e1 * t417;
t1555 = t749 + t751 + t809 - t810;
t461 = t1250 + t1555;
t1181 = t1527 * t1659 + t639 * t1540 + t461 * t1544;
t320 = t1631 * t994 + t1633 * t995;
t1546 = 0.2e1 * t320;
t1249 = t1225 * t946;
t397 = t1267 * t995 + t1637 * t994 + t1249;
t35 = t1528 * t379 + t1529 * t380 + t1546 * t397 + t1181;
t551 = (t995 * t701 - t705 * t994) * t996;
t1537 = 0.2e1 * t551;
t1566 = t750 + t752;
t534 = t1250 + t1566;
t1178 = t589 * t1520 + t717 * t1532 + t534 * t1537;
t440 = t1084 * t994 + t1085 * t995;
t1543 = 0.2e1 * t440;
t508 = t994 * t701 + t995 * t705 + t1249;
t90 = t1521 * t513 + t1522 * t514 + t1543 * t508 + t1178;
t1439 = t1503 * t35 + t1505 * t90;
t988 = t995 * pkin(7);
t1260 = -t994 * (pkin(2) * t994 + t1229 - t988) + t995 * (-t994 * pkin(7) - t978 + (-pkin(2) + t990) * t995);
t350 = t397 + t1260;
t174 = t350 * t1546;
t403 = t508 + t1260;
t273 = t403 * t1543;
t334 = t379 * t1530;
t335 = t380 * t1531;
t438 = t513 * t1524;
t439 = t514 * t1525;
t1440 = (t334 + t335 + t174 + t1181) * t1503 + (t438 + t439 + t273 + t1178) * t1505;
t1584 = t1439 - t1440;
t374 = t1049 * t997 + (t1001 * t760 - t1004 * t768 + t684) * t996;
t376 = -t1047 * t997 + (t1001 * t766 - t1004 * t770 + t681) * t996;
t1582 = t374 + t376 + t1599;
t1412 = Icges(5,4) * t996;
t938 = Icges(5,2) * t997 + t1412;
t941 = Icges(5,1) * t997 - t1412;
t1579 = t941 * t1623 + t938 * t1624 + t1320 / 0.2e1 + t1321 / 0.2e1 + t1636 * t1622 + (t862 + t860) * t1481 + (t867 + t865) * t1211 / 0.2e1;
t998 = Icges(4,4) * t1005;
t961 = -Icges(4,2) * t1002 + t998;
t962 = Icges(4,1) * t1002 + t998;
t646 = t994 * t830 + t995 * (rSges(5,1) * t1322 + t1131);
t1558 = t1225 * t944 * t943 + t646 * t707;
t1557 = m(6) * (t440 * t551 + t513 * t589 - t514 * t591) + m(7) * (t1659 * t379 + t320 * t417 - t380 * t512);
t1112 = 0.2e1 * t350 + 0.2e1 * t397;
t470 = t1266 * t994 + t1125;
t1556 = ((t403 + t508) * t575 + ((t698 + t718) * t995 + (t696 + t716) * t994) * t927) * t1644 + (t470 * t1112 + (-t613 - t640) * t1519 - 0.2e1 * (-t611 - t638) * t720) * t1503;
t1102 = t403 * t534 - t696 * t717 - t698 * t719;
t1103 = t350 * t461 - t611 * t639 - t613 * t641;
t1509 = m(4) / 0.4e1;
t1415 = rSges(4,1) * t1005;
t1138 = pkin(2) + t1415;
t1214 = t1002 * t994;
t1226 = rSges(4,2) * t1214 + t995 * rSges(4,3);
t753 = -t1138 * t994 + t1226 - t1426 + t988;
t1213 = t1002 * t995;
t974 = rSges(4,2) * t1213;
t754 = t1425 - t974 + t1138 * t995 + (rSges(4,3) + pkin(7)) * t994;
t966 = rSges(4,1) * t1002 + rSges(4,2) * t1005;
t917 = t966 * t994;
t918 = t966 * t995;
t1223 = Icges(4,4) * t1002;
t960 = Icges(4,2) * t1005 + t1223;
t963 = Icges(4,1) * t1005 - t1223;
t1554 = -(t963 / 0.2e1 - t960 / 0.2e1) * t1002 - 0.4e1 * (-t1563 * t710 + t1564 * t709) * t1507 - 0.4e1 * (t753 * t917 - t754 * t918) * t1509;
t1115 = t1485 * t1641 + t1486 * t1640;
t855 = Icges(4,5) * t994 + t963 * t995;
t1244 = -t960 * t995 + t855;
t1210 = t1005 * t994;
t972 = Icges(4,4) * t1214;
t854 = Icges(4,1) * t1210 - Icges(4,5) * t995 - t972;
t1245 = -Icges(4,2) * t1210 + t854 - t972;
t853 = Icges(4,6) * t994 + t961 * t995;
t1246 = -t962 * t995 - t853;
t852 = Icges(4,4) * t1210 - Icges(4,2) * t1214 - Icges(4,6) * t995;
t1247 = t962 * t994 + t852;
t1549 = (-t1244 * t994 + t1245 * t995) * t1002 + (t1246 * t994 + t1247 * t995) * t1005;
t395 = -t727 * t997 + (t1001 * t1271 + t1004 * t1275) * t996;
t1393 = t395 * t994;
t394 = -t726 * t997 + (t1001 * t1272 + t1004 * t1276) * t996;
t1394 = t394 * t995;
t393 = -t729 * t997 + (t1001 * t1273 + t1004 * t1277) * t996;
t1395 = t393 * t994;
t392 = -t728 * t997 + (t1001 * t1274 + t1004 * t1278) * t996;
t1396 = t392 * t995;
t1038 = -t1396 / 0.4e1 + t1395 / 0.4e1 - t1394 / 0.4e1 + t1393 / 0.4e1 + t1653 * t994 / 0.4e1 - t1654 * t1626;
t1547 = (t1626 + t995 / 0.4e1) * t1638;
t824 = Icges(5,5) * t1342 - Icges(5,6) * t1343 - Icges(5,3) * t995;
t827 = Icges(5,6) * t994 + t939 * t995;
t1120 = t827 * t996 - t824;
t829 = Icges(5,5) * t994 + t941 * t995;
t757 = t829 * t1342;
t937 = Icges(5,5) * t997 - Icges(5,6) * t996;
t1359 = t937 * t995;
t825 = Icges(5,3) * t994 + t1359;
t1122 = t995 * t825 - t757;
t1264 = t829 * t1322 + t994 * t825;
t949 = Icges(5,4) * t1343;
t828 = Icges(5,1) * t1342 - Icges(5,5) * t995 - t949;
t1265 = -t828 * t1322 - t994 * t824;
t826 = Icges(5,4) * t1342 - Icges(5,2) * t1343 - Icges(5,6) * t995;
t1370 = t826 * t996;
t568 = -t1343 * t827 - t1122;
t569 = -t1323 * t826 - t1265;
t570 = -t1323 * t827 + t1264;
t1020 = ((t568 - t757 + (t825 + t1370) * t995 + t1265) * t995 + t1264 * t994 + t1645) * t1485 + (-t569 * t995 + t570 * t994 + t1645) * t1484 + ((t1120 * t994 + t1122 + t568 + t569) * t994 + (-t1264 + t570 + (-t828 * t997 + t1370) * t994 + (t1120 + t824) * t995) * t995) * t1486;
t450 = (t1266 * t995 - t604 * t994) * t996;
t1542 = 0.2e1 * t450;
t571 = (t738 * t995 - t743 * t994) * t996;
t1536 = 0.2e1 * t571;
t1535 = 0.4e1 * t575;
t1361 = t901 * t994;
t1365 = t899 * t995;
t665 = (-t1361 + t1365) * t996;
t1526 = 0.2e1 * t665;
t1523 = 0.4e1 * t708;
t993 = t996 ^ 2;
t755 = t1217 * t993 + t899 * t997;
t1516 = 0.2e1 * t755;
t756 = -t1201 * t993 - t901 * t997;
t1515 = -0.2e1 * t756;
t1514 = 0.2e1 * t899;
t1513 = -0.2e1 * t900;
t1512 = 0.2e1 * t901;
t1511 = 0.2e1 * t902;
t1117 = -0.2e1 * t1169;
t1118 = 0.2e1 * t1170;
t1292 = t1117 * t1659 + t512 * t1118;
t1493 = m(7) * (t379 * t1512 + t380 * t1514 + 0.2e1 * (t320 * t996 + t417 * t997) * t1001 + t1292);
t191 = 0.4e1 * t1659 * t755 + 0.4e1 * t417 * t665 - 0.4e1 * t512 * t756;
t1488 = t191 / 0.4e1;
t1487 = -t994 / 0.2e1;
t1473 = 0.4e1 * m(5) * (t709 * t894 - t710 * t895);
t1179 = t708 * t1544 + t1292;
t1462 = m(7) * (t1515 * t611 - t1516 * t613 + t1526 * t350 + t1179);
t1461 = m(7) * (t1515 * t638 - t1516 * t640 + t1526 * t397 + t1179);
t1281 = -t640 * t1117 + t638 * t1118;
t601 = t611 * t1118;
t602 = t613 * t1117;
t1282 = t601 - t602;
t1460 = m(7) * (t1112 * t708 + t1281 + t1282);
t1364 = t901 * t512;
t1296 = -0.2e1 * t1659 * t899 - 0.2e1 * t1364;
t1459 = m(7) * (t1514 * t1659 + t1296 + 0.2e1 * t1364);
t1191 = 0.2e1 * t1216;
t1176 = t461 * t1191 + t641 * t1512 + t639 * t1514;
t1189 = 0.2e1 * t1200;
t346 = t350 * t1189;
t1458 = m(7) * (t346 + t1176 + t1282);
t1175 = t470 * t1191 - t721 * t1512 - t720 * t1514;
t1190 = 0.2e1 * t1212;
t1457 = m(7) * (t1190 * t350 - t1511 * t613 + t1513 * t611 + t1175);
t161 = t1189 * t397 + t1176 + t1281;
t1456 = m(7) * t161;
t1455 = m(7) * (t1190 * t397 - t1511 * t640 + t1513 * t638 + t1175);
t1363 = t901 * t611;
t1367 = t899 * t613;
t1286 = -0.2e1 * t1363 + 0.2e1 * t1367;
t1454 = m(7) * (t1286 + 0.2e1 * t1363 - 0.2e1 * t1367);
t1451 = m(7) * (t1516 * t527 + t1539 * t756 + t1296);
t1362 = t901 * t638;
t1366 = t899 * t640;
t1283 = -0.2e1 * t1362 + 0.2e1 * t1366;
t1448 = m(7) * (t1283 + 0.2e1 * t1362 - 0.2e1 * t1366);
t1447 = (t527 * t902 + t528 * t900 + t603 * t901 + t604 * t899) * t1197;
t1290 = t527 * t1117 - 0.2e1 * t528 * t1170;
t1446 = m(7) * (t1512 * t594 + t1514 * t593 + t1290);
t1445 = m(7) * (t1512 * t606 + t1514 * t605 + t1290);
t1444 = m(7) * (t1286 + t1290);
t1443 = m(7) * (t1283 + t1290);
t1438 = m(5) * qJD(3);
t1437 = m(5) * qJD(4);
t1435 = m(6) * qJD(3);
t1434 = m(6) * qJD(4);
t1432 = m(7) * qJD(3);
t1431 = m(7) * qJD(4);
t1430 = m(7) * qJD(5);
t1391 = t496 * t994;
t1087 = t498 * t995 + t1391;
t112 = (t1087 - t503) * t997 + (t374 * t994 + t375 * t995 + t595) * t996;
t1390 = t499 * t994;
t1086 = t501 * t995 + t1390;
t113 = (t1086 - t504) * t997 + (t376 * t994 + t377 * t995 + t596) * t996;
t269 = t1087 * t996 - t595 * t997;
t270 = t1086 * t996 - t596 * t997;
t12 = (-t113 / 0.2e1 - t112 / 0.2e1 + (t1658 / 0.2e1 + t1657 / 0.2e1) * t995 + (-t79 / 0.2e1 - t80 / 0.2e1) * t994) * t997 + (t270 / 0.2e1 + t269 / 0.2e1 + (t68 / 0.2e1 + t67 / 0.2e1) * t995 + (t65 / 0.2e1 + t66 / 0.2e1) * t994) * t996 + t1557;
t1208 = t137 * qJD(2);
t1416 = t12 * qJD(5) + t1208;
t1400 = t374 * t995;
t1399 = t375 * t994;
t1398 = t376 * t995;
t1397 = t377 * t994;
t532 = t646 + t1260;
t1387 = t532 * t707;
t1203 = t660 * qJD(2);
t1119 = 0.4e1 * t1216;
t617 = (-t708 + t1200) * t1119;
t1297 = -t617 * t1429 / 0.4e1 - t1203;
t1288 = t527 - t1019 - t740 - t1234;
t1185 = t1429 / 0.4e1;
t1280 = t617 * t1185 + t1203;
t1209 = t1005 * t995;
t850 = Icges(4,5) * t1210 - Icges(4,6) * t1214 - Icges(4,3) * t995;
t1258 = -t854 * t1209 - t994 * t850;
t1074 = Icges(4,5) * t1005 - Icges(4,6) * t1002;
t851 = Icges(4,3) * t994 + t1074 * t995;
t1257 = t855 * t1209 + t994 * t851;
t1256 = t940 * t994 + t826;
t1255 = -t940 * t995 - t827;
t1254 = -Icges(5,2) * t1342 + t828 - t949;
t1253 = -t938 * t995 + t829;
t1227 = -0.2e1 * t1108;
t1218 = qJD(5) * t996;
t1215 = t1002 * t852;
t632 = (t1004 / 0.2e1 - t1365 / 0.2e1 + t1361 / 0.2e1) * t996 * m(7);
t1204 = t632 * qJD(2);
t1194 = 0.4e1 * t927;
t1180 = t1519 * t1659 - t720 * t1540 + t470 * t1544;
t1177 = t589 * t1193 + 0.2e1 * t591 * t1360 + t575 * t1537;
t1165 = t1343 / 0.2e1;
t1157 = t1323 / 0.2e1;
t1146 = -t1216 / 0.2e1;
t1145 = t1216 / 0.2e1;
t1143 = t1212 / 0.2e1;
t1140 = -t1200 / 0.2e1;
t1128 = -t944 - t1423;
t1127 = -t946 - t1423;
t806 = t855 * t1210;
t1121 = t995 * t851 - t806;
t1116 = t1002 * t853 - t850;
t1105 = t1127 - t872;
t1101 = t1387 + t1564 * t1358 - t1563 * t1192 / 0.2e1;
t1093 = Icges(5,5) * t996 + Icges(5,6) * t997;
t1083 = t826 * t997 + t828 * t996;
t1080 = -t917 * t994 - t918 * t995;
t1073 = -Icges(4,5) * t1002 - Icges(4,6) * t1005;
t1062 = t1127 - t1237;
t1060 = t1115 + t1556;
t1035 = -t1253 * t996 + t1255 * t997;
t1036 = t1254 * t996 + t1256 * t997;
t880 = t1093 * t994;
t881 = t995 * t1093;
t1057 = (-t991 * t881 + (t1036 * t995 + (t880 + t1035) * t994) * t995 + t1602) * t1486 + (-t992 * t880 + (t1035 * t994 + (t881 + t1036) * t995) * t994 + t1603) * t1485;
t1043 = -t1108 + t1250;
t1040 = t397 * t461 - t638 * t639 - t640 * t641;
t1039 = t508 * t534 - t716 * t717 - t718 * t719;
t1034 = (-t938 + t941) * t997 - t1636 * t996;
t1016 = t857 * t1139 + t863 * t1140 + t1143 * t1595 + t858 * t1145 + t864 * t1146 + t1579;
t1015 = t1547 + t1655;
t1014 = -t1115 + t1618 * t1486 + t1619 * t1485 + (t1596 * t994 + t1663 * t995) * t1623 + (t1399 - t1400 + t1397 - t1398) * t1481 + t1603 * t1165 + (-(-t445 - t443) * t994 - (t442 + t444) * t995) * t1342 / 0.2e1 + t1602 * t1157 + t1645 * t1322 / 0.2e1;
t1013 = t1164 * t1582 + t1676;
t1012 = -t1557 + t1604 * t1486 + t1605 * t1485 + (t269 + t270) * t1624 + (t1395 - t1396 + t1393 - t1394) * t1481 + (t113 + t112) * t1622 - t1619 * t1343 / 0.2e1 + t1641 * t1165 - (-t80 - t79) * t1342 / 0.2e1 - t1618 * t1323 / 0.2e1 + t1640 * t1157 - t1638 * t1322 / 0.2e1;
t1011 = t1397 / 0.2e1 - t1398 / 0.2e1 + t1399 / 0.2e1 - t1400 / 0.2e1 - t1020 + (t1034 * t995 + t1253 * t997 + t1255 * t996 + t937 * t994 + t1598) * t1486 + (t1034 * t994 + t1254 * t997 - t1256 * t996 - t1359 + t1599) * t1485;
t1010 = t858 * t1146 + t1083 * t1487 + t863 * t1139 + t857 * t1140 + t864 * t1145 + t1390 / 0.2e1 + t1391 / 0.2e1 - t1595 * t1212 / 0.2e1 + (t1083 + t1663) * t1486 - t1579;
t1009 = t1013 + t1038 - t1547 + t1655;
t1008 = t1166 * t1582 + t1015 + t1038 - t1676;
t1007 = t1015 + t1013 - t1038;
t968 = -rSges(4,2) * t1002 + t1415;
t912 = t1073 * t995;
t911 = t1073 * t994;
t838 = t1128 * t995;
t836 = t1128 * t994;
t699 = t1105 * t995;
t697 = t1105 * t994;
t634 = -t1108 + t707;
t633 = m(7) * t1143 + t1133 * t665;
t616 = -t1323 * t927 - t743 * t997;
t615 = t1343 * t927 + t738 * t997;
t614 = t1062 * t995;
t612 = t1062 * t994;
t592 = 0.4e1 * t1001 * t1004 * t993 + 0.4e1 * t899 * t900 + 0.4e1 * t901 * t902;
t588 = -t1213 * t853 + t1257;
t587 = -t1213 * t852 - t1258;
t586 = -t1214 * t853 - t1121;
t536 = -t997 * t920 + (t1001 * t1240 + t1004 * t1242) * t996;
t535 = -t997 * t921 + (t1001 * t1241 + t1004 * t1243) * t996;
t520 = -t604 * t997 - t721 * t996;
t519 = t1233 * t1343 + t1266 * t997;
t516 = t1043 + t1566;
t502 = 0.4e1 * t1558;
t460 = -t587 * t995 + t588 * t994;
t459 = -(-(-t1005 * t854 + t1215) * t994 - t995 * t850) * t995 + t586 * t994;
t441 = t1043 + t1555;
t391 = 0.4e1 * t1387 + 0.4e1 * (t1563 * t995 + t1564 * t994) * t944;
t390 = -0.4e1 * t579 * t738 + 0.4e1 * t580 * t743;
t378 = 0.4e1 * t579 * t676 + 0.4e1 * t580 * t677;
t369 = -0.4e1 * t527 * t899 + 0.4e1 * t528 * t901;
t360 = 0.4e1 * t579 * t648 + 0.4e1 * t580 * t649;
t317 = t508 * t1535 + (t716 * t994 + t718 * t995) * t1194;
t296 = 0.4e1 * t527 * t606 + 0.4e1 * t528 * t605;
t275 = 0.4e1 * t527 * t603 + 0.4e1 * t528 * t604;
t274 = 0.4e1 * t527 * t594 + 0.4e1 * t528 * t593;
t271 = t397 * t1523 + (t638 * t994 + t640 * t995) * t1119;
t268 = t1443 / 0.4e1;
t265 = 0.4e1 * t1039;
t263 = t403 * t1535 + (t696 * t994 + t698 * t995) * t1194;
t260 = t1444 / 0.4e1;
t256 = t1445 / 0.4e1;
t250 = t1446 / 0.4e1;
t247 = t1447 / 0.4e1;
t245 = t350 * t1523 + (t611 * t994 + t613 * t995) * t1119;
t238 = 0.4e1 * t1102;
t229 = (t586 - t806 + (t851 + t1215) * t995 + t1258) * t995 + t1257 * t994;
t228 = (t1116 * t995 - t1257 + t588) * t995 + (t1116 * t994 + t1121 + t587) * t994;
t221 = t1448 / 0.4e1;
t194 = 0.4e1 * t397 * t470 + 0.4e1 * t638 * t720 + 0.4e1 * t640 * t721;
t184 = t1451 / 0.4e1;
t169 = t1454 / 0.4e1;
t168 = 0.4e1 * t1040;
t167 = 0.4e1 * t350 * t470 + 0.4e1 * t611 * t720 + 0.4e1 * t613 * t721;
t163 = t1455 / 0.4e1;
t160 = t1456 / 0.4e1;
t157 = 0.4e1 * t1103;
t151 = t1457 / 0.4e1;
t139 = t1458 / 0.4e1;
t130 = t1459 / 0.4e1;
t129 = t275 * t1503 + t390 * t1505 - t1656;
t127 = t1460 / 0.4e1;
t122 = t1461 / 0.4e1;
t118 = t1462 / 0.4e1;
t91 = t1016 + t1473 / 0.4e1 + t378 * t1505 + t296 * t1503;
t56 = (t962 / 0.2e1 + t961 / 0.2e1) * t1005 + t1016 + t360 * t1505 + t274 * t1503 - t1554;
t53 = t256 + t221 - t1443 / 0.4e1;
t52 = t268 + t256 - t1448 / 0.4e1;
t51 = t268 + t221 - t1445 / 0.4e1;
t49 = -t1493 / 0.4e1;
t48 = t1493 / 0.4e1;
t47 = t250 + t169 - t1444 / 0.4e1;
t46 = t260 + t250 - t1454 / 0.4e1;
t45 = t260 + t169 - t1446 / 0.4e1;
t38 = t247 + t130 - t1451 / 0.4e1;
t37 = t184 + t247 - t1459 / 0.4e1;
t36 = t184 + t130 - t1447 / 0.4e1;
t31 = t139 + t160 - t1460 / 0.4e1;
t30 = t127 + t160 - t1458 / 0.4e1;
t29 = t127 + t139 - t1456 / 0.4e1;
t28 = t163 + t48 - t1461 / 0.4e1;
t27 = t122 + t163 + t49;
t26 = t122 + t48 - t1455 / 0.4e1;
t25 = t151 + t48 - t1462 / 0.4e1;
t24 = t118 + t151 + t49;
t23 = t118 + t48 - t1457 / 0.4e1;
t22 = t1503 * t194 + t1505 * t317 + t1115;
t21 = t1503 * t167 + t1505 * t263 + t1115;
t20 = t1503 * t168 + t1505 * t265 + t1507 * t502 + t1057;
t19 = t20 * qJD(4);
t18 = t1503 * t157 + t1505 * t238 + t1507 * t391 + t1057;
t15 = t1060 + t1584;
t14 = t1060 - t1584;
t13 = (t228 / 0.2e1 + t459 / 0.2e1) * t994 + (t460 / 0.2e1 - t229 / 0.2e1) * t995 + t1020;
t10 = t1020 - t1585;
t9 = t1020 + t1585;
t8 = t1011 + t1182 + t1183;
t7 = t1014 + t1439 + t1440 - t1556;
t6 = t1009 + t1567 + t1573;
t5 = t1007 + t1678;
t4 = t1008 - t1678;
t3 = t1007 + t1677;
t2 = t1008 - t1677;
t1 = t1009 + t1568 + t1574;
t11 = [-t1288 * t1433 * t528 + t56 * qJD(3) + t91 * qJD(4) + t129 * qJD(5) + t369 * t1185, 0, t56 * qJD(1) + t8 * qJD(4) + t1 * qJD(5) + t46 * qJD(6) + (t527 * t614 + t528 * t612 - t593 * t611 - t594 * t613) * t1432 + (t579 * t699 + t580 * t697 - t648 * t698 - t649 * t696) * t1435 + (t709 * t838 + t710 * t836) * t1438 + (t229 * t1484 + t1011 + (t1002 * t1246 + t1005 * t1244) * t1486 + (t228 + t459) * t1487 + (-t1002 * t1247 + t1005 * t1245 + t460) * t1485 + (t992 / 0.2e1 + t991 / 0.2e1) * t1074 + ((-t753 * t995 - t754 * t994) * t968 + (-t917 * t995 + t918 * t994) * t966) * m(4)) * qJD(3), t91 * qJD(1) + t8 * qJD(3) + t1011 * qJD(4) + t6 * qJD(5) + t52 * qJD(6) + (-t605 * t638 - t606 * t640 + t453 / 0.2e1 + t454 / 0.2e1) * t1431 + (-t676 * t718 - t677 * t716 + t523 / 0.2e1 + t524 / 0.2e1) * t1434 + (t642 / 0.2e1 + t643 / 0.2e1 + (-t894 * t995 + t895 * t994) * t943) * t1437, t129 * qJD(1) + t1 * qJD(3) + t6 * qJD(4) + t37 * qJD(6) + (t1659 * t603 - t512 * t604 + t519 * t527 + t520 * t528) * t1430 + ((t395 / 0.2e1 + t393 / 0.2e1 + t483 / 0.2e1 + t482 / 0.2e1) * t995 + (t481 / 0.2e1 + t480 / 0.2e1 + t394 / 0.2e1 + t392 / 0.2e1) * t994) * t1218 + ((-t536 - t535) * t997 + (t579 * t615 + t580 * t616 - t589 * t738 - t591 * t743) * m(6)) * qJD(5), t46 * qJD(3) + t52 * qJD(4) + t37 * qJD(5) + t369 * t1433 / 0.4e1; 0, 0 (m(4) * t1080 + (t1227 + t1252) * t1508 + (t1113 + t1227) * t1644 + (t1056 + t1227) * t1504) * qJD(3) + t1648, qJD(3) * t321 + t1648 ((m(6) * t738 + t1109 * t1504) * t995 + (-m(6) * t743 - t604 * t1197 / 0.2e1) * t994) * t1218 + t633 * qJD(6) + t1196 * t138, qJD(5) * t633 + t1196 * t661; t13 * qJD(3) + t10 * qJD(4) + t2 * qJD(5) + t45 * qJD(6) + (-t274 / 0.4e1 + t1288 * t611) * t1433 + t360 * t1643 + (t1010 - (t961 + t962) * t1005 / 0.2e1 + t1554) * qJD(1), t1646, t13 * qJD(1) + (m(5) * (-t1563 * t838 - t1564 * t836 + t532 * t634) + 0.4e1 * (t1225 * t966 * t968 + t1080 * (t994 * (rSges(4,1) * t1210 - t1226) + t995 * (rSges(4,1) * t1209 + t994 * rSges(4,3) - t974))) * t1509 + (t991 * t912 + (-t994 * t911 + t1549) * t995) * t1486 + (t992 * t911 + (-t995 * t912 + t1549) * t994) * t1485 + (t350 * t441 - t611 * t612 - t613 * t614) * m(7) + (t403 * t516 - t696 * t697 - t698 * t699) * m(6) + t1057) * qJD(3) + t18 * qJD(4) + t21 * qJD(5) + t245 * t1185, t10 * qJD(1) + t18 * qJD(3) + t1057 * qJD(4) + t14 * qJD(5) + t29 * qJD(6) + (-t168 / 0.4e1 + t1040 + t1103) * t1431 + (-t265 / 0.4e1 + t1039 + t1102) * t1434 + (-t502 / 0.4e1 + t1101 + t1558) * t1437, t2 * qJD(1) - t1208 + t21 * qJD(3) + t14 * qJD(4) + (t1012 + (t1524 * t615 + t1525 * t616 + t1536 * t403 + t1177) * t1644 + (t1530 * t519 + t1531 * t520 + t1542 * t350 + t1180) * t1504) * qJD(5) + t24 * qJD(6), t45 * qJD(1) + t245 * t1432 / 0.4e1 + t29 * qJD(4) + t24 * qJD(5) + t1297; (t1010 - t1473 / 0.4e1) * qJD(1) + t9 * qJD(3) + t1020 * qJD(4) + t4 * qJD(5) + t51 * qJD(6) + (-t296 / 0.4e1 + t1288 * t638) * t1433 + t378 * t1643, t1646, t9 * qJD(1) + t1057 * qJD(3) + t19 + t15 * qJD(5) + t30 * qJD(6) + (-t157 / 0.4e1 + t397 * t441 - t612 * t638 - t614 * t640 + t1103) * t1432 + (-t238 / 0.4e1 + t508 * t516 - t697 * t716 - t699 * t718 + t1102) * t1435 + (t634 * t646 - t391 / 0.4e1 + (-t836 * t994 - t838 * t995) * t943 + t1101) * t1438, qJD(1) * t1020 + t20 * qJD(3) + t22 * qJD(5) + t1185 * t271 + t19, t4 * qJD(1) - t1208 + t15 * qJD(3) + t22 * qJD(4) + (t1012 + (t1521 * t615 + t1522 * t616 + t1536 * t508 + t1177) * t1644 + (t1528 * t519 + t1529 * t520 + t1542 * t397 + t1180) * t1504) * qJD(5) + t27 * qJD(6), t51 * qJD(1) + t30 * qJD(3) + t271 * t1431 / 0.4e1 + t27 * qJD(5) + t1297; t3 * qJD(3) + t5 * qJD(4) + t36 * qJD(6) + (-t275 / 0.4e1 + t1288 * t512) * t1433 + t390 * t1643 + t1656 * qJD(1), -qJD(6) * t632 + t1196 * t137, t3 * qJD(1) + t1014 * qJD(3) + t7 * qJD(4) + t23 * qJD(6) + (-t167 / 0.4e1 + t417 * t441 + t1659 * t614 - t512 * t612 + t174 / 0.2e1 + t334 / 0.2e1 + t335 / 0.2e1) * t1432 + (-t263 / 0.4e1 + t516 * t551 + t589 * t699 - t591 * t697 + t273 / 0.2e1 + t438 / 0.2e1 + t439 / 0.2e1) * t1435 + t1416, t5 * qJD(1) + t7 * qJD(3) + ((-t194 / 0.4e1 + t35 / 0.2e1) * m(7) + (-t317 / 0.4e1 + t90 / 0.2e1) * m(6) + t1014) * qJD(4) + t26 * qJD(6) + t1416, t1196 * t12 + (t1605 * t1486 + t1604 * t1484 + ((t393 + t395) * t995 + (t392 + t394) * t994) * t1481) * t1218 + t1488 * t1429 + ((t551 * t571 + t589 * t615 - t591 * t616) * m(6) + (t536 / 0.2e1 + t535 / 0.2e1) * t997 ^ 2 + (t1659 * t519 + t417 * t450 - t512 * t520) * m(7)) * qJD(5), t36 * qJD(1) - t1204 + t23 * qJD(3) + t26 * qJD(4) + t1430 * t1488 + (t665 * t1216 + t901 * t755 + t899 * t756 - t592 / 0.4e1) * t1429; (-t369 / 0.4e1 - t1288 * t899) * t1433 + t47 * qJD(3) + t53 * qJD(4) + t38 * qJD(5), qJD(5) * t632 + t1196 * t660, t47 * qJD(1) + (t441 * t1216 + t899 * t612 + t901 * t614 + t346 / 0.2e1 + t601 / 0.2e1 - t602 / 0.2e1 - t245 / 0.4e1) * t1432 + t31 * qJD(4) + t25 * qJD(5) + t1280, t53 * qJD(1) + t31 * qJD(3) + (t161 / 0.2e1 - t271 / 0.4e1) * t1431 + t28 * qJD(5) + t1280, t38 * qJD(1) + t1204 + t25 * qJD(3) + t28 * qJD(4) + (t902 * t1659 - t900 * t512 + t901 * t519 + t899 * t520 - t191 / 0.4e1 + (t1001 * t450 + t1004 * t417) * t996) * t1430 + t592 * t1185 (t592 * qJD(5) / 0.4e1 + (qJD(3) / 0.4e1 + qJD(4) / 0.4e1) * t617) * m(7);];
Cq  = t11;
