% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:51:21
% EndTime: 2019-03-09 03:52:36
% DurationCPUTime: 64.19s
% Computational Cost: add. (218588->1320), mult. (166592->1869), div. (0->0), fcn. (180383->11), ass. (0->806)
t1080 = 2 * m(7);
t887 = -pkin(8) - qJ(4);
t888 = sin(qJ(1));
t1188 = t887 * t888;
t884 = sin(pkin(11));
t889 = cos(qJ(1));
t1190 = t884 * t889;
t881 = pkin(10) + qJ(3);
t870 = sin(t881);
t1096 = pkin(4) * t1190 + t870 * t1188;
t885 = cos(pkin(11));
t865 = t885 * pkin(4) + pkin(3);
t872 = cos(t881);
t1209 = t865 * t872;
t1192 = t872 * t888;
t880 = pkin(11) + qJ(5);
t871 = cos(t880);
t835 = pkin(5) * t871 + t865;
t799 = t835 * t1192;
t869 = sin(t880);
t1300 = pkin(5) * t869;
t1301 = pkin(4) * t884;
t838 = t1300 + t1301;
t879 = -pkin(9) + t887;
t563 = t889 * t838 - t799 + (t870 * t879 + t1209) * t888 - t1096;
t1205 = t870 * t888;
t873 = qJ(6) + t880;
t863 = sin(t873);
t1183 = t888 * t863;
t864 = cos(t873);
t767 = t1183 * t872 + t864 * t889;
t1182 = t888 * t864;
t768 = t1182 * t872 - t863 * t889;
t1423 = -t768 * rSges(7,1) + t767 * rSges(7,2);
t604 = rSges(7,3) * t1205 - t1423;
t1136 = -t563 + t604;
t1094 = t879 - t887;
t1098 = t835 - t865;
t667 = t1094 * t872 + t1098 * t870;
t627 = t667 * t1205;
t1283 = rSges(7,2) * t863;
t1285 = rSges(7,1) * t864;
t980 = -t1283 + t1285;
t731 = -rSges(7,3) * t872 + t870 * t980;
t691 = t731 * t1205;
t370 = t1136 * t872 + t627 + t691;
t1230 = t563 * t872;
t1477 = t872 * t604 + t691;
t371 = t1477 + t627 - t1230;
t1152 = t370 - t371;
t1485 = t1152 * t1080;
t1484 = -m(6) / 0.4e1;
t1041 = t1205 / 0.4e1;
t1043 = -t1205 / 0.4e1;
t1483 = t1041 + t1043;
t1292 = m(7) * qJD(1);
t1481 = -t1292 / 0.4e1;
t1480 = qJD(1) * t1484;
t1204 = t870 * t889;
t1044 = -t1205 / 0.2e1;
t967 = Icges(7,5) * t864 - Icges(7,6) * t863;
t725 = -Icges(7,3) * t872 + t870 * t967;
t1263 = Icges(7,4) * t864;
t971 = -Icges(7,2) * t863 + t1263;
t727 = -Icges(7,6) * t872 + t870 * t971;
t1264 = Icges(7,4) * t863;
t974 = Icges(7,1) * t864 - t1264;
t729 = -Icges(7,5) * t872 + t870 * t974;
t1191 = t872 * t889;
t769 = -t1191 * t863 + t1182;
t770 = t1191 * t864 + t1183;
t456 = t1204 * t725 + t769 * t727 + t770 * t729;
t595 = Icges(7,5) * t768 - Icges(7,6) * t767 + Icges(7,3) * t1205;
t760 = Icges(7,4) * t768;
t598 = -Icges(7,2) * t767 + Icges(7,6) * t1205 + t760;
t759 = Icges(7,4) * t767;
t602 = -Icges(7,1) * t768 - Icges(7,5) * t1205 + t759;
t375 = t595 * t1204 + t769 * t598 - t770 * t602;
t597 = Icges(7,5) * t770 + Icges(7,6) * t769 + Icges(7,3) * t1204;
t1265 = Icges(7,4) * t770;
t600 = Icges(7,2) * t769 + Icges(7,6) * t1204 + t1265;
t761 = Icges(7,4) * t769;
t603 = Icges(7,1) * t770 + Icges(7,5) * t1204 + t761;
t376 = t597 * t1204 + t769 * t600 + t770 * t603;
t965 = t888 * t375 + t376 * t889;
t1474 = -t456 * t872 + t870 * t965;
t1479 = t1474 * t1044;
t968 = Icges(6,5) * t871 - Icges(6,6) * t869;
t736 = -Icges(6,3) * t872 + t870 * t968;
t1266 = Icges(6,4) * t871;
t972 = -Icges(6,2) * t869 + t1266;
t738 = -Icges(6,6) * t872 + t870 * t972;
t1267 = Icges(6,4) * t869;
t975 = Icges(6,1) * t871 - t1267;
t740 = -Icges(6,5) * t872 + t870 * t975;
t1181 = t888 * t869;
t1201 = t871 * t889;
t795 = t1181 * t872 + t1201;
t1180 = t888 * t871;
t796 = t1180 * t872 - t869 * t889;
t464 = t1205 * t736 - t738 * t795 + t740 * t796;
t614 = Icges(6,5) * t796 - Icges(6,6) * t795 + Icges(6,3) * t1205;
t781 = Icges(6,4) * t796;
t617 = -Icges(6,2) * t795 + Icges(6,6) * t1205 + t781;
t780 = Icges(6,4) * t795;
t621 = -Icges(6,1) * t796 - Icges(6,5) * t1205 + t780;
t389 = t1205 * t614 - t617 * t795 - t621 * t796;
t797 = -t869 * t1191 + t1180;
t798 = t1191 * t871 + t1181;
t616 = Icges(6,5) * t798 + Icges(6,6) * t797 + Icges(6,3) * t1204;
t1268 = Icges(6,4) * t798;
t619 = Icges(6,2) * t797 + Icges(6,6) * t1204 + t1268;
t782 = Icges(6,4) * t797;
t622 = Icges(6,1) * t798 + Icges(6,5) * t1204 + t782;
t390 = t616 * t1205 - t795 * t619 + t796 * t622;
t964 = t389 * t888 + t390 * t889;
t193 = -t464 * t872 + t870 * t964;
t454 = t1205 * t725 - t727 * t767 + t729 * t768;
t373 = t1205 * t595 - t598 * t767 - t602 * t768;
t374 = t597 * t1205 - t767 * t600 + t768 * t603;
t966 = t373 * t888 + t374 * t889;
t187 = -t454 * t872 + t870 * t966;
t1227 = t614 * t872;
t1466 = t617 * t869 + t621 * t871;
t416 = t1466 * t870 + t1227;
t1229 = t595 * t872;
t1467 = t598 * t863 + t602 * t864;
t408 = t1467 * t870 + t1229;
t1472 = -t375 * t889 + t376 * t888;
t1478 = t1483 * t1472;
t1422 = -t796 * rSges(6,1) + t795 * rSges(6,2);
t623 = rSges(6,3) * t1205 - t1422;
t1284 = rSges(6,2) * t869;
t1286 = rSges(6,1) * t871;
t983 = -t1284 + t1286;
t742 = -rSges(6,3) * t872 + t870 * t983;
t1476 = t742 * t1205 + t623 * t872;
t466 = t1204 * t736 + t797 * t738 + t798 * t740;
t391 = t614 * t1204 + t797 * t617 - t798 * t621;
t392 = t616 * t1204 + t797 * t619 + t798 * t622;
t963 = t888 * t391 + t392 * t889;
t1475 = -t466 * t872 + t870 * t963;
t1473 = -t391 * t889 + t392 * t888;
t1471 = -t870 / 0.2e1;
t1362 = -t884 / 0.2e1;
t1361 = t885 / 0.2e1;
t1463 = t1080 / 0.4e1;
t1034 = t1192 / 0.4e1;
t1465 = 0.2e1 * m(6);
t1470 = -t1465 / 0.4e1;
t1469 = t1465 / 0.4e1;
t882 = t888 ^ 2;
t883 = t889 ^ 2;
t1093 = t882 + t883;
t1116 = t729 + (-Icges(7,2) * t864 - t1264) * t870;
t1468 = t1116 * t863;
t1382 = m(5) / 0.2e1;
t1381 = m(5) / 0.4e1;
t1380 = m(6) / 0.4e1;
t1379 = m(7) / 0.4e1;
t1364 = t870 / 0.2e1;
t1363 = -t872 / 0.2e1;
t1360 = t888 / 0.2e1;
t1359 = t888 / 0.4e1;
t1358 = -t889 / 0.2e1;
t1357 = -t889 / 0.4e1;
t687 = t727 * t888;
t689 = t729 * t888;
t936 = -t725 * t888 + t1467;
t293 = -t936 * t872 + (t687 * t863 - t689 * t864 + t595) * t870;
t728 = Icges(7,6) * t870 + t872 * t971;
t730 = Icges(7,5) * t870 + t872 * t974;
t1219 = t725 * t872;
t726 = Icges(7,3) * t870 + t872 * t967;
t1212 = t864 * t729;
t1215 = t863 * t727;
t947 = t1212 - t1215;
t931 = t726 - t947;
t902 = t870 * t931 + t1219;
t335 = -t728 * t767 + t730 * t768 + t888 * t902;
t1454 = t293 + t335;
t688 = t727 * t889;
t690 = t729 * t889;
t955 = -t600 * t863 + t603 * t864;
t935 = -t725 * t889 - t955;
t294 = -t935 * t872 + (t688 * t863 - t690 * t864 + t597) * t870;
t336 = t769 * t728 + t770 * t730 + t889 * t902;
t1453 = t294 + t336;
t1451 = -t408 + t454;
t1228 = t597 * t872;
t409 = t870 * t955 - t1228;
t1450 = t409 + t456;
t1269 = Icges(4,4) * t870;
t828 = Icges(4,1) * t872 - t1269;
t776 = Icges(4,5) * t888 + t828 * t889;
t733 = t776 * t1192;
t824 = Icges(4,5) * t872 - Icges(4,6) * t870;
t772 = Icges(4,3) * t888 + t824 * t889;
t1011 = t889 * t772 - t733;
t771 = Icges(4,5) * t1192 - Icges(4,6) * t1205 - Icges(4,3) * t889;
t847 = Icges(4,4) * t1205;
t775 = Icges(4,1) * t1192 - Icges(4,5) * t889 - t847;
t1115 = -t775 * t1191 - t888 * t771;
t773 = Icges(4,4) * t1192 - Icges(4,2) * t1205 - Icges(4,6) * t889;
t1259 = Icges(4,2) * t870;
t861 = Icges(4,4) * t872;
t774 = Icges(4,6) * t888 + (t861 - t1259) * t889;
t1449 = -t1204 * t773 - t1205 * t774 - t1011 - t1115;
t643 = t769 * rSges(7,1) - t770 * rSges(7,2);
t784 = t797 * pkin(5);
t1122 = -t643 - t784;
t1448 = t1122 * t888;
t1287 = rSges(5,1) * t885;
t986 = -rSges(5,2) * t884 + t1287;
t927 = t986 * t870;
t1443 = t872 * rSges(5,3) - t927;
t1297 = -pkin(7) - qJ(2);
t1027 = t889 * t1297;
t1282 = rSges(6,3) * t870;
t866 = cos(pkin(10)) * pkin(2) + pkin(1);
t941 = t1209 + t866 + t1282;
t527 = -t888 * t941 - t1027 + t1096 + t1422;
t1273 = rSges(5,3) + qJ(4);
t1299 = t872 * pkin(3);
t1407 = t1273 * t870 + t1299 + t866;
t1179 = t888 * t884;
t1189 = t885 * t889;
t817 = t1179 * t872 + t1189;
t1178 = t888 * t885;
t818 = t1178 * t872 - t1190;
t1421 = -t818 * rSges(5,1) + t817 * rSges(5,2);
t561 = -t1407 * t888 - t1027 + t1421;
t1272 = -rSges(7,3) + t879;
t488 = (t1272 * t870 - t866) * t888 + (t838 - t1297) * t889 - t799 + t1423;
t831 = pkin(3) * t870 - qJ(4) * t872;
t1111 = t831 - t1443;
t673 = t1111 * t889;
t1220 = t673 * t889;
t1169 = qJ(4) + t887;
t1298 = pkin(3) - t865;
t1418 = t1298 * t870;
t1119 = t1169 * t872 - t1418 + t831;
t1050 = t742 + t1119;
t549 = t1050 * t889;
t1231 = t549 * t889;
t1121 = t667 + t731;
t1002 = t1119 + t1121;
t470 = t1002 * t889;
t1238 = t470 * t889;
t468 = t1002 * t888;
t1239 = t468 * t888;
t467 = -0.2e1 * t1239;
t547 = t1050 * t888;
t1232 = t547 * t888;
t540 = -0.2e1 * t1232;
t671 = t1111 * t888;
t1221 = t671 * t888;
t646 = -0.2e1 * t1221;
t832 = rSges(4,1) * t870 + rSges(4,2) * t872;
t997 = 0.2e1 * t1093;
t1006 = (t467 - 0.2e1 * t1238) * t1379 + (t540 - 0.2e1 * t1231) * t1380 + (t646 - 0.2e1 * t1220) * t1381 - m(4) * t832 * t997 / 0.4e1;
t1052 = (t467 + 0.2e1 * t1239) * t1379 + (t540 + 0.2e1 * t1232) * t1380 + (t646 + 0.2e1 * t1221) * t1381;
t1442 = t1006 - t1052;
t1075 = -0.2e1 * t1205;
t1071 = 0.2e1 * t1191;
t1072 = 0.2e1 * t1192;
t867 = t888 * t1297;
t819 = -t1190 * t872 + t1178;
t820 = t1189 * t872 + t1179;
t987 = t820 * rSges(5,1) + t819 * rSges(5,2);
t562 = t1407 * t889 - t867 + t987;
t1138 = t561 * t1071 + t562 * t1072;
t984 = t798 * rSges(6,1) + t797 * rSges(6,2);
t841 = t887 * t1204;
t996 = pkin(4) * t1179 - t841;
t528 = t889 * t941 - t867 + t984 + t996;
t1143 = t527 * t1071 + t528 * t1072;
t1280 = rSges(7,3) * t870;
t837 = t879 * t1204;
t981 = t770 * rSges(7,1) + t769 * rSges(7,2);
t489 = t888 * t838 - t837 - t867 + (t835 * t872 + t1280 + t866) * t889 + t981;
t1146 = t488 * t1071 + t489 * t1072;
t1064 = t468 * t1204;
t463 = -0.2e1 * t1064;
t1063 = t547 * t1204;
t537 = -0.2e1 * t1063;
t1062 = t671 * t1204;
t626 = -0.2e1 * t1062;
t1054 = (-t1075 * t470 + t1146 + t463) * t1379 + (-t1075 * t549 + t1143 + t537) * t1380 + (-t1075 * t673 + t1138 + t626) * t1381;
t1056 = (0.2e1 * t1064 + t463) * t1379 + (0.2e1 * t1063 + t537) * t1380 + (0.2e1 * t1062 + t626) * t1381;
t1441 = t1054 - t1056;
t1012 = t1098 * t872;
t564 = -t837 + t841 + (t838 - t1301) * t888 + t889 * t1012;
t606 = rSges(7,3) * t1204 + t981;
t1135 = t564 + t606;
t1420 = t1135 * t872;
t372 = t1121 * t1204 + t1420;
t1244 = t372 * t888;
t1413 = t371 * t889 - t1244;
t222 = 0.2e1 * t1413;
t625 = rSges(6,3) * t1204 + t984;
t535 = t1204 * t742 + t872 * t625;
t1233 = t535 * t888;
t1415 = t1476 * t889 - t1233;
t412 = 0.2e1 * t1415;
t1163 = t1379 * t222 + t1380 * t412;
t1384 = -0.2e1 * t889;
t1167 = (t1384 * t370 + 0.2e1 * t1244 + t222) * t1379 + (t1384 * t1476 + 0.2e1 * t1233 + t412) * t1380;
t1440 = t1163 - t1167;
t1074 = -0.2e1 * t1204;
t214 = t372 * t1074 + t371 * t1075;
t404 = t535 * t1074 + t1075 * t1476;
t1164 = t1379 * t214 + t1380 * t404;
t1078 = 0.2e1 * t870;
t1168 = ((t370 * t888 + t372 * t889) * t1078 + t214) * t1379 + ((t1476 * t888 + t535 * t889) * t1078 + t404) * t1380;
t1439 = t1164 - t1168;
t1073 = -0.2e1 * t1192;
t1142 = t1071 * t1476 + t535 * t1073;
t1153 = t371 * t1071 + t372 * t1073;
t822 = t1093 * t870;
t1432 = 0.2e1 * t822;
t578 = t604 * t1204;
t319 = t578 + (-t1135 * t888 - t563 * t889) * t870;
t952 = t623 * t889 - t625 * t888;
t483 = t952 * t870;
t1166 = (t1432 * t319 + t1153) * t1379 + (t1432 * t483 + t1142) * t1380;
t1385 = -0.2e1 * t872;
t1101 = rSges(7,3) * t1191 + t1204 * t1283;
t695 = -t1204 * t1285 + t1101;
t1128 = -t667 * t889 + t695;
t694 = t731 * t888;
t1134 = t604 * t1191 - t694 * t1204;
t1100 = t872 * t1188 + t865 * t1205;
t1194 = t872 * t879;
t608 = (-t835 * t870 - t1194) * t888 + t1100;
t209 = (t608 * t870 - t1230) * t889 + (-t1128 * t870 - t1420) * t888 + t1134;
t732 = t872 * t980 + t1280;
t1051 = t731 * t1192 + t732 * t1205 - t872 * t694;
t668 = -t1094 * t870 + t1012;
t239 = (t667 * t888 + t608) * t872 + (t668 * t888 - t1136) * t870 + t1051;
t1120 = -t668 - t732;
t586 = t870 * t606;
t240 = t586 + (t1120 * t889 + t564) * t870 + (-t1121 * t889 - t1128) * t872;
t1217 = t742 * t888;
t1099 = rSges(6,3) * t1191 + t1204 * t1284;
t710 = -rSges(6,1) * t1201 * t870 + t1099;
t401 = t952 * t872 + (-t1217 * t889 - t710 * t888) * t870;
t743 = t872 * t983 + t1282;
t436 = (t743 * t888 - t623) * t870;
t437 = (-t742 * t889 - t710) * t872 + (-t743 * t889 + t625) * t870;
t1271 = (t401 * t1385 + (t436 * t889 + t437 * t888 + t483) * t1078 + t1142) * t1380 + (t209 * t1385 + (t239 * t889 + t240 * t888 + t319) * t1078 + t1153) * t1379;
t1438 = t1166 - t1271;
t1026 = t1298 * t872;
t1252 = qJ(4) * t870;
t1437 = t1252 + t1026;
t642 = -t767 * rSges(7,1) - t768 * rSges(7,2);
t783 = t795 * pkin(5);
t571 = -t642 + t783;
t669 = -rSges(6,1) * t795 - rSges(6,2) * t796;
t670 = rSges(6,1) * t797 - rSges(6,2) * t798;
t762 = (-rSges(7,1) * t863 - rSges(7,2) * t864) * t870;
t995 = t1300 * t870 - t762;
t683 = t995 * t888;
t684 = t995 * t889;
t785 = (-rSges(6,1) * t869 - rSges(6,2) * t871) * t870;
t958 = t527 * t889 + t528 * t888;
t1436 = (t1122 * t468 - t470 * t571 + t488 * t684 + t489 * t683) * t1463 + (-t547 * t670 + t549 * t669 - t785 * t958) * t1469;
t573 = (t1272 * t872 + (t835 + t980) * t870) * t888;
t574 = (-t1194 + (-t835 - t1285) * t870) * t889 + t1101;
t592 = t1100 + t1217;
t1193 = t872 * t887;
t593 = (-t1193 + (-t865 - t1286) * t870) * t889 + t1099;
t1435 = (t1476 * t592 + t436 * t527 + t437 * t528 - t535 * t593) * t1470 - (t239 * t488 + t240 * t489 + t371 * t573 - t372 * t574) * t1080 / 0.4e1;
t1434 = t468 * t1485 / 0.4e1;
t651 = Icges(5,4) * t818 - Icges(5,2) * t817 + Icges(5,6) * t1205;
t654 = Icges(5,1) * t818 - Icges(5,4) * t817 + Icges(5,5) * t1205;
t951 = -t651 * t817 + t654 * t818;
t756 = (-Icges(7,5) * t863 - Icges(7,6) * t864) * t870;
t1196 = t872 * t756;
t1433 = -t1196 / 0.2e1 + t1468 * t1471;
t1356 = t889 / 0.2e1;
t1430 = m(7) * t1078;
t969 = Icges(5,5) * t885 - Icges(5,6) * t884;
t1427 = t870 * t969;
t1424 = -t888 * t642 - t889 * t643;
t1009 = t1381 + t1380 + t1379;
t1076 = 0.4e1 * t872;
t494 = t1009 * (-0.1e1 + t1093) * t870 * t1076;
t1197 = t872 * t606;
t520 = t1204 * t731 + t1197;
t1235 = t520 * t888;
t1414 = t1477 * t889 - t1235;
t1131 = -Icges(7,2) * t768 - t602 - t759;
t1133 = -Icges(7,1) * t767 - t598 - t760;
t636 = -Icges(7,5) * t767 - Icges(7,6) * t768;
t276 = -t1131 * t767 + t1133 * t768 + t1205 * t636;
t1130 = -Icges(7,2) * t770 + t603 + t761;
t1132 = Icges(7,1) * t769 - t1265 - t600;
t637 = Icges(7,5) * t769 - Icges(7,6) * t770;
t277 = -t1130 * t767 + t1132 * t768 + t1205 * t637;
t158 = -t276 * t889 + t277 * t888;
t278 = t1131 * t769 + t1133 * t770 + t1204 * t636;
t279 = t1130 * t769 + t1132 * t770 + t1204 * t637;
t159 = -t278 * t889 + t279 * t888;
t1165 = t1358 * t158 + t1360 * t159;
t1211 = t864 * t730;
t1214 = t863 * t728;
t364 = -t931 * t872 + (t1211 + t725 - t1214) * t870;
t484 = t870 * t947 - t1219;
t962 = -t408 * t888 + t409 * t889;
t1411 = (-t484 * t872 + t870 * t962) * t1364 + ((-t364 + t962) * t872 + (t293 * t888 + t294 * t889 + t484) * t870) * t1363;
t1410 = t1363 * t364 + t484 * t1364;
t529 = t888 * t669 + t670 * t889;
t1147 = (t1122 * t889 + t888 * t571) * t1463 + t529 * t1470;
t950 = t669 * t889 - t670 * t888;
t1148 = (t571 * t889 - t1448) * t1430 / 0.4e1 + t950 * t1078 * t1484;
t451 = -t888 * t783 + t784 * t889 - t1424;
t978 = t870 * t997;
t1154 = (t451 * t1385 + (t683 * t888 + t684 * t889) * t1078) * t1379 + (t1385 * t529 - t785 * t978) * t1380;
t607 = t642 * t1204;
t490 = -t1205 * t643 + t607;
t538 = t762 * t1205 + t872 * t642;
t539 = -t1204 * t762 - t872 * t643;
t990 = t319 * t490 + t371 * t538 - t372 * t539;
t973 = Icges(5,4) * t885 - Icges(5,2) * t884;
t1409 = -Icges(5,6) * t872 + t870 * t973;
t976 = Icges(5,1) * t885 - Icges(5,4) * t884;
t1408 = -Icges(5,5) * t872 + t870 * t976;
t1203 = t871 * t740;
t1208 = t869 * t738;
t1270 = Icges(4,1) * t870;
t1288 = rSges(4,1) * t872;
t1016 = t866 + t1288;
t1095 = rSges(4,2) * t1205 + t889 * rSges(4,3);
t700 = -t1016 * t888 - t1027 + t1095;
t1015 = -rSges(4,2) * t1204 + t888 * rSges(4,3);
t701 = t1016 * t889 + t1015 - t867;
t737 = Icges(6,3) * t870 + t872 * t968;
t814 = t832 * t888;
t816 = t832 * t889;
t1406 = -m(4) * (t700 * t814 - t701 * t816) - (t861 - t1259 / 0.2e1 + t1270 / 0.2e1 + Icges(5,3) * t1471 + t969 * t1363 - t737 / 0.2e1 - t726 / 0.2e1 + t1408 * t1361 + t1409 * t1362 + t1203 / 0.2e1 - t1208 / 0.2e1 + t1212 / 0.2e1 - t1215 / 0.2e1) * t872;
t854 = pkin(3) * t1205;
t644 = t854 + (-t1273 * t872 + t927) * t888;
t1097 = t870 * rSges(5,2) * t1190 + rSges(5,3) * t1191;
t844 = qJ(4) * t1191;
t645 = t844 + (-pkin(3) - t1287) * t1204 + t1097;
t1053 = ((t573 * t889 + t574 * t888) * t1078 + t1146) * t1379 + ((t592 * t889 + t593 * t888) * t1078 + t1143) * t1380 + ((t644 * t889 + t645 * t888) * t1078 + t1138) * t1381;
t977 = -t861 - t1270;
t812 = t977 * t888;
t813 = t977 * t889;
t1405 = ((t773 - t812) * t889 + (-t774 + t813) * t888) * t872;
t778 = (-Icges(6,2) * t871 - t1267) * t870;
t779 = (-Icges(6,1) * t869 - t1266) * t870;
t1404 = (t779 / 0.2e1 - t738 / 0.2e1) * t871 - (t740 / 0.2e1 + t778 / 0.2e1) * t869;
t1125 = -Icges(6,2) * t796 - t621 - t780;
t1127 = -Icges(6,1) * t795 - t617 - t781;
t661 = -Icges(6,5) * t795 - Icges(6,6) * t796;
t323 = -t661 * t872 + (-t1125 * t869 + t1127 * t871) * t870;
t1124 = -Icges(6,2) * t798 + t622 + t782;
t1126 = Icges(6,1) * t797 - t1268 - t619;
t662 = Icges(6,5) * t797 - Icges(6,6) * t798;
t324 = -t662 * t872 + (-t1124 * t869 + t1126 * t871) * t870;
t1112 = t740 + t778;
t1113 = -t738 + t779;
t777 = (-Icges(6,5) * t869 - Icges(6,6) * t871) * t870;
t386 = -t1112 * t795 + t1113 * t796 + t1205 * t777;
t387 = t1112 * t797 + t1113 * t798 + t1204 * t777;
t1403 = t1436 + (t324 + t387) * t1359 + (t323 + t386) * t1357;
t315 = -t637 * t872 + (-t1130 * t863 + t1132 * t864) * t870;
t1250 = t315 * t888;
t314 = -t636 * t872 + (-t1131 * t863 + t1133 * t864) * t870;
t1251 = t314 * t889;
t758 = (-Icges(7,1) * t863 - t1263) * t870;
t1117 = -t727 + t758;
t368 = -t1116 * t767 + t1117 * t768 + t1205 * t756;
t369 = t1116 * t769 + t1117 * t770 + t1204 * t756;
t1007 = -t1251 / 0.4e1 + t1250 / 0.4e1 + t369 * t1359 + t368 * t1357;
t944 = t888 * t814 + t816 * t889;
t1005 = (t888 * t573 - t574 * t889) * t1463 + (t888 * t592 - t593 * t889) * t1469 + (t888 * t644 - t645 * t889) * t1382 + m(4) * t944 / 0.2e1;
t1176 = t889 * t1474;
t1186 = t888 * t187;
t1402 = t1176 / 0.4e1 + t1186 / 0.4e1 + t1474 * t1357 - t187 * t1359;
t1175 = t889 * t1475;
t1185 = t888 * t193;
t1400 = t1185 / 0.4e1 + t1175 / 0.4e1 - t193 * t1359 + t1475 * t1357 - t1434 + t1483 * t1473;
t1031 = t1191 / 0.4e1;
t1037 = t1204 / 0.4e1;
t704 = t738 * t888;
t706 = t740 * t888;
t934 = -t736 * t888 + t1466;
t307 = -t934 * t872 + (t704 * t869 - t706 * t871 + t614) * t870;
t705 = t738 * t889;
t707 = t740 * t889;
t953 = -t619 * t869 + t622 * t871;
t933 = -t736 * t889 - t953;
t308 = -t933 * t872 + (t705 * t869 - t707 * t871 + t616) * t870;
t739 = Icges(6,6) * t870 + t872 * t972;
t741 = Icges(6,5) * t870 + t872 * t975;
t1218 = t736 * t872;
t946 = t1203 - t1208;
t930 = t737 - t946;
t901 = t870 * t930 + t1218;
t347 = -t739 * t795 + t741 * t796 + t888 * t901;
t348 = t797 * t739 + t798 * t741 + t889 * t901;
t1202 = t871 * t741;
t1207 = t869 * t739;
t380 = -t930 * t872 + (t1202 + t736 - t1207) * t870;
t1226 = t616 * t872;
t417 = t870 * t953 - t1226;
t493 = t870 * t946 - t1218;
t1399 = t380 * t1363 + t493 * t1364 - t1435 + (t307 + t347) * t1041 + (t308 + t348) * t1037 + (-t416 + t464) * t1034 + (t417 + t466) * t1031;
t1392 = 0.2e1 * t1477;
t1391 = -0.2e1 * t520;
t1390 = 0.2e1 * t538;
t1389 = -0.2e1 * t642;
t1388 = -0.2e1 * t643;
t1386 = 0.4e1 * t822;
t378 = (-t695 * t870 - t1197) * t888 + t1134;
t424 = -t604 * t870 + t1051;
t425 = -t732 * t1204 + t586 + (-t731 * t889 - t695) * t872;
t473 = -t1205 * t606 + t578;
t1377 = (t1477 * t239 + t209 * t473 - t240 * t520 + t319 * t378 + t371 * t424 - t372 * t425) * t1080;
t1375 = t520 * t1485;
t833 = t1252 + t1299;
t1104 = t1093 * t833;
t1003 = -t888 * (t1437 * t888 + t1096) + t1104 + t889 * (-t1437 * t889 + t996);
t246 = t1135 * t889 + t1136 * t888 + t1003;
t203 = 0.2e1 * t490 * t246;
t383 = t470 * t1390;
t384 = -0.2e1 * t539 * t468;
t1055 = t203 - t383 + t384;
t1079 = 0.2e1 * t762;
t1372 = m(7) * (-t1079 * t1413 - 0.2e1 * t1424 * t319 + t1055);
t1370 = m(7) * (t1391 * t683 + t1392 * t684 + 0.2e1 * t451 * t473 + t1055);
t1369 = -t1475 / 0.2e1;
t1368 = t323 / 0.2e1;
t648 = Icges(5,5) * t818 - Icges(5,6) * t817 + Icges(5,3) * t1205;
t1367 = -t648 / 0.2e1;
t751 = Icges(5,6) * t870 + t872 * t973;
t1366 = t751 / 0.2e1;
t405 = t488 * t1390;
t406 = 0.2e1 * t539 * t489;
t1151 = t405 + t406;
t1328 = m(7) * (t1388 * t372 + t1389 * t371 + t1151);
t1144 = t1071 * t1477 + t520 * t1073;
t1327 = m(7) * (t378 * t1385 + (t424 * t889 + t425 * t888 + t473) * t1078 + t1144);
t1325 = (t1477 * t573 + t424 * t488 + t425 * t489 - t520 * t574) * t1080;
t1323 = m(7) * (-t1122 * t1391 + t1392 * t571 + t1151);
t959 = t488 * t889 + t489 * t888;
t1320 = m(7) * (-t1079 * t959 + t1388 * t468 - t1389 * t470);
t385 = t520 * t1074 + t1075 * t1477;
t1319 = m(7) * ((t1477 * t888 + t520 * t889) * t1078 + t385);
t400 = 0.2e1 * t1414;
t1315 = m(7) * (t1384 * t1477 + 0.2e1 * t1235 + t400);
t1313 = m(7) * (t1432 * t473 + t1144);
t1309 = m(7) * t385;
t1308 = m(7) * t400;
t1307 = m(7) * (-t1385 * t1424 - t762 * t978);
t1303 = (-t642 * t889 + t643 * t888) * t1430;
t1302 = t1424 * t1080;
t1296 = m(5) * qJD(3);
t1294 = m(6) * qJD(3);
t1293 = m(6) * qJD(5);
t1291 = m(7) * qJD(3);
t1290 = m(7) * qJD(5);
t1289 = m(7) * qJD(6);
t1255 = Icges(5,3) * t872;
t418 = -t1196 + (t1117 * t864 - t1468) * t870;
t1243 = t418 * t872;
t1225 = t651 * t884;
t653 = Icges(5,4) * t820 + Icges(5,2) * t819 + Icges(5,6) * t1204;
t1224 = t653 * t884;
t1223 = t654 * t885;
t656 = Icges(5,1) * t820 + Icges(5,4) * t819 + Icges(5,5) * t1204;
t1222 = t656 * t885;
t1216 = t773 * t870;
t1210 = t864 * t870;
t1195 = t872 * t777;
t922 = (t436 * t888 - t437 * t889) * t1465;
t99 = -t922 / 0.4e1 + ((-t683 / 0.2e1 + t240 / 0.2e1) * t889 + (t684 / 0.2e1 - t239 / 0.2e1) * t888) * m(7);
t1170 = t99 * qJD(2);
t1118 = t1169 * t870 + t1026 - t833;
t1114 = t776 * t1191 + t888 * t772;
t1110 = -rSges(5,3) * t870 - t872 * t986 - t833;
t810 = -Icges(4,2) * t1192 - t847;
t1107 = t775 + t810;
t991 = qJ(4) * t1192 - t854;
t1105 = t888 * t991 + t889 * (-pkin(3) * t1204 + t844);
t948 = t700 * t889 + t701 * t888;
t148 = 0.4e1 * t959 * t1379 + 0.4e1 * t958 * t1380 + 0.4e1 * (t561 * t889 + t562 * t888) * t1381 + t948 * m(4) + m(3) * t1093 * (rSges(3,3) + qJ(2));
t1092 = qJD(1) * t148;
t174 = 0.4e1 * ((-t488 * t888 + t489 * t889) * t1379 + (-t527 * t888 + t528 * t889) * t1380 + (-t561 * t888 + t562 * t889) * t1381) * t870;
t1091 = qJD(1) * t174;
t1090 = qJD(1) * t870;
t1088 = qJD(3) * t870;
t1087 = qJD(3) * t872;
t1086 = qJD(5) * t870;
t1085 = qJD(5) * t872;
t272 = 0.2e1 * t424 * t888 - 0.2e1 * t425 * t889;
t1084 = t272 * qJD(2);
t1008 = m(7) / 0.2e1 + m(6) / 0.2e1 + t1382;
t501 = t1008 * t1432;
t1083 = t501 * qJD(1);
t1070 = t870 ^ 2 * t1300;
t1066 = t272 * t1379;
t1061 = t1290 / 0.2e1;
t1060 = t1289 / 0.2e1;
t1059 = t1369 + t1475 / 0.2e1;
t1038 = t1204 / 0.2e1;
t1042 = t1205 / 0.2e1;
t114 = -t368 * t872 + (t276 * t888 + t277 * t889) * t870;
t115 = -t369 * t872 + (t278 * t888 + t279 * t889) * t870;
t1057 = t115 * t1038 + t114 * t1042 + (-t1243 + (t314 * t888 + t315 * t889) * t870) * t1363;
t650 = Icges(5,5) * t820 + Icges(5,6) * t819 + Icges(5,3) * t1204;
t423 = t650 * t1204 + t819 * t653 + t820 * t656;
t1049 = -t743 + t1118;
t1046 = -t1210 / 0.2e1;
t1045 = t1210 / 0.2e1;
t1035 = t1192 / 0.2e1;
t1032 = t1191 / 0.2e1;
t825 = Icges(4,2) * t872 + t1269;
t811 = t825 * t889;
t1013 = (-t776 + t811) * t888;
t1010 = t774 * t870 - t771;
t1001 = t1118 + t1120;
t1000 = t888 * (-t991 - t1100) + t889 * (-t844 + (-t1193 + t1418) * t889) + t1105;
t475 = -t1070 * t888 - t783 * t872 + t538;
t557 = t1205 * t785 + t669 * t872;
t993 = m(6) * t557 + m(7) * t475;
t476 = (-t762 * t870 + t1070) * t889 + t1122 * t872;
t558 = -t1204 * t785 - t872 * t670;
t992 = m(6) * t558 + m(7) * t476;
t989 = t1474 * t1042 + t1479;
t979 = t758 * t1045 + t727 * t1046 + t1433;
t970 = -Icges(4,5) * t870 - Icges(4,6) * t872;
t961 = -t416 * t888 + t417 * t889;
t960 = t1238 + t1239;
t957 = t1231 + t1232;
t949 = -t1220 - t1221;
t943 = -t1375 / 0.4e1 + t989;
t906 = t870 * t936 + t1229;
t264 = t687 * t767 - t689 * t768 + t888 * t906;
t905 = t870 * t935 + t1228;
t265 = t688 * t767 - t690 * t768 + t888 * t905;
t57 = (-t335 + t966) * t872 + (t264 * t888 + t265 * t889 + t454) * t870;
t266 = -t769 * t687 - t770 * t689 + t889 * t906;
t267 = -t769 * t688 - t770 * t690 + t889 * t905;
t58 = (-t336 + t965) * t872 + (t266 * t888 + t267 * t889 + t456) * t870;
t942 = t1032 * t1474 + t187 * t1035 + t58 * t1038 + t57 * t1042 + t1411;
t911 = t1255 - t1427;
t714 = t911 * t888;
t932 = -t1223 + t714 + t1225;
t295 = -t1125 * t795 + t1127 * t796 + t1205 * t661;
t296 = -t1124 * t795 + t1126 * t796 + t1205 * t662;
t168 = -t295 * t889 + t296 * t888;
t297 = t1125 * t797 + t1127 * t798 + t1204 * t661;
t298 = t1124 * t797 + t1126 * t798 + t1204 * t662;
t169 = -t297 * t889 + t298 * t888;
t929 = t1358 * t168 + t1360 * t169;
t926 = t1377 / 0.4e1 + t942;
t925 = (-t648 * t889 + t650 * t888) * t872;
t715 = t911 * t889;
t924 = (-t1222 + t715 + t1224) * t888;
t923 = t1477 * t538 + t473 * t490 - t520 * t539;
t921 = t1402 + t1478;
t920 = t1479 + (t315 + t369) * t1038 + (t1474 + t314 + t368) * t1042;
t919 = t727 * t1045 + t758 * t1046 - t1433;
t144 = -t264 * t889 + t265 * t888;
t145 = -t266 * t889 + t267 * t888;
t918 = t145 * t1038 + t144 * t1042 + t57 * t1358 + t58 * t1360 + (-t293 * t889 + t294 * t888) * t1363 + (-t373 * t889 + t374 * t888) * t1035 + t1472 * t1032 + (t408 * t889 + t409 * t888) * t1364 - t1165;
t908 = t1031 * t1450 + t1034 * t1451 + t1037 * t1453 + t1041 * t1454 + t1410;
t907 = t159 * t1038 - t58 * t1204 / 0.2e1 + t158 * t1042 + t57 * t1044 + t114 * t1358 + t115 * t1360 - t1411 + (t1250 - t1251 + t1186 + t1176) * t1363;
t904 = t870 * t934 + t1227;
t903 = t870 * t933 + t1226;
t900 = t920 - t1243;
t544 = -t1204 * t774 + t1114;
t899 = -(t1114 + t423) * t888 / 0.2e1 + (t544 + t423) * t1360 + (-t733 + (t772 + t1216) * t889 + t1115 + t1449) * t1358;
t898 = (-t888 * (-t775 * t872 + t1216) - t889 * t771 + t1205 * t648 + t951) * t1358 + (t1010 * t889 - t1114 + t544 + t951) * t1356 + (t1010 * t888 + t648 * t1204 + t1011 + t1449) * t1360;
t894 = t908 + t921 - t1007;
t893 = t1007 - t1410 + t921 + t1454 * t1043 - t1453 * t1204 / 0.4e1 - t1451 * t1192 / 0.4e1 - t1450 * t1191 / 0.4e1;
t892 = t1007 - t1402 + t1478 + t908;
t753 = Icges(5,5) * t870 + t872 * t976;
t891 = t736 / 0.2e1 + t725 / 0.2e1 - t825 / 0.2e1 + t828 / 0.2e1 + t1427 / 0.2e1 - t1255 / 0.2e1 + t753 * t1361 + t751 * t1362 + t1202 / 0.2e1 - t1207 / 0.2e1 + t1211 / 0.2e1 - t1214 / 0.2e1;
t834 = -rSges(4,2) * t870 + t1288;
t809 = t970 * t889;
t808 = t970 * t888;
t719 = t1408 * t889;
t718 = t1408 * t888;
t717 = t1409 * t889;
t716 = t1409 * t888;
t674 = t1110 * t889;
t672 = t1110 * t888;
t550 = t1049 * t889;
t548 = t1049 * t888;
t509 = t950 * t870;
t502 = t1302 / 0.4e1;
t500 = -t1009 * t1432 + (m(5) + m(6) + m(7)) * t1432 / 0.4e1;
t492 = t889 * (-rSges(5,1) * t1189 * t870 + t1097) + t1443 * t882 + t1105;
t486 = t1303 / 0.4e1;
t471 = t1001 * t889;
t469 = t1001 * t888;
t461 = t888 * (rSges(5,3) * t1205 - t1421) + t889 * (rSges(5,3) * t1204 + t987) + t1104;
t439 = -t1195 + (-t1112 * t869 + t1113 * t871) * t870;
t433 = t607 + (-t783 * t889 + t1448) * t870;
t429 = t1307 / 0.4e1;
t426 = 0.2e1 * t538 * t888 - 0.2e1 * t539 * t889;
t419 = t426 * t1060;
t403 = -t1217 * t888 + t710 * t889 + t1000;
t388 = t1308 / 0.4e1;
t379 = t1309 / 0.4e1;
t346 = 0.4e1 * t561 * t644 + 0.4e1 * t562 * t645;
t344 = t888 * t623 + t625 * t889 + t1003;
t341 = t1076 * t949 + t1386 * t461;
t325 = -0.4e1 * t527 * t669 + 0.4e1 * t528 * t670;
t317 = t490 * t1385 + (t538 * t889 + t539 * t888) * t1078;
t316 = t317 * t1060;
t313 = 0.4e1 * t527 * t592 + 0.4e1 * t528 * t593;
t302 = t1128 * t889 + (t608 - t694) * t888 + t1000;
t299 = -0.4e1 * t488 * t642 + 0.4e1 * t489 * t643;
t288 = -t797 * t705 - t798 * t707 + t889 * t903;
t287 = -t797 * t704 - t798 * t706 + t889 * t904;
t286 = t705 * t795 - t707 * t796 + t888 * t903;
t285 = t704 * t795 - t706 * t796 + t888 * t904;
t274 = t1313 / 0.4e1;
t273 = 0.4e1 * t488 * t573 + 0.4e1 * t489 * t574;
t269 = -0.4e1 * t1122 * t489 + 0.4e1 * t488 * t571;
t268 = qJD(3) * t1066;
t245 = t1315 / 0.4e1;
t230 = -t1076 * t957 + t1386 * t344;
t210 = -t493 * t872 + t870 * t961;
t208 = t1379 * t299 + t979;
t197 = 0.4e1 * t344 * t529 + 0.4e1 * t785 * t957;
t195 = t1319 / 0.4e1;
t190 = 0.4e1 * t923;
t171 = t1320 / 0.4e1;
t165 = -t287 * t889 + t288 * t888;
t164 = -t285 * t889 + t286 * t888;
t162 = -t1076 * t960 + t1386 * t246;
t157 = t388 + t245 - t1302 / 0.4e1;
t156 = t502 + t388 - t1315 / 0.4e1;
t155 = t502 + t245 - t1308 / 0.4e1;
t146 = t1323 / 0.4e1;
t137 = -0.4e1 * t1424 * t246 + 0.4e1 * t762 * t960;
t133 = 0.4e1 * t1476 * t436 + 0.4e1 * t401 * t483 - 0.4e1 * t437 * t535;
t130 = -t387 * t872 + (t297 * t888 + t298 * t889) * t870;
t129 = -t386 * t872 + (t295 * t888 + t296 * t889) * t870;
t128 = t379 + t195 - t1303 / 0.4e1;
t127 = t486 + t379 - t1319 / 0.4e1;
t126 = t486 + t195 - t1309 / 0.4e1;
t125 = 0.4e1 * t1477 * t424 + 0.4e1 * t378 * t473 - 0.4e1 * t425 * t520;
t124 = 0.4e1 * t246 * t451 - 0.4e1 * t468 * t683 - 0.4e1 * t470 * t684;
t122 = t1325 / 0.4e1;
t118 = t1327 / 0.4e1;
t116 = t1328 / 0.4e1;
t109 = t269 * t1379 - t1195 / 0.2e1 + t325 * t1380 + t1404 * t870 + t979;
t108 = 0.4e1 * t990;
t98 = t922 / 0.4e1 + ((-t240 - t683) * t889 + (t239 + t684) * t888) * t1463;
t96 = (-t380 + t961) * t872 + (t307 * t888 + t308 * t889 + t493) * t870;
t87 = t1379 * t162 + t1380 * t230 + t1381 * t341;
t79 = t274 + t118 - t1307 / 0.4e1;
t78 = t429 + t274 - t1327 / 0.4e1;
t77 = t429 + t118 - t1313 / 0.4e1;
t76 = 0.4e1 * t1152 * t372;
t75 = t1006 + t1052 - t1005;
t74 = t1005 - t1442;
t73 = t1005 + t1442;
t71 = t1370 / 0.4e1;
t70 = (-t348 + t963) * t872 + (t287 * t888 + t288 * t889 + t466) * t870;
t69 = (-t347 + t964) * t872 + (t285 * t888 + t286 * t889 + t464) * t870;
t48 = t273 * t1379 + t313 * t1380 + t346 * t1381 + t891 * t870 - t1406;
t46 = t1372 / 0.4e1;
t43 = 0.4e1 * t209 * t319 + 0.4e1 * t239 * t371 - 0.4e1 * t240 * t372;
t38 = t1163 + t1167 - t1147;
t37 = t1147 - t1440;
t36 = t1147 + t1440;
t35 = t137 * t1379 + t1165;
t34 = t1164 + t1168 - t1148;
t33 = t1148 - t1439;
t32 = t1148 + t1439;
t27 = t1379 * t190 + t1057;
t26 = t27 * qJD(6);
t25 = t1054 + t1056 - t1053;
t24 = t1053 + t1441;
t23 = t1053 - t1441;
t22 = t108 * t1379 + t1057;
t21 = t124 * t1379 + t1380 * t197 + t1165 + t929;
t20 = t1166 + t1271 - t1154;
t19 = t1154 - t1438;
t18 = t1154 + t1438;
t16 = t125 * t1379 + t942;
t15 = t146 - t1328 / 0.4e1 + t943;
t14 = t116 - t1323 / 0.4e1 + t943;
t13 = t116 + t146 + t1375 / 0.4e1 + t900;
t12 = t1059 * t1205 - t76 * t1379 + t989;
t11 = t71 - t1372 / 0.4e1 + t926;
t10 = t46 - t1370 / 0.4e1 + t926;
t9 = t888 * t898 + t889 * t899;
t8 = t892 + t122 + t171;
t7 = t894 + t122 - t1320 / 0.4e1;
t6 = t893 + t171 - t1325 / 0.4e1;
t5 = t43 * t1379 + t133 * t1380 + (-t96 / 0.2e1 + t1175 / 0.2e1 + t1185 / 0.2e1) * t872 + (t210 / 0.2e1 + t70 * t1356 + t69 * t1360) * t870 + t942;
t4 = t907 + t71 + t46 - t1377 / 0.4e1;
t3 = t892 + t1399 + t1403 + t1434;
t2 = t894 + (-t387 / 0.4e1 - t324 / 0.4e1) * t888 + (t386 / 0.4e1 + t323 / 0.4e1) * t889 + t1400 + t1399 - t1436;
t1 = t893 + (-t493 / 0.2e1 + (-t348 / 0.4e1 - t308 / 0.4e1) * t889 + (-t347 / 0.4e1 - t307 / 0.4e1) * t888) * t870 + (t380 / 0.2e1 + (-t466 / 0.4e1 - t417 / 0.4e1) * t889 + (-t464 / 0.4e1 + t416 / 0.4e1) * t888) * t872 + t1400 + t1403 + t1435;
t17 = [t148 * qJD(2) + t48 * qJD(3) + t174 * qJD(4) + t109 * qJD(5) + t208 * qJD(6), qJD(3) * t73 + qJD(4) * t500 + qJD(5) * t36 + qJD(6) * t156 + t1092, t48 * qJD(1) + t73 * qJD(2) + t24 * qJD(4) + t3 * qJD(5) + t8 * qJD(6) + ((-t307 / 0.2e1 - t347 / 0.2e1 - t335 / 0.2e1 - t293 / 0.2e1 + t817 * t1366 - t753 * t818 / 0.2e1 + t824 * t1356 - t899) * qJD(3) + (t1225 / 0.2e1 - t1223 / 0.2e1 + t714 / 0.2e1 - t775 / 0.2e1 - t810 / 0.2e1) * t1087 + (t716 * t1362 + t718 * t1361 + t1367 + t773 / 0.2e1 - t812 / 0.2e1) * t1088) * t889 + ((t294 / 0.2e1 + t348 / 0.2e1 + t336 / 0.2e1 + t819 * t1366 + t820 * t753 / 0.2e1 + t308 / 0.2e1 + t824 * t1360 - t898) * qJD(3) + (-t1224 / 0.2e1 + t1222 / 0.2e1 - t715 / 0.2e1 + t776 / 0.2e1 - t811 / 0.2e1) * t1087 + (-t717 * t1362 - t719 * t1361 + t650 / 0.2e1 - t774 / 0.2e1 + t813 / 0.2e1) * t1088) * t888 + (-t468 * t574 + t469 * t489 - t470 * t573 + t471 * t488) * t1291 + (t527 * t550 + t528 * t548 - t547 * t593 - t549 * t592) * t1294 + (t561 * t674 + t562 * t672 - t644 * t673 - t645 * t671) * t1296 + (-t948 * t834 + (-t814 * t889 + t816 * t888) * t832) * qJD(3) * m(4), qJD(2) * t500 + qJD(3) * t24 + qJD(5) * t32 + qJD(6) * t127 + t1091, t109 * qJD(1) + t36 * qJD(2) + t3 * qJD(3) + t32 * qJD(4) + t920 * qJD(5) + t13 * qJD(6) + (-t418 - t439) * t1085 + (t76 / 0.4e1 + t371 * t571 + t372 * t1122 + t475 * t488 + t476 * t489) * t1290 + (-t1476 * t669 + t527 * t557 + t528 * t558 - t535 * t670) * t1293 + ((t324 / 0.2e1 + t387 / 0.2e1) * t889 + (t386 / 0.2e1 + t1368 - t1059) * t888) * t1086, t208 * qJD(1) + t156 * qJD(2) + t8 * qJD(3) + t127 * qJD(4) + t13 * qJD(5) + t900 * qJD(6) + (-t1477 * t642 - t520 * t643 + t405 / 0.2e1 + t406 / 0.2e1) * t1289; t74 * qJD(3) - t501 * qJD(4) + t37 * qJD(5) + t155 * qJD(6) - t1092, 0, t74 * qJD(1) + ((-m(5) * t672 - m(6) * t548 - m(7) * t469) * t889 + (m(5) * t674 + m(6) * t550 + m(7) * t471) * t888) * qJD(3) + t98 * qJD(5) + qJD(6) * t1066, -t1083, t37 * qJD(1) + t98 * qJD(3) + (t888 * t993 - t889 * t992) * qJD(5) + t419, t155 * qJD(1) + t1061 * t426 + t268 + t419; t75 * qJD(2) + t9 * qJD(3) + t25 * qJD(4) + t1 * qJD(5) + t6 * qJD(6) + t273 * t1481 + t313 * t1480 - t891 * t1090 + ((t1367 + t648 / 0.2e1) * t1192 - m(5) * t346 / 0.4e1 + t1406) * qJD(1), t75 * qJD(1) + t99 * qJD(5) - t272 * t1289 / 0.4e1, t9 * qJD(1) + t87 * qJD(4) + t21 * qJD(5) + t35 * qJD(6) + ((t1093 * t832 * t834 - (t888 * (rSges(4,1) * t1192 - t1095) + t889 * (rSges(4,1) * t1191 + t1015)) * t944) * m(4) + (t246 * t302 - t468 * t469 - t470 * t471) * m(7) + (t344 * t403 - t547 * t548 - t549 * t550) * m(6) + (t461 * t492 - t671 * t672 - t673 * t674) * m(5) + (t882 * t809 + (-t888 * t808 + t1405 + (t1107 * t889 + t1013) * t870) * t889 + (-t819 * t717 - t820 * t719) * t888 + (t819 * t716 + t820 * t718 + t925 + (-t889 * t932 + t924) * t870) * t889 + t145 + t165) * t1360 + (t883 * t808 - (t817 * t716 - t818 * t718) * t889 + t144 + t164 + (-t889 * t809 + t1405 + t817 * t717 - t818 * t719 + t925 + (t1013 + t924 + (t1107 - t932) * t889) * t870) * t888) * t1358) * qJD(3), t25 * qJD(1) + t87 * qJD(3) + t18 * qJD(5) + t78 * qJD(6) + (-t494 + t1008 * (-0.2e1 * t822 + t1432) * t872) * qJD(4), t1 * qJD(1) + t1170 + t21 * qJD(3) + t18 * qJD(4) + (t129 * t1358 + t130 * t1360 + t907) * qJD(5) + t4 * qJD(6) + (t96 / 0.2e1 + (t1368 + t1369) * t889 + (-t324 / 0.2e1 - t193 / 0.2e1) * t888) * t1085 + (-t210 / 0.2e1 + (t169 / 0.2e1 - t70 / 0.2e1) * t889 + (t168 / 0.2e1 - t69 / 0.2e1) * t888) * t1086 + (t246 * t433 + t319 * t451 + t371 * t684 - t372 * t683 - t468 * t476 - t470 * t475 - t43 / 0.4e1) * t1290 + (t509 * t344 + t483 * t529 - t558 * t547 - t557 * t549 - t133 / 0.4e1 - t1415 * t785) * t1293, t6 * qJD(1) + t35 * qJD(3) + t78 * qJD(4) + t4 * qJD(5) + t907 * qJD(6) + (-t1084 / 0.4e1 + (-t473 * t1424 + t203 / 0.2e1 - t383 / 0.2e1 + t384 / 0.2e1 - t125 / 0.4e1 - t1414 * t762) * qJD(6)) * m(7); t501 * qJD(2) + t23 * qJD(3) + t33 * qJD(5) + t126 * qJD(6) - t1091, t1083, t23 * qJD(1) + t494 * qJD(4) + t19 * qJD(5) + t77 * qJD(6) + (-t162 / 0.4e1 + (-t302 - t960) * t872 + (t469 * t888 + t471 * t889 + t246) * t870) * t1291 + (-t230 / 0.4e1 + (-t403 - t957) * t872 + (t548 * t888 + t550 * t889 + t344) * t870) * t1294 + (-t341 / 0.4e1 + (-t492 + t949) * t872 + (t672 * t888 + t674 * t889 + t461) * t870) * t1296, t494 * qJD(3), t33 * qJD(1) + t19 * qJD(3) + ((-m(6) * t509 - m(7) * t433) * t872 + (t888 * t992 + t889 * t993) * t870) * qJD(5) + t316, t126 * qJD(1) + t77 * qJD(3) + t1061 * t317 + t316; (t1195 / 0.2e1 + t919) * qJD(1) + t38 * qJD(2) + t2 * qJD(3) + t34 * qJD(4) + t12 * qJD(5) + t14 * qJD(6) + (-t269 / 0.4e1 + t1152 * t489) * t1292 + t325 * t1480 - t1404 * t1090, qJD(1) * t38 - qJD(3) * t99, t2 * qJD(1) - t1170 + (t918 + t69 * t1358 + (-t307 * t889 + t308 * t888) * t1363 + (t416 * t889 + t417 * t888) * t1364 + t70 * t1360 + t1473 * t1032 + t165 * t1038 + (-t389 * t889 + t390 * t888) * t1035 + t164 * t1042 - t929) * qJD(3) + t20 * qJD(4) + t5 * qJD(5) + t10 * qJD(6) + (-t124 / 0.4e1 + t209 * t246 - t239 * t470 - t240 * t468 + t302 * t319 + t371 * t471 - t372 * t469) * t1291 + (t344 * t401 + t403 * t483 - t436 * t549 - t437 * t547 + t1476 * t550 - t535 * t548 - t197 / 0.4e1) * t1294, qJD(1) * t34 + qJD(3) * t20, t12 * qJD(1) + t5 * qJD(3) + ((t319 * t433 + t371 * t475 - t372 * t476) * m(7) + (t1476 * t557 + t483 * t509 - t535 * t558) * m(6) + t872 ^ 2 * t439 / 0.2e1 + (t130 * t1356 + t129 * t1360 + (t323 * t888 + t324 * t889) * t1363) * t870 + t1057) * qJD(5) + t22 * qJD(6), t14 * qJD(1) + t10 * qJD(3) + t22 * qJD(5) + t1057 * qJD(6) + (-t190 / 0.4e1 + t923 + t990) * t1289; t919 * qJD(1) + t157 * qJD(2) + t7 * qJD(3) + t128 * qJD(4) + t15 * qJD(5) + qJD(6) * t989 + t1481 * t299, qJD(1) * t157 + t268, t7 * qJD(1) + t918 * qJD(3) + t79 * qJD(4) + t11 * qJD(5) + t16 * qJD(6) + (t1084 / 0.4e1 + (t246 * t378 + t302 * t473 - t424 * t470 - t425 * t468 - t469 * t520 + t471 * t1477 - t137 / 0.4e1) * qJD(3)) * m(7), qJD(1) * t128 + qJD(3) * t79, t15 * qJD(1) + t11 * qJD(3) + t1057 * qJD(5) + t26 + (t433 * t473 + t475 * t1477 - t476 * t520 - t108 / 0.4e1 + t990) * t1290, qJD(1) * t989 + qJD(3) * t16 + qJD(5) * t27 + t26;];
Cq  = t17;
