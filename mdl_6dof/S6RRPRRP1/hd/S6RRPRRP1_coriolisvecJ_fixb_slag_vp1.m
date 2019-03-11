% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:39:00
% EndTime: 2019-03-09 11:41:02
% DurationCPUTime: 105.39s
% Computational Cost: add. (60891->1421), mult. (61500->1809), div. (0->0), fcn. (55983->10), ass. (0->725)
t1411 = -Icges(6,4) - Icges(7,4);
t1349 = Icges(6,1) + Icges(7,1);
t1374 = Icges(6,5) + Icges(7,5);
t1373 = -Icges(6,2) - Icges(7,2);
t1372 = Icges(6,6) + Icges(7,6);
t1410 = Icges(3,3) + Icges(4,3);
t707 = cos(qJ(5));
t1409 = t1411 * t707;
t704 = sin(qJ(5));
t1408 = t1411 * t704;
t699 = qJ(2) + pkin(10);
t670 = sin(t699);
t671 = cos(t699);
t705 = sin(qJ(2));
t708 = cos(qJ(2));
t1380 = Icges(3,5) * t708 + Icges(4,5) * t671 - Icges(3,6) * t705 - Icges(4,6) * t670;
t1407 = -t1372 * t704 + t1374 * t707;
t1406 = t1373 * t704 - t1409;
t1405 = t1349 * t707 + t1408;
t1403 = Icges(6,3) + Icges(7,3);
t709 = cos(qJ(1));
t1404 = t1410 * t709;
t672 = qJ(4) + t699;
t654 = sin(t672);
t655 = cos(t672);
t1382 = -t1372 * t655 + t1406 * t654;
t1358 = -t1374 * t655 + t1405 * t654;
t706 = sin(qJ(1));
t1083 = t706 * t708;
t1086 = t705 * t706;
t1090 = t671 * t706;
t1092 = t670 * t706;
t1368 = -Icges(3,5) * t1083 - Icges(4,5) * t1090 + Icges(3,6) * t1086 + Icges(4,6) * t1092 + t1404;
t1381 = t1380 * t709 + t1410 * t706;
t1394 = Icges(3,5) * t705 + Icges(4,5) * t670 + Icges(3,6) * t708 + Icges(4,6) * t671;
t1402 = t1407 * t655;
t1401 = t1406 * t655;
t1400 = t1405 * t655;
t1399 = -t1372 * t707 - t1374 * t704;
t1398 = t1373 * t707 + t1408;
t1397 = -t1349 * t704 + t1409;
t1145 = Icges(4,6) * t709;
t470 = Icges(4,4) * t1090 - Icges(4,2) * t1092 - t1145;
t1146 = Icges(3,6) * t709;
t525 = Icges(3,4) * t1083 - Icges(3,2) * t1086 - t1146;
t1396 = t470 * t670 + t525 * t705;
t1096 = t654 * t709;
t1360 = -t1403 * t655 + t1407 * t654;
t1078 = t709 * t704;
t1084 = t706 * t707;
t541 = -t1078 * t655 + t1084;
t1082 = t707 * t709;
t1087 = t704 * t706;
t542 = t1082 * t655 + t1087;
t1251 = t1096 * t1360 + t1358 * t542 + t1382 * t541;
t286 = Icges(7,5) * t542 + Icges(7,6) * t541 + Icges(7,3) * t1096;
t289 = Icges(6,5) * t542 + Icges(6,6) * t541 + Icges(6,3) * t1096;
t1361 = t286 + t289;
t513 = Icges(7,4) * t541;
t298 = Icges(7,1) * t542 + Icges(7,5) * t1096 + t513;
t516 = Icges(6,4) * t541;
t301 = Icges(6,1) * t542 + Icges(6,5) * t1096 + t516;
t1383 = t298 + t301;
t1157 = Icges(7,4) * t542;
t292 = Icges(7,2) * t541 + Icges(7,6) * t1096 + t1157;
t1160 = Icges(6,4) * t542;
t295 = Icges(6,2) * t541 + Icges(6,6) * t1096 + t1160;
t1384 = t292 + t295;
t1307 = t1096 * t1361 + t1383 * t542 + t1384 * t541;
t1097 = t654 * t706;
t539 = t1087 * t655 + t1082;
t511 = Icges(7,4) * t539;
t540 = t1084 * t655 - t1078;
t297 = -Icges(7,1) * t540 - Icges(7,5) * t1097 + t511;
t514 = Icges(6,4) * t539;
t300 = -Icges(6,1) * t540 - Icges(6,5) * t1097 + t514;
t1298 = t297 + t300;
t512 = Icges(7,4) * t540;
t290 = -Icges(7,2) * t539 + Icges(7,6) * t1097 + t512;
t515 = Icges(6,4) * t540;
t293 = -Icges(6,2) * t539 + Icges(6,6) * t1097 + t515;
t1299 = t290 + t293;
t284 = Icges(7,5) * t540 - Icges(7,6) * t539 + Icges(7,3) * t1097;
t287 = Icges(6,5) * t540 - Icges(6,6) * t539 + Icges(6,3) * t1097;
t1300 = t284 + t287;
t1308 = t1096 * t1300 - t1298 * t542 + t1299 * t541;
t1005 = qJD(5) * t709;
t676 = qJD(2) * t706;
t615 = qJD(4) * t706 + t676;
t508 = t1005 * t654 + t615;
t1006 = qJD(5) * t706;
t698 = qJD(2) + qJD(4);
t616 = t698 * t709;
t509 = -t1006 * t654 + t616;
t1007 = qJD(5) * t655;
t612 = qJD(1) - t1007;
t1316 = t1251 * t612 + t1307 * t508 - t1308 * t509;
t1252 = t1097 * t1360 + t1358 * t540 - t1382 * t539;
t1309 = t1097 * t1361 + t1383 * t540 - t1384 * t539;
t1310 = t1097 * t1300 - t1298 * t540 - t1299 * t539;
t1317 = t1252 * t612 + t1309 * t508 - t1310 * t509;
t1162 = Icges(4,4) * t670;
t575 = Icges(4,1) * t671 - t1162;
t473 = Icges(4,5) * t706 + t575 * t709;
t1163 = Icges(3,4) * t705;
t622 = Icges(3,1) * t708 - t1163;
t528 = Icges(3,5) * t706 + t622 * t709;
t1395 = -t528 * t1083 - t473 * t1090;
t572 = Icges(4,2) * t671 + t1162;
t653 = Icges(4,4) * t671;
t574 = Icges(4,1) * t670 + t653;
t619 = Icges(3,2) * t708 + t1163;
t686 = Icges(3,4) * t708;
t621 = Icges(3,1) * t705 + t686;
t1379 = t572 * t670 - t574 * t671 + t619 * t705 - t621 * t708;
t1153 = Icges(4,5) * t709;
t630 = Icges(4,4) * t1092;
t472 = Icges(4,1) * t1090 - t1153 - t630;
t1154 = Icges(3,5) * t709;
t648 = Icges(3,4) * t1086;
t527 = Icges(3,1) * t1083 - t1154 - t648;
t1364 = -t472 * t671 - t527 * t708 + t1396;
t912 = qJD(1) * t655 - qJD(5);
t989 = t654 * t616;
t1229 = t706 * t912 + t989;
t832 = t707 * t612;
t269 = t1229 * t704 + t709 * t832;
t831 = t612 * t704;
t270 = -t1229 * t707 + t709 * t831;
t1011 = qJD(1) * t706;
t971 = t654 * t1011;
t988 = t655 * t616;
t784 = -t971 + t988;
t149 = Icges(7,5) * t270 + Icges(7,6) * t269 + Icges(7,3) * t784;
t151 = Icges(6,5) * t270 + Icges(6,6) * t269 + Icges(6,3) * t784;
t1393 = t149 + t151;
t1088 = t698 * t706;
t990 = t654 * t1088;
t271 = t706 * t832 + (-t709 * t912 + t990) * t704;
t1099 = t654 * t698;
t272 = t912 * t1082 + (-t1099 * t707 + t831) * t706;
t1010 = qJD(1) * t709;
t785 = t1010 * t654 + t1088 * t655;
t150 = Icges(7,5) * t272 + Icges(7,6) * t271 + Icges(7,3) * t785;
t152 = Icges(6,5) * t272 + Icges(6,6) * t271 + Icges(6,3) * t785;
t1392 = t150 + t152;
t153 = Icges(7,4) * t270 + Icges(7,2) * t269 + Icges(7,6) * t784;
t155 = Icges(6,4) * t270 + Icges(6,2) * t269 + Icges(6,6) * t784;
t1391 = t153 + t155;
t154 = Icges(7,4) * t272 + Icges(7,2) * t271 + Icges(7,6) * t785;
t156 = Icges(6,4) * t272 + Icges(6,2) * t271 + Icges(6,6) * t785;
t1390 = t154 + t156;
t157 = Icges(7,1) * t270 + Icges(7,4) * t269 + Icges(7,5) * t784;
t159 = Icges(6,1) * t270 + Icges(6,4) * t269 + Icges(6,5) * t784;
t1389 = t157 + t159;
t158 = Icges(7,1) * t272 + Icges(7,4) * t271 + Icges(7,5) * t785;
t160 = Icges(6,1) * t272 + Icges(6,4) * t271 + Icges(6,5) * t785;
t1388 = t158 + t160;
t1387 = t1402 * t698 + (qJD(5) * t1399 + t1403 * t698) * t654;
t1386 = t1401 * t698 + (qJD(5) * t1398 + t1372 * t698) * t654;
t1385 = t1400 * t698 + (t1397 * qJD(5) + t1374 * t698) * t654;
t1378 = t1358 * t707 - t1382 * t704;
t1377 = t1381 * t709 + t1395;
t1081 = t708 * t709;
t1089 = t671 * t709;
t1376 = -t527 * t1081 - t472 * t1089 + t1368 * t706;
t1266 = t528 * t1081 + t473 * t1089 + t1381 * t706;
t1265 = t1394 * t709;
t1264 = t1394 * t706;
t874 = -Icges(4,2) * t670 + t653;
t471 = Icges(4,6) * t706 + t709 * t874;
t875 = -Icges(3,2) * t705 + t686;
t526 = Icges(3,6) * t706 + t709 * t875;
t1375 = t471 * t670 + t526 * t705;
t1371 = -t1364 * t706 + t1368 * t709;
t1348 = -t1086 * t526 - t1092 * t471 - t1377;
t1085 = t705 * t709;
t1091 = t670 * t709;
t1347 = -t1085 * t525 - t1091 * t470 - t1376;
t1346 = -t1085 * t526 - t1091 * t471 + t1266;
t1370 = t1379 * t706 + t1265;
t1369 = -t1379 * t709 + t1264;
t1302 = t470 * t671 + t472 * t670 + t525 * t708 + t527 * t705;
t1301 = t471 * t671 + t473 * t670 + t526 * t708 + t528 * t705;
t1366 = t1394 * qJD(2);
t1365 = t473 * t671 + t528 * t708 - t1375;
t1323 = t1096 * t1392 - t1298 * t270 + t1299 * t269 + t1300 * t784 + t1388 * t542 + t1390 * t541;
t1322 = t1096 * t1393 + t1361 * t784 + t1383 * t270 + t1384 * t269 + t1389 * t542 + t1391 * t541;
t1321 = t1097 * t1392 - t1298 * t272 + t1299 * t271 + t1300 * t785 + t1388 * t540 - t1390 * t539;
t1320 = t1097 * t1393 + t1361 * t785 + t1383 * t272 + t1384 * t271 + t1389 * t540 - t1391 * t539;
t1314 = t1096 * t1387 + t1358 * t270 + t1360 * t784 + t1382 * t269 + t1385 * t542 + t1386 * t541;
t1313 = t1097 * t1387 + t1358 * t272 + t1360 * t785 + t1382 * t271 + t1385 * t540 - t1386 * t539;
t850 = -t290 * t704 - t297 * t707;
t127 = -t284 * t655 + t654 * t850;
t848 = -t293 * t704 - t300 * t707;
t129 = -t287 * t655 + t654 * t848;
t1363 = t127 + t129;
t849 = -t292 * t704 + t298 * t707;
t128 = -t286 * t655 + t654 * t849;
t847 = -t295 * t704 + t301 * t707;
t130 = -t289 * t655 + t654 * t847;
t1362 = t128 + t130;
t1250 = -t1360 * t655 + t1378 * t654;
t702 = -qJ(6) - pkin(9);
t1185 = pkin(9) + t702;
t664 = pkin(5) * t707 + pkin(4);
t1186 = pkin(4) - t664;
t885 = rSges(7,1) * t707 - rSges(7,2) * t704;
t1062 = (-rSges(7,3) + t1185) * t655 + (-t1186 + t885) * t654;
t1094 = t655 * t706;
t1357 = -rSges(7,1) * t540 + rSges(7,2) * t539 - t664 * t1094;
t545 = t874 * qJD(2);
t546 = t575 * qJD(2);
t588 = t875 * qJD(2);
t589 = t622 * qJD(2);
t1356 = -t545 * t670 + t546 * t671 - t588 * t705 + t589 * t708 + (-t572 * t671 - t574 * t670 - t619 * t708 - t621 * t705) * qJD(2) + t1394 * qJD(1);
t1002 = qJD(6) * t709;
t599 = t654 * t1002;
t1176 = pkin(5) * qJD(5);
t991 = t707 * t1176;
t1190 = pkin(5) * t704;
t993 = qJD(1) * t1190;
t1355 = t270 * rSges(7,1) + t269 * rSges(7,2) + rSges(7,3) * t988 + t702 * t971 + t706 * t991 + t709 * t993 + t599;
t960 = t655 * t1006;
t1354 = t1011 - t960;
t1353 = t1381 * qJD(1);
t1352 = (t1403 * t654 - t1378 + t1402) * t612 + (t1360 * t706 + t848 + t850) * t509 + (-t1360 * t709 - t847 - t849) * t508;
t1351 = t1379 * qJD(1) + qJD(2) * t1380;
t1247 = (t1373 * t540 - t1298 - t511 - t514) * t509 - (t1373 * t542 + t1383 + t513 + t516) * t508 - (t1398 * t654 + t1358) * t612;
t1350 = (t1378 * t698 - t1387) * t655 + (t1385 * t707 - t1386 * t704 + t1360 * t698 + (-t1358 * t704 - t1382 * t707) * qJD(5)) * t654;
t583 = pkin(9) * t988;
t992 = t704 * t1176;
t1077 = -rSges(7,3) * t971 - t583 + (pkin(9) * t1011 + t1186 * t616) * t654 + ((-t698 * t702 - t992) * t709 + t1186 * t1011) * t655 + t1355;
t1268 = t993 - (t664 * t698 - qJD(6)) * t654;
t1238 = t1185 * t654;
t381 = -t1186 * t655 - t1238;
t582 = pkin(4) * t990;
t887 = rSges(7,1) * t272 + rSges(7,2) * t271;
t1076 = -rSges(7,3) * t785 - t887 - t582 - (qJD(1) * t381 - t991) * t709 - ((-t1185 * t698 - t992) * t655 + t1268) * t706;
t1345 = t1369 * qJD(1);
t1191 = pkin(4) * t655;
t1073 = pkin(5) * t1078 + (t1191 + t1238) * t706 - rSges(7,3) * t1097 + t1357;
t1344 = t885 * t655 + t381;
t797 = qJD(2) * t572;
t314 = -t709 * t797 + (-t706 * t874 + t1145) * qJD(1);
t799 = qJD(2) * t574;
t316 = -t709 * t799 + (-t575 * t706 + t1153) * qJD(1);
t798 = qJD(2) * t619;
t367 = -t709 * t798 + (-t706 * t875 + t1146) * qJD(1);
t800 = qJD(2) * t621;
t369 = -t709 * t800 + (-t622 * t706 + t1154) * qJD(1);
t1343 = -qJD(2) * t1301 - t314 * t670 + t316 * t671 - t367 * t705 + t369 * t708 + t1353;
t315 = qJD(1) * t471 - t706 * t797;
t317 = qJD(1) * t473 - t706 * t799;
t368 = qJD(1) * t526 - t706 * t798;
t370 = qJD(1) * t528 - t706 * t800;
t1342 = qJD(1) * t1368 + qJD(2) * t1302 + t315 * t670 - t317 * t671 + t368 * t705 - t370 * t708;
t1341 = (-t1372 * t540 - t1374 * t539) * t509 + (t1372 * t542 - t1374 * t541) * t508 - t1399 * t612 * t654;
t1340 = (t1346 * t706 - t1347 * t709) * qJD(2);
t1339 = (t1348 * t706 - t1371 * t709) * qJD(2);
t1338 = t1370 * qJD(1);
t1093 = t655 * t709;
t1337 = t542 * rSges(7,1) + t541 * rSges(7,2) + rSges(7,3) * t1096 + pkin(5) * t1087 + t664 * t1093;
t1333 = qJD(1) * t1364 - t1366 * t706 + t1353;
t1332 = -t1366 * t709 + (-t1380 * t706 - t1365 + t1404) * qJD(1);
t1329 = t1073 * t612 + t599;
t1003 = qJD(6) * t706;
t1008 = qJD(5) * t654;
t614 = pkin(4) * t1093;
t507 = pkin(9) * t1096 + t614;
t1072 = -t1096 * t702 + t1337 - t507;
t1328 = t655 * t1003 + t1008 * t1072;
t1327 = 0.2e1 * qJD(2);
t998 = qJD(1) * qJD(2);
t660 = t706 * t998;
t996 = qJD(1) * qJD(4);
t585 = t706 * t996 + t660;
t376 = qJD(5) * t785 + t585;
t661 = t709 * t998;
t586 = t709 * t996 + t661;
t377 = qJD(5) * t784 + t586;
t961 = t698 * t1008;
t1326 = t1251 * t961 + t1307 * t377 + t1308 * t376 + t1314 * t612 + t1322 * t508 - t1323 * t509;
t1325 = t1252 * t961 + t1309 * t377 + t1310 * t376 + t1313 * t612 + t1320 * t508 - t1321 * t509;
t1324 = rSges(7,1) + pkin(5);
t37 = (t698 * t850 - t150) * t655 + (-t154 * t704 + t158 * t707 + t284 * t698 + (-t290 * t707 + t297 * t704) * qJD(5)) * t654;
t39 = (t698 * t848 - t152) * t655 + (-t156 * t704 + t160 * t707 + t287 * t698 + (-t293 * t707 + t300 * t704) * qJD(5)) * t654;
t1319 = t37 + t39;
t38 = (t698 * t849 - t149) * t655 + (-t153 * t704 + t157 * t707 + t286 * t698 + (-t292 * t707 - t298 * t704) * qJD(5)) * t654;
t40 = (t698 * t847 - t151) * t655 + (-t155 * t704 + t159 * t707 + t289 * t698 + (-t295 * t707 - t301 * t704) * qJD(5)) * t654;
t1318 = t38 + t40;
t1315 = t1250 * t612 + t1362 * t508 - t1363 * t509;
t1312 = -t1338 + t1339;
t1311 = t1340 + t1345;
t1306 = t1351 * t706 + t1356 * t709;
t1305 = -t1351 * t709 + t1356 * t706;
t1304 = qJD(2) * t1364 - t315 * t671 - t317 * t670 - t368 * t708 - t370 * t705;
t1303 = qJD(2) * t1365 + t314 * t671 + t316 * t670 + t367 * t708 + t369 * t705;
t1297 = t1062 * t706;
t1296 = t1062 * t709;
t1178 = rSges(7,3) * t654;
t1295 = -t1178 - t1344;
t1294 = t1382 * t706;
t1293 = t1382 * t709;
t1292 = t1358 * t706;
t1291 = t1358 * t709;
t890 = rSges(6,1) * t540 - rSges(6,2) * t539;
t303 = rSges(6,3) * t1097 + t890;
t889 = rSges(6,1) * t707 - rSges(6,2) * t704;
t430 = -rSges(6,3) * t655 + t654 * t889;
t395 = t430 * t706;
t1180 = rSges(6,3) * t654;
t811 = t889 * t655;
t432 = t811 + t1180;
t556 = pkin(4) * t654 - pkin(9) * t655;
t504 = t556 * t706;
t557 = pkin(9) * t654 + t1191;
t916 = qJD(1) * t504 - t616 * t557;
t1290 = t1008 * t303 + t1354 * t430 - t395 * t612 + t509 * t432 - t916;
t1289 = t1372 * t654 + t1401;
t1288 = t1374 * t654 + t1400;
t1144 = Icges(5,6) * t709;
t452 = Icges(5,4) * t1094 - Icges(5,2) * t1097 - t1144;
t644 = Icges(5,4) * t655;
t552 = Icges(5,1) * t654 + t644;
t1287 = -t552 * t706 - t452;
t873 = -Icges(5,2) * t654 + t644;
t453 = Icges(5,6) * t706 + t709 * t873;
t1286 = -t552 * t709 - t453;
t1161 = Icges(5,4) * t654;
t553 = Icges(5,1) * t655 - t1161;
t455 = Icges(5,5) * t706 + t553 * t709;
t550 = Icges(5,2) * t655 + t1161;
t1285 = -t550 * t709 + t455;
t1284 = -t550 + t553;
t1283 = t552 + t873;
t1282 = t1352 * t654;
t1218 = (-t572 * t709 + t473) * t706 - (-Icges(4,2) * t1090 + t472 - t630) * t709;
t1219 = (-t619 * t709 + t528) * t706 - (-Icges(3,2) * t1083 + t527 - t648) * t709;
t769 = t525 * t709 - t526 * t706;
t770 = t470 * t709 - t471 * t706;
t1281 = -t1218 * t670 - t1219 * t705 + t770 * t671 + t769 * t708;
t1025 = t621 + t875;
t1026 = -t619 + t622;
t1031 = t574 + t874;
t1032 = -t572 + t575;
t1280 = (-t1025 * t705 + t1026 * t708 - t1031 * t670 + t1032 * t671) * qJD(1);
t1279 = -t1300 * t509 + t1360 * t612 + t1361 * t508;
t1244 = t1362 * t709 + t1363 * t706;
t1278 = -t1362 * t706 + t1363 * t709;
t1243 = t1307 * t709 + t1308 * t706;
t1277 = -t1307 * t706 + t1308 * t709;
t1242 = t1309 * t709 + t1310 * t706;
t1276 = -t1309 * t706 + t1310 * t709;
t1275 = t1368 + t1375;
t969 = t705 * t1011;
t641 = pkin(2) * t969;
t970 = t670 * t1011;
t1027 = pkin(3) * t970 + t641;
t510 = t556 * t1011;
t1194 = pkin(2) * t705;
t1193 = pkin(3) * t670;
t596 = -t1193 - t1194;
t941 = t596 + t1194;
t486 = t941 * t706;
t915 = -qJD(1) * t486 + t641;
t1274 = t510 + t1027 - t915;
t640 = t676 * t1194;
t675 = qJD(3) * t709;
t1020 = t640 + t675;
t606 = t676 * t1193;
t973 = -t606 - t1020;
t693 = t706 * pkin(7);
t634 = t709 * pkin(1) + t693;
t694 = t708 * pkin(2);
t665 = t694 + pkin(1);
t636 = t709 * t665;
t703 = -qJ(3) - pkin(7);
t918 = -t703 * t706 + t636;
t466 = t918 - t634;
t1042 = t466 + t634;
t697 = -pkin(8) + t703;
t1018 = -t697 + t703;
t1192 = pkin(3) * t671;
t590 = t665 + t1192;
t559 = t709 * t590;
t379 = t1018 * t706 + t559 - t636;
t976 = t379 + t1042;
t752 = (t507 + t976) * qJD(1) - t615 * t556 + t973;
t958 = t654 * t1003;
t79 = -t1062 * t508 + t1072 * t612 + t752 + t958;
t1135 = qJD(1) * t79;
t1004 = qJD(6) * t655;
t884 = -rSges(7,1) * t704 - rSges(7,2) * t707;
t1075 = -t1004 + t1344 * t698 + (rSges(7,3) * t698 + qJD(5) * t884 - t992) * t654;
t600 = t655 * t1002;
t483 = t557 * t698;
t997 = qJD(1) * qJD(3);
t1024 = qJD(1) * t640 + t709 * t997;
t710 = qJD(2) ^ 2;
t822 = (-t694 - t1192) * t710;
t751 = qJD(1) * t606 + t709 * t822 + t1024;
t740 = -t616 * t483 + t585 * t556 + t751;
t309 = pkin(9) * t785 + qJD(1) * t614 - t582;
t1187 = pkin(1) - t665;
t645 = t703 * t1011;
t972 = t645 + t1020;
t346 = (-t1187 * t709 - t693) * qJD(1) - t972;
t593 = t634 * qJD(1);
t1066 = -t346 - t593;
t1030 = t590 - t665;
t576 = t596 * qJD(2);
t638 = t697 * t1011;
t919 = t576 * t706 - t638;
t219 = t1010 * t1030 + t640 + t645 + t919;
t986 = -t219 + t1066;
t910 = -t309 + t986;
t937 = t654 * t1073;
t26 = t1076 * t612 - t1075 * t509 + t1062 * t376 + (qJD(5) * t937 + t600) * t698 + (t910 - t958) * qJD(1) + t740;
t1273 = (t1135 + t26) * t709;
t307 = t542 * rSges(6,1) + t541 * rSges(6,2) + rSges(6,3) * t1096;
t105 = t307 * t612 - t430 * t508 + t752;
t891 = rSges(6,1) * t272 + rSges(6,2) * t271;
t164 = rSges(6,3) * t785 + t891;
t888 = -rSges(6,1) * t704 - rSges(6,2) * t707;
t231 = t698 * t811 + (rSges(6,3) * t698 + qJD(5) * t888) * t654;
t42 = qJD(1) * t910 - t164 * t612 - t231 * t509 - t303 * t961 + t376 * t430 + t740;
t1272 = (qJD(1) * t105 + t42) * t709;
t1009 = qJD(2) * t709;
t662 = t709 * t703;
t1023 = t706 * t665 + t662;
t695 = t709 * pkin(7);
t633 = pkin(1) * t706 - t695;
t465 = t633 - t1023;
t1052 = t466 * t1009 - t465 * t676;
t478 = rSges(4,1) * t1090 - rSges(4,2) * t1092 - t709 * rSges(4,3);
t688 = t706 * rSges(4,3);
t479 = rSges(4,1) * t1089 - rSges(4,2) * t1091 + t688;
t838 = t478 * t706 + t479 * t709;
t177 = qJD(2) * t838 + t1052;
t1051 = -t706 * t465 + t709 * t466;
t788 = t838 + t1051;
t1271 = qJD(2) * t788 + t177;
t1270 = t1062 * t1354 - t1295 * t509 - t1297 * t612 - t600 - t916;
t505 = t557 * t706;
t1034 = -t706 * t590 - t709 * t697;
t378 = t1023 + t1034;
t908 = t379 * t1009 - t378 * t676 + t1052;
t781 = t615 * t505 + t507 * t616 + t908;
t61 = t1072 * t509 - t1073 * t508 - t1004 + t781;
t939 = t79 * t1062;
t1269 = t1073 * t61 + t939;
t1267 = t1380 * qJD(1);
t1263 = -t303 * t612 - t509 * t430;
t1260 = -rSges(6,3) - pkin(9);
t786 = -t1011 * t655 - t989;
t308 = pkin(4) * t786 - pkin(9) * t971 + t583;
t1230 = qJD(1) * t308 - t615 * t483 - t586 * t556;
t643 = qJD(6) * t654;
t1079 = t709 * t576;
t966 = t705 * t1009;
t913 = pkin(2) * t966;
t218 = t913 + t1079 + (t1018 * t709 - t1030 * t706) * qJD(1);
t669 = pkin(7) * t1010;
t674 = qJD(3) * t706;
t345 = -t913 - t669 + t674 + (t1187 * t706 - t662) * qJD(1);
t581 = qJD(1) * (-pkin(1) * t1011 + t669);
t979 = qJD(1) * t345 + t706 * t997 + t581;
t911 = qJD(1) * t218 + t979;
t750 = t706 * t822 + t911;
t25 = t1077 * t612 - t1075 * t508 - t1062 * t377 + t1328 * t698 + (t576 + t643) * t1010 + t750 + t1230;
t1259 = t25 * t706;
t597 = qJD(1) * t633;
t1248 = qJD(1) * t465 - t597;
t1246 = (t1397 * t654 - t1382) * t612 + (t1349 * t539 + t1299 + t512 + t515) * t509 + (t1349 * t541 - t1157 - t1160 - t1384) * t508;
t1245 = t1341 * t654;
t1237 = t1079 + t674;
t1232 = t1250 * t961 + t1350 * t612;
t983 = t270 * rSges(6,1) + t269 * rSges(6,2) + rSges(6,3) * t988;
t162 = -rSges(6,3) * t971 + t983;
t1231 = t303 * t1010 + t709 * t162 + t706 * t164;
t548 = Icges(5,5) * t654 + Icges(5,6) * t655;
t1107 = t548 * t709;
t1116 = t453 * t654;
t1139 = Icges(5,3) * t709;
t549 = Icges(5,5) * t655 - Icges(5,6) * t654;
t1228 = -t698 * t1107 + (-t455 * t655 - t549 * t706 + t1116 + t1139) * qJD(1);
t451 = Icges(5,3) * t706 + t549 * t709;
t1016 = qJD(1) * t451;
t1108 = t548 * t706;
t1152 = Icges(5,5) * t709;
t605 = Icges(5,4) * t1097;
t454 = Icges(5,1) * t1094 - t1152 - t605;
t842 = t452 * t654 - t454 * t655;
t1227 = qJD(1) * t842 - t1108 * t698 + t1016;
t835 = t550 * t654 - t552 * t655;
t1222 = qJD(1) * t835 + t549 * t698;
t1217 = -t1010 * t1073 - t1076 * t706 + t1077 * t709;
t1214 = qJD(1) * t1283 + t1285 * t615 - (-Icges(5,2) * t1094 + t454 - t605) * t616;
t1213 = -m(7) / 0.2e1;
t1212 = m(7) / 0.2e1;
t1211 = t376 / 0.2e1;
t1210 = t377 / 0.2e1;
t1209 = -t508 / 0.2e1;
t1208 = t508 / 0.2e1;
t1207 = -t509 / 0.2e1;
t1206 = t509 / 0.2e1;
t1205 = t585 / 0.2e1;
t1204 = t586 / 0.2e1;
t1203 = -t612 / 0.2e1;
t1202 = t612 / 0.2e1;
t1201 = -t615 / 0.2e1;
t1200 = t615 / 0.2e1;
t1199 = -t616 / 0.2e1;
t1198 = t616 / 0.2e1;
t1196 = t706 / 0.2e1;
t1195 = -t709 / 0.2e1;
t1189 = -qJD(1) / 0.2e1;
t1188 = qJD(1) / 0.2e1;
t1184 = rSges(3,1) * t708;
t1183 = rSges(4,1) * t671;
t1182 = rSges(5,1) * t655;
t1181 = rSges(4,2) * t671;
t1175 = t37 * t509;
t1174 = t38 * t508;
t1173 = t39 * t509;
t1172 = t40 * t508;
t1080 = t708 * t710;
t41 = t307 * t961 + t162 * t612 - t231 * t508 - t377 * t430 + (-t1090 * t710 - t661 * t670) * pkin(3) + (-t1080 * t706 - t661 * t705) * pkin(2) + t911 + t1230;
t1171 = t41 * t709;
t1170 = t42 * t706;
t689 = t706 * rSges(3,3);
t687 = t706 * rSges(5,3);
t1169 = -rSges(7,3) + t702;
t1136 = qJD(1) * t61;
t87 = t303 * t508 + t307 * t509 + t781;
t1134 = qJD(1) * t87;
t1133 = t105 * t706;
t1128 = t127 * t376;
t1127 = t128 * t377;
t1126 = t129 * t376;
t1125 = t130 * t377;
t554 = rSges(5,1) * t654 + rSges(5,2) * t655;
t1102 = t616 * t554;
t457 = rSges(5,1) * t1094 - rSges(5,2) * t1097 - t709 * rSges(5,3);
t1044 = t465 - t633;
t978 = t378 + t1044;
t141 = -t1102 + (-t457 + t978) * qJD(1) + t1237;
t1124 = t141 * t706;
t577 = rSges(4,1) * t670 + t1181;
t943 = -t577 - t1194;
t894 = t709 * t943;
t827 = qJD(2) * t894;
t783 = t674 + t827;
t184 = (-t478 + t1044) * qJD(1) + t783;
t1121 = t184 * t577;
t1019 = rSges(3,2) * t1086 + t709 * rSges(3,3);
t537 = rSges(3,1) * t1083 - t1019;
t624 = rSges(3,1) * t705 + rSges(3,2) * t708;
t962 = t624 * t1009;
t329 = -t962 + (-t537 - t633) * qJD(1);
t1120 = t329 * t706;
t1119 = t329 * t709;
t538 = rSges(3,1) * t1081 - rSges(3,2) * t1085 + t689;
t1037 = t538 + t634;
t967 = t624 * t676;
t330 = qJD(1) * t1037 - t967;
t569 = t624 * t709;
t1118 = t330 * t569;
t450 = Icges(5,5) * t1094 - Icges(5,6) * t1097 - t1139;
t1117 = t450 * t709;
t448 = t616 * t556;
t1098 = t654 * t702;
t1095 = t655 * t698;
t1074 = -t231 - t483;
t1067 = -t307 - t507;
t1065 = -rSges(7,2) * t540 - t1324 * t539;
t1064 = -rSges(7,2) * t542 + t1324 * t541;
t1063 = t379 + t466;
t1061 = -t454 * t1093 - t706 * t450;
t1060 = t455 * t1093 + t706 * t451;
t458 = rSges(5,1) * t1093 - rSges(5,2) * t1096 + t687;
t1054 = t706 * t457 + t709 * t458;
t487 = t941 * t709;
t1050 = t487 * t1009 + t486 * t676;
t1045 = t706 * t505 + t709 * t507;
t1029 = rSges(5,2) * t971 + rSges(5,3) * t1010;
t1028 = rSges(4,2) * t970 + rSges(4,3) * t1010;
t1022 = rSges(3,2) * t969 + rSges(3,3) * t1010;
t1021 = t638 + t675;
t199 = -t706 * t835 - t1107;
t1001 = t199 * qJD(1);
t995 = pkin(2) * t1085;
t994 = qJD(2) * t694;
t987 = -t483 - t1075;
t267 = rSges(5,1) * t786 - rSges(5,2) * t988 + t1029;
t502 = t554 * t706;
t555 = -rSges(5,2) * t654 + t1182;
t268 = -t698 * t502 + (t555 * t709 + t687) * qJD(1);
t985 = t457 * t1010 + t709 * t267 + t706 * t268;
t982 = t505 * t1010 + t709 * t308 + t706 * t309;
t981 = t345 * t1009 + t346 * t676 - t465 * t661;
t980 = -t465 * t1010 + t709 * t345 + t706 * t346;
t977 = t507 + t1063;
t968 = t577 * t676;
t965 = t708 * t1009;
t959 = t655 * t1005;
t953 = t1095 / 0.2e1;
t952 = -pkin(1) - t1184;
t951 = t1011 / 0.2e1;
t950 = t1010 / 0.2e1;
t949 = -t676 / 0.2e1;
t946 = t1009 / 0.2e1;
t578 = -rSges(4,2) * t670 + t1183;
t942 = -t578 - t694;
t940 = t61 * t1072;
t936 = t705 * (-t706 ^ 2 - t709 ^ 2);
t935 = (-t706 * t873 + t1144) * qJD(1) + t1285 * t698;
t934 = qJD(1) * t453 - t1088 * t550 + t454 * t698;
t933 = (-t553 * t706 + t1152) * qJD(1) + t1286 * t698;
t932 = qJD(1) * t455 + t1287 * t698;
t398 = t455 * t1094;
t931 = t451 * t709 - t398;
t503 = t554 * t709;
t930 = -t615 * t502 - t503 * t616;
t506 = t556 * t709;
t928 = -t615 * t504 - t506 * t616;
t927 = -t450 + t1116;
t924 = t1283 * t698;
t923 = t1284 * t698;
t922 = -qJD(1) * t503 - t555 * t615;
t921 = -qJD(1) * t506 - t557 * t615;
t917 = qJD(1) * t502 - t616 * t555;
t907 = -t706 * t378 + t709 * t379 + t1051;
t906 = t706 * t303 + t709 * t307 + t1045;
t897 = (t1190 - t884) * t654;
t892 = -rSges(3,2) * t705 + t1184;
t720 = (-t505 + t978) * qJD(1) - t448 + t1237;
t78 = -t1062 * t509 + t1329 + t720;
t883 = -t706 * t79 - t709 * t78;
t104 = t1263 + t720;
t864 = t104 * t709 + t1133;
t142 = -t554 * t615 + (t458 + t976) * qJD(1) + t973;
t851 = -t141 * t709 - t142 * t706;
t846 = t303 * t709 - t307 * t706;
t845 = -t330 * t706 - t1119;
t220 = t452 * t655 + t454 * t654;
t830 = (t487 - t995) * qJD(1);
t829 = -t554 + t596;
t828 = -t556 + t596;
t826 = -t378 * t1010 + t709 * t218 + t706 * t219 + t980;
t825 = t1072 * t709 - t1073 * t706 + t1045;
t547 = t578 * qJD(2);
t824 = -pkin(2) * t1080 - qJD(2) * t547;
t823 = -t655 * t664 - t1178 - t590;
t568 = t624 * t706;
t535 = t577 * t706;
t808 = -t430 + t828;
t801 = t842 * t706;
t325 = (t537 * t706 + t538 * t709) * qJD(2);
t790 = -qJD(2) * t1192 - t994;
t789 = -t557 - t1180;
t782 = t828 - t1062;
t480 = t555 * t698;
t774 = -t480 + t790;
t773 = -t483 + t790;
t772 = qJD(1) * t549 - t1107 * t615 + t1108 * t616;
t771 = t1062 * t78 - t940;
t767 = -t231 + t773;
t766 = t826 + t982;
t755 = t773 - t1075;
t749 = qJD(1) * t378 + t1237 + t1248;
t743 = t218 * t1009 - t1063 * t660 + t219 * t676 - t378 * t661 + t981;
t741 = qJD(1) * t1284 + t1286 * t615 - t1287 * t616;
t397 = t430 * t709;
t739 = t303 * t959 - t307 * t960 - t508 * t395 - t397 * t509 + t928;
t738 = t307 * t1008 - t612 * t397 - t430 * t959 - t432 * t508 + t921;
t737 = -t1073 * t959 - t1296 * t509 - t1297 * t508 + t643 + t928;
t736 = t1295 * t508 - t1296 * t612 + t1328 + t921;
t732 = qJD(1) * t450 - t654 * t934 + t655 * t932;
t731 = -t654 * t935 + t655 * t933 + t1016;
t730 = qJD(1) * t548 - t654 * t924 + t655 * t923;
t721 = -qJD(1) * t505 - t448 + t749;
t719 = (t78 * t937 + (-t706 * t940 - t709 * t939) * t655) * qJD(5);
t718 = t616 * t308 + t615 * t309 + t586 * t505 - t585 * t507 + t743;
t713 = -t1214 * t654 + t741 * t655;
t108 = t1222 * t706 + t730 * t709;
t109 = -t1222 * t709 + t730 * t706;
t123 = t654 * t932 + t655 * t934;
t124 = t654 * t933 + t655 * t935;
t178 = -t801 - t1117;
t179 = -t1097 * t453 - t931;
t180 = -t1096 * t452 - t1061;
t181 = -t1096 * t453 + t1060;
t221 = t453 * t655 + t455 * t654;
t82 = t1227 * t706 + t732 * t709;
t83 = t1228 * t706 + t731 * t709;
t84 = -t1227 * t709 + t732 * t706;
t85 = -t1228 * t709 + t731 * t706;
t94 = -t178 * t616 + t179 * t615 + t1001;
t200 = -t709 * t835 + t1108;
t192 = t200 * qJD(1);
t95 = -t180 * t616 + t181 * t615 + t192;
t712 = (t1214 * t655 + t741 * t654) * t1189 + (-t123 * t709 + t124 * t706 + (t220 * t706 + t221 * t709) * qJD(1)) * t1188 + (t706 * t713 - t709 * t772) * t1198 + (t706 * t85 - t709 * t84 + (t178 * t706 + t179 * t709) * qJD(1)) * t1199 + (t706 * t83 - t709 * t82 + (t180 * t706 + t181 * t709) * qJD(1)) * t1200 + (t706 * t772 + t709 * t713) * t1201 + (-t180 * t709 + t181 * t706) * t1204 + (-t178 * t709 + t179 * t706) * t1205 - t1276 * t376 / 0.2e1 - t1277 * t377 / 0.2e1 + ((t1094 * t1308 + t1251 * t654) * qJD(5) + ((qJD(5) * t1307 + t1279) * t655 + t1282) * t709 + (t1288 * t542 + t1289 * t541) * t612 + (t1292 * t542 + t1294 * t541) * t509 + (-t1291 * t542 - t1293 * t541) * t508) * t1209 + (qJD(1) * t1243 + t1322 * t706 - t1323 * t709) * t1208 + (qJD(1) * t1242 + t1320 * t706 - t1321 * t709) * t1207 + ((t1093 * t1309 + t1252 * t654) * qJD(5) + ((qJD(5) * t1310 + t1279) * t655 + t1282) * t706 + (t1288 * t540 - t1289 * t539) * t612 + (t1292 * t540 - t1294 * t539) * t509 + (-t1291 * t540 + t1293 * t539) * t508) * t1206 + ((qJD(5) * t1244 - t1352) * t655 + ((t1288 * t707 - t1289 * t704 + t1360) * t612 + (t1292 * t707 - t1294 * t704 - t1300) * t509 + (-t1291 * t707 + t1293 * t704 + t1361) * t508 + t1250 * qJD(5)) * t654) * t1203 + (qJD(1) * t1244 + t1318 * t706 - t1319 * t709) * t1202 + (t94 + t1317) * t951 + (t95 + t1316) * t950 + (qJD(1) * t108 + t180 * t585 + t181 * t586 + t615 * t83 - t616 * t82 + t1326) * t1196 + (qJD(1) * t109 + t178 * t585 + t179 * t586 + t615 * t85 - t616 * t84 + t1325) * t1195 - (t1278 * t698 + t1315) * t1008 / 0.2e1 - (t1316 * t709 + t1317 * t706) * t1007 / 0.2e1;
t591 = t892 * qJD(2);
t536 = t577 * t709;
t501 = t888 * t654;
t375 = -qJD(2) * t568 + (t709 * t892 + t689) * qJD(1);
t374 = -rSges(3,2) * t965 + (-t1011 * t708 - t966) * rSges(3,1) + t1022;
t364 = rSges(6,1) * t541 - rSges(6,2) * t542;
t362 = -rSges(6,1) * t539 - rSges(6,2) * t540;
t319 = -qJD(2) * t535 + (t578 * t709 + t688) * qJD(1);
t318 = -t1009 * t1181 + (-t1009 * t670 - t1011 * t671) * rSges(4,1) + t1028;
t191 = -t591 * t1009 + (-t375 - t593 + t967) * qJD(1);
t190 = -t591 * t676 + t581 + (t374 - t962) * qJD(1);
t185 = -t968 + (t479 + t1042) * qJD(1) - t1020;
t132 = t824 * t709 + (-t319 + t968 + t1066) * qJD(1) + t1024;
t131 = t824 * t706 + (t318 + t827) * qJD(1) + t979;
t112 = t457 * t615 + t458 * t616 + t908;
t97 = -t480 * t616 + t554 * t585 + (-t268 + t986) * qJD(1) + t751;
t96 = -t480 * t615 - t554 * t586 + (t267 + t1079) * qJD(1) + t750;
t51 = t267 * t616 + t268 * t615 + t457 * t586 - t458 * t585 + t743;
t24 = t162 * t509 + t164 * t508 + t303 * t377 - t307 * t376 + t718;
t15 = -t1072 * t376 - t1073 * t377 - t1076 * t508 + t1077 * t509 + t643 * t698 + t718;
t1 = [(-(-qJD(1) * t478 + t1248 - t184 + t783) * t185 + t132 * (-t478 - t1023) + t184 * t972 + t131 * (t479 + t918) + t185 * (t674 + t1028) + (t1121 * t706 + t185 * t894) * qJD(2) + ((-t184 * rSges(4,3) + t185 * (-t665 - t1183)) * t706 + (t184 * (-t578 - t665) - t185 * t703) * t709) * qJD(1)) * m(4) + (t26 * (t1034 + t1357) + t78 * (-t887 + t1021) + t25 * (t559 + t1337) + t79 * (t1237 + t1355) + (-t25 * t1098 + t79 * (-t654 * t664 - t655 * t702) * t698 + (t26 * t704 + (-t655 * t704 * t79 + t707 * t78) * qJD(5)) * pkin(5) + (t78 * (t823 + t1098) - t79 * t697) * qJD(1)) * t709 + (t823 * t1135 + t26 * t1169 * t654 - t25 * t697 + (-t576 + (t1169 * t698 + t992) * t655 - t1268) * t78) * t706 - (t1329 + t721 - t78) * t79 + t939 * t509) * m(7) + ((t1301 + t1369) * t709 + (t1302 - t1370) * t706) * t998 / 0.2e1 + (t200 + t221) * t1204 + (-(-qJD(1) * t537 - t329 - t597 - t962) * t330 + t191 * (t706 * t952 + t1019 + t695) + t190 * t1037 + t330 * (t669 + t1022) + (t1120 * t624 - t1118) * qJD(2) + ((-pkin(1) - t892) * t1119 + (t329 * (-rSges(3,3) - pkin(7)) + t330 * t952) * t706) * qJD(1)) * m(3) + t1316 * t1206 + (t1300 * t655 + (t1298 * t707 + t1299 * t704) * t654 + t1363) * t508 * t1203 + (t199 + t220) * t1205 + (t42 * (-t890 + t1034) + t104 * (t582 - t891 + t1021) + t41 * (t559 - t1067) + t105 * (-pkin(4) * t989 + t1237 + t583 + t983) + (t42 * t789 + t104 * (t1095 * t1260 - t576) - t41 * t697) * t706 + ((t1260 * t654 - t1191 - t590) * t1133 + (t104 * (-t590 + t789) - t105 * t697) * t709) * qJD(1) - (-t104 + t1263 + t721) * t105) * m(6) + t1125 / 0.2e1 + t1126 / 0.2e1 + (((t1275 * t709 - t1266 + t1346) * t709 + (t1275 * t706 + t1347 + t1377) * t706) * qJD(2) + t1312 + t1338) * t949 + (-t1379 * qJD(2) + t545 * t671 + t546 * t670 + t588 * t708 + t589 * t705 + t654 * t923 + t655 * t924) * qJD(1) - (-t1304 + t1305 + t1311) * t1009 / 0.2e1 + t1127 / 0.2e1 + t1128 / 0.2e1 + (t109 + t123 + t95) * t1199 + (t108 + t124) * t1200 + t1251 * t1210 + t1252 * t1211 + (t192 + (t179 + (t452 * t709 + t453 * t706) * t654 + t931 + t1061) * t616 + (-t454 * t1094 + t1117 + t178 + (t452 * t706 - t453 * t709) * t654 + t1060) * t615) * t1198 + t1232 + ((t1266 * t706 + ((t1381 + t1396) * t709 + t1348 + t1376 + t1395) * t709) * qJD(2) + t1345) * t946 + (-(-qJD(1) * t457 - t1102 - t141 + t749) * t142 + t97 * (-t457 + t1034) + t141 * (t675 - t919) + t96 * (-t697 * t706 + t458 + t559) + t142 * (t1029 + t1237) + (t1124 * t554 - t142 * t503) * t698 + ((-t141 * rSges(5,3) + t142 * (-t590 - t1182)) * t706 + (t141 * (-t555 - t590) - t142 * t697) * t709) * qJD(1)) * m(5) - t1175 / 0.2e1 + (t94 - t1001 + (t181 - t801 - t1060) * t616 + (t706 * t927 + t180 - t398) * t615 + ((t451 + t842) * t615 + t927 * t616) * t709) * t1201 + t1172 / 0.2e1 + t1314 * t1208 + (t1313 + t1316) * t1207 + (t1303 + t1306) * t676 / 0.2e1 - t1173 / 0.2e1 + t1174 / 0.2e1; t712 + ((-t1265 * t676 + t1267) * t706 + ((t1264 * t706 + t1281) * qJD(2) + t1280) * t709) * t949 + ((-t1009 * t1264 - t1267) * t709 + ((t1265 * t709 + t1281) * qJD(2) + t1280) * t706) * t946 + ((t1218 * t671 + t1219 * t708 + t770 * t670 + t769 * t705) * qJD(2) + (t1025 * t708 + t1026 * t705 + t1031 * t671 + t1032 * t670) * qJD(1)) * t1189 + (t1304 * t709 + t1303 * t706 + (t1301 * t709 + t1302 * t706) * qJD(1)) * t1188 + (t15 * (t825 + t907) + t61 * (t766 + t1217) + t782 * t1273 + (t25 * t782 + t79 * t755 + (-t977 - t1072) * t1136) * t706 - t79 * (t736 + t830) - t61 * (t737 + t1050) - t719 - (t883 * t1192 + (t61 * t936 + t708 * t883) * pkin(2)) * qJD(2) + (t709 * t755 + t1270 + t1274) * t78) * m(7) + (t24 * (t906 + t907) + t87 * (t766 + t1231) + t808 * t1272 + (t41 * t808 + t105 * t767 + (-t307 - t977) * t1134) * t706 - t105 * (t738 + t830) - t87 * (t739 + t1050) - (-t864 * t1192 + (-t708 * t864 + t87 * t936) * pkin(2)) * qJD(2) + (t709 * t767 + t1274 + t1290) * t104) * m(6) + (t141 * t1027 + t51 * (t907 + t1054) + t112 * (t826 + t985) + (t141 * t774 + (qJD(1) * t142 + t97) * t829) * t709 + (t96 * t829 + t142 * t774 + (t141 * t554 + t112 * (-t458 - t1063)) * qJD(1)) * t706 - t141 * (t915 + t917) - t142 * (t830 + t922) - t112 * (t930 + t1050) - (t851 * t1192 + (t112 * t936 + t708 * t851) * pkin(2)) * qJD(2)) * m(5) + (t132 * t894 - t184 * pkin(2) * t965 + (t1009 * t318 + t319 * t676 + t981) * t788 + t177 * t980 + (-t184 * t547 + t177 * t318 + (t1271 * t478 + t185 * t943) * qJD(1)) * t709 + (t131 * t943 + t185 * (-t547 - t994) + t177 * t319 + (t1121 + t1271 * (-t466 - t479)) * qJD(1)) * t706 - (t184 * t535 + t185 * (-t536 - t995)) * qJD(1) - (t177 * pkin(2) * t936 + (-t177 * t536 + t184 * t942) * t709 + (-t177 * t535 + t185 * t942) * t706) * qJD(2)) * m(4) + (-(t329 * t568 - t1118) * qJD(1) - (t325 * (-t568 * t706 - t569 * t709) + t845 * t892) * qJD(2) + 0.2e1 * t325 * (t374 * t709 + t375 * t706 + (t537 * t709 - t538 * t706) * qJD(1)) + t845 * t591 + (-t190 * t706 - t191 * t709 + (-t330 * t709 + t1120) * qJD(1)) * t624) * m(3) + (t1306 * qJD(1) + ((t1346 * qJD(1) + t1342 * t709) * t709 + (t1332 * t706 + t1347 * qJD(1) + (-t1333 + t1343) * t709) * t706) * t1327) * t1196 + (t1305 * qJD(1) + ((t1348 * qJD(1) + t1333 * t709) * t709 + (t1343 * t706 + t1371 * qJD(1) + (-t1332 + t1342) * t709) * t706) * t1327) * t1195 + (t1312 + t1339) * t951 + (t1311 + t1340) * t950; 0.2e1 * (t1195 * t25 + t1196 * t26) * m(7) + 0.2e1 * (-t1171 / 0.2e1 + t1170 / 0.2e1) * m(6) + 0.2e1 * (t1195 * t96 + t1196 * t97) * m(5) + 0.2e1 * (t1195 * t131 + t1196 * t132) * m(4); t712 + (t15 * t825 + (t79 * t987 + (-t507 - t1072) * t1136) * t706 - t736 * t79 - t719 + (t1259 + t1273) * (-t556 - t1062) + (t709 * t987 + t1270 + t510) * t78 + (t982 - t737 + t1217) * t61) * m(7) + (-t105 * t738 + t24 * t906 + (t105 * t1074 + t1067 * t1134) * t706 + (-t739 + t982 + t1231) * t87 + (t41 * t706 + t1272) * (-t430 - t556) + (t1074 * t709 + t1290 + t510) * t104) * m(6) + (-t141 * t917 - t142 * t922 + t51 * t1054 + t851 * t480 + (-t96 * t706 - t97 * t709 + (-t142 * t709 + t1124) * qJD(1)) * t554 + (-t1011 * t458 - t930 + t985) * t112) * m(5); t1325 * t1097 / 0.2e1 + t1326 * t1096 / 0.2e1 + (t1242 * t654 - t1252 * t655) * t1211 + (t1243 * t654 - t1251 * t655) * t1210 + (-t1245 * t709 + t1246 * t542 - t1247 * t541) * t1209 + ((t1243 * t698 - t1314) * t655 + (qJD(1) * t1277 + t1251 * t698 + t1322 * t709 + t1323 * t706) * t654) * t1208 + ((t1242 * t698 - t1313) * t655 + (qJD(1) * t1276 + t1252 * t698 + t1320 * t709 + t1321 * t706) * t654) * t1207 + (-t1245 * t706 + t1246 * t540 + t1247 * t539) * t1206 + (t1341 * t655 + (t1246 * t707 + t1247 * t704) * t654) * t1203 + ((t1244 * t698 - t1350) * t655 + (qJD(1) * t1278 + t1250 * t698 + t1318 * t709 + t1319 * t706) * t654) * t1202 - (t1127 + t1128 + t1174 - t1175 + t1125 + t1126 + t1172 - t1173 + t1232) * t655 / 0.2e1 + (-(t1064 * t79 - t1065 * t78) * t612 - (t1064 * t61 + t78 * t897) * t509 - (t1065 * t61 + t79 * t897) * t508 + (-t26 * t1073 - t78 * t1076 - t25 * t1072 - t79 * t1077 + (-t1269 * t709 + t771 * t706) * t698) * t655 + ((t1072 * t79 + t1073 * t78) * t698 + (qJD(1) * t771 - t1062 * t25 - t1073 * t15 - t1075 * t79 - t1076 * t61) * t709 + (qJD(1) * t1269 + t26 * t1062 - t15 * t1072 + t78 * t1075 - t61 * t1077) * t706) * t654) * m(7) + ((t104 * t164 - t105 * t162 + t42 * t303 - t41 * t307 + (t87 * t846 + (t104 * t706 - t105 * t709) * t430) * t698) * t655 + (t104 * (t231 * t706 - t303 * t698) + t105 * (-t231 * t709 + t307 * t698) + t24 * t846 + t87 * (-t1010 * t307 - t1011 * t303 - t162 * t706 + t164 * t709) + (qJD(1) * t864 + t1170 - t1171) * t430) * t654 - t104 * (-t362 * t612 - t501 * t509) - t105 * (t364 * t612 - t501 * t508) - t87 * (t362 * t508 + t364 * t509)) * m(6) + ((t1244 * t654 - t1250 * t655) * qJD(5) + t1315) * t1099 / 0.2e1 + t1317 * (t654 * t950 + t706 * t953) + t1316 * (-t971 / 0.2e1 + t709 * t953); 0.2e1 * ((t1088 * t79 + t616 * t78 - t15) * t1212 + (t508 * t79 + t509 * t78) * t1213) * t655 + 0.2e1 * ((t1010 * t79 - t1011 * t78 + t26 * t709 + t61 * t698 + t1259) * t1212 + (t61 * (t508 * t706 + t509 * t709) + (-t706 * t78 + t709 * t79) * t612) * t1213) * t654;];
tauc  = t1(:);
