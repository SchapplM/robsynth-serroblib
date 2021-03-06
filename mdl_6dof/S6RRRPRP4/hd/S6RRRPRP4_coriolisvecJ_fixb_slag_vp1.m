% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:36
% EndTime: 2019-03-09 16:45:39
% DurationCPUTime: 107.72s
% Computational Cost: add. (46535->1414), mult. (64685->1786), div. (0->0), fcn. (59124->8), ass. (0->716)
t1441 = Icges(6,4) - Icges(7,5);
t1382 = Icges(6,1) + Icges(7,1);
t1381 = Icges(7,4) + Icges(6,5);
t1380 = Icges(6,2) + Icges(7,3);
t1435 = Icges(6,6) - Icges(7,6);
t1440 = Icges(5,4) - Icges(4,5);
t1439 = Icges(5,5) - Icges(4,6);
t713 = cos(qJ(5));
t1438 = t1441 * t713;
t710 = sin(qJ(5));
t1437 = t1441 * t710;
t1436 = Icges(5,1) + Icges(4,3);
t709 = qJ(2) + qJ(3);
t689 = cos(t709);
t674 = Icges(4,4) * t689;
t688 = sin(t709);
t590 = Icges(4,1) * t688 + t674;
t1152 = Icges(5,6) * t689;
t873 = Icges(5,2) * t688 + t1152;
t1434 = t590 + t873;
t1433 = -t1380 * t713 - t1437;
t1432 = t1381 * t710 + t1435 * t713;
t1431 = t1382 * t710 + t1438;
t1426 = Icges(7,2) + Icges(6,3);
t712 = sin(qJ(1));
t1111 = t689 * t712;
t1114 = t688 * t712;
t715 = cos(qJ(1));
t1154 = Icges(4,6) * t715;
t464 = Icges(4,4) * t1111 - Icges(4,2) * t1114 - t1154;
t1163 = Icges(5,5) * t715;
t469 = Icges(5,6) * t1111 - Icges(5,3) * t1114 + t1163;
t1430 = t464 + t469;
t885 = -Icges(4,2) * t688 + t674;
t465 = Icges(4,6) * t712 + t715 * t885;
t872 = -Icges(5,3) * t688 + t1152;
t468 = Icges(5,5) * t712 - t715 * t872;
t1425 = t465 - t468;
t1164 = Icges(4,5) * t715;
t650 = Icges(4,4) * t1114;
t466 = Icges(4,1) * t1111 - t1164 - t650;
t1170 = Icges(5,4) * t715;
t645 = Icges(5,6) * t1114;
t471 = Icges(5,2) * t1111 + t1170 - t645;
t1395 = t466 + t471;
t1172 = Icges(4,4) * t688;
t591 = Icges(4,1) * t689 - t1172;
t467 = Icges(4,5) * t712 + t591 * t715;
t1109 = t689 * t715;
t1171 = Icges(5,4) * t712;
t1113 = t688 * t715;
t646 = Icges(5,6) * t1113;
t470 = -Icges(5,2) * t1109 + t1171 + t646;
t1429 = t467 - t470;
t1394 = t1439 * t688 - t1440 * t689;
t1153 = Icges(5,6) * t688;
t874 = Icges(5,2) * t689 - t1153;
t1428 = -t591 - t874;
t1427 = t872 + t885;
t584 = Icges(4,5) * t688 + Icges(4,6) * t689;
t883 = Icges(5,4) * t688 + Icges(5,5) * t689;
t1424 = t584 - t883;
t1357 = t1433 * t689 + t1435 * t688;
t1356 = -t1381 * t688 + t1431 * t689;
t1423 = t1433 * t688;
t1422 = t1432 * t688;
t1421 = t1431 * t688;
t1420 = -t1380 * t710 + t1438;
t1419 = -t1381 * t713 + t1435 * t710;
t1418 = -t1382 * t713 + t1437;
t588 = Icges(4,2) * t689 + t1172;
t1417 = t588 * t715 - t1429;
t1416 = t1436 * t715;
t1415 = -t1427 - t1434;
t871 = Icges(5,3) * t689 + t1153;
t1414 = -t871 - t588 - t1428;
t1413 = t1434 * t712 + t1430;
t1412 = -t1434 * t715 - t1425;
t1377 = t1426 * t688 - t1432 * t689;
t1101 = t713 * t715;
t1107 = t710 * t712;
t557 = -t1101 * t688 + t1107;
t1103 = t712 * t713;
t1106 = t710 * t715;
t558 = t1106 * t688 + t1103;
t1257 = t1109 * t1377 - t1356 * t558 - t1357 * t557;
t559 = t1103 * t688 + t1106;
t512 = Icges(7,5) * t559;
t560 = -t1107 * t688 + t1101;
t300 = -Icges(7,1) * t560 + Icges(7,4) * t1111 - t512;
t515 = Icges(6,4) * t559;
t302 = Icges(6,1) * t560 - Icges(6,5) * t1111 - t515;
t1297 = t300 - t302;
t513 = Icges(7,5) * t560;
t287 = -Icges(7,6) * t1111 + Icges(7,3) * t559 + t513;
t516 = Icges(6,4) * t560;
t297 = Icges(6,2) * t559 + Icges(6,6) * t1111 - t516;
t1299 = t287 + t297;
t291 = -Icges(6,5) * t560 + Icges(6,6) * t559 + Icges(6,3) * t1111;
t294 = -Icges(7,4) * t560 + Icges(7,2) * t1111 - Icges(7,6) * t559;
t1364 = t291 + t294;
t1370 = t1109 * t1364 + t1297 * t558 - t1299 * t557;
t1161 = Icges(7,5) * t557;
t298 = Icges(7,1) * t558 + Icges(7,4) * t1109 + t1161;
t514 = Icges(6,4) * t557;
t301 = Icges(6,1) * t558 + Icges(6,5) * t1109 - t514;
t1398 = t298 + t301;
t289 = Icges(6,5) * t558 - Icges(6,6) * t557 + Icges(6,3) * t1109;
t292 = Icges(7,4) * t558 + Icges(7,2) * t1109 + Icges(7,6) * t557;
t1399 = t289 + t292;
t511 = Icges(7,5) * t558;
t286 = Icges(7,6) * t1109 + Icges(7,3) * t557 + t511;
t1169 = Icges(6,4) * t558;
t295 = -Icges(6,2) * t557 + Icges(6,6) * t1109 + t1169;
t1400 = t286 - t295;
t1371 = t1109 * t1399 + t1398 * t558 + t1400 * t557;
t1016 = qJD(5) * t715;
t687 = qJD(2) * t712;
t620 = qJD(3) * t712 + t687;
t504 = t1016 * t689 + t620;
t1017 = qJD(5) * t712;
t706 = qJD(2) + qJD(3);
t621 = t706 * t715;
t505 = -t1017 * t689 + t621;
t1020 = qJD(5) * t688;
t659 = qJD(1) + t1020;
t1325 = t1257 * t659 - t1370 * t505 + t1371 * t504;
t1256 = t1111 * t1377 + t1356 * t560 + t1357 * t559;
t1368 = t1111 * t1364 - t1297 * t560 + t1299 * t559;
t1369 = t1111 * t1399 - t1398 * t560 - t1400 * t559;
t1324 = t1256 * t659 - t1368 * t505 + t1369 * t504;
t1396 = t1111 * t1440 - t1439 * t1114 + t1416;
t1411 = t1394 * t715 + t1436 * t712;
t841 = t588 * t688 - t590 * t689;
t1410 = -t688 * t871 + t689 * t873 - t841;
t1393 = -t1395 * t689 + t1430 * t688;
t1108 = t706 * t712;
t1001 = t689 * t1108;
t1024 = qJD(1) * t715;
t988 = t688 * t1024;
t1242 = t988 + t1001;
t257 = t1107 * t659 + (-t1016 - t1242) * t713;
t835 = t713 * t659;
t931 = qJD(1) * t688 + qJD(5);
t258 = t712 * t835 + (t715 * t931 + t1001) * t710;
t1004 = t688 * t1108;
t986 = t689 * t1024;
t797 = t986 - t1004;
t154 = Icges(7,5) * t258 + Icges(7,6) * t797 + Icges(7,3) * t257;
t160 = Icges(6,4) * t258 - Icges(6,2) * t257 + Icges(6,6) * t797;
t1409 = t154 - t160;
t1000 = t621 * t689;
t259 = qJD(1) * t559 + qJD(5) * t558 - t1000 * t713;
t260 = t715 * t835 + (-t712 * t931 + t1000) * t710;
t1002 = t688 * t621;
t1025 = qJD(1) * t712;
t987 = t689 * t1025;
t796 = -t987 - t1002;
t155 = Icges(7,5) * t260 + Icges(7,6) * t796 + Icges(7,3) * t259;
t161 = Icges(6,4) * t260 - Icges(6,2) * t259 + Icges(6,6) * t796;
t1408 = t155 - t161;
t156 = Icges(6,5) * t258 - Icges(6,6) * t257 + Icges(6,3) * t797;
t158 = Icges(7,4) * t258 + Icges(7,2) * t797 + Icges(7,6) * t257;
t1407 = t156 + t158;
t157 = Icges(6,5) * t260 - Icges(6,6) * t259 + Icges(6,3) * t796;
t159 = Icges(7,4) * t260 + Icges(7,2) * t796 + Icges(7,6) * t259;
t1406 = t157 + t159;
t162 = Icges(7,1) * t258 + Icges(7,4) * t797 + Icges(7,5) * t257;
t164 = Icges(6,1) * t258 - Icges(6,4) * t257 + Icges(6,5) * t797;
t1405 = t162 + t164;
t163 = Icges(7,1) * t260 + Icges(7,4) * t796 + Icges(7,5) * t259;
t165 = Icges(6,1) * t260 - Icges(6,4) * t259 + Icges(6,5) * t796;
t1404 = t163 + t165;
t1403 = t1423 * t706 + (qJD(5) * t1420 - t1435 * t706) * t689;
t1402 = t1422 * t706 + (qJD(5) * t1419 + t1426 * t706) * t689;
t1401 = t1421 * t706 + (t1418 * qJD(5) + t1381 * t706) * t689;
t1397 = -t467 * t1111 - t468 * t1114;
t1392 = -t1356 * t710 + t1357 * t713;
t1391 = t1414 * t706;
t1390 = t1415 * t706;
t1389 = t1412 * t706 + (t1428 * t712 + t1164 - t1170) * qJD(1);
t1388 = t1413 * t706 + (-t715 * t874 + t1171 - t467) * qJD(1);
t814 = t706 * t871;
t1387 = t715 * t814 + t1417 * t706 + (t1427 * t712 - t1154 + t1163) * qJD(1);
t1386 = qJD(1) * t1425 - t1108 * t588 + t1395 * t706 - t712 * t814;
t1385 = t1424 * t712;
t1119 = t584 * t715;
t528 = t883 * t715;
t1384 = -t528 + t1119;
t1383 = t465 * t688 + t470 * t689;
t1332 = rSges(7,1) + pkin(5);
t1316 = rSges(7,3) + qJ(6);
t1378 = t1410 * t715 + t1385;
t1376 = t1393 * t712;
t1375 = t1411 * t715 + t1397;
t1374 = t1411 * qJD(1);
t1269 = -t467 * t1109 - t468 * t1113 - t1411 * t712;
t1373 = -t466 * t1109 + t469 * t1113 + t1396 * t712;
t1372 = t1396 * t715;
t1331 = t1111 * t1406 + t1398 * t258 + t1399 * t797 + t1400 * t257 - t1404 * t560 - t1408 * t559;
t1330 = t1407 * t1111 + t1297 * t258 - t1299 * t257 + t1364 * t797 - t1405 * t560 - t1409 * t559;
t1329 = t1109 * t1406 + t1398 * t260 + t1399 * t796 + t1400 * t259 + t1404 * t558 + t1408 * t557;
t1328 = t1407 * t1109 + t1297 * t260 - t1299 * t259 + t1364 * t796 + t1405 * t558 + t1409 * t557;
t1322 = t1111 * t1402 - t1356 * t258 - t1357 * t257 + t1377 * t797 - t1401 * t560 - t1403 * t559;
t1321 = t1109 * t1402 - t1356 * t260 - t1357 * t259 + t1377 * t796 + t1401 * t558 + t1403 * t557;
t856 = t286 * t713 - t298 * t710;
t140 = t292 * t688 + t689 * t856;
t854 = t295 * t713 + t301 * t710;
t142 = t289 * t688 - t689 * t854;
t1367 = t140 + t142;
t855 = -t287 * t713 - t300 * t710;
t141 = t294 * t688 + t689 * t855;
t853 = t297 * t713 - t302 * t710;
t143 = t291 * t688 - t689 * t853;
t1366 = t141 + t143;
t1255 = t1377 * t688 - t1392 * t689;
t1305 = t1372 - t1376;
t1304 = -t1111 * t470 - t1114 * t465 - t1375;
t1303 = t1109 * t471 - t1113 * t464 - t1373;
t1302 = -t1109 * t470 - t1113 * t465 - t1269;
t1363 = (t1426 * t689 + t1392 + t1422) * t659 + (t1377 * t712 - t853 + t855) * t505 + (-t1377 * t715 + t854 - t856) * t504;
t1362 = qJD(1) * t1424 + t1390 * t688 + t1391 * t689;
t1361 = qJD(1) * t1396 + t1386 * t688 + t1388 * t689;
t1360 = t1387 * t688 + t1389 * t689 + t1374;
t1359 = qJD(1) * t1414 + t1412 * t620 + t1413 * t621;
t895 = pkin(5) * t710 - qJ(6) * t713;
t897 = rSges(7,1) * t710 - rSges(7,3) * t713;
t1063 = rSges(7,2) * t688 + (-t895 - t897) * t689;
t1355 = (-t645 - t650 + (-Icges(4,2) - Icges(5,3)) * t1111 + t1395) * t621 + (Icges(5,3) * t1109 + t1417 + t646) * t620 + t1415 * qJD(1);
t1354 = qJD(1) * t1410 - t1394 * t706;
t1353 = t1384 * t706 + (t1394 * t712 + t467 * t689 + t468 * t688 - t1383 - t1416) * qJD(1);
t1352 = qJD(1) * t1393 - t1385 * t706 + t1374;
t1346 = t1316 * t559 + t1332 * t560;
t1083 = rSges(7,2) * t1111 - t1346;
t507 = qJD(6) * t557;
t1351 = -t1063 * t505 - t1083 * t659 + t507;
t1252 = (t1382 * t557 + t1169 - t1400 - t511) * t504 - (-t1382 * t559 + t1299 + t513 - t516) * t505 + (-t1418 * t689 + t1357) * t659;
t1350 = (-t1380 * t558 + t1161 + t1398 - t514) * t504 - (t1380 * t560 + t1297 - t512 + t515) * t505 + (-t1420 * t689 - t1356) * t659;
t1349 = (t1403 * t713 - t1401 * t710 + t1377 * t706 + (t1356 * t713 + t1357 * t710) * qJD(5)) * t689 + (t1392 * t706 + t1402) * t688;
t506 = t559 * qJD(6);
t1348 = -t1316 * t257 - t1332 * t258 + t506;
t1347 = t1378 * qJD(1);
t1345 = (-t1381 * t559 - t1435 * t560) * t505 + (-t1381 * t557 - t1435 * t558) * t504 + t1419 * t659 * t689;
t1344 = t1316 * t259 + t1332 * t260 + t507;
t216 = -t712 * t841 - t1119;
t219 = -t1111 * t873 + t1114 * t871 - t528;
t1343 = (-t216 + t219) * qJD(1);
t1010 = qJD(1) * qJD(3);
t1011 = qJD(1) * qJD(2);
t677 = t712 * t1011;
t598 = t712 * t1010 + t677;
t367 = qJD(5) * t797 + t598;
t678 = t715 * t1011;
t599 = t715 * t1010 + t678;
t368 = qJD(5) * t796 + t599;
t1019 = qJD(5) * t689;
t977 = t706 * t1019;
t1334 = t1256 * t977 + t1322 * t659 - t1330 * t505 + t1331 * t504 + t1368 * t367 + t1369 * t368;
t1333 = t1257 * t977 + t1321 * t659 - t1328 * t505 + t1329 * t504 + t1370 * t367 + t1371 * t368;
t41 = (-t706 * t856 + t159) * t688 + (t155 * t713 - t163 * t710 + t292 * t706 + (-t286 * t710 - t298 * t713) * qJD(5)) * t689;
t43 = (t706 * t854 + t157) * t688 + (-t161 * t713 - t165 * t710 + t289 * t706 + (t295 * t710 - t301 * t713) * qJD(5)) * t689;
t1327 = t41 + t43;
t42 = (-t706 * t855 + t158) * t688 + (t154 * t713 - t162 * t710 + t294 * t706 + (t287 * t710 - t300 * t713) * qJD(5)) * t689;
t44 = (t706 * t853 + t156) * t688 + (-t160 * t713 - t164 * t710 + t291 * t706 + (t297 * t710 + t302 * t713) * qJD(5)) * t689;
t1326 = t42 + t44;
t1323 = t1255 * t659 - t1366 * t505 + t1367 * t504;
t1320 = -t1352 * t712 + t1361 * t715;
t1319 = -t1353 * t712 + t1360 * t715;
t1318 = t1352 * t715 + t1361 * t712;
t1317 = t1353 * t715 + t1360 * t712;
t1115 = t688 * t706;
t1315 = t1304 * t620 - t1305 * t621 - t1343;
t1314 = t1302 * t620 - t1303 * t621 + t1347;
t1313 = -t1354 * t712 + t1362 * t715;
t1312 = t1354 * t715 + t1362 * t712;
t1307 = -t1386 * t689 + t1388 * t688;
t1306 = -t1387 * t689 + t1389 * t688;
t1098 = rSges(7,2) * t797 - t1348;
t1097 = rSges(7,2) * t796 + t1344;
t1301 = t1395 * t688 + t1430 * t689;
t1300 = t1425 * t689 + t1429 * t688;
t1085 = rSges(7,2) * t1109 + t1316 * t557 + t1332 * t558;
t1296 = t1357 * t712;
t1295 = t1357 * t715;
t1294 = t1356 * t712;
t1293 = t1356 * t715;
t1292 = t1063 * t712;
t1291 = t1063 * t715;
t1290 = -t1435 * t689 + t1423;
t1289 = t1381 * t689 + t1421;
t1187 = rSges(7,2) * t689;
t819 = t897 * t688;
t1288 = -t895 * t688 - t1187 - t819;
t1280 = t1355 * t688 + t1359 * t689;
t1279 = t1363 * t689;
t1278 = qJD(1) * t1394 - t1384 * t620 + t1385 * t621;
t1277 = t1364 * t505 - t1377 * t659 - t1399 * t504;
t1250 = t1366 * t712 + t1367 * t715;
t1276 = -t1366 * t715 + t1367 * t712;
t1275 = -t1368 * t715 + t1369 * t712;
t1249 = t1368 * t712 + t1369 * t715;
t1248 = t1370 * t712 + t1371 * t715;
t1274 = -t1370 * t715 + t1371 * t712;
t1271 = t1383 + t1396;
t1015 = qJD(6) * t713;
t635 = t689 * t1015;
t643 = qJ(4) * t1113;
t553 = pkin(3) * t1109 + t643;
t703 = t712 * pkin(4);
t578 = pkin(9) * t1109 + t703;
t1048 = t553 + t578;
t579 = pkin(4) * t715 - pkin(9) * t1111;
t1023 = qJD(2) * t715;
t714 = cos(qJ(2));
t1198 = pkin(2) * t714;
t679 = pkin(1) + t1198;
t716 = -pkin(8) - pkin(7);
t683 = t715 * t716;
t1037 = -t712 * t679 - t683;
t704 = t715 * pkin(7);
t641 = pkin(1) * t712 - t704;
t458 = t641 + t1037;
t702 = t712 * pkin(7);
t642 = t715 * pkin(1) + t702;
t652 = t715 * t679;
t936 = -t712 * t716 + t652;
t459 = t936 - t642;
t1073 = t459 * t1023 - t458 * t687;
t595 = pkin(3) * t689 + qJ(4) * t688;
t545 = t595 * t712;
t830 = -qJD(4) * t689 + t620 * t545 + t1073;
t759 = t1048 * t621 - t620 * t579 + t830;
t70 = t1083 * t504 + t1085 * t505 + t635 + t759;
t1197 = pkin(9) * t688;
t564 = t621 * t1197;
t592 = pkin(3) * t688 - qJ(4) * t689;
t1021 = qJD(4) * t715;
t639 = t688 * t1021;
t711 = sin(qJ(2));
t982 = t711 * t1023;
t932 = pkin(2) * t982;
t834 = t639 - t932;
t818 = -t621 * t592 + t834;
t1062 = t458 - t641;
t993 = -t545 + t1062;
t746 = (t579 + t993) * qJD(1) - t564 + t818;
t94 = t1351 + t746;
t959 = t94 * t1063;
t1270 = -t1085 * t70 + t959;
t903 = rSges(6,1) * t560 - rSges(6,2) * t559;
t310 = rSges(6,3) * t1111 - t903;
t901 = rSges(6,1) * t710 + rSges(6,2) * t713;
t771 = -rSges(6,3) * t688 + t689 * t901;
t1268 = -t310 * t659 + t505 * t771;
t563 = t620 * t1197;
t1022 = qJD(4) * t712;
t637 = t688 * t1022;
t1199 = pkin(2) * t711;
t664 = t687 * t1199;
t991 = -t620 * t592 + t637 - t664;
t1061 = t459 + t642;
t992 = t553 + t1061;
t761 = qJD(1) * (t578 + t992) - t563 + t991;
t95 = -t1063 * t504 + t1085 * t659 - t506 + t761;
t1267 = 0.2e1 * qJD(2);
t1266 = rSges(3,2) * t711;
t795 = -t1019 * t710 - t713 * t1115;
t898 = -rSges(7,1) * t713 - rSges(7,3) * t710;
t1096 = t706 * t819 + (rSges(7,2) * t706 + qJD(5) * t898) * t689 + t635 + t795 * qJ(6) + (-t1019 * t713 + t1115 * t710) * pkin(5);
t1112 = t689 * t706;
t685 = pkin(4) * t1024;
t407 = pkin(9) * t796 + t685;
t419 = qJ(4) * t1115 + (pkin(3) * t706 - qJD(4)) * t689;
t605 = qJ(4) * t1000;
t989 = t688 * t1025;
t244 = pkin(3) * t796 - qJ(4) * t989 + t605 + t639;
t638 = t689 * t1022;
t1099 = t714 * qJD(2) ^ 2;
t1194 = pkin(1) - t679;
t684 = pkin(7) * t1024;
t363 = -t932 - t684 + (t1194 * t712 - t683) * qJD(1);
t577 = qJD(1) * (-pkin(1) * t1025 + t684);
t755 = qJD(1) * t363 + t577 + (-t1099 * t712 - t678 * t711) * pkin(2);
t734 = t706 * t638 + t755 + (t244 + t639) * qJD(1);
t731 = qJD(1) * t407 - t620 * t419 + t734;
t961 = -t592 - t1197;
t26 = qJD(6) * t257 + t1097 * t659 + t961 * t599 - t1096 * t504 - t1063 * t368 + (-pkin(9) * t620 + qJD(5) * t1085) * t1112 + t731;
t1265 = t26 * t712;
t1094 = t260 * rSges(6,1) - t259 * rSges(6,2);
t169 = rSges(6,3) * t796 + t1094;
t820 = t901 * t688;
t902 = -rSges(6,1) * t713 + rSges(6,2) * t710;
t236 = t706 * t820 + (rSges(6,3) * t706 + qJD(5) * t902) * t689;
t306 = t558 * rSges(6,1) - t557 * rSges(6,2) + rSges(6,3) * t1109;
t45 = t306 * t977 + t169 * t659 - t236 * t504 + t368 * t771 - t599 * t592 + (-t1112 * t620 - t599 * t688) * pkin(9) + t731;
t1264 = t45 * t712;
t1190 = rSges(5,2) * t688;
t896 = rSges(5,3) * t689 + t1190;
t1043 = -t592 + t896;
t1186 = rSges(5,3) * t688;
t1189 = rSges(5,2) * t689;
t596 = t1186 - t1189;
t495 = t596 * t706;
t1071 = -t419 - t495;
t924 = rSges(5,1) * t1024 - rSges(5,2) * t796 + rSges(5,3) * t1000;
t280 = -rSges(5,3) * t989 + t924;
t96 = qJD(1) * t280 + t1043 * t599 + t1071 * t620 + t734;
t1263 = t712 * t96;
t1110 = t689 * t713;
t699 = t712 * rSges(5,1);
t478 = -rSges(5,2) * t1109 + rSges(5,3) * t1113 + t699;
t145 = t896 * t620 + (t478 + t992) * qJD(1) + t991;
t543 = t896 * t712;
t279 = t706 * t543 + (t596 * t715 + t699) * qJD(1);
t640 = t689 * t1021;
t839 = -pkin(2) * t1099 * t715 + qJD(1) * t664;
t780 = -t621 * t419 + t598 * t592 + t706 * t640 + t839;
t1035 = t716 * t1025 + t664;
t364 = (-t1194 * t715 - t702) * qJD(1) - t1035;
t612 = t642 * qJD(1);
t1080 = -t364 - t612;
t1039 = -pkin(3) * t1004 + t637;
t245 = pkin(3) * t986 + qJ(4) * t1242 + t1039;
t805 = -t245 - t637 + t1080;
t97 = -t495 * t621 - t896 * t598 + (-t279 + t805) * qJD(1) + t780;
t1259 = qJD(1) * t145 + t97;
t610 = pkin(9) * t1004;
t406 = qJD(1) * t578 - t610;
t565 = pkin(9) * t1000;
t730 = t598 * t1197 + (-t406 + t805) * qJD(1) - t706 * t565 + t780;
t955 = t1083 * t689;
t27 = -qJD(5) * t706 * t955 + qJD(6) * t259 + t1063 * t367 - t1096 * t505 - t1098 * t659 + t730;
t1258 = qJD(1) * t95 + t27;
t614 = qJD(1) * t641;
t1254 = qJD(1) * t458 - t614;
t1251 = t1345 * t689;
t555 = t592 * t1025;
t622 = pkin(9) * t989;
t1046 = t555 + t622;
t385 = t771 * t712;
t1184 = rSges(6,3) * t689;
t448 = t820 + t1184;
t1038 = -t621 * t595 + t640;
t542 = t592 * t712;
t907 = qJD(1) * t542 + t1038;
t804 = -t565 + t622 + t907;
t1245 = -t1025 * t771 + t385 * t659 + t505 * t448 + t1046 - t804;
t696 = Icges(3,4) * t714;
t886 = -Icges(3,2) * t711 + t696;
t628 = Icges(3,1) * t711 + t696;
t1239 = t1255 * t977 + t1349 * t659;
t476 = rSges(4,1) * t1111 - rSges(4,2) * t1114 - t715 * rSges(4,3);
t697 = t712 * rSges(4,3);
t477 = rSges(4,1) * t1109 - rSges(4,2) * t1113 + t697;
t176 = t476 * t620 + t477 * t621 + t1073;
t594 = rSges(4,1) * t688 + rSges(4,2) * t689;
t798 = -t594 * t621 - t932;
t185 = (-t476 + t1062) * qJD(1) + t798;
t186 = -t594 * t620 - t664 + (t477 + t1061) * qJD(1);
t544 = t594 * t712;
t552 = t594 * t715;
t1192 = rSges(4,1) * t689;
t597 = -rSges(4,2) * t688 + t1192;
t1238 = -t185 * (qJD(1) * t544 - t621 * t597) - t176 * (-t620 * t544 - t552 * t621) - t186 * (-qJD(1) * t552 - t597 * t620);
t778 = -qJD(1) * t545 + t1254 + t834;
t1237 = qJD(1) * t579 + t621 * t961 + t778;
t1149 = Icges(3,3) * t715;
t625 = Icges(3,5) * t714 - Icges(3,6) * t711;
t624 = Icges(3,5) * t711 + Icges(3,6) * t714;
t801 = qJD(2) * t624;
t501 = Icges(3,6) * t712 + t715 * t886;
t1122 = t501 * t711;
t1173 = Icges(3,4) * t711;
t629 = Icges(3,1) * t714 - t1173;
t503 = Icges(3,5) * t712 + t629 * t715;
t843 = -t503 * t714 + t1122;
t1232 = -t715 * t801 + (-t625 * t712 + t1149 + t843) * qJD(1);
t499 = Icges(3,3) * t712 + t625 * t715;
t1027 = qJD(1) * t499;
t1102 = t712 * t714;
t1105 = t711 * t712;
t1155 = Icges(3,6) * t715;
t500 = Icges(3,4) * t1102 - Icges(3,2) * t1105 - t1155;
t1123 = t500 * t711;
t1165 = Icges(3,5) * t715;
t670 = Icges(3,4) * t1105;
t502 = Icges(3,1) * t1102 - t1165 - t670;
t844 = -t502 * t714 + t1123;
t1231 = qJD(1) * t844 - t712 * t801 + t1027;
t626 = Icges(3,2) * t714 + t1173;
t840 = t626 * t711 - t628 * t714;
t1228 = t840 * qJD(1) + t625 * qJD(2);
t498 = Icges(3,5) * t1102 - Icges(3,6) * t1105 - t1149;
t197 = -t498 * t715 - t712 * t844;
t1224 = (-t626 * t715 + t503) * t712 - (-Icges(3,2) * t1102 + t502 - t670) * t715;
t114 = t306 * t659 + t504 * t771 + t761;
t904 = rSges(6,1) * t258 - rSges(6,2) * t257;
t167 = rSges(6,3) * t797 + t904;
t46 = -t167 * t659 - t236 * t505 - t310 * t977 - t367 * t771 + t730;
t1223 = (qJD(1) * t114 + t46) * t715 + t1264;
t1220 = m(5) / 0.2e1;
t1219 = m(6) / 0.2e1;
t1218 = m(7) / 0.2e1;
t1217 = -pkin(3) - pkin(9);
t1216 = t367 / 0.2e1;
t1215 = t368 / 0.2e1;
t1214 = -t504 / 0.2e1;
t1213 = t504 / 0.2e1;
t1212 = -t505 / 0.2e1;
t1211 = t505 / 0.2e1;
t1210 = t598 / 0.2e1;
t1209 = t599 / 0.2e1;
t1208 = -t620 / 0.2e1;
t1207 = t620 / 0.2e1;
t1206 = -t621 / 0.2e1;
t1205 = t621 / 0.2e1;
t1204 = -t659 / 0.2e1;
t1203 = t659 / 0.2e1;
t1201 = t712 / 0.2e1;
t1200 = -t715 / 0.2e1;
t1196 = -qJD(1) / 0.2e1;
t1195 = qJD(1) / 0.2e1;
t1193 = rSges(3,1) * t714;
t1191 = rSges(3,2) * t714;
t1183 = t41 * t504;
t1182 = t42 * t505;
t1181 = t43 * t504;
t1180 = t44 * t505;
t698 = t712 * rSges(3,3);
t1179 = t712 * t95;
t1146 = qJD(1) * t70;
t113 = t1268 + t746;
t1144 = t113 * t712;
t1143 = t114 * t712;
t1138 = t140 * t368;
t1137 = t141 * t367;
t1136 = t142 * t368;
t1135 = t143 * t367;
t1132 = t185 * t712;
t1034 = rSges(3,2) * t1105 + t715 * rSges(3,3);
t508 = rSges(3,1) * t1102 - t1034;
t631 = rSges(3,1) * t711 + t1191;
t980 = t631 * t1023;
t312 = -t980 + (-t508 - t641) * qJD(1);
t1131 = t312 * t712;
t1130 = t312 * t715;
t1100 = t714 * t715;
t1104 = t711 * t715;
t509 = rSges(3,1) * t1100 - rSges(3,2) * t1104 + t698;
t1051 = t509 + t642;
t983 = t631 * t687;
t313 = qJD(1) * t1051 - t983;
t574 = t631 * t715;
t1129 = t313 * t574;
t1118 = t620 * t712;
t1117 = t624 * t712;
t1116 = t624 * t715;
t1082 = t1316 * t558 - t1332 * t557;
t1081 = -t1316 * t560 + t1332 * t559;
t673 = qJD(4) * t688;
t1075 = -t620 * t542 + t673;
t1072 = -t712 * t458 + t715 * t459;
t1070 = t712 * t476 + t715 * t477;
t1069 = -t502 * t1100 - t712 * t498;
t1068 = t503 * t1100 + t712 * t499;
t1056 = t478 + t553;
t1055 = t712 * t545 + t715 * t553;
t550 = t592 * t715;
t1054 = -qJD(1) * t550 + t638;
t1049 = (pkin(5) * t713 + qJ(6) * t710 - t898) * t689;
t1047 = -t1025 * t896 + t555;
t1042 = rSges(4,2) * t989 + rSges(4,3) * t1024;
t1041 = -t626 + t629;
t1040 = t628 + t886;
t1036 = rSges(3,3) * t1024 + t1025 * t1266;
t100 = t306 * t505 + t310 * t504 + t759;
t1033 = qJD(1) * t100;
t957 = rSges(5,1) * t715 - rSges(5,3) * t1114;
t479 = rSges(5,2) * t1111 + t957;
t125 = t1056 * t621 - t479 * t620 + t830;
t1031 = qJD(1) * t125;
t1026 = qJD(1) * t625;
t316 = -t712 * t840 - t1116;
t1012 = t316 * qJD(1);
t1009 = -rSges(7,2) + t1217;
t1008 = -rSges(6,3) + t1217;
t1007 = pkin(9) * t1112;
t1006 = qJD(2) * t1198;
t1005 = qJ(4) * t1111;
t999 = t545 * t1024 + t715 * t244 + t712 * t245;
t277 = rSges(4,1) * t796 - rSges(4,2) * t1000 + t1042;
t278 = -t706 * t544 + (t597 * t715 + t697) * qJD(1);
t998 = t476 * t1024 + t715 * t277 + t712 * t278;
t997 = -t306 - t1048;
t996 = -t458 * t1024 + t715 * t363 + t712 * t364;
t990 = t652 + t553;
t984 = t711 * t1024;
t979 = t688 * t1017;
t978 = t688 * t1016;
t975 = -t1115 / 0.2e1;
t974 = t1112 / 0.2e1;
t971 = -pkin(1) - t1193;
t970 = t1025 / 0.2e1;
t969 = t1024 / 0.2e1;
t968 = -t687 / 0.2e1;
t965 = t1023 / 0.2e1;
t962 = -t594 - t1199;
t960 = t70 * t1083;
t956 = t711 * (-t712 ^ 2 - t715 ^ 2);
t430 = t503 * t1102;
t944 = t499 * t715 - t430;
t937 = -t498 + t1122;
t934 = t715 * t1009;
t933 = t715 * t1008;
t930 = -t1048 - t1085;
t929 = t1025 * t1063 + t1046;
t926 = t715 * t478 - t712 * t479 + t1055;
t925 = t715 * t578 - t712 * t579 + t1055;
t923 = t1035 - t1039;
t915 = qJD(5) * t974;
t914 = -t419 - t1007;
t913 = t1043 - t1199;
t909 = -t419 - t1006;
t496 = t597 * t706;
t908 = -t496 - t1006;
t905 = t1193 - t1266;
t894 = t610 + t923;
t893 = t990 + t578;
t870 = t113 * t715 + t1143;
t857 = -t185 * t715 - t186 * t712;
t852 = t306 * t712 - t310 * t715;
t851 = -t313 * t712 - t1130;
t314 = t500 * t714 + t502 * t711;
t315 = t501 * t714 + t503 * t711;
t838 = -t236 + t914;
t837 = t579 + t1037;
t833 = -t495 + t909;
t832 = t961 - t1199;
t831 = -t550 * t621 - t715 * t564 + t1075;
t829 = -t479 * t1024 + t712 * t279 + t715 * t280 + t999;
t828 = -t579 * t1024 + t712 * t406 + t715 * t407 + t999;
t827 = t715 * t306 + t712 * t310 + t925;
t445 = t620 * t689 + t988;
t822 = t914 - t1096;
t573 = t631 * t712;
t817 = t605 + t834;
t803 = qJD(2) * t628;
t802 = qJD(2) * t626;
t198 = -t1105 * t501 - t944;
t800 = (-t197 * t715 + t198 * t712) * qJD(2);
t199 = -t1104 * t500 - t1069;
t200 = -t1104 * t501 + t1068;
t799 = (-t199 * t715 + t200 * t712) * qJD(2);
t285 = (t508 * t712 + t509 * t715) * qJD(2);
t794 = t685 + t817;
t793 = t832 - t1063;
t788 = t909 - t1007;
t785 = t363 * t1023 + t364 * t687 - t458 * t678 - t459 * t677;
t784 = t95 * t1063 - t960;
t783 = t1083 * t712 + t1085 * t715 + t925;
t782 = t500 * t715 - t501 * t712;
t781 = -t595 - t679 - t1186;
t779 = -t643 - t652 - t703;
t777 = -qJ(4) * t1114 + t1037;
t776 = -t236 + t788;
t775 = t310 * t1024 + t712 * t167 + t715 * t169 + t828;
t769 = t788 - t1096;
t762 = (-t1040 * t711 + t1041 * t714) * qJD(1);
t758 = t620 * t245 + t599 * t545 + t706 * t673 + t785;
t756 = t1024 * t1083 + t1097 * t715 + t1098 * t712 + t828;
t747 = -pkin(9) * t445 - t595 * t620 + t1054;
t334 = qJD(1) * t501 - t712 * t802;
t336 = qJD(1) * t503 - t712 * t803;
t737 = qJD(1) * t498 - qJD(2) * t314 - t334 * t711 + t336 * t714;
t333 = -t715 * t802 + (-t712 * t886 + t1155) * qJD(1);
t335 = -t715 * t803 + (-t629 * t712 + t1165) * qJD(1);
t736 = -qJD(2) * t315 - t333 * t711 + t335 * t714 + t1027;
t603 = t886 * qJD(2);
t604 = t629 * qJD(2);
t735 = qJD(1) * t624 - t603 * t711 + t604 * t714 + (-t626 * t714 - t628 * t711) * qJD(2);
t551 = t896 * t715;
t733 = t125 * (t620 * t543 + (-t550 + t551) * t621 + t1075) + t145 * (qJD(1) * t551 + (-t595 - t596) * t620 + t1054);
t732 = t1288 * t505 + t1292 * t659 - t635 * t715 + t804;
t729 = -t1224 * t711 + t782 * t714;
t728 = t620 * t406 + (t244 + t407) * t621 - t1048 * t598 - t599 * t579 + t758;
t387 = t771 * t715;
t721 = t100 * (t387 * t505 + t306 * t979 + t504 * t385 + (-pkin(9) * t1118 - t1016 * t310) * t688 + t831) + t114 * (t306 * t1019 + t659 * t387 - t448 * t504 - t771 * t978 + t747);
t720 = -t1323 * t1019 / 0.2e1 + t1276 * t915 + t1275 * t1216 + t1274 * t1215 + ((-t1114 * t1370 + t1257 * t689) * qJD(5) + ((-qJD(5) * t1371 + t1277) * t688 + t1279) * t715 + (t1289 * t558 + t1290 * t557) * t659 + (-t1294 * t558 - t1296 * t557) * t505 + (t1293 * t558 + t1295 * t557) * t504) * t1214 + (qJD(1) * t1248 - t1328 * t715 + t1329 * t712) * t1213 + (qJD(1) * t1249 - t1330 * t715 + t1331 * t712) * t1212 + ((-t1113 * t1369 + t1256 * t689) * qJD(5) + ((-qJD(5) * t1368 + t1277) * t688 + t1279) * t712 + (-t1289 * t560 - t1290 * t559) * t659 + (t1294 * t560 + t1296 * t559) * t505 + (-t1293 * t560 - t1295 * t559) * t504) * t1211 + (t1304 * t712 - t1305 * t715) * t1210 + (t1302 * t712 - t1303 * t715) * t1209 + (t1278 * t712 + t1280 * t715) * t1208 + (t1320 * t715 + t1319 * t712 + (t1302 * t715 + t1303 * t712) * qJD(1)) * t1207 + (t1318 * t715 + t1317 * t712 + (t1304 * t715 + t1305 * t712) * qJD(1)) * t1206 + (-t1278 * t715 + t1280 * t712) * t1205 + (((-t1289 * t710 + t1290 * t713 + t1377) * t659 + (t1294 * t710 - t1296 * t713 - t1364) * t505 + (-t1293 * t710 + t1295 * t713 + t1399) * t504 + t1255 * qJD(5)) * t689 + (-qJD(5) * t1250 + t1363) * t688) * t1204 + (qJD(1) * t1250 - t1326 * t715 + t1327 * t712) * t1203 + (-t1355 * t689 + t1359 * t688) * t1196 + (t1307 * t715 + t1306 * t712 + (t1300 * t715 + t1301 * t712) * qJD(1)) * t1195 + (t1315 + t1324) * t970 + (t1314 + t1325) * t969 + (qJD(1) * t1313 + t1302 * t599 + t1303 * t598 + t1319 * t620 + t1320 * t621 + t1333) * t1201 + (qJD(1) * t1312 + t1304 * t599 + t1305 * t598 + t1317 * t620 + t1318 * t621 + t1334) * t1200 + (t1324 * t712 + t1325 * t715) * t1020 / 0.2e1;
t719 = t70 * (-t1015 * t688 + t1085 * t979 - t1291 * t505 - t1292 * t504 - t563 * t712 + t831) + t95 * (t1019 * t1085 + t1063 * t978 + t1288 * t504 - t1291 * t659 - t712 * t635 + t747) + (-t94 * t955 + (-t712 * t959 - t715 * t960) * t688) * qJD(5);
t608 = t905 * qJD(2);
t549 = t902 * t689;
t456 = t621 * t596;
t446 = -t989 + t1000;
t366 = (t621 * t715 + t1118) * t688;
t360 = rSges(6,1) * t559 + rSges(6,2) * t560;
t356 = -rSges(6,1) * t557 - rSges(6,2) * t558;
t338 = -qJD(2) * t573 + (t715 * t905 + t698) * qJD(1);
t337 = -t1023 * t1191 + (-t1025 * t714 - t982) * rSges(3,1) + t1036;
t317 = -t715 * t840 + t1117;
t311 = t317 * qJD(1);
t196 = -t608 * t1023 + (-t338 - t612 + t983) * qJD(1);
t195 = -t608 * t687 + t577 + (t337 - t980) * qJD(1);
t149 = -t1228 * t715 + t735 * t712;
t148 = t1228 * t712 + t735 * t715;
t147 = -qJD(2) * t843 + t333 * t714 + t335 * t711;
t146 = -qJD(2) * t844 + t334 * t714 + t336 * t711;
t144 = t896 * t621 + (t479 + t993) * qJD(1) + t818;
t135 = -t496 * t621 + t594 * t598 + (-t278 + t1080) * qJD(1) + t839;
t134 = qJD(1) * t277 - t496 * t620 - t594 * t599 + t755;
t116 = t311 + t799;
t115 = t800 + t1012;
t85 = t277 * t621 + t278 * t620 + t476 * t599 - t477 * t598 + t785;
t55 = t279 * t620 - t479 * t599 + (t244 + t280) * t621 - t1056 * t598 + t758;
t28 = t167 * t504 + t169 * t505 - t306 * t367 + t310 * t368 + t728;
t19 = qJD(6) * t795 + t1083 * t368 - t1085 * t367 + t1097 * t505 + t1098 * t504 + t728;
t1 = [(t1306 + t1313) * t1207 + (-t1307 + t1312 + t1314) * t1206 + t1321 * t1213 + (t1322 + t1325) * t1212 + ((t1297 * t710 + t1299 * t713) * t689 - t1364 * t688 + t1366) * t504 * t1204 + (-t1012 + ((t715 * t937 - t1068 + t200) * t715 + (t712 * t937 + t199 + t944) * t712) * qJD(2) + t115) * t968 + (t196 * (t712 * t971 + t1034 + t704) + t195 * t1051 + t313 * (t684 + t1036) + (t1131 * t631 - t1129) * qJD(2) + ((-pkin(1) - t905) * t1130 + (t312 * (-rSges(3,3) - pkin(7)) + t313 * t971) * t712) * qJD(1) - (-qJD(1) * t508 - t312 - t614 - t980) * t313) * m(3) + t1256 * t1216 + t1257 * t1215 + (t216 + t1301) * t1210 + (t97 * (t957 + t1037) + t144 * t923 + t96 * (t478 + t990) + t145 * (-pkin(3) * t1002 + t817 + t924) + (t97 * (-t595 + t1189) - t96 * t716 + t144 * (-t1190 + (-rSges(5,3) - qJ(4)) * t689) * t706) * t712 + ((-t144 * rSges(5,1) + t145 * t781) * t712 + (t144 * (t781 + t1189) - t145 * t716) * t715) * qJD(1) - (qJD(1) * t479 + t1043 * t621 - t144 + t778) * t145) * m(5) - (t146 + t149 + t116) * t1023 / 0.2e1 + (t147 + t148) * t687 / 0.2e1 + ((t1271 * t715 + t1269 + t1302 - t1376) * t621 + ((t1393 + t1411) * t715 + t1271 * t712 + t1303 + t1397) * t620 + t1315 + t1343) * t1208 + ((t315 + t317) * t715 + (t314 + t316) * t712) * t1011 / 0.2e1 + (t46 * (t837 + t903) + t113 * (t894 - t904) + t45 * (t893 + t306) + t114 * (t794 + t1094) + (t46 * (-t595 - t1184) - t45 * t716) * t712 + (-t113 * t1005 + (rSges(6,3) * t1144 + t114 * t933) * t688) * t706 - (-t113 + t1237 + t1268) * t114) * m(6) + t1239 + t1183 / 0.2e1 + (t311 + ((t198 - t430 + (t499 + t1123) * t715 + t1069) * t715 + t1068 * t712) * qJD(2)) * t965 + (-qJD(2) * t840 + t603 * t714 + t604 * t711 + m(6) * (t113 * t779 + t114 * t777) + m(7) * (t95 * t777 + t94 * t779) + (m(6) * (t1008 * t1143 + t113 * t933) + m(7) * (t1009 * t1179 + t934 * t94) - t1390) * t689 + t1391 * t688) * qJD(1) + t1137 / 0.2e1 + t1138 / 0.2e1 + t1135 / 0.2e1 + t1136 / 0.2e1 + (t27 * (t837 + t1346) + t26 * (t893 + t1085) + (t27 * (-t595 - t1187) - t26 * t716) * t712 + (t1115 * t934 - t1237 + t1344 - t1351 + t794) * t95 + (t894 + (rSges(7,2) * t1114 - t1005) * t706 + t95 + t1348) * t94) * m(7) + t1325 * t1211 + t1181 / 0.2e1 - t1182 / 0.2e1 + (t135 * (-t476 + t1037) + t185 * t1035 + t134 * (t477 + t936) + t186 * (-t932 + t1042) + (t1132 * t594 - t186 * t552) * t706 + ((-t185 * rSges(4,3) + t186 * (-t679 - t1192)) * t712 + (t185 * (-t597 - t679) - t186 * t716) * t715) * qJD(1) - (-qJD(1) * t476 + t1254 - t185 + t798) * t186) * m(4) - t1180 / 0.2e1 - t598 * t219 / 0.2e1 + (t1300 + t1378) * t1209 + (((t464 * t715 + t465 * t712) * t688 + (t470 * t712 - t471 * t715) * t689 + t1304 + t1373 + t1375) * t621 + (-t466 * t1111 + (t464 * t712 - t465 * t715) * t688 + t469 * t1114 + (-t470 * t715 - t471 * t712) * t689 - t1269 + t1305 - t1372) * t620 + t1347) * t1205; t720 + (-t146 * t715 + t147 * t712 + (t314 * t712 + t315 * t715) * qJD(1)) * t1195 + ((-t1116 * t687 + t1026) * t712 + (t762 + (t1117 * t712 + t729) * qJD(2)) * t715) * t968 + ((t1040 * t714 + t1041 * t711) * qJD(1) + (t1224 * t714 + t782 * t711) * qJD(2)) * t1196 + ((-t1023 * t1117 - t1026) * t715 + (t762 + (t1116 * t715 + t729) * qJD(2)) * t712) * t965 + (qJD(1) * t148 + (-(t1231 * t712 + t737 * t715) * t715 + (t1232 * t712 + t736 * t715) * t712 + (t199 * t712 + t200 * t715) * qJD(1)) * t1267) * t1201 + (qJD(1) * t149 + (-(-t1231 * t715 + t737 * t712) * t715 + (-t1232 * t715 + t736 * t712) * t712 + (t197 * t712 + t198 * t715) * qJD(1)) * t1267) * t1200 + (t115 + t800) * t970 + (t799 + t116) * t969 + (t94 * t929 + t19 * (t783 + t1072) + t70 * (t756 + t996) + (t1258 * t793 + t769 * t94) * t715 + (t26 * t793 + t95 * t769 + (-t459 + t930) * t1146) * t712 - t732 * t94 - (-t95 * t984 + ((-t715 * t94 - t1179) * t714 + t70 * t956) * qJD(2)) * pkin(2) - t719) * m(7) + (t28 * (t827 + t1072) + t100 * (t775 + t996) + (t114 * t776 + (-t459 + t997) * t1033) * t712 - (-t114 * t984 + (t100 * t956 - t714 * t870) * qJD(2)) * pkin(2) - t721 + t1223 * (t771 + t832) + (t1019 * t310 + t715 * t776 - t771 * t979 + t1245) * t113) * m(6) + (-t144 * (-qJD(1) * t543 - t456 + t907) - (-t145 * t984 + ((-t144 * t715 - t145 * t712) * t714 + t125 * t956) * qJD(2)) * pkin(2) - t733 + t144 * t1047 + t55 * (t926 + t1072) + t125 * (t829 + t996) + (t1259 * t913 + t144 * t833) * t715 + (t96 * t913 + t145 * t833 + (-t459 - t1056) * t1031) * t712) * m(5) + (-(-t186 * t984 + (t176 * t956 + t714 * t857) * qJD(2)) * pkin(2) + t85 * (t1070 + t1072) + t176 * (t996 + t998) + (t185 * t908 + (qJD(1) * t186 + t135) * t962) * t715 + (t134 * t962 + t186 * t908 + (t185 * t594 + t176 * (-t459 - t477)) * qJD(1)) * t712 + t1238) * m(4) + (0.2e1 * t285 * (t337 * t715 + t338 * t712 + (t508 * t715 - t509 * t712) * qJD(1)) + t851 * t608 + (-t195 * t712 - t196 * t715 + (-t313 * t715 + t1131) * qJD(1)) * t631 - (t312 * t573 - t1129) * qJD(1) - (t285 * (-t573 * t712 - t574 * t715) + t851 * t905) * qJD(2)) * m(3); t720 + (t19 * t783 + t70 * t756 + (t1146 * t930 + t822 * t95) * t712 - t719 + (t1258 * t715 + t1265) * (t961 - t1063) + (t715 * t822 - t732 + t929) * t94) * m(7) + (t28 * t827 + t100 * t775 + (t1033 * t997 + t114 * t838) * t712 - t721 + t1223 * (t771 + t961) + (t838 * t715 - (t1114 * t771 - t310 * t689) * qJD(5) + t1245) * t113) * m(6) + (-t733 + t55 * t926 + t125 * t829 + (-t1031 * t1056 + t1071 * t145) * t712 + (t1259 * t715 + t1263) * t1043 + (t456 - (t542 - t543) * qJD(1) - t1038 + t1047 + t1071 * t715) * t144) * m(5) + (t85 * t1070 + t176 * (-t1025 * t477 + t998) + t857 * t496 + (-t134 * t712 - t135 * t715 + (-t186 * t715 + t1132) * qJD(1)) * t594 + t1238) * m(4); -m(5) * (t125 * t366 + t144 * t446 + t145 * t445) - m(6) * (t100 * t366 + t113 * t446 + t114 * t445) - m(7) * (t366 * t70 + t445 * t95 + t446 * t94) + 0.2e1 * ((t1108 * t145 + t144 * t621 - t55) * t1220 + (t1108 * t114 + t113 * t621 - t28) * t1219 + (t1108 * t95 + t621 * t94 - t19) * t1218) * t689 + 0.2e1 * ((t1024 * t145 - t1025 * t144 + t125 * t706 + t715 * t97 + t1263) * t1220 + (t100 * t706 + t1024 * t114 - t1025 * t113 + t46 * t715 + t1264) * t1219 + (t1024 * t95 - t1025 * t94 + t27 * t715 + t70 * t706 + t1265) * t1218) * t688; t1323 * t974 + t1334 * t1111 / 0.2e1 + t1333 * t1109 / 0.2e1 + (t1250 * t689 + t1255 * t688) * t915 + (t1249 * t689 + t1256 * t688) * t1216 + (t1248 * t689 + t1257 * t688) * t1215 + (t1251 * t715 - t1252 * t558 - t1350 * t557) * t1214 + ((-qJD(1) * t1274 + t1257 * t706 + t1328 * t712 + t1329 * t715) * t689 + (-t1248 * t706 + t1321) * t688) * t1213 + ((-qJD(1) * t1275 + t1256 * t706 + t1330 * t712 + t1331 * t715) * t689 + (-t1249 * t706 + t1322) * t688) * t1212 + (t1251 * t712 + t1252 * t560 + t1350 * t559) * t1211 + ((t1252 * t710 - t1350 * t713) * t689 + t1345 * t688) * t1204 + ((-qJD(1) * t1276 + t1255 * t706 + t1326 * t712 + t1327 * t715) * t689 + (-t1250 * t706 + t1349) * t688) * t1203 + (t1137 + t1138 - t1182 + t1183 + t1135 + t1136 - t1180 + t1181 + t1239) * t688 / 0.2e1 + (-(-t689 * t70 * t710 + t558 * t94 - t560 * t95) * qJD(6) - (-t1081 * t94 + t1082 * t95) * t659 - (t1049 * t94 + t1082 * t70) * t505 - (t1049 * t95 + t1081 * t70) * t504 + (-t27 * t1083 - t94 * t1098 + t26 * t1085 + t95 * t1097 + (-t1270 * t712 + t784 * t715) * t706) * t688 + ((-t1083 * t94 + t1085 * t95) * t706 + (qJD(1) * t1270 - t26 * t1063 + t19 * t1083 - t95 * t1096 + t70 * t1098) * t715 + (qJD(1) * t784 + t1063 * t27 - t1085 * t19 + t1096 * t94 - t1097 * t70) * t712) * t689) * m(7) + (-t113 * (-t360 * t659 - t505 * t549) - t114 * (t356 * t659 - t504 * t549) - t100 * (t356 * t505 + t360 * t504) + (-t113 * t167 + t114 * t169 + t45 * t306 - t46 * t310 + (t100 * t852 - (t114 * t715 - t1144) * t771) * t706) * t688 + (t113 * (t236 * t712 - t310 * t706) + t114 * (-t236 * t715 + t306 * t706) - t28 * t852 + t100 * (-t1024 * t306 - t1025 * t310 + t167 * t715 - t169 * t712) - (qJD(1) * t870 - t45 * t715 + t46 * t712) * t771) * t689) * m(6) + t1324 * (t689 * t969 + t712 * t975) + t1325 * (-t987 / 0.2e1 + t715 * t975); (t19 * t1110 - t26 * t559 + t27 * t557 + (t1110 * t504 - t557 * t659 + t257) * t95 + (t1110 * t505 - t559 * t659 + t259) * t94 + (t504 * t559 - t505 * t557 + t795) * t70) * m(7);];
tauc  = t1(:);
