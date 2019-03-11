% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:15
% EndTime: 2019-03-09 09:56:52
% DurationCPUTime: 142.32s
% Computational Cost: add. (44167->1520), mult. (70651->1879), div. (0->0), fcn. (66796->8), ass. (0->695)
t1296 = Icges(6,4) - Icges(5,5) - Icges(7,5);
t1295 = Icges(7,4) + Icges(6,5) - Icges(5,6);
t687 = pkin(9) + qJ(4);
t676 = cos(t687);
t1346 = (-Icges(5,4) - Icges(6,6)) * t676;
t675 = sin(t687);
t1345 = (Icges(5,4) - Icges(7,6)) * t675;
t1344 = -Icges(5,1) - Icges(6,2);
t1341 = Icges(5,2) + Icges(7,2);
t1340 = -t1345 + (Icges(5,1) + Icges(7,3)) * t676;
t1339 = -t1346 + (-Icges(5,2) - Icges(6,3)) * t675;
t691 = sin(qJ(2));
t693 = cos(qJ(2));
t1338 = t1295 * t693 + t1339 * t691;
t1337 = t1296 * t693 + t1340 * t691;
t694 = cos(qJ(1));
t1062 = t691 * t694;
t1335 = t1295 * t675 - t1296 * t676;
t1336 = Icges(6,1) + Icges(7,1) + Icges(5,3);
t1283 = t1335 * t691 - t1336 * t693;
t1067 = t676 * t691;
t1068 = t675 * t691;
t634 = Icges(6,6) * t1068;
t1314 = -Icges(6,2) * t1067 - t1337 + t634;
t633 = Icges(7,6) * t1067;
t1332 = Icges(7,2) * t1068 - t1338 + t633;
t1059 = t694 * t675;
t692 = sin(qJ(1));
t522 = t1059 * t693 - t692 * t676;
t1060 = t693 * t694;
t523 = t1060 * t676 + t675 * t692;
t1189 = t1062 * t1283 - t1314 * t523 + t1332 * t522;
t256 = Icges(7,1) * t1062 + Icges(7,4) * t522 + Icges(7,5) * t523;
t259 = Icges(6,1) * t1062 - Icges(6,4) * t523 + Icges(6,5) * t522;
t262 = Icges(5,5) * t523 - Icges(5,6) * t522 + Icges(5,3) * t1062;
t1284 = t256 + t259 + t262;
t1092 = Icges(6,6) * t523;
t247 = Icges(6,5) * t1062 + Icges(6,3) * t522 - t1092;
t497 = Icges(7,6) * t523;
t250 = Icges(7,4) * t1062 + Icges(7,2) * t522 + t497;
t1106 = Icges(5,4) * t523;
t265 = -Icges(5,2) * t522 + Icges(5,6) * t1062 + t1106;
t1304 = t247 + t250 - t265;
t1089 = Icges(7,6) * t522;
t244 = Icges(7,5) * t1062 + Icges(7,3) * t523 + t1089;
t500 = Icges(6,6) * t522;
t253 = Icges(6,4) * t1062 - Icges(6,2) * t523 + t500;
t503 = Icges(5,4) * t522;
t268 = Icges(5,1) * t523 + Icges(5,5) * t1062 - t503;
t1306 = t244 - t253 + t268;
t1261 = t1062 * t1284 + t1304 * t522 + t1306 * t523;
t1063 = t691 * t692;
t1061 = t692 * t693;
t520 = t1061 * t675 + t676 * t694;
t521 = t1061 * t676 - t1059;
t254 = Icges(7,1) * t1063 + Icges(7,4) * t520 + Icges(7,5) * t521;
t257 = Icges(6,1) * t1063 - Icges(6,4) * t521 + Icges(6,5) * t520;
t260 = Icges(5,5) * t521 - Icges(5,6) * t520 + Icges(5,3) * t1063;
t1303 = t254 + t257 + t260;
t499 = Icges(6,6) * t521;
t245 = Icges(6,5) * t1063 + Icges(6,3) * t520 - t499;
t496 = Icges(7,6) * t521;
t249 = -Icges(7,4) * t1063 - Icges(7,2) * t520 - t496;
t502 = Icges(5,4) * t521;
t263 = -Icges(5,2) * t520 + Icges(5,6) * t1063 + t502;
t1305 = t245 - t249 - t263;
t495 = Icges(7,6) * t520;
t242 = Icges(7,5) * t1063 + Icges(7,3) * t521 + t495;
t498 = Icges(6,6) * t520;
t252 = -Icges(6,4) * t1063 + Icges(6,2) * t521 - t498;
t501 = Icges(5,4) * t520;
t267 = -Icges(5,1) * t521 - Icges(5,5) * t1063 + t501;
t1307 = t242 + t252 - t267;
t1262 = t1303 * t1062 + t1305 * t522 + t1307 * t523;
t991 = qJD(4) * t694;
t999 = qJD(2) * t692;
t590 = t691 * t991 + t999;
t993 = qJD(4) * t692;
t997 = qJD(2) * t694;
t591 = -t691 * t993 + t997;
t992 = qJD(4) * t693;
t667 = qJD(1) - t992;
t1226 = t1189 * t667 + t1261 * t590 - t1262 * t591;
t1190 = t1063 * t1283 - t1314 * t521 + t1332 * t520;
t1263 = t1063 * t1284 + t1304 * t520 + t1306 * t521;
t1264 = t1303 * t1063 + t1305 * t520 + t1307 * t521;
t1227 = t1190 * t667 + t1263 * t590 - t1264 * t591;
t1331 = (t1341 * t676 + t1345) * t691;
t1330 = (t1344 * t675 + t1346) * t691;
t1325 = t1335 * t693 + t1336 * t691;
t1087 = Icges(7,6) * t676;
t836 = Icges(7,2) * t675 + t1087;
t1220 = (t836 - t1339) * t693 + t1295 * t691;
t1091 = Icges(6,6) * t675;
t837 = Icges(6,2) * t676 - t1091;
t1219 = (t837 + t1340) * t693 - t1296 * t691;
t1322 = (t1295 * t676 + t1296 * t675) * t691;
t688 = sin(pkin(9));
t689 = cos(pkin(9));
t564 = t1061 * t688 + t689 * t694;
t1058 = t694 * t688;
t565 = t1061 * t689 - t1058;
t361 = Icges(4,5) * t565 - Icges(4,6) * t564 + Icges(4,3) * t1063;
t1086 = Icges(3,3) * t694;
t524 = Icges(3,5) * t1061 - Icges(3,6) * t1063 - t1086;
t1095 = Icges(3,6) * t694;
t526 = Icges(3,4) * t1061 - Icges(3,2) * t1063 - t1095;
t1072 = t526 * t691;
t1101 = Icges(3,5) * t694;
t660 = Icges(3,4) * t1063;
t528 = Icges(3,1) * t1061 - t1101 - t660;
t818 = -t528 * t693 + t1072;
t364 = Icges(4,4) * t565 - Icges(4,2) * t564 + Icges(4,6) * t1063;
t367 = Icges(4,1) * t565 - Icges(4,4) * t564 + Icges(4,5) * t1063;
t825 = -t364 * t564 + t367 * t565;
t1319 = t1063 * t361 - t524 * t694 - t692 * t818 + t825;
t566 = -t1058 * t693 + t689 * t692;
t1066 = t688 * t692;
t567 = t1060 * t689 + t1066;
t363 = Icges(4,5) * t567 + Icges(4,6) * t566 + Icges(4,3) * t1062;
t366 = Icges(4,4) * t567 + Icges(4,2) * t566 + Icges(4,6) * t1062;
t369 = Icges(4,1) * t567 + Icges(4,4) * t566 + Icges(4,5) * t1062;
t102 = t363 * t1063 - t564 * t366 + t565 * t369;
t681 = Icges(3,4) * t693;
t844 = -Icges(3,2) * t691 + t681;
t527 = Icges(3,6) * t692 + t694 * t844;
t1107 = Icges(3,4) * t691;
t616 = Icges(3,1) * t693 - t1107;
t529 = Icges(3,5) * t692 + t616 * t694;
t446 = t529 * t1061;
t612 = Icges(3,5) * t693 - Icges(3,6) * t691;
t525 = Icges(3,3) * t692 + t612 * t694;
t911 = t525 * t694 - t446;
t184 = -t1063 * t527 - t911;
t1318 = t102 + t184;
t611 = Icges(3,5) * t691 + Icges(3,6) * t693;
t1069 = t611 * t694;
t839 = Icges(4,5) * t689 - Icges(4,6) * t688;
t512 = -Icges(4,3) * t693 + t691 * t839;
t843 = Icges(4,4) * t689 - Icges(4,2) * t688;
t514 = -Icges(4,6) * t693 + t691 * t843;
t846 = Icges(4,1) * t689 - Icges(4,4) * t688;
t516 = -Icges(4,5) * t693 + t691 * t846;
t613 = Icges(3,2) * t693 + t1107;
t615 = Icges(3,1) * t691 + t681;
t816 = t613 * t691 - t615 * t693;
t1317 = t1063 * t512 - t514 * t564 + t516 * t565 - t692 * t816 - t1069;
t1070 = t611 * t692;
t1316 = t1062 * t512 + t514 * t566 + t516 * t567 - t694 * t816 + t1070;
t944 = t676 * t991;
t946 = t675 * t993;
t954 = t691 * t997;
t237 = qJD(1) * t520 + t675 * t954 - t693 * t944 - t946;
t1001 = qJD(1) * t694;
t768 = t1001 * t675 + t676 * t993;
t1002 = qJD(1) * t693;
t769 = t1002 * t692 + t954;
t945 = t675 * t991;
t238 = t676 * t769 + t693 * t945 - t768;
t952 = t693 * t997;
t1003 = qJD(1) * t692;
t958 = t691 * t1003;
t571 = t952 - t958;
t120 = Icges(7,5) * t571 - Icges(7,6) * t237 - Icges(7,3) * t238;
t126 = Icges(6,4) * t571 + Icges(6,2) * t238 - Icges(6,6) * t237;
t136 = -Icges(5,1) * t238 + Icges(5,4) * t237 + Icges(5,5) * t571;
t1313 = t120 - t126 + t136;
t955 = t691 * t999;
t239 = -t1003 * t676 - t675 * t955 + t693 * t768 - t945;
t240 = qJD(1) * t523 - t676 * t955 - t693 * t946 - t944;
t998 = qJD(2) * t693;
t953 = t692 * t998;
t570 = t691 * t1001 + t953;
t121 = Icges(7,5) * t570 + Icges(7,6) * t239 + Icges(7,3) * t240;
t127 = Icges(6,4) * t570 - Icges(6,2) * t240 + Icges(6,6) * t239;
t137 = Icges(5,1) * t240 - Icges(5,4) * t239 + Icges(5,5) * t570;
t1312 = t121 - t127 + t137;
t122 = Icges(6,5) * t571 + Icges(6,6) * t238 - Icges(6,3) * t237;
t124 = Icges(7,4) * t571 - Icges(7,2) * t237 - Icges(7,6) * t238;
t134 = -Icges(5,4) * t238 + Icges(5,2) * t237 + Icges(5,6) * t571;
t1311 = t122 + t124 - t134;
t123 = Icges(6,5) * t570 - Icges(6,6) * t240 + Icges(6,3) * t239;
t125 = Icges(7,4) * t570 + Icges(7,2) * t239 + Icges(7,6) * t240;
t135 = Icges(5,4) * t240 - Icges(5,2) * t239 + Icges(5,6) * t570;
t1310 = t123 + t125 - t135;
t128 = Icges(7,1) * t571 - Icges(7,4) * t237 - Icges(7,5) * t238;
t130 = Icges(6,1) * t571 + Icges(6,4) * t238 - Icges(6,5) * t237;
t132 = -Icges(5,5) * t238 + Icges(5,6) * t237 + Icges(5,3) * t571;
t1309 = t128 + t130 + t132;
t129 = Icges(7,1) * t570 + Icges(7,4) * t239 + Icges(7,5) * t240;
t131 = Icges(6,1) * t570 - Icges(6,4) * t240 + Icges(6,5) * t239;
t133 = Icges(5,5) * t240 - Icges(5,6) * t239 + Icges(5,3) * t570;
t1308 = t129 + t131 + t133;
t1302 = qJD(2) * t1325 + qJD(4) * t1322;
t994 = qJD(4) * t691;
t1301 = (Icges(6,3) * t676 + t1091) * t994 + t1331 * qJD(4) + t1220 * qJD(2);
t1300 = (-Icges(7,3) * t675 + t1087) * t994 + t1330 * qJD(4) + t1219 * qJD(2);
t1297 = -t1314 * t676 + t1332 * t675;
t1254 = rSges(7,1) + pkin(5);
t1294 = t1316 * qJD(1);
t1209 = t566 * t364 + t567 * t367;
t103 = t361 * t1062 + t1209;
t104 = t363 * t1062 + t566 * t366 + t567 * t369;
t1213 = -t103 * t694 + t104 * t692;
t1033 = -t528 * t1060 - t692 * t524;
t185 = -t1062 * t526 - t1033;
t1032 = t529 * t1060 + t692 * t525;
t186 = -t1062 * t527 + t1032;
t1291 = (-t185 * t694 + t186 * t692 + t1213) * qJD(2);
t1290 = (t1318 * t692 - t1319 * t694) * qJD(2);
t776 = qJD(2) * t611;
t1289 = (-t692 * t776 + (t525 + t818) * qJD(1)) * t694;
t1288 = t1317 * qJD(1);
t1233 = t1062 * t1308 + t1303 * t571 - t1305 * t237 - t1307 * t238 + t1310 * t522 + t1312 * t523;
t1232 = t1062 * t1309 + t1284 * t571 - t1304 * t237 - t1306 * t238 + t1311 * t522 + t1313 * t523;
t1231 = t1063 * t1308 + t1303 * t570 + t1305 * t239 + t1307 * t240 + t1310 * t520 + t1312 * t521;
t1230 = t1063 * t1309 + t1284 * t570 + t1304 * t239 + t1306 * t240 + t1311 * t520 + t1313 * t521;
t1267 = t1063 * t1302 + t1283 * t570 + t1300 * t521 + t1301 * t520 - t1314 * t240 + t1332 * t239;
t1266 = t1062 * t1302 + t1283 * t571 + t1300 * t523 + t1301 * t522 + t1314 * t238 - t1332 * t237;
t1250 = t260 * t693;
t829 = -t263 * t675 - t267 * t676;
t95 = t691 * t829 - t1250;
t1252 = t254 * t693;
t833 = t242 * t676 - t249 * t675;
t97 = t691 * t833 - t1252;
t1251 = t257 * t693;
t831 = t245 * t675 + t252 * t676;
t99 = t691 * t831 - t1251;
t1287 = t95 + t97 + t99;
t830 = t247 * t675 - t253 * t676;
t100 = -t259 * t693 + t691 * t830;
t828 = -t265 * t675 + t268 * t676;
t96 = -t262 * t693 + t691 * t828;
t832 = t244 * t676 + t250 * t675;
t98 = -t256 * t693 + t691 * t832;
t1286 = t100 + t96 + t98;
t1188 = -t1283 * t693 + t1297 * t691;
t1282 = -t691 * t836 + t1338;
t1281 = t691 * t837 + t1337;
t1280 = (-t1297 + t1325) * t667 + (t1283 * t692 + t829 + t831 + t833) * t591 + (-t1283 * t694 - t828 - t830 - t832) * t590;
t1277 = t1288 + t1290;
t1276 = t1291 + t1294;
t1253 = rSges(7,3) + qJ(6);
t1208 = t364 * t688 - t367 * t689;
t432 = qJD(1) * t566 + t688 * t955;
t433 = qJD(1) * t567 - t689 * t955;
t176 = Icges(4,5) * t433 + Icges(4,6) * t432 + Icges(4,3) * t570;
t178 = Icges(4,4) * t433 + Icges(4,2) * t432 + Icges(4,6) * t570;
t180 = Icges(4,1) * t433 + Icges(4,4) * t432 + Icges(4,5) * t570;
t777 = qJD(2) * t613;
t778 = qJD(2) * t615;
t1275 = (-qJD(1) * t527 + t692 * t777 + t176) * t693 + (-qJD(1) * t529 + t178 * t688 - t180 * t689 + t692 * t778) * t691 + (t1208 * t693 - t361 * t691 + t818) * qJD(2);
t430 = qJD(1) * t564 + t688 * t954;
t431 = -qJD(1) * t565 - t689 * t954;
t175 = Icges(4,5) * t431 + Icges(4,6) * t430 + Icges(4,3) * t571;
t177 = Icges(4,4) * t431 + Icges(4,2) * t430 + Icges(4,6) * t571;
t179 = Icges(4,1) * t431 + Icges(4,4) * t430 + Icges(4,5) * t571;
t1071 = t527 * t691;
t817 = -t529 * t693 + t1071;
t823 = -t366 * t688 + t369 * t689;
t1274 = (-t694 * t777 + (-t692 * t844 + t1095) * qJD(1) - t175) * t693 + (-t177 * t688 + t179 * t689 - t694 * t778 + (-t616 * t692 + t1101) * qJD(1)) * t691 + (t363 * t691 + t693 * t823 - t817) * qJD(2);
t1177 = qJD(1) * t816 + t612 * qJD(2);
t513 = Icges(4,3) * t691 + t693 * t839;
t452 = t513 * qJD(2);
t515 = Icges(4,6) * t691 + t693 * t843;
t453 = t515 * qJD(2);
t517 = Icges(4,5) * t691 + t693 * t846;
t454 = t517 * qJD(2);
t593 = t844 * qJD(2);
t594 = t616 * qJD(2);
t703 = qJD(1) * t611 - t593 * t691 + t594 * t693 + (-t613 * t693 - t615 * t691) * qJD(2);
t1273 = t1062 * t452 + t1177 * t692 + t430 * t514 + t431 * t516 + t453 * t566 + t454 * t567 + t512 * t571 + t703 * t694;
t1272 = t1063 * t452 - t1177 * t694 + t432 * t514 + t433 * t516 - t453 * t564 + t454 * t565 + t512 * t570 + t703 * t692;
t1271 = t103 + t185;
t1270 = (-t361 + t526) * t693 + (-t1208 + t528) * t691;
t1269 = (-t363 + t527) * t693 + (t529 + t823) * t691;
t873 = rSges(7,2) * t675 + rSges(7,3) * t676;
t1248 = -qJ(6) * t1067 + t1254 * t693 - t691 * t873;
t1007 = t694 * pkin(1) + t692 * pkin(7);
t1234 = qJD(1) * t1007;
t628 = rSges(3,1) * t691 + rSges(3,2) * t693;
t1268 = t628 * t999 - t1234;
t1265 = (qJD(2) * t1297 - t1302) * t693 + (t1300 * t676 + t1301 * t675 + (t1314 * t675 + t1332 * t676) * qJD(4) + t1283 * qJD(2)) * t691;
t1260 = Icges(7,3) - t1344;
t1259 = Icges(6,3) + t1341;
t1256 = -t1322 * t667 + (t1295 * t521 + t1296 * t520) * t591 + (-t1295 * t523 - t1296 * t522) * t590;
t1045 = t522 * rSges(7,2) + t1062 * t1254 + t1253 * t523;
t872 = pkin(4) * t676 + qJ(5) * t675;
t1207 = t872 * t691;
t358 = t523 * pkin(4) + qJ(5) * t522;
t669 = pkin(2) * t1060;
t583 = qJ(3) * t1062 + t669;
t656 = pkin(3) * t1066;
t671 = pkin(3) * t689 + pkin(2);
t1015 = t671 * t1060 + t656;
t690 = -pkin(8) - qJ(3);
t979 = t690 * t1062;
t814 = -t979 + t1015;
t374 = t814 - t583;
t1057 = qJ(3) + t690;
t1147 = pkin(2) - t671;
t938 = t1147 * t691;
t757 = -t1057 * t693 + t938;
t627 = pkin(2) * t691 - qJ(3) * t693;
t677 = qJD(3) * t691;
t652 = t692 * t677;
t901 = qJD(1) * t583 - t627 * t999 + t1234 + t652;
t809 = qJD(1) * t374 + t757 * t999 + t901;
t990 = qJD(5) * t520;
t1245 = t1207 * t590 - t667 * t358 - t809 - t990;
t488 = qJD(6) * t521;
t57 = t1045 * t667 + t1248 * t590 - t1245 + t488;
t1255 = 0.2e1 * qJD(2);
t1012 = t615 + t844;
t1013 = -t613 + t616;
t819 = -t514 * t688 + t516 * t689;
t1210 = qJD(2) * (-t1208 * t694 - t692 * t823) + (t513 - t819) * qJD(1);
t1247 = t512 * t1002 + t1210 * t691 + (-t1012 * t691 + t1013 * t693) * qJD(1);
t1174 = (-t613 * t694 + t529) * t692 - (-Icges(3,2) * t1061 + t528 - t660) * t694;
t753 = t526 * t694 - t527 * t692;
t1246 = -t1174 * t691 + (-t361 * t694 + t363 * t692 + t753) * t693;
t587 = t627 * t1003;
t1034 = -t1003 * t757 + t587;
t455 = t692 * t1207;
t478 = t1207 * t1003;
t553 = t872 * t693;
t390 = t757 * t692;
t578 = t627 * t692;
t995 = qJD(3) * t694;
t909 = qJD(1) * t578 + t693 * t995;
t782 = -qJD(1) * t390 + t909;
t1244 = -t455 * t992 + t591 * t553 + t1034 + t478 - t782;
t1083 = qJ(3) * t691;
t939 = t1147 * t693;
t1237 = -t939 - t1083;
t982 = qJD(2) * qJD(4);
t933 = t693 * t982;
t457 = qJD(1) * t590 + t692 * t933;
t458 = qJD(1) * t591 + t694 * t933;
t934 = t691 * t982;
t1236 = t1189 * t934 + t1232 * t590 - t1233 * t591 + t1261 * t458 + t1262 * t457 + t1266 * t667;
t1235 = t1190 * t934 + t1230 * t590 - t1231 * t591 + t1263 * t458 + t1264 * t457 + t1267 * t667;
t26 = (qJD(2) * t829 - t133) * t693 + (qJD(2) * t260 - t135 * t675 + t137 * t676 + (-t263 * t676 + t267 * t675) * qJD(4)) * t691;
t28 = (qJD(2) * t833 - t129) * t693 + (qJD(2) * t254 + t121 * t676 + t125 * t675 + (-t242 * t675 - t249 * t676) * qJD(4)) * t691;
t30 = (qJD(2) * t831 - t131) * t693 + (qJD(2) * t257 + t123 * t675 - t127 * t676 + (t245 * t676 - t252 * t675) * qJD(4)) * t691;
t1229 = t26 + t28 + t30;
t27 = (qJD(2) * t828 - t132) * t693 + (qJD(2) * t262 - t134 * t675 + t136 * t676 + (-t265 * t676 - t268 * t675) * qJD(4)) * t691;
t29 = (qJD(2) * t832 - t128) * t693 + (qJD(2) * t256 + t120 * t676 + t124 * t675 + (-t244 * t675 + t250 * t676) * qJD(4)) * t691;
t31 = (qJD(2) * t830 - t130) * t693 + (qJD(2) * t259 + t122 * t675 - t126 * t676 + (t247 * t676 + t253 * t675) * qJD(4)) * t691;
t1228 = t27 + t29 + t31;
t1225 = t1188 * t667 + t1286 * t590 - t1287 * t591;
t1224 = t1282 * t692;
t1223 = t1282 * t694;
t1222 = t1281 * t692;
t1221 = t1281 * t694;
t1218 = t1280 * t691;
t1217 = t1283 * t667 + t1284 * t590 - t1303 * t591;
t1183 = t1286 * t694 + t1287 * t692;
t1216 = t1286 * t692 - t1287 * t694;
t1182 = t1261 * t694 + t1262 * t692;
t1181 = t1263 * t694 + t1264 * t692;
t1215 = t1261 * t692 - t1262 * t694;
t1214 = t1263 * t692 - t1264 * t694;
t878 = rSges(5,1) * t521 - rSges(5,2) * t520;
t280 = rSges(5,3) * t1063 + t878;
t877 = rSges(5,1) * t676 - rSges(5,2) * t675;
t479 = -rSges(5,3) * t693 + t691 * t877;
t1212 = -t280 * t667 - t479 * t591;
t801 = t1068 * t591 + t520 * t667;
t1206 = -t237 + t801;
t800 = t1068 * t590 - t522 * t667;
t1205 = t239 + t800;
t1140 = rSges(7,2) * t520;
t1048 = t1063 * t1254 + t1253 * t521 + t1140;
t1031 = t757 - t627;
t1154 = pkin(2) * t693;
t630 = t1083 + t1154;
t580 = t630 * t692;
t685 = t694 * pkin(7);
t635 = pkin(1) * t692 - t685;
t604 = qJD(1) * t635;
t1204 = -qJD(1) * t580 - t604;
t489 = qJD(6) * t523;
t1203 = -t237 * rSges(7,2) - t1253 * t238 + t1254 * t952 + t489;
t279 = rSges(6,1) * t1062 - t523 * rSges(6,2) + t522 * rSges(6,3);
t874 = rSges(6,2) * t676 - rSges(6,3) * t675;
t744 = rSges(6,1) * t693 + t691 * t874;
t69 = t279 * t667 + t590 * t744 - t1245;
t1202 = t1057 * t691;
t1199 = -rSges(7,2) * t239 - t488;
t494 = t520 * qJ(5);
t352 = pkin(4) * t521 + t494;
t1186 = (Icges(6,3) * t1067 + t1314 + t1331 + t634) * t667 + (-t1259 * t521 + t1307 + t495 - t498 - t501) * t591 + (t1259 * t523 - t1089 - t1306 + t500 + t503) * t590;
t1185 = (-Icges(7,3) * t1068 + t1330 + t1332 + t633) * t667 + (t1260 * t520 - t1305 - t496 + t499 + t502) * t591 + (-t1260 * t522 - t1092 - t1106 + t1304 + t497) * t590;
t1184 = t1256 * t691;
t1180 = t1188 * t934 + t1265 * t667;
t1179 = -t694 * t776 + (-t612 * t692 + t1086 + t817) * qJD(1);
t880 = rSges(4,1) * t689 - rSges(4,2) * t688;
t518 = -rSges(4,3) * t693 + t691 * t880;
t1052 = t239 * qJ(5) + t990;
t1152 = pkin(4) * t240;
t115 = t1052 + t1152;
t988 = qJD(5) * t691;
t623 = t676 * t988;
t989 = qJD(5) * t675;
t624 = t693 * t989;
t1176 = qJD(2) * t624 + qJD(4) * t623 + t590 * t115 + t458 * t352;
t1170 = m(4) / 0.2e1;
t1169 = m(5) / 0.2e1;
t1168 = m(6) / 0.2e1;
t1167 = m(7) / 0.2e1;
t1166 = t457 / 0.2e1;
t1165 = t458 / 0.2e1;
t1164 = -t590 / 0.2e1;
t1163 = t590 / 0.2e1;
t1162 = -t591 / 0.2e1;
t1161 = t591 / 0.2e1;
t1160 = -t667 / 0.2e1;
t1159 = t667 / 0.2e1;
t1153 = pkin(3) * t688;
t1150 = pkin(5) * t691;
t1146 = rSges(3,1) * t693;
t1145 = rSges(6,1) * t691;
t1143 = rSges(7,1) * t691;
t1139 = rSges(4,3) * t691;
t1137 = rSges(5,3) * t691;
t206 = t520 * t590 + t522 * t591;
t986 = qJD(6) * t691;
t621 = t676 * t986;
t622 = t675 * t988;
t601 = t671 * t1061;
t980 = pkin(3) * t1058;
t892 = -t601 + t980;
t373 = (t1154 + t1202) * t692 + t892;
t996 = qJD(3) * t693;
t885 = t580 * t999 + t583 * t997 - t996;
t779 = -t373 * t999 + t374 * t997 + t885;
t751 = t590 * t352 + t622 + t779;
t972 = t358 + t1045;
t44 = t1048 * t590 + t591 * t972 + t621 + t751;
t1135 = t206 * t44;
t1134 = t26 * t591;
t1133 = t27 * t590;
t1132 = t28 * t591;
t1131 = t29 * t590;
t1130 = t30 * t591;
t1129 = t31 * t590;
t875 = rSges(6,2) * t521 - rSges(6,3) * t520;
t276 = rSges(6,1) * t1063 - t875;
t1044 = t279 + t358;
t61 = t1044 * t591 + t276 * t590 + t751;
t1122 = t61 * t276;
t1046 = -t276 - t352;
t492 = qJD(5) * t522;
t1016 = -t580 - t635;
t654 = t691 * t995;
t913 = t1031 * t694;
t787 = qJD(2) * t913 + t654;
t726 = (t373 + t1016) * qJD(1) + t787;
t712 = -t1207 * t591 + t492 + t726;
t68 = t1046 * t667 + t591 * t744 + t712;
t1121 = t68 * t744;
t682 = t692 * rSges(3,3);
t79 = t1212 + t726;
t1120 = t694 * t79;
t1119 = t95 * t457;
t1118 = t96 * t458;
t1117 = t97 * t457;
t1116 = t98 * t458;
t1115 = t99 * t457;
t1114 = -rSges(4,3) - qJ(3);
t1082 = t100 * t458;
t1008 = rSges(3,2) * t1063 + t694 * rSges(3,3);
t554 = rSges(3,1) * t1061 - t1008;
t949 = t628 * t997;
t313 = -t949 + (-t554 - t635) * qJD(1);
t1077 = t313 * t692;
t1076 = t313 * t694;
t555 = rSges(3,1) * t1060 - rSges(3,2) * t1062 + t682;
t314 = qJD(1) * t555 - t1268;
t582 = t628 * t694;
t1075 = t314 * t582;
t1065 = t690 * t691;
t1064 = t690 * t693;
t316 = t693 * t352;
t1056 = t115 * t1062 + t352 * t952;
t114 = -t238 * pkin(4) - qJ(5) * t237 + t492;
t976 = rSges(6,1) * t952 + t238 * rSges(6,2) - t237 * rSges(6,3);
t139 = -rSges(6,1) * t958 + t976;
t1055 = t114 + t139;
t1054 = -t1254 * t958 + t1203;
t1053 = t1253 * t240 + t1254 * t570 - t1199;
t350 = -pkin(4) * t520 + qJ(5) * t521;
t1051 = t590 * t350 + t623;
t766 = t675 * t998 + t676 * t994;
t767 = -t675 * t994 + t676 * t998;
t241 = pkin(4) * t767 + qJ(5) * t766 + t622;
t482 = -t693 * t874 + t1145;
t549 = (rSges(6,2) * t675 + rSges(6,3) * t676) * t691;
t299 = qJD(2) * t482 + qJD(4) * t549;
t1050 = -t241 - t299;
t356 = -pkin(4) * t522 + qJ(5) * t523;
t1049 = qJD(5) * t521 + t667 * t356;
t1000 = qJD(2) * t691;
t1043 = t358 * t1000 + t691 * t478;
t1042 = t1063 * t1207 + t316;
t372 = t567 * rSges(4,1) + t566 * rSges(4,2) + rSges(4,3) * t1062;
t1041 = -t372 - t583;
t1040 = -t374 - t583;
t647 = pkin(2) * t955;
t1010 = t652 - t647;
t376 = qJ(3) * t570 + qJD(1) * t669 + t1010;
t1039 = -t376 - t1234;
t481 = t693 * t873 + t1143;
t552 = (rSges(7,2) * t676 - rSges(7,3) * t675) * t691;
t1038 = -pkin(5) * t1000 - qJ(6) * t767 - qJD(2) * t481 - qJD(4) * t552 - t621;
t1037 = t1248 * t692;
t1036 = t1248 * t694;
t451 = -t939 - t1202;
t560 = qJD(2) * t630 - t996;
t1035 = -qJD(2) * t451 - t560;
t1030 = -t451 - t630;
t519 = t693 * t880 + t1139;
t1029 = -t519 * qJD(2) - t560;
t1028 = -qJ(6) * t676 * t693 - t1150 - t481;
t1026 = t744 - t1207;
t548 = (-pkin(4) * t675 + qJ(5) * t676) * t691;
t1025 = qJD(5) * t523 - t591 * t548;
t1024 = -t518 - t627;
t1023 = -t519 - t630;
t1019 = t692 * t580 + t694 * t583;
t581 = t627 * t694;
t1018 = -qJD(1) * t581 + t692 * t996;
t983 = qJD(2) * qJD(3);
t935 = t693 * t983;
t984 = qJD(1) * qJD(2);
t937 = t692 * t984;
t1017 = t627 * t937 + t694 * t935;
t1014 = qJD(1) * t980 + t690 * t958;
t1011 = rSges(3,2) * t958 + rSges(3,3) * t1001;
t674 = pkin(7) * t1001;
t1009 = t654 + t674;
t1004 = qJD(1) * t612;
t987 = qJD(6) * t676;
t981 = -pkin(4) - t1253;
t978 = t114 + t1054;
t977 = -t238 * rSges(5,1) + t237 * rSges(5,2) + rSges(5,3) * t952;
t975 = -t241 + t1038;
t974 = -t241 + t1035;
t973 = -t352 - t1048;
t480 = t693 * t877 + t1137;
t550 = (-rSges(5,1) * t675 - rSges(5,2) * t676) * t691;
t297 = qJD(2) * t480 + qJD(4) * t550;
t971 = -t297 + t1035;
t640 = qJ(3) * t952;
t375 = -pkin(2) * t769 - qJ(3) * t958 + t640 + t654;
t970 = t580 * t1001 + t694 * t375 + t692 * t376;
t969 = -t358 + t1040;
t391 = t757 * t694;
t968 = qJD(1) * t391 + t1018;
t967 = t431 * rSges(4,1) + t430 * rSges(4,2) + rSges(4,3) * t952;
t965 = -t479 + t1031;
t964 = -t1207 + t1248;
t282 = t523 * rSges(5,1) - t522 * rSges(5,2) + rSges(5,3) * t1062;
t963 = -t578 * t999 - t581 * t997 + t677;
t962 = -t570 * t690 - t671 * t955;
t960 = t654 + t1204;
t959 = -pkin(7) - t1153;
t948 = t675 * t986;
t947 = t693 * t987;
t940 = -pkin(1) - t1146;
t936 = t694 * t984;
t931 = t1001 / 0.2e1;
t930 = t1000 / 0.2e1;
t929 = -t999 / 0.2e1;
t928 = t999 / 0.2e1;
t926 = t997 / 0.2e1;
t923 = -t671 * t693 - pkin(1);
t922 = -pkin(1) + t1065;
t921 = t44 * t1048;
t56 = t1248 * t591 + t667 * t973 + t489 + t712;
t920 = t56 * t1248;
t919 = t69 * t1026;
t154 = qJD(1) * t372 - t518 * t999 + t901;
t917 = t154 * t1024;
t916 = t692 * t1030;
t915 = t694 * t1030;
t914 = t1035 * t694;
t912 = t1024 * t694;
t910 = -t524 + t1071;
t908 = t241 * t1063 + t693 * t115 + t1207 * t570;
t907 = -t299 + t974;
t905 = t375 * t997 + t376 * t999 + t580 * t936 + t691 * t983;
t904 = -t692 * t373 + t694 * t374 + t1019;
t588 = qJD(1) * (-pkin(1) * t1003 + t674);
t903 = t692 * t935 + t588 + (t375 + t654) * qJD(1);
t902 = t1026 + t1031;
t900 = -t652 - t962;
t899 = t1009 + t1014;
t891 = t57 * t964;
t888 = qJD(4) * t930;
t887 = t1040 * t1003;
t876 = -rSges(6,2) * t240 + rSges(6,3) * t239;
t141 = rSges(6,1) * t570 + t876;
t171 = -t640 + (t938 - t1064) * t997 - t1237 * t1003 + t1014;
t172 = -qJ(3) * t953 + t647 + (t1237 * t694 + t656) * qJD(1) + t962;
t770 = t171 * t997 + t172 * t999 - t373 * t936 + t905;
t715 = qJD(2) * t887 + t770;
t8 = -t1044 * t457 + t1055 * t591 + t141 * t590 + t276 * t458 + t1176 + t715;
t886 = t61 * t141 + t8 * t276;
t883 = -rSges(3,2) * t691 + t1146;
t882 = rSges(4,1) * t433 + rSges(4,2) * t432;
t881 = -rSges(4,1) * t565 + rSges(4,2) * t564;
t879 = rSges(5,1) * t240 - rSges(5,2) * t239;
t855 = qJD(1) * t171 + t903;
t854 = t591 * t352;
t853 = t974 + t1038;
t456 = t694 * t1207;
t852 = t358 * t994 - t667 * t456 + t968;
t851 = t390 * t999 + t391 * t997 + t963;
t850 = t964 + t1031;
t827 = t280 * t694 - t282 * t692;
t826 = -t314 * t692 - t1076;
t815 = t685 + t892;
t813 = t923 - t1145;
t812 = t923 - t1137;
t811 = -t373 * t1001 + t694 * t171 + t692 * t172 + t970;
t810 = t692 * t352 + t694 * t358 + t904;
t808 = t900 - t1052;
t579 = t628 * t692;
t784 = -t1044 * t61 - t1121;
t783 = t1046 * t68 + t69 * t279;
t306 = (t554 * t692 + t555 * t694) * qJD(2);
t771 = t1114 * t691 - pkin(1) - t1154;
t765 = t358 + t1007 + t1015;
t7 = t978 * t591 + t770 - t972 * t457 - qJD(4) * t948 + t1053 * t590 + (t887 + t947) * qJD(2) + t1048 * t458 + t1176;
t758 = t1048 * t7 + t1053 * t44;
t756 = qJD(5) * t239 + t667 * t114 + t358 * t934 + t855;
t754 = t991 * t316 - t590 * t455 + t624 + t851;
t370 = rSges(4,3) * t1063 - t881;
t750 = -t44 * t972 - t920;
t749 = t1045 * t57 + t56 * t973;
t748 = t352 * t1001 + t694 * t114 + t692 * t115 + t811;
t740 = qJD(1) * t913 + t1035 * t692;
t725 = (-t987 - t989) * t691 + t1030 * qJD(2);
t724 = -t757 * t937 + (-t172 - t652 + t1039) * qJD(1) + t1017;
t714 = -t671 * t954 - t690 * t952 + t899;
t328 = qJD(1) * t373;
t713 = t1031 * t997 - t352 * t667 + t328 + t492 + t960;
t76 = t280 * t590 + t282 * t591 + t779;
t80 = t282 * t667 - t479 * t590 + t809;
t707 = t76 * t827 + (t692 * t79 - t694 * t80) * t479;
t706 = -qJD(5) * t237 + t1207 * t457 - t591 * t241 + t724;
t702 = (t919 + t1122) * t694 + t784 * t692;
t700 = (t891 + t921) * t694 + t750 * t692;
t595 = t883 * qJD(2);
t569 = (t692 ^ 2 + t694 ^ 2) * t1000;
t441 = t518 * t694;
t440 = t518 * t692;
t439 = t516 * t694;
t438 = t516 * t692;
t437 = t514 * t694;
t436 = t514 * t692;
t420 = t744 * t694;
t418 = t744 * t692;
t416 = t479 * t694;
t415 = t479 * t692;
t386 = -qJD(2) * t579 + (t694 * t883 + t682) * qJD(1);
t385 = -rSges(3,1) * t769 - rSges(3,2) * t952 + t1011;
t359 = rSges(7,2) * t523 - rSges(7,3) * t522;
t357 = -rSges(5,1) * t522 - rSges(5,2) * t523;
t355 = rSges(6,2) * t522 + rSges(6,3) * t523;
t353 = rSges(7,2) * t521 - rSges(7,3) * t520;
t351 = -rSges(5,1) * t520 - rSges(5,2) * t521;
t349 = rSges(6,2) * t520 + rSges(6,3) * t521;
t293 = t352 * t1062;
t182 = rSges(4,3) * t570 + t882;
t181 = -rSges(4,3) * t958 + t967;
t174 = -t595 * t997 + (-t386 + t1268) * qJD(1);
t173 = -t595 * t999 + t588 + (t385 - t949) * qJD(1);
t153 = t654 + qJD(2) * t912 + (-t370 + t1016) * qJD(1);
t143 = rSges(5,3) * t570 + t879;
t142 = -rSges(5,3) * t958 + t977;
t113 = (t370 * t692 + t372 * t694) * qJD(2) + t885;
t78 = t1029 * t997 + (-t182 + (qJD(2) * t518 - t677) * t692 + t1039) * qJD(1) + t1017;
t77 = qJD(1) * t181 + (qJD(1) * t912 + t1029 * t692) * qJD(2) + t903;
t58 = (t181 * t694 + t182 * t692 + (t1041 * t692 + t370 * t694) * qJD(1)) * qJD(2) + t905;
t33 = -t143 * t667 - t297 * t591 + t457 * t479 + (-t280 * t994 + t914) * qJD(2) + t724;
t32 = t142 * t667 - t297 * t590 - t458 * t479 + (t282 * t994 + t740) * qJD(2) + t855;
t13 = t142 * t591 + t143 * t590 + t280 * t458 - t282 * t457 + t715;
t12 = t139 * t667 + t1050 * t590 + t1026 * t458 + (t279 * t994 + t740) * qJD(2) + t756;
t11 = -t299 * t591 - t457 * t744 + (-t115 - t141) * t667 + (t1046 * t994 + t914) * qJD(2) + t706;
t10 = qJD(6) * t240 + t1054 * t667 + t975 * t590 + t964 * t458 + (t1045 * t994 + t740) * qJD(2) + t756;
t9 = -qJD(6) * t238 + t1038 * t591 - t1248 * t457 + (-t115 - t1053) * t667 + (t973 * t994 + t914) * qJD(2) + t706;
t1 = [(((t694 * t910 - t1032 + t186 + t825) * t694 + (t692 * t910 - t102 - t1209 + t1271 + t911) * t692) * qJD(2) + t1277 - t1288) * t929 + (-t917 * t997 + t77 * (t1007 - t1041) + (-t882 - t1010 + t771 * t1001 + (-pkin(7) * qJD(1) + t1114 * t998) * t692) * t153 + (t692 * t771 + t685 + t881) * t78 + (-pkin(2) * t954 + t1009 + t153 + t640 - t960 + t967 + (t370 + (-pkin(1) - t630 - t1139) * t692) * qJD(1)) * t154) * m(4) + (((t184 - t446 + (t525 + t1072) * t694 + t1033) * t694 + t1032 * t692 + t1213) * qJD(2) + t1294) * t926 + (t174 * (t692 * t940 + t1008 + t685) + t173 * (t555 + t1007) + t314 * (t674 + t1011) + (t1077 * t628 - t1075) * qJD(2) + ((-pkin(1) - t883) * t1076 + (t313 * (-rSges(3,3) - pkin(7)) + t314 * t940) * t692) * qJD(1) - (-qJD(1) * t554 - t313 - t604 - t949) * t314) * m(3) + (-t61 * t854 - (-t352 * t61 + t919) * t591 - (-t276 * t667 - t68 + t713) * t69 + t11 * (-t352 + t815 + t875) + t68 * (t808 - t876 - t1152) + t12 * (t279 + t765 - t979) + t69 * (t114 + t714 + t976) + (t11 * (t922 - t1145) - t68 * rSges(6,1) * t998) * t692 + (t68 * t813 * t694 + (t68 * t959 + t69 * t813) * t692) * qJD(1)) * m(6) + (-(t1204 + t1212 + t328 + t787 - t79) * t80 + t33 * (t815 - t878) + t79 * (-t879 + t900) + t32 * (t814 + t282 + t1007) + t80 * (t714 + t977) + (t33 * (t922 - t1137) - t79 * rSges(5,3) * t998) * t692 + (t812 * t1120 + (t79 * t959 + t80 * t812) * t692) * qJD(1)) * m(5) + (t1250 + t1251 + t1252 + (-t1305 * t675 - t1307 * t676) * t691 + t1287) * t590 * t1160 + (t1273 + t1274) * t928 + t1189 * t1165 + t1190 * t1166 + t1266 * t1163 + (t1267 + t1226) * t1162 + t1082 / 0.2e1 - (t1272 - t1275 + t1276) * t997 / 0.2e1 + t1115 / 0.2e1 + t1129 / 0.2e1 + t1180 + t1133 / 0.2e1 - t1134 / 0.2e1 - t1132 / 0.2e1 - t1130 / 0.2e1 + t1131 / 0.2e1 + ((t593 - t452) * t693 + (-t453 * t688 + t454 * t689 + t594) * t691 + (t512 * t691 + t693 * t819 - t816) * qJD(2)) * qJD(1) + ((-t1065 * t694 + t1045 + t765) * t10 + (t685 - t494 - t601 - t1140 + t981 * t521 + t1153 * t694 + (-pkin(1) + (t690 - t1254) * t691) * t692) * t9 - t44 * t854 - (-t352 * t44 + t891) * t591 + (t114 + t899 + (-t671 * t691 - t1064) * t997 + (-t1143 + t923 - t1150) * t1003 + t1203 - t489 - t713 + t1048 * t667) * t57 + (t808 + t981 * t240 + (-t1254 * t691 + t923) * t1001 + (qJD(1) * t959 - t1254 * t998) * t692 + t1199 + t57) * t56) * m(7) + t1226 * t1161 + t1118 / 0.2e1 + t1119 / 0.2e1 + t1116 / 0.2e1 + t1117 / 0.2e1 + ((t1269 + t1316) * t694 + (t1270 + t1317) * t692) * t984 / 0.2e1; (t1275 * t694 + t1274 * t692 + (t1269 * t694 + t1270 * t692) * qJD(1)) * qJD(1) / 0.2e1 + (-(t313 * t579 - t1075) * qJD(1) - (t306 * (-t579 * t692 - t582 * t694) + t826 * t883) * qJD(2) + 0.2e1 * t306 * (t385 * t694 + t386 * t692 + (t554 * t694 - t555 * t692) * qJD(1)) + t826 * t595 + (-t173 * t692 - t174 * t694 + (-t314 * t694 + t1077) * qJD(1)) * t628) * m(3) - (t1272 * qJD(1) + t1235 + ((-t1063 * t176 + t178 * t564 - t180 * t565 - t361 * t570 - t364 * t432 - t367 * t433 + t1289) * t694 + (t1063 * t175 - t1179 * t694 - t177 * t564 + t179 * t565 + t363 * t570 + t366 * t432 + t369 * t433) * t692 + (t1318 * t694 + t1319 * t692) * qJD(1)) * t1255) * t694 / 0.2e1 + ((t1061 * t1262 + t1189 * t691) * qJD(4) + ((qJD(4) * t1261 + t1217) * t693 + t1218) * t694 + (t1219 * t523 + t1220 * t522) * t667 + (t1222 * t523 - t1224 * t522) * t591 + (-t1221 * t523 + t1223 * t522) * t590) * t1164 + ((t1060 * t1263 + t1190 * t691) * qJD(4) + ((qJD(4) * t1264 + t1217) * t693 + t1218) * t692 + (t1219 * t521 + t1220 * t520) * t667 + (t1222 * t521 - t1224 * t520) * t591 + (-t1221 * t521 + t1223 * t520) * t590) * t1161 + (-t69 * (t420 * t667 - t622 * t692 + t852 + (-t482 - t553) * t590) - t61 * (t418 * t590 + t754 + (t420 - t456) * t591) - t69 * t916 * qJD(2) - (t691 * t783 + t693 * t702) * qJD(4) + t8 * t810 + t61 * t748 + (t11 * t902 + t8 * t279 + t61 * t139 + (t69 * t902 + t1122) * qJD(1)) * t694 + (t12 * t902 + t69 * t907 + (-t1121 + t61 * (-t279 + t969)) * qJD(1) + t886) * t692 + (t482 * t591 - (-t418 + t455) * t667 - t915 * qJD(2) + (t622 + t907) * t694 + t1244) * t68) * m(6) + (t7 * t810 + t44 * t748 + (t9 * t850 + t7 * t1045 + t44 * t1054 + (t57 * t850 + t921) * qJD(1)) * t694 + (t10 * t850 + t57 * t853 + (-t920 + t44 * (t969 - t1045)) * qJD(1) + t758) * t692 - (t691 * t749 + t693 * t700) * qJD(4) - (t754 + t947 + (-t456 + t1036) * t591 + t1037 * t590) * t44 - (t852 + t725 * t692 + t1036 * t667 + (-t553 + t1028) * t590) * t57 + (-(t455 - t1037) * t667 - t1028 * t591 + (-t725 + t853) * t694 + t1244) * t56) * m(7) + (t1226 + t1276 + t1291) * t931 + (t1227 + t1277 + t1290) * t1003 / 0.2e1 + ((qJD(4) * t1183 - t1280) * t693 + ((t1219 * t676 + t1220 * t675 + t1283) * t667 + (t1222 * t676 - t1224 * t675 - t1303) * t591 + (-t1221 * t676 + t1223 * t675 + t1284) * t590 + t1188 * qJD(4)) * t691) * t1160 + (t79 * t1034 + t13 * t904 + (t13 * t282 + t79 * t971 + (qJD(1) * t80 + t33) * t965) * t694 + (qJD(1) * t479 * t79 + t13 * t280 + t32 * t965 + t80 * t971) * t692 - t79 * (t415 * t667 - t480 * t591 + t782) - t80 * (-t416 * t667 - t480 * t590 + t968) - (t79 * t915 + t80 * t916) * qJD(2) - ((-t280 * t79 + t282 * t80) * t691 + t707 * t693) * qJD(4) + (t811 + (qJD(1) * t280 + t142) * t694 + (t143 + (-t282 + t1040) * qJD(1)) * t692 + t415 * t590 + t416 * t591 - t851) * t76) * m(5) + t1215 * t1165 + t1216 * t888 + (t1273 * qJD(1) + t1236 + ((-t1062 * t176 - t178 * t566 - t180 * t567 - t361 * t571 - t364 * t430 - t367 * t431) * t694 + (t1062 * t175 + t1179 * t692 + t177 * t566 + t179 * t567 + t363 * t571 + t366 * t430 + t369 * t431 - t1289) * t692 + ((t104 + t186) * t694 + t1271 * t692) * qJD(1)) * t1255) * t692 / 0.2e1 + t1214 * t1166 + (qJD(1) * t1181 + t1230 * t692 - t1231 * t694) * t1162 + (qJD(1) * t1182 + t1232 * t692 - t1233 * t694) * t1163 - t1225 * t994 / 0.2e1 + (qJD(1) * t1183 + t1228 * t692 - t1229 * t694) * t1159 + ((-t1069 * t999 + t1004) * t692 + (-t437 * t566 - t439 * t567) * t999 + (t515 * t566 + t517 * t567) * qJD(1) + ((t1070 * t692 + t436 * t566 + t438 * t567 + t1246) * qJD(2) + t1247) * t694) * t929 + ((-t1070 * t997 - t1004) * t694 - (t436 * t564 - t438 * t565) * t997 + (-t515 * t564 + t517 * t565) * qJD(1) + ((t1069 * t694 + t437 * t564 - t439 * t565 + t1246) * qJD(2) + t1247) * t692) * t926 - (-t1210 * t693 + ((-t515 * t688 + t517 * t689 + t512) * qJD(1) + ((t437 * t688 - t439 * t689 + t363) * t692 - (t436 * t688 - t438 * t689 + t361) * t694) * qJD(2)) * t691 + (t1012 * t693 + t1013 * t691) * qJD(1) + (t1174 * t693 + t753 * t691) * qJD(2)) * qJD(1) / 0.2e1 + (-t153 * (qJD(1) * t440 + t909) - t154 * (-qJD(1) * t441 + t1018) - t113 * t963 - ((t1023 * t153 - t113 * t441) * t694 + (t1023 * t154 - t113 * t440) * t692) * qJD(2) + t153 * t587 + t58 * t1019 + t113 * t970 + (t78 * t1024 + t153 * t1029 + t58 * t372 + t113 * t181 + (t113 * t370 + t917) * qJD(1)) * t694 + (t77 * t1024 + t154 * t1029 + t58 * t370 + t113 * t182 + (t1041 * t113 + t153 * t518) * qJD(1)) * t692) * m(4) - (t1226 * t694 + t1227 * t692) * t992 / 0.2e1; -m(4) * (t113 * t569 + t153 * t571 + t154 * t570) - m(5) * (t569 * t76 + t570 * t80 + t571 * t79) - m(6) * (t569 * t61 + t570 * t69 + t571 * t68) - m(7) * (t44 * t569 + t56 * t571 + t57 * t570) + 0.2e1 * ((t153 * t997 + t154 * t999 - t58) * t1170 + (t79 * t997 + t80 * t999 - t13) * t1169 + (t68 * t997 + t69 * t999 - t8) * t1168 + (t56 * t997 + t57 * t999 - t7) * t1167) * t693 + 0.2e1 * ((qJD(2) * t113 + t1001 * t154 - t1003 * t153 + t692 * t77 + t694 * t78) * t1170 + (qJD(2) * t76 + t1001 * t80 - t1003 * t79 + t32 * t692 + t33 * t694) * t1169 + (qJD(2) * t61 + t1001 * t69 - t1003 * t68 + t11 * t694 + t12 * t692) * t1168 + (qJD(2) * t44 + t10 * t692 + t1001 * t57 - t1003 * t56 + t694 * t9) * t1167) * t691; (-t56 * (-qJD(6) * t522 - t552 * t591 + t1025 + (-t350 - t353) * t667) - t57 * (-qJD(6) * t520 + t359 * t667 + t1049 + (-t548 - t552) * t590) - t44 * (t353 * t590 + t1051 - t948 + (t356 + t359) * t591) - (t56 * t801 + t57 * t800 - t1135) * qJ(6) + t9 * t1042 + t56 * t908 + t57 * t1043 + t7 * t293 + t44 * t1056 + (qJD(2) * t700 - t10 * t972 + t1048 * t9 + t1053 * t56 - t57 * t978) * t693 + (t749 * qJD(2) + (qJD(1) * t750 + t10 * t964 + t57 * t975 + t758) * t694 + (-t9 * t1248 - t56 * t1038 - t7 * t972 - t44 * t978 + (-t1248 * t57 + t44 * t973) * qJD(1)) * t692) * t691) * m(7) + (-t68 * (-t549 * t591 + (-t349 - t350) * t667 + t1025) - t69 * (t355 * t667 + (-t548 - t549) * t590 + t1049) - t61 * (t349 * t590 + (t355 + t356) * t591 + t1051) + t11 * t1042 + t68 * t908 + t69 * t1043 + t8 * t293 + t61 * t1056 + (qJD(2) * t702 - t1044 * t12 - t1055 * t69 + t11 * t276 + t68 * t141) * t693 + (t783 * qJD(2) + (qJD(1) * t784 + t1026 * t12 + t1050 * t69 + t886) * t694 + (-t11 * t744 + t68 * t299 - t8 * t1044 - t61 * t1055 + (t1046 * t61 - t69 * t744) * qJD(1)) * t692) * t691) * m(6) + ((qJD(2) * t707 - t80 * t142 + t79 * t143 + t33 * t280 - t32 * t282) * t693 + (t79 * (-qJD(2) * t280 + t297 * t692) + t80 * (qJD(2) * t282 - t297 * t694) + t13 * t827 + t76 * (-t1001 * t282 - t1003 * t280 - t142 * t692 + t143 * t694) + (-t32 * t694 + t33 * t692 + (t692 * t80 + t1120) * qJD(1)) * t479) * t691 - t79 * (-t351 * t667 - t550 * t591) - t80 * (t357 * t667 - t550 * t590) - t76 * (t351 * t590 + t357 * t591)) * m(5) + t1235 * t1063 / 0.2e1 + t1236 * t1062 / 0.2e1 + t1225 * t930 + (t1183 * t691 - t1188 * t693) * t888 + (t1181 * t691 - t1190 * t693) * t1166 + (t1182 * t691 - t1189 * t693) * t1165 + (-t1184 * t694 + t1185 * t523 + t1186 * t522) * t1164 + ((qJD(2) * t1182 - t1266) * t693 + (-qJD(1) * t1215 + t1189 * qJD(2) + t1232 * t694 + t1233 * t692) * t691) * t1163 + ((qJD(2) * t1181 - t1267) * t693 + (-qJD(1) * t1214 + t1190 * qJD(2) + t1230 * t694 + t1231 * t692) * t691) * t1162 + (-t1184 * t692 + t1185 * t521 + t1186 * t520) * t1161 + (t1256 * t693 + (t1185 * t676 + t1186 * t675) * t691) * t1160 + ((qJD(2) * t1183 - t1265) * t693 + (-qJD(1) * t1216 + t1188 * qJD(2) + t1228 * t694 + t1229 * t692) * t691) * t1159 - (t1118 + t1119 + t1133 - t1134 + t1116 + t1117 + t1131 - t1132 + t1082 + t1115 + t1129 - t1130 + t1180) * t693 / 0.2e1 + (t691 * t931 + t693 * t928) * t1227 + (-t958 / 0.2e1 + t693 * t926) * t1226; (t10 * t520 + t1068 * t7 + t1205 * t57 + t1206 * t56 + t44 * t766 + t522 * t9 - t1135) * m(7) + (t1068 * t8 + t11 * t522 + t12 * t520 + t1205 * t69 + t1206 * t68 + (-t206 + t766) * t61) * m(6); (t10 * t521 + t1067 * t7 + t523 * t9 + (t1067 * t590 - t523 * t667 + t240) * t57 + (t1067 * t591 + t521 * t667 - t238) * t56 + (-t521 * t590 - t523 * t591 + t767) * t44) * m(7);];
tauc  = t1(:);
