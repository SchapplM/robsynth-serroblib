% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:33
% EndTime: 2019-12-31 21:08:59
% DurationCPUTime: 77.52s
% Computational Cost: add. (58709->1328), mult. (147294->1784), div. (0->0), fcn. (163265->6), ass. (0->730)
t804 = sin(qJ(2));
t805 = sin(qJ(1));
t1072 = t804 * t805;
t806 = cos(qJ(3));
t808 = cos(qJ(1));
t1065 = t806 * t808;
t807 = cos(qJ(2));
t1067 = t805 * t807;
t803 = sin(qJ(3));
t734 = t1067 * t803 + t1065;
t713 = Icges(6,6) * t734;
t1062 = t808 * t803;
t1066 = t806 * t807;
t735 = t1066 * t805 - t1062;
t527 = Icges(6,5) * t1072 + Icges(6,3) * t735 + t713;
t717 = Icges(5,6) * t735;
t530 = Icges(5,5) * t1072 + Icges(5,3) * t734 - t717;
t714 = Icges(6,6) * t735;
t534 = -Icges(6,4) * t1072 - Icges(6,2) * t734 - t714;
t716 = Icges(5,6) * t734;
t537 = -Icges(5,4) * t1072 + Icges(5,2) * t735 - t716;
t539 = Icges(6,1) * t1072 + Icges(6,4) * t734 + Icges(6,5) * t735;
t542 = Icges(5,1) * t1072 - Icges(5,4) * t735 + Icges(5,5) * t734;
t545 = Icges(4,5) * t735 - Icges(4,6) * t734 + Icges(4,3) * t1072;
t720 = Icges(4,4) * t735;
t548 = -Icges(4,2) * t734 + Icges(4,6) * t1072 + t720;
t719 = Icges(4,4) * t734;
t552 = -Icges(4,1) * t735 - Icges(4,5) * t1072 + t719;
t1325 = (t527 + t537 - t552) * t735 + (t530 - t534 - t548) * t734 + (t539 + t542 + t545) * t1072;
t1070 = t804 * t808;
t1068 = t805 * t806;
t736 = t1062 * t807 - t1068;
t1098 = Icges(6,6) * t736;
t1064 = t807 * t808;
t1069 = t805 * t803;
t737 = t1064 * t806 + t1069;
t529 = Icges(6,5) * t1070 + Icges(6,3) * t737 + t1098;
t1100 = Icges(5,6) * t737;
t532 = Icges(5,5) * t1070 + Icges(5,3) * t736 - t1100;
t715 = Icges(6,6) * t737;
t535 = Icges(6,4) * t1070 + Icges(6,2) * t736 + t715;
t718 = Icges(5,6) * t736;
t538 = Icges(5,4) * t1070 - Icges(5,2) * t737 + t718;
t541 = Icges(6,1) * t1070 + Icges(6,4) * t736 + Icges(6,5) * t737;
t544 = Icges(5,1) * t1070 - Icges(5,4) * t737 + Icges(5,5) * t736;
t547 = Icges(4,5) * t737 - Icges(4,6) * t736 + Icges(4,3) * t1070;
t1110 = Icges(4,4) * t737;
t550 = -Icges(4,2) * t736 + Icges(4,6) * t1070 + t1110;
t721 = Icges(4,4) * t736;
t553 = Icges(4,1) * t737 + Icges(4,5) * t1070 - t721;
t1324 = (t529 - t538 + t553) * t735 + (t532 + t535 - t550) * t734 + (t541 + t544 + t547) * t1072;
t1097 = Icges(6,6) * t803;
t896 = Icges(6,3) * t806 + t1097;
t682 = -Icges(6,5) * t807 + t804 * t896;
t1073 = t803 * t804;
t1106 = Icges(6,4) * t807;
t1071 = t804 * t806;
t779 = Icges(6,6) * t1071;
t686 = Icges(6,2) * t1073 - t1106 + t779;
t902 = Icges(6,4) * t803 + Icges(6,5) * t806;
t690 = -Icges(6,1) * t807 + t804 * t902;
t423 = t1072 * t690 + t682 * t735 + t686 * t734;
t1107 = Icges(5,4) * t807;
t780 = Icges(5,6) * t1073;
t688 = -Icges(5,2) * t1071 - t1107 + t780;
t1099 = Icges(5,6) * t806;
t897 = -Icges(5,3) * t803 + t1099;
t829 = Icges(5,5) * t807 + t804 * t897;
t903 = Icges(5,4) * t806 - Icges(5,5) * t803;
t834 = Icges(5,1) * t807 + t804 * t903;
t425 = -t1072 * t834 - t688 * t735 - t734 * t829;
t900 = Icges(4,5) * t806 - Icges(4,6) * t803;
t669 = -Icges(4,3) * t807 + t804 * t900;
t1108 = Icges(4,4) * t806;
t904 = -Icges(4,2) * t803 + t1108;
t673 = -Icges(4,6) * t807 + t804 * t904;
t1109 = Icges(4,4) * t803;
t905 = Icges(4,1) * t806 - t1109;
t677 = -Icges(4,5) * t807 + t804 * t905;
t429 = t1072 * t669 - t673 * t734 + t677 * t735;
t1262 = t423 + t425 + t429;
t1323 = t1324 * t808 + t1325 * t805;
t1191 = t804 / 0.2e1;
t1190 = t805 / 0.2e1;
t1188 = -t807 / 0.2e1;
t1255 = -t808 / 0.2e1;
t1214 = m(6) / 0.2e1;
t1216 = m(5) / 0.2e1;
t777 = pkin(2) * t804 - pkin(7) * t807;
t1135 = pkin(3) * t806;
t907 = qJ(4) * t803 + t1135;
t756 = t907 * t804;
t1297 = rSges(6,1) + pkin(4);
t1121 = rSges(6,2) * t803;
t908 = rSges(6,3) * t806 + t1121;
t852 = t908 * t804;
t988 = qJ(5) * t1071 - t1297 * t807 + t852;
t954 = t756 + t988;
t933 = t777 + t954;
t468 = t933 * t805;
t470 = t933 * t808;
t909 = rSges(5,2) * t806 - rSges(5,3) * t803;
t853 = t909 * t804;
t699 = -rSges(5,1) * t807 - t853;
t987 = t699 + t756;
t953 = t777 + t987;
t507 = t953 * t805;
t509 = t953 * t808;
t1181 = rSges(5,2) - pkin(3);
t709 = t734 * qJ(4);
t723 = t734 * rSges(5,3);
t1182 = rSges(5,1) + pkin(7);
t1136 = pkin(2) * t807;
t951 = pkin(1) + t1136;
t1236 = t1182 * t804 + t951;
t799 = t808 * pkin(6);
t836 = -t1236 * t805 + t799;
t420 = t1181 * t735 - t709 - t723 + t836;
t1117 = t736 * rSges(5,3);
t1134 = t805 * pkin(6);
t711 = t736 * qJ(4);
t944 = t711 + t1134;
t421 = -t1181 * t737 + t1236 * t808 + t1117 + t944;
t845 = (-t420 * t808 - t421 * t805) * t1073;
t724 = t734 * rSges(6,2);
t969 = pkin(7) + t1297;
t1238 = t804 * t969 + t951;
t825 = -t1238 * t805 + t799;
t1116 = rSges(6,3) + qJ(5);
t967 = pkin(3) + t1116;
t401 = -t735 * t967 - t709 - t724 + t825;
t1239 = t736 * rSges(6,2) + t737 * qJ(5);
t402 = (rSges(6,3) + pkin(3)) * t737 + t1238 * t808 + t944 + t1239;
t851 = t804 * (-t401 * t808 - t402 * t805);
t847 = t803 * t851;
t1057 = (-t736 * t468 + t734 * t470 + t847) * t1214 + (-t736 * t507 + t734 * t509 + t845) * t1216;
t667 = (-pkin(3) * t1068 - qJ(4) * t1069) * t804;
t789 = pkin(2) * t1072;
t952 = t789 - t667;
t452 = (-t969 * t807 + (t1116 * t806 + t1121) * t804) * t805 + t952;
t1241 = t1297 * t1064;
t790 = pkin(7) * t1064;
t453 = t790 + (-pkin(2) + (-rSges(6,2) - qJ(4)) * t803 - t967 * t806) * t1070 + t1241;
t961 = t804 * t1065;
t977 = rSges(5,1) * t1064 + rSges(5,2) * t961;
t483 = t790 + (-t1135 - pkin(2) + (-rSges(5,3) - qJ(4)) * t803) * t1070 + t977;
t484 = (-t1182 * t807 - t853) * t805 + t952;
t1058 = (t736 * t452 + t734 * t453 + t847) * t1214 + (t734 * t483 + t736 * t484 + t845) * t1216;
t19 = t1058 - t1057;
t1322 = t19 * qJD(1);
t1014 = t737 * rSges(6,3) + t1297 * t1070 + t1239;
t1305 = t1116 * t735 + t724;
t1274 = t1297 * t1072 + t1305;
t608 = -t735 * pkin(3) - t709;
t612 = t737 * pkin(3) + t711;
t1225 = t808 ^ 2;
t1226 = t805 ^ 2;
t1302 = t1225 + t1226;
t778 = pkin(7) * t804 + t1136;
t985 = t1302 * t778;
t937 = -t608 * t805 + t808 * t612 + t985;
t265 = t1014 * t808 + t1274 * t805 + t937;
t981 = t735 * rSges(5,2) - t723;
t556 = rSges(5,1) * t1072 - t981;
t559 = rSges(5,1) * t1070 - t737 * rSges(5,2) + t1117;
t298 = t805 * t556 + t559 * t808 + t937;
t710 = t735 * qJ(4);
t604 = -pkin(3) * t734 + t710;
t610 = -t736 * pkin(3) + qJ(4) * t737;
t1010 = t805 * t604 + t808 * t610;
t577 = t805 * t734 + t736 * t808;
t1122 = rSges(6,2) * t735;
t607 = -rSges(6,3) * t734 + t1122;
t613 = t737 * rSges(6,2) - t736 * rSges(6,3);
t301 = -qJ(5) * t577 + t805 * t607 + t613 * t808 + t1010;
t1119 = rSges(5,3) * t735;
t603 = rSges(5,2) * t734 + t1119;
t609 = t736 * rSges(5,2) + t737 * rSges(5,3);
t355 = t805 * t603 + t609 * t808 + t1010;
t1242 = -t735 * rSges(4,1) + t734 * rSges(4,2);
t561 = -rSges(4,3) * t1072 + t1242;
t911 = t737 * rSges(4,1) - t736 * rSges(4,2);
t562 = rSges(4,3) * t1070 + t911;
t396 = -t561 * t805 + t562 * t808 + t985;
t605 = -rSges(4,1) * t734 - rSges(4,2) * t735;
t611 = -rSges(4,1) * t736 - rSges(4,2) * t737;
t451 = t805 * t605 + t611 * t808;
t753 = (-pkin(3) * t803 + qJ(4) * t806) * t804;
t757 = (rSges(6,2) * t806 - rSges(6,3) * t803) * t804;
t979 = -t753 - t757;
t871 = qJ(5) * t1073 + t979;
t520 = t871 * t805;
t521 = t871 * t808;
t754 = (rSges(5,2) * t803 + rSges(5,3) * t806) * t804;
t980 = -t753 - t754;
t579 = t980 * t805;
t580 = t980 * t808;
t1125 = rSges(4,1) * t806;
t910 = -rSges(4,2) * t803 + t1125;
t854 = t910 * t804;
t693 = -rSges(4,3) * t807 + t854;
t992 = t693 + t777;
t581 = t992 * t805;
t583 = t992 * t808;
t755 = (-rSges(4,1) * t803 - rSges(4,2) * t806) * t804;
t1321 = -m(6) * (t265 * t301 - t468 * t520 - t470 * t521) - m(5) * (t298 * t355 - t507 * t579 - t509 * t580) - m(4) * (t396 * t451 + (t581 * t805 + t583 * t808) * t755);
t926 = t1262 * t1188 + t1323 * t1191;
t664 = t756 * t1072;
t1001 = t699 * t1072 + t664;
t1084 = t556 * t807;
t572 = t807 * t608;
t391 = t1001 - t572 + t1084;
t1315 = t556 - t608;
t392 = t1315 * t807 + t1001;
t1037 = t391 - t392;
t943 = t1274 * t807;
t958 = t988 * t1072 + t664;
t321 = -t572 + t943 + t958;
t1304 = -t608 + t1274;
t322 = t1304 * t807 + t958;
t1042 = t321 - t322;
t1299 = -m(5) / 0.2e1;
t1298 = -m(6) / 0.2e1;
t1319 = t1298 * t468;
t1320 = -t1037 * t507 * t1299 - t1042 * t1319;
t1086 = t545 * t807;
t1301 = t548 * t803 + t552 * t806;
t372 = t1301 * t804 + t1086;
t427 = t1070 * t690 + t737 * t682 + t736 * t686;
t328 = t539 * t1070 + t737 * t527 - t736 * t534;
t329 = t541 * t1070 + t737 * t529 + t736 * t535;
t892 = t805 * t328 + t329 * t808;
t1310 = -t427 * t807 + t804 * t892;
t428 = -t1070 * t834 - t737 * t688 - t736 * t829;
t330 = t542 * t1070 + t736 * t530 + t737 * t537;
t331 = t544 * t1070 + t736 * t532 - t737 * t538;
t891 = t805 * t330 + t331 * t808;
t1311 = -t428 * t807 + t804 * t891;
t431 = t1070 * t669 - t736 * t673 + t737 * t677;
t334 = t545 * t1070 - t736 * t548 - t737 * t552;
t335 = t547 * t1070 - t736 * t550 + t737 * t553;
t889 = t805 * t334 + t335 * t808;
t1312 = -t431 * t807 + t804 * t889;
t925 = t1310 / 0.2e1 + t1311 / 0.2e1 + t1312 / 0.2e1;
t1083 = t577 * t806;
t1254 = m(6) * t804;
t578 = t805 * t735 + t737 * t808;
t1303 = (t578 * t803 - t1083) * t1254;
t433 = -t1303 / 0.2e1;
t432 = t1303 / 0.2e1;
t1306 = m(6) * (-t734 * t737 + t735 * t736);
t474 = -t1306 / 0.2e1;
t475 = t1306 / 0.2e1;
t1090 = t539 * t807;
t882 = t527 * t806 - t534 * t803;
t374 = t804 * t882 - t1090;
t1088 = t542 * t807;
t880 = t530 * t803 + t537 * t806;
t377 = t804 * t880 - t1088;
t1314 = -t372 + t374 + t377;
t1085 = t547 * t807;
t877 = -t550 * t803 + t553 * t806;
t373 = t804 * t877 - t1085;
t1089 = t541 * t807;
t881 = t529 * t806 + t535 * t803;
t376 = t804 * t881 - t1089;
t1087 = t544 * t807;
t879 = t532 * t803 - t538 * t806;
t379 = t804 * t879 - t1087;
t1313 = t373 + t376 + t379;
t463 = t693 * t1072 - t561 * t807;
t1309 = -t334 * t808 + t335 * t805;
t1308 = -t330 * t808 + t331 * t805;
t1307 = -t328 * t808 + t329 * t805;
t573 = t808 * t608;
t1300 = -m(4) / 0.2e1;
t1284 = t693 * t805;
t1283 = t699 * t805;
t1282 = t805 * t988;
t671 = Icges(3,5) * t1067 - Icges(3,6) * t1072 - Icges(3,3) * t808;
t783 = Icges(3,4) * t1072;
t679 = Icges(3,1) * t1067 - Icges(3,5) * t808 - t783;
t1003 = -t679 * t1064 - t805 * t671;
t675 = Icges(3,4) * t1067 - Icges(3,2) * t1072 - Icges(3,6) * t808;
t1102 = Icges(3,2) * t804;
t797 = Icges(3,4) * t807;
t676 = Icges(3,6) * t805 + (t797 - t1102) * t808;
t1111 = Icges(3,4) * t804;
t768 = Icges(3,1) * t807 - t1111;
t680 = Icges(3,5) * t805 + t768 * t808;
t649 = t680 * t1067;
t764 = Icges(3,5) * t807 - Icges(3,6) * t804;
t672 = Icges(3,3) * t805 + t764 * t808;
t940 = t672 * t808 - t649;
t1281 = -t1070 * t675 - t1072 * t676 - t1003 - t940;
t1180 = rSges(4,3) + pkin(7);
t1237 = t1180 * t804 + t951;
t477 = -t1237 * t805 + t1242 + t799;
t1013 = t559 + t612;
t1240 = t1013 * t807;
t393 = t1070 * t987 + t1240;
t369 = t736 * t393;
t1041 = -t391 * t734 - t369;
t960 = t612 + t1014;
t1246 = t807 * t960;
t323 = t1070 * t954 + t1246;
t304 = t736 * t323;
t1044 = -t321 * t734 - t304;
t801 = t804 ^ 2;
t618 = t1069 * t801 + t734 * t807;
t620 = -t1062 * t801 - t807 * t736;
t1061 = (t401 * t618 + t402 * t620 + t1044) * t1214 + (t420 * t618 + t421 * t620 + t1041) * t1216;
t1131 = (t322 * t734 + t1044 + t304) * t1214 + (t392 * t734 + t1041 + t369) * t1216;
t1272 = t1061 - t1131;
t566 = t804 * t573;
t280 = -t566 + (t1274 * t808 - t805 * t960) * t804;
t320 = -t566 + (-t1013 * t805 + t556 * t808) * t804;
t525 = (t734 * t808 - t736 * t805) * t804;
t895 = -t321 * t808 + t323 * t805;
t850 = t895 * t804;
t885 = -t391 * t808 + t393 * t805;
t1132 = (t525 * t265 + t280 * t577 - t620 * t468 - t618 * t470 + t803 * t850) * t1214 + (t1073 * t885 + t525 * t298 + t320 * t577 - t620 * t507 - t618 * t509) * t1216;
t759 = t907 * t807;
t957 = t756 * t1067 + t759 * t1072 + t807 * t667;
t990 = qJ(5) * t1066 + t1297 * t804 + t807 * t908;
t218 = (t805 * t990 - t1304) * t804 + t957;
t568 = t804 * t612;
t956 = -t759 - t990;
t1005 = -qJ(5) * t961 - t808 * t852 + t1241;
t668 = t907 * t1070;
t959 = t668 - t1005;
t219 = t568 + (t808 * t956 + t1014) * t804 + (-t808 * t954 + t959) * t807;
t962 = t804 * t1062;
t647 = -rSges(5,3) * t962 + t977;
t1004 = -t647 + t668;
t1011 = t667 * t1070 - t807 * t573;
t232 = (-t1283 * t804 + t1084) * t808 + (t1004 * t804 - t1240) * t805 + t1011;
t697 = rSges(5,1) * t804 - t807 * t909;
t281 = (t697 * t805 - t1315) * t804 + t957;
t989 = -t697 - t759;
t282 = t568 + (t808 * t989 + t559) * t804 + (-t808 * t987 + t1004) * t807;
t176 = (-t1282 * t804 + t943) * t808 + (t804 * t959 - t1246) * t805 + t1011;
t824 = t280 * t807 + (t176 + t895) * t804;
t1133 = (t736 * t218 + t734 * t219 + t803 * t824) * t1214 + (t736 * t281 + t734 * t282 + (t320 * t807 + (t232 + t885) * t804) * t803) * t1216;
t1270 = t1132 - t1133;
t626 = t673 * t805;
t628 = t677 * t805;
t859 = -t669 * t805 + t1301;
t273 = -t859 * t807 + (t626 * t803 - t628 * t806 + t545) * t804;
t630 = t682 * t805;
t898 = Icges(6,2) * t803 + Icges(6,6) * t806;
t831 = -t804 * t898 + t1106;
t634 = t831 * t805;
t863 = t690 * t805 + t882;
t275 = t863 * t807 + (-t630 * t806 + t634 * t803 + t539) * t804;
t632 = t829 * t805;
t899 = Icges(5,2) * t806 - Icges(5,6) * t803;
t832 = t804 * t899 + t1107;
t636 = t832 * t805;
t861 = -t834 * t805 + t880;
t277 = t861 * t807 + (t632 * t803 - t636 * t806 + t542) * t804;
t1269 = t273 + t275 + t277;
t627 = t673 * t808;
t629 = t677 * t808;
t858 = -t669 * t808 - t877;
t274 = -t858 * t807 + (t627 * t803 - t629 * t806 + t547) * t804;
t631 = t682 * t808;
t635 = t831 * t808;
t862 = t690 * t808 + t881;
t276 = t862 * t807 + (-t631 * t806 + t635 * t803 + t541) * t804;
t633 = t829 * t808;
t637 = t832 * t808;
t860 = -t834 * t808 + t879;
t278 = t860 * t807 + (t633 * t803 - t637 * t806 + t544) * t804;
t1268 = t274 + t276 + t278;
t1019 = Icges(4,2) * t735 + t552 + t719;
t1021 = -Icges(4,1) * t734 - t548 - t720;
t593 = -Icges(4,5) * t734 - Icges(4,6) * t735;
t290 = -t593 * t807 + (t1019 * t803 + t1021 * t806) * t804;
t1025 = -Icges(6,3) * t734 - t534 + t714;
t1029 = Icges(6,2) * t735 - t527 - t713;
t595 = Icges(6,4) * t735 - Icges(6,5) * t734;
t292 = -t595 * t807 + (t1025 * t806 + t1029 * t803) * t804;
t1023 = Icges(5,3) * t735 - t537 + t716;
t1027 = -Icges(5,2) * t734 + t530 - t717;
t597 = Icges(5,4) * t734 + Icges(5,5) * t735;
t294 = -t597 * t807 + (t1023 * t803 + t1027 * t806) * t804;
t1267 = t290 + t292 + t294;
t1018 = Icges(4,2) * t737 - t553 + t721;
t1020 = -Icges(4,1) * t736 - t1110 - t550;
t594 = -Icges(4,5) * t736 - Icges(4,6) * t737;
t291 = -t594 * t807 + (t1018 * t803 + t1020 * t806) * t804;
t1024 = -Icges(6,3) * t736 + t535 + t715;
t1028 = Icges(6,2) * t737 - t1098 - t529;
t596 = Icges(6,4) * t737 - Icges(6,5) * t736;
t293 = -t596 * t807 + (t1024 * t806 + t1028 * t803) * t804;
t1022 = Icges(5,3) * t737 + t538 + t718;
t1026 = -Icges(5,2) * t736 - t1100 + t532;
t598 = Icges(5,4) * t736 + Icges(5,5) * t737;
t295 = -t598 * t807 + (t1022 * t803 + t1026 * t806) * t804;
t1266 = t291 + t293 + t295;
t674 = Icges(4,6) * t804 + t807 * t904;
t678 = Icges(4,5) * t804 + t807 * t905;
t1081 = t669 * t807;
t670 = Icges(4,3) * t804 + t807 * t900;
t875 = -t673 * t803 + t806 * t677;
t857 = t670 - t875;
t817 = t804 * t857 + t1081;
t313 = -t674 * t734 + t678 * t735 + t805 * t817;
t681 = Icges(6,5) * t804 + t807 * t896;
t685 = Icges(6,4) * t804 + t807 * t898;
t1078 = t690 * t807;
t689 = Icges(6,1) * t804 + t807 * t902;
t873 = t806 * t682 + t686 * t803;
t856 = -t689 + t873;
t816 = -t804 * t856 + t1078;
t315 = t681 * t735 + t685 * t734 + t805 * t816;
t683 = Icges(5,5) * t804 - t807 * t897;
t687 = Icges(5,4) * t804 - t807 * t899;
t1077 = t834 * t807;
t691 = Icges(5,1) * t804 - t807 * t903;
t872 = -t806 * t688 - t803 * t829;
t855 = -t691 + t872;
t815 = -t804 * t855 - t1077;
t316 = t683 * t734 - t687 * t735 + t805 * t815;
t1265 = -t313 - t315 - t316;
t314 = -t736 * t674 + t737 * t678 + t808 * t817;
t317 = t737 * t681 + t736 * t685 + t808 * t816;
t318 = t736 * t683 - t737 * t687 + t808 * t815;
t1264 = -t314 - t317 - t318;
t363 = -t857 * t807 + (-t674 * t803 + t678 * t806 + t669) * t804;
t364 = t856 * t807 + (t681 * t806 + t685 * t803 + t690) * t804;
t365 = t855 * t807 + (t683 * t803 - t687 * t806 - t834) * t804;
t1263 = t365 + t364 + t363;
t1261 = t427 + t428 + t431;
t454 = t804 * t875 - t1081;
t459 = t804 * t873 - t1078;
t460 = t804 * t872 + t1077;
t1260 = t454 + t459 + t460;
t1259 = t1313 * t808 + t1314 * t805;
t1258 = -t804 / 0.2e1;
t1257 = -t805 / 0.2e1;
t1256 = t807 / 0.2e1;
t1184 = t808 / 0.2e1;
t1074 = t801 * t803;
t461 = t1074 * t806 + t734 * t735 + t736 * t737;
t1213 = m(6) / 0.4e1;
t1215 = m(5) / 0.4e1;
t968 = t1215 + t1213;
t1245 = t968 * t461;
t800 = t803 ^ 2;
t1244 = t968 * (-t577 * t803 + t800 * t807) * t804;
t972 = qJD(2) * t808;
t442 = -qJ(5) * t736 + t610 + t613;
t1007 = -t609 - t610;
t441 = t734 * t967 - t1122 - t710;
t472 = -t1181 * t734 - t1119 - t710;
t1059 = (t401 * t737 + t402 * t735 + t441 * t736 + t442 * t734) * t1214 + (-t1007 * t734 + t420 * t737 + t421 * t735 + t472 * t736) * t1216;
t1115 = (-t507 * t735 - t509 * t737 + t579 * t734 + t580 * t736 + (t298 * t806 + t355 * t803) * t804) * t1216 + (-t468 * t735 - t470 * t737 + t520 * t734 + t521 * t736 + (t265 * t806 + t301 * t803) * t804) * t1214;
t1235 = t273 / 0.2e1 + t275 / 0.2e1 + t277 / 0.2e1;
t1234 = -t274 / 0.2e1 - t276 / 0.2e1 - t278 / 0.2e1;
t478 = t1237 * t808 + t1134 + t911;
t1233 = (t401 * t521 + t402 * t520 - t441 * t470 - t442 * t468) * t1298 + (t1007 * t507 + t420 * t580 + t421 * t579 - t472 * t509) * t1299 + (-t581 * t611 + t583 * t605 + (-t477 * t808 - t478 * t805) * t755) * t1300;
t695 = rSges(4,3) * t804 + t807 * t910;
t394 = (t695 * t805 + t561) * t804;
t978 = rSges(4,2) * t962 + rSges(4,3) * t1064;
t643 = -rSges(4,1) * t961 + t978;
t395 = (-t693 * t808 - t643) * t807 + (-t695 * t808 + t562) * t804;
t465 = t1070 * t693 + t807 * t562;
t563 = t789 + (-t1180 * t807 + t854) * t805;
t564 = t790 + (-pkin(2) - t1125) * t1070 + t978;
t1232 = (t394 * t477 + t395 * t478 + t463 * t563 - t465 * t564) * t1300 + (t218 * t401 + t219 * t402 + t321 * t452 - t323 * t453) * t1298 + (t281 * t420 + t282 * t421 + t391 * t484 - t393 * t483) * t1299;
t738 = -Icges(6,3) * t1073 + t779;
t739 = Icges(5,3) * t1071 + t780;
t740 = (Icges(6,2) * t806 - t1097) * t804;
t741 = (Icges(5,2) * t803 + t1099) * t804;
t747 = (-Icges(4,2) * t806 - t1109) * t804;
t750 = (-Icges(4,1) * t803 - t1108) * t804;
t914 = t688 / 0.2e1 - t682 / 0.2e1 - t677 / 0.2e1;
t916 = -t829 / 0.2e1 + t686 / 0.2e1 - t673 / 0.2e1;
t1230 = t803 * (t739 / 0.2e1 + t740 / 0.2e1 - t747 / 0.2e1 + t914) - t806 * (t741 / 0.2e1 - t738 / 0.2e1 - t750 / 0.2e1 - t916);
t1114 = Icges(3,1) * t804;
t1229 = t803 * t916 - t806 * t914 - t670 / 0.2e1 - t689 / 0.2e1 - t691 / 0.2e1 + t797 - t1102 / 0.2e1 + t1114 / 0.2e1;
t765 = Icges(3,2) * t807 + t1111;
t1228 = t803 * (t683 / 0.2e1 + t685 / 0.2e1 - t674 / 0.2e1) - t806 * (t687 / 0.2e1 - t681 / 0.2e1 - t678 / 0.2e1) + t669 / 0.2e1 + t690 / 0.2e1 - t834 / 0.2e1 - t765 / 0.2e1 + t768 / 0.2e1;
t748 = -Icges(3,2) * t1067 - t783;
t749 = t765 * t808;
t906 = -t797 - t1114;
t751 = t906 * t805;
t752 = t906 * t808;
t1227 = (t808 * (t675 - t751) + (-t676 + t752) * t805) * t807 + ((t679 + t748) * t808 + (-t680 + t749) * t805) * t804;
t1224 = 0.2e1 * qJD(1);
t1223 = 0.4e1 * qJD(1);
t1222 = 0.2e1 * qJD(2);
t1220 = 2 * qJD(3);
t1219 = 4 * qJD(3);
t1218 = m(4) / 0.2e1;
t1209 = m(5) * (t232 * t320 + t281 * t391 - t282 * t393);
t1206 = m(5) * t1037 * t393;
t1202 = m(6) * (t737 * t218 + t735 * t219 + t806 * t824);
t1201 = m(6) * (t176 * t280 + t218 * t321 - t219 * t323);
t1199 = m(6) * t1042 * t323;
t526 = (t735 * t808 - t737 * t805) * t804;
t619 = t1068 * t801 + t735 * t807;
t621 = -t1065 * t801 - t807 * t737;
t1197 = m(6) * (t526 * t265 + t280 * t578 - t621 * t468 - t619 * t470 + t806 * t850);
t305 = t737 * t323;
t1043 = -t321 * t735 - t305;
t1195 = m(6) * (t322 * t735 + t1043 + t305);
t1193 = m(6) * (t468 * t734 + t470 * t736 + t520 * t735 + t521 * t737 + (-t265 * t803 + t301 * t806) * t804);
t1185 = -t808 / 0.4e1;
t1126 = rSges(3,1) * t807;
t945 = pkin(1) + t1126;
t976 = rSges(3,2) * t1072 + t808 * rSges(3,3);
t622 = -t805 * t945 + t799 + t976;
t785 = rSges(3,2) * t1070;
t623 = -t785 + t945 * t808 + (rSges(3,3) + pkin(6)) * t805;
t772 = rSges(3,1) * t804 + rSges(3,2) * t807;
t758 = t772 * t805;
t760 = t772 * t808;
t1179 = m(3) * (t622 * t758 - t623 * t760);
t876 = -t561 * t808 - t562 * t805;
t337 = t876 * t807 + (-t1284 * t808 - t643 * t805) * t804;
t436 = t876 * t804;
t1178 = m(4) * (t337 * t436 + t394 * t463 - t395 * t465);
t1172 = m(4) * (t477 * t563 + t478 * t564);
t1171 = m(4) * (-t477 * t605 + t478 * t611);
t1161 = m(5) * (-t1007 * t421 + t420 * t472);
t1160 = m(5) * (t420 * t484 + t421 * t483);
t1157 = m(6) * (t401 * t619 + t402 * t621 + t1043);
t1153 = m(6) * (-t401 * t736 - t402 * t734 + t441 * t737 + t442 * t735);
t846 = t806 * t851;
t1150 = m(6) * (t737 * t452 + t735 * t453 + t846);
t1148 = m(6) * (-t737 * t468 + t735 * t470 + t846);
t1147 = m(6) * (t401 * t441 + t402 * t442);
t1146 = m(6) * (t401 * t452 + t402 * t453);
t1145 = m(6) * (t1073 * t526 + t619 * t736 + t621 * t734);
t1144 = m(6) * (t1071 * t525 + t618 * t737 + t620 * t735);
t1142 = (-t1083 + (-t578 + 0.2e1 * t1066) * t803) * t1254;
t802 = t806 ^ 2;
t1141 = m(6) * (-t734 ^ 2 + t735 ^ 2 - t736 ^ 2 + t737 ^ 2 + (-t800 + t802) * t801);
t1130 = m(6) * qJD(1);
t1129 = m(6) * qJD(2);
t1128 = m(6) * qJD(3);
t1127 = m(6) * qJD(5);
t384 = t402 * t736;
t385 = t402 * t737;
t406 = t421 * t736;
t1079 = t675 * t804;
t236 = -t401 * t734 + t384;
t237 = -t401 * t735 + t385;
t400 = t608 + t825 - t1305;
t1036 = t400 - t401;
t283 = -t420 * t734 + t406;
t422 = t608 + t836 + t981;
t1035 = t420 - t422;
t1009 = t753 * t1072 + t807 * t604;
t1002 = t680 * t1064 + t805 * t672;
t1000 = -t673 + t750;
t998 = -t677 - t747;
t996 = -t682 + t740;
t995 = -t829 - t741;
t994 = t686 + t738;
t993 = t688 + t739;
t991 = -t695 - t778;
t986 = t805 * (pkin(7) * t1067 - t789) + t808 * (-pkin(2) * t1070 + t790);
t975 = qJD(1) * t804;
t974 = qJD(1) * t807;
t973 = qJD(2) * t805;
t971 = qJD(3) * t804;
t970 = qJD(3) * t807;
t506 = (-t578 * t806 + t802 * t807) * t804;
t966 = t506 * t1127;
t955 = -t778 + t989;
t949 = t1072 / 0.4e1;
t939 = t676 * t804 - t671;
t935 = t805 * t667 - t808 * t668 + t986;
t934 = -t778 + t956;
t819 = t804 * t859 + t1086;
t238 = t626 * t734 - t628 * t735 + t805 * t819;
t818 = t804 * t858 + t1085;
t239 = t627 * t734 - t629 * t735 + t805 * t818;
t823 = -t804 * t863 + t1090;
t242 = -t630 * t735 + t634 * t734 + t805 * t823;
t822 = -t804 * t862 + t1089;
t243 = -t631 * t735 + t635 * t734 + t805 * t822;
t821 = -t804 * t861 + t1088;
t244 = t632 * t734 - t636 * t735 + t805 * t821;
t820 = -t804 * t860 + t1087;
t245 = t633 * t734 - t637 * t735 + t805 * t820;
t932 = (t1265 + t1323) * t1256 + ((t239 + t243 + t245) * t808 + (t238 + t242 + t244) * t805 + t1262) * t1191;
t240 = t736 * t626 - t737 * t628 + t808 * t819;
t241 = t736 * t627 - t737 * t629 + t808 * t818;
t246 = -t737 * t630 + t736 * t634 + t808 * t823;
t247 = -t737 * t631 + t736 * t635 + t808 * t822;
t248 = t736 * t632 - t737 * t636 + t808 * t821;
t249 = t736 * t633 - t737 * t637 + t808 * t820;
t931 = (t892 + t891 + t889 + t1264) * t1256 + ((t241 + t247 + t249) * t808 + (t240 + t246 + t248) * t805 + t1261) * t1191;
t930 = (t1268 * t808 + t1269 * t805 + t1260) * t1258 + (t1259 - t1263) * t1188;
t41 = m(5) * (t320 * t525 + t391 * t618 - t393 * t620) + m(6) * (t280 * t525 + t321 * t618 - t323 * t620);
t884 = t468 * t805 + t470 * t808;
t849 = t884 * t804;
t883 = t507 * t805 + t509 * t808;
t92 = m(5) * (t1073 * t883 + t298 * t577) + m(6) * (t265 * t577 + t803 * t849);
t143 = m(5) * t283 + m(6) * t236;
t253 = t1025 * t735 + t1029 * t734 + t1072 * t595;
t254 = t1024 * t735 + t1028 * t734 + t1072 * t596;
t255 = t1023 * t734 + t1027 * t735 + t1072 * t597;
t256 = t1022 * t734 + t1026 * t735 + t1072 * t598;
t261 = t1019 * t734 + t1021 * t735 + t1072 * t593;
t262 = t1018 * t734 + t1020 * t735 + t1072 * t594;
t928 = (t254 + t256 + t262) * t1257 + (t253 + t255 + t261) * t1184;
t257 = t1025 * t737 + t1029 * t736 + t1070 * t595;
t258 = t1024 * t737 + t1028 * t736 + t1070 * t596;
t259 = t1023 * t736 + t1027 * t737 + t1070 * t597;
t260 = t1022 * t736 + t1026 * t737 + t1070 * t598;
t263 = t1019 * t736 + t1021 * t737 + t1070 * t593;
t264 = t1018 * t736 + t1020 * t737 + t1070 * t594;
t927 = (t257 + t259 + t263) * t1255 + (t258 + t260 + t264) * t1190;
t924 = t1260 * t1256 + t1259 * t1258;
t920 = t290 / 0.2e1 + t292 / 0.2e1 + t294 / 0.2e1;
t919 = -t291 / 0.2e1 - t293 / 0.2e1 - t295 / 0.2e1;
t742 = (-Icges(4,5) * t803 - Icges(4,6) * t806) * t804;
t745 = (Icges(6,4) * t806 - Icges(6,5) * t803) * t804;
t746 = (Icges(5,4) * t803 + Icges(5,5) * t806) * t804;
t913 = t746 / 0.2e1 + t745 / 0.2e1 + t742 / 0.2e1;
t901 = -Icges(3,5) * t804 - Icges(3,6) * t807;
t312 = t1005 * t808 - t1282 * t805 + t935;
t864 = t312 + t884;
t842 = -t928 - t932;
t840 = t927 - t931;
t458 = -t1070 * t676 + t1002;
t814 = t458 * t1190 + t1002 * t1257 + (-t649 + (t672 + t1079) * t808 + t1003 + t1281) * t1255;
t813 = (t808 * t939 - t1002 + t458) * t1184 + (-t805 * (-t679 * t807 + t1079) - t671 * t808) * t1255 + (t805 * t939 + t1281 + t940) * t1190;
t346 = t1072 * t745 + t734 * t996 + t735 * t994;
t347 = t1072 * t746 + t734 * t993 + t735 * t995;
t348 = t1070 * t745 + t736 * t996 + t737 * t994;
t349 = t1070 * t746 + t736 * t993 + t737 * t995;
t350 = t1000 * t735 + t1072 * t742 + t734 * t998;
t351 = t1000 * t737 + t1070 * t742 + t736 * t998;
t812 = -t1233 + (t348 + t349 + t351 + t1266) * t805 / 0.4e1 + (t346 + t347 + t350 + t1267) * t1185;
t810 = (-t1072 / 0.4e1 + t949) * (t1309 + t1308 + t1307) + (t1185 + t808 / 0.4e1) * (t1312 + t1311 + t1310) + t1320;
t809 = -t1232 + t1260 * t1191 + t1263 * t1188 + (-t1265 + t1269) * t949 + (-t1264 + t1268) * t1070 / 0.4e1 + (t1262 + t1314) * t1067 / 0.4e1 + (t1261 + t1313) * t1064 / 0.4e1;
t774 = -rSges(3,2) * t804 + t1126;
t744 = t901 * t808;
t743 = t901 * t805;
t584 = t991 * t808;
t582 = t991 * t805;
t565 = t604 * t1070;
t510 = t955 * t808;
t508 = t955 * t805;
t486 = -t1070 * t755 - t807 * t611;
t485 = t1072 * t755 + t605 * t807;
t471 = t934 * t808;
t469 = t934 * t805;
t447 = (t605 * t808 - t611 * t805) * t804;
t434 = -t1284 * t805 + t643 * t808 + t986;
t411 = t1141 / 0.2e1;
t410 = t1142 / 0.2e1;
t409 = -t807 * t746 + (t803 * t993 + t806 * t995) * t804;
t408 = -t807 * t745 + (t803 * t996 + t806 * t994) * t804;
t407 = -t807 * t742 + (t1000 * t806 + t803 * t998) * t804;
t404 = t1007 * t807 + t580 * t804;
t403 = t1072 * t754 + t603 * t807 + t1009;
t399 = 0.4e1 * t1244;
t381 = t1144 / 0.2e1;
t380 = t1145 / 0.2e1;
t370 = -t1283 * t805 + t647 * t808 + t935;
t367 = (qJ(5) * t1074 + t804 * t979) * t808 - t442 * t807;
t366 = -qJ(5) * t618 + t1072 * t757 + t607 * t807 + t1009;
t356 = 0.4e1 * t1245;
t336 = t565 + (t1007 * t805 + t603 * t808) * t804;
t297 = t565 + ((-qJ(5) * t734 + t607) * t808 - t442 * t805) * t804;
t252 = 0.2e1 * t475;
t251 = t474 + t475;
t250 = 0.2e1 * t474;
t225 = 0.2e1 * t433 + t410;
t224 = t432 + t433 - t1142 / 0.2e1;
t223 = 0.2e1 * t432 + t410;
t183 = t381 + t411 - t1145 / 0.2e1;
t182 = t380 + t381 - t1141 / 0.2e1;
t181 = t380 + t411 - t1144 / 0.2e1;
t172 = t1148 / 0.2e1;
t156 = t1150 / 0.2e1;
t154 = t265 * t578 + t806 * t849;
t146 = t1153 / 0.2e1;
t130 = -t248 * t808 + t249 * t805;
t129 = -t246 * t808 + t247 * t805;
t128 = -t244 * t808 + t245 * t805;
t127 = -t242 * t808 + t243 * t805;
t126 = -t240 * t808 + t241 * t805;
t125 = -t238 * t808 + t239 * t805;
t109 = t1157 / 0.2e1;
t107 = -t351 * t807 + (t263 * t805 + t264 * t808) * t804;
t106 = -t350 * t807 + (t261 * t805 + t262 * t808) * t804;
t105 = -t349 * t807 + (t259 * t805 + t260 * t808) * t804;
t104 = -t348 * t807 + (t257 * t805 + t258 * t808) * t804;
t103 = -t347 * t807 + (t255 * t805 + t256 * t808) * t804;
t102 = -t346 * t807 + (t253 * t805 + t254 * t808) * t804;
t101 = t280 * t526 + t321 * t619 - t323 * t621;
t96 = t1193 / 0.2e1;
t72 = t1195 / 0.2e1;
t51 = t1197 / 0.2e1;
t38 = t1230 * t804 - t913 * t807 + t1147 + t1161 + t1171;
t35 = t156 - t1148 / 0.2e1;
t34 = t172 + t156;
t33 = t172 - t1150 / 0.2e1;
t30 = t1228 * t804 + t1229 * t807 + t1146 + t1160 + t1172 + t1179;
t26 = t1202 / 0.2e1;
t22 = t146 + t72 - t1157 / 0.2e1;
t21 = t109 + t146 - t1195 / 0.2e1;
t20 = t109 + t72 - t1153 / 0.2e1;
t18 = t1057 + t1058;
t16 = t96 + t26 - t1197 / 0.2e1;
t15 = t51 + t96 - t1202 / 0.2e1;
t14 = t51 + t26 - t1193 / 0.2e1;
t13 = t1059 - t1272;
t12 = t1059 + t1272;
t11 = t1061 + t1131 - t1059;
t10 = t805 * t927 + t808 * t928 - t1321;
t9 = t1115 - t1270;
t8 = t1115 + t1270;
t7 = t1132 + t1133 - t1115;
t6 = t805 * t813 + t808 * t814;
t5 = t1199 + t1206;
t4 = t1178 + t1209 + t1201 + (t805 * t926 + t808 * t925 + t930) * t807 + (t805 * t932 + t808 * t931 - t924) * t804;
t3 = t812 + t809 - t1320;
t2 = (-t351 / 0.4e1 - t349 / 0.4e1 - t348 / 0.4e1 - t295 / 0.4e1 - t293 / 0.4e1 - t291 / 0.4e1) * t805 + t810 + (t350 / 0.4e1 + t347 / 0.4e1 + t346 / 0.4e1 + t294 / 0.4e1 + t292 / 0.4e1 + t290 / 0.4e1) * t808 + t809 + t1233;
t1 = t812 + t810 + (-t460 / 0.2e1 - t459 / 0.2e1 - t454 / 0.2e1 + (-t318 / 0.4e1 - t317 / 0.4e1 - t314 / 0.4e1 - t278 / 0.4e1 - t276 / 0.4e1 - t274 / 0.4e1) * t808 + (-t316 / 0.4e1 - t315 / 0.4e1 - t313 / 0.4e1 - t277 / 0.4e1 - t275 / 0.4e1 - t273 / 0.4e1) * t805) * t804 + (t365 / 0.2e1 + t364 / 0.2e1 + t363 / 0.2e1 + (-t379 / 0.4e1 - t376 / 0.4e1 - t373 / 0.4e1 - t431 / 0.4e1 - t428 / 0.4e1 - t427 / 0.4e1) * t808 + (-t377 / 0.4e1 - t374 / 0.4e1 + t372 / 0.4e1 - t429 / 0.4e1 - t425 / 0.4e1 - t423 / 0.4e1) * t805) * t807 + t1232;
t17 = [(-m(5) * t1035 * t421 / 0.4e1 + t1036 * t402 * t1213) * t1223 + t30 * qJD(2) + t38 * qJD(3) + t143 * qJD(4) + t237 * t1127, t30 * qJD(1) + t3 * qJD(3) + t18 * qJD(4) + t34 * qJD(5) + ((t477 * t584 + t478 * t582 - t563 * t583 - t564 * t581) * t1218 + (t420 * t510 + t421 * t508 - t483 * t507 - t484 * t509) * t1216 + (t401 * t471 + t402 * t469 - t452 * t470 - t453 * t468) * t1214) * t1222 + (-t313 / 0.2e1 + t764 * t1184 + (t675 / 0.2e1 - t751 / 0.2e1) * t804 - t315 / 0.2e1 + (-t679 / 0.2e1 - t748 / 0.2e1) * t807 - t316 / 0.2e1 + m(3) * (-t622 * t774 - t758 * t772) - t814 - t1235) * t972 + (t317 / 0.2e1 + t314 / 0.2e1 + (t680 / 0.2e1 - t749 / 0.2e1) * t807 + (-t676 / 0.2e1 + t752 / 0.2e1) * t804 + t318 / 0.2e1 + m(3) * (-t623 * t774 + t760 * t772) - t813 + t764 * t1190 - t1234) * t973, t38 * qJD(1) + t3 * qJD(2) + t12 * qJD(4) + t21 * qJD(5) + (-t1199 / 0.4e1 - t1206 / 0.4e1) * t1219 + ((-t463 * t605 - t465 * t611 + t477 * t485 + t478 * t486) * t1218 + (t1007 * t393 + t391 * t472 + t403 * t420 + t404 * t421) * t1216 + (t321 * t441 - t323 * t442 + t366 * t401 + t367 * t402) * t1214) * t1220 + (-t409 - t407 - t408) * t970 + ((t348 / 0.2e1 + t349 / 0.2e1 + t351 / 0.2e1 - t919) * t808 + (t346 / 0.2e1 + t347 / 0.2e1 + t350 / 0.2e1 + t920) * t805) * t971, t143 * qJD(1) + t18 * qJD(2) + t12 * qJD(3) + t251 * qJD(5), t34 * qJD(2) + t21 * qJD(3) + t251 * qJD(4) + t1130 * t237; t6 * qJD(2) + t1 * qJD(3) - t19 * qJD(4) + t33 * qJD(5) + (t1035 * t507 * t1216 + t1036 * t1319) * t1224 + (-t1179 / 0.4e1 - t1172 / 0.4e1 - t1146 / 0.4e1 - t1160 / 0.4e1) * t1223 - t1229 * t974 - t1228 * t975, t6 * qJD(1) + t10 * qJD(3) + t92 * qJD(4) + t154 * t1127 - (t128 + t127 + t125 + t1225 * t743 + (-t808 * t744 + t1227) * t805) * t972 / 0.2e1 + (m(6) * (t265 * t312 - t468 * t469 - t470 * t471) + m(5) * (t298 * t370 - t507 * t508 - t509 * t510) + m(4) * (t396 * t434 - t581 * t582 - t583 * t584) + m(3) * ((t805 * (rSges(3,1) * t1067 - t976) + t808 * (rSges(3,1) * t1064 + t805 * rSges(3,3) - t785)) * (-t805 * t758 - t760 * t808) + t1302 * t774 * t772) + (t130 + t129 + t126 + t1226 * t744 + (-t805 * t743 + t1227) * t808) * t1190) * qJD(2), t1 * qJD(1) + t10 * qJD(2) + t8 * qJD(4) + t15 * qJD(5) + (-t1201 / 0.4e1 - t1209 / 0.4e1 - t1178 / 0.4e1) * t1219 + ((t447 * t396 + t436 * t451 - t485 * t583 - t486 * t581 + (-t463 * t808 + t465 * t805) * t755) * t1218 + (t298 * t336 + t320 * t355 + t391 * t580 - t393 * t579 - t403 * t509 - t404 * t507) * t1216 + (t265 * t297 + t280 * t301 + t321 * t521 - t323 * t520 - t366 * t470 - t367 * t468) * t1214) * t1220 + ((t920 - t925) * t808 + (t919 - t926) * t805 - t930) * t970 + (t805 * t842 + t808 * t840 + t924) * t971 + ((-t106 / 0.2e1 - t102 / 0.2e1 - t103 / 0.2e1) * t808 + (t104 / 0.2e1 + t105 / 0.2e1 + t107 / 0.2e1) * t805) * qJD(3), t92 * qJD(2) + t8 * qJD(3) - 0.4e1 * qJD(4) * t1244 + t224 * qJD(5) - t1322, t33 * qJD(1) + t15 * qJD(3) + t224 * qJD(4) + t1129 * t154 - t966; t2 * qJD(2) + t5 * qJD(3) + t11 * qJD(4) + t20 * qJD(5) + t913 * t974 + (-t1171 / 0.4e1 - t1147 / 0.4e1 - t1161 / 0.4e1) * t1223 + ((t1035 * t393 - t1037 * t421) * t1216 + (-t1036 * t323 - t1042 * t402) * t1214) * t1224 - t1230 * t975, t2 * qJD(1) + t4 * qJD(3) + t7 * qJD(4) + t14 * qJD(5) + ((t337 * t396 - t394 * t583 - t395 * t581 + t434 * t436 + t463 * t584 - t465 * t582) * t1218 + (t232 * t298 - t281 * t509 - t282 * t507 + t320 * t370 + t391 * t510 - t393 * t508) * t1216 + (t176 * t265 - t218 * t470 - t219 * t468 + t280 * t312 + t321 * t471 - t323 * t469) * t1214) * t1222 + ((t1307 / 0.2e1 + t1308 / 0.2e1 + t1309 / 0.2e1 + t1235) * t807 + t842) * t972 + (-t840 + (t1324 * t1190 + t1325 * t1255 + t1234) * t807) * t973 + (((t126 / 0.2e1 + t129 / 0.2e1 + t130 / 0.2e1 + t372 / 0.2e1 - t374 / 0.2e1 - t377 / 0.2e1) * t808 + (t125 / 0.2e1 + t127 / 0.2e1 + t128 / 0.2e1 + t379 / 0.2e1 + t376 / 0.2e1 + t373 / 0.2e1) * t805) * t804 + t1321) * qJD(2), t5 * qJD(1) + t4 * qJD(2) + t41 * qJD(4) + (t409 / 0.2e1 + t407 / 0.2e1 + t408 / 0.2e1) * qJD(3) * t807 ^ 2 + ((t280 * t297 + t321 * t366 - t323 * t367) * t1213 + (t320 * t336 + t391 * t403 - t393 * t404) * t1215 + (t436 * t447 + t463 * t485 - t465 * t486) * m(4) / 0.4e1) * t1219 + t101 * t1127 + ((t103 + t106 + t102) * t1190 + (t1266 * t808 + t1267 * t805) * t1188 + (t107 + t104 + t105) * t1184) * t971, t11 * qJD(1) + t7 * qJD(2) + t41 * qJD(3) + (-0.4e1 * t1245 + 0.2e1 * (t1214 + t1216) * (t1073 * t525 + t618 * t736 + t620 * t734)) * qJD(4) + t182 * qJD(5), t20 * qJD(1) + t14 * qJD(2) + t101 * t1128 + t182 * qJD(4) + (t1071 * t526 + t619 * t737 + t621 * t735 + t461) * t1127; (m(6) * (t400 * t734 + t236 - t384) + m(5) * (t422 * t734 + t283 - t406) - t143) * qJD(1) + t19 * qJD(2) + t13 * qJD(3) + t250 * qJD(5), t1322 + (m(6) * (t734 * t469 + t736 * t471) + m(5) * (t734 * t508 + t736 * t510) + 0.2e1 * ((t265 * t1214 + t1216 * t298) * t807 + (t864 * t1214 + (t370 + t883) * t1216) * t804) * t803 - t92) * qJD(2) + t9 * qJD(3) + t399 * qJD(4) + t223 * qJD(5), t13 * qJD(1) + t9 * qJD(2) + (m(6) * (t321 * t737 - t323 * t735 + t366 * t736 + t367 * t734 + (t280 * t806 + t297 * t803) * t804) + m(5) * (t391 * t737 - t393 * t735 + t403 * t736 + t404 * t734 + (t320 * t806 + t336 * t803) * t804) - t41) * qJD(3) + t356 * qJD(4) + t181 * qJD(5), qJD(2) * t399 + qJD(3) * t356, qJD(1) * t250 + qJD(2) * t223 + qJD(3) * t181; (t400 * t735 - t385) * t1130 + t35 * qJD(2) + t22 * qJD(3) + t252 * qJD(4), t35 * qJD(1) + (t735 * t469 + t737 * t471 - t154 + (t265 * t807 + t804 * t864) * t806) * t1129 + t16 * qJD(3) + t225 * qJD(4) + t966, t22 * qJD(1) + t16 * qJD(2) + (-t321 * t736 + t323 * t734 + t366 * t737 + t367 * t735 + (-t280 * t803 + t297 * t806) * t804 - t101) * t1128 + t183 * qJD(4) - t461 * t1127, qJD(1) * t252 + qJD(2) * t225 + qJD(3) * t183, 0.4e1 * (t506 * qJD(2) / 0.4e1 - t461 * qJD(3) / 0.4e1) * m(6);];
Cq = t17;
