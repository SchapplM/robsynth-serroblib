% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR10_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:40
% EndTime: 2019-03-09 04:08:51
% DurationCPUTime: 59.48s
% Computational Cost: add. (159312->1358), mult. (166444->1890), div. (0->0), fcn. (180317->10), ass. (0->821)
t874 = sin(qJ(3));
t877 = cos(qJ(1));
t1160 = t874 * t877;
t876 = cos(qJ(3));
t1147 = t876 * t877;
t875 = sin(qJ(1));
t1149 = t875 * t876;
t867 = pkin(10) + qJ(5);
t858 = qJ(6) + t867;
t853 = sin(t858);
t1155 = t875 * t853;
t854 = cos(t858);
t749 = -t1155 * t874 + t854 * t877;
t1154 = t875 * t854;
t750 = t1154 * t874 + t853 * t877;
t588 = Icges(7,5) * t750 + Icges(7,6) * t749 - Icges(7,3) * t1149;
t1227 = Icges(7,4) * t750;
t591 = Icges(7,2) * t749 - Icges(7,6) * t1149 + t1227;
t736 = Icges(7,4) * t749;
t594 = Icges(7,1) * t750 - Icges(7,5) * t1149 + t736;
t751 = t1160 * t853 + t1154;
t752 = t1160 * t854 - t1155;
t359 = t588 * t1147 + t751 * t591 - t752 * t594;
t1452 = t359 * t877;
t856 = sin(t867);
t1153 = t875 * t856;
t857 = cos(t867);
t772 = -t874 * t1153 + t857 * t877;
t1152 = t875 * t857;
t773 = t1152 * t874 + t856 * t877;
t605 = Icges(6,5) * t773 + Icges(6,6) * t772 - Icges(6,3) * t1149;
t1230 = Icges(6,4) * t773;
t608 = Icges(6,2) * t772 - Icges(6,6) * t1149 + t1230;
t757 = Icges(6,4) * t772;
t611 = Icges(6,1) * t773 - Icges(6,5) * t1149 + t757;
t774 = t1160 * t856 + t1152;
t775 = t1160 * t857 - t1153;
t381 = t605 * t1147 + t774 * t608 - t775 * t611;
t1451 = t381 * t877;
t1450 = t875 * t359;
t1449 = t875 * t381;
t607 = -Icges(6,5) * t775 + Icges(6,6) * t774 + Icges(6,3) * t1147;
t1191 = t607 * t874;
t759 = Icges(6,4) * t775;
t610 = Icges(6,2) * t774 + Icges(6,6) * t1147 - t759;
t758 = Icges(6,4) * t774;
t612 = Icges(6,1) * t775 - Icges(6,5) * t1147 - t758;
t1433 = t610 * t856 + t612 * t857;
t417 = t1433 * t876 - t1191;
t590 = -Icges(7,5) * t752 + Icges(7,6) * t751 + Icges(7,3) * t1147;
t1194 = t590 * t874;
t738 = Icges(7,4) * t752;
t593 = Icges(7,2) * t751 + Icges(7,6) * t1147 - t738;
t737 = Icges(7,4) * t751;
t595 = Icges(7,1) * t752 - Icges(7,5) * t1147 - t737;
t1434 = t593 * t853 + t595 * t854;
t402 = t1434 * t876 - t1194;
t1248 = m(7) * qJD(1);
t1448 = -t1248 / 0.4e1;
t1447 = -m(6) * qJD(1) / 0.4e1;
t871 = sin(pkin(10));
t1151 = t875 * t871;
t872 = cos(pkin(10));
t855 = t872 * pkin(4) + pkin(3);
t873 = -pkin(8) - qJ(4);
t1019 = pkin(4) * t1151 - t873 * t1147 - t855 * t1160;
t813 = pkin(5) * t857 + t855;
t866 = -pkin(9) + t873;
t1070 = t866 * t1147 + t813 * t1160;
t1254 = pkin(5) * t856;
t1255 = pkin(4) * t871;
t820 = t1254 + t1255;
t557 = t875 * t820 - t1019 - t1070;
t598 = t752 * rSges(7,1) - t751 * rSges(7,2) - rSges(7,3) * t1147;
t1108 = t557 - t598;
t1067 = t813 - t855;
t1064 = t866 - t873;
t989 = t1064 * t874;
t662 = t1067 * t876 - t989;
t634 = t662 * t1147;
t1242 = rSges(7,1) * t854;
t972 = -rSges(7,2) * t853 + t1242;
t723 = rSges(7,3) * t874 + t876 * t972;
t700 = t723 * t1147;
t373 = -t1108 * t874 + t634 + t700;
t1438 = t874 * t598 + t700;
t374 = -t557 * t874 + t1438 + t634;
t1125 = t373 - t374;
t1446 = m(7) * t1125;
t1161 = t874 * t875;
t795 = -t1151 * t874 + t872 * t877;
t1150 = t875 * t872;
t1169 = t871 * t877;
t796 = t1150 * t874 + t1169;
t647 = Icges(5,4) * t796 + Icges(5,2) * t795 - Icges(5,6) * t1149;
t650 = Icges(5,1) * t796 + Icges(5,4) * t795 - Icges(5,5) * t1149;
t1232 = Icges(4,4) * t874;
t965 = Icges(4,2) * t876 + t1232;
t778 = Icges(4,6) * t877 + t875 * t965;
t843 = Icges(4,4) * t1149;
t780 = Icges(4,1) * t1161 + Icges(4,5) * t877 + t843;
t960 = Icges(4,5) * t874 + Icges(4,6) * t876;
t1445 = t795 * t647 + t796 * t650 + t778 * t1149 + t780 * t1161 + t877 * (Icges(4,3) * t877 + t875 * t960);
t1111 = t749 * t591 + t750 * t594;
t938 = -t751 * t593 - t752 * t595;
t1444 = t1111 + (-t588 * t875 - t590 * t877) * t876 + t938;
t1107 = t772 * t608 + t773 * t611;
t934 = -t774 * t610 - t775 * t612;
t1443 = t1107 + (-t605 * t875 - t607 * t877) * t876 + t934;
t1051 = 0.2e1 * m(7);
t1430 = -t1051 / 0.4e1;
t992 = t1051 / 0.4e1;
t1010 = -t1160 / 0.4e1;
t1432 = 0.2e1 * m(6);
t1442 = -t1432 / 0.4e1;
t1441 = t1432 / 0.4e1;
t1362 = t1064 * t876 + t1067 * t874;
t1440 = t1362 * t875;
t797 = t1160 * t871 + t1150;
t798 = t1160 * t872 - t1151;
t649 = -Icges(5,4) * t798 + Icges(5,2) * t797 + Icges(5,6) * t1147;
t652 = -Icges(5,1) * t798 + Icges(5,4) * t797 + Icges(5,5) * t1147;
t1103 = t795 * t649 + t796 * t652;
t922 = -t778 * t876 - t780 * t874;
t1439 = t877 * t922 - t1103;
t615 = t775 * rSges(6,1) - t774 * rSges(6,2) - rSges(6,3) * t1147;
t1243 = rSges(6,1) * t857;
t973 = -rSges(6,2) * t856 + t1243;
t900 = t973 * t876;
t745 = rSges(6,3) * t874 + t900;
t1437 = t745 * t1147 + t615 * t874;
t1187 = t652 * t872;
t1188 = t650 * t872;
t1189 = t649 * t871;
t1190 = t647 * t871;
t959 = Icges(5,5) * t872 - Icges(5,6) * t871;
t765 = Icges(5,3) * t874 + t876 * t959;
t706 = t765 * t875;
t707 = t765 * t877;
t779 = -Icges(4,6) * t875 + t877 * t965;
t1231 = Icges(4,4) * t876;
t969 = Icges(4,1) * t874 + t1231;
t781 = -Icges(4,5) * t875 + t877 * t969;
t801 = -Icges(4,2) * t1161 + t843;
t824 = -Icges(4,2) * t874 + t1231;
t802 = t824 * t877;
t826 = Icges(4,1) * t876 - t1232;
t803 = t826 * t875;
t804 = t826 * t877;
t644 = Icges(5,5) * t796 + Icges(5,6) * t795 - Icges(5,3) * t1149;
t646 = -Icges(5,5) * t798 + Icges(5,6) * t797 + Icges(5,3) * t1147;
t929 = t644 * t877 + t646 * t875;
t1436 = ((-t780 - t801 - t1188 + t706 + t1190) * t877 + (t781 + t802 - t1187 - t707 + t1189) * t875) * t876 + (-(t779 - t804) * t875 + (t778 - t803) * t877 - t929) * t874;
t627 = t749 * rSges(7,1) - t750 * rSges(7,2);
t760 = t772 * pkin(5);
t1090 = -t627 - t760;
t628 = t751 * rSges(7,1) + t752 * rSges(7,2);
t761 = t774 * pkin(5);
t1417 = -t628 - t761;
t1435 = t1090 * t877 + t1417 * t875;
t380 = -t1149 * t607 + t772 * t610 - t612 * t773;
t358 = -t1149 * t590 + t749 * t593 - t595 * t750;
t1337 = m(5) / 0.4e1;
t1336 = m(6) / 0.4e1;
t1335 = m(7) / 0.4e1;
t1316 = t874 / 0.2e1;
t1314 = t875 / 0.2e1;
t1313 = t875 / 0.4e1;
t1312 = t876 / 0.2e1;
t1311 = t877 / 0.2e1;
t1310 = t877 / 0.4e1;
t957 = Icges(7,5) * t854 - Icges(7,6) * t853;
t717 = Icges(7,3) * t874 + t876 * t957;
t1225 = Icges(7,4) * t854;
t962 = -Icges(7,2) * t853 + t1225;
t719 = Icges(7,6) * t874 + t876 * t962;
t1226 = Icges(7,4) * t853;
t966 = Icges(7,1) * t854 - t1226;
t721 = Icges(7,5) * t874 + t876 * t966;
t459 = t1147 * t717 + t751 * t719 - t752 * t721;
t1426 = t459 * t874;
t360 = t1147 * t590 - t938;
t954 = t360 * t877 - t1450;
t188 = t876 * t954 + t1426;
t357 = -t1149 * t588 + t1111;
t1134 = -t357 + t1444;
t68 = -t1426 + (t1134 * t877 + t1450) * t876;
t1428 = t188 + t68;
t220 = t360 * t875 + t1452;
t89 = t1134 * t875 - t1452;
t1427 = t220 + t89;
t958 = Icges(6,5) * t857 - Icges(6,6) * t856;
t731 = Icges(6,3) * t874 + t876 * t958;
t1228 = Icges(6,4) * t857;
t963 = -Icges(6,2) * t856 + t1228;
t733 = Icges(6,6) * t874 + t876 * t963;
t1229 = Icges(6,4) * t856;
t967 = Icges(6,1) * t857 - t1229;
t735 = Icges(6,5) * t874 + t876 * t967;
t471 = t1147 * t731 + t774 * t733 - t775 * t735;
t1425 = t471 * t874;
t1424 = t745 * t877;
t679 = t719 * t875;
t681 = t721 * t875;
t939 = -t591 * t853 + t594 * t854;
t910 = t717 * t875 - t939;
t294 = (-t679 * t853 + t681 * t854 + t588) * t876 + t910 * t874;
t1185 = t717 * t874;
t716 = Icges(7,3) * t876 - t874 * t957;
t1177 = t854 * t721;
t1180 = t853 * t719;
t924 = t1177 - t1180;
t904 = t716 - t924;
t1371 = t876 * t904 - t1185;
t718 = Icges(7,6) * t876 - t874 * t962;
t720 = Icges(7,5) * t876 - t874 * t966;
t334 = -t1371 * t875 + t718 * t749 + t720 * t750;
t1423 = t294 + t334;
t680 = t719 * t877;
t682 = t721 * t877;
t909 = -t717 * t877 + t1434;
t295 = (t680 * t853 - t682 * t854 + t590) * t876 + t909 * t874;
t335 = t1371 * t877 + t751 * t718 - t752 * t720;
t1422 = t295 + t335;
t1196 = t588 * t874;
t401 = t876 * t939 + t1196;
t457 = -t1149 * t717 + t719 * t749 + t721 * t750;
t1421 = t401 + t457;
t1420 = -t402 + t459;
t777 = -Icges(4,3) * t875 + t877 * t960;
t550 = -t779 * t1149 - t781 * t1161 - t877 * t777;
t1419 = -t1149 * t646 + t1103 + t550;
t1079 = t721 + (-Icges(7,2) * t854 - t1226) * t876;
t1240 = rSges(5,3) * t876;
t1256 = pkin(3) * t874;
t1383 = t798 * rSges(5,1) - t797 * rSges(5,2);
t842 = qJ(4) * t1147;
t1334 = pkin(1) + pkin(7);
t860 = t877 * qJ(2);
t981 = -t1334 * t875 + t860;
t566 = (-t1240 + t1256) * t877 - t842 + t981 + t1383;
t1407 = -rSges(5,1) * t872 + rSges(5,2) * t871;
t771 = t874 * rSges(5,3) - t1407 * t876;
t831 = pkin(3) * t876 + qJ(4) * t874;
t675 = (t771 + t831) * t877;
t1186 = t675 * t875;
t1142 = -qJ(4) - t873;
t1253 = pkin(3) - t855;
t729 = t1142 * t874 - t1253 * t876;
t1077 = t729 + t831;
t573 = (t745 + t1077) * t877;
t1197 = t573 * t875;
t1085 = t662 + t723;
t483 = (t1077 + t1085) * t877;
t1203 = t483 * t875;
t474 = -0.2e1 * t1203;
t555 = -0.2e1 * t1197;
t659 = -0.2e1 * t1186;
t1021 = (t474 + 0.2e1 * t1203) * t1335 + (t555 + 0.2e1 * t1197) * t1336 + (t659 + 0.2e1 * t1186) * t1337;
t1340 = 0.2e1 * t877;
t812 = t875 * t831;
t1081 = t875 * t729 + t812;
t990 = t1085 * t875;
t481 = t990 + t1081;
t1182 = t745 * t875;
t571 = t1081 + t1182;
t673 = t771 * t875 + t812;
t1022 = (t1340 * t481 + t474) * t1335 + (t1340 * t571 + t555) * t1336 + (t1340 * t673 + t659) * t1337;
t1412 = t1021 - t1022;
t1042 = -0.2e1 * t1147;
t1044 = -0.2e1 * t1160;
t1045 = 0.2e1 * t1161;
t840 = qJ(4) * t1149;
t915 = t796 * rSges(5,1) + t795 * rSges(5,2) - rSges(5,3) * t1149;
t567 = -t840 + t1334 * t877 + (qJ(2) + t1256) * t875 + t915;
t1110 = t567 * t1044 + t566 * t1045;
t517 = -t1019 + t981 + t615;
t1168 = t873 * t876;
t614 = t773 * rSges(6,1) + t772 * rSges(6,2) - rSges(6,3) * t1149;
t518 = (t1255 + t1334) * t877 + (t855 * t874 + qJ(2) + t1168) * t875 + t614;
t1117 = t518 * t1044 + t517 * t1045;
t1050 = t820 + t1334;
t487 = -t1050 * t875 + t1070 + t598 + t860;
t597 = t750 * rSges(7,1) + t749 * rSges(7,2) - rSges(7,3) * t1149;
t488 = t1050 * t877 + (t813 * t874 + t866 * t876 + qJ(2)) * t875 + t597;
t1119 = t488 * t1044 + t487 * t1045;
t1034 = t483 * t1149;
t468 = 0.2e1 * t1034;
t1033 = t573 * t1149;
t548 = 0.2e1 * t1033;
t1032 = t675 * t1149;
t653 = 0.2e1 * t1032;
t1024 = (t1042 * t481 + t1119 + t468) * t1335 + (t1042 * t571 + t1117 + t548) * t1336 + (t1042 * t673 + t1110 + t653) * t1337;
t1026 = (-0.2e1 * t1034 + t468) * t1335 + (-0.2e1 * t1033 + t548) * t1336 + (-0.2e1 * t1032 + t653) * t1337;
t1411 = t1024 - t1026;
t556 = (t820 - t1255) * t877 + t1440;
t1109 = t556 + t597;
t372 = t1109 * t874 + t876 * t990;
t1211 = t372 * t875;
t228 = t373 * t1340 + 0.2e1 * t1211;
t530 = t1149 * t745 + t614 * t874;
t1199 = t530 * t875;
t415 = t1340 * t1437 + 0.2e1 * t1199;
t1132 = t1335 * t228 + t1336 * t415;
t1198 = t1437 * t877;
t1209 = t374 * t877;
t1140 = (-0.2e1 * t1209 + t228 - 0.2e1 * t1211) * t1335 + (-0.2e1 * t1198 + t415 - 0.2e1 * t1199) * t1336;
t1410 = t1132 - t1140;
t1043 = -0.2e1 * t1149;
t226 = t373 * t1042 + t372 * t1043;
t413 = t1042 * t1437 + t530 * t1043;
t1133 = t1335 * t226 + t1336 * t413;
t1047 = 0.2e1 * t876;
t1141 = ((t1209 + t1211) * t1047 + t226) * t1335 + ((t1198 + t1199) * t1047 + t413) * t1336;
t1409 = t1133 - t1141;
t1113 = t530 * t1044 + t1045 * t1437;
t1126 = t372 * t1044 + t373 * t1045;
t868 = t875 ^ 2;
t870 = t877 ^ 2;
t1063 = t868 + t870;
t815 = t1063 * t876;
t1343 = 0.2e1 * t815;
t319 = (-t1108 * t875 - t1109 * t877) * t876;
t932 = t614 * t877 - t615 * t875;
t478 = t932 * t876;
t1139 = (t1343 * t319 + t1126) * t1335 + (-t1343 * t478 + t1113) * t1336;
t1341 = 0.2e1 * t874;
t1388 = t723 * t877;
t818 = t855 * t1147;
t1091 = t818 + (-t813 * t876 + t989) * t877 - t1388;
t787 = t813 * t1149;
t816 = t855 * t1149;
t617 = -t875 * t989 + t787 - t816;
t1179 = t853 * t876;
t1039 = rSges(7,2) * t1179;
t1068 = rSges(7,3) * t1161 + t1149 * t1242;
t683 = -t1039 * t875 + t1068;
t1092 = -t617 - t683;
t1105 = t597 * t1160 - t1161 * t598;
t207 = (t556 * t877 + t557 * t875) * t874 + (-t1091 * t875 + t1092 * t877) * t876 + t1105;
t1173 = t856 * t876;
t1040 = rSges(6,2) * t1173;
t1066 = rSges(6,3) * t1161 + t1149 * t1243;
t701 = -t1040 * t875 + t1066;
t407 = (t1424 * t875 - t701 * t877) * t876 + t932 * t874;
t744 = rSges(6,3) * t876 - t874 * t973;
t441 = (t744 * t877 + t615) * t876;
t1204 = t441 * t875;
t1183 = t744 * t875;
t440 = (t614 + t1183) * t876 + (t701 - t1182) * t874;
t1205 = t440 * t877;
t948 = t1204 - t1205;
t722 = rSges(7,3) * t876 - t874 * t972;
t1020 = t722 * t1149 + t876 * t597 + t874 * t683;
t255 = (t556 - t1440) * t876 + (t617 - t990) * t874 + t1020;
t699 = t722 * t1147;
t256 = t699 + (-t1362 * t877 - t1108) * t876 + (-t1085 * t877 - t1091) * t874;
t956 = -t255 * t877 + t256 * t875;
t1233 = (t407 * t1341 + (-t478 - t948) * t1047 + t1113) * t1336 + (t207 * t1341 + (t319 - t956) * t1047 + t1126) * t1335;
t1408 = t1139 - t1233;
t642 = rSges(6,1) * t772 - rSges(6,2) * t773;
t643 = rSges(6,1) * t774 + rSges(6,2) * t775;
t756 = (-rSges(7,1) * t853 - rSges(7,2) * t854) * t876;
t983 = pkin(5) * t1173 - t756;
t664 = t983 * t875;
t665 = t983 * t877;
t785 = (-rSges(6,1) * t856 - rSges(6,2) * t857) * t876;
t1405 = (t1090 * t483 + t1417 * t481 - t487 * t664 + t488 * t665) * t992 + (-t571 * t643 - t573 * t642 + (t517 * t875 - t518 * t877) * t785) * t1441;
t574 = t787 + (-t866 * t874 - t1039) * t875 + t1068;
t575 = ((rSges(7,3) - t866) * t874 + (t813 + t972) * t876) * t877;
t601 = t816 + (-t873 * t874 - t1040) * t875 + t1066;
t602 = t818 + (t900 + (rSges(6,3) - t873) * t874) * t877;
t1404 = (t1437 * t602 + t440 * t518 + t441 * t517 + t530 * t601) * t1442 + (t255 * t488 + t256 * t487 + t372 * t574 + t373 * t575) * t1430;
t1402 = -t483 * t1446 / 0.2e1;
t753 = (-Icges(7,5) * t853 - Icges(7,6) * t854) * t876;
t1163 = t874 * t753;
t1401 = t1163 / 0.2e1 - t1079 * t1179 / 0.2e1;
t1397 = t868 / 0.2e1;
t1396 = -t874 / 0.2e1;
t1315 = -t875 / 0.2e1;
t930 = t642 * t877 + t875 * t643;
t1394 = m(6) * t930;
t1393 = m(7) * t1047;
t1392 = t358 * t875;
t1391 = t358 * t877;
t1390 = t380 * t875;
t1389 = t380 * t877;
t1387 = (t779 * t876 + t781 * t874) * t877;
t1046 = 0.4e1 * t876;
t545 = (t1337 + t1336 + t1335) * (0.1e1 - t1063) * t874 * t1046;
t657 = (rSges(5,3) + qJ(4)) * t1161 + (pkin(3) - t1407) * t1149;
t1381 = t1253 * t874;
t1380 = (t566 * t877 + t567 * t875) * t1337;
t1379 = (t517 * t877 + t518 * t875) * t1336;
t1378 = (t487 * t877 + t488 * t875) * t1335;
t1178 = t854 * t720;
t1181 = t853 * t718;
t377 = (t1178 + t717 - t1181) * t876 + t904 * t874;
t497 = t876 * t924 + t1185;
t951 = -t401 * t875 - t402 * t877;
t1374 = (t497 * t874 + t876 * t951) * t1312 + ((-t294 * t875 + t295 * t877 + t497) * t876 + (t377 - t951) * t874) * t1316;
t1100 = -Icges(7,2) * t750 + t594 + t736;
t1102 = -Icges(7,1) * t749 + t1227 + t591;
t620 = Icges(7,5) * t749 - Icges(7,6) * t750;
t262 = t1100 * t749 - t1102 * t750 - t1149 * t620;
t1099 = Icges(7,2) * t752 - t595 + t737;
t1101 = -Icges(7,1) * t751 + t593 - t738;
t621 = Icges(7,5) * t751 + Icges(7,6) * t752;
t263 = t1099 * t749 - t1101 * t750 - t1149 * t621;
t150 = t262 * t877 + t263 * t875;
t264 = t1100 * t751 + t1102 * t752 + t1147 * t620;
t265 = t1099 * t751 + t1101 * t752 + t1147 * t621;
t151 = t264 * t877 + t265 * t875;
t1138 = t1311 * t150 + t1314 * t151;
t1373 = t497 * t1312 + t1316 * t377;
t1120 = -t1435 * t1430 - t1394 / 0.2e1;
t1121 = -t1435 * t1393 / 0.4e1 + t1047 * t1394 / 0.4e1;
t604 = t877 * t628;
t434 = t1090 * t875 + t761 * t877 + t604;
t513 = -t875 * t642 + t643 * t877;
t926 = -t664 * t875 - t665 * t877;
t970 = t1063 * t1047;
t1127 = (-t1047 * t926 + t1341 * t434) * t1335 + (t1341 * t513 - t785 * t970) * t1336;
t931 = t627 * t877 + t875 * t628;
t489 = t931 * t876;
t537 = t756 * t1149 + t874 * t627;
t713 = t756 * t1147;
t538 = -t874 * t628 + t713;
t978 = -t319 * t489 + t372 * t537 + t373 * t538;
t1184 = t731 * t874;
t730 = Icges(6,3) * t876 - t874 * t958;
t1171 = t857 * t735;
t1174 = t856 * t733;
t923 = t1171 - t1174;
t903 = t730 - t923;
t1372 = t876 * t903 - t1184;
t907 = -t731 * t877 + t1433;
t1370 = t876 * t907 - t1191;
t1193 = t605 * t874;
t935 = -t608 * t856 + t611 * t857;
t908 = t731 * t875 - t935;
t1369 = t876 * t908 - t1193;
t1368 = t876 * t909 - t1194;
t1367 = t876 * t910 - t1196;
t964 = Icges(5,4) * t872 - Icges(5,2) * t871;
t767 = Icges(5,6) * t874 + t876 * t964;
t968 = Icges(5,1) * t872 - Icges(5,4) * t871;
t769 = Icges(5,5) * t874 + t876 * t968;
t1318 = -t872 / 0.2e1;
t1319 = t871 / 0.2e1;
t976 = rSges(4,1) * t874 + rSges(4,2) * t876;
t894 = -t875 * rSges(4,3) + t877 * t976;
t690 = t894 + t981;
t691 = (rSges(4,3) + t1334) * t877 + (qJ(2) + t976) * t875;
t832 = rSges(4,1) * t876 - rSges(4,2) * t874;
t807 = t832 * t875;
t809 = t832 * t877;
t1366 = -m(4) * (t690 * t809 + t691 * t807) - (t767 * t1319 + t769 * t1318 - t1171 / 0.2e1 + t1174 / 0.2e1 - t1177 / 0.2e1 + t1180 / 0.2e1 + t730 / 0.2e1 + t716 / 0.2e1 + t965 / 0.2e1 - t826 / 0.2e1 + Icges(5,3) * t1312 + t959 * t1396) * t874;
t1065 = -pkin(3) * t1147 - qJ(4) * t1160;
t658 = t771 * t877 - t1065;
t927 = t657 * t877 - t875 * t658;
t936 = t601 * t877 - t875 * t602;
t940 = t574 * t877 - t875 * t575;
t1023 = (t1047 * t940 + t1119) * t1335 + (t1047 * t936 + t1117) * t1336 + (t1047 * t927 + t1110) * t1337;
t1041 = t876 ^ 2 * t1254;
t475 = -t1041 * t875 + t760 * t874 + t537;
t476 = -t1041 * t877 + t1417 * t874 + t713;
t553 = t1149 * t785 + t642 * t874;
t554 = t1147 * t785 - t874 * t643;
t1364 = t875 * (m(6) * t554 + m(7) * t476) - t877 * (m(6) * t553 + m(7) * t475);
t783 = (-Icges(6,2) * t857 - t1229) * t876;
t784 = (-Icges(6,1) * t856 - t1228) * t876;
t1363 = -(t735 / 0.2e1 + t783 / 0.2e1) * t856 + t857 * (t784 / 0.2e1 - t733 / 0.2e1);
t1095 = -Icges(6,2) * t773 + t611 + t757;
t1097 = -Icges(6,1) * t772 + t1230 + t608;
t635 = Icges(6,5) * t772 - Icges(6,6) * t773;
t324 = t635 * t874 + (-t1095 * t856 - t1097 * t857) * t876;
t1094 = Icges(6,2) * t775 - t612 + t758;
t1096 = -Icges(6,1) * t774 + t610 - t759;
t636 = Icges(6,5) * t774 + Icges(6,6) * t775;
t325 = t636 * t874 + (-t1094 * t856 - t1096 * t857) * t876;
t1075 = t735 + t783;
t1076 = t733 - t784;
t782 = (-Icges(6,5) * t856 - Icges(6,6) * t857) * t876;
t394 = t1075 * t772 - t1076 * t773 - t1149 * t782;
t395 = t1075 * t774 + t1076 * t775 + t1147 * t782;
t1359 = t1405 + (t325 + t395) * t1313 + (t324 + t394) * t1310;
t308 = t621 * t874 + (-t1099 * t853 - t1101 * t854) * t876;
t1214 = t308 * t875;
t307 = t620 * t874 + (-t1100 * t853 - t1102 * t854) * t876;
t1215 = t307 * t877;
t755 = (-Icges(7,1) * t853 - t1225) * t876;
t1080 = t719 - t755;
t368 = t1079 * t749 - t1080 * t750 - t1149 * t753;
t369 = t1079 * t751 + t1080 * t752 + t1147 * t753;
t987 = t1215 / 0.4e1 + t1214 / 0.4e1 + t369 * t1313 + t368 * t1310;
t920 = -t807 * t877 + t875 * t809;
t986 = t940 * t1430 + t936 * t1442 - m(5) * t927 / 0.2e1 + m(4) * t920 / 0.2e1;
t1145 = t877 * t188;
t455 = t457 * t874;
t955 = -t357 * t875 + t1391;
t187 = t876 * t955 + t455;
t1158 = t875 * t187;
t1137 = t360 + t1444;
t67 = t455 + (-t1137 * t875 + t1391) * t876;
t1357 = t1145 / 0.4e1 - t1158 / 0.4e1 + t68 * t1310 + t67 * t1313;
t382 = t1147 * t607 - t934;
t1131 = t382 + t1443;
t100 = t1131 * t877 + t1390;
t1001 = t1147 / 0.4e1;
t1003 = -t1147 / 0.4e1;
t1007 = -t1149 / 0.4e1;
t379 = -t1149 * t605 + t1107;
t1128 = -t379 + t1443;
t101 = t1128 * t875 - t1451;
t952 = t382 * t877 - t1449;
t194 = t876 * t952 + t1425;
t1144 = t877 * t194;
t469 = -t1149 * t731 + t733 * t772 + t735 * t773;
t467 = t469 * t874;
t953 = -t379 * t875 + t1389;
t193 = t876 * t953 + t467;
t1157 = t875 * t193;
t231 = t379 * t877 + t1390;
t232 = t382 * t875 + t1451;
t82 = t467 + (-t1131 * t875 + t1389) * t876;
t83 = -t1425 + (t1128 * t877 + t1449) * t876;
t1356 = -t1157 / 0.4e1 + t1144 / 0.4e1 + t231 * t1003 + t82 * t1313 + t83 * t1310 + t100 * t1001 - t1402 + (t232 + t101) * t1007;
t1012 = t1161 / 0.4e1;
t694 = t733 * t875;
t696 = t735 * t875;
t315 = (-t694 * t856 + t696 * t857 + t605) * t876 + t908 * t874;
t695 = t733 * t877;
t697 = t735 * t877;
t316 = (t695 * t856 - t697 * t857 + t607) * t876 + t907 * t874;
t732 = Icges(6,6) * t876 - t874 * t963;
t734 = Icges(6,5) * t876 - t874 * t967;
t361 = -t1372 * t875 + t732 * t772 + t734 * t773;
t362 = t1372 * t877 + t774 * t732 - t775 * t734;
t1172 = t857 * t734;
t1175 = t856 * t732;
t408 = (t1172 + t731 - t1175) * t876 + t903 * t874;
t416 = t876 * t935 + t1193;
t515 = t876 * t923 + t1184;
t1355 = t515 * t1312 + t408 * t1316 - t1404 + (t416 + t469) * t1012 + (-t417 + t471) * t1010 + (t315 + t361) * t1007 + (t316 + t362) * t1001;
t1348 = -0.2e1 * t483;
t519 = t1149 * t723 + t597 * t874;
t1347 = 0.2e1 * t519;
t1346 = 0.2e1 * t1438;
t1345 = 0.2e1 * t538;
t1344 = -0.2e1 * t628;
t1342 = 0.4e1 * t815;
t378 = (t1388 * t875 - t683 * t877) * t876 + t1105;
t429 = -t1161 * t723 + t1020;
t430 = t598 * t876 + t699;
t465 = (-t597 * t877 + t598 * t875) * t876;
t1332 = (t1438 * t256 + t207 * t465 + t255 * t519 + t319 * t378 + t372 * t429 + t373 * t430) * t1051;
t1330 = -0.2e1 * t519 * t1446;
t805 = pkin(3) * t1161 - t840;
t1087 = -pkin(4) * t1169 - t840 - (t1168 - t1381) * t875 - t805;
t984 = pkin(3) * t1160 - t842;
t791 = t877 * t984;
t1088 = t877 * (t984 + t1019) - t791;
t242 = t1108 * t877 + (t1087 - t1109) * t875 + t1088;
t198 = 0.2e1 * t489 * t242;
t396 = t537 * t1348;
t397 = t481 * t1345;
t1025 = -t198 + t396 + t397;
t1048 = 0.2e1 * t756;
t495 = -t875 * t627 + t604;
t1327 = m(7) * (0.2e1 * t495 * t319 + (-t372 * t877 + t373 * t875) * t1048 + t1025);
t1325 = m(7) * (-t1346 * t664 + t1347 * t665 + 0.2e1 * t434 * t465 + t1025);
t1324 = -t193 / 0.2e1;
t1323 = t325 / 0.2e1;
t1322 = t646 / 0.2e1;
t766 = Icges(5,6) * t876 - t874 * t964;
t1321 = t766 / 0.2e1;
t1320 = -t871 / 0.2e1;
t1317 = t872 / 0.2e1;
t399 = 0.2e1 * t537 * t488;
t400 = t487 * t1345;
t1124 = t399 + t400;
t1283 = m(7) * (t1344 * t373 + 0.2e1 * t372 * t627 + t1124);
t1116 = t519 * t1044 + t1045 * t1438;
t1206 = t430 * t875;
t1207 = t429 * t877;
t949 = t1206 - t1207;
t1282 = m(7) * (t378 * t1341 + (t465 - t949) * t1047 + t1116);
t1280 = (t1438 * t575 + t429 * t488 + t430 * t487 + t519 * t574) * t1051;
t1278 = m(7) * (-t1090 * t1347 + t1346 * t1417 + t1124);
t1275 = m(7) * (t481 * t1344 + t627 * t1348 + (t487 * t875 - t488 * t877) * t1048);
t1200 = t1438 * t877;
t1201 = t519 * t875;
t393 = t1042 * t1438 + t519 * t1043;
t1274 = m(7) * ((t1200 + t1201) * t1047 + t393);
t405 = t1340 * t1438 + 0.2e1 * t1201;
t1270 = m(7) * (-0.2e1 * t1200 + t405 - 0.2e1 * t1201);
t1268 = m(7) * (t1343 * t465 + t1116);
t1264 = m(7) * t393;
t1263 = m(7) * t405;
t1262 = m(7) * (t1341 * t495 - t756 * t970);
t1258 = t931 * t1393;
t1257 = t931 * t1051;
t1252 = m(5) * qJD(3);
t1250 = m(6) * qJD(3);
t1249 = m(6) * qJD(5);
t1247 = m(7) * qJD(3);
t1246 = m(7) * qJD(5);
t1245 = m(7) * qJD(6);
t1176 = t854 * t876;
t1162 = t874 * t782;
t1086 = -t1362 + t722;
t792 = t877 * t1065;
t1084 = t877 * (t1160 * t873 - t1065 - t818) + t792;
t1083 = t1161 * t873 - t816;
t728 = t1142 * t876 + t1381;
t828 = qJ(4) * t876 - t1256;
t811 = t875 * t828;
t1082 = t875 * t728 + t811;
t1078 = -t728 - t828;
t156 = 0.4e1 * t1378 + 0.4e1 * t1379 + 0.4e1 * t1380 + (t690 * t877 + t691 * t875) * m(4) + m(3) * ((rSges(3,3) * t877 + t860) * t877 + (rSges(3,3) + qJ(2)) * t868);
t1062 = qJD(1) * t156;
t186 = (-t1378 - t1379 - t1380) * t1046;
t1061 = qJD(1) * t186;
t1060 = qJD(1) * t876;
t1058 = qJD(3) * t874;
t1057 = qJD(3) * t876;
t1056 = qJD(5) * t876;
t998 = t1397 + t870 / 0.2e1;
t103 = ((-t665 / 0.2e1 + t255 / 0.2e1) * t877 + (-t664 / 0.2e1 - t256 / 0.2e1) * t875) * m(7) + (t1205 / 0.2e1 - t1204 / 0.2e1 + t998 * t785) * m(6);
t1055 = t103 * qJD(2);
t252 = (t1207 / 0.2e1 - t1206 / 0.2e1 + t998 * t756) * m(7);
t1054 = t252 * qJD(2);
t1031 = t1246 / 0.2e1;
t1030 = t1245 / 0.2e1;
t1029 = t1324 + t82 / 0.2e1;
t1028 = -t83 / 0.2e1 - t194 / 0.2e1;
t1002 = t1147 / 0.2e1;
t1008 = -t1149 / 0.2e1;
t114 = t368 * t874 + (-t262 * t875 + t263 * t877) * t876;
t115 = t369 * t874 + (-t264 * t875 + t265 * t877) * t876;
t427 = (t1163 + (-t1079 * t853 - t1080 * t854) * t876) * t874;
t1027 = t115 * t1002 + t114 * t1008 + (t427 + (-t307 * t875 + t308 * t877) * t876) * t1316;
t1016 = -t1176 / 0.2e1;
t1015 = t1176 / 0.2e1;
t1013 = t1161 / 0.2e1;
t1011 = -t1160 / 0.2e1;
t1006 = t1149 / 0.2e1;
t1005 = t1149 / 0.4e1;
t1004 = -t1147 / 0.2e1;
t985 = 0.2e1 * t1063;
t977 = t67 * t1002 + t187 * t1004 + t1428 * t1008;
t971 = t755 * t1015 + t719 * t1016 + t1401;
t961 = Icges(4,5) * t876 - Icges(4,6) * t874;
t950 = -t416 * t875 - t417 * t877;
t947 = t481 * t875 + t483 * t877;
t944 = t537 * t877 - t538 * t875;
t941 = t571 * t875 + t573 * t877;
t925 = t673 * t875 + t675 * t877;
t919 = t985 / 0.4e1;
t918 = t874 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1);
t267 = -t1367 * t875 + t679 * t749 + t681 * t750;
t268 = -t1368 * t875 - t680 * t749 - t682 * t750;
t55 = (-t267 * t875 + t268 * t877 + t457) * t876 + (t334 - t955) * t874;
t269 = t1367 * t877 + t751 * t679 - t752 * t681;
t270 = t1368 * t877 - t751 * t680 + t752 * t682;
t56 = (-t269 * t875 + t270 * t877 + t459) * t876 + (t335 - t954) * t874;
t917 = t56 * t1002 + t55 * t1008 + t188 * t1011 + t187 * t1013 + t1374;
t916 = t1330 / 0.4e1 + t977;
t289 = t1095 * t772 - t1097 * t773 - t1149 * t635;
t290 = t1094 * t772 - t1096 * t773 - t1149 * t636;
t164 = t289 * t877 + t290 * t875;
t291 = t1095 * t774 + t1097 * t775 + t1147 * t635;
t292 = t1094 * t774 + t1096 * t775 + t1147 * t636;
t165 = t291 * t877 + t292 * t875;
t902 = t1311 * t164 + t1314 * t165;
t899 = t1332 / 0.4e1 + t917;
t897 = 0.4e1 * t947;
t896 = 0.4e1 * t941;
t895 = t1438 * t538 - t465 * t489 + t519 * t537;
t219 = t357 * t877 + t1392;
t88 = t1137 * t877 + t1392;
t893 = t88 * t1001 + t219 * t1003 + t1007 * t1427 + t1357;
t892 = t719 * t1015 + t755 * t1016 - t1401;
t891 = t67 * t1004 + t427 + (t307 + t368) * t1008 + t1428 * t1006 + (t187 + t308 + t369) * t1002;
t152 = t267 * t877 + t268 * t875;
t153 = t269 * t877 + t270 * t875;
t890 = t55 * t1311 + t56 * t1314 + t152 * t1008 + t153 * t1002 + (t294 * t877 + t295 * t875) * t1316 + t219 * t1013 + t220 * t1011 + (t401 * t877 - t402 * t875) * t1312 - t1138;
t889 = t1001 * t1422 + t1007 * t1423 + t1010 * t1420 + t1012 * t1421 + t1373;
t888 = t151 * t1002 + t56 * t1004 + t55 * t1006 + t150 * t1008 + t114 * t1311 + t115 * t1314 + t1158 * t1396 - t1374 + (t1214 + t1215 + t1145) * t1316;
t552 = -t875 * t777 + t1387;
t884 = -t219 / 0.2e1 + t88 / 0.2e1 - t231 / 0.2e1 + t100 / 0.2e1 - (-t1149 * t644 + t1445) * t877 / 0.2e1 + t1419 * t1315 + (-t876 * t929 - t1439 + t550) * t1314 + (t552 - t1387 + (t777 + t922) * t875 + t1445) * t1311;
t883 = t89 / 0.2e1 + t220 / 0.2e1 + t101 / 0.2e1 + t232 / 0.2e1 + t777 * t1397 + (t1147 * t646 + t552) * t1314 + ((t777 - t922) * t877 + t1419 + t1439) * t1311;
t882 = t219 * t1001 + t88 * t1003 + t1005 * t1427 - t1357 + t889 + t987;
t881 = -t1373 + t893 + t987 - t1421 * t1161 / 0.4e1 + t1420 * t1160 / 0.4e1 + t1423 * t1005 + t1422 * t1003;
t880 = t889 + t893 - t987;
t768 = Icges(5,5) * t876 - t874 * t968;
t879 = t766 * t1319 + t768 * t1318 - t1172 / 0.2e1 + t1175 / 0.2e1 - t1178 / 0.2e1 + t1181 / 0.2e1 + t824 / 0.2e1 + t969 / 0.2e1 - t765 / 0.2e1 - t731 / 0.2e1 - t717 / 0.2e1;
t800 = t961 * t877;
t799 = t875 * t961;
t770 = t1407 * t874 + t1240;
t711 = t769 * t877;
t710 = t769 * t875;
t709 = t767 * t877;
t708 = t767 * t875;
t674 = (-t770 - t828) * t877;
t672 = t770 * t875 + t811;
t633 = t918 * t985;
t572 = (-t744 + t1078) * t877;
t570 = t1082 + t1183;
t500 = -t657 * t875 - t771 * t870 + t792;
t498 = t930 * t876;
t491 = -t1257 / 0.4e1;
t484 = t1258 / 0.4e1;
t482 = (t1078 - t1086) * t877;
t480 = t1086 * t875 + t1082;
t462 = t877 * (rSges(5,3) * t1147 - t1383) - t791 + (-t805 - t915) * t875;
t445 = (t1162 + (-t1075 * t856 - t1076 * t857) * t876) * t874;
t426 = t1435 * t876;
t425 = t1262 / 0.4e1;
t424 = 0.2e1 * t944;
t419 = t424 * t1030;
t410 = -t1424 * t877 + (-t701 + t1083) * t875 + t1084;
t398 = t1263 / 0.4e1;
t385 = t1264 / 0.4e1;
t349 = 0.4e1 * t566 * t658 + 0.4e1 * t567 * t657;
t346 = t1342 * t462 + 0.4e1 * t874 * t925;
t341 = -t615 * t877 + (-t614 + t1087) * t875 + t1088;
t323 = -0.4e1 * t517 * t643 + 0.4e1 * t518 * t642;
t318 = t1047 * t944 - t1341 * t489;
t317 = 0.4e1 * t517 * t602 + 0.4e1 * t518 * t601;
t314 = t318 * t1030;
t311 = t1091 * t877 + (t1083 + t1092) * t875 + t1084;
t293 = -0.4e1 * t487 * t628 + 0.4e1 * t488 * t627;
t284 = t1370 * t877 - t774 * t695 + t775 * t697;
t283 = t1369 * t877 + t774 * t694 - t775 * t696;
t282 = -t1370 * t875 - t695 * t772 - t697 * t773;
t281 = -t1369 * t875 + t694 * t772 + t696 * t773;
t278 = t1268 / 0.4e1;
t273 = 0.4e1 * t487 * t575 + 0.4e1 * t488 * t574;
t259 = -0.4e1 * t1090 * t488 + 0.4e1 * t1417 * t487;
t253 = t1270 / 0.4e1;
t251 = m(7) * t756 * t919 + t949 * t992;
t238 = t1342 * t341 + t874 * t896;
t211 = t515 * t874 + t876 * t950;
t202 = t1335 * t293 + t971;
t197 = 0.4e1 * t341 * t513 + t785 * t896;
t195 = t1274 / 0.4e1;
t189 = 0.4e1 * t895;
t172 = t1275 / 0.4e1;
t166 = t1342 * t242 + t874 * t897;
t163 = t283 * t877 + t284 * t875;
t162 = t281 * t877 + t282 * t875;
t159 = t398 + t253 + t1257 / 0.4e1;
t158 = t491 + t398 - t1270 / 0.4e1;
t157 = t491 + t253 - t1263 / 0.4e1;
t154 = t1278 / 0.4e1;
t135 = 0.4e1 * t1437 * t441 - 0.4e1 * t407 * t478 + 0.4e1 * t440 * t530;
t134 = 0.4e1 * t242 * t495 + t756 * t897;
t130 = t385 + t195 - t1258 / 0.4e1;
t129 = t484 + t385 - t1274 / 0.4e1;
t128 = t484 + t195 - t1264 / 0.4e1;
t127 = 0.4e1 * t1438 * t430 + 0.4e1 * t378 * t465 + 0.4e1 * t429 * t519;
t126 = 0.4e1 * t242 * t434 - 0.4e1 * t481 * t664 - 0.4e1 * t483 * t665;
t124 = t395 * t874 + (-t291 * t875 + t292 * t877) * t876;
t123 = t394 * t874 + (-t289 * t875 + t290 * t877) * t876;
t122 = t1280 / 0.4e1;
t118 = t1282 / 0.4e1;
t116 = t1283 / 0.4e1;
t111 = t259 * t1335 + t1162 / 0.2e1 + t323 * t1336 + t1363 * t876 + t971;
t108 = 0.4e1 * t978;
t102 = t948 * t1441 + m(6) * t785 * t919 + (t926 + t956) * t992;
t97 = (-t315 * t875 + t316 * t877 + t515) * t876 + (t408 - t950) * t874;
t95 = t1335 * t166 + t1336 * t238 + t1337 * t346;
t94 = t1021 + t1022 - t986;
t93 = t986 + t1412;
t92 = t986 - t1412;
t76 = t278 + t118 - t1262 / 0.4e1;
t75 = t425 + t278 - t1282 / 0.4e1;
t74 = t425 + t118 - t1268 / 0.4e1;
t73 = 0.4e1 * t1125 * t372;
t72 = (-t283 * t875 + t284 * t877 + t471) * t876 + (t362 - t952) * t874;
t71 = (-t281 * t875 + t282 * t877 + t469) * t876 + (t361 - t953) * t874;
t69 = t1325 / 0.4e1;
t48 = t273 * t1335 + t317 * t1336 + t349 * t1337 - t879 * t876 - t1366;
t46 = t1327 / 0.4e1;
t43 = 0.4e1 * t207 * t319 + 0.4e1 * t255 * t372 + 0.4e1 * t256 * t373;
t38 = t1132 + t1140 - t1120;
t37 = t1120 - t1410;
t36 = t1120 + t1410;
t35 = t1335 * t134 + t1138;
t34 = t1133 + t1141 - t1121;
t33 = t1121 - t1409;
t32 = t1121 + t1409;
t27 = t1024 + t1026 - t1023;
t26 = t1023 + t1411;
t25 = t1023 - t1411;
t24 = t1335 * t189 + t1027;
t23 = t24 * qJD(6);
t22 = t108 * t1335 + t1027;
t21 = t126 * t1335 + t1336 * t197 + t1138 + t902;
t20 = t1139 + t1233 - t1127;
t19 = t1127 - t1408;
t18 = t1127 + t1408;
t16 = t127 * t1335 + t917;
t15 = t154 - t1283 / 0.4e1 + t916;
t14 = t116 - t1278 / 0.4e1 + t916;
t13 = t116 + t154 - t1330 / 0.4e1 + t891;
t12 = -t73 * t1335 + (t1028 * t875 + t1029 * t877) * t876 + t977;
t11 = t69 - t1327 / 0.4e1 + t899;
t10 = t46 - t1325 / 0.4e1 + t899;
t9 = t875 * t884 + t877 * t883;
t8 = t43 * t1335 + t135 * t1336 + (t211 / 0.2e1 + t71 * t1315 + t72 * t1311) * t876 + (t97 / 0.2e1 + t1157 / 0.2e1 - t1144 / 0.2e1) * t874 + t917;
t7 = t882 + t172 + t122;
t6 = t880 + t122 - t1275 / 0.4e1;
t5 = t881 + t172 - t1280 / 0.4e1;
t4 = t46 + t888 - t1332 / 0.4e1 + t69;
t3 = t882 + (-t83 / 0.4e1 - t194 / 0.4e1 + (-t100 / 0.4e1 + t231 / 0.4e1) * t876) * t877 + (-t82 / 0.4e1 + t193 / 0.4e1 + (t101 / 0.4e1 + t232 / 0.4e1) * t876) * t875 + t1355 + t1359 + t1402;
t2 = t881 + (-t515 / 0.2e1 + (-t362 / 0.4e1 - t316 / 0.4e1) * t877 + (t315 / 0.4e1 + t361 / 0.4e1) * t875) * t876 + (-t408 / 0.2e1 + (t471 / 0.4e1 - t417 / 0.4e1) * t877 + (-t416 / 0.4e1 - t469 / 0.4e1) * t875) * t874 + t1356 + t1359 + t1404;
t1 = t880 + (-t395 / 0.4e1 - t325 / 0.4e1) * t875 + (-t394 / 0.4e1 - t324 / 0.4e1) * t877 + t1356 + t1355 - t1405;
t17 = [t156 * qJD(2) + t48 * qJD(3) + t186 * qJD(4) + t111 * qJD(5) + t202 * qJD(6), qJD(3) * t92 + qJD(5) * t36 + qJD(6) * t158 + t1062, t48 * qJD(1) + t92 * qJD(2) + t26 * qJD(4) + t3 * qJD(5) + t7 * qJD(6) + ((t795 * t1321 + t768 * t796 / 0.2e1 + t361 / 0.2e1 + t334 / 0.2e1 + t315 / 0.2e1 + t294 / 0.2e1 - t960 * t1311 - t883) * qJD(3) + (t708 * t1320 + t710 * t1317 + t644 / 0.2e1 - t778 / 0.2e1 + t803 / 0.2e1) * t1057 + (t1190 / 0.2e1 - t1188 / 0.2e1 + t706 / 0.2e1 - t780 / 0.2e1 - t801 / 0.2e1) * t1058) * t877 + ((t335 / 0.2e1 + t295 / 0.2e1 + t362 / 0.2e1 + t316 / 0.2e1 + t797 * t1321 - t798 * t768 / 0.2e1 - t960 * t1314 - t884) * qJD(3) + (t779 / 0.2e1 - t804 / 0.2e1 + t1322 - t709 * t1320 - t711 * t1317) * t1057 + (t781 / 0.2e1 + t802 / 0.2e1 - t707 / 0.2e1 + t1189 / 0.2e1 - t1187 / 0.2e1) * t1058) * t875 + (t480 * t487 + t481 * t575 + t482 * t488 - t483 * t574) * t1247 + (t517 * t570 + t518 * t572 + t571 * t602 - t573 * t601) * t1250 + (t566 * t672 + t567 * t674 - t657 * t675 + t658 * t673) * t1252 + (t920 * t832 - (t690 * t875 - t691 * t877) * t976) * qJD(3) * m(4), qJD(3) * t26 + qJD(5) * t32 + qJD(6) * t129 + t1061, t111 * qJD(1) + t36 * qJD(2) + t3 * qJD(3) + t32 * qJD(4) + (t445 + t891) * qJD(5) + t13 * qJD(6) + (t73 / 0.4e1 - t372 * t1090 + t373 * t1417 + t475 * t488 + t476 * t487) * t1246 + (-t1437 * t643 + t517 * t554 + t518 * t553 + t530 * t642) * t1249 + ((t395 / 0.2e1 + t1323 - t1029) * t877 + (-t394 / 0.2e1 - t324 / 0.2e1 - t1028) * t875) * t1056, t202 * qJD(1) + t158 * qJD(2) + t7 * qJD(3) + t129 * qJD(4) + t13 * qJD(5) + t891 * qJD(6) + (t519 * t627 - t1438 * t628 + t399 / 0.2e1 + t400 / 0.2e1) * t1245; t93 * qJD(3) + t37 * qJD(5) + t157 * qJD(6) - t1062, 0, t93 * qJD(1) + (-m(4) * t976 * t985 / 0.2e1 + (-m(5) * t674 - m(6) * t572 - m(7) * t482) * t877 + (m(5) * t672 + m(6) * t570 + m(7) * t480) * t875) * qJD(3) + t633 * qJD(4) + t102 * qJD(5) + t251 * qJD(6), t633 * qJD(3), t37 * qJD(1) + t102 * qJD(3) + qJD(5) * t1364 - t419, t157 * qJD(1) + t251 * qJD(3) - t1031 * t424 - t419; t94 * qJD(2) + t9 * qJD(3) + t27 * qJD(4) + t2 * qJD(5) + t5 * qJD(6) + t879 * t1060 + t317 * t1447 + t273 * t1448 + ((-t646 / 0.2e1 + t1322) * t1160 - m(5) * t349 / 0.4e1 + t1366) * qJD(1), qJD(1) * t94 + qJD(5) * t103 + qJD(6) * t252, t9 * qJD(1) + t95 * qJD(4) + t21 * qJD(5) + t35 * qJD(6) + ((t242 * t311 + t480 * t481 - t482 * t483) * m(7) + (t341 * t410 + t570 * t571 - t572 * t573) * m(6) + (t462 * t500 + t672 * t673 - t674 * t675) * m(5) + ((-t877 * t894 + (-t877 * rSges(4,3) - t875 * t976) * t875) * (-t875 * t807 - t809 * t877) - t1063 * t832 * t976) * m(4) + (t153 + t163 - t868 * t800 + (-t797 * t709 + t798 * t711) * t875 + (t797 * t708 - t798 * t710 + t875 * t799 + t1436) * t877) * t1314 + (t152 + t162 + (t708 * t795 + t710 * t796) * t877 + t870 * t799 + (-t795 * t709 - t796 * t711 - t877 * t800 - t1436) * t875) * t1311) * qJD(3), t27 * qJD(1) + t95 * qJD(3) + t18 * qJD(5) + t75 * qJD(6) + (-t545 + (t1343 - t970) * t918) * qJD(4), t2 * qJD(1) + t1055 + t21 * qJD(3) + t18 * qJD(4) + t4 * qJD(6) + (-t211 / 0.2e1 + (-t72 / 0.2e1 + t165 / 0.2e1) * t877 + (t71 / 0.2e1 - t164 / 0.2e1) * t875) * t1056 + (-t43 / 0.4e1 + t242 * t426 + t319 * t434 + t372 * t665 - t373 * t664 - t475 * t483 + t476 * t481) * t1246 + (-t135 / 0.4e1 - t498 * t341 - t478 * t513 - t553 * t573 + t554 * t571 + (t1437 * t875 - t530 * t877) * t785) * t1249 + (t123 * t1311 + t124 * t1314 + t888 + (-t97 / 0.2e1 + (t194 / 0.2e1 + t324 / 0.2e1) * t877 + (t1324 + t1323) * t875) * t874) * qJD(5), t5 * qJD(1) + t1054 + t35 * qJD(3) + t75 * qJD(4) + t4 * qJD(5) + t888 * qJD(6) + (t465 * t495 - t198 / 0.2e1 + t396 / 0.2e1 + t397 / 0.2e1 - t127 / 0.4e1 + (t1438 * t875 - t519 * t877) * t756) * t1245; t25 * qJD(3) + t33 * qJD(5) + t128 * qJD(6) - t1061, 0, t25 * qJD(1) + t545 * qJD(4) + t19 * qJD(5) + t74 * qJD(6) + (-t166 / 0.4e1 + (-t480 * t875 + t482 * t877 + t242) * t876 + (t311 + t947) * t874) * t1247 + (-t238 / 0.4e1 + (-t570 * t875 + t572 * t877 + t341) * t876 + (t410 + t941) * t874) * t1250 + (-t346 / 0.4e1 + (-t672 * t875 + t674 * t877 + t462) * t876 + (t500 + t925) * t874) * t1252, t545 * qJD(3), t33 * qJD(1) + t19 * qJD(3) + ((-m(6) * t498 + m(7) * t426) * t874 - t1364 * t876) * qJD(5) + t314, t128 * qJD(1) + t74 * qJD(3) + t1031 * t318 + t314; (-t1162 / 0.2e1 + t892) * qJD(1) + t38 * qJD(2) + t1 * qJD(3) + t34 * qJD(4) + t12 * qJD(5) + t14 * qJD(6) + (-t259 / 0.4e1 - t1125 * t488) * t1248 + t323 * t1447 - t1363 * t1060, qJD(1) * t38 - qJD(3) * t103, t1 * qJD(1) - t1055 + (t890 + t163 * t1002 + t232 * t1011 + t162 * t1008 + t231 * t1013 + t71 * t1311 + (t416 * t877 - t417 * t875) * t1312 + (t315 * t877 + t316 * t875) * t1316 + t72 * t1314 - t902) * qJD(3) + t20 * qJD(4) + t8 * qJD(5) + t10 * qJD(6) + (-t126 / 0.4e1 + t207 * t242 - t255 * t483 + t256 * t481 + t311 * t319 + t372 * t482 + t373 * t480) * t1247 + (-t197 / 0.4e1 + t341 * t407 - t410 * t478 - t440 * t573 + t441 * t571 + t530 * t572 + t1437 * t570) * t1250, qJD(1) * t34 + qJD(3) * t20, t12 * qJD(1) + t8 * qJD(3) + ((t319 * t426 + t372 * t475 + t373 * t476) * m(7) + (t1437 * t554 + t478 * t498 + t530 * t553) * m(6) + t445 * t1316 + (t123 * t1315 + t124 * t1311 + (-t324 * t875 + t325 * t877) * t1316) * t876 + t1027) * qJD(5) + t22 * qJD(6), t14 * qJD(1) + t10 * qJD(3) + t22 * qJD(5) + t1027 * qJD(6) + (-t189 / 0.4e1 + t895 + t978) * t1245; t892 * qJD(1) + t159 * qJD(2) + t6 * qJD(3) + t130 * qJD(4) + t15 * qJD(5) + qJD(6) * t977 + t1448 * t293, qJD(1) * t159 - qJD(3) * t252, t6 * qJD(1) - t1054 + t890 * qJD(3) + t76 * qJD(4) + t11 * qJD(5) + t16 * qJD(6) + (t242 * t378 + t311 * t465 - t429 * t483 + t430 * t481 + t480 * t1438 + t482 * t519 - t134 / 0.4e1) * t1247, qJD(1) * t130 + qJD(3) * t76, t15 * qJD(1) + t11 * qJD(3) + t1027 * qJD(5) + t23 + (t426 * t465 + t475 * t519 + t476 * t1438 - t108 / 0.4e1 + t978) * t1246, qJD(1) * t977 + qJD(3) * t16 + qJD(5) * t24 + t23;];
Cq  = t17;
