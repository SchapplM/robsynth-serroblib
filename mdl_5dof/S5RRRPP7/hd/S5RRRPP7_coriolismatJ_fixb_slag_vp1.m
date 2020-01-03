% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP7
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:04:00
% EndTime: 2019-12-31 21:05:29
% DurationCPUTime: 78.94s
% Computational Cost: add. (56857->1346), mult. (142738->1790), div. (0->0), fcn. (158003->6), ass. (0->751)
t787 = sin(qJ(2));
t788 = sin(qJ(1));
t1060 = t787 * t788;
t786 = sin(qJ(3));
t791 = cos(qJ(1));
t1048 = t791 * t786;
t789 = cos(qJ(3));
t790 = cos(qJ(2));
t1053 = t789 * t790;
t712 = t1053 * t788 - t1048;
t690 = Icges(5,5) * t712;
t1052 = t789 * t791;
t1054 = t788 * t790;
t711 = t1054 * t786 + t1052;
t513 = -Icges(5,6) * t1060 - Icges(5,3) * t711 - t690;
t515 = Icges(4,5) * t712 - Icges(4,6) * t711 + Icges(4,3) * t1060;
t521 = Icges(5,4) * t712 + Icges(5,2) * t1060 + Icges(5,6) * t711;
t696 = Icges(4,4) * t712;
t524 = -Icges(4,2) * t711 + Icges(4,6) * t1060 + t696;
t689 = Icges(5,5) * t711;
t530 = Icges(5,1) * t712 + Icges(5,4) * t1060 + t689;
t695 = Icges(4,4) * t711;
t534 = -Icges(4,1) * t712 - Icges(4,5) * t1060 + t695;
t1300 = (t530 - t534) * t712 + (-t513 - t524) * t711 + (t515 + t521) * t1060;
t1057 = t787 * t791;
t1051 = t790 * t791;
t1056 = t788 * t786;
t714 = t1051 * t789 + t1056;
t691 = Icges(5,5) * t714;
t1055 = t788 * t789;
t713 = t1048 * t790 - t1055;
t514 = Icges(5,6) * t1057 + Icges(5,3) * t713 + t691;
t517 = Icges(4,5) * t714 - Icges(4,6) * t713 + Icges(4,3) * t1057;
t523 = Icges(5,4) * t714 + Icges(5,2) * t1057 + Icges(5,6) * t713;
t1096 = Icges(4,4) * t714;
t526 = -Icges(4,2) * t713 + Icges(4,6) * t1057 + t1096;
t1089 = Icges(5,5) * t713;
t532 = Icges(5,1) * t714 + Icges(5,4) * t1057 + t1089;
t697 = Icges(4,4) * t713;
t535 = Icges(4,1) * t714 + Icges(4,5) * t1057 - t697;
t1299 = (t532 + t535) * t712 + (t514 - t526) * t711 + (t517 + t523) * t1060;
t1061 = t786 * t787;
t1083 = Icges(5,6) * t790;
t1059 = t787 * t789;
t764 = Icges(5,5) * t1059;
t648 = Icges(5,3) * t1061 - t1083 + t764;
t887 = Icges(5,4) * t789 + Icges(5,6) * t786;
t656 = -Icges(5,2) * t790 + t787 * t887;
t1088 = Icges(5,5) * t786;
t890 = Icges(5,1) * t789 + t1088;
t664 = -Icges(5,4) * t790 + t787 * t890;
t411 = t1060 * t656 + t648 * t711 + t664 * t712;
t884 = Icges(4,5) * t789 - Icges(4,6) * t786;
t650 = -Icges(4,3) * t790 + t787 * t884;
t1094 = Icges(4,4) * t789;
t888 = -Icges(4,2) * t786 + t1094;
t658 = -Icges(4,6) * t790 + t787 * t888;
t1095 = Icges(4,4) * t786;
t891 = Icges(4,1) * t789 - t1095;
t666 = -Icges(4,5) * t790 + t787 * t891;
t412 = t1060 * t650 - t658 * t711 + t666 * t712;
t1298 = t411 + t412;
t1297 = t1299 * t791 + t1300 * t788;
t1173 = t787 / 0.2e1;
t1172 = t788 / 0.2e1;
t1170 = -t790 / 0.2e1;
t1233 = -t791 / 0.2e1;
t1194 = m(6) / 0.2e1;
t1196 = m(5) / 0.2e1;
t754 = pkin(2) * t787 - pkin(7) * t790;
t1100 = -rSges(6,3) - qJ(5);
t1292 = t1100 * t790;
t1105 = rSges(6,2) * t786;
t894 = rSges(6,1) * t789 + t1105;
t1257 = -pkin(4) * t1059 - t787 * t894 + t1292;
t893 = pkin(3) * t789 + qJ(4) * t786;
t734 = t893 * t787;
t943 = t734 - t1257;
t920 = t754 + t943;
t450 = t920 * t788;
t452 = t920 * t791;
t895 = rSges(5,1) * t789 + rSges(5,3) * t786;
t1223 = t895 * t787;
t671 = -rSges(5,2) * t790 + t1223;
t976 = t671 + t734;
t942 = t754 + t976;
t488 = t942 * t788;
t490 = t942 * t791;
t1164 = rSges(5,1) + pkin(3);
t686 = t711 * qJ(4);
t699 = t711 * rSges(5,3);
t1162 = rSges(5,2) + pkin(7);
t1118 = pkin(2) * t790;
t938 = pkin(1) + t1118;
t1215 = t1162 * t787 + t938;
t782 = t791 * pkin(6);
t818 = -t1215 * t788 + t782;
t407 = -t1164 * t712 - t686 - t699 + t818;
t1101 = t713 * rSges(5,3);
t1116 = t788 * pkin(6);
t688 = t713 * qJ(4);
t931 = t688 + t1116;
t408 = t1164 * t714 + t1215 * t791 + t1101 + t931;
t830 = (-t407 * t791 - t408 * t788) * t1061;
t700 = t711 * rSges(6,2);
t1219 = t1100 * t1060;
t755 = pkin(7) * t787 + t1118;
t854 = pkin(1) + t755;
t807 = -t788 * t854 - t1219 + t782;
t1163 = rSges(6,1) + pkin(4);
t952 = pkin(3) + t1163;
t390 = -t712 * t952 - t686 - t700 + t807;
t1218 = t713 * rSges(6,2) + t1057 * t1100;
t391 = t714 * t952 + t791 * t854 + t1218 + t931;
t870 = -t390 * t791 - t391 * t788;
t831 = t870 * t1061;
t1043 = (-t713 * t450 + t711 * t452 + t831) * t1194 + (-t713 * t488 + t711 * t490 + t830) * t1196;
t774 = pkin(7) * t1051;
t434 = t774 + (t1292 + (-pkin(2) + (-rSges(6,2) - qJ(4)) * t786 - t952 * t789) * t787) * t791;
t644 = (-pkin(3) * t1055 - qJ(4) * t1056) * t787;
t773 = pkin(2) * t1060;
t939 = t773 - t644;
t435 = ((-pkin(7) - t1100) * t790 + (t1163 * t789 + t1105) * t787) * t788 + t939;
t772 = rSges(5,2) * t1051;
t463 = t772 + t774 + (-pkin(2) - t1164 * t789 + (-rSges(5,3) - qJ(4)) * t786) * t1057;
t464 = (-t1162 * t790 + t1223) * t788 + t939;
t1044 = (t711 * t434 + t713 * t435 + t831) * t1194 + (t711 * t463 + t713 * t464 + t830) * t1196;
t16 = t1044 - t1043;
t1296 = t16 * qJD(1);
t1217 = -t1163 * t712 - t700;
t1212 = -t1217 + t1219;
t589 = -t712 * pkin(3) - t686;
t594 = t714 * pkin(3) + t688;
t784 = t788 ^ 2;
t785 = t791 ^ 2;
t960 = t784 + t785;
t970 = t960 * t755;
t923 = -t589 * t788 + t791 * t594 + t970;
t998 = t1163 * t714 + t1218;
t252 = t1212 * t788 + t791 * t998 + t923;
t967 = -t712 * rSges(5,1) - t699;
t537 = rSges(5,2) * t1060 - t967;
t543 = t714 * rSges(5,1) + rSges(5,2) * t1057 + t1101;
t289 = t788 * t537 + t543 * t791 + t923;
t559 = t788 * t711 + t713 * t791;
t1106 = rSges(6,2) * t712;
t585 = -rSges(6,1) * t711 + t1106;
t591 = -t713 * rSges(6,1) + t714 * rSges(6,2);
t687 = t712 * qJ(4);
t584 = -pkin(3) * t711 + t687;
t590 = -t713 * pkin(3) + qJ(4) * t714;
t995 = t788 * t584 + t791 * t590;
t292 = -pkin(4) * t559 + t788 * t585 + t591 * t791 + t995;
t1103 = rSges(5,3) * t712;
t586 = -rSges(5,1) * t711 + t1103;
t592 = -t713 * rSges(5,1) + t714 * rSges(5,3);
t349 = t788 * t586 + t592 * t791 + t995;
t1220 = -t712 * rSges(4,1) + t711 * rSges(4,2);
t541 = -rSges(4,3) * t1060 + t1220;
t897 = t714 * rSges(4,1) - t713 * rSges(4,2);
t544 = rSges(4,3) * t1057 + t897;
t387 = -t541 * t788 + t544 * t791 + t970;
t587 = -rSges(4,1) * t711 - rSges(4,2) * t712;
t593 = -rSges(4,1) * t713 - rSges(4,2) * t714;
t433 = t788 * t587 + t593 * t791;
t730 = (-pkin(3) * t786 + qJ(4) * t789) * t787;
t731 = (-rSges(6,1) * t786 + rSges(6,2) * t789) * t787;
t966 = -t730 - t731;
t853 = pkin(4) * t1061 + t966;
t501 = t853 * t788;
t502 = t853 * t791;
t732 = (-rSges(5,1) * t786 + rSges(5,3) * t789) * t787;
t965 = -t730 - t732;
t560 = t965 * t788;
t561 = t965 * t791;
t1108 = rSges(4,1) * t789;
t896 = -rSges(4,2) * t786 + t1108;
t835 = t896 * t787;
t672 = -rSges(4,3) * t790 + t835;
t975 = t672 + t754;
t562 = t975 * t788;
t564 = t975 * t791;
t733 = (-rSges(4,1) * t786 - rSges(4,2) * t789) * t787;
t1295 = -m(6) * (t252 * t292 - t450 * t501 - t452 * t502) - m(5) * (t289 * t349 - t488 * t560 - t490 * t561) - m(4) * (t387 * t433 + (t562 * t788 + t564 * t791) * t733);
t1294 = t1298 * t1170 + t1297 * t1173;
t1293 = t1299 * t1172 + t1300 * t1233;
t991 = t1257 * t788;
t326 = t521 * t1057 - t713 * t513 + t714 * t530;
t327 = t523 * t1057 + t713 * t514 + t714 * t532;
t1279 = -t326 * t791 + t327 * t788;
t328 = t515 * t1057 - t713 * t524 - t714 * t534;
t329 = t517 * t1057 - t713 * t526 + t714 * t535;
t1280 = -t328 * t791 + t329 * t788;
t1291 = t1279 + t1280;
t417 = t1057 * t656 + t713 * t648 + t714 * t664;
t876 = t788 * t326 + t327 * t791;
t1281 = -t417 * t790 + t787 * t876;
t418 = t1057 * t650 - t713 * t658 + t714 * t666;
t875 = t788 * t328 + t329 * t791;
t1282 = -t418 * t790 + t787 * t875;
t1290 = t1281 + t1282;
t1289 = t1279 / 0.2e1 + t1280 / 0.2e1;
t1288 = t1281 / 0.2e1 + t1282 / 0.2e1;
t1070 = t521 * t790;
t1275 = t513 * t786 - t530 * t789;
t369 = t1275 * t787 + t1070;
t1072 = t515 * t790;
t1273 = t524 * t786 + t534 * t789;
t372 = t1273 * t787 + t1072;
t510 = -Icges(6,5) * t712 - Icges(6,6) * t711 + Icges(6,3) * t1060;
t1076 = t510 * t790;
t693 = Icges(6,4) * t712;
t519 = -Icges(6,2) * t711 + Icges(6,6) * t1060 - t693;
t692 = Icges(6,4) * t711;
t527 = Icges(6,1) * t712 - Icges(6,5) * t1060 + t692;
t1274 = t519 * t786 - t527 * t789;
t366 = t1274 * t787 + t1076;
t1287 = t537 - t589;
t694 = Icges(6,4) * t714;
t520 = Icges(6,2) * t713 - Icges(6,6) * t1057 + t694;
t1092 = Icges(6,4) * t713;
t529 = Icges(6,1) * t714 - Icges(6,5) * t1057 + t1092;
t1015 = t713 * t520 + t714 * t529;
t511 = Icges(6,5) * t714 + Icges(6,6) * t713 - Icges(6,3) * t1057;
t865 = t711 * t519 - t712 * t527;
t1286 = t1015 + (-t510 * t788 - t511 * t791) * t787 + t865;
t1016 = -t713 * t519 + t714 * t527;
t1017 = t711 * t520 + t712 * t529;
t1285 = -t1016 - (t510 * t791 - t511 * t788) * t787 - t1017;
t1284 = -t366 - t369 - t372;
t1074 = t511 * t790;
t863 = t520 * t786 + t529 * t789;
t367 = t787 * t863 + t1074;
t1069 = t523 * t790;
t866 = t514 * t786 + t532 * t789;
t370 = t787 * t866 - t1069;
t1071 = t517 * t790;
t861 = -t526 * t786 + t535 * t789;
t373 = t787 * t861 - t1071;
t1283 = t367 + t370 + t373;
t445 = t672 * t1060 - t541 * t790;
t425 = -pkin(4) * t713 + t590 + t591;
t1278 = t425 * t788;
t555 = t791 * t589;
t1277 = -t589 + t1212;
t1272 = -m(4) / 0.2e1;
t1271 = -m(5) / 0.2e1;
t1270 = -m(6) / 0.2e1;
t882 = Icges(6,5) * t789 + Icges(6,6) * t786;
t646 = Icges(6,3) * t790 + t787 * t882;
t1082 = Icges(6,6) * t790;
t765 = Icges(6,4) * t1059;
t654 = Icges(6,2) * t1061 + t1082 + t765;
t1091 = Icges(6,4) * t786;
t889 = Icges(6,1) * t789 + t1091;
t662 = Icges(6,5) * t790 + t787 * t889;
t413 = t1060 * t646 - t654 * t711 - t662 * t712;
t1261 = t413 * t790;
t1260 = t671 * t788;
t1259 = t672 * t788;
t660 = Icges(3,4) * t1054 - Icges(3,2) * t1060 - Icges(3,6) * t791;
t1086 = Icges(3,2) * t787;
t780 = Icges(3,4) * t790;
t661 = Icges(3,6) * t788 + (t780 - t1086) * t791;
t1097 = Icges(3,4) * t787;
t746 = Icges(3,1) * t790 - t1097;
t669 = Icges(3,5) * t788 + t746 * t791;
t624 = t669 * t1054;
t742 = Icges(3,5) * t790 - Icges(3,6) * t787;
t653 = Icges(3,3) * t788 + t742 * t791;
t926 = t791 * t653 - t624;
t652 = Icges(3,5) * t1054 - Icges(3,6) * t1060 - Icges(3,3) * t791;
t766 = Icges(3,4) * t1060;
t668 = Icges(3,1) * t1054 - Icges(3,5) * t791 - t766;
t988 = -t668 * t1051 - t788 * t652;
t1258 = -t1057 * t660 - t1060 * t661 - t926 - t988;
t1161 = rSges(4,3) + pkin(7);
t1216 = t1161 * t787 + t938;
t459 = -t1216 * t788 + t1220 + t782;
t997 = t543 + t594;
t1224 = t790 * t997;
t383 = t1057 * t976 + t1224;
t363 = t713 * t383;
t1068 = t537 * t790;
t554 = t790 * t589;
t641 = t734 * t1060;
t986 = t671 * t1060 + t641;
t381 = -t554 + t986 + t1068;
t1028 = -t381 * t711 - t363;
t947 = t594 + t998;
t1225 = t790 * t947;
t315 = t1057 * t943 + t1225;
t294 = t713 * t315;
t929 = t1212 * t790;
t945 = -t1060 * t1257 + t641;
t313 = -t554 + t929 + t945;
t1030 = -t313 * t711 - t294;
t783 = t787 ^ 2;
t595 = t1056 * t783 + t711 * t790;
t596 = -t1048 * t783 - t790 * t713;
t1047 = (t390 * t595 + t391 * t596 + t1030) * t1194 + (t407 * t595 + t408 * t596 + t1028) * t1196;
t314 = t1277 * t790 + t945;
t382 = t1287 * t790 + t986;
t1113 = (t314 * t711 + t1030 + t294) * t1194 + (t382 * t711 + t1028 + t363) * t1196;
t1251 = t1047 - t1113;
t548 = t787 * t555;
t268 = -t548 + (t1212 * t791 - t788 * t947) * t787;
t312 = -t548 + (t537 * t791 - t788 * t997) * t787;
t506 = -t711 * t1057 + t1060 * t713;
t871 = -t381 * t791 + t383 * t788;
t881 = -t313 * t791 + t315 * t788;
t1114 = (t1061 * t881 - t506 * t252 + t268 * t559 - t596 * t450 - t595 * t452) * t1194 + (t1061 * t871 - t506 * t289 + t312 * t559 - t596 * t488 - t595 * t490) * t1196;
t736 = t893 * t790;
t944 = t734 * t1054 + t736 * t1060 + t790 * t644;
t974 = pkin(4) * t1053 + t1100 * t787 + t790 * t894;
t210 = (t788 * t974 - t1277) * t787 + t944;
t550 = t787 * t594;
t941 = -t736 - t974;
t645 = t893 * t1057;
t990 = t1257 * t791;
t946 = t645 - t990;
t211 = t550 + (t791 * t941 + t998) * t787 + (-t791 * t943 + t946) * t790;
t621 = -t1223 * t791 + t772;
t989 = -t621 + t645;
t996 = t644 * t1057 - t790 * t555;
t223 = (-t1260 * t787 + t1068) * t791 + (t787 * t989 - t1224) * t788 + t996;
t675 = rSges(5,2) * t787 + t790 * t895;
t269 = (t675 * t788 - t1287) * t787 + t944;
t973 = -t675 - t736;
t270 = t550 + (t791 * t973 + t543) * t787 + (-t791 * t976 + t989) * t790;
t169 = (t787 * t991 + t929) * t791 + (t787 * t946 - t1225) * t788 + t996;
t846 = t169 + t881;
t1115 = (t713 * t210 + t711 * t211 + (t268 * t790 + t787 * t846) * t786) * t1194 + (t713 * t269 + t711 * t270 + (t312 * t790 + (t223 + t871) * t787) * t786) * t1196;
t1250 = t1114 - t1115;
t886 = Icges(6,4) * t789 + Icges(6,2) * t786;
t810 = -t787 * t886 - t1082;
t605 = t810 * t788;
t611 = t662 * t788;
t842 = -t646 * t788 - t1274;
t261 = t842 * t790 + (t605 * t786 - t611 * t789 + t510) * t787;
t883 = Icges(5,5) * t789 + Icges(5,3) * t786;
t811 = -t787 * t883 + t1083;
t601 = t811 * t788;
t613 = t664 * t788;
t844 = t656 * t788 - t1275;
t263 = t844 * t790 + (t601 * t786 - t613 * t789 + t521) * t787;
t609 = t658 * t788;
t615 = t666 * t788;
t840 = -t650 * t788 + t1273;
t265 = -t840 * t790 + (t609 * t786 - t615 * t789 + t515) * t787;
t1249 = t261 + t263 + t265;
t606 = t810 * t791;
t612 = t662 * t791;
t841 = -t646 * t791 + t863;
t262 = t841 * t790 + (t606 * t786 - t612 * t789 - t511) * t787;
t602 = t811 * t791;
t614 = t664 * t791;
t843 = t656 * t791 + t866;
t264 = t843 * t790 + (t602 * t786 - t614 * t789 + t523) * t787;
t610 = t658 * t791;
t616 = t666 * t791;
t839 = -t650 * t791 - t861;
t266 = -t839 * t790 + (t610 * t786 - t616 * t789 + t517) * t787;
t1248 = t262 + t264 + t266;
t1008 = Icges(6,2) * t712 - t527 - t692;
t1012 = -Icges(6,1) * t711 - t519 + t693;
t566 = -Icges(6,5) * t711 + Icges(6,6) * t712;
t278 = t566 * t790 + (t1008 * t786 + t1012 * t789) * t787;
t1006 = Icges(5,3) * t712 - t530 - t689;
t1014 = -Icges(5,1) * t711 - t513 + t690;
t574 = -Icges(5,4) * t711 + Icges(5,6) * t712;
t280 = -t574 * t790 + (t1006 * t786 + t1014 * t789) * t787;
t1004 = Icges(4,2) * t712 + t534 + t695;
t1010 = -Icges(4,1) * t711 - t524 - t696;
t570 = -Icges(4,5) * t711 - Icges(4,6) * t712;
t282 = -t570 * t790 + (t1004 * t786 + t1010 * t789) * t787;
t1247 = t278 + t280 + t282;
t1007 = Icges(6,2) * t714 - t1092 - t529;
t1011 = -Icges(6,1) * t713 + t520 + t694;
t567 = -Icges(6,5) * t713 + Icges(6,6) * t714;
t279 = t567 * t790 + (t1007 * t786 + t1011 * t789) * t787;
t1005 = Icges(5,3) * t714 - t1089 - t532;
t1013 = -Icges(5,1) * t713 + t514 + t691;
t575 = -Icges(5,4) * t713 + Icges(5,6) * t714;
t281 = -t575 * t790 + (t1005 * t786 + t1013 * t789) * t787;
t1003 = Icges(4,2) * t714 - t535 + t697;
t1009 = -Icges(4,1) * t713 - t1096 - t526;
t571 = -Icges(4,5) * t713 - Icges(4,6) * t714;
t283 = -t571 * t790 + (t1003 * t786 + t1009 * t789) * t787;
t1246 = t279 + t281 + t283;
t655 = -Icges(6,6) * t787 + t790 * t886;
t663 = -Icges(6,5) * t787 + t790 * t889;
t1066 = t646 * t790;
t647 = -Icges(6,3) * t787 + t790 * t882;
t857 = t786 * t654 + t789 * t662;
t838 = t647 + t857;
t800 = -t787 * t838 - t1066;
t305 = t655 * t711 + t663 * t712 + t788 * t800;
t649 = Icges(5,6) * t787 + t790 * t883;
t665 = Icges(5,4) * t787 + t790 * t890;
t1064 = t656 * t790;
t657 = Icges(5,2) * t787 + t790 * t887;
t858 = t786 * t648 + t789 * t664;
t837 = -t657 + t858;
t798 = -t787 * t837 + t1064;
t306 = t649 * t711 + t665 * t712 + t788 * t798;
t659 = Icges(4,6) * t787 + t790 * t888;
t667 = Icges(4,5) * t787 + t790 * t891;
t1065 = t650 * t790;
t651 = Icges(4,3) * t787 + t790 * t884;
t856 = -t786 * t658 + t789 * t666;
t836 = t651 - t856;
t799 = t787 * t836 + t1065;
t307 = -t659 * t711 + t667 * t712 + t788 * t799;
t1245 = -t305 - t306 - t307;
t308 = t713 * t655 + t714 * t663 + t791 * t800;
t309 = t713 * t649 + t714 * t665 + t791 * t798;
t310 = -t713 * t659 + t714 * t667 + t791 * t799;
t1244 = -t308 - t309 - t310;
t357 = t838 * t790 + (t786 * t655 + t789 * t663 - t646) * t787;
t358 = t837 * t790 + (t786 * t649 + t789 * t665 + t656) * t787;
t359 = -t836 * t790 + (-t786 * t659 + t789 * t667 + t650) * t787;
t1243 = t359 + t358 + t357;
t1242 = -t413 + t1298;
t416 = -t1057 * t646 + t713 * t654 + t714 * t662;
t1241 = t416 + t417 + t418;
t436 = t787 * t857 + t1066;
t437 = t787 * t858 - t1064;
t438 = t787 * t856 - t1065;
t1240 = t436 + t437 + t438;
t1239 = t1283 * t791 + t1284 * t788;
t1236 = -t787 / 0.2e1;
t1235 = -t788 / 0.2e1;
t1234 = t790 / 0.2e1;
t1166 = t791 / 0.2e1;
t1231 = m(6) * t787;
t1062 = t783 * t786;
t1193 = m(6) / 0.4e1;
t1195 = m(5) / 0.4e1;
t951 = t1195 + t1193;
t1222 = t951 * (t1062 * t789 + t711 * t712 + t713 * t714);
t829 = t786 * t790 - t559;
t1221 = t951 * t829 * t1061;
t956 = qJD(2) * t791;
t424 = t711 * t952 - t1106 - t687;
t454 = t1164 * t711 - t1103 - t687;
t992 = -t590 - t592;
t1045 = (t390 * t714 + t391 * t712 + t424 * t713 + t425 * t711) * t1194 + (t407 * t714 + t408 * t712 + t454 * t713 - t711 * t992) * t1196;
t1099 = (-t488 * t712 - t490 * t714 + t560 * t711 + t561 * t713 + (t289 * t789 + t349 * t786) * t787) * t1196 + (-t450 * t712 - t452 * t714 + t501 * t711 + t502 * t713 + (t252 * t789 + t292 * t786) * t787) * t1194;
t1214 = t261 / 0.2e1 + t263 / 0.2e1 + t265 / 0.2e1;
t1213 = -t262 / 0.2e1 - t264 / 0.2e1 - t266 / 0.2e1;
t460 = t1216 * t791 + t1116 + t897;
t1211 = (t390 * t502 + t391 * t501 - t424 * t452 - t425 * t450) * t1270 + (t407 * t561 + t408 * t560 - t454 * t490 + t488 * t992) * t1271 + (-t562 * t593 + t564 * t587 + (-t459 * t791 - t460 * t788) * t733) * t1272;
t676 = rSges(4,3) * t787 + t790 * t896;
t385 = (t676 * t788 + t541) * t787;
t964 = t787 * rSges(4,2) * t1048 + rSges(4,3) * t1051;
t622 = -rSges(4,1) * t1052 * t787 + t964;
t386 = (-t672 * t791 - t622) * t790 + (-t676 * t791 + t544) * t787;
t447 = t1057 * t672 + t790 * t544;
t545 = t773 + (-t1161 * t790 + t835) * t788;
t546 = t774 + (-pkin(2) - t1108) * t1057 + t964;
t1210 = (t385 * t459 + t386 * t460 + t445 * t545 - t447 * t546) * t1272 + (t210 * t390 + t211 * t391 + t313 * t435 - t315 * t434) * t1270 + (t269 * t407 + t270 * t408 + t381 * t464 - t383 * t463) * t1271;
t1024 = t381 - t382;
t1029 = t313 - t314;
t1209 = t1024 * t488 * t1271 + t1029 * t450 * t1270;
t716 = (Icges(5,3) * t789 - t1088) * t787;
t720 = (Icges(6,2) * t789 - t1091) * t787;
t722 = (-Icges(4,2) * t789 - t1095) * t787;
t725 = -Icges(6,1) * t1061 + t765;
t726 = -Icges(5,1) * t1061 + t764;
t727 = (-Icges(4,1) * t786 - t1094) * t787;
t901 = t666 / 0.2e1 + t664 / 0.2e1 + t662 / 0.2e1;
t903 = t658 / 0.2e1 - t648 / 0.2e1 - t654 / 0.2e1;
t1208 = -t786 * (t722 / 0.2e1 - t716 / 0.2e1 - t720 / 0.2e1 + t901) + t789 * (t727 / 0.2e1 + t726 / 0.2e1 + t725 / 0.2e1 - t903);
t1098 = Icges(3,1) * t787;
t1207 = t786 * t903 - t789 * t901 - t647 / 0.2e1 + t651 / 0.2e1 + t657 / 0.2e1 - t780 + t1086 / 0.2e1 - t1098 / 0.2e1;
t743 = Icges(3,2) * t790 + t1097;
t1206 = t786 * (t659 / 0.2e1 - t649 / 0.2e1 - t655 / 0.2e1) - t789 * (t667 / 0.2e1 + t665 / 0.2e1 + t663 / 0.2e1) + t646 / 0.2e1 - t650 / 0.2e1 - t656 / 0.2e1 + t743 / 0.2e1 - t746 / 0.2e1;
t723 = -Icges(3,2) * t1054 - t766;
t724 = t743 * t791;
t892 = -t780 - t1098;
t728 = t892 * t788;
t729 = t892 * t791;
t1205 = (t791 * (t660 - t728) + (-t661 + t729) * t788) * t790 + (t791 * (t668 + t723) + (-t669 + t724) * t788) * t787;
t1204 = 0.2e1 * qJD(1);
t1203 = 0.4e1 * qJD(1);
t1202 = 0.2e1 * qJD(2);
t1200 = 2 * qJD(3);
t1199 = 4 * qJD(3);
t1198 = m(4) / 0.2e1;
t1189 = m(5) * (t223 * t312 + t269 * t381 - t270 * t383);
t1186 = m(5) * t1024 * t383;
t1182 = m(6) * (t846 * t790 + (-t210 * t791 - t211 * t788 - t268) * t787);
t1181 = m(6) * (t169 * t268 + t210 * t313 - t211 * t315);
t1179 = m(6) * t1029 * t315;
t304 = t313 * t1060;
t1176 = m(6) * (-t1060 * t314 + t304);
t1171 = t788 / 0.4e1;
t1167 = -t791 / 0.4e1;
t1109 = rSges(3,1) * t790;
t932 = pkin(1) + t1109;
t961 = rSges(3,2) * t1060 + t791 * rSges(3,3);
t597 = -t788 * t932 + t782 + t961;
t770 = rSges(3,2) * t1057;
t598 = -t770 + t932 * t791 + (rSges(3,3) + pkin(6)) * t788;
t749 = rSges(3,1) * t787 + rSges(3,2) * t790;
t735 = t749 * t788;
t737 = t749 * t791;
t1160 = m(3) * (t597 * t735 - t598 * t737);
t860 = -t541 * t791 - t544 * t788;
t331 = t860 * t790 + (-t1259 * t791 - t622 * t788) * t787;
t421 = t860 * t787;
t1159 = m(4) * (t331 * t421 + t385 * t445 - t386 * t447);
t1153 = m(4) * (t459 * t545 + t460 * t546);
t1152 = m(4) * (-t459 * t587 + t460 * t593);
t1142 = m(5) * (t407 * t454 - t408 * t992);
t1141 = m(5) * (t407 * t464 + t408 * t463);
t739 = t960 * t787;
t1137 = m(6) * (-t268 * t739 + t790 * t881);
t832 = t870 * t790;
t1132 = m(6) * (t832 + (-t434 * t788 - t435 * t791) * t787);
t1131 = m(6) * (t390 * t424 + t391 * t425);
t1129 = m(6) * (t832 + (t450 * t791 - t452 * t788) * t787);
t1128 = m(6) * (t390 * t435 + t391 * t434);
t1127 = m(6) * (t1057 * t315 + t304);
t1126 = m(6) * (t790 * t292 + (-t501 * t788 - t502 * t791) * t787);
t1125 = (-t424 * t791 - t1278) * t1231;
t963 = t960 * t1062;
t1124 = m(6) * (t790 * t829 - t1062 + t963);
t1067 = t559 * t790;
t1123 = m(6) * (-t1061 * t739 - t1067);
t1122 = m(6) * (t963 + t1067);
t1121 = m(6) * t506;
t1120 = (t711 * t791 - t713 * t788) * t1231;
t1112 = m(6) * qJD(1);
t1111 = m(6) * qJD(2);
t1110 = m(6) * qJD(5);
t378 = t391 * t713;
t397 = t408 * t713;
t1079 = t416 * t790;
t1063 = t660 * t787;
t1058 = t787 * t790;
t325 = -t1057 * t511 + t1015;
t1042 = -t325 + t1286;
t324 = t1057 * t510 + t1016;
t1041 = t324 + t1285;
t319 = -t1060 * t511 + t1017;
t1036 = -t319 - t1285;
t318 = t1060 * t510 - t865;
t1035 = t318 + t1286;
t227 = -t390 * t711 + t378;
t392 = t589 + t807 + t1217;
t1023 = t390 - t392;
t271 = -t407 * t711 + t397;
t409 = t589 + t818 + t967;
t1022 = t407 - t409;
t994 = t730 * t1060 + t790 * t584;
t987 = t669 * t1051 + t788 * t653;
t985 = t648 + t726;
t984 = t654 + t725;
t983 = -t658 + t727;
t981 = -t662 + t720;
t980 = -t664 + t716;
t979 = -t666 - t722;
t972 = -t676 - t755;
t971 = t788 * (pkin(7) * t1054 - t773) + t791 * (-pkin(2) * t1057 + t774);
t962 = t960 * t1058;
t959 = qJD(1) * t787;
t958 = qJD(1) * t790;
t957 = qJD(2) * t788;
t955 = qJD(3) * t787;
t954 = qJD(3) * t790;
t826 = m(6) * (-t712 * t788 - t714 * t791 + t1053);
t859 = -t595 * t791 - t596 * t788;
t948 = -t1121 / 0.2e1;
t317 = t790 * t948 + 0.2e1 * (t859 * t1193 - t826 / 0.4e1) * t787;
t953 = t317 * qJD(5);
t940 = -t755 + t973;
t936 = t1060 / 0.4e1;
t935 = t1057 / 0.4e1;
t925 = t661 * t787 - t652;
t921 = t788 * t644 - t791 * t645 + t971;
t919 = -t755 + t941;
t806 = -t787 * t842 + t1076;
t228 = t605 * t711 - t611 * t712 + t788 * t806;
t805 = -t787 * t841 - t1074;
t229 = t606 * t711 - t612 * t712 + t788 * t805;
t802 = -t787 * t844 + t1070;
t230 = t601 * t711 - t613 * t712 + t788 * t802;
t801 = -t787 * t843 + t1069;
t231 = t602 * t711 - t614 * t712 + t788 * t801;
t804 = t787 * t840 + t1072;
t232 = t609 * t711 - t615 * t712 + t788 * t804;
t803 = t787 * t839 + t1071;
t233 = t610 * t711 - t616 * t712 + t788 * t803;
t880 = t318 * t788 + t319 * t791;
t918 = (t880 + t1245 + t1297) * t1234 + ((t229 + t231 + t233) * t791 + (t228 + t230 + t232) * t788 + t1242) * t1173;
t234 = t713 * t605 - t714 * t611 + t791 * t806;
t235 = t713 * t606 - t714 * t612 + t791 * t805;
t236 = t713 * t601 - t714 * t613 + t791 * t802;
t237 = t713 * t602 - t714 * t614 + t791 * t801;
t238 = t713 * t609 - t714 * t615 + t791 * t804;
t239 = t713 * t610 - t714 * t616 + t791 * t803;
t877 = t788 * t324 + t325 * t791;
t917 = (t876 + t875 + t877 + t1244) * t1234 + ((t235 + t237 + t239) * t791 + (t234 + t236 + t238) * t788 + t1241) * t1173;
t916 = (t1248 * t791 + t1249 * t788 + t1240) * t1236 + (t1239 - t1243) * t1170;
t41 = m(5) * (-t312 * t506 + t381 * t595 - t383 * t596) + m(6) * (-t268 * t506 + t313 * t595 - t315 * t596);
t868 = t488 * t788 + t490 * t791;
t869 = t450 * t788 + t452 * t791;
t90 = m(5) * (t1061 * t868 + t289 * t559) + m(6) * (t1061 * t869 + t252 * t559);
t136 = m(5) * t271 + m(6) * t227;
t240 = t1008 * t711 + t1012 * t712 - t1060 * t566;
t241 = t1007 * t711 + t1011 * t712 - t1060 * t567;
t242 = t1006 * t711 + t1014 * t712 + t1060 * t574;
t243 = t1005 * t711 + t1013 * t712 + t1060 * t575;
t244 = t1004 * t711 + t1010 * t712 + t1060 * t570;
t245 = t1003 * t711 + t1009 * t712 + t1060 * t571;
t914 = (t241 + t243 + t245) * t1235 + (t240 + t242 + t244) * t1166;
t246 = t1008 * t713 + t1012 * t714 - t1057 * t566;
t247 = t1007 * t713 + t1011 * t714 - t1057 * t567;
t248 = t1006 * t713 + t1014 * t714 + t1057 * t574;
t249 = t1005 * t713 + t1013 * t714 + t1057 * t575;
t250 = t1004 * t713 + t1010 * t714 + t1057 * t570;
t251 = t1003 * t713 + t1009 * t714 + t1057 * t571;
t913 = (t246 + t248 + t250) * t1233 + (t247 + t249 + t251) * t1172;
t153 = t787 * t880 + t1261;
t912 = t153 / 0.2e1 + t1294;
t156 = t787 * t877 - t1079;
t911 = t156 / 0.2e1 + t1288;
t910 = t1234 * t1240 + t1236 * t1239;
t190 = -t318 * t791 + t319 * t788;
t909 = t190 / 0.2e1 + t1293;
t193 = -t324 * t791 + t325 * t788;
t908 = t193 / 0.2e1 + t1289;
t906 = t278 / 0.2e1 + t280 / 0.2e1 + t282 / 0.2e1;
t905 = -t279 / 0.2e1 - t281 / 0.2e1 - t283 / 0.2e1;
t715 = (-Icges(6,5) * t786 + Icges(6,6) * t789) * t787;
t717 = (-Icges(4,5) * t786 - Icges(4,6) * t789) * t787;
t721 = (-Icges(5,4) * t786 + Icges(5,6) * t789) * t787;
t899 = t717 / 0.2e1 + t721 / 0.2e1 - t715 / 0.2e1;
t885 = -Icges(3,5) * t787 - Icges(3,6) * t790;
t302 = t788 * t991 + t791 * t990 + t921;
t845 = t302 + t869;
t61 = -t1079 + (t1035 * t791 + t1036 * t788) * t787;
t825 = t61 / 0.2e1 - t911 + t1288;
t824 = -t914 - t918;
t823 = t913 - t917;
t58 = -t1261 + (t1041 * t791 + t1042 * t788) * t787;
t822 = t58 / 0.2e1 + t912 - t1294;
t442 = -t1057 * t661 + t987;
t80 = t1035 * t788 - t1036 * t791;
t797 = t987 * t1235 + t442 * t1172 - t80 / 0.2e1 + t908 + (-t624 + (t653 + t1063) * t791 + t988 + t1258) * t1233 - t1289;
t77 = t1041 * t788 - t1042 * t791;
t796 = t77 / 0.2e1 + (-t788 * (-t668 * t790 + t1063) - t791 * t652) * t1233 + (t791 * t925 + t442 - t987) * t1166 + t909 + (t788 * t925 + t1258 + t926) * t1172 - t1293;
t340 = -t1060 * t715 + t711 * t981 + t712 * t984;
t341 = t1060 * t721 + t711 * t980 + t712 * t985;
t342 = t1060 * t717 + t711 * t979 + t712 * t983;
t343 = -t1057 * t715 + t713 * t981 + t714 * t984;
t344 = t1057 * t721 + t713 * t980 + t714 * t985;
t345 = t1057 * t717 + t713 * t979 + t714 * t983;
t795 = -t1211 + (t343 + t344 + t345 + t1246) * t1171 + (t340 + t341 + t342 + t1247) * t1167;
t793 = -t1209 - (t193 + t1291) * t1060 / 0.4e1 + (t80 + t1291) * t936 + (t61 + t1290) * t1167 + (t156 + t1290) * t791 / 0.4e1 + (t190 + t77) * t935 + (t153 + t58) * t1171;
t792 = -t1210 + t1240 * t1173 + t1243 * t1170 + (-t1245 + t1249) * t936 + (-t1244 + t1248) * t935 + (t1242 + t1284) * t1054 / 0.4e1 + (t1241 + t1283) * t1051 / 0.4e1;
t751 = -rSges(3,2) * t787 + t1109;
t719 = t885 * t791;
t718 = t885 * t788;
t633 = t962 - t1058;
t565 = t972 * t791;
t563 = t972 * t788;
t547 = t584 * t1057;
t499 = t1120 / 0.2e1;
t498 = t1121 / 0.2e1;
t491 = t940 * t791;
t489 = t940 * t788;
t466 = -t1057 * t733 - t790 * t593;
t465 = t1060 * t733 + t587 * t790;
t457 = t1122 / 0.2e1;
t456 = t1123 / 0.2e1;
t453 = t919 * t791;
t451 = t919 * t788;
t432 = t1124 / 0.2e1;
t428 = (t587 * t791 - t593 * t788) * t787;
t419 = -t1259 * t788 + t622 * t791 + t971;
t400 = -t790 * t717 + (t786 * t979 + t789 * t983) * t787;
t399 = -t790 * t721 + (t786 * t980 + t789 * t985) * t787;
t398 = t790 * t715 + (t786 * t981 + t789 * t984) * t787;
t395 = t561 * t787 + t790 * t992;
t394 = t1060 * t732 + t586 * t790 + t994;
t389 = 0.4e1 * t1221;
t384 = t390 * t1060;
t376 = t499 + t948;
t375 = t498 - t1120 / 0.2e1;
t374 = t498 + t499;
t364 = -t1260 * t788 + t621 * t791 + t921;
t361 = (pkin(4) * t1062 + t787 * t966) * t791 - t425 * t790;
t360 = -pkin(4) * t595 + t1060 * t731 + t585 * t790 + t994;
t350 = 0.4e1 * t1222;
t330 = t547 + (t586 * t791 + t788 * t992) * t787;
t316 = (-t790 * t506 + t787 * t859) * t1194 + t826 * t1173;
t295 = t1125 / 0.2e1;
t288 = t456 + t457 - t1124 / 0.2e1;
t287 = t456 + t432 - t1122 / 0.2e1;
t286 = t457 + t432 - t1123 / 0.2e1;
t285 = t547 + ((-pkin(4) * t711 + t585) * t791 - t1278) * t787;
t254 = -t1057 * t391 + t384;
t215 = t1126 / 0.2e1;
t178 = t1127 / 0.2e1;
t167 = t1129 / 0.2e1;
t162 = -t252 * t739 + t790 * t869;
t161 = t1132 / 0.2e1;
t123 = -t238 * t791 + t239 * t788;
t122 = -t236 * t791 + t237 * t788;
t121 = -t234 * t791 + t235 * t788;
t120 = -t232 * t791 + t233 * t788;
t119 = -t230 * t791 + t231 * t788;
t118 = -t228 * t791 + t229 * t788;
t110 = t1137 / 0.2e1;
t102 = -t345 * t790 + (t250 * t788 + t251 * t791) * t787;
t101 = -t344 * t790 + (t248 * t788 + t249 * t791) * t787;
t100 = -t343 * t790 + (t246 * t788 + t247 * t791) * t787;
t99 = -t342 * t790 + (t244 * t788 + t245 * t791) * t787;
t98 = -t341 * t790 + (t242 * t788 + t243 * t791) * t787;
t97 = -t340 * t790 + (t240 * t788 + t241 * t791) * t787;
t83 = t1176 / 0.2e1;
t39 = t167 - t1132 / 0.2e1;
t38 = t167 + t161;
t37 = t161 - t1129 / 0.2e1;
t35 = t1208 * t787 - t899 * t790 + t1131 + t1142 + t1152;
t30 = t178 + t295 - t1176 / 0.2e1;
t29 = t178 + t83 - t1125 / 0.2e1;
t28 = t295 + t83 - t1127 / 0.2e1;
t27 = -t1206 * t787 - t1207 * t790 + t1128 + t1141 + t1153 + t1160;
t24 = t1182 / 0.2e1;
t19 = t110 + t215 - t1182 / 0.2e1;
t18 = t110 + t24 - t1126 / 0.2e1;
t17 = t215 + t24 - t1137 / 0.2e1;
t15 = t1043 + t1044;
t13 = t1045 - t1251;
t12 = t1047 + t1113 - t1045;
t11 = t1045 + t1251;
t10 = t788 * t913 + t791 * t914 - t1295;
t9 = t1099 - t1250;
t8 = t1114 + t1115 - t1099;
t7 = t1099 + t1250;
t6 = t788 * t796 + t791 * t797;
t5 = t1186 + t1179 + (t788 * t825 + t791 * t822) * t787;
t4 = t1159 + t1189 + t1181 + (t788 * t912 + t791 * t911 + t916) * t790 + (t788 * t918 + t791 * t917 - t910) * t787;
t3 = t795 + (-t153 / 0.4e1 - t58 / 0.4e1) * t788 + ((-t190 / 0.4e1 - t77 / 0.4e1) * t791 + (-t80 / 0.4e1 + t193 / 0.4e1) * t788) * t787 + (-t156 / 0.4e1 + t61 / 0.4e1) * t791 + t792 + t1209;
t2 = t793 + (-t345 / 0.4e1 - t344 / 0.4e1 - t343 / 0.4e1 - t283 / 0.4e1 - t281 / 0.4e1 - t279 / 0.4e1) * t788 + t792 + (t342 / 0.4e1 + t341 / 0.4e1 + t340 / 0.4e1 + t282 / 0.4e1 + t280 / 0.4e1 + t278 / 0.4e1) * t791 + t1211;
t1 = t795 + (-t438 / 0.2e1 - t437 / 0.2e1 - t436 / 0.2e1 + (-t310 / 0.4e1 - t309 / 0.4e1 - t308 / 0.4e1 - t266 / 0.4e1 - t264 / 0.4e1 - t262 / 0.4e1) * t791 + (-t307 / 0.4e1 - t306 / 0.4e1 - t305 / 0.4e1 - t265 / 0.4e1 - t263 / 0.4e1 - t261 / 0.4e1) * t788) * t787 + t793 + (t359 / 0.2e1 + t358 / 0.2e1 + t357 / 0.2e1 + (-t418 / 0.4e1 - t417 / 0.4e1 - t416 / 0.4e1 - t373 / 0.4e1 - t370 / 0.4e1 - t367 / 0.4e1) * t791 + (t372 / 0.4e1 + t369 / 0.4e1 + t366 / 0.4e1 - t412 / 0.4e1 - t411 / 0.4e1 + t413 / 0.4e1) * t788) * t790 + t1210;
t14 = [(-m(5) * t1022 * t408 / 0.4e1 - m(6) * t1023 * t391 / 0.4e1) * t1203 + t27 * qJD(2) + t35 * qJD(3) + t136 * qJD(4) + t254 * t1110, t27 * qJD(1) + t3 * qJD(3) + t15 * qJD(4) + t38 * qJD(5) + ((t390 * t453 + t391 * t451 - t434 * t450 - t435 * t452) * t1194 + (t407 * t491 + t408 * t489 - t463 * t488 - t464 * t490) * t1196 + (t459 * t565 + t460 * t563 - t545 * t564 - t546 * t562) * t1198) * t1202 + (-t305 / 0.2e1 + m(3) * (-t597 * t751 - t735 * t749) + (-t668 / 0.2e1 - t723 / 0.2e1) * t790 + (t660 / 0.2e1 - t728 / 0.2e1) * t787 - t306 / 0.2e1 + t742 * t1166 - t797 - t307 / 0.2e1 - t1214) * t956 + ((-t661 / 0.2e1 + t729 / 0.2e1) * t787 + t310 / 0.2e1 + m(3) * (-t598 * t751 + t737 * t749) + (t669 / 0.2e1 - t724 / 0.2e1) * t790 + t308 / 0.2e1 + t309 / 0.2e1 - t796 + t742 * t1172 - t1213) * t957, t35 * qJD(1) + t3 * qJD(2) + t11 * qJD(4) + t30 * qJD(5) + (-t1186 / 0.4e1 - t1179 / 0.4e1) * t1199 + ((t313 * t424 - t315 * t425 + t360 * t390 + t361 * t391) * t1194 + (t381 * t454 + t383 * t992 + t394 * t407 + t395 * t408) * t1196 + (-t445 * t587 - t447 * t593 + t459 * t465 + t460 * t466) * t1198) * t1200 + (-t399 - t398 - t400) * t954 + ((t343 / 0.2e1 + t344 / 0.2e1 + t345 / 0.2e1 - t822 - t905) * t791 + (t340 / 0.2e1 + t341 / 0.2e1 + t342 / 0.2e1 - t825 + t906) * t788) * t955, t136 * qJD(1) + t15 * qJD(2) + t11 * qJD(3) + t374 * qJD(5), t38 * qJD(2) + t30 * qJD(3) + t374 * qJD(4) + t1112 * t254; t6 * qJD(2) + t1 * qJD(3) - t16 * qJD(4) + t39 * qJD(5) + (t1022 * t488 * t1196 + t1023 * t450 * t1194) * t1204 + (-t1153 / 0.4e1 - t1141 / 0.4e1 - t1128 / 0.4e1 - t1160 / 0.4e1) * t1203 + t1207 * t958 + t1206 * t959, t6 * qJD(1) + t10 * qJD(3) + t90 * qJD(4) + t162 * t1110 - (t785 * t718 + (-t791 * t719 + t1205) * t788 + t118 + t120 + t119) * t956 / 0.2e1 + (m(4) * (t387 * t419 - t562 * t563 - t564 * t565) + m(3) * ((t788 * (rSges(3,1) * t1054 - t961) + t791 * (rSges(3,1) * t1051 + t788 * rSges(3,3) - t770)) * (-t788 * t735 - t737 * t791) + t960 * t751 * t749) + m(6) * (t252 * t302 - t450 * t451 - t452 * t453) + m(5) * (t289 * t364 - t488 * t489 - t490 * t491) + (t784 * t719 + (-t788 * t718 + t1205) * t791 + t122 + t121 + t123) * t1172) * qJD(2), t1 * qJD(1) + t10 * qJD(2) + t7 * qJD(4) + t19 * qJD(5) + (-t1159 / 0.4e1 - t1189 / 0.4e1 - t1181 / 0.4e1) * t1199 + ((t428 * t387 + t421 * t433 - t465 * t564 - t466 * t562 + (-t445 * t791 + t447 * t788) * t733) * t1198 + (t252 * t285 + t268 * t292 + t313 * t502 - t315 * t501 - t360 * t452 - t361 * t450) * t1194 + (t289 * t330 + t312 * t349 + t381 * t561 - t383 * t560 - t394 * t490 - t395 * t488) * t1196) * t1200 + ((t906 - t911) * t791 + (t905 - t912) * t788 - t916) * t954 + (t788 * t824 + t791 * t823 + t910) * t955 + ((-t97 / 0.2e1 - t98 / 0.2e1 - t99 / 0.2e1) * t791 + (t100 / 0.2e1 + t101 / 0.2e1 + t102 / 0.2e1) * t788) * qJD(3), t90 * qJD(2) + t7 * qJD(3) - 0.4e1 * qJD(4) * t1221 + t288 * qJD(5) - t1296, t39 * qJD(1) + t162 * t1111 + t19 * qJD(3) + t288 * qJD(4) + (-t739 * t790 - t633 + t962) * t1110; t2 * qJD(2) + t5 * qJD(3) + t12 * qJD(4) + t29 * qJD(5) + t899 * t958 + (-t1152 / 0.4e1 - t1142 / 0.4e1 - t1131 / 0.4e1) * t1203 + ((t1022 * t383 - t1024 * t408) * t1196 + (t1023 * t315 - t1029 * t391) * t1194) * t1204 - t1208 * t959, t2 * qJD(1) + t4 * qJD(3) + t8 * qJD(4) + t18 * qJD(5) + ((t169 * t252 - t210 * t452 - t211 * t450 + t268 * t302 + t313 * t453 - t315 * t451) * t1194 + (t223 * t289 - t269 * t490 - t270 * t488 + t312 * t364 + t381 * t491 - t383 * t489) * t1196 + (t331 * t387 - t385 * t564 - t386 * t562 + t419 * t421 + t445 * t565 - t447 * t563) * t1198) * t1202 + ((t908 + t1214) * t790 + t824) * t956 + ((t909 + t1213) * t790 - t823) * t957 + (((t121 / 0.2e1 + t122 / 0.2e1 + t123 / 0.2e1 + t366 / 0.2e1 + t369 / 0.2e1 + t372 / 0.2e1) * t791 + (t118 / 0.2e1 + t119 / 0.2e1 + t120 / 0.2e1 + t373 / 0.2e1 + t367 / 0.2e1 + t370 / 0.2e1) * t788) * t787 + t1295) * qJD(2), t5 * qJD(1) + t4 * qJD(2) + t41 * qJD(4) + ((t312 * t330 + t381 * t394 - t383 * t395) * t1195 + (t268 * t285 + t313 * t360 - t315 * t361) * t1193 + (t421 * t428 + t445 * t465 - t447 * t466) * m(4) / 0.4e1) * t1199 + (t399 / 0.2e1 + t398 / 0.2e1 + t400 / 0.2e1) * qJD(3) * t790 ^ 2 + ((t97 + t98 + t99) * t1172 + (t1246 * t791 + t1247 * t788) * t1170 + (t101 + t102 + t100) * t1166) * t955, t12 * qJD(1) + t8 * qJD(2) + t41 * qJD(3) + (-0.4e1 * t1222 + 0.2e1 * (t1196 + t1194) * (-t1061 * t506 + t595 * t713 + t596 * t711)) * qJD(4) + t953, qJD(1) * t29 + qJD(2) * t18 + qJD(4) * t317; (m(6) * (t392 * t711 + t227 - t378) + m(5) * (t409 * t711 + t271 - t397) - t136) * qJD(1) + t16 * qJD(2) + t13 * qJD(3) + t375 * qJD(5), t1296 + (m(6) * (t711 * t451 + t713 * t453) + m(5) * (t711 * t489 + t713 * t491) + 0.2e1 * ((t252 * t1194 + t289 * t1196) * t790 + (t845 * t1194 + (t364 + t868) * t1196) * t787) * t786 - t90) * qJD(2) + t9 * qJD(3) + t389 * qJD(4) + t287 * qJD(5), t13 * qJD(1) + t9 * qJD(2) + (m(5) * (t381 * t714 - t383 * t712 + t394 * t713 + t395 * t711 + (t312 * t789 + t330 * t786) * t787) + m(6) * (t313 * t714 - t315 * t712 + t360 * t713 + t361 * t711 + (t268 * t789 + t285 * t786) * t787) - t41) * qJD(3) + t350 * qJD(4) - t953, qJD(2) * t389 + qJD(3) * t350, qJD(1) * t375 + qJD(2) * t287 - qJD(3) * t317; (-t1060 * t392 - t254 + t384) * t1112 + t37 * qJD(2) + t28 * qJD(3) + t376 * qJD(4), t37 * qJD(1) + (t845 * t790 + (-t451 * t788 - t453 * t791 - t252) * t787 - t162) * t1111 + t17 * qJD(3) + t286 * qJD(4) + t633 * t1110, t28 * qJD(1) + t17 * qJD(2) + m(6) * (t790 * t285 + (-t360 * t791 - t361 * t788) * t787) * qJD(3) + t316 * qJD(4), qJD(1) * t376 + qJD(2) * t286 + qJD(3) * t316, t633 * t1111;];
Cq = t14;