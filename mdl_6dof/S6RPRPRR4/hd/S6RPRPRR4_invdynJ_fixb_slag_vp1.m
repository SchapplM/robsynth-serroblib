% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:28
% EndTime: 2019-03-09 03:46:25
% DurationCPUTime: 111.34s
% Computational Cost: add. (45519->1531), mult. (49856->1946), div. (0->0), fcn. (46088->10), ass. (0->698)
t1169 = Icges(5,1) + Icges(4,3);
t613 = sin(qJ(3));
t616 = cos(qJ(3));
t514 = Icges(4,5) * t616 - Icges(4,6) * t613;
t735 = Icges(5,4) * t616 - Icges(5,5) * t613;
t1168 = t514 - t735;
t610 = qJ(1) + pkin(10);
t596 = cos(t610);
t1167 = t1169 * t596;
t595 = sin(t610);
t957 = t595 * t616;
t959 = t595 * t613;
t1142 = t1167 + (-Icges(5,5) + Icges(4,6)) * t959 + (Icges(5,4) - Icges(4,5)) * t957;
t1148 = t1168 * t596 + t1169 * t595;
t982 = Icges(4,6) * t596;
t319 = Icges(4,4) * t957 - Icges(4,2) * t959 - t982;
t543 = Icges(5,6) * t959;
t994 = Icges(5,4) * t596;
t326 = Icges(5,2) * t957 - t543 + t994;
t1166 = t319 * t613 - t326 * t616;
t996 = Icges(4,4) * t613;
t520 = Icges(4,1) * t616 - t996;
t322 = Icges(4,5) * t595 + t520 * t596;
t980 = Icges(5,6) * t616;
t727 = -Icges(5,3) * t613 + t980;
t323 = Icges(5,5) * t595 - t596 * t727;
t1165 = -t322 * t957 - t323 * t959;
t517 = Icges(4,2) * t616 + t996;
t981 = Icges(5,6) * t613;
t726 = Icges(5,3) * t616 + t981;
t1164 = -t517 - t726;
t603 = Icges(4,4) * t616;
t519 = Icges(4,1) * t613 + t603;
t728 = Icges(5,2) * t613 + t980;
t1163 = t519 + t728;
t729 = Icges(5,2) * t616 - t981;
t1162 = t520 + t729;
t548 = Icges(4,4) * t959;
t987 = Icges(4,5) * t596;
t321 = Icges(4,1) * t957 - t548 - t987;
t986 = Icges(5,5) * t596;
t324 = Icges(5,6) * t957 - Icges(5,3) * t959 + t986;
t1139 = -t321 * t616 + t324 * t613 + t1166;
t736 = -Icges(4,2) * t613 + t603;
t1161 = t727 + t736;
t703 = t517 * t613 - t519 * t616;
t1151 = -t613 * t726 + t616 * t728 - t703;
t1160 = t1148 * t596 + t1165;
t954 = t596 * t616;
t955 = t596 * t613;
t1113 = t1148 * t595 + t322 * t954 + t323 * t955;
t1159 = t1142 * t595 - t321 * t954 + t324 * t955;
t1100 = -t1139 * t595 + t1142 * t596;
t320 = Icges(4,6) * t595 + t596 * t736;
t544 = Icges(5,6) * t955;
t995 = Icges(5,4) * t595;
t325 = -Icges(5,2) * t954 + t544 + t995;
t1099 = -t320 * t959 - t325 * t957 - t1160;
t1098 = -t319 * t955 + t326 * t954 - t1159;
t1097 = -t320 * t955 - t325 * t954 + t1113;
t1158 = t320 - t323;
t1157 = t1161 * qJD(3);
t1156 = t1162 * qJD(3);
t513 = Icges(4,5) * t613 + Icges(4,6) * t616;
t734 = Icges(5,4) * t613 + Icges(5,5) * t616;
t1155 = t513 - t734;
t1153 = t1164 * qJD(3);
t1152 = t1163 * qJD(3);
t1150 = -t1163 * t613 + t1164 * t616;
t961 = t513 * t596;
t237 = -t595 * t703 - t961;
t412 = t734 * t596;
t240 = t726 * t959 - t728 * t957 - t412;
t1149 = t237 - t240;
t411 = t734 * t595;
t962 = t513 * t595;
t1117 = t1151 * t596 - t411 + t962;
t1147 = t320 * t613 + t325 * t616;
t1096 = (t319 + t324) * t616 + (t321 + t326) * t613;
t1095 = t1158 * t616 + (t322 - t325) * t613;
t1146 = t1153 * t596 + (-t1161 * t595 + t982 - t986) * qJD(1);
t1145 = qJD(1) * t1158 + t1153 * t595;
t1144 = -t1152 * t596 + (-t1162 * t595 + t987 - t994) * qJD(1);
t1143 = t1152 * t595 + (-t596 * t729 - t322 + t995) * qJD(1);
t1141 = qJD(1) * t1155 + qJD(3) * t1150 + t1156 * t616 - t1157 * t613;
t1140 = t1155 * qJD(3);
t1138 = t322 * t616 + t323 * t613 - t1147;
t1137 = t1097 * t595 - t1098 * t596;
t1136 = t1099 * t595 - t1100 * t596;
t1034 = -rSges(7,3) - pkin(3);
t612 = sin(qJ(5));
t805 = -pkin(5) * t612 - qJ(4);
t688 = t613 * t805 - pkin(2);
t1135 = t1034 * t616 + t688;
t1134 = t1151 * qJD(1) - qJD(3) * t1168;
t617 = cos(qJ(1));
t608 = t617 * pkin(1);
t1133 = t1117 * qJD(1);
t1132 = t1149 * qJD(1);
t1131 = t1148 * qJD(1);
t614 = sin(qJ(1));
t1027 = t614 * pkin(1);
t455 = rSges(3,1) * t595 + rSges(3,2) * t596;
t434 = -t455 - t1027;
t602 = t613 * qJ(4);
t1081 = t616 * pkin(3) + t602;
t420 = t1081 * t595;
t587 = t596 * pkin(7);
t461 = t595 * pkin(2) - t587;
t450 = qJD(1) * t461;
t876 = qJD(1) * t596;
t570 = pkin(7) * t876;
t599 = qJD(4) * t613;
t529 = t596 * t599;
t871 = qJD(3) * t616;
t827 = t596 * t871;
t891 = qJ(4) * t827 + t529;
t1129 = qJD(1) * t420 + t450 + t570 + t891;
t598 = qJD(5) * t613;
t1077 = qJD(1) + t598;
t611 = qJ(5) + qJ(6);
t600 = sin(t611);
t601 = cos(t611);
t951 = t601 * t613;
t357 = t595 * t951 + t596 * t600;
t953 = t600 * t613;
t358 = -t595 * t953 + t596 * t601;
t910 = t358 * rSges(7,1) - t357 * rSges(7,2);
t189 = rSges(7,3) * t957 - t910;
t446 = t596 * pkin(4) - pkin(8) * t957;
t948 = t612 * t613;
t847 = t595 * t948;
t615 = cos(qJ(5));
t1033 = pkin(5) * t615;
t591 = pkin(4) + t1033;
t618 = -pkin(9) - pkin(8);
t943 = t616 * t618;
t892 = -t596 * t591 - t595 * t943;
t255 = pkin(5) * t847 + t446 + t892;
t868 = qJD(5) * t616;
t873 = qJD(3) * t596;
t440 = -t595 * t868 + t873;
t867 = qJD(6) * t616;
t316 = -t595 * t867 + t440;
t761 = rSges(7,1) * t600 + rSges(7,2) * t601;
t372 = rSges(7,3) * t613 - t616 * t761;
t1024 = -pkin(8) - t618;
t947 = t612 * t616;
t437 = -pkin(5) * t947 + t1024 * t613;
t489 = qJD(6) * t613 + t1077;
t1032 = pkin(8) * t613;
t531 = pkin(3) * t613 - qJ(4) * t616;
t804 = -t531 - t1032;
t667 = t804 * t873 + t529;
t1128 = t1077 * t255 + t189 * t489 + t316 * t372 + t437 * t440 - t667;
t946 = t613 * t615;
t397 = t595 * t946 + t596 * t612;
t398 = t596 * t615 - t847;
t905 = t398 * rSges(6,1) - t397 * rSges(6,2);
t221 = rSges(6,3) * t957 - t905;
t764 = rSges(6,1) * t612 + rSges(6,2) * t615;
t426 = rSges(6,3) * t613 - t616 * t764;
t1127 = t1077 * t221 + t426 * t440 - t667;
t1126 = t596 ^ 2;
t620 = qJD(1) ^ 2;
t857 = t620 * t608;
t1125 = qJD(3) * t1136 + t1132;
t1124 = qJD(3) * t1137 + t1133;
t1123 = qJD(3) * t1139 + t1143 * t613 - t1145 * t616;
t1122 = qJD(3) * t1138 + t1144 * t613 + t1146 * t616;
t872 = qJD(3) * t613;
t831 = t595 * t872;
t874 = qJD(1) * t616;
t833 = t596 * t874;
t669 = -t831 + t833;
t653 = -t1077 * t612 + t615 * t871;
t875 = qJD(1) * t613;
t785 = qJD(5) + t875;
t696 = t596 * t785;
t201 = t595 * t653 + t615 * t696;
t826 = t612 * t871;
t654 = t1077 * t615 + t826;
t202 = t595 * t654 + t612 * t696;
t766 = rSges(6,1) * t202 + rSges(6,2) * t201;
t115 = rSges(6,3) * t669 + t766;
t527 = t596 * t868;
t860 = qJDD(5) * t616;
t1083 = qJD(1) * t527 + t595 * t860;
t863 = qJD(1) * qJD(3);
t566 = t595 * t863;
t448 = -qJDD(3) * t596 + t566;
t257 = -qJD(5) * t831 + t1083 + t448;
t454 = (-rSges(6,1) * t615 + rSges(6,2) * t612) * t616;
t606 = t616 * rSges(6,3);
t272 = qJD(5) * t454 + (t613 * t764 + t606) * qJD(3);
t463 = qJD(3) * t868 + qJDD(5) * t613 + qJDD(1);
t445 = t595 * pkin(4) + pkin(8) * t954;
t505 = pkin(8) * t831;
t313 = qJD(1) * t445 - t505;
t869 = qJD(4) * t616;
t441 = qJD(3) * t1081 - t869;
t862 = qJD(3) * qJD(4);
t1090 = qJDD(4) * t613 + t616 * t862;
t682 = t1090 * t596 + t448 * t531 - t857;
t393 = t595 * t871 + t596 * t875;
t506 = pkin(3) * t831;
t824 = t595 * t599;
t218 = pkin(3) * t833 + qJ(4) * t393 - t506 + t824;
t462 = t596 * pkin(2) + t595 * pkin(7);
t442 = t462 * qJD(1);
t691 = -t218 - t442 - t824;
t809 = -t461 - t1027;
t775 = -t420 + t809;
t699 = t446 + t775;
t942 = t616 * qJD(3) ^ 2;
t622 = t448 * t1032 + (-pkin(8) * t942 - qJD(3) * t441) * t596 + t699 * qJDD(1) + (-t313 + t691) * qJD(1) + t682;
t38 = -t1077 * t115 - t463 * t221 + t257 * t426 - t440 * t272 + t622;
t1121 = t38 * t595;
t828 = t596 * t872;
t1119 = -t1134 * t595 + t1141 * t596;
t1118 = t1134 * t596 + t1141 * t595;
t1116 = t1142 * qJD(1) + qJD(3) * t1096 + t1143 * t616 + t1145 * t613;
t1115 = -qJD(3) * t1095 + t1144 * t616 - t1146 * t613 + t1131;
t1114 = t1142 + t1147;
t1018 = pkin(1) * qJD(1);
t851 = t614 * t1018;
t1112 = -qJD(1) * t446 + t1129 + t851;
t1111 = -t1140 * t596 + (-t1168 * t595 - t1138 + t1167) * qJD(1);
t1110 = qJD(1) * t1139 - t1140 * t595 + t1131;
t1080 = t608 + t462;
t424 = pkin(3) * t954 + t596 * t602;
t1087 = t424 + t1080;
t582 = t595 * rSges(5,1);
t339 = -rSges(5,2) * t954 + rSges(5,3) * t955 + t582;
t177 = t339 + t1087;
t1107 = t445 + t1087;
t730 = Icges(7,5) * t600 + Icges(7,6) * t601;
t647 = -Icges(7,3) * t613 + t616 * t730;
t989 = Icges(7,4) * t600;
t732 = Icges(7,2) * t601 + t989;
t649 = -Icges(7,6) * t613 + t616 * t732;
t988 = Icges(7,4) * t601;
t737 = Icges(7,1) * t600 + t988;
t651 = -Icges(7,5) * t613 + t616 * t737;
t119 = -t357 * t649 + t358 * t651 - t647 * t957;
t575 = qJD(3) * t595;
t439 = t575 + t527;
t315 = t596 * t867 + t439;
t355 = -t595 * t600 + t596 * t951;
t356 = t595 * t601 + t596 * t953;
t178 = Icges(7,5) * t356 + Icges(7,6) * t355 + Icges(7,3) * t954;
t990 = Icges(7,4) * t356;
t181 = Icges(7,2) * t355 + Icges(7,6) * t954 + t990;
t341 = Icges(7,4) * t355;
t184 = Icges(7,1) * t356 + Icges(7,5) * t954 + t341;
t75 = t178 * t957 + t357 * t181 - t358 * t184;
t180 = -Icges(7,5) * t358 + Icges(7,6) * t357 + Icges(7,3) * t957;
t343 = Icges(7,4) * t358;
t183 = Icges(7,2) * t357 + Icges(7,6) * t957 - t343;
t342 = Icges(7,4) * t357;
t185 = Icges(7,1) * t358 - Icges(7,5) * t957 - t342;
t76 = t180 * t957 + t183 * t357 + t185 * t358;
t37 = t119 * t489 + t315 * t75 - t316 * t76;
t731 = Icges(6,5) * t612 + Icges(6,6) * t615;
t648 = -Icges(6,3) * t613 + t616 * t731;
t992 = Icges(6,4) * t612;
t733 = Icges(6,2) * t615 + t992;
t650 = -Icges(6,6) * t613 + t616 * t733;
t991 = Icges(6,4) * t615;
t738 = Icges(6,1) * t612 + t991;
t652 = -Icges(6,5) * t613 + t616 * t738;
t128 = -t397 * t650 + t398 * t652 - t648 * t957;
t843 = t596 * t946;
t960 = t595 * t612;
t395 = t843 - t960;
t845 = t596 * t948;
t958 = t595 * t615;
t396 = t845 + t958;
t208 = Icges(6,5) * t396 + Icges(6,6) * t395 + Icges(6,3) * t954;
t993 = Icges(6,4) * t396;
t211 = Icges(6,2) * t395 + Icges(6,6) * t954 + t993;
t377 = Icges(6,4) * t395;
t214 = Icges(6,1) * t396 + Icges(6,5) * t954 + t377;
t79 = t208 * t957 + t397 * t211 - t398 * t214;
t210 = -Icges(6,5) * t398 + Icges(6,6) * t397 + Icges(6,3) * t957;
t379 = Icges(6,4) * t398;
t213 = Icges(6,2) * t397 + Icges(6,6) * t957 - t379;
t378 = Icges(6,4) * t397;
t215 = Icges(6,1) * t398 - Icges(6,5) * t957 - t378;
t80 = t210 * t957 + t213 * t397 + t215 * t398;
t41 = t1077 * t128 + t439 * t79 - t440 * t80;
t127 = -t395 * t650 - t396 * t652 - t648 * t954;
t77 = t208 * t954 + t395 * t211 + t396 * t214;
t78 = t210 * t954 + t395 * t213 - t215 * t396;
t40 = t1077 * t127 + t439 * t77 - t440 * t78;
t118 = -t355 * t649 - t356 * t651 - t647 * t954;
t73 = t178 * t954 + t355 * t181 + t356 * t184;
t74 = t180 * t954 + t355 * t183 - t185 * t356;
t36 = t118 * t489 + t315 * t73 - t316 * t74;
t717 = t213 * t615 - t215 * t612;
t87 = t210 * t613 - t616 * t717;
t719 = t183 * t601 - t185 * t600;
t83 = t180 * t613 - t616 * t719;
t1094 = t1116 * t1126 + (t1111 * t595 + (-t1110 + t1115) * t596) * t595;
t1093 = t1110 * t1126 + (t1115 * t595 + (-t1111 + t1116) * t596) * t595;
t812 = -t872 / 0.2e1;
t834 = t595 * t874;
t1092 = t596 * t812 - t834 / 0.2e1;
t817 = t876 / 0.2e1;
t1091 = t595 * t812 + t616 * t817;
t1089 = t828 + t834;
t950 = t601 * t616;
t952 = t600 * t616;
t1088 = rSges(7,1) * t952 + rSges(7,2) * t950;
t877 = qJD(1) * t595;
t581 = t595 * rSges(4,3);
t338 = rSges(4,1) * t954 - rSges(4,2) * t955 + t581;
t268 = t338 + t1080;
t187 = t356 * rSges(7,1) + t355 * rSges(7,2) + rSges(7,3) * t954;
t692 = pkin(5) * t845 + t595 * t591 - t596 * t943;
t254 = t692 - t445;
t690 = t420 * t575 + t424 * t873 + qJD(2) - t869;
t663 = t445 * t873 - t446 * t575 + t690;
t48 = t187 * t316 + t189 * t315 + t254 * t440 + t255 * t439 + t663;
t932 = t189 + t255;
t1086 = t48 * t932;
t671 = t699 * qJD(1);
t64 = -t1128 + t671;
t906 = t372 + t437;
t1085 = t64 * t906;
t456 = t596 * rSges(3,1) - rSges(3,2) * t595;
t435 = t456 + t608;
t577 = pkin(5) * t948;
t1082 = t1024 * t616 + t577;
t806 = -pkin(2) - t602;
t1035 = -rSges(6,3) - pkin(3);
t859 = -pkin(8) + t1035;
t1079 = t859 * t616 + t806;
t803 = -pkin(8) * t616 - t1081;
t604 = t613 * rSges(5,3);
t760 = -rSges(5,2) * t616 + t604;
t1028 = g(2) * t595;
t1066 = (g(1) * t596 + t1028) * t613;
t609 = qJD(5) + qJD(6);
t949 = t609 * t613;
t791 = qJD(1) + t949;
t655 = -t600 * t791 + t601 * t871;
t790 = t609 + t875;
t702 = t595 * t790;
t160 = t596 * t655 - t601 * t702;
t656 = t600 * t871 + t601 * t791;
t161 = t596 * t656 - t600 * t702;
t938 = t161 * rSges(7,1) + t160 * rSges(7,2);
t103 = -rSges(7,3) * t1089 + t938;
t701 = t596 * t790;
t158 = t595 * t655 + t601 * t701;
t159 = t595 * t656 + t600 * t701;
t763 = rSges(7,1) * t159 + rSges(7,2) * t158;
t102 = rSges(7,3) * t669 + t763;
t870 = qJD(3) * t618;
t825 = t613 * t870;
t121 = t595 * t825 + t505 + (qJD(5) * t397 + t595 * t826) * pkin(5) + ((-pkin(4) + t591) * t595 + t1082 * t596) * qJD(1);
t571 = pkin(4) * t876;
t695 = t612 * t785;
t487 = pkin(5) * t843;
t844 = t596 * t947;
t488 = pkin(5) * t844;
t741 = qJD(3) * t488 + qJD(5) * t487 + t591 * t876 + t596 * t825 + t618 * t834;
t122 = pkin(8) * t828 - t571 + (-pkin(5) * t695 + pkin(8) * t874) * t595 + t741;
t659 = -qJD(3) * t949 + qJDD(6) * t616;
t447 = qJDD(3) * t595 + t596 * t863;
t837 = t596 * t860 + t447;
t146 = t596 * t659 - t609 * t834 + t837;
t147 = t566 + (qJD(1) * t867 - qJDD(3)) * t596 + t659 * t595 + t1083;
t256 = -qJD(5) * t1089 + t837;
t314 = -pkin(8) * t1089 + t571;
t835 = t595 * t875;
t217 = -pkin(3) * t1089 - qJ(4) * t835 + t891;
t662 = -qJDD(4) * t616 + t217 * t873 + t218 * t575 + t447 * t420 + t613 * t862 + qJDD(2);
t901 = -t424 - t445;
t627 = t313 * t575 + t314 * t873 - t447 * t446 + t448 * t901 + t662;
t11 = t102 * t315 + t103 * t316 + t121 * t439 + t122 * t440 + t146 * t189 - t147 * t187 - t254 * t257 + t255 * t256 + t627;
t933 = t187 + t254;
t1064 = t11 * t933 + t48 * (t103 + t122);
t913 = -Icges(5,3) * t957 + t326 - t543;
t915 = t728 * t595 + t324;
t1063 = -t613 * t913 - t616 * t915;
t918 = -Icges(4,2) * t957 + t321 - t548;
t920 = t519 * t595 + t319;
t1062 = -t613 * t918 - t616 * t920;
t428 = (Icges(7,2) * t600 - t988) * t616;
t640 = t315 * (-Icges(7,2) * t356 + t184 + t341) - t316 * (Icges(7,2) * t358 - t185 + t342) + t489 * (-t651 + t428);
t429 = (-Icges(7,1) * t601 + t989) * t616;
t641 = t315 * (-Icges(7,1) * t355 + t181 + t990) - t316 * (-Icges(7,1) * t357 + t183 - t343) + t489 * (-t649 - t429);
t452 = (Icges(6,2) * t612 - t991) * t616;
t638 = t439 * (-Icges(6,2) * t396 + t214 + t377) - t440 * (Icges(6,2) * t398 - t215 + t378) + t1077 * (-t652 + t452);
t453 = (-Icges(6,1) * t615 + t992) * t616;
t639 = t439 * (-Icges(6,1) * t395 + t211 + t993) - t440 * (-Icges(6,1) * t397 + t213 - t379) + t1077 * (-t650 - t453);
t1061 = m(5) / 0.2e1;
t1060 = m(6) / 0.2e1;
t1059 = m(7) / 0.2e1;
t1058 = t146 / 0.2e1;
t1057 = t147 / 0.2e1;
t1056 = t256 / 0.2e1;
t1055 = t257 / 0.2e1;
t1054 = -t315 / 0.2e1;
t1053 = t315 / 0.2e1;
t1052 = -t316 / 0.2e1;
t1051 = t316 / 0.2e1;
t354 = qJD(3) * t867 + qJDD(6) * t613 + t463;
t1050 = t354 / 0.2e1;
t1049 = -t439 / 0.2e1;
t1048 = t439 / 0.2e1;
t1047 = -t440 / 0.2e1;
t1046 = t440 / 0.2e1;
t1045 = t447 / 0.2e1;
t1044 = t448 / 0.2e1;
t1043 = t463 / 0.2e1;
t1042 = -t489 / 0.2e1;
t1041 = t489 / 0.2e1;
t1040 = -t1077 / 0.2e1;
t1039 = t1077 / 0.2e1;
t1036 = t613 / 0.2e1;
t1030 = g(1) * t595;
t1023 = rSges(4,1) * t616;
t1022 = rSges(7,1) * t601;
t1020 = rSges(5,2) * t613;
t1017 = pkin(5) * qJD(5);
t101 = Icges(7,1) * t161 + Icges(7,4) * t160 - Icges(7,5) * t1089;
t720 = t181 * t601 + t184 * t600;
t97 = Icges(7,5) * t161 + Icges(7,6) * t160 - Icges(7,3) * t1089;
t99 = Icges(7,4) * t161 + Icges(7,2) * t160 - Icges(7,6) * t1089;
t20 = (qJD(3) * t720 + t97) * t613 + (qJD(3) * t178 + (-t184 * t609 - t99) * t601 + (t181 * t609 - t101) * t600) * t616;
t1016 = t20 * t315;
t100 = Icges(7,1) * t159 + Icges(7,4) * t158 + Icges(7,5) * t669;
t96 = Icges(7,5) * t159 + Icges(7,6) * t158 + Icges(7,3) * t669;
t98 = Icges(7,4) * t159 + Icges(7,2) * t158 + Icges(7,6) * t669;
t21 = (qJD(3) * t719 + t96) * t613 + (qJD(3) * t180 + (t185 * t609 - t98) * t601 + (t183 * t609 - t100) * t600) * t616;
t1015 = t21 * t316;
t203 = t596 * t653 - t785 * t958;
t204 = -t595 * t695 + t596 * t654;
t110 = Icges(6,5) * t204 + Icges(6,6) * t203 - Icges(6,3) * t1089;
t112 = Icges(6,4) * t204 + Icges(6,2) * t203 - Icges(6,6) * t1089;
t114 = Icges(6,1) * t204 + Icges(6,4) * t203 - Icges(6,5) * t1089;
t718 = t211 * t615 + t214 * t612;
t28 = (qJD(3) * t718 + t110) * t613 + (qJD(3) * t208 - t112 * t615 - t114 * t612 + (t211 * t612 - t214 * t615) * qJD(5)) * t616;
t1014 = t28 * t439;
t109 = Icges(6,5) * t202 + Icges(6,6) * t201 + Icges(6,3) * t669;
t111 = Icges(6,4) * t202 + Icges(6,2) * t201 + Icges(6,6) * t669;
t113 = Icges(6,1) * t202 + Icges(6,4) * t201 + Icges(6,5) * t669;
t29 = (qJD(3) * t717 + t109) * t613 + (qJD(3) * t210 - t111 * t615 - t113 * t612 + (t213 * t612 + t215 * t615) * qJD(5)) * t616;
t1013 = t29 * t440;
t84 = -t1127 + t671;
t1008 = t596 * t84;
t605 = t616 * rSges(7,3);
t82 = t178 * t613 - t616 * t720;
t1007 = t82 * t146;
t1006 = t83 * t147;
t86 = t208 * t613 - t616 * t718;
t1005 = t86 * t256;
t1004 = t87 * t257;
t1002 = -rSges(5,3) - qJ(4);
t708 = -t600 * t651 - t601 * t649;
t138 = -t613 * t647 - t616 * t708;
t359 = Icges(7,3) * t616 + t613 * t730;
t427 = (-Icges(7,5) * t601 + Icges(7,6) * t600) * t616;
t222 = qJD(3) * t359 + t427 * t609;
t361 = Icges(7,6) * t616 + t613 * t732;
t223 = qJD(3) * t361 + t428 * t609;
t363 = Icges(7,5) * t616 + t613 * t737;
t224 = qJD(3) * t363 + t429 * t609;
t61 = (qJD(3) * t708 + t222) * t613 + (-qJD(3) * t647 + (t609 * t651 - t223) * t601 + (-t609 * t649 - t224) * t600) * t616;
t1001 = t138 * t354 + t61 * t489;
t707 = -t612 * t652 - t615 * t650;
t149 = -t613 * t648 - t616 * t707;
t399 = Icges(6,3) * t616 + t613 * t731;
t451 = (-Icges(6,5) * t615 + Icges(6,6) * t612) * t616;
t269 = qJD(3) * t399 + qJD(5) * t451;
t401 = Icges(6,6) * t616 + t613 * t733;
t270 = qJD(3) * t401 + qJD(5) * t452;
t403 = Icges(6,5) * t616 + t613 * t738;
t271 = qJD(3) * t403 + qJD(5) * t453;
t71 = (qJD(3) * t707 + t269) * t613 + (-qJD(3) * t648 - t270 * t615 - t271 * t612 + (-t612 * t650 + t615 * t652) * qJD(5)) * t616;
t1000 = t1077 * t71 + t149 * t463;
t999 = t102 * t954 + t187 * t831;
t974 = qJD(3) * t64;
t883 = rSges(4,2) * t959 + t596 * rSges(4,3);
t337 = rSges(4,1) * t957 - t883;
t776 = -t337 + t809;
t533 = rSges(4,1) * t613 + rSges(4,2) * t616;
t832 = t533 * t873;
t168 = qJD(1) * t776 - t832;
t971 = t168 * t595;
t970 = t168 * t596;
t169 = qJD(1) * t268 - t533 * t575;
t423 = t533 * t596;
t969 = t169 * t423;
t508 = t609 * t616;
t945 = t613 * t618;
t944 = t615 * t616;
t940 = t102 + t121;
t241 = (rSges(7,2) * t600 - t1022) * t508 + (t613 * t761 + t605) * qJD(3);
t931 = t241 * t957 + t372 * t833;
t930 = t204 * rSges(6,1) + t203 * rSges(6,2);
t849 = t615 * t1017;
t312 = qJD(3) * t1082 - t616 * t849;
t925 = -t241 - t312;
t919 = -t519 * t596 - t320;
t917 = -t517 * t596 + t322;
t916 = -t728 * t596 + t323;
t914 = Icges(5,3) * t954 + t325 + t544;
t911 = -t339 - t424;
t252 = t355 * rSges(7,1) - t356 * rSges(7,2);
t253 = t357 * rSges(7,1) + t358 * rSges(7,2);
t909 = t595 * t420 + t596 * t424;
t541 = qJ(4) * t954;
t421 = -pkin(3) * t955 + t541;
t904 = qJD(1) * t421 + t595 * t869;
t900 = -t760 * qJD(3) - t441;
t444 = t531 * t877;
t507 = pkin(8) * t835;
t899 = t444 + t507;
t898 = t1088 * t595;
t897 = t1088 * t596;
t848 = t595 * t947;
t855 = rSges(6,2) * t944;
t896 = rSges(6,1) * t848 + t595 * t855;
t895 = rSges(6,1) * t844 + t596 * t855;
t894 = pkin(5) * t848 + t595 * t945;
t893 = t596 * t945 + t488;
t890 = rSges(4,2) * t835 + rSges(4,3) * t876;
t889 = -t726 + t729;
t888 = -t727 - t728;
t887 = -t517 + t520;
t886 = t519 + t736;
t759 = rSges(5,3) * t616 + t1020;
t885 = -t531 + t759;
t884 = -t1081 - t760;
t418 = rSges(5,2) * t959 + rSges(5,3) * t957;
t422 = rSges(5,2) * t955 + rSges(5,3) * t954;
t381 = t397 * pkin(5);
t882 = t595 ^ 2 + t1126;
t879 = qJD(1) * t514;
t878 = qJD(1) * t735;
t864 = -m(5) - m(6) - m(7);
t858 = pkin(5) * t944;
t850 = t612 * t1017;
t443 = t531 * t575;
t637 = qJD(1) * t1107 - t443 - t505 + t824;
t65 = t1077 * t254 + t187 * t489 - t315 * t372 - t437 * t439 + t637;
t842 = t65 * t876;
t219 = t396 * rSges(6,1) + t395 * rSges(6,2) + rSges(6,3) * t954;
t85 = t1077 * t219 - t426 * t439 + t637;
t841 = t85 * t876;
t840 = t596 * t217 + t595 * t218 + t420 * t876;
t539 = qJ(4) * t957;
t417 = -pkin(3) * t959 + t539;
t839 = t417 * t575 + t421 * t873 + t599;
t371 = rSges(7,1) * t953 + rSges(7,2) * t951 + t605;
t425 = rSges(6,1) * t948 + rSges(6,2) * t946 + t606;
t823 = t957 / 0.2e1;
t822 = t954 / 0.2e1;
t821 = -pkin(2) - t1023;
t816 = -t575 / 0.2e1;
t815 = t575 / 0.2e1;
t814 = -t873 / 0.2e1;
t813 = t873 / 0.2e1;
t811 = t871 / 0.2e1;
t807 = t587 - t1027;
t802 = rSges(5,1) * t596 - rSges(5,3) * t959;
t801 = t252 * t316 + t315 * t253;
t565 = rSges(7,2) * t952;
t430 = -rSges(7,1) * t950 + t565;
t800 = t489 * t252 - t315 * t430;
t799 = -t253 * t489 - t316 * t430;
t794 = qJD(3) * t900;
t793 = t882 * qJD(3);
t792 = -qJD(1) * t417 + t596 * t869;
t789 = t596 * t859;
t783 = t596 * t445 - t595 * t446 + t909;
t782 = rSges(5,1) * t876 + rSges(5,2) * t1089 + rSges(5,3) * t827;
t777 = qJDD(1) * t608 - t1027 * t620;
t380 = -pkin(5) * t960 + t487;
t773 = -t426 + t804;
t772 = -pkin(2) + (rSges(5,2) - pkin(3)) * t616;
t771 = -pkin(8) * t871 - t441;
t770 = t507 + t792;
t769 = t595 * t803;
t768 = t803 * t596;
t538 = rSges(2,1) * t617 - rSges(2,2) * t614;
t534 = rSges(2,1) * t614 + rSges(2,2) * t617;
t537 = -rSges(4,2) * t613 + t1023;
t753 = t595 * t74 + t596 * t73;
t752 = t595 * t73 - t596 * t74;
t751 = t595 * t76 + t596 * t75;
t750 = t595 * t75 - t596 * t76;
t749 = t595 * t78 + t596 * t77;
t748 = t595 * t77 - t596 * t78;
t747 = t595 * t80 + t596 * t79;
t746 = t595 * t79 - t596 * t80;
t745 = t595 * t83 + t596 * t82;
t744 = t595 * t82 - t596 * t83;
t743 = t595 * t87 + t596 * t86;
t742 = t595 * t86 - t596 * t87;
t721 = -t169 * t595 - t970;
t716 = t219 * t595 - t221 * t596;
t248 = -rSges(4,1) * t1089 - rSges(4,2) * t827 + t890;
t419 = t533 * t595;
t249 = -qJD(3) * t419 + (t537 * t596 + t581) * qJD(1);
t715 = t248 * t596 + t249 * t595;
t709 = t337 * t595 + t338 * t596;
t340 = rSges(5,2) * t957 + t802;
t700 = t340 + t775;
t697 = t804 - t906;
t694 = -t272 + t771;
t693 = -pkin(2) - t1081;
t689 = t595 * t313 + t596 * t314 - t446 * t876 + t840;
t683 = t873 * t885 + t529;
t681 = qJD(1) * (-pkin(2) * t877 + t570) + qJDD(1) * t462 + t777;
t680 = t65 * t437 - t1086;
t679 = t771 + t925;
t666 = t1077 * t648 - t208 * t439 + t210 * t440;
t665 = (Icges(7,5) * t355 - Icges(7,6) * t356) * t315 - (Icges(7,5) * t357 + Icges(7,6) * t358) * t316 + t427 * t489;
t664 = (Icges(6,5) * t395 - Icges(6,6) * t396) * t439 - (Icges(6,5) * t397 + Icges(6,6) * t398) * t440 + t451 * t1077;
t661 = -t613 * t917 + t616 * t919;
t660 = t613 * t914 + t616 * t916;
t658 = t616 * t665;
t657 = t616 * t664;
t646 = (t613 * t888 + t616 * t889) * qJD(1);
t645 = (-t613 * t886 + t616 * t887) * qJD(1);
t642 = qJDD(1) * t424 + t681 + t1090 * t595 + (t217 + t529) * qJD(1);
t623 = qJD(1) * t314 + qJDD(1) * t445 + (-t447 * t613 - t595 * t942) * pkin(8) - t447 * t531 - t441 * t575 + t642;
t14 = t489 * t103 + t1077 * t122 - t146 * t372 + t354 * t187 - t315 * t241 + t463 * t254 - t256 * t437 - t439 * t312 + t623;
t643 = t11 * t189 * t954 + t14 * t613 * t187 + t65 * (t613 * t103 + t1089 * t372 + t187 * t871);
t16 = -t101 * t358 + t158 * t181 + t159 * t184 + t178 * t669 + t357 * t99 + t957 * t97;
t17 = -t100 * t358 + t158 * t183 - t159 * t185 + t180 * t669 + t357 * t98 + t957 * t96;
t18 = t101 * t356 - t1089 * t178 + t160 * t181 + t161 * t184 + t355 * t99 + t954 * t97;
t19 = t100 * t356 - t1089 * t180 + t160 * t183 - t161 * t185 + t355 * t98 + t954 * t96;
t43 = t138 * t489 + t315 * t82 - t316 * t83;
t46 = -t158 * t649 - t159 * t651 + t222 * t957 + t223 * t357 - t224 * t358 - t647 * t669;
t47 = t1089 * t647 - t160 * t649 - t161 * t651 + t222 * t954 + t223 * t355 + t224 * t356;
t5 = t119 * t354 + t146 * t75 + t147 * t76 + t16 * t315 - t17 * t316 + t46 * t489;
t6 = t118 * t354 + t146 * t73 + t147 * t74 + t18 * t315 - t19 * t316 + t47 * t489;
t636 = ((-qJD(3) * t753 + t47) * t613 + (-qJD(1) * t752 + qJD(3) * t118 + t18 * t596 + t19 * t595) * t616) * t1053 + (t640 * t355 - t356 * t641 + t596 * t658) * t1054 + (t357 * t640 + t358 * t641 + t595 * t658) * t1051 + ((-qJD(3) * t751 + t46) * t613 + (-qJD(1) * t750 + qJD(3) * t119 + t16 * t596 + t17 * t595) * t616) * t1052 + (t665 * t613 + (t641 * t600 - t601 * t640) * t616) * t1042 + t5 * t823 + (t118 * t613 + t616 * t753) * t1058 + (t119 * t613 + t616 * t751) * t1057 + t6 * t822 + t43 * t811 + (t138 * t613 + t616 * t745) * t1050 + ((-qJD(3) * t745 + t61) * t613 + (-qJD(1) * t744 + qJD(3) * t138 + t20 * t596 + t21 * t595) * t616) * t1041 + (t1001 + t1006 + t1007 - t1015 + t1016) * t1036 + t1091 * t37 + t1092 * t36;
t72 = t219 * t440 + t221 * t439 + t663;
t635 = t72 * t716 + (-t595 * t84 + t596 * t85) * t426;
t626 = (t647 * t596 + t720) * t315 - (t647 * t595 + t719) * t316 + (t359 + t708) * t489;
t625 = (t648 * t596 + t718) * t439 - (t648 * t595 + t717) * t440 + (t399 + t707) * t1077;
t624 = t625 * t616;
t621 = (-t178 * t315 + t180 * t316 + t489 * t647) * t613 + t626 * t616;
t475 = t537 * qJD(3);
t394 = t827 - t835;
t392 = t882 * t872;
t391 = t596 * t949;
t390 = t595 * t949;
t311 = t372 * t957;
t305 = pkin(8) * t955 + t893;
t304 = pkin(8) * t959 + t894;
t299 = -rSges(6,3) * t955 + t895;
t298 = -rSges(6,3) * t959 + t896;
t297 = t652 * t596;
t296 = t652 * t595;
t295 = t650 * t596;
t294 = t650 * t595;
t280 = -rSges(7,3) * t955 + t897;
t279 = -rSges(7,3) * t959 + t898;
t278 = t651 * t596;
t277 = t651 * t595;
t276 = t649 * t596;
t275 = t649 * t595;
t266 = rSges(6,1) * t397 + rSges(6,2) * t398;
t265 = rSges(6,1) * t395 - rSges(6,2) * t396;
t251 = -rSges(5,3) * t835 + t782;
t250 = t759 * t575 + (t596 * t760 + t582) * qJD(1);
t167 = qJD(3) * t709 + qJD(2);
t126 = -t443 + (qJD(3) * t759 + t599) * t595 + t177 * qJD(1);
t125 = qJD(1) * t700 + t683;
t120 = (t339 * t596 - t340 * t595) * qJD(3) + t690;
t116 = -rSges(6,3) * t1089 + t930;
t94 = qJD(1) * t248 + qJDD(1) * t338 - t447 * t533 - t475 * t575 + t681;
t93 = -t857 - t475 * t873 + t448 * t533 + (-t249 - t442) * qJD(1) + t776 * qJDD(1);
t81 = qJD(3) * t715 + t337 * t447 - t338 * t448 + qJDD(2);
t63 = qJD(1) * t251 + qJDD(1) * t339 + t447 * t885 + t595 * t794 + t642;
t62 = -t448 * t759 + t596 * t794 + t700 * qJDD(1) + (-t250 + t691) * qJD(1) + t682;
t51 = t1089 * t648 - t203 * t650 - t204 * t652 + t269 * t954 + t270 * t395 + t271 * t396;
t50 = -t201 * t650 - t202 * t652 + t269 * t957 + t270 * t397 - t271 * t398 - t648 * t669;
t49 = -t340 * t447 + t911 * t448 + (t250 * t595 + t251 * t596) * qJD(3) + t662;
t45 = t1077 * t149 + t439 * t86 - t440 * t87;
t39 = t1077 * t116 + t219 * t463 - t256 * t426 - t272 * t439 + t623;
t27 = -t1089 * t210 + t109 * t954 + t111 * t395 + t113 * t396 + t203 * t213 - t204 * t215;
t26 = -t1089 * t208 + t110 * t954 + t112 * t395 + t114 * t396 + t203 * t211 + t204 * t214;
t25 = t109 * t957 + t111 * t397 - t113 * t398 + t201 * t213 - t202 * t215 + t210 * t669;
t24 = t110 * t957 + t112 * t397 - t114 * t398 + t201 * t211 + t202 * t214 + t208 * t669;
t22 = t115 * t439 + t116 * t440 - t219 * t257 + t221 * t256 + t627;
t15 = -t489 * t102 - t1077 * t121 + t147 * t372 - t354 * t189 - t316 * t241 - t463 * t255 + t257 * t437 - t440 * t312 + t622;
t10 = t1077 * t51 + t127 * t463 + t256 * t77 + t257 * t78 + t26 * t439 - t27 * t440;
t9 = t1077 * t50 + t128 * t463 + t24 * t439 - t25 * t440 + t256 * t79 + t257 * t80;
t1 = [(((t1114 * t596 + t1097 - t1113) * t596 + (t1114 * t595 + t1098 + t1160) * t595) * qJD(3) + t1125 - t1132) * t816 + ((t693 - t606) * t1121 + t1079 * t1008 * qJD(1) - (t1035 * t616 + t806) * t1030 + (t505 + t506 - t766 + (rSges(6,3) * t872 - qJ(4) * t871 - t599) * t595 + (-t608 + (-pkin(4) - pkin(7)) * t595) * qJD(1)) * t84 + (t84 + t571 + t930 + t789 * t872 + (t1079 * t595 - t1027) * qJD(1) + t1112 + t1127) * t85 + (t39 - g(2)) * (t219 + t1107) + (-g(1) + t38) * (t446 + t807 + t905)) * m(6) + (t46 + t36) * t1052 + (t50 + t40) * t1047 + (t1151 * qJD(3) + t1156 * t613 + t1157 * t616) * qJD(1) + (m(2) * (t534 ^ 2 + t538 ^ 2) + m(3) * (t434 ^ 2 + t456 * t435) + Icges(2,3) + Icges(3,3) - t1150) * qJDD(1) + ((t805 * t974 * t616 + t1135 * t15) * t595 - t1135 * t1030 + (-t1018 * t617 + t506 - t763 + (-t850 + ((t618 + t1034) * t616 + t688) * qJD(1)) * t596 + ((rSges(7,3) * qJD(3) - qJD(4) - t849 - t870) * t613 + (-pkin(7) - t591) * qJD(1)) * t595) * t64 + (t64 + t741 - t851 + t938 + t1034 * t828 + (-t850 + (-t577 + t693 - t605) * qJD(1)) * t595 + t1112 + t1128) * t65 + (t14 - g(2)) * (t692 + t1087 + t187) + (-g(1) + t15) * (t807 - t892 + t910)) * m(7) + (t1117 + t1095) * t1045 + (t237 + t1096) * t1044 + ((t1113 * t595 + ((t1148 + t1166) * t596 + t1099 + t1159 + t1165) * t596) * qJD(3) + t1133) * t813 + ((-t455 * t620 - g(2) + t777) * t435 + (-g(1) - t857 + (-0.2e1 * t456 - t608 + t435) * t620) * t434) * m(3) + (t1119 + t1122) * t815 + (t1118 - t1123 + t1124) * t814 - m(2) * (-g(1) * t534 + g(2) * t538) + t118 * t1058 - t448 * t240 / 0.2e1 + (-(-t832 - t168 - t450 + (-t337 - t1027) * qJD(1)) * t169 + t169 * (t570 + t890) + (t533 * t971 - t969) * qJD(3) + ((-t168 * t617 - t169 * t614) * pkin(1) + (-pkin(2) - t537) * t970 + (t168 * (-rSges(4,3) - pkin(7)) + t169 * t821) * t595) * qJD(1) + (t94 - g(2)) * t268 + (t93 - g(1)) * (t595 * t821 + t807 + t883)) * m(4) - t1015 / 0.2e1 + t1016 / 0.2e1 + t51 * t1048 + t47 * t1053 + t128 * t1055 + t127 * t1056 + t119 * t1057 + t1000 + t1001 + t36 * t1051 + t40 * t1046 - t1013 / 0.2e1 + t1014 / 0.2e1 + t1006 / 0.2e1 + t1007 / 0.2e1 + t1004 / 0.2e1 + t1005 / 0.2e1 + ((t506 + (-t599 + (t1002 * t616 - t1020) * qJD(3)) * t595 + (-t608 + (t1002 * t613 + t772) * t596 + (-rSges(5,1) - pkin(7)) * t595) * qJD(1)) * t125 + (t63 - g(2)) * t177 + (t62 - g(1)) * ((t772 - t602) * t595 + t802 + t807) + (-pkin(3) * t828 + t125 - t683 + t782 + (-t340 + (t693 - t604) * t595) * qJD(1) + t1129) * t126) * m(5); m(3) * qJDD(2) + (-m(3) - m(4) + t864) * g(3) + m(4) * t81 + m(5) * t49 + m(6) * t22 + m(7) * t11; (t40 * t596 + t41 * t595) * t598 / 0.2e1 + (t138 * t508 - t83 * t390 - t82 * t391 + ((-t276 * t601 - t278 * t600 + t178) * t315 - (-t275 * t601 - t277 * t600 + t180) * t316 + (-t361 * t601 - t363 * t600 - t647) * t489) * t616 + t626 * t613) * t1042 + (((-t295 * t615 - t297 * t612 + t208) * t439 - (-t294 * t615 - t296 * t612 + t210) * t440 + (-t401 * t615 - t403 * t612 - t648) * t1077 + t149 * qJD(5)) * t616 + (-qJD(5) * t743 + t625) * t613) * t1040 + t1136 * t1044 + t1137 * t1045 + (g(1) * t423 + g(2) * t419 - g(3) * t537 - (t168 * t419 - t969) * qJD(1) - (t167 * (-t419 * t595 - t423 * t596) + t721 * t537) * qJD(3) + t81 * t709 + t167 * ((t337 * t596 - t338 * t595) * qJD(1) + t715) + t721 * t475 + (-t94 * t595 - t93 * t596 + (-t169 * t596 + t971) * qJD(1)) * t533) * m(4) + ((t295 * t397 - t297 * t398) * t439 - (t294 * t397 - t296 * t398) * t440 + (t397 * t401 - t398 * t403) * t1077 + (t128 * t616 - t79 * t955) * qJD(5) + ((-qJD(5) * t80 + t666) * t613 + t624) * t595) * t1046 + ((t295 * t395 + t297 * t396) * t439 - (t294 * t395 + t296 * t396) * t440 + (t395 * t401 + t396 * t403) * t1077 + (t127 * t616 - t78 * t959) * qJD(5) + ((-qJD(5) * t77 + t666) * t613 + t624) * t596) * t1049 + (qJD(1) * t1119 + t1094 * qJD(3) + qJDD(1) * t1117 + t1097 * t447 + t1098 * t448 + t10 + t6) * t595 / 0.2e1 - (t1118 * qJD(1) + t1093 * qJD(3) + t1149 * qJDD(1) + t1099 * t447 + t1100 * t448 + t5 + t9) * t596 / 0.2e1 + ((t1099 * t596 + t1100 * t595) * qJD(1) + t1093) * t814 + (t1095 * t595 - t1096 * t596) * qJDD(1) / 0.2e1 + ((t1097 * t596 + t1098 * t595) * qJD(1) + t1094) * t815 + (-g(1) * (t541 + t895) - g(2) * (t539 + t896) - g(3) * (t425 - t803) + t22 * t783 + (t22 * t219 + t38 * t773) * t596 + (t22 * t221 + t39 * t773) * t595 + (-qJD(3) * t769 - t1077 * t299 - t219 * t868 + t425 * t439 + t595 * t694 + t773 * t876 - t904) * t85 + (-qJD(3) * t768 + t1077 * t298 + t221 * t868 + t425 * t440 + t426 * t877 + t596 * t694 - t770 + t899) * t84 + (t689 + (qJD(1) * t221 + t116) * t596 + (t115 + (-t219 + t901) * qJD(1)) * t595 - t298 * t439 - t299 * t440 - t839) * t72 + (-g(1) * t789 - t1028 * t859 - (-t72 * t793 - t841) * pkin(8) - t635 * qJD(5)) * t613) * m(6) - ((((-t913 - t918) * t596 + (-t914 + t917) * t595) * t616 + ((t915 + t920) * t596 + (t916 + t919) * t595) * t613) * qJD(3) + ((t886 - t888) * t616 + (t887 + t889) * t613) * qJD(1)) * qJD(1) / 0.2e1 + ((-t873 * t962 - t879) * t596 + (t645 + (t661 * t595 + (-t1062 + t961) * t596) * qJD(3)) * t595 + (t411 * t873 + t878) * t596 + (t646 + (t660 * t595 + (-t1063 - t412) * t596) * qJD(3)) * t595) * t813 + ((-t575 * t961 + t879) * t595 + (t645 + (-t1062 * t596 + (t962 + t661) * t595) * qJD(3)) * t596 + (t412 * t575 - t878) * t595 + (t646 + (-t1063 * t596 + (-t411 + t660) * t595) * qJD(3)) * t596) * t816 + (-g(1) * (t421 + t422) - g(2) * (t417 + t418) + g(3) * t884 - t125 * (-qJD(1) * t418 + t792) - t126 * (qJD(1) * t422 + t904) - t120 * t839 - ((t120 * t422 + t125 * t884) * t596 + (t120 * t418 + t126 * t884) * t595) * qJD(3) + t125 * t444 + t49 * t909 + t120 * t840 + (t62 * t885 + t125 * t900 + t49 * t339 + t120 * t251 + (-t120 * t340 + t126 * t885) * qJD(1)) * t596 + (t63 * t885 + t126 * t900 - t49 * t340 + t120 * t250 + (t120 * t911 - t125 * t759) * qJD(1)) * t595) * m(5) + (t1123 * t596 + t1122 * t595 + (t1095 * t596 + t1096 * t595) * qJD(1)) * qJD(1) / 0.2e1 + (t40 + t36 + t1124) * t817 + (t41 + t37 + t1125) * t877 / 0.2e1 + (qJD(1) * t743 + t28 * t595 - t29 * t596) * t1039 + (qJD(1) * t745 + t20 * t595 - t21 * t596) * t1041 + t742 * t1043 + (qJD(1) * t747 + t24 * t595 - t25 * t596) * t1047 - t508 * t43 / 0.2e1 + t752 * t1058 + (qJD(1) * t749 + t26 * t595 - t27 * t596) * t1048 + t744 * t1050 + (qJD(1) * t751 + t16 * t595 - t17 * t596) * t1052 + (qJD(1) * t753 + t18 * t595 - t19 * t596) * t1053 + t746 * t1055 + t748 * t1056 + t750 * t1057 + t390 * t37 / 0.2e1 + t391 * t36 / 0.2e1 + (-g(1) * (t541 + t893 + t897) - g(2) * (t539 + t894 + t898) - g(3) * (t371 + t1081 + t577 - t943) - t1034 * t1066 - t64 * (-t1077 * t304 - t1082 * t440 - t189 * t508 - t255 * t868 - t279 * t489 - t316 * t371 - t372 * t390 + t770) - t65 * (t1077 * t305 - t1082 * t439 + t187 * t508 + t254 * t868 + t280 * t489 - t315 * t371 + t372 * t391 + t904) - t48 * (t187 * t390 - t189 * t391 + t279 * t315 + t280 * t316 + t304 * t439 + t305 * t440 + t839) - (t64 * t768 + t65 * t769) * qJD(3) - ((-t48 * t793 - t842) * pkin(8) + (t48 * (t254 * t595 - t255 * t596) + (-t595 * t64 + t596 * t65) * t437) * qJD(5)) * t613 + t64 * t899 + t11 * t783 + t48 * t689 + (t15 * t697 + t64 * t679 + (t65 * t697 + t1086) * qJD(1) + t1064) * t596 + (t14 * t697 + t65 * t679 + t11 * t932 + t48 * t940 + (t1085 + t48 * (t901 - t933)) * qJD(1)) * t595) * m(7) + ((t276 * t357 - t278 * t358) * t315 - t75 * t391 - (t275 * t357 - t277 * t358) * t316 - t76 * t390 + (t357 * t361 - t358 * t363) * t489 + t119 * t508 + t621 * t595) * t1051 + ((t276 * t355 + t278 * t356) * t315 - t73 * t391 - (t275 * t355 + t277 * t356) * t316 - t74 * t390 + (t355 * t361 + t356 * t363) * t489 + t118 * t508 + t621 * t596) * t1054 - t45 * t868 / 0.2e1; t864 * (-g(3) * t616 + t1066) - m(5) * (t120 * t392 + t125 * t394 + t126 * t393) - m(6) * (t392 * t72 + t393 * t85 + t394 * t84) - m(7) * (t392 * t48 + t393 * t65 + t394 * t64) + 0.2e1 * ((t125 * t873 + t126 * t575 - t49) * t1061 + (t575 * t85 + t84 * t873 - t22) * t1060 + (t575 * t65 + t64 * t873 - t11) * t1059) * t616 + 0.2e1 * ((qJD(3) * t120 - t125 * t877 + t126 * t876 + t595 * t63 + t596 * t62) * t1061 + (qJD(3) * t72 + t38 * t596 + t39 * t595 - t84 * t877 + t841) * t1060 + (qJD(3) * t48 + t14 * t595 + t15 * t596 - t64 * t877 + t842) * t1059) * t613; t636 + ((-qJD(3) * t743 + t71) * t613 + (-qJD(1) * t742 + qJD(3) * t149 + t28 * t596 + t29 * t595) * t616) * t1039 + (t149 * t613 + t616 * t743) * t1043 + (t128 * t613 + t616 * t747) * t1055 + (t127 * t613 + t616 * t749) * t1056 + ((-qJD(3) * t747 + t50) * t613 + (-qJD(1) * t746 + qJD(3) * t128 + t24 * t596 + t25 * t595) * t616) * t1047 + ((-qJD(3) * t749 + t51) * t613 + (-qJD(1) * t748 + qJD(3) * t127 + t26 * t596 + t27 * t595) * t616) * t1048 + (t664 * t613 + (t639 * t612 - t615 * t638) * t616) * t1040 + (t397 * t638 + t398 * t639 + t595 * t657) * t1046 + (t1000 + t1004 + t1005 - t1013 + t1014) * t1036 + t10 * t822 + t9 * t823 + t45 * t811 + (t638 * t395 - t396 * t639 + t596 * t657) * t1049 + t1091 * t41 + t1092 * t40 + (-t64 * (-t1077 * t381 + t440 * t858 + t799) - t65 * (t1077 * t380 + t439 * t858 + t800) - t48 * (t380 * t440 + t381 * t439 + t801) + t15 * t311 + t64 * t931 + t48 * t999 + (-t15 * t932 - t64 * t940 + t14 * t254 + t65 * t122 + (t680 * t596 + (t48 * t254 - t1085) * t595) * qJD(3)) * t613 + ((t65 * t254 - t64 * t932) * qJD(3) + (-t14 * t906 + t65 * t925 + t11 * t255 + t48 * t121 + (t64 * t437 - t48 * t933) * qJD(1)) * t596 + (qJD(1) * t680 + t15 * t437 + t64 * t312 - t1064) * t595) * t616 + t643 - g(1) * (t380 + t252) - g(2) * (t381 + t253) - g(3) * (t565 + (-t1022 - t1033) * t616)) * m(7) + (-t84 * (-t1077 * t266 - t440 * t454) - t85 * (t1077 * t265 - t439 * t454) - t72 * (t265 * t440 + t266 * t439) + (qJD(3) * t635 - t84 * t115 + t85 * t116 + t39 * t219 - t38 * t221) * t613 + (t84 * (-qJD(3) * t221 + t272 * t595) + t85 * (qJD(3) * t219 - t272 * t596) - t22 * t716 + t72 * (t115 * t596 - t116 * t595 - t219 * t876 - t221 * t877) + (t1121 - t39 * t596 + (t595 * t85 + t1008) * qJD(1)) * t426) * t616 - g(1) * t265 - g(2) * t266 - g(3) * t454) * m(6); t636 + (t15 * (-t189 * t613 + t311) + t48 * (-t189 * t828 + t999) + (-t189 * t974 + (-qJD(1) * t187 * t48 - t14 * t372 - t241 * t65) * t596 + (-t11 * t187 + t48 * (-qJD(1) * t189 - t103)) * t595) * t616 + t643 - t48 * t801 - t65 * t800 - g(1) * t252 - g(2) * t253 - g(3) * t430 + (-t102 * t613 - t372 * t831 - t799 + t931) * t64) * m(7);];
tau  = t1;
