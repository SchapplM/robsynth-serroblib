% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:22
% EndTime: 2019-03-09 04:48:26
% DurationCPUTime: 113.63s
% Computational Cost: add. (34503->1485), mult. (60159->1802), div. (0->0), fcn. (56581->8), ass. (0->669)
t1209 = Icges(6,1) + Icges(7,1);
t1190 = Icges(7,4) + Icges(6,5);
t1189 = Icges(6,6) - Icges(7,6);
t630 = qJ(4) + pkin(9);
t606 = sin(t630);
t1240 = (Icges(6,4) - Icges(7,5)) * t606;
t607 = cos(t630);
t1239 = t1209 * t607 - t1240;
t1208 = Icges(6,2) + Icges(7,3);
t633 = sin(qJ(3));
t636 = cos(qJ(3));
t1004 = Icges(6,4) * t607;
t746 = -Icges(6,2) * t606 + t1004;
t1238 = -t1189 * t633 - t636 * t746;
t637 = cos(qJ(1));
t972 = t633 * t637;
t634 = sin(qJ(1));
t978 = t607 * t634;
t424 = t606 * t972 + t978;
t396 = Icges(7,5) * t424;
t965 = t637 * t607;
t971 = t634 * t606;
t425 = t633 * t965 - t971;
t966 = t636 * t637;
t210 = -Icges(7,1) * t425 + Icges(7,4) * t966 - t396;
t399 = Icges(6,4) * t424;
t212 = Icges(6,1) * t425 - Icges(6,5) * t966 - t399;
t1121 = -t212 + t210;
t397 = Icges(7,5) * t425;
t197 = -Icges(7,6) * t966 + Icges(7,3) * t424 + t397;
t400 = Icges(6,4) * t425;
t207 = Icges(6,2) * t424 + Icges(6,6) * t966 - t400;
t1122 = -t207 - t197;
t635 = cos(qJ(4));
t967 = t635 * t637;
t632 = sin(qJ(4));
t970 = t634 * t632;
t478 = t633 * t967 - t970;
t463 = Icges(5,4) * t478;
t969 = t634 * t635;
t477 = t632 * t972 + t969;
t283 = Icges(5,2) * t477 + Icges(5,6) * t966 - t463;
t462 = Icges(5,4) * t477;
t285 = Icges(5,1) * t478 - Icges(5,5) * t966 - t462;
t1237 = t1121 * t425 + t1122 * t424 - t283 * t477 - t285 * t478;
t1236 = Icges(7,2) + Icges(5,3) + Icges(6,3);
t1235 = Icges(5,5) * t635 - Icges(5,6) * t632 - t1189 * t606 + t1190 * t607;
t977 = t607 * t636;
t553 = Icges(7,5) * t977;
t979 = t606 * t636;
t1207 = Icges(7,3) * t979 + t1238 + t553;
t1202 = t1190 * t633 + t1239 * t636;
t1230 = t1235 * t636 + t1236 * t633;
t1007 = Icges(5,4) * t635;
t747 = -Icges(5,2) * t632 + t1007;
t431 = Icges(5,6) * t633 + t636 * t747;
t1008 = Icges(5,4) * t632;
t751 = Icges(5,1) * t635 - t1008;
t435 = Icges(5,5) * t633 + t636 * t751;
t1098 = -t1202 * t425 - t1207 * t424 + t1230 * t966 + t431 * t477 - t435 * t478;
t422 = t633 * t971 - t965;
t973 = t633 * t634;
t423 = t606 * t637 + t607 * t973;
t968 = t634 * t636;
t199 = Icges(6,5) * t423 - Icges(6,6) * t422 - Icges(6,3) * t968;
t202 = Icges(7,4) * t423 - Icges(7,2) * t968 + Icges(7,6) * t422;
t475 = -t633 * t970 + t967;
t974 = t632 * t637;
t476 = t633 * t969 + t974;
t278 = Icges(5,5) * t476 + Icges(5,6) * t475 - Icges(5,3) * t968;
t1212 = t199 + t202 + t278;
t999 = Icges(7,5) * t422;
t208 = Icges(7,1) * t423 - Icges(7,4) * t968 + t999;
t398 = Icges(6,4) * t422;
t211 = Icges(6,1) * t423 - Icges(6,5) * t968 - t398;
t1219 = t208 + t211;
t395 = Icges(7,5) * t423;
t196 = -Icges(7,6) * t968 + Icges(7,3) * t422 + t395;
t1006 = Icges(6,4) * t423;
t205 = -Icges(6,2) * t422 - Icges(6,6) * t968 + t1006;
t1220 = t196 - t205;
t1009 = Icges(5,4) * t476;
t281 = Icges(5,2) * t475 - Icges(5,6) * t968 + t1009;
t461 = Icges(5,4) * t475;
t284 = Icges(5,1) * t476 - Icges(5,5) * t968 + t461;
t1139 = t1212 * t966 - t1219 * t425 - t1220 * t424 + t477 * t281 - t478 * t284;
t886 = qJD(4) * t636;
t888 = qJD(3) * t637;
t508 = -t634 * t886 + t888;
t887 = qJD(4) * t633;
t588 = qJD(1) + t887;
t1234 = -t1098 * t588 - t1139 * t508;
t201 = -Icges(6,5) * t425 + Icges(6,6) * t424 + Icges(6,3) * t966;
t204 = -Icges(7,4) * t425 + Icges(7,2) * t966 - Icges(7,6) * t424;
t280 = -Icges(5,5) * t478 + Icges(5,6) * t477 + Icges(5,3) * t966;
t1096 = t201 + t280 + t204;
t1203 = t1219 * t423 + t1220 * t422 + t475 * t281 + t476 * t284;
t1233 = (-t1096 * t637 - t1212 * t634) * t636 + t1203 + t1237;
t997 = Icges(7,5) * t607;
t741 = Icges(7,3) * t606 + t997;
t1168 = -t1189 * t636 + (-t741 + t746) * t633;
t1167 = t1190 * t636 - t1239 * t633;
t1232 = (t1208 * t607 + t1240) * t636;
t1229 = -t1235 * t633 + t1236 * t636;
t1228 = (-Icges(5,5) * t632 - Icges(5,6) * t635 - t1189 * t607 - t1190 * t606) * t636;
t1198 = -t1096 * t968 + t1121 * t423 + t1122 * t422 + t475 * t283 - t285 * t476;
t1210 = -t1202 * t607 - t1207 * t606 + t431 * t632 - t435 * t635;
t838 = t633 * t888;
t893 = qJD(1) * t636;
t1162 = t634 * t893 + t838;
t809 = qJD(1) * t633 + qJD(4);
t835 = t636 * t888;
t1165 = -t809 * t634 + t835;
t885 = qJD(4) * t637;
t833 = t633 * t885;
t892 = qJD(1) * t637;
t187 = (-t833 - t892) * t607 - t1165 * t606;
t714 = t588 * t637;
t188 = -t1165 * t607 + t606 * t714;
t100 = Icges(7,5) * t188 - Icges(7,6) * t1162 + Icges(7,3) * t187;
t106 = Icges(6,4) * t188 - Icges(6,2) * t187 - Icges(6,6) * t1162;
t1224 = t100 - t106;
t712 = t809 * t637;
t889 = qJD(3) * t636;
t836 = t634 * t889;
t1164 = -t836 - t712;
t189 = -t1164 * t606 + t588 * t978;
t716 = t634 * t588;
t190 = -t1164 * t607 - t606 * t716;
t890 = qJD(3) * t634;
t839 = t633 * t890;
t840 = t636 * t892;
t685 = t839 - t840;
t101 = Icges(7,5) * t190 + Icges(7,6) * t685 + Icges(7,3) * t189;
t107 = Icges(6,4) * t190 - Icges(6,2) * t189 + Icges(6,6) * t685;
t1223 = t101 - t107;
t108 = Icges(7,1) * t188 - Icges(7,4) * t1162 + Icges(7,5) * t187;
t110 = Icges(6,1) * t188 - Icges(6,4) * t187 - Icges(6,5) * t1162;
t1222 = t108 + t110;
t109 = Icges(7,1) * t190 + Icges(7,4) * t685 + Icges(7,5) * t189;
t111 = Icges(6,1) * t190 - Icges(6,4) * t189 + Icges(6,5) * t685;
t1221 = t109 + t111;
t1218 = qJD(3) * t1168 + t1232 * qJD(4);
t443 = (-Icges(6,1) * t606 - t1004) * t636;
t1217 = (-Icges(7,1) * t606 + t997) * t886 + qJD(4) * t443 + t1167 * qJD(3);
t102 = Icges(6,5) * t188 - Icges(6,6) * t187 - Icges(6,3) * t1162;
t104 = Icges(7,4) * t188 - Icges(7,2) * t1162 + Icges(7,6) * t187;
t247 = t1165 * t632 + t635 * t714;
t248 = -t1165 * t635 + t632 * t714;
t131 = Icges(5,5) * t248 + Icges(5,6) * t247 - Icges(5,3) * t1162;
t1214 = t102 + t104 + t131;
t103 = Icges(6,5) * t190 - Icges(6,6) * t189 + Icges(6,3) * t685;
t105 = Icges(7,4) * t190 + Icges(7,2) * t685 + Icges(7,6) * t189;
t249 = t1164 * t632 - t635 * t716;
t715 = t588 * t632;
t250 = t635 * t712 + (t635 * t889 - t715) * t634;
t132 = Icges(5,5) * t250 + Icges(5,6) * t249 + Icges(5,3) * t685;
t1213 = t103 + t105 + t132;
t1211 = qJD(3) * t1229 + qJD(4) * t1228;
t832 = t636 * t885;
t507 = t832 + t890;
t729 = t283 * t632 + t285 * t635;
t731 = t281 * t632 - t284 * t635;
t732 = t207 * t606 + t212 * t607;
t734 = t205 * t606 - t211 * t607;
t735 = -t197 * t606 + t210 * t607;
t737 = t196 * t606 + t208 * t607;
t1205 = (t1210 + t1229) * t588 + (t1230 * t634 + t731 + t734 - t737) * t508 + (-t1230 * t637 + t729 + t732 - t735) * t507;
t1099 = t1202 * t423 + t1207 * t422 - t1230 * t968 + t431 * t475 + t435 * t476;
t1174 = rSges(7,1) + pkin(5);
t1200 = t636 * t741 + t1238;
t1109 = -t424 * rSges(7,3) - t1174 * t425;
t1120 = -rSges(7,2) * t966 - t1109;
t990 = qJ(6) * t424;
t947 = t990 + t1120;
t134 = Icges(5,4) * t250 + Icges(5,2) * t249 + Icges(5,6) * t685;
t136 = Icges(5,1) * t250 + Icges(5,4) * t249 + Icges(5,5) * t685;
t1148 = -t1162 * t1212 + t1213 * t966 + t1219 * t188 + t1220 * t187 - t1221 * t425 - t1223 * t424 + t134 * t477 - t136 * t478 + t247 * t281 + t248 * t284;
t133 = Icges(5,4) * t248 + Icges(5,2) * t247 - Icges(5,6) * t1162;
t135 = Icges(5,1) * t248 + Icges(5,4) * t247 - Icges(5,5) * t1162;
t1147 = -t1096 * t1162 + t1121 * t188 + t1122 * t187 + t1214 * t966 - t1222 * t425 - t1224 * t424 + t133 * t477 - t135 * t478 + t247 * t283 - t248 * t285;
t1146 = t1212 * t685 - t1213 * t968 + t1219 * t190 + t1220 * t189 + t1221 * t423 + t1223 * t422 + t134 * t475 + t136 * t476 + t249 * t281 + t250 * t284;
t1145 = t1096 * t685 + t1121 * t190 + t1122 * t189 - t1214 * t968 + t1222 * t423 + t1224 * t422 + t133 * t475 + t135 * t476 + t249 * t283 - t250 * t285;
t430 = Icges(5,6) * t636 - t633 * t747;
t483 = (-Icges(5,2) * t635 - t1008) * t636;
t295 = qJD(3) * t430 + qJD(4) * t483;
t434 = Icges(5,5) * t636 - t633 * t751;
t486 = (-Icges(5,1) * t632 - t1007) * t636;
t298 = qJD(3) * t434 + qJD(4) * t486;
t1101 = -t1162 * t1230 + t1202 * t188 + t1207 * t187 + t1211 * t966 - t1217 * t425 - t1218 * t424 + t247 * t431 + t248 * t435 + t295 * t477 - t298 * t478;
t1100 = t1202 * t190 + t1207 * t189 - t1211 * t968 + t1217 * t423 + t1218 * t422 + t1230 * t685 + t249 * t431 + t250 * t435 + t295 * t475 + t298 * t476;
t1141 = -t1212 * t968 + t1203;
t1197 = t1096 * t966 - t1237;
t79 = t202 * t633 + t636 * t737;
t81 = t199 * t633 - t636 * t734;
t90 = t278 * t633 - t636 * t731;
t1196 = t79 + t81 + t90;
t80 = t204 * t633 + t636 * t735;
t82 = t201 * t633 - t636 * t732;
t91 = t280 * t633 - t636 * t729;
t1195 = t80 + t82 + t91;
t1097 = -t1210 * t636 + t1230 * t633;
t1116 = t478 * rSges(5,1) - t477 * rSges(5,2);
t290 = -rSges(5,3) * t966 + t1116;
t1040 = rSges(5,3) * t633;
t1042 = rSges(5,2) * t632;
t1045 = rSges(5,1) * t635;
t783 = -t1042 + t1045;
t452 = t636 * t783 + t1040;
t1193 = t290 * t588 + t452 * t507;
t1192 = (t1209 * t422 + t1006 - t1220 - t395) * t508 + (-t1209 * t424 - t1122 + t397 - t400) * t507 + (Icges(7,1) * t979 - t1207 - t443 - t553) * t588;
t1191 = (-t1208 * t423 + t1219 - t398 + t999) * t508 + (t1208 * t425 + t1121 - t396 + t399) * t507 + (t1202 - t1232) * t588;
t1173 = rSges(7,3) + qJ(6);
t590 = pkin(3) * t973;
t845 = pkin(3) * t835 + pkin(8) * t1162;
t323 = qJD(1) * t590 - t845;
t493 = pkin(3) * t968 + pkin(8) * t973;
t1054 = pkin(8) * t633;
t557 = pkin(3) * t636 + t1054;
t495 = t557 * t637;
t1188 = t637 * t323 + t493 * t890 + t495 * t888;
t1181 = (-t295 * t632 + t298 * t635 + t1217 * t607 + t1218 * t606 + (-t1202 * t606 + t1207 * t607 - t431 * t635 - t435 * t632) * qJD(4) + t1230 * qJD(3)) * t636 + (qJD(3) * t1210 + t1211) * t633;
t1178 = t1099 * t588 + t1198 * t507;
t1177 = t1205 * t636;
t1172 = t1200 * t634;
t1171 = t1200 * t637;
t1170 = t1202 * t634;
t1169 = t1202 * t637;
t631 = -qJ(5) - pkin(8);
t1046 = -pkin(8) - t631;
t612 = qJD(5) * t633;
t602 = pkin(4) * t635 + pkin(3);
t1047 = pkin(3) - t602;
t825 = t1047 * t633;
t1039 = pkin(4) * qJD(4);
t869 = t632 * t1039;
t276 = -t636 * t869 + t612 + (t1046 * t636 + t825) * qJD(3);
t1108 = t1047 * t636;
t373 = t1046 * t633 - t1108;
t352 = t373 * t892;
t976 = t631 * t636;
t1160 = t976 - t825;
t627 = t636 * pkin(8);
t372 = -t627 - t1160;
t579 = t634 * t612;
t1056 = pkin(3) * t633;
t556 = t627 - t1056;
t811 = qJD(1) * t495 + t556 * t890;
t519 = qJD(3) * t556;
t910 = t634 * t519 + t557 * t892;
t1166 = t634 * t276 - t507 * t372 + t352 - t579 - t811 + t910;
t1163 = t633 * t892 + t836;
t781 = -t425 * rSges(6,1) + t424 * rSges(6,2);
t217 = -rSges(6,3) * t966 - t781;
t504 = t602 * t835;
t884 = qJD(5) * t637;
t580 = t636 * t884;
t593 = pkin(4) * t974;
t129 = t631 * t838 - t504 + t580 + t477 * t1039 + (t1160 * t634 + t593) * qJD(1) + t845;
t492 = -pkin(8) * t968 + t590;
t849 = t602 * t973 + t631 * t968 + t593;
t287 = -t492 + t849;
t595 = pkin(8) * t966;
t494 = pkin(3) * t972 - t595;
t908 = t602 * t972 + t631 * t966;
t288 = pkin(4) * t970 + t494 - t908;
t561 = t631 * t972;
t320 = t561 + (t1054 + t1108) * t637;
t613 = qJD(5) * t636;
t834 = t634 * t887;
t1161 = t637 * t129 - t287 * t833 - t288 * t834 - t508 * t320 + t1188 - t613;
t1159 = qJD(1) * t452;
t548 = t637 * pkin(1) + t634 * qJ(2);
t628 = t637 * pkin(7);
t1112 = t628 + t548;
t626 = t637 * rSges(4,3);
t451 = rSges(4,1) * t973 + rSges(4,2) * t968 + t626;
t1157 = t451 + t1112;
t528 = t602 * t968;
t319 = -t631 * t973 - t493 + t528;
t894 = qJD(1) * t634;
t500 = t557 * t894;
t1155 = -t287 * t886 - t588 * t319 + t373 * t894 + t500;
t880 = qJD(1) * qJD(3);
t524 = qJDD(3) * t634 + t637 * t880;
t879 = qJDD(4) * t636;
t310 = -qJD(4) * t1162 + t637 * t879 + t524;
t610 = qJDD(3) * t637;
t311 = -qJD(1) * t832 + t610 + (-t879 + (-qJD(1) + t887) * qJD(3)) * t634;
t502 = qJD(3) * t886 + qJDD(4) * t633 + qJDD(1);
t1154 = t1098 * t502 + t1101 * t588 + t1139 * t311 + t1147 * t507 + t1148 * t508 + t1197 * t310;
t1153 = t1099 * t502 + t1100 * t588 + t1141 * t311 + t1145 * t507 + t1146 * t508 + t1198 * t310;
t1010 = Icges(4,4) * t636;
t536 = -Icges(4,2) * t633 + t1010;
t752 = Icges(4,1) * t633 + t1010;
t906 = t536 + t752;
t1011 = Icges(4,4) * t633;
t538 = Icges(4,1) * t636 - t1011;
t748 = Icges(4,2) * t636 + t1011;
t907 = -t748 + t538;
t1152 = (t633 * t906 - t636 * t907) * qJD(1);
t393 = t424 * qJD(6);
t949 = -rSges(7,2) * t968 + t1173 * t422 + t1174 * t423;
t858 = t287 + t949;
t1027 = t633 * rSges(7,2);
t776 = rSges(7,1) * t607 + rSges(7,3) * t606;
t388 = t636 * t776 + t1027;
t775 = pkin(5) * t607 + qJ(6) * t606;
t920 = t775 * t636 + t388;
t1151 = -t508 * t920 + t858 * t588 - t393;
t1002 = Icges(4,5) * t637;
t581 = Icges(4,4) * t968;
t436 = Icges(4,1) * t973 + t1002 + t581;
t437 = -Icges(4,5) * t634 + t637 * t752;
t485 = t536 * t637;
t670 = t634 * (t437 + t485) - t637 * (-Icges(4,2) * t973 + t436 + t581);
t432 = Icges(4,6) * t637 + t634 * t748;
t433 = -Icges(4,6) * t634 + t637 * t748;
t487 = t538 * t634;
t488 = t538 * t637;
t671 = t634 * (t433 - t488) - t637 * (t432 - t487);
t1150 = -t671 * t633 + t670 * t636;
t21 = (-qJD(3) * t737 + t105) * t633 + (qJD(3) * t202 + t101 * t606 + t109 * t607 + (t196 * t607 - t208 * t606) * qJD(4)) * t636;
t23 = (qJD(3) * t734 + t103) * t633 + (qJD(3) * t199 - t107 * t606 + t111 * t607 + (-t205 * t607 - t211 * t606) * qJD(4)) * t636;
t29 = (qJD(3) * t731 + t132) * t633 + (qJD(3) * t278 - t134 * t632 + t136 * t635 + (-t281 * t635 - t284 * t632) * qJD(4)) * t636;
t1144 = t21 + t23 + t29;
t22 = (-qJD(3) * t735 + t104) * t633 + (qJD(3) * t204 + t100 * t606 + t108 * t607 + (-t197 * t607 - t210 * t606) * qJD(4)) * t636;
t24 = (qJD(3) * t732 + t102) * t633 + (qJD(3) * t201 - t106 * t606 + t110 * t607 + (-t207 * t607 + t212 * t606) * qJD(4)) * t636;
t30 = (qJD(3) * t729 + t131) * t633 + (qJD(3) * t280 - t133 * t632 + t135 * t635 + (-t283 * t635 + t285 * t632) * qJD(4)) * t636;
t1143 = t22 + t24 + t30;
t1103 = t1141 * t508 + t1178;
t1102 = t1197 * t507 - t1234;
t1142 = t1097 * t588 + t1195 * t507 + t1196 * t508;
t1137 = -t1096 * t507 - t1212 * t508 - t1230 * t588;
t1094 = -t1195 * t637 + t1196 * t634;
t1093 = t1139 * t634 - t1197 * t637;
t1092 = t1141 * t634 - t1198 * t637;
t1136 = t1195 * t634 + t1196 * t637;
t1135 = t1139 * t637 + t1197 * t634;
t1134 = t1141 * t637 + t1198 * t634;
t394 = qJD(6) * t422;
t498 = t557 * t890;
t614 = qJD(2) * t634;
t1050 = t634 * pkin(7);
t617 = t637 * qJ(2);
t544 = t634 * pkin(1) - t617;
t815 = -t544 - t1050;
t794 = t494 + t815;
t666 = qJD(1) * t794 + t498 + t614;
t831 = t634 * t613;
t648 = t507 * t373 + t666 - t831;
t856 = -t288 + t947;
t50 = t507 * t920 + t588 * t856 + t394 + t648;
t308 = t508 * t373;
t1086 = (t492 + t1112) * qJD(1) - t557 * t888;
t615 = qJD(2) * t637;
t665 = t1086 - t615;
t657 = t580 - t308 + t665;
t51 = t657 + t1151;
t638 = qJD(1) ^ 2;
t881 = qJD(1) * qJD(2);
t902 = qJDD(2) * t634 + t637 * t881;
t717 = -t628 * t638 + t902;
t647 = qJDD(1) * t794 + t519 * t890 + t524 * t557 + t717;
t878 = qJDD(5) * t636;
t464 = qJD(1) * t548 - t615;
t932 = -t323 - t464;
t640 = t310 * t373 + t507 * t276 + qJD(3) * t579 + (-t580 + t932) * qJD(1) - t634 * t878 + t647;
t1130 = -t1173 * t187 - t1174 * t188 + t393;
t1016 = -rSges(7,2) * t1162 - t1130;
t866 = -t129 - t1016;
t625 = t636 * rSges(7,2);
t1113 = -t776 * t633 + t625;
t883 = qJD(6) * t636;
t541 = t606 * t883;
t891 = qJD(3) * t633;
t683 = -t606 * t891 + t607 * t886;
t958 = t541 + t683 * qJ(6) + (-t606 * t886 - t607 * t891) * pkin(5) + (-rSges(7,1) * t606 + rSges(7,3) * t607) * t886 + t1113 * qJD(3);
t9 = qJD(6) * t189 + qJDD(6) * t422 + t310 * t920 + t502 * t856 + t507 * t958 + t588 * t866 + t640;
t1133 = t50 * t958 + (qJD(1) * t51 + t9) * t920;
t1129 = rSges(7,2) * t839 + t1173 * t189 + t1174 * t190 + t394;
t1015 = -rSges(7,2) * t840 + t1129;
t1115 = -t519 - t612;
t846 = pkin(3) * t1163 + pkin(8) * t839;
t682 = pkin(8) * t840 - t846;
t687 = -t631 * t891 - t613;
t868 = t635 * t1039;
t803 = t1163 * t602 + t631 * t840 + t637 * t868;
t130 = (-pkin(4) * t715 + t687) * t634 + t682 + t803;
t525 = -t634 * t880 + t610;
t901 = qJ(2) * t892 + t614;
t851 = qJD(1) * (-pkin(1) * t894 + t901) + qJDD(1) * t548 + t634 * t881;
t805 = qJDD(1) * t628 + t851;
t663 = -qJD(1) * t682 + qJDD(1) * t492 - t525 * t557 + t805;
t639 = t588 * t130 + t502 * t287 + (-pkin(7) * t638 - qJD(1) * t613) * t634 + t663 + (qJD(3) * t1115 - qJDD(2) + t878) * t637;
t852 = t373 + t920;
t853 = -t276 - t958;
t8 = qJD(6) * t187 - qJDD(6) * t424 + t1015 * t588 - t311 * t852 + t502 * t949 + t508 * t853 + t639;
t1132 = t8 - g(2);
t1026 = t633 * rSges(6,3);
t1044 = rSges(6,1) * t607;
t780 = -rSges(6,2) * t606 + t1044;
t389 = t636 * t780 + t1026;
t215 = t423 * rSges(6,1) - t422 * rSges(6,2) - rSges(6,3) * t968;
t948 = t215 + t287;
t68 = -t389 * t508 + t588 * t948 + t657;
t926 = t373 + t389;
t1131 = t68 * t926;
t721 = t433 * t636 + t437 * t633;
t1123 = t721 * t637;
t465 = t475 * pkin(4);
t942 = t1173 * t423 - t1174 * t422;
t1118 = t465 + t942;
t466 = t477 * pkin(4);
t940 = -t1173 * t425 + t1174 * t424;
t1117 = t466 + t940;
t1114 = t1173 * t977;
t549 = -rSges(3,2) * t637 + t634 * rSges(3,3);
t547 = rSges(4,1) * t636 - rSges(4,2) * t633;
t491 = t547 * t637;
t786 = rSges(4,1) * t633 + rSges(4,2) * t636;
t302 = -qJD(3) * t491 + (t634 * t786 + t626) * qJD(1);
t513 = t786 * qJD(3);
t453 = -t634 * rSges(4,3) + t637 * t786;
t796 = t453 + t815;
t87 = -t513 * t890 + t524 * t547 + (-t302 - t464) * qJD(1) + t796 * qJDD(1) + t717;
t847 = rSges(4,1) * t1163 + rSges(4,2) * t840;
t303 = (-rSges(4,2) * t891 - rSges(4,3) * qJD(1)) * t634 + t847;
t875 = t638 * t1050;
t88 = -t875 + qJD(1) * t303 + qJDD(1) * t451 - t525 * t547 + (qJD(3) * t513 - qJDD(2)) * t637 + t805;
t1104 = t87 * t634 - t88 * t637;
t1095 = t1228 * t588 + (Icges(5,5) * t475 - Icges(5,6) * t476 - t1189 * t423 - t1190 * t422) * t508 + (Icges(5,5) * t477 + Icges(5,6) * t478 + t1189 * t425 + t1190 * t424) * t507;
t1089 = t1097 * t502 + t1181 * t588;
t243 = t433 * t633 - t437 * t636;
t296 = qJD(1) * t432 - qJD(3) * t485;
t299 = -qJD(3) * t488 + (t634 * t752 + t1002) * qJD(1);
t744 = Icges(4,5) * t633 + Icges(4,6) * t636;
t429 = -Icges(4,3) * t634 + t637 * t744;
t896 = qJD(1) * t429;
t1085 = qJD(3) * t243 + t296 * t636 + t299 * t633 + t896;
t510 = t748 * qJD(3);
t511 = t752 * qJD(3);
t534 = Icges(4,5) * t636 - Icges(4,6) * t633;
t718 = t536 * t633 - t538 * t636;
t1084 = qJD(1) * t534 + qJD(3) * t718 + t510 * t636 + t511 * t633;
t297 = qJD(1) * t433 + t536 * t890;
t300 = qJD(1) * t437 + qJD(3) * t487;
t722 = t432 * t633 - t436 * t636;
t428 = Icges(4,3) * t637 + t634 * t744;
t897 = qJD(1) * t428;
t1083 = qJD(3) * t722 - t297 * t636 - t300 * t633 + t897;
t1055 = pkin(4) * t632;
t1075 = -pkin(1) - pkin(7);
t810 = -t1055 + t1075;
t1082 = -t625 * t637 + t634 * t810;
t650 = t507 * (Icges(5,2) * t478 - t285 + t462) + t508 * (-Icges(5,2) * t476 + t284 + t461) + t588 * (t435 + t483);
t1079 = t507 * (-Icges(5,1) * t477 + t283 - t463) + t508 * (-Icges(5,1) * t475 + t1009 + t281) + t588 * (t431 - t486);
t1078 = m(6) / 0.2e1;
t1077 = m(7) / 0.2e1;
t1076 = -m(6) - m(7);
t1074 = t310 / 0.2e1;
t1073 = t311 / 0.2e1;
t1072 = t502 / 0.2e1;
t1071 = -t507 / 0.2e1;
t1070 = t507 / 0.2e1;
t1069 = -t508 / 0.2e1;
t1068 = t508 / 0.2e1;
t1067 = t524 / 0.2e1;
t1066 = t525 / 0.2e1;
t1065 = -t588 / 0.2e1;
t1064 = t588 / 0.2e1;
t1062 = t634 / 0.2e1;
t1061 = -t637 / 0.2e1;
t1058 = rSges(3,2) - pkin(1);
t1057 = -rSges(5,3) - pkin(8);
t1053 = g(1) * t634;
t1052 = g(2) * t637;
t1051 = g(3) * t636;
t1049 = t8 * t637;
t1048 = t9 * t634;
t1041 = rSges(3,3) * t637;
t779 = -rSges(6,1) * t606 - rSges(6,2) * t607;
t448 = t779 * t636;
t623 = t636 * rSges(6,3);
t227 = qJD(4) * t448 + (-t633 * t780 + t623) * qJD(3);
t945 = t217 - t288;
t782 = rSges(6,1) * t188 - rSges(6,2) * t187;
t115 = -rSges(6,3) * t1162 + t782;
t964 = -t115 - t129;
t11 = t227 * t507 + t310 * t389 + t502 * t945 + t588 * t964 + t640;
t1038 = t11 * t634;
t859 = t190 * rSges(6,1) - t189 * rSges(6,2) + rSges(6,3) * t839;
t117 = -rSges(6,3) * t840 + t859;
t938 = -t276 - t227;
t12 = t588 * t117 + t502 * t215 - t311 * t926 + t508 * t938 + t639;
t1037 = t12 * t637;
t1036 = t21 * t508;
t1035 = t22 * t507;
t1034 = t23 * t508;
t1033 = t24 * t507;
t1032 = t29 * t508;
t1031 = t30 * t507;
t624 = t636 * rSges(5,3);
t1025 = t68 * t389;
t1024 = t79 * t311;
t1023 = t80 * t310;
t1022 = t81 * t311;
t1021 = t82 * t310;
t1018 = t90 * t311;
t1017 = t91 * t310;
t496 = t547 * t890;
t183 = qJD(1) * t796 + t496 + t614;
t986 = t183 * t637;
t981 = t534 * t637;
t980 = t606 * t633;
t975 = t632 * t636;
t481 = t634 * t534;
t963 = -t117 - t130;
t944 = t633 * t287 + t373 * t968;
t460 = t637 * t494;
t943 = t637 * t288 - t460;
t269 = -t422 * rSges(6,1) - t423 * rSges(6,2);
t941 = -t269 - t465;
t273 = t424 * rSges(6,1) + t425 * rSges(6,2);
t939 = -t273 - t466;
t912 = t476 * rSges(5,1) + t475 * rSges(5,2);
t289 = -rSges(5,3) * t968 + t912;
t933 = -t289 - t492;
t863 = t607 * t968;
t864 = t606 * t968;
t931 = rSges(7,2) * t973 + t1173 * t864 + t1174 * t863;
t930 = -t388 * t637 - t775 * t966;
t512 = t634 * t557;
t928 = t634 * t373 + t512;
t876 = pkin(4) * t975;
t927 = t588 * t465 + t508 * t876;
t921 = -t775 * t633 + t1113;
t913 = t1174 * t979 - t1114;
t909 = rSges(6,1) * t863 + rSges(6,3) * t973;
t545 = rSges(3,2) * t634 + t1041;
t905 = -t544 + t545;
t459 = t548 + t549;
t904 = rSges(6,2) * t980 + t623;
t903 = t633 * t1042 + t624;
t900 = rSges(3,2) * t894 + rSges(3,3) * t892;
t899 = -qJD(1) * t544 + t614;
t898 = t615 - t580;
t895 = qJD(1) * t744;
t719 = t536 * t636 + t538 * t633;
t252 = t637 * t719 - t481;
t882 = t252 * qJD(1);
t877 = -rSges(4,3) + t1075;
t874 = t636 * t1045;
t873 = rSges(5,2) * t975;
t872 = rSges(6,2) * t979;
t865 = -t130 - t1015;
t857 = -t492 - t948;
t855 = t1162 * t287 + t288 * t839;
t854 = t250 * rSges(5,1) + t249 * rSges(5,2) + rSges(5,3) * t839;
t149 = t637 * t428 + t432 * t968 + t436 * t973;
t150 = -t637 * t429 - t433 * t968 - t437 * t973;
t850 = t504 + t898;
t848 = t617 + t908;
t843 = t636 * t1057;
t837 = t682 * t890;
t830 = qJD(6) * t980;
t826 = -pkin(3) - t1045;
t823 = -t893 / 0.2e1;
t821 = -t890 / 0.2e1;
t820 = t890 / 0.2e1;
t818 = -t888 / 0.2e1;
t817 = t888 / 0.2e1;
t812 = -t602 - t1044;
t808 = t633 * t130 + t276 * t968 + t287 * t889 + t636 * t352;
t807 = -t492 - t858;
t804 = qJD(1) * t494 + t498 + t899;
t793 = t1075 * t634 + t617;
t792 = -t1052 + t1053;
t791 = t323 * t888 - t524 * t492 - t525 * t494;
t790 = -t624 + t1056;
t789 = -t492 * t890 - t494 * t888;
t473 = qJD(1) * t493;
t787 = -t556 * t888 + t473;
t550 = rSges(2,1) * t637 - rSges(2,2) * t634;
t546 = rSges(2,1) * t634 + rSges(2,2) * t637;
t785 = rSges(5,1) * t248 + rSges(5,2) * t247;
t67 = t389 * t507 + t588 * t945 + t648;
t774 = t11 * t389 + t67 * t227;
t723 = t432 * t636 + t436 * t633;
t668 = qJD(1) * t723 + qJD(3) * t481 + t896;
t669 = -qJD(1) * t721 - qJD(3) * t981 + t897;
t773 = (t1083 * t637 + t668 * t634) * t637 + (-t1085 * t637 + t669 * t634) * t634;
t772 = (-t1083 * t634 + t668 * t637) * t637 + t634 * (t1085 * t634 + t669 * t637);
t753 = -t588 * t288 + t804;
t740 = t149 * t637 + t150 * t634;
t390 = t634 * t428;
t152 = -t637 * t723 + t390;
t153 = -t429 * t634 + t1123;
t739 = t152 * t637 + t153 * t634;
t184 = qJD(1) * t1157 - t547 * t888 - t615;
t738 = t183 * t634 - t184 * t637;
t728 = t289 * t637 - t290 * t634;
t727 = t302 * t637 - t303 * t634;
t720 = -t451 * t634 - t453 * t637;
t711 = t848 + t990;
t709 = t810 * t1053;
t708 = t803 + t901;
t707 = t1112 + t849;
t700 = -t602 * t633 - qJ(2) - t976;
t490 = (-rSges(5,1) * t632 - rSges(5,2) * t635) * t636;
t689 = t68 * t215 + t67 * t945;
t686 = t508 * t288 + t612 + t789;
t58 = -t217 * t508 - t507 * t948 + t686;
t688 = t58 * t215 - t67 * t926;
t362 = rSges(5,3) * t973 + (-t873 + t874) * t634;
t667 = t719 * qJD(1) - t744 * qJD(3);
t664 = qJD(3) * t613 + qJDD(5) * t633 + t508 * t129 + t311 * t288 + t791;
t661 = t50 * t856 + t51 * t949;
t46 = -t507 * t858 - t508 * t947 + t541 + t686;
t660 = t46 * t949 - t50 * t852;
t656 = -t633 * t869 + t687;
t112 = t1193 + t666;
t113 = t289 * t588 - t452 * t508 + t665;
t89 = -t289 * t507 - t290 * t508 + t789;
t646 = t89 * t728 + (-t112 * t637 - t113 * t634) * t452;
t645 = t688 * t637 + (-t217 * t58 - t1131) * t634;
t644 = t660 * t637 + (-t46 * t947 - t51 * t852) * t634;
t540 = t637 * t873;
t523 = t637 * t872;
t489 = t547 * t634;
t450 = -t1045 * t633 + t903;
t387 = -t1044 * t633 + t904;
t363 = t540 + (-t874 - t1040) * t637;
t361 = t435 * t637;
t360 = t435 * t634;
t359 = t431 * t637;
t358 = t431 * t634;
t354 = qJD(1) * t459 - t615;
t353 = qJD(1) * t905 + t614;
t350 = -t508 * t633 - t588 * t968;
t349 = t507 * t633 - t588 * t966;
t348 = t373 * t966;
t346 = t508 * t466;
t345 = t523 + (-rSges(6,1) * t977 - t1026) * t637;
t343 = -rSges(6,2) * t864 + t909;
t322 = rSges(5,1) * t477 + rSges(5,2) * t478;
t321 = rSges(5,1) * t475 - rSges(5,2) * t476;
t305 = (t507 * t634 + t508 * t637) * t636;
t301 = qJD(4) * t490 + (-t633 * t783 + t624) * qJD(3);
t251 = t634 * t719 + t981;
t239 = t251 * qJD(1);
t236 = t720 * qJD(3);
t235 = t276 * t966;
t167 = qJD(1) * t900 + qJDD(1) * t549 - qJDD(2) * t637 + t851;
t166 = t905 * qJDD(1) + (-qJD(1) * t549 - t464) * qJD(1) + t902;
t138 = -rSges(5,3) * t840 + t854;
t137 = -rSges(5,3) * t1162 + t785;
t99 = -t1084 * t634 + t667 * t637;
t98 = t1084 * t637 + t667 * t634;
t97 = qJD(3) * t721 - t296 * t633 + t299 * t636;
t96 = -qJD(3) * t723 - t297 * t633 + t300 * t636;
t70 = qJD(3) * t739 - t882;
t69 = qJD(3) * t740 + t239;
t41 = -t875 + t138 * t588 + t289 * t502 - t301 * t508 - t311 * t452 + (-qJD(3) * t519 - qJDD(2)) * t637 + t663;
t40 = qJD(1) * t932 - t137 * t588 + t290 * t502 + t301 * t507 + t310 * t452 + t647;
t35 = t137 * t508 - t138 * t507 - t289 * t310 - t290 * t311 + t791 + t837;
t10 = t115 * t508 - t217 * t311 - t310 * t948 + t507 * t963 + t664 + t837;
t7 = (qJD(4) * qJD(6) * t607 + qJDD(6) * t606) * t636 + t1016 * t508 - t947 * t311 + (t634 * t682 - t830) * qJD(3) + t865 * t507 - t858 * t310 + t664;
t1 = [(t98 + t69 + t97) * t820 + (t239 + ((-t152 + t390 + t150) * t634 + (t153 - t1123 + (t429 - t723) * t634 + t149) * t637) * qJD(3)) * t821 + t1098 * t1074 + t1099 * t1073 + t1100 * t1068 + (-t718 + m(2) * (t546 ^ 2 + t550 ^ 2) + Icges(2,3) + Icges(3,1)) * qJDD(1) + t1033 / 0.2e1 + t1034 / 0.2e1 + (t1101 + t1103) * t1070 + t1035 / 0.2e1 + t1036 / 0.2e1 + ((t1197 + t1233) * t508 + t1178) * t1071 + ((-t1141 + t1233) * t507 + t1102 + t1234) * t1069 + t1031 / 0.2e1 + t1032 / 0.2e1 + (-qJD(3) * t719 + t510 * t633 - t511 * t636) * qJD(1) + t1089 + (t251 - t722) * t1066 + (-g(1) * (t711 + t1120) - t709 + (t1082 + t711 - t1109) * t9 + (-t852 * t507 - t947 * t588 + t656 * t634 - t394 + t708 - t753 + t831 + (t1050 + t1082) * qJD(1) + t1129) * t51 + (-t308 - t898 + t850 + ((-t869 + (rSges(7,2) - t631) * qJD(3)) * t633 + t810 * qJD(1)) * t637 + (-t868 + (t700 + t625) * qJD(1)) * t634 + t1086 + t1130 + t1151) * t50 + t1132 * (t707 + t949)) * m(7) + (t96 + t99) * t817 + (t353 * t615 + t354 * (t900 + t901) + (t353 * t1058 * t637 + (t353 * (-rSges(3,3) - qJ(2)) - t354 * pkin(1)) * t634) * qJD(1) - (qJD(1) * t545 - t353 + t899) * t354 + (t167 - g(2)) * t459 + (t166 - g(1)) * (t1058 * t634 + t1041 + t617)) * m(3) + (t112 * (rSges(5,3) * t838 + t615 - t785 + t845) + t113 * (t846 + t854 + t901) + ((t1075 * t112 + t113 * t843) * t637 + (t112 * (-qJ(2) - t790) + t113 * t1075) * t634) * qJD(1) - (-pkin(7) * t894 - t112 + t1193 + t804) * t113 + (t41 - g(2)) * (t634 * t843 + t1112 + t590 + t912) + (t40 - g(1)) * (t637 * t790 + t1116 - t595 + t793)) * m(5) + (t882 + (t429 * t634 ^ 2 + (-t390 + t150 + (t429 + t723) * t637) * t637) * qJD(3) + t70) * t818 + ((-t1121 * t607 - t1122 * t606 + t729) * t636 - t1096 * t633 + t1195) * t508 * t1065 + t1024 / 0.2e1 + t1022 / 0.2e1 + t1023 / 0.2e1 + t1021 / 0.2e1 + t243 * t1067 + t1017 / 0.2e1 + t1018 / 0.2e1 - t524 * t252 / 0.2e1 + (t183 * (rSges(4,1) * t835 - rSges(4,2) * t838 + t615) + t184 * (-rSges(4,2) * t839 + t847 + t901) + (t877 * t986 + (t183 * (-qJ(2) - t786) + t184 * t877) * t634) * qJD(1) - (-t183 + t496 + (t453 - t1050) * qJD(1) + t899) * t184 + (t88 - g(2)) * t1157 + (t87 - g(1)) * (t453 + t793)) * m(4) + (-(t588 * t217 - t67 + (-pkin(7) * qJD(1) - t613) * t634 + t753) * t68 - t1131 * t507 + t11 * (-t781 + t848) + t67 * (-t782 + t850) + t68 * (t708 + t859) + (-t11 * t623 + t67 * (-t869 + (rSges(6,3) - t631) * qJD(3)) * t633 + (-t623 * t68 + t67 * t810) * qJD(1)) * t637 + (t11 * t810 - t67 * t868 + t68 * t656 + (t67 * (t700 + t623) + t68 * t810) * qJD(1)) * t634 - g(1) * (t217 + t848) - t709 + (t12 - g(2)) * (t215 + t707)) * m(6) - m(2) * (-g(1) * t546 + g(2) * t550); (-m(3) - m(5) + t1076) * t792 + 0.2e1 * (t1048 / 0.2e1 - t1049 / 0.2e1) * m(7) + 0.2e1 * (t1038 / 0.2e1 - t1037 / 0.2e1) * m(6) + 0.2e1 * (t1061 * t41 + t1062 * t40) * m(5) + 0.2e1 * (t1061 * t167 + t1062 * t166) * m(3) + (-t792 + t1104) * m(4); (-(t633 * t644 + t636 * t661) * qJD(4) - g(1) * (t528 + t931) - g(2) * t561 - g(3) * (t625 - t976) - (-t631 * t1053 + g(3) * (-t602 - t775 - t776)) * t633 - (-t1027 + (-t1173 * t606 - t1174 * t607 - t602) * t636) * t1052 + t9 * t928 + t7 * t943 + (t7 * t807 + t1133) * t634 + (t8 * (-t557 - t852) - t7 * t947) * t637 + (-t473 - t931 * t588 - (-t372 - t921) * t508 + (t541 - t1115 - t519 + t853) * t637 + t1155) * t51 + (-t541 * t634 - (-t320 - t930) * t588 - t921 * t507 + t920 * t892 + t1166) * t50 + (t830 - t930 * t508 - (-t319 - t931) * t507 + (t682 + t865 + (t494 + t856) * qJD(1)) * t634 + (qJD(1) * t807 + t1016) * t637 + t1161) * t46) * m(7) + qJD(1) * (t634 * t97 + t637 * t96 + (t243 * t637 + t634 * t722) * qJD(1)) / 0.2e1 + qJDD(1) * (t243 * t634 - t637 * t722) / 0.2e1 + ((qJD(3) * t727 - t451 * t524 - t453 * t525) * t720 + t236 * ((-t451 * t637 + t453 * t634) * qJD(1) + t727) - t738 * t513 + ((t184 * t634 + t986) * qJD(1) + t1104) * t547 - (t183 * t491 + t184 * t489) * qJD(1) - (t236 * (-t489 * t634 - t491 * t637) - t738 * t786) * qJD(3) - g(1) * t489 + g(2) * t491 + g(3) * t786) * m(4) - qJD(1) * ((-t633 * t907 - t636 * t906) * qJD(1) + (t633 * t670 + t636 * t671) * qJD(3)) / 0.2e1 - t1142 * t886 / 0.2e1 + (t70 + t1102) * t892 / 0.2e1 + t1102 * t833 / 0.2e1 - (t69 + t1103) * t894 / 0.2e1 - t1103 * t834 / 0.2e1 + (-qJD(1) * t1094 + t1143 * t634 + t1144 * t637) * t1064 + ((-t152 * t634 + t153 * t637) * qJD(1) + t773) * t820 + t740 * t1066 + t739 * t1067 + (((t1167 * t607 + t1168 * t606 - t430 * t632 + t434 * t635 + t1230) * t588 + (t1170 * t607 + t1172 * t606 - t358 * t632 + t360 * t635 + t1212) * t508 + (-t1169 * t607 - t1171 * t606 + t359 * t632 - t361 * t635 + t1096) * t507 + t1097 * qJD(4)) * t636 + (qJD(4) * t1094 + t1205) * t633) * t1065 + (-qJD(1) * t1092 + t1145 * t634 + t1146 * t637) * t1068 + (-qJD(1) * t1093 + t1147 * t634 + t1148 * t637) * t1070 + (t11 * t928 + t10 * t943 + (t1025 * qJD(1) + t10 * t857 + t774) * t634 + (t12 * (-t557 - t926) - t10 * t217) * t637 - g(1) * (t528 + (-t631 * t633 - t872) * t634 + t909) - g(2) * (t523 + t561) - g(3) * (t633 * t812 + t904 - t976) - (t636 * t812 - t1026) * t1052 - (t633 * t645 + t636 * t689) * qJD(4) + (-t343 * t588 + t633 * t884 - t787 - (-t372 - t387) * t508 + (-t519 + t938) * t637 + t1155) * t68 + (-t387 * t507 - (-t320 - t345) * t588 + t389 * t892 + t1166) * t67 + (-t345 * t508 - (-t319 - t343) * t507 + (t682 + t963 + (t494 + t945) * qJD(1)) * t634 + (qJD(1) * t857 + t115) * t637 + t1161) * t58) * m(6) + (t40 * t512 + t112 * t910 + t113 * t500 - t35 * t460 + (t112 * t301 + t113 * t1159 + t35 * t933 + t40 * t452) * t634 + (t41 * (-t452 - t557) + t113 * (-t301 - t519) - t35 * t290 + t112 * t1159) * t637 - g(1) * (t362 + t493) - g(2) * t540 - g(3) * (t633 * t826 + t627 + t903) - (t1057 * t633 + t636 * t826) * t1052 - t112 * (-t363 * t588 + t450 * t507 + t811) - t113 * (t362 * t588 - t450 * t508 + t787) - ((t112 * t290 + t113 * t289) * t636 + t646 * t633) * qJD(4) + ((-t138 + t682 + (t290 + t494) * qJD(1)) * t634 + (qJD(1) * t933 + t137) * t637 + t362 * t507 - t363 * t508 + t1188) * t89) * m(5) + ((t481 * t888 - t895) * t637 + (-t1152 + (-t637 * t981 - t1150) * qJD(3)) * t634) * t818 + ((-t890 * t981 - t895) * t634 + (t1152 + (t634 * t481 + t1150) * qJD(3)) * t637) * t821 + (qJD(1) * t99 + qJD(3) * t772 + qJDD(1) * t251 + t149 * t525 + t150 * t524 + t1153) * t637 / 0.2e1 + (qJD(1) * t98 + qJD(3) * t773 - qJDD(1) * t252 + t152 * t525 + t153 * t524 + t1154) * t1062 + ((-t149 * t634 + t150 * t637) * qJD(1) + t772) * t817 + ((t1098 * t636 + t1139 * t973) * qJD(4) + ((-qJD(4) * t1197 + t1137) * t633 + t1177) * t637 + (-t1167 * t425 - t1168 * t424 + t430 * t477 - t434 * t478) * t588 + (-t1170 * t425 - t1172 * t424 + t358 * t477 - t360 * t478) * t508 + (t1169 * t425 + t1171 * t424 - t359 * t477 + t361 * t478) * t507) * t1071 + ((t1099 * t636 - t1198 * t972) * qJD(4) + ((qJD(4) * t1141 - t1137) * t633 - t1177) * t634 + (t1167 * t423 + t1168 * t422 + t430 * t475 + t434 * t476) * t588 + (t1170 * t423 + t1172 * t422 + t358 * t475 + t360 * t476) * t508 + (-t1169 * t423 - t1171 * t422 - t359 * t475 - t361 * t476) * t507) * t1069 + t1134 * t1073 + t1135 * t1074 + t1136 * t1072; -t1153 * t968 / 0.2e1 + t1154 * t966 / 0.2e1 + t1142 * t889 / 0.2e1 + (-t1093 * t636 + t1098 * t633) * t1074 + (-t1092 * t636 + t1099 * t633) * t1073 + (-t1094 * t636 + t1097 * t633) * t1072 + (t1079 * t478 + t1095 * t966 + t1191 * t424 + t1192 * t425 + t650 * t477) * t1071 + ((-qJD(1) * t1135 + t1098 * qJD(3) + t1147 * t637 - t1148 * t634) * t636 + (qJD(3) * t1093 + t1101) * t633) * t1070 + (-t1079 * t476 - t1095 * t968 - t1191 * t422 - t1192 * t423 + t475 * t650) * t1069 + ((-qJD(1) * t1134 + t1099 * qJD(3) + t1145 * t637 - t1146 * t634) * t636 + (qJD(3) * t1092 + t1100) * t633) * t1068 + ((-t1079 * t635 - t1191 * t606 - t1192 * t607 - t632 * t650) * t636 + t1095 * t633) * t1065 + ((-qJD(1) * t1136 + t1097 * qJD(3) + t1143 * t637 - t1144 * t634) * t636 + (qJD(3) * t1094 + t1181) * t633) * t1064 + (t1023 + t1024 + t1035 + t1036 + t1021 + t1022 + t1033 + t1034 + t1017 + t1018 + t1031 + t1032 + t1089) * t633 / 0.2e1 + (-g(1) * t1118 - g(2) * t1117 - g(3) * t1114 - (-t1174 * t606 - t1055) * t1051 + t9 * t348 + t8 * t944 + (qJD(3) * t644 + t8 * t949 + t856 * t9) * t633 + (t661 * qJD(3) + (-t7 * t858 + t1133) * t637 + (qJD(1) * t660 + t7 * t856 + t8 * t920) * t634) * t636 + (qJD(6) * t425 + t1015 * t633 - t508 * t913 - t588 * t942 + t958 * t968 + t808 - t927) * t51 + (-t607 * t883 - t346 - t940 * t508 + t1118 * t507 + t855 + ((qJD(1) * t856 + t865) * t637 + t866 * t634) * t636) * t46 + (-qJD(6) * t423 + t1117 * t588 - (-t876 - t913) * t507 + t235 + t866 * t633) * t50) * m(7) + (t11 * t348 + t67 * t235 + t12 * t944 + t68 * t808 + t58 * t855 + (qJD(3) * t645 + t11 * t945 + t68 * t117 + t12 * t215 + t67 * t964) * t633 + (t689 * qJD(3) + (-t10 * t948 + t58 * t963 + (t58 * t945 + t1025) * qJD(1) + t774) * t637 + (qJD(1) * t688 + t10 * t945 + t12 * t389 + t68 * t227 + t58 * t964) * t634) * t636 - t68 * (-t448 * t508 + t927) - t58 * (t273 * t508 + t346) - (t68 * t269 + t67 * t939) * t588 - (t67 * (t448 - t876) + t58 * t941) * t507 + g(1) * t941 + g(2) * t939 - (t779 - t1055) * t1051) * m(6) + (-t112 * (-t322 * t588 + t490 * t507) - t113 * (t321 * t588 - t490 * t508) - t89 * (-t321 * t507 + t322 * t508) + (qJD(3) * t646 - t112 * t137 + t113 * t138 + t41 * t289 + t290 * t40) * t633 + (t112 * (qJD(3) * t290 + t301 * t637) + t113 * (qJD(3) * t289 + t301 * t634) - t35 * t728 + t89 * (-t137 * t634 - t138 * t637 + t289 * t894 + t290 * t892) + (t40 * t637 + t41 * t634 + (-t112 * t634 + t113 * t637) * qJD(1)) * t452) * t636 - g(1) * t321 - g(2) * t322 - g(3) * t490) * m(5) + (t633 * t820 + t637 * t823) * t1103 + (t633 * t818 + t634 * t823) * t1102; t1076 * (g(3) * t633 - t636 * t792) - m(6) * (t305 * t58 + t349 * t67 + t350 * t68) - m(7) * (t305 * t46 + t349 * t50 + t350 * t51) + 0.2e1 * ((t67 * t890 - t68 * t888 + t10) * t1078 + (t50 * t890 - t51 * t888 + t7) * t1077) * t633 + 0.2e1 * ((qJD(3) * t58 - t67 * t892 - t68 * t894 + t1037 - t1038) * t1078 + (qJD(3) * t46 - t50 * t892 - t51 * t894 - t1048 + t1049) * t1077) * t636; (t187 * t51 + t189 * t50 + t683 * t46 + (-t50 * t507 + t508 * t51 - g(3) + t7) * t979 + (t46 * t508 - t50 * t588 - t1132) * t424 + (t46 * t507 - t51 * t588 - g(1) + t9) * t422) * m(7);];
tau  = t1;
