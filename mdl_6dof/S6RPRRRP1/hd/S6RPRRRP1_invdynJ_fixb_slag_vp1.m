% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP1
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:32
% EndTime: 2019-03-09 05:57:05
% DurationCPUTime: 83.61s
% Computational Cost: add. (62045->1338), mult. (56053->1687), div. (0->0), fcn. (52471->10), ass. (0->666)
t1183 = Icges(6,1) + Icges(7,1);
t1208 = Icges(7,4) + Icges(6,5);
t1207 = Icges(6,6) - Icges(7,6);
t638 = sin(qJ(5));
t1235 = (Icges(6,4) - Icges(7,5)) * t638;
t1182 = Icges(6,2) + Icges(7,3);
t1231 = Icges(7,2) + Icges(6,3);
t641 = cos(qJ(5));
t1234 = -t1207 * t638 + t1208 * t641;
t1233 = t1183 * t641 - t1235;
t637 = qJ(3) + qJ(4);
t630 = sin(t637);
t631 = cos(t637);
t1034 = Icges(6,4) * t641;
t787 = -Icges(6,2) * t638 + t1034;
t1232 = -t1207 * t631 + t630 * t787;
t636 = qJ(1) + pkin(10);
t628 = sin(t636);
t629 = cos(t636);
t971 = t631 * t641;
t470 = t628 * t971 - t629 * t638;
t434 = Icges(7,5) * t470;
t972 = t631 * t638;
t978 = t629 * t641;
t469 = t628 * t972 + t978;
t986 = t628 * t630;
t238 = -Icges(7,6) * t986 - Icges(7,3) * t469 - t434;
t437 = Icges(6,4) * t470;
t246 = -Icges(6,2) * t469 + Icges(6,6) * t986 + t437;
t1219 = t238 + t246;
t472 = t628 * t638 + t629 * t971;
t435 = Icges(7,5) * t472;
t983 = t628 * t641;
t471 = t629 * t972 - t983;
t981 = t629 * t630;
t239 = Icges(7,6) * t981 + Icges(7,3) * t471 + t435;
t1036 = Icges(6,4) * t472;
t248 = -Icges(6,2) * t471 + Icges(6,6) * t981 + t1036;
t1218 = t239 - t248;
t240 = Icges(6,5) * t470 - Icges(6,6) * t469 + Icges(6,3) * t986;
t243 = Icges(7,4) * t470 + Icges(7,2) * t986 + Icges(7,6) * t469;
t1217 = t240 + t243;
t242 = Icges(6,5) * t472 - Icges(6,6) * t471 + Icges(6,3) * t981;
t245 = Icges(7,4) * t472 + Icges(7,2) * t981 + Icges(7,6) * t471;
t1192 = t242 + t245;
t433 = Icges(7,5) * t469;
t249 = Icges(7,1) * t470 + Icges(7,4) * t986 + t433;
t436 = Icges(6,4) * t469;
t253 = -Icges(6,1) * t470 - Icges(6,5) * t986 + t436;
t1216 = t249 - t253;
t1029 = Icges(7,5) * t471;
t251 = Icges(7,1) * t472 + Icges(7,4) * t981 + t1029;
t438 = Icges(6,4) * t471;
t254 = Icges(6,1) * t472 + Icges(6,5) * t981 - t438;
t1215 = t251 + t254;
t974 = t630 * t641;
t583 = Icges(7,5) * t974;
t975 = t630 * t638;
t1206 = Icges(7,3) * t975 - t1232 + t583;
t1190 = -t1231 * t631 + t1234 * t630;
t1188 = -t1208 * t631 + t1233 * t630;
t1027 = Icges(7,5) * t641;
t781 = Icges(7,3) * t638 + t1027;
t1230 = (t781 - t787) * t631;
t1229 = t1234 * t631;
t1228 = t1233 * t631;
t1227 = t1182 * t641 + t1235;
t1226 = -t1207 * t641 - t1208 * t638;
t640 = sin(qJ(1));
t1057 = t640 * pkin(1);
t644 = -pkin(8) - pkin(7);
t598 = t629 * t644;
t642 = cos(qJ(3));
t633 = t642 * pkin(3);
t624 = t633 + pkin(2);
t930 = -t628 * t624 - t598;
t813 = t930 - t1057;
t1172 = rSges(7,1) + pkin(5);
t1162 = rSges(7,3) + qJ(6);
t1197 = t1216 * t472 + t1217 * t981 - t1219 * t471;
t1196 = t1192 * t981 + t1215 * t472 + t1218 * t471;
t908 = qJD(5) * t631;
t864 = t629 * t908;
t635 = qJD(3) + qJD(4);
t970 = t635 * t638;
t885 = t630 * t970;
t906 = qJD(5) * t638;
t217 = qJD(1) * t469 - t628 * t906 + t629 * t885 - t641 * t864;
t588 = qJD(1) - t908;
t969 = t635 * t641;
t684 = t588 * t638 - t630 * t969;
t914 = qJD(1) * t631;
t825 = -qJD(5) + t914;
t218 = t629 * t684 - t825 * t983;
t915 = qJD(1) * t630;
t871 = t628 * t915;
t973 = t631 * t635;
t886 = t629 * t973;
t700 = -t871 + t886;
t130 = Icges(7,5) * t218 + Icges(7,6) * t700 - Icges(7,3) * t217;
t136 = Icges(6,4) * t218 + Icges(6,2) * t217 + Icges(6,6) * t700;
t1225 = t130 - t136;
t905 = qJD(5) * t641;
t916 = qJD(1) * t629;
t917 = qJD(1) * t628;
t219 = -t628 * t885 - t629 * t906 - t641 * t917 + (t628 * t905 + t638 * t916) * t631;
t220 = t628 * t684 + t825 * t978;
t888 = t628 * t973;
t702 = t629 * t915 + t888;
t131 = Icges(7,5) * t220 + Icges(7,6) * t702 + Icges(7,3) * t219;
t137 = Icges(6,4) * t220 - Icges(6,2) * t219 + Icges(6,6) * t702;
t1224 = t131 - t137;
t132 = Icges(6,5) * t218 + Icges(6,6) * t217 + Icges(6,3) * t700;
t134 = Icges(7,4) * t218 + Icges(7,2) * t700 - Icges(7,6) * t217;
t1223 = t132 + t134;
t133 = Icges(6,5) * t220 - Icges(6,6) * t219 + Icges(6,3) * t702;
t135 = Icges(7,4) * t220 + Icges(7,2) * t702 + Icges(7,6) * t219;
t1222 = t133 + t135;
t138 = Icges(7,1) * t218 + Icges(7,4) * t700 - Icges(7,5) * t217;
t140 = Icges(6,1) * t218 + Icges(6,4) * t217 + Icges(6,5) * t700;
t1221 = t138 + t140;
t139 = Icges(7,1) * t220 + Icges(7,4) * t702 + Icges(7,5) * t219;
t141 = Icges(6,1) * t220 - Icges(6,4) * t219 + Icges(6,5) * t702;
t1220 = t139 + t141;
t1120 = t1188 * t472 + t1190 * t981 + t1206 * t471;
t1214 = t1230 * t635 + (qJD(5) * t1227 - t1207 * t635) * t630;
t1213 = t1229 * t635 + (qJD(5) * t1226 + t1231 * t635) * t630;
t791 = -Icges(6,1) * t638 - t1034;
t1212 = t1228 * t635 + (t1208 * t635 + (-Icges(7,1) * t638 + t1027 + t791) * qJD(5)) * t630;
t985 = t628 * t631;
t896 = rSges(5,1) * t985;
t1211 = -t896 + t813;
t1210 = t1188 * t641 + t1206 * t638;
t1175 = t246 * t638 + t253 * t641;
t1176 = t238 * t638 - t249 * t641;
t1209 = t1175 + t1176;
t643 = cos(qJ(1));
t634 = t643 * pkin(1);
t1199 = t1216 * t470 + t1217 * t986 - t1219 * t469;
t1198 = t1192 * t986 + t1215 * t470 + t1218 * t469;
t1121 = t1188 * t470 + t1190 * t986 + t1206 * t469;
t1063 = -rSges(7,2) - pkin(9);
t623 = t631 * pkin(4);
t738 = t1063 * t630 - t623;
t516 = rSges(3,1) * t628 + rSges(3,2) * t629;
t501 = -t516 - t1057;
t803 = t470 * rSges(6,1) - t469 * rSges(6,2);
t259 = -rSges(6,3) * t986 - t803;
t510 = t629 * t635;
t909 = qJD(5) * t630;
t407 = -t628 * t909 + t510;
t1052 = rSges(6,1) * t641;
t802 = -rSges(6,2) * t638 + t1052;
t426 = -rSges(6,3) * t631 + t630 * t802;
t1205 = t259 * t588 - t407 * t426;
t432 = qJD(6) * t471;
t1191 = -t1162 * t469 - t1172 * t470;
t963 = -rSges(7,2) * t986 + t1191;
t1204 = t588 * t963 + t432;
t901 = qJD(1) * qJD(3);
t504 = qJDD(3) * t628 + t629 * t901;
t900 = qJD(1) * qJD(4);
t378 = qJDD(4) * t628 + t629 * t900 + t504;
t899 = qJDD(5) * t630;
t187 = qJD(5) * t700 + t629 * t899 + t378;
t603 = qJD(3) * t628;
t509 = qJD(4) * t628 + t603;
t406 = t629 * t909 + t509;
t907 = qJD(5) * t635;
t486 = -qJDD(5) * t631 + t630 * t907 + qJDD(1);
t514 = pkin(9) * t886;
t976 = t630 * t635;
t887 = t629 * t976;
t701 = -t628 * t914 - t887;
t278 = pkin(4) * t701 - pkin(9) * t871 + t514;
t980 = t629 * t631;
t556 = pkin(4) * t980;
t466 = pkin(9) * t981 + t556;
t1113 = t630 * pkin(9) + t623;
t478 = t1113 * t635;
t530 = pkin(4) * t630 - pkin(9) * t631;
t1054 = pkin(2) - t624;
t599 = pkin(7) * t916;
t639 = sin(qJ(3));
t911 = qJD(3) * t639;
t866 = t629 * t911;
t826 = pkin(3) * t866;
t301 = -t826 - t599 + (t1054 * t628 - t598) * qJD(1);
t615 = t628 * pkin(7);
t521 = t629 * pkin(2) + t615;
t829 = t629 * t624 - t628 * t644;
t375 = t829 - t521;
t646 = qJD(1) ^ 2;
t817 = qJDD(1) * t634 - t1057 * t646;
t724 = qJD(1) * (-pkin(2) * t917 + t599) + qJDD(1) * t521 + t817;
t968 = t642 * qJD(3) ^ 2;
t657 = qJD(1) * t301 + qJDD(1) * t375 + (-t504 * t639 - t628 * t968) * pkin(3) + t724;
t655 = qJD(1) * t278 + qJDD(1) * t466 - t378 * t530 - t509 * t478 + t657;
t798 = rSges(7,1) * t641 + rSges(7,3) * t638;
t1187 = (-pkin(5) * t641 - qJ(6) * t638 - t798) * t630;
t939 = -rSges(7,2) * t631 - t1187;
t904 = qJD(6) * t638;
t579 = t630 * t904;
t884 = t631 * t970;
t699 = t630 * t905 + t884;
t960 = t798 * t973 + (rSges(7,2) * t635 + (-rSges(7,1) * t638 + rSges(7,3) * t641) * qJD(5)) * t630 + t579 + t699 * qJ(6) + (-t630 * t906 + t631 * t969) * pkin(5);
t962 = rSges(7,2) * t981 + t1162 * t471 + t1172 * t472;
t1181 = rSges(7,2) * t886 - t1162 * t217 + t1172 * t218 + t432;
t966 = -rSges(7,2) * t871 + t1181;
t24 = qJD(6) * t219 + qJDD(6) * t469 - t187 * t939 - t406 * t960 + t486 * t962 + t588 * t966 + t655;
t1203 = t24 - g(2);
t898 = -qJDD(3) - qJDD(4);
t592 = t628 * t901;
t922 = t628 * t900 + t592;
t188 = (t631 * t907 + t899) * t628 + (qJD(1) * t909 + t898) * t629 + t922;
t889 = t628 * t976;
t513 = pkin(4) * t889;
t279 = pkin(9) * t702 + qJD(1) * t556 - t513;
t379 = t629 * t898 + t922;
t1061 = pkin(3) * t639;
t505 = -qJDD(3) * t629 + t592;
t897 = t646 * t634;
t688 = -pkin(3) * t629 * t968 + t505 * t1061 - t897;
t464 = t1113 * t628;
t616 = t629 * pkin(7);
t520 = pkin(2) * t628 - t616;
t374 = t520 + t930;
t847 = -t520 - t1057;
t816 = t374 + t847;
t745 = -t464 + t816;
t891 = pkin(3) * t911;
t558 = t628 * t891;
t926 = t644 * t917 + t558;
t302 = (-t1054 * t629 - t615) * qJD(1) - t926;
t503 = t521 * qJD(1);
t958 = -t302 - t503;
t654 = t379 * t530 + (-t279 + t958) * qJD(1) + t745 * qJDD(1) - t510 * t478 + t688;
t431 = qJD(6) * t469;
t1180 = -t1162 * t219 - t1172 * t220 - t431;
t965 = -rSges(7,2) * t702 + t1180;
t25 = -qJD(6) * t217 + qJDD(6) * t471 + t188 * t939 - t407 * t960 + t486 * t963 + t588 * t965 + t654;
t1202 = t25 - g(1);
t1171 = t1216 * t218 + t1217 * t700 + t1219 * t217 + t1220 * t472 + t1222 * t981 + t1224 * t471;
t1170 = t1192 * t700 + t1215 * t218 - t1218 * t217 + t1221 * t472 + t1223 * t981 + t1225 * t471;
t1169 = t1216 * t220 + t1217 * t702 - t1219 * t219 + t1220 * t470 + t1222 * t986 + t1224 * t469;
t1168 = t1192 * t702 + t1215 * t220 + t1218 * t219 + t1221 * t470 + t1223 * t986 + t1225 * t469;
t1164 = t1120 * t588 + t1196 * t406 - t1197 * t407;
t1201 = t1188 * t218 + t1190 * t700 - t1206 * t217 + t1212 * t472 + t1213 * t981 + t1214 * t471;
t1200 = t1188 * t220 + t1190 * t702 + t1206 * t219 + t1212 * t470 + t1213 * t986 + t1214 * t469;
t1158 = t243 * t631;
t115 = -t1176 * t630 - t1158;
t1161 = t240 * t631;
t117 = -t1175 * t630 - t1161;
t1195 = t115 + t117;
t761 = t239 * t638 + t251 * t641;
t116 = -t245 * t631 + t630 * t761;
t759 = -t248 * t638 + t254 * t641;
t118 = -t242 * t631 + t630 * t759;
t1194 = t116 + t118;
t1119 = -t1190 * t631 + t1210 * t630;
t1189 = -t630 * t781 + t1232;
t1186 = (t1231 * t630 - t1210 + t1229) * t588 + (t1190 * t628 - t1209) * t407 + (-t1190 * t629 - t759 - t761) * t406;
t1185 = -t406 * t939 + t962 * t588 + t431;
t935 = -t1162 * t974 + t1172 * t975;
t1184 = (t1210 * t635 - t1213) * t631 + (t1212 * t641 + t1214 * t638 + t1190 * t635 + (-t1188 * t638 + t1206 * t641) * qJD(5)) * t630;
t1179 = (-t1207 * t470 - t1208 * t469) * t407 + (t1207 * t472 + t1208 * t471) * t406 - t1226 * t588 * t630;
t1174 = t1120 * t486 + t1170 * t406 - t1171 * t407 + t1196 * t187 + t1197 * t188 + t1201 * t588;
t1173 = t1121 * t486 + t1168 * t406 - t1169 * t407 + t1198 * t187 + t1199 * t188 + t1200 * t588;
t37 = (-t1176 * t635 - t135) * t631 + (t131 * t638 + t139 * t641 + t243 * t635 + (-t238 * t641 - t249 * t638) * qJD(5)) * t630;
t39 = (-t1175 * t635 - t133) * t631 + (-t137 * t638 + t141 * t641 + t240 * t635 + (-t246 * t641 + t253 * t638) * qJD(5)) * t630;
t1167 = t37 + t39;
t38 = (t635 * t761 - t134) * t631 + (t130 * t638 + t138 * t641 + t245 * t635 + (t239 * t641 - t251 * t638) * qJD(5)) * t630;
t40 = (t635 * t759 - t132) * t631 + (-t136 * t638 + t140 * t641 + t242 * t635 + (-t248 * t641 - t254 * t638) * qJD(5)) * t630;
t1166 = t38 + t40;
t1165 = t1121 * t588 + t1198 * t406 - t1199 * t407;
t1163 = t1119 * t588 + t1194 * t406 - t1195 * t407;
t1155 = t1189 * t628;
t1154 = t1189 * t629;
t1153 = t1188 * t628;
t1152 = t1188 * t629;
t547 = rSges(7,2) * t985;
t1151 = t1187 * t628 + t547;
t551 = rSges(7,2) * t980;
t1150 = t1187 * t629 + t551;
t1022 = Icges(5,6) * t629;
t384 = Icges(5,4) * t985 - Icges(5,2) * t986 - t1022;
t618 = Icges(5,4) * t631;
t526 = Icges(5,1) * t630 + t618;
t1149 = -t526 * t628 - t384;
t788 = -Icges(5,2) * t630 + t618;
t385 = Icges(5,6) * t628 + t629 * t788;
t1148 = -t526 * t629 - t385;
t1037 = Icges(5,4) * t630;
t527 = Icges(5,1) * t631 - t1037;
t387 = Icges(5,5) * t628 + t527 * t629;
t524 = Icges(5,2) * t631 + t1037;
t1147 = -t524 * t629 + t387;
t1146 = -t1207 * t630 + t1230;
t1145 = t1208 * t630 + t1228;
t620 = t630 * rSges(7,2);
t1116 = -t1162 * t972 - t1172 * t971 - t620;
t1144 = t526 + t788;
t1096 = g(1) * t629 + g(2) * t628;
t1089 = t630 * t1096;
t912 = qJD(3) * t629;
t858 = -t374 * t603 + t375 * t912 + qJD(2);
t711 = t509 * t464 + t466 * t510 + t858;
t67 = -t406 * t963 + t407 * t962 + t579 + t711;
t555 = pkin(9) * t980;
t465 = -pkin(4) * t981 + t555;
t832 = qJD(1) * t465 - t1113 * t509;
t553 = pkin(9) * t985;
t463 = -pkin(4) * t986 + t553;
t836 = t509 * t463 + t465 * t510;
t846 = t521 + t634;
t815 = t375 + t846;
t744 = t466 + t815;
t944 = -t509 * t530 - t558;
t683 = t744 * qJD(1) + t944;
t84 = t1185 + t683;
t842 = t84 * t939;
t1117 = -t510 * t530 - t826;
t665 = qJD(1) * t745 + t1117;
t83 = -t407 * t939 + t1204 + t665;
t843 = t83 * t963;
t844 = t67 * t962;
t1143 = -t67 * (t1150 * t407 + t1151 * t406 + t631 * t904 - t864 * t963 + t836) - t84 * (t1116 * t406 + t1150 * t588 - t579 * t628 + t909 * t962 + t832) - (t630 * t843 + (-t628 * t844 - t629 * t842) * t631) * qJD(5) - (-t1162 * t638 - t1172 * t641 - pkin(4)) * t1089;
t261 = t472 * rSges(6,1) - t471 * rSges(6,2) + rSges(6,3) * t981;
t101 = t261 * t588 - t406 * t426 + t683;
t895 = rSges(6,1) * t974;
t894 = rSges(6,2) * t975;
t934 = rSges(6,3) * t985 + t628 * t894;
t351 = -t628 * t895 + t934;
t933 = rSges(6,3) * t980 + t629 * t894;
t353 = -t629 * t895 + t933;
t619 = t630 * rSges(6,3);
t428 = rSges(6,1) * t971 - rSges(6,2) * t972 + t619;
t865 = t628 * t908;
t95 = -t259 * t406 + t261 * t407 + t711;
t1142 = -t101 * (t261 * t909 + t588 * t353 - t406 * t428 - t426 * t864 + t832) - t95 * (-t259 * t864 - t261 * t865 + t406 * t351 + t353 * t407 + t836) - g(1) * (t555 + t933) - g(2) * (t553 + t934);
t1141 = t1186 * t630;
t1140 = t1190 * t588 + t1192 * t406 - t1217 * t407;
t1031 = Icges(5,5) * t629;
t542 = Icges(5,4) * t986;
t386 = Icges(5,1) * t985 - t1031 - t542;
t997 = t384 * t630;
t756 = -t386 * t631 + t997;
t712 = t756 * t628;
t523 = Icges(5,5) * t631 - Icges(5,6) * t630;
t383 = Icges(5,3) * t628 + t523 * t629;
t954 = t628 * t383 + t387 * t980;
t1139 = -t712 - t954;
t1108 = t1194 * t629 + t1195 * t628;
t1138 = -t1194 * t628 + t1195 * t629;
t1107 = t1196 * t629 + t1197 * t628;
t1137 = -t1196 * t628 + t1197 * t629;
t1106 = t1198 * t629 + t1199 * t628;
t1136 = -t1198 * t628 + t1199 * t629;
t479 = t530 * t917;
t827 = -qJD(1) * t463 - t1113 * t510;
t1135 = -t259 * t909 + t351 * t588 + t407 * t428 + t479 - t827 + (-t865 + t917) * t426;
t804 = rSges(6,1) * t220 - rSges(6,2) * t219;
t145 = rSges(6,3) * t702 + t804;
t801 = -rSges(6,1) * t638 - rSges(6,2) * t641;
t272 = t802 * t973 + (rSges(6,3) * t635 + qJD(5) * t801) * t630;
t41 = -t145 * t588 + t188 * t426 + t259 * t486 - t272 * t407 + t654;
t1134 = (qJD(1) * t101 + t41) * t629;
t1132 = t67 * t963 + t842;
t1100 = t631 * rSges(5,1) - rSges(5,2) * t630;
t1115 = -rSges(5,2) * t986 - t629 * rSges(5,3);
t389 = t896 + t1115;
t610 = t628 * rSges(5,3);
t390 = rSges(5,1) * t980 - rSges(5,2) * t981 + t610;
t123 = t389 * t509 + t390 * t510 + t858;
t528 = rSges(5,1) * t630 + rSges(5,2) * t631;
t704 = -t510 * t528 - t826;
t746 = -t389 + t816;
t158 = qJD(1) * t746 + t704;
t159 = -t509 * t528 - t558 + (t390 + t815) * qJD(1);
t461 = t528 * t628;
t462 = t528 * t629;
t1131 = -t158 * (qJD(1) * t461 - t1100 * t510) - t123 * (-t509 * t461 - t462 * t510) - t159 * (-qJD(1) * t462 - t1100 * t509);
t1122 = qJD(1) * t84 + t25;
t508 = qJD(1) * t520;
t1118 = qJD(1) * t374 - t508;
t611 = t628 * rSges(4,3);
t977 = t629 * t642;
t979 = t629 * t639;
t409 = rSges(4,1) * t977 - rSges(4,2) * t979 + t611;
t329 = t409 + t846;
t517 = t629 * rSges(3,1) - rSges(3,2) * t628;
t502 = t517 + t634;
t1112 = (t1227 * t630 - t1188) * t588 + (-t1182 * t470 + t1216 + t433 - t436) * t407 + (t1182 * t472 - t1029 - t1215 + t438) * t406;
t1111 = (-Icges(7,1) * t975 + t791 * t630 + t1206 + t583) * t588 + (t1183 * t469 + t1219 - t434 + t437) * t407 + (-t1183 * t471 - t1036 + t1218 + t435) * t406;
t1110 = t1179 * t630;
t632 = Icges(4,4) * t642;
t789 = -Icges(4,2) * t639 + t632;
t564 = Icges(4,1) * t639 + t632;
t1097 = t1119 * t486 + t1184 * t588;
t1019 = Icges(4,3) * t629;
t982 = t628 * t642;
t984 = t628 * t639;
t400 = Icges(4,5) * t982 - Icges(4,6) * t984 - t1019;
t1032 = Icges(4,5) * t629;
t576 = Icges(4,4) * t984;
t404 = Icges(4,1) * t982 - t1032 - t576;
t1023 = Icges(4,6) * t629;
t402 = Icges(4,4) * t982 - Icges(4,2) * t984 - t1023;
t994 = t402 * t639;
t754 = -t404 * t642 + t994;
t167 = -t400 * t629 - t628 * t754;
t1018 = Icges(5,3) * t629;
t522 = Icges(5,5) * t630 + Icges(5,6) * t631;
t990 = t522 * t629;
t996 = t385 * t630;
t1095 = -t635 * t990 + (-t387 * t631 - t523 * t628 + t1018 + t996) * qJD(1);
t920 = qJD(1) * t383;
t991 = t522 * t628;
t1094 = qJD(1) * t756 - t635 * t991 + t920;
t561 = Icges(4,5) * t642 - Icges(4,6) * t639;
t560 = Icges(4,5) * t639 + Icges(4,6) * t642;
t708 = qJD(3) * t560;
t1038 = Icges(4,4) * t639;
t565 = Icges(4,1) * t642 - t1038;
t405 = Icges(4,5) * t628 + t565 * t629;
t403 = Icges(4,6) * t628 + t629 * t789;
t993 = t403 * t639;
t753 = -t405 * t642 + t993;
t1093 = -t629 * t708 + (-t561 * t628 + t1019 + t753) * qJD(1);
t401 = Icges(4,3) * t628 + t561 * t629;
t919 = qJD(1) * t401;
t1092 = qJD(1) * t754 - t628 * t708 + t919;
t749 = t524 * t630 - t526 * t631;
t1091 = qJD(1) * t749 + t523 * t635;
t562 = Icges(4,2) * t642 + t1038;
t747 = t562 * t639 - t564 * t642;
t1090 = t747 * qJD(1) + t561 * qJD(3);
t941 = -Icges(4,2) * t982 + t404 - t576;
t943 = t564 * t628 + t402;
t1088 = -t639 * t941 - t642 * t943;
t1085 = qJD(1) * t1144 + t509 * t1147 - t510 * (-Icges(5,2) * t985 + t386 - t542);
t1084 = t187 / 0.2e1;
t1083 = t188 / 0.2e1;
t1082 = t378 / 0.2e1;
t1081 = t379 / 0.2e1;
t1080 = -t406 / 0.2e1;
t1079 = t406 / 0.2e1;
t1078 = -t407 / 0.2e1;
t1077 = t407 / 0.2e1;
t1075 = t504 / 0.2e1;
t1074 = t505 / 0.2e1;
t1073 = -t509 / 0.2e1;
t1072 = t509 / 0.2e1;
t1071 = -t510 / 0.2e1;
t1070 = t510 / 0.2e1;
t1069 = -t588 / 0.2e1;
t1068 = t588 / 0.2e1;
t1067 = t628 / 0.2e1;
t1066 = -t629 / 0.2e1;
t1062 = -rSges(6,3) - pkin(9);
t1060 = g(1) * t628;
t1056 = -qJD(1) / 0.2e1;
t1055 = qJD(1) / 0.2e1;
t1053 = rSges(4,1) * t642;
t1049 = t37 * t407;
t1048 = t38 * t406;
t1047 = t39 * t407;
t1046 = t40 * t406;
t1045 = t628 * t84;
t1043 = qJDD(1) / 0.2e1;
t1015 = qJD(1) * t67;
t1013 = qJD(1) * t95;
t1012 = t101 * t628;
t1007 = t115 * t188;
t1006 = t116 * t187;
t1005 = t117 * t188;
t1004 = t118 * t187;
t1001 = t159 * t629;
t923 = rSges(4,2) * t984 + t629 * rSges(4,3);
t408 = rSges(4,1) * t982 - t923;
t814 = -t408 + t847;
t569 = rSges(4,1) * t639 + rSges(4,2) * t642;
t867 = t569 * t912;
t204 = qJD(1) * t814 - t867;
t1000 = t204 * t628;
t999 = t204 * t629;
t205 = qJD(1) * t329 - t569 * t603;
t488 = t569 * t629;
t998 = t205 * t488;
t989 = t524 * t635;
t988 = t560 * t628;
t987 = t560 * t629;
t961 = -t261 - t466;
t959 = -t272 - t478;
t957 = t1162 * t470 - t1172 * t469;
t956 = t1162 * t472 - t1172 * t471;
t382 = Icges(5,5) * t985 - Icges(5,6) * t986 - t1018;
t955 = -t628 * t382 - t386 * t980;
t952 = -t628 * t374 + t629 * t375;
t951 = t628 * t389 + t629 * t390;
t950 = -t628 * t400 - t404 * t977;
t949 = t628 * t401 + t405 * t977;
t945 = t628 * t464 + t629 * t466;
t942 = -t564 * t629 - t403;
t940 = -t562 * t629 + t405;
t938 = -t426 - t530;
t931 = rSges(5,2) * t871 + rSges(5,3) * t916;
t929 = t547 + t553;
t928 = t551 + t555;
t913 = qJD(1) * t639;
t927 = rSges(4,2) * t628 * t913 + rSges(4,3) * t916;
t925 = -t562 + t565;
t924 = t564 + t789;
t918 = qJD(1) * t561;
t910 = qJD(3) * t642;
t202 = -t628 * t749 - t990;
t903 = t202 * qJD(1);
t293 = -t628 * t747 - t987;
t902 = t293 * qJD(1);
t890 = pkin(3) * t910;
t231 = rSges(5,1) * t701 - rSges(5,2) * t886 + t931;
t721 = t528 * t635;
t232 = -t628 * t721 + (t1100 * t629 + t610) * qJD(1);
t883 = t629 * t231 + t628 * t232 + t389 * t916;
t882 = t218 * rSges(6,1) + t217 * rSges(6,2) + rSges(6,3) * t886;
t880 = t629 * t278 + t628 * t279 + t464 * t916;
t879 = -t466 - t962;
t878 = -t478 - t960;
t877 = t629 * t301 + t628 * t302 - t374 * t916;
t876 = t917 * t939 + t479;
t875 = -t530 - t939;
t874 = t513 + t926;
t872 = t1062 * t630;
t869 = t629 * t913;
t859 = t973 / 0.2e1;
t857 = -pkin(2) - t1053;
t856 = -pkin(4) - t1052;
t855 = t917 / 0.2e1;
t854 = t916 / 0.2e1;
t853 = -t603 / 0.2e1;
t852 = t603 / 0.2e1;
t851 = -t912 / 0.2e1;
t850 = t912 / 0.2e1;
t703 = -t528 - t1061;
t845 = -t530 - t1061;
t841 = t639 * (-t628 ^ 2 - t629 ^ 2);
t840 = (-t628 * t788 + t1022) * qJD(1) + t1147 * t635;
t839 = qJD(1) * t385 + t386 * t635 - t628 * t989;
t838 = (-t527 * t628 + t1031) * qJD(1) + t1148 * t635;
t837 = qJD(1) * t387 + t1149 * t635;
t371 = t405 * t982;
t835 = t401 * t629 - t371;
t834 = -t382 + t996;
t833 = -t400 + t993;
t831 = t1144 * t635;
t830 = t527 * t635 - t989;
t824 = -t259 * t628 + t629 * t261 + t945;
t812 = -t426 + t845;
t477 = t1100 * t635;
t811 = -t477 - t890;
t810 = -t478 - t890;
t330 = t387 * t985;
t808 = t385 * t986 - t330;
t807 = t634 + t829;
t572 = rSges(2,1) * t643 - rSges(2,2) * t640;
t570 = rSges(2,1) * t640 + rSges(2,2) * t643;
t571 = -rSges(4,2) * t639 + t1053;
t234 = t403 * t642 + t405 * t639;
t709 = qJD(3) * t562;
t289 = -t629 * t709 + (-t628 * t789 + t1023) * qJD(1);
t710 = qJD(3) * t564;
t291 = -t629 * t710 + (-t565 * t628 + t1032) * qJD(1);
t659 = -qJD(3) * t234 - t289 * t639 + t291 * t642 + t919;
t233 = t402 * t642 + t404 * t639;
t290 = qJD(1) * t403 - t628 * t709;
t292 = qJD(1) * t405 - t628 * t710;
t660 = qJD(1) * t400 - qJD(3) * t233 - t290 * t639 + t292 * t642;
t796 = t628 * (t1093 * t628 + t659 * t629) - t629 * (t1092 * t628 + t660 * t629);
t795 = t628 * (-t1093 * t629 + t659 * t628) - t629 * (-t1092 * t629 + t660 * t628);
t100 = t1205 + t665;
t779 = t100 * t629 + t1012;
t766 = -t158 * t629 - t159 * t628;
t168 = -t403 * t984 - t835;
t765 = -t167 * t629 + t168 * t628;
t169 = -t402 * t979 - t950;
t170 = -t403 * t979 + t949;
t764 = -t169 * t629 + t170 * t628;
t763 = -t205 * t628 - t999;
t758 = -t259 * t629 - t261 * t628;
t297 = -rSges(4,2) * t629 * t910 + (-t642 * t917 - t866) * rSges(4,1) + t927;
t487 = t569 * t628;
t298 = -qJD(3) * t487 + (t571 * t629 + t611) * qJD(1);
t757 = t297 * t629 + t298 * t628;
t185 = t384 * t631 + t386 * t630;
t752 = t408 * t628 + t409 * t629;
t748 = t562 * t642 + t564 * t639;
t743 = t845 - t939;
t742 = -t272 + t810;
t741 = -t624 - t1113;
t737 = t872 - t623;
t143 = -rSges(6,3) * t871 + t882;
t736 = t629 * t143 + t628 * t145 - t259 * t916 + t880;
t735 = -t628 * t963 + t629 * t962 + t945;
t720 = t810 - t960;
t719 = t807 + t466;
t706 = t301 * t912 + t302 * t603 - t504 * t374 - t505 * t375 + qJDD(2);
t705 = t1113 - t1116;
t698 = t428 + t1113;
t692 = qJD(1) * t523 - t509 * t990 + t510 * t991;
t691 = t83 * t939 - t844;
t690 = -t639 * t940 + t642 * t942;
t687 = -t628 * t965 + t629 * t966 - t916 * t963 + t880;
t681 = -pkin(4) * t887 + t514 - t826;
t674 = (-t639 * t924 + t642 * t925) * qJD(1);
t673 = -qJD(1) * t464 + t1117 + t1118;
t666 = t1148 * t509 - t1149 * t510 + (-t524 + t527) * qJD(1);
t664 = t510 * t278 + t509 * t279 + t378 * t464 - t379 * t466 + t706;
t663 = qJD(1) * t382 - t630 * t839 + t631 * t837;
t662 = -t630 * t840 + t631 * t838 + t920;
t661 = qJD(1) * t522 - t630 * t831 + t631 * t830;
t533 = t789 * qJD(3);
t534 = t565 * qJD(3);
t658 = qJD(1) * t560 - qJD(3) * t748 - t533 * t639 + t534 * t642;
t656 = t1116 * t407 - t1151 * t588 - t579 * t629 + t865 * t939 + t827;
t650 = -t1085 * t630 + t666 * t631;
t110 = t630 * t837 + t631 * t839;
t111 = t630 * t838 + t631 * t840;
t112 = t1091 * t628 + t661 * t629;
t113 = -t1091 * t629 + t661 * t628;
t162 = -t382 * t629 - t712;
t163 = -t383 * t629 - t808;
t164 = -t384 * t981 - t955;
t165 = -t385 * t981 + t954;
t186 = t385 * t631 + t387 * t630;
t203 = -t629 * t749 + t991;
t74 = t1094 * t628 + t663 * t629;
t75 = t1095 * t628 + t662 * t629;
t76 = -t1094 * t629 + t663 * t628;
t77 = -t1095 * t629 + t662 * t628;
t93 = -t162 * t510 + t163 * t509 + t903;
t192 = t203 * qJD(1);
t94 = -t164 * t510 + t165 * t509 + t192;
t648 = (qJD(1) * t1106 + t1168 * t628 - t1169 * t629) * t1078 + (qJD(1) * t1107 + t1170 * t628 - t1171 * t629) * t1079 - t1138 * t486 / 0.2e1 + (qJD(1) * t113 + qJDD(1) * t202 + t162 * t379 + t163 * t378 + t509 * t77 - t510 * t76 + t1173) * t1066 - t1163 * t909 / 0.2e1 + (t94 + t1164) * t854 + (t93 + t1165) * t855 + (qJD(1) * t1108 + t1166 * t628 - t1167 * t629) * t1068 + (-t110 * t629 + t111 * t628 + (t185 * t628 + t186 * t629) * qJD(1)) * t1055 - t1136 * t188 / 0.2e1 - t1137 * t187 / 0.2e1 + (-t162 * t629 + t163 * t628) * t1081 + (-t164 * t629 + t165 * t628) * t1082 + ((qJD(5) * t1108 - t1186) * t631 + ((t1145 * t641 + t1146 * t638 + t1190) * t588 + (t1153 * t641 - t1155 * t638 - t1217) * t407 + (-t1152 * t641 + t1154 * t638 + t1192) * t406 + t1119 * qJD(5)) * t630) * t1069 + ((t1120 * t630 + t1197 * t985) * qJD(5) + ((qJD(5) * t1196 + t1140) * t631 + t1141) * t629 + (t1145 * t472 + t1146 * t471) * t588 + (t1153 * t472 - t1155 * t471) * t407 + (-t1152 * t472 + t1154 * t471) * t406) * t1080 + ((t1121 * t630 + t1198 * t980) * qJD(5) + ((qJD(5) * t1199 + t1140) * t631 + t1141) * t628 + (t1145 * t470 + t1146 * t469) * t588 + (t1153 * t470 - t1155 * t469) * t407 + (-t1152 * t470 + t1154 * t469) * t406) * t1077 + (-t185 * t629 + t186 * t628) * t1043 + (qJD(1) * t112 + qJDD(1) * t203 + t164 * t379 + t165 * t378 + t509 * t75 - t510 * t74 + t1174) * t1067 + (t628 * t650 - t629 * t692) * t1070 + (t628 * t77 - t629 * t76 + (t162 * t628 + t163 * t629) * qJD(1)) * t1071 + (t628 * t75 - t629 * t74 + (t164 * t628 + t165 * t629) * qJD(1)) * t1072 + (t628 * t692 + t629 * t650) * t1073 + (t1085 * t631 + t666 * t630) * t1056 - (t1164 * t629 + t1165 * t628) * t908 / 0.2e1;
t536 = t571 * qJD(3);
t497 = t801 * t630;
t324 = -rSges(6,1) * t471 - rSges(6,2) * t472;
t319 = -rSges(6,1) * t469 - rSges(6,2) * t470;
t294 = -t629 * t747 + t988;
t282 = t294 * qJD(1);
t201 = qJD(3) * t752 + qJD(2);
t150 = -t1090 * t629 + t658 * t628;
t149 = t1090 * t628 + t658 * t629;
t122 = qJD(1) * t297 + qJDD(1) * t409 - t504 * t569 - t536 * t603 + t724;
t121 = -t897 - t536 * t912 + t505 * t569 + (-t298 - t503) * qJD(1) + t814 * qJDD(1);
t120 = -qJD(3) * t753 + t289 * t642 + t291 * t639;
t119 = -qJD(3) * t754 + t290 * t642 + t292 * t639;
t114 = qJD(3) * t757 + t408 * t504 - t409 * t505 + qJDD(2);
t99 = qJD(3) * t764 + t282;
t98 = qJD(3) * t765 + t902;
t92 = qJD(1) * t231 + qJDD(1) * t390 - t378 * t528 - t477 * t509 + t657;
t91 = t379 * t528 - t477 * t510 + (-t232 + t958) * qJD(1) + t746 * qJDD(1) + t688;
t66 = t231 * t510 + t232 * t509 + t378 * t389 - t379 * t390 + t706;
t42 = t143 * t588 - t187 * t426 + t261 * t486 - t272 * t406 + t655;
t28 = t143 * t407 + t145 * t406 - t187 * t259 - t188 * t261 + t664;
t19 = qJD(6) * t699 + qJDD(6) * t975 - t187 * t963 - t188 * t962 - t406 * t965 + t407 * t966 + t664;
t1 = [(t113 + t110 + t94) * t1071 + t1097 + (t149 + t120) * t852 + (t112 + t111) * t1072 + (t119 + t150 + t99) * t851 + (t192 + (t163 + (t383 + t997) * t629 + t808 + t955) * t510 + (-t629 * t834 - t1139 + t162) * t509) * t1070 + (t1209 * t630 + t1158 + t1161 + t1195) * t406 * t1069 + ((-t721 - t891) * t1001 + (t92 - g(2)) * (t390 + t807) + (t91 - g(1)) * (-t1115 + t1211) + (rSges(5,1) * t889 + rSges(5,2) * t888 + t926 + (-t610 - t634 + (-t624 - t1100) * t629) * qJD(1)) * t158 + (t158 - t704 - t1118 + t931 + (t1057 + t389 + t1211) * qJD(1)) * t159) * m(5) + t1048 / 0.2e1 - t1049 / 0.2e1 + t1006 / 0.2e1 + t1007 / 0.2e1 + t1004 / 0.2e1 + t1005 / 0.2e1 + (-(-t867 - t204 - t508 + (-t408 - t1057) * qJD(1)) * t205 + t205 * (t599 + t927) + (t1000 * t569 - t998) * qJD(3) + ((-t204 * t643 - t205 * t640) * pkin(1) + (-pkin(2) - t571) * t999 + (t204 * (-rSges(4,3) - pkin(7)) + t205 * t857) * t628) * qJD(1) + (-g(2) + t122) * t329 + (-g(1) + t121) * (t628 * t857 - t1057 + t616 + t923)) * m(4) - m(2) * (-g(1) * t570 + g(2) * t572) + ((-t516 * t646 - g(2) + t817) * t502 + (-t897 + (-0.2e1 * t517 - t634 + t502) * t646 - g(1)) * t501) * m(3) + (t203 + t186) * t1082 + (t294 + t234) * t1075 + t1120 * t1084 + t1121 * t1083 + (-t737 * t1060 + t100 * (-t804 + t874) + t101 * (t681 + t882) + (t41 * t872 + (t100 * t1062 * t635 - t41 * pkin(4)) * t631) * t628 + ((-t100 * t643 - t101 * t640) * pkin(1) + (t100 * (t741 - t619) - t101 * t644) * t629 + (-t624 + t737) * t1012) * qJD(1) - (-qJD(1) * t1057 - t100 + t1205 + t673) * t101 + (-g(2) + t42) * (t719 + t261) + (-g(1) + t41) * (-t803 + t813)) * m(6) + t1201 * t1079 + (t1200 + t1164) * t1078 + t1164 * t1077 + t1046 / 0.2e1 - t1047 / 0.2e1 + (t282 + ((t168 - t371 + (t401 + t994) * t629 + t950) * t629 + t949 * t628) * qJD(3)) * t850 + (-t738 * t1060 + t842 * t407 + t738 * t628 * t25 + (t1063 * t888 + t1180 + t1185 + t874 + t944) * t83 + (t1181 - t1204 - t673 + t681) * t84 + t1203 * (t719 + t962) + t1202 * (t813 + t1191) + (-t84 * t598 + (-t624 + t738) * t1045 + (t744 - t634 + (t741 - t620) * t629) * t83) * qJD(1)) * m(7) + (-qJD(3) * t747 + t533 * t642 + t534 * t639 + t630 * t830 + t631 * t831) * qJD(1) + (m(3) * (t501 ^ 2 + t517 * t502) + m(2) * (t570 ^ 2 + t572 ^ 2) + t748 + t524 * t631 + t526 * t630 + Icges(2,3) + Icges(3,3)) * qJDD(1) + (t293 + t233) * t1074 + (t202 + t185) * t1081 + (t98 - t902 + ((t629 * t833 + t170 - t949) * t629 + (t628 * t833 + t169 + t835) * t628) * qJD(3)) * t853 + (-t903 + (t165 + t1139) * t510 + (t628 * t834 + t164 - t330) * t509 + ((t383 + t756) * t509 + t834 * t510) * t629 + t93) * t1073; m(3) * qJDD(2) + (-m(3) - m(4) - m(5) - m(6) - m(7)) * g(3) + m(4) * t114 + m(5) * t66 + m(6) * t28 + m(7) * t19; t99 * t854 + t98 * t855 + (-t119 * t629 + t120 * t628 + (t233 * t628 + t234 * t629) * qJD(1)) * t1055 + ((t639 * t925 + t642 * t924) * qJD(1) + ((t628 * t940 - t629 * t941) * t642 + (t628 * t942 + t629 * t943) * t639) * qJD(3)) * t1056 + (qJD(1) * t150 + qJD(3) * t795 + qJDD(1) * t293 + t167 * t505 + t168 * t504) * t1066 + (-t233 * t629 + t234 * t628) * t1043 + t648 + (qJD(1) * t149 + qJD(3) * t796 + qJDD(1) * t294 + t169 * t505 + t170 * t504) * t1067 + t765 * t1074 + t764 * t1075 + ((t167 * t628 + t168 * t629) * qJD(1) + t795) * t851 + ((t169 * t628 + t170 * t629) * qJD(1) + t796) * t852 + ((-t912 * t988 - t918) * t629 + (t674 + (t690 * t628 + (-t1088 + t987) * t629) * qJD(3)) * t628) * t850 + ((-t603 * t987 + t918) * t628 + (t674 + (-t1088 * t629 + (t988 + t690) * t628) * qJD(3)) * t629) * t853 + (-t656 * t83 - (-t84 * t869 + ((-t629 * t83 - t1045) * t642 + t67 * t841) * qJD(3)) * pkin(3) - g(1) * (-pkin(3) * t979 + t928) - g(2) * (-pkin(3) * t984 + t929) - g(3) * (t633 + t705) + t83 * t876 + t19 * (t735 + t952) + t67 * (t687 + t877) + (t1122 * t743 + t720 * t83) * t629 + (t24 * t743 + t84 * t720 + (-t375 + t879) * t1015) * t628 + t1143) * m(7) + (-(-t101 * t869 + (-t642 * t779 + t841 * t95) * qJD(3)) * pkin(3) - g(3) * (t633 + t698) - t1096 * (t630 * t856 - t1061) + t28 * (t824 + t952) + t95 * (t736 + t877) + t812 * t1134 + (t42 * t812 + t101 * t742 + (-t375 + t961) * t1013) * t628 + (t629 * t742 + t1135) * t100 + t1142) * m(6) + (-(-t159 * t869 + (t123 * t841 + t642 * t766) * qJD(3)) * pkin(3) - g(3) * (t1100 + t633) - t1096 * t703 + t66 * (t951 + t952) + t123 * (t877 + t883) + (t158 * t811 + (qJD(1) * t159 + t91) * t703) * t629 + (t92 * t703 + t159 * t811 + (t158 * t528 + t123 * (-t375 - t390)) * qJD(1)) * t628 + t1131) * m(5) + (t114 * t752 + t201 * ((t408 * t629 - t409 * t628) * qJD(1) + t757) + t763 * t536 + (-t121 * t629 - t122 * t628 + (-t205 * t629 + t1000) * qJD(1)) * t569 + g(1) * t488 + g(2) * t487 - g(3) * t571 - (t204 * t487 - t998) * qJD(1) - (t201 * (-t487 * t628 - t488 * t629) + t763 * t571) * qJD(3)) * m(4); t648 + (t19 * t735 + t67 * t687 + t1122 * t875 * t629 + (t1015 * t879 + t24 * t875 + t84 * t878) * t628 - g(1) * t928 - g(2) * t929 - g(3) * t705 + (t629 * t878 - t656 + t876) * t83 + t1143) * m(7) + (t28 * t824 + t95 * t736 + t938 * t1134 + (t101 * t959 + t1013 * t961 + t42 * t938) * t628 - g(3) * t698 - t856 * t1089 + (t629 * t959 + t1135) * t100 + t1142) * m(6) + (t66 * t951 + t123 * (-t390 * t917 + t883) + t766 * t477 + (-t92 * t628 - t91 * t629 + (t158 * t628 - t1001) * qJD(1)) * t528 + g(1) * t462 + g(2) * t461 - g(3) * t1100 + t1131) * m(5); t1173 * t986 / 0.2e1 + t1174 * t981 / 0.2e1 + t1163 * t976 / 0.2e1 + (t1107 * t630 - t1120 * t631) * t1084 + (t1106 * t630 - t1121 * t631) * t1083 + (-t1110 * t629 + t1111 * t472 + t1112 * t471) * t1080 + ((t1107 * t635 - t1201) * t631 + (qJD(1) * t1137 + t1120 * t635 + t1170 * t629 + t1171 * t628) * t630) * t1079 + ((t1106 * t635 - t1200) * t631 + (qJD(1) * t1136 + t1121 * t635 + t1168 * t629 + t1169 * t628) * t630) * t1078 + (-t1110 * t628 + t1111 * t470 + t1112 * t469) * t1077 + (t1108 * t630 - t1119 * t631) * t486 / 0.2e1 + (t1179 * t631 + (t1111 * t641 + t1112 * t638) * t630) * t1069 + ((t1108 * t635 - t1184) * t631 + (qJD(1) * t1138 + t1119 * t635 + t1166 * t629 + t1167 * t628) * t630) * t1068 - (t1006 + t1007 + t1048 - t1049 + t1004 + t1005 + t1046 - t1047 + t1097) * t631 / 0.2e1 + ((-t25 * t963 - t83 * t965 - t24 * t962 - t84 * t966 + (-t1132 * t629 + t691 * t628) * t635) * t631 + ((t84 * t962 + t843) * t635 + (qJD(1) * t691 - t19 * t963 - t24 * t939 - t67 * t965 - t84 * t960) * t629 + (qJD(1) * t1132 - t19 * t962 + t25 * t939 - t67 * t966 + t83 * t960) * t628) * t630 - (t470 * t84 + t472 * t83 + t67 * t974) * qJD(6) - (-t83 * t957 + t84 * t956) * t588 - (t67 * t956 + t83 * t935) * t407 - (t67 * t957 + t84 * t935) * t406 - g(1) * t956 - g(2) * t957 + g(3) * t935) * m(7) + (-t100 * (-t319 * t588 - t407 * t497) - t101 * (t324 * t588 - t406 * t497) - t95 * (t319 * t406 + t324 * t407) - g(1) * t324 - g(2) * t319 - g(3) * t497 + (t100 * t145 - t101 * t143 - t41 * t259 - t42 * t261 + (t95 * t758 + (t100 * t628 - t101 * t629) * t426) * t635) * t631 + (t100 * (t259 * t635 + t272 * t628) + t101 * (t261 * t635 - t272 * t629) + t28 * t758 + t95 * (-t143 * t628 + t145 * t629 + t259 * t917 - t261 * t916) + (qJD(1) * t779 + t41 * t628 - t42 * t629) * t426) * t630) * m(6) + t1165 * (t628 * t859 + t630 * t854) + t1164 * (-t871 / 0.2e1 + t629 * t859); (t67 * t884 - t217 * t83 + t219 * t84 + (t19 * t638 + t67 * t905) * t630 + (t406 * t84 + t407 * t83 - g(3)) * t975 + (-t407 * t67 - t588 * t84 + t1202) * t471 + (-t406 * t67 + t588 * t83 + t1203) * t469) * m(7);];
tau  = t1;
