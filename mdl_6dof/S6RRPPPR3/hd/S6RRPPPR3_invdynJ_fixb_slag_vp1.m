% Calculate vector of inverse dynamics joint torques for
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:53
% EndTime: 2019-03-09 08:15:37
% DurationCPUTime: 96.74s
% Computational Cost: add. (21512->1512), mult. (42438->1790), div. (0->0), fcn. (37253->8), ass. (0->711)
t646 = sin(qJ(2));
t1053 = Icges(5,4) * t646;
t621 = Icges(4,5) * t646;
t648 = cos(qJ(2));
t501 = -Icges(4,3) * t648 + t621;
t1213 = -Icges(5,1) * t648 - t1053 + t501;
t1055 = Icges(3,4) * t646;
t509 = Icges(3,2) * t648 + t1055;
t1249 = -t509 + t1213;
t1047 = Icges(4,5) * t648;
t502 = Icges(4,3) * t646 + t1047;
t647 = sin(qJ(1));
t649 = cos(qJ(1));
t349 = -Icges(4,6) * t649 + t502 * t647;
t1010 = t646 * t647;
t1046 = Icges(5,5) * t649;
t1007 = t647 * t648;
t574 = Icges(5,4) * t1007;
t359 = Icges(5,1) * t1010 + t1046 - t574;
t1248 = -t349 - t359;
t1009 = t646 * t649;
t1040 = Icges(4,6) * t647;
t1006 = t648 * t649;
t573 = Icges(4,5) * t1006;
t350 = Icges(4,3) * t1009 + t1040 + t573;
t1052 = Icges(5,4) * t648;
t511 = Icges(5,1) * t646 - t1052;
t360 = -Icges(5,5) * t647 + t511 * t649;
t1247 = -t350 - t360;
t504 = Icges(3,5) * t648 - Icges(3,6) * t646;
t352 = Icges(3,3) * t647 + t504 * t649;
t508 = Icges(4,4) * t648 + Icges(4,6) * t646;
t356 = Icges(4,2) * t647 + t508 * t649;
t1246 = -t352 - t356;
t792 = Icges(4,1) * t648 + t621;
t362 = Icges(4,4) * t647 + t649 * t792;
t516 = Icges(3,1) * t648 - t1055;
t364 = Icges(3,5) * t647 + t516 * t649;
t1236 = -t362 - t364;
t624 = Icges(3,4) * t648;
t788 = -Icges(3,2) * t646 + t624;
t1234 = -t511 + t788;
t1245 = t502 - t1234;
t505 = -Icges(5,2) * t648 + t1053;
t1235 = t505 - t516;
t1243 = t792 - t1235;
t1037 = Icges(5,6) * t649;
t353 = Icges(5,4) * t1010 - Icges(5,2) * t1007 + t1037;
t361 = -Icges(4,4) * t649 + t647 * t792;
t1048 = Icges(3,5) * t649;
t576 = Icges(3,4) * t1010;
t363 = Icges(3,1) * t1007 - t1048 - t576;
t1242 = t353 - t361 - t363;
t513 = Icges(4,1) * t646 - t1047;
t787 = Icges(5,2) * t646 + t1052;
t1240 = Icges(3,1) * t646 + t513 + t624 + t787;
t499 = Icges(5,5) * t646 - Icges(5,6) * t648;
t348 = -Icges(5,3) * t647 + t499 * t649;
t325 = t649 * t348;
t1239 = -t364 * t1007 - t360 * t1010 - t325;
t1226 = (Icges(5,5) + Icges(3,6) - Icges(4,6)) * t648 + (Icges(4,4) + Icges(3,5) + Icges(5,6)) * t646;
t358 = Icges(3,6) * t647 + t649 * t788;
t1237 = t360 - t358;
t1233 = -t352 * t649 - t1239;
t1214 = t499 - t504;
t1232 = -t508 + t1214;
t1038 = Icges(5,6) * t647;
t575 = Icges(5,4) * t1009;
t354 = -Icges(5,2) * t1006 - t1038 + t575;
t1177 = -t1007 * t354 - t1010 * t358 + t1233;
t643 = sin(pkin(9));
t644 = cos(pkin(9));
t414 = -t1009 * t643 - t644 * t647;
t1004 = t649 * t644;
t1012 = t643 * t647;
t415 = t1004 * t646 - t1012;
t214 = Icges(6,5) * t415 + Icges(6,6) * t414 + Icges(6,3) * t1006;
t217 = Icges(6,4) * t415 + Icges(6,2) * t414 + Icges(6,6) * t1006;
t220 = Icges(6,1) * t415 + Icges(6,4) * t414 + Icges(6,5) * t1006;
t412 = t1010 * t643 - t1004;
t1011 = t643 * t649;
t413 = t1010 * t644 + t1011;
t1181 = t350 * t1010 - t412 * t217 + t413 * t220 - t356 * t649 + (t214 + t362) * t1007;
t1124 = t1177 + t1181;
t355 = -Icges(4,2) * t649 + t508 * t647;
t1021 = t355 * t649;
t763 = t349 * t646 + t361 * t648;
t1158 = t647 * t763;
t1032 = Icges(5,3) * t649;
t347 = Icges(5,5) * t1010 - Icges(5,6) * t1007 + t1032;
t760 = t353 * t648 - t359 * t646;
t1195 = t347 * t649 - t647 * t760 - t1021 + t1158;
t212 = Icges(6,5) * t413 - Icges(6,6) * t412 + Icges(6,3) * t1007;
t1034 = Icges(3,3) * t649;
t351 = Icges(3,5) * t1007 - Icges(3,6) * t1010 - t1034;
t1041 = Icges(3,6) * t649;
t357 = Icges(3,4) * t1007 - Icges(3,2) * t1010 - t1041;
t1020 = t357 * t646;
t756 = -t363 * t648 + t1020;
t215 = Icges(6,4) * t413 - Icges(6,2) * t412 + Icges(6,6) * t1007;
t218 = Icges(6,1) * t413 - Icges(6,4) * t412 + Icges(6,5) * t1007;
t770 = -t215 * t412 + t218 * t413;
t1125 = t1007 * t212 - t351 * t649 - t647 * t756 + t1195 + t770;
t1231 = t1006 * t1236 + t1009 * t1247 + t1246 * t647;
t1230 = t357 + t1248;
t1229 = t350 + t1237;
t1227 = t354 + t1236;
t1225 = t1249 * qJD(2);
t1224 = t1240 * qJD(2);
t1217 = t1240 * t648 + t1249 * t646;
t323 = t647 * t355;
t1223 = t1242 * t1006 + t1009 * t1248 - t647 * t351 - t323;
t1169 = -t1006 * t354 - t1009 * t358 - t348 * t647 - t1231;
t71 = t214 * t1006 + t414 * t217 + t415 * t220;
t1122 = t71 + t1169;
t1221 = t354 * t648 + t358 * t646;
t1220 = t1009 * t357 + t347 * t647 + t1223;
t1219 = t1245 * qJD(2);
t1218 = t1243 * qJD(2);
t1216 = -t1240 * t646 + t1249 * t648;
t1141 = t1226 * t649;
t1142 = t1226 * t647;
t1215 = t347 - t351;
t1212 = t1230 * t648 - t1242 * t646;
t1211 = t1227 * t646 + t1229 * t648;
t1210 = -t1225 * t649 + (t1234 * t647 - t1041 - t1046 - t349) * qJD(1);
t1209 = t1225 * t647 + (-t502 * t649 - t1040 - t1237) * qJD(1);
t1208 = -t1224 * t649 + (t1235 * t647 + t1037 + t1048 - t361) * qJD(1);
t1207 = t1224 * t647 + (t505 * t649 - t1038 + t1236) * qJD(1);
t1206 = t1226 * qJD(2);
t1205 = t1236 * t648 + t1247 * t646 + t1221;
t1204 = t756 + t760 - t763;
t640 = pkin(9) + qJ(6);
t608 = cos(t640);
t1005 = t649 * t608;
t607 = sin(t640);
t1008 = t647 * t607;
t343 = t1008 * t646 - t1005;
t344 = t1010 * t608 + t607 * t649;
t816 = t344 * rSges(7,1) - t343 * rSges(7,2);
t167 = -rSges(7,3) * t1007 - t816;
t1072 = rSges(7,3) * t646;
t1075 = rSges(7,2) * t607;
t1079 = rSges(7,1) * t608;
t815 = -t1075 + t1079;
t320 = -t648 * t815 + t1072;
t928 = qJD(6) * t648;
t936 = qJD(2) * t649;
t467 = -t647 * t928 + t936;
t929 = qJD(6) * t646;
t591 = qJD(1) + t929;
t645 = -pkin(8) - qJ(5);
t1003 = qJ(5) + t645;
t598 = pkin(5) * t644 + pkin(4);
t1083 = pkin(4) - t598;
t1146 = t1083 * t648;
t307 = -t1003 * t646 + t1146;
t523 = -pkin(4) * t648 + qJ(5) * t646;
t1091 = pkin(3) * t646;
t518 = pkin(2) * t646 - qJ(3) * t648;
t866 = -t518 - t1091;
t842 = -t523 + t866;
t746 = -t307 + t842;
t1203 = -t167 * t591 + t320 * t467 - t746 * t936;
t782 = Icges(6,5) * t644 - Icges(6,6) * t643;
t336 = Icges(6,3) * t646 - t648 * t782;
t786 = Icges(6,4) * t644 - Icges(6,2) * t643;
t685 = -Icges(6,6) * t646 + t648 * t786;
t790 = Icges(6,1) * t644 - Icges(6,4) * t643;
t687 = -Icges(6,5) * t646 + t648 * t790;
t1161 = t1007 * t336 + t1217 * t647 + t412 * t685 - t413 * t687 - t1141;
t1160 = t1006 * t336 + t1217 * t649 - t414 * t685 - t415 * t687 + t1142;
t1080 = rSges(6,1) * t644;
t630 = t648 * rSges(6,3);
t616 = t648 * qJ(5);
t634 = t646 * pkin(4);
t1155 = t634 + t616;
t615 = t646 * qJ(3);
t637 = t648 * pkin(2);
t1154 = t637 + t615;
t636 = t648 * pkin(3);
t893 = t636 + t1154;
t843 = -t1155 - t893;
t1202 = t646 * t1080 + t630 - t843;
t1157 = -t413 * rSges(6,1) + t412 * rSges(6,2);
t222 = -rSges(6,3) * t1007 + t1157;
t1201 = qJD(1) * t222;
t1200 = qJD(1) * t1226 + t1216 * qJD(2) + t1218 * t648 + t1219 * t646;
t1199 = (-t348 - t1246) * qJD(1);
t1198 = t1217 * qJD(1) + qJD(2) * t1232;
t1159 = -t414 * t215 - t415 * t218;
t70 = t212 * t1006 - t1159;
t1065 = t649 * t70;
t1197 = t1122 * t647 + t1220 * t649 - t1065;
t1196 = t1124 * t647 - t1125 * t649;
t565 = pkin(5) * t1011;
t970 = -t598 * t1010 - t565;
t1194 = -t816 + t970;
t626 = t646 * rSges(5,1);
t1145 = -rSges(5,2) * t648 + t626;
t1193 = t1145 + t893;
t1192 = qJD(2) * t1211 + t1208 * t648 + t1210 * t646 + t1199;
t1150 = qJD(1) * t355;
t1191 = qJD(1) * t1215 + t1212 * qJD(2) + t1207 * t648 + t1209 * t646 - t1150;
t1190 = (Icges(5,2) * t1010 + t1230 + t574) * t649 + (-Icges(4,1) * t1009 + t573 + (t513 - t787) * t649 + t1229) * t647;
t1189 = -t1240 + t1245;
t1188 = t1243 + t1249;
t1187 = qJD(1) * t1204 - t1206 * t647 + t1199;
t1186 = -t1150 - t1206 * t649 + (t1214 * t647 + t1032 + t1034 + t1205) * qJD(1);
t938 = qJD(2) * t647;
t466 = t649 * t928 + t938;
t157 = Icges(7,5) * t344 - Icges(7,6) * t343 + Icges(7,3) * t1007;
t329 = Icges(7,4) * t344;
t160 = -Icges(7,2) * t343 + Icges(7,6) * t1007 + t329;
t328 = Icges(7,4) * t343;
t164 = -Icges(7,1) * t344 - Icges(7,5) * t1007 + t328;
t60 = t1007 * t157 - t160 * t343 - t164 * t344;
t345 = -t1009 * t607 - t608 * t647;
t346 = t1005 * t646 - t1008;
t159 = Icges(7,5) * t346 + Icges(7,6) * t345 + Icges(7,3) * t1006;
t1051 = Icges(7,4) * t346;
t162 = Icges(7,2) * t345 + Icges(7,6) * t1006 + t1051;
t330 = Icges(7,4) * t345;
t165 = Icges(7,1) * t346 + Icges(7,5) * t1006 + t330;
t61 = t159 * t1007 - t343 * t162 + t344 * t165;
t1013 = t608 * t648;
t1044 = Icges(7,5) * t646;
t1014 = t607 * t648;
t529 = Icges(7,4) * t1014;
t317 = -Icges(7,1) * t1013 + t1044 + t529;
t781 = Icges(7,5) * t608 - Icges(7,6) * t607;
t682 = -Icges(7,3) * t646 + t648 * t781;
t1049 = Icges(7,4) * t608;
t785 = -Icges(7,2) * t607 + t1049;
t684 = -Icges(7,6) * t646 + t648 * t785;
t98 = -t1007 * t682 + t317 * t344 + t343 * t684;
t13 = t466 * t61 - t467 * t60 + t591 * t98;
t62 = t157 * t1006 + t345 * t160 - t164 * t346;
t63 = t159 * t1006 + t345 * t162 + t346 * t165;
t99 = -t1006 * t682 + t317 * t346 - t345 * t684;
t14 = t466 * t63 - t467 * t62 + t99 * t591;
t774 = t160 * t607 + t164 * t608;
t66 = t157 * t646 + t648 * t774;
t1185 = t1160 * qJD(1);
t1184 = -Icges(3,2) * t1007 + t1213 * t647 - t1242 + t212 - t576;
t1183 = Icges(5,1) * t1006 - t214 + t575 + (-t501 + t509) * t649 + t1227;
t1182 = t1161 * qJD(1);
t885 = t646 * t936;
t941 = qJD(1) * t648;
t889 = t647 * t941;
t1173 = t885 + t889;
t933 = qJD(4) * t647;
t940 = qJD(1) * t649;
t286 = -pkin(3) * t1173 - qJ(4) * t940 - t933;
t886 = t646 * t938;
t548 = pkin(3) * t886;
t942 = qJD(1) * t647;
t955 = -qJ(4) * t942 - t548;
t287 = (pkin(3) * t941 + qJD(4)) * t649 + t955;
t1028 = qJ(4) * t649;
t482 = pkin(3) * t1007 + t1028;
t564 = qJ(3) * t1007;
t444 = -pkin(2) * t1010 + t564;
t568 = qJ(3) * t1006;
t451 = -pkin(2) * t1009 + t568;
t614 = qJD(3) * t646;
t904 = t444 * t938 + t451 * t936 + t614;
t891 = t646 * t942;
t883 = t648 * t936;
t542 = qJ(3) * t883;
t934 = qJD(3) * t649;
t562 = t646 * t934;
t958 = t542 + t562;
t227 = -pkin(2) * t1173 - qJ(3) * t891 + t958;
t937 = qJD(2) * t648;
t884 = t647 * t937;
t890 = t646 * t940;
t420 = t884 + t890;
t549 = pkin(2) * t886;
t879 = t647 * t614;
t830 = -t549 + t879;
t888 = t648 * t940;
t229 = pkin(2) * t888 + qJ(3) * t420 + t830;
t448 = t1154 * t647;
t906 = t649 * t227 + t647 * t229 + t448 * t940;
t1176 = t649 * t286 + t647 * t287 + t482 * t940 - t904 + t906;
t1175 = t1215 + t1221;
t922 = qJD(2) * qJD(3);
t1174 = qJDD(3) * t646 + t648 * t922;
t460 = t518 * t942;
t563 = t648 * t934;
t1172 = t460 - t563;
t823 = t648 * rSges(4,1) + t646 * rSges(4,3);
t960 = -t1154 - t823;
t1168 = -t1021 + t1231;
t1167 = qJD(2) * t1196 + t1182;
t1166 = qJD(2) * t1197 + t1185;
t278 = qJD(1) * t414 - t643 * t884;
t279 = qJD(1) * t415 + t644 * t884;
t419 = -t886 + t888;
t115 = Icges(6,5) * t279 + Icges(6,6) * t278 + Icges(6,3) * t419;
t117 = Icges(6,4) * t279 + Icges(6,2) * t278 + Icges(6,6) * t419;
t119 = Icges(6,1) * t279 + Icges(6,4) * t278 + Icges(6,5) * t419;
t769 = t215 * t643 - t218 * t644;
t1165 = (-t117 * t643 + t119 * t644 - t1209) * t648 + (-t115 + t1207) * t646 + (-t212 * t648 + t646 * t769 + t1204) * qJD(2);
t276 = qJD(1) * t412 - t643 * t883;
t277 = -qJD(1) * t413 + t644 * t883;
t114 = Icges(6,5) * t277 + Icges(6,6) * t276 - Icges(6,3) * t1173;
t116 = Icges(6,4) * t277 + Icges(6,2) * t276 - Icges(6,6) * t1173;
t118 = Icges(6,1) * t277 + Icges(6,4) * t276 - Icges(6,5) * t1173;
t768 = t217 * t643 - t220 * t644;
t1164 = (t116 * t643 - t118 * t644 - t1210) * t648 + (t114 + t1208) * t646 + (t214 * t648 - t646 * t768 - t1205) * qJD(2);
t335 = Icges(6,3) * t648 + t646 * t782;
t308 = t335 * qJD(2);
t337 = Icges(6,6) * t648 + t646 * t786;
t309 = t337 * qJD(2);
t339 = Icges(6,5) * t648 + t646 * t790;
t310 = t339 * qJD(2);
t1163 = t1006 * t308 - t336 * t1173 - t1198 * t647 + t1200 * t649 - t276 * t685 - t277 * t687 + t309 * t414 + t310 * t415;
t1162 = t1007 * t308 + t1198 * t649 + t1200 * t647 - t278 * t685 - t279 * t687 - t309 * t412 + t310 * t413 + t336 * t419;
t1123 = t70 - t1220;
t1120 = t214 * t646 + t648 * t768 - t1211;
t1121 = -t212 * t646 - t648 * t769 - t1212;
t638 = t649 * pkin(7);
t532 = pkin(1) * t647 - t638;
t494 = qJD(1) * t532;
t1156 = -qJD(1) * t448 - t494;
t533 = t649 * pkin(1) + t647 * pkin(7);
t550 = pkin(4) * t883;
t930 = qJD(5) * t649;
t560 = t648 * t930;
t226 = -pkin(4) * t891 - qJ(5) * t1173 + t550 + t560;
t1071 = pkin(4) * qJD(2);
t541 = qJ(5) * t886;
t228 = pkin(4) * t890 - t541 + (qJ(5) * t940 + (qJD(5) + t1071) * t647) * t648;
t443 = t1155 * t647;
t592 = pkin(4) * t1007;
t447 = -qJ(5) * t1010 + t592;
t594 = pkin(4) * t1006;
t454 = -qJ(5) * t1009 + t594;
t613 = qJD(5) * t648;
t1153 = t649 * t226 + t647 * t228 + t443 * t940 - t447 * t938 - t454 * t936 + t1176 - t613;
t628 = t647 * rSges(4,2);
t387 = rSges(4,1) * t1006 + rSges(4,3) * t1009 + t628;
t455 = pkin(2) * t1006 + qJ(3) * t1009;
t845 = t455 + t533;
t231 = t845 + t387;
t488 = t533 * qJD(1);
t996 = -t229 - t488;
t1151 = t996 - t287;
t580 = rSges(5,2) * t1007;
t735 = -rSges(5,1) * t1010 - rSges(5,3) * t649;
t383 = -t580 - t735;
t1149 = qJD(1) * t383;
t1148 = t1003 * t648;
t1147 = t1083 * t646;
t1087 = g(2) * t647;
t836 = g(1) * t649 + t1087;
t1140 = t1232 * qJD(1);
t943 = qJD(1) * t646;
t853 = qJD(6) + t943;
t1139 = t647 * t853 - t883;
t1138 = t649 * t853 + t884;
t382 = (rSges(7,1) * t607 + rSges(7,2) * t608) * t648;
t629 = t648 * rSges(7,3);
t175 = qJD(6) * t382 + (t646 * t815 + t629) * qJD(2);
t698 = -t1147 - t1148;
t274 = t698 * qJD(2);
t612 = qJD(5) * t646;
t399 = qJD(2) * t1155 + t612;
t935 = qJD(3) * t648;
t400 = qJD(2) * t1154 - t935;
t837 = -pkin(3) * t937 - t400;
t740 = -t399 + t837;
t530 = t646 * t598;
t859 = -t645 * t648 + t530;
t1128 = -qJD(2) * (-t859 + t1155 + t843) - t175 - t274 + t740;
t1076 = rSges(6,2) * t643;
t818 = -t1076 + t1080;
t318 = (t646 * t818 + t630) * qJD(2);
t1127 = -qJD(2) * (t1076 * t646 - t1202) - t318 + t740;
t477 = t1145 * qJD(2);
t1126 = qJD(2) * t1193 - t477 + t837;
t765 = -t643 * t685 + t644 * t687;
t654 = (t335 - t765) * qJD(1) + (-t647 * t768 + t649 * t769) * qJD(2);
t1119 = -t336 * t943 + t648 * t654 + (t1188 * t648 + t1189 * t646) * qJD(1);
t1118 = t1190 * t648 + (t1183 * t647 + t1184 * t649) * t646;
t1117 = (-t1007 * t115 + t117 * t412 + t1187 * t649 - t119 * t413 - t212 * t419 - t215 * t278 - t218 * t279) * t649 + (t1007 * t114 - t116 * t412 + t118 * t413 + t214 * t419 + t217 * t278 + t220 * t279 + t1192 * t647 + (-t1186 + t1191) * t649) * t647;
t1116 = (-t1006 * t115 - t117 * t414 + t1173 * t212 - t119 * t415 + t1191 * t649 - t215 * t276 - t218 * t277) * t649 + (t1006 * t114 + t116 * t414 + t118 * t415 - t214 * t1173 + t217 * t276 + t220 * t277 + t1186 * t647 + (-t1187 + t1192) * t649) * t647;
t1115 = t523 * t942 + t646 * t930 + t1172;
t673 = t466 * (-Icges(7,2) * t346 + t165 + t330) - t467 * (-Icges(7,2) * t344 - t164 - t328) + t591 * (Icges(7,2) * t1013 + t317 + t529);
t367 = (Icges(7,1) * t607 + t1049) * t648;
t1114 = (-Icges(7,1) * t345 + t1051 + t162) * t466 - (Icges(7,1) * t343 + t160 + t329) * t467 + t591 * (-t684 - t367);
t1113 = 0.2e1 * t646;
t1112 = 0.2e1 * t648;
t1111 = m(4) / 0.2e1;
t1110 = m(5) / 0.2e1;
t1109 = m(6) / 0.2e1;
t1108 = m(7) / 0.2e1;
t1107 = -m(6) - m(7);
t1106 = -pkin(2) - pkin(3);
t923 = qJD(1) * qJD(2);
t491 = qJDD(2) * t647 + t649 * t923;
t918 = qJDD(6) * t648;
t256 = -qJD(6) * t1173 + t649 * t918 + t491;
t1105 = t256 / 0.2e1;
t492 = -qJDD(2) * t649 + t647 * t923;
t257 = qJD(6) * t419 + t647 * t918 + t492;
t1104 = t257 / 0.2e1;
t465 = qJD(2) * t928 + qJDD(6) * t646 + qJDD(1);
t1103 = t465 / 0.2e1;
t1102 = -t466 / 0.2e1;
t1101 = t466 / 0.2e1;
t1100 = -t467 / 0.2e1;
t1099 = t467 / 0.2e1;
t1098 = t491 / 0.2e1;
t1097 = t492 / 0.2e1;
t1096 = -t591 / 0.2e1;
t1095 = t591 / 0.2e1;
t1092 = -rSges(4,1) - pkin(2);
t1090 = pkin(5) * t643;
t1089 = g(1) * t647;
t1084 = -pkin(4) - qJ(3);
t1082 = rSges(3,1) * t648;
t1081 = rSges(4,1) * t646;
t1078 = rSges(5,2) * t646;
t1074 = rSges(5,3) * t647;
t1073 = rSges(6,3) * t646;
t743 = t647 * t591;
t153 = -t1138 * t607 - t608 * t743;
t154 = t1138 * t608 - t607 * t743;
t87 = Icges(7,5) * t154 + Icges(7,6) * t153 + Icges(7,3) * t419;
t89 = Icges(7,4) * t154 + Icges(7,2) * t153 + Icges(7,6) * t419;
t91 = Icges(7,1) * t154 + Icges(7,4) * t153 + Icges(7,5) * t419;
t10 = (-qJD(2) * t774 + t87) * t646 + (qJD(2) * t157 + t607 * t89 - t608 * t91 + (t160 * t608 - t164 * t607) * qJD(6)) * t648;
t1070 = t10 * t467;
t773 = t162 * t607 - t165 * t608;
t742 = t649 * t591;
t151 = t1139 * t607 - t608 * t742;
t152 = -t1139 * t608 - t607 * t742;
t86 = Icges(7,5) * t152 + Icges(7,6) * t151 - Icges(7,3) * t1173;
t88 = Icges(7,4) * t152 + Icges(7,2) * t151 - Icges(7,6) * t1173;
t90 = Icges(7,1) * t152 + Icges(7,4) * t151 - Icges(7,5) * t1173;
t11 = (-qJD(2) * t773 + t86) * t646 + (qJD(2) * t159 + t607 * t88 - t608 * t90 + (t162 * t608 + t165 * t607) * qJD(6)) * t648;
t1069 = t11 * t466;
t627 = t647 * rSges(3,3);
t1064 = t66 * t257;
t67 = t159 * t646 + t648 * t773;
t1063 = t67 * t256;
t1061 = -rSges(5,1) - qJ(3);
t1060 = -rSges(4,3) - qJ(3);
t1059 = -rSges(5,3) - qJ(4);
t766 = -t317 * t608 - t607 * t684;
t107 = -t646 * t682 + t648 * t766;
t312 = Icges(7,3) * t648 + t646 * t781;
t365 = (Icges(7,5) * t607 + Icges(7,6) * t608) * t648;
t169 = qJD(2) * t312 + qJD(6) * t365;
t1050 = Icges(7,4) * t607;
t314 = Icges(7,6) * t648 + t646 * t785;
t170 = (Icges(7,2) * t608 + t1050) * t928 + t314 * qJD(2);
t789 = Icges(7,1) * t608 - t1050;
t316 = Icges(7,5) * t648 + t646 * t789;
t171 = qJD(2) * t316 + qJD(6) * t367;
t33 = (-qJD(2) * t766 + t169) * t646 + (-qJD(2) * t682 + t170 * t607 - t171 * t608 + (t317 * t607 - t608 * t684) * qJD(6)) * t648;
t1058 = t107 * t465 + t33 * t591;
t1029 = qJ(4) * t647;
t521 = rSges(3,1) * t646 + rSges(3,2) * t648;
t887 = t521 * t936;
t952 = rSges(3,2) * t1010 + t649 * rSges(3,3);
t385 = rSges(3,1) * t1007 - t952;
t981 = -t385 - t532;
t180 = qJD(1) * t981 - t887;
t1027 = t180 * t647;
t1026 = t180 * t649;
t388 = rSges(3,1) * t1006 - rSges(3,2) * t1009 + t627;
t291 = t388 + t533;
t181 = qJD(1) * t291 - t521 * t938;
t453 = t521 * t649;
t1025 = t181 * t453;
t1002 = t152 * rSges(7,1) + t151 * rSges(7,2);
t224 = (t634 + t1148) * t647 + t970;
t997 = -t167 - t224;
t995 = t277 * rSges(6,1) + t276 * rSges(6,2);
t980 = -t387 - t455;
t979 = t647 * t448 + t649 * t455;
t977 = -t823 * qJD(2) - t400;
t976 = qJD(1) * t451 + t647 * t935;
t975 = -t444 - t447;
t974 = -t448 - t532;
t595 = pkin(3) * t1006;
t483 = t595 - t1029;
t973 = t455 + t483;
t606 = pkin(7) * t940;
t971 = qJD(1) * (-pkin(1) * t942 + t606) + qJDD(1) * t533;
t969 = t598 * t1007 + t645 * t1010;
t968 = t598 * t1006 + t645 * t1009;
t520 = -rSges(4,3) * t648 + t1081;
t961 = -t518 - t520;
t959 = t646 * t1079 + t629;
t957 = rSges(4,2) * t940 + rSges(4,3) * t883;
t956 = rSges(3,2) * t891 + rSges(3,3) * t940;
t449 = rSges(5,1) * t1007 + rSges(5,2) * t1010;
t456 = rSges(5,1) * t1006 + rSges(5,2) * t1009;
t450 = pkin(4) * t1009 + qJ(5) * t1006;
t951 = t647 ^ 2 + t649 ^ 2;
t939 = qJD(2) * t646;
t932 = qJD(4) * t649;
t931 = qJD(5) * t647;
t924 = -m(5) + t1107;
t921 = qJD(4) * qJD(1);
t919 = qJDD(5) * t648;
t917 = -rSges(7,3) + t1106;
t916 = pkin(3) * t1009;
t915 = pkin(5) * t1012;
t914 = qJD(2) ^ 2 * t636;
t913 = t648 * t1080;
t912 = rSges(7,1) * t1013;
t911 = t648 * t1076;
t910 = rSges(7,2) * t1014;
t907 = t649 * t1106;
t905 = t1174 * t649 + t492 * t518;
t168 = t346 * rSges(7,1) + t345 * rSges(7,2) + rSges(7,3) * t1006;
t633 = t649 * rSges(4,2);
t384 = t647 * t823 - t633;
t903 = -t384 + t974;
t838 = rSges(5,1) * t1009 - rSges(5,2) * t1006;
t386 = t838 - t1074;
t902 = -t386 - t973;
t223 = t415 * rSges(6,1) + t414 * rSges(6,2) + rSges(6,3) * t1006;
t901 = qJD(1) * t454 + t976;
t900 = -t482 + t974;
t899 = -t450 - t973;
t897 = t1173 * t645 + t598 * t883;
t896 = t606 + t958;
t895 = rSges(5,1) * t883 + rSges(5,2) * t1173;
t894 = t549 - t955;
t892 = t649 * t1092;
t878 = -pkin(1) - t1082;
t875 = t940 / 0.2e1;
t874 = -t938 / 0.2e1;
t873 = t938 / 0.2e1;
t872 = -t936 / 0.2e1;
t871 = t936 / 0.2e1;
t869 = -rSges(6,3) - qJ(5) + t1106;
t868 = t645 + t917;
t867 = -pkin(1) - t615;
t864 = -qJ(4) - t1090;
t863 = t638 - t1028;
t858 = qJD(2) * t977;
t857 = qJD(2) * t951;
t458 = t518 * t938;
t828 = -t458 - t548 + t932;
t736 = -t523 * t938 + t828;
t741 = t613 + t614;
t833 = -t1006 * t645 + t598 * t1009;
t225 = t833 - t915 - t450;
t847 = -t225 + t899;
t36 = t168 * t591 - t320 * t466 + (-qJD(2) * t307 + t741) * t647 + (t533 - t847) * qJD(1) + t736;
t856 = t36 * t917;
t855 = t562 - t933;
t854 = -qJ(3) * qJD(2) - qJD(5);
t852 = -t223 + t899;
t850 = -t383 + t900;
t849 = t647 * t482 + t649 * t483 + t979;
t848 = -t443 + t900;
t846 = t560 + t896;
t844 = t951 * t1091;
t822 = rSges(5,1) * t648 + t1078;
t841 = t822 + t866;
t110 = qJ(5) * t885 - t550 + (-t565 + (t616 + t1147) * t647) * qJD(1) + t897;
t705 = qJDD(1) * t455 + t971 + t1174 * t647 + (t227 + t562) * qJD(1);
t677 = qJD(1) * t286 + qJDD(1) * t483 + qJDD(4) * t649 + t705;
t672 = qJDD(1) * t450 + t647 * t919 + t677 + (t226 + t560) * qJD(1);
t826 = -t399 - t400 - t612;
t676 = -t914 + (-t274 + t826) * qJD(2);
t92 = -rSges(7,3) * t1173 + t1002;
t4 = t746 * t491 + t672 + qJD(1) * t110 + (t676 - t921) * t647 + t591 * t92 + qJDD(1) * t225 - t256 * t320 + t465 * t168 - t466 * t175;
t111 = t541 + (t645 * t646 - t1146) * t938 + (t649 * t698 - t915) * qJD(1);
t670 = -t647 * t741 + t1151 - t228 - t932;
t707 = -qJDD(4) * t647 + t492 * t1091 + t905;
t691 = t492 * t523 + t649 * t919 + t707;
t795 = t224 + t848;
t817 = rSges(7,1) * t154 + rSges(7,2) * t153;
t93 = rSges(7,3) * t419 + t817;
t5 = t465 * t167 - t467 * t175 + t257 * t320 + t492 * t307 - t591 * t93 + t676 * t649 + t795 * qJDD(1) + (-t111 + t670) * qJD(1) + t691;
t834 = t4 * t649 - t5 * t647;
t342 = -t648 * t818 + t1073;
t65 = (-qJD(2) * t342 + t741) * t647 + (t533 - t852) * qJD(1) + t736;
t831 = t65 * t869;
t829 = t448 * t938 + t455 * t936 - t935;
t827 = t560 + t855;
t825 = t869 * t648;
t824 = t868 * t648;
t528 = rSges(2,1) * t649 - rSges(2,2) * t647;
t522 = rSges(2,1) * t647 + rSges(2,2) * t649;
t527 = -rSges(3,2) * t646 + t1082;
t820 = rSges(6,1) * t279 + rSges(6,2) * t278;
t804 = t60 * t649 - t61 * t647;
t803 = t60 * t647 + t61 * t649;
t802 = t62 * t649 - t63 * t647;
t801 = t62 * t647 + t63 * t649;
t800 = t647 * t67 - t649 * t66;
t799 = t647 * t66 + t649 * t67;
t796 = t222 + t848;
t794 = t595 + t845;
t772 = -t167 * t649 - t168 * t647;
t771 = -t181 * t647 - t1026;
t252 = -rSges(3,1) * t1173 - rSges(3,2) * t883 + t956;
t446 = t521 * t647;
t255 = -qJD(2) * t446 + (t527 * t649 + t627) * qJD(1);
t767 = t252 * t649 + t255 * t647;
t754 = t385 * t647 + t388 * t649;
t744 = -t342 + t842;
t737 = t867 - t626;
t734 = -t911 - t1073;
t733 = -t910 - t1072;
t732 = t867 - t530;
t728 = t647 * t443 + t649 * t450 + t849;
t725 = t936 * t961 + t562;
t724 = -t320 + t746;
t723 = -t879 - t932;
t721 = t65 * t744;
t719 = t482 * t938 + t483 * t936 + t829;
t718 = -qJD(1) * t482 + t1156 + t855;
t708 = -qJDD(3) * t648 + t227 * t936 + t229 * t938 + t491 * t448 + t646 * t922;
t706 = t841 * t936;
t703 = -t914 + (-t400 - t477) * qJD(2);
t702 = t157 * t467 - t159 * t466 + t591 * t682;
t701 = -(-Icges(7,5) * t343 - Icges(7,6) * t344) * t467 + (Icges(7,5) * t345 - Icges(7,6) * t346) * t466 + t365 * t591;
t681 = t286 * t936 + t287 * t938 + t491 * t482 + t708;
t669 = qJD(2) * t613 + qJDD(5) * t646 + t226 * t936 + t228 * t938 + t491 * t443 + t681;
t3 = -t491 * t224 + t669 + (t110 * t649 + t111 * t647) * qJD(2) + t847 * t492 - t256 * t167 - t257 * t168 + t466 * t93 + t467 * t92;
t35 = qJD(1) * t795 - t1203 + t827;
t700 = t35 * t936 + t36 * t938 - t3;
t120 = -rSges(6,3) * t1173 + t995;
t121 = rSges(6,3) * t419 + t820;
t12 = -t222 * t491 + (t120 * t649 + t121 * t647) * qJD(2) + t852 * t492 + t669;
t64 = qJD(1) * t796 + t744 * t936 + t827;
t697 = t64 * t936 + t65 * t938 - t12;
t693 = t1060 * t646 + t1092 * t648 - pkin(1);
t692 = -qJD(1) * t443 + t560 + t718;
t689 = t648 * t701;
t688 = t443 * t938 + t450 * t936 + t612 + t719;
t686 = t648 * t789 - t1044;
t675 = -t914 + (-t318 + t826) * qJD(2);
t29 = -t167 * t466 + t168 * t467 + (-t224 * t647 + t225 * t649) * qJD(2) + t688;
t671 = -t29 * t772 + (-t35 * t647 + t36 * t649) * t320;
t656 = (t682 * t649 - t773) * t466 - (t682 * t647 - t774) * t467 + (t312 - t766) * t591;
t655 = t656 * t648;
t15 = (t675 - t921) * t647 + t744 * t491 + t672 + qJD(1) * t120 + qJDD(1) * t223;
t16 = t492 * t342 + t675 * t649 + t796 * qJDD(1) + (-t121 + t670) * qJD(1) + t691;
t53 = (-t222 * t647 + t223 * t649) * qJD(2) + t688;
t652 = (qJD(2) * t53 + t15 * t647 + t16 * t649 - t64 * t942 + t65 * t940) * t1109 + (qJD(2) * t29 - t35 * t942 + t36 * t940 + t4 * t647 + t5 * t649) * t1108;
t587 = rSges(4,3) * t1006;
t579 = rSges(4,3) * t1007;
t498 = t649 * t913;
t497 = t647 * t913;
t490 = t649 * t912;
t489 = t647 * t912;
t479 = t527 * qJD(2);
t452 = -rSges(4,1) * t1009 + t587;
t445 = -rSges(4,1) * t1010 + t579;
t422 = t883 - t891;
t418 = t648 * t857;
t417 = t646 * t857;
t319 = -t1075 * t646 + t959;
t289 = t649 * t734 + t498;
t288 = t647 * t734 + t497;
t285 = t687 * t649;
t284 = t687 * t647;
t283 = t685 * t649;
t282 = t685 * t647;
t269 = t649 * t733 + t490;
t268 = t647 * t733 + t489;
t267 = t686 * t649;
t266 = t686 * t647;
t265 = t684 * t649;
t264 = t684 * t647;
t259 = -t454 + t968;
t258 = -t447 + t969;
t254 = -t520 * t938 + (t649 * t823 + t628) * qJD(1);
t253 = t822 * t938 + (t1145 * t649 - t1074) * qJD(1);
t251 = -rSges(4,1) * t1173 - rSges(4,3) * t891 + t957;
t250 = qJD(1) * t735 + t895;
t209 = rSges(7,1) * t345 - rSges(7,2) * t346;
t208 = -rSges(7,1) * t343 - rSges(7,2) * t344;
t176 = t754 * qJD(2);
t113 = -t458 + (-qJD(2) * t520 + t614) * t647 + t231 * qJD(1);
t112 = qJD(1) * t903 + t725;
t108 = (t384 * t647 + t387 * t649) * qJD(2) + t829;
t101 = (qJD(2) * t822 + t614) * t647 + (t533 - t902) * qJD(1) + t828;
t100 = qJD(1) * t850 + t706 + t855;
t97 = (t383 * t647 + t386 * t649) * qJD(2) + t719;
t95 = qJD(1) * t252 + qJDD(1) * t388 - t479 * t938 - t491 * t521 + t971;
t94 = -t479 * t936 + t492 * t521 + t981 * qJDD(1) + (-t255 - t488) * qJD(1);
t52 = qJD(1) * t251 + qJDD(1) * t387 + t491 * t961 + t647 * t858 + t705;
t51 = t492 * t520 + t649 * t858 + t903 * qJDD(1) + (-t254 - t879 + t996) * qJD(1) + t905;
t34 = t384 * t491 + t980 * t492 + (t251 * t649 + t254 * t647) * qJD(2) + t708;
t32 = (t703 - t921) * t647 + qJD(1) * t250 + qJDD(1) * t386 + t677 + t841 * t491;
t31 = -t492 * t822 + t703 * t649 + t850 * qJDD(1) + (-t253 + t723 + t1151) * qJD(1) + t707;
t26 = t1007 * t169 - t153 * t684 + t154 * t317 - t170 * t343 + t171 * t344 - t419 * t682;
t25 = t1006 * t169 + t1173 * t682 - t151 * t684 + t152 * t317 + t170 * t345 + t171 * t346;
t24 = t383 * t491 + (t250 * t649 + t253 * t647) * qJD(2) + t902 * t492 + t681;
t21 = t107 * t591 + t466 * t67 - t467 * t66;
t9 = t1007 * t86 + t153 * t162 + t154 * t165 + t159 * t419 - t343 * t88 + t344 * t90;
t8 = t1007 * t87 + t153 * t160 - t154 * t164 + t157 * t419 - t343 * t89 + t344 * t91;
t7 = t1006 * t86 - t1173 * t159 + t151 * t162 + t152 * t165 + t345 * t88 + t346 * t90;
t6 = t1006 * t87 - t1173 * t157 + t151 * t160 - t152 * t164 + t345 * t89 + t346 * t91;
t2 = t256 * t61 + t257 * t60 + t26 * t591 + t465 * t98 + t466 * t9 - t467 * t8;
t1 = t25 * t591 + t256 * t63 + t257 * t62 + t465 * t99 + t466 * t7 - t467 * t6;
t17 = [(t26 + t14) * t1100 + t14 * t1099 + (t5 * (t638 + t1194) + t35 * (-t817 + t894) + (-t5 * qJ(4) + t856 * t939 + (-qJD(4) + (t732 + t824) * qJD(1)) * t35) * t649 + (-t5 * pkin(1) + (-t5 * qJ(3) + t35 * (-qJD(3) + (rSges(7,3) - t645) * qJD(2))) * t646 + t35 * (-pkin(7) + t1090) * qJD(1) + (t5 * t868 + t35 * (-qJD(2) * t598 + t854) + qJD(1) * t856) * t648) * t647 - g(1) * (t863 + t1194) - (t824 + t867) * t1089 + (-t933 + t1002 + t35 - t692 + t846 + t897 + (t647 * t732 + t649 * t864 - t224) * qJD(1) + t1203) * t36 + (t4 - g(2)) * (t647 * t864 + t168 + t794 + t833)) * m(7) + (t16 * (t638 + t1157) + t64 * (t541 - t820 + t894) + t65 * (t550 + t846 + t995) + (-t16 * qJ(4) - t64 * qJD(4) + t831 * t939 + (-t65 * qJ(4) + (t867 - t634 + t825) * t64) * qJD(1)) * t649 + (-t16 * pkin(1) - t64 * pkin(7) * qJD(1) + t65 * (-pkin(1) * qJD(1) - qJD(4)) + (t64 * (rSges(6,3) * qJD(2) - qJD(3)) + (qJD(1) * t65 + t16) * t1084) * t646 + (t16 * t869 + t64 * (t854 - t1071) + qJD(1) * t831) * t648) * t647 - (-t64 + t692 + t1201) * t65 - t721 * t936 - g(1) * (t863 + t1157) - (t1084 * t646 - pkin(1) + t825) * t1089 + (t15 - g(2)) * (t223 + t450 + t794 - t1029)) * m(6) + ((t309 * t643 - t310 * t644 - t1219) * t648 + (t308 + t1218) * t646 + (t336 * t648 - t646 * t765 + t1217) * qJD(2)) * qJD(1) + (t1161 - t1121) * t1097 + (t1160 + t1120) * t1098 + t1069 / 0.2e1 - t1070 / 0.2e1 + t1063 / 0.2e1 + t1064 / 0.2e1 + t1058 + (-(-t100 + t706 + t718 - t1149) * t101 + t100 * (t723 + t894) + t101 * (t542 + t606 + t855 + t895) + (t101 * t646 * t907 + t100 * (t1061 * t648 - t1078) * t647) * qJD(2) + ((t100 * (rSges(5,3) - pkin(7)) + t101 * (-t636 + t737 - t637)) * t647 + (t101 * t1059 + (t737 + (rSges(5,2) + t1106) * t648) * t100) * t649) * qJD(1) + (-g(2) + t32) * (t1059 * t647 + t794 + t838) + (-g(1) + t31) * (t580 + t638 + t1059 * t649 + (t1061 * t646 + t1106 * t648 - pkin(1)) * t647)) * m(5) + (t181 * (t606 + t956) + (t1027 * t521 - t1025) * qJD(2) + ((-pkin(1) - t527) * t1026 + (t180 * (-rSges(3,3) - pkin(7)) + t181 * t878) * t647) * qJD(1) - (-qJD(1) * t385 - t180 - t494 - t887) * t181 + (t95 - g(2)) * t291 + (t94 - g(1)) * (t647 * t878 + t638 + t952)) * m(3) + t25 * t1101 + t98 * t1104 + t99 * t1105 - m(2) * (-g(1) * t522 + g(2) * t528) + (t336 * t646 + t648 * t765 + m(2) * (t522 ^ 2 + t528 ^ 2) + Icges(2,3) - t1216) * qJDD(1) + (((t1175 * t649 + t1168 + t1169 + t770) * t649 + (t1175 * t647 + t1123 + t1159 - t1181 - t1233 - t323 + t325) * t647) * qJD(2) + t1167 - t1182) * t874 + ((-t1065 + ((t352 + t1020) * t649 + t1177 + t1223 + t1239) * t649 + (-t1158 + t71 + (-t348 + t760) * t647 - t1168 + t1195) * t647) * qJD(2) + t1185) * t871 + (t1163 + t1164) * t873 + (t1162 - t1165 + t1166) * t872 + (-(-qJD(1) * t384 - t112 + t1156 + t725) * t113 - t112 * t830 + t113 * (t896 + t957) + (t113 * t646 * t892 + t112 * (t1060 * t648 + t1081) * t647) * qJD(2) + (t112 * t693 * t649 + (t112 * (-rSges(4,2) - pkin(7)) + t113 * (-pkin(1) + t960)) * t647) * qJD(1) + (t52 - g(2)) * t231 + (t51 - g(1)) * (t647 * t693 + t633 + t638)) * m(4); (t3 * t728 + (t5 * t724 + t3 * (t168 + t225)) * t649 + (t3 * t997 + t4 * t724) * t647 - g(1) * (-t649 * t910 + t490 + t568 + t968) - g(2) * (-t647 * t910 + t489 + t564 + t969) - g(3) * (t859 + t893 + t959) + (g(3) * t1075 - t671 * qJD(6) - t836 * t917) * t646 + (-t168 * t928 - t269 * t591 + t319 * t466 - t901 - (-pkin(3) * t940 - t931) * t646 + t1128 * t647 + (t649 * t724 - t259) * qJD(1)) * t36 + (-t167 * t928 + t268 * t591 + t319 * t467 + t1128 * t649 + ((t307 + t320) * t647 + t258 - t975) * qJD(1) + t1115) * t35 + ((qJD(1) * t997 + t110 + t92) * t649 + (t111 + t93 + (-t168 + t847) * qJD(1)) * t647 - t268 * t466 - t269 * t467 - (t258 * t647 + t259 * t649 - t844) * qJD(2) + t1153) * t29) * m(7) + (g(1) * t453 + g(2) * t446 - g(3) * t527 - (t180 * t446 - t1025) * qJD(1) - (t176 * (-t446 * t647 - t453 * t649) + t771 * t527) * qJD(2) + (qJD(2) * t767 + t385 * t491 - t388 * t492) * t754 + t176 * ((t385 * t649 - t388 * t647) * qJD(1) + t767) + t771 * t479 + (-t95 * t647 - t94 * t649 + (-t181 * t649 + t1027) * qJD(1)) * t521) * m(3) - t21 * t928 / 0.2e1 + (((t265 * t607 - t267 * t608 + t159) * t466 - (t264 * t607 - t266 * t608 + t157) * t467 + (t314 * t607 - t316 * t608 - t682) * t591 + t107 * qJD(6)) * t648 + (-qJD(6) * t799 + t656) * t646) * t1096 + (t13 * t647 + t14 * t649) * t929 / 0.2e1 + (t24 * t849 + (t24 * t386 + t31 * t841) * t649 + (t24 * t383 + t32 * t841) * t647 - g(1) * (t568 + t456) - g(2) * (t564 + t449) - g(3) * t1193 - (g(1) * t907 + t1087 * t1106) * t646 + (-t976 + t1126 * t647 + (t649 * t841 - t456 + t916) * qJD(1)) * t101 + (t1126 * t649 + (-t647 * t822 + t444 + t449) * qJD(1) + t1172) * t100 + ((t250 + t1149) * t649 + (qJD(1) * t902 + t253) * t647 - (t449 * t647 + t456 * t649 - t844) * qJD(2) + t1176) * t97) * m(5) + t1197 * t1098 + ((t283 * t414 + t285 * t415) * t938 + (t337 * t414 + t339 * t415) * qJD(1) + (-t1141 * t938 - t1140) * t647 + ((t1142 * t647 - t414 * t282 - t415 * t284 + t1118) * qJD(2) + t1119) * t649) * t874 + (-(-t282 * t412 + t284 * t413) * t936 + (-t337 * t412 + t339 * t413) * qJD(1) + (-t1142 * t936 + t1140) * t649 + ((t1141 * t649 - t412 * t283 + t413 * t285 + t1118) * qJD(2) + t1119) * t647) * t871 - t256 * t802 / 0.2e1 - t257 * t804 / 0.2e1 + (t12 * t728 + (qJD(1) * t721 + t12 * t223 + t16 * t744) * t649 + (-t12 * t222 + t15 * t744) * t647 - g(1) * (-t649 * t911 + t498 + t568 + t594) - g(2) * (-t647 * t911 + t497 + t564 + t592) - g(3) * t1202 - (-g(3) * t1076 + t836 * t869) * t646 + (t646 * t931 - t901 - (t289 - t916) * qJD(1) + t1127 * t647) * t65 + (t1127 * t649 + (t342 * t647 + t288 - t975) * qJD(1) + t1115) * t64 + ((t120 - t1201) * t649 + (qJD(1) * t852 + t121) * t647 - (t288 * t647 + t289 * t649 - t844) * qJD(2) + t1153) * t53) * m(6) + ((t1124 * t649 + t1125 * t647) * qJD(1) + t1117) * t872 + ((t1122 * t649 + t1123 * t647) * qJD(1) + t1116) * t873 + (t1120 * t647 + t1121 * t649) * qJDD(1) / 0.2e1 + ((-t265 * t343 + t267 * t344) * t466 - (-t264 * t343 + t266 * t344) * t467 + (-t314 * t343 + t316 * t344) * t591 + (-t1009 * t61 + t648 * t98) * qJD(6) + ((-qJD(6) * t60 + t702) * t646 + t655) * t647) * t1099 + (qJD(1) * t803 + t647 * t9 - t649 * t8) * t1100 + (qJD(1) * t801 - t6 * t649 + t647 * t7) * t1101 + ((t265 * t345 + t267 * t346) * t466 - (t264 * t345 + t266 * t346) * t467 + (t314 * t345 + t316 * t346) * t591 + (-t1010 * t62 + t648 * t99) * qJD(6) + ((-qJD(6) * t63 + t702) * t646 + t655) * t649) * t1102 + t800 * t1103 + (qJD(1) * t799 - t10 * t649 + t11 * t647) * t1095 + (-t112 * (t563 + (-t444 - t445) * qJD(1)) - t113 * (qJD(1) * t452 + t976) - t108 * t904 - ((t108 * t452 + t112 * t960) * t649 + (t108 * t445 + t113 * t960) * t647) * qJD(2) + t112 * t460 + t34 * t979 + t108 * t906 + (t51 * t961 + t112 * t977 + t34 * t387 + t108 * t251 + (t108 * t384 + t113 * t961) * qJD(1)) * t649 + (t52 * t961 + t113 * t977 + t34 * t384 + t108 * t254 + (t108 * t980 + t112 * t520) * qJD(1)) * t647 - g(1) * (t568 + t587) - g(2) * (t564 + t579) + g(3) * t960 - (g(1) * t892 + t1087 * t1092) * t646) * m(4) + t1196 * t1097 - (t654 * t646 + (t1190 * t646 + ((-t282 * t643 + t284 * t644 - t1184) * t649 + (t283 * t643 - t285 * t644 - t1183) * t647) * t648) * qJD(2) + (t1188 * t646 + (t337 * t643 - t339 * t644 - t1189 + t336) * t648) * qJD(1)) * qJD(1) / 0.2e1 - (qJD(1) * t1162 + t1117 * qJD(2) + qJDD(1) * t1161 + t1124 * t491 + t1125 * t492 + t2) * t649 / 0.2e1 + (qJD(1) * t1163 + t1116 * qJD(2) + qJDD(1) * t1160 + t1122 * t491 + t1123 * t492 + t1) * t647 / 0.2e1 + (t1165 * t649 + t1164 * t647 + (t1120 * t649 - t1121 * t647) * qJD(1)) * qJD(1) / 0.2e1 + (t14 + t1166) * t875 + (t13 + t1167) * t942 / 0.2e1; (-m(4) + t924) * (-g(3) * t648 + t646 * t836) - m(4) * (t108 * t417 + t112 * t422 + t113 * t420) - m(5) * (t100 * t422 + t101 * t420 + t417 * t97) - m(6) * (t417 * t53 + t420 * t65 + t422 * t64) - m(7) * (t29 * t417 + t35 * t422 + t36 * t420) + ((t112 * t936 + t113 * t938 - t34) * t1111 + (t100 * t936 + t101 * t938 - t24) * t1110 + t697 * t1109 + t700 * t1108) * t1112 + ((qJD(2) * t108 - t112 * t942 + t113 * t940 + t51 * t649 + t52 * t647) * t1111 + (qJD(2) * t97 - t100 * t942 + t101 * t940 + t31 * t649 + t32 * t647) * t1110 + t652) * t1113; t924 * (g(2) * t649 - t1089) + m(5) * (-t31 * t647 + t32 * t649) + m(6) * (t15 * t649 - t16 * t647) + m(7) * t834; t1107 * (g(3) * t646 + t648 * t836) - m(6) * (-t1173 * t64 + t418 * t53 + t419 * t65) - m(7) * (-t1173 * t35 + t29 * t418 + t36 * t419) + (-m(6) * t697 / 0.2e1 - m(7) * t700 / 0.2e1) * t1113 + t652 * t1112; t1 * t1006 / 0.2e1 + (t646 * t99 + t648 * t801) * t1105 + ((-qJD(2) * t801 + t25) * t646 + (qJD(1) * t802 + qJD(2) * t99 + t6 * t647 + t649 * t7) * t648) * t1101 + t2 * t1007 / 0.2e1 + (t646 * t98 + t648 * t803) * t1104 + ((-qJD(2) * t803 + t26) * t646 + (qJD(1) * t804 + qJD(2) * t98 + t647 * t8 + t649 * t9) * t648) * t1100 + t21 * t937 / 0.2e1 + t646 * (t1058 + t1063 + t1064 + t1069 - t1070) / 0.2e1 + (t107 * t646 + t648 * t799) * t1103 + ((-qJD(2) * t799 + t33) * t646 + (-qJD(1) * t800 + qJD(2) * t107 + t10 * t647 + t11 * t649) * t648) * t1095 + (-t1114 * t346 + t345 * t673 + t649 * t689) * t1102 + (-t1114 * t344 - t343 * t673 + t647 * t689) * t1099 + (t701 * t646 + (t1114 * t608 + t673 * t607) * t648) * t1096 + (-t889 / 0.2e1 + t646 * t872) * t14 + (t646 * t874 + t648 * t875) * t13 + ((qJD(2) * t671 + t167 * t5 + t4 * t168 - t35 * t93 + t36 * t92) * t646 + (t35 * (qJD(2) * t167 + t175 * t647) + t36 * (qJD(2) * t168 - t175 * t649) + t3 * t772 + t29 * (t167 * t942 - t168 * t940 - t647 * t92 + t649 * t93) + ((t35 * t649 + t36 * t647) * qJD(1) - t834) * t320) * t648 - t35 * (-t208 * t591 - t382 * t467) - t36 * (t209 * t591 - t382 * t466) - t29 * (t208 * t466 + t209 * t467) - g(1) * t209 - g(2) * t208 - g(3) * t382) * m(7);];
tau  = t17;
