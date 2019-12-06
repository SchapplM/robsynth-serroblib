% Calculate vector of inverse dynamics joint torques for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:03
% EndTime: 2019-12-05 19:00:14
% DurationCPUTime: 5.00s
% Computational Cost: add. (8737->482), mult. (13204->639), div. (0->0), fcn. (8699->16), ass. (0->253)
t277 = sin(qJ(4));
t278 = sin(qJ(3));
t282 = cos(qJ(4));
t283 = cos(qJ(3));
t206 = -t277 * t278 + t282 * t283;
t286 = -pkin(8) - pkin(7);
t315 = qJD(3) * t286;
t216 = t278 * t315;
t217 = t283 * t315;
t230 = t286 * t278;
t266 = t283 * pkin(8);
t231 = pkin(7) * t283 + t266;
t284 = cos(qJ(2));
t344 = qJD(1) * pkin(1);
t318 = t284 * t344;
t322 = qJD(4) * t282;
t323 = qJD(4) * t277;
t378 = -t206 * t318 + t282 * t216 + t277 * t217 + t230 * t322 - t231 * t323;
t160 = t277 * t230 + t282 * t231;
t207 = t277 * t283 + t278 * t282;
t377 = -qJD(4) * t160 + t207 * t318 - t216 * t277 + t282 * t217;
t297 = t206 * qJD(4);
t141 = t206 * qJD(3) + t297;
t363 = pkin(9) * t141;
t395 = -t363 + t377;
t298 = t207 * qJD(4);
t142 = -t207 * qJD(3) - t298;
t139 = t142 * pkin(9);
t394 = -t139 - t378;
t369 = t278 / 0.2e1;
t270 = qJD(1) + qJD(2);
t180 = t206 * t270;
t181 = t207 * t270;
t276 = sin(qJ(5));
t281 = cos(qJ(5));
t309 = t281 * t180 - t181 * t276;
t113 = Ifges(6,4) * t309;
t117 = t180 * t276 + t181 * t281;
t269 = qJD(3) + qJD(4);
t260 = qJD(5) + t269;
t54 = Ifges(6,1) * t117 + Ifges(6,5) * t260 + t113;
t393 = t113 + t54;
t251 = pkin(3) * t282 + pkin(4);
t320 = qJD(5) * t281;
t321 = qJD(5) * t276;
t327 = t277 * t281;
t279 = sin(qJ(2));
t319 = t279 * t344;
t221 = pkin(7) * t270 + t319;
t313 = pkin(8) * t270 + t221;
t166 = t313 * t283;
t155 = t282 * t166;
t165 = t313 * t278;
t104 = t165 * t277 - t155;
t362 = pkin(9) * t180;
t79 = t104 - t362;
t153 = t277 * t166;
t105 = -t282 * t165 - t153;
t174 = t181 * pkin(9);
t80 = -t174 + t105;
t392 = t276 * t80 - t281 * t79 - t251 * t321 + (-t277 * t320 + (-t276 * t282 - t327) * qJD(4)) * pkin(3);
t328 = t276 * t277;
t391 = -t276 * t79 - t281 * t80 + t251 * t320 + (-t277 * t321 + (t281 * t282 - t328) * qJD(4)) * pkin(3);
t222 = -pkin(2) * t270 - t318;
t355 = mrSges(4,2) * t283;
t307 = mrSges(4,1) * t278 + t355;
t390 = t222 * t307 + qJD(3) * (Ifges(4,5) * t283 - Ifges(4,6) * t278) / 0.2e1;
t389 = -m(4) * pkin(7) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3) + mrSges(3,2) + m(6) * (-pkin(9) + t286) + m(5) * t286;
t347 = Ifges(6,4) * t117;
t53 = Ifges(6,2) * t309 + Ifges(6,6) * t260 + t347;
t387 = t53 / 0.2e1;
t375 = -t309 / 0.2e1;
t384 = mrSges(6,2) * t309;
t383 = Ifges(6,1) * t309;
t382 = Ifges(6,5) * t309;
t343 = qJD(3) * pkin(3);
t156 = -t165 + t343;
t101 = t156 * t277 + t155;
t77 = t101 + t362;
t340 = t276 * t77;
t100 = t282 * t156 - t153;
t76 = t100 - t174;
t72 = pkin(4) * t269 + t76;
t35 = t281 * t72 - t340;
t381 = t309 * t35;
t159 = t282 * t230 - t231 * t277;
t361 = pkin(9) * t207;
t125 = t159 - t361;
t199 = t206 * pkin(9);
t126 = t199 + t160;
t67 = t125 * t281 - t126 * t276;
t380 = t67 * qJD(5) + t395 * t276 - t394 * t281;
t68 = t125 * t276 + t126 * t281;
t379 = -t68 * qJD(5) + t394 * t276 + t395 * t281;
t250 = pkin(1) * t279 + pkin(7);
t358 = -pkin(8) - t250;
t200 = t358 * t278;
t333 = t250 * t283;
t201 = t266 + t333;
t135 = t277 * t200 + t282 * t201;
t330 = t270 * t278;
t219 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t330;
t329 = t270 * t283;
t220 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t329;
t376 = (mrSges(3,2) * t270 + t219 * t278 - t220 * t283) * t284 - m(4) * (t222 * t279 + (t278 ^ 2 + t283 ^ 2) * t284 * t221);
t374 = -t117 / 0.2e1;
t373 = t117 / 0.2e1;
t371 = t181 / 0.2e1;
t370 = -t260 / 0.2e1;
t274 = qJ(3) + qJ(4);
t261 = sin(t274);
t364 = pkin(3) * t278;
t367 = m(6) * (-pkin(4) * t261 - t364);
t339 = t281 * t77;
t36 = t276 * t72 + t339;
t366 = mrSges(6,3) * t36;
t365 = pkin(1) * t284;
t360 = g(1) * t283;
t275 = qJ(1) + qJ(2);
t264 = cos(t275);
t359 = g(3) * t264;
t357 = mrSges(4,1) * t283;
t356 = mrSges(4,2) * t278;
t354 = mrSges(5,2) * t261;
t263 = cos(t274);
t353 = mrSges(5,2) * t263;
t265 = qJ(5) + t274;
t249 = cos(t265);
t352 = mrSges(6,2) * t249;
t351 = mrSges(5,3) * t180;
t350 = Ifges(4,4) * t278;
t349 = Ifges(4,4) * t283;
t348 = Ifges(5,4) * t181;
t346 = Ifges(4,2) * t283;
t345 = pkin(1) * qJD(2);
t342 = t181 * mrSges(5,3);
t248 = sin(t265);
t341 = t248 * mrSges(6,2);
t237 = t249 * mrSges(6,1);
t242 = t263 * mrSges(5,1);
t338 = qJDD(1) * pkin(1);
t316 = qJD(2) * t344;
t204 = t279 * t338 + t284 * t316;
t268 = qJDD(1) + qJDD(2);
t184 = pkin(7) * t268 + t204;
t325 = qJD(3) * t278;
t127 = t283 * t184 - t221 * t325;
t337 = t127 * t283;
t324 = qJD(3) * t283;
t128 = -t184 * t278 - t221 * t324;
t336 = t128 * t278;
t262 = sin(t275);
t334 = t248 * t262;
t332 = t261 * t262;
t331 = t262 * t278;
t326 = mrSges(6,1) * t334 + t262 * t352;
t267 = qJDD(3) + qJDD(4);
t317 = t284 * t345;
t255 = pkin(3) * t325;
t252 = pkin(3) * t283 + pkin(2);
t129 = -pkin(4) * t142 + t255;
t312 = qJD(3) * t358;
t229 = t356 - t357;
t311 = (-mrSges(3,1) + t229) * t270;
t310 = t237 - t341;
t134 = t282 * t200 - t201 * t277;
t308 = -mrSges(5,1) * t332 - t262 * t353 - t326;
t203 = -t279 * t316 + t284 * t338;
t306 = -mrSges(6,1) * t248 - t352;
t305 = t346 + t350;
t108 = t134 - t361;
t109 = t199 + t135;
t50 = t108 * t281 - t109 * t276;
t51 = t108 * t276 + t109 * t281;
t137 = t206 * t281 - t207 * t276;
t138 = t206 * t276 + t207 * t281;
t303 = t219 * t283 + t220 * t278;
t177 = -pkin(4) * t206 - t252;
t302 = mrSges(2,1) + (m(3) + m(4) + m(5) + m(6)) * pkin(1);
t259 = qJDD(5) + t267;
t193 = t268 * t278 + t270 * t324;
t96 = qJDD(3) * pkin(3) - pkin(8) * t193 + t128;
t192 = t268 * t283 - t270 * t325;
t99 = pkin(8) * t192 + t127;
t29 = -t101 * qJD(4) - t277 * t99 + t282 * t96;
t93 = t192 * t277 + t193 * t282 + t270 * t297;
t15 = pkin(4) * t267 - pkin(9) * t93 + t29;
t28 = t156 * t322 - t166 * t323 + t277 * t96 + t282 * t99;
t94 = t192 * t282 - t193 * t277 - t270 * t298;
t19 = pkin(9) * t94 + t28;
t3 = t35 * qJD(5) + t15 * t276 + t19 * t281;
t32 = qJD(5) * t309 + t276 * t94 + t281 * t93;
t33 = -t117 * qJD(5) - t276 * t93 + t281 * t94;
t4 = -t36 * qJD(5) + t15 * t281 - t19 * t276;
t301 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t32 + Ifges(6,6) * t33 + Ifges(6,3) * t259;
t299 = t278 * (Ifges(4,1) * t283 - t350);
t161 = t278 * t312 + t283 * t317;
t162 = -t278 * t317 + t283 * t312;
t60 = t282 * t161 + t277 * t162 + t200 * t322 - t201 * t323;
t183 = -pkin(2) * t268 - t203;
t296 = -t242 - t310 + t354;
t182 = -t252 * t270 - t318;
t130 = -pkin(3) * t192 + t183;
t293 = mrSges(5,1) * t261 - t306 + t353;
t61 = -qJD(4) * t135 - t161 * t277 + t282 * t162;
t292 = -t303 * qJD(3) + t283 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t192) - t278 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t193);
t247 = pkin(4) * t263;
t291 = m(4) * pkin(2) + m(5) * t252 + m(6) * (t247 + t252) + mrSges(3,1) + t237 + t242 + t357;
t111 = Ifges(5,2) * t180 + Ifges(5,6) * t269 + t348;
t173 = Ifges(5,4) * t180;
t112 = Ifges(5,1) * t181 + Ifges(5,5) * t269 + t173;
t124 = -pkin(4) * t180 + t182;
t290 = mrSges(6,3) * t381 + t29 * mrSges(5,1) - t28 * mrSges(5,2) - t182 * (mrSges(5,1) * t181 + mrSges(5,2) * t180) - t124 * t384 - t269 * (Ifges(5,5) * t180 - Ifges(5,6) * t181) / 0.2e1 + Ifges(5,3) * t267 + t301 + t383 * t374 + t382 * t370 + t100 * t351 + t111 * t371 - t181 * (Ifges(5,1) * t180 - t348) / 0.2e1 + Ifges(5,6) * t94 + Ifges(5,5) * t93 + t393 * t375 - (-Ifges(5,2) * t181 + t112 + t173) * t180 / 0.2e1 + (-t124 * mrSges(6,1) - Ifges(6,4) * t374 - Ifges(6,2) * t375 - Ifges(6,6) * t370 + t366 + t387) * t117;
t289 = (-t341 - t354 - t356 + t291) * t264 - t389 * t262;
t288 = -mrSges(4,2) * t331 - mrSges(5,2) * t332 - mrSges(6,2) * t334 + t291 * t262 + t389 * t264;
t178 = Ifges(4,6) * qJD(3) + t305 * t270;
t235 = Ifges(4,4) * t329;
t179 = Ifges(4,1) * t330 + Ifges(4,5) * qJD(3) + t235;
t55 = t137 * qJD(5) + t141 * t281 + t142 * t276;
t56 = -t138 * qJD(5) - t141 * t276 + t142 * t281;
t57 = -pkin(4) * t94 + t130;
t287 = (t270 * (-Ifges(4,2) * t278 + t349) + t179) * t324 / 0.2e1 + (Ifges(4,1) * t193 + Ifges(4,4) * t192) * t369 + t309 * (Ifges(6,4) * t55 + Ifges(6,2) * t56) / 0.2e1 + (t299 * t270 / 0.2e1 + t390) * qJD(3) + (-mrSges(6,1) * t57 + mrSges(6,3) * t3 + Ifges(6,4) * t32 + Ifges(6,2) * t33 + Ifges(6,6) * t259) * t137 - t35 * t55 * mrSges(6,3) + t56 * t387 + (-t100 * t141 + t101 * t142) * mrSges(5,3) + t283 * (Ifges(4,4) * t193 + Ifges(4,2) * t192) / 0.2e1 + (mrSges(6,2) * t57 - mrSges(6,3) * t4 + Ifges(6,1) * t32 + Ifges(6,4) * t33 + Ifges(6,5) * t259) * t138 + (t130 * mrSges(5,2) - t29 * mrSges(5,3) + Ifges(5,1) * t93 + Ifges(5,4) * t94 + Ifges(5,5) * t267) * t207 + (-t130 * mrSges(5,1) + t28 * mrSges(5,3) + Ifges(5,4) * t93 + Ifges(5,2) * t94 + Ifges(5,6) * t267) * t206 + t56 * t366 + (Ifges(5,1) * t141 + Ifges(5,4) * t142) * t371 + (Ifges(6,1) * t55 + Ifges(6,4) * t56) * t373 + mrSges(4,3) * t337 + t192 * t305 / 0.2e1 + (0.2e1 * Ifges(4,5) * t369 + Ifges(4,6) * t283) * qJDD(3) - t178 * t325 / 0.2e1 + t193 * (t278 * Ifges(4,1) + t349) / 0.2e1 + t55 * t54 / 0.2e1 + t124 * (-mrSges(6,1) * t56 + mrSges(6,2) * t55) + t141 * t112 / 0.2e1 + t142 * t111 / 0.2e1 + t180 * (Ifges(5,4) * t141 + Ifges(5,2) * t142) / 0.2e1 + t182 * (-mrSges(5,1) * t142 + mrSges(5,2) * t141) + t203 * mrSges(3,1) - t204 * mrSges(3,2) + t183 * t229 + t260 * (Ifges(6,5) * t55 + Ifges(6,6) * t56) / 0.2e1 + Ifges(3,3) * t268 + t269 * (Ifges(5,5) * t141 + Ifges(5,6) * t142) / 0.2e1;
t285 = cos(qJ(1));
t280 = sin(qJ(1));
t256 = t279 * t345;
t253 = -pkin(2) - t365;
t228 = -t252 - t365;
t218 = t256 + t255;
t190 = pkin(3) * t327 + t251 * t276;
t189 = -pkin(3) * t328 + t251 * t281;
t164 = t177 - t365;
t152 = mrSges(5,1) * t269 - t342;
t151 = -mrSges(5,2) * t269 + t351;
t143 = pkin(3) * t330 + pkin(4) * t181;
t131 = -mrSges(4,1) * t192 + mrSges(4,2) * t193;
t122 = -mrSges(5,1) * t180 + mrSges(5,2) * t181;
t120 = t129 + t256;
t103 = mrSges(6,1) * t260 - mrSges(6,3) * t117;
t102 = -mrSges(6,2) * t260 + mrSges(6,3) * t309;
t85 = -mrSges(5,2) * t267 + mrSges(5,3) * t94;
t84 = mrSges(5,1) * t267 - mrSges(5,3) * t93;
t59 = -mrSges(6,1) * t309 + mrSges(6,2) * t117;
t49 = t61 - t363;
t48 = t139 + t60;
t47 = -mrSges(5,1) * t94 + mrSges(5,2) * t93;
t38 = t281 * t76 - t340;
t37 = -t276 * t76 - t339;
t26 = -mrSges(6,2) * t259 + mrSges(6,3) * t33;
t25 = mrSges(6,1) * t259 - mrSges(6,3) * t32;
t9 = -mrSges(6,1) * t33 + mrSges(6,2) * t32;
t8 = -t51 * qJD(5) - t276 * t48 + t281 * t49;
t7 = t50 * qJD(5) + t276 * t49 + t281 * t48;
t1 = [(mrSges(2,2) * t285 + t302 * t280 + t288) * g(3) + t292 * t250 + t287 + m(4) * (t127 * t333 + t183 * t253 - t250 * t336) + m(6) * (t120 * t124 + t164 * t57 + t3 * t51 + t35 * t8 + t36 * t7 + t4 * t50) + m(5) * (t100 * t61 + t101 * t60 + t130 * t228 + t134 * t29 + t135 * t28 + t182 * t218) + (-mrSges(2,2) * t280 + t302 * t285 + t289) * g(2) - mrSges(4,3) * t336 + t50 * t25 + t51 * t26 + t7 * t102 + t8 * t103 + t120 * t59 + t134 * t84 + t135 * t85 + t60 * t151 + t61 * t152 + t164 * t9 + t218 * t122 + Ifges(2,3) * qJDD(1) + t228 * t47 + t253 * t131 + ((mrSges(3,1) * t284 - mrSges(3,2) * t279) * t268 + m(3) * (t203 * t284 + t204 * t279) + (t311 * t279 - t376) * qJD(2)) * pkin(1); t288 * g(3) + t292 * pkin(7) + t289 * g(2) + (-t128 * mrSges(4,3) + t122 * t343) * t278 + t287 + t380 * t102 + t379 * t103 + m(4) * (-pkin(2) * t183 + (-t336 + t337) * pkin(7)) + t67 * t25 + t68 * t26 + t129 * t59 - pkin(2) * t131 + t159 * t84 + t160 * t85 + t177 * t9 - t252 * t47 + t378 * t151 + t377 * t152 + ((-t122 - t311 - t59) * t279 + t376) * t344 + (t177 * t57 + t3 * t68 + t4 * t67 + t380 * t36 + t379 * t35 + (t129 - t319) * t124) * m(6) + (-t130 * t252 + t159 * t29 + t160 * t28 + (t255 - t319) * t182 + t378 * t101 + t377 * t100) * m(5); (-t367 + t355 + (m(5) * pkin(3) + mrSges(4,1)) * t278 + t293) * t359 + t101 * t342 + (t277 * t85 + t282 * t84 + (t282 * t151 - t277 * t152) * qJD(4) + (-g(2) * t331 - t100 * t323 + t101 * t322 + t277 * t28 + t282 * t29 - t360) * m(5)) * pkin(3) + (t178 * t369 + (t346 * t369 - t299 / 0.2e1) * t270 + (-m(5) * t182 - t122) * t364 - (t179 + t235) * t283 / 0.2e1 - t390) * t270 + t303 * t221 + ((-t307 + t367) * t262 + t308) * g(2) + (t229 + t296) * g(1) - m(5) * (t100 * t104 + t101 * t105) + t290 + Ifges(4,3) * qJDD(3) - t127 * mrSges(4,2) + t128 * mrSges(4,1) - t143 * t59 - t105 * t151 - t104 * t152 + t189 * t25 + t190 * t26 + Ifges(4,6) * t192 + Ifges(4,5) * t193 + t391 * t102 + t392 * t103 + (-t360 * pkin(3) - t247 * g(1) - t124 * t143 + t189 * t4 + t190 * t3 + t392 * t35 + t391 * t36) * m(6); t293 * t359 + (-t181 * t59 + t281 * t25 + t276 * t26 + (t102 * t281 - t103 * t276) * qJD(5) + (-t124 * t181 - t35 * t321 + t36 * t320 + t276 * t3 + t281 * t4 - g(1) * t263 + (-g(2) * t262 + t359) * t261) * m(6)) * pkin(4) - m(6) * (t35 * t37 + t36 * t38) + t308 * g(2) + t296 * g(1) + (t152 + t342) * t101 + t290 - t38 * t102 - t37 * t103 - t100 * t151; -t124 * (mrSges(6,1) * t117 + t384) + (-t347 + t383) * t374 + t53 * t373 + (-Ifges(6,6) * t117 + t382) * t370 - t35 * t102 + t36 * t103 - g(1) * t310 - g(2) * t326 - t306 * t359 + (t117 * t36 + t381) * mrSges(6,3) + t301 + (-Ifges(6,2) * t117 + t393) * t375;];
tau = t1;
