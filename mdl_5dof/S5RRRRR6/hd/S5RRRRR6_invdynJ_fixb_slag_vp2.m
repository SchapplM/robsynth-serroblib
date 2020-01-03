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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:14:47
% EndTime: 2020-01-03 12:14:58
% DurationCPUTime: 5.27s
% Computational Cost: add. (8737->494), mult. (13204->659), div. (0->0), fcn. (8699->16), ass. (0->261)
t287 = sin(qJ(4));
t288 = sin(qJ(3));
t292 = cos(qJ(4));
t293 = cos(qJ(3));
t208 = -t287 * t288 + t292 * t293;
t296 = -pkin(8) - pkin(7);
t323 = qJD(3) * t296;
t218 = t288 * t323;
t219 = t293 * t323;
t234 = t296 * t288;
t275 = t293 * pkin(8);
t235 = pkin(7) * t293 + t275;
t294 = cos(qJ(2));
t359 = qJD(1) * pkin(1);
t326 = t294 * t359;
t330 = qJD(4) * t292;
t331 = qJD(4) * t287;
t394 = -t208 * t326 + t292 * t218 + t287 * t219 + t234 * t330 - t235 * t331;
t160 = t287 * t234 + t292 * t235;
t209 = t287 * t293 + t288 * t292;
t393 = -qJD(4) * t160 + t209 * t326 - t218 * t287 + t292 * t219;
t233 = -mrSges(4,1) * t293 + mrSges(4,2) * t288;
t412 = -mrSges(3,1) + t233;
t305 = t208 * qJD(4);
t141 = qJD(3) * t208 + t305;
t377 = pkin(9) * t141;
t411 = -t377 + t393;
t306 = t209 * qJD(4);
t142 = -qJD(3) * t209 - t306;
t139 = t142 * pkin(9);
t410 = -t139 - t394;
t383 = t288 / 0.2e1;
t280 = qJD(1) + qJD(2);
t180 = t208 * t280;
t181 = t209 * t280;
t286 = sin(qJ(5));
t291 = cos(qJ(5));
t319 = t291 * t180 - t181 * t286;
t113 = Ifges(6,4) * t319;
t117 = t180 * t286 + t181 * t291;
t279 = qJD(3) + qJD(4);
t268 = qJD(5) + t279;
t54 = Ifges(6,1) * t117 + Ifges(6,5) * t268 + t113;
t409 = t113 + t54;
t259 = pkin(3) * t292 + pkin(4);
t328 = qJD(5) * t291;
t329 = qJD(5) * t286;
t340 = t287 * t291;
t289 = sin(qJ(2));
t327 = t289 * t359;
t223 = pkin(7) * t280 + t327;
t321 = pkin(8) * t280 + t223;
t166 = t321 * t293;
t155 = t292 * t166;
t165 = t321 * t288;
t104 = t165 * t287 - t155;
t376 = pkin(9) * t180;
t79 = t104 - t376;
t153 = t287 * t166;
t105 = -t292 * t165 - t153;
t174 = t181 * pkin(9);
t80 = -t174 + t105;
t408 = t286 * t80 - t291 * t79 - t259 * t329 + (-t287 * t328 + (-t286 * t292 - t340) * qJD(4)) * pkin(3);
t341 = t286 * t287;
t407 = -t286 * t79 - t291 * t80 + t259 * t328 + (-t287 * t329 + (t291 * t292 - t341) * qJD(4)) * pkin(3);
t343 = t280 * t288;
t221 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t343;
t342 = t280 * t293;
t222 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t342;
t406 = (mrSges(3,2) * t280 + t288 * t221 - t293 * t222) * t294;
t284 = qJ(3) + qJ(4);
t273 = qJ(5) + t284;
t257 = cos(t273);
t241 = t257 * mrSges(6,1);
t271 = cos(t284);
t405 = -t271 * mrSges(5,1) - t241;
t269 = sin(t284);
t256 = sin(t273);
t356 = t256 * mrSges(6,2);
t313 = mrSges(5,2) * t269 + t356;
t404 = t313 + t412;
t403 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t224 = -pkin(2) * t280 - t326;
t366 = mrSges(4,2) * t293;
t315 = mrSges(4,1) * t288 + t366;
t402 = t224 * t315 + qJD(3) * (Ifges(4,5) * t293 - Ifges(4,6) * t288) / 0.2e1;
t361 = Ifges(6,4) * t117;
t53 = Ifges(6,2) * t319 + Ifges(6,6) * t268 + t361;
t401 = t53 / 0.2e1;
t389 = -t319 / 0.2e1;
t400 = mrSges(6,2) * t319;
t399 = Ifges(6,1) * t319;
t398 = Ifges(6,5) * t319;
t358 = qJD(3) * pkin(3);
t156 = -t165 + t358;
t101 = t156 * t287 + t155;
t77 = t101 + t376;
t355 = t286 * t77;
t100 = t292 * t156 - t153;
t76 = t100 - t174;
t72 = pkin(4) * t279 + t76;
t35 = t291 * t72 - t355;
t397 = t319 * t35;
t159 = t292 * t234 - t235 * t287;
t375 = pkin(9) * t209;
t125 = t159 - t375;
t201 = t208 * pkin(9);
t126 = t201 + t160;
t67 = t125 * t291 - t126 * t286;
t396 = qJD(5) * t67 + t411 * t286 - t410 * t291;
t68 = t125 * t286 + t126 * t291;
t395 = -qJD(5) * t68 + t410 * t286 + t411 * t291;
t258 = pkin(1) * t289 + pkin(7);
t371 = -pkin(8) - t258;
t202 = t371 * t288;
t345 = t258 * t293;
t203 = t275 + t345;
t135 = t287 * t202 + t292 * t203;
t392 = t412 * t280;
t390 = t224 * t289 + (t288 ^ 2 + t293 ^ 2) * t223 * t294;
t388 = -t117 / 0.2e1;
t387 = t117 / 0.2e1;
t385 = t181 / 0.2e1;
t384 = -t268 / 0.2e1;
t378 = pkin(3) * t288;
t381 = m(6) * (-pkin(4) * t269 - t378);
t354 = t291 * t77;
t36 = t286 * t72 + t354;
t380 = mrSges(6,3) * t36;
t379 = pkin(1) * t294;
t374 = g(1) * t293;
t285 = qJ(1) + qJ(2);
t270 = sin(t285);
t373 = g(2) * t270;
t272 = cos(t285);
t372 = g(3) * t272;
t369 = mrSges(5,1) * t269;
t368 = mrSges(6,1) * t256;
t365 = mrSges(5,3) * t180;
t364 = Ifges(4,4) * t288;
t363 = Ifges(4,4) * t293;
t362 = Ifges(5,4) * t181;
t360 = pkin(1) * qJD(2);
t357 = t181 * mrSges(5,3);
t353 = t293 * Ifges(4,2);
t352 = qJDD(1) * pkin(1);
t324 = qJD(2) * t359;
t206 = t289 * t352 + t294 * t324;
t278 = qJDD(1) + qJDD(2);
t184 = pkin(7) * t278 + t206;
t333 = qJD(3) * t288;
t127 = t293 * t184 - t223 * t333;
t351 = t127 * t293;
t332 = qJD(3) * t293;
t128 = -t184 * t288 - t223 * t332;
t350 = t128 * t288;
t346 = t257 * t272;
t344 = t271 * t272;
t255 = pkin(4) * t271;
t260 = pkin(3) * t293 + pkin(2);
t217 = t255 + t260;
t281 = -pkin(9) + t296;
t337 = t270 * t217 + t272 * t281;
t336 = mrSges(6,2) * t346 + t272 * t368;
t335 = t270 * t260 + t272 * t296;
t334 = t272 * pkin(2) + t270 * pkin(7);
t277 = qJDD(3) + qJDD(4);
t325 = t294 * t360;
t263 = pkin(3) * t333;
t129 = -pkin(4) * t142 + t263;
t320 = qJD(3) * t371;
t134 = t292 * t202 - t203 * t287;
t318 = t272 * t217 - t270 * t281;
t317 = t272 * t260 - t270 * t296;
t316 = -mrSges(5,2) * t344 - t272 * t369 - t336;
t205 = -t289 * t324 + t294 * t352;
t314 = -mrSges(6,2) * t257 - t368;
t312 = t353 + t364;
t108 = t134 - t375;
t109 = t201 + t135;
t50 = t108 * t291 - t109 * t286;
t51 = t108 * t286 + t109 * t291;
t137 = t208 * t291 - t209 * t286;
t138 = t208 * t286 + t209 * t291;
t310 = t221 * t293 + t222 * t288;
t177 = -pkin(4) * t208 - t260;
t267 = qJDD(5) + t277;
t195 = t278 * t288 + t280 * t332;
t96 = qJDD(3) * pkin(3) - pkin(8) * t195 + t128;
t194 = t278 * t293 - t280 * t333;
t99 = pkin(8) * t194 + t127;
t29 = -qJD(4) * t101 - t287 * t99 + t292 * t96;
t93 = t194 * t287 + t195 * t292 + t280 * t305;
t15 = pkin(4) * t277 - pkin(9) * t93 + t29;
t28 = t156 * t330 - t166 * t331 + t287 * t96 + t292 * t99;
t94 = t194 * t292 - t195 * t287 - t280 * t306;
t19 = pkin(9) * t94 + t28;
t3 = qJD(5) * t35 + t15 * t286 + t19 * t291;
t32 = qJD(5) * t319 + t286 * t94 + t291 * t93;
t33 = -qJD(5) * t117 - t286 * t93 + t291 * t94;
t4 = -qJD(5) * t36 + t15 * t291 - t19 * t286;
t309 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t32 + Ifges(6,6) * t33 + Ifges(6,3) * t267;
t307 = t288 * (Ifges(4,1) * t293 - t364);
t161 = t288 * t320 + t293 * t325;
t162 = -t288 * t325 + t293 * t320;
t60 = t292 * t161 + t287 * t162 + t202 * t330 - t203 * t331;
t183 = -pkin(2) * t278 - t205;
t304 = t313 + t405;
t182 = -t260 * t280 - t326;
t130 = -pkin(3) * t194 + t183;
t302 = mrSges(5,2) * t271 - t314 + t369;
t61 = -qJD(4) * t135 - t161 * t287 + t292 * t162;
t301 = -qJD(3) * t310 + t293 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t194) - t288 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t195);
t300 = -mrSges(5,1) * t344 - mrSges(6,1) * t346 + t403 * t270 + t404 * t272;
t299 = (m(4) * pkin(7) - t403) * t272 + (t404 + t405) * t270;
t111 = Ifges(5,2) * t180 + Ifges(5,6) * t279 + t362;
t173 = Ifges(5,4) * t180;
t112 = Ifges(5,1) * t181 + Ifges(5,5) * t279 + t173;
t124 = -pkin(4) * t180 + t182;
t298 = mrSges(6,3) * t397 + t29 * mrSges(5,1) - t28 * mrSges(5,2) - t182 * (mrSges(5,1) * t181 + mrSges(5,2) * t180) - t124 * t400 - t279 * (Ifges(5,5) * t180 - Ifges(5,6) * t181) / 0.2e1 + Ifges(5,3) * t277 + t309 + t399 * t388 + t398 * t384 + t100 * t365 + t111 * t385 - t181 * (Ifges(5,1) * t180 - t362) / 0.2e1 + Ifges(5,6) * t94 + Ifges(5,5) * t93 + t409 * t389 - (-Ifges(5,2) * t181 + t112 + t173) * t180 / 0.2e1 + (-t124 * mrSges(6,1) - Ifges(6,4) * t388 - Ifges(6,2) * t389 - Ifges(6,6) * t384 + t380 + t401) * t117;
t178 = Ifges(4,6) * qJD(3) + t312 * t280;
t239 = Ifges(4,4) * t342;
t179 = Ifges(4,1) * t343 + Ifges(4,5) * qJD(3) + t239;
t55 = qJD(5) * t137 + t141 * t291 + t142 * t286;
t56 = -qJD(5) * t138 - t141 * t286 + t142 * t291;
t57 = -pkin(4) * t94 + t130;
t297 = (t130 * mrSges(5,2) - t29 * mrSges(5,3) + Ifges(5,1) * t93 + Ifges(5,4) * t94 + Ifges(5,5) * t277) * t209 + (-t130 * mrSges(5,1) + t28 * mrSges(5,3) + Ifges(5,4) * t93 + Ifges(5,2) * t94 + Ifges(5,6) * t277) * t208 + (t179 + t280 * (-Ifges(4,2) * t288 + t363)) * t332 / 0.2e1 + t56 * t401 + (mrSges(6,2) * t57 - mrSges(6,3) * t4 + Ifges(6,1) * t32 + Ifges(6,4) * t33 + Ifges(6,5) * t267) * t138 + t319 * (Ifges(6,4) * t55 + Ifges(6,2) * t56) / 0.2e1 + (-t100 * t141 + t101 * t142) * mrSges(5,3) + (0.2e1 * Ifges(4,5) * t383 + Ifges(4,6) * t293) * qJDD(3) + (t307 * t280 / 0.2e1 + t402) * qJD(3) + t56 * t380 + (Ifges(5,1) * t141 + Ifges(5,4) * t142) * t385 + (Ifges(6,1) * t55 + Ifges(6,4) * t56) * t387 + mrSges(4,3) * t351 + t194 * t312 / 0.2e1 + (-mrSges(6,1) * t57 + mrSges(6,3) * t3 + Ifges(6,4) * t32 + Ifges(6,2) * t33 + Ifges(6,6) * t267) * t137 + (Ifges(4,1) * t195 + Ifges(4,4) * t194) * t383 - t35 * t55 * mrSges(6,3) - t178 * t333 / 0.2e1 + t195 * (Ifges(4,1) * t288 + t363) / 0.2e1 + t293 * (Ifges(4,4) * t195 + Ifges(4,2) * t194) / 0.2e1 + t55 * t54 / 0.2e1 + t124 * (-mrSges(6,1) * t56 + mrSges(6,2) * t55) + t141 * t112 / 0.2e1 + t142 * t111 / 0.2e1 + t180 * (Ifges(5,4) * t141 + Ifges(5,2) * t142) / 0.2e1 + t182 * (-mrSges(5,1) * t142 + mrSges(5,2) * t141) + t205 * mrSges(3,1) - t206 * mrSges(3,2) + t183 * t233 + t268 * (Ifges(6,5) * t55 + Ifges(6,6) * t56) / 0.2e1 + Ifges(3,3) * t278 + t279 * (Ifges(5,5) * t141 + Ifges(5,6) * t142) / 0.2e1;
t295 = cos(qJ(1));
t290 = sin(qJ(1));
t276 = t295 * pkin(1);
t274 = t290 * pkin(1);
t264 = t289 * t360;
t261 = -pkin(2) - t379;
t253 = t270 * pkin(2);
t232 = -t260 - t379;
t220 = t264 + t263;
t192 = pkin(3) * t340 + t259 * t286;
t191 = -pkin(3) * t341 + t259 * t291;
t164 = t177 - t379;
t152 = mrSges(5,1) * t279 - t357;
t151 = -mrSges(5,2) * t279 + t365;
t143 = pkin(3) * t343 + pkin(4) * t181;
t131 = -mrSges(4,1) * t194 + mrSges(4,2) * t195;
t122 = -mrSges(5,1) * t180 + mrSges(5,2) * t181;
t120 = t129 + t264;
t103 = mrSges(6,1) * t268 - mrSges(6,3) * t117;
t102 = -mrSges(6,2) * t268 + mrSges(6,3) * t319;
t85 = -mrSges(5,2) * t277 + mrSges(5,3) * t94;
t84 = mrSges(5,1) * t277 - mrSges(5,3) * t93;
t59 = -mrSges(6,1) * t319 + mrSges(6,2) * t117;
t49 = t61 - t377;
t48 = t139 + t60;
t47 = -mrSges(5,1) * t94 + mrSges(5,2) * t93;
t38 = t291 * t76 - t355;
t37 = -t286 * t76 - t354;
t26 = -mrSges(6,2) * t267 + mrSges(6,3) * t33;
t25 = mrSges(6,1) * t267 - mrSges(6,3) * t32;
t9 = -mrSges(6,1) * t33 + mrSges(6,2) * t32;
t8 = -qJD(5) * t51 - t286 * t48 + t291 * t49;
t7 = qJD(5) * t50 + t286 * t49 + t291 * t48;
t1 = [m(6) * (t120 * t124 + t164 * t57 + t3 * t51 + t35 * t8 + t36 * t7 + t4 * t50) + m(5) * (t100 * t61 + t101 * t60 + t130 * t232 + t134 * t29 + t135 * t28 + t182 * t220) + m(4) * (t127 * t345 + t183 * t261 - t258 * t350) + t301 * t258 + t297 + (-m(6) * (t276 + t318) - m(4) * (t276 + t334) - mrSges(2,1) * t295 + mrSges(2,2) * t290 - m(5) * (t276 + t317) + t300) * g(2) + (-m(6) * (t274 + t337) - m(5) * (t274 + t335) - m(4) * (t253 + t274) - mrSges(2,1) * t290 - mrSges(2,2) * t295 + t299) * g(3) - mrSges(4,3) * t350 + t50 * t25 + t51 * t26 + t7 * t102 + t8 * t103 + t120 * t59 + t134 * t84 + t135 * t85 + t60 * t151 + t61 * t152 + t164 * t9 + t220 * t122 + Ifges(2,3) * qJDD(1) + t232 * t47 + t261 * t131 + ((mrSges(3,1) * t294 - mrSges(3,2) * t289) * t278 + (-g(2) * t295 - g(3) * t290 + t205 * t294 + t206 * t289) * m(3) + (m(4) * t390 + t289 * t392 - t406) * qJD(2)) * pkin(1); t393 * t152 + (-t128 * mrSges(4,3) + t122 * t358) * t288 + t394 * t151 + t297 + (t406 + (-t122 - t392 - t59) * t289) * t359 + t300 * g(2) + t301 * pkin(7) + t299 * g(3) + t67 * t25 + t68 * t26 + t129 * t59 - pkin(2) * t131 + t159 * t84 + t160 * t85 + t177 * t9 - t260 * t47 + t396 * t102 + t395 * t103 + (-t318 * g(2) - t337 * g(3) + t177 * t57 + t3 * t68 + t4 * t67 + t396 * t36 + t395 * t35 + (t129 - t327) * t124) * m(6) + (-t317 * g(2) - t335 * g(3) - t130 * t260 + t159 * t29 + t160 * t28 + (t263 - t327) * t182 + t394 * t101 + t393 * t100) * m(5) + (-pkin(2) * t183 + (-t350 + t351) * pkin(7) - t334 * g(2) - t253 * g(3) - t390 * t359) * m(4); t310 * t223 + (t233 + t304) * g(1) - m(5) * (t100 * t104 + t101 * t105) + t298 + ((-t315 + t381) * t272 + t316) * g(3) + (-t381 + t366 + (m(5) * pkin(3) + mrSges(4,1)) * t288 + t302) * t373 + t101 * t357 + Ifges(4,3) * qJDD(3) - t127 * mrSges(4,2) + t128 * mrSges(4,1) - t143 * t59 - t105 * t151 - t104 * t152 + t191 * t25 + t192 * t26 + Ifges(4,6) * t194 + Ifges(4,5) * t195 + t407 * t102 + t408 * t103 + (t178 * t383 + (t353 * t383 - t307 / 0.2e1) * t280 + (-m(5) * t182 - t122) * t378 - (t179 + t239) * t293 / 0.2e1 - t402) * t280 + (t287 * t85 + t292 * t84 + (t292 * t151 - t287 * t152) * qJD(4) + (-t100 * t331 + t101 * t330 + t28 * t287 - t288 * t372 + t29 * t292 - t374) * m(5)) * pkin(3) + (-t374 * pkin(3) - g(1) * t255 - t124 * t143 + t191 * t4 + t192 * t3 + t408 * t35 + t407 * t36) * m(6); t302 * t373 + t316 * g(3) + t304 * g(1) + (t152 + t357) * t101 + t298 - t38 * t102 - t37 * t103 - t100 * t151 + (-t181 * t59 + t291 * t25 + t286 * t26 + (t291 * t102 - t286 * t103) * qJD(5) + (-t124 * t181 - t35 * t329 + t36 * t328 + t286 * t3 + t291 * t4 - g(1) * t271 + (-t372 + t373) * t269) * m(6)) * pkin(4) - m(6) * (t35 * t37 + t36 * t38); -t124 * (mrSges(6,1) * t117 + t400) + (-t361 + t399) * t388 + t53 * t387 + (-Ifges(6,6) * t117 + t398) * t384 - t35 * t102 + t36 * t103 - g(1) * (t241 - t356) - g(3) * t336 - t314 * t373 + (t117 * t36 + t397) * mrSges(6,3) + t309 + (-Ifges(6,2) * t117 + t409) * t389;];
tau = t1;
