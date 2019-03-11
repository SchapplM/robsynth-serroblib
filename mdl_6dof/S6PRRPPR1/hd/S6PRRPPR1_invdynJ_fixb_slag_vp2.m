% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:34
% EndTime: 2019-03-08 20:59:12
% DurationCPUTime: 23.35s
% Computational Cost: add. (8779->717), mult. (20496->986), div. (0->0), fcn. (16002->18), ass. (0->325)
t254 = pkin(12) + qJ(6);
t250 = sin(t254);
t252 = cos(t254);
t256 = sin(pkin(12));
t259 = cos(pkin(12));
t293 = -mrSges(6,1) * t259 + mrSges(6,2) * t256;
t375 = t259 * pkin(5);
t410 = mrSges(5,1) + m(7) * (pkin(4) + t375) + mrSges(7,1) * t252 - mrSges(7,2) * t250 + m(6) * pkin(4) - t293;
t265 = sin(qJ(2));
t258 = sin(pkin(6));
t343 = qJD(1) * t258;
t320 = t265 * t343;
t264 = sin(qJ(3));
t337 = qJD(3) * t264;
t416 = pkin(3) * t337 - t320;
t405 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(9) - qJ(5));
t257 = sin(pkin(11));
t267 = cos(qJ(3));
t362 = cos(pkin(11));
t214 = t257 * t267 + t264 * t362;
t200 = t214 * qJD(3);
t304 = t362 * t267;
t280 = -t257 * t264 + t304;
t202 = t280 * qJD(3);
t459 = pkin(4) * t200 - qJ(5) * t202 - qJD(5) * t214 + t416;
t261 = -qJ(4) - pkin(8);
t310 = qJD(3) * t261;
t193 = qJD(4) * t267 + t264 * t310;
t194 = -qJD(4) * t264 + t267 * t310;
t268 = cos(qJ(2));
t319 = t268 * t343;
t419 = t193 * t362 + t257 * t194 - t280 * t319;
t255 = qJ(3) + pkin(11);
t251 = sin(t255);
t253 = cos(t255);
t458 = (-mrSges(7,3) + t405) * t251 - t410 * t253;
t427 = -t256 * t419 + t259 * t459;
t426 = t256 * t459 + t259 * t419;
t443 = -m(7) - m(6);
t412 = -m(5) + t443;
t455 = pkin(3) * t412 - mrSges(4,1);
t351 = t202 * t259;
t454 = pkin(5) * t200 - pkin(9) * t351 + t427;
t352 = t202 * t256;
t453 = pkin(9) * t352 - t426;
t449 = Ifges(5,5) * qJD(3);
t448 = Ifges(5,6) * qJD(3);
t341 = qJD(2) * t258;
t313 = qJD(1) * t341;
t234 = t268 * t313;
t333 = qJDD(1) * t258;
t189 = t265 * t333 + t234;
t179 = qJDD(2) * pkin(8) + t189;
t260 = cos(pkin(6));
t342 = qJD(1) * t260;
t447 = qJD(3) * t342 + t179;
t292 = mrSges(6,1) * t256 + mrSges(6,2) * t259;
t446 = -m(7) * pkin(5) * t256 - t250 * mrSges(7,1) - t252 * mrSges(7,2) - mrSges(5,3) - t292;
t348 = t258 * t265;
t205 = t260 * t267 - t264 * t348;
t201 = t214 * qJD(2);
t169 = qJD(3) * t256 + t201 * t259;
t263 = sin(qJ(6));
t266 = cos(qJ(6));
t300 = t259 * qJD(3) - t201 * t256;
t444 = -t169 * t263 + t266 * t300;
t98 = t169 * t266 + t263 * t300;
t335 = qJD(2) * qJD(3);
t219 = qJDD(2) * t267 - t264 * t335;
t220 = qJDD(2) * t264 + t267 * t335;
t156 = t257 * t219 + t220 * t362;
t125 = qJDD(3) * t259 - t156 * t256;
t126 = qJDD(3) * t256 + t156 * t259;
t34 = qJD(6) * t444 + t125 * t263 + t126 * t266;
t401 = t34 / 0.2e1;
t35 = -qJD(6) * t98 + t125 * t266 - t126 * t263;
t400 = t35 / 0.2e1;
t340 = qJD(2) * t264;
t199 = -qJD(2) * t304 + t257 * t340;
t247 = pkin(3) * t267 + pkin(2);
t186 = -qJD(2) * t247 + qJD(4) - t319;
t108 = pkin(4) * t199 - qJ(5) * t201 + t186;
t241 = t267 * t342;
t223 = qJD(2) * pkin(8) + t320;
t299 = qJ(4) * qJD(2) + t223;
t160 = -t264 * t299 + t241;
t154 = qJD(3) * pkin(3) + t160;
t318 = t264 * t342;
t161 = t267 * t299 + t318;
t305 = t362 * t161;
t86 = t257 * t154 + t305;
t81 = qJD(3) * qJ(5) + t86;
t41 = t259 * t108 - t256 * t81;
t23 = pkin(5) * t199 - pkin(9) * t169 + t41;
t42 = t256 * t108 + t259 * t81;
t29 = pkin(9) * t300 + t42;
t6 = t23 * t266 - t263 * t29;
t442 = t6 * mrSges(7,1);
t7 = t23 * t263 + t266 * t29;
t441 = t7 * mrSges(7,2);
t392 = t125 / 0.2e1;
t391 = t126 / 0.2e1;
t155 = -t362 * t219 + t220 * t257;
t149 = qJDD(6) + t155;
t390 = t149 / 0.2e1;
t389 = t155 / 0.2e1;
t440 = t219 / 0.2e1;
t349 = t214 * t259;
t140 = -pkin(4) * t280 - qJ(5) * t214 - t247;
t229 = t261 * t264;
t230 = t261 * t267;
t166 = t257 * t229 - t230 * t362;
t82 = t259 * t140 - t166 * t256;
t56 = -pkin(5) * t280 - pkin(9) * t349 + t82;
t350 = t214 * t256;
t83 = t256 * t140 + t259 * t166;
t60 = -pkin(9) * t350 + t83;
t19 = -t263 * t60 + t266 * t56;
t439 = qJD(6) * t19 + t263 * t454 - t266 * t453;
t20 = t263 * t56 + t266 * t60;
t438 = -qJD(6) * t20 + t263 * t453 + t266 * t454;
t436 = t41 * mrSges(6,1);
t435 = t42 * mrSges(6,2);
t190 = qJD(6) + t199;
t432 = t300 * Ifges(6,6);
t433 = t169 * Ifges(6,5);
t434 = t98 * Ifges(7,5) + Ifges(7,6) * t444 + t199 * Ifges(6,3) + t190 * Ifges(7,3) + t432 + t433;
t150 = t257 * t161;
t85 = t154 * t362 - t150;
t77 = -qJD(3) * pkin(4) + qJD(5) - t85;
t431 = t77 * t292;
t364 = qJDD(3) / 0.2e1;
t378 = pkin(3) * t257;
t243 = qJ(5) + t378;
t374 = pkin(9) + t243;
t209 = t374 * t256;
t210 = t374 * t259;
t136 = -t209 * t263 + t210 * t266;
t215 = t256 * t266 + t259 * t263;
t353 = t199 * t259;
t330 = pkin(3) * t340;
t122 = pkin(4) * t201 + qJ(5) * t199 + t330;
t91 = t160 * t362 - t150;
t50 = t259 * t122 - t256 * t91;
t30 = pkin(5) * t201 + pkin(9) * t353 + t50;
t354 = t199 * t256;
t51 = t256 * t122 + t259 * t91;
t40 = pkin(9) * t354 + t51;
t429 = -qJD(5) * t215 - qJD(6) * t136 + t263 * t40 - t266 * t30;
t135 = -t209 * t266 - t210 * t263;
t284 = t256 * t263 - t259 * t266;
t428 = -qJD(5) * t284 + qJD(6) * t135 - t263 * t30 - t266 * t40;
t138 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t156;
t64 = -t125 * mrSges(6,1) + t126 * mrSges(6,2);
t425 = t64 - t138;
t372 = mrSges(5,3) * t201;
t422 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t300 + mrSges(6,2) * t169 + t372;
t116 = t215 * t199;
t204 = t215 * qJD(6);
t421 = t116 + t204;
t117 = t284 * t199;
t203 = t284 * qJD(6);
t420 = t117 + t203;
t111 = -mrSges(6,2) * t199 + mrSges(6,3) * t300;
t112 = mrSges(6,1) * t199 - mrSges(6,3) * t169;
t415 = t111 * t259 - t112 * t256;
t332 = qJDD(1) * t260;
t100 = -t223 * t337 + t264 * t332 + t267 * t447;
t173 = t223 * t267 + t318;
t239 = t267 * t332;
t101 = -t173 * qJD(3) - t179 * t264 + t239;
t414 = t100 * t267 - t101 * t264;
t334 = qJD(2) * qJD(4);
t336 = qJD(3) * t267;
t70 = -t223 * t336 + qJDD(3) * pkin(3) - qJ(4) * t220 + t239 + (-t334 - t447) * t264;
t73 = qJ(4) * t219 + t267 * t334 + t100;
t28 = t257 * t70 + t362 * t73;
t25 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t28;
t233 = t265 * t313;
t188 = t268 * t333 - t233;
t178 = -qJDD(2) * pkin(2) - t188;
t141 = -pkin(3) * t219 + qJDD(4) + t178;
t47 = pkin(4) * t155 - qJ(5) * t156 - qJD(5) * t201 + t141;
t12 = -t25 * t256 + t259 * t47;
t13 = t259 * t25 + t256 * t47;
t413 = -t12 * t256 + t13 * t259;
t411 = 0.2e1 * t364;
t409 = m(5) * t85 - t422;
t5 = pkin(5) * t155 - pkin(9) * t126 + t12;
t8 = pkin(9) * t125 + t13;
t1 = qJD(6) * t6 + t263 * t5 + t266 * t8;
t2 = -qJD(6) * t7 - t263 * t8 + t266 * t5;
t408 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t327 = m(4) * pkin(8) + mrSges(4,3);
t407 = mrSges(3,2) - t327 + t446;
t295 = -mrSges(4,1) * t267 + mrSges(4,2) * t264;
t279 = m(4) * pkin(2) - t295;
t406 = mrSges(3,1) + t279 - t458;
t269 = qJD(2) ^ 2;
t404 = Ifges(7,4) * t401 + Ifges(7,2) * t400 + Ifges(7,6) * t390;
t403 = m(5) * pkin(3);
t402 = Ifges(7,1) * t401 + Ifges(7,4) * t400 + Ifges(7,5) * t390;
t379 = Ifges(7,4) * t98;
t38 = Ifges(7,2) * t444 + Ifges(7,6) * t190 + t379;
t399 = t38 / 0.2e1;
t95 = Ifges(7,4) * t444;
t39 = Ifges(7,1) * t98 + Ifges(7,5) * t190 + t95;
t398 = t39 / 0.2e1;
t397 = Ifges(6,1) * t391 + Ifges(6,4) * t392 + Ifges(6,5) * t389;
t396 = -t444 / 0.2e1;
t395 = t444 / 0.2e1;
t394 = -t98 / 0.2e1;
t393 = t98 / 0.2e1;
t388 = -t190 / 0.2e1;
t387 = t190 / 0.2e1;
t386 = -t199 / 0.2e1;
t385 = t199 / 0.2e1;
t382 = t201 / 0.2e1;
t380 = t259 / 0.2e1;
t376 = g(3) * t258;
t373 = mrSges(5,3) * t199;
t371 = Ifges(4,4) * t264;
t370 = Ifges(4,4) * t267;
t369 = Ifges(5,4) * t201;
t368 = Ifges(6,4) * t256;
t367 = Ifges(6,4) * t259;
t363 = cos(pkin(10));
t361 = sin(pkin(10));
t347 = t258 * t268;
t339 = qJD(2) * t265;
t338 = qJD(2) * t267;
t331 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t149;
t46 = -mrSges(7,1) * t444 + mrSges(7,2) * t98;
t328 = t46 + t422;
t326 = mrSges(4,3) * t340;
t325 = mrSges(4,3) * t338;
t323 = -t256 * (t169 * Ifges(6,4) + Ifges(6,2) * t300 + t199 * Ifges(6,6)) / 0.2e1;
t322 = (t169 * Ifges(6,1) + Ifges(6,4) * t300 + t199 * Ifges(6,5)) * t380;
t321 = t362 * pkin(3);
t317 = t258 * t339;
t316 = t268 * t341;
t11 = -t35 * mrSges(7,1) + t34 * mrSges(7,2);
t309 = t258 * t363;
t308 = t258 * t361;
t307 = t363 * t265;
t306 = t363 * t268;
t303 = t361 * t265;
t302 = t361 * t268;
t87 = t155 * mrSges(5,1) + t156 * mrSges(5,2);
t89 = t160 * t257 + t305;
t123 = t193 * t257 - t362 * t194;
t165 = -t362 * t229 - t230 * t257;
t246 = -t321 - pkin(4);
t27 = -t257 * t73 + t362 * t70;
t290 = Ifges(6,1) * t259 - t368;
t289 = t267 * Ifges(4,2) + t371;
t288 = -Ifges(6,2) * t256 + t367;
t287 = Ifges(4,5) * t267 - Ifges(4,6) * t264;
t286 = Ifges(6,5) * t259 - Ifges(6,6) * t256;
t285 = t256 * t41 - t259 * t42;
t206 = t260 * t264 + t267 * t348;
t130 = t257 * t205 + t206 * t362;
t109 = -t130 * t256 - t259 * t347;
t110 = t130 * t259 - t256 * t347;
t52 = t109 * t266 - t110 * t263;
t53 = t109 * t263 + t110 * t266;
t224 = -qJD(2) * pkin(2) - t319;
t282 = t224 * (mrSges(4,1) * t264 + mrSges(4,2) * t267);
t281 = t264 * (Ifges(4,1) * t267 - t371);
t196 = t260 * t307 + t302;
t144 = t196 * t251 + t253 * t309;
t198 = -t260 * t303 + t306;
t146 = t198 * t251 - t253 * t308;
t181 = t251 * t348 - t260 * t253;
t277 = -g(1) * t146 - g(2) * t144 - g(3) * t181;
t195 = -t260 * t306 + t303;
t197 = t260 * t302 + t307;
t273 = -g(1) * t197 - g(2) * t195 + g(3) * t347;
t26 = -qJDD(3) * pkin(4) + qJDD(5) - t27;
t248 = Ifges(4,4) * t338;
t228 = -qJD(3) * mrSges(4,2) + t325;
t227 = qJD(3) * mrSges(4,1) - t326;
t225 = t246 - t375;
t217 = t295 * qJD(2);
t208 = Ifges(4,1) * t340 + Ifges(4,5) * qJD(3) + t248;
t207 = Ifges(4,6) * qJD(3) + qJD(2) * t289;
t192 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t220;
t191 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t219;
t187 = Ifges(5,4) * t199;
t182 = t251 * t260 + t253 * t348;
t176 = -qJD(3) * mrSges(5,2) - t373;
t172 = -t223 * t264 + t241;
t162 = -mrSges(4,1) * t219 + mrSges(4,2) * t220;
t159 = qJD(3) * t205 + t267 * t316;
t158 = -qJD(3) * t206 - t264 * t316;
t147 = t198 * t253 + t251 * t308;
t145 = t196 * t253 - t251 * t309;
t139 = mrSges(5,1) * t199 + mrSges(5,2) * t201;
t137 = -qJDD(3) * mrSges(5,2) - mrSges(5,3) * t155;
t132 = t201 * Ifges(5,1) - t187 + t449;
t131 = -t199 * Ifges(5,2) + t369 + t448;
t129 = -t205 * t362 + t206 * t257;
t128 = t284 * t214;
t127 = t215 * t214;
t118 = pkin(5) * t350 + t165;
t92 = pkin(5) * t352 + t123;
t90 = t257 * t158 + t159 * t362;
t88 = -t158 * t362 + t159 * t257;
t75 = mrSges(6,1) * t155 - mrSges(6,3) * t126;
t74 = -mrSges(6,2) * t155 + mrSges(6,3) * t125;
t72 = mrSges(7,1) * t190 - mrSges(7,3) * t98;
t71 = -mrSges(7,2) * t190 + mrSges(7,3) * t444;
t69 = t256 * t317 + t259 * t90;
t68 = -t256 * t90 + t259 * t317;
t63 = -pkin(5) * t354 + t89;
t62 = -t202 * t215 + t203 * t214;
t61 = -t202 * t284 - t204 * t214;
t57 = -pkin(5) * t300 + t77;
t48 = t126 * Ifges(6,4) + t125 * Ifges(6,2) + t155 * Ifges(6,6);
t22 = -mrSges(7,2) * t149 + mrSges(7,3) * t35;
t21 = mrSges(7,1) * t149 - mrSges(7,3) * t34;
t18 = -t125 * pkin(5) + t26;
t17 = -qJD(6) * t53 - t263 * t69 + t266 * t68;
t16 = qJD(6) * t52 + t263 * t68 + t266 * t69;
t3 = [m(2) * qJDD(1) + t109 * t75 + t110 * t74 + t69 * t111 + t68 * t112 + t130 * t137 + t158 * t227 + t159 * t228 + t16 * t71 + t17 * t72 + t90 * t176 + t206 * t191 + t205 * t192 + t52 * t21 + t53 * t22 + t328 * t88 + (t11 + t425) * t129 + (-m(2) - m(3) - m(4) + t412) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t269 - t162 - t87) * t268 + (-mrSges(3,1) * t269 - mrSges(3,2) * qJDD(2) + (t139 + t217) * qJD(2)) * t265) * t258 + m(3) * (qJDD(1) * t260 ^ 2 + (t188 * t268 + t189 * t265) * t258) + m(5) * (-t129 * t27 + t130 * t28 - t85 * t88 + t86 * t90 + (-t141 * t268 + t186 * t339) * t258) + m(4) * (t100 * t206 + t101 * t205 + t158 * t172 + t159 * t173 + (-t178 * t268 + t224 * t339) * t258) + m(7) * (t1 * t53 + t129 * t18 + t16 * t7 + t17 * t6 + t2 * t52 + t57 * t88) + m(6) * (t109 * t12 + t110 * t13 + t129 * t26 + t41 * t68 + t42 * t69 + t77 * t88); ((t458 * t268 + (-t261 * t412 + t446) * t265) * t258 + t412 * t247 * t347) * g(3) + (m(4) * ((-t172 * t267 - t173 * t264) * qJD(3) + t414) - t227 * t336 - t228 * t337 + t267 * t191 - t264 * t192) * pkin(8) + (-t172 * t336 - t173 * t337 + t414) * mrSges(4,3) + (t412 * (-t197 * t247 - t198 * t261) + t407 * t198 + t406 * t197) * g(1) + (t412 * (-t195 * t247 - t196 * t261) + t407 * t196 + t406 * t195) * g(2) + (t141 * mrSges(5,2) - t27 * mrSges(5,3) + Ifges(5,1) * t156 - Ifges(5,4) * t155 + t411 * Ifges(5,5) + t26 * t292 + t286 * t389 + t288 * t392 + t290 * t391 + (-m(6) * t77 - m(7) * t57 + t409 - t46) * t319) * t214 + (t434 / 0.2e1 - t86 * mrSges(5,3) - t435 + t436 + t432 / 0.2e1 + t433 / 0.2e1 + t442 - t441 - Ifges(5,4) * t382 + Ifges(6,3) * t385 - Ifges(5,2) * t386 + Ifges(7,3) * t387 + Ifges(7,5) * t393 + Ifges(7,6) * t395 - t131 / 0.2e1 + t186 * mrSges(5,1) - t448 / 0.2e1) * t200 + (-t85 * mrSges(5,3) + t322 + t323 + t431 + t300 * t288 / 0.2e1 + t169 * t290 / 0.2e1 + Ifges(5,1) * t382 + t286 * t385 + Ifges(5,4) * t386 + t132 / 0.2e1 + t186 * mrSges(5,2) + t449 / 0.2e1) * t202 + (-pkin(2) * t178 - (t224 * t265 + (-t172 * t264 + t173 * t267) * t268) * t343) * m(4) + (-t268 * t376 + t188 + t233) * mrSges(3,1) + (-Ifges(7,1) * t128 - Ifges(7,4) * t127) * t401 + (-t1 * t127 + t128 * t2 - t6 * t61 + t62 * t7) * mrSges(7,3) + (-Ifges(7,4) * t128 - Ifges(7,2) * t127) * t400 + (-Ifges(7,5) * t128 - Ifges(7,6) * t127) * t390 + t18 * (mrSges(7,1) * t127 - mrSges(7,2) * t128) + (t281 + t267 * (-Ifges(4,2) * t264 + t370)) * t335 / 0.2e1 - (t141 * mrSges(5,1) + t12 * mrSges(6,1) - t13 * mrSges(6,2) - t28 * mrSges(5,3) - Ifges(5,4) * t156 + Ifges(6,5) * t391 + Ifges(7,5) * t401 + Ifges(5,2) * t155 - t411 * Ifges(5,6) + Ifges(6,6) * t392 + Ifges(7,6) * t400 + Ifges(6,3) * t389 + Ifges(7,3) * t390 + t408) * t280 - (Ifges(6,5) * t126 + Ifges(6,6) * t125 + Ifges(6,3) * t155 + t331) * t280 / 0.2e1 + (Ifges(7,1) * t61 + Ifges(7,4) * t62) * t393 + t425 * t165 + t426 * t111 + (t12 * t82 + t123 * t77 + t13 * t83 + t165 * t26 + t41 * t427 + t42 * t426) * m(6) + t427 * t112 + t20 * t22 + t19 * t21 + (t265 * t376 - t189 + t234) * mrSges(3,2) + t220 * t370 / 0.2e1 + (-t12 * t349 - t13 * t350 - t351 * t41 - t352 * t42) * mrSges(6,3) + t438 * t72 + t439 * t71 + (t1 * t20 + t118 * t18 + t19 * t2 + t438 * t6 + t439 * t7 + t57 * t92) * m(7) + (Ifges(4,1) * t220 + Ifges(4,4) * t440 + t411 * Ifges(4,5) + t227 * t319) * t264 - t217 * t320 + Ifges(4,6) * t267 * t364 + (Ifges(7,4) * t61 + Ifges(7,2) * t62) * t395 - t48 * t350 / 0.2e1 + t416 * t139 + (-t123 * t85 - t141 * t247 - t165 * t27 + t166 * t28 + t186 * t416 + t419 * t86) * m(5) + t419 * t176 + t422 * t123 + t208 * t336 / 0.2e1 - t207 * t337 / 0.2e1 + t289 * t440 - t267 * t228 * t319 - (t265 * t327 + t268 * t279) * t376 + (t282 + t287 * qJD(3) / 0.2e1) * qJD(3) + t178 * t295 + Ifges(3,3) * qJDD(2) + (Ifges(7,5) * t61 + Ifges(7,6) * t62) * t387 + t57 * (-mrSges(7,1) * t62 + mrSges(7,2) * t61) + t82 * t75 + t83 * t74 + t92 * t46 + t118 * t11 + t349 * t397 + t61 * t398 + t62 * t399 - t128 * t402 - t127 * t404 - pkin(2) * t162 + t166 * t137 - t247 * t87 + t267 * (Ifges(4,4) * t220 + Ifges(4,2) * t219 + Ifges(4,6) * qJDD(3)) / 0.2e1; (-t353 * t41 - t354 * t42 + t413) * mrSges(6,3) + (-t285 * qJD(5) + t243 * t413 + t246 * t26 - t41 * t50 - t42 * t51 - t77 * t89) * m(6) + t415 * qJD(5) - t201 * t436 - (-Ifges(5,1) * t199 - t369 + t434) * t201 / 0.2e1 + (-Ifges(7,5) * t203 - Ifges(7,6) * t204) * t387 + (-Ifges(7,1) * t203 - Ifges(7,4) * t204) * t393 + (-Ifges(7,4) * t203 - Ifges(7,2) * t204) * t395 - (-Ifges(4,2) * t340 + t208 + t248) * t338 / 0.2e1 - t300 * (Ifges(6,6) * t201 - t199 * t288) / 0.2e1 + (-(-t196 * t267 + t264 * t309) * mrSges(4,2) + t405 * t145 + t410 * t144 + (pkin(3) * t443 - mrSges(4,1) - t403) * (-t196 * t264 - t267 * t309)) * g(2) + (Ifges(7,5) * t215 - Ifges(7,6) * t284) * t390 + (Ifges(7,4) * t215 - Ifges(7,2) * t284) * t400 + (Ifges(7,1) * t215 - Ifges(7,4) * t284) * t401 + t18 * (mrSges(7,1) * t284 + mrSges(7,2) * t215) + (-g(1) * t147 - g(2) * t145 - g(3) * t182 - t1 * t284 - t2 * t215 + t420 * t6 - t421 * t7) * mrSges(7,3) - t284 * t404 + t428 * t71 + (t1 * t136 + t135 * t2 + t18 * t225 + t428 * t7 + t429 * t6 - t57 * t63) * m(7) + t429 * t72 - t28 * mrSges(5,2) + t27 * mrSges(5,1) + t409 * t89 - qJD(2) * t282 - t287 * t335 / 0.2e1 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (-t256 * t75 + t259 * t74) * t243 - t85 * t373 - t139 * t330 + t138 * t321 + (-Ifges(5,2) * t201 + t132 - t187) * t385 - t201 * t442 + (t326 + t227) * t173 + (mrSges(7,1) * t421 - mrSges(7,2) * t420) * t57 + t199 * t322 + t199 * t323 + (-(-t198 * t267 - t264 * t308) * mrSges(4,2) + t405 * t147 + t410 * t146 + t455 * (-t198 * t264 + t267 * t308)) * g(1) + (mrSges(4,2) * t206 + t410 * t181 + t405 * t182 + t205 * t455) * g(3) + t207 * t340 / 0.2e1 + t199 * t431 + t201 * t435 - t269 * t281 / 0.2e1 + t201 * t441 + (t325 - t228) * t172 + t26 * t293 - t169 * (Ifges(6,5) * t201 - t199 * t290) / 0.2e1 - m(5) * (t186 * t330 + t86 * t91) - t63 * t46 - t100 * mrSges(4,2) + t101 * mrSges(4,1) - t51 * t111 - t50 * t112 - t116 * t38 / 0.2e1 - t117 * t39 / 0.2e1 + t135 * t21 + t136 * t22 + t86 * t372 + t137 * t378 + t48 * t380 + t131 * t382 + (Ifges(6,3) * t201 - t199 * t286) * t386 + (Ifges(7,5) * t117 + Ifges(7,6) * t116 + Ifges(7,3) * t201) * t388 + (Ifges(6,5) * t256 + Ifges(6,6) * t259) * t389 + (Ifges(6,1) * t256 + t367) * t391 + (Ifges(6,2) * t259 + t368) * t392 + (Ifges(7,1) * t117 + Ifges(7,4) * t116 + Ifges(7,5) * t201) * t394 + (Ifges(7,4) * t117 + Ifges(7,2) * t116 + Ifges(7,6) * t201) * t396 + t256 * t397 - t203 * t398 - t204 * t399 + t215 * t402 + (t257 * t28 + t27 * t362) * t403 - Ifges(5,6) * t155 + Ifges(5,5) * t156 - t91 * t176 - qJD(3) * (-Ifges(5,5) * t199 - Ifges(5,6) * t201) / 0.2e1 - t186 * (mrSges(5,1) * t201 - mrSges(5,2) * t199) + Ifges(4,6) * t219 + Ifges(4,5) * t220 + t225 * t11 + t246 * t64; -t284 * t21 + t215 * t22 + t256 * t74 + t259 * t75 - t421 * t72 - t420 * t71 - t328 * t201 + (t176 + t415) * t199 + t87 + (t1 * t215 - t2 * t284 - t201 * t57 - t420 * t7 - t421 * t6 + t273) * m(7) + (t12 * t259 + t13 * t256 - t199 * t285 - t201 * t77 + t273) * m(6) + (t199 * t86 + t201 * t85 + t141 + t273) * m(5); -t300 * t111 + t169 * t112 - t444 * t71 + t98 * t72 + t11 + t64 + (-t444 * t7 + t6 * t98 + t18 + t277) * m(7) + (t169 * t41 - t300 * t42 + t26 + t277) * m(6); -t57 * (mrSges(7,1) * t98 + mrSges(7,2) * t444) + (Ifges(7,1) * t444 - t379) * t394 + t38 * t393 + (Ifges(7,5) * t444 - Ifges(7,6) * t98) * t388 - t6 * t71 + t7 * t72 - g(1) * ((-t147 * t250 + t197 * t252) * mrSges(7,1) + (-t147 * t252 - t197 * t250) * mrSges(7,2)) - g(2) * ((-t145 * t250 + t195 * t252) * mrSges(7,1) + (-t145 * t252 - t195 * t250) * mrSges(7,2)) - g(3) * ((-t182 * t250 - t252 * t347) * mrSges(7,1) + (-t182 * t252 + t250 * t347) * mrSges(7,2)) + (t444 * t6 + t7 * t98) * mrSges(7,3) + t331 + (-Ifges(7,2) * t98 + t39 + t95) * t396 + t408;];
tau  = t3;
