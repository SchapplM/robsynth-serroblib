% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP2
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
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:58:53
% EndTime: 2019-03-09 05:59:36
% DurationCPUTime: 26.52s
% Computational Cost: add. (9616->722), mult. (20289->945), div. (0->0), fcn. (13396->14), ass. (0->313)
t477 = -mrSges(7,1) - mrSges(6,1);
t476 = -mrSges(7,2) - mrSges(6,2);
t474 = Ifges(6,4) + Ifges(7,4);
t264 = sin(pkin(10));
t245 = pkin(1) * t264 + pkin(7);
t228 = t245 * qJD(1);
t268 = sin(qJ(3));
t272 = cos(qJ(3));
t174 = qJD(2) * t272 - t268 * t228;
t312 = pkin(3) * t268 - pkin(8) * t272;
t217 = t312 * qJD(1);
t267 = sin(qJ(4));
t271 = cos(qJ(4));
t122 = t271 * t174 + t267 * t217;
t350 = qJD(1) * t272;
t328 = t267 * t350;
t274 = -pkin(9) - pkin(8);
t329 = qJD(4) * t274;
t501 = pkin(9) * t328 + t267 * t329 - t122;
t121 = -t174 * t267 + t271 * t217;
t357 = t271 * t272;
t294 = pkin(4) * t268 - pkin(9) * t357;
t500 = -qJD(1) * t294 + t271 * t329 - t121;
t499 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t263 = qJ(4) + qJ(5);
t259 = cos(t263);
t398 = pkin(4) * t271;
t224 = pkin(5) * t259 + t398;
t216 = pkin(3) + t224;
t252 = pkin(3) + t398;
t258 = sin(t263);
t309 = -mrSges(5,1) * t271 + mrSges(5,2) * t267;
t498 = m(5) * pkin(3) + m(6) * t252 + m(7) * t216 + t476 * t258 - t259 * t477 - t309;
t475 = Ifges(6,1) + Ifges(7,1);
t473 = Ifges(6,5) + Ifges(7,5);
t472 = Ifges(6,2) + Ifges(7,2);
t471 = Ifges(6,6) + Ifges(7,6);
t470 = Ifges(6,3) + Ifges(7,3);
t266 = sin(qJ(5));
t270 = cos(qJ(5));
t296 = t266 * t267 - t270 * t271;
t446 = qJD(4) + qJD(5);
t139 = t446 * t296;
t289 = t296 * t272;
t167 = qJD(1) * t289;
t497 = t139 - t167;
t211 = t266 * t271 + t267 * t270;
t140 = t446 * t211;
t290 = t272 * t211;
t166 = qJD(1) * t290;
t496 = t140 - t166;
t163 = -qJD(3) * pkin(3) - t174;
t348 = qJD(3) * t271;
t351 = qJD(1) * t268;
t208 = -t267 * t351 + t348;
t128 = -pkin(4) * t208 + t163;
t240 = qJD(4) - t350;
t209 = qJD(3) * t267 + t271 * t351;
t175 = t268 * qJD(2) + t272 * t228;
t164 = qJD(3) * pkin(8) + t175;
t313 = pkin(3) * t272 + pkin(8) * t268;
t295 = -pkin(2) - t313;
t265 = cos(pkin(10));
t404 = pkin(1) * t265;
t198 = t295 - t404;
t165 = t198 * qJD(1);
t93 = -t164 * t267 + t271 * t165;
t80 = -pkin(9) * t209 + t93;
t73 = pkin(4) * t240 + t80;
t94 = t164 * t271 + t165 * t267;
t81 = pkin(9) * t208 + t94;
t77 = t266 * t81;
t28 = t270 * t73 - t77;
t136 = t208 * t266 + t209 * t270;
t488 = qJ(6) * t136;
t23 = t28 - t488;
t236 = qJD(5) + t240;
t21 = pkin(5) * t236 + t23;
t314 = t270 * t208 - t209 * t266;
t74 = -pkin(5) * t314 + qJD(6) + t128;
t440 = t74 * mrSges(7,2) - mrSges(6,3) * t28 - mrSges(7,3) * t21;
t495 = t128 * mrSges(6,2) + t440;
t79 = t270 * t81;
t29 = t266 * t73 + t79;
t457 = qJ(6) * t314;
t24 = t29 + t457;
t441 = -t74 * mrSges(7,1) + mrSges(6,3) * t29 + mrSges(7,3) * t24;
t485 = t474 * t136;
t462 = t236 * t471 + t314 * t472 + t485;
t494 = t441 + t462 / 0.2e1 - t128 * mrSges(6,1);
t261 = -qJ(6) + t274;
t493 = -m(5) * pkin(8) + m(6) * t274 + m(7) * t261 - t499;
t341 = qJD(1) * qJD(3);
t221 = qJDD(1) * t272 - t268 * t341;
t203 = qJDD(4) - t221;
t196 = qJDD(5) + t203;
t222 = qJDD(1) * t268 + t272 * t341;
t129 = qJD(4) * t208 + qJDD(3) * t267 + t222 * t271;
t130 = -qJD(4) * t209 + qJDD(3) * t271 - t222 * t267;
t45 = qJD(5) * t314 + t129 * t270 + t130 * t266;
t46 = -qJD(5) * t136 - t129 * t266 + t130 * t270;
t492 = -t472 * t46 / 0.2e1 - t474 * t45 / 0.2e1 - t471 * t196 / 0.2e1;
t349 = qJD(3) * t268;
t486 = qJD(2) * qJD(3) + t245 * qJDD(1);
t115 = t268 * qJDD(2) - t228 * t349 + t272 * t486;
t107 = qJDD(3) * pkin(8) + t115;
t246 = -pkin(2) - t404;
t227 = t246 * qJDD(1);
t126 = -pkin(3) * t221 - pkin(8) * t222 + t227;
t344 = qJD(4) * t271;
t346 = qJD(4) * t267;
t30 = t271 * t107 + t267 * t126 - t164 * t346 + t165 * t344;
t491 = t30 * mrSges(5,2);
t31 = -qJD(4) * t94 - t107 * t267 + t271 * t126;
t490 = t31 * mrSges(5,1);
t234 = t274 * t267;
t235 = t274 * t271;
t145 = t266 * t234 - t270 * t235;
t464 = -qJD(5) * t145 - t266 * t501 + t270 * t500;
t342 = qJD(5) * t270;
t343 = qJD(5) * t266;
t463 = t234 * t342 + t235 * t343 + t266 * t500 + t270 * t501;
t487 = t474 * t314;
t347 = qJD(3) * t272;
t284 = t267 * t347 + t268 * t344;
t455 = -t175 + (-t328 + t346) * pkin(4);
t262 = qJ(1) + pkin(10);
t254 = sin(t262);
t255 = cos(t262);
t484 = g(1) * t255 + g(2) * t254;
t410 = t236 / 0.2e1;
t424 = t136 / 0.2e1;
t427 = t314 / 0.2e1;
t482 = t410 * t471 + t424 * t474 + t427 * t472;
t481 = t410 * t473 + t424 * t475 + t427 * t474;
t384 = Ifges(4,4) * t268;
t304 = t272 * Ifges(4,2) + t384;
t480 = t24 * mrSges(7,2) + t29 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t304 / 0.2e1 - t21 * mrSges(7,1) - t28 * mrSges(6,1);
t438 = m(6) * pkin(4);
t431 = t129 / 0.2e1;
t430 = t130 / 0.2e1;
t416 = t203 / 0.2e1;
t428 = -t314 / 0.2e1;
t478 = -qJD(1) / 0.2e1;
t468 = t196 * t473 + t45 * t475 + t46 * t474;
t467 = -pkin(5) * t351 + qJ(6) * t497 - qJD(6) * t211 + t464;
t37 = -mrSges(7,2) * t196 + mrSges(7,3) * t46;
t38 = -mrSges(6,2) * t196 + mrSges(6,3) * t46;
t466 = t37 + t38;
t465 = -qJ(6) * t496 - qJD(6) * t296 + t463;
t461 = t136 * t475 + t236 * t473 + t487;
t460 = t438 + mrSges(5,1);
t72 = -mrSges(5,1) * t130 + mrSges(5,2) * t129;
t459 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t222 + t72;
t458 = pkin(5) * t496 + t455;
t177 = t296 * t268;
t179 = t271 * t198;
t359 = t268 * t271;
t109 = -pkin(9) * t359 + t179 + (-t245 * t267 - pkin(4)) * t272;
t215 = t245 * t357;
t138 = t267 * t198 + t215;
t361 = t267 * t268;
t123 = -pkin(9) * t361 + t138;
t59 = t266 * t109 + t270 * t123;
t333 = mrSges(4,3) * t351;
t456 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t208 + mrSges(5,2) * t209 + t333;
t454 = t196 * t470 + t45 * t473 + t46 * t471;
t363 = t258 * t272;
t160 = t254 * t259 - t255 * t363;
t362 = t259 * t272;
t161 = t254 * t258 + t255 * t362;
t453 = t160 * t477 - t476 * t161;
t158 = t254 * t363 + t255 * t259;
t159 = -t254 * t362 + t255 * t258;
t452 = -t158 * t477 + t476 * t159;
t199 = Ifges(5,4) * t208;
t119 = t209 * Ifges(5,1) + t240 * Ifges(5,5) + t199;
t253 = Ifges(4,4) * t350;
t451 = Ifges(4,1) * t351 + Ifges(4,5) * qJD(3) + t271 * t119 + t253;
t116 = qJDD(2) * t272 - t228 * t347 - t268 * t486;
t450 = t115 * t272 - t116 * t268;
t90 = mrSges(5,1) * t203 - mrSges(5,3) * t129;
t91 = -mrSges(5,2) * t203 + mrSges(5,3) * t130;
t449 = -t267 * t90 + t271 * t91;
t448 = -t267 * t31 + t271 * t30;
t447 = mrSges(6,1) * t258 - t259 * t476;
t445 = t209 * Ifges(5,5) + t208 * Ifges(5,6) + t240 * Ifges(5,3) + t136 * t473 + t236 * t470 + t314 * t471;
t444 = -m(7) - m(6) - m(5) - m(4);
t400 = pkin(4) * t267;
t223 = pkin(5) * t258 + t400;
t442 = -m(7) * t223 + mrSges(3,2) - mrSges(4,3);
t311 = t272 * mrSges(4,1) - mrSges(4,2) * t268;
t439 = t268 * t499 + mrSges(3,1) + t311;
t437 = m(7) * pkin(5);
t436 = t45 / 0.2e1;
t435 = t46 / 0.2e1;
t434 = Ifges(5,1) * t431 + Ifges(5,4) * t430 + Ifges(5,5) * t416;
t425 = -t136 / 0.2e1;
t417 = t196 / 0.2e1;
t414 = t209 / 0.2e1;
t411 = -t236 / 0.2e1;
t269 = sin(qJ(1));
t403 = pkin(1) * t269;
t402 = pkin(4) * t209;
t399 = pkin(4) * t270;
t394 = g(3) * t268;
t273 = cos(qJ(1));
t260 = t273 * pkin(1);
t33 = t270 * t80 - t77;
t390 = mrSges(5,3) * t208;
t389 = mrSges(5,3) * t209;
t388 = mrSges(6,3) * t314;
t387 = mrSges(6,3) * t136;
t386 = mrSges(7,3) * t314;
t385 = mrSges(7,3) * t136;
t383 = Ifges(4,4) * t272;
t382 = Ifges(5,4) * t267;
t381 = Ifges(5,4) * t271;
t378 = t209 * Ifges(5,4);
t366 = t223 * t272;
t365 = t254 * t267;
t364 = t255 * t267;
t360 = t267 * t272;
t233 = t268 * t245;
t220 = t312 * qJD(3);
t327 = t245 * t349;
t352 = t271 * t220 + t267 * t327;
t182 = pkin(4) * t361 + t233;
t345 = qJD(4) * t268;
t336 = pkin(4) * t343;
t335 = pkin(4) * t342;
t332 = mrSges(4,3) * t350;
t331 = Ifges(5,5) * t129 + Ifges(5,6) * t130 + Ifges(5,3) * t203;
t225 = t245 * t347;
t142 = pkin(4) * t284 + t225;
t118 = t208 * Ifges(5,2) + t240 * Ifges(5,6) + t378;
t324 = -t267 * t118 / 0.2e1;
t15 = -t46 * mrSges(7,1) + t45 * mrSges(7,2);
t32 = -t266 * t80 - t79;
t58 = t270 * t109 - t123 * t266;
t144 = t270 * t234 + t235 * t266;
t310 = mrSges(4,1) * t268 + mrSges(4,2) * t272;
t308 = mrSges(5,1) * t267 + mrSges(5,2) * t271;
t306 = Ifges(5,1) * t271 - t382;
t305 = Ifges(5,1) * t267 + t381;
t303 = -Ifges(5,2) * t267 + t381;
t302 = Ifges(5,2) * t271 + t382;
t301 = Ifges(4,5) * t272 - Ifges(4,6) * t268;
t300 = Ifges(5,5) * t271 - Ifges(5,6) * t267;
t299 = Ifges(5,5) * t267 + Ifges(5,6) * t271;
t298 = t216 * t272 - t261 * t268;
t297 = t252 * t272 - t268 * t274;
t171 = t254 * t271 - t255 * t360;
t169 = t254 * t360 + t255 * t271;
t20 = pkin(4) * t203 - pkin(9) * t129 + t31;
t25 = pkin(9) * t130 + t30;
t5 = t266 * t20 + t270 * t25 + t73 * t342 - t343 * t81;
t60 = t294 * qJD(3) + (-t215 + (pkin(9) * t268 - t198) * t267) * qJD(4) + t352;
t82 = t198 * t344 + t267 * t220 + (-t268 * t348 - t272 * t346) * t245;
t69 = -pkin(9) * t284 + t82;
t9 = t109 * t342 - t123 * t343 + t266 * t60 + t270 * t69;
t293 = t163 * t308;
t292 = t246 * qJD(1) * t310;
t291 = t268 * (Ifges(4,1) * t272 - t384);
t285 = -t267 * t345 + t271 * t347;
t108 = -qJDD(3) * pkin(3) - t116;
t283 = Ifges(5,5) * t268 + t272 * t306;
t282 = Ifges(5,6) * t268 + t272 * t303;
t281 = Ifges(5,3) * t268 + t272 * t300;
t6 = -qJD(5) * t29 + t270 * t20 - t25 * t266;
t10 = -qJD(5) * t59 - t266 * t69 + t270 * t60;
t68 = -pkin(4) * t130 + t108;
t2 = pkin(5) * t196 - qJ(6) * t45 - qJD(6) * t136 + t6;
t3 = qJ(6) * t46 + qJD(6) * t314 + t5;
t280 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2) + t454;
t279 = (-t267 * t94 - t271 * t93) * qJD(4) + t448;
t251 = pkin(5) + t399;
t231 = -qJD(3) * mrSges(4,2) + t332;
t197 = t308 * t268;
t180 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t221;
t176 = t211 * t268;
t172 = t255 * t357 + t365;
t170 = -t254 * t357 + t364;
t168 = pkin(5) * t296 - t252;
t157 = mrSges(5,1) * t240 - t389;
t156 = -mrSges(5,2) * t240 + t390;
t137 = -t245 * t360 + t179;
t127 = pkin(5) * t176 + t182;
t112 = -qJ(6) * t296 + t145;
t111 = -qJ(6) * t211 + t144;
t104 = mrSges(6,1) * t236 - t387;
t103 = mrSges(7,1) * t236 - t385;
t102 = -mrSges(6,2) * t236 + t388;
t101 = -mrSges(7,2) * t236 + t386;
t95 = pkin(5) * t136 + t402;
t85 = -qJD(3) * t290 + t177 * t446;
t84 = -qJD(3) * t289 - t140 * t268;
t83 = -qJD(4) * t138 + t352;
t76 = -mrSges(6,1) * t314 + mrSges(6,2) * t136;
t75 = -mrSges(7,1) * t314 + mrSges(7,2) * t136;
t55 = -pkin(5) * t85 + t142;
t53 = t129 * Ifges(5,4) + t130 * Ifges(5,2) + t203 * Ifges(5,6);
t48 = -qJ(6) * t176 + t59;
t47 = -pkin(5) * t272 + qJ(6) * t177 + t58;
t36 = mrSges(6,1) * t196 - mrSges(6,3) * t45;
t35 = mrSges(7,1) * t196 - mrSges(7,3) * t45;
t27 = t33 - t488;
t26 = t32 - t457;
t17 = -pkin(5) * t46 + qJDD(6) + t68;
t16 = -mrSges(6,1) * t46 + mrSges(6,2) * t45;
t8 = qJ(6) * t85 - qJD(6) * t176 + t9;
t7 = pkin(5) * t349 - qJ(6) * t84 + qJD(6) * t177 + t10;
t1 = [(t384 + t304) * t221 / 0.2e1 - t227 * t311 + t272 * t491 + t176 * t492 + (-t174 * t347 + t450) * mrSges(4,3) + (t272 * (-Ifges(4,2) * t268 + t383) + t291) * t341 / 0.2e1 - (t271 * t118 + t267 * t119) * t345 / 0.2e1 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t265 - 0.2e1 * mrSges(3,2) * t264 + m(3) * (t264 ^ 2 + t265 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t222 * t268 * Ifges(4,1) + t459 * t233 + t5 * (mrSges(6,2) * t272 - mrSges(6,3) * t176) + t3 * (mrSges(7,2) * t272 - mrSges(7,3) * t176) + t272 * (Ifges(4,4) * t222 + Ifges(4,2) * t221) / 0.2e1 - t272 * t490 + (-t284 * t94 - t285 * t93 - t30 * t361 - t31 * t359) * mrSges(5,3) + qJD(3) ^ 2 * t301 / 0.2e1 + m(4) * (t227 * t246 + ((-t174 * t272 - t175 * t268) * qJD(3) + t450) * t245) + t451 * t347 / 0.2e1 - (t331 + t454) * t272 / 0.2e1 + t456 * t225 + t17 * (mrSges(7,1) * t176 - mrSges(7,2) * t177) + t68 * (mrSges(6,1) * t176 - mrSges(6,2) * t177) + t6 * (-mrSges(6,1) * t272 + mrSges(6,3) * t177) + t2 * (-mrSges(7,1) * t272 + mrSges(7,3) * t177) + (-Ifges(5,6) * t272 + t268 * t303) * t430 + (-Ifges(5,5) * t272 + t268 * t306) * t431 + t359 * t434 + (qJD(3) * t283 - t305 * t345) * t414 + (-Ifges(5,3) * t272 + t268 * t300) * t416 + t222 * t383 / 0.2e1 + (t445 / 0.2e1 - t175 * mrSges(4,3) - t94 * mrSges(5,2) + t93 * mrSges(5,1) + t470 * t410 + t471 * t427 + t473 * t424 - t480) * t349 + t324 * t347 + qJD(3) * t292 + m(6) * (t10 * t28 + t128 * t142 + t182 * t68 + t29 * t9 + t5 * t59 + t58 * t6) + m(7) * (t127 * t17 + t2 * t47 + t21 * t7 + t24 * t8 + t3 * t48 + t55 * t74) + t163 * (mrSges(5,1) * t284 + mrSges(5,2) * t285) - t53 * t361 / 0.2e1 + m(5) * (t137 * t31 + t138 * t30 + t82 * t94 + t83 * t93 + (t108 * t268 + t163 * t347) * t245) + t208 * (qJD(3) * t282 - t302 * t345) / 0.2e1 + t240 * (qJD(3) * t281 - t299 * t345) / 0.2e1 - t231 * t327 + (t482 + t494) * t85 + (t461 / 0.2e1 + t481 + t495) * t84 + t47 * t35 + t48 * t37 - t468 * t177 / 0.2e1 + t272 * t245 * t180 + t58 * t36 + t59 * t38 + (-t176 * t471 - t177 * t473 - t272 * t470) * t417 + (-t176 * t472 - t177 * t474 - t272 * t471) * t435 + (-t176 * t474 - t177 * t475 - t272 * t473) * t436 + t55 * t75 + (-t365 * t438 - m(3) * t260 - mrSges(2,1) * t273 - t172 * mrSges(5,1) + mrSges(2,2) * t269 - t171 * mrSges(5,2) + t444 * (t255 * pkin(2) + t254 * pkin(7) + t260) + t442 * t254 + t477 * t161 + t476 * t160 + (-m(5) * t313 - m(6) * t297 - m(7) * t298 - t439) * t255) * g(2) + (-t364 * t438 + m(3) * t403 + mrSges(2,1) * t269 - t170 * mrSges(5,1) + mrSges(2,2) * t273 - t169 * mrSges(5,2) + t444 * (t255 * pkin(7) - t403) + t442 * t255 + t477 * t159 + t476 * t158 + (-m(7) * (-pkin(2) - t298) - m(6) * (-pkin(2) - t297) - m(5) * t295 + m(4) * pkin(2) + t439) * t254) * g(1) + t8 * t101 + t9 * t102 + t7 * t103 + t10 * t104 + t127 * t15 + t137 * t90 + t138 * t91 + t142 * t76 + t82 * t156 + t83 * t157 + t182 * t16 + t108 * t197 + t246 * (-mrSges(4,1) * t221 + mrSges(4,2) * t222) + qJDD(3) * (Ifges(4,5) * t268 + Ifges(4,6) * t272); m(3) * qJDD(2) + (t103 + t104) * t85 + (t101 + t102) * t84 - t466 * t177 - (t35 + t36) * t176 + (-m(3) + t444) * g(3) + (-t15 - t16 + (t156 * t271 - t157 * t267 + t231) * qJD(3) - t459) * t272 + (t180 + (-t267 * t156 - t271 * t157) * qJD(4) + (t75 + t76 + t456) * qJD(3) + t449) * t268 + m(6) * (t128 * t349 - t176 * t6 - t177 * t5 - t272 * t68 + t28 * t85 + t29 * t84) + m(7) * (-t17 * t272 - t176 * t2 - t177 * t3 + t21 * t85 + t24 * t84 + t349 * t74) + m(4) * (t115 * t268 + t116 * t272 + (-t174 * t268 + t175 * t272) * qJD(3)) + m(5) * ((-t108 + (-t267 * t93 + t271 * t94) * qJD(3)) * t272 + (qJD(3) * t163 + t279) * t268); t17 * (mrSges(7,1) * t296 + mrSges(7,2) * t211) + t68 * (mrSges(6,1) * t296 + mrSges(6,2) * t211) + (t211 * t473 - t296 * t471) * t417 + (t211 * t474 - t296 * t472) * t435 + (t211 * t475 - t296 * t474) * t436 + (t166 / 0.2e1 - t140 / 0.2e1) * t462 + t296 * t492 + (t293 + t324) * qJD(4) - t445 * t351 / 0.2e1 + (-t344 * t93 - t346 * t94 + t448) * mrSges(5,3) + (m(5) * t279 - t156 * t346 - t157 * t344 + t449) * pkin(8) - (t440 + t481) * t139 + (t167 / 0.2e1 - t139 / 0.2e1) * t461 + (t166 * t24 - t167 * t21 - t2 * t211 - t296 * t3) * mrSges(7,3) + (t166 * t29 - t167 * t28 - t211 * t6 - t296 * t5) * mrSges(6,3) + (-t166 * t471 - t167 * t473) * t411 + (-t166 * t472 - t167 * t474) * t428 + (-t166 * t474 - t167 * t475) * t425 - t74 * (mrSges(7,1) * t166 - mrSges(7,2) * t167) + (t208 * t303 + t209 * t306 + t240 * t300) * qJD(4) / 0.2e1 + (-pkin(3) * t108 - t121 * t93 - t122 * t94) * m(5) + (t332 - t231) * t174 - (t441 + t482) * t140 + t458 * t75 + t108 * t309 - (-Ifges(4,2) * t351 + t253 + t451) * t350 / 0.2e1 + t455 * t76 + (-m(5) * t163 + t333 - t456) * t175 + t484 * (t268 * t498 + t310) + (-g(3) * t498 + t484 * t493) * t272 + t302 * t430 + t305 * t431 + t267 * t434 + t299 * t416 + (t470 * t411 + t473 * t425 + t471 * t428 + t480) * t351 + (t208 * t282 + t209 * t283 + t240 * t281) * t478 - t293 * t350 + t119 * t344 / 0.2e1 - t301 * t341 / 0.2e1 + t118 * t328 / 0.2e1 + (mrSges(6,1) * t496 - mrSges(6,2) * t497) * t128 + (t493 * t268 - t311) * g(3) + Ifges(4,3) * qJDD(3) + t463 * t102 + t464 * t104 + (t128 * t455 + t144 * t6 + t145 * t5 - t252 * t68 + t28 * t464 + t29 * t463) * m(6) + t465 * t101 + t467 * t103 + (t111 * t2 + t112 * t3 + t168 * t17 + t21 * t467 + t24 * t465 + t458 * t74) * m(7) + t468 * t211 / 0.2e1 + (-t292 - t93 * (mrSges(5,1) * t268 - mrSges(5,3) * t357) - t94 * (-mrSges(5,2) * t268 - mrSges(5,3) * t360) + t291 * t478) * qJD(1) - pkin(3) * t72 + t111 * t35 + t112 * t37 - t115 * mrSges(4,2) + t116 * mrSges(4,1) + t144 * t36 + t145 * t38 - t122 * t156 - t121 * t157 + t168 * t15 + Ifges(4,6) * t221 + Ifges(4,5) * t222 - t252 * t16 + t271 * t53 / 0.2e1; (t389 + t157) * t94 + (t390 - t156) * t93 + t490 - t491 + (-t336 - t32) * t104 + (-t336 - t26) * t103 + (t335 - t33) * t102 + (t335 - t27) * t101 + (m(6) * t400 + mrSges(7,1) * t258 + t447) * t394 - (-Ifges(5,2) * t209 + t119 + t199) * t208 / 0.2e1 + (-m(7) * (-t224 * t255 - t254 * t366) - mrSges(5,2) * t170 + t460 * t169 + t452) * g(2) + (-m(7) * (t224 * t254 - t255 * t366) + mrSges(5,2) * t172 - t460 * t171 + t453) * g(1) + (t2 * t251 + (t266 * t3 + (-t21 * t266 + t24 * t270) * qJD(5)) * pkin(4) + t223 * t394 - t21 * t26 - t24 * t27 - t74 * t95) * m(7) - t209 * (Ifges(5,1) * t208 - t378) / 0.2e1 + (t266 * t5 + t270 * t6 + (-t266 * t28 + t270 * t29) * qJD(5)) * t438 + t118 * t414 - t76 * t402 - m(6) * (t128 * t402 + t28 * t32 + t29 * t33) + t331 + t280 + t461 * t428 + t36 * t399 + (-t471 * t411 - t474 * t425 - t472 * t428 + t494) * t136 + (t411 * t473 + t425 * t475 + t428 * t474 - t495) * t314 + t466 * pkin(4) * t266 - t95 * t75 + g(3) * t197 - t163 * (mrSges(5,1) * t209 + mrSges(5,2) * t208) - t240 * (Ifges(5,5) * t208 - Ifges(5,6) * t209) / 0.2e1 + t251 * t35; t2 * t437 + t21 * t386 + t280 - t23 * t101 - t74 * (mrSges(7,1) * t136 + mrSges(7,2) * t314) - t128 * (mrSges(6,1) * t136 + mrSges(6,2) * t314) + (t314 * t475 - t485) * t425 + t462 * t424 + (-t136 * t471 + t314 * t473) * t411 + (-(-mrSges(7,1) - t437) * t258 + t447) * t394 + (t387 + t104) * t29 + (t388 - t102) * t28 + (-m(7) * (-t21 + t23) + t385 + t103) * t24 + (t158 * t437 + t452) * g(2) + (-t160 * t437 + t453) * g(1) + (-t136 * t472 + t461 + t487) * t428 + (t35 + (-m(7) * t74 - t75) * t136) * pkin(5); -t314 * t101 + t136 * t103 + (g(3) * t272 + t21 * t136 - t24 * t314 - t268 * t484 + t17) * m(7) + t15;];
tau  = t1;
