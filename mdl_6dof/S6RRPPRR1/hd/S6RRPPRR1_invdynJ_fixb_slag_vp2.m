% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:53
% EndTime: 2019-03-09 08:46:33
% DurationCPUTime: 25.80s
% Computational Cost: add. (11282->713), mult. (25362->914), div. (0->0), fcn. (18424->12), ass. (0->325)
t470 = -m(7) - m(6);
t469 = mrSges(4,1) + mrSges(5,1);
t468 = -mrSges(5,2) - mrSges(4,3);
t458 = Ifges(4,1) + Ifges(5,1);
t457 = Ifges(5,4) + Ifges(4,5);
t247 = sin(qJ(6));
t251 = cos(qJ(6));
t209 = -mrSges(7,1) * t251 + mrSges(7,2) * t247;
t467 = -m(7) * pkin(5) - mrSges(6,1) + t209;
t239 = -qJD(2) + qJD(5);
t248 = sin(qJ(5));
t252 = cos(qJ(5));
t246 = -qJ(3) - pkin(7);
t249 = sin(qJ(2));
t210 = t246 * t249;
t192 = qJD(1) * t210;
t182 = qJD(2) * pkin(2) + t192;
t253 = cos(qJ(2));
t212 = t246 * t253;
t193 = qJD(1) * t212;
t244 = sin(pkin(10));
t343 = t244 * t193;
t354 = cos(pkin(10));
t116 = t182 * t354 + t343;
t272 = qJD(4) - t116;
t188 = t244 * t253 + t249 * t354;
t171 = t188 * qJD(1);
t386 = t171 * pkin(8);
t405 = -pkin(3) - pkin(4);
t73 = qJD(2) * t405 + t272 - t386;
t309 = t354 * t193;
t117 = t244 * t182 - t309;
t105 = qJD(2) * qJ(4) + t117;
t308 = t354 * t253;
t340 = qJD(1) * t249;
t169 = -qJD(1) * t308 + t244 * t340;
t389 = pkin(8) * t169;
t77 = t105 + t389;
t38 = -t248 * t77 + t252 * t73;
t36 = -pkin(5) * t239 - t38;
t421 = t169 * t248 + t252 * t171;
t376 = mrSges(6,3) * t421;
t81 = t239 * t251 - t247 * t421;
t82 = t239 * t247 + t251 * t421;
t453 = mrSges(6,1) * t239 + mrSges(7,1) * t81 - mrSges(7,2) * t82 - t376;
t466 = -m(7) * t36 + t453;
t456 = Ifges(5,5) - Ifges(4,4);
t315 = t354 * pkin(2);
t230 = -t315 - pkin(3);
t224 = -pkin(4) + t230;
t392 = pkin(2) * t244;
t226 = qJ(4) + t392;
t147 = t248 * t224 + t252 * t226;
t124 = t192 * t244 - t309;
t274 = t124 + t389;
t125 = t354 * t192 + t343;
t83 = t125 + t386;
t447 = qJD(5) * t147 + t252 * t274 + (qJD(4) - t83) * t248;
t211 = -mrSges(3,1) * t253 + mrSges(3,2) * t249;
t240 = qJ(2) + pkin(10);
t235 = sin(t240);
t236 = cos(t240);
t465 = t211 - t469 * t236 - (-mrSges(4,2) + mrSges(5,3)) * t235;
t331 = qJD(1) * qJD(2);
t313 = t249 * t331;
t330 = qJDD(1) * t253;
t198 = -t313 + t330;
t199 = qJDD(1) * t249 + t253 * t331;
t128 = -t198 * t354 + t199 * t244;
t129 = t244 * t198 + t199 * t354;
t353 = qJDD(1) * pkin(1);
t269 = pkin(2) * t198 - qJDD(3) + t353;
t257 = qJ(4) * t129 + qJD(4) * t171 + t269;
t35 = t128 * t405 + t257;
t102 = t169 * t252 - t248 * t171;
t49 = qJD(5) * t102 + t128 * t248 + t129 * t252;
t50 = -qJD(5) * t421 + t128 * t252 - t129 * t248;
t11 = -pkin(5) * t50 - pkin(9) * t49 + t35;
t339 = qJD(1) * t253;
t204 = -qJD(1) * pkin(1) - pkin(2) * t339 + qJD(3);
t86 = t169 * pkin(3) - t171 * qJ(4) + t204;
t71 = -pkin(4) * t169 - t86;
t33 = -pkin(5) * t102 - pkin(9) * t421 + t71;
t39 = t248 * t73 + t252 * t77;
t37 = pkin(9) * t239 + t39;
t15 = -t247 * t37 + t251 * t33;
t238 = -qJDD(2) + qJDD(5);
t186 = t199 * pkin(7);
t336 = qJD(3) * t249;
t112 = qJDD(2) * pkin(2) - qJ(3) * t199 - qJD(1) * t336 - t186;
t232 = pkin(7) * t330;
t338 = qJD(2) * t249;
t325 = pkin(7) * t338;
t335 = qJD(3) * t253;
t120 = qJ(3) * t198 + t232 + (-t325 + t335) * qJD(1);
t65 = t112 * t354 - t244 * t120;
t267 = qJDD(4) - t65;
t45 = -t129 * pkin(8) + qJDD(2) * t405 + t267;
t66 = t244 * t112 + t354 * t120;
t61 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t66;
t46 = pkin(8) * t128 + t61;
t9 = qJD(5) * t38 + t248 * t45 + t252 * t46;
t7 = pkin(9) * t238 + t9;
t1 = qJD(6) * t15 + t11 * t247 + t251 * t7;
t16 = t247 * t33 + t251 * t37;
t2 = -qJD(6) * t16 + t11 * t251 - t247 * t7;
t256 = t1 * t251 - t2 * t247 + (-t15 * t251 - t16 * t247) * qJD(6);
t333 = qJD(6) * t251;
t334 = qJD(6) * t247;
t28 = qJD(6) * t81 + t238 * t247 + t251 * t49;
t48 = qJDD(6) - t50;
t13 = mrSges(7,1) * t48 - mrSges(7,3) * t28;
t29 = -qJD(6) * t82 + t238 * t251 - t247 * t49;
t14 = -mrSges(7,2) * t48 + mrSges(7,3) * t29;
t417 = -t247 * t13 + t251 * t14;
t94 = qJD(6) - t102;
t55 = -mrSges(7,2) * t94 + mrSges(7,3) * t81;
t56 = mrSges(7,1) * t94 - mrSges(7,3) * t82;
t464 = m(7) * t256 - t56 * t333 - t55 * t334 + t417;
t451 = t239 * Ifges(6,5);
t93 = Ifges(6,4) * t102;
t60 = Ifges(6,1) * t421 + t451 + t93;
t463 = -t38 * mrSges(6,3) + t60 / 0.2e1;
t30 = t82 * Ifges(7,5) + t81 * Ifges(7,6) + t94 * Ifges(7,3);
t370 = Ifges(6,4) * t421;
t450 = t239 * Ifges(6,6);
t452 = t102 * Ifges(6,2);
t59 = t370 + t450 + t452;
t462 = -t39 * mrSges(6,3) - t59 / 0.2e1 + t30 / 0.2e1;
t395 = t188 / 0.2e1;
t460 = t198 / 0.2e1;
t459 = t71 * mrSges(6,2);
t381 = qJD(2) / 0.2e1;
t455 = Ifges(5,6) - Ifges(4,6);
t12 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t43 = mrSges(6,1) * t238 - mrSges(6,3) * t49;
t454 = t12 - t43;
t160 = Ifges(4,4) * t169;
t367 = Ifges(5,5) * t169;
t448 = qJD(2) * t457 + t458 * t171 - t160 + t367;
t285 = t247 * t56 - t251 * t55;
t377 = mrSges(6,3) * t102;
t84 = -mrSges(6,2) * t239 + t377;
t277 = t285 - t84;
t446 = qJD(2) * t469 + t171 * t468;
t378 = mrSges(5,2) * t169;
t145 = qJD(2) * mrSges(5,3) - t378;
t445 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t169 + t145;
t164 = t235 * t248 + t236 * t252;
t346 = t235 * t252;
t165 = -t236 * t248 + t346;
t444 = -t165 * mrSges(6,2) + t164 * t467;
t254 = cos(qJ(1));
t341 = t252 * t254;
t342 = t248 * t254;
t138 = -t235 * t342 - t236 * t341;
t139 = -t235 * t341 + t236 * t342;
t443 = t138 * mrSges(6,2) + t139 * t467;
t250 = sin(qJ(1));
t345 = t236 * t250;
t136 = t248 * t345 - t250 * t346;
t137 = t164 * t250;
t442 = -t137 * mrSges(6,2) + t136 * t467;
t265 = -t244 * t249 + t308;
t119 = t188 * t252 - t248 * t265;
t170 = t188 * qJD(2);
t172 = t265 * qJD(2);
t278 = -t188 * t248 - t252 * t265;
t67 = qJD(5) * t278 + t170 * t248 + t172 * t252;
t271 = t119 * t333 + t247 * t67;
t185 = -pkin(7) * t313 + t232;
t440 = t185 * t253 + t186 * t249;
t439 = g(1) * t254 + g(2) * t250;
t322 = m(7) * pkin(9) + mrSges(7,3);
t437 = -mrSges(6,2) + t322;
t436 = -m(3) * pkin(1) - mrSges(2,1) + t465;
t435 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t63 = pkin(5) * t421 - pkin(9) * t102;
t433 = t204 * mrSges(4,2) - t116 * mrSges(4,3) - t86 * mrSges(5,3);
t432 = -t204 * mrSges(4,1) - t86 * mrSges(5,1) + t105 * mrSges(5,2) + t117 * mrSges(4,3);
t431 = t71 * mrSges(6,1) + t15 * mrSges(7,1) - t16 * mrSges(7,2);
t430 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) + mrSges(6,3) + t470 * (-pkin(8) - t246) + t468;
t394 = -t239 / 0.2e1;
t404 = -t102 / 0.2e1;
t407 = -t94 / 0.2e1;
t409 = -t82 / 0.2e1;
t411 = -t81 / 0.2e1;
t429 = Ifges(7,5) * t409 - Ifges(6,2) * t404 - Ifges(6,6) * t394 + Ifges(7,6) * t411 + Ifges(7,3) * t407 - t431;
t292 = t247 * mrSges(7,1) + t251 * mrSges(7,2);
t268 = t36 * t292;
t287 = Ifges(7,5) * t251 - Ifges(7,6) * t247;
t368 = Ifges(7,4) * t251;
t289 = -Ifges(7,2) * t247 + t368;
t369 = Ifges(7,4) * t247;
t291 = Ifges(7,1) * t251 - t369;
t383 = t82 * Ifges(7,4);
t31 = t81 * Ifges(7,2) + t94 * Ifges(7,6) + t383;
t79 = Ifges(7,4) * t81;
t32 = t82 * Ifges(7,1) + t94 * Ifges(7,5) + t79;
t358 = t251 * t32;
t374 = mrSges(7,3) * t251;
t375 = mrSges(7,3) * t247;
t393 = t247 / 0.2e1;
t403 = -t421 / 0.2e1;
t428 = Ifges(6,1) * t403 + Ifges(6,5) * t394 + t15 * t374 + t16 * t375 + t287 * t407 + t289 * t411 + t291 * t409 - t268 - t358 / 0.2e1 + t31 * t393 - t459;
t216 = qJ(4) * t345;
t391 = pkin(2) * t249;
t276 = t235 * t405 - t391;
t420 = t250 * t276 + t216;
t225 = t235 * qJ(4);
t344 = t236 * t254;
t419 = pkin(3) * t344 + t254 * t225;
t218 = qJ(4) * t344;
t418 = t254 * t276 + t218;
t10 = -qJD(5) * t39 - t248 * t46 + t252 * t45;
t415 = t28 / 0.2e1;
t414 = t29 / 0.2e1;
t412 = t48 / 0.2e1;
t410 = t81 / 0.2e1;
t408 = t82 / 0.2e1;
t406 = t94 / 0.2e1;
t402 = t421 / 0.2e1;
t400 = -t169 / 0.2e1;
t399 = t169 / 0.2e1;
t397 = t171 / 0.2e1;
t390 = pkin(7) * t253;
t228 = t236 * pkin(3);
t237 = t253 * pkin(2);
t382 = -qJD(2) / 0.2e1;
t231 = t237 + pkin(1);
t373 = Ifges(3,4) * t249;
t372 = Ifges(3,4) * t253;
t371 = Ifges(4,4) * t171;
t352 = t119 * t247;
t351 = t119 * t251;
t310 = qJD(2) * t246;
t166 = t249 * t310 + t335;
t167 = t253 * t310 - t336;
t90 = t354 * t166 + t244 * t167;
t135 = t244 * t210 - t354 * t212;
t337 = qJD(2) * t252;
t332 = m(5) - t470;
t329 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t48;
t328 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t340) * t390;
t327 = pkin(2) * t340;
t326 = pkin(2) * t338;
t317 = t358 / 0.2e1;
t316 = t228 + t225 + t237;
t311 = -t334 / 0.2e1;
t306 = t128 * mrSges(4,1) + t129 * mrSges(4,2);
t305 = t128 * mrSges(5,1) - t129 * mrSges(5,3);
t304 = -t231 - t225;
t89 = t166 * t244 - t354 * t167;
t134 = -t354 * t210 - t212 * t244;
t221 = t254 * t231;
t303 = -t246 * t250 + t221;
t115 = -pkin(3) * t265 - t188 * qJ(4) - t231;
t108 = -qJDD(2) * mrSges(5,1) + t129 * mrSges(5,2);
t302 = pkin(4) * t344 + t221 + t419;
t301 = t236 * pkin(4) + t316;
t299 = -t50 * mrSges(6,1) + t49 * mrSges(6,2);
t298 = mrSges(3,1) * t249 + mrSges(3,2) * t253;
t290 = t253 * Ifges(3,2) + t373;
t288 = Ifges(3,5) * t253 - Ifges(3,6) * t249;
t286 = t15 * t247 - t16 * t251;
t80 = pkin(4) * t265 - t115;
t41 = -pkin(5) * t278 - pkin(9) * t119 + t80;
t91 = -pkin(8) * t188 + t134;
t92 = -pkin(8) * t265 + t135;
t58 = t248 * t91 + t252 * t92;
t20 = -t247 * t58 + t251 * t41;
t21 = t247 * t41 + t251 * t58;
t57 = t248 * t92 - t252 * t91;
t281 = -t137 * t251 - t247 * t254;
t280 = t137 * t247 - t251 * t254;
t146 = t224 * t252 - t226 * t248;
t88 = t171 * pkin(3) + t169 * qJ(4) + t327;
t275 = -pkin(8) * t172 + t89;
t273 = pkin(1) * t298;
t270 = t119 * t334 - t251 * t67;
t266 = t249 * (Ifges(3,1) * t253 - t373);
t76 = t170 * pkin(3) - t172 * qJ(4) - t188 * qJD(4) + t326;
t74 = -pkin(4) * t171 - t88;
t259 = t250 * (t236 * t405 + t304);
t69 = -pkin(4) * t170 - t76;
t258 = m(5) * (-pkin(3) * t235 - t391) - t235 * mrSges(5,1) + t236 * mrSges(5,3);
t5 = t28 * Ifges(7,4) + t29 * Ifges(7,2) + t48 * Ifges(7,6);
t6 = t28 * Ifges(7,1) + t29 * Ifges(7,4) + t48 * Ifges(7,5);
t8 = -pkin(5) * t238 - t10;
t255 = Ifges(6,5) * t49 + Ifges(6,6) * t50 + Ifges(6,3) * t238 + t10 * mrSges(6,1) + t8 * t209 + t1 * t374 + t6 * t393 + t251 * t5 / 0.2e1 + (Ifges(7,1) * t247 + t368) * t415 + (Ifges(7,2) * t251 + t369) * t414 + (Ifges(7,5) * t247 + Ifges(7,6) * t251) * t412 - t9 * mrSges(6,2) + t31 * t311 - t2 * t375 + (-t15 * t333 - t16 * t334) * mrSges(7,3) + (t287 * t406 + t289 * t410 + t291 * t408 + t268 + t317) * qJD(6);
t233 = Ifges(3,4) * t339;
t208 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t339;
t179 = Ifges(3,1) * t340 + Ifges(3,5) * qJD(2) + t233;
t178 = Ifges(3,6) * qJD(2) + qJD(1) * t290;
t159 = Ifges(5,5) * t171;
t140 = pkin(5) - t146;
t133 = t171 * t247 + t251 * t337;
t132 = t171 * t251 - t247 * t337;
t127 = -t138 * t251 - t247 * t250;
t126 = t138 * t247 - t250 * t251;
t114 = mrSges(4,1) * t169 + mrSges(4,2) * t171;
t113 = mrSges(5,1) * t169 - mrSges(5,3) * t171;
t109 = -mrSges(5,2) * t128 + qJDD(2) * mrSges(5,3);
t107 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t129;
t106 = -qJDD(2) * mrSges(4,2) - mrSges(4,3) * t128;
t98 = -t169 * Ifges(4,2) + Ifges(4,6) * qJD(2) + t371;
t97 = Ifges(5,6) * qJD(2) + t169 * Ifges(5,3) + t159;
t96 = -qJD(2) * pkin(3) + t272;
t75 = pkin(8) * t170 + t90;
t68 = qJD(5) * t119 - t252 * t170 + t172 * t248;
t64 = -qJDD(2) * pkin(3) + t267;
t62 = -mrSges(6,1) * t102 + mrSges(6,2) * t421;
t54 = pkin(3) * t128 - t257;
t52 = t248 * t274 + t252 * t83;
t44 = -mrSges(6,2) * t238 + mrSges(6,3) * t50;
t34 = -t63 + t74;
t25 = t247 * t63 + t251 * t38;
t24 = -t247 * t38 + t251 * t63;
t22 = -qJD(5) * t57 + t248 * t275 + t252 * t75;
t19 = pkin(5) * t68 - pkin(9) * t67 + t69;
t18 = t247 * t34 + t251 * t52;
t17 = -t247 * t52 + t251 * t34;
t4 = -qJD(6) * t21 + t19 * t251 - t22 * t247;
t3 = qJD(6) * t20 + t19 * t247 + t22 * t251;
t23 = [(Ifges(3,1) * t199 + Ifges(3,4) * t460 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t199)) * t249 + (m(4) * t66 + m(5) * t61 + t106 + t109) * t135 + (t198 * t390 + t440) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t440) - t271 * t31 / 0.2e1 + (Ifges(6,1) * t402 + t451 / 0.2e1 + t93 / 0.2e1 + t459 + t317 + t463) * t67 + m(4) * (t204 * t326 + t231 * t269) - t269 * (-mrSges(4,1) * t265 + mrSges(4,2) * t188) + (t188 * t64 + t265 * t61) * mrSges(5,2) + t265 * (Ifges(4,4) * t129 + Ifges(4,6) * qJDD(2)) / 0.2e1 + (-t188 * t65 + t265 * t66) * mrSges(4,3) - t265 * (Ifges(5,5) * t129 + Ifges(5,6) * qJDD(2)) / 0.2e1 + t54 * (-mrSges(5,1) * t265 - mrSges(5,3) * t188) + t199 * t372 / 0.2e1 - qJDD(2) * mrSges(3,2) * t390 + (-m(4) * t65 + m(5) * t64 - t107 + t108) * t134 + (m(4) * t117 + m(5) * t105 + t445) * t90 + (-m(4) * t116 + m(5) * t96 - t446) * t89 + (t35 * mrSges(6,2) - t10 * mrSges(6,3) + Ifges(6,1) * t49 + Ifges(6,4) * t50 + Ifges(6,5) * t238 + t287 * t412 + t289 * t414 + t291 * t415 + t292 * t8 + t311 * t32) * t119 + (-m(6) * t38 - t466) * (qJD(5) * t58 + t248 * t75 - t252 * t275) + t290 * t460 + (Ifges(7,3) * t406 + Ifges(7,5) * t408 + Ifges(7,6) * t410 - Ifges(6,4) * t402 - t450 / 0.2e1 - t452 / 0.2e1 + t431 + t462) * t68 + t114 * t326 + t253 * t179 * t381 + (-Ifges(7,3) * t412 - Ifges(7,6) * t414 - Ifges(7,5) * t415 + Ifges(6,6) * t238 - t329 / 0.2e1 + Ifges(6,4) * t49 + Ifges(6,2) * t50 - t35 * mrSges(6,1) + t9 * mrSges(6,3) + t435) * t278 + (-m(6) * t259 + t137 * mrSges(6,1) - t281 * mrSges(7,1) - t280 * mrSges(7,2) - (-t137 * pkin(5) + t259) * m(7) + t437 * t136 + (-m(5) * (t304 - t228) + m(4) * t231 - t436) * t250 + ((m(4) + m(5)) * t246 + t430) * t254) * g(1) + (-m(5) * (t303 + t419) - m(7) * (-pkin(5) * t138 + t302) - t127 * mrSges(7,1) - t126 * mrSges(7,2) - m(4) * t303 - m(6) * t302 + t138 * mrSges(6,1) - t437 * t139 + t436 * t254 + t430 * t250) * g(2) + t253 * (Ifges(3,4) * t199 + Ifges(3,2) * t198) / 0.2e1 + (-Ifges(7,5) * t270 - Ifges(7,6) * t271) * t406 + (-Ifges(7,1) * t270 - Ifges(7,4) * t271) * t408 - t178 * t338 / 0.2e1 + (-t1 * t352 + t15 * t270 - t16 * t271 - t2 * t351) * mrSges(7,3) + t6 * t351 / 0.2e1 - t5 * t352 / 0.2e1 - t211 * t353 + (t253 * (-Ifges(3,2) * t249 + t372) + t266) * t331 / 0.2e1 + t253 * qJDD(2) * Ifges(3,6) - t273 * t331 + (t448 / 0.2e1 + mrSges(5,2) * t96 + Ifges(4,4) * t400 + Ifges(5,5) * t399 + t381 * t457 + t397 * t458 + t433) * t172 - t208 * t325 + Ifges(2,3) * qJDD(1) + t115 * t305 - t231 * t306 - pkin(1) * (-mrSges(3,1) * t198 + mrSges(3,2) * t199) + t76 * t113 + t22 * t84 + t69 * t62 + t3 * t55 + t4 * t56 + t58 * t44 + t20 * t13 + t21 * t14 + t80 * t299 + m(5) * (t115 * t54 + t76 * t86) + t36 * (mrSges(7,1) * t271 - mrSges(7,2) * t270) + m(7) * (t1 * t21 + t15 * t4 + t16 * t3 + t2 * t20) + m(6) * (t22 * t39 + t35 * t80 + t58 * t9 + t69 * t71) + (-m(6) * t10 + m(7) * t8 + t454) * t57 + (Ifges(5,3) * t399 - Ifges(4,2) * t400 - t98 / 0.2e1 + t97 / 0.2e1 + t456 * t397 + t455 * t381 - t432) * t170 + (t188 * t458 - t265 * t456) * t129 / 0.2e1 + (qJDD(2) * t457 + t129 * t458) * t395 + (t288 * t381 - t328) * qJD(2) + ((-Ifges(5,3) - Ifges(4,2)) * t265 + 0.2e1 * t456 * t395) * t128 + (-Ifges(7,4) * t270 - Ifges(7,2) * t271) * t410 + (0.2e1 * Ifges(3,5) * t249 + t188 * t457 - t265 * t455) * qJDD(2) / 0.2e1; (m(6) * t39 - m(7) * t286 - t277) * (qJD(4) * t252 + qJD(5) * t146) + (t140 * t8 - t15 * t17 - t16 * t18 + t447 * t36) * m(7) + (t10 * t146 + t147 * t9 - t447 * t38 - t39 * t52 - t71 * t74) * m(6) + (-Ifges(4,2) * t171 - t160 + t448) * t399 + (-t420 * m(6) + t137 * mrSges(7,3) - (-pkin(9) * t137 + t420) * m(7) - m(5) * t216 - t250 * t258 + t442) * g(2) + (-t138 * mrSges(7,3) - (t138 * pkin(9) + t418) * m(7) - m(5) * t218 - t254 * t258 - t418 * m(6) + t443) * g(1) + (-Ifges(6,4) * t404 - t428 + t463) * t102 + t464 * (-pkin(9) + t147) + (-m(7) * (-pkin(9) * t165 + t301) + t165 * mrSges(7,3) - m(4) * t237 - m(5) * t316 - m(6) * t301 + t444 + t465) * g(3) + (Ifges(3,3) + Ifges(5,2) + Ifges(4,3)) * qJDD(2) - t445 * t125 + t446 * t124 + (t116 * t124 - t117 * t125 - t204 * t327 + (t244 * t66 + t354 * t65) * pkin(2)) * m(4) + (Ifges(6,4) * t403 - t429 + t462) * t421 + t106 * t392 + t98 * t397 + (-t124 * t96 + t226 * t61 + t230 * t64 - t86 * t88 + (qJD(4) - t125) * t105) * m(5) + (t328 + (-t266 / 0.2e1 + t273) * qJD(1)) * qJD(1) + (m(4) * t391 + mrSges(4,1) * t235 + mrSges(4,2) * t236 + t298) * t439 - t255 - (-Ifges(3,2) * t340 + t179 + t233) * t339 / 0.2e1 + (t178 / 0.2e1 + pkin(7) * t208) * t340 - t288 * t331 / 0.2e1 + t96 * t378 - t114 * t327 + t230 * t108 + t226 * t109 + Ifges(3,5) * t199 + Ifges(3,6) * t198 - t185 * mrSges(3,2) - t186 * mrSges(3,1) + qJD(4) * t145 + t146 * t43 + t147 * t44 + t140 * t12 - t88 * t113 - t52 * t84 - t74 * t62 - t66 * mrSges(4,2) + t61 * mrSges(5,3) - t64 * mrSges(5,1) + t65 * mrSges(4,1) - t18 * t55 - t17 * t56 - t447 * t453 + (Ifges(5,3) * t400 + t455 * t382 + t432) * t171 + t455 * t128 + t457 * t129 + (-t457 * t382 + t433) * t169 - (-t458 * t169 + t159 - t371 + t97) * t171 / 0.2e1 - t367 * t400 + t107 * t315; -t453 * t421 - t299 + t446 * t171 + t445 * t169 + t305 - t251 * t13 - t247 * t14 + t285 * qJD(6) - t277 * t102 + t306 + (-g(1) * t250 + g(2) * t254) * (m(4) + t332) + (-t1 * t247 - t2 * t251 + t286 * t94 + t421 * t36) * m(7) + (t102 * t39 - t38 * t421 - t35) * m(6) + (t105 * t169 - t171 * t96 + t54) * m(5) + (t116 * t171 + t117 * t169 - t269) * m(4); -qJD(2) * t145 - t132 * t56 - t133 * t55 + (t113 - t62) * t171 + (-qJD(2) * t84 - qJD(5) * t277 - t454) * t252 + (t44 + (-t247 * t55 - t251 * t56) * qJD(6) - t239 * t453 + t417) * t248 + t108 + (t236 * g(3) - t235 * t439) * t332 + (-t132 * t15 - t133 * t16 + (-qJD(5) * t286 - t8) * t252 + (t239 * t36 + t256) * t248) * m(7) + (t10 * t252 - t171 * t71 + t248 * t9 + t239 * (-t248 * t38 + t252 * t39)) * m(6) + (-qJD(2) * t105 + t171 * t86 + t64) * m(5); -pkin(5) * t12 - t24 * t56 - t25 * t55 + t59 * t402 + t255 + (t93 + t60) * t404 + (t30 - t370) * t403 + (t377 - t84) * t38 + (-t165 * t322 - t444) * g(3) + (-t137 * t322 - t442) * g(2) + (t138 * t322 - t443) * g(1) + (-pkin(5) * t8 - t15 * t24 - t16 * t25) * m(7) + (t376 + t466) * t39 + t429 * t421 + t428 * t102 + t464 * pkin(9); -t36 * (mrSges(7,1) * t82 + mrSges(7,2) * t81) + (Ifges(7,1) * t81 - t383) * t409 + t31 * t408 + (Ifges(7,5) * t81 - Ifges(7,6) * t82) * t407 - t15 * t55 + t16 * t56 - g(1) * (mrSges(7,1) * t126 - mrSges(7,2) * t127) - g(2) * (-mrSges(7,1) * t280 + mrSges(7,2) * t281) + g(3) * t292 * t165 + (t15 * t81 + t16 * t82) * mrSges(7,3) + t329 + (-Ifges(7,2) * t82 + t32 + t79) * t411 - t435;];
tau  = t23;
