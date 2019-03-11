% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:51:14
% EndTime: 2019-03-09 03:51:57
% DurationCPUTime: 29.17s
% Computational Cost: add. (18609->790), mult. (45071->1026), div. (0->0), fcn. (35794->18), ass. (0->332)
t293 = sin(pkin(10));
t300 = sin(qJ(3));
t295 = cos(pkin(10));
t403 = cos(qJ(3));
t346 = t403 * t295;
t315 = -t300 * t293 + t346;
t244 = t315 * qJD(1);
t292 = sin(pkin(11));
t294 = cos(pkin(11));
t299 = sin(qJ(5));
t303 = cos(qJ(5));
t317 = t292 * t299 - t294 * t303;
t173 = t317 * t244;
t246 = t317 * qJD(5);
t492 = t246 - t173;
t257 = t292 * t303 + t294 * t299;
t172 = t257 * t244;
t247 = t257 * qJD(5);
t450 = -t247 + t172;
t240 = qJD(5) - t244;
t258 = t293 * t403 + t300 * t295;
t245 = t258 * qJD(1);
t214 = qJD(3) * t292 + t245 * t294;
t332 = t294 * qJD(3) - t245 * t292;
t153 = t214 * t303 + t299 * t332;
t384 = Ifges(6,4) * t153;
t476 = -t214 * t299 + t303 * t332;
t79 = Ifges(6,2) * t476 + t240 * Ifges(6,6) + t384;
t430 = t79 / 0.2e1;
t150 = Ifges(6,4) * t476;
t80 = t153 * Ifges(6,1) + t240 * Ifges(6,5) + t150;
t429 = t80 / 0.2e1;
t491 = mrSges(6,3) + mrSges(7,3);
t193 = pkin(3) * t245 - qJ(4) * t244;
t391 = pkin(7) + qJ(2);
t268 = t391 * t295;
t260 = qJD(1) * t268;
t242 = t300 * t260;
t266 = t391 * t293;
t259 = qJD(1) * t266;
t201 = -t259 * t403 - t242;
t136 = t294 * t193 - t201 * t292;
t371 = t244 * t294;
t100 = pkin(4) * t245 - pkin(8) * t371 + t136;
t137 = t292 * t193 + t294 * t201;
t372 = t244 * t292;
t112 = -pkin(8) * t372 + t137;
t390 = pkin(8) + qJ(4);
t265 = t390 * t292;
t267 = t390 * t294;
t355 = qJD(5) * t303;
t357 = qJD(4) * t294;
t358 = qJD(4) * t292;
t462 = -t265 * t355 + (-t112 + t357) * t303 + (-qJD(5) * t267 - t100 - t358) * t299;
t206 = -t299 * t265 + t303 * t267;
t461 = -t257 * qJD(4) - qJD(5) * t206 - t303 * t100 + t112 * t299;
t352 = qJD(1) * qJD(2);
t271 = qJ(2) * qJDD(1) + t352;
t298 = sin(qJ(6));
t302 = cos(qJ(6));
t490 = -t153 * t298 + t302 * t476;
t94 = t153 * t302 + t298 * t476;
t489 = pkin(9) * t450 + t462;
t488 = -pkin(5) * t245 + pkin(9) * t492 + t461;
t360 = t293 ^ 2 + t295 ^ 2;
t278 = pkin(4) * t294 + pkin(3);
t290 = pkin(11) + qJ(5);
t284 = cos(t290);
t261 = pkin(5) * t284 + t278;
t286 = qJ(6) + t290;
t276 = sin(t286);
t277 = cos(t286);
t282 = sin(t290);
t485 = -m(6) * t278 - m(7) * t261 - t284 * mrSges(6,1) - t277 * mrSges(7,1) + t282 * mrSges(6,2) + t276 * mrSges(7,2);
t287 = -pkin(9) - t390;
t484 = -m(6) * t390 + m(7) * t287 - t491;
t190 = -qJD(3) * pkin(3) + qJD(4) - t201;
t279 = pkin(2) * t295 + pkin(1);
t264 = -qJD(1) * t279 + qJD(2);
t385 = Ifges(5,4) * t294;
t322 = -Ifges(5,2) * t292 + t385;
t386 = Ifges(5,4) * t292;
t323 = Ifges(5,1) * t294 - t386;
t325 = mrSges(5,1) * t292 + mrSges(5,2) * t294;
t404 = t294 / 0.2e1;
t460 = Ifges(4,5) * qJD(3);
t468 = t332 / 0.2e1;
t469 = t214 / 0.2e1;
t483 = t264 * mrSges(4,2) + t322 * t468 + t323 * t469 + (t214 * Ifges(5,1) + Ifges(5,4) * t332 - Ifges(5,5) * t244) * t404 - t292 * (t214 * Ifges(5,4) + Ifges(5,2) * t332 - Ifges(5,6) * t244) / 0.2e1 - t201 * mrSges(4,3) + t190 * t325 + t460 / 0.2e1;
t163 = -pkin(3) * t244 - qJ(4) * t245 + t264;
t202 = -t300 * t259 + t403 * t260;
t192 = qJD(3) * qJ(4) + t202;
t119 = t294 * t163 - t192 * t292;
t120 = t292 * t163 + t294 * t192;
t82 = -pkin(4) * t244 - pkin(8) * t214 + t119;
t99 = pkin(8) * t332 + t120;
t49 = -t299 * t99 + t303 * t82;
t35 = -pkin(9) * t153 + t49;
t34 = pkin(5) * t240 + t35;
t50 = t299 * t82 + t303 * t99;
t36 = pkin(9) * t476 + t50;
t375 = t298 * t36;
t13 = t302 * t34 - t375;
t374 = t302 * t36;
t14 = t298 * t34 + t374;
t459 = Ifges(4,6) * qJD(3);
t482 = -t264 * mrSges(4,1) - t119 * mrSges(5,1) - t49 * mrSges(6,1) - t13 * mrSges(7,1) + t120 * mrSges(5,2) + t50 * mrSges(6,2) + t14 * mrSges(7,2) + t202 * mrSges(4,3) + t459 / 0.2e1;
t481 = m(6) + m(7);
t452 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t332 - mrSges(5,2) * t214 - mrSges(4,3) * t245;
t301 = sin(qJ(1));
t304 = cos(qJ(1));
t479 = g(1) * t304 + g(2) * t301;
t232 = qJD(6) + t240;
t464 = t332 * Ifges(5,6);
t465 = t214 * Ifges(5,5);
t478 = t153 * Ifges(6,5) + t94 * Ifges(7,5) + Ifges(6,6) * t476 + Ifges(7,6) * t490 - t244 * Ifges(5,3) + t240 * Ifges(6,3) + t232 * Ifges(7,3) + t464 + t465;
t341 = m(3) * qJ(2) + mrSges(3,3);
t475 = -m(5) * t391 + mrSges(2,2) - mrSges(4,3) - t325 - t341;
t248 = t315 * qJD(3);
t197 = qJD(1) * t248 + qJDD(1) * t258;
t175 = qJDD(3) * t292 + t197 * t294;
t249 = t258 * qJD(3);
t351 = qJDD(1) * t293;
t198 = qJD(1) * t249 - qJDD(1) * t346 + t300 * t351;
t263 = -qJDD(1) * t279 + qJDD(2);
t109 = pkin(3) * t198 - qJ(4) * t197 - qJD(4) * t245 + t263;
t334 = pkin(7) * qJDD(1) + t271;
t233 = t334 * t293;
t234 = t334 * t295;
t337 = qJD(3) * t403;
t347 = -t300 * t233 + t403 * t234 - t259 * t337;
t128 = qJDD(3) * qJ(4) + (qJD(4) - t242) * qJD(3) + t347;
t64 = t294 * t109 - t128 * t292;
t42 = pkin(4) * t198 - pkin(8) * t175 + t64;
t174 = qJDD(3) * t294 - t197 * t292;
t65 = t292 * t109 + t294 * t128;
t55 = pkin(8) * t174 + t65;
t12 = -qJD(5) * t50 - t299 * t55 + t303 * t42;
t195 = qJDD(5) + t198;
t69 = qJD(5) * t476 + t174 * t299 + t175 * t303;
t6 = pkin(5) * t195 - pkin(9) * t69 + t12;
t356 = qJD(5) * t299;
t11 = t299 * t42 + t303 * t55 + t82 * t355 - t356 * t99;
t70 = -qJD(5) * t153 + t174 * t303 - t175 * t299;
t7 = pkin(9) * t70 + t11;
t2 = qJD(6) * t13 + t298 * t6 + t302 * t7;
t3 = -qJD(6) * t14 - t298 * t7 + t302 * t6;
t474 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t473 = -t12 * mrSges(6,1) + t11 * mrSges(6,2);
t291 = pkin(10) + qJ(3);
t283 = sin(t291);
t285 = cos(t291);
t328 = mrSges(4,1) * t285 - mrSges(4,2) * t283;
t329 = -mrSges(3,1) * t295 + mrSges(3,2) * t293;
t472 = m(3) * pkin(1) + t283 * t491 + mrSges(2,1) + t328 - t329;
t30 = qJD(6) * t490 + t298 * t70 + t302 * t69;
t440 = t30 / 0.2e1;
t31 = -qJD(6) * t94 - t298 * t69 + t302 * t70;
t439 = t31 / 0.2e1;
t432 = t69 / 0.2e1;
t431 = t70 / 0.2e1;
t419 = t174 / 0.2e1;
t418 = t175 / 0.2e1;
t189 = qJDD(6) + t195;
t417 = t189 / 0.2e1;
t416 = t195 / 0.2e1;
t415 = t198 / 0.2e1;
t204 = -t303 * t265 - t267 * t299;
t176 = -pkin(9) * t257 + t204;
t177 = -pkin(9) * t317 + t206;
t118 = t176 * t298 + t177 * t302;
t467 = -qJD(6) * t118 - t298 * t489 + t488 * t302;
t117 = t176 * t302 - t177 * t298;
t466 = qJD(6) * t117 + t488 * t298 + t302 * t489;
t441 = m(7) * pkin(5);
t463 = mrSges(6,1) + t441;
t326 = -mrSges(5,1) * t294 + mrSges(5,2) * t292;
t313 = m(5) * pkin(3) - t326;
t455 = t313 * t285;
t196 = -pkin(3) * t315 - qJ(4) * t258 - t279;
t207 = -t300 * t266 + t268 * t403;
t143 = t294 * t196 - t207 * t292;
t367 = t258 * t294;
t108 = -pkin(4) * t315 - pkin(8) * t367 + t143;
t144 = t292 * t196 + t294 * t207;
t368 = t258 * t292;
t121 = -pkin(8) * t368 + t144;
t61 = t299 * t108 + t303 * t121;
t114 = -t172 * t298 - t173 * t302;
t199 = -t257 * t298 - t302 * t317;
t138 = qJD(6) * t199 - t246 * t302 - t247 * t298;
t454 = t138 - t114;
t113 = -t172 * t302 + t173 * t298;
t200 = t257 * t302 - t298 * t317;
t139 = -qJD(6) * t200 + t246 * t298 - t247 * t302;
t453 = t139 - t113;
t449 = -t403 * t266 - t268 * t300;
t448 = -t292 * t64 + t294 * t65;
t158 = pkin(4) * t372 + t202;
t447 = -pkin(5) * t450 - t158;
t146 = -pkin(4) * t332 + t190;
t89 = -pkin(5) * t476 + t146;
t446 = -mrSges(7,1) * t89 + t14 * mrSges(7,3);
t445 = mrSges(7,2) * t89 - mrSges(7,3) * t13;
t443 = Ifges(7,4) * t440 + Ifges(7,2) * t439 + Ifges(7,6) * t417;
t442 = Ifges(7,1) * t440 + Ifges(7,4) * t439 + Ifges(7,5) * t417;
t438 = Ifges(6,4) * t432 + Ifges(6,2) * t431 + Ifges(6,6) * t416;
t437 = Ifges(6,1) * t432 + Ifges(6,4) * t431 + Ifges(6,5) * t416;
t400 = Ifges(7,4) * t94;
t44 = Ifges(7,2) * t490 + t232 * Ifges(7,6) + t400;
t436 = -t44 / 0.2e1;
t435 = t44 / 0.2e1;
t88 = Ifges(7,4) * t490;
t45 = t94 * Ifges(7,1) + t232 * Ifges(7,5) + t88;
t434 = -t45 / 0.2e1;
t433 = t45 / 0.2e1;
t428 = Ifges(5,1) * t418 + Ifges(5,4) * t419 + Ifges(5,5) * t415;
t427 = -t490 / 0.2e1;
t426 = t490 / 0.2e1;
t425 = -t94 / 0.2e1;
t424 = t94 / 0.2e1;
t423 = -t476 / 0.2e1;
t422 = t476 / 0.2e1;
t421 = -t153 / 0.2e1;
t420 = t153 / 0.2e1;
t414 = -t232 / 0.2e1;
t413 = t232 / 0.2e1;
t412 = -t240 / 0.2e1;
t411 = t240 / 0.2e1;
t410 = t244 / 0.2e1;
t409 = -t244 / 0.2e1;
t407 = t245 / 0.2e1;
t399 = pkin(4) * t292;
t398 = pkin(5) * t153;
t396 = pkin(5) * t282;
t393 = g(3) * t283;
t388 = mrSges(6,3) * t476;
t387 = mrSges(6,3) * t153;
t380 = t245 * Ifges(4,4);
t370 = t248 * t292;
t369 = t248 * t294;
t365 = t285 * t301;
t364 = t285 * t304;
t159 = pkin(3) * t249 - qJ(4) * t248 - qJD(4) * t258;
t167 = t315 * qJD(2) + qJD(3) * t449;
t107 = t292 * t159 + t294 * t167;
t216 = t276 * t365 + t277 * t304;
t217 = t276 * t304 - t277 * t365;
t362 = -t216 * mrSges(7,1) + t217 * mrSges(7,2);
t218 = -t276 * t364 + t277 * t301;
t219 = t276 * t301 + t277 * t364;
t361 = t218 * mrSges(7,1) - t219 * mrSges(7,2);
t359 = qJD(3) * t300;
t354 = m(5) + t481;
t350 = qJDD(1) * t295;
t349 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t189;
t348 = Ifges(6,5) * t69 + Ifges(6,6) * t70 + Ifges(6,3) * t195;
t340 = m(5) * qJ(4) + mrSges(5,3);
t37 = -t70 * mrSges(6,1) + t69 * mrSges(6,2);
t10 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t336 = t198 * mrSges(4,1) + t197 * mrSges(4,2);
t116 = -t174 * mrSges(5,1) + t175 * mrSges(5,2);
t60 = t303 * t108 - t121 * t299;
t106 = t294 * t159 - t167 * t292;
t171 = pkin(4) * t368 - t449;
t330 = -mrSges(3,1) * t350 + mrSges(3,2) * t351;
t324 = -mrSges(7,1) * t276 - mrSges(7,2) * t277;
t321 = Ifges(5,5) * t294 - Ifges(5,6) * t292;
t182 = t317 * t258;
t51 = -pkin(5) * t315 + pkin(9) * t182 + t60;
t181 = t257 * t258;
t54 = -pkin(9) * t181 + t61;
t21 = -t298 * t54 + t302 * t51;
t22 = t298 * t51 + t302 * t54;
t320 = t119 * t292 - t120 * t294;
t125 = -t181 * t302 + t182 * t298;
t126 = -t181 * t298 - t182 * t302;
t319 = t261 * t285 - t283 * t287;
t318 = t278 * t285 + t283 * t390;
t316 = t349 - t474;
t222 = -t282 * t364 + t284 * t301;
t220 = t282 * t365 + t284 * t304;
t76 = pkin(4) * t249 - pkin(8) * t369 + t106;
t87 = -pkin(8) * t370 + t107;
t25 = t108 * t355 - t121 * t356 + t299 * t76 + t303 * t87;
t135 = -t233 * t403 - t300 * t234 + t259 * t359 - t260 * t337;
t129 = -qJDD(3) * pkin(3) + qJDD(4) - t135;
t26 = -qJD(5) * t61 - t299 * t87 + t303 * t76;
t83 = -t174 * pkin(4) + t129;
t309 = t283 * t340 + t455;
t168 = qJD(2) * t258 + qJD(3) * t207;
t145 = pkin(4) * t370 + t168;
t280 = -qJDD(1) * pkin(1) + qJDD(2);
t269 = t304 * t279;
t262 = t396 + t399;
t236 = Ifges(4,4) * t244;
t229 = pkin(5) * t317 - t278;
t224 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t244;
t223 = t282 * t301 + t284 * t364;
t221 = t282 * t304 - t284 * t365;
t186 = t245 * Ifges(4,1) + t236 + t460;
t185 = t244 * Ifges(4,2) + t380 + t459;
t170 = -mrSges(5,1) * t244 - mrSges(5,3) * t214;
t169 = mrSges(5,2) * t244 + mrSges(5,3) * t332;
t134 = -t260 * t359 + t347;
t133 = mrSges(6,1) * t240 - t387;
t132 = -mrSges(6,2) * t240 + t388;
t131 = mrSges(5,1) * t198 - mrSges(5,3) * t175;
t130 = -mrSges(5,2) * t198 + mrSges(5,3) * t174;
t127 = pkin(5) * t181 + t171;
t124 = t246 * t258 - t248 * t257;
t123 = -t247 * t258 - t248 * t317;
t96 = -mrSges(6,1) * t476 + mrSges(6,2) * t153;
t85 = t175 * Ifges(5,4) + t174 * Ifges(5,2) + t198 * Ifges(5,6);
t73 = mrSges(7,1) * t232 - mrSges(7,3) * t94;
t72 = -mrSges(7,2) * t232 + mrSges(7,3) * t490;
t71 = -t124 * pkin(5) + t145;
t63 = -mrSges(6,2) * t195 + mrSges(6,3) * t70;
t62 = mrSges(6,1) * t195 - mrSges(6,3) * t69;
t52 = -mrSges(7,1) * t490 + mrSges(7,2) * t94;
t47 = -qJD(6) * t126 - t123 * t298 + t124 * t302;
t46 = qJD(6) * t125 + t123 * t302 + t124 * t298;
t38 = -t70 * pkin(5) + t83;
t24 = -mrSges(7,2) * t189 + mrSges(7,3) * t31;
t23 = mrSges(7,1) * t189 - mrSges(7,3) * t30;
t18 = pkin(9) * t124 + t25;
t17 = pkin(5) * t249 - pkin(9) * t123 + t26;
t16 = t302 * t35 - t375;
t15 = -t298 * t35 - t374;
t5 = -qJD(6) * t22 + t17 * t302 - t18 * t298;
t4 = qJD(6) * t21 + t17 * t298 + t18 * t302;
t1 = [(Ifges(7,1) * t46 + Ifges(7,4) * t47) * t424 + (Ifges(7,1) * t126 + Ifges(7,4) * t125) * t440 + m(4) * (t134 * t207 + t167 * t202 - t263 * t279) + m(5) * (t106 * t119 + t107 * t120 + t143 * t64 + t144 * t65) - (Ifges(5,5) * t175 + Ifges(5,6) * t174 + Ifges(5,3) * t198 + t348 + t349) * t315 / 0.2e1 + t83 * (mrSges(6,1) * t181 - mrSges(6,2) * t182) + (-t119 * t369 - t120 * t370 - t367 * t64 - t368 * t65) * mrSges(5,3) + (-m(4) * t201 + m(5) * t190 - t452) * t168 + 0.2e1 * t360 * t271 * mrSges(3,3) + (Ifges(7,4) * t46 + Ifges(7,2) * t47) * t426 + (Ifges(7,4) * t126 + Ifges(7,2) * t125) * t439 + (Ifges(5,6) * t468 + Ifges(5,5) * t469 - Ifges(4,4) * t407 + Ifges(5,3) * t409 - Ifges(4,2) * t410 + Ifges(6,3) * t411 + Ifges(7,3) * t413 + Ifges(6,5) * t420 + Ifges(6,6) * t422 + Ifges(7,5) * t424 + Ifges(7,6) * t426 - t185 / 0.2e1 + t478 / 0.2e1 - t482) * t249 + (mrSges(4,2) * t263 - mrSges(4,3) * t135 + Ifges(4,1) * t197 - Ifges(4,4) * t198 + Ifges(4,5) * qJDD(3) + t129 * t325 + t321 * t415 + t322 * t419 + t323 * t418) * t258 + (-Ifges(6,4) * t182 - Ifges(6,2) * t181) * t431 + m(7) * (t127 * t38 + t13 * t5 + t14 * t4 + t2 * t22 + t21 * t3 + t71 * t89) + m(6) * (t11 * t61 + t12 * t60 + t145 * t146 + t171 * t83 + t25 * t50 + t26 * t49) - (-m(4) * t135 + m(5) * t129 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t197 + t116) * t449 + (-t11 * t181 + t12 * t182 - t123 * t49 + t124 * t50) * mrSges(6,3) + (-Ifges(6,1) * t182 - Ifges(6,4) * t181) * t432 + (Ifges(6,5) * t123 + Ifges(6,6) * t124) * t411 + (Ifges(3,4) * t293 + Ifges(3,2) * t295) * t350 + (Ifges(3,1) * t293 + Ifges(3,4) * t295) * t351 + (-Ifges(6,5) * t182 - Ifges(6,6) * t181) * t416 + (Ifges(7,5) * t46 + Ifges(7,6) * t47) * t413 + (Ifges(7,5) * t126 + Ifges(7,6) * t125) * t417 + (-Ifges(6,5) * t432 - Ifges(6,3) * t416 - Ifges(6,6) * t431 - t263 * mrSges(4,1) + Ifges(4,6) * qJDD(3) + Ifges(4,4) * t197 - Ifges(4,2) * t198 - Ifges(5,3) * t415 - Ifges(7,3) * t417 - Ifges(5,5) * t418 - Ifges(5,6) * t419 - Ifges(7,6) * t439 - Ifges(7,5) * t440 - t64 * mrSges(5,1) + t65 * mrSges(5,2) + t134 * mrSges(4,3) + t473 + t474) * t315 + (-t221 * mrSges(6,1) - t217 * mrSges(7,1) - t220 * mrSges(6,2) - t216 * mrSges(7,2) + (-m(4) * t391 - m(6) * (t391 + t399) - m(7) * (t262 + t391) + t475) * t304 + (m(4) * t279 - m(6) * (-t279 - t318) - m(7) * (-t279 - t319) - m(5) * (-qJ(4) * t283 - t279) + t283 * mrSges(5,3) + t455 + t472) * t301) * g(1) + (Ifges(6,1) * t123 + Ifges(6,4) * t124) * t420 + (t125 * t2 - t126 * t3 - t13 * t46 + t14 * t47) * mrSges(7,3) + (t186 / 0.2e1 + Ifges(4,1) * t407 + t321 * t409 + Ifges(4,4) * t410 + t483) * t248 + (Ifges(6,4) * t123 + Ifges(6,2) * t124) * t422 + Ifges(2,3) * qJDD(1) + m(3) * (-pkin(1) * t280 + (t271 + t352) * qJ(2) * t360) + t167 * t224 + t207 * (-qJDD(3) * mrSges(4,2) - mrSges(4,3) * t198) + t107 * t169 + t106 * t170 + t171 * t37 + t145 * t96 + t146 * (-mrSges(6,1) * t124 + mrSges(6,2) * t123) + t143 * t131 + t144 * t130 + t25 * t132 + t26 * t133 + t38 * (-mrSges(7,1) * t125 + mrSges(7,2) * t126) + t127 * t10 + t21 * t23 + t22 * t24 - t279 * t336 + (-m(5) * t269 - t223 * mrSges(6,1) - t219 * mrSges(7,1) - t222 * mrSges(6,2) - t218 * mrSges(7,2) + (-m(4) - t481) * (t301 * t391 + t269) + (-m(6) * t399 - m(7) * t262 + t475) * t301 + (-m(6) * t318 - m(7) * t319 - t309 - t472) * t304) * g(2) + t280 * t329 - pkin(1) * t330 + t367 * t428 + t123 * t429 + t124 * t430 + t46 * t433 + t47 * t435 - t182 * t437 - t181 * t438 + t126 * t442 + t125 * t443 + t60 * t62 + t61 * t63 + t71 * t52 + t4 * t72 + t5 * t73 + t89 * (-mrSges(7,1) * t47 + mrSges(7,2) * t46) - t85 * t368 / 0.2e1; (-t52 - t96 + t452) * t245 - (t169 * t294 - t170 * t292 + t224) * t244 + m(3) * t280 + t453 * t73 + t454 * t72 + t450 * t133 - t492 * t132 + t336 + t294 * t131 + t292 * t130 + t257 * t63 - t317 * t62 + t200 * t24 + t199 * t23 + t330 + (-g(1) * t301 + g(2) * t304) * (m(3) + m(4) + t354) - t341 * t360 * qJD(1) ^ 2 + (t13 * t453 + t14 * t454 + t199 * t3 + t2 * t200 - t245 * t89) * m(7) + (t11 * t257 - t12 * t317 - t146 * t245 + t450 * t49 - t492 * t50) * m(6) + (-t190 * t245 + t244 * t320 + t292 * t65 + t294 * t64) * m(5) + (t201 * t245 - t202 * t244 + t263) * m(4); t83 * (mrSges(6,1) * t317 + mrSges(6,2) * t257) + (Ifges(6,5) * t257 - Ifges(6,6) * t317) * t416 + (Ifges(6,4) * t257 - Ifges(6,2) * t317) * t431 + (Ifges(6,1) * t257 - Ifges(6,4) * t317) * t432 - t317 * t438 + (-Ifges(6,5) * t246 - Ifges(6,6) * t247) * t411 + (-Ifges(6,1) * t246 - Ifges(6,4) * t247) * t420 + (-Ifges(6,4) * t246 - Ifges(6,2) * t247) * t422 + (-t309 - t328) * g(3) + (Ifges(7,4) * t114 + Ifges(7,2) * t113) * t427 + (-Ifges(6,1) * t173 - Ifges(6,4) * t172) * t421 + (-Ifges(6,4) * t173 - Ifges(6,2) * t172) * t423 + (-Ifges(6,5) * t173 - Ifges(6,6) * t172) * t412 + (t236 + t186) * t409 + (t130 * t294 - t131 * t292) * qJ(4) + (-Ifges(4,2) * t409 + Ifges(5,3) * t410 + Ifges(6,3) * t412 + Ifges(7,3) * t414 + Ifges(6,5) * t421 + Ifges(6,6) * t423 + Ifges(7,5) * t425 + Ifges(7,6) * t427 - t465 / 0.2e1 - t464 / 0.2e1 + t482) * t245 + (t484 * g(3) + t479 * (mrSges(4,1) + t313 - t485)) * t283 + (t485 * g(3) + t479 * (mrSges(4,2) - t340 + t484)) * t285 + (t357 - t137) * t169 - (Ifges(4,1) * t244 - t380 + t478) * t245 / 0.2e1 + (t119 * t371 + t120 * t372 + t448) * mrSges(5,3) + (Ifges(7,1) * t424 + Ifges(7,4) * t426 + Ifges(7,5) * t413 + t433 + t445) * t138 + (Ifges(7,4) * t424 + Ifges(7,2) * t426 + Ifges(7,6) * t413 + t435 + t446) * t139 + (-t136 - t358) * t170 + (-t11 * t317 - t12 * t257 + t450 * t50 + t49 * t492) * mrSges(6,3) + (-mrSges(6,1) * t450 - mrSges(6,2) * t492) * t146 + (Ifges(7,1) * t114 + Ifges(7,4) * t113) * t425 + (Ifges(7,5) * t114 + Ifges(7,6) * t113) * t414 + t450 * t430 - t492 * t429 + t129 * t326 + (-t113 * t14 + t114 * t13 + t199 * t2 - t200 * t3) * mrSges(7,3) - (-t321 * t410 + t483) * t244 + Ifges(4,3) * qJDD(3) - t278 * t37 + t229 * t10 - t201 * t224 + t206 * t63 + t38 * (-mrSges(7,1) * t199 + mrSges(7,2) * t200) + t204 * t62 + Ifges(4,5) * t197 - Ifges(4,6) * t198 - t158 * t96 - t134 * mrSges(4,2) + t135 * mrSges(4,1) + t447 * t52 + (-pkin(3) * t129 + qJ(4) * t448 - t320 * qJD(4) - t119 * t136 - t120 * t137 - t190 * t202) * m(5) + t452 * t202 + t461 * t133 + t462 * t132 + (t11 * t206 + t12 * t204 - t146 * t158 - t278 * t83 + t461 * t49 + t462 * t50) * m(6) + t85 * t404 + t185 * t407 + (Ifges(5,5) * t292 + Ifges(5,6) * t294) * t415 + (Ifges(7,5) * t200 + Ifges(7,6) * t199) * t417 + (Ifges(5,1) * t292 + t385) * t418 + (Ifges(5,2) * t294 + t386) * t419 + t292 * t428 + t114 * t434 + t113 * t436 + t257 * t437 + (Ifges(7,4) * t200 + Ifges(7,2) * t199) * t439 + (Ifges(7,1) * t200 + Ifges(7,4) * t199) * t440 + t200 * t442 + t199 * t443 + t466 * t72 + t467 * t73 + (t117 * t3 + t118 * t2 + t13 * t467 + t14 * t466 + t229 * t38 + t447 * t89) * m(7) - t89 * (-mrSges(7,1) * t113 + mrSges(7,2) * t114) - pkin(3) * t116 + t117 * t23 + t118 * t24; -t476 * t132 + t153 * t133 - t332 * t169 + t214 * t170 - t490 * t72 + t94 * t73 + t10 + t116 + t37 + (t13 * t94 - t14 * t490 + t38) * m(7) + (t153 * t49 - t476 * t50 + t83) * m(6) + (t119 * t214 - t120 * t332 + t129) * m(5) + (t285 * g(3) - t283 * t479) * t354; t316 + (t387 + t133) * t50 + (t388 - t132) * t49 + (-Ifges(6,2) * t153 + t150 + t80) * t423 + (m(7) * t396 + mrSges(6,1) * t282 + mrSges(6,2) * t284 - t324) * t393 + (Ifges(7,1) * t425 + Ifges(7,4) * t427 + Ifges(7,5) * t414 + t434 - t445) * t490 - (Ifges(7,4) * t425 + Ifges(7,2) * t427 + Ifges(7,6) * t414 + t436 - t446) * t94 + t348 + ((-t298 * t73 + t302 * t72) * qJD(6) + t23 * t302 + t24 * t298) * pkin(5) - t52 * t398 - m(7) * (t13 * t15 + t14 * t16 + t398 * t89) - t146 * (mrSges(6,1) * t153 + mrSges(6,2) * t476) + (-mrSges(6,2) * t221 + t220 * t463 - t362) * g(2) + (mrSges(6,2) * t223 - t222 * t463 - t361) * g(1) + (Ifges(6,5) * t476 - Ifges(6,6) * t153) * t412 + t79 * t420 + (Ifges(6,1) * t476 - t384) * t421 + (t2 * t298 + t3 * t302 + (-t13 * t298 + t14 * t302) * qJD(6)) * t441 - t16 * t72 - t15 * t73 - t473; -t89 * (mrSges(7,1) * t94 + mrSges(7,2) * t490) + (Ifges(7,1) * t490 - t400) * t425 + t44 * t424 + (Ifges(7,5) * t490 - Ifges(7,6) * t94) * t414 - t13 * t72 + t14 * t73 - g(1) * t361 - g(2) * t362 - t324 * t393 + (t13 * t490 + t14 * t94) * mrSges(7,3) + t316 + (-Ifges(7,2) * t94 + t45 + t88) * t427;];
tau  = t1;
