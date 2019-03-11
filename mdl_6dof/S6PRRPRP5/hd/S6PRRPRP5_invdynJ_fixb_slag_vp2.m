% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:33
% EndTime: 2019-03-08 21:46:05
% DurationCPUTime: 20.63s
% Computational Cost: add. (4547->665), mult. (10080->870), div. (0->0), fcn. (6884->10), ass. (0->290)
t476 = Ifges(6,1) + Ifges(7,1);
t465 = Ifges(7,4) + Ifges(6,5);
t477 = Ifges(6,6) - Ifges(7,6);
t467 = -mrSges(5,2) + mrSges(4,1);
t464 = Ifges(4,5) - Ifges(5,4);
t463 = Ifges(5,5) - Ifges(4,6);
t475 = Ifges(7,2) + Ifges(6,3);
t218 = sin(qJ(5));
t221 = cos(qJ(5));
t440 = t218 * t465 + t221 * t477;
t388 = Ifges(7,5) * t221;
t391 = Ifges(6,4) * t221;
t437 = t218 * t476 - t388 + t391;
t219 = sin(qJ(3));
t222 = cos(qJ(3));
t257 = pkin(9) * t219 - qJ(4) * t222;
t354 = qJD(3) * t219;
t291 = pkin(3) * t354 - qJD(4) * t219;
t108 = qJD(3) * t257 + t291;
t216 = sin(pkin(6));
t220 = sin(qJ(2));
t223 = cos(qJ(2));
t366 = t219 * t223;
t112 = (t218 * t366 + t220 * t221) * t216;
t375 = qJ(4) * t219;
t304 = -pkin(2) - t375;
t416 = pkin(3) + pkin(9);
t154 = -t222 * t416 + t304;
t352 = qJD(3) * t222;
t415 = pkin(4) + pkin(8);
t167 = t415 * t352;
t188 = t415 * t219;
t349 = qJD(5) * t221;
t350 = qJD(5) * t218;
t460 = -qJD(1) * t112 + t221 * t108 - t154 * t350 + t218 * t167 + t188 * t349;
t360 = qJD(1) * t216;
t324 = t223 * t360;
t287 = t219 * t324;
t325 = t220 * t360;
t452 = t221 * t154 + t218 * t188;
t459 = -qJD(5) * t452 + (t167 - t287) * t221 + (-t108 + t325) * t218;
t357 = qJD(2) * t216;
t316 = qJD(1) * t357;
t191 = t223 * t316;
t343 = qJDD(1) * t216;
t131 = t220 * t343 + t191;
t217 = cos(pkin(6));
t355 = qJD(3) * t217;
t473 = qJDD(2) * pkin(8) + qJD(1) * t355 + t131;
t346 = t219 * qJD(2);
t332 = mrSges(4,3) * t346;
t334 = mrSges(5,1) * t346;
t447 = -qJD(3) * t467 + t332 + t334;
t387 = Ifges(5,6) * t219;
t260 = -t222 * Ifges(5,3) - t387;
t345 = t222 * qJD(2);
t160 = qJD(3) * t218 + t221 * t345;
t156 = Ifges(6,4) * t160;
t318 = t218 * t345;
t353 = qJD(3) * t221;
t161 = -t318 + t353;
t203 = qJD(5) + t346;
t390 = Ifges(7,5) * t160;
t457 = t476 * t161 + t465 * t203 - t156 + t390;
t385 = t161 * Ifges(6,4);
t55 = -t160 * Ifges(6,2) + t203 * Ifges(6,6) + t385;
t472 = Ifges(5,5) * qJD(3) + qJD(2) * t260 + t218 * t457 + t221 * t55;
t344 = qJD(2) * qJD(3);
t168 = -t222 * qJDD(2) + t219 * t344;
t351 = qJD(5) * t160;
t64 = qJDD(3) * t221 + t168 * t218 - t351;
t419 = t64 / 0.2e1;
t65 = qJD(3) * t349 - qJD(5) * t318 + qJDD(3) * t218 - t221 * t168;
t417 = t65 / 0.2e1;
t471 = -m(6) - m(7);
t169 = qJDD(2) * t219 + t222 * t344;
t158 = qJDD(5) + t169;
t470 = (-Ifges(6,4) + Ifges(7,5)) * t65 + t476 * t64 + t465 * t158;
t171 = qJD(2) * pkin(8) + t325;
t342 = qJDD(1) * t217;
t33 = -t171 * t352 - t219 * t473 + t222 * t342;
t236 = qJDD(4) - t33;
t15 = pkin(4) * t169 - qJDD(3) * t416 + t236;
t190 = t220 * t316;
t130 = t223 * t343 - t190;
t114 = -qJDD(2) * pkin(2) - t130;
t226 = -qJ(4) * t169 - qJD(4) * t346 + t114;
t26 = t168 * t416 + t226;
t359 = qJD(1) * t217;
t198 = t222 * t359;
t88 = -t219 * (pkin(4) * qJD(2) + t171) + t198;
t454 = qJD(4) - t88;
t63 = -qJD(3) * t416 + t454;
t90 = qJD(2) * t154 - t324;
t3 = t218 * t15 + t221 * t26 + t63 * t349 - t350 * t90;
t1 = qJ(6) * t158 + qJD(6) * t203 + t3;
t469 = t1 * mrSges(7,3);
t468 = t3 * mrSges(6,2);
t414 = t158 / 0.2e1;
t466 = mrSges(6,3) + mrSges(7,2);
t462 = qJ(6) * t352 + qJD(6) * t219 + t460;
t461 = -pkin(5) * t352 - t459;
t41 = mrSges(6,1) * t158 - mrSges(6,3) * t64;
t42 = -t158 * mrSges(7,1) + t64 * mrSges(7,2);
t398 = t41 - t42;
t43 = -mrSges(6,2) * t158 - mrSges(6,3) * t65;
t44 = -mrSges(7,2) * t65 + mrSges(7,3) * t158;
t458 = t43 + t44;
t259 = pkin(5) * t221 + qJ(6) * t218;
t249 = -pkin(4) - t259;
t456 = -t198 - (qJD(2) * t249 - t171) * t219 + qJD(5) * t259 - qJD(6) * t221 + qJD(4);
t333 = mrSges(5,1) * t345;
t182 = -qJD(3) * mrSges(5,3) - t333;
t79 = mrSges(6,1) * t160 + mrSges(6,2) * t161;
t455 = t79 - t182;
t23 = -t218 * t90 + t221 * t63;
t453 = qJD(6) - t23;
t134 = mrSges(5,1) * t168 - qJDD(3) * mrSges(5,3);
t451 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t168 - t134;
t135 = t169 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t450 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t169 + t135;
t394 = Ifges(4,4) * t219;
t272 = t222 * Ifges(4,2) + t394;
t155 = Ifges(7,5) * t161;
t52 = t203 * Ifges(7,6) + t160 * Ifges(7,3) + t155;
t449 = Ifges(4,6) * qJD(3) + qJD(2) * t272 + t221 * t52;
t331 = mrSges(4,3) * t345;
t181 = -qJD(3) * mrSges(4,2) + t331;
t446 = t181 - t182;
t78 = mrSges(7,1) * t160 - mrSges(7,3) * t161;
t339 = -t78 - t455;
t445 = t181 - t339;
t444 = t219 * t440 + t222 * t475;
t443 = t219 * t437 + t222 * t465;
t386 = Ifges(5,6) * t222;
t442 = t219 * (-Ifges(5,2) * t222 + t387) + t222 * (Ifges(5,3) * t219 - t386);
t174 = -pkin(3) * t222 + t304;
t107 = qJD(2) * t174 - t324;
t172 = -qJD(2) * pkin(2) - t324;
t441 = t172 * (mrSges(4,1) * t219 + mrSges(4,2) * t222) + t107 * (-mrSges(5,2) * t219 - mrSges(5,3) * t222);
t439 = -t218 * t477 + t221 * t465;
t438 = t219 * t463 + t222 * t464;
t389 = Ifges(7,5) * t218;
t392 = Ifges(6,4) * t218;
t436 = t221 * t476 + t389 - t392;
t434 = t158 * t475 + t465 * t64 - t477 * t65;
t32 = -t171 * t354 + t219 * t342 + t222 * t473;
t433 = -t219 * t33 + t222 * t32;
t27 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4) - t32;
t28 = -qJDD(3) * pkin(3) + t236;
t432 = t219 * t28 - t222 * t27;
t24 = t218 * t63 + t221 * t90;
t4 = -qJD(5) * t24 + t15 * t221 - t218 * t26;
t431 = t218 * t3 + t221 * t4;
t2 = -pkin(5) * t158 + qJDD(6) - t4;
t430 = t1 * t218 - t2 * t221;
t429 = -t64 * Ifges(6,4) / 0.2e1 - t158 * Ifges(6,6) / 0.2e1 + Ifges(7,5) * t419 + Ifges(7,6) * t414 + (Ifges(6,2) + Ifges(7,3)) * t417;
t428 = -mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t205 = Ifges(4,4) * t345;
t427 = Ifges(4,1) * t346 + Ifges(4,5) * qJD(3) - t160 * t477 + t161 * t465 + t203 * t475 + t205;
t277 = mrSges(5,2) * t222 - mrSges(5,3) * t219;
t282 = mrSges(4,1) * t222 - mrSges(4,2) * t219;
t426 = t282 - t277 + mrSges(3,1);
t290 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t288 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t258 = -pkin(5) * t218 + qJ(6) * t221;
t173 = qJ(4) - t258;
t278 = t218 * mrSges(7,1) - t221 * mrSges(7,3);
t280 = mrSges(6,1) * t218 + mrSges(6,2) * t221;
t423 = -m(7) * t173 + mrSges(4,2) - mrSges(5,3) - t278 - t280 + (-m(5) - m(6)) * qJ(4);
t372 = t216 * t220;
t330 = t219 * t372;
t143 = -t217 * t222 + t330;
t376 = sin(pkin(10));
t295 = t376 * t223;
t377 = cos(pkin(10));
t298 = t377 * t220;
t139 = t217 * t298 + t295;
t300 = t216 * t377;
t80 = t139 * t219 + t222 * t300;
t296 = t376 * t220;
t297 = t377 * t223;
t141 = -t217 * t296 + t297;
t299 = t216 * t376;
t82 = t141 * t219 - t222 * t299;
t422 = g(1) * t82 + g(2) * t80 + g(3) * t143;
t225 = qJD(2) ^ 2;
t418 = -t65 / 0.2e1;
t413 = -t160 / 0.2e1;
t412 = t160 / 0.2e1;
t410 = t161 / 0.2e1;
t408 = t203 / 0.2e1;
t210 = t222 * pkin(8);
t396 = mrSges(6,3) * t160;
t395 = mrSges(6,3) * t161;
t393 = Ifges(4,4) * t222;
t206 = pkin(3) * t346;
t124 = qJD(2) * t257 + t206;
t106 = t222 * t171 + t219 * t359;
t89 = pkin(4) * t345 + t106;
t36 = t221 * t124 + t218 * t89;
t138 = -t217 * t297 + t296;
t374 = t138 * t222;
t140 = t217 * t295 + t298;
t373 = t140 * t222;
t371 = t216 * t223;
t370 = t218 * t219;
t369 = t218 * t222;
t367 = t219 * t221;
t100 = -mrSges(6,2) * t203 - t396;
t101 = -mrSges(7,2) * t160 + mrSges(7,3) * t203;
t363 = t100 + t101;
t102 = mrSges(6,1) * t203 - t395;
t103 = -mrSges(7,1) * t203 + mrSges(7,2) * t161;
t362 = -t102 + t103;
t361 = pkin(2) * t371 + pkin(8) * t372;
t189 = t222 * pkin(4) + t210;
t356 = qJD(2) * t220;
t348 = qJD(5) * t222;
t347 = qJD(5) * t416;
t329 = t221 * t371;
t328 = t222 * t371;
t323 = t216 * t356;
t322 = t223 * t357;
t321 = t218 * t348;
t320 = t218 * t347;
t319 = t221 * t347;
t303 = -t138 * pkin(2) + pkin(8) * t139;
t302 = -t140 * pkin(2) + pkin(8) * t141;
t294 = -t344 / 0.2e1;
t105 = t171 * t219 - t198;
t215 = qJD(3) * qJ(4);
t71 = t215 + t89;
t289 = t216 * qJ(4) * t366 + pkin(3) * t328 + t361;
t281 = mrSges(6,1) * t221 - mrSges(6,2) * t218;
t279 = mrSges(7,1) * t221 + mrSges(7,3) * t218;
t270 = -Ifges(6,2) * t218 + t391;
t269 = Ifges(6,2) * t221 + t392;
t263 = Ifges(7,3) * t218 + t388;
t262 = -Ifges(7,3) * t221 + t389;
t261 = -Ifges(5,2) * t219 - t386;
t12 = -pkin(5) * t203 + t453;
t13 = qJ(6) * t203 + t24;
t256 = t12 * t218 + t13 * t221;
t254 = t23 * t218 - t24 * t221;
t35 = -t124 * t218 + t221 * t89;
t75 = -t154 * t218 + t188 * t221;
t251 = -pkin(3) * t374 - t138 * t375 + t303;
t250 = -pkin(3) * t373 - t140 * t375 + t302;
t86 = t143 * t221 + t218 * t371;
t144 = t217 * t219 + t222 * t372;
t245 = t219 * (Ifges(4,1) * t222 - t394);
t232 = Ifges(6,6) * t222 + t219 * t269;
t231 = Ifges(7,6) * t222 + t219 * t262;
t16 = -pkin(4) * t168 - t27;
t229 = qJD(5) * t256 + t430;
t228 = -qJD(5) * t254 + t431;
t166 = t415 * t354;
t165 = -qJ(4) * t345 + t206;
t164 = t282 * qJD(2);
t163 = t277 * qJD(2);
t148 = Ifges(5,4) * qJD(3) + qJD(2) * t261;
t137 = -qJ(4) * t352 + t291;
t136 = t143 * pkin(3);
t104 = t222 * t259 + t189;
t94 = -t215 - t106;
t93 = -qJD(3) * pkin(3) + qJD(4) + t105;
t92 = -mrSges(5,2) * t168 - mrSges(5,3) * t169;
t91 = mrSges(4,1) * t168 + mrSges(4,2) * t169;
t87 = -t143 * t218 + t329;
t85 = -qJD(3) * t330 + (t322 + t355) * t222;
t84 = qJD(3) * t144 + t219 * t322;
t77 = pkin(5) * t161 + qJ(6) * t160;
t74 = t82 * pkin(3);
t73 = t80 * pkin(3);
t67 = -pkin(5) * t219 - t75;
t66 = qJ(6) * t219 + t452;
t45 = (qJD(5) * t258 + qJD(6) * t218) * t222 + (-pkin(8) + t249) * t354;
t39 = t140 * t218 - t82 * t221;
t37 = t138 * t218 - t80 * t221;
t34 = pkin(3) * t168 + t226;
t31 = -pkin(5) * t345 - t35;
t30 = qJ(6) * t345 + t36;
t29 = pkin(5) * t160 - qJ(6) * t161 + t71;
t22 = qJD(5) * t86 + t218 * t84 + t221 * t323;
t21 = -qJD(5) * t329 - t84 * t221 + (qJD(5) * t143 + t323) * t218;
t18 = mrSges(6,1) * t65 + mrSges(6,2) * t64;
t17 = mrSges(7,1) * t65 - mrSges(7,3) * t64;
t5 = pkin(5) * t65 - qJ(6) * t64 - qJD(6) * t161 + t16;
t6 = [m(2) * qJDD(1) - t458 * t87 + t398 * t86 + t447 * t84 + t363 * t22 + t362 * t21 + t450 * t143 + t445 * t85 + (t17 + t18 + t451) * t144 + (-m(2) - m(3) - m(4) - m(5) + t471) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t225 - t91 - t92) * t223 + (-mrSges(3,1) * t225 - mrSges(3,2) * qJDD(2) + (t163 - t164) * qJD(2)) * t220) * t216 + m(6) * (t144 * t16 - t21 * t23 + t22 * t24 - t3 * t87 + t4 * t86 + t71 * t85) + m(7) * (-t1 * t87 + t12 * t21 + t13 * t22 + t144 * t5 - t2 * t86 + t29 * t85) + m(5) * (t143 * t28 - t144 * t27 + t84 * t93 - t85 * t94 + (t107 * t356 - t223 * t34) * t216) + m(4) * (t105 * t84 + t106 * t85 - t143 * t33 + t144 * t32 + (-t114 * t223 + t172 * t356) * t216) + m(3) * (qJDD(1) * t217 ^ 2 + (t130 * t223 + t131 * t220) * t216); (t354 * t94 + t432) * mrSges(5,1) + (t231 * t412 + t232 * t413 + t408 * t444 + t410 * t443 + t441 + t438 * qJD(3) / 0.2e1) * qJD(3) + (-t131 + t191) * mrSges(3,2) + t2 * (-mrSges(7,1) * t219 - mrSges(7,2) * t369) + t4 * (mrSges(6,1) * t219 + mrSges(6,3) * t369) + t459 * t102 + t460 * t100 + t461 * t103 + (t1 * t66 + t104 * t5 + t12 * t461 + t13 * t462 + t2 * t67 + t29 * t45) * m(7) + t462 * t101 + (-m(6) * t71 - m(7) * t29 - t445) * t222 * t324 + (-mrSges(7,2) * t1 - mrSges(6,3) * t3 + t429) * t221 * t222 + (-t325 + t137) * t163 + t5 * t279 * t222 + t16 * t281 * t222 - t219 * t468 - t470 * t369 / 0.2e1 + (-t148 / 0.2e1 - t24 * mrSges(6,2) - t12 * mrSges(7,1) + t13 * mrSges(7,3) + t23 * mrSges(6,1) + t427 / 0.2e1 + t93 * mrSges(5,1) + t105 * mrSges(4,3)) * t352 + (mrSges(6,2) * t71 + t12 * mrSges(7,2) - t23 * mrSges(6,3) - mrSges(7,3) * t29) * (t218 * t354 - t221 * t348) + (-mrSges(6,1) * t71 - mrSges(7,1) * t29 + mrSges(7,2) * t13 + mrSges(6,3) * t24) * (t219 * t353 + t321) - (t218 * t52 + t221 * t457) * t348 / 0.2e1 + (t16 * t189 - t166 * t71 + t23 * t459 + t24 * t460 + t3 * t452 + t4 * t75) * m(6) + t452 * t43 + t164 * t325 + t55 * t321 / 0.2e1 + (-t263 * t412 - t270 * t413 - t408 * t439 - t410 * t436) * t348 + t442 * t294 - t449 * t354 / 0.2e1 + t451 * t210 + t169 * t219 * Ifges(4,1) + (t130 + t190) * mrSges(3,1) + t472 * t354 / 0.2e1 + t447 * (pkin(8) * t352 - t287) + t169 * t393 / 0.2e1 + (-Ifges(4,4) * t168 + Ifges(4,5) * qJDD(3) + t434) * t219 / 0.2e1 + (-m(4) * t302 - m(5) * t250 + t471 * (t141 * pkin(4) - pkin(9) * t373 + t250) - t290 * (-t140 * t370 + t141 * t221) - t288 * (t140 * t367 + t141 * t218) + t466 * t373 + t428 * t141 + t426 * t140) * g(1) + (-m(4) * t303 - m(5) * t251 + t471 * (t139 * pkin(4) - pkin(9) * t374 + t251) - t290 * (-t138 * t370 + t139 * t221) - t288 * (t138 * t367 + t139 * t218) + t466 * t374 + t428 * t139 + t426 * t138) * g(2) + (-m(4) * t361 - m(5) * t289 - t466 * t328 + t471 * (pkin(4) * t372 + pkin(9) * t328 + t289) - t290 * t112 - t288 * (t218 * t372 - t219 * t329) + (t220 * t428 - t223 * t426) * t216) * g(3) + (t219 * t475 - t440 * t222) * t414 + (t245 + t222 * (-Ifges(4,2) * t219 + t393)) * t344 / 0.2e1 + (t219 * t464 - t222 * t463) * qJDD(3) / 0.2e1 + (t219 * t465 - t222 * t437) * t419 + (-t106 * t354 + t433) * mrSges(4,3) + Ifges(3,3) * qJDD(2) + t66 * t44 + t67 * t42 + t75 * t41 + t45 * t78 - pkin(2) * t91 + t168 * t260 / 0.2e1 - t169 * t261 / 0.2e1 + t104 * t17 - t168 * t272 / 0.2e1 + t34 * t277 - t114 * t282 + t219 * t469 - t166 * t79 + t174 * t92 + t189 * t18 + (t107 * t137 + t174 * t34 - (t107 * t220 + (t219 * t93 - t222 * t94) * t223) * t360) * m(5) + (-pkin(2) * t114 - (t172 * t220 + (t105 * t219 + t106 * t222) * t223) * t360) * m(4) + (-t446 * t354 + t450 * t219 + ((t219 * t94 + t222 * t93) * qJD(3) + t432) * m(5) + ((t105 * t222 - t106 * t219) * qJD(3) + t433) * m(4)) * pkin(8) - t219 * (Ifges(5,4) * qJDD(3) - Ifges(5,2) * t169 + Ifges(5,6) * t168) / 0.2e1 + t222 * (Ifges(4,4) * t169 - Ifges(4,2) * t168 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t222 * (Ifges(5,5) * qJDD(3) - Ifges(5,6) * t169 + Ifges(5,3) * t168) / 0.2e1 + (Ifges(7,6) * t219 - t222 * t262) * t417 + (Ifges(6,6) * t219 - t222 * t269) * t418; ((-t231 / 0.2e1 + t232 / 0.2e1) * t160 - t24 * (-mrSges(6,2) * t222 + mrSges(6,3) * t367) - t13 * (mrSges(7,2) * t367 + mrSges(7,3) * t222) - t23 * (mrSges(6,1) * t222 - mrSges(6,3) * t370) - t12 * (-mrSges(7,1) * t222 + mrSges(7,2) * t370) - t441) * qJD(2) - t457 * t350 / 0.2e1 + (t269 / 0.2e1 - t262 / 0.2e1) * t351 + (-t55 / 0.2e1 + t52 / 0.2e1) * t349 + t148 * t345 / 0.2e1 + (t442 / 0.2e1 - t245 / 0.2e1) * t225 - t93 * t333 - t94 * t334 - (t161 * t443 + t203 * t444) * qJD(2) / 0.2e1 + (qJ(4) * t16 - t228 * t416 - t23 * t35 - t24 * t36 + t454 * t71) * m(6) + (-t12 * t31 - t13 * t30 + t173 * t5 - t229 * t416 + t29 * t456) * m(7) + (-t398 * t416 + t470 / 0.2e1) * t221 + (-t416 * t458 + t429) * t218 + (-t319 - t36) * t100 + (-t319 - t30) * t101 + (t320 - t35) * t102 + (-t320 - t31) * t103 - (t161 * t437 + t203 * t440) * qJD(5) / 0.2e1 + (-m(5) * t94 - t331 + t446) * t105 + (-m(5) * t93 + t332 - t447) * t106 + t449 * t346 / 0.2e1 + t455 * qJD(4) + t456 * t78 - t27 * mrSges(5,3) + t28 * mrSges(5,2) - t32 * mrSges(4,2) + t33 * mrSges(4,1) - t472 * t346 / 0.2e1 + t203 * (t279 * t29 + t281 * t71) + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) - (-Ifges(4,2) * t346 + t205 + t427) * t345 / 0.2e1 + (-pkin(3) * t28 - qJ(4) * t27 - qJD(4) * t94 - t107 * t165) * m(5) + (t18 - t134) * qJ(4) + (-t12 * t350 - t13 * t349 + t422 - t430) * mrSges(7,2) + (t23 * t350 - t24 * t349 + t422 - t431) * mrSges(6,3) + t436 * t419 + t438 * t294 + t439 * t414 + (m(5) * t73 + t467 * t80 + t471 * (-pkin(9) * t80 - t73) + t423 * (t139 * t222 - t219 * t300)) * g(2) + (m(5) * t136 + t471 * (-pkin(9) * t143 - t136) + t467 * t143 + t423 * t144) * g(3) + (m(5) * t74 + t467 * t82 + t471 * (-pkin(9) * t82 - t74) + t423 * (t141 * t222 + t219 * t299)) * g(1) + t463 * t168 + t464 * t169 - t88 * t79 + t5 * t278 + t16 * t280 - pkin(3) * t135 - t165 * t163 + t173 * t17 + t263 * t417 + t270 * t418; t163 * t346 + t339 * qJD(3) + (t203 * t363 + t398) * t221 + (t203 * t362 + t458) * t218 + t135 + (-qJD(3) * t29 + t256 * t346 + t229 - t422) * m(7) + (-qJD(3) * t71 - t254 * t346 + t228 - t422) * m(6) + (qJD(3) * t94 + t107 * t346 + t28 - t422) * m(5); (-Ifges(6,2) * t161 - t156 + t457) * t412 - (-t465 * t160 - t161 * t477) * t203 / 0.2e1 - (-t160 * t476 + t155 - t385 + t52) * t161 / 0.2e1 + (t12 * t160 + t13 * t161) * mrSges(7,2) + (-pkin(5) * t2 + qJ(6) * t1 - t12 * t24 + t13 * t453 - t29 * t77) * m(7) - pkin(5) * t42 + qJ(6) * t44 + t4 * mrSges(6,1) + t469 - t2 * mrSges(7,1) - t468 + t434 + (-t288 * (t138 * t221 + t218 * t80) + t290 * t37) * g(2) + (-t288 * (t140 * t221 + t218 * t82) + t290 * t39) * g(1) + (t288 * t87 - t290 * t86) * g(3) + (-t363 - t396) * t23 + (-t362 + t395) * t24 - t77 * t78 + qJD(6) * t101 - t29 * (mrSges(7,1) * t161 + mrSges(7,3) * t160) - t71 * (mrSges(6,1) * t161 - mrSges(6,2) * t160) + t55 * t410 + (Ifges(7,3) * t161 - t390) * t413; -t203 * t101 + t161 * t78 + (-g(1) * t39 - g(2) * t37 + g(3) * t86 - t13 * t203 + t29 * t161 + t2) * m(7) + t42;];
tau  = t6;
