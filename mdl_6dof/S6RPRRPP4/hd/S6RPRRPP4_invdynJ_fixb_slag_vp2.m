% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:47
% EndTime: 2019-03-09 04:39:25
% DurationCPUTime: 26.16s
% Computational Cost: add. (11530->675), mult. (27788->834), div. (0->0), fcn. (21000->14), ass. (0->306)
t502 = mrSges(6,1) + mrSges(7,1);
t501 = -mrSges(6,2) + mrSges(7,3);
t488 = -mrSges(6,3) - mrSges(7,2);
t253 = sin(pkin(9));
t370 = pkin(7) + qJ(2);
t226 = t370 * t253;
t220 = qJD(1) * t226;
t254 = cos(pkin(9));
t227 = t370 * t254;
t221 = qJD(1) * t227;
t258 = sin(qJ(3));
t261 = cos(qJ(3));
t163 = -t220 * t261 - t258 * t221;
t155 = -qJD(3) * pkin(3) - t163;
t218 = t253 * t261 + t254 * t258;
t210 = t218 * qJD(1);
t257 = sin(qJ(4));
t260 = cos(qJ(4));
t173 = qJD(3) * t260 - t210 * t257;
t117 = -pkin(4) * t173 + qJD(5) + t155;
t216 = t253 * t258 - t261 * t254;
t208 = t216 * qJD(1);
t197 = t208 + qJD(4);
t354 = cos(pkin(10));
t252 = sin(pkin(10));
t241 = pkin(2) * t254 + pkin(1);
t224 = -qJD(1) * t241 + qJD(2);
t130 = pkin(3) * t208 - pkin(8) * t210 + t224;
t164 = -t258 * t220 + t261 * t221;
t156 = qJD(3) * pkin(8) + t164;
t88 = t130 * t257 + t156 * t260;
t78 = qJ(5) * t173 + t88;
t355 = t252 * t78;
t174 = qJD(3) * t257 + t210 * t260;
t87 = t260 * t130 - t156 * t257;
t77 = -qJ(5) * t174 + t87;
t70 = pkin(4) * t197 + t77;
t25 = t354 * t70 - t355;
t21 = -t197 * pkin(5) + qJD(6) - t25;
t119 = -t354 * t173 + t174 * t252;
t273 = t252 * t173 + t174 * t354;
t47 = pkin(5) * t119 - qJ(6) * t273 + t117;
t500 = -mrSges(6,2) * t117 - mrSges(7,2) * t21 + mrSges(6,3) * t25 + mrSges(7,3) * t47;
t72 = t354 * t78;
t26 = t252 * t70 + t72;
t22 = qJ(6) * t197 + t26;
t391 = t197 / 0.2e1;
t392 = -t197 / 0.2e1;
t399 = t273 / 0.2e1;
t400 = -t273 / 0.2e1;
t403 = -t119 / 0.2e1;
t499 = -t117 * mrSges(6,1) - t47 * mrSges(7,1) + Ifges(6,4) * t399 + Ifges(7,5) * t400 + Ifges(6,6) * t391 + Ifges(7,6) * t392 + (Ifges(6,2) + Ifges(7,3)) * t403 + mrSges(7,2) * t22 + mrSges(6,3) * t26;
t211 = t216 * qJD(3);
t161 = -qJD(1) * t211 + qJDD(1) * t218;
t115 = qJD(4) * t173 + qJDD(3) * t257 + t161 * t260;
t116 = -qJD(4) * t174 + qJDD(3) * t260 - t161 * t257;
t58 = t115 * t252 - t116 * t354;
t412 = -t58 / 0.2e1;
t59 = t115 * t354 + t252 * t116;
t410 = t59 / 0.2e1;
t212 = t218 * qJD(3);
t162 = -qJD(1) * t212 - qJDD(1) * t216;
t154 = qJDD(4) - t162;
t396 = t154 / 0.2e1;
t447 = Ifges(6,1) + Ifges(7,1);
t473 = Ifges(6,4) - Ifges(7,5);
t445 = Ifges(7,4) + Ifges(6,5);
t472 = Ifges(6,6) - Ifges(7,6);
t157 = pkin(3) * t210 + pkin(8) * t208;
t106 = t257 * t157 + t260 * t163;
t255 = -qJ(5) - pkin(8);
t305 = qJD(4) * t255;
t353 = t208 * t257;
t498 = -qJ(5) * t353 + qJD(5) * t260 + t257 * t305 - t106;
t105 = t260 * t157 - t163 * t257;
t352 = t208 * t260;
t497 = -pkin(4) * t210 - qJ(5) * t352 - qJD(5) * t257 + t260 * t305 - t105;
t323 = qJD(1) * qJD(2);
t232 = qJ(2) * qJDD(1) + t323;
t326 = qJD(4) * t257;
t496 = t326 + t353;
t402 = t119 / 0.2e1;
t441 = -t473 * t119 + t445 * t197 + t447 * t273;
t476 = t441 / 0.2e1;
t494 = Ifges(6,4) * t403 + Ifges(7,5) * t402 + t445 * t391 + t447 * t399 + t476 - t500;
t222 = -qJDD(1) * t241 + qJDD(2);
t96 = -pkin(3) * t162 - pkin(8) * t161 + t222;
t300 = pkin(7) * qJDD(1) + t232;
t191 = t300 * t253;
t192 = t300 * t254;
t328 = qJD(3) * t261;
t329 = qJD(3) * t258;
t103 = -t258 * t191 + t261 * t192 - t220 * t328 - t221 * t329;
t97 = qJDD(3) * pkin(8) + t103;
t24 = -t88 * qJD(4) - t257 * t97 + t260 * t96;
t11 = pkin(4) * t154 - qJ(5) * t115 - qJD(5) * t174 + t24;
t325 = qJD(4) * t260;
t23 = t130 * t325 - t156 * t326 + t257 * t96 + t260 * t97;
t13 = qJ(5) * t116 + qJD(5) * t173 + t23;
t3 = t11 * t354 - t252 * t13;
t2 = -t154 * pkin(5) + qJDD(6) - t3;
t411 = t58 / 0.2e1;
t104 = -t191 * t261 - t258 * t192 + t220 * t329 - t221 * t328;
t98 = -qJDD(3) * pkin(3) - t104;
t50 = -pkin(4) * t116 + qJDD(5) + t98;
t5 = pkin(5) * t58 - qJ(6) * t59 - qJD(6) * t273 + t50;
t493 = mrSges(6,2) * t50 + mrSges(7,2) * t2 - mrSges(6,3) * t3 - mrSges(7,3) * t5 + Ifges(7,5) * t411 + 0.2e1 * t396 * t445 + 0.2e1 * t410 * t447 + (Ifges(6,4) + t473) * t412;
t464 = Ifges(5,3) + Ifges(6,3) + Ifges(7,2);
t4 = t252 * t11 + t354 * t13;
t1 = qJ(6) * t154 + qJD(6) * t197 + t4;
t492 = mrSges(6,1) * t50 + mrSges(7,1) * t5 - mrSges(7,2) * t1 - mrSges(6,3) * t4 - t59 * Ifges(6,4) / 0.2e1 - t154 * Ifges(6,6) / 0.2e1 + 0.2e1 * Ifges(7,3) * t411 + (-t412 + t411) * Ifges(6,2) + (-t473 + Ifges(7,5)) * t410 + (-t472 + Ifges(7,6)) * t396;
t491 = Ifges(6,2) * t403 - Ifges(7,3) * t402 + t391 * t472 + t399 * t473 + t499;
t490 = m(7) + m(6);
t474 = pkin(4) * t490;
t487 = mrSges(5,1) + t474;
t304 = t354 * t257;
t217 = t252 * t260 + t304;
t133 = t217 * t208;
t207 = t217 * qJD(4);
t428 = t207 + t133;
t303 = t354 * t260;
t343 = t252 * t257;
t272 = t303 - t343;
t134 = t272 * t208;
t209 = t272 * qJD(4);
t427 = t209 + t134;
t330 = t253 ^ 2 + t254 ^ 2;
t250 = pkin(9) + qJ(3);
t244 = sin(t250);
t486 = t488 * t244;
t362 = t210 * mrSges(4,3);
t429 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t173 - mrSges(5,2) * t174 - t362;
t485 = -m(5) * t155 + t429;
t434 = Ifges(4,5) * qJD(3);
t483 = -t224 * mrSges(4,2) + t163 * mrSges(4,3) - t434 / 0.2e1;
t433 = Ifges(4,6) * qJD(3);
t482 = -t224 * mrSges(4,1) - t87 * mrSges(5,1) - t25 * mrSges(6,1) + t21 * mrSges(7,1) + t88 * mrSges(5,2) + t26 * mrSges(6,2) - t22 * mrSges(7,3) + t433 / 0.2e1;
t251 = qJ(4) + pkin(10);
t245 = sin(t251);
t247 = cos(t251);
t288 = -mrSges(5,1) * t260 + mrSges(5,2) * t257;
t481 = m(5) * pkin(3) + t501 * t245 + t247 * t502 - t288;
t478 = -m(5) - m(4);
t469 = t173 * Ifges(5,6);
t438 = -t498 * t252 + t497 * t354;
t435 = t497 * t252 + t498 * t354;
t430 = t496 * pkin(4) - t164;
t309 = t218 * t325;
t276 = -t211 * t257 + t309;
t246 = cos(t250);
t259 = sin(qJ(1));
t336 = t259 * t260;
t262 = cos(qJ(1));
t338 = t257 * t262;
t200 = -t246 * t338 + t336;
t468 = g(1) * t262 + g(2) * t259;
t307 = m(3) * qJ(2) + mrSges(3,3);
t463 = mrSges(2,2) - mrSges(4,3) - t307;
t462 = t174 * Ifges(5,5) - t472 * t119 + t464 * t197 + t445 * t273 + t469;
t461 = Ifges(5,5) * t115 + Ifges(5,6) * t116 + t464 * t154 + t445 * t59 - t472 * t58;
t459 = -Ifges(6,2) * t402 + Ifges(7,3) * t403 - t400 * t473 + t499;
t422 = m(7) * pkin(5) + t502;
t290 = t246 * mrSges(4,1) - t244 * mrSges(4,2);
t291 = -mrSges(3,1) * t254 + mrSges(3,2) * t253;
t457 = m(3) * pkin(1) + t244 * mrSges(5,3) + mrSges(2,1) + t290 - t291;
t456 = m(7) * qJ(6) + t501;
t453 = Ifges(6,4) * t402 + Ifges(7,5) * t403 + t445 * t392 + t400 * t447 + t500;
t452 = t24 * mrSges(5,1) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t23 * mrSges(5,2) - t4 * mrSges(6,2) + t1 * mrSges(7,3);
t405 = t115 / 0.2e1;
t404 = t116 / 0.2e1;
t41 = -mrSges(7,2) * t58 + mrSges(7,3) * t154;
t42 = -mrSges(6,2) * t154 - mrSges(6,3) * t58;
t443 = t41 + t42;
t43 = mrSges(6,1) * t154 - mrSges(6,3) * t59;
t44 = -t154 * mrSges(7,1) + t59 * mrSges(7,2);
t442 = t44 - t43;
t100 = -mrSges(6,2) * t197 - mrSges(6,3) * t119;
t99 = -mrSges(7,2) * t119 + mrSges(7,3) * t197;
t440 = t100 + t99;
t439 = pkin(5) * t428 - qJ(6) * t427 - qJD(6) * t217 + t430;
t437 = t210 * pkin(5) - t438;
t436 = -qJ(6) * t210 + t435;
t432 = t244 * t468;
t101 = mrSges(6,1) * t197 - mrSges(6,3) * t273;
t102 = -mrSges(7,1) * t197 + mrSges(7,2) * t273;
t431 = t102 - t101;
t160 = pkin(3) * t216 - pkin(8) * t218 - t241;
t169 = -t226 * t258 + t227 * t261;
t165 = t260 * t169;
t114 = t257 * t160 + t165;
t426 = -t261 * t226 - t227 * t258;
t425 = t23 * t260 - t24 * t257;
t413 = Ifges(5,1) * t405 + Ifges(5,4) * t404 + Ifges(5,5) * t396;
t395 = -t173 / 0.2e1;
t394 = -t174 / 0.2e1;
t393 = t174 / 0.2e1;
t389 = -t208 / 0.2e1;
t386 = t210 / 0.2e1;
t377 = pkin(4) * t174;
t376 = pkin(4) * t252;
t372 = g(3) * t244;
t278 = qJ(5) * t211 - qJD(5) * t218;
t131 = -t216 * qJD(2) + qJD(3) * t426;
t158 = pkin(3) * t212 + pkin(8) * t211;
t299 = -t131 * t257 + t260 * t158;
t32 = pkin(4) * t212 + t278 * t260 + (-t165 + (qJ(5) * t218 - t160) * t257) * qJD(4) + t299;
t314 = t260 * t131 + t257 * t158 + t160 * t325;
t36 = -qJ(5) * t309 + (-qJD(4) * t169 + t278) * t257 + t314;
t9 = t252 * t32 + t354 * t36;
t113 = t260 * t160 - t169 * t257;
t349 = t218 * t260;
t83 = pkin(4) * t216 - qJ(5) * t349 + t113;
t350 = t218 * t257;
t89 = -qJ(5) * t350 + t114;
t40 = t252 * t83 + t354 * t89;
t368 = mrSges(5,3) * t173;
t367 = mrSges(5,3) * t174;
t366 = Ifges(5,4) * t257;
t365 = Ifges(5,4) * t260;
t363 = t174 * Ifges(5,4);
t361 = t210 * Ifges(4,4);
t347 = t244 * t255;
t346 = t244 * t262;
t345 = t245 * t259;
t242 = pkin(4) * t260 + pkin(3);
t223 = t246 * t242;
t344 = t247 * t262;
t341 = t370 * t262;
t109 = t173 * Ifges(5,2) + t197 * Ifges(5,6) + t363;
t340 = t257 * t109;
t339 = t257 * t259;
t337 = t259 * t247;
t172 = Ifges(5,4) * t173;
t110 = t174 * Ifges(5,1) + t197 * Ifges(5,5) + t172;
t335 = t260 * t110;
t334 = t260 * t262;
t333 = t262 * t245;
t327 = qJD(4) * t218;
t322 = qJDD(1) * t253;
t321 = qJDD(1) * t254;
t317 = m(5) * pkin(8) + mrSges(5,3);
t313 = t354 * pkin(4);
t310 = t218 * t326;
t308 = t335 / 0.2e1;
t20 = t58 * mrSges(6,1) + t59 * mrSges(6,2);
t19 = t58 * mrSges(7,1) - t59 * mrSges(7,3);
t306 = -t326 / 0.2e1;
t302 = -t162 * mrSges(4,1) + t161 * mrSges(4,2);
t297 = t262 * t241 + t259 * t370;
t295 = pkin(3) * t246 + pkin(8) * t244;
t135 = pkin(4) * t350 - t426;
t292 = -mrSges(3,1) * t321 + mrSges(3,2) * t322;
t287 = mrSges(5,1) * t257 + mrSges(5,2) * t260;
t284 = Ifges(5,1) * t260 - t366;
t283 = -Ifges(5,2) * t257 + t365;
t282 = Ifges(5,5) * t260 - Ifges(5,6) * t257;
t281 = pkin(5) * t247 + qJ(6) * t245;
t128 = -mrSges(5,2) * t197 + t368;
t129 = mrSges(5,1) * t197 - t367;
t279 = t128 * t260 - t129 * t257;
t198 = t246 * t339 + t334;
t8 = -t252 * t36 + t32 * t354;
t39 = -t252 * t89 + t354 * t83;
t275 = t211 * t260 + t310;
t274 = t155 * t287;
t132 = qJD(2) * t218 + qJD(3) * t169;
t95 = pkin(4) * t276 + t132;
t243 = -qJDD(1) * pkin(1) + qJDD(2);
t240 = -t313 - pkin(5);
t236 = qJ(6) + t376;
t228 = t255 * t260;
t201 = t246 * t334 + t339;
t199 = -t246 * t336 + t338;
t193 = Ifges(4,4) * t208;
t181 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t208;
t180 = t246 * t344 + t345;
t179 = t246 * t333 - t337;
t178 = t246 * t337 - t333;
t177 = t246 * t345 + t344;
t171 = -t228 * t354 + t255 * t343;
t170 = -t228 * t252 - t255 * t304;
t159 = -pkin(5) * t272 - qJ(6) * t217 - t242;
t142 = t210 * Ifges(4,1) - t193 + t434;
t141 = -t208 * Ifges(4,2) + t361 + t433;
t140 = t272 * t218;
t139 = t217 * t218;
t91 = -t211 * t272 - t217 * t327;
t90 = t211 * t217 + t252 * t310 - t303 * t327;
t81 = -mrSges(5,2) * t154 + mrSges(5,3) * t116;
t80 = mrSges(5,1) * t154 - mrSges(5,3) * t115;
t75 = mrSges(6,1) * t119 + mrSges(6,2) * t273;
t74 = mrSges(7,1) * t119 - mrSges(7,3) * t273;
t71 = pkin(5) * t139 - qJ(6) * t140 + t135;
t62 = -mrSges(5,1) * t116 + mrSges(5,2) * t115;
t61 = pkin(5) * t273 + qJ(6) * t119 + t377;
t49 = -qJD(4) * t114 + t299;
t48 = -t169 * t326 + t314;
t45 = t115 * Ifges(5,4) + t116 * Ifges(5,2) + t154 * Ifges(5,6);
t38 = -t216 * pkin(5) - t39;
t37 = qJ(6) * t216 + t40;
t28 = t354 * t77 - t355;
t27 = t252 * t77 + t72;
t18 = -pkin(5) * t90 - qJ(6) * t91 - qJD(6) * t140 + t95;
t7 = -t212 * pkin(5) - t8;
t6 = qJ(6) * t212 + qJD(6) * t216 + t9;
t10 = [(Ifges(3,4) * t253 + Ifges(3,2) * t254) * t321 + (Ifges(3,1) * t253 + Ifges(3,4) * t254) * t322 + m(5) * (t113 * t24 + t114 * t23 + t48 * t88 + t49 * t87) + m(4) * (t103 * t169 + t131 * t164 - t222 * t241) + (-t201 * mrSges(5,1) - t200 * mrSges(5,2) + t488 * t346 + t478 * t297 - t490 * (pkin(4) * t339 + t262 * t223 - t255 * t346 + t297) - t422 * t180 - t456 * t179 + t463 * t259 + (-m(5) * t295 - t457) * t262) * g(2) + (t222 * mrSges(4,1) - t103 * mrSges(4,3) - Ifges(4,4) * t161 + Ifges(5,5) * t405 - Ifges(4,2) * t162 - Ifges(4,6) * qJDD(3) + Ifges(5,6) * t404 + Ifges(6,6) * t412 + Ifges(7,6) * t411 + t396 * t464 + t410 * t445 + t452 + t461 / 0.2e1) * t216 + (-Ifges(4,1) * t386 - Ifges(4,4) * t389 - t308 - t142 / 0.2e1 + t340 / 0.2e1 + t483) * t211 + (-m(4) * t163 - t485) * t132 + (-Ifges(5,5) * t275 - Ifges(5,6) * t276) * t391 + (-Ifges(5,1) * t275 - Ifges(5,4) * t276) * t393 + (t445 * t399 + t464 * t391 + t469 / 0.2e1 - t164 * mrSges(4,3) - t141 / 0.2e1 - Ifges(4,4) * t386 - Ifges(4,2) * t389 + Ifges(5,5) * t393 + Ifges(7,6) * t402 + Ifges(6,6) * t403 + t462 / 0.2e1 - t482) * t212 - (-m(4) * t104 + m(5) * t98 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t161 + t62) * t426 + (mrSges(4,2) * t222 - mrSges(4,3) * t104 + Ifges(4,1) * t161 + Ifges(4,4) * t162 + Ifges(4,5) * qJDD(3) + t110 * t306 + t282 * t396 + t283 * t404 + t284 * t405 + t287 * t98) * t218 + t349 * t413 + m(7) * (t1 * t37 + t18 * t47 + t2 * t38 + t21 * t7 + t22 * t6 + t5 * t71) + m(6) * (t117 * t95 + t135 * t50 + t25 * t8 + t26 * t9 + t3 * t39 + t4 * t40) + t173 * (-Ifges(5,4) * t275 - Ifges(5,2) * t276) / 0.2e1 + (-t23 * t350 - t24 * t349 + t275 * t87 - t276 * t88) * mrSges(5,3) + (-m(5) * t341 - t199 * mrSges(5,1) - t198 * mrSges(5,2) - t490 * (pkin(4) * t338 + t259 * t347 + t341) + t422 * t178 + t456 * t177 + (-m(4) * t370 + t463) * t262 + (-m(5) * (-t241 - t295) + m(4) * t241 - t490 * (-t241 - t223) + t457 - t486) * t259) * g(1) + Ifges(2,3) * qJDD(1) + t131 * t181 + m(3) * (-pkin(1) * t243 + (t232 + t323) * qJ(2) * t330) + t169 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t162) + t135 * t20 + t48 * t128 + t49 * t129 + t113 * t80 + t114 * t81 + t9 * t100 + t8 * t101 + t7 * t102 + t95 * t75 + t6 * t99 + t71 * t19 + t18 * t74 + t39 * t43 + t38 * t44 + t37 * t41 + t40 * t42 - t45 * t350 / 0.2e1 + 0.2e1 * t330 * t232 * mrSges(3,3) + t155 * (mrSges(5,1) * t276 - mrSges(5,2) * t275) + t243 * t291 - pkin(1) * t292 + t491 * t90 + t492 * t139 + t493 * t140 + t494 * t91 - t241 * t302 - t109 * t309 / 0.2e1; t302 + t292 + t443 * t217 - t442 * t272 + (-t74 - t75 + t429) * t210 + t260 * t80 + m(3) * t243 + t257 * t81 - (-t181 - t279) * t208 + t279 * qJD(4) + t431 * t428 + t440 * t427 + (-g(1) * t259 + g(2) * t262) * (m(3) + t490 - t478) - t307 * t330 * qJD(1) ^ 2 + (t1 * t217 - t2 * t272 + t21 * t428 - t210 * t47 + t22 * t427) * m(7) + (-t117 * t210 + t217 * t4 - t25 * t428 + t26 * t427 + t272 * t3) * m(6) + (-t155 * t210 + t23 * t257 + t24 * t260 + t197 * (-t257 * t87 + t260 * t88)) * m(5) + (t163 * t210 + t164 * t208 + t222) * m(4); (t274 + t308) * qJD(4) - (-t441 / 0.2e1 + t453) * t134 + (-t129 * t325 - t128 * t326 + m(5) * ((-t88 * t257 - t87 * t260) * qJD(4) + t425) - t257 * t80 + t260 * t81) * pkin(8) + t468 * ((t255 * t490 + mrSges(4,2) - t317 + t488) * t246 + (mrSges(4,1) + m(6) * t242 - m(7) * (-t242 - t281) + t481) * t244) + t435 * t100 + t436 * t99 + t437 * t102 + (t117 * t430 - t170 * t3 + t171 * t4 - t242 * t50 + t25 * t438 + t26 * t435) * m(6) + t438 * t101 + (t1 * t171 + t159 * t5 + t170 * t2 + t21 * t437 + t22 * t436 + t439 * t47) * m(7) + t439 * t74 - (t283 * t395 + t284 * t394 - t274 + t483) * t208 + (t362 + t485) * t164 + t141 * t386 + t340 * t389 + t442 * t170 + t443 * t171 + (Ifges(5,5) * t394 + Ifges(5,6) * t395 + Ifges(6,6) * t402 + Ifges(7,6) * t403 + t445 * t400 + t482) * t210 + t257 * t413 + (Ifges(5,5) * t257 + Ifges(5,6) * t260) * t396 + (Ifges(5,2) * t260 + t366) * t404 + (Ifges(5,1) * t257 + t365) * t405 + (-Ifges(4,2) * t210 + t142 - t193 + t335) * t208 / 0.2e1 + t260 * t45 / 0.2e1 + t109 * t306 + (t173 * t283 + t174 * t284 + t197 * t282) * qJD(4) / 0.2e1 + t430 * t75 + (-pkin(3) * t98 - t105 * t87 - t106 * t88) * m(5) + (-t244 * t317 - t290 - t490 * (t223 - t347) + (-m(7) * t281 - t481) * t246 + t486) * g(3) - t242 * t20 - t163 * t181 + t159 * t19 + Ifges(4,5) * t161 + Ifges(4,6) * t162 - t106 * t128 - t105 * t129 - t103 * mrSges(4,2) + t104 * mrSges(4,1) - pkin(3) * t62 + Ifges(4,3) * qJDD(3) - t459 * t133 - (-Ifges(4,1) * t208 - t361 + t462) * t210 / 0.2e1 + (t133 * t472 - t208 * t282 + t210 * t464) * t392 + t98 * t288 - t491 * t207 - t492 * t272 + t493 * t217 + t494 * t209 + (-t496 * t88 + (-t325 - t352) * t87 + t425) * mrSges(5,3); (t1 * t236 + t2 * t240 - t21 * t27 - t47 * t61 + (qJD(6) - t28) * t22) * m(7) + (t368 - t128) * t87 + t43 * t313 + t452 - t440 * t28 + t42 * t376 + (Ifges(5,5) * t173 - Ifges(5,6) * t174) * t392 + t109 * t393 + (Ifges(5,1) * t173 - t363) * t394 + ((t252 * t4 + t3 * t354) * pkin(4) - t117 * t377 + t25 * t27 - t26 * t28) * m(6) + (t367 + t129) * t88 + t461 - t431 * t27 + (-Ifges(5,2) * t174 + t110 + t172) * t395 + t240 * t44 + t236 * t41 - t155 * (mrSges(5,1) * t174 + mrSges(5,2) * t173) - t75 * t377 + qJD(6) * t99 - t61 * t74 + (-mrSges(5,2) * t199 + t177 * t422 - t178 * t456 + t198 * t487) * g(2) + (mrSges(5,2) * t201 + t422 * t179 - t180 * t456 - t200 * t487) * g(1) + (-t392 * t472 + t459) * t273 + (t245 * t422 - t247 * t456 + t257 * t474 + t287) * t372 + (-t453 + t476) * t119; t440 * t119 - t431 * t273 + t490 * t246 * g(3) + t19 + t20 + (t119 * t22 - t21 * t273 - t432 + t5) * m(7) + (t119 * t26 + t25 * t273 - t432 + t50) * m(6); t273 * t74 - t197 * t99 + (-g(1) * t179 - g(2) * t177 - t197 * t22 - t245 * t372 + t273 * t47 + t2) * m(7) + t44;];
tau  = t10;
