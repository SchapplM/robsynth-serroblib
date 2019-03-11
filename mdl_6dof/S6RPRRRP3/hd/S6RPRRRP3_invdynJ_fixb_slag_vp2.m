% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP3
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:03:00
% EndTime: 2019-03-09 06:03:44
% DurationCPUTime: 25.59s
% Computational Cost: add. (9625->730), mult. (19966->937), div. (0->0), fcn. (13034->14), ass. (0->326)
t502 = mrSges(6,1) + mrSges(7,1);
t501 = mrSges(6,2) - mrSges(7,3);
t265 = cos(qJ(4));
t247 = pkin(4) * t265 + pkin(3);
t266 = cos(qJ(3));
t227 = t266 * t247;
t262 = sin(qJ(3));
t268 = -pkin(9) - pkin(8);
t366 = t262 * t268;
t500 = -mrSges(6,3) - mrSges(7,2);
t458 = t500 * t262;
t505 = m(6) + m(7);
t507 = t458 - t505 * (t227 - t366);
t357 = qJD(1) * t266;
t234 = qJD(4) - t357;
t229 = qJD(5) + t234;
t418 = t229 / 0.2e1;
t261 = sin(qJ(4));
t358 = qJD(1) * t262;
t197 = qJD(3) * t261 + t265 * t358;
t260 = sin(qJ(5));
t264 = cos(qJ(5));
t355 = qJD(3) * t265;
t287 = t261 * t358 - t355;
t275 = t264 * t197 - t260 * t287;
t427 = t275 / 0.2e1;
t123 = t260 * t197 + t264 * t287;
t430 = t123 / 0.2e1;
t431 = -t123 / 0.2e1;
t476 = Ifges(7,4) + Ifges(6,5);
t478 = Ifges(6,1) + Ifges(7,1);
t258 = sin(pkin(10));
t239 = pkin(1) * t258 + pkin(7);
t215 = t239 * qJD(1);
t164 = t266 * qJD(2) - t262 * t215;
t154 = -qJD(3) * pkin(3) - t164;
t116 = pkin(4) * t287 + t154;
t165 = t262 * qJD(2) + t266 * t215;
t155 = qJD(3) * pkin(8) + t165;
t315 = pkin(3) * t266 + pkin(8) * t262;
t293 = -pkin(2) - t315;
t259 = cos(pkin(10));
t414 = pkin(1) * t259;
t187 = t293 - t414;
t156 = t187 * qJD(1);
t88 = t265 * t155 + t261 * t156;
t75 = -pkin(9) * t287 + t88;
t390 = t260 * t75;
t87 = -t155 * t261 + t265 * t156;
t74 = -pkin(9) * t197 + t87;
t68 = pkin(4) * t234 + t74;
t25 = t264 * t68 - t390;
t493 = qJD(6) - t25;
t23 = -pkin(5) * t229 + t493;
t44 = t123 * pkin(5) - qJ(6) * t275 + t116;
t120 = Ifges(6,4) * t123;
t391 = Ifges(7,5) * t123;
t469 = t229 * t476 + t275 * t478 - t120 + t391;
t499 = mrSges(6,2) * t116 + mrSges(7,2) * t23 - mrSges(6,3) * t25 - t44 * mrSges(7,3) + t469 / 0.2e1;
t506 = Ifges(6,4) * t431 + Ifges(7,5) * t430 + t476 * t418 + t478 * t427 + t499;
t348 = qJD(1) * qJD(3);
t205 = qJDD(1) * t266 - t262 * t348;
t206 = qJDD(1) * t262 + t266 * t348;
t241 = -pkin(2) - t414;
t212 = t241 * qJDD(1);
t115 = -pkin(3) * t205 - pkin(8) * t206 + t212;
t351 = qJD(4) * t265;
t353 = qJD(4) * t261;
t356 = qJD(3) * t262;
t492 = qJD(2) * qJD(3) + t239 * qJDD(1);
t104 = t262 * qJDD(2) - t215 * t356 + t266 * t492;
t98 = qJDD(3) * pkin(8) + t104;
t27 = t261 * t115 - t155 * t353 + t156 * t351 + t265 * t98;
t504 = t27 * mrSges(5,2);
t28 = -qJD(4) * t88 + t265 * t115 - t261 * t98;
t503 = t28 * mrSges(5,1);
t477 = -Ifges(6,4) + Ifges(7,5);
t475 = Ifges(7,2) + Ifges(6,3);
t474 = -Ifges(6,6) + Ifges(7,6);
t384 = t264 * t75;
t26 = t260 * t68 + t384;
t24 = qJ(6) * t229 + t26;
t119 = Ifges(7,5) * t275;
t58 = Ifges(7,6) * t229 + Ifges(7,3) * t123 + t119;
t392 = Ifges(6,4) * t275;
t61 = -Ifges(6,2) * t123 + Ifges(6,6) * t229 + t392;
t498 = -t24 * mrSges(7,2) - t26 * mrSges(6,3) + t116 * mrSges(6,1) + t44 * mrSges(7,1) + t58 / 0.2e1 - t61 / 0.2e1;
t257 = qJ(4) + qJ(5);
t253 = sin(t257);
t254 = cos(t257);
t311 = -mrSges(5,1) * t265 + mrSges(5,2) * t261;
t497 = m(5) * pkin(3) - t501 * t253 + t254 * t502 - t311;
t397 = Ifges(4,4) * t262;
t305 = t266 * Ifges(4,2) + t397;
t495 = t23 * mrSges(7,1) + t26 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t305 / 0.2e1 - t24 * mrSges(7,3) - t25 * mrSges(6,1);
t354 = qJD(3) * t266;
t332 = t261 * t354;
t282 = t262 * t351 + t332;
t256 = qJ(1) + pkin(10);
t249 = sin(t256);
t250 = cos(t256);
t368 = t261 * t266;
t161 = t249 * t265 - t250 * t368;
t334 = t261 * t357;
t465 = -t165 + (-t334 + t353) * pkin(4);
t490 = -Ifges(6,2) * t431 + Ifges(7,3) * t430 + t418 * t474 + t427 * t477 + t498;
t419 = -t229 / 0.2e1;
t428 = -t275 / 0.2e1;
t489 = -Ifges(6,2) * t430 + Ifges(7,3) * t431 + t419 * t474 + t428 * t477 - t498;
t70 = pkin(5) * t275 + qJ(6) * t123;
t485 = Ifges(6,4) * t430 + Ifges(7,5) * t431 + t419 * t476 + t428 * t478 - t499;
t280 = t287 * qJD(4);
t117 = qJDD(3) * t261 + t206 * t265 - t280;
t118 = -qJD(4) * t197 + qJDD(3) * t265 - t206 * t261;
t40 = -qJD(5) * t123 + t264 * t117 + t260 * t118;
t442 = t40 / 0.2e1;
t41 = qJD(5) * t275 + t260 * t117 - t264 * t118;
t440 = t41 / 0.2e1;
t484 = -m(5) - m(4);
t433 = t117 / 0.2e1;
t432 = t118 / 0.2e1;
t192 = qJDD(4) - t205;
t185 = qJDD(5) + t192;
t423 = t185 / 0.2e1;
t422 = t192 / 0.2e1;
t481 = t287 / 0.2e1;
t405 = -qJD(1) / 0.2e1;
t479 = mrSges(3,2) - mrSges(4,3);
t473 = t185 * t476 + t40 * t478 + t41 * t477;
t31 = mrSges(6,1) * t185 - mrSges(6,3) * t40;
t32 = -t185 * mrSges(7,1) + t40 * mrSges(7,2);
t472 = t32 - t31;
t33 = -mrSges(6,2) * t185 - mrSges(6,3) * t41;
t34 = -mrSges(7,2) * t41 + mrSges(7,3) * t185;
t471 = t33 + t34;
t198 = t260 * t261 - t264 * t265;
t457 = qJD(4) + qJD(5);
t128 = t457 * t198;
t199 = t260 * t265 + t261 * t264;
t129 = t457 * t199;
t157 = t199 * t357;
t286 = t198 * t266;
t158 = qJD(1) * t286;
t470 = -qJD(6) * t199 + t465 + (t128 - t158) * qJ(6) + (t129 - t157) * pkin(5);
t93 = -mrSges(7,2) * t123 + mrSges(7,3) * t229;
t399 = mrSges(6,3) * t123;
t94 = -mrSges(6,2) * t229 - t399;
t403 = t93 + t94;
t398 = mrSges(6,3) * t275;
t95 = mrSges(6,1) * t229 - t398;
t96 = -mrSges(7,1) * t229 + mrSges(7,2) * t275;
t402 = -t95 + t96;
t67 = -mrSges(5,1) * t118 + mrSges(5,2) * t117;
t468 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t206 + t67;
t166 = t199 * t262;
t340 = mrSges(4,3) * t358;
t467 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t287 + t197 * mrSges(5,2) + t340;
t188 = Ifges(5,4) * t287;
t109 = t197 * Ifges(5,1) + t234 * Ifges(5,5) - t188;
t248 = Ifges(4,4) * t357;
t466 = Ifges(4,1) * t358 + Ifges(4,5) * qJD(3) + t265 * t109 + t248;
t364 = t265 * t266;
t201 = t239 * t364;
t127 = t261 * t187 + t201;
t352 = qJD(4) * t262;
t283 = -t261 * t352 + t265 * t354;
t464 = t185 * t475 + t40 * t476 + t41 * t474;
t372 = t253 * t266;
t151 = -t249 * t254 + t250 * t372;
t370 = t254 * t266;
t152 = t249 * t253 + t250 * t370;
t463 = t151 * t502 + t152 * t501;
t149 = t249 * t372 + t250 * t254;
t150 = t249 * t370 - t250 * t253;
t462 = t149 * t502 + t150 * t501;
t105 = qJDD(2) * t266 - t215 * t354 - t262 * t492;
t461 = t104 * t266 - t105 * t262;
t84 = mrSges(5,1) * t192 - mrSges(5,3) * t117;
t85 = -mrSges(5,2) * t192 + mrSges(5,3) * t118;
t460 = -t261 * t84 + t265 * t85;
t459 = -t261 * t28 + t265 * t27;
t456 = Ifges(5,5) * t197 - Ifges(5,6) * t287 + Ifges(5,3) * t234 + t123 * t474 + t229 * t475 + t275 * t476;
t313 = t266 * mrSges(4,1) - t262 * mrSges(4,2);
t455 = t262 * mrSges(5,3) + mrSges(3,1) + t313;
t19 = pkin(4) * t192 - pkin(9) * t117 + t28;
t21 = pkin(9) * t118 + t27;
t6 = -qJD(5) * t26 + t19 * t264 - t21 * t260;
t169 = t265 * t187;
t367 = t262 * t265;
t100 = -pkin(9) * t367 + t169 + (-t239 * t261 - pkin(4)) * t266;
t369 = t261 * t262;
t112 = -pkin(9) * t369 + t127;
t381 = t260 * t100 + t264 * t112;
t291 = pkin(4) * t262 - pkin(9) * t364;
t314 = pkin(3) * t262 - pkin(8) * t266;
t204 = t314 * qJD(3);
t333 = t239 * t356;
t359 = t265 * t204 + t261 * t333;
t56 = t291 * qJD(3) + (-t201 + (pkin(9) * t262 - t187) * t261) * qJD(4) + t359;
t76 = t187 * t351 + t261 * t204 + (-t262 * t355 - t266 * t353) * t239;
t65 = -pkin(9) * t282 + t76;
t11 = -qJD(5) * t381 - t260 * t65 + t264 * t56;
t454 = m(7) * pkin(5) + t502;
t452 = m(7) * qJ(6) - t501;
t349 = qJD(5) * t264;
t350 = qJD(5) * t260;
t5 = t260 * t19 + t264 * t21 + t68 * t349 - t350 * t75;
t2 = qJ(6) * t185 + qJD(6) * t229 + t5;
t441 = -t41 / 0.2e1;
t99 = -qJDD(3) * pkin(3) - t105;
t64 = -pkin(4) * t118 + t99;
t7 = pkin(5) * t41 - qJ(6) * t40 - qJD(6) * t275 + t64;
t446 = mrSges(6,1) * t64 + mrSges(7,1) * t7 - mrSges(7,2) * t2 - mrSges(6,3) * t5 + 0.2e1 * Ifges(7,3) * t440 - t40 * Ifges(6,4) / 0.2e1 - t185 * Ifges(6,6) / 0.2e1 + (t477 + Ifges(7,5)) * t442 + (t474 + Ifges(7,6)) * t423 + (-t441 + t440) * Ifges(6,2);
t445 = m(6) * pkin(4);
t439 = Ifges(5,1) * t433 + Ifges(5,4) * t432 + Ifges(5,5) * t422;
t421 = t197 / 0.2e1;
t263 = sin(qJ(1));
t413 = pkin(1) * t263;
t412 = pkin(4) * t197;
t411 = pkin(4) * t260;
t410 = pkin(4) * t264;
t267 = cos(qJ(1));
t255 = t267 * pkin(1);
t404 = qJD(4) / 0.2e1;
t202 = t314 * qJD(1);
t110 = -t164 * t261 + t265 * t202;
t86 = qJD(1) * t291 + t110;
t111 = t265 * t164 + t261 * t202;
t91 = -pkin(9) * t334 + t111;
t46 = t260 * t86 + t264 * t91;
t401 = mrSges(6,2) * t254;
t400 = mrSges(5,3) * t197;
t396 = Ifges(4,4) * t266;
t395 = Ifges(5,4) * t197;
t394 = Ifges(5,4) * t261;
t393 = Ifges(5,4) * t265;
t377 = t249 * t261;
t375 = t250 * t261;
t373 = t253 * t262;
t371 = t254 * t262;
t220 = t262 * t239;
t235 = pkin(4) * t369;
t172 = t220 + t235;
t343 = pkin(4) * t350;
t342 = pkin(4) * t349;
t341 = m(5) * pkin(8) + mrSges(5,3);
t339 = mrSges(4,3) * t357;
t337 = Ifges(5,5) * t117 + Ifges(5,6) * t118 + Ifges(5,3) * t192;
t207 = t239 * t354;
t131 = pkin(4) * t282 + t207;
t336 = t250 * pkin(2) + t249 * pkin(7) + t255;
t335 = qJD(4) * t268;
t108 = -Ifges(5,2) * t287 + Ifges(5,6) * t234 + t395;
t328 = -t261 * t108 / 0.2e1;
t322 = t250 * pkin(7) - t413;
t319 = -t149 * pkin(5) + qJ(6) * t150;
t318 = -t151 * pkin(5) + qJ(6) * t152;
t316 = t265 * t335;
t312 = mrSges(4,1) * t262 + mrSges(4,2) * t266;
t310 = mrSges(5,1) * t261 + mrSges(5,2) * t265;
t307 = Ifges(5,1) * t265 - t394;
t306 = Ifges(5,1) * t261 + t393;
t304 = -Ifges(5,2) * t261 + t393;
t303 = Ifges(5,2) * t265 + t394;
t302 = Ifges(4,5) * t266 - Ifges(4,6) * t262;
t301 = Ifges(5,5) * t265 - Ifges(5,6) * t261;
t300 = Ifges(5,5) * t261 + Ifges(5,6) * t265;
t299 = pkin(5) * t254 + qJ(6) * t253;
t45 = -t260 * t91 + t264 * t86;
t52 = t100 * t264 - t112 * t260;
t221 = t268 * t261;
t222 = t268 * t265;
t295 = t264 * t221 + t222 * t260;
t134 = t221 * t260 - t222 * t264;
t294 = t161 * pkin(4);
t159 = t249 * t368 + t250 * t265;
t10 = t100 * t349 - t112 * t350 + t260 * t56 + t264 * t65;
t290 = t154 * t310;
t289 = t241 * qJD(1) * t312;
t288 = t262 * (Ifges(4,1) * t266 - t397);
t284 = t159 * pkin(4);
t281 = t287 * mrSges(5,3);
t278 = Ifges(5,5) * t262 + t266 * t307;
t277 = Ifges(5,6) * t262 + t266 * t304;
t276 = Ifges(5,3) * t262 + t266 * t301;
t3 = -pkin(5) * t185 + qJDD(6) - t6;
t274 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) + t464;
t273 = (-t261 * t88 - t265 * t87) * qJD(4) + t459;
t246 = -pkin(5) - t410;
t240 = qJ(6) + t411;
t233 = qJD(6) + t342;
t228 = mrSges(7,3) * t371;
t225 = qJ(6) * t371;
t218 = -qJD(3) * mrSges(4,2) + t339;
t203 = t261 * t335;
t186 = t310 * t262;
t170 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t205;
t167 = t198 * t262;
t162 = t250 * t364 + t377;
t160 = -t249 * t364 + t375;
t148 = mrSges(5,1) * t234 - t400;
t147 = -t234 * mrSges(5,2) - t281;
t126 = -t239 * t368 + t169;
t121 = pkin(5) * t198 - qJ(6) * t199 - t247;
t82 = pkin(5) * t166 + qJ(6) * t167 + t172;
t81 = qJD(5) * t134 + t203 * t260 - t264 * t316;
t80 = qJD(5) * t295 + t264 * t203 + t260 * t316;
t79 = -t350 * t369 + (t367 * t457 + t332) * t264 + t283 * t260;
t78 = -qJD(3) * t286 - t166 * t457;
t77 = -qJD(4) * t127 + t359;
t72 = mrSges(6,1) * t123 + mrSges(6,2) * t275;
t71 = mrSges(7,1) * t123 - mrSges(7,3) * t275;
t55 = t412 + t70;
t49 = t117 * Ifges(5,4) + t118 * Ifges(5,2) + t192 * Ifges(5,6);
t48 = pkin(5) * t266 - t52;
t47 = -qJ(6) * t266 + t381;
t43 = -pkin(5) * t358 - t45;
t42 = qJ(6) * t358 + t46;
t30 = t264 * t74 - t390;
t29 = t260 * t74 + t384;
t22 = pkin(5) * t79 - qJ(6) * t78 + qJD(6) * t167 + t131;
t17 = mrSges(6,1) * t41 + mrSges(6,2) * t40;
t16 = mrSges(7,1) * t41 - mrSges(7,3) * t40;
t9 = -pkin(5) * t356 - t11;
t8 = qJ(6) * t356 - qJD(6) * t266 + t10;
t1 = [(-t27 * t369 - t28 * t367 - t282 * t88 - t283 * t87) * mrSges(5,3) + (-Ifges(7,5) * t167 - Ifges(7,6) * t266) * t440 + (-Ifges(6,4) * t167 - Ifges(6,6) * t266) * t441 + t3 * (mrSges(7,1) * t266 - mrSges(7,2) * t167) + t6 * (-mrSges(6,1) * t266 + mrSges(6,3) * t167) + m(6) * (t10 * t26 + t11 * t25 + t116 * t131 + t172 * t64 + t381 * t5 + t52 * t6) + t381 * t33 + (t266 * (-Ifges(4,2) * t262 + t396) + t288) * t348 / 0.2e1 - (t108 * t265 + t109 * t261) * t352 / 0.2e1 - (t337 + t464) * t266 / 0.2e1 + t466 * t354 / 0.2e1 + t467 * t207 + t206 * t396 / 0.2e1 + t446 * t166 + t266 * t504 - t266 * t503 + m(4) * (t212 * t241 + ((-t164 * t266 - t165 * t262) * qJD(3) + t461) * t239) + (-t167 * t478 - t266 * t476) * t442 + (t397 + t305) * t205 / 0.2e1 - t212 * t313 + (-t164 * t354 + t461) * mrSges(4,3) + t468 * t220 + t490 * t79 + (-t167 * t64 + t266 * t5) * mrSges(6,2) + (t167 * t7 - t2 * t266) * mrSges(7,3) + (m(3) * t413 + mrSges(2,1) * t263 - t160 * mrSges(5,1) + mrSges(2,2) * t267 - t159 * mrSges(5,2) + t484 * t322 - t505 * (pkin(4) * t375 + t249 * t366 + t322) + t479 * t250 + t454 * t150 + t452 * t149 + (m(4) * pkin(2) - m(5) * t293 - t505 * (-pkin(2) - t227) + t455 - t458) * t249) * g(1) + (-m(3) * t255 - mrSges(2,1) * t267 - t162 * mrSges(5,1) + mrSges(2,2) * t263 - t161 * mrSges(5,2) + t484 * t336 - t505 * (pkin(4) * t377 + t336) + t479 * t249 - t454 * t152 - t452 * t151 + (-m(5) * t315 - t455 + t507) * t250) * g(2) + t206 * Ifges(4,1) * t262 + m(7) * (t2 * t47 + t22 * t44 + t23 * t9 + t24 * t8 + t3 * t48 + t7 * t82) + t154 * (mrSges(5,1) * t282 + mrSges(5,2) * t283) + (-t165 * mrSges(4,3) + t476 * t427 + Ifges(7,6) * t430 + Ifges(6,6) * t431 + t456 / 0.2e1 + t475 * t418 - t88 * mrSges(5,2) + t87 * mrSges(5,1) - t495) * t356 + (m(3) * (t258 ^ 2 + t259 ^ 2) * pkin(1) ^ 2 + Ifges(3,3) + Ifges(2,3) + 0.2e1 * (mrSges(3,1) * t259 - mrSges(3,2) * t258) * pkin(1)) * qJDD(1) - t473 * t167 / 0.2e1 + (-t167 * t476 - t266 * t475) * t423 + qJD(3) ^ 2 * t302 / 0.2e1 + (-Ifges(5,6) * t266 + t262 * t304) * t432 + (-Ifges(5,5) * t266 + t262 * t307) * t433 + t367 * t439 - t49 * t369 / 0.2e1 + qJDD(3) * (Ifges(4,5) * t262 + Ifges(4,6) * t266) + m(5) * (t126 * t28 + t127 * t27 + t76 * t88 + t77 * t87 + (t154 * t354 + t262 * t99) * t239) - t287 * (qJD(3) * t277 - t303 * t352) / 0.2e1 + t234 * (qJD(3) * t276 - t300 * t352) / 0.2e1 + t266 * t239 * t170 + t47 * t34 + t48 * t32 + t52 * t31 + qJD(3) * t289 + t328 * t354 + (qJD(3) * t278 - t306 * t352) * t421 + (-Ifges(5,3) * t266 + t262 * t301) * t422 + t22 * t71 + t82 * t16 + t8 * t93 + t10 * t94 + t11 * t95 + t9 * t96 + t506 * t78 + t126 * t84 + t127 * t85 + t131 * t72 - t218 * t333 + t76 * t147 + t77 * t148 + t266 * (Ifges(4,4) * t206 + Ifges(4,2) * t205) / 0.2e1 + t172 * t17 + t99 * t186 + t241 * (-mrSges(4,1) * t205 + mrSges(4,2) * t206); m(3) * qJDD(2) + t402 * t79 + t403 * t78 - t471 * t167 + t472 * t166 + (-m(3) - t505 + t484) * g(3) + (-t16 - t17 + (t147 * t265 - t148 * t261 + t218) * qJD(3) - t468) * t266 + (t170 + (-t261 * t147 - t265 * t148) * qJD(4) + (t71 + t72 + t467) * qJD(3) + t460) * t262 + m(7) * (t166 * t3 - t167 * t2 + t23 * t79 + t24 * t78 - t266 * t7 + t356 * t44) + m(6) * (t116 * t356 - t166 * t6 - t167 * t5 - t25 * t79 + t26 * t78 - t266 * t64) + m(4) * (t104 * t262 + t105 * t266 + (-t164 * t262 + t165 * t266) * qJD(3)) + m(5) * ((-t99 + (-t261 * t87 + t265 * t88) * qJD(3)) * t266 + (qJD(3) * t154 + t273) * t262); t465 * t72 - (-Ifges(4,2) * t358 + t248 + t466) * t357 / 0.2e1 + (t278 * t405 + t307 * t404) * t197 + (-m(5) * t154 + t340 - t467) * t165 + (-pkin(3) * t99 - t110 * t87 - t111 * t88) * m(5) + (t473 / 0.2e1 + mrSges(6,2) * t64 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t7 + Ifges(6,4) * t441 + Ifges(7,5) * t440 + t476 * t423 + t478 * t442) * t199 + t446 * t198 + (t276 * t405 + t301 * t404) * t234 + (-t351 * t87 - t353 * t88 + t459) * mrSges(5,3) + (m(5) * t273 - t147 * t353 - t148 * t351 + t460) * pkin(8) + (t295 * t6 + t134 * t5 - t247 * t64 + (-t46 + t80) * t26 + (-t45 - t81) * t25 + t465 * t116) * m(6) - t472 * t295 + (t121 * t7 - t295 * t3 + t134 * t2 + t470 * t44 + (-t42 + t80) * t24 + (-t43 + t81) * t23) * m(7) + t489 * t157 + t490 * t129 + (g(1) * t250 + g(2) * t249) * (t312 + (t268 * t505 - t341 + t500) * t266 + (m(6) * t247 - m(7) * (-t247 - t299) + t497) * t262) - t485 * t158 + (t328 + t290) * qJD(4) + (Ifges(6,6) * t430 + Ifges(7,6) * t431 + t475 * t419 + t476 * t428 + t495) * t358 - t456 * t358 / 0.2e1 + (t339 - t218) * t164 + t99 * t311 + (-t87 * (mrSges(5,1) * t262 - mrSges(5,3) * t364) - t88 * (-mrSges(5,2) * t262 - mrSges(5,3) * t368) - t289 + t277 * t481 + t288 * t405) * qJD(1) + t471 * t134 - t304 * t280 / 0.2e1 + t402 * t81 + t403 * t80 + t470 * t71 + t303 * t432 + t306 * t433 + t261 * t439 + Ifges(4,3) * qJDD(3) + t265 * t49 / 0.2e1 - t290 * t357 + t109 * t351 / 0.2e1 - t302 * t348 / 0.2e1 + t300 * t422 - pkin(3) * t67 - t42 * t93 - t46 * t94 - t45 * t95 - t43 * t96 - t104 * mrSges(4,2) + t105 * mrSges(4,1) + t121 * t16 - t506 * t128 + (-t262 * t341 - t313 + (-m(7) * t299 - t497) * t266 + t507) * g(3) - t111 * t147 - t110 * t148 + t108 * t334 / 0.2e1 + Ifges(4,6) * t205 + Ifges(4,5) * t206 - t247 * t17; (-(-mrSges(6,1) * t253 - t261 * t445 - t401) * t262 - m(7) * (-pkin(5) * t373 + t225 - t235) + mrSges(7,1) * t373 - t228 + t186) * g(3) + (-t281 - t147) * t87 + (m(6) * t284 - m(7) * (-t284 + t319) + mrSges(5,1) * t159 - mrSges(5,2) * t160 + t462) * g(2) + (-m(7) * (t294 + t318) - m(6) * t294 - mrSges(5,1) * t161 + mrSges(5,2) * t162 + t463) * g(1) + t489 * t275 - t485 * t123 + (-Ifges(5,2) * t197 + t109 - t188) * t481 + t94 * t342 + (t2 * t240 + t246 * t3 - t44 * t55 + (-t30 + t233) * t24 + (t343 - t29) * t23) * m(7) - t154 * (t197 * mrSges(5,1) - mrSges(5,2) * t287) - t234 * (-Ifges(5,5) * t287 - Ifges(5,6) * t197) / 0.2e1 + t274 + (t400 + t148) * t88 + t337 + (m(6) * t25 - t402) * t29 + t402 * t343 + (-m(6) * t26 - t403) * t30 + (t260 * t5 + t264 * t6 + (-t25 * t260 + t26 * t264) * qJD(5)) * t445 + (-m(6) * t116 - t72) * t412 - t197 * (-Ifges(5,1) * t287 - t395) / 0.2e1 + t503 - t504 + t31 * t410 + t33 * t411 + t108 * t421 - t55 * t71 + t233 * t93 + t240 * t34 + t246 * t32; (t123 * t23 + t24 * t275) * mrSges(7,2) + t274 + (Ifges(7,3) * t275 - t391) * t431 + t61 * t427 + (t398 - t402) * t26 + (-t399 - t403) * t25 + (-t228 + (t253 * t454 + t401) * t262) * g(3) + t463 * g(1) + t462 * g(2) - pkin(5) * t32 + qJ(6) * t34 - t70 * t71 + qJD(6) * t93 - t44 * (mrSges(7,1) * t275 + mrSges(7,3) * t123) - t116 * (mrSges(6,1) * t275 - mrSges(6,2) * t123) + (-t123 * t476 + t275 * t474) * t419 + (-Ifges(6,2) * t275 - t120 + t469) * t430 + (-t123 * t478 + t119 - t392 + t58) * t428 + (-pkin(5) * t3 - t318 * g(1) - t319 * g(2) - t225 * g(3) + qJ(6) * t2 - t23 * t26 + t24 * t493 - t44 * t70) * m(7); t275 * t71 - t229 * t93 + (-g(1) * t151 - g(2) * t149 - g(3) * t373 - t24 * t229 + t275 * t44 + t3) * m(7) + t32;];
tau  = t1;
