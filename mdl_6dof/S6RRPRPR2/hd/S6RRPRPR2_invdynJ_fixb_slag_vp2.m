% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:12:05
% EndTime: 2019-03-09 10:12:40
% DurationCPUTime: 20.54s
% Computational Cost: add. (12963->723), mult. (30021->897), div. (0->0), fcn. (22587->14), ass. (0->356)
t522 = -mrSges(6,2) + mrSges(5,1);
t277 = sin(qJ(6));
t281 = cos(qJ(6));
t515 = mrSges(7,1) * t277 + mrSges(7,2) * t281;
t276 = -qJ(3) - pkin(7);
t279 = sin(qJ(2));
t244 = t276 * t279;
t225 = qJD(1) * t244;
t282 = cos(qJ(2));
t246 = t276 * t282;
t226 = qJD(1) * t246;
t397 = cos(pkin(10));
t343 = t397 * t226;
t396 = sin(pkin(10));
t162 = -t225 * t396 + t343;
t300 = t279 * t396 - t282 * t397;
t202 = t300 * qJD(1);
t432 = t202 * pkin(8);
t125 = t162 + t432;
t204 = t396 * t226;
t163 = t397 * t225 + t204;
t299 = -t397 * t279 - t396 * t282;
t203 = t299 * qJD(1);
t431 = t203 * pkin(8);
t126 = t163 + t431;
t278 = sin(qJ(4));
t352 = pkin(2) * t396;
t255 = t278 * t352;
t351 = t397 * pkin(2);
t259 = t351 + pkin(3);
t440 = cos(qJ(4));
t348 = qJD(4) * t440;
t495 = -qJD(4) * t255 - t278 * t125 - t126 * t440 + t259 * t348;
t146 = -t440 * t202 + t203 * t278;
t274 = qJD(2) + qJD(4);
t124 = -t146 * t277 + t274 * t281;
t138 = Ifges(5,4) * t146;
t310 = -t278 * t202 - t203 * t440;
t139 = qJD(6) + t310;
t498 = t139 * Ifges(7,3);
t123 = -t146 * t281 - t274 * t277;
t499 = t123 * Ifges(7,6);
t500 = Ifges(5,1) * t310 + t274 * Ifges(5,5) + t124 * Ifges(7,5) + t138 + t498 + t499;
t211 = qJD(2) * pkin(2) + t225;
t153 = t397 * t211 + t204;
t117 = qJD(2) * pkin(3) + t153 + t431;
t154 = t396 * t211 - t343;
t120 = t154 - t432;
t70 = -t440 * t117 + t120 * t278;
t63 = -pkin(4) * t274 + qJD(5) + t70;
t137 = Ifges(6,6) * t146;
t97 = t274 * Ifges(6,4) - Ifges(6,2) * t310 - t137;
t521 = -mrSges(6,1) * t63 - t70 * mrSges(5,3) + t97 / 0.2e1 - t500 / 0.2e1;
t520 = -t202 / 0.2e1;
t433 = pkin(5) * t146;
t518 = qJD(2) / 0.2e1;
t503 = Ifges(6,4) - Ifges(5,5);
t502 = Ifges(6,5) - Ifges(5,6);
t494 = -qJD(5) - t495;
t395 = qJ(5) * t146;
t245 = -mrSges(3,1) * t282 + mrSges(3,2) * t279;
t275 = qJ(2) + pkin(10);
t268 = sin(t275);
t269 = cos(t275);
t517 = -mrSges(4,1) * t269 + mrSges(4,2) * t268 + t245;
t326 = mrSges(7,1) * t281 - mrSges(7,2) * t277;
t71 = t278 * t117 + t120 * t440;
t64 = -t274 * qJ(5) - t71;
t46 = -t64 + t433;
t412 = Ifges(7,4) * t124;
t56 = Ifges(7,2) * t123 + Ifges(7,6) * t139 + t412;
t516 = t46 * t326 - t281 * t56 / 0.2e1;
t280 = sin(qJ(1));
t283 = cos(qJ(1));
t482 = g(1) * t283 + g(2) * t280;
t270 = qJ(4) + t275;
t257 = sin(t270);
t258 = cos(t270);
t514 = -t522 * t258 + (mrSges(5,2) - mrSges(6,3)) * t257;
t504 = pkin(5) * t310;
t313 = t70 + t504;
t513 = t313 + qJD(5);
t457 = pkin(4) + pkin(9);
t334 = -m(7) * t457 - mrSges(7,3);
t512 = t257 * t334 + t258 * t515;
t370 = qJD(4) * t278;
t367 = qJD(1) * qJD(2);
t228 = qJDD(1) * t279 + t282 * t367;
t215 = t228 * pkin(7);
t372 = qJD(3) * t279;
t151 = qJDD(2) * pkin(2) - qJ(3) * t228 - qJD(1) * t372 - t215;
t347 = t279 * t367;
t366 = qJDD(1) * t282;
t227 = -t347 + t366;
t261 = pkin(7) * t366;
t373 = qJD(2) * t279;
t360 = pkin(7) * t373;
t371 = qJD(3) * t282;
t160 = qJ(3) * t227 + t261 + (-t360 + t371) * qJD(1);
t103 = t397 * t151 - t160 * t396;
t164 = t227 * t396 + t228 * t397;
t75 = qJDD(2) * pkin(3) - t164 * pkin(8) + t103;
t104 = t396 * t151 + t397 * t160;
t301 = -t227 * t397 + t228 * t396;
t86 = -pkin(8) * t301 + t104;
t18 = t117 * t348 - t120 * t370 + t278 * t75 + t440 * t86;
t272 = qJDD(2) + qJDD(4);
t84 = qJD(4) * t310 + t278 * t164 + t301 * t440;
t510 = m(5) * t18 - mrSges(5,2) * t272 - mrSges(5,3) * t84;
t443 = -t274 / 0.2e1;
t447 = t310 / 0.2e1;
t448 = -t310 / 0.2e1;
t452 = -t139 / 0.2e1;
t454 = -t124 / 0.2e1;
t455 = -t123 / 0.2e1;
t271 = t282 * pkin(2);
t260 = t271 + pkin(1);
t231 = -qJD(1) * t260 + qJD(3);
t165 = pkin(3) * t202 + t231;
t45 = -t274 * t457 + t513;
t290 = -qJ(5) * t310 + t165;
t51 = -t146 * t457 + t290;
t21 = -t277 * t51 + t281 * t45;
t22 = t277 * t45 + t281 * t51;
t85 = -pkin(4) * t146 + t290;
t467 = t21 * mrSges(7,1) + t165 * mrSges(5,2) - t22 * mrSges(7,2) - t85 * mrSges(6,3);
t509 = -Ifges(5,1) * t448 - Ifges(7,5) * t454 + Ifges(6,2) * t447 - Ifges(7,6) * t455 - Ifges(7,3) * t452 + t503 * t443 + t467;
t408 = t310 * Ifges(6,6);
t96 = t274 * Ifges(6,5) - Ifges(6,3) * t146 - t408;
t409 = t310 * Ifges(5,4);
t98 = Ifges(5,2) * t146 + t274 * Ifges(5,6) + t409;
t508 = t71 * mrSges(5,3) - t64 * mrSges(6,1) - t96 / 0.2e1 + t98 / 0.2e1;
t43 = qJD(6) * t123 + t272 * t281 + t277 * t84;
t462 = t43 / 0.2e1;
t44 = -qJD(6) * t124 - t272 * t277 + t281 * t84;
t461 = t44 / 0.2e1;
t83 = -t440 * t164 + t202 * t348 - t203 * t370 + t278 * t301;
t82 = qJDD(6) - t83;
t460 = t82 / 0.2e1;
t507 = -t164 / 0.2e1;
t506 = t227 / 0.2e1;
t444 = t272 / 0.2e1;
t505 = pkin(4) * t310;
t16 = -mrSges(7,1) * t44 + mrSges(7,2) * t43;
t67 = mrSges(6,1) * t84 - mrSges(6,3) * t272;
t501 = t16 - t67;
t497 = t231 * mrSges(4,1);
t399 = qJDD(2) / 0.2e1;
t496 = t504 - t494;
t423 = mrSges(6,1) * t146;
t129 = -mrSges(6,3) * t274 - t423;
t90 = -mrSges(7,1) * t123 + mrSges(7,2) * t124;
t493 = t90 - t129;
t492 = qJD(2) * mrSges(4,2);
t427 = -qJD(2) / 0.2e1;
t491 = t203 * t427;
t490 = t257 * t482;
t489 = t310 * t457;
t418 = mrSges(5,3) * t310;
t422 = mrSges(6,1) * t310;
t488 = t274 * t522 - t418 - t422;
t384 = t258 * t283;
t386 = t257 * t283;
t487 = pkin(4) * t384 + qJ(5) * t386;
t214 = -pkin(7) * t347 + t261;
t484 = t214 * t282 + t215 * t279;
t91 = -mrSges(7,2) * t139 + mrSges(7,3) * t123;
t92 = mrSges(7,1) * t139 - mrSges(7,3) * t124;
t316 = -t277 * t91 - t281 * t92;
t28 = mrSges(7,1) * t82 - mrSges(7,3) * t43;
t29 = -mrSges(7,2) * t82 + mrSges(7,3) * t44;
t483 = t277 * t29 + t281 * t28;
t413 = Ifges(4,4) * t203;
t481 = Ifges(4,2) * t520 + Ifges(4,6) * t518 - t413 / 0.2e1 + t154 * mrSges(4,3);
t480 = 0.2e1 * t444;
t479 = 0.2e1 * t399;
t477 = -t258 * mrSges(7,3) - t257 * t515 + t514;
t419 = mrSges(5,3) * t146;
t127 = -mrSges(5,2) * t274 + t419;
t476 = -m(6) * t64 + t127 - t129;
t475 = -m(6) * t63 + t488;
t394 = qJDD(1) * pkin(1);
t188 = -t227 * pkin(2) + qJDD(3) - t394;
t122 = pkin(3) * t301 + t188;
t286 = t83 * qJ(5) - qJD(5) * t310 + t122;
t12 = t457 * t84 + t286;
t19 = -t117 * t370 - t120 * t348 - t278 * t86 + t440 * t75;
t303 = qJDD(5) - t19;
t6 = -t83 * pkin(5) - t272 * t457 + t303;
t1 = qJD(6) * t21 + t12 * t281 + t277 * t6;
t2 = -qJD(6) * t22 - t12 * t277 + t281 * t6;
t474 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t473 = t165 * mrSges(5,1) - t85 * mrSges(6,2);
t472 = -m(3) * pkin(1) - m(4) * t260 - mrSges(2,1) + t514 + t517;
t471 = m(5) * t70 - t475;
t318 = t21 * t277 - t22 * t281;
t289 = -qJD(6) * t318 + t1 * t277 + t2 * t281;
t368 = qJD(6) * t281;
t369 = qJD(6) * t277;
t470 = m(7) * t289 + t91 * t368 - t92 * t369 + t483;
t273 = -pkin(8) + t276;
t469 = -m(3) * pkin(7) + m(4) * t276 - m(7) * (pkin(5) - t273) - mrSges(6,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t10 = Ifges(7,4) * t43 + Ifges(7,2) * t44 + Ifges(7,6) * t82;
t13 = -qJ(5) * t272 - qJD(5) * t274 - t18;
t15 = -t272 * pkin(4) + t303;
t320 = Ifges(7,5) * t277 + Ifges(7,6) * t281;
t411 = Ifges(7,4) * t277;
t322 = Ifges(7,2) * t281 + t411;
t410 = Ifges(7,4) * t281;
t324 = Ifges(7,1) * t277 + t410;
t345 = -t369 / 0.2e1;
t385 = t258 * t280;
t387 = t257 * t280;
t416 = mrSges(7,3) * t281;
t417 = mrSges(7,3) * t277;
t441 = -t277 / 0.2e1;
t463 = Ifges(7,1) * t462 + Ifges(7,4) * t461 + Ifges(7,5) * t460;
t121 = Ifges(7,4) * t123;
t57 = Ifges(7,1) * t124 + Ifges(7,5) * t139 + t121;
t7 = -pkin(5) * t84 - t13;
t468 = t7 * t515 - t2 * t416 - t1 * t417 + t19 * mrSges(5,1) - t18 * mrSges(5,2) - t13 * mrSges(6,3) + t15 * mrSges(6,2) + t57 * t345 + t10 * t441 + (Ifges(7,5) * t281 - Ifges(7,6) * t277) * t460 + (-Ifges(7,2) * t277 + t410) * t461 + (Ifges(7,1) * t281 - t411) * t462 + t281 * t463 + t502 * t84 + t503 * t83 + (Ifges(6,1) + Ifges(5,3)) * t272 + t516 * qJD(6) + (t21 * t369 - t22 * t368) * mrSges(7,3) - (t123 * t322 + t124 * t324 + t139 * t320) * qJD(6) / 0.2e1 + (-mrSges(6,2) * t386 - mrSges(6,3) * t384 - t283 * t512) * g(1) + (-mrSges(6,2) * t387 - mrSges(6,3) * t385 - t280 * t512) * g(2);
t450 = -t146 / 0.2e1;
t451 = t146 / 0.2e1;
t465 = -Ifges(5,2) * t450 + Ifges(6,3) * t451 + t21 * t417 - t22 * t416 + t320 * t452 + t322 * t455 + t324 * t454 + t57 * t441 + t443 * t502 - t473 + t516;
t458 = m(6) + m(7);
t453 = t124 / 0.2e1;
t445 = t202 / 0.2e1;
t442 = t274 / 0.2e1;
t438 = pkin(2) * t279;
t437 = pkin(7) * t282;
t434 = g(3) * t258;
t254 = t258 * pkin(4);
t415 = Ifges(3,4) * t279;
t414 = Ifges(3,4) * t282;
t407 = t153 * mrSges(4,3);
t293 = qJD(2) * t299;
t294 = qJD(2) * t300;
t295 = t278 * t300;
t106 = -qJD(4) * t295 - t278 * t294 - t293 * t440 - t299 * t348;
t393 = t106 * t277;
t392 = t106 * t281;
t292 = t440 * t300;
t158 = -t278 * t299 + t292;
t391 = t158 * t277;
t390 = t158 * t281;
t247 = t257 * qJ(5);
t383 = t277 * t280;
t382 = t277 * t283;
t381 = t280 * t281;
t380 = t281 * t283;
t344 = qJD(2) * t276;
t199 = t279 * t344 + t371;
t200 = t282 * t344 - t372;
t134 = t397 * t199 + t396 * t200;
t167 = t396 * t244 - t397 * t246;
t377 = t254 + t247;
t376 = pkin(3) * t269 + t271;
t375 = qJD(1) * t279;
t374 = qJD(1) * t282;
t365 = qJD(2) ^ 2 / 0.2e1;
t364 = Ifges(7,5) * t43 + Ifges(7,6) * t44 + Ifges(7,3) * t82;
t363 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t375) * t437;
t264 = pkin(2) * t373;
t263 = pkin(2) * t375;
t68 = -t83 * mrSges(6,1) + t272 * mrSges(6,2);
t168 = -pkin(3) * t203 + t263;
t224 = pkin(1) + t376;
t337 = -t224 - t247;
t87 = -t440 * t125 + t126 * t278;
t166 = t397 * t244 + t246 * t396;
t135 = pkin(8) * t299 + t166;
t136 = -pkin(8) * t300 + t167;
t94 = -t440 * t135 + t136 * t278;
t208 = t283 * t224;
t336 = -t273 * t280 + t208;
t335 = t376 + t377;
t232 = qJ(5) * t385;
t333 = -pkin(4) * t387 + t232;
t234 = qJ(5) * t384;
t332 = -pkin(4) * t386 + t234;
t197 = t259 * t440 - t255;
t330 = mrSges(3,1) * t279 + mrSges(3,2) * t282;
t327 = mrSges(5,1) * t257 + mrSges(5,2) * t258;
t323 = t282 * Ifges(3,2) + t415;
t321 = Ifges(3,5) * t282 - Ifges(3,6) * t279;
t319 = t21 * t281 + t22 * t277;
t159 = -t299 * t440 - t295;
t174 = pkin(3) * t300 - t260;
t288 = -t159 * qJ(5) + t174;
t58 = t158 * t457 + t288;
t59 = pkin(5) * t159 + t94;
t34 = t277 * t59 + t281 * t58;
t33 = -t277 * t58 + t281 * t59;
t317 = -t277 * t92 + t281 * t91;
t133 = -t199 * t396 + t397 * t200;
t194 = -pkin(4) - t197;
t314 = t168 - t395;
t312 = pkin(1) * t330;
t95 = t278 * t135 + t136 * t440;
t309 = t158 * t368 + t393;
t308 = t158 * t369 - t392;
t307 = t279 * (Ifges(3,1) * t282 - t415);
t114 = pkin(8) * t293 + t134;
t287 = pkin(8) * t294 + t133;
t35 = -t440 * t114 - t135 * t348 + t136 * t370 - t278 * t287;
t198 = t278 * t259 + t352 * t440;
t296 = Ifges(4,4) * t300;
t291 = mrSges(4,1) * t301 + t164 * mrSges(4,2);
t169 = -pkin(3) * t293 + t264;
t36 = qJD(4) * t95 + t278 * t114 - t440 * t287;
t105 = t274 * t292 - t278 * t293 - t299 * t370;
t37 = t106 * pkin(4) + t105 * qJ(5) - t159 * qJD(5) + t169;
t262 = Ifges(3,4) * t374;
t253 = t258 * pkin(9);
t243 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t374;
t229 = -pkin(3) * t268 - t438;
t213 = t283 * t229;
t212 = t280 * t229;
t210 = Ifges(3,1) * t375 + Ifges(3,5) * qJD(2) + t262;
t209 = Ifges(3,6) * qJD(2) + qJD(1) * t323;
t196 = Ifges(4,4) * t202;
t193 = qJ(5) + t198;
t192 = -t257 * t383 + t380;
t191 = t257 * t381 + t382;
t190 = t257 * t382 + t381;
t189 = t257 * t380 - t383;
t171 = qJD(2) * mrSges(4,1) + mrSges(4,3) * t203;
t170 = -mrSges(4,3) * t202 - t492;
t152 = mrSges(4,1) * t202 - mrSges(4,2) * t203;
t149 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t164;
t148 = -qJDD(2) * mrSges(4,2) - mrSges(4,3) * t301;
t143 = -t203 * Ifges(4,1) + Ifges(4,5) * qJD(2) - t196;
t102 = mrSges(6,2) * t146 - mrSges(6,3) * t310;
t101 = -mrSges(5,1) * t146 + mrSges(5,2) * t310;
t100 = -t395 + t505;
t93 = t158 * pkin(4) + t288;
t89 = t314 + t505;
t77 = t83 * mrSges(5,2);
t76 = t83 * mrSges(6,3);
t65 = mrSges(5,1) * t272 + mrSges(5,3) * t83;
t62 = -t395 + t489;
t60 = -t158 * pkin(5) + t95;
t54 = t314 + t489;
t52 = t87 + t433;
t48 = t71 + t433;
t32 = t277 * t48 + t281 * t62;
t31 = -t277 * t62 + t281 * t48;
t30 = t106 * pkin(9) + t37;
t27 = -t105 * pkin(5) + t36;
t26 = -pkin(5) * t106 - t35;
t25 = t277 * t52 + t281 * t54;
t24 = -t277 * t54 + t281 * t52;
t23 = t84 * pkin(4) + t286;
t4 = -qJD(6) * t34 + t27 * t281 - t277 * t30;
t3 = qJD(6) * t33 + t27 * t277 + t281 * t30;
t5 = [(Ifges(7,1) * t309 - Ifges(7,4) * t308) * t453 + (-m(6) * t13 + t510 - t67) * t95 + (t1 * t390 - t2 * t391 - t21 * t309 - t22 * t308) * mrSges(7,3) - qJDD(2) * mrSges(3,2) * t437 + (-t299 * t497 - t363) * qJD(2) + (Ifges(3,1) * t228 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t228) + Ifges(3,4) * t506 + t479 * Ifges(3,5)) * t279 + (t188 * mrSges(4,1) - t104 * mrSges(4,3) - Ifges(4,1) * t491 + Ifges(4,4) * t507 - Ifges(4,5) * t365 + Ifges(4,2) * t301 - Ifges(4,6) * t479 - t231 * t492) * t300 + Ifges(3,6) * t282 * t399 + (-Ifges(5,4) * t447 - Ifges(5,2) * t451 + Ifges(6,6) * t448 + Ifges(6,3) * t450 + t502 * t442 + t473 - t508) * t106 + (-m(7) * (pkin(9) * t384 + t208 + t487) - t190 * mrSges(7,1) - t189 * mrSges(7,2) - mrSges(7,3) * t384 - m(6) * (t336 + t487) - m(5) * t336 + t472 * t283 + t469 * t280) * g(2) + (t227 * t437 + t484) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t484) + (-t202 * (Ifges(4,2) * t299 - t296) + t282 * t210) * t518 + (t122 * mrSges(5,1) + t13 * mrSges(6,1) - t23 * mrSges(6,2) - t18 * mrSges(5,3) + t320 * t460 + t322 * t461 + t324 * t462 - t7 * t326 + t56 * t345 + (Ifges(5,2) + Ifges(6,3)) * t84 + (Ifges(5,4) + Ifges(6,6)) * t83 + t502 * t480) * t158 + t471 * t36 + t228 * t414 / 0.2e1 - (t188 * mrSges(4,2) - t103 * mrSges(4,3) + Ifges(4,1) * t164 - Ifges(4,4) * t301 + Ifges(4,5) * t479) * t299 + m(5) * (t122 * t174 + t165 * t169) + m(6) * (t23 * t93 + t37 * t85) + (t15 * mrSges(6,1) + t122 * mrSges(5,2) - t19 * mrSges(5,3) - t23 * mrSges(6,3) + Ifges(5,5) * t444 + Ifges(7,5) * t462 + Ifges(7,6) * t461 + Ifges(7,3) * t460 + (-Ifges(5,4) / 0.2e1 - Ifges(6,6)) * t84 - Ifges(6,2) * t83 - t480 * Ifges(6,4) + t474) * t159 + (-Ifges(5,4) * t84 + Ifges(5,5) * t272 + t364) * t159 / 0.2e1 + (-Ifges(5,1) * t447 + Ifges(6,2) * t448 + Ifges(6,6) * t450 - Ifges(5,4) * t451 - Ifges(7,5) * t453 - t498 / 0.2e1 - t499 / 0.2e1 + t503 * t442 - t467 + t521) * t105 - t245 * t394 + t56 * t392 / 0.2e1 + t57 * t393 / 0.2e1 + t60 * t16 + (-t192 * mrSges(7,1) + t191 * mrSges(7,2) + ((m(5) + m(6)) * t273 + t469) * t283 + (-m(6) * (t337 - t254) - m(7) * t337 - t258 * t334 + m(5) * t224 - t472) * t280) * g(1) + t323 * t506 + t296 * t507 + t123 * (Ifges(7,4) * t309 - Ifges(7,2) * t308) / 0.2e1 + (-m(5) * t71 - t476) * t35 + t174 * (t84 * mrSges(5,1) - t77) + t139 * (Ifges(7,5) * t309 - Ifges(7,6) * t308) / 0.2e1 + m(7) * (t1 * t34 + t2 * t33 + t21 * t4 + t22 * t3 + t26 * t46 + t60 * t7) + (Ifges(4,6) * t299 + t321) * t365 + (-m(5) * t19 + m(6) * t15 - t65 + t68) * t94 - t83 * Ifges(5,1) * t159 + t481 * t293 - (-t407 + t143 / 0.2e1) * t294 + (t282 * (-Ifges(3,2) * t279 + t414) + t307) * t367 / 0.2e1 + (qJD(6) * t57 + t10) * t390 / 0.2e1 + t152 * t264 - t209 * t373 / 0.2e1 - t312 * t367 + t33 * t28 + t34 * t29 + m(4) * (t103 * t166 + t104 * t167 + t133 * t153 + t134 * t154 - t188 * t260 + t231 * t264) - t243 * t360 + Ifges(2,3) * qJDD(1) - t260 * t291 + t282 * (Ifges(3,4) * t228 + Ifges(3,2) * t227 + Ifges(3,6) * qJDD(2)) / 0.2e1 + t26 * t90 + t3 * t91 + t4 * t92 + t93 * (-t84 * mrSges(6,2) + t76) + t37 * t102 + Ifges(4,4) * t299 * t491 + t391 * t463 + t166 * t149 + t167 * t148 + t169 * t101 + t134 * t170 + t133 * t171 + t46 * (mrSges(7,1) * t308 + mrSges(7,2) * t309) - pkin(1) * (-mrSges(3,1) * t227 + mrSges(3,2) * t228); ((m(7) * t319 - t316 + t471) * qJD(4) + t510) * t198 + (-t63 * t87 - t85 * t89 - g(1) * (t213 + t332) - g(2) * (t212 + t333) - t13 * t193 + t15 * t194 + t494 * t64) * m(6) + t494 * t129 + t495 * t127 + (-t21 * t24 - t22 * t25 - g(1) * (t213 + t234) - g(2) * (t212 + t232) + t193 * t7 + t496 * t46) * m(7) + t496 * t90 + t482 * (m(4) * t438 - m(5) * t229 + mrSges(4,1) * t268 + mrSges(4,2) * t269 + t327 + t330) + (t209 / 0.2e1 + pkin(7) * t243) * t375 + t488 * t87 + (-Ifges(5,4) * t448 + Ifges(6,6) * t447 + t465 + t508) * t310 + t470 * (-pkin(9) + t194) + (-t196 + t143) * t445 - (-mrSges(4,2) * t231 + Ifges(4,5) * t427 + t407) * t202 + t501 * t193 + (Ifges(5,4) * t450 - Ifges(6,6) * t451 - t509 + t521) * t146 + t149 * t351 + t148 * t352 - Ifges(4,6) * t301 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t468 + ((t103 * t397 + t104 * t396) * pkin(2) - t153 * t162 - t154 * t163 - t231 * t263) * m(4) + (-m(4) * t271 - m(5) * t376 - m(7) * (t253 + t335) - m(6) * t335 + t477 + t517) * g(3) + (t363 + (t312 - t307 / 0.2e1) * qJD(1)) * qJD(1) - (-Ifges(3,2) * t375 + t210 + t262) * t374 / 0.2e1 + (-t165 * t168 + t19 * t197 + t495 * t71 - t70 * t87) * m(5) - t321 * t367 / 0.2e1 - t152 * t263 + (Ifges(4,1) * t520 + t413 / 0.2e1 + Ifges(4,6) * t427 + Ifges(4,2) * t445 + t497 - t481) * t203 - t25 * t91 - t24 * t92 - t89 * t102 + t103 * mrSges(4,1) - t104 * mrSges(4,2) + Ifges(4,5) * t164 - t168 * t101 - t163 * t170 - t162 * t171 + t194 * t68 + t197 * t65 - t214 * mrSges(3,2) - t215 * mrSges(3,1) + Ifges(3,6) * t227 + Ifges(3,5) * t228; -t77 + t76 + t522 * t84 + (-t127 - t493) * t146 + (t316 + t488) * t310 + t316 * qJD(6) + t291 + t281 * t29 + t202 * t170 - t203 * t171 - t277 * t28 + (-g(1) * t280 + g(2) * t283) * (m(4) + m(5) + t458) + (t1 * t281 - t139 * t319 - t146 * t46 - t2 * t277) * m(7) + (t146 * t64 - t310 * t63 + t23) * m(6) + (-t146 * t71 - t310 * t70 + t122) * m(5) + (-t153 * t203 + t154 * t202 + t188) * m(4); t493 * qJD(5) + t501 * qJ(5) + (t408 + t98) * t447 + (-t409 + t96) * t448 + (-t137 + t97) * t451 + (t138 + t500) * t450 + t465 * t310 - t509 * t146 + (-pkin(4) * t15 - g(1) * t332 - g(2) * t333 - qJ(5) * t13 - qJD(5) * t64 - t100 * t85) * m(6) - pkin(4) * t68 + t468 - t64 * t422 + t482 * t327 + (-m(7) * (t253 + t377) - m(6) * t377 + t477) * g(3) - t470 * t457 + (-g(1) * t234 - g(2) * t232 + t7 * qJ(5) - t21 * t31 - t22 * t32 + t513 * t46) * m(7) + (t418 + t475) * t71 + (-t419 + t476) * t70 + t313 * t90 - t32 * t91 - t31 * t92 - t100 * t102 - t63 * t423; -t493 * t274 + t317 * qJD(6) + t458 * t434 + (t102 + t317) * t310 + t68 + (-t274 * t46 - t310 * t318 + t289 - t490) * m(7) + (t274 * t64 + t310 * t85 + t15 - t490) * m(6) + t483; -t46 * (mrSges(7,1) * t124 + mrSges(7,2) * t123) + (Ifges(7,1) * t123 - t412) * t454 + t56 * t453 + (Ifges(7,5) * t123 - Ifges(7,6) * t124) * t452 - t21 * t91 + t22 * t92 - g(1) * (mrSges(7,1) * t189 - mrSges(7,2) * t190) - g(2) * (mrSges(7,1) * t191 + mrSges(7,2) * t192) + t326 * t434 + (t123 * t21 + t124 * t22) * mrSges(7,3) + t364 + (-Ifges(7,2) * t124 + t121 + t57) * t455 + t474;];
tau  = t5;
