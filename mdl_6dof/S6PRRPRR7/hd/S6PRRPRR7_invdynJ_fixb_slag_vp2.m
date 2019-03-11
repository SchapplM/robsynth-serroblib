% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:17
% EndTime: 2019-03-08 22:30:54
% DurationCPUTime: 24.63s
% Computational Cost: add. (7577->756), mult. (16320->1032), div. (0->0), fcn. (11425->14), ass. (0->356)
t533 = -mrSges(4,1) + mrSges(5,2);
t260 = sin(qJ(3));
t264 = cos(qJ(3));
t300 = pkin(9) * t260 - qJ(4) * t264;
t377 = qJD(3) * t260;
t321 = pkin(3) * t377 - qJD(4) * t260;
t138 = qJD(3) * t300 + t321;
t400 = qJ(4) * t260;
t336 = -pkin(2) - t400;
t449 = pkin(3) + pkin(9);
t182 = -t264 * t449 + t336;
t375 = qJD(3) * t264;
t448 = pkin(4) + pkin(8);
t198 = t448 * t375;
t218 = t448 * t260;
t259 = sin(qJ(5));
t263 = cos(qJ(5));
t261 = sin(qJ(2));
t265 = cos(qJ(2));
t389 = t260 * t265;
t290 = t259 * t389 + t261 * t263;
t373 = qJD(5) * t263;
t374 = qJD(5) * t259;
t256 = sin(pkin(6));
t382 = qJD(1) * t256;
t497 = t263 * t138 - t182 * t374 + t259 * t198 + t218 * t373 - t290 * t382;
t385 = t263 * t265;
t291 = -t259 * t261 + t260 * t385;
t532 = -t138 * t259 + t263 * t198 - t291 * t382;
t531 = -m(5) - m(6);
t262 = cos(qJ(6));
t246 = t260 * qJD(2);
t235 = t246 + qJD(5);
t369 = t264 * qJD(2);
t376 = qJD(3) * t263;
t284 = t259 * t369 - t376;
t349 = t265 * t382;
t122 = qJD(2) * t182 - t349;
t350 = t261 * t382;
t201 = qJD(2) * pkin(8) + t350;
t257 = cos(pkin(6));
t381 = qJD(1) * t257;
t229 = t264 * t381;
t118 = -t260 * (pkin(4) * qJD(2) + t201) + t229;
t478 = -t118 + qJD(4);
t90 = -qJD(3) * t449 + t478;
t49 = -t122 * t259 + t263 * t90;
t35 = pkin(10) * t284 + t49;
t29 = pkin(5) * t235 + t35;
t258 = sin(qJ(6));
t188 = -qJD(3) * t259 - t263 * t369;
t50 = t122 * t263 + t259 * t90;
t36 = pkin(10) * t188 + t50;
t413 = t258 * t36;
t15 = t262 * t29 - t413;
t368 = qJD(2) * qJD(3);
t200 = qJDD(2) * t260 + t264 * t368;
t366 = qJDD(1) * t257;
t380 = qJD(2) * t256;
t343 = qJD(1) * t380;
t221 = t265 * t343;
t367 = qJDD(1) * t256;
t156 = t261 * t367 + t221;
t378 = qJD(3) * t257;
t518 = qJDD(2) * pkin(8) + qJD(1) * t378 + t156;
t69 = -t201 * t375 - t260 * t518 + t264 * t366;
t277 = qJDD(4) - t69;
t42 = pkin(4) * t200 - qJDD(3) * t449 + t277;
t199 = -t264 * qJDD(2) + t260 * t368;
t220 = t261 * t343;
t155 = t265 * t367 - t220;
t140 = -qJDD(2) * pkin(2) - t155;
t270 = -qJ(4) * t200 - qJD(4) * t246 + t140;
t54 = t199 * t449 + t270;
t12 = -qJD(5) * t50 - t259 * t54 + t263 * t42;
t187 = qJDD(5) + t200;
t91 = qJD(5) * t188 + qJDD(3) * t263 + t199 * t259;
t6 = pkin(5) * t187 - pkin(10) * t91 + t12;
t11 = -t122 * t374 + t259 * t42 + t263 * t54 + t90 * t373;
t92 = qJD(5) * t284 - qJDD(3) * t259 + t199 * t263;
t7 = pkin(10) * t92 + t11;
t2 = qJD(6) * t15 + t258 * t6 + t262 * t7;
t530 = t2 * mrSges(7,2);
t408 = t262 * t36;
t16 = t258 * t29 + t408;
t3 = -qJD(6) * t16 - t258 * t7 + t262 * t6;
t529 = t3 * mrSges(7,1);
t528 = t11 * mrSges(6,2);
t527 = t12 * mrSges(6,1);
t507 = Ifges(4,5) - Ifges(5,4);
t506 = Ifges(5,5) - Ifges(4,6);
t192 = t259 * t218;
t392 = t259 * t260;
t293 = pkin(5) * t264 - pkin(10) * t392;
t333 = pkin(10) * t264 - t182;
t526 = t293 * qJD(3) + (t263 * t333 - t192) * qJD(5) + t532;
t372 = qJD(5) * t264;
t345 = t259 * t372;
t280 = t260 * t376 + t345;
t525 = -pkin(10) * t280 - t497;
t430 = pkin(10) + t449;
t136 = t264 * t201 + t260 * t381;
t119 = pkin(4) * t369 + t136;
t242 = pkin(3) * t246;
t149 = qJD(2) * t300 + t242;
t73 = t263 * t119 - t149 * t259;
t524 = -qJD(2) * t293 + t430 * t374 - t73;
t204 = t430 * t263;
t346 = t263 * t246;
t74 = t259 * t119 + t263 * t149;
t523 = pkin(10) * t346 + qJD(5) * t204 + t74;
t363 = qJD(5) + qJD(6);
t370 = qJD(6) * t258;
t388 = t262 * t263;
t107 = -t258 * t374 - t259 * t370 + t363 * t388;
t394 = t258 * t259;
t144 = -t246 * t394 + t262 * t346;
t493 = t144 + t107;
t295 = t258 * t263 + t262 * t259;
t108 = t363 * t295;
t283 = t295 * t260;
t145 = qJD(2) * t283;
t492 = -t145 - t108;
t355 = mrSges(4,3) * t246;
t357 = mrSges(5,1) * t246;
t489 = qJD(3) * t533 + t355 + t357;
t426 = mrSges(6,3) * t188;
t133 = -mrSges(6,2) * t235 + t426;
t425 = mrSges(6,3) * t284;
t134 = mrSges(6,1) * t235 + t425;
t296 = t263 * t133 - t259 * t134;
t75 = mrSges(6,1) * t187 - mrSges(6,3) * t91;
t76 = -mrSges(6,2) * t187 + mrSges(6,3) * t92;
t522 = t296 * qJD(5) + t259 * t76 + t263 * t75;
t322 = t262 * t188 + t258 * t284;
t180 = qJDD(6) + t187;
t26 = qJD(6) * t322 + t258 * t92 + t262 * t91;
t99 = t188 * t258 - t262 * t284;
t27 = -qJD(6) * t99 - t258 * t91 + t262 * t92;
t362 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t180;
t440 = Ifges(7,4) * t99;
t224 = t246 + t363;
t443 = -t224 / 0.2e1;
t451 = -t99 / 0.2e1;
t254 = qJD(3) * qJ(4);
t100 = t254 + t119;
t77 = -pkin(5) * t188 + t100;
t521 = t529 - t530 + t362 + (Ifges(7,5) * t322 - Ifges(7,6) * t99) * t443 + (t15 * t322 + t16 * t99) * mrSges(7,3) - t77 * (mrSges(7,1) * t99 + mrSges(7,2) * t322) + (Ifges(7,1) * t322 - t440) * t451;
t418 = Ifges(5,6) * t264;
t302 = -t260 * Ifges(5,2) - t418;
t520 = t16 * mrSges(7,2) + Ifges(5,4) * qJD(3) / 0.2e1 + qJD(2) * t302 / 0.2e1 - t15 * mrSges(7,1);
t477 = m(7) - t531;
t519 = -m(4) - t477;
t312 = mrSges(5,2) * t264 - mrSges(5,3) * t260;
t255 = qJ(5) + qJ(6);
t247 = sin(t255);
t248 = cos(t255);
t313 = t247 * mrSges(7,1) + t248 * mrSges(7,2);
t317 = mrSges(4,1) * t264 - mrSges(4,2) * t260;
t438 = pkin(5) * t259;
t508 = m(7) * (-pkin(10) - pkin(9));
t517 = -t260 * (m(7) * t438 + t313) + t264 * (-mrSges(7,3) + t508) + t312 - mrSges(3,1) - t317;
t401 = sin(pkin(11));
t327 = t401 * t265;
t402 = cos(pkin(11));
t330 = t402 * t261;
t165 = t257 * t330 + t327;
t332 = t256 * t402;
t110 = t165 * t260 + t264 * t332;
t328 = t401 * t261;
t329 = t402 * t265;
t167 = -t257 * t328 + t329;
t331 = t256 * t401;
t112 = t167 * t260 - t264 * t331;
t396 = t256 * t261;
t353 = t260 * t396;
t171 = -t257 * t264 + t353;
t282 = -g(1) * t112 - g(2) * t110 - g(3) * t171;
t294 = -t388 + t394;
t516 = -t15 * t492 - t16 * t493 - t2 * t295 + t294 * t3 - t282;
t473 = -m(6) * pkin(9) - mrSges(6,3);
t515 = pkin(3) * t477 - t473 - t508 - t533;
t514 = -t248 * mrSges(7,1) + t247 * mrSges(7,2) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t93 = Ifges(7,4) * t322;
t512 = -Ifges(7,2) * t99 + t93;
t461 = t26 / 0.2e1;
t460 = t27 / 0.2e1;
t447 = t180 / 0.2e1;
t193 = t263 * t218;
t80 = pkin(5) * t260 + t259 * t333 + t193;
t106 = t263 * t182 + t192;
t386 = t263 * t264;
t86 = -pkin(10) * t386 + t106;
t34 = t258 * t80 + t262 * t86;
t510 = -qJD(6) * t34 + t258 * t525 + t262 * t526;
t33 = -t258 * t86 + t262 * t80;
t509 = qJD(6) * t33 + t258 * t526 - t262 * t525;
t203 = t430 * t259;
t120 = t203 * t258 - t204 * t262;
t505 = qJD(6) * t120 + t258 * t524 - t262 * t523;
t121 = -t203 * t262 - t204 * t258;
t504 = -qJD(6) * t121 + t258 * t523 + t262 * t524;
t462 = m(7) * pkin(5);
t499 = -mrSges(6,1) - t462;
t498 = -qJD(5) * t106 + t532;
t240 = pkin(5) * t263 + pkin(4);
t496 = pkin(5) * t373 + qJD(4) - t229 - (-qJD(2) * t240 - t201) * t260;
t164 = -t257 * t329 + t328;
t399 = t164 * t264;
t495 = -pkin(3) * t399 - t164 * t400;
t166 = t257 * t327 + t330;
t398 = t166 * t264;
t494 = -pkin(3) * t398 - t166 * t400;
t160 = t200 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t491 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t200 + t160;
t159 = mrSges(5,1) * t199 - qJDD(3) * mrSges(5,3);
t490 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t199 - t159;
t354 = mrSges(4,3) * t369;
t211 = -qJD(3) * mrSges(4,2) + t354;
t109 = -mrSges(6,1) * t188 - mrSges(6,2) * t284;
t356 = mrSges(5,1) * t369;
t212 = -qJD(3) * mrSges(5,3) - t356;
t487 = -t212 + t109;
t52 = -mrSges(7,1) * t322 + mrSges(7,2) * t99;
t358 = -t52 - t487;
t488 = t211 - t358;
t486 = t212 - t211;
t419 = Ifges(5,6) * t260;
t485 = t260 * (-Ifges(5,2) * t264 + t419) + t264 * (Ifges(5,3) * t260 - t418);
t484 = t260 * t506 + t264 * t507;
t68 = -t201 * t377 + t260 * t366 + t264 * t518;
t482 = -t260 * t69 + t264 * t68;
t59 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4) - t68;
t60 = -qJDD(3) * pkin(3) + t277;
t481 = t260 * t60 - t264 * t59;
t479 = t11 * t259 + t12 * t263;
t241 = Ifges(4,4) * t369;
t475 = Ifges(4,1) * t246 + Ifges(4,5) * qJD(3) - Ifges(6,5) * t284 + t99 * Ifges(7,5) + t188 * Ifges(6,6) + Ifges(7,6) * t322 + t235 * Ifges(6,3) + t224 * Ifges(7,3) + t241;
t301 = -t264 * Ifges(5,3) - t419;
t422 = Ifges(6,4) * t284;
t84 = t188 * Ifges(6,2) + t235 * Ifges(6,6) - t422;
t183 = Ifges(6,4) * t188;
t85 = -Ifges(6,1) * t284 + t235 * Ifges(6,5) + t183;
t472 = Ifges(5,5) * qJD(3) + qJD(2) * t301 + t259 * t85 + t263 * t84;
t470 = t264 * t363;
t316 = t263 * mrSges(6,1) - t259 * mrSges(6,2);
t237 = qJ(4) + t438;
t315 = mrSges(6,1) * t259 + mrSges(6,2) * t263;
t469 = -m(7) * t237 + qJ(4) * t531 + mrSges(4,2) - mrSges(5,3) - t313 - t315;
t390 = t260 * t263;
t467 = mrSges(6,1) * t392 + mrSges(6,2) * t390 - t517;
t431 = pkin(8) + t240;
t465 = -m(6) * t448 - m(7) * t431 - t316 + t514;
t268 = qJD(2) ^ 2;
t464 = Ifges(7,4) * t461 + Ifges(7,2) * t460 + Ifges(7,6) * t447;
t463 = Ifges(7,1) * t461 + Ifges(7,4) * t460 + Ifges(7,5) * t447;
t459 = -t91 * Ifges(6,4) / 0.2e1 - t92 * Ifges(6,2) / 0.2e1 - t187 * Ifges(6,6) / 0.2e1;
t40 = Ifges(7,2) * t322 + Ifges(7,6) * t224 + t440;
t458 = -t40 / 0.2e1;
t41 = Ifges(7,1) * t99 + Ifges(7,5) * t224 + t93;
t457 = -t41 / 0.2e1;
t456 = t41 / 0.2e1;
t455 = t91 / 0.2e1;
t454 = t92 / 0.2e1;
t453 = -t322 / 0.2e1;
t452 = t322 / 0.2e1;
t450 = t99 / 0.2e1;
t446 = t187 / 0.2e1;
t444 = -t284 / 0.2e1;
t442 = t224 / 0.2e1;
t439 = pkin(5) * t284;
t249 = t264 * pkin(8);
t429 = (t110 * t248 - t164 * t247) * mrSges(7,1) + (-t110 * t247 - t164 * t248) * mrSges(7,2);
t428 = (t112 * t248 - t166 * t247) * mrSges(7,1) + (-t112 * t247 - t166 * t248) * mrSges(7,2);
t395 = t256 * t265;
t427 = (t171 * t248 + t247 * t395) * mrSges(7,1) + (-t171 * t247 + t248 * t395) * mrSges(7,2);
t424 = Ifges(4,4) * t260;
t423 = Ifges(4,4) * t264;
t421 = Ifges(6,4) * t259;
t420 = Ifges(6,4) * t263;
t391 = t259 * t264;
t384 = t264 * t265;
t219 = t264 * pkin(4) + t249;
t379 = qJD(2) * t261;
t361 = Ifges(6,5) * t91 + Ifges(6,6) * t92 + Ifges(6,3) * t187;
t153 = t164 * pkin(2);
t352 = -t153 + t495;
t154 = t166 * pkin(2);
t351 = -t154 + t494;
t348 = t256 * t379;
t347 = t265 * t380;
t339 = -t373 / 0.2e1;
t335 = pkin(8) * t165 - t153;
t334 = pkin(8) * t167 - t154;
t326 = -t368 / 0.2e1;
t135 = t201 * t260 - t229;
t311 = Ifges(6,1) * t263 - t421;
t310 = Ifges(6,1) * t259 + t420;
t309 = t264 * Ifges(4,2) + t424;
t307 = -Ifges(6,2) * t259 + t420;
t306 = Ifges(6,2) * t263 + t421;
t304 = Ifges(6,5) * t263 - Ifges(6,6) * t259;
t303 = Ifges(6,5) * t259 + Ifges(6,6) * t263;
t299 = t49 * t259 - t50 * t263;
t116 = t171 * t263 + t259 * t395;
t117 = -t171 * t259 + t256 * t385;
t57 = t116 * t262 + t117 * t258;
t58 = t116 * t258 - t117 * t262;
t205 = -pkin(3) * t264 + t336;
t172 = t257 * t260 + t264 * t396;
t137 = qJD(2) * t205 - t349;
t289 = t137 * (-mrSges(5,2) * t260 - mrSges(5,3) * t264);
t202 = -qJD(2) * pkin(2) - t349;
t288 = t202 * (mrSges(4,1) * t260 + mrSges(4,2) * t264);
t287 = t260 * (Ifges(4,1) * t264 - t424);
t279 = t259 * t377 - t263 * t372;
t275 = Ifges(6,5) * t264 + t260 * t310;
t274 = Ifges(6,6) * t264 + t260 * t306;
t273 = Ifges(6,3) * t264 + t260 * t303;
t43 = -pkin(4) * t199 - t59;
t271 = -qJD(5) * t299 + t479;
t197 = t448 * t377;
t196 = -qJ(4) * t369 + t242;
t195 = t317 * qJD(2);
t194 = t312 * qJD(2);
t173 = Ifges(4,6) * qJD(3) + qJD(2) * t309;
t170 = pkin(5) * t386 + t219;
t162 = -qJ(4) * t375 + t321;
t151 = t295 * t264;
t150 = t294 * t264;
t127 = -t254 - t136;
t126 = -qJD(3) * pkin(3) + qJD(4) + t135;
t125 = -mrSges(5,2) * t199 - mrSges(5,3) * t200;
t124 = mrSges(4,1) * t199 + mrSges(4,2) * t200;
t123 = -pkin(5) * t345 - t377 * t431;
t115 = -qJD(3) * t353 + (t347 + t378) * t264;
t114 = qJD(3) * t172 + t260 * t347;
t105 = -t182 * t259 + t193;
t79 = mrSges(7,1) * t224 - mrSges(7,3) * t99;
t78 = -mrSges(7,2) * t224 + mrSges(7,3) * t322;
t72 = pkin(3) * t199 + t270;
t63 = -t294 * t377 + t295 * t470;
t62 = qJD(3) * t283 + t294 * t470;
t48 = qJD(5) * t116 + t114 * t259 + t263 * t348;
t47 = qJD(5) * t117 + t114 * t263 - t259 * t348;
t44 = -mrSges(6,1) * t92 + mrSges(6,2) * t91;
t32 = t91 * Ifges(6,1) + t92 * Ifges(6,4) + t187 * Ifges(6,5);
t23 = -pkin(5) * t92 + t43;
t22 = -mrSges(7,2) * t180 + mrSges(7,3) * t27;
t21 = mrSges(7,1) * t180 - mrSges(7,3) * t26;
t18 = t262 * t35 - t413;
t17 = -t258 * t35 - t408;
t14 = -qJD(6) * t58 - t258 * t48 + t262 * t47;
t13 = qJD(6) * t57 + t258 * t47 + t262 * t48;
t10 = -mrSges(7,1) * t27 + mrSges(7,2) * t26;
t1 = [m(2) * qJDD(1) + t116 * t75 - t117 * t76 + t13 * t78 + t48 * t133 + t47 * t134 + t14 * t79 + t57 * t21 + t58 * t22 + t491 * t171 + t489 * t114 + (t44 + t10 + t490) * t172 + t488 * t115 + (-m(2) - m(3) + t519) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t268 - t124 - t125) * t265 + (-mrSges(3,1) * t268 - mrSges(3,2) * qJDD(2) + (t194 - t195) * qJD(2)) * t261) * t256 + m(3) * (qJDD(1) * t257 ^ 2 + (t155 * t265 + t156 * t261) * t256) + m(4) * (t114 * t135 + t115 * t136 - t171 * t69 + t172 * t68 + (-t140 * t265 + t202 * t379) * t256) + m(5) * (t114 * t126 - t115 * t127 + t171 * t60 - t172 * t59 + (t137 * t379 - t265 * t72) * t256) + m(7) * (t115 * t77 + t13 * t16 + t14 * t15 + t172 * t23 + t2 * t58 + t3 * t57) + m(6) * (t100 * t115 - t11 * t117 + t116 * t12 + t172 * t43 + t47 * t49 + t48 * t50); (g(1) * t398 + g(2) * t399 - t11 * t386 + t12 * t391 - t279 * t49 + t280 * t50) * mrSges(6,3) + (-t15 * t62 + t150 * t2 + t151 * t3 + t16 * t63) * mrSges(7,3) + t260 * t529 + (Ifges(7,4) * t62 + Ifges(7,2) * t63) * t452 + t481 * mrSges(5,1) + t482 * mrSges(4,3) + t200 * t423 / 0.2e1 + t489 * (pkin(8) * t375 - t260 * t349) + (t221 - t156) * mrSges(3,2) + (t264 * (-Ifges(4,2) * t260 + t423) + t287) * t368 / 0.2e1 + (-Ifges(7,4) * t151 + Ifges(7,2) * t150 + Ifges(7,6) * t260) * t460 + (-Ifges(7,1) * t151 + Ifges(7,4) * t150 + Ifges(7,5) * t260) * t461 + t23 * (-mrSges(7,1) * t150 - mrSges(7,2) * t151) + (-Ifges(7,5) * t151 + Ifges(7,6) * t150 + Ifges(7,3) * t260) * t447 - t260 * t528 + (-t50 * mrSges(6,2) + t49 * mrSges(6,1) + Ifges(7,3) * t442 + Ifges(7,5) * t450 + Ifges(7,6) * t452 + t475 / 0.2e1 + t126 * mrSges(5,1) + t135 * mrSges(4,3) - t520) * t375 + (Ifges(7,5) * t62 + Ifges(7,6) * t63) * t442 + (-t350 + t162) * t194 + t491 * pkin(8) * t260 + (-Ifges(4,4) * t199 + Ifges(4,5) * qJDD(3) + t361 + t362) * t260 / 0.2e1 + t260 * t527 + t43 * t316 * t264 + (-m(7) * t351 - m(5) * (t334 + t494) - m(4) * t334 - m(6) * (-pkin(9) * t398 + t351) + t465 * t167 + t467 * t166) * g(1) + (-m(7) * t352 - m(5) * (t335 + t495) - m(4) * t335 - m(6) * (-pkin(9) * t399 + t352) + t465 * t165 + t467 * t164) * g(2) + (-m(6) * t100 - m(7) * t77 - t488) * t264 * t349 + (t486 * pkin(8) - t173 / 0.2e1 + t472 / 0.2e1 + t127 * mrSges(5,1) - t136 * mrSges(4,3)) * t377 + (t220 + t155) * mrSges(3,1) - t32 * t391 / 0.2e1 + t188 * (qJD(3) * t274 - t307 * t372) / 0.2e1 + t235 * (qJD(3) * t273 - t304 * t372) / 0.2e1 + (t519 * (pkin(2) * t395 + pkin(8) * t396) + (-t477 * (pkin(3) * t384 + qJ(4) * t389) - t290 * mrSges(6,1) - t291 * mrSges(6,2) + t473 * t384 + t517 * t265 + (-m(6) * pkin(4) - m(7) * t240 + t514) * t261) * t256) * g(3) + t195 * t350 + t84 * t345 / 0.2e1 - t260 * t530 + t200 * t260 * Ifges(4,1) + t62 * t456 + t386 * t459 - t151 * t463 + t150 * t464 + (Ifges(7,1) * t62 + Ifges(7,4) * t63) * t450 + t219 * t44 + t205 * t125 - t197 * t109 + t170 * t10 + Ifges(3,3) * qJDD(2) + t123 * t52 - pkin(2) * t124 + t105 * t75 + t106 * t76 + t77 * (-mrSges(7,1) * t63 + mrSges(7,2) * t62) + t63 * t40 / 0.2e1 + t34 * t22 + t33 * t21 - t140 * t317 + t72 * t312 - t199 * t309 / 0.2e1 + t199 * t301 / 0.2e1 - t200 * t302 / 0.2e1 + (qJD(3) * t275 - t311 * t372) * t444 + (Ifges(6,3) * t260 - t264 * t303) * t446 + (Ifges(6,6) * t260 - t264 * t306) * t454 + (Ifges(6,5) * t260 - t264 * t310) * t455 + qJD(3) * t288 + qJD(3) * t289 + (-(t137 * t261 + (t126 * t260 - t127 * t264) * t265) * t382 + t137 * t162 + t205 * t72 + ((t126 * t264 + t127 * t260) * qJD(3) + t481) * pkin(8)) * m(5) + (-(t202 * t261 + (t135 * t260 + t136 * t264) * t265) * t382 - pkin(2) * t140 + ((t135 * t264 - t136 * t260) * qJD(3) + t482) * pkin(8)) * m(4) + t484 * qJD(3) ^ 2 / 0.2e1 + t485 * t326 + t490 * t249 + t497 * t133 + t498 * t134 + (-t100 * t197 + t105 * t12 + t106 * t11 + t219 * t43 + t49 * t498 + t497 * t50) * m(6) - t260 * (Ifges(5,4) * qJDD(3) - Ifges(5,2) * t200 + Ifges(5,6) * t199) / 0.2e1 + t264 * (Ifges(4,4) * t200 - Ifges(4,2) * t199 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t264 * (Ifges(5,5) * qJDD(3) - Ifges(5,6) * t200 + Ifges(5,3) * t199) / 0.2e1 + (t260 * t507 - t264 * t506) * qJDD(3) / 0.2e1 + t509 * t78 + t510 * t79 + (t123 * t77 + t15 * t510 + t16 * t509 + t170 * t23 + t2 * t34 + t3 * t33) * m(7) + t100 * (-mrSges(6,1) * t280 + mrSges(6,2) * t279) + t264 * t85 * t339; (-t159 + t44) * qJ(4) + (Ifges(7,4) * t145 + Ifges(7,2) * t144) * t453 + (-pkin(3) * t60 - qJ(4) * t59 - qJD(4) * t127 - t137 * t196) * m(5) + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + (Ifges(7,5) * t451 + Ifges(7,6) * t453 + Ifges(7,3) * t443 + t520) * t369 + (t171 * t515 + t172 * t469) * g(3) + (t469 * (t167 * t264 + t260 * t331) + t515 * t112) * g(1) + (t469 * (t165 * t264 - t260 * t332) + t515 * t110) * g(2) + t235 * t316 * t100 + (Ifges(7,5) * t145 + Ifges(7,6) * t144) * t443 - t85 * t374 / 0.2e1 + t173 * t246 / 0.2e1 + t516 * mrSges(7,3) - t126 * t356 - t127 * t357 + t237 * t10 + (qJ(4) * t43 + t100 * t478 - t49 * t73 - t50 * t74) * m(6) - (t188 * t306 + t235 * t303 - t284 * t310) * qJD(5) / 0.2e1 - (t188 * t274 + t235 * t273 - t275 * t284) * qJD(2) / 0.2e1 - t295 * t464 + (-Ifges(7,4) * t294 - Ifges(7,2) * t295) * t460 + (-Ifges(7,1) * t294 - Ifges(7,4) * t295) * t461 + t23 * (mrSges(7,1) * t295 - mrSges(7,2) * t294) + (-Ifges(7,5) * t294 - Ifges(7,6) * t295) * t447 - t294 * t463 + (-Ifges(7,5) * t108 - Ifges(7,6) * t107) * t442 + (-Ifges(7,1) * t108 - Ifges(7,4) * t107) * t450 + (-Ifges(7,4) * t108 - Ifges(7,2) * t107) * t452 + t84 * t339 + t145 * t457 + t259 * t459 + (-t49 * (mrSges(6,1) * t264 - mrSges(6,3) * t392) - t50 * (-mrSges(6,2) * t264 + mrSges(6,3) * t390) - t289 - t288) * qJD(2) + (Ifges(7,1) * t145 + Ifges(7,4) * t144) * t451 + (-m(6) * t271 - t522) * t449 - t196 * t194 - pkin(3) * t160 - t74 * t133 - t73 * t134 - t118 * t109 + t120 * t21 + t121 * t22 - t68 * mrSges(4,2) + t69 * mrSges(4,1) - t59 * mrSges(5,3) + t60 * mrSges(5,2) + t43 * t315 + t304 * t446 + t307 * t454 + t311 * t455 - t108 * t456 - t472 * t246 / 0.2e1 - (-Ifges(4,2) * t246 + t241 + t475) * t369 / 0.2e1 + (-t373 * t50 + t374 * t49 - t479) * mrSges(6,3) + t484 * t326 + (-m(5) * t127 - t354 - t486) * t135 + t487 * qJD(4) + (-m(5) * t126 + t355 - t489) * t136 + t493 * t458 + (mrSges(7,1) * t493 + mrSges(7,2) * t492) * t77 + t496 * t52 + (t485 / 0.2e1 - t287 / 0.2e1) * t268 + t263 * t32 / 0.2e1 + t504 * t79 + t505 * t78 + (t120 * t3 + t121 * t2 + t15 * t504 + t16 * t505 + t23 * t237 + t496 * t77) * m(7) + t506 * t199 + t507 * t200; t295 * t22 - t294 * t21 + t492 * t79 + t493 * t78 + t358 * qJD(3) + (t194 + t296) * t246 + t160 + (-qJD(3) * t77 - t516) * m(7) + (-qJD(3) * t100 - t246 * t299 + t271 + t282) * m(6) + (qJD(3) * t127 + t137 * t246 + t282 + t60) * m(5) + t522; t322 * t457 + (t134 - t425) * t50 + t361 + t512 * t453 + (t22 * t258 - t370 * t79 + (qJD(6) * t78 + t21) * t262) * pkin(5) + t521 + (-t133 + t426) * t49 - (Ifges(6,2) * t284 + t183 + t85) * t188 / 0.2e1 - t235 * (Ifges(6,5) * t188 + Ifges(6,6) * t284) / 0.2e1 - t100 * (-mrSges(6,1) * t284 + mrSges(6,2) * t188) + t284 * (Ifges(6,1) * t188 + t422) / 0.2e1 + t52 * t439 + (t2 * t258 + t262 * t3 + (-t15 * t258 + t16 * t262) * qJD(6)) * t462 - t99 * t458 - t18 * t78 - t17 * t79 - t528 + t527 + t84 * t444 - m(7) * (t15 * t17 + t16 * t18 - t439 * t77) + (-mrSges(6,2) * t117 + t116 * t499 - t427) * g(3) + (-(-t112 * t259 - t166 * t263) * mrSges(6,2) - t428 + t499 * (t112 * t263 - t166 * t259)) * g(1) + (-(-t110 * t259 - t164 * t263) * mrSges(6,2) - t429 + t499 * (t110 * t263 - t164 * t259)) * g(2); t40 * t450 - t15 * t78 + t16 * t79 - g(1) * t428 - g(2) * t429 - g(3) * t427 + (t41 + t512) * t453 + t521;];
tau  = t1;
