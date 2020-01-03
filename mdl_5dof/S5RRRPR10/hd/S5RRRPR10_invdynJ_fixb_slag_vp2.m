% Calculate vector of inverse dynamics joint torques for
% S5RRRPR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:40
% EndTime: 2019-12-31 21:27:35
% DurationCPUTime: 29.66s
% Computational Cost: add. (12197->785), mult. (30206->1088), div. (0->0), fcn. (23851->14), ass. (0->350)
t296 = cos(qJ(2));
t292 = sin(qJ(2));
t286 = sin(pkin(5));
t381 = qJD(1) * t286;
t351 = t292 * t381;
t288 = cos(pkin(5));
t380 = qJD(1) * t288;
t368 = pkin(1) * t380;
t222 = -pkin(7) * t351 + t296 * t368;
t310 = (pkin(2) * t292 - pkin(8) * t296) * t286;
t223 = qJD(1) * t310;
t291 = sin(qJ(3));
t295 = cos(qJ(3));
t147 = -t291 * t222 + t295 * t223;
t289 = -qJ(4) - pkin(8);
t341 = qJD(3) * t289;
t385 = t295 * t296;
t525 = -(pkin(3) * t292 - qJ(4) * t385) * t381 - t147 - qJD(4) * t291 + t295 * t341;
t148 = t295 * t222 + t291 * t223;
t350 = t296 * t381;
t332 = t291 * t350;
t524 = -qJ(4) * t332 - qJD(4) * t295 - t291 * t341 + t148;
t253 = qJD(3) - t350;
t437 = -t253 / 0.2e1;
t270 = qJD(2) + t380;
t196 = t270 * t295 - t291 * t351;
t197 = t270 * t291 + t295 * t351;
t285 = sin(pkin(10));
t287 = cos(pkin(10));
t312 = t196 * t285 + t287 * t197;
t442 = -t312 / 0.2e1;
t340 = t287 * t196 - t197 * t285;
t444 = -t340 / 0.2e1;
t523 = -Ifges(5,4) * t442 - Ifges(5,2) * t444 - Ifges(5,6) * t437;
t131 = qJD(5) - t340;
t445 = t131 / 0.2e1;
t290 = sin(qJ(5));
t294 = cos(qJ(5));
t107 = t253 * t290 + t294 * t312;
t451 = t107 / 0.2e1;
t106 = t253 * t294 - t290 * t312;
t453 = t106 / 0.2e1;
t522 = Ifges(6,5) * t451 + Ifges(6,6) * t453 + Ifges(6,3) * t445;
t521 = Ifges(4,3) + Ifges(5,3);
t324 = -mrSges(6,1) * t294 + mrSges(6,2) * t290;
t498 = m(6) * pkin(4) + mrSges(5,1) - t324;
t338 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t493 = t525 * t285 - t524 * t287;
t225 = pkin(7) * t350 + t292 * t368;
t377 = qJD(3) * t291;
t489 = -t225 + (-t332 + t377) * pkin(3);
t183 = pkin(8) * t270 + t225;
t211 = (-pkin(2) * t296 - pkin(8) * t292 - pkin(1)) * t286;
t187 = qJD(1) * t211;
t122 = -t183 * t291 + t295 * t187;
t97 = -qJ(4) * t197 + t122;
t91 = pkin(3) * t253 + t97;
t123 = t183 * t295 + t187 * t291;
t98 = qJ(4) * t196 + t123;
t95 = t287 * t98;
t46 = t285 * t91 + t95;
t43 = pkin(9) * t253 + t46;
t182 = -t270 * pkin(2) - t222;
t132 = -t196 * pkin(3) + qJD(4) + t182;
t53 = -pkin(4) * t340 - pkin(9) * t312 + t132;
t13 = -t290 * t43 + t294 * t53;
t14 = t290 * t53 + t294 * t43;
t520 = mrSges(5,1) * t132 + mrSges(6,1) * t13 - mrSges(6,2) * t14 - mrSges(5,3) * t46 + t522 - t523;
t519 = -pkin(9) * t351 + t493;
t243 = t285 * t295 + t287 * t291;
t175 = t243 * t350;
t242 = t285 * t291 - t287 * t295;
t176 = t242 * t350;
t231 = t243 * qJD(3);
t232 = t242 * qJD(3);
t518 = t489 + (-t176 + t232) * pkin(9) + (-t175 + t231) * pkin(4);
t373 = qJD(1) * qJD(2);
t229 = (qJDD(1) * t292 + t296 * t373) * t286;
t371 = qJDD(1) * t288;
t269 = qJDD(2) + t371;
t124 = qJD(3) * t196 + t229 * t295 + t269 * t291;
t228 = (-qJDD(1) * t296 + t292 * t373) * t286;
t215 = qJDD(3) + t228;
t372 = qJDD(1) * t286;
t509 = pkin(7) * t372 + qJD(2) * t368;
t510 = -pkin(7) * t286 * t373 + pkin(1) * t371;
t156 = t292 * t510 + t296 * t509;
t141 = pkin(8) * t269 + t156;
t146 = -pkin(1) * t372 + pkin(2) * t228 - pkin(8) * t229;
t52 = -qJD(3) * t123 - t141 * t291 + t295 * t146;
t29 = pkin(3) * t215 - qJ(4) * t124 - qJD(4) * t197 + t52;
t125 = -qJD(3) * t197 - t229 * t291 + t269 * t295;
t376 = qJD(3) * t295;
t51 = t295 * t141 + t291 * t146 - t183 * t377 + t187 * t376;
t32 = qJ(4) * t125 + qJD(4) * t196 + t51;
t11 = t285 * t29 + t287 * t32;
t438 = t215 / 0.2e1;
t82 = t124 * t287 + t125 * t285;
t461 = t82 / 0.2e1;
t81 = -t124 * t285 + t125 * t287;
t462 = t81 / 0.2e1;
t80 = qJDD(5) - t81;
t463 = t80 / 0.2e1;
t38 = -qJD(5) * t107 + t215 * t294 - t290 * t82;
t471 = t38 / 0.2e1;
t37 = qJD(5) * t106 + t215 * t290 + t294 * t82;
t472 = t37 / 0.2e1;
t157 = -t292 * t509 + t296 * t510;
t142 = -t269 * pkin(2) - t157;
t88 = -t125 * pkin(3) + qJDD(4) + t142;
t15 = -t81 * pkin(4) - t82 * pkin(9) + t88;
t9 = pkin(9) * t215 + t11;
t1 = qJD(5) * t13 + t15 * t290 + t294 * t9;
t2 = -qJD(5) * t14 + t15 * t294 - t290 * t9;
t481 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t5 = Ifges(6,5) * t37 + Ifges(6,6) * t38 + Ifges(6,3) * t80;
t516 = t481 + mrSges(5,1) * t88 - mrSges(5,3) * t11 + Ifges(6,5) * t472 + Ifges(6,6) * t471 + Ifges(6,3) * t463 + t5 / 0.2e1 + (-t438 - t215 / 0.2e1) * Ifges(5,6) + (-t462 - t81 / 0.2e1) * Ifges(5,2) + (-t461 - t82 / 0.2e1) * Ifges(5,4);
t436 = t253 / 0.2e1;
t441 = t312 / 0.2e1;
t443 = t340 / 0.2e1;
t402 = t285 * t98;
t45 = t287 * t91 - t402;
t484 = mrSges(5,2) * t132 - t45 * mrSges(5,3);
t86 = Ifges(5,1) * t312 + Ifges(5,4) * t340 + t253 * Ifges(5,5);
t515 = t484 + Ifges(5,1) * t441 + Ifges(5,4) * t443 + Ifges(5,5) * t436 + t86 / 0.2e1;
t513 = -Ifges(5,4) * t441 - Ifges(5,2) * t443 - Ifges(5,6) * t436 + t520 + t522;
t393 = t286 * t292;
t233 = t288 * t295 - t291 * t393;
t311 = t233 * pkin(3);
t512 = t156 * mrSges(3,2);
t495 = t524 * t285 + t525 * t287;
t284 = qJ(3) + pkin(10);
t281 = sin(t284);
t282 = cos(t284);
t326 = -mrSges(4,1) * t295 + mrSges(4,2) * t291;
t511 = -m(4) * pkin(2) + t338 * t281 - t498 * t282 + t326;
t508 = Ifges(4,5) * t124 + Ifges(5,5) * t82 + Ifges(4,6) * t125 + Ifges(5,6) * t81 + t521 * t215;
t433 = cos(qJ(1));
t352 = t433 * t296;
t293 = sin(qJ(1));
t387 = t292 * t293;
t238 = -t288 * t387 + t352;
t391 = t286 * t295;
t172 = -t238 * t291 + t293 * t391;
t446 = -t131 / 0.2e1;
t452 = -t107 / 0.2e1;
t454 = -t106 / 0.2e1;
t507 = Ifges(6,5) * t452 + Ifges(6,6) * t454 + Ifges(6,3) * t446 - t520 + t523;
t10 = -t285 * t32 + t287 * t29;
t506 = -t52 * mrSges(4,1) - t10 * mrSges(5,1) + t51 * mrSges(4,2) + t11 * mrSges(5,2);
t504 = t88 * mrSges(5,2) - mrSges(5,3) * t10 + 0.2e1 * Ifges(5,1) * t461 + 0.2e1 * Ifges(5,4) * t462 + 0.2e1 * Ifges(5,5) * t438;
t501 = -m(6) - m(5);
t448 = t124 / 0.2e1;
t447 = t125 / 0.2e1;
t280 = pkin(3) * t295 + pkin(2);
t159 = pkin(4) * t242 - pkin(9) * t243 - t280;
t257 = t289 * t295;
t345 = t289 * t291;
t181 = -t287 * t257 + t285 * t345;
t104 = t159 * t294 - t181 * t290;
t500 = qJD(5) * t104 + t518 * t290 + t294 * t519;
t105 = t159 * t290 + t181 * t294;
t499 = -qJD(5) * t105 - t290 * t519 + t518 * t294;
t497 = -m(4) * pkin(8) - mrSges(4,3) - mrSges(5,3);
t113 = mrSges(5,1) * t253 - mrSges(5,3) * t312;
t64 = -mrSges(6,1) * t106 + mrSges(6,2) * t107;
t496 = t113 - t64;
t494 = pkin(4) * t351 - t495;
t337 = mrSges(3,3) * t351;
t492 = -mrSges(3,1) * t270 - mrSges(4,1) * t196 + mrSges(4,2) * t197 + t337;
t152 = t176 * t290 + t294 * t351;
t374 = qJD(5) * t294;
t307 = -t290 * t232 + t243 * t374;
t491 = t152 + t307;
t153 = -t176 * t294 + t290 * t351;
t375 = qJD(5) * t290;
t306 = t294 * t232 + t243 * t375;
t490 = t153 + t306;
t390 = t286 * t296;
t241 = t288 * t292 * pkin(1) + pkin(7) * t390;
t210 = pkin(8) * t288 + t241;
t145 = t295 * t210 + t291 * t211;
t271 = pkin(7) * t393;
t431 = pkin(1) * t296;
t240 = t288 * t431 - t271;
t353 = t433 * t292;
t386 = t293 * t296;
t236 = t288 * t353 + t386;
t354 = t286 * t433;
t300 = t236 * t291 + t295 * t354;
t487 = -t291 * t52 + t295 * t51;
t486 = t1 * t294 - t2 * t290;
t323 = mrSges(6,1) * t290 + mrSges(6,2) * t294;
t485 = -t323 + t497;
t483 = mrSges(3,2) + t485;
t482 = mrSges(3,1) - t511;
t478 = t286 ^ 2;
t476 = Ifges(6,1) * t472 + Ifges(6,4) * t471 + Ifges(6,5) * t463;
t475 = m(5) * pkin(3);
t414 = Ifges(6,4) * t107;
t40 = Ifges(6,2) * t106 + Ifges(6,6) * t131 + t414;
t467 = t40 / 0.2e1;
t103 = Ifges(6,4) * t106;
t41 = Ifges(6,1) * t107 + Ifges(6,5) * t131 + t103;
t466 = -t41 / 0.2e1;
t465 = Ifges(4,4) * t448 + Ifges(4,2) * t447 + Ifges(4,6) * t438;
t464 = Ifges(4,1) * t448 + Ifges(4,4) * t447 + Ifges(4,5) * t438;
t458 = -t86 / 0.2e1;
t456 = pkin(1) * mrSges(3,1);
t455 = pkin(1) * mrSges(3,2);
t405 = t197 * Ifges(4,4);
t117 = t196 * Ifges(4,2) + t253 * Ifges(4,6) + t405;
t450 = t117 / 0.2e1;
t188 = Ifges(4,4) * t196;
t118 = t197 * Ifges(4,1) + t253 * Ifges(4,5) + t188;
t449 = t118 / 0.2e1;
t439 = t197 / 0.2e1;
t435 = t288 / 0.2e1;
t434 = t294 / 0.2e1;
t430 = pkin(3) * t197;
t429 = pkin(3) * t285;
t428 = pkin(3) * t287;
t378 = qJD(2) * t286;
t348 = t296 * t378;
t170 = qJD(3) * t233 + t295 * t348;
t234 = t288 * t291 + t292 * t391;
t349 = t292 * t378;
t224 = qJD(2) * t310;
t226 = t240 * qJD(2);
t93 = -qJD(3) * t145 + t295 * t224 - t226 * t291;
t57 = pkin(3) * t349 - qJ(4) * t170 - qJD(4) * t234 + t93;
t169 = -qJD(3) * t234 - t291 * t348;
t92 = -t210 * t377 + t211 * t376 + t291 * t224 + t295 * t226;
t63 = qJ(4) * t169 + qJD(4) * t233 + t92;
t23 = t285 * t57 + t287 * t63;
t418 = Ifges(3,4) * t292;
t417 = Ifges(3,4) * t296;
t416 = Ifges(4,4) * t291;
t415 = Ifges(4,4) * t295;
t413 = Ifges(6,4) * t290;
t412 = Ifges(6,4) * t294;
t411 = Ifges(3,6) * t270;
t410 = t122 * mrSges(4,3);
t409 = t123 * mrSges(4,3);
t408 = t340 * Ifges(5,6);
t407 = t312 * Ifges(5,5);
t406 = t196 * Ifges(4,6);
t404 = t197 * Ifges(4,5);
t403 = t270 * Ifges(3,5);
t399 = t340 * t290;
t398 = t340 * t294;
t395 = t243 * t290;
t394 = t243 * t294;
t392 = t286 * t293;
t144 = -t291 * t210 + t295 * t211;
t102 = -pkin(3) * t390 - t234 * qJ(4) + t144;
t114 = qJ(4) * t233 + t145;
t60 = t285 * t102 + t287 * t114;
t382 = t433 * pkin(1) + pkin(7) * t392;
t362 = t290 * t390;
t360 = t291 * t392;
t358 = t294 * t390;
t357 = t41 * t434;
t355 = Ifges(3,5) * t229 - Ifges(3,6) * t228 + Ifges(3,3) * t269;
t31 = -t81 * mrSges(5,1) + t82 * mrSges(5,2);
t343 = -t375 / 0.2e1;
t342 = -pkin(1) * t293 + pkin(7) * t354;
t163 = t236 * t282 - t281 * t354;
t162 = -t236 * t281 - t282 * t354;
t260 = t291 * t354;
t339 = -t236 * t295 + t260;
t336 = mrSges(3,3) * t350;
t329 = t172 * pkin(3);
t328 = mrSges(3,2) + t497;
t327 = t233 * mrSges(4,1) - t234 * mrSges(4,2);
t322 = Ifges(4,1) * t295 - t416;
t321 = Ifges(6,1) * t294 - t413;
t320 = Ifges(3,2) * t296 + t418;
t319 = -Ifges(4,2) * t291 + t415;
t318 = -Ifges(6,2) * t290 + t412;
t317 = Ifges(4,5) * t295 - Ifges(4,6) * t291;
t316 = Ifges(6,5) * t294 - Ifges(6,6) * t290;
t22 = -t285 * t63 + t287 * t57;
t56 = -pkin(9) * t390 + t60;
t154 = -t287 * t233 + t234 * t285;
t155 = t233 * t285 + t234 * t287;
t209 = t271 + (-pkin(2) - t431) * t288;
t158 = t209 - t311;
t83 = t154 * pkin(4) - t155 * pkin(9) + t158;
t25 = t290 * t83 + t294 * t56;
t24 = -t290 * t56 + t294 * t83;
t74 = -mrSges(6,2) * t131 + mrSges(6,3) * t106;
t75 = mrSges(6,1) * t131 - mrSges(6,3) * t107;
t314 = -t290 * t75 + t294 * t74;
t237 = t288 * t386 + t353;
t313 = pkin(3) * t360 - t237 * t289 + t238 * t280 + t382;
t59 = t287 * t102 - t285 * t114;
t129 = -t290 * t155 - t358;
t309 = -t294 * t155 + t362;
t42 = -pkin(4) * t253 - t45;
t308 = t42 * t323;
t235 = -t288 * t352 + t387;
t298 = -g(1) * t237 - g(2) * t235 + g(3) * t390;
t227 = t241 * qJD(2);
t140 = -t169 * pkin(3) + t227;
t279 = -pkin(4) - t428;
t265 = Ifges(3,4) * t350;
t239 = (-mrSges(3,1) * t296 + mrSges(3,2) * t292) * t286;
t221 = -t270 * mrSges(3,2) + t336;
t208 = t281 * t288 + t282 * t393;
t180 = -t257 * t285 - t287 * t345;
t178 = Ifges(3,1) * t351 + t265 + t403;
t177 = t320 * t381 + t411;
t173 = t238 * t295 + t360;
t167 = t238 * t282 + t281 * t392;
t166 = t238 * t281 - t282 * t392;
t161 = mrSges(4,1) * t253 - mrSges(4,3) * t197;
t160 = -mrSges(4,2) * t253 + mrSges(4,3) * t196;
t127 = t167 * t294 + t237 * t290;
t126 = -t167 * t290 + t237 * t294;
t116 = t253 * Ifges(4,3) + t404 + t406;
t112 = -mrSges(5,2) * t253 + mrSges(5,3) * t340;
t111 = t169 * t285 + t170 * t287;
t110 = -t287 * t169 + t170 * t285;
t100 = -mrSges(4,2) * t215 + mrSges(4,3) * t125;
t99 = mrSges(4,1) * t215 - mrSges(4,3) * t124;
t89 = -mrSges(5,1) * t340 + mrSges(5,2) * t312;
t87 = -mrSges(4,1) * t125 + mrSges(4,2) * t124;
t84 = t253 * Ifges(5,3) + t407 + t408;
t71 = pkin(4) * t312 - pkin(9) * t340 + t430;
t68 = qJD(5) * t309 - t290 * t111 + t294 * t349;
t67 = qJD(5) * t129 + t294 * t111 + t290 * t349;
t62 = mrSges(5,1) * t215 - mrSges(5,3) * t82;
t61 = -mrSges(5,2) * t215 + mrSges(5,3) * t81;
t55 = pkin(4) * t390 - t59;
t48 = t287 * t97 - t402;
t47 = t285 * t97 + t95;
t44 = t110 * pkin(4) - t111 * pkin(9) + t140;
t21 = pkin(9) * t349 + t23;
t20 = -pkin(4) * t349 - t22;
t19 = t290 * t71 + t294 * t48;
t18 = -t290 * t48 + t294 * t71;
t17 = -mrSges(6,2) * t80 + mrSges(6,3) * t38;
t16 = mrSges(6,1) * t80 - mrSges(6,3) * t37;
t12 = -mrSges(6,1) * t38 + mrSges(6,2) * t37;
t8 = -pkin(4) * t215 - t10;
t6 = t37 * Ifges(6,4) + t38 * Ifges(6,2) + t80 * Ifges(6,6);
t4 = -qJD(5) * t25 - t21 * t290 + t294 * t44;
t3 = qJD(5) * t24 + t21 * t294 + t290 * t44;
t7 = [(Ifges(4,1) * t234 + Ifges(4,4) * t233) * t448 + (Ifges(6,1) * t67 + Ifges(6,4) * t68) * t451 + (Ifges(4,4) * t234 + Ifges(4,2) * t233) * t447 + (Ifges(6,4) * t67 + Ifges(6,2) * t68) * t453 + (Ifges(4,5) * t234 + Ifges(4,6) * t233) * t438 + t492 * t227 + m(6) * (t1 * t25 + t13 * t4 + t14 * t3 + t2 * t24 + t20 * t42 + t55 * t8) + m(5) * (t10 * t59 + t11 * t60 + t132 * t140 + t158 * t88 + t22 * t45 + t23 * t46) + m(4) * (t122 * t93 + t123 * t92 + t142 * t209 + t144 * t52 + t145 * t51 + t182 * t227) + (-pkin(1) * t239 * t286 + Ifges(2,3)) * qJDD(1) + t504 * t155 - t288 * t512 + (-t508 / 0.2e1 + mrSges(3,3) * t156 - Ifges(4,5) * t448 - Ifges(5,5) * t461 - Ifges(4,6) * t447 - Ifges(5,6) * t462 - t521 * t438 + t506) * t390 + t182 * (-mrSges(4,1) * t169 + mrSges(4,2) * t170) + (-t122 * t170 + t123 * t169 + t233 * t51 - t234 * t52) * mrSges(4,3) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t478 + t156 * t241 + t157 * t240 - t222 * t227 + t225 * t226) + (-m(3) * t342 + t236 * mrSges(3,1) - mrSges(3,3) * t354 + t293 * mrSges(2,1) + t433 * mrSges(2,2) - m(4) * (-pkin(2) * t236 + t342) - t339 * mrSges(4,1) - t300 * mrSges(4,2) + t338 * t162 + t498 * t163 + (t323 - t328) * t235 - t501 * (-pkin(3) * t260 - t235 * t289 + t236 * t280 - t342)) * g(1) + (-t240 * mrSges(3,3) + Ifges(3,5) * t435 + (t292 * Ifges(3,1) + t417 - t455) * t286) * t229 + ((-t222 * mrSges(3,3) + t178 / 0.2e1 + t403 / 0.2e1 + (-t455 + t417 / 0.2e1) * t381) * t296 + (-t225 * mrSges(3,3) + t84 / 0.2e1 + t116 / 0.2e1 - t177 / 0.2e1 + t406 / 0.2e1 + t404 / 0.2e1 - t123 * mrSges(4,2) + t122 * mrSges(4,1) + t408 / 0.2e1 + t407 / 0.2e1 + t45 * mrSges(5,1) - t46 * mrSges(5,2) - t411 / 0.2e1 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t253 + (-t456 - t418 / 0.2e1 + (Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t296) * t381) * t292) * t378 - (t241 * mrSges(3,3) + Ifges(3,6) * t435 + (t320 + t456) * t286) * t228 + t209 * t87 + t196 * (Ifges(4,4) * t170 + Ifges(4,2) * t169) / 0.2e1 + (-t433 * mrSges(2,1) - m(4) * (pkin(2) * t238 + t382) - t173 * mrSges(4,1) - t172 * mrSges(4,2) - m(6) * (pkin(4) * t167 + t313) - t127 * mrSges(6,1) - t126 * mrSges(6,2) - m(5) * t313 - t167 * mrSges(5,1) - m(3) * t382 - t238 * mrSges(3,1) + (-mrSges(3,3) * t286 + mrSges(2,2)) * t293 + t338 * t166 + t328 * t237) * g(2) + (t240 * mrSges(3,1) - t241 * mrSges(3,2) + Ifges(3,3) * t435 + (Ifges(3,5) * t292 + Ifges(3,6) * t296) * t286) * t269 + (t1 * t129 - t13 * t67 + t14 * t68 + t2 * t309) * mrSges(6,3) + (-Ifges(6,4) * t309 + Ifges(6,2) * t129) * t471 + (-Ifges(6,5) * t309 + Ifges(6,6) * t129) * t463 + t8 * (-mrSges(6,1) * t129 - mrSges(6,2) * t309) - t309 * t476 + (-Ifges(6,1) * t309 + Ifges(6,4) * t129) * t472 + t513 * t110 + t92 * t160 + t93 * t161 + t158 * t31 + t140 * t89 + t144 * t99 + t145 * t100 - t142 * t327 + t129 * t6 / 0.2e1 + t23 * t112 + t22 * t113 + t3 * t74 + t4 * t75 + t67 * t41 / 0.2e1 + t42 * (-mrSges(6,1) * t68 + mrSges(6,2) * t67) + t60 * t61 + t59 * t62 + t20 * t64 + t55 * t12 + t24 * t16 + t25 * t17 + t226 * t221 + t515 * t111 + t516 * t154 + (Ifges(6,5) * t67 + Ifges(6,6) * t68) * t445 + t355 * t435 + (Ifges(4,5) * t170 + Ifges(4,6) * t169) * t436 + (Ifges(4,1) * t170 + Ifges(4,4) * t169) * t439 + t170 * t449 + t169 * t450 + t234 * t464 + t233 * t465 + t68 * t467 + t157 * (mrSges(3,1) * t288 - mrSges(3,3) * t393); (-pkin(2) * t142 - t122 * t147 - t123 * t148) * m(4) + (-t123 * (-mrSges(4,3) * t291 * t296 - mrSges(4,2) * t292) - t122 * (mrSges(4,1) * t292 - mrSges(4,3) * t385)) * t381 + (t336 - t221) * t222 + (-m(4) * t182 + t337 - t492) * t225 + (t239 + t501 * t280 * t390 + (t511 * t296 + (-t289 * t501 + t485) * t292) * t286) * g(3) + (pkin(1) * (mrSges(3,1) * t292 + mrSges(3,2) * t296) - t292 * (Ifges(3,1) * t296 - t418) / 0.2e1) * qJD(1) ^ 2 * t478 + (-t161 * t376 - t160 * t377 + m(4) * ((-t122 * t295 - t123 * t291) * qJD(3) + t487) - t291 * t99 + t295 * t100) * pkin(8) + t487 * mrSges(4,3) + t253 * t182 * (mrSges(4,1) * t291 + mrSges(4,2) * t295) + t355 + (t316 * t463 + t318 * t471 + t321 * t472 + t323 * t8 + t343 * t41 + t504) * t243 + t181 * t61 - t512 + (-Ifges(6,5) * t306 - Ifges(6,6) * t307) * t445 + (Ifges(6,5) * t153 + Ifges(6,6) * t152) * t446 - t6 * t395 / 0.2e1 + (t501 * (-t235 * t280 - t236 * t289) + t483 * t236 + t482 * t235) * g(2) + (t501 * (-t237 * t280 - t238 * t289) + t483 * t238 + t482 * t237) * g(1) + (-Ifges(6,4) * t306 - Ifges(6,2) * t307) * t453 + (Ifges(6,4) * t153 + Ifges(6,2) * t152) * t454 + (t12 - t62) * t180 + t177 * t351 / 0.2e1 + t489 * t89 - t491 * t40 / 0.2e1 + (mrSges(6,1) * t491 - mrSges(6,2) * t490) * t42 + (-t1 * t395 + t13 * t490 - t14 * t491 - t2 * t394) * mrSges(6,3) - (t253 * (Ifges(4,3) * t292 + t296 * t317) + t197 * (Ifges(4,5) * t292 + t296 * t322) + t196 * (Ifges(4,6) * t292 + t296 * t319) + t270 * (Ifges(3,5) * t296 - Ifges(3,6) * t292) + (-Ifges(3,2) * t351 + t295 * t118 + t178 + t265) * t296 + (t116 + t84) * t292) * t381 / 0.2e1 + t499 * t75 + (t1 * t105 + t104 * t2 + t13 * t499 + t14 * t500 + t180 * t8 + t42 * t494) * m(6) + t500 * t74 + t513 * t231 - t148 * t160 - t147 * t161 + t157 * mrSges(3,1) + t142 * t326 + t104 * t16 + t105 * t17 - pkin(2) * t87 + (-Ifges(6,1) * t306 - Ifges(6,4) * t307) * t451 + (Ifges(6,1) * t153 + Ifges(6,4) * t152) * t452 + t493 * t112 + t494 * t64 + t495 * t113 + (-t10 * t180 + t11 * t181 + t132 * t489 - t280 * t88 + t45 * t495 + t46 * t493) * m(5) - (t357 + t515) * t232 + t516 * t242 + (t196 * t319 + t197 * t322 + t253 * t317) * qJD(3) / 0.2e1 - t280 * t31 + (t132 * t176 + t351 * t46) * mrSges(5,2) + (-Ifges(5,5) * t176 + Ifges(5,3) * t351) * t437 - t45 * (mrSges(5,1) * t351 + mrSges(5,3) * t176) + (-Ifges(5,4) * t176 + Ifges(5,6) * t351) * t444 + (-Ifges(5,1) * t176 + Ifges(5,5) * t351) * t442 + (Ifges(4,5) * t291 + Ifges(4,6) * t295) * t438 + (Ifges(4,2) * t295 + t416) * t447 + (Ifges(4,1) * t291 + t415) * t448 + t332 * t450 - t176 * t458 + t291 * t464 + t295 * t465 + t153 * t466 + t394 * t476 + (-t410 + t449) * t376 + t507 * t175 + (-t117 / 0.2e1 - t409) * t377; (t308 + t357) * qJD(5) - m(5) * (t132 * t430 + t46 * t48) - t506 - t89 * t430 + (-g(1) * t167 - g(2) * t163 - g(3) * t208 + (-t375 + t399) * t14 + (-t374 + t398) * t13 + t486) * mrSges(6,3) + (-t75 * t374 - t74 * t375 + m(6) * ((-t13 * t294 - t14 * t290) * qJD(5) + t486) - t290 * t16 + t294 * t17) * (pkin(9) + t429) + t508 - t182 * (mrSges(4,1) * t197 + mrSges(4,2) * t196) + (Ifges(5,1) * t442 + Ifges(5,4) * t444 + Ifges(5,5) * t437 + t316 * t446 + t318 * t454 + t321 * t452 - t308 + t458 - t484) * t340 + (-t327 - m(5) * t311 + t208 * mrSges(5,2) - m(6) * (pkin(9) * t208 + t311) - t498 * (-t281 * t393 + t282 * t288)) * g(3) + (-t13 * t18 - t14 * t19 + t279 * t8 - t42 * t47) * m(6) - t122 * t160 + t123 * t161 + t8 * t324 - t48 * t112 - t19 * t74 - t18 * t75 - t197 * (Ifges(4,1) * t196 - t405) / 0.2e1 + (m(5) * t45 + t496) * t47 + (-t172 * mrSges(4,1) + t173 * mrSges(4,2) - m(5) * t329 + t167 * mrSges(5,2) - m(6) * (pkin(9) * t167 + t329) + t498 * t166) * g(1) + (t163 * mrSges(5,2) - m(6) * (-pkin(3) * t300 + t163 * pkin(9)) - t339 * mrSges(4,2) + (t475 + mrSges(4,1)) * t300 - t498 * t162) * g(2) + t279 * t12 + (t106 * t318 + t107 * t321 + t131 * t316) * qJD(5) / 0.2e1 - (-Ifges(4,2) * t197 + t118 + t188) * t196 / 0.2e1 + t40 * t343 + t197 * t409 + t196 * t410 + t62 * t428 + t61 * t429 + t6 * t434 + (Ifges(4,5) * t196 - Ifges(4,6) * t197) * t437 + t117 * t439 + (Ifges(6,5) * t290 + Ifges(6,6) * t294) * t463 + t398 * t466 + t399 * t467 + (Ifges(6,2) * t294 + t413) * t471 + (Ifges(6,1) * t290 + t412) * t472 + (t10 * t287 + t11 * t285) * t475 + t290 * t476 + t507 * t312; t294 * t16 + t290 * t17 + t496 * t312 + t314 * qJD(5) + (-t112 - t314) * t340 + t31 + (t1 * t290 - t312 * t42 + t2 * t294 + t298 + t131 * (-t13 * t290 + t14 * t294)) * m(6) + (t312 * t45 - t340 * t46 + t298 + t88) * m(5); -t42 * (mrSges(6,1) * t107 + mrSges(6,2) * t106) + (Ifges(6,1) * t106 - t414) * t452 + t40 * t451 + (Ifges(6,5) * t106 - Ifges(6,6) * t107) * t446 - t13 * t74 + t14 * t75 - g(1) * (mrSges(6,1) * t126 - mrSges(6,2) * t127) - g(2) * ((-t163 * t290 + t235 * t294) * mrSges(6,1) + (-t163 * t294 - t235 * t290) * mrSges(6,2)) - g(3) * ((-t208 * t290 - t358) * mrSges(6,1) + (-t208 * t294 + t362) * mrSges(6,2)) + (t106 * t13 + t107 * t14) * mrSges(6,3) + t5 + (-Ifges(6,2) * t107 + t103 + t41) * t454 + t481;];
tau = t7;
