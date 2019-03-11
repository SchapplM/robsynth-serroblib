% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:36:04
% EndTime: 2019-03-08 22:36:51
% DurationCPUTime: 27.21s
% Computational Cost: add. (10367->867), mult. (25867->1179), div. (0->0), fcn. (21068->14), ass. (0->402)
t290 = sin(pkin(7));
t297 = sin(qJ(3));
t418 = qJD(3) * t297;
t388 = t290 * t418;
t274 = pkin(3) * t388;
t301 = cos(qJ(3));
t336 = pkin(10) * t297 - qJ(4) * t301;
t416 = qJD(4) * t297;
t143 = t274 + (qJD(3) * t336 - t416) * t290;
t448 = t290 * t297;
t283 = pkin(9) * t448;
t293 = cos(pkin(7));
t394 = -pkin(2) * t301 - pkin(3);
t148 = pkin(4) * t448 + t283 + (-pkin(10) + t394) * t293;
t452 = qJ(4) * t297;
t495 = pkin(3) + pkin(10);
t328 = -t301 * t495 - t452;
t166 = (-pkin(2) + t328) * t290;
t440 = t293 * t297;
t287 = pkin(2) * t440;
t446 = t290 * t301;
t494 = pkin(4) + pkin(9);
t175 = (t446 * t494 + t287) * qJD(3);
t291 = sin(pkin(6));
t298 = sin(qJ(2));
t433 = t298 * t301;
t302 = cos(qJ(2));
t434 = t297 * t302;
t327 = t293 * t433 + t434;
t196 = t327 * t291;
t178 = qJD(1) * t196;
t296 = sin(qJ(5));
t300 = cos(qJ(5));
t423 = qJD(1) * t291;
t392 = t298 * t423;
t363 = t290 * t392;
t414 = qJD(5) * t300;
t415 = qJD(5) * t296;
t537 = t148 * t414 - t166 * t415 + (t143 - t363) * t300 + (t175 - t178) * t296;
t417 = qJD(3) * t301;
t386 = t293 * t417;
t275 = pkin(2) * t386;
t288 = t293 * qJD(4);
t147 = -t388 * t494 + t275 + t288;
t444 = t291 * t298;
t398 = t297 * t444;
t366 = t293 * t398;
t391 = t302 * t423;
t179 = -qJD(1) * t366 + t301 * t391;
t552 = t147 - t179;
t387 = t290 * t417;
t551 = -pkin(11) * t387 - t537;
t236 = t293 * t296 + t300 * t446;
t160 = -qJD(5) * t236 + t296 * t388;
t399 = t296 * t446;
t161 = -qJD(5) * t399 + t293 * t414 - t300 * t388;
t550 = pkin(5) * t161 - pkin(11) * t160 + t552;
t413 = qJD(5) * t495;
t250 = qJD(2) * pkin(2) + t391;
t228 = t250 * t440;
t421 = qJD(2) * t290;
t240 = pkin(9) * t421 + t392;
t294 = cos(pkin(6));
t422 = qJD(1) * t294;
t393 = t290 * t422;
t121 = t301 * t240 + t297 * t393 + t228;
t389 = t301 * t421;
t108 = pkin(4) * t389 + t121;
t390 = t297 * t421;
t272 = pkin(3) * t390;
t174 = t336 * t421 + t272;
t66 = t296 * t108 + t300 * t174;
t549 = -t300 * t413 - t66;
t295 = sin(qJ(6));
t299 = cos(qJ(6));
t347 = t295 * mrSges(7,1) + t299 * mrSges(7,2);
t496 = -m(7) - m(6);
t548 = -mrSges(4,1) + mrSges(5,2);
t309 = pkin(10) * t496 - mrSges(6,3) - t347 + t548;
t351 = pkin(5) * t300 + pkin(11) * t296;
t439 = t293 * t301;
t425 = t250 * t439 + t301 * t393;
t547 = qJD(5) * t351 + qJD(4) - (-t240 + (-pkin(4) - t351) * t421) * t297 - t425;
t546 = pkin(11) * t389 - t549;
t523 = t296 * t148 + t300 * t166;
t79 = pkin(11) * t448 + t523;
t239 = pkin(9) * t446 + t287;
t202 = -t293 * qJ(4) - t239;
t165 = pkin(4) * t446 - t202;
t237 = t293 * t300 - t399;
t89 = pkin(5) * t236 - pkin(11) * t237 + t165;
t32 = t295 * t89 + t299 * t79;
t545 = -qJD(6) * t32 + t295 * t551 + t550 * t299;
t31 = -t295 * t79 + t299 * t89;
t544 = qJD(6) * t31 + t550 * t295 - t299 * t551;
t261 = qJD(5) + t390;
t277 = t293 * t422;
t117 = t277 + (qJD(2) * t328 - t250) * t290;
t282 = qJD(2) * t293 + qJD(3);
t107 = -t297 * (pkin(4) * t421 + t240) + t425;
t516 = qJD(4) - t107;
t75 = -t282 * t495 + t516;
t41 = t117 * t300 + t296 * t75;
t30 = pkin(11) * t261 + t41;
t198 = -t282 * t296 - t300 * t389;
t316 = -t282 * t300 + t296 * t389;
t270 = t282 * qJ(4);
t88 = t270 + t108;
t50 = -pkin(5) * t198 + pkin(11) * t316 + t88;
t15 = -t295 * t30 + t299 * t50;
t543 = t15 * mrSges(7,1);
t16 = t295 * t50 + t299 * t30;
t542 = t16 * mrSges(7,2);
t541 = -m(5) + t496;
t186 = Ifges(6,4) * t198;
t540 = Ifges(6,2) * t198;
t539 = Ifges(6,6) * t261;
t538 = t261 * Ifges(6,5);
t420 = qJD(2) * t291;
t382 = qJD(1) * t420;
t269 = t302 * t382;
t408 = qJDD(1) * t291;
t227 = t298 * t408 + t269;
t406 = qJDD(2) * t290;
t536 = pkin(9) * t406 + qJD(3) * t393 + t227;
t134 = t261 * t299 + t295 * t316;
t135 = t261 * t295 - t299 * t316;
t189 = qJD(6) - t198;
t51 = Ifges(7,5) * t135 + Ifges(7,6) * t134 + Ifges(7,3) * t189;
t467 = Ifges(6,4) * t316;
t93 = -t467 + t539 + t540;
t534 = -t93 / 0.2e1 + t51 / 0.2e1;
t533 = -pkin(3) * t541 - t309;
t224 = qJD(2) * t388 - t301 * t406;
t281 = qJDD(2) * t293 + qJDD(3);
t104 = qJD(5) * t198 + t224 * t296 + t281 * t300;
t105 = qJD(5) * t316 + t224 * t300 - t281 * t296;
t360 = t298 * t382;
t226 = t302 * t408 - t360;
t201 = qJDD(2) * pkin(2) + t226;
t407 = qJDD(1) * t294;
t381 = t290 * t407;
t48 = t201 * t440 - t240 * t418 + t250 * t386 + t297 * t381 + t301 * t536;
t37 = -t281 * qJ(4) - t282 * qJD(4) - t48;
t26 = -pkin(4) * t224 - t37;
t14 = -pkin(5) * t105 - pkin(11) * t104 + t26;
t225 = (qJD(2) * t417 + qJDD(2) * t297) * t290;
t209 = qJDD(5) + t225;
t49 = t301 * (t201 * t293 + t381) - qJD(3) * t228 - t240 * t417 - t536 * t297;
t306 = qJDD(4) - t49;
t25 = pkin(4) * t225 - t281 * t495 + t306;
t276 = t293 * t407;
t310 = -qJ(4) * t225 + t276 + (-qJD(2) * t416 - t201) * t290;
t55 = t224 * t495 + t310;
t5 = -t117 * t415 + t296 * t25 + t300 * t55 + t75 * t414;
t3 = pkin(11) * t209 + t5;
t1 = qJD(6) * t15 + t14 * t295 + t299 * t3;
t2 = -qJD(6) * t16 + t14 * t299 - t295 * t3;
t532 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t481 = t209 / 0.2e1;
t491 = t105 / 0.2e1;
t492 = t104 / 0.2e1;
t505 = Ifges(6,1) * t492 + Ifges(6,4) * t491 + Ifges(6,5) * t481;
t45 = qJD(6) * t134 + t104 * t299 + t209 * t295;
t504 = t45 / 0.2e1;
t46 = -qJD(6) * t135 - t104 * t295 + t209 * t299;
t503 = t46 / 0.2e1;
t100 = qJDD(6) - t105;
t493 = t100 / 0.2e1;
t531 = mrSges(4,2) - mrSges(5,3);
t530 = Ifges(4,5) - Ifges(5,4);
t529 = Ifges(5,5) - Ifges(4,6);
t17 = -mrSges(7,1) * t46 + mrSges(7,2) * t45;
t69 = mrSges(6,1) * t209 - mrSges(6,3) * t104;
t472 = t17 - t69;
t348 = -mrSges(7,1) * t299 + mrSges(7,2) * t295;
t320 = m(7) * pkin(5) - t348;
t528 = mrSges(6,1) + t320;
t527 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t375 = -pkin(11) * t300 + qJ(4);
t254 = pkin(5) * t296 + t375;
t436 = t296 * t495;
t190 = t254 * t299 + t295 * t436;
t526 = qJD(6) * t190 + t295 * t547 - t299 * t546;
t191 = t254 * t295 - t299 * t436;
t525 = -qJD(6) * t191 + t295 * t546 + t299 * t547;
t470 = mrSges(6,3) * t316;
t141 = mrSges(6,1) * t261 + t470;
t67 = -mrSges(7,1) * t134 + mrSges(7,2) * t135;
t454 = t141 - t67;
t158 = mrSges(5,1) * t224 - mrSges(5,3) * t281;
t47 = -mrSges(6,1) * t105 + mrSges(6,2) * t104;
t524 = -t158 + t47;
t437 = t296 * t297;
t182 = (-t295 * t437 + t299 * t301) * t421;
t385 = t295 * t415;
t410 = qJD(6) * t300;
t522 = t299 * t410 + t182 - t385;
t435 = t297 * t299;
t183 = (t295 * t301 + t296 * t435) * t421;
t521 = t295 * t410 + t299 * t415 + t183;
t118 = -mrSges(6,1) * t198 - mrSges(6,2) * t316;
t369 = mrSges(5,1) * t389;
t215 = -mrSges(5,3) * t282 - t369;
t520 = -t215 + t118;
t368 = mrSges(4,3) * t390;
t370 = mrSges(5,1) * t390;
t427 = t282 * t548 + t368 + t370;
t361 = t300 * t390;
t519 = t361 + t414;
t20 = mrSges(7,1) * t100 - mrSges(7,3) * t45;
t21 = -mrSges(7,2) * t100 + mrSges(7,3) * t46;
t518 = -t295 * t20 + t299 * t21;
t517 = t1 * t299 - t2 * t295;
t6 = -qJD(5) * t41 + t25 * t300 - t296 * t55;
t512 = t6 * mrSges(6,1) - t5 * mrSges(6,2) + Ifges(6,5) * t104 + Ifges(6,6) * t105 + Ifges(6,3) * t209;
t36 = -qJD(5) * t523 - t143 * t296 + t175 * t300;
t349 = mrSges(6,1) * t296 + mrSges(6,2) * t300;
t510 = -m(7) * t375 + t300 * mrSges(7,3) - t296 * t320 - t349 + t531 + (-m(5) - m(6)) * qJ(4);
t304 = qJD(2) ^ 2;
t9 = Ifges(7,5) * t45 + Ifges(7,6) * t46 + Ifges(7,3) * t100;
t509 = t9 / 0.2e1;
t10 = t45 * Ifges(7,4) + t46 * Ifges(7,2) + t100 * Ifges(7,6);
t508 = t10 / 0.2e1;
t507 = Ifges(7,1) * t504 + Ifges(7,4) * t503 + Ifges(7,5) * t493;
t506 = -t104 * Ifges(6,4) / 0.2e1 - t105 * Ifges(6,2) / 0.2e1 - t209 * Ifges(6,6) / 0.2e1;
t464 = Ifges(7,4) * t135;
t52 = t134 * Ifges(7,2) + t189 * Ifges(7,6) + t464;
t501 = -t52 / 0.2e1;
t500 = t52 / 0.2e1;
t131 = Ifges(7,4) * t134;
t53 = t135 * Ifges(7,1) + t189 * Ifges(7,5) + t131;
t499 = -t53 / 0.2e1;
t498 = t53 / 0.2e1;
t490 = -t134 / 0.2e1;
t489 = t134 / 0.2e1;
t488 = -t135 / 0.2e1;
t487 = t135 / 0.2e1;
t486 = -t189 / 0.2e1;
t485 = t189 / 0.2e1;
t482 = -t316 / 0.2e1;
t478 = t296 * t5;
t4 = -pkin(5) * t209 - t6;
t477 = t300 * t4;
t474 = -mrSges(4,3) - mrSges(5,1);
t473 = Ifges(4,4) + Ifges(5,6);
t471 = mrSges(6,3) * t198;
t469 = mrSges(6,3) * t300;
t468 = Ifges(4,4) * t297;
t466 = Ifges(6,4) * t296;
t465 = Ifges(6,4) * t300;
t463 = Ifges(7,4) * t295;
t462 = Ifges(7,4) * t299;
t461 = Ifges(5,6) * t297;
t460 = Ifges(5,6) * t301;
t459 = t198 * Ifges(6,6);
t458 = t316 * Ifges(6,5);
t457 = t261 * Ifges(6,3);
t453 = sin(pkin(12));
t451 = t198 * t295;
t450 = t198 * t299;
t449 = t290 * t296;
t447 = t290 * t300;
t292 = cos(pkin(12));
t445 = t291 * t292;
t443 = t291 * t302;
t442 = t292 * t298;
t441 = t292 * t302;
t438 = t295 * t300;
t432 = t299 * t300;
t217 = (mrSges(5,2) * t301 - mrSges(5,3) * t297) * t421;
t426 = t217 + (-mrSges(4,1) * t301 + mrSges(4,2) * t297) * t421;
t400 = t290 * t444;
t424 = pkin(2) * t443 + pkin(9) * t400;
t419 = t304 * t290 ^ 2;
t412 = qJD(6) * t295;
t411 = qJD(6) * t299;
t405 = Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1;
t404 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t403 = -m(4) + t541;
t401 = t290 * t445;
t397 = t301 * t443;
t197 = -t366 + t397;
t396 = t197 * pkin(3) + t424;
t367 = mrSges(4,3) * t389;
t214 = -mrSges(4,2) * t282 + t367;
t395 = t214 + t520;
t384 = t296 * t413;
t377 = -t415 / 0.2e1;
t376 = -t410 / 0.2e1;
t374 = t291 * t453;
t373 = t453 * t298;
t372 = t453 * t302;
t159 = t225 * mrSges(5,1) + t281 * mrSges(5,2);
t362 = qJD(2) * t400;
t353 = t290 * t374;
t120 = t240 * t297 - t425;
t346 = Ifges(6,1) * t296 + t465;
t345 = Ifges(7,1) * t299 - t463;
t344 = Ifges(7,1) * t295 + t462;
t343 = Ifges(6,2) * t300 + t466;
t342 = -Ifges(7,2) * t295 + t462;
t341 = Ifges(7,2) * t299 + t463;
t340 = Ifges(6,5) * t296 + Ifges(6,6) * t300;
t339 = Ifges(7,5) * t299 - Ifges(7,6) * t295;
t338 = Ifges(7,5) * t295 + Ifges(7,6) * t299;
t337 = -pkin(3) * t301 - t452;
t40 = -t117 * t296 + t300 * t75;
t334 = t296 * t40 - t300 * t41;
t65 = t108 * t300 - t174 * t296;
t149 = -t293 * t397 - t294 * t446 + t398;
t231 = -t290 * t443 + t293 * t294;
t114 = t149 * t296 + t231 * t300;
t326 = t293 * t434 + t433;
t150 = t291 * t326 + t294 * t448;
t59 = t114 * t299 + t150 * t295;
t58 = -t114 * t295 + t150 * t299;
t84 = t148 * t300 - t166 * t296;
t331 = t149 * t300 - t231 * t296;
t220 = -pkin(9) * t388 + t275;
t329 = qJ(4) * t541 + t531;
t162 = -t237 * t295 + t290 * t435;
t163 = t237 * t299 + t295 * t448;
t324 = t297 * (Ifges(4,1) * t301 - t468);
t323 = t297 * (-Ifges(5,2) * t301 + t461);
t322 = t301 * (Ifges(5,3) * t297 - t460);
t232 = t294 * t441 - t373;
t233 = t294 * t442 + t372;
t109 = -t232 * t439 + t233 * t297 + t301 * t401;
t234 = -t294 * t372 - t442;
t235 = -t294 * t373 + t441;
t111 = -t234 * t439 + t235 * t297 - t301 * t353;
t318 = -g(1) * t111 - g(2) * t109 - g(3) * t149;
t42 = -pkin(3) * t281 + t306;
t312 = t49 * mrSges(4,1) - t48 * mrSges(4,2) + t42 * mrSges(5,2) - t37 * mrSges(5,3);
t308 = (-t15 * t299 - t16 * t295) * qJD(6) + t517;
t307 = -qJD(5) * t334 + t300 * t6 + t478;
t271 = Ifges(4,4) * t389;
t268 = Ifges(5,1) * t281;
t267 = Ifges(4,3) * t281;
t238 = pkin(2) * t439 - t283;
t223 = t234 * pkin(2);
t222 = t232 * pkin(2);
t221 = t239 * qJD(3);
t219 = -qJ(4) * t389 + t272;
t208 = Ifges(5,4) * t225;
t207 = Ifges(4,5) * t225;
t206 = Ifges(5,5) * t224;
t205 = Ifges(4,6) * t224;
t204 = t293 * t394 + t283;
t203 = (-pkin(2) + t337) * t290;
t188 = -t220 - t288;
t185 = -t250 * t290 + t277;
t181 = t282 * t295 - t299 * t361;
t180 = t282 * t299 + t295 * t361;
t176 = t274 + (-qJ(4) * t417 - t416) * t290;
t170 = t282 * Ifges(5,4) + (-Ifges(5,2) * t297 - t460) * t421;
t169 = t282 * Ifges(5,5) + (-t301 * Ifges(5,3) - t461) * t421;
t168 = Ifges(4,1) * t390 + t282 * Ifges(4,5) + t271;
t167 = t282 * Ifges(4,6) + (t301 * Ifges(4,2) + t468) * t421;
t157 = mrSges(4,1) * t281 - mrSges(4,3) * t225;
t156 = -mrSges(4,2) * t281 - mrSges(4,3) * t224;
t153 = -t234 * t290 + t293 * t374;
t152 = -t232 * t290 - t293 * t445;
t145 = -t201 * t290 + t276;
t140 = -mrSges(6,2) * t261 + t471;
t133 = -mrSges(5,2) * t224 - mrSges(5,3) * t225;
t132 = mrSges(4,1) * t224 + mrSges(4,2) * t225;
t130 = t277 + (qJD(2) * t337 - t250) * t290;
t128 = -t300 * t178 + t296 * t363;
t127 = t234 * t301 - t235 * t440;
t126 = t234 * t297 + t235 * t439;
t125 = t232 * t301 - t233 * t440;
t124 = t232 * t297 + t233 * t439;
t119 = -pkin(5) * t316 - pkin(11) * t198;
t112 = t235 * t301 + (t234 * t293 + t353) * t297;
t110 = t232 * t440 + t233 * t301 - t297 * t401;
t106 = -t270 - t121;
t101 = -pkin(3) * t282 + qJD(4) + t120;
t96 = -qJD(2) * t366 - qJD(3) * t398 + (t302 * t420 + (t290 * t294 + t293 * t443) * qJD(3)) * t301;
t95 = t294 * t388 + (qJD(2) * t327 + qJD(3) * t326) * t291;
t94 = -Ifges(6,1) * t316 + t186 + t538;
t92 = t457 - t458 + t459;
t87 = mrSges(7,1) * t189 - mrSges(7,3) * t135;
t86 = -mrSges(7,2) * t189 + mrSges(7,3) * t134;
t78 = -pkin(5) * t448 - t84;
t74 = -qJD(6) * t163 - t160 * t295 + t299 * t387;
t73 = qJD(6) * t162 + t160 * t299 + t295 * t387;
t70 = -mrSges(6,2) * t209 + mrSges(6,3) * t105;
t68 = pkin(3) * t224 + t310;
t64 = t111 * t296 + t153 * t300;
t62 = t109 * t296 + t152 * t300;
t56 = -pkin(5) * t389 - t65;
t39 = qJD(5) * t331 + t296 * t95 + t300 * t362;
t38 = qJD(5) * t114 + t296 * t362 - t95 * t300;
t29 = -pkin(5) * t261 - t40;
t28 = -pkin(5) * t387 - t36;
t19 = t119 * t295 + t299 * t40;
t18 = t119 * t299 - t295 * t40;
t13 = qJD(6) * t58 + t295 * t96 + t299 * t39;
t12 = -qJD(6) * t59 - t295 * t39 + t299 * t96;
t7 = [m(2) * qJDD(1) + t114 * t70 + t12 * t87 + t13 * t86 + t39 * t140 + t58 * t20 + t59 * t21 + t427 * t95 - t454 * t38 + (t132 + t133) * t231 + (-t157 + t159) * t149 - t472 * t331 + t395 * t96 + (t156 + t524) * t150 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t304) * t302 + (-mrSges(3,1) * t304 - mrSges(3,2) * qJDD(2) + t421 * t426) * t298) * t291 + (-m(2) - m(3) + t403) * g(3) + m(3) * (qJDD(1) * t294 ^ 2 + (t226 * t302 + t227 * t298) * t291) + m(5) * (t101 * t95 - t106 * t96 + t130 * t362 + t149 * t42 - t150 * t37 + t231 * t68) + m(4) * (t120 * t95 + t121 * t96 + t145 * t231 - t149 * t49 + t150 * t48 + t185 * t362) + m(7) * (t1 * t59 + t12 * t15 + t13 * t16 + t2 * t58 + t29 * t38 - t331 * t4) + m(6) * (t114 * t5 + t150 * t26 + t331 * t6 - t38 * t40 + t39 * t41 + t88 * t96); (mrSges(6,1) * t26 - mrSges(6,3) * t5 - Ifges(6,4) * t492 + Ifges(7,5) * t504 - Ifges(6,2) * t491 - Ifges(6,6) * t481 + Ifges(7,6) * t503 + Ifges(7,3) * t493 + t506 + t509 + t532) * t236 + (mrSges(6,2) * t26 - t6 * mrSges(6,3) + 0.2e1 * t505) * t237 + (t1 * t162 - t15 * t73 + t16 * t74 - t163 * t2) * mrSges(7,3) + (t165 * t26 + t5 * t523 + t6 * t84 + t552 * t88 + t537 * t41 + (t128 + t36) * t40) * m(6) + (Ifges(7,5) * t73 + Ifges(7,6) * t74) * t485 + (Ifges(7,5) * t163 + Ifges(7,6) * t162) * t493 + t523 * t70 + (Ifges(7,4) * t163 + Ifges(7,2) * t162) * t503 + (Ifges(7,4) * t73 + Ifges(7,2) * t74) * t489 + (Ifges(7,1) * t163 + Ifges(7,4) * t162) * t504 + (Ifges(7,1) * t73 + Ifges(7,4) * t74) * t487 + (-m(5) * t396 - m(4) * t424 + (-mrSges(3,1) * t302 + mrSges(3,2) * t298) * t291 + t329 * t196 + t527 * (-t196 * t300 + t296 * t400) - t528 * (t196 * t296 + t300 * t400) + t309 * t197 + t496 * (pkin(4) * t400 + t396)) * g(3) + (-m(4) * t223 - mrSges(3,1) * t234 + mrSges(3,2) * t235 - t527 * (t126 * t300 - t235 * t449) + t329 * t126 - t528 * (t126 * t296 + t235 * t447) + t309 * t127 + t541 * (t127 * pkin(3) + t223)) * g(1) + (-m(4) * t222 - mrSges(3,1) * t232 + mrSges(3,2) * t233 - t527 * (t124 * t300 - t233 * t449) + t329 * t124 - t528 * (t124 * t296 + t233 * t447) + t309 * t125 + t541 * (t125 * pkin(3) + t222)) * g(2) + (-t227 + t269) * mrSges(3,2) + (t226 + t360) * mrSges(3,1) + m(5) * (t101 * t221 + t106 * t188 + t130 * t176 + t202 * t37 + t203 * t68 + t204 * t42) + (t221 - t178) * t427 + t454 * t128 + Ifges(3,3) * qJDD(2) + t163 * t507 + t162 * t508 + t73 * t498 + t74 * t500 + (t207 / 0.2e1 - t205 / 0.2e1 + t267 / 0.2e1 + t268 / 0.2e1 - t208 / 0.2e1 + t206 / 0.2e1 + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t281 + t405 * t225 + t404 * t224 + t312) * t293 - t395 * t179 + t238 * t157 + t239 * t156 - m(5) * (t101 * t178 - t106 * t179) - m(4) * (t120 * t178 + t121 * t179) + m(4) * (t120 * t221 + t121 * t220 + t238 * t49 + t239 * t48) + t176 * t217 + t220 * t214 + t188 * t215 + t202 * t158 + t203 * t133 + t204 * t159 + t4 * (-mrSges(7,1) * t162 + mrSges(7,2) * t163) + t165 * t47 + t36 * t141 + t147 * t118 + t84 * t69 + t78 * t17 + t29 * (-mrSges(7,1) * t74 + mrSges(7,2) * t73) + t28 * t67 + t31 * t20 + t32 * t21 + t537 * t140 + (t186 / 0.2e1 - t40 * mrSges(6,3) + t538 / 0.2e1 + Ifges(6,1) * t482 + t88 * mrSges(6,2) + t94 / 0.2e1) * t160 + (-t540 / 0.2e1 - t41 * mrSges(6,3) - t539 / 0.2e1 + Ifges(7,3) * t485 + Ifges(7,5) * t487 + Ifges(7,6) * t489 - Ifges(6,4) * t482 + t88 * mrSges(6,1) + t543 - t542 + t534) * t161 + t544 * t86 + t545 * t87 + (t1 * t32 + t2 * t31 + t4 * t78 + (-t128 + t28) * t29 + t544 * t16 + t545 * t15) * m(7) + ((-m(4) * t145 - t132) * pkin(2) + (t474 * g(3) + (-m(4) * t185 - m(5) * t130 - t426) * qJD(1)) * t444 + (-t145 * mrSges(4,1) - t37 * mrSges(5,1) + t68 * mrSges(5,2) + t48 * mrSges(4,3) - t529 * t281 + t473 * t225 + (-Ifges(4,2) - Ifges(5,3)) * t224) * t301 + (t42 * mrSges(5,1) - t49 * mrSges(4,3) + t145 * mrSges(4,2) - t68 * mrSges(5,3) + t530 * t281 + (Ifges(4,1) + Ifges(5,2)) * t225 - t473 * t224 + t512) * t297 + ((-t167 / 0.2e1 + t169 / 0.2e1 - t121 * mrSges(4,3) + t106 * mrSges(5,1) - t130 * mrSges(5,2) + t185 * mrSges(4,1) + t404 * t282) * t297 + (t168 / 0.2e1 - t170 / 0.2e1 + t92 / 0.2e1 + t120 * mrSges(4,3) + t101 * mrSges(5,1) + t40 * mrSges(6,1) + t459 / 0.2e1 - t458 / 0.2e1 - t41 * mrSges(6,2) + t457 / 0.2e1 - t130 * mrSges(5,3) + t185 * mrSges(4,2) + t405 * t282) * t301 + (-t323 / 0.2e1 - t322 / 0.2e1 + t324 / 0.2e1 + t301 * (Ifges(4,4) * t301 - Ifges(4,2) * t297) / 0.2e1) * t421) * qJD(3) + (g(1) * t235 + g(2) * t233) * (pkin(4) * t496 + pkin(9) * t403 + t474)) * t290; t472 * t300 * t495 + t94 * t377 + ((t300 * t51 + t167) * t297 + t301 * t170) * t421 / 0.2e1 - ((t296 * t94 + t300 * t93 + t169) * t297 + (-Ifges(4,2) * t390 + t168 + t271 + t92) * t301 + t261 * (Ifges(6,3) * t301 + t297 * t340) - t316 * (Ifges(6,5) * t301 + t297 * t346) + t198 * (Ifges(6,6) * t301 + t297 * t343) + (t297 * t529 + t301 * t530) * t282) * t421 / 0.2e1 - (t198 * t343 + t261 * t340 - t316 * t346) * qJD(5) / 0.2e1 + (qJ(4) * t26 - t307 * t495 - t40 * t65 - t41 * t66 + t516 * t88) * m(6) + (-t29 * t56 + t1 * t191 + t190 * t2 - (t29 * t415 - t477) * t495 + t526 * t16 + t525 * t15) * m(7) + (-m(5) * t106 + t214 - t215 - t367) * t120 + (t111 * t533 + t112 * t510) * g(1) + (t149 * t533 + t150 * t510) * g(3) + (t109 * t533 + t110 * t510) * g(2) + t534 * t414 + t261 * t88 * (mrSges(6,1) * t300 - mrSges(6,2) * t296) + t267 + t268 + (t384 - t65) * t141 + (-t384 - t56) * t67 + (t40 * t415 - t41 * t414 - t478) * mrSges(6,3) + t312 + (-pkin(3) * t42 - qJ(4) * t37 - qJD(4) * t106 - t130 * t219) * m(5) - t208 - t205 + t206 + t207 - t6 * t469 + (t323 + t322) * t419 / 0.2e1 + (t376 * t52 + t377 * t53) * t299 + t295 * t53 * t376 + t549 * t140 - t10 * t438 / 0.2e1 + t1 * (-mrSges(7,2) * t296 - mrSges(7,3) * t438) + (Ifges(7,5) * t296 + t300 * t345) * t504 + t300 * t505 + t296 * t506 + t432 * t507 + t296 * t509 + t183 * t499 + t385 * t500 + t182 * t501 + (Ifges(7,6) * t296 + t300 * t342) * t503 + (-t338 * t410 + (Ifges(7,3) * t300 - t296 * t339) * qJD(5)) * t485 + (Ifges(7,5) * t183 + Ifges(7,6) * t182 - Ifges(7,3) * t361) * t486 + (-t344 * t410 + (Ifges(7,5) * t300 - t296 * t345) * qJD(5)) * t487 + (Ifges(7,1) * t183 + Ifges(7,4) * t182 - Ifges(7,5) * t361) * t488 + (-t341 * t410 + (Ifges(7,6) * t300 - t296 * t342) * qJD(5)) * t489 + (Ifges(7,4) * t183 + Ifges(7,2) * t182 - Ifges(7,6) * t361) * t490 + (-Ifges(6,2) * t296 + t465) * t491 + (Ifges(6,1) * t300 - t466) * t492 + (Ifges(7,3) * t296 + t300 * t339) * t493 + t347 * t477 + (Ifges(6,5) * t300 - Ifges(6,6) * t296) * t481 + t2 * (mrSges(7,1) * t296 - mrSges(7,3) * t432) - t324 * t419 / 0.2e1 - t219 * t217 + t190 * t20 + t191 * t21 - pkin(3) * t159 - t107 * t118 - t106 * t370 - t101 * t369 - t70 * t436 + (-m(5) * t101 + t368 - t427) * t121 + t520 * qJD(4) + (mrSges(7,1) * t519 + mrSges(7,3) * t521) * t15 + (mrSges(7,1) * t522 - mrSges(7,2) * t521) * t29 + (-mrSges(7,2) * t519 - mrSges(7,3) * t522) * t16 + t26 * t349 + t524 * qJ(4) + t525 * t87 + t526 * t86 + (-t185 * (mrSges(4,1) * t297 + mrSges(4,2) * t301) - t130 * (-mrSges(5,2) * t297 - mrSges(5,3) * t301) - t41 * (-mrSges(6,2) * t301 + t297 * t469) - t40 * (mrSges(6,1) * t301 - mrSges(6,3) * t437)) * t421; t217 * t390 - t180 * t87 - t181 * t86 - t520 * t282 + (t140 * t390 + (-t295 * t87 + t299 * t86 + t140) * qJD(5) - t472) * t300 + (t70 + (-t295 * t86 - t299 * t87) * qJD(6) - t261 * t454 + t518) * t296 + t159 + (t318 + (-t4 + (-t15 * t295 + t16 * t299) * qJD(5)) * t300 - t15 * t180 - t16 * t181 + (t261 * t29 + t308) * t296) * m(7) + (-t282 * t88 - t334 * t390 + t307 + t318) * m(6) + (t106 * t282 + t130 * t390 + t318 + t42) * m(5); (-pkin(5) * t4 - t15 * t18 - t16 * t19) * m(7) + (t114 * t527 - t331 * t528) * g(3) + (Ifges(6,1) * t198 + t467 + t51) * t316 / 0.2e1 - (Ifges(6,2) * t316 + t186 + t94) * t198 / 0.2e1 + (-Ifges(7,3) * t316 + t198 * t339) * t486 + (-Ifges(7,5) * t316 + t198 * t345) * t488 + (-Ifges(7,6) * t316 + t198 * t342) * t490 - t261 * (Ifges(6,5) * t198 + Ifges(6,6) * t316) / 0.2e1 - t88 * (-mrSges(6,1) * t316 + mrSges(6,2) * t198) + t189 * t29 * t347 + (t471 - t140) * t40 + t512 + (t134 * t342 + t135 * t345 + t189 * t339) * qJD(6) / 0.2e1 + t344 * t504 + t295 * t507 + t299 * t508 + t411 * t498 + t450 * t499 + t451 * t500 + t412 * t501 + t341 * t503 + t338 * t493 + t93 * t482 - t19 * t86 - t18 * t87 - pkin(5) * t17 + t316 * t543 - t316 * t542 + ((-t412 + t451) * t16 + (-t411 + t450) * t15 + t517) * mrSges(7,3) + (m(7) * t308 - t411 * t87 - t412 * t86 + t518) * pkin(11) + t4 * t348 + (-m(7) * t29 + t454 - t470) * t41 + (t527 * t62 - t528 * (t109 * t300 - t152 * t296)) * g(2) + (t527 * t64 - t528 * (t111 * t300 - t153 * t296)) * g(1); -t29 * (mrSges(7,1) * t135 + mrSges(7,2) * t134) + (Ifges(7,1) * t134 - t464) * t488 + t52 * t487 + (Ifges(7,5) * t134 - Ifges(7,6) * t135) * t486 - t15 * t86 + t16 * t87 - g(1) * ((t112 * t299 - t295 * t64) * mrSges(7,1) + (-t112 * t295 - t299 * t64) * mrSges(7,2)) - g(2) * ((t110 * t299 - t295 * t62) * mrSges(7,1) + (-t110 * t295 - t299 * t62) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t58 - mrSges(7,2) * t59) + (t134 * t15 + t135 * t16) * mrSges(7,3) + t9 + (-Ifges(7,2) * t135 + t131 + t53) * t490 + t532;];
tau  = t7;
