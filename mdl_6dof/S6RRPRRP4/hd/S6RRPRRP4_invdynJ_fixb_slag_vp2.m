% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:52:35
% EndTime: 2019-03-09 11:53:29
% DurationCPUTime: 32.17s
% Computational Cost: add. (15634->811), mult. (35196->1029), div. (0->0), fcn. (25929->14), ass. (0->378)
t568 = mrSges(6,1) + mrSges(7,1);
t566 = mrSges(6,2) - mrSges(7,3);
t540 = -Ifges(6,4) + Ifges(7,5);
t577 = t540 + Ifges(7,5);
t293 = sin(qJ(2));
t297 = cos(qJ(2));
t261 = -mrSges(3,1) * t297 + mrSges(3,2) * t293;
t288 = qJ(2) + pkin(10);
t283 = sin(t288);
t284 = cos(t288);
t575 = -mrSges(6,3) - mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t576 = -t284 * mrSges(4,1) + t283 * t575 + t261;
t377 = qJD(1) * qJD(2);
t357 = t293 * t377;
t376 = qJDD(1) * t297;
t250 = -t357 + t376;
t251 = qJDD(1) * t293 + t297 * t377;
t417 = sin(pkin(10));
t418 = cos(pkin(10));
t180 = t250 * t417 + t251 * t418;
t292 = sin(qJ(4));
t296 = cos(qJ(4));
t347 = qJD(1) * t417;
t348 = qJD(1) * t418;
t228 = -t293 * t348 - t297 * t347;
t318 = qJD(2) * t296 + t228 * t292;
t308 = t318 * qJD(4);
t119 = qJDD(2) * t292 + t180 * t296 + t308;
t189 = qJD(2) * t292 - t228 * t296;
t120 = -qJD(4) * t189 + qJDD(2) * t296 - t180 * t292;
t291 = sin(qJ(5));
t295 = cos(qJ(5));
t125 = t189 * t291 - t295 * t318;
t50 = -qJD(5) * t125 + t295 * t119 + t291 * t120;
t492 = t50 / 0.2e1;
t304 = t295 * t189 + t291 * t318;
t51 = qJD(5) * t304 + t291 * t119 - t295 * t120;
t490 = t51 / 0.2e1;
t179 = t250 * t418 - t417 * t251;
t176 = qJDD(4) - t179;
t173 = qJDD(5) + t176;
t474 = t173 / 0.2e1;
t361 = t417 * pkin(2);
t274 = t361 + pkin(8);
t442 = pkin(9) + t274;
t352 = qJD(4) * t442;
t227 = -t293 * t347 + t297 * t348;
t415 = t227 * t292;
t386 = qJD(1) * t293;
t372 = pkin(2) * t386;
t145 = -pkin(3) * t228 - pkin(8) * t227 + t372;
t290 = -qJ(3) - pkin(7);
t262 = t290 * t297;
t249 = qJD(1) * t262;
t231 = t417 * t249;
t260 = t290 * t293;
t248 = qJD(1) * t260;
t178 = t248 * t418 + t231;
t99 = t292 * t145 + t296 * t178;
t574 = -pkin(9) * t415 + t292 * t352 + t99;
t414 = t227 * t296;
t98 = t296 * t145 - t178 * t292;
t573 = pkin(4) * t228 + pkin(9) * t414 - t296 * t352 - t98;
t381 = qJD(4) * t292;
t572 = t381 - t415;
t237 = qJD(2) * pkin(2) + t248;
t170 = t237 * t418 + t231;
t156 = -qJD(2) * pkin(3) - t170;
t121 = -pkin(4) * t318 + t156;
t218 = qJD(4) - t227;
t212 = qJD(5) + t218;
t287 = t297 * pkin(2);
t280 = t287 + pkin(1);
t255 = -qJD(1) * t280 + qJD(3);
t134 = -pkin(3) * t227 + pkin(8) * t228 + t255;
t351 = t418 * t249;
t171 = t417 * t237 - t351;
t157 = qJD(2) * pkin(8) + t171;
t90 = t292 * t134 + t296 * t157;
t78 = pkin(9) * t318 + t90;
t421 = t291 * t78;
t89 = t296 * t134 - t157 * t292;
t77 = -pkin(9) * t189 + t89;
t70 = pkin(4) * t218 + t77;
t27 = t295 * t70 - t421;
t557 = qJD(6) - t27;
t23 = -pkin(5) * t212 + t557;
t54 = t125 * pkin(5) - qJ(6) * t304 + t121;
t123 = Ifges(6,4) * t125;
t427 = Ifges(7,5) * t125;
t539 = Ifges(7,4) + Ifges(6,5);
t541 = Ifges(6,1) + Ifges(7,1);
t533 = t212 * t539 + t304 * t541 - t123 + t427;
t560 = t533 / 0.2e1;
t571 = t121 * mrSges(6,2) + mrSges(7,2) * t23 - mrSges(6,3) * t27 - t54 * mrSges(7,3) + t560;
t570 = m(3) * pkin(1) + mrSges(2,1) - t576;
t569 = m(5) + m(4);
t548 = -m(6) - m(7);
t538 = -Ifges(6,6) + Ifges(7,6);
t537 = Ifges(6,3) + Ifges(7,2);
t246 = t291 * t292 - t295 * t296;
t141 = t246 * t227;
t512 = qJD(4) + qJD(5);
t182 = t512 * t246;
t522 = -t182 + t141;
t401 = t291 * t296;
t247 = t292 * t295 + t401;
t140 = t247 * t227;
t183 = t512 * t247;
t521 = t183 - t140;
t468 = t212 / 0.2e1;
t477 = t304 / 0.2e1;
t480 = t125 / 0.2e1;
t481 = -t125 / 0.2e1;
t564 = Ifges(6,4) * t481 + Ifges(7,5) * t480 + t468 * t539 + t477 * t541 + t571;
t420 = t295 * t78;
t28 = t291 * t70 + t420;
t24 = qJ(6) * t212 + t28;
t122 = Ifges(7,5) * t304;
t63 = t212 * Ifges(7,6) + t125 * Ifges(7,3) + t122;
t428 = Ifges(6,4) * t304;
t66 = -t125 * Ifges(6,2) + t212 * Ifges(6,6) + t428;
t563 = -mrSges(7,2) * t24 - mrSges(6,3) * t28 + t121 * mrSges(6,1) + t54 * mrSges(7,1) + t63 / 0.2e1 - t66 / 0.2e1;
t242 = t251 * pkin(7);
t383 = qJD(3) * t293;
t162 = qJDD(2) * pkin(2) - qJ(3) * t251 - qJD(1) * t383 - t242;
t281 = pkin(7) * t376;
t384 = qJD(2) * t293;
t367 = pkin(7) * t384;
t382 = qJD(3) * t297;
t174 = qJ(3) * t250 + t281 + (-t367 + t382) * qJD(1);
t109 = t162 * t418 - t417 * t174;
t101 = -qJDD(2) * pkin(3) - t109;
t59 = -t120 * pkin(4) + t101;
t11 = t51 * pkin(5) - t50 * qJ(6) - qJD(6) * t304 + t59;
t110 = t417 * t162 + t418 * t174;
t102 = qJDD(2) * pkin(8) + t110;
t416 = qJDD(1) * pkin(1);
t215 = -pkin(2) * t250 + qJDD(3) - t416;
t95 = -pkin(3) * t179 - pkin(8) * t180 + t215;
t26 = -t90 * qJD(4) - t102 * t292 + t296 * t95;
t20 = pkin(4) * t176 - pkin(9) * t119 + t26;
t380 = qJD(4) * t296;
t25 = t296 * t102 + t134 * t380 - t157 * t381 + t292 * t95;
t22 = pkin(9) * t120 + t25;
t6 = -qJD(5) * t28 + t20 * t295 - t22 * t291;
t3 = -pkin(5) * t173 + qJDD(6) - t6;
t491 = -t51 / 0.2e1;
t562 = mrSges(6,2) * t59 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t11 + Ifges(6,4) * t491 + 0.2e1 * t474 * t539 + t490 * t577 + 0.2e1 * t492 * t541;
t289 = qJ(4) + qJ(5);
t285 = sin(t289);
t286 = cos(t289);
t335 = -mrSges(5,1) * t296 + mrSges(5,2) * t292;
t561 = m(5) * pkin(3) - t566 * t285 + t286 * t568 - t335;
t558 = -m(7) * qJ(6) - mrSges(7,3);
t238 = t442 * t292;
t239 = t442 * t296;
t167 = -t238 * t291 + t239 * t295;
t529 = -qJD(5) * t167 + t291 * t574 + t295 * t573;
t319 = -t295 * t238 - t239 * t291;
t525 = qJD(5) * t319 + t291 * t573 - t295 * t574;
t435 = mrSges(6,3) * t304;
t105 = mrSges(6,1) * t212 - t435;
t106 = -mrSges(7,1) * t212 + mrSges(7,2) * t304;
t523 = t106 - t105;
t445 = t296 * pkin(4);
t279 = pkin(3) + t445;
t253 = t284 * t279;
t299 = -pkin(9) - pkin(8);
t405 = t283 * t299;
t556 = t253 - t405;
t244 = t293 * t418 + t297 * t417;
t359 = t244 * t380;
t243 = t293 * t417 - t297 * t418;
t230 = t243 * qJD(2);
t413 = t230 * t292;
t315 = t359 - t413;
t294 = sin(qJ(1));
t396 = t294 * t296;
t298 = cos(qJ(1));
t398 = t292 * t298;
t221 = -t284 * t398 + t396;
t177 = t248 * t417 - t351;
t518 = pkin(4) * t572 - t177;
t553 = -Ifges(6,2) * t481 + Ifges(7,3) * t480 + t468 * t538 + t477 * t540 + t563;
t469 = -t212 / 0.2e1;
t478 = -t304 / 0.2e1;
t552 = -Ifges(6,2) * t480 + Ifges(7,3) * t481 + t469 * t538 + t478 * t540 - t563;
t79 = pkin(5) * t304 + qJ(6) * t125;
t483 = t119 / 0.2e1;
t482 = t120 / 0.2e1;
t473 = t176 / 0.2e1;
t547 = t250 / 0.2e1;
t546 = -t318 / 0.2e1;
t543 = t89 * mrSges(5,1);
t542 = t90 * mrSges(5,2);
t444 = qJD(2) / 0.2e1;
t33 = mrSges(6,1) * t173 - mrSges(6,3) * t50;
t34 = -t173 * mrSges(7,1) + t50 * mrSges(7,2);
t535 = t34 - t33;
t35 = -mrSges(6,2) * t173 - mrSges(6,3) * t51;
t36 = -mrSges(7,2) * t51 + mrSges(7,3) * t173;
t534 = t36 + t35;
t532 = pkin(5) * t521 - qJ(6) * t522 - qJD(6) * t247 + t518;
t531 = Ifges(5,3) * t218;
t530 = t318 * Ifges(5,6);
t419 = qJDD(2) / 0.2e1;
t528 = -pkin(5) * t228 - t529;
t527 = qJ(6) * t228 + t525;
t151 = t247 * t244;
t103 = -mrSges(7,2) * t125 + mrSges(7,3) * t212;
t436 = mrSges(6,3) * t125;
t104 = -mrSges(6,2) * t212 - t436;
t524 = t104 + t103;
t169 = pkin(3) * t243 - pkin(8) * t244 - t280;
t187 = t260 * t417 - t262 * t418;
t181 = t296 * t187;
t116 = t292 * t169 + t181;
t438 = mrSges(4,3) * t228;
t520 = qJD(2) * mrSges(4,1) + mrSges(5,1) * t318 - t189 * mrSges(5,2) + t438;
t517 = t173 * t537 + t50 * t539 + t51 * t538;
t393 = t298 * t285;
t397 = t294 * t286;
t204 = t284 * t393 - t397;
t403 = t286 * t298;
t404 = t285 * t294;
t205 = t284 * t403 + t404;
t516 = t204 * t568 + t205 * t566;
t202 = t284 * t404 + t403;
t203 = t284 * t397 - t393;
t515 = t202 * t568 + t203 * t566;
t241 = -pkin(7) * t357 + t281;
t514 = t241 * t297 + t242 * t293;
t513 = t25 * t296 - t26 * t292;
t511 = Ifges(5,5) * t189 + t125 * t538 + t212 * t537 + t304 * t539 + t530 + t531;
t510 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t509 = 0.2e1 * t419;
t229 = t244 * qJD(2);
t353 = qJD(2) * t290;
t225 = t293 * t353 + t382;
t226 = t297 * t353 - t383;
t144 = t225 * t418 + t226 * t417;
t371 = pkin(2) * t384;
t146 = pkin(3) * t229 + pkin(8) * t230 + t371;
t344 = -t144 * t292 + t296 * t146;
t412 = t230 * t296;
t40 = pkin(9) * t412 + pkin(4) * t229 + (-t181 + (pkin(9) * t244 - t169) * t292) * qJD(4) + t344;
t115 = t296 * t169 - t187 * t292;
t408 = t244 * t296;
t87 = pkin(4) * t243 - pkin(9) * t408 + t115;
t409 = t244 * t292;
t94 = -pkin(9) * t409 + t116;
t441 = t291 * t87 + t295 * t94;
t57 = t296 * t144 + t292 * t146 + t169 * t380 - t187 * t381;
t49 = -pkin(9) * t315 + t57;
t10 = -qJD(5) * t441 - t291 * t49 + t295 * t40;
t508 = m(7) * pkin(5) + t568;
t506 = -mrSges(6,2) - t558;
t505 = t26 * mrSges(5,1) - t25 * mrSges(5,2);
t378 = qJD(5) * t295;
t379 = qJD(5) * t291;
t5 = t291 * t20 + t295 * t22 + t70 * t378 - t379 * t78;
t2 = qJ(6) * t173 + qJD(6) * t212 + t5;
t501 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3);
t496 = mrSges(6,1) * t59 + mrSges(7,1) * t11 - mrSges(7,2) * t2 - mrSges(6,3) * t5 + 0.2e1 * Ifges(7,3) * t490 - t50 * Ifges(6,4) / 0.2e1 - t173 * Ifges(6,6) / 0.2e1 + t577 * t492 + (t538 + Ifges(7,6)) * t474 + (-t491 + t490) * Ifges(6,2);
t489 = Ifges(5,1) * t483 + Ifges(5,4) * t482 + Ifges(5,5) * t473;
t471 = -t189 / 0.2e1;
t470 = t189 / 0.2e1;
t467 = -t218 / 0.2e1;
t465 = t227 / 0.2e1;
t464 = -t228 / 0.2e1;
t455 = pkin(2) * t293;
t454 = pkin(4) * t189;
t453 = pkin(4) * t291;
t452 = pkin(4) * t292;
t451 = pkin(4) * t295;
t450 = pkin(7) * t297;
t449 = pkin(8) * t283;
t446 = g(3) * t283;
t440 = mrSges(6,2) * t286;
t439 = mrSges(4,3) * t227;
t437 = mrSges(5,3) * t189;
t434 = Ifges(3,4) * t293;
t433 = Ifges(3,4) * t297;
t432 = Ifges(4,4) * t228;
t431 = Ifges(5,4) * t189;
t430 = Ifges(5,4) * t292;
t429 = Ifges(5,4) * t296;
t402 = t290 * t298;
t112 = Ifges(5,2) * t318 + Ifges(5,6) * t218 + t431;
t400 = t292 * t112;
t399 = t292 * t294;
t188 = Ifges(5,4) * t318;
t113 = t189 * Ifges(5,1) + t218 * Ifges(5,5) + t188;
t395 = t296 * t113;
t394 = t296 * t298;
t385 = qJD(1) * t297;
t375 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t386) * t450;
t369 = pkin(4) * t379;
t368 = pkin(4) * t378;
t363 = Ifges(5,5) * t119 + Ifges(5,6) * t120 + Ifges(5,3) * t176;
t362 = t418 * pkin(2);
t360 = t244 * t381;
t358 = t395 / 0.2e1;
t355 = -t381 / 0.2e1;
t354 = t558 * t283 * t286;
t349 = -t179 * mrSges(4,1) + t180 * mrSges(4,2);
t346 = -t202 * pkin(5) + qJ(6) * t203;
t345 = -t204 * pkin(5) + qJ(6) * t205;
t343 = t298 * t280 - t294 * t290;
t276 = -t362 - pkin(3);
t340 = pkin(3) * t284 + t449;
t337 = mrSges(3,1) * t293 + mrSges(3,2) * t297;
t334 = mrSges(5,1) * t292 + t296 * mrSges(5,2);
t331 = Ifges(5,1) * t296 - t430;
t330 = t297 * Ifges(3,2) + t434;
t329 = -Ifges(5,2) * t292 + t429;
t328 = Ifges(3,5) * t297 - Ifges(3,6) * t293;
t327 = Ifges(5,5) * t296 - Ifges(5,6) * t292;
t325 = pkin(5) * t286 + qJ(6) * t285;
t52 = -t291 * t94 + t295 * t87;
t143 = t225 * t417 - t418 * t226;
t186 = -t418 * t260 - t262 * t417;
t311 = t318 * mrSges(5,3);
t132 = -t218 * mrSges(5,2) + t311;
t133 = mrSges(5,1) * t218 - t437;
t320 = t132 * t296 - t133 * t292;
t317 = t221 * pkin(4);
t142 = pkin(4) * t409 + t186;
t316 = pkin(1) * t337;
t219 = t284 * t399 + t394;
t9 = t291 * t40 + t295 * t49 + t87 * t378 - t379 * t94;
t314 = t360 + t412;
t313 = t156 * t334;
t312 = t293 * (Ifges(3,1) * t297 - t434);
t257 = t276 - t445;
t100 = pkin(4) * t315 + t143;
t307 = t219 * pkin(4);
t303 = t501 + t517;
t282 = Ifges(3,4) * t385;
t278 = -pkin(5) - t451;
t275 = qJ(6) + t453;
t266 = qJD(6) + t368;
t259 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t385;
t235 = Ifges(3,1) * t386 + Ifges(3,5) * qJD(2) + t282;
t234 = Ifges(3,6) * qJD(2) + qJD(1) * t330;
t222 = t284 * t394 + t399;
t220 = -t284 * t396 + t398;
t216 = Ifges(4,4) * t227;
t200 = -qJD(2) * mrSges(4,2) + t439;
t163 = -mrSges(4,1) * t227 - mrSges(4,2) * t228;
t159 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t180;
t158 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t179;
t153 = t246 * pkin(5) - t247 * qJ(6) + t257;
t152 = t246 * t244;
t150 = -t228 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t216;
t149 = t227 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t432;
t86 = -mrSges(5,2) * t176 + mrSges(5,3) * t120;
t85 = mrSges(5,1) * t176 - mrSges(5,3) * t119;
t81 = mrSges(6,1) * t125 + mrSges(6,2) * t304;
t80 = mrSges(7,1) * t125 - mrSges(7,3) * t304;
t73 = t151 * pkin(5) + t152 * qJ(6) + t142;
t72 = -t230 * t401 - t291 * t360 - t379 * t409 + (t408 * t512 - t413) * t295;
t71 = -t151 * t512 + t246 * t230;
t69 = -mrSges(5,1) * t120 + mrSges(5,2) * t119;
t61 = t454 + t79;
t58 = -qJD(4) * t116 + t344;
t55 = t119 * Ifges(5,4) + t120 * Ifges(5,2) + t176 * Ifges(5,6);
t42 = -pkin(5) * t243 - t52;
t41 = qJ(6) * t243 + t441;
t32 = t295 * t77 - t421;
t31 = t291 * t77 + t420;
t18 = t72 * pkin(5) - t71 * qJ(6) + t152 * qJD(6) + t100;
t17 = mrSges(6,1) * t51 + mrSges(6,2) * t50;
t16 = mrSges(7,1) * t51 - mrSges(7,3) * t50;
t8 = -pkin(5) * t229 - t10;
t7 = qJ(6) * t229 + qJD(6) * t243 + t9;
t1 = [(t215 * mrSges(4,2) - t109 * mrSges(4,3) + Ifges(4,1) * t180 + Ifges(4,4) * t179 + Ifges(4,5) * t509 + t101 * t334 + t113 * t355 + t327 * t473 + t329 * t482 + t331 * t483) * t244 - t562 * t152 + m(5) * (t115 * t26 + t116 * t25 + t57 * t90 + t58 * t89) + m(4) * (t110 * t187 + t144 * t171 - t215 * t280 + t255 * t371) + m(6) * (t10 * t27 + t100 * t121 + t142 * t59 + t28 * t9 + t441 * t5 + t52 * t6) + t441 * t35 + (t297 * (-Ifges(3,2) * t293 + t433) + t312) * t377 / 0.2e1 + (t215 * mrSges(4,1) - t110 * mrSges(4,3) - Ifges(4,4) * t180 + Ifges(5,5) * t483 - Ifges(4,2) * t179 - t509 * Ifges(4,6) + Ifges(5,6) * t482 + Ifges(6,6) * t491 + Ifges(7,6) * t490 + Ifges(5,3) * t473 + t537 * t474 + t539 * t492 + t501 + t505) * t243 + (t328 * t444 - t375) * qJD(2) - (-t400 / 0.2e1 + t255 * mrSges(4,2) - mrSges(4,3) * t170 + t150 / 0.2e1 + t358 + Ifges(4,5) * t444 + Ifges(4,1) * t464 + Ifges(4,4) * t465) * t230 + m(7) * (t11 * t73 + t18 * t54 + t2 * t41 + t23 * t8 + t24 * t7 + t3 * t42) - qJDD(2) * mrSges(3,2) * t450 + (-pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t251) + Ifges(3,1) * t251 + Ifges(3,4) * t547 + t509 * Ifges(3,5)) * t293 + t496 * t151 + t318 * (-Ifges(5,4) * t314 - Ifges(5,2) * t315) / 0.2e1 + t564 * t71 + (-m(4) * t109 + m(5) * t101 - t159 + t69) * t186 + t553 * t72 + (t511 / 0.2e1 + t543 - t542 + t531 / 0.2e1 + t255 * mrSges(4,1) - mrSges(4,3) * t171 + Ifges(7,6) * t480 + Ifges(6,6) * t481 - t28 * mrSges(6,2) + t27 * mrSges(6,1) - t23 * mrSges(7,1) + t24 * mrSges(7,3) - t149 / 0.2e1 + t530 / 0.2e1 - Ifges(4,6) * t444 - Ifges(4,4) * t464 - Ifges(4,2) * t465 + Ifges(5,5) * t470 + t539 * t477 + t537 * t468) * t229 + Ifges(3,6) * t297 * t419 + t156 * (mrSges(5,1) * t315 - mrSges(5,2) * t314) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t514) + (t363 + t517) * t243 / 0.2e1 + (-m(4) * t170 + m(5) * t156 - t520) * t143 + t251 * t433 / 0.2e1 + t330 * t547 + (-t25 * t409 - t26 * t408 + t314 * t89 - t315 * t90) * mrSges(5,3) + t297 * t235 * t444 + (-t222 * mrSges(5,1) - t221 * mrSges(5,2) - t569 * t343 + t548 * (pkin(4) * t399 + t343) - t508 * t205 - t506 * t204 + t510 * t294 + (-m(5) * t340 + t548 * t556 - t570) * t298) * g(2) - t55 * t409 / 0.2e1 + t163 * t371 + (-Ifges(5,1) * t314 - Ifges(5,4) * t315) * t470 + t218 * (-Ifges(5,5) * t314 - Ifges(5,6) * t315) / 0.2e1 - t234 * t384 / 0.2e1 - t316 * t377 + Ifges(2,3) * qJDD(1) - t259 * t367 - t112 * t359 / 0.2e1 + t297 * (Ifges(3,4) * t251 + Ifges(3,2) * t250 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t280 * t349 - pkin(1) * (-mrSges(3,1) * t250 + mrSges(3,2) * t251) + t408 * t489 + t144 * t200 + t187 * t158 + t142 * t17 + t57 * t132 + t58 * t133 + t115 * t85 + t116 * t86 + t100 * t81 + t7 * t103 + t9 * t104 + t10 * t105 + t8 * t106 + t18 * t80 + t73 * t16 + t52 * t33 + t42 * t34 + t41 * t36 - t261 * t416 + (m(5) * t402 - t220 * mrSges(5,1) - t219 * mrSges(5,2) + t548 * (pkin(4) * t398 + t294 * t405 - t402) + t508 * t203 + t506 * t202 + (m(4) * t290 + t510) * t298 + (-m(5) * (-t280 - t340) + m(4) * t280 + t548 * (-t280 - t253) + t570) * t294) * g(1) + (t250 * t450 + t514) * mrSges(3,3); (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (t101 * t276 - t156 * t177 - t89 * t98 - t90 * t99) * m(5) + (t121 * t518 + t167 * t5 + t257 * t59 + t27 * t529 + t28 * t525 + t319 * t6) * m(6) + (t11 * t153 + t167 * t2 + t23 * t528 + t24 * t527 - t3 * t319 + t532 * t54) * m(7) - t535 * t319 + (Ifges(4,1) * t227 + t432 + t511) * t228 / 0.2e1 + t562 * t247 + (-Ifges(6,4) * t141 - Ifges(6,6) * t228) * t480 + (t121 * t141 - t228 * t28) * mrSges(6,2) + (-Ifges(7,5) * t141 - Ifges(7,6) * t228) * t481 + (-t141 * t54 + t228 * t24) * mrSges(7,3) - t27 * (-mrSges(6,1) * t228 + mrSges(6,3) * t141) - t23 * (mrSges(7,1) * t228 - mrSges(7,2) * t141) - (-Ifges(3,2) * t386 + t235 + t282) * t385 / 0.2e1 + (t189 * t331 + t218 * t327) * qJD(4) / 0.2e1 - (Ifges(4,2) * t228 + t150 + t216 + t395) * t227 / 0.2e1 + (-t141 * t539 - t228 * t537) * t469 + (-t141 * t541 - t228 * t539) * t478 + ((t109 * t418 + t110 * t417) * pkin(2) + t170 * t177 - t171 * t178 - t255 * t372) * m(4) + t525 * t104 + t527 * t103 + t528 * t106 + t529 * t105 + (t313 + t358) * qJD(4) + t496 * t246 - t564 * t182 + (t234 / 0.2e1 + pkin(7) * t259) * t386 + t552 * t140 + t553 * t183 - t227 * t313 + (t375 + (-t312 / 0.2e1 + t316) * qJD(1)) * qJD(1) + t518 * t81 + t520 * t177 + t228 * t543 + (-Ifges(5,6) * t228 + t227 * t329) * t546 - t228 * t542 + t532 * t80 + t534 * t167 - t171 * t438 - t328 * t377 / 0.2e1 - t163 * t372 + t296 * t55 / 0.2e1 + t276 * t69 + Ifges(3,5) * t251 - t255 * (-mrSges(4,1) * t228 + mrSges(4,2) * t227) + t257 * t17 + Ifges(3,6) * t250 - t241 * mrSges(3,2) - t242 * mrSges(3,1) + t292 * t489 + (Ifges(5,2) * t296 + t430) * t482 + (Ifges(5,1) * t292 + t429) * t483 - qJD(2) * (Ifges(4,5) * t227 + Ifges(4,6) * t228) / 0.2e1 - t178 * t200 + Ifges(4,6) * t179 + Ifges(4,5) * t180 + t153 * t16 - t99 * t132 - t98 * t133 + t109 * mrSges(4,1) - t110 * mrSges(4,2) + (g(1) * t298 + g(2) * t294) * (t337 + t569 * t455 + t548 * (-t284 * t299 - t455) + (-m(5) * pkin(8) + t575) * t284 + (m(6) * t279 - m(7) * (-t279 - t325) + mrSges(4,1) + t561) * t283) + (-m(5) * (t287 + t449) - m(4) * t287 + t548 * (t287 + t556) + (-m(7) * t325 - t561) * t284 + t576) * g(3) + t329 * t308 / 0.2e1 + t101 * t335 + t112 * t355 + t158 * t361 + t159 * t362 + t141 * t560 + t170 * t439 + t149 * t464 + t400 * t465 + (-Ifges(5,3) * t228 + t227 * t327) * t467 + (-Ifges(5,5) * t228 + t227 * t331) * t471 + (Ifges(5,5) * t292 + Ifges(5,6) * t296) * t473 + (-t572 * t90 + (-t380 + t414) * t89 + t513) * mrSges(5,3) + (-t133 * t380 - t132 * t381 - t292 * t85 + m(5) * ((-t292 * t90 - t296 * t89) * qJD(4) + t513) + t296 * t86) * t274; t349 + t534 * t247 + t535 * t246 + (t80 + t81 - t520) * t228 + t296 * t85 + t292 * t86 + (-t200 - t320) * t227 + t320 * qJD(4) + t524 * t522 + t523 * t521 + (-g(1) * t294 + g(2) * t298) * (t569 - t548) + (t2 * t247 + t228 * t54 + t23 * t521 + t24 * t522 + t246 * t3) * m(7) + (t121 * t228 - t246 * t6 + t247 * t5 - t27 * t521 + t28 * t522) * m(6) + (t156 * t228 + t25 * t292 + t26 * t296 + t218 * (-t292 * t89 + t296 * t90)) * m(5) + (-t170 * t228 - t171 * t227 + t215) * m(4); (mrSges(6,1) * t285 + t334 + t440) * t446 + ((t291 * t5 + t295 * t6 + (-t27 * t291 + t28 * t295) * qJD(5)) * pkin(4) + t452 * t446 - t121 * t454 + t27 * t31 - t28 * t32) * m(6) + (t133 + t437) * t90 + (t311 - t132) * t89 + t505 + t303 + t523 * (t369 - t31) - m(7) * (t23 * t31 + t24 * t32 + t54 * t61) + t552 * t304 + (-Ifges(5,2) * t189 + t113 + t188) * t546 + (t368 - t32) * t104 + (t266 - t32) * t103 + (-m(7) * (-t307 + t346) + m(6) * t307 + mrSges(5,1) * t219 - mrSges(5,2) * t220 + t515) * g(2) + (-m(7) * (t317 + t345) - m(6) * t317 - mrSges(5,1) * t221 + mrSges(5,2) * t222 + t516) * g(1) - g(3) * ((m(7) * (-pkin(5) * t285 - t452) - t285 * mrSges(7,1)) * t283 - t354) - t81 * t454 + t363 + m(7) * (t2 * t275 + t23 * t369 + t24 * t266 + t278 * t3) + t278 * t34 + t275 * t36 - t61 * t80 - t156 * (t189 * mrSges(5,1) + mrSges(5,2) * t318) + t33 * t451 + t35 * t453 + (Ifges(5,5) * t318 - Ifges(5,6) * t189) * t467 + t112 * t470 + (Ifges(5,1) * t318 - t431) * t471 + (-Ifges(6,4) * t480 - Ifges(7,5) * t481 - t539 * t469 - t541 * t478 + t571) * t125; t303 + ((t285 * t508 + t440) * t283 + t354) * g(3) + (-t436 - t524) * t27 + (t125 * t23 + t24 * t304) * mrSges(7,2) + t516 * g(1) + t515 * g(2) + (Ifges(7,3) * t304 - t427) * t481 - t121 * (mrSges(6,1) * t304 - mrSges(6,2) * t125) - t54 * (mrSges(7,1) * t304 + mrSges(7,3) * t125) + qJD(6) * t103 - t79 * t80 + qJ(6) * t36 - pkin(5) * t34 + (t435 - t523) * t28 + t66 * t477 + (-t125 * t539 + t304 * t538) * t469 + (-Ifges(6,2) * t304 - t123 + t533) * t480 + (-t125 * t541 + t122 - t428 + t63) * t478 + (-pkin(5) * t3 - t345 * g(1) - t346 * g(2) + qJ(6) * t2 - t23 * t28 + t24 * t557 - t54 * t79) * m(7); -t212 * t103 + t304 * t80 + (-g(1) * t204 - g(2) * t202 - t212 * t24 - t285 * t446 + t304 * t54 + t3) * m(7) + t34;];
tau  = t1;
