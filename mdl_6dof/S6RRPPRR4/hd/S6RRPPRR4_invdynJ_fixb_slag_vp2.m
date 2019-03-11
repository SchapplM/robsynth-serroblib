% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:14
% EndTime: 2019-03-09 09:02:05
% DurationCPUTime: 30.44s
% Computational Cost: add. (14451->881), mult. (39977->1137), div. (0->0), fcn. (32146->12), ass. (0->424)
t499 = m(6) + m(7);
t404 = m(5) + t499;
t440 = cos(pkin(11));
t388 = t440 * pkin(2);
t286 = -t388 - pkin(3);
t283 = -pkin(9) + t286;
t298 = cos(qJ(5));
t408 = qJD(5) * t298;
t290 = sin(pkin(11));
t295 = sin(qJ(2));
t299 = cos(qJ(2));
t252 = -t299 * t290 - t295 * t440;
t291 = sin(pkin(6));
t231 = t252 * t291;
t225 = qJD(1) * t231;
t374 = t299 * t440;
t360 = t291 * t374;
t413 = qJD(1) * t291;
t387 = t295 * t413;
t222 = -qJD(1) * t360 + t290 * t387;
t271 = pkin(2) * t387;
t372 = qJ(4) * t222 + t271;
t498 = pkin(3) + pkin(9);
t100 = -t225 * t498 + t372;
t294 = sin(qJ(5));
t292 = cos(pkin(6));
t477 = pkin(1) * t292;
t282 = t299 * t477;
t273 = qJD(1) * t282;
t462 = pkin(8) + qJ(3);
t379 = t462 * t295;
t365 = t291 * t379;
t209 = -qJD(1) * t365 + t273;
t429 = t292 * t295;
t281 = pkin(1) * t429;
t431 = t291 * t299;
t210 = (t431 * t462 + t281) * qJD(1);
t373 = t440 * t210;
t127 = t209 * t290 + t373;
t475 = pkin(4) * t222;
t95 = t127 - t475;
t54 = t298 * t100 + t294 * t95;
t566 = t283 * t408 - t54;
t200 = t290 * t210;
t128 = t209 * t440 - t200;
t565 = t128 - qJD(4);
t465 = -mrSges(5,2) + mrSges(4,1);
t359 = pkin(5) * t298 + pkin(10) * t294;
t564 = qJD(5) * t359 + (-pkin(4) - t359) * t225 - t565;
t563 = pkin(10) * t222 + t566;
t403 = qJD(1) * qJD(2);
t238 = (qJDD(1) * t299 - t295 * t403) * t291;
t239 = (qJDD(1) * t295 + t299 * t403) * t291;
t164 = t290 * t238 + t239 * t440;
t162 = qJDD(5) + t164;
t493 = t162 / 0.2e1;
t163 = -t238 * t440 + t239 * t290;
t276 = qJD(1) * t292 + qJD(2);
t189 = t222 * t294 + t276 * t298;
t401 = qJDD(1) * t292;
t275 = qJDD(2) + t401;
t89 = -qJD(5) * t189 + t163 * t298 - t275 * t294;
t500 = t89 / 0.2e1;
t87 = qJDD(6) - t89;
t502 = t87 / 0.2e1;
t216 = qJD(5) - t225;
t293 = sin(qJ(6));
t297 = cos(qJ(6));
t118 = t189 * t297 + t216 * t293;
t188 = t222 * t298 - t276 * t294;
t88 = qJD(5) * t188 + t163 * t294 + t275 * t298;
t37 = -qJD(6) * t118 + t162 * t297 - t293 * t88;
t508 = t37 / 0.2e1;
t117 = -t189 * t293 + t216 * t297;
t36 = qJD(6) * t117 + t162 * t293 + t297 * t88;
t509 = t36 / 0.2e1;
t190 = pkin(2) * t276 + t209;
t114 = t190 * t440 - t200;
t327 = qJD(4) - t114;
t471 = t225 * pkin(4);
t70 = -t276 * t498 + t327 - t471;
t287 = pkin(2) * t299 + pkin(1);
t241 = -t287 * t413 + qJD(3);
t308 = qJ(4) * t225 + t241;
t93 = t222 * t498 + t308;
t41 = t294 * t70 + t298 * t93;
t31 = pkin(10) * t216 + t41;
t115 = t290 * t190 + t373;
t103 = -t276 * qJ(4) - t115;
t76 = -t103 - t475;
t49 = -pkin(5) * t188 - pkin(10) * t189 + t76;
t13 = -t293 * t31 + t297 * t49;
t414 = pkin(8) * t431 + t281;
t237 = t414 * qJD(2);
t395 = pkin(1) * t401;
t269 = t299 * t395;
t433 = t291 * t295;
t384 = qJD(3) * t433;
t402 = qJDD(1) * t291;
t394 = pkin(8) * t402;
t110 = -t295 * t394 + pkin(2) * t275 - qJ(3) * t239 + t269 + (-t237 - t384) * qJD(1);
t399 = qJD(2) * t477;
t366 = qJD(1) * t399;
t389 = t295 * t395 + (t366 + t394) * t299;
t410 = qJD(3) * t299;
t411 = qJD(2) * t295;
t119 = qJ(3) * t238 + (-pkin(8) * t411 + t410) * t413 + t389;
t60 = t290 * t110 + t440 * t119;
t57 = -t275 * qJ(4) - t276 * qJD(4) - t60;
t39 = -pkin(4) * t163 - t57;
t15 = -pkin(5) * t89 - pkin(10) * t88 + t39;
t59 = t110 * t440 - t290 * t119;
t325 = qJDD(4) - t59;
t38 = t164 * pkin(4) - t275 * t498 + t325;
t409 = qJD(5) * t294;
t205 = -pkin(1) * t402 - pkin(2) * t238 + qJDD(3);
t305 = -qJ(4) * t164 + qJD(4) * t225 + t205;
t48 = t163 * t498 + t305;
t7 = t294 * t38 + t298 * t48 + t70 * t408 - t409 * t93;
t5 = pkin(10) * t162 + t7;
t1 = qJD(6) * t13 + t15 * t293 + t297 * t5;
t14 = t293 * t49 + t297 * t31;
t2 = -qJD(6) * t14 + t15 * t297 - t293 * t5;
t517 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t9 = Ifges(7,5) * t36 + Ifges(7,6) * t37 + Ifges(7,3) * t87;
t561 = t517 + Ifges(7,5) * t509 + Ifges(7,6) * t508 + Ifges(7,3) * t502 - t88 * Ifges(6,4) / 0.2e1 + t9 / 0.2e1 + (-t493 - t162 / 0.2e1) * Ifges(6,6) + (-t500 - t89 / 0.2e1) * Ifges(6,2);
t296 = sin(qJ(1));
t421 = t296 * t299;
t300 = cos(qJ(1));
t424 = t295 * t300;
t244 = -t292 * t421 - t424;
t315 = t244 * pkin(2);
t560 = t76 * mrSges(6,2);
t542 = Ifges(4,5) - Ifges(5,4);
t541 = Ifges(5,5) - Ifges(4,6);
t181 = Ifges(6,4) * t188;
t559 = t188 * Ifges(6,2);
t558 = t188 * Ifges(6,6);
t557 = t216 * Ifges(6,5);
t556 = t216 * Ifges(6,6);
t555 = t216 * Ifges(6,3);
t554 = -mrSges(6,3) - t465;
t215 = Ifges(4,4) * t222;
t553 = -t225 * Ifges(4,1) + t276 * Ifges(4,5) + t189 * Ifges(6,5) - t215 + t555 + t558;
t351 = -mrSges(7,1) * t297 + mrSges(7,2) * t293;
t319 = m(7) * pkin(5) - t351;
t353 = mrSges(6,1) * t294 + mrSges(6,2) * t298;
t375 = pkin(10) * t298 - qJ(4);
t464 = mrSges(5,3) - mrSges(4,2);
t552 = -m(6) * qJ(4) + m(7) * t375 + t298 * mrSges(7,3) - t294 * t319 - t353 - t464;
t419 = t299 * t300;
t422 = t296 * t295;
t551 = t292 * t419 - t422;
t425 = t295 * t290;
t322 = t374 - t425;
t534 = t292 * t252;
t177 = t296 * t534 + t322 * t300;
t550 = -t296 * t322 + t300 * t534;
t122 = pkin(3) * t222 + t308;
t549 = t241 * mrSges(4,1) - t122 * mrSges(5,2);
t501 = t88 / 0.2e1;
t548 = Ifges(6,1) * t501 + Ifges(6,5) * t493;
t547 = mrSges(6,1) * t76 + mrSges(7,1) * t13 - mrSges(7,2) * t14;
t40 = -t294 * t93 + t298 * t70;
t546 = t40 * mrSges(6,1) + t241 * mrSges(4,2) - t41 * mrSges(6,2) - t122 * mrSges(5,3);
t12 = -mrSges(7,1) * t37 + mrSges(7,2) * t36;
t61 = mrSges(6,1) * t162 - mrSges(6,3) * t88;
t540 = -t12 + t61;
t539 = mrSges(6,1) + t319;
t396 = m(7) * pkin(10) + mrSges(7,3);
t370 = mrSges(6,2) - t396;
t230 = t291 * t425 - t360;
t278 = pkin(2) * t431;
t415 = -t230 * pkin(3) + t278;
t358 = -qJ(4) * t231 + t415;
t478 = pkin(1) * t291;
t144 = -t358 - t478;
t474 = pkin(9) * t230;
t111 = t144 + t474;
t206 = pkin(2) * t292 + t282 - t365;
t217 = qJ(3) * t431 + t414;
t335 = -t206 * t440 + t290 * t217;
t94 = -t231 * pkin(4) - t292 * t498 + t335;
t538 = t298 * t111 + t294 * t94;
t476 = pkin(2) * t290;
t249 = pkin(5) * t294 - t375 + t476;
t428 = t293 * t294;
t197 = t249 * t297 - t283 * t428;
t537 = qJD(6) * t197 + t293 * t564 + t297 * t563;
t426 = t294 * t297;
t198 = t249 * t293 + t283 * t426;
t536 = -qJD(6) * t198 - t293 * t563 + t297 * t564;
t456 = mrSges(6,3) * t189;
t124 = mrSges(6,1) * t216 - t456;
t67 = -mrSges(7,1) * t117 + mrSges(7,2) * t118;
t535 = -t124 + t67;
t145 = -t222 * t297 + t225 * t428;
t381 = t293 * t409;
t405 = qJD(6) * t298;
t533 = t297 * t405 + t145 - t381;
t146 = -t222 * t293 - t225 * t426;
t532 = t293 * t405 + t297 * t409 + t146;
t458 = mrSges(4,3) * t225;
t460 = mrSges(5,1) * t225;
t416 = t276 * t465 + t458 + t460;
t106 = -mrSges(6,1) * t188 + mrSges(6,2) * t189;
t461 = mrSges(5,1) * t222;
t195 = -mrSges(5,3) * t276 + t461;
t531 = -t195 + t106;
t436 = t225 * t298;
t530 = t408 - t436;
t24 = mrSges(7,1) * t87 - mrSges(7,3) * t36;
t25 = -mrSges(7,2) * t87 + mrSges(7,3) * t37;
t529 = -t293 * t24 + t297 * t25;
t8 = -qJD(5) * t41 - t294 * t48 + t298 * t38;
t528 = -t294 * t7 - t298 * t8;
t527 = t1 * t297 - t2 * t293;
t455 = Ifges(3,4) * t295;
t526 = -t295 * (Ifges(3,1) * t299 - t455) / 0.2e1 + pkin(1) * (mrSges(3,1) * t295 + mrSges(3,2) * t299);
t385 = t291 * t411;
t369 = pkin(8) * t385;
t182 = -qJD(1) * t369 + t389;
t183 = -pkin(8) * t239 - t295 * t366 + t269;
t58 = -t275 * pkin(3) + t325;
t524 = t183 * mrSges(3,1) + t59 * mrSges(4,1) - t182 * mrSges(3,2) - t60 * mrSges(4,2) + t58 * mrSges(5,2) - t57 * mrSges(5,3) + Ifges(3,5) * t239 + Ifges(3,6) * t238;
t523 = -pkin(9) * t499 + t554;
t6 = -pkin(5) * t162 - t8;
t522 = qJD(5) * (t13 * t293 - t14 * t297) + t6;
t350 = t293 * mrSges(7,1) + t297 * mrSges(7,2);
t521 = m(7) * pkin(9) + t350 - t554;
t520 = qJ(4) * t404 + t464;
t519 = t8 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,5) * t88 + Ifges(6,6) * t89 + Ifges(6,3) * t162;
t457 = mrSges(6,3) * t188;
t123 = -mrSges(6,2) * t216 + t457;
t187 = qJD(6) - t188;
t71 = -mrSges(7,2) * t187 + mrSges(7,3) * t117;
t72 = mrSges(7,1) * t187 - mrSges(7,3) * t118;
t518 = qJD(5) * (t293 * t72 - t297 * t71 - t123) + t225 * t123 - t540;
t223 = qJD(2) * t231;
t224 = qJD(2) * t360 - t290 * t385;
t324 = pkin(2) * t385 - qJ(4) * t224 + qJD(4) * t231;
t80 = -t223 * t498 + t324;
t274 = t299 * t399;
t191 = t274 + (-qJD(2) * t379 + t410) * t291;
t380 = t462 * t291;
t192 = -t384 + (-t299 * t380 - t281) * qJD(2);
t107 = t191 * t290 - t440 * t192;
t81 = pkin(4) * t224 + t107;
t19 = -qJD(5) * t538 - t294 * t80 + t298 * t81;
t479 = t298 / 0.2e1;
t50 = t118 * Ifges(7,5) + t117 * Ifges(7,6) + t187 * Ifges(7,3);
t453 = Ifges(6,4) * t189;
t78 = t453 + t556 + t559;
t503 = -t78 / 0.2e1;
t516 = t76 * (mrSges(6,1) * t298 - mrSges(6,2) * t294) + t298 * t503 + t50 * t479;
t515 = -m(6) * pkin(9) - t521;
t514 = m(5) * qJ(4) - t552;
t10 = t36 * Ifges(7,4) + t37 * Ifges(7,2) + t87 * Ifges(7,6);
t512 = t10 / 0.2e1;
t511 = Ifges(7,1) * t509 + Ifges(7,4) * t508 + Ifges(7,5) * t502;
t450 = Ifges(7,4) * t118;
t51 = Ifges(7,2) * t117 + Ifges(7,6) * t187 + t450;
t507 = -t51 / 0.2e1;
t506 = t51 / 0.2e1;
t116 = Ifges(7,4) * t117;
t52 = Ifges(7,1) * t118 + Ifges(7,5) * t187 + t116;
t505 = -t52 / 0.2e1;
t504 = t52 / 0.2e1;
t497 = -t117 / 0.2e1;
t496 = t117 / 0.2e1;
t495 = -t118 / 0.2e1;
t494 = t118 / 0.2e1;
t492 = -t187 / 0.2e1;
t491 = t187 / 0.2e1;
t490 = -t188 / 0.2e1;
t489 = -t189 / 0.2e1;
t488 = t189 / 0.2e1;
t487 = -t216 / 0.2e1;
t486 = -t222 / 0.2e1;
t485 = t222 / 0.2e1;
t483 = t225 / 0.2e1;
t482 = -t225 / 0.2e1;
t481 = -t276 / 0.2e1;
t480 = t276 / 0.2e1;
t469 = t298 * t6;
t463 = -Ifges(4,4) - Ifges(5,6);
t459 = mrSges(4,3) * t222;
t454 = Ifges(4,4) * t225;
t452 = Ifges(6,4) * t294;
t451 = Ifges(6,4) * t298;
t449 = Ifges(7,4) * t293;
t448 = Ifges(7,4) * t297;
t447 = Ifges(5,6) * t225;
t446 = t276 * Ifges(3,5);
t445 = t276 * Ifges(3,6);
t439 = t188 * t293;
t438 = t188 * t297;
t437 = t225 * t294;
t432 = t291 * t296;
t430 = t291 * t300;
t427 = t293 * t298;
t420 = t297 * t298;
t108 = t440 * t191 + t290 * t192;
t132 = t290 * t206 + t440 * t217;
t407 = qJD(6) * t293;
t406 = qJD(6) * t297;
t398 = m(4) + t404;
t397 = -m(3) * pkin(1) - mrSges(2,1);
t99 = -t292 * qJD(4) - t108;
t125 = -t292 * qJ(4) - t132;
t386 = t299 * t413;
t383 = t283 * t409;
t377 = -t409 / 0.2e1;
t376 = -t405 / 0.2e1;
t143 = t164 * mrSges(5,1) + t275 * mrSges(5,2);
t240 = pkin(2) * t429 - t380;
t371 = -t240 * t296 + t300 * t287;
t368 = mrSges(3,3) * t387;
t367 = mrSges(3,3) * t386;
t362 = t551 * pkin(2);
t357 = t177 * pkin(3) + t371;
t204 = t230 * t294 + t292 * t298;
t332 = t230 * t298 - t292 * t294;
t355 = -mrSges(6,1) * t332 + mrSges(6,2) * t204;
t138 = -t204 * t293 - t231 * t297;
t139 = t204 * t297 - t231 * t293;
t352 = mrSges(7,1) * t138 - mrSges(7,2) * t139;
t349 = Ifges(6,1) * t294 + t451;
t348 = Ifges(7,1) * t297 - t449;
t347 = Ifges(7,1) * t293 + t448;
t346 = Ifges(6,2) * t298 + t452;
t345 = -Ifges(7,2) * t293 + t448;
t344 = Ifges(7,2) * t297 + t449;
t343 = Ifges(6,5) * t294 + Ifges(6,6) * t298;
t342 = Ifges(7,5) * t297 - Ifges(7,6) * t293;
t341 = Ifges(7,5) * t293 + Ifges(7,6) * t297;
t46 = -pkin(10) * t231 + t538;
t98 = -pkin(4) * t230 - t125;
t63 = -pkin(5) * t332 - pkin(10) * t204 + t98;
t21 = t293 * t63 + t297 * t46;
t20 = -t293 * t46 + t297 * t63;
t337 = t294 * t40 - t298 * t41;
t53 = -t100 * t294 + t298 * t95;
t55 = -t111 * t294 + t298 * t94;
t73 = pkin(4) * t223 - t99;
t330 = pkin(4) * t432 + t357;
t310 = t292 * t322;
t173 = t296 * t252 + t300 * t310;
t152 = t173 * t294 + t298 * t430;
t150 = -t173 * t298 + t294 * t430;
t30 = -pkin(5) * t216 - t40;
t326 = t30 * t350;
t18 = -t111 * t409 + t294 * t81 + t298 * t80 + t94 * t408;
t176 = t252 * t300 - t296 * t310;
t317 = g(1) * t176 + g(2) * t173 - g(3) * t230;
t307 = (-t13 * t297 - t14 * t293) * qJD(6) + t527;
t306 = -qJD(5) * t337 - t528;
t303 = qJD(5) * t30 + t307;
t62 = -mrSges(6,2) * t162 + mrSges(6,3) * t89;
t301 = t62 + (-t293 * t71 - t297 * t72) * qJD(6) + t216 * t535 + t529;
t284 = qJ(4) + t476;
t270 = Ifges(3,4) * t386;
t264 = Ifges(5,1) * t275;
t263 = Ifges(3,3) * t275;
t262 = Ifges(4,3) * t275;
t253 = -t278 - t478;
t247 = -pkin(8) * t433 + t282;
t246 = (-mrSges(3,1) * t299 + mrSges(3,2) * t295) * t291;
t245 = -t292 * t422 + t419;
t243 = -t292 * t424 - t421;
t236 = t274 - t369;
t235 = t414 * qJD(1);
t234 = -pkin(8) * t387 + t273;
t233 = -mrSges(3,2) * t276 + t367;
t232 = mrSges(3,1) * t276 - t368;
t214 = Ifges(5,6) * t222;
t208 = Ifges(3,1) * t387 + t270 + t446;
t207 = t445 + (t299 * Ifges(3,2) + t455) * t413;
t193 = -mrSges(4,2) * t276 - t459;
t179 = t225 * t420 + t276 * t293;
t178 = -t225 * t427 + t276 * t297;
t169 = t550 * pkin(3);
t161 = Ifges(5,4) * t164;
t160 = Ifges(4,5) * t164;
t159 = Ifges(5,5) * t163;
t158 = Ifges(4,6) * t163;
t157 = t164 * mrSges(5,3);
t156 = t164 * mrSges(4,2);
t154 = -mrSges(5,2) * t222 + mrSges(5,3) * t225;
t153 = mrSges(4,1) * t222 - mrSges(4,2) * t225;
t148 = -t176 * t294 + t298 * t432;
t147 = t176 * t298 + t294 * t432;
t142 = mrSges(5,1) * t163 - mrSges(5,3) * t275;
t141 = mrSges(4,1) * t275 - mrSges(4,3) * t164;
t140 = -mrSges(4,2) * t275 - mrSges(4,3) * t163;
t137 = -pkin(3) * t225 + t372;
t135 = -t222 * Ifges(4,2) + t276 * Ifges(4,6) - t454;
t134 = t276 * Ifges(5,4) + t225 * Ifges(5,2) + t214;
t133 = t276 * Ifges(5,5) + t222 * Ifges(5,3) + t447;
t130 = qJD(5) * t204 + t223 * t298;
t129 = qJD(5) * t332 - t223 * t294;
t126 = -t292 * pkin(3) + t335;
t109 = pkin(5) * t189 - pkin(10) * t188;
t105 = -pkin(3) * t223 + t324;
t101 = -t276 * pkin(3) + t327;
t96 = t128 + t471;
t84 = t148 * t297 + t177 * t293;
t83 = -t148 * t293 + t177 * t297;
t79 = t189 * Ifges(6,1) + t181 + t557;
t66 = pkin(3) * t163 + t305;
t65 = qJD(6) * t138 + t129 * t297 + t224 * t293;
t64 = -qJD(6) * t139 - t129 * t293 + t224 * t297;
t45 = pkin(5) * t231 - t55;
t43 = pkin(5) * t222 - t53;
t42 = -mrSges(6,1) * t89 + mrSges(6,2) * t88;
t32 = pkin(5) * t130 - pkin(10) * t129 + t73;
t29 = t88 * Ifges(6,1) + t89 * Ifges(6,4) + t162 * Ifges(6,5);
t27 = t109 * t293 + t297 * t40;
t26 = t109 * t297 - t293 * t40;
t17 = -pkin(5) * t224 - t19;
t16 = pkin(10) * t224 + t18;
t4 = -qJD(6) * t21 - t16 * t293 + t297 * t32;
t3 = qJD(6) * t20 + t16 * t297 + t293 * t32;
t11 = [(-m(5) * t169 - t243 * mrSges(3,1) + t551 * mrSges(3,2) + (t240 * t398 + mrSges(2,2)) * t300 + (t287 * t398 - t397) * t296 - t520 * t173 + t370 * t150 - t539 * t152 - (t350 - t523) * t550 + t499 * (-pkin(4) * t430 - t169)) * g(1) + ((mrSges(3,1) * t238 - mrSges(3,2) * t239 + (m(3) * t478 - t246) * qJDD(1)) * pkin(1) + (mrSges(3,3) * t182 + Ifges(3,4) * t239 + Ifges(3,2) * t238 + Ifges(3,6) * t275) * t299 + (-mrSges(3,3) * t183 + Ifges(3,1) * t239 + Ifges(3,4) * t238 + Ifges(3,5) * t275) * t295 + ((t208 / 0.2e1 - t234 * mrSges(3,3) + t446 / 0.2e1) * t299 + (-t207 / 0.2e1 - t235 * mrSges(3,3) - t445 / 0.2e1 + (m(4) * t241 + t153) * pkin(2)) * t295 + (t299 * (Ifges(3,4) * t299 - Ifges(3,2) * t295) / 0.2e1 - t526) * t413) * qJD(2) + (g(1) * t300 + g(2) * t296) * (-m(3) * pkin(8) - mrSges(5,1) - mrSges(3,3) - mrSges(4,3))) * t291 + (t1 * t138 - t13 * t65 - t139 * t2 + t14 * t64) * mrSges(7,3) + (Ifges(7,5) * t139 + Ifges(7,6) * t138) * t502 + (Ifges(7,5) * t65 + Ifges(7,6) * t64) * t491 + (Ifges(7,4) * t65 + Ifges(7,2) * t64) * t496 + (Ifges(7,4) * t139 + Ifges(7,2) * t138) * t508 + (-m(5) * t357 - t245 * mrSges(3,1) - t244 * mrSges(3,2) - m(6) * t330 - t148 * mrSges(6,1) - m(7) * (pkin(5) * t148 + t330) - t84 * mrSges(7,1) - t83 * mrSges(7,2) + t296 * mrSges(2,2) - m(4) * t371 + t397 * t300 + t520 * t176 + t370 * t147 + t523 * t177) * g(2) + (Ifges(7,1) * t65 + Ifges(7,4) * t64) * t494 + (Ifges(7,1) * t139 + Ifges(7,4) * t138) * t509 + t414 * (-mrSges(3,2) * t275 + mrSges(3,3) * t238) - (-t115 * mrSges(4,3) - Ifges(4,4) * t482 + Ifges(5,6) * t483 + Ifges(5,3) * t485 - Ifges(4,2) * t486 + t103 * mrSges(5,1) - t135 / 0.2e1 + t133 / 0.2e1 + t541 * t480 + t549) * t223 + t139 * t511 + t138 * t512 - t416 * t107 + ((Ifges(3,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t275 + (-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t164 + t262 / 0.2e1 + t263 / 0.2e1 + t264 / 0.2e1 + (Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1) * t163 + t160 / 0.2e1 - t161 / 0.2e1 + t159 / 0.2e1 - t158 / 0.2e1 + t524) * t292 + m(7) * (t1 * t21 + t13 * t4 + t14 * t3 + t17 * t30 + t2 * t20 + t45 * t6) + m(5) * (t101 * t107 + t103 * t99 + t105 * t122 + t125 * t57 + t126 * t58 + t144 * t66) + t65 * t504 + t64 * t506 + m(6) * (t18 * t41 + t19 * t40 + t39 * t98 + t538 * t7 + t55 * t8 + t73 * t76) + t538 * t62 + m(4) * (-t107 * t114 + t108 * t115 + t132 * t60 + t205 * t253 - t335 * t59) - t335 * t141 + (mrSges(4,1) * t205 + mrSges(5,1) * t57 - mrSges(5,2) * t66 - mrSges(4,3) * t60 + t541 * t275 + t463 * t164 + (Ifges(4,2) + Ifges(5,3)) * t163) * t230 - (-t66 * mrSges(5,3) + t205 * mrSges(4,2) - t59 * mrSges(4,3) + t58 * mrSges(5,1) + t542 * t275 + (Ifges(4,1) + Ifges(5,2)) * t164 + t463 * t163 + t519) * t231 + (Ifges(6,1) * t488 + t557 / 0.2e1 + t181 / 0.2e1 + t79 / 0.2e1 + t560 - t40 * mrSges(6,3)) * t129 + t247 * (mrSges(3,1) * t275 - mrSges(3,3) * t239) + (Ifges(7,5) * t494 + Ifges(7,6) * t496 + t503 - Ifges(6,4) * t488 + Ifges(7,3) * t491 - t556 / 0.2e1 - t559 / 0.2e1 + t50 / 0.2e1 - t41 * mrSges(6,3) + t547) * t130 + (Ifges(6,4) * t500 + t29 / 0.2e1 - t8 * mrSges(6,3) + t548) * t204 + t253 * (t163 * mrSges(4,1) + t156) + t236 * t233 - t237 * t232 + Ifges(2,3) * qJDD(1) + t108 * t193 + t99 * t195 + t144 * (-t163 * mrSges(5,2) - t157) + t105 * t154 + t132 * t140 + t125 * t142 + t126 * t143 + t18 * t123 + t19 * t124 + t73 * t106 + t98 * t42 + t3 * t71 + t4 * t72 + t30 * (-mrSges(7,1) * t64 + mrSges(7,2) * t65) + t17 * t67 + t55 * t61 + t45 * t12 + t21 * t25 + t20 * t24 + (-t114 * mrSges(4,3) + Ifges(4,1) * t482 - Ifges(5,2) * t483 - Ifges(5,6) * t485 + Ifges(4,4) * t486 + Ifges(6,5) * t488 + t101 * mrSges(5,1) - t134 / 0.2e1 + t555 / 0.2e1 + t558 / 0.2e1 + t542 * t480 + t546 + t553 / 0.2e1) * t224 + (mrSges(6,3) * t7 + Ifges(6,4) * t501 - t561) * t332 - t6 * t352 + t39 * t355 + m(3) * (t182 * t414 + t183 * t247 - t234 * t237 + t235 * t236); (-t101 * t127 + t103 * t565 - t122 * t137 - t284 * t57 + t286 * t58) * m(5) + t526 * qJD(1) ^ 2 * t291 ^ 2 + (-t193 + t195) * t128 + t540 * t283 * t298 + ((Ifges(7,5) * t298 - t294 * t348) * t494 + (Ifges(7,6) * t298 - t294 * t345) * t496 + (Ifges(7,3) * t298 - t294 * t342) * t491 + t516) * qJD(5) + (t42 - t142) * t284 + t451 * t500 + (t367 - t233) * t234 - t452 * t501 + (t383 - t43) * t67 + (-m(4) * t362 - mrSges(3,1) * t551 - mrSges(3,2) * t243 - t404 * (t173 * pkin(3) + t362) + t515 * t173 + t514 * t550) * g(2) + (t376 * t51 + t377 * t52) * t297 + (t368 + t232) * t235 + (-t383 - t53) * t124 + t207 * t387 / 0.2e1 - (t188 * t346 + t189 * t349 + t216 * t343) * qJD(5) / 0.2e1 + (t214 + t134) * t486 - t10 * t427 / 0.2e1 + t420 * t511 + t103 * t460 + t524 + (Ifges(7,1) * t146 + Ifges(7,4) * t145 + Ifges(7,5) * t436) * t495 + (Ifges(7,4) * t146 + Ifges(7,2) * t145 + Ifges(7,6) * t436) * t497 + (Ifges(7,5) * t146 + Ifges(7,6) * t145 + Ifges(7,3) * t436) * t492 + (-t530 * t41 + (t409 - t437) * t40 + t528) * mrSges(6,3) + (t437 / 0.2e1 + t377) * t79 + (t135 - t447) * t482 + (t133 + t454) * t483 - t153 * t271 + t146 * t505 + t381 * t506 + t145 * t507 + t350 * t469 + t140 * t476 + t29 * t479 + t101 * t461 - t115 * t458 + (-t215 + t553) * t485 - ((-Ifges(3,2) * t387 + t208 + t270) * t299 + t276 * (Ifges(3,5) * t299 - Ifges(3,6) * t295)) * t413 / 0.2e1 + t566 * t123 + t531 * qJD(4) + t416 * t127 + (-t14 * t530 - t30 * t532) * mrSges(7,2) + (t13 * t530 + t30 * t533) * mrSges(7,1) + t262 + t263 + t264 + (-t341 * t491 - t344 * t496 - t347 * t494) * t405 + (t283 * t306 + t284 * t39 - t40 * t53 - t41 * t54 + (-t96 + qJD(4)) * t76) * m(6) + (t114 * t127 - t115 * t128 - t241 * t271 + (t290 * t60 + t440 * t59) * pkin(2)) * m(4) + (-t1 * t427 + t13 * t532 - t14 * t533 - t2 * t420) * mrSges(7,3) + t293 * t52 * t376 - t114 * t459 + (-m(4) * t315 - mrSges(3,1) * t244 + mrSges(3,2) * t245 + t515 * t176 - t514 * t177 - t404 * (t176 * pkin(3) + t315)) * g(1) + (t342 * t502 + t345 * t508 + t348 * t509 + t548) * t298 + t160 - t161 + t159 - t158 - (-Ifges(4,2) * t485 + Ifges(5,3) * t486 + t343 * t487 + t346 * t490 + t349 * t489 + t481 * t541 + t516 - t549) * t225 + t197 * t24 + t198 * t25 - t137 * t154 - t96 * t106 + t141 * t388 + t286 * t143 + (-m(7) * t415 - m(6) * (t415 - t474) - m(4) * t278 - m(5) * t358 + t246 - t552 * t231 + t521 * t230) * g(3) + (t283 * t62 + t561) * t294 + t536 * t72 + t537 * t71 + (t1 * t198 + t197 * t2 + (t30 * t409 - t469) * t283 - t30 * t43 + t537 * t14 + t536 * t13) * m(7) + (-Ifges(4,1) * t483 - Ifges(6,5) * t489 + Ifges(5,2) * t482 - Ifges(6,6) * t490 - Ifges(6,3) * t487 - t481 * t542 + t546) * t222 + t39 * t353; -t145 * t72 - t146 * t71 + t156 - t157 - t416 * t225 + t465 * t163 + (t193 + t531) * t222 + t518 * t294 + t301 * t298 + (-t13 * t145 - t14 * t146 + t294 * t522 + t303 * t298 - t30 * t436) * m(7) + (t222 * t76 - t294 * t8 + t298 * t7 - t216 * (t294 * t41 + t298 * t40)) * m(6) + (t101 * t225 - t103 * t222 + t66) * m(5) + (-t114 * t225 + t115 * t222 + t205) * m(4) + (-t292 * g(3) + (-g(1) * t296 + g(2) * t300) * t291) * t398; -t225 * t154 - t178 * t72 - t179 * t71 - t531 * t276 - t518 * t298 + t301 * t294 + t143 + (-t13 * t178 - t14 * t179 + t303 * t294 - t298 * t522 - t30 * t437 + t317) * m(7) + (t225 * t337 - t276 * t76 + t306 + t317) * m(6) + (t103 * t276 - t122 * t225 + t317 + t58) * m(5); t519 + ((-t407 + t439) * t14 + (-t406 + t438) * t13 + t527) * mrSges(7,3) + (-pkin(5) * t6 - t13 * t26 - t14 * t27) * m(7) + (t117 * t345 + t118 * t348 + t187 * t342) * qJD(6) / 0.2e1 + (t50 - t453) * t489 + (Ifges(7,5) * t495 - Ifges(6,2) * t490 - Ifges(6,6) * t487 + Ifges(7,6) * t497 + Ifges(7,3) * t492 - t547) * t189 + t293 * t511 + t297 * t512 + t341 * t502 + t406 * t504 + t438 * t505 + t439 * t506 + t407 * t507 + t344 * t508 + t347 * t509 + t78 * t488 + (-t204 * t396 - t319 * t332 + t355) * g(3) + (m(7) * t307 - t406 * t72 - t407 * t71 + t529) * pkin(10) + (Ifges(6,1) * t489 + Ifges(6,5) * t487 + t342 * t492 + t345 * t497 + t348 * t495 - t326 - t560) * t188 + (t181 + t79) * t490 + (t457 - t123) * t40 + (t147 * t539 + t148 * t370) * g(1) + (-t150 * t539 - t152 * t370) * g(2) - t27 * t71 - t26 * t72 - pkin(5) * t12 + qJD(6) * t326 + (-m(7) * t30 + t456 - t535) * t41 + t6 * t351; -t30 * (mrSges(7,1) * t118 + mrSges(7,2) * t117) + (Ifges(7,1) * t117 - t450) * t495 + t51 * t494 + (Ifges(7,5) * t117 - Ifges(7,6) * t118) * t492 - t13 * t71 + t14 * t72 - g(1) * (mrSges(7,1) * t83 - mrSges(7,2) * t84) - g(2) * ((t152 * t293 - t297 * t550) * mrSges(7,1) + (t152 * t297 + t293 * t550) * mrSges(7,2)) - g(3) * t352 + (t117 * t13 + t118 * t14) * mrSges(7,3) + t9 + (-Ifges(7,2) * t118 + t116 + t52) * t497 + t517;];
tau  = t11;
