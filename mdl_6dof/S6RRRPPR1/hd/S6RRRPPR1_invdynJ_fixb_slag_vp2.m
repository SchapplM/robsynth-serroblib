% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:21:14
% EndTime: 2019-03-09 15:21:52
% DurationCPUTime: 21.59s
% Computational Cost: add. (21080->807), mult. (49288->1055), div. (0->0), fcn. (37347->18), ass. (0->382)
t385 = qJD(2) + qJD(3);
t532 = t385 / 0.2e1;
t393 = sin(qJ(3));
t394 = sin(qJ(2));
t397 = cos(qJ(3));
t398 = cos(qJ(2));
t323 = -t393 * t394 + t397 * t398;
t298 = t323 * qJD(1);
t324 = t393 * t398 + t394 * t397;
t299 = t324 * qJD(1);
t388 = sin(pkin(10));
t390 = cos(pkin(10));
t243 = -t390 * t298 + t299 * t388;
t540 = -t243 / 0.2e1;
t421 = t298 * t388 + t390 * t299;
t592 = t421 * Ifges(5,4);
t613 = -t592 / 0.2e1 - Ifges(5,2) * t540 - Ifges(5,6) * t532;
t612 = -mrSges(6,3) - mrSges(7,3);
t386 = qJ(2) + qJ(3);
t377 = pkin(10) + t386;
t362 = sin(t377);
t384 = pkin(11) + qJ(6);
t375 = sin(t384);
t512 = mrSges(7,2) * t375;
t387 = sin(pkin(11));
t513 = mrSges(6,2) * t387;
t611 = (-t512 - t513) * t362;
t400 = -pkin(8) - pkin(7);
t350 = t400 * t398;
t330 = qJD(1) * t350;
t304 = t393 * t330;
t349 = t400 * t394;
t329 = qJD(1) * t349;
t311 = qJD(2) * pkin(2) + t329;
t251 = t397 * t311 + t304;
t290 = t299 * qJ(4);
t218 = t251 - t290;
t206 = pkin(3) * t385 + t218;
t307 = t397 * t330;
t252 = t311 * t393 - t307;
t497 = qJ(4) * t298;
t219 = t252 + t497;
t474 = t390 * t219;
t138 = t388 * t206 + t474;
t392 = sin(qJ(6));
t396 = cos(qJ(6));
t389 = cos(pkin(11));
t223 = t385 * t387 + t389 * t421;
t130 = qJ(5) * t385 + t138;
t381 = t398 * pkin(2);
t368 = t381 + pkin(1);
t348 = t368 * qJD(1);
t264 = -pkin(3) * t298 + qJD(4) - t348;
t148 = pkin(4) * t243 - qJ(5) * t421 + t264;
t73 = -t130 * t387 + t389 * t148;
t55 = pkin(5) * t243 - pkin(9) * t223 + t73;
t439 = t389 * t385 - t387 * t421;
t74 = t389 * t130 + t387 * t148;
t63 = pkin(9) * t439 + t74;
t17 = -t392 * t63 + t396 * t55;
t18 = t392 * t55 + t396 * t63;
t591 = t439 * Ifges(6,6);
t595 = t223 * Ifges(6,5);
t610 = t138 * mrSges(5,3) + t18 * mrSges(7,2) + t74 * mrSges(6,2) - t17 * mrSges(7,1) - t264 * mrSges(5,1) - t73 * mrSges(6,1) - t591 / 0.2e1 - t595 / 0.2e1 - t613;
t539 = t243 / 0.2e1;
t152 = t223 * t396 + t392 * t439;
t237 = qJD(6) + t243;
t602 = -t223 * t392 + t396 * t439;
t609 = t152 * Ifges(7,5) + Ifges(7,6) * t602 + t243 * Ifges(6,3) + t237 * Ifges(7,3) + t591 + t595;
t320 = t387 * t396 + t389 * t392;
t180 = t320 * t243;
t295 = t320 * qJD(6);
t608 = t180 + t295;
t418 = t387 * t392 - t389 * t396;
t181 = t418 * t243;
t294 = t418 * qJD(6);
t607 = t181 + t294;
t363 = cos(t377);
t378 = sin(t386);
t379 = cos(t386);
t606 = mrSges(4,1) * t378 + mrSges(5,1) * t362 + mrSges(4,2) * t379 + mrSges(5,2) * t363;
t461 = qJD(1) * qJD(2);
t333 = qJDD(1) * t398 - t394 * t461;
t376 = cos(t384);
t515 = mrSges(7,1) * t376;
t516 = mrSges(6,1) * t389;
t605 = (t515 + t516) * t362;
t593 = t421 * Ifges(5,1);
t604 = -t243 * Ifges(5,4) + t385 * Ifges(5,5) + t389 * (t223 * Ifges(6,1) + Ifges(6,4) * t439 + t243 * Ifges(6,5)) + t593;
t382 = qJDD(2) + qJDD(3);
t334 = qJDD(1) * t394 + t398 * t461;
t316 = t334 * pkin(7);
t263 = qJDD(2) * pkin(2) - pkin(8) * t334 - t316;
t315 = t333 * pkin(7);
t267 = pkin(8) * t333 + t315;
t173 = -qJD(3) * t252 + t397 * t263 - t267 * t393;
t412 = t323 * qJD(3);
t229 = qJD(1) * t412 + t333 * t393 + t334 * t397;
t103 = pkin(3) * t382 - qJ(4) * t229 - qJD(4) * t299 + t173;
t463 = qJD(3) * t397;
t464 = qJD(3) * t393;
t172 = t393 * t263 + t397 * t267 + t311 * t463 + t330 * t464;
t413 = t324 * qJD(3);
t230 = -qJD(1) * t413 + t333 * t397 - t334 * t393;
t107 = qJ(4) * t230 + qJD(4) * t298 + t172;
t48 = t388 * t103 + t390 * t107;
t37 = qJ(5) * t382 + qJD(5) * t385 + t48;
t166 = t229 * t388 - t390 * t230;
t167 = t229 * t390 + t230 * t388;
t496 = qJDD(1) * pkin(1);
t291 = -pkin(2) * t333 - t496;
t196 = -pkin(3) * t230 + qJDD(4) + t291;
t58 = pkin(4) * t166 - qJ(5) * t167 - qJD(5) * t421 + t196;
t14 = -t37 * t387 + t389 * t58;
t15 = t389 * t37 + t387 * t58;
t423 = -t14 * t387 + t15 * t389;
t430 = t513 - t516;
t395 = sin(qJ(1));
t399 = cos(qJ(1));
t570 = g(1) * t399 + g(2) * t395;
t603 = -t379 * mrSges(4,1) - t363 * mrSges(5,1) + mrSges(4,2) * t378 + (mrSges(5,2) - mrSges(7,3)) * t362;
t211 = t388 * t219;
t137 = t206 * t390 - t211;
t129 = -pkin(4) * t385 + qJD(5) - t137;
t508 = Ifges(6,4) * t389;
t426 = -Ifges(6,2) * t387 + t508;
t509 = Ifges(6,4) * t387;
t428 = Ifges(6,1) * t389 - t509;
t429 = mrSges(6,1) * t387 + mrSges(6,2) * t389;
t601 = t129 * t429 - t387 * (t223 * Ifges(6,4) + Ifges(6,2) * t439 + t243 * Ifges(6,6)) / 0.2e1 + t264 * mrSges(5,2) - t137 * mrSges(5,3) + t223 * t428 / 0.2e1 + t439 * t426 / 0.2e1;
t133 = -t167 * t387 + t382 * t389;
t134 = t167 * t389 + t382 * t387;
t45 = qJD(6) * t602 + t133 * t392 + t134 * t396;
t559 = t45 / 0.2e1;
t46 = -qJD(6) * t152 + t133 * t396 - t134 * t392;
t558 = t46 / 0.2e1;
t503 = t152 * Ifges(7,4);
t68 = Ifges(7,2) * t602 + t237 * Ifges(7,6) + t503;
t556 = t68 / 0.2e1;
t146 = Ifges(7,4) * t602;
t69 = Ifges(7,1) * t152 + Ifges(7,5) * t237 + t146;
t555 = t69 / 0.2e1;
t553 = t133 / 0.2e1;
t552 = t134 / 0.2e1;
t160 = qJDD(6) + t166;
t547 = t160 / 0.2e1;
t546 = t166 / 0.2e1;
t599 = t333 / 0.2e1;
t534 = t382 / 0.2e1;
t531 = t398 / 0.2e1;
t597 = Ifges(4,5) * t324;
t596 = Ifges(4,6) * t323;
t594 = t398 * Ifges(3,2);
t528 = pkin(2) * t397;
t367 = pkin(3) + t528;
t473 = t390 * t393;
t289 = pkin(2) * t473 + t388 * t367;
t283 = qJ(5) + t289;
t265 = (-pkin(9) - t283) * t387;
t380 = t389 * pkin(9);
t489 = t283 * t389;
t266 = t380 + t489;
t199 = t265 * t396 - t266 * t392;
t452 = pkin(2) * t463;
t453 = pkin(2) * t464;
t287 = -t388 * t453 + t390 * t452;
t280 = qJD(5) + t287;
t492 = t243 * t389;
t437 = pkin(5) * t421 + pkin(9) * t492;
t258 = -t329 * t393 + t307;
t227 = t258 - t497;
t259 = t397 * t329 + t304;
t228 = -t290 + t259;
t164 = t227 * t388 + t228 * t390;
t527 = pkin(3) * t299;
t174 = pkin(4) * t421 + qJ(5) * t243 + t527;
t467 = qJD(1) * t394;
t371 = pkin(2) * t467;
t168 = t174 + t371;
t87 = -t164 * t387 + t389 * t168;
t61 = t437 + t87;
t493 = t243 * t387;
t459 = pkin(9) * t493;
t88 = t389 * t164 + t387 * t168;
t71 = t459 + t88;
t590 = qJD(6) * t199 - t280 * t418 - t392 * t61 - t396 * t71;
t200 = t265 * t392 + t266 * t396;
t589 = -qJD(6) * t200 - t280 * t320 + t392 * t71 - t396 * t61;
t525 = pkin(3) * t388;
t360 = qJ(5) + t525;
t300 = (-pkin(9) - t360) * t387;
t486 = t360 * t389;
t301 = t380 + t486;
t247 = t300 * t392 + t301 * t396;
t145 = t218 * t390 - t211;
t85 = -t145 * t387 + t389 * t174;
t60 = t437 + t85;
t86 = t389 * t145 + t387 * t174;
t70 = t459 + t86;
t588 = -qJD(5) * t320 - qJD(6) * t247 + t392 * t70 - t396 * t60;
t246 = t300 * t396 - t301 * t392;
t587 = -qJD(5) * t418 + qJD(6) * t246 - t392 * t60 - t396 * t70;
t585 = t362 * t570;
t582 = -t164 + t287;
t581 = mrSges(5,1) * t385 + mrSges(6,1) * t439 - mrSges(6,2) * t223 - mrSges(5,3) * t421;
t521 = pkin(5) * t389;
t364 = pkin(4) + t521;
t391 = -pkin(9) - qJ(5);
t420 = -t362 * t391 + t363 * t364;
t269 = t393 * t349 - t397 * t350;
t419 = -t362 * t364 - t363 * t391;
t523 = pkin(4) * t362;
t526 = pkin(3) * t378;
t577 = -m(7) * (t419 - t526) - m(6) * (-t523 - t526) + t605;
t576 = -m(7) * t419 + t605;
t466 = qJD(1) * t398;
t519 = pkin(7) * t398;
t520 = pkin(7) * t394;
t575 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t467) * t519 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t466) * t520;
t483 = t363 * t399;
t574 = t399 * t611 + t483 * t612;
t484 = t363 * t395;
t573 = t395 * t611 + t484 * t612;
t177 = -mrSges(6,2) * t243 + mrSges(6,3) * t439;
t178 = mrSges(6,1) * t243 - mrSges(6,3) * t223;
t572 = t389 * t177 - t387 * t178;
t571 = t315 * t398 + t316 * t394;
t569 = 0.2e1 * t534;
t354 = t362 * mrSges(6,3);
t567 = -t354 + t603 + (t512 - t515 + t430) * t363;
t11 = pkin(9) * t133 + t15;
t6 = pkin(5) * t166 - pkin(9) * t134 + t14;
t2 = qJD(6) * t17 + t11 * t396 + t392 * t6;
t3 = -qJD(6) * t18 - t11 * t392 + t396 * t6;
t566 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t383 = -qJ(4) + t400;
t565 = -m(3) * pkin(7) + m(4) * t400 + m(6) * t383 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - t429;
t564 = -m(5) * t137 + m(6) * t129 - t581;
t346 = -mrSges(3,1) * t398 + mrSges(3,2) * t394;
t563 = -m(4) * t368 - (m(6) * pkin(4) - t430) * t363 - mrSges(2,1) - m(3) * pkin(1) + t346 + t603;
t561 = Ifges(7,4) * t559 + Ifges(7,2) * t558 + Ifges(7,6) * t547;
t560 = Ifges(7,1) * t559 + Ifges(7,4) * t558 + Ifges(7,5) * t547;
t557 = Ifges(6,1) * t552 + Ifges(6,4) * t553 + Ifges(6,5) * t546;
t554 = m(6) + m(7);
t551 = -t602 / 0.2e1;
t550 = t602 / 0.2e1;
t549 = -t152 / 0.2e1;
t548 = t152 / 0.2e1;
t543 = -t237 / 0.2e1;
t542 = t237 / 0.2e1;
t535 = t299 / 0.2e1;
t533 = -t385 / 0.2e1;
t530 = pkin(2) * t393;
t529 = pkin(2) * t394;
t366 = pkin(3) * t379;
t524 = pkin(3) * t390;
t522 = pkin(5) * t387;
t447 = qJD(2) * t400;
t331 = t394 * t447;
t332 = t398 * t447;
t202 = t397 * t331 + t393 * t332 + t349 * t463 + t350 * t464;
t261 = -qJD(2) * t324 - t413;
t161 = qJ(4) * t261 + qJD(4) * t323 + t202;
t203 = -qJD(3) * t269 - t331 * t393 + t397 * t332;
t260 = qJD(2) * t323 + t412;
t162 = -qJ(4) * t260 - qJD(4) * t324 + t203;
t84 = t161 * t390 + t162 * t388;
t197 = t260 * t388 - t390 * t261;
t198 = t260 * t390 + t261 * t388;
t465 = qJD(2) * t394;
t372 = pkin(2) * t465;
t248 = -pkin(3) * t261 + t372;
t255 = t323 * t388 + t324 * t390;
t91 = pkin(4) * t197 - qJ(5) * t198 - qJD(5) * t255 + t248;
t32 = t387 * t91 + t389 * t84;
t511 = Ifges(3,4) * t394;
t510 = Ifges(3,4) * t398;
t502 = t251 * mrSges(4,3);
t501 = t299 * mrSges(4,3);
t500 = t299 * Ifges(4,4);
t79 = mrSges(6,1) * t166 - mrSges(6,3) * t134;
t498 = t387 * t79;
t495 = t198 * t387;
t494 = t198 * t389;
t491 = t255 * t387;
t490 = t255 * t389;
t352 = t362 * qJ(5);
t482 = t375 * t395;
t481 = t375 * t399;
t480 = t376 * t395;
t479 = t376 * t399;
t254 = -t390 * t323 + t324 * t388;
t278 = -pkin(3) * t323 - t368;
t182 = pkin(4) * t254 - qJ(5) * t255 + t278;
t268 = t397 * t349 + t350 * t393;
t240 = -qJ(4) * t324 + t268;
t241 = qJ(4) * t323 + t269;
t184 = t240 * t388 + t241 * t390;
t97 = t387 * t182 + t389 * t184;
t468 = t366 + t381;
t456 = Ifges(7,5) * t45 + Ifges(7,6) * t46 + Ifges(7,3) * t160;
t448 = t363 * pkin(4) + t352 + t366;
t12 = -t46 * mrSges(7,1) + t45 * mrSges(7,2);
t31 = -t387 * t84 + t389 * t91;
t443 = t461 / 0.2e1;
t442 = t166 * mrSges(5,1) + t167 * mrSges(5,2);
t75 = -t133 * mrSges(6,1) + t134 * mrSges(6,2);
t47 = t103 * t390 - t388 * t107;
t83 = t161 * t388 - t390 * t162;
t96 = t389 * t182 - t184 * t387;
t144 = t218 * t388 + t474;
t163 = -t390 * t227 + t228 * t388;
t183 = -t390 * t240 + t241 * t388;
t288 = t367 * t390 - t388 * t530;
t284 = -pkin(4) - t288;
t435 = t366 + t420;
t434 = mrSges(3,1) * t394 + mrSges(3,2) * t398;
t427 = t511 + t594;
t425 = Ifges(3,5) * t398 - Ifges(3,6) * t394;
t424 = Ifges(6,5) * t389 - Ifges(6,6) * t387;
t422 = -t387 * t73 + t389 * t74;
t72 = pkin(5) * t254 - pkin(9) * t490 + t96;
t80 = -pkin(9) * t491 + t97;
t26 = -t392 * t80 + t396 * t72;
t27 = t392 * t72 + t396 * t80;
t415 = pkin(1) * t434;
t414 = t394 * (Ifges(3,1) * t398 - t511);
t38 = -pkin(4) * t382 + qJDD(5) - t47;
t102 = -pkin(5) * t439 + t129;
t238 = t298 * Ifges(4,2) + t385 * Ifges(4,6) + t500;
t292 = Ifges(4,4) * t298;
t239 = t299 * Ifges(4,1) + t385 * Ifges(4,5) + t292;
t28 = -pkin(5) * t133 + t38;
t51 = Ifges(6,4) * t134 + Ifges(6,2) * t133 + Ifges(6,6) * t166;
t401 = (Ifges(4,3) + Ifges(5,3)) * t382 + (-Ifges(5,4) * t539 - Ifges(5,5) * t533 - t424 * t540 + t601) * t243 - (-Ifges(4,2) * t299 + t239 + t292) * t298 / 0.2e1 + (Ifges(7,4) * t320 - Ifges(7,2) * t418) * t558 + (Ifges(7,1) * t320 - Ifges(7,4) * t418) * t559 + (Ifges(7,5) * t320 - Ifges(7,6) * t418) * t547 + t28 * (mrSges(7,1) * t418 + mrSges(7,2) * t320) - t418 * t561 - (-Ifges(5,1) * t243 - t592 + t609) * t421 / 0.2e1 + (-Ifges(7,1) * t294 - Ifges(7,4) * t295) * t548 + (-Ifges(7,4) * t294 - Ifges(7,2) * t295) * t550 + (-Ifges(7,5) * t294 - Ifges(7,6) * t295) * t542 + t348 * (mrSges(4,1) * t299 + mrSges(4,2) * t298) + (-t492 * t73 - t493 * t74 + t423) * mrSges(6,3) - t607 * t555 + (mrSges(7,1) * t608 - mrSges(7,2) * t607) * t102 + (t17 * t607 - t18 * t608 - t2 * t418 - t3 * t320) * mrSges(7,3) - t608 * t556 + (Ifges(7,5) * t549 - Ifges(5,2) * t539 - Ifges(5,6) * t533 + Ifges(7,6) * t551 + Ifges(6,3) * t540 + Ifges(7,3) * t543 + t610) * t421 - Ifges(5,6) * t166 + Ifges(5,5) * t167 + (Ifges(7,1) * t181 + Ifges(7,4) * t180) * t549 + (Ifges(7,4) * t181 + Ifges(7,2) * t180) * t551 + (Ifges(7,5) * t181 + Ifges(7,6) * t180) * t543 + (Ifges(6,1) * t387 + t508) * t552 + (Ifges(6,2) * t389 + t509) * t553 + t387 * t557 + t320 * t560 + (Ifges(4,5) * t298 - Ifges(4,6) * t299) * t533 + t238 * t535 + (Ifges(6,5) * t387 + Ifges(6,6) * t389) * t546 + t252 * t501 + t298 * t502 + t604 * t539 - t299 * (Ifges(4,1) * t298 - t500) / 0.2e1 - t48 * mrSges(5,2) + t47 * mrSges(5,1) + t389 * t51 / 0.2e1 - t172 * mrSges(4,2) + t173 * mrSges(4,1) + Ifges(4,5) * t229 + Ifges(4,6) * t230 + t38 * t430;
t370 = Ifges(3,4) * t466;
t365 = -pkin(4) - t524;
t339 = qJ(5) * t483;
t338 = qJ(5) * t484;
t337 = -t364 - t524;
t335 = -t526 - t529;
t328 = pkin(1) + t468;
t313 = t399 * t335;
t312 = t395 * t335;
t310 = t399 * t328;
t297 = Ifges(3,1) * t467 + Ifges(3,5) * qJD(2) + t370;
t296 = Ifges(3,6) * qJD(2) + qJD(1) * t427;
t277 = t284 - t521;
t276 = mrSges(4,1) * t385 - t501;
t275 = -mrSges(4,2) * t385 + mrSges(4,3) * t298;
t274 = t363 * t479 + t482;
t273 = -t363 * t481 + t480;
t272 = -t363 * t480 + t481;
t271 = t363 * t482 + t479;
t270 = t371 + t527;
t250 = -mrSges(4,1) * t298 + mrSges(4,2) * t299;
t234 = pkin(5) * t493;
t231 = -mrSges(5,2) * t385 - mrSges(5,3) * t243;
t210 = -mrSges(4,2) * t382 + mrSges(4,3) * t230;
t209 = mrSges(4,1) * t382 - mrSges(4,3) * t229;
t194 = t418 * t255;
t193 = t320 * t255;
t192 = mrSges(5,1) * t243 + mrSges(5,2) * t421;
t141 = mrSges(5,1) * t382 - mrSges(5,3) * t167;
t140 = -mrSges(5,2) * t382 - mrSges(5,3) * t166;
t125 = pkin(5) * t491 + t183;
t112 = t163 - t234;
t111 = mrSges(7,1) * t237 - mrSges(7,3) * t152;
t110 = -mrSges(7,2) * t237 + mrSges(7,3) * t602;
t109 = t144 - t234;
t82 = -mrSges(7,1) * t602 + mrSges(7,2) * t152;
t78 = -mrSges(6,2) * t166 + mrSges(6,3) * t133;
t77 = -t198 * t320 + t255 * t294;
t76 = -t198 * t418 - t255 * t295;
t62 = pkin(5) * t495 + t83;
t30 = -mrSges(7,2) * t160 + mrSges(7,3) * t46;
t29 = mrSges(7,1) * t160 - mrSges(7,3) * t45;
t24 = -pkin(9) * t495 + t32;
t19 = pkin(5) * t197 - pkin(9) * t494 + t31;
t5 = -qJD(6) * t27 + t19 * t396 - t24 * t392;
t4 = qJD(6) * t26 + t19 * t392 + t24 * t396;
t1 = [t564 * t83 + (-t272 * mrSges(7,1) - t271 * mrSges(7,2) + (-m(7) * (-t383 + t522) + m(5) * t383 + t565) * t399 + (-m(7) * (-t328 - t420) - m(6) * (-t328 - t352) + t354 + m(5) * t328 - t563) * t395) * g(1) + (-m(6) * t310 - t274 * mrSges(7,1) - t273 * mrSges(7,2) + (-m(7) - m(5)) * (-t383 * t395 + t310) + (-m(7) * t522 + t565) * t395 + (-m(7) * t420 - (m(6) * qJ(5) + mrSges(6,3)) * t362 + t563) * t399) * g(2) + (t596 / 0.2e1 + t597 / 0.2e1) * t382 + (t596 + t597) * t534 + m(6) * (t14 * t96 + t15 * t97 + t31 * t73 + t32 * t74) + m(5) * (t138 * t84 + t184 * t48 + t196 * t278 + t248 * t264) + (t297 * t531 + t425 * qJD(2) / 0.2e1 - t575) * qJD(2) + (Ifges(5,5) * t532 + t424 * t539 + Ifges(5,4) * t540 + t593 / 0.2e1 + t601 + t604 / 0.2e1) * t198 + (-mrSges(3,1) * t520 - mrSges(3,2) * t519 + 0.2e1 * Ifges(3,6) * t531) * qJDD(2) + (Ifges(3,1) * t334 + Ifges(3,4) * t599 + Ifges(3,5) * qJDD(2) - t443 * t594) * t394 + (-t14 * t490 - t15 * t491 - t494 * t73 - t495 * t74) * mrSges(6,3) + (-t17 * t76 + t18 * t77 - t193 * t2 + t194 * t3) * mrSges(7,3) + t28 * (mrSges(7,1) * t193 - mrSges(7,2) * t194) + (-Ifges(7,5) * t194 - Ifges(7,6) * t193) * t547 + (-Ifges(7,1) * t194 - Ifges(7,4) * t193) * t559 + (-Ifges(7,4) * t194 - Ifges(7,2) * t193) * t558 + (Ifges(6,5) * t134 + Ifges(6,6) * t133 + Ifges(6,3) * t166 + t456) * t254 / 0.2e1 + (t398 * t510 + t414) * t443 + (t196 * mrSges(5,2) - t47 * mrSges(5,3) + Ifges(5,1) * t167 - Ifges(5,4) * t166 + Ifges(5,5) * t569 + t38 * t429 + t424 * t546 + t426 * t553 + t428 * t552) * t255 + (t196 * mrSges(5,1) + t14 * mrSges(6,1) - t15 * mrSges(6,2) - t48 * mrSges(5,3) - Ifges(5,4) * t167 + Ifges(6,5) * t552 + Ifges(7,5) * t559 + Ifges(5,2) * t166 - Ifges(5,6) * t569 + Ifges(6,6) * t553 + Ifges(7,6) * t558 + Ifges(6,3) * t546 + Ifges(7,3) * t547 + t566) * t254 + (t333 * t519 + t334 * t520 + t571) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t571) + (-m(5) * t47 + m(6) * t38 - t141 + t75) * t183 + m(4) * (t172 * t269 + t173 * t268 + t202 * t252 + t203 * t251 - t291 * t368 - t348 * t372) + m(7) * (t102 * t62 + t125 * t28 + t17 * t5 + t18 * t4 + t2 * t27 + t26 * t3) + t427 * t599 + t76 * t555 + t77 * t556 + t490 * t557 - t194 * t560 - t193 * t561 + (-mrSges(4,2) * t368 + Ifges(4,1) * t324 + Ifges(4,4) * t323) * t229 + (Ifges(7,5) * t76 + Ifges(7,6) * t77) * t542 + t334 * t510 / 0.2e1 + (mrSges(4,1) * t368 + Ifges(4,4) * t324 + Ifges(4,2) * t323) * t230 + (Ifges(3,4) * t334 + Ifges(3,2) * t333) * t531 + (Ifges(7,1) * t76 + Ifges(7,4) * t77) * t548 + (Ifges(4,5) * t260 + Ifges(4,6) * t261) * t532 + (Ifges(4,1) * t260 + Ifges(4,4) * t261) * t535 + (Ifges(7,4) * t76 + Ifges(7,2) * t77) * t550 + (t172 * t323 - t173 * t324 + t252 * t261) * mrSges(4,3) + Ifges(2,3) * qJDD(1) + t125 * t12 + t4 * t110 + t5 * t111 - t260 * t502 + t102 * (-mrSges(7,1) * t77 + mrSges(7,2) * t76) + t96 * t79 + t97 * t78 + (Ifges(7,6) * t550 + Ifges(6,3) * t539 + Ifges(7,3) * t542 + Ifges(7,5) * t548 + t609 / 0.2e1 - t610 + t613) * t197 - t346 * t496 + t62 * t82 - t51 * t491 / 0.2e1 - t415 * t461 - t296 * t465 / 0.2e1 + t27 * t30 + t26 * t29 + t278 * t442 - t348 * (-mrSges(4,1) * t261 + mrSges(4,2) * t260) + t32 * t177 + t31 * t178 + t184 * t140 + t250 * t372 + t84 * t231 + t248 * t192 + t260 * t239 / 0.2e1 + t261 * t238 / 0.2e1 + t268 * t209 + t269 * t210 + t202 * t275 + t203 * t276 + t298 * (Ifges(4,4) * t260 + Ifges(4,2) * t261) / 0.2e1 + t291 * (-mrSges(4,1) * t323 + mrSges(4,2) * t324) - pkin(1) * (-mrSges(3,1) * t333 + mrSges(3,2) * t334); (t576 * t395 + t573) * g(2) + (t576 * t399 + t574) * g(1) + t581 * t163 + (t137 * t163 + t138 * t582 - t264 * t270 + t288 * t47 + t289 * t48) * m(5) + t582 * t231 + (-m(4) * t381 - m(5) * t468 - m(6) * (t381 + t448) + t346 - m(7) * (t381 + t435) + t567) * g(3) - (-Ifges(3,2) * t467 + t297 + t370) * t466 / 0.2e1 + t570 * (m(4) * t529 - m(5) * t335 + t434 + t606) + t572 * t280 + (t575 + (-t414 / 0.2e1 + t415) * qJD(1)) * qJD(1) + (t452 - t259) * t275 + t401 + ((t172 * t393 + t173 * t397 + (-t251 * t393 + t252 * t397) * qJD(3)) * pkin(2) - t251 * t258 - t252 * t259 + t348 * t371) * m(4) + (-t453 - t258) * t276 + t589 * t111 + (-g(1) * t313 - g(2) * t312 - t102 * t112 + t589 * t17 + t590 * t18 + t199 * t3 + t2 * t200 + t277 * t28) * m(7) + t590 * t110 + (m(7) * t102 + t564 + t82) * (t388 * t397 + t473) * qJD(3) * pkin(2) + t210 * t530 + t209 * t528 + t78 * t489 + (t280 * t422 + t283 * t423 + t284 * t38 - t129 * t163 - t73 * t87 - t74 * t88 - g(2) * (-t395 * t523 + t312 + t338) - g(1) * (-t399 * t523 + t313 + t339)) * m(6) - t112 * t82 - t283 * t498 + Ifges(3,3) * qJDD(2) + t296 * t467 / 0.2e1 - t425 * t461 / 0.2e1 - t250 * t371 - t88 * t177 - t87 * t178 + t199 * t29 + t200 * t30 - t270 * t192 + t277 * t12 + t284 * t75 + t288 * t141 + t289 * t140 - t315 * mrSges(3,2) - t316 * mrSges(3,1) + Ifges(3,6) * t333 + Ifges(3,5) * t334; (t395 * t577 + t573) * g(2) + (t399 * t577 + t574) * g(1) + t581 * t144 + (-m(5) * t366 - m(6) * t448 - m(7) * t435 + t567) * g(3) + (m(5) * t526 + t606) * t570 + t572 * qJD(5) + t401 + (-g(1) * t339 - g(2) * t338 + qJD(5) * t422 - t129 * t144 + t360 * t423 + t365 * t38 - t73 * t85 - t74 * t86) * m(6) + t587 * t110 + t588 * t111 + (-t102 * t109 + t17 * t588 + t18 * t587 + t2 * t247 + t246 * t3 + t28 * t337) * m(7) + t141 * t524 + t140 * t525 + ((t388 * t48 + t390 * t47) * pkin(3) + t137 * t144 - t138 * t145 - t264 * t527) * m(5) + t78 * t486 - t192 * t527 - t109 * t82 - t360 * t498 + t365 * t75 - t86 * t177 - t85 * t178 - t145 * t231 + t246 * t29 + t247 * t30 - t251 * t275 + t252 * t276 + t337 * t12; -t418 * t29 + t320 * t30 + t387 * t78 + t389 * t79 - t608 * t111 - t607 * t110 - (-t231 - t572) * t243 + (-t82 + t581) * t421 + t442 + (-g(1) * t395 + g(2) * t399) * (m(5) + t554) + (-t102 * t421 - t17 * t608 - t18 * t607 + t2 * t320 - t3 * t418) * m(7) + (-t129 * t421 + t14 * t389 + t15 * t387 + t243 * t422) * m(6) + (t137 * t421 + t138 * t243 + t196) * m(5); t554 * t363 * g(3) - t602 * t110 + t152 * t111 - t439 * t177 + t223 * t178 + t12 + t75 + (t152 * t17 - t18 * t602 + t28 - t585) * m(7) + (t223 * t73 - t439 * t74 + t38 - t585) * m(6); -t102 * (mrSges(7,1) * t152 + mrSges(7,2) * t602) + (Ifges(7,1) * t602 - t503) * t549 + t68 * t548 + (Ifges(7,5) * t602 - Ifges(7,6) * t152) * t543 - t17 * t110 + t18 * t111 - g(1) * (mrSges(7,1) * t273 - mrSges(7,2) * t274) - g(2) * (-mrSges(7,1) * t271 + mrSges(7,2) * t272) - g(3) * (-mrSges(7,1) * t375 - mrSges(7,2) * t376) * t362 + (t152 * t18 + t17 * t602) * mrSges(7,3) + t456 + (-Ifges(7,2) * t152 + t146 + t69) * t551 + t566;];
tau  = t1;
