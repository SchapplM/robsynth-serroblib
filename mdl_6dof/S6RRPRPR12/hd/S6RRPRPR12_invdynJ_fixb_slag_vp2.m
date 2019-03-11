% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR12_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR12_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR12_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:17:18
% EndTime: 2019-03-09 11:18:18
% DurationCPUTime: 38.89s
% Computational Cost: add. (16837->944), mult. (39880->1245), div. (0->0), fcn. (30608->14), ass. (0->430)
t407 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t346 = sin(qJ(2));
t341 = sin(pkin(6));
t458 = qJD(1) * t341;
t423 = t346 * t458;
t307 = pkin(2) * t423;
t350 = cos(qJ(2));
t381 = pkin(9) * t346 - qJ(3) * t350;
t211 = t381 * t458 + t307;
t422 = t350 * t458;
t310 = pkin(8) * t422;
t342 = cos(pkin(6));
t457 = qJD(1) * t342;
t445 = pkin(1) * t457;
t258 = t346 * t445 + t310;
t213 = pkin(3) * t422 + t258;
t345 = sin(qJ(4));
t349 = cos(qJ(4));
t129 = -t211 * t345 + t349 * t213;
t452 = qJD(4) * t345;
t544 = pkin(2) + pkin(9);
t463 = qJ(5) + t544;
t470 = t345 * t346;
t625 = -(pkin(4) * t350 - qJ(5) * t470) * t458 - t129 - t349 * qJD(5) + t452 * t463;
t130 = t349 * t211 + t345 * t213;
t400 = t349 * t423;
t410 = t349 * t463;
t624 = qJ(5) * t400 + qJD(4) * t410 + t345 * qJD(5) + t130;
t344 = sin(qJ(6));
t348 = cos(qJ(6));
t290 = qJD(4) + t423;
t340 = sin(pkin(11));
t319 = qJD(2) + t457;
t543 = pkin(3) + pkin(8);
t312 = t350 * t445;
t619 = qJD(3) - t312;
t154 = -t319 * t544 + t423 * t543 + t619;
t413 = -qJ(3) * t346 - pkin(1);
t183 = (-t350 * t544 + t413) * t458;
t102 = t349 * t154 - t183 * t345;
t364 = -t319 * t349 + t345 * t422;
t95 = qJ(5) * t364 + t102;
t88 = pkin(4) * t290 + t95;
t486 = cos(pkin(11));
t103 = t154 * t345 + t183 * t349;
t222 = -t319 * t345 - t349 * t422;
t96 = qJ(5) * t222 + t103;
t93 = t486 * t96;
t40 = t340 * t88 + t93;
t35 = pkin(10) * t290 + t40;
t298 = t319 * qJ(3);
t171 = t298 + t213;
t131 = -pkin(4) * t222 + qJD(5) + t171;
t368 = t340 * t222 - t364 * t486;
t409 = t486 * t222 + t340 * t364;
t63 = -pkin(5) * t409 - pkin(10) * t368 + t131;
t13 = -t344 * t35 + t348 * t63;
t14 = t344 * t63 + t348 * t35;
t145 = qJD(6) - t409;
t536 = -t145 / 0.2e1;
t121 = t290 * t344 + t348 * t368;
t540 = -t121 / 0.2e1;
t120 = t290 * t348 - t344 * t368;
t542 = -t120 / 0.2e1;
t616 = Ifges(7,5) * t540 + Ifges(7,6) * t542 + Ifges(7,3) * t536;
t526 = t290 / 0.2e1;
t531 = t368 / 0.2e1;
t533 = t409 / 0.2e1;
t617 = Ifges(6,4) * t531 + Ifges(6,2) * t533 + Ifges(6,6) * t526;
t623 = mrSges(6,1) * t131 + mrSges(7,1) * t13 - mrSges(7,2) * t14 - t616 - t617;
t592 = -m(7) - m(6);
t593 = m(5) + m(4);
t615 = -t593 + t592;
t570 = qJ(3) * t615 + mrSges(3,2) - mrSges(4,3);
t478 = t340 * t345;
t202 = t400 * t486 - t423 * t478;
t411 = t486 * t349;
t267 = -qJD(4) * t411 + t340 * t452;
t622 = t202 - t267;
t280 = t340 * t349 + t345 * t486;
t203 = t280 * t423;
t268 = t280 * qJD(4);
t577 = t203 + t268;
t527 = -t290 / 0.2e1;
t532 = -t368 / 0.2e1;
t534 = -t409 / 0.2e1;
t621 = Ifges(6,4) * t532 + Ifges(6,2) * t534 + Ifges(6,6) * t527 - t616 + t623;
t535 = t145 / 0.2e1;
t539 = t121 / 0.2e1;
t541 = t120 / 0.2e1;
t620 = -Ifges(7,5) * t539 - Ifges(7,6) * t541 - Ifges(7,3) * t535 + t617 - t623;
t339 = qJ(4) + pkin(11);
t334 = sin(t339);
t335 = cos(t339);
t392 = t345 * mrSges(5,1) + t349 * mrSges(5,2);
t390 = -mrSges(7,1) * t348 + mrSges(7,2) * t344;
t581 = m(7) * pkin(5) + mrSges(6,1) - t390;
t618 = -t334 * t581 - t335 * t407 - t392;
t582 = t340 * t625 - t624 * t486;
t331 = pkin(4) * t349 + pkin(3);
t451 = qJD(4) * t349;
t575 = pkin(4) * t451 - (-pkin(8) - t331) * t423 + t619;
t477 = t341 * t346;
t421 = qJD(2) * t477;
t475 = t341 * t350;
t261 = qJD(1) * t421 - qJDD(1) * t475;
t447 = qJDD(1) * t342;
t317 = qJDD(2) + t447;
t138 = qJD(4) * t222 + t261 * t345 + t317 * t349;
t139 = qJD(4) * t364 + t261 * t349 - t317 * t345;
t86 = -t340 * t138 + t139 * t486;
t549 = t86 / 0.2e1;
t87 = t138 * t486 + t340 * t139;
t548 = t87 / 0.2e1;
t454 = qJD(2) * t350;
t262 = (qJD(1) * t454 + qJDD(1) * t346) * t341;
t245 = qJDD(4) + t262;
t528 = t245 / 0.2e1;
t614 = -mrSges(3,1) + mrSges(4,2);
t613 = pkin(4) * t592 - mrSges(5,1);
t612 = -m(5) * pkin(9) - mrSges(5,3) - mrSges(6,3);
t611 = -pkin(10) * t422 + t582;
t610 = pkin(5) * t622 + t577 * pkin(10) + t575;
t85 = qJDD(6) - t86;
t550 = t85 / 0.2e1;
t46 = -qJD(6) * t121 + t245 * t348 - t344 * t87;
t556 = t46 / 0.2e1;
t45 = qJD(6) * t120 + t245 * t344 + t348 * t87;
t557 = t45 / 0.2e1;
t523 = pkin(1) * t342;
t444 = qJD(2) * t523;
t401 = qJD(1) * t444;
t436 = pkin(1) * t447;
t165 = -pkin(8) * t261 + t346 * t436 + t350 * t401;
t137 = -t317 * qJ(3) - t319 * qJD(3) - t165;
t113 = -pkin(3) * t261 - t137;
t73 = -pkin(4) * t139 + qJDD(5) + t113;
t19 = -pkin(5) * t86 - pkin(10) * t87 + t73;
t320 = pkin(8) * t477;
t166 = -qJD(2) * t310 - qJDD(1) * t320 - t346 * t401 + t350 * t436;
t358 = qJDD(3) - t166;
t112 = pkin(3) * t262 - t317 * t544 + t358;
t453 = qJD(3) * t346;
t356 = -qJ(3) * t262 + (-pkin(1) * qJDD(1) - qJD(1) * t453) * t341;
t115 = t261 * t544 + t356;
t37 = -qJD(4) * t103 + t349 * t112 - t115 * t345;
t23 = pkin(4) * t245 - qJ(5) * t138 + qJD(5) * t364 + t37;
t36 = t345 * t112 + t349 * t115 + t154 * t451 - t183 * t452;
t27 = qJ(5) * t139 + qJD(5) * t222 + t36;
t8 = t340 * t23 + t486 * t27;
t6 = pkin(10) * t245 + t8;
t1 = qJD(6) * t13 + t19 * t344 + t348 * t6;
t2 = -qJD(6) * t14 + t19 * t348 - t344 * t6;
t568 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t9 = Ifges(7,5) * t45 + Ifges(7,6) * t46 + Ifges(7,3) * t85;
t608 = t568 + mrSges(6,1) * t73 + Ifges(7,5) * t557 + Ifges(7,6) * t556 + Ifges(7,3) * t550 + t9 / 0.2e1 + (-t528 - t245 / 0.2e1) * Ifges(6,6) + (-t549 - t86 / 0.2e1) * Ifges(6,2) + (-t548 - t87 / 0.2e1) * Ifges(6,4);
t606 = mrSges(6,3) * t40;
t491 = t340 * t96;
t39 = t486 * t88 - t491;
t604 = t39 * mrSges(6,3);
t603 = t131 * mrSges(6,2);
t584 = t624 * t340 + t486 * t625;
t375 = t612 + t614;
t389 = mrSges(7,1) * t344 + mrSges(7,2) * t348;
t565 = -t375 + t389;
t347 = sin(qJ(1));
t467 = t347 * t350;
t351 = cos(qJ(1));
t468 = t346 * t351;
t273 = t342 * t467 + t468;
t476 = t341 * t347;
t199 = t273 * t349 - t345 * t476;
t109 = mrSges(5,1) * t245 - mrSges(5,3) * t138;
t110 = -mrSges(5,2) * t245 + mrSges(5,3) * t139;
t172 = -mrSges(5,2) * t290 + mrSges(5,3) * t222;
t173 = mrSges(5,1) * t290 + mrSges(5,3) * t364;
t377 = t349 * t172 - t345 * t173;
t598 = t377 * qJD(4) + t349 * t109 + t345 * t110;
t597 = Ifges(6,1) * t531 + Ifges(6,4) * t533 + Ifges(6,5) * t526;
t596 = t73 * mrSges(6,2) + 0.2e1 * Ifges(6,1) * t548 + 0.2e1 * Ifges(6,4) * t549 + 0.2e1 * Ifges(6,5) * t528;
t538 = t138 / 0.2e1;
t537 = t139 / 0.2e1;
t269 = -t342 * t345 - t349 * t475;
t520 = pkin(4) * t269;
t591 = t39 * mrSges(6,1);
t590 = -Ifges(4,4) + Ifges(3,5);
t589 = Ifges(4,5) - Ifges(3,6);
t12 = -mrSges(7,1) * t46 + mrSges(7,2) * t45;
t67 = mrSges(6,1) * t245 - mrSges(6,3) * t87;
t588 = t12 - t67;
t279 = -t411 + t478;
t337 = t345 * pkin(4);
t330 = qJ(3) + t337;
t170 = pkin(5) * t280 + pkin(10) * t279 + t330;
t285 = t463 * t345;
t189 = -t285 * t486 - t340 * t410;
t116 = t170 * t348 - t189 * t344;
t587 = qJD(6) * t116 + t344 * t610 + t348 * t611;
t117 = t170 * t344 + t189 * t348;
t586 = -qJD(6) * t117 - t344 * t611 + t348 * t610;
t127 = mrSges(6,1) * t290 - mrSges(6,3) * t368;
t68 = -mrSges(7,1) * t120 + mrSges(7,2) * t121;
t585 = t127 - t68;
t583 = pkin(5) * t422 - t584;
t153 = -mrSges(5,1) * t222 - mrSges(5,2) * t364;
t404 = mrSges(4,1) * t422;
t253 = -mrSges(4,3) * t319 - t404;
t580 = t153 - t253;
t161 = -t203 * t344 + t348 * t422;
t448 = qJD(6) * t348;
t370 = t344 * t268 + t279 * t448;
t579 = t161 - t370;
t162 = t203 * t348 + t344 * t422;
t449 = qJD(6) * t344;
t369 = -t348 * t268 + t279 * t449;
t578 = t162 - t369;
t522 = pkin(1) * t350;
t425 = -pkin(2) - t522;
t186 = pkin(3) * t477 + t320 + (-pkin(9) + t425) * t342;
t460 = pkin(2) * t475 + qJ(3) * t477;
t205 = (-pkin(9) * t350 - pkin(1)) * t341 - t460;
t125 = t345 * t186 + t349 * t205;
t403 = mrSges(3,3) * t423;
t405 = mrSges(4,1) * t423;
t576 = t319 * t614 + t403 + t405;
t372 = -t342 * t349 + t345 * t475;
t164 = t340 * t269 - t372 * t486;
t300 = t348 * t477;
t143 = -t164 * t344 + t300;
t573 = -t345 * t36 - t349 * t37;
t24 = mrSges(7,1) * t85 - mrSges(7,3) * t45;
t25 = -mrSges(7,2) * t85 + mrSges(7,3) * t46;
t572 = -t344 * t24 + t348 * t25;
t7 = t23 * t486 - t340 * t27;
t567 = t37 * mrSges(5,1) + t7 * mrSges(6,1) - t36 * mrSges(5,2) - t8 * mrSges(6,2) + Ifges(5,5) * t138 + Ifges(6,5) * t87 + Ifges(5,6) * t139 + Ifges(6,6) * t86;
t505 = Ifges(7,4) * t121;
t50 = Ifges(7,2) * t120 + Ifges(7,6) * t145 + t505;
t525 = t344 / 0.2e1;
t91 = Ifges(6,1) * t368 + Ifges(6,4) * t409 + t290 * Ifges(6,5);
t545 = -t91 / 0.2e1;
t119 = Ifges(7,4) * t120;
t51 = t121 * Ifges(7,1) + t145 * Ifges(7,5) + t119;
t552 = -t51 / 0.2e1;
t566 = t348 * t552 + t50 * t525 + t545 - t603;
t563 = t570 + t618;
t562 = t279 * t7 - t280 * t8 + t577 * t39 - t40 * t622;
t10 = t45 * Ifges(7,4) + t46 * Ifges(7,2) + t85 * Ifges(7,6);
t559 = t10 / 0.2e1;
t553 = -t50 / 0.2e1;
t551 = Ifges(5,1) * t538 + Ifges(5,4) * t537 + Ifges(5,5) * t528;
t529 = -t364 / 0.2e1;
t524 = pkin(1) * t341;
t521 = pkin(4) * t364;
t519 = pkin(4) * t340;
t5 = -t245 * pkin(5) - t7;
t518 = t279 * t5;
t512 = -Ifges(4,6) - Ifges(3,4);
t197 = qJD(4) * t269 + t345 * t421;
t420 = t341 * t454;
t309 = pkin(2) * t421;
t176 = t309 + (qJD(2) * t381 - t453) * t341;
t326 = t346 * t523;
t214 = (t475 * t543 + t326) * qJD(2);
t72 = -qJD(4) * t125 - t176 * t345 + t349 * t214;
t52 = pkin(4) * t420 - qJ(5) * t197 + qJD(5) * t372 + t72;
t198 = qJD(4) * t372 + t349 * t421;
t71 = t349 * t176 + t186 * t451 - t205 * t452 + t345 * t214;
t56 = qJ(5) * t198 + qJD(5) * t269 + t71;
t18 = t340 * t52 + t486 * t56;
t511 = mrSges(5,2) * t345;
t510 = mrSges(7,3) * t344;
t509 = mrSges(7,3) * t348;
t508 = Ifges(3,4) * t346;
t507 = Ifges(5,4) * t345;
t506 = Ifges(5,4) * t349;
t504 = Ifges(7,4) * t344;
t503 = Ifges(7,4) * t348;
t502 = Ifges(4,6) * t346;
t501 = Ifges(4,6) * t350;
t500 = t102 * mrSges(5,3);
t499 = t103 * mrSges(5,3);
t498 = t409 * Ifges(6,6);
t497 = t368 * Ifges(6,5);
t496 = t222 * Ifges(5,6);
t495 = t364 * Ifges(5,4);
t494 = t364 * Ifges(5,5);
t124 = t349 * t186 - t205 * t345;
t446 = pkin(4) * t477;
t100 = qJ(5) * t372 + t124 + t446;
t105 = qJ(5) * t269 + t125;
t58 = t340 * t100 + t486 * t105;
t464 = t350 * t351;
t469 = t346 * t347;
t271 = -t342 * t464 + t469;
t484 = t271 * t345;
t482 = t273 * t345;
t480 = t279 * t344;
t479 = t279 * t348;
t474 = t341 * t351;
t277 = pkin(8) * t475 + t326;
t459 = t351 * pkin(1) + pkin(8) * t476;
t455 = qJD(1) ^ 2 * t341 ^ 2;
t443 = Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1;
t442 = Ifges(3,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t441 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t440 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t435 = t344 * t477;
t274 = -t342 * t469 + t464;
t427 = t274 * pkin(2) + t459;
t231 = -t342 * qJ(3) - t277;
t424 = t486 * pkin(4);
t38 = -t86 * mrSges(6,1) + t87 * mrSges(6,2);
t414 = t448 / 0.2e1;
t412 = -pkin(1) * t347 + pkin(8) * t474;
t196 = t262 * mrSges(4,1) + t317 * mrSges(4,2);
t408 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(3,3);
t406 = t543 * t477;
t402 = mrSges(3,3) * t422;
t204 = pkin(3) * t475 - t231;
t272 = t342 * t468 + t467;
t395 = t272 * pkin(2) - t412;
t394 = mrSges(5,1) * t269 + mrSges(5,2) * t372;
t388 = mrSges(4,2) * t350 - mrSges(4,3) * t346;
t387 = Ifges(5,1) * t345 + t506;
t386 = Ifges(7,1) * t348 - t504;
t385 = Ifges(5,2) * t349 + t507;
t384 = -Ifges(7,2) * t344 + t503;
t383 = Ifges(5,5) * t345 + Ifges(5,6) * t349;
t382 = Ifges(7,5) * t348 - Ifges(7,6) * t344;
t380 = t13 * t344 - t14 * t348;
t54 = pkin(10) * t477 + t58;
t157 = t204 - t520;
t163 = -t269 * t486 - t340 * t372;
t77 = pkin(5) * t163 - pkin(10) * t164 + t157;
t29 = t344 * t77 + t348 * t54;
t28 = -t344 * t54 + t348 * t77;
t78 = -mrSges(7,2) * t145 + mrSges(7,3) * t120;
t79 = mrSges(7,1) * t145 - mrSges(7,3) * t121;
t379 = -t344 * t79 + t348 * t78;
t378 = t102 * t345 - t103 * t349;
t257 = pkin(8) * t423 - t312;
t313 = t350 * t444;
t259 = -pkin(8) * t421 + t313;
t126 = -mrSges(6,2) * t290 + mrSges(6,3) * t409;
t374 = -t126 - t379;
t343 = -qJ(5) - pkin(9);
t373 = pkin(4) * t482 - t274 * t343 + t331 * t476 + t427;
t144 = t164 * t348 + t435;
t181 = -t271 * t334 + t335 * t474;
t179 = t271 * t335 + t334 * t474;
t34 = -t290 * pkin(5) - t39;
t371 = t34 * t389;
t17 = -t340 * t56 + t486 * t52;
t57 = t100 * t486 - t340 * t105;
t333 = t342 * qJD(3);
t185 = -qJD(2) * t406 + t313 + t333;
t362 = -g(1) * t274 - g(2) * t272 - g(3) * t477;
t151 = -pkin(2) * t317 + t358;
t357 = t166 * mrSges(3,1) - t165 * mrSges(3,2) + t151 * mrSges(4,2) - t137 * mrSges(4,3);
t128 = -pkin(4) * t198 + t185;
t355 = t1 * t348 - t2 * t344 + (-t13 * t348 - t14 * t344) * qJD(6);
t354 = -qJD(4) * t378 - t573;
t329 = -t424 - pkin(5);
t306 = Ifges(3,4) * t422;
t301 = t349 * t474;
t297 = Ifges(4,1) * t317;
t296 = Ifges(3,3) * t317;
t276 = t342 * t522 - t320;
t275 = (-mrSges(3,1) * t350 + mrSges(3,2) * t346) * t341;
t265 = t273 * pkin(2);
t263 = t271 * pkin(2);
t260 = t277 * qJD(2);
t256 = -qJ(3) * t422 + t307;
t255 = t388 * t458;
t252 = -mrSges(3,2) * t319 + t402;
t244 = Ifges(4,4) * t262;
t243 = Ifges(3,5) * t262;
t242 = Ifges(4,5) * t261;
t241 = Ifges(3,6) * t261;
t239 = t342 * t425 + t320;
t232 = -t460 - t524;
t230 = -t334 * t475 + t335 * t342;
t228 = Ifges(5,3) * t245;
t227 = Ifges(6,3) * t245;
t220 = -t259 - t333;
t219 = Ifges(5,4) * t222;
t218 = (-pkin(2) * t350 + t413) * t458;
t215 = t309 + (-qJ(3) * t454 - t453) * t341;
t212 = -qJD(1) * t406 + t312;
t210 = -t298 - t258;
t209 = t319 * Ifges(4,4) + (-Ifges(4,2) * t346 - t501) * t458;
t208 = t319 * Ifges(4,5) + (-t350 * Ifges(4,3) - t502) * t458;
t207 = Ifges(3,1) * t423 + t319 * Ifges(3,5) + t306;
t206 = t319 * Ifges(3,6) + (t350 * Ifges(3,2) + t508) * t458;
t201 = -pkin(2) * t319 + qJD(3) + t257;
t200 = t349 * t476 + t482;
t195 = mrSges(4,1) * t261 - mrSges(4,3) * t317;
t188 = -t285 * t340 + t410 * t486;
t178 = t273 * t334 + t335 * t476;
t177 = -t273 * t335 + t334 * t476;
t160 = -t202 * t348 + t319 * t344;
t159 = t202 * t344 + t319 * t348;
t142 = pkin(2) * t261 + t356;
t141 = t178 * t348 + t274 * t344;
t140 = -t178 * t344 + t274 * t348;
t134 = -Ifges(5,1) * t364 + t290 * Ifges(5,5) + t219;
t133 = t222 * Ifges(5,2) + t290 * Ifges(5,6) - t495;
t132 = t290 * Ifges(5,3) - t494 + t496;
t123 = t197 * t486 + t340 * t198;
t122 = t197 * t340 - t198 * t486;
t97 = -mrSges(6,1) * t409 + mrSges(6,2) * t368;
t92 = -mrSges(5,1) * t139 + mrSges(5,2) * t138;
t89 = t290 * Ifges(6,3) + t497 + t498;
t76 = pkin(5) * t368 - pkin(10) * t409 - t521;
t75 = -qJD(6) * t144 - t123 * t344 + t348 * t420;
t74 = qJD(6) * t143 + t123 * t348 + t344 * t420;
t69 = t138 * Ifges(5,4) + t139 * Ifges(5,2) + t245 * Ifges(5,6);
t66 = -mrSges(6,2) * t245 + mrSges(6,3) * t86;
t53 = -pkin(5) * t477 - t57;
t47 = pkin(5) * t122 - pkin(10) * t123 + t128;
t42 = t486 * t95 - t491;
t41 = t340 * t95 + t93;
t21 = t344 * t76 + t348 * t42;
t20 = -t344 * t42 + t348 * t76;
t16 = pkin(10) * t420 + t18;
t15 = -pkin(5) * t420 - t17;
t11 = t45 * Ifges(7,1) + t46 * Ifges(7,4) + t85 * Ifges(7,5);
t4 = -qJD(6) * t29 - t16 * t344 + t348 * t47;
t3 = qJD(6) * t28 + t16 * t348 + t344 * t47;
t22 = [(t1 * t143 - t13 * t74 + t14 * t75 - t144 * t2) * mrSges(7,3) + t576 * t260 + (-m(3) * t459 - t200 * mrSges(5,1) - t199 * mrSges(5,2) - mrSges(2,1) * t351 + mrSges(2,2) * t347 - m(7) * (pkin(5) * t178 + t373) - t141 * mrSges(7,1) - t140 * mrSges(7,2) - m(6) * t373 - t178 * mrSges(6,1) + t570 * t273 + t407 * t177 + t375 * t274 - t593 * t427) * g(2) + (Ifges(7,4) * t144 + Ifges(7,2) * t143) * t556 + (Ifges(7,4) * t74 + Ifges(7,2) * t75) * t541 + (-t606 - t620) * t122 + (-t102 * t197 + t103 * t198 + t269 * t36 + t37 * t372) * mrSges(5,3) + (-Ifges(5,4) * t372 + Ifges(5,2) * t269) * t537 + (-Ifges(5,1) * t372 + Ifges(5,4) * t269) * t538 + (-Ifges(5,5) * t372 + Ifges(5,6) * t269) * t528 - t372 * t551 + (Ifges(7,5) * t144 + Ifges(7,6) * t143) * t550 + (Ifges(7,5) * t74 + Ifges(7,6) * t75) * t535 + (t408 * t347 * g(2) + (-mrSges(3,1) * t261 - mrSges(3,2) * t262 + (m(3) * t524 - t275) * qJDD(1)) * pkin(1) + (t408 + t511) * g(1) * t351 + (-t137 * mrSges(4,1) + t142 * mrSges(4,2) + t165 * mrSges(3,3) - t589 * t317 - t512 * t262 + (-Ifges(3,2) - Ifges(4,3)) * t261) * t350 + (t227 / 0.2e1 + t228 / 0.2e1 + (Ifges(4,2) + Ifges(3,1)) * t262 + t512 * t261 + t440 * t245 - t142 * mrSges(4,3) + t151 * mrSges(4,1) - t166 * mrSges(3,3) + t590 * t317 + t567) * t346 + ((-t206 / 0.2e1 + t208 / 0.2e1 - t218 * mrSges(4,2) - t258 * mrSges(3,3) + t210 * mrSges(4,1) + t441 * t319 + (-pkin(1) * mrSges(3,1) - t346 * t443) * t458) * t346 + (t89 / 0.2e1 + t207 / 0.2e1 - t209 / 0.2e1 + t132 / 0.2e1 - t218 * mrSges(4,3) + t591 - t40 * mrSges(6,2) + t102 * mrSges(5,1) - t103 * mrSges(5,2) + t496 / 0.2e1 - t494 / 0.2e1 + t257 * mrSges(3,3) + t201 * mrSges(4,1) + t497 / 0.2e1 + t498 / 0.2e1 + t442 * t319 + t440 * t290 + ((-pkin(1) * mrSges(3,2) + t350 * t443) * t341 + (Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(4,3) / 0.2e1) * t477) * qJD(1)) * t350) * qJD(2)) * t341 + (-t604 + t603 + t91 / 0.2e1 + t597) * t123 + m(3) * (t165 * t277 + t166 * t276 + t257 * t260 + t258 * t259) + m(7) * (t1 * t29 + t13 * t4 + t14 * t3 + t15 * t34 + t2 * t28 + t5 * t53) + m(6) * (t128 * t131 + t157 * t73 + t17 * t39 + t18 * t40 + t57 * t7 + t58 * t8) + m(5) * (t102 * t72 + t103 * t71 + t113 * t204 + t124 * t37 + t125 * t36 + t171 * t185) + m(4) * (t137 * t231 + t142 * t232 + t151 * t239 + t201 * t260 + t210 * t220 + t215 * t218) - t113 * t394 + (t243 / 0.2e1 - t241 / 0.2e1 + t296 / 0.2e1 + t297 / 0.2e1 - t244 / 0.2e1 + t242 / 0.2e1 + (Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1) * t317 + t442 * t262 + t441 * t261 + t357) * t342 + (-mrSges(6,3) * t8 + t608) * t163 + (-t7 * mrSges(6,3) + t596) * t164 + (Ifges(7,1) * t144 + Ifges(7,4) * t143) * t557 + (Ifges(7,1) * t74 + Ifges(7,4) * t75) * t539 + (-t301 * mrSges(5,1) - m(3) * t412 + mrSges(2,1) * t347 + mrSges(2,2) * t351 + t407 * t179 + (t392 - t570) * t271 - t581 * t181 + t565 * t272 + t593 * t395 - t592 * (pkin(4) * t484 - t272 * t343 - t331 * t474 + t395)) * g(1) + t157 * t38 + t5 * (-mrSges(7,1) * t143 + mrSges(7,2) * t144) + t144 * t11 / 0.2e1 + t124 * t109 + t125 * t110 + t18 * t126 + t17 * t127 + t128 * t97 + t222 * (Ifges(5,4) * t197 + Ifges(5,2) * t198) / 0.2e1 + Ifges(2,3) * qJDD(1) + t3 * t78 + t4 * t79 + t74 * t51 / 0.2e1 + t34 * (-mrSges(7,1) * t75 + mrSges(7,2) * t74) + t75 * t50 / 0.2e1 + t58 * t66 + t57 * t67 + t15 * t68 + t53 * t12 + t204 * t92 + t197 * t134 / 0.2e1 + t198 * t133 / 0.2e1 + t171 * (-mrSges(5,1) * t198 + mrSges(5,2) * t197) + t29 * t25 + t28 * t24 + t143 * t559 + t185 * t153 + (Ifges(5,5) * t197 + Ifges(5,6) * t198) * t526 + t231 * t195 + t239 * t196 + t220 * t253 + t215 * t255 + t259 * t252 + t232 * (-mrSges(4,2) * t261 - mrSges(4,3) * t262) + t269 * t69 / 0.2e1 + (Ifges(5,1) * t197 + Ifges(5,4) * t198) * t529 + t276 * (mrSges(3,1) * t317 - mrSges(3,3) * t262) + t277 * (-mrSges(3,2) * t317 - mrSges(3,3) * t261) + t71 * t172 + t72 * t173; (mrSges(7,1) * t579 - mrSges(7,2) * t578) * t34 + (t1 * t480 + t13 * t578 - t14 * t579 + t2 * t479) * mrSges(7,3) + t580 * qJD(3) + t573 * mrSges(5,3) + t575 * t97 + (-m(4) * t201 + t403 - t576) * t258 + t586 * t79 + (t1 * t117 + t116 * t2 + t13 * t586 + t14 * t587 + t188 * t5 + t34 * t583) * m(7) + t587 * t78 + t588 * t188 + t562 * mrSges(6,3) + (t592 * (t273 * t343 + t274 * t337 - t265) + t593 * t265 + t563 * t274 + t565 * t273) * g(1) + (t592 * (t271 * t343 + t272 * t337 - t263) + t593 * t263 + t563 * t272 + t565 * t271) * g(2) + (-pkin(2) * t151 - qJ(3) * t137 - qJD(3) * t210 - t218 * t256) * m(4) + (Ifges(7,1) * t369 + Ifges(7,4) * t370) * t539 + (Ifges(7,1) * t162 + Ifges(7,4) * t161) * t540 + t290 * (mrSges(5,1) * t349 - t511) * t171 + (-t102 * (mrSges(5,1) * t350 - mrSges(5,3) * t470) - t103 * (mrSges(5,3) * t346 * t349 - mrSges(5,2) * t350) - t218 * (-mrSges(4,2) * t346 - mrSges(4,3) * t350)) * t458 + (t275 - t593 * t460 + t592 * (t345 * t446 + t460) + (t388 + (-t343 * t592 - t389 + t612) * t350 + t618 * t346) * t341) * g(3) + (-t346 * (Ifges(3,1) * t350 - t508) / 0.2e1 + pkin(1) * (mrSges(3,1) * t346 + mrSges(3,2) * t350)) * t455 + (-t499 - t133 / 0.2e1) * t451 + (-t131 * t203 + t40 * t422) * mrSges(6,2) + (Ifges(7,5) * t369 + Ifges(7,6) * t370) * t535 + (Ifges(7,5) * t162 + Ifges(7,6) * t161) * t536 + (Ifges(7,4) * t369 + Ifges(7,2) * t370) * t541 + (Ifges(7,4) * t162 + Ifges(7,2) * t161) * t542 + t620 * t267 + t621 * t202 + (Ifges(6,4) * t203 + Ifges(6,6) * t422) * t534 - (t222 * t385 + t290 * t383 - t364 * t387) * qJD(4) / 0.2e1 - ((t349 * t133 + t345 * t134 + t208) * t346 + (-Ifges(3,2) * t423 + t132 + t207 + t306 + t89) * t350 + t290 * (Ifges(5,3) * t350 + t346 * t383) - t364 * (Ifges(5,5) * t350 + t346 * t387) + t222 * (Ifges(5,6) * t350 + t346 * t385) + (t346 * t589 + t350 * t590) * t319) * t458 / 0.2e1 + (t113 * qJ(3) - t102 * t129 - t103 * t130 + (qJD(3) - t212) * t171) * m(5) + (t350 * (Ifges(4,3) * t346 - t501) + t346 * (-Ifges(4,2) * t350 + t502)) * t455 / 0.2e1 + (t346 * t206 + t350 * t209) * t458 / 0.2e1 + (qJD(6) * t51 + t10) * t480 / 0.2e1 + (-t134 / 0.2e1 + t500) * t452 + (Ifges(6,5) * t203 + Ifges(6,3) * t422) * t527 + (-m(5) * t354 - t598) * t544 + t582 * t126 + t583 * t68 + t584 * t127 + (t131 * t575 - t188 * t7 + t189 * t8 + t330 * t73 + t39 * t584 + t40 * t582) * m(6) + t357 - t422 * t591 - t210 * t405 - t244 - t241 + t242 + t243 - t389 * t518 - t11 * t479 / 0.2e1 + t113 * t392 + t608 * t280 + (-t382 * t550 - t384 * t556 - t386 * t557 + t50 * t414 - t596) * t279 + (-m(4) * t210 + t252 - t253 - t402) * t257 + t296 + t297 + (t92 - t195) * qJ(3) + t116 * t24 + t117 * t25 - t212 * t153 - pkin(2) * t196 - t345 * t69 / 0.2e1 + t349 * t551 + t162 * t552 + t161 * t553 + (-Ifges(5,2) * t345 + t506) * t537 + (Ifges(5,1) * t349 - t507) * t538 + t203 * t545 - t201 * t404 + t189 * t66 - t256 * t255 + (Ifges(5,5) * t349 - Ifges(5,6) * t345) * t528 + t330 * t38 + (t566 - t597) * t268 + (Ifges(6,1) * t203 + Ifges(6,5) * t422) * t532 - t130 * t172 - t129 * t173; t196 + (t255 + t377) * t423 + (t66 + (-t344 * t78 - t348 * t79) * qJD(6) + t572) * t280 + t374 * t267 - t159 * t79 - t160 * t78 + t202 * t126 + (-t97 - t580) * t319 + t588 * t279 - t585 * t577 - (-g(1) * t273 - g(2) * t271 + g(3) * t475) * t615 + (-t13 * t159 - t14 * t160 + t380 * t267 + t355 * t280 + t34 * t577 + t518) * m(7) + (-t131 * t319 - t562) * m(6) + (-t171 * t319 - t378 * t423 + t354) * m(5) + (t210 * t319 + t218 * t423 + t151) * m(4) + t598; (m(7) * t355 - t448 * t79 - t449 * t78 + t572) * (pkin(10) + t519) + t567 + (t131 * t521 + t39 * t41 - t40 * t42 + (t340 * t8 + t486 * t7) * pkin(4)) * m(6) + (t606 - t621) * t368 + (mrSges(5,2) * t200 + t581 * t177 + t407 * t178 + t199 * t613) * g(1) + (-(t301 - t484) * mrSges(5,2) - t407 * t181 - t581 * t179 + t613 * (t271 * t349 + t345 * t474)) * g(2) + t364 * (Ifges(5,1) * t222 + t495) / 0.2e1 - (Ifges(5,2) * t364 + t134 + t219) * t222 / 0.2e1 - t171 * (-mrSges(5,1) * t364 + mrSges(5,2) * t222) + (Ifges(5,5) * t222 + Ifges(5,6) * t364) * t527 - t364 * t499 + (t120 * t384 + t121 * t386 + t145 * t382) * qJD(6) / 0.2e1 + (-t13 * t448 - t14 * t449) * mrSges(7,3) + (-t13 * t20 - t14 * t21 + t329 * t5 - t34 * t41) * m(7) + t585 * t41 + (Ifges(6,1) * t532 + Ifges(6,4) * t534 + Ifges(6,5) * t527 + t13 * t509 + t14 * t510 + t382 * t536 + t384 * t542 + t386 * t540 - t371 + t566 + t604) * t409 + t97 * t521 + t227 + t228 + qJD(6) * t371 - t2 * t510 + (-t394 + t407 * t230 - t581 * (-t334 * t342 - t335 * t475) + t592 * t520) * g(3) + t51 * t414 + t5 * t390 + t222 * t500 + t1 * t509 - t42 * t126 - t21 * t78 - t20 * t79 + (Ifges(7,2) * t348 + t504) * t556 + (Ifges(7,1) * t344 + t503) * t557 + t348 * t559 + (Ifges(7,5) * t344 + Ifges(7,6) * t348) * t550 + t449 * t553 + t66 * t519 + t11 * t525 + t133 * t529 + t329 * t12 + t67 * t424 - t102 * t172 + t103 * t173; t348 * t24 + t344 * t25 + t585 * t368 + t379 * qJD(6) + t374 * t409 + t38 + (t1 * t344 - t145 * t380 + t2 * t348 - t368 * t34 + t362) * m(7) + (t368 * t39 - t40 * t409 + t362 + t73) * m(6); -t34 * (mrSges(7,1) * t121 + mrSges(7,2) * t120) + (Ifges(7,1) * t120 - t505) * t540 + t50 * t539 + (Ifges(7,5) * t120 - Ifges(7,6) * t121) * t536 - t13 * t78 + t14 * t79 - g(1) * (mrSges(7,1) * t140 - mrSges(7,2) * t141) - g(2) * ((t181 * t344 + t272 * t348) * mrSges(7,1) + (t181 * t348 - t272 * t344) * mrSges(7,2)) - g(3) * ((-t230 * t344 + t300) * mrSges(7,1) + (-t230 * t348 - t435) * mrSges(7,2)) + (t120 * t13 + t121 * t14) * mrSges(7,3) + t9 + (-Ifges(7,2) * t121 + t119 + t51) * t542 + t568;];
tau  = t22;
