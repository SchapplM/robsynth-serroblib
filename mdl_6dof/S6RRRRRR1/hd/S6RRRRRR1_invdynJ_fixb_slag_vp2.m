% Calculate vector of inverse dynamics joint torques for
% S6RRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:27:52
% EndTime: 2019-03-10 03:28:35
% DurationCPUTime: 24.08s
% Computational Cost: add. (37593->841), mult. (89400->1105), div. (0->0), fcn. (68431->18), ass. (0->407)
t375 = sin(qJ(2));
t381 = cos(qJ(2));
t470 = qJD(1) * qJD(2);
t316 = qJDD(1) * t381 - t375 * t470;
t317 = qJDD(1) * t375 + t381 * t470;
t374 = sin(qJ(3));
t380 = cos(qJ(3));
t306 = -t374 * t375 + t380 * t381;
t403 = t306 * qJD(3);
t210 = qJD(1) * t403 + t316 * t374 + t317 * t380;
t307 = t374 * t381 + t375 * t380;
t404 = t307 * qJD(3);
t211 = -qJD(1) * t404 + t316 * t380 - t317 * t374;
t373 = sin(qJ(4));
t379 = cos(qJ(4));
t291 = t306 * qJD(1);
t292 = t307 * qJD(1);
t438 = t379 * t291 - t292 * t373;
t128 = qJD(4) * t438 + t210 * t379 + t211 * t373;
t238 = t291 * t373 + t292 * t379;
t129 = -qJD(4) * t238 - t210 * t373 + t211 * t379;
t371 = sin(qJ(6));
t377 = cos(qJ(6));
t368 = qJD(2) + qJD(3);
t361 = qJD(4) + t368;
t344 = qJD(5) + t361;
t383 = -pkin(8) - pkin(7);
t329 = t383 * t381;
t313 = qJD(1) * t329;
t293 = t374 * t313;
t328 = t383 * t375;
t312 = qJD(1) * t328;
t299 = qJD(2) * pkin(2) + t312;
t247 = t380 * t299 + t293;
t284 = t292 * pkin(9);
t206 = t247 - t284;
t195 = pkin(3) * t368 + t206;
t296 = t380 * t313;
t248 = t299 * t374 - t296;
t542 = pkin(9) * t291;
t207 = t248 + t542;
t196 = t373 * t207;
t141 = t379 * t195 - t196;
t233 = pkin(10) * t238;
t117 = t141 - t233;
t111 = pkin(4) * t361 + t117;
t372 = sin(qJ(5));
t198 = t379 * t207;
t142 = t195 * t373 + t198;
t541 = pkin(10) * t438;
t118 = t142 + t541;
t378 = cos(qJ(5));
t484 = t378 * t118;
t70 = t111 * t372 + t484;
t68 = pkin(11) * t344 + t70;
t366 = t381 * pkin(2);
t353 = t366 + pkin(1);
t327 = t353 * qJD(1);
t260 = -pkin(3) * t291 - t327;
t184 = -pkin(4) * t438 + t260;
t619 = t378 * t238 + t372 * t438;
t620 = -t238 * t372 + t378 * t438;
t92 = -pkin(5) * t620 - pkin(11) * t619 + t184;
t27 = -t371 * t68 + t377 * t92;
t28 = t371 * t92 + t377 * t68;
t57 = qJD(5) * t620 + t128 * t378 + t129 * t372;
t58 = -qJD(5) * t619 - t128 * t372 + t129 * t378;
t503 = qJDD(1) * pkin(1);
t282 = -pkin(2) * t316 - t503;
t183 = -pkin(3) * t211 + t282;
t95 = -pkin(4) * t129 + t183;
t19 = -pkin(5) * t58 - pkin(11) * t57 + t95;
t367 = qJDD(2) + qJDD(3);
t360 = qJDD(4) + t367;
t305 = t317 * pkin(7);
t258 = qJDD(2) * pkin(2) - pkin(8) * t317 - t305;
t304 = t316 * pkin(7);
t259 = pkin(8) * t316 + t304;
t154 = -qJD(3) * t248 + t380 * t258 - t259 * t374;
t115 = pkin(3) * t367 - pkin(9) * t210 + t154;
t477 = qJD(3) * t380;
t478 = qJD(3) * t374;
t153 = t374 * t258 + t380 * t259 + t299 * t477 + t313 * t478;
t123 = pkin(9) * t211 + t153;
t47 = -qJD(4) * t142 + t379 * t115 - t123 * t373;
t24 = pkin(4) * t360 - pkin(10) * t128 + t47;
t475 = qJD(4) * t379;
t476 = qJD(4) * t373;
t46 = t373 * t115 + t379 * t123 + t195 * t475 - t207 * t476;
t26 = pkin(10) * t129 + t46;
t473 = qJD(5) * t378;
t474 = qJD(5) * t372;
t10 = t111 * t473 - t118 * t474 + t372 * t24 + t378 * t26;
t342 = qJDD(5) + t360;
t7 = pkin(11) * t342 + t10;
t3 = -qJD(6) * t28 + t19 * t377 - t371 * t7;
t11 = -qJD(5) * t70 + t24 * t378 - t26 * t372;
t157 = t344 * t377 - t371 * t619;
t41 = qJD(6) * t157 + t342 * t371 + t377 * t57;
t158 = t344 * t371 + t377 * t619;
t42 = -qJD(6) * t158 + t342 * t377 - t371 * t57;
t56 = qJDD(6) - t58;
t14 = t41 * Ifges(7,4) + t42 * Ifges(7,2) + t56 * Ifges(7,6);
t2 = qJD(6) * t27 + t19 * t371 + t377 * t7;
t172 = qJD(6) - t620;
t423 = Ifges(7,5) * t377 - Ifges(7,6) * t371;
t406 = t172 * t423;
t522 = Ifges(7,4) * t371;
t427 = Ifges(7,1) * t377 - t522;
t407 = t158 * t427;
t521 = Ifges(7,4) * t377;
t425 = -Ifges(7,2) * t371 + t521;
t408 = t157 * t425;
t428 = mrSges(7,1) * t371 + mrSges(7,2) * t377;
t501 = t118 * t372;
t69 = t111 * t378 - t501;
t67 = -pkin(5) * t344 - t69;
t409 = t67 * t428;
t472 = qJD(6) * t371;
t445 = -t472 / 0.2e1;
t156 = Ifges(7,4) * t157;
t83 = Ifges(7,1) * t158 + Ifges(7,5) * t172 + t156;
t505 = t377 * t83;
t455 = t505 / 0.2e1;
t527 = mrSges(7,3) * t377;
t572 = t56 / 0.2e1;
t573 = t42 / 0.2e1;
t574 = t41 / 0.2e1;
t575 = Ifges(7,1) * t574 + Ifges(7,4) * t573 + Ifges(7,5) * t572;
t530 = mrSges(7,2) * t371;
t533 = mrSges(7,1) * t377;
t624 = t530 - t533;
t8 = -pkin(5) * t342 - t11;
t513 = t158 * Ifges(7,4);
t82 = t157 * Ifges(7,2) + t172 * Ifges(7,6) + t513;
t386 = -t10 * mrSges(6,2) + t2 * t527 + t371 * t575 + t377 * t14 / 0.2e1 + Ifges(6,3) * t342 + (Ifges(7,1) * t371 + t521) * t574 + (Ifges(7,2) * t377 + t522) * t573 + (Ifges(7,5) * t371 + Ifges(7,6) * t377) * t572 + Ifges(6,6) * t58 + Ifges(6,5) * t57 + t8 * t624 + t82 * t445 + t11 * mrSges(6,1) + (t409 + t455) * qJD(6) + (t408 + t407 + t406) * qJD(6) / 0.2e1;
t515 = t141 * mrSges(5,3);
t528 = mrSges(7,3) * t371;
t557 = -t344 / 0.2e1;
t563 = -t619 / 0.2e1;
t564 = -t620 / 0.2e1;
t565 = -t172 / 0.2e1;
t567 = -t158 / 0.2e1;
t568 = -t157 / 0.2e1;
t518 = Ifges(6,6) * t344;
t523 = Ifges(6,4) * t619;
t102 = Ifges(6,2) * t620 + t518 + t523;
t535 = t70 * mrSges(6,3);
t609 = t28 * mrSges(7,2);
t610 = t27 * mrSges(7,1);
t516 = Ifges(7,3) * t172;
t517 = Ifges(7,6) * t157;
t519 = Ifges(7,5) * t158;
t81 = t516 + t517 + t519;
t578 = t609 - t610 - t184 * mrSges(6,1) + t535 + t102 / 0.2e1 - t81 / 0.2e1;
t171 = Ifges(6,4) * t620;
t520 = Ifges(6,5) * t344;
t103 = Ifges(6,1) * t619 + t171 + t520;
t507 = t371 * t82;
t536 = t69 * mrSges(6,3);
t579 = -t184 * mrSges(6,2) - t409 - t505 / 0.2e1 + t507 / 0.2e1 + t536 - t103 / 0.2e1;
t471 = qJD(6) * t377;
t595 = -t27 * t471 - t28 * t472;
t633 = (-Ifges(6,4) * t563 + Ifges(7,5) * t567 - Ifges(6,2) * t564 - Ifges(6,6) * t557 + Ifges(7,6) * t568 + Ifges(7,3) * t565 + t578) * t619 + t47 * mrSges(5,1) - t46 * mrSges(5,2) + mrSges(7,3) * t595 + Ifges(5,5) * t128 + Ifges(5,6) * t129 + Ifges(5,3) * t360 - t3 * t528 + t438 * t515 + t386 + (Ifges(6,1) * t563 + Ifges(6,4) * t564 + Ifges(6,5) * t557 + t27 * t527 + t28 * t528 + t423 * t565 + t425 * t568 + t427 * t567 + t579) * t620;
t370 = qJ(2) + qJ(3);
t364 = qJ(4) + t370;
t356 = qJ(5) + t364;
t340 = sin(t356);
t341 = cos(t356);
t632 = t340 * t530 + t341 * (m(7) * pkin(11) + mrSges(7,3));
t603 = mrSges(6,1) * t344 + mrSges(7,1) * t157 - mrSges(7,2) * t158 - mrSges(6,3) * t619;
t630 = -m(6) * t69 + m(7) * t67 - t603;
t109 = pkin(5) * t619 - pkin(11) * t620;
t597 = t341 * pkin(5) + t340 * pkin(11);
t628 = m(7) * t597;
t549 = pkin(4) * t238;
t253 = -t312 * t374 + t296;
t212 = t253 - t542;
t254 = t380 * t312 + t293;
t213 = -t284 + t254;
t149 = t373 * t212 + t379 * t213;
t130 = -t233 + t149;
t553 = pkin(2) * t380;
t352 = pkin(3) + t553;
t489 = t373 * t374;
t286 = -pkin(2) * t489 + t379 * t352;
t280 = pkin(4) + t286;
t487 = t374 * t379;
t288 = pkin(2) * t487 + t352 * t373;
t223 = t372 * t280 + t378 * t288;
t242 = t352 * t475 + (-t374 * t476 + (t379 * t380 - t489) * qJD(3)) * pkin(2);
t243 = -t352 * t476 + (-t374 * t475 + (-t373 * t380 - t487) * qJD(3)) * pkin(2);
t148 = t379 * t212 - t213 * t373;
t413 = t148 - t541;
t604 = qJD(5) * t223 + (-t243 + t413) * t378 + (-t130 + t242) * t372;
t230 = Ifges(5,4) * t438;
t168 = t238 * Ifges(5,1) + t361 * Ifges(5,5) + t230;
t627 = t168 + t230;
t626 = -t341 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t340;
t347 = sin(t364);
t348 = cos(t364);
t531 = mrSges(6,2) * t341;
t625 = mrSges(5,1) * t347 + mrSges(6,1) * t340 + mrSges(5,2) * t348 + t531;
t376 = sin(qJ(1));
t382 = cos(qJ(1));
t623 = g(1) * t382 + g(2) * t376;
t362 = sin(t370);
t363 = cos(t370);
t622 = mrSges(4,1) * t362 + mrSges(4,2) * t363 + t625;
t613 = -m(7) - m(6);
t524 = Ifges(5,4) * t238;
t167 = Ifges(5,2) * t438 + t361 * Ifges(5,6) + t524;
t612 = t167 / 0.2e1;
t611 = t316 / 0.2e1;
t555 = t381 / 0.2e1;
t562 = -t438 / 0.2e1;
t608 = mrSges(5,2) * t438;
t607 = Ifges(5,1) * t438;
t606 = Ifges(5,5) * t438;
t605 = t381 * Ifges(3,2);
t261 = t380 * t328 + t329 * t374;
t231 = -pkin(9) * t307 + t261;
t262 = t374 * t328 - t380 * t329;
t232 = pkin(9) * t306 + t262;
t170 = t373 * t231 + t379 * t232;
t599 = t242 - t149;
t598 = t243 - t148;
t442 = t348 * mrSges(5,1) - t347 * mrSges(5,2);
t443 = t363 * mrSges(4,1) - mrSges(4,2) * t362;
t596 = t341 * t624 + t626;
t480 = qJD(1) * t381;
t481 = qJD(1) * t375;
t543 = pkin(7) * t381;
t544 = pkin(7) * t375;
t594 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t481) * t543 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t480) * t544;
t593 = t304 * t381 + t305 * t375;
t590 = -t442 + t596;
t100 = mrSges(7,1) * t172 - mrSges(7,3) * t158;
t159 = -mrSges(6,2) * t344 + mrSges(6,3) * t620;
t99 = -mrSges(7,2) * t172 + mrSges(7,3) * t157;
t589 = -t100 * t371 + t377 * t99 + t159;
t586 = -t443 + t590;
t20 = mrSges(7,1) * t56 - mrSges(7,3) * t41;
t21 = -mrSges(7,2) * t56 + mrSges(7,3) * t42;
t585 = -t100 * t471 - t371 * t20 + t377 * t21 - t99 * t472;
t584 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t369 = -pkin(9) + t383;
t583 = -m(3) * pkin(7) + m(4) * t383 + m(5) * t369 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t326 = -mrSges(3,1) * t381 + mrSges(3,2) * t375;
t346 = pkin(3) * t363;
t482 = t346 + t366;
t582 = mrSges(2,1) + m(5) * (pkin(1) + t482) + t442 + m(4) * t353 + t443 + m(3) * pkin(1) - t326 - t626;
t421 = t27 * t377 + t28 * t371;
t537 = t3 * t371;
t391 = -qJD(6) * t421 - t537;
t538 = t2 * t377;
t581 = m(7) * (t391 + t538) + t585;
t420 = -t27 * t371 + t28 * t377;
t580 = m(6) * t70 + m(7) * t420 + t589;
t576 = m(7) * pkin(5);
t566 = t158 / 0.2e1;
t561 = -t238 / 0.2e1;
t560 = t238 / 0.2e1;
t558 = t292 / 0.2e1;
t556 = -t361 / 0.2e1;
t554 = pkin(2) * t375;
t552 = pkin(3) * t292;
t551 = pkin(3) * t362;
t550 = pkin(3) * t379;
t548 = pkin(4) * t347;
t338 = pkin(4) * t348;
t547 = pkin(4) * t372;
t546 = pkin(4) * t378;
t545 = pkin(5) * t340;
t529 = mrSges(4,3) * t291;
t526 = Ifges(3,4) * t375;
t525 = Ifges(3,4) * t381;
t514 = t142 * mrSges(5,3);
t512 = t248 * mrSges(4,3);
t511 = t292 * Ifges(4,4);
t251 = t306 * t379 - t307 * t373;
t252 = t306 * t373 + t307 * t379;
t182 = t251 * t372 + t252 * t378;
t500 = t182 * t371;
t499 = t182 * t377;
t492 = t371 * t376;
t491 = t371 * t382;
t490 = t372 * t373;
t488 = t373 * t378;
t486 = t376 * t377;
t485 = t377 * t382;
t146 = t379 * t206 - t196;
t351 = pkin(4) + t550;
t287 = pkin(3) * t488 + t372 * t351;
t479 = qJD(2) * t375;
t468 = Ifges(7,5) * t41 + Ifges(7,6) * t42 + Ifges(7,3) * t56;
t358 = pkin(2) * t479;
t462 = t340 * t533;
t357 = pkin(2) * t481;
t454 = t338 + t597;
t453 = t338 + t482;
t452 = qJD(2) * t383;
t451 = t182 * t471;
t256 = -qJD(2) * t307 - t404;
t239 = -pkin(3) * t256 + t358;
t444 = t470 / 0.2e1;
t145 = -t206 * t373 - t198;
t169 = t379 * t231 - t232 * t373;
t437 = t346 + t454;
t268 = -pkin(3) * t306 - t353;
t436 = t632 * t376;
t435 = t632 * t382;
t192 = t552 + t549;
t300 = -t548 - t551;
t433 = mrSges(3,1) * t375 + mrSges(3,2) * t381;
t426 = t526 + t605;
t424 = Ifges(3,5) * t381 - Ifges(3,6) * t375;
t137 = -pkin(10) * t252 + t169;
t138 = pkin(10) * t251 + t170;
t91 = t137 * t372 + t138 * t378;
t199 = -pkin(4) * t251 + t268;
t416 = t378 * t251 - t252 * t372;
t97 = -pkin(5) * t416 - pkin(11) * t182 + t199;
t44 = t371 * t97 + t377 * t91;
t43 = -t371 * t91 + t377 * t97;
t418 = t378 * t137 - t138 * t372;
t222 = t280 * t378 - t288 * t372;
t255 = qJD(2) * t306 + t403;
t151 = -qJD(4) * t252 - t255 * t373 + t256 * t379;
t131 = -pkin(4) * t151 + t239;
t285 = -pkin(3) * t490 + t351 * t378;
t414 = t145 - t541;
t412 = pkin(1) * t433;
t150 = qJD(4) * t251 + t255 * t379 + t256 * t373;
t79 = qJD(5) * t416 + t150 * t378 + t151 * t372;
t411 = t371 * t79 + t451;
t410 = t182 * t472 - t377 * t79;
t405 = t375 * (Ifges(3,1) * t381 - t526);
t314 = t375 * t452;
t315 = t381 * t452;
t188 = t380 * t314 + t374 * t315 + t328 * t477 + t329 * t478;
t164 = pkin(9) * t256 + t188;
t189 = -qJD(3) * t262 - t314 * t374 + t380 * t315;
t165 = -pkin(9) * t255 + t189;
t84 = t379 * t164 + t373 * t165 + t231 * t475 - t232 * t476;
t277 = t300 - t554;
t393 = m(7) * (t277 - t545) - t462;
t392 = m(7) * (t300 - t545) - t462;
t94 = t109 + t192;
t390 = m(7) * (-t545 - t548) - t462;
t389 = t531 + (mrSges(6,1) + t533 + t576) * t340;
t85 = -qJD(4) * t170 - t164 * t373 + t379 * t165;
t387 = -pkin(10) * t150 + t85;
t228 = t291 * Ifges(4,2) + t368 * Ifges(4,6) + t511;
t283 = Ifges(4,4) * t291;
t229 = t292 * Ifges(4,1) + t368 * Ifges(4,5) + t283;
t384 = -t292 * (Ifges(4,1) * t291 - t511) / 0.2e1 + t228 * t558 + t292 * t512 + t327 * (mrSges(4,1) * t292 + mrSges(4,2) * t291) + t247 * t529 + Ifges(4,5) * t210 + Ifges(4,6) * t211 - t153 * mrSges(4,2) + t154 * mrSges(4,1) + t606 * t556 + t607 * t561 - t260 * t608 + Ifges(4,3) * t367 - t368 * (Ifges(4,5) * t291 - Ifges(4,6) * t292) / 0.2e1 + t627 * t562 - (-Ifges(4,2) * t292 + t229 + t283) * t291 / 0.2e1 + (-t260 * mrSges(5,1) - Ifges(5,4) * t561 - Ifges(5,2) * t562 - Ifges(5,6) * t556 + t514 + t612) * t238 + t633;
t365 = -pkin(10) + t369;
t355 = Ifges(3,4) * t480;
t350 = -pkin(5) - t546;
t290 = Ifges(3,1) * t481 + Ifges(3,5) * qJD(2) + t355;
t289 = Ifges(3,6) * qJD(2) + qJD(1) * t426;
t279 = -pkin(5) - t285;
t275 = pkin(1) + t453;
t273 = t341 * t485 + t492;
t272 = -t341 * t491 + t486;
t271 = -t341 * t486 + t491;
t270 = t341 * t492 + t485;
t265 = mrSges(4,1) * t368 - mrSges(4,3) * t292;
t264 = -mrSges(4,2) * t368 + t529;
t263 = t357 + t552;
t246 = -mrSges(4,1) * t291 + mrSges(4,2) * t292;
t220 = -pkin(5) - t222;
t215 = mrSges(5,1) * t361 - mrSges(5,3) * t238;
t214 = -mrSges(5,2) * t361 + mrSges(5,3) * t438;
t194 = -mrSges(4,2) * t367 + mrSges(4,3) * t211;
t193 = mrSges(4,1) * t367 - mrSges(4,3) * t210;
t185 = t192 + t357;
t180 = -mrSges(5,1) * t438 + mrSges(5,2) * t238;
t124 = -t233 + t146;
t120 = -mrSges(5,2) * t360 + mrSges(5,3) * t129;
t119 = mrSges(5,1) * t360 - mrSges(5,3) * t128;
t108 = -mrSges(6,1) * t620 + mrSges(6,2) * t619;
t96 = t109 + t549;
t93 = t357 + t94;
t80 = qJD(5) * t182 + t150 * t372 - t378 * t151;
t76 = t378 * t130 + t372 * t413;
t74 = t378 * t124 + t372 * t414;
t72 = t117 * t378 - t501;
t71 = t117 * t372 + t484;
t59 = pkin(10) * t151 + t84;
t51 = -mrSges(6,2) * t342 + mrSges(6,3) * t58;
t50 = mrSges(6,1) * t342 - mrSges(6,3) * t57;
t38 = t109 * t371 + t377 * t69;
t37 = t109 * t377 - t371 * t69;
t34 = t371 * t96 + t377 * t72;
t33 = -t371 * t72 + t377 * t96;
t32 = t371 * t93 + t377 * t76;
t31 = -t371 * t76 + t377 * t93;
t30 = t371 * t94 + t377 * t74;
t29 = -t371 * t74 + t377 * t94;
t22 = pkin(5) * t80 - pkin(11) * t79 + t131;
t18 = -mrSges(7,1) * t42 + mrSges(7,2) * t41;
t16 = qJD(5) * t418 + t372 * t387 + t378 * t59;
t5 = -qJD(6) * t44 - t16 * t371 + t22 * t377;
t4 = qJD(6) * t43 + t16 * t377 + t22 * t371;
t1 = [m(4) * (t153 * t262 + t154 * t261 + t188 * t248 + t189 * t247 - t282 * t353 - t327 * t358) - (t468 / 0.2e1 + Ifges(7,3) * t572 + Ifges(7,6) * t573 + Ifges(7,5) * t574 - Ifges(6,4) * t57 + t95 * mrSges(6,1) - Ifges(6,2) * t58 - Ifges(6,6) * t342 - t10 * mrSges(6,3) + t584) * t416 + (-t271 * mrSges(7,1) - t270 * mrSges(7,2) + (-t365 * t613 + t583) * t382 + (-m(7) * (-t275 - t597) + m(6) * t275 + t582) * t376) * g(1) + (-t273 * mrSges(7,1) - t272 * mrSges(7,2) + t613 * (t382 * t275 - t365 * t376) + t583 * t376 + (-t582 - t628) * t382) * g(2) - t82 * t451 / 0.2e1 + (Ifges(3,4) * t317 + Ifges(3,2) * t316) * t555 - t80 * t535 - t79 * t536 - t150 * t515 + (Ifges(4,1) * t255 + Ifges(4,4) * t256) * t558 + (Ifges(5,1) * t150 + Ifges(5,4) * t151) * t560 + (-Ifges(7,1) * t410 - Ifges(7,4) * t411 + Ifges(7,5) * t80) * t566 + t79 * t455 + t67 * (mrSges(7,1) * t411 - mrSges(7,2) * t410) + t157 * (-Ifges(7,4) * t410 - Ifges(7,2) * t411 + Ifges(7,6) * t80) / 0.2e1 + t172 * (-Ifges(7,5) * t410 - Ifges(7,6) * t411 + Ifges(7,3) * t80) / 0.2e1 + t256 * t512 + t255 * t229 / 0.2e1 + t256 * t228 / 0.2e1 + t260 * (-mrSges(5,1) * t151 + mrSges(5,2) * t150) + t261 * t193 + t262 * t194 - t79 * t507 / 0.2e1 + (-mrSges(4,1) * t282 + mrSges(4,3) * t153 + Ifges(4,4) * t210 + Ifges(4,2) * t211 + Ifges(4,6) * t367) * t306 + (t381 * t525 + t405) * t444 - t326 * t503 + m(5) * (t141 * t85 + t142 * t84 + t169 * t47 + t170 * t46 + t183 * t268 + t239 * t260) + (-t2 * t500 + t27 * t410 - t28 * t411 - t3 * t499) * mrSges(7,3) - (-m(6) * t11 + m(7) * t8 + t18 - t50) * t418 - t412 * t470 + (mrSges(4,2) * t282 - mrSges(4,3) * t154 + Ifges(4,1) * t210 + Ifges(4,4) * t211 + Ifges(4,5) * t367) * t307 + t317 * t525 / 0.2e1 + (t95 * mrSges(6,2) - t11 * mrSges(6,3) + Ifges(6,1) * t57 + Ifges(6,4) * t58 + Ifges(6,5) * t342 + t423 * t572 + t425 * t573 + t427 * t574 + t428 * t8 + t445 * t83) * t182 + t619 * (Ifges(6,1) * t79 - Ifges(6,4) * t80) / 0.2e1 + t620 * (Ifges(6,4) * t79 - Ifges(6,2) * t80) / 0.2e1 + t151 * t514 + (mrSges(5,2) * t183 - mrSges(5,3) * t47 + Ifges(5,1) * t128 + Ifges(5,4) * t129 + Ifges(5,5) * t360) * t252 + m(6) * (t10 * t91 + t131 * t184 + t16 * t70 + t199 * t95) + m(7) * (t2 * t44 + t27 * t5 + t28 * t4 + t3 * t43) + (-mrSges(5,1) * t183 + mrSges(5,3) * t46 + Ifges(5,4) * t128 + Ifges(5,2) * t129 + Ifges(5,6) * t360) * t251 + t499 * t575 - t289 * t479 / 0.2e1 + t239 * t180 + t84 * t214 + t85 * t215 + Ifges(2,3) * qJDD(1) + t199 * (-mrSges(6,1) * t58 + mrSges(6,2) * t57) + t184 * (mrSges(6,1) * t80 + mrSges(6,2) * t79) + t150 * t168 / 0.2e1 + t169 * t119 + t170 * t120 + t16 * t159 + t131 * t108 + t438 * (Ifges(5,4) * t150 + Ifges(5,2) * t151) / 0.2e1 - t80 * t102 / 0.2e1 + t79 * t103 / 0.2e1 + t4 * t99 + t5 * t100 + t91 * t51 + t80 * t81 / 0.2e1 + t43 * t20 + t44 * t21 - t247 * t255 * mrSges(4,3) + t188 * t264 + t189 * t265 + t268 * (-mrSges(5,1) * t129 + mrSges(5,2) * t128) + t291 * (Ifges(4,4) * t255 + Ifges(4,2) * t256) / 0.2e1 + t630 * (qJD(5) * t91 + t372 * t59 - t378 * t387) + t246 * t358 - pkin(1) * (-mrSges(3,1) * t316 + mrSges(3,2) * t317) - t327 * (-mrSges(4,1) * t256 + mrSges(4,2) * t255) - t80 * t609 + (-mrSges(3,1) * t544 - mrSges(3,2) * t543 + 0.2e1 * Ifges(3,6) * t555) * qJDD(2) + t344 * (Ifges(6,5) * t79 - Ifges(6,6) * t80) / 0.2e1 + (Ifges(3,1) * t317 + Ifges(3,4) * t611 + Ifges(3,5) * qJDD(2) - t444 * t605) * t375 - t353 * (-mrSges(4,1) * t211 + mrSges(4,2) * t210) + t361 * (Ifges(5,5) * t150 + Ifges(5,6) * t151) / 0.2e1 + t80 * t610 + t426 * t611 + t151 * t612 + t368 * (Ifges(4,5) * t255 + Ifges(4,6) * t256) / 0.2e1 - t14 * t500 / 0.2e1 + (t316 * t543 + t317 * t544 + t593) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t593) + (t290 * t555 + t424 * qJD(2) / 0.2e1 - t594) * qJD(2); -m(4) * (t247 * t253 + t248 * t254 - t327 * t357) - t604 * t603 + t623 * (m(4) * t554 - m(6) * t277 - m(5) * (-t551 - t554) + t433 + t622) + (t594 + (-t405 / 0.2e1 + t412) * qJD(1)) * qJD(1) + t580 * (qJD(5) * t222 + t242 * t378 + t243 * t372) + t581 * (pkin(11) + t223) - t246 * t357 + (t264 * t477 + m(4) * (t153 * t374 + t154 * t380 + (-t247 * t374 + t248 * t380) * qJD(3)) + t374 * t194 - t265 * t478) * pkin(2) + t384 - t424 * t470 / 0.2e1 + (-m(7) * (t366 + t437) - m(4) * t366 - m(5) * t482 - m(6) * t453 + t326 + t586) * g(3) + t193 * t553 + t598 * t215 + t599 * t214 + (t141 * t598 + t142 * t599 - t260 * t263 + t286 * t47 + t288 * t46) * m(5) + t289 * t481 / 0.2e1 + t220 * t18 + t222 * t50 + t223 * t51 - t185 * t108 - t76 * t159 - (-Ifges(3,2) * t481 + t290 + t355) * t480 / 0.2e1 - t32 * t99 - t31 * t100 + Ifges(3,3) * qJDD(2) - g(1) * (t382 * t393 + t435) - g(2) * (t376 * t393 + t436) - t263 * t180 - t254 * t264 - t253 * t265 + (t220 * t8 - t27 * t31 - t28 * t32 + t604 * t67) * m(7) + (t10 * t223 + t11 * t222 - t184 * t185 - t604 * t69 - t70 * t76) * m(6) + t286 * t119 + t288 * t120 - t304 * mrSges(3,2) - t305 * mrSges(3,1) + Ifges(3,6) * t316 + Ifges(3,5) * t317; t623 * (m(5) * t551 - m(6) * t300 + t622) + t580 * (t351 * t473 + (-t373 * t474 + (t378 * t379 - t490) * qJD(4)) * pkin(3)) + t581 * (pkin(11) + t287) + (-t27 * t29 + t279 * t8 - t28 * t30) * m(7) + (t10 * t287 + t11 * t285 - t184 * t192 - t70 * t74) * m(6) + (m(5) * (t373 * t46 + t379 * t47 + (-t141 * t373 + t142 * t379) * qJD(4)) + t373 * t120 + t214 * t475 - t215 * t476) * pkin(3) - t180 * t552 - m(5) * (t141 * t145 + t142 * t146 + t260 * t552) + t384 + (-m(7) * t437 - m(5) * t346 - m(6) * (t338 + t346) + t586) * g(3) + t119 * t550 - t146 * t214 - t145 * t215 - t192 * t108 - t74 * t159 - t30 * t99 - t29 * t100 - g(1) * (t382 * t392 + t435) - g(2) * (t376 * t392 + t436) - t247 * t264 + t248 * t265 + t279 * t18 + t285 * t50 + t287 * t51 + t630 * (-t124 * t372 + t378 * t414 + t351 * t474 + (t373 * t473 + (t372 * t379 + t488) * qJD(4)) * pkin(3)); -t603 * (pkin(4) * t474 - t71) + (-Ifges(5,2) * t238 + t627) * t562 + (m(6) * t548 + t625) * t623 + t581 * (pkin(11) + t547) - t108 * t549 + t167 * t560 + t589 * pkin(4) * t473 + t633 + (-m(6) * t338 - m(7) * t454 + t590) * g(3) + t50 * t546 + t51 * t547 + t238 * t514 - t141 * t214 + t142 * t215 - t72 * t159 - t34 * t99 - t33 * t100 + (-t27 * t33 - t28 * t34 - t67 * t71 + t350 * t8 + (t372 * t67 + t378 * t420) * qJD(5) * pkin(4)) * m(7) + (-t184 * t549 + t69 * t71 - t70 * t72 + (t10 * t372 + t11 * t378 + (-t372 * t69 + t378 * t70) * qJD(5)) * pkin(4)) * m(6) - g(1) * (t382 * t390 + t435) - g(2) * (t376 * t390 + t436) + (-Ifges(5,6) * t238 + t606) * t556 + (-t524 + t607) * t561 - t260 * (mrSges(5,1) * t238 + t608) + t350 * t18; (t596 - t628) * g(3) - t8 * t576 - m(7) * (t27 * t37 + t28 * t38 + t67 * t70) + t603 * t70 + (m(7) * (-t537 + t538 + t595) + t585) * pkin(11) + (-t171 / 0.2e1 - t520 / 0.2e1 - t406 / 0.2e1 - t408 / 0.2e1 - t407 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t619 + t421 * mrSges(7,3) + t579) * t620 + (t382 * t389 - t435) * g(1) + (t518 / 0.2e1 - t516 / 0.2e1 - t517 / 0.2e1 - t519 / 0.2e1 + t523 / 0.2e1 + t578) * t619 + t391 * mrSges(7,3) - t69 * t159 - t38 * t99 - t37 * t100 - pkin(5) * t18 + t386 + (t376 * t389 - t436) * g(2); -t67 * (mrSges(7,1) * t158 + mrSges(7,2) * t157) + (Ifges(7,1) * t157 - t513) * t567 + t82 * t566 + (Ifges(7,5) * t157 - Ifges(7,6) * t158) * t565 - t27 * t99 + t28 * t100 - g(1) * (mrSges(7,1) * t272 - mrSges(7,2) * t273) - g(2) * (-mrSges(7,1) * t270 + mrSges(7,2) * t271) + g(3) * t428 * t340 + (t157 * t27 + t158 * t28) * mrSges(7,3) + t468 + (-Ifges(7,2) * t158 + t156 + t83) * t568 + t584;];
tau  = t1;
