% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:06:13
% EndTime: 2019-03-09 22:06:44
% DurationCPUTime: 17.44s
% Computational Cost: add. (24733->772), mult. (59246->1044), div. (0->0), fcn. (44107->10), ass. (0->375)
t419 = cos(qJ(4));
t416 = sin(qJ(3));
t420 = cos(qJ(3));
t421 = cos(qJ(2));
t466 = qJD(1) * t421;
t417 = sin(qJ(2));
t467 = qJD(1) * t417;
t367 = -t416 * t467 + t420 * t466;
t415 = sin(qJ(4));
t478 = t367 * t415;
t628 = -qJ(5) * t478 - t419 * qJD(5);
t384 = t416 * t421 + t420 * t417;
t368 = t384 * qJD(1);
t307 = pkin(3) * t368 - pkin(9) * t367;
t280 = pkin(2) * t467 + t307;
t553 = -pkin(8) - pkin(7);
t398 = t553 * t421;
t387 = qJD(1) * t398;
t369 = t416 * t387;
t397 = t553 * t417;
t386 = qJD(1) * t397;
t319 = t386 * t420 + t369;
t211 = t419 * t280 - t319 * t415;
t410 = t419 * qJ(5);
t443 = t368 * pkin(4) - t367 * t410;
t402 = pkin(2) * t416 + pkin(9);
t469 = -qJ(5) - t402;
t445 = qJD(4) * t469;
t501 = pkin(2) * qJD(3);
t460 = t420 * t501;
t627 = -t211 - t443 + (-qJD(5) - t460) * t415 + t419 * t445;
t373 = qJD(2) * pkin(2) + t386;
t311 = t373 * t420 + t369;
t215 = t419 * t307 - t311 * t415;
t508 = -qJ(5) - pkin(9);
t448 = qJD(4) * t508;
t626 = -qJD(5) * t415 + t419 * t448 - t215 - t443;
t212 = t415 * t280 + t419 * t319;
t625 = -t415 * t445 - t419 * t460 + t212 + t628;
t216 = t415 * t307 + t419 * t311;
t624 = -t415 * t448 + t216 + t628;
t412 = sin(pkin(11));
t413 = cos(pkin(11));
t379 = t412 * t419 + t413 * t415;
t267 = t379 * t367;
t363 = t379 * qJD(4);
t599 = t267 - t363;
t432 = t412 * t415 - t413 * t419;
t268 = t432 * t367;
t364 = t432 * qJD(4);
t598 = t268 - t364;
t411 = qJD(2) + qJD(3);
t284 = -pkin(3) * t411 - t311;
t331 = -t368 * t415 + t411 * t419;
t231 = -pkin(4) * t331 + qJD(5) + t284;
t332 = t368 * t419 + t411 * t415;
t446 = t413 * t331 - t332 * t412;
t139 = -pkin(5) * t446 + t231;
t244 = t331 * t412 + t332 * t413;
t414 = sin(qJ(6));
t418 = cos(qJ(6));
t144 = t244 * t418 + t414 * t446;
t361 = qJD(4) - t367;
t605 = pkin(10) * t244;
t405 = -pkin(2) * t421 - pkin(1);
t396 = qJD(1) * t405;
t271 = -t367 * pkin(3) - t368 * pkin(9) + t396;
t470 = t420 * t387;
t312 = t416 * t373 - t470;
t285 = pkin(9) * t411 + t312;
t190 = t419 * t271 - t285 * t415;
t153 = -qJ(5) * t332 + t190;
t138 = pkin(4) * t361 + t153;
t191 = t271 * t415 + t285 * t419;
t154 = qJ(5) * t331 + t191;
t147 = t412 * t154;
t84 = t413 * t138 - t147;
t55 = pkin(5) * t361 - t605 + t84;
t590 = pkin(10) * t446;
t473 = t413 * t154;
t85 = t412 * t138 + t473;
t58 = t85 + t590;
t20 = -t414 * t58 + t418 * t55;
t21 = t414 * t55 + t418 * t58;
t324 = t411 * t384;
t304 = t324 * qJD(1);
t383 = t416 * t417 - t420 * t421;
t323 = t411 * t383;
t303 = t323 * qJD(1);
t234 = qJD(4) * t331 - t303 * t419;
t235 = -qJD(4) * t332 + t303 * t415;
t129 = -t234 * t412 + t235 * t413;
t130 = t234 * t413 + t235 * t412;
t593 = -t244 * t414 + t418 * t446;
t39 = qJD(6) * t593 + t129 * t414 + t130 * t418;
t40 = -qJD(6) * t144 + t129 * t418 - t130 * t414;
t462 = Ifges(7,5) * t39 + Ifges(7,6) * t40 + Ifges(7,3) * t304;
t502 = Ifges(7,4) * t144;
t346 = qJD(6) + t361;
t529 = -t346 / 0.2e1;
t541 = -t144 / 0.2e1;
t465 = qJD(2) * t417;
t461 = pkin(2) * t465;
t594 = qJD(1) * t461;
t204 = pkin(3) * t304 + pkin(9) * t303 + t594;
t455 = qJD(2) * t553;
t444 = qJD(1) * t455;
t375 = t417 * t444;
t431 = t421 * t444;
t223 = qJD(3) * t311 + t420 * t375 + t416 * t431;
t89 = -qJD(4) * t191 + t419 * t204 - t223 * t415;
t43 = pkin(4) * t304 - qJ(5) * t234 - qJD(5) * t332 + t89;
t463 = qJD(4) * t419;
t464 = qJD(4) * t415;
t88 = t415 * t204 + t419 * t223 + t271 * t463 - t285 * t464;
t50 = qJ(5) * t235 + qJD(5) * t331 + t88;
t17 = t412 * t43 + t413 * t50;
t10 = pkin(10) * t129 + t17;
t16 = -t412 * t50 + t413 * t43;
t7 = pkin(5) * t304 - pkin(10) * t130 + t16;
t3 = qJD(6) * t20 + t10 * t418 + t414 * t7;
t4 = -qJD(6) * t21 - t10 * t414 + t418 * t7;
t589 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t623 = t462 + t589 + (Ifges(7,5) * t593 - Ifges(7,6) * t144) * t529 + (t144 * t21 + t20 * t593) * mrSges(7,3) - t139 * (mrSges(7,1) * t144 + mrSges(7,2) * t593) + (Ifges(7,1) * t593 - t502) * t541;
t522 = -t415 / 0.2e1;
t580 = t625 * t412 + t413 * t627;
t579 = t412 * t627 - t625 * t413;
t576 = t412 * t624 + t413 * t626;
t575 = t412 * t626 - t413 * t624;
t622 = t599 * pkin(10);
t621 = -t368 * pkin(5) - pkin(10) * t598;
t140 = Ifges(7,4) * t593;
t620 = -Ifges(7,2) * t144 + t140;
t618 = t411 * Ifges(4,6) / 0.2e1;
t617 = Ifges(5,3) + Ifges(6,3);
t612 = t621 + t580;
t611 = t622 + t579;
t610 = t621 + t576;
t609 = t622 + t575;
t505 = Ifges(5,4) * t332;
t220 = Ifges(5,2) * t331 + Ifges(5,6) * t361 + t505;
t441 = mrSges(5,1) * t415 + mrSges(5,2) * t419;
t608 = t220 * t522 + t284 * t441;
t489 = t411 * Ifges(4,5);
t607 = t396 * mrSges(4,2) + t489 / 0.2e1;
t490 = t368 * Ifges(4,4);
t606 = t618 + t490 / 0.2e1 + t367 * Ifges(4,2) / 0.2e1;
t604 = Ifges(6,4) * t244;
t602 = t190 * t419;
t181 = -t267 * t418 + t268 * t414;
t314 = t379 * t418 - t414 * t432;
t218 = -qJD(6) * t314 - t363 * t418 + t364 * t414;
t601 = t181 - t218;
t182 = -t267 * t414 - t268 * t418;
t313 = -t379 * t414 - t418 * t432;
t217 = qJD(6) * t313 - t363 * t414 - t364 * t418;
t600 = t182 - t217;
t407 = pkin(4) * t464;
t597 = -pkin(5) * t599 + t407;
t318 = t386 * t416 - t470;
t345 = pkin(4) * t478;
t596 = t416 * t501 - t318 - t345;
t595 = -t89 * t415 + t419 * t88;
t436 = Ifges(5,5) * t419 - Ifges(5,6) * t415;
t503 = Ifges(5,4) * t419;
t438 = -Ifges(5,2) * t415 + t503;
t504 = Ifges(5,4) * t415;
t440 = Ifges(5,1) * t419 - t504;
t328 = Ifges(5,4) * t331;
t221 = Ifges(5,1) * t332 + Ifges(5,5) * t361 + t328;
t471 = t419 * t221;
t526 = t361 / 0.2e1;
t530 = t332 / 0.2e1;
t592 = t436 * t526 + t440 * t530 + t331 * t438 / 0.2e1 + t471 / 0.2e1 + t608;
t591 = -t396 * mrSges(4,1) - t190 * mrSges(5,1) - t84 * mrSges(6,1) - t20 * mrSges(7,1) + t191 * mrSges(5,2) + t85 * mrSges(6,2) + t21 * mrSges(7,2) + t606;
t561 = t39 / 0.2e1;
t560 = t40 / 0.2e1;
t549 = t129 / 0.2e1;
t548 = t130 / 0.2e1;
t533 = t304 / 0.2e1;
t376 = t469 * t415;
t377 = t402 * t419 + t410;
t301 = t413 * t376 - t377 * t412;
t513 = pkin(10) * t379;
t257 = t301 - t513;
t302 = t412 * t376 + t413 * t377;
t374 = t432 * pkin(10);
t258 = -t374 + t302;
t165 = t257 * t414 + t258 * t418;
t588 = -qJD(6) * t165 - t414 * t611 + t418 * t612;
t164 = t257 * t418 - t258 * t414;
t587 = qJD(6) * t164 + t414 * t612 + t418 * t611;
t394 = t508 * t415;
t395 = pkin(9) * t419 + t410;
t329 = t413 * t394 - t395 * t412;
t275 = t329 - t513;
t330 = t412 * t394 + t413 * t395;
t276 = -t374 + t330;
t188 = t275 * t414 + t276 * t418;
t586 = -qJD(6) * t188 - t414 * t609 + t418 * t610;
t187 = t275 * t418 - t276 * t414;
t585 = qJD(6) * t187 + t414 * t610 + t418 * t609;
t584 = mrSges(4,2) * t368;
t583 = Ifges(6,4) * t446;
t400 = pkin(4) * t413 + pkin(5);
t516 = pkin(4) * t412;
t354 = t400 * t418 - t414 * t516;
t91 = -t153 * t412 - t473;
t64 = t91 - t590;
t92 = t413 * t153 - t147;
t65 = t92 - t605;
t582 = t354 * qJD(6) - t414 * t64 - t418 * t65;
t355 = t400 * t414 + t418 * t516;
t581 = -t355 * qJD(6) + t414 * t65 - t418 * t64;
t574 = Ifges(5,5) * t234 + Ifges(5,6) * t235;
t573 = t596 + t597;
t310 = t383 * pkin(3) - t384 * pkin(9) + t405;
t334 = t397 * t416 - t398 * t420;
t327 = t419 * t334;
t237 = t415 * t310 + t327;
t253 = t345 + t312;
t572 = -t253 + t597;
t571 = t420 * t397 + t398 * t416;
t570 = t407 + t596;
t569 = Ifges(6,5) * t130 + Ifges(6,6) * t129 + t304 * t617 + t574;
t568 = -t190 * t415 + t191 * t419;
t135 = -mrSges(5,1) * t235 + mrSges(5,2) * t234;
t224 = qJD(3) * t312 + t375 * t416 - t420 * t431;
t567 = m(5) * t224 + t135;
t566 = t89 * mrSges(5,1) + t16 * mrSges(6,1) - t88 * mrSges(5,2) - t17 * mrSges(6,2);
t565 = m(4) / 0.2e1;
t564 = Ifges(5,3) / 0.2e1;
t563 = Ifges(7,4) * t561 + Ifges(7,2) * t560 + Ifges(7,6) * t533;
t562 = Ifges(7,1) * t561 + Ifges(7,4) * t560 + Ifges(7,5) * t533;
t559 = Ifges(6,4) * t548 + Ifges(6,2) * t549 + Ifges(6,6) * t533;
t558 = Ifges(6,1) * t548 + Ifges(6,4) * t549 + Ifges(6,5) * t533;
t76 = Ifges(7,2) * t593 + Ifges(7,6) * t346 + t502;
t557 = -t76 / 0.2e1;
t556 = t76 / 0.2e1;
t77 = Ifges(7,1) * t144 + Ifges(7,5) * t346 + t140;
t555 = -t77 / 0.2e1;
t554 = t77 / 0.2e1;
t552 = pkin(1) * mrSges(3,1);
t551 = pkin(1) * mrSges(3,2);
t132 = Ifges(6,2) * t446 + Ifges(6,6) * t361 + t604;
t547 = -t132 / 0.2e1;
t546 = t132 / 0.2e1;
t133 = Ifges(6,1) * t244 + Ifges(6,5) * t361 + t583;
t545 = -t133 / 0.2e1;
t544 = t133 / 0.2e1;
t543 = -t593 / 0.2e1;
t542 = t593 / 0.2e1;
t540 = t144 / 0.2e1;
t539 = t234 / 0.2e1;
t538 = t235 / 0.2e1;
t537 = -t446 / 0.2e1;
t536 = t446 / 0.2e1;
t535 = -t244 / 0.2e1;
t534 = t244 / 0.2e1;
t532 = -t331 / 0.2e1;
t531 = -t332 / 0.2e1;
t528 = t346 / 0.2e1;
t527 = -t361 / 0.2e1;
t525 = -t367 / 0.2e1;
t521 = t419 / 0.2e1;
t520 = m(4) * t396;
t518 = pkin(2) * t420;
t517 = pkin(4) * t332;
t430 = qJ(5) * t323 - qJD(5) * t384;
t229 = pkin(3) * t324 + pkin(9) * t323 + t461;
t389 = t417 * t455;
t390 = t421 * t455;
t247 = qJD(3) * t571 + t389 * t420 + t390 * t416;
t447 = t419 * t229 - t247 * t415;
t63 = pkin(4) * t324 + t430 * t419 + (-t327 + (qJ(5) * t384 - t310) * t415) * qJD(4) + t447;
t454 = t384 * t463;
t457 = t415 * t229 + t419 * t247 + t310 * t463;
t81 = -qJ(5) * t454 + (-qJD(4) * t334 + t430) * t415 + t457;
t25 = t412 * t63 + t413 * t81;
t507 = mrSges(4,3) * t367;
t506 = Ifges(3,4) * t417;
t356 = Ifges(4,4) * t367;
t500 = t593 * Ifges(7,6);
t499 = t144 * Ifges(7,5);
t498 = t446 * Ifges(6,6);
t497 = t244 * Ifges(6,5);
t496 = t331 * Ifges(5,6);
t495 = t332 * Ifges(5,5);
t494 = t346 * Ifges(7,3);
t492 = t368 * mrSges(4,3);
t491 = t368 * Ifges(4,1);
t485 = Ifges(3,5) * qJD(2);
t484 = Ifges(3,6) * qJD(2);
t483 = qJD(2) * mrSges(3,1);
t482 = qJD(2) * mrSges(3,2);
t479 = t224 * t571;
t475 = t384 * t415;
t236 = t419 * t310 - t334 * t415;
t178 = pkin(4) * t383 - t384 * t410 + t236;
t201 = -qJ(5) * t475 + t237;
t106 = t412 * t178 + t413 * t201;
t468 = mrSges(4,1) * t411 + mrSges(5,1) * t331 - mrSges(5,2) * t332 - t492;
t404 = -pkin(4) * t419 - pkin(3);
t453 = t485 / 0.2e1;
t452 = -t484 / 0.2e1;
t13 = -t40 * mrSges(7,1) + t39 * mrSges(7,2);
t24 = -t412 * t81 + t413 * t63;
t67 = -t129 * mrSges(6,1) + t130 * mrSges(6,2);
t105 = t413 * t178 - t201 * t412;
t272 = pkin(4) * t475 - t571;
t442 = mrSges(5,1) * t419 - mrSges(5,2) * t415;
t439 = Ifges(5,1) * t415 + t503;
t437 = Ifges(5,2) * t419 + t504;
t435 = Ifges(5,5) * t415 + Ifges(5,6) * t419;
t288 = t432 * t384;
t90 = pkin(5) * t383 + pkin(10) * t288 + t105;
t287 = t379 * t384;
t94 = -pkin(10) * t287 + t106;
t32 = -t414 * t94 + t418 * t90;
t33 = t414 * t90 + t418 * t94;
t167 = mrSges(5,1) * t304 - mrSges(5,3) * t234;
t168 = -mrSges(5,2) * t304 + mrSges(5,3) * t235;
t434 = -t415 * t167 + t419 * t168;
t433 = -t191 * t415 - t602;
t199 = -t287 * t418 + t288 * t414;
t200 = -t287 * t414 - t288 * t418;
t339 = pkin(5) * t432 + t404;
t427 = t433 * mrSges(5,3);
t248 = qJD(3) * t334 + t389 * t416 - t420 * t390;
t163 = t248 + (-t323 * t415 + t454) * pkin(4);
t134 = -pkin(4) * t235 + t224;
t423 = m(5) * (qJD(4) * t433 + t595);
t113 = t234 * Ifges(5,4) + t235 * Ifges(5,2) + t304 * Ifges(5,6);
t114 = Ifges(5,1) * t234 + Ifges(5,4) * t235 + Ifges(5,5) * t304;
t131 = t361 * Ifges(6,3) + t497 + t498;
t219 = t361 * Ifges(5,3) + t495 + t496;
t282 = t356 + t489 + t491;
t70 = -pkin(5) * t129 + t134;
t75 = t494 + t499 + t500;
t422 = t311 * t507 + t113 * t521 + (Ifges(7,5) * t182 + Ifges(7,6) * t181) * t529 + (-Ifges(6,5) * t268 - Ifges(6,6) * t267) * t527 + (Ifges(7,5) * t217 + Ifges(7,6) * t218) * t528 + (t356 + t471 + t282) * t525 + (-Ifges(6,5) * t364 - Ifges(6,6) * t363) * t526 + (-Ifges(6,1) * t364 - Ifges(6,4) * t363) * t534 + (-Ifges(6,4) * t364 - Ifges(6,2) * t363) * t536 + t592 * qJD(4) + t182 * t555 + t218 * t556 + t181 * t557 + t379 * t558 + (Ifges(7,4) * t314 + Ifges(7,2) * t313) * t560 + (Ifges(7,1) * t314 + Ifges(7,4) * t313) * t561 + t314 * t562 + t313 * t563 - (Ifges(4,1) * t367 + t131 + t219 - t490 + t75) * t368 / 0.2e1 + (Ifges(6,1) * t379 - Ifges(6,4) * t432) * t548 + (Ifges(6,4) * t379 - Ifges(6,2) * t432) * t549 + t134 * (mrSges(6,1) * t432 + mrSges(6,2) * t379) - t432 * t559 + (Ifges(6,5) * t379 + Ifges(7,5) * t314 - Ifges(6,6) * t432 + Ifges(7,6) * t313 + t435) * t533 + (-Ifges(6,4) * t268 - Ifges(6,2) * t267) * t537 + (-t16 * t379 - t17 * t432 - t598 * t84 + t599 * t85) * mrSges(6,3) + (-mrSges(6,1) * t599 + mrSges(6,2) * t598) * t231 + (mrSges(7,1) * t601 - mrSges(7,2) * t600) * t139 + (t20 * t600 - t21 * t601 + t3 * t313 - t314 * t4) * mrSges(7,3) + (t191 * t478 + t595) * mrSges(5,3) + t437 * t538 + t439 * t539 + (Ifges(7,1) * t217 + Ifges(7,4) * t218) * t540 + (Ifges(7,4) * t217 + Ifges(7,2) * t218) * t542 - t364 * t544 - t268 * t545 - t363 * t546 - t267 * t547 + t217 * t554 + (-t442 - mrSges(4,1)) * t224 + (Ifges(7,1) * t182 + Ifges(7,4) * t181) * t541 + (-Ifges(6,1) * t268 - Ifges(6,4) * t267) * t535 + t415 * t114 / 0.2e1 + t70 * (-mrSges(7,1) * t313 + mrSges(7,2) * t314) - Ifges(4,5) * t303 - Ifges(4,6) * t304 + (mrSges(5,3) * t602 + t436 * t527 + t438 * t532 + t440 * t531 - t607 - t608) * t367 + (Ifges(7,4) * t182 + Ifges(7,2) * t181) * t543 + (Ifges(5,5) * t531 + Ifges(6,5) * t535 + Ifges(7,5) * t541 - Ifges(4,2) * t525 + Ifges(5,6) * t532 + Ifges(6,6) * t537 + Ifges(7,6) * t543 + Ifges(7,3) * t529 + t617 * t527 + t591 + t618) * t368 - t223 * mrSges(4,2);
t406 = Ifges(3,4) * t466;
t393 = mrSges(3,3) * t466 - t482;
t392 = -mrSges(3,3) * t467 + t483;
t391 = t404 - t518;
t366 = Ifges(3,1) * t467 + t406 + t485;
t365 = t484 + (Ifges(3,2) * t421 + t506) * qJD(1);
t337 = -mrSges(4,2) * t411 + t507;
t335 = t339 - t518;
t306 = -mrSges(4,1) * t367 + t584;
t263 = mrSges(5,1) * t361 - mrSges(5,3) * t332;
t262 = -mrSges(5,2) * t361 + mrSges(5,3) * t331;
t206 = mrSges(6,1) * t361 - mrSges(6,3) * t244;
t205 = -mrSges(6,2) * t361 + mrSges(6,3) * t446;
t203 = pkin(5) * t287 + t272;
t189 = pkin(5) * t244 + t517;
t161 = t323 * t432 - t363 * t384;
t160 = t323 * t379 + t364 * t384;
t146 = -mrSges(6,1) * t446 + mrSges(6,2) * t244;
t125 = mrSges(7,1) * t346 - mrSges(7,3) * t144;
t124 = -mrSges(7,2) * t346 + mrSges(7,3) * t593;
t110 = mrSges(6,1) * t304 - mrSges(6,3) * t130;
t109 = -mrSges(6,2) * t304 + mrSges(6,3) * t129;
t103 = -qJD(4) * t237 + t447;
t102 = -t334 * t464 + t457;
t95 = -pkin(5) * t160 + t163;
t86 = -mrSges(7,1) * t593 + mrSges(7,2) * t144;
t60 = -qJD(6) * t200 + t160 * t418 - t161 * t414;
t59 = qJD(6) * t199 + t160 * t414 + t161 * t418;
t35 = -mrSges(7,2) * t304 + mrSges(7,3) * t40;
t34 = mrSges(7,1) * t304 - mrSges(7,3) * t39;
t19 = pkin(10) * t160 + t25;
t18 = pkin(5) * t324 - pkin(10) * t161 + t24;
t6 = -qJD(6) * t33 + t18 * t418 - t19 * t414;
t5 = qJD(6) * t32 + t18 * t414 + t19 * t418;
t1 = [(Ifges(6,5) * t161 + Ifges(6,6) * t160) * t526 + (Ifges(7,5) * t59 + Ifges(7,6) * t60) * t528 + (Ifges(6,1) * t161 + Ifges(6,4) * t160) * t534 + (t462 + t569 + t574) * t383 / 0.2e1 + m(7) * (t139 * t95 + t20 * t6 + t203 * t70 + t21 * t5 + t3 * t33 + t32 * t4) + m(6) * (t105 * t16 + t106 * t17 + t134 * t272 + t163 * t231 + t24 * t84 + t25 * t85) + (t405 * mrSges(4,1) + (t564 + Ifges(4,2)) * t383 - Ifges(4,4) * t384 - t334 * mrSges(4,3)) * t304 + (mrSges(4,1) * t594 - mrSges(4,3) * t223 + Ifges(4,4) * t303 + Ifges(6,5) * t548 + Ifges(7,5) * t561 + Ifges(6,6) * t549 + Ifges(7,6) * t560 + (Ifges(6,3) + Ifges(7,3)) * t533 + t566 + t589) * t383 + (t440 * t539 + t438 * t538 - Ifges(4,1) * t303 + t113 * t522 + t114 * t521 + (mrSges(4,3) + t441) * t224 + (-t415 * t88 - t419 * t89) * mrSges(5,3) + (-t419 * t220 / 0.2e1 + t221 * t522 + t284 * t442 + t435 * t527 + t437 * t532 + t439 * t531 - t568 * mrSges(5,3)) * qJD(4)) * t384 + t59 * t554 + t60 * t556 - t288 * t558 - t287 * t559 + t200 * t562 + t199 * t563 + (t75 / 0.2e1 - t591 + t219 / 0.2e1 + t131 / 0.2e1 + t497 / 0.2e1 + t498 / 0.2e1 + t494 / 0.2e1 + t495 / 0.2e1 + t496 / 0.2e1 + (t564 + Ifges(6,3) / 0.2e1) * t361 + t500 / 0.2e1 + t499 / 0.2e1 - t606) * t324 + (t311 * t323 - t312 * t324) * mrSges(4,3) + m(4) * (t223 * t334 + t247 * t312 - t248 * t311 - t479) + m(5) * (t102 * t191 + t103 * t190 + t236 * t89 + t237 * t88 + t248 * t284 - t479) + (t199 * t3 - t20 * t59 - t200 * t4 + t21 * t60) * mrSges(7,3) + (Ifges(6,4) * t161 + Ifges(6,2) * t160) * t536 + (Ifges(7,1) * t59 + Ifges(7,4) * t60) * t540 + (Ifges(7,4) * t59 + Ifges(7,2) * t60) * t542 + t161 * t544 + t160 * t546 + (t16 * t288 + t160 * t85 - t161 * t84 - t17 * t287) * mrSges(6,3) + (-Ifges(6,5) * t288 + Ifges(7,5) * t200 - Ifges(6,6) * t287 + Ifges(7,6) * t199 + t384 * t436) * t533 + (-Ifges(6,4) * t288 - Ifges(6,2) * t287) * t549 + t134 * (mrSges(6,1) * t287 - mrSges(6,2) * t288) + (-Ifges(6,1) * t288 - Ifges(6,4) * t287) * t548 + (Ifges(7,4) * t200 + Ifges(7,2) * t199) * t560 + t247 * t337 - (mrSges(4,2) * t405 - mrSges(4,3) * t571) * t303 - t571 * t135 + (-pkin(7) * t392 + t366 / 0.2e1 + t453 + (-0.2e1 * t551 + 0.3e1 / 0.2e1 * Ifges(3,4) * t421) * qJD(1)) * t421 * qJD(2) - (t356 / 0.2e1 + t491 / 0.2e1 + t282 / 0.2e1 + t427 + t592 + t607) * t323 + (-pkin(7) * t393 - t365 / 0.2e1 + t452 + (-0.2e1 * t552 - 0.3e1 / 0.2e1 * t506 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t421) * qJD(1) + (t306 + 0.2e1 * t520 + t584) * pkin(2)) * t465 + t32 * t34 + t33 * t35 + t95 * t86 + t106 * t109 + t105 * t110 + t5 * t124 + t6 * t125 + t139 * (-mrSges(7,1) * t60 + mrSges(7,2) * t59) - t468 * t248 + t163 * t146 + t70 * (-mrSges(7,1) * t199 + mrSges(7,2) * t200) + t203 * t13 + t25 * t205 + t24 * t206 + (Ifges(7,1) * t200 + Ifges(7,4) * t199) * t561 + t231 * (-mrSges(6,1) * t160 + mrSges(6,2) * t161) + t236 * t167 + t237 * t168 + t102 * t262 + t103 * t263 + t272 * t67; t570 * t146 + t573 * t86 + t422 + t427 * qJD(4) + (t423 + (-t262 * t415 - t263 * t419) * qJD(4) + t434) * t402 + ((t303 * t420 - t304 * t416) * mrSges(4,3) + (-t468 * t416 + (t262 * t419 - t263 * t415 + t337) * t420) * qJD(3)) * pkin(2) + t579 * t205 + (t134 * t391 + t16 * t301 + t17 * t302 + t231 * t570 + t579 * t85 + t580 * t84) * m(6) + t580 * t206 + t567 * (-pkin(3) - t518) + 0.2e1 * ((t223 * t416 - t224 * t420) * t565 + ((-t311 * t416 + t312 * t420) * t565 + m(5) * (t284 * t416 + t420 * t568) / 0.2e1) * qJD(3)) * pkin(2) - m(5) * (t190 * t211 + t191 * t212 + t284 * t318) + t312 * t492 - m(4) * (-t311 * t318 + t312 * t319) + t391 * t67 + t335 * t13 - t319 * t337 + t301 * t110 + t302 * t109 + ((t453 - t366 / 0.2e1 - t406 / 0.2e1 + qJD(1) * t551 + (t392 - t483) * pkin(7)) * t421 + (t452 + t365 / 0.2e1 + (t552 + t506 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t421) * qJD(1) + (t393 + t482) * pkin(7) + (-t306 - t520) * pkin(2)) * t417) * qJD(1) + t587 * t124 + (t139 * t573 + t164 * t4 + t165 * t3 + t20 * t588 + t21 * t587 + t335 * t70) * m(7) + t588 * t125 + t468 * t318 + t164 * t34 + t165 * t35 - t212 * t262 - t211 * t263; -m(5) * (t190 * t215 + t191 * t216 + t284 * t312) + t572 * t86 + t576 * t206 + t575 * t205 + t422 + t586 * t125 + pkin(9) * t423 + t585 * t124 + t404 * t67 + ((-mrSges(5,3) * t190 - pkin(9) * t263) * t419 + (-t191 * mrSges(5,3) + pkin(4) * t146 - pkin(9) * t262) * t415) * qJD(4) + t339 * t13 - t311 * t337 + t329 * t110 + t330 * t109 + t434 * pkin(9) + t187 * t34 + t188 * t35 + (t468 + t492) * t312 - t253 * t146 - t216 * t262 - t215 * t263 - t567 * pkin(3) + (t139 * t572 + t187 * t4 + t188 * t3 + t20 * t586 + t21 * t585 + t339 * t70) * m(7) + (t134 * t404 + t16 * t329 + t17 * t330 + t575 * t85 + t576 * t84 + (-t253 + t407) * t231) * m(6); t220 * t530 + (Ifges(5,1) * t331 - t505) * t531 + t620 * t543 + t446 * t545 + (-Ifges(5,2) * t332 + t221 + t328) * t532 + t581 * t125 + t582 * t124 + (-t139 * t189 + t20 * t581 + t21 * t582 + t3 * t355 + t354 * t4) * m(7) + ((t16 * t413 + t17 * t412) * pkin(4) - t231 * t517 - t84 * t91 - t85 * t92) * m(6) + (Ifges(5,5) * t331 + Ifges(6,5) * t446 - Ifges(5,6) * t332 - Ifges(6,6) * t244) * t527 + (t244 * t85 + t446 * t84) * mrSges(6,3) - t231 * (mrSges(6,1) * t244 + mrSges(6,2) * t446) - t244 * t547 + (-Ifges(6,2) * t244 + t583) * t537 + (t109 * t412 + t110 * t413 - t146 * t332) * pkin(4) + t623 + (t190 * t331 + t191 * t332) * mrSges(5,3) + t354 * t34 + t355 * t35 - t284 * (mrSges(5,1) * t332 + mrSges(5,2) * t331) + t593 * t555 + t569 + t566 + (Ifges(6,1) * t446 - t604) * t535 - t189 * t86 - t144 * t557 - t92 * t205 - t91 * t206 - t190 * t262 + t191 * t263; -t593 * t124 + t144 * t125 - t446 * t205 + t244 * t206 + t13 + t67 + (t144 * t20 - t21 * t593 + t70) * m(7) + (t244 * t84 - t446 * t85 + t134) * m(6); t76 * t540 - t20 * t124 + t21 * t125 + (t620 + t77) * t543 + t623;];
tauc  = t1(:);
