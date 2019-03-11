% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:33
% EndTime: 2019-03-09 12:29:01
% DurationCPUTime: 55.78s
% Computational Cost: add. (21912->924), mult. (53473->1218), div. (0->0), fcn. (43873->14), ass. (0->406)
t647 = Ifges(6,4) + Ifges(7,4);
t366 = sin(pkin(11));
t368 = cos(pkin(11));
t373 = sin(qJ(4));
t377 = cos(qJ(4));
t322 = t366 * t377 + t368 * t373;
t367 = sin(pkin(6));
t378 = cos(qJ(2));
t492 = t367 * t378;
t383 = t322 * t492;
t259 = qJD(1) * t383;
t312 = t322 * qJD(4);
t666 = t259 - t312;
t321 = t366 * t373 - t377 * t368;
t382 = t321 * t492;
t260 = qJD(1) * t382;
t311 = t321 * qJD(4);
t679 = -t260 + t311;
t648 = Ifges(6,1) + Ifges(7,1);
t646 = Ifges(6,5) + Ifges(7,5);
t645 = Ifges(6,2) + Ifges(7,2);
t644 = Ifges(6,6) + Ifges(7,6);
t374 = sin(qJ(2));
t406 = pkin(2) * t374 - qJ(3) * t378;
t480 = qJD(1) * t367;
t300 = t406 * t480;
t455 = t374 * t480;
t369 = cos(pkin(6));
t479 = qJD(1) * t369;
t469 = pkin(1) * t479;
t301 = -pkin(8) * t455 + t378 * t469;
t215 = t368 * t300 - t366 * t301;
t490 = t368 * t378;
t385 = (pkin(3) * t374 - pkin(9) * t490) * t367;
t174 = qJD(1) * t385 + t215;
t216 = t366 * t300 + t368 * t301;
t454 = t378 * t480;
t429 = t366 * t454;
t194 = -pkin(9) * t429 + t216;
t539 = pkin(9) + qJ(3);
t328 = t539 * t366;
t329 = t539 * t368;
t618 = -t377 * t328 - t329 * t373;
t623 = -t321 * qJD(3) + qJD(4) * t618 - t373 * t174 - t377 * t194;
t335 = qJD(4) - t454;
t352 = qJD(2) + t479;
t279 = t352 * t368 - t366 * t455;
t280 = t352 * t366 + t368 * t455;
t434 = t377 * t279 - t280 * t373;
t401 = t279 * t373 + t377 * t280;
t528 = Ifges(5,4) * t401;
t128 = Ifges(5,2) * t434 + t335 * Ifges(5,6) + t528;
t556 = t335 / 0.2e1;
t569 = t401 / 0.2e1;
t571 = t434 / 0.2e1;
t201 = qJD(5) - t434;
t573 = t201 / 0.2e1;
t372 = sin(qJ(5));
t376 = cos(qJ(5));
t168 = t335 * t372 + t376 * t401;
t577 = t168 / 0.2e1;
t167 = t335 * t376 - t372 * t401;
t579 = t167 / 0.2e1;
t643 = Ifges(6,3) + Ifges(7,3);
t254 = -t352 * pkin(2) + qJD(3) - t301;
t202 = -t279 * pkin(3) + t254;
t346 = pkin(8) * t454;
t302 = t374 * t469 + t346;
t263 = qJ(3) * t352 + t302;
t288 = (-pkin(2) * t378 - qJ(3) * t374 - pkin(1)) * t367;
t268 = qJD(1) * t288;
t176 = -t366 * t263 + t368 * t268;
t137 = -pkin(3) * t454 - t280 * pkin(9) + t176;
t177 = t368 * t263 + t366 * t268;
t148 = pkin(9) * t279 + t177;
t78 = t137 * t373 + t148 * t377;
t75 = pkin(10) * t335 + t78;
t89 = -pkin(4) * t434 - pkin(10) * t401 + t202;
t29 = -t372 * t75 + t376 * t89;
t25 = -qJ(6) * t168 + t29;
t21 = pkin(5) * t201 + t25;
t30 = t372 * t89 + t376 * t75;
t26 = qJ(6) * t167 + t30;
t663 = t202 * mrSges(5,1) + t29 * mrSges(6,1) + t21 * mrSges(7,1) - t30 * mrSges(6,2) - t26 * mrSges(7,2);
t678 = t663 - Ifges(5,4) * t569 - Ifges(5,2) * t571 - Ifges(5,6) * t556 + t573 * t643 + t577 * t646 + t579 * t644 - t128 / 0.2e1;
t677 = pkin(10) * t455 - t623;
t252 = pkin(3) * t429 + t302;
t676 = -pkin(4) * t666 + pkin(10) * t679 - t252;
t675 = t647 * t167;
t674 = t647 * t168;
t476 = qJD(2) * t378;
t306 = (qJD(1) * t476 + qJDD(1) * t374) * t367;
t470 = qJDD(1) * t369;
t351 = qJDD(2) + t470;
t240 = -t306 * t366 + t351 * t368;
t241 = t306 * t368 + t351 * t366;
t114 = -qJD(4) * t401 + t240 * t377 - t241 * t373;
t112 = qJDD(5) - t114;
t473 = qJD(4) * t377;
t474 = qJD(4) * t373;
t550 = pkin(1) * t369;
t468 = qJD(2) * t550;
t431 = qJD(1) * t468;
t464 = pkin(1) * t470;
t477 = qJD(2) * t367;
t453 = t374 * t477;
t656 = -qJD(1) * t453 + qJDD(1) * t492;
t230 = pkin(8) * t656 + t374 * t464 + t378 * t431;
t185 = qJ(3) * t351 + qJD(3) * t352 + t230;
t475 = qJD(3) * t374;
t197 = -pkin(2) * t656 - qJ(3) * t306 + (-pkin(1) * qJDD(1) - qJD(1) * t475) * t367;
t125 = -t185 * t366 + t368 * t197;
t83 = -pkin(3) * t656 - pkin(9) * t241 + t125;
t126 = t368 * t185 + t366 * t197;
t93 = pkin(9) * t240 + t126;
t19 = t137 * t473 - t148 * t474 + t373 * t83 + t377 * t93;
t293 = qJDD(4) - t656;
t17 = pkin(10) * t293 + t19;
t113 = qJD(4) * t434 + t240 * t373 + t241 * t377;
t494 = t367 * t374;
t353 = pkin(8) * t494;
t231 = -qJD(2) * t346 - qJDD(1) * t353 - t374 * t431 + t378 * t464;
t208 = -t351 * pkin(2) + qJDD(3) - t231;
t147 = -t240 * pkin(3) + t208;
t38 = -t114 * pkin(4) - t113 * pkin(10) + t147;
t4 = -qJD(5) * t30 - t17 * t372 + t376 * t38;
t64 = qJD(5) * t167 + t113 * t376 + t293 * t372;
t1 = pkin(5) * t112 - qJ(6) * t64 - qJD(6) * t168 + t4;
t562 = t293 / 0.2e1;
t584 = t114 / 0.2e1;
t585 = t113 / 0.2e1;
t586 = t112 / 0.2e1;
t65 = -qJD(5) * t168 - t113 * t372 + t293 * t376;
t589 = t65 / 0.2e1;
t590 = t64 / 0.2e1;
t471 = qJD(5) * t376;
t472 = qJD(5) * t372;
t3 = t376 * t17 + t372 * t38 + t89 * t471 - t472 * t75;
t2 = qJ(6) * t65 + qJD(6) * t167 + t3;
t598 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t642 = t112 * t643 + t64 * t646 + t644 * t65;
t673 = t598 + t1 * mrSges(7,1) + t147 * mrSges(5,1) + t586 * t643 + t589 * t644 + t590 * t646 + t642 / 0.2e1 - t19 * mrSges(5,3) + (-t293 / 0.2e1 - t562) * Ifges(5,6) + (-t114 / 0.2e1 - t584) * Ifges(5,2) + (-t113 / 0.2e1 - t585) * Ifges(5,4);
t200 = Ifges(5,4) * t434;
t129 = Ifges(5,1) * t401 + t335 * Ifges(5,5) + t200;
t632 = t202 * mrSges(5,2);
t672 = t632 + Ifges(5,1) * t569 + Ifges(5,4) * t571 + Ifges(5,5) * t556 + t129 / 0.2e1;
t671 = t147 * mrSges(5,2) + 0.2e1 * Ifges(5,1) * t585 + 0.2e1 * Ifges(5,4) * t584 + 0.2e1 * Ifges(5,5) * t562;
t361 = pkin(5) * t376 + pkin(4);
t670 = -m(6) * pkin(4) - m(7) * t361 - mrSges(5,1);
t370 = -qJ(6) - pkin(10);
t386 = -m(6) * pkin(10) + m(7) * t370 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t641 = t112 * t644 + t64 * t647 + t645 * t65;
t668 = t641 / 0.2e1;
t634 = t167 * t645 + t201 * t644 + t674;
t633 = t168 * t648 + t646 * t201 + t675;
t258 = -t328 * t373 + t329 * t377;
t622 = -qJD(3) * t322 - qJD(4) * t258 - t174 * t377 + t373 * t194;
t221 = t260 * t372 + t376 * t455;
t452 = t322 * t471;
t506 = t311 * t372;
t393 = t452 - t506;
t620 = t221 + t393;
t509 = t434 * t372;
t665 = t472 - t509;
t593 = m(7) * pkin(5);
t640 = t646 * t112 + t64 * t648 + t647 * t65;
t660 = t640 / 0.2e1;
t635 = t167 * t644 + t168 * t646 + t201 * t643;
t659 = t230 * mrSges(3,2);
t658 = t372 * t677 + t376 * t676;
t360 = pkin(3) * t368 + pkin(2);
t234 = pkin(4) * t321 - pkin(10) * t322 - t360;
t657 = t234 * t471 + t372 * t676 - t376 * t677;
t621 = pkin(4) * t455 - t622;
t655 = t647 * t376;
t654 = t647 * t372;
t20 = -t137 * t474 - t148 * t473 - t373 * t93 + t377 * t83;
t653 = -t20 * mrSges(5,1) + t19 * mrSges(5,2);
t365 = pkin(11) + qJ(4);
t362 = sin(t365);
t363 = cos(t365);
t421 = -mrSges(4,1) * t368 + mrSges(4,2) * t366;
t388 = m(4) * pkin(2) - t421;
t652 = t362 * t386 + t363 * t670 - t388;
t551 = cos(qJ(1));
t457 = t551 * t374;
t375 = sin(qJ(1));
t487 = t375 * t378;
t314 = t369 * t457 + t487;
t458 = t367 * t551;
t243 = t314 * t363 - t362 * t458;
t456 = t551 * t378;
t488 = t374 * t375;
t313 = -t369 * t456 + t488;
t651 = t243 * t372 - t313 * t376;
t650 = -t243 * t376 - t313 * t372;
t541 = -mrSges(7,1) - mrSges(6,1);
t649 = mrSges(6,2) + mrSges(7,2);
t222 = -t260 * t376 + t372 * t455;
t248 = t376 * t258;
t398 = qJ(6) * t311 - qJD(6) * t322;
t639 = qJ(6) * t222 + t398 * t376 + (-t248 + (qJ(6) * t322 - t234) * t372) * qJD(5) + t658 - t666 * pkin(5);
t638 = (-qJD(5) * t258 + t398) * t372 + t657 + (-t221 - t452) * qJ(6);
t158 = t372 * t234 + t248;
t637 = -qJD(5) * t158 + t658;
t636 = -t258 * t472 + t657;
t631 = -m(4) * qJ(3) - mrSges(4,3) - mrSges(5,3);
t436 = qJD(5) * t370;
t133 = pkin(4) * t401 - pkin(10) * t434;
t77 = t137 * t377 - t373 * t148;
t47 = t372 * t133 + t376 * t77;
t630 = qJ(6) * t509 + qJD(6) * t376 + t372 * t436 - t47;
t46 = t376 * t133 - t372 * t77;
t508 = t434 * t376;
t629 = -pkin(5) * t401 + qJ(6) * t508 - qJD(6) * t372 + t376 * t436 - t46;
t628 = pkin(5) * t665 - t78;
t627 = t593 + mrSges(7,1);
t626 = pkin(5) * t620 + t621;
t309 = -t366 * t494 + t368 * t369;
t310 = t366 * t369 + t368 * t494;
t224 = t309 * t373 + t310 * t377;
t549 = pkin(1) * t378;
t290 = t353 + (-pkin(2) - t549) * t369;
t232 = -t309 * pkin(3) + t290;
t400 = t377 * t309 - t310 * t373;
t124 = -pkin(4) * t400 - t224 * pkin(10) + t232;
t319 = pkin(8) * t492 + t374 * t550;
t287 = qJ(3) * t369 + t319;
t210 = -t366 * t287 + t368 * t288;
t154 = -pkin(3) * t492 - t310 * pkin(9) + t210;
t211 = t368 * t287 + t366 * t288;
t171 = pkin(9) * t309 + t211;
t99 = t373 * t154 + t377 * t171;
t95 = -pkin(10) * t492 + t99;
t49 = t372 * t124 + t376 * t95;
t101 = -mrSges(6,1) * t167 + mrSges(6,2) * t168;
t537 = mrSges(5,3) * t401;
t173 = mrSges(5,1) * t335 - t537;
t625 = t173 - t101;
t433 = mrSges(3,3) * t455;
t624 = -mrSges(3,1) * t352 - mrSges(4,1) * t279 + mrSges(4,2) * t280 + t433;
t505 = t311 * t376;
t392 = t322 * t472 + t505;
t619 = t222 + t392;
t415 = mrSges(7,1) * t372 + mrSges(7,2) * t376;
t417 = mrSges(6,1) * t372 + mrSges(6,2) * t376;
t74 = -pkin(4) * t335 - t77;
t53 = -pkin(5) * t167 + qJD(6) + t74;
t616 = t53 * t415 + t74 * t417;
t615 = -t372 * t644 + t376 * t646;
t614 = -t372 * t645 + t655;
t613 = t376 * t648 - t654;
t612 = -t471 + t508;
t610 = -t125 * t366 + t126 * t368;
t609 = t3 * t376 - t372 * t4;
t608 = m(5) + m(6) + m(7);
t532 = Ifges(3,4) * t374;
t607 = pkin(1) * (mrSges(3,1) * t374 + mrSges(3,2) * t378) - t374 * (Ifges(3,1) * t378 - t532) / 0.2e1;
t606 = -mrSges(6,1) - t627;
t605 = mrSges(3,2) + t631;
t416 = -mrSges(7,1) * t376 + mrSges(7,2) * t372;
t418 = -mrSges(6,1) * t376 + mrSges(6,2) * t372;
t604 = -t416 - t418 - t670;
t601 = -t372 * t593 + t631;
t600 = mrSges(3,1) - t652;
t522 = Ifges(3,6) * t352;
t599 = t78 * mrSges(5,2) + t522 / 0.2e1 + (t378 * Ifges(3,2) + t532) * t480 / 0.2e1 - t77 * mrSges(5,1);
t557 = -t335 / 0.2e1;
t570 = -t401 / 0.2e1;
t597 = Ifges(5,1) * t570 + Ifges(5,5) * t557 - t632;
t572 = -t434 / 0.2e1;
t574 = -t201 / 0.2e1;
t578 = -t168 / 0.2e1;
t580 = -t167 / 0.2e1;
t595 = -Ifges(5,2) * t572 - Ifges(5,6) * t557 + t574 * t643 + t578 * t646 + t580 * t644 - t663;
t594 = t367 ^ 2;
t565 = t240 / 0.2e1;
t564 = t241 / 0.2e1;
t561 = t309 / 0.2e1;
t560 = t310 / 0.2e1;
t555 = t369 / 0.2e1;
t552 = t378 / 0.2e1;
t548 = pkin(5) * t168;
t538 = mrSges(5,3) * t434;
t536 = mrSges(6,3) * t167;
t535 = mrSges(6,3) * t168;
t534 = mrSges(7,3) * t167;
t533 = mrSges(7,3) * t168;
t531 = Ifges(3,4) * t378;
t530 = Ifges(4,4) * t366;
t529 = Ifges(4,4) * t368;
t523 = Ifges(4,5) * t280;
t521 = Ifges(4,6) * t279;
t520 = Ifges(4,3) * t374;
t517 = t434 * Ifges(5,6);
t516 = t401 * Ifges(5,5);
t515 = t335 * Ifges(5,3);
t514 = t352 * Ifges(3,5);
t502 = t314 * t372;
t316 = -t369 * t488 + t456;
t501 = t316 * t372;
t500 = t322 * t372;
t499 = t322 * t376;
t497 = t363 * t372;
t496 = t363 * t376;
t495 = t366 * t378;
t493 = t367 * t375;
t491 = t368 * (Ifges(4,1) * t280 + Ifges(4,4) * t279 - Ifges(4,5) * t454);
t489 = t372 * t378;
t486 = t376 * t378;
t265 = (qJD(2) * t406 - t475) * t367;
t303 = -pkin(8) * t453 + t378 * t468;
t274 = qJD(3) * t369 + t303;
t187 = t366 * t265 + t368 * t274;
t481 = t551 * pkin(1) + pkin(8) * t493;
t463 = t367 * t489;
t462 = t367 * t486;
t461 = Ifges(5,5) * t113 + Ifges(5,6) * t114 + Ifges(5,3) * t293;
t459 = Ifges(3,5) * t306 + Ifges(3,6) * t656 + Ifges(3,3) * t351;
t22 = -t65 * mrSges(7,1) + t64 * mrSges(7,2);
t439 = -t472 / 0.2e1;
t437 = -pkin(1) * t375 + pkin(8) * t458;
t160 = -t240 * mrSges(4,1) + t241 * mrSges(4,2);
t54 = -t114 * mrSges(5,1) + t113 * mrSges(5,2);
t48 = t376 * t124 - t372 * t95;
t98 = t154 * t377 - t373 * t171;
t157 = t376 * t234 - t258 * t372;
t186 = t368 * t265 - t366 * t274;
t432 = mrSges(3,3) * t454;
t430 = t366 * t458;
t428 = t366 * t367 * t476;
t94 = pkin(4) * t492 - t98;
t420 = mrSges(4,1) * t366 + mrSges(4,2) * t368;
t414 = Ifges(4,1) * t368 - t530;
t411 = -Ifges(4,2) * t366 + t529;
t315 = t369 * t487 + t457;
t403 = t366 * pkin(3) * t493 + t315 * t539 + t316 * t360 + t481;
t247 = t316 * t363 + t362 * t493;
t192 = -t247 * t372 + t315 * t376;
t152 = qJD(2) * t385 + t186;
t169 = -pkin(9) * t428 + t187;
t45 = t152 * t377 - t154 * t474 - t373 * t169 - t171 * t473;
t198 = -t372 * t224 - t462;
t396 = -t376 * t224 + t463;
t44 = t373 * t152 + t154 * t473 + t377 * t169 - t171 * t474;
t41 = pkin(10) * t453 + t44;
t163 = -qJD(2) * t382 + qJD(4) * t400;
t164 = qJD(2) * t383 + qJD(4) * t224;
t304 = t319 * qJD(2);
t253 = pkin(3) * t428 + t304;
t84 = t164 * pkin(4) - t163 * pkin(10) + t253;
t7 = t124 * t471 + t372 * t84 + t376 * t41 - t472 * t95;
t391 = pkin(3) * t430 - t313 * t539 - t314 * t360 + t437;
t390 = t254 * t420;
t242 = t314 * t362 + t363 * t458;
t18 = -pkin(4) * t293 - t20;
t381 = -mrSges(3,2) - t601;
t42 = -pkin(4) * t453 - t45;
t8 = -qJD(5) * t49 - t372 * t41 + t376 * t84;
t344 = Ifges(3,4) * t454;
t337 = t370 * t376;
t336 = t370 * t372;
t318 = t369 * t549 - t353;
t317 = (-mrSges(3,1) * t378 + mrSges(3,2) * t374) * t367;
t299 = -t352 * mrSges(3,2) + t432;
t286 = t362 * t369 + t363 * t494;
t285 = t362 * t494 - t369 * t363;
t262 = Ifges(3,1) * t455 + t344 + t514;
t246 = t316 * t362 - t363 * t493;
t237 = -mrSges(4,1) * t454 - t280 * mrSges(4,3);
t236 = mrSges(4,2) * t454 + t279 * mrSges(4,3);
t214 = pkin(5) * t500 - t618;
t196 = -mrSges(4,1) * t656 - mrSges(4,3) * t241;
t195 = mrSges(4,2) * t656 + mrSges(4,3) * t240;
t193 = t247 * t376 + t315 * t372;
t183 = Ifges(4,4) * t280 + Ifges(4,2) * t279 - Ifges(4,6) * t454;
t182 = -Ifges(4,3) * t454 + t521 + t523;
t172 = -mrSges(5,2) * t335 + t538;
t142 = -qJ(6) * t500 + t158;
t141 = t241 * Ifges(4,1) + t240 * Ifges(4,4) - Ifges(4,5) * t656;
t140 = t241 * Ifges(4,4) + t240 * Ifges(4,2) - Ifges(4,6) * t656;
t134 = pkin(5) * t321 - qJ(6) * t499 + t157;
t132 = -mrSges(5,1) * t434 + mrSges(5,2) * t401;
t127 = t515 + t516 + t517;
t121 = mrSges(6,1) * t201 - t535;
t120 = mrSges(7,1) * t201 - t533;
t119 = -mrSges(6,2) * t201 + t536;
t118 = -mrSges(7,2) * t201 + t534;
t103 = qJD(5) * t396 - t372 * t163 + t376 * t453;
t102 = qJD(5) * t198 + t376 * t163 + t372 * t453;
t100 = -mrSges(7,1) * t167 + mrSges(7,2) * t168;
t97 = -mrSges(5,2) * t293 + mrSges(5,3) * t114;
t96 = mrSges(5,1) * t293 - mrSges(5,3) * t113;
t67 = -pkin(5) * t198 + t94;
t39 = qJ(6) * t198 + t49;
t34 = -mrSges(6,2) * t112 + mrSges(6,3) * t65;
t33 = -mrSges(7,2) * t112 + mrSges(7,3) * t65;
t32 = mrSges(6,1) * t112 - mrSges(6,3) * t64;
t31 = mrSges(7,1) * t112 - mrSges(7,3) * t64;
t28 = -pkin(5) * t400 + qJ(6) * t396 + t48;
t24 = -pkin(5) * t103 + t42;
t23 = -mrSges(6,1) * t65 + mrSges(6,2) * t64;
t9 = -pkin(5) * t65 + qJDD(6) + t18;
t6 = qJ(6) * t103 + qJD(6) * t198 + t7;
t5 = pkin(5) * t164 - qJ(6) * t102 + qJD(6) * t396 + t8;
t10 = [m(4) * (t125 * t210 + t126 * t211 + t176 * t186 + t177 * t187 + t208 * t290 + t254 * t304) + m(5) * (t147 * t232 + t19 * t99 + t20 * t98 + t202 * t253 + t44 * t78 + t45 * t77) + m(7) * (t1 * t28 + t2 * t39 + t21 * t5 + t24 * t53 + t26 * t6 + t67 * t9) + m(6) * (t18 * t94 + t29 * t8 + t3 * t49 + t30 * t7 + t4 * t48 + t42 * t74) + (-Ifges(5,6) * t584 - Ifges(5,5) * t585 - Ifges(5,3) * t562 - t461 / 0.2e1 - t125 * mrSges(4,1) + t230 * mrSges(3,3) + t126 * mrSges(4,2) - Ifges(4,6) * t240 - t241 * Ifges(4,5) + t653) * t492 + t459 * t555 + (t633 / 0.2e1 + t646 * t573 + t647 * t579 + t648 * t577 + mrSges(7,2) * t53 + mrSges(6,2) * t74 - t29 * mrSges(6,3) - t21 * mrSges(7,3)) * t102 + (-pkin(1) * t317 * t367 + Ifges(2,3)) * qJDD(1) + t126 * t309 * mrSges(4,3) + (t634 / 0.2e1 + t644 * t573 + t645 * t579 + t647 * t577 - mrSges(7,1) * t53 - mrSges(6,1) * t74 + t30 * mrSges(6,3) + t26 * mrSges(7,3)) * t103 + (-t640 / 0.2e1 - t646 * t586 - t647 * t589 - t648 * t590 - t18 * mrSges(6,2) - t9 * mrSges(7,2) + t4 * mrSges(6,3) + t1 * mrSges(7,3)) * t396 - (-Ifges(3,6) * t369 / 0.2e1 + Ifges(4,5) * t560 + Ifges(4,6) * t561 - t319 * mrSges(3,3) + (-pkin(1) * mrSges(3,1) - t532 + (-Ifges(3,2) - Ifges(4,3)) * t378) * t367) * t656 + t231 * (mrSges(3,1) * t369 - mrSges(3,3) * t494) + t141 * t560 + t140 * t561 + (Ifges(4,1) * t310 + Ifges(4,4) * t309) * t564 + (Ifges(4,4) * t310 + Ifges(4,2) * t309) * t565 + (-t18 * mrSges(6,1) - t9 * mrSges(7,1) + t3 * mrSges(6,3) + t2 * mrSges(7,3) + t644 * t586 + t645 * t589 + t647 * t590 + t668) * t198 - t125 * t310 * mrSges(4,3) + (t635 / 0.2e1 - t78 * mrSges(5,3) + t678) * t164 + (-m(4) * (-pkin(2) * t314 + t437) - (-t314 * t368 + t430) * mrSges(4,1) - (t314 * t366 + t368 * t458) * mrSges(4,2) - m(3) * t437 + t314 * mrSges(3,1) - mrSges(3,3) * t458 + t375 * mrSges(2,1) + t551 * mrSges(2,2) - m(7) * (-t243 * t361 + t391) - m(5) * t391 + t243 * mrSges(5,1) - m(6) * (-pkin(4) * t243 + t391) + t541 * t650 - t649 * t651 + t381 * t313 - t386 * t242) * g(1) + (-t551 * mrSges(2,1) - m(7) * (t247 * t361 + t403) - m(6) * (pkin(4) * t247 + t403) - m(5) * t403 - t247 * mrSges(5,1) + (-mrSges(3,1) - t388) * t316 + (mrSges(2,2) + (-mrSges(3,3) - t420) * t367) * t375 + t541 * t193 - t649 * t192 - t381 * t315 + t386 * t246 + (-m(4) - m(3)) * t481) * g(2) + t39 * t33 + t28 * t31 + t208 * (-mrSges(4,1) * t309 + mrSges(4,2) * t310) + t303 * t299 + t290 * t160 + t253 * t132 + t187 * t236 + t186 * t237 + t232 * t54 + t210 * t196 + t211 * t195 + t45 * t173 + t44 * t172 + (Ifges(3,3) * t555 + t318 * mrSges(3,1) - t319 * mrSges(3,2) + (Ifges(3,5) * t374 + Ifges(3,6) * t378) * t367) * t351 + (Ifges(3,5) * t555 - t318 * mrSges(3,3) + (-pkin(1) * mrSges(3,2) + t374 * Ifges(3,1) + t531) * t367) * t306 + t48 * t32 + t49 * t34 - t369 * t659 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t594 + t230 * t319 + t231 * t318 - t301 * t304 + t302 * t303) + t67 * t22 + (-mrSges(5,3) * t20 + t671) * t224 + (-mrSges(5,3) * t77 + t672) * t163 - t673 * t400 + t94 * t23 + ((t279 * t411 / 0.2e1 + t280 * t414 / 0.2e1 + t514 / 0.2e1 + t390 + t262 / 0.2e1 - t366 * t183 / 0.2e1 - t301 * mrSges(3,3) + t491 / 0.2e1 + (-t176 * t368 - t177 * t366) * mrSges(4,3)) * t378 + (t521 / 0.2e1 + t523 / 0.2e1 + t176 * mrSges(4,1) - t177 * mrSges(4,2) - t522 / 0.2e1 + t515 / 0.2e1 + t517 / 0.2e1 + t516 / 0.2e1 - t302 * mrSges(3,3) + t182 / 0.2e1 + t127 / 0.2e1 - t599) * t374 + ((-Ifges(3,2) * t374 + t531) * t552 - t378 * (Ifges(4,5) * t490 - Ifges(4,6) * t495 + t520) / 0.2e1 - t607) * t480) * t477 + t98 * t96 + t99 * t97 + t24 * t100 + t42 * t101 + t624 * t304 + t6 * t118 + t7 * t119 + t5 * t120 + t8 * t121; (t312 / 0.2e1 - t259 / 0.2e1) * t635 + t183 * t429 / 0.2e1 + t208 * t421 + ((t520 + (Ifges(4,5) * t368 - Ifges(4,6) * t366) * t378) * t552 + t607) * qJD(1) ^ 2 * t594 + (Ifges(5,5) * t570 + Ifges(5,6) * t572 + Ifges(5,3) * t557 + t599) * t455 - (Ifges(5,4) * t572 - t129 / 0.2e1 + t597) * t260 + (t236 * t368 - t237 * t366) * qJD(3) - t656 * (Ifges(4,5) * t366 + Ifges(4,6) * t368) / 0.2e1 + (-t20 * t322 + t666 * t78 + t679 * t77) * mrSges(5,3) + (-t608 * t360 * t492 + t317 + (t541 * (t363 * t486 + t372 * t374) - t649 * (-t363 * t489 + t374 * t376) + t652 * t378 + (-t539 * t608 + t601) * t374) * t367) * g(3) + (-Ifges(5,4) * t570 + t128 / 0.2e1 + t595) * t259 + (-t176 * (mrSges(4,1) * t374 - mrSges(4,3) * t490) - t177 * (-mrSges(4,2) * t374 - mrSges(4,3) * t495)) * t480 + (t195 * t368 - t196 * t366) * qJ(3) + (Ifges(4,1) * t366 + t529) * t564 + (Ifges(4,2) * t368 + t530) * t565 + t499 * t660 + t678 * t312 + t633 * (t322 * t439 - t505 / 0.2e1 - t222 / 0.2e1) + t634 * (-t221 / 0.2e1 - t452 / 0.2e1 + t506 / 0.2e1) - t659 - t390 * t454 + (t432 - t299) * t301 - ((-Ifges(3,2) * t455 + t262 + t344 + t491) * t378 + (t182 + t127) * t374 + t280 * (Ifges(4,5) * t374 + t378 * t414) + t279 * (Ifges(4,6) * t374 + t378 * t411) + t352 * (Ifges(3,5) * t378 - Ifges(3,6) * t374)) * t480 / 0.2e1 - (t23 - t96) * t618 + (t157 * t4 + t158 * t3 - t18 * t618 + t29 * t637 + t30 * t636 + t621 * t74) * m(6) + (-t147 * t360 + t19 * t258 + t20 * t618 - t202 * t252 + t622 * t77 + t623 * t78) * m(5) + t368 * t140 / 0.2e1 + t366 * t141 / 0.2e1 - t360 * t54 - t252 * t132 + t258 * t97 - t216 * t236 - t215 * t237 + t231 * mrSges(3,1) + t214 * t22 + (t221 * t645 + t222 * t647) * t580 + (-t392 * t647 - t393 * t645) * t579 + (t221 * t647 + t222 * t648) * t578 + (-t392 * t648 - t393 * t647) * t577 + (t221 * t644 + t222 * t646) * t574 + (-t392 * t646 - t393 * t644) * t573 + t636 * t119 + t637 * t121 + t638 * t118 + t639 * t120 + (t1 * t134 + t142 * t2 + t21 * t639 + t214 * t9 + t26 * t638 + t53 * t626) * m(7) - t641 * t500 / 0.2e1 + (-t501 * t593 - t608 * (-t315 * t360 + t316 * t539) + t541 * (-t315 * t496 + t501) - t649 * (t315 * t497 + t316 * t376) + t605 * t316 + t600 * t315) * g(1) + (-t502 * t593 - t608 * (-t313 * t360 + t314 * t539) + t541 * (-t313 * t496 + t502) - t649 * (t313 * t497 + t314 * t376) + t605 * t314 + t600 * t313) * g(2) + (t18 * t417 + t415 * t9 + t586 * t615 + t589 * t614 + t590 * t613 + t671) * t322 - t672 * t311 + t673 * t321 + t610 * mrSges(4,3) + (-pkin(2) * t208 + (-t176 * t366 + t177 * t368) * qJD(3) + t610 * qJ(3) - t176 * t215 - t177 * t216) * m(4) + (mrSges(6,1) * t620 - mrSges(6,2) * t619) * t74 + (mrSges(7,1) * t620 - mrSges(7,2) * t619) * t53 + (-t1 * t499 - t2 * t500 + t21 * t619 - t26 * t620) * mrSges(7,3) + (t29 * t619 - t3 * t500 - t30 * t620 - t4 * t499) * mrSges(6,3) + t621 * t101 + t622 * t173 + t623 * t172 + (-m(4) * t254 + t433 - t624) * t302 + t134 * t31 + t142 * t33 + t157 * t32 + t626 * t100 + t158 * t34 - pkin(2) * t160 + t459; -t434 * t172 - t279 * t236 + t280 * t237 + (-t100 + t625) * t401 + (t31 + t32 + t201 * (t118 + t119)) * t376 + (t33 + t34 - t201 * (t120 + t121)) * t372 + t54 + t160 + (t1 * t376 + t2 * t372 - t401 * t53 + t201 * (-t21 * t372 + t26 * t376)) * m(7) + (-t401 * t74 + t3 * t372 + t4 * t376 + t201 * (-t29 * t372 + t30 * t376)) * m(6) + (t401 * t77 - t434 * t78 + t147) * m(5) + (t176 * t280 - t177 * t279 + t208) * m(4) + (-g(1) * t315 - g(2) * t313 + g(3) * t492) * (m(4) + t608); (t376 * t645 + t654) * t589 + (t372 * t648 + t655) * t590 + t9 * t416 + t18 * t418 + (-pkin(4) * t18 - t29 * t46 - t30 * t47) * m(6) + (-t1 * t372 + t2 * t376 + t21 * t612 - t26 * t665) * mrSges(7,3) + (t29 * t612 - t30 * t665 + t609) * mrSges(6,3) - t653 + t595 * t401 + (t200 + t129) * t572 + t128 * t569 + t372 * t660 + (-t508 / 0.2e1 + t471 / 0.2e1) * t633 - t361 * t22 - pkin(4) * t23 + t336 * t31 - t337 * t33 + t376 * t668 + (t372 * t646 + t376 * t644) * t586 + (-t528 + t635) * t570 + (t538 - t172) * t77 + (t167 * t614 + t168 * t613 + t201 * t615) * qJD(5) / 0.2e1 + (t242 * t604 + t243 * t386) * g(2) + (t246 * t604 + t247 * t386) * g(1) + (t285 * t604 + t286 * t386) * g(3) + (m(6) * ((-t29 * t376 - t30 * t372) * qJD(5) + t609) - t121 * t471 - t119 * t472 - t372 * t32 + t376 * t34) * pkin(10) + (t574 * t615 + t578 * t613 + t580 * t614 + t597 - t616) * t434 + t616 * qJD(5) + (-m(6) * t74 + t537 + t625) * t78 - t47 * t119 - t46 * t121 + (t509 / 0.2e1 + t439) * t634 + t628 * t100 + t629 * t120 + t630 * t118 + (t1 * t336 - t2 * t337 + t21 * t629 + t26 * t630 - t361 * t9 + t53 * t628) * m(7) + t461; (t167 * t648 - t674) * t578 + (-m(7) * t548 - mrSges(7,1) * t168 - mrSges(7,2) * t167) * t53 + t21 * t534 + (t536 - t119) * t29 + (-t606 * t651 - t649 * t650) * g(2) + (-t649 * (-t286 * t376 + t463) + t606 * (-t286 * t372 - t462)) * g(3) - t100 * t548 + t634 * t577 + (t535 + t121) * t30 + (t167 * t646 - t168 * t644) * t574 + (t533 - m(7) * (-t21 + t25) + t120) * t26 + (-t168 * t645 + t633 + t675) * t580 + t627 * t1 + pkin(5) * t31 - t74 * (mrSges(6,1) * t168 + mrSges(6,2) * t167) - t25 * t118 + t598 + (t192 * t606 + t193 * t649) * g(1) + t642; -t167 * t118 + t168 * t120 + (-g(1) * t246 - g(2) * t242 - g(3) * t285 - t167 * t26 + t168 * t21 + t9) * m(7) + t22;];
tau  = t10;
