% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:50:07
% EndTime: 2019-03-09 08:50:29
% DurationCPUTime: 14.21s
% Computational Cost: add. (43145->643), mult. (84971->887), div. (0->0), fcn. (101064->10), ass. (0->332)
t385 = sin(pkin(11));
t386 = cos(pkin(11));
t563 = sin(qJ(5));
t565 = cos(qJ(5));
t365 = -t385 * t565 - t386 * t563;
t387 = sin(qJ(6));
t388 = cos(qJ(6));
t622 = t385 * t563 - t386 * t565;
t324 = -t365 * t388 - t387 * t622;
t442 = t365 * t387 - t388 * t622;
t680 = mrSges(7,1) * t324 + mrSges(7,2) * t442;
t677 = qJD(6) * t680;
t507 = sin(pkin(10));
t508 = cos(pkin(10));
t564 = sin(qJ(2));
t566 = cos(qJ(2));
t363 = -t507 * t566 - t508 * t564;
t404 = t365 * t363;
t666 = t622 * t363;
t220 = -t387 * t404 + t388 * t666;
t656 = t220 / 0.2e1;
t473 = Ifges(7,5) * t442 - Ifges(7,6) * t324;
t454 = t507 * pkin(2);
t373 = t454 + qJ(4);
t554 = pkin(8) + t373;
t355 = t554 * t385;
t356 = t554 * t386;
t302 = -t355 * t563 + t356 * t565;
t262 = -pkin(9) * t622 + t302;
t301 = -t355 * t565 - t356 * t563;
t409 = pkin(9) * t365 + t301;
t160 = t262 * t388 + t387 * t409;
t619 = -t262 * t387 + t388 * t409;
t673 = -t160 * mrSges(7,1) - mrSges(7,2) * t619;
t32 = t473 + t673;
t681 = t32 * qJD(6);
t580 = -t324 / 0.2e1;
t581 = t442 / 0.2e1;
t670 = -t387 * t666 - t388 * t404;
t627 = t670 * t581;
t679 = (-t220 * t580 + t627) * mrSges(7,3);
t640 = t670 / 0.2e1;
t361 = t507 * t564 - t508 * t566;
t531 = t220 * Ifges(7,4);
t109 = Ifges(7,2) * t670 + t361 * Ifges(7,6) + t531;
t170 = -mrSges(7,2) * t361 + mrSges(7,3) * t670;
t172 = mrSges(7,1) * t361 - mrSges(7,3) * t220;
t467 = t564 * pkin(7);
t371 = -qJ(3) * t564 - t467;
t469 = t566 * pkin(7);
t372 = qJ(3) * t566 + t469;
t335 = -t371 * t508 + t372 * t507;
t487 = t363 * t385;
t269 = -pkin(4) * t487 + t335;
t222 = pkin(5) * t404 + t269;
t456 = t508 * pkin(2);
t380 = -t456 - pkin(3);
t370 = -pkin(4) * t386 + t380;
t338 = pkin(5) * t622 + t370;
t604 = -t160 / 0.2e1;
t634 = t670 * mrSges(7,2);
t653 = t220 * mrSges(7,1);
t661 = t653 + t634;
t665 = Ifges(7,1) * t670 - t531;
t675 = t619 / 0.2e1;
t676 = (t665 / 0.4e1 - t109 / 0.4e1) * t324 + t170 * t675 + t172 * t604 + t222 * t680 / 0.2e1 + t338 * t661 / 0.2e1;
t584 = t666 / 0.2e1;
t674 = t222 * t661;
t546 = Ifges(7,4) * t324;
t242 = Ifges(7,2) * t442 + t546;
t664 = Ifges(7,1) * t442 - t546;
t672 = t242 / 0.4e1 - t664 / 0.4e1;
t612 = -m(7) / 0.2e1;
t468 = t564 * pkin(2);
t668 = m(4) * t468;
t636 = Ifges(7,5) * t670;
t655 = Ifges(7,6) * t220;
t479 = t636 - t655;
t207 = Ifges(7,4) * t670;
t111 = Ifges(7,1) * t220 + Ifges(7,5) * t361 + t207;
t644 = -Ifges(7,2) * t220 + t207;
t663 = t644 + t111;
t657 = -mrSges(7,2) / 0.2e1;
t658 = mrSges(7,1) / 0.2e1;
t66 = t653 / 0.2e1 - t220 * t658 + t670 * t657 + t634 / 0.2e1;
t451 = t655 / 0.2e1 - t636 / 0.2e1;
t660 = -mrSges(4,1) * t363 - t361 * mrSges(4,2);
t65 = 0.2e1 * mrSges(7,1) * t656 + 0.2e1 * mrSges(7,2) * t640;
t601 = -t220 / 0.2e1;
t576 = -t363 / 0.2e1;
t574 = -t622 / 0.2e1;
t571 = t365 / 0.2e1;
t582 = -t442 / 0.2e1;
t598 = -t670 / 0.2e1;
t320 = -pkin(3) * t363 + qJ(4) * t361 + t468;
t247 = t320 * t386 + t335 * t385;
t248 = t320 * t385 - t335 * t386;
t648 = -t247 * t385 + t248 * t386;
t287 = t361 * t365;
t290 = t361 * t622;
t219 = -t287 * t387 + t290 * t388;
t594 = t219 / 0.2e1;
t215 = -t287 * t388 - t290 * t387;
t599 = t215 / 0.2e1;
t414 = Ifges(7,5) * t594 + Ifges(7,6) * t599;
t488 = t361 * t386;
t174 = -pkin(4) * t363 + pkin(8) * t488 + t247;
t489 = t361 * t385;
t209 = pkin(8) * t489 + t248;
t105 = t174 * t565 - t209 * t563;
t68 = -pkin(5) * t363 - pkin(9) * t290 + t105;
t106 = t174 * t563 + t209 * t565;
t73 = -pkin(9) * t287 + t106;
t55 = -t387 * t73 + t388 * t68;
t56 = t387 * t68 + t388 * t73;
t437 = Ifges(7,3) * t576 + t55 * t658 + t56 * t657 + t414;
t381 = -pkin(2) * t566 - pkin(1);
t315 = pkin(3) * t361 + qJ(4) * t363 + t381;
t623 = t371 * t507 + t372 * t508;
t245 = t315 * t386 - t385 * t623;
t486 = t363 * t386;
t173 = pkin(4) * t361 + pkin(8) * t486 + t245;
t246 = t315 * t385 + t386 * t623;
t208 = pkin(8) * t487 + t246;
t104 = t173 * t563 + t208 * t565;
t71 = -pkin(9) * t404 + t104;
t514 = t387 * t71;
t103 = t173 * t565 - t208 * t563;
t70 = -pkin(9) * t666 + t103;
t67 = pkin(5) * t361 + t70;
t51 = t388 * t67 - t514;
t513 = t388 * t71;
t52 = t387 * t67 + t513;
t423 = -t324 * t51 + t442 * t52;
t528 = t404 * mrSges(6,3);
t251 = -mrSges(6,2) * t361 - t528;
t253 = mrSges(6,1) * t361 - mrSges(6,3) * t666;
t440 = t170 * t581 + t172 * t580 + t251 * t574 + t253 * t571;
t59 = -t387 * t70 - t513;
t60 = t388 * t70 - t514;
t647 = (t324 * t60 + t442 * t59 + t423) * t612 - t440;
t319 = Ifges(7,4) * t442;
t643 = -Ifges(7,2) * t324 + t319;
t587 = -t404 / 0.2e1;
t384 = t386 ^ 2;
t472 = t385 ^ 2 + t384;
t633 = t472 * mrSges(5,3);
t632 = t365 * t404;
t629 = t404 * t574;
t524 = t442 * mrSges(7,3);
t628 = t524 * t598;
t626 = -Ifges(6,5) * t622 + Ifges(6,6) * t365 + t473;
t359 = Ifges(6,4) * t622;
t331 = -Ifges(6,1) * t365 - t359;
t625 = Ifges(6,2) * t365 + t331 - t359;
t282 = Ifges(6,4) * t404;
t180 = Ifges(6,1) * t666 + Ifges(6,5) * t361 - t282;
t624 = -Ifges(6,2) * t666 + t180 - t282;
t621 = -Ifges(6,5) * t404 - Ifges(6,6) * t666 + t479;
t620 = t51 * t582 + t52 * t580;
t617 = 0.2e1 * t363;
t616 = -m(5) / 0.2e1;
t615 = m(5) / 0.2e1;
t614 = -m(6) / 0.2e1;
t613 = m(6) / 0.2e1;
t611 = m(7) / 0.2e1;
t610 = m(4) * pkin(2);
t609 = m(7) * pkin(5);
t607 = -mrSges(6,2) / 0.2e1;
t605 = -t619 / 0.2e1;
t591 = t242 / 0.2e1;
t244 = Ifges(7,1) * t324 + t319;
t590 = t244 / 0.2e1;
t589 = -t287 / 0.2e1;
t586 = t290 / 0.2e1;
t585 = -t666 / 0.2e1;
t583 = t666 / 0.4e1;
t579 = t324 / 0.2e1;
t577 = t361 / 0.2e1;
t575 = t622 / 0.4e1;
t573 = -t365 / 0.2e1;
t572 = -t365 / 0.4e1;
t570 = t385 / 0.2e1;
t569 = t386 / 0.2e1;
t568 = -t387 / 0.2e1;
t567 = t388 / 0.2e1;
t562 = pkin(5) * t666;
t561 = pkin(5) * t365;
t560 = t51 * mrSges(7,2);
t559 = t52 * mrSges(7,1);
t556 = t59 * mrSges(7,1);
t555 = t60 * mrSges(7,2);
t553 = mrSges(4,3) * t363;
t552 = mrSges(6,3) * t365;
t551 = mrSges(7,3) * t324;
t550 = Ifges(5,4) * t385;
t549 = Ifges(5,4) * t386;
t548 = Ifges(6,4) * t666;
t547 = Ifges(6,4) * t365;
t542 = pkin(5) * qJD(5);
t536 = t215 * mrSges(7,1);
t533 = t219 * mrSges(7,2);
t530 = t287 * mrSges(6,1);
t529 = t404 * mrSges(6,2);
t527 = t290 * mrSges(6,2);
t108 = Ifges(7,4) * t219 + Ifges(7,2) * t215 - Ifges(7,6) * t363;
t110 = Ifges(7,1) * t219 + Ifges(7,4) * t215 - Ifges(7,5) * t363;
t116 = t533 - t536;
t117 = -mrSges(7,1) * t670 + mrSges(7,2) * t220;
t169 = mrSges(7,2) * t363 + mrSges(7,3) * t215;
t171 = -mrSges(7,1) * t363 - mrSges(7,3) * t219;
t177 = Ifges(6,4) * t290 - Ifges(6,2) * t287 - Ifges(6,6) * t363;
t178 = -Ifges(6,2) * t404 + Ifges(6,6) * t361 + t548;
t179 = Ifges(6,1) * t290 - Ifges(6,4) * t287 - Ifges(6,5) * t363;
t268 = -pkin(4) * t489 + t623;
t221 = pkin(5) * t287 + t268;
t223 = t527 + t530;
t224 = mrSges(6,1) * t404 + mrSges(6,2) * t666;
t250 = mrSges(6,2) * t363 - mrSges(6,3) * t287;
t252 = -mrSges(6,1) * t363 - mrSges(6,3) * t290;
t518 = t385 * Ifges(5,2);
t263 = -Ifges(5,6) * t363 + (t518 - t549) * t361;
t264 = -Ifges(5,5) * t363 + (-Ifges(5,1) * t386 + t550) * t361;
t516 = t386 * mrSges(5,2);
t519 = t385 * mrSges(5,1);
t435 = -t516 - t519;
t313 = t435 * t361;
t314 = t435 * t363;
t325 = mrSges(5,2) * t363 + mrSges(5,3) * t489;
t326 = -mrSges(5,2) * t361 + mrSges(5,3) * t487;
t327 = -mrSges(5,1) * t363 + mrSges(5,3) * t488;
t328 = mrSges(5,1) * t361 + mrSges(5,3) * t486;
t415 = Ifges(6,5) * t586 + Ifges(6,6) * t589;
t515 = t386 * Ifges(5,5);
t517 = t385 * Ifges(5,6);
t3 = (Ifges(3,1) - Ifges(3,2)) * t566 * t564 + t179 * t584 + (-t564 ^ 2 + t566 ^ 2) * Ifges(3,4) + t110 * t656 + m(5) * (t245 * t247 + t246 * t248 + t335 * t623) + (t660 + t668) * t381 + t335 * t313 + t246 * t325 + t248 * t326 + t245 * t327 + t247 * t328 + t268 * t224 + t269 * t223 + t104 * t250 + t106 * t251 + t103 * t252 + t105 * t253 + t221 * t117 + t222 * t116 + t52 * t169 + t56 * t170 + t51 * t171 + t55 * t172 - pkin(1) * (mrSges(3,1) * t564 + mrSges(3,2) * t566) + t623 * t314 + t108 * t640 + (mrSges(4,1) * t468 + (-Ifges(5,3) + t384 * Ifges(5,1) / 0.2e1 - Ifges(6,3) - Ifges(7,3) + Ifges(4,1) - Ifges(4,2) + (-t549 + t518 / 0.2e1) * t385) * t363 + t414 + t415 + (Ifges(4,4) - t515 + t517) * t361) * t361 + (Ifges(7,5) * t601 + Ifges(7,6) * t598 + Ifges(6,5) * t585 + Ifges(6,6) * t404 / 0.2e1 - t386 * t264 / 0.2e1 + t263 * t570 - mrSges(4,2) * t468 + (t515 / 0.2e1 - t517 / 0.2e1 - Ifges(4,4)) * t363) * t363 + t180 * t586 + t177 * t587 + t178 * t589 + t111 * t594 + t109 * t599 + m(6) * (t103 * t105 + t104 * t106 + t268 * t269) + m(7) * (t221 * t222 + t51 * t55 + t52 * t56);
t526 = t3 * qJD(1);
t521 = t361 * mrSges(4,3);
t520 = t622 * mrSges(6,3);
t431 = -Ifges(6,1) * t404 - t548;
t434 = mrSges(6,1) * t666 - t529;
t4 = t59 * t172 + t117 * t562 + t674 + t665 * t656 + t109 * t601 + t60 * t170 + m(7) * (t222 * t562 + t51 * t59 + t52 * t60) + t178 * t585 + t103 * t251 + t269 * t434 - t104 * t253 + t431 * t584 + (-t220 * t52 - t51 * t670) * mrSges(7,3) + (t103 * t404 - t104 * t666) * mrSges(6,3) + t624 * t587 + t621 * t577 + t663 * t640;
t512 = t4 * qJD(1);
t7 = -t52 * t172 + t674 + t479 * t577 + t51 * t170 + (t665 / 0.2e1 - t109 / 0.2e1 - t52 * mrSges(7,3)) * t220 + (-t51 * mrSges(7,3) + t111 / 0.2e1 + t644 / 0.2e1) * t670;
t509 = t7 * qJD(1);
t240 = -mrSges(7,1) * t442 + mrSges(7,2) * t324;
t433 = mrSges(6,1) * t622 - mrSges(6,2) * t365;
t436 = -mrSges(5,1) * t386 + mrSges(5,2) * t385;
t458 = t524 / 0.2e1;
t463 = -t551 / 0.2e1;
t464 = t552 / 0.2e1;
t390 = (-t361 * t373 * t472 - t363 * t380) * t615 + (-t287 * t301 + t290 * t302 - t363 * t370) * t613 + (t160 * t219 + t215 * t619 - t338 * t363) * t611 + (-t361 * t507 + t363 * t508) * t610 / 0.2e1 + t215 * t463 + t219 * t458 - t287 * t464 - t520 * t586 - t361 * t633 / 0.2e1 + (t240 + t436 + t433) * t576;
t392 = (t247 * t386 + t248 * t385) * t615 + (-t105 * t622 - t106 * t365) * t613 + (t324 * t56 + t442 * t55) * t611 + t171 * t581 + t169 * t579 + t252 * t574 + t250 * t573 + t325 * t570 + t327 * t569 + t668 / 0.2e1;
t14 = -t390 + t392 + t660;
t506 = qJD(1) * t14;
t20 = t670 * t170 - t220 * t172 - t404 * t251 - t666 * t253 + m(7) * (-t220 * t51 + t52 * t670) + m(6) * (-t103 * t666 - t104 * t404) + (t386 * t328 + t385 * t326 + m(5) * (t245 * t386 + t246 * t385)) * t363;
t505 = qJD(1) * t20;
t418 = t245 * t385 - t246 * t386;
t482 = t386 * t326;
t483 = t385 * t328;
t490 = t335 * t363;
t16 = t219 * t170 + t215 * t172 + t290 * t251 - t287 * t253 + (-t482 + t483 + t521) * t361 + (-t117 - t224 - t314 + t553) * t363 + m(7) * (t215 * t51 + t219 * t52 - t222 * t363) + m(6) * (-t103 * t287 + t104 * t290 - t269 * t363) + m(5) * (t361 * t418 - t490) + m(4) * (-t361 * t623 - t490);
t504 = t16 * qJD(1);
t411 = t170 * t582 + t172 * t579;
t412 = t536 / 0.2e1 - t533 / 0.2e1;
t17 = (t220 * t579 + t627) * mrSges(7,3) + t411 + t412;
t503 = t17 * qJD(1);
t502 = t220 * t388;
t501 = t670 * t387;
t497 = t666 * t240;
t496 = t302 * t666;
t494 = t324 * t388;
t492 = t442 * t387;
t471 = t609 / 0.2e1;
t466 = -t561 / 0.2e1;
t462 = mrSges(7,3) * t605;
t461 = mrSges(7,3) * t601;
t460 = t528 / 0.2e1;
t453 = t666 * t573;
t452 = -t488 / 0.2e1;
t449 = m(5) * t472;
t447 = -mrSges(6,1) * t365 - mrSges(6,2) * t622;
t441 = -m(7) / 0.4e1 - m(6) / 0.4e1 - m(5) / 0.4e1;
t430 = -Ifges(6,1) * t622 + t547;
t26 = t338 * t680 + (t664 / 0.2e1 - t242 / 0.2e1) * t324 + (t590 + t643 / 0.2e1) * t442;
t393 = (t111 / 0.4e1 + t644 / 0.4e1) * t442 + (mrSges(7,3) * t604 - t672) * t220 + (t462 + t244 / 0.4e1 + t643 / 0.4e1) * t670 + t361 * t473 / 0.4e1 + t676;
t6 = t393 - t437;
t422 = qJD(1) * t6 + qJD(2) * t26;
t398 = (t622 * t666 + t632) * t614 + (-t220 * t442 + t324 * t670) * t612;
t61 = (-t449 / 0.4e1 + t441) * t617 + t398;
t421 = qJD(1) * t61;
t419 = -t222 * t365 + t338 * t666;
t410 = -t680 - t447;
t102 = (t573 + t494 / 0.2e1 - t492 / 0.2e1) * t609 - t410;
t36 = 0.2e1 * t587 * mrSges(6,2) + 0.2e1 * t584 * mrSges(6,1) + (t584 + t502 / 0.2e1 - t501 / 0.2e1) * t609 + t65;
t417 = qJD(1) * t36 + qJD(2) * t102;
t416 = qJD(1) * t65 + qJD(2) * t680;
t408 = (t492 - t494) * t609;
t330 = -Ifges(6,2) * t622 - t547;
t389 = -t269 * t447 / 0.2e1 + t330 * t583 - t301 * t251 / 0.2e1 + t302 * t253 / 0.2e1 + t178 * t572 - t370 * t434 / 0.2e1 + t625 * t404 / 0.4e1 - t626 * t361 / 0.4e1 + t624 * t575 - t663 * t442 / 0.4e1 - (t643 + t244) * t670 / 0.4e1 + t672 * t220 - t676;
t396 = Ifges(6,3) * t576 + t105 * mrSges(6,1) / 0.2e1 + t106 * t607 + t415 + t437;
t399 = t171 * t567 + (t387 * t56 + t388 * t55) * t611 + t387 * t169 / 0.2e1;
t403 = (t52 + t59) * t619 + (-t51 + t60) * t160;
t1 = (t160 * t656 + t59 * t579 + t60 * t582 + t619 * t640 - t620) * mrSges(7,3) + (t496 / 0.2e1 + t301 * t587) * mrSges(6,3) + (t117 * t571 - t497 / 0.2e1 + t419 * t612 + t399) * pkin(5) + t403 * t612 + Ifges(6,4) * t453 + t396 + t389 + (-t632 / 0.4e1 + t666 * t575) * Ifges(6,1);
t22 = -t240 * t561 - t324 * t591 + t330 * t571 + t370 * t447 + t664 * t579 + t430 * t573 + t442 * t590 + t625 * t574 + t581 * t643 + (-m(7) * t561 + t680) * t338;
t407 = -qJD(1) * t1 + qJD(2) * t22;
t400 = -t530 / 0.2e1 - t527 / 0.2e1 + t412;
t395 = (t215 * t388 + t219 * t387) * t471 + t400;
t8 = t679 + (t453 - t629) * mrSges(6,3) + t395 + t647;
t406 = t8 * qJD(1);
t391 = t679 + (-t571 * t666 - t629) * mrSges(6,3) + t418 * t616 + (t103 * t365 - t104 * t622 - t301 * t666 - t302 * t404) * t613 + (t160 * t670 - t220 * t619 + t423) * t611 - t483 / 0.2e1 + t482 / 0.2e1 + t440;
t394 = t221 * t612 + t268 * t614 + t616 * t623 + t400;
t13 = (t516 / 0.2e1 + t519 / 0.2e1) * t361 + t391 + t394;
t38 = (t324 ^ 2 + t442 ^ 2) * mrSges(7,3) + (t365 ^ 2 + t622 ^ 2) * mrSges(6,3) + t633 + m(7) * (t160 * t442 - t324 * t619) + m(6) * (t301 * t365 - t302 * t622) + t373 * t449;
t405 = -qJD(1) * t13 - qJD(2) * t38;
t397 = (t172 * t568 + t170 * t567 + (t220 * t568 + t388 * t598) * mrSges(7,3)) * pkin(5) - t451;
t11 = (-t51 / 0.2e1 + t60 / 0.2e1) * mrSges(7,2) + (-t52 / 0.2e1 - t59 / 0.2e1) * mrSges(7,1) + t397 + t451;
t31 = (t605 + t675) * mrSges(7,2) + (t604 + t160 / 0.2e1) * mrSges(7,1);
t369 = (mrSges(7,1) * t387 + mrSges(7,2) * t388) * pkin(5);
t402 = -qJD(1) * t11 - qJD(2) * t31 + qJD(5) * t369;
t366 = t369 * qJD(6);
t181 = m(7) * t466 + t408 / 0.2e1;
t62 = t363 * t449 / 0.2e1 + t441 * t617 - t398;
t37 = -t529 / 0.2e1 - t404 * t607 + t66 + (t501 - t502 + t666) * t471;
t18 = t220 * t463 - t411 + t412 + t628;
t15 = t390 + t392;
t12 = mrSges(5,2) * t452 - mrSges(5,1) * t489 / 0.2e1 + t391 - t394;
t10 = -t560 / 0.2e1 - t559 / 0.2e1 - t555 / 0.2e1 + t556 / 0.2e1 + t397 - t451;
t9 = t324 * t461 - t460 * t622 + t464 * t666 + t395 + t628 - t647;
t5 = t393 + t437;
t2 = t403 * t611 + t117 * t466 + t431 * t572 + t430 * t583 + t60 * t458 + t59 * t463 - mrSges(6,3) * t496 / 0.2e1 + t301 * t460 + t160 * t461 + t670 * t462 + t396 - t389 + (t399 + t419 * t611 + t497 / 0.2e1) * pkin(5) + t620 * mrSges(7,3);
t19 = [qJD(2) * t3 + qJD(3) * t16 + qJD(4) * t20 + qJD(5) * t4 + qJD(6) * t7, t15 * qJD(3) + t12 * qJD(4) + t2 * qJD(5) + t5 * qJD(6) + t526 + ((m(5) * t648 + t386 * t325 - t385 * t327) * t373 + t648 * mrSges(5,3) + t110 * t579 + t108 * t581 + t263 * t569 + t264 * t570 + t179 * t573 + t177 * t574 + t105 * t552 + t454 * t553 - t106 * t520 + t456 * t521 + t56 * t524 + m(6) * (t105 * t301 + t106 * t302 + t268 * t370) + (Ifges(5,5) * t385 - Ifges(6,5) * t365 + Ifges(7,5) * t324 + Ifges(5,6) * t386 - Ifges(6,6) * t622 + Ifges(7,6) * t442) * t576 + t370 * t223 + t380 * t313 - Ifges(4,5) * t361 + Ifges(4,6) * t363 + mrSges(3,2) * t467 + t338 * t116 + t301 * t252 + t302 * t250 + t221 * t240 + t160 * t169 - Ifges(3,6) * t564 + Ifges(3,5) * t566 + (Ifges(5,1) * t385 + t549) * t452 + (Ifges(5,2) * t386 + t550) * t489 / 0.2e1 - t55 * t551 - mrSges(3,1) * t469 + (m(5) * t380 - t508 * t610 - mrSges(4,1) + t436) * t623 + (-t507 * t610 + mrSges(4,2)) * t335 + m(7) * (t160 * t56 + t221 * t338 + t55 * t619) + t619 * t171 + t331 * t586 + t330 * t589 + t219 * t590 + t215 * t591 + t268 * t433) * qJD(2), t504 + t15 * qJD(2) + 0.2e1 * ((t215 * t442 + t219 * t324) * t611 + (t287 * t622 - t290 * t365) * t613) * qJD(3) + t62 * qJD(4) + t9 * qJD(5) + t18 * qJD(6), qJD(2) * t12 + qJD(3) * t62 + qJD(5) * t37 + qJD(6) * t66 + t505, t512 + t2 * qJD(2) + t9 * qJD(3) + t37 * qJD(4) + (-t104 * mrSges(6,1) - t103 * mrSges(6,2) - t555 + t556 + t621) * qJD(5) + t10 * qJD(6) + (m(7) * (t387 * t60 + t388 * t59) + (-t220 * t387 - t388 * t670) * mrSges(7,3)) * t542, t509 + t5 * qJD(2) + t18 * qJD(3) + t66 * qJD(4) + t10 * qJD(5) + (t479 - t559 - t560) * qJD(6); -qJD(3) * t14 + qJD(4) * t13 - qJD(5) * t1 + qJD(6) * t6 - t526, qJD(4) * t38 + qJD(5) * t22 + qJD(6) * t26, -t506, qJD(5) * t181 - t405, t181 * qJD(4) + (-mrSges(6,1) * t302 - mrSges(6,2) * t301 + t626 + t673) * qJD(5) + t681 + (m(7) * (-t160 * t388 + t387 * t619) + (-t324 * t387 - t388 * t442) * mrSges(7,3)) * t542 + t407, qJD(5) * t32 + t422 + t681; qJD(2) * t14 - qJD(4) * t61 - qJD(5) * t8 - qJD(6) * t17 - t504, t506, 0, -t421 (t408 + t410) * qJD(5) - t677 - t406, -qJD(5) * t680 - t503 - t677; -qJD(2) * t13 + qJD(3) * t61 + qJD(5) * t36 + qJD(6) * t65 - t505, qJD(5) * t102 + t405 + t677, t421, 0, t417, t416; qJD(2) * t1 + qJD(3) * t8 - qJD(4) * t36 + qJD(6) * t11 - t512, -qJD(4) * t102 + qJD(6) * t31 - t407, t406, -t417, -t366, -t366 - t402; -qJD(2) * t6 + qJD(3) * t17 - qJD(4) * t65 - qJD(5) * t11 - t509, -qJD(4) * t680 - qJD(5) * t31 - t422, t503, -t416, t402, 0;];
Cq  = t19;
