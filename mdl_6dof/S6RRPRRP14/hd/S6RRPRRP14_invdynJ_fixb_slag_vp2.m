% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP14_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP14_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP14_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:47
% EndTime: 2019-03-09 13:04:56
% DurationCPUTime: 44.64s
% Computational Cost: add. (13175->965), mult. (31295->1233), div. (0->0), fcn. (23471->10), ass. (0->427)
t331 = cos(qJ(2));
t324 = cos(pkin(6));
t483 = qJD(1) * t324;
t465 = pkin(1) * t483;
t303 = t331 * t465;
t664 = qJD(3) - t303;
t327 = sin(qJ(2));
t323 = sin(pkin(6));
t484 = qJD(1) * t323;
t439 = t327 * t484;
t505 = t323 * t331;
t249 = qJD(2) * t439 - qJDD(1) * t505;
t468 = qJDD(1) * t324;
t307 = qJDD(2) + t468;
t326 = sin(qJ(4));
t330 = cos(qJ(4));
t309 = qJD(2) + t483;
t446 = t331 * t484;
t343 = -t309 * t330 + t326 * t446;
t122 = qJD(4) * t343 + t249 * t330 - t307 * t326;
t118 = qJDD(5) - t122;
t567 = t118 / 0.2e1;
t215 = -t309 * t326 - t330 * t446;
t121 = qJD(4) * t215 + t249 * t326 + t307 * t330;
t325 = sin(qJ(5));
t329 = cos(qJ(5));
t364 = qJD(4) + t439;
t149 = t325 * t364 - t329 * t343;
t480 = qJD(2) * t331;
t250 = (qJD(1) * t480 + qJDD(1) * t327) * t323;
t230 = qJDD(4) + t250;
t62 = qJD(5) * t149 + t325 * t121 - t329 * t230;
t577 = t62 / 0.2e1;
t578 = -t62 / 0.2e1;
t148 = -t325 * t343 - t329 * t364;
t473 = qJD(5) * t148;
t61 = t329 * t121 + t325 * t230 - t473;
t579 = t61 / 0.2e1;
t639 = Ifges(7,4) + Ifges(6,5);
t640 = Ifges(6,1) + Ifges(7,1);
t655 = Ifges(6,4) * t578 + Ifges(7,5) * t577 + t567 * t639 + t579 * t640;
t636 = Ifges(7,2) + Ifges(6,3);
t653 = Ifges(6,6) - Ifges(7,6);
t395 = pkin(4) * t330 + pkin(10) * t326;
t569 = pkin(3) + pkin(8);
t663 = -(-t395 - t569) * t439 + qJD(4) * t395 + t664;
t642 = m(6) + m(7);
t634 = t118 * t636 + t61 * t639 - t62 * t653;
t661 = t634 / 0.2e1;
t660 = -mrSges(3,1) + mrSges(4,2);
t641 = mrSges(6,3) + mrSges(7,2);
t407 = t330 * t439;
t475 = qJD(4) * t330;
t608 = t407 + t475;
t570 = pkin(2) + pkin(9);
t141 = -t309 * t570 + t439 * t569 + t664;
t426 = -qJ(3) * t327 - pkin(1);
t163 = (-t331 * t570 + t426) * t484;
t477 = qJD(4) * t326;
t301 = pkin(8) * t446;
t507 = t323 * t327;
t310 = pkin(8) * t507;
t547 = pkin(1) * t324;
t464 = qJD(2) * t547;
t411 = qJD(1) * t464;
t458 = pkin(1) * t468;
t151 = -qJD(2) * t301 - qJDD(1) * t310 - t327 * t411 + t331 * t458;
t337 = qJDD(3) - t151;
t88 = pkin(3) * t250 - t307 * t570 + t337;
t479 = qJD(3) * t327;
t335 = -qJ(3) * t250 + (-pkin(1) * qJDD(1) - qJD(1) * t479) * t323;
t93 = t249 * t570 + t335;
t22 = t141 * t475 - t163 * t477 + t326 * t88 + t330 * t93;
t18 = pkin(10) * t230 + t22;
t150 = -pkin(8) * t249 + t327 * t458 + t331 * t411;
t120 = -t307 * qJ(3) - t309 * qJD(3) - t150;
t90 = -pkin(3) * t249 - t120;
t27 = -pkin(4) * t122 - pkin(10) * t121 + t90;
t78 = t326 * t141 + t330 * t163;
t72 = pkin(10) * t364 + t78;
t246 = t327 * t465 + t301;
t201 = pkin(3) * t446 + t246;
t288 = t309 * qJ(3);
t158 = t288 + t201;
t79 = -pkin(4) * t215 + pkin(10) * t343 + t158;
t29 = t325 * t79 + t329 * t72;
t4 = -qJD(5) * t29 - t18 * t325 + t27 * t329;
t2 = -pkin(5) * t118 + qJDD(6) - t4;
t659 = mrSges(7,2) * t2 - mrSges(6,3) * t4 + t655;
t211 = qJD(5) - t215;
t658 = Ifges(6,4) * t579 + Ifges(6,6) * t567 - t61 * Ifges(7,5) / 0.2e1 - t118 * Ifges(7,6) / 0.2e1 + (Ifges(6,2) + Ifges(7,3)) * t578;
t552 = t230 / 0.2e1;
t565 = t122 / 0.2e1;
t590 = t4 * mrSges(6,1) - t2 * mrSges(7,1);
t657 = t567 * t636 + t579 * t639 + t590 + Ifges(6,6) * t578 + Ifges(7,6) * t577 - t121 * Ifges(5,4) / 0.2e1 + (-t552 - t230 / 0.2e1) * Ifges(5,6) + (-t565 - t122 / 0.2e1) * Ifges(5,2);
t631 = -t148 * t653 + t149 * t639 + t211 * t636;
t656 = t631 / 0.2e1;
t654 = Ifges(6,4) - Ifges(7,5);
t147 = Ifges(6,4) * t148;
t523 = Ifges(7,5) * t148;
t630 = t149 * t640 + t211 * t639 - t147 + t523;
t421 = -pkin(10) * t330 + qJ(3);
t545 = pkin(4) * t326;
t271 = t421 + t545;
t474 = qJD(4) * t570;
t441 = t330 * t474;
t471 = qJD(5) * t329;
t298 = pkin(2) * t439;
t365 = pkin(9) * t327 - qJ(3) * t331;
t199 = t365 * t484 + t298;
t110 = t330 * t199 + t326 * t201;
t98 = pkin(10) * t446 + t110;
t652 = t271 * t471 + (-t441 - t98) * t329 + t663 * t325;
t499 = t326 * t570;
t615 = t325 * t271 - t329 * t499;
t651 = -qJD(5) * t615 + t325 * t98 + t329 * t663;
t28 = -t325 * t72 + t329 * t79;
t472 = qJD(5) * t325;
t3 = t329 * t18 + t325 * t27 + t471 * t79 - t472 * t72;
t392 = t3 * t329 - t325 * t4;
t650 = -t28 * t471 - t29 * t472 + t392;
t619 = qJD(6) - t28;
t24 = -pkin(5) * t211 + t619;
t25 = qJ(6) * t211 + t29;
t1 = qJ(6) * t118 + qJD(6) * t211 + t3;
t393 = t1 * t329 + t2 * t325;
t649 = t24 * t471 - t25 * t472 + t393;
t146 = Ifges(7,5) * t149;
t65 = Ifges(7,6) * t211 + Ifges(7,3) * t148 + t146;
t527 = Ifges(6,4) * t149;
t68 = -Ifges(6,2) * t148 + Ifges(6,6) * t211 + t527;
t648 = -t68 / 0.2e1 + t65 / 0.2e1;
t647 = -t3 * mrSges(6,2) + t1 * mrSges(7,3);
t566 = t121 / 0.2e1;
t645 = Ifges(5,1) * t566 + Ifges(5,5) * t552;
t644 = mrSges(7,2) * t1 + mrSges(6,3) * t3 + t658;
t580 = Ifges(5,4) * t565 + t645;
t77 = t330 * t141 - t326 * t163;
t71 = -pkin(4) * t364 - t77;
t34 = t148 * pkin(5) - t149 * qJ(6) + t71;
t587 = -mrSges(6,1) * t71 - mrSges(7,1) * t34 + mrSges(7,2) * t25 + mrSges(6,3) * t29;
t535 = mrSges(4,3) - mrSges(3,2);
t638 = Ifges(3,5) - Ifges(4,4);
t637 = Ifges(4,5) - Ifges(3,6);
t21 = mrSges(6,1) * t62 + mrSges(6,2) * t61;
t84 = mrSges(5,1) * t230 - mrSges(5,3) * t121;
t632 = -t21 + t84;
t629 = t364 * Ifges(5,5);
t628 = t364 * Ifges(5,6);
t440 = t570 * t472;
t627 = (qJD(6) + t440) * t326 + t652 + t608 * qJ(6);
t424 = -t325 * t570 - pkin(5);
t626 = -pkin(5) * t407 + t424 * t475 - t651;
t625 = t326 * t440 + t652;
t624 = t325 * t441 + t651;
t501 = t326 * t327;
t453 = t323 * t501;
t410 = t325 * t453;
t207 = qJD(1) * t410 - t329 * t446;
t497 = t327 * t329;
t218 = (t325 * t331 + t326 * t497) * t323;
t208 = qJD(1) * t218;
t366 = pkin(5) * t325 - qJ(6) * t329;
t355 = t570 + t366;
t367 = pkin(5) * t329 + qJ(6) * t325;
t109 = -t326 * t199 + t201 * t330;
t97 = -pkin(4) * t446 - t109;
t623 = -pkin(5) * t207 + qJ(6) * t208 + (qJD(5) * t367 - qJD(6) * t329) * t330 - t355 * t477 - t97;
t622 = -qJD(6) * t325 + t211 * t366 - t78;
t514 = t343 * mrSges(5,3);
t160 = mrSges(5,1) * t364 + t514;
t82 = mrSges(6,1) * t148 + mrSges(6,2) * t149;
t621 = t82 - t160;
t319 = t327 * t547;
t266 = pkin(8) * t505 + t319;
t221 = -t324 * qJ(3) - t266;
t191 = pkin(3) * t505 - t221;
t258 = t324 * t326 + t330 * t505;
t451 = t326 * t505;
t259 = t324 * t330 - t451;
t418 = -t258 * pkin(4) + pkin(10) * t259;
t108 = t191 - t418;
t546 = pkin(1) * t331;
t447 = -pkin(2) - t546;
t166 = pkin(3) * t507 + t310 + (-pkin(9) + t447) * t324;
t486 = pkin(2) * t505 + qJ(3) * t507;
t548 = pkin(1) * t323;
t222 = -t486 - t548;
t313 = pkin(9) * t505;
t192 = t222 - t313;
t101 = t326 * t166 + t330 * t192;
t95 = pkin(10) * t507 + t101;
t620 = t325 * t108 + t329 * t95;
t139 = -mrSges(5,1) * t215 - mrSges(5,2) * t343;
t414 = mrSges(4,1) * t446;
t241 = -mrSges(4,3) * t309 - t414;
t618 = t139 - t241;
t470 = qJD(5) * t330;
t476 = qJD(4) * t329;
t617 = t325 * t470 + t326 * t476 + t208;
t413 = mrSges(3,3) * t439;
t415 = mrSges(4,1) * t439;
t616 = t309 * t660 + t413 + t415;
t386 = -t329 * mrSges(7,1) - t325 * mrSges(7,3);
t388 = mrSges(6,1) * t329 - mrSges(6,2) * t325;
t614 = m(7) * t367 - t386 + t388;
t385 = mrSges(7,1) * t325 - mrSges(7,3) * t329;
t387 = mrSges(6,1) * t325 + mrSges(6,2) * t329;
t613 = t34 * t385 + t387 * t71;
t612 = t325 * t639 + t329 * t653;
t611 = -t325 * t653 + t329 * t639;
t521 = Ifges(7,5) * t329;
t525 = Ifges(6,4) * t329;
t610 = t325 * t640 - t521 + t525;
t522 = Ifges(7,5) * t325;
t526 = Ifges(6,4) * t325;
t609 = t329 * t640 + t522 - t526;
t23 = -t141 * t477 - t163 * t475 - t326 * t93 + t330 * t88;
t605 = t22 * t326 + t23 * t330;
t573 = t68 / 0.2e1;
t604 = t573 - t65 / 0.2e1;
t564 = -t148 / 0.2e1;
t603 = qJD(4) * t564;
t531 = Ifges(3,4) * t327;
t602 = pkin(1) * (mrSges(3,1) * t327 + mrSges(3,2) * t331) - t327 * (Ifges(3,1) * t331 - t531) / 0.2e1;
t601 = mrSges(5,3) - t660;
t600 = mrSges(5,1) + t614;
t30 = -mrSges(7,2) * t62 + mrSges(7,3) * t118;
t33 = -mrSges(6,2) * t118 - mrSges(6,3) * t62;
t598 = (t30 + t33) * t329;
t31 = mrSges(6,1) * t118 - mrSges(6,3) * t61;
t32 = -t118 * mrSges(7,1) + t61 * mrSges(7,2);
t597 = (-t31 + t32) * t325;
t469 = -m(5) - t642;
t596 = pkin(9) * t469 - t601;
t459 = m(4) - t469;
t595 = qJ(3) * t459 + t535;
t328 = sin(qJ(1));
t494 = t328 * t331;
t332 = cos(qJ(1));
t495 = t327 * t332;
t261 = t324 * t495 + t494;
t489 = t331 * t332;
t498 = t327 * t328;
t260 = -t324 * t489 + t498;
t504 = t323 * t332;
t352 = -t260 * t326 + t330 * t504;
t594 = t261 * t329 + t325 * t352;
t593 = -t261 * t325 + t329 * t352;
t444 = t323 * t480;
t445 = qJD(2) * t507;
t300 = pkin(2) * t445;
t162 = t300 + (qJD(2) * t365 - t479) * t323;
t202 = (t505 * t569 + t319) * qJD(2);
t51 = t330 * t162 + t166 * t475 - t192 * t477 + t326 * t202;
t40 = pkin(10) * t444 + t51;
t304 = t331 * t464;
t320 = t324 * qJD(3);
t416 = t569 * t507;
t165 = -qJD(2) * t416 + t304 + t320;
t179 = -qJD(4) * t258 + t326 * t445;
t180 = -qJD(4) * t451 + t324 * t475 - t330 * t445;
t74 = pkin(4) * t180 - pkin(10) * t179 + t165;
t9 = -qJD(5) * t620 - t325 * t40 + t329 * t74;
t417 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t591 = t23 * mrSges(5,1) - t22 * mrSges(5,2) + Ifges(5,5) * t121 + Ifges(5,6) * t122 + Ifges(5,3) * t230;
t408 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t262 = t324 * t494 + t495;
t506 = t323 * t328;
t184 = t262 * t326 + t330 * t506;
t589 = -g(1) * t184 + g(2) * t352 - g(3) * t259;
t389 = mrSges(5,1) * t326 + mrSges(5,2) * t330;
t588 = -t389 - t535 - t642 * t421 + t641 * t330 + (-m(4) - m(5)) * qJ(3);
t586 = m(5) / 0.2e1;
t585 = m(6) / 0.2e1;
t584 = m(7) / 0.2e1;
t530 = Ifges(5,4) * t343;
t114 = Ifges(5,2) * t215 - t530 + t628;
t568 = -t114 / 0.2e1;
t563 = t148 / 0.2e1;
t562 = -t149 / 0.2e1;
t561 = t149 / 0.2e1;
t557 = -t211 / 0.2e1;
t556 = t211 / 0.2e1;
t555 = -t215 / 0.2e1;
t554 = t343 / 0.2e1;
t553 = -t343 / 0.2e1;
t537 = qJD(4) / 0.2e1;
t534 = Ifges(3,4) + Ifges(4,6);
t533 = mrSges(6,3) * t148;
t532 = mrSges(6,3) * t149;
t529 = Ifges(5,4) * t326;
t528 = Ifges(5,4) * t330;
t524 = Ifges(5,5) * t343;
t520 = Ifges(4,6) * t327;
t519 = Ifges(4,6) * t331;
t518 = Ifges(5,6) * t215;
t517 = Ifges(5,3) * t331;
t19 = -pkin(4) * t230 - t23;
t516 = t19 * t330;
t515 = t215 * mrSges(5,3);
t140 = -pkin(4) * t343 - pkin(10) * t215;
t45 = t325 * t140 + t329 * t77;
t511 = t215 * t325;
t510 = t215 * t329;
t503 = t325 * t326;
t500 = t326 * t329;
t496 = t327 * t330;
t492 = t329 * t271;
t102 = -mrSges(7,2) * t148 + mrSges(7,3) * t211;
t103 = -mrSges(6,2) * t211 - t533;
t488 = t102 + t103;
t104 = mrSges(6,1) * t211 - t532;
t105 = -mrSges(7,1) * t211 + mrSges(7,2) * t149;
t487 = -t104 + t105;
t485 = t332 * pkin(1) + pkin(8) * t506;
t481 = qJD(1) ^ 2 * t323 ^ 2;
t478 = qJD(4) * t325;
t463 = pkin(10) * t472;
t462 = pkin(10) * t471;
t461 = Ifges(3,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t460 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t452 = t323 * t496;
t263 = -t324 * t498 + t489;
t450 = t263 * pkin(2) + t485;
t448 = t313 + t486;
t443 = t325 * t477;
t442 = t326 * t474;
t434 = t484 / 0.2e1;
t432 = -t477 / 0.2e1;
t428 = t471 / 0.2e1;
t427 = -t470 / 0.2e1;
t425 = -pkin(1) * t328 + pkin(8) * t504;
t252 = t260 * pkin(2);
t423 = -pkin(9) * t260 - t252;
t254 = t262 * pkin(2);
t422 = -pkin(9) * t262 - t254;
t178 = t250 * mrSges(4,1) + t307 * mrSges(4,2);
t100 = t166 * t330 - t326 * t192;
t412 = mrSges(3,3) * t446;
t409 = pkin(3) * t506 + t450;
t400 = t327 * t434;
t399 = -t261 * pkin(2) + t425;
t391 = mrSges(5,1) * t258 + mrSges(5,2) * t259;
t384 = mrSges(4,2) * t331 - mrSges(4,3) * t327;
t383 = Ifges(5,1) * t326 + t528;
t378 = Ifges(5,2) * t330 + t529;
t377 = -Ifges(6,2) * t325 + t525;
t376 = Ifges(6,2) * t329 + t526;
t373 = Ifges(5,5) * t179 - Ifges(5,6) * t180;
t372 = Ifges(5,5) * t326 + Ifges(5,6) * t330;
t369 = Ifges(7,3) * t325 + t521;
t368 = -Ifges(7,3) * t329 + t522;
t42 = t108 * t329 - t325 * t95;
t44 = t140 * t329 - t325 * t77;
t357 = pkin(3) * t504 + t399;
t245 = pkin(8) * t439 - t303;
t247 = -pkin(8) * t445 + t304;
t52 = -t326 * t162 - t166 * t477 - t192 * t475 + t202 * t330;
t356 = -pkin(10) * t642 + mrSges(5,2) - t641;
t353 = -t259 * t325 + t323 * t497;
t182 = t259 * t329 + t325 * t507;
t185 = t260 * t330 + t326 * t504;
t8 = t108 * t471 + t325 * t74 + t329 * t40 - t472 * t95;
t346 = t327 * (-Ifges(4,2) * t331 + t520);
t345 = t331 * (Ifges(4,3) * t327 - t519);
t94 = -pkin(4) * t507 - t100;
t136 = -pkin(2) * t307 + t337;
t336 = t151 * mrSges(3,1) - t150 * mrSges(3,2) + t136 * mrSges(4,2) - t120 * mrSges(4,3);
t41 = -pkin(4) * t444 - t52;
t297 = Ifges(3,4) * t446;
t287 = Ifges(4,1) * t307;
t286 = Ifges(3,3) * t307;
t272 = -pkin(4) - t367;
t265 = t324 * t546 - t310;
t264 = (-mrSges(3,1) * t331 + mrSges(3,2) * t327) * t323;
t248 = t266 * qJD(2);
t244 = -qJ(3) * t446 + t298;
t243 = t384 * t484;
t240 = -mrSges(3,2) * t309 + t412;
t229 = Ifges(4,4) * t250;
t228 = Ifges(3,5) * t250;
t227 = Ifges(4,5) * t249;
t226 = Ifges(3,6) * t249;
t225 = t324 * t447 + t310;
t220 = t355 * t330;
t212 = t325 * t499 + t492;
t210 = -t247 - t320;
t209 = Ifges(5,4) * t215;
t206 = t309 * t325 - t329 * t407;
t205 = t309 * t329 + t325 * t407;
t204 = (-pkin(2) * t331 + t426) * t484;
t203 = t300 + (-qJ(3) * t480 - t479) * t323;
t200 = -qJD(1) * t416 + t303;
t198 = t326 * t424 - t492;
t197 = -t288 - t246;
t196 = t309 * Ifges(4,4) + (-Ifges(4,2) * t327 - t519) * t484;
t195 = t309 * Ifges(4,5) + (-t331 * Ifges(4,3) - t520) * t484;
t194 = Ifges(3,1) * t439 + t309 * Ifges(3,5) + t297;
t193 = t309 * Ifges(3,6) + (t331 * Ifges(3,2) + t531) * t484;
t190 = qJ(6) * t326 + t615;
t189 = -pkin(2) * t309 + qJD(3) + t245;
t183 = -t262 * t330 + t326 * t506;
t177 = mrSges(4,1) * t249 - mrSges(4,3) * t307;
t159 = -mrSges(5,2) * t364 + t515;
t128 = t184 * t329 + t263 * t325;
t127 = t184 * t325 - t263 * t329;
t123 = pkin(2) * t249 + t335;
t115 = -Ifges(5,1) * t343 + t209 + t629;
t113 = Ifges(5,3) * t364 + t518 - t524;
t92 = qJD(5) * t353 + t179 * t329 + t325 * t444;
t91 = qJD(5) * t182 + t179 * t325 - t329 * t444;
t85 = -mrSges(5,2) * t230 + mrSges(5,3) * t122;
t81 = mrSges(7,1) * t148 - mrSges(7,3) * t149;
t80 = pkin(5) * t149 + qJ(6) * t148;
t64 = -mrSges(5,1) * t122 + mrSges(5,2) * t121;
t53 = -pkin(5) * t353 - qJ(6) * t182 + t94;
t38 = -pkin(5) * t258 - t42;
t37 = qJ(6) * t258 + t620;
t36 = pkin(5) * t343 - t44;
t35 = -qJ(6) * t343 + t45;
t20 = mrSges(7,1) * t62 - mrSges(7,3) * t61;
t10 = pkin(5) * t91 - qJ(6) * t92 - qJD(6) * t182 + t41;
t7 = -pkin(5) * t180 - t9;
t6 = qJ(6) * t180 + qJD(6) * t258 + t8;
t5 = pkin(5) * t62 - qJ(6) * t61 - qJD(6) * t149 + t19;
t11 = [(-mrSges(5,3) * t22 - Ifges(5,4) * t566 + t647 + t657 + t661) * t258 + (mrSges(6,2) * t19 - mrSges(7,3) * t5 + t655 + t659) * t182 + t90 * t391 + t215 * (Ifges(5,4) * t179 - Ifges(5,2) * t180) / 0.2e1 + t165 * t139 + t51 * t159 + t52 * t160 + m(6) * (t19 * t94 + t28 * t9 + t29 * t8 + t3 * t620 + t4 * t42 + t41 * t71) + t620 * t33 + m(3) * (t150 * t266 + t151 * t265 + t245 * t248 + t246 * t247) + (Ifges(7,5) * t92 + Ifges(7,6) * t180) * t563 + (-t179 * t77 - t180 * t78) * mrSges(5,3) + (-mrSges(5,3) * t23 + 0.2e1 * t580) * t259 + (-t180 * t29 + t71 * t92) * mrSges(6,2) + (t228 / 0.2e1 - t226 / 0.2e1 + t286 / 0.2e1 + t287 / 0.2e1 - t229 / 0.2e1 + t227 / 0.2e1 + (Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1) * t307 + t461 * t250 + t460 * t249 + t336) * t324 + t180 * t568 + t100 * t84 + t101 * t85 + t6 * t102 + t8 * t103 + t9 * t104 + t7 * t105 + t94 * t21 + t10 * t81 + t41 * t82 + Ifges(2,3) * qJDD(1) + (-m(3) * t425 + mrSges(2,1) * t328 + mrSges(2,2) * t332 - m(5) * t357 - t352 * mrSges(5,1) - m(4) * t399 + t595 * t260 - t417 * t593 + t408 * t594 - t596 * t261 + t356 * t185 + t642 * (-pkin(4) * t352 - t357)) * g(1) + (-m(3) * t485 - mrSges(2,1) * t332 + mrSges(2,2) * t328 - m(5) * t409 - t184 * mrSges(5,1) - m(4) * t450 - t595 * t262 - t417 * t128 + t408 * t127 + t596 * t263 + t356 * t183 - t642 * (t184 * pkin(4) + t409)) * g(2) + t53 * t20 + ((-mrSges(3,1) * t249 - mrSges(3,2) * t250 + (m(3) * t548 - t264) * qJDD(1)) * pkin(1) + (-t120 * mrSges(4,1) + t123 * mrSges(4,2) + t150 * mrSges(3,3) - t637 * t307 + t534 * t250 + (-Ifges(3,2) - Ifges(4,3)) * t249) * t331 + (-t151 * mrSges(3,3) + t136 * mrSges(4,1) + qJD(1) * t373 / 0.2e1 - t123 * mrSges(4,3) + t638 * t307 + (Ifges(3,1) + Ifges(4,2)) * t250 - t534 * t249 + t591) * t327 + ((-t193 / 0.2e1 + t195 / 0.2e1 - t204 * mrSges(4,2) + t197 * mrSges(4,1) - t246 * mrSges(3,3) + t460 * t309) * t327 + (t194 / 0.2e1 - t196 / 0.2e1 + t113 / 0.2e1 - t204 * mrSges(4,3) + t77 * mrSges(5,1) - t78 * mrSges(5,2) + t518 / 0.2e1 - t524 / 0.2e1 + Ifges(5,3) * t537 + t189 * mrSges(4,1) + t245 * mrSges(3,3) + t461 * t309) * t331 + (-t345 / 0.2e1 + t331 * (Ifges(3,4) * t331 - Ifges(3,2) * t327) / 0.2e1 - t346 / 0.2e1 + t327 * t517 / 0.2e1 - t602) * t484) * qJD(2) + (mrSges(3,3) + mrSges(4,1)) * (-g(1) * t332 - g(2) * t328)) * t323 + m(4) * (t120 * t221 + t123 * t222 + t136 * t225 + t189 * t248 + t197 * t210 + t203 * t204) + m(5) * (t100 * t23 + t101 * t22 + t158 * t165 + t191 * t90 + t51 * t78 + t52 * t77) + m(7) * (t1 * t37 + t10 * t34 + t2 * t38 + t24 * t7 + t25 * t6 + t5 * t53) + t42 * t31 + t37 * t30 + t38 * t32 + (Ifges(6,4) * t92 + Ifges(6,6) * t180) * t564 + t616 * t248 + (t180 * t25 - t34 * t92) * mrSges(7,3) + (Ifges(5,1) * t179 - Ifges(5,4) * t180) * t553 + t373 * t537 + t180 * t656 + t179 * t115 / 0.2e1 + t28 * (mrSges(6,1) * t180 - mrSges(6,3) * t92) + t24 * (-mrSges(7,1) * t180 + mrSges(7,2) * t92) + t158 * (mrSges(5,1) * t180 + mrSges(5,2) * t179) + (-Ifges(6,2) * t564 + Ifges(7,3) * t563 - t587 + t648) * t91 + t191 * t64 + t221 * t177 + t225 * t178 + t210 * t241 + t203 * t243 + t247 * t240 + t222 * (-mrSges(4,2) * t249 - mrSges(4,3) * t250) + t630 * t92 / 0.2e1 + (t180 * t636 + t639 * t92 - t653 * t91) * t556 + (-mrSges(6,1) * t19 - mrSges(7,1) * t5 + Ifges(6,2) * t578 - Ifges(7,3) * t577 + t567 * t653 + t579 * t654 + t644) * t353 + (t180 * t639 + t640 * t92 - t654 * t91) * t561 + t265 * (mrSges(3,1) * t307 - mrSges(3,3) * t250) + t266 * (-mrSges(3,2) * t307 - mrSges(3,3) * t249); (Ifges(7,5) * t208 - Ifges(7,6) * t407 - t376 * t470) * t564 + (t477 * t77 - t605) * mrSges(5,3) - t644 * t325 * t330 + t659 * t329 * t330 + (-pkin(2) * t136 - qJ(3) * t120 - qJD(3) * t197 - t204 * t244) * m(4) + (t325 * t65 + t115) * t432 + t326 * t661 + (Ifges(6,6) * t603 + t369 * t577 + t377 * t578 + t5 * t385 + t65 * t428 + t567 * t611 + t579 * t609 + t580 + t645) * t330 + (t346 + t345) * t481 / 0.2e1 + (-t377 * t603 + t657) * t326 - t529 * t566 + (-Ifges(6,2) * t563 + Ifges(7,3) * t564 - t557 * t653 - t562 * t654 + t587 + t604) * t207 - t632 * t330 * t570 - t197 * t415 + t90 * t389 + t364 * (mrSges(5,1) * t330 - mrSges(5,2) * t326) * t158 + t615 * t33 - ((t327 * t637 + t331 * t638) * t309 + t215 * (Ifges(5,6) * t331 + t327 * t378) - t343 * (Ifges(5,5) * t331 + t327 * t383) + t364 * (t327 * t372 + t517) + (-Ifges(3,2) * t439 + t113 + t194 + t297) * t331 + (t114 * t330 + t115 * t326 + t195) * t327) * t484 / 0.2e1 - (t215 * t378 - t343 * t383 + t364 * t372) * qJD(4) / 0.2e1 + (-t71 * t97 + t4 * t212 + t3 * t615 - (t477 * t71 - t516) * t570 + t625 * t29 + t624 * t28) * m(6) + (-t109 * t77 - t110 * t78 + t90 * qJ(3) - ((-t326 * t77 + t330 * t78) * qJD(4) + t605) * t570 + (-t200 + qJD(3)) * t158) * m(5) + (-t204 * (-mrSges(4,2) * t327 - mrSges(4,3) * t331) - t78 * (-mrSges(5,2) * t331 + mrSges(5,3) * t496) - t77 * (mrSges(5,1) * t331 - mrSges(5,3) * t501)) * t484 + (Ifges(6,4) * t208 - Ifges(6,6) * t407 - t368 * t470 + (Ifges(7,6) * t330 - t326 * t369) * qJD(4)) * t563 + t228 - t229 - t226 + t227 + (-m(4) * t197 + t240 - t241 - t412) * t245 + (t330 * t631 + t193) * t400 + (-t110 - t441) * t159 + (-t442 - t97) * t82 + (-t109 + t442) * t160 + t286 + t287 + t329 * t68 * t427 + t331 * t196 * t434 + t336 - t189 * t414 + t443 * t573 + t528 * t565 + (-m(4) * t486 - m(5) * t448 + t264 + t641 * t452 - t642 * (pkin(4) * t453 - pkin(10) * t452 + t448) + (-t331 * mrSges(5,3) - t327 * t389 + t384) * t323 - t417 * t218 + t408 * (-t329 * t505 + t410)) * g(3) + (m(4) * t252 - m(5) * t423 - t642 * (t261 * t545 + t423) - t417 * (-t260 * t325 + t261 * t500) + t408 * (t260 * t329 + t261 * t503) + t601 * t260 + t588 * t261) * g(2) + (m(4) * t254 - m(5) * t422 - t642 * (t263 * t545 + t422) - t417 * (-t262 * t325 + t263 * t500) + t408 * (t262 * t329 + t263 * t503) + t601 * t262 + t588 * t263) * g(1) + (-t610 * t470 + (-t326 * t609 + t330 * t639) * qJD(4)) * t561 + (t208 * t640 - t407 * t639) * t562 + (-t612 * t470 + (-t326 * t611 + t330 * t636) * qJD(4)) * t556 + (t208 * t639 - t407 * t636) * t557 + t625 * t103 + t626 * t105 + (t1 * t190 + t198 * t2 + t220 * t5 + t24 * t626 + t25 * t627 + t34 * t623) * m(7) + t627 * t102 + t623 * t81 + t624 * t104 + (mrSges(6,1) * t608 + mrSges(6,3) * t617) * t28 + (t1 * t326 + t25 * t608 + t34 * t617) * mrSges(7,3) + (-t29 * t608 - t3 * t326 - t617 * t71) * mrSges(6,2) + (-mrSges(7,1) * t608 - mrSges(7,2) * t617) * t24 + t618 * qJD(3) + (-m(4) * t189 + t413 - t616) * t246 + t602 * t481 + t587 * (-t329 * t470 + t443) + t387 * t516 - pkin(2) * t178 + t190 * t30 + t198 * t32 - t200 * t139 + t212 * t31 + t220 * t20 - t244 * t243 - t85 * t499 + t630 * (t325 * t427 + t329 * t432 - t208 / 0.2e1) + (-t78 * mrSges(5,3) + t568 + t656) * t475 + (t64 - t177) * qJ(3); -t488 * t206 + t487 * t205 + t243 * t439 + (t159 * t439 - t20 + (t325 * t487 + t329 * t488 + t159) * qJD(4) + t632) * t330 + (t85 + t598 + t597 + (-t325 * t488 + t329 * t487) * qJD(5) + t364 * (t81 + t621)) * t326 - m(6) * (t205 * t28 + t206 * t29) - m(7) * (-t205 * t24 + t206 * t25) + 0.2e1 * (m(5) * t78 * t400 + (qJD(4) * t78 + t23) * t586 + (t24 * t478 + t25 * t476 - t5) * t584 + (-t28 * t478 + t29 * t476 - t19) * t585) * t330 + 0.2e1 * ((-qJD(4) * t77 + t22) * t586 + (qJD(4) * t34 + t649) * t584 + (qJD(4) * t71 + t650) * t585 + (-m(5) * t77 / 0.2e1 + t34 * t584 + t71 * t585) * t439) * t326 + t178 + (-g(1) * t262 - g(2) * t260 + g(3) * t505) * t459 + (-m(5) * t158 - t618) * t309 + (t197 * t309 + t204 * t439 + t136) * m(4); t658 * t329 + (-t510 / 0.2e1 + t428) * t630 + t5 * t386 - t19 * t388 + (-t463 - t35) * t102 + (t149 * t609 + t211 * t611) * qJD(5) / 0.2e1 - (-t28 * mrSges(6,1) + t24 * mrSges(7,1) + t628 / 0.2e1 + t29 * mrSges(6,2) - t25 * mrSges(7,3) + Ifges(6,6) * t563 + Ifges(7,6) * t564 - Ifges(5,2) * t555 - t158 * mrSges(5,1) + t639 * t562 + t636 * t557) * t343 + (-pkin(4) * t19 + ((-t28 * t329 - t29 * t325) * qJD(5) + t392) * pkin(10) - t28 * t44 - t29 * t45) * m(6) + t591 + (-t463 - t45) * t103 + (t209 + t115) * t555 + (t462 - t36) * t105 + (-t462 - t44) * t104 + (-t159 + t515) * t77 + t368 * t577 + t376 * t578 + (-t377 / 0.2e1 + t369 / 0.2e1) * t473 + (t258 * t614 - t418 * t642 + t391) * g(3) + (mrSges(5,2) * t184 - t642 * (-t183 * pkin(4) + pkin(10) * t184) + t600 * t183) * g(1) + (-mrSges(5,2) * t352 - t642 * (t185 * pkin(4) - pkin(10) * t352) - t600 * t185) * g(2) + (t530 + t631) * t554 + (-t629 / 0.2e1 + t377 * t563 + t369 * t564 + Ifges(5,1) * t554 - t158 * mrSges(5,2) + t609 * t562 + t611 * t557 - t613) * t215 + (-m(6) * t71 - t514 - t621) * t78 + t622 * t81 + (t272 * t5 + ((t24 * t329 - t25 * t325) * qJD(5) + t393) * pkin(10) - t24 * t36 - t25 * t35 + t622 * t34) * m(7) - pkin(4) * t21 + t612 * t567 + t613 * qJD(5) + t610 * t579 + t604 * t511 + t114 * t553 + pkin(10) * t598 + t325 * t655 + pkin(10) * t597 + t648 * t472 + (-t24 * t510 + t25 * t511 + t589 + t649) * mrSges(7,2) + (t28 * t510 + t29 * t511 + t589 + t650) * mrSges(6,3) + t272 * t20; t647 + t634 + (t182 * t408 - t353 * t417) * g(3) + (-t408 * t593 - t417 * t594) * g(2) + t590 + (-t487 + t532) * t29 + (-t488 - t533) * t28 + (t148 * t24 + t149 * t25) * mrSges(7,2) + (t127 * t417 + t128 * t408) * g(1) - t34 * (mrSges(7,1) * t149 + mrSges(7,3) * t148) - t71 * (mrSges(6,1) * t149 - mrSges(6,2) * t148) + qJD(6) * t102 - t80 * t81 + (-t148 * t639 - t149 * t653) * t557 + (-t148 * t640 + t146 - t527 + t65) * t562 + (-Ifges(6,2) * t149 - t147 + t630) * t563 - pkin(5) * t32 + qJ(6) * t30 + (-pkin(5) * t2 + qJ(6) * t1 - t24 * t29 + t25 * t619 - t34 * t80) * m(7) + (Ifges(7,3) * t149 - t523) * t564 + t68 * t561; -t211 * t102 + t149 * t81 + (-g(1) * t127 + g(2) * t594 + g(3) * t353 + t34 * t149 - t25 * t211 + t2) * m(7) + t32;];
tau  = t11;
