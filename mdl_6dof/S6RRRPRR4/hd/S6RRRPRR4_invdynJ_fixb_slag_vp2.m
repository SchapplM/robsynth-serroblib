% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:15:25
% EndTime: 2019-03-09 18:16:12
% DurationCPUTime: 28.08s
% Computational Cost: add. (28205->921), mult. (60984->1190), div. (0->0), fcn. (46328->18), ass. (0->421)
t463 = sin(qJ(3));
t468 = cos(qJ(3));
t469 = cos(qJ(2));
t548 = qJD(1) * t469;
t464 = sin(qJ(2));
t549 = qJD(1) * t464;
t367 = t463 * t549 - t468 * t548;
t458 = sin(pkin(11));
t459 = cos(pkin(11));
t462 = sin(qJ(5));
t467 = cos(qJ(5));
t390 = t458 * t467 + t459 * t462;
t275 = t390 * t367;
t364 = t390 * qJD(5);
t711 = t275 + t364;
t560 = t459 * t467;
t490 = t458 * t462 - t560;
t276 = t490 * t367;
t363 = t490 * qJD(5);
t710 = t276 + t363;
t455 = pkin(11) + qJ(5);
t448 = qJ(6) + t455;
t435 = sin(t448);
t446 = sin(t455);
t730 = -mrSges(6,2) * t446 - mrSges(7,2) * t435;
t729 = -mrSges(6,3) - mrSges(7,3);
t393 = t463 * t469 + t464 * t468;
t368 = t393 * qJD(1);
t300 = pkin(3) * t368 + qJ(4) * t367;
t471 = -pkin(8) - pkin(7);
t416 = t471 * t469;
t397 = qJD(1) * t416;
t369 = t463 * t397;
t415 = t471 * t464;
t396 = qJD(1) * t415;
t378 = qJD(2) * pkin(2) + t396;
t306 = t378 * t468 + t369;
t213 = t459 * t300 - t306 * t458;
t573 = t367 * t459;
t508 = t368 * pkin(4) + pkin(9) * t573;
t162 = t213 + t508;
t214 = t458 * t300 + t459 * t306;
t574 = t367 * t458;
t537 = pkin(9) * t574;
t183 = t537 + t214;
t460 = -pkin(9) - qJ(4);
t408 = t460 * t458;
t451 = t459 * pkin(9);
t409 = qJ(4) * t459 + t451;
t541 = qJD(5) * t467;
t543 = qJD(4) * t459;
t544 = qJD(4) * t458;
t689 = t408 * t541 + (-t183 + t543) * t467 + (-qJD(5) * t409 - t162 - t544) * t462;
t320 = t462 * t408 + t467 * t409;
t688 = -t390 * qJD(4) - qJD(5) * t320 - t467 * t162 + t183 * t462;
t530 = pkin(2) * t549;
t284 = t300 + t530;
t312 = t396 * t468 + t369;
t211 = t459 * t284 - t312 * t458;
t156 = t211 + t508;
t212 = t458 * t284 + t459 * t312;
t179 = t537 + t212;
t617 = pkin(2) * t463;
t437 = qJ(4) + t617;
t373 = (-pkin(9) - t437) * t458;
t374 = t437 * t459 + t451;
t299 = t462 * t373 + t467 * t374;
t545 = qJD(3) * t468;
t527 = pkin(2) * t545;
t426 = qJD(4) + t527;
t687 = -qJD(5) * t299 - t467 * t156 + t179 * t462 - t390 * t426;
t567 = t426 * t458;
t686 = -t467 * t179 + t373 * t541 + t426 * t560 + (-qJD(5) * t374 - t156 - t567) * t462;
t728 = t711 * pkin(10);
t436 = cos(t448);
t447 = cos(t455);
t727 = mrSges(6,1) * t447 + mrSges(7,1) * t436;
t726 = -t368 * pkin(5) + t710 * pkin(10);
t456 = qJD(2) + qJD(3);
t323 = t368 * t459 + t456 * t458;
t511 = -t368 * t458 + t459 * t456;
t238 = t323 * t467 + t462 * t511;
t461 = sin(qJ(6));
t466 = cos(qJ(6));
t705 = -t323 * t462 + t467 * t511;
t725 = -t238 * t461 + t466 * t705;
t144 = t238 * t466 + t461 * t705;
t724 = m(5) * qJ(4) + mrSges(5,3);
t723 = -t728 + t686;
t722 = t726 + t687;
t721 = -t728 + t689;
t720 = t726 + t688;
t283 = -pkin(3) * t456 + qJD(4) - t306;
t452 = t469 * pkin(2);
t442 = t452 + pkin(1);
t413 = t442 * qJD(1);
t589 = Ifges(5,4) * t459;
t499 = -Ifges(5,2) * t458 + t589;
t590 = Ifges(5,4) * t458;
t501 = Ifges(5,1) * t459 - t590;
t503 = mrSges(5,1) * t458 + mrSges(5,2) * t459;
t621 = t459 / 0.2e1;
t694 = t456 * Ifges(4,5);
t717 = (t323 * Ifges(5,1) + Ifges(5,4) * t511 + t367 * Ifges(5,5)) * t621 - t458 * (Ifges(5,4) * t323 + Ifges(5,2) * t511 + Ifges(5,6) * t367) / 0.2e1 + t283 * t503 - t306 * mrSges(4,3) + t323 * t501 / 0.2e1 - t413 * mrSges(4,2) + t511 * t499 / 0.2e1 + t694 / 0.2e1;
t457 = qJ(2) + qJ(3);
t449 = sin(t457);
t597 = mrSges(5,2) * t458;
t716 = (-t597 + t730) * t449;
t715 = -t724 + t729;
t277 = pkin(3) * t367 - qJ(4) * t368 - t413;
t559 = t468 * t397;
t307 = t463 * t378 - t559;
t288 = qJ(4) * t456 + t307;
t194 = t459 * t277 - t288 * t458;
t195 = t458 * t277 + t459 * t288;
t361 = qJD(5) + t367;
t136 = pkin(4) * t367 - pkin(9) * t323 + t194;
t152 = pkin(9) * t511 + t195;
t79 = t467 * t136 - t152 * t462;
t58 = -pkin(10) * t238 + t79;
t55 = pkin(5) * t361 + t58;
t80 = t136 * t462 + t152 * t467;
t59 = pkin(10) * t705 + t80;
t580 = t461 * t59;
t20 = t466 * t55 - t580;
t579 = t466 * t59;
t21 = t461 * t55 + t579;
t691 = t511 * Ifges(5,6);
t693 = t456 * Ifges(4,6);
t695 = t323 * Ifges(5,5);
t714 = t413 * mrSges(4,1) - t194 * mrSges(5,1) - t79 * mrSges(6,1) - t20 * mrSges(7,1) + t195 * mrSges(5,2) + t80 * mrSges(6,2) + t21 * mrSges(7,2) - t691 / 0.2e1 + t693 / 0.2e1 - t695 / 0.2e1;
t539 = qJD(1) * qJD(2);
t402 = qJDD(1) * t469 - t464 * t539;
t403 = qJDD(1) * t464 + t469 * t539;
t392 = t463 * t464 - t468 * t469;
t482 = t392 * qJD(3);
t264 = -qJD(1) * t482 + t402 * t463 + t403 * t468;
t483 = t393 * qJD(3);
t265 = qJD(1) * t483 - t468 * t402 + t403 * t463;
t578 = qJDD(1) * pkin(1);
t356 = -pkin(2) * t402 - t578;
t137 = pkin(3) * t265 - qJ(4) * t264 - qJD(4) * t368 + t356;
t385 = t403 * pkin(7);
t317 = qJDD(2) * pkin(2) - pkin(8) * t403 - t385;
t384 = t402 * pkin(7);
t325 = pkin(8) * t402 + t384;
t546 = qJD(3) * t463;
t176 = t463 * t317 + t468 * t325 + t378 * t545 + t397 * t546;
t453 = qJDD(2) + qJDD(3);
t167 = qJ(4) * t453 + qJD(4) * t456 + t176;
t82 = t459 * t137 - t167 * t458;
t83 = t458 * t137 + t459 * t167;
t496 = -t458 * t82 + t459 * t83;
t601 = mrSges(5,1) * t459;
t504 = t597 - t601;
t709 = t711 * pkin(5);
t352 = qJD(6) + t361;
t708 = t238 * Ifges(6,5) + t144 * Ifges(7,5) + Ifges(6,6) * t705 + Ifges(7,6) * t725 + t367 * Ifges(5,3) + t361 * Ifges(6,3) + t352 * Ifges(7,3) + t691 + t695;
t450 = cos(t457);
t707 = -t450 * mrSges(4,1) + (mrSges(4,2) + t729) * t449;
t706 = (t601 + t727) * t449;
t438 = pkin(4) * t459 + pkin(3);
t678 = t450 * t438 - t449 * t460;
t391 = pkin(5) * t447 + t438;
t454 = -pkin(10) + t460;
t680 = t450 * t391 - t449 * t454;
t704 = -m(6) * t678 - m(7) * t680;
t241 = -t264 * t458 + t453 * t459;
t242 = t264 * t459 + t453 * t458;
t110 = qJD(5) * t705 + t241 * t462 + t242 * t467;
t111 = -qJD(5) * t238 + t241 * t467 - t242 * t462;
t37 = qJD(6) * t725 + t110 * t466 + t111 * t461;
t657 = t37 / 0.2e1;
t38 = -qJD(6) * t144 - t110 * t461 + t111 * t466;
t656 = t38 / 0.2e1;
t701 = -m(6) - m(7);
t649 = t110 / 0.2e1;
t648 = t111 / 0.2e1;
t636 = t241 / 0.2e1;
t635 = t242 / 0.2e1;
t263 = qJDD(5) + t265;
t258 = qJDD(6) + t263;
t634 = t258 / 0.2e1;
t633 = t263 / 0.2e1;
t632 = t265 / 0.2e1;
t700 = t402 / 0.2e1;
t620 = t469 / 0.2e1;
t298 = t467 * t373 - t374 * t462;
t605 = pkin(10) * t390;
t255 = t298 - t605;
t380 = t490 * pkin(10);
t256 = -t380 + t299;
t165 = t255 * t466 - t256 * t461;
t699 = qJD(6) * t165 + t461 * t722 + t466 * t723;
t166 = t255 * t461 + t256 * t466;
t698 = -qJD(6) * t166 - t461 * t723 + t466 * t722;
t319 = t467 * t408 - t409 * t462;
t280 = t319 - t605;
t281 = -t380 + t320;
t191 = t280 * t466 - t281 * t461;
t697 = qJD(6) * t191 + t461 * t720 + t466 * t721;
t192 = t280 * t461 + t281 * t466;
t696 = -qJD(6) * t192 - t461 * t721 + t466 * t720;
t692 = t469 * Ifges(3,2);
t660 = m(7) * pkin(5);
t690 = mrSges(6,1) + t660;
t305 = pkin(3) * t392 - qJ(4) * t393 - t442;
t331 = t415 * t463 - t416 * t468;
t226 = t459 * t305 - t331 * t458;
t569 = t393 * t459;
t182 = pkin(4) * t392 - pkin(9) * t569 + t226;
t227 = t458 * t305 + t459 * t331;
t570 = t393 * t458;
t203 = -pkin(9) * t570 + t227;
t107 = t462 * t182 + t467 * t203;
t311 = t396 * t463 - t559;
t345 = pkin(4) * t574;
t260 = t311 - t345;
t528 = pkin(2) * t546;
t682 = t528 - t260 + t709;
t584 = t368 * mrSges(4,3);
t681 = mrSges(4,1) * t456 + mrSges(5,1) * t511 - mrSges(5,2) * t323 - t584;
t679 = t468 * t415 + t416 * t463;
t607 = pkin(7) * t469;
t608 = pkin(7) * t464;
t676 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t549) * t607 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t548) * t608;
t163 = -mrSges(5,2) * t265 + mrSges(5,3) * t241;
t164 = mrSges(5,1) * t265 - mrSges(5,3) * t242;
t675 = t459 * t163 - t458 * t164;
t674 = t384 * t469 + t385 * t464;
t465 = sin(qJ(1));
t470 = cos(qJ(1));
t673 = g(1) * t470 + g(2) * t465;
t251 = -t345 + t307;
t672 = -t251 + t709;
t491 = -t438 * t449 - t450 * t460;
t493 = -t391 * t449 - t450 * t454;
t614 = pkin(3) * t449;
t616 = pkin(2) * t464;
t671 = -m(7) * (t493 - t616) - m(6) * (t491 - t616) - m(5) * (-t614 - t616) + t706;
t670 = -m(3) * pkin(7) + m(5) * t471 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - t503;
t433 = t449 * mrSges(5,3);
t669 = -t433 + t707 + (t504 - t727 - t730) * t450;
t51 = pkin(4) * t265 - pkin(9) * t242 + t82;
t63 = pkin(9) * t241 + t83;
t16 = -qJD(5) * t80 - t462 * t63 + t467 * t51;
t7 = pkin(5) * t263 - pkin(10) * t110 + t16;
t542 = qJD(5) * t462;
t15 = t136 * t541 - t152 * t542 + t462 * t51 + t467 * t63;
t8 = pkin(10) * t111 + t15;
t3 = qJD(6) * t20 + t461 * t7 + t466 * t8;
t4 = -qJD(6) * t21 - t461 * t8 + t466 * t7;
t668 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t667 = m(5) * t283 - t681;
t666 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t564 = t450 * t465;
t665 = t465 * t716 + t564 * t715;
t563 = t450 * t470;
t664 = t470 * t716 + t563 * t715;
t663 = m(5) * t614 - m(6) * t491 - m(7) * t493 + t706;
t412 = -mrSges(3,1) * t469 + mrSges(3,2) * t464;
t662 = -(m(5) * pkin(3) - t504) * t450 - mrSges(2,1) - m(3) * pkin(1) + t412 + t707;
t659 = Ifges(7,4) * t657 + Ifges(7,2) * t656 + Ifges(7,6) * t634;
t658 = Ifges(7,1) * t657 + Ifges(7,4) * t656 + Ifges(7,5) * t634;
t655 = Ifges(6,4) * t649 + Ifges(6,2) * t648 + Ifges(6,6) * t633;
t654 = Ifges(6,1) * t649 + Ifges(6,4) * t648 + Ifges(6,5) * t633;
t587 = Ifges(7,4) * t144;
t72 = Ifges(7,2) * t725 + Ifges(7,6) * t352 + t587;
t653 = -t72 / 0.2e1;
t652 = t72 / 0.2e1;
t139 = Ifges(7,4) * t725;
t73 = Ifges(7,1) * t144 + Ifges(7,5) * t352 + t139;
t651 = -t73 / 0.2e1;
t650 = t73 / 0.2e1;
t647 = Ifges(5,1) * t635 + Ifges(5,4) * t636 + Ifges(5,5) * t632;
t588 = Ifges(6,4) * t238;
t130 = Ifges(6,2) * t705 + Ifges(6,6) * t361 + t588;
t646 = t130 / 0.2e1;
t233 = Ifges(6,4) * t705;
t131 = Ifges(6,1) * t238 + Ifges(6,5) * t361 + t233;
t645 = t131 / 0.2e1;
t644 = -t725 / 0.2e1;
t643 = t725 / 0.2e1;
t642 = -t144 / 0.2e1;
t641 = t144 / 0.2e1;
t640 = -t705 / 0.2e1;
t639 = t705 / 0.2e1;
t638 = -t238 / 0.2e1;
t637 = t238 / 0.2e1;
t630 = -t352 / 0.2e1;
t629 = t352 / 0.2e1;
t628 = -t361 / 0.2e1;
t627 = t361 / 0.2e1;
t626 = -t367 / 0.2e1;
t625 = t367 / 0.2e1;
t623 = t368 / 0.2e1;
t619 = mrSges(7,3) * t20;
t618 = mrSges(7,3) * t21;
t615 = pkin(2) * t468;
t613 = pkin(4) * t458;
t612 = pkin(5) * t238;
t609 = pkin(5) * t446;
t602 = g(3) * t449;
t594 = mrSges(6,3) * t705;
t593 = mrSges(6,3) * t238;
t592 = Ifges(3,4) * t464;
t591 = Ifges(3,4) * t469;
t586 = pkin(5) * qJD(6);
t583 = t368 * Ifges(4,4);
t313 = -qJD(2) * t392 - t482;
t576 = t313 * t458;
t575 = t313 * t459;
t430 = t449 * qJ(4);
t314 = qJD(2) * t393 + t483;
t547 = qJD(2) * t464;
t529 = pkin(2) * t547;
t204 = pkin(3) * t314 - qJ(4) * t313 - qJD(4) * t393 + t529;
t521 = qJD(2) * t471;
t398 = t464 * t521;
t399 = t469 * t521;
t243 = qJD(3) * t679 + t398 * t468 + t399 * t463;
t128 = t458 * t204 + t459 * t243;
t333 = t435 * t564 + t436 * t470;
t334 = t435 * t470 - t436 * t564;
t558 = -t333 * mrSges(7,1) + t334 * mrSges(7,2);
t335 = -t435 * t563 + t436 * t465;
t336 = t435 * t465 + t436 * t563;
t557 = t335 * mrSges(7,1) - t336 * mrSges(7,2);
t550 = t450 * pkin(3) + t430;
t534 = Ifges(7,5) * t37 + Ifges(7,6) * t38 + Ifges(7,3) * t258;
t522 = Ifges(6,5) * t110 + Ifges(6,6) * t111 + Ifges(6,3) * t263;
t13 = -t38 * mrSges(7,1) + t37 * mrSges(7,2);
t515 = t539 / 0.2e1;
t149 = -t241 * mrSges(5,1) + t242 * mrSges(5,2);
t46 = -t111 * mrSges(6,1) + t110 * mrSges(6,2);
t106 = t467 * t182 - t203 * t462;
t127 = t459 * t204 - t243 * t458;
t279 = pkin(4) * t570 - t679;
t507 = mrSges(3,1) * t464 + mrSges(3,2) * t469;
t505 = mrSges(4,1) * t449 + mrSges(4,2) * t450;
t502 = -mrSges(7,1) * t435 - mrSges(7,2) * t436;
t500 = t592 + t692;
t498 = Ifges(3,5) * t469 - Ifges(3,6) * t464;
t497 = Ifges(5,5) * t459 - Ifges(5,6) * t458;
t292 = t490 * t393;
t84 = pkin(5) * t392 + pkin(10) * t292 + t106;
t291 = t390 * t393;
t88 = -pkin(10) * t291 + t107;
t39 = -t461 * t88 + t466 * t84;
t40 = t461 * t84 + t466 * t88;
t495 = -t194 * t458 + t195 * t459;
t205 = -t291 * t466 + t292 * t461;
t206 = -t291 * t461 - t292 * t466;
t308 = -t390 * t461 - t466 * t490;
t309 = t390 * t466 - t461 * t490;
t344 = pkin(5) * t490 - t438;
t177 = t317 * t468 - t463 * t325 - t378 * t546 + t397 * t545;
t486 = t534 + t668;
t485 = pkin(1) * t507;
t342 = -t446 * t563 + t447 * t465;
t340 = t446 * t564 + t447 * t470;
t484 = t464 * (Ifges(3,1) * t469 - t592);
t113 = -pkin(9) * t576 + t128;
t99 = pkin(4) * t314 - pkin(9) * t575 + t127;
t24 = t467 * t113 + t182 * t541 - t203 * t542 + t462 * t99;
t168 = -pkin(3) * t453 + qJDD(4) - t177;
t223 = -pkin(4) * t511 + t283;
t25 = -qJD(5) * t107 - t113 * t462 + t467 * t99;
t244 = qJD(3) * t331 + t398 * t463 - t468 * t399;
t121 = -pkin(4) * t241 + t168;
t189 = pkin(4) * t576 + t244;
t117 = Ifges(5,4) * t242 + Ifges(5,2) * t241 + Ifges(5,6) * t265;
t138 = -pkin(5) * t705 + t223;
t187 = t275 * t466 - t276 * t461;
t188 = t275 * t461 + t276 * t466;
t215 = qJD(6) * t308 - t363 * t466 - t364 * t461;
t216 = -qJD(6) * t309 + t363 * t461 - t364 * t466;
t286 = -t367 * Ifges(4,2) + t583 + t693;
t358 = Ifges(4,4) * t367;
t287 = t368 * Ifges(4,1) - t358 + t694;
t52 = -pkin(5) * t111 + t121;
t472 = (Ifges(7,4) * t188 + Ifges(7,2) * t187) * t644 + t121 * (mrSges(6,1) * t490 + mrSges(6,2) * t390) + (Ifges(6,5) * t390 - Ifges(6,6) * t490) * t633 + (Ifges(6,4) * t390 - Ifges(6,2) * t490) * t648 + (Ifges(6,1) * t390 - Ifges(6,4) * t490) * t649 - t490 * t655 + (-t187 * t21 + t188 * t20 + t3 * t308 - t309 * t4) * mrSges(7,3) + (Ifges(7,1) * t641 + Ifges(7,4) * t643 + Ifges(7,5) * t629 - t619 + t650) * t215 + (Ifges(7,4) * t641 + Ifges(7,2) * t643 + Ifges(7,6) * t629 + t618 + t652) * t216 + (-Ifges(6,1) * t363 - Ifges(6,4) * t364) * t637 + (-Ifges(6,4) * t363 - Ifges(6,2) * t364) * t639 + (-Ifges(6,5) * t363 - Ifges(6,6) * t364) * t627 + (Ifges(6,1) * t276 + Ifges(6,4) * t275) * t638 + (Ifges(6,5) * t276 + Ifges(6,6) * t275) * t628 + (Ifges(7,5) * t188 + Ifges(7,6) * t187) * t630 + ((-t188 + t215) * mrSges(7,2) + (t187 - t216) * mrSges(7,1)) * t138 + (Ifges(7,1) * t188 + Ifges(7,4) * t187) * t642 + t168 * t504 + (Ifges(6,4) * t276 + Ifges(6,2) * t275) * t640 + (Ifges(6,5) * t638 + Ifges(7,5) * t642 - Ifges(4,2) * t625 + Ifges(6,6) * t640 + Ifges(7,6) * t644 + Ifges(5,3) * t626 + Ifges(6,3) * t628 + Ifges(7,3) * t630 + t714) * t368 + (-t358 + t287) * t625 + Ifges(4,3) * t453 + t52 * (-mrSges(7,1) * t308 + mrSges(7,2) * t309) - t275 * t130 / 0.2e1 - t276 * t131 / 0.2e1 + (mrSges(6,1) * t711 - mrSges(6,2) * t710) * t223 + (-t15 * t490 - t16 * t390 + t710 * t79 - t711 * t80) * mrSges(6,3) + (-t194 * t573 - t195 * t574 + t496) * mrSges(5,3) + Ifges(4,5) * t264 - Ifges(4,6) * t265 - (-Ifges(4,1) * t367 - t583 + t708) * t368 / 0.2e1 + t307 * t584 + (-t497 * t626 + t717) * t367 - t176 * mrSges(4,2) + t177 * mrSges(4,1) + t117 * t621 + t286 * t623 + (Ifges(5,5) * t458 + Ifges(5,6) * t459) * t632 + (Ifges(7,5) * t309 + Ifges(7,6) * t308) * t634 + (Ifges(5,1) * t458 + t589) * t635 + (Ifges(5,2) * t459 + t590) * t636 - t363 * t645 - t364 * t646 + t458 * t647 + t188 * t651 + t187 * t653 + t390 * t654 + (Ifges(7,4) * t309 + Ifges(7,2) * t308) * t656 + (Ifges(7,1) * t309 + Ifges(7,4) * t308) * t657 + t309 * t658 + t308 * t659;
t444 = Ifges(3,4) * t548;
t441 = -pkin(3) - t615;
t419 = t470 * t442;
t407 = -t438 - t615;
t404 = t609 + t613;
t366 = Ifges(3,1) * t549 + Ifges(3,5) * qJD(2) + t444;
t365 = Ifges(3,6) * qJD(2) + qJD(1) * t500;
t343 = t446 * t465 + t447 * t563;
t341 = t446 * t470 - t447 * t564;
t338 = -mrSges(4,2) * t456 - mrSges(4,3) * t367;
t332 = t344 - t615;
t302 = mrSges(4,1) * t367 + mrSges(4,2) * t368;
t269 = mrSges(5,1) * t367 - mrSges(5,3) * t323;
t268 = -mrSges(5,2) * t367 + mrSges(5,3) * t511;
t249 = -mrSges(4,2) * t453 - mrSges(4,3) * t265;
t248 = mrSges(4,1) * t453 - mrSges(4,3) * t264;
t210 = mrSges(6,1) * t361 - t593;
t209 = -mrSges(6,2) * t361 + t594;
t207 = pkin(5) * t291 + t279;
t161 = -t313 * t390 + t363 * t393;
t160 = -t313 * t490 - t364 * t393;
t147 = -mrSges(6,1) * t705 + mrSges(6,2) * t238;
t126 = mrSges(7,1) * t352 - mrSges(7,3) * t144;
t125 = -mrSges(7,2) * t352 + mrSges(7,3) * t725;
t100 = -pkin(5) * t161 + t189;
t86 = -mrSges(6,2) * t263 + mrSges(6,3) * t111;
t85 = mrSges(6,1) * t263 - mrSges(6,3) * t110;
t78 = -mrSges(7,1) * t725 + mrSges(7,2) * t144;
t61 = -qJD(6) * t206 - t160 * t461 + t161 * t466;
t60 = qJD(6) * t205 + t160 * t466 + t161 * t461;
t29 = -mrSges(7,2) * t258 + mrSges(7,3) * t38;
t28 = mrSges(7,1) * t258 - mrSges(7,3) * t37;
t23 = t466 * t58 - t580;
t22 = -t461 * t58 - t579;
t19 = pkin(10) * t161 + t24;
t18 = pkin(5) * t314 - pkin(10) * t160 + t25;
t6 = -qJD(6) * t40 + t18 * t466 - t19 * t461;
t5 = qJD(6) * t39 + t18 * t461 + t19 * t466;
t1 = [(t402 * t607 + t403 * t608 + t674) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t674) + (t366 * t620 + t498 * qJD(2) / 0.2e1 - t676) * qJD(2) + (-t341 * mrSges(6,1) - t334 * mrSges(7,1) - t340 * mrSges(6,2) - t333 * mrSges(7,2) + (-m(6) * (-t471 + t613) + m(4) * t471 - m(7) * (t404 - t471) + t670) * t470 + (-m(6) * (-t442 - t678) + m(4) * t442 - m(7) * (-t442 - t680) - m(5) * (-t442 - t430) + t433 - t662) * t465) * g(1) + m(4) * (t176 * t331 + t243 * t307 - t356 * t442 - t413 * t529) + (t82 * mrSges(5,1) - t83 * mrSges(5,2) - Ifges(4,4) * t264 + Ifges(4,2) * t265 - Ifges(4,6) * t453 + t356 * mrSges(4,1) - t176 * mrSges(4,3) + Ifges(5,3) * t632 + Ifges(6,3) * t633 + Ifges(7,3) * t634 + Ifges(5,5) * t635 + Ifges(5,6) * t636 + Ifges(6,6) * t648 + Ifges(6,5) * t649 + Ifges(7,6) * t656 + Ifges(7,5) * t657 + t666 + t668) * t392 + m(5) * (t127 * t194 + t128 * t195 + t226 * t82 + t227 * t83) + (-m(4) * t306 + t667) * t244 + (Ifges(6,5) * t160 + Ifges(6,6) * t161) * t627 - (-m(4) * t177 + m(5) * t168 + t149 - t248) * t679 + (Ifges(7,1) * t60 + Ifges(7,4) * t61) * t641 + (Ifges(7,1) * t206 + Ifges(7,4) * t205) * t657 + (t356 * mrSges(4,2) - t177 * mrSges(4,3) + Ifges(4,1) * t264 - Ifges(4,4) * t265 + Ifges(4,5) * t453 + t168 * t503 + t497 * t632 + t499 * t636 + t501 * t635) * t393 + (Ifges(5,5) * t242 + Ifges(5,6) * t241 + Ifges(5,3) * t265 + t522 + t534) * t392 / 0.2e1 + (-Ifges(6,5) * t292 - Ifges(6,6) * t291) * t633 + (-Ifges(6,4) * t292 - Ifges(6,2) * t291) * t648 + (-t15 * t291 + t16 * t292 - t160 * t79 + t161 * t80) * mrSges(6,3) + t121 * (mrSges(6,1) * t291 - mrSges(6,2) * t292) + (-Ifges(6,1) * t292 - Ifges(6,4) * t291) * t649 + (t469 * t591 + t484) * t515 + (-t20 * t60 + t205 * t3 - t206 * t4 + t21 * t61) * mrSges(7,3) + t500 * t700 + (Ifges(7,4) * t60 + Ifges(7,2) * t61) * t643 + (Ifges(7,4) * t206 + Ifges(7,2) * t205) * t656 - t117 * t570 / 0.2e1 + m(7) * (t100 * t138 + t20 * t6 + t207 * t52 + t21 * t5 + t3 * t40 + t39 * t4) + m(6) * (t106 * t16 + t107 * t15 + t121 * t279 + t189 * t223 + t24 * t80 + t25 * t79) + (Ifges(3,4) * t403 + Ifges(3,2) * t402) * t620 + (Ifges(6,4) * t160 + Ifges(6,2) * t161) * t639 - t365 * t547 / 0.2e1 - t485 * t539 + t403 * t591 / 0.2e1 - t412 * t578 + (Ifges(7,5) * t60 + Ifges(7,6) * t61) * t629 + (Ifges(7,5) * t206 + Ifges(7,6) * t205) * t634 + (-t194 * t575 - t195 * t576 - t569 * t82 - t570 * t83) * mrSges(5,3) + t302 * t529 + (t708 / 0.2e1 - t307 * mrSges(4,3) - Ifges(4,4) * t623 + Ifges(5,3) * t625 - Ifges(4,2) * t626 - t286 / 0.2e1 + Ifges(6,3) * t627 + Ifges(7,3) * t629 + Ifges(6,5) * t637 + Ifges(6,6) * t639 + Ifges(7,5) * t641 + Ifges(7,6) * t643 - t714) * t314 + Ifges(2,3) * qJDD(1) - t442 * (mrSges(4,1) * t265 + mrSges(4,2) * t264) - pkin(1) * (-mrSges(3,1) * t402 + mrSges(3,2) * t403) + t243 * t338 + t331 * t249 + t39 * t28 + t40 * t29 + t279 * t46 + t128 * t268 + t127 * t269 + t100 * t78 + (-mrSges(3,1) * t608 - mrSges(3,2) * t607 + 0.2e1 * Ifges(3,6) * t620) * qJDD(2) + (Ifges(3,1) * t403 + Ifges(3,4) * t700 + Ifges(3,5) * qJDD(2) - t515 * t692) * t464 + t106 * t85 + t107 * t86 + (t287 / 0.2e1 + Ifges(4,1) * t623 + t497 * t625 + Ifges(4,4) * t626 + t717) * t313 + (Ifges(6,1) * t160 + Ifges(6,4) * t161) * t637 + t5 * t125 + t6 * t126 + t138 * (-mrSges(7,1) * t61 + mrSges(7,2) * t60) + (-m(5) * t419 - t343 * mrSges(6,1) - t336 * mrSges(7,1) - t342 * mrSges(6,2) - t335 * mrSges(7,2) + (-m(4) + t701) * (-t465 * t471 + t419) + (-m(6) * t613 - m(7) * t404 + t670) * t465 + (-t449 * t724 + t662 + t704) * t470) * g(2) + t189 * t147 + t52 * (-mrSges(7,1) * t205 + mrSges(7,2) * t206) + t207 * t13 + t24 * t209 + t25 * t210 + t223 * (-mrSges(6,1) * t161 + mrSges(6,2) * t160) + t226 * t164 + t227 * t163 + t160 * t645 + t161 * t646 + t569 * t647 + t60 * t650 + t61 * t652 - t292 * t654 - t291 * t655 + t206 * t658 + t205 * t659; t675 * t437 + (t465 * t671 + t665) * g(2) + (t470 * t671 + t664) * g(1) + (-m(5) * (t452 + t550) - m(4) * t452 - m(6) * (t452 + t678) - m(7) * (t452 + t680) + t412 + t669) * g(3) + ((t176 * t463 + t177 * t468 + (-t306 * t463 + t307 * t468) * qJD(3)) * pkin(2) + t306 * t311 - t307 * t312 + t413 * t530) * m(4) + (m(6) * t223 + t147 + t667) * t528 + (t168 * t441 - t194 * t211 - t195 * t212 - t283 * t311 + t426 * t495 + t437 * t496) * m(5) + t472 + (m(4) * t616 + t505 + t507) * t673 - (-Ifges(3,2) * t549 + t366 + t444) * t548 / 0.2e1 + t686 * t209 + (t121 * t407 + t15 * t299 + t16 * t298 - t223 * t260 + t686 * t80 + t687 * t79) * m(6) + t687 * t210 + (-t567 - t211) * t269 + (t426 * t459 - t212) * t268 + (t527 - t312) * t338 + t365 * t549 / 0.2e1 - t498 * t539 / 0.2e1 + (t676 + (-t484 / 0.2e1 + t485) * qJD(1)) * qJD(1) - t302 * t530 + t698 * t126 + t699 * t125 + (t682 * t138 + t165 * t4 + t166 * t3 + t698 * t20 + t699 * t21 + t332 * t52) * m(7) + t681 * t311 + t682 * t78 + t441 * t149 + Ifges(3,5) * t403 + t407 * t46 + Ifges(3,6) * t402 - t384 * mrSges(3,2) - t385 * mrSges(3,1) + t332 * t13 + Ifges(3,3) * qJDD(2) + t298 * t85 + t299 * t86 - t260 * t147 + t165 * t28 + t166 * t29 + t248 * t615 + t249 * t617; (-pkin(3) * t168 + qJ(4) * t496 + qJD(4) * t495 - t194 * t213 - t195 * t214 - t283 * t307) * m(5) + t675 * qJ(4) + t672 * t78 + t673 * t505 + (t543 - t214) * t268 + t472 + (-t544 - t213) * t269 + t688 * t210 + (-t121 * t438 + t15 * t320 + t16 * t319 - t223 * t251 + t688 * t79 + t689 * t80) * m(6) + t689 * t209 + (t465 * t663 + t665) * g(2) + (t470 * t663 + t664) * g(1) + t696 * t126 + t697 * t125 + (t138 * t672 + t191 * t4 + t192 * t3 + t20 * t696 + t21 * t697 + t344 * t52) * m(7) + t681 * t307 - t438 * t46 + t344 * t13 - t306 * t338 + t319 * t85 + t320 * t86 - t251 * t147 - pkin(3) * t149 + t191 * t28 + t192 * t29 + (-m(5) * t550 + t669 + t704) * g(3); -t725 * t125 + t144 * t126 - t705 * t209 + t238 * t210 - t511 * t268 + t323 * t269 + t13 + t149 + t46 + (t144 * t20 - t21 * t725 + t52) * m(7) + (t238 * t79 - t705 * t80 + t121) * m(6) + (t194 * t323 - t195 * t511 + t168) * m(5) + (t450 * g(3) - t449 * t673) * (m(5) - t701); (mrSges(6,2) * t343 - t342 * t690 - t557) * g(1) + (-mrSges(6,2) * t341 + t340 * t690 - t558) * g(2) - (mrSges(7,1) * t138 + Ifges(7,4) * t642 + Ifges(7,2) * t644 + Ifges(7,6) * t630 - t618 + t653) * t144 + (-mrSges(7,2) * t138 + Ifges(7,1) * t642 + Ifges(7,4) * t644 + Ifges(7,5) * t630 + t619 + t651) * t725 + (t466 * t586 - t23) * t125 + (t466 * t28 + t461 * t29) * pkin(5) + (m(7) * t609 + mrSges(6,1) * t446 + mrSges(6,2) * t447 - t502) * t602 + (t593 + t210) * t80 + (t594 - t209) * t79 - t78 * t612 - m(7) * (t138 * t612 + t20 * t22 + t21 * t23) + (-Ifges(6,2) * t238 + t131 + t233) * t640 + t666 + (-t461 * t586 - t22) * t126 + t522 + t486 - t223 * (mrSges(6,1) * t238 + mrSges(6,2) * t705) + (Ifges(6,5) * t705 - Ifges(6,6) * t238) * t628 + t130 * t637 + (Ifges(6,1) * t705 - t588) * t638 + (t3 * t461 + t4 * t466 + (-t20 * t461 + t21 * t466) * qJD(6)) * t660; -t138 * (mrSges(7,1) * t144 + mrSges(7,2) * t725) + (Ifges(7,1) * t725 - t587) * t642 + t72 * t641 + (Ifges(7,5) * t725 - Ifges(7,6) * t144) * t630 - t20 * t125 + t21 * t126 - g(1) * t557 - g(2) * t558 - t502 * t602 + (t144 * t21 + t20 * t725) * mrSges(7,3) + t486 + (-Ifges(7,2) * t144 + t139 + t73) * t644;];
tau  = t1;
