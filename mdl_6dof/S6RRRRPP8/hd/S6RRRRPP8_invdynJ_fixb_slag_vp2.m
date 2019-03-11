% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:20
% EndTime: 2019-03-09 21:31:50
% DurationCPUTime: 57.21s
% Computational Cost: add. (15534->1042), mult. (37754->1312), div. (0->0), fcn. (29551->10), ass. (0->447)
t706 = Ifges(7,4) + Ifges(6,5);
t707 = Ifges(6,4) + Ifges(5,5);
t655 = -Ifges(7,5) + t707;
t563 = cos(pkin(6));
t481 = t563 * qJD(1);
t360 = t481 + qJD(2);
t377 = sin(qJ(3));
t380 = cos(qJ(3));
t378 = sin(qJ(2));
t375 = sin(pkin(6));
t539 = qJD(1) * t375;
t508 = t378 * t539;
t261 = t380 * t360 - t377 * t508;
t381 = cos(qJ(2));
t526 = qJD(1) * qJD(2);
t305 = (qJDD(1) * t378 + t381 * t526) * t375;
t478 = t563 * qJDD(1);
t359 = t478 + qJDD(2);
t154 = t261 * qJD(3) + t305 * t380 + t359 * t377;
t304 = (-qJDD(1) * t381 + t378 * t526) * t375;
t287 = qJDD(3) + t304;
t376 = sin(qJ(4));
t379 = cos(qJ(4));
t262 = t360 * t377 + t380 * t508;
t504 = t381 * t539;
t425 = -qJD(3) + t504;
t202 = t376 * t262 + t379 * t425;
t535 = qJD(4) * t202;
t74 = t379 * t154 + t376 * t287 - t535;
t636 = t74 / 0.2e1;
t203 = t379 * t262 - t376 * t425;
t75 = qJD(4) * t203 + t376 * t154 - t379 * t287;
t635 = -t75 / 0.2e1;
t155 = -qJD(3) * t262 - t305 * t377 + t380 * t359;
t149 = qJDD(4) - t155;
t626 = t149 / 0.2e1;
t705 = Ifges(7,2) + Ifges(6,3);
t703 = -Ifges(7,6) + Ifges(6,6);
t656 = -Ifges(5,4) + t706;
t657 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t741 = t657 * t636;
t470 = t377 * t504;
t537 = qJD(3) * t377;
t666 = t470 - t537;
t608 = t287 / 0.2e1;
t624 = t155 / 0.2e1;
t627 = -t149 / 0.2e1;
t634 = t75 / 0.2e1;
t702 = Ifges(5,3) + Ifges(6,2);
t469 = pkin(1) * t481;
t297 = -pkin(8) * t508 + t381 * t469;
t238 = -t360 * pkin(2) - t297;
t121 = -t261 * pkin(3) - t262 * pkin(10) + t238;
t511 = pkin(1) * t563;
t368 = t378 * t511;
t552 = t375 * t381;
t541 = pkin(8) * t552 + t368;
t300 = t541 * qJD(1);
t239 = t360 * pkin(9) + t300;
t250 = (-pkin(2) * t381 - pkin(9) * t378 - pkin(1)) * t539;
t151 = t380 * t239 + t377 * t250;
t125 = -pkin(10) * t425 + t151;
t722 = -pkin(8) * t375 * t526 + pkin(1) * t478;
t525 = qJDD(1) * t375;
t723 = pkin(8) * t525 + qJD(2) * t469;
t204 = t378 * t722 + t381 * t723;
t180 = pkin(9) * t359 + t204;
t518 = pkin(1) * t525;
t188 = pkin(2) * t304 - pkin(9) * t305 - t518;
t536 = qJD(3) * t380;
t52 = t380 * t180 + t377 * t188 - t239 * t537 + t250 * t536;
t44 = pkin(10) * t287 + t52;
t205 = -t378 * t723 + t381 * t722;
t181 = -t359 * pkin(2) - t205;
t50 = -t155 * pkin(3) - t154 * pkin(10) + t181;
t531 = qJD(4) * t379;
t533 = qJD(4) * t376;
t7 = -t121 * t533 - t125 * t531 - t376 * t44 + t379 * t50;
t405 = qJDD(5) - t7;
t630 = pkin(4) + pkin(5);
t1 = -qJ(6) * t74 - qJD(6) * t203 - t149 * t630 + t405;
t253 = qJD(4) - t261;
t6 = t121 * t531 - t125 * t533 + t376 * t50 + t379 * t44;
t4 = t149 * qJ(5) + t253 * qJD(5) + t6;
t2 = qJ(6) * t75 + qJD(6) * t202 + t4;
t5 = -pkin(4) * t149 + t405;
t713 = t7 * mrSges(5,1) - t5 * mrSges(6,1) - t1 * mrSges(7,1) - t6 * mrSges(5,2) + t2 * mrSges(7,2) + t4 * mrSges(6,3);
t740 = -Ifges(4,2) * t624 - Ifges(4,6) * t608 + Ifges(5,6) * t635 - Ifges(7,3) * t627 + t626 * t702 + t634 * t703 + t713;
t739 = Ifges(5,6) - Ifges(6,6);
t416 = (pkin(2) * t378 - pkin(9) * t381) * t375;
t298 = qJD(1) * t416;
t191 = t380 * t297 + t377 * t298;
t167 = pkin(10) * t508 + t191;
t458 = pkin(3) * t377 - pkin(10) * t380;
t201 = (t368 + (pkin(8) + t458) * t552) * qJD(1);
t598 = pkin(3) * t380;
t337 = -pkin(10) * t377 - pkin(2) - t598;
t547 = t379 * t380;
t370 = pkin(9) * t547;
t738 = qJD(4) * t370 - t376 * t167 + t201 * t379 + t337 * t533;
t737 = t706 * t203;
t561 = qJ(5) * t376;
t597 = pkin(4) * t379;
t427 = t561 + t597;
t736 = t706 * t202;
t735 = -t706 * t74 / 0.2e1 + t703 * t627 - mrSges(7,3) * t2 + Ifges(5,4) * t636 + Ifges(5,6) * t626 + (t705 + Ifges(5,2)) * t635;
t55 = t379 * t121 - t376 * t125;
t32 = qJ(6) * t203 + t55;
t727 = qJD(5) - t32;
t29 = -t253 * t630 + t727;
t244 = t253 * qJ(5);
t56 = t376 * t121 + t379 * t125;
t33 = qJ(6) * t202 + t56;
t30 = t244 + t33;
t726 = qJD(5) - t55;
t46 = -pkin(4) * t253 + t726;
t47 = t244 + t56;
t695 = Ifges(4,6) * t425;
t734 = -t238 * mrSges(4,1) - t55 * mrSges(5,1) + t46 * mrSges(6,1) + t29 * mrSges(7,1) + t56 * mrSges(5,2) - t30 * mrSges(7,2) - t47 * mrSges(6,3) - t695 / 0.2e1;
t53 = -t377 * t180 + t380 * t188 - t239 * t536 - t250 * t537;
t411 = pkin(3) * t287 + t53;
t385 = qJ(5) * t74 + qJD(5) * t203 + t411;
t3 = -t630 * t75 + qJDD(6) + t385;
t733 = t3 * mrSges(7,1);
t732 = t626 * t655 + t634 * t656 + t741;
t697 = -t202 * t739 + t203 * t707 + t253 * t702;
t731 = t697 / 0.2e1;
t698 = t202 * t705 + t253 * t703 + t737;
t730 = Ifges(4,2) * t261;
t451 = mrSges(6,1) * t379 + mrSges(6,3) * t376;
t453 = mrSges(5,1) * t379 - mrSges(5,2) * t376;
t519 = m(7) * pkin(5) + mrSges(7,1);
t567 = t376 * mrSges(7,2);
t649 = t379 * t519 + t451 + t453 + t567;
t729 = mrSges(4,1) + t649;
t581 = Ifges(4,4) * t262;
t690 = t203 * Ifges(7,5) + t202 * Ifges(7,6) - t253 * Ifges(7,3) + t581 - t695 + t730;
t101 = t379 * t167 + t376 * t201;
t728 = -qJ(5) * t666 - qJD(5) * t380 - t101;
t329 = t458 * qJD(3);
t725 = t329 * t379 - t738;
t724 = qJD(5) * t376 + t151;
t721 = t706 * t379;
t720 = t706 * t376;
t578 = Ifges(5,4) * t203;
t91 = -Ifges(5,2) * t202 + Ifges(5,6) * t253 + t578;
t633 = -t91 / 0.2e1;
t719 = t698 / 0.2e1 + t633;
t200 = Ifges(5,4) * t202;
t653 = t657 * t203 + t253 * t655 - t200 + t736;
t601 = cos(qJ(1));
t460 = t563 * t601;
t600 = sin(qJ(1));
t317 = t378 * t460 + t381 * t600;
t510 = t375 * t601;
t228 = t317 * t380 - t377 * t510;
t316 = t378 * t600 - t381 * t460;
t168 = t228 * t376 - t316 * t379;
t718 = t228 * t379 + t316 * t376;
t625 = t154 / 0.2e1;
t716 = Ifges(4,1) * t625 + Ifges(4,5) * t608;
t715 = mrSges(6,2) * t5 - mrSges(5,3) * t7 - mrSges(7,3) * t1 + t732;
t714 = -mrSges(6,2) * t4 - mrSges(5,3) * t6 - t735;
t631 = m(6) + m(7);
t710 = mrSges(5,2) - mrSges(6,3);
t709 = -mrSges(6,2) - mrSges(5,3);
t708 = mrSges(4,3) - mrSges(3,2);
t700 = t149 * t702 + t707 * t74 - t739 * t75;
t14 = Ifges(7,5) * t74 + Ifges(7,6) * t75 - Ifges(7,3) * t149;
t699 = t154 * Ifges(4,4) + t155 * Ifges(4,2) + t287 * Ifges(4,6) + t14;
t696 = Ifges(4,5) * t425;
t694 = Ifges(4,3) * t425;
t546 = t380 * t381;
t264 = (t376 * t378 + t379 * t546) * t375;
t249 = qJD(1) * t264;
t512 = -pkin(9) * t376 - pkin(4);
t527 = qJD(6) * t379;
t693 = qJ(6) * t249 + t470 * t630 + (-qJ(6) * t536 - t329) * t379 + (qJ(6) * t533 - t527 + (-pkin(5) + t512) * qJD(3)) * t377 + t738;
t515 = t376 * t552;
t474 = t380 * t515;
t248 = qJD(1) * t474 - t379 * t508;
t544 = t376 * t329 + t337 * t531;
t549 = t377 * t379;
t692 = -qJ(6) * t248 + (-pkin(9) * qJD(3) + qJ(6) * qJD(4)) * t549 + (qJD(6) * t377 + (-pkin(9) * qJD(4) + qJ(6) * qJD(3)) * t380) * t376 + t544 + t728;
t190 = -t377 * t297 + t380 * t298;
t406 = pkin(3) * t508 + t190;
t388 = qJ(5) * t249 + t406;
t560 = qJ(5) * t379;
t415 = -t376 * t630 + t560;
t403 = -pkin(9) + t415;
t529 = qJD(5) * t379;
t646 = -t379 * t630 - t561;
t691 = t248 * t630 + (qJD(4) * t646 + t529) * t377 + t403 * t536 - t388;
t183 = (-t379 * t537 - t380 * t533) * pkin(9) + t544;
t689 = t183 + t728;
t688 = pkin(4) * t470 + t512 * t537 - t725;
t426 = pkin(4) * t376 - t560;
t417 = pkin(9) + t426;
t687 = -pkin(4) * t248 + (qJD(4) * t427 - t529) * t377 + t417 * t536 + t388;
t686 = t253 * t415 + t724;
t559 = t261 * t376;
t586 = pkin(10) - qJ(6);
t150 = -t377 * t239 + t380 * t250;
t185 = pkin(3) * t262 - pkin(10) * t261;
t82 = t379 * t150 + t376 * t185;
t62 = t262 * qJ(5) + t82;
t685 = -qJ(6) * t559 - t533 * t586 - t527 - t62;
t684 = t253 * t426 - t724;
t133 = t376 * t150;
t344 = t586 * t379;
t683 = qJD(4) * t344 - qJD(6) * t376 - t133 - (-qJ(6) * t261 - t185) * t379 + t630 * t262;
t454 = mrSges(4,1) * t380 - mrSges(4,2) * t377;
t682 = t454 + mrSges(3,1);
t462 = m(7) * t586 - mrSges(7,3);
t681 = -t462 + mrSges(4,2);
t680 = t425 * (Ifges(4,5) * t380 - Ifges(4,6) * t377);
t679 = t183 - t101;
t523 = pkin(9) * t537;
t678 = t376 * t523 + t725;
t479 = -t317 * t377 - t380 * t510;
t677 = t427 * t479;
t459 = t563 * t600;
t319 = -t378 * t459 + t381 * t601;
t509 = t375 * t600;
t231 = t319 * t377 - t380 * t509;
t676 = t427 * t231;
t505 = t376 * t536;
t675 = -t377 * t531 + t248 - t505;
t532 = qJD(4) * t377;
t674 = t376 * t532 - t379 * t536 + t249;
t553 = t375 * t378;
t314 = t377 * t553 - t380 * t563;
t303 = t314 * qJ(5);
t673 = -t376 * t303 - t314 * t597;
t322 = -pkin(8) * t553 + t381 * t511;
t670 = -t379 * t705 + t720;
t669 = t376 * t705 + t721;
t668 = t376 * t707 + t379 * t739;
t667 = -t376 * t739 + t379 * t707;
t456 = mrSges(5,1) + mrSges(6,1) + t519;
t665 = pkin(4) * t631 + t456;
t558 = t261 * t379;
t664 = t531 - t558;
t663 = t533 - t559;
t662 = -t377 * t53 + t380 * t52;
t661 = -t376 * t7 + t379 * t6;
t660 = t376 * t5 + t379 * t4;
t583 = Ifges(3,4) * t378;
t640 = t375 ^ 2;
t659 = (pkin(1) * (mrSges(3,1) * t378 + mrSges(3,2) * t381) - t378 * (Ifges(3,1) * t381 - t583) / 0.2e1) * t640;
t658 = -mrSges(7,2) + t710;
t124 = pkin(3) * t425 - t150;
t386 = -t203 * qJ(5) + t124;
t31 = -t202 * t630 + qJD(6) - t386;
t449 = -mrSges(7,1) * t376 + mrSges(7,2) * t379;
t450 = mrSges(6,1) * t376 - mrSges(6,3) * t379;
t452 = mrSges(5,1) * t376 + mrSges(5,2) * t379;
t54 = t202 * pkin(4) + t386;
t652 = t124 * t452 + t31 * t449 + t54 * t450;
t576 = Ifges(5,4) * t379;
t651 = t376 * t657 + t576 - t721;
t577 = Ifges(5,4) * t376;
t650 = t379 * t657 - t577 + t720;
t647 = t681 + t709;
t477 = mrSges(3,3) * t508;
t645 = -m(4) * t238 + mrSges(3,1) * t360 + mrSges(4,1) * t261 - mrSges(4,2) * t262 - t477;
t643 = m(7) * qJ(6) + mrSges(7,3) + t709;
t232 = t319 * t380 + t377 * t509;
t315 = t377 * t563 + t380 * t553;
t642 = -g(1) * t232 - g(2) * t228 - g(3) * t315;
t637 = Ifges(4,4) * t624 + t716;
t632 = t91 / 0.2e1;
t623 = -t202 / 0.2e1;
t622 = t202 / 0.2e1;
t621 = -t203 / 0.2e1;
t620 = t203 / 0.2e1;
t613 = -t253 / 0.2e1;
t612 = t253 / 0.2e1;
t611 = -t261 / 0.2e1;
t610 = -t262 / 0.2e1;
t609 = t262 / 0.2e1;
t599 = pkin(1) * t375;
t595 = pkin(10) * t231;
t592 = t479 * pkin(10);
t585 = mrSges(5,3) * t202;
t584 = mrSges(5,3) * t203;
t582 = Ifges(3,4) * t381;
t580 = Ifges(4,4) * t377;
t579 = Ifges(4,4) * t380;
t569 = t261 * mrSges(4,3);
t568 = t262 * mrSges(4,3);
t566 = t377 * t411;
t562 = qJ(5) * t202;
t556 = t316 * t377;
t318 = t378 * t601 + t381 * t459;
t554 = t318 * t377;
t551 = t376 * t377;
t550 = t376 * t380;
t548 = t377 * t381;
t126 = -mrSges(6,2) * t202 + mrSges(6,3) * t253;
t127 = mrSges(7,2) * t253 + mrSges(7,3) * t202;
t545 = t126 + t127;
t282 = -pkin(2) * t563 - t322;
t306 = t314 * pkin(3);
t482 = t315 * pkin(10) - t306;
t158 = t282 - t482;
t283 = pkin(9) * t563 + t541;
t542 = pkin(2) * t552 + pkin(9) * t553;
t284 = -t542 - t599;
t187 = t380 * t283 + t377 * t284;
t160 = -pkin(10) * t552 + t187;
t79 = t376 * t158 + t379 * t160;
t186 = -t377 * t283 + t380 * t284;
t275 = t376 * t337 + t370;
t540 = t601 * pkin(1) + pkin(8) * t509;
t538 = qJD(2) * t381;
t364 = pkin(3) * t552;
t522 = pkin(9) * t536;
t521 = pkin(10) * t533;
t520 = pkin(10) * t531;
t516 = t375 * t548;
t514 = Ifges(4,5) * t154 + Ifges(4,6) * t155 + Ifges(4,3) * t287;
t65 = t303 + t79;
t159 = t364 - t186;
t513 = Ifges(3,5) * t305 - Ifges(3,6) * t304 + Ifges(3,3) * t359;
t507 = qJD(2) * t553;
t506 = t375 * t538;
t501 = t553 / 0.2e1;
t27 = -t75 * mrSges(7,1) + t74 * mrSges(7,2);
t497 = t539 / 0.2e1;
t493 = t536 / 0.2e1;
t488 = -t532 / 0.2e1;
t487 = t531 / 0.2e1;
t37 = -t149 * mrSges(6,1) + t74 * mrSges(6,2);
t35 = -t149 * mrSges(7,1) - t74 * mrSges(7,3);
t486 = -t316 * pkin(2) + pkin(9) * t317;
t485 = -t318 * pkin(2) + pkin(9) * t319;
t218 = t479 * pkin(3);
t484 = pkin(10) * t228 + t218;
t220 = t231 * pkin(3);
t483 = pkin(10) * t232 - t220;
t81 = t185 * t379 - t133;
t78 = t158 * t379 - t376 * t160;
t369 = pkin(9) * t550;
t274 = t337 * t379 - t369;
t476 = mrSges(3,3) * t504;
t299 = qJD(2) * t416;
t301 = t322 * qJD(2);
t103 = -t283 * t536 - t284 * t537 + t380 * t299 - t377 * t301;
t473 = pkin(10) * t516 + t380 * t364 + t542;
t461 = -pkin(1) * t600 + pkin(8) * t510;
t241 = -qJ(5) * t380 + t275;
t455 = mrSges(4,1) * t314 + mrSges(4,2) * t315;
t448 = Ifges(4,1) * t380 - t580;
t441 = -Ifges(4,2) * t377 + t579;
t440 = -Ifges(5,2) * t376 + t576;
t439 = Ifges(5,2) * t379 + t577;
t429 = Ifges(7,5) * t379 + Ifges(7,6) * t376;
t428 = Ifges(7,5) * t376 - Ifges(7,6) * t379;
t422 = -qJ(5) * t631 + t710;
t421 = -pkin(10) * t556 - t316 * t598 + t486;
t420 = -pkin(10) * t554 - t318 * t598 + t485;
t419 = t319 * pkin(2) + pkin(9) * t318 + t540;
t223 = qJD(3) * t315 + t377 * t506;
t224 = -qJD(3) * t314 + t380 * t506;
t302 = t541 * qJD(2);
t113 = t223 * pkin(3) - t224 * pkin(10) + t302;
t102 = -t283 * t537 + t284 * t536 + t377 * t299 + t380 * t301;
t98 = pkin(10) * t507 + t102;
t24 = t113 * t379 - t158 * t533 - t160 * t531 - t376 * t98;
t414 = -mrSges(7,2) + t422;
t225 = t315 * t376 + t379 * t552;
t410 = t232 * pkin(3) + t419;
t408 = t238 * (mrSges(4,1) * t377 + mrSges(4,2) * t380);
t407 = (Ifges(3,2) * t381 + t583) * t375;
t23 = t376 * t113 + t158 * t531 - t160 * t533 + t379 * t98;
t217 = t225 * pkin(4);
t226 = t315 * t379 - t515;
t80 = -qJ(5) * t226 + t159 + t217;
t172 = t232 * t376 - t318 * t379;
t402 = -g(1) * t172 - g(2) * t168 - g(3) * t225;
t400 = t360 * t375 * (Ifges(3,5) * t381 - Ifges(3,6) * t378);
t394 = -t317 * pkin(2) - t316 * pkin(9) + t461;
t392 = pkin(3) * t507 + t103;
t12 = t223 * qJ(5) + t314 * qJD(5) + t23;
t391 = -pkin(3) * t228 + t394;
t173 = t232 * t379 + t318 * t376;
t387 = t173 * pkin(4) + qJ(5) * t172 + t410;
t384 = -pkin(4) * t718 - qJ(5) * t168 + t391;
t119 = -qJD(4) * t225 + t224 * t379 + t376 * t507;
t383 = qJ(5) * t119 + qJD(5) * t226 + t392;
t373 = t380 * pkin(4);
t356 = Ifges(3,4) * t504;
t343 = t586 * t376;
t333 = -pkin(3) - t427;
t326 = pkin(3) - t646;
t320 = (-mrSges(3,1) * t381 + mrSges(3,2) * t378) * t375;
t295 = -mrSges(3,2) * t360 + t476;
t288 = t417 * t377;
t263 = -t379 * t553 + t474;
t252 = Ifges(4,4) * t261;
t242 = -t274 + t373;
t240 = t403 * t377;
t236 = Ifges(3,1) * t508 + Ifges(3,5) * t360 + t356;
t235 = t360 * Ifges(3,6) + qJD(1) * t407;
t216 = t226 * mrSges(7,2);
t213 = qJ(6) * t551 + t241;
t208 = -mrSges(4,1) * t425 - t568;
t207 = mrSges(4,2) * t425 + t569;
t206 = pkin(5) * t380 + t369 + t373 + (-qJ(6) * t377 - t337) * t379;
t197 = -t318 * t547 + t319 * t376;
t196 = -t318 * t550 - t319 * t379;
t195 = -t316 * t547 + t317 * t376;
t194 = -t316 * t550 - t317 * t379;
t143 = Ifges(4,1) * t262 + t252 - t696;
t141 = Ifges(4,5) * t262 + Ifges(4,6) * t261 - t694;
t131 = -mrSges(6,1) * t253 + mrSges(6,2) * t203;
t130 = mrSges(5,1) * t253 - t584;
t129 = -mrSges(7,1) * t253 - mrSges(7,3) * t203;
t128 = -mrSges(5,2) * t253 - t585;
t118 = -qJD(4) * t515 + t224 * t376 + t315 * t531 - t379 * t507;
t112 = -mrSges(4,2) * t287 + mrSges(4,3) * t155;
t111 = mrSges(4,1) * t287 - mrSges(4,3) * t154;
t110 = mrSges(5,1) * t202 + mrSges(5,2) * t203;
t109 = -mrSges(7,1) * t202 + mrSges(7,2) * t203;
t108 = mrSges(6,1) * t202 - mrSges(6,3) * t203;
t107 = pkin(4) * t203 + t562;
t77 = -mrSges(4,1) * t155 + mrSges(4,2) * t154;
t76 = -t203 * t630 - t562;
t66 = -pkin(4) * t314 - t78;
t64 = -pkin(4) * t262 - t81;
t58 = -pkin(5) * t225 - t80;
t49 = qJ(6) * t225 + t65;
t43 = -qJ(6) * t226 - t314 * t630 - t78;
t39 = -mrSges(5,2) * t149 - mrSges(5,3) * t75;
t38 = mrSges(7,2) * t149 + mrSges(7,3) * t75;
t36 = mrSges(5,1) * t149 - mrSges(5,3) * t74;
t34 = -mrSges(6,2) * t75 + mrSges(6,3) * t149;
t28 = mrSges(5,1) * t75 + mrSges(5,2) * t74;
t26 = mrSges(6,1) * t75 - mrSges(6,3) * t74;
t25 = pkin(4) * t118 - t383;
t13 = -pkin(4) * t223 - t24;
t11 = -t118 * t630 + t383;
t10 = qJ(6) * t118 + qJD(6) * t225 + t12;
t9 = pkin(4) * t75 - t385;
t8 = -qJ(6) * t119 - qJD(6) * t226 - t223 * t630 - t24;
t15 = [(-t690 / 0.2e1 + t731 - t151 * mrSges(4,3) - Ifges(7,3) * t613 - Ifges(4,4) * t609 + t655 * t620 + t703 * t622 + t702 * t612 + Ifges(5,6) * t623 - t730 / 0.2e1 - t734) * t223 + m(4) * (t102 * t151 + t103 * t150 + t181 * t282 + t186 * t53 + t187 * t52) + (Ifges(4,5) * t315 - Ifges(4,3) * t552) * t608 - t425 * (Ifges(4,5) * t224 + Ifges(4,3) * t507) / 0.2e1 + (t204 * t552 - t205 * t553 - t297 * t506 - t300 * t507 - t304 * t541 - t305 * t322) * mrSges(3,3) + (Ifges(3,4) * t305 - Ifges(3,2) * t304 + Ifges(3,6) * t359) * t552 / 0.2e1 + (t205 * t563 - t304 * t599 + t322 * t359) * mrSges(3,1) + (Ifges(3,1) * t305 - Ifges(3,4) * t304 + Ifges(3,5) * t359) * t501 + (-m(5) * (t391 + t592) - m(6) * (t384 + t592) - m(4) * t394 + t228 * mrSges(4,1) - m(3) * t461 + t317 * mrSges(3,1) - mrSges(3,3) * t510 - m(7) * t384 + mrSges(2,1) * t600 + mrSges(2,2) * t601 + t708 * t316 + t456 * t718 - t658 * t168 + t647 * t479) * g(1) + (-m(3) * t540 - t319 * mrSges(3,1) - mrSges(3,3) * t509 - m(5) * (t410 + t595) - m(6) * (t387 + t595) - m(7) * t387 - mrSges(2,1) * t601 + mrSges(2,2) * t600 - m(4) * t419 - t232 * mrSges(4,1) - t708 * t318 - t456 * t173 + t658 * t172 + t647 * t231) * g(2) + m(7) * (t1 * t43 + t10 * t30 + t11 * t31 + t2 * t49 + t29 * t8 + t3 * t58) + m(6) * (t12 * t47 + t13 * t46 + t25 * t54 + t4 * t65 + t5 * t66 + t80 * t9) + t53 * (-mrSges(4,1) * t552 - mrSges(4,3) * t315) + t305 * (Ifges(3,5) * t563 + (Ifges(3,1) * t378 + t582) * t375) / 0.2e1 + (-t411 * mrSges(5,1) + t9 * mrSges(6,1) - Ifges(5,2) * t635 + Ifges(7,6) * t627 - t626 * t739 + t705 * t634 + t656 * t636 + t714 - t733) * t225 + (t124 * mrSges(5,1) + t54 * mrSges(6,1) - t31 * mrSges(7,1) - t47 * mrSges(6,2) - t56 * mrSges(5,3) + t30 * mrSges(7,3) - Ifges(5,2) * t623 + Ifges(7,6) * t613 - t612 * t739 + t656 * t620 + t705 * t622 + t719) * t118 + (-t151 * t507 + t224 * t238 + t52 * t552) * mrSges(4,2) + (t640 * qJD(1) * (-Ifges(3,2) * t378 + t582) + t375 * t236) * t538 / 0.2e1 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t640 + t204 * t541 + t205 * t322 + t300 * t301) + (-t204 * t563 - t305 * t599 - t359 * t541) * mrSges(3,2) + (Ifges(7,5) * t613 + t653 / 0.2e1 + t657 * t620 + t706 * t622 + t707 * t612 + Ifges(5,4) * t623 - t29 * mrSges(7,3) + t46 * mrSges(6,2) - t55 * mrSges(5,3) + t124 * mrSges(5,2) - t54 * mrSges(6,3) + t31 * mrSges(7,2)) * t119 + t3 * t216 - t659 * t526 - t514 * t552 / 0.2e1 + (t400 / 0.2e1 + t141 * t501) * qJD(2) - t235 * t507 / 0.2e1 + t150 * (mrSges(4,1) * t507 - mrSges(4,3) * t224) + t315 * t637 - t320 * t518 - t304 * (Ifges(3,6) * t563 + t407) / 0.2e1 + t563 * t513 / 0.2e1 + t359 * (Ifges(3,3) * t563 + (Ifges(3,5) * t378 + Ifges(3,6) * t381) * t375) / 0.2e1 + t301 * t295 + t282 * t77 + Ifges(2,3) * qJDD(1) + t224 * t143 / 0.2e1 + t103 * t208 + t102 * t207 + t187 * t112 + t186 * t111 + t159 * t28 + t13 * t131 + t8 * t129 + t24 * t130 + t12 * t126 + t10 * t127 + t23 * t128 + t25 * t108 + t11 * t109 + t78 * t36 + t79 * t39 + t80 * t26 + t65 * t34 + t66 * t37 + t58 * t27 + t49 * t38 + t43 * t35 + t181 * t455 + (Ifges(4,1) * t224 + Ifges(4,5) * t507) * t609 + (Ifges(4,1) * t315 - Ifges(4,5) * t552) * t625 + t261 * (Ifges(4,4) * t224 + Ifges(4,6) * t507) / 0.2e1 + (Ifges(4,4) * t315 - Ifges(4,6) * t552) * t624 + (-t52 * mrSges(4,3) + t655 * t636 - t699 / 0.2e1 + t700 / 0.2e1 - Ifges(4,4) * t625 + t740) * t314 + (-t411 * mrSges(5,2) - t9 * mrSges(6,3) + Ifges(5,4) * t635 + Ifges(7,5) * t627 + t707 * t626 + t706 * t634 + t715 + t741) * t226 - t392 * t110 + m(5) * (-t124 * t392 - t159 * t411 + t23 * t56 + t24 * t55 + t6 * t79 + t7 * t78) + (-m(3) * t297 - t645) * t302; (-t400 / 0.2e1 + t659 * qJD(1)) * qJD(1) - (t262 * (Ifges(4,5) * t378 + t381 * t448) + t261 * (Ifges(4,6) * t378 + t381 * t441) + t378 * t141 + (-Ifges(3,2) * t508 + t143 * t380 + t377 * t697 + t236 + t356) * t381) * t539 / 0.2e1 + (-m(5) * t473 - m(4) * t542 - (t378 * mrSges(4,3) + t381 * t454) * t375 + t320 - t631 * (t264 * pkin(4) + qJ(5) * t263 + t473) + t643 * t516 - t456 * t264 + t658 * t263) * g(3) + (-m(4) * t486 - m(5) * t421 - t631 * (t195 * pkin(4) + qJ(5) * t194 + t421) - t708 * t317 + t682 * t316 - t643 * t556 - t456 * t195 + t658 * t194) * g(2) + (-m(4) * t485 - m(5) * t420 - t631 * (t197 * pkin(4) + qJ(5) * t196 + t420) - t708 * t319 + t682 * t318 - t643 * t554 - t456 * t197 + t658 * t196) * g(1) + (-t522 - t190) * t208 + t580 * t624 + (-t249 / 0.2e1 + t376 * t488 + t379 * t493) * t653 - t408 * t504 + (-t523 - t191) * t207 + t381 * t497 * t680 + ((-mrSges(4,1) * t150 + mrSges(4,2) * t151) * t539 + (t235 + t694) * t497) * t378 + (-t428 * t532 + (-Ifges(7,3) * t377 + t380 * t429) * qJD(3) + t702 * t470 + t707 * t249 - t739 * t248) * t613 + t579 * t625 + (t261 * t441 + t262 * t448) * qJD(3) / 0.2e1 + ((t539 * t548 - t537) * t151 + (t539 * t546 - t536) * t150 + t662) * mrSges(4,3) + (-pkin(2) * t181 + ((-t150 * t380 - t151 * t377) * qJD(3) + t662) * pkin(9) - t150 * t190 - t151 * t191) * m(4) + (-t651 * t532 + (t377 * t655 + t380 * t650) * qJD(3)) * t620 + (t377 * t650 - t380 * t655) * t636 + (t248 * t656 + t249 * t657 + t470 * t655) * t621 + (-t680 / 0.2e1 + t408) * qJD(3) + (mrSges(7,1) * t675 - mrSges(7,2) * t674) * t31 + (mrSges(5,2) * t666 + mrSges(5,3) * t675) * t56 + (mrSges(6,2) * t675 - mrSges(6,3) * t666) * t47 + (-mrSges(6,1) * t675 + mrSges(6,3) * t674) * t54 + (-mrSges(5,1) * t675 - mrSges(5,2) * t674) * t124 + (-mrSges(7,2) * t666 - mrSges(7,3) * t675) * t30 + t678 * t130 + t679 * t128 + (mrSges(7,1) * t666 + mrSges(7,3) * t674) * t29 + (-mrSges(5,1) * t666 + mrSges(5,3) * t674) * t55 + (mrSges(6,1) * t666 - mrSges(6,2) * t674) * t46 + t699 * t380 / 0.2e1 - t452 * t566 + (-t111 + t28) * pkin(9) * t377 + t248 * t632 + t505 * t633 - t700 * t380 / 0.2e1 + (-t295 + t476) * t297 + (Ifges(7,5) * t249 + Ifges(7,6) * t248 - Ifges(7,3) * t470 - t668 * t532 + (t377 * t702 + t380 * t667) * qJD(3)) * t612 + (Ifges(5,4) * t249 - Ifges(5,2) * t248 + Ifges(5,6) * t470 - t670 * t532 + (t377 * t703 + t380 * t669) * qJD(3)) * t622 + (-t439 * t532 + (Ifges(5,6) * t377 + t380 * t440) * qJD(3) + t703 * t470 + t706 * t249 + t705 * t248) * t623 + t537 * t731 + t288 * t26 + t274 * t36 + t275 * t39 + t241 * t34 + t242 * t37 + t240 * t27 + t698 * (t376 * t493 + t377 * t487 - t248 / 0.2e1) + t213 * t38 + t690 * (t470 / 0.2e1 - t537 / 0.2e1) + t206 * t35 - t204 * mrSges(3,2) + t205 * mrSges(3,1) + (t3 * t449 + t429 * t627 + t440 * t635 + t450 * t9 + t626 * t667 + t634 * t669 + t637 + t716) * t377 - pkin(2) * t77 + t715 * t549 + t714 * t551 - t181 * t454 + t379 * t91 * t488 + (pkin(9) * t112 - t740) * t380 + t143 * t493 + (t522 + t406) * t110 + (t274 * t7 + t275 * t6 + (t124 * t536 - t566) * pkin(9) + t124 * t406 + t679 * t56 + t678 * t55) * m(5) + (t477 + t645) * t300 + t513 + t687 * t108 + t688 * t131 + t689 * t126 + (t241 * t4 + t242 * t5 + t288 * t9 + t46 * t688 + t47 * t689 + t54 * t687) * m(6) + t691 * t109 + t692 * t127 + t693 * t129 + (t1 * t206 + t2 * t213 + t240 * t3 + t29 * t693 + t30 * t692 + t31 * t691) * m(7); (-Ifges(4,2) * t611 + Ifges(5,6) * t622 - Ifges(7,3) * t612 + t702 * t613 + t655 * t621 + t703 * t623 + t734) * t262 + (t733 + t735) * t379 + (-t1 * t376 - t29 * t664 + t30 * t663) * mrSges(7,3) + t3 * t567 + (t667 / 0.2e1 - t429 / 0.2e1) * qJD(4) * t253 + (t669 / 0.2e1 - t440 / 0.2e1) * t535 + (t143 + t252) * t611 + (-m(6) * (t483 - t676) - m(7) * (-t220 - t676) - m(5) * t483 + t681 * t232 + t729 * t231) * g(1) + (t520 - t64) * t131 + (-t698 / 0.2e1 + t632) * t559 + t343 * t35 + t344 * t38 + (Ifges(4,1) * t610 + t429 * t612 + t440 * t622 - t238 * mrSges(4,2) + t696 / 0.2e1 + t669 * t623 + t667 * t613 + t650 * t621 - t652) * t261 + (t620 * t650 + t652) * qJD(4) + (t333 * t9 - t46 * t64 - t47 * t62 + t684 * t54) * m(6) + ((t37 - t36) * t376 + ((-t376 * t56 - t379 * t55) * qJD(4) + t661) * m(5) + (t39 + t34) * t379 + ((-t376 * t47 + t379 * t46) * qJD(4) + t660) * m(6)) * pkin(10) + (-m(5) * t124 - t110 + t208 + t568) * t151 + (-t521 - t62) * t126 + t326 * t27 + t333 * t26 + t668 * t626 + t670 * t634 + (-m(6) * (t482 + t673) - m(7) * (-t306 + t673) - t462 * t315 + t455 - m(5) * t482 + t649 * t314) * g(3) + (-t581 + t697) * t610 + (t46 * t664 - t47 * t663 + t642 + t660) * mrSges(6,2) + (-t55 * t664 - t56 * t663 + t642 + t661) * mrSges(5,3) + t439 * t635 + (-t520 - t81) * t130 + t376 * t732 + (t569 - t207) * t150 + (-t521 - t82) * t128 + t428 * t627 + (-m(7) * (t218 + t677) - m(6) * (t484 + t677) - m(5) * t484 + t681 * t228 - t729 * t479) * g(2) + t719 * t533 - t52 * mrSges(4,2) + t53 * mrSges(4,1) - pkin(3) * t28 - t9 * t451 + (t487 - t558 / 0.2e1) * t653 + (pkin(3) * t411 - t55 * t81 - t56 * t82) * m(5) + t411 * t453 + t514 + t683 * t129 + t684 * t108 + t685 * t127 + (t1 * t343 + t2 * t344 + t29 * t683 + t3 * t326 + t30 * t685 + t31 * t686) * m(7) + t686 * t109 + t690 * t609 + t651 * t636; t700 + (t203 * t705 - t736) * t623 + t713 + t545 * qJD(5) + t91 * t620 + (-Ifges(7,5) * t202 + Ifges(7,6) * t203) * t612 + (t202 * t46 + t203 * t47) * mrSges(6,2) + (-t202 * t657 - t578 + t698 + t737) * t621 + (-t202 * t29 - t203 * t30) * mrSges(7,3) + (-t202 * t707 - t203 * t739) * t613 + (t130 - t131 + t584) * t56 + (-t126 - t128 - t585) * t55 + (t172 * t665 + t173 * t414) * g(1) + (t217 * t631 + t225 * t456 + t226 * t422 - t216) * g(3) + (t38 + t34) * qJ(5) - t14 - t124 * (mrSges(5,1) * t203 - mrSges(5,2) * t202) - t54 * (mrSges(6,1) * t203 + mrSges(6,3) * t202) - t31 * (-mrSges(7,1) * t203 - mrSges(7,2) * t202) + (-pkin(4) * t5 + qJ(5) * t4 - t107 * t54 - t46 * t56 + t47 * t726) * m(6) + (t2 * qJ(5) - t1 * t630 - t29 * t33 + t30 * t727 - t31 * t76) * m(7) + (-Ifges(5,2) * t203 - t200 + t653) * t622 + (t665 * t168 + t414 * t718) * g(2) - t32 * t127 - t33 * t129 - t107 * t108 - t76 * t109 - pkin(4) * t37 - t630 * t35; -t545 * t253 + (t108 - t109) * t203 + t35 + t37 + (-t203 * t31 - t253 * t30 + t1 + t402) * m(7) + (t203 * t54 - t253 * t47 + t402 + t5) * m(6); -t202 * t127 + t203 * t129 + (g(1) * t231 - g(2) * t479 + g(3) * t314 - t202 * t30 + t203 * t29 + t3) * m(7) + t27;];
tau  = t15;
