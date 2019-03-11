% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:46
% EndTime: 2019-03-09 05:05:09
% DurationCPUTime: 12.41s
% Computational Cost: add. (14448->657), mult. (29104->875), div. (0->0), fcn. (26719->8), ass. (0->344)
t464 = sin(qJ(4));
t667 = pkin(8) * t464;
t414 = -t464 * pkin(9) + t667;
t467 = cos(qJ(4));
t415 = (pkin(8) - pkin(9)) * t467;
t463 = sin(qJ(6));
t466 = cos(qJ(6));
t253 = t414 * t466 - t415 * t463;
t510 = t414 * t463 + t415 * t466;
t599 = t464 * t466;
t385 = -t463 * t467 + t599;
t509 = t463 * t464 + t466 * t467;
t587 = -Ifges(7,5) * t509 - Ifges(7,6) * t385;
t761 = -t510 * mrSges(7,1) - t253 * mrSges(7,2) + t587;
t772 = qJD(6) * t761;
t764 = -Ifges(5,6) + Ifges(6,6);
t771 = mrSges(6,3) / 0.2e1;
t465 = sin(qJ(3));
t468 = cos(qJ(3));
t598 = t464 * t468;
t390 = -mrSges(6,2) * t598 + mrSges(6,3) * t465;
t770 = t390 / 0.2e1;
t552 = t598 / 0.2e1;
t651 = Ifges(5,4) * t464;
t404 = Ifges(5,2) * t467 + t651;
t453 = Ifges(6,5) * t464;
t521 = Ifges(6,3) * t467 - t453;
t769 = t404 + t521;
t597 = t465 * t467;
t551 = -t597 / 0.2e1;
t665 = pkin(8) * t468;
t669 = pkin(3) * t465;
t416 = -t665 + t669;
t441 = sin(pkin(10)) * pkin(1) + pkin(7);
t257 = t464 * t416 - t441 * t597;
t220 = t465 * qJ(5) + t257;
t540 = -t441 * t464 - pkin(4);
t607 = t416 * t467;
t222 = t465 * t540 - t607;
t342 = t385 * t468;
t260 = mrSges(7,2) * t465 + mrSges(7,3) * t342;
t702 = pkin(4) + pkin(5);
t394 = -t463 * qJ(5) - t466 * t702;
t395 = t466 * qJ(5) - t463 * t702;
t679 = -t465 / 0.2e1;
t343 = t509 * t468;
t687 = t343 / 0.2e1;
t688 = t342 / 0.2e1;
t733 = Ifges(7,5) * t687 + Ifges(7,6) * t688;
t754 = mrSges(7,1) / 0.2e1;
t172 = (-pkin(9) * t468 - t416) * t467 + (-pkin(5) + t540) * t465;
t190 = pkin(9) * t598 + t220;
t81 = t172 * t466 - t190 * t463;
t82 = t172 * t463 + t190 * t466;
t532 = Ifges(7,3) * t679 - t82 * mrSges(7,2) / 0.2e1 + t81 * t754 + t733;
t592 = t467 * t468;
t637 = t257 * mrSges(5,2);
t600 = t464 * t465;
t256 = t441 * t600 + t607;
t638 = t256 * mrSges(5,1);
t643 = t222 * mrSges(6,1);
t684 = -t395 / 0.2e1;
t436 = mrSges(6,2) * t592;
t629 = t465 * mrSges(6,1);
t389 = t436 - t629;
t685 = t389 / 0.2e1;
t261 = -mrSges(7,1) * t465 - mrSges(7,3) * t343;
t694 = -t261 / 0.2e1;
t709 = -m(7) / 0.2e1;
t711 = -m(6) / 0.2e1;
t741 = Ifges(6,2) + Ifges(5,3);
t742 = Ifges(6,4) + Ifges(5,5);
t768 = -(-pkin(4) * t222 + qJ(5) * t220) * t711 - (t394 * t81 + t395 * t82) * t709 - pkin(4) * t685 + qJ(5) * t770 + t220 * t771 - t643 / 0.2e1 + t638 / 0.2e1 - t637 / 0.2e1 - t394 * t694 - t260 * t684 - t741 * t679 + t742 * t592 / 0.2e1 - t532 + t764 * t552;
t340 = t463 * t597 - t465 * t599;
t301 = Ifges(7,4) * t340;
t341 = t509 * t465;
t161 = t341 * Ifges(7,1) + t468 * Ifges(7,5) - t301;
t174 = t341 * mrSges(7,1) - t340 * mrSges(7,2);
t563 = t702 * t464;
t622 = qJ(5) * t467;
t382 = -t563 + t622;
t214 = (-t441 + t382) * t465;
t601 = t464 * qJ(5);
t718 = -t467 * t702 - t601;
t372 = pkin(3) - t718;
t523 = Ifges(7,2) * t341 + t301;
t636 = t340 * mrSges(7,3);
t258 = -mrSges(7,2) * t468 - t636;
t695 = t258 / 0.2e1;
t223 = mrSges(7,1) * t385 - mrSges(7,2) * t509;
t765 = t223 / 0.2e1;
t767 = -(t161 / 0.4e1 - t523 / 0.4e1) * t509 + t587 * t468 / 0.4e1 + t214 * t765 + t253 * t695 + t372 * t174 / 0.2e1;
t377 = Ifges(7,4) * t385;
t228 = -Ifges(7,2) * t509 + t377;
t376 = Ifges(7,4) * t509;
t229 = Ifges(7,2) * t385 + t376;
t231 = t385 * Ifges(7,1) - t376;
t232 = Ifges(7,1) * t509 + t377;
t766 = -t372 * t223 + (t228 / 0.2e1 + t232 / 0.2e1) * t385 - (t229 / 0.2e1 - t231 / 0.2e1) * t509;
t558 = -cos(pkin(10)) * pkin(1) - pkin(2);
t666 = pkin(8) * t465;
t371 = -pkin(3) * t468 + t558 - t666;
t586 = -t467 * t371 + t441 * t598;
t187 = pkin(9) * t597 - t586;
t459 = t468 * pkin(4);
t143 = pkin(5) * t468 - t187 + t459;
t217 = t464 * t371 + t441 * t592;
t621 = qJ(5) * t468;
t198 = t217 - t621;
t439 = pkin(9) * t600;
t169 = t198 + t439;
t74 = t143 * t463 + t169 * t466;
t188 = t439 + t217;
t83 = -t187 * t463 + t188 * t466;
t763 = t74 - t83;
t581 = t464 ^ 2 + t467 ^ 2;
t756 = (mrSges(6,2) + mrSges(5,3)) * t581;
t758 = -t229 / 0.4e1 + t231 / 0.4e1;
t753 = -t463 / 0.2e1;
t678 = t465 / 0.2e1;
t676 = -t468 / 0.2e1;
t674 = t468 / 0.2e1;
t734 = t340 * t463 + t466 * t341;
t670 = m(7) * (t597 - t734) * t468;
t568 = t670 / 0.2e1;
t456 = Ifges(5,4) * t467;
t729 = -Ifges(5,2) * t464 + t456;
t752 = t729 / 0.2e1;
t751 = m(7) * t382;
t73 = t143 * t466 - t169 * t463;
t84 = t187 * t466 + t188 * t463;
t740 = t84 + t73;
t545 = t223 * t674;
t602 = t463 * t468;
t628 = t466 * mrSges(7,2);
t584 = t602 * t754 + t628 * t674;
t635 = t341 * mrSges(7,3);
t259 = mrSges(7,1) * t468 - t635;
t595 = t466 * t258;
t591 = t259 * t753 + t595 / 0.2e1;
t747 = t584 - t591;
t528 = t463 * mrSges(7,1) + t628;
t746 = qJD(6) * t528;
t512 = t253 * t463 - t466 * t510;
t698 = t510 / 0.2e1;
t744 = -mrSges(5,1) - mrSges(6,1);
t743 = mrSges(5,2) - mrSges(6,3);
t632 = t509 * mrSges(7,3);
t594 = t466 * t340;
t604 = t463 * t341;
t739 = (t604 / 0.2e1 - t594 / 0.2e1) * mrSges(7,3);
t199 = t459 + t586;
t735 = t199 - t586;
t299 = Ifges(7,6) * t341;
t300 = Ifges(7,5) * t340;
t520 = t300 + t299;
t302 = Ifges(7,4) * t341;
t539 = Ifges(7,1) * t340 + t302;
t730 = t467 * Ifges(6,1) + t453;
t332 = -t468 * Ifges(6,4) + t465 * t730;
t409 = t467 * Ifges(5,1) - t651;
t334 = -t468 * Ifges(5,5) + t409 * t465;
t732 = t334 + t332;
t400 = t464 * mrSges(6,1) - t467 * mrSges(6,3);
t668 = pkin(4) * t464;
t399 = -t622 + t668;
t671 = m(6) * t399;
t731 = t400 + t671;
t408 = Ifges(5,1) * t464 + t456;
t574 = mrSges(5,3) * t600;
t624 = t468 * mrSges(5,2);
t504 = -t574 + t624;
t573 = mrSges(5,3) * t597;
t506 = -t468 * mrSges(5,1) - t573;
t728 = t464 * t504 + t467 * t506;
t437 = Ifges(6,5) * t597;
t328 = -t468 * Ifges(6,6) + Ifges(6,3) * t600 + t437;
t727 = -Ifges(6,1) * t600 + t328 + t437;
t159 = -t340 * Ifges(7,2) + t468 * Ifges(7,6) + t302;
t725 = t539 + t159;
t723 = t464 * t764 + t467 * t742;
t511 = -t256 * t464 + t257 * t467;
t513 = t220 * t467 + t222 * t464;
t722 = t74 / 0.2e1 - t83 / 0.2e1;
t503 = -t300 / 0.2e1 - t299 / 0.2e1;
t625 = t467 * Ifges(6,5);
t403 = Ifges(6,3) * t464 + t625;
t721 = Ifges(5,6) * t678 + Ifges(6,6) * t679 + t403 * t676 + t468 * t752;
t451 = t468 * mrSges(6,3);
t571 = mrSges(6,2) * t600;
t391 = -t451 - t571;
t720 = m(6) * t198 + t391;
t691 = -t341 / 0.2e1;
t693 = -t340 / 0.2e1;
t14 = t174 * t674 + (-t341 ^ 2 / 0.2e1 - t340 ^ 2 / 0.2e1) * mrSges(7,3) + t258 * t693 + t259 * t691;
t719 = t14 * qJD(1);
t716 = (-qJD(4) + qJD(6)) * t174;
t715 = t465 ^ 2;
t714 = 2 * qJD(3);
t713 = 0.2e1 * qJD(4);
t712 = m(5) / 0.2e1;
t710 = m(6) / 0.2e1;
t708 = m(7) / 0.2e1;
t706 = mrSges(6,1) / 0.2e1;
t705 = mrSges(7,3) / 0.2e1;
t175 = mrSges(7,1) * t340 + mrSges(7,2) * t341;
t701 = t175 / 0.2e1;
t224 = mrSges(7,1) * t509 + mrSges(7,2) * t385;
t699 = t224 / 0.2e1;
t697 = -t253 / 0.2e1;
t696 = -t510 / 0.2e1;
t692 = t340 / 0.2e1;
t690 = t341 / 0.2e1;
t683 = t400 / 0.2e1;
t401 = t464 * mrSges(5,1) + t467 * mrSges(5,2);
t682 = t401 / 0.2e1;
t518 = pkin(4) * t467 + t601;
t396 = -pkin(3) - t518;
t672 = m(6) * t396;
t664 = t73 * mrSges(7,2);
t663 = t74 * mrSges(7,1);
t660 = t83 * mrSges(7,1);
t659 = t84 * mrSges(7,2);
t656 = mrSges(5,1) * t465;
t655 = mrSges(5,2) * t465;
t634 = t342 * mrSges(7,1);
t633 = t343 * mrSges(7,2);
t631 = t385 * mrSges(7,3);
t627 = t467 * mrSges(6,2);
t530 = t467 * mrSges(5,1) - t464 * mrSges(5,2);
t623 = -t530 - mrSges(4,1);
t87 = -mrSges(7,1) * t395 - t394 * mrSges(7,2);
t620 = qJD(6) * t87;
t618 = t214 * t174;
t612 = t342 * t466;
t611 = t343 * t463;
t610 = t394 * t340;
t609 = t395 * t341;
t608 = t396 * t465;
t606 = t441 * t465;
t603 = t463 * t385;
t596 = t465 * t468;
t593 = t466 * t509;
t529 = t467 * mrSges(6,1) + t464 * mrSges(6,3);
t585 = -t529 - t224;
t583 = t581 * t665;
t580 = qJD(3) * t468;
t579 = qJD(4) * t465;
t569 = mrSges(6,2) * t597;
t566 = mrSges(5,1) / 0.2e1 + t706;
t565 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t564 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t554 = t600 / 0.2e1;
t550 = t597 / 0.2e1;
t368 = t400 * t465;
t544 = t701 - t368 / 0.2e1;
t542 = -t521 / 0.2e1 - t404 / 0.2e1;
t406 = Ifges(6,1) * t464 - t625;
t541 = t406 / 0.2e1 + t408 / 0.2e1;
t533 = (t699 + t529 / 0.2e1) * t467;
t524 = Ifges(6,4) * t464 - Ifges(6,6) * t467;
t522 = Ifges(5,5) * t464 + Ifges(5,6) * t467;
t160 = Ifges(7,4) * t343 + Ifges(7,2) * t342 - Ifges(7,6) * t465;
t162 = Ifges(7,1) * t343 + Ifges(7,4) * t342 - Ifges(7,5) * t465;
t176 = t633 - t634;
t435 = qJ(5) * t592;
t215 = t435 + (-t441 - t563) * t468;
t262 = (t399 + t441) * t465;
t263 = -t435 + (t441 + t668) * t468;
t330 = -t468 * Ifges(5,6) + t465 * t729;
t333 = Ifges(6,4) * t465 + t468 * t730;
t335 = Ifges(5,5) * t465 + t409 * t468;
t369 = t468 * t400;
t370 = t468 * t401;
t387 = -mrSges(5,3) * t598 - t655;
t388 = -mrSges(5,3) * t592 + t656;
t4 = t217 * t387 - t586 * t388 + t199 * t389 + t198 * t390 + t220 * t391 + t263 * t368 + t262 * t369 + t161 * t687 + t159 * t688 + t160 * t693 + t162 * t690 + t82 * t258 + t81 * t259 + t74 * t260 + t73 * t261 + t214 * t176 + t215 * t175 + m(6) * (t198 * t220 + t199 * t222 + t262 * t263) + m(7) * (t214 * t215 + t73 * t81 + t74 * t82) + m(5) * (t217 * t257 - t256 * t586) + (-t638 + t637 + Ifges(4,4) * t468 + t643 + t558 * mrSges(4,2) + (t334 / 0.2e1 + t332 / 0.2e1 - t565 * t468) * t467 + (t328 / 0.2e1 - t330 / 0.2e1 - t564 * t468) * t464 + t733) * t468 + (t441 * t370 + t558 * mrSges(4,1) + Ifges(7,5) * t691 + Ifges(7,6) * t692 - Ifges(4,4) * t465 + (-t257 * mrSges(5,3) + t564 * t465 - t721) * t464 + (t335 / 0.2e1 + t333 / 0.2e1 - t256 * mrSges(5,3) + t222 * mrSges(6,2) + t565 * t465) * t467 + (-Ifges(4,2) + Ifges(4,1) - Ifges(7,3) + (m(5) * t441 + t401) * t441 - t741) * t468) * t465;
t8 = t258 * t687 + t260 * t690 + t259 * t688 + t261 * t693 + ((t387 / 0.2e1 + t770 + t655 / 0.2e1) * t467 + (-t388 / 0.2e1 + t685 + t656 / 0.2e1) * t464 - t544) * t465 + (t176 / 0.2e1 - t370 / 0.2e1 - t369 / 0.2e1 + (t391 / 0.2e1 + t624 / 0.2e1) * t467 + (mrSges(6,2) * t550 + t468 * t566) * t464) * t468 + ((t198 * t467 + t199 * t464 - t263) * t468 + (t262 + t513) * t465) * t710 + (-t214 * t465 + t215 * t468 - t340 * t81 + t341 * t82 + t342 * t73 + t343 * t74) * t708 + ((t217 * t467 - t441 * t468 + t464 * t586) * t468 + (t511 + t606) * t465) * t712;
t517 = t4 * qJD(1) + t8 * qJD(2);
t367 = t518 * t465;
t483 = t735 * t467 + (-t198 + t217) * t464;
t487 = t465 * t718;
t496 = t529 * t465;
t490 = t496 / 0.2e1;
t505 = t468 * mrSges(6,1) + t569;
t492 = t467 * t505;
t497 = t530 * t465;
t10 = t492 * t679 + t391 * t554 + (-t367 * t468 + t465 * t483) * t711 + (t340 * t763 + t740 * t341 + t468 * t487) * t709 + t497 * t674 + t468 * t490 + t728 * t678 + t14 + t715 * t756 / 0.2e1;
t5 = -t175 * t487 - m(7) * (t214 * t487 + t73 * t83 + t74 * t84) - t74 * t635 + t330 * t550 + t198 * t569 + t199 * t571 + t618 + t523 * t692 + t161 * t693 - t715 * t441 * t530 + t520 * t676 - t367 * t368 - t84 * t258 - t83 * t259 + t73 * t636 + t725 * t691 + (-t524 / 0.2e1 - t522 / 0.2e1) * t596 + t732 * t554 + t727 * t551 + (-m(6) * t367 - t496) * t262 - (-t467 * t408 + t769 * t464) * t715 / 0.2e1 + (-m(6) * t199 - t505 + t506 + t573) * t217 - (-t504 - t574 - t720) * t586;
t516 = -t5 * qJD(1) - t10 * qJD(2);
t9 = t73 * t258 - t520 * t674 - t74 * t259 + t618 + (-t74 * mrSges(7,3) - t159 / 0.2e1 - t539 / 0.2e1) * t341 + (t523 / 0.2e1 - t161 / 0.2e1 + t73 * mrSges(7,3)) * t340;
t515 = t9 * qJD(1) + t14 * qJD(2);
t54 = m(7) * (-t340 * t342 + t341 * t343 - t596) + 0.4e1 * (m(6) / 0.4e1 + m(5) / 0.4e1) * (-0.1e1 + t581) * t596;
t514 = t8 * qJD(1) + t54 * qJD(2);
t507 = t10 * qJD(1);
t501 = t634 / 0.2e1 - t633 / 0.2e1;
t498 = m(7) * t512;
t24 = -t259 * t602 + (t595 - m(7) * (t463 * t73 - t466 * t74) + t720) * t468 + (m(6) * t262 - m(7) * t214 - t175 + t368) * t597;
t495 = -t24 * qJD(1) + qJD(2) * t568;
t491 = t512 * t708;
t13 = -t399 * t529 - pkin(3) * t401 + t382 * t224 + t372 * t751 + t731 * t396 + (-t403 / 0.2e1 + t752 + t541) * t467 + (t730 / 0.2e1 + t409 / 0.2e1 + t542) * t464 + t766;
t470 = -t767 - t391 * t667 / 0.2e1 - t367 * t529 / 0.2e1 - (t408 + t406 + t729) * t600 / 0.4e1 + t769 * t551 + (-t408 + t403) * t600 / 0.4e1 + (t232 + t228) * t341 / 0.4e1 + (t635 + t259) * t698 + (t382 * t214 - t253 * t763 + t372 * t487 + t740 * t510) * t708 + t758 * t340 - t666 * t756 / 0.2e1 + ((t217 / 0.2e1 - t198 / 0.2e1) * mrSges(6,2) + t727 / 0.4e1 - t330 / 0.4e1) * t464 + (t492 / 0.2e1 + t483 * t710 - t728 / 0.2e1) * pkin(8) + t636 * t697 + t487 * t699 + t382 * t701 - t723 * t468 / 0.4e1 + t725 * t385 / 0.4e1 + t722 * t631 + t396 * t490 + t606 * t682 + t262 * t683 - t740 * t632 / 0.2e1 + t399 * t368 / 0.2e1 + t732 * t467 / 0.4e1 + t735 * t627 / 0.2e1 + (t262 * t399 + t367 * t396) * t710 - pkin(3) * t497 / 0.2e1 + (t409 + t730) * t597 / 0.4e1;
t2 = t470 - t768;
t475 = (-pkin(4) * t598 + t435) * t710 + (t342 * t394 + t343 * t395) * t708 - t501;
t478 = t468 * t751;
t481 = (-mrSges(5,2) / 0.2e1 + t771) * t467 - t566 * t464;
t21 = (t671 / 0.2e1 + t765 + t682 + t683 + t481) * t468 - t478 / 0.2e1 + t475;
t489 = t2 * qJD(1) - t21 * qJD(2) + t13 * qJD(3);
t32 = t545 - t501;
t472 = (-t539 / 0.4e1 - t159 / 0.4e1) * t385 + (t253 * t705 - t758) * t340 + (-t228 / 0.4e1 - t232 / 0.4e1 + mrSges(7,3) * t696) * t341 + t259 * t696 + t767;
t6 = t472 - t532;
t488 = t6 * qJD(1) + t32 * qJD(2) - qJD(3) * t766;
t115 = (t552 - t612 / 0.2e1 - t611 / 0.2e1) * m(7);
t473 = t544 * t464 + (t593 / 0.2e1 - t603 / 0.2e1) * t468 * mrSges(7,3) + (-t262 * t464 + (-t608 - t665) * t467) * t710 + (t214 * t464 + t372 * t597 + t468 * t512) * t708;
t477 = t222 * t711 + (t463 * t82 + t466 * t81) * t709 + t260 * t753 + t466 * t694;
t15 = -t436 + (t706 + t533) * t465 + t473 + t477;
t85 = (m(7) * t372 - t585 - t672) * t464;
t486 = qJD(1) * t15 + qJD(2) * t115 + qJD(3) * t85;
t138 = mrSges(6,3) + m(6) * qJ(5) + m(7) * (-t394 * t463 + t395 * t466) + t528;
t474 = -t451 + (t217 - 0.2e1 * t621) * t710 + ((-t395 * t468 + t74) * t466 + (t394 * t468 - t73) * t463) * t708 - t747;
t480 = t217 * t711 + (t463 * t84 + t466 * t83) * t709;
t18 = t474 + t480 - t739;
t51 = -t491 + t498 / 0.2e1;
t485 = t18 * qJD(1) - t51 * qJD(3) + t138 * qJD(4);
t484 = -mrSges(7,3) * t604 / 0.2e1 + t594 * t705;
t34 = t739 + t747;
t482 = t34 * qJD(1) - qJD(4) * t528;
t476 = (-t609 / 0.2e1 + t610 / 0.2e1) * mrSges(7,3) + t394 * t695 + t259 * t684 - t503;
t11 = (t73 / 0.2e1 + t84 / 0.2e1) * mrSges(7,2) + t722 * mrSges(7,1) + t476 + t503;
t30 = (t253 / 0.2e1 + t697) * mrSges(7,2) + (t698 + t696) * mrSges(7,1);
t479 = t11 * qJD(1) - t30 * qJD(3) - t87 * qJD(4);
t102 = (t611 + t612) * t708 + (m(6) + t708) * t598;
t65 = m(6) * t597 + m(7) * t734;
t45 = -t498 / 0.2e1 - t491 + (m(6) * pkin(8) + mrSges(6,2)) * t467 + (-t593 + t603) * mrSges(7,3);
t35 = t484 + t584 + t591;
t33 = t545 + t501;
t22 = t478 / 0.2e1 - t545 + t481 * t468 + t475 + (t401 + t731) * t676;
t17 = -t629 / 0.2e1 + t465 * t533 + t473 - t477;
t16 = t474 - t480 - t484 - t571;
t12 = t664 / 0.2e1 + t663 / 0.2e1 - t659 / 0.2e1 + t660 / 0.2e1 + t476 - t503;
t7 = t472 + t532;
t3 = t8 * qJD(3) - t10 * qJD(4) + qJD(5) * t568 + t14 * qJD(6);
t1 = t470 + t768;
t19 = [qJD(3) * t4 - qJD(4) * t5 - qJD(5) * t24 + qJD(6) * t9, t3, t1 * qJD(4) + t17 * qJD(5) + t7 * qJD(6) + ((t215 * t372 + t253 * t81 + t510 * t82) * t708 + t263 * t672 / 0.2e1) * t714 + (Ifges(4,5) + t541 * t467 + t542 * t464 + (-m(5) * pkin(3) + t623) * t441) * t580 + t517 + (-t81 * t631 + t231 * t687 + t228 * t688 + (Ifges(7,5) * t385 - Ifges(7,6) * t509) * t679 - Ifges(4,6) * t465 + t396 * t369 - t263 * t529 - t509 * t160 / 0.2e1 + t385 * t162 / 0.2e1 - pkin(3) * t370 + t372 * t176 + t510 * t260 + t253 * t261 + t215 * t224 + mrSges(4,2) * t606 - t82 * t632 + ((-t388 + t389) * t464 + m(5) * t511 + m(6) * t513) * pkin(8) + (t335 + t333) * t464 / 0.2e1 + (t524 + t522) * t678 + ((t387 + t390) * pkin(8) + t721) * t467 + t511 * mrSges(5,3) + t513 * mrSges(6,2)) * qJD(3), t1 * qJD(3) + (-t520 + t659 - t660 + t744 * t217 + t743 * t586 + (t609 - t610) * mrSges(7,3)) * qJD(4) + t16 * qJD(5) + t12 * qJD(6) + ((t394 * t83 + t395 * t84) * t708 + (-pkin(4) * t217 - qJ(5) * t586) * t710) * t713 + ((-mrSges(6,2) * qJ(5) + t764) * t467 + (mrSges(6,2) * pkin(4) - t742) * t464) * t579 + t516, t17 * qJD(3) + t16 * qJD(4) + t35 * qJD(6) + t495, t7 * qJD(3) + t12 * qJD(4) + t35 * qJD(5) + (-t520 - t663 - t664) * qJD(6) + t515; t3, t54 * qJD(3), t22 * qJD(4) + t102 * qJD(5) + t33 * qJD(6) + ((t583 + t608) * t710 + (t253 * t342 + t343 * t510 - t372 * t465) * t708 + (t583 - t669) * t712) * t714 + (-mrSges(4,2) + t756) * t580 + t514 + ((-t342 * t385 - t343 * t509) * mrSges(7,3) + (t585 + t623) * t465) * qJD(3), t22 * qJD(3) + t65 * qJD(5) + (t367 * t711 + (t340 * t395 + t394 * t341) * t708) * t713 + (t464 * t743 + t467 * t744) * t579 - t507 + t716, qJD(1) * t568 + t102 * qJD(3) + t65 * qJD(4), t33 * qJD(3) - t716 + t719; qJD(4) * t2 + qJD(5) * t15 + qJD(6) * t6 - t517, -qJD(4) * t21 + qJD(5) * t115 + qJD(6) * t32 - t514, qJD(4) * t13 + qJD(5) * t85 - qJD(6) * t766, t45 * qJD(5) - t772 + t489 + (m(7) * (-t253 * t395 + t394 * t510) + t395 * t631 - t394 * t632 - pkin(4) * t627 - mrSges(6,2) * t601 + (-m(6) * t518 - t529 - t530) * pkin(8) + t723 + t761) * qJD(4), qJD(4) * t45 + t486, -qJD(4) * t761 + t488 + t772; -qJD(3) * t2 + qJD(5) * t18 + qJD(6) * t11 - t516, t21 * qJD(3) + t507, -qJD(5) * t51 - qJD(6) * t30 - t489, qJD(5) * t138 - t620, t485, t479 + t620; -t15 * qJD(3) - t18 * qJD(4) - t34 * qJD(6) - t495, -qJD(1) * t670 / 0.2e1 - t115 * qJD(3), qJD(4) * t51 - t486, -t485 + t746, 0, -t482 - t746; -qJD(3) * t6 - qJD(4) * t11 + qJD(5) * t34 - t515, -t32 * qJD(3) - t719, qJD(4) * t30 - t488, -qJD(5) * t528 - t479, t482, 0;];
Cq  = t19;
