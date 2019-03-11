% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:54:18
% EndTime: 2019-03-09 06:54:39
% DurationCPUTime: 14.28s
% Computational Cost: add. (35468->570), mult. (66621->738), div. (0->0), fcn. (75313->10), ass. (0->364)
t392 = sin(qJ(6));
t394 = cos(qJ(6));
t653 = sin(qJ(4));
t654 = sin(qJ(3));
t656 = cos(qJ(4));
t657 = cos(qJ(3));
t365 = -t653 * t657 - t654 * t656;
t393 = sin(qJ(5));
t436 = t653 * t654 - t656 * t657;
t655 = cos(qJ(5));
t322 = t365 * t393 - t655 * t436;
t530 = -cos(pkin(11)) * pkin(1) - pkin(2);
t367 = -pkin(3) * t657 + t530;
t328 = pkin(4) * t436 + t367;
t421 = -t365 * t655 - t393 * t436;
t181 = -pkin(5) * t322 - pkin(10) * t421 + t328;
t491 = sin(pkin(11)) * pkin(1) + pkin(7);
t455 = t654 * t491;
t428 = -pkin(8) * t654 - t455;
t456 = t657 * t491;
t429 = pkin(8) * t657 + t456;
t305 = t653 * t428 + t429 * t656;
t403 = t436 * pkin(9) - t305;
t682 = t656 * t428 - t653 * t429;
t728 = t365 * pkin(9) + t682;
t770 = t393 * t728 - t403 * t655;
t68 = t392 * t181 + t394 * t770;
t613 = t394 * t68;
t67 = t181 * t394 - t392 * t770;
t483 = -t392 * t67 + t613;
t790 = t770 - t483;
t615 = t394 * mrSges(7,2);
t618 = t392 * mrSges(7,1);
t486 = t615 + t618;
t215 = t486 * t421;
t614 = t394 * mrSges(7,3);
t717 = mrSges(7,1) * t421 - t322 * t614;
t731 = t486 * t322;
t771 = t393 * t403 + t655 * t728;
t789 = t770 * t215 + t67 * t717 - t771 * t731;
t318 = Ifges(6,5) * t322;
t643 = Ifges(7,2) * t394;
t645 = Ifges(7,4) * t392;
t370 = t643 + t645;
t572 = t392 * t370;
t360 = -t572 / 0.2e1;
t386 = Ifges(7,4) * t394;
t647 = Ifges(7,1) * t392;
t372 = t386 + t647;
t561 = t394 * t372;
t729 = t360 + t561 / 0.2e1;
t369 = Ifges(7,5) * t392 + Ifges(7,6) * t394;
t699 = t421 * t369;
t705 = Ifges(6,6) * t421;
t693 = -Ifges(7,2) * t392 + t386;
t716 = Ifges(7,6) * t421 + t322 * t693;
t746 = t394 * t716;
t373 = Ifges(7,1) * t394 - t645;
t718 = Ifges(7,5) * t421 + t322 * t373;
t748 = t392 * t718;
t755 = t699 / 0.2e1 - t705 + t746 / 0.2e1 + t748 / 0.2e1;
t439 = -Ifges(5,5) * t436 + Ifges(5,6) * t365 + t729 * t322 + t318 + t755;
t733 = t682 * mrSges(5,2);
t754 = t305 * mrSges(5,1);
t368 = -mrSges(7,1) * t394 + t392 * mrSges(7,2);
t776 = t770 * t368;
t779 = t771 * mrSges(6,2);
t780 = t770 * mrSges(6,1);
t773 = t776 - t780 - t779;
t787 = t439 - t733 - t754 + t773;
t786 = m(6) * t328 - mrSges(6,1) * t322 + mrSges(6,2) * t421;
t676 = -m(7) / 0.2e1;
t596 = t771 * t421;
t784 = -t790 * t322 - t596;
t674 = m(6) * pkin(4);
t766 = t771 * t393;
t783 = (-t655 * t770 + t766) * t674;
t573 = t392 * t322;
t708 = mrSges(7,2) * t421;
t217 = -mrSges(7,3) * t573 - t708;
t703 = t421 * mrSges(6,1);
t734 = t322 * mrSges(6,2);
t223 = t703 + t734;
t324 = -t365 * mrSges(5,1) - mrSges(5,2) * t436;
t665 = -t322 / 0.2e1;
t156 = -Ifges(7,5) * t322 + t373 * t421;
t566 = t394 * t156;
t153 = -Ifges(7,6) * t322 + t421 * t693;
t579 = t392 * t153;
t712 = t421 / 0.2e1;
t686 = Ifges(6,1) * t712 + Ifges(6,4) * t322 - t579 / 0.2e1 + t566 / 0.2e1;
t385 = Ifges(7,5) * t394;
t641 = Ifges(7,6) * t392;
t694 = t385 - t641;
t658 = t394 / 0.2e1;
t660 = -t392 / 0.2e1;
t735 = t718 * t658 + t716 * t660 + t694 * t712 - Ifges(6,4) * t421 + (Ifges(6,2) + Ifges(7,3)) * t665;
t670 = Ifges(7,3) / 0.2e1;
t738 = Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1 - t670;
t782 = t328 * t223 + t367 * t324 + t68 * t217 + t436 ^ 2 * Ifges(5,4) + (t665 * t694 + t686) * t322 + (t322 * t738 + t735) * t421 + t789;
t681 = t779 / 0.2e1 - t776 / 0.2e1 + t780 / 0.2e1;
t781 = pkin(5) * t770;
t626 = t322 * mrSges(6,3);
t551 = t656 * pkin(3);
t384 = t551 + pkin(4);
t526 = t653 * t393;
t349 = -pkin(3) * t526 + t384 * t655;
t343 = -pkin(5) - t349;
t751 = t343 * t770;
t550 = t655 * pkin(4);
t383 = -t550 - pkin(5);
t750 = t383 * t770;
t764 = t392 * t771;
t775 = t770 * t771;
t493 = t653 * t655;
t350 = pkin(3) * t493 + t393 * t384;
t767 = t771 * t350;
t765 = t771 * t394;
t265 = -t731 / 0.2e1;
t564 = t394 * t217;
t749 = t392 * t717;
t769 = -t749 / 0.2e1;
t724 = t564 / 0.2e1 + t769;
t726 = t561 / 0.4e1 - t572 / 0.4e1;
t756 = t699 / 0.4e1 - t705 / 0.2e1 + t746 / 0.4e1 + t748 / 0.4e1;
t435 = pkin(10) * t724 + t322 * t726 + t318 / 0.2e1 + pkin(5) * t265 + t756;
t710 = pkin(5) * t421;
t227 = -pkin(10) * t322 + t710;
t652 = pkin(4) * t365;
t187 = t227 - t652;
t82 = t392 * t187 + t765;
t608 = t82 * t394;
t81 = t187 * t394 - t764;
t609 = t81 * t392;
t481 = t608 - t609;
t774 = (pkin(10) * t481 - t781) * t676 - t435;
t772 = (Ifges(6,5) / 0.2e1 + t726) * t322;
t651 = pkin(4) * t393;
t382 = pkin(10) + t651;
t763 = t382 * t749;
t389 = t392 ^ 2;
t390 = t394 ^ 2;
t557 = t389 + t390;
t700 = t421 * t368;
t762 = t700 - t223;
t467 = t615 / 0.2e1 + t618 / 0.2e1;
t450 = t467 * t322;
t741 = t265 - t450;
t759 = t741 * qJD(6);
t758 = -t733 / 0.2e1 - t754 / 0.2e1;
t659 = -t394 / 0.2e1;
t667 = t731 / 0.2e1;
t707 = mrSges(7,3) * t557;
t742 = -t763 / 0.2e1;
t619 = t390 * mrSges(7,3);
t620 = t389 * mrSges(7,3);
t739 = -t620 / 0.2e1 - t619 / 0.2e1;
t737 = t700 / 0.2e1 - t703 / 0.2e1 - t734 / 0.2e1;
t730 = t557 * t322;
t736 = pkin(10) * t730 - t710;
t675 = m(7) / 0.2e1;
t732 = t322 * t393;
t593 = t322 * t421;
t698 = (-t385 / 0.2e1 + t641 / 0.2e1) * t322;
t344 = pkin(10) + t350;
t723 = t557 * t344;
t721 = t373 * t660 + t693 * t659;
t490 = (t390 / 0.2e1 + t389 / 0.2e1) * mrSges(7,3);
t701 = t383 * t421;
t502 = t557 * t349;
t457 = t383 * t486;
t461 = t572 / 0.2e1 - t561 / 0.2e1;
t684 = t461 + t721;
t414 = -t457 / 0.2e1 + t684;
t458 = t343 * t486;
t447 = -t458 / 0.2e1;
t355 = (t655 * t656 - t526) * pkin(3);
t448 = t467 * t355;
t100 = t447 - t448 + t414;
t587 = t421 * t394;
t588 = t421 * t392;
t212 = -mrSges(7,1) * t587 + mrSges(7,2) * t588;
t422 = -t771 * t486 / 0.2e1 - t579 / 0.4e1 + t566 / 0.4e1;
t464 = -t322 * t694 / 0.4e1;
t415 = t464 + t422;
t416 = (t373 / 0.4e1 - t370 / 0.4e1 - t643 / 0.4e1) * t394 + (-t372 / 0.4e1 - t693 / 0.4e1 - t386 / 0.2e1 - t647 / 0.4e1) * t392;
t549 = mrSges(7,3) * t588;
t218 = mrSges(7,2) * t322 - t549;
t221 = -mrSges(7,1) * t322 - mrSges(7,3) * t587;
t462 = t218 * t660 + t221 * t659;
t400 = t462 * t382 + (-t382 * t490 + t416) * t421 - t383 * t212 / 0.2e1 + t415;
t671 = mrSges(7,2) / 0.2e1;
t672 = -mrSges(7,1) / 0.2e1;
t470 = t671 * t82 + t672 * t81;
t640 = Ifges(7,3) * t421;
t691 = t698 - t640 / 0.2e1;
t17 = t400 + t470 + t691;
t495 = -t721 + t729;
t222 = t457 + t495;
t110 = t667 - t450;
t556 = t110 * qJD(2);
t692 = t17 * qJD(1) - t100 * qJD(3) + t222 * qJD(4) - t556;
t664 = -t343 / 0.2e1;
t401 = t462 * t344 + (-t344 * t490 + t416) * t421 + t212 * t664 + t415;
t388 = t654 * pkin(3);
t185 = t187 + t388;
t77 = t392 * t185 + t765;
t668 = t77 / 0.2e1;
t76 = t185 * t394 - t764;
t669 = -t76 / 0.2e1;
t471 = mrSges(7,1) * t669 + mrSges(7,2) * t668;
t15 = t401 + t471 + t691;
t197 = t458 + t495;
t690 = t15 * qJD(1) + t197 * qJD(3) - t556;
t563 = t394 * t218;
t513 = t563 / 0.2e1;
t574 = t392 * t221;
t460 = t513 - t574 / 0.2e1;
t617 = t392 * mrSges(7,3);
t216 = -t322 * t617 - t708;
t565 = t394 * t216;
t688 = t565 / 0.2e1 + t769;
t610 = t77 * t394;
t611 = t76 * t392;
t482 = t610 - t611;
t87 = t227 * t394 - t764;
t88 = t392 * t227 + t765;
t480 = -t392 * t87 + t394 * t88;
t521 = -t573 / 0.2e1;
t540 = t385 / 0.2e1;
t685 = Ifges(7,6) * t521 + t322 * t540 + t640 / 0.2e1;
t683 = (-t653 * mrSges(5,1) - t656 * mrSges(5,2)) * pkin(3);
t680 = (t610 / 0.2e1 - t611 / 0.2e1) * mrSges(7,3) - t681;
t299 = t322 * t620;
t300 = t322 * t619;
t466 = t736 * t675 + t299 / 0.2e1 + t300 / 0.2e1 + t737;
t677 = m(6) / 0.2e1;
t673 = m(7) * pkin(4);
t666 = t215 / 0.2e1;
t663 = t349 / 0.2e1;
t354 = (t393 * t656 + t493) * pkin(3);
t662 = -t354 / 0.2e1;
t661 = t355 / 0.2e1;
t629 = t421 * mrSges(6,3);
t624 = t349 * mrSges(6,2);
t623 = t350 * mrSges(6,1);
t622 = t354 * mrSges(6,1);
t621 = t355 * mrSges(6,2);
t51 = m(7) * (t421 * t730 - t593);
t559 = t51 * qJD(2);
t607 = -t110 * qJD(6) - t559;
t606 = t559 + t759;
t30 = t212 * t665 + (t421 * t490 - t462) * t421;
t605 = qJD(1) * t30;
t463 = t421 * t666 + t665 * t731;
t12 = t322 * t513 + t221 * t521 + t784 * t675 - t596 * t677 + t463 + (t482 * t675 + t677 * t771 + t724) * t421;
t604 = t12 * qJD(1);
t13 = t724 * t421 + t460 * t322 + (t421 * t481 + t784) * t675 + t463;
t603 = t13 * qJD(1);
t591 = t322 * t354;
t586 = t343 * t731;
t583 = t350 * t368;
t582 = t354 * t771;
t581 = t354 * t368;
t580 = t383 * t731;
t560 = t394 * t382;
t558 = t350 * t322 - t349 * t421;
t555 = qJD(3) + qJD(4);
t554 = mrSges(6,1) * t651;
t553 = mrSges(7,3) * t609;
t552 = mrSges(7,3) * t608;
t548 = t349 * t626;
t547 = t350 * t629;
t546 = t651 / 0.2e1;
t535 = -t617 / 0.2e1;
t534 = t614 / 0.2e1;
t531 = t322 * t723 + t343 * t421;
t529 = t392 * t655;
t528 = t394 * t655;
t509 = t560 / 0.2e1;
t501 = t557 * t355;
t500 = t557 * t382;
t499 = t629 * t651;
t497 = mrSges(6,2) * t550;
t496 = t550 / 0.2e1;
t492 = -t529 / 0.2e1;
t489 = t550 * t626;
t488 = t82 * t534 + t81 * t535 - t681;
t465 = pkin(5) * t486;
t479 = -t465 / 0.2e1 + t495;
t409 = -Ifges(5,4) * t365 + (Ifges(5,1) - Ifges(5,2)) * t436;
t431 = pkin(3) * t436;
t454 = mrSges(4,1) * t654 + mrSges(4,2) * t657;
t3 = (-mrSges(5,2) * t388 + t409) * t365 + m(7) * (t67 * t76 + t68 * t77 - t775) + m(5) * t367 * t388 + t77 * t218 + t76 * t221 + t530 * t454 + t657 ^ 2 * Ifges(4,4) + (mrSges(5,1) * t431 + (-Ifges(4,2) + Ifges(4,1)) * t657 - Ifges(4,4) * t654) * t654 + t782 + t786 * (t388 - t652);
t478 = t3 * qJD(1) + t12 * qJD(2);
t4 = (-t786 * pkin(4) + t409) * t365 + m(7) * (t67 * t81 + t68 * t82 - t775) + t82 * t218 + t81 * t221 + t782;
t477 = t4 * qJD(1) + t13 * qJD(2);
t10 = m(7) * (t67 * t87 + t68 * t88 - t775) + t88 * t218 + t68 * t216 + t87 * t221 + (t328 * mrSges(6,1) + t735) * t421 - (-t328 * mrSges(6,2) - t421 * t738 - t686 - t698) * t322 + t789;
t18 = ((-t771 + t480) * t675 + t666 + t688) * t421 - (t790 * t675 - t460 + t667) * t322;
t476 = t10 * qJD(1) + t18 * qJD(2);
t11 = -t771 * t212 + t68 * t221 + (t392 * t156 / 0.2e1 + t153 * t658 + mrSges(7,3) * t613 + t369 * t665 + (-t372 * t659 + t360) * t421) * t421 + (-t218 - t549) * t67;
t475 = -t11 * qJD(1) - t30 * qJD(2);
t52 = m(7) * (0.1e1 - t557) * t593;
t474 = t18 * qJD(1) - t52 * qJD(2);
t469 = t671 * t88 + t672 * t87;
t453 = t465 / 0.2e1;
t452 = t557 * t655;
t449 = t467 * t349;
t405 = ((t343 + t502) * t675 + t368 / 0.2e1 - mrSges(6,1) / 0.2e1) * t421 - ((t350 - t723) * t675 + mrSges(6,2) / 0.2e1 - t490) * t322;
t28 = t405 - t466;
t425 = mrSges(7,3) * t502 + t583 - t623 - t624;
t55 = m(7) * (t343 * t350 + t344 * t502) + t425;
t407 = (t344 * t480 + t349 * t483 + t751 - t767) * t676 + t731 * t664 - t350 * t215 / 0.2e1;
t412 = (pkin(10) * t482 - t781) * t675 + t435;
t6 = -t772 + (Ifges(6,6) / 0.2e1 - t369 / 0.4e1) * t421 + t407 + t412 + (t221 * t663 + t344 * t717 / 0.2e1 - t718 / 0.4e1) * t392 + (-t349 * t218 / 0.2e1 - t344 * t216 / 0.2e1 - t716 / 0.4e1) * t394 + ((-t88 / 0.2e1 + t668) * t394 + (t87 / 0.2e1 + t669) * t392) * mrSges(7,3);
t446 = -t6 * qJD(1) + t28 * qJD(2) + t55 * qJD(3);
t406 = (t355 * t421 + t558 - t591) * t677 + (t421 * t501 + t531 - t591) * t675;
t440 = m(7) * (t322 * t500 + t701);
t442 = (-t421 * t655 + t732) * t674;
t417 = -t440 / 0.2e1 - t442 / 0.2e1;
t34 = t406 + t417;
t56 = (-mrSges(6,1) + t368) * t354 + t683 + (-mrSges(6,2) + t707) * t355 + m(7) * (t343 * t354 + t344 * t501) + m(6) * (-t349 * t354 + t350 * t355);
t395 = -m(6) * (t767 - t582 + (-t349 + t355) * t770) / 0.2e1 + (t355 * t483 - t582 + t751) * t676 - t586 / 0.2e1 + t215 * t662 + t548 / 0.2e1 + t547 / 0.2e1 + t574 * t661 - t355 * t563 / 0.2e1 + (t481 * t676 - t724) * t344 + (-t322 * t661 + t421 * t662) * mrSges(6,3) - t758;
t397 = (t382 * t482 + t750) * t675 + t580 / 0.2e1 + t783 / 0.2e1 + t742 + t217 * t509 - t499 / 0.2e1 - t489 / 0.2e1 + t680 + t758;
t420 = t553 / 0.2e1 - t552 / 0.2e1 + t681;
t7 = t397 + t395 + t420;
t445 = -t7 * qJD(1) + t34 * qJD(2) + t56 * qJD(3);
t413 = t368 * t651 + t550 * t707 - t497 - t554;
t198 = (t382 * t452 + t383 * t393) * t673 + t413;
t402 = (t701 + t382 * t730 + (t421 * t452 - t732) * pkin(4)) * t675 - t739 * t322 + t737;
t33 = t402 - t466;
t398 = (t383 * t350 + t349 * t500 + (t343 * t393 + t344 * t452) * pkin(4)) * t675 - t624 / 0.2e1 - t623 / 0.2e1 + t583 / 0.2e1 - t554 / 0.2e1 + t368 * t546 - t497 / 0.2e1 + (t496 + t663) * t707;
t408 = (-pkin(5) * t354 + pkin(10) * t501) * t676 + t622 / 0.2e1 - t581 / 0.2e1 + t621 / 0.2e1 + t739 * t355;
t43 = t398 + t408;
t426 = t88 * t534 + t87 * t535 - t681 + t756 + t772;
t399 = t215 * t546 + t216 * t509 + t742 + t426 + (t750 + t480 * t382 + (t528 * t68 - t529 * t67 - t766) * pkin(4)) * t675 + pkin(4) * t221 * t492 + t496 * t563 + t383 * t667;
t9 = t399 + t420 + t774;
t433 = t9 * qJD(1) + t33 * qJD(2) + t43 * qJD(3) + t198 * qJD(4);
t102 = t447 + t453 - t449 + t684;
t427 = (-mrSges(7,2) * t528 / 0.2e1 + mrSges(7,1) * t492) * pkin(4);
t167 = t453 + t427 + t414;
t410 = -pkin(10) * t490 + t416;
t411 = t462 * pkin(10) + pkin(5) * t212 / 0.2e1 + t422;
t20 = -(-0.3e1 / 0.4e1 * t641 + t385 / 0.4e1 + t540) * t322 + (-Ifges(7,3) / 0.2e1 + t410) * t421 + t411 + t469;
t228 = -t465 + t495;
t432 = t20 * qJD(1) - t102 * qJD(3) - t167 * qJD(4) + t228 * qJD(5);
t418 = t299 + t300 - t324 + t762;
t347 = t457 / 0.2e1;
t325 = t458 / 0.2e1;
t168 = t347 + t427 + t479;
t103 = t325 - t449 + t479;
t101 = t325 + t347 - t448 + t495;
t42 = t398 - t408;
t32 = t402 + t466;
t29 = t406 - t417 + t418;
t27 = t405 + t466;
t19 = t411 + t464 - t469 + (t410 + t670) * t421 - t698;
t16 = t400 - t470 + t685;
t14 = t401 - t471 + t685;
t8 = t399 + t488 - t774;
t5 = t344 * t688 + t349 * t460 - t407 + t412 + t426 + t680;
t2 = t439 + t397 - t395 + t488;
t1 = qJD(3) * t12 + qJD(4) * t13 + qJD(5) * t18 - qJD(6) * t30;
t21 = [qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t10 - qJD(6) * t11, t1 (mrSges(4,2) * t455 - mrSges(4,1) * t456 + Ifges(4,5) * t657 - Ifges(4,6) * t654 + m(5) * (-t305 * t656 + t653 * t682) * pkin(3) - t547 - t548 + m(6) * (-t349 * t770 + t767) + m(7) * t751 + t586 + (m(7) * t482 + t564 - t749) * t344 + t482 * mrSges(7,3) + (pkin(3) * t365 * t653 + t431 * t656) * mrSges(5,3) + t787) * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t14 * qJD(6) + t478, t2 * qJD(3) + (t783 + t552 + t217 * t560 - t553 - t763 + m(7) * (t382 * t481 + t750) - t499 - t489 + t580 + t787) * qJD(4) + t8 * qJD(5) + t16 * qJD(6) + t477, t5 * qJD(3) + t8 * qJD(4) + t19 * qJD(6) + t476 + (-(-Ifges(6,5) + t461) * t322 + (-m(7) * t770 - t731) * pkin(5) + (m(7) * t480 + t565 - t749) * pkin(10) + t480 * mrSges(7,3) + t755 + t773) * qJD(5), t14 * qJD(3) + t16 * qJD(4) + t19 * qJD(5) + (-t68 * mrSges(7,1) - t67 * mrSges(7,2) - t699) * qJD(6) + t475; t1, -qJD(5) * t52 + t51 * t555, t29 * qJD(4) + t27 * qJD(5) + t604 + t606 + (t418 - t454 + 0.2e1 * t531 * t675 + 0.2e1 * t558 * t677 + (t365 * t551 - t431 * t653) * m(5)) * qJD(3), t603 + t29 * qJD(3) + (t440 + t442 + t418) * qJD(4) + t32 * qJD(5) + t606, t27 * qJD(3) + t32 * qJD(4) + (m(7) * t736 + mrSges(7,3) * t730 + t762) * qJD(5) + t759 + t474, qJD(6) * t212 - t605 + (qJD(5) + t555) * t741; -qJD(4) * t7 - qJD(5) * t6 + qJD(6) * t15 - t478, qJD(4) * t34 + qJD(5) * t28 - t604 + t607, qJD(4) * t56 + qJD(5) * t55 + qJD(6) * t197 (t581 - t621 - t622 + (m(7) * t383 - t655 * t674) * t354 + (m(7) * t500 + t393 * t674 + t619 + t620) * t355 + t683) * qJD(4) + t42 * qJD(5) + t101 * qJD(6) + t445, t42 * qJD(4) + (m(7) * (-pkin(5) * t350 + pkin(10) * t502) + t425) * qJD(5) + t103 * qJD(6) + t446, t101 * qJD(4) + t103 * qJD(5) + (t344 * t368 + t694) * qJD(6) + t690; qJD(3) * t7 + qJD(5) * t9 + qJD(6) * t17 - t477, -qJD(3) * t34 + qJD(5) * t33 - t603 + t607, qJD(5) * t43 - qJD(6) * t100 - t445, qJD(5) * t198 + qJD(6) * t222 ((-pkin(5) * t393 + pkin(10) * t452) * t673 + t413) * qJD(5) + t168 * qJD(6) + t433, t168 * qJD(5) + (t368 * t382 + t694) * qJD(6) + t692; qJD(3) * t6 - qJD(4) * t9 + qJD(6) * t20 - t476, -qJD(3) * t28 - qJD(4) * t33 - t474, -qJD(4) * t43 - qJD(6) * t102 - t446, -qJD(6) * t167 - t433, t228 * qJD(6) (pkin(10) * t368 + t694) * qJD(6) + t432; -qJD(3) * t15 - qJD(4) * t17 - qJD(5) * t20 - t475, t110 * t555 + t605, qJD(4) * t100 + qJD(5) * t102 - t690, qJD(5) * t167 - t692, -t432, 0;];
Cq  = t21;
