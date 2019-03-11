% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:53
% EndTime: 2019-03-08 23:16:26
% DurationCPUTime: 20.29s
% Computational Cost: add. (30230->866), mult. (70383->1224), div. (0->0), fcn. (78633->12), ass. (0->438)
t519 = sin(qJ(6));
t523 = cos(qJ(6));
t516 = sin(pkin(12));
t518 = cos(pkin(12));
t520 = sin(qJ(4));
t524 = cos(qJ(4));
t561 = t516 * t520 - t518 * t524;
t562 = t516 * t524 + t518 * t520;
t352 = -t519 * t561 + t523 * t562;
t591 = -t519 * t562 - t523 * t561;
t621 = Ifges(7,5) * t591 - Ifges(7,6) * t352;
t716 = -qJ(5) - pkin(9);
t480 = t716 * t520;
t482 = t716 * t524;
t386 = t516 * t480 - t518 * t482;
t313 = pkin(10) * t561 - t386;
t813 = t518 * t480 + t516 * t482;
t835 = -pkin(10) * t562 + t813;
t186 = t313 * t523 - t519 * t835;
t847 = t313 * t519 + t523 * t835;
t877 = t186 * mrSges(7,1) - t847 * mrSges(7,2);
t31 = t621 + t877;
t879 = t31 * qJD(6);
t521 = sin(qJ(3));
t525 = cos(qJ(3));
t478 = -pkin(3) * t525 - pkin(9) * t521 - pkin(2);
t460 = t524 * t478;
t629 = t521 * t524;
t586 = -qJ(5) * t629 + t460;
t630 = t520 * t525;
t615 = pkin(8) * t630;
t378 = t586 - t615;
t627 = t524 * t525;
t415 = pkin(8) * t627 + t478 * t520;
t631 = t520 * t521;
t379 = -qJ(5) * t631 + t415;
t634 = t518 * t379;
t243 = -t378 * t516 - t634;
t434 = t562 * t521;
t722 = pkin(10) * t434;
t201 = t243 + t722;
t364 = t516 * t379;
t244 = t518 * t378 - t364;
t433 = t516 * t631 - t518 * t629;
t723 = pkin(10) * t433;
t202 = t244 + t723;
t102 = t201 * t519 + t202 * t523;
t361 = (-pkin(8) * t520 - pkin(4)) * t525 + t586;
t230 = t518 * t361 - t364;
t185 = -pkin(5) * t525 + t230 + t723;
t231 = t516 * t361 + t634;
t193 = t231 - t722;
t81 = t185 * t523 - t193 * t519;
t878 = t102 - t81;
t101 = t201 * t523 - t202 * t519;
t82 = t185 * t519 + t193 * t523;
t875 = t101 + t82;
t517 = sin(pkin(6));
t522 = sin(qJ(2));
t636 = t517 * t522;
t657 = cos(pkin(6));
t450 = t521 * t657 + t525 * t636;
t526 = cos(qJ(2));
t635 = t517 * t526;
t376 = -t450 * t520 - t524 * t635;
t377 = t450 * t524 - t520 * t635;
t240 = t376 * t516 + t377 * t518;
t593 = t518 * t376 - t377 * t516;
t123 = t240 * t523 + t519 * t593;
t831 = -t240 * t519 + t523 * t593;
t28 = -t123 * mrSges(7,1) - t831 * mrSges(7,2);
t874 = t28 * qJD(6);
t592 = t433 * t519 - t523 * t434;
t303 = t433 * t523 + t434 * t519;
t679 = Ifges(7,4) * t303;
t167 = Ifges(7,2) * t592 - t525 * Ifges(7,6) - t679;
t510 = t521 * pkin(8);
t472 = pkin(4) * t631 + t510;
t362 = pkin(5) * t434 + t472;
t726 = pkin(4) * t524;
t504 = -pkin(3) - t726;
t418 = pkin(5) * t561 + t504;
t266 = -mrSges(7,1) * t525 + mrSges(7,3) * t303;
t779 = -t266 / 0.2e1;
t851 = Ifges(7,1) * t592 + t679;
t820 = t591 * mrSges(7,2);
t849 = t352 * mrSges(7,1) + t820;
t870 = t849 / 0.2e1;
t840 = t303 * mrSges(7,1);
t848 = t592 * mrSges(7,2) - t840;
t871 = t848 / 0.2e1;
t873 = -t186 * t779 + (t851 / 0.4e1 - t167 / 0.4e1) * t352 + t362 * t870 + t418 * t871;
t872 = t123 / 0.2e1;
t673 = t352 * mrSges(7,3);
t867 = t362 * t848;
t866 = t418 * t849;
t864 = -t230 + t244;
t341 = Ifges(7,4) * t591;
t222 = t352 * Ifges(7,1) + t341;
t832 = -Ifges(7,2) * t352 + t341;
t863 = t832 + t222;
t672 = t352 * Ifges(7,4);
t220 = Ifges(7,2) * t591 + t672;
t852 = Ifges(7,1) * t591 - t672;
t862 = -t220 / 0.4e1 + t852 / 0.4e1;
t775 = t840 / 0.2e1;
t512 = t520 ^ 2;
t514 = t524 ^ 2;
t811 = t514 + t512;
t859 = pkin(9) * t811;
t858 = t434 * mrSges(6,3);
t508 = Ifges(5,4) * t524;
t580 = -t520 * Ifges(5,1) - t508;
t732 = t524 / 0.2e1;
t856 = t580 * t732;
t264 = mrSges(7,2) * t525 + mrSges(7,3) * t592;
t780 = t264 / 0.2e1;
t855 = t831 * t780;
t837 = t243 + t231;
t626 = t525 * t526;
t412 = (-t520 * t626 + t522 * t524) * t517;
t413 = (t520 * t522 + t524 * t626) * t517;
t268 = t412 * t518 - t413 * t516;
t269 = t412 * t516 + t413 * t518;
t157 = t268 * t523 - t269 * t519;
t158 = t268 * t519 + t269 * t523;
t503 = pkin(4) * t518 + pkin(5);
t728 = pkin(4) * t516;
t440 = t503 * t523 - t519 * t728;
t441 = t503 * t519 + t523 * t728;
t795 = m(6) * pkin(4);
t616 = t795 / 0.2e1;
t845 = -mrSges(7,2) / 0.2e1;
t846 = mrSges(7,1) / 0.2e1;
t625 = t157 * t846 + t158 * t845;
t776 = -t269 / 0.2e1;
t792 = -mrSges(5,2) / 0.2e1;
t793 = mrSges(6,1) / 0.2e1;
t794 = mrSges(5,1) / 0.2e1;
t796 = m(7) / 0.2e1;
t854 = -(t157 * t440 + t158 * t441) * t796 - t268 * t793 - mrSges(6,2) * t776 - t412 * t794 - t413 * t792 - (t268 * t518 + t269 * t516) * t616 - t625;
t449 = t521 * t636 - t525 * t657;
t300 = t562 * t449;
t301 = t561 * t449;
t165 = t300 * t523 - t301 * t519;
t166 = t300 * t519 + t301 * t523;
t617 = -t795 / 0.2e1;
t624 = t165 * t846 + t166 * t845;
t797 = -m(7) / 0.2e1;
t853 = (t165 * t440 + t166 * t441) * t797 - t300 * mrSges(6,1) / 0.2e1 + t301 * mrSges(6,2) / 0.2e1 + (t300 * t518 + t301 * t516) * t617 - t624;
t821 = Ifges(7,5) * t592;
t842 = Ifges(7,6) * t303;
t623 = t821 + t842;
t294 = Ifges(7,4) * t592;
t169 = -Ifges(7,1) * t303 - t525 * Ifges(7,5) + t294;
t833 = Ifges(7,2) * t303 + t294;
t850 = t833 + t169;
t599 = t842 / 0.2e1 + t821 / 0.2e1;
t844 = -t352 / 0.2e1;
t761 = t352 / 0.2e1;
t393 = -mrSges(6,1) * t525 + t433 * mrSges(6,3);
t757 = t393 / 0.2e1;
t753 = t433 / 0.2e1;
t843 = t521 / 0.2e1;
t608 = t820 / 0.2e1;
t669 = t561 * mrSges(6,3);
t668 = t562 * mrSges(6,3);
t827 = -t592 / 0.2e1;
t838 = t831 * t827;
t435 = t562 * t525;
t436 = t561 * t525;
t305 = -t435 * t523 + t436 * t519;
t308 = -t435 * t519 - t436 * t523;
t557 = -Ifges(7,5) * t308 / 0.2e1 - Ifges(7,6) * t305 / 0.2e1;
t489 = pkin(3) * t521 - pkin(9) * t525;
t423 = pkin(8) * t631 + t524 * t489;
t368 = pkin(4) * t521 - qJ(5) * t627 + t423;
t424 = -pkin(8) * t629 + t520 * t489;
t382 = -qJ(5) * t630 + t424;
t234 = t518 * t368 - t382 * t516;
t190 = pkin(5) * t521 + pkin(10) * t436 + t234;
t235 = t516 * t368 + t518 * t382;
t198 = -pkin(10) * t435 + t235;
t89 = t190 * t523 - t198 * t519;
t90 = t190 * t519 + t198 * t523;
t583 = Ifges(7,3) * t843 + t90 * t845 + t846 * t89 - t557;
t610 = t521 * t635;
t643 = t413 * t524;
t644 = t412 * t520;
t799 = -m(6) / 0.2e1;
t834 = -(t643 / 0.2e1 - t644 / 0.2e1) * mrSges(5,3) - m(5) * (-pkin(3) * t610 + (t643 - t644) * pkin(9)) / 0.2e1 + (t268 * t813 + t386 * t269 + t504 * t610) * t799 + (t157 * t847 - t158 * t186 + t418 * t610) * t797;
t830 = t230 / 0.2e1;
t754 = -t433 / 0.4e1;
t742 = t562 / 0.4e1;
t828 = t591 / 0.2e1;
t825 = t592 / 0.2e1;
t702 = Ifges(5,6) * t520;
t706 = Ifges(5,5) * t524;
t559 = t706 / 0.2e1 - t702 / 0.2e1;
t818 = Ifges(4,4) - t559;
t218 = -mrSges(7,1) * t591 + mrSges(7,2) * t352;
t372 = mrSges(6,1) * t561 + mrSges(6,2) * t562;
t815 = t218 + t372;
t453 = Ifges(6,4) * t561;
t374 = Ifges(6,1) * t562 - t453;
t814 = -Ifges(6,2) * t562 + t374 - t453;
t812 = -Ifges(5,2) * t520 + t508;
t810 = -Ifges(6,5) * t561 - Ifges(6,6) * t562 + t621;
t422 = Ifges(6,4) * t434;
t298 = -t433 * Ifges(6,1) - t525 * Ifges(6,5) - t422;
t809 = Ifges(6,2) * t433 + t298 - t422;
t808 = -Ifges(6,5) * t434 + Ifges(6,6) * t433 + t623;
t807 = -t423 * t520 + t424 * t524;
t663 = t525 * mrSges(4,2);
t484 = t521 * mrSges(4,1) + t663;
t806 = -mrSges(5,1) * t524 + mrSges(5,2) * t520;
t804 = t386 * t753 + t813 * t434 / 0.2e1;
t802 = 0.2e1 * m(7);
t801 = 2 * qJD(3);
t800 = m(5) / 0.2e1;
t798 = m(6) / 0.2e1;
t791 = -mrSges(7,3) / 0.2e1;
t789 = -t831 / 0.2e1;
t788 = -t123 / 0.2e1;
t786 = -t847 / 0.2e1;
t785 = t186 / 0.2e1;
t784 = -t218 / 0.2e1;
t782 = -t593 / 0.2e1;
t781 = -t264 / 0.2e1;
t778 = t266 / 0.2e1;
t267 = mrSges(7,1) * t521 - mrSges(7,3) * t308;
t777 = t267 / 0.2e1;
t774 = t303 / 0.2e1;
t771 = t305 / 0.2e1;
t768 = -t303 / 0.2e1;
t767 = t308 / 0.2e1;
t765 = -t591 / 0.2e1;
t760 = -t372 / 0.2e1;
t391 = mrSges(6,2) * t525 - t858;
t759 = t391 / 0.2e1;
t392 = -mrSges(6,2) * t521 - mrSges(6,3) * t435;
t758 = t392 / 0.2e1;
t394 = mrSges(6,1) * t521 + mrSges(6,3) * t436;
t756 = t394 / 0.2e1;
t755 = -t433 / 0.2e1;
t751 = -t434 / 0.2e1;
t750 = -t435 / 0.2e1;
t749 = -t436 / 0.2e1;
t748 = t441 / 0.2e1;
t747 = t449 / 0.2e1;
t454 = t806 * t521;
t746 = -t454 / 0.2e1;
t744 = -t561 / 0.2e1;
t743 = t562 / 0.2e1;
t741 = -t562 / 0.2e1;
t468 = mrSges(5,2) * t525 - mrSges(5,3) * t631;
t740 = t468 / 0.2e1;
t470 = -mrSges(5,1) * t525 - mrSges(5,3) * t629;
t739 = -t470 / 0.2e1;
t664 = t524 * mrSges(5,2);
t667 = t520 * mrSges(5,1);
t483 = t664 + t667;
t738 = t483 / 0.2e1;
t737 = -t504 / 0.2e1;
t736 = -t520 / 0.2e1;
t735 = t520 / 0.2e1;
t733 = -t524 / 0.2e1;
t731 = -t525 / 0.2e1;
t730 = -t525 / 0.4e1;
t727 = pkin(4) * t520;
t725 = pkin(9) * t520;
t724 = pkin(9) * t524;
t511 = t525 * pkin(8);
t720 = t81 * mrSges(7,2);
t719 = t82 * mrSges(7,1);
t714 = mrSges(4,2) * t521;
t712 = mrSges(5,3) * t521;
t711 = mrSges(7,3) * t440;
t710 = mrSges(7,3) * t441;
t709 = Ifges(5,4) * t520;
t708 = Ifges(6,4) * t433;
t707 = Ifges(6,4) * t562;
t698 = pkin(4) * qJD(4);
t697 = t101 * mrSges(7,1);
t696 = t102 * mrSges(7,2);
t682 = t305 * mrSges(7,1);
t678 = t308 * mrSges(7,2);
t675 = t591 * mrSges(7,3);
t671 = t435 * mrSges(6,1);
t670 = t436 * mrSges(6,2);
t662 = t525 * Ifges(5,5);
t661 = t525 * Ifges(5,6);
t660 = t81 * t591;
t659 = t82 * t352;
t658 = -mrSges(4,1) + t806;
t656 = t186 * t303;
t655 = t303 * t440;
t654 = t592 * t441;
t395 = t449 * t610;
t33 = m(7) * (t123 * t158 + t157 * t831 + t395) + m(6) * (t240 * t269 + t268 * t593 + t395) + m(5) * (t376 * t412 + t377 * t413 + t395) + m(4) * (t449 * t521 + t450 * t525 - t636) * t635;
t653 = t33 * qJD(1);
t360 = t449 * t450;
t647 = t377 * t524;
t648 = t376 * t520;
t34 = m(7) * (t123 * t166 + t165 * t831 + t360) + m(6) * (t240 * t301 + t300 * t593 + t360) + m(5) * (t360 + (-t647 + t648) * t449);
t652 = t34 * qJD(1);
t651 = t352 * t440;
t650 = t591 * t441;
t640 = t440 * t592;
t639 = t441 * t303;
t638 = t449 * t520;
t637 = t472 * t520;
t429 = t521 * t812 - t661;
t632 = t520 * t429;
t581 = Ifges(5,1) * t524 - t709;
t550 = t581 * t521;
t431 = t550 - t662;
t628 = t524 * t431;
t473 = pkin(4) * t630 + t511;
t149 = -mrSges(7,1) * t441 - mrSges(7,2) * t440;
t618 = t149 * qJD(6);
t614 = pkin(4) * t629;
t613 = t727 / 0.2e1;
t612 = -Ifges(5,2) / 0.4e1 + Ifges(5,1) / 0.4e1;
t611 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t609 = mrSges(7,3) * t786;
t607 = t675 / 0.2e1;
t606 = -t673 / 0.2e1;
t605 = -t669 / 0.2e1;
t604 = -t668 / 0.2e1;
t603 = t433 * t741;
t602 = t849 * t747;
t371 = mrSges(6,1) * t562 - mrSges(6,2) * t561;
t317 = -t433 * mrSges(6,1) - t434 * mrSges(6,2);
t589 = t614 / 0.2e1;
t587 = m(6) * t589;
t585 = 0.2e1 * (m(7) / 0.4e1 + m(6) / 0.4e1) * t450;
t584 = t870 + t738 + t371 / 0.2e1;
t579 = -Ifges(6,1) * t434 + t708;
t578 = -Ifges(6,1) * t561 - t707;
t485 = t524 * Ifges(5,2) + t709;
t572 = -t702 + t706;
t571 = Ifges(5,5) * t520 + Ifges(5,6) * t524;
t181 = -mrSges(7,1) * t592 - mrSges(7,2) * t303;
t182 = t678 - t682;
t265 = -mrSges(7,2) * t521 + mrSges(7,3) * t305;
t318 = mrSges(6,1) * t434 - mrSges(6,2) * t433;
t319 = -t670 + t671;
t363 = pkin(5) * t435 + t473;
t455 = t483 * t521;
t456 = t483 * t525;
t469 = -mrSges(5,2) * t521 - mrSges(5,3) * t630;
t471 = mrSges(5,1) * t521 - mrSges(5,3) * t627;
t414 = t460 - t615;
t563 = t414 * t520 - t415 * t524;
t528 = (t318 / 0.2e1 + t181 / 0.2e1 + t455 / 0.2e1) * t450 + (t468 * t733 + t470 * t735 + t456 / 0.2e1 + t319 / 0.2e1 + t182 / 0.2e1) * t449 + (t450 * t510 + t376 * t423 + t377 * t424 + (t563 + t511) * t449) * t800 + (t230 * t300 + t231 * t301 + t234 * t593 + t235 * t240 + t449 * t473 + t450 * t472) * t798 + (t123 * t90 + t165 * t81 + t166 * t82 + t362 * t450 + t363 * t449 + t831 * t89) * t796 + t831 * t777 + t265 * t872 + t165 * t778 + t166 * t780 + t593 * t756 + t240 * t758 + t300 * t757 + t301 * t759 + t376 * t471 / 0.2e1 + t377 * t469 / 0.2e1;
t3 = t528 + (t157 * t761 + t158 * t765) * mrSges(7,3) + (t268 * t743 - t561 * t776) * mrSges(6,3) + (-t484 / 0.2e1 + t663 / 0.2e1 + (mrSges(4,1) / 0.2e1 - t806 / 0.2e1 + t760 + t784) * t521) * t635 + t834;
t168 = Ifges(7,4) * t308 + Ifges(7,2) * t305 + Ifges(7,6) * t521;
t170 = Ifges(7,1) * t308 + Ifges(7,4) * t305 + Ifges(7,5) * t521;
t296 = -t434 * Ifges(6,2) - t525 * Ifges(6,6) - t708;
t297 = -Ifges(6,4) * t436 - Ifges(6,2) * t435 + t521 * Ifges(6,6);
t299 = -Ifges(6,1) * t436 - Ifges(6,4) * t435 + Ifges(6,5) * t521;
t430 = Ifges(5,6) * t521 + t525 * t812;
t432 = Ifges(5,5) * t521 + t525 * t581;
t558 = Ifges(6,5) * t436 / 0.2e1 + Ifges(6,6) * t435 / 0.2e1;
t5 = -pkin(2) * t484 + t415 * t469 + t423 * t470 + t414 * t471 + t472 * t319 + t473 * t318 + t424 * t468 + t235 * t391 + t231 * t392 + t234 * t393 + t230 * t394 + t362 * t182 + t363 * t181 + t89 * t266 + t81 * t267 + t90 * t264 + t82 * t265 + t169 * t767 + t170 * t768 + t167 * t771 + t296 * t750 + t297 * t751 + t299 * t755 + t298 * t749 + (t628 / 0.2e1 - t632 / 0.2e1 + pkin(8) * t455 + t557 + t558 + t818 * t525) * t525 + t168 * t825 + (t432 * t732 + t430 * t736 + pkin(8) * t456 + Ifges(6,5) * t755 + Ifges(6,6) * t751 + Ifges(7,5) * t768 + Ifges(7,6) * t825 - t818 * t521 + (m(5) * pkin(8) ^ 2 + Ifges(4,1) - Ifges(4,2) - Ifges(5,3) - Ifges(6,3) - Ifges(7,3)) * t525) * t521 + m(5) * (t414 * t423 + t415 * t424) + m(6) * (t230 * t234 + t231 * t235 + t472 * t473) + m(7) * (t362 * t363 + t81 * t89 + t82 * t90);
t568 = t3 * qJD(1) + t5 * qJD(2);
t380 = -pkin(5) * t433 + t614;
t529 = (t123 * t774 + t838) * mrSges(7,3) - t782 * t858 + (-t647 / 0.2e1 + t648 / 0.2e1) * t712 + (t123 * t878 + t875 * t831) * t796 - t123 * t778 + t855 + t376 * t740 + t377 * t739 + (t798 * t837 + t759) * t593 + (t746 + t317 / 0.2e1 + t871 + t614 * t798 + t380 * t796) * t449 + (t753 * mrSges(6,3) + t798 * t864 - t757) * t240;
t7 = t529 + t854;
t8 = -t415 * t470 + t472 * t317 + t414 * t468 + t579 * t755 + t296 * t753 + t244 * t391 + t243 * t393 + t380 * t181 + t867 + t851 * t768 + t167 * t774 + t101 * t266 + t102 * t264 + (t303 * t82 - t592 * t81) * mrSges(7,3) + (t230 * t434 + t231 * t433) * mrSges(6,3) + m(7) * (t101 * t81 + t102 * t82 + t362 * t380) + m(6) * (t230 * t243 + t231 * t244) + (t429 * t733 + t431 * t736 - pkin(8) * t454 + t525 * t571 / 0.2e1 + (t485 * t735 + t856) * t521 + (m(6) * t472 + t318) * t726 + t563 * mrSges(5,3)) * t521 + t809 * t751 + t808 * t731 + t850 * t825;
t567 = t7 * qJD(1) + t8 * qJD(2);
t13 = t81 * t264 - t82 * t266 + t867 + t623 * t731 - (-t82 * mrSges(7,3) - t167 / 0.2e1 + t851 / 0.2e1) * t303 + (-t81 * mrSges(7,3) + t833 / 0.2e1 + t169 / 0.2e1) * t592;
t538 = (-t303 * t788 + t838) * mrSges(7,3) + t855 + t123 * t779 + t848 * t747;
t14 = t538 - t625;
t566 = t14 * qJD(1) + t13 * qJD(2);
t37 = t592 * t264 + t303 * t266 - t434 * t391 + t433 * t393 + m(7) * (t303 * t81 + t592 * t82) + m(6) * (t230 * t433 - t231 * t434);
t540 = (-t240 * t434 + t433 * t593) * t798 + (t123 * t592 + t303 * t831) * t796;
t553 = (t796 + t798) * t610;
t42 = t553 - t540;
t565 = -qJD(1) * t42 + qJD(2) * t37;
t549 = (t433 * t518 - t434 * t516) * t795;
t56 = t587 + (t380 / 0.4e1 - t655 / 0.4e1 - t654 / 0.4e1) * t802 - t549 / 0.2e1 + t848 + t317;
t425 = pkin(5) * t562 + t727;
t544 = (-t516 * t561 - t518 * t562) * t616;
t58 = t520 * t617 + (-t651 / 0.4e1 + t650 / 0.4e1 - t425 / 0.4e1) * t802 + t544 - t849 - t371;
t564 = qJD(2) * t56 - qJD(3) * t58;
t108 = 0.2e1 * mrSges(7,1) * t761 + 0.2e1 * t608;
t79 = 0.2e1 * t827 * mrSges(7,2) + 0.2e1 * t775;
t560 = qJD(2) * t79 - qJD(3) * t108;
t552 = t516 * t758 + t518 * t756;
t551 = t485 * t736 - t856;
t548 = (t234 * t518 + t235 * t516) * t798;
t530 = -(t593 / 0.2e1 + t782) * t669 + (t591 * t789 + t828 * t831) * mrSges(7,3) + pkin(4) * t638 * t798 + t425 * t449 * t796;
t12 = (-t664 / 0.2e1 - t667 / 0.2e1 + t584) * t449 + t530 + t853;
t373 = -Ifges(6,2) * t561 + t707;
t16 = t372 * t727 + t581 * t735 + t812 * t732 - pkin(3) * t483 + t578 * t743 + t373 * t741 + t866 + t852 * t761 + t220 * t844 + t551 + t814 * t744 + t863 * t828 + (m(6) * t727 + t371) * t504 + (m(7) * t418 + t218) * t425;
t527 = -t472 * t371 / 0.2e1 - t425 * t181 / 0.2e1 - t813 * t391 / 0.2e1 + t386 * t757 + (-t186 * t878 + t362 * t425 + t380 * t418) * t797 + t380 * t784 + t373 * t754 + t317 * t737 + t296 * t742 + pkin(3) * t746 + t814 * t434 / 0.4e1 + t809 * t561 / 0.4e1 + t810 * t525 / 0.4e1 - t850 * t591 / 0.4e1 - t863 * t592 / 0.4e1 + t862 * t303 + (t797 * t875 + t781) * t847 - t873;
t534 = (t440 * t89 + t441 * t90) * t796 + t234 * t793 - t235 * mrSges(6,2) / 0.2e1 + t423 * t794 + t424 * t792 + t440 * t777 + t265 * t748 - t558 + t583;
t542 = t386 * t864 + t837 * t813;
t2 = t527 + (0.3e1 / 0.4e1 * t662 + pkin(9) * t470 / 0.2e1 - t431 / 0.4e1 + (0.3e1 / 0.4e1 * t709 + t485 / 0.4e1 - t612 * t524 + (m(6) * t737 + t760) * pkin(4)) * t521) * t524 + Ifges(6,4) * t603 + pkin(4) * t548 + (t434 * t742 + t561 * t754) * Ifges(6,1) + (t659 / 0.2e1 + t102 * t765 + t660 / 0.2e1 + t101 * t761 + t847 * t825 + t656 / 0.2e1) * mrSges(7,3) + (-pkin(8) * t483 / 0.2e1 + (t514 / 0.2e1 + t512 / 0.2e1) * pkin(9) * mrSges(5,3) + (-t580 / 0.4e1 + t508 / 0.4e1 + t612 * t520) * t520 + t611) * t521 + (-(-t231 / 0.2e1 - t243 / 0.2e1) * t562 - (-t244 / 0.2e1 + t830) * t561 - t804) * mrSges(6,3) + (pkin(4) * t637 + t542) * t799 + t552 * pkin(4) + (-0.3e1 / 0.4e1 * t661 + pkin(9) * t740 - pkin(4) * t318 / 0.2e1 + t429 / 0.4e1) * t520 + t534;
t547 = t12 * qJD(1) - t2 * qJD(2) + t16 * qJD(3);
t532 = (t169 / 0.4e1 + t833 / 0.4e1) * t591 + (t609 + t832 / 0.4e1 + t222 / 0.4e1) * t592 - (mrSges(7,3) * t785 + t862) * t303 + t847 * t780 + t621 * t730 + t873;
t10 = t532 - t583;
t21 = t602 - t624;
t29 = t866 + (t852 / 0.2e1 - t220 / 0.2e1) * t352 + (t222 / 0.2e1 + t832 / 0.2e1) * t591;
t546 = t21 * qJD(1) + t10 * qJD(2) + t29 * qJD(3);
t531 = (t303 * t844 + t592 * t828) * mrSges(7,3) + (-t434 * t744 + t603) * mrSges(6,3) + (-t230 * t562 - t231 * t561 - t386 * t434 + t433 * t813) * t798 + (-t186 * t592 + t303 * t847 - t352 * t81 + t591 * t82) * t796 + t266 * t844 + t264 * t828 + t391 * t744 + t393 * t741;
t537 = t473 * t798 + t363 * t796 - t682 / 0.2e1 + t678 / 0.2e1 + t671 / 0.2e1 - t670 / 0.2e1;
t20 = t531 - t537;
t539 = (-t240 * t561 - t562 * t593) * t799 + (t123 * t591 - t352 * t831) * t797;
t46 = t585 + t539;
t49 = (t352 ^ 2 + t591 ^ 2) * mrSges(7,3) + (t561 ^ 2 + t562 ^ 2) * mrSges(6,3) + m(7) * (-t186 * t591 - t352 * t847) + m(6) * (-t386 * t561 - t562 * t813);
t545 = qJD(1) * t46 - qJD(2) * t20 - qJD(3) * t49;
t543 = t266 * t748 + t440 * t781 - t599;
t17 = (-t639 / 0.2e1 + t640 / 0.2e1) * mrSges(7,3) + (-t102 / 0.2e1 + t81 / 0.2e1) * mrSges(7,2) + (t101 / 0.2e1 + t82 / 0.2e1) * mrSges(7,1) + t543 + t599;
t27 = (t789 + t831 / 0.2e1) * mrSges(7,2) + (t788 + t872) * mrSges(7,1);
t32 = (t786 + t847 / 0.2e1) * mrSges(7,2) + (t785 - t186 / 0.2e1) * mrSges(7,1);
t541 = t27 * qJD(1) - t17 * qJD(2) + t32 * qJD(3) + t149 * qJD(4);
t515 = t525 ^ 2;
t513 = t521 ^ 2;
t479 = t513 * pkin(8) * t635;
t109 = t608 - t820 / 0.2e1;
t91 = m(6) * t613 + t544 + (t425 + t650 - t651) * t796;
t80 = -t840 / 0.2e1 + t775;
t67 = t549 / 0.2e1 + t587 + (t654 + t655 + t380) * t796;
t47 = t585 - t539;
t43 = t553 + t540;
t22 = t602 + t624;
t19 = t531 + t537;
t18 = -t719 / 0.2e1 - t720 / 0.2e1 + t697 / 0.2e1 - t696 / 0.2e1 - t543 + t599 + (-t639 + t640) * t791;
t15 = t538 + t625;
t11 = t449 * t584 + t638 * t794 + t664 * t747 + t530 - t853;
t9 = t532 + t583;
t6 = t529 - t854;
t4 = t157 * t606 + t158 * t607 + t268 * t604 + t269 * t605 + t528 - t484 * t635 + (t806 + t815) * t610 / 0.2e1 - t834;
t1 = -t527 + t628 / 0.4e1 + t101 * t606 + t102 * t607 + t318 * t613 + t244 * t605 - t632 / 0.4e1 - t712 * t859 / 0.2e1 - t468 * t725 / 0.2e1 - t485 * t629 / 0.2e1 + t611 * t521 + t592 * t609 + ((t504 * t629 + t637) * pkin(4) + t542) * t798 + t660 * t791 + t578 * t754 + t510 * t738 + t724 * t739 + t579 * t742 + t572 * t730 + t372 * t589 + t804 * mrSges(6,3) + (-t659 - t656) * mrSges(7,3) / 0.2e1 + (t580 / 0.2e1 - t812 / 0.4e1) * t631 + t559 * t525 + (t548 + t552) * pkin(4) + t524 * t550 / 0.4e1 + t534 + t837 * t604 + t669 * t830;
t23 = [t33 * qJD(2) + t34 * qJD(3), t4 * qJD(3) + t6 * qJD(4) + t43 * qJD(5) + t15 * qJD(6) + t653 + (t157 * t266 + t158 * t264 + t268 * t393 + t269 * t391 + t412 * t470 + t413 * t468 + ((-mrSges(4,1) * t525 - mrSges(3,1) + t714) * t522 + (-mrSges(3,2) + (t513 + t515) * mrSges(4,3) + (t181 + t318 + t455) * t521) * t526) * t517 + 0.2e1 * (t81 * t157 + t82 * t158 + t362 * t610) * t796 + 0.2e1 * (t230 * t268 + t231 * t269 + t472 * t610) * t798 + 0.2e1 * (t412 * t414 + t413 * t415 + t479) * t800 + m(4) * (t479 + (pkin(8) * t515 * t526 - pkin(2) * t522) * t517)) * qJD(2), t652 + t4 * qJD(2) + t11 * qJD(4) + t47 * qJD(5) + t22 * qJD(6) + ((t165 * t847 - t166 * t186 + t418 * t450) * t796 + (t300 * t813 + t301 * t386 + t450 * t504) * t798 + (-pkin(3) * t450 - t449 * t859) * t800) * t801 + (-t165 * t673 + t166 * t675 - t300 * t668 - t301 * t669 + (-mrSges(5,3) * t811 + mrSges(4,2)) * t449 + (t658 + t815) * t450) * qJD(3), t6 * qJD(2) + t11 * qJD(3) + (-t377 * mrSges(5,1) - mrSges(6,1) * t240 - t376 * mrSges(5,2) - mrSges(6,2) * t593 + t28) * qJD(4) + t874 + 0.2e1 * ((-t123 * t440 + t441 * t831) * t796 + (-t240 * t518 + t516 * t593) * t616) * qJD(4), qJD(2) * t43 + qJD(3) * t47, t15 * qJD(2) + t22 * qJD(3) + t28 * qJD(4) + t874; qJD(3) * t3 + qJD(4) * t7 - qJD(5) * t42 + qJD(6) * t14 - t653, qJD(3) * t5 + qJD(4) * t8 + qJD(5) * t37 + qJD(6) * t13, t1 * qJD(4) + t19 * qJD(5) + t9 * qJD(6) + ((-pkin(3) * t511 + pkin(9) * t807) * t800 + (t234 * t813 + t235 * t386 + t473 * t504) * t798 + (-t186 * t90 + t363 * t418 + t847 * t89) * t796) * t801 + t568 + ((Ifges(6,5) * t562 + Ifges(7,5) * t352 - Ifges(6,6) * t561 + Ifges(7,6) * t591 + t571) * t843 - Ifges(4,6) * t521 + t504 * t319 + t473 * t372 - pkin(3) * t456 + (pkin(8) * t658 + Ifges(4,5) + t551) * t525 + t418 * t182 + t386 * t392 + t363 * t218 - t471 * t725 - t89 * t673 + t222 * t767 + t220 * t771 + t170 * t761 + t432 * t735 + t299 * t743 + t297 * t744 + t374 * t749 + t373 * t750 + t469 * t724 + t430 * t732 + pkin(8) * t714 + t90 * t675 + t807 * mrSges(5,3) - t186 * t265 + t168 * t828 + t813 * t394 + t847 * t267 - t234 * t668 - t235 * t669) * qJD(3), t1 * qJD(3) + (-Ifges(5,5) * t631 - Ifges(5,6) * t629 + m(7) * (t101 * t440 + t102 * t441) + t303 * t710 - t592 * t711 + t697 - t696 - t244 * mrSges(6,2) + t243 * mrSges(6,1) - t414 * mrSges(5,2) - t415 * mrSges(5,1) + t808) * qJD(4) + t67 * qJD(5) + t18 * qJD(6) + (m(6) * (t243 * t518 + t244 * t516) + (t433 * t516 + t434 * t518) * mrSges(6,3)) * t698 + t567, qJD(3) * t19 + qJD(4) * t67 + qJD(6) * t80 + t565, t9 * qJD(3) + t18 * qJD(4) + t80 * qJD(5) + (t623 - t719 - t720) * qJD(6) + t566; -qJD(2) * t3 + qJD(4) * t12 - qJD(5) * t46 + qJD(6) * t21 - t652, -qJD(4) * t2 + qJD(5) * t20 + qJD(6) * t10 - t568, qJD(4) * t16 + qJD(5) * t49 + qJD(6) * t29 (m(7) * (t186 * t440 + t441 * t847) - t386 * mrSges(6,1) - t813 * mrSges(6,2) - t352 * t710 - t591 * t711 + t572 + t806 * pkin(9) + t810 + t877) * qJD(4) + t91 * qJD(5) + t879 + (m(6) * (-t386 * t518 + t516 * t813) + (-t516 * t562 + t518 * t561) * mrSges(6,3)) * t698 + t547, qJD(4) * t91 + qJD(6) * t109 - t545, t31 * qJD(4) + t109 * qJD(5) + t546 + t879; -qJD(2) * t7 - qJD(3) * t12 + qJD(6) * t27, qJD(3) * t2 - qJD(5) * t56 - qJD(6) * t17 - t567, qJD(5) * t58 + qJD(6) * t32 - t547, t618, -t564, t541 + t618; qJD(2) * t42 + qJD(3) * t46, -qJD(3) * t20 + qJD(4) * t56 - qJD(6) * t79 - t565, -qJD(4) * t58 + qJD(6) * t108 + t545, t564, 0, -t560; -t14 * qJD(2) - t21 * qJD(3) - t27 * qJD(4), -qJD(3) * t10 + qJD(4) * t17 + qJD(5) * t79 - t566, -qJD(4) * t32 - qJD(5) * t108 - t546, -t541, t560, 0;];
Cq  = t23;
