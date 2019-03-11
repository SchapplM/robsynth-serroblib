% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRPR2
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
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:05:35
% EndTime: 2019-03-08 23:05:58
% DurationCPUTime: 13.47s
% Computational Cost: add. (27354->675), mult. (61143->921), div. (0->0), fcn. (69875->12), ass. (0->389)
t803 = qJD(3) + qJD(4);
t463 = cos(qJ(3));
t461 = sin(qJ(3));
t672 = sin(qJ(4));
t551 = t672 * t461;
t674 = cos(qJ(4));
t420 = -t463 * t674 + t551;
t460 = cos(pkin(12));
t641 = t460 * mrSges(6,2);
t458 = sin(pkin(12));
t643 = t458 * mrSges(6,1);
t462 = cos(qJ(6));
t671 = sin(qJ(6));
t502 = t671 * t458 - t462 * t460;
t321 = t502 * t420;
t655 = t321 * mrSges(7,2);
t501 = t462 * t458 + t460 * t671;
t318 = t501 * t420;
t657 = t318 * mrSges(7,1);
t711 = m(7) / 0.2e1;
t713 = m(6) / 0.2e1;
t704 = -pkin(9) - pkin(8);
t431 = t704 * t463;
t729 = -t674 * t431 + t704 * t551;
t601 = t420 * t458;
t750 = -pkin(5) * t601 + t729;
t805 = (-t641 / 0.2e1 - t643 / 0.2e1) * t420 - t657 / 0.2e1 + t655 / 0.2e1 + t711 * t750 + t713 * t729;
t459 = sin(pkin(6));
t673 = sin(qJ(2));
t555 = t459 * t673;
t634 = cos(pkin(6));
t403 = t461 * t634 + t463 * t555;
t494 = t461 * t555 - t463 * t634;
t311 = t403 * t672 + t674 * t494;
t353 = mrSges(7,1) * t501 - mrSges(7,2) * t502;
t683 = t353 / 0.2e1;
t545 = t311 * t683;
t706 = mrSges(7,2) / 0.2e1;
t707 = -mrSges(7,1) / 0.2e1;
t754 = t502 * t311;
t755 = t501 * t311;
t773 = t706 * t754 + t707 * t755;
t793 = t545 + t773;
t804 = qJD(1) * t793;
t484 = t403 * t674 - t494 * t672;
t464 = cos(qJ(2));
t592 = t459 * t464;
t244 = -t458 * t484 - t460 * t592;
t245 = -t458 * t592 + t460 * t484;
t120 = t462 * t244 - t245 * t671;
t121 = t244 * t671 + t462 * t245;
t514 = t244 * t458 - t245 * t460;
t737 = t311 * t484;
t769 = m(6) * (t311 * t514 + t737) + m(7) * (t120 * t755 + t121 * t754 + t737);
t790 = t769 * qJD(1);
t802 = qJD(6) * t793 - t790;
t794 = t545 - t773;
t801 = t794 * qJD(6) + t790;
t552 = t674 * t461;
t381 = -t431 * t672 - t704 * t552;
t421 = -t463 * t672 - t552;
t444 = -pkin(3) * t463 - pkin(2);
t343 = pkin(4) * t420 + qJ(5) * t421 + t444;
t213 = t460 * t343 - t458 * t729;
t214 = t458 * t343 + t460 * t729;
t517 = -t458 * t213 + t214 * t460;
t800 = t484 * t381 + (t729 - t517) * t311;
t354 = mrSges(7,1) * t502 + mrSges(7,2) * t501;
t428 = -mrSges(6,1) * t460 + mrSges(6,2) * t458;
t799 = t428 + t354;
t495 = t501 * t421;
t679 = t501 / 0.2e1;
t680 = t502 / 0.2e1;
t732 = t502 * t421;
t230 = -mrSges(7,2) * t420 + mrSges(7,3) * t495;
t232 = mrSges(7,1) * t420 - mrSges(7,3) * t732;
t599 = t421 * t458;
t350 = -mrSges(6,2) * t420 + mrSges(6,3) * t599;
t598 = t421 * t460;
t352 = t420 * mrSges(6,1) + mrSges(6,3) * t598;
t675 = t460 / 0.2e1;
t503 = t350 * t675 - t458 * t352 / 0.2e1;
t714 = -m(6) / 0.2e1;
t747 = t230 * t680 + t232 * t679 + t517 * t714 - t503;
t798 = (t495 * t680 - t679 * t732) * mrSges(7,3) + t747 + t805;
t646 = t501 * mrSges(7,3);
t564 = -t646 / 0.2e1;
t647 = t502 * mrSges(7,3);
t566 = -t647 / 0.2e1;
t797 = t495 * t566 - t564 * t732 - t747 + t805;
t688 = t321 / 0.2e1;
t691 = t318 / 0.2e1;
t508 = Ifges(7,5) * t688 + Ifges(7,6) * t691;
t726 = Ifges(7,3) * t421 / 0.2e1 - t508;
t361 = -pkin(4) * t421 + qJ(5) * t420;
t759 = t381 * t458;
t225 = t460 * t361 + t759;
t600 = t420 * t460;
t526 = -t421 * pkin(5) + pkin(10) * t600;
t151 = t225 + t526;
t756 = t460 * t381;
t226 = t458 * t361 - t756;
t582 = pkin(10) * t601;
t179 = t582 + t226;
t85 = t462 * t151 - t179 * t671;
t86 = t151 * t671 + t462 * t179;
t796 = -t706 * t86 - t707 * t85 - t726;
t670 = pkin(3) * t461;
t344 = t361 + t670;
t219 = t460 * t344 + t759;
t136 = t219 + t526;
t220 = t458 * t344 - t756;
t171 = t582 + t220;
t80 = t462 * t136 - t171 * t671;
t81 = t136 * t671 + t462 * t171;
t795 = -t706 * t81 - t707 * t80 - t726;
t524 = -t641 - t643;
t338 = t524 * t420;
t792 = t338 / 0.2e1 - t503;
t141 = Ifges(7,1) * t321 + Ifges(7,4) * t318 - t421 * Ifges(7,5);
t661 = Ifges(6,2) * t458;
t664 = Ifges(6,4) * t460;
t249 = -Ifges(6,6) * t421 + (t661 - t664) * t420;
t665 = Ifges(6,4) * t458;
t250 = -Ifges(6,5) * t421 + (-Ifges(6,1) * t460 + t665) * t420;
t662 = Ifges(7,4) * t501;
t356 = -Ifges(7,2) * t502 + t662;
t411 = Ifges(7,4) * t502;
t358 = Ifges(7,1) * t501 - t411;
t676 = t458 / 0.2e1;
t678 = -t421 / 0.2e1;
t703 = Ifges(7,4) * t688 + Ifges(7,2) * t691 + Ifges(7,6) * t678;
t488 = t141 * t679 + t249 * t675 + t250 * t676 + t356 * t691 + t358 * t688 - t502 * t703 + (Ifges(6,2) * t460 + t665) * t601 / 0.2e1 - (Ifges(6,1) * t458 + t664) * t600 / 0.2e1 + Ifges(5,6) * t421 - Ifges(5,5) * t420 + (Ifges(6,5) * t458 + Ifges(7,5) * t501 + Ifges(6,6) * t460 - Ifges(7,6) * t502) * t678;
t752 = t729 * t428;
t761 = t729 * mrSges(5,1);
t762 = t381 * mrSges(5,2);
t781 = t750 * t354;
t791 = t488 + t752 + t762 - t761 + t781;
t789 = t762 / 0.2e1 + t781 / 0.2e1;
t387 = t420 * t592;
t330 = -t460 * t387 + t458 * t555;
t618 = t330 * t460;
t328 = t458 * t387 + t460 * t555;
t620 = t328 * t458;
t788 = -t387 * mrSges(5,2) / 0.2e1 - (-t620 / 0.2e1 + t618 / 0.2e1) * mrSges(6,3);
t541 = -t592 / 0.2e1;
t767 = -t732 / 0.2e1;
t786 = mrSges(7,1) * t767;
t453 = t458 ^ 2;
t455 = t460 ^ 2;
t585 = t453 + t455;
t785 = -mrSges(6,3) * t585 + mrSges(5,2);
t286 = -pkin(5) * t599 + t381;
t784 = t286 * t750;
t580 = t674 * pkin(3);
t443 = -t580 - pkin(4);
t668 = t460 * pkin(5);
t426 = t443 - t668;
t783 = t426 * t750;
t439 = -pkin(4) - t668;
t782 = t439 * t750;
t586 = t585 * qJ(5);
t779 = t752 / 0.2e1 - t761 / 0.2e1;
t175 = t655 - t657;
t778 = (t175 / 0.2e1 + t792) * t311;
t663 = Ifges(7,4) * t732;
t140 = Ifges(7,2) * t495 + t420 * Ifges(7,6) + t663;
t309 = Ifges(7,4) * t495;
t142 = Ifges(7,1) * t732 + t420 * Ifges(7,5) + t309;
t176 = -mrSges(7,1) * t495 + mrSges(7,2) * t732;
t229 = mrSges(7,2) * t421 + mrSges(7,3) * t318;
t231 = -mrSges(7,1) * t421 - mrSges(7,3) * t321;
t339 = t524 * t421;
t349 = mrSges(6,2) * t421 + mrSges(6,3) * t601;
t351 = -mrSges(6,1) * t421 + mrSges(6,3) * t600;
t362 = -mrSges(5,1) * t421 - mrSges(5,2) * t420;
t640 = t460 * Ifges(6,5);
t642 = t458 * Ifges(6,6);
t741 = -t495 / 0.2e1;
t133 = pkin(5) * t420 + pkin(10) * t598 + t213;
t166 = pkin(10) * t599 + t214;
t77 = t462 * t133 - t166 * t671;
t78 = t133 * t671 + t462 * t166;
t777 = t750 * t176 + t729 * t339 + ((t455 * Ifges(6,1) / 0.2e1 - Ifges(6,3) + Ifges(5,1) - Ifges(5,2) - Ifges(7,3) + (-t664 + t661 / 0.2e1) * t458) * t421 + t508 + (Ifges(5,4) - t640 + t642) * t420) * t420 + (Ifges(7,5) * t767 + Ifges(7,6) * t741 - t460 * t250 / 0.2e1 + t249 * t676 + (t640 / 0.2e1 - t642 / 0.2e1 - Ifges(5,4)) * t421) * t421 + t213 * t351 + t214 * t349 + t286 * t175 + t381 * t338 + t444 * t362 + t77 * t231 + t78 * t229 + t140 * t691 + t495 * t703 + t142 * t688 + t732 * t141 / 0.2e1;
t579 = t672 * pkin(3);
t438 = t579 + qJ(5);
t405 = (-pkin(10) - t438) * t458;
t449 = t460 * pkin(10);
t406 = t438 * t460 + t449;
t329 = t462 * t405 - t406 * t671;
t331 = t405 * t671 + t462 * t406;
t774 = t329 * t755 + t331 * t754;
t427 = (-pkin(10) - qJ(5)) * t458;
t429 = qJ(5) * t460 + t449;
t376 = t462 * t427 - t429 * t671;
t377 = t427 * t671 + t462 * t429;
t772 = t376 * t755 + t377 * t754 + t439 * t484;
t698 = t232 / 0.2e1;
t699 = t230 / 0.2e1;
t771 = t755 * t698 + t754 * t699 + t120 * t231 / 0.2e1 + t121 * t229 / 0.2e1 + t244 * t351 / 0.2e1 + t245 * t349 / 0.2e1;
t770 = t286 * t484 + t311 * t750 + t754 * t78 + t755 * t77;
t694 = t311 / 0.2e1;
t766 = m(5) * t444;
t765 = pkin(4) * t729;
t764 = mrSges(7,2) * t741;
t758 = t381 * t672;
t757 = t443 * t729;
t613 = t729 * t381;
t742 = m(6) + m(7);
t751 = t742 * t579;
t749 = t461 ^ 2 + t463 ^ 2;
t731 = (t339 / 0.2e1 + t176 / 0.2e1) * t484;
t748 = t362 * t541 + t731;
t746 = -mrSges(5,1) + t799;
t386 = t421 * t592;
t745 = (mrSges(5,1) / 0.2e1 - t799 / 0.2e1) * t386;
t534 = t585 * t438;
t740 = m(6) * t534;
t739 = pkin(4) * t484;
t430 = t461 * mrSges(4,1) + t463 * mrSges(4,2);
t736 = t430 * t541;
t583 = m(6) / 0.4e1 + m(7) / 0.4e1;
t733 = t484 * t583;
t730 = t176 + t339;
t447 = t453 * mrSges(6,3);
t448 = t455 * mrSges(6,3);
t728 = t447 + t448;
t727 = -t447 / 0.2e1 - t448 / 0.2e1;
t629 = t220 * t460;
t630 = t219 * t458;
t516 = t629 - t630;
t723 = -mrSges(4,1) * t463 + mrSges(4,2) * t461;
t558 = t428 / 0.2e1 - mrSges(5,1) / 0.2e1;
t693 = t484 / 0.2e1;
t722 = t693 - t484 / 0.2e1;
t92 = 0.2e1 * t764 + 0.2e1 * t786;
t720 = -qJD(2) * t92 + t353 * t803;
t355 = -Ifges(7,2) * t501 - t411;
t357 = -Ifges(7,1) * t502 - t662;
t530 = -t501 * t356 / 0.2e1 + t357 * t679 - (t355 + t358) * t502 / 0.2e1;
t172 = t462 * t328 - t330 * t671;
t173 = t328 * t671 + t462 * t330;
t718 = t172 * t564 + t173 * t566 + t745 - t788;
t717 = (t172 * t679 + t173 * t680) * mrSges(7,3) + t788;
t591 = t460 * t349;
t593 = t458 * t351;
t596 = t439 * t175;
t614 = t377 * t229;
t615 = t376 * t231;
t636 = t81 * t502;
t637 = t80 * t501;
t669 = pkin(4) * t338;
t716 = (qJ(5) * t516 - t765) * t713 + (t376 * t80 + t377 * t81 + t782) * t711 - t669 / 0.2e1 + t615 / 0.2e1 + t614 / 0.2e1 + t596 / 0.2e1 + (t591 / 0.2e1 - t593 / 0.2e1) * qJ(5) + (t629 / 0.2e1 - t630 / 0.2e1) * mrSges(6,3) + (-t636 / 0.2e1 - t637 / 0.2e1) * mrSges(7,3) + t789;
t715 = 2 * qJD(4);
t712 = -m(7) / 0.2e1;
t710 = m(5) * pkin(3);
t709 = m(6) * pkin(3);
t705 = -mrSges(7,3) / 0.2e1;
t654 = t732 * mrSges(7,1);
t656 = t495 * mrSges(7,2);
t174 = t654 + t656;
t702 = t174 / 0.2e1;
t700 = -t230 / 0.2e1;
t696 = t484 / 0.4e1;
t695 = -t311 / 0.2e1;
t686 = -t331 / 0.2e1;
t682 = -t377 / 0.2e1;
t527 = t674 * t671;
t553 = t462 * t674;
t392 = (-t458 * t553 - t460 * t527) * pkin(3);
t681 = -t392 / 0.2e1;
t659 = t484 * mrSges(5,1);
t658 = t311 * mrSges(5,2);
t645 = t420 * mrSges(5,3);
t644 = t421 * mrSges(5,3);
t68 = -t120 * t501 - t121 * t502;
t635 = -t514 * t713 + t68 * t711;
t632 = t120 * t495;
t631 = t121 * t732;
t628 = t225 * t458;
t627 = t226 * t460;
t622 = t484 * t354;
t621 = t484 * t428;
t224 = t311 * t386;
t619 = t329 * t231;
t617 = t331 * t229;
t557 = t459 ^ 2 * t673;
t34 = m(6) * (t244 * t328 + t245 * t330 - t224) + m(7) * (t120 * t172 + t121 * t173 - t224) + m(5) * (-t387 * t484 - t464 * t557 - t224) + m(4) * (-t557 + (t403 * t463 + t461 * t494) * t459) * t464;
t616 = t34 * qJD(1);
t611 = t381 * t386;
t597 = t426 * t175;
t595 = t439 * t353;
t594 = t443 * t338;
t589 = -t376 * t501 - t377 * t502;
t588 = Ifges(7,5) * t495 - Ifges(7,6) * t732;
t587 = -Ifges(7,5) * t502 - Ifges(7,6) * t501;
t577 = mrSges(6,3) * t628;
t576 = mrSges(6,3) * t627;
t571 = t438 * t593;
t570 = t438 * t591;
t567 = t647 / 0.2e1;
t565 = t646 / 0.2e1;
t556 = t458 * t674;
t554 = t460 * t674;
t537 = t355 / 0.4e1 + t358 / 0.4e1;
t536 = -t356 / 0.4e1 + t357 / 0.4e1;
t535 = t585 * t311;
t529 = t728 + (t501 ^ 2 + t502 ^ 2) * mrSges(7,3);
t525 = -0.2e1 * t583 * t386;
t467 = m(5) * t670 * t541 + (t219 * t244 + t220 * t245 + t800) * t713 + (t120 * t80 + t121 * t81 + t770) * t711 + t771;
t512 = t618 - t620;
t472 = (-t386 * t443 + t438 * t512) * t713 + (t172 * t329 + t173 * t331 - t386 * t426) * t711 + (t386 * t674 - t387 * t672) * t710 / 0.2e1 + t736;
t3 = t175 * t694 + t311 * t792 + t722 * t644 + t467 - t472 + t717 + t736 - t745 + t748;
t363 = mrSges(5,1) * t420 - mrSges(5,2) * t421;
t4 = t670 * t766 + m(6) * (t213 * t219 + t214 * t220 + t613) + m(7) * (t77 * t80 + t78 * t81 + t784) + (Ifges(4,4) * t463 + (Ifges(4,1) - Ifges(4,2)) * t461) * t463 - pkin(2) * t430 + t220 * t350 + t219 * t352 + t81 * t230 + t80 * t232 + (-Ifges(4,4) * t461 + pkin(3) * t363) * t461 + t777;
t523 = t3 * qJD(1) + t4 * qJD(2);
t466 = t778 + (t225 * t244 + t226 * t245 + t800) * t713 + (t120 * t85 + t121 * t86 + t770) * t711 + t748 + t771;
t478 = (pkin(4) * t386 + qJ(5) * t512) * t714 + (t172 * t376 + t173 * t377 - t386 * t439) * t712;
t7 = -(-t354 / 0.2e1 - t558) * t386 + t466 + t478 + t717;
t8 = m(6) * (t213 * t225 + t214 * t226 + t613) + m(7) * (t77 * t85 + t78 * t86 + t784) + t226 * t350 + t225 * t352 + t86 * t230 + t85 * t232 + t777;
t522 = t7 * qJD(1) + t8 * qJD(2);
t521 = -t501 * t77 - t502 * t78;
t177 = -Ifges(7,2) * t732 + t309;
t178 = Ifges(7,1) * t495 - t663;
t15 = -t78 * t232 + t420 * t588 / 0.2e1 + t286 * t174 + t77 * t230 + (t178 / 0.2e1 - t140 / 0.2e1 - t78 * mrSges(7,3)) * t732 + (-t77 * mrSges(7,3) + t142 / 0.2e1 + t177 / 0.2e1) * t495;
t493 = t120 * t700 + t121 * t698 + t174 * t695;
t505 = t172 * mrSges(7,1) / 0.2e1 - t173 * mrSges(7,2) / 0.2e1;
t16 = (t632 / 0.2e1 + t631 / 0.2e1) * mrSges(7,3) + t493 + t505;
t519 = -t16 * qJD(1) + t15 * qJD(2);
t32 = -t732 * t232 + t495 * t230 + m(7) * (t495 * t78 - t732 * t77) + (t460 * t352 + t458 * t350 + m(6) * (t213 * t460 + t214 * t458)) * t421;
t483 = (-t120 * t732 + t121 * t495) * t712 + m(6) * (t244 * t460 + t245 * t458) * t678;
t45 = t525 + t483;
t518 = qJD(1) * t45 - qJD(2) * t32;
t515 = t627 - t628;
t511 = -t329 * t501 - t331 * t502;
t393 = (-t458 * t527 + t460 * t553) * pkin(3);
t504 = mrSges(7,1) * t681 + t393 * t706;
t499 = t585 * t674;
t282 = t311 * t579;
t468 = (t282 - t438 * t535 + (-t244 * t556 + t245 * t554) * pkin(3)) * t713 + (t120 * t392 + t121 * t393 + t282 + t774) * t711 + mrSges(5,2) * t694 + t354 * t693 + t755 * t564 + t754 * t566 + (t426 * t711 + t443 * t713 + t558) * t484 + t727 * t311;
t469 = -t714 * t739 + t772 * t712 + t659 / 0.2e1 - t622 / 0.2e1 - t621 / 0.2e1 - t658 / 0.2e1 + t755 * t565 + t754 * t567 - (t586 * t714 + t727) * t311;
t14 = t468 + t469;
t465 = (t515 * t438 + t757) * t714 + (t286 * t579 + t329 * t85 + t331 * t86 + t392 * t77 + t393 * t78 + t783) * t712 - t619 / 0.2e1 - t617 / 0.2e1 + t232 * t681 + t393 * t700 - t597 / 0.2e1 - t594 / 0.2e1 + t577 / 0.2e1 - t576 / 0.2e1 + t571 / 0.2e1 - t570 / 0.2e1 + t85 * t565 + t86 * t567 - t730 * t579 / 0.2e1 + ((-t213 * t556 + t214 * t554 + t758) * t714 + t352 * t556 / 0.2e1 - t350 * t554 / 0.2e1) * pkin(3) - t779 - t789;
t5 = t465 + t716 + t779;
t482 = t392 * t646 + t393 * t647 - t746 * t579 + t580 * t785;
t56 = -m(7) * (t329 * t392 + t331 * t393 + t426 * t579) - (t438 * t499 + t443 * t672) * t709 + t482;
t498 = t14 * qJD(1) - t5 * qJD(2) - t56 * qJD(3);
t481 = -(t142 / 0.4e1 + t177 / 0.4e1) * t502 - (-t178 / 0.4e1 + t140 / 0.4e1) * t501 + t286 * t683 + t420 * t587 / 0.4e1;
t471 = (mrSges(7,3) * t686 + t536) * t732 + (t329 * t705 + t537) * t495 + t329 * t699 + t232 * t686 + t426 * t702 + t481;
t10 = t471 - t795;
t293 = t426 * t353;
t48 = -t293 - t530;
t497 = t10 * qJD(2) - t48 * qJD(3) + t804;
t490 = m(7) * (-t329 * t732 + t331 * t495 + t521);
t20 = -t490 / 0.2e1 + t798;
t49 = 0.2e1 * (t696 - t68 / 0.4e1) * m(7) + 0.2e1 * (t696 + t514 / 0.4e1) * m(6);
t82 = -m(7) * t511 - t529 - t740;
t496 = -qJD(1) * t49 - qJD(2) * t20 - qJD(3) * t82;
t470 = (mrSges(7,3) * t682 + t536) * t732 + (t376 * t705 + t537) * t495 + t376 * t699 + t232 * t682 + t439 * t702 + t481;
t12 = t470 - t796;
t485 = t293 / 0.2e1 + t595 / 0.2e1 + t530 + (t564 + t565) * t331;
t37 = t485 + t504;
t54 = t530 + t595;
t492 = t12 * qJD(2) + t37 * qJD(3) + t54 * qJD(4) + t804;
t489 = m(7) * (-t376 * t732 + t377 * t495 + t521);
t23 = -t489 / 0.2e1 + t798;
t52 = t635 - 0.2e1 * t733;
t477 = (t534 + t586) * t713 + (t511 + t589) * t711 + t529;
t61 = t477 - t751 / 0.2e1;
t95 = m(6) * t586 + m(7) * t589 + t529;
t491 = qJD(1) * t52 - qJD(2) * t23 + qJD(3) * t61 + qJD(4) * t95;
t348 = t353 * qJD(5);
t347 = t353 * qJD(6);
t93 = t656 / 0.2e1 + t654 / 0.2e1 + t786 + t764;
t60 = t477 + t751 / 0.2e1;
t51 = t693 * t742 + t635;
t50 = t635 + 0.2e1 * t733;
t46 = t525 - t483;
t36 = t485 - t504;
t24 = t489 / 0.2e1 + t797;
t21 = t490 / 0.2e1 + t797;
t17 = -t493 + t505 + (t631 + t632) * t705;
t13 = t468 - t469;
t11 = t470 + t796;
t9 = t471 + t795;
t6 = t466 - t478 + t718;
t2 = (t722 * t421 + (t695 + t694) * t420) * mrSges(5,3) + t731 + t467 + (-t430 / 0.2e1 - t362 / 0.2e1) * t592 + t472 + t718 + t778;
t1 = t558 * t729 - t465 + t488 + t716;
t18 = [t34 * qJD(2) + t769 * t803, t2 * qJD(3) + t6 * qJD(4) + t46 * qJD(5) + t17 * qJD(6) + t616 + (m(6) * (t213 * t328 + t214 * t330 - t611) + t172 * t232 + t173 * t230 + m(7) * (t172 * t77 + t173 * t78) + t330 * t350 + t328 * t352 + m(5) * (-t387 * t729 - t611) + t387 * t645 + m(4) * (pkin(8) * t464 * t749 - t673 * pkin(2)) * t459 + (-mrSges(3,1) + t363 + t723 + t766) * t555 - (m(7) * t286 - t644 + t730) * t386 + (mrSges(4,3) * t749 - mrSges(3,2)) * t592) * qJD(2), t2 * qJD(2) + (-t755 * t646 - t754 * t647 + t622 + m(7) * t774 + t621 + t658 - t659 + t494 * mrSges(4,2) - t403 * mrSges(4,1) + (m(6) * t443 + m(7) * t426 - t674 * t710) * t484 - (t672 * t710 + t728 + t740) * t311) * qJD(3) + t13 * qJD(4) + t50 * qJD(5) + t801, t6 * qJD(2) + t13 * qJD(3) + t51 * qJD(5) + ((-qJ(5) * t535 - t739) * t713 + t772 * t711) * t715 + ((-t501 * t755 - t502 * t754) * mrSges(7,3) + t746 * t484 + t785 * t311) * qJD(4) + t801, qJD(2) * t46 + qJD(3) * t50 + qJD(4) * t51, t17 * qJD(2) + (-mrSges(7,1) * t121 - mrSges(7,2) * t120) * qJD(6) + t803 * t794; qJD(3) * t3 + qJD(4) * t7 - qJD(5) * t45 - qJD(6) * t16 - t616, qJD(3) * t4 + qJD(4) * t8 + qJD(5) * t32 + qJD(6) * t15 (t723 * pkin(8) + t516 * mrSges(6,3) + t619 + t617 + t597 + t594 + (-t636 - t637) * mrSges(7,3) + t579 * t644 + t580 * t645 + Ifges(4,5) * t463 - Ifges(4,6) * t461 - t571 + t570 + m(7) * (t329 * t80 + t331 * t81 + t783) + m(6) * (t438 * t516 + t757) + (-t674 * t729 - t758) * t710 + t791) * qJD(3) + t1 * qJD(4) + t21 * qJD(5) + t9 * qJD(6) + t523, t1 * qJD(3) + (-t646 * t85 - t647 * t86 + t576 - t577 + t596 + t614 + t615 - t669 + (t591 - t593) * qJ(5) + t791) * qJD(4) + t24 * qJD(5) + t11 * qJD(6) + ((qJ(5) * t515 - t765) * t713 + (t376 * t85 + t377 * t86 + t782) * t711) * t715 + t522, qJD(3) * t21 + qJD(4) * t24 + qJD(6) * t93 - t518, t9 * qJD(3) + t11 * qJD(4) + t93 * qJD(5) + (-mrSges(7,1) * t78 - mrSges(7,2) * t77 + t588) * qJD(6) + t519; -qJD(2) * t3 + qJD(4) * t14 - qJD(5) * t49 + t802, -qJD(4) * t5 - qJD(5) * t20 + qJD(6) * t10 - t523, -qJD(4) * t56 - qJD(5) * t82 - qJD(6) * t48 (m(7) * (t376 * t392 + t377 * t393 + t439 * t579) + (-pkin(4) * t672 + qJ(5) * t499) * t709 - t482) * qJD(4) + t60 * qJD(5) + t36 * qJD(6) + t498, qJD(4) * t60 + t496, t36 * qJD(4) + (-mrSges(7,1) * t331 - mrSges(7,2) * t329 + t587) * qJD(6) + t497; -qJD(2) * t7 - qJD(3) * t14 + qJD(5) * t52 + t802, qJD(3) * t5 - qJD(5) * t23 + qJD(6) * t12 - t522, qJD(5) * t61 + qJD(6) * t37 - t498, qJD(5) * t95 + qJD(6) * t54, t491 (-mrSges(7,1) * t377 - mrSges(7,2) * t376 + t587) * qJD(6) + t492; qJD(2) * t45 + qJD(3) * t49 - qJD(4) * t52, qJD(3) * t20 + qJD(4) * t23 - qJD(6) * t92 + t518, -qJD(4) * t61 + t347 - t496, t347 - t491, 0, t720; t16 * qJD(2) - t793 * t803, -qJD(3) * t10 - qJD(4) * t12 + qJD(5) * t92 - t519, -qJD(4) * t37 - t348 - t497, -t348 - t492, -t720, 0;];
Cq  = t18;
