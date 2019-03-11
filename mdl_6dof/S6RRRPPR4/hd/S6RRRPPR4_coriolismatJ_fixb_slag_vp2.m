% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:56
% EndTime: 2019-03-09 15:32:30
% DurationCPUTime: 20.65s
% Computational Cost: add. (27001->811), mult. (56854->1074), div. (0->0), fcn. (59560->8), ass. (0->425)
t496 = sin(qJ(3));
t649 = sin(pkin(10));
t592 = t649 * t496;
t650 = cos(pkin(10));
t730 = cos(qJ(3));
t442 = -t650 * t730 + t592;
t495 = sin(qJ(6));
t498 = cos(qJ(6));
t594 = t650 * t496;
t515 = t649 * t730 + t594;
t330 = t442 * t498 - t495 * t515;
t543 = t442 * t495 + t498 * t515;
t635 = Ifges(7,5) * t330 - Ifges(7,6) * t543;
t609 = t730 * qJ(4);
t625 = t730 * pkin(8);
t460 = t625 + t609;
t716 = -qJ(4) - pkin(8);
t359 = -t649 * t460 + t716 * t594;
t288 = t515 * pkin(9) + t359;
t812 = t650 * t460 + t716 * t592;
t813 = t442 * pkin(9) + t812;
t136 = t288 * t498 + t495 * t813;
t860 = -t288 * t495 + t498 * t813;
t876 = -t860 * mrSges(7,1) + t136 * mrSges(7,2);
t22 = -t635 - t876;
t878 = t22 * qJD(6);
t175 = -mrSges(7,1) * t330 + mrSges(7,2) * t543;
t483 = -pkin(3) * t730 - pkin(2);
t526 = -qJ(5) * t515 + t483;
t780 = -pkin(4) - pkin(5);
t237 = t442 * t780 - t526;
t877 = -m(7) * t237 - t175;
t873 = t136 * t495 + t498 * t860;
t737 = t496 / 0.2e1;
t843 = t330 * mrSges(7,3);
t824 = t543 * mrSges(7,1);
t174 = t330 * mrSges(7,2) + t824;
t871 = t237 * t174;
t855 = Ifges(5,6) - Ifges(6,6);
t856 = Ifges(6,4) + Ifges(5,5);
t869 = -t442 * t856 - t515 * t855 + t635;
t320 = Ifges(7,4) * t330;
t179 = Ifges(7,1) * t543 + t320;
t558 = Ifges(7,2) * t543 - t320;
t868 = -t558 / 0.4e1 + t179 / 0.4e1;
t497 = sin(qJ(2));
t416 = t442 * t497;
t417 = t515 * t497;
t279 = t416 * t495 + t417 * t498;
t544 = -t416 * t498 + t417 * t495;
t822 = t544 * mrSges(7,1);
t126 = t279 * mrSges(7,2) + t822;
t493 = t497 * pkin(7);
t643 = t496 * t497;
t457 = pkin(3) * t643 + t493;
t566 = -qJ(5) * t416 - t457;
t171 = t417 * t780 + t566;
t825 = Ifges(7,6) * t544;
t848 = Ifges(7,5) * t279;
t636 = t848 - t825;
t499 = cos(qJ(2));
t732 = t499 / 0.2e1;
t867 = t171 * t126 + t636 * t732;
t707 = Ifges(7,4) * t544;
t116 = Ifges(7,2) * t279 + t499 * Ifges(7,6) + t707;
t268 = Ifges(7,4) * t279;
t118 = Ifges(7,1) * t544 + t499 * Ifges(7,5) + t268;
t130 = Ifges(7,1) * t279 - t707;
t671 = t543 * Ifges(7,4);
t177 = t330 * Ifges(7,2) + t671;
t178 = Ifges(7,1) * t330 - t671;
t821 = t544 * mrSges(7,3);
t235 = mrSges(7,1) * t499 - t821;
t559 = Ifges(7,2) * t544 - t268;
t784 = -mrSges(7,3) / 0.2e1;
t866 = (t178 / 0.4e1 + t860 * t784 - t177 / 0.4e1) * t544 + (-t116 / 0.4e1 + t130 / 0.4e1) * t543 + t237 * t126 / 0.2e1 + t171 * t174 / 0.2e1 - t860 * t235 / 0.2e1 + (-t559 / 0.4e1 + t118 / 0.4e1) * t330;
t731 = t499 / 0.4e1;
t315 = t442 * pkin(4) + t526;
t346 = mrSges(6,1) * t442 - mrSges(6,3) * t515;
t801 = m(6) * t315 + t346;
t347 = mrSges(5,1) * t442 + mrSges(5,2) * t515;
t863 = m(5) * t483 + t347;
t711 = Ifges(4,4) * t496;
t533 = -Ifges(4,2) * t730 - t711;
t491 = Ifges(4,4) * t730;
t462 = Ifges(4,1) * t496 + t491;
t626 = t730 / 0.2e1;
t568 = t462 * t626;
t862 = t533 * t737 + t568;
t722 = t497 * pkin(8);
t459 = -pkin(2) * t499 - pkin(1) - t722;
t446 = t730 * t459;
t540 = -t497 * t609 + t446;
t335 = (-pkin(7) * t496 - pkin(3)) * t499 + t540;
t610 = t499 * t730;
t575 = pkin(7) * t610;
t355 = t575 + (-qJ(4) * t497 + t459) * t496;
t593 = t649 * t355;
t188 = t335 * t650 - t593;
t173 = t499 * pkin(4) - t188;
t726 = t416 * pkin(9);
t120 = t499 * pkin(5) + t173 + t726;
t596 = t650 * t355;
t189 = t649 * t335 + t596;
t170 = -qJ(5) * t499 + t189;
t725 = t417 * pkin(9);
t131 = t170 + t725;
t57 = t120 * t498 - t131 * t495;
t642 = t496 * t499;
t628 = pkin(7) * t642;
t354 = t540 - t628;
t192 = t354 * t649 + t596;
t146 = t192 + t725;
t193 = t650 * t354 - t593;
t147 = t193 - t726;
t71 = t146 * t495 + t147 * t498;
t861 = -t71 / 0.2e1 - t57 / 0.2e1;
t858 = 0.2e1 * mrSges(7,2);
t857 = Ifges(5,1) + Ifges(6,1);
t565 = t495 * mrSges(7,1) + t498 * mrSges(7,2);
t58 = t120 * t495 + t131 * t498;
t70 = t146 * t498 - t147 * t495;
t854 = t70 / 0.2e1 - t58 / 0.2e1;
t535 = -t848 / 0.2e1 + t825 / 0.2e1;
t770 = -t279 / 0.2e1;
t853 = -t330 / 0.2e1;
t852 = t330 / 0.2e1;
t611 = t497 * t730;
t571 = -t611 / 0.2e1;
t615 = -t822 / 0.2e1;
t613 = -t824 / 0.2e1;
t831 = mrSges(5,3) + mrSges(6,2);
t845 = t279 * mrSges(7,3);
t632 = t565 * t732;
t233 = -mrSges(7,2) * t499 + t845;
t640 = t498 * t233;
t644 = t495 * t235;
t637 = -t644 / 0.2e1 + t640 / 0.2e1;
t839 = t632 - t637;
t838 = qJD(6) * t565;
t837 = -mrSges(5,1) / 0.2e1;
t833 = t543 / 0.2e1;
t769 = t544 / 0.2e1;
t823 = t543 * mrSges(7,3);
t404 = Ifges(5,4) * t417;
t705 = Ifges(6,5) * t417;
t820 = -t416 * t857 - t499 * t856 - t404 + t705;
t418 = t515 * t499;
t419 = t499 * t442;
t819 = t856 * t497 - t857 * t419 + (-Ifges(5,4) + Ifges(6,5)) * t418;
t438 = Ifges(5,4) * t442;
t658 = t442 * Ifges(6,5);
t818 = t515 * t857 - t438 + t658;
t490 = t499 * mrSges(6,3);
t667 = t417 * mrSges(6,2);
t364 = -t490 - t667;
t666 = t417 * mrSges(5,3);
t365 = mrSges(5,2) * t499 - t666;
t817 = -t364 - t365;
t668 = t416 * mrSges(5,3);
t367 = -mrSges(5,1) * t499 + t668;
t669 = t416 * mrSges(6,2);
t368 = mrSges(6,1) * t499 - t669;
t816 = t367 - t368;
t815 = -Ifges(4,2) * t496 + t491;
t531 = t730 * mrSges(4,1) - t496 * mrSges(4,2);
t810 = t416 * t855 - t417 * t856;
t809 = -t368 / 0.2e1 + t367 / 0.2e1;
t808 = -t365 / 0.2e1 - t364 / 0.2e1;
t806 = -Ifges(5,2) * t515 - t438 + t818;
t805 = Ifges(5,2) * t416 - t404 + t820;
t437 = Ifges(6,5) * t515;
t348 = t442 * Ifges(6,3) + t437;
t709 = Ifges(5,4) * t515;
t804 = -t442 * t857 + t348 + t437 - t709;
t401 = Ifges(6,5) * t416;
t269 = -t499 * Ifges(6,6) + Ifges(6,3) * t417 - t401;
t710 = Ifges(5,4) * t416;
t803 = -t417 * t857 + t269 - t401 + t710;
t606 = t649 * pkin(3);
t477 = t606 + qJ(5);
t802 = m(6) * t477 + mrSges(6,3);
t607 = t650 * pkin(3);
t482 = -t607 - pkin(4);
t788 = m(5) * pkin(3);
t799 = m(6) * t482 - t650 * t788 - mrSges(5,1) - mrSges(6,1);
t798 = t649 * t788 - mrSges(5,2) + t802;
t796 = 0.2e1 * t416;
t795 = -0.2e1 * t515;
t794 = -m(5) / 0.2e1;
t793 = -m(6) / 0.2e1;
t792 = m(6) / 0.2e1;
t791 = -m(7) / 0.2e1;
t790 = m(7) / 0.2e1;
t789 = m(7) / 0.4e1;
t787 = -mrSges(6,1) / 0.2e1;
t786 = -mrSges(5,2) / 0.2e1;
t785 = mrSges(6,3) / 0.2e1;
t778 = -t136 / 0.2e1;
t576 = pkin(3) * t611;
t249 = -t416 * pkin(4) + t417 * qJ(5) + t576;
t187 = -t416 * pkin(5) + t249;
t776 = t187 / 0.2e1;
t775 = -t233 / 0.2e1;
t729 = pkin(3) * t496;
t323 = pkin(4) * t515 + t442 * qJ(5) + t729;
t242 = pkin(5) * t515 + t323;
t773 = t242 / 0.2e1;
t772 = t279 / 0.2e1;
t280 = t418 * t498 + t419 * t495;
t771 = t280 / 0.2e1;
t283 = t418 * t495 - t419 * t498;
t768 = t283 / 0.2e1;
t762 = -t543 / 0.2e1;
t474 = -pkin(5) + t482;
t393 = t474 * t498 - t477 * t495;
t756 = t393 / 0.2e1;
t394 = t474 * t495 + t477 * t498;
t755 = t394 / 0.2e1;
t754 = -t416 / 0.2e1;
t752 = -t417 / 0.2e1;
t751 = t417 / 0.2e1;
t749 = -t418 / 0.2e1;
t748 = t418 / 0.2e1;
t747 = -t419 / 0.2e1;
t746 = -t442 / 0.2e1;
t745 = t442 / 0.2e1;
t743 = t515 / 0.2e1;
t742 = -t515 / 0.2e1;
t740 = -t495 / 0.2e1;
t739 = t495 / 0.2e1;
t736 = -t497 / 0.2e1;
t735 = t497 / 0.2e1;
t734 = t498 / 0.2e1;
t733 = -t499 / 0.2e1;
t727 = pkin(8) * t496;
t723 = t497 * pkin(2);
t721 = t57 * mrSges(7,2);
t720 = t58 * mrSges(7,1);
t719 = t70 * mrSges(7,1);
t718 = t71 * mrSges(7,2);
t712 = Ifges(3,4) * t497;
t708 = Ifges(6,4) * t419;
t706 = Ifges(5,5) * t419;
t703 = Ifges(7,5) * t283;
t700 = Ifges(6,2) * t497;
t699 = Ifges(5,6) * t418;
t698 = Ifges(6,6) * t418;
t696 = Ifges(7,6) * t280;
t694 = Ifges(4,3) * t497;
t693 = Ifges(5,3) * t497;
t692 = Ifges(7,3) * t497;
t689 = t136 * mrSges(7,3);
t683 = t280 * mrSges(7,1);
t679 = t283 * mrSges(7,2);
t117 = Ifges(7,4) * t283 + Ifges(7,2) * t280 - Ifges(7,6) * t497;
t119 = Ifges(7,1) * t283 + Ifges(7,4) * t280 - Ifges(7,5) * t497;
t127 = -mrSges(7,1) * t279 + mrSges(7,2) * t544;
t128 = t679 - t683;
t458 = (pkin(7) + t729) * t499;
t222 = t418 * pkin(4) + t419 * qJ(5) + t458;
t172 = t418 * pkin(5) + t222;
t465 = -pkin(8) * t499 + t723;
t407 = pkin(7) * t643 + t730 * t465;
t343 = t497 * pkin(3) - t499 * t609 + t407;
t408 = -pkin(7) * t611 + t496 * t465;
t356 = -qJ(4) * t642 + t408;
t191 = t649 * t343 + t650 * t356;
t180 = t497 * qJ(5) + t191;
t551 = -t343 * t650 + t649 * t356;
t181 = -t497 * pkin(4) + t551;
t221 = pkin(4) * t417 - t566;
t234 = mrSges(7,2) * t497 + mrSges(7,3) * t280;
t236 = -mrSges(7,1) * t497 - mrSges(7,3) * t283;
t270 = -Ifges(6,5) * t419 + Ifges(6,6) * t497 + Ifges(6,3) * t418;
t271 = -Ifges(5,2) * t417 - t499 * Ifges(5,6) - t710;
t272 = -Ifges(5,4) * t419 - Ifges(5,2) * t418 + Ifges(5,6) * t497;
t290 = mrSges(6,1) * t417 + mrSges(6,3) * t416;
t291 = mrSges(5,1) * t417 - mrSges(5,2) * t416;
t661 = t419 * mrSges(6,3);
t664 = t418 * mrSges(6,1);
t292 = t661 + t664;
t663 = t419 * mrSges(5,2);
t665 = t418 * mrSges(5,1);
t293 = -t663 + t665;
t363 = -mrSges(6,2) * t418 + mrSges(6,3) * t497;
t366 = -mrSges(5,2) * t497 - mrSges(5,3) * t418;
t369 = mrSges(5,1) * t497 + mrSges(5,3) * t419;
t653 = t497 * mrSges(6,1);
t662 = t419 * mrSges(6,2);
t370 = -t653 - t662;
t391 = t446 - t628;
t392 = t496 * t459 + t575;
t411 = -Ifges(4,6) * t499 + t497 * t815;
t412 = Ifges(4,6) * t497 + t499 * t815;
t539 = Ifges(4,1) * t730 - t711;
t524 = t539 * t497;
t413 = -Ifges(4,5) * t499 + t524;
t414 = Ifges(4,5) * t497 + t499 * t539;
t530 = t496 * mrSges(4,1) + mrSges(4,2) * t730;
t440 = t530 * t499;
t623 = mrSges(4,3) * t643;
t453 = mrSges(4,2) * t499 - t623;
t454 = -t497 * mrSges(4,2) - mrSges(4,3) * t642;
t572 = mrSges(4,3) * t611;
t455 = -t499 * mrSges(4,1) - t572;
t456 = t497 * mrSges(4,1) - mrSges(4,3) * t610;
t518 = t530 * t493;
t536 = Ifges(4,5) * t730 - Ifges(4,6) * t496;
t525 = t499 * t536;
t569 = t610 / 0.2e1;
t570 = t611 / 0.2e1;
t598 = -t642 / 0.2e1;
t600 = -t643 / 0.2e1;
t121 = t419 * pkin(9) + t497 * t780 + t551;
t133 = pkin(9) * t418 + t180;
t64 = t121 * t498 - t133 * t495;
t641 = t497 * t499;
t65 = t121 * t495 + t133 * t498;
t3 = m(6) * (t170 * t180 + t173 * t181 + t221 * t222) + m(7) * (-t171 * t172 + t57 * t64 + t58 * t65) + m(4) * (pkin(7) ^ 2 * t641 + t391 * t407 + t392 * t408) + t499 * t518 + (t698 + t700 - t708 + t693 - t699 - t706 + t525 + t694) * t733 + (Ifges(7,5) * t544 + Ifges(7,6) * t279 + t712 + (Ifges(3,2) + Ifges(7,3)) * t499) * t736 - pkin(1) * (t497 * mrSges(3,1) + mrSges(3,2) * t499) + t408 * t453 + t392 * t454 + t407 * t455 + t391 * t456 + t457 * t293 + t458 * t291 + (t497 * t536 - t712 - t855 * t417 - t856 * t416 + (Ifges(3,1) - Ifges(6,2) - Ifges(4,3) - Ifges(5,3)) * t499) * t735 + t180 * t364 + t191 * t365 + t189 * t366 + t181 * t368 + t188 * t369 + t173 * t370 + t170 * t363 + t222 * t290 + t221 * t292 + t65 * t233 + t58 * t234 + t64 * t235 + t57 * t236 + t171 * t128 - t172 * t127 + (0.2e1 * Ifges(3,4) * t499 - t692 + t696 + t703 + (Ifges(3,1) - Ifges(3,2)) * t497) * t732 + t118 * t768 + t119 * t769 + t116 * t771 + t117 * t772 + t271 * t749 + t270 * t751 + t272 * t752 + t269 * t748 + t440 * t493 + t819 * t754 + t820 * t747 + t413 * t569 + t414 * t570 + t411 * t598 + t412 * t600 + m(5) * (-t188 * t551 + t189 * t191 + t457 * t458) - t551 * t367;
t678 = t3 * qJD(1);
t532 = -Ifges(4,5) * t496 - Ifges(4,6) * t730;
t556 = -Ifges(6,3) * t416 - t705;
t587 = -t416 * mrSges(6,1) + t417 * mrSges(6,3);
t591 = -t416 * mrSges(5,1) - t417 * mrSges(5,2);
t4 = (pkin(7) * t531 - t862) * t497 ^ 2 + m(7) * (-t171 * t187 + t57 * t70 + t58 * t71) - t532 * t641 / 0.2e1 + t57 * t845 + t291 * t576 - t173 * t667 - t867 + t416 * t271 / 0.2e1 + t249 * t290 + t71 * t233 + t70 * t235 - t187 * t127 + (-t130 + t116) * t769 + (m(5) * t576 + t591) * t457 + (-t455 - t572) * t392 + (t623 + t453) * t391 + t118 * t770 + t559 * t772 + t556 * t751 + (m(6) * t249 + t587) * t221 + t189 * t668 + t170 * t669 + (-m(5) * t188 + m(6) * t173 - t816) * t192 + (m(5) * t189 + m(6) * t170 - t817) * t193 + t810 * t733 + t411 * t571 + t413 * t600 + t188 * t666 + t58 * t821 + t803 * t754 + t805 * t752;
t670 = t4 * qJD(1);
t660 = t442 * mrSges(6,2);
t659 = t442 * mrSges(5,3);
t657 = t515 * mrSges(6,2);
t656 = t515 * mrSges(5,3);
t8 = t57 * t233 - t58 * t235 + (-t58 * mrSges(7,3) + t130 / 0.2e1 - t116 / 0.2e1) * t544 + (-t57 * mrSges(7,3) + t118 / 0.2e1 - t559 / 0.2e1) * t279 + t867;
t651 = t8 * qJD(1);
t18 = -t279 * t233 + t544 * t235 + t817 * t417 + t816 * t416 + m(7) * (-t279 * t58 + t544 * t57) + m(6) * (-t170 * t417 - t173 * t416) + m(5) * (t188 * t416 - t189 * t417);
t648 = qJD(1) * t18;
t23 = (-t127 + t290) * t416 + (-t364 - t640 + t644) * t499 + m(7) * (-t171 * t416 + (t495 * t57 - t498 * t58) * t499) + m(6) * (-t170 * t499 + t221 * t416);
t647 = qJD(1) * t23;
t94 = t394 * mrSges(7,1) + mrSges(7,2) * t393;
t646 = qJD(6) * t94;
t639 = t498 * t279;
t638 = t498 * t330;
t631 = t788 / 0.2e1;
t630 = t789 + t792;
t629 = t789 + m(6) / 0.4e1;
t620 = mrSges(5,3) / 0.2e1 + mrSges(6,2) / 0.2e1;
t618 = mrSges(7,3) * t739;
t616 = -t845 / 0.2e1;
t614 = -t821 / 0.2e1;
t612 = t660 / 0.2e1;
t601 = t543 * t740;
t590 = mrSges(5,1) * t515 - t442 * mrSges(5,2);
t586 = mrSges(6,1) * t515 + t442 * mrSges(6,3);
t579 = t359 * t416 - t417 * t812;
t574 = mrSges(5,3) * t607;
t573 = mrSges(5,3) * t606;
t555 = Ifges(6,3) * t515 - t658;
t349 = -t442 * Ifges(5,2) + t709;
t545 = -t192 * t359 + t193 * t812;
t500 = t854 * t823 + t861 * t843 + t866 + t453 * t727 / 0.2e1 - t291 * t729 / 0.2e1 + (-t349 / 0.4e1 + t804 / 0.4e1) * t416 + t531 * t723 / 0.2e1 + (t815 + 0.2e1 * t462) * t643 / 0.4e1 + (t496 ^ 2 + t730 ^ 2) * mrSges(4,3) * t722 / 0.2e1 + t525 / 0.4e1 + t455 * t625 / 0.2e1 - t483 * t591 / 0.2e1 - t221 * t586 / 0.2e1 - t315 * t587 / 0.2e1 - t457 * t590 / 0.2e1 + (-t556 / 0.4e1 + t805 / 0.4e1) * t442 - t188 * t659 / 0.2e1 + (-t555 / 0.4e1 + t806 / 0.4e1) * t417 + t496 * t411 / 0.4e1 - t249 * t346 / 0.2e1 - t323 * t290 / 0.2e1 + (pkin(3) * t347 + t533) * t571 - (t524 + t413) * t730 / 0.4e1 - t689 * t770 + t127 * t773 + t136 * t775 + t175 * t776 + ((t457 * t496 + t483 * t611) * pkin(3) + t545) * t794 + (t221 * t323 + t249 * t315 + t545) * t793 - t518 / 0.2e1 + t831 * (t192 * t742 + t193 * t745 - t359 * t751 + t754 * t812) + t189 * t656 / 0.2e1 + t170 * t657 / 0.2e1 + (-t171 * t242 - t187 * t237 + (t57 + t71) * t860 + (-t70 + t58) * t136) * t791 + (t170 * t793 + t189 * t794 + t808) * t359 + t868 * t279 + t869 * t731 + t173 * t612 + (t173 * t793 - t188 * t794 + t809) * t812 + (t271 / 0.4e1 - t803 / 0.4e1) * t515;
t509 = t703 / 0.2e1 + t696 / 0.2e1 - t692 / 0.2e1 + t64 * mrSges(7,1) / 0.2e1 - t65 * mrSges(7,2) / 0.2e1;
t501 = t366 * t606 / 0.2e1 + t369 * t607 / 0.2e1 + t698 / 0.2e1 + t700 / 0.2e1 - t699 / 0.2e1 + t477 * t363 / 0.2e1 + t482 * t370 / 0.2e1 + t407 * mrSges(4,1) / 0.2e1 - t408 * mrSges(4,2) / 0.2e1 + t180 * t785 + t191 * t786 + t181 * t787 + (t393 * t64 + t394 * t65) * t790 + (t180 * t477 + t181 * t482) * t792 + t234 * t755 + t236 * t756 - t708 / 0.2e1 - t706 / 0.2e1 + t693 / 0.2e1 + t694 / 0.2e1 + Ifges(4,5) * t569 - t509 + Ifges(4,6) * t598 + (t191 * t649 - t551 * t650) * t631 + t551 * t837;
t1 = t501 + t500;
t7 = t315 * t586 - t871 + t483 * t590 + t179 * t853 + t558 * t852 + t349 * t742 + t555 * t745 + t539 * t737 - pkin(2) * t530 + t815 * t626 + t806 * t746 + t804 * t743 + (t177 - t178) * t833 + t863 * t729 + t801 * t323 + t862 + t877 * t242;
t552 = -t1 * qJD(1) + t7 * qJD(2);
t502 = (-t279 * t852 + t544 * t762) * mrSges(7,3) + (t417 * t620 + t808) * t442 - (t416 * t620 + t809) * t515 + m(5) * (-t188 * t515 - t189 * t442 + t579) / 0.2e1 + (-t170 * t442 + t173 * t515 + t579) * t792 + (-t136 * t544 - t279 * t860 - t330 * t58 + t543 * t57) * t790 + t235 * t833 + t233 * t853;
t510 = t458 * t794 + t222 * t793 + t172 * t791 - t683 / 0.2e1 + t679 / 0.2e1;
t10 = -(t786 + t785) * t419 + (t837 + t787) * t418 + t502 + t510;
t26 = (-t330 ^ 2 - t543 ^ 2) * mrSges(7,3) + m(7) * (-t136 * t543 - t330 * t860) + (t442 ^ 2 + t515 ^ 2) * t831 + (m(6) + m(5)) * (-t359 * t515 - t442 * t812);
t550 = -qJD(1) * t10 - qJD(2) * t26;
t503 = -(-t127 / 0.2e1 + t290 / 0.2e1) * t515 + (-t175 / 0.2e1 + t346 / 0.2e1) * t416 + (t612 + (-t638 / 0.2e1 + t601) * mrSges(7,3)) * t499 + (-t221 * t515 + t315 * t416 - t499 * t812) * t792 + (t171 * t515 - t237 * t416 - t499 * t873) * t790;
t508 = t181 * t792 + (t495 * t65 + t498 * t64) * t790 - t662 / 0.2e1 + t234 * t739 - t653 / 0.2e1 + t236 * t734;
t14 = t503 - t508;
t48 = (t801 + t877) * t515;
t549 = -qJD(1) * t14 + qJD(2) * t48;
t507 = (-t416 * t482 - t417 * t477) * t792 + (-t279 * t394 + t393 * t544) * t790 + (t416 * t650 - t417 * t649) * t631;
t511 = m(7) * t776 + t249 * t792 + t570 * t788;
t32 = t507 - t511 - t587 - t126 - t591;
t506 = (-t442 * t477 + t482 * t515) * t792 + (-t330 * t394 + t393 * t543) * t790 + (-t442 * t649 - t515 * t650) * t631;
t514 = m(7) * t773 + t323 * t792 + t496 * t631;
t36 = t506 - t514 - t586 - t174 - t590;
t548 = qJD(1) * t32 + qJD(2) * t36;
t66 = t770 * t858 + 0.2e1 * t615;
t79 = t853 * t858 + 0.2e1 * t613;
t547 = qJD(1) * t66 + qJD(2) * t79;
t527 = m(7) * (-t330 * t495 + t498 * t543);
t100 = -t527 / 0.2e1 + t630 * t795;
t528 = m(7) * (-t279 * t495 + t498 * t544);
t73 = -t528 / 0.2e1 + t630 * t796;
t542 = qJD(1) * t73 + qJD(2) * t100;
t529 = m(7) * t873;
t523 = t873 * t790;
t19 = t871 + (t178 / 0.2e1 - t177 / 0.2e1) * t543 + (t179 / 0.2e1 - t558 / 0.2e1) * t330;
t302 = t543 * t618;
t45 = mrSges(7,3) * t601 + t302;
t504 = (t689 / 0.2e1 + t868) * t279 + t233 * t778 + t635 * t731 + t866;
t5 = t504 - t509;
t522 = t5 * qJD(1) + t19 * qJD(2) + t45 * qJD(5);
t513 = t235 * t755 + t393 * t775 - t535;
t11 = (t279 * t756 + t544 * t755) * mrSges(7,3) + t861 * mrSges(7,2) + t854 * mrSges(7,1) + t513 + t535;
t21 = (t778 + t136 / 0.2e1) * mrSges(7,2);
t521 = -t11 * qJD(1) + t21 * qJD(2) + t94 * qJD(3);
t162 = m(7) * (-t393 * t495 + t394 * t498) + t565 + t802;
t505 = -t490 + ((-qJ(5) - t477) * t499 + t189) * t792 + ((-t394 * t499 + t58) * t498 + (t393 * t499 - t57) * t495) * t790 - t839;
t512 = t192 * t793 + (t495 * t71 + t498 * t70) * t791;
t17 = (t544 * t740 - t639 / 0.2e1) * mrSges(7,3) + t505 + t512;
t27 = ((t852 + t853) * t498 + (t833 + t762) * t495) * mrSges(7,3) + t523 - t529 / 0.2e1;
t520 = t17 * qJD(1) - t27 * qJD(2) + t162 * qJD(3);
t29 = (t279 * t734 + t544 * t739) * mrSges(7,3) + t839;
t519 = t29 * qJD(1) - t45 * qJD(2) - qJD(3) * t565;
t99 = t527 / 0.2e1 + m(6) * t743 + t629 * t795;
t80 = t613 + t824 / 0.2e1;
t72 = t528 / 0.2e1 + m(6) * t754 + t629 * t796;
t67 = t615 + t822 / 0.2e1;
t44 = t45 * qJD(6);
t43 = t506 + t514;
t40 = t507 + t511;
t30 = t495 * t614 + t498 * t616 + t632 + t637;
t28 = -t660 + t302 + t529 / 0.2e1 + t523 + 0.2e1 * t812 * t792 + (t638 / 0.2e1 + t543 * t739 + t498 * t852) * mrSges(7,3);
t16 = t544 * t618 - t639 * t784 + t505 - t512 - t667;
t13 = t503 + t508;
t12 = t721 / 0.2e1 + t720 / 0.2e1 + t394 * t614 + t393 * t616 - t718 / 0.2e1 + t719 / 0.2e1 - t513 + t535;
t9 = t502 + t665 / 0.2e1 - t663 / 0.2e1 + t661 / 0.2e1 + t664 / 0.2e1 - t510;
t6 = t504 + t509;
t2 = t501 - t500;
t15 = [qJD(2) * t3 + qJD(3) * t4 + qJD(4) * t18 + qJD(5) * t23 + qJD(6) * t8, t678 + (t863 * t458 - t456 * t727 + t65 * t843 + t117 * t852 + (Ifges(7,5) * t543 + Ifges(7,6) * t330 + t532) * t736 - t191 * t659 - t180 * t660 - Ifges(3,6) * t497 + t483 * t293 - pkin(2) * t440 + t860 * t234 + (-t442 * t855 + t515 * t856) * t735 + t315 * t292 + t237 * t128 - t172 * t175 + m(7) * (-t136 * t64 - t172 * t237 + t65 * t860) - t136 * t236 + ((-m(4) * pkin(2) - mrSges(3,1) - t531) * pkin(7) + t568 + Ifges(3,5)) * t499 + t179 * t768 + t177 * t771 + t270 * t745 + t272 * t746 + t348 * t748 + t349 * t749 + t414 * t737 + mrSges(3,2) * t493 + t181 * t657 + t818 * t747 + t819 * t743 + t454 * t625 + t412 * t626 + (m(4) * pkin(8) + mrSges(4,3)) * (-t407 * t496 + t730 * t408) + t551 * t656 - (m(5) * t551 + m(6) * t181 - t369 + t370) * t359 + (m(5) * t191 + m(6) * t180 + t363 + t366) * t812 - t533 * t598 - t64 * t823 + t119 * t833 + t801 * t222) * qJD(2) + t2 * qJD(3) + t9 * qJD(4) + t13 * qJD(5) + t6 * qJD(6), t670 + t2 * qJD(2) + (m(7) * (t393 * t70 + t394 * t71) - Ifges(4,6) * t611 + t416 * t573 + t417 * t574 - t391 * mrSges(4,2) - t392 * mrSges(4,1) + t718 - t719 + t636 - Ifges(4,5) * t643 + t477 * t669 - t482 * t667 + t393 * t845 + t394 * t821 + t798 * t193 + t799 * t192 + t810) * qJD(3) + t40 * qJD(4) + t16 * qJD(5) + t12 * qJD(6), qJD(2) * t9 + qJD(3) * t40 + qJD(5) * t72 + qJD(6) * t67 + t648, qJD(2) * t13 + qJD(3) * t16 + qJD(4) * t72 + qJD(6) * t30 + t647, t651 + t6 * qJD(2) + t12 * qJD(3) + t67 * qJD(4) + t30 * qJD(5) + (t636 - t720 - t721) * qJD(6); -qJD(3) * t1 + qJD(4) * t10 + qJD(5) * t14 + qJD(6) * t5 - t678, qJD(3) * t7 + qJD(4) * t26 - qJD(5) * t48 + qJD(6) * t19 (m(7) * (t136 * t394 + t393 * t860) + t536 + t442 * t574 - t515 * t573 - t477 * t657 - t482 * t660 + t394 * t823 + t393 * t843 + t798 * t359 + t799 * t812 - t531 * pkin(8) + t869 + t876) * qJD(3) + t43 * qJD(4) + t28 * qJD(5) + t878 + t552, qJD(3) * t43 + qJD(5) * t99 + qJD(6) * t80 - t550, qJD(3) * t28 + qJD(4) * t99 + t44 - t549, t22 * qJD(3) + t80 * qJD(4) + t522 - t878; qJD(2) * t1 + qJD(4) * t32 + qJD(5) * t17 - qJD(6) * t11 - t670, qJD(4) * t36 - qJD(5) * t27 + qJD(6) * t21 - t552, qJD(5) * t162 + t646, t548, t520, t521 - t646; -qJD(2) * t10 - qJD(3) * t32 + qJD(5) * t73 + qJD(6) * t66 - t648, -qJD(3) * t36 + qJD(5) * t100 + qJD(6) * t79 + t550, -t548, 0, t542, t547; -qJD(2) * t14 - qJD(3) * t17 - qJD(4) * t73 - qJD(6) * t29 - t647, qJD(3) * t27 - qJD(4) * t100 + t44 + t549, -t520 + t838, -t542, 0, -t519 - t838; -qJD(2) * t5 + qJD(3) * t11 - qJD(4) * t66 + qJD(5) * t29 - t651, -qJD(3) * t21 - qJD(4) * t79 - t522, -qJD(5) * t565 - t521, -t547, t519, 0;];
Cq  = t15;
