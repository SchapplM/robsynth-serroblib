% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:40
% EndTime: 2019-03-09 15:57:22
% DurationCPUTime: 23.97s
% Computational Cost: add. (27340->900), mult. (55691->1207), div. (0->0), fcn. (57632->8), ass. (0->432)
t518 = sin(qJ(6));
t521 = cos(qJ(6));
t517 = sin(pkin(10));
t519 = sin(qJ(3));
t522 = cos(qJ(3));
t685 = cos(pkin(10));
t610 = t685 * t522;
t547 = t517 * t519 + t610;
t611 = t685 * t519;
t548 = -t517 * t522 + t611;
t562 = -t518 * t547 + t521 * t548;
t603 = -t518 * t548 - t521 * t547;
t653 = Ifges(7,5) * t603 - Ifges(7,6) * t562;
t748 = pkin(8) * t519;
t463 = -t519 * qJ(5) + t748;
t466 = (pkin(8) - qJ(5)) * t522;
t336 = -t685 * t463 + t517 * t466;
t261 = -pkin(9) * t548 - t336;
t839 = t517 * t463 + t685 * t466;
t854 = -pkin(9) * t547 + t839;
t142 = t261 * t521 - t518 * t854;
t862 = t261 * t518 + t521 * t854;
t892 = t862 * mrSges(7,1) + t142 * mrSges(7,2);
t895 = t653 - t892;
t896 = qJD(6) * t895;
t520 = sin(qJ(2));
t747 = pkin(8) * t520;
t523 = cos(qJ(2));
t752 = pkin(2) * t523;
t462 = -pkin(1) - t747 - t752;
t663 = t519 * t523;
t651 = pkin(7) * t663 - t522 * t462;
t662 = t520 * t522;
t329 = qJ(5) * t662 - t651;
t516 = t523 * pkin(3);
t289 = pkin(4) * t523 - t329 + t516;
t667 = t462 * t519;
t749 = pkin(7) * t522;
t348 = t667 + (-qJ(4) + t749) * t523;
t664 = t519 * t520;
t487 = qJ(5) * t664;
t307 = t487 + t348;
t162 = t685 * t289 - t307 * t517;
t401 = -t517 * t664 - t520 * t610;
t744 = pkin(9) * t401;
t121 = pkin(5) * t523 + t162 + t744;
t163 = t517 * t289 + t685 * t307;
t402 = t517 * t662 - t520 * t611;
t743 = pkin(9) * t402;
t128 = t163 - t743;
t58 = t121 * t521 - t128 * t518;
t661 = t522 * t523;
t379 = pkin(7) * t661 + t667;
t330 = t379 + t487;
t183 = -t329 * t517 + t685 * t330;
t151 = t183 - t743;
t184 = t685 * t329 + t517 * t330;
t152 = t184 - t744;
t78 = t151 * t518 + t152 * t521;
t886 = t78 + t58;
t510 = t522 * qJ(4);
t488 = t520 * t510;
t524 = -pkin(3) - pkin(4);
t641 = t524 * t519;
t599 = -pkin(7) + t641;
t345 = t520 * t599 + t488;
t244 = pkin(5) * t402 + t345;
t253 = t401 * t521 + t402 * t518;
t856 = t253 * mrSges(7,1);
t256 = -t401 * t518 + t521 * t402;
t870 = t256 * mrSges(7,2);
t589 = t870 + t856;
t894 = t244 * t589;
t508 = t519 * qJ(4);
t573 = t522 * pkin(3) + t508;
t461 = -pkin(2) - t573;
t435 = t522 * pkin(4) - t461;
t326 = pkin(5) * t547 + t435;
t846 = t562 * mrSges(7,1);
t869 = t603 * mrSges(7,2);
t588 = -t869 - t846;
t893 = t326 * t588;
t293 = Ifges(7,4) * t603;
t170 = -Ifges(7,2) * t562 + t293;
t173 = Ifges(7,1) * t562 + t293;
t891 = -t170 / 0.4e1 - t173 / 0.4e1;
t715 = t253 * Ifges(7,4);
t124 = -Ifges(7,2) * t256 + t523 * Ifges(7,6) - t715;
t242 = Ifges(7,4) * t256;
t126 = -Ifges(7,1) * t253 + t523 * Ifges(7,5) - t242;
t580 = -Ifges(7,2) * t253 + t242;
t585 = Ifges(7,1) * t256 - t715;
t211 = mrSges(7,1) * t523 + mrSges(7,3) * t253;
t798 = -t211 / 0.2e1;
t890 = t862 * t798 + (-t124 / 0.4e1 - t585 / 0.4e1) * t562 - t244 * t588 / 0.2e1 - t326 * t589 / 0.2e1 + (-t580 / 0.4e1 + t126 / 0.4e1) * t603;
t887 = m(7) * t326;
t551 = t521 * t517 + t518 * t685;
t400 = t551 * t523;
t550 = t518 * t517 - t521 * t685;
t403 = t550 * t523;
t553 = t400 * mrSges(7,1) / 0.2e1 - t403 * mrSges(7,2) / 0.2e1;
t209 = -mrSges(7,2) * t523 - mrSges(7,3) * t256;
t769 = -t550 / 0.2e1;
t660 = t209 * t769 + t551 * t798;
t884 = t553 - t660;
t704 = t562 * Ifges(7,4);
t171 = Ifges(7,2) * t603 + t704;
t172 = Ifges(7,1) * t603 - t704;
t883 = t171 / 0.4e1 - t172 / 0.4e1;
t882 = -t142 * t551 - t550 * t862;
t868 = t603 * mrSges(7,3);
t880 = t868 / 0.2e1;
t459 = -qJ(4) * t517 + t685 * t524;
t449 = -pkin(5) + t459;
t460 = qJ(4) * t685 + t517 * t524;
t309 = t449 * t521 - t460 * t518;
t310 = t449 * t518 + t460 * t521;
t698 = t401 * mrSges(6,3);
t350 = mrSges(6,1) * t523 + t698;
t59 = t121 * t518 + t128 * t521;
t615 = t523 * t685;
t598 = -t615 / 0.2e1;
t697 = t402 * mrSges(6,3);
t349 = -mrSges(6,2) * t523 - t697;
t613 = t685 * t349;
t614 = t685 * t163;
t665 = t517 * t523;
t689 = t523 * mrSges(5,3);
t764 = -t517 / 0.2e1;
t811 = -mrSges(6,1) / 0.2e1;
t814 = m(7) / 0.2e1;
t816 = m(6) / 0.2e1;
t818 = m(5) / 0.2e1;
t876 = -(t667 + (-0.2e1 * qJ(4) + t749) * t523) * t818 - (-t460 * t615 + t614 + (t459 * t523 - t162) * t517) * t816 - (t309 * t400 + t310 * t403 - t550 * t59 - t551 * t58) * t814 - t350 * t764 + t689 - t613 / 0.2e1 - t665 * t811 - mrSges(6,2) * t598 + t884;
t857 = Ifges(7,6) * t253;
t872 = Ifges(7,5) * t256;
t657 = -t872 + t857;
t659 = t856 / 0.2e1 + t870 / 0.2e1;
t864 = t869 / 0.2e1 + t846 / 0.2e1;
t558 = -t872 / 0.2e1 + t857 / 0.2e1;
t875 = t256 / 0.2e1;
t787 = -t603 / 0.2e1;
t873 = m(6) * t435;
t867 = t184 + t162;
t866 = (-t459 * t548 - t460 * t547) * t816 + (-t309 * t562 + t310 * t603) * t814 + t864;
t865 = (t401 * t459 - t402 * t460) * t816 + (t253 * t309 - t256 * t310) * t814 - t659;
t609 = mrSges(7,1) * t551 - mrSges(7,2) * t550;
t863 = qJD(6) * t609;
t549 = t517 * t336 + t685 * t839;
t593 = t522 * mrSges(5,1) + t519 * mrSges(5,3);
t766 = t593 / 0.2e1;
t694 = t547 * mrSges(6,3);
t808 = mrSges(7,3) / 0.2e1;
t853 = -t124 / 0.2e1;
t852 = t163 / 0.2e1;
t851 = t562 / 0.2e1;
t850 = Ifges(5,4) + Ifges(4,5);
t849 = Ifges(4,6) - Ifges(5,6);
t845 = t562 * mrSges(7,3);
t137 = mrSges(7,1) * t256 - mrSges(7,2) * t253;
t265 = mrSges(6,1) * t402 - mrSges(6,2) * t401;
t842 = t137 + t265;
t169 = -mrSges(7,1) * t603 + mrSges(7,2) * t562;
t315 = mrSges(6,1) * t547 + mrSges(6,2) * t548;
t841 = t169 + t315;
t353 = t516 + t651;
t840 = t353 - t651;
t451 = mrSges(4,2) * t523 - mrSges(4,3) * t664;
t648 = mrSges(5,2) * t664;
t458 = -t648 - t689;
t838 = t458 + t451;
t511 = Ifges(5,5) * t519;
t837 = Ifges(5,1) * t522 + t511;
t512 = Ifges(4,4) * t522;
t836 = -Ifges(4,2) * t519 + t512;
t470 = Ifges(4,1) * t519 + t512;
t835 = Ifges(6,5) * t547 + Ifges(6,6) * t548 - t653;
t834 = Ifges(6,5) * t402 - Ifges(6,6) * t401 - t657;
t833 = -(t519 ^ 2 + t522 ^ 2) * t747 / 0.2e1;
t430 = Ifges(6,4) * t548;
t316 = -Ifges(6,2) * t547 + t430;
t832 = Ifges(6,1) * t547 + t316 + t430;
t380 = Ifges(6,4) * t401;
t246 = -t402 * Ifges(6,2) + t523 * Ifges(6,6) - t380;
t831 = Ifges(6,1) * t402 + t246 - t380;
t825 = t522 * Ifges(5,3) - t511;
t830 = t837 - t825;
t745 = pkin(8) * t523;
t753 = pkin(2) * t520;
t475 = -t745 + t753;
t632 = -pkin(7) * t519 - pkin(3);
t294 = (-qJ(5) * t523 - t475) * t522 + (-pkin(4) + t632) * t520;
t382 = -pkin(7) * t662 + t519 * t475;
t356 = t520 * qJ(4) + t382;
t308 = qJ(5) * t663 + t356;
t166 = t685 * t294 - t308 * t517;
t405 = t547 * t523;
t122 = -pkin(5) * t520 - pkin(9) * t405 + t166;
t167 = t517 * t294 + t685 * t308;
t404 = t548 * t523;
t129 = pkin(9) * t404 + t167;
t63 = t122 * t521 - t129 * t518;
t64 = t122 * t518 + t129 * t521;
t828 = t64 * mrSges(7,2) / 0.2e1 - t63 * mrSges(7,1) / 0.2e1;
t666 = t475 * t522;
t381 = pkin(7) * t664 + t666;
t827 = -t381 * t519 + t382 * t522;
t357 = t520 * t632 - t666;
t826 = t356 * t522 + t357 * t519;
t732 = Ifges(5,5) * t522;
t469 = Ifges(5,1) * t519 - t732;
t824 = t836 + t470 + t469;
t453 = -mrSges(4,1) * t523 - mrSges(4,3) * t662;
t454 = mrSges(5,1) * t523 + mrSges(5,2) * t662;
t823 = t454 / 0.2e1 - t453 / 0.2e1;
t394 = -t523 * Ifges(5,4) + t520 * t837;
t735 = Ifges(4,4) * t519;
t587 = Ifges(4,1) * t522 - t735;
t552 = t587 * t520;
t396 = -t523 * Ifges(4,5) + t552;
t643 = Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1;
t822 = -t643 * t523 + t394 / 0.2e1 + t396 / 0.2e1;
t768 = t551 / 0.2e1;
t77 = t151 * t521 - t152 * t518;
t780 = -t402 / 0.2e1;
t783 = t379 / 0.2e1;
t821 = m(5) * t783 + (t183 * t685 + t517 * t184) * t816 + (-t550 * t77 + t551 * t78) * t814 + (t401 * t764 + t685 * t780) * mrSges(6,3) + (-t253 * t768 + t550 * t875) * mrSges(7,3);
t820 = 2 * qJD(3);
t819 = -m(5) / 0.2e1;
t817 = -m(6) / 0.2e1;
t815 = -m(7) / 0.2e1;
t812 = mrSges(4,1) / 0.2e1;
t810 = -mrSges(4,2) / 0.2e1;
t807 = -Ifges(5,5) / 0.2e1;
t806 = pkin(7) * mrSges(4,1);
t805 = pkin(7) * mrSges(4,2);
t803 = -t142 / 0.2e1;
t802 = t142 / 0.2e1;
t801 = -t862 / 0.2e1;
t799 = t209 / 0.2e1;
t796 = -t256 / 0.2e1;
t255 = t404 * t521 - t405 * t518;
t795 = t255 / 0.2e1;
t793 = -t253 / 0.2e1;
t259 = t404 * t518 + t405 * t521;
t792 = t259 / 0.2e1;
t791 = t603 / 0.2e1;
t789 = -t562 / 0.2e1;
t786 = t309 / 0.2e1;
t785 = -t310 / 0.2e1;
t351 = mrSges(6,2) * t520 + mrSges(6,3) * t404;
t784 = t351 / 0.2e1;
t782 = -t401 / 0.2e1;
t779 = -t402 / 0.4e1;
t778 = t402 / 0.2e1;
t777 = t404 / 0.2e1;
t776 = t405 / 0.2e1;
t775 = t547 / 0.2e1;
t774 = -t547 / 0.2e1;
t773 = -t547 / 0.4e1;
t772 = t548 / 0.2e1;
t771 = -t548 / 0.2e1;
t767 = -t461 / 0.2e1;
t468 = t522 * Ifges(4,2) + t735;
t765 = -t468 / 0.2e1;
t763 = -t519 / 0.2e1;
t762 = t519 / 0.2e1;
t759 = t520 / 0.2e1;
t758 = -t522 / 0.2e1;
t757 = t522 / 0.2e1;
t755 = -t523 / 0.4e1;
t754 = t523 / 0.2e1;
t751 = pkin(3) * t519;
t750 = pkin(7) * t520;
t746 = pkin(8) * t522;
t741 = t58 * mrSges(7,2);
t740 = t59 * mrSges(7,1);
t739 = t77 * mrSges(7,1);
t738 = t78 * mrSges(7,2);
t734 = Ifges(6,4) * t402;
t733 = Ifges(6,4) * t547;
t719 = t255 * mrSges(7,1);
t714 = t259 * mrSges(7,2);
t125 = Ifges(7,4) * t259 + Ifges(7,2) * t255 - Ifges(7,6) * t520;
t127 = Ifges(7,1) * t259 + Ifges(7,4) * t255 - t520 * Ifges(7,5);
t138 = t714 - t719;
t210 = mrSges(7,2) * t520 + mrSges(7,3) * t255;
t212 = -mrSges(7,1) * t520 - mrSges(7,3) * t259;
t489 = qJ(4) * t661;
t346 = t523 * t599 + t489;
t245 = -pkin(5) * t404 + t346;
t247 = Ifges(6,4) * t405 + Ifges(6,2) * t404 - Ifges(6,6) * t520;
t248 = -t401 * Ifges(6,1) + t523 * Ifges(6,5) - t734;
t249 = Ifges(6,1) * t405 + Ifges(6,4) * t404 - Ifges(6,5) * t520;
t695 = t405 * mrSges(6,2);
t696 = t404 * mrSges(6,1);
t266 = t695 - t696;
t352 = -mrSges(6,1) * t520 - mrSges(6,3) * t405;
t578 = Ifges(5,3) * t519 + t732;
t391 = Ifges(5,6) * t520 + t523 * t578;
t393 = Ifges(4,6) * t520 + t523 * t836;
t395 = Ifges(5,4) * t520 + t523 * t837;
t397 = Ifges(4,5) * t520 + t523 * t587;
t631 = pkin(7) + t751;
t398 = t520 * t631 - t488;
t399 = t523 * t631 - t489;
t592 = mrSges(5,1) * t519 - mrSges(5,3) * t522;
t432 = t592 * t520;
t433 = t592 * t523;
t594 = mrSges(4,1) * t519 + mrSges(4,2) * t522;
t434 = t594 * t523;
t452 = -mrSges(4,2) * t520 - mrSges(4,3) * t663;
t455 = mrSges(4,1) * t520 - mrSges(4,3) * t661;
t490 = mrSges(5,2) * t661;
t691 = t520 * mrSges(5,1);
t456 = t490 - t691;
t457 = -mrSges(5,2) * t663 + mrSges(5,3) * t520;
t557 = Ifges(7,5) * t792 + Ifges(7,6) * t795;
t542 = Ifges(6,5) * t776 + Ifges(6,6) * t777 + t557;
t492 = Ifges(5,5) * t662;
t390 = -Ifges(5,6) * t523 + Ifges(5,3) * t664 + t492;
t688 = t523 * Ifges(4,6);
t392 = t520 * t836 - t688;
t619 = t390 / 0.2e1 - t392 / 0.2e1;
t642 = Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t3 = t126 * t792 + t127 * t793 + t124 * t795 + t125 * t796 + t248 * t776 + t246 * t777 + t247 * t780 + t249 * t782 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t523 + t822 * t522 + (-t642 * t523 + t619) * t519 + t542) * t523 + (-pkin(1) * mrSges(3,1) + Ifges(6,5) * t401 / 0.2e1 + Ifges(6,6) * t778 + Ifges(7,5) * t253 / 0.2e1 + Ifges(7,6) * t875 + pkin(7) * t434 - Ifges(3,4) * t520 + (t395 / 0.2e1 + t397 / 0.2e1 + t643 * t520) * t522 + (t391 / 0.2e1 - t393 / 0.2e1 + t642 * t520) * t519 + (-Ifges(5,2) - Ifges(4,3) - Ifges(3,2) + Ifges(3,1) - Ifges(6,3) - Ifges(7,3) + (m(4) * pkin(7) + t594) * pkin(7)) * t523) * t520 + m(4) * (t379 * t382 - t381 * t651) - t651 * t455 + t353 * t456 + t348 * t457 + t356 * t458 + t382 * t451 + t379 * t452 + t381 * t453 + t357 * t454 + t399 * t432 + t398 * t433 + t163 * t351 + t162 * t352 + t345 * t266 + t346 * t265 + t167 * t349 + t166 * t350 + t244 * t138 + t245 * t137 + t64 * t209 + t59 * t210 + t63 * t211 + t58 * t212 + m(7) * (t244 * t245 + t58 * t63 + t59 * t64) + m(6) * (t162 * t166 + t163 * t167 + t345 * t346) + m(5) * (t348 * t356 + t353 * t357 + t398 * t399);
t712 = t3 * qJD(1);
t701 = t398 * mrSges(5,1);
t700 = t398 * mrSges(5,3);
t377 = (t522 * t524 - t508) * t520;
t264 = pkin(5) * t401 + t377;
t431 = t573 * t520;
t491 = Ifges(5,6) * t662;
t582 = -Ifges(6,2) * t401 + t734;
t591 = t401 * mrSges(6,1) + t402 * mrSges(6,2);
t4 = t831 * t782 + t834 * t754 + t585 * t793 + t580 * t796 + t248 * t778 + t582 * t780 + (-t162 * t402 - t163 * t401) * mrSges(6,3) + (-t253 * t59 - t256 * t58) * mrSges(7,3) + t253 * t853 + t126 * t875 + m(5) * (-t348 * t651 + t353 * t379 + t398 * t431) - t838 * t651 + ((t688 / 0.2e1 + t701 - t348 * mrSges(5,2) - t379 * mrSges(4,3) + t492 / 0.2e1 + (t806 - t512 / 0.2e1) * t520 + t619) * t522 + (t700 - t353 * mrSges(5,2) - t651 * mrSges(4,3) + (-t805 + (Ifges(4,4) / 0.2e1 + t807) * t519) * t520 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(5,1) / 0.2e1) * t662 - t822) * t519) * t520 + t431 * t432 + t377 * t265 + t184 * t349 + t183 * t350 + t264 * t137 + t78 * t209 + t77 * t211 + (-t453 + t454) * t379 + t345 * t591 + m(7) * (t244 * t264 + t58 * t77 + t59 * t78) + m(6) * (t162 * t183 + t163 * t184 + t345 * t377) - t523 * t491 / 0.2e1 + t894;
t699 = t4 * qJD(1);
t693 = t548 * mrSges(6,3);
t687 = t59 * t562;
t7 = t657 * t754 + t58 * t209 - t59 * t211 - t894 - (-t59 * mrSges(7,3) - t585 / 0.2e1 + t853) * t253 - (-t58 * mrSges(7,3) + t126 / 0.2e1 - t580 / 0.2e1) * t256;
t686 = t7 * qJD(1);
t25 = -t256 * t209 + t253 * t211 - t402 * t349 + t401 * t350 + m(7) * (t253 * t58 - t256 * t59) + m(6) * (t162 * t401 - t163 * t402);
t684 = qJD(1) * t25;
t65 = t310 * mrSges(7,1) + mrSges(7,2) * t309;
t683 = qJD(6) * t65;
t682 = t862 * t253;
t22 = t403 * t209 + t400 * t211 + m(7) * (t400 * t58 + t403 * t59) + t350 * t665 + (-t613 + m(6) * (t517 * t162 - t614) - t458 - m(5) * t348) * t523 + (-m(5) * t398 + m(6) * t345 + m(7) * t244 - t432 + t842) * t662;
t681 = t22 * qJD(1);
t680 = t336 * t402;
t679 = t839 * t401;
t671 = t550 * t256;
t669 = t551 * t253;
t668 = t461 * t520;
t652 = -t348 + t379;
t650 = t814 + t816;
t640 = mrSges(7,3) * t803;
t638 = -t845 / 0.2e1;
t635 = -t693 / 0.2e1;
t633 = t685 / 0.2e1;
t630 = t402 * t775;
t617 = -mrSges(5,2) * qJ(4) - Ifges(4,6);
t600 = mrSges(5,2) * pkin(3) - t850;
t438 = t510 + t641;
t595 = t522 * mrSges(4,1) - t519 * mrSges(4,2);
t590 = -mrSges(6,1) * t548 + mrSges(6,2) * t547;
t581 = Ifges(6,2) * t548 + t733;
t317 = Ifges(6,1) * t548 - t733;
t328 = -pkin(5) * t548 + t438;
t465 = -t510 + t751;
t525 = t832 * t401 / 0.4e1 + t835 * t755 + t209 * t802 + t248 * t773 + t317 * t779 + t431 * t766 + t890 + t891 * t256 + (-t350 / 0.2e1 + t867 * t817) * t839 + t883 * t253 - t465 * t432 / 0.2e1 - t438 * t265 / 0.2e1 - t377 * t315 / 0.2e1 - t328 * t137 / 0.2e1 - t264 * t169 / 0.2e1 + (t244 * t328 + t264 * t326 + t886 * t862 + (-t59 + t77) * t142) * t815 - t831 * t548 / 0.4e1 + (t345 * t438 + t377 * t435) * t817 + ((t163 - t183) * t817 - t349 / 0.2e1) * t336 - t345 * t590 / 0.2e1 - t435 * t591 / 0.2e1;
t527 = (t309 * t63 + t310 * t64) * t814 + (-pkin(3) * t357 + qJ(4) * t356) * t818 + (t166 * t459 + t167 * t460) * t816 - pkin(3) * t456 / 0.2e1 + qJ(4) * t457 / 0.2e1 + t459 * t352 / 0.2e1 + t460 * t784 + t381 * t812 + t382 * t810 + t356 * mrSges(5,3) / 0.2e1 - t357 * mrSges(5,1) / 0.2e1 + t212 * t786 + t310 * t210 / 0.2e1 + t167 * mrSges(6,2) / 0.2e1 + t166 * t811 - t542 + t828;
t559 = Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t560 = Ifges(4,1) / 0.4e1 + Ifges(5,1) / 0.4e1 - Ifges(4,2) / 0.4e1 - Ifges(5,3) / 0.4e1;
t561 = (mrSges(4,3) + mrSges(5,2)) * pkin(8) / 0.2e1;
t563 = t398 * t465 + t431 * t461;
t2 = t527 + t525 + (t746 * t840 + t563) * t819 + ((pkin(2) * t812 + mrSges(5,1) * t767 - t805 / 0.2e1 - t511 / 0.4e1 + t468 / 0.4e1 + t825 / 0.4e1 + (-t560 + t561) * t522) * t522 + (pkin(2) * t810 + mrSges(5,3) * t767 + t470 / 0.4e1 + t469 / 0.4e1 + t512 / 0.4e1 - t806 / 0.2e1 + (t560 + t561) * t519 + (0.3e1 / 0.4e1 * Ifges(4,4) + t807) * t522) * t519 + t559) * t520 + (t401 * t773 - t548 * t779) * Ifges(6,2) + (t682 / 0.2e1 + t142 * t875 - t687 / 0.2e1 + t77 * t851 + t886 * t787) * mrSges(7,3) + Ifges(6,4) * t630 + (t679 / 0.2e1 - t680 / 0.2e1 - (-t183 / 0.2e1 + t852) * t548 - (-t184 / 0.2e1 - t162 / 0.2e1) * t547) * mrSges(6,3) + (-t492 / 0.4e1 - t701 / 0.2e1 + t392 / 0.4e1 - t390 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(5,6) - 0.3e1 / 0.4e1 * Ifges(4,6)) * t523 + (-t379 / 0.2e1 + t348 / 0.2e1) * mrSges(5,2) + (t458 / 0.2e1 + t451 / 0.2e1 + t652 * t819) * pkin(8)) * t519 + (t700 / 0.2e1 - t396 / 0.4e1 - t394 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(5,4) + 0.3e1 / 0.4e1 * Ifges(4,5)) * t523 - t823 * pkin(8) + (t651 / 0.2e1 - t353 / 0.2e1) * mrSges(5,2)) * t522;
t8 = -t170 * t791 + t581 * t774 + t317 * t775 + t173 * t787 + t468 * t763 + t578 * t758 - t465 * t593 + t435 * t590 + (m(5) * t465 + t592) * t461 - pkin(2) * t594 + t893 + t832 * t772 + (t587 + t830) * t762 + t824 * t757 + (t171 - t172) * t851 + (t315 + t873) * t438 + (t169 + t887) * t328;
t572 = -t2 * qJD(1) + t8 * qJD(2);
t21 = -t893 + (-t171 / 0.2e1 + t172 / 0.2e1) * t562 + (t170 / 0.2e1 + t173 / 0.2e1) * t603;
t530 = -(t640 - t891) * t256 - (mrSges(7,3) * t801 - t883) * t253 + t142 * t799 + t523 * t653 / 0.4e1 + t890;
t541 = Ifges(7,3) * t759 - t557 + t828;
t6 = t530 + t541;
t571 = t6 * qJD(1) + t21 * qJD(2);
t570 = t650 * t662;
t526 = (-t398 * t519 + (-t668 - t745) * t522) * t818 + (t519 * t345 - t523 * t549) * t816 + (t142 * t400 + t244 * t519 + t403 * t862) * t814 + t432 * t763 + t400 * t638 + t403 * t880 + t635 * t665 - t598 * t694 + t842 * t762 + (t435 * t816 + t326 * t814 + t766 + t841 / 0.2e1) * t662;
t531 = t357 * t818 + (t166 * t685 + t517 * t167) * t816 + (-t550 * t63 + t551 * t64) * t814 + t212 * t769 + t210 * t768 + t517 * t784 - t691 / 0.2e1 + t352 * t633;
t11 = t526 - t490 - t531;
t66 = (-m(5) * t461 + t593 + t841 + t873 + t887) * t519;
t569 = qJD(1) * t11 + qJD(2) * t66;
t529 = (t253 * t789 - t256 * t791) * mrSges(7,3) + (t401 * t771 + t630) * mrSges(6,3) + (-t162 * t548 - t163 * t547 - t336 * t401 - t402 * t839) * t816 + (t142 * t253 - t256 * t862 - t562 * t58 + t59 * t603) * t814 + t562 * t798 + t603 * t799 + t349 * t774 + t350 * t771;
t533 = t346 * t816 + t245 * t814 - t719 / 0.2e1 + t714 / 0.2e1 - t696 / 0.2e1 + t695 / 0.2e1;
t16 = t529 - t533;
t35 = (t562 ^ 2 + t603 ^ 2) * mrSges(7,3) + (t547 ^ 2 + t548 ^ 2) * mrSges(6,3) + m(7) * (-t142 * t562 + t603 * t862) + m(6) * (t336 * t548 - t547 * t839);
t568 = -qJD(1) * t16 - qJD(2) * t35;
t544 = t264 * t814 + t377 * t816 + t659;
t31 = t544 + t591 - t865;
t543 = t328 * t814 + t438 * t816 - t864;
t40 = t543 + t590 - t866;
t567 = qJD(1) * t31 + qJD(2) * t40;
t67 = -0.2e1 * t659;
t84 = 0.2e1 * t864;
t566 = qJD(1) * t67 + qJD(2) * t84;
t537 = (t401 * t685 - t517 * t402) * t816 + (-t253 * t550 - t256 * t551) * t814;
t70 = t570 - t537;
t536 = (-t517 * t547 - t548 * t685) * t816 + (t550 * t562 + t551 * t603) * t814;
t83 = -t519 * t650 + t536;
t565 = qJD(1) * t70 - qJD(2) * t83;
t19 = (-t669 / 0.2e1 + t671 / 0.2e1) * mrSges(7,3) + t884;
t564 = t19 * qJD(1) - qJD(3) * t609;
t534 = (-t253 * t785 + t309 * t875) * mrSges(7,3) + t209 * t786 + t211 * t785 - t558;
t10 = (t58 / 0.2e1 + t78 / 0.2e1) * mrSges(7,2) + (t59 / 0.2e1 - t77 / 0.2e1) * mrSges(7,1) + t534 + t558;
t23 = (t802 + t803) * mrSges(7,2) + (t862 / 0.2e1 + t801) * mrSges(7,1);
t546 = t10 * qJD(1) - t23 * qJD(2) + t65 * qJD(3);
t14 = t821 + t876;
t535 = t549 * t816 + t814 * t882;
t538 = t549 * t817 + t815 * t882;
t29 = (-(t789 + t851) * t551 - (t791 + t787) * t550) * mrSges(7,3) + t535 + t538;
t76 = m(5) * qJ(4) + t517 * mrSges(6,1) + t685 * mrSges(6,2) + mrSges(5,3) + m(6) * (-t459 * t517 + t460 * t685) + m(7) * (-t309 * t551 - t310 * t550) + t609;
t545 = -t14 * qJD(1) + t29 * qJD(2) + t76 * qJD(3);
t82 = t536 + (m(6) + m(7)) * t762;
t71 = t570 + t537;
t41 = t543 + t866;
t32 = t544 + t865;
t28 = -t550 * t880 + t845 * t768 + (m(5) * pkin(8) + mrSges(5,2)) * t522 + (t562 * t768 + t603 * t769) * mrSges(7,3) + t535 - t538 + (t517 * t548 - 0.2e1 * t547 * t633) * mrSges(6,3);
t20 = t553 + t660 + (t669 - t671) * t808;
t15 = t529 + t533;
t13 = t526 + t531;
t12 = -t648 + t821 - t876;
t9 = t741 / 0.2e1 + t740 / 0.2e1 - t738 / 0.2e1 + t739 / 0.2e1 + t534 - t558;
t5 = t530 - t541;
t1 = (-Ifges(5,1) * t664 + t390 + t492) * t519 / 0.4e1 + (t552 + t396 + t394) * t522 / 0.4e1 + (t578 - t470) * t664 / 0.4e1 - t838 * t748 / 0.2e1 + ((t652 * t519 + t522 * t840) * pkin(8) + t563) * t818 + t833 * mrSges(4,3) + (-t519 * t849 + t522 * t850) * t755 + t582 * t773 + t581 * t779 + t668 * t766 - t824 * t664 / 0.4e1 + t823 * t746 + t886 * t880 + t527 - t525 - t867 * t694 / 0.2e1 - (t679 - t680) * mrSges(6,3) / 0.2e1 - t595 * t753 / 0.2e1 + t594 * t750 / 0.2e1 + (t830 / 0.4e1 - t825 / 0.4e1 + t765) * t662 + (t833 + (t783 - t348 / 0.2e1) * t519 + t840 * t757) * mrSges(5,2) + (t519 * t642 + t522 * t643) * t523 - t519 * t392 / 0.4e1 + t398 * t592 / 0.2e1 + t559 * t520 + t693 * t852 + (t687 - t682) * t808 + t183 * t635 + t77 * t638 + t256 * t640;
t17 = [qJD(2) * t3 + qJD(3) * t4 + qJD(4) * t22 + qJD(5) * t25 + qJD(6) * t7, t1 * qJD(3) + t13 * qJD(4) + t15 * qJD(5) + t5 * qJD(6) + t712 + (m(4) * (-pkin(7) * t752 + pkin(8) * t827) + t827 * mrSges(4,3) + 0.2e1 * (pkin(8) * t826 + t399 * t461) * t818 + t826 * mrSges(5,2) + (t519 * t850 + t522 * t849) * t759 + t125 * t791 + t173 * t792 + t171 * t795 + t247 * t774 + t317 * t776 + t316 * t777 + t391 * t758 + t249 * t772 + mrSges(3,2) * t750 + t393 * t757 - t166 * t693 - t167 * t694 - t399 * t593 - t63 * t845 + t127 * t851 - t336 * t352 + 0.2e1 * (-t166 * t336 + t167 * t839 + t346 * t435) * t816 + (Ifges(3,5) + (t469 / 0.2e1 + t470 / 0.2e1) * t522 + (-t825 / 0.2e1 + t765) * t519 + (-mrSges(3,1) - t595) * pkin(7)) * t523 + 0.2e1 * (t142 * t63 + t245 * t326 + t64 * t862) * t814 + t862 * t210 + t64 * t868 + t839 * t351 - Ifges(3,6) * t520 + t461 * t433 - pkin(2) * t434 + t435 * t266 + t346 * t315 + t326 * t138 + t245 * t169 + t142 * t212 + ((t452 + t457) * t522 + (-t455 + t456) * t519) * pkin(8) - (Ifges(6,5) * t548 + Ifges(7,5) * t562 - Ifges(6,6) * t547 + Ifges(7,6) * t603) * t520 / 0.2e1 + (t397 + t395) * t762) * qJD(2), t699 + t1 * qJD(2) + t12 * qJD(4) + t32 * qJD(5) + t9 * qJD(6) + ((t309 * t77 + t310 * t78) * t814 + (-pkin(3) * t379 - qJ(4) * t651) * t818 + (t183 * t459 + t184 * t460) * t816) * t820 + (-t183 * mrSges(6,1) + t184 * mrSges(6,2) - t459 * t697 - t460 * t698 + t491 + t738 - t739 + (t519 * t600 + t522 * t617) * t520 + (-mrSges(4,1) - mrSges(5,1)) * t379 - (-mrSges(4,2) + mrSges(5,3)) * t651 + (-t253 * t310 - t256 * t309) * mrSges(7,3) - t834) * qJD(3), t681 + t13 * qJD(2) + t12 * qJD(3) + m(7) * (-t400 * t550 + t403 * t551) * qJD(4) + t71 * qJD(5) + t20 * qJD(6), qJD(2) * t15 + qJD(3) * t32 + qJD(4) * t71 + t684, t686 + t5 * qJD(2) + t9 * qJD(3) + t20 * qJD(4) + (t657 - t740 - t741) * qJD(6); -qJD(3) * t2 + qJD(4) * t11 + qJD(5) * t16 + qJD(6) * t6 - t712, qJD(3) * t8 + qJD(4) * t66 + qJD(5) * t35 + qJD(6) * t21, t28 * qJD(4) + t41 * qJD(5) - t896 + ((-t142 * t310 + t309 * t862) * t814 + (t336 * t460 + t459 * t839) * t816) * t820 + t572 + (-t839 * mrSges(6,1) + t336 * mrSges(6,2) + t309 * t868 + t310 * t845 - t459 * t694 + t460 * t693 - t600 * t522 + (Ifges(5,6) + t617) * t519 + (-m(5) * t573 - t593 - t595) * pkin(8) - t835 - t892) * qJD(3), qJD(3) * t28 + qJD(5) * t82 + t569, qJD(3) * t41 + qJD(4) * t82 - t568, -qJD(3) * t895 + t571 + t896; qJD(2) * t2 - qJD(4) * t14 - qJD(5) * t31 + qJD(6) * t10 - t699, qJD(4) * t29 - qJD(5) * t40 - qJD(6) * t23 - t572, qJD(4) * t76 + t683, t545, -t567, t546 - t683; -qJD(2) * t11 + qJD(3) * t14 - qJD(5) * t70 - qJD(6) * t19 - t681, -qJD(3) * t29 + qJD(5) * t83 - t569, -t545 + t863, 0, -t565, -t564 - t863; -qJD(2) * t16 + qJD(3) * t31 + qJD(4) * t70 + qJD(6) * t67 - t684, qJD(3) * t40 - qJD(4) * t83 + qJD(6) * t84 + t568, t567, t565, 0, t566; -qJD(2) * t6 - qJD(3) * t10 + qJD(4) * t19 - qJD(5) * t67 - t686, qJD(3) * t23 - qJD(5) * t84 - t571, -qJD(4) * t609 - t546, t564, -t566, 0;];
Cq  = t17;
