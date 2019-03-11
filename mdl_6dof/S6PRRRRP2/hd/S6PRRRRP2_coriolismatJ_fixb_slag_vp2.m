% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:30
% EndTime: 2019-03-09 00:01:53
% DurationCPUTime: 15.03s
% Computational Cost: add. (20000->726), mult. (45519->960), div. (0->0), fcn. (48774->10), ass. (0->411)
t477 = sin(qJ(5));
t479 = cos(qJ(5));
t428 = -t479 * mrSges(7,1) - t477 * mrSges(7,3);
t711 = qJ(6) * t479;
t744 = pkin(5) * t477;
t430 = -t711 + t744;
t732 = Ifges(7,5) * t479;
t434 = Ifges(7,3) * t477 + t732;
t469 = Ifges(6,4) * t479;
t552 = Ifges(6,2) * t477 - t469;
t759 = -t479 / 0.2e1;
t761 = t477 / 0.2e1;
t466 = Ifges(7,5) * t477;
t433 = -Ifges(7,3) * t479 + t466;
t733 = Ifges(6,4) * t477;
t437 = Ifges(6,2) * t479 + t733;
t762 = -t477 / 0.2e1;
t439 = Ifges(7,1) * t477 - t732;
t441 = Ifges(6,1) * t477 + t469;
t885 = t439 + t441;
t830 = t433 * t762 + t437 * t761 + t759 * t885;
t440 = Ifges(7,1) * t479 + t466;
t442 = Ifges(6,1) * t479 - t733;
t881 = t440 + t442;
t514 = t430 * t428 + t881 * t761 + (t434 + t552) * t759 - t830;
t472 = t477 ^ 2;
t474 = t479 ^ 2;
t635 = t472 + t474;
t478 = sin(qJ(3));
t480 = cos(qJ(3));
t476 = sin(pkin(6));
t755 = sin(qJ(2));
t607 = t476 * t755;
t713 = cos(pkin(6));
t399 = t478 * t713 + t480 * t607;
t516 = t478 * t607 - t480 * t713;
t754 = sin(qJ(4));
t756 = cos(qJ(4));
t288 = t399 * t754 + t756 * t516;
t601 = t754 * t478;
t419 = -t480 * t756 + t601;
t462 = -pkin(3) * t480 - pkin(2);
t602 = t756 * t478;
t420 = -t480 * t754 - t602;
t741 = pkin(10) * t420;
t324 = pkin(4) * t419 + t462 + t741;
t775 = -pkin(9) - pkin(8);
t447 = t775 * t480;
t810 = -t756 * t447 + t775 * t601;
t149 = t324 * t479 - t477 * t810;
t646 = t479 * t810;
t654 = t477 * t324;
t150 = t646 + t654;
t528 = t149 * t477 - t150 * t479 + t810;
t886 = t288 * t528;
t656 = t477 * qJ(6);
t551 = -t479 * pkin(5) - t656;
t426 = -pkin(4) + t551;
t634 = t756 * pkin(3);
t404 = -t634 + t426;
t637 = t404 + t426;
t715 = t479 * mrSges(6,2);
t718 = t477 * mrSges(6,1);
t432 = t715 + t718;
t461 = -t634 - pkin(4);
t658 = t461 * t432;
t747 = pkin(4) * t432;
t714 = t479 * mrSges(7,3);
t717 = t477 * mrSges(7,1);
t431 = -t714 + t717;
t751 = m(7) * t430;
t626 = t751 / 0.2e1;
t834 = t626 + t431 / 0.2e1;
t884 = t747 / 0.2e1 - t658 / 0.2e1 - t834 * t637 - t514;
t791 = 0.2e1 * m(7);
t481 = cos(qJ(2));
t657 = t476 * t481;
t883 = -t657 / 0.2e1;
t500 = t399 * t756 - t516 * t754;
t319 = t432 * t420;
t318 = t431 * t420;
t862 = t318 / 0.2e1;
t832 = t319 / 0.2e1 + t862;
t882 = t832 * t500;
t668 = t419 * t477;
t328 = mrSges(6,2) * t420 + mrSges(6,3) * t668;
t335 = mrSges(7,2) * t668 - mrSges(7,3) * t420;
t813 = t328 + t335;
t880 = t635 * t756;
t869 = t635 * pkin(10) * t288;
t879 = t426 * t500 - t869;
t878 = -pkin(4) * t500 - t869;
t195 = t477 * t500 + t479 * t657;
t196 = -t477 * t657 + t479 * t500;
t825 = m(6) + m(7);
t847 = t288 * t479;
t848 = t288 * t477;
t864 = t825 * (-t195 * t848 - t196 * t847 + t288 * t500);
t877 = t864 * qJD(1);
t665 = t420 * t479;
t332 = mrSges(6,1) * t419 + mrSges(6,3) * t665;
t333 = -mrSges(7,1) * t419 - mrSges(7,2) * t665;
t666 = t420 * t477;
t329 = -mrSges(6,2) * t419 + mrSges(6,3) * t666;
t722 = t419 * mrSges(7,3);
t334 = mrSges(7,2) * t666 + t722;
t812 = t329 + t334;
t863 = t288 / 0.2e1;
t872 = -t848 / 0.2e1;
t876 = t333 * t872 + t477 * t332 * t863 - t847 * t812 / 0.2e1;
t875 = t196 / 0.2e1;
t874 = -t288 / 0.2e1;
t667 = t419 * t479;
t628 = mrSges(7,2) * t667;
t720 = t420 * mrSges(7,1);
t331 = -t628 + t720;
t873 = t331 / 0.2e1;
t821 = mrSges(7,2) + mrSges(6,3);
t670 = t419 * qJ(6);
t104 = t150 + t670;
t105 = -pkin(5) * t419 - t149;
t811 = -pkin(5) * t668 + qJ(6) * t667 + t810;
t529 = -t104 * t479 - t105 * t477 + t811;
t871 = t288 * t529;
t870 = t635 * t821;
t374 = t419 * t657;
t307 = -t479 * t374 + t477 * t607;
t679 = t307 * t479;
t306 = -t374 * t477 - t479 * t607;
t680 = t306 * t477;
t826 = mrSges(5,2) / 0.2e1;
t868 = t374 * t826 + t821 * (t679 / 0.2e1 + t680 / 0.2e1);
t207 = -Ifges(7,6) * t420 - t419 * t434;
t209 = -Ifges(6,6) * t420 + t419 * t552;
t211 = -Ifges(7,4) * t420 - t419 * t440;
t213 = -Ifges(6,5) * t420 - t419 * t442;
t758 = t479 / 0.2e1;
t766 = -t420 / 0.2e1;
t854 = Ifges(7,4) + Ifges(6,5);
t506 = Ifges(5,6) * t420 + t207 * t759 + t209 * t758 + ((Ifges(6,6) - Ifges(7,6)) * t479 + t854 * t477) * t766 + (t211 + t213) * t761 + (-Ifges(5,5) + t830) * t419;
t838 = t811 * t428;
t429 = -t479 * mrSges(6,1) + t477 * mrSges(6,2);
t839 = t810 * t429;
t851 = t810 * mrSges(5,1);
t362 = -t447 * t754 - t775 * t602;
t852 = t362 * mrSges(5,2);
t867 = t506 + t838 + t839 + t852 - t851;
t865 = -t838 / 0.2e1 - t839 / 0.2e1 + t851 / 0.2e1 - t852 / 0.2e1;
t861 = -t335 / 0.2e1;
t860 = -t431 / 0.2e1;
t859 = -t432 / 0.2e1;
t857 = m(5) * t462;
t856 = pkin(4) * t810;
t853 = Ifges(6,3) + Ifges(7,2);
t809 = t428 + t429;
t850 = -mrSges(5,1) + t809;
t173 = -t420 * t430 + t362;
t849 = t173 * t811;
t845 = t362 * t477;
t820 = t362 * t500;
t844 = t362 * t754;
t843 = t404 * t811;
t842 = t426 * t811;
t841 = t461 * t810;
t840 = t479 * t362;
t676 = t810 * t362;
t640 = t105 + t149;
t435 = Ifges(6,5) * t479 - Ifges(6,6) * t477;
t436 = Ifges(7,4) * t479 + Ifges(7,6) * t477;
t837 = t436 + t435;
t836 = t478 ^ 2 + t480 ^ 2;
t835 = t880 * pkin(3);
t633 = t754 * pkin(3);
t460 = t633 + pkin(10);
t833 = t635 * t460;
t571 = t472 / 0.2e1 + t474 / 0.2e1;
t511 = t821 * t571;
t555 = t429 / 0.2e1 + t428 / 0.2e1 - mrSges(5,1) / 0.2e1;
t831 = (t826 - t511) * t288 + t555 * t500;
t316 = t431 * t419;
t317 = t432 * t419;
t330 = -mrSges(6,1) * t420 + mrSges(6,3) * t667;
t828 = t813 * t875 + (t317 + t316) * t874 + (t873 - t330 / 0.2e1) * t195;
t397 = Ifges(7,5) * t665;
t208 = Ifges(7,6) * t419 - Ifges(7,3) * t666 - t397;
t721 = t419 * Ifges(6,6);
t210 = t420 * t552 + t721;
t576 = -t208 / 0.2e1 + t210 / 0.2e1;
t620 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t621 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t212 = t419 * Ifges(7,4) - t420 * t440;
t214 = t419 * Ifges(6,5) - t420 * t442;
t797 = -t621 * t419 - t212 / 0.2e1 - t214 / 0.2e1;
t827 = -t811 * t318 - t810 * t319 + (-t462 * mrSges(5,2) + Ifges(5,4) * t419 + t797 * t479 + (-t419 * t620 + t576) * t477) * t419 + (-t462 * mrSges(5,1) - Ifges(5,4) * t420 + (-t213 / 0.2e1 - t211 / 0.2e1 + t621 * t420) * t479 + (t209 / 0.2e1 - t207 / 0.2e1 + t620 * t420) * t477 + (Ifges(5,1) - Ifges(5,2) - t853) * t419) * t420 + t104 * t335 + t105 * t331 + t149 * t330 + t150 * t328 - t173 * t316 - t362 * t317;
t789 = -m(6) / 0.2e1;
t787 = -m(7) / 0.2e1;
t785 = m(5) * pkin(3);
t824 = -t785 / 0.2e1;
t822 = mrSges(6,1) + mrSges(7,1);
t819 = t404 * t500;
t817 = t461 * t500;
t816 = t479 * t640;
t814 = -t318 - t319;
t579 = t428 * t758;
t781 = -mrSges(7,1) / 0.2e1;
t808 = (t579 + t781) * t420 + t628;
t742 = pkin(10) * t419;
t337 = -pkin(4) * t420 + t742;
t749 = pkin(3) * t478;
t325 = t337 + t749;
t160 = t477 * t325 - t840;
t405 = t420 * qJ(6);
t107 = -t405 + t160;
t159 = t325 * t479 + t845;
t745 = pkin(5) * t420;
t108 = -t159 + t745;
t542 = t107 * t479 + t108 * t477;
t803 = -mrSges(4,1) * t480 + mrSges(4,2) * t478;
t801 = pkin(5) * t873 + qJ(6) * t861;
t463 = m(7) * qJ(6) + mrSges(7,3);
t800 = t420 * t821;
t660 = t460 * t479;
t799 = t813 * t660;
t798 = (mrSges(4,1) * t478 + mrSges(4,2) * t480) * t883;
t592 = -t667 / 0.2e1;
t464 = t479 * mrSges(7,2);
t609 = t464 / 0.2e1;
t796 = mrSges(7,2) * t592 + t720 / 0.2e1 + t420 * t579 + t419 * t609;
t795 = t332 * t759 + t333 * t758 + t762 * t812;
t313 = t551 * t420;
t794 = t430 * t862 + t362 * t859 - t313 * t428 / 0.2e1 + t173 * t860;
t525 = (-mrSges(5,1) * t420 - mrSges(5,2) * t419) * t657;
t719 = t420 * mrSges(5,3);
t793 = t500 * t719 / 0.2e1 - t525 / 0.2e1 + t828;
t790 = 2 * qJD(4);
t788 = m(6) / 0.2e1;
t786 = m(7) / 0.2e1;
t784 = m(7) * pkin(3);
t782 = mrSges(6,1) / 0.2e1;
t780 = mrSges(7,1) / 0.2e1;
t779 = -mrSges(6,2) / 0.2e1;
t778 = mrSges(6,2) / 0.2e1;
t777 = -mrSges(7,3) / 0.2e1;
t776 = mrSges(7,3) / 0.2e1;
t315 = t429 * t420;
t771 = t315 / 0.2e1;
t768 = t404 / 0.2e1;
t753 = m(7) * t108;
t169 = t337 * t479 + t845;
t114 = -t169 + t745;
t752 = m(7) * t114;
t750 = m(7) * t477;
t748 = pkin(4) * t317;
t740 = pkin(10) * t479;
t729 = t500 * mrSges(5,1);
t728 = t288 * mrSges(5,2);
t723 = t419 * mrSges(5,3);
t716 = t477 * mrSges(6,3);
t170 = t477 * t337 - t840;
t113 = -t405 + t170;
t708 = t113 * t479;
t707 = t114 * t477;
t706 = t150 * t477;
t705 = t159 * t477;
t704 = t160 * t479;
t703 = t169 * t477;
t702 = t170 * t479;
t700 = t173 * t430;
t698 = t173 * t477;
t373 = t420 * t657;
t168 = t288 * t373;
t608 = t476 ^ 2 * t755;
t26 = m(5) * (-t374 * t500 - t481 * t608 - t168) + m(4) * (-t608 + (t399 * t480 + t478 * t516) * t476) * t481 + t825 * (t195 * t306 + t196 * t307 - t168);
t690 = t26 * qJD(1);
t689 = t500 * t428;
t688 = t500 * t429;
t683 = t288 * t472;
t682 = t288 * t474;
t674 = t362 * t373;
t671 = t404 * t316;
t669 = t419 * t196;
t664 = t426 * t316;
t659 = t461 * t317;
t653 = t477 * t330;
t652 = t477 * t331;
t641 = -t104 + t150;
t639 = t833 * t288;
t636 = t835 * pkin(10);
t631 = mrSges(7,2) * t707;
t630 = mrSges(6,3) * t704;
t629 = mrSges(6,3) * t702;
t627 = -t751 / 0.2e1;
t623 = t782 + t780;
t622 = t778 + t777;
t619 = t460 * t653;
t618 = t460 * t652;
t615 = t721 / 0.2e1;
t614 = -t719 / 0.2e1;
t613 = mrSges(7,2) * t762;
t612 = -t716 / 0.2e1;
t611 = t716 / 0.2e1;
t610 = -t464 / 0.2e1;
t606 = t477 * t756;
t605 = t479 * t756;
t590 = -t665 / 0.4e1;
t586 = t318 * t762;
t574 = -t332 / 0.2e1 + t333 / 0.2e1;
t573 = t334 / 0.2e1 + t329 / 0.2e1;
t572 = t860 + t859;
t570 = t431 + t751;
t564 = -t634 / 0.2e1;
t563 = t634 / 0.2e1;
t561 = t500 * t614;
t557 = (t679 + t680) * pkin(10);
t483 = m(5) * t749 * t883 + (-t159 * t195 + t160 * t196 + t820 + t886) * t788 + (t107 * t196 + t108 * t195 + t871) * t786 + t798 + (t173 * t786 + t614 - t832) * t500 + t876;
t486 = t754 * t374 * t824 + t798 + (t788 + t786) * (t307 * t660 + t460 * t680) - (t404 * t786 + t461 * t788 + t756 * t824 + t555) * t373 + t868;
t517 = t525 / 0.2e1;
t2 = -t483 + t486 + t517 + t561 - t828;
t336 = mrSges(5,1) * t419 - mrSges(5,2) * t420;
t4 = (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t478 + pkin(3) * t336) * t478 + t749 * t857 + m(6) * (t149 * t159 + t150 * t160 + t676) + m(7) * (t104 * t107 + t105 * t108 + t849) + t107 * t334 + t160 * t329 + t159 * t332 + t108 * t333 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) * t480 + (Ifges(4,1) - Ifges(4,2)) * t478) * t480 + t827;
t550 = -t2 * qJD(1) + t4 * qJD(2);
t489 = -t555 * t373 + (pkin(4) * t373 + t557) * t788 + (-t373 * t426 + t557) * t786 + t868;
t521 = -t169 * t195 + t170 * t196 + t820;
t522 = t113 * t196 + t114 * t195 + t173 * t500;
t6 = t517 + t882 + (t861 - t328 / 0.2e1) * t196 + (t330 / 0.2e1 - t331 / 0.2e1) * t195 + t521 * t789 + t522 * t787 + t489 + (t316 / 0.2e1 + t317 / 0.2e1 + t573 * t479 + t574 * t477 + t528 * t789 + t529 * t787) * t288;
t7 = t113 * t334 + t114 * t333 + t169 * t332 + t170 * t329 + m(6) * (t149 * t169 + t150 * t170 + t676) + m(7) * (t104 * t113 + t105 * t114 + t849) + t827;
t549 = -t6 * qJD(1) + t7 * qJD(2);
t548 = t477 * t564;
t547 = t479 * t563;
t314 = t428 * t420;
t320 = t420 * t433;
t321 = t420 * t437;
t322 = Ifges(7,1) * t666 - t397;
t323 = t441 * t420;
t396 = Ifges(7,6) * t665;
t13 = m(7) * (t104 * t149 + t105 * t150 + t173 * t313) - t313 * t318 + t173 * t314 + t149 * t334 + t150 * t333 + t362 * t315 + t149 * t329 - t150 * t332 - t419 * t396 / 0.2e1 + ((t104 * mrSges(7,2) + t150 * mrSges(6,3) + t615 - t323 / 0.2e1 - t322 / 0.2e1 + t576) * t479 + (t105 * mrSges(7,2) - t149 * mrSges(6,3) + t321 / 0.2e1 - t320 / 0.2e1 - t797) * t477) * t420;
t487 = (t771 + t314 / 0.2e1 + t313 * t786) * t288 + (t641 * t786 + t761 * t800 - t573) * t195 + (t640 * t786 + t758 * t800 + t574) * t196;
t531 = m(7) * (-pkin(5) * t306 + qJ(6) * t307);
t14 = t622 * t307 + t623 * t306 - t531 / 0.2e1 + t487;
t545 = t14 * qJD(1) + t13 * qJD(2);
t45 = m(7) * (t104 * t419 + t173 * t665) - t318 * t665 + t419 * t334;
t56 = (t306 / 0.4e1 + t288 * t590 - t669 / 0.4e1) * t791;
t544 = -qJD(1) * t56 + qJD(2) * t45;
t541 = t707 + t708;
t539 = t704 - t705;
t538 = t702 - t703;
t537 = t313 * t426 + t700;
t53 = t722 + (t670 / 0.2e1 + t654 / 0.4e1 + t646 / 0.4e1 - t150 / 0.4e1) * t791;
t535 = qJD(2) * t53 + qJD(5) * t463;
t532 = pkin(4) * t771 - t426 * t314 / 0.2e1;
t513 = (t195 * t606 + t196 * t605) * pkin(3) + (t633 - t833) * t288;
t488 = (t513 + t817) * t788 + (t513 + t819) * t786 + t831;
t496 = t787 * t879 + t789 * t878;
t17 = t488 + t496 - t831;
t504 = -mrSges(5,2) * t634 + t850 * t633 + t821 * t835;
t515 = t880 * t460;
t59 = -(t404 * t754 + t515) * t784 - m(6) * (t461 * t754 + t515) * pkin(3) - t504;
t482 = (t841 + t538 * t460 + (-t149 * t606 + t150 * t605 + t844) * pkin(3)) * t788 + (t843 + t541 * t460 + (t104 * t605 + t105 * t606 + t173 * t754) * pkin(3)) * t786 - t671 / 0.2e1 - t659 / 0.2e1 + t113 * t609 + t631 / 0.2e1 + t169 * t612 + t629 / 0.2e1 - t619 / 0.2e1 + t618 / 0.2e1 + t332 * t548 + t477 * t333 * t563 + t799 / 0.2e1 + t814 * t633 / 0.2e1 + t812 * t547 - t865;
t485 = -t789 * t856 + t787 * t842 - t748 / 0.2e1 + t664 / 0.2e1 + t107 * t610 + t108 * t613 + t159 * t611 - t630 / 0.2e1 - t813 * t740 / 0.2e1 + (t539 * t789 + t542 * t787 + t653 / 0.2e1 - t652 / 0.2e1) * pkin(10) + t865;
t8 = t482 + t485;
t524 = t17 * qJD(1) + t8 * qJD(2) - t59 * qJD(3);
t492 = t150 * t612 + t104 * t613 - t794 - t479 * t320 / 0.4e1 - t477 * t210 / 0.4e1 + t437 * t665 / 0.4e1 - t434 * t666 / 0.4e1 + t837 * t419 / 0.4e1 + t640 * t609 + t821 * t706 / 0.2e1 + (t323 + t322 + t208) * t477 / 0.4e1 + (t321 + t214 + t212) * t479 / 0.4e1 + (-t552 + t885) * t666 / 0.4e1 + (t433 + t881) * t590;
t484 = t461 * t771 + t314 * t768 + ((-t104 * t477 + t706 + t816) * t786 + t571 * t800 + t795) * t460 + (t313 * t404 + t700) * t786 + t492;
t497 = (-pkin(5) * t108 + qJ(6) * t107) * t787 + t107 * t777 + t108 * t780 - t159 * mrSges(6,1) / 0.2e1 + t160 * t778;
t10 = t484 + (t477 * t620 + t479 * t621) * t419 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t420 + t497 + t801;
t507 = -t718 / 0.2e1 - t717 / 0.2e1 - t715 / 0.2e1 + t714 / 0.2e1;
t495 = (t627 + t507) * t288;
t18 = (t627 + t572) * t288 - t495;
t48 = t404 * t570 + t514 + t658;
t523 = -t18 * qJD(1) + t10 * qJD(2) + t48 * qJD(3);
t327 = (m(7) * t404 + t428) * t477;
t502 = (-t698 + (t404 * t420 + t419 * t460) * t479) * t786 - t586;
t37 = -t753 / 0.2e1 + t502 + t808;
t84 = (t874 + t863) * t750;
t520 = -qJD(1) * t84 - qJD(2) * t37 + qJD(3) * t327;
t519 = -Ifges(7,6) * t668 / 0.2e1 + t477 * t615 - t801 + t853 * t766 + t854 * t592;
t176 = (t428 + (t563 + t768 + t426 / 0.2e1) * m(7)) * t477;
t344 = (m(7) * t426 + t428) * t477;
t503 = (-t698 + (t420 * t426 + t742) * t479) * t786 - t586;
t43 = -t752 / 0.2e1 + t503 + t808;
t518 = -qJD(2) * t43 + qJD(3) * t176 + qJD(4) * t344;
t512 = t432 * t863 + t834 * t288 + (t611 + t612) * t196 + (t609 + t610) * t195;
t493 = (-pkin(5) * t114 + qJ(6) * t113) * t786 + t113 * t776 + t114 * t781 + t169 * t782 + t170 * t779 + t519;
t11 = (-t321 / 0.4e1 + t320 / 0.4e1 - t214 / 0.4e1 - t212 / 0.4e1 + (-t149 / 0.2e1 - t105 / 0.2e1) * mrSges(7,2) + (t640 * t787 - t574) * pkin(10)) * t479 + (-t323 / 0.4e1 - t322 / 0.4e1 + t210 / 0.4e1 - t208 / 0.4e1 + (t104 / 0.2e1 - t150 / 0.2e1) * mrSges(7,2) + (t641 * t787 + t573) * pkin(10)) * t477 + ((t440 / 0.4e1 - t437 / 0.4e1 + t433 / 0.4e1 + t442 / 0.4e1) * t479 + (-t441 / 0.4e1 - t439 / 0.4e1 + t552 / 0.4e1 + t434 / 0.4e1) * t477 - t511 * pkin(10)) * t420 + t537 * t787 + (-t435 / 0.4e1 - t436 / 0.4e1) * t419 + t493 + t532 + t794;
t20 = (t622 * t479 + t623 * t477 + (t744 / 0.4e1 - t711 / 0.4e1 - t430 / 0.4e1) * t791 + t572) * t288;
t490 = (-pkin(5) * t606 + qJ(6) * t605) * t784 / 0.2e1 + t564 * t715 + mrSges(7,3) * t547 + t822 * t548;
t40 = t490 + t884;
t50 = t426 * t570 + t514 - t747;
t509 = t20 * qJD(1) + t11 * qJD(2) + t40 * qJD(3) - t50 * qJD(4);
t505 = (-mrSges(7,2) * t656 - pkin(5) * t464 + t837) * qJD(5);
t499 = qJD(5) * (m(7) * t551 + t809);
t427 = m(7) * t740 + t464;
t403 = m(7) * t660 + t464;
t177 = -t637 * t750 / 0.2e1 + (m(7) * t563 - t428) * t477;
t83 = t872 * t791;
t81 = m(7) * t848;
t57 = (t288 * t665 + t306 + t669) * t786;
t49 = 0.2e1 * t104 * t786 + t334;
t42 = t752 / 0.2e1 + t503 + t796;
t41 = t490 - t884;
t36 = t753 / 0.2e1 + t502 + t796;
t21 = (t626 - t507) * t288 + t512;
t19 = -t495 + t512;
t16 = t688 / 0.2e1 + t689 / 0.2e1 + t728 / 0.2e1 - t729 / 0.2e1 + t488 - t496 + t821 * (-t683 / 0.2e1 - t682 / 0.2e1);
t15 = t531 / 0.2e1 + t487 - t822 * t306 / 0.2e1 + (t779 + t776) * t307;
t12 = t537 * t786 + t492 + t493 - t532 + ((t477 * t641 + t816) * t786 + t795) * pkin(10) + t741 * t870 / 0.2e1;
t9 = t484 - t497 + t519;
t5 = t793 + t561 - t882 + t489 + (t522 + t871) * t786 + (t521 + t886) * t788 + t876;
t3 = t506 + t482 - t485;
t1 = t483 + t486 + t793;
t22 = [qJD(2) * t26 + (qJD(3) + qJD(4)) * t864, t1 * qJD(3) + t5 * qJD(4) + t15 * qJD(5) + t57 * qJD(6) + t690 + (-m(6) * t674 + m(5) * (-t374 * t810 - t674) + t374 * t723 + m(4) * (pkin(8) * t481 * t836 - t755 * pkin(2)) * t476 + (-mrSges(3,1) + t336 + t803 + t857) * t607 - (m(7) * t173 - t719 + t814) * t373 + (m(6) * t150 + m(7) * t104 + t812) * t307 + (-m(6) * t149 + m(7) * t105 - t332 + t333) * t306 + (mrSges(4,3) * t836 - mrSges(3,2)) * t657) * qJD(2), t877 + t1 * qJD(2) + (m(6) * (-t639 + t817) + t688 + m(7) * (-t639 + t819) + t689 + (-t288 * t754 - t500 * t756) * t785 + t728 - t729 + t516 * mrSges(4,2) - t399 * mrSges(4,1) + t821 * (-t683 - t682)) * qJD(3) + t16 * qJD(4) + t19 * qJD(5) + t83 * qJD(6), t877 + t5 * qJD(2) + t16 * qJD(3) + t21 * qJD(5) - t81 * qJD(6) + (t786 * t879 + t788 * t878) * t790 + (t850 * t500 + (mrSges(5,2) - t870) * t288) * qJD(4), t15 * qJD(2) + t19 * qJD(3) + t21 * qJD(4) + (-t822 * t196 + (mrSges(6,2) - mrSges(7,3)) * t195) * qJD(5) + ((-pkin(5) * t196 - qJ(6) * t195) * qJD(5) / 0.2e1 + qJD(6) * t875) * t791, m(7) * t196 * qJD(5) + t57 * qJD(2) + t83 * qJD(3) - t81 * qJD(4); -qJD(3) * t2 - qJD(4) * t6 + qJD(5) * t14 - qJD(6) * t56 - t690, qJD(3) * t4 + qJD(4) * t7 + qJD(5) * t13 + qJD(6) * t45 (-t659 + t799 + t630 + (-t756 * t810 - t844) * t785 + m(7) * (t460 * t542 + t843) + m(6) * (t460 * t539 + t841) + t803 * pkin(8) + t542 * mrSges(7,2) + Ifges(4,5) * t480 - Ifges(4,6) * t478 + t633 * t719 + t634 * t723 + t618 - t619 - t671 + t867 - mrSges(6,3) * t705) * qJD(3) + t3 * qJD(4) + t9 * qJD(5) + t36 * qJD(6) + t550, t3 * qJD(3) + t12 * qJD(5) + t42 * qJD(6) + ((pkin(10) * t541 + t842) * t786 + (pkin(10) * t538 - t856) * t788) * t790 + t549 + (mrSges(7,2) * t708 - mrSges(6,3) * t703 + t629 + t631 - t664 + t748 + (t813 * t479 + (-t330 + t331) * t477) * pkin(10) + t867) * qJD(4), t9 * qJD(3) + t12 * qJD(4) + t49 * qJD(6) + t545 + (-t396 + (-m(7) * pkin(5) - t822) * t150 + (-mrSges(6,2) + t463) * t149 + ((mrSges(7,2) * qJ(6) + Ifges(6,6)) * t479 + (-mrSges(7,2) * pkin(5) + t854) * t477) * t420) * qJD(5), qJD(3) * t36 + qJD(4) * t42 + qJD(5) * t49 + t544; qJD(2) * t2 + qJD(4) * t17 - qJD(5) * t18 + qJD(6) * t84 - t877, qJD(4) * t8 + qJD(5) * t10 + qJD(6) * t37 - t550, -qJD(4) * t59 + qJD(5) * t48 - qJD(6) * t327, t504 * qJD(4) + t41 * qJD(5) + t177 * qJD(6) + ((t426 * t633 + t636) * t786 + (-pkin(4) * t633 + t636) * t788) * t790 + t524, t41 * qJD(4) + t403 * qJD(6) + t460 * t499 + t505 + t523, qJD(4) * t177 + qJD(5) * t403 - t520; qJD(2) * t6 - qJD(3) * t17 - qJD(5) * t20 - t877, -qJD(3) * t8 - qJD(5) * t11 + qJD(6) * t43 - t549, -qJD(5) * t40 - qJD(6) * t176 - t524, qJD(5) * t50 - qJD(6) * t344, pkin(10) * t499 + t427 * qJD(6) + t505 - t509, qJD(5) * t427 - t518; -qJD(2) * t14 + qJD(3) * t18 + qJD(4) * t20, -qJD(3) * t10 + qJD(4) * t11 + qJD(6) * t53 - t545, qJD(4) * t40 - t523, t509, t463 * qJD(6), t535; qJD(2) * t56 - qJD(3) * t84, -qJD(3) * t37 - qJD(4) * t43 - qJD(5) * t53 - t544, qJD(4) * t176 + t520, t518, -t535, 0;];
Cq  = t22;
