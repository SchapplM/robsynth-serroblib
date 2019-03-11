% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:53:08
% EndTime: 2019-03-09 21:53:44
% DurationCPUTime: 23.01s
% Computational Cost: add. (65812->759), mult. (126780->981), div. (0->0), fcn. (154799->10), ass. (0->429)
t469 = sin(qJ(4));
t471 = sin(qJ(2));
t473 = cos(qJ(4));
t470 = sin(qJ(3));
t474 = cos(qJ(3));
t475 = cos(qJ(2));
t444 = -t470 * t471 + t474 * t475;
t782 = t444 * pkin(9);
t464 = t475 * pkin(7);
t644 = t475 * pkin(8) + t464;
t844 = t474 * t644;
t525 = t844 + t782;
t445 = -t470 * t475 - t474 * t471;
t781 = t445 * pkin(9);
t845 = t470 * t644;
t526 = -t845 + t781;
t804 = pkin(8) + pkin(7);
t631 = t474 * t804;
t632 = t470 * t804;
t282 = -t469 * t525 + t473 * t526 + (t469 * t632 - t473 * t631) * t471;
t630 = t804 * t471;
t854 = -t474 * t630 - t845;
t371 = t854 + t781;
t405 = -t470 * t630 + t844;
t499 = t405 + t782;
t920 = t473 * t371 - t469 * t499;
t951 = -t282 + t920;
t941 = t920 * mrSges(5,2);
t626 = -t941 / 0.2e1;
t467 = sin(pkin(11));
t395 = t444 * t469 - t445 * t473;
t388 = t395 * qJ(5);
t937 = -t388 + t920;
t950 = t467 * t937;
t472 = cos(qJ(6));
t570 = t473 * t444 + t445 * t469;
t708 = cos(pkin(11));
t318 = t395 * t467 - t570 * t708;
t468 = sin(qJ(6));
t720 = t468 * mrSges(7,3);
t851 = t395 * t708 + t467 * t570;
t888 = -mrSges(7,2) * t851 + t318 * t720;
t923 = t472 * t888;
t949 = t708 * t937;
t780 = t473 * pkin(3);
t458 = pkin(4) + t780;
t581 = t708 * t469;
t428 = pkin(3) * t581 + t467 * t458;
t420 = pkin(10) + t428;
t717 = t472 * mrSges(7,3);
t889 = mrSges(7,1) * t851 + t318 * t717;
t925 = t468 * t889;
t943 = -t925 / 0.2e1;
t857 = t923 / 0.2e1 + t943;
t948 = t857 * t420;
t767 = Ifges(7,4) * t472;
t449 = Ifges(7,1) * t468 + t767;
t647 = t472 * t449;
t768 = Ifges(7,4) * t468;
t448 = Ifges(7,2) * t472 + t768;
t660 = t468 * t448;
t838 = -t660 / 0.2e1 + t647 / 0.2e1;
t847 = Ifges(5,5) * t570;
t863 = Ifges(5,6) * t395;
t877 = t851 * (Ifges(7,5) * t468 + Ifges(7,6) * t472);
t882 = Ifges(6,6) * t851;
t908 = Ifges(6,5) * t318;
t761 = Ifges(7,2) * t468;
t552 = -t761 + t767;
t892 = Ifges(7,6) * t851 - t318 * t552;
t922 = t472 * t892;
t553 = Ifges(7,1) * t472 - t768;
t891 = Ifges(7,5) * t851 - t318 * t553;
t924 = t468 * t891;
t931 = t847 - t863 + t877 / 0.2e1 - t882 - t908 + t922 / 0.2e1 + t924 / 0.2e1;
t918 = -t405 * mrSges(4,1) - t854 * mrSges(4,2) + Ifges(4,5) * t444 + Ifges(4,6) * t445 - t838 * t318 + t931;
t947 = t918 - t941;
t946 = t847 / 0.2e1 - t863 / 0.2e1 + t877 / 0.4e1 - t882 / 0.2e1 - t908 / 0.2e1 + t922 / 0.4e1 + t924 / 0.4e1 - t318 * (t647 / 0.4e1 - t660 / 0.4e1);
t786 = pkin(2) * t474;
t459 = pkin(3) + t786;
t657 = t470 * t473;
t430 = pkin(2) * t657 + t459 * t469;
t411 = t708 * t430;
t658 = t469 * t470;
t429 = -pkin(2) * t658 + t473 * t459;
t424 = pkin(4) + t429;
t363 = t467 * t424 + t411;
t352 = pkin(10) + t363;
t740 = t851 * mrSges(6,3);
t899 = t363 * t740;
t670 = t467 * t430;
t362 = t424 * t708 - t670;
t743 = t318 * mrSges(6,3);
t927 = t362 * t743;
t938 = t352 * t923;
t945 = t352 * t943 + t938 / 0.2e1 + t927 / 0.2e1 - t899 / 0.2e1;
t284 = (-t469 * t631 - t473 * t632) * t471 + t473 * t525 + t469 * t526;
t929 = t284 * mrSges(5,1);
t942 = -t929 / 0.2e1;
t940 = t430 * t920;
t939 = t469 * t920;
t555 = mrSges(7,1) * t472 - mrSges(7,2) * t468;
t229 = t282 - t388;
t280 = -t371 * t469 - t473 * t499;
t680 = t570 * qJ(5);
t538 = t280 - t680;
t901 = t538 * t708;
t827 = t229 * t467 - t901;
t704 = t827 * t555;
t747 = t282 * mrSges(5,2);
t903 = t467 * t538;
t126 = -t708 * t229 - t903;
t752 = t126 * mrSges(6,2);
t755 = t827 * mrSges(6,1);
t936 = -t704 / 0.2e1 + t942 - t747 / 0.2e1 + t752 / 0.2e1 - t755 / 0.2e1 + t946;
t935 = t923 - t925;
t610 = t708 * pkin(4);
t569 = mrSges(6,3) * t610;
t783 = pkin(4) * t467;
t642 = mrSges(6,3) * t783;
t802 = -t851 / 0.2e1;
t455 = -t610 - pkin(5);
t869 = t455 / 0.2e1;
t871 = t318 / 0.2e1;
t718 = t472 * mrSges(7,2);
t721 = t468 * mrSges(7,1);
t554 = t718 + t721;
t900 = t554 * t318;
t934 = -t569 * t871 - t642 * t802 + t900 * t869 - t946;
t933 = t951 * t570;
t654 = t472 * t891;
t930 = -t654 / 0.2e1;
t928 = t284 * t395;
t466 = t472 ^ 2;
t722 = t466 * mrSges(7,3);
t465 = t468 ^ 2;
t723 = t465 * mrSges(7,3);
t897 = -t722 / 0.2e1 - t723 / 0.2e1;
t921 = t897 * t318;
t880 = t851 * mrSges(6,1);
t907 = t318 * mrSges(6,2);
t580 = -t907 + t880;
t531 = t455 * t554;
t787 = t468 / 0.2e1;
t868 = t472 / 0.2e1;
t563 = t552 * t868 + t553 * t787 + t838;
t544 = t531 / 0.2e1 + t563;
t534 = t126 * t554;
t757 = Ifges(7,6) * t468;
t763 = Ifges(7,5) * t472;
t551 = -t757 + t763;
t911 = -t318 / 0.2e1;
t916 = 0.2e1 * t551 * t911 - t534 / 0.2e1;
t914 = t888 * t787 + t889 * t868;
t460 = -pkin(2) * t475 - pkin(1);
t413 = -t444 * pkin(3) + t460;
t328 = -pkin(4) * t570 + t413;
t159 = t318 * pkin(5) - pkin(10) * t851 + t328;
t69 = t159 * t472 - t468 * t827;
t70 = t159 * t468 + t472 * t827;
t913 = t126 * t900 - t69 * t889 - t70 * t888;
t912 = t318 ^ 2;
t910 = t880 / 0.2e1;
t840 = t318 * t551;
t896 = pkin(5) * t851 + pkin(10) * t318;
t895 = t405 * mrSges(4,3) + Ifges(4,4) * t445;
t528 = t553 * t851;
t764 = Ifges(7,5) * t318;
t500 = t472 * (t528 + t764);
t527 = t552 * t851;
t758 = Ifges(7,6) * t318;
t501 = t468 * (t527 + t758);
t616 = t720 / 0.2e1;
t617 = -t720 / 0.2e1;
t684 = t851 * t472;
t685 = t851 * t468;
t893 = -t318 * (-t757 / 0.2e1 + t763 / 0.2e1) + t534 / 0.2e1 - t501 / 0.4e1 + t500 / 0.4e1 - t449 * t685 / 0.2e1 - t448 * t684 / 0.2e1 - t468 * t527 / 0.4e1 + t472 * t528 / 0.4e1 + t840 / 0.4e1 + (t616 + t617) * t70;
t799 = t851 / 0.2e1;
t864 = Ifges(6,4) * t851;
t870 = -t395 / 0.2e1;
t867 = t570 / 0.2e1;
t784 = pkin(4) * t395;
t865 = Ifges(4,2) - Ifges(4,1);
t769 = Ifges(5,4) * t395;
t365 = t429 * t708 - t670;
t669 = t467 * t469;
t434 = (t473 * t708 - t669) * pkin(3);
t860 = t365 + t434;
t462 = t471 * pkin(2);
t785 = pkin(3) * t445;
t558 = -t785 + t784;
t330 = t462 + t558;
t564 = t910 - t907 / 0.2e1 + t914;
t481 = -t284 - t680;
t129 = t467 * t481 + t949;
t174 = t558 + t896;
t160 = t174 + t462;
t71 = -t129 * t468 + t160 * t472;
t72 = t129 * t472 + t160 * t468;
t813 = m(7) / 0.2e1;
t815 = m(6) / 0.2e1;
t859 = -t330 * t815 - (t468 * t72 + t472 * t71) * t813 - t564;
t181 = t784 + t896;
t81 = t126 * t468 + t181 * t472;
t82 = -t126 * t472 + t181 * t468;
t547 = -t468 * t81 + t472 * t82;
t212 = t554 * t851;
t427 = -pkin(3) * t669 + t458 * t708;
t419 = -pkin(5) - t427;
t614 = t717 / 0.2e1;
t494 = t82 * t614 + t81 * t617 + t936;
t215 = -mrSges(7,2) * t318 - mrSges(7,3) * t685;
t651 = t472 * t215;
t218 = mrSges(7,1) * t318 - mrSges(7,3) * t684;
t661 = t468 * t218;
t536 = t651 / 0.2e1 - t661 / 0.2e1;
t550 = -t468 * t69 + t472 * t70;
t433 = (t467 * t473 + t581) * pkin(3);
t695 = t126 * t433;
t789 = -t434 / 0.2e1;
t790 = t433 / 0.2e1;
t791 = -t428 / 0.2e1;
t792 = -t427 / 0.2e1;
t793 = t419 / 0.2e1;
t855 = t212 * t790 - t900 * t793 + (-t126 * t428 + t695 + (-t427 + t434) * t827) * t815 + (t419 * t827 + t420 * t547 + t434 * t550 + t695) * t813 + t494 + t948 + ((t790 + t791) * t851 + (t789 - t792) * t318) * mrSges(6,3) + t536 * t434;
t537 = t215 * t787 + t218 * t868;
t435 = (-t469 * t474 - t657) * pkin(2);
t436 = (t473 * t474 - t658) * pkin(2);
t377 = -t435 * t708 + t467 * t436;
t344 = t377 * t555;
t378 = t467 * t435 + t436 * t708;
t366 = t378 * t723;
t367 = t378 * t722;
t376 = t377 * mrSges(6,1);
t431 = t435 * mrSges(5,1);
t725 = t436 * mrSges(5,2);
t735 = t378 * mrSges(6,2);
t852 = -(mrSges(4,1) * t470 + mrSges(4,2) * t474) * pkin(2) - t344 + t366 + t367 - t376 + t431 - t725 - t735;
t646 = t473 * t570;
t659 = t469 * t395;
t849 = (-t646 / 0.2e1 - t659 / 0.2e1) * pkin(3);
t560 = mrSges(7,3) * (t466 / 0.2e1 + t465 / 0.2e1);
t770 = Ifges(5,4) * t570;
t709 = -mrSges(6,1) - t555;
t842 = t851 * t551;
t128 = t903 + t949;
t76 = t128 * t472 + t174 * t468;
t805 = -t76 / 0.2e1;
t75 = -t128 * t468 + t174 * t472;
t806 = t75 / 0.2e1;
t835 = mrSges(7,1) * t806 + mrSges(7,2) * t805;
t808 = -mrSges(7,2) / 0.2e1;
t809 = mrSges(7,1) / 0.2e1;
t834 = t71 * t809 + t72 * t808;
t548 = -t75 * t468 + t76 * t472;
t713 = t72 * t472;
t714 = t71 * t468;
t549 = t713 - t714;
t831 = -t743 / 0.2e1 + t536;
t830 = t435 * t870 + t436 * t867;
t829 = (-mrSges(5,1) * t469 - mrSges(5,2) * t473) * pkin(3);
t125 = -t481 * t708 + t950;
t817 = m(5) / 0.2e1;
t828 = pkin(3) * (-t284 * t473 + t939) * t817 + (-t125 * t427 + t129 * t428) * t815 + (t125 * t419 + t420 * t549) * t813 + t626;
t825 = t851 * t897 - t537;
t811 = m(6) * pkin(4);
t824 = m(7) * t455 - t708 * t811;
t750 = t129 * mrSges(6,2);
t823 = t942 - t750 / 0.2e1 + (t713 / 0.2e1 - t714 / 0.2e1) * mrSges(7,3);
t124 = -t901 + t950;
t703 = t124 * t555;
t749 = t280 * mrSges(5,1);
t751 = t128 * mrSges(6,2);
t754 = t124 * mrSges(6,1);
t822 = -t703 / 0.2e1 + t749 / 0.2e1 - t751 / 0.2e1 - t754 / 0.2e1;
t454 = pkin(10) + t783;
t643 = t465 + t466;
t572 = t643 * t454;
t821 = m(7) * t572 + t467 * t811 + t722 + t723;
t816 = -m(6) / 0.2e1;
t814 = -m(7) / 0.2e1;
t812 = m(5) * pkin(3);
t810 = -mrSges(7,1) / 0.2e1;
t807 = -Ifges(7,1) / 0.2e1;
t351 = -pkin(5) - t362;
t798 = -t351 / 0.2e1;
t797 = t351 / 0.2e1;
t364 = t429 * t467 + t411;
t796 = t364 / 0.2e1;
t795 = -t377 / 0.2e1;
t794 = t378 / 0.2e1;
t779 = -Ifges(6,2) - Ifges(7,3);
t772 = Ifges(5,1) * t395;
t762 = Ifges(5,2) * t570;
t756 = Ifges(7,3) * t851;
t753 = t125 * mrSges(6,1);
t737 = t364 * mrSges(6,1);
t736 = t365 * mrSges(6,2);
t219 = mrSges(6,1) * t318 + mrSges(6,2) * t851;
t323 = -mrSges(5,1) * t570 + mrSges(5,2) * t395;
t418 = t462 - t785;
t439 = Ifges(4,4) * t444;
t724 = t445 * mrSges(4,3);
t508 = t827 * t740 - t328 * t580 + mrSges(5,3) * t928 - t413 * (mrSges(5,1) * t395 + mrSges(5,2) * t570) - t405 * t724 - t460 * (-mrSges(4,1) * t445 + mrSges(4,2) * t444) + t913;
t668 = t468 * t892;
t698 = t126 * t318;
t702 = t125 * t126;
t4 = (-mrSges(4,2) * t462 + t444 * t865 - t895) * t445 + (-pkin(2) * mrSges(4,1) * t444 - pkin(1) * mrSges(3,1) - Ifges(3,4) * t471) * t471 + (t125 * t851 - t129 * t318 - t698) * mrSges(6,3) + t444 * t439 + m(6) * (t129 * t827 + t328 * t330 + t702) + (t928 + t933) * mrSges(5,3) + m(5) * (t284 * t951 + t413 * t418) + (t654 + t842 - t864 + (-Ifges(6,1) + Ifges(7,3)) * t318) * t799 + (-Ifges(5,2) * t395 + 0.2e1 * t770 + t772) * t867 + (t762 + t769) * t870 + t418 * t323 + t330 * t219 + t125 * t212 + t72 * t215 + t71 * t218 + (-Ifges(6,2) * t318 + t668 + t864) * t802 - t508 + m(7) * (t69 * t71 + t70 * t72 + t702) + t395 * (Ifges(5,1) * t570 - t769) / 0.2e1 + m(4) * t460 * t462 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t475 + (Ifges(3,1) - Ifges(3,2)) * t471) * t475 + (t756 - t840 + t501) * t871 + (t500 - 0.2e1 * Ifges(6,4) * t318 + (Ifges(6,1) - Ifges(6,2)) * t851) * t911;
t734 = t4 * qJD(1);
t729 = t429 * mrSges(5,2);
t728 = t430 * mrSges(5,1);
t727 = t433 * mrSges(6,1);
t726 = t434 * mrSges(6,2);
t497 = t466 * t807 - Ifges(6,1) + (t767 - t761 / 0.2e1) * t468;
t700 = t126 * t124;
t7 = (t851 ^ 2 - t912) * Ifges(6,4) + (t395 ^ 2 - t570 ^ 2) * Ifges(5,4) + mrSges(6,3) * t698 + t128 * t743 - t558 * t219 - t76 * t215 - t75 * t218 - t124 * t212 + (-t445 * t865 - t439) * t444 + (t779 * t851 + t840) * t318 + (pkin(3) * t323 + t895) * t445 + t508 - t570 * t772 + t395 * t762 - m(5) * (t282 * t280 + t284 * t920 - t413 * t785) - m(6) * (t128 * t827 + t328 * t558 + t700) - m(7) * (t69 * t75 + t70 * t76 + t700) + (t930 + t668 / 0.2e1 - t124 * mrSges(6,3) + t551 * t802 - t497 * t318) * t851 + (t280 * t395 - t933) * mrSges(5,3);
t715 = t7 * qJD(1);
t9 = -t81 * t218 - t82 * t215 - t827 * t212 - m(7) * (t827 * t126 + t69 * t81 + t70 * t82) + (-t413 * mrSges(5,2) - t770) * t570 + (-t413 * mrSges(5,1) - pkin(4) * t219 + t769 + (-Ifges(5,1) + Ifges(5,2)) * t570) * t395 + (-Ifges(6,4) + t551) * t912 + (t930 + t892 * t787 - t842 / 0.2e1 + t864 + (-t497 + t779) * t318) * t851 + (-m(6) * t784 - t580) * t328 + t913;
t710 = t9 * qJD(1);
t699 = t126 * t851;
t17 = (t212 + t740) * t851 - (t651 - t661 - t743) * t318 + m(7) * (-t318 * t550 + t699) + m(6) * (-t318 * t827 + t699);
t707 = qJD(1) * t17;
t10 = t69 * t215 - t70 * t218 + ((mrSges(7,1) * t126 - mrSges(7,3) * t70 - Ifges(7,4) * t684 - t758) * t472 + (-t126 * mrSges(7,2) + t69 * mrSges(7,3) - t764 + (t768 + (-Ifges(7,1) + Ifges(7,2)) * t472) * t851) * t468) * t851;
t706 = t10 * qJD(1);
t701 = t125 * t555;
t697 = t126 * t364;
t696 = t126 * t377;
t688 = t851 * t555;
t611 = -t688 / 0.2e1 + t921;
t545 = -t880 / 0.2e1 + t907 / 0.2e1 + t611;
t575 = t643 * t352;
t487 = (-t318 * t363 - t362 * t851) * t815 + (-t318 * t575 + t351 * t851) * t813 + t545;
t21 = t487 + t859;
t694 = t21 * qJD(1);
t491 = t558 * t815 + (t468 * t76 + t472 * t75) * t813 + t564;
t574 = t643 * t420;
t492 = (-t318 * t428 - t427 * t851) * t816 + (-t318 * t574 + t419 * t851) * t814;
t609 = t688 / 0.2e1;
t25 = t609 + t910 - (mrSges(6,2) / 0.2e1 - t560) * t318 + t491 + t492;
t692 = t25 * qJD(1);
t641 = t811 / 0.2e1;
t490 = (-t318 * t572 + t455 * t851) * t813 + (-t318 * t467 - t708 * t851) * t641;
t493 = (t468 * t82 + t472 * t81) * t813 + t395 * t641 + t914;
t27 = -t490 + t493 + t580 + t609 - t921;
t691 = t27 * qJD(1);
t539 = -t718 / 0.2e1 - t721 / 0.2e1;
t521 = t539 * t318;
t31 = -t521 - t536;
t689 = t31 * qJD(1);
t683 = t351 * t900;
t682 = t364 * t555;
t679 = t419 * t900;
t678 = t427 * t318;
t677 = t428 * t851;
t676 = t429 * t570;
t675 = t430 * t395;
t674 = t433 * t555;
t640 = mrSges(5,3) * t676;
t639 = mrSges(5,3) * t675;
t638 = -pkin(3) * t469 / 0.2e1;
t637 = -t780 / 0.2e1;
t629 = t756 / 0.2e1;
t624 = t740 / 0.2e1;
t615 = -t717 / 0.2e1;
t586 = t420 / 0.2e1 - t352 / 0.2e1;
t573 = t643 * t434;
t559 = t941 / 0.2e1;
t557 = mrSges(7,3) * t643 - mrSges(6,2);
t556 = -t728 - t729;
t476 = (-t126 * t363 - t362 * t827 + t697) * t815 + (t351 * t827 + t352 * t547 + t697) * t813 - t900 * t797 + t212 * t796 + t364 * t624 + (t550 * t813 + t815 * t827 + t831) * t365 + t945;
t495 = t454 * t857 - t934;
t477 = t495 + (t125 * t455 + t454 * t549) * t813 + (-t125 * t708 + t129 * t467) * t641 - t701 / 0.2e1 + t626 - t753 / 0.2e1 + t823;
t1 = t82 * t615 + t81 * t616 - t476 + t477 - t936;
t48 = t709 * t364 + t557 * t365 + m(7) * (t351 * t364 + t365 * t575) + m(6) * (-t362 * t364 + t363 * t365) + t556;
t543 = -t1 * qJD(1) + t48 * qJD(2);
t49 = -m(7) * (t351 * t377 + t378 * t575) - m(6) * (-t362 * t377 + t363 * t378) - m(5) * (t429 * t435 + t430 * t436) - t852;
t479 = -m(5) * (t429 * t280 + t435 * t282 + t436 * t284 + t940) / 0.2e1 + (-t124 * t362 + t128 * t363 + t378 * t827 + t696) * t816 + (t124 * t351 + t352 * t548 + t378 * t550 + t696) * t814 + t212 * t795;
t8 = t479 + ((t795 + t791 + t363 / 0.2e1) * t851 + (t794 - t792 - t362 / 0.2e1) * t318) * mrSges(6,3) - (t793 + t798) * t900 + (t218 * t794 - t586 * t889 + (t806 - t71 / 0.2e1) * mrSges(7,3)) * t468 + (t128 / 0.2e1 - t129 / 0.2e1) * mrSges(6,2) + t559 + (-t280 / 0.2e1 - t284 / 0.2e1) * mrSges(5,1) + (-t378 * t215 / 0.2e1 + t586 * t888 + (t805 + t72 / 0.2e1) * mrSges(7,3)) * t472 + (t676 / 0.2e1 + t675 / 0.2e1 + t849 - t830) * mrSges(5,3) + t709 * (t125 / 0.2e1 - t124 / 0.2e1) + t828;
t542 = -t8 * qJD(1) - t49 * qJD(2);
t486 = t629 + t916;
t505 = (Ifges(7,1) / 0.4e1 - Ifges(7,2) / 0.2e1) * t468 + t449 / 0.4e1 + 0.3e1 / 0.2e1 * t767;
t524 = (t807 + Ifges(7,2) / 0.4e1) * t472 + t448 / 0.4e1;
t11 = t537 * t352 + (t352 * t560 + (mrSges(7,1) * t798 + t524) * t472 + (mrSges(7,2) * t797 + t505) * t468) * t851 + t486 + t834;
t533 = t351 * t554;
t244 = t533 + t563;
t541 = -t11 * qJD(1) + t244 * qJD(2);
t532 = t419 * t554;
t522 = t76 * t614 + t75 * t617 + t626 + t822;
t520 = t539 * t365;
t519 = t539 * t378;
t518 = t539 * t434;
t515 = -t533 / 0.2e1;
t514 = -t532 / 0.2e1;
t511 = t555 * t799;
t109 = t709 * t433 + t829 + t557 * t434 + m(7) * (t419 * t433 + t420 * t573) + m(6) * (-t427 * t433 + t428 * t434);
t480 = -t344 / 0.2e1 + t366 / 0.2e1 + t367 / 0.2e1 - t376 / 0.2e1 + t431 / 0.2e1 + (t377 * t455 + t378 * t572) * t813 - t735 / 0.2e1 - t725 / 0.2e1 + (-t377 * t708 + t378 * t467) * t641;
t482 = (-t362 * t433 + t363 * t434 - t364 * t427 + t365 * t428) * t815 + (t351 * t433 + t352 * t573 + t364 * t419 + t365 * t574) * t813;
t30 = -t480 + mrSges(5,2) * t637 + mrSges(5,1) * t638 - t674 / 0.2e1 - t682 / 0.2e1 - t727 / 0.2e1 - t726 / 0.2e1 - t729 / 0.2e1 - t728 / 0.2e1 - t737 / 0.2e1 - t736 / 0.2e1 + t482 - t860 * t897;
t489 = (t124 * t455 + t454 * t548) * t813 + (-t124 * t708 + t128 * t467) * t641;
t6 = -t489 + t75 * t616 + t559 + t76 * t615 + (t925 / 0.2e1 - t923 / 0.2e1) * t454 - t822 + t855 + t934;
t510 = t6 * qJD(1) + t30 * qJD(2) + t109 * qJD(3);
t13 = t537 * t420 + (t420 * t560 + (t419 * t810 + t524) * t472 + (mrSges(7,2) * t793 + t505) * t468) * t851 + t486 + t835;
t153 = t515 + t514 + t519 - t563;
t292 = t532 + t563;
t509 = -t13 * qJD(1) - t153 * qJD(2) + t292 * qJD(3);
t507 = Ifges(7,3) * t799 + t82 * t808 + t81 * t809;
t15 = t537 * t454 + (t454 * t560 + (t455 * t810 + t524) * t472 + (mrSges(7,2) * t869 + t505) * t468) * t851 + t507 + t916;
t155 = t515 + t520 - t544;
t221 = t514 + t518 - t544;
t322 = t531 + t563;
t502 = t15 * qJD(1) + t155 * qJD(2) + t221 * qJD(3) - t322 * qJD(4);
t485 = t629 + t893;
t396 = t532 / 0.2e1;
t327 = t533 / 0.2e1;
t222 = t396 + t518 + t544;
t156 = t327 + t520 + t544;
t154 = t327 + t396 + t519 + t563;
t32 = -t521 + t536;
t29 = t480 + (t637 - t429 / 0.2e1) * mrSges(5,2) + (t789 - t365 / 0.2e1) * mrSges(6,2) + (t638 - t430 / 0.2e1) * mrSges(5,1) + t482 + t709 * (t790 + t796) + t860 * t560;
t28 = t490 + t493 + t611;
t26 = t491 - t492 + t545;
t20 = t487 - t859;
t16 = t454 * t825 + t455 * t511 + t507 + t893;
t14 = t419 * t511 + t420 * t825 + t485 + t835;
t12 = t351 * t511 + t352 * t825 + t485 + t834;
t5 = t489 + t495 + t522 + t855;
t3 = t828 + t831 * t378 + t823 - t639 / 0.2e1 - t640 / 0.2e1 - t479 - t679 / 0.2e1 + t522 + (-t677 / 0.2e1 + t678 / 0.2e1) * mrSges(6,3) + (-t555 / 0.2e1 - mrSges(6,1) / 0.2e1) * t125 - t683 / 0.2e1 + t377 * t624 + t948 + (t830 + t849) * mrSges(5,3) + t918 + t945;
t2 = t477 + t476 + t494;
t18 = [qJD(2) * t4 - qJD(3) * t7 - qJD(4) * t9 + qJD(5) * t17 + qJD(6) * t10, t734 + (-t444 * mrSges(4,3) * t786 + t470 * pkin(2) * t724 - mrSges(3,1) * t464 + Ifges(3,5) * t475 + t938 - t899 + t927 - t352 * t925 - t639 - t640 - t683 - t701 - t929 - t750 - t753 + (mrSges(3,2) * pkin(7) - Ifges(3,6)) * t471 + t549 * mrSges(7,3) + t947) * qJD(2) + t3 * qJD(3) + t2 * qJD(4) + t20 * qJD(5) + t12 * qJD(6) + 0.2e1 * ((t125 * t351 + t352 * t549) * t813 + (-t284 * t429 + t940) * t817 + (-t125 * t362 + t129 * t363) * t815 + m(4) * (-t405 * t474 + t470 * t854) * pkin(2) / 0.2e1) * qJD(2), -t715 + t3 * qJD(2) + ((t473 * t280 + t939) * t812 - t703 - t679 + t749 - t751 - t754 + m(6) * (-t124 * t427 + t128 * t428) + m(7) * t124 * t419 + (-t646 - t659) * mrSges(5,3) * pkin(3) + (m(7) * t548 + t935) * t420 + t548 * mrSges(7,3) + (-t677 + t678) * mrSges(6,3) + t947) * qJD(3) + t5 * qJD(4) + t26 * qJD(5) + t14 * qJD(6), t2 * qJD(2) + t5 * qJD(3) + t28 * qJD(5) + t16 * qJD(6) - t710 + (-t704 + t752 - t755 - t929 - t747 + (-t126 * t467 - t708 * t827) * t811 - t851 * t642 - (-t569 + t838) * t318 + (m(7) * t827 - t900) * t455 + (m(7) * t547 + t935) * t454 + t547 * mrSges(7,3) + t931) * qJD(4), qJD(2) * t20 + qJD(3) * t26 + qJD(4) * t28 + qJD(6) * t32 + t707, t706 + t12 * qJD(2) + t14 * qJD(3) + t16 * qJD(4) + t32 * qJD(5) + (-mrSges(7,1) * t70 - mrSges(7,2) * t69 - t877) * qJD(6); -qJD(3) * t8 - qJD(4) * t1 + qJD(5) * t21 - qJD(6) * t11 - t734, -qJD(3) * t49 + qJD(4) * t48 + qJD(6) * t244, t29 * qJD(4) + t154 * qJD(6) + t542 + (0.2e1 * (t377 * t419 + t378 * t574) * t813 + 0.2e1 * (-t377 * t427 + t378 * t428) * t815 + (t435 * t473 + t436 * t469) * t812 + t852) * qJD(3), t29 * qJD(3) + (t364 * t824 + t365 * t821 + t556 - t682 - t736 - t737) * qJD(4) + t156 * qJD(6) + t543, t694, t154 * qJD(3) + t156 * qJD(4) + (-t352 * t555 + t551) * qJD(6) + t541; qJD(2) * t8 + qJD(4) * t6 - qJD(5) * t25 - qJD(6) * t13 + t715, qJD(4) * t30 - qJD(6) * t153 - t542, qJD(4) * t109 + qJD(6) * t292 (t433 * t824 + t434 * t821 - t674 - t726 - t727 + t829) * qJD(4) + t222 * qJD(6) + t510, -t692, t222 * qJD(4) + (-t420 * t555 + t551) * qJD(6) + t509; qJD(2) * t1 - qJD(3) * t6 - qJD(5) * t27 - qJD(6) * t15 + t710, -qJD(3) * t30 - qJD(6) * t155 - t543, -qJD(6) * t221 - t510, t322 * qJD(6), -t691 (-t454 * t555 + t551) * qJD(6) - t502; -qJD(2) * t21 + qJD(3) * t25 + qJD(4) * t27 - qJD(6) * t31 - t707, -t694, t692, t691, 0, -qJD(6) * t554 - t689; qJD(2) * t11 + qJD(3) * t13 + qJD(4) * t15 + qJD(5) * t31 - t706, qJD(3) * t153 + qJD(4) * t155 - t541, qJD(4) * t221 - t509, t502, t689, 0;];
Cq  = t18;
