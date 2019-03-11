% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRP8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:22
% EndTime: 2019-03-09 17:14:55
% DurationCPUTime: 17.47s
% Computational Cost: add. (20296->935), mult. (41577->1209), div. (0->0), fcn. (39575->6), ass. (0->449)
t530 = sin(qJ(5));
t533 = cos(qJ(3));
t531 = sin(qJ(3));
t774 = cos(qJ(5));
t634 = t774 * t531;
t454 = -t530 * t533 + t634;
t563 = t530 * t531 + t533 * t774;
t259 = mrSges(7,1) * t563 + mrSges(7,2) * t454;
t519 = t531 * qJ(4);
t579 = t533 * pkin(3) + t519;
t470 = -pkin(2) - t579;
t434 = t533 * pkin(4) - t470;
t278 = pkin(5) * t563 + t434;
t772 = m(7) * t278;
t613 = t259 + t772;
t260 = mrSges(6,1) * t563 + mrSges(6,2) * t454;
t898 = m(6) * t434 + t260;
t765 = pkin(8) * t531;
t482 = -t531 * pkin(9) + t765;
t483 = (pkin(8) - pkin(9)) * t533;
t307 = -t774 * t482 + t530 * t483;
t564 = t530 * t482 + t483 * t774;
t550 = -qJ(6) * t563 + t564;
t591 = -qJ(6) * t454 - t307;
t897 = -t564 * mrSges(6,1) - t550 * mrSges(7,1) + t307 * mrSges(6,2) - t591 * mrSges(7,2);
t817 = -t307 / 0.2e1;
t565 = t530 * t307 + t564 * t774;
t566 = -t530 * t591 + t550 * t774;
t532 = sin(qJ(2));
t681 = t533 * t532;
t500 = qJ(4) * t681;
t535 = -pkin(3) - pkin(4);
t653 = t535 * t531;
t592 = -pkin(7) + t653;
t302 = t532 * t592 + t500;
t396 = t530 * t681 - t532 * t634;
t181 = pkin(5) * t396 + t302;
t397 = t563 * t532;
t208 = mrSges(6,1) * t397 - mrSges(6,2) * t396;
t534 = cos(qJ(2));
t364 = t396 * mrSges(7,2);
t586 = -t397 * mrSges(7,1) + t364;
t365 = Ifges(7,6) * t397;
t366 = Ifges(6,6) * t397;
t367 = Ifges(7,5) * t396;
t368 = Ifges(6,5) * t396;
t856 = -t367 - t365 - t368 - t366;
t893 = -t181 * t586 + t302 * t208 + t856 * t534 / 0.2e1;
t258 = mrSges(6,1) * t454 - mrSges(6,2) * t563;
t314 = mrSges(6,1) * t534 - t397 * mrSges(6,3);
t813 = -t314 / 0.2e1;
t313 = mrSges(7,1) * t534 - t397 * mrSges(7,3);
t814 = -t313 / 0.2e1;
t722 = t396 * mrSges(7,3);
t311 = -mrSges(7,2) * t534 - t722;
t816 = t311 / 0.2e1;
t436 = t563 * mrSges(7,2);
t717 = t454 * mrSges(7,1);
t852 = -t717 + t436;
t892 = -t278 * t586 / 0.2e1 + t591 * t816 - t181 * t852 / 0.2e1 + t550 * t814 + t564 * t813 + t434 * t208 / 0.2e1 + t302 * t258 / 0.2e1;
t588 = t533 * mrSges(5,1) + t531 * mrSges(5,3);
t789 = t588 / 0.2e1;
t756 = mrSges(6,1) + mrSges(7,1);
t755 = mrSges(6,2) + mrSges(7,2);
t883 = mrSges(6,3) + mrSges(7,3);
t764 = pkin(8) * t532;
t770 = pkin(2) * t534;
t471 = -pkin(1) - t764 - t770;
t682 = t531 * t534;
t670 = pkin(7) * t682 - t533 * t471;
t284 = pkin(9) * t681 - t670;
t680 = t533 * t534;
t361 = pkin(7) * t680 + t531 * t471;
t683 = t531 * t532;
t505 = pkin(9) * t683;
t285 = t505 + t361;
t118 = t774 * t284 + t530 * t285;
t529 = t534 * pkin(3);
t233 = pkin(4) * t534 - t284 + t529;
t705 = qJ(4) * t534;
t309 = t361 - t705;
t245 = t309 + t505;
t98 = t774 * t233 - t245 * t530;
t889 = t118 + t98;
t886 = t313 + t314;
t468 = -qJ(4) * t530 + t774 * t535;
t467 = -pkin(5) + t468;
t671 = -t467 + t468;
t438 = Ifges(7,6) * t454;
t439 = Ifges(6,6) * t454;
t440 = Ifges(7,5) * t563;
t441 = Ifges(6,5) * t563;
t857 = -t440 - t438 - t441 - t439;
t807 = -t396 / 0.2e1;
t818 = -t550 / 0.2e1;
t773 = m(7) * t181;
t882 = Ifges(5,4) + Ifges(4,5);
t881 = Ifges(6,5) + Ifges(7,5);
t880 = Ifges(4,6) - Ifges(5,6);
t879 = Ifges(6,6) + Ifges(7,6);
t878 = Ifges(6,3) + Ifges(7,3);
t704 = qJ(6) * t397;
t75 = t98 - t704;
t767 = pkin(5) * t534;
t68 = t75 + t767;
t731 = t550 * t68;
t835 = m(7) / 0.2e1;
t872 = t671 * t835;
t369 = Ifges(7,4) * t396;
t190 = t397 * Ifges(7,1) + t534 * Ifges(7,5) - t369;
t371 = Ifges(6,4) * t396;
t192 = t397 * Ifges(6,1) + t534 * Ifges(6,5) - t371;
t871 = t192 + t190;
t209 = mrSges(7,1) * t396 + mrSges(7,2) * t397;
t210 = mrSges(6,1) * t396 + mrSges(6,2) * t397;
t870 = t209 + t210;
t443 = Ifges(7,4) * t454;
t262 = -Ifges(7,2) * t563 + t443;
t445 = Ifges(6,4) * t454;
t264 = -Ifges(6,2) * t563 + t445;
t869 = t264 + t262;
t442 = Ifges(7,4) * t563;
t266 = t454 * Ifges(7,1) - t442;
t444 = Ifges(6,4) * t563;
t268 = t454 * Ifges(6,1) - t444;
t868 = t268 + t266;
t310 = t529 + t670;
t867 = t310 - t670;
t370 = Ifges(7,4) * t397;
t609 = Ifges(7,1) * t396 + t370;
t372 = Ifges(6,4) * t397;
t610 = Ifges(6,1) * t396 + t372;
t261 = -Ifges(7,2) * t454 - t442;
t265 = -Ifges(7,1) * t563 - t443;
t263 = -Ifges(6,2) * t454 - t444;
t267 = -Ifges(6,1) * t563 - t445;
t459 = mrSges(4,2) * t534 - mrSges(4,3) * t683;
t522 = t534 * mrSges(5,3);
t666 = mrSges(5,2) * t683;
t466 = -t522 - t666;
t866 = t466 + t459;
t523 = Ifges(5,5) * t531;
t865 = Ifges(5,1) * t533 + t523;
t524 = Ifges(4,4) * t533;
t477 = t531 * Ifges(4,1) + t524;
t864 = -Ifges(4,2) * t531 + t524;
t581 = Ifges(7,2) * t397 + t369;
t582 = Ifges(6,2) * t397 + t371;
t863 = t582 + t581;
t859 = -(t531 ^ 2 + t533 ^ 2) * t764 / 0.2e1;
t853 = t533 * Ifges(5,3) - t523;
t858 = t865 - t853;
t762 = pkin(8) * t534;
t771 = pkin(2) * t532;
t484 = -t762 + t771;
t687 = t484 * t533;
t362 = pkin(7) * t683 + t687;
t363 = -pkin(7) * t681 + t531 * t484;
t855 = -t362 * t531 + t363 * t533;
t319 = t532 * qJ(4) + t363;
t643 = -pkin(7) * t531 - pkin(3);
t320 = t532 * t643 - t687;
t854 = t319 * t533 + t320 * t531;
t117 = -t284 * t530 + t774 * t285;
t99 = t530 * t233 + t245 * t774;
t645 = t117 / 0.2e1 - t99 / 0.2e1;
t745 = Ifges(5,5) * t533;
t476 = t531 * Ifges(5,1) - t745;
t851 = t864 + t476 + t477;
t461 = -mrSges(4,1) * t534 - mrSges(4,3) * t681;
t462 = mrSges(5,1) * t534 + mrSges(5,2) * t681;
t850 = t462 / 0.2e1 - t461 / 0.2e1;
t723 = t396 * mrSges(6,3);
t312 = -mrSges(6,2) * t534 - t723;
t638 = t774 * t312;
t639 = t774 * t311;
t849 = t638 / 0.2e1 + t639 / 0.2e1;
t618 = t813 + t814;
t826 = mrSges(7,3) / 0.2e1;
t827 = mrSges(6,3) / 0.2e1;
t663 = t826 + t827;
t848 = t663 * t397 - t618;
t469 = qJ(4) * t774 + t530 * t535;
t847 = t755 * t468 + t469 * t756;
t382 = -t534 * Ifges(5,4) + t532 * t865;
t746 = Ifges(4,4) * t531;
t585 = Ifges(4,1) * t533 - t746;
t560 = t585 * t532;
t384 = -t534 * Ifges(4,5) + t560;
t662 = Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1;
t846 = -t662 * t534 + t382 / 0.2e1 + t384 / 0.2e1;
t186 = -t396 * Ifges(7,2) + t534 * Ifges(7,6) + t370;
t188 = -t396 * Ifges(6,2) + t534 * Ifges(6,6) + t372;
t845 = t610 + t609 + t188 + t186;
t561 = t438 / 0.2e1 + t439 / 0.2e1 + t440 / 0.2e1 + t441 / 0.2e1;
t562 = t365 / 0.2e1 + t366 / 0.2e1 + t367 / 0.2e1 + t368 / 0.2e1;
t398 = t534 * t454;
t399 = t563 * t534;
t659 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t661 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t844 = t659 * t398 + t661 * t399;
t842 = 0.2e1 * m(7);
t841 = 2 * qJD(3);
t840 = -m(5) / 0.2e1;
t839 = m(5) / 0.2e1;
t838 = -m(6) / 0.2e1;
t837 = m(6) / 0.2e1;
t836 = -m(7) / 0.2e1;
t834 = -pkin(5) / 0.2e1;
t832 = m(7) * pkin(5);
t831 = mrSges(4,1) / 0.2e1;
t830 = mrSges(7,1) / 0.2e1;
t829 = -mrSges(4,2) / 0.2e1;
t828 = -mrSges(6,3) / 0.2e1;
t825 = -Ifges(5,5) / 0.2e1;
t824 = -t68 / 0.2e1;
t823 = t75 / 0.2e1;
t693 = t396 * qJ(6);
t76 = t99 - t693;
t822 = t76 / 0.2e1;
t88 = t117 - t693;
t821 = -t88 / 0.2e1;
t820 = pkin(7) * mrSges(4,1);
t819 = pkin(7) * mrSges(4,2);
t815 = t312 / 0.2e1;
t812 = -t361 / 0.2e1;
t806 = -t396 / 0.4e1;
t805 = t396 / 0.2e1;
t804 = -t397 / 0.4e1;
t802 = t398 / 0.2e1;
t801 = t399 / 0.2e1;
t648 = t436 / 0.2e1;
t796 = -t563 / 0.2e1;
t795 = -t563 / 0.4e1;
t793 = t454 / 0.2e1;
t792 = t467 / 0.2e1;
t791 = t468 / 0.2e1;
t790 = -t470 / 0.2e1;
t475 = t533 * Ifges(4,2) + t746;
t788 = -t475 / 0.2e1;
t787 = -t530 / 0.2e1;
t785 = -t531 / 0.2e1;
t784 = t531 / 0.2e1;
t782 = -t532 / 0.2e1;
t780 = -t533 / 0.2e1;
t779 = t533 / 0.2e1;
t777 = -t534 / 0.4e1;
t769 = pkin(3) * t531;
t317 = -mrSges(7,1) * t532 - mrSges(7,3) * t399;
t768 = pkin(5) * t317;
t766 = pkin(7) * t532;
t763 = pkin(8) * t533;
t761 = t75 * mrSges(7,2);
t760 = t88 * mrSges(7,1);
t89 = t118 + t704;
t759 = t89 * mrSges(7,2);
t758 = t98 * mrSges(6,2);
t757 = t99 * mrSges(6,1);
t754 = -t68 + t75;
t753 = mrSges(7,3) * t563;
t752 = mrSges(7,3) * t454;
t739 = t117 * mrSges(6,1);
t738 = t118 * mrSges(6,2);
t730 = t550 * t75;
t642 = pkin(7) + t769;
t386 = t532 * t642 - t500;
t725 = t386 * mrSges(5,1);
t724 = t386 * mrSges(5,3);
t720 = t398 * mrSges(7,1);
t719 = t399 * mrSges(7,2);
t718 = t563 * mrSges(6,3);
t716 = t454 * mrSges(6,3);
t715 = t467 * mrSges(7,3);
t714 = t468 * mrSges(6,3);
t234 = (-pkin(9) * t534 - t484) * t533 + (-pkin(4) + t643) * t532;
t249 = pkin(9) * t682 + t319;
t100 = t774 * t234 - t249 * t530;
t101 = t530 * t234 + t774 * t249;
t501 = qJ(4) * t680;
t303 = t534 * t592 + t501;
t182 = -pkin(5) * t398 + t303;
t187 = Ifges(7,4) * t399 + Ifges(7,2) * t398 - t532 * Ifges(7,6);
t189 = Ifges(6,4) * t399 + Ifges(6,2) * t398 - t532 * Ifges(6,6);
t191 = Ifges(7,1) * t399 + Ifges(7,4) * t398 - t532 * Ifges(7,5);
t193 = Ifges(6,1) * t399 + Ifges(6,4) * t398 - t532 * Ifges(6,5);
t211 = t719 - t720;
t212 = -mrSges(6,1) * t398 + mrSges(6,2) * t399;
t315 = mrSges(7,2) * t532 + mrSges(7,3) * t398;
t316 = mrSges(6,2) * t532 + mrSges(6,3) * t398;
t318 = -mrSges(6,1) * t532 - mrSges(6,3) * t399;
t580 = Ifges(5,3) * t531 + t745;
t379 = t532 * Ifges(5,6) + t534 * t580;
t381 = t532 * Ifges(4,6) + t534 * t864;
t383 = t532 * Ifges(5,4) + t534 * t865;
t385 = t532 * Ifges(4,5) + t534 * t585;
t387 = t534 * t642 - t501;
t587 = mrSges(5,1) * t531 - mrSges(5,3) * t533;
t431 = t587 * t532;
t432 = t587 * t534;
t589 = mrSges(4,1) * t531 + mrSges(4,2) * t533;
t433 = t534 * t589;
t460 = -mrSges(4,2) * t532 - mrSges(4,3) * t682;
t463 = mrSges(4,1) * t532 - mrSges(4,3) * t680;
t502 = mrSges(5,2) * t680;
t710 = t532 * mrSges(5,1);
t464 = t502 - t710;
t465 = -mrSges(5,2) * t682 + mrSges(5,3) * t532;
t504 = Ifges(5,5) * t681;
t378 = -Ifges(5,6) * t534 + Ifges(5,3) * t683 + t504;
t708 = t534 * Ifges(4,6);
t380 = t532 * t864 - t708;
t617 = t378 / 0.2e1 - t380 / 0.2e1;
t622 = t192 / 0.2e1 + t190 / 0.2e1;
t623 = t186 / 0.2e1 + t188 / 0.2e1;
t660 = Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t73 = -pkin(5) * t532 - qJ(6) * t399 + t100;
t77 = qJ(6) * t398 + t101;
t5 = t623 * t398 + m(4) * (t361 * t363 - t362 * t670) - t670 * t463 + t622 * t399 + (t193 / 0.2e1 + t191 / 0.2e1) * t397 + (-t187 / 0.2e1 - t189 / 0.2e1) * t396 + m(5) * (t309 * t319 + t310 * t320 + t386 * t387) + m(6) * (t100 * t98 + t101 * t99 + t302 * t303) + m(7) * (t181 * t182 + t68 * t73 + t76 * t77) + t310 * t464 + t309 * t465 + t319 * t466 + t363 * t459 + t361 * t460 + t362 * t461 + t320 * t462 + t387 * t431 + t386 * t432 + t100 * t314 + t76 * t315 + t99 * t316 + t68 * t317 + t98 * t318 + t77 * t311 + t101 * t312 + t73 * t313 + t302 * t212 + t303 * t210 + t182 * t209 + t181 * t211 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t532 + pkin(7) * t433 - t661 * t397 + t659 * t396 + (t383 / 0.2e1 + t385 / 0.2e1 + t662 * t532) * t533 + (t379 / 0.2e1 - t381 / 0.2e1 + t660 * t532) * t531 + (-Ifges(5,2) - Ifges(4,3) - Ifges(3,2) + Ifges(3,1) + (m(4) * pkin(7) + t589) * pkin(7) - t878) * t534) * t532 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t534 + t846 * t533 + (-t534 * t660 + t617) * t531 + t844) * t534;
t713 = t5 * qJD(1);
t337 = (t533 * t535 - t519) * t532;
t201 = -pkin(5) * t397 + t337;
t430 = t579 * t532;
t503 = Ifges(5,6) * t681;
t6 = -t893 - t866 * t670 + m(5) * (-t309 * t670 + t310 * t361 + t386 * t430) + ((t708 / 0.2e1 + t725 + t504 / 0.2e1 - t309 * mrSges(5,2) - t361 * mrSges(4,3) + (t820 - t524 / 0.2e1) * t532 + t617) * t533 + (t724 - t310 * mrSges(5,2) - t670 * mrSges(4,3) + (-t819 + (Ifges(4,4) / 0.2e1 + t825) * t531) * t532 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(5,3) / 0.2e1) * t681 - t846) * t531) * t532 + (-t396 * t68 + t397 * t76) * mrSges(7,3) + (-t396 * t98 + t397 * t99) * mrSges(6,3) - t534 * t503 / 0.2e1 + m(6) * (t117 * t98 + t118 * t99 + t302 * t337) + m(7) * (t181 * t201 + t68 * t88 + t76 * t89) + t430 * t431 + (-t461 + t462) * t361 + t337 * t210 + t117 * t314 + t89 * t311 + t118 * t312 + t88 * t313 + t201 * t209 + t871 * t805 + t863 * t807 + t845 * t397 / 0.2e1;
t707 = t6 * qJD(1);
t654 = m(7) * t754;
t7 = t75 * t311 - t99 * t314 + t98 * t312 + (t654 - t313) * t76 + (-t609 / 0.2e1 - t76 * mrSges(7,3) - t99 * mrSges(6,3) - t610 / 0.2e1 + (t209 + t773) * pkin(5) - t623) * t397 + (t68 * mrSges(7,3) + t98 * mrSges(6,3) + t582 / 0.2e1 + t581 / 0.2e1 - t622) * t396 + t893;
t706 = t7 * qJD(1);
t651 = t774 * t99;
t652 = t774 * t76;
t684 = t530 * t534;
t20 = t886 * t684 + (-m(5) * t386 + m(6) * t302 - t431 + t773 + t870) * t681 + (m(7) * (t530 * t68 - t652) - t639 - t638 + m(6) * (t530 * t98 - t651) - m(5) * t309 - t466) * t534;
t703 = qJD(1) * t20;
t33 = m(7) * (-t396 * t76 - t397 * t68) - t397 * t313 - t396 * t311;
t702 = qJD(1) * t33;
t701 = t550 * t397;
t699 = t307 * t396;
t698 = t564 * t397;
t692 = t467 * t397;
t691 = t467 * t454;
t690 = t469 * t396;
t689 = t469 * t563;
t688 = t470 * t532;
t686 = t530 * t396;
t685 = t530 * t563;
t678 = -t309 + t361;
t669 = t832 / 0.2e1;
t667 = t774 / 0.2e1;
t658 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t657 = t823 + t824;
t656 = t88 / 0.2e1 - t76 / 0.2e1;
t655 = -mrSges(7,1) - t832;
t650 = -t753 / 0.2e1;
t649 = -t722 / 0.2e1;
t644 = t118 / 0.2e1 + t98 / 0.2e1;
t641 = t469 * t774;
t640 = t534 * t774;
t637 = t774 * t397;
t636 = t774 * t454;
t632 = t563 * t807;
t629 = -t684 / 0.2e1;
t628 = t684 / 0.2e1;
t625 = t681 / 0.2e1;
t621 = -t262 / 0.2e1 - t264 / 0.2e1;
t620 = t266 / 0.2e1 + t268 / 0.2e1;
t619 = t815 + t816;
t615 = -t467 / 0.2e1 + t791;
t612 = -mrSges(5,2) * qJ(4) - Ifges(4,6);
t606 = t671 * t550;
t602 = mrSges(5,2) * pkin(3) - t882;
t601 = t787 * t886 + t849;
t600 = t756 * t628 + t755 * t640 / 0.2e1;
t599 = t469 * t640;
t598 = t654 / 0.2e1;
t595 = -t640 / 0.2e1;
t521 = t533 * qJ(4);
t448 = t521 + t653;
t590 = t533 * mrSges(4,1) - t531 * mrSges(4,2);
t283 = -pkin(5) * t454 + t448;
t473 = -t521 + t769;
t536 = t868 * t806 + t871 * t795 - t857 * t777 + (-t267 - t265 + t869) * t804 + (t181 * t283 + t201 * t278 + t550 * t89 + t731 + (-t76 + t88) * t591) * t836 + (t302 * t448 + t337 * t434 + t564 * t889 + (-t117 + t99) * t307) * t838 - t845 * t454 / 0.4e1 + t430 * t789 - t473 * t431 / 0.2e1 - t448 * t210 / 0.2e1 - t337 * t260 / 0.2e1 + t312 * t817 - t283 * t209 / 0.2e1 - t201 * t259 / 0.2e1 + t892;
t549 = t100 * mrSges(6,1) / 0.2e1 - t101 * mrSges(6,2) / 0.2e1 + t73 * t830 - t77 * mrSges(7,2) / 0.2e1;
t544 = -t549 - t844;
t539 = (t316 / 0.2e1 + t315 / 0.2e1) * t469 + t544 + (t100 * t468 + t101 * t469) * t837 + (t467 * t73 + t469 * t77) * t835 + (-pkin(3) * t320 + qJ(4) * t319) * t839 - pkin(3) * t464 / 0.2e1 + qJ(4) * t465 / 0.2e1 + t317 * t792 + t318 * t791 + t362 * t831 + t363 * t829 + t319 * mrSges(5,3) / 0.2e1 - t320 * mrSges(5,1) / 0.2e1;
t568 = Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + t658;
t569 = Ifges(4,1) / 0.4e1 + Ifges(5,1) / 0.4e1 - Ifges(5,3) / 0.4e1 - Ifges(4,2) / 0.4e1;
t572 = (mrSges(4,3) + mrSges(5,2)) * pkin(8) / 0.2e1;
t574 = t386 * t473 + t430 * t470;
t1 = (t724 / 0.2e1 - t384 / 0.4e1 - t382 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(5,4) + 0.3e1 / 0.4e1 * Ifges(4,5)) * t534 - t850 * pkin(8) + (t670 / 0.2e1 - t310 / 0.2e1) * mrSges(5,2)) * t533 + (-t698 / 0.2e1 - t699 / 0.2e1 + t645 * t454 + t644 * t563) * mrSges(6,3) + (t763 * t867 + t574) * t840 + ((mrSges(5,1) * t790 + pkin(2) * t831 + t475 / 0.4e1 + t853 / 0.4e1 - t819 / 0.2e1 - t523 / 0.4e1 + (-t569 + t572) * t533) * t533 + (mrSges(5,3) * t790 + pkin(2) * t829 - t820 / 0.2e1 + t477 / 0.4e1 + t476 / 0.4e1 + t524 / 0.4e1 + (t569 + t572) * t531 + (0.3e1 / 0.4e1 * Ifges(4,4) + t825) * t533) * t531 + t568) * t532 + t539 + (-t701 / 0.2e1 + t591 * t805 + t656 * t454 - (-t89 / 0.2e1 + t824) * t563) * mrSges(7,3) + (-t504 / 0.4e1 - t725 / 0.2e1 + t380 / 0.4e1 - t378 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(5,6) - 0.3e1 / 0.4e1 * Ifges(4,6)) * t534 + (t812 + t309 / 0.2e1) * mrSges(5,2) + (t678 * t840 + t466 / 0.2e1 + t459 / 0.2e1) * pkin(8)) * t531 + t536 - (Ifges(7,4) + Ifges(6,4)) * t632 + (Ifges(7,2) + Ifges(6,2)) * (t454 * t396 / 0.4e1 - t563 * t804);
t10 = t580 * t780 - pkin(2) * t589 + t475 * t785 - t473 * t588 + (m(5) * t473 + t587) * t470 + (t443 / 0.2e1 + t445 / 0.2e1 - t278 * mrSges(7,1) - t434 * mrSges(6,1) - t621) * t454 - (-t434 * mrSges(6,2) - t278 * mrSges(7,2) - (-Ifges(6,4) / 0.2e1 - Ifges(7,4) / 0.2e1) * t563 + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1 - Ifges(7,1) / 0.2e1 - Ifges(6,1) / 0.2e1) * t454 - t620) * t563 + (t585 + t858) * t784 + t851 * t779 + t898 * t448 + t613 * t283;
t578 = -t1 * qJD(1) + t10 * qJD(2);
t538 = (-t386 * t531 + (-t688 - t762) * t533) * t839 + (t531 * t302 + t434 * t681 - t534 * t565) * t837 + (t531 * t181 + t278 * t681 - t534 * t566) * t835 + t431 * t785 + t681 * t789 + t870 * t784 + (t259 + t260) * t625 + t883 * (t454 * t629 - t563 * t595);
t542 = t320 * t839 + (t100 * t774 + t530 * t101) * t837 + (t530 * t77 + t73 * t774) * t835 - t710 / 0.2e1 + (t315 + t316) * t530 / 0.2e1 + (t317 + t318) * t667;
t11 = -t502 + t538 - t542;
t55 = (-m(5) * t470 + t588 + t613 + t898) * t531;
t577 = qJD(1) * t11 + qJD(2) * t55;
t543 = (t397 * t793 - t632) * mrSges(7,3) + (-t396 * t550 - t397 * t591 - t454 * t68 - t563 * t76) * t835 + t311 * t796 + t454 * t814;
t551 = t182 * t835 - t720 / 0.2e1 + t719 / 0.2e1;
t24 = t543 - t551;
t51 = m(7) * (-t454 * t591 - t550 * t563) + (t454 ^ 2 + t563 ^ 2) * mrSges(7,3);
t576 = qJD(1) * t24 + qJD(2) * t51;
t42 = (-t690 / 0.4e1 - t692 / 0.4e1 - t201 / 0.4e1) * t842 - t586;
t53 = 0.2e1 * t648 - t717 + (t283 / 0.4e1 + t689 / 0.4e1 + t691 / 0.4e1) * t842;
t575 = qJD(1) * t42 - qJD(2) * t53;
t123 = t397 * t655 + t364;
t185 = t454 * t655 + t436;
t571 = qJD(1) * t123 + qJD(2) * t185;
t139 = (t625 + t686 / 0.2e1 + t637 / 0.2e1) * m(7);
t200 = (-t685 / 0.2e1 - t636 / 0.2e1 + t785) * m(7);
t570 = qJD(1) * t139 - qJD(2) * t200;
t559 = t714 + t715 + t881;
t557 = t755 * t774;
t556 = qJD(3) * (t469 * t883 - t879);
t16 = -t659 * t454 + (t817 + t307 / 0.2e1) * mrSges(6,2) + (t550 / 0.2e1 + t818) * mrSges(7,1) + (t606 / 0.2e1 + pkin(5) * t818) * m(7) - ((t834 + t615) * mrSges(7,3) + t661) * t563 + t561;
t52 = m(7) * t469 * t671 + t847;
t8 = t618 * t469 + t619 * t468 + (t89 / 0.2e1 + t823) * mrSges(7,2) + t644 * mrSges(6,2) - t656 * mrSges(7,1) - t645 * mrSges(6,1) + (t671 * t822 + t754 * t469 / 0.2e1 + pkin(5) * t821) * m(7) + (-t469 * t663 - t659) * t397 + (t714 / 0.2e1 + (t792 + pkin(5) / 0.2e1) * mrSges(7,3) - t661) * t396 + t562;
t555 = t8 * qJD(1) + t16 * qJD(2) + t52 * qJD(3);
t13 = -t278 * t852 + t434 * t258 + (t267 / 0.2e1 + t265 / 0.2e1 + t613 * pkin(5) + t621) * t454 - (t263 / 0.2e1 + t261 / 0.2e1 + t620) * t563;
t537 = -(-t582 / 0.4e1 - t581 / 0.4e1 + t192 / 0.4e1 + t190 / 0.4e1 + t657 * mrSges(7,3)) * t563 + (-t610 / 0.4e1 - t609 / 0.4e1 - t188 / 0.4e1 - t186 / 0.4e1 + (t209 / 0.2e1 + t773 / 0.2e1) * pkin(5)) * t454 + (-t261 / 0.4e1 - t263 / 0.4e1 - t268 / 0.4e1 - t266 / 0.4e1 + mrSges(6,3) * t817 + t591 * t826) * t396 + (-t264 / 0.4e1 - t262 / 0.4e1 + t267 / 0.4e1 + t265 / 0.4e1 + t564 * t828 + mrSges(7,3) * t818 + (t772 / 0.2e1 + t259 / 0.2e1) * pkin(5)) * t397 - t307 * t815 + t857 * t534 / 0.4e1 + t892;
t4 = (-t731 / 0.2e1 + t730 / 0.2e1 + t73 * t834) * m(7) - t768 / 0.2e1 + t658 * t532 + t537 + t544;
t554 = t4 * qJD(1) + t13 * qJD(2);
t540 = -t522 + (t361 - 0.2e1 * t705) * t839 + (-t599 + t651 + (t468 * t534 - t98) * t530) * t837 + (-t599 + t652 + (t467 * t534 - t68) * t530) * t835 + t601 + t756 * t629 + t755 * t595;
t545 = t883 * (t396 * t667 + t397 * t787);
t541 = m(5) * t812 + (t117 * t774 + t530 * t118) * t838 + (t530 * t89 + t774 * t88) * t836 + t545;
t15 = t540 + t541;
t546 = t565 * t837 + t566 * t835;
t547 = t565 * t838 + t566 * t836;
t30 = t546 + t547;
t84 = m(5) * qJ(4) + mrSges(5,3) + t756 * t530 + m(6) * (-t468 * t530 + t641) + m(7) * (-t467 * t530 + t641) + t557;
t553 = -qJD(1) * t15 - qJD(2) * t30 - qJD(3) * t84;
t18 = ((t767 / 0.2e1 - t657) * m(7) + t848) * t530 + t600 + t883 * t774 * t807 - t849;
t548 = t530 * t872;
t57 = (t669 + t756) * t530 + t548 + t557;
t552 = t18 * qJD(1) - t57 * qJD(3);
t199 = (-t636 - t685) * t835 + m(7) * t784;
t140 = (-t637 - t686) * t835 + m(7) * t625;
t65 = t787 * t832 + t548;
t54 = -t436 / 0.2e1 + t648 + (-t689 - t691 + t283) * t835;
t41 = (-t690 - t692 + t201) * t835;
t27 = (m(5) * pkin(8) + mrSges(5,2)) * t533 + t546 - t547 + t883 * (t530 * t454 - 0.2e1 * t563 * t667);
t23 = t543 + t551;
t19 = t530 * t598 + t628 * t832 + t545 + t600 + t601;
t17 = pkin(5) * t650 + t550 * t669 + t606 * t835 - t615 * t753 + 0.2e1 * t561 - t897;
t14 = t540 - t541 - t666;
t12 = t538 + t542;
t9 = 0.2e1 * t562 + (t830 + t872) * t76 + (t598 - t848) * t469 + (t723 / 0.2e1 + t619) * t468 + t715 * t805 + t739 / 0.2e1 - t738 / 0.2e1 + t757 / 0.2e1 + pkin(5) * t649 + t758 / 0.2e1 + t760 / 0.2e1 - t759 / 0.2e1 + t761 / 0.2e1 + t88 * t669;
t3 = (t730 - t731) * t835 + t768 / 0.2e1 + t73 * t669 + t537 + t549 + t879 * t802 + t881 * t801 + t878 * t782;
t2 = (t821 + t822) * t752 + (t89 + t68) * t650 + (t788 - t853 / 0.4e1 + t858 / 0.4e1) * t681 + ((t361 / 0.2e1 - t309 / 0.2e1) * t531 + t859 + t867 * t779) * mrSges(5,2) + t591 * t649 + (t560 + t384 + t382) * t533 / 0.4e1 + (-Ifges(5,1) * t683 + t378 + t504) * t531 / 0.4e1 + (t531 * t660 + t533 * t662) * t534 + t539 + t701 * t826 + t698 * t827 - t699 * t828 - t531 * t380 / 0.4e1 + t386 * t587 / 0.2e1 + (-t531 * t880 + t533 * t882) * t777 + (-t263 - t261) * t806 + ((t678 * t531 + t533 * t867) * pkin(8) + t574) * t839 + t863 * t795 + (t580 - t477) * t683 / 0.4e1 - t866 * t765 / 0.2e1 + t859 * mrSges(4,3) - t851 * t683 / 0.4e1 - t645 * t716 + t850 * t763 - t536 - t889 * t718 / 0.2e1 + t688 * t789 + t568 * t532 + t589 * t766 / 0.2e1 - t590 * t771 / 0.2e1;
t21 = [qJD(2) * t5 + qJD(3) * t6 + qJD(4) * t20 + qJD(5) * t7 + qJD(6) * t33, t2 * qJD(3) + t12 * qJD(4) + t3 * qJD(5) + t23 * qJD(6) + t713 + (0.2e1 * (-t100 * t307 + t101 * t564 + t303 * t434) * t837 - t307 * t318 + (Ifges(3,5) + (t476 / 0.2e1 + t477 / 0.2e1) * t533 + (-t853 / 0.2e1 + t788) * t531 + (-mrSges(3,1) - t590) * pkin(7)) * t534 + t591 * t317 + 0.2e1 * (t182 * t278 + t550 * t77 + t591 * t73) * t835 + t550 * t315 + t564 * t316 + (t454 * t881 - t563 * t879) * t782 - t387 * t588 + (t189 + t187) * t796 + (t193 + t191) * t793 + (t385 + t383) * t784 + t381 * t779 + t379 * t780 - t77 * t753 + mrSges(3,2) * t766 - t101 * t718 - Ifges(3,6) * t532 + t470 * t432 - pkin(2) * t433 + t434 * t212 + t303 * t260 + ((t460 + t465) * t533 + (-t463 + t464) * t531) * pkin(8) + t278 * t211 + t182 * t259 + (t882 * t531 + t880 * t533) * t532 / 0.2e1 + t868 * t801 + t869 * t802 + 0.2e1 * (pkin(8) * t854 + t387 * t470) * t839 + t854 * mrSges(5,2) + t855 * mrSges(4,3) + m(4) * (-pkin(7) * t770 + pkin(8) * t855) - t100 * t716 - t73 * t752) * qJD(2), t707 + t2 * qJD(2) + t14 * qJD(4) + t9 * qJD(5) + t41 * qJD(6) + ((t117 * t468 + t118 * t469) * t837 + (t467 * t88 + t469 * t89) * t835 + (-pkin(3) * t361 - qJ(4) * t670) * t839) * t841 + t397 * t556 + (t503 + t738 - t739 + t759 - t760 - t559 * t396 + (t531 * t602 + t533 * t612) * t532 + (-mrSges(4,1) - mrSges(5,1)) * t361 - (-mrSges(4,2) + mrSges(5,3)) * t670) * qJD(3), qJD(2) * t12 + qJD(3) * t14 + qJD(5) * t19 + qJD(6) * t140 + t703, t706 + t3 * qJD(2) + t9 * qJD(3) + t19 * qJD(4) + (-t757 - t76 * mrSges(7,1) - t758 - t761 + (-m(7) * t76 + t722) * pkin(5) + t856) * qJD(5), qJD(2) * t23 + qJD(3) * t41 + qJD(4) * t140 + t702; -qJD(3) * t1 + qJD(4) * t11 + qJD(5) * t4 + qJD(6) * t24 - t713, qJD(3) * t10 + qJD(4) * t55 + qJD(5) * t13 + qJD(6) * t51, t27 * qJD(4) + t17 * qJD(5) + t54 * qJD(6) + ((t307 * t469 + t468 * t564) * t837 + (t467 * t550 - t469 * t591) * t835) * t841 + t454 * t556 + t578 + (-t602 * t533 + (Ifges(5,6) + t612) * t531 - t559 * t563 + (-m(5) * t579 - t588 - t590) * pkin(8) + t897) * qJD(3), qJD(3) * t27 + qJD(6) * t199 + t577, t17 * qJD(3) + ((-m(7) * t550 + t753) * pkin(5) + t857 + t897) * qJD(5) + t554, qJD(3) * t54 + qJD(4) * t199 + t576; qJD(2) * t1 + qJD(4) * t15 + qJD(5) * t8 + qJD(6) * t42 - t707, qJD(4) * t30 + qJD(5) * t16 - qJD(6) * t53 - t578, qJD(4) * t84 + qJD(5) * t52, qJD(5) * t65 - t553, t65 * qJD(4) + (-t469 * t832 - t847) * qJD(5) + t555, t575; -qJD(2) * t11 - qJD(3) * t15 - qJD(5) * t18 - qJD(6) * t139 - t703, -qJD(3) * t30 + qJD(6) * t200 - t577, qJD(5) * t57 + t553, 0, -t552 + (-t557 + (-mrSges(6,1) + t655) * t530) * qJD(5), -t570; -qJD(2) * t4 - qJD(3) * t8 + qJD(4) * t18 + qJD(6) * t123 - t706, -qJD(3) * t16 + qJD(6) * t185 - t554, -qJD(4) * t57 - t555, t552, 0, t571; -qJD(2) * t24 - qJD(3) * t42 + qJD(4) * t139 - qJD(5) * t123 - t702, qJD(3) * t53 - qJD(4) * t200 - qJD(5) * t185 - t576, -t575, t570, -t571, 0;];
Cq  = t21;
