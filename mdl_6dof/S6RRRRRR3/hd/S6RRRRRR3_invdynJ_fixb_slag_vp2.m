% Calculate vector of inverse dynamics joint torques for
% S6RRRRRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:39:14
% EndTime: 2019-03-10 03:40:17
% DurationCPUTime: 34.94s
% Computational Cost: add. (37074->1014), mult. (79033->1317), div. (0->0), fcn. (59272->18), ass. (0->481)
t509 = sin(qJ(3));
t515 = cos(qJ(3));
t516 = cos(qJ(2));
t606 = qJD(1) * t516;
t510 = sin(qJ(2));
t607 = qJD(1) * t510;
t406 = -t509 * t607 + t515 * t606;
t430 = t509 * t516 + t510 * t515;
t407 = t430 * qJD(1);
t324 = pkin(3) * t407 - pkin(9) * t406;
t586 = pkin(2) * t607;
t299 = t324 + t586;
t519 = -pkin(8) - pkin(7);
t460 = t519 * t516;
t438 = qJD(1) * t460;
t409 = t509 * t438;
t458 = t519 * t510;
t437 = qJD(1) * t458;
t335 = t437 * t515 + t409;
t508 = sin(qJ(4));
t514 = cos(qJ(4));
t234 = t508 * t299 + t514 * t335;
t603 = qJD(3) * t515;
t584 = pkin(2) * t603;
t836 = t514 * t584 - t234;
t233 = t514 * t299 - t335 * t508;
t835 = -t508 * t584 - t233;
t642 = t406 * t514;
t556 = t407 * pkin(4) - pkin(10) * t642;
t693 = pkin(2) * t509;
t484 = pkin(9) + t693;
t674 = -pkin(10) - t484;
t567 = qJD(4) * t674;
t834 = t514 * t567 - t556 + t835;
t416 = qJD(2) * pkin(2) + t437;
t328 = t416 * t515 + t409;
t236 = t514 * t324 - t328 * t508;
t518 = -pkin(10) - pkin(9);
t576 = qJD(4) * t518;
t833 = t514 * t576 - t236 - t556;
t643 = t406 * t508;
t593 = pkin(10) * t643;
t832 = -t508 * t567 - t593 - t836;
t237 = t508 * t324 + t514 * t328;
t831 = -t508 * t576 + t237 - t593;
t507 = sin(qJ(5));
t513 = cos(qJ(5));
t429 = t507 * t514 + t508 * t513;
t292 = t429 * t406;
t761 = qJD(4) + qJD(5);
t338 = t761 * t429;
t830 = t292 - t338;
t539 = t507 * t508 - t513 * t514;
t293 = t539 * t406;
t337 = t761 * t539;
t829 = t293 - t337;
t504 = qJ(4) + qJ(5);
t498 = qJ(6) + t504;
t482 = sin(t498);
t494 = sin(t504);
t828 = -mrSges(6,2) * t494 - mrSges(7,2) * t482;
t420 = t674 * t508;
t499 = t514 * pkin(10);
t421 = t484 * t514 + t499;
t321 = t507 * t420 + t513 * t421;
t782 = -qJD(5) * t321 + t832 * t507 + t513 * t834;
t599 = qJD(5) * t513;
t600 = qJD(5) * t507;
t781 = t420 * t599 - t421 * t600 + t507 * t834 - t832 * t513;
t457 = t518 * t508;
t459 = pkin(9) * t514 + t499;
t355 = t507 * t457 + t513 * t459;
t780 = -qJD(5) * t355 + t507 * t831 + t513 * t833;
t779 = t457 * t599 - t459 * t600 + t507 * t833 - t513 * t831;
t827 = t830 * pkin(11);
t483 = cos(t498);
t496 = cos(t504);
t825 = mrSges(6,1) * t496 + mrSges(7,1) * t483;
t824 = -t407 * pkin(5) - pkin(11) * t829;
t823 = -mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t822 = t827 + t781;
t821 = t824 + t782;
t820 = t827 + t779;
t819 = t824 + t780;
t619 = t515 * t438;
t334 = t437 * t509 - t619;
t604 = qJD(3) * t509;
t818 = pkin(2) * t604 - t334;
t500 = t516 * pkin(2);
t488 = t500 + pkin(1);
t456 = t488 * qJD(1);
t294 = -pkin(3) * t406 - pkin(9) * t407 - t456;
t329 = t509 * t416 - t619;
t502 = qJD(2) + qJD(3);
t306 = pkin(9) * t502 + t329;
t212 = t514 * t294 - t306 * t508;
t213 = t294 * t508 + t306 * t514;
t512 = cos(qJ(6));
t399 = qJD(4) - t406;
t384 = qJD(5) + t399;
t351 = -t407 * t508 + t502 * t514;
t352 = t407 * t514 + t502 * t508;
t256 = t351 * t507 + t352 * t513;
t811 = pkin(11) * t256;
t168 = -pkin(10) * t352 + t212;
t144 = pkin(4) * t399 + t168;
t169 = pkin(10) * t351 + t213;
t161 = t507 * t169;
t94 = t513 * t144 - t161;
t70 = t94 - t811;
t63 = pkin(5) * t384 + t70;
t506 = sin(qJ(6));
t562 = t513 * t351 - t352 * t507;
t797 = pkin(11) * t562;
t163 = t513 * t169;
t95 = t144 * t507 + t163;
t71 = t95 + t797;
t652 = t506 * t71;
t27 = t512 * t63 - t652;
t650 = t512 * t71;
t28 = t506 * t63 + t650;
t789 = t502 * Ifges(4,6);
t817 = t456 * mrSges(4,1) - t212 * mrSges(5,1) - t94 * mrSges(6,1) - t27 * mrSges(7,1) + t213 * mrSges(5,2) + t95 * mrSges(6,2) + t28 * mrSges(7,2) + t329 * mrSges(4,3) + t789 / 0.2e1;
t790 = t502 * Ifges(4,5);
t816 = t456 * mrSges(4,2) + t328 * mrSges(4,3) - t790 / 0.2e1;
t305 = -pkin(3) * t502 - t328;
t244 = -pkin(4) * t351 + t305;
t147 = -pkin(5) * t562 + t244;
t154 = t256 * t512 + t506 * t562;
t596 = qJD(1) * qJD(2);
t444 = qJDD(1) * t516 - t510 * t596;
t445 = qJDD(1) * t510 + t516 * t596;
t283 = -qJD(3) * t407 + t444 * t515 - t509 * t445;
t279 = qJDD(4) - t283;
t274 = qJDD(5) + t279;
t269 = qJDD(6) + t274;
t428 = t509 * t510 - t515 * t516;
t530 = t428 * qJD(3);
t282 = -qJD(1) * t530 + t444 * t509 + t445 * t515;
t501 = qJDD(2) + qJDD(3);
t216 = qJD(4) * t351 + t282 * t514 + t501 * t508;
t217 = -qJD(4) * t352 - t282 * t508 + t501 * t514;
t102 = qJD(5) * t562 + t216 * t513 + t217 * t507;
t103 = -qJD(5) * t256 - t216 * t507 + t217 * t513;
t563 = -t256 * t506 + t512 * t562;
t33 = qJD(6) * t563 + t102 * t512 + t103 * t506;
t34 = -qJD(6) * t154 - t102 * t506 + t103 * t512;
t590 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t269;
t648 = qJDD(1) * pkin(1);
t394 = -pkin(2) * t444 - t648;
t170 = -pkin(3) * t283 - pkin(9) * t282 + t394;
t424 = t445 * pkin(7);
t343 = qJDD(2) * pkin(2) - pkin(8) * t445 - t424;
t423 = t444 * pkin(7);
t350 = pkin(8) * t444 + t423;
t191 = t509 * t343 + t515 * t350 + t416 * t603 + t438 * t604;
t182 = pkin(9) * t501 + t191;
t73 = -qJD(4) * t213 + t514 * t170 - t182 * t508;
t51 = pkin(4) * t279 - pkin(10) * t216 + t73;
t601 = qJD(4) * t514;
t602 = qJD(4) * t508;
t72 = t508 * t170 + t514 * t182 + t294 * t601 - t306 * t602;
t54 = pkin(10) * t217 + t72;
t17 = -qJD(5) * t95 - t507 * t54 + t513 * t51;
t7 = pkin(5) * t274 - pkin(11) * t102 + t17;
t16 = t144 * t599 - t169 * t600 + t507 * t51 + t513 * t54;
t8 = pkin(11) * t103 + t16;
t3 = qJD(6) * t27 + t506 * t7 + t512 * t8;
t4 = -qJD(6) * t28 - t506 * t8 + t512 * t7;
t756 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t535 = t590 + t756;
t579 = Ifges(6,5) * t102 + Ifges(6,6) * t103 + Ifges(6,3) * t274;
t694 = mrSges(7,3) * t28;
t695 = mrSges(7,3) * t27;
t379 = qJD(6) + t384;
t708 = -t379 / 0.2e1;
t723 = -t154 / 0.2e1;
t725 = -t563 / 0.2e1;
t146 = Ifges(7,4) * t563;
t86 = Ifges(7,1) * t154 + Ifges(7,5) * t379 + t146;
t734 = -t86 / 0.2e1;
t656 = Ifges(7,4) * t154;
t85 = Ifges(7,2) * t563 + Ifges(7,6) * t379 + t656;
t736 = -t85 / 0.2e1;
t753 = t17 * mrSges(6,1) - t16 * mrSges(6,2);
t815 = t535 + t579 + t753 + (-mrSges(7,2) * t147 + Ifges(7,1) * t723 + Ifges(7,4) * t725 + Ifges(7,5) * t708 + t695 + t734) * t563 - (mrSges(7,1) * t147 + Ifges(7,4) * t723 + Ifges(7,2) * t725 + Ifges(7,6) * t708 - t694 + t736) * t154;
t505 = qJ(2) + qJ(3);
t495 = sin(t505);
t669 = mrSges(5,2) * t508;
t814 = (-t669 + t828) * t495;
t813 = -m(5) * pkin(9) + t823;
t685 = pkin(5) * t256;
t491 = pkin(4) * t602;
t808 = -pkin(5) * t830 + t491;
t776 = mrSges(4,1) * t502 + mrSges(5,1) * t351 - mrSges(5,2) * t352 - mrSges(4,3) * t407;
t380 = pkin(4) * t643;
t807 = -t380 + t818;
t339 = -qJD(2) * t428 - t530;
t575 = t430 * t601;
t533 = t339 * t508 + t575;
t806 = -t508 * t73 + t514 * t72;
t673 = mrSges(5,1) * t514;
t805 = t669 - t673;
t791 = t399 * Ifges(5,3);
t792 = t351 * Ifges(5,6);
t804 = t352 * Ifges(5,5) + t256 * Ifges(6,5) + t154 * Ifges(7,5) + Ifges(6,6) * t562 + Ifges(7,6) * t563 + t384 * Ifges(6,3) + t379 * Ifges(7,3) + t791 + t792;
t803 = (t673 + t825) * t495;
t497 = cos(t505);
t801 = -t497 * mrSges(4,1) + (mrSges(4,2) + t823) * t495;
t770 = t497 * pkin(3) + t495 * pkin(9);
t686 = pkin(4) * t514;
t486 = pkin(3) + t686;
t772 = t497 * t486 - t495 * t518;
t449 = pkin(5) * t496 + t686;
t436 = pkin(3) + t449;
t503 = -pkin(11) + t518;
t775 = t497 * t436 - t495 * t503;
t800 = -m(5) * t770 - m(6) * t772 - m(7) * t775;
t744 = m(6) * pkin(4);
t740 = t33 / 0.2e1;
t739 = t34 / 0.2e1;
t732 = t102 / 0.2e1;
t731 = t103 / 0.2e1;
t721 = t216 / 0.2e1;
t720 = t217 / 0.2e1;
t715 = t269 / 0.2e1;
t714 = t274 / 0.2e1;
t713 = t279 / 0.2e1;
t798 = t444 / 0.2e1;
t698 = t516 / 0.2e1;
t320 = t513 * t420 - t421 * t507;
t679 = pkin(11) * t429;
t280 = t320 - t679;
t419 = t539 * pkin(11);
t281 = -t419 + t321;
t180 = t280 * t512 - t281 * t506;
t796 = qJD(6) * t180 + t821 * t506 + t512 * t822;
t181 = t280 * t506 + t281 * t512;
t795 = -qJD(6) * t181 - t506 * t822 + t821 * t512;
t353 = t513 * t457 - t459 * t507;
t303 = t353 - t679;
t304 = -t419 + t355;
t219 = t303 * t512 - t304 * t506;
t794 = qJD(6) * t219 + t506 * t819 + t512 * t820;
t220 = t303 * t506 + t304 * t512;
t793 = -qJD(6) * t220 - t506 * t820 + t512 * t819;
t788 = t516 * Ifges(3,2);
t687 = pkin(4) * t513;
t485 = pkin(5) + t687;
t597 = qJD(6) * t512;
t598 = qJD(6) * t506;
t631 = t507 * t512;
t105 = -t168 * t507 - t163;
t74 = t105 - t797;
t106 = t513 * t168 - t161;
t75 = t106 - t811;
t787 = t506 * t75 - t512 * t74 - t485 * t598 + (-t507 * t597 + (-t506 * t513 - t631) * qJD(5)) * pkin(4);
t632 = t506 * t507;
t786 = -t506 * t74 - t512 * t75 + t485 * t597 + (-t507 * t598 + (t512 * t513 - t632) * qJD(5)) * pkin(4);
t785 = t744 + mrSges(5,1);
t551 = mrSges(5,1) * t508 + mrSges(5,2) * t514;
t784 = t305 * t551;
t312 = t539 * t430;
t267 = t380 + t329;
t778 = -t267 + t808;
t777 = t807 + t808;
t327 = pkin(3) * t428 - pkin(9) * t430 - t488;
t356 = t458 * t509 - t460 * t515;
t245 = t514 * t327 - t356 * t508;
t638 = t430 * t514;
t198 = pkin(4) * t428 - pkin(10) * t638 + t245;
t345 = t514 * t356;
t246 = t508 * t327 + t345;
t639 = t430 * t508;
t224 = -pkin(10) * t639 + t246;
t126 = t507 * t198 + t513 * t224;
t774 = t515 * t458 + t460 * t509;
t773 = t491 + t807;
t769 = t491 - t267;
t549 = -mrSges(7,1) * t482 - mrSges(7,2) * t483;
t767 = mrSges(6,1) * t494 + mrSges(6,2) * t496 - t549;
t681 = pkin(7) * t516;
t682 = pkin(7) * t510;
t766 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t607) * t681 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t606) * t682;
t511 = sin(qJ(1));
t517 = cos(qJ(1));
t633 = t497 * t517;
t376 = -t494 * t633 + t496 * t511;
t377 = t494 * t511 + t496 * t633;
t366 = -t482 * t633 + t483 * t511;
t367 = t482 * t511 + t483 * t633;
t617 = t366 * mrSges(7,1) - t367 * mrSges(7,2);
t765 = -t376 * mrSges(6,1) + t377 * mrSges(6,2) - t617;
t634 = t497 * t511;
t374 = t494 * t634 + t496 * t517;
t375 = t494 * t517 - t496 * t634;
t364 = t482 * t634 + t483 * t517;
t365 = t482 * t517 - t483 * t634;
t618 = -t364 * mrSges(7,1) + t365 * mrSges(7,2);
t764 = t374 * mrSges(6,1) - t375 * mrSges(6,2) - t618;
t763 = t423 * t516 + t424 * t510;
t762 = g(1) * t517 + g(2) * t511;
t760 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t759 = m(7) + m(5) + m(4) + m(6);
t540 = -t486 * t495 - t497 * t518;
t542 = -t436 * t495 - t497 * t503;
t690 = pkin(3) * t495;
t692 = pkin(2) * t510;
t758 = -m(7) * (t542 - t692) - m(6) * (t540 - t692) - m(5) * (-t690 - t692) + t803;
t757 = t801 + (t805 - t825 - t828) * t497;
t683 = pkin(5) * t494;
t688 = pkin(4) * t508;
t448 = t683 + t688;
t755 = -m(6) * t688 - m(7) * t448;
t754 = t73 * mrSges(5,1) - t72 * mrSges(5,2);
t752 = t511 * t814 + t634 * t813;
t751 = t517 * t814 + t633 * t813;
t750 = m(5) * t690 - m(6) * t540 - m(7) * t542 + t803;
t455 = -mrSges(3,1) * t516 + mrSges(3,2) * t510;
t749 = m(3) * pkin(1) + mrSges(2,1) - t455 - t801;
t140 = mrSges(5,1) * t279 - mrSges(5,3) * t216;
t141 = -mrSges(5,2) * t279 + mrSges(5,3) * t217;
t666 = mrSges(5,3) * t351;
t285 = -mrSges(5,2) * t399 + t666;
t665 = mrSges(5,3) * t352;
t286 = mrSges(5,1) * t399 - t665;
t748 = m(5) * ((-t212 * t514 - t213 * t508) * qJD(4) + t806) - t286 * t601 - t285 * t602 + t514 * t141 - t508 * t140;
t743 = m(7) * pkin(5);
t742 = Ifges(7,4) * t740 + Ifges(7,2) * t739 + Ifges(7,6) * t715;
t741 = Ifges(7,1) * t740 + Ifges(7,4) * t739 + Ifges(7,5) * t715;
t738 = Ifges(6,4) * t732 + Ifges(6,2) * t731 + Ifges(6,6) * t714;
t737 = Ifges(6,1) * t732 + Ifges(6,4) * t731 + Ifges(6,5) * t714;
t735 = t85 / 0.2e1;
t733 = t86 / 0.2e1;
t730 = Ifges(5,1) * t721 + Ifges(5,4) * t720 + Ifges(5,5) * t713;
t657 = Ifges(6,4) * t256;
t138 = Ifges(6,2) * t562 + Ifges(6,6) * t384 + t657;
t729 = -t138 / 0.2e1;
t728 = t138 / 0.2e1;
t251 = Ifges(6,4) * t562;
t139 = Ifges(6,1) * t256 + Ifges(6,5) * t384 + t251;
t727 = -t139 / 0.2e1;
t726 = t139 / 0.2e1;
t724 = t563 / 0.2e1;
t722 = t154 / 0.2e1;
t719 = -t562 / 0.2e1;
t718 = t562 / 0.2e1;
t717 = -t256 / 0.2e1;
t716 = t256 / 0.2e1;
t711 = -t351 / 0.2e1;
t710 = -t352 / 0.2e1;
t709 = t352 / 0.2e1;
t707 = t379 / 0.2e1;
t706 = -t384 / 0.2e1;
t705 = t384 / 0.2e1;
t704 = -t399 / 0.2e1;
t702 = t406 / 0.2e1;
t700 = t407 / 0.2e1;
t697 = mrSges(6,3) * t94;
t696 = mrSges(6,3) * t95;
t691 = pkin(2) * t515;
t689 = pkin(4) * t352;
t676 = g(3) * t495;
t664 = mrSges(6,3) * t562;
t663 = mrSges(6,3) * t256;
t662 = Ifges(3,4) * t510;
t661 = Ifges(3,4) * t516;
t660 = Ifges(5,4) * t352;
t659 = Ifges(5,4) * t508;
t658 = Ifges(5,4) * t514;
t653 = t407 * Ifges(4,4);
t645 = t339 * t514;
t239 = Ifges(5,2) * t351 + Ifges(5,6) * t399 + t660;
t629 = t508 * t239;
t628 = t508 * t511;
t627 = t508 * t517;
t624 = t511 * t448;
t623 = t511 * t514;
t349 = Ifges(5,4) * t351;
t240 = Ifges(5,1) * t352 + Ifges(5,5) * t399 + t349;
t621 = t514 * t240;
t620 = t514 * t517;
t605 = qJD(2) * t510;
t585 = pkin(2) * t605;
t578 = Ifges(5,5) * t216 + Ifges(5,6) * t217 + Ifges(5,3) * t279;
t577 = qJD(2) * t519;
t572 = t621 / 0.2e1;
t569 = -t602 / 0.2e1;
t566 = t596 / 0.2e1;
t125 = t513 * t198 - t224 * t507;
t340 = t502 * t430;
t243 = pkin(3) * t340 - pkin(9) * t339 + t585;
t442 = t510 * t577;
t443 = t516 * t577;
t259 = qJD(3) * t774 + t442 * t515 + t443 * t509;
t564 = t514 * t243 - t259 * t508;
t295 = pkin(4) * t639 - t774;
t554 = mrSges(3,1) * t510 + mrSges(3,2) * t516;
t552 = mrSges(4,1) * t495 + mrSges(4,2) * t497;
t548 = Ifges(5,1) * t514 - t659;
t547 = t662 + t788;
t546 = -Ifges(5,2) * t508 + t658;
t545 = Ifges(3,5) * t516 - Ifges(3,6) * t510;
t544 = Ifges(5,5) * t514 - Ifges(5,6) * t508;
t104 = pkin(5) * t428 + pkin(11) * t312 + t125;
t311 = t429 * t430;
t113 = -pkin(11) * t311 + t126;
t47 = t104 * t512 - t113 * t506;
t48 = t104 * t506 + t113 * t512;
t225 = -t311 * t512 + t312 * t506;
t226 = -t311 * t506 - t312 * t512;
t331 = -t429 * t506 - t512 * t539;
t332 = t429 * t512 - t506 * t539;
t378 = pkin(5) * t539 - t486;
t192 = t343 * t515 - t509 * t350 - t416 * t604 + t438 * t603;
t534 = pkin(1) * t554;
t392 = -t497 * t627 + t623;
t390 = t497 * t628 + t620;
t532 = t430 * t602 - t645;
t531 = t510 * (Ifges(3,1) * t516 - t662);
t82 = -pkin(10) * t645 + pkin(4) * t340 + (-t345 + (pkin(10) * t430 - t327) * t508) * qJD(4) + t564;
t119 = t508 * t243 + t514 * t259 + t327 * t601 - t356 * t602;
t99 = -pkin(10) * t533 + t119;
t21 = t198 * t599 - t224 * t600 + t507 * t82 + t513 * t99;
t183 = -pkin(3) * t501 - t192;
t121 = -pkin(4) * t217 + t183;
t22 = -qJD(5) * t126 - t507 * t99 + t513 * t82;
t260 = qJD(3) * t356 + t442 * t509 - t515 * t443;
t178 = pkin(4) * t533 + t260;
t109 = Ifges(5,4) * t216 + Ifges(5,2) * t217 + Ifges(5,6) * t279;
t185 = qJD(6) * t331 - t337 * t512 - t338 * t506;
t186 = -qJD(6) * t332 + t337 * t506 - t338 * t512;
t205 = -t292 * t512 + t293 * t506;
t206 = -t292 * t506 - t293 * t512;
t300 = t406 * Ifges(4,2) + t653 + t789;
t396 = Ifges(4,4) * t406;
t301 = t407 * Ifges(4,1) + t396 + t790;
t52 = -pkin(5) * t103 + t121;
t520 = (-t205 * t28 + t206 * t27 + t3 * t331 - t332 * t4) * mrSges(7,3) + ((t185 - t206) * mrSges(7,2) + (-t186 + t205) * mrSges(7,1)) * t147 - (Ifges(6,4) * t716 + Ifges(6,2) * t718 + Ifges(6,6) * t705 + t696 + t728) * t338 + t183 * t805 + ((-t602 + t643) * t213 + (-t601 + t642) * t212 + t806) * mrSges(5,3) - (Ifges(4,1) * t406 - t653 + t804) * t407 / 0.2e1 - (Ifges(6,1) * t716 + Ifges(6,4) * t718 + Ifges(6,5) * t705 - t697 + t726) * t337 + (Ifges(7,1) * t722 + Ifges(7,4) * t724 + Ifges(7,5) * t707 - t695 + t733) * t185 + (Ifges(7,4) * t722 + Ifges(7,2) * t724 + Ifges(7,6) * t707 + t694 + t735) * t186 + (Ifges(7,1) * t332 + Ifges(7,4) * t331) * t740 + t332 * t741 + t331 * t742 + (t784 + t572) * qJD(4) + (-mrSges(6,1) * t830 + mrSges(6,2) * t829) * t244 + (-Ifges(6,4) * t293 - Ifges(6,2) * t292) * t719 - t293 * t727 - t292 * t729 + t508 * t730 + t206 * t734 + t205 * t736 + t429 * t737 + (Ifges(7,4) * t332 + Ifges(7,2) * t331) * t739 + (Ifges(7,4) * t206 + Ifges(7,2) * t205) * t725 + t239 * t569 + (-Ifges(6,5) * t293 - Ifges(6,6) * t292) * t706 + (Ifges(5,5) * t508 + Ifges(5,6) * t514) * t713 + (Ifges(7,5) * t332 + Ifges(7,6) * t331) * t715 + (Ifges(5,2) * t514 + t659) * t720 + (Ifges(5,1) * t508 + t658) * t721 + t300 * t700 + t629 * t702 + (Ifges(7,5) * t206 + Ifges(7,6) * t205) * t708 + (Ifges(7,1) * t206 + Ifges(7,4) * t205) * t723 + (-Ifges(6,1) * t293 - Ifges(6,4) * t292) * t717 + t514 * t109 / 0.2e1 + Ifges(4,3) * t501 + (t544 * t704 + t546 * t711 + t548 * t710 - t784 + t816) * t406 + (Ifges(5,5) * t710 + Ifges(6,5) * t717 + Ifges(7,5) * t723 + Ifges(5,6) * t711 + Ifges(6,6) * t719 + Ifges(7,6) * t725 + Ifges(5,3) * t704 + Ifges(6,3) * t706 + Ifges(7,3) * t708 + t817) * t407 - t191 * mrSges(4,2) + t192 * mrSges(4,1) + Ifges(4,5) * t282 + Ifges(4,6) * t283 + t52 * (-mrSges(7,1) * t331 + mrSges(7,2) * t332) + (t351 * t546 + t352 * t548 + t399 * t544) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t407 + t301 + t396 + t621) * t406 / 0.2e1 + (Ifges(6,4) * t429 - Ifges(6,2) * t539) * t731 + (-t16 * t539 - t17 * t429 + t292 * t95 - t293 * t94) * mrSges(6,3) + (Ifges(6,1) * t429 - Ifges(6,4) * t539) * t732 + (Ifges(6,5) * t429 - Ifges(6,6) * t539) * t714 + t121 * (mrSges(6,1) * t539 + mrSges(6,2) * t429) - t539 * t738;
t490 = Ifges(3,4) * t606;
t487 = -pkin(3) - t691;
t452 = -t486 - t691;
t405 = Ifges(3,1) * t607 + Ifges(3,5) * qJD(2) + t490;
t404 = Ifges(3,6) * qJD(2) + qJD(1) * t547;
t402 = pkin(4) * t631 + t485 * t506;
t401 = -pkin(4) * t632 + t485 * t512;
t393 = t497 * t620 + t628;
t391 = -t497 * t623 + t627;
t368 = -mrSges(4,2) * t502 + mrSges(4,3) * t406;
t363 = t378 - t691;
t323 = -mrSges(4,1) * t406 + mrSges(4,2) * t407;
t265 = -mrSges(4,2) * t501 + mrSges(4,3) * t283;
t264 = mrSges(4,1) * t501 - mrSges(4,3) * t282;
t230 = pkin(5) * t311 + t295;
t229 = mrSges(6,1) * t384 - t663;
t228 = -mrSges(6,2) * t384 + t664;
t211 = t689 + t685;
t160 = -mrSges(6,1) * t562 + mrSges(6,2) * t256;
t136 = mrSges(7,1) * t379 - mrSges(7,3) * t154;
t135 = -mrSges(7,2) * t379 + mrSges(7,3) * t563;
t134 = t312 * t761 - t429 * t339;
t133 = -t338 * t430 - t339 * t539;
t127 = -mrSges(5,1) * t217 + mrSges(5,2) * t216;
t120 = -qJD(4) * t246 + t564;
t98 = -mrSges(7,1) * t563 + mrSges(7,2) * t154;
t93 = -pkin(5) * t134 + t178;
t77 = -mrSges(6,2) * t274 + mrSges(6,3) * t103;
t76 = mrSges(6,1) * t274 - mrSges(6,3) * t102;
t56 = -qJD(6) * t226 - t133 * t506 + t134 * t512;
t55 = qJD(6) * t225 + t133 * t512 + t134 * t506;
t46 = -mrSges(6,1) * t103 + mrSges(6,2) * t102;
t30 = t512 * t70 - t652;
t29 = -t506 * t70 - t650;
t26 = -mrSges(7,2) * t269 + mrSges(7,3) * t34;
t25 = mrSges(7,1) * t269 - mrSges(7,3) * t33;
t19 = pkin(11) * t134 + t21;
t18 = pkin(5) * t340 - pkin(11) * t133 + t22;
t13 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t6 = -qJD(6) * t48 + t18 * t512 - t19 * t506;
t5 = qJD(6) * t47 + t18 * t506 + t19 * t512;
t1 = [(-t133 * t94 + t134 * t95 - t16 * t311 + t17 * t312) * mrSges(6,3) + (-Ifges(6,5) * t312 - Ifges(6,6) * t311) * t714 + (-Ifges(6,4) * t312 - Ifges(6,2) * t311) * t731 + (-Ifges(6,1) * t312 - Ifges(6,4) * t311) * t732 + t121 * (mrSges(6,1) * t311 - mrSges(6,2) * t312) + (t590 + t579 + t578) * t428 / 0.2e1 + t445 * t661 / 0.2e1 + (t225 * t3 - t226 * t4 - t27 * t55 + t28 * t56) * mrSges(7,3) + (Ifges(7,4) * t226 + Ifges(7,2) * t225) * t739 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t724 + (t212 * t532 - t213 * t533 - t638 * t73 - t639 * t72) * mrSges(5,3) + (-t628 * t744 - m(7) * t624 - t393 * mrSges(5,1) - t377 * mrSges(6,1) - t367 * mrSges(7,1) - t392 * mrSges(5,2) - t376 * mrSges(6,2) - t366 * mrSges(7,2) - t759 * (t517 * t488 - t511 * t519) + t760 * t511 + (-t749 + t800) * t517) * g(2) + (t516 * t661 + t531) * t566 - (-m(4) * t192 + m(5) * t183 + t127 - t264) * t774 + (Ifges(7,1) * t226 + Ifges(7,4) * t225) * t740 + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t722 + t399 * (-Ifges(5,5) * t532 - Ifges(5,6) * t533) / 0.2e1 + (Ifges(6,5) * t133 + Ifges(6,6) * t134) * t705 - t109 * t639 / 0.2e1 + t226 * t741 + t225 * t742 + m(4) * (t191 * t356 + t259 * t329 - t394 * t488 - t456 * t585) + (-Ifges(5,1) * t532 - Ifges(5,4) * t533) * t709 - t455 * t648 + t133 * t726 + t134 * t728 + t638 * t730 + t55 * t733 + t56 * t735 - t312 * t737 - t311 * t738 + (Ifges(6,4) * t133 + Ifges(6,2) * t134) * t718 - t534 * t596 - t404 * t605 / 0.2e1 - t239 * t575 / 0.2e1 + m(5) * (t119 * t213 + t120 * t212 + t245 * t73 + t246 * t72) + (Ifges(7,5) * t226 + Ifges(7,6) * t225) * t715 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t707 + (Ifges(6,1) * t133 + Ifges(6,4) * t134) * t716 + t147 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + (Ifges(3,4) * t445 + Ifges(3,2) * t444) * t698 + t5 * t135 + t6 * t136 + t125 * t76 + t126 * t77 + t93 * t98 + Ifges(2,3) * qJDD(1) + t48 * t26 + t47 * t25 + m(7) * (t147 * t93 + t230 * t52 + t27 * t6 + t28 * t5 + t3 * t48 + t4 * t47) + m(6) * (t121 * t295 + t125 * t17 + t126 * t16 + t178 * t244 + t21 * t95 + t22 * t94) + t305 * (mrSges(5,1) * t533 - mrSges(5,2) * t532) - t488 * (-mrSges(4,1) * t283 + mrSges(4,2) * t282) - pkin(1) * (-mrSges(3,1) * t444 + mrSges(3,2) * t445) + (-mrSges(3,1) * t682 - mrSges(3,2) * t681 + 0.2e1 * Ifges(3,6) * t698) * qJDD(2) + (Ifges(3,1) * t445 + Ifges(3,4) * t798 + Ifges(3,5) * qJDD(2) - t566 * t788) * t510 + t259 * t368 + t356 * t265 + t547 * t798 + t323 * t585 + (Ifges(4,1) * t700 + Ifges(4,4) * t702 - t629 / 0.2e1 + t301 / 0.2e1 + t572 - t816) * t339 + (t791 / 0.2e1 + Ifges(7,6) * t724 - Ifges(4,4) * t700 - Ifges(4,2) * t702 + Ifges(6,3) * t705 + Ifges(7,3) * t707 + Ifges(5,5) * t709 + t792 / 0.2e1 - t300 / 0.2e1 + Ifges(6,5) * t716 + Ifges(6,6) * t718 + Ifges(7,5) * t722 + t804 / 0.2e1 - t817) * t340 + t178 * t160 + (t394 * mrSges(4,1) - t191 * mrSges(4,3) - Ifges(4,4) * t282 + Ifges(5,5) * t721 + Ifges(6,5) * t732 + Ifges(7,5) * t740 - Ifges(4,2) * t283 - Ifges(4,6) * t501 + Ifges(5,6) * t720 + Ifges(6,6) * t731 + Ifges(7,6) * t739 + Ifges(5,3) * t713 + Ifges(6,3) * t714 + Ifges(7,3) * t715 + t753 + t754 + t756) * t428 + t351 * (-Ifges(5,4) * t532 - Ifges(5,2) * t533) / 0.2e1 + t52 * (-mrSges(7,1) * t225 + mrSges(7,2) * t226) + t21 * t228 + t22 * t229 + t230 * t13 + t244 * (-mrSges(6,1) * t134 + mrSges(6,2) * t133) + t245 * t140 + t246 * t141 + (t394 * mrSges(4,2) - t192 * mrSges(4,3) + Ifges(4,1) * t282 + Ifges(4,4) * t283 + Ifges(4,5) * t501 + t183 * t551 + t240 * t569 + t544 * t713 + t546 * t720 + t548 * t721) * t430 + (t444 * t681 + t445 * t682 + t763) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t763) + (t405 * t698 + t545 * qJD(2) / 0.2e1 - t766) * qJD(2) + (-m(4) * t328 + m(5) * t305 - t776) * t260 + t119 * t285 + t120 * t286 + t295 * t46 + (-t391 * mrSges(5,1) - t375 * mrSges(6,1) - t365 * mrSges(7,1) - t390 * mrSges(5,2) - t374 * mrSges(6,2) - t364 * mrSges(7,2) + (t519 * t759 + t755 + t760) * t517 + (-m(7) * (-t488 - t775) - m(5) * (-t488 - t770) + m(4) * t488 - m(6) * (-t488 - t772) + t749) * t511) * g(1); -(-Ifges(3,2) * t607 + t405 + t490) * t606 / 0.2e1 + (-t335 + t584) * t368 + t748 * t484 + (m(4) * t692 + t552 + t554) * t762 + (t328 * t334 - t329 * t335 + t456 * t586 + (t191 * t509 + t192 * t515 + (-t328 * t509 + t329 * t515) * qJD(3)) * pkin(2)) * m(4) - t323 * t586 + (t183 * t487 + (t305 * t509 + (-t212 * t508 + t213 * t514) * t515) * qJD(3) * pkin(2) - t212 * t233 - t213 * t234 - t305 * t334) * m(5) + t835 * t286 + t836 * t285 - t545 * t596 / 0.2e1 + (t766 + (-t531 / 0.2e1 + t534) * qJD(1)) * qJD(1) + t520 + t404 * t607 / 0.2e1 + t264 * t691 + t265 * t693 + Ifges(3,3) * qJDD(2) + t487 * t127 + Ifges(3,5) * t445 + t452 * t46 + Ifges(3,6) * t444 - t423 * mrSges(3,2) - t424 * mrSges(3,1) + t795 * t136 + t796 * t135 + (t147 * t777 + t180 * t4 + t181 * t3 + t27 * t795 + t28 * t796 + t363 * t52) * m(7) + t363 * t13 - t776 * t818 + t180 * t25 + t181 * t26 + (t517 * t758 + t751) * g(1) + (t511 * t758 + t752) * g(2) + t773 * t160 + t777 * t98 + t781 * t228 + t782 * t229 + (t121 * t452 + t16 * t321 + t17 * t320 + t773 * t244 + t781 * t95 + t782 * t94) * m(6) + t320 * t76 + t321 * t77 + (-m(4) * t500 - m(6) * (t500 + t772) - m(7) * (t500 + t775) + t455 - m(5) * (t500 + t770) + t757) * g(3); (-pkin(3) * t183 - t212 * t236 - t213 * t237 - t305 * t329) * m(5) + (t517 * t750 + t751) * g(1) + (t757 + t800) * g(3) + t748 * pkin(9) + t520 - pkin(3) * t127 - t486 * t46 + t793 * t136 + t794 * t135 + (t147 * t778 + t219 * t4 + t220 * t3 + t27 * t793 + t28 * t794 + t378 * t52) * m(7) + t378 * t13 - t328 * t368 + t355 * t77 + t353 * t76 + (t511 * t750 + t752) * g(2) + t219 * t25 + t220 * t26 + t762 * t552 + t769 * t160 + t776 * t329 + t778 * t98 + t779 * t228 + t780 * t229 + (-t121 * t486 + t16 * t355 + t17 * t353 + t244 * t769 + t779 * t95 + t780 * t94) * m(6) - t237 * t285 - t236 * t286; (-mrSges(6,2) * t244 + Ifges(6,1) * t717 + Ifges(6,4) * t719 + Ifges(6,5) * t706 + t697 + t727) * t562 + (t228 * t599 - t229 * t600 + t507 * t77) * pkin(4) + (t666 - t285) * t212 + t754 + (t16 * t507 + t17 * t513 + (-t507 * t94 + t513 * t95) * qJD(5)) * t744 + t815 + (t665 + t286) * t213 + (-Ifges(5,2) * t352 + t240 + t349) * t711 + t578 + (Ifges(5,5) * t351 - Ifges(5,6) * t352) * t704 + t239 * t709 + (Ifges(5,1) * t351 - t660) * t710 + t76 * t687 + t401 * t25 + t402 * t26 - t305 * (mrSges(5,1) * t352 + mrSges(5,2) * t351) - t211 * t98 - t106 * t228 - t105 * t229 + (t551 - t755 + t767) * t676 + (-m(7) * (-t449 * t517 - t497 * t624) - mrSges(5,2) * t391 + t785 * t390 + t764) * g(2) + (-m(7) * (-t448 * t633 + t449 * t511) + mrSges(5,2) * t393 - t785 * t392 + t765) * g(1) + t786 * t135 + t787 * t136 + (-t147 * t211 + t27 * t787 + t28 * t786 + t3 * t402 + t4 * t401) * m(7) - t160 * t689 - m(6) * (t105 * t94 + t106 * t95 + t244 * t689) - (mrSges(6,1) * t244 + Ifges(6,4) * t717 + Ifges(6,2) * t719 + Ifges(6,6) * t706 - t696 + t729) * t256; (t3 * t506 + t4 * t512 + (-t27 * t506 + t28 * t512) * qJD(6)) * t743 - t98 * t685 - m(7) * (t147 * t685 + t27 * t29 + t28 * t30) + t138 * t716 + (Ifges(6,1) * t562 - t657) * t717 + (Ifges(6,5) * t562 - Ifges(6,6) * t256) * t706 - t30 * t135 - t29 * t136 - t244 * (mrSges(6,1) * t256 + mrSges(6,2) * t562) + (t663 + t229) * t95 + (t664 - t228) * t94 + (-Ifges(6,2) * t256 + t139 + t251) * t719 + (m(7) * t683 + t767) * t676 + (t374 * t743 + t764) * g(2) + (-t376 * t743 + t765) * g(1) + (t135 * t597 - t136 * t598 + t25 * t512 + t26 * t506) * pkin(5) + t815; -t147 * (mrSges(7,1) * t154 + mrSges(7,2) * t563) + (Ifges(7,1) * t563 - t656) * t723 + t85 * t722 + (Ifges(7,5) * t563 - Ifges(7,6) * t154) * t708 - t27 * t135 + t28 * t136 - g(1) * t617 - g(2) * t618 - t549 * t676 + (t154 * t28 + t27 * t563) * mrSges(7,3) + t535 + (-Ifges(7,2) * t154 + t146 + t86) * t725;];
tau  = t1;
