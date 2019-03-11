% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRRP3
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:34
% EndTime: 2019-03-09 00:07:06
% DurationCPUTime: 17.46s
% Computational Cost: add. (22022->803), mult. (52083->1090), div. (0->0), fcn. (54828->10), ass. (0->414)
t505 = sin(qJ(5));
t506 = sin(qJ(4));
t509 = cos(qJ(4));
t766 = cos(qJ(5));
t602 = t766 * t509;
t541 = t505 * t506 - t602;
t542 = t505 * t509 + t506 * t766;
t330 = mrSges(7,1) * t541 + mrSges(7,2) * t542;
t762 = pkin(4) * t509;
t492 = -pkin(3) - t762;
t400 = pkin(5) * t541 + t492;
t884 = m(7) * t400 + t330;
t866 = Ifges(6,6) + Ifges(7,6);
t867 = Ifges(6,5) + Ifges(7,5);
t567 = -t541 * t867 - t542 * t866;
t790 = -pkin(10) - pkin(9);
t477 = t790 * t506;
t478 = t790 * t509;
t818 = t766 * t477 + t505 * t478;
t841 = t818 * mrSges(6,2);
t833 = -t542 * qJ(6) + t818;
t863 = t833 * mrSges(7,2);
t362 = t505 * t477 - t766 * t478;
t864 = t362 * mrSges(6,1);
t276 = -t541 * qJ(6) + t362;
t878 = t276 * mrSges(7,1);
t883 = t567 - t841 - t863 - t864 - t878;
t882 = t841 / 0.2e1 + t863 / 0.2e1 + t864 / 0.2e1 + t878 / 0.2e1;
t507 = sin(qJ(3));
t510 = cos(qJ(3));
t468 = -pkin(3) * t510 - pkin(9) * t507 - pkin(2);
t642 = t509 * t510;
t398 = pkin(8) * t642 + t468 * t506;
t646 = t506 * t507;
t349 = -pkin(10) * t646 + t398;
t317 = t505 * t349;
t451 = t509 * t468;
t644 = t507 * t509;
t565 = -pkin(10) * t644 + t451;
t645 = t506 * t510;
t627 = pkin(8) * t645;
t348 = t565 - t627;
t176 = t766 * t348 - t317;
t424 = t505 * t646 - t507 * t602;
t404 = t424 * qJ(6);
t118 = t404 + t176;
t797 = -m(7) / 0.2e1;
t312 = (-pkin(8) * t506 - pkin(4)) * t510 + t565;
t150 = t766 * t312 - t317;
t102 = t150 + t404;
t97 = -pkin(5) * t510 + t102;
t881 = (t118 - t97) * t797;
t504 = sin(pkin(6));
t508 = sin(qJ(2));
t650 = t504 * t508;
t701 = cos(pkin(6));
t438 = t507 * t701 + t510 * t650;
t511 = cos(qJ(2));
t649 = t504 * t511;
t342 = -t438 * t506 - t509 * t649;
t343 = t438 * t509 - t506 * t649;
t578 = t766 * t342 - t343 * t505;
t880 = t578 / 0.2e1;
t168 = t505 * t342 + t343 * t766;
t425 = t542 * t507;
t713 = t425 * mrSges(6,3);
t365 = mrSges(6,2) * t510 - t713;
t785 = t365 / 0.2e1;
t712 = t425 * mrSges(7,3);
t364 = mrSges(7,2) * t510 - t712;
t786 = t364 / 0.2e1;
t582 = t786 + t785;
t825 = mrSges(6,3) + mrSges(7,3);
t877 = t582 * t578 + t825 * (t168 * t424 / 0.2e1 + t425 * t880);
t794 = m(7) * pkin(5);
t633 = t794 / 0.2e1;
t319 = t766 * t349;
t151 = t505 * t312 + t319;
t671 = t425 * qJ(6);
t103 = t151 - t671;
t175 = -t348 * t505 - t319;
t117 = t175 + t671;
t876 = t103 + t117;
t875 = -t150 + t176;
t702 = t102 - t97;
t714 = t424 * mrSges(7,3);
t368 = -mrSges(7,1) * t510 + t714;
t784 = -t368 / 0.2e1;
t796 = m(7) / 0.2e1;
t874 = t702 * t796 + t784;
t479 = pkin(3) * t507 - pkin(9) * t510;
t401 = pkin(8) * t646 + t509 * t479;
t316 = pkin(4) * t507 - pkin(10) * t642 + t401;
t402 = -pkin(8) * t644 + t506 * t479;
t355 = -pkin(10) * t645 + t402;
t163 = t766 * t316 - t355 * t505;
t427 = t541 * t510;
t100 = pkin(5) * t507 + qJ(6) * t427 + t163;
t870 = -t100 / 0.2e1;
t775 = t542 / 0.2e1;
t816 = t506 ^ 2 + t509 ^ 2;
t869 = pkin(9) * t816;
t827 = mrSges(6,1) + mrSges(7,1);
t868 = mrSges(7,2) + mrSges(6,2);
t865 = Ifges(6,3) + Ifges(7,3);
t715 = t424 * mrSges(6,3);
t369 = -mrSges(6,1) * t510 + t715;
t782 = -t369 / 0.2e1;
t862 = t362 * t782;
t437 = t507 * t650 - t510 * t701;
t861 = t437 * t796;
t496 = Ifges(5,4) * t509;
t560 = -t506 * Ifges(5,1) - t496;
t768 = t509 / 0.2e1;
t860 = t560 * t768;
t840 = t168 * t542;
t859 = t825 * t840;
t858 = t151 + t175;
t747 = Ifges(7,4) * t542;
t335 = -Ifges(7,2) * t541 + t747;
t749 = Ifges(6,4) * t542;
t337 = -Ifges(6,2) * t541 + t749;
t857 = t335 + t337;
t448 = Ifges(7,4) * t541;
t339 = Ifges(7,1) * t542 - t448;
t449 = Ifges(6,4) * t541;
t341 = Ifges(6,1) * t542 - t449;
t856 = t341 + t339;
t626 = t766 * pkin(4);
t491 = t626 + pkin(5);
t853 = t626 - t491;
t639 = t510 * t511;
t389 = (-t506 * t639 + t508 * t509) * t504;
t390 = (t506 * t508 + t509 * t639) * t504;
t217 = t389 * t766 - t505 * t390;
t218 = t505 * t389 + t390 * t766;
t610 = t507 * t649;
t680 = t390 * t509;
t681 = t389 * t506;
t799 = -m(6) / 0.2e1;
t852 = t825 * (t217 * t775 + t218 * t541 / 0.2e1) - (t680 / 0.2e1 - t681 / 0.2e1) * mrSges(5,3) - m(5) * (-pkin(3) * t610 + (t680 - t681) * pkin(9)) / 0.2e1 + (t217 * t818 + t362 * t218 + t492 * t610) * t799 + (t217 * t833 + t276 * t218 + t400 * t610) * t797;
t443 = t541 * mrSges(7,2);
t328 = mrSges(7,1) * t542 - t443;
t329 = mrSges(6,1) * t542 - mrSges(6,2) * t541;
t334 = -Ifges(7,2) * t542 - t448;
t336 = -Ifges(6,2) * t542 - t449;
t753 = mrSges(7,3) * t541;
t851 = t400 * t328 + t492 * t329 + t753 * t833 - (mrSges(7,3) * t833 + t334 / 0.2e1 + t336 / 0.2e1 + t339 / 0.2e1 + t341 / 0.2e1) * t541;
t242 = t542 * t437;
t850 = t242 / 0.2e1;
t777 = t437 / 0.2e1;
t764 = pkin(4) * t505;
t838 = t764 * t825;
t629 = pkin(4) * t644;
t761 = pkin(5) * t424;
t352 = t629 - t761;
t459 = mrSges(5,2) * t510 - mrSges(5,3) * t646;
t581 = t368 / 0.2e1 + t369 / 0.2e1;
t461 = -mrSges(5,1) * t510 - mrSges(5,3) * t644;
t773 = t461 / 0.2e1;
t835 = -t342 * t459 / 0.2e1 + t343 * t773 + (t876 * t797 + t858 * t799) * t578 + (t352 * t797 + t629 * t799) * t437 + (t875 * t799 + t581 + t881) * t168;
t243 = t541 * t437;
t569 = t827 * t850 - t868 * t243 / 0.2e1;
t795 = m(6) * pkin(4);
t634 = t795 / 0.2e1;
t666 = t437 * t506;
t689 = t243 * t505;
t754 = mrSges(5,2) * t509;
t793 = mrSges(5,1) / 0.2e1;
t834 = -(pkin(4) * t689 + t242 * t491) * t796 - (t242 * t766 + t689) * t634 - t666 * t793 - t754 * t777 - t569;
t606 = t712 / 0.2e1;
t608 = t753 / 0.2e1;
t700 = t103 * t542;
t708 = t542 * mrSges(6,3);
t791 = mrSges(7,3) / 0.2e1;
t748 = Ifges(7,4) * t424;
t244 = -Ifges(7,2) * t425 - t510 * Ifges(7,6) - t748;
t750 = Ifges(6,4) * t424;
t246 = -Ifges(6,2) * t425 - t510 * Ifges(6,6) - t750;
t410 = Ifges(7,4) * t425;
t248 = -Ifges(7,1) * t424 - Ifges(7,5) * t510 - t410;
t411 = Ifges(6,4) * t425;
t250 = -Ifges(6,1) * t424 - Ifges(6,5) * t510 - t411;
t405 = t425 * mrSges(7,2);
t277 = -t424 * mrSges(7,1) - t405;
t278 = -mrSges(6,1) * t424 - mrSges(6,2) * t425;
t283 = Ifges(7,2) * t424 - t410;
t284 = Ifges(6,2) * t424 - t411;
t285 = -Ifges(7,1) * t425 + t748;
t286 = -Ifges(6,1) * t425 + t750;
t498 = t507 * pkin(8);
t463 = pkin(4) * t646 + t498;
t314 = pkin(5) * t425 + t463;
t338 = -Ifges(7,1) * t541 - t747;
t340 = -Ifges(6,1) * t541 - t749;
t803 = t276 * t714 / 0.2e1 + t362 * t715 / 0.2e1 + t818 * t713 / 0.2e1 + t492 * t278 / 0.2e1 + t463 * t329 / 0.2e1 + t400 * t277 / 0.2e1 + t314 * t328 / 0.2e1 - (t340 + t338) * t424 / 0.4e1 + t857 * t424 / 0.4e1 - t567 * t510 / 0.4e1 - (t246 + t244) * t542 / 0.4e1 + (t286 + t285) * t542 / 0.4e1 - (t336 + t334 + t856) * t425 / 0.4e1 - (t284 + t283 + t250 + t248) * t541 / 0.4e1;
t522 = t803 + t833 * t606 + t97 * t608 - t700 * t791 - t151 * t708 / 0.2e1;
t604 = t708 / 0.2e1;
t571 = t151 * t604;
t601 = t700 / 0.2e1;
t776 = -t541 / 0.2e1;
t830 = (t102 * t776 + t601) * mrSges(7,3) + t862 + t833 * t786 + t818 * t785 + t571 + t522 + t874 * t276;
t166 = t505 * t316 + t766 * t355;
t426 = t542 * t510;
t112 = -qJ(6) * t426 + t166;
t829 = t166 * mrSges(6,2) / 0.2e1 - t163 * mrSges(6,1) / 0.2e1 + t112 * mrSges(7,2) / 0.2e1 + mrSges(7,1) * t870;
t828 = m(6) + m(7);
t744 = Ifges(5,6) * t506;
t746 = Ifges(5,5) * t509;
t824 = Ifges(4,4) - t746 / 0.2e1 + t744 / 0.2e1;
t331 = mrSges(6,1) * t541 + mrSges(6,2) * t542;
t820 = t330 + t331;
t817 = -Ifges(5,2) * t506 + t496;
t279 = mrSges(7,1) * t425 - mrSges(7,2) * t424;
t787 = -t330 / 0.2e1;
t815 = t279 * t775 + t424 * t787;
t813 = t868 * t766;
t568 = t424 * t866 - t425 * t867;
t812 = -t401 * t506 + t402 * t509;
t704 = t510 * mrSges(4,2);
t472 = t507 * mrSges(4,1) + t704;
t810 = -mrSges(5,1) * t509 + mrSges(5,2) * t506;
t807 = t424 * t838;
t728 = t176 * mrSges(6,2);
t729 = t175 * mrSges(6,1);
t738 = t118 * mrSges(7,2);
t739 = t117 * mrSges(7,1);
t806 = t739 / 0.2e1 - t738 / 0.2e1 + t729 / 0.2e1 - t728 / 0.2e1 + (t117 * t796 + t606) * pkin(5);
t770 = t507 / 0.2e1;
t778 = -t427 / 0.2e1;
t779 = -t426 / 0.2e1;
t805 = t770 * t865 + t778 * t867 + t779 * t866 - t829;
t804 = (t277 + t278) * t777 + t877;
t802 = 0.2e1 * m(7);
t801 = 2 * qJD(3);
t800 = m(5) / 0.2e1;
t798 = m(6) / 0.2e1;
t792 = -mrSges(5,2) / 0.2e1;
t789 = -t578 / 0.2e1;
t370 = mrSges(7,1) * t507 + mrSges(7,3) * t427;
t780 = t370 / 0.2e1;
t774 = -t542 / 0.2e1;
t772 = -t506 / 0.2e1;
t771 = t506 / 0.2e1;
t769 = -t509 / 0.2e1;
t763 = pkin(4) * t506;
t760 = pkin(5) * t542;
t759 = pkin(9) * t506;
t758 = pkin(9) * t509;
t499 = t510 * pkin(8);
t756 = mrSges(4,2) * t507;
t752 = mrSges(7,3) * t542;
t751 = Ifges(5,4) * t506;
t742 = t102 * mrSges(7,2);
t741 = t103 * mrSges(7,1);
t737 = t150 * mrSges(6,2);
t736 = t151 * mrSges(6,1);
t734 = t168 * mrSges(6,1);
t733 = t168 * mrSges(7,1);
t731 = t578 * mrSges(6,2);
t730 = t578 * mrSges(7,2);
t711 = t426 * mrSges(7,1);
t710 = t427 * mrSges(7,2);
t709 = t541 * mrSges(6,3);
t707 = t491 * mrSges(7,3);
t703 = -mrSges(4,1) + t810;
t698 = t578 * t424;
t695 = t168 * t425;
t310 = t437 * t438;
t682 = t343 * t509;
t683 = t342 * t506;
t25 = m(5) * (t310 + (-t682 + t683) * t437) + (t168 * t243 + t242 * t578 + t310) * t828;
t688 = t25 * qJD(1);
t363 = t437 * t610;
t26 = m(5) * (t342 * t389 + t343 * t390 + t363) + (t168 * t218 + t217 * t578 + t363) * t828 + m(4) * (t437 * t507 + t438 * t510 - t650) * t649;
t687 = t26 * qJD(1);
t684 = t314 * t542;
t678 = t400 * t424;
t664 = t541 * t168;
t659 = t542 * t578;
t653 = t491 * t424;
t652 = t491 * t542;
t419 = -t510 * Ifges(5,6) + t507 * t817;
t647 = t506 * t419;
t561 = Ifges(5,1) * t509 - t751;
t538 = t561 * t507;
t421 = -t510 * Ifges(5,5) + t538;
t643 = t509 * t421;
t464 = pkin(4) * t645 + t499;
t628 = t437 * t760;
t621 = t764 / 0.2e1;
t620 = t764 / 0.4e1;
t619 = pkin(5) * t850;
t618 = mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1;
t617 = -mrSges(7,2) / 0.2e1 - mrSges(6,2) / 0.2e1;
t616 = t791 + mrSges(6,3) / 0.2e1;
t615 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t614 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t613 = mrSges(7,1) + t794;
t609 = mrSges(5,3) * t770;
t607 = -t753 / 0.2e1;
t605 = -t712 / 0.2e1;
t603 = m(7) * t702;
t587 = -t644 / 0.2e1;
t585 = t244 / 0.2e1 + t246 / 0.2e1;
t584 = t248 / 0.2e1 + t250 / 0.2e1;
t583 = -t278 / 0.2e1 - t277 / 0.2e1;
t577 = mrSges(6,3) * t626;
t575 = -mrSges(6,1) - t613;
t572 = t626 / 0.2e1;
t570 = (t328 + t329) * t777 - t859 / 0.2e1;
t566 = t610 / 0.2e1;
t564 = t617 * t218;
t563 = t425 * t577;
t471 = mrSges(5,1) * t506 + t754;
t473 = t509 * Ifges(5,2) + t751;
t558 = -t744 + t746;
t557 = Ifges(5,5) * t506 + Ifges(5,6) * t509;
t280 = mrSges(6,1) * t425 - mrSges(6,2) * t424;
t281 = -t710 + t711;
t282 = mrSges(6,1) * t426 - mrSges(6,2) * t427;
t315 = pkin(5) * t426 + t464;
t366 = -mrSges(7,2) * t507 - mrSges(7,3) * t426;
t367 = -mrSges(6,2) * t507 - mrSges(6,3) * t426;
t371 = mrSges(6,1) * t507 + mrSges(6,3) * t427;
t440 = t471 * t507;
t441 = t471 * t510;
t460 = -mrSges(5,2) * t507 - mrSges(5,3) * t645;
t462 = mrSges(5,1) * t507 - mrSges(5,3) * t642;
t397 = t451 - t627;
t551 = t397 * t506 - t398 * t509;
t513 = t582 * t243 + (t367 / 0.2e1 + t366 / 0.2e1) * t168 + t581 * t242 + (t780 + t371 / 0.2e1) * t578 + (t280 / 0.2e1 + t279 / 0.2e1 + t440 / 0.2e1) * t438 + (t282 / 0.2e1 + t281 / 0.2e1 + t441 / 0.2e1 + t459 * t769 + t461 * t771) * t437 + (t438 * t498 + t342 * t401 + t343 * t402 + (t551 + t499) * t437) * t800 + (t150 * t242 + t151 * t243 + t163 * t578 + t166 * t168 + t437 * t464 + t438 * t463) * t798 + (t100 * t578 + t103 * t243 + t112 * t168 + t242 * t97 + t314 * t438 + t315 * t437) * t796 + t342 * t462 / 0.2e1 + t343 * t460 / 0.2e1;
t6 = t513 + (t704 / 0.2e1 - t472 / 0.2e1 + (mrSges(4,1) / 0.2e1 - t810 / 0.2e1 - t331 / 0.2e1 + t787) * t507) * t649 + t852;
t245 = -Ifges(7,4) * t427 - Ifges(7,2) * t426 + Ifges(7,6) * t507;
t247 = -Ifges(6,4) * t427 - Ifges(6,2) * t426 + Ifges(6,6) * t507;
t249 = -Ifges(7,1) * t427 - Ifges(7,4) * t426 + Ifges(7,5) * t507;
t251 = -Ifges(6,1) * t427 - Ifges(6,4) * t426 + Ifges(6,5) * t507;
t420 = Ifges(5,6) * t507 + t510 * t817;
t422 = Ifges(5,5) * t507 + t510 * t561;
t532 = t426 * t614 - t427 * t615;
t7 = (t643 / 0.2e1 - t647 / 0.2e1 + pkin(8) * t440 + t532 + t824 * t510) * t510 + m(5) * (t397 * t401 + t398 * t402) + t397 * t462 + t463 * t282 + t464 * t280 - pkin(2) * t472 + t402 * t459 + t398 * t460 + t401 * t461 + m(6) * (t150 * t163 + t151 * t166 + t463 * t464) + m(7) * (t100 * t97 + t103 * t112 + t314 * t315) + t151 * t367 + t100 * t368 + t163 * t369 + t97 * t370 + t150 * t371 + t112 * t364 + t166 * t365 + t103 * t366 + t314 * t281 + t315 * t279 + (-t251 / 0.2e1 - t249 / 0.2e1) * t424 + (t422 * t768 + t420 * t772 + pkin(8) * t441 - t824 * t507 - t614 * t425 + t615 * t424 + (m(5) * pkin(8) ^ 2 + Ifges(4,1) - Ifges(4,2) - Ifges(5,3) - t865) * t510) * t507 - (t245 / 0.2e1 + t247 / 0.2e1) * t425 - t584 * t427 - t585 * t426;
t556 = t6 * qJD(1) + t7 * qJD(2);
t525 = t314 * t277 + t463 * t278 - (t284 / 0.2e1 + t283 / 0.2e1 + t584) * t425 + t150 * t713 + t97 * t712 - t568 * t510 / 0.2e1;
t530 = t103 * mrSges(7,3) - t286 / 0.2e1 - t285 / 0.2e1 + t151 * mrSges(6,3) + t585;
t10 = m(6) * (t150 * t175 + t151 * t176) + m(7) * (t103 * t118 + t117 * t97 + t314 * t352) + t530 * t424 + t397 * t459 - t398 * t461 + t117 * t368 + t175 * t369 + t118 * t364 + t176 * t365 + t352 * t279 + (t510 * t557 / 0.2e1 + t419 * t769 + t421 * t772 + (-pkin(8) * t810 + t473 * t771 + t860) * t507 + (m(6) * t463 + t280) * t762 + t551 * mrSges(5,3)) * t507 + t525;
t205 = t218 * t764;
t521 = t618 * t217 + t564 + (t217 * t626 + t205) * t798 + (t217 * t491 + t205) * t796 + t389 * t793 + t390 * t792;
t9 = t583 * t437 + (t810 * t777 + (t682 / 0.2e1 - t683 / 0.2e1) * mrSges(5,3)) * t507 + t521 + t835 - t877;
t555 = -t9 * qJD(1) + t10 * qJD(2);
t550 = t633 + t618;
t528 = t217 * t550 + t564;
t11 = (t424 * t633 + t583) * t437 + (-t603 / 0.2e1 - t616 * t424 + t581) * t168 + (-t425 * t616 - t582) * t578 + t528;
t13 = -t151 * t369 + t102 * t364 + t150 * t365 + (t603 - t368) * t103 + ((-m(7) * t314 - t279) * pkin(5) + t530) * t424 + t525;
t554 = -t11 * qJD(1) + t13 * qJD(2);
t48 = m(7) * (-t103 * t425 + t424 * t97) + t424 * t368 - t425 * t364;
t58 = (t566 - t698 / 0.2e1 + t695 / 0.2e1) * m(7);
t553 = -qJD(1) * t58 + qJD(2) * t48;
t70 = (-t425 * t620 + t653 / 0.4e1 - t352 / 0.4e1) * t802 - t277;
t423 = t760 + t763;
t82 = (-t541 * t620 - t652 / 0.4e1 - t423 / 0.4e1) * t802 - t328;
t552 = qJD(2) * t70 + qJD(3) * t82;
t182 = t424 * t613 + t405;
t241 = -t542 * t613 + t443;
t549 = qJD(2) * t182 + qJD(3) * t241;
t548 = m(7) * t853;
t543 = t473 * t772 - t860;
t539 = t335 / 0.2e1 + t337 / 0.2e1 - t338 / 0.2e1 - t340 / 0.2e1;
t537 = t810 * t770;
t517 = pkin(4) * t666 * t798 + t423 * t861 + t471 * t777 + t570;
t15 = t517 + t825 * (-t541 * t789 + t578 * t776 + t840 / 0.2e1) + t834;
t512 = -t862 + t609 * t869 + pkin(4) * t331 * t587 - pkin(3) * t537 - (-t176 / 0.2e1 + t150 / 0.2e1) * t709 + t118 * t608 + t510 * t558 / 0.4e1 - t833 * t364 / 0.2e1 + t352 * t787 + t758 * t773 - t509 * t538 / 0.4e1 - t818 * t365 / 0.2e1 - t280 * t763 / 0.2e1 + t459 * t759 / 0.2e1 - t471 * t498 / 0.2e1 + t117 * t752 / 0.2e1 + t473 * t644 / 0.2e1 - t643 / 0.4e1 - t423 * t279 / 0.2e1 + (t817 / 0.4e1 - t560 / 0.2e1) * t646 + ((t463 * t506 + t492 * t644) * pkin(4) + t858 * t818 + t875 * t362) * t799 + (-t784 + t881) * t276 + t175 * t604 + (t314 * t423 + t352 * t400 + t833 * t876) * t797 + t647 / 0.4e1;
t515 = (t366 + t367) * t621 + (t100 * t491 + t112 * t764) * t796 + t491 * t780 + t402 * t792 + t401 * t793 + Ifges(5,3) * t770 + t805 + t371 * t572 + (t163 * t766 + t166 * t505) * t634 - Ifges(5,6) * t645 / 0.2e1 + Ifges(5,5) * t642 / 0.2e1;
t2 = mrSges(7,3) * t601 + t605 * t833 + t97 * t607 + t512 + t515 + t571 - t803;
t20 = t561 * t771 + t817 * t768 - pkin(3) * t471 - t539 * t542 + t543 + t851 + (m(6) * t492 + t331) * t763 + t884 * t423;
t536 = t15 * qJD(1) - t2 * qJD(2) + t20 * qJD(3);
t17 = (-t328 / 0.2e1 - t329 / 0.2e1) * t437 + (t619 - t628 / 0.2e1) * m(7) + t569;
t21 = -(-pkin(5) * t884 + t539) * t542 + t851;
t4 = (-t370 / 0.2e1 + (t684 / 0.2e1 - t678 / 0.2e1 + t870) * m(7) + t815) * pkin(5) + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1) * t507 + t532 + t829 + t830;
t535 = -t17 * qJD(1) + t4 * qJD(2) + t21 * qJD(3);
t523 = (t424 * t774 - t425 * t776) * mrSges(7,3) + (-t103 * t541 - t276 * t425 + t424 * t833 - t542 * t97) * t796 + t364 * t776 + t368 * t774;
t533 = t315 * t796 + t711 / 0.2e1 - t710 / 0.2e1;
t37 = t523 - t533;
t60 = (t438 / 0.4e1 + t659 / 0.4e1 + t664 / 0.4e1) * t802;
t65 = m(7) * (-t276 * t541 - t542 * t833) + (t541 ^ 2 + t542 ^ 2) * mrSges(7,3);
t534 = -qJD(1) * t60 + qJD(2) * t37 + qJD(3) * t65;
t514 = (-t491 * t103 + (t103 * t766 + t505 * t702) * pkin(4)) * t797 + t742 / 0.2e1 + t741 / 0.2e1 + t737 / 0.2e1 + t736 / 0.2e1 + t491 * t605 - t563 / 0.2e1 + (t368 + t369) * t621 - t807 / 0.2e1 - (t364 + t365) * t626 / 0.2e1;
t19 = t514 + t806;
t526 = t868 * t789 + (t548 / 0.2e1 - t827 / 0.2e1) * t168;
t22 = -t168 * t633 - t730 / 0.2e1 - t733 / 0.2e1 - t731 / 0.2e1 - t734 / 0.2e1 - t526;
t261 = t813 * pkin(4) + (-t548 + t827) * t764;
t520 = t853 * t796 * t276 + t491 * t608 - t572 * t753 - t882;
t527 = pkin(5) * t607 + t276 * t633 + t882;
t27 = t520 + t527;
t531 = t22 * qJD(1) + t19 * qJD(2) - t27 * qJD(3) + t261 * qJD(4);
t503 = t510 ^ 2;
t501 = t507 ^ 2;
t469 = t501 * pkin(8) * t649;
t170 = (-t541 * t764 + t423 - t652) * t796;
t138 = t578 * t764;
t113 = (-t425 * t764 + t352 + t653) * t796;
t61 = (-t659 - t664 + t438) * t796;
t59 = (-t695 + t698) * t796 + m(7) * t566;
t36 = t523 + t533;
t24 = t520 - t527 + t567;
t23 = -t168 * t550 + t578 * t617 + t526;
t18 = t628 * t796 + m(7) * t619 + t569 + t570 + t859 / 0.2e1;
t16 = -t514 + t568 + t806;
t14 = t616 * t840 + t517 - t825 * t541 * (t880 + t789) - t834;
t12 = -t761 * t861 + t528 + (t782 + t874) * t168 + t804;
t8 = t343 * mrSges(5,3) * t587 - t437 * t537 + t609 * t683 + t521 + t804 - t835;
t5 = t513 - t472 * t649 + (t810 + t820) * t566 - t852;
t3 = (t780 + (-t678 + t684) * t796 + t815) * pkin(5) + t100 * t633 + t805 + t830;
t1 = -t512 + t515 + t522;
t28 = [qJD(2) * t26 + qJD(3) * t25, t5 * qJD(3) + t8 * qJD(4) + t12 * qJD(5) + t59 * qJD(6) + t687 + (t217 * t368 + t217 * t369 + t218 * t364 + t218 * t365 + t389 * t461 + t390 * t459 + ((-mrSges(4,1) * t510 - mrSges(3,1) + t756) * t508 + (-mrSges(3,2) + (t501 + t503) * mrSges(4,3) + (t279 + t280 + t440) * t507) * t511) * t504 + 0.2e1 * (t150 * t217 + t151 * t218 + t463 * t610) * t798 + 0.2e1 * (t103 * t218 + t97 * t217 + t314 * t610) * t796 + 0.2e1 * (t389 * t397 + t390 * t398 + t469) * t800 + m(4) * (t469 + (pkin(8) * t503 * t511 - pkin(2) * t508) * t504)) * qJD(2), t688 + t5 * qJD(2) + t14 * qJD(4) + t18 * qJD(5) + t61 * qJD(6) + ((t242 * t818 + t243 * t362 + t438 * t492) * t798 + (t242 * t833 + t243 * t276 + t400 * t438) * t796 + (-pkin(3) * t438 - t437 * t869) * t800) * t801 + ((-mrSges(5,3) * t816 + mrSges(4,2)) * t437 + (t703 + t820) * t438 + t825 * (-t242 * t542 - t243 * t541)) * qJD(3), t8 * qJD(2) + t14 * qJD(3) + (-t343 * mrSges(5,1) - t342 * mrSges(5,2) - t730 - t731 - t733 - t734) * qJD(4) + t23 * qJD(5) + 0.2e1 * ((-t168 * t626 + t138) * t798 + (-t168 * t491 + t138) * t796) * qJD(4), t12 * qJD(2) + t18 * qJD(3) + t23 * qJD(4) + (t575 * t168 - t578 * t868) * qJD(5), qJD(2) * t59 + qJD(3) * t61; qJD(3) * t6 - qJD(4) * t9 - qJD(5) * t11 - qJD(6) * t58 - t687, qJD(3) * t7 + qJD(4) * t10 + qJD(5) * t13 + qJD(6) * t48, t1 * qJD(4) + t3 * qJD(5) + t36 * qJD(6) + ((-pkin(3) * t499 + pkin(9) * t812) * t800 + (t163 * t818 + t166 * t362 + t464 * t492) * t798 + (t100 * t833 + t112 * t276 + t315 * t400) * t796) * t801 + t556 + (t856 * t778 + t857 * t779 + t812 * mrSges(5,3) + t833 * t370 + t422 * t771 + t420 * t768 - t100 * t752 - t112 * t753 + pkin(8) * t756 + t460 * t758 + (pkin(8) * t703 + Ifges(4,5) + t543) * t510 - t163 * t708 - t166 * t709 + t818 * t371 - t462 * t759 - Ifges(4,6) * t507 + t492 * t282 + t464 * t331 - pkin(3) * t441 + t400 * t281 + t362 * t367 + t276 * t366 + t315 * t330 + (-t541 * t866 + t542 * t867 + t557) * t770 + (t251 + t249) * t775 + (t247 + t245) * t776) * qJD(3), t1 * qJD(3) + (-Ifges(5,5) * t646 - Ifges(5,6) * t644 + m(7) * (t117 * t491 + t118 * t764) + t425 * t707 + t739 - t738 - t728 + t729 + (t175 * t766 + t176 * t505) * t795 - t398 * mrSges(5,1) - t397 * mrSges(5,2) + t563 + t568 + t807) * qJD(4) + t16 * qJD(5) + t113 * qJD(6) + t555, t3 * qJD(3) + t16 * qJD(4) + (-t736 - t741 - t737 - t742 + (-m(7) * t103 + t712) * pkin(5) + t568) * qJD(5) + t554, qJD(3) * t36 + qJD(4) * t113 + t553; -qJD(2) * t6 + qJD(4) * t15 - qJD(5) * t17 - qJD(6) * t60 - t688, -qJD(4) * t2 + qJD(5) * t4 + qJD(6) * t37 - t556, qJD(4) * t20 + qJD(5) * t21 + qJD(6) * t65 ((-t362 * t766 + t505 * t818) * t795 + m(7) * (-t491 * t276 + t764 * t833) + t558 - t542 * t838 - (-t577 - t707) * t541 + t810 * pkin(9) + t883) * qJD(4) + t24 * qJD(5) + t170 * qJD(6) + t536, t24 * qJD(4) + ((-m(7) * t276 + t753) * pkin(5) + t883) * qJD(5) + t535, qJD(4) * t170 + t534; qJD(2) * t9 - qJD(3) * t15 - qJD(5) * t22, qJD(3) * t2 - qJD(5) * t19 + qJD(6) * t70 - t555, qJD(5) * t27 + qJD(6) * t82 - t536, -t261 * qJD(5) (t505 * t575 - t813) * qJD(5) * pkin(4) - t531, t552; qJD(2) * t11 + qJD(3) * t17 + qJD(4) * t22, -qJD(3) * t4 + qJD(4) * t19 + qJD(6) * t182 - t554, -qJD(4) * t27 + qJD(6) * t241 - t535, t531, 0, t549; qJD(2) * t58 + qJD(3) * t60, -qJD(3) * t37 - qJD(4) * t70 - qJD(5) * t182 - t553, -qJD(4) * t82 - qJD(5) * t241 - t534, -t552, -t549, 0;];
Cq  = t28;
