% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:45
% EndTime: 2019-03-09 05:23:12
% DurationCPUTime: 16.63s
% Computational Cost: add. (27395->789), mult. (56362->1065), div. (0->0), fcn. (60106->8), ass. (0->388)
t482 = sin(qJ(6));
t485 = cos(qJ(6));
t483 = sin(qJ(4));
t486 = cos(qJ(4));
t615 = sin(pkin(10));
t616 = cos(pkin(10));
t508 = t483 * t616 + t486 * t615;
t756 = t615 * t483 - t616 * t486;
t326 = -t482 * t756 + t485 * t508;
t532 = -t482 * t508 - t485 * t756;
t582 = Ifges(7,5) * t532 - Ifges(7,6) * t326;
t666 = -qJ(5) - pkin(8);
t444 = t666 * t483;
t446 = t666 * t486;
t353 = t615 * t444 - t616 * t446;
t279 = pkin(9) * t756 - t353;
t576 = t616 * t444 + t615 * t446;
t753 = -pkin(9) * t508 + t576;
t156 = t279 * t485 - t482 * t753;
t788 = t279 * t482 + t485 * t753;
t811 = t156 * mrSges(7,1) - t788 * mrSges(7,2);
t24 = t582 + t811;
t813 = t24 * qJD(6);
t487 = cos(qJ(3));
t676 = pkin(8) * t487;
t484 = sin(qJ(3));
t679 = pkin(3) * t484;
t443 = qJ(2) - t676 + t679;
t425 = t486 * t443;
t488 = -pkin(1) - pkin(7);
t542 = -t483 * t488 + pkin(4);
t588 = t486 * t487;
t569 = qJ(5) * t588;
t322 = t484 * t542 + t425 - t569;
t475 = t484 * t488;
t364 = t483 * t443 + t475 * t486;
t594 = t483 * t487;
t345 = -qJ(5) * t594 + t364;
t327 = t615 * t345;
t185 = t616 * t322 - t327;
t397 = t756 * t487;
t675 = t397 * pkin(9);
t137 = t484 * pkin(5) + t185 + t675;
t540 = t616 * t345;
t186 = t615 * t322 + t540;
t399 = t508 * t487;
t674 = t399 * pkin(9);
t151 = t186 - t674;
t66 = t137 * t485 - t151 * t482;
t363 = -t475 * t483 + t425;
t344 = t363 - t569;
t203 = -t344 * t615 - t540;
t163 = t203 + t674;
t204 = t616 * t344 - t327;
t164 = t204 + t675;
t78 = t163 * t482 + t164 * t485;
t812 = t78 - t66;
t67 = t137 * t482 + t151 * t485;
t77 = t163 * t485 - t164 * t482;
t809 = t67 + t77;
t464 = pkin(4) * t594;
t587 = t487 * t488;
t429 = t464 - t587;
t332 = t399 * pkin(5) + t429;
t270 = t397 * t485 + t399 * t482;
t534 = t397 * t482 - t485 * t399;
t765 = t534 * mrSges(7,2);
t790 = -t270 * mrSges(7,1) + t765;
t807 = t332 * t790;
t677 = pkin(4) * t486;
t469 = -pkin(3) - t677;
t376 = pkin(5) * t756 + t469;
t766 = t532 * mrSges(7,2);
t791 = t326 * mrSges(7,1) + t766;
t806 = t376 * t791;
t398 = t508 * t484;
t400 = t756 * t484;
t272 = t485 * t398 - t400 * t482;
t276 = -t398 * t482 - t485 * t400;
t39 = -t276 * mrSges(7,1) + t272 * mrSges(7,2);
t805 = t39 * qJD(6);
t803 = -t185 + t204;
t454 = pkin(3) * t487 + pkin(8) * t484;
t436 = t486 * t454;
t593 = t484 * t486;
t333 = qJ(5) * t593 + t487 * t542 + t436;
t375 = t483 * t454 + t486 * t587;
t595 = t483 * t484;
t346 = qJ(5) * t595 + t375;
t201 = t616 * t333 - t346 * t615;
t202 = t615 * t333 + t616 * t346;
t223 = -mrSges(7,2) * t487 + mrSges(7,3) * t272;
t795 = mrSges(7,3) * t276;
t225 = mrSges(7,1) * t487 + t795;
t354 = -mrSges(6,2) * t487 + mrSges(6,3) * t398;
t356 = mrSges(6,1) * t487 - mrSges(6,3) * t400;
t374 = -t483 * t587 + t436;
t557 = t616 * pkin(4);
t468 = t557 + pkin(5);
t558 = pkin(4) * t615;
t405 = t485 * t468 - t482 * t558;
t406 = t482 * t468 + t485 * t558;
t696 = t400 / 0.2e1;
t699 = t398 / 0.2e1;
t517 = Ifges(6,5) * t696 + Ifges(6,6) * t699;
t797 = -t276 / 0.2e1;
t798 = t272 / 0.2e1;
t516 = Ifges(7,5) * t797 + Ifges(7,6) * t798;
t681 = -t487 / 0.2e1;
t142 = t487 * pkin(5) - t400 * pkin(9) + t201;
t160 = pkin(9) * t398 + t202;
t70 = t142 * t485 - t160 * t482;
t71 = t142 * t482 + t160 * t485;
t526 = -Ifges(7,3) * t681 - t71 * mrSges(7,2) / 0.2e1 + t70 * mrSges(7,1) / 0.2e1 + t516;
t621 = t484 * Ifges(5,5);
t563 = -t621 / 0.2e1;
t738 = m(6) * pkin(4);
t572 = t738 / 0.2e1;
t680 = t487 / 0.2e1;
t695 = t405 / 0.2e1;
t737 = -mrSges(5,2) / 0.2e1;
t740 = m(7) / 0.2e1;
t771 = Ifges(5,3) + Ifges(6,3);
t802 = (t405 * t70 + t406 * t71) * t740 + t201 * mrSges(6,1) / 0.2e1 - t202 * mrSges(6,2) / 0.2e1 + t374 * mrSges(5,1) / 0.2e1 + t375 * t737 + t225 * t695 + t406 * t223 / 0.2e1 + (t201 * t616 + t202 * t615) * t572 + t486 * t563 + Ifges(5,6) * t595 / 0.2e1 + t356 * t557 / 0.2e1 + t354 * t558 / 0.2e1 + t517 + t771 * t680 + t526;
t658 = Ifges(7,4) * t326;
t192 = Ifges(7,2) * t532 + t658;
t793 = Ifges(7,1) * t532 - t658;
t801 = t192 / 0.4e1 - t793 / 0.4e1;
t659 = Ifges(7,4) * t270;
t134 = Ifges(7,2) * t534 + t484 * Ifges(7,6) - t659;
t255 = Ifges(7,4) * t534;
t136 = -Ifges(7,1) * t270 + Ifges(7,5) * t484 + t255;
t224 = -mrSges(7,2) * t484 + mrSges(7,3) * t534;
t728 = t224 / 0.2e1;
t778 = Ifges(7,2) * t270 + t255;
t792 = Ifges(7,1) * t534 + t659;
t800 = (t136 / 0.4e1 + t778 / 0.4e1) * t532 + t788 * t728 + t790 * t376 / 0.2e1 + t791 * t332 / 0.2e1 + (t792 / 0.4e1 - t134 / 0.4e1) * t326;
t799 = 0.2e1 * mrSges(7,1);
t480 = t483 ^ 2;
t481 = t486 ^ 2;
t575 = t480 + t481;
t796 = mrSges(5,3) * t575;
t767 = Ifges(7,5) * t534;
t784 = Ifges(7,6) * t270;
t585 = t767 + t784;
t315 = Ifges(7,4) * t532;
t195 = Ifges(7,1) * t326 + t315;
t777 = -Ifges(7,2) * t326 + t315;
t789 = t777 / 0.4e1 + t195 / 0.4e1;
t545 = -t784 / 0.2e1 - t767 / 0.2e1;
t786 = t270 / 0.2e1;
t785 = -t326 / 0.2e1;
t708 = t326 / 0.2e1;
t725 = -t765 / 0.2e1;
t565 = t766 / 0.2e1;
t758 = t203 + t186;
t518 = -t374 * t483 + t375 * t486;
t776 = t270 * t67 - t534 * t66;
t357 = mrSges(6,1) * t484 + t397 * mrSges(6,3);
t773 = t357 / 0.2e1;
t772 = t532 / 0.2e1;
t623 = t508 * mrSges(6,3);
t447 = t483 * mrSges(5,1) + t486 * mrSges(5,2);
t761 = t484 * t447;
t283 = -t397 * mrSges(6,1) - t399 * mrSges(6,2);
t760 = t790 + t283;
t338 = mrSges(6,1) * t508 - mrSges(6,2) * t756;
t759 = t791 + t338;
t479 = Ifges(5,4) * t486;
t757 = -Ifges(5,2) * t483 + t479;
t450 = Ifges(5,1) * t483 + t479;
t607 = t353 * t397;
t608 = t576 * t399;
t755 = t607 / 0.2e1 + t608 / 0.2e1;
t422 = t487 * t450;
t438 = -mrSges(5,2) * t484 - mrSges(5,3) * t594;
t739 = -pkin(8) / 0.2e1;
t754 = t438 * t739 - t422 / 0.4e1;
t752 = -Ifges(6,5) * t399 + Ifges(6,6) * t397 + t585;
t751 = -Ifges(6,5) * t756 - Ifges(6,6) * t508 + t582;
t750 = -t400 * t397 - t398 * t399;
t619 = t486 * mrSges(5,1);
t622 = t483 * mrSges(5,2);
t749 = t622 - t619;
t440 = mrSges(5,1) * t484 - mrSges(5,3) * t588;
t589 = t486 * t440;
t748 = -t589 / 0.2e1 + t681 * t796;
t633 = t326 * mrSges(7,2);
t637 = t532 * mrSges(7,1);
t583 = t637 / 0.2e1 - t633 / 0.2e1;
t640 = t270 * mrSges(7,2);
t646 = t534 * mrSges(7,1);
t586 = t646 / 0.2e1 + t640 / 0.2e1;
t146 = -t640 - t646;
t189 = t633 - t637;
t661 = Ifges(6,4) * t397;
t266 = -Ifges(6,2) * t399 + t484 * Ifges(6,6) - t661;
t287 = -Ifges(6,1) * t399 + t661;
t660 = Ifges(6,4) * t508;
t341 = -Ifges(6,2) * t756 + t660;
t342 = -Ifges(6,1) * t756 - t660;
t349 = pkin(4) * t588 - pkin(5) * t397;
t678 = pkin(4) * t483;
t381 = pkin(5) * t508 + t678;
t419 = t749 * t487;
t684 = t484 / 0.4e1;
t687 = t469 / 0.2e1;
t701 = t381 / 0.2e1;
t355 = -mrSges(6,2) * t484 - t399 * mrSges(6,3);
t703 = t355 / 0.2e1;
t704 = t349 / 0.2e1;
t226 = mrSges(7,1) * t484 + t270 * mrSges(7,3);
t726 = t226 / 0.2e1;
t732 = -t788 / 0.2e1;
t735 = -t66 / 0.2e1;
t746 = t800 + t751 * t684 + t576 * t703 + t429 * t338 / 0.2e1 + pkin(3) * t419 / 0.2e1 - (t266 / 0.4e1 - t287 / 0.4e1) * t508 + (t341 / 0.4e1 - t342 / 0.4e1) * t397 + t801 * t270 + (t332 * t381 + t349 * t376 + t809 * t788) * t740 + (t532 * t735 + t534 * t732 + t772 * t78 + t809 * t785) * mrSges(7,3) + t789 * t534 - t353 * t773 + t283 * t687 + t146 * t701 + t189 * t704 + (-t786 * mrSges(7,3) - t740 * t812 + t726) * t156;
t745 = 2 * qJD(3);
t744 = m(5) / 0.2e1;
t743 = -m(6) / 0.2e1;
t742 = m(6) / 0.2e1;
t741 = -m(7) / 0.2e1;
t734 = t136 / 0.2e1;
t731 = t156 / 0.2e1;
t729 = t204 / 0.2e1;
t727 = -t226 / 0.2e1;
t723 = t276 / 0.2e1;
t719 = t534 / 0.2e1;
t718 = -t272 / 0.2e1;
t714 = -t270 / 0.2e1;
t700 = -t397 / 0.2e1;
t698 = -t399 / 0.2e1;
t662 = Ifges(5,4) * t483;
t448 = Ifges(5,2) * t486 + t662;
t421 = t487 * t448;
t694 = -t421 / 0.4e1;
t692 = -t756 / 0.2e1;
t690 = -t508 / 0.2e1;
t689 = -t448 / 0.4e1;
t451 = Ifges(5,1) * t486 - t662;
t688 = t451 / 0.4e1;
t686 = t483 / 0.2e1;
t685 = t484 / 0.2e1;
t683 = -t486 / 0.2e1;
t682 = t486 / 0.2e1;
t672 = t66 * mrSges(7,2);
t671 = t67 * mrSges(7,1);
t668 = t77 * mrSges(7,1);
t667 = t78 * mrSges(7,2);
t665 = mrSges(7,3) * t405;
t664 = mrSges(7,3) * t406;
t478 = Ifges(5,5) * t486;
t656 = Ifges(5,6) * t483;
t647 = t272 * mrSges(7,1);
t642 = t276 * mrSges(7,2);
t133 = -Ifges(7,4) * t276 + Ifges(7,2) * t272 + t487 * Ifges(7,6);
t135 = -Ifges(7,1) * t276 + Ifges(7,4) * t272 + t487 * Ifges(7,5);
t145 = -t642 - t647;
t265 = Ifges(6,4) * t400 + Ifges(6,2) * t398 + t487 * Ifges(6,6);
t267 = Ifges(6,1) * t400 + Ifges(6,4) * t398 + Ifges(6,5) * t487;
t380 = Ifges(6,4) * t399;
t268 = -Ifges(6,1) * t397 + Ifges(6,5) * t484 - t380;
t627 = t400 * mrSges(6,2);
t630 = t398 * mrSges(6,1);
t284 = t627 - t630;
t629 = t399 * mrSges(6,1);
t631 = t397 * mrSges(6,2);
t285 = t629 - t631;
t428 = -pkin(4) * t595 + t475;
t331 = -pkin(5) * t398 + t428;
t393 = t487 * Ifges(5,6) - t484 * t757;
t395 = t487 * Ifges(5,5) - t451 * t484;
t437 = -mrSges(5,2) * t487 + mrSges(5,3) * t595;
t439 = mrSges(5,1) * t487 + mrSges(5,3) * t593;
t513 = t478 / 0.2e1 - t656 / 0.2e1 - Ifges(4,4);
t396 = t451 * t487 + t621;
t591 = t486 * t396;
t620 = t484 * Ifges(5,6);
t394 = t487 * t757 + t620;
t597 = t483 * t394;
t3 = t429 * t284 + t364 * t437 + t375 * t438 + t363 * t439 + t374 * t440 + t428 * t285 + t201 * t357 + t186 * t354 + t202 * t355 + t185 * t356 + t331 * t146 + t332 * t145 + t67 * t223 + t71 * t224 + t66 * t225 + t70 * t226 + (Ifges(7,5) * t714 + Ifges(7,6) * t719 + qJ(2) * mrSges(4,1) + Ifges(6,5) * t700 + Ifges(6,6) * t698 + t488 * t761 + t395 * t682 - t483 * t393 / 0.2e1 + t513 * t487 + (Ifges(7,3) - Ifges(4,1) + Ifges(4,2) + (-m(5) * t488 + t447) * t488 + t771) * t484) * t487 + t136 * t797 + t134 * t798 + (-t591 / 0.2e1 + t597 / 0.2e1 - qJ(2) * mrSges(4,2) - t513 * t484 + t516 + t517) * t484 + m(5) * (t363 * t374 + t364 * t375) + t268 * t696 + t265 * t698 + t266 * t699 + t267 * t700 + t135 * t714 + t133 * t719 + m(6) * (t185 * t201 + t186 * t202 + t428 * t429) + m(7) * (t331 * t332 + t66 * t70 + t67 * t71);
t639 = t3 * qJD(1);
t635 = t532 * mrSges(7,3);
t632 = t326 * mrSges(7,3);
t286 = Ifges(6,2) * t397 - t380;
t4 = t363 * t438 - t364 * t440 + t429 * t283 + t203 * t357 + t204 * t355 + t349 * t146 + t807 + t792 * t714 + t534 * t734 + t778 * t719 + t134 * t786 + t78 * t224 + t77 * t226 + t776 * mrSges(7,3) + m(6) * (t185 * t203 + t186 * t204) + m(7) * (t332 * t349 + t66 * t77 + t67 * t78) - (t286 / 0.2e1 + t268 / 0.2e1 - t185 * mrSges(6,3)) * t399 + (t266 / 0.2e1 - t287 / 0.2e1 + t186 * mrSges(6,3)) * t397 + (t488 * t419 + (t563 + t363 * mrSges(5,3) + t421 / 0.2e1 - t396 / 0.2e1) * t483 + (-t620 / 0.2e1 - t364 * mrSges(5,3) - t422 / 0.2e1 - t394 / 0.2e1 + (m(6) * t429 + t285) * pkin(4)) * t486) * t487 + t752 * t685;
t628 = t4 * qJD(1);
t626 = t756 * mrSges(6,1);
t625 = t756 * mrSges(6,3);
t624 = t508 * mrSges(6,2);
t9 = t807 + t585 * t685 + t66 * t224 - t67 * t226 - (-t67 * mrSges(7,3) - t134 / 0.2e1 + t792 / 0.2e1) * t270 + (-t66 * mrSges(7,3) + t778 / 0.2e1 + t734) * t534;
t618 = t9 * qJD(1);
t617 = -mrSges(4,1) + t749;
t51 = (-t270 * t272 + t276 * t534) * t740 + (-t397 * t398 + t399 * t400) * t742;
t614 = qJD(1) * t51;
t613 = t185 * t756;
t609 = t276 * t270;
t610 = t272 * t534;
t502 = (t609 / 0.2e1 + t610 / 0.2e1) * mrSges(7,3) - t272 * t728 + t226 * t797 + t790 * t681;
t19 = t502 - t583;
t612 = t19 * qJD(1);
t27 = t484 * mrSges(4,1) + t487 * mrSges(4,2) + t326 * t224 + t532 * t226 + t508 * t355 - t756 * t357 + t483 * t438 + t589 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(7) * (t326 * t67 + t532 * t66) + m(6) * (t186 * t508 - t613) + m(5) * (t363 * t486 + t364 * t483);
t611 = t27 * qJD(1);
t600 = t405 * t534;
t599 = t406 * t270;
t598 = t429 * t483;
t596 = t483 * t439;
t592 = t484 * t487;
t590 = t486 * t437;
t574 = qJD(3) * t484;
t115 = -mrSges(7,1) * t406 - mrSges(7,2) * t405;
t573 = t115 * qJD(6);
t571 = t678 / 0.2e1;
t567 = mrSges(7,3) * t798;
t566 = t795 / 0.2e1;
t552 = -t595 / 0.2e1;
t551 = -t594 / 0.2e1;
t548 = t588 / 0.2e1;
t547 = t791 * t681;
t546 = -t488 * t447 / 0.2e1;
t536 = t478 - t656;
t531 = pkin(4) * t548;
t530 = mrSges(6,3) * t557;
t529 = mrSges(6,3) * t558;
t46 = m(7) * (-t592 - t609 - t610) + m(6) * (-t592 - t750) + m(5) * (-0.1e1 + t575) * t592;
t491 = -m(5) * ((-t363 * t483 + t364 * t486) * t487 + (t518 - 0.2e1 * t587) * t484) / 0.2e1 + (-t185 * t399 - t186 * t397 - t201 * t398 - t202 * t400 - t428 * t487 + t429 * t484) * t743 + (-t272 * t70 + t276 * t71 - t331 * t487 + t332 * t484 - t776) * t741 + t225 * t798 + t534 * t727 + t223 * t797 + t224 * t786 + t397 * t703 + t356 * t699 + t399 * t773 + t354 * t696;
t504 = (-t156 * t326 + t532 * t788) * t740 + (t353 * t508 - t576 * t756) * t742;
t7 = t491 + (t145 / 0.2e1 + t284 / 0.2e1 + t440 * t686 - t761 + t438 * t683) * t487 + (-t146 / 0.2e1 - t285 / 0.2e1 + t596 / 0.2e1 - t590 / 0.2e1) * t484 + t504;
t523 = -t7 * qJD(1) + t46 * qJD(2);
t490 = (-t487 ^ 2 * t677 - t758 * t398 - t803 * t400) * t742 + (-t809 * t272 + t276 * t812 - t349 * t487) * t740 - t276 * t726 + t224 * t718 - t398 * t355 / 0.2e1 + t357 * t696 + t534 * t567 + t270 * t566 + t438 * t552 + t750 * mrSges(6,3) / 0.2e1 + (-t419 + t760) * t681 + t748 * t484;
t497 = (t326 * t406 + t405 * t532) * t740 - t626 / 0.2e1 - t624 / 0.2e1 - t622 / 0.2e1 + t619 / 0.2e1 + (t508 * t615 - t616 * t756) * t572 + t583;
t10 = -t490 + t497;
t522 = t10 * qJD(1);
t30 = t534 * t224 + t270 * t226 - t399 * t355 + t397 * t357 + m(7) * (t270 * t66 + t534 * t67) + m(6) * (t185 * t397 - t186 * t399);
t521 = qJD(1) * t30 + qJD(2) * t51;
t499 = (t270 * t405 + t406 * t534) * t740 + (t397 * t616 - t399 * t615) * t572;
t512 = m(6) * t531 + m(7) * t704;
t45 = t499 - t512 - t760;
t498 = (-t326 * t405 + t406 * t532) * t740 + (-t508 * t616 - t615 * t756) * t572;
t514 = m(6) * t571 + m(7) * t701;
t47 = t498 - t514 - t759;
t520 = qJD(1) * t45 + qJD(3) * t47;
t72 = t786 * t799 + 0.2e1 * t725;
t83 = t708 * t799 + 0.2e1 * t565;
t519 = qJD(1) * t72 - qJD(3) * t83;
t339 = t624 + t626;
t418 = Ifges(6,4) * t756;
t340 = -Ifges(6,2) * t508 - t418;
t343 = Ifges(6,1) * t508 - t418;
t507 = t803 * t353 + t758 * t576;
t1 = t339 * t531 + t285 * t571 + t487 * t546 + t746 + t591 / 0.4e1 - (t450 + t757) * t594 / 0.4e1 - t758 * t623 / 0.2e1 - (t286 + t268) * t756 / 0.4e1 - (t343 + t340) * t399 / 0.4e1 + t754 * t483 + t748 * pkin(8) + (t613 / 0.2e1 + t755) * mrSges(6,3) - t597 / 0.4e1 + t536 * t684 + t588 * t688 + t588 * t689 + t486 * t694 - t625 * t729 + ((t469 * t588 + t598) * pkin(4) + t507) * t742 - t802;
t12 = -pkin(3) * t447 + t806 + t793 * t708 + t192 * t785 + (t757 / 0.2e1 + t450 / 0.2e1) * t486 - (t341 / 0.2e1 - t342 / 0.2e1) * t508 - (t343 / 0.2e1 + t340 / 0.2e1) * t756 + (-t448 / 0.2e1 + t451 / 0.2e1 + pkin(4) * t339) * t483 + (t777 + t195) * t772 + (m(6) * t678 + t338) * t469 + (m(7) * t376 + t189) * t381;
t494 = -t464 * t742 - t381 * t487 * t740 + t632 * t723 + t532 * t567 + t635 * t718 - t326 * t566 + (t447 + t759) * t681;
t496 = (-t599 + t600) * t740 + t631 / 0.2e1 - t629 / 0.2e1 + (-t397 * t615 - t399 * t616) * t572 + mrSges(5,1) * t551 + t588 * t737 + t586;
t17 = -t494 + t496;
t511 = t1 * qJD(1) - t17 * qJD(2) + t12 * qJD(3);
t23 = t806 + (-t192 / 0.2e1 + t793 / 0.2e1) * t326 + (t777 / 0.2e1 + t195 / 0.2e1) * t532;
t32 = t547 - t586;
t495 = (mrSges(7,3) * t732 + t789) * t534 - (mrSges(7,3) * t731 - t801) * t270 - t156 * t727 + t582 * t684 + t800;
t6 = t495 - t526;
t510 = t6 * qJD(1) + t32 * qJD(2) + t23 * qJD(3);
t493 = (t270 * t785 + t534 * t772) * mrSges(7,3) + (t397 * t690 - t399 * t692) * mrSges(6,3) + (-t185 * t508 - t186 * t756 - t353 * t399 + t397 * t576) * t742 + (-t156 * t534 + t270 * t788 - t326 * t66 + t532 * t67) * t740 + t226 * t785 + t224 * t772 + t355 * t692 + t357 * t690;
t500 = t428 * t743 + t331 * t741 + t647 / 0.2e1 + t642 / 0.2e1 + t630 / 0.2e1 - t627 / 0.2e1;
t16 = t493 + t500;
t37 = (t326 ^ 2 + t532 ^ 2) * mrSges(7,3) + (t508 ^ 2 + t756 ^ 2) * mrSges(6,3) + m(7) * (-t156 * t532 - t326 * t788) + m(6) * (-t353 * t756 - t508 * t576);
t503 = (t272 * t326 + t276 * t532) * t740 + (t398 * t508 + t400 * t756) * t742;
t55 = (t741 + t743) * t484 + t503;
t509 = -qJD(1) * t16 - qJD(2) * t55 - qJD(3) * t37;
t501 = (t599 / 0.2e1 - t600 / 0.2e1) * mrSges(7,3) + t224 * t695 + t406 * t727 - t545;
t14 = (t735 + t78 / 0.2e1) * mrSges(7,2) + (-t67 / 0.2e1 - t77 / 0.2e1) * mrSges(7,1) + t501 + t545;
t25 = (t732 + t788 / 0.2e1) * mrSges(7,2) + (t731 - t156 / 0.2e1) * mrSges(7,1);
t38 = (t798 + t718) * mrSges(7,2) + (t797 + t723) * mrSges(7,1);
t506 = t14 * qJD(1) + t38 * qJD(2) + t25 * qJD(3) + t115 * qJD(4);
t84 = t565 - t766 / 0.2e1;
t79 = t498 + t514;
t73 = t765 / 0.2e1 + t725;
t60 = t499 + t512;
t54 = t503 + (m(6) + m(7)) * t685;
t49 = t51 * qJD(5);
t31 = t547 + t586;
t20 = t502 + t583;
t18 = t494 + t496;
t15 = t493 - t500;
t13 = -t672 / 0.2e1 - t671 / 0.2e1 - t667 / 0.2e1 + t668 / 0.2e1 + t501 - t545;
t11 = t490 + t497;
t8 = t438 * t548 + t439 * t552 + t440 * t551 + t761 * t680 - t491 + t504 + (t146 + t285 + t590) * t685 + (t145 + t284 - t761) * t681;
t5 = t495 + t526;
t2 = (pkin(4) * t598 + t507) * t742 + (t694 + t396 / 0.4e1 + t440 * t739) * t486 - (t286 / 0.4e1 + t268 / 0.4e1) * t756 - (t343 / 0.4e1 + t340 / 0.4e1) * t399 + (-t620 / 0.4e1 - t394 / 0.4e1 + pkin(4) * t285 / 0.2e1 + t754) * t483 + t478 * t684 + (t546 + (-t450 / 0.4e1 - t757 / 0.4e1) * t483 + (-t481 / 0.2e1 - t480 / 0.2e1) * pkin(8) * mrSges(5,3) + (t688 + t689 + (t339 / 0.2e1 + m(6) * t687) * pkin(4)) * t486) * t487 + (-(t186 / 0.2e1 + t203 / 0.2e1) * t508 - (t729 - t185 / 0.2e1) * t756 + t755) * mrSges(6,3) + t746 + t802;
t21 = [qJD(2) * t27 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t30 + qJD(6) * t9, t611 + 0.2e1 * ((-t272 * t532 + t276 * t326) * t740 + (t398 * t756 - t400 * t508) * t742) * qJD(2) + t8 * qJD(3) + t11 * qJD(4) + t49 + t20 * qJD(6), t639 + t8 * qJD(2) + t2 * qJD(4) + t15 * qJD(5) + t5 * qJD(6) + (t448 * t686 + t450 * t683 + t488 * t617 - Ifges(4,5)) * t574 + ((-pkin(3) * t475 + pkin(8) * t518) * t744 + (t201 * t576 + t202 * t353 + t428 * t469) * t742 + (-t156 * t71 + t331 * t376 + t70 * t788) * t740) * t745 + (-t70 * t632 + t576 * t356 - Ifges(4,6) * t487 + t518 * mrSges(5,3) + t469 * t284 + t428 * t339 + t376 * t145 + t353 * t354 + t331 * t189 + (-t596 + t590) * pkin(8) + t133 * t772 + pkin(3) * t761 + t508 * t267 / 0.2e1 + t71 * t635 - t201 * t623 - t202 * t625 + t788 * t225 + t195 * t797 + t192 * t798 - t156 * t223 + (Ifges(5,5) * t483 + Ifges(6,5) * t508 + Ifges(7,5) * t326 + Ifges(5,6) * t486 - Ifges(6,6) * t756 + Ifges(7,6) * t532) * t680 + t393 * t682 + t395 * t686 + t265 * t692 + t343 * t696 + t341 * t699 + t135 * t708 - mrSges(4,2) * t587) * qJD(3), t628 + t11 * qJD(2) + t2 * qJD(3) + (-Ifges(5,5) * t594 - Ifges(5,6) * t588 + m(7) * (t405 * t77 + t406 * t78) - t667 + t668 + t270 * t664 - t534 * t665 - t204 * mrSges(6,2) + t203 * mrSges(6,1) + (t203 * t616 + t204 * t615) * t738 - t364 * mrSges(5,1) - t363 * mrSges(5,2) + t397 * t529 + t399 * t530 + t752) * qJD(4) + t60 * qJD(5) + t13 * qJD(6), qJD(3) * t15 + qJD(4) * t60 + qJD(6) * t73 + t521, t618 + t20 * qJD(2) + t5 * qJD(3) + t13 * qJD(4) + t73 * qJD(5) + (t585 - t671 - t672) * qJD(6); -qJD(3) * t7 - qJD(4) * t10 + qJD(6) * t19 + t49 - t611, t46 * qJD(3), t18 * qJD(4) + t54 * qJD(5) + t31 * qJD(6) + (t189 + t339 + t617) * t574 + ((t156 * t270 + t376 * t484 + t534 * t788) * t740 + (t469 * t484 - t607 - t608) * t742 + (t575 * t676 - t679) * t744) * t745 + t523 + (-t534 * t632 - t270 * t635 + t397 * t625 + t399 * t623 + (-mrSges(4,2) + t796) * t487) * qJD(3), t18 * qJD(3) + (m(7) * (-t272 * t406 - t276 * t405) + (-t398 * t615 + t400 * t616) * t738 + t398 * mrSges(6,2) + t400 * mrSges(6,1) + mrSges(5,2) * t595 - mrSges(5,1) * t593 + t39) * qJD(4) + t805 - t522, qJD(3) * t54 + t614, t31 * qJD(3) + t39 * qJD(4) + t612 + t805; qJD(2) * t7 + qJD(4) * t1 + qJD(5) * t16 + qJD(6) * t6 - t639, -qJD(4) * t17 + qJD(5) * t55 + qJD(6) * t32 - t523, qJD(4) * t12 + qJD(5) * t37 + qJD(6) * t23 (m(7) * (t156 * t405 + t406 * t788) + (-t353 * t616 + t576 * t615) * t738 - t576 * mrSges(6,2) - t353 * mrSges(6,1) - t326 * t664 - t532 * t665 + t756 * t530 - t508 * t529 + t536 + t749 * pkin(8) + t751 + t811) * qJD(4) + t79 * qJD(5) + t813 + t511, qJD(4) * t79 + qJD(6) * t84 - t509, t24 * qJD(4) + t84 * qJD(5) + t510 + t813; qJD(2) * t10 - qJD(3) * t1 + qJD(5) * t45 + qJD(6) * t14 - t628, qJD(3) * t17 + qJD(6) * t38 + t522, qJD(5) * t47 + qJD(6) * t25 - t511, t573, t520, t506 + t573; -qJD(3) * t16 - qJD(4) * t45 - qJD(6) * t72 - t521, -qJD(3) * t55 - t614, -qJD(4) * t47 + qJD(6) * t83 + t509, -t520, 0, -t519; -qJD(2) * t19 - qJD(3) * t6 - qJD(4) * t14 + qJD(5) * t72 - t618, -t32 * qJD(3) - t38 * qJD(4) - t612, -qJD(4) * t25 - qJD(5) * t83 - t510, -t506, t519, 0;];
Cq  = t21;
