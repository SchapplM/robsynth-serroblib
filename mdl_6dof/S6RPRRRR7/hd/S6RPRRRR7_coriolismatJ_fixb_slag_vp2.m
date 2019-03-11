% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:17:03
% EndTime: 2019-03-09 07:17:26
% DurationCPUTime: 15.31s
% Computational Cost: add. (35880->602), mult. (63533->759), div. (0->0), fcn. (72534->8), ass. (0->368)
t397 = sin(qJ(6));
t399 = cos(qJ(6));
t663 = sin(qJ(4));
t664 = sin(qJ(3));
t666 = cos(qJ(4));
t667 = cos(qJ(3));
t368 = t663 * t664 - t666 * t667;
t398 = sin(qJ(5));
t563 = t663 * t667 + t666 * t664;
t665 = cos(qJ(5));
t311 = t665 * t368 + t398 * t563;
t382 = t664 * pkin(3) + qJ(2);
t334 = t563 * pkin(4) + t382;
t314 = -t398 * t368 + t665 * t563;
t658 = pkin(5) * t314;
t178 = pkin(10) * t311 + t334 + t658;
t682 = pkin(1) + pkin(7);
t496 = t664 * t682;
t444 = -t664 * pkin(8) - t496;
t497 = t667 * t682;
t445 = -t667 * pkin(8) - t497;
t324 = t666 * t444 + t663 * t445;
t414 = t563 * pkin(9) - t324;
t698 = -t663 * t444 + t666 * t445;
t747 = t368 * pkin(9) + t698;
t803 = t398 * t747 - t665 * t414;
t73 = t178 * t399 - t397 * t803;
t74 = t397 * t178 + t399 * t803;
t482 = t397 * t73 - t399 * t74;
t828 = t803 + t482;
t781 = t314 * mrSges(6,1);
t782 = t311 * mrSges(6,2);
t222 = t781 - t782;
t662 = m(6) * t334;
t827 = t222 + t662;
t623 = t399 * mrSges(7,3);
t736 = -mrSges(7,1) * t311 + t314 * t623;
t624 = t399 * mrSges(7,2);
t629 = t397 * mrSges(7,1);
t485 = t624 + t629;
t773 = t311 * t485;
t826 = t73 * t736 - t803 * t773;
t308 = Ifges(6,5) * t314;
t392 = Ifges(7,4) * t399;
t374 = Ifges(7,1) * t397 + t392;
t567 = t399 * t374;
t655 = Ifges(7,4) * t397;
t372 = Ifges(7,2) * t399 + t655;
t578 = t397 * t372;
t462 = t578 / 0.2e1 - t567 / 0.2e1;
t371 = Ifges(7,5) * t397 + Ifges(7,6) * t399;
t719 = t311 * t371;
t723 = Ifges(6,6) * t311;
t714 = -Ifges(7,2) * t397 + t392;
t735 = -Ifges(7,6) * t311 - t314 * t714;
t765 = t399 * t735;
t375 = Ifges(7,1) * t399 - t655;
t737 = -Ifges(7,5) * t311 - t314 * t375;
t767 = t397 * t737;
t785 = t723 - t719 / 0.2e1 + t765 / 0.2e1 + t767 / 0.2e1;
t441 = -Ifges(5,5) * t563 + Ifges(5,6) * t368 + t314 * t462 - t308 + t785;
t752 = t698 * mrSges(5,2);
t780 = t324 * mrSges(5,1);
t625 = t399 * mrSges(7,1);
t628 = t397 * mrSges(7,2);
t486 = t625 - t628;
t813 = t803 * t486;
t804 = t398 * t414 + t665 * t747;
t817 = t804 * mrSges(6,2);
t818 = t803 * mrSges(6,1);
t807 = -t813 - t818 - t817;
t825 = t441 - t752 - t780 + t807;
t579 = t397 * t314;
t725 = mrSges(7,2) * t311;
t212 = mrSges(7,3) * t579 + t725;
t683 = Ifges(7,3) / 0.2e1;
t554 = t683 + Ifges(6,2) / 0.2e1;
t676 = t314 / 0.2e1;
t391 = Ifges(7,5) * t399;
t652 = Ifges(7,6) * t397;
t715 = t391 - t652;
t753 = Ifges(7,5) * t314;
t140 = -t311 * t375 + t753;
t573 = t399 * t140;
t137 = Ifges(7,6) * t314 - t311 * t714;
t585 = t397 * t137;
t730 = -t311 / 0.2e1;
t738 = -t334 * mrSges(6,2) - Ifges(6,1) * t730 + Ifges(6,4) * t314 + t585 / 0.2e1 - t573 / 0.2e1;
t642 = t314 * mrSges(6,3);
t751 = t485 * t314;
t750 = -t751 - t642;
t754 = -t314 / 0.2e1;
t668 = t399 / 0.2e1;
t669 = -t397 / 0.2e1;
t756 = t668 * t737 + t669 * t735 + t334 * mrSges(6,1) + t715 * t730 + Ifges(6,4) * t311 + (Ifges(6,2) + Ifges(7,3)) * t676;
t799 = t804 * t314;
t824 = -t799 * mrSges(6,3) + t382 * (-t368 * mrSges(5,1) - t563 * mrSges(5,2)) + t74 * t212 + t563 ^ 2 * Ifges(5,4) - (t676 * t715 - t738) * t314 - (Ifges(6,1) * t754 + t314 * t554 + t756) * t311 - t750 * t804 + t826;
t689 = m(7) / 0.2e1;
t822 = t311 * t828 - t799;
t688 = m(6) * pkin(4);
t797 = t804 * t398;
t821 = (-t665 * t803 + t797) * t688;
t820 = -t813 / 0.2e1 - t817 / 0.2e1 - t818 / 0.2e1;
t819 = pkin(5) * t803;
t557 = t666 * pkin(3);
t388 = t557 + pkin(4);
t534 = t663 * t398;
t350 = -pkin(3) * t534 + t665 * t388;
t344 = -pkin(5) - t350;
t770 = t344 * t803;
t556 = t665 * pkin(4);
t387 = -t556 - pkin(5);
t815 = t387 * t803;
t795 = t397 * t804;
t494 = t665 * t663;
t351 = pkin(3) * t494 + t398 * t388;
t798 = t804 * t351;
t796 = t804 * t399;
t812 = t804 * t803;
t757 = -t578 / 0.4e1 + t567 / 0.4e1;
t810 = -t757 * t314 - t719 / 0.4e1 + t723 / 0.2e1 + t765 / 0.4e1 + t767 / 0.4e1;
t571 = t399 * t212;
t768 = t397 * t736;
t802 = -t768 / 0.2e1;
t710 = t571 / 0.2e1 + t802;
t755 = -pkin(5) / 0.2e1;
t437 = t710 * pkin(10) - t308 / 0.2e1 - t751 * t755 + t810;
t225 = -pkin(5) * t311 + pkin(10) * t314;
t661 = pkin(4) * t368;
t180 = t225 - t661;
t79 = t180 * t399 - t795;
t80 = t397 * t180 + t796;
t480 = -t79 * t397 + t80 * t399;
t809 = t437 + (t480 * pkin(10) - t819) * t689;
t345 = pkin(10) + t351;
t395 = t397 ^ 2;
t396 = t399 ^ 2;
t562 = t395 + t396;
t743 = t562 * t345;
t808 = t351 - t743;
t784 = -t399 / 0.2e1;
t498 = -t375 * t669 - t714 * t784 - t462;
t263 = t773 / 0.2e1;
t467 = t624 / 0.2e1 + t629 / 0.2e1;
t454 = t467 * t311;
t794 = t263 + t454;
t806 = t794 * qJD(6);
t771 = t314 * t486;
t697 = t771 / 0.2e1 + t781 / 0.2e1 - t782 / 0.2e1;
t630 = t396 * mrSges(7,3);
t631 = t395 * mrSges(7,3);
t744 = -t631 / 0.2e1 - t630 / 0.2e1;
t805 = -t744 * t311 + t697;
t801 = -t773 / 0.2e1;
t790 = t752 / 0.2e1 + t780 / 0.2e1;
t764 = t562 * t311;
t788 = -pkin(10) * t764 - t658;
t787 = -t563 * mrSges(5,1) + t368 * mrSges(5,2);
t783 = m(5) * t382;
t775 = t311 * t314;
t774 = t311 * t398;
t769 = t387 * t314;
t592 = t311 * t351;
t593 = t350 * t314;
t758 = t593 + t592;
t507 = t562 * t350;
t749 = t344 + t507;
t356 = (t665 * t666 - t534) * pkin(3);
t216 = mrSges(7,1) * t314 + t311 * t623;
t627 = t397 * mrSges(7,3);
t213 = -mrSges(7,2) * t314 + t311 * t627;
t570 = t399 * t213;
t461 = t570 / 0.2e1 + t216 * t669;
t540 = t623 / 0.2e1;
t543 = -t627 / 0.2e1;
t489 = t80 * t540 + t79 * t543 + t820;
t355 = (t666 * t398 + t494) * pkin(3);
t590 = t355 * t804;
t594 = t344 * t751;
t673 = t355 / 0.2e1;
t691 = m(6) / 0.2e1;
t742 = (t798 - t590 + (-t350 + t356) * t803) * t691 + (t480 * t345 - t482 * t356 - t590 + t770) * t689 - t594 / 0.2e1 - t773 * t673 + t489 + t710 * t345 + (t356 * t754 - t311 * t673 + t592 / 0.2e1 + t593 / 0.2e1) * mrSges(6,3) + t461 * t356 - t790;
t731 = t344 / 0.2e1;
t728 = m(7) * t387;
t721 = (t396 / 0.2e1 + t395 / 0.2e1) * mrSges(7,3);
t701 = t562 * mrSges(7,3);
t205 = t486 * t311;
t660 = pkin(4) * t398;
t386 = pkin(10) + t660;
t433 = (-t374 / 0.4e1 - t714 / 0.4e1) * t397 + (t375 / 0.4e1 - t372 / 0.4e1) * t399;
t671 = t387 / 0.2e1;
t417 = -(-t386 * t721 + t433) * t311 - t205 * t671;
t210 = t374 * t311;
t679 = t210 / 0.4e1;
t681 = -t804 / 0.2e1;
t429 = (-t314 / 0.4e1 + t754) * Ifges(7,6) - t137 / 0.4e1 + t679 + mrSges(7,1) * t681;
t209 = t311 * t372;
t680 = t209 / 0.4e1;
t440 = t140 / 0.4e1 + t680 + t753 / 0.2e1 + mrSges(7,2) * t681;
t651 = Ifges(7,3) * t311;
t670 = t391 / 0.4e1;
t464 = t651 / 0.2e1 + t314 * t670;
t685 = mrSges(7,2) / 0.2e1;
t686 = -mrSges(7,1) / 0.2e1;
t470 = t80 * t685 + t79 * t686;
t672 = -t386 / 0.2e1;
t20 = (t213 * t672 + t429) * t397 + (t216 * t672 + t440) * t399 + t417 + t464 + t470;
t458 = t387 * t485;
t220 = t458 + t498;
t106 = t801 + t454;
t561 = t106 * qJD(2);
t420 = -t458 / 0.2e1 - t498;
t459 = t344 * t485;
t451 = -t459 / 0.2e1;
t452 = t467 * t356;
t96 = t451 - t452 + t420;
t713 = t20 * qJD(1) - t96 * qJD(3) + t220 * qJD(4) - t561;
t418 = -(-t345 * t721 + t433) * t311 - t205 * t731;
t394 = t667 * pkin(3);
t179 = t180 + t394;
t75 = t179 * t399 - t795;
t76 = t397 * t179 + t796;
t471 = t76 * t685 + t75 * t686;
t674 = -t345 / 0.2e1;
t18 = (t213 * t674 + t429) * t397 + (t216 * t674 + t440) * t399 + t418 + t464 + t471;
t190 = t459 + t498;
t712 = t18 * qJD(1) + t190 * qJD(3) - t561;
t711 = t663 * pkin(3) * t368 + t563 * t557;
t211 = t314 * t627 + t725;
t572 = t399 * t211;
t708 = t572 / 0.2e1 + t802;
t569 = t399 * t216;
t582 = t397 * t213;
t463 = -t582 / 0.2e1 - t569 / 0.2e1;
t706 = t571 - t768;
t481 = -t75 * t397 + t76 * t399;
t87 = t225 * t399 - t795;
t88 = t397 * t225 + t796;
t479 = -t397 * t87 + t399 * t88;
t294 = t311 * t631;
t295 = t311 * t630;
t704 = -t771 - t294 - t295;
t702 = (-mrSges(5,1) * t663 - mrSges(5,2) * t666) * pkin(3);
t700 = (-t391 / 0.2e1 + t652 / 0.2e1) * t314;
t690 = -m(7) / 0.2e1;
t687 = m(7) * pkin(4);
t684 = -mrSges(7,3) / 0.2e1;
t644 = t311 * mrSges(6,3);
t635 = t350 * mrSges(6,2);
t634 = t351 * mrSges(6,1);
t633 = t355 * mrSges(6,1);
t632 = t356 * mrSges(6,2);
t53 = m(7) * (-t314 * t764 + t775);
t566 = t53 * qJD(2);
t617 = -t106 * qJD(6) - t566;
t616 = t566 + t806;
t457 = -t222 + t787;
t425 = -t664 * mrSges(4,1) - t667 * mrSges(4,2) + t457;
t37 = t582 + t569 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(7) * (t74 * t397 + t399 * t73) + t662 + t783 - t425;
t615 = qJD(1) * t37;
t405 = -t773 * t676 + (t644 / 0.2e1 + t710) * t314 + (t750 / 0.2e1 - t461) * t311;
t11 = (t481 * t314 + t822) * t689 + t405;
t614 = t11 * qJD(1);
t13 = (t480 * t314 + t822) * t689 + t405;
t613 = t13 * qJD(1);
t14 = t804 * t205 + t73 * t213 - t74 * t216 - (t371 * t754 + t210 * t668 + t137 * t784 + t482 * mrSges(7,3) + (t140 + t209) * t669) * t311;
t612 = t14 * qJD(1);
t601 = t311 * t205;
t599 = t311 * t355;
t468 = -t628 / 0.2e1 + t625 / 0.2e1;
t32 = t601 / 0.2e1 + (-t311 * t721 - t463) * t314 + t468;
t595 = t32 * qJD(1);
t591 = t351 * t486;
t589 = t355 * t486;
t588 = t387 * t751;
t560 = qJD(3) + qJD(4);
t559 = mrSges(6,1) * t660;
t16 = (t801 + (-t804 + t479) * t689 + t708) * t314 + (-t751 / 0.2e1 + t828 * t689 - t461) * t311;
t558 = t11 * qJD(3) + t13 * qJD(4) + t16 * qJD(5);
t555 = t660 / 0.2e1;
t550 = t391 / 0.2e1;
t542 = t627 / 0.2e1;
t541 = -t623 / 0.2e1;
t539 = t486 / 0.2e1 + mrSges(6,1) / 0.2e1;
t538 = -t311 * t743 + t314 * t344;
t537 = t397 * t665;
t536 = t399 * t665;
t506 = t562 * t356;
t505 = t562 * t386;
t504 = t644 * t660;
t502 = mrSges(6,2) * t556;
t500 = t556 / 0.2e1;
t493 = -t537 / 0.2e1;
t490 = t556 * t642;
t488 = t76 * t540 + t75 * t543 + t820;
t466 = pkin(5) * t485;
t478 = -t466 / 0.2e1 + t498;
t434 = -Ifges(5,4) * t368 + (Ifges(5,1) - Ifges(5,2)) * t563;
t1 = qJ(2) * (t667 * mrSges(4,1) - t664 * mrSges(4,2)) + m(7) * (t73 * t75 + t74 * t76 - t812) + t434 * t368 + (t664 ^ 2 - t667 ^ 2) * Ifges(4,4) + t75 * t216 + t76 * t213 + (Ifges(4,2) - Ifges(4,1)) * t667 * t664 + (t783 - t787) * t394 + t827 * (t394 - t661) + t824;
t477 = t1 * qJD(1) + t11 * qJD(2);
t3 = m(7) * (t73 * t79 + t74 * t80 - t812) + t79 * t216 + t80 * t213 + (-pkin(4) * t827 + t434) * t368 + t824;
t476 = t3 * qJD(1) + t13 * qJD(2);
t9 = t804 * t751 + m(7) * (t73 * t87 + t74 * t88 - t812) + t88 * t213 + t74 * t211 + t87 * t216 - t756 * t311 + (t700 - (-Ifges(6,1) / 0.2e1 + t554) * t311 + t738) * t314 + t826;
t475 = t9 * qJD(1) + t16 * qJD(2);
t55 = m(7) * (0.1e1 - t562) * t775;
t474 = t16 * qJD(1) + t55 * qJD(2);
t469 = t88 * t685 + t87 * t686;
t465 = t314 * t715 / 0.4e1;
t456 = t466 / 0.2e1;
t455 = t562 * t665;
t453 = t467 * t350;
t423 = -t294 / 0.2e1 - t295 / 0.2e1 + t788 * t689 - t697;
t28 = (t749 * t690 + t539) * t314 + (t808 * t690 - mrSges(6,2) / 0.2e1 + t721) * t311 + t423;
t409 = (t479 * t345 - t482 * t350 + t770 - t798) * t690 + t751 * t731 + t351 * t263;
t411 = t437 + (t481 * pkin(10) - t819) * t689 + t488;
t4 = t411 + t409 + (t88 * t684 - t350 * t213 / 0.2e1 + t211 * t674 - t735 / 0.4e1) * t399 + (t350 * t216 / 0.2e1 + t345 * t736 / 0.2e1 + t87 * mrSges(7,3) / 0.2e1 - t737 / 0.4e1) * t397 + t539 * t803 + t817 / 0.2e1 + (Ifges(6,5) / 0.2e1 + t757) * t314 - (-t371 / 0.4e1 + Ifges(6,6) / 0.2e1) * t311;
t426 = mrSges(7,3) * t507 - t591 - t634 - t635;
t58 = m(7) * (t344 * t351 + t345 * t507) + t426;
t450 = -t4 * qJD(1) - t28 * qJD(2) + t58 * qJD(3);
t406 = (t314 * t356 + t599 - t758) * t691 + (t314 * t506 + t538 + t599) * t689;
t442 = m(7) * (-t311 * t505 + t769);
t443 = (-t314 * t665 - t774) * t688;
t422 = -t442 / 0.2e1 - t443 / 0.2e1;
t36 = t406 + t422;
t59 = (-t486 - mrSges(6,1)) * t355 + t702 + (-mrSges(6,2) + t701) * t356 + m(7) * (t344 * t355 + t345 * t506) + m(6) * (-t350 * t355 + t351 * t356);
t404 = t690 * t815 + t588 / 0.2e1 - t821 / 0.2e1 - t504 / 0.2e1 - t490 / 0.2e1 + (t481 * t690 - t710) * t386 + t790;
t6 = t76 * t541 + t75 * t542 + t404 + t742 - t820;
t449 = t6 * qJD(1) + t36 * qJD(2) + t59 * qJD(3);
t438 = t457 + t704;
t419 = -t486 * t660 + t556 * t701 - t502 - t559;
t191 = (t455 * t386 + t387 * t398) * t687 + t419;
t416 = m(7) * (t769 - t386 * t764 + (t455 * t314 + t774) * pkin(4));
t34 = -t416 / 0.2e1 + t423 + t805;
t400 = (t387 * t351 + (t344 * t398 + t455 * t345) * pkin(4)) * t689 - t635 / 0.2e1 - t634 / 0.2e1 - t591 / 0.2e1 - t559 / 0.2e1 - t486 * t555 - t502 / 0.2e1 + t500 * t701 + (t505 * t689 - t744) * t350;
t410 = (-pkin(5) * t355 + pkin(10) * t506) * t690 + t633 / 0.2e1 + t589 / 0.2e1 + t632 / 0.2e1 + t744 * t356;
t42 = t400 + t410;
t427 = Ifges(6,5) * t754 + t88 * t540 + t87 * t543 + t810 + t820;
t401 = -t773 * t555 + t427 - t751 * t671 + pkin(4) * t216 * t493 + t500 * t570 + (t815 + (t74 * t536 - t73 * t537 - t797) * pkin(4)) * t689 + (t479 * t689 + t708) * t386;
t8 = t80 * t541 + t79 * t542 + t401 - t809 - t820;
t435 = t8 * qJD(1) - t34 * qJD(2) + t42 * qJD(3) + t191 * qJD(4);
t428 = (mrSges(7,1) * t493 - mrSges(7,2) * t536 / 0.2e1) * pkin(4);
t149 = t456 + t428 + t420;
t421 = t485 * t681 - t585 / 0.4e1 + t397 * t679 + t573 / 0.4e1 + t399 * t680;
t407 = t463 * pkin(10) - t205 * t755 + t421;
t424 = -pkin(10) * t721 + t433;
t22 = (-0.3e1 / 0.4e1 * t652 + t670 + t550) * t314 - (-Ifges(7,3) / 0.2e1 + t424) * t311 + t407 + t469;
t230 = -t466 + t498;
t98 = t451 + t456 - t453 - t498;
t432 = t22 * qJD(1) - t98 * qJD(3) - t149 * qJD(4) + t230 * qJD(5);
t413 = t423 - t805;
t408 = Ifges(7,6) * t579 / 0.2e1 - t314 * t550 - t651 / 0.2e1 + t465 + t421;
t348 = t458 / 0.2e1;
t321 = t459 / 0.2e1;
t150 = t348 + t428 + t478;
t99 = t321 - t453 + t478;
t97 = t321 + t348 - t452 + t498;
t41 = t400 - t410;
t35 = t416 / 0.2e1 + t413;
t33 = -t601 / 0.2e1 + t468 - t562 * t775 * t684 + t463 * t314;
t30 = t406 - t422 + t438;
t29 = (t808 * t311 + t749 * t314) * t689 + t413;
t21 = t407 + t465 - t469 - (t424 + t683) * t311 + t700;
t19 = t463 * t386 + t408 + t417 - t470;
t17 = t463 * t345 + t408 + t418 - t471;
t7 = t401 + t489 + t809;
t5 = t345 * t708 + t461 * t350 - t409 + t411 + t427;
t2 = -t404 + t441 + t488 + t742;
t10 = [qJD(2) * t37 + qJD(3) * t1 + qJD(4) * t3 + qJD(5) * t9 + qJD(6) * t14, qJD(6) * t33 + t558 + t615 (m(5) * (-t324 * t666 + t663 * t698) * pkin(3) + mrSges(4,1) * t496 + mrSges(4,2) * t497 - Ifges(4,5) * t664 - Ifges(4,6) * t667 + m(6) * (-t350 * t803 + t798) - t594 + m(7) * t770 + (m(7) * t481 + t706) * t345 + t481 * mrSges(7,3) + t758 * mrSges(6,3) + t711 * mrSges(5,3) + t825) * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t17 * qJD(6) + t477, t2 * qJD(3) + (t821 - t588 + t490 + t504 + t803 * t728 + (m(7) * t480 + t706) * t386 + t480 * mrSges(7,3) + t825) * qJD(4) + t7 * qJD(5) + t19 * qJD(6) + t476, t5 * qJD(3) + t7 * qJD(4) + t21 * qJD(6) + t475 + ((-Ifges(6,5) + t462) * t314 + (-m(7) * t803 + t751) * pkin(5) + (m(7) * t479 + t572 - t768) * pkin(10) + t479 * mrSges(7,3) + t785 + t807) * qJD(5), t612 + t33 * qJD(2) + t17 * qJD(3) + t19 * qJD(4) + t21 * qJD(5) + (-t74 * mrSges(7,1) - t73 * mrSges(7,2) + t719) * qJD(6); -qJD(6) * t32 + t558 - t615, qJD(5) * t55 + t560 * t53, t30 * qJD(4) + t29 * qJD(5) + t614 + t616 + (-m(5) * t711 + 0.2e1 * t538 * t689 - 0.2e1 * t691 * t758 + t425 + t704) * qJD(3), t613 + t30 * qJD(3) + (t442 + t443 + t438) * qJD(4) + t35 * qJD(5) + t616, t29 * qJD(3) + t35 * qJD(4) + (m(7) * t788 - mrSges(7,3) * t764 - t222 - t771) * qJD(5) + t806 + t474, -qJD(6) * t771 - t595 + (qJD(5) + t560) * t794; qJD(4) * t6 - qJD(5) * t4 + qJD(6) * t18 - t477, qJD(4) * t36 - qJD(5) * t28 - t614 + t617, qJD(4) * t59 + qJD(5) * t58 + qJD(6) * t190 (-t589 - t632 - t633 + (-t665 * t688 + t728) * t355 + (m(7) * t505 + t398 * t688 + t630 + t631) * t356 + t702) * qJD(4) + t41 * qJD(5) + t97 * qJD(6) + t449, t41 * qJD(4) + (m(7) * (-pkin(5) * t351 + pkin(10) * t507) + t426) * qJD(5) + t99 * qJD(6) + t450, t97 * qJD(4) + t99 * qJD(5) + (-t345 * t486 + t715) * qJD(6) + t712; -qJD(3) * t6 + qJD(5) * t8 + qJD(6) * t20 - t476, -qJD(3) * t36 - qJD(5) * t34 - t613 + t617, qJD(5) * t42 - qJD(6) * t96 - t449, qJD(5) * t191 + qJD(6) * t220 ((-pkin(5) * t398 + pkin(10) * t455) * t687 + t419) * qJD(5) + t150 * qJD(6) + t435, t150 * qJD(5) + (-t386 * t486 + t715) * qJD(6) + t713; qJD(3) * t4 - qJD(4) * t8 + qJD(6) * t22 - t475, qJD(3) * t28 + qJD(4) * t34 - t474, -qJD(4) * t42 - qJD(6) * t98 - t450, -qJD(6) * t149 - t435, t230 * qJD(6) (-pkin(10) * t486 + t715) * qJD(6) + t432; qJD(2) * t32 - qJD(3) * t18 - qJD(4) * t20 - qJD(5) * t22 - t612, t560 * t106 + t595, qJD(4) * t96 + qJD(5) * t98 - t712, qJD(5) * t149 - t713, -t432, 0;];
Cq  = t10;
