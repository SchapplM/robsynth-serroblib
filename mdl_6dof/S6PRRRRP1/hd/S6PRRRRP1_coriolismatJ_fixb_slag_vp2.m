% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRRP1
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
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:55:58
% EndTime: 2019-03-08 23:56:25
% DurationCPUTime: 16.49s
% Computational Cost: add. (20460->732), mult. (46505->964), div. (0->0), fcn. (49870->10), ass. (0->415)
t463 = cos(qJ(3));
t461 = sin(qJ(3));
t724 = sin(qJ(4));
t597 = t724 * t461;
t726 = cos(qJ(4));
t413 = -t463 * t726 + t597;
t459 = sin(pkin(6));
t464 = cos(qJ(2));
t646 = t459 * t464;
t362 = t413 * t646;
t460 = sin(qJ(5));
t462 = cos(qJ(5));
t725 = sin(qJ(2));
t603 = t459 * t725;
t291 = t460 * t362 + t462 * t603;
t292 = -t462 * t362 + t460 * t603;
t598 = t726 * t461;
t414 = -t463 * t724 - t598;
t361 = t414 * t646;
t627 = t724 * pkin(3);
t440 = t627 + pkin(10);
t641 = qJ(6) + t440;
t398 = t641 * t460;
t399 = t641 * t462;
t628 = t726 * pkin(3);
t441 = -t628 - pkin(4);
t714 = t462 * pkin(5);
t420 = t441 - t714;
t667 = t292 * t462;
t668 = t291 * t460;
t524 = t667 - t668;
t772 = m(5) * pkin(3);
t773 = m(7) / 0.2e1;
t775 = m(6) / 0.2e1;
t808 = mrSges(6,3) + mrSges(7,3);
t873 = -mrSges(5,2) / 0.2e1;
t844 = -(t667 / 0.2e1 - t668 / 0.2e1) * t808 + t362 * t873;
t423 = -mrSges(7,1) * t462 + mrSges(7,2) * t460;
t858 = -t423 / 0.2e1;
t708 = mrSges(6,1) * t462;
t424 = mrSges(6,2) * t460 - t708;
t861 = mrSges(5,1) / 0.2e1 - t424 / 0.2e1;
t778 = (t858 + t861) * t361 - t844;
t427 = mrSges(4,1) * t461 + mrSges(4,2) * t463;
t580 = -t646 / 0.2e1;
t804 = t427 * t580;
t874 = -(-t361 * t441 + t440 * t524) * t775 - (-t291 * t398 + t292 * t399 - t361 * t420) * t773 - (t361 * t726 - t362 * t724) * t772 / 0.2e1 - t804 - t778;
t685 = cos(pkin(6));
t390 = t461 * t685 + t463 * t603;
t498 = t461 * t603 - t463 * t685;
t277 = t390 * t724 + t726 * t498;
t426 = t460 * mrSges(6,1) + t462 * mrSges(6,2);
t769 = m(7) * pkin(5);
t632 = t769 / 0.2e1;
t766 = mrSges(7,1) / 0.2e1;
t555 = t766 + t632;
t767 = mrSges(6,1) / 0.2e1;
t520 = t767 + t555;
t764 = mrSges(7,2) / 0.2e1;
t765 = mrSges(6,2) / 0.2e1;
t620 = t765 + t764;
t551 = t620 * t462;
t483 = (t460 * t520 + t551) * t277;
t716 = pkin(5) * t460;
t621 = t716 / 0.2e1;
t446 = t460 * mrSges(7,1);
t448 = t462 * mrSges(7,2);
t639 = t448 + t446;
t758 = t277 / 0.2e1;
t865 = t483 + m(7) * t277 * t621 + (t426 + t639) * t758;
t480 = t390 * t726 - t498 * t724;
t186 = -t460 * t480 - t462 * t646;
t187 = -t460 * t646 + t462 * t480;
t79 = -t186 * t460 + t187 * t462;
t814 = m(7) + m(6);
t859 = (t480 - t79) * t277 * t814;
t868 = t859 * qJD(1);
t872 = t865 * qJD(5) + t868;
t871 = qJD(3) + qJD(4);
t455 = t460 ^ 2;
t457 = t462 ^ 2;
t638 = t455 + t457;
t788 = -t448 / 0.2e1 - t446 / 0.2e1;
t655 = t413 * t460;
t761 = -pkin(9) - pkin(8);
t434 = t761 * t463;
t797 = -t726 * t434 + t761 * t597;
t823 = -pkin(5) * t655 + t797;
t870 = t788 * t413 + t823 * t773;
t443 = -pkin(3) * t463 - pkin(2);
t715 = pkin(10) * t414;
t310 = pkin(4) * t413 + t443 + t715;
t141 = t462 * t310 - t460 * t797;
t142 = t310 * t460 + t462 * t797;
t869 = t141 * t460 - t142 * t462 + t797;
t653 = t414 * t460;
t318 = -mrSges(7,2) * t413 + mrSges(7,3) * t653;
t652 = t414 * t462;
t322 = mrSges(7,1) * t413 + mrSges(7,3) * t652;
t728 = t462 / 0.2e1;
t732 = -t460 / 0.2e1;
t799 = t318 * t728 + t322 * t732;
t866 = -t799 + t870;
t864 = t799 + t870;
t701 = Ifges(7,4) * t462;
t538 = -Ifges(7,2) * t460 + t701;
t195 = -Ifges(7,6) * t414 - t413 * t538;
t703 = Ifges(6,4) * t462;
t540 = -Ifges(6,2) * t460 + t703;
t197 = -Ifges(6,6) * t414 - t413 * t540;
t702 = Ifges(7,4) * t460;
t542 = Ifges(7,1) * t462 - t702;
t199 = -Ifges(7,5) * t414 - t413 * t542;
t704 = Ifges(6,4) * t460;
t544 = Ifges(6,1) * t462 - t704;
t201 = -Ifges(6,5) * t414 - t413 * t544;
t541 = Ifges(7,1) * t460 + t701;
t513 = t462 * t541;
t543 = Ifges(6,1) * t460 + t703;
t514 = t462 * t543;
t537 = Ifges(7,2) * t462 + t702;
t515 = t460 * t537;
t539 = Ifges(6,2) * t462 + t704;
t516 = t460 * t539;
t730 = t460 / 0.2e1;
t737 = -t414 / 0.2e1;
t838 = Ifges(7,5) + Ifges(6,5);
t712 = Ifges(6,6) + Ifges(7,6);
t852 = t712 * t462;
t489 = -Ifges(5,5) * t413 + Ifges(5,6) * t414 + (t515 + t516) * t413 / 0.2e1 - (t513 + t514) * t413 / 0.2e1 + (t838 * t460 + t852) * t737 + (t199 + t201) * t730 + (t195 + t197) * t728;
t826 = t797 * t424;
t835 = t797 * mrSges(5,1);
t352 = -t434 * t724 - t761 * t598;
t836 = t352 * mrSges(5,2);
t851 = t823 * t423;
t863 = t489 + t826 + t836 - t835 + t851;
t862 = -t826 / 0.2e1 + t835 / 0.2e1 - t836 / 0.2e1 - t851 / 0.2e1;
t253 = -pkin(5) * t653 + t352;
t857 = t253 * t823;
t856 = t420 * t823;
t827 = t638 * t277;
t855 = t440 * t827;
t442 = -pkin(4) - t714;
t854 = t442 * t823;
t853 = t638 * t808;
t324 = -mrSges(5,1) * t414 - mrSges(5,2) * t413;
t304 = t639 * t414;
t305 = t426 * t414;
t316 = mrSges(7,2) * t414 + mrSges(7,3) * t655;
t317 = mrSges(6,2) * t414 + mrSges(6,3) * t655;
t654 = t413 * t462;
t320 = -mrSges(7,1) * t414 + mrSges(7,3) * t654;
t321 = -mrSges(6,1) * t414 + mrSges(6,3) * t654;
t781 = (t316 / 0.2e1 + t317 / 0.2e1) * t187 + (t320 / 0.2e1 + t321 / 0.2e1) * t186 + (-t304 / 0.2e1 - t305 / 0.2e1) * t480;
t849 = t324 * t580 + t781;
t507 = t638 * t726;
t319 = -mrSges(6,2) * t413 + mrSges(6,3) * t653;
t746 = t319 / 0.2e1;
t747 = t318 / 0.2e1;
t570 = t746 + t747;
t302 = t639 * t413;
t303 = t426 * t413;
t572 = -t302 / 0.2e1 - t303 / 0.2e1;
t323 = mrSges(6,1) * t413 + mrSges(6,3) * t652;
t740 = t323 / 0.2e1;
t742 = t322 / 0.2e1;
t847 = t460 * (t740 + t742) - t462 * t570 + t572;
t106 = qJ(6) * t653 + t142;
t105 = qJ(6) * t652 + t141;
t85 = pkin(5) * t413 + t105;
t532 = t106 * t462 - t460 * t85;
t846 = t823 - t532;
t196 = t413 * Ifges(7,6) - t414 * t538;
t198 = t413 * Ifges(6,6) - t414 * t540;
t617 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t501 = t198 / 0.2e1 + t196 / 0.2e1 + t617 * t413;
t618 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t200 = t413 * Ifges(7,5) - t414 * t542;
t202 = t413 * Ifges(6,5) - t414 * t544;
t782 = -t618 * t413 - t200 / 0.2e1 - t202 / 0.2e1;
t837 = Ifges(7,3) + Ifges(6,3);
t845 = -t823 * t304 - t797 * t305 + (Ifges(5,4) * t413 + t501 * t460 + t462 * t782) * t413 + (-Ifges(5,4) * t414 + (-t201 / 0.2e1 - t199 / 0.2e1 + t618 * t414) * t462 + (t197 / 0.2e1 + t195 / 0.2e1 - t617 * t414) * t460 + (-Ifges(5,2) + Ifges(5,1) - t837) * t413) * t414 + t106 * t316 + t141 * t321 + t142 * t317 - t253 * t302 - t352 * t303 + t85 * t320 + t443 * t324;
t711 = -qJ(6) - pkin(10);
t422 = t711 * t460;
t425 = t711 * t462;
t522 = t422 * t460 + t425 * t462;
t477 = (-pkin(4) * t480 - pkin(10) * t827) * t775 + (t277 * t522 + t442 * t480) * t773;
t841 = -t426 / 0.2e1;
t840 = m(5) * t443;
t839 = pkin(4) * t797;
t692 = t413 * mrSges(5,3);
t656 = t399 * t462;
t658 = t398 * t460;
t523 = -t656 - t658;
t834 = t277 * t523;
t833 = t277 * t724;
t831 = t352 * t460;
t805 = t352 * t480;
t830 = t352 * t724;
t829 = t441 * t797;
t828 = t462 * t352;
t666 = t797 * t352;
t825 = t318 + t319;
t824 = t322 + t323;
t822 = t461 ^ 2 + t463 ^ 2;
t821 = t538 + t540;
t820 = t544 + t542;
t723 = m(7) * t253;
t809 = mrSges(6,2) + mrSges(7,2);
t807 = -mrSges(7,1) - t769;
t806 = t253 * t480;
t713 = Ifges(7,4) + Ifges(6,4);
t800 = t713 * t460;
t798 = -t304 - t305;
t741 = -t323 / 0.2e1;
t794 = t319 * t732 + t462 * t741;
t789 = -mrSges(4,1) * t463 + mrSges(4,2) * t461;
t786 = mrSges(6,1) - t807;
t785 = t424 + t423 - mrSges(5,1);
t736 = t423 / 0.2e1;
t549 = t736 - t861;
t784 = -t253 * t639 / 0.2e1 + t352 * t841;
t300 = -mrSges(7,1) * t652 + mrSges(7,2) * t653;
t630 = pkin(5) * t652;
t566 = m(7) * t630;
t185 = -t300 + t566;
t634 = m(7) * t716;
t397 = -t634 - t639;
t783 = qJD(2) * t185 + t397 * t871;
t326 = -pkin(4) * t414 + pkin(10) * t413;
t162 = t460 * t326 - t828;
t616 = qJ(6) * t655;
t111 = t616 + t162;
t161 = t462 * t326 + t831;
t548 = -t414 * pkin(5) + qJ(6) * t654;
t89 = t161 + t548;
t779 = (t161 * t186 + t162 * t187 + t805) * t775 + (t111 * t187 + t186 * t89 + t806) * t773 + (t846 * t773 + t775 * t869 + t847) * t277 + t849;
t777 = 2 * qJD(4);
t776 = -m(6) / 0.2e1;
t774 = -m(7) / 0.2e1;
t771 = m(6) * pkin(3);
t770 = m(7) * pkin(3);
t763 = mrSges(6,3) / 0.2e1;
t720 = pkin(3) * t461;
t311 = t326 + t720;
t152 = t462 * t311 + t831;
t86 = t152 + t548;
t762 = -t86 / 0.2e1;
t759 = -t277 / 0.2e1;
t756 = -t300 / 0.2e1;
t755 = t300 / 0.2e1;
t301 = t424 * t414;
t754 = t301 / 0.2e1;
t743 = -t322 / 0.2e1;
t727 = t462 / 0.4e1;
t722 = m(7) * t480;
t721 = m(7) * t442;
t719 = pkin(4) * t303;
t718 = pkin(4) * t426;
t717 = pkin(5) * t320;
t705 = mrSges(7,3) * t462;
t700 = t480 * mrSges(5,1);
t699 = t277 * mrSges(5,2);
t691 = t414 * mrSges(5,3);
t690 = t422 * mrSges(7,3);
t689 = t425 * mrSges(7,3);
t688 = t460 * mrSges(6,3);
t687 = t460 * mrSges(7,3);
t686 = t105 - t85;
t684 = t152 * t460;
t153 = t460 * t311 - t828;
t683 = t153 * t462;
t682 = t162 * t462;
t677 = t253 * t460;
t604 = t459 ^ 2 * t725;
t671 = t277 * t361;
t26 = m(5) * (-t362 * t480 - t464 * t604 - t671) + m(4) * (-t604 + (t390 * t463 + t461 * t498) * t459) * t464 + t814 * (t186 * t291 + t187 * t292 - t671);
t675 = t26 * qJD(1);
t673 = t480 * t423;
t672 = t480 * t424;
t670 = t277 * t455;
t669 = t277 * t457;
t664 = t352 * t361;
t659 = t398 * t320;
t657 = t399 * t316;
t651 = t420 * t302;
t650 = t422 * t320;
t649 = t425 * t316;
t648 = t441 * t303;
t647 = t442 * t302;
t645 = t460 * t321;
t643 = t462 * t317;
t640 = t638 * mrSges(7,3);
t636 = qJD(5) * t460;
t633 = mrSges(6,3) * t715;
t629 = -t726 / 0.2e1;
t625 = mrSges(6,3) * t683;
t622 = t361 * t773;
t619 = t763 + mrSges(7,3) / 0.2e1;
t615 = t440 * t645;
t614 = t440 * t643;
t613 = -t705 / 0.2e1;
t612 = t705 / 0.2e1;
t609 = -t688 / 0.2e1;
t608 = t688 / 0.2e1;
t607 = -t687 / 0.2e1;
t606 = t687 / 0.2e1;
t605 = m(7) * t686;
t602 = t460 * t726;
t601 = t462 * t726;
t584 = -t653 / 0.2e1;
t581 = t652 / 0.2e1;
t568 = t457 / 0.2e1 + t455 / 0.2e1;
t567 = t686 * t399;
t564 = (t398 - t422) * t460;
t563 = (t399 - t425) * t462;
t557 = (t420 + t442) * t769;
t556 = t627 / 0.2e1;
t550 = t568 * mrSges(6,3);
t547 = t841 + t788;
t412 = t423 * t716;
t546 = -t457 * t713 - t412;
t107 = t616 + t153;
t468 = m(5) * t720 * t580 + (t152 * t186 + t153 * t187 + t277 * t869 + t805) * t775 + (t107 * t187 + t186 * t86 + t846 * t277 + t806) * t773;
t2 = t468 + t804 + (-t319 * t728 - t323 * t732 + t572 - t799) * t277 + t849 + t874;
t325 = mrSges(5,1) * t413 - mrSges(5,2) * t414;
t4 = -pkin(2) * t427 + t107 * t318 + t153 * t319 + t86 * t322 + t152 * t323 + (Ifges(4,4) * t463 + (Ifges(4,1) - Ifges(4,2)) * t461) * t463 + (-Ifges(4,4) * t461 + pkin(3) * t325) * t461 + t720 * t840 + m(6) * (t141 * t152 + t142 * t153 + t666) + m(7) * (t106 * t107 + t85 * t86 + t857) + t845;
t536 = t2 * qJD(1) + t4 * qJD(2);
t476 = (pkin(4) * t361 + pkin(10) * t524) * t776 + (t291 * t422 - t292 * t425 - t361 * t442) * t774;
t6 = t549 * t361 + t476 + t779 + t844;
t7 = t111 * t318 + t162 * t319 + t89 * t322 + t161 * t323 + m(6) * (t141 * t161 + t142 * t162 + t666) + m(7) * (t106 * t111 + t85 * t89 + t857) + t845;
t535 = t6 * qJD(1) + t7 * qJD(2);
t306 = t537 * t414;
t307 = t539 * t414;
t308 = t541 * t414;
t309 = t543 * t414;
t13 = t105 * t318 + t141 * t319 - t142 * t323 + t253 * t300 + t352 * t301 + (t605 - t322) * t106 + ((-t85 * mrSges(7,3) - t141 * mrSges(6,3) + t306 / 0.2e1 + t307 / 0.2e1 - t782) * t460 + (t106 * mrSges(7,3) + t142 * mrSges(6,3) - t308 / 0.2e1 - t309 / 0.2e1 + (t304 - t723) * pkin(5) + t501) * t462) * t414;
t488 = t291 * t520 - t292 * t620;
t518 = -t605 / 0.2e1 + t742;
t14 = (-t301 / 0.2e1 + t756 + t566 / 0.2e1) * t277 + (-t619 * t652 + t518 + t740) * t187 + (t619 * t653 - t570) * t186 + t488;
t530 = -t14 * qJD(1) + t13 * qJD(2);
t44 = (m(7) * (t106 * t460 + t85 * t462) + t462 * t322 + t460 * t318) * t414;
t508 = m(7) * (t186 * t462 + t187 * t460) * t414;
t60 = -t622 - t508 / 0.2e1;
t529 = -qJD(1) * t60 + qJD(2) * t44;
t527 = t683 - t684;
t526 = -t161 * t460 + t682;
t449 = Ifges(7,5) * t462;
t450 = Ifges(6,5) * t462;
t521 = -pkin(5) * t705 + t449 + t450;
t519 = Ifges(6,1) / 0.4e1 + Ifges(7,1) / 0.4e1 - Ifges(6,2) / 0.4e1 - Ifges(7,2) / 0.4e1;
t512 = t717 / 0.2e1 + t837 * t737 + t712 * t655 / 0.2e1 - t838 * t654 / 0.2e1;
t509 = (-t420 * t652 + t677) * pkin(5);
t482 = (-t186 * t602 + t187 * t601 + t833) * pkin(3);
t466 = (t482 - t855) * t775 + (t482 + t834) * t773 + mrSges(5,2) * t758 + (t420 * t773 + t441 * t775 + t549) * t480 + t759 * t853;
t16 = -t466 + t673 / 0.2e1 + t672 / 0.2e1 - t700 / 0.2e1 + t699 / 0.2e1 + t477 + t808 * (-t670 / 0.2e1 - t669 / 0.2e1);
t472 = t808 * t507 * pkin(3) - mrSges(5,2) * t628 + t785 * t627;
t59 = (t398 * t602 + t399 * t601 + t420 * t724) * t770 + (t440 * t507 + t441 * t724) * t771 + t472;
t465 = (t526 * t440 + t829) * t775 + (t399 * t111 - t398 * t89 + t856) * t773 - t659 / 0.2e1 + t657 / 0.2e1 - t651 / 0.2e1 - t648 / 0.2e1 + t111 * t612 + t161 * t609 + t682 * t763 - t615 / 0.2e1 + t614 / 0.2e1 + t89 * t607 + t798 * t556 + ((-t141 * t602 + t142 * t601 + t830) * t775 + (t106 * t601 + t253 * t724 - t602 * t85) * t773 - t824 * t602 / 0.2e1 + t825 * t601 / 0.2e1) * pkin(3) - t862;
t467 = -t776 * t839 + (-t107 * t425 + t422 * t86 + t854) * t774 - t719 / 0.2e1 - t650 / 0.2e1 + t649 / 0.2e1 + t647 / 0.2e1 + t107 * t613 + t152 * t608 - t625 / 0.2e1 + t86 * t606 + (t527 * t776 + t645 / 0.2e1 - t643 / 0.2e1) * pkin(10) + t862;
t8 = t467 + t465;
t505 = -t16 * qJD(1) + t8 * qJD(2) + t59 * qJD(3);
t504 = (t656 / 0.2e1 + t658 / 0.2e1) * mrSges(7,3);
t428 = -Ifges(7,6) * t460 + t449;
t429 = -Ifges(6,6) * t460 + t450;
t474 = t630 * t858 - t784 + t105 * t612 + t85 * t613 - t304 * t621 + 0.2e1 * (t307 + t306) * t727 + (t309 + t308) * t730 + (t429 + t428) * t413 / 0.4e1 - (t198 + t196) * t460 / 0.4e1 + (t202 + t200) * t727 + t821 * t653 / 0.4e1 - t820 * t652 / 0.4e1 + (t608 + t609) * t142 + (t606 + t607) * t106;
t469 = t441 * t754 + t420 * t755 - t398 * t747 + t399 * t743 + (t414 * t550 + t794) * t440 + t474;
t490 = t107 * t764 - t152 * mrSges(6,1) / 0.2e1 + t153 * t765 + mrSges(7,1) * t762;
t10 = t469 - t717 / 0.2e1 + (t567 / 0.2e1 + t509 / 0.2e1 + pkin(5) * t762) * m(7) + (-t460 * t617 + t462 * t618) * t413 + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + t504) * t414 + t490;
t18 = (-t634 / 0.2e1 + t547) * t277 + t483;
t368 = t420 * t639;
t391 = t441 * t426;
t496 = t800 + (-Ifges(7,1) + Ifges(7,2) - Ifges(6,1) + Ifges(6,2)) * t462;
t48 = -t368 - t391 + (-t420 * t769 + t496) * t460 + t546;
t503 = -t18 * qJD(1) + t10 * qJD(2) - t48 * qJD(3);
t194 = -m(7) * t523 + t640;
t486 = m(7) * ((-t398 * t462 + t399 * t460) * t414 + t532);
t37 = -t486 / 0.2e1 + t866;
t62 = 0.2e1 * (t480 / 0.4e1 - t79 / 0.4e1) * m(7);
t502 = -qJD(1) * t62 - qJD(2) * t37 + qJD(3) * t194;
t500 = t809 * t629;
t392 = t442 * t639;
t499 = -t368 / 0.2e1 - t391 / 0.2e1 - t392 / 0.2e1 - t412 + t718 / 0.2e1;
t497 = pkin(4) * t754 - t422 * t318 / 0.2e1 + t442 * t756;
t475 = t555 * t89 - t111 * mrSges(7,2) / 0.2e1 + t161 * t767 - t162 * mrSges(6,2) / 0.2e1 + t512;
t11 = -t568 * t633 + t475 - t518 * t425 + (-t309 / 0.4e1 - t308 / 0.4e1 + t198 / 0.4e1 + t196 / 0.4e1 + pkin(10) * t746 + (-t723 / 0.2e1 + t304 / 0.2e1) * pkin(5) + (t690 / 0.2e1 - t519 * t460) * t414) * t460 + (-t307 / 0.4e1 - t306 / 0.4e1 - t202 / 0.4e1 - t200 / 0.4e1 + pkin(10) * t740 + (-t105 / 0.2e1 + t85 / 0.2e1) * mrSges(7,3) + (t689 / 0.2e1 + t519 * t462 - t800 + (t721 / 0.2e1 + t736) * pkin(5)) * t414) * t462 + (-t428 / 0.4e1 - t429 / 0.4e1) * t413 + t497 + t784;
t20 = (t551 + (t767 + t766) * t460 + t547) * t277;
t487 = t786 * t629;
t35 = (pkin(3) * t500 - t462 * t713) * t462 + (-t557 / 0.2e1 + t487 * pkin(3) + t496) * t460 + t499;
t56 = t718 - t392 + (-pkin(5) * t721 + t496) * t460 + t546;
t494 = -t20 * qJD(1) - t11 * qJD(2) - t35 * qJD(3) - t56 * qJD(4);
t109 = (t556 - t563 / 0.2e1 - t564 / 0.2e1) * m(7) - t640;
t267 = -m(7) * t522 + t640;
t485 = m(7) * ((t422 * t462 - t425 * t460) * t414 + t532);
t39 = -t485 / 0.2e1 + t866;
t77 = t79 * t773;
t65 = t77 - t722 / 0.2e1;
t491 = qJD(1) * t65 - qJD(2) * t39 - qJD(3) * t109 + qJD(4) * t267;
t389 = t397 * qJD(5);
t388 = t397 * qJD(6);
t110 = (t563 + t564) * t773 + m(7) * t556 + t640;
t64 = t77 + t722 / 0.2e1;
t63 = t480 * t773 + t77;
t61 = t508 / 0.2e1 - t622;
t40 = t485 / 0.2e1 + t864;
t38 = t486 / 0.2e1 + t864;
t36 = t513 / 0.2e1 - t515 / 0.2e1 + t514 / 0.2e1 - t516 / 0.2e1 + (t460 * t487 + t462 * t500) * pkin(3) - t499 + t821 * t728 + (t557 + t820) * t730;
t17 = t466 - (t568 * mrSges(7,3) + t550 + t873) * t277 + t549 * t480 + t477;
t15 = t488 + (-t630 * t773 + t754 + t755) * t277 + (t686 * t773 + t741 + t743) * t187 + t570 * t186 + t808 * (t186 * t584 + t187 * t581);
t12 = -t425 * t743 + t475 + (-t686 * t425 + (-t442 * t652 + t677) * pkin(5)) * t773 - t581 * t689 + t584 * t690 + t474 - t497 + t638 * t633 / 0.2e1 + t794 * pkin(10);
t9 = t469 + (t567 + t509) * t773 + t86 * t632 + t414 * t504 - t490 + t512;
t5 = -t476 + t778 + t779;
t3 = -t467 + t465 + t489;
t1 = t468 + (t758 + t759) * t692 + (-t427 / 0.2e1 - t324 / 0.2e1) * t646 + t781 + t847 * t277 - t874;
t19 = [qJD(2) * t26 + t859 * t871, t1 * qJD(3) + t5 * qJD(4) + t15 * qJD(5) + t61 * qJD(6) + t675 + (-m(6) * t664 + m(5) * (-t362 * t797 - t664) + t362 * t692 + m(4) * (pkin(8) * t464 * t822 - t725 * pkin(2)) * t459 + (-mrSges(3,1) + t325 + t789 + t840) * t603 - (-t691 + t798 + t723) * t361 + (m(6) * t142 + m(7) * t106 + t825) * t292 + (m(6) * t141 + m(7) * t85 + t824) * t291 + (t822 * mrSges(4,3) - mrSges(3,2)) * t646) * qJD(2), t1 * qJD(2) + (t672 + t673 + m(7) * (t420 * t480 + t834) + m(6) * (t441 * t480 - t855) + (-t480 * t726 - t833) * t772 + t699 - t700 + t498 * mrSges(4,2) - t390 * mrSges(4,1) + t808 * (-t669 - t670)) * qJD(3) + t17 * qJD(4) + t63 * qJD(6) + t872, t5 * qJD(2) + t17 * qJD(3) + t64 * qJD(6) + t477 * t777 + (t785 * t480 + (mrSges(5,2) - t853) * t277) * qJD(4) + t872, t15 * qJD(2) + (-t186 * t809 - t187 * t786) * qJD(5) + t871 * t865, qJD(2) * t61 + qJD(3) * t63 + qJD(4) * t64; qJD(3) * t2 + qJD(4) * t6 - qJD(5) * t14 - qJD(6) * t60 - t675, qJD(3) * t4 + qJD(4) * t7 + qJD(5) * t13 + qJD(6) * t44 (t614 - t615 - t659 + t657 + t625 + Ifges(4,5) * t463 - Ifges(4,6) * t461 + t107 * t705 + (-t726 * t797 - t830) * t772 + m(6) * (t440 * t527 + t829) + t627 * t691 + t628 * t692 - t648 + t789 * pkin(8) + m(7) * (t107 * t399 - t398 * t86 + t856) - t86 * t687 - mrSges(6,3) * t684 - t651 + t863) * qJD(3) + t3 * qJD(4) + t9 * qJD(5) + t38 * qJD(6) + t536, t3 * qJD(3) + t12 * qJD(5) + t40 * qJD(6) + ((pkin(10) * t526 - t839) * t775 + (-t111 * t425 + t422 * t89 + t854) * t773) * t777 + t535 + (-t647 - t649 + t650 + t719 + (t162 * mrSges(6,3) + t111 * mrSges(7,3) + pkin(10) * t317) * t462 + (-t161 * mrSges(6,3) - t89 * mrSges(7,3) - pkin(10) * t321) * t460 + t863) * qJD(4), t9 * qJD(3) + t12 * qJD(4) + t530 + (-mrSges(6,1) * t142 - mrSges(6,2) * t141 - mrSges(7,2) * t105 + (t852 + (-mrSges(7,3) * pkin(5) + t838) * t460) * t414 + t807 * t106) * qJD(5), qJD(3) * t38 + qJD(4) * t40 + t529; -qJD(2) * t2 - qJD(4) * t16 - qJD(5) * t18 - qJD(6) * t62 - t868, qJD(4) * t8 + qJD(5) * t10 - qJD(6) * t37 - t536, qJD(4) * t59 - qJD(5) * t48 + qJD(6) * t194 ((-t422 * t602 - t425 * t601 + t442 * t724) * t770 + (-pkin(4) * t724 + pkin(10) * t507) * t771 + t472) * qJD(4) + t36 * qJD(5) + t110 * qJD(6) + t505, t36 * qJD(4) + (mrSges(7,2) * t398 + t399 * t807 - t440 * t708 + t521) * qJD(5) + (mrSges(6,2) * t440 - t712) * t636 + t503, qJD(4) * t110 + t502; -qJD(2) * t6 + qJD(3) * t16 - qJD(5) * t20 + qJD(6) * t65 - t868, -qJD(3) * t8 - qJD(5) * t11 - qJD(6) * t39 - t535, -qJD(5) * t35 - qJD(6) * t109 - t505, -qJD(5) * t56 + qJD(6) * t267 (-mrSges(7,2) * t422 - pkin(10) * t708 - t425 * t807 + t521) * qJD(5) + (mrSges(6,2) * pkin(10) - t712) * t636 + t494, t491; qJD(2) * t14 + qJD(3) * t18 + qJD(4) * t20, -qJD(3) * t10 + qJD(4) * t11 + qJD(6) * t185 - t530, qJD(4) * t35 + t388 - t503, t388 - t494, 0, t783; qJD(2) * t60 + qJD(3) * t62 - qJD(4) * t65, qJD(3) * t37 + qJD(4) * t39 - qJD(5) * t185 - t529, qJD(4) * t109 - t389 - t502, -t389 - t491, -t783, 0;];
Cq  = t19;
