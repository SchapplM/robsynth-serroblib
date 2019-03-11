% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:48
% EndTime: 2019-03-09 08:55:24
% DurationCPUTime: 19.35s
% Computational Cost: add. (51475->896), mult. (129099->1238), div. (0->0), fcn. (148283->12), ass. (0->453)
t520 = sin(pkin(11));
t523 = sin(qJ(2));
t697 = sin(pkin(6));
t698 = cos(pkin(11));
t603 = t698 * t697;
t764 = cos(qJ(2));
t607 = t764 * t697;
t472 = t520 * t607 + t523 * t603;
t519 = sin(pkin(12));
t521 = cos(pkin(12));
t699 = cos(pkin(6));
t437 = -t472 * t521 - t519 * t699;
t559 = t472 * t519 - t521 * t699;
t762 = sin(qJ(5));
t763 = cos(qJ(5));
t342 = t437 * t762 - t559 * t763;
t847 = -t342 / 0.2e1;
t846 = pkin(10) * t342;
t840 = -t437 * t763 - t559 * t762;
t844 = -t840 / 0.2e1;
t845 = mrSges(6,1) * t844;
t843 = pkin(5) * t840;
t616 = t523 * t697;
t471 = t520 * t616 - t603 * t764;
t842 = Ifges(3,5) * t607 - Ifges(4,5) * t471 - Ifges(3,6) * t616 - Ifges(4,6) * t472;
t637 = pkin(1) * t699;
t486 = pkin(8) * t607 + t523 * t637;
t507 = t764 * t637;
t841 = -t486 * mrSges(3,1) - (-pkin(8) * t616 + t507) * mrSges(3,2);
t464 = t507 + (-pkin(8) - qJ(3)) * t616;
t448 = pkin(2) * t699 + t464;
t465 = qJ(3) * t607 + t486;
t449 = t520 * t465;
t360 = t448 * t698 - t449;
t548 = -pkin(3) * t699 - t360;
t708 = t472 * mrSges(4,3);
t839 = m(4) * t360 - m(5) * t548 + mrSges(4,1) * t699 - mrSges(5,1) * t559 + mrSges(5,2) * t437 - t708;
t838 = -mrSges(5,3) / 0.2e1;
t716 = t840 * mrSges(6,3);
t522 = sin(qJ(6));
t517 = t522 ^ 2;
t524 = cos(qJ(6));
t518 = t524 ^ 2;
t837 = (t518 / 0.2e1 + t517 / 0.2e1) * mrSges(7,3);
t497 = -mrSges(7,1) * t524 + t522 * mrSges(7,2);
t836 = -t497 + mrSges(6,1);
t279 = t471 * t524 - t522 * t840;
t280 = t471 * t522 + t524 * t840;
t165 = -mrSges(7,1) * t279 + mrSges(7,2) * t280;
t286 = mrSges(6,1) * t471 - t716;
t835 = t165 - t286;
t514 = Ifges(7,4) * t524;
t503 = Ifges(7,1) * t522 + t514;
t651 = t519 ^ 2 + t521 ^ 2;
t493 = t519 * t762 - t521 * t763;
t382 = t471 * t493;
t314 = -t522 * t382 + t472 * t524;
t495 = -t519 * t763 - t521 * t762;
t381 = t471 * t495;
t230 = -mrSges(7,2) * t381 + mrSges(7,3) * t314;
t315 = t382 * t524 + t472 * t522;
t231 = mrSges(7,1) * t381 - t315 * mrSges(7,3);
t766 = -t524 / 0.2e1;
t767 = t522 / 0.2e1;
t834 = t230 * t766 + t231 * t767;
t702 = t524 * mrSges(7,2);
t704 = t522 * mrSges(7,1);
t498 = t702 + t704;
t422 = t498 * t495;
t706 = t495 * mrSges(6,3);
t833 = -t706 / 0.2e1 - t422 / 0.2e1;
t150 = Ifges(7,1) * t315 + Ifges(7,4) * t314 + t381 * Ifges(7,5);
t832 = -pkin(10) * t231 / 0.2e1 + t150 / 0.4e1;
t149 = Ifges(7,4) * t315 + Ifges(7,2) * t314 + t381 * Ifges(7,6);
t831 = pkin(10) * t230 / 0.2e1 + t149 / 0.4e1;
t760 = pkin(2) * t520;
t508 = qJ(4) + t760;
t755 = pkin(9) + t508;
t484 = t755 * t521;
t617 = t755 * t519;
t412 = t484 * t762 + t617 * t763;
t758 = pkin(5) * t495;
t436 = pkin(10) * t493 - t758;
t317 = t522 * t412 + t436 * t524;
t318 = -t412 * t524 + t522 * t436;
t586 = -t317 * t522 + t318 * t524;
t377 = t464 * t698 - t449;
t612 = pkin(2) * t616;
t400 = t472 * pkin(3) + t471 * qJ(4) + t612;
t249 = -t377 * t519 + t521 * t400;
t250 = t521 * t377 + t519 * t400;
t830 = -t249 * t519 + t250 * t521;
t513 = Ifges(7,5) * t524;
t732 = Ifges(7,6) * t522;
t596 = t732 - t513;
t598 = Ifges(7,2) * t522 - t514;
t727 = t279 * mrSges(7,3);
t192 = mrSges(7,2) * t342 + t727;
t726 = t280 * mrSges(7,3);
t193 = -mrSges(7,1) * t342 - t726;
t765 = t524 / 0.2e1;
t768 = -t522 / 0.2e1;
t829 = mrSges(7,3) * (t279 * t768 + t280 * t765) - t192 * t768 - t193 * t766;
t828 = -m(7) * pkin(5) - t836;
t576 = t513 / 0.2e1 - t732 / 0.2e1;
t744 = Ifges(6,5) * t471;
t827 = -t576 * t342 + t744 / 0.2e1;
t615 = t698 * t465;
t652 = t520 * t448 + t615;
t351 = qJ(4) * t699 + t652;
t635 = t697 * pkin(1);
t543 = -pkin(2) * t607 - t635;
t371 = t471 * pkin(3) - t472 * qJ(4) + t543;
t237 = -t351 * t519 + t521 * t371;
t180 = pkin(4) * t471 + pkin(9) * t437 + t237;
t238 = t521 * t351 + t519 * t371;
t203 = -pkin(9) * t559 + t238;
t100 = t180 * t763 - t203 * t762;
t92 = -t471 * pkin(5) - t100;
t826 = m(7) * t92 + t835;
t825 = m(5) / 0.2e1;
t824 = m(6) / 0.2e1;
t823 = -m(7) / 0.2e1;
t822 = m(7) / 0.2e1;
t821 = -pkin(5) / 0.2e1;
t820 = m(4) * pkin(2);
t819 = -mrSges(6,2) / 0.2e1;
t818 = mrSges(7,2) / 0.2e1;
t817 = Ifges(6,1) / 0.2e1;
t816 = -t92 / 0.2e1;
t734 = Ifges(7,6) * t342;
t747 = Ifges(7,4) * t280;
t122 = Ifges(7,2) * t279 - t734 + t747;
t815 = -t122 / 0.2e1;
t812 = -t192 / 0.2e1;
t281 = pkin(4) * t559 + t548;
t811 = t281 / 0.2e1;
t810 = t286 / 0.2e1;
t636 = t698 * pkin(2);
t512 = -t636 - pkin(3);
t496 = -t521 * pkin(4) + t512;
t408 = t493 * pkin(5) + t495 * pkin(10) + t496;
t413 = t484 * t763 - t617 * t762;
t295 = t522 * t408 + t413 * t524;
t809 = -t295 / 0.2e1;
t808 = -t317 / 0.2e1;
t806 = t342 / 0.2e1;
t803 = -t342 / 0.4e1;
t800 = t840 / 0.2e1;
t798 = -t381 / 0.2e1;
t797 = t381 / 0.2e1;
t733 = Ifges(7,6) * t493;
t394 = t495 * t598 + t733;
t796 = t394 / 0.2e1;
t746 = Ifges(7,4) * t522;
t504 = Ifges(7,1) * t524 - t746;
t740 = Ifges(7,5) * t493;
t396 = -t495 * t504 + t740;
t795 = t396 / 0.2e1;
t794 = t412 / 0.2e1;
t793 = t413 / 0.2e1;
t792 = t422 / 0.2e1;
t672 = t495 * t522;
t429 = -mrSges(7,2) * t493 + mrSges(7,3) * t672;
t790 = -t429 / 0.2e1;
t789 = t429 / 0.2e1;
t671 = t495 * t524;
t431 = t493 * mrSges(7,1) + mrSges(7,3) * t671;
t788 = -t431 / 0.2e1;
t787 = -t471 / 0.2e1;
t786 = t471 / 0.2e1;
t784 = t472 / 0.2e1;
t783 = -t493 / 0.2e1;
t781 = t493 / 0.2e1;
t780 = t493 / 0.4e1;
t779 = -t495 / 0.2e1;
t776 = t496 / 0.2e1;
t775 = -t497 / 0.2e1;
t774 = t497 / 0.2e1;
t773 = t498 / 0.2e1;
t772 = -t596 / 0.4e1;
t771 = t503 / 0.4e1;
t770 = t519 / 0.2e1;
t769 = t521 / 0.2e1;
t650 = t517 + t518;
t613 = t650 * t493;
t399 = -pkin(10) * t613 + t758;
t761 = m(7) * t399;
t759 = pkin(5) * t381;
t133 = -pkin(5) * t342 - pkin(10) * t840 + t281;
t101 = t180 * t762 + t203 * t763;
t93 = t471 * pkin(10) + t101;
t61 = t133 * t524 - t522 * t93;
t757 = t61 * mrSges(7,3);
t62 = t522 * t133 + t524 * t93;
t756 = t62 * mrSges(7,3);
t752 = Ifges(4,4) * t472;
t751 = Ifges(5,4) * t519;
t750 = Ifges(5,4) * t521;
t749 = Ifges(6,4) * t840;
t748 = Ifges(6,4) * t495;
t745 = Ifges(6,5) * t382;
t743 = Ifges(7,5) * t280;
t742 = Ifges(7,5) * t315;
t741 = Ifges(7,5) * t342;
t738 = Ifges(6,6) * t471;
t737 = Ifges(6,6) * t495;
t736 = Ifges(7,6) * t279;
t735 = Ifges(7,6) * t314;
t731 = Ifges(6,3) * t472;
t730 = Ifges(7,3) * t840;
t729 = Ifges(7,3) * t381;
t728 = Ifges(7,3) * t495;
t725 = t281 * mrSges(6,1);
t294 = t408 * t524 - t522 * t413;
t724 = t294 * mrSges(7,3);
t723 = t295 * mrSges(7,3);
t676 = t471 * t521;
t204 = pkin(4) * t472 + pkin(9) * t676 + t249;
t677 = t471 * t519;
t229 = pkin(9) * t677 + t250;
t110 = t204 * t763 - t229 * t762;
t103 = -t472 * pkin(5) - t110;
t111 = t762 * t204 + t763 * t229;
t278 = Ifges(7,4) * t279;
t123 = Ifges(7,1) * t280 + t278 - t741;
t720 = t315 * mrSges(7,2);
t721 = t314 * mrSges(7,1);
t198 = t720 - t721;
t718 = t342 * mrSges(6,3);
t285 = -mrSges(6,2) * t471 + t718;
t376 = t464 * t520 + t615;
t316 = -pkin(4) * t677 + t376;
t319 = -t472 * mrSges(6,2) - t381 * mrSges(6,3);
t320 = mrSges(6,1) * t472 - mrSges(6,3) * t382;
t710 = t471 * mrSges(5,2);
t372 = -mrSges(5,3) * t559 - t710;
t711 = t471 * mrSges(5,1);
t373 = t437 * mrSges(5,3) + t711;
t404 = -t472 * mrSges(5,2) + mrSges(5,3) * t677;
t405 = t472 * mrSges(5,1) + mrSges(5,3) * t676;
t544 = Ifges(5,6) * t472 + (Ifges(5,2) * t519 - t750) * t471;
t545 = Ifges(5,5) * t472 + (-Ifges(5,1) * t521 + t751) * t471;
t554 = t729 + t735 + t742;
t555 = -Ifges(7,3) * t342 + t736 + t743;
t556 = Ifges(6,4) * t382 - Ifges(6,2) * t381 + Ifges(6,6) * t472;
t557 = Ifges(6,2) * t342 + t738 + t749;
t558 = Ifges(6,1) * t382 - Ifges(6,4) * t381 + Ifges(6,5) * t472;
t560 = (-mrSges(5,1) * t519 - mrSges(5,2) * t521) * t471;
t709 = t471 * mrSges(4,3);
t569 = mrSges(4,2) * t699 + t709;
t335 = Ifges(6,4) * t342;
t578 = Ifges(6,1) * t840 + t335 + t744;
t580 = -t607 / 0.2e1;
t714 = t382 * mrSges(6,2);
t715 = t381 * mrSges(6,1);
t600 = t714 + t715;
t601 = -mrSges(6,1) * t342 + mrSges(6,2) * t840;
t606 = t616 / 0.2e1;
t611 = Ifges(3,4) * t616;
t614 = t472 * mrSges(4,1) - t471 * mrSges(4,2);
t634 = -t677 / 0.2e1;
t104 = pkin(10) * t472 + t111;
t176 = -pkin(10) * t382 + t316 + t759;
t73 = -t522 * t104 + t176 * t524;
t74 = t104 * t524 + t522 * t176;
t3 = (0.2e1 * Ifges(3,4) * t607 + (Ifges(3,1) - Ifges(3,2)) * t616) * t580 + (-0.2e1 * Ifges(4,4) * t471 + (Ifges(4,1) - Ifges(4,2)) * t472) * t786 + (-Ifges(5,4) * t437 - Ifges(5,2) * t559 + Ifges(5,6) * t471) * t634 + (-Ifges(5,1) * t437 - Ifges(5,4) * t559 + Ifges(5,5) * t471) * t676 / 0.2e1 + t437 * t545 / 0.2e1 + (Ifges(3,2) * t607 + t611) * t606 + (-Ifges(4,2) * t471 + t752) * t784 + (Ifges(5,3) * t472 + (-Ifges(5,5) * t521 + Ifges(5,6) * t519) * t471 - Ifges(6,6) * t381 + t731 + t745) * t787 + t652 * t708 - m(6) * (t100 * t110 + t101 * t111 + t281 * t316) - m(7) * (t103 * t92 + t61 * t73 + t62 * t74) + t839 * t376 - m(4) * (t652 * t377 + (-pkin(2) * t764 - pkin(1)) * t523 * pkin(2) * t697 ^ 2) - t360 * t709 + (-t842 / 0.2e1 + Ifges(3,5) * t580 + Ifges(4,5) * t786 + Ifges(3,6) * t606 + Ifges(4,6) * t784 - t841) * t699 - (-Ifges(5,5) * t437 + Ifges(6,5) * t840 - Ifges(5,6) * t559 + Ifges(6,6) * t342 - t752 + (-Ifges(4,1) + Ifges(5,3) + Ifges(6,3)) * t471) * t472 / 0.2e1 + t558 * t844 + t556 * t847 - (Ifges(3,1) * t607 - t611) * t616 / 0.2e1 - t543 * t614 + t314 * t815 - (mrSges(4,1) * t471 + mrSges(4,2) * t472) * t612 - t281 * t600 - t316 * t601 - t382 * t578 / 0.2e1 + t377 * t569 - t548 * t560 + t559 * t544 / 0.2e1 + (mrSges(3,1) * t616 + mrSges(3,2) * t607) * t635 + t557 * t797 + t555 * t798 + t554 * t806 - m(5) * (t237 * t249 + t238 * t250) - t103 * t165 - t74 * t192 - t73 * t193 - t92 * t198 - t62 * t230 - t61 * t231 - t279 * t149 / 0.2e1 - t280 * t150 / 0.2e1 - t111 * t285 - t110 * t286 - t315 * t123 / 0.2e1 - t101 * t319 - t100 * t320 - t250 * t372 - t249 * t373 - t238 * t404 - t237 * t405;
t722 = t3 * qJD(1);
t719 = t342 * mrSges(6,2);
t717 = t840 * mrSges(6,1);
t145 = Ifges(7,6) * t840 - t342 * t598;
t146 = Ifges(7,5) * t840 + t342 * t504;
t220 = t498 * t342;
t685 = t342 * t522;
t221 = -mrSges(7,2) * t840 - mrSges(7,3) * t685;
t684 = t342 * t524;
t222 = mrSges(7,1) * t840 - mrSges(7,3) * t684;
t643 = -t738 / 0.2e1;
t653 = Ifges(6,5) * t342 - Ifges(6,6) * t840;
t660 = t524 * t123;
t668 = t522 * t122;
t225 = t843 - t846;
t78 = -t522 * t100 + t225 * t524;
t79 = t100 * t524 + t522 * t225;
t4 = t279 * t145 / 0.2e1 + t79 * t192 + t62 * t221 + t78 * t193 + t61 * t222 + t280 * t146 / 0.2e1 + t92 * t220 + m(7) * (t61 * t78 + t62 * t79) + t653 * t786 + t100 * t285 + (t743 / 0.2e1 + t736 / 0.2e1 + t725 + t643 - t749) * t840 - (-t281 * mrSges(6,2) - t335 + t668 / 0.2e1 - t660 / 0.2e1 + t100 * mrSges(6,3) + (Ifges(7,3) + Ifges(6,2) - Ifges(6,1)) * t840 - t827) * t342 + (t826 - t716) * t101;
t713 = t4 * qJD(1);
t712 = t412 * mrSges(6,3);
t707 = t493 * mrSges(6,3);
t705 = t496 * mrSges(6,1);
t703 = t522 * t73;
t701 = t524 * t74;
t164 = mrSges(7,1) * t280 + mrSges(7,2) * t279;
t166 = Ifges(7,5) * t279 - Ifges(7,6) * t280;
t167 = -Ifges(7,2) * t280 + t278;
t168 = Ifges(7,1) * t279 - t747;
t7 = t61 * t192 - t62 * t193 + t166 * t847 + t92 * t164 + (t815 - t756 + t168 / 0.2e1) * t280 + (t167 / 0.2e1 - t757 + t123 / 0.2e1) * t279;
t700 = t7 * qJD(1);
t599 = mrSges(6,1) * t493 - mrSges(6,2) * t495;
t602 = -t521 * mrSges(5,1) + t519 * mrSges(5,2);
t640 = -t707 / 0.2e1;
t682 = t381 * t412;
t525 = (-t471 * t508 * t651 + t472 * t512) * t825 + (t382 * t413 + t472 * t496 + t682) * t824 + (t294 * t314 + t295 * t315 + t682) * t822 + t314 * t431 / 0.2e1 + t315 * t789 + (-t471 * t520 - t472 * t698) * t820 / 0.2e1 + t382 * t640 + (t602 + t599) * t784 + t651 * mrSges(5,3) * t787 + t833 * t381;
t594 = -t701 + t703;
t530 = (t249 * t521 + t250 * t519) * t825 + (-t110 * t493 - t111 * t495) * t824 + (t493 * t103 + t495 * t594) * t822 + t404 * t770 + t405 * t769 + t606 * t820;
t15 = t198 * t781 + t319 * t779 + t320 * t783 + t495 * t834 - t525 + t530 + t614;
t696 = qJD(1) * t15;
t595 = t522 * t61 - t524 * t62;
t658 = t524 * t192;
t666 = t522 * t193;
t17 = -t559 * t372 + t437 * t373 + t835 * t840 + (t285 + t658 - t666) * t342 + m(7) * (-t342 * t595 + t840 * t92) + m(6) * (-t100 * t840 + t101 * t342) + m(5) * (t237 * t437 - t238 * t559);
t695 = qJD(1) * t17;
t589 = t237 * t519 - t238 * t521;
t669 = t521 * t372;
t670 = t519 * t373;
t14 = t315 * t192 + t314 * t193 + m(7) * (t314 * t61 + t315 * t62) + (m(6) * t101 + t285) * t382 + (-m(6) * t100 + t826) * t381 + (m(6) * t281 + t601 - t839) * t472 + (-m(4) * t652 + m(5) * t589 + t569 - t669 + t670) * t471;
t694 = t14 * qJD(1);
t691 = t314 * t522;
t690 = t315 * t524;
t688 = t317 * t524;
t686 = t840 * t493;
t683 = t342 * t495;
t681 = t412 * t840;
t680 = t412 * t495;
t679 = t437 * t519;
t678 = t559 * t521;
t675 = t493 * t381;
t674 = t493 * t522;
t673 = t493 * t524;
t667 = t522 * t145;
t665 = t522 * t222;
t664 = t522 * t318;
t663 = t522 * t394;
t430 = -t495 * mrSges(7,1) + mrSges(7,3) * t673;
t662 = t522 * t430;
t661 = t522 * t431;
t659 = t524 * t146;
t657 = t524 * t221;
t656 = t524 * t396;
t428 = mrSges(7,2) * t495 + mrSges(7,3) * t674;
t655 = t524 * t428;
t654 = t524 * t429;
t212 = m(7) * (-0.1e1 + t650) * t495 * t493;
t211 = t212 / 0.2e1;
t421 = t493 * t498;
t564 = t661 / 0.2e1 - t654 / 0.2e1;
t588 = t294 * t522 - t295 * t524;
t43 = ((-t412 - t586) * t822 - t655 / 0.2e1 + t662 / 0.2e1 + t792) * t495 + ((t413 + t588) * t822 - t421 / 0.2e1 + t564) * t493;
t420 = mrSges(7,1) * t671 - mrSges(7,2) * t672;
t86 = t420 * t783 + (t429 * t767 + t431 * t765 - t495 * t837) * t495;
t648 = t211 * qJD(4) + t43 * qJD(5) + t86 * qJD(6);
t645 = -Ifges(6,2) / 0.2e1 - Ifges(7,3) / 0.2e1;
t642 = -t726 / 0.2e1;
t641 = -t716 / 0.2e1;
t633 = -t676 / 0.2e1;
t632 = t674 / 0.2e1;
t631 = -t673 / 0.2e1;
t630 = t672 / 0.2e1;
t629 = -t671 / 0.2e1;
t626 = t122 / 0.4e1 - t168 / 0.4e1;
t625 = -t165 / 0.2e1 + t810;
t624 = t167 / 0.4e1 + t123 / 0.4e1;
t623 = -t220 / 0.2e1 + t285 / 0.2e1;
t425 = t503 * t495;
t622 = -t394 / 0.4e1 + t425 / 0.4e1;
t501 = Ifges(7,2) * t524 + t746;
t424 = t501 * t495;
t621 = t396 / 0.4e1 + t424 / 0.4e1;
t620 = t771 - t598 / 0.4e1;
t619 = t504 / 0.4e1 - t501 / 0.4e1;
t605 = t812 + t727 / 0.2e1;
t604 = t642 - t193 / 0.2e1;
t499 = Ifges(7,5) * t522 + Ifges(7,6) * t524;
t593 = -t522 * t78 + t524 * t79;
t592 = t165 * t779 + t192 * t631 + t193 * t632 + t285 * t783 + t495 * t810;
t534 = t164 * t781 + t495 * t829;
t575 = t721 / 0.2e1 - t720 / 0.2e1;
t22 = t534 - t575;
t591 = -qJD(1) * t22 - qJD(2) * t86;
t536 = (t522 * t79 + t524 * t78) * t823 - t719 / 0.2e1 + t221 * t768 + t222 * t766;
t539 = (t650 * t846 - t843) * t822 + t840 * t774;
t25 = 0.2e1 * t845 + (t819 + t837) * t342 + t536 + t539;
t487 = t493 * mrSges(6,2);
t565 = t428 * t768 + t430 * t766;
t581 = t493 * t837;
t89 = t487 + (t775 + mrSges(6,1)) * t495 - t581 + 0.2e1 * (t399 / 0.4e1 - t688 / 0.4e1 - t664 / 0.4e1) * m(7) + t565;
t590 = qJD(1) * t25 + qJD(2) * t89;
t587 = -t690 + t691;
t572 = t702 / 0.2e1 + t704 / 0.2e1;
t553 = t572 * t493;
t112 = t553 + t564;
t28 = (mrSges(7,2) * t847 + t605) * t524 + (mrSges(7,1) * t847 - t604) * t522;
t585 = qJD(1) * t28 + qJD(2) * t112;
t533 = (-t650 * t683 + t686) * t822 + (-t683 + t686) * t824 + (t437 * t521 - t519 * t559) * t825;
t561 = m(7) * (t314 * t524 + t522 * t315);
t49 = -t561 / 0.2e1 + 0.2e1 * (-m(6) / 0.4e1 - m(5) / 0.4e1) * t472 + t533;
t584 = qJD(1) * t49 + qJD(2) * t211;
t579 = -t593 - t92;
t577 = t101 + t595;
t574 = mrSges(7,1) * t808 + t318 * t818;
t573 = -t715 / 0.2e1 - t714 / 0.2e1;
t571 = pkin(10) * t788 + t621;
t570 = pkin(10) * t790 + t622;
t566 = -t666 / 0.2e1 + t658 / 0.2e1;
t562 = t501 * t767 + t503 * t766;
t423 = t495 * t499;
t392 = Ifges(7,3) * t493 + t495 * t596;
t433 = -Ifges(6,2) * t493 - t748;
t488 = Ifges(6,4) * t493;
t434 = -Ifges(6,1) * t495 - t488;
t393 = -Ifges(7,6) * t495 + t493 * t598;
t395 = -Ifges(7,5) * t495 - t493 * t504;
t528 = (t101 * t412 + t294 * t78 + t295 * t79 + t317 * t61 + t318 * t62 + t413 * t92) * t823 + t101 * t792 - t279 * t393 / 0.4e1 - t280 * t395 / 0.4e1 - t294 * t222 / 0.2e1 + t221 * t809 + t193 * t808 + t318 * t812 - t61 * t430 / 0.2e1 - t62 * t428 / 0.2e1 + t78 * t788 + t79 * t790 - t421 * t816;
t529 = (t499 / 0.4e1 - Ifges(6,6) / 0.2e1) * t381 + (-pkin(5) * t103 - pkin(10) * t594) * t822 + t745 / 0.2e1 + t731 / 0.2e1 + t198 * t821 + t103 * t774 + t110 * mrSges(6,1) / 0.2e1 + t111 * t819 + t314 * t501 / 0.4e1 + t315 * t771;
t2 = (t660 / 0.4e1 - t668 / 0.4e1 + t335 / 0.2e1 + (-Ifges(7,3) / 0.4e1 + t817 - Ifges(6,2) / 0.4e1) * t840 + t827) * t493 - (-t488 / 0.4e1 + t656 / 0.4e1 - t663 / 0.4e1 + t712 / 0.2e1 + mrSges(6,2) * t776 + t434 / 0.4e1) * t342 + t487 * t811 + t529 + t528 + t623 * t412 + (t74 * mrSges(7,3) / 0.2e1 + t831) * t524 + (-t73 * mrSges(7,3) / 0.2e1 + t832) * t522 + (mrSges(6,3) * t793 - t705 / 0.2e1 - t392 / 0.4e1 + t433 / 0.4e1) * t840 + t625 * t413 + (t659 / 0.4e1 - t667 / 0.4e1 + t743 / 0.4e1 + t736 / 0.4e1 + t643 - 0.3e1 / 0.4e1 * t749 + t725 / 0.2e1 - (-Ifges(6,1) / 0.4e1 - t645) * t342) * t495;
t27 = t318 * t429 + t295 * t428 + m(7) * (t294 * t317 + t295 * t318 + t412 * t413) - t413 * t422 - t412 * t421 + t317 * t431 + t294 * t430 - t496 * t487 + (-t705 - t392 / 0.2e1 + t433 / 0.2e1 + t395 * t766 + t393 * t767 - t748 / 0.2e1) * t495 + (-t434 / 0.2e1 + t488 / 0.2e1 - t656 / 0.2e1 + t663 / 0.2e1 - t576 * t493 + (t817 + t645) * t495) * t493;
t552 = -t2 * qJD(1) + t27 * qJD(2) + t43 * qJD(3);
t38 = t294 * t429 - t412 * t420 - t295 * t431 + t423 * t781 + ((t723 - t425 / 0.2e1 + t796) * t524 + (-t724 + t795 + t424 / 0.2e1) * t522) * t495;
t527 = ((t756 / 0.2e1 + t626) * t524 + (-t757 / 0.2e1 + t624) * t522) * t495 + (-t724 / 0.2e1 + t621) * t279 + (-t723 / 0.2e1 + t622) * t280 + t294 * t192 / 0.2e1 + t193 * t809 + t423 * t803 + t164 * t794 + t166 * t780 + t61 * t789 + t62 * t788 + t420 * t816;
t538 = t742 / 0.2e1 + t735 / 0.2e1 + t729 / 0.2e1 + t73 * mrSges(7,1) / 0.2e1 - t74 * mrSges(7,2) / 0.2e1;
t5 = t527 - t538;
t551 = t5 * qJD(1) + t38 * qJD(2) + t86 * qJD(3);
t531 = (t690 / 0.2e1 - t691 / 0.2e1) * mrSges(7,3) + (-pkin(10) * t587 - t759) * t822 + t381 * t774 + t573;
t8 = (t657 / 0.2e1 - t665 / 0.2e1 + t579 * t823 + t641 - t625) * t495 + (t577 * t823 - t718 / 0.2e1 + t566 + t623) * t493 + t531;
t550 = -t8 * qJD(1) + t43 * qJD(2) + t212 * qJD(3);
t526 = t833 * t840 + (t640 - t564) * t342 + ((-t678 - t679) * t508 - t589) * t825 + (t100 * t495 - t493 * t101 + t342 * t413 + t681) * t824 + (-t342 * t588 + t493 * t595 - t495 * t92 + t681) * t822 + t592;
t532 = -m(5) * t376 / 0.2e1 - m(6) * t316 / 0.2e1 + (t522 * t74 + t524 * t73) * t823 + t230 * t768 + t231 * t766 + t573;
t11 = t526 + (t710 / 0.2e1 + t559 * t838 + t372 / 0.2e1) * t521 + (t711 / 0.2e1 + t437 * t838 - t373 / 0.2e1) * t519 + t532;
t63 = (t422 + t706) * t495 + (-t654 + t661 + t707) * t493 + m(7) * (t493 * t588 - t680) + m(6) * (-t493 * t413 - t680) + (m(5) * t508 + mrSges(5,3)) * t651;
t549 = -qJD(1) * t11 - qJD(2) * t63 - qJD(3) * t211;
t547 = -t730 / 0.2e1 - t78 * mrSges(7,1) / 0.2e1 + t79 * t818;
t546 = pkin(5) * t420 / 0.2e1 + t412 * t773 + t493 * t772;
t535 = t164 * t821 + t279 * t620 + t280 * t619 - t342 * t772 + t773 * t92;
t12 = (-t741 / 0.2e1 + t604 * pkin(10) + t624) * t524 + (t734 / 0.2e1 + t605 * pkin(10) - t626) * t522 + t535 + t547;
t346 = pkin(5) * t498 + (-t503 / 0.2e1 + t598 / 0.2e1) * t524 + (-t504 / 0.2e1 + t501 / 0.2e1) * t522;
t352 = (-t498 / 0.2e1 + t572) * t493;
t537 = pkin(10) * t837 + t522 * t620 - t524 * t619;
t40 = (t740 / 0.2e1 + t571) * t524 + (-t733 / 0.2e1 + t570) * t522 + (Ifges(7,3) / 0.2e1 + t537) * t495 + t546 + t574;
t542 = t12 * qJD(1) + t40 * qJD(2) - t352 * qJD(3) - t346 * qJD(5);
t353 = t493 * t773 + t553;
t113 = t553 - t564;
t96 = (t664 + t688) * t822 + t761 / 0.2e1 + t495 * t775 - t581 - t565;
t48 = t561 / 0.2e1 + t533 + (m(6) + m(5)) * t784;
t39 = Ifges(7,5) * t631 + Ifges(7,6) * t632 - t728 / 0.2e1 + t571 * t524 + t570 * t522 + t537 * t495 + t546 - t574;
t29 = -t342 * t572 + t522 * t642 + t727 * t766 + t566;
t26 = t717 / 0.2e1 + mrSges(6,2) * t847 + t845 + t342 * t837 - t536 + t539;
t23 = t534 + t575;
t16 = t525 + (-t319 / 0.2e1 + t834) * t495 + (-t320 / 0.2e1 + t198 / 0.2e1) * t493 + t530;
t13 = Ifges(7,5) * t684 / 0.2e1 - Ifges(7,6) * t685 / 0.2e1 + t624 * t524 - t626 * t522 - t829 * pkin(10) + t535 - t547;
t10 = t526 + mrSges(5,2) * t633 + mrSges(5,1) * t634 + t669 / 0.2e1 - t670 / 0.2e1 + (-t678 / 0.2e1 - t679 / 0.2e1) * mrSges(5,3) - t532;
t9 = t220 * t781 + t221 * t629 + t222 * t630 + (t493 * t577 + t495 * t579) * t822 + t706 * t800 - t342 * t640 + t531 + t592;
t6 = t527 + t538;
t1 = t529 - t528 + t471 * (-Ifges(6,5) * t493 + t737) / 0.4e1 + (t701 / 0.2e1 - t703 / 0.2e1) * mrSges(7,3) + (-t495 * mrSges(6,1) - t487) * t811 + (t717 + t719) * t776 + t165 * t793 + t220 * t794 + t712 * t806 - t412 * t285 / 0.2e1 - t840 * t433 / 0.4e1 + (t493 * t596 + t663 - t728) * t803 + (-Ifges(6,1) * t493 + t392 + t748) * t840 / 0.4e1 + (-t342 * t596 + t668 + t730) * t780 + (t667 + t557) * t495 / 0.4e1 + t831 * t524 + t832 * t522 + (t641 - t286 / 0.2e1) * t413 + (Ifges(6,2) * t495 + t434 - t488 + t656) * t342 / 0.4e1 - (-Ifges(6,2) * t840 + t335 + t578 + t660) * t493 / 0.4e1 - (Ifges(6,1) * t342 + t555 + t659 - t749) * t495 / 0.4e1;
t18 = [-qJD(2) * t3 + qJD(3) * t14 + qJD(4) * t17 + qJD(5) * t4 + qJD(6) * t7, -t722 + (t841 + t110 * t706 + t636 * t709 + t150 * t629 + t149 * t630 - t708 * t760 + (Ifges(5,1) * t519 + t750) * t633 + (Ifges(5,2) * t521 + t751) * t677 / 0.2e1 - t111 * t707 + t512 * t560 + t842 + (m(5) * t830 + t521 * t404 - t519 * t405) * t508 + t830 * mrSges(5,3) + t316 * t599 + t496 * t600 + m(6) * (t111 * t413 + t316 * t496) + m(7) * (t294 * t73 + t295 * t74) + (Ifges(5,5) * t519 - Ifges(6,5) * t495 + Ifges(5,6) * t521 - Ifges(6,6) * t493) * t784 + (-m(6) * t110 + m(7) * t103 + t198 - t320) * t412 + t544 * t769 + t545 * t770 + t558 * t779 + t554 * t781 + t556 * t783 + t315 * t795 + t314 * t796 + t392 * t797 + t433 * t798 + t294 * t231 + t295 * t230 + t413 * t319 - t103 * t422 + t74 * t429 + t73 * t431 + t382 * t434 / 0.2e1 + (m(5) * t512 - t698 * t820 - mrSges(4,1) + t602) * t376 + (t520 * t820 - mrSges(4,2)) * t377) * qJD(2) + t16 * qJD(3) + t10 * qJD(4) + t1 * qJD(5) + t6 * qJD(6), t694 + t16 * qJD(2) + 0.2e1 * ((t495 * t587 + t675) * t822 + (-t382 * t495 + t675) * t824) * qJD(3) + t48 * qJD(4) + t9 * qJD(5) + t23 * qJD(6), qJD(2) * t10 + qJD(3) * t48 + qJD(5) * t26 + qJD(6) * t29 + t695, t713 + t1 * qJD(2) + t9 * qJD(3) + t26 * qJD(4) + (-t100 * mrSges(6,2) + t593 * mrSges(7,3) - pkin(5) * t220 + t145 * t765 + t146 * t767 - t562 * t342 + t499 * t800 + t653 + t828 * t101 + (m(7) * t593 + t657 - t665) * pkin(10)) * qJD(5) + t13 * qJD(6), t700 + t6 * qJD(2) + t23 * qJD(3) + t29 * qJD(4) + t13 * qJD(5) + (-mrSges(7,1) * t62 - mrSges(7,2) * t61 + t166) * qJD(6); -qJD(3) * t15 + qJD(4) * t11 - qJD(5) * t2 + qJD(6) * t5 + t722, qJD(4) * t63 + qJD(5) * t27 + qJD(6) * t38, t648 - t696, qJD(5) * t96 + qJD(6) * t113 - t549, t96 * qJD(4) + t39 * qJD(6) + t552 + (t395 * t767 + t393 * t765 + pkin(5) * t421 + t499 * t779 + t737 + t412 * mrSges(6,2) + (-Ifges(6,5) + t562) * t493 + t828 * t413 + (m(7) * t586 + t655 - t662) * pkin(10) + t586 * mrSges(7,3)) * qJD(5), t113 * qJD(4) + t39 * qJD(5) + (-t295 * mrSges(7,1) - t294 * mrSges(7,2) + t423) * qJD(6) + t551; qJD(2) * t15 + qJD(4) * t49 - qJD(5) * t8 + qJD(6) * t22 - t694, t648 + t696, t212 * qJD(5), t584 (-mrSges(7,3) * t613 + t495 * t836 + t487 + t761) * qJD(5) + t353 * qJD(6) + t550, qJD(5) * t353 + qJD(6) * t420 - t591; -qJD(2) * t11 - qJD(3) * t49 - qJD(5) * t25 - qJD(6) * t28 - t695, -qJD(5) * t89 - qJD(6) * t112 + t549, -t584, 0, -t590, -qJD(6) * t498 - t585; qJD(2) * t2 + qJD(3) * t8 + qJD(4) * t25 + qJD(6) * t12 - t713, qJD(4) * t89 + qJD(6) * t40 - t552, -qJD(6) * t352 - t550, t590, -t346 * qJD(6) (pkin(10) * t497 - t596) * qJD(6) + t542; -qJD(2) * t5 - qJD(3) * t22 + qJD(4) * t28 - qJD(5) * t12 - t700, qJD(4) * t112 - qJD(5) * t40 - t551, qJD(5) * t352 + t591, t585, -t542, 0;];
Cq  = t18;
