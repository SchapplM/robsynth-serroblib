% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:45:23
% EndTime: 2019-05-06 09:45:53
% DurationCPUTime: 27.06s
% Computational Cost: add. (389355->386), mult. (924976->489), div. (0->0), fcn. (674204->12), ass. (0->148)
t758 = -2 * qJD(3);
t729 = sin(qJ(2));
t733 = cos(qJ(2));
t750 = qJD(1) * qJD(2);
t714 = qJDD(1) * t729 + t733 * t750;
t730 = sin(qJ(1));
t734 = cos(qJ(1));
t720 = -g(1) * t734 - g(2) * t730;
t736 = qJD(1) ^ 2;
t709 = -pkin(1) * t736 + qJDD(1) * pkin(7) + t720;
t755 = t729 * t709;
t757 = pkin(2) * t736;
t666 = qJDD(2) * pkin(2) - t714 * qJ(3) - t755 + (qJ(3) * t750 + t729 * t757 - g(3)) * t733;
t694 = -g(3) * t729 + t733 * t709;
t715 = qJDD(1) * t733 - t729 * t750;
t753 = qJD(1) * t729;
t716 = qJD(2) * pkin(2) - qJ(3) * t753;
t723 = t733 ^ 2;
t668 = qJ(3) * t715 - qJD(2) * t716 - t723 * t757 + t694;
t725 = sin(pkin(10));
t756 = cos(pkin(10));
t703 = (t725 * t733 + t729 * t756) * qJD(1);
t647 = t666 * t756 - t725 * t668 + t703 * t758;
t752 = qJD(1) * t733;
t702 = t725 * t753 - t756 * t752;
t648 = t725 * t666 + t756 * t668 + t702 * t758;
t682 = mrSges(4,1) * t702 + mrSges(4,2) * t703;
t686 = t714 * t725 - t756 * t715;
t696 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t703;
t681 = pkin(3) * t702 - qJ(4) * t703;
t735 = qJD(2) ^ 2;
t634 = -pkin(3) * t735 + qJDD(2) * qJ(4) - t681 * t702 + t648;
t719 = t730 * g(1) - t734 * g(2);
t741 = -qJDD(1) * pkin(1) - t719;
t670 = -t715 * pkin(2) + qJDD(3) + t716 * t753 + (-qJ(3) * t723 - pkin(7)) * t736 + t741;
t687 = t714 * t756 + t725 * t715;
t637 = (qJD(2) * t702 - t687) * qJ(4) + (qJD(2) * t703 + t686) * pkin(3) + t670;
t724 = sin(pkin(11));
t726 = cos(pkin(11));
t692 = qJD(2) * t724 + t703 * t726;
t623 = -0.2e1 * qJD(4) * t692 - t724 * t634 + t726 * t637;
t676 = qJDD(2) * t724 + t687 * t726;
t691 = qJD(2) * t726 - t703 * t724;
t616 = (t691 * t702 - t676) * pkin(8) + (t691 * t692 + t686) * pkin(4) + t623;
t624 = 0.2e1 * qJD(4) * t691 + t726 * t634 + t724 * t637;
t673 = pkin(4) * t702 - pkin(8) * t692;
t675 = qJDD(2) * t726 - t687 * t724;
t690 = t691 ^ 2;
t618 = -pkin(4) * t690 + pkin(8) * t675 - t673 * t702 + t624;
t728 = sin(qJ(5));
t732 = cos(qJ(5));
t610 = t732 * t616 - t728 * t618;
t662 = t691 * t732 - t692 * t728;
t640 = qJD(5) * t662 + t675 * t728 + t676 * t732;
t663 = t691 * t728 + t692 * t732;
t685 = qJDD(5) + t686;
t701 = qJD(5) + t702;
t608 = (t662 * t701 - t640) * pkin(9) + (t662 * t663 + t685) * pkin(5) + t610;
t611 = t728 * t616 + t732 * t618;
t639 = -qJD(5) * t663 + t675 * t732 - t676 * t728;
t654 = pkin(5) * t701 - pkin(9) * t663;
t661 = t662 ^ 2;
t609 = -pkin(5) * t661 + pkin(9) * t639 - t654 * t701 + t611;
t727 = sin(qJ(6));
t731 = cos(qJ(6));
t606 = t608 * t731 - t609 * t727;
t649 = t662 * t731 - t663 * t727;
t622 = qJD(6) * t649 + t639 * t727 + t640 * t731;
t650 = t662 * t727 + t663 * t731;
t631 = -mrSges(7,1) * t649 + mrSges(7,2) * t650;
t697 = qJD(6) + t701;
t641 = -mrSges(7,2) * t697 + mrSges(7,3) * t649;
t683 = qJDD(6) + t685;
t602 = m(7) * t606 + mrSges(7,1) * t683 - mrSges(7,3) * t622 - t631 * t650 + t641 * t697;
t607 = t608 * t727 + t609 * t731;
t621 = -qJD(6) * t650 + t639 * t731 - t640 * t727;
t642 = mrSges(7,1) * t697 - mrSges(7,3) * t650;
t603 = m(7) * t607 - mrSges(7,2) * t683 + mrSges(7,3) * t621 + t631 * t649 - t642 * t697;
t596 = t731 * t602 + t727 * t603;
t651 = -mrSges(6,1) * t662 + mrSges(6,2) * t663;
t652 = -mrSges(6,2) * t701 + mrSges(6,3) * t662;
t594 = m(6) * t610 + mrSges(6,1) * t685 - mrSges(6,3) * t640 - t651 * t663 + t652 * t701 + t596;
t653 = mrSges(6,1) * t701 - mrSges(6,3) * t663;
t744 = -t602 * t727 + t731 * t603;
t595 = m(6) * t611 - mrSges(6,2) * t685 + mrSges(6,3) * t639 + t651 * t662 - t653 * t701 + t744;
t590 = t732 * t594 + t728 * t595;
t667 = -mrSges(5,1) * t691 + mrSges(5,2) * t692;
t671 = -mrSges(5,2) * t702 + mrSges(5,3) * t691;
t588 = m(5) * t623 + mrSges(5,1) * t686 - mrSges(5,3) * t676 - t667 * t692 + t671 * t702 + t590;
t672 = mrSges(5,1) * t702 - mrSges(5,3) * t692;
t745 = -t594 * t728 + t732 * t595;
t589 = m(5) * t624 - mrSges(5,2) * t686 + mrSges(5,3) * t675 + t667 * t691 - t672 * t702 + t745;
t746 = -t588 * t724 + t726 * t589;
t581 = m(4) * t648 - qJDD(2) * mrSges(4,2) - mrSges(4,3) * t686 - qJD(2) * t696 - t682 * t702 + t746;
t695 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t702;
t633 = -qJDD(2) * pkin(3) - t735 * qJ(4) + t703 * t681 + qJDD(4) - t647;
t625 = -t675 * pkin(4) - t690 * pkin(8) + t692 * t673 + t633;
t613 = -t639 * pkin(5) - t661 * pkin(9) + t663 * t654 + t625;
t742 = m(7) * t613 - t621 * mrSges(7,1) + t622 * mrSges(7,2) - t649 * t641 + t650 * t642;
t739 = m(6) * t625 - t639 * mrSges(6,1) + mrSges(6,2) * t640 - t662 * t652 + t653 * t663 + t742;
t737 = -m(5) * t633 + t675 * mrSges(5,1) - mrSges(5,2) * t676 + t691 * t671 - t672 * t692 - t739;
t605 = m(4) * t647 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t687 + qJD(2) * t695 - t682 * t703 + t737;
t576 = t725 * t581 + t756 * t605;
t693 = -t733 * g(3) - t755;
t713 = (-mrSges(3,1) * t733 + mrSges(3,2) * t729) * qJD(1);
t718 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t752;
t574 = m(3) * t693 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t714 + qJD(2) * t718 - t713 * t753 + t576;
t717 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t753;
t747 = t756 * t581 - t605 * t725;
t575 = m(3) * t694 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t715 - qJD(2) * t717 + t713 * t752 + t747;
t748 = -t574 * t729 + t733 * t575;
t566 = m(2) * t720 - mrSges(2,1) * t736 - qJDD(1) * mrSges(2,2) + t748;
t708 = -t736 * pkin(7) + t741;
t582 = t726 * t588 + t724 * t589;
t740 = m(4) * t670 + t686 * mrSges(4,1) + mrSges(4,2) * t687 + t702 * t695 + t696 * t703 + t582;
t738 = -m(3) * t708 + t715 * mrSges(3,1) - mrSges(3,2) * t714 - t717 * t753 + t718 * t752 - t740;
t578 = m(2) * t719 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t736 + t738;
t754 = t730 * t566 + t734 * t578;
t567 = t733 * t574 + t729 * t575;
t749 = t734 * t566 - t578 * t730;
t706 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t729 + Ifges(3,4) * t733) * qJD(1);
t705 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t729 + Ifges(3,2) * t733) * qJD(1);
t704 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t729 + Ifges(3,6) * t733) * qJD(1);
t679 = Ifges(4,1) * t703 - Ifges(4,4) * t702 + Ifges(4,5) * qJD(2);
t678 = Ifges(4,4) * t703 - Ifges(4,2) * t702 + Ifges(4,6) * qJD(2);
t677 = Ifges(4,5) * t703 - Ifges(4,6) * t702 + Ifges(4,3) * qJD(2);
t657 = Ifges(5,1) * t692 + Ifges(5,4) * t691 + Ifges(5,5) * t702;
t656 = Ifges(5,4) * t692 + Ifges(5,2) * t691 + Ifges(5,6) * t702;
t655 = Ifges(5,5) * t692 + Ifges(5,6) * t691 + Ifges(5,3) * t702;
t645 = Ifges(6,1) * t663 + Ifges(6,4) * t662 + Ifges(6,5) * t701;
t644 = Ifges(6,4) * t663 + Ifges(6,2) * t662 + Ifges(6,6) * t701;
t643 = Ifges(6,5) * t663 + Ifges(6,6) * t662 + Ifges(6,3) * t701;
t628 = Ifges(7,1) * t650 + Ifges(7,4) * t649 + Ifges(7,5) * t697;
t627 = Ifges(7,4) * t650 + Ifges(7,2) * t649 + Ifges(7,6) * t697;
t626 = Ifges(7,5) * t650 + Ifges(7,6) * t649 + Ifges(7,3) * t697;
t598 = mrSges(7,2) * t613 - mrSges(7,3) * t606 + Ifges(7,1) * t622 + Ifges(7,4) * t621 + Ifges(7,5) * t683 + t626 * t649 - t627 * t697;
t597 = -mrSges(7,1) * t613 + mrSges(7,3) * t607 + Ifges(7,4) * t622 + Ifges(7,2) * t621 + Ifges(7,6) * t683 - t626 * t650 + t628 * t697;
t584 = mrSges(6,2) * t625 - mrSges(6,3) * t610 + Ifges(6,1) * t640 + Ifges(6,4) * t639 + Ifges(6,5) * t685 - pkin(9) * t596 - t597 * t727 + t598 * t731 + t643 * t662 - t644 * t701;
t583 = -mrSges(6,1) * t625 + mrSges(6,3) * t611 + Ifges(6,4) * t640 + Ifges(6,2) * t639 + Ifges(6,6) * t685 - pkin(5) * t742 + pkin(9) * t744 + t731 * t597 + t727 * t598 - t663 * t643 + t701 * t645;
t570 = mrSges(5,2) * t633 - mrSges(5,3) * t623 + Ifges(5,1) * t676 + Ifges(5,4) * t675 + Ifges(5,5) * t686 - pkin(8) * t590 - t583 * t728 + t584 * t732 + t655 * t691 - t656 * t702;
t569 = -mrSges(5,1) * t633 + mrSges(5,3) * t624 + Ifges(5,4) * t676 + Ifges(5,2) * t675 + Ifges(5,6) * t686 - pkin(4) * t739 + pkin(8) * t745 + t732 * t583 + t728 * t584 - t692 * t655 + t702 * t657;
t568 = (-Ifges(5,3) - Ifges(4,2)) * t686 - t703 * t677 + Ifges(4,4) * t687 + t691 * t657 - t692 * t656 + qJD(2) * t679 - Ifges(7,3) * t683 - Ifges(6,3) * t685 - mrSges(4,1) * t670 - Ifges(5,6) * t675 - Ifges(5,5) * t676 + t662 * t645 - t663 * t644 + mrSges(4,3) * t648 + t649 * t628 - t650 * t627 - Ifges(6,6) * t639 - Ifges(6,5) * t640 - Ifges(7,6) * t621 - Ifges(7,5) * t622 - mrSges(5,1) * t623 + mrSges(5,2) * t624 - mrSges(6,1) * t610 + mrSges(6,2) * t611 + mrSges(7,2) * t607 - mrSges(7,1) * t606 - pkin(5) * t596 - pkin(4) * t590 - pkin(3) * t582 + Ifges(4,6) * qJDD(2);
t563 = mrSges(4,2) * t670 - mrSges(4,3) * t647 + Ifges(4,1) * t687 - Ifges(4,4) * t686 + Ifges(4,5) * qJDD(2) - qJ(4) * t582 - qJD(2) * t678 - t569 * t724 + t570 * t726 - t677 * t702;
t562 = mrSges(3,2) * t708 - mrSges(3,3) * t693 + Ifges(3,1) * t714 + Ifges(3,4) * t715 + Ifges(3,5) * qJDD(2) - qJ(3) * t576 - qJD(2) * t705 + t563 * t756 - t725 * t568 + t704 * t752;
t561 = mrSges(2,1) * g(3) - t702 * t679 - t703 * t678 + t736 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - Ifges(4,5) * t687 + Ifges(4,6) * t686 - mrSges(4,1) * t647 + mrSges(4,2) * t648 - t724 * t570 - t726 * t569 - pkin(3) * t737 - qJ(4) * t746 - Ifges(3,5) * t714 - Ifges(3,6) * t715 - mrSges(3,1) * t693 + mrSges(3,2) * t694 - pkin(2) * t576 - pkin(1) * t567 + mrSges(2,3) * t720 + (-Ifges(4,3) - Ifges(3,3)) * qJDD(2) + (-t705 * t729 + t706 * t733) * qJD(1);
t560 = -mrSges(3,1) * t708 + mrSges(3,3) * t694 + Ifges(3,4) * t714 + Ifges(3,2) * t715 + Ifges(3,6) * qJDD(2) - pkin(2) * t740 + qJ(3) * t747 + qJD(2) * t706 + t725 * t563 + t568 * t756 - t704 * t753;
t559 = -mrSges(2,2) * g(3) - mrSges(2,3) * t719 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t736 - pkin(7) * t567 - t560 * t729 + t562 * t733;
t1 = [-m(1) * g(1) + t749; -m(1) * g(2) + t754; (-m(1) - m(2)) * g(3) + t567; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t754 + t734 * t559 - t730 * t561; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t749 + t730 * t559 + t734 * t561; -mrSges(1,1) * g(2) + mrSges(2,1) * t719 + mrSges(1,2) * g(1) - mrSges(2,2) * t720 + Ifges(2,3) * qJDD(1) + pkin(1) * t738 + pkin(7) * t748 + t733 * t560 + t729 * t562;];
tauB  = t1;
