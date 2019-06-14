% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 21:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 21:19:04
% EndTime: 2019-05-06 21:19:55
% DurationCPUTime: 49.64s
% Computational Cost: add. (760752->397), mult. (1989524->516), div. (0->0), fcn. (1595520->14), ass. (0->163)
t805 = -2 * qJD(3);
t764 = sin(pkin(12));
t766 = cos(pkin(12));
t771 = sin(qJ(2));
t776 = cos(qJ(2));
t765 = sin(pkin(6));
t797 = qJD(1) * t765;
t745 = (t764 * t771 - t766 * t776) * t797;
t795 = qJD(1) * qJD(2);
t754 = (qJDD(1) * t771 + t776 * t795) * t765;
t767 = cos(pkin(6));
t759 = qJDD(1) * t767 + qJDD(2);
t760 = qJD(1) * t767 + qJD(2);
t772 = sin(qJ(1));
t777 = cos(qJ(1));
t756 = t772 * g(1) - g(2) * t777;
t778 = qJD(1) ^ 2;
t804 = pkin(8) * t765;
t751 = qJDD(1) * pkin(1) + t778 * t804 + t756;
t757 = -g(1) * t777 - g(2) * t772;
t752 = -pkin(1) * t778 + qJDD(1) * t804 + t757;
t799 = t767 * t776;
t785 = t751 * t799 - t771 * t752;
t803 = t765 ^ 2 * t778;
t689 = t759 * pkin(2) - t754 * qJ(3) + (pkin(2) * t771 * t803 + (qJ(3) * qJD(1) * t760 - g(3)) * t765) * t776 + t785;
t800 = t767 * t771;
t802 = t765 * t771;
t718 = -g(3) * t802 + t751 * t800 + t776 * t752;
t793 = t771 * t797;
t748 = pkin(2) * t760 - qJ(3) * t793;
t755 = (qJDD(1) * t776 - t771 * t795) * t765;
t794 = t776 ^ 2 * t803;
t693 = -pkin(2) * t794 + qJ(3) * t755 - t748 * t760 + t718;
t746 = (t764 * t776 + t766 * t771) * t797;
t669 = t766 * t689 - t764 * t693 + t746 * t805;
t801 = t765 * t776;
t670 = t764 * t689 + t766 * t693 + t745 * t805;
t719 = mrSges(4,1) * t745 + mrSges(4,2) * t746;
t724 = -t754 * t764 + t755 * t766;
t732 = mrSges(4,1) * t760 - mrSges(4,3) * t746;
t720 = pkin(3) * t745 - pkin(9) * t746;
t758 = t760 ^ 2;
t663 = -pkin(3) * t758 + pkin(9) * t759 - t720 * t745 + t670;
t736 = -t767 * g(3) - t765 * t751;
t704 = -t755 * pkin(2) - qJ(3) * t794 + t748 * t793 + qJDD(3) + t736;
t725 = t754 * t766 + t755 * t764;
t672 = (t745 * t760 - t725) * pkin(9) + (t746 * t760 - t724) * pkin(3) + t704;
t770 = sin(qJ(4));
t775 = cos(qJ(4));
t659 = t775 * t663 + t770 * t672;
t730 = t746 * t775 + t760 * t770;
t701 = -t730 * qJD(4) - t770 * t725 + t759 * t775;
t729 = -t770 * t746 + t760 * t775;
t705 = -mrSges(5,1) * t729 + mrSges(5,2) * t730;
t744 = qJD(4) + t745;
t712 = mrSges(5,1) * t744 - mrSges(5,3) * t730;
t723 = qJDD(4) - t724;
t706 = -pkin(4) * t729 - pkin(10) * t730;
t743 = t744 ^ 2;
t651 = -pkin(4) * t743 + pkin(10) * t723 + t706 * t729 + t659;
t662 = -t759 * pkin(3) - t758 * pkin(9) + t746 * t720 - t669;
t702 = qJD(4) * t729 + t725 * t775 + t759 * t770;
t654 = (-t729 * t744 - t702) * pkin(10) + (t730 * t744 - t701) * pkin(4) + t662;
t769 = sin(qJ(5));
t774 = cos(qJ(5));
t646 = -t769 * t651 + t774 * t654;
t709 = -t730 * t769 + t744 * t774;
t675 = qJD(5) * t709 + t702 * t774 + t723 * t769;
t700 = qJDD(5) - t701;
t710 = t730 * t774 + t744 * t769;
t728 = qJD(5) - t729;
t644 = (t709 * t728 - t675) * pkin(11) + (t709 * t710 + t700) * pkin(5) + t646;
t647 = t774 * t651 + t769 * t654;
t674 = -qJD(5) * t710 - t702 * t769 + t723 * t774;
t692 = pkin(5) * t728 - pkin(11) * t710;
t708 = t709 ^ 2;
t645 = -pkin(5) * t708 + pkin(11) * t674 - t692 * t728 + t647;
t768 = sin(qJ(6));
t773 = cos(qJ(6));
t642 = t644 * t773 - t645 * t768;
t682 = t709 * t773 - t710 * t768;
t657 = qJD(6) * t682 + t674 * t768 + t675 * t773;
t683 = t709 * t768 + t710 * t773;
t668 = -mrSges(7,1) * t682 + mrSges(7,2) * t683;
t726 = qJD(6) + t728;
t676 = -mrSges(7,2) * t726 + mrSges(7,3) * t682;
t698 = qJDD(6) + t700;
t640 = m(7) * t642 + mrSges(7,1) * t698 - mrSges(7,3) * t657 - t668 * t683 + t676 * t726;
t643 = t644 * t768 + t645 * t773;
t656 = -qJD(6) * t683 + t674 * t773 - t675 * t768;
t677 = mrSges(7,1) * t726 - mrSges(7,3) * t683;
t641 = m(7) * t643 - mrSges(7,2) * t698 + mrSges(7,3) * t656 + t668 * t682 - t677 * t726;
t632 = t773 * t640 + t768 * t641;
t684 = -mrSges(6,1) * t709 + mrSges(6,2) * t710;
t690 = -mrSges(6,2) * t728 + mrSges(6,3) * t709;
t630 = m(6) * t646 + mrSges(6,1) * t700 - mrSges(6,3) * t675 - t684 * t710 + t690 * t728 + t632;
t691 = mrSges(6,1) * t728 - mrSges(6,3) * t710;
t787 = -t640 * t768 + t773 * t641;
t631 = m(6) * t647 - mrSges(6,2) * t700 + mrSges(6,3) * t674 + t684 * t709 - t691 * t728 + t787;
t788 = -t630 * t769 + t774 * t631;
t627 = m(5) * t659 - mrSges(5,2) * t723 + mrSges(5,3) * t701 + t705 * t729 - t712 * t744 + t788;
t658 = -t770 * t663 + t672 * t775;
t711 = -mrSges(5,2) * t744 + mrSges(5,3) * t729;
t650 = -pkin(4) * t723 - pkin(10) * t743 + t730 * t706 - t658;
t648 = -pkin(5) * t674 - pkin(11) * t708 + t692 * t710 + t650;
t782 = m(7) * t648 - t656 * mrSges(7,1) + mrSges(7,2) * t657 - t682 * t676 + t677 * t683;
t779 = -m(6) * t650 + t674 * mrSges(6,1) - mrSges(6,2) * t675 + t709 * t690 - t691 * t710 - t782;
t636 = m(5) * t658 + mrSges(5,1) * t723 - mrSges(5,3) * t702 - t705 * t730 + t711 * t744 + t779;
t789 = t775 * t627 - t636 * t770;
t617 = m(4) * t670 - mrSges(4,2) * t759 + mrSges(4,3) * t724 - t719 * t745 - t732 * t760 + t789;
t731 = -mrSges(4,2) * t760 - mrSges(4,3) * t745;
t628 = t630 * t774 + t631 * t769;
t780 = -m(5) * t662 + t701 * mrSges(5,1) - mrSges(5,2) * t702 + t729 * t711 - t712 * t730 - t628;
t624 = m(4) * t669 + mrSges(4,1) * t759 - mrSges(4,3) * t725 - t719 * t746 + t731 * t760 + t780;
t613 = t764 * t617 + t766 * t624;
t717 = -g(3) * t801 + t785;
t792 = t776 * t797;
t750 = -mrSges(3,2) * t760 + mrSges(3,3) * t792;
t753 = (-mrSges(3,1) * t776 + mrSges(3,2) * t771) * t797;
t611 = m(3) * t717 + mrSges(3,1) * t759 - mrSges(3,3) * t754 + t750 * t760 - t753 * t793 + t613;
t749 = mrSges(3,1) * t760 - mrSges(3,3) * t793;
t790 = t766 * t617 - t624 * t764;
t612 = m(3) * t718 - mrSges(3,2) * t759 + mrSges(3,3) * t755 - t749 * t760 + t753 * t792 + t790;
t620 = t770 * t627 + t775 * t636;
t781 = m(4) * t704 - t724 * mrSges(4,1) + t725 * mrSges(4,2) + t745 * t731 + t746 * t732 + t620;
t619 = m(3) * t736 - t755 * mrSges(3,1) + t754 * mrSges(3,2) + (t749 * t771 - t750 * t776) * t797 + t781;
t599 = t611 * t799 + t612 * t800 - t619 * t765;
t597 = m(2) * t756 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t778 + t599;
t604 = -t611 * t771 + t776 * t612;
t603 = m(2) * t757 - mrSges(2,1) * t778 - qJDD(1) * mrSges(2,2) + t604;
t798 = t777 * t597 + t772 * t603;
t598 = t611 * t801 + t612 * t802 + t767 * t619;
t791 = -t597 * t772 + t777 * t603;
t664 = Ifges(7,5) * t683 + Ifges(7,6) * t682 + Ifges(7,3) * t726;
t666 = Ifges(7,1) * t683 + Ifges(7,4) * t682 + Ifges(7,5) * t726;
t633 = -mrSges(7,1) * t648 + mrSges(7,3) * t643 + Ifges(7,4) * t657 + Ifges(7,2) * t656 + Ifges(7,6) * t698 - t664 * t683 + t666 * t726;
t665 = Ifges(7,4) * t683 + Ifges(7,2) * t682 + Ifges(7,6) * t726;
t634 = mrSges(7,2) * t648 - mrSges(7,3) * t642 + Ifges(7,1) * t657 + Ifges(7,4) * t656 + Ifges(7,5) * t698 + t664 * t682 - t665 * t726;
t678 = Ifges(6,5) * t710 + Ifges(6,6) * t709 + Ifges(6,3) * t728;
t680 = Ifges(6,1) * t710 + Ifges(6,4) * t709 + Ifges(6,5) * t728;
t621 = -mrSges(6,1) * t650 + mrSges(6,3) * t647 + Ifges(6,4) * t675 + Ifges(6,2) * t674 + Ifges(6,6) * t700 - pkin(5) * t782 + pkin(11) * t787 + t773 * t633 + t768 * t634 - t710 * t678 + t728 * t680;
t679 = Ifges(6,4) * t710 + Ifges(6,2) * t709 + Ifges(6,6) * t728;
t622 = mrSges(6,2) * t650 - mrSges(6,3) * t646 + Ifges(6,1) * t675 + Ifges(6,4) * t674 + Ifges(6,5) * t700 - pkin(11) * t632 - t633 * t768 + t634 * t773 + t678 * t709 - t679 * t728;
t694 = Ifges(5,5) * t730 + Ifges(5,6) * t729 + Ifges(5,3) * t744;
t695 = Ifges(5,4) * t730 + Ifges(5,2) * t729 + Ifges(5,6) * t744;
t605 = mrSges(5,2) * t662 - mrSges(5,3) * t658 + Ifges(5,1) * t702 + Ifges(5,4) * t701 + Ifges(5,5) * t723 - pkin(10) * t628 - t621 * t769 + t622 * t774 + t694 * t729 - t695 * t744;
t696 = Ifges(5,1) * t730 + Ifges(5,4) * t729 + Ifges(5,5) * t744;
t614 = Ifges(5,4) * t702 + Ifges(5,2) * t701 + Ifges(5,6) * t723 - t730 * t694 + t744 * t696 - mrSges(5,1) * t662 + mrSges(5,3) * t659 - Ifges(6,5) * t675 - Ifges(6,6) * t674 - Ifges(6,3) * t700 - t710 * t679 + t709 * t680 - mrSges(6,1) * t646 + mrSges(6,2) * t647 - Ifges(7,5) * t657 - Ifges(7,6) * t656 - Ifges(7,3) * t698 - t683 * t665 + t682 * t666 - mrSges(7,1) * t642 + mrSges(7,2) * t643 - pkin(5) * t632 - pkin(4) * t628;
t713 = Ifges(4,5) * t746 - Ifges(4,6) * t745 + Ifges(4,3) * t760;
t714 = Ifges(4,4) * t746 - Ifges(4,2) * t745 + Ifges(4,6) * t760;
t595 = mrSges(4,2) * t704 - mrSges(4,3) * t669 + Ifges(4,1) * t725 + Ifges(4,4) * t724 + Ifges(4,5) * t759 - pkin(9) * t620 + t605 * t775 - t614 * t770 - t713 * t745 - t714 * t760;
t715 = Ifges(4,1) * t746 - Ifges(4,4) * t745 + Ifges(4,5) * t760;
t600 = Ifges(4,4) * t725 + Ifges(4,2) * t724 + Ifges(4,6) * t759 - t746 * t713 + t760 * t715 - mrSges(4,1) * t704 + mrSges(4,3) * t670 - Ifges(5,5) * t702 - Ifges(5,6) * t701 - Ifges(5,3) * t723 - t730 * t695 + t729 * t696 - mrSges(5,1) * t658 + mrSges(5,2) * t659 - t769 * t622 - t774 * t621 - pkin(4) * t779 - pkin(10) * t788 - pkin(3) * t620;
t733 = Ifges(3,3) * t760 + (Ifges(3,5) * t771 + Ifges(3,6) * t776) * t797;
t735 = Ifges(3,5) * t760 + (Ifges(3,1) * t771 + Ifges(3,4) * t776) * t797;
t592 = -mrSges(3,1) * t736 + mrSges(3,3) * t718 + Ifges(3,4) * t754 + Ifges(3,2) * t755 + Ifges(3,6) * t759 - pkin(2) * t781 + qJ(3) * t790 + t764 * t595 + t766 * t600 - t733 * t793 + t760 * t735;
t734 = Ifges(3,6) * t760 + (Ifges(3,4) * t771 + Ifges(3,2) * t776) * t797;
t593 = mrSges(3,2) * t736 - mrSges(3,3) * t717 + Ifges(3,1) * t754 + Ifges(3,4) * t755 + Ifges(3,5) * t759 - qJ(3) * t613 + t595 * t766 - t600 * t764 + t733 * t792 - t734 * t760;
t783 = pkin(8) * t604 + t592 * t776 + t593 * t771;
t594 = Ifges(3,5) * t754 + Ifges(3,6) * t755 + mrSges(3,1) * t717 - mrSges(3,2) * t718 + Ifges(4,5) * t725 + Ifges(4,6) * t724 + t746 * t714 + t745 * t715 + mrSges(4,1) * t669 - mrSges(4,2) * t670 + t770 * t605 + t775 * t614 + pkin(3) * t780 + pkin(9) * t789 + pkin(2) * t613 + (Ifges(3,3) + Ifges(4,3)) * t759 + (t734 * t771 - t735 * t776) * t797;
t591 = -mrSges(2,2) * g(3) - mrSges(2,3) * t756 + Ifges(2,5) * qJDD(1) - t778 * Ifges(2,6) - t771 * t592 + t776 * t593 + (-t598 * t765 - t599 * t767) * pkin(8);
t590 = mrSges(2,1) * g(3) + mrSges(2,3) * t757 + t778 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t598 - t765 * t594 + t783 * t767;
t1 = [-m(1) * g(1) + t791; -m(1) * g(2) + t798; (-m(1) - m(2)) * g(3) + t598; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t798 - t772 * t590 + t777 * t591; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t791 + t777 * t590 + t772 * t591; -mrSges(1,1) * g(2) + mrSges(2,1) * t756 + mrSges(1,2) * g(1) - mrSges(2,2) * t757 + Ifges(2,3) * qJDD(1) + pkin(1) * t599 + t767 * t594 + t783 * t765;];
tauB  = t1;
