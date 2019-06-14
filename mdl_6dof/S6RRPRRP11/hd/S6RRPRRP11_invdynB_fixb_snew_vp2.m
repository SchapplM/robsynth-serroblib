% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP11
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 18:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP11_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:53:26
% EndTime: 2019-05-06 18:53:38
% DurationCPUTime: 6.37s
% Computational Cost: add. (66801->354), mult. (137174->405), div. (0->0), fcn. (84263->8), ass. (0->137)
t810 = -2 * qJD(3);
t809 = Ifges(3,1) + Ifges(4,2);
t808 = Ifges(6,1) + Ifges(7,1);
t800 = Ifges(3,4) + Ifges(4,6);
t799 = Ifges(6,4) + Ifges(7,4);
t798 = Ifges(3,5) - Ifges(4,4);
t797 = Ifges(6,5) + Ifges(7,5);
t807 = Ifges(3,2) + Ifges(4,3);
t806 = Ifges(6,2) + Ifges(7,2);
t796 = Ifges(3,6) - Ifges(4,5);
t795 = Ifges(6,6) + Ifges(7,6);
t805 = Ifges(3,3) + Ifges(4,1);
t804 = Ifges(6,3) + Ifges(7,3);
t754 = sin(qJ(4));
t758 = cos(qJ(4));
t759 = cos(qJ(2));
t784 = qJD(1) * t759;
t726 = t758 * qJD(2) - t754 * t784;
t755 = sin(qJ(2));
t783 = qJD(1) * qJD(2);
t777 = t755 * t783;
t731 = t759 * qJDD(1) - t777;
t692 = -t726 * qJD(4) - t754 * qJDD(2) - t758 * t731;
t725 = -t754 * qJD(2) - t758 * t784;
t693 = t725 * qJD(4) + t758 * qJDD(2) - t754 * t731;
t753 = sin(qJ(5));
t757 = cos(qJ(5));
t696 = t753 * t725 + t757 * t726;
t657 = -t696 * qJD(5) + t757 * t692 - t753 * t693;
t778 = (mrSges(6,1) + mrSges(7,1)) * t657;
t695 = t757 * t725 - t753 * t726;
t746 = t755 * qJD(1);
t743 = t746 + qJD(4);
t741 = qJD(5) + t743;
t681 = -t741 * mrSges(7,2) + t695 * mrSges(7,3);
t682 = -t741 * mrSges(6,2) + t695 * mrSges(6,3);
t803 = -(t681 + t682) * t695 - t778;
t756 = sin(qJ(1));
t760 = cos(qJ(1));
t740 = -t760 * g(1) - t756 * g(2);
t762 = qJD(1) ^ 2;
t716 = -t762 * pkin(1) + qJDD(1) * pkin(7) + t740;
t701 = -t755 * g(3) + t759 * t716;
t727 = (-pkin(2) * t759 - qJ(3) * t755) * qJD(1);
t761 = qJD(2) ^ 2;
t679 = t761 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t810 - t727 * t784 - t701;
t802 = t762 * pkin(7);
t794 = t695 * t681;
t700 = -t759 * g(3) - t755 * t716;
t728 = (mrSges(4,2) * t759 - mrSges(4,3) * t755) * qJD(1);
t729 = (-mrSges(3,1) * t759 + mrSges(3,2) * t755) * qJD(1);
t776 = t759 * t783;
t730 = t755 * qJDD(1) + t776;
t735 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t784;
t736 = -mrSges(4,1) * t784 - qJD(2) * mrSges(4,3);
t738 = pkin(3) * t746 - qJD(2) * pkin(8);
t752 = t759 ^ 2;
t739 = t756 * g(1) - t760 * g(2);
t771 = -qJDD(1) * pkin(1) - t739;
t765 = pkin(2) * t777 + t746 * t810 + (-t730 - t776) * qJ(3) + t771;
t661 = -t738 * t746 + (-pkin(3) * t752 - pkin(7)) * t762 + (-pkin(2) - pkin(8)) * t731 + t765;
t680 = -qJDD(2) * pkin(2) - t761 * qJ(3) + t727 * t746 + qJDD(3) - t700;
t666 = (-t755 * t759 * t762 - qJDD(2)) * pkin(8) + (t730 - t776) * pkin(3) + t680;
t650 = -t754 * t661 + t758 * t666;
t724 = qJDD(4) + t730;
t647 = (t725 * t743 - t693) * pkin(9) + (t725 * t726 + t724) * pkin(4) + t650;
t651 = t758 * t661 + t754 * t666;
t702 = t743 * pkin(4) - t726 * pkin(9);
t723 = t725 ^ 2;
t649 = -t723 * pkin(4) + t692 * pkin(9) - t743 * t702 + t651;
t641 = t757 * t647 - t753 * t649;
t658 = t695 * qJD(5) + t753 * t692 + t757 * t693;
t675 = -t695 * mrSges(7,1) + t696 * mrSges(7,2);
t676 = -t695 * mrSges(6,1) + t696 * mrSges(6,2);
t717 = qJDD(5) + t724;
t638 = -0.2e1 * qJD(6) * t696 + (t695 * t741 - t658) * qJ(6) + (t695 * t696 + t717) * pkin(5) + t641;
t781 = m(7) * t638 + t717 * mrSges(7,1) + t741 * t681;
t629 = m(6) * t641 + t717 * mrSges(6,1) + t741 * t682 + (-t675 - t676) * t696 + (-mrSges(6,3) - mrSges(7,3)) * t658 + t781;
t642 = t753 * t647 + t757 * t649;
t684 = t741 * mrSges(7,1) - t696 * mrSges(7,3);
t685 = t741 * mrSges(6,1) - t696 * mrSges(6,3);
t683 = t741 * pkin(5) - t696 * qJ(6);
t694 = t695 ^ 2;
t640 = -t694 * pkin(5) + t657 * qJ(6) + 0.2e1 * qJD(6) * t695 - t741 * t683 + t642;
t780 = m(7) * t640 + t657 * mrSges(7,3) + t695 * t675;
t632 = m(6) * t642 + t657 * mrSges(6,3) + t695 * t676 + (-t684 - t685) * t741 + (-mrSges(6,2) - mrSges(7,2)) * t717 + t780;
t627 = t757 * t629 + t753 * t632;
t697 = -t725 * mrSges(5,1) + t726 * mrSges(5,2);
t698 = -t743 * mrSges(5,2) + t725 * mrSges(5,3);
t624 = m(5) * t650 + t724 * mrSges(5,1) - t693 * mrSges(5,3) - t726 * t697 + t743 * t698 + t627;
t699 = t743 * mrSges(5,1) - t726 * mrSges(5,3);
t772 = -t753 * t629 + t757 * t632;
t625 = m(5) * t651 - t724 * mrSges(5,2) + t692 * mrSges(5,3) + t725 * t697 - t743 * t699 + t772;
t620 = t758 * t624 + t754 * t625;
t767 = -m(4) * t680 - t730 * mrSges(4,1) - t620;
t618 = m(3) * t700 - t730 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t735 - t736) * qJD(2) + (-t728 - t729) * t746 + t767;
t734 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t746;
t737 = mrSges(4,1) * t746 + qJD(2) * mrSges(4,2);
t665 = -t752 * t762 * pkin(8) + t731 * pkin(3) + qJD(2) * t738 - t679;
t653 = -t692 * pkin(4) - t723 * pkin(9) + t726 * t702 + t665;
t644 = -t657 * pkin(5) - t694 * qJ(6) + t696 * t683 + qJDD(6) + t653;
t779 = m(7) * t644 + t658 * mrSges(7,2) + t696 * t684;
t770 = m(6) * t653 + t658 * mrSges(6,2) + t696 * t685 + t779;
t766 = -m(5) * t665 + t692 * mrSges(5,1) - t693 * mrSges(5,2) + t725 * t698 - t726 * t699 - t770;
t764 = -m(4) * t679 + qJDD(2) * mrSges(4,3) + qJD(2) * t737 + t728 * t784 - t766;
t635 = t764 + m(3) * t701 - qJD(2) * t734 - qJDD(2) * mrSges(3,2) + (mrSges(4,1) + mrSges(3,3)) * t731 + t729 * t784 + t803;
t773 = -t755 * t618 + t759 * t635;
t613 = m(2) * t740 - t762 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t773;
t715 = t771 - t802;
t677 = -t731 * pkin(2) + t765 - t802;
t792 = -t754 * t624 + t758 * t625;
t769 = -m(4) * t677 - t731 * mrSges(4,2) + t737 * t746 - t792;
t763 = -m(3) * t715 + t735 * t784 + t731 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t730 + (-t734 * t755 - t736 * t759) * qJD(1) + t769;
t616 = m(2) * t739 + qJDD(1) * mrSges(2,1) - t762 * mrSges(2,2) + t763;
t793 = t756 * t613 + t760 * t616;
t614 = t759 * t618 + t755 * t635;
t791 = -t795 * t695 - t797 * t696 - t804 * t741;
t790 = t806 * t695 + t799 * t696 + t795 * t741;
t789 = -t799 * t695 - t808 * t696 - t797 * t741;
t787 = t805 * qJD(2) + (t798 * t755 + t796 * t759) * qJD(1);
t786 = -t796 * qJD(2) + (-t800 * t755 - t807 * t759) * qJD(1);
t785 = t798 * qJD(2) + (t809 * t755 + t800 * t759) * qJD(1);
t774 = t760 * t613 - t756 * t616;
t688 = Ifges(5,1) * t726 + Ifges(5,4) * t725 + Ifges(5,5) * t743;
t687 = Ifges(5,4) * t726 + Ifges(5,2) * t725 + Ifges(5,6) * t743;
t686 = Ifges(5,5) * t726 + Ifges(5,6) * t725 + Ifges(5,3) * t743;
t636 = -t658 * mrSges(7,3) - t696 * t675 + t781;
t626 = mrSges(6,2) * t653 + mrSges(7,2) * t644 - mrSges(6,3) * t641 - mrSges(7,3) * t638 - qJ(6) * t636 + t799 * t657 + t808 * t658 - t791 * t695 + t797 * t717 - t790 * t741;
t621 = -mrSges(6,1) * t653 + mrSges(6,3) * t642 - mrSges(7,1) * t644 + mrSges(7,3) * t640 - pkin(5) * (t779 - t794) + qJ(6) * t780 + (-qJ(6) * t684 - t789) * t741 + (-qJ(6) * mrSges(7,2) + t795) * t717 + t791 * t696 + t799 * t658 + (pkin(5) * mrSges(7,1) + t806) * t657;
t619 = -t730 * mrSges(4,3) + t736 * t784 - t769;
t610 = mrSges(5,2) * t665 - mrSges(5,3) * t650 + Ifges(5,1) * t693 + Ifges(5,4) * t692 + Ifges(5,5) * t724 - pkin(9) * t627 - t753 * t621 + t757 * t626 + t725 * t686 - t743 * t687;
t609 = Ifges(5,4) * t693 + Ifges(5,2) * t692 + Ifges(5,6) * t724 - t726 * t686 + t743 * t688 - mrSges(5,1) * t665 + mrSges(5,3) * t651 + t753 * t626 + t757 * t621 - pkin(4) * (t770 + t803) + pkin(9) * t772;
t608 = t786 * qJD(2) + t787 * t784 + mrSges(6,1) * t641 - mrSges(6,2) * t642 + mrSges(7,1) * t638 - mrSges(7,2) * t640 + pkin(5) * t636 + pkin(4) * t627 + pkin(3) * t620 - qJ(3) * t619 + mrSges(3,2) * t715 - mrSges(3,3) * t700 + Ifges(5,6) * t692 + Ifges(5,5) * t693 + mrSges(4,1) * t680 - mrSges(4,3) * t677 + mrSges(5,1) * t650 - mrSges(5,2) * t651 + Ifges(5,3) * t724 - t725 * t688 + t726 * t687 + t789 * t695 + t790 * t696 + t809 * t730 + t804 * t717 + t795 * t657 + t797 * t658 + t798 * qJDD(2) + t800 * t731;
t607 = -mrSges(3,1) * t715 + mrSges(3,3) * t701 - mrSges(4,1) * t679 + mrSges(4,2) * t677 - t754 * t610 - t758 * t609 - pkin(3) * (t766 - t803) - pkin(8) * t792 - pkin(2) * t619 + t807 * t731 + t800 * t730 + t796 * qJDD(2) + t785 * qJD(2) - t787 * t746;
t606 = -pkin(1) * t614 + mrSges(2,3) * t740 - pkin(2) * (-qJD(2) * t736 + t767) - qJ(3) * (-t695 * t682 + t764 - t778 - t794) - t758 * t610 + t754 * t609 + pkin(8) * t620 - mrSges(3,1) * t700 + mrSges(3,2) * t701 - mrSges(4,2) * t680 + mrSges(4,3) * t679 + mrSges(2,1) * g(3) + t762 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-qJ(3) * mrSges(4,1) - t796) * t731 - t798 * t730 + (pkin(2) * mrSges(4,2) - t805) * qJDD(2) + (t785 * t759 + (pkin(2) * t728 + t786) * t755) * qJD(1);
t605 = -mrSges(2,2) * g(3) - mrSges(2,3) * t739 + Ifges(2,5) * qJDD(1) - t762 * Ifges(2,6) - pkin(7) * t614 - t755 * t607 + t759 * t608;
t1 = [-m(1) * g(1) + t774; -m(1) * g(2) + t793; (-m(1) - m(2)) * g(3) + t614; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t793 + t760 * t605 - t756 * t606; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t774 + t756 * t605 + t760 * t606; -mrSges(1,1) * g(2) + mrSges(2,1) * t739 + mrSges(1,2) * g(1) - mrSges(2,2) * t740 + Ifges(2,3) * qJDD(1) + pkin(1) * t763 + pkin(7) * t773 + t759 * t607 + t755 * t608;];
tauB  = t1;
