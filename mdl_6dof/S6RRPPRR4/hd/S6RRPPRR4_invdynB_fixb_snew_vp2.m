% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 10:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:20:40
% EndTime: 2019-05-06 10:20:57
% DurationCPUTime: 16.21s
% Computational Cost: add. (228567->378), mult. (614723->473), div. (0->0), fcn. (470114->12), ass. (0->160)
t837 = -2 * qJD(4);
t836 = Ifges(4,1) + Ifges(5,2);
t835 = -Ifges(5,1) - Ifges(4,3);
t830 = Ifges(4,5) - Ifges(5,4);
t834 = -Ifges(4,2) - Ifges(5,3);
t829 = Ifges(4,6) - Ifges(5,5);
t828 = -Ifges(5,6) - Ifges(4,4);
t780 = sin(pkin(11));
t789 = cos(qJ(2));
t781 = sin(pkin(6));
t814 = qJD(1) * t781;
t807 = t789 * t814;
t785 = sin(qJ(2));
t808 = t785 * t814;
t827 = cos(pkin(11));
t756 = t780 * t808 - t827 * t807;
t757 = (t780 * t789 + t827 * t785) * t814;
t723 = t756 * pkin(3) - t757 * qJ(4);
t782 = cos(pkin(6));
t776 = t782 * qJD(1) + qJD(2);
t774 = t776 ^ 2;
t775 = t782 * qJDD(1) + qJDD(2);
t810 = qJD(1) * qJD(2);
t766 = (qJDD(1) * t785 + t789 * t810) * t781;
t786 = sin(qJ(1));
t790 = cos(qJ(1));
t772 = t786 * g(1) - t790 * g(2);
t791 = qJD(1) ^ 2;
t832 = pkin(8) * t781;
t763 = qJDD(1) * pkin(1) + t791 * t832 + t772;
t773 = -t790 * g(1) - t786 * g(2);
t764 = -t791 * pkin(1) + qJDD(1) * t832 + t773;
t821 = t782 * t789;
t802 = t763 * t821 - t785 * t764;
t825 = t781 ^ 2 * t791;
t690 = t775 * pkin(2) - t766 * qJ(3) + (pkin(2) * t785 * t825 + (qJ(3) * qJD(1) * t776 - g(3)) * t781) * t789 + t802;
t822 = t782 * t785;
t824 = t781 * t785;
t722 = -g(3) * t824 + t763 * t822 + t789 * t764;
t760 = t776 * pkin(2) - qJ(3) * t808;
t767 = (qJDD(1) * t789 - t785 * t810) * t781;
t809 = t789 ^ 2 * t825;
t693 = -pkin(2) * t809 + t767 * qJ(3) - t776 * t760 + t722;
t819 = t780 * t690 + t827 * t693;
t833 = t774 * pkin(3) - t775 * qJ(4) + t756 * t723 + t776 * t837 - t819;
t677 = -0.2e1 * qJD(3) * t757 + t827 * t690 - t780 * t693;
t831 = mrSges(4,1) - mrSges(5,2);
t826 = t756 * t776;
t823 = t781 * t789;
t724 = t756 * mrSges(4,1) + t757 * mrSges(4,2);
t731 = t827 * t766 + t780 * t767;
t674 = -t775 * pkin(3) - t774 * qJ(4) + t757 * t723 + qJDD(4) - t677;
t668 = (t756 * t757 - t775) * pkin(9) + (t731 + t826) * pkin(4) + t674;
t730 = t780 * t766 - t827 * t767;
t741 = t757 * pkin(4) - t776 * pkin(9);
t755 = t756 ^ 2;
t746 = -t782 * g(3) - t781 * t763;
t703 = -t767 * pkin(2) - qJ(3) * t809 + t760 * t808 + qJDD(3) + t746;
t792 = (-t731 + t826) * qJ(4) + t703 + (t776 * pkin(3) + t837) * t757;
t672 = -t755 * pkin(4) - t757 * t741 + (pkin(3) + pkin(9)) * t730 + t792;
t784 = sin(qJ(5));
t788 = cos(qJ(5));
t665 = t784 * t668 + t788 * t672;
t736 = t784 * t756 + t788 * t776;
t700 = -t736 * qJD(5) + t788 * t730 - t784 * t775;
t735 = t788 * t756 - t784 * t776;
t704 = -t735 * mrSges(6,1) + t736 * mrSges(6,2);
t754 = qJD(5) + t757;
t712 = t754 * mrSges(6,1) - t736 * mrSges(6,3);
t729 = qJDD(5) + t731;
t705 = -t735 * pkin(5) - t736 * pkin(10);
t753 = t754 ^ 2;
t663 = -t753 * pkin(5) + t729 * pkin(10) + t735 * t705 + t665;
t813 = qJD(3) * t756;
t749 = -0.2e1 * t813;
t670 = -t730 * pkin(4) - t755 * pkin(9) + t776 * t741 + t749 - t833;
t701 = t735 * qJD(5) + t784 * t730 + t788 * t775;
t666 = (-t735 * t754 - t701) * pkin(10) + (t736 * t754 - t700) * pkin(5) + t670;
t783 = sin(qJ(6));
t787 = cos(qJ(6));
t660 = -t783 * t663 + t787 * t666;
t709 = -t783 * t736 + t787 * t754;
t681 = t709 * qJD(6) + t787 * t701 + t783 * t729;
t710 = t787 * t736 + t783 * t754;
t686 = -t709 * mrSges(7,1) + t710 * mrSges(7,2);
t734 = qJD(6) - t735;
t691 = -t734 * mrSges(7,2) + t709 * mrSges(7,3);
t699 = qJDD(6) - t700;
t658 = m(7) * t660 + t699 * mrSges(7,1) - t681 * mrSges(7,3) - t710 * t686 + t734 * t691;
t661 = t787 * t663 + t783 * t666;
t680 = -t710 * qJD(6) - t783 * t701 + t787 * t729;
t692 = t734 * mrSges(7,1) - t710 * mrSges(7,3);
t659 = m(7) * t661 - t699 * mrSges(7,2) + t680 * mrSges(7,3) + t709 * t686 - t734 * t692;
t803 = -t783 * t658 + t787 * t659;
t650 = m(6) * t665 - t729 * mrSges(6,2) + t700 * mrSges(6,3) + t735 * t704 - t754 * t712 + t803;
t664 = t788 * t668 - t784 * t672;
t711 = -t754 * mrSges(6,2) + t735 * mrSges(6,3);
t662 = -t729 * pkin(5) - t753 * pkin(10) + t736 * t705 - t664;
t797 = -m(7) * t662 + t680 * mrSges(7,1) - t681 * mrSges(7,2) + t709 * t691 - t710 * t692;
t654 = m(6) * t664 + t729 * mrSges(6,1) - t701 * mrSges(6,3) - t736 * t704 + t754 * t711 + t797;
t645 = t784 * t650 + t788 * t654;
t725 = -t756 * mrSges(5,2) - t757 * mrSges(5,3);
t796 = -m(5) * t674 - t731 * mrSges(5,1) - t757 * t725 - t645;
t739 = t756 * mrSges(5,1) - t776 * mrSges(5,3);
t815 = -t776 * mrSges(4,2) - t756 * mrSges(4,3) - t739;
t641 = m(4) * t677 - t731 * mrSges(4,3) - t757 * t724 + t831 * t775 + t815 * t776 + t796;
t678 = t749 + t819;
t738 = t776 * mrSges(4,1) - t757 * mrSges(4,3);
t673 = 0.2e1 * t813 + t833;
t740 = t757 * mrSges(5,1) + t776 * mrSges(5,2);
t651 = t787 * t658 + t783 * t659;
t795 = -m(6) * t670 + t700 * mrSges(6,1) - t701 * mrSges(6,2) + t735 * t711 - t736 * t712 - t651;
t794 = -m(5) * t673 + t775 * mrSges(5,3) + t776 * t740 - t795;
t648 = t794 + (-mrSges(4,3) - mrSges(5,1)) * t730 + (-t724 - t725) * t756 - t775 * mrSges(4,2) - t776 * t738 + m(4) * t678;
t637 = t827 * t641 + t780 * t648;
t721 = -g(3) * t823 + t802;
t762 = -t776 * mrSges(3,2) + mrSges(3,3) * t807;
t765 = (-mrSges(3,1) * t789 + mrSges(3,2) * t785) * t814;
t635 = m(3) * t721 + t775 * mrSges(3,1) - t766 * mrSges(3,3) + t776 * t762 - t765 * t808 + t637;
t761 = t776 * mrSges(3,1) - mrSges(3,3) * t808;
t805 = -t780 * t641 + t827 * t648;
t636 = m(3) * t722 - t775 * mrSges(3,2) + t767 * mrSges(3,3) - t776 * t761 + t765 * t807 + t805;
t676 = t730 * pkin(3) + t792;
t804 = t788 * t650 - t784 * t654;
t799 = m(5) * t676 - t731 * mrSges(5,3) - t757 * t740 + t804;
t793 = m(4) * t703 + t731 * mrSges(4,2) + t831 * t730 + t757 * t738 + t815 * t756 + t799;
t643 = t793 + (t761 * t785 - t762 * t789) * t814 + t766 * mrSges(3,2) - t767 * mrSges(3,1) + m(3) * t746;
t624 = t635 * t821 + t636 * t822 - t781 * t643;
t622 = m(2) * t772 + qJDD(1) * mrSges(2,1) - t791 * mrSges(2,2) + t624;
t629 = -t785 * t635 + t789 * t636;
t628 = m(2) * t773 - t791 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t629;
t820 = t790 * t622 + t786 * t628;
t818 = t829 * t756 - t830 * t757 + t835 * t776;
t817 = t834 * t756 - t828 * t757 + t829 * t776;
t816 = t828 * t756 + t836 * t757 + t830 * t776;
t623 = t635 * t823 + t636 * t824 + t782 * t643;
t806 = -t786 * t622 + t790 * t628;
t682 = Ifges(7,5) * t710 + Ifges(7,6) * t709 + Ifges(7,3) * t734;
t684 = Ifges(7,1) * t710 + Ifges(7,4) * t709 + Ifges(7,5) * t734;
t652 = -mrSges(7,1) * t662 + mrSges(7,3) * t661 + Ifges(7,4) * t681 + Ifges(7,2) * t680 + Ifges(7,6) * t699 - t710 * t682 + t734 * t684;
t683 = Ifges(7,4) * t710 + Ifges(7,2) * t709 + Ifges(7,6) * t734;
t653 = mrSges(7,2) * t662 - mrSges(7,3) * t660 + Ifges(7,1) * t681 + Ifges(7,4) * t680 + Ifges(7,5) * t699 + t709 * t682 - t734 * t683;
t694 = Ifges(6,5) * t736 + Ifges(6,6) * t735 + Ifges(6,3) * t754;
t695 = Ifges(6,4) * t736 + Ifges(6,2) * t735 + Ifges(6,6) * t754;
t638 = mrSges(6,2) * t670 - mrSges(6,3) * t664 + Ifges(6,1) * t701 + Ifges(6,4) * t700 + Ifges(6,5) * t729 - pkin(10) * t651 - t783 * t652 + t787 * t653 + t735 * t694 - t754 * t695;
t696 = Ifges(6,1) * t736 + Ifges(6,4) * t735 + Ifges(6,5) * t754;
t639 = -mrSges(6,1) * t670 - mrSges(7,1) * t660 + mrSges(7,2) * t661 + mrSges(6,3) * t665 + Ifges(6,4) * t701 - Ifges(7,5) * t681 + Ifges(6,2) * t700 + Ifges(6,6) * t729 - Ifges(7,6) * t680 - Ifges(7,3) * t699 - pkin(5) * t651 - t710 * t683 + t709 * t684 - t736 * t694 + t754 * t696;
t644 = -t730 * mrSges(5,2) - t756 * t739 + t799;
t620 = -mrSges(4,1) * t703 - mrSges(5,1) * t673 + mrSges(5,2) * t676 + mrSges(4,3) * t678 - pkin(3) * t644 - pkin(4) * t795 - pkin(9) * t804 - t784 * t638 - t788 * t639 + t834 * t730 - t828 * t731 + t818 * t757 + t829 * t775 + t816 * t776;
t625 = -qJ(4) * t644 + pkin(4) * t645 + t787 * t652 + pkin(10) * t803 + t783 * t653 - t735 * t696 + t736 * t695 + Ifges(6,3) * t729 + pkin(5) * t797 + mrSges(4,2) * t703 + Ifges(6,6) * t700 + Ifges(6,5) * t701 - mrSges(4,3) * t677 + mrSges(5,1) * t674 - mrSges(5,3) * t676 - mrSges(6,2) * t665 + mrSges(6,1) * t664 - t817 * t776 + t830 * t775 + t818 * t756 + t836 * t731 + t828 * t730;
t743 = Ifges(3,3) * t776 + (Ifges(3,5) * t785 + Ifges(3,6) * t789) * t814;
t745 = Ifges(3,5) * t776 + (Ifges(3,1) * t785 + Ifges(3,4) * t789) * t814;
t617 = -mrSges(3,1) * t746 + mrSges(3,3) * t722 + Ifges(3,4) * t766 + Ifges(3,2) * t767 + Ifges(3,6) * t775 - pkin(2) * t793 + qJ(3) * t805 + t827 * t620 + t780 * t625 - t743 * t808 + t776 * t745;
t744 = Ifges(3,6) * t776 + (Ifges(3,4) * t785 + Ifges(3,2) * t789) * t814;
t618 = mrSges(3,2) * t746 - mrSges(3,3) * t721 + Ifges(3,1) * t766 + Ifges(3,4) * t767 + Ifges(3,5) * t775 - qJ(3) * t637 - t780 * t620 + t827 * t625 + t743 * t807 - t776 * t744;
t798 = pkin(8) * t629 + t617 * t789 + t618 * t785;
t619 = -pkin(9) * t645 + t788 * t638 + pkin(3) * (-t776 * t739 + t796) - t784 * t639 + Ifges(3,5) * t766 + Ifges(3,6) * t767 + qJ(4) * t794 + pkin(2) * t637 + mrSges(3,1) * t721 - mrSges(3,2) * t722 + mrSges(4,1) * t677 - mrSges(4,2) * t678 + mrSges(5,2) * t674 - mrSges(5,3) * t673 + t817 * t757 + (-qJ(4) * t725 + t816) * t756 + t830 * t731 + (-qJ(4) * mrSges(5,1) - t829) * t730 + (t744 * t785 - t745 * t789) * t814 + (-pkin(3) * mrSges(5,2) + Ifges(3,3) - t835) * t775;
t616 = -mrSges(2,2) * g(3) - mrSges(2,3) * t772 + Ifges(2,5) * qJDD(1) - t791 * Ifges(2,6) - t785 * t617 + t789 * t618 + (-t623 * t781 - t624 * t782) * pkin(8);
t615 = mrSges(2,1) * g(3) + mrSges(2,3) * t773 + t791 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t623 - t781 * t619 + t798 * t782;
t1 = [-m(1) * g(1) + t806; -m(1) * g(2) + t820; (-m(1) - m(2)) * g(3) + t623; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t820 - t786 * t615 + t790 * t616; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t806 + t790 * t615 + t786 * t616; -mrSges(1,1) * g(2) + mrSges(2,1) * t772 + mrSges(1,2) * g(1) - mrSges(2,2) * t773 + Ifges(2,3) * qJDD(1) + pkin(1) * t624 + t782 * t619 + t798 * t781;];
tauB  = t1;
