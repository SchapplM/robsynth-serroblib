% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 16:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR13_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR13_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:33:14
% EndTime: 2019-05-06 16:33:33
% DurationCPUTime: 17.94s
% Computational Cost: add. (278854->382), mult. (632358->479), div. (0->0), fcn. (465831->12), ass. (0->159)
t835 = -2 * qJD(3);
t834 = Ifges(3,1) + Ifges(4,2);
t827 = Ifges(3,4) + Ifges(4,6);
t826 = Ifges(3,5) - Ifges(4,4);
t833 = Ifges(3,2) + Ifges(4,3);
t825 = Ifges(3,6) - Ifges(4,5);
t832 = Ifges(3,3) + Ifges(4,1);
t784 = cos(pkin(6));
t776 = t784 * qJD(1) + qJD(2);
t787 = sin(qJ(2));
t782 = sin(pkin(6));
t814 = qJD(1) * t782;
t807 = t787 * t814;
t831 = (pkin(2) * t776 + t835) * t807;
t788 = sin(qJ(1));
t792 = cos(qJ(1));
t771 = t788 * g(1) - t792 * g(2);
t793 = qJD(1) ^ 2;
t756 = t793 * t782 * pkin(8) + qJDD(1) * pkin(1) + t771;
t772 = -t792 * g(1) - t788 * g(2);
t811 = qJDD(1) * t782;
t757 = -t793 * pkin(1) + pkin(8) * t811 + t772;
t791 = cos(qJ(2));
t821 = t784 * t787;
t823 = t782 * t787;
t719 = -g(3) * t823 + t756 * t821 + t791 * t757;
t758 = (-pkin(2) * t791 - qJ(3) * t787) * t814;
t774 = t776 ^ 2;
t775 = t784 * qJDD(1) + qJDD(2);
t813 = qJD(1) * t791;
t808 = t782 * t813;
t693 = t774 * pkin(2) - t775 * qJ(3) - t758 * t808 + t776 * t835 - t719;
t830 = -pkin(2) - pkin(9);
t829 = t784 * g(3);
t828 = mrSges(3,1) - mrSges(4,2);
t824 = t782 ^ 2 * t793;
t822 = t782 * t791;
t820 = t784 * t791;
t735 = -t782 * t756 - t829;
t752 = t776 * mrSges(3,1) - mrSges(3,3) * t807;
t753 = -t776 * mrSges(3,2) + mrSges(3,3) * t808;
t755 = mrSges(4,1) * t807 + t776 * mrSges(4,2);
t762 = (qJD(2) * t813 + qJDD(1) * t787) * t782;
t763 = -qJD(2) * t807 + t791 * t811;
t694 = -t763 * pkin(2) + (-t776 * t808 - t762) * qJ(3) + t735 + t831;
t754 = -mrSges(4,1) * t808 - t776 * mrSges(4,3);
t761 = pkin(3) * t807 - t776 * pkin(9);
t810 = t791 ^ 2 * t824;
t687 = -pkin(3) * t810 - t829 - t762 * qJ(3) + t830 * t763 + (-t756 + (-qJ(3) * t776 * t791 - t761 * t787) * qJD(1)) * t782 + t831;
t815 = g(3) * t822 + t787 * t757;
t801 = -t774 * qJ(3) + t758 * t807 + qJDD(3) + t815;
t689 = t762 * pkin(3) + t830 * t775 + (-pkin(3) * t776 * t814 - pkin(9) * t787 * t824 - t756 * t784) * t791 + t801;
t786 = sin(qJ(4));
t790 = cos(qJ(4));
t674 = t790 * t687 + t786 * t689;
t745 = t790 * t776 - t786 * t808;
t716 = t745 * qJD(4) + t790 * t763 + t786 * t775;
t744 = t786 * t776 + t790 * t808;
t721 = t744 * mrSges(5,1) + t745 * mrSges(5,2);
t767 = qJD(4) + t807;
t728 = t767 * mrSges(5,1) - t745 * mrSges(5,3);
t751 = qJDD(4) + t762;
t720 = t744 * pkin(4) - t745 * qJ(5);
t765 = t767 ^ 2;
t669 = -t765 * pkin(4) + t751 * qJ(5) - t744 * t720 + t674;
t686 = t763 * pkin(3) - pkin(9) * t810 + t776 * t761 - t693;
t717 = -t744 * qJD(4) - t786 * t763 + t790 * t775;
t672 = (t744 * t767 - t717) * qJ(5) + (t745 * t767 + t716) * pkin(4) + t686;
t781 = sin(pkin(11));
t783 = cos(pkin(11));
t726 = t783 * t745 + t781 * t767;
t664 = -0.2e1 * qJD(5) * t726 - t781 * t669 + t783 * t672;
t704 = t783 * t717 + t781 * t751;
t725 = -t781 * t745 + t783 * t767;
t662 = (t725 * t744 - t704) * pkin(10) + (t725 * t726 + t716) * pkin(5) + t664;
t665 = 0.2e1 * qJD(5) * t725 + t783 * t669 + t781 * t672;
t703 = -t781 * t717 + t783 * t751;
t709 = t744 * pkin(5) - t726 * pkin(10);
t724 = t725 ^ 2;
t663 = -t724 * pkin(5) + t703 * pkin(10) - t744 * t709 + t665;
t785 = sin(qJ(6));
t789 = cos(qJ(6));
t660 = t789 * t662 - t785 * t663;
t700 = t789 * t725 - t785 * t726;
t677 = t700 * qJD(6) + t785 * t703 + t789 * t704;
t701 = t785 * t725 + t789 * t726;
t682 = -t700 * mrSges(7,1) + t701 * mrSges(7,2);
t743 = qJD(6) + t744;
t690 = -t743 * mrSges(7,2) + t700 * mrSges(7,3);
t714 = qJDD(6) + t716;
t658 = m(7) * t660 + t714 * mrSges(7,1) - t677 * mrSges(7,3) - t701 * t682 + t743 * t690;
t661 = t785 * t662 + t789 * t663;
t676 = -t701 * qJD(6) + t789 * t703 - t785 * t704;
t691 = t743 * mrSges(7,1) - t701 * mrSges(7,3);
t659 = m(7) * t661 - t714 * mrSges(7,2) + t676 * mrSges(7,3) + t700 * t682 - t743 * t691;
t651 = t789 * t658 + t785 * t659;
t705 = -t725 * mrSges(6,1) + t726 * mrSges(6,2);
t707 = -t744 * mrSges(6,2) + t725 * mrSges(6,3);
t649 = m(6) * t664 + t716 * mrSges(6,1) - t704 * mrSges(6,3) - t726 * t705 + t744 * t707 + t651;
t708 = t744 * mrSges(6,1) - t726 * mrSges(6,3);
t803 = -t785 * t658 + t789 * t659;
t650 = m(6) * t665 - t716 * mrSges(6,2) + t703 * mrSges(6,3) + t725 * t705 - t744 * t708 + t803;
t804 = -t781 * t649 + t783 * t650;
t644 = m(5) * t674 - t751 * mrSges(5,2) - t716 * mrSges(5,3) - t744 * t721 - t767 * t728 + t804;
t673 = -t786 * t687 + t790 * t689;
t727 = -t767 * mrSges(5,2) - t744 * mrSges(5,3);
t668 = -t751 * pkin(4) - t765 * qJ(5) + t745 * t720 + qJDD(5) - t673;
t666 = -t703 * pkin(5) - t724 * pkin(10) + t726 * t709 + t668;
t798 = m(7) * t666 - t676 * mrSges(7,1) + t677 * mrSges(7,2) - t700 * t690 + t701 * t691;
t794 = -m(6) * t668 + t703 * mrSges(6,1) - t704 * mrSges(6,2) + t725 * t707 - t726 * t708 - t798;
t654 = m(5) * t673 + t751 * mrSges(5,1) - t717 * mrSges(5,3) - t745 * t721 + t767 * t727 + t794;
t805 = t790 * t644 - t786 * t654;
t802 = m(4) * t694 - t762 * mrSges(4,3) + t754 * t808 + t805;
t633 = m(3) * t735 + t762 * mrSges(3,2) - t828 * t763 + (-t753 * t791 + (t752 - t755) * t787) * t814 + t802;
t809 = t756 * t820;
t718 = t809 - t815;
t759 = (mrSges(4,2) * t791 - mrSges(4,3) * t787) * t814;
t760 = (-mrSges(3,1) * t791 + mrSges(3,2) * t787) * t814;
t636 = t786 * t644 + t790 * t654;
t699 = -t775 * pkin(2) + t801 - t809;
t799 = -m(4) * t699 - t762 * mrSges(4,1) - t636;
t634 = m(3) * t718 - t762 * mrSges(3,3) + (t753 - t754) * t776 + t828 * t775 + (-t759 - t760) * t807 + t799;
t645 = t783 * t649 + t781 * t650;
t796 = m(5) * t686 + t716 * mrSges(5,1) + t717 * mrSges(5,2) + t744 * t727 + t745 * t728 + t645;
t795 = -m(4) * t693 + t775 * mrSges(4,3) + t776 * t755 + t759 * t808 + t796;
t642 = t760 * t808 + t795 + (mrSges(3,3) + mrSges(4,1)) * t763 - t775 * mrSges(3,2) - t776 * t752 + m(3) * t719;
t623 = -t782 * t633 + t634 * t820 + t642 * t821;
t621 = m(2) * t771 + qJDD(1) * mrSges(2,1) - t793 * mrSges(2,2) + t623;
t628 = -t787 * t634 + t791 * t642;
t627 = m(2) * t772 - t793 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t628;
t819 = t792 * t621 + t788 * t627;
t818 = (t826 * t787 + t825 * t791) * t814 + t832 * t776;
t817 = (-t827 * t787 - t833 * t791) * t814 - t825 * t776;
t816 = (t834 * t787 + t827 * t791) * t814 + t826 * t776;
t622 = t784 * t633 + t634 * t822 + t642 * t823;
t806 = -t788 * t621 + t792 * t627;
t678 = Ifges(7,5) * t701 + Ifges(7,6) * t700 + Ifges(7,3) * t743;
t680 = Ifges(7,1) * t701 + Ifges(7,4) * t700 + Ifges(7,5) * t743;
t652 = -mrSges(7,1) * t666 + mrSges(7,3) * t661 + Ifges(7,4) * t677 + Ifges(7,2) * t676 + Ifges(7,6) * t714 - t701 * t678 + t743 * t680;
t679 = Ifges(7,4) * t701 + Ifges(7,2) * t700 + Ifges(7,6) * t743;
t653 = mrSges(7,2) * t666 - mrSges(7,3) * t660 + Ifges(7,1) * t677 + Ifges(7,4) * t676 + Ifges(7,5) * t714 + t700 * t678 - t743 * t679;
t695 = Ifges(6,5) * t726 + Ifges(6,6) * t725 + Ifges(6,3) * t744;
t697 = Ifges(6,1) * t726 + Ifges(6,4) * t725 + Ifges(6,5) * t744;
t637 = -mrSges(6,1) * t668 + mrSges(6,3) * t665 + Ifges(6,4) * t704 + Ifges(6,2) * t703 + Ifges(6,6) * t716 - pkin(5) * t798 + pkin(10) * t803 + t789 * t652 + t785 * t653 - t726 * t695 + t744 * t697;
t696 = Ifges(6,4) * t726 + Ifges(6,2) * t725 + Ifges(6,6) * t744;
t638 = mrSges(6,2) * t668 - mrSges(6,3) * t664 + Ifges(6,1) * t704 + Ifges(6,4) * t703 + Ifges(6,5) * t716 - pkin(10) * t651 - t785 * t652 + t789 * t653 + t725 * t695 - t744 * t696;
t710 = Ifges(5,5) * t745 - Ifges(5,6) * t744 + Ifges(5,3) * t767;
t711 = Ifges(5,4) * t745 - Ifges(5,2) * t744 + Ifges(5,6) * t767;
t624 = mrSges(5,2) * t686 - mrSges(5,3) * t673 + Ifges(5,1) * t717 - Ifges(5,4) * t716 + Ifges(5,5) * t751 - qJ(5) * t645 - t781 * t637 + t783 * t638 - t744 * t710 - t767 * t711;
t712 = Ifges(5,1) * t745 - Ifges(5,4) * t744 + Ifges(5,5) * t767;
t629 = Ifges(5,4) * t717 + Ifges(5,6) * t751 - t745 * t710 + t767 * t712 - mrSges(5,1) * t686 + mrSges(5,3) * t674 - Ifges(6,5) * t704 - Ifges(6,6) * t703 - t726 * t696 + t725 * t697 - mrSges(6,1) * t664 + mrSges(6,2) * t665 - Ifges(7,5) * t677 - Ifges(7,6) * t676 - Ifges(7,3) * t714 - t701 * t679 + t700 * t680 - mrSges(7,1) * t660 + mrSges(7,2) * t661 - pkin(5) * t651 - pkin(4) * t645 + (-Ifges(5,2) - Ifges(6,3)) * t716;
t635 = t763 * mrSges(4,2) - t755 * t807 + t802;
t618 = -mrSges(3,1) * t735 - mrSges(4,1) * t693 + mrSges(4,2) * t694 + mrSges(3,3) * t719 - pkin(2) * t635 + pkin(3) * t796 - pkin(9) * t805 - t786 * t624 - t790 * t629 + t827 * t762 + t833 * t763 + t825 * t775 + t816 * t776 - t818 * t807;
t619 = pkin(3) * t636 - qJ(3) * t635 + qJ(5) * t804 + t781 * t638 + t783 * t637 + Ifges(5,3) * t751 + t744 * t712 + t745 * t711 + mrSges(3,2) * t735 + pkin(4) * t794 - Ifges(5,6) * t716 + Ifges(5,5) * t717 - mrSges(3,3) * t718 - mrSges(4,3) * t694 + mrSges(4,1) * t699 + mrSges(5,1) * t673 - mrSges(5,2) * t674 + t817 * t776 + t826 * t775 + t827 * t763 + t834 * t762 + t818 * t808;
t800 = pkin(8) * t628 + t618 * t791 + t619 * t787;
t617 = mrSges(3,1) * t718 - mrSges(3,2) * t719 + mrSges(4,2) * t699 - mrSges(4,3) * t693 + t790 * t624 - t786 * t629 - pkin(9) * t636 + pkin(2) * (-t776 * t754 + t799) + qJ(3) * t795 + (-pkin(2) * mrSges(4,2) + t832) * t775 + (qJ(3) * mrSges(4,1) + t825) * t763 + t826 * t762 + (-t816 * t791 + (-pkin(2) * t759 - t817) * t787) * t814;
t616 = -mrSges(2,2) * g(3) - mrSges(2,3) * t771 + Ifges(2,5) * qJDD(1) - t793 * Ifges(2,6) - t787 * t618 + t791 * t619 + (-t622 * t782 - t623 * t784) * pkin(8);
t615 = mrSges(2,1) * g(3) + mrSges(2,3) * t772 + t793 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t622 - t782 * t617 + t800 * t784;
t1 = [-m(1) * g(1) + t806; -m(1) * g(2) + t819; (-m(1) - m(2)) * g(3) + t622; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t819 - t788 * t615 + t792 * t616; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t806 + t792 * t615 + t788 * t616; -mrSges(1,1) * g(2) + mrSges(2,1) * t771 + mrSges(1,2) * g(1) - mrSges(2,2) * t772 + Ifges(2,3) * qJDD(1) + pkin(1) * t623 + t784 * t617 + t800 * t782;];
tauB  = t1;
