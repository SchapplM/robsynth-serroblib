% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 06:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPPR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:33:11
% EndTime: 2019-05-07 06:33:35
% DurationCPUTime: 19.12s
% Computational Cost: add. (299795->377), mult. (655953->474), div. (0->0), fcn. (513908->12), ass. (0->157)
t826 = -2 * qJD(4);
t825 = Ifges(5,1) + Ifges(6,1);
t815 = Ifges(5,4) - Ifges(6,5);
t824 = -Ifges(5,5) - Ifges(6,4);
t823 = Ifges(5,2) + Ifges(6,3);
t822 = -Ifges(6,2) - Ifges(5,3);
t813 = Ifges(5,6) - Ifges(6,6);
t771 = sin(pkin(6));
t775 = sin(qJ(2));
t778 = cos(qJ(2));
t797 = qJD(1) * qJD(2);
t757 = (-qJDD(1) * t778 + t775 * t797) * t771;
t800 = qJD(1) * t771;
t755 = (-pkin(2) * t778 - pkin(9) * t775) * t800;
t772 = cos(pkin(6));
t767 = qJD(1) * t772 + qJD(2);
t765 = t767 ^ 2;
t766 = qJDD(1) * t772 + qJDD(2);
t799 = qJD(1) * t778;
t776 = sin(qJ(1));
t779 = cos(qJ(1));
t763 = t776 * g(1) - g(2) * t779;
t780 = qJD(1) ^ 2;
t818 = pkin(8) * t771;
t752 = qJDD(1) * pkin(1) + t780 * t818 + t763;
t764 = -g(1) * t779 - g(2) * t776;
t753 = -pkin(1) * t780 + qJDD(1) * t818 + t764;
t808 = t772 * t775;
t801 = t752 * t808 + t778 * t753;
t694 = -t765 * pkin(2) + t766 * pkin(9) + (-g(3) * t775 + t755 * t799) * t771 + t801;
t756 = (qJDD(1) * t775 + t778 * t797) * t771;
t817 = t772 * g(3);
t695 = t757 * pkin(2) - t756 * pkin(9) - t817 + (-t752 + (pkin(2) * t775 - pkin(9) * t778) * t767 * qJD(1)) * t771;
t774 = sin(qJ(3));
t819 = cos(qJ(3));
t670 = t819 * t694 + t774 * t695;
t796 = t775 * t800;
t745 = -t819 * t767 + t774 * t796;
t746 = t774 * t767 + t796 * t819;
t726 = pkin(3) * t745 - qJ(4) * t746;
t749 = qJDD(3) + t757;
t795 = t771 * t799;
t762 = qJD(3) - t795;
t761 = t762 ^ 2;
t662 = -pkin(3) * t761 + qJ(4) * t749 - t726 * t745 + t670;
t807 = t772 * t778;
t809 = t771 * t778;
t724 = -g(3) * t809 + t752 * t807 - t775 * t753;
t693 = -t766 * pkin(2) - t765 * pkin(9) + t755 * t796 - t724;
t722 = qJD(3) * t746 + t756 * t774 - t819 * t766;
t723 = -t745 * qJD(3) + t756 * t819 + t774 * t766;
t664 = (t745 * t762 - t723) * qJ(4) + (t746 * t762 + t722) * pkin(3) + t693;
t770 = sin(pkin(11));
t812 = cos(pkin(11));
t732 = t746 * t812 + t770 * t762;
t657 = -t770 * t662 + t664 * t812 + t732 * t826;
t703 = t723 * t812 + t770 * t749;
t669 = -t774 * t694 + t819 * t695;
t784 = t749 * pkin(3) + t761 * qJ(4) - t746 * t726 - qJDD(4) + t669;
t731 = t770 * t746 - t762 * t812;
t811 = t731 * t745;
t821 = (-t703 + t811) * qJ(5) - t784;
t820 = 2 * qJD(5);
t816 = -mrSges(5,3) - mrSges(6,2);
t810 = t771 * t775;
t725 = -g(3) * t810 + t801;
t750 = mrSges(3,1) * t767 - mrSges(3,3) * t796;
t754 = (-mrSges(3,1) * t778 + mrSges(3,2) * t775) * t800;
t727 = mrSges(4,1) * t745 + mrSges(4,2) * t746;
t734 = mrSges(4,1) * t762 - mrSges(4,3) * t746;
t658 = t812 * t662 + t770 * t664 + t731 * t826;
t702 = t770 * t723 - t749 * t812;
t710 = mrSges(5,1) * t745 - mrSges(5,3) * t732;
t704 = pkin(4) * t731 - qJ(5) * t732;
t744 = t745 ^ 2;
t654 = -pkin(4) * t744 + t722 * qJ(5) - t731 * t704 + t745 * t820 + t658;
t711 = -mrSges(6,1) * t745 + mrSges(6,2) * t732;
t655 = -t722 * pkin(4) - t744 * qJ(5) + t732 * t704 + qJDD(5) - t657;
t649 = (-t703 - t811) * pkin(10) + (t731 * t732 - t722) * pkin(5) + t655;
t712 = -pkin(5) * t745 - pkin(10) * t732;
t730 = t731 ^ 2;
t650 = -pkin(5) * t730 + pkin(10) * t702 + t712 * t745 + t654;
t773 = sin(qJ(6));
t777 = cos(qJ(6));
t647 = t649 * t777 - t650 * t773;
t698 = t731 * t777 - t732 * t773;
t668 = qJD(6) * t698 + t702 * t773 + t703 * t777;
t699 = t731 * t773 + t732 * t777;
t676 = -mrSges(7,1) * t698 + mrSges(7,2) * t699;
t742 = qJD(6) - t745;
t679 = -mrSges(7,2) * t742 + mrSges(7,3) * t698;
t720 = qJDD(6) - t722;
t644 = m(7) * t647 + mrSges(7,1) * t720 - mrSges(7,3) * t668 - t676 * t699 + t679 * t742;
t648 = t649 * t773 + t650 * t777;
t667 = -qJD(6) * t699 + t702 * t777 - t703 * t773;
t680 = mrSges(7,1) * t742 - mrSges(7,3) * t699;
t645 = m(7) * t648 - mrSges(7,2) * t720 + mrSges(7,3) * t667 + t676 * t698 - t680 * t742;
t791 = -t773 * t644 + t777 * t645;
t787 = m(6) * t654 + t722 * mrSges(6,3) + t745 * t711 + t791;
t705 = mrSges(6,1) * t731 - mrSges(6,3) * t732;
t802 = -mrSges(5,1) * t731 - mrSges(5,2) * t732 - t705;
t636 = m(5) * t658 - t722 * mrSges(5,2) + t702 * t816 - t745 * t710 + t731 * t802 + t787;
t709 = -mrSges(5,2) * t745 - mrSges(5,3) * t731;
t638 = t777 * t644 + t773 * t645;
t708 = -mrSges(6,2) * t731 + mrSges(6,3) * t745;
t783 = -m(6) * t655 + t722 * mrSges(6,1) + t745 * t708 - t638;
t637 = m(5) * t657 + t722 * mrSges(5,1) + t703 * t816 + t745 * t709 + t732 * t802 + t783;
t792 = t812 * t636 - t637 * t770;
t633 = m(4) * t670 - mrSges(4,2) * t749 - mrSges(4,3) * t722 - t727 * t745 - t734 * t762 + t792;
t733 = -mrSges(4,2) * t762 - mrSges(4,3) * t745;
t656 = -0.2e1 * qJD(5) * t732 + (t732 * t745 + t702) * pkin(4) + t821;
t652 = -t730 * pkin(10) + (-pkin(4) - pkin(5)) * t702 + (-pkin(4) * t745 + t712 + t820) * t732 - t821;
t788 = -m(7) * t652 + t667 * mrSges(7,1) - t668 * mrSges(7,2) + t698 * t679 - t699 * t680;
t646 = m(6) * t656 + mrSges(6,1) * t702 - t703 * mrSges(6,3) + t708 * t731 - t732 * t711 + t788;
t781 = m(5) * t784 - t702 * mrSges(5,1) - mrSges(5,2) * t703 - t731 * t709 - t710 * t732 - t646;
t642 = m(4) * t669 + mrSges(4,1) * t749 - mrSges(4,3) * t723 - t727 * t746 + t733 * t762 + t781;
t793 = t819 * t633 - t642 * t774;
t623 = m(3) * t725 - mrSges(3,2) * t766 - mrSges(3,3) * t757 - t750 * t767 + t754 * t795 + t793;
t627 = t774 * t633 + t819 * t642;
t738 = -t771 * t752 - t817;
t751 = -mrSges(3,2) * t767 + mrSges(3,3) * t795;
t625 = m(3) * t738 + t757 * mrSges(3,1) + t756 * mrSges(3,2) + (t750 * t775 - t751 * t778) * t800 + t627;
t634 = t770 * t636 + t637 * t812;
t782 = -m(4) * t693 - t722 * mrSges(4,1) - t723 * mrSges(4,2) - t745 * t733 - t746 * t734 - t634;
t630 = m(3) * t724 + t766 * mrSges(3,1) - t756 * mrSges(3,3) + t767 * t751 - t754 * t796 + t782;
t613 = t623 * t808 - t625 * t771 + t630 * t807;
t611 = m(2) * t763 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t780 + t613;
t618 = t778 * t623 - t630 * t775;
t617 = m(2) * t764 - mrSges(2,1) * t780 - qJDD(1) * mrSges(2,2) + t618;
t806 = t779 * t611 + t776 * t617;
t805 = t823 * t731 - t815 * t732 - t813 * t745;
t804 = t813 * t731 + t824 * t732 + t822 * t745;
t803 = -t815 * t731 + t825 * t732 - t824 * t745;
t612 = t623 * t810 + t772 * t625 + t630 * t809;
t794 = -t611 * t776 + t779 * t617;
t671 = Ifges(7,5) * t699 + Ifges(7,6) * t698 + Ifges(7,3) * t742;
t673 = Ifges(7,1) * t699 + Ifges(7,4) * t698 + Ifges(7,5) * t742;
t639 = -mrSges(7,1) * t652 + mrSges(7,3) * t648 + Ifges(7,4) * t668 + Ifges(7,2) * t667 + Ifges(7,6) * t720 - t671 * t699 + t673 * t742;
t672 = Ifges(7,4) * t699 + Ifges(7,2) * t698 + Ifges(7,6) * t742;
t640 = mrSges(7,2) * t652 - mrSges(7,3) * t647 + Ifges(7,1) * t668 + Ifges(7,4) * t667 + Ifges(7,5) * t720 + t671 * t698 - t672 * t742;
t619 = mrSges(5,1) * t784 - mrSges(6,1) * t656 + mrSges(6,2) * t654 + mrSges(5,3) * t658 - pkin(4) * t646 - pkin(5) * t788 - pkin(10) * t791 - t777 * t639 - t773 * t640 - t823 * t702 + t815 * t703 + t813 * t722 + t804 * t732 + t803 * t745;
t626 = -mrSges(5,2) * t784 + mrSges(6,2) * t655 - mrSges(5,3) * t657 - mrSges(6,3) * t656 - pkin(10) * t638 - qJ(5) * t646 - t773 * t639 + t777 * t640 - t815 * t702 + t825 * t703 - t722 * t824 + t804 * t731 + t805 * t745;
t713 = Ifges(4,5) * t746 - Ifges(4,6) * t745 + Ifges(4,3) * t762;
t714 = Ifges(4,4) * t746 - Ifges(4,2) * t745 + Ifges(4,6) * t762;
t609 = mrSges(4,2) * t693 - mrSges(4,3) * t669 + Ifges(4,1) * t723 - Ifges(4,4) * t722 + Ifges(4,5) * t749 - qJ(4) * t634 - t770 * t619 + t626 * t812 - t745 * t713 - t762 * t714;
t715 = Ifges(4,1) * t746 - Ifges(4,4) * t745 + Ifges(4,5) * t762;
t614 = -qJ(5) * t787 + (qJ(5) * t705 - t803) * t731 + (pkin(4) * t705 + t805) * t732 - pkin(4) * t783 + (-Ifges(4,2) + t822) * t722 + (mrSges(6,2) * qJ(5) + t813) * t702 + t762 * t715 - t746 * t713 + Ifges(4,6) * t749 + Ifges(7,3) * t720 + Ifges(4,4) * t723 - t698 * t673 + t699 * t672 - mrSges(4,1) * t693 + Ifges(7,6) * t667 + Ifges(7,5) * t668 + mrSges(4,3) * t670 - mrSges(5,1) * t657 + mrSges(5,2) * t658 - mrSges(6,3) * t654 + mrSges(6,1) * t655 - mrSges(7,2) * t648 + mrSges(7,1) * t647 + pkin(5) * t638 - pkin(3) * t634 + (mrSges(6,2) * pkin(4) + t824) * t703;
t735 = Ifges(3,3) * t767 + (Ifges(3,5) * t775 + Ifges(3,6) * t778) * t800;
t736 = Ifges(3,6) * t767 + (Ifges(3,4) * t775 + Ifges(3,2) * t778) * t800;
t607 = mrSges(3,2) * t738 - mrSges(3,3) * t724 + Ifges(3,1) * t756 - Ifges(3,4) * t757 + Ifges(3,5) * t766 - pkin(9) * t627 + t609 * t819 - t774 * t614 + t735 * t795 - t767 * t736;
t737 = Ifges(3,5) * t767 + (Ifges(3,1) * t775 + Ifges(3,4) * t778) * t800;
t608 = Ifges(3,4) * t756 - Ifges(3,2) * t757 + Ifges(3,6) * t766 - t735 * t796 + t767 * t737 - mrSges(3,1) * t738 + mrSges(3,3) * t725 - Ifges(4,5) * t723 + Ifges(4,6) * t722 - Ifges(4,3) * t749 - t746 * t714 - t745 * t715 - mrSges(4,1) * t669 + mrSges(4,2) * t670 - t770 * t626 - t812 * t619 - pkin(3) * t781 - qJ(4) * t792 - pkin(2) * t627;
t785 = pkin(8) * t618 + t607 * t775 + t608 * t778;
t606 = Ifges(3,5) * t756 - Ifges(3,6) * t757 + Ifges(3,3) * t766 + mrSges(3,1) * t724 - mrSges(3,2) * t725 + t774 * t609 + t819 * t614 + pkin(2) * t782 + pkin(9) * t793 + (t736 * t775 - t737 * t778) * t800;
t605 = -mrSges(2,2) * g(3) - mrSges(2,3) * t763 + Ifges(2,5) * qJDD(1) - t780 * Ifges(2,6) + t778 * t607 - t775 * t608 + (-t612 * t771 - t613 * t772) * pkin(8);
t604 = mrSges(2,1) * g(3) + mrSges(2,3) * t764 + t780 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t612 - t771 * t606 + t772 * t785;
t1 = [-m(1) * g(1) + t794; -m(1) * g(2) + t806; (-m(1) - m(2)) * g(3) + t612; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t806 - t776 * t604 + t779 * t605; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t794 + t779 * t604 + t776 * t605; -mrSges(1,1) * g(2) + mrSges(2,1) * t763 + mrSges(1,2) * g(1) - mrSges(2,2) * t764 + Ifges(2,3) * qJDD(1) + pkin(1) * t613 + t772 * t606 + t771 * t785;];
tauB  = t1;
