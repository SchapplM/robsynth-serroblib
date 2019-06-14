% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR12
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-08 00:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR12_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR12_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 00:16:56
% EndTime: 2019-05-08 00:20:11
% DurationCPUTime: 141.95s
% Computational Cost: add. (2275610->417), mult. (5639277->556), div. (0->0), fcn. (4788522->16), ass. (0->181)
t799 = cos(pkin(6));
t792 = qJD(1) * t799 + qJD(2);
t795 = sin(pkin(7));
t798 = cos(pkin(7));
t796 = sin(pkin(6));
t808 = cos(qJ(2));
t830 = qJD(1) * t808;
t826 = t796 * t830;
t777 = (t792 * t795 + t798 * t826) * pkin(10);
t803 = sin(qJ(2));
t832 = qJD(1) * t796;
t845 = pkin(10) * t795;
t781 = (-pkin(2) * t808 - t803 * t845) * t832;
t829 = qJD(1) * qJD(2);
t787 = (qJDD(1) * t803 + t808 * t829) * t796;
t791 = qJDD(1) * t799 + qJDD(2);
t804 = sin(qJ(1));
t809 = cos(qJ(1));
t789 = t804 * g(1) - g(2) * t809;
t810 = qJD(1) ^ 2;
t846 = pkin(9) * t796;
t784 = qJDD(1) * pkin(1) + t810 * t846 + t789;
t790 = -g(1) * t809 - g(2) * t804;
t785 = -pkin(1) * t810 + qJDD(1) * t846 + t790;
t835 = t799 * t808;
t821 = t784 * t835 - t803 * t785;
t831 = qJD(1) * t803;
t844 = pkin(10) * t798;
t738 = -t787 * t844 + t791 * pkin(2) + t792 * t777 + (-g(3) * t808 - t781 * t831) * t796 + t821;
t827 = t796 * t831;
t780 = pkin(2) * t792 - t827 * t844;
t788 = (qJDD(1) * t808 - t803 * t829) * t796;
t818 = t788 * t798 + t791 * t795;
t836 = t799 * t803;
t833 = t784 * t836 + t808 * t785;
t739 = -t792 * t780 + (-g(3) * t803 + t781 * t830) * t796 + t818 * pkin(10) + t833;
t843 = t799 * g(3);
t745 = -t787 * t845 - t788 * pkin(2) - t843 + (-t784 + (-t777 * t808 + t780 * t803) * qJD(1)) * t796;
t802 = sin(qJ(3));
t807 = cos(qJ(3));
t705 = -t802 * t739 + (t738 * t798 + t745 * t795) * t807;
t837 = t798 * t808;
t842 = t795 * t802;
t768 = t792 * t842 + (t802 * t837 + t803 * t807) * t832;
t754 = -t768 * qJD(3) - t802 * t787 + t807 * t818;
t841 = t795 * t807;
t767 = (-t802 * t803 + t807 * t837) * t832 + t792 * t841;
t847 = 2 * qJD(5);
t840 = t796 * t803;
t839 = t796 * t808;
t838 = t798 * t802;
t706 = t738 * t838 + t807 * t739 + t745 * t842;
t756 = -mrSges(4,1) * t767 + mrSges(4,2) * t768;
t778 = t792 * t798 - t795 * t826 + qJD(3);
t762 = mrSges(4,1) * t778 - mrSges(4,3) * t768;
t769 = -t788 * t795 + t791 * t798 + qJDD(3);
t757 = -pkin(3) * t767 - pkin(11) * t768;
t776 = t778 ^ 2;
t697 = -pkin(3) * t776 + pkin(11) * t769 + t757 * t767 + t706;
t720 = -t795 * t738 + t798 * t745;
t755 = t767 * qJD(3) + t807 * t787 + t802 * t818;
t700 = (-t767 * t778 - t755) * pkin(11) + (t768 * t778 - t754) * pkin(3) + t720;
t801 = sin(qJ(4));
t806 = cos(qJ(4));
t689 = -t801 * t697 + t806 * t700;
t759 = -t768 * t801 + t778 * t806;
t723 = qJD(4) * t759 + t755 * t806 + t769 * t801;
t753 = qJDD(4) - t754;
t760 = t768 * t806 + t778 * t801;
t766 = qJD(4) - t767;
t686 = (t759 * t766 - t723) * qJ(5) + (t759 * t760 + t753) * pkin(4) + t689;
t690 = t806 * t697 + t801 * t700;
t722 = -qJD(4) * t760 - t755 * t801 + t769 * t806;
t748 = pkin(4) * t766 - qJ(5) * t760;
t758 = t759 ^ 2;
t688 = -pkin(4) * t758 + qJ(5) * t722 - t748 * t766 + t690;
t794 = sin(pkin(13));
t797 = cos(pkin(13));
t740 = t759 * t797 - t760 * t794;
t683 = t794 * t686 + t797 * t688 + t740 * t847;
t709 = t722 * t797 - t723 * t794;
t741 = t759 * t794 + t760 * t797;
t718 = -mrSges(6,1) * t740 + mrSges(6,2) * t741;
t727 = mrSges(6,1) * t766 - mrSges(6,3) * t741;
t719 = -pkin(5) * t740 - pkin(12) * t741;
t765 = t766 ^ 2;
t681 = -pkin(5) * t765 + pkin(12) * t753 + t719 * t740 + t683;
t696 = -t769 * pkin(3) - t776 * pkin(11) + t768 * t757 - t705;
t691 = -t722 * pkin(4) - t758 * qJ(5) + t760 * t748 + qJDD(5) + t696;
t710 = t722 * t794 + t723 * t797;
t684 = (-t740 * t766 - t710) * pkin(12) + (t741 * t766 - t709) * pkin(5) + t691;
t800 = sin(qJ(6));
t805 = cos(qJ(6));
t678 = -t681 * t800 + t684 * t805;
t724 = -t741 * t800 + t766 * t805;
t694 = qJD(6) * t724 + t710 * t805 + t753 * t800;
t708 = qJDD(6) - t709;
t725 = t741 * t805 + t766 * t800;
t711 = -mrSges(7,1) * t724 + mrSges(7,2) * t725;
t737 = qJD(6) - t740;
t712 = -mrSges(7,2) * t737 + mrSges(7,3) * t724;
t676 = m(7) * t678 + mrSges(7,1) * t708 - mrSges(7,3) * t694 - t711 * t725 + t712 * t737;
t679 = t681 * t805 + t684 * t800;
t693 = -qJD(6) * t725 - t710 * t800 + t753 * t805;
t713 = mrSges(7,1) * t737 - mrSges(7,3) * t725;
t677 = m(7) * t679 - mrSges(7,2) * t708 + mrSges(7,3) * t693 + t711 * t724 - t713 * t737;
t822 = -t676 * t800 + t805 * t677;
t664 = m(6) * t683 - mrSges(6,2) * t753 + mrSges(6,3) * t709 + t718 * t740 - t727 * t766 + t822;
t820 = -t797 * t686 + t794 * t688;
t682 = -0.2e1 * qJD(5) * t741 - t820;
t726 = -mrSges(6,2) * t766 + mrSges(6,3) * t740;
t680 = -t753 * pkin(5) - t765 * pkin(12) + (t847 + t719) * t741 + t820;
t813 = -m(7) * t680 + t693 * mrSges(7,1) - mrSges(7,2) * t694 + t724 * t712 - t713 * t725;
t672 = m(6) * t682 + mrSges(6,1) * t753 - mrSges(6,3) * t710 - t718 * t741 + t726 * t766 + t813;
t661 = t794 * t664 + t797 * t672;
t742 = -mrSges(5,1) * t759 + mrSges(5,2) * t760;
t747 = -mrSges(5,2) * t766 + mrSges(5,3) * t759;
t659 = m(5) * t689 + mrSges(5,1) * t753 - mrSges(5,3) * t723 - t742 * t760 + t747 * t766 + t661;
t749 = mrSges(5,1) * t766 - mrSges(5,3) * t760;
t823 = t797 * t664 - t672 * t794;
t660 = m(5) * t690 - mrSges(5,2) * t753 + mrSges(5,3) * t722 + t742 * t759 - t749 * t766 + t823;
t824 = -t659 * t801 + t806 * t660;
t650 = m(4) * t706 - mrSges(4,2) * t769 + mrSges(4,3) * t754 + t756 * t767 - t762 * t778 + t824;
t653 = t806 * t659 + t801 * t660;
t761 = -mrSges(4,2) * t778 + mrSges(4,3) * t767;
t652 = m(4) * t720 - mrSges(4,1) * t754 + mrSges(4,2) * t755 - t761 * t767 + t762 * t768 + t653;
t668 = t805 * t676 + t800 * t677;
t812 = m(6) * t691 - t709 * mrSges(6,1) + mrSges(6,2) * t710 - t740 * t726 + t727 * t741 + t668;
t811 = -m(5) * t696 + t722 * mrSges(5,1) - mrSges(5,2) * t723 + t759 * t747 - t749 * t760 - t812;
t667 = m(4) * t705 + mrSges(4,1) * t769 - mrSges(4,3) * t755 - t756 * t768 + t761 * t778 + t811;
t639 = t798 * t807 * t667 + t650 * t838 - t652 * t795;
t763 = -g(3) * t839 + t821;
t783 = -mrSges(3,2) * t792 + mrSges(3,3) * t826;
t786 = (-mrSges(3,1) * t808 + mrSges(3,2) * t803) * t832;
t635 = m(3) * t763 + mrSges(3,1) * t791 - mrSges(3,3) * t787 + t783 * t792 - t786 * t827 + t639;
t638 = t650 * t842 + t798 * t652 + t667 * t841;
t773 = -t796 * t784 - t843;
t782 = mrSges(3,1) * t792 - mrSges(3,3) * t827;
t637 = m(3) * t773 - t788 * mrSges(3,1) + t787 * mrSges(3,2) + (t782 * t803 - t783 * t808) * t832 + t638;
t646 = t807 * t650 - t667 * t802;
t764 = -g(3) * t840 + t833;
t645 = m(3) * t764 - mrSges(3,2) * t791 + mrSges(3,3) * t788 - t782 * t792 + t786 * t826 + t646;
t625 = t635 * t835 - t637 * t796 + t645 * t836;
t623 = m(2) * t789 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t810 + t625;
t631 = -t635 * t803 + t808 * t645;
t630 = m(2) * t790 - mrSges(2,1) * t810 - qJDD(1) * mrSges(2,2) + t631;
t834 = t809 * t623 + t804 * t630;
t624 = t635 * t839 + t799 * t637 + t645 * t840;
t825 = -t623 * t804 + t809 * t630;
t701 = Ifges(7,5) * t725 + Ifges(7,6) * t724 + Ifges(7,3) * t737;
t703 = Ifges(7,1) * t725 + Ifges(7,4) * t724 + Ifges(7,5) * t737;
t669 = -mrSges(7,1) * t680 + mrSges(7,3) * t679 + Ifges(7,4) * t694 + Ifges(7,2) * t693 + Ifges(7,6) * t708 - t701 * t725 + t703 * t737;
t702 = Ifges(7,4) * t725 + Ifges(7,2) * t724 + Ifges(7,6) * t737;
t670 = mrSges(7,2) * t680 - mrSges(7,3) * t678 + Ifges(7,1) * t694 + Ifges(7,4) * t693 + Ifges(7,5) * t708 + t701 * t724 - t702 * t737;
t714 = Ifges(6,5) * t741 + Ifges(6,6) * t740 + Ifges(6,3) * t766;
t715 = Ifges(6,4) * t741 + Ifges(6,2) * t740 + Ifges(6,6) * t766;
t654 = mrSges(6,2) * t691 - mrSges(6,3) * t682 + Ifges(6,1) * t710 + Ifges(6,4) * t709 + Ifges(6,5) * t753 - pkin(12) * t668 - t669 * t800 + t670 * t805 + t714 * t740 - t715 * t766;
t716 = Ifges(6,1) * t741 + Ifges(6,4) * t740 + Ifges(6,5) * t766;
t655 = -mrSges(6,1) * t691 - mrSges(7,1) * t678 + mrSges(7,2) * t679 + mrSges(6,3) * t683 + Ifges(6,4) * t710 - Ifges(7,5) * t694 + Ifges(6,2) * t709 + Ifges(6,6) * t753 - Ifges(7,6) * t693 - Ifges(7,3) * t708 - pkin(5) * t668 - t702 * t725 + t703 * t724 - t714 * t741 + t716 * t766;
t728 = Ifges(5,5) * t760 + Ifges(5,6) * t759 + Ifges(5,3) * t766;
t730 = Ifges(5,1) * t760 + Ifges(5,4) * t759 + Ifges(5,5) * t766;
t640 = -mrSges(5,1) * t696 + mrSges(5,3) * t690 + Ifges(5,4) * t723 + Ifges(5,2) * t722 + Ifges(5,6) * t753 - pkin(4) * t812 + qJ(5) * t823 + t794 * t654 + t797 * t655 - t760 * t728 + t766 * t730;
t729 = Ifges(5,4) * t760 + Ifges(5,2) * t759 + Ifges(5,6) * t766;
t641 = mrSges(5,2) * t696 - mrSges(5,3) * t689 + Ifges(5,1) * t723 + Ifges(5,4) * t722 + Ifges(5,5) * t753 - qJ(5) * t661 + t654 * t797 - t655 * t794 + t728 * t759 - t729 * t766;
t751 = Ifges(4,4) * t768 + Ifges(4,2) * t767 + Ifges(4,6) * t778;
t752 = Ifges(4,1) * t768 + Ifges(4,4) * t767 + Ifges(4,5) * t778;
t626 = mrSges(4,1) * t705 - mrSges(4,2) * t706 + Ifges(4,5) * t755 + Ifges(4,6) * t754 + Ifges(4,3) * t769 + pkin(3) * t811 + pkin(11) * t824 + t806 * t640 + t801 * t641 + t768 * t751 - t767 * t752;
t770 = Ifges(3,3) * t792 + (Ifges(3,5) * t803 + Ifges(3,6) * t808) * t832;
t772 = Ifges(3,5) * t792 + (Ifges(3,1) * t803 + Ifges(3,4) * t808) * t832;
t750 = Ifges(4,5) * t768 + Ifges(4,6) * t767 + Ifges(4,3) * t778;
t627 = mrSges(4,2) * t720 - mrSges(4,3) * t705 + Ifges(4,1) * t755 + Ifges(4,4) * t754 + Ifges(4,5) * t769 - pkin(11) * t653 - t640 * t801 + t641 * t806 + t750 * t767 - t751 * t778;
t632 = -pkin(3) * t653 - pkin(12) * t822 - pkin(5) * t813 - pkin(4) * t661 + (-Ifges(5,3) - Ifges(6,3)) * t753 - t805 * t669 - t800 * t670 + t778 * t752 - t768 * t750 + Ifges(4,6) * t769 + t759 * t730 - t760 * t729 + Ifges(4,2) * t754 + Ifges(4,4) * t755 - t741 * t715 + t740 * t716 - mrSges(6,1) * t682 + mrSges(6,2) * t683 - mrSges(5,1) * t689 + mrSges(5,2) * t690 + mrSges(4,3) * t706 - Ifges(6,6) * t709 - Ifges(6,5) * t710 - mrSges(4,1) * t720 - Ifges(5,6) * t722 - Ifges(5,5) * t723;
t814 = pkin(10) * t646 + t627 * t802 + t632 * t807;
t620 = -mrSges(3,1) * t773 + mrSges(3,3) * t764 + Ifges(3,4) * t787 + Ifges(3,2) * t788 + Ifges(3,6) * t791 - pkin(2) * t638 - t795 * t626 - t770 * t827 + t792 * t772 + t798 * t814;
t771 = Ifges(3,6) * t792 + (Ifges(3,4) * t803 + Ifges(3,2) * t808) * t832;
t621 = t770 * t826 + mrSges(3,2) * t773 - mrSges(3,3) * t763 + Ifges(3,1) * t787 + Ifges(3,4) * t788 + Ifges(3,5) * t791 + t807 * t627 - t802 * t632 - t792 * t771 + (-t638 * t795 - t639 * t798) * pkin(10);
t815 = pkin(9) * t631 + t620 * t808 + t621 * t803;
t619 = mrSges(3,1) * t763 - mrSges(3,2) * t764 + Ifges(3,5) * t787 + Ifges(3,6) * t788 + Ifges(3,3) * t791 + pkin(2) * t639 + t798 * t626 + (t771 * t803 - t772 * t808) * t832 + t814 * t795;
t618 = -mrSges(2,2) * g(3) - mrSges(2,3) * t789 + Ifges(2,5) * qJDD(1) - t810 * Ifges(2,6) - t803 * t620 + t808 * t621 + (-t624 * t796 - t625 * t799) * pkin(9);
t617 = mrSges(2,1) * g(3) + mrSges(2,3) * t790 + t810 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t624 - t796 * t619 + t799 * t815;
t1 = [-m(1) * g(1) + t825; -m(1) * g(2) + t834; (-m(1) - m(2)) * g(3) + t624; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t834 - t804 * t617 + t809 * t618; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t825 + t809 * t617 + t804 * t618; -mrSges(1,1) * g(2) + mrSges(2,1) * t789 + mrSges(1,2) * g(1) - mrSges(2,2) * t790 + Ifges(2,3) * qJDD(1) + pkin(1) * t625 + t799 * t619 + t796 * t815;];
tauB  = t1;
