% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 05:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR11_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR11_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 05:34:50
% EndTime: 2019-05-06 05:36:41
% DurationCPUTime: 115.61s
% Computational Cost: add. (1807094->403), mult. (5603668->550), div. (0->0), fcn. (4831991->16), ass. (0->184)
t809 = sin(pkin(13));
t811 = sin(pkin(6));
t812 = cos(pkin(13));
t814 = cos(pkin(6));
t818 = sin(qJ(3));
t813 = cos(pkin(7));
t823 = cos(qJ(3));
t857 = t813 * t823;
t810 = sin(pkin(7));
t862 = t810 * t823;
t828 = t811 * (-t809 * t818 + t812 * t857) + t814 * t862;
t784 = t828 * qJD(1);
t858 = t813 * t818;
t863 = t810 * t818;
t830 = t814 * t863 + (t809 * t823 + t812 * t858) * t811;
t785 = t830 * qJD(1);
t770 = -t785 * qJD(3) + qJDD(1) * t828;
t860 = t811 * t813;
t796 = (t810 * t814 + t812 * t860) * qJD(1) * pkin(9);
t819 = sin(qJ(1));
t824 = cos(qJ(1));
t807 = -g(1) * t824 - g(2) * t819;
t825 = qJD(1) ^ 2;
t866 = qJ(2) * t811;
t800 = -pkin(1) * t825 + qJDD(1) * t866 + t807;
t869 = pkin(9) * t809;
t841 = -pkin(2) * t812 - t810 * t869;
t855 = qJD(1) * t811;
t867 = pkin(9) * qJDD(1);
t836 = qJD(1) * t841 * t855 + t813 * t867;
t806 = t819 * g(1) - g(2) * t824;
t799 = qJDD(1) * pkin(1) + t825 * t866 + t806;
t851 = qJD(2) * t855;
t859 = t812 * t814;
t861 = t811 * t812;
t842 = -g(3) * t861 + t799 * t859 - 0.2e1 * t809 * t851;
t750 = (pkin(2) * qJDD(1) + qJD(1) * t796) * t814 + (-t811 * t836 - t800) * t809 + t842;
t801 = (pkin(2) * t814 - t860 * t869) * qJD(1);
t864 = t809 * t814;
t852 = t799 * t864 + (t800 + 0.2e1 * t851) * t812;
t751 = (-qJD(1) * t801 + t810 * t867) * t814 + (-g(3) * t809 + t812 * t836) * t811 + t852;
t850 = -t814 * g(3) + qJDD(2);
t760 = (-t799 + t841 * qJDD(1) + (-t796 * t812 + t801 * t809) * qJD(1)) * t811 + t850;
t719 = -t818 * t751 + (t750 * t813 + t760 * t810) * t823;
t868 = Ifges(3,3) * t814;
t865 = t809 * t811;
t720 = t750 * t858 + t823 * t751 + t760 * t863;
t768 = -mrSges(4,1) * t784 + mrSges(4,2) * t785;
t837 = -t810 * t861 + t813 * t814;
t797 = qJD(1) * t837 + qJD(3);
t780 = mrSges(4,1) * t797 - mrSges(4,3) * t785;
t794 = qJDD(1) * t837 + qJDD(3);
t769 = -pkin(3) * t784 - pkin(10) * t785;
t793 = t797 ^ 2;
t711 = -pkin(3) * t793 + pkin(10) * t794 + t769 * t784 + t720;
t730 = -t810 * t750 + t813 * t760;
t771 = t784 * qJD(3) + qJDD(1) * t830;
t716 = (-t784 * t797 - t771) * pkin(10) + (t785 * t797 - t770) * pkin(3) + t730;
t817 = sin(qJ(4));
t822 = cos(qJ(4));
t704 = t822 * t711 + t817 * t716;
t778 = t785 * t822 + t797 * t817;
t745 = -t778 * qJD(4) - t817 * t771 + t794 * t822;
t777 = -t817 * t785 + t797 * t822;
t752 = -mrSges(5,1) * t777 + mrSges(5,2) * t778;
t783 = qJD(4) - t784;
t762 = mrSges(5,1) * t783 - mrSges(5,3) * t778;
t767 = qJDD(4) - t770;
t753 = -pkin(4) * t777 - pkin(11) * t778;
t782 = t783 ^ 2;
t699 = -pkin(4) * t782 + pkin(11) * t767 + t753 * t777 + t704;
t710 = -t794 * pkin(3) - t793 * pkin(10) + t785 * t769 - t719;
t746 = qJD(4) * t777 + t771 * t822 + t794 * t817;
t702 = (-t777 * t783 - t746) * pkin(11) + (t778 * t783 - t745) * pkin(4) + t710;
t816 = sin(qJ(5));
t821 = cos(qJ(5));
t694 = -t816 * t699 + t821 * t702;
t758 = -t778 * t816 + t783 * t821;
t723 = qJD(5) * t758 + t746 * t821 + t767 * t816;
t743 = qJDD(5) - t745;
t759 = t778 * t821 + t783 * t816;
t774 = qJD(5) - t777;
t692 = (t758 * t774 - t723) * pkin(12) + (t758 * t759 + t743) * pkin(5) + t694;
t695 = t821 * t699 + t816 * t702;
t722 = -qJD(5) * t759 - t746 * t816 + t767 * t821;
t737 = pkin(5) * t774 - pkin(12) * t759;
t757 = t758 ^ 2;
t693 = -pkin(5) * t757 + pkin(12) * t722 - t737 * t774 + t695;
t815 = sin(qJ(6));
t820 = cos(qJ(6));
t690 = t692 * t820 - t693 * t815;
t731 = t758 * t820 - t759 * t815;
t707 = qJD(6) * t731 + t722 * t815 + t723 * t820;
t732 = t758 * t815 + t759 * t820;
t718 = -mrSges(7,1) * t731 + mrSges(7,2) * t732;
t772 = qJD(6) + t774;
t724 = -mrSges(7,2) * t772 + mrSges(7,3) * t731;
t738 = qJDD(6) + t743;
t688 = m(7) * t690 + mrSges(7,1) * t738 - mrSges(7,3) * t707 - t718 * t732 + t724 * t772;
t691 = t692 * t815 + t693 * t820;
t706 = -qJD(6) * t732 + t722 * t820 - t723 * t815;
t725 = mrSges(7,1) * t772 - mrSges(7,3) * t732;
t689 = m(7) * t691 - mrSges(7,2) * t738 + mrSges(7,3) * t706 + t718 * t731 - t725 * t772;
t680 = t820 * t688 + t815 * t689;
t733 = -mrSges(6,1) * t758 + mrSges(6,2) * t759;
t735 = -mrSges(6,2) * t774 + mrSges(6,3) * t758;
t678 = m(6) * t694 + mrSges(6,1) * t743 - mrSges(6,3) * t723 - t733 * t759 + t735 * t774 + t680;
t736 = mrSges(6,1) * t774 - mrSges(6,3) * t759;
t846 = -t688 * t815 + t820 * t689;
t679 = m(6) * t695 - mrSges(6,2) * t743 + mrSges(6,3) * t722 + t733 * t758 - t736 * t774 + t846;
t847 = -t678 * t816 + t821 * t679;
t675 = m(5) * t704 - mrSges(5,2) * t767 + mrSges(5,3) * t745 + t752 * t777 - t762 * t783 + t847;
t703 = -t817 * t711 + t716 * t822;
t761 = -mrSges(5,2) * t783 + mrSges(5,3) * t777;
t698 = -pkin(4) * t767 - pkin(11) * t782 + t778 * t753 - t703;
t696 = -pkin(5) * t722 - pkin(12) * t757 + t737 * t759 + t698;
t833 = m(7) * t696 - t706 * mrSges(7,1) + mrSges(7,2) * t707 - t731 * t724 + t725 * t732;
t826 = -m(6) * t698 + t722 * mrSges(6,1) - mrSges(6,2) * t723 + t758 * t735 - t736 * t759 - t833;
t684 = m(5) * t703 + mrSges(5,1) * t767 - mrSges(5,3) * t746 - t752 * t778 + t761 * t783 + t826;
t848 = t822 * t675 - t684 * t817;
t664 = m(4) * t720 - mrSges(4,2) * t794 + mrSges(4,3) * t770 + t768 * t784 - t780 * t797 + t848;
t667 = t817 * t675 + t822 * t684;
t779 = -mrSges(4,2) * t797 + mrSges(4,3) * t784;
t666 = m(4) * t730 - mrSges(4,1) * t770 + mrSges(4,2) * t771 - t779 * t784 + t780 * t785 + t667;
t676 = t678 * t821 + t679 * t816;
t827 = -m(5) * t710 + t745 * mrSges(5,1) - mrSges(5,2) * t746 + t777 * t761 - t762 * t778 - t676;
t672 = m(4) * t719 + mrSges(4,1) * t794 - mrSges(4,3) * t771 - t768 * t785 + t779 * t797 + t827;
t653 = t664 * t858 - t810 * t666 + t672 * t857;
t775 = -t809 * t800 + t842;
t845 = -mrSges(3,1) * t812 + mrSges(3,2) * t809;
t798 = t845 * t855;
t839 = -mrSges(3,2) * t814 + mrSges(3,3) * t861;
t803 = t839 * qJD(1);
t840 = mrSges(3,1) * t814 - mrSges(3,3) * t865;
t649 = m(3) * t775 + t840 * qJDD(1) + (-t798 * t865 + t803 * t814) * qJD(1) + t653;
t652 = t664 * t863 + t813 * t666 + t672 * t862;
t786 = -t811 * t799 + t850;
t802 = t840 * qJD(1);
t651 = m(3) * t786 + (t845 * qJDD(1) + (t802 * t809 - t803 * t812) * qJD(1)) * t811 + t652;
t659 = t823 * t664 - t818 * t672;
t776 = -g(3) * t865 + t852;
t658 = m(3) * t776 + t839 * qJDD(1) + (t798 * t861 - t802 * t814) * qJD(1) + t659;
t639 = t649 * t859 - t651 * t811 + t658 * t864;
t637 = m(2) * t806 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t825 + t639;
t645 = -t649 * t809 + t812 * t658;
t644 = m(2) * t807 - mrSges(2,1) * t825 - qJDD(1) * mrSges(2,2) + t645;
t856 = t824 * t637 + t819 * t644;
t638 = t649 * t861 + t814 * t651 + t658 * t865;
t849 = -t637 * t819 + t824 * t644;
t844 = Ifges(3,5) * t809 + Ifges(3,6) * t812;
t713 = Ifges(7,5) * t732 + Ifges(7,6) * t731 + Ifges(7,3) * t772;
t715 = Ifges(7,1) * t732 + Ifges(7,4) * t731 + Ifges(7,5) * t772;
t681 = -mrSges(7,1) * t696 + mrSges(7,3) * t691 + Ifges(7,4) * t707 + Ifges(7,2) * t706 + Ifges(7,6) * t738 - t713 * t732 + t715 * t772;
t714 = Ifges(7,4) * t732 + Ifges(7,2) * t731 + Ifges(7,6) * t772;
t682 = mrSges(7,2) * t696 - mrSges(7,3) * t690 + Ifges(7,1) * t707 + Ifges(7,4) * t706 + Ifges(7,5) * t738 + t713 * t731 - t714 * t772;
t726 = Ifges(6,5) * t759 + Ifges(6,6) * t758 + Ifges(6,3) * t774;
t728 = Ifges(6,1) * t759 + Ifges(6,4) * t758 + Ifges(6,5) * t774;
t668 = -mrSges(6,1) * t698 + mrSges(6,3) * t695 + Ifges(6,4) * t723 + Ifges(6,2) * t722 + Ifges(6,6) * t743 - pkin(5) * t833 + pkin(12) * t846 + t820 * t681 + t815 * t682 - t759 * t726 + t774 * t728;
t727 = Ifges(6,4) * t759 + Ifges(6,2) * t758 + Ifges(6,6) * t774;
t669 = mrSges(6,2) * t698 - mrSges(6,3) * t694 + Ifges(6,1) * t723 + Ifges(6,4) * t722 + Ifges(6,5) * t743 - pkin(12) * t680 - t681 * t815 + t682 * t820 + t726 * t758 - t727 * t774;
t739 = Ifges(5,5) * t778 + Ifges(5,6) * t777 + Ifges(5,3) * t783;
t740 = Ifges(5,4) * t778 + Ifges(5,2) * t777 + Ifges(5,6) * t783;
t654 = mrSges(5,2) * t710 - mrSges(5,3) * t703 + Ifges(5,1) * t746 + Ifges(5,4) * t745 + Ifges(5,5) * t767 - pkin(11) * t676 - t668 * t816 + t669 * t821 + t739 * t777 - t740 * t783;
t741 = Ifges(5,1) * t778 + Ifges(5,4) * t777 + Ifges(5,5) * t783;
t660 = Ifges(5,4) * t746 + Ifges(5,2) * t745 + Ifges(5,6) * t767 - t778 * t739 + t783 * t741 - mrSges(5,1) * t710 + mrSges(5,3) * t704 - Ifges(6,5) * t723 - Ifges(6,6) * t722 - Ifges(6,3) * t743 - t759 * t727 + t758 * t728 - mrSges(6,1) * t694 + mrSges(6,2) * t695 - Ifges(7,5) * t707 - Ifges(7,6) * t706 - Ifges(7,3) * t738 - t732 * t714 + t731 * t715 - mrSges(7,1) * t690 + mrSges(7,2) * t691 - pkin(5) * t680 - pkin(4) * t676;
t763 = Ifges(4,5) * t785 + Ifges(4,6) * t784 + Ifges(4,3) * t797;
t764 = Ifges(4,4) * t785 + Ifges(4,2) * t784 + Ifges(4,6) * t797;
t641 = mrSges(4,2) * t730 - mrSges(4,3) * t719 + Ifges(4,1) * t771 + Ifges(4,4) * t770 + Ifges(4,5) * t794 - pkin(10) * t667 + t654 * t822 - t660 * t817 + t763 * t784 - t764 * t797;
t765 = Ifges(4,1) * t785 + Ifges(4,4) * t784 + Ifges(4,5) * t797;
t646 = Ifges(4,4) * t771 + Ifges(4,2) * t770 + Ifges(4,6) * t794 - t785 * t763 + t797 * t765 - mrSges(4,1) * t730 + mrSges(4,3) * t720 - Ifges(5,5) * t746 - Ifges(5,6) * t745 - Ifges(5,3) * t767 - t778 * t740 + t777 * t741 - mrSges(5,1) * t703 + mrSges(5,2) * t704 - t816 * t669 - t821 * t668 - pkin(4) * t826 - pkin(11) * t847 - pkin(3) * t667;
t835 = pkin(9) * t659 + t641 * t818 + t646 * t823;
t640 = mrSges(4,1) * t719 - mrSges(4,2) * t720 + Ifges(4,5) * t771 + Ifges(4,6) * t770 + Ifges(4,3) * t794 + pkin(3) * t827 + pkin(10) * t848 + t817 * t654 + t822 * t660 + t785 * t764 - t784 * t765;
t789 = (t811 * t844 + t868) * qJD(1);
t832 = Ifges(3,5) * t814 + (Ifges(3,1) * t809 + Ifges(3,4) * t812) * t811;
t791 = t832 * qJD(1);
t831 = Ifges(3,6) * t814 + (Ifges(3,4) * t809 + Ifges(3,2) * t812) * t811;
t634 = -mrSges(3,1) * t786 + mrSges(3,3) * t776 - pkin(2) * t652 - t810 * t640 + (-t789 * t865 + t791 * t814) * qJD(1) + t835 * t813 + t831 * qJDD(1);
t790 = t831 * qJD(1);
t635 = mrSges(3,2) * t786 - mrSges(3,3) * t775 + t823 * t641 - t818 * t646 + (t789 * t861 - t790 * t814) * qJD(1) + (-t652 * t810 - t653 * t813) * pkin(9) + t832 * qJDD(1);
t834 = qJ(2) * t645 + t634 * t812 + t635 * t809;
t633 = qJDD(1) * t868 + mrSges(3,1) * t775 - mrSges(3,2) * t776 + pkin(2) * t653 + t813 * t640 + t835 * t810 + (t844 * qJDD(1) + (t790 * t809 - t791 * t812) * qJD(1)) * t811;
t632 = -mrSges(2,2) * g(3) - mrSges(2,3) * t806 + Ifges(2,5) * qJDD(1) - t825 * Ifges(2,6) - t809 * t634 + t812 * t635 + (-t638 * t811 - t639 * t814) * qJ(2);
t631 = mrSges(2,1) * g(3) + mrSges(2,3) * t807 + t825 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t638 - t811 * t633 + t814 * t834;
t1 = [-m(1) * g(1) + t849; -m(1) * g(2) + t856; (-m(1) - m(2)) * g(3) + t638; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t856 - t819 * t631 + t824 * t632; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t849 + t824 * t631 + t819 * t632; -mrSges(1,1) * g(2) + mrSges(2,1) * t806 + mrSges(1,2) * g(1) - mrSges(2,2) * t807 + Ifges(2,3) * qJDD(1) + pkin(1) * t639 + t814 * t633 + t811 * t834;];
tauB  = t1;
