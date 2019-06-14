% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 20:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR11_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR11_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:14:51
% EndTime: 2019-05-05 20:16:35
% DurationCPUTime: 107.58s
% Computational Cost: add. (1612205->402), mult. (5116504->549), div. (0->0), fcn. (4394391->16), ass. (0->182)
t823 = cos(pkin(6));
t826 = sin(qJ(3));
t822 = cos(pkin(7));
t873 = cos(qJ(3));
t856 = t822 * t873;
t818 = sin(pkin(7));
t857 = t818 * t873;
t819 = sin(pkin(6));
t821 = cos(pkin(12));
t865 = t819 * t821;
t817 = sin(pkin(12));
t868 = t817 * t819;
t874 = -t823 * t857 + t826 * t868 - t856 * t865;
t872 = pkin(9) * t817;
t871 = Ifges(3,3) * t823;
t870 = pkin(9) * qJDD(1);
t869 = qJ(2) * t819;
t867 = t817 * t823;
t866 = t818 * t826;
t864 = t819 * t822;
t863 = t821 * t823;
t862 = t822 * t826;
t799 = (t818 * t823 + t821 * t864) * qJD(1) * pkin(9);
t827 = sin(qJ(1));
t830 = cos(qJ(1));
t814 = -g(1) * t830 - g(2) * t827;
t831 = qJD(1) ^ 2;
t803 = -pkin(1) * t831 + qJDD(1) * t869 + t814;
t844 = -pkin(2) * t821 - t818 * t872;
t860 = qJD(1) * t819;
t840 = qJD(1) * t844 * t860 + t822 * t870;
t813 = t827 * g(1) - g(2) * t830;
t802 = qJDD(1) * pkin(1) + t831 * t869 + t813;
t855 = qJD(2) * t860;
t845 = -g(3) * t865 + t802 * t863 - 0.2e1 * t817 * t855;
t752 = (pkin(2) * qJDD(1) + qJD(1) * t799) * t823 + (-t819 * t840 - t803) * t817 + t845;
t804 = (pkin(2) * t823 - t864 * t872) * qJD(1);
t858 = t802 * t867 + (t803 + 0.2e1 * t855) * t821;
t753 = (-qJD(1) * t804 + t818 * t870) * t823 + (-g(3) * t817 + t821 * t840) * t819 + t858;
t854 = -t823 * g(3) + qJDD(2);
t761 = (-t802 + t844 * qJDD(1) + (-t799 * t821 + t804 * t817) * qJD(1)) * t819 + t854;
t725 = t752 * t862 + t873 * t753 + t761 * t866;
t787 = t874 * qJD(1);
t833 = t823 * t866 + (t817 * t873 + t821 * t862) * t819;
t788 = t833 * qJD(1);
t773 = mrSges(4,1) * t787 + mrSges(4,2) * t788;
t774 = qJD(3) * t788 + qJDD(1) * t874;
t841 = -t818 * t865 + t822 * t823;
t800 = qJD(1) * t841 + qJD(3);
t784 = mrSges(4,1) * t800 - mrSges(4,3) * t788;
t797 = qJDD(1) * t841 + qJDD(3);
t772 = pkin(3) * t787 - qJ(4) * t788;
t796 = t800 ^ 2;
t716 = -pkin(3) * t796 + qJ(4) * t797 - t772 * t787 + t725;
t739 = -t818 * t752 + t822 * t761;
t775 = -t787 * qJD(3) + qJDD(1) * t833;
t719 = (t787 * t800 - t775) * qJ(4) + (t788 * t800 + t774) * pkin(3) + t739;
t816 = sin(pkin(13));
t820 = cos(pkin(13));
t782 = t788 * t820 + t800 * t816;
t708 = -0.2e1 * qJD(4) * t782 - t816 * t716 + t820 * t719;
t767 = t775 * t820 + t797 * t816;
t781 = -t788 * t816 + t800 * t820;
t705 = (t781 * t787 - t767) * pkin(10) + (t781 * t782 + t774) * pkin(4) + t708;
t709 = 0.2e1 * qJD(4) * t781 + t820 * t716 + t816 * t719;
t765 = pkin(4) * t787 - pkin(10) * t782;
t766 = -t775 * t816 + t797 * t820;
t778 = t781 ^ 2;
t707 = -pkin(4) * t778 + pkin(10) * t766 - t765 * t787 + t709;
t825 = sin(qJ(5));
t829 = cos(qJ(5));
t702 = t825 * t705 + t829 * t707;
t756 = t781 * t825 + t782 * t829;
t729 = -qJD(5) * t756 + t766 * t829 - t767 * t825;
t755 = t781 * t829 - t782 * t825;
t737 = -mrSges(6,1) * t755 + mrSges(6,2) * t756;
t786 = qJD(5) + t787;
t743 = mrSges(6,1) * t786 - mrSges(6,3) * t756;
t771 = qJDD(5) + t774;
t738 = -pkin(5) * t755 - pkin(11) * t756;
t785 = t786 ^ 2;
t700 = -pkin(5) * t785 + pkin(11) * t771 + t738 * t755 + t702;
t724 = t752 * t856 - t826 * t753 + t761 * t857;
t715 = -t797 * pkin(3) - t796 * qJ(4) + t788 * t772 + qJDD(4) - t724;
t710 = -t766 * pkin(4) - t778 * pkin(10) + t782 * t765 + t715;
t730 = qJD(5) * t755 + t766 * t825 + t767 * t829;
t703 = (-t755 * t786 - t730) * pkin(11) + (t756 * t786 - t729) * pkin(5) + t710;
t824 = sin(qJ(6));
t828 = cos(qJ(6));
t697 = -t700 * t824 + t703 * t828;
t740 = -t756 * t824 + t786 * t828;
t713 = qJD(6) * t740 + t730 * t828 + t771 * t824;
t741 = t756 * t828 + t786 * t824;
t726 = -mrSges(7,1) * t740 + mrSges(7,2) * t741;
t728 = qJDD(6) - t729;
t754 = qJD(6) - t755;
t731 = -mrSges(7,2) * t754 + mrSges(7,3) * t740;
t695 = m(7) * t697 + mrSges(7,1) * t728 - mrSges(7,3) * t713 - t726 * t741 + t731 * t754;
t698 = t700 * t828 + t703 * t824;
t712 = -qJD(6) * t741 - t730 * t824 + t771 * t828;
t732 = mrSges(7,1) * t754 - mrSges(7,3) * t741;
t696 = m(7) * t698 - mrSges(7,2) * t728 + mrSges(7,3) * t712 + t726 * t740 - t732 * t754;
t850 = -t695 * t824 + t828 * t696;
t683 = m(6) * t702 - mrSges(6,2) * t771 + mrSges(6,3) * t729 + t737 * t755 - t743 * t786 + t850;
t701 = t705 * t829 - t707 * t825;
t742 = -mrSges(6,2) * t786 + mrSges(6,3) * t755;
t699 = -pkin(5) * t771 - pkin(11) * t785 + t738 * t756 - t701;
t837 = -m(7) * t699 + t712 * mrSges(7,1) - mrSges(7,2) * t713 + t740 * t731 - t732 * t741;
t691 = m(6) * t701 + mrSges(6,1) * t771 - mrSges(6,3) * t730 - t737 * t756 + t742 * t786 + t837;
t680 = t825 * t683 + t829 * t691;
t757 = -mrSges(5,1) * t781 + mrSges(5,2) * t782;
t763 = -mrSges(5,2) * t787 + mrSges(5,3) * t781;
t678 = m(5) * t708 + mrSges(5,1) * t774 - mrSges(5,3) * t767 - t757 * t782 + t763 * t787 + t680;
t764 = mrSges(5,1) * t787 - mrSges(5,3) * t782;
t851 = t829 * t683 - t691 * t825;
t679 = m(5) * t709 - mrSges(5,2) * t774 + mrSges(5,3) * t766 + t757 * t781 - t764 * t787 + t851;
t852 = -t678 * t816 + t820 * t679;
t669 = m(4) * t725 - mrSges(4,2) * t797 - mrSges(4,3) * t774 - t773 * t787 - t784 * t800 + t852;
t672 = t820 * t678 + t816 * t679;
t783 = -mrSges(4,2) * t800 - mrSges(4,3) * t787;
t671 = m(4) * t739 + mrSges(4,1) * t774 + mrSges(4,2) * t775 + t783 * t787 + t784 * t788 + t672;
t687 = t828 * t695 + t824 * t696;
t834 = m(6) * t710 - t729 * mrSges(6,1) + mrSges(6,2) * t730 - t755 * t742 + t743 * t756 + t687;
t832 = -m(5) * t715 + t766 * mrSges(5,1) - mrSges(5,2) * t767 + t781 * t763 - t764 * t782 - t834;
t686 = m(4) * t724 + mrSges(4,1) * t797 - mrSges(4,3) * t775 - t773 * t788 + t783 * t800 + t832;
t658 = t669 * t862 - t818 * t671 + t686 * t856;
t779 = -t817 * t803 + t845;
t848 = -mrSges(3,1) * t821 + mrSges(3,2) * t817;
t801 = t848 * t860;
t842 = -mrSges(3,2) * t823 + mrSges(3,3) * t865;
t806 = t842 * qJD(1);
t843 = mrSges(3,1) * t823 - mrSges(3,3) * t868;
t654 = m(3) * t779 + t843 * qJDD(1) + (-t801 * t868 + t806 * t823) * qJD(1) + t658;
t657 = t669 * t866 + t822 * t671 + t686 * t857;
t789 = -t819 * t802 + t854;
t805 = t843 * qJD(1);
t656 = m(3) * t789 + (t848 * qJDD(1) + (t805 * t817 - t806 * t821) * qJD(1)) * t819 + t657;
t665 = t873 * t669 - t826 * t686;
t780 = -g(3) * t868 + t858;
t664 = m(3) * t780 + t842 * qJDD(1) + (t801 * t865 - t805 * t823) * qJD(1) + t665;
t644 = t654 * t863 - t656 * t819 + t664 * t867;
t642 = m(2) * t813 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t831 + t644;
t650 = -t654 * t817 + t821 * t664;
t649 = m(2) * t814 - mrSges(2,1) * t831 - qJDD(1) * mrSges(2,2) + t650;
t861 = t830 * t642 + t827 * t649;
t643 = t654 * t865 + t823 * t656 + t664 * t868;
t853 = -t642 * t827 + t830 * t649;
t847 = Ifges(3,5) * t817 + Ifges(3,6) * t821;
t720 = Ifges(7,5) * t741 + Ifges(7,6) * t740 + Ifges(7,3) * t754;
t722 = Ifges(7,1) * t741 + Ifges(7,4) * t740 + Ifges(7,5) * t754;
t688 = -mrSges(7,1) * t699 + mrSges(7,3) * t698 + Ifges(7,4) * t713 + Ifges(7,2) * t712 + Ifges(7,6) * t728 - t720 * t741 + t722 * t754;
t721 = Ifges(7,4) * t741 + Ifges(7,2) * t740 + Ifges(7,6) * t754;
t689 = mrSges(7,2) * t699 - mrSges(7,3) * t697 + Ifges(7,1) * t713 + Ifges(7,4) * t712 + Ifges(7,5) * t728 + t720 * t740 - t721 * t754;
t733 = Ifges(6,5) * t756 + Ifges(6,6) * t755 + Ifges(6,3) * t786;
t734 = Ifges(6,4) * t756 + Ifges(6,2) * t755 + Ifges(6,6) * t786;
t673 = mrSges(6,2) * t710 - mrSges(6,3) * t701 + Ifges(6,1) * t730 + Ifges(6,4) * t729 + Ifges(6,5) * t771 - pkin(11) * t687 - t688 * t824 + t689 * t828 + t733 * t755 - t734 * t786;
t735 = Ifges(6,1) * t756 + Ifges(6,4) * t755 + Ifges(6,5) * t786;
t674 = -mrSges(6,1) * t710 - mrSges(7,1) * t697 + mrSges(7,2) * t698 + mrSges(6,3) * t702 + Ifges(6,4) * t730 - Ifges(7,5) * t713 + Ifges(6,2) * t729 + Ifges(6,6) * t771 - Ifges(7,6) * t712 - Ifges(7,3) * t728 - pkin(5) * t687 - t721 * t741 + t722 * t740 - t733 * t756 + t735 * t786;
t744 = Ifges(5,5) * t782 + Ifges(5,6) * t781 + Ifges(5,3) * t787;
t746 = Ifges(5,1) * t782 + Ifges(5,4) * t781 + Ifges(5,5) * t787;
t659 = -mrSges(5,1) * t715 + mrSges(5,3) * t709 + Ifges(5,4) * t767 + Ifges(5,2) * t766 + Ifges(5,6) * t774 - pkin(4) * t834 + pkin(10) * t851 + t825 * t673 + t829 * t674 - t782 * t744 + t787 * t746;
t745 = Ifges(5,4) * t782 + Ifges(5,2) * t781 + Ifges(5,6) * t787;
t660 = mrSges(5,2) * t715 - mrSges(5,3) * t708 + Ifges(5,1) * t767 + Ifges(5,4) * t766 + Ifges(5,5) * t774 - pkin(10) * t680 + t673 * t829 - t674 * t825 + t744 * t781 - t745 * t787;
t769 = Ifges(4,4) * t788 - Ifges(4,2) * t787 + Ifges(4,6) * t800;
t770 = Ifges(4,1) * t788 - Ifges(4,4) * t787 + Ifges(4,5) * t800;
t645 = mrSges(4,1) * t724 - mrSges(4,2) * t725 + Ifges(4,5) * t775 - Ifges(4,6) * t774 + Ifges(4,3) * t797 + pkin(3) * t832 + qJ(4) * t852 + t820 * t659 + t816 * t660 + t788 * t769 + t787 * t770;
t792 = (t819 * t847 + t871) * qJD(1);
t836 = Ifges(3,5) * t823 + (Ifges(3,1) * t817 + Ifges(3,4) * t821) * t819;
t794 = t836 * qJD(1);
t835 = Ifges(3,6) * t823 + (Ifges(3,4) * t817 + Ifges(3,2) * t821) * t819;
t768 = Ifges(4,5) * t788 - Ifges(4,6) * t787 + Ifges(4,3) * t800;
t646 = mrSges(4,2) * t739 - mrSges(4,3) * t724 + Ifges(4,1) * t775 - Ifges(4,4) * t774 + Ifges(4,5) * t797 - qJ(4) * t672 - t659 * t816 + t660 * t820 - t768 * t787 - t769 * t800;
t651 = -pkin(11) * t850 - mrSges(6,1) * t701 + mrSges(6,2) * t702 - mrSges(5,1) * t708 + mrSges(5,2) * t709 - pkin(4) * t680 + (-Ifges(5,3) - Ifges(4,2)) * t774 - pkin(5) * t837 - pkin(3) * t672 + mrSges(4,3) * t725 - Ifges(6,6) * t729 - Ifges(6,5) * t730 - mrSges(4,1) * t739 + t755 * t735 - t756 * t734 - Ifges(5,6) * t766 - Ifges(5,5) * t767 - Ifges(6,3) * t771 + Ifges(4,4) * t775 + t781 * t746 - t782 * t745 - t788 * t768 + Ifges(4,6) * t797 + t800 * t770 - t824 * t689 - t828 * t688;
t838 = pkin(9) * t665 + t646 * t826 + t651 * t873;
t639 = -mrSges(3,1) * t789 + mrSges(3,3) * t780 - pkin(2) * t657 - t818 * t645 + (-t792 * t868 + t794 * t823) * qJD(1) + t838 * t822 + t835 * qJDD(1);
t793 = t835 * qJD(1);
t640 = mrSges(3,2) * t789 - mrSges(3,3) * t779 + t873 * t646 - t826 * t651 + (t792 * t865 - t793 * t823) * qJD(1) + (-t657 * t818 - t658 * t822) * pkin(9) + t836 * qJDD(1);
t839 = qJ(2) * t650 + t639 * t821 + t640 * t817;
t638 = qJDD(1) * t871 + mrSges(3,1) * t779 - mrSges(3,2) * t780 + pkin(2) * t658 + t822 * t645 + t838 * t818 + (t847 * qJDD(1) + (t793 * t817 - t794 * t821) * qJD(1)) * t819;
t637 = -mrSges(2,2) * g(3) - mrSges(2,3) * t813 + Ifges(2,5) * qJDD(1) - t831 * Ifges(2,6) - t817 * t639 + t821 * t640 + (-t643 * t819 - t644 * t823) * qJ(2);
t636 = mrSges(2,1) * g(3) + mrSges(2,3) * t814 + t831 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t643 - t819 * t638 + t823 * t839;
t1 = [-m(1) * g(1) + t853; -m(1) * g(2) + t861; (-m(1) - m(2)) * g(3) + t643; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t861 - t827 * t636 + t830 * t637; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t853 + t830 * t636 + t827 * t637; -mrSges(1,1) * g(2) + mrSges(2,1) * t813 + mrSges(1,2) * g(1) - mrSges(2,2) * t814 + Ifges(2,3) * qJDD(1) + pkin(1) * t644 + t823 * t638 + t819 * t839;];
tauB  = t1;
