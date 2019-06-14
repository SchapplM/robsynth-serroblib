% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-05-06 10:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:38:55
% EndTime: 2019-05-06 10:39:13
% DurationCPUTime: 8.46s
% Computational Cost: add. (84240->366), mult. (194217->447), div. (0->0), fcn. (131356->10), ass. (0->150)
t894 = Ifges(3,1) + Ifges(4,1) + Ifges(5,1);
t872 = Ifges(3,4) - Ifges(4,5) - Ifges(5,4);
t893 = Ifges(5,5) - Ifges(3,5) - Ifges(4,4);
t870 = Ifges(3,6) + Ifges(5,6) - Ifges(4,6);
t892 = Ifges(4,3) + Ifges(5,2) + Ifges(3,2);
t891 = Ifges(5,3) + Ifges(3,3) + Ifges(4,2);
t834 = cos(pkin(6));
t826 = t834 * qJD(1) + qJD(2);
t890 = t826 ^ 2;
t889 = -2 * qJD(4);
t888 = pkin(2) * t826;
t887 = -mrSges(5,2) - mrSges(4,3);
t886 = mrSges(3,3) + mrSges(4,2);
t833 = sin(pkin(6));
t841 = cos(qJ(2));
t874 = qJD(1) * t841;
t863 = t833 * t874;
t794 = t826 * mrSges(5,2) - mrSges(5,3) * t863;
t885 = t826 * t794;
t843 = qJD(1) ^ 2;
t884 = t833 ^ 2 * t843;
t837 = sin(qJ(2));
t883 = t833 * t837;
t882 = t833 * t841;
t881 = t834 * t837;
t880 = t834 * t841;
t869 = g(3) * t883;
t838 = sin(qJ(1));
t842 = cos(qJ(1));
t817 = t838 * g(1) - t842 * g(2);
t797 = t843 * t833 * pkin(8) + qJDD(1) * pkin(1) + t817;
t818 = -t842 * g(1) - t838 * g(2);
t873 = qJDD(1) * t833;
t798 = -t843 * pkin(1) + pkin(8) * t873 + t818;
t878 = t797 * t881 + t841 * t798;
t747 = -t869 + t878;
t875 = qJD(1) * t833;
t801 = (mrSges(5,1) * t841 + mrSges(5,2) * t837) * t875;
t802 = (-mrSges(3,1) * t841 + mrSges(3,2) * t837) * t875;
t864 = t837 * t875;
t805 = -qJD(2) * t864 + t841 * t873;
t825 = t834 * qJDD(1) + qJDD(2);
t799 = (-pkin(2) * t841 - qJ(3) * t837) * t875;
t854 = -pkin(2) * t890 + t825 * qJ(3) + 0.2e1 * qJD(3) * t826 + t799 * t863 + t878;
t726 = t854 - t869;
t793 = -t826 * mrSges(4,1) + mrSges(4,2) * t864;
t800 = (-mrSges(4,1) * t841 - mrSges(4,3) * t837) * t875;
t790 = -t826 * pkin(3) - qJ(4) * t864;
t868 = t841 ^ 2 * t884;
t845 = -pkin(3) * t868 - t805 * qJ(4) + t826 * t790 + t863 * t889 + t854;
t718 = t845 - t869;
t804 = (qJD(2) * t874 + qJDD(1) * t837) * t833;
t763 = -t834 * g(3) - t833 * t797;
t858 = t826 * t863;
t857 = -t805 * pkin(2) + t763 + (-t804 - t858) * qJ(3);
t859 = 0.2e1 * t864;
t846 = -qJ(4) * t868 + qJD(3) * t859 + t790 * t864 + qJDD(4) - t857;
t712 = -t804 * pkin(9) + (pkin(3) + pkin(4)) * t805 + (-pkin(9) * t841 + (-pkin(2) - pkin(4)) * t837) * t826 * t875 + t846;
t803 = (pkin(4) * t841 - pkin(9) * t837) * t875;
t714 = -t890 * pkin(4) - t825 * pkin(9) + (-g(3) * t837 - t803 * t874) * t833 + t845;
t836 = sin(qJ(5));
t840 = cos(qJ(5));
t709 = t836 * t712 + t840 * t714;
t778 = -t836 * t826 + t840 * t864;
t744 = -t778 * qJD(5) - t836 * t804 - t840 * t825;
t777 = -t840 * t826 - t836 * t864;
t748 = -t777 * mrSges(6,1) + t778 * mrSges(6,2);
t810 = qJD(5) + t863;
t753 = t810 * mrSges(6,1) - t778 * mrSges(6,3);
t788 = qJDD(5) + t805;
t749 = -t777 * pkin(5) - t778 * pkin(10);
t808 = t810 ^ 2;
t707 = -t808 * pkin(5) + t788 * pkin(10) + t777 * t749 + t709;
t746 = -g(3) * t882 + t797 * t880 - t837 * t798;
t732 = -t825 * pkin(2) - qJ(3) * t890 + t799 * t864 + qJDD(3) - t746;
t849 = -t804 * qJ(4) + t732 + (-t837 * t841 * t884 - t825) * pkin(3);
t716 = t825 * pkin(4) - pkin(9) * t890 - qJ(4) * t858 + qJD(4) * t859 + t803 * t864 - t849;
t745 = t777 * qJD(5) + t840 * t804 - t836 * t825;
t710 = (-t777 * t810 - t745) * pkin(10) + (t778 * t810 - t744) * pkin(5) + t716;
t835 = sin(qJ(6));
t839 = cos(qJ(6));
t704 = -t835 * t707 + t839 * t710;
t750 = -t835 * t778 + t839 * t810;
t723 = t750 * qJD(6) + t839 * t745 + t835 * t788;
t751 = t839 * t778 + t835 * t810;
t733 = -t750 * mrSges(7,1) + t751 * mrSges(7,2);
t775 = qJD(6) - t777;
t734 = -t775 * mrSges(7,2) + t750 * mrSges(7,3);
t741 = qJDD(6) - t744;
t702 = m(7) * t704 + t741 * mrSges(7,1) - t723 * mrSges(7,3) - t751 * t733 + t775 * t734;
t705 = t839 * t707 + t835 * t710;
t722 = -t751 * qJD(6) - t835 * t745 + t839 * t788;
t735 = t775 * mrSges(7,1) - t751 * mrSges(7,3);
t703 = m(7) * t705 - t741 * mrSges(7,2) + t722 * mrSges(7,3) + t750 * t733 - t775 * t735;
t860 = -t835 * t702 + t839 * t703;
t694 = m(6) * t709 - t788 * mrSges(6,2) + t744 * mrSges(6,3) + t777 * t748 - t810 * t753 + t860;
t708 = t840 * t712 - t836 * t714;
t752 = -t810 * mrSges(6,2) + t777 * mrSges(6,3);
t706 = -t788 * pkin(5) - t808 * pkin(10) + t778 * t749 - t708;
t851 = -m(7) * t706 + t722 * mrSges(7,1) - t723 * mrSges(7,2) + t750 * t734 - t751 * t735;
t698 = m(6) * t708 + t788 * mrSges(6,1) - t745 * mrSges(6,3) - t778 * t748 + t810 * t752 + t851;
t861 = t840 * t694 - t836 * t698;
t856 = m(5) * t718 - t805 * mrSges(5,3) + t861;
t850 = m(4) * t726 + t825 * mrSges(4,3) + t826 * t793 + t800 * t863 + t856;
t791 = -t826 * mrSges(5,1) - mrSges(5,3) * t864;
t877 = -t826 * mrSges(3,1) + mrSges(3,3) * t864 + t791;
t684 = m(3) * t747 + t877 * t826 + (-mrSges(3,2) + mrSges(5,2)) * t825 + t886 * t805 + (-t801 + t802) * t863 + t850;
t796 = mrSges(4,2) * t863 + t826 * mrSges(4,3);
t727 = (-0.2e1 * qJD(3) + t888) * t864 + t857;
t688 = t836 * t694 + t840 * t698;
t720 = t805 * pkin(3) - t864 * t888 + t846;
t852 = -m(5) * t720 - t805 * mrSges(5,1) - t688;
t848 = m(4) * t727 - t805 * mrSges(4,1) + t852;
t876 = -t826 * mrSges(3,2) + mrSges(3,3) * t863 + t794;
t686 = m(3) * t763 - t805 * mrSges(3,1) + (mrSges(3,2) + t887) * t804 + ((-t796 - t876) * t841 + (-t793 - t877) * t837) * t875 + t848;
t719 = (qJ(4) * t826 * t841 + t837 * t889) * t875 + t849;
t695 = t839 * t702 + t835 * t703;
t855 = -m(6) * t716 + t744 * mrSges(6,1) - t745 * mrSges(6,2) + t777 * t752 - t778 * t753 - t695;
t847 = -m(5) * t719 + t804 * mrSges(5,3) + t801 * t864 - t855;
t844 = -m(4) * t732 + t825 * mrSges(4,1) + t826 * t796 + t847;
t691 = m(3) * t746 + (-t800 - t802) * t864 + t876 * t826 + t844 + (mrSges(3,1) + mrSges(5,1)) * t825 - t886 * t804;
t674 = t684 * t881 - t833 * t686 + t691 * t880;
t672 = m(2) * t817 + qJDD(1) * mrSges(2,1) - t843 * mrSges(2,2) + t674;
t678 = t841 * t684 - t837 * t691;
t677 = m(2) * t818 - t843 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t678;
t879 = t842 * t672 + t838 * t677;
t673 = t684 * t883 + t834 * t686 + t691 * t882;
t867 = (-t837 * t893 + t841 * t870) * t875 + t891 * t826;
t866 = (-t837 * t872 - t841 * t892) * t875 - t870 * t826;
t865 = (-t837 * t894 - t872 * t841) * t875 + t893 * t826;
t862 = -t838 * t672 + t842 * t677;
t687 = t887 * t804 + ((-t794 - t796) * t841 + (-t791 - t793) * t837) * t875 + t848;
t728 = Ifges(7,5) * t751 + Ifges(7,6) * t750 + Ifges(7,3) * t775;
t730 = Ifges(7,1) * t751 + Ifges(7,4) * t750 + Ifges(7,5) * t775;
t696 = -mrSges(7,1) * t706 + mrSges(7,3) * t705 + Ifges(7,4) * t723 + Ifges(7,2) * t722 + Ifges(7,6) * t741 - t751 * t728 + t775 * t730;
t729 = Ifges(7,4) * t751 + Ifges(7,2) * t750 + Ifges(7,6) * t775;
t697 = mrSges(7,2) * t706 - mrSges(7,3) * t704 + Ifges(7,1) * t723 + Ifges(7,4) * t722 + Ifges(7,5) * t741 + t750 * t728 - t775 * t729;
t737 = Ifges(6,4) * t778 + Ifges(6,2) * t777 + Ifges(6,6) * t810;
t738 = Ifges(6,1) * t778 + Ifges(6,4) * t777 + Ifges(6,5) * t810;
t669 = ((pkin(3) * t794 + qJ(4) * t801) * t841 + (pkin(3) * t791 - t867) * t837) * t875 + t839 * t696 + t835 * t697 + Ifges(6,3) * t788 - t777 * t738 + t778 * t737 - mrSges(3,1) * t763 + Ifges(6,6) * t744 + Ifges(6,5) * t745 + mrSges(3,3) * t747 - pkin(3) * t852 + pkin(5) * t851 - qJ(4) * t856 + t892 * t805 + (-qJ(4) * t791 - t865) * t826 + (-qJ(4) * mrSges(5,2) + t870) * t825 + (pkin(3) * mrSges(5,2) + t872) * t804 + pkin(10) * t860 + mrSges(5,1) * t720 + mrSges(4,2) * t726 - mrSges(4,1) * t727 - mrSges(5,3) * t718 - pkin(2) * t687 + pkin(4) * t688 + mrSges(6,1) * t708 - mrSges(6,2) * t709;
t736 = Ifges(6,5) * t778 + Ifges(6,6) * t777 + Ifges(6,3) * t810;
t679 = mrSges(6,2) * t716 - mrSges(6,3) * t708 + Ifges(6,1) * t745 + Ifges(6,4) * t744 + Ifges(6,5) * t788 - pkin(10) * t695 - t835 * t696 + t839 * t697 + t777 * t736 - t810 * t737;
t680 = -mrSges(6,1) * t716 - mrSges(7,1) * t704 + mrSges(7,2) * t705 + mrSges(6,3) * t709 + Ifges(6,4) * t745 - Ifges(7,5) * t723 + Ifges(6,2) * t744 + Ifges(6,6) * t788 - Ifges(7,6) * t722 - Ifges(7,3) * t741 - pkin(5) * t695 - t751 * t729 + t750 * t730 - t778 * t736 + t810 * t738;
t692 = -t825 * mrSges(5,1) - t847 - t885;
t670 = mrSges(3,2) * t763 + mrSges(4,2) * t732 + mrSges(5,2) * t720 - mrSges(3,3) * t746 - mrSges(4,3) * t727 - mrSges(5,3) * t719 - pkin(9) * t688 - qJ(3) * t687 - qJ(4) * t692 + t840 * t679 - t836 * t680 + t866 * t826 - t893 * t825 + t872 * t805 + t894 * t804 + t867 * t863;
t853 = pkin(8) * t678 + t669 * t841 + t670 * t837;
t668 = -t836 * t679 - pkin(9) * t861 + pkin(2) * (t844 + t885) - t840 * t680 + qJ(3) * (t826 * t791 + t850) + mrSges(3,1) * t746 - mrSges(3,2) * t747 - mrSges(4,1) * t732 - pkin(4) * t855 + mrSges(4,3) * t726 + mrSges(5,2) * t718 - mrSges(5,1) * t719 - pkin(3) * t692 + (pkin(2) * mrSges(5,1) + qJ(3) * mrSges(5,2) + t891) * t825 + (qJ(3) * mrSges(4,2) + t870) * t805 + (-pkin(2) * mrSges(4,2) - t893) * t804 + ((-qJ(3) * t801 + t865) * t841 + (-pkin(2) * t800 - t866) * t837) * t875;
t667 = -mrSges(2,2) * g(3) - mrSges(2,3) * t817 + Ifges(2,5) * qJDD(1) - t843 * Ifges(2,6) - t837 * t669 + t841 * t670 + (-t673 * t833 - t674 * t834) * pkin(8);
t666 = mrSges(2,1) * g(3) + mrSges(2,3) * t818 + t843 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t673 - t833 * t668 + t853 * t834;
t1 = [-m(1) * g(1) + t862; -m(1) * g(2) + t879; (-m(1) - m(2)) * g(3) + t673; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t879 - t838 * t666 + t842 * t667; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t862 + t842 * t666 + t838 * t667; -mrSges(1,1) * g(2) + mrSges(2,1) * t817 + mrSges(1,2) * g(1) - mrSges(2,2) * t818 + Ifges(2,3) * qJDD(1) + pkin(1) * t674 + t834 * t668 + t853 * t833;];
tauB  = t1;
