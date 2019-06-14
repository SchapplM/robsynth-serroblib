% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 07:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:43:00
% EndTime: 2019-05-05 07:43:27
% DurationCPUTime: 26.46s
% Computational Cost: add. (468440->341), mult. (943339->438), div. (0->0), fcn. (688360->14), ass. (0->145)
t856 = sin(pkin(11));
t859 = cos(pkin(11));
t845 = g(1) * t856 - g(2) * t859;
t846 = -g(1) * t859 - g(2) * t856;
t854 = -g(3) + qJDD(1);
t868 = cos(qJ(2));
t860 = cos(pkin(6));
t864 = sin(qJ(2));
t889 = t860 * t864;
t857 = sin(pkin(6));
t890 = t857 * t864;
t806 = t845 * t889 + t868 * t846 + t854 * t890;
t870 = qJD(2) ^ 2;
t801 = -pkin(2) * t870 + qJDD(2) * pkin(8) + t806;
t825 = -t845 * t857 + t854 * t860;
t863 = sin(qJ(3));
t867 = cos(qJ(3));
t796 = t867 * t801 + t863 * t825;
t842 = (-pkin(3) * t867 - pkin(9) * t863) * qJD(2);
t869 = qJD(3) ^ 2;
t886 = qJD(2) * t867;
t780 = -pkin(3) * t869 + qJDD(3) * pkin(9) + t842 * t886 + t796;
t805 = -t864 * t846 + (t845 * t860 + t854 * t857) * t868;
t800 = -qJDD(2) * pkin(2) - t870 * pkin(8) - t805;
t885 = qJD(2) * qJD(3);
t884 = t867 * t885;
t843 = qJDD(2) * t863 + t884;
t852 = t863 * t885;
t844 = qJDD(2) * t867 - t852;
t785 = (-t843 - t884) * pkin(9) + (-t844 + t852) * pkin(3) + t800;
t862 = sin(qJ(4));
t866 = cos(qJ(4));
t769 = -t862 * t780 + t866 * t785;
t887 = qJD(2) * t863;
t839 = qJD(3) * t866 - t862 * t887;
t816 = qJD(4) * t839 + qJDD(3) * t862 + t843 * t866;
t836 = qJDD(4) - t844;
t840 = qJD(3) * t862 + t866 * t887;
t851 = qJD(4) - t886;
t759 = (t839 * t851 - t816) * qJ(5) + (t839 * t840 + t836) * pkin(4) + t769;
t770 = t866 * t780 + t862 * t785;
t815 = -qJD(4) * t840 + qJDD(3) * t866 - t843 * t862;
t823 = pkin(4) * t851 - qJ(5) * t840;
t835 = t839 ^ 2;
t761 = -pkin(4) * t835 + qJ(5) * t815 - t823 * t851 + t770;
t855 = sin(pkin(12));
t858 = cos(pkin(12));
t819 = t839 * t855 + t840 * t858;
t753 = -0.2e1 * qJD(5) * t819 + t858 * t759 - t855 * t761;
t791 = t815 * t855 + t816 * t858;
t818 = t839 * t858 - t840 * t855;
t751 = (t818 * t851 - t791) * pkin(10) + (t818 * t819 + t836) * pkin(5) + t753;
t754 = 0.2e1 * qJD(5) * t818 + t855 * t759 + t858 * t761;
t790 = t815 * t858 - t816 * t855;
t804 = pkin(5) * t851 - pkin(10) * t819;
t817 = t818 ^ 2;
t752 = -pkin(5) * t817 + pkin(10) * t790 - t804 * t851 + t754;
t861 = sin(qJ(6));
t865 = cos(qJ(6));
t749 = t751 * t865 - t752 * t861;
t793 = t818 * t865 - t819 * t861;
t767 = qJD(6) * t793 + t790 * t861 + t791 * t865;
t794 = t818 * t861 + t819 * t865;
t777 = -mrSges(7,1) * t793 + mrSges(7,2) * t794;
t850 = qJD(6) + t851;
t783 = -mrSges(7,2) * t850 + mrSges(7,3) * t793;
t832 = qJDD(6) + t836;
t742 = m(7) * t749 + mrSges(7,1) * t832 - mrSges(7,3) * t767 - t777 * t794 + t783 * t850;
t750 = t751 * t861 + t752 * t865;
t766 = -qJD(6) * t794 + t790 * t865 - t791 * t861;
t784 = mrSges(7,1) * t850 - mrSges(7,3) * t794;
t743 = m(7) * t750 - mrSges(7,2) * t832 + mrSges(7,3) * t766 + t777 * t793 - t784 * t850;
t736 = t865 * t742 + t861 * t743;
t797 = -mrSges(6,1) * t818 + mrSges(6,2) * t819;
t802 = -mrSges(6,2) * t851 + mrSges(6,3) * t818;
t734 = m(6) * t753 + mrSges(6,1) * t836 - mrSges(6,3) * t791 - t797 * t819 + t802 * t851 + t736;
t803 = mrSges(6,1) * t851 - mrSges(6,3) * t819;
t880 = -t742 * t861 + t865 * t743;
t735 = m(6) * t754 - mrSges(6,2) * t836 + mrSges(6,3) * t790 + t797 * t818 - t803 * t851 + t880;
t730 = t858 * t734 + t855 * t735;
t788 = Ifges(6,4) * t819 + Ifges(6,2) * t818 + Ifges(6,6) * t851;
t789 = Ifges(6,1) * t819 + Ifges(6,4) * t818 + Ifges(6,5) * t851;
t808 = Ifges(5,4) * t840 + Ifges(5,2) * t839 + Ifges(5,6) * t851;
t809 = Ifges(5,1) * t840 + Ifges(5,4) * t839 + Ifges(5,5) * t851;
t773 = Ifges(7,4) * t794 + Ifges(7,2) * t793 + Ifges(7,6) * t850;
t774 = Ifges(7,1) * t794 + Ifges(7,4) * t793 + Ifges(7,5) * t850;
t874 = -mrSges(7,1) * t749 + mrSges(7,2) * t750 - Ifges(7,5) * t767 - Ifges(7,6) * t766 - Ifges(7,3) * t832 - t794 * t773 + t793 * t774;
t894 = mrSges(5,1) * t769 + mrSges(6,1) * t753 - mrSges(5,2) * t770 - mrSges(6,2) * t754 + Ifges(5,5) * t816 + Ifges(6,5) * t791 + Ifges(5,6) * t815 + Ifges(6,6) * t790 + pkin(4) * t730 + pkin(5) * t736 + t819 * t788 - t818 * t789 + t840 * t808 - t839 * t809 + (Ifges(5,3) + Ifges(6,3)) * t836 - t874;
t795 = -t863 * t801 + t825 * t867;
t779 = -qJDD(3) * pkin(3) - pkin(9) * t869 + t842 * t887 - t795;
t771 = -pkin(4) * t815 - qJ(5) * t835 + t840 * t823 + qJDD(5) + t779;
t756 = -pkin(5) * t790 - pkin(10) * t817 + t804 * t819 + t771;
t772 = Ifges(7,5) * t794 + Ifges(7,6) * t793 + Ifges(7,3) * t850;
t737 = -mrSges(7,1) * t756 + mrSges(7,3) * t750 + Ifges(7,4) * t767 + Ifges(7,2) * t766 + Ifges(7,6) * t832 - t772 * t794 + t774 * t850;
t738 = mrSges(7,2) * t756 - mrSges(7,3) * t749 + Ifges(7,1) * t767 + Ifges(7,4) * t766 + Ifges(7,5) * t832 + t772 * t793 - t773 * t850;
t787 = Ifges(6,5) * t819 + Ifges(6,6) * t818 + Ifges(6,3) * t851;
t878 = m(7) * t756 - t766 * mrSges(7,1) + t767 * mrSges(7,2) - t793 * t783 + t794 * t784;
t725 = -mrSges(6,1) * t771 + mrSges(6,3) * t754 + Ifges(6,4) * t791 + Ifges(6,2) * t790 + Ifges(6,6) * t836 - pkin(5) * t878 + pkin(10) * t880 + t865 * t737 + t861 * t738 - t819 * t787 + t851 * t789;
t726 = mrSges(6,2) * t771 - mrSges(6,3) * t753 + Ifges(6,1) * t791 + Ifges(6,4) * t790 + Ifges(6,5) * t836 - pkin(10) * t736 - t737 * t861 + t738 * t865 + t787 * t818 - t788 * t851;
t747 = m(6) * t771 - t790 * mrSges(6,1) + t791 * mrSges(6,2) - t818 * t802 + t819 * t803 + t878;
t807 = Ifges(5,5) * t840 + Ifges(5,6) * t839 + Ifges(5,3) * t851;
t881 = -t734 * t855 + t858 * t735;
t708 = -mrSges(5,1) * t779 + mrSges(5,3) * t770 + Ifges(5,4) * t816 + Ifges(5,2) * t815 + Ifges(5,6) * t836 - pkin(4) * t747 + qJ(5) * t881 + t858 * t725 + t855 * t726 - t840 * t807 + t851 * t809;
t709 = mrSges(5,2) * t779 - mrSges(5,3) * t769 + Ifges(5,1) * t816 + Ifges(5,4) * t815 + Ifges(5,5) * t836 - qJ(5) * t730 - t725 * t855 + t726 * t858 + t807 * t839 - t808 * t851;
t820 = -mrSges(5,1) * t839 + mrSges(5,2) * t840;
t822 = -mrSges(5,2) * t851 + mrSges(5,3) * t839;
t728 = m(5) * t769 + mrSges(5,1) * t836 - mrSges(5,3) * t816 - t820 * t840 + t822 * t851 + t730;
t824 = mrSges(5,1) * t851 - mrSges(5,3) * t840;
t729 = m(5) * t770 - mrSges(5,2) * t836 + mrSges(5,3) * t815 + t820 * t839 - t824 * t851 + t881;
t724 = -t728 * t862 + t866 * t729;
t746 = -m(5) * t779 + t815 * mrSges(5,1) - t816 * mrSges(5,2) + t839 * t822 - t840 * t824 - t747;
t830 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t863 + Ifges(4,2) * t867) * qJD(2);
t831 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t863 + Ifges(4,4) * t867) * qJD(2);
t893 = mrSges(4,1) * t795 - mrSges(4,2) * t796 + Ifges(4,5) * t843 + Ifges(4,6) * t844 + Ifges(4,3) * qJDD(3) + pkin(3) * t746 + pkin(9) * t724 + t866 * t708 + t862 * t709 + (t830 * t863 - t831 * t867) * qJD(2);
t723 = t728 * t866 + t729 * t862;
t847 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t887;
t848 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t886;
t873 = -m(4) * t800 + t844 * mrSges(4,1) - mrSges(4,2) * t843 - t847 * t887 + t848 * t886 - t723;
t719 = m(3) * t805 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t870 + t873;
t891 = t719 * t868;
t841 = (-mrSges(4,1) * t867 + mrSges(4,2) * t863) * qJD(2);
t722 = m(4) * t796 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t844 - qJD(3) * t847 + t841 * t886 + t724;
t745 = m(4) * t795 + qJDD(3) * mrSges(4,1) - t843 * mrSges(4,3) + qJD(3) * t848 - t841 * t887 + t746;
t882 = t867 * t722 - t745 * t863;
t713 = m(3) * t806 - mrSges(3,1) * t870 - qJDD(2) * mrSges(3,2) + t882;
t716 = t863 * t722 + t867 * t745;
t715 = m(3) * t825 + t716;
t702 = t713 * t889 - t715 * t857 + t860 * t891;
t700 = m(2) * t845 + t702;
t706 = t868 * t713 - t719 * t864;
t705 = m(2) * t846 + t706;
t888 = t859 * t700 + t856 * t705;
t701 = t713 * t890 + t860 * t715 + t857 * t891;
t883 = -t700 * t856 + t859 * t705;
t879 = m(2) * t854 + t701;
t829 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t863 + Ifges(4,6) * t867) * qJD(2);
t698 = mrSges(4,2) * t800 - mrSges(4,3) * t795 + Ifges(4,1) * t843 + Ifges(4,4) * t844 + Ifges(4,5) * qJDD(3) - pkin(9) * t723 - qJD(3) * t830 - t708 * t862 + t709 * t866 + t829 * t886;
t707 = -mrSges(4,1) * t800 + mrSges(4,3) * t796 + Ifges(4,4) * t843 + Ifges(4,2) * t844 + Ifges(4,6) * qJDD(3) - pkin(3) * t723 + qJD(3) * t831 - t829 * t887 - t894;
t696 = mrSges(3,2) * t825 - mrSges(3,3) * t805 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t870 - pkin(8) * t716 + t698 * t867 - t707 * t863;
t697 = -mrSges(3,1) * t825 + mrSges(3,3) * t806 + t870 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t716 - t893;
t875 = pkin(7) * t706 + t696 * t864 + t697 * t868;
t695 = mrSges(3,1) * t805 - mrSges(3,2) * t806 + Ifges(3,3) * qJDD(2) + pkin(2) * t873 + pkin(8) * t882 + t863 * t698 + t867 * t707;
t694 = mrSges(2,2) * t854 - mrSges(2,3) * t845 + t868 * t696 - t864 * t697 + (-t701 * t857 - t702 * t860) * pkin(7);
t693 = -mrSges(2,1) * t854 + mrSges(2,3) * t846 - pkin(1) * t701 - t857 * t695 + t860 * t875;
t1 = [-m(1) * g(1) + t883; -m(1) * g(2) + t888; -m(1) * g(3) + t879; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t888 - t856 * t693 + t859 * t694; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t883 + t859 * t693 + t856 * t694; -mrSges(1,1) * g(2) + mrSges(2,1) * t845 + mrSges(1,2) * g(1) - mrSges(2,2) * t846 + pkin(1) * t702 + t860 * t695 + t857 * t875; t879; t695; t893; t894; t747; -t874;];
tauJB  = t1;
