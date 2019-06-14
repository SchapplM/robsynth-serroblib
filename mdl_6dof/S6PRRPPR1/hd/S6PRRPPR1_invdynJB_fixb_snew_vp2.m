% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-05-05 02:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:31:26
% EndTime: 2019-05-05 02:31:47
% DurationCPUTime: 20.98s
% Computational Cost: add. (334743->340), mult. (745017->442), div. (0->0), fcn. (537529->14), ass. (0->143)
t907 = -2 * qJD(4);
t867 = sin(pkin(10));
t870 = cos(pkin(10));
t855 = g(1) * t867 - g(2) * t870;
t856 = -g(1) * t870 - g(2) * t867;
t864 = -g(3) + qJDD(1);
t877 = cos(qJ(2));
t871 = cos(pkin(6));
t874 = sin(qJ(2));
t901 = t871 * t874;
t868 = sin(pkin(6));
t902 = t868 * t874;
t819 = t855 * t901 + t877 * t856 + t864 * t902;
t879 = qJD(2) ^ 2;
t810 = -pkin(2) * t879 + qJDD(2) * pkin(8) + t819;
t836 = -t855 * t868 + t864 * t871;
t873 = sin(qJ(3));
t876 = cos(qJ(3));
t797 = -t873 * t810 + t876 * t836;
t896 = qJD(2) * qJD(3);
t895 = t876 * t896;
t853 = qJDD(2) * t873 + t895;
t792 = (-t853 + t895) * qJ(4) + (t873 * t876 * t879 + qJDD(3)) * pkin(3) + t797;
t798 = t876 * t810 + t873 * t836;
t854 = qJDD(2) * t876 - t873 * t896;
t899 = qJD(2) * t873;
t857 = qJD(3) * pkin(3) - qJ(4) * t899;
t863 = t876 ^ 2;
t793 = -pkin(3) * t863 * t879 + qJ(4) * t854 - qJD(3) * t857 + t798;
t866 = sin(pkin(11));
t904 = cos(pkin(11));
t841 = (t866 * t876 + t873 * t904) * qJD(2);
t776 = t792 * t904 - t866 * t793 + t841 * t907;
t818 = -t874 * t856 + (t855 * t871 + t864 * t868) * t877;
t898 = qJD(2) * t876;
t840 = t866 * t899 - t904 * t898;
t777 = t866 * t792 + t904 * t793 + t840 * t907;
t821 = pkin(4) * t840 - qJ(5) * t841;
t878 = qJD(3) ^ 2;
t775 = -pkin(4) * t878 + qJDD(3) * qJ(5) - t821 * t840 + t777;
t883 = -qJDD(2) * pkin(2) - t818;
t794 = -t854 * pkin(3) + qJDD(4) + t857 * t899 + (-qJ(4) * t863 - pkin(8)) * t879 + t883;
t825 = t853 * t866 - t904 * t854;
t826 = t853 * t904 + t866 * t854;
t780 = (qJD(3) * t840 - t826) * qJ(5) + (qJD(3) * t841 + t825) * pkin(4) + t794;
t865 = sin(pkin(12));
t869 = cos(pkin(12));
t831 = qJD(3) * t865 + t841 * t869;
t770 = -0.2e1 * qJD(5) * t831 - t865 * t775 + t869 * t780;
t814 = qJDD(3) * t865 + t826 * t869;
t830 = qJD(3) * t869 - t841 * t865;
t768 = (t830 * t840 - t814) * pkin(9) + (t830 * t831 + t825) * pkin(5) + t770;
t771 = 0.2e1 * qJD(5) * t830 + t869 * t775 + t865 * t780;
t811 = pkin(5) * t840 - pkin(9) * t831;
t813 = qJDD(3) * t869 - t826 * t865;
t829 = t830 ^ 2;
t769 = -pkin(5) * t829 + pkin(9) * t813 - t811 * t840 + t771;
t872 = sin(qJ(6));
t875 = cos(qJ(6));
t766 = t768 * t875 - t769 * t872;
t803 = t830 * t875 - t831 * t872;
t783 = qJD(6) * t803 + t813 * t872 + t814 * t875;
t804 = t830 * t872 + t831 * t875;
t788 = -mrSges(7,1) * t803 + mrSges(7,2) * t804;
t839 = qJD(6) + t840;
t795 = -mrSges(7,2) * t839 + mrSges(7,3) * t803;
t824 = qJDD(6) + t825;
t763 = m(7) * t766 + mrSges(7,1) * t824 - mrSges(7,3) * t783 - t788 * t804 + t795 * t839;
t767 = t768 * t872 + t769 * t875;
t782 = -qJD(6) * t804 + t813 * t875 - t814 * t872;
t796 = mrSges(7,1) * t839 - mrSges(7,3) * t804;
t764 = m(7) * t767 - mrSges(7,2) * t824 + mrSges(7,3) * t782 + t788 * t803 - t796 * t839;
t755 = t875 * t763 + t872 * t764;
t805 = -mrSges(6,1) * t830 + mrSges(6,2) * t831;
t889 = -mrSges(6,2) * t840 + mrSges(6,3) * t830;
t753 = m(6) * t770 + t825 * mrSges(6,1) - t814 * mrSges(6,3) - t831 * t805 + t840 * t889 + t755;
t808 = mrSges(6,1) * t840 - mrSges(6,3) * t831;
t891 = -t763 * t872 + t875 * t764;
t754 = m(6) * t771 - mrSges(6,2) * t825 + mrSges(6,3) * t813 + t805 * t830 - t808 * t840 + t891;
t749 = -t753 * t865 + t869 * t754;
t822 = mrSges(5,1) * t840 + mrSges(5,2) * t841;
t835 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t841;
t746 = m(5) * t777 - qJDD(3) * mrSges(5,2) - mrSges(5,3) * t825 - qJD(3) * t835 - t822 * t840 + t749;
t774 = -qJDD(3) * pkin(4) - t878 * qJ(5) + t841 * t821 + qJDD(5) - t776;
t772 = -t813 * pkin(5) - t829 * pkin(9) + t831 * t811 + t774;
t884 = m(7) * t772 - t782 * mrSges(7,1) + mrSges(7,2) * t783 - t803 * t795 + t796 * t804;
t765 = m(6) * t774 - t813 * mrSges(6,1) + mrSges(6,2) * t814 + t808 * t831 - t830 * t889 + t884;
t834 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t840;
t759 = m(5) * t776 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t826 + qJD(3) * t834 - t822 * t841 - t765;
t738 = t866 * t746 + t904 * t759;
t784 = Ifges(7,5) * t804 + Ifges(7,6) * t803 + Ifges(7,3) * t839;
t786 = Ifges(7,1) * t804 + Ifges(7,4) * t803 + Ifges(7,5) * t839;
t756 = -mrSges(7,1) * t772 + mrSges(7,3) * t767 + Ifges(7,4) * t783 + Ifges(7,2) * t782 + Ifges(7,6) * t824 - t784 * t804 + t786 * t839;
t785 = Ifges(7,4) * t804 + Ifges(7,2) * t803 + Ifges(7,6) * t839;
t757 = mrSges(7,2) * t772 - mrSges(7,3) * t766 + Ifges(7,1) * t783 + Ifges(7,4) * t782 + Ifges(7,5) * t824 + t784 * t803 - t785 * t839;
t799 = Ifges(6,5) * t831 + Ifges(6,6) * t830 + Ifges(6,3) * t840;
t801 = Ifges(6,1) * t831 + Ifges(6,4) * t830 + Ifges(6,5) * t840;
t739 = -mrSges(6,1) * t774 + mrSges(6,3) * t771 + Ifges(6,4) * t814 + Ifges(6,2) * t813 + Ifges(6,6) * t825 - pkin(5) * t884 + pkin(9) * t891 + t875 * t756 + t872 * t757 - t831 * t799 + t840 * t801;
t800 = Ifges(6,4) * t831 + Ifges(6,2) * t830 + Ifges(6,6) * t840;
t740 = mrSges(6,2) * t774 - mrSges(6,3) * t770 + Ifges(6,1) * t814 + Ifges(6,4) * t813 + Ifges(6,5) * t825 - pkin(9) * t755 - t756 * t872 + t757 * t875 + t799 * t830 - t800 * t840;
t816 = Ifges(5,4) * t841 - Ifges(5,2) * t840 + Ifges(5,6) * qJD(3);
t817 = Ifges(5,1) * t841 - Ifges(5,4) * t840 + Ifges(5,5) * qJD(3);
t844 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t873 + Ifges(4,2) * t876) * qJD(2);
t845 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t873 + Ifges(4,4) * t876) * qJD(2);
t906 = (t844 * t873 - t845 * t876) * qJD(2) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + mrSges(4,1) * t797 + mrSges(5,1) * t776 - mrSges(4,2) * t798 - mrSges(5,2) * t777 + Ifges(4,5) * t853 + Ifges(5,5) * t826 + Ifges(4,6) * t854 - Ifges(5,6) * t825 + pkin(3) * t738 - pkin(4) * t765 + qJ(5) * t749 + t869 * t739 + t865 * t740 + t841 * t816 + t840 * t817;
t748 = t869 * t753 + t865 * t754;
t747 = m(5) * t794 + t825 * mrSges(5,1) + mrSges(5,2) * t826 + t840 * t834 + t835 * t841 + t748;
t809 = -t879 * pkin(8) + t883;
t858 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t899;
t859 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t898;
t881 = -m(4) * t809 + t854 * mrSges(4,1) - mrSges(4,2) * t853 - t858 * t899 + t859 * t898 - t747;
t743 = m(3) * t818 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t879 + t881;
t903 = t743 * t877;
t852 = (-mrSges(4,1) * t876 + mrSges(4,2) * t873) * qJD(2);
t736 = m(4) * t797 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t853 + qJD(3) * t859 - t852 * t899 + t738;
t892 = t904 * t746 - t759 * t866;
t737 = m(4) * t798 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t854 - qJD(3) * t858 + t852 * t898 + t892;
t893 = -t736 * t873 + t876 * t737;
t728 = m(3) * t819 - mrSges(3,1) * t879 - qJDD(2) * mrSges(3,2) + t893;
t731 = t876 * t736 + t873 * t737;
t730 = m(3) * t836 + t731;
t719 = t728 * t901 - t730 * t868 + t871 * t903;
t717 = m(2) * t855 + t719;
t723 = t877 * t728 - t743 * t874;
t722 = m(2) * t856 + t723;
t900 = t870 * t717 + t867 * t722;
t718 = t728 * t902 + t871 * t730 + t868 * t903;
t894 = -t717 * t867 + t870 * t722;
t890 = m(2) * t864 + t718;
t815 = Ifges(5,5) * t841 - Ifges(5,6) * t840 + Ifges(5,3) * qJD(3);
t724 = mrSges(5,2) * t794 - mrSges(5,3) * t776 + Ifges(5,1) * t826 - Ifges(5,4) * t825 + Ifges(5,5) * qJDD(3) - qJ(5) * t748 - qJD(3) * t816 - t739 * t865 + t740 * t869 - t815 * t840;
t882 = mrSges(7,1) * t766 - mrSges(7,2) * t767 + Ifges(7,5) * t783 + Ifges(7,6) * t782 + Ifges(7,3) * t824 + t804 * t785 - t803 * t786;
t732 = Ifges(5,6) * qJDD(3) - t882 - pkin(5) * t755 - pkin(4) * t748 + (-Ifges(5,2) - Ifges(6,3)) * t825 - mrSges(6,1) * t770 + mrSges(6,2) * t771 + mrSges(5,3) * t777 - mrSges(5,1) * t794 - Ifges(6,6) * t813 - Ifges(6,5) * t814 + qJD(3) * t817 + Ifges(5,4) * t826 + t830 * t801 - t831 * t800 - t841 * t815;
t843 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t873 + Ifges(4,6) * t876) * qJD(2);
t713 = -mrSges(4,1) * t809 + mrSges(4,3) * t798 + Ifges(4,4) * t853 + Ifges(4,2) * t854 + Ifges(4,6) * qJDD(3) - pkin(3) * t747 + qJ(4) * t892 + qJD(3) * t845 + t866 * t724 + t732 * t904 - t843 * t899;
t715 = mrSges(4,2) * t809 - mrSges(4,3) * t797 + Ifges(4,1) * t853 + Ifges(4,4) * t854 + Ifges(4,5) * qJDD(3) - qJ(4) * t738 - qJD(3) * t844 + t724 * t904 - t866 * t732 + t843 * t898;
t712 = mrSges(3,2) * t836 - mrSges(3,3) * t818 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t879 - pkin(8) * t731 - t713 * t873 + t715 * t876;
t714 = -mrSges(3,1) * t836 + mrSges(3,3) * t819 + t879 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t731 - t906;
t885 = pkin(7) * t723 + t712 * t874 + t714 * t877;
t711 = mrSges(3,1) * t818 - mrSges(3,2) * t819 + Ifges(3,3) * qJDD(2) + pkin(2) * t881 + pkin(8) * t893 + t876 * t713 + t873 * t715;
t710 = mrSges(2,2) * t864 - mrSges(2,3) * t855 + t877 * t712 - t874 * t714 + (-t718 * t868 - t719 * t871) * pkin(7);
t709 = -mrSges(2,1) * t864 + mrSges(2,3) * t856 - pkin(1) * t718 - t868 * t711 + t871 * t885;
t1 = [-m(1) * g(1) + t894; -m(1) * g(2) + t900; -m(1) * g(3) + t890; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t900 - t867 * t709 + t870 * t710; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t894 + t870 * t709 + t867 * t710; -mrSges(1,1) * g(2) + mrSges(2,1) * t855 + mrSges(1,2) * g(1) - mrSges(2,2) * t856 + pkin(1) * t719 + t871 * t711 + t868 * t885; t890; t711; t906; t747; t765; t882;];
tauJB  = t1;
