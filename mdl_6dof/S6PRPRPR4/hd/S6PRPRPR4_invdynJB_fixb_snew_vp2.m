% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:46:42
% EndTime: 2019-05-04 22:47:02
% DurationCPUTime: 19.45s
% Computational Cost: add. (308191->320), mult. (708386->414), div. (0->0), fcn. (538212->14), ass. (0->148)
t882 = qJD(2) ^ 2;
t870 = sin(pkin(10));
t874 = cos(pkin(10));
t854 = g(1) * t870 - g(2) * t874;
t855 = -g(1) * t874 - g(2) * t870;
t867 = -g(3) + qJDD(1);
t871 = sin(pkin(6));
t875 = cos(pkin(6));
t878 = sin(qJ(2));
t880 = cos(qJ(2));
t827 = -t878 * t855 + (t854 * t875 + t867 * t871) * t880;
t917 = cos(qJ(4));
t873 = cos(pkin(11));
t916 = pkin(3) * t873;
t869 = sin(pkin(11));
t915 = mrSges(4,2) * t869;
t887 = qJDD(3) - t827;
t817 = -qJDD(2) * pkin(2) - t882 * qJ(3) + t887;
t865 = t869 ^ 2;
t911 = t875 * t878;
t912 = t871 * t878;
t828 = t854 * t911 + t880 * t855 + t867 * t912;
t819 = -pkin(2) * t882 + qJDD(2) * qJ(3) + t828;
t845 = -t854 * t871 + t867 * t875;
t905 = qJD(2) * qJD(3);
t909 = t873 * t845 - 0.2e1 * t869 * t905;
t801 = (-pkin(8) * qJDD(2) + t882 * t916 - t819) * t869 + t909;
t804 = t869 * t845 + (t819 + 0.2e1 * t905) * t873;
t903 = qJDD(2) * t873;
t866 = t873 ^ 2;
t913 = t866 * t882;
t802 = -pkin(3) * t913 + pkin(8) * t903 + t804;
t877 = sin(qJ(4));
t786 = t877 * t801 + t917 * t802;
t902 = t873 * t917;
t908 = qJD(2) * t869;
t847 = -qJD(2) * t902 + t877 * t908;
t890 = t869 * t917 + t873 * t877;
t848 = t890 * qJD(2);
t830 = pkin(4) * t847 - qJ(5) * t848;
t881 = qJD(4) ^ 2;
t784 = -pkin(4) * t881 + qJDD(4) * qJ(5) - t830 * t847 + t786;
t811 = (-pkin(2) - t916) * qJDD(2) + (-qJ(3) + (-t865 - t866) * pkin(8)) * t882 + t887;
t904 = qJDD(2) * t869;
t907 = qJD(4) * t848;
t834 = -qJDD(2) * t902 + t877 * t904 + t907;
t906 = t847 * qJD(4);
t835 = qJDD(2) * t890 - t906;
t789 = (-t835 + t906) * qJ(5) + (t834 + t907) * pkin(4) + t811;
t868 = sin(pkin(12));
t872 = cos(pkin(12));
t840 = qJD(4) * t868 + t848 * t872;
t779 = -0.2e1 * qJD(5) * t840 - t868 * t784 + t872 * t789;
t823 = qJDD(4) * t868 + t835 * t872;
t839 = qJD(4) * t872 - t848 * t868;
t777 = (t839 * t847 - t823) * pkin(9) + (t839 * t840 + t834) * pkin(5) + t779;
t780 = 0.2e1 * qJD(5) * t839 + t872 * t784 + t868 * t789;
t821 = pkin(5) * t847 - pkin(9) * t840;
t822 = qJDD(4) * t872 - t835 * t868;
t838 = t839 ^ 2;
t778 = -pkin(5) * t838 + pkin(9) * t822 - t821 * t847 + t780;
t876 = sin(qJ(6));
t879 = cos(qJ(6));
t775 = t777 * t879 - t778 * t876;
t812 = t839 * t879 - t840 * t876;
t792 = qJD(6) * t812 + t822 * t876 + t823 * t879;
t813 = t839 * t876 + t840 * t879;
t797 = -mrSges(7,1) * t812 + mrSges(7,2) * t813;
t846 = qJD(6) + t847;
t805 = -mrSges(7,2) * t846 + mrSges(7,3) * t812;
t833 = qJDD(6) + t834;
t772 = m(7) * t775 + mrSges(7,1) * t833 - t792 * mrSges(7,3) - t797 * t813 + t805 * t846;
t776 = t777 * t876 + t778 * t879;
t791 = -qJD(6) * t813 + t822 * t879 - t823 * t876;
t806 = mrSges(7,1) * t846 - mrSges(7,3) * t813;
t773 = m(7) * t776 - mrSges(7,2) * t833 + t791 * mrSges(7,3) + t797 * t812 - t806 * t846;
t764 = t879 * t772 + t876 * t773;
t814 = -mrSges(6,1) * t839 + mrSges(6,2) * t840;
t896 = -mrSges(6,2) * t847 + mrSges(6,3) * t839;
t762 = m(6) * t779 + t834 * mrSges(6,1) - t823 * mrSges(6,3) - t840 * t814 + t847 * t896 + t764;
t820 = mrSges(6,1) * t847 - mrSges(6,3) * t840;
t898 = -t772 * t876 + t879 * t773;
t763 = m(6) * t780 - mrSges(6,2) * t834 + mrSges(6,3) * t822 + t814 * t839 - t820 * t847 + t898;
t757 = t872 * t762 + t868 * t763;
t843 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t847;
t844 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t848;
t886 = m(5) * t811 + t834 * mrSges(5,1) + t835 * mrSges(5,2) + t847 * t843 + t848 * t844 + t757;
t885 = -m(4) * t817 + mrSges(4,1) * t903 - t886 + (t865 * t882 + t913) * mrSges(4,3);
t752 = (mrSges(3,1) - t915) * qJDD(2) + t885 - t882 * mrSges(3,2) + m(3) * t827;
t914 = t752 * t880;
t758 = -t762 * t868 + t872 * t763;
t831 = mrSges(5,1) * t847 + mrSges(5,2) * t848;
t755 = m(5) * t786 - qJDD(4) * mrSges(5,2) - mrSges(5,3) * t834 - qJD(4) * t844 - t831 * t847 + t758;
t785 = t801 * t917 - t877 * t802;
t783 = -qJDD(4) * pkin(4) - t881 * qJ(5) + t848 * t830 + qJDD(5) - t785;
t781 = -t822 * pkin(5) - t838 * pkin(9) + t840 * t821 + t783;
t888 = m(7) * t781 - t791 * mrSges(7,1) + t792 * mrSges(7,2) - t812 * t805 + t806 * t813;
t774 = m(6) * t783 - t822 * mrSges(6,1) + mrSges(6,2) * t823 + t820 * t840 - t839 * t896 + t888;
t768 = m(5) * t785 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t835 + qJD(4) * t843 - t831 * t848 - t774;
t747 = t877 * t755 + t917 * t768;
t803 = -t819 * t869 + t909;
t891 = mrSges(4,3) * qJDD(2) + t882 * (-mrSges(4,1) * t873 + t915);
t745 = m(4) * t803 - t869 * t891 + t747;
t899 = t917 * t755 - t877 * t768;
t746 = m(4) * t804 + t873 * t891 + t899;
t900 = -t745 * t869 + t873 * t746;
t737 = m(3) * t828 - mrSges(3,1) * t882 - qJDD(2) * mrSges(3,2) + t900;
t740 = t873 * t745 + t869 * t746;
t739 = m(3) * t845 + t740;
t728 = t737 * t911 - t739 * t871 + t875 * t914;
t726 = m(2) * t854 + t728;
t732 = t880 * t737 - t752 * t878;
t731 = m(2) * t855 + t732;
t910 = t874 * t726 + t870 * t731;
t727 = t737 * t912 + t875 * t739 + t871 * t914;
t901 = -t726 * t870 + t874 * t731;
t897 = m(2) * t867 + t727;
t895 = Ifges(4,1) * t869 + Ifges(4,4) * t873;
t894 = Ifges(4,4) * t869 + Ifges(4,2) * t873;
t893 = Ifges(4,5) * t869 + Ifges(4,6) * t873;
t793 = Ifges(7,5) * t813 + Ifges(7,6) * t812 + Ifges(7,3) * t846;
t795 = Ifges(7,1) * t813 + Ifges(7,4) * t812 + Ifges(7,5) * t846;
t765 = -mrSges(7,1) * t781 + mrSges(7,3) * t776 + Ifges(7,4) * t792 + Ifges(7,2) * t791 + Ifges(7,6) * t833 - t793 * t813 + t795 * t846;
t794 = Ifges(7,4) * t813 + Ifges(7,2) * t812 + Ifges(7,6) * t846;
t766 = mrSges(7,2) * t781 - mrSges(7,3) * t775 + Ifges(7,1) * t792 + Ifges(7,4) * t791 + Ifges(7,5) * t833 + t793 * t812 - t794 * t846;
t807 = Ifges(6,5) * t840 + Ifges(6,6) * t839 + Ifges(6,3) * t847;
t809 = Ifges(6,1) * t840 + Ifges(6,4) * t839 + Ifges(6,5) * t847;
t748 = -mrSges(6,1) * t783 + mrSges(6,3) * t780 + Ifges(6,4) * t823 + Ifges(6,2) * t822 + Ifges(6,6) * t834 - pkin(5) * t888 + pkin(9) * t898 + t879 * t765 + t876 * t766 - t840 * t807 + t847 * t809;
t808 = Ifges(6,4) * t840 + Ifges(6,2) * t839 + Ifges(6,6) * t847;
t749 = mrSges(6,2) * t783 - mrSges(6,3) * t779 + Ifges(6,1) * t823 + Ifges(6,4) * t822 + Ifges(6,5) * t834 - pkin(9) * t764 - t765 * t876 + t766 * t879 + t807 * t839 - t808 * t847;
t824 = Ifges(5,5) * t848 - Ifges(5,6) * t847 + Ifges(5,3) * qJD(4);
t825 = Ifges(5,4) * t848 - Ifges(5,2) * t847 + Ifges(5,6) * qJD(4);
t733 = mrSges(5,2) * t811 - mrSges(5,3) * t785 + Ifges(5,1) * t835 - Ifges(5,4) * t834 + Ifges(5,5) * qJDD(4) - qJ(5) * t757 - qJD(4) * t825 - t748 * t868 + t749 * t872 - t824 * t847;
t826 = Ifges(5,1) * t848 - Ifges(5,4) * t847 + Ifges(5,5) * qJD(4);
t884 = mrSges(7,1) * t775 - mrSges(7,2) * t776 + Ifges(7,5) * t792 + Ifges(7,6) * t791 + Ifges(7,3) * t833 + t813 * t794 - t812 * t795;
t741 = Ifges(5,6) * qJDD(4) + (-Ifges(5,2) - Ifges(6,3)) * t834 - t884 - t848 * t824 + t839 * t809 - t840 * t808 + Ifges(5,4) * t835 - Ifges(6,6) * t822 - Ifges(6,5) * t823 + qJD(4) * t826 - mrSges(5,1) * t811 + mrSges(5,3) * t786 - mrSges(6,1) * t779 + mrSges(6,2) * t780 - pkin(5) * t764 - pkin(4) * t757;
t853 = t893 * qJD(2);
t722 = -mrSges(4,1) * t817 + mrSges(4,3) * t804 - pkin(3) * t886 + pkin(8) * t899 + qJDD(2) * t894 + t877 * t733 + t741 * t917 - t853 * t908;
t724 = t873 * qJD(2) * t853 + mrSges(4,2) * t817 - mrSges(4,3) * t803 - pkin(8) * t747 + qJDD(2) * t895 + t733 * t917 - t877 * t741;
t721 = mrSges(3,2) * t845 - mrSges(3,3) * t827 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t882 - qJ(3) * t740 - t722 * t869 + t724 * t873;
t883 = mrSges(5,1) * t785 - mrSges(5,2) * t786 + Ifges(5,5) * t835 - Ifges(5,6) * t834 + Ifges(5,3) * qJDD(4) - pkin(4) * t774 + qJ(5) * t758 + t872 * t748 + t868 * t749 + t848 * t825 + t847 * t826;
t723 = -t883 + (Ifges(3,6) - t893) * qJDD(2) - mrSges(3,1) * t845 + mrSges(3,3) * t828 - mrSges(4,1) * t803 + mrSges(4,2) * t804 - pkin(3) * t747 - pkin(2) * t740 + (-t869 * t894 + t873 * t895 + Ifges(3,5)) * t882;
t889 = pkin(7) * t732 + t721 * t878 + t723 * t880;
t756 = mrSges(4,2) * t904 - t885;
t720 = mrSges(3,1) * t827 - mrSges(3,2) * t828 + Ifges(3,3) * qJDD(2) - pkin(2) * t756 + qJ(3) * t900 + t873 * t722 + t869 * t724;
t719 = mrSges(2,2) * t867 - mrSges(2,3) * t854 + t880 * t721 - t878 * t723 + (-t727 * t871 - t728 * t875) * pkin(7);
t718 = -mrSges(2,1) * t867 + mrSges(2,3) * t855 - pkin(1) * t727 - t871 * t720 + t875 * t889;
t1 = [-m(1) * g(1) + t901; -m(1) * g(2) + t910; -m(1) * g(3) + t897; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t910 - t870 * t718 + t874 * t719; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t901 + t874 * t718 + t870 * t719; -mrSges(1,1) * g(2) + mrSges(2,1) * t854 + mrSges(1,2) * g(1) - mrSges(2,2) * t855 + pkin(1) * t728 + t875 * t720 + t871 * t889; t897; t720; t756; t883; t774; t884;];
tauJB  = t1;
