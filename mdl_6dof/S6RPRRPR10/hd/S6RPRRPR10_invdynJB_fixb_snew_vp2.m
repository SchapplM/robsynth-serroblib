% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-05-05 23:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:53:33
% EndTime: 2019-05-05 23:53:40
% DurationCPUTime: 4.54s
% Computational Cost: add. (52283->318), mult. (99870->369), div. (0->0), fcn. (60417->8), ass. (0->132)
t915 = Ifges(5,1) + Ifges(6,1);
t903 = Ifges(5,4) - Ifges(6,5);
t901 = Ifges(5,5) + Ifges(6,4);
t914 = -Ifges(5,2) - Ifges(6,3);
t900 = Ifges(5,6) - Ifges(6,6);
t913 = Ifges(5,3) + Ifges(6,2);
t862 = sin(qJ(1));
t865 = cos(qJ(1));
t843 = -t865 * g(1) - t862 * g(2);
t912 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t843;
t860 = sin(qJ(4));
t864 = cos(qJ(3));
t892 = qJD(1) * t864;
t907 = cos(qJ(4));
t832 = -qJD(3) * t907 + t860 * t892;
t861 = sin(qJ(3));
t890 = qJD(1) * qJD(3);
t888 = t861 * t890;
t837 = t864 * qJDD(1) - t888;
t790 = -t832 * qJD(4) + t860 * qJDD(3) + t837 * t907;
t842 = t862 * g(1) - t865 * g(2);
t867 = qJD(1) ^ 2;
t878 = -t867 * qJ(2) + qJDD(2) - t842;
t908 = -pkin(1) - pkin(7);
t810 = qJDD(1) * t908 + t878;
t799 = t861 * g(3) + t864 * t810;
t835 = (pkin(3) * t861 - pkin(8) * t864) * qJD(1);
t866 = qJD(3) ^ 2;
t877 = qJDD(3) * pkin(3) + t866 * pkin(8) - t835 * t892 + t799;
t891 = t861 * qJD(1);
t846 = qJD(4) + t891;
t898 = t832 * t846;
t911 = (-t790 + t898) * qJ(5) - t877;
t833 = t860 * qJD(3) + t892 * t907;
t797 = t832 * mrSges(6,1) - t833 * mrSges(6,3);
t809 = t867 * t908 - t912;
t887 = t864 * t890;
t836 = -t861 * qJDD(1) - t887;
t770 = (-t837 + t888) * pkin(8) + (-t836 + t887) * pkin(3) + t809;
t800 = -t864 * g(3) + t861 * t810;
t774 = -t866 * pkin(3) + qJDD(3) * pkin(8) - t835 * t891 + t800;
t754 = t770 * t907 - t860 * t774;
t796 = t832 * pkin(4) - t833 * qJ(5);
t831 = qJDD(4) - t836;
t845 = t846 ^ 2;
t752 = -t831 * pkin(4) - t845 * qJ(5) + t833 * t796 + qJDD(5) - t754;
t746 = (-t790 - t898) * pkin(9) + (t832 * t833 - t831) * pkin(5) + t752;
t755 = t860 * t770 + t907 * t774;
t909 = 2 * qJD(5);
t751 = -t845 * pkin(4) + t831 * qJ(5) - t832 * t796 + t846 * t909 + t755;
t789 = t833 * qJD(4) - qJDD(3) * t907 + t860 * t837;
t808 = -t846 * pkin(5) - t833 * pkin(9);
t830 = t832 ^ 2;
t747 = -t830 * pkin(5) + t789 * pkin(9) + t846 * t808 + t751;
t859 = sin(qJ(6));
t863 = cos(qJ(6));
t744 = t863 * t746 - t859 * t747;
t791 = t863 * t832 - t859 * t833;
t762 = t791 * qJD(6) + t859 * t789 + t863 * t790;
t792 = t859 * t832 + t863 * t833;
t768 = -t791 * mrSges(7,1) + t792 * mrSges(7,2);
t844 = qJD(6) - t846;
t775 = -t844 * mrSges(7,2) + t791 * mrSges(7,3);
t824 = qJDD(6) - t831;
t741 = m(7) * t744 + t824 * mrSges(7,1) - t762 * mrSges(7,3) - t792 * t768 + t844 * t775;
t745 = t859 * t746 + t863 * t747;
t761 = -t792 * qJD(6) + t863 * t789 - t859 * t790;
t776 = t844 * mrSges(7,1) - t792 * mrSges(7,3);
t742 = m(7) * t745 - t824 * mrSges(7,2) + t761 * mrSges(7,3) + t791 * t768 - t844 * t776;
t734 = t863 * t741 + t859 * t742;
t804 = -t832 * mrSges(6,2) + t846 * mrSges(6,3);
t874 = -m(6) * t752 + t831 * mrSges(6,1) + t846 * t804 - t734;
t733 = t790 * mrSges(6,2) + t833 * t797 - t874;
t764 = Ifges(7,4) * t792 + Ifges(7,2) * t791 + Ifges(7,6) * t844;
t765 = Ifges(7,1) * t792 + Ifges(7,4) * t791 + Ifges(7,5) * t844;
t873 = -mrSges(7,1) * t744 + mrSges(7,2) * t745 - Ifges(7,5) * t762 - Ifges(7,6) * t761 - Ifges(7,3) * t824 - t792 * t764 + t791 * t765;
t803 = -t846 * mrSges(6,1) + t833 * mrSges(6,2);
t883 = -t859 * t741 + t863 * t742;
t879 = m(6) * t751 + t831 * mrSges(6,3) + t846 * t803 + t883;
t894 = -t903 * t832 + t915 * t833 + t901 * t846;
t895 = t914 * t832 + t903 * t833 + t900 * t846;
t910 = -t789 * t900 + t790 * t901 + t913 * t831 + t832 * t894 + t833 * t895 + mrSges(5,1) * t754 - mrSges(6,1) * t752 - mrSges(5,2) * t755 + mrSges(6,3) * t751 - pkin(4) * t733 - pkin(5) * t734 + qJ(5) * (-t789 * mrSges(6,2) - t832 * t797 + t879) + t873;
t906 = mrSges(2,1) - mrSges(3,2);
t905 = -mrSges(5,3) - mrSges(6,2);
t904 = -Ifges(3,4) + Ifges(2,5);
t902 = (Ifges(3,5) - Ifges(2,6));
t834 = (mrSges(4,1) * t861 + mrSges(4,2) * t864) * qJD(1);
t840 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t892;
t802 = t846 * mrSges(5,1) - t833 * mrSges(5,3);
t893 = -t832 * mrSges(5,1) - t833 * mrSges(5,2) - t797;
t729 = m(5) * t755 - t831 * mrSges(5,2) + t789 * t905 - t846 * t802 + t832 * t893 + t879;
t801 = -t846 * mrSges(5,2) - t832 * mrSges(5,3);
t731 = m(5) * t754 + t831 * mrSges(5,1) + t790 * t905 + t846 * t801 + t833 * t893 + t874;
t884 = t907 * t729 - t860 * t731;
t724 = m(4) * t800 - qJDD(3) * mrSges(4,2) + t836 * mrSges(4,3) - qJD(3) * t840 - t834 * t891 + t884;
t839 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t891;
t753 = -0.2e1 * qJD(5) * t833 + (t833 * t846 + t789) * pkin(4) + t911;
t749 = -t830 * pkin(9) + (-pkin(4) - pkin(5)) * t789 + (-pkin(4) * t846 + t808 + t909) * t833 - t911;
t881 = -m(7) * t749 + t761 * mrSges(7,1) - t762 * mrSges(7,2) + t791 * t775 - t792 * t776;
t739 = m(6) * t753 + t789 * mrSges(6,1) - t790 * mrSges(6,3) - t833 * t803 + t832 * t804 + t881;
t869 = m(5) * t877 - t789 * mrSges(5,1) - t790 * mrSges(5,2) - t832 * t801 - t833 * t802 - t739;
t737 = m(4) * t799 + qJDD(3) * mrSges(4,1) - t837 * mrSges(4,3) + qJD(3) * t839 - t834 * t892 + t869;
t718 = t861 * t724 + t864 * t737;
t815 = -qJDD(1) * pkin(1) + t878;
t876 = -m(3) * t815 + (t867 * mrSges(3,3)) - t718;
t714 = m(2) * t842 - (t867 * mrSges(2,2)) + qJDD(1) * t906 + t876;
t813 = t867 * pkin(1) + t912;
t726 = t860 * t729 + t907 * t731;
t875 = -m(4) * t809 + t836 * mrSges(4,1) - t837 * mrSges(4,2) - t839 * t891 - t840 * t892 - t726;
t871 = -m(3) * t813 + (t867 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t875;
t721 = m(2) * t843 - (t867 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t871;
t897 = t865 * t714 + t862 * t721;
t896 = t900 * t832 - t901 * t833 - t913 * t846;
t886 = -t862 * t714 + t865 * t721;
t885 = t864 * t724 - t861 * t737;
t763 = Ifges(7,5) * t792 + Ifges(7,6) * t791 + Ifges(7,3) * t844;
t735 = -mrSges(7,1) * t749 + mrSges(7,3) * t745 + Ifges(7,4) * t762 + Ifges(7,2) * t761 + Ifges(7,6) * t824 - t792 * t763 + t844 * t765;
t736 = mrSges(7,2) * t749 - mrSges(7,3) * t744 + Ifges(7,1) * t762 + Ifges(7,4) * t761 + Ifges(7,5) * t824 + t791 * t763 - t844 * t764;
t710 = mrSges(5,1) * t877 - mrSges(6,1) * t753 + mrSges(6,2) * t751 + mrSges(5,3) * t755 - pkin(4) * t739 - pkin(5) * t881 - pkin(9) * t883 - t863 * t735 - t859 * t736 + t914 * t789 + t903 * t790 + t900 * t831 + t896 * t833 + t894 * t846;
t712 = -mrSges(5,2) * t877 + mrSges(6,2) * t752 - mrSges(5,3) * t754 - mrSges(6,3) * t753 - pkin(9) * t734 - qJ(5) * t739 - t859 * t735 + t863 * t736 - t903 * t789 + t915 * t790 + t901 * t831 + t896 * t832 - t895 * t846;
t819 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t864 - Ifges(4,2) * t861) * qJD(1);
t820 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t864 - Ifges(4,4) * t861) * qJD(1);
t872 = mrSges(4,1) * t799 - mrSges(4,2) * t800 + Ifges(4,5) * t837 + Ifges(4,6) * t836 + Ifges(4,3) * qJDD(3) + pkin(3) * t869 + pkin(8) * t884 + t907 * t710 + t860 * t712 + t819 * t892 + t820 * t891;
t818 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t864 - Ifges(4,6) * t861) * qJD(1);
t707 = mrSges(4,2) * t809 - mrSges(4,3) * t799 + Ifges(4,1) * t837 + Ifges(4,4) * t836 + Ifges(4,5) * qJDD(3) - pkin(8) * t726 - qJD(3) * t819 - t860 * t710 + t712 * t907 - t818 * t891;
t708 = -mrSges(4,1) * t809 + mrSges(4,3) * t800 + Ifges(4,4) * t837 + Ifges(4,2) * t836 + Ifges(4,6) * qJDD(3) - pkin(3) * t726 + qJD(3) * t820 - t818 * t892 - t910;
t716 = qJDD(1) * mrSges(3,2) - t876;
t870 = mrSges(2,1) * t842 - mrSges(2,2) * t843 + mrSges(3,2) * t815 - mrSges(3,3) * t813 - pkin(1) * t716 - pkin(7) * t718 + qJ(2) * t871 + t864 * t707 - t861 * t708 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t717 = -m(3) * g(3) + t885;
t705 = t872 - mrSges(2,3) * t842 + mrSges(3,1) * t815 + t904 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + pkin(2) * t718 - qJ(2) * t717 + (t902 * t867);
t704 = -mrSges(3,1) * t813 + mrSges(2,3) * t843 - pkin(1) * t717 - pkin(2) * t875 - pkin(7) * t885 + g(3) * t906 - qJDD(1) * t902 - t861 * t707 - t864 * t708 + t867 * t904;
t1 = [-m(1) * g(1) + t886; -m(1) * g(2) + t897; (-m(1) - m(2) - m(3)) * g(3) + t885; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t897 - t862 * t704 + t865 * t705; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t886 + t865 * t704 + t862 * t705; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t870; t870; t716; t872; t910; t733; -t873;];
tauJB  = t1;
