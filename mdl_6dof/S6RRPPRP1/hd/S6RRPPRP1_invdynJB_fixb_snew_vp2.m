% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-05-06 09:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:03:16
% EndTime: 2019-05-06 09:03:30
% DurationCPUTime: 13.64s
% Computational Cost: add. (203051->363), mult. (478956->447), div. (0->0), fcn. (338104->10), ass. (0->142)
t969 = -2 * qJD(3);
t968 = Ifges(6,1) + Ifges(7,1);
t960 = Ifges(6,4) - Ifges(7,5);
t959 = -Ifges(6,5) - Ifges(7,4);
t967 = Ifges(6,2) + Ifges(7,3);
t958 = Ifges(6,6) - Ifges(7,6);
t966 = -Ifges(6,3) - Ifges(7,2);
t924 = sin(qJ(2));
t926 = cos(qJ(2));
t945 = qJD(1) * qJD(2);
t908 = qJDD(1) * t924 + t926 * t945;
t925 = sin(qJ(1));
t927 = cos(qJ(1));
t915 = -g(1) * t927 - g(2) * t925;
t929 = qJD(1) ^ 2;
t903 = -pkin(1) * t929 + qJDD(1) * pkin(7) + t915;
t954 = t924 * t903;
t962 = pkin(2) * t929;
t858 = qJDD(2) * pkin(2) - t908 * qJ(3) - t954 + (qJ(3) * t945 + t924 * t962 - g(3)) * t926;
t888 = -g(3) * t924 + t926 * t903;
t909 = qJDD(1) * t926 - t924 * t945;
t948 = qJD(1) * t924;
t911 = qJD(2) * pkin(2) - qJ(3) * t948;
t919 = t926 ^ 2;
t860 = qJ(3) * t909 - qJD(2) * t911 - t919 * t962 + t888;
t921 = sin(pkin(9));
t955 = cos(pkin(9));
t897 = (t921 * t926 + t924 * t955) * qJD(1);
t837 = t858 * t955 - t921 * t860 + t897 * t969;
t947 = qJD(1) * t926;
t896 = t921 * t948 - t955 * t947;
t873 = pkin(3) * t896 - qJ(4) * t897;
t928 = qJD(2) ^ 2;
t820 = -qJDD(2) * pkin(3) - t928 * qJ(4) + t897 * t873 + qJDD(4) - t837;
t920 = sin(pkin(10));
t922 = cos(pkin(10));
t886 = qJD(2) * t920 + t897 * t922;
t864 = pkin(4) * t896 - pkin(8) * t886;
t880 = t908 * t955 + t921 * t909;
t867 = qJDD(2) * t922 - t880 * t920;
t885 = qJD(2) * t922 - t897 * t920;
t884 = t885 ^ 2;
t818 = -t867 * pkin(4) - t884 * pkin(8) + t886 * t864 + t820;
t923 = sin(qJ(5));
t963 = cos(qJ(5));
t855 = t923 * t885 + t886 * t963;
t868 = qJDD(2) * t920 + t880 * t922;
t827 = t855 * qJD(5) - t867 * t963 + t923 * t868;
t854 = -t885 * t963 + t923 * t886;
t828 = -t854 * qJD(5) + t923 * t867 + t868 * t963;
t895 = qJD(5) + t896;
t811 = -0.2e1 * qJD(6) * t855 + (t854 * t895 - t828) * qJ(6) + (t855 * t895 + t827) * pkin(5) + t818;
t844 = -mrSges(7,2) * t854 + mrSges(7,3) * t895;
t847 = -mrSges(7,1) * t895 + mrSges(7,2) * t855;
t805 = m(7) * t811 + t827 * mrSges(7,1) - t828 * mrSges(7,3) + t854 * t844 - t855 * t847;
t838 = t921 * t858 + t955 * t860 + t896 * t969;
t821 = -pkin(3) * t928 + qJDD(2) * qJ(4) - t873 * t896 + t838;
t914 = t925 * g(1) - t927 * g(2);
t935 = -qJDD(1) * pkin(1) - t914;
t862 = -t909 * pkin(2) + qJDD(3) + t911 * t948 + (-qJ(3) * t919 - pkin(7)) * t929 + t935;
t879 = t908 * t921 - t955 * t909;
t824 = (qJD(2) * t896 - t880) * qJ(4) + (qJD(2) * t897 + t879) * pkin(3) + t862;
t816 = -0.2e1 * qJD(4) * t886 - t920 * t821 + t922 * t824;
t813 = (t885 * t896 - t868) * pkin(8) + (t885 * t886 + t879) * pkin(4) + t816;
t817 = 0.2e1 * qJD(4) * t885 + t922 * t821 + t920 * t824;
t815 = -pkin(4) * t884 + pkin(8) * t867 - t864 * t896 + t817;
t810 = t923 * t813 + t815 * t963;
t839 = pkin(5) * t854 - qJ(6) * t855;
t878 = qJDD(5) + t879;
t894 = t895 ^ 2;
t807 = -pkin(5) * t894 + qJ(6) * t878 + 0.2e1 * qJD(6) * t895 - t839 * t854 + t810;
t950 = t960 * t854 - t968 * t855 + t959 * t895;
t951 = t958 * t854 + t959 * t855 + t966 * t895;
t793 = -mrSges(6,1) * t818 - mrSges(7,1) * t811 + mrSges(7,2) * t807 + mrSges(6,3) * t810 - pkin(5) * t805 - t967 * t827 + t960 * t828 + t951 * t855 + t958 * t878 - t950 * t895;
t809 = t813 * t963 - t923 * t815;
t808 = -t878 * pkin(5) - t894 * qJ(6) + t855 * t839 + qJDD(6) - t809;
t952 = t967 * t854 - t960 * t855 - t958 * t895;
t794 = mrSges(6,2) * t818 + mrSges(7,2) * t808 - mrSges(6,3) * t809 - mrSges(7,3) * t811 - qJ(6) * t805 - t960 * t827 + t968 * t828 + t951 * t854 - t959 * t878 + t952 * t895;
t848 = Ifges(5,5) * t886 + Ifges(5,6) * t885 + Ifges(5,3) * t896;
t850 = Ifges(5,1) * t886 + Ifges(5,4) * t885 + Ifges(5,5) * t896;
t845 = -mrSges(6,2) * t895 - mrSges(6,3) * t854;
t846 = mrSges(6,1) * t895 - mrSges(6,3) * t855;
t932 = m(6) * t818 + t827 * mrSges(6,1) + t828 * mrSges(6,2) + t854 * t845 + t855 * t846 + t805;
t944 = m(7) * t807 + t878 * mrSges(7,3) + t895 * t847;
t840 = mrSges(7,1) * t854 - mrSges(7,3) * t855;
t949 = -mrSges(6,1) * t854 - mrSges(6,2) * t855 - t840;
t961 = -mrSges(6,3) - mrSges(7,2);
t797 = m(6) * t810 - t878 * mrSges(6,2) + t827 * t961 - t895 * t846 + t854 * t949 + t944;
t939 = -m(7) * t808 + t878 * mrSges(7,1) + t895 * t844;
t799 = m(6) * t809 + t878 * mrSges(6,1) + t828 * t961 + t895 * t845 + t855 * t949 + t939;
t940 = t797 * t963 - t799 * t923;
t770 = -mrSges(5,1) * t820 + mrSges(5,3) * t817 + Ifges(5,4) * t868 + Ifges(5,2) * t867 + Ifges(5,6) * t879 - pkin(4) * t932 + pkin(8) * t940 + t793 * t963 + t923 * t794 - t886 * t848 + t896 * t850;
t792 = t923 * t797 + t799 * t963;
t849 = Ifges(5,4) * t886 + Ifges(5,2) * t885 + Ifges(5,6) * t896;
t771 = mrSges(5,2) * t820 - mrSges(5,3) * t816 + Ifges(5,1) * t868 + Ifges(5,4) * t867 + Ifges(5,5) * t879 - pkin(8) * t792 - t923 * t793 + t794 * t963 + t885 * t848 - t896 * t849;
t859 = -mrSges(5,1) * t885 + mrSges(5,2) * t886;
t938 = -mrSges(5,2) * t896 + mrSges(5,3) * t885;
t790 = m(5) * t816 + t879 * mrSges(5,1) - t868 * mrSges(5,3) - t886 * t859 + t896 * t938 + t792;
t863 = mrSges(5,1) * t896 - mrSges(5,3) * t886;
t791 = m(5) * t817 - mrSges(5,2) * t879 + mrSges(5,3) * t867 + t859 * t885 - t863 * t896 + t940;
t786 = -t790 * t920 + t922 * t791;
t874 = mrSges(4,1) * t896 + mrSges(4,2) * t897;
t890 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t897;
t783 = m(4) * t838 - qJDD(2) * mrSges(4,2) - mrSges(4,3) * t879 - qJD(2) * t890 - t874 * t896 + t786;
t802 = m(5) * t820 - t867 * mrSges(5,1) + t868 * mrSges(5,2) + t886 * t863 - t885 * t938 + t932;
t889 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t896;
t801 = m(4) * t837 + qJDD(2) * mrSges(4,1) - t880 * mrSges(4,3) + qJD(2) * t889 - t897 * t874 - t802;
t777 = t921 * t783 + t955 * t801;
t870 = Ifges(4,4) * t897 - Ifges(4,2) * t896 + Ifges(4,6) * qJD(2);
t871 = Ifges(4,1) * t897 - Ifges(4,4) * t896 + Ifges(4,5) * qJD(2);
t887 = -t926 * g(3) - t954;
t899 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t924 + Ifges(3,2) * t926) * qJD(1);
t900 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t924 + Ifges(3,4) * t926) * qJD(1);
t965 = mrSges(3,1) * t887 + mrSges(4,1) * t837 - mrSges(3,2) * t888 - mrSges(4,2) * t838 + Ifges(3,5) * t908 + Ifges(4,5) * t880 + Ifges(3,6) * t909 - Ifges(4,6) * t879 + pkin(2) * t777 - pkin(3) * t802 + qJ(4) * t786 + t922 * t770 + t920 * t771 + t897 * t870 + t896 * t871 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t899 * t924 - t900 * t926) * qJD(1);
t804 = t828 * mrSges(7,2) + t855 * t840 - t939;
t964 = t958 * t827 + t959 * t828 + t950 * t854 + t952 * t855 + t966 * t878 - mrSges(6,1) * t809 + mrSges(7,1) * t808 + mrSges(6,2) * t810 - mrSges(7,3) * t807 + pkin(5) * t804 - qJ(6) * (-t827 * mrSges(7,2) - t854 * t840 + t944);
t907 = (-mrSges(3,1) * t926 + mrSges(3,2) * t924) * qJD(1);
t913 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t947;
t775 = m(3) * t887 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t908 + qJD(2) * t913 - t907 * t948 + t777;
t912 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t948;
t941 = t955 * t783 - t801 * t921;
t776 = m(3) * t888 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t909 - qJD(2) * t912 + t907 * t947 + t941;
t942 = -t775 * t924 + t926 * t776;
t766 = m(2) * t915 - mrSges(2,1) * t929 - qJDD(1) * mrSges(2,2) + t942;
t785 = t922 * t790 + t920 * t791;
t784 = m(4) * t862 + t879 * mrSges(4,1) + mrSges(4,2) * t880 + t896 * t889 + t890 * t897 + t785;
t902 = -t929 * pkin(7) + t935;
t931 = -m(3) * t902 + t909 * mrSges(3,1) - mrSges(3,2) * t908 - t912 * t948 + t913 * t947 - t784;
t779 = m(2) * t914 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t929 + t931;
t953 = t925 * t766 + t927 * t779;
t768 = t926 * t775 + t924 * t776;
t943 = t927 * t766 - t779 * t925;
t869 = Ifges(4,5) * t897 - Ifges(4,6) * t896 + Ifges(4,3) * qJD(2);
t763 = mrSges(4,2) * t862 - mrSges(4,3) * t837 + Ifges(4,1) * t880 - Ifges(4,4) * t879 + Ifges(4,5) * qJDD(2) - qJ(4) * t785 - qJD(2) * t870 - t770 * t920 + t771 * t922 - t869 * t896;
t769 = t964 + Ifges(4,6) * qJDD(2) + (-Ifges(5,3) - Ifges(4,2)) * t879 - t897 * t869 + t885 * t850 - t886 * t849 + Ifges(4,4) * t880 - mrSges(4,1) * t862 - Ifges(5,6) * t867 - Ifges(5,5) * t868 + qJD(2) * t871 + mrSges(4,3) * t838 - mrSges(5,1) * t816 + mrSges(5,2) * t817 - pkin(4) * t792 - pkin(3) * t785;
t898 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t924 + Ifges(3,6) * t926) * qJD(1);
t759 = -mrSges(3,1) * t902 + mrSges(3,3) * t888 + Ifges(3,4) * t908 + Ifges(3,2) * t909 + Ifges(3,6) * qJDD(2) - pkin(2) * t784 + qJ(3) * t941 + qJD(2) * t900 + t921 * t763 + t769 * t955 - t898 * t948;
t762 = mrSges(3,2) * t902 - mrSges(3,3) * t887 + Ifges(3,1) * t908 + Ifges(3,4) * t909 + Ifges(3,5) * qJDD(2) - qJ(3) * t777 - qJD(2) * t899 + t763 * t955 - t921 * t769 + t898 * t947;
t934 = mrSges(2,1) * t914 - mrSges(2,2) * t915 + Ifges(2,3) * qJDD(1) + pkin(1) * t931 + pkin(7) * t942 + t926 * t759 + t924 * t762;
t760 = mrSges(2,1) * g(3) + mrSges(2,3) * t915 + t929 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t768 - t965;
t757 = -mrSges(2,2) * g(3) - mrSges(2,3) * t914 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t929 - pkin(7) * t768 - t759 * t924 + t762 * t926;
t1 = [-m(1) * g(1) + t943; -m(1) * g(2) + t953; (-m(1) - m(2)) * g(3) + t768; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t953 + t927 * t757 - t925 * t760; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t943 + t925 * t757 + t927 * t760; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t934; t934; t965; t784; t802; -t964; t804;];
tauJB  = t1;
