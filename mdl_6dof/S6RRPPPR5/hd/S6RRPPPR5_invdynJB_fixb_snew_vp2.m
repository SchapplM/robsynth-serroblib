% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-05-06 08:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:52:12
% EndTime: 2019-05-06 08:52:20
% DurationCPUTime: 5.53s
% Computational Cost: add. (53037->348), mult. (116896->404), div. (0->0), fcn. (71863->8), ass. (0->142)
t963 = sin(qJ(2));
t1004 = qJD(1) * t963;
t1012 = cos(pkin(9));
t961 = sin(pkin(9));
t938 = t961 * qJD(2) + t1004 * t1012;
t1024 = -0.2e1 * t938;
t1023 = 2 * qJD(4);
t1013 = -mrSges(5,2) + mrSges(6,3);
t1022 = Ifges(4,1) + Ifges(6,3) + Ifges(5,2);
t1021 = Ifges(5,1) + Ifges(6,2) + Ifges(4,3);
t1020 = Ifges(6,1) + Ifges(4,2) + Ifges(5,3);
t999 = -Ifges(4,5) + Ifges(6,6) + Ifges(5,4);
t998 = -Ifges(4,6) + Ifges(5,5) - Ifges(6,4);
t997 = -Ifges(5,6) + Ifges(6,5) - Ifges(4,4);
t964 = sin(qJ(1));
t967 = cos(qJ(1));
t950 = -t967 * g(1) - t964 * g(2);
t969 = qJD(1) ^ 2;
t932 = -t969 * pkin(1) + qJDD(1) * pkin(7) + t950;
t966 = cos(qJ(2));
t904 = -t966 * g(3) - t963 * t932;
t942 = (-pkin(2) * t966 - qJ(3) * t963) * qJD(1);
t968 = qJD(2) ^ 2;
t979 = qJDD(2) * pkin(2) + t968 * qJ(3) - t942 * t1004 - qJDD(3) + t904;
t1001 = t966 * qJD(1);
t937 = -qJD(2) * t1012 + t1004 * t961;
t993 = t937 * t1001;
t1019 = -qJ(4) * t993 - t979;
t949 = t964 * g(1) - t967 * g(2);
t931 = -qJDD(1) * pkin(1) - t969 * pkin(7) - t949;
t1000 = qJD(1) * qJD(2);
t990 = t966 * t1000;
t944 = t963 * qJDD(1) + t990;
t954 = t963 * t1000;
t945 = t966 * qJDD(1) - t954;
t873 = (-t944 - t990) * qJ(3) + (-t945 + t954) * pkin(2) + t931;
t905 = -t963 * g(3) + t966 * t932;
t877 = -t968 * pkin(2) + qJDD(2) * qJ(3) + t1001 * t942 + t905;
t1006 = t1012 * t877 + t961 * t873;
t1009 = t966 ^ 2 * t969;
t900 = t937 * pkin(3) - t938 * qJ(4);
t1018 = pkin(3) * t1009 + t945 * qJ(4) + t1001 * t1023 + t937 * t900 - t1006;
t862 = qJD(3) * t1024 + t1012 * t873 - t961 * t877;
t913 = -t938 * mrSges(6,2) + mrSges(6,3) * t1001;
t1010 = t937 * t938;
t1016 = 2 * qJD(5);
t860 = t945 * pkin(3) - qJ(4) * t1009 + t938 * t900 + qJDD(4) - t862;
t919 = t961 * qJDD(2) + t1012 * t944;
t854 = t1001 * t1016 + (t945 + t1010) * qJ(5) + (t919 - t993) * pkin(4) + t860;
t920 = -pkin(5) * t1001 - t937 * pkin(8);
t935 = t938 ^ 2;
t851 = -t935 * pkin(5) + t919 * pkin(8) + t1001 * t920 + t854;
t912 = t938 * pkin(4) + qJ(5) * t1001;
t918 = -qJDD(2) * t1012 + t961 * t944;
t1003 = qJD(3) * t937;
t926 = -0.2e1 * t1003;
t934 = t937 ^ 2;
t856 = -t918 * pkin(4) - t934 * qJ(5) - t1001 * t912 + qJDD(5) - t1018 + t926;
t992 = t938 * t1001;
t852 = t856 + (-t945 + t1010) * pkin(5) + (-t918 - t992) * pkin(8);
t962 = sin(qJ(6));
t965 = cos(qJ(6));
t849 = -t962 * t851 + t965 * t852;
t897 = -t962 * t937 + t965 * t938;
t866 = t897 * qJD(6) + t965 * t918 + t962 * t919;
t898 = t965 * t937 + t962 * t938;
t871 = -t897 * mrSges(7,1) + t898 * mrSges(7,2);
t951 = qJD(6) - t1001;
t878 = -t951 * mrSges(7,2) + t897 * mrSges(7,3);
t941 = qJDD(6) - t945;
t847 = m(7) * t849 + t941 * mrSges(7,1) - t866 * mrSges(7,3) - t898 * t871 + t951 * t878;
t850 = t965 * t851 + t962 * t852;
t865 = -t898 * qJD(6) - t962 * t918 + t965 * t919;
t879 = t951 * mrSges(7,1) - t898 * mrSges(7,3);
t848 = m(7) * t850 - t941 * mrSges(7,2) + t865 * mrSges(7,3) + t897 * t871 - t951 * t879;
t836 = t965 * t847 + t962 * t848;
t899 = t938 * mrSges(6,1) - t937 * mrSges(6,3);
t986 = m(6) * t856 + t918 * mrSges(6,2) + t937 * t899 + t836;
t835 = -t945 * mrSges(6,1) - t1001 * t913 + t986;
t1015 = t934 * pkin(4);
t973 = qJD(4) * t1024 + (t918 - t992) * pkin(3) + t1019;
t853 = -t1015 - t935 * pkin(8) + t973 + (t1016 + t920) * t937 + (-pkin(5) - qJ(4)) * t919 + t918 * qJ(5) - t938 * t912;
t867 = Ifges(7,5) * t898 + Ifges(7,6) * t897 + Ifges(7,3) * t951;
t869 = Ifges(7,1) * t898 + Ifges(7,4) * t897 + Ifges(7,5) * t951;
t840 = -mrSges(7,1) * t853 + mrSges(7,3) * t850 + Ifges(7,4) * t866 + Ifges(7,2) * t865 + Ifges(7,6) * t941 - t898 * t867 + t951 * t869;
t868 = Ifges(7,4) * t898 + Ifges(7,2) * t897 + Ifges(7,6) * t951;
t841 = mrSges(7,2) * t853 - mrSges(7,3) * t849 + Ifges(7,1) * t866 + Ifges(7,4) * t865 + Ifges(7,5) * t941 + t897 * t867 - t951 * t868;
t915 = t938 * mrSges(5,1) - mrSges(5,2) * t1001;
t1011 = t919 * qJ(4);
t861 = t973 - t1011;
t914 = t937 * mrSges(5,1) + mrSges(5,3) * t1001;
t858 = t1015 + t1011 - 0.2e1 * qJD(5) * t937 + (-pkin(3) - qJ(5)) * t918 + (pkin(3) * t1001 + t1023 + t912) * t938 - t1019;
t916 = -mrSges(6,1) * t1001 + t937 * mrSges(6,2);
t980 = m(7) * t853 - t865 * mrSges(7,1) + t866 * mrSges(7,2) - t897 * t878 + t898 * t879;
t977 = -m(6) * t858 - t919 * mrSges(6,1) - t938 * t913 + t937 * t916 + t980;
t972 = m(5) * t861 + t1013 * t918 - t937 * t914 + t977;
t842 = -t919 * mrSges(5,3) - t938 * t915 + t972;
t859 = 0.2e1 * t1003 + t1018;
t863 = t926 + t1006;
t995 = t1001 * t1021 - t937 * t998 + t938 * t999;
t996 = -t1001 * t999 - t1022 * t938 - t937 * t997;
t815 = pkin(4) * t835 + pkin(8) * t836 - pkin(3) * t842 - qJ(5) * t977 - mrSges(6,2) * t856 + t962 * t840 - t965 * t841 + mrSges(6,3) * t858 - mrSges(5,1) * t859 + mrSges(5,2) * t861 + mrSges(4,3) * t863 + mrSges(4,1) * t979 + t998 * t945 + t995 * t938 - t997 * t919 + (-qJ(5) * mrSges(6,3) - t1020) * t918 + t996 * t1001;
t1007 = -t962 * t847 + t965 * t848;
t978 = m(6) * t854 - t919 * mrSges(6,2) + t1001 * t916 - t938 * t899 + t1007;
t834 = t945 * mrSges(6,3) + t978;
t994 = -t1001 * t998 + t1020 * t937 + t938 * t997;
t816 = -qJ(4) * t842 + pkin(4) * t834 - mrSges(6,2) * t854 + t962 * t841 + t965 * t840 + mrSges(6,1) * t858 + mrSges(5,1) * t860 - mrSges(5,3) * t861 - mrSges(4,3) * t862 - mrSges(4,2) * t979 - pkin(5) * t980 + pkin(8) * t1007 + t999 * t945 + t995 * t937 + t1022 * t919 + t997 * t918 - t994 * t1001;
t901 = t937 * mrSges(4,1) + t938 * mrSges(4,2);
t902 = -t937 * mrSges(5,2) - t938 * mrSges(5,3);
t976 = m(5) * t860 + t919 * mrSges(5,1) - t914 * t1001 + t938 * t902 + t978;
t984 = mrSges(4,2) * t1001 - t937 * mrSges(4,3);
t830 = -t984 * t1001 - t976 + (-mrSges(4,1) - t1013) * t945 + m(4) * t862 - t919 * mrSges(4,3) - t938 * t901;
t1005 = mrSges(4,1) * t1001 + t938 * mrSges(4,3) + t915;
t982 = -m(5) * t859 - t945 * mrSges(5,3) + t986;
t832 = m(4) * t863 + (-mrSges(6,1) + mrSges(4,2)) * t945 + (-t901 - t902) * t937 + (-mrSges(4,3) - mrSges(5,1)) * t918 + (-t913 - t1005) * t1001 + t982;
t829 = t1012 * t832 - t961 * t830;
t839 = -m(4) * t979 + t918 * mrSges(4,1) - t1005 * t938 + (mrSges(4,2) - mrSges(5,3)) * t919 + t937 * t984 + t972;
t929 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t963 + Ifges(3,2) * t966) * qJD(1);
t930 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t963 + Ifges(3,4) * t966) * qJD(1);
t1017 = mrSges(3,1) * t904 - mrSges(3,2) * t905 + Ifges(3,5) * t944 + Ifges(3,6) * t945 + Ifges(3,3) * qJDD(2) - pkin(2) * t839 + qJ(3) * t829 + (t929 * t963 - t930 * t966) * qJD(1) + t1012 * t815 + t961 * t816;
t943 = (-mrSges(3,1) * t966 + mrSges(3,2) * t963) * qJD(1);
t947 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1004;
t827 = m(3) * t905 - qJDD(2) * mrSges(3,2) + t945 * mrSges(3,3) - qJD(2) * t947 + t1001 * t943 + t829;
t948 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1001;
t838 = m(3) * t904 + qJDD(2) * mrSges(3,1) - t944 * mrSges(3,3) + qJD(2) * t948 - t1004 * t943 - t839;
t988 = t966 * t827 - t963 * t838;
t819 = m(2) * t950 - t969 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t988;
t828 = t1012 * t830 + t961 * t832;
t975 = -m(3) * t931 + t945 * mrSges(3,1) - t944 * mrSges(3,2) + t948 * t1001 - t1004 * t947 - t828;
t823 = m(2) * t949 + qJDD(1) * mrSges(2,1) - t969 * mrSges(2,2) + t975;
t1008 = t964 * t819 + t967 * t823;
t821 = t963 * t827 + t966 * t838;
t989 = t967 * t819 - t964 * t823;
t928 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t963 + Ifges(3,6) * t966) * qJD(1);
t812 = mrSges(3,2) * t931 - mrSges(3,3) * t904 + Ifges(3,1) * t944 + Ifges(3,4) * t945 + Ifges(3,5) * qJDD(2) - qJ(3) * t828 - qJD(2) * t929 + t1001 * t928 + t1012 * t816 - t961 * t815;
t833 = t1013 * t945 + t976;
t974 = mrSges(7,1) * t849 - mrSges(7,2) * t850 + Ifges(7,5) * t866 + Ifges(7,6) * t865 + Ifges(7,3) * t941 + t898 * t868 - t897 * t869;
t814 = -pkin(2) * t828 - pkin(5) * t836 + qJ(5) * t834 + qJD(2) * t930 - mrSges(3,1) * t931 + pkin(3) * t833 + mrSges(6,3) * t854 - mrSges(6,1) * t856 - t974 + mrSges(5,3) * t859 - mrSges(5,2) * t860 - mrSges(4,1) * t862 + mrSges(4,2) * t863 + Ifges(3,6) * qJDD(2) - qJ(4) * t982 + Ifges(3,4) * t944 + (-t963 * t928 - qJ(4) * (-t913 - t915) * t966) * qJD(1) + mrSges(3,3) * t905 + (qJ(4) * mrSges(5,1) - t998) * t918 + t999 * t919 + (qJ(4) * t902 + t996) * t937 + t994 * t938 + (qJ(4) * mrSges(6,1) + Ifges(3,2) + t1021) * t945;
t981 = mrSges(2,1) * t949 - mrSges(2,2) * t950 + Ifges(2,3) * qJDD(1) + pkin(1) * t975 + pkin(7) * t988 + t963 * t812 + t966 * t814;
t810 = mrSges(2,1) * g(3) + mrSges(2,3) * t950 + t969 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t821 - t1017;
t809 = -mrSges(2,2) * g(3) - mrSges(2,3) * t949 + Ifges(2,5) * qJDD(1) - t969 * Ifges(2,6) - pkin(7) * t821 + t966 * t812 - t963 * t814;
t1 = [-m(1) * g(1) + t989; -m(1) * g(2) + t1008; (-m(1) - m(2)) * g(3) + t821; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1008 + t967 * t809 - t964 * t810; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t989 + t964 * t809 + t967 * t810; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t981; t981; t1017; t839; t833; t835; t974;];
tauJB  = t1;
