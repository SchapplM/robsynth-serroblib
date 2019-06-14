% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-07 05:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:03:10
% EndTime: 2019-05-07 05:04:18
% DurationCPUTime: 64.99s
% Computational Cost: add. (1107670->399), mult. (2468148->519), div. (0->0), fcn. (1997364->14), ass. (0->164)
t1028 = -2 * qJD(4);
t1013 = qJD(1) * qJD(2);
t987 = sin(pkin(6));
t992 = sin(qJ(2));
t996 = cos(qJ(2));
t973 = (-qJDD(1) * t996 + t1013 * t992) * t987;
t1023 = cos(pkin(11));
t1015 = qJD(1) * t996;
t989 = cos(pkin(6));
t1020 = t989 * t992;
t1026 = pkin(8) * t987;
t993 = sin(qJ(1));
t997 = cos(qJ(1));
t977 = t993 * g(1) - g(2) * t997;
t998 = qJD(1) ^ 2;
t968 = qJDD(1) * pkin(1) + t1026 * t998 + t977;
t978 = -g(1) * t997 - g(2) * t993;
t969 = -pkin(1) * t998 + qJDD(1) * t1026 + t978;
t1017 = t968 * t1020 + t996 * t969;
t1016 = qJD(1) * t987;
t971 = (-pkin(2) * t996 - pkin(9) * t992) * t1016;
t981 = qJD(1) * t989 + qJD(2);
t979 = t981 ^ 2;
t980 = qJDD(1) * t989 + qJDD(2);
t924 = -t979 * pkin(2) + t980 * pkin(9) + (-g(3) * t992 + t1015 * t971) * t987 + t1017;
t1025 = t989 * g(3);
t972 = (qJDD(1) * t992 + t1013 * t996) * t987;
t925 = t973 * pkin(2) - t972 * pkin(9) - t1025 + (-t968 + (pkin(2) * t992 - pkin(9) * t996) * t981 * qJD(1)) * t987;
t991 = sin(qJ(3));
t995 = cos(qJ(3));
t895 = -t991 * t924 + t995 * t925;
t1012 = t992 * t1016;
t961 = -t1012 * t991 + t981 * t995;
t941 = qJD(3) * t961 + t972 * t995 + t980 * t991;
t962 = t1012 * t995 + t981 * t991;
t965 = qJDD(3) + t973;
t1011 = t987 * t1015;
t976 = qJD(3) - t1011;
t883 = (t961 * t976 - t941) * qJ(4) + (t961 * t962 + t965) * pkin(3) + t895;
t896 = t995 * t924 + t991 * t925;
t940 = -qJD(3) * t962 - t972 * t991 + t980 * t995;
t951 = pkin(3) * t976 - qJ(4) * t962;
t960 = t961 ^ 2;
t890 = -pkin(3) * t960 + qJ(4) * t940 - t951 * t976 + t896;
t986 = sin(pkin(11));
t948 = t1023 * t962 + t986 * t961;
t874 = t1023 * t883 + t1028 * t948 - t986 * t890;
t985 = sin(pkin(12));
t988 = cos(pkin(12));
t929 = -t948 * t985 + t976 * t988;
t947 = -t1023 * t961 + t962 * t986;
t1006 = -mrSges(6,2) * t947 + mrSges(6,3) * t929;
t875 = t1023 * t890 + t1028 * t947 + t986 * t883;
t918 = pkin(4) * t947 - qJ(5) * t948;
t975 = t976 ^ 2;
t873 = -pkin(4) * t975 + qJ(5) * t965 - t918 * t947 + t875;
t1019 = t989 * t996;
t1021 = t987 * t996;
t942 = -g(3) * t1021 + t1019 * t968 - t992 * t969;
t923 = -t980 * pkin(2) - t979 * pkin(9) + t971 * t1012 - t942;
t892 = -t940 * pkin(3) - t960 * qJ(4) + t962 * t951 + qJDD(4) + t923;
t912 = -t1023 * t940 + t941 * t986;
t913 = t1023 * t941 + t986 * t940;
t878 = (t947 * t976 - t913) * qJ(5) + (t948 * t976 + t912) * pkin(4) + t892;
t930 = t948 * t988 + t976 * t985;
t868 = -0.2e1 * qJD(5) * t930 - t985 * t873 + t988 * t878;
t904 = t913 * t988 + t965 * t985;
t866 = (t929 * t947 - t904) * pkin(10) + (t929 * t930 + t912) * pkin(5) + t868;
t869 = 0.2e1 * qJD(5) * t929 + t988 * t873 + t985 * t878;
t903 = -t913 * t985 + t965 * t988;
t909 = pkin(5) * t947 - pkin(10) * t930;
t928 = t929 ^ 2;
t867 = -pkin(5) * t928 + pkin(10) * t903 - t909 * t947 + t869;
t990 = sin(qJ(6));
t994 = cos(qJ(6));
t864 = t866 * t994 - t867 * t990;
t905 = t929 * t994 - t930 * t990;
t881 = qJD(6) * t905 + t903 * t990 + t904 * t994;
t906 = t929 * t990 + t930 * t994;
t891 = -mrSges(7,1) * t905 + mrSges(7,2) * t906;
t946 = qJD(6) + t947;
t893 = -mrSges(7,2) * t946 + mrSges(7,3) * t905;
t911 = qJDD(6) + t912;
t861 = m(7) * t864 + mrSges(7,1) * t911 - mrSges(7,3) * t881 - t891 * t906 + t893 * t946;
t865 = t866 * t990 + t867 * t994;
t880 = -qJD(6) * t906 + t903 * t994 - t904 * t990;
t894 = mrSges(7,1) * t946 - mrSges(7,3) * t906;
t862 = m(7) * t865 - mrSges(7,2) * t911 + mrSges(7,3) * t880 + t891 * t905 - t894 * t946;
t853 = t994 * t861 + t990 * t862;
t907 = -mrSges(6,1) * t929 + mrSges(6,2) * t930;
t851 = m(6) * t868 + t912 * mrSges(6,1) - t904 * mrSges(6,3) + t1006 * t947 - t930 * t907 + t853;
t1007 = -t861 * t990 + t994 * t862;
t908 = mrSges(6,1) * t947 - mrSges(6,3) * t930;
t852 = m(6) * t869 - mrSges(6,2) * t912 + mrSges(6,3) * t903 + t907 * t929 - t908 * t947 + t1007;
t847 = -t851 * t985 + t988 * t852;
t919 = mrSges(5,1) * t947 + mrSges(5,2) * t948;
t932 = mrSges(5,1) * t976 - mrSges(5,3) * t948;
t844 = m(5) * t875 - mrSges(5,2) * t965 - mrSges(5,3) * t912 - t919 * t947 - t932 * t976 + t847;
t872 = -t965 * pkin(4) - t975 * qJ(5) + t948 * t918 + qJDD(5) - t874;
t870 = -t903 * pkin(5) - t928 * pkin(10) + t930 * t909 + t872;
t1003 = m(7) * t870 - t880 * mrSges(7,1) + mrSges(7,2) * t881 - t905 * t893 + t894 * t906;
t863 = m(6) * t872 - t903 * mrSges(6,1) + mrSges(6,2) * t904 - t929 * t1006 + t908 * t930 + t1003;
t931 = -mrSges(5,2) * t976 - mrSges(5,3) * t947;
t857 = m(5) * t874 + mrSges(5,1) * t965 - mrSges(5,3) * t913 - t919 * t948 + t931 * t976 - t863;
t836 = t1023 * t857 + t986 * t844;
t884 = Ifges(7,5) * t906 + Ifges(7,6) * t905 + Ifges(7,3) * t946;
t886 = Ifges(7,1) * t906 + Ifges(7,4) * t905 + Ifges(7,5) * t946;
t854 = -mrSges(7,1) * t870 + mrSges(7,3) * t865 + Ifges(7,4) * t881 + Ifges(7,2) * t880 + Ifges(7,6) * t911 - t884 * t906 + t886 * t946;
t885 = Ifges(7,4) * t906 + Ifges(7,2) * t905 + Ifges(7,6) * t946;
t855 = mrSges(7,2) * t870 - mrSges(7,3) * t864 + Ifges(7,1) * t881 + Ifges(7,4) * t880 + Ifges(7,5) * t911 + t884 * t905 - t885 * t946;
t897 = Ifges(6,5) * t930 + Ifges(6,6) * t929 + Ifges(6,3) * t947;
t899 = Ifges(6,1) * t930 + Ifges(6,4) * t929 + Ifges(6,5) * t947;
t837 = -mrSges(6,1) * t872 + mrSges(6,3) * t869 + Ifges(6,4) * t904 + Ifges(6,2) * t903 + Ifges(6,6) * t912 - pkin(5) * t1003 + pkin(10) * t1007 + t994 * t854 + t990 * t855 - t930 * t897 + t947 * t899;
t898 = Ifges(6,4) * t930 + Ifges(6,2) * t929 + Ifges(6,6) * t947;
t838 = mrSges(6,2) * t872 - mrSges(6,3) * t868 + Ifges(6,1) * t904 + Ifges(6,4) * t903 + Ifges(6,5) * t912 - pkin(10) * t853 - t854 * t990 + t855 * t994 + t897 * t929 - t898 * t947;
t915 = Ifges(5,4) * t948 - Ifges(5,2) * t947 + Ifges(5,6) * t976;
t916 = Ifges(5,1) * t948 - Ifges(5,4) * t947 + Ifges(5,5) * t976;
t935 = Ifges(4,4) * t962 + Ifges(4,2) * t961 + Ifges(4,6) * t976;
t936 = Ifges(4,1) * t962 + Ifges(4,4) * t961 + Ifges(4,5) * t976;
t1027 = Ifges(4,5) * t941 + Ifges(4,6) * t940 + t962 * t935 - t961 * t936 + mrSges(4,1) * t895 - mrSges(4,2) * t896 + Ifges(5,5) * t913 - Ifges(5,6) * t912 + t948 * t915 + t947 * t916 + mrSges(5,1) * t874 - mrSges(5,2) * t875 + t985 * t838 + t988 * t837 - pkin(4) * t863 + qJ(5) * t847 + pkin(3) * t836 + (Ifges(4,3) + Ifges(5,3)) * t965;
t1022 = t987 * t992;
t949 = -mrSges(4,1) * t961 + mrSges(4,2) * t962;
t950 = -mrSges(4,2) * t976 + mrSges(4,3) * t961;
t834 = m(4) * t895 + mrSges(4,1) * t965 - mrSges(4,3) * t941 - t949 * t962 + t950 * t976 + t836;
t1008 = t1023 * t844 - t857 * t986;
t952 = mrSges(4,1) * t976 - mrSges(4,3) * t962;
t835 = m(4) * t896 - mrSges(4,2) * t965 + mrSges(4,3) * t940 + t949 * t961 - t952 * t976 + t1008;
t1009 = -t834 * t991 + t995 * t835;
t943 = -g(3) * t1022 + t1017;
t966 = mrSges(3,1) * t981 - mrSges(3,3) * t1012;
t970 = (-mrSges(3,1) * t996 + mrSges(3,2) * t992) * t1016;
t826 = m(3) * t943 - mrSges(3,2) * t980 - mrSges(3,3) * t973 + t1011 * t970 - t966 * t981 + t1009;
t829 = t995 * t834 + t991 * t835;
t956 = -t987 * t968 - t1025;
t967 = -mrSges(3,2) * t981 + mrSges(3,3) * t1011;
t828 = m(3) * t956 + t973 * mrSges(3,1) + t972 * mrSges(3,2) + (t966 * t992 - t967 * t996) * t1016 + t829;
t846 = t988 * t851 + t985 * t852;
t845 = m(5) * t892 + t912 * mrSges(5,1) + mrSges(5,2) * t913 + t947 * t931 + t932 * t948 + t846;
t1000 = -m(4) * t923 + t940 * mrSges(4,1) - mrSges(4,2) * t941 + t961 * t950 - t952 * t962 - t845;
t841 = m(3) * t942 + mrSges(3,1) * t980 - mrSges(3,3) * t972 - t1012 * t970 + t967 * t981 + t1000;
t816 = t841 * t1019 + t826 * t1020 - t828 * t987;
t813 = m(2) * t977 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t998 + t816;
t821 = t996 * t826 - t841 * t992;
t819 = m(2) * t978 - mrSges(2,1) * t998 - qJDD(1) * mrSges(2,2) + t821;
t1018 = t997 * t813 + t993 * t819;
t815 = t841 * t1021 + t826 * t1022 + t989 * t828;
t1010 = -t813 * t993 + t997 * t819;
t914 = Ifges(5,5) * t948 - Ifges(5,6) * t947 + Ifges(5,3) * t976;
t822 = mrSges(5,2) * t892 - mrSges(5,3) * t874 + Ifges(5,1) * t913 - Ifges(5,4) * t912 + Ifges(5,5) * t965 - qJ(5) * t846 - t837 * t985 + t838 * t988 - t914 * t947 - t915 * t976;
t1001 = mrSges(7,1) * t864 - mrSges(7,2) * t865 + Ifges(7,5) * t881 + Ifges(7,6) * t880 + Ifges(7,3) * t911 + t906 * t885 - t905 * t886;
t830 = -t1001 + (-Ifges(5,2) - Ifges(6,3)) * t912 + t976 * t916 + Ifges(5,6) * t965 - t948 * t914 + t929 * t899 - t930 * t898 + Ifges(5,4) * t913 - pkin(4) * t846 - Ifges(6,6) * t903 - Ifges(6,5) * t904 - mrSges(5,1) * t892 + mrSges(5,3) * t875 - mrSges(6,1) * t868 + mrSges(6,2) * t869 - pkin(5) * t853;
t934 = Ifges(4,5) * t962 + Ifges(4,6) * t961 + Ifges(4,3) * t976;
t808 = -mrSges(4,1) * t923 + mrSges(4,3) * t896 + Ifges(4,4) * t941 + Ifges(4,2) * t940 + Ifges(4,6) * t965 - pkin(3) * t845 + qJ(4) * t1008 + t1023 * t830 + t986 * t822 - t962 * t934 + t976 * t936;
t811 = mrSges(4,2) * t923 - mrSges(4,3) * t895 + Ifges(4,1) * t941 + Ifges(4,4) * t940 + Ifges(4,5) * t965 - qJ(4) * t836 + t1023 * t822 - t986 * t830 + t961 * t934 - t976 * t935;
t954 = Ifges(3,6) * t981 + (Ifges(3,4) * t992 + Ifges(3,2) * t996) * t1016;
t955 = Ifges(3,5) * t981 + (Ifges(3,1) * t992 + Ifges(3,4) * t996) * t1016;
t805 = Ifges(3,5) * t972 - Ifges(3,6) * t973 + Ifges(3,3) * t980 + mrSges(3,1) * t942 - mrSges(3,2) * t943 + t991 * t811 + t995 * t808 + pkin(2) * t1000 + pkin(9) * t1009 + (t954 * t992 - t955 * t996) * t1016;
t953 = Ifges(3,3) * t981 + (Ifges(3,5) * t992 + Ifges(3,6) * t996) * t1016;
t807 = mrSges(3,2) * t956 - mrSges(3,3) * t942 + Ifges(3,1) * t972 - Ifges(3,4) * t973 + Ifges(3,5) * t980 - pkin(9) * t829 + t1011 * t953 - t808 * t991 + t811 * t995 - t954 * t981;
t810 = -mrSges(3,1) * t956 + mrSges(3,3) * t943 + Ifges(3,4) * t972 - Ifges(3,2) * t973 + Ifges(3,6) * t980 - pkin(2) * t829 - t1012 * t953 + t981 * t955 - t1027;
t1002 = mrSges(2,1) * t977 - mrSges(2,2) * t978 + Ifges(2,3) * qJDD(1) + pkin(1) * t816 + t810 * t1021 + t807 * t1022 + t1026 * t821 + t989 * t805;
t803 = -mrSges(2,2) * g(3) - mrSges(2,3) * t977 + Ifges(2,5) * qJDD(1) - t998 * Ifges(2,6) + t996 * t807 - t992 * t810 + (-t815 * t987 - t816 * t989) * pkin(8);
t802 = mrSges(2,1) * g(3) + mrSges(2,3) * t978 + t998 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t815 - t987 * t805 + (pkin(8) * t821 + t807 * t992 + t810 * t996) * t989;
t1 = [-m(1) * g(1) + t1010; -m(1) * g(2) + t1018; (-m(1) - m(2)) * g(3) + t815; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1018 - t993 * t802 + t997 * t803; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1010 + t997 * t802 + t993 * t803; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1002; t1002; t805; t1027; t845; t863; t1001;];
tauJB  = t1;
