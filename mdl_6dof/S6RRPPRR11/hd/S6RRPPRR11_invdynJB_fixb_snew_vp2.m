% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 12:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR11_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:00:17
% EndTime: 2019-05-06 12:00:40
% DurationCPUTime: 21.03s
% Computational Cost: add. (325866->384), mult. (769250->481), div. (0->0), fcn. (570437->12), ass. (0->162)
t1070 = -2 * qJD(3);
t1069 = Ifges(3,1) + Ifges(4,2);
t1064 = Ifges(3,4) + Ifges(4,6);
t1063 = Ifges(3,5) - Ifges(4,4);
t1068 = Ifges(3,2) + Ifges(4,3);
t1062 = Ifges(3,6) - Ifges(4,5);
t1067 = Ifges(3,3) + Ifges(4,1);
t1018 = cos(pkin(6));
t1009 = qJD(1) * t1018 + qJD(2);
t1016 = sin(pkin(6));
t1021 = sin(qJ(2));
t1050 = t1016 * t1021;
t1044 = qJD(1) * t1050;
t1066 = (pkin(2) * t1009 + t1070) * t1044;
t1007 = t1009 ^ 2;
t1008 = qJDD(1) * t1018 + qJDD(2);
t1025 = cos(qJ(2));
t1049 = t1016 * t1025;
t1043 = qJD(1) * t1049;
t1048 = t1018 * t1021;
t1022 = sin(qJ(1));
t1026 = cos(qJ(1));
t1004 = t1022 * g(1) - g(2) * t1026;
t1027 = qJD(1) ^ 2;
t1060 = pkin(8) * t1016;
t989 = qJDD(1) * pkin(1) + t1027 * t1060 + t1004;
t1005 = -g(1) * t1026 - g(2) * t1022;
t1046 = qJDD(1) * t1016;
t990 = -pkin(1) * t1027 + pkin(8) * t1046 + t1005;
t948 = -g(3) * t1050 + t1025 * t990 + t989 * t1048;
t1053 = qJD(1) * t1016;
t991 = (-pkin(2) * t1025 - qJ(3) * t1021) * t1053;
t935 = pkin(2) * t1007 - t1008 * qJ(3) + t1009 * t1070 - t991 * t1043 - t948;
t1065 = mrSges(3,1) - mrSges(4,2);
t1061 = -pkin(2) - qJ(4);
t1059 = g(3) * t1018;
t1047 = t1018 * t1025;
t1015 = sin(pkin(11));
t1017 = cos(pkin(11));
t1020 = sin(qJ(5));
t1024 = cos(qJ(5));
t1000 = qJD(5) + t1044;
t1019 = sin(qJ(6));
t1023 = cos(qJ(6));
t1051 = t1016 ^ 2 * t1027;
t1042 = t1025 ^ 2 * t1051;
t986 = pkin(3) * t1044 - qJ(4) * t1009;
t994 = (qJD(1) * qJD(2) * t1025 + qJDD(1) * t1021) * t1016;
t995 = -qJD(2) * t1044 + t1025 * t1046;
t915 = -pkin(3) * t1042 - t1059 - qJ(3) * t994 + t1061 * t995 + (-t989 + (-qJ(3) * t1009 * t1025 - t1021 * t986) * qJD(1)) * t1016 + t1066;
t1054 = g(3) * t1049 + t1021 * t990;
t1036 = -qJ(3) * t1007 + t991 * t1044 + qJDD(3) + t1054;
t918 = pkin(3) * t994 + t1061 * t1008 + (-pkin(3) * t1009 * t1053 - qJ(4) * t1021 * t1051 - t1018 * t989) * t1025 + t1036;
t977 = t1009 * t1017 - t1015 * t1043;
t900 = -0.2e1 * qJD(4) * t977 - t1015 * t915 + t1017 * t918;
t957 = t1008 * t1017 - t1015 * t995;
t976 = -t1009 * t1015 - t1017 * t1043;
t897 = (t1044 * t976 - t957) * pkin(9) + (t976 * t977 + t994) * pkin(4) + t900;
t901 = 0.2e1 * qJD(4) * t976 + t1015 * t918 + t1017 * t915;
t956 = -t1008 * t1015 - t1017 * t995;
t958 = pkin(4) * t1044 - pkin(9) * t977;
t975 = t976 ^ 2;
t899 = -pkin(4) * t975 + pkin(9) * t956 - t1044 * t958 + t901;
t894 = t1020 * t897 + t1024 * t899;
t950 = -t1020 * t977 + t1024 * t976;
t951 = t1020 * t976 + t1024 * t977;
t934 = -pkin(5) * t950 - pkin(10) * t951;
t983 = qJDD(5) + t994;
t998 = t1000 ^ 2;
t891 = -pkin(5) * t998 + pkin(10) * t983 + t934 * t950 + t894;
t914 = pkin(3) * t995 - qJ(4) * t1042 + t1009 * t986 + qJDD(4) - t935;
t903 = -pkin(4) * t956 - pkin(9) * t975 + t977 * t958 + t914;
t924 = -qJD(5) * t951 - t1020 * t957 + t1024 * t956;
t925 = qJD(5) * t950 + t1020 * t956 + t1024 * t957;
t895 = (-t1000 * t950 - t925) * pkin(10) + t903 + (t1000 * t951 - t924) * pkin(5);
t888 = -t1019 * t891 + t1023 * t895;
t938 = t1000 * t1023 - t1019 * t951;
t906 = qJD(6) * t938 + t1019 * t983 + t1023 * t925;
t939 = t1000 * t1019 + t1023 * t951;
t919 = -mrSges(7,1) * t938 + mrSges(7,2) * t939;
t923 = qJDD(6) - t924;
t949 = qJD(6) - t950;
t926 = -mrSges(7,2) * t949 + mrSges(7,3) * t938;
t885 = m(7) * t888 + mrSges(7,1) * t923 - mrSges(7,3) * t906 - t919 * t939 + t926 * t949;
t889 = t1019 * t895 + t1023 * t891;
t905 = -qJD(6) * t939 - t1019 * t925 + t1023 * t983;
t927 = mrSges(7,1) * t949 - mrSges(7,3) * t939;
t886 = m(7) * t889 - mrSges(7,2) * t923 + mrSges(7,3) * t905 + t919 * t938 - t927 * t949;
t1040 = -t1019 * t885 + t1023 * t886;
t933 = -mrSges(6,1) * t950 + mrSges(6,2) * t951;
t941 = mrSges(6,1) * t1000 - mrSges(6,3) * t951;
t872 = m(6) * t894 - mrSges(6,2) * t983 + mrSges(6,3) * t924 - t1000 * t941 + t933 * t950 + t1040;
t893 = -t1020 * t899 + t1024 * t897;
t890 = -pkin(5) * t983 - pkin(10) * t998 + t934 * t951 - t893;
t1034 = -m(7) * t890 + t905 * mrSges(7,1) - mrSges(7,2) * t906 + t938 * t926 - t927 * t939;
t940 = -mrSges(6,2) * t1000 + mrSges(6,3) * t950;
t881 = m(6) * t893 + mrSges(6,1) * t983 - mrSges(6,3) * t925 + t1000 * t940 - t933 * t951 + t1034;
t865 = t1020 * t872 + t1024 * t881;
t952 = -mrSges(5,1) * t976 + mrSges(5,2) * t977;
t954 = -mrSges(5,2) * t1044 + mrSges(5,3) * t976;
t863 = m(5) * t900 + mrSges(5,1) * t994 - mrSges(5,3) * t957 + t1044 * t954 - t952 * t977 + t865;
t1039 = -t1020 * t881 + t1024 * t872;
t955 = mrSges(5,1) * t1044 - mrSges(5,3) * t977;
t864 = m(5) * t901 - mrSges(5,2) * t994 + mrSges(5,3) * t956 - t1044 * t955 + t952 * t976 + t1039;
t1041 = -t1015 * t863 + t1017 * t864;
t965 = -t1016 * t989 - t1059;
t936 = -pkin(2) * t995 + (-t1009 * t1043 - t994) * qJ(3) + t965 + t1066;
t987 = -mrSges(4,1) * t1043 - mrSges(4,3) * t1009;
t1037 = m(4) * t936 - t994 * mrSges(4,3) + t987 * t1043 + t1041;
t984 = mrSges(3,1) * t1009 - mrSges(3,3) * t1044;
t985 = -mrSges(3,2) * t1009 + mrSges(3,3) * t1043;
t988 = mrSges(4,1) * t1044 + mrSges(4,2) * t1009;
t855 = m(3) * t965 + mrSges(3,2) * t994 - t1065 * t995 + (-t1025 * t985 + (t984 - t988) * t1021) * t1053 + t1037;
t859 = t1015 * t864 + t1017 * t863;
t1045 = t989 * t1047;
t937 = -pkin(2) * t1008 + t1036 - t1045;
t1035 = -m(4) * t937 - t994 * mrSges(4,1) - t859;
t947 = t1045 - t1054;
t992 = (mrSges(4,2) * t1025 - mrSges(4,3) * t1021) * t1053;
t993 = (-mrSges(3,1) * t1025 + mrSges(3,2) * t1021) * t1053;
t856 = m(3) * t947 - mrSges(3,3) * t994 + (t985 - t987) * t1009 + t1065 * t1008 + (-t992 - t993) * t1044 + t1035;
t875 = t1019 * t886 + t1023 * t885;
t1032 = m(6) * t903 - t924 * mrSges(6,1) + t925 * mrSges(6,2) - t950 * t940 + t951 * t941 + t875;
t873 = m(5) * t914 - t956 * mrSges(5,1) + t957 * mrSges(5,2) - t976 * t954 + t977 * t955 + t1032;
t1028 = -m(4) * t935 + t1008 * mrSges(4,3) + t1009 * t988 + t992 * t1043 + t873;
t869 = t993 * t1043 + (mrSges(3,3) + mrSges(4,1)) * t995 - t1009 * t984 - t1008 * mrSges(3,2) + m(3) * t948 + t1028;
t844 = -t1016 * t855 + t856 * t1047 + t869 * t1048;
t841 = m(2) * t1004 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1027 + t844;
t851 = -t1021 * t856 + t1025 * t869;
t849 = m(2) * t1005 - mrSges(2,1) * t1027 - qJDD(1) * mrSges(2,2) + t851;
t1058 = t1022 * t849 + t1026 * t841;
t1057 = (t1063 * t1021 + t1062 * t1025) * t1053 + t1067 * t1009;
t1056 = (t1064 * t1021 + t1068 * t1025) * t1053 + t1062 * t1009;
t1055 = (t1069 * t1021 + t1064 * t1025) * t1053 + t1063 * t1009;
t843 = t1018 * t855 + t856 * t1049 + t869 * t1050;
t1038 = -t1022 * t841 + t1026 * t849;
t907 = Ifges(7,5) * t939 + Ifges(7,6) * t938 + Ifges(7,3) * t949;
t909 = Ifges(7,1) * t939 + Ifges(7,4) * t938 + Ifges(7,5) * t949;
t878 = -mrSges(7,1) * t890 + mrSges(7,3) * t889 + Ifges(7,4) * t906 + Ifges(7,2) * t905 + Ifges(7,6) * t923 - t907 * t939 + t909 * t949;
t908 = Ifges(7,4) * t939 + Ifges(7,2) * t938 + Ifges(7,6) * t949;
t879 = mrSges(7,2) * t890 - mrSges(7,3) * t888 + Ifges(7,1) * t906 + Ifges(7,4) * t905 + Ifges(7,5) * t923 + t907 * t938 - t908 * t949;
t928 = Ifges(6,5) * t951 + Ifges(6,6) * t950 + Ifges(6,3) * t1000;
t929 = Ifges(6,4) * t951 + Ifges(6,2) * t950 + Ifges(6,6) * t1000;
t860 = mrSges(6,2) * t903 - mrSges(6,3) * t893 + Ifges(6,1) * t925 + Ifges(6,4) * t924 + Ifges(6,5) * t983 - pkin(10) * t875 - t1000 * t929 - t1019 * t878 + t1023 * t879 + t928 * t950;
t1029 = mrSges(7,1) * t888 - mrSges(7,2) * t889 + Ifges(7,5) * t906 + Ifges(7,6) * t905 + Ifges(7,3) * t923 + t908 * t939 - t909 * t938;
t930 = Ifges(6,1) * t951 + Ifges(6,4) * t950 + Ifges(6,5) * t1000;
t861 = -mrSges(6,1) * t903 + mrSges(6,3) * t894 + Ifges(6,4) * t925 + Ifges(6,2) * t924 + Ifges(6,6) * t983 - pkin(5) * t875 + t1000 * t930 - t928 * t951 - t1029;
t942 = Ifges(5,5) * t977 + Ifges(5,6) * t976 + Ifges(5,3) * t1044;
t944 = Ifges(5,1) * t977 + Ifges(5,4) * t976 + Ifges(5,5) * t1044;
t845 = -mrSges(5,1) * t914 + mrSges(5,3) * t901 + Ifges(5,4) * t957 + Ifges(5,2) * t956 + Ifges(5,6) * t994 - pkin(4) * t1032 + pkin(9) * t1039 + t1020 * t860 + t1024 * t861 + t1044 * t944 - t977 * t942;
t943 = Ifges(5,4) * t977 + Ifges(5,2) * t976 + Ifges(5,6) * t1044;
t846 = mrSges(5,2) * t914 - mrSges(5,3) * t900 + Ifges(5,1) * t957 + Ifges(5,4) * t956 + Ifges(5,5) * t994 - pkin(9) * t865 - t1020 * t861 + t1024 * t860 - t1044 * t943 + t942 * t976;
t858 = mrSges(4,2) * t1008 + t1009 * t987 + t1044 * t992 - t1035;
t835 = mrSges(3,1) * t947 - mrSges(3,2) * t948 + mrSges(4,2) * t937 - mrSges(4,3) * t935 + t1017 * t846 - t1015 * t845 - qJ(4) * t859 - pkin(2) * t858 + qJ(3) * t1028 + (qJ(3) * mrSges(4,1) + t1062) * t995 + t1063 * t994 + t1067 * t1008 + (t1021 * t1056 - t1025 * t1055) * t1053;
t857 = mrSges(4,2) * t995 - t1044 * t988 + t1037;
t837 = -mrSges(3,1) * t965 - mrSges(4,1) * t935 + mrSges(4,2) * t936 + mrSges(3,3) * t948 - pkin(2) * t857 + pkin(3) * t873 - qJ(4) * t1041 + t1062 * t1008 + t1055 * t1009 - t1015 * t846 - t1017 * t845 - t1057 * t1044 + t1064 * t994 + t1068 * t995;
t1030 = mrSges(6,1) * t893 - mrSges(6,2) * t894 + Ifges(6,5) * t925 + Ifges(6,6) * t924 + Ifges(6,3) * t983 + pkin(5) * t1034 + pkin(10) * t1040 + t1019 * t879 + t1023 * t878 + t951 * t929 - t950 * t930;
t839 = pkin(3) * t859 + t1064 * t995 + (Ifges(5,3) + t1069) * t994 + t1030 - t1056 * t1009 + t1063 * t1008 - qJ(3) * t857 - t976 * t944 + t977 * t943 + Ifges(5,6) * t956 + Ifges(5,5) * t957 + mrSges(3,2) * t965 - mrSges(3,3) * t947 - mrSges(4,3) * t936 + mrSges(4,1) * t937 + mrSges(5,1) * t900 - mrSges(5,2) * t901 + t1057 * t1043 + pkin(4) * t865;
t1033 = mrSges(2,1) * t1004 - mrSges(2,2) * t1005 + Ifges(2,3) * qJDD(1) + pkin(1) * t844 + t1018 * t835 + t837 * t1049 + t839 * t1050 + t851 * t1060;
t833 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1004 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t1027 - t1021 * t837 + t1025 * t839 + (-t1016 * t843 - t1018 * t844) * pkin(8);
t832 = mrSges(2,1) * g(3) + mrSges(2,3) * t1005 + Ifges(2,5) * t1027 + Ifges(2,6) * qJDD(1) - pkin(1) * t843 - t1016 * t835 + (pkin(8) * t851 + t1021 * t839 + t1025 * t837) * t1018;
t1 = [-m(1) * g(1) + t1038; -m(1) * g(2) + t1058; (-m(1) - m(2)) * g(3) + t843; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1058 - t1022 * t832 + t1026 * t833; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1038 + t1022 * t833 + t1026 * t832; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1033; t1033; t835; t858; t873; t1030; t1029;];
tauJB  = t1;
