% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 13:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 13:03:01
% EndTime: 2019-05-07 13:06:04
% DurationCPUTime: 169.30s
% Computational Cost: add. (2637238->418), mult. (6918048->556), div. (0->0), fcn. (5874879->16), ass. (0->181)
t1024 = sin(pkin(13));
t1027 = cos(pkin(13));
t1029 = cos(pkin(6));
t1021 = qJD(1) * t1029 + qJD(2);
t1025 = sin(pkin(7));
t1028 = cos(pkin(7));
t1026 = sin(pkin(6));
t1038 = cos(qJ(2));
t1060 = t1026 * t1038;
t1052 = qJD(1) * t1060;
t1006 = t1021 * t1028 - t1025 * t1052 + qJD(3);
t1032 = sin(qJ(3));
t1037 = cos(qJ(3));
t1058 = t1028 * t1037;
t1062 = t1025 * t1037;
t1005 = (t1021 * t1025 + t1028 * t1052) * pkin(10);
t1033 = sin(qJ(2));
t1066 = qJD(1) * t1026;
t1071 = pkin(10) * t1025;
t1009 = (-pkin(2) * t1038 - t1033 * t1071) * t1066;
t1064 = qJD(1) * t1038;
t1015 = (qJD(2) * t1064 + qJDD(1) * t1033) * t1026;
t1020 = qJDD(1) * t1029 + qJDD(2);
t1034 = sin(qJ(1));
t1039 = cos(qJ(1));
t1018 = t1034 * g(1) - g(2) * t1039;
t1040 = qJD(1) ^ 2;
t1072 = pkin(9) * t1026;
t1012 = qJDD(1) * pkin(1) + t1040 * t1072 + t1018;
t1019 = -g(1) * t1039 - g(2) * t1034;
t1013 = -pkin(1) * t1040 + qJDD(1) * t1072 + t1019;
t1055 = t1029 * t1038;
t1047 = t1012 * t1055 - t1033 * t1013;
t1065 = qJD(1) * t1033;
t1070 = pkin(10) * t1028;
t965 = -t1015 * t1070 + t1020 * pkin(2) + t1021 * t1005 + (-g(3) * t1038 - t1009 * t1065) * t1026 + t1047;
t1061 = t1026 * t1033;
t1053 = qJD(1) * t1061;
t1008 = pkin(2) * t1021 - t1053 * t1070;
t1016 = (-qJD(2) * t1065 + qJDD(1) * t1038) * t1026;
t1046 = t1016 * t1028 + t1020 * t1025;
t1056 = t1029 * t1033;
t1054 = t1012 * t1056 + t1038 * t1013;
t966 = -t1021 * t1008 + (-g(3) * t1033 + t1009 * t1064) * t1026 + t1046 * pkin(10) + t1054;
t1069 = t1029 * g(3);
t970 = -t1015 * t1071 - t1016 * pkin(2) - t1069 + (-t1012 + (-t1005 * t1038 + t1008 * t1033) * qJD(1)) * t1026;
t931 = -t1032 * t966 + t965 * t1058 + t970 * t1062;
t1057 = t1028 * t1038;
t995 = t1021 * t1062 + (-t1032 * t1033 + t1037 * t1057) * t1066;
t981 = t995 * qJD(3) + t1037 * t1015 + t1032 * t1046;
t1063 = t1025 * t1032;
t996 = t1021 * t1063 + (t1032 * t1057 + t1033 * t1037) * t1066;
t997 = -t1016 * t1025 + t1020 * t1028 + qJDD(3);
t921 = (t1006 * t995 - t981) * qJ(4) + (t995 * t996 + t997) * pkin(3) + t931;
t1059 = t1028 * t1032;
t932 = t1037 * t966 + t965 * t1059 + t970 * t1063;
t980 = -t996 * qJD(3) - t1032 * t1015 + t1037 * t1046;
t990 = pkin(3) * t1006 - qJ(4) * t996;
t994 = t995 ^ 2;
t924 = -pkin(3) * t994 + qJ(4) * t980 - t1006 * t990 + t932;
t987 = t1024 * t995 + t1027 * t996;
t913 = -0.2e1 * qJD(4) * t987 - t1024 * t924 + t1027 * t921;
t1011 = -mrSges(3,2) * t1021 + mrSges(3,3) * t1052;
t1014 = (-mrSges(3,1) * t1038 + mrSges(3,2) * t1033) * t1066;
t1031 = sin(qJ(5));
t1036 = cos(qJ(5));
t1030 = sin(qJ(6));
t1035 = cos(qJ(6));
t1004 = t1006 ^ 2;
t986 = -t1024 * t996 + t1027 * t995;
t914 = 0.2e1 * qJD(4) * t986 + t1024 * t921 + t1027 * t924;
t961 = -pkin(4) * t986 - pkin(11) * t987;
t912 = -pkin(4) * t1004 + pkin(11) * t997 + t961 * t986 + t914;
t944 = -t1025 * t965 + t1028 * t970;
t930 = -t980 * pkin(3) - t994 * qJ(4) + t996 * t990 + qJDD(4) + t944;
t952 = -t1024 * t981 + t1027 * t980;
t953 = t1024 * t980 + t1027 * t981;
t916 = (-t1006 * t986 - t953) * pkin(11) + (t1006 * t987 - t952) * pkin(4) + t930;
t908 = t1031 * t916 + t1036 * t912;
t972 = t1006 * t1036 - t1031 * t987;
t973 = t1006 * t1031 + t1036 * t987;
t947 = -pkin(5) * t972 - pkin(12) * t973;
t949 = qJDD(5) - t952;
t985 = qJD(5) - t986;
t984 = t985 ^ 2;
t906 = -pkin(5) * t984 + pkin(12) * t949 + t947 * t972 + t908;
t911 = -t997 * pkin(4) - t1004 * pkin(11) + t987 * t961 - t913;
t936 = -qJD(5) * t973 - t1031 * t953 + t1036 * t997;
t937 = qJD(5) * t972 + t1031 * t997 + t1036 * t953;
t909 = (-t972 * t985 - t937) * pkin(12) + (t973 * t985 - t936) * pkin(5) + t911;
t902 = -t1030 * t906 + t1035 * t909;
t950 = -t1030 * t973 + t1035 * t985;
t919 = qJD(6) * t950 + t1030 * t949 + t1035 * t937;
t951 = t1030 * t985 + t1035 * t973;
t933 = -mrSges(7,1) * t950 + mrSges(7,2) * t951;
t935 = qJDD(6) - t936;
t971 = qJD(6) - t972;
t938 = -mrSges(7,2) * t971 + mrSges(7,3) * t950;
t900 = m(7) * t902 + mrSges(7,1) * t935 - mrSges(7,3) * t919 - t933 * t951 + t938 * t971;
t903 = t1030 * t909 + t1035 * t906;
t918 = -qJD(6) * t951 - t1030 * t937 + t1035 * t949;
t939 = mrSges(7,1) * t971 - mrSges(7,3) * t951;
t901 = m(7) * t903 - mrSges(7,2) * t935 + mrSges(7,3) * t918 + t933 * t950 - t939 * t971;
t894 = -t1030 * t900 + t1035 * t901;
t946 = -mrSges(6,1) * t972 + mrSges(6,2) * t973;
t955 = mrSges(6,1) * t985 - mrSges(6,3) * t973;
t892 = m(6) * t908 - mrSges(6,2) * t949 + mrSges(6,3) * t936 + t946 * t972 - t955 * t985 + t894;
t907 = -t1031 * t912 + t1036 * t916;
t905 = -pkin(5) * t949 - pkin(12) * t984 + t947 * t973 - t907;
t904 = -m(7) * t905 + t918 * mrSges(7,1) - mrSges(7,2) * t919 + t950 * t938 - t939 * t951;
t954 = -mrSges(6,2) * t985 + mrSges(6,3) * t972;
t898 = m(6) * t907 + mrSges(6,1) * t949 - mrSges(6,3) * t937 - t946 * t973 + t954 * t985 + t904;
t1049 = -t1031 * t898 + t1036 * t892;
t960 = -mrSges(5,1) * t986 + mrSges(5,2) * t987;
t975 = mrSges(5,1) * t1006 - mrSges(5,3) * t987;
t883 = m(5) * t914 - mrSges(5,2) * t997 + mrSges(5,3) * t952 - t1006 * t975 + t960 * t986 + t1049;
t893 = t1030 * t901 + t1035 * t900;
t1043 = -m(6) * t911 + t936 * mrSges(6,1) - mrSges(6,2) * t937 + t972 * t954 - t955 * t973 - t893;
t974 = -mrSges(5,2) * t1006 + mrSges(5,3) * t986;
t889 = m(5) * t913 + mrSges(5,1) * t997 - mrSges(5,3) * t953 + t1006 * t974 - t960 * t987 + t1043;
t878 = t1024 * t883 + t1027 * t889;
t988 = -mrSges(4,1) * t995 + mrSges(4,2) * t996;
t989 = -mrSges(4,2) * t1006 + mrSges(4,3) * t995;
t876 = m(4) * t931 + mrSges(4,1) * t997 - mrSges(4,3) * t981 + t1006 * t989 - t988 * t996 + t878;
t1051 = -t1024 * t889 + t1027 * t883;
t991 = mrSges(4,1) * t1006 - mrSges(4,3) * t996;
t877 = m(4) * t932 - mrSges(4,2) * t997 + mrSges(4,3) * t980 - t1006 * t991 + t988 * t995 + t1051;
t887 = t1031 * t892 + t1036 * t898;
t886 = m(5) * t930 - mrSges(5,1) * t952 + t953 * mrSges(5,2) - t974 * t986 + t987 * t975 + t887;
t885 = m(4) * t944 - mrSges(4,1) * t980 + mrSges(4,2) * t981 - t989 * t995 + t991 * t996 + t886;
t863 = -t1025 * t885 + t876 * t1058 + t877 * t1059;
t992 = -g(3) * t1060 + t1047;
t859 = m(3) * t992 + mrSges(3,1) * t1020 - mrSges(3,3) * t1015 + t1011 * t1021 - t1014 * t1053 + t863;
t1001 = -t1026 * t1012 - t1069;
t1010 = mrSges(3,1) * t1021 - mrSges(3,3) * t1053;
t862 = t1028 * t885 + t876 * t1062 + t877 * t1063;
t861 = m(3) * t1001 - t1016 * mrSges(3,1) + t1015 * mrSges(3,2) + (t1010 * t1033 - t1011 * t1038) * t1066 + t862;
t869 = -t1032 * t876 + t1037 * t877;
t993 = -g(3) * t1061 + t1054;
t868 = m(3) * t993 - mrSges(3,2) * t1020 + mrSges(3,3) * t1016 - t1010 * t1021 + t1014 * t1052 + t869;
t848 = -t1026 * t861 + t859 * t1055 + t868 * t1056;
t845 = m(2) * t1018 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1040 + t848;
t855 = -t1033 * t859 + t1038 * t868;
t853 = m(2) * t1019 - mrSges(2,1) * t1040 - qJDD(1) * mrSges(2,2) + t855;
t1068 = t1034 * t853 + t1039 * t845;
t847 = t1029 * t861 + t859 * t1060 + t868 * t1061;
t1048 = -t1034 * t845 + t1039 * t853;
t925 = Ifges(7,5) * t951 + Ifges(7,6) * t950 + Ifges(7,3) * t971;
t927 = Ifges(7,1) * t951 + Ifges(7,4) * t950 + Ifges(7,5) * t971;
t895 = -mrSges(7,1) * t905 + mrSges(7,3) * t903 + Ifges(7,4) * t919 + Ifges(7,2) * t918 + Ifges(7,6) * t935 - t925 * t951 + t927 * t971;
t926 = Ifges(7,4) * t951 + Ifges(7,2) * t950 + Ifges(7,6) * t971;
t896 = mrSges(7,2) * t905 - mrSges(7,3) * t902 + Ifges(7,1) * t919 + Ifges(7,4) * t918 + Ifges(7,5) * t935 + t925 * t950 - t926 * t971;
t940 = Ifges(6,5) * t973 + Ifges(6,6) * t972 + Ifges(6,3) * t985;
t941 = Ifges(6,4) * t973 + Ifges(6,2) * t972 + Ifges(6,6) * t985;
t879 = mrSges(6,2) * t911 - mrSges(6,3) * t907 + Ifges(6,1) * t937 + Ifges(6,4) * t936 + Ifges(6,5) * t949 - pkin(12) * t893 - t1030 * t895 + t1035 * t896 + t940 * t972 - t941 * t985;
t1042 = mrSges(7,1) * t902 - mrSges(7,2) * t903 + Ifges(7,5) * t919 + Ifges(7,6) * t918 + Ifges(7,3) * t935 + t926 * t951 - t927 * t950;
t942 = Ifges(6,1) * t973 + Ifges(6,4) * t972 + Ifges(6,5) * t985;
t880 = -mrSges(6,1) * t911 + mrSges(6,3) * t908 + Ifges(6,4) * t937 + Ifges(6,2) * t936 + Ifges(6,6) * t949 - pkin(5) * t893 - t940 * t973 + t942 * t985 - t1042;
t956 = Ifges(5,5) * t987 + Ifges(5,6) * t986 + Ifges(5,3) * t1006;
t957 = Ifges(5,4) * t987 + Ifges(5,2) * t986 + Ifges(5,6) * t1006;
t864 = mrSges(5,2) * t930 - mrSges(5,3) * t913 + Ifges(5,1) * t953 + Ifges(5,4) * t952 + Ifges(5,5) * t997 - pkin(11) * t887 - t1006 * t957 - t1031 * t880 + t1036 * t879 + t956 * t986;
t1041 = mrSges(6,1) * t907 - mrSges(6,2) * t908 + Ifges(6,5) * t937 + Ifges(6,6) * t936 + Ifges(6,3) * t949 + pkin(5) * t904 + pkin(12) * t894 + t1030 * t896 + t1035 * t895 + t973 * t941 - t972 * t942;
t958 = Ifges(5,1) * t987 + Ifges(5,4) * t986 + Ifges(5,5) * t1006;
t870 = -mrSges(5,1) * t930 + mrSges(5,3) * t914 + Ifges(5,4) * t953 + Ifges(5,2) * t952 + Ifges(5,6) * t997 - pkin(4) * t887 + t1006 * t958 - t987 * t956 - t1041;
t976 = Ifges(4,5) * t996 + Ifges(4,6) * t995 + Ifges(4,3) * t1006;
t978 = Ifges(4,1) * t996 + Ifges(4,4) * t995 + Ifges(4,5) * t1006;
t849 = -mrSges(4,1) * t944 + mrSges(4,3) * t932 + Ifges(4,4) * t981 + Ifges(4,2) * t980 + Ifges(4,6) * t997 - pkin(3) * t886 + qJ(4) * t1051 + t1006 * t978 + t1024 * t864 + t1027 * t870 - t996 * t976;
t977 = Ifges(4,4) * t996 + Ifges(4,2) * t995 + Ifges(4,6) * t1006;
t850 = mrSges(4,2) * t944 - mrSges(4,3) * t931 + Ifges(4,1) * t981 + Ifges(4,4) * t980 + Ifges(4,5) * t997 - qJ(4) * t878 - t1006 * t977 - t1024 * t870 + t1027 * t864 + t976 * t995;
t1045 = pkin(10) * t869 + t1032 * t850 + t1037 * t849;
t1000 = Ifges(3,5) * t1021 + (Ifges(3,1) * t1033 + Ifges(3,4) * t1038) * t1066;
t856 = Ifges(4,5) * t981 + Ifges(4,6) * t980 + t996 * t977 - t995 * t978 + mrSges(4,1) * t931 - mrSges(4,2) * t932 + Ifges(5,5) * t953 + Ifges(5,6) * t952 + t987 * t957 - t986 * t958 + mrSges(5,1) * t913 - mrSges(5,2) * t914 + t1031 * t879 + t1036 * t880 + pkin(4) * t1043 + pkin(11) * t1049 + pkin(3) * t878 + (Ifges(4,3) + Ifges(5,3)) * t997;
t999 = Ifges(3,6) * t1021 + (Ifges(3,4) * t1033 + Ifges(3,2) * t1038) * t1066;
t839 = mrSges(3,1) * t992 - mrSges(3,2) * t993 + Ifges(3,5) * t1015 + Ifges(3,6) * t1016 + Ifges(3,3) * t1020 + pkin(2) * t863 + t1028 * t856 + (-t1000 * t1038 + t1033 * t999) * t1066 + t1045 * t1025;
t998 = Ifges(3,3) * t1021 + (Ifges(3,5) * t1033 + Ifges(3,6) * t1038) * t1066;
t841 = -mrSges(3,1) * t1001 + mrSges(3,3) * t993 + Ifges(3,4) * t1015 + Ifges(3,2) * t1016 + Ifges(3,6) * t1020 - pkin(2) * t862 + t1021 * t1000 - t1025 * t856 + t1028 * t1045 - t1053 * t998;
t843 = t998 * t1052 + mrSges(3,2) * t1001 - mrSges(3,3) * t992 + Ifges(3,1) * t1015 + Ifges(3,4) * t1016 + Ifges(3,5) * t1020 - t1021 * t999 - t1032 * t849 + t1037 * t850 + (-t1025 * t862 - t1028 * t863) * pkin(10);
t1044 = mrSges(2,1) * t1018 - mrSges(2,2) * t1019 + Ifges(2,3) * qJDD(1) + pkin(1) * t848 + t1029 * t839 + t841 * t1060 + t843 * t1061 + t1072 * t855;
t837 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1018 + Ifges(2,5) * qJDD(1) - t1040 * Ifges(2,6) - t1033 * t841 + t1038 * t843 + (-t1026 * t847 - t1029 * t848) * pkin(9);
t836 = mrSges(2,1) * g(3) + mrSges(2,3) * t1019 + t1040 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t847 - t1026 * t839 + (pkin(9) * t855 + t1033 * t843 + t1038 * t841) * t1029;
t1 = [-m(1) * g(1) + t1048; -m(1) * g(2) + t1068; (-m(1) - m(2)) * g(3) + t847; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1068 - t1034 * t836 + t1039 * t837; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1048 + t1034 * t837 + t1039 * t836; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1044; t1044; t839; t856; t886; t1041; t1042;];
tauJB  = t1;
