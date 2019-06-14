% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-08 02:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR14_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR14_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR14_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 02:18:15
% EndTime: 2019-05-08 02:21:18
% DurationCPUTime: 179.34s
% Computational Cost: add. (3007278->418), mult. (7456738->557), div. (0->0), fcn. (6323127->16), ass. (0->180)
t1028 = sin(pkin(7));
t1031 = cos(pkin(7));
t1035 = sin(qJ(3));
t1039 = cos(qJ(3));
t1032 = cos(pkin(6));
t1024 = qJD(1) * t1032 + qJD(2);
t1029 = sin(pkin(6));
t1040 = cos(qJ(2));
t1066 = t1029 * t1040;
t1059 = qJD(1) * t1066;
t1008 = (t1024 * t1028 + t1031 * t1059) * pkin(10);
t1036 = sin(qJ(2));
t1072 = qJD(1) * t1029;
t1076 = pkin(10) * t1028;
t1012 = (-pkin(2) * t1040 - t1036 * t1076) * t1072;
t1070 = qJD(1) * t1040;
t1018 = (qJD(2) * t1070 + qJDD(1) * t1036) * t1029;
t1023 = qJDD(1) * t1032 + qJDD(2);
t1037 = sin(qJ(1));
t1041 = cos(qJ(1));
t1021 = t1037 * g(1) - g(2) * t1041;
t1042 = qJD(1) ^ 2;
t1077 = pkin(9) * t1029;
t1015 = qJDD(1) * pkin(1) + t1042 * t1077 + t1021;
t1022 = -g(1) * t1041 - g(2) * t1037;
t1016 = -pkin(1) * t1042 + qJDD(1) * t1077 + t1022;
t1062 = t1032 * t1040;
t1054 = t1015 * t1062 - t1036 * t1016;
t1071 = qJD(1) * t1036;
t1075 = pkin(10) * t1031;
t964 = -t1018 * t1075 + t1023 * pkin(2) + t1024 * t1008 + (-g(3) * t1040 - t1012 * t1071) * t1029 + t1054;
t1067 = t1029 * t1036;
t1060 = qJD(1) * t1067;
t1011 = pkin(2) * t1024 - t1060 * t1075;
t1019 = (-qJD(2) * t1071 + qJDD(1) * t1040) * t1029;
t1050 = t1019 * t1031 + t1023 * t1028;
t1063 = t1032 * t1036;
t1061 = t1015 * t1063 + t1040 * t1016;
t965 = -t1024 * t1011 + (-g(3) * t1036 + t1012 * t1070) * t1029 + t1050 * pkin(10) + t1061;
t1074 = t1032 * g(3);
t970 = -t1018 * t1076 - t1019 * pkin(2) - t1074 + (-t1015 + (-t1008 * t1040 + t1011 * t1036) * qJD(1)) * t1029;
t935 = -t1035 * t965 + (t1028 * t970 + t1031 * t964) * t1039;
t1064 = t1031 * t1040;
t1069 = t1028 * t1035;
t998 = t1024 * t1069 + (t1035 * t1064 + t1036 * t1039) * t1072;
t983 = -t998 * qJD(3) - t1035 * t1018 + t1050 * t1039;
t1068 = t1028 * t1039;
t997 = (-t1035 * t1036 + t1039 * t1064) * t1072 + t1024 * t1068;
t1078 = cos(qJ(4));
t1014 = -mrSges(3,2) * t1024 + mrSges(3,3) * t1059;
t1017 = (-mrSges(3,1) * t1040 + mrSges(3,2) * t1036) * t1072;
t1065 = t1031 * t1035;
t1009 = t1024 * t1031 - t1028 * t1059 + qJD(3);
t1034 = sin(qJ(4));
t1027 = sin(pkin(13));
t1030 = cos(pkin(13));
t989 = t1034 * t1009 + t1078 * t998;
t995 = qJD(4) - t997;
t975 = -t1027 * t989 + t1030 * t995;
t988 = -t1078 * t1009 + t1034 * t998;
t1053 = -mrSges(6,2) * t988 + mrSges(6,3) * t975;
t1033 = sin(qJ(6));
t1038 = cos(qJ(6));
t1007 = t1009 ^ 2;
t936 = t1039 * t965 + t964 * t1065 + t970 * t1069;
t986 = -pkin(3) * t997 - pkin(11) * t998;
t999 = -t1019 * t1028 + t1023 * t1031 + qJDD(3);
t927 = -pkin(3) * t1007 + pkin(11) * t999 + t986 * t997 + t936;
t946 = -t1028 * t964 + t1031 * t970;
t984 = t997 * qJD(3) + t1039 * t1018 + t1050 * t1035;
t929 = (-t1009 * t997 - t984) * pkin(11) + (t1009 * t998 - t983) * pkin(3) + t946;
t920 = t1034 * t929 + t1078 * t927;
t966 = pkin(4) * t988 - qJ(5) * t989;
t982 = qJDD(4) - t983;
t994 = t995 ^ 2;
t915 = -pkin(4) * t994 + qJ(5) * t982 - t966 * t988 + t920;
t926 = -t999 * pkin(3) - t1007 * pkin(11) + t998 * t986 - t935;
t952 = qJD(4) * t989 + t1034 * t984 - t1078 * t999;
t953 = -t988 * qJD(4) + t1034 * t999 + t1078 * t984;
t918 = (t988 * t995 - t953) * qJ(5) + (t989 * t995 + t952) * pkin(4) + t926;
t976 = t1027 * t995 + t1030 * t989;
t910 = -0.2e1 * qJD(5) * t976 - t1027 * t915 + t1030 * t918;
t939 = t1027 * t982 + t1030 * t953;
t908 = (t975 * t988 - t939) * pkin(12) + (t975 * t976 + t952) * pkin(5) + t910;
t911 = 0.2e1 * qJD(5) * t975 + t1027 * t918 + t1030 * t915;
t938 = -t1027 * t953 + t1030 * t982;
t956 = pkin(5) * t988 - pkin(12) * t976;
t974 = t975 ^ 2;
t909 = -pkin(5) * t974 + pkin(12) * t938 - t956 * t988 + t911;
t906 = -t1033 * t909 + t1038 * t908;
t947 = -t1033 * t976 + t1038 * t975;
t923 = qJD(6) * t947 + t1033 * t938 + t1038 * t939;
t948 = t1033 * t975 + t1038 * t976;
t934 = -mrSges(7,1) * t947 + mrSges(7,2) * t948;
t987 = qJD(6) + t988;
t940 = -mrSges(7,2) * t987 + mrSges(7,3) * t947;
t951 = qJDD(6) + t952;
t903 = m(7) * t906 + mrSges(7,1) * t951 - mrSges(7,3) * t923 - t934 * t948 + t940 * t987;
t907 = t1033 * t908 + t1038 * t909;
t922 = -qJD(6) * t948 - t1033 * t939 + t1038 * t938;
t941 = mrSges(7,1) * t987 - mrSges(7,3) * t948;
t904 = m(7) * t907 - mrSges(7,2) * t951 + mrSges(7,3) * t922 + t934 * t947 - t941 * t987;
t895 = t1033 * t904 + t1038 * t903;
t949 = -mrSges(6,1) * t975 + mrSges(6,2) * t976;
t893 = m(6) * t910 + t952 * mrSges(6,1) - t939 * mrSges(6,3) + t988 * t1053 - t976 * t949 + t895;
t1057 = -t1033 * t903 + t1038 * t904;
t955 = mrSges(6,1) * t988 - mrSges(6,3) * t976;
t894 = m(6) * t911 - mrSges(6,2) * t952 + mrSges(6,3) * t938 + t949 * t975 - t955 * t988 + t1057;
t891 = -t1027 * t893 + t1030 * t894;
t967 = mrSges(5,1) * t988 + mrSges(5,2) * t989;
t978 = mrSges(5,1) * t995 - mrSges(5,3) * t989;
t889 = m(5) * t920 - mrSges(5,2) * t982 - mrSges(5,3) * t952 - t967 * t988 - t978 * t995 + t891;
t919 = -t1034 * t927 + t1078 * t929;
t914 = -t982 * pkin(4) - t994 * qJ(5) + t989 * t966 + qJDD(5) - t919;
t912 = -t938 * pkin(5) - t974 * pkin(12) + t976 * t956 + t914;
t1047 = m(7) * t912 - t922 * mrSges(7,1) + mrSges(7,2) * t923 - t947 * t940 + t941 * t948;
t905 = m(6) * t914 - t938 * mrSges(6,1) + mrSges(6,2) * t939 - t975 * t1053 + t955 * t976 + t1047;
t977 = -mrSges(5,2) * t995 - mrSges(5,3) * t988;
t899 = m(5) * t919 + mrSges(5,1) * t982 - mrSges(5,3) * t953 - t967 * t989 + t977 * t995 - t905;
t1056 = -t1034 * t899 + t1078 * t889;
t985 = -mrSges(4,1) * t997 + mrSges(4,2) * t998;
t991 = mrSges(4,1) * t1009 - mrSges(4,3) * t998;
t878 = m(4) * t936 - mrSges(4,2) * t999 + mrSges(4,3) * t983 - t1009 * t991 + t985 * t997 + t1056;
t881 = t1034 * t889 + t1078 * t899;
t990 = -mrSges(4,2) * t1009 + mrSges(4,3) * t997;
t880 = m(4) * t946 - mrSges(4,1) * t983 + mrSges(4,2) * t984 - t990 * t997 + t991 * t998 + t881;
t890 = t1027 * t894 + t1030 * t893;
t1045 = -m(5) * t926 - t952 * mrSges(5,1) - mrSges(5,2) * t953 - t988 * t977 - t978 * t989 - t890;
t886 = m(4) * t935 + mrSges(4,1) * t999 - mrSges(4,3) * t984 + t1009 * t990 - t985 * t998 + t1045;
t867 = t1031 * t1039 * t886 - t1028 * t880 + t878 * t1065;
t992 = -g(3) * t1066 + t1054;
t863 = m(3) * t992 + mrSges(3,1) * t1023 - mrSges(3,3) * t1018 + t1014 * t1024 - t1017 * t1060 + t867;
t1003 = -t1029 * t1015 - t1074;
t1013 = mrSges(3,1) * t1024 - mrSges(3,3) * t1060;
t866 = t1031 * t880 + t886 * t1068 + t878 * t1069;
t865 = m(3) * t1003 - t1019 * mrSges(3,1) + t1018 * mrSges(3,2) + (t1013 * t1036 - t1014 * t1040) * t1072 + t866;
t873 = -t1035 * t886 + t1039 * t878;
t993 = -g(3) * t1067 + t1061;
t872 = m(3) * t993 - mrSges(3,2) * t1023 + mrSges(3,3) * t1019 - t1013 * t1024 + t1017 * t1059 + t873;
t852 = -t1029 * t865 + t863 * t1062 + t872 * t1063;
t849 = m(2) * t1021 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1042 + t852;
t859 = -t1036 * t863 + t1040 * t872;
t857 = m(2) * t1022 - mrSges(2,1) * t1042 - qJDD(1) * mrSges(2,2) + t859;
t1073 = t1037 * t857 + t1041 * t849;
t851 = t1032 * t865 + t863 * t1066 + t872 * t1067;
t1055 = -t1037 * t849 + t1041 * t857;
t930 = Ifges(7,5) * t948 + Ifges(7,6) * t947 + Ifges(7,3) * t987;
t932 = Ifges(7,1) * t948 + Ifges(7,4) * t947 + Ifges(7,5) * t987;
t896 = -mrSges(7,1) * t912 + mrSges(7,3) * t907 + Ifges(7,4) * t923 + Ifges(7,2) * t922 + Ifges(7,6) * t951 - t930 * t948 + t932 * t987;
t931 = Ifges(7,4) * t948 + Ifges(7,2) * t947 + Ifges(7,6) * t987;
t897 = mrSges(7,2) * t912 - mrSges(7,3) * t906 + Ifges(7,1) * t923 + Ifges(7,4) * t922 + Ifges(7,5) * t951 + t930 * t947 - t931 * t987;
t942 = Ifges(6,5) * t976 + Ifges(6,6) * t975 + Ifges(6,3) * t988;
t944 = Ifges(6,1) * t976 + Ifges(6,4) * t975 + Ifges(6,5) * t988;
t882 = -mrSges(6,1) * t914 + mrSges(6,3) * t911 + Ifges(6,4) * t939 + Ifges(6,2) * t938 + Ifges(6,6) * t952 - pkin(5) * t1047 + pkin(12) * t1057 + t1033 * t897 + t1038 * t896 - t976 * t942 + t988 * t944;
t943 = Ifges(6,4) * t976 + Ifges(6,2) * t975 + Ifges(6,6) * t988;
t883 = mrSges(6,2) * t914 - mrSges(6,3) * t910 + Ifges(6,1) * t939 + Ifges(6,4) * t938 + Ifges(6,5) * t952 - pkin(12) * t895 - t1033 * t896 + t1038 * t897 + t942 * t975 - t943 * t988;
t957 = Ifges(5,5) * t989 - Ifges(5,6) * t988 + Ifges(5,3) * t995;
t958 = Ifges(5,4) * t989 - Ifges(5,2) * t988 + Ifges(5,6) * t995;
t868 = mrSges(5,2) * t926 - mrSges(5,3) * t919 + Ifges(5,1) * t953 - Ifges(5,4) * t952 + Ifges(5,5) * t982 - qJ(5) * t890 - t1027 * t882 + t1030 * t883 - t957 * t988 - t958 * t995;
t1044 = mrSges(7,1) * t906 - mrSges(7,2) * t907 + Ifges(7,5) * t923 + Ifges(7,6) * t922 + Ifges(7,3) * t951 + t948 * t931 - t947 * t932;
t959 = Ifges(5,1) * t989 - Ifges(5,4) * t988 + Ifges(5,5) * t995;
t874 = -t1044 - pkin(5) * t895 - t989 * t957 - Ifges(6,6) * t938 - mrSges(5,1) * t926 + mrSges(5,3) * t920 + mrSges(6,2) * t911 - mrSges(6,1) * t910 + t975 * t944 - t976 * t943 + t995 * t959 - pkin(4) * t890 + Ifges(5,4) * t953 - Ifges(6,5) * t939 + Ifges(5,6) * t982 + (-Ifges(5,2) - Ifges(6,3)) * t952;
t979 = Ifges(4,5) * t998 + Ifges(4,6) * t997 + Ifges(4,3) * t1009;
t980 = Ifges(4,4) * t998 + Ifges(4,2) * t997 + Ifges(4,6) * t1009;
t854 = mrSges(4,2) * t946 - mrSges(4,3) * t935 + Ifges(4,1) * t984 + Ifges(4,4) * t983 + Ifges(4,5) * t999 - pkin(11) * t881 - t1009 * t980 - t1034 * t874 + t1078 * t868 + t997 * t979;
t1043 = mrSges(5,1) * t919 - mrSges(5,2) * t920 + Ifges(5,5) * t953 - Ifges(5,6) * t952 + Ifges(5,3) * t982 - pkin(4) * t905 + qJ(5) * t891 + t1027 * t883 + t1030 * t882 + t989 * t958 + t988 * t959;
t981 = Ifges(4,1) * t998 + Ifges(4,4) * t997 + Ifges(4,5) * t1009;
t860 = -mrSges(4,1) * t946 + mrSges(4,3) * t936 + Ifges(4,4) * t984 + Ifges(4,2) * t983 + Ifges(4,6) * t999 - pkin(3) * t881 + t1009 * t981 - t998 * t979 - t1043;
t1048 = pkin(10) * t873 + t1035 * t854 + t1039 * t860;
t1001 = Ifges(3,6) * t1024 + (Ifges(3,4) * t1036 + Ifges(3,2) * t1040) * t1072;
t1002 = Ifges(3,5) * t1024 + (Ifges(3,1) * t1036 + Ifges(3,4) * t1040) * t1072;
t853 = mrSges(4,1) * t935 - mrSges(4,2) * t936 + Ifges(4,5) * t984 + Ifges(4,6) * t983 + Ifges(4,3) * t999 + pkin(3) * t1045 + pkin(11) * t1056 + t1034 * t868 + t1078 * t874 + t998 * t980 - t997 * t981;
t843 = mrSges(3,1) * t992 - mrSges(3,2) * t993 + Ifges(3,5) * t1018 + Ifges(3,6) * t1019 + Ifges(3,3) * t1023 + pkin(2) * t867 + t1031 * t853 + (t1001 * t1036 - t1002 * t1040) * t1072 + t1048 * t1028;
t1000 = Ifges(3,3) * t1024 + (Ifges(3,5) * t1036 + Ifges(3,6) * t1040) * t1072;
t845 = -mrSges(3,1) * t1003 + mrSges(3,3) * t993 + Ifges(3,4) * t1018 + Ifges(3,2) * t1019 + Ifges(3,6) * t1023 - pkin(2) * t866 - t1000 * t1060 + t1024 * t1002 - t1028 * t853 + t1048 * t1031;
t847 = t1000 * t1059 + mrSges(3,2) * t1003 - mrSges(3,3) * t992 + Ifges(3,1) * t1018 + Ifges(3,4) * t1019 + Ifges(3,5) * t1023 - t1024 * t1001 - t1035 * t860 + t1039 * t854 + (-t1028 * t866 - t1031 * t867) * pkin(10);
t1046 = mrSges(2,1) * t1021 - mrSges(2,2) * t1022 + Ifges(2,3) * qJDD(1) + pkin(1) * t852 + t1032 * t843 + t845 * t1066 + t847 * t1067 + t859 * t1077;
t841 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1021 + Ifges(2,5) * qJDD(1) - t1042 * Ifges(2,6) - t1036 * t845 + t1040 * t847 + (-t1029 * t851 - t1032 * t852) * pkin(9);
t840 = mrSges(2,1) * g(3) + mrSges(2,3) * t1022 + t1042 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t851 - t1029 * t843 + (pkin(9) * t859 + t1036 * t847 + t1040 * t845) * t1032;
t1 = [-m(1) * g(1) + t1055; -m(1) * g(2) + t1073; (-m(1) - m(2)) * g(3) + t851; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1073 - t1037 * t840 + t1041 * t841; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1055 + t1037 * t841 + t1041 * t840; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1046; t1046; t843; t853; t1043; t905; t1044;];
tauJB  = t1;
