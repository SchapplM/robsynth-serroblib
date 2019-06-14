% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 20:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 20:42:56
% EndTime: 2019-05-06 20:44:03
% DurationCPUTime: 56.52s
% Computational Cost: add. (873446->400), mult. (2279607->518), div. (0->0), fcn. (1833761->14), ass. (0->165)
t1004 = sin(pkin(12));
t1006 = cos(pkin(12));
t1005 = sin(pkin(6));
t1011 = sin(qJ(2));
t1016 = cos(qJ(2));
t1007 = cos(pkin(6));
t1036 = t1007 * t1016;
t1018 = qJD(1) ^ 2;
t1045 = pkin(8) * t1005;
t1012 = sin(qJ(1));
t1017 = cos(qJ(1));
t995 = t1012 * g(1) - t1017 * g(2);
t988 = qJDD(1) * pkin(1) + t1018 * t1045 + t995;
t996 = -t1017 * g(1) - t1012 * g(2);
t989 = -t1018 * pkin(1) + qJDD(1) * t1045 + t996;
t1027 = -t1011 * t989 + t988 * t1036;
t1040 = t1005 ^ 2 * t1018;
t1042 = qJD(1) * qJD(2);
t991 = (qJDD(1) * t1011 + t1016 * t1042) * t1005;
t998 = t1007 * qJDD(1) + qJDD(2);
t999 = t1007 * qJD(1) + qJD(2);
t928 = t998 * pkin(2) - t991 * qJ(3) + (pkin(2) * t1011 * t1040 + (qJ(3) * qJD(1) * t999 - g(3)) * t1005) * t1016 + t1027;
t1033 = t1016 ^ 2 * t1040;
t1037 = t1007 * t1011;
t1039 = t1005 * t1011;
t955 = -g(3) * t1039 + t1016 * t989 + t988 * t1037;
t1035 = qJD(1) * t1039;
t985 = t999 * pkin(2) - qJ(3) * t1035;
t992 = (qJDD(1) * t1016 - t1011 * t1042) * t1005;
t931 = -pkin(2) * t1033 + t992 * qJ(3) - t999 * t985 + t955;
t1041 = qJD(1) * t1005;
t982 = (t1004 * t1016 + t1006 * t1011) * t1041;
t909 = -0.2e1 * qJD(3) * t982 - t1004 * t931 + t1006 * t928;
t1010 = sin(qJ(4));
t1015 = cos(qJ(4));
t1009 = sin(qJ(5));
t1014 = cos(qJ(5));
t1008 = sin(qJ(6));
t1013 = cos(qJ(6));
t1038 = t1005 * t1016;
t1034 = qJD(1) * t1038;
t981 = -t1004 * t1035 + t1006 * t1034;
t910 = 0.2e1 * qJD(3) * t981 + t1004 * t928 + t1006 * t931;
t958 = -t981 * pkin(3) - t982 * pkin(9);
t997 = t999 ^ 2;
t897 = -t997 * pkin(3) + t998 * pkin(9) + t981 * t958 + t910;
t972 = -t1007 * g(3) - t1005 * t988;
t941 = -t992 * pkin(2) - qJ(3) * t1033 + t985 * t1035 + qJDD(3) + t972;
t962 = -t1004 * t991 + t1006 * t992;
t963 = t1004 * t992 + t1006 * t991;
t913 = (-t981 * t999 - t963) * pkin(9) + (t982 * t999 - t962) * pkin(3) + t941;
t889 = -t1010 * t897 + t1015 * t913;
t965 = -t1010 * t982 + t1015 * t999;
t939 = t965 * qJD(4) + t1010 * t998 + t1015 * t963;
t961 = qJDD(4) - t962;
t966 = t1010 * t999 + t1015 * t982;
t980 = qJD(4) - t981;
t886 = (t965 * t980 - t939) * pkin(10) + (t965 * t966 + t961) * pkin(4) + t889;
t890 = t1010 * t913 + t1015 * t897;
t938 = -t966 * qJD(4) - t1010 * t963 + t1015 * t998;
t949 = t980 * pkin(4) - t966 * pkin(10);
t964 = t965 ^ 2;
t888 = -t964 * pkin(4) + t938 * pkin(10) - t980 * t949 + t890;
t883 = t1009 * t886 + t1014 * t888;
t943 = -t1009 * t966 + t1014 * t965;
t944 = t1009 * t965 + t1014 * t966;
t922 = -t943 * pkin(5) - t944 * pkin(11);
t959 = qJDD(5) + t961;
t975 = qJD(5) + t980;
t974 = t975 ^ 2;
t880 = -t974 * pkin(5) + t959 * pkin(11) + t943 * t922 + t883;
t896 = -t998 * pkin(3) - t997 * pkin(9) + t982 * t958 - t909;
t891 = -t938 * pkin(4) - t964 * pkin(10) + t966 * t949 + t896;
t906 = -t944 * qJD(5) - t1009 * t939 + t1014 * t938;
t907 = t943 * qJD(5) + t1009 * t938 + t1014 * t939;
t884 = (-t943 * t975 - t907) * pkin(11) + (t944 * t975 - t906) * pkin(5) + t891;
t877 = -t1008 * t880 + t1013 * t884;
t924 = -t1008 * t944 + t1013 * t975;
t894 = t924 * qJD(6) + t1008 * t959 + t1013 * t907;
t905 = qJDD(6) - t906;
t925 = t1008 * t975 + t1013 * t944;
t914 = -mrSges(7,1) * t924 + mrSges(7,2) * t925;
t942 = qJD(6) - t943;
t915 = -t942 * mrSges(7,2) + t924 * mrSges(7,3);
t873 = m(7) * t877 + t905 * mrSges(7,1) - t894 * mrSges(7,3) - t925 * t914 + t942 * t915;
t878 = t1008 * t884 + t1013 * t880;
t893 = -t925 * qJD(6) - t1008 * t907 + t1013 * t959;
t916 = t942 * mrSges(7,1) - t925 * mrSges(7,3);
t874 = m(7) * t878 - t905 * mrSges(7,2) + t893 * mrSges(7,3) + t924 * t914 - t942 * t916;
t1030 = -t1008 * t873 + t1013 * t874;
t921 = -t943 * mrSges(6,1) + t944 * mrSges(6,2);
t930 = t975 * mrSges(6,1) - t944 * mrSges(6,3);
t860 = m(6) * t883 - t959 * mrSges(6,2) + t906 * mrSges(6,3) + t943 * t921 - t975 * t930 + t1030;
t882 = -t1009 * t888 + t1014 * t886;
t879 = -t959 * pkin(5) - t974 * pkin(11) + t944 * t922 - t882;
t1025 = -m(7) * t879 + t893 * mrSges(7,1) - t894 * mrSges(7,2) + t924 * t915 - t925 * t916;
t929 = -t975 * mrSges(6,2) + t943 * mrSges(6,3);
t869 = m(6) * t882 + t959 * mrSges(6,1) - t907 * mrSges(6,3) - t944 * t921 + t975 * t929 + t1025;
t855 = t1009 * t860 + t1014 * t869;
t945 = -t965 * mrSges(5,1) + t966 * mrSges(5,2);
t947 = -t980 * mrSges(5,2) + t965 * mrSges(5,3);
t853 = m(5) * t889 + t961 * mrSges(5,1) - t939 * mrSges(5,3) - t966 * t945 + t980 * t947 + t855;
t1029 = -t1009 * t869 + t1014 * t860;
t948 = t980 * mrSges(5,1) - t966 * mrSges(5,3);
t854 = m(5) * t890 - t961 * mrSges(5,2) + t938 * mrSges(5,3) + t965 * t945 - t980 * t948 + t1029;
t1028 = -t1010 * t853 + t1015 * t854;
t956 = -t981 * mrSges(4,1) + t982 * mrSges(4,2);
t968 = t999 * mrSges(4,1) - t982 * mrSges(4,3);
t843 = m(4) * t910 - t998 * mrSges(4,2) + t962 * mrSges(4,3) + t981 * t956 - t999 * t968 + t1028;
t862 = t1008 * t874 + t1013 * t873;
t1023 = m(6) * t891 - t906 * mrSges(6,1) + t907 * mrSges(6,2) - t943 * t929 + t944 * t930 + t862;
t1020 = -m(5) * t896 + t938 * mrSges(5,1) - t939 * mrSges(5,2) + t965 * t947 - t966 * t948 - t1023;
t967 = -t999 * mrSges(4,2) + t981 * mrSges(4,3);
t857 = m(4) * t909 + t998 * mrSges(4,1) - t963 * mrSges(4,3) - t982 * t956 + t999 * t967 + t1020;
t840 = t1004 * t843 + t1006 * t857;
t954 = -g(3) * t1038 + t1027;
t987 = -t999 * mrSges(3,2) + mrSges(3,3) * t1034;
t990 = (-mrSges(3,1) * t1016 + mrSges(3,2) * t1011) * t1041;
t838 = m(3) * t954 + t998 * mrSges(3,1) - t991 * mrSges(3,3) - t990 * t1035 + t999 * t987 + t840;
t1032 = -t1004 * t857 + t1006 * t843;
t986 = t999 * mrSges(3,1) - mrSges(3,3) * t1035;
t839 = m(3) * t955 - t998 * mrSges(3,2) + t992 * mrSges(3,3) + t990 * t1034 - t999 * t986 + t1032;
t847 = t1010 * t854 + t1015 * t853;
t846 = m(4) * t941 - t962 * mrSges(4,1) + t963 * mrSges(4,2) - t981 * t967 + t982 * t968 + t847;
t845 = m(3) * t972 - t992 * mrSges(3,1) + t991 * mrSges(3,2) + (t1011 * t986 - t1016 * t987) * t1041 + t846;
t824 = -t1005 * t845 + t838 * t1036 + t839 * t1037;
t821 = m(2) * t995 + qJDD(1) * mrSges(2,1) - t1018 * mrSges(2,2) + t824;
t829 = -t1011 * t838 + t1016 * t839;
t827 = m(2) * t996 - t1018 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t829;
t1044 = t1012 * t827 + t1017 * t821;
t823 = t1007 * t845 + t838 * t1038 + t839 * t1039;
t1026 = -t1012 * t821 + t1017 * t827;
t898 = Ifges(7,5) * t925 + Ifges(7,6) * t924 + Ifges(7,3) * t942;
t900 = Ifges(7,1) * t925 + Ifges(7,4) * t924 + Ifges(7,5) * t942;
t866 = -mrSges(7,1) * t879 + mrSges(7,3) * t878 + Ifges(7,4) * t894 + Ifges(7,2) * t893 + Ifges(7,6) * t905 - t925 * t898 + t942 * t900;
t899 = Ifges(7,4) * t925 + Ifges(7,2) * t924 + Ifges(7,6) * t942;
t867 = mrSges(7,2) * t879 - mrSges(7,3) * t877 + Ifges(7,1) * t894 + Ifges(7,4) * t893 + Ifges(7,5) * t905 + t924 * t898 - t942 * t899;
t917 = Ifges(6,5) * t944 + Ifges(6,6) * t943 + Ifges(6,3) * t975;
t918 = Ifges(6,4) * t944 + Ifges(6,2) * t943 + Ifges(6,6) * t975;
t848 = mrSges(6,2) * t891 - mrSges(6,3) * t882 + Ifges(6,1) * t907 + Ifges(6,4) * t906 + Ifges(6,5) * t959 - pkin(11) * t862 - t1008 * t866 + t1013 * t867 + t943 * t917 - t975 * t918;
t1021 = mrSges(7,1) * t877 - mrSges(7,2) * t878 + Ifges(7,5) * t894 + Ifges(7,6) * t893 + Ifges(7,3) * t905 + t925 * t899 - t924 * t900;
t919 = Ifges(6,1) * t944 + Ifges(6,4) * t943 + Ifges(6,5) * t975;
t849 = -mrSges(6,1) * t891 + mrSges(6,3) * t883 + Ifges(6,4) * t907 + Ifges(6,2) * t906 + Ifges(6,6) * t959 - pkin(5) * t862 - t944 * t917 + t975 * t919 - t1021;
t932 = Ifges(5,5) * t966 + Ifges(5,6) * t965 + Ifges(5,3) * t980;
t934 = Ifges(5,1) * t966 + Ifges(5,4) * t965 + Ifges(5,5) * t980;
t831 = -mrSges(5,1) * t896 + mrSges(5,3) * t890 + Ifges(5,4) * t939 + Ifges(5,2) * t938 + Ifges(5,6) * t961 - pkin(4) * t1023 + pkin(10) * t1029 + t1009 * t848 + t1014 * t849 - t966 * t932 + t980 * t934;
t933 = Ifges(5,4) * t966 + Ifges(5,2) * t965 + Ifges(5,6) * t980;
t832 = mrSges(5,2) * t896 - mrSges(5,3) * t889 + Ifges(5,1) * t939 + Ifges(5,4) * t938 + Ifges(5,5) * t961 - pkin(10) * t855 - t1009 * t849 + t1014 * t848 + t965 * t932 - t980 * t933;
t950 = Ifges(4,5) * t982 + Ifges(4,6) * t981 + Ifges(4,3) * t999;
t951 = Ifges(4,4) * t982 + Ifges(4,2) * t981 + Ifges(4,6) * t999;
t819 = mrSges(4,2) * t941 - mrSges(4,3) * t909 + Ifges(4,1) * t963 + Ifges(4,4) * t962 + Ifges(4,5) * t998 - pkin(9) * t847 - t1010 * t831 + t1015 * t832 + t981 * t950 - t999 * t951;
t1022 = -mrSges(6,1) * t882 + mrSges(6,2) * t883 - Ifges(6,5) * t907 - Ifges(6,6) * t906 - Ifges(6,3) * t959 - pkin(5) * t1025 - pkin(11) * t1030 - t1008 * t867 - t1013 * t866 - t944 * t918 + t943 * t919;
t1019 = mrSges(5,1) * t889 - mrSges(5,2) * t890 + Ifges(5,5) * t939 + Ifges(5,6) * t938 + Ifges(5,3) * t961 + pkin(4) * t855 + t966 * t933 - t965 * t934 - t1022;
t952 = Ifges(4,1) * t982 + Ifges(4,4) * t981 + Ifges(4,5) * t999;
t830 = -mrSges(4,1) * t941 + mrSges(4,3) * t910 + Ifges(4,4) * t963 + Ifges(4,2) * t962 + Ifges(4,6) * t998 - pkin(3) * t847 - t982 * t950 + t999 * t952 - t1019;
t969 = Ifges(3,3) * t999 + (Ifges(3,5) * t1011 + Ifges(3,6) * t1016) * t1041;
t971 = Ifges(3,5) * t999 + (Ifges(3,1) * t1011 + Ifges(3,4) * t1016) * t1041;
t814 = -mrSges(3,1) * t972 + mrSges(3,3) * t955 + Ifges(3,4) * t991 + Ifges(3,2) * t992 + Ifges(3,6) * t998 - pkin(2) * t846 + qJ(3) * t1032 + t1004 * t819 + t1006 * t830 - t969 * t1035 + t999 * t971;
t970 = Ifges(3,6) * t999 + (Ifges(3,4) * t1011 + Ifges(3,2) * t1016) * t1041;
t816 = mrSges(3,2) * t972 - mrSges(3,3) * t954 + Ifges(3,1) * t991 + Ifges(3,4) * t992 + Ifges(3,5) * t998 - qJ(3) * t840 - t1004 * t830 + t1006 * t819 + t969 * t1034 - t999 * t970;
t818 = Ifges(3,5) * t991 + Ifges(3,6) * t992 + mrSges(3,1) * t954 - mrSges(3,2) * t955 + Ifges(4,5) * t963 + Ifges(4,6) * t962 + t982 * t951 - t981 * t952 + mrSges(4,1) * t909 - mrSges(4,2) * t910 + t1010 * t832 + t1015 * t831 + pkin(3) * t1020 + pkin(9) * t1028 + pkin(2) * t840 + (Ifges(3,3) + Ifges(4,3)) * t998 + (t1011 * t970 - t1016 * t971) * t1041;
t1024 = mrSges(2,1) * t995 - mrSges(2,2) * t996 + Ifges(2,3) * qJDD(1) + pkin(1) * t824 + t1007 * t818 + t814 * t1038 + t816 * t1039 + t829 * t1045;
t812 = -mrSges(2,2) * g(3) - mrSges(2,3) * t995 + Ifges(2,5) * qJDD(1) - t1018 * Ifges(2,6) - t1011 * t814 + t1016 * t816 + (-t1005 * t823 - t1007 * t824) * pkin(8);
t811 = mrSges(2,1) * g(3) + mrSges(2,3) * t996 + t1018 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t823 - t1005 * t818 + (pkin(8) * t829 + t1011 * t816 + t1016 * t814) * t1007;
t1 = [-m(1) * g(1) + t1026; -m(1) * g(2) + t1044; (-m(1) - m(2)) * g(3) + t823; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1044 - t1012 * t811 + t1017 * t812; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1026 + t1012 * t812 + t1017 * t811; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1024; t1024; t818; t846; t1019; -t1022; t1021;];
tauJB  = t1;
