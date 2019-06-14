% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 21:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 21:15:06
% EndTime: 2019-05-07 21:16:50
% DurationCPUTime: 78.33s
% Computational Cost: add. (1333477->402), mult. (2922338->518), div. (0->0), fcn. (2396547->14), ass. (0->168)
t1046 = 2 * qJD(5);
t1006 = sin(pkin(6));
t1045 = pkin(8) * t1006;
t1008 = cos(pkin(6));
t1044 = g(3) * t1008;
t1013 = sin(qJ(1));
t1018 = cos(qJ(1));
t1019 = qJD(1) ^ 2;
t1017 = cos(qJ(2));
t1036 = t1008 * t1017;
t1012 = sin(qJ(2));
t1037 = t1008 * t1012;
t1000 = qJDD(1) * t1008 + qJDD(2);
t1001 = qJD(1) * t1008 + qJD(2);
t1011 = sin(qJ(3));
t1016 = cos(qJ(3));
t1010 = sin(qJ(4));
t1015 = cos(qJ(4));
t1005 = sin(pkin(12));
t1007 = cos(pkin(12));
t1009 = sin(qJ(6));
t1014 = cos(qJ(6));
t1040 = qJD(1) * t1017;
t996 = t1013 * g(1) - g(2) * t1018;
t986 = qJDD(1) * pkin(1) + t1019 * t1045 + t996;
t1035 = qJDD(1) * t1006;
t997 = -g(1) * t1018 - g(2) * t1013;
t987 = -pkin(1) * t1019 + pkin(8) * t1035 + t997;
t1042 = t1017 * t987 + t986 * t1037;
t1041 = qJD(1) * t1006;
t989 = (-pkin(2) * t1017 - pkin(9) * t1012) * t1041;
t999 = t1001 ^ 2;
t947 = -pkin(2) * t999 + pkin(9) * t1000 + (-g(3) * t1012 + t989 * t1040) * t1006 + t1042;
t990 = (qJD(2) * t1040 + qJDD(1) * t1012) * t1006;
t1039 = t1006 * t1012;
t1034 = qJD(1) * t1039;
t991 = -qJD(2) * t1034 + t1017 * t1035;
t948 = -pkin(2) * t991 - pkin(9) * t990 - t1044 + (-t986 + (pkin(2) * t1012 - pkin(9) * t1017) * t1001 * qJD(1)) * t1006;
t921 = -t1011 * t947 + t1016 * t948;
t978 = t1001 * t1016 - t1011 * t1034;
t959 = qJD(3) * t978 + t1000 * t1011 + t1016 * t990;
t979 = t1001 * t1011 + t1016 * t1034;
t983 = qJDD(3) - t991;
t1038 = t1006 * t1017;
t1033 = qJD(1) * t1038;
t995 = qJD(3) - t1033;
t899 = (t978 * t995 - t959) * pkin(10) + (t978 * t979 + t983) * pkin(3) + t921;
t922 = t1011 * t948 + t1016 * t947;
t958 = -qJD(3) * t979 + t1000 * t1016 - t1011 * t990;
t968 = pkin(3) * t995 - pkin(10) * t979;
t977 = t978 ^ 2;
t909 = -pkin(3) * t977 + pkin(10) * t958 - t968 * t995 + t922;
t886 = -t1010 * t909 + t1015 * t899;
t963 = -t1010 * t979 + t1015 * t978;
t927 = qJD(4) * t963 + t1010 * t958 + t1015 * t959;
t964 = t1010 * t978 + t1015 * t979;
t982 = qJDD(4) + t983;
t993 = qJD(4) + t995;
t882 = (t963 * t993 - t927) * qJ(5) + (t963 * t964 + t982) * pkin(4) + t886;
t887 = t1010 * t899 + t1015 * t909;
t926 = -qJD(4) * t964 - t1010 * t959 + t1015 * t958;
t950 = pkin(4) * t993 - qJ(5) * t964;
t962 = t963 ^ 2;
t884 = -pkin(4) * t962 + qJ(5) * t926 - t950 * t993 + t887;
t940 = -t1005 * t964 + t1007 * t963;
t879 = t1005 * t882 + t1007 * t884 + t940 * t1046;
t941 = t1005 * t963 + t1007 * t964;
t920 = -pkin(5) * t940 - pkin(11) * t941;
t992 = t993 ^ 2;
t876 = -pkin(5) * t992 + pkin(11) * t982 + t920 * t940 + t879;
t960 = -g(3) * t1038 - t1012 * t987 + t986 * t1036;
t946 = -pkin(2) * t1000 - pkin(9) * t999 + t989 * t1034 - t960;
t913 = -pkin(3) * t958 - pkin(10) * t977 + t979 * t968 + t946;
t889 = -pkin(4) * t926 - qJ(5) * t962 + t964 * t950 + qJDD(5) + t913;
t905 = -t1005 * t927 + t1007 * t926;
t906 = t1005 * t926 + t1007 * t927;
t880 = t889 + (t941 * t993 - t905) * pkin(5) + (-t940 * t993 - t906) * pkin(11);
t873 = -t1009 * t876 + t1014 * t880;
t929 = -t1009 * t941 + t1014 * t993;
t892 = qJD(6) * t929 + t1009 * t982 + t1014 * t906;
t904 = qJDD(6) - t905;
t930 = t1009 * t993 + t1014 * t941;
t910 = -mrSges(7,1) * t929 + mrSges(7,2) * t930;
t939 = qJD(6) - t940;
t911 = -mrSges(7,2) * t939 + mrSges(7,3) * t929;
t869 = m(7) * t873 + mrSges(7,1) * t904 - mrSges(7,3) * t892 - t910 * t930 + t911 * t939;
t874 = t1009 * t880 + t1014 * t876;
t891 = -qJD(6) * t930 - t1009 * t906 + t1014 * t982;
t912 = mrSges(7,1) * t939 - mrSges(7,3) * t930;
t870 = m(7) * t874 - mrSges(7,2) * t904 + mrSges(7,3) * t891 + t910 * t929 - t912 * t939;
t1031 = -t1009 * t869 + t1014 * t870;
t919 = -mrSges(6,1) * t940 + mrSges(6,2) * t941;
t932 = mrSges(6,1) * t993 - mrSges(6,3) * t941;
t855 = m(6) * t879 - mrSges(6,2) * t982 + mrSges(6,3) * t905 + t919 * t940 - t932 * t993 + t1031;
t1027 = t1005 * t884 - t1007 * t882;
t875 = -pkin(5) * t982 - pkin(11) * t992 + (t1046 + t920) * t941 + t1027;
t1026 = -m(7) * t875 + t891 * mrSges(7,1) - mrSges(7,2) * t892 + t929 * t911 - t912 * t930;
t878 = -0.2e1 * qJD(5) * t941 - t1027;
t931 = -mrSges(6,2) * t993 + mrSges(6,3) * t940;
t865 = m(6) * t878 + mrSges(6,1) * t982 - mrSges(6,3) * t906 - t919 * t941 + t931 * t993 + t1026;
t849 = t1005 * t855 + t1007 * t865;
t942 = -mrSges(5,1) * t963 + mrSges(5,2) * t964;
t949 = -mrSges(5,2) * t993 + mrSges(5,3) * t963;
t846 = m(5) * t886 + mrSges(5,1) * t982 - mrSges(5,3) * t927 - t942 * t964 + t949 * t993 + t849;
t1032 = -t1005 * t865 + t1007 * t855;
t951 = mrSges(5,1) * t993 - mrSges(5,3) * t964;
t847 = m(5) * t887 - mrSges(5,2) * t982 + mrSges(5,3) * t926 + t942 * t963 - t951 * t993 + t1032;
t840 = t1010 * t847 + t1015 * t846;
t965 = -mrSges(4,1) * t978 + mrSges(4,2) * t979;
t966 = -mrSges(4,2) * t995 + mrSges(4,3) * t978;
t838 = m(4) * t921 + mrSges(4,1) * t983 - mrSges(4,3) * t959 - t965 * t979 + t966 * t995 + t840;
t1030 = -t1010 * t846 + t1015 * t847;
t967 = mrSges(4,1) * t995 - mrSges(4,3) * t979;
t839 = m(4) * t922 - mrSges(4,2) * t983 + mrSges(4,3) * t958 + t965 * t978 - t967 * t995 + t1030;
t1029 = -t1011 * t838 + t1016 * t839;
t961 = -g(3) * t1039 + t1042;
t984 = mrSges(3,1) * t1001 - mrSges(3,3) * t1034;
t988 = (-mrSges(3,1) * t1017 + mrSges(3,2) * t1012) * t1041;
t830 = m(3) * t961 - mrSges(3,2) * t1000 + mrSges(3,3) * t991 - t1001 * t984 + t988 * t1033 + t1029;
t833 = t1011 * t839 + t1016 * t838;
t972 = -t1006 * t986 - t1044;
t985 = -mrSges(3,2) * t1001 + mrSges(3,3) * t1033;
t832 = m(3) * t972 - mrSges(3,1) * t991 + mrSges(3,2) * t990 + (t1012 * t984 - t1017 * t985) * t1041 + t833;
t858 = t1009 * t870 + t1014 * t869;
t856 = m(6) * t889 - t905 * mrSges(6,1) + t906 * mrSges(6,2) - t940 * t931 + t941 * t932 + t858;
t1024 = m(5) * t913 - t926 * mrSges(5,1) + t927 * mrSges(5,2) - t963 * t949 + t964 * t951 + t856;
t1021 = -m(4) * t946 + t958 * mrSges(4,1) - t959 * mrSges(4,2) + t978 * t966 - t979 * t967 - t1024;
t852 = m(3) * t960 + t1000 * mrSges(3,1) - t990 * mrSges(3,3) + t1001 * t985 - t988 * t1034 + t1021;
t820 = -t1006 * t832 + t852 * t1036 + t830 * t1037;
t817 = m(2) * t996 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1019 + t820;
t825 = -t1012 * t852 + t1017 * t830;
t823 = m(2) * t997 - mrSges(2,1) * t1019 - qJDD(1) * mrSges(2,2) + t825;
t1043 = t1013 * t823 + t1018 * t817;
t819 = t1008 * t832 + t852 * t1038 + t830 * t1039;
t1028 = -t1013 * t817 + t1018 * t823;
t893 = Ifges(7,5) * t930 + Ifges(7,6) * t929 + Ifges(7,3) * t939;
t895 = Ifges(7,1) * t930 + Ifges(7,4) * t929 + Ifges(7,5) * t939;
t862 = -mrSges(7,1) * t875 + mrSges(7,3) * t874 + Ifges(7,4) * t892 + Ifges(7,2) * t891 + Ifges(7,6) * t904 - t893 * t930 + t895 * t939;
t894 = Ifges(7,4) * t930 + Ifges(7,2) * t929 + Ifges(7,6) * t939;
t863 = mrSges(7,2) * t875 - mrSges(7,3) * t873 + Ifges(7,1) * t892 + Ifges(7,4) * t891 + Ifges(7,5) * t904 + t893 * t929 - t894 * t939;
t914 = Ifges(6,5) * t941 + Ifges(6,6) * t940 + Ifges(6,3) * t993;
t915 = Ifges(6,4) * t941 + Ifges(6,2) * t940 + Ifges(6,6) * t993;
t841 = mrSges(6,2) * t889 - mrSges(6,3) * t878 + Ifges(6,1) * t906 + Ifges(6,4) * t905 + Ifges(6,5) * t982 - pkin(11) * t858 - t1009 * t862 + t1014 * t863 + t914 * t940 - t915 * t993;
t1023 = mrSges(7,1) * t873 - mrSges(7,2) * t874 + Ifges(7,5) * t892 + Ifges(7,6) * t891 + Ifges(7,3) * t904 + t894 * t930 - t895 * t929;
t916 = Ifges(6,1) * t941 + Ifges(6,4) * t940 + Ifges(6,5) * t993;
t842 = -mrSges(6,1) * t889 + mrSges(6,3) * t879 + Ifges(6,4) * t906 + Ifges(6,2) * t905 + Ifges(6,6) * t982 - pkin(5) * t858 - t914 * t941 + t916 * t993 - t1023;
t933 = Ifges(5,5) * t964 + Ifges(5,6) * t963 + Ifges(5,3) * t993;
t935 = Ifges(5,1) * t964 + Ifges(5,4) * t963 + Ifges(5,5) * t993;
t826 = -mrSges(5,1) * t913 + mrSges(5,3) * t887 + Ifges(5,4) * t927 + Ifges(5,2) * t926 + Ifges(5,6) * t982 - pkin(4) * t856 + qJ(5) * t1032 + t1005 * t841 + t1007 * t842 - t964 * t933 + t993 * t935;
t934 = Ifges(5,4) * t964 + Ifges(5,2) * t963 + Ifges(5,6) * t993;
t834 = mrSges(5,2) * t913 - mrSges(5,3) * t886 + Ifges(5,1) * t927 + Ifges(5,4) * t926 + Ifges(5,5) * t982 - qJ(5) * t849 - t1005 * t842 + t1007 * t841 + t933 * t963 - t934 * t993;
t952 = Ifges(4,5) * t979 + Ifges(4,6) * t978 + Ifges(4,3) * t995;
t954 = Ifges(4,1) * t979 + Ifges(4,4) * t978 + Ifges(4,5) * t995;
t812 = -mrSges(4,1) * t946 + mrSges(4,3) * t922 + Ifges(4,4) * t959 + Ifges(4,2) * t958 + Ifges(4,6) * t983 - pkin(3) * t1024 + pkin(10) * t1030 + t1010 * t834 + t1015 * t826 - t979 * t952 + t995 * t954;
t953 = Ifges(4,4) * t979 + Ifges(4,2) * t978 + Ifges(4,6) * t995;
t813 = mrSges(4,2) * t946 - mrSges(4,3) * t921 + Ifges(4,1) * t959 + Ifges(4,4) * t958 + Ifges(4,5) * t983 - pkin(10) * t840 - t1010 * t826 + t1015 * t834 + t952 * t978 - t953 * t995;
t970 = Ifges(3,6) * t1001 + (Ifges(3,4) * t1012 + Ifges(3,2) * t1017) * t1041;
t971 = Ifges(3,5) * t1001 + (Ifges(3,1) * t1012 + Ifges(3,4) * t1017) * t1041;
t809 = Ifges(3,5) * t990 + Ifges(3,6) * t991 + Ifges(3,3) * t1000 + mrSges(3,1) * t960 - mrSges(3,2) * t961 + t1011 * t813 + t1016 * t812 + pkin(2) * t1021 + pkin(9) * t1029 + (t1012 * t970 - t1017 * t971) * t1041;
t969 = Ifges(3,3) * t1001 + (Ifges(3,5) * t1012 + Ifges(3,6) * t1017) * t1041;
t811 = mrSges(3,2) * t972 - mrSges(3,3) * t960 + Ifges(3,1) * t990 + Ifges(3,4) * t991 + Ifges(3,5) * t1000 - pkin(9) * t833 - t1001 * t970 - t1011 * t812 + t1016 * t813 + t969 * t1033;
t1022 = -mrSges(5,1) * t886 - mrSges(6,1) * t878 + mrSges(5,2) * t887 + mrSges(6,2) * t879 - Ifges(6,6) * t905 - pkin(4) * t849 - pkin(5) * t1026 - pkin(11) * t1031 - t1009 * t863 - t1014 * t862 + t940 * t916 + t963 * t935 - Ifges(6,5) * t906 - t941 * t915 - Ifges(5,6) * t926 - Ifges(5,5) * t927 - t964 * t934 + (-Ifges(6,3) - Ifges(5,3)) * t982;
t1020 = mrSges(4,1) * t921 - mrSges(4,2) * t922 + Ifges(4,5) * t959 + Ifges(4,6) * t958 + Ifges(4,3) * t983 + pkin(3) * t840 + t979 * t953 - t978 * t954 - t1022;
t815 = -mrSges(3,1) * t972 + mrSges(3,3) * t961 + Ifges(3,4) * t990 + Ifges(3,2) * t991 + Ifges(3,6) * t1000 - pkin(2) * t833 + t1001 * t971 - t969 * t1034 - t1020;
t1025 = mrSges(2,1) * t996 - mrSges(2,2) * t997 + Ifges(2,3) * qJDD(1) + pkin(1) * t820 + t1008 * t809 + t815 * t1038 + t811 * t1039 + t825 * t1045;
t807 = -mrSges(2,2) * g(3) - mrSges(2,3) * t996 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t1019 - t1012 * t815 + t1017 * t811 + (-t1006 * t819 - t1008 * t820) * pkin(8);
t806 = mrSges(2,1) * g(3) + mrSges(2,3) * t997 + Ifges(2,5) * t1019 + Ifges(2,6) * qJDD(1) - pkin(1) * t819 - t1006 * t809 + (pkin(8) * t825 + t1012 * t811 + t1017 * t815) * t1008;
t1 = [-m(1) * g(1) + t1028; -m(1) * g(2) + t1043; (-m(1) - m(2)) * g(3) + t819; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1043 - t1013 * t806 + t1018 * t807; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1028 + t1013 * t807 + t1018 * t806; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1025; t1025; t809; t1020; -t1022; t856; t1023;];
tauJB  = t1;
