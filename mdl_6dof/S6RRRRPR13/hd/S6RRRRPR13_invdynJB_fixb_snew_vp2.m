% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-08 01:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR13_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR13_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 01:21:44
% EndTime: 2019-05-08 01:22:24
% DurationCPUTime: 24.94s
% Computational Cost: add. (412567->381), mult. (872809->476), div. (0->0), fcn. (689120->12), ass. (0->161)
t1058 = Ifges(5,1) + Ifges(6,1);
t1050 = Ifges(5,4) - Ifges(6,5);
t1049 = Ifges(5,5) + Ifges(6,4);
t1057 = -Ifges(5,2) - Ifges(6,3);
t1048 = Ifges(5,6) - Ifges(6,6);
t1056 = Ifges(5,3) + Ifges(6,2);
t1004 = sin(pkin(6));
t1009 = sin(qJ(2));
t1013 = cos(qJ(2));
t989 = (qJD(1) * qJD(2) * t1009 - qJDD(1) * t1013) * t1004;
t1008 = sin(qJ(3));
t1012 = cos(qJ(3));
t1036 = qJD(1) * t1013;
t1005 = cos(pkin(6));
t1033 = t1005 * t1009;
t1015 = qJD(1) ^ 2;
t1046 = pkin(8) * t1004;
t1010 = sin(qJ(1));
t1014 = cos(qJ(1));
t996 = t1010 * g(1) - g(2) * t1014;
t984 = qJDD(1) * pkin(1) + t1015 * t1046 + t996;
t997 = -g(1) * t1014 - g(2) * t1010;
t985 = -pkin(1) * t1015 + qJDD(1) * t1046 + t997;
t1038 = t1013 * t985 + t984 * t1033;
t1037 = qJD(1) * t1004;
t987 = (-pkin(2) * t1013 - pkin(9) * t1009) * t1037;
t1000 = qJD(1) * t1005 + qJD(2);
t998 = t1000 ^ 2;
t999 = qJDD(1) * t1005 + qJDD(2);
t929 = -t998 * pkin(2) + t999 * pkin(9) + (-g(3) * t1009 + t1036 * t987) * t1004 + t1038;
t1045 = t1005 * g(3);
t988 = (qJD(2) * t1036 + qJDD(1) * t1009) * t1004;
t930 = t989 * pkin(2) - t988 * pkin(9) - t1045 + (-t984 + (pkin(2) * t1009 - pkin(9) * t1013) * t1000 * qJD(1)) * t1004;
t900 = -t1008 * t929 + t1012 * t930;
t1035 = t1004 * t1009;
t1031 = qJD(1) * t1035;
t977 = t1012 * t1000 - t1008 * t1031;
t978 = t1000 * t1008 + t1012 * t1031;
t960 = -pkin(3) * t977 - pkin(10) * t978;
t981 = qJDD(3) + t989;
t1034 = t1004 * t1013;
t1030 = qJD(1) * t1034;
t995 = qJD(3) - t1030;
t994 = t995 ^ 2;
t1022 = t981 * pkin(3) + t994 * pkin(10) - t978 * t960 + t900;
t1007 = sin(qJ(4));
t1052 = cos(qJ(4));
t962 = t1007 * t978 - t1052 * t995;
t975 = qJD(4) - t977;
t1044 = t962 * t975;
t956 = qJD(3) * t977 + t1008 * t999 + t1012 * t988;
t911 = -t962 * qJD(4) + t1007 * t981 + t1052 * t956;
t1055 = (-t911 + t1044) * qJ(5) - t1022;
t1006 = sin(qJ(6));
t1011 = cos(qJ(6));
t901 = t1008 * t930 + t1012 * t929;
t896 = -pkin(3) * t994 + pkin(10) * t981 + t960 * t977 + t901;
t1032 = t1005 * t1013;
t957 = -g(3) * t1034 - t1009 * t985 + t1032 * t984;
t928 = -t999 * pkin(2) - t998 * pkin(9) + t987 * t1031 - t957;
t955 = -qJD(3) * t978 - t1008 * t988 + t1012 * t999;
t899 = (-t977 * t995 - t956) * pkin(10) + (t978 * t995 - t955) * pkin(3) + t928;
t885 = -t1007 * t896 + t1052 * t899;
t963 = t1007 * t995 + t1052 * t978;
t935 = pkin(4) * t962 - qJ(5) * t963;
t953 = qJDD(4) - t955;
t974 = t975 ^ 2;
t883 = -t953 * pkin(4) - t974 * qJ(5) + t963 * t935 + qJDD(5) - t885;
t877 = (-t911 - t1044) * pkin(11) + (t962 * t963 - t953) * pkin(5) + t883;
t1053 = 2 * qJD(5);
t886 = t1007 * t899 + t1052 * t896;
t882 = -pkin(4) * t974 + t953 * qJ(5) + t1053 * t975 - t962 * t935 + t886;
t910 = qJD(4) * t963 + t1007 * t956 - t1052 * t981;
t943 = -pkin(5) * t975 - pkin(11) * t963;
t961 = t962 ^ 2;
t878 = -pkin(5) * t961 + pkin(11) * t910 + t943 * t975 + t882;
t875 = -t1006 * t878 + t1011 * t877;
t876 = t1006 * t877 + t1011 * t878;
t934 = t1006 * t962 + t1011 * t963;
t891 = -qJD(6) * t934 - t1006 * t911 + t1011 * t910;
t933 = -t1006 * t963 + t1011 * t962;
t892 = qJD(6) * t933 + t1006 * t910 + t1011 * t911;
t973 = qJD(6) - t975;
t903 = Ifges(7,4) * t934 + Ifges(7,2) * t933 + Ifges(7,6) * t973;
t904 = Ifges(7,1) * t934 + Ifges(7,4) * t933 + Ifges(7,5) * t973;
t948 = qJDD(6) - t953;
t1019 = -mrSges(7,1) * t875 + mrSges(7,2) * t876 - Ifges(7,5) * t892 - Ifges(7,6) * t891 - Ifges(7,3) * t948 - t934 * t903 + t933 * t904;
t907 = -mrSges(7,1) * t933 + mrSges(7,2) * t934;
t914 = -mrSges(7,2) * t973 + mrSges(7,3) * t933;
t872 = m(7) * t875 + mrSges(7,1) * t948 - mrSges(7,3) * t892 - t907 * t934 + t914 * t973;
t915 = mrSges(7,1) * t973 - mrSges(7,3) * t934;
t873 = m(7) * t876 - mrSges(7,2) * t948 + mrSges(7,3) * t891 + t907 * t933 - t915 * t973;
t1028 = -t1006 * t872 + t1011 * t873;
t942 = -mrSges(6,1) * t975 + mrSges(6,2) * t963;
t1024 = m(6) * t882 + t953 * mrSges(6,3) + t975 * t942 + t1028;
t1040 = t1049 * t975 - t1050 * t962 + t1058 * t963;
t1041 = t1048 * t975 + t1050 * t963 + t1057 * t962;
t864 = t1006 * t873 + t1011 * t872;
t939 = -mrSges(6,2) * t962 + mrSges(6,3) * t975;
t1020 = -m(6) * t883 + t953 * mrSges(6,1) + t975 * t939 - t864;
t936 = mrSges(6,1) * t962 - mrSges(6,3) * t963;
t863 = t911 * mrSges(6,2) + t963 * t936 - t1020;
t1054 = t1040 * t962 + t1041 * t963 + t1056 * t953 - t1048 * t910 + t1049 * t911 + mrSges(5,1) * t885 - mrSges(6,1) * t883 - mrSges(5,2) * t886 + mrSges(6,3) * t882 - pkin(4) * t863 - pkin(5) * t864 + qJ(5) * (-t910 * mrSges(6,2) - t962 * t936 + t1024) + t1019;
t1051 = -mrSges(5,3) - mrSges(6,2);
t1039 = -mrSges(5,1) * t962 - mrSges(5,2) * t963 - t936;
t941 = mrSges(5,1) * t975 - mrSges(5,3) * t963;
t860 = m(5) * t886 - t953 * mrSges(5,2) + t1039 * t962 + t1051 * t910 - t975 * t941 + t1024;
t940 = -mrSges(5,2) * t975 - mrSges(5,3) * t962;
t861 = m(5) * t885 + t953 * mrSges(5,1) + t1039 * t963 + t1051 * t911 + t975 * t940 + t1020;
t858 = -t1007 * t861 + t1052 * t860;
t959 = -mrSges(4,1) * t977 + mrSges(4,2) * t978;
t965 = mrSges(4,1) * t995 - mrSges(4,3) * t978;
t856 = m(4) * t901 - mrSges(4,2) * t981 + mrSges(4,3) * t955 + t959 * t977 - t965 * t995 + t858;
t880 = -t961 * pkin(11) + (-pkin(4) - pkin(5)) * t910 + (-pkin(4) * t975 + t1053 + t943) * t963 - t1055;
t1025 = -m(7) * t880 + t891 * mrSges(7,1) - t892 * mrSges(7,2) + t933 * t914 - t934 * t915;
t884 = -0.2e1 * qJD(5) * t963 + (t963 * t975 + t910) * pkin(4) + t1055;
t870 = m(6) * t884 + mrSges(6,1) * t910 - t911 * mrSges(6,3) + t939 * t962 - t963 * t942 + t1025;
t869 = m(5) * t1022 - t910 * mrSges(5,1) - mrSges(5,2) * t911 - t962 * t940 - t941 * t963 - t870;
t964 = -mrSges(4,2) * t995 + mrSges(4,3) * t977;
t868 = m(4) * t900 + mrSges(4,1) * t981 - mrSges(4,3) * t956 - t959 * t978 + t964 * t995 + t869;
t1027 = -t1008 * t868 + t1012 * t856;
t958 = -g(3) * t1035 + t1038;
t982 = mrSges(3,1) * t1000 - mrSges(3,3) * t1031;
t986 = (-mrSges(3,1) * t1013 + mrSges(3,2) * t1009) * t1037;
t847 = m(3) * t958 - mrSges(3,2) * t999 - mrSges(3,3) * t989 - t1000 * t982 + t1030 * t986 + t1027;
t850 = t1008 * t856 + t1012 * t868;
t969 = -t1004 * t984 - t1045;
t983 = -mrSges(3,2) * t1000 + mrSges(3,3) * t1030;
t849 = m(3) * t969 + t989 * mrSges(3,1) + t988 * mrSges(3,2) + (t1009 * t982 - t1013 * t983) * t1037 + t850;
t857 = t1007 * t860 + t1052 * t861;
t1018 = -m(4) * t928 + t955 * mrSges(4,1) - t956 * mrSges(4,2) + t977 * t964 - t978 * t965 - t857;
t853 = m(3) * t957 + t999 * mrSges(3,1) - t988 * mrSges(3,3) + t1000 * t983 - t1031 * t986 + t1018;
t835 = -t1004 * t849 + t853 * t1032 + t847 * t1033;
t832 = m(2) * t996 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1015 + t835;
t841 = -t1009 * t853 + t1013 * t847;
t839 = m(2) * t997 - mrSges(2,1) * t1015 - qJDD(1) * mrSges(2,2) + t841;
t1043 = t1010 * t839 + t1014 * t832;
t1042 = t1048 * t962 - t1049 * t963 - t1056 * t975;
t834 = t1005 * t849 + t853 * t1034 + t847 * t1035;
t1026 = -t1010 * t832 + t1014 * t839;
t902 = Ifges(7,5) * t934 + Ifges(7,6) * t933 + Ifges(7,3) * t973;
t865 = -mrSges(7,1) * t880 + mrSges(7,3) * t876 + Ifges(7,4) * t892 + Ifges(7,2) * t891 + Ifges(7,6) * t948 - t902 * t934 + t904 * t973;
t866 = mrSges(7,2) * t880 - mrSges(7,3) * t875 + Ifges(7,1) * t892 + Ifges(7,4) * t891 + Ifges(7,5) * t948 + t902 * t933 - t903 * t973;
t842 = mrSges(5,1) * t1022 - mrSges(6,1) * t884 + mrSges(6,2) * t882 + mrSges(5,3) * t886 - pkin(4) * t870 - pkin(5) * t1025 - pkin(11) * t1028 - t1006 * t866 - t1011 * t865 + t1040 * t975 + t1042 * t963 + t1048 * t953 + t1050 * t911 + t1057 * t910;
t843 = -mrSges(5,2) * t1022 + mrSges(6,2) * t883 - mrSges(5,3) * t885 - mrSges(6,3) * t884 - pkin(11) * t864 - qJ(5) * t870 - t1006 * t865 + t1011 * t866 - t1041 * t975 + t1042 * t962 + t1049 * t953 - t1050 * t910 + t1058 * t911;
t949 = Ifges(4,5) * t978 + Ifges(4,6) * t977 + Ifges(4,3) * t995;
t950 = Ifges(4,4) * t978 + Ifges(4,2) * t977 + Ifges(4,6) * t995;
t830 = mrSges(4,2) * t928 - mrSges(4,3) * t900 + Ifges(4,1) * t956 + Ifges(4,4) * t955 + Ifges(4,5) * t981 - pkin(10) * t857 - t1007 * t842 + t1052 * t843 + t977 * t949 - t995 * t950;
t951 = Ifges(4,1) * t978 + Ifges(4,4) * t977 + Ifges(4,5) * t995;
t836 = -mrSges(4,1) * t928 + mrSges(4,3) * t901 + Ifges(4,4) * t956 + Ifges(4,2) * t955 + Ifges(4,6) * t981 - pkin(3) * t857 - t978 * t949 + t995 * t951 - t1054;
t967 = Ifges(3,6) * t1000 + (Ifges(3,4) * t1009 + Ifges(3,2) * t1013) * t1037;
t968 = Ifges(3,5) * t1000 + (Ifges(3,1) * t1009 + Ifges(3,4) * t1013) * t1037;
t825 = Ifges(3,5) * t988 - Ifges(3,6) * t989 + Ifges(3,3) * t999 + mrSges(3,1) * t957 - mrSges(3,2) * t958 + t1008 * t830 + t1012 * t836 + pkin(2) * t1018 + pkin(9) * t1027 + (t1009 * t967 - t1013 * t968) * t1037;
t966 = Ifges(3,3) * t1000 + (Ifges(3,5) * t1009 + Ifges(3,6) * t1013) * t1037;
t827 = mrSges(3,2) * t969 - mrSges(3,3) * t957 + Ifges(3,1) * t988 - Ifges(3,4) * t989 + Ifges(3,5) * t999 - pkin(9) * t850 - t1000 * t967 - t1008 * t836 + t1012 * t830 + t1030 * t966;
t1016 = mrSges(4,1) * t900 - mrSges(4,2) * t901 + Ifges(4,5) * t956 + Ifges(4,6) * t955 + Ifges(4,3) * t981 + pkin(3) * t869 + pkin(10) * t858 + t1007 * t843 + t1052 * t842 + t978 * t950 - t977 * t951;
t829 = -mrSges(3,1) * t969 + mrSges(3,3) * t958 + Ifges(3,4) * t988 - Ifges(3,2) * t989 + Ifges(3,6) * t999 - pkin(2) * t850 + t1000 * t968 - t1031 * t966 - t1016;
t1021 = mrSges(2,1) * t996 - mrSges(2,2) * t997 + Ifges(2,3) * qJDD(1) + pkin(1) * t835 + t1005 * t825 + t829 * t1034 + t827 * t1035 + t841 * t1046;
t823 = -mrSges(2,2) * g(3) - mrSges(2,3) * t996 + Ifges(2,5) * qJDD(1) - t1015 * Ifges(2,6) - t1009 * t829 + t1013 * t827 + (-t1004 * t834 - t1005 * t835) * pkin(8);
t822 = mrSges(2,1) * g(3) + mrSges(2,3) * t997 + t1015 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t834 - t1004 * t825 + (pkin(8) * t841 + t1009 * t827 + t1013 * t829) * t1005;
t1 = [-m(1) * g(1) + t1026; -m(1) * g(2) + t1043; (-m(1) - m(2)) * g(3) + t834; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1043 - t1010 * t822 + t1014 * t823; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1026 + t1010 * t823 + t1014 * t822; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1021; t1021; t825; t1016; t1054; t863; -t1019;];
tauJB  = t1;
