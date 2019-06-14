% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 06:39
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRP10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:21:29
% EndTime: 2019-05-08 06:22:20
% DurationCPUTime: 34.50s
% Computational Cost: add. (581559->378), mult. (1234017->476), div. (0->0), fcn. (985988->12), ass. (0->159)
t1040 = Ifges(6,1) + Ifges(7,1);
t1029 = Ifges(6,4) - Ifges(7,5);
t1038 = Ifges(7,4) + Ifges(6,5);
t1039 = Ifges(6,2) + Ifges(7,3);
t1037 = Ifges(6,6) - Ifges(7,6);
t1036 = -Ifges(6,3) - Ifges(7,2);
t1033 = cos(qJ(5));
t1000 = cos(qJ(3));
t997 = sin(qJ(2));
t1021 = qJD(1) * t997;
t992 = sin(pkin(6));
t1016 = t992 * t1021;
t993 = cos(pkin(6));
t988 = t993 * qJD(1) + qJD(2);
t996 = sin(qJ(3));
t967 = t1000 * t1016 + t996 * t988;
t1001 = cos(qJ(2));
t1018 = qJD(1) * t1001;
t1015 = t992 * t1018;
t983 = qJD(3) - t1015;
t995 = sin(qJ(4));
t999 = cos(qJ(4));
t952 = -t995 * t967 + t999 * t983;
t953 = t999 * t967 + t995 * t983;
t994 = sin(qJ(5));
t926 = -t1033 * t952 + t994 * t953;
t927 = t1033 * t953 + t994 * t952;
t966 = t1000 * t988 - t996 * t1016;
t965 = qJD(4) - t966;
t963 = qJD(5) + t965;
t1035 = -t1029 * t927 - t1037 * t963 + t1039 * t926;
t1034 = -t1029 * t926 + t1038 * t963 + t1040 * t927;
t979 = (qJD(2) * t1021 - qJDD(1) * t1001) * t992;
t1032 = pkin(8) * t992;
t1031 = t993 * g(3);
t1030 = -mrSges(6,3) - mrSges(7,2);
t1028 = t992 * t997;
t1027 = t993 * t997;
t1002 = cos(qJ(1));
t1003 = qJD(1) ^ 2;
t1019 = t1001 * t993;
t998 = sin(qJ(1));
t984 = t998 * g(1) - t1002 * g(2);
t974 = qJDD(1) * pkin(1) + t1003 * t1032 + t984;
t985 = -t1002 * g(1) - t998 * g(2);
t975 = -t1003 * pkin(1) + qJDD(1) * t1032 + t985;
t1023 = t1001 * t975 + t974 * t1027;
t1022 = qJD(1) * t992;
t977 = (-pkin(2) * t1001 - pkin(9) * t997) * t1022;
t986 = t988 ^ 2;
t987 = t993 * qJDD(1) + qJDD(2);
t924 = -t986 * pkin(2) + t987 * pkin(9) + (-g(3) * t997 + t1018 * t977) * t992 + t1023;
t978 = (qJD(2) * t1018 + qJDD(1) * t997) * t992;
t925 = t979 * pkin(2) - t978 * pkin(9) - t1031 + (-t974 + (pkin(2) * t997 - pkin(9) * t1001) * t988 * qJD(1)) * t992;
t892 = t1000 * t924 + t996 * t925;
t949 = -t966 * pkin(3) - t967 * pkin(10);
t971 = qJDD(3) + t979;
t981 = t983 ^ 2;
t885 = -t981 * pkin(3) + t971 * pkin(10) + t966 * t949 + t892;
t1020 = t1001 * t992;
t946 = -g(3) * t1020 + t1019 * t974 - t997 * t975;
t923 = -t987 * pkin(2) - t986 * pkin(9) + t977 * t1016 - t946;
t944 = -t967 * qJD(3) + t1000 * t987 - t996 * t978;
t945 = t966 * qJD(3) + t1000 * t978 + t996 * t987;
t890 = (-t966 * t983 - t945) * pkin(10) + (t967 * t983 - t944) * pkin(3) + t923;
t871 = -t995 * t885 + t999 * t890;
t910 = t952 * qJD(4) + t999 * t945 + t995 * t971;
t942 = qJDD(4) - t944;
t868 = (t952 * t965 - t910) * pkin(11) + (t952 * t953 + t942) * pkin(4) + t871;
t872 = t999 * t885 + t995 * t890;
t909 = -t953 * qJD(4) - t995 * t945 + t999 * t971;
t932 = t965 * pkin(4) - t953 * pkin(11);
t951 = t952 ^ 2;
t870 = -t951 * pkin(4) + t909 * pkin(11) - t965 * t932 + t872;
t864 = t1033 * t870 + t994 * t868;
t902 = t926 * pkin(5) - t927 * qJ(6);
t937 = qJDD(5) + t942;
t961 = t963 ^ 2;
t860 = -t961 * pkin(5) + t937 * qJ(6) + 0.2e1 * qJD(6) * t963 - t926 * t902 + t864;
t914 = -t963 * mrSges(7,1) + t927 * mrSges(7,2);
t1017 = m(7) * t860 + t937 * mrSges(7,3) + t963 * t914;
t903 = t926 * mrSges(7,1) - t927 * mrSges(7,3);
t1024 = -t926 * mrSges(6,1) - t927 * mrSges(6,2) - t903;
t881 = t927 * qJD(5) - t1033 * t909 + t994 * t910;
t913 = t963 * mrSges(6,1) - t927 * mrSges(6,3);
t847 = m(6) * t864 - t937 * mrSges(6,2) + t1024 * t926 + t1030 * t881 - t963 * t913 + t1017;
t863 = t1033 * t868 - t994 * t870;
t861 = -t937 * pkin(5) - t961 * qJ(6) + t927 * t902 + qJDD(6) - t863;
t911 = -t926 * mrSges(7,2) + t963 * mrSges(7,3);
t1011 = -m(7) * t861 + t937 * mrSges(7,1) + t963 * t911;
t882 = -t926 * qJD(5) + t1033 * t910 + t994 * t909;
t912 = -t963 * mrSges(6,2) - t926 * mrSges(6,3);
t849 = m(6) * t863 + t937 * mrSges(6,1) + t1024 * t927 + t1030 * t882 + t963 * t912 + t1011;
t844 = t1033 * t849 + t994 * t847;
t928 = -t952 * mrSges(5,1) + t953 * mrSges(5,2);
t930 = -t965 * mrSges(5,2) + t952 * mrSges(5,3);
t840 = m(5) * t871 + t942 * mrSges(5,1) - t910 * mrSges(5,3) - t953 * t928 + t965 * t930 + t844;
t1012 = t1033 * t847 - t994 * t849;
t931 = t965 * mrSges(5,1) - t953 * mrSges(5,3);
t841 = m(5) * t872 - t942 * mrSges(5,2) + t909 * mrSges(5,3) + t952 * t928 - t965 * t931 + t1012;
t838 = -t995 * t840 + t999 * t841;
t948 = -t966 * mrSges(4,1) + t967 * mrSges(4,2);
t955 = t983 * mrSges(4,1) - t967 * mrSges(4,3);
t836 = m(4) * t892 - t971 * mrSges(4,2) + t944 * mrSges(4,3) + t966 * t948 - t983 * t955 + t838;
t891 = t1000 * t925 - t996 * t924;
t884 = -t971 * pkin(3) - t981 * pkin(10) + t967 * t949 - t891;
t873 = -t909 * pkin(4) - t951 * pkin(11) + t953 * t932 + t884;
t866 = -0.2e1 * qJD(6) * t927 + (t926 * t963 - t882) * qJ(6) + (t927 * t963 + t881) * pkin(5) + t873;
t857 = m(7) * t866 + t881 * mrSges(7,1) - t882 * mrSges(7,3) + t926 * t911 - t927 * t914;
t1008 = m(6) * t873 + t881 * mrSges(6,1) + t882 * mrSges(6,2) + t926 * t912 + t927 * t913 + t857;
t852 = -m(5) * t884 + t909 * mrSges(5,1) - t910 * mrSges(5,2) + t952 * t930 - t953 * t931 - t1008;
t954 = -t983 * mrSges(4,2) + t966 * mrSges(4,3);
t851 = m(4) * t891 + t971 * mrSges(4,1) - t945 * mrSges(4,3) - t967 * t948 + t983 * t954 + t852;
t1013 = t1000 * t836 - t996 * t851;
t947 = -g(3) * t1028 + t1023;
t972 = t988 * mrSges(3,1) - mrSges(3,3) * t1016;
t976 = (-mrSges(3,1) * t1001 + mrSges(3,2) * t997) * t1022;
t827 = m(3) * t947 - t987 * mrSges(3,2) - t979 * mrSges(3,3) + t1015 * t976 - t988 * t972 + t1013;
t830 = t1000 * t851 + t996 * t836;
t959 = -t992 * t974 - t1031;
t973 = -t988 * mrSges(3,2) + mrSges(3,3) * t1015;
t829 = m(3) * t959 + t979 * mrSges(3,1) + t978 * mrSges(3,2) + (-t1001 * t973 + t972 * t997) * t1022 + t830;
t837 = t999 * t840 + t995 * t841;
t1007 = -m(4) * t923 + t944 * mrSges(4,1) - t945 * mrSges(4,2) + t966 * t954 - t967 * t955 - t837;
t833 = m(3) * t946 + t987 * mrSges(3,1) - t978 * mrSges(3,3) - t1016 * t976 + t988 * t973 + t1007;
t815 = t833 * t1019 + t827 * t1027 - t992 * t829;
t812 = m(2) * t984 + qJDD(1) * mrSges(2,1) - t1003 * mrSges(2,2) + t815;
t820 = t1001 * t827 - t997 * t833;
t818 = m(2) * t985 - t1003 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t820;
t1026 = t1002 * t812 + t998 * t818;
t1025 = t1036 * t963 + t1037 * t926 - t1038 * t927;
t814 = t833 * t1020 + t827 * t1028 + t993 * t829;
t1014 = t1002 * t818 - t998 * t812;
t842 = -mrSges(6,1) * t873 - mrSges(7,1) * t866 + mrSges(7,2) * t860 + mrSges(6,3) * t864 - pkin(5) * t857 + t1025 * t927 + t1029 * t882 + t1034 * t963 + t1037 * t937 - t1039 * t881;
t843 = mrSges(6,2) * t873 + mrSges(7,2) * t861 - mrSges(6,3) * t863 - mrSges(7,3) * t866 - qJ(6) * t857 + t1025 * t926 - t1029 * t881 + t1035 * t963 + t1038 * t937 + t1040 * t882;
t915 = Ifges(5,5) * t953 + Ifges(5,6) * t952 + Ifges(5,3) * t965;
t917 = Ifges(5,1) * t953 + Ifges(5,4) * t952 + Ifges(5,5) * t965;
t822 = -mrSges(5,1) * t884 + mrSges(5,3) * t872 + Ifges(5,4) * t910 + Ifges(5,2) * t909 + Ifges(5,6) * t942 - pkin(4) * t1008 + pkin(11) * t1012 + t1033 * t842 + t994 * t843 - t953 * t915 + t965 * t917;
t916 = Ifges(5,4) * t953 + Ifges(5,2) * t952 + Ifges(5,6) * t965;
t823 = mrSges(5,2) * t884 - mrSges(5,3) * t871 + Ifges(5,1) * t910 + Ifges(5,4) * t909 + Ifges(5,5) * t942 - pkin(11) * t844 + t1033 * t843 - t994 * t842 + t952 * t915 - t965 * t916;
t938 = Ifges(4,5) * t967 + Ifges(4,6) * t966 + Ifges(4,3) * t983;
t939 = Ifges(4,4) * t967 + Ifges(4,2) * t966 + Ifges(4,6) * t983;
t810 = mrSges(4,2) * t923 - mrSges(4,3) * t891 + Ifges(4,1) * t945 + Ifges(4,4) * t944 + Ifges(4,5) * t971 - pkin(10) * t837 - t995 * t822 + t999 * t823 + t966 * t938 - t983 * t939;
t856 = t882 * mrSges(7,2) + t927 * t903 - t1011;
t1006 = -mrSges(6,1) * t863 + mrSges(7,1) * t861 + mrSges(6,2) * t864 - mrSges(7,3) * t860 + pkin(5) * t856 - qJ(6) * t1017 + t1036 * t937 + t1035 * t927 + (qJ(6) * t903 - t1034) * t926 - t1038 * t882 + (qJ(6) * mrSges(7,2) + t1037) * t881;
t1004 = mrSges(5,1) * t871 - mrSges(5,2) * t872 + Ifges(5,5) * t910 + Ifges(5,6) * t909 + Ifges(5,3) * t942 + pkin(4) * t844 + t953 * t916 - t952 * t917 - t1006;
t940 = Ifges(4,1) * t967 + Ifges(4,4) * t966 + Ifges(4,5) * t983;
t821 = -mrSges(4,1) * t923 + mrSges(4,3) * t892 + Ifges(4,4) * t945 + Ifges(4,2) * t944 + Ifges(4,6) * t971 - pkin(3) * t837 - t967 * t938 + t983 * t940 - t1004;
t957 = Ifges(3,6) * t988 + (Ifges(3,4) * t997 + Ifges(3,2) * t1001) * t1022;
t958 = Ifges(3,5) * t988 + (Ifges(3,1) * t997 + Ifges(3,4) * t1001) * t1022;
t805 = Ifges(3,5) * t978 - Ifges(3,6) * t979 + Ifges(3,3) * t987 + mrSges(3,1) * t946 - mrSges(3,2) * t947 + t996 * t810 + t1000 * t821 + pkin(2) * t1007 + pkin(9) * t1013 + (-t1001 * t958 + t957 * t997) * t1022;
t956 = Ifges(3,3) * t988 + (Ifges(3,5) * t997 + Ifges(3,6) * t1001) * t1022;
t807 = mrSges(3,2) * t959 - mrSges(3,3) * t946 + Ifges(3,1) * t978 - Ifges(3,4) * t979 + Ifges(3,5) * t987 - pkin(9) * t830 + t1000 * t810 + t1015 * t956 - t996 * t821 - t988 * t957;
t1005 = mrSges(4,1) * t891 - mrSges(4,2) * t892 + Ifges(4,5) * t945 + Ifges(4,6) * t944 + Ifges(4,3) * t971 + pkin(3) * t852 + pkin(10) * t838 + t999 * t822 + t995 * t823 + t967 * t939 - t966 * t940;
t809 = -mrSges(3,1) * t959 + mrSges(3,3) * t947 + Ifges(3,4) * t978 - Ifges(3,2) * t979 + Ifges(3,6) * t987 - pkin(2) * t830 - t1016 * t956 + t988 * t958 - t1005;
t1009 = mrSges(2,1) * t984 - mrSges(2,2) * t985 + Ifges(2,3) * qJDD(1) + pkin(1) * t815 + t809 * t1020 + t807 * t1028 + t820 * t1032 + t993 * t805;
t803 = -mrSges(2,2) * g(3) - mrSges(2,3) * t984 + Ifges(2,5) * qJDD(1) - t1003 * Ifges(2,6) + t1001 * t807 - t997 * t809 + (-t814 * t992 - t815 * t993) * pkin(8);
t802 = mrSges(2,1) * g(3) + mrSges(2,3) * t985 + t1003 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t814 - t992 * t805 + (pkin(8) * t820 + t1001 * t809 + t807 * t997) * t993;
t1 = [-m(1) * g(1) + t1014; -m(1) * g(2) + t1026; (-m(1) - m(2)) * g(3) + t814; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1026 + t1002 * t803 - t998 * t802; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1014 + t1002 * t802 + t998 * t803; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1009; t1009; t805; t1005; t1004; -t1006; t856;];
tauJB  = t1;
