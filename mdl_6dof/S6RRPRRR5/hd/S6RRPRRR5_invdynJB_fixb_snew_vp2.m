% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRR5
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
% Datum: 2019-05-06 21:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 21:19:55
% EndTime: 2019-05-06 21:20:59
% DurationCPUTime: 60.46s
% Computational Cost: add. (940615->400), mult. (2459925->519), div. (0->0), fcn. (1972947->14), ass. (0->167)
t1037 = -2 * qJD(3);
t1002 = sin(qJ(2));
t1007 = cos(qJ(2));
t996 = sin(pkin(6));
t1034 = qJD(1) * t996;
t995 = sin(pkin(12));
t997 = cos(pkin(12));
t974 = (t1002 * t995 - t1007 * t997) * t1034;
t998 = cos(pkin(6));
t1029 = t1007 * t998;
t1009 = qJD(1) ^ 2;
t1036 = pkin(8) * t996;
t1003 = sin(qJ(1));
t1008 = cos(qJ(1));
t986 = t1003 * g(1) - g(2) * t1008;
t980 = qJDD(1) * pkin(1) + t1009 * t1036 + t986;
t987 = -g(1) * t1008 - g(2) * t1003;
t981 = -pkin(1) * t1009 + qJDD(1) * t1036 + t987;
t1018 = -t1002 * t981 + t980 * t1029;
t1028 = t1009 * t996 ^ 2;
t1026 = qJD(1) * t1007;
t983 = (qJD(2) * t1026 + qJDD(1) * t1002) * t996;
t989 = qJDD(1) * t998 + qJDD(2);
t990 = qJD(1) * t998 + qJD(2);
t917 = t989 * pkin(2) - t983 * qJ(3) + (pkin(2) * t1002 * t1028 + (qJ(3) * qJD(1) * t990 - g(3)) * t996) * t1007 + t1018;
t1025 = t1007 ^ 2 * t1028;
t1031 = t1002 * t998;
t1032 = t1002 * t996;
t947 = -g(3) * t1032 + t1007 * t981 + t980 * t1031;
t1027 = qJD(1) * t1002;
t1024 = t996 * t1027;
t977 = pkin(2) * t990 - qJ(3) * t1024;
t984 = (-qJD(2) * t1027 + qJDD(1) * t1007) * t996;
t921 = -pkin(2) * t1025 + qJ(3) * t984 - t977 * t990 + t947;
t975 = (t1002 * t997 + t1007 * t995) * t1034;
t897 = t1037 * t975 + t997 * t917 - t995 * t921;
t1001 = sin(qJ(4));
t1006 = cos(qJ(4));
t1000 = sin(qJ(5));
t1005 = cos(qJ(5));
t1004 = cos(qJ(6));
t898 = t1037 * t974 + t995 * t917 + t997 * t921;
t949 = pkin(3) * t974 - pkin(9) * t975;
t988 = t990 ^ 2;
t891 = -pkin(3) * t988 + pkin(9) * t989 - t949 * t974 + t898;
t965 = -t998 * g(3) - t996 * t980;
t933 = -t984 * pkin(2) - qJ(3) * t1025 + t977 * t1024 + qJDD(3) + t965;
t953 = -t983 * t995 + t984 * t997;
t954 = t983 * t997 + t984 * t995;
t900 = (t974 * t990 - t954) * pkin(9) + (t975 * t990 - t953) * pkin(3) + t933;
t886 = t1001 * t900 + t1006 * t891;
t958 = -t1001 * t975 + t1006 * t990;
t959 = t1001 * t990 + t1006 * t975;
t935 = -pkin(4) * t958 - pkin(10) * t959;
t952 = qJDD(4) - t953;
t973 = qJD(4) + t974;
t972 = t973 ^ 2;
t876 = -pkin(4) * t972 + pkin(10) * t952 + t935 * t958 + t886;
t890 = -t989 * pkin(3) - t988 * pkin(9) + t975 * t949 - t897;
t930 = -t959 * qJD(4) - t1001 * t954 + t1006 * t989;
t931 = qJD(4) * t958 + t1001 * t989 + t1006 * t954;
t879 = (-t958 * t973 - t931) * pkin(10) + (t959 * t973 - t930) * pkin(4) + t890;
t871 = -t1000 * t876 + t1005 * t879;
t938 = -t1000 * t959 + t1005 * t973;
t903 = qJD(5) * t938 + t1000 * t952 + t1005 * t931;
t929 = qJDD(5) - t930;
t939 = t1000 * t973 + t1005 * t959;
t957 = qJD(5) - t958;
t869 = (t938 * t957 - t903) * pkin(11) + (t938 * t939 + t929) * pkin(5) + t871;
t872 = t1000 * t879 + t1005 * t876;
t902 = -qJD(5) * t939 - t1000 * t931 + t1005 * t952;
t920 = pkin(5) * t957 - pkin(11) * t939;
t937 = t938 ^ 2;
t870 = -pkin(5) * t937 + pkin(11) * t902 - t920 * t957 + t872;
t999 = sin(qJ(6));
t867 = t1004 * t869 - t870 * t999;
t910 = t1004 * t938 - t939 * t999;
t884 = qJD(6) * t910 + t1004 * t903 + t902 * t999;
t911 = t1004 * t939 + t938 * t999;
t896 = -mrSges(7,1) * t910 + mrSges(7,2) * t911;
t955 = qJD(6) + t957;
t904 = -mrSges(7,2) * t955 + mrSges(7,3) * t910;
t927 = qJDD(6) + t929;
t863 = m(7) * t867 + mrSges(7,1) * t927 - mrSges(7,3) * t884 - t896 * t911 + t904 * t955;
t868 = t1004 * t870 + t869 * t999;
t883 = -qJD(6) * t911 + t1004 * t902 - t903 * t999;
t905 = mrSges(7,1) * t955 - mrSges(7,3) * t911;
t864 = m(7) * t868 - mrSges(7,2) * t927 + mrSges(7,3) * t883 + t896 * t910 - t905 * t955;
t855 = t1004 * t863 + t999 * t864;
t912 = -mrSges(6,1) * t938 + mrSges(6,2) * t939;
t918 = -mrSges(6,2) * t957 + mrSges(6,3) * t938;
t853 = m(6) * t871 + mrSges(6,1) * t929 - mrSges(6,3) * t903 - t912 * t939 + t918 * t957 + t855;
t1021 = t1004 * t864 - t863 * t999;
t919 = mrSges(6,1) * t957 - mrSges(6,3) * t939;
t854 = m(6) * t872 - mrSges(6,2) * t929 + mrSges(6,3) * t902 + t912 * t938 - t919 * t957 + t1021;
t851 = -t1000 * t853 + t1005 * t854;
t934 = -mrSges(5,1) * t958 + mrSges(5,2) * t959;
t941 = mrSges(5,1) * t973 - mrSges(5,3) * t959;
t849 = m(5) * t886 - mrSges(5,2) * t952 + mrSges(5,3) * t930 + t934 * t958 - t941 * t973 + t851;
t885 = -t1001 * t891 + t1006 * t900;
t875 = -pkin(4) * t952 - pkin(10) * t972 + t959 * t935 - t885;
t873 = -pkin(5) * t902 - pkin(11) * t937 + t920 * t939 + t875;
t1015 = m(7) * t873 - t883 * mrSges(7,1) + mrSges(7,2) * t884 - t910 * t904 + t905 * t911;
t865 = -m(6) * t875 + t902 * mrSges(6,1) - mrSges(6,2) * t903 + t938 * t918 - t919 * t939 - t1015;
t940 = -mrSges(5,2) * t973 + mrSges(5,3) * t958;
t859 = m(5) * t885 + mrSges(5,1) * t952 - mrSges(5,3) * t931 - t934 * t959 + t940 * t973 + t865;
t1019 = -t1001 * t859 + t1006 * t849;
t948 = mrSges(4,1) * t974 + mrSges(4,2) * t975;
t961 = mrSges(4,1) * t990 - mrSges(4,3) * t975;
t838 = m(4) * t898 - mrSges(4,2) * t989 + mrSges(4,3) * t953 - t948 * t974 - t961 * t990 + t1019;
t850 = t1000 * t854 + t1005 * t853;
t1012 = -m(5) * t890 + t930 * mrSges(5,1) - mrSges(5,2) * t931 + t958 * t940 - t941 * t959 - t850;
t960 = -mrSges(4,2) * t990 - mrSges(4,3) * t974;
t846 = m(4) * t897 + mrSges(4,1) * t989 - mrSges(4,3) * t954 - t948 * t975 + t960 * t990 + t1012;
t834 = t995 * t838 + t997 * t846;
t1030 = t1007 * t996;
t946 = -g(3) * t1030 + t1018;
t1023 = t996 * t1026;
t979 = -mrSges(3,2) * t990 + mrSges(3,3) * t1023;
t982 = (-mrSges(3,1) * t1007 + mrSges(3,2) * t1002) * t1034;
t832 = m(3) * t946 + mrSges(3,1) * t989 - mrSges(3,3) * t983 - t982 * t1024 + t979 * t990 + t834;
t1022 = t997 * t838 - t846 * t995;
t978 = mrSges(3,1) * t990 - mrSges(3,3) * t1024;
t833 = m(3) * t947 - mrSges(3,2) * t989 + mrSges(3,3) * t984 + t982 * t1023 - t978 * t990 + t1022;
t842 = t1001 * t849 + t1006 * t859;
t841 = m(4) * t933 - t953 * mrSges(4,1) + t954 * mrSges(4,2) + t974 * t960 + t975 * t961 + t842;
t840 = m(3) * t965 - t984 * mrSges(3,1) + t983 * mrSges(3,2) + (t1002 * t978 - t1007 * t979) * t1034 + t841;
t819 = t832 * t1029 + t833 * t1031 - t840 * t996;
t816 = m(2) * t986 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1009 + t819;
t825 = -t1002 * t832 + t1007 * t833;
t823 = m(2) * t987 - mrSges(2,1) * t1009 - qJDD(1) * mrSges(2,2) + t825;
t1035 = t1003 * t823 + t1008 * t816;
t818 = t832 * t1030 + t833 * t1032 + t998 * t840;
t1017 = -t1003 * t816 + t1008 * t823;
t892 = Ifges(7,5) * t911 + Ifges(7,6) * t910 + Ifges(7,3) * t955;
t894 = Ifges(7,1) * t911 + Ifges(7,4) * t910 + Ifges(7,5) * t955;
t856 = -mrSges(7,1) * t873 + mrSges(7,3) * t868 + Ifges(7,4) * t884 + Ifges(7,2) * t883 + Ifges(7,6) * t927 - t892 * t911 + t894 * t955;
t893 = Ifges(7,4) * t911 + Ifges(7,2) * t910 + Ifges(7,6) * t955;
t857 = mrSges(7,2) * t873 - mrSges(7,3) * t867 + Ifges(7,1) * t884 + Ifges(7,4) * t883 + Ifges(7,5) * t927 + t892 * t910 - t893 * t955;
t906 = Ifges(6,5) * t939 + Ifges(6,6) * t938 + Ifges(6,3) * t957;
t908 = Ifges(6,1) * t939 + Ifges(6,4) * t938 + Ifges(6,5) * t957;
t843 = -mrSges(6,1) * t875 + mrSges(6,3) * t872 + Ifges(6,4) * t903 + Ifges(6,2) * t902 + Ifges(6,6) * t929 - pkin(5) * t1015 + pkin(11) * t1021 + t1004 * t856 + t999 * t857 - t939 * t906 + t957 * t908;
t907 = Ifges(6,4) * t939 + Ifges(6,2) * t938 + Ifges(6,6) * t957;
t844 = mrSges(6,2) * t875 - mrSges(6,3) * t871 + Ifges(6,1) * t903 + Ifges(6,4) * t902 + Ifges(6,5) * t929 - pkin(11) * t855 + t1004 * t857 - t856 * t999 + t906 * t938 - t907 * t957;
t922 = Ifges(5,5) * t959 + Ifges(5,6) * t958 + Ifges(5,3) * t973;
t923 = Ifges(5,4) * t959 + Ifges(5,2) * t958 + Ifges(5,6) * t973;
t826 = mrSges(5,2) * t890 - mrSges(5,3) * t885 + Ifges(5,1) * t931 + Ifges(5,4) * t930 + Ifges(5,5) * t952 - pkin(10) * t850 - t1000 * t843 + t1005 * t844 + t922 * t958 - t923 * t973;
t1013 = -mrSges(7,1) * t867 + mrSges(7,2) * t868 - Ifges(7,5) * t884 - Ifges(7,6) * t883 - Ifges(7,3) * t927 - t911 * t893 + t910 * t894;
t1010 = mrSges(6,1) * t871 - mrSges(6,2) * t872 + Ifges(6,5) * t903 + Ifges(6,6) * t902 + Ifges(6,3) * t929 + pkin(5) * t855 + t939 * t907 - t938 * t908 - t1013;
t924 = Ifges(5,1) * t959 + Ifges(5,4) * t958 + Ifges(5,5) * t973;
t835 = -mrSges(5,1) * t890 + mrSges(5,3) * t886 + Ifges(5,4) * t931 + Ifges(5,2) * t930 + Ifges(5,6) * t952 - pkin(4) * t850 - t959 * t922 + t973 * t924 - t1010;
t942 = Ifges(4,5) * t975 - Ifges(4,6) * t974 + Ifges(4,3) * t990;
t943 = Ifges(4,4) * t975 - Ifges(4,2) * t974 + Ifges(4,6) * t990;
t814 = mrSges(4,2) * t933 - mrSges(4,3) * t897 + Ifges(4,1) * t954 + Ifges(4,4) * t953 + Ifges(4,5) * t989 - pkin(9) * t842 - t1001 * t835 + t1006 * t826 - t942 * t974 - t943 * t990;
t1011 = mrSges(5,1) * t885 - mrSges(5,2) * t886 + Ifges(5,5) * t931 + Ifges(5,6) * t930 + Ifges(5,3) * t952 + pkin(4) * t865 + pkin(10) * t851 + t1000 * t844 + t1005 * t843 + t959 * t923 - t958 * t924;
t944 = Ifges(4,1) * t975 - Ifges(4,4) * t974 + Ifges(4,5) * t990;
t820 = -mrSges(4,1) * t933 + mrSges(4,3) * t898 + Ifges(4,4) * t954 + Ifges(4,2) * t953 + Ifges(4,6) * t989 - pkin(3) * t842 - t975 * t942 + t990 * t944 - t1011;
t962 = Ifges(3,3) * t990 + (Ifges(3,5) * t1002 + Ifges(3,6) * t1007) * t1034;
t964 = Ifges(3,5) * t990 + (Ifges(3,1) * t1002 + Ifges(3,4) * t1007) * t1034;
t809 = -mrSges(3,1) * t965 + mrSges(3,3) * t947 + Ifges(3,4) * t983 + Ifges(3,2) * t984 + Ifges(3,6) * t989 - pkin(2) * t841 + qJ(3) * t1022 - t962 * t1024 + t995 * t814 + t997 * t820 + t990 * t964;
t963 = Ifges(3,6) * t990 + (Ifges(3,4) * t1002 + Ifges(3,2) * t1007) * t1034;
t811 = mrSges(3,2) * t965 - mrSges(3,3) * t946 + Ifges(3,1) * t983 + Ifges(3,4) * t984 + Ifges(3,5) * t989 - qJ(3) * t834 + t962 * t1023 + t814 * t997 - t820 * t995 - t963 * t990;
t813 = Ifges(3,5) * t983 + Ifges(3,6) * t984 + mrSges(3,1) * t946 - mrSges(3,2) * t947 + Ifges(4,5) * t954 + Ifges(4,6) * t953 + t975 * t943 + t974 * t944 + mrSges(4,1) * t897 - mrSges(4,2) * t898 + t1001 * t826 + t1006 * t835 + pkin(3) * t1012 + pkin(9) * t1019 + pkin(2) * t834 + (Ifges(3,3) + Ifges(4,3)) * t989 + (t1002 * t963 - t1007 * t964) * t1034;
t1014 = mrSges(2,1) * t986 - mrSges(2,2) * t987 + Ifges(2,3) * qJDD(1) + pkin(1) * t819 + t809 * t1030 + t811 * t1032 + t825 * t1036 + t998 * t813;
t807 = -mrSges(2,2) * g(3) - mrSges(2,3) * t986 + Ifges(2,5) * qJDD(1) - t1009 * Ifges(2,6) - t1002 * t809 + t1007 * t811 + (-t818 * t996 - t819 * t998) * pkin(8);
t806 = mrSges(2,1) * g(3) + mrSges(2,3) * t987 + t1009 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t818 - t996 * t813 + (pkin(8) * t825 + t1002 * t811 + t1007 * t809) * t998;
t1 = [-m(1) * g(1) + t1017; -m(1) * g(2) + t1035; (-m(1) - m(2)) * g(3) + t818; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1035 - t1003 * t806 + t1008 * t807; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1017 + t1003 * t807 + t1008 * t806; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1014; t1014; t813; t841; t1011; t1010; -t1013;];
tauJB  = t1;
