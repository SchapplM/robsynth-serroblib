% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 15:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:39:24
% EndTime: 2019-05-06 15:39:57
% DurationCPUTime: 24.55s
% Computational Cost: add. (388041->380), mult. (884863->475), div. (0->0), fcn. (703512->12), ass. (0->159)
t1045 = Ifges(5,1) + Ifges(6,2);
t1037 = Ifges(5,4) + Ifges(6,6);
t1036 = Ifges(5,5) - Ifges(6,4);
t1044 = -Ifges(5,2) - Ifges(6,3);
t1035 = Ifges(5,6) - Ifges(6,5);
t1043 = Ifges(5,3) + Ifges(6,1);
t1002 = cos(qJ(6));
t1040 = cos(qJ(4));
t1003 = cos(qJ(2));
t1021 = qJD(1) * t1003;
t995 = sin(pkin(6));
t1018 = t995 * t1021;
t1000 = sin(qJ(2));
t997 = cos(pkin(6));
t1024 = t1000 * t997;
t1005 = qJD(1) ^ 2;
t1039 = pkin(8) * t995;
t1001 = sin(qJ(1));
t1004 = cos(qJ(1));
t985 = t1001 * g(1) - g(2) * t1004;
t975 = qJDD(1) * pkin(1) + t1005 * t1039 + t985;
t1020 = qJDD(1) * t995;
t986 = -g(1) * t1004 - g(2) * t1001;
t976 = -pkin(1) * t1005 + pkin(8) * t1020 + t986;
t1027 = t1003 * t976 + t975 * t1024;
t1026 = qJD(1) * t995;
t977 = (-pkin(2) * t1003 - qJ(3) * t1000) * t1026;
t990 = qJD(1) * t997 + qJD(2);
t988 = t990 ^ 2;
t989 = qJDD(1) * t997 + qJDD(2);
t926 = -t988 * pkin(2) + t989 * qJ(3) + (-g(3) * t1000 + t977 * t1021) * t995 + t1027;
t1038 = t997 * g(3);
t979 = (qJD(2) * t1021 + qJDD(1) * t1000) * t995;
t1025 = t1000 * t995;
t1019 = qJD(1) * t1025;
t980 = -qJD(2) * t1019 + t1003 * t1020;
t927 = -t980 * pkin(2) - t1038 - t979 * qJ(3) + (-t975 + (pkin(2) * t1000 - qJ(3) * t1003) * t990 * qJD(1)) * t995;
t994 = sin(pkin(11));
t996 = cos(pkin(11));
t968 = t996 * t1019 + t990 * t994;
t885 = -0.2e1 * qJD(3) * t968 - t994 * t926 + t996 * t927;
t955 = t979 * t996 + t989 * t994;
t967 = -t994 * t1019 + t990 * t996;
t880 = (-t967 * t1018 - t955) * pkin(9) + (t967 * t968 - t980) * pkin(3) + t885;
t886 = 0.2e1 * qJD(3) * t967 + t996 * t926 + t994 * t927;
t954 = -t979 * t994 + t989 * t996;
t956 = -pkin(3) * t1018 - pkin(9) * t968;
t966 = t967 ^ 2;
t883 = -pkin(3) * t966 + pkin(9) * t954 + t956 * t1018 + t886;
t999 = sin(qJ(4));
t876 = t1040 * t883 + t999 * t880;
t947 = -t1040 * t967 + t968 * t999;
t948 = t1040 * t968 + t999 * t967;
t919 = pkin(4) * t947 - qJ(5) * t948;
t972 = qJDD(4) - t980;
t983 = -qJD(4) + t1018;
t982 = t983 ^ 2;
t1013 = -pkin(4) * t982 + qJ(5) * t972 - t919 * t947 + t876;
t1041 = -2 * qJD(5);
t906 = qJD(4) * t948 - t1040 * t954 + t955 * t999;
t935 = pkin(5) * t948 + pkin(10) * t983;
t946 = t947 ^ 2;
t871 = -pkin(5) * t906 - pkin(10) * t946 + (t1041 - t935) * t983 + t1013;
t998 = sin(qJ(6));
t930 = -t1002 * t983 + t947 * t998;
t890 = -qJD(6) * t930 + t1002 * t906 - t972 * t998;
t929 = t1002 * t947 + t983 * t998;
t891 = qJD(6) * t929 + t1002 * t972 + t906 * t998;
t945 = qJD(6) + t948;
t908 = -mrSges(7,2) * t945 + mrSges(7,3) * t929;
t909 = mrSges(7,1) * t945 - mrSges(7,3) * t930;
t1014 = -m(7) * t871 + mrSges(7,1) * t890 - t891 * mrSges(7,2) + t908 * t929 - t930 * t909;
t873 = 0.2e1 * qJD(5) * t983 - t1013;
t932 = mrSges(6,1) * t948 - mrSges(6,2) * t983;
t1009 = -m(6) * t873 + t972 * mrSges(6,3) - t983 * t932 - t1014;
t1028 = -t1036 * t983 - t1037 * t947 + t1045 * t948;
t1029 = -t1035 * t983 + t1037 * t948 + t1044 * t947;
t1033 = t947 * t983;
t875 = t1040 * t880 - t999 * t883;
t874 = -t972 * pkin(4) - t982 * qJ(5) + t948 * t919 + qJDD(5) - t875;
t907 = -t947 * qJD(4) + t1040 * t955 + t999 * t954;
t869 = (t947 * t948 - t972) * pkin(10) + (t907 - t1033) * pkin(5) + t874;
t1022 = t1003 * t997;
t1023 = t1003 * t995;
t943 = -g(3) * t1023 - t1000 * t976 + t975 * t1022;
t925 = -pkin(2) * t989 - qJ(3) * t988 + t977 * t1019 + qJDD(3) - t943;
t892 = -pkin(3) * t954 - pkin(9) * t966 + t968 * t956 + t925;
t1006 = (-t907 - t1033) * qJ(5) + t892 + (-t983 * pkin(4) + t1041) * t948;
t872 = t1006 + (pkin(4) + pkin(10)) * t906 - pkin(5) * t946 - t935 * t948;
t867 = t1002 * t869 - t872 * t998;
t899 = -mrSges(7,1) * t929 + mrSges(7,2) * t930;
t905 = qJDD(6) + t907;
t864 = m(7) * t867 + mrSges(7,1) * t905 - mrSges(7,3) * t891 - t899 * t930 + t908 * t945;
t868 = t1002 * t872 + t869 * t998;
t865 = m(7) * t868 - mrSges(7,2) * t905 + mrSges(7,3) * t890 + t899 * t929 - t909 * t945;
t855 = t1002 * t864 + t865 * t998;
t921 = -mrSges(6,2) * t947 - mrSges(6,3) * t948;
t1011 = -m(6) * t874 - t907 * mrSges(6,1) - t948 * t921 - t855;
t931 = mrSges(6,1) * t947 + mrSges(6,3) * t983;
t853 = mrSges(6,2) * t972 - t931 * t983 - t1011;
t893 = Ifges(7,5) * t930 + Ifges(7,6) * t929 + Ifges(7,3) * t945;
t895 = Ifges(7,1) * t930 + Ifges(7,4) * t929 + Ifges(7,5) * t945;
t856 = -mrSges(7,1) * t871 + mrSges(7,3) * t868 + Ifges(7,4) * t891 + Ifges(7,2) * t890 + Ifges(7,6) * t905 - t893 * t930 + t895 * t945;
t894 = Ifges(7,4) * t930 + Ifges(7,2) * t929 + Ifges(7,6) * t945;
t857 = mrSges(7,2) * t871 - mrSges(7,3) * t867 + Ifges(7,1) * t891 + Ifges(7,4) * t890 + Ifges(7,5) * t905 + t893 * t929 - t894 * t945;
t1042 = t1028 * t947 + t1029 * t948 + t1043 * t972 - t1035 * t906 + t1036 * t907 + mrSges(5,1) * t875 - mrSges(5,2) * t876 + mrSges(6,2) * t874 - mrSges(6,3) * t873 - pkin(4) * t853 - pkin(10) * t855 + qJ(5) * (-mrSges(6,1) * t906 - t921 * t947 + t1009) + t1002 * t857 - t998 * t856;
t920 = mrSges(5,1) * t947 + mrSges(5,2) * t948;
t933 = mrSges(5,2) * t983 - mrSges(5,3) * t947;
t848 = m(5) * t875 - mrSges(5,3) * t907 - t920 * t948 + (t931 - t933) * t983 + (mrSges(5,1) - mrSges(6,2)) * t972 + t1011;
t934 = -mrSges(5,1) * t983 - mrSges(5,3) * t948;
t860 = m(5) * t876 - mrSges(5,2) * t972 + t934 * t983 + (-t920 - t921) * t947 + (-mrSges(5,3) - mrSges(6,1)) * t906 + t1009;
t846 = t1040 * t848 + t999 * t860;
t949 = -mrSges(4,1) * t967 + mrSges(4,2) * t968;
t952 = mrSges(4,2) * t1018 + mrSges(4,3) * t967;
t844 = m(4) * t885 - mrSges(4,1) * t980 - mrSges(4,3) * t955 - t952 * t1018 - t949 * t968 + t846;
t1016 = t1040 * t860 - t848 * t999;
t953 = -mrSges(4,1) * t1018 - mrSges(4,3) * t968;
t845 = m(4) * t886 + mrSges(4,2) * t980 + mrSges(4,3) * t954 + t953 * t1018 + t949 * t967 + t1016;
t1017 = -t844 * t994 + t996 * t845;
t944 = -g(3) * t1025 + t1027;
t973 = mrSges(3,1) * t990 - mrSges(3,3) * t1019;
t978 = (-mrSges(3,1) * t1003 + mrSges(3,2) * t1000) * t1026;
t836 = m(3) * t944 - mrSges(3,2) * t989 + mrSges(3,3) * t980 + t978 * t1018 - t973 * t990 + t1017;
t839 = t996 * t844 + t994 * t845;
t960 = -t995 * t975 - t1038;
t974 = -mrSges(3,2) * t990 + mrSges(3,3) * t1018;
t838 = m(3) * t960 - mrSges(3,1) * t980 + mrSges(3,2) * t979 + (t1000 * t973 - t1003 * t974) * t1026 + t839;
t1031 = t1002 * t865 - t998 * t864;
t878 = pkin(4) * t906 + t1006;
t854 = m(6) * t878 - t906 * mrSges(6,2) - t907 * mrSges(6,3) - t947 * t931 - t948 * t932 + t1031;
t1008 = m(5) * t892 + t906 * mrSges(5,1) + t907 * mrSges(5,2) + t947 * t933 + t948 * t934 + t854;
t852 = m(4) * t925 - t954 * mrSges(4,1) + t955 * mrSges(4,2) - t967 * t952 + t968 * t953 + t1008;
t851 = m(3) * t943 + t989 * mrSges(3,1) - t979 * mrSges(3,3) - t978 * t1019 + t990 * t974 - t852;
t826 = t851 * t1022 + t836 * t1024 - t838 * t995;
t823 = m(2) * t985 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1005 + t826;
t831 = -t1000 * t851 + t1003 * t836;
t829 = m(2) * t986 - mrSges(2,1) * t1005 - qJDD(1) * mrSges(2,2) + t831;
t1032 = t1001 * t829 + t1004 * t823;
t1030 = t1035 * t947 - t1036 * t948 + t1043 * t983;
t825 = t851 * t1023 + t836 * t1025 + t997 * t838;
t1015 = -t1001 * t823 + t1004 * t829;
t832 = -mrSges(5,1) * t892 - mrSges(6,1) * t873 + mrSges(6,2) * t878 + mrSges(5,3) * t876 - pkin(4) * t854 - pkin(5) * t1014 - pkin(10) * t1031 - t1002 * t856 - t1028 * t983 + t1030 * t948 + t1035 * t972 + t1037 * t907 + t1044 * t906 - t998 * t857;
t1010 = mrSges(7,1) * t867 - mrSges(7,2) * t868 + Ifges(7,5) * t891 + Ifges(7,6) * t890 + Ifges(7,3) * t905 + t930 * t894 - t929 * t895;
t840 = mrSges(6,1) * t874 + mrSges(5,2) * t892 - mrSges(5,3) * t875 - mrSges(6,3) * t878 + pkin(5) * t855 - qJ(5) * t854 + t1029 * t983 + t1030 * t947 + t1036 * t972 - t1037 * t906 + t1045 * t907 + t1010;
t937 = Ifges(4,5) * t968 + Ifges(4,6) * t967 - Ifges(4,3) * t1018;
t939 = Ifges(4,1) * t968 + Ifges(4,4) * t967 - Ifges(4,5) * t1018;
t818 = -mrSges(4,1) * t925 + mrSges(4,3) * t886 + Ifges(4,4) * t955 + Ifges(4,2) * t954 - Ifges(4,6) * t980 - pkin(3) * t1008 + pkin(9) * t1016 - t939 * t1018 + t1040 * t832 + t999 * t840 - t968 * t937;
t938 = Ifges(4,4) * t968 + Ifges(4,2) * t967 - Ifges(4,6) * t1018;
t821 = mrSges(4,2) * t925 - mrSges(4,3) * t885 + Ifges(4,1) * t955 + Ifges(4,4) * t954 - Ifges(4,5) * t980 - pkin(9) * t846 + t938 * t1018 + t1040 * t840 - t999 * t832 + t967 * t937;
t958 = Ifges(3,6) * t990 + (Ifges(3,4) * t1000 + Ifges(3,2) * t1003) * t1026;
t959 = Ifges(3,5) * t990 + (Ifges(3,1) * t1000 + Ifges(3,4) * t1003) * t1026;
t815 = Ifges(3,5) * t979 + Ifges(3,6) * t980 + Ifges(3,3) * t989 + mrSges(3,1) * t943 - mrSges(3,2) * t944 + t994 * t821 + t996 * t818 - pkin(2) * t852 + qJ(3) * t1017 + (t1000 * t958 - t1003 * t959) * t1026;
t957 = Ifges(3,3) * t990 + (Ifges(3,5) * t1000 + Ifges(3,6) * t1003) * t1026;
t817 = mrSges(3,2) * t960 - mrSges(3,3) * t943 + Ifges(3,1) * t979 + Ifges(3,4) * t980 + Ifges(3,5) * t989 - qJ(3) * t839 + t957 * t1018 - t818 * t994 + t821 * t996 - t958 * t990;
t820 = (Ifges(3,2) + Ifges(4,3)) * t980 - pkin(3) * t846 - t1042 - pkin(2) * t839 + Ifges(3,4) * t979 - mrSges(3,1) * t960 + t967 * t939 - t968 * t938 + Ifges(3,6) * t989 + t990 * t959 - Ifges(4,6) * t954 - Ifges(4,5) * t955 + mrSges(3,3) * t944 - mrSges(4,1) * t885 + mrSges(4,2) * t886 - t957 * t1019;
t1012 = mrSges(2,1) * t985 - mrSges(2,2) * t986 + Ifges(2,3) * qJDD(1) + pkin(1) * t826 + t820 * t1023 + t817 * t1025 + t831 * t1039 + t997 * t815;
t813 = -mrSges(2,2) * g(3) - mrSges(2,3) * t985 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t1005 - t1000 * t820 + t1003 * t817 + (-t825 * t995 - t826 * t997) * pkin(8);
t812 = mrSges(2,1) * g(3) + mrSges(2,3) * t986 + Ifges(2,5) * t1005 + Ifges(2,6) * qJDD(1) - pkin(1) * t825 - t815 * t995 + (pkin(8) * t831 + t1000 * t817 + t1003 * t820) * t997;
t1 = [-m(1) * g(1) + t1015; -m(1) * g(2) + t1032; (-m(1) - m(2)) * g(3) + t825; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1032 - t1001 * t812 + t1004 * t813; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1015 + t1001 * t813 + t1004 * t812; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1012; t1012; t815; t852; t1042; t853; t1010;];
tauJB  = t1;
