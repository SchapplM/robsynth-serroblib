% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 13:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:34:10
% EndTime: 2019-05-06 13:35:21
% DurationCPUTime: 57.28s
% Computational Cost: add. (828770->399), mult. (2201658->519), div. (0->0), fcn. (1758536->14), ass. (0->165)
t1029 = -2 * qJD(3);
t988 = sin(pkin(6));
t1019 = qJD(1) * t988;
t987 = sin(pkin(11));
t990 = cos(pkin(11));
t994 = sin(qJ(2));
t998 = cos(qJ(2));
t965 = (t987 * t994 - t990 * t998) * t1019;
t991 = cos(pkin(6));
t1021 = t991 * t998;
t1000 = qJD(1) ^ 2;
t1026 = pkin(8) * t988;
t995 = sin(qJ(1));
t999 = cos(qJ(1));
t977 = t995 * g(1) - g(2) * t999;
t971 = qJDD(1) * pkin(1) + t1000 * t1026 + t977;
t978 = -g(1) * t999 - g(2) * t995;
t972 = -pkin(1) * t1000 + qJDD(1) * t1026 + t978;
t1007 = t971 * t1021 - t994 * t972;
t1017 = t1000 * t988 ^ 2;
t1016 = qJD(1) * qJD(2);
t974 = (qJDD(1) * t994 + t998 * t1016) * t988;
t980 = qJDD(1) * t991 + qJDD(2);
t981 = qJD(1) * t991 + qJD(2);
t914 = t980 * pkin(2) - t974 * qJ(3) + (pkin(2) * t994 * t1017 + (qJ(3) * qJD(1) * t981 - g(3)) * t988) * t998 + t1007;
t1013 = t998 ^ 2 * t1017;
t1022 = t991 * t994;
t1024 = t988 * t994;
t942 = -g(3) * t1024 + t971 * t1022 + t998 * t972;
t1015 = t994 * t1019;
t968 = pkin(2) * t981 - qJ(3) * t1015;
t975 = (qJDD(1) * t998 - t994 * t1016) * t988;
t917 = -pkin(2) * t1013 + qJ(3) * t975 - t968 * t981 + t942;
t966 = (t987 * t998 + t990 * t994) * t1019;
t891 = t966 * t1029 + t990 * t914 - t987 * t917;
t1027 = 2 * qJD(5);
t892 = t965 * t1029 + t987 * t914 + t990 * t917;
t944 = pkin(3) * t965 - pkin(9) * t966;
t979 = t981 ^ 2;
t886 = -pkin(3) * t979 + pkin(9) * t980 - t944 * t965 + t892;
t957 = -t991 * g(3) - t988 * t971;
t928 = -t975 * pkin(2) - qJ(3) * t1013 + t968 * t1015 + qJDD(3) + t957;
t947 = -t974 * t987 + t975 * t990;
t948 = t974 * t990 + t975 * t987;
t895 = (t965 * t981 - t948) * pkin(9) + (t966 * t981 - t947) * pkin(3) + t928;
t993 = sin(qJ(4));
t997 = cos(qJ(4));
t878 = -t993 * t886 + t997 * t895;
t950 = -t966 * t993 + t981 * t997;
t925 = qJD(4) * t950 + t948 * t997 + t980 * t993;
t946 = qJDD(4) - t947;
t951 = t966 * t997 + t981 * t993;
t964 = qJD(4) + t965;
t875 = (t950 * t964 - t925) * qJ(5) + (t950 * t951 + t946) * pkin(4) + t878;
t879 = t997 * t886 + t993 * t895;
t924 = -qJD(4) * t951 - t948 * t993 + t980 * t997;
t935 = pkin(4) * t964 - qJ(5) * t951;
t949 = t950 ^ 2;
t877 = -pkin(4) * t949 + qJ(5) * t924 - t935 * t964 + t879;
t986 = sin(pkin(12));
t989 = cos(pkin(12));
t930 = t950 * t989 - t951 * t986;
t872 = t930 * t1027 + t986 * t875 + t989 * t877;
t931 = t950 * t986 + t951 * t989;
t908 = -pkin(5) * t930 - pkin(10) * t931;
t963 = t964 ^ 2;
t870 = -pkin(5) * t963 + pkin(10) * t946 + t908 * t930 + t872;
t885 = -t980 * pkin(3) - t979 * pkin(9) + t966 * t944 - t891;
t880 = -t924 * pkin(4) - t949 * qJ(5) + t951 * t935 + qJDD(5) + t885;
t899 = t924 * t989 - t925 * t986;
t900 = t924 * t986 + t925 * t989;
t873 = (-t930 * t964 - t900) * pkin(10) + (t931 * t964 - t899) * pkin(5) + t880;
t992 = sin(qJ(6));
t996 = cos(qJ(6));
t867 = -t870 * t992 + t873 * t996;
t910 = -t931 * t992 + t964 * t996;
t883 = qJD(6) * t910 + t900 * t996 + t946 * t992;
t911 = t931 * t996 + t964 * t992;
t896 = -mrSges(7,1) * t910 + mrSges(7,2) * t911;
t898 = qJDD(6) - t899;
t929 = qJD(6) - t930;
t901 = -mrSges(7,2) * t929 + mrSges(7,3) * t910;
t864 = m(7) * t867 + mrSges(7,1) * t898 - mrSges(7,3) * t883 - t896 * t911 + t901 * t929;
t868 = t870 * t996 + t873 * t992;
t882 = -qJD(6) * t911 - t900 * t992 + t946 * t996;
t902 = mrSges(7,1) * t929 - mrSges(7,3) * t911;
t865 = m(7) * t868 - mrSges(7,2) * t898 + mrSges(7,3) * t882 + t896 * t910 - t902 * t929;
t856 = -t864 * t992 + t996 * t865;
t907 = -mrSges(6,1) * t930 + mrSges(6,2) * t931;
t916 = mrSges(6,1) * t964 - mrSges(6,3) * t931;
t853 = m(6) * t872 - mrSges(6,2) * t946 + mrSges(6,3) * t899 + t907 * t930 - t916 * t964 + t856;
t1006 = -t989 * t875 + t986 * t877;
t869 = -t946 * pkin(5) - t963 * pkin(10) + (t1027 + t908) * t931 + t1006;
t866 = -m(7) * t869 + t882 * mrSges(7,1) - mrSges(7,2) * t883 + t910 * t901 - t902 * t911;
t871 = -0.2e1 * qJD(5) * t931 - t1006;
t915 = -mrSges(6,2) * t964 + mrSges(6,3) * t930;
t860 = m(6) * t871 + mrSges(6,1) * t946 - mrSges(6,3) * t900 - t907 * t931 + t915 * t964 + t866;
t848 = t986 * t853 + t989 * t860;
t887 = Ifges(7,5) * t911 + Ifges(7,6) * t910 + Ifges(7,3) * t929;
t889 = Ifges(7,1) * t911 + Ifges(7,4) * t910 + Ifges(7,5) * t929;
t857 = -mrSges(7,1) * t869 + mrSges(7,3) * t868 + Ifges(7,4) * t883 + Ifges(7,2) * t882 + Ifges(7,6) * t898 - t887 * t911 + t889 * t929;
t888 = Ifges(7,4) * t911 + Ifges(7,2) * t910 + Ifges(7,6) * t929;
t858 = mrSges(7,2) * t869 - mrSges(7,3) * t867 + Ifges(7,1) * t883 + Ifges(7,4) * t882 + Ifges(7,5) * t898 + t887 * t910 - t888 * t929;
t904 = Ifges(6,4) * t931 + Ifges(6,2) * t930 + Ifges(6,6) * t964;
t905 = Ifges(6,1) * t931 + Ifges(6,4) * t930 + Ifges(6,5) * t964;
t919 = Ifges(5,4) * t951 + Ifges(5,2) * t950 + Ifges(5,6) * t964;
t920 = Ifges(5,1) * t951 + Ifges(5,4) * t950 + Ifges(5,5) * t964;
t1028 = Ifges(5,5) * t925 + Ifges(5,6) * t924 + t951 * t919 - t950 * t920 + mrSges(5,1) * t878 - mrSges(5,2) * t879 + Ifges(6,5) * t900 + Ifges(6,6) * t899 + t931 * t904 - t930 * t905 + mrSges(6,1) * t871 - mrSges(6,2) * t872 + t992 * t858 + t996 * t857 + pkin(5) * t866 + pkin(10) * t856 + pkin(4) * t848 + (Ifges(5,3) + Ifges(6,3)) * t946;
t1023 = t988 * t998;
t932 = -mrSges(5,1) * t950 + mrSges(5,2) * t951;
t934 = -mrSges(5,2) * t964 + mrSges(5,3) * t950;
t846 = m(5) * t878 + mrSges(5,1) * t946 - mrSges(5,3) * t925 - t932 * t951 + t934 * t964 + t848;
t1009 = t989 * t853 - t860 * t986;
t936 = mrSges(5,1) * t964 - mrSges(5,3) * t951;
t847 = m(5) * t879 - mrSges(5,2) * t946 + mrSges(5,3) * t924 + t932 * t950 - t936 * t964 + t1009;
t1010 = -t846 * t993 + t997 * t847;
t943 = mrSges(4,1) * t965 + mrSges(4,2) * t966;
t953 = mrSges(4,1) * t981 - mrSges(4,3) * t966;
t836 = m(4) * t892 - mrSges(4,2) * t980 + mrSges(4,3) * t947 - t943 * t965 - t953 * t981 + t1010;
t855 = t996 * t864 + t992 * t865;
t854 = m(6) * t880 - t899 * mrSges(6,1) + mrSges(6,2) * t900 - t930 * t915 + t916 * t931 + t855;
t1002 = -m(5) * t885 + t924 * mrSges(5,1) - mrSges(5,2) * t925 + t950 * t934 - t936 * t951 - t854;
t952 = -mrSges(4,2) * t981 - mrSges(4,3) * t965;
t850 = m(4) * t891 + mrSges(4,1) * t980 - mrSges(4,3) * t948 - t943 * t966 + t952 * t981 + t1002;
t833 = t987 * t836 + t990 * t850;
t941 = -g(3) * t1023 + t1007;
t1014 = t998 * t1019;
t970 = -mrSges(3,2) * t981 + mrSges(3,3) * t1014;
t973 = (-mrSges(3,1) * t998 + mrSges(3,2) * t994) * t1019;
t831 = m(3) * t941 + mrSges(3,1) * t980 - mrSges(3,3) * t974 - t973 * t1015 + t970 * t981 + t833;
t1011 = t990 * t836 - t850 * t987;
t969 = mrSges(3,1) * t981 - mrSges(3,3) * t1015;
t832 = m(3) * t942 - mrSges(3,2) * t980 + mrSges(3,3) * t975 + t973 * t1014 - t969 * t981 + t1011;
t840 = t997 * t846 + t993 * t847;
t839 = m(4) * t928 - t947 * mrSges(4,1) + t948 * mrSges(4,2) + t965 * t952 + t966 * t953 + t840;
t838 = m(3) * t957 - t975 * mrSges(3,1) + t974 * mrSges(3,2) + (t969 * t994 - t970 * t998) * t1019 + t839;
t817 = t831 * t1021 + t832 * t1022 - t838 * t988;
t814 = m(2) * t977 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1000 + t817;
t822 = -t831 * t994 + t998 * t832;
t820 = m(2) * t978 - mrSges(2,1) * t1000 - qJDD(1) * mrSges(2,2) + t822;
t1020 = t999 * t814 + t995 * t820;
t816 = t831 * t1023 + t832 * t1024 + t991 * t838;
t1012 = -t814 * t995 + t999 * t820;
t903 = Ifges(6,5) * t931 + Ifges(6,6) * t930 + Ifges(6,3) * t964;
t841 = mrSges(6,2) * t880 - mrSges(6,3) * t871 + Ifges(6,1) * t900 + Ifges(6,4) * t899 + Ifges(6,5) * t946 - pkin(10) * t855 - t857 * t992 + t858 * t996 + t903 * t930 - t904 * t964;
t1003 = mrSges(7,1) * t867 - mrSges(7,2) * t868 + Ifges(7,5) * t883 + Ifges(7,6) * t882 + Ifges(7,3) * t898 + t888 * t911 - t889 * t910;
t842 = -mrSges(6,1) * t880 + mrSges(6,3) * t872 + Ifges(6,4) * t900 + Ifges(6,2) * t899 + Ifges(6,6) * t946 - pkin(5) * t855 - t903 * t931 + t905 * t964 - t1003;
t918 = Ifges(5,5) * t951 + Ifges(5,6) * t950 + Ifges(5,3) * t964;
t824 = -mrSges(5,1) * t885 + mrSges(5,3) * t879 + Ifges(5,4) * t925 + Ifges(5,2) * t924 + Ifges(5,6) * t946 - pkin(4) * t854 + qJ(5) * t1009 + t986 * t841 + t989 * t842 - t951 * t918 + t964 * t920;
t825 = mrSges(5,2) * t885 - mrSges(5,3) * t878 + Ifges(5,1) * t925 + Ifges(5,4) * t924 + Ifges(5,5) * t946 - qJ(5) * t848 + t841 * t989 - t842 * t986 + t918 * t950 - t919 * t964;
t937 = Ifges(4,5) * t966 - Ifges(4,6) * t965 + Ifges(4,3) * t981;
t938 = Ifges(4,4) * t966 - Ifges(4,2) * t965 + Ifges(4,6) * t981;
t812 = mrSges(4,2) * t928 - mrSges(4,3) * t891 + Ifges(4,1) * t948 + Ifges(4,4) * t947 + Ifges(4,5) * t980 - pkin(9) * t840 - t824 * t993 + t825 * t997 - t937 * t965 - t938 * t981;
t939 = Ifges(4,1) * t966 - Ifges(4,4) * t965 + Ifges(4,5) * t981;
t823 = -mrSges(4,1) * t928 + mrSges(4,3) * t892 + Ifges(4,4) * t948 + Ifges(4,2) * t947 + Ifges(4,6) * t980 - pkin(3) * t840 - t966 * t937 + t981 * t939 - t1028;
t954 = Ifges(3,3) * t981 + (Ifges(3,5) * t994 + Ifges(3,6) * t998) * t1019;
t956 = Ifges(3,5) * t981 + (Ifges(3,1) * t994 + Ifges(3,4) * t998) * t1019;
t807 = -mrSges(3,1) * t957 + mrSges(3,3) * t942 + Ifges(3,4) * t974 + Ifges(3,2) * t975 + Ifges(3,6) * t980 - pkin(2) * t839 + qJ(3) * t1011 - t954 * t1015 + t987 * t812 + t990 * t823 + t981 * t956;
t955 = Ifges(3,6) * t981 + (Ifges(3,4) * t994 + Ifges(3,2) * t998) * t1019;
t809 = mrSges(3,2) * t957 - mrSges(3,3) * t941 + Ifges(3,1) * t974 + Ifges(3,4) * t975 + Ifges(3,5) * t980 - qJ(3) * t833 + t954 * t1014 + t812 * t990 - t823 * t987 - t955 * t981;
t811 = Ifges(3,5) * t974 + Ifges(3,6) * t975 + mrSges(3,1) * t941 - mrSges(3,2) * t942 + Ifges(4,5) * t948 + Ifges(4,6) * t947 + t966 * t938 + t965 * t939 + mrSges(4,1) * t891 - mrSges(4,2) * t892 + t993 * t825 + t997 * t824 + pkin(3) * t1002 + pkin(9) * t1010 + pkin(2) * t833 + (Ifges(3,3) + Ifges(4,3)) * t980 + (t955 * t994 - t956 * t998) * t1019;
t1004 = mrSges(2,1) * t977 - mrSges(2,2) * t978 + Ifges(2,3) * qJDD(1) + pkin(1) * t817 + t807 * t1023 + t809 * t1024 + t822 * t1026 + t991 * t811;
t805 = -mrSges(2,2) * g(3) - mrSges(2,3) * t977 + Ifges(2,5) * qJDD(1) - t1000 * Ifges(2,6) - t994 * t807 + t998 * t809 + (-t816 * t988 - t817 * t991) * pkin(8);
t804 = mrSges(2,1) * g(3) + mrSges(2,3) * t978 + t1000 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t816 - t988 * t811 + (pkin(8) * t822 + t807 * t998 + t809 * t994) * t991;
t1 = [-m(1) * g(1) + t1012; -m(1) * g(2) + t1020; (-m(1) - m(2)) * g(3) + t816; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1020 - t995 * t804 + t999 * t805; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1012 + t999 * t804 + t995 * t805; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1004; t1004; t811; t839; t1028; t854; t1003;];
tauJB  = t1;
