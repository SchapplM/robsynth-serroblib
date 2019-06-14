% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 08:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:15:37
% EndTime: 2019-05-07 08:16:19
% DurationCPUTime: 31.00s
% Computational Cost: add. (514717->376), mult. (1134214->477), div. (0->0), fcn. (905306->12), ass. (0->156)
t1025 = Ifges(6,1) + Ifges(7,1);
t1016 = Ifges(6,4) - Ifges(7,5);
t1015 = -Ifges(6,5) - Ifges(7,4);
t1024 = Ifges(6,2) + Ifges(7,3);
t1014 = Ifges(6,6) - Ifges(7,6);
t1023 = -Ifges(6,3) - Ifges(7,2);
t974 = sin(pkin(6));
t979 = sin(qJ(2));
t982 = cos(qJ(2));
t998 = qJD(1) * qJD(2);
t960 = (-qJDD(1) * t982 + t979 * t998) * t974;
t1000 = qJD(1) * t982;
t976 = cos(pkin(6));
t1009 = t976 * t979;
t1019 = pkin(8) * t974;
t980 = sin(qJ(1));
t983 = cos(qJ(1));
t965 = t980 * g(1) - g(2) * t983;
t984 = qJD(1) ^ 2;
t955 = qJDD(1) * pkin(1) + t1019 * t984 + t965;
t966 = -g(1) * t983 - g(2) * t980;
t956 = -pkin(1) * t984 + qJDD(1) * t1019 + t966;
t1002 = t955 * t1009 + t982 * t956;
t1001 = qJD(1) * t974;
t958 = (-pkin(2) * t982 - pkin(9) * t979) * t1001;
t969 = qJD(1) * t976 + qJD(2);
t967 = t969 ^ 2;
t968 = qJDD(1) * t976 + qJDD(2);
t914 = -t967 * pkin(2) + t968 * pkin(9) + (-g(3) * t979 + t1000 * t958) * t974 + t1002;
t1018 = t976 * g(3);
t959 = (qJDD(1) * t979 + t982 * t998) * t974;
t915 = t960 * pkin(2) - t959 * pkin(9) - t1018 + (-t955 + (pkin(2) * t979 - pkin(9) * t982) * t969 * qJD(1)) * t974;
t978 = sin(qJ(3));
t981 = cos(qJ(3));
t878 = -t978 * t914 + t981 * t915;
t996 = t979 * t1001;
t948 = t969 * t981 - t978 * t996;
t927 = qJD(3) * t948 + t959 * t981 + t968 * t978;
t949 = t969 * t978 + t981 * t996;
t952 = qJDD(3) + t960;
t995 = t974 * t1000;
t964 = qJD(3) - t995;
t869 = (t948 * t964 - t927) * qJ(4) + (t948 * t949 + t952) * pkin(3) + t878;
t879 = t981 * t914 + t978 * t915;
t926 = -qJD(3) * t949 - t959 * t978 + t968 * t981;
t938 = pkin(3) * t964 - qJ(4) * t949;
t947 = t948 ^ 2;
t872 = -pkin(3) * t947 + qJ(4) * t926 - t938 * t964 + t879;
t973 = sin(pkin(11));
t975 = cos(pkin(11));
t935 = t948 * t973 + t949 * t975;
t864 = -0.2e1 * qJD(4) * t935 + t975 * t869 - t973 * t872;
t1020 = cos(qJ(5));
t977 = sin(qJ(5));
t916 = -t1020 * t964 + t935 * t977;
t917 = t1020 * t935 + t977 * t964;
t892 = mrSges(7,1) * t916 - mrSges(7,3) * t917;
t1003 = -mrSges(6,1) * t916 - mrSges(6,2) * t917 - t892;
t1017 = -mrSges(6,3) - mrSges(7,2);
t934 = t948 * t975 - t949 * t973;
t865 = 0.2e1 * qJD(4) * t934 + t973 * t869 + t975 * t872;
t909 = -pkin(4) * t934 - pkin(10) * t935;
t963 = t964 ^ 2;
t863 = -pkin(4) * t963 + pkin(10) * t952 + t909 * t934 + t865;
t1008 = t976 * t982;
t1010 = t974 * t982;
t928 = -g(3) * t1010 + t1008 * t955 - t979 * t956;
t913 = -t968 * pkin(2) - t967 * pkin(9) + t958 * t996 - t928;
t873 = -t926 * pkin(3) - t947 * qJ(4) + t949 * t938 + qJDD(4) + t913;
t902 = t926 * t975 - t927 * t973;
t903 = t926 * t973 + t927 * t975;
t867 = (-t934 * t964 - t903) * pkin(10) + (t935 * t964 - t902) * pkin(4) + t873;
t860 = t1020 * t863 + t977 * t867;
t876 = qJD(5) * t917 - t1020 * t952 + t903 * t977;
t933 = qJD(5) - t934;
t896 = mrSges(6,1) * t933 - mrSges(6,3) * t917;
t901 = qJDD(5) - t902;
t891 = pkin(5) * t916 - qJ(6) * t917;
t932 = t933 ^ 2;
t856 = -pkin(5) * t932 + qJ(6) * t901 + 0.2e1 * qJD(6) * t933 - t891 * t916 + t860;
t897 = -mrSges(7,1) * t933 + mrSges(7,2) * t917;
t997 = m(7) * t856 + t901 * mrSges(7,3) + t933 * t897;
t848 = m(6) * t860 - t901 * mrSges(6,2) + t1003 * t916 + t1017 * t876 - t933 * t896 + t997;
t859 = t1020 * t867 - t977 * t863;
t877 = -t916 * qJD(5) + t1020 * t903 + t977 * t952;
t895 = -mrSges(6,2) * t933 - mrSges(6,3) * t916;
t857 = -t901 * pkin(5) - t932 * qJ(6) + t917 * t891 + qJDD(6) - t859;
t894 = -mrSges(7,2) * t916 + mrSges(7,3) * t933;
t990 = -m(7) * t857 + t901 * mrSges(7,1) + t933 * t894;
t850 = m(6) * t859 + t901 * mrSges(6,1) + t1003 * t917 + t1017 * t877 + t933 * t895 + t990;
t843 = t1020 * t848 - t850 * t977;
t908 = -mrSges(5,1) * t934 + mrSges(5,2) * t935;
t919 = mrSges(5,1) * t964 - mrSges(5,3) * t935;
t838 = m(5) * t865 - mrSges(5,2) * t952 + mrSges(5,3) * t902 + t908 * t934 - t919 * t964 + t843;
t862 = -t952 * pkin(4) - t963 * pkin(10) + t935 * t909 - t864;
t858 = -0.2e1 * qJD(6) * t917 + (t916 * t933 - t877) * qJ(6) + (t917 * t933 + t876) * pkin(5) + t862;
t854 = m(7) * t858 + mrSges(7,1) * t876 - t877 * mrSges(7,3) + t894 * t916 - t917 * t897;
t851 = -m(6) * t862 - t876 * mrSges(6,1) - mrSges(6,2) * t877 - t916 * t895 - t896 * t917 - t854;
t918 = -mrSges(5,2) * t964 + mrSges(5,3) * t934;
t845 = m(5) * t864 + mrSges(5,1) * t952 - mrSges(5,3) * t903 - t908 * t935 + t918 * t964 + t851;
t832 = t973 * t838 + t975 * t845;
t1004 = t1015 * t933 + t1016 * t916 - t1025 * t917;
t1006 = t1014 * t916 + t1015 * t917 + t1023 * t933;
t839 = -mrSges(6,1) * t862 - mrSges(7,1) * t858 + mrSges(7,2) * t856 + mrSges(6,3) * t860 - pkin(5) * t854 - t1004 * t933 + t1006 * t917 + t1014 * t901 + t1016 * t877 - t1024 * t876;
t1005 = -t1014 * t933 - t1016 * t917 + t1024 * t916;
t840 = mrSges(6,2) * t862 + mrSges(7,2) * t857 - mrSges(6,3) * t859 - mrSges(7,3) * t858 - qJ(6) * t854 + t1005 * t933 + t1006 * t916 - t1015 * t901 - t1016 * t876 + t1025 * t877;
t905 = Ifges(5,4) * t935 + Ifges(5,2) * t934 + Ifges(5,6) * t964;
t906 = Ifges(5,1) * t935 + Ifges(5,4) * t934 + Ifges(5,5) * t964;
t921 = Ifges(4,4) * t949 + Ifges(4,2) * t948 + Ifges(4,6) * t964;
t922 = Ifges(4,1) * t949 + Ifges(4,4) * t948 + Ifges(4,5) * t964;
t1022 = (Ifges(4,3) + Ifges(5,3)) * t952 + mrSges(4,1) * t878 + mrSges(5,1) * t864 + Ifges(4,5) * t927 + Ifges(5,5) * t903 + Ifges(4,6) * t926 + Ifges(5,6) * t902 + pkin(3) * t832 + pkin(4) * t851 + pkin(10) * t843 + t935 * t905 + t949 * t921 + t977 * t840 + t1020 * t839 - mrSges(4,2) * t879 - mrSges(5,2) * t865 - t934 * t906 - t948 * t922;
t853 = t877 * mrSges(7,2) + t917 * t892 - t990;
t1021 = -t1004 * t916 - t1005 * t917 - t1023 * t901 - t1014 * t876 - t1015 * t877 + mrSges(6,1) * t859 - mrSges(7,1) * t857 - mrSges(6,2) * t860 + mrSges(7,3) * t856 - pkin(5) * t853 + qJ(6) * (-t876 * mrSges(7,2) - t916 * t892 + t997);
t1011 = t974 * t979;
t929 = -g(3) * t1011 + t1002;
t953 = mrSges(3,1) * t969 - mrSges(3,3) * t996;
t957 = (-mrSges(3,1) * t982 + mrSges(3,2) * t979) * t1001;
t936 = -mrSges(4,1) * t948 + mrSges(4,2) * t949;
t937 = -mrSges(4,2) * t964 + mrSges(4,3) * t948;
t830 = m(4) * t878 + mrSges(4,1) * t952 - mrSges(4,3) * t927 - t936 * t949 + t937 * t964 + t832;
t939 = mrSges(4,1) * t964 - mrSges(4,3) * t949;
t992 = t975 * t838 - t845 * t973;
t831 = m(4) * t879 - mrSges(4,2) * t952 + mrSges(4,3) * t926 + t936 * t948 - t939 * t964 + t992;
t993 = -t830 * t978 + t981 * t831;
t821 = m(3) * t929 - mrSges(3,2) * t968 - mrSges(3,3) * t960 - t953 * t969 + t957 * t995 + t993;
t824 = t981 * t830 + t978 * t831;
t943 = -t974 * t955 - t1018;
t954 = -mrSges(3,2) * t969 + mrSges(3,3) * t995;
t823 = m(3) * t943 + t960 * mrSges(3,1) + t959 * mrSges(3,2) + (t953 * t979 - t954 * t982) * t1001 + t824;
t842 = t1020 * t850 + t977 * t848;
t841 = m(5) * t873 - t902 * mrSges(5,1) + mrSges(5,2) * t903 - t934 * t918 + t919 * t935 + t842;
t986 = -m(4) * t913 + t926 * mrSges(4,1) - mrSges(4,2) * t927 + t948 * t937 - t939 * t949 - t841;
t835 = m(3) * t928 + mrSges(3,1) * t968 - mrSges(3,3) * t959 + t954 * t969 - t957 * t996 + t986;
t811 = t835 * t1008 + t821 * t1009 - t823 * t974;
t808 = m(2) * t965 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t984 + t811;
t817 = t982 * t821 - t835 * t979;
t815 = m(2) * t966 - mrSges(2,1) * t984 - qJDD(1) * mrSges(2,2) + t817;
t1007 = t983 * t808 + t980 * t815;
t810 = t835 * t1010 + t821 * t1011 + t976 * t823;
t994 = -t808 * t980 + t983 * t815;
t904 = Ifges(5,5) * t935 + Ifges(5,6) * t934 + Ifges(5,3) * t964;
t825 = mrSges(5,2) * t873 - mrSges(5,3) * t864 + Ifges(5,1) * t903 + Ifges(5,4) * t902 + Ifges(5,5) * t952 - pkin(10) * t842 + t1020 * t840 - t977 * t839 + t934 * t904 - t964 * t905;
t826 = -mrSges(5,1) * t873 + mrSges(5,3) * t865 + Ifges(5,4) * t903 + Ifges(5,2) * t902 + Ifges(5,6) * t952 - pkin(4) * t842 - t935 * t904 + t964 * t906 - t1021;
t920 = Ifges(4,5) * t949 + Ifges(4,6) * t948 + Ifges(4,3) * t964;
t806 = -mrSges(4,1) * t913 + mrSges(4,3) * t879 + Ifges(4,4) * t927 + Ifges(4,2) * t926 + Ifges(4,6) * t952 - pkin(3) * t841 + qJ(4) * t992 + t973 * t825 + t975 * t826 - t949 * t920 + t964 * t922;
t812 = mrSges(4,2) * t913 - mrSges(4,3) * t878 + Ifges(4,1) * t927 + Ifges(4,4) * t926 + Ifges(4,5) * t952 - qJ(4) * t832 + t825 * t975 - t826 * t973 + t920 * t948 - t921 * t964;
t941 = Ifges(3,6) * t969 + (Ifges(3,4) * t979 + Ifges(3,2) * t982) * t1001;
t942 = Ifges(3,5) * t969 + (Ifges(3,1) * t979 + Ifges(3,4) * t982) * t1001;
t801 = Ifges(3,5) * t959 - Ifges(3,6) * t960 + Ifges(3,3) * t968 + mrSges(3,1) * t928 - mrSges(3,2) * t929 + t978 * t812 + t981 * t806 + pkin(2) * t986 + pkin(9) * t993 + (t941 * t979 - t942 * t982) * t1001;
t940 = Ifges(3,3) * t969 + (Ifges(3,5) * t979 + Ifges(3,6) * t982) * t1001;
t803 = mrSges(3,2) * t943 - mrSges(3,3) * t928 + Ifges(3,1) * t959 - Ifges(3,4) * t960 + Ifges(3,5) * t968 - pkin(9) * t824 - t806 * t978 + t812 * t981 + t940 * t995 - t941 * t969;
t805 = -mrSges(3,1) * t943 + mrSges(3,3) * t929 + Ifges(3,4) * t959 - Ifges(3,2) * t960 + Ifges(3,6) * t968 - pkin(2) * t824 - t940 * t996 + t969 * t942 - t1022;
t988 = mrSges(2,1) * t965 - mrSges(2,2) * t966 + Ifges(2,3) * qJDD(1) + pkin(1) * t811 + t805 * t1010 + t803 * t1011 + t817 * t1019 + t976 * t801;
t799 = -mrSges(2,2) * g(3) - mrSges(2,3) * t965 + Ifges(2,5) * qJDD(1) - t984 * Ifges(2,6) + t982 * t803 - t979 * t805 + (-t810 * t974 - t811 * t976) * pkin(8);
t798 = mrSges(2,1) * g(3) + mrSges(2,3) * t966 + t984 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t810 - t974 * t801 + (pkin(8) * t817 + t803 * t979 + t805 * t982) * t976;
t1 = [-m(1) * g(1) + t994; -m(1) * g(2) + t1007; (-m(1) - m(2)) * g(3) + t810; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1007 - t980 * t798 + t983 * t799; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t994 + t983 * t798 + t980 * t799; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t988; t988; t801; t1022; t841; t1021; t853;];
tauJB  = t1;
