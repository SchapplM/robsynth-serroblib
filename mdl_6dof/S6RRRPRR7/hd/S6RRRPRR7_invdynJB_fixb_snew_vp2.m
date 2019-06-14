% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 11:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:33:27
% EndTime: 2019-05-07 11:35:21
% DurationCPUTime: 76.84s
% Computational Cost: add. (1283442->401), mult. (2845044->518), div. (0->0), fcn. (2327637->14), ass. (0->166)
t1002 = cos(qJ(6));
t1003 = cos(qJ(5));
t1004 = cos(qJ(3));
t1000 = sin(qJ(2));
t1005 = cos(qJ(2));
t1023 = qJD(1) * t1005;
t996 = cos(pkin(6));
t1026 = t1000 * t996;
t1007 = qJD(1) ^ 2;
t994 = sin(pkin(6));
t1033 = pkin(8) * t994;
t1001 = sin(qJ(1));
t1006 = cos(qJ(1));
t984 = t1001 * g(1) - g(2) * t1006;
t974 = qJDD(1) * pkin(1) + t1007 * t1033 + t984;
t1022 = qJDD(1) * t994;
t985 = -g(1) * t1006 - g(2) * t1001;
t975 = -pkin(1) * t1007 + pkin(8) * t1022 + t985;
t1029 = t1005 * t975 + t974 * t1026;
t1028 = qJD(1) * t994;
t977 = (-pkin(2) * t1005 - pkin(9) * t1000) * t1028;
t989 = qJD(1) * t996 + qJD(2);
t987 = t989 ^ 2;
t988 = qJDD(1) * t996 + qJDD(2);
t934 = -t987 * pkin(2) + t988 * pkin(9) + (-g(3) * t1000 + t977 * t1023) * t994 + t1029;
t1032 = t996 * g(3);
t978 = (qJD(2) * t1023 + qJDD(1) * t1000) * t994;
t1027 = t1000 * t994;
t1021 = qJD(1) * t1027;
t979 = -qJD(2) * t1021 + t1005 * t1022;
t935 = -t979 * pkin(2) - t978 * pkin(9) - t1032 + (-t974 + (pkin(2) * t1000 - pkin(9) * t1005) * t989 * qJD(1)) * t994;
t999 = sin(qJ(3));
t912 = t1004 * t935 - t999 * t934;
t966 = t1004 * t989 - t999 * t1021;
t946 = qJD(3) * t966 + t1004 * t978 + t988 * t999;
t967 = t1004 * t1021 + t989 * t999;
t971 = qJDD(3) - t979;
t1020 = t994 * t1023;
t983 = qJD(3) - t1020;
t897 = (t966 * t983 - t946) * qJ(4) + (t966 * t967 + t971) * pkin(3) + t912;
t913 = t1004 * t934 + t999 * t935;
t945 = -qJD(3) * t967 + t1004 * t988 - t978 * t999;
t956 = pkin(3) * t983 - qJ(4) * t967;
t965 = t966 ^ 2;
t900 = -pkin(3) * t965 + qJ(4) * t945 - t956 * t983 + t913;
t993 = sin(pkin(12));
t995 = cos(pkin(12));
t953 = t966 * t993 + t967 * t995;
t877 = -0.2e1 * qJD(4) * t953 + t995 * t897 - t993 * t900;
t920 = t945 * t993 + t946 * t995;
t952 = t966 * t995 - t967 * t993;
t874 = (t952 * t983 - t920) * pkin(10) + (t952 * t953 + t971) * pkin(4) + t877;
t878 = 0.2e1 * qJD(4) * t952 + t993 * t897 + t995 * t900;
t919 = t945 * t995 - t946 * t993;
t938 = pkin(4) * t983 - pkin(10) * t953;
t951 = t952 ^ 2;
t876 = -pkin(4) * t951 + pkin(10) * t919 - t938 * t983 + t878;
t998 = sin(qJ(5));
t870 = t1003 * t874 - t876 * t998;
t927 = t1003 * t952 - t953 * t998;
t928 = t1003 * t953 + t952 * t998;
t911 = -pkin(5) * t927 - pkin(11) * t928;
t970 = qJDD(5) + t971;
t981 = qJD(5) + t983;
t980 = t981 ^ 2;
t867 = -pkin(5) * t970 - pkin(11) * t980 + t911 * t928 - t870;
t896 = qJD(5) * t927 + t1003 * t920 + t919 * t998;
t997 = sin(qJ(6));
t915 = t1002 * t928 + t981 * t997;
t882 = -qJD(6) * t915 + t1002 * t970 - t896 * t997;
t914 = t1002 * t981 - t928 * t997;
t883 = qJD(6) * t914 + t1002 * t896 + t970 * t997;
t926 = qJD(6) - t927;
t902 = -mrSges(7,2) * t926 + mrSges(7,3) * t914;
t903 = mrSges(7,1) * t926 - mrSges(7,3) * t915;
t1013 = -m(7) * t867 + t882 * mrSges(7,1) - mrSges(7,2) * t883 + t914 * t902 - t903 * t915;
t871 = t1003 * t876 + t998 * t874;
t868 = -pkin(5) * t980 + pkin(11) * t970 + t911 * t927 + t871;
t1024 = t1005 * t996;
t1025 = t1005 * t994;
t947 = -g(3) * t1025 - t1000 * t975 + t974 * t1024;
t933 = -t988 * pkin(2) - t987 * pkin(9) + t977 * t1021 - t947;
t904 = -t945 * pkin(3) - t965 * qJ(4) + t967 * t956 + qJDD(4) + t933;
t880 = -t919 * pkin(4) - t951 * pkin(10) + t953 * t938 + t904;
t895 = -qJD(5) * t928 + t1003 * t919 - t920 * t998;
t872 = t880 + (-t927 * t981 - t896) * pkin(11) + (t928 * t981 - t895) * pkin(5);
t865 = t1002 * t872 - t868 * t997;
t894 = qJDD(6) - t895;
t901 = -mrSges(7,1) * t914 + mrSges(7,2) * t915;
t861 = m(7) * t865 + mrSges(7,1) * t894 - mrSges(7,3) * t883 - t901 * t915 + t902 * t926;
t866 = t1002 * t868 + t872 * t997;
t862 = m(7) * t866 - mrSges(7,2) * t894 + mrSges(7,3) * t882 + t901 * t914 - t903 * t926;
t1016 = t1002 * t862 - t861 * t997;
t884 = Ifges(7,5) * t915 + Ifges(7,6) * t914 + Ifges(7,3) * t926;
t886 = Ifges(7,1) * t915 + Ifges(7,4) * t914 + Ifges(7,5) * t926;
t854 = -mrSges(7,1) * t867 + mrSges(7,3) * t866 + Ifges(7,4) * t883 + Ifges(7,2) * t882 + Ifges(7,6) * t894 - t884 * t915 + t886 * t926;
t885 = Ifges(7,4) * t915 + Ifges(7,2) * t914 + Ifges(7,6) * t926;
t855 = mrSges(7,2) * t867 - mrSges(7,3) * t865 + Ifges(7,1) * t883 + Ifges(7,4) * t882 + Ifges(7,5) * t894 + t884 * t914 - t885 * t926;
t906 = Ifges(6,4) * t928 + Ifges(6,2) * t927 + Ifges(6,6) * t981;
t907 = Ifges(6,1) * t928 + Ifges(6,4) * t927 + Ifges(6,5) * t981;
t1011 = -mrSges(6,1) * t870 + mrSges(6,2) * t871 - Ifges(6,5) * t896 - Ifges(6,6) * t895 - Ifges(6,3) * t970 - pkin(5) * t1013 - pkin(11) * t1016 - t1002 * t854 - t997 * t855 - t928 * t906 + t927 * t907;
t910 = -mrSges(6,1) * t927 + mrSges(6,2) * t928;
t917 = mrSges(6,1) * t981 - mrSges(6,3) * t928;
t844 = m(6) * t871 - mrSges(6,2) * t970 + mrSges(6,3) * t895 + t910 * t927 - t917 * t981 + t1016;
t916 = -mrSges(6,2) * t981 + mrSges(6,3) * t927;
t857 = m(6) * t870 + mrSges(6,1) * t970 - mrSges(6,3) * t896 - t910 * t928 + t916 * t981 + t1013;
t841 = t1003 * t857 + t998 * t844;
t929 = -mrSges(5,1) * t952 + mrSges(5,2) * t953;
t936 = -mrSges(5,2) * t983 + mrSges(5,3) * t952;
t839 = m(5) * t877 + mrSges(5,1) * t971 - mrSges(5,3) * t920 - t929 * t953 + t936 * t983 + t841;
t1017 = t1003 * t844 - t857 * t998;
t937 = mrSges(5,1) * t983 - mrSges(5,3) * t953;
t840 = m(5) * t878 - mrSges(5,2) * t971 + mrSges(5,3) * t919 + t929 * t952 - t937 * t983 + t1017;
t833 = t995 * t839 + t993 * t840;
t922 = Ifges(5,4) * t953 + Ifges(5,2) * t952 + Ifges(5,6) * t983;
t923 = Ifges(5,1) * t953 + Ifges(5,4) * t952 + Ifges(5,5) * t983;
t940 = Ifges(4,4) * t967 + Ifges(4,2) * t966 + Ifges(4,6) * t983;
t941 = Ifges(4,1) * t967 + Ifges(4,4) * t966 + Ifges(4,5) * t983;
t1034 = mrSges(4,1) * t912 + mrSges(5,1) * t877 - mrSges(4,2) * t913 - mrSges(5,2) * t878 + Ifges(4,5) * t946 + Ifges(5,5) * t920 + Ifges(4,6) * t945 + Ifges(5,6) * t919 + pkin(3) * t833 + pkin(4) * t841 + (Ifges(4,3) + Ifges(5,3)) * t971 + t953 * t922 - t952 * t923 + t967 * t940 - t966 * t941 - t1011;
t954 = -mrSges(4,1) * t966 + mrSges(4,2) * t967;
t955 = -mrSges(4,2) * t983 + mrSges(4,3) * t966;
t831 = m(4) * t912 + mrSges(4,1) * t971 - mrSges(4,3) * t946 - t954 * t967 + t955 * t983 + t833;
t1018 = -t839 * t993 + t995 * t840;
t957 = mrSges(4,1) * t983 - mrSges(4,3) * t967;
t832 = m(4) * t913 - mrSges(4,2) * t971 + mrSges(4,3) * t945 + t954 * t966 - t957 * t983 + t1018;
t1019 = t1004 * t832 - t831 * t999;
t948 = -g(3) * t1027 + t1029;
t972 = mrSges(3,1) * t989 - mrSges(3,3) * t1021;
t976 = (-mrSges(3,1) * t1005 + mrSges(3,2) * t1000) * t1028;
t823 = m(3) * t948 - mrSges(3,2) * t988 + mrSges(3,3) * t979 + t976 * t1020 - t972 * t989 + t1019;
t826 = t1004 * t831 + t999 * t832;
t961 = -t994 * t974 - t1032;
t973 = -mrSges(3,2) * t989 + mrSges(3,3) * t1020;
t825 = m(3) * t961 - t979 * mrSges(3,1) + t978 * mrSges(3,2) + (t1000 * t972 - t1005 * t973) * t1028 + t826;
t850 = t1002 * t861 + t997 * t862;
t1014 = m(6) * t880 - t895 * mrSges(6,1) + t896 * mrSges(6,2) - t927 * t916 + t928 * t917 + t850;
t848 = m(5) * t904 - t919 * mrSges(5,1) + t920 * mrSges(5,2) - t952 * t936 + t953 * t937 + t1014;
t1009 = -m(4) * t933 + t945 * mrSges(4,1) - t946 * mrSges(4,2) + t966 * t955 - t967 * t957 - t848;
t847 = m(3) * t947 + t988 * mrSges(3,1) - t978 * mrSges(3,3) - t976 * t1021 + t989 * t973 + t1009;
t813 = t847 * t1024 + t823 * t1026 - t825 * t994;
t810 = m(2) * t984 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1007 + t813;
t818 = -t1000 * t847 + t1005 * t823;
t816 = m(2) * t985 - mrSges(2,1) * t1007 - qJDD(1) * mrSges(2,2) + t818;
t1030 = t1001 * t816 + t1006 * t810;
t812 = t847 * t1025 + t823 * t1027 + t996 * t825;
t1015 = -t1001 * t810 + t1006 * t816;
t905 = Ifges(6,5) * t928 + Ifges(6,6) * t927 + Ifges(6,3) * t981;
t834 = mrSges(6,2) * t880 - mrSges(6,3) * t870 + Ifges(6,1) * t896 + Ifges(6,4) * t895 + Ifges(6,5) * t970 - pkin(11) * t850 + t1002 * t855 - t854 * t997 + t905 * t927 - t906 * t981;
t1010 = mrSges(7,1) * t865 - mrSges(7,2) * t866 + Ifges(7,5) * t883 + Ifges(7,6) * t882 + Ifges(7,3) * t894 + t885 * t915 - t886 * t914;
t835 = -mrSges(6,1) * t880 + mrSges(6,3) * t871 + Ifges(6,4) * t896 + Ifges(6,2) * t895 + Ifges(6,6) * t970 - pkin(5) * t850 - t905 * t928 + t907 * t981 - t1010;
t921 = Ifges(5,5) * t953 + Ifges(5,6) * t952 + Ifges(5,3) * t983;
t819 = -mrSges(5,1) * t904 + mrSges(5,3) * t878 + Ifges(5,4) * t920 + Ifges(5,2) * t919 + Ifges(5,6) * t971 - pkin(4) * t1014 + pkin(10) * t1017 + t1003 * t835 + t998 * t834 - t953 * t921 + t983 * t923;
t827 = mrSges(5,2) * t904 - mrSges(5,3) * t877 + Ifges(5,1) * t920 + Ifges(5,4) * t919 + Ifges(5,5) * t971 - pkin(10) * t841 + t1003 * t834 - t835 * t998 + t921 * t952 - t922 * t983;
t939 = Ifges(4,5) * t967 + Ifges(4,6) * t966 + Ifges(4,3) * t983;
t805 = -mrSges(4,1) * t933 + mrSges(4,3) * t913 + Ifges(4,4) * t946 + Ifges(4,2) * t945 + Ifges(4,6) * t971 - pkin(3) * t848 + qJ(4) * t1018 + t995 * t819 + t993 * t827 - t967 * t939 + t983 * t941;
t806 = mrSges(4,2) * t933 - mrSges(4,3) * t912 + Ifges(4,1) * t946 + Ifges(4,4) * t945 + Ifges(4,5) * t971 - qJ(4) * t833 - t819 * t993 + t827 * t995 + t939 * t966 - t940 * t983;
t959 = Ifges(3,6) * t989 + (Ifges(3,4) * t1000 + Ifges(3,2) * t1005) * t1028;
t960 = Ifges(3,5) * t989 + (Ifges(3,1) * t1000 + Ifges(3,4) * t1005) * t1028;
t802 = Ifges(3,5) * t978 + Ifges(3,6) * t979 + Ifges(3,3) * t988 + mrSges(3,1) * t947 - mrSges(3,2) * t948 + t999 * t806 + t1004 * t805 + pkin(2) * t1009 + pkin(9) * t1019 + (t1000 * t959 - t1005 * t960) * t1028;
t958 = Ifges(3,3) * t989 + (Ifges(3,5) * t1000 + Ifges(3,6) * t1005) * t1028;
t804 = mrSges(3,2) * t961 - mrSges(3,3) * t947 + Ifges(3,1) * t978 + Ifges(3,4) * t979 + Ifges(3,5) * t988 - pkin(9) * t826 + t1004 * t806 + t958 * t1020 - t805 * t999 - t959 * t989;
t808 = -mrSges(3,1) * t961 + mrSges(3,3) * t948 + Ifges(3,4) * t978 + Ifges(3,2) * t979 + Ifges(3,6) * t988 - pkin(2) * t826 - t958 * t1021 + t989 * t960 - t1034;
t1012 = mrSges(2,1) * t984 - mrSges(2,2) * t985 + Ifges(2,3) * qJDD(1) + pkin(1) * t813 + t808 * t1025 + t804 * t1027 + t818 * t1033 + t996 * t802;
t800 = -mrSges(2,2) * g(3) - mrSges(2,3) * t984 + Ifges(2,5) * qJDD(1) - t1007 * Ifges(2,6) - t1000 * t808 + t1005 * t804 + (-t812 * t994 - t813 * t996) * pkin(8);
t799 = mrSges(2,1) * g(3) + mrSges(2,3) * t985 + t1007 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t812 - t994 * t802 + (pkin(8) * t818 + t1000 * t804 + t1005 * t808) * t996;
t1 = [-m(1) * g(1) + t1015; -m(1) * g(2) + t1030; (-m(1) - m(2)) * g(3) + t812; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1030 - t1001 * t799 + t1006 * t800; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1015 + t1001 * t800 + t1006 * t799; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1012; t1012; t802; t1034; t848; -t1011; t1010;];
tauJB  = t1;
