% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:47:24
% EndTime: 2019-05-06 17:47:53
% DurationCPUTime: 27.33s
% Computational Cost: add. (417004->377), mult. (1090889->476), div. (0->0), fcn. (858536->12), ass. (0->157)
t1024 = -2 * qJD(3);
t1023 = Ifges(6,1) + Ifges(7,1);
t1017 = Ifges(6,4) + Ifges(7,4);
t1016 = Ifges(6,5) + Ifges(7,5);
t1022 = Ifges(6,2) + Ifges(7,2);
t1015 = Ifges(6,6) + Ifges(7,6);
t1021 = Ifges(6,3) + Ifges(7,3);
t973 = sin(pkin(6));
t1003 = qJD(1) * t973;
t972 = sin(pkin(11));
t974 = cos(pkin(11));
t978 = sin(qJ(2));
t982 = cos(qJ(2));
t951 = (t972 * t978 - t974 * t982) * t1003;
t984 = qJD(1) ^ 2;
t1013 = t973 ^ 2 * t984;
t1001 = qJD(1) * qJD(2);
t960 = (qJDD(1) * t978 + t1001 * t982) * t973;
t975 = cos(pkin(6));
t966 = qJDD(1) * t975 + qJDD(2);
t967 = qJD(1) * t975 + qJD(2);
t1009 = t975 * t982;
t1019 = pkin(8) * t973;
t979 = sin(qJ(1));
t983 = cos(qJ(1));
t963 = t979 * g(1) - g(2) * t983;
t957 = qJDD(1) * pkin(1) + t1019 * t984 + t963;
t964 = -g(1) * t983 - g(2) * t979;
t958 = -pkin(1) * t984 + qJDD(1) * t1019 + t964;
t991 = t957 * t1009 - t978 * t958;
t897 = t966 * pkin(2) - t960 * qJ(3) + (pkin(2) * t978 * t1013 + (qJ(3) * qJD(1) * t967 - g(3)) * t973) * t982 + t991;
t1000 = t982 ^ 2 * t1013;
t1010 = t975 * t978;
t1012 = t973 * t978;
t928 = -g(3) * t1012 + t957 * t1010 + t982 * t958;
t997 = t978 * t1003;
t954 = pkin(2) * t967 - qJ(3) * t997;
t961 = (qJDD(1) * t982 - t1001 * t978) * t973;
t903 = -pkin(2) * t1000 + qJ(3) * t961 - t954 * t967 + t928;
t952 = (t972 * t982 + t974 * t978) * t1003;
t872 = t1024 * t952 + t974 * t897 - t972 * t903;
t977 = sin(qJ(4));
t981 = cos(qJ(4));
t937 = t952 * t981 + t967 * t977;
t950 = qJD(4) + t951;
t976 = sin(qJ(5));
t980 = cos(qJ(5));
t919 = -t937 * t976 + t950 * t980;
t920 = t937 * t980 + t950 * t976;
t936 = -t952 * t977 + t967 * t981;
t935 = qJD(5) - t936;
t1005 = t1016 * t935 + t1017 * t919 + t1023 * t920;
t1006 = -t1015 * t935 - t1017 * t920 - t1022 * t919;
t934 = t960 * t974 + t961 * t972;
t912 = qJD(4) * t936 + t934 * t981 + t966 * t977;
t933 = -t960 * t972 + t961 * t974;
t932 = qJDD(4) - t933;
t880 = qJD(5) * t919 + t912 * t980 + t932 * t976;
t891 = -mrSges(7,1) * t919 + mrSges(7,2) * t920;
t873 = t1024 * t951 + t972 * t897 + t974 * t903;
t930 = pkin(3) * t951 - pkin(9) * t952;
t965 = t967 ^ 2;
t871 = -pkin(3) * t965 + pkin(9) * t966 - t930 * t951 + t873;
t943 = -t975 * g(3) - t973 * t957;
t914 = -t961 * pkin(2) - qJ(3) * t1000 + t954 * t997 + qJDD(3) + t943;
t875 = (t951 * t967 - t934) * pkin(9) + (t952 * t967 - t933) * pkin(3) + t914;
t867 = t981 * t871 + t977 * t875;
t916 = -pkin(4) * t936 - pkin(10) * t937;
t949 = t950 ^ 2;
t862 = -pkin(4) * t949 + pkin(10) * t932 + t916 * t936 + t867;
t870 = -t966 * pkin(3) - t965 * pkin(9) + t952 * t930 - t872;
t911 = -qJD(4) * t937 - t934 * t977 + t966 * t981;
t865 = (-t936 * t950 - t912) * pkin(10) + (t937 * t950 - t911) * pkin(4) + t870;
t857 = -t976 * t862 + t980 * t865;
t910 = qJDD(5) - t911;
t853 = -0.2e1 * qJD(6) * t920 + (t919 * t935 - t880) * qJ(6) + (t919 * t920 + t910) * pkin(5) + t857;
t898 = -mrSges(7,2) * t935 + mrSges(7,3) * t919;
t999 = m(7) * t853 + t910 * mrSges(7,1) + t935 * t898;
t851 = -t880 * mrSges(7,3) - t920 * t891 + t999;
t858 = t980 * t862 + t976 * t865;
t879 = -qJD(5) * t920 - t912 * t976 + t932 * t980;
t900 = pkin(5) * t935 - qJ(6) * t920;
t918 = t919 ^ 2;
t856 = -pkin(5) * t918 + qJ(6) * t879 + 0.2e1 * qJD(6) * t919 - t900 * t935 + t858;
t1020 = mrSges(6,1) * t857 + mrSges(7,1) * t853 - mrSges(6,2) * t858 - mrSges(7,2) * t856 + pkin(5) * t851 - t1005 * t919 - t1006 * t920 + t1015 * t879 + t1016 * t880 + t1021 * t910;
t1018 = -mrSges(6,2) - mrSges(7,2);
t1011 = t973 * t982;
t929 = mrSges(4,1) * t951 + mrSges(4,2) * t952;
t939 = mrSges(4,1) * t967 - mrSges(4,3) * t952;
t892 = -mrSges(6,1) * t919 + mrSges(6,2) * t920;
t899 = -mrSges(6,2) * t935 + mrSges(6,3) * t919;
t845 = m(6) * t857 + t910 * mrSges(6,1) + t935 * t899 + (-t891 - t892) * t920 + (-mrSges(6,3) - mrSges(7,3)) * t880 + t999;
t901 = mrSges(7,1) * t935 - mrSges(7,3) * t920;
t1004 = -mrSges(6,1) * t935 + mrSges(6,3) * t920 - t901;
t998 = m(7) * t856 + t879 * mrSges(7,3) + t919 * t891;
t847 = m(6) * t858 + t879 * mrSges(6,3) + t1004 * t935 + t1018 * t910 + t919 * t892 + t998;
t844 = -t845 * t976 + t980 * t847;
t915 = -mrSges(5,1) * t936 + mrSges(5,2) * t937;
t922 = mrSges(5,1) * t950 - mrSges(5,3) * t937;
t841 = m(5) * t867 - mrSges(5,2) * t932 + mrSges(5,3) * t911 + t915 * t936 - t922 * t950 + t844;
t866 = -t977 * t871 + t875 * t981;
t861 = -pkin(4) * t932 - pkin(10) * t949 + t937 * t916 - t866;
t859 = -pkin(5) * t879 - qJ(6) * t918 + t900 * t920 + qJDD(6) + t861;
t990 = -m(7) * t859 + t879 * mrSges(7,1) + t919 * t898;
t850 = -m(6) * t861 + t879 * mrSges(6,1) + t1004 * t920 + t1018 * t880 + t919 * t899 + t990;
t921 = -mrSges(5,2) * t950 + mrSges(5,3) * t936;
t849 = m(5) * t866 + t932 * mrSges(5,1) - t912 * mrSges(5,3) - t937 * t915 + t950 * t921 + t850;
t993 = t981 * t841 - t849 * t977;
t831 = m(4) * t873 - mrSges(4,2) * t966 + mrSges(4,3) * t933 - t929 * t951 - t939 * t967 + t993;
t938 = -mrSges(4,2) * t967 - mrSges(4,3) * t951;
t843 = t845 * t980 + t847 * t976;
t986 = -m(5) * t870 + t911 * mrSges(5,1) - mrSges(5,2) * t912 + t936 * t921 - t922 * t937 - t843;
t838 = m(4) * t872 + mrSges(4,1) * t966 - mrSges(4,3) * t934 - t929 * t952 + t938 * t967 + t986;
t827 = t972 * t831 + t974 * t838;
t927 = -g(3) * t1011 + t991;
t996 = t982 * t1003;
t956 = -mrSges(3,2) * t967 + mrSges(3,3) * t996;
t959 = (-mrSges(3,1) * t982 + mrSges(3,2) * t978) * t1003;
t825 = m(3) * t927 + mrSges(3,1) * t966 - mrSges(3,3) * t960 + t956 * t967 - t959 * t997 + t827;
t955 = mrSges(3,1) * t967 - mrSges(3,3) * t997;
t994 = t974 * t831 - t838 * t972;
t826 = m(3) * t928 - mrSges(3,2) * t966 + mrSges(3,3) * t961 - t955 * t967 + t959 * t996 + t994;
t835 = t977 * t841 + t981 * t849;
t834 = m(4) * t914 - t933 * mrSges(4,1) + t934 * mrSges(4,2) + t951 * t938 + t952 * t939 + t835;
t833 = m(3) * t943 - t961 * mrSges(3,1) + t960 * mrSges(3,2) + (t955 * t978 - t956 * t982) * t1003 + t834;
t812 = t825 * t1009 + t826 * t1010 - t833 * t973;
t809 = m(2) * t963 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t984 + t812;
t818 = -t825 * t978 + t982 * t826;
t816 = m(2) * t964 - mrSges(2,1) * t984 - qJDD(1) * mrSges(2,2) + t818;
t1008 = t983 * t809 + t979 * t816;
t1007 = -t1015 * t919 - t1016 * t920 - t1021 * t935;
t811 = t825 * t1011 + t826 * t1012 + t975 * t833;
t995 = -t809 * t979 + t983 * t816;
t854 = t880 * mrSges(7,2) + t920 * t901 - t990;
t836 = -mrSges(6,1) * t861 + mrSges(6,3) * t858 - mrSges(7,1) * t859 + mrSges(7,3) * t856 - pkin(5) * t854 + qJ(6) * t998 + (-qJ(6) * t901 + t1005) * t935 + t1007 * t920 + (-mrSges(7,2) * qJ(6) + t1015) * t910 + t1017 * t880 + t1022 * t879;
t842 = mrSges(6,2) * t861 + mrSges(7,2) * t859 - mrSges(6,3) * t857 - mrSges(7,3) * t853 - qJ(6) * t851 + t1006 * t935 - t1007 * t919 + t1016 * t910 + t1017 * t879 + t1023 * t880;
t904 = Ifges(5,5) * t937 + Ifges(5,6) * t936 + Ifges(5,3) * t950;
t905 = Ifges(5,4) * t937 + Ifges(5,2) * t936 + Ifges(5,6) * t950;
t819 = mrSges(5,2) * t870 - mrSges(5,3) * t866 + Ifges(5,1) * t912 + Ifges(5,4) * t911 + Ifges(5,5) * t932 - pkin(10) * t843 - t836 * t976 + t842 * t980 + t904 * t936 - t905 * t950;
t906 = Ifges(5,1) * t937 + Ifges(5,4) * t936 + Ifges(5,5) * t950;
t828 = -mrSges(5,1) * t870 + mrSges(5,3) * t867 + Ifges(5,4) * t912 + Ifges(5,2) * t911 + Ifges(5,6) * t932 - pkin(4) * t843 - t937 * t904 + t950 * t906 - t1020;
t923 = Ifges(4,5) * t952 - Ifges(4,6) * t951 + Ifges(4,3) * t967;
t924 = Ifges(4,4) * t952 - Ifges(4,2) * t951 + Ifges(4,6) * t967;
t807 = mrSges(4,2) * t914 - mrSges(4,3) * t872 + Ifges(4,1) * t934 + Ifges(4,4) * t933 + Ifges(4,5) * t966 - pkin(9) * t835 + t819 * t981 - t828 * t977 - t923 * t951 - t924 * t967;
t925 = Ifges(4,1) * t952 - Ifges(4,4) * t951 + Ifges(4,5) * t967;
t985 = mrSges(5,1) * t866 - mrSges(5,2) * t867 + Ifges(5,5) * t912 + Ifges(5,6) * t911 + Ifges(5,3) * t932 + pkin(4) * t850 + pkin(10) * t844 + t980 * t836 + t976 * t842 + t937 * t905 - t936 * t906;
t813 = -mrSges(4,1) * t914 + mrSges(4,3) * t873 + Ifges(4,4) * t934 + Ifges(4,2) * t933 + Ifges(4,6) * t966 - pkin(3) * t835 - t952 * t923 + t967 * t925 - t985;
t940 = Ifges(3,3) * t967 + (Ifges(3,5) * t978 + Ifges(3,6) * t982) * t1003;
t942 = Ifges(3,5) * t967 + (Ifges(3,1) * t978 + Ifges(3,4) * t982) * t1003;
t802 = -mrSges(3,1) * t943 + mrSges(3,3) * t928 + Ifges(3,4) * t960 + Ifges(3,2) * t961 + Ifges(3,6) * t966 - pkin(2) * t834 + qJ(3) * t994 + t972 * t807 + t974 * t813 - t940 * t997 + t967 * t942;
t941 = Ifges(3,6) * t967 + (Ifges(3,4) * t978 + Ifges(3,2) * t982) * t1003;
t804 = mrSges(3,2) * t943 - mrSges(3,3) * t927 + Ifges(3,1) * t960 + Ifges(3,4) * t961 + Ifges(3,5) * t966 - qJ(3) * t827 + t807 * t974 - t813 * t972 + t940 * t996 - t941 * t967;
t806 = Ifges(3,5) * t960 + Ifges(3,6) * t961 + mrSges(3,1) * t927 - mrSges(3,2) * t928 + Ifges(4,5) * t934 + Ifges(4,6) * t933 + t952 * t924 + t951 * t925 + mrSges(4,1) * t872 - mrSges(4,2) * t873 + t977 * t819 + t981 * t828 + pkin(3) * t986 + pkin(9) * t993 + pkin(2) * t827 + (Ifges(3,3) + Ifges(4,3)) * t966 + (t941 * t978 - t942 * t982) * t1003;
t988 = mrSges(2,1) * t963 - mrSges(2,2) * t964 + Ifges(2,3) * qJDD(1) + pkin(1) * t812 + t802 * t1011 + t804 * t1012 + t818 * t1019 + t975 * t806;
t800 = -mrSges(2,2) * g(3) - mrSges(2,3) * t963 + Ifges(2,5) * qJDD(1) - t984 * Ifges(2,6) - t978 * t802 + t982 * t804 + (-t811 * t973 - t812 * t975) * pkin(8);
t799 = mrSges(2,1) * g(3) + mrSges(2,3) * t964 + t984 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t811 - t973 * t806 + (pkin(8) * t818 + t802 * t982 + t804 * t978) * t975;
t1 = [-m(1) * g(1) + t995; -m(1) * g(2) + t1008; (-m(1) - m(2)) * g(3) + t811; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1008 - t979 * t799 + t983 * t800; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t995 + t983 * t799 + t979 * t800; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t988; t988; t806; t834; t985; t1020; t854;];
tauJB  = t1;
