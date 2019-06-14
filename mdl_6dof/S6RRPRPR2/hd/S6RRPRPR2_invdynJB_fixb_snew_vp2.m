% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 13:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:02:32
% EndTime: 2019-05-06 13:02:51
% DurationCPUTime: 14.14s
% Computational Cost: add. (189898->367), mult. (446187->447), div. (0->0), fcn. (324723->10), ass. (0->147)
t1001 = Ifges(5,4) + Ifges(6,6);
t1013 = -Ifges(5,2) - Ifges(6,3);
t1009 = Ifges(5,6) - Ifges(6,5);
t1003 = cos(qJ(4));
t963 = sin(pkin(10));
t964 = cos(pkin(10));
t967 = sin(qJ(2));
t970 = cos(qJ(2));
t935 = (-t963 * t967 + t964 * t970) * qJD(1);
t936 = (t963 * t970 + t964 * t967) * qJD(1);
t966 = sin(qJ(4));
t916 = -t1003 * t935 + t966 * t936;
t917 = t1003 * t936 + t966 * t935;
t960 = qJD(2) + qJD(4);
t1007 = t1001 * t917 + t1009 * t960 + t1013 * t916;
t1008 = Ifges(5,3) + Ifges(6,1);
t1010 = Ifges(5,5) - Ifges(6,4);
t906 = mrSges(6,1) * t916 - mrSges(6,3) * t960;
t959 = qJDD(2) + qJDD(4);
t972 = qJD(1) ^ 2;
t1002 = pkin(2) * t972;
t991 = qJD(1) * qJD(2);
t946 = qJDD(1) * t967 + t970 * t991;
t968 = sin(qJ(1));
t971 = cos(qJ(1));
t953 = -g(1) * t971 - g(2) * t968;
t941 = -pkin(1) * t972 + qJDD(1) * pkin(7) + t953;
t998 = t967 * t941;
t898 = qJDD(2) * pkin(2) - t946 * qJ(3) - t998 + (qJ(3) * t991 + t967 * t1002 - g(3)) * t970;
t926 = -g(3) * t967 + t970 * t941;
t947 = qJDD(1) * t970 - t967 * t991;
t993 = qJD(1) * t967;
t949 = qJD(2) * pkin(2) - qJ(3) * t993;
t962 = t970 ^ 2;
t899 = qJ(3) * t947 - qJD(2) * t949 - t962 * t1002 + t926;
t860 = -0.2e1 * qJD(3) * t936 + t964 * t898 - t963 * t899;
t924 = t946 * t964 + t947 * t963;
t844 = (qJD(2) * t935 - t924) * pkin(8) + (t935 * t936 + qJDD(2)) * pkin(3) + t860;
t861 = 0.2e1 * qJD(3) * t935 + t963 * t898 + t964 * t899;
t923 = -t946 * t963 + t947 * t964;
t929 = qJD(2) * pkin(3) - pkin(8) * t936;
t934 = t935 ^ 2;
t847 = -pkin(3) * t934 + pkin(8) * t923 - qJD(2) * t929 + t861;
t841 = t1003 * t844 - t966 * t847;
t892 = pkin(4) * t916 - qJ(5) * t917;
t958 = t960 ^ 2;
t837 = -t959 * pkin(4) - t958 * qJ(5) + t917 * t892 + qJDD(5) - t841;
t878 = -t916 * qJD(4) + t1003 * t924 + t966 * t923;
t999 = t916 * t960;
t831 = (t916 * t917 - t959) * pkin(9) + (t878 + t999) * pkin(5) + t837;
t877 = t917 * qJD(4) - t1003 * t923 + t966 * t924;
t908 = pkin(5) * t917 - pkin(9) * t960;
t912 = t916 ^ 2;
t1004 = -2 * qJD(5);
t952 = t968 * g(1) - t971 * g(2);
t984 = -qJDD(1) * pkin(1) - t952;
t901 = -t947 * pkin(2) + qJDD(3) + t949 * t993 + (-qJ(3) * t962 - pkin(7)) * t972 + t984;
t859 = -t923 * pkin(3) - t934 * pkin(8) + t936 * t929 + t901;
t975 = (-t878 + t999) * qJ(5) + t859 + (t960 * pkin(4) + t1004) * t917;
t834 = t975 - t917 * t908 - t912 * pkin(5) + (pkin(4) + pkin(9)) * t877;
t965 = sin(qJ(6));
t969 = cos(qJ(6));
t829 = t831 * t969 - t834 * t965;
t902 = t916 * t969 - t960 * t965;
t853 = qJD(6) * t902 + t877 * t965 + t959 * t969;
t876 = qJDD(6) + t878;
t903 = t916 * t965 + t960 * t969;
t879 = -mrSges(7,1) * t902 + mrSges(7,2) * t903;
t911 = qJD(6) + t917;
t880 = -mrSges(7,2) * t911 + mrSges(7,3) * t902;
t826 = m(7) * t829 + mrSges(7,1) * t876 - mrSges(7,3) * t853 - t879 * t903 + t880 * t911;
t830 = t831 * t965 + t834 * t969;
t852 = -qJD(6) * t903 + t877 * t969 - t959 * t965;
t881 = mrSges(7,1) * t911 - mrSges(7,3) * t903;
t827 = m(7) * t830 - mrSges(7,2) * t876 + mrSges(7,3) * t852 + t879 * t902 - t881 * t911;
t815 = t969 * t826 + t965 * t827;
t894 = -mrSges(6,2) * t916 - mrSges(6,3) * t917;
t981 = -m(6) * t837 - t878 * mrSges(6,1) - t917 * t894 - t815;
t814 = t959 * mrSges(6,2) + t960 * t906 - t981;
t842 = t1003 * t847 + t966 * t844;
t980 = -t958 * pkin(4) + t959 * qJ(5) - t916 * t892 + t842;
t833 = -t877 * pkin(5) - t912 * pkin(9) + ((2 * qJD(5)) + t908) * t960 + t980;
t854 = Ifges(7,5) * t903 + Ifges(7,6) * t902 + Ifges(7,3) * t911;
t856 = Ifges(7,1) * t903 + Ifges(7,4) * t902 + Ifges(7,5) * t911;
t817 = -mrSges(7,1) * t833 + mrSges(7,3) * t830 + Ifges(7,4) * t853 + Ifges(7,2) * t852 + Ifges(7,6) * t876 - t854 * t903 + t856 * t911;
t855 = Ifges(7,4) * t903 + Ifges(7,2) * t902 + Ifges(7,6) * t911;
t818 = mrSges(7,2) * t833 - mrSges(7,3) * t829 + Ifges(7,1) * t853 + Ifges(7,4) * t852 + Ifges(7,5) * t876 + t854 * t902 - t855 * t911;
t835 = t960 * t1004 - t980;
t907 = mrSges(6,1) * t917 + mrSges(6,2) * t960;
t983 = -m(7) * t833 + t852 * mrSges(7,1) - t853 * mrSges(7,2) + t902 * t880 - t903 * t881;
t977 = -m(6) * t835 + t959 * mrSges(6,3) + t960 * t907 - t983;
t1011 = Ifges(5,1) + Ifges(6,2);
t994 = -t1001 * t916 + t1010 * t960 + t1011 * t917;
t1012 = -mrSges(5,2) * t842 - mrSges(6,3) * t835 - pkin(4) * t814 - pkin(9) * t815 - t965 * t817 + t969 * t818 + qJ(5) * (-t916 * t894 + t977) + mrSges(6,2) * t837 + mrSges(5,1) * t841 + t1008 * t959 + t1007 * t917 + t1010 * t878 + (-qJ(5) * mrSges(6,1) - t1009) * t877 + t994 * t916;
t893 = mrSges(5,1) * t916 + mrSges(5,2) * t917;
t904 = -mrSges(5,2) * t960 - mrSges(5,3) * t916;
t811 = m(5) * t841 - t878 * mrSges(5,3) - t917 * t893 + (t904 - t906) * t960 + (mrSges(5,1) - mrSges(6,2)) * t959 + t981;
t905 = mrSges(5,1) * t960 - mrSges(5,3) * t917;
t821 = m(5) * t842 - t959 * mrSges(5,2) - t960 * t905 + (-t893 - t894) * t916 + (-mrSges(5,3) - mrSges(6,1)) * t877 + t977;
t805 = t1003 * t811 + t966 * t821;
t920 = -mrSges(4,1) * t935 + mrSges(4,2) * t936;
t927 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t935;
t803 = m(4) * t860 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t924 + qJD(2) * t927 - t920 * t936 + t805;
t928 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t936;
t986 = t1003 * t821 - t811 * t966;
t804 = m(4) * t861 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t923 - qJD(2) * t928 + t920 * t935 + t986;
t798 = t964 * t803 + t963 * t804;
t914 = Ifges(4,4) * t936 + Ifges(4,2) * t935 + Ifges(4,6) * qJD(2);
t915 = Ifges(4,1) * t936 + Ifges(4,4) * t935 + Ifges(4,5) * qJD(2);
t925 = -t970 * g(3) - t998;
t938 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t967 + Ifges(3,2) * t970) * qJD(1);
t939 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t967 + Ifges(3,4) * t970) * qJD(1);
t1006 = mrSges(3,1) * t925 + mrSges(4,1) * t860 - mrSges(3,2) * t926 - mrSges(4,2) * t861 + Ifges(3,5) * t946 + Ifges(4,5) * t924 + Ifges(3,6) * t947 + Ifges(4,6) * t923 + pkin(2) * t798 + pkin(3) * t805 + (t938 * t967 - t939 * t970) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t936 * t914 - t935 * t915 + t1012;
t945 = (-mrSges(3,1) * t970 + mrSges(3,2) * t967) * qJD(1);
t992 = qJD(1) * t970;
t951 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t992;
t796 = m(3) * t925 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t946 + qJD(2) * t951 - t945 * t993 + t798;
t950 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t993;
t987 = -t803 * t963 + t964 * t804;
t797 = m(3) * t926 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t947 - qJD(2) * t950 + t945 * t992 + t987;
t988 = -t796 * t967 + t970 * t797;
t789 = m(2) * t953 - mrSges(2,1) * t972 - qJDD(1) * mrSges(2,2) + t988;
t839 = t877 * pkin(4) + t975;
t996 = -t965 * t826 + t969 * t827;
t812 = m(6) * t839 - t877 * mrSges(6,2) - t878 * mrSges(6,3) - t916 * t906 - t917 * t907 + t996;
t978 = m(5) * t859 + t877 * mrSges(5,1) + t878 * mrSges(5,2) + t916 * t904 + t917 * t905 + t812;
t809 = m(4) * t901 - t923 * mrSges(4,1) + t924 * mrSges(4,2) - t935 * t927 + t936 * t928 + t978;
t940 = -t972 * pkin(7) + t984;
t974 = -m(3) * t940 + t947 * mrSges(3,1) - t946 * mrSges(3,2) - t950 * t993 + t951 * t992 - t809;
t807 = m(2) * t952 + qJDD(1) * mrSges(2,1) - t972 * mrSges(2,2) + t974;
t997 = t968 * t789 + t971 * t807;
t791 = t970 * t796 + t967 * t797;
t995 = -t1008 * t960 + t1009 * t916 - t1010 * t917;
t989 = t971 * t789 - t807 * t968;
t792 = -mrSges(5,1) * t859 - mrSges(6,1) * t835 + mrSges(6,2) * t839 + mrSges(5,3) * t842 - pkin(4) * t812 - pkin(5) * t983 - pkin(9) * t996 + t1001 * t878 + t1009 * t959 + t1013 * t877 - t969 * t817 - t965 * t818 + t995 * t917 + t994 * t960;
t979 = mrSges(7,1) * t829 - mrSges(7,2) * t830 + Ifges(7,5) * t853 + Ifges(7,6) * t852 + Ifges(7,3) * t876 + t903 * t855 - t902 * t856;
t799 = mrSges(6,1) * t837 + mrSges(5,2) * t859 - mrSges(5,3) * t841 - mrSges(6,3) * t839 + pkin(5) * t815 - qJ(5) * t812 - t1001 * t877 - t1007 * t960 + t1010 * t959 + t1011 * t878 + t995 * t916 + t979;
t913 = Ifges(4,5) * t936 + Ifges(4,6) * t935 + Ifges(4,3) * qJD(2);
t785 = -mrSges(4,1) * t901 + mrSges(4,3) * t861 + Ifges(4,4) * t924 + Ifges(4,2) * t923 + Ifges(4,6) * qJDD(2) - pkin(3) * t978 + pkin(8) * t986 + qJD(2) * t915 + t1003 * t792 + t966 * t799 - t936 * t913;
t786 = mrSges(4,2) * t901 - mrSges(4,3) * t860 + Ifges(4,1) * t924 + Ifges(4,4) * t923 + Ifges(4,5) * qJDD(2) - pkin(8) * t805 - qJD(2) * t914 + t1003 * t799 - t966 * t792 + t935 * t913;
t937 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t967 + Ifges(3,6) * t970) * qJD(1);
t781 = -mrSges(3,1) * t940 + mrSges(3,3) * t926 + Ifges(3,4) * t946 + Ifges(3,2) * t947 + Ifges(3,6) * qJDD(2) - pkin(2) * t809 + qJ(3) * t987 + qJD(2) * t939 + t964 * t785 + t963 * t786 - t937 * t993;
t783 = mrSges(3,2) * t940 - mrSges(3,3) * t925 + Ifges(3,1) * t946 + Ifges(3,4) * t947 + Ifges(3,5) * qJDD(2) - qJ(3) * t798 - qJD(2) * t938 - t785 * t963 + t786 * t964 + t937 * t992;
t982 = mrSges(2,1) * t952 - mrSges(2,2) * t953 + Ifges(2,3) * qJDD(1) + pkin(1) * t974 + pkin(7) * t988 + t970 * t781 + t967 * t783;
t784 = mrSges(2,1) * g(3) + mrSges(2,3) * t953 + t972 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t791 - t1006;
t779 = -mrSges(2,2) * g(3) - mrSges(2,3) * t952 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t972 - pkin(7) * t791 - t781 * t967 + t783 * t970;
t1 = [-m(1) * g(1) + t989; -m(1) * g(2) + t997; (-m(1) - m(2)) * g(3) + t791; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t997 + t971 * t779 - t968 * t784; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t989 + t968 * t779 + t971 * t784; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t982; t982; t1006; t809; t1012; t814; t979;];
tauJB  = t1;
