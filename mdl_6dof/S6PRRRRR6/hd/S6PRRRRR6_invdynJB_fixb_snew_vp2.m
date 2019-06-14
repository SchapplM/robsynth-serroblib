% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRR6
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 12:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_invdynJB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 12:37:23
% EndTime: 2019-05-05 12:39:13
% DurationCPUTime: 113.53s
% Computational Cost: add. (2155472->371), mult. (5071134->504), div. (0->0), fcn. (4283042->18), ass. (0->175)
t963 = cos(pkin(8));
t1009 = pkin(11) * t963;
t960 = sin(pkin(7));
t1011 = pkin(10) * t960;
t961 = sin(pkin(6));
t970 = sin(qJ(2));
t1004 = t961 * t970;
t958 = sin(pkin(14));
t962 = cos(pkin(14));
t950 = t958 * g(1) - t962 * g(2);
t951 = -t962 * g(1) - t958 * g(2);
t957 = -g(3) + qJDD(1);
t975 = cos(qJ(2));
t965 = cos(pkin(6));
t998 = t965 * t970;
t925 = t957 * t1004 + t950 * t998 + t975 * t951;
t976 = qJD(2) ^ 2;
t922 = -t976 * pkin(2) + qJDD(2) * t1011 + t925;
t964 = cos(pkin(7));
t955 = t964 * qJD(2) + qJD(3);
t959 = sin(pkin(8));
t974 = cos(qJ(3));
t994 = qJD(2) * t960;
t990 = t974 * t994;
t937 = (t955 * t959 + t963 * t990) * pkin(11);
t1010 = pkin(11) * t959;
t969 = sin(qJ(3));
t941 = (-pkin(3) * t974 - t969 * t1010) * t994;
t993 = qJD(2) * qJD(3);
t945 = (qJDD(2) * t969 + t974 * t993) * t960;
t954 = t964 * qJDD(2) + qJDD(3);
t1005 = t960 * t974;
t1003 = t961 * t975;
t997 = t965 * t975;
t924 = t957 * t1003 + t950 * t997 - t970 * t951;
t921 = qJDD(2) * pkin(2) + t976 * t1011 + t924;
t940 = -t961 * t950 + t965 * t957;
t999 = t964 * t974;
t995 = t940 * t1005 + t921 * t999;
t880 = -t945 * t1009 + t954 * pkin(3) + t955 * t937 + (-t941 * t994 - t922) * t969 + t995;
t1000 = t964 * t969;
t1006 = t960 * t969;
t898 = t921 * t1000 + t940 * t1006 + t974 * t922;
t991 = t969 * t994;
t939 = t955 * pkin(3) - t991 * t1009;
t946 = (qJDD(2) * t974 - t969 * t993) * t960;
t985 = t946 * t963 + t954 * t959;
t881 = t985 * pkin(11) - t955 * t939 + t941 * t990 + t898;
t935 = t964 * t940;
t894 = -t945 * t1010 - t946 * pkin(3) + t935 + (-t921 + (-t937 * t974 + t939 * t969) * qJD(2)) * t960;
t968 = sin(qJ(4));
t973 = cos(qJ(4));
t866 = -t968 * t881 + (t880 * t963 + t894 * t959) * t973;
t1001 = t963 * t974;
t1008 = t959 * t968;
t928 = t955 * t1008 + (t968 * t1001 + t969 * t973) * t994;
t908 = -t928 * qJD(4) - t968 * t945 + t985 * t973;
t1007 = t959 * t973;
t927 = (t973 * t1001 - t968 * t969) * t994 + t955 * t1007;
t1002 = t963 * t968;
t867 = t880 * t1002 + t894 * t1008 + t973 * t881;
t910 = -t927 * mrSges(5,1) + t928 * mrSges(5,2);
t938 = t963 * t955 - t959 * t990 + qJD(4);
t917 = t938 * mrSges(5,1) - t928 * mrSges(5,3);
t929 = -t959 * t946 + t963 * t954 + qJDD(4);
t911 = -t927 * pkin(4) - t928 * pkin(12);
t936 = t938 ^ 2;
t863 = -t936 * pkin(4) + t929 * pkin(12) + t927 * t911 + t867;
t871 = -t959 * t880 + t963 * t894;
t909 = t927 * qJD(4) + t973 * t945 + t985 * t968;
t865 = (-t927 * t938 - t909) * pkin(12) + (t928 * t938 - t908) * pkin(4) + t871;
t967 = sin(qJ(5));
t972 = cos(qJ(5));
t859 = t972 * t863 + t967 * t865;
t914 = -t967 * t928 + t972 * t938;
t915 = t972 * t928 + t967 * t938;
t896 = -t914 * pkin(5) - t915 * pkin(13);
t907 = qJDD(5) - t908;
t926 = qJD(5) - t927;
t923 = t926 ^ 2;
t857 = -t923 * pkin(5) + t907 * pkin(13) + t914 * t896 + t859;
t862 = -t929 * pkin(4) - t936 * pkin(12) + t928 * t911 - t866;
t884 = -t915 * qJD(5) - t967 * t909 + t972 * t929;
t885 = t914 * qJD(5) + t972 * t909 + t967 * t929;
t860 = (-t914 * t926 - t885) * pkin(13) + (t915 * t926 - t884) * pkin(5) + t862;
t966 = sin(qJ(6));
t971 = cos(qJ(6));
t854 = -t966 * t857 + t971 * t860;
t900 = -t966 * t915 + t971 * t926;
t870 = t900 * qJD(6) + t971 * t885 + t966 * t907;
t901 = t971 * t915 + t966 * t926;
t876 = -t900 * mrSges(7,1) + t901 * mrSges(7,2);
t883 = qJDD(6) - t884;
t913 = qJD(6) - t914;
t886 = -t913 * mrSges(7,2) + t900 * mrSges(7,3);
t851 = m(7) * t854 + t883 * mrSges(7,1) - t870 * mrSges(7,3) - t901 * t876 + t913 * t886;
t855 = t971 * t857 + t966 * t860;
t869 = -t901 * qJD(6) - t966 * t885 + t971 * t907;
t887 = t913 * mrSges(7,1) - t901 * mrSges(7,3);
t852 = m(7) * t855 - t883 * mrSges(7,2) + t869 * mrSges(7,3) + t900 * t876 - t913 * t887;
t845 = -t966 * t851 + t971 * t852;
t895 = -t914 * mrSges(6,1) + t915 * mrSges(6,2);
t903 = t926 * mrSges(6,1) - t915 * mrSges(6,3);
t843 = m(6) * t859 - t907 * mrSges(6,2) + t884 * mrSges(6,3) + t914 * t895 - t926 * t903 + t845;
t858 = -t967 * t863 + t972 * t865;
t856 = -t907 * pkin(5) - t923 * pkin(13) + t915 * t896 - t858;
t853 = -m(7) * t856 + t869 * mrSges(7,1) - t870 * mrSges(7,2) + t900 * t886 - t901 * t887;
t902 = -t926 * mrSges(6,2) + t914 * mrSges(6,3);
t849 = m(6) * t858 + t907 * mrSges(6,1) - t885 * mrSges(6,3) - t915 * t895 + t926 * t902 + t853;
t988 = t972 * t843 - t967 * t849;
t834 = m(5) * t867 - t929 * mrSges(5,2) + t908 * mrSges(5,3) + t927 * t910 - t938 * t917 + t988;
t837 = t967 * t843 + t972 * t849;
t916 = -t938 * mrSges(5,2) + t927 * mrSges(5,3);
t836 = m(5) * t871 - t908 * mrSges(5,1) + t909 * mrSges(5,2) - t927 * t916 + t928 * t917 + t837;
t844 = t971 * t851 + t966 * t852;
t979 = -m(6) * t862 + t884 * mrSges(6,1) - t885 * mrSges(6,2) + t914 * t902 - t915 * t903 - t844;
t840 = m(5) * t866 + t929 * mrSges(5,1) - t909 * mrSges(5,3) - t928 * t910 + t938 * t916 + t979;
t823 = t963 * t973 * t840 + t834 * t1002 - t959 * t836;
t897 = -t969 * t922 + t995;
t943 = -t955 * mrSges(4,2) + mrSges(4,3) * t990;
t944 = (-mrSges(4,1) * t974 + mrSges(4,2) * t969) * t994;
t819 = m(4) * t897 + t954 * mrSges(4,1) - t945 * mrSges(4,3) + t955 * t943 - t944 * t991 + t823;
t822 = t840 * t1007 + t834 * t1008 + t963 * t836;
t912 = -t960 * t921 + t935;
t942 = t955 * mrSges(4,1) - mrSges(4,3) * t991;
t821 = m(4) * t912 - t946 * mrSges(4,1) + t945 * mrSges(4,2) + (t942 * t969 - t943 * t974) * t994 + t822;
t828 = t973 * t834 - t968 * t840;
t827 = m(4) * t898 - t954 * mrSges(4,2) + t946 * mrSges(4,3) - t955 * t942 + t944 * t990 + t828;
t808 = t827 * t1000 + t819 * t999 - t960 * t821;
t804 = m(3) * t924 + qJDD(2) * mrSges(3,1) - t976 * mrSges(3,2) + t808;
t807 = t819 * t1005 + t827 * t1006 + t964 * t821;
t806 = m(3) * t940 + t807;
t813 = -t969 * t819 + t974 * t827;
t812 = m(3) * t925 - t976 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t813;
t794 = t804 * t997 - t961 * t806 + t812 * t998;
t792 = m(2) * t950 + t794;
t801 = -t970 * t804 + t975 * t812;
t800 = m(2) * t951 + t801;
t996 = t962 * t792 + t958 * t800;
t793 = t804 * t1003 + t812 * t1004 + t965 * t806;
t989 = -t958 * t792 + t962 * t800;
t987 = m(2) * t957 + t793;
t872 = Ifges(7,5) * t901 + Ifges(7,6) * t900 + Ifges(7,3) * t913;
t874 = Ifges(7,1) * t901 + Ifges(7,4) * t900 + Ifges(7,5) * t913;
t846 = -mrSges(7,1) * t856 + mrSges(7,3) * t855 + Ifges(7,4) * t870 + Ifges(7,2) * t869 + Ifges(7,6) * t883 - t901 * t872 + t913 * t874;
t873 = Ifges(7,4) * t901 + Ifges(7,2) * t900 + Ifges(7,6) * t913;
t847 = mrSges(7,2) * t856 - mrSges(7,3) * t854 + Ifges(7,1) * t870 + Ifges(7,4) * t869 + Ifges(7,5) * t883 + t900 * t872 - t913 * t873;
t890 = Ifges(6,5) * t915 + Ifges(6,6) * t914 + Ifges(6,3) * t926;
t891 = Ifges(6,4) * t915 + Ifges(6,2) * t914 + Ifges(6,6) * t926;
t829 = mrSges(6,2) * t862 - mrSges(6,3) * t858 + Ifges(6,1) * t885 + Ifges(6,4) * t884 + Ifges(6,5) * t907 - pkin(13) * t844 - t966 * t846 + t971 * t847 + t914 * t890 - t926 * t891;
t892 = Ifges(6,1) * t915 + Ifges(6,4) * t914 + Ifges(6,5) * t926;
t978 = mrSges(7,1) * t854 - mrSges(7,2) * t855 + Ifges(7,5) * t870 + Ifges(7,6) * t869 + Ifges(7,3) * t883 + t901 * t873 - t900 * t874;
t830 = -mrSges(6,1) * t862 + mrSges(6,3) * t859 + Ifges(6,4) * t885 + Ifges(6,2) * t884 + Ifges(6,6) * t907 - pkin(5) * t844 - t915 * t890 + t926 * t892 - t978;
t905 = Ifges(5,4) * t928 + Ifges(5,2) * t927 + Ifges(5,6) * t938;
t906 = Ifges(5,1) * t928 + Ifges(5,4) * t927 + Ifges(5,5) * t938;
t814 = mrSges(5,1) * t866 - mrSges(5,2) * t867 + Ifges(5,5) * t909 + Ifges(5,6) * t908 + Ifges(5,3) * t929 + pkin(4) * t979 + pkin(12) * t988 + t967 * t829 + t972 * t830 + t928 * t905 - t927 * t906;
t933 = Ifges(4,6) * t955 + (Ifges(4,4) * t969 + Ifges(4,2) * t974) * t994;
t934 = Ifges(4,5) * t955 + (Ifges(4,1) * t969 + Ifges(4,4) * t974) * t994;
t904 = Ifges(5,5) * t928 + Ifges(5,6) * t927 + Ifges(5,3) * t938;
t815 = mrSges(5,2) * t871 - mrSges(5,3) * t866 + Ifges(5,1) * t909 + Ifges(5,4) * t908 + Ifges(5,5) * t929 - pkin(12) * t837 + t972 * t829 - t967 * t830 + t927 * t904 - t938 * t905;
t977 = mrSges(6,1) * t858 - mrSges(6,2) * t859 + Ifges(6,5) * t885 + Ifges(6,6) * t884 + Ifges(6,3) * t907 + pkin(5) * t853 + pkin(13) * t845 + t971 * t846 + t966 * t847 + t915 * t891 - t914 * t892;
t816 = -mrSges(5,1) * t871 + mrSges(5,3) * t867 + Ifges(5,4) * t909 + Ifges(5,2) * t908 + Ifges(5,6) * t929 - pkin(4) * t837 - t928 * t904 + t938 * t906 - t977;
t980 = pkin(11) * t828 + t815 * t968 + t816 * t973;
t795 = mrSges(4,1) * t897 - mrSges(4,2) * t898 + Ifges(4,5) * t945 + Ifges(4,6) * t946 + Ifges(4,3) * t954 + pkin(3) * t823 + t963 * t814 + (t933 * t969 - t934 * t974) * t994 + t980 * t959;
t932 = Ifges(4,3) * t955 + (Ifges(4,5) * t969 + Ifges(4,6) * t974) * t994;
t796 = -mrSges(4,1) * t912 + mrSges(4,3) * t898 + Ifges(4,4) * t945 + Ifges(4,2) * t946 + Ifges(4,6) * t954 - pkin(3) * t822 - t959 * t814 - t932 * t991 + t955 * t934 + t980 * t963;
t797 = t932 * t990 + mrSges(4,2) * t912 - mrSges(4,3) * t897 + Ifges(4,1) * t945 + Ifges(4,4) * t946 + Ifges(4,5) * t954 + t973 * t815 - t968 * t816 - t955 * t933 + (-t822 * t959 - t823 * t963) * pkin(11);
t981 = pkin(10) * t813 + t796 * t974 + t797 * t969;
t789 = -mrSges(3,1) * t940 + mrSges(3,3) * t925 + t976 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t807 - t960 * t795 + t981 * t964;
t790 = mrSges(3,2) * t940 - mrSges(3,3) * t924 + Ifges(3,5) * qJDD(2) - t976 * Ifges(3,6) - t969 * t796 + t974 * t797 + (-t807 * t960 - t808 * t964) * pkin(10);
t982 = pkin(9) * t801 + t789 * t975 + t790 * t970;
t788 = mrSges(3,1) * t924 - mrSges(3,2) * t925 + Ifges(3,3) * qJDD(2) + pkin(2) * t808 + t964 * t795 + t981 * t960;
t787 = mrSges(2,2) * t957 - mrSges(2,3) * t950 - t970 * t789 + t975 * t790 + (-t793 * t961 - t794 * t965) * pkin(9);
t786 = -mrSges(2,1) * t957 + mrSges(2,3) * t951 - pkin(1) * t793 - t961 * t788 + t982 * t965;
t1 = [-m(1) * g(1) + t989; -m(1) * g(2) + t996; -m(1) * g(3) + t987; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t996 - t958 * t786 + t962 * t787; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t989 + t962 * t786 + t958 * t787; -mrSges(1,1) * g(2) + mrSges(2,1) * t950 + mrSges(1,2) * g(1) - mrSges(2,2) * t951 + pkin(1) * t794 + t965 * t788 + t982 * t961; t987; t788; t795; t814; t977; t978;];
tauJB  = t1;
