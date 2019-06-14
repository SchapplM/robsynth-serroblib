% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 09:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:55:09
% EndTime: 2019-05-08 08:56:06
% DurationCPUTime: 42.57s
% Computational Cost: add. (716509->391), mult. (1426971->488), div. (0->0), fcn. (1064065->12), ass. (0->157)
t975 = qJD(1) ^ 2;
t1000 = pkin(2) * t975;
t967 = sin(qJ(2));
t973 = cos(qJ(2));
t995 = qJD(1) * qJD(2);
t947 = qJDD(1) * t967 + t973 * t995;
t968 = sin(qJ(1));
t974 = cos(qJ(1));
t954 = -g(1) * t974 - g(2) * t968;
t941 = -pkin(1) * t975 + qJDD(1) * pkin(7) + t954;
t999 = t967 * t941;
t900 = qJDD(2) * pkin(2) - t947 * pkin(8) - t999 + (pkin(8) * t995 + t1000 * t967 - g(3)) * t973;
t927 = -g(3) * t967 + t973 * t941;
t948 = qJDD(1) * t973 - t967 * t995;
t997 = qJD(1) * t967;
t952 = qJD(2) * pkin(2) - pkin(8) * t997;
t962 = t973 ^ 2;
t902 = pkin(8) * t948 - qJD(2) * t952 - t1000 * t962 + t927;
t966 = sin(qJ(3));
t972 = cos(qJ(3));
t880 = t966 * t900 + t972 * t902;
t939 = (t966 * t973 + t967 * t972) * qJD(1);
t911 = -t939 * qJD(3) - t966 * t947 + t948 * t972;
t996 = qJD(1) * t973;
t938 = -t966 * t997 + t972 * t996;
t921 = -mrSges(4,1) * t938 + mrSges(4,2) * t939;
t960 = qJD(2) + qJD(3);
t929 = mrSges(4,1) * t960 - mrSges(4,3) * t939;
t959 = qJDD(2) + qJDD(3);
t912 = qJD(3) * t938 + t947 * t972 + t948 * t966;
t953 = t968 * g(1) - t974 * g(2);
t986 = -qJDD(1) * pkin(1) - t953;
t913 = -t948 * pkin(2) + t952 * t997 + (-pkin(8) * t962 - pkin(7)) * t975 + t986;
t865 = (-t938 * t960 - t912) * pkin(9) + (t939 * t960 - t911) * pkin(3) + t913;
t922 = -pkin(3) * t938 - pkin(9) * t939;
t958 = t960 ^ 2;
t868 = -pkin(3) * t958 + pkin(9) * t959 + t922 * t938 + t880;
t965 = sin(qJ(4));
t971 = cos(qJ(4));
t848 = t971 * t865 - t965 * t868;
t924 = -t939 * t965 + t960 * t971;
t883 = qJD(4) * t924 + t912 * t971 + t959 * t965;
t910 = qJDD(4) - t911;
t925 = t939 * t971 + t960 * t965;
t934 = qJD(4) - t938;
t844 = (t924 * t934 - t883) * pkin(10) + (t924 * t925 + t910) * pkin(4) + t848;
t849 = t965 * t865 + t971 * t868;
t882 = -qJD(4) * t925 - t912 * t965 + t959 * t971;
t916 = pkin(4) * t934 - pkin(10) * t925;
t923 = t924 ^ 2;
t846 = -pkin(4) * t923 + pkin(10) * t882 - t916 * t934 + t849;
t964 = sin(qJ(5));
t970 = cos(qJ(5));
t832 = t970 * t844 - t964 * t846;
t893 = t924 * t970 - t925 * t964;
t861 = qJD(5) * t893 + t882 * t964 + t883 * t970;
t894 = t924 * t964 + t925 * t970;
t905 = qJDD(5) + t910;
t932 = qJD(5) + t934;
t829 = (t893 * t932 - t861) * pkin(11) + (t893 * t894 + t905) * pkin(5) + t832;
t833 = t964 * t844 + t970 * t846;
t860 = -qJD(5) * t894 + t882 * t970 - t883 * t964;
t886 = pkin(5) * t932 - pkin(11) * t894;
t892 = t893 ^ 2;
t830 = -pkin(5) * t892 + pkin(11) * t860 - t886 * t932 + t833;
t963 = sin(qJ(6));
t969 = cos(qJ(6));
t827 = t829 * t969 - t830 * t963;
t875 = t893 * t969 - t894 * t963;
t841 = qJD(6) * t875 + t860 * t963 + t861 * t969;
t876 = t893 * t963 + t894 * t969;
t856 = -mrSges(7,1) * t875 + mrSges(7,2) * t876;
t930 = qJD(6) + t932;
t869 = -mrSges(7,2) * t930 + mrSges(7,3) * t875;
t904 = qJDD(6) + t905;
t822 = m(7) * t827 + mrSges(7,1) * t904 - mrSges(7,3) * t841 - t856 * t876 + t869 * t930;
t828 = t829 * t963 + t830 * t969;
t840 = -qJD(6) * t876 + t860 * t969 - t861 * t963;
t870 = mrSges(7,1) * t930 - mrSges(7,3) * t876;
t823 = m(7) * t828 - mrSges(7,2) * t904 + mrSges(7,3) * t840 + t856 * t875 - t870 * t930;
t814 = t969 * t822 + t963 * t823;
t877 = -mrSges(6,1) * t893 + mrSges(6,2) * t894;
t884 = -mrSges(6,2) * t932 + mrSges(6,3) * t893;
t811 = m(6) * t832 + mrSges(6,1) * t905 - mrSges(6,3) * t861 - t877 * t894 + t884 * t932 + t814;
t885 = mrSges(6,1) * t932 - mrSges(6,3) * t894;
t989 = -t822 * t963 + t969 * t823;
t812 = m(6) * t833 - mrSges(6,2) * t905 + mrSges(6,3) * t860 + t877 * t893 - t885 * t932 + t989;
t807 = t970 * t811 + t964 * t812;
t898 = -mrSges(5,1) * t924 + mrSges(5,2) * t925;
t914 = -mrSges(5,2) * t934 + mrSges(5,3) * t924;
t805 = m(5) * t848 + mrSges(5,1) * t910 - mrSges(5,3) * t883 - t898 * t925 + t914 * t934 + t807;
t915 = mrSges(5,1) * t934 - mrSges(5,3) * t925;
t990 = -t811 * t964 + t970 * t812;
t806 = m(5) * t849 - mrSges(5,2) * t910 + mrSges(5,3) * t882 + t898 * t924 - t915 * t934 + t990;
t991 = -t805 * t965 + t971 * t806;
t796 = m(4) * t880 - mrSges(4,2) * t959 + mrSges(4,3) * t911 + t921 * t938 - t929 * t960 + t991;
t879 = t900 * t972 - t966 * t902;
t928 = -mrSges(4,2) * t960 + mrSges(4,3) * t938;
t867 = -pkin(3) * t959 - pkin(9) * t958 + t939 * t922 - t879;
t850 = -pkin(4) * t882 - pkin(10) * t923 + t925 * t916 + t867;
t835 = -pkin(5) * t860 - pkin(11) * t892 + t886 * t894 + t850;
t988 = m(7) * t835 - t840 * mrSges(7,1) + t841 * mrSges(7,2) - t875 * t869 + t876 * t870;
t981 = m(6) * t850 - t860 * mrSges(6,1) + mrSges(6,2) * t861 - t893 * t884 + t885 * t894 + t988;
t978 = -m(5) * t867 + t882 * mrSges(5,1) - mrSges(5,2) * t883 + t924 * t914 - t915 * t925 - t981;
t818 = m(4) * t879 + mrSges(4,1) * t959 - mrSges(4,3) * t912 - t921 * t939 + t928 * t960 + t978;
t790 = t966 * t796 + t972 * t818;
t926 = -t973 * g(3) - t999;
t936 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t967 + Ifges(3,2) * t973) * qJD(1);
t937 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t967 + Ifges(3,4) * t973) * qJD(1);
t851 = Ifges(7,5) * t876 + Ifges(7,6) * t875 + Ifges(7,3) * t930;
t853 = Ifges(7,1) * t876 + Ifges(7,4) * t875 + Ifges(7,5) * t930;
t815 = -mrSges(7,1) * t835 + mrSges(7,3) * t828 + Ifges(7,4) * t841 + Ifges(7,2) * t840 + Ifges(7,6) * t904 - t851 * t876 + t853 * t930;
t852 = Ifges(7,4) * t876 + Ifges(7,2) * t875 + Ifges(7,6) * t930;
t816 = mrSges(7,2) * t835 - mrSges(7,3) * t827 + Ifges(7,1) * t841 + Ifges(7,4) * t840 + Ifges(7,5) * t904 + t851 * t875 - t852 * t930;
t871 = Ifges(6,5) * t894 + Ifges(6,6) * t893 + Ifges(6,3) * t932;
t873 = Ifges(6,1) * t894 + Ifges(6,4) * t893 + Ifges(6,5) * t932;
t800 = -mrSges(6,1) * t850 + mrSges(6,3) * t833 + Ifges(6,4) * t861 + Ifges(6,2) * t860 + Ifges(6,6) * t905 - pkin(5) * t988 + pkin(11) * t989 + t969 * t815 + t963 * t816 - t894 * t871 + t932 * t873;
t872 = Ifges(6,4) * t894 + Ifges(6,2) * t893 + Ifges(6,6) * t932;
t801 = mrSges(6,2) * t850 - mrSges(6,3) * t832 + Ifges(6,1) * t861 + Ifges(6,4) * t860 + Ifges(6,5) * t905 - pkin(11) * t814 - t815 * t963 + t816 * t969 + t871 * t893 - t872 * t932;
t887 = Ifges(5,5) * t925 + Ifges(5,6) * t924 + Ifges(5,3) * t934;
t889 = Ifges(5,1) * t925 + Ifges(5,4) * t924 + Ifges(5,5) * t934;
t782 = -mrSges(5,1) * t867 + mrSges(5,3) * t849 + Ifges(5,4) * t883 + Ifges(5,2) * t882 + Ifges(5,6) * t910 - pkin(4) * t981 + pkin(10) * t990 + t970 * t800 + t964 * t801 - t925 * t887 + t934 * t889;
t888 = Ifges(5,4) * t925 + Ifges(5,2) * t924 + Ifges(5,6) * t934;
t784 = mrSges(5,2) * t867 - mrSges(5,3) * t848 + Ifges(5,1) * t883 + Ifges(5,4) * t882 + Ifges(5,5) * t910 - pkin(10) * t807 - t800 * t964 + t801 * t970 + t887 * t924 - t888 * t934;
t918 = Ifges(4,4) * t939 + Ifges(4,2) * t938 + Ifges(4,6) * t960;
t919 = Ifges(4,1) * t939 + Ifges(4,4) * t938 + Ifges(4,5) * t960;
t982 = -mrSges(4,1) * t879 + mrSges(4,2) * t880 - Ifges(4,5) * t912 - Ifges(4,6) * t911 - Ifges(4,3) * t959 - pkin(3) * t978 - pkin(9) * t991 - t971 * t782 - t965 * t784 - t939 * t918 + t938 * t919;
t1001 = mrSges(3,1) * t926 - mrSges(3,2) * t927 + Ifges(3,5) * t947 + Ifges(3,6) * t948 + Ifges(3,3) * qJDD(2) + pkin(2) * t790 + (t936 * t967 - t937 * t973) * qJD(1) - t982;
t946 = (-mrSges(3,1) * t973 + mrSges(3,2) * t967) * qJD(1);
t951 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t996;
t788 = m(3) * t926 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t947 + qJD(2) * t951 - t946 * t997 + t790;
t950 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t997;
t992 = t972 * t796 - t818 * t966;
t789 = m(3) * t927 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t948 - qJD(2) * t950 + t946 * t996 + t992;
t993 = -t788 * t967 + t973 * t789;
t777 = m(2) * t954 - mrSges(2,1) * t975 - qJDD(1) * mrSges(2,2) + t993;
t940 = -t975 * pkin(7) + t986;
t798 = t971 * t805 + t965 * t806;
t983 = m(4) * t913 - t911 * mrSges(4,1) + mrSges(4,2) * t912 - t938 * t928 + t929 * t939 + t798;
t979 = -m(3) * t940 + t948 * mrSges(3,1) - mrSges(3,2) * t947 - t950 * t997 + t951 * t996 - t983;
t792 = m(2) * t953 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t975 + t979;
t998 = t968 * t777 + t974 * t792;
t779 = t973 * t788 + t967 * t789;
t994 = t974 * t777 - t792 * t968;
t917 = Ifges(4,5) * t939 + Ifges(4,6) * t938 + Ifges(4,3) * t960;
t774 = mrSges(4,2) * t913 - mrSges(4,3) * t879 + Ifges(4,1) * t912 + Ifges(4,4) * t911 + Ifges(4,5) * t959 - pkin(9) * t798 - t782 * t965 + t784 * t971 + t917 * t938 - t918 * t960;
t984 = -mrSges(7,1) * t827 + mrSges(7,2) * t828 - Ifges(7,5) * t841 - Ifges(7,6) * t840 - Ifges(7,3) * t904 - t876 * t852 + t875 * t853;
t980 = -mrSges(6,1) * t832 + mrSges(6,2) * t833 - Ifges(6,5) * t861 - Ifges(6,6) * t860 - Ifges(6,3) * t905 - pkin(5) * t814 - t894 * t872 + t893 * t873 + t984;
t976 = mrSges(5,1) * t848 - mrSges(5,2) * t849 + Ifges(5,5) * t883 + Ifges(5,6) * t882 + Ifges(5,3) * t910 + pkin(4) * t807 + t925 * t888 - t924 * t889 - t980;
t780 = -mrSges(4,1) * t913 + mrSges(4,3) * t880 + Ifges(4,4) * t912 + Ifges(4,2) * t911 + Ifges(4,6) * t959 - pkin(3) * t798 - t939 * t917 + t960 * t919 - t976;
t935 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t967 + Ifges(3,6) * t973) * qJD(1);
t770 = -mrSges(3,1) * t940 + mrSges(3,3) * t927 + Ifges(3,4) * t947 + Ifges(3,2) * t948 + Ifges(3,6) * qJDD(2) - pkin(2) * t983 + pkin(8) * t992 + qJD(2) * t937 + t966 * t774 + t972 * t780 - t935 * t997;
t773 = mrSges(3,2) * t940 - mrSges(3,3) * t926 + Ifges(3,1) * t947 + Ifges(3,4) * t948 + Ifges(3,5) * qJDD(2) - pkin(8) * t790 - qJD(2) * t936 + t774 * t972 - t780 * t966 + t935 * t996;
t985 = mrSges(2,1) * t953 - mrSges(2,2) * t954 + Ifges(2,3) * qJDD(1) + pkin(1) * t979 + pkin(7) * t993 + t973 * t770 + t967 * t773;
t771 = mrSges(2,1) * g(3) + mrSges(2,3) * t954 + t975 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t779 - t1001;
t768 = -mrSges(2,2) * g(3) - mrSges(2,3) * t953 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t975 - pkin(7) * t779 - t770 * t967 + t773 * t973;
t1 = [-m(1) * g(1) + t994; -m(1) * g(2) + t998; (-m(1) - m(2)) * g(3) + t779; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t998 + t974 * t768 - t968 * t771; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t994 + t968 * t768 + t974 * t771; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t985; t985; t1001; -t982; t976; -t980; -t984;];
tauJB  = t1;
