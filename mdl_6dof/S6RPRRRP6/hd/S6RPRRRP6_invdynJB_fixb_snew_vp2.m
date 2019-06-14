% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:35:09
% EndTime: 2019-05-06 01:35:24
% DurationCPUTime: 14.96s
% Computational Cost: add. (213037->344), mult. (503343->416), div. (0->0), fcn. (381729->10), ass. (0->148)
t992 = Ifges(6,4) + Ifges(7,4);
t1001 = Ifges(6,2) + Ifges(7,2);
t997 = Ifges(6,6) + Ifges(7,6);
t948 = sin(pkin(10));
t949 = cos(pkin(10));
t952 = sin(qJ(3));
t956 = cos(qJ(3));
t968 = t948 * t956 + t949 * t952;
t984 = t949 * qJD(1);
t985 = t948 * qJD(1);
t929 = -t952 * t985 + t956 * t984;
t986 = t929 * qJD(3);
t916 = t968 * qJDD(1) + t986;
t930 = t968 * qJD(1);
t951 = sin(qJ(4));
t955 = cos(qJ(4));
t921 = t951 * qJD(3) + t955 * t930;
t887 = -t921 * qJD(4) + t955 * qJDD(3) - t951 * t916;
t920 = t955 * qJD(3) - t951 * t930;
t888 = t920 * qJD(4) + t951 * qJDD(3) + t955 * t916;
t950 = sin(qJ(5));
t954 = cos(qJ(5));
t890 = t954 * t920 - t950 * t921;
t852 = t890 * qJD(5) + t950 * t887 + t954 * t888;
t891 = t950 * t920 + t954 * t921;
t871 = -t890 * mrSges(7,1) + t891 * mrSges(7,2);
t953 = sin(qJ(1));
t957 = cos(qJ(1));
t936 = -t957 * g(1) - t953 * g(2);
t959 = qJD(1) ^ 2;
t931 = -t959 * pkin(1) + qJDD(1) * qJ(2) + t936;
t983 = qJD(1) * qJD(2);
t978 = -t949 * g(3) - 0.2e1 * t948 * t983;
t993 = pkin(2) * t949;
t900 = (-pkin(7) * qJDD(1) + t959 * t993 - t931) * t948 + t978;
t918 = -t948 * g(3) + (t931 + 0.2e1 * t983) * t949;
t981 = qJDD(1) * t949;
t946 = t949 ^ 2;
t990 = t946 * t959;
t901 = -pkin(2) * t990 + pkin(7) * t981 + t918;
t875 = t952 * t900 + t956 * t901;
t913 = -t929 * pkin(3) - t930 * pkin(8);
t958 = qJD(3) ^ 2;
t858 = -t958 * pkin(3) + qJDD(3) * pkin(8) + t929 * t913 + t875;
t945 = t948 ^ 2;
t935 = t953 * g(1) - t957 * g(2);
t972 = qJDD(2) - t935;
t914 = (-pkin(1) - t993) * qJDD(1) + (-qJ(2) + (-t945 - t946) * pkin(7)) * t959 + t972;
t926 = t930 * qJD(3);
t982 = qJDD(1) * t948;
t915 = -t952 * t982 + t956 * t981 - t926;
t861 = (-t916 - t986) * pkin(8) + (-t915 + t926) * pkin(3) + t914;
t841 = -t951 * t858 + t955 * t861;
t912 = qJDD(4) - t915;
t927 = qJD(4) - t929;
t837 = (t920 * t927 - t888) * pkin(9) + (t920 * t921 + t912) * pkin(4) + t841;
t842 = t955 * t858 + t951 * t861;
t899 = t927 * pkin(4) - t921 * pkin(9);
t919 = t920 ^ 2;
t839 = -t919 * pkin(4) + t887 * pkin(9) - t927 * t899 + t842;
t831 = t954 * t837 - t950 * t839;
t909 = qJDD(5) + t912;
t925 = qJD(5) + t927;
t826 = -0.2e1 * qJD(6) * t891 + (t890 * t925 - t852) * qJ(6) + (t890 * t891 + t909) * pkin(5) + t831;
t876 = -t925 * mrSges(7,2) + t890 * mrSges(7,3);
t980 = m(7) * t826 + t909 * mrSges(7,1) + t925 * t876;
t823 = -t852 * mrSges(7,3) - t891 * t871 + t980;
t832 = t950 * t837 + t954 * t839;
t851 = -t891 * qJD(5) + t954 * t887 - t950 * t888;
t878 = t925 * pkin(5) - t891 * qJ(6);
t889 = t890 ^ 2;
t829 = -t889 * pkin(5) + t851 * qJ(6) + 0.2e1 * qJD(6) * t890 - t925 * t878 + t832;
t998 = Ifges(6,5) + Ifges(7,5);
t999 = Ifges(6,1) + Ifges(7,1);
t987 = -t992 * t890 - t891 * t999 - t998 * t925;
t995 = t1001 * t890 + t992 * t891 + t997 * t925;
t996 = Ifges(6,3) + Ifges(7,3);
t1000 = mrSges(6,1) * t831 + mrSges(7,1) * t826 - mrSges(6,2) * t832 - mrSges(7,2) * t829 + pkin(5) * t823 + t851 * t997 + t852 * t998 + t987 * t890 + t891 * t995 + t909 * t996;
t872 = -t890 * mrSges(6,1) + t891 * mrSges(6,2);
t877 = -t925 * mrSges(6,2) + t890 * mrSges(6,3);
t815 = m(6) * t831 + t909 * mrSges(6,1) + t925 * t877 + (-t871 - t872) * t891 + (-mrSges(6,3) - mrSges(7,3)) * t852 + t980;
t879 = t925 * mrSges(7,1) - t891 * mrSges(7,3);
t880 = t925 * mrSges(6,1) - t891 * mrSges(6,3);
t979 = m(7) * t829 + t851 * mrSges(7,3) + t890 * t871;
t818 = m(6) * t832 + t851 * mrSges(6,3) + t890 * t872 + (-t879 - t880) * t925 + (-mrSges(6,2) - mrSges(7,2)) * t909 + t979;
t813 = t954 * t815 + t950 * t818;
t882 = Ifges(5,4) * t921 + Ifges(5,2) * t920 + Ifges(5,6) * t927;
t883 = Ifges(5,1) * t921 + Ifges(5,4) * t920 + Ifges(5,5) * t927;
t994 = mrSges(5,1) * t841 - mrSges(5,2) * t842 + Ifges(5,5) * t888 + Ifges(5,6) * t887 + Ifges(5,3) * t912 + pkin(4) * t813 + t921 * t882 - t920 * t883 + t1000;
t991 = mrSges(3,2) * t948;
t892 = -t920 * mrSges(5,1) + t921 * mrSges(5,2);
t895 = -t927 * mrSges(5,2) + t920 * mrSges(5,3);
t810 = m(5) * t841 + t912 * mrSges(5,1) - t888 * mrSges(5,3) - t921 * t892 + t927 * t895 + t813;
t896 = t927 * mrSges(5,1) - t921 * mrSges(5,3);
t973 = -t950 * t815 + t954 * t818;
t811 = m(5) * t842 - t912 * mrSges(5,2) + t887 * mrSges(5,3) + t920 * t892 - t927 * t896 + t973;
t805 = -t951 * t810 + t955 * t811;
t910 = -t929 * mrSges(4,1) + t930 * mrSges(4,2);
t923 = qJD(3) * mrSges(4,1) - t930 * mrSges(4,3);
t803 = m(4) * t875 - qJDD(3) * mrSges(4,2) + t915 * mrSges(4,3) - qJD(3) * t923 + t929 * t910 + t805;
t874 = t956 * t900 - t952 * t901;
t857 = -qJDD(3) * pkin(3) - t958 * pkin(8) + t930 * t913 - t874;
t840 = -t887 * pkin(4) - t919 * pkin(9) + t921 * t899 + t857;
t834 = -t851 * pkin(5) - t889 * qJ(6) + t891 * t878 + qJDD(6) + t840;
t827 = m(7) * t834 - t851 * mrSges(7,1) + t852 * mrSges(7,2) - t890 * t876 + t891 * t879;
t964 = m(6) * t840 - t851 * mrSges(6,1) + t852 * mrSges(6,2) - t890 * t877 + t891 * t880 + t827;
t821 = -m(5) * t857 + t887 * mrSges(5,1) - t888 * mrSges(5,2) + t920 * t895 - t921 * t896 - t964;
t922 = -qJD(3) * mrSges(4,2) + t929 * mrSges(4,3);
t820 = m(4) * t874 + qJDD(3) * mrSges(4,1) - t916 * mrSges(4,3) + qJD(3) * t922 - t930 * t910 + t821;
t796 = t952 * t803 + t956 * t820;
t917 = -t948 * t931 + t978;
t967 = mrSges(3,3) * qJDD(1) + t959 * (-mrSges(3,1) * t949 + t991);
t794 = m(3) * t917 - t967 * t948 + t796;
t974 = t956 * t803 - t952 * t820;
t795 = m(3) * t918 + t967 * t949 + t974;
t975 = -t948 * t794 + t949 * t795;
t785 = m(2) * t936 - t959 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t975;
t928 = -qJDD(1) * pkin(1) - t959 * qJ(2) + t972;
t804 = t955 * t810 + t951 * t811;
t965 = m(4) * t914 - t915 * mrSges(4,1) + t916 * mrSges(4,2) - t929 * t922 + t930 * t923 + t804;
t962 = -m(3) * t928 + mrSges(3,1) * t981 - t965 + (t945 * t959 + t990) * mrSges(3,3);
t798 = -t959 * mrSges(2,2) + m(2) * t935 + (mrSges(2,1) - t991) * qJDD(1) + t962;
t989 = t953 * t785 + t957 * t798;
t787 = t949 * t794 + t948 * t795;
t988 = -t890 * t997 - t891 * t998 - t925 * t996;
t976 = t957 * t785 - t953 * t798;
t971 = Ifges(3,1) * t948 + Ifges(3,4) * t949;
t970 = Ifges(3,4) * t948 + Ifges(3,2) * t949;
t969 = Ifges(3,5) * t948 + Ifges(3,6) * t949;
t806 = -mrSges(6,1) * t840 + mrSges(6,3) * t832 - mrSges(7,1) * t834 + mrSges(7,3) * t829 - pkin(5) * t827 + qJ(6) * t979 + (-qJ(6) * t879 - t987) * t925 + (-qJ(6) * mrSges(7,2) + t997) * t909 + t988 * t891 + t992 * t852 + t1001 * t851;
t812 = mrSges(6,2) * t840 + mrSges(7,2) * t834 - mrSges(6,3) * t831 - mrSges(7,3) * t826 - qJ(6) * t823 + t992 * t851 + t852 * t999 - t988 * t890 + t998 * t909 - t995 * t925;
t881 = Ifges(5,5) * t921 + Ifges(5,6) * t920 + Ifges(5,3) * t927;
t789 = -mrSges(5,1) * t857 + mrSges(5,3) * t842 + Ifges(5,4) * t888 + Ifges(5,2) * t887 + Ifges(5,6) * t912 - pkin(4) * t964 + pkin(9) * t973 + t954 * t806 + t950 * t812 - t921 * t881 + t927 * t883;
t790 = mrSges(5,2) * t857 - mrSges(5,3) * t841 + Ifges(5,1) * t888 + Ifges(5,4) * t887 + Ifges(5,5) * t912 - pkin(9) * t813 - t950 * t806 + t954 * t812 + t920 * t881 - t927 * t882;
t902 = Ifges(4,5) * t930 + Ifges(4,6) * t929 + Ifges(4,3) * qJD(3);
t903 = Ifges(4,4) * t930 + Ifges(4,2) * t929 + Ifges(4,6) * qJD(3);
t782 = mrSges(4,2) * t914 - mrSges(4,3) * t874 + Ifges(4,1) * t916 + Ifges(4,4) * t915 + Ifges(4,5) * qJDD(3) - pkin(8) * t804 - qJD(3) * t903 - t951 * t789 + t955 * t790 + t929 * t902;
t904 = Ifges(4,1) * t930 + Ifges(4,4) * t929 + Ifges(4,5) * qJD(3);
t788 = -mrSges(4,1) * t914 + mrSges(4,3) * t875 + Ifges(4,4) * t916 + Ifges(4,2) * t915 + Ifges(4,6) * qJDD(3) - pkin(3) * t804 + qJD(3) * t904 - t930 * t902 - t994;
t933 = t969 * qJD(1);
t778 = -mrSges(3,1) * t928 + mrSges(3,3) * t918 - pkin(2) * t965 + pkin(7) * t974 + t970 * qJDD(1) + t952 * t782 + t956 * t788 - t933 * t985;
t781 = mrSges(3,2) * t928 - mrSges(3,3) * t917 - pkin(7) * t796 + t971 * qJDD(1) + t956 * t782 - t952 * t788 + t933 * t984;
t800 = mrSges(3,2) * t982 - t962;
t966 = mrSges(2,1) * t935 - mrSges(2,2) * t936 + Ifges(2,3) * qJDD(1) - pkin(1) * t800 + qJ(2) * t975 + t949 * t778 + t948 * t781;
t961 = mrSges(4,1) * t874 - mrSges(4,2) * t875 + Ifges(4,5) * t916 + Ifges(4,6) * t915 + Ifges(4,3) * qJDD(3) + pkin(3) * t821 + pkin(8) * t805 + t955 * t789 + t951 * t790 + t930 * t903 - t929 * t904;
t779 = mrSges(2,3) * t936 - mrSges(3,1) * t917 + mrSges(3,2) * t918 + mrSges(2,1) * g(3) - pkin(1) * t787 - t961 - pkin(2) * t796 + (Ifges(2,6) - t969) * qJDD(1) + (-t948 * t970 + t949 * t971 + Ifges(2,5)) * t959;
t776 = -mrSges(2,2) * g(3) - mrSges(2,3) * t935 + Ifges(2,5) * qJDD(1) - t959 * Ifges(2,6) - qJ(2) * t787 - t948 * t778 + t949 * t781;
t1 = [-m(1) * g(1) + t976; -m(1) * g(2) + t989; (-m(1) - m(2)) * g(3) + t787; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t989 + t957 * t776 - t953 * t779; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t976 + t953 * t776 + t957 * t779; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t966; t966; t800; t961; t994; t1000; t827;];
tauJB  = t1;
