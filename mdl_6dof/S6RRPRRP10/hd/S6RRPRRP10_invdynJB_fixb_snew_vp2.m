% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP10
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
% Datum: 2019-05-06 18:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:41:14
% EndTime: 2019-05-06 18:41:50
% DurationCPUTime: 30.42s
% Computational Cost: add. (488198->376), mult. (1106110->476), div. (0->0), fcn. (886072->12), ass. (0->156)
t1009 = Ifges(6,1) + Ifges(7,1);
t1001 = Ifges(6,4) - Ifges(7,5);
t1000 = -Ifges(6,5) - Ifges(7,4);
t1008 = Ifges(6,2) + Ifges(7,3);
t999 = Ifges(6,6) - Ifges(7,6);
t1007 = -Ifges(6,3) - Ifges(7,2);
t1005 = cos(qJ(5));
t965 = cos(pkin(6));
t958 = t965 * qJD(1) + qJD(2);
t962 = sin(pkin(11));
t964 = cos(pkin(11));
t968 = sin(qJ(2));
t963 = sin(pkin(6));
t987 = qJD(1) * t963;
t983 = t968 * t987;
t936 = t964 * t958 - t962 * t983;
t937 = t962 * t958 + t964 * t983;
t967 = sin(qJ(4));
t970 = cos(qJ(4));
t917 = t970 * t936 - t967 * t937;
t971 = cos(qJ(2));
t986 = qJD(1) * t971;
t947 = (qJD(2) * t986 + qJDD(1) * t968) * t963;
t957 = t965 * qJDD(1) + qJDD(2);
t923 = -t962 * t947 + t964 * t957;
t924 = t964 * t947 + t962 * t957;
t887 = t917 * qJD(4) + t967 * t923 + t970 * t924;
t918 = t967 * t936 + t970 * t937;
t982 = t963 * t986;
t952 = qJD(4) - t982;
t966 = sin(qJ(5));
t904 = -t1005 * t952 + t966 * t918;
t985 = qJDD(1) * t963;
t948 = -qJD(2) * t983 + t971 * t985;
t940 = qJDD(4) - t948;
t866 = -t904 * qJD(5) + t1005 * t887 + t966 * t940;
t905 = t1005 * t918 + t966 * t952;
t880 = t904 * mrSges(7,1) - t905 * mrSges(7,3);
t945 = (-pkin(2) * t971 - qJ(3) * t968) * t987;
t956 = t958 ^ 2;
t1004 = pkin(8) * t963;
t969 = sin(qJ(1));
t972 = cos(qJ(1));
t953 = t969 * g(1) - t972 * g(2);
t973 = qJD(1) ^ 2;
t943 = qJDD(1) * pkin(1) + t973 * t1004 + t953;
t954 = -t972 * g(1) - t969 * g(2);
t944 = -t973 * pkin(1) + pkin(8) * t985 + t954;
t995 = t965 * t968;
t988 = t943 * t995 + t971 * t944;
t902 = -t956 * pkin(2) + t957 * qJ(3) + (-g(3) * t968 + t945 * t986) * t963 + t988;
t1003 = t965 * g(3);
t903 = -t948 * pkin(2) - t1003 - t947 * qJ(3) + (-t943 + (pkin(2) * t968 - qJ(3) * t971) * t958 * qJD(1)) * t963;
t861 = -0.2e1 * qJD(3) * t937 - t962 * t902 + t964 * t903;
t857 = (-t936 * t982 - t924) * pkin(9) + (t936 * t937 - t948) * pkin(3) + t861;
t862 = 0.2e1 * qJD(3) * t936 + t964 * t902 + t962 * t903;
t925 = -pkin(3) * t982 - t937 * pkin(9);
t934 = t936 ^ 2;
t860 = -t934 * pkin(3) + t923 * pkin(9) + t925 * t982 + t862;
t853 = t967 * t857 + t970 * t860;
t897 = -t917 * pkin(4) - t918 * pkin(10);
t951 = t952 ^ 2;
t851 = -t951 * pkin(4) + t940 * pkin(10) + t917 * t897 + t853;
t994 = t965 * t971;
t996 = t963 * t971;
t913 = -g(3) * t996 + t943 * t994 - t968 * t944;
t901 = -t957 * pkin(2) - t956 * qJ(3) + t945 * t983 + qJDD(3) - t913;
t867 = -t923 * pkin(3) - t934 * pkin(9) + t937 * t925 + t901;
t886 = -t918 * qJD(4) + t970 * t923 - t967 * t924;
t855 = (-t917 * t952 - t887) * pkin(10) + (t918 * t952 - t886) * pkin(4) + t867;
t847 = t1005 * t855 - t966 * t851;
t879 = t904 * pkin(5) - t905 * qJ(6);
t885 = qJDD(5) - t886;
t916 = qJD(5) - t917;
t915 = t916 ^ 2;
t845 = -t885 * pkin(5) - t915 * qJ(6) + t905 * t879 + qJDD(6) - t847;
t888 = -t904 * mrSges(7,2) + t916 * mrSges(7,3);
t978 = -m(7) * t845 + t885 * mrSges(7,1) + t916 * t888;
t841 = t866 * mrSges(7,2) + t905 * t880 - t978;
t848 = t1005 * t851 + t966 * t855;
t844 = -t915 * pkin(5) + t885 * qJ(6) + 0.2e1 * qJD(6) * t916 - t904 * t879 + t848;
t865 = t905 * qJD(5) - t1005 * t940 + t966 * t887;
t891 = -t916 * mrSges(7,1) + t905 * mrSges(7,2);
t984 = m(7) * t844 + t885 * mrSges(7,3) + t916 * t891;
t990 = t1000 * t916 + t1001 * t904 - t1009 * t905;
t991 = -t1001 * t905 + t1008 * t904 - t999 * t916;
t1006 = -t1000 * t866 - t999 * t865 - t1007 * t885 - t990 * t904 - t991 * t905 + mrSges(6,1) * t847 - mrSges(7,1) * t845 - mrSges(6,2) * t848 + mrSges(7,3) * t844 - pkin(5) * t841 + qJ(6) * (-t865 * mrSges(7,2) - t904 * t880 + t984);
t1002 = -mrSges(6,3) - mrSges(7,2);
t997 = t963 * t968;
t914 = -g(3) * t997 + t988;
t941 = t958 * mrSges(3,1) - mrSges(3,3) * t983;
t946 = (-mrSges(3,1) * t971 + mrSges(3,2) * t968) * t987;
t890 = t916 * mrSges(6,1) - t905 * mrSges(6,3);
t989 = -t904 * mrSges(6,1) - t905 * mrSges(6,2) - t880;
t836 = m(6) * t848 - t885 * mrSges(6,2) + t1002 * t865 - t916 * t890 + t989 * t904 + t984;
t889 = -t916 * mrSges(6,2) - t904 * mrSges(6,3);
t838 = m(6) * t847 + t885 * mrSges(6,1) + t1002 * t866 + t916 * t889 + t989 * t905 + t978;
t831 = t1005 * t836 - t966 * t838;
t896 = -t917 * mrSges(5,1) + t918 * mrSges(5,2);
t907 = t952 * mrSges(5,1) - t918 * mrSges(5,3);
t826 = m(5) * t853 - t940 * mrSges(5,2) + t886 * mrSges(5,3) + t917 * t896 - t952 * t907 + t831;
t852 = t970 * t857 - t967 * t860;
t850 = -t940 * pkin(4) - t951 * pkin(10) + t918 * t897 - t852;
t846 = -0.2e1 * qJD(6) * t905 + (t904 * t916 - t866) * qJ(6) + (t905 * t916 + t865) * pkin(5) + t850;
t842 = m(7) * t846 + t865 * mrSges(7,1) - t866 * mrSges(7,3) + t904 * t888 - t905 * t891;
t839 = -m(6) * t850 - t865 * mrSges(6,1) - t866 * mrSges(6,2) - t904 * t889 - t905 * t890 - t842;
t906 = -t952 * mrSges(5,2) + t917 * mrSges(5,3);
t833 = m(5) * t852 + t940 * mrSges(5,1) - t887 * mrSges(5,3) - t918 * t896 + t952 * t906 + t839;
t820 = t967 * t826 + t970 * t833;
t919 = -t936 * mrSges(4,1) + t937 * mrSges(4,2);
t921 = mrSges(4,2) * t982 + t936 * mrSges(4,3);
t818 = m(4) * t861 - t948 * mrSges(4,1) - t924 * mrSges(4,3) - t937 * t919 - t921 * t982 + t820;
t922 = -mrSges(4,1) * t982 - t937 * mrSges(4,3);
t979 = t970 * t826 - t967 * t833;
t819 = m(4) * t862 + t948 * mrSges(4,2) + t923 * mrSges(4,3) + t936 * t919 + t922 * t982 + t979;
t980 = -t962 * t818 + t964 * t819;
t809 = m(3) * t914 - t957 * mrSges(3,2) + t948 * mrSges(3,3) - t958 * t941 + t946 * t982 + t980;
t812 = t964 * t818 + t962 * t819;
t929 = -t963 * t943 - t1003;
t942 = -t958 * mrSges(3,2) + mrSges(3,3) * t982;
t811 = m(3) * t929 - t948 * mrSges(3,1) + t947 * mrSges(3,2) + (t941 * t968 - t942 * t971) * t987 + t812;
t830 = t1005 * t838 + t966 * t836;
t976 = m(5) * t867 - t886 * mrSges(5,1) + t887 * mrSges(5,2) - t917 * t906 + t918 * t907 + t830;
t827 = m(4) * t901 - t923 * mrSges(4,1) + t924 * mrSges(4,2) - t936 * t921 + t937 * t922 + t976;
t823 = m(3) * t913 + t957 * mrSges(3,1) - t947 * mrSges(3,3) + t958 * t942 - t946 * t983 - t827;
t799 = t809 * t995 - t963 * t811 + t823 * t994;
t796 = m(2) * t953 + qJDD(1) * mrSges(2,1) - t973 * mrSges(2,2) + t799;
t805 = t971 * t809 - t968 * t823;
t803 = m(2) * t954 - t973 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t805;
t993 = t972 * t796 + t969 * t803;
t992 = t1000 * t905 + t1007 * t916 + t999 * t904;
t798 = t809 * t997 + t965 * t811 + t823 * t996;
t981 = -t969 * t796 + t972 * t803;
t828 = -mrSges(6,1) * t850 - mrSges(7,1) * t846 + mrSges(7,2) * t844 + mrSges(6,3) * t848 - pkin(5) * t842 + t1001 * t866 - t1008 * t865 + t999 * t885 + t992 * t905 - t990 * t916;
t829 = mrSges(6,2) * t850 + mrSges(7,2) * t845 - mrSges(6,3) * t847 - mrSges(7,3) * t846 - qJ(6) * t842 - t1000 * t885 - t1001 * t865 + t1009 * t866 + t992 * t904 + t991 * t916;
t892 = Ifges(5,5) * t918 + Ifges(5,6) * t917 + Ifges(5,3) * t952;
t893 = Ifges(5,4) * t918 + Ifges(5,2) * t917 + Ifges(5,6) * t952;
t813 = mrSges(5,2) * t867 - mrSges(5,3) * t852 + Ifges(5,1) * t887 + Ifges(5,4) * t886 + Ifges(5,5) * t940 - pkin(10) * t830 + t1005 * t829 - t966 * t828 + t917 * t892 - t952 * t893;
t894 = Ifges(5,1) * t918 + Ifges(5,4) * t917 + Ifges(5,5) * t952;
t814 = -mrSges(5,1) * t867 + mrSges(5,3) * t853 + Ifges(5,4) * t887 + Ifges(5,2) * t886 + Ifges(5,6) * t940 - pkin(4) * t830 - t918 * t892 + t952 * t894 - t1006;
t908 = Ifges(4,5) * t937 + Ifges(4,6) * t936 - Ifges(4,3) * t982;
t910 = Ifges(4,1) * t937 + Ifges(4,4) * t936 - Ifges(4,5) * t982;
t794 = -mrSges(4,1) * t901 + mrSges(4,3) * t862 + Ifges(4,4) * t924 + Ifges(4,2) * t923 - Ifges(4,6) * t948 - pkin(3) * t976 + pkin(9) * t979 + t967 * t813 + t970 * t814 - t937 * t908 - t910 * t982;
t909 = Ifges(4,4) * t937 + Ifges(4,2) * t936 - Ifges(4,6) * t982;
t800 = mrSges(4,2) * t901 - mrSges(4,3) * t861 + Ifges(4,1) * t924 + Ifges(4,4) * t923 - Ifges(4,5) * t948 - pkin(9) * t820 + t970 * t813 - t967 * t814 + t936 * t908 + t909 * t982;
t927 = Ifges(3,6) * t958 + (Ifges(3,4) * t968 + Ifges(3,2) * t971) * t987;
t928 = Ifges(3,5) * t958 + (Ifges(3,1) * t968 + Ifges(3,4) * t971) * t987;
t789 = Ifges(3,5) * t947 + Ifges(3,6) * t948 + Ifges(3,3) * t957 + mrSges(3,1) * t913 - mrSges(3,2) * t914 + t962 * t800 + t964 * t794 - pkin(2) * t827 + qJ(3) * t980 + (t927 * t968 - t928 * t971) * t987;
t926 = Ifges(3,3) * t958 + (Ifges(3,5) * t968 + Ifges(3,6) * t971) * t987;
t791 = mrSges(3,2) * t929 - mrSges(3,3) * t913 + Ifges(3,1) * t947 + Ifges(3,4) * t948 + Ifges(3,5) * t957 - qJ(3) * t812 - t962 * t794 + t964 * t800 + t926 * t982 - t958 * t927;
t974 = mrSges(5,1) * t852 - mrSges(5,2) * t853 + Ifges(5,5) * t887 + Ifges(5,6) * t886 + Ifges(5,3) * t940 + pkin(4) * t839 + pkin(10) * t831 + t1005 * t828 + t966 * t829 + t918 * t893 - t917 * t894;
t793 = -t926 * t983 + Ifges(3,6) * t957 + t958 * t928 + Ifges(3,4) * t947 + t936 * t910 - t937 * t909 - Ifges(4,5) * t924 - mrSges(3,1) * t929 - Ifges(4,6) * t923 + mrSges(3,3) * t914 - t974 - mrSges(4,1) * t861 + mrSges(4,2) * t862 - pkin(3) * t820 + (Ifges(3,2) + Ifges(4,3)) * t948 - pkin(2) * t812;
t977 = mrSges(2,1) * t953 - mrSges(2,2) * t954 + Ifges(2,3) * qJDD(1) + pkin(1) * t799 + t805 * t1004 + t965 * t789 + t791 * t997 + t793 * t996;
t787 = -mrSges(2,2) * g(3) - mrSges(2,3) * t953 + Ifges(2,5) * qJDD(1) - t973 * Ifges(2,6) + t971 * t791 - t968 * t793 + (-t798 * t963 - t799 * t965) * pkin(8);
t786 = mrSges(2,1) * g(3) + mrSges(2,3) * t954 + t973 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t798 - t963 * t789 + (pkin(8) * t805 + t791 * t968 + t793 * t971) * t965;
t1 = [-m(1) * g(1) + t981; -m(1) * g(2) + t993; (-m(1) - m(2)) * g(3) + t798; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t993 - t969 * t786 + t972 * t787; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t981 + t972 * t786 + t969 * t787; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t977; t977; t789; t827; t974; t1006; t841;];
tauJB  = t1;
