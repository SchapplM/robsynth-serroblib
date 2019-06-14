% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 06:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:16:29
% EndTime: 2019-05-05 06:16:45
% DurationCPUTime: 15.71s
% Computational Cost: add. (266701->339), mult. (572665->432), div. (0->0), fcn. (423392->14), ass. (0->159)
t1015 = -2 * qJD(4);
t1014 = Ifges(4,1) + Ifges(5,2);
t1007 = Ifges(4,4) + Ifges(5,6);
t1006 = Ifges(4,5) - Ifges(5,4);
t1013 = Ifges(4,2) + Ifges(5,3);
t1005 = Ifges(4,6) - Ifges(5,5);
t1012 = Ifges(4,3) + Ifges(5,1);
t954 = sin(pkin(12));
t957 = cos(pkin(12));
t939 = g(1) * t954 - g(2) * t957;
t940 = -g(1) * t957 - g(2) * t954;
t953 = -g(3) + qJDD(1);
t963 = sin(qJ(2));
t959 = cos(pkin(6));
t967 = cos(qJ(2));
t996 = t959 * t967;
t956 = sin(pkin(6));
t999 = t956 * t967;
t894 = t939 * t996 - t940 * t963 + t953 * t999;
t955 = sin(pkin(7));
t968 = qJD(2) ^ 2;
t890 = pkin(9) * t955 * t968 + qJDD(2) * pkin(2) + t894;
t958 = cos(pkin(7));
t1003 = t890 * t958;
t917 = -t939 * t956 + t953 * t959;
t966 = cos(qJ(3));
t1011 = (t917 * t955 + t1003) * t966;
t949 = qJD(2) * t958 + qJD(3);
t962 = sin(qJ(3));
t991 = qJD(2) * t955;
t986 = t962 * t991;
t1010 = (pkin(3) * t949 + t1015) * t986;
t1001 = t955 * t962;
t1000 = t956 * t963;
t997 = t959 * t963;
t895 = t953 * t1000 + t939 * t997 + t967 * t940;
t988 = qJDD(2) * t955;
t891 = -pkin(2) * t968 + pkin(9) * t988 + t895;
t998 = t958 * t962;
t864 = t917 * t1001 + t890 * t998 + t966 * t891;
t926 = (-pkin(3) * t966 - qJ(4) * t962) * t991;
t947 = t949 ^ 2;
t948 = qJDD(2) * t958 + qJDD(3);
t990 = qJD(2) * t966;
t985 = t955 * t990;
t859 = t947 * pkin(3) - t948 * qJ(4) + t949 * t1015 - t926 * t985 - t864;
t1009 = -pkin(3) - pkin(10);
t1008 = mrSges(4,1) - mrSges(5,2);
t888 = t962 * t891;
t863 = -t888 + t1011;
t923 = -mrSges(4,2) * t949 + mrSges(4,3) * t985;
t924 = -mrSges(5,1) * t985 - mrSges(5,3) * t949;
t927 = (mrSges(5,2) * t966 - mrSges(5,3) * t962) * t991;
t928 = (-mrSges(4,1) * t966 + mrSges(4,2) * t962) * t991;
t930 = (qJD(3) * t990 + qJDD(2) * t962) * t955;
t1002 = t955 ^ 2 * t968;
t979 = -t947 * qJ(4) + t926 * t986 + qJDD(4) + t888;
t856 = t930 * pkin(4) + t1009 * t948 + (-pkin(10) * t962 * t1002 - t1003 + (-pkin(4) * qJD(2) * t949 - t917) * t955) * t966 + t979;
t913 = t958 * t917;
t929 = pkin(4) * t986 - pkin(10) * t949;
t931 = -qJD(3) * t986 + t966 * t988;
t987 = t966 ^ 2 * t1002;
t858 = -pkin(4) * t987 - t930 * qJ(4) + t913 + t1009 * t931 + (-t890 + (-qJ(4) * t949 * t966 - t929 * t962) * qJD(2)) * t955 + t1010;
t961 = sin(qJ(5));
t965 = cos(qJ(5));
t851 = t961 * t856 + t965 * t858;
t916 = t949 * t965 - t961 * t985;
t885 = -qJD(5) * t916 - t931 * t965 - t948 * t961;
t915 = -t949 * t961 - t965 * t985;
t892 = -mrSges(6,1) * t915 + mrSges(6,2) * t916;
t938 = qJD(5) + t986;
t899 = mrSges(6,1) * t938 - mrSges(6,3) * t916;
t921 = qJDD(5) + t930;
t893 = -pkin(5) * t915 - pkin(11) * t916;
t936 = t938 ^ 2;
t848 = -pkin(5) * t936 + pkin(11) * t921 + t893 * t915 + t851;
t854 = t931 * pkin(4) - pkin(10) * t987 + t949 * t929 - t859;
t886 = qJD(5) * t915 - t931 * t961 + t948 * t965;
t852 = (-t915 * t938 - t886) * pkin(11) + (t916 * t938 - t885) * pkin(5) + t854;
t960 = sin(qJ(6));
t964 = cos(qJ(6));
t845 = -t848 * t960 + t852 * t964;
t896 = -t916 * t960 + t938 * t964;
t867 = qJD(6) * t896 + t886 * t964 + t921 * t960;
t897 = t916 * t964 + t938 * t960;
t872 = -mrSges(7,1) * t896 + mrSges(7,2) * t897;
t914 = qJD(6) - t915;
t875 = -mrSges(7,2) * t914 + mrSges(7,3) * t896;
t883 = qJDD(6) - t885;
t842 = m(7) * t845 + mrSges(7,1) * t883 - mrSges(7,3) * t867 - t872 * t897 + t875 * t914;
t846 = t848 * t964 + t852 * t960;
t866 = -qJD(6) * t897 - t886 * t960 + t921 * t964;
t876 = mrSges(7,1) * t914 - mrSges(7,3) * t897;
t843 = m(7) * t846 - mrSges(7,2) * t883 + mrSges(7,3) * t866 + t872 * t896 - t876 * t914;
t982 = -t842 * t960 + t964 * t843;
t831 = m(6) * t851 - mrSges(6,2) * t921 + mrSges(6,3) * t885 + t892 * t915 - t899 * t938 + t982;
t850 = t856 * t965 - t858 * t961;
t898 = -mrSges(6,2) * t938 + mrSges(6,3) * t915;
t847 = -pkin(5) * t921 - pkin(11) * t936 + t893 * t916 - t850;
t973 = -m(7) * t847 + t866 * mrSges(7,1) - mrSges(7,2) * t867 + t896 * t875 - t876 * t897;
t838 = m(6) * t850 + mrSges(6,1) * t921 - mrSges(6,3) * t886 - t892 * t916 + t898 * t938 + t973;
t825 = t961 * t831 + t965 * t838;
t860 = -t948 * pkin(3) - t1011 + t979;
t974 = -m(5) * t860 - t930 * mrSges(5,1) - t825;
t820 = m(4) * t863 - t930 * mrSges(4,3) + (t923 - t924) * t949 + t1008 * t948 + (-t927 - t928) * t986 + t974;
t1004 = t820 * t966;
t874 = -t955 * t890 + t913;
t922 = mrSges(4,1) * t949 - mrSges(4,3) * t986;
t925 = mrSges(5,1) * t986 + mrSges(5,2) * t949;
t862 = -t931 * pkin(3) + (-t949 * t985 - t930) * qJ(4) + t874 + t1010;
t983 = t965 * t831 - t961 * t838;
t978 = m(5) * t862 - t930 * mrSges(5,3) + t924 * t985 + t983;
t822 = m(4) * t874 + t930 * mrSges(4,2) - t1008 * t931 + (-t923 * t966 + (t922 - t925) * t962) * t991 + t978;
t833 = t964 * t842 + t960 * t843;
t972 = -m(6) * t854 + t885 * mrSges(6,1) - t886 * mrSges(6,2) + t915 * t898 - t916 * t899 - t833;
t970 = -m(5) * t859 + t948 * mrSges(5,3) + t949 * t925 + t927 * t985 - t972;
t829 = t928 * t985 - t948 * mrSges(4,2) - t949 * t922 + m(4) * t864 + t970 + (mrSges(4,3) + mrSges(5,1)) * t931;
t810 = t958 * t1004 - t822 * t955 + t829 * t998;
t806 = m(3) * t894 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t968 + t810;
t809 = t829 * t1001 + t955 * t1004 + t958 * t822;
t808 = m(3) * t917 + t809;
t815 = -t820 * t962 + t966 * t829;
t814 = m(3) * t895 - mrSges(3,1) * t968 - qJDD(2) * mrSges(3,2) + t815;
t796 = t806 * t996 - t808 * t956 + t814 * t997;
t794 = m(2) * t939 + t796;
t802 = -t806 * t963 + t967 * t814;
t801 = m(2) * t940 + t802;
t995 = t957 * t794 + t954 * t801;
t994 = (t1005 * t966 + t1006 * t962) * t991 + t1012 * t949;
t993 = (t1007 * t966 + t1014 * t962) * t991 + t1006 * t949;
t992 = (-t1007 * t962 - t1013 * t966) * t991 - t1005 * t949;
t795 = t814 * t1000 + t806 * t999 + t959 * t808;
t984 = -t794 * t954 + t957 * t801;
t981 = m(2) * t953 + t795;
t868 = Ifges(7,5) * t897 + Ifges(7,6) * t896 + Ifges(7,3) * t914;
t870 = Ifges(7,1) * t897 + Ifges(7,4) * t896 + Ifges(7,5) * t914;
t836 = -mrSges(7,1) * t847 + mrSges(7,3) * t846 + Ifges(7,4) * t867 + Ifges(7,2) * t866 + Ifges(7,6) * t883 - t868 * t897 + t870 * t914;
t869 = Ifges(7,4) * t897 + Ifges(7,2) * t896 + Ifges(7,6) * t914;
t837 = mrSges(7,2) * t847 - mrSges(7,3) * t845 + Ifges(7,1) * t867 + Ifges(7,4) * t866 + Ifges(7,5) * t883 + t868 * t896 - t869 * t914;
t877 = Ifges(6,5) * t916 + Ifges(6,6) * t915 + Ifges(6,3) * t938;
t878 = Ifges(6,4) * t916 + Ifges(6,2) * t915 + Ifges(6,6) * t938;
t816 = mrSges(6,2) * t854 - mrSges(6,3) * t850 + Ifges(6,1) * t886 + Ifges(6,4) * t885 + Ifges(6,5) * t921 - pkin(11) * t833 - t836 * t960 + t837 * t964 + t877 * t915 - t878 * t938;
t879 = Ifges(6,1) * t916 + Ifges(6,4) * t915 + Ifges(6,5) * t938;
t969 = mrSges(7,1) * t845 - mrSges(7,2) * t846 + Ifges(7,5) * t867 + Ifges(7,6) * t866 + Ifges(7,3) * t883 + t869 * t897 - t870 * t896;
t817 = -mrSges(6,1) * t854 + mrSges(6,3) * t851 + Ifges(6,4) * t886 + Ifges(6,2) * t885 + Ifges(6,6) * t921 - pkin(5) * t833 - t877 * t916 + t879 * t938 - t969;
t823 = t948 * mrSges(5,2) + t949 * t924 + t927 * t986 - t974;
t797 = mrSges(4,1) * t863 - mrSges(4,2) * t864 + mrSges(5,2) * t860 - mrSges(5,3) * t859 + t965 * t816 - t961 * t817 - pkin(10) * t825 - pkin(3) * t823 + qJ(4) * t970 + t1012 * t948 + (mrSges(5,1) * qJ(4) + t1005) * t931 + t1006 * t930 + (-t992 * t962 - t993 * t966) * t991;
t824 = t931 * mrSges(5,2) - t925 * t986 + t978;
t798 = -mrSges(4,1) * t874 - mrSges(5,1) * t859 + mrSges(5,2) * t862 + mrSges(4,3) * t864 - pkin(3) * t824 - pkin(4) * t972 - pkin(10) * t983 + t1005 * t948 + t1007 * t930 + t1013 * t931 - t961 * t816 - t965 * t817 + t993 * t949 - t994 * t986;
t971 = mrSges(6,1) * t850 - mrSges(6,2) * t851 + Ifges(6,5) * t886 + Ifges(6,6) * t885 + Ifges(6,3) * t921 + pkin(5) * t973 + pkin(11) * t982 + t964 * t836 + t960 * t837 + t916 * t878 - t915 * t879;
t803 = mrSges(5,1) * t860 + mrSges(4,2) * t874 - mrSges(4,3) * t863 - mrSges(5,3) * t862 + pkin(4) * t825 - qJ(4) * t824 + t1006 * t948 + t1007 * t931 + t1014 * t930 + t992 * t949 + t994 * t985 + t971;
t975 = pkin(9) * t815 + t798 * t966 + t803 * t962;
t791 = -mrSges(3,1) * t917 + mrSges(3,3) * t895 + t968 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t809 - t955 * t797 + t975 * t958;
t792 = mrSges(3,2) * t917 - mrSges(3,3) * t894 + Ifges(3,5) * qJDD(2) - t968 * Ifges(3,6) - t962 * t798 + t966 * t803 + (-t809 * t955 - t810 * t958) * pkin(9);
t976 = pkin(8) * t802 + t791 * t967 + t792 * t963;
t790 = mrSges(3,1) * t894 - mrSges(3,2) * t895 + Ifges(3,3) * qJDD(2) + pkin(2) * t810 + t958 * t797 + t975 * t955;
t789 = mrSges(2,2) * t953 - mrSges(2,3) * t939 - t963 * t791 + t967 * t792 + (-t795 * t956 - t796 * t959) * pkin(8);
t788 = -mrSges(2,1) * t953 + mrSges(2,3) * t940 - pkin(1) * t795 - t956 * t790 + t976 * t959;
t1 = [-m(1) * g(1) + t984; -m(1) * g(2) + t995; -m(1) * g(3) + t981; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t995 - t954 * t788 + t957 * t789; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t984 + t957 * t788 + t954 * t789; -mrSges(1,1) * g(2) + mrSges(2,1) * t939 + mrSges(1,2) * g(1) - mrSges(2,2) * t940 + pkin(1) * t796 + t959 * t790 + t976 * t956; t981; t790; t797; t823; t971; t969;];
tauJB  = t1;
