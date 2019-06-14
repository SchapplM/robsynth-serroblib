% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:30:39
% EndTime: 2019-05-06 17:30:55
% DurationCPUTime: 15.62s
% Computational Cost: add. (231850->365), mult. (528785->447), div. (0->0), fcn. (380008->10), ass. (0->143)
t985 = Ifges(6,4) + Ifges(7,4);
t995 = Ifges(6,2) + Ifges(7,2);
t991 = Ifges(6,6) + Ifges(7,6);
t952 = sin(qJ(2));
t956 = cos(qJ(2));
t976 = qJD(1) * qJD(2);
t936 = t952 * qJDD(1) + t956 * t976;
t937 = t956 * qJDD(1) - t952 * t976;
t948 = sin(pkin(10));
t949 = cos(pkin(10));
t911 = t949 * t936 + t948 * t937;
t925 = (t948 * t956 + t949 * t952) * qJD(1);
t951 = sin(qJ(4));
t955 = cos(qJ(4));
t914 = t951 * qJD(2) + t955 * t925;
t882 = -t914 * qJD(4) + t955 * qJDD(2) - t951 * t911;
t913 = t955 * qJD(2) - t951 * t925;
t883 = t913 * qJD(4) + t951 * qJDD(2) + t955 * t911;
t950 = sin(qJ(5));
t954 = cos(qJ(5));
t888 = t954 * t913 - t950 * t914;
t847 = t888 * qJD(5) + t950 * t882 + t954 * t883;
t889 = t950 * t913 + t954 * t914;
t868 = -t888 * mrSges(7,1) + t889 * mrSges(7,2);
t953 = sin(qJ(1));
t957 = cos(qJ(1));
t943 = -t957 * g(1) - t953 * g(2);
t959 = qJD(1) ^ 2;
t931 = -t959 * pkin(1) + qJDD(1) * pkin(7) + t943;
t983 = t952 * t931;
t986 = pkin(2) * t959;
t890 = qJDD(2) * pkin(2) - t936 * qJ(3) - t983 + (qJ(3) * t976 + t952 * t986 - g(3)) * t956;
t916 = -t952 * g(3) + t956 * t931;
t979 = qJD(1) * t952;
t939 = qJD(2) * pkin(2) - qJ(3) * t979;
t947 = t956 ^ 2;
t891 = t937 * qJ(3) - qJD(2) * t939 - t947 * t986 + t916;
t978 = qJD(1) * t956;
t924 = -t948 * t979 + t949 * t978;
t867 = 0.2e1 * qJD(3) * t924 + t948 * t890 + t949 * t891;
t906 = -t924 * pkin(3) - t925 * pkin(8);
t958 = qJD(2) ^ 2;
t852 = -t958 * pkin(3) + qJDD(2) * pkin(8) + t924 * t906 + t867;
t942 = t953 * g(1) - t957 * g(2);
t966 = -qJDD(1) * pkin(1) - t942;
t894 = -t937 * pkin(2) + qJDD(3) + t939 * t979 + (-qJ(3) * t947 - pkin(7)) * t959 + t966;
t910 = -t948 * t936 + t949 * t937;
t856 = (-qJD(2) * t924 - t911) * pkin(8) + (qJD(2) * t925 - t910) * pkin(3) + t894;
t836 = -t951 * t852 + t955 * t856;
t909 = qJDD(4) - t910;
t923 = qJD(4) - t924;
t832 = (t913 * t923 - t883) * pkin(9) + (t913 * t914 + t909) * pkin(4) + t836;
t837 = t955 * t852 + t951 * t856;
t897 = t923 * pkin(4) - t914 * pkin(9);
t912 = t913 ^ 2;
t834 = -t912 * pkin(4) + t882 * pkin(9) - t923 * t897 + t837;
t826 = t954 * t832 - t950 * t834;
t907 = qJDD(5) + t909;
t919 = qJD(5) + t923;
t821 = -0.2e1 * qJD(6) * t889 + (t888 * t919 - t847) * qJ(6) + (t888 * t889 + t907) * pkin(5) + t826;
t871 = -t919 * mrSges(7,2) + t888 * mrSges(7,3);
t975 = m(7) * t821 + t907 * mrSges(7,1) + t919 * t871;
t818 = -t847 * mrSges(7,3) - t889 * t868 + t975;
t827 = t950 * t832 + t954 * t834;
t846 = -t889 * qJD(5) + t954 * t882 - t950 * t883;
t873 = t919 * pkin(5) - t889 * qJ(6);
t887 = t888 ^ 2;
t823 = -t887 * pkin(5) + t846 * qJ(6) + 0.2e1 * qJD(6) * t888 - t919 * t873 + t827;
t992 = Ifges(6,5) + Ifges(7,5);
t993 = Ifges(6,1) + Ifges(7,1);
t980 = -t985 * t888 - t889 * t993 - t992 * t919;
t989 = t995 * t888 + t889 * t985 + t991 * t919;
t990 = Ifges(6,3) + Ifges(7,3);
t994 = mrSges(6,1) * t826 + mrSges(7,1) * t821 - mrSges(6,2) * t827 - mrSges(7,2) * t823 + pkin(5) * t818 + t846 * t991 + t847 * t992 + t980 * t888 + t889 * t989 + t907 * t990;
t869 = -t888 * mrSges(6,1) + t889 * mrSges(6,2);
t872 = -t919 * mrSges(6,2) + t888 * mrSges(6,3);
t810 = m(6) * t826 + t907 * mrSges(6,1) + t919 * t872 + (-t868 - t869) * t889 + (-mrSges(6,3) - mrSges(7,3)) * t847 + t975;
t874 = t919 * mrSges(7,1) - t889 * mrSges(7,3);
t875 = t919 * mrSges(6,1) - t889 * mrSges(6,3);
t974 = m(7) * t823 + t846 * mrSges(7,3) + t888 * t868;
t813 = m(6) * t827 + t846 * mrSges(6,3) + t888 * t869 + (-t874 - t875) * t919 + (-mrSges(6,2) - mrSges(7,2)) * t907 + t974;
t808 = t954 * t810 + t950 * t813;
t877 = Ifges(5,4) * t914 + Ifges(5,2) * t913 + Ifges(5,6) * t923;
t878 = Ifges(5,1) * t914 + Ifges(5,4) * t913 + Ifges(5,5) * t923;
t988 = mrSges(5,1) * t836 - mrSges(5,2) * t837 + Ifges(5,5) * t883 + Ifges(5,6) * t882 + Ifges(5,3) * t909 + pkin(4) * t808 + t914 * t877 - t913 * t878 + t994;
t866 = -0.2e1 * qJD(3) * t925 + t949 * t890 - t948 * t891;
t851 = -qJDD(2) * pkin(3) - t958 * pkin(8) + t925 * t906 - t866;
t835 = -t882 * pkin(4) - t912 * pkin(9) + t914 * t897 + t851;
t829 = -t846 * pkin(5) - t887 * qJ(6) + t889 * t873 + qJDD(6) + t835;
t824 = m(7) * t829 - t846 * mrSges(7,1) + t847 * mrSges(7,2) - t888 * t871 + t889 * t874;
t981 = -t888 * t991 - t889 * t992 - t919 * t990;
t801 = -mrSges(6,1) * t835 + mrSges(6,3) * t827 - mrSges(7,1) * t829 + mrSges(7,3) * t823 - pkin(5) * t824 + qJ(6) * t974 + (-qJ(6) * t874 - t980) * t919 + (-qJ(6) * mrSges(7,2) + t991) * t907 + t981 * t889 + t985 * t847 + t995 * t846;
t807 = mrSges(6,2) * t835 + mrSges(7,2) * t829 - mrSges(6,3) * t826 - mrSges(7,3) * t821 - qJ(6) * t818 + t985 * t846 + t847 * t993 - t981 * t888 + t992 * t907 - t989 * t919;
t876 = Ifges(5,5) * t914 + Ifges(5,6) * t913 + Ifges(5,3) * t923;
t964 = m(6) * t835 - t846 * mrSges(6,1) + t847 * mrSges(6,2) - t888 * t872 + t889 * t875 + t824;
t969 = -t950 * t810 + t954 * t813;
t784 = -mrSges(5,1) * t851 + mrSges(5,3) * t837 + Ifges(5,4) * t883 + Ifges(5,2) * t882 + Ifges(5,6) * t909 - pkin(4) * t964 + pkin(9) * t969 + t954 * t801 + t950 * t807 - t914 * t876 + t923 * t878;
t785 = mrSges(5,2) * t851 - mrSges(5,3) * t836 + Ifges(5,1) * t883 + Ifges(5,4) * t882 + Ifges(5,5) * t909 - pkin(9) * t808 - t950 * t801 + t954 * t807 + t913 * t876 - t923 * t877;
t892 = -t913 * mrSges(5,1) + t914 * mrSges(5,2);
t895 = -t923 * mrSges(5,2) + t913 * mrSges(5,3);
t805 = m(5) * t836 + t909 * mrSges(5,1) - t883 * mrSges(5,3) - t914 * t892 + t923 * t895 + t808;
t896 = t923 * mrSges(5,1) - t914 * mrSges(5,3);
t806 = m(5) * t837 - t909 * mrSges(5,2) + t882 * mrSges(5,3) + t913 * t892 - t923 * t896 + t969;
t800 = -t951 * t805 + t955 * t806;
t902 = -t924 * mrSges(4,1) + t925 * mrSges(4,2);
t918 = qJD(2) * mrSges(4,1) - t925 * mrSges(4,3);
t797 = m(4) * t867 - qJDD(2) * mrSges(4,2) + t910 * mrSges(4,3) - qJD(2) * t918 + t924 * t902 + t800;
t816 = -m(5) * t851 + t882 * mrSges(5,1) - t883 * mrSges(5,2) + t913 * t895 - t914 * t896 - t964;
t917 = -qJD(2) * mrSges(4,2) + t924 * mrSges(4,3);
t815 = m(4) * t866 + qJDD(2) * mrSges(4,1) - t911 * mrSges(4,3) + qJD(2) * t917 - t925 * t902 + t816;
t791 = t948 * t797 + t949 * t815;
t899 = Ifges(4,4) * t925 + Ifges(4,2) * t924 + Ifges(4,6) * qJD(2);
t900 = Ifges(4,1) * t925 + Ifges(4,4) * t924 + Ifges(4,5) * qJD(2);
t915 = -t956 * g(3) - t983;
t927 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t952 + Ifges(3,2) * t956) * qJD(1);
t928 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t952 + Ifges(3,4) * t956) * qJD(1);
t987 = (t952 * t927 - t956 * t928) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + mrSges(3,1) * t915 + mrSges(4,1) * t866 - mrSges(3,2) * t916 - mrSges(4,2) * t867 + Ifges(3,5) * t936 + Ifges(4,5) * t911 + Ifges(3,6) * t937 + Ifges(4,6) * t910 + pkin(2) * t791 + pkin(3) * t816 + pkin(8) * t800 + t955 * t784 + t951 * t785 + t925 * t899 - t924 * t900;
t935 = (-mrSges(3,1) * t956 + mrSges(3,2) * t952) * qJD(1);
t941 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t978;
t789 = m(3) * t915 + qJDD(2) * mrSges(3,1) - t936 * mrSges(3,3) + qJD(2) * t941 - t935 * t979 + t791;
t940 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t979;
t970 = t949 * t797 - t948 * t815;
t790 = m(3) * t916 - qJDD(2) * mrSges(3,2) + t937 * mrSges(3,3) - qJD(2) * t940 + t935 * t978 + t970;
t971 = -t952 * t789 + t956 * t790;
t780 = m(2) * t943 - t959 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t971;
t799 = t955 * t805 + t951 * t806;
t798 = m(4) * t894 - t910 * mrSges(4,1) + t911 * mrSges(4,2) - t924 * t917 + t925 * t918 + t799;
t930 = -t959 * pkin(7) + t966;
t962 = -m(3) * t930 + t937 * mrSges(3,1) - t936 * mrSges(3,2) - t940 * t979 + t941 * t978 - t798;
t793 = m(2) * t942 + qJDD(1) * mrSges(2,1) - t959 * mrSges(2,2) + t962;
t982 = t953 * t780 + t957 * t793;
t782 = t956 * t789 + t952 * t790;
t972 = t957 * t780 - t953 * t793;
t898 = Ifges(4,5) * t925 + Ifges(4,6) * t924 + Ifges(4,3) * qJD(2);
t777 = mrSges(4,2) * t894 - mrSges(4,3) * t866 + Ifges(4,1) * t911 + Ifges(4,4) * t910 + Ifges(4,5) * qJDD(2) - pkin(8) * t799 - qJD(2) * t899 - t951 * t784 + t955 * t785 + t924 * t898;
t783 = -mrSges(4,1) * t894 + mrSges(4,3) * t867 + Ifges(4,4) * t911 + Ifges(4,2) * t910 + Ifges(4,6) * qJDD(2) - pkin(3) * t799 + qJD(2) * t900 - t925 * t898 - t988;
t926 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t952 + Ifges(3,6) * t956) * qJD(1);
t773 = -mrSges(3,1) * t930 + mrSges(3,3) * t916 + Ifges(3,4) * t936 + Ifges(3,2) * t937 + Ifges(3,6) * qJDD(2) - pkin(2) * t798 + qJ(3) * t970 + qJD(2) * t928 + t948 * t777 + t949 * t783 - t926 * t979;
t776 = mrSges(3,2) * t930 - mrSges(3,3) * t915 + Ifges(3,1) * t936 + Ifges(3,4) * t937 + Ifges(3,5) * qJDD(2) - qJ(3) * t791 - qJD(2) * t927 + t949 * t777 - t948 * t783 + t926 * t978;
t965 = mrSges(2,1) * t942 - mrSges(2,2) * t943 + Ifges(2,3) * qJDD(1) + pkin(1) * t962 + pkin(7) * t971 + t956 * t773 + t952 * t776;
t774 = mrSges(2,1) * g(3) + mrSges(2,3) * t943 + t959 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t782 - t987;
t771 = -mrSges(2,2) * g(3) - mrSges(2,3) * t942 + Ifges(2,5) * qJDD(1) - t959 * Ifges(2,6) - pkin(7) * t782 - t952 * t773 + t956 * t776;
t1 = [-m(1) * g(1) + t972; -m(1) * g(2) + t982; (-m(1) - m(2)) * g(3) + t782; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t982 + t957 * t771 - t953 * t774; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t972 + t953 * t771 + t957 * t774; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t965; t965; t987; t798; t988; t994; t824;];
tauJB  = t1;
