% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP6
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
% Datum: 2019-05-06 18:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:00:12
% EndTime: 2019-05-06 18:00:43
% DurationCPUTime: 27.02s
% Computational Cost: add. (406444->375), mult. (1061497->476), div. (0->0), fcn. (833626->12), ass. (0->156)
t1021 = -2 * qJD(3);
t1020 = Ifges(6,1) + Ifges(7,1);
t1013 = Ifges(6,4) - Ifges(7,5);
t1012 = -Ifges(6,5) - Ifges(7,4);
t1019 = Ifges(6,2) + Ifges(7,3);
t1011 = Ifges(6,6) - Ifges(7,6);
t1018 = -Ifges(6,3) - Ifges(7,2);
t970 = sin(pkin(11));
t972 = cos(pkin(11));
t976 = sin(qJ(2));
t979 = cos(qJ(2));
t971 = sin(pkin(6));
t999 = qJD(1) * t971;
t949 = (t970 * t976 - t972 * t979) * t999;
t981 = qJD(1) ^ 2;
t1009 = t971 ^ 2 * t981;
t997 = qJD(1) * qJD(2);
t958 = (qJDD(1) * t976 + t979 * t997) * t971;
t973 = cos(pkin(6));
t964 = qJDD(1) * t973 + qJDD(2);
t965 = qJD(1) * t973 + qJD(2);
t1005 = t973 * t979;
t1015 = pkin(8) * t971;
t977 = sin(qJ(1));
t980 = cos(qJ(1));
t961 = t977 * g(1) - g(2) * t980;
t955 = qJDD(1) * pkin(1) + t1015 * t981 + t961;
t962 = -g(1) * t980 - g(2) * t977;
t956 = -pkin(1) * t981 + qJDD(1) * t1015 + t962;
t988 = t955 * t1005 - t976 * t956;
t893 = t964 * pkin(2) - t958 * qJ(3) + (pkin(2) * t976 * t1009 + (qJ(3) * qJD(1) * t965 - g(3)) * t971) * t979 + t988;
t1006 = t973 * t976;
t1008 = t971 * t976;
t923 = -g(3) * t1008 + t955 * t1006 + t979 * t956;
t994 = t976 * t999;
t952 = pkin(2) * t965 - qJ(3) * t994;
t959 = (qJDD(1) * t979 - t976 * t997) * t971;
t996 = t979 ^ 2 * t1009;
t898 = -pkin(2) * t996 + qJ(3) * t959 - t952 * t965 + t923;
t950 = (t970 * t979 + t972 * t976) * t999;
t868 = t1021 * t950 + t972 * t893 - t970 * t898;
t1016 = cos(qJ(5));
t975 = sin(qJ(4));
t978 = cos(qJ(4));
t934 = t950 * t978 + t965 * t975;
t948 = qJD(4) + t949;
t974 = sin(qJ(5));
t914 = -t1016 * t948 + t974 * t934;
t915 = t1016 * t934 + t974 * t948;
t933 = -t950 * t975 + t965 * t978;
t932 = qJD(5) - t933;
t1001 = t1012 * t932 + t1013 * t914 - t1020 * t915;
t1002 = -t1011 * t932 - t1013 * t915 + t1019 * t914;
t930 = t958 * t972 + t959 * t970;
t908 = qJD(4) * t933 + t930 * t978 + t964 * t975;
t929 = -t958 * t970 + t959 * t972;
t928 = qJDD(4) - t929;
t875 = -t914 * qJD(5) + t1016 * t908 + t974 * t928;
t887 = mrSges(7,1) * t914 - mrSges(7,3) * t915;
t869 = t1021 * t949 + t970 * t893 + t972 * t898;
t926 = pkin(3) * t949 - pkin(9) * t950;
t963 = t965 ^ 2;
t867 = -pkin(3) * t963 + pkin(9) * t964 - t926 * t949 + t869;
t940 = -t973 * g(3) - t971 * t955;
t910 = -t959 * pkin(2) - qJ(3) * t996 + t952 * t994 + qJDD(3) + t940;
t871 = (t949 * t965 - t930) * pkin(9) + (t950 * t965 - t929) * pkin(3) + t910;
t863 = t978 * t867 + t975 * t871;
t912 = -pkin(4) * t933 - pkin(10) * t934;
t947 = t948 ^ 2;
t859 = -pkin(4) * t947 + pkin(10) * t928 + t912 * t933 + t863;
t866 = -t964 * pkin(3) - t963 * pkin(9) + t950 * t926 - t868;
t907 = -qJD(4) * t934 - t930 * t975 + t964 * t978;
t861 = (-t933 * t948 - t908) * pkin(10) + (t934 * t948 - t907) * pkin(4) + t866;
t855 = t1016 * t861 - t974 * t859;
t886 = pkin(5) * t914 - qJ(6) * t915;
t906 = qJDD(5) - t907;
t931 = t932 ^ 2;
t853 = -t906 * pkin(5) - t931 * qJ(6) + t915 * t886 + qJDD(6) - t855;
t894 = -mrSges(7,2) * t914 + mrSges(7,3) * t932;
t987 = -m(7) * t853 + t906 * mrSges(7,1) + t932 * t894;
t849 = t875 * mrSges(7,2) + t915 * t887 - t987;
t856 = t1016 * t859 + t974 * t861;
t852 = -pkin(5) * t931 + qJ(6) * t906 + 0.2e1 * qJD(6) * t932 - t886 * t914 + t856;
t874 = t915 * qJD(5) - t1016 * t928 + t974 * t908;
t897 = -mrSges(7,1) * t932 + mrSges(7,2) * t915;
t995 = m(7) * t852 + t906 * mrSges(7,3) + t932 * t897;
t1017 = -t1001 * t914 - t1002 * t915 - t1018 * t906 - t1011 * t874 - t1012 * t875 + mrSges(6,1) * t855 - mrSges(7,1) * t853 - mrSges(6,2) * t856 + mrSges(7,3) * t852 - pkin(5) * t849 + qJ(6) * (-t874 * mrSges(7,2) - t914 * t887 + t995);
t1014 = -mrSges(6,3) - mrSges(7,2);
t1007 = t971 * t979;
t924 = mrSges(4,1) * t949 + mrSges(4,2) * t950;
t936 = mrSges(4,1) * t965 - mrSges(4,3) * t950;
t1000 = -mrSges(6,1) * t914 - mrSges(6,2) * t915 - t887;
t896 = mrSges(6,1) * t932 - mrSges(6,3) * t915;
t845 = m(6) * t856 - t906 * mrSges(6,2) + t1000 * t914 + t1014 * t874 - t932 * t896 + t995;
t895 = -mrSges(6,2) * t932 - mrSges(6,3) * t914;
t846 = m(6) * t855 + t906 * mrSges(6,1) + t1000 * t915 + t1014 * t875 + t932 * t895 + t987;
t841 = t1016 * t845 - t846 * t974;
t911 = -mrSges(5,1) * t933 + mrSges(5,2) * t934;
t917 = mrSges(5,1) * t948 - mrSges(5,3) * t934;
t837 = m(5) * t863 - mrSges(5,2) * t928 + mrSges(5,3) * t907 + t911 * t933 - t917 * t948 + t841;
t862 = -t975 * t867 + t978 * t871;
t858 = -t928 * pkin(4) - t947 * pkin(10) + t934 * t912 - t862;
t854 = -0.2e1 * qJD(6) * t915 + (t914 * t932 - t875) * qJ(6) + (t915 * t932 + t874) * pkin(5) + t858;
t850 = m(7) * t854 + mrSges(7,1) * t874 - t875 * mrSges(7,3) + t894 * t914 - t915 * t897;
t847 = -m(6) * t858 - t874 * mrSges(6,1) - mrSges(6,2) * t875 - t914 * t895 - t896 * t915 - t850;
t916 = -mrSges(5,2) * t948 + mrSges(5,3) * t933;
t843 = m(5) * t862 + mrSges(5,1) * t928 - mrSges(5,3) * t908 - t911 * t934 + t916 * t948 + t847;
t990 = t978 * t837 - t843 * t975;
t828 = m(4) * t869 - mrSges(4,2) * t964 + mrSges(4,3) * t929 - t924 * t949 - t936 * t965 + t990;
t935 = -mrSges(4,2) * t965 - mrSges(4,3) * t949;
t840 = t1016 * t846 + t974 * t845;
t983 = -m(5) * t866 + t907 * mrSges(5,1) - t908 * mrSges(5,2) + t933 * t916 - t934 * t917 - t840;
t834 = m(4) * t868 + t964 * mrSges(4,1) - t930 * mrSges(4,3) - t950 * t924 + t965 * t935 + t983;
t824 = t970 * t828 + t972 * t834;
t922 = -g(3) * t1007 + t988;
t993 = t979 * t999;
t954 = -mrSges(3,2) * t965 + mrSges(3,3) * t993;
t957 = (-mrSges(3,1) * t979 + mrSges(3,2) * t976) * t999;
t822 = m(3) * t922 + mrSges(3,1) * t964 - mrSges(3,3) * t958 + t954 * t965 - t957 * t994 + t824;
t953 = mrSges(3,1) * t965 - mrSges(3,3) * t994;
t991 = t972 * t828 - t834 * t970;
t823 = m(3) * t923 - mrSges(3,2) * t964 + mrSges(3,3) * t959 - t953 * t965 + t957 * t993 + t991;
t832 = t975 * t837 + t978 * t843;
t831 = m(4) * t910 - t929 * mrSges(4,1) + t930 * mrSges(4,2) + t949 * t935 + t950 * t936 + t832;
t830 = m(3) * t940 - t959 * mrSges(3,1) + t958 * mrSges(3,2) + (t953 * t976 - t954 * t979) * t999 + t831;
t809 = t822 * t1005 + t823 * t1006 - t830 * t971;
t806 = m(2) * t961 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t981 + t809;
t814 = -t822 * t976 + t979 * t823;
t812 = m(2) * t962 - mrSges(2,1) * t981 - qJDD(1) * mrSges(2,2) + t814;
t1004 = t980 * t806 + t977 * t812;
t1003 = t1011 * t914 + t1012 * t915 + t1018 * t932;
t808 = t822 * t1007 + t823 * t1008 + t973 * t830;
t992 = -t806 * t977 + t980 * t812;
t838 = -mrSges(6,1) * t858 - mrSges(7,1) * t854 + mrSges(7,2) * t852 + mrSges(6,3) * t856 - pkin(5) * t850 - t1001 * t932 + t1003 * t915 + t1011 * t906 + t1013 * t875 - t1019 * t874;
t839 = mrSges(6,2) * t858 + mrSges(7,2) * t853 - mrSges(6,3) * t855 - mrSges(7,3) * t854 - qJ(6) * t850 + t1002 * t932 + t1003 * t914 - t1012 * t906 - t1013 * t874 + t1020 * t875;
t899 = Ifges(5,5) * t934 + Ifges(5,6) * t933 + Ifges(5,3) * t948;
t900 = Ifges(5,4) * t934 + Ifges(5,2) * t933 + Ifges(5,6) * t948;
t816 = mrSges(5,2) * t866 - mrSges(5,3) * t862 + Ifges(5,1) * t908 + Ifges(5,4) * t907 + Ifges(5,5) * t928 - pkin(10) * t840 + t1016 * t839 - t974 * t838 + t933 * t899 - t948 * t900;
t901 = Ifges(5,1) * t934 + Ifges(5,4) * t933 + Ifges(5,5) * t948;
t825 = -mrSges(5,1) * t866 + mrSges(5,3) * t863 + Ifges(5,4) * t908 + Ifges(5,2) * t907 + Ifges(5,6) * t928 - pkin(4) * t840 - t934 * t899 + t948 * t901 - t1017;
t918 = Ifges(4,5) * t950 - Ifges(4,6) * t949 + Ifges(4,3) * t965;
t919 = Ifges(4,4) * t950 - Ifges(4,2) * t949 + Ifges(4,6) * t965;
t804 = mrSges(4,2) * t910 - mrSges(4,3) * t868 + Ifges(4,1) * t930 + Ifges(4,4) * t929 + Ifges(4,5) * t964 - pkin(9) * t832 + t816 * t978 - t825 * t975 - t918 * t949 - t919 * t965;
t920 = Ifges(4,1) * t950 - Ifges(4,4) * t949 + Ifges(4,5) * t965;
t982 = mrSges(5,1) * t862 - mrSges(5,2) * t863 + Ifges(5,5) * t908 + Ifges(5,6) * t907 + Ifges(5,3) * t928 + pkin(4) * t847 + pkin(10) * t841 + t1016 * t838 + t974 * t839 + t934 * t900 - t933 * t901;
t815 = -mrSges(4,1) * t910 + mrSges(4,3) * t869 + Ifges(4,4) * t930 + Ifges(4,2) * t929 + Ifges(4,6) * t964 - pkin(3) * t832 - t950 * t918 + t965 * t920 - t982;
t937 = Ifges(3,3) * t965 + (Ifges(3,5) * t976 + Ifges(3,6) * t979) * t999;
t939 = Ifges(3,5) * t965 + (Ifges(3,1) * t976 + Ifges(3,4) * t979) * t999;
t799 = -mrSges(3,1) * t940 + mrSges(3,3) * t923 + Ifges(3,4) * t958 + Ifges(3,2) * t959 + Ifges(3,6) * t964 - pkin(2) * t831 + qJ(3) * t991 + t970 * t804 + t972 * t815 - t937 * t994 + t965 * t939;
t938 = Ifges(3,6) * t965 + (Ifges(3,4) * t976 + Ifges(3,2) * t979) * t999;
t801 = mrSges(3,2) * t940 - mrSges(3,3) * t922 + Ifges(3,1) * t958 + Ifges(3,4) * t959 + Ifges(3,5) * t964 - qJ(3) * t824 + t804 * t972 - t815 * t970 + t937 * t993 - t938 * t965;
t803 = Ifges(3,5) * t958 + Ifges(3,6) * t959 + mrSges(3,1) * t922 - mrSges(3,2) * t923 + Ifges(4,5) * t930 + Ifges(4,6) * t929 + t950 * t919 + t949 * t920 + mrSges(4,1) * t868 - mrSges(4,2) * t869 + t975 * t816 + t978 * t825 + pkin(3) * t983 + pkin(9) * t990 + pkin(2) * t824 + (Ifges(3,3) + Ifges(4,3)) * t964 + (t938 * t976 - t939 * t979) * t999;
t985 = mrSges(2,1) * t961 - mrSges(2,2) * t962 + Ifges(2,3) * qJDD(1) + pkin(1) * t809 + t799 * t1007 + t801 * t1008 + t814 * t1015 + t973 * t803;
t797 = -mrSges(2,2) * g(3) - mrSges(2,3) * t961 + Ifges(2,5) * qJDD(1) - t981 * Ifges(2,6) - t976 * t799 + t979 * t801 + (-t808 * t971 - t809 * t973) * pkin(8);
t796 = mrSges(2,1) * g(3) + mrSges(2,3) * t962 + t981 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t808 - t971 * t803 + (pkin(8) * t814 + t799 * t979 + t801 * t976) * t973;
t1 = [-m(1) * g(1) + t992; -m(1) * g(2) + t1004; (-m(1) - m(2)) * g(3) + t808; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1004 - t977 * t796 + t980 * t797; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t992 + t980 * t796 + t977 * t797; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t985; t985; t803; t831; t982; t1017; t849;];
tauJB  = t1;
