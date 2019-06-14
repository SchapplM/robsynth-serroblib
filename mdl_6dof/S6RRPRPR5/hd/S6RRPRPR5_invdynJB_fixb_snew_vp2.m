% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 14:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:58:41
% EndTime: 2019-05-06 13:59:41
% DurationCPUTime: 57.95s
% Computational Cost: add. (901552->398), mult. (2388875->519), div. (0->0), fcn. (1903976->14), ass. (0->163)
t1019 = -2 * qJD(3);
t980 = sin(pkin(6));
t1010 = qJD(1) * t980;
t979 = sin(pkin(11));
t982 = cos(pkin(11));
t986 = sin(qJ(2));
t989 = cos(qJ(2));
t955 = (t979 * t986 - t982 * t989) * t1010;
t991 = qJD(1) ^ 2;
t1016 = t980 ^ 2 * t991;
t1008 = qJD(1) * qJD(2);
t964 = (qJDD(1) * t986 + t1008 * t989) * t980;
t983 = cos(pkin(6));
t972 = qJDD(1) * t983 + qJDD(2);
t973 = qJD(1) * t983 + qJD(2);
t1012 = t983 * t989;
t1017 = pkin(8) * t980;
t987 = sin(qJ(1));
t990 = cos(qJ(1));
t969 = t987 * g(1) - g(2) * t990;
t961 = qJDD(1) * pkin(1) + t1017 * t991 + t969;
t970 = -g(1) * t990 - g(2) * t987;
t962 = -pkin(1) * t991 + qJDD(1) * t1017 + t970;
t999 = t961 * t1012 - t986 * t962;
t903 = t972 * pkin(2) - t964 * qJ(3) + (pkin(2) * t986 * t1016 + (qJ(3) * qJD(1) * t973 - g(3)) * t980) * t989 + t999;
t1007 = t989 ^ 2 * t1016;
t1013 = t983 * t986;
t1015 = t980 * t986;
t932 = -g(3) * t1015 + t961 * t1013 + t989 * t962;
t1006 = t986 * t1010;
t958 = pkin(2) * t973 - qJ(3) * t1006;
t965 = (qJDD(1) * t989 - t1008 * t986) * t980;
t906 = -pkin(2) * t1007 + qJ(3) * t965 - t958 * t973 + t932;
t956 = (t979 * t989 + t982 * t986) * t1010;
t883 = t1019 * t956 + t982 * t903 - t979 * t906;
t1018 = cos(qJ(4));
t1014 = t980 * t989;
t884 = t1019 * t955 + t979 * t903 + t982 * t906;
t934 = pkin(3) * t955 - pkin(9) * t956;
t971 = t973 ^ 2;
t877 = -pkin(3) * t971 + pkin(9) * t972 - t934 * t955 + t884;
t947 = -t983 * g(3) - t980 * t961;
t916 = -t965 * pkin(2) - qJ(3) * t1007 + t958 * t1006 + qJDD(3) + t947;
t937 = -t964 * t979 + t965 * t982;
t938 = t964 * t982 + t965 * t979;
t886 = (t955 * t973 - t938) * pkin(9) + (t956 * t973 - t937) * pkin(3) + t916;
t985 = sin(qJ(4));
t870 = t1018 * t877 + t985 * t886;
t940 = -t1018 * t973 + t956 * t985;
t941 = t1018 * t956 + t985 * t973;
t917 = pkin(4) * t940 - qJ(5) * t941;
t936 = qJDD(4) - t937;
t954 = qJD(4) + t955;
t953 = t954 ^ 2;
t865 = -pkin(4) * t953 + qJ(5) * t936 - t917 * t940 + t870;
t876 = -t972 * pkin(3) - t971 * pkin(9) + t956 * t934 - t883;
t913 = qJD(4) * t941 - t1018 * t972 + t938 * t985;
t914 = -t940 * qJD(4) + t1018 * t938 + t985 * t972;
t868 = (t940 * t954 - t914) * qJ(5) + (t941 * t954 + t913) * pkin(4) + t876;
t978 = sin(pkin(12));
t981 = cos(pkin(12));
t924 = t941 * t981 + t954 * t978;
t860 = -0.2e1 * qJD(5) * t924 - t978 * t865 + t981 * t868;
t895 = t914 * t981 + t936 * t978;
t923 = -t941 * t978 + t954 * t981;
t858 = (t923 * t940 - t895) * pkin(10) + (t923 * t924 + t913) * pkin(5) + t860;
t861 = 0.2e1 * qJD(5) * t923 + t981 * t865 + t978 * t868;
t894 = -t914 * t978 + t936 * t981;
t905 = pkin(5) * t940 - pkin(10) * t924;
t922 = t923 ^ 2;
t859 = -pkin(5) * t922 + pkin(10) * t894 - t905 * t940 + t861;
t984 = sin(qJ(6));
t988 = cos(qJ(6));
t856 = t858 * t988 - t859 * t984;
t896 = t923 * t988 - t924 * t984;
t873 = qJD(6) * t896 + t894 * t984 + t895 * t988;
t897 = t923 * t984 + t924 * t988;
t882 = -mrSges(7,1) * t896 + mrSges(7,2) * t897;
t939 = qJD(6) + t940;
t887 = -mrSges(7,2) * t939 + mrSges(7,3) * t896;
t912 = qJDD(6) + t913;
t853 = m(7) * t856 + mrSges(7,1) * t912 - mrSges(7,3) * t873 - t882 * t897 + t887 * t939;
t857 = t858 * t984 + t859 * t988;
t872 = -qJD(6) * t897 + t894 * t988 - t895 * t984;
t888 = mrSges(7,1) * t939 - mrSges(7,3) * t897;
t854 = m(7) * t857 - mrSges(7,2) * t912 + mrSges(7,3) * t872 + t882 * t896 - t888 * t939;
t845 = t988 * t853 + t984 * t854;
t898 = -mrSges(6,1) * t923 + mrSges(6,2) * t924;
t998 = -mrSges(6,2) * t940 + mrSges(6,3) * t923;
t843 = m(6) * t860 + t913 * mrSges(6,1) - t895 * mrSges(6,3) - t924 * t898 + t940 * t998 + t845;
t1001 = -t853 * t984 + t988 * t854;
t904 = mrSges(6,1) * t940 - mrSges(6,3) * t924;
t844 = m(6) * t861 - mrSges(6,2) * t913 + mrSges(6,3) * t894 + t898 * t923 - t904 * t940 + t1001;
t841 = -t843 * t978 + t981 * t844;
t918 = mrSges(5,1) * t940 + mrSges(5,2) * t941;
t926 = mrSges(5,1) * t954 - mrSges(5,3) * t941;
t839 = m(5) * t870 - mrSges(5,2) * t936 - mrSges(5,3) * t913 - t918 * t940 - t926 * t954 + t841;
t869 = t1018 * t886 - t985 * t877;
t864 = -t936 * pkin(4) - t953 * qJ(5) + t941 * t917 + qJDD(5) - t869;
t862 = -t894 * pkin(5) - t922 * pkin(10) + t924 * t905 + t864;
t996 = m(7) * t862 - t872 * mrSges(7,1) + mrSges(7,2) * t873 - t896 * t887 + t888 * t897;
t855 = m(6) * t864 - t894 * mrSges(6,1) + mrSges(6,2) * t895 + t904 * t924 - t923 * t998 + t996;
t925 = -mrSges(5,2) * t954 - mrSges(5,3) * t940;
t852 = m(5) * t869 + mrSges(5,1) * t936 - mrSges(5,3) * t914 - t918 * t941 + t925 * t954 - t855;
t1002 = t1018 * t839 - t852 * t985;
t933 = mrSges(4,1) * t955 + mrSges(4,2) * t956;
t943 = mrSges(4,1) * t973 - mrSges(4,3) * t956;
t828 = m(4) * t884 - mrSges(4,2) * t972 + mrSges(4,3) * t937 - t933 * t955 - t943 * t973 + t1002;
t942 = -mrSges(4,2) * t973 - mrSges(4,3) * t955;
t840 = t843 * t981 + t844 * t978;
t994 = -m(5) * t876 - t913 * mrSges(5,1) - mrSges(5,2) * t914 - t940 * t925 - t926 * t941 - t840;
t836 = m(4) * t883 + mrSges(4,1) * t972 - mrSges(4,3) * t938 - t933 * t956 + t942 * t973 + t994;
t824 = t979 * t828 + t982 * t836;
t931 = -g(3) * t1014 + t999;
t1005 = t989 * t1010;
t960 = -mrSges(3,2) * t973 + mrSges(3,3) * t1005;
t963 = (-mrSges(3,1) * t989 + mrSges(3,2) * t986) * t1010;
t822 = m(3) * t931 + mrSges(3,1) * t972 - mrSges(3,3) * t964 - t1006 * t963 + t960 * t973 + t824;
t1003 = t982 * t828 - t836 * t979;
t959 = mrSges(3,1) * t973 - mrSges(3,3) * t1006;
t823 = m(3) * t932 - mrSges(3,2) * t972 + mrSges(3,3) * t965 + t1005 * t963 - t959 * t973 + t1003;
t832 = t1018 * t852 + t985 * t839;
t831 = m(4) * t916 - t937 * mrSges(4,1) + t938 * mrSges(4,2) + t955 * t942 + t956 * t943 + t832;
t830 = m(3) * t947 - t965 * mrSges(3,1) + t964 * mrSges(3,2) + (t959 * t986 - t960 * t989) * t1010 + t831;
t809 = t822 * t1012 + t823 * t1013 - t830 * t980;
t806 = m(2) * t969 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t991 + t809;
t815 = -t822 * t986 + t989 * t823;
t813 = m(2) * t970 - mrSges(2,1) * t991 - qJDD(1) * mrSges(2,2) + t815;
t1011 = t990 * t806 + t987 * t813;
t808 = t822 * t1014 + t823 * t1015 + t983 * t830;
t1004 = -t806 * t987 + t990 * t813;
t878 = Ifges(7,5) * t897 + Ifges(7,6) * t896 + Ifges(7,3) * t939;
t880 = Ifges(7,1) * t897 + Ifges(7,4) * t896 + Ifges(7,5) * t939;
t846 = -mrSges(7,1) * t862 + mrSges(7,3) * t857 + Ifges(7,4) * t873 + Ifges(7,2) * t872 + Ifges(7,6) * t912 - t878 * t897 + t880 * t939;
t879 = Ifges(7,4) * t897 + Ifges(7,2) * t896 + Ifges(7,6) * t939;
t847 = mrSges(7,2) * t862 - mrSges(7,3) * t856 + Ifges(7,1) * t873 + Ifges(7,4) * t872 + Ifges(7,5) * t912 + t878 * t896 - t879 * t939;
t889 = Ifges(6,5) * t924 + Ifges(6,6) * t923 + Ifges(6,3) * t940;
t891 = Ifges(6,1) * t924 + Ifges(6,4) * t923 + Ifges(6,5) * t940;
t833 = -mrSges(6,1) * t864 + mrSges(6,3) * t861 + Ifges(6,4) * t895 + Ifges(6,2) * t894 + Ifges(6,6) * t913 - pkin(5) * t996 + pkin(10) * t1001 + t988 * t846 + t984 * t847 - t924 * t889 + t940 * t891;
t890 = Ifges(6,4) * t924 + Ifges(6,2) * t923 + Ifges(6,6) * t940;
t834 = mrSges(6,2) * t864 - mrSges(6,3) * t860 + Ifges(6,1) * t895 + Ifges(6,4) * t894 + Ifges(6,5) * t913 - pkin(10) * t845 - t846 * t984 + t847 * t988 + t889 * t923 - t890 * t940;
t907 = Ifges(5,5) * t941 - Ifges(5,6) * t940 + Ifges(5,3) * t954;
t908 = Ifges(5,4) * t941 - Ifges(5,2) * t940 + Ifges(5,6) * t954;
t816 = mrSges(5,2) * t876 - mrSges(5,3) * t869 + Ifges(5,1) * t914 - Ifges(5,4) * t913 + Ifges(5,5) * t936 - qJ(5) * t840 - t833 * t978 + t834 * t981 - t907 * t940 - t908 * t954;
t909 = Ifges(5,1) * t941 - Ifges(5,4) * t940 + Ifges(5,5) * t954;
t993 = mrSges(7,1) * t856 - mrSges(7,2) * t857 + Ifges(7,5) * t873 + Ifges(7,6) * t872 + Ifges(7,3) * t912 + t897 * t879 - t896 * t880;
t825 = -Ifges(6,6) * t894 - Ifges(6,5) * t895 + t923 * t891 - t924 * t890 - t993 + Ifges(5,4) * t914 - mrSges(5,1) * t876 - pkin(5) * t845 - pkin(4) * t840 + mrSges(6,2) * t861 - mrSges(6,1) * t860 + (-Ifges(5,2) - Ifges(6,3)) * t913 + mrSges(5,3) * t870 + Ifges(5,6) * t936 - t941 * t907 + t954 * t909;
t927 = Ifges(4,5) * t956 - Ifges(4,6) * t955 + Ifges(4,3) * t973;
t928 = Ifges(4,4) * t956 - Ifges(4,2) * t955 + Ifges(4,6) * t973;
t804 = mrSges(4,2) * t916 - mrSges(4,3) * t883 + Ifges(4,1) * t938 + Ifges(4,4) * t937 + Ifges(4,5) * t972 - pkin(9) * t832 + t1018 * t816 - t985 * t825 - t955 * t927 - t973 * t928;
t929 = Ifges(4,1) * t956 - Ifges(4,4) * t955 + Ifges(4,5) * t973;
t992 = mrSges(5,1) * t869 - mrSges(5,2) * t870 + Ifges(5,5) * t914 - Ifges(5,6) * t913 + Ifges(5,3) * t936 - pkin(4) * t855 + qJ(5) * t841 + t981 * t833 + t978 * t834 + t941 * t908 + t940 * t909;
t810 = -mrSges(4,1) * t916 + mrSges(4,3) * t884 + Ifges(4,4) * t938 + Ifges(4,2) * t937 + Ifges(4,6) * t972 - pkin(3) * t832 - t956 * t927 + t973 * t929 - t992;
t944 = Ifges(3,3) * t973 + (Ifges(3,5) * t986 + Ifges(3,6) * t989) * t1010;
t946 = Ifges(3,5) * t973 + (Ifges(3,1) * t986 + Ifges(3,4) * t989) * t1010;
t799 = -mrSges(3,1) * t947 + mrSges(3,3) * t932 + Ifges(3,4) * t964 + Ifges(3,2) * t965 + Ifges(3,6) * t972 - pkin(2) * t831 + qJ(3) * t1003 - t1006 * t944 + t979 * t804 + t982 * t810 + t973 * t946;
t945 = Ifges(3,6) * t973 + (Ifges(3,4) * t986 + Ifges(3,2) * t989) * t1010;
t801 = mrSges(3,2) * t947 - mrSges(3,3) * t931 + Ifges(3,1) * t964 + Ifges(3,4) * t965 + Ifges(3,5) * t972 - qJ(3) * t824 + t1005 * t944 + t804 * t982 - t810 * t979 - t945 * t973;
t803 = Ifges(3,5) * t964 + Ifges(3,6) * t965 + mrSges(3,1) * t931 - mrSges(3,2) * t932 + Ifges(4,5) * t938 + Ifges(4,6) * t937 + t956 * t928 + t955 * t929 + mrSges(4,1) * t883 - mrSges(4,2) * t884 + t985 * t816 + t1018 * t825 + pkin(3) * t994 + pkin(9) * t1002 + pkin(2) * t824 + (Ifges(3,3) + Ifges(4,3)) * t972 + (t945 * t986 - t946 * t989) * t1010;
t995 = mrSges(2,1) * t969 - mrSges(2,2) * t970 + Ifges(2,3) * qJDD(1) + pkin(1) * t809 + t799 * t1014 + t801 * t1015 + t815 * t1017 + t983 * t803;
t797 = -mrSges(2,2) * g(3) - mrSges(2,3) * t969 + Ifges(2,5) * qJDD(1) - t991 * Ifges(2,6) - t986 * t799 + t989 * t801 + (-t808 * t980 - t809 * t983) * pkin(8);
t796 = mrSges(2,1) * g(3) + mrSges(2,3) * t970 + t991 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t808 - t980 * t803 + (pkin(8) * t815 + t799 * t989 + t801 * t986) * t983;
t1 = [-m(1) * g(1) + t1004; -m(1) * g(2) + t1011; (-m(1) - m(2)) * g(3) + t808; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1011 - t987 * t796 + t990 * t797; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1004 + t990 * t796 + t987 * t797; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t995; t995; t803; t831; t992; t855; t993;];
tauJB  = t1;
