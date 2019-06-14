% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 15:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR12_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR12_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 15:01:22
% EndTime: 2019-05-07 15:02:59
% DurationCPUTime: 80.42s
% Computational Cost: add. (1370745->401), mult. (2994183->519), div. (0->0), fcn. (2420187->14), ass. (0->166)
t1001 = qJD(1) * qJD(2);
t973 = sin(pkin(6));
t979 = sin(qJ(2));
t983 = cos(qJ(2));
t958 = (-qJDD(1) * t983 + t1001 * t979) * t973;
t1012 = cos(qJ(3));
t1011 = pkin(8) * t973;
t975 = cos(pkin(6));
t1010 = t975 * g(3);
t1009 = t973 * t979;
t1008 = t973 * t983;
t1007 = t975 * t979;
t1006 = t975 * t983;
t980 = sin(qJ(1));
t984 = cos(qJ(1));
t964 = t980 * g(1) - g(2) * t984;
t985 = qJD(1) ^ 2;
t953 = qJDD(1) * pkin(1) + t1011 * t985 + t964;
t965 = -g(1) * t984 - g(2) * t980;
t954 = -pkin(1) * t985 + qJDD(1) * t1011 + t965;
t1004 = t953 * t1007 + t983 * t954;
t928 = -g(3) * t1009 + t1004;
t1003 = qJD(1) * t973;
t1000 = t979 * t1003;
t968 = qJD(1) * t975 + qJD(2);
t951 = mrSges(3,1) * t968 - mrSges(3,3) * t1000;
t955 = (-mrSges(3,1) * t983 + mrSges(3,2) * t979) * t1003;
t967 = qJDD(1) * t975 + qJDD(2);
t1002 = qJD(1) * t983;
t956 = (-pkin(2) * t983 - pkin(9) * t979) * t1003;
t966 = t968 ^ 2;
t905 = -t966 * pkin(2) + t967 * pkin(9) + (-g(3) * t979 + t1002 * t956) * t973 + t1004;
t957 = (qJDD(1) * t979 + t1001 * t983) * t973;
t906 = t958 * pkin(2) - t957 * pkin(9) - t1010 + (-t953 + (pkin(2) * t979 - pkin(9) * t983) * t968 * qJD(1)) * t973;
t978 = sin(qJ(3));
t885 = t1012 * t905 + t978 * t906;
t946 = t1000 * t978 - t1012 * t968;
t947 = t1000 * t1012 + t978 * t968;
t929 = pkin(3) * t946 - qJ(4) * t947;
t950 = qJDD(3) + t958;
t999 = t973 * t1002;
t963 = qJD(3) - t999;
t962 = t963 ^ 2;
t875 = -pkin(3) * t962 + qJ(4) * t950 - t929 * t946 + t885;
t927 = -g(3) * t1008 + t1006 * t953 - t979 * t954;
t904 = -t967 * pkin(2) - t966 * pkin(9) + t956 * t1000 - t927;
t925 = qJD(3) * t947 - t1012 * t967 + t957 * t978;
t926 = -t946 * qJD(3) + t1012 * t957 + t978 * t967;
t878 = (t946 * t963 - t926) * qJ(4) + (t947 * t963 + t925) * pkin(3) + t904;
t972 = sin(pkin(12));
t974 = cos(pkin(12));
t935 = t947 * t974 + t963 * t972;
t864 = -0.2e1 * qJD(4) * t935 - t972 * t875 + t974 * t878;
t912 = t926 * t974 + t950 * t972;
t934 = -t947 * t972 + t963 * t974;
t856 = (t934 * t946 - t912) * pkin(10) + (t934 * t935 + t925) * pkin(4) + t864;
t865 = 0.2e1 * qJD(4) * t934 + t974 * t875 + t972 * t878;
t911 = -t926 * t972 + t950 * t974;
t916 = pkin(4) * t946 - pkin(10) * t935;
t933 = t934 ^ 2;
t858 = -pkin(4) * t933 + pkin(10) * t911 - t916 * t946 + t865;
t977 = sin(qJ(5));
t982 = cos(qJ(5));
t851 = t982 * t856 - t977 * t858;
t908 = t934 * t982 - t935 * t977;
t881 = qJD(5) * t908 + t911 * t977 + t912 * t982;
t909 = t934 * t977 + t935 * t982;
t923 = qJDD(5) + t925;
t945 = qJD(5) + t946;
t849 = (t908 * t945 - t881) * pkin(11) + (t908 * t909 + t923) * pkin(5) + t851;
t852 = t977 * t856 + t982 * t858;
t880 = -qJD(5) * t909 + t911 * t982 - t912 * t977;
t895 = pkin(5) * t945 - pkin(11) * t909;
t907 = t908 ^ 2;
t850 = -pkin(5) * t907 + pkin(11) * t880 - t895 * t945 + t852;
t976 = sin(qJ(6));
t981 = cos(qJ(6));
t847 = t849 * t981 - t850 * t976;
t890 = t908 * t981 - t909 * t976;
t863 = qJD(6) * t890 + t880 * t976 + t881 * t981;
t891 = t908 * t976 + t909 * t981;
t872 = -mrSges(7,1) * t890 + mrSges(7,2) * t891;
t944 = qJD(6) + t945;
t882 = -mrSges(7,2) * t944 + mrSges(7,3) * t890;
t918 = qJDD(6) + t923;
t841 = m(7) * t847 + mrSges(7,1) * t918 - mrSges(7,3) * t863 - t872 * t891 + t882 * t944;
t848 = t849 * t976 + t850 * t981;
t862 = -qJD(6) * t891 + t880 * t981 - t881 * t976;
t883 = mrSges(7,1) * t944 - mrSges(7,3) * t891;
t842 = m(7) * t848 - mrSges(7,2) * t918 + mrSges(7,3) * t862 + t872 * t890 - t883 * t944;
t835 = t981 * t841 + t976 * t842;
t892 = -mrSges(6,1) * t908 + mrSges(6,2) * t909;
t893 = -mrSges(6,2) * t945 + mrSges(6,3) * t908;
t833 = m(6) * t851 + mrSges(6,1) * t923 - mrSges(6,3) * t881 - t892 * t909 + t893 * t945 + t835;
t894 = mrSges(6,1) * t945 - mrSges(6,3) * t909;
t995 = -t841 * t976 + t981 * t842;
t834 = m(6) * t852 - mrSges(6,2) * t923 + mrSges(6,3) * t880 + t892 * t908 - t894 * t945 + t995;
t829 = t982 * t833 + t977 * t834;
t913 = -mrSges(5,1) * t934 + mrSges(5,2) * t935;
t994 = -mrSges(5,2) * t946 + mrSges(5,3) * t934;
t827 = m(5) * t864 + t925 * mrSges(5,1) - t912 * mrSges(5,3) - t935 * t913 + t946 * t994 + t829;
t915 = mrSges(5,1) * t946 - mrSges(5,3) * t935;
t996 = -t833 * t977 + t982 * t834;
t828 = m(5) * t865 - mrSges(5,2) * t925 + mrSges(5,3) * t911 + t913 * t934 - t915 * t946 + t996;
t823 = -t827 * t972 + t974 * t828;
t930 = mrSges(4,1) * t946 + mrSges(4,2) * t947;
t937 = mrSges(4,1) * t963 - mrSges(4,3) * t947;
t821 = m(4) * t885 - mrSges(4,2) * t950 - mrSges(4,3) * t925 - t930 * t946 - t937 * t963 + t823;
t884 = t1012 * t906 - t978 * t905;
t874 = -t950 * pkin(3) - t962 * qJ(4) + t947 * t929 + qJDD(4) - t884;
t867 = -t911 * pkin(4) - t933 * pkin(10) + t935 * t916 + t874;
t853 = -t880 * pkin(5) - t907 * pkin(11) + t909 * t895 + t867;
t992 = m(7) * t853 - t862 * mrSges(7,1) + t863 * mrSges(7,2) - t890 * t882 + t891 * t883;
t988 = m(6) * t867 - t880 * mrSges(6,1) + t881 * mrSges(6,2) - t908 * t893 + t909 * t894 + t992;
t845 = m(5) * t874 - t911 * mrSges(5,1) + t912 * mrSges(5,2) + t935 * t915 - t934 * t994 + t988;
t936 = -mrSges(4,2) * t963 - mrSges(4,3) * t946;
t844 = m(4) * t884 + t950 * mrSges(4,1) - t926 * mrSges(4,3) - t947 * t930 + t963 * t936 - t845;
t997 = t1012 * t821 - t844 * t978;
t812 = m(3) * t928 - mrSges(3,2) * t967 - mrSges(3,3) * t958 - t951 * t968 + t955 * t999 + t997;
t815 = t1012 * t844 + t978 * t821;
t941 = -t973 * t953 - t1010;
t952 = -mrSges(3,2) * t968 + mrSges(3,3) * t999;
t814 = m(3) * t941 + t958 * mrSges(3,1) + t957 * mrSges(3,2) + (t951 * t979 - t952 * t983) * t1003 + t815;
t822 = t827 * t974 + t828 * t972;
t989 = -m(4) * t904 - t925 * mrSges(4,1) - mrSges(4,2) * t926 - t946 * t936 - t937 * t947 - t822;
t818 = m(3) * t927 + mrSges(3,1) * t967 - mrSges(3,3) * t957 - t1000 * t955 + t952 * t968 + t989;
t800 = t818 * t1006 + t812 * t1007 - t814 * t973;
t797 = m(2) * t964 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t985 + t800;
t805 = t983 * t812 - t818 * t979;
t803 = m(2) * t965 - mrSges(2,1) * t985 - qJDD(1) * mrSges(2,2) + t805;
t1005 = t984 * t797 + t980 * t803;
t799 = t818 * t1008 + t812 * t1009 + t975 * t814;
t998 = -t797 * t980 + t984 * t803;
t868 = Ifges(7,5) * t891 + Ifges(7,6) * t890 + Ifges(7,3) * t944;
t870 = Ifges(7,1) * t891 + Ifges(7,4) * t890 + Ifges(7,5) * t944;
t836 = -mrSges(7,1) * t853 + mrSges(7,3) * t848 + Ifges(7,4) * t863 + Ifges(7,2) * t862 + Ifges(7,6) * t918 - t868 * t891 + t870 * t944;
t869 = Ifges(7,4) * t891 + Ifges(7,2) * t890 + Ifges(7,6) * t944;
t837 = mrSges(7,2) * t853 - mrSges(7,3) * t847 + Ifges(7,1) * t863 + Ifges(7,4) * t862 + Ifges(7,5) * t918 + t868 * t890 - t869 * t944;
t886 = Ifges(6,5) * t909 + Ifges(6,6) * t908 + Ifges(6,3) * t945;
t888 = Ifges(6,1) * t909 + Ifges(6,4) * t908 + Ifges(6,5) * t945;
t824 = -mrSges(6,1) * t867 + mrSges(6,3) * t852 + Ifges(6,4) * t881 + Ifges(6,2) * t880 + Ifges(6,6) * t923 - pkin(5) * t992 + pkin(11) * t995 + t981 * t836 + t976 * t837 - t909 * t886 + t945 * t888;
t887 = Ifges(6,4) * t909 + Ifges(6,2) * t908 + Ifges(6,6) * t945;
t825 = mrSges(6,2) * t867 - mrSges(6,3) * t851 + Ifges(6,1) * t881 + Ifges(6,4) * t880 + Ifges(6,5) * t923 - pkin(11) * t835 - t836 * t976 + t837 * t981 + t886 * t908 - t887 * t945;
t896 = Ifges(5,5) * t935 + Ifges(5,6) * t934 + Ifges(5,3) * t946;
t898 = Ifges(5,1) * t935 + Ifges(5,4) * t934 + Ifges(5,5) * t946;
t807 = -mrSges(5,1) * t874 + mrSges(5,3) * t865 + Ifges(5,4) * t912 + Ifges(5,2) * t911 + Ifges(5,6) * t925 - pkin(4) * t988 + pkin(10) * t996 + t982 * t824 + t977 * t825 - t935 * t896 + t946 * t898;
t897 = Ifges(5,4) * t935 + Ifges(5,2) * t934 + Ifges(5,6) * t946;
t808 = mrSges(5,2) * t874 - mrSges(5,3) * t864 + Ifges(5,1) * t912 + Ifges(5,4) * t911 + Ifges(5,5) * t925 - pkin(10) * t829 - t824 * t977 + t825 * t982 + t896 * t934 - t897 * t946;
t919 = Ifges(4,5) * t947 - Ifges(4,6) * t946 + Ifges(4,3) * t963;
t920 = Ifges(4,4) * t947 - Ifges(4,2) * t946 + Ifges(4,6) * t963;
t795 = mrSges(4,2) * t904 - mrSges(4,3) * t884 + Ifges(4,1) * t926 - Ifges(4,4) * t925 + Ifges(4,5) * t950 - qJ(4) * t822 - t807 * t972 + t808 * t974 - t919 * t946 - t920 * t963;
t921 = Ifges(4,1) * t947 - Ifges(4,4) * t946 + Ifges(4,5) * t963;
t990 = -mrSges(7,1) * t847 + mrSges(7,2) * t848 - Ifges(7,5) * t863 - Ifges(7,6) * t862 - Ifges(7,3) * t918 - t891 * t869 + t890 * t870;
t986 = mrSges(6,1) * t851 - mrSges(6,2) * t852 + Ifges(6,5) * t881 + Ifges(6,6) * t880 + Ifges(6,3) * t923 + pkin(5) * t835 + t909 * t887 - t908 * t888 - t990;
t806 = -t986 - pkin(4) * t829 - pkin(3) * t822 + (-Ifges(5,3) - Ifges(4,2)) * t925 + t963 * t921 - t947 * t919 + Ifges(4,6) * t950 + t934 * t898 - t935 * t897 + Ifges(4,4) * t926 - Ifges(5,6) * t911 - Ifges(5,5) * t912 - mrSges(4,1) * t904 + mrSges(4,3) * t885 - mrSges(5,1) * t864 + mrSges(5,2) * t865;
t939 = Ifges(3,6) * t968 + (Ifges(3,4) * t979 + Ifges(3,2) * t983) * t1003;
t940 = Ifges(3,5) * t968 + (Ifges(3,1) * t979 + Ifges(3,4) * t983) * t1003;
t790 = Ifges(3,5) * t957 - Ifges(3,6) * t958 + Ifges(3,3) * t967 + mrSges(3,1) * t927 - mrSges(3,2) * t928 + t978 * t795 + t1012 * t806 + pkin(2) * t989 + pkin(9) * t997 + (t939 * t979 - t940 * t983) * t1003;
t938 = Ifges(3,3) * t968 + (Ifges(3,5) * t979 + Ifges(3,6) * t983) * t1003;
t792 = mrSges(3,2) * t941 - mrSges(3,3) * t927 + Ifges(3,1) * t957 - Ifges(3,4) * t958 + Ifges(3,5) * t967 - pkin(9) * t815 + t1012 * t795 - t978 * t806 + t938 * t999 - t968 * t939;
t987 = mrSges(4,1) * t884 - mrSges(4,2) * t885 + Ifges(4,5) * t926 - Ifges(4,6) * t925 + Ifges(4,3) * t950 - pkin(3) * t845 + qJ(4) * t823 + t974 * t807 + t972 * t808 + t947 * t920 + t946 * t921;
t794 = -mrSges(3,1) * t941 + mrSges(3,3) * t928 + Ifges(3,4) * t957 - Ifges(3,2) * t958 + Ifges(3,6) * t967 - pkin(2) * t815 - t1000 * t938 + t968 * t940 - t987;
t991 = mrSges(2,1) * t964 - mrSges(2,2) * t965 + Ifges(2,3) * qJDD(1) + pkin(1) * t800 + t794 * t1008 + t792 * t1009 + t805 * t1011 + t975 * t790;
t788 = -mrSges(2,2) * g(3) - mrSges(2,3) * t964 + Ifges(2,5) * qJDD(1) - t985 * Ifges(2,6) + t983 * t792 - t979 * t794 + (-t799 * t973 - t800 * t975) * pkin(8);
t787 = mrSges(2,1) * g(3) + mrSges(2,3) * t965 + t985 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t799 - t973 * t790 + (pkin(8) * t805 + t792 * t979 + t794 * t983) * t975;
t1 = [-m(1) * g(1) + t998; -m(1) * g(2) + t1005; (-m(1) - m(2)) * g(3) + t799; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1005 - t980 * t787 + t984 * t788; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t998 + t984 * t787 + t980 * t788; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t991; t991; t790; t987; t845; t986; -t990;];
tauJB  = t1;
