% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR9
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
% Datum: 2019-05-06 15:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:15:49
% EndTime: 2019-05-06 15:16:57
% DurationCPUTime: 64.69s
% Computational Cost: add. (1053613->399), mult. (2412030->519), div. (0->0), fcn. (1961973->14), ass. (0->163)
t1007 = cos(qJ(4));
t972 = sin(pkin(6));
t1006 = pkin(8) * t972;
t975 = cos(pkin(6));
t1005 = g(3) * t975;
t978 = sin(qJ(2));
t1004 = t972 * t978;
t981 = cos(qJ(2));
t1003 = t972 * t981;
t1002 = t975 * t978;
t1001 = t975 * t981;
t979 = sin(qJ(1));
t982 = cos(qJ(1));
t961 = t979 * g(1) - g(2) * t982;
t983 = qJD(1) ^ 2;
t952 = qJDD(1) * pkin(1) + t1006 * t983 + t961;
t962 = -g(1) * t982 - g(2) * t979;
t996 = qJDD(1) * t972;
t953 = -pkin(1) * t983 + pkin(8) * t996 + t962;
t999 = t952 * t1002 + t981 * t953;
t923 = -g(3) * t1004 + t999;
t966 = qJD(1) * t975 + qJD(2);
t998 = qJD(1) * t972;
t995 = t978 * t998;
t950 = mrSges(3,1) * t966 - mrSges(3,3) * t995;
t955 = (-mrSges(3,1) * t981 + mrSges(3,2) * t978) * t998;
t957 = -qJD(2) * t995 + t981 * t996;
t965 = qJDD(1) * t975 + qJDD(2);
t954 = (-pkin(2) * t981 - qJ(3) * t978) * t998;
t964 = t966 ^ 2;
t997 = qJD(1) * t981;
t908 = -pkin(2) * t964 + qJ(3) * t965 + (-g(3) * t978 + t954 * t997) * t972 + t999;
t956 = (qJD(2) * t997 + qJDD(1) * t978) * t972;
t909 = -pkin(2) * t957 - t1005 - qJ(3) * t956 + (-t952 + (pkin(2) * t978 - qJ(3) * t981) * t966 * qJD(1)) * t972;
t971 = sin(pkin(11));
t974 = cos(pkin(11));
t946 = t966 * t971 + t974 * t995;
t876 = -0.2e1 * qJD(3) * t946 - t908 * t971 + t974 * t909;
t933 = t956 * t974 + t965 * t971;
t945 = t966 * t974 - t971 * t995;
t994 = t972 * t997;
t867 = (-t945 * t994 - t933) * pkin(9) + (t945 * t946 - t957) * pkin(3) + t876;
t877 = 0.2e1 * qJD(3) * t945 + t974 * t908 + t971 * t909;
t932 = -t956 * t971 + t965 * t974;
t934 = -pkin(3) * t994 - pkin(9) * t946;
t944 = t945 ^ 2;
t874 = -pkin(3) * t944 + pkin(9) * t932 + t934 * t994 + t877;
t977 = sin(qJ(4));
t859 = t1007 * t874 + t977 * t867;
t925 = -t1007 * t945 + t946 * t977;
t926 = t1007 * t946 + t977 * t945;
t902 = pkin(4) * t925 - qJ(5) * t926;
t949 = qJDD(4) - t957;
t960 = qJD(4) - t994;
t959 = t960 ^ 2;
t857 = -pkin(4) * t959 + qJ(5) * t949 - t902 * t925 + t859;
t922 = -g(3) * t1003 + t1001 * t952 - t978 * t953;
t907 = -pkin(2) * t965 - qJ(3) * t964 + t954 * t995 + qJDD(3) - t922;
t878 = -pkin(3) * t932 - pkin(9) * t944 + t946 * t934 + t907;
t894 = qJD(4) * t926 - t1007 * t932 + t933 * t977;
t895 = -t925 * qJD(4) + t1007 * t933 + t977 * t932;
t862 = (t925 * t960 - t895) * qJ(5) + (t926 * t960 + t894) * pkin(4) + t878;
t970 = sin(pkin(12));
t973 = cos(pkin(12));
t914 = t926 * t973 + t960 * t970;
t852 = -0.2e1 * qJD(5) * t914 - t857 * t970 + t973 * t862;
t888 = t895 * t973 + t949 * t970;
t913 = -t926 * t970 + t960 * t973;
t850 = (t913 * t925 - t888) * pkin(10) + (t913 * t914 + t894) * pkin(5) + t852;
t853 = 0.2e1 * qJD(5) * t913 + t973 * t857 + t970 * t862;
t887 = -t895 * t970 + t949 * t973;
t897 = pkin(5) * t925 - pkin(10) * t914;
t912 = t913 ^ 2;
t851 = -pkin(5) * t912 + pkin(10) * t887 - t897 * t925 + t853;
t976 = sin(qJ(6));
t980 = cos(qJ(6));
t848 = t850 * t980 - t851 * t976;
t889 = t913 * t980 - t914 * t976;
t865 = qJD(6) * t889 + t887 * t976 + t888 * t980;
t890 = t913 * t976 + t914 * t980;
t875 = -mrSges(7,1) * t889 + mrSges(7,2) * t890;
t924 = qJD(6) + t925;
t879 = -mrSges(7,2) * t924 + mrSges(7,3) * t889;
t893 = qJDD(6) + t894;
t845 = m(7) * t848 + mrSges(7,1) * t893 - mrSges(7,3) * t865 - t875 * t890 + t879 * t924;
t849 = t850 * t976 + t851 * t980;
t864 = -qJD(6) * t890 + t887 * t980 - t888 * t976;
t880 = mrSges(7,1) * t924 - mrSges(7,3) * t890;
t846 = m(7) * t849 - mrSges(7,2) * t893 + mrSges(7,3) * t864 + t875 * t889 - t880 * t924;
t837 = t980 * t845 + t976 * t846;
t891 = -mrSges(6,1) * t913 + mrSges(6,2) * t914;
t989 = -mrSges(6,2) * t925 + mrSges(6,3) * t913;
t835 = m(6) * t852 + t894 * mrSges(6,1) - t888 * mrSges(6,3) - t914 * t891 + t925 * t989 + t837;
t896 = mrSges(6,1) * t925 - mrSges(6,3) * t914;
t990 = -t845 * t976 + t980 * t846;
t836 = m(6) * t853 - mrSges(6,2) * t894 + mrSges(6,3) * t887 + t891 * t913 - t896 * t925 + t990;
t831 = -t835 * t970 + t973 * t836;
t903 = mrSges(5,1) * t925 + mrSges(5,2) * t926;
t916 = mrSges(5,1) * t960 - mrSges(5,3) * t926;
t828 = m(5) * t859 - mrSges(5,2) * t949 - mrSges(5,3) * t894 - t903 * t925 - t916 * t960 + t831;
t858 = t1007 * t867 - t977 * t874;
t856 = -t949 * pkin(4) - t959 * qJ(5) + t926 * t902 + qJDD(5) - t858;
t854 = -t887 * pkin(5) - t912 * pkin(10) + t914 * t897 + t856;
t988 = m(7) * t854 - t864 * mrSges(7,1) + mrSges(7,2) * t865 - t889 * t879 + t880 * t890;
t847 = m(6) * t856 - t887 * mrSges(6,1) + mrSges(6,2) * t888 + t896 * t914 - t913 * t989 + t988;
t915 = -mrSges(5,2) * t960 - mrSges(5,3) * t925;
t841 = m(5) * t858 + mrSges(5,1) * t949 - mrSges(5,3) * t895 - t903 * t926 + t915 * t960 - t847;
t820 = t1007 * t841 + t977 * t828;
t927 = -mrSges(4,1) * t945 + mrSges(4,2) * t946;
t930 = mrSges(4,2) * t994 + mrSges(4,3) * t945;
t818 = m(4) * t876 - mrSges(4,1) * t957 - mrSges(4,3) * t933 - t927 * t946 - t930 * t994 + t820;
t931 = -mrSges(4,1) * t994 - mrSges(4,3) * t946;
t991 = t1007 * t828 - t841 * t977;
t819 = m(4) * t877 + mrSges(4,2) * t957 + mrSges(4,3) * t932 + t927 * t945 + t931 * t994 + t991;
t992 = -t818 * t971 + t974 * t819;
t810 = m(3) * t923 - mrSges(3,2) * t965 + mrSges(3,3) * t957 - t950 * t966 + t955 * t994 + t992;
t813 = t974 * t818 + t971 * t819;
t938 = -t952 * t972 - t1005;
t951 = -mrSges(3,2) * t966 + mrSges(3,3) * t994;
t812 = m(3) * t938 - mrSges(3,1) * t957 + mrSges(3,2) * t956 + (t950 * t978 - t951 * t981) * t998 + t813;
t830 = t973 * t835 + t970 * t836;
t986 = m(5) * t878 + t894 * mrSges(5,1) + mrSges(5,2) * t895 + t925 * t915 + t916 * t926 + t830;
t829 = m(4) * t907 - t932 * mrSges(4,1) + mrSges(4,2) * t933 - t945 * t930 + t931 * t946 + t986;
t825 = m(3) * t922 + mrSges(3,1) * t965 - mrSges(3,3) * t956 + t951 * t966 - t955 * t995 - t829;
t800 = t825 * t1001 + t810 * t1002 - t812 * t972;
t797 = m(2) * t961 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t983 + t800;
t805 = t981 * t810 - t825 * t978;
t803 = m(2) * t962 - mrSges(2,1) * t983 - qJDD(1) * mrSges(2,2) + t805;
t1000 = t982 * t797 + t979 * t803;
t799 = t825 * t1003 + t810 * t1004 + t975 * t812;
t993 = -t797 * t979 + t982 * t803;
t868 = Ifges(7,5) * t890 + Ifges(7,6) * t889 + Ifges(7,3) * t924;
t870 = Ifges(7,1) * t890 + Ifges(7,4) * t889 + Ifges(7,5) * t924;
t838 = -mrSges(7,1) * t854 + mrSges(7,3) * t849 + Ifges(7,4) * t865 + Ifges(7,2) * t864 + Ifges(7,6) * t893 - t868 * t890 + t870 * t924;
t869 = Ifges(7,4) * t890 + Ifges(7,2) * t889 + Ifges(7,6) * t924;
t839 = mrSges(7,2) * t854 - mrSges(7,3) * t848 + Ifges(7,1) * t865 + Ifges(7,4) * t864 + Ifges(7,5) * t893 + t868 * t889 - t869 * t924;
t881 = Ifges(6,5) * t914 + Ifges(6,6) * t913 + Ifges(6,3) * t925;
t883 = Ifges(6,1) * t914 + Ifges(6,4) * t913 + Ifges(6,5) * t925;
t821 = -mrSges(6,1) * t856 + mrSges(6,3) * t853 + Ifges(6,4) * t888 + Ifges(6,2) * t887 + Ifges(6,6) * t894 - pkin(5) * t988 + pkin(10) * t990 + t980 * t838 + t976 * t839 - t914 * t881 + t925 * t883;
t882 = Ifges(6,4) * t914 + Ifges(6,2) * t913 + Ifges(6,6) * t925;
t822 = mrSges(6,2) * t856 - mrSges(6,3) * t852 + Ifges(6,1) * t888 + Ifges(6,4) * t887 + Ifges(6,5) * t894 - pkin(10) * t837 - t838 * t976 + t839 * t980 + t881 * t913 - t882 * t925;
t898 = Ifges(5,5) * t926 - Ifges(5,6) * t925 + Ifges(5,3) * t960;
t899 = Ifges(5,4) * t926 - Ifges(5,2) * t925 + Ifges(5,6) * t960;
t806 = mrSges(5,2) * t878 - mrSges(5,3) * t858 + Ifges(5,1) * t895 - Ifges(5,4) * t894 + Ifges(5,5) * t949 - qJ(5) * t830 - t821 * t970 + t822 * t973 - t898 * t925 - t899 * t960;
t900 = Ifges(5,1) * t926 - Ifges(5,4) * t925 + Ifges(5,5) * t960;
t985 = mrSges(7,1) * t848 - mrSges(7,2) * t849 + Ifges(7,5) * t865 + Ifges(7,6) * t864 + Ifges(7,3) * t893 + t890 * t869 - t889 * t870;
t814 = -t985 + (-Ifges(5,2) - Ifges(6,3)) * t894 - pkin(4) * t830 + t960 * t900 + Ifges(5,6) * t949 - t926 * t898 - t914 * t882 + t913 * t883 + Ifges(5,4) * t895 - Ifges(6,6) * t887 - Ifges(6,5) * t888 - mrSges(5,1) * t878 + mrSges(5,3) * t859 + mrSges(6,2) * t853 - mrSges(6,1) * t852 - pkin(5) * t837;
t917 = Ifges(4,5) * t946 + Ifges(4,6) * t945 - Ifges(4,3) * t994;
t919 = Ifges(4,1) * t946 + Ifges(4,4) * t945 - Ifges(4,5) * t994;
t792 = -mrSges(4,1) * t907 + mrSges(4,3) * t877 + Ifges(4,4) * t933 + Ifges(4,2) * t932 - Ifges(4,6) * t957 - pkin(3) * t986 + pkin(9) * t991 + t1007 * t814 + t977 * t806 - t946 * t917 - t919 * t994;
t918 = Ifges(4,4) * t946 + Ifges(4,2) * t945 - Ifges(4,6) * t994;
t795 = mrSges(4,2) * t907 - mrSges(4,3) * t876 + Ifges(4,1) * t933 + Ifges(4,4) * t932 - Ifges(4,5) * t957 - pkin(9) * t820 + t1007 * t806 - t977 * t814 + t945 * t917 + t918 * t994;
t936 = Ifges(3,6) * t966 + (Ifges(3,4) * t978 + Ifges(3,2) * t981) * t998;
t937 = Ifges(3,5) * t966 + (Ifges(3,1) * t978 + Ifges(3,4) * t981) * t998;
t789 = Ifges(3,5) * t956 + Ifges(3,6) * t957 + Ifges(3,3) * t965 + mrSges(3,1) * t922 - mrSges(3,2) * t923 + t971 * t795 + t974 * t792 - pkin(2) * t829 + qJ(3) * t992 + (t936 * t978 - t937 * t981) * t998;
t935 = Ifges(3,3) * t966 + (Ifges(3,5) * t978 + Ifges(3,6) * t981) * t998;
t791 = mrSges(3,2) * t938 - mrSges(3,3) * t922 + Ifges(3,1) * t956 + Ifges(3,4) * t957 + Ifges(3,5) * t965 - qJ(3) * t813 - t792 * t971 + t795 * t974 + t935 * t994 - t936 * t966;
t984 = mrSges(5,1) * t858 - mrSges(5,2) * t859 + Ifges(5,5) * t895 - Ifges(5,6) * t894 + Ifges(5,3) * t949 - pkin(4) * t847 + qJ(5) * t831 + t973 * t821 + t970 * t822 + t926 * t899 + t925 * t900;
t794 = -pkin(3) * t820 - t984 + (Ifges(3,2) + Ifges(4,3)) * t957 + Ifges(3,6) * t965 + t966 * t937 + Ifges(3,4) * t956 + t945 * t919 - t946 * t918 - Ifges(4,6) * t932 - Ifges(4,5) * t933 - mrSges(3,1) * t938 + mrSges(3,3) * t923 - mrSges(4,1) * t876 + mrSges(4,2) * t877 - pkin(2) * t813 - t935 * t995;
t987 = mrSges(2,1) * t961 - mrSges(2,2) * t962 + Ifges(2,3) * qJDD(1) + pkin(1) * t800 + t794 * t1003 + t791 * t1004 + t805 * t1006 + t975 * t789;
t787 = -mrSges(2,2) * g(3) - mrSges(2,3) * t961 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t983 + t791 * t981 - t794 * t978 + (-t799 * t972 - t800 * t975) * pkin(8);
t786 = mrSges(2,1) * g(3) + mrSges(2,3) * t962 + t983 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t799 - t972 * t789 + (pkin(8) * t805 + t791 * t978 + t794 * t981) * t975;
t1 = [-m(1) * g(1) + t993; -m(1) * g(2) + t1000; (-m(1) - m(2)) * g(3) + t799; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1000 - t979 * t786 + t982 * t787; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t993 + t982 * t786 + t979 * t787; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t987; t987; t789; t829; t984; t847; t985;];
tauJB  = t1;
