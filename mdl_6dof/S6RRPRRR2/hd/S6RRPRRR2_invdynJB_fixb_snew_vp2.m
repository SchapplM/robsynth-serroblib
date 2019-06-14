% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 19:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:54:35
% EndTime: 2019-05-06 19:55:21
% DurationCPUTime: 36.16s
% Computational Cost: add. (540069->389), mult. (1247917->490), div. (0->0), fcn. (944500->12), ass. (0->154)
t954 = sin(qJ(2));
t959 = cos(qJ(2));
t979 = qJD(1) * qJD(2);
t934 = qJDD(1) * t954 + t959 * t979;
t955 = sin(qJ(1));
t960 = cos(qJ(1));
t941 = -g(1) * t960 - g(2) * t955;
t961 = qJD(1) ^ 2;
t929 = -pkin(1) * t961 + qJDD(1) * pkin(7) + t941;
t983 = t954 * t929;
t985 = pkin(2) * t961;
t890 = qJDD(2) * pkin(2) - t934 * qJ(3) - t983 + (qJ(3) * t979 + t954 * t985 - g(3)) * t959;
t914 = -g(3) * t954 + t959 * t929;
t935 = qJDD(1) * t959 - t954 * t979;
t981 = qJD(1) * t954;
t937 = qJD(2) * pkin(2) - qJ(3) * t981;
t948 = t959 ^ 2;
t891 = qJ(3) * t935 - qJD(2) * t937 - t948 * t985 + t914;
t949 = sin(pkin(11));
t950 = cos(pkin(11));
t924 = (t949 * t959 + t950 * t954) * qJD(1);
t861 = -0.2e1 * qJD(3) * t924 + t950 * t890 - t949 * t891;
t912 = t934 * t950 + t935 * t949;
t923 = (-t949 * t954 + t950 * t959) * qJD(1);
t845 = (qJD(2) * t923 - t912) * pkin(8) + (t923 * t924 + qJDD(2)) * pkin(3) + t861;
t862 = 0.2e1 * qJD(3) * t923 + t949 * t890 + t950 * t891;
t911 = -t934 * t949 + t935 * t950;
t917 = qJD(2) * pkin(3) - pkin(8) * t924;
t922 = t923 ^ 2;
t849 = -pkin(3) * t922 + pkin(8) * t911 - qJD(2) * t917 + t862;
t953 = sin(qJ(4));
t958 = cos(qJ(4));
t838 = t953 * t845 + t958 * t849;
t905 = t923 * t953 + t924 * t958;
t872 = -t905 * qJD(4) + t911 * t958 - t953 * t912;
t904 = t923 * t958 - t953 * t924;
t885 = -mrSges(5,1) * t904 + mrSges(5,2) * t905;
t946 = qJD(2) + qJD(4);
t897 = mrSges(5,1) * t946 - mrSges(5,3) * t905;
t945 = qJDD(2) + qJDD(4);
t886 = -pkin(4) * t904 - pkin(9) * t905;
t944 = t946 ^ 2;
t827 = -pkin(4) * t944 + pkin(9) * t945 + t886 * t904 + t838;
t940 = t955 * g(1) - t960 * g(2);
t971 = -qJDD(1) * pkin(1) - t940;
t893 = -t935 * pkin(2) + qJDD(3) + t937 * t981 + (-qJ(3) * t948 - pkin(7)) * t961 + t971;
t860 = -t911 * pkin(3) - t922 * pkin(8) + t924 * t917 + t893;
t873 = qJD(4) * t904 + t911 * t953 + t912 * t958;
t835 = (-t904 * t946 - t873) * pkin(9) + (t905 * t946 - t872) * pkin(4) + t860;
t952 = sin(qJ(5));
t957 = cos(qJ(5));
t822 = -t952 * t827 + t957 * t835;
t894 = -t905 * t952 + t946 * t957;
t852 = qJD(5) * t894 + t873 * t957 + t945 * t952;
t871 = qJDD(5) - t872;
t895 = t905 * t957 + t946 * t952;
t900 = qJD(5) - t904;
t820 = (t894 * t900 - t852) * pkin(10) + (t894 * t895 + t871) * pkin(5) + t822;
t823 = t957 * t827 + t952 * t835;
t851 = -qJD(5) * t895 - t873 * t952 + t945 * t957;
t879 = pkin(5) * t900 - pkin(10) * t895;
t892 = t894 ^ 2;
t821 = -pkin(5) * t892 + pkin(10) * t851 - t879 * t900 + t823;
t951 = sin(qJ(6));
t956 = cos(qJ(6));
t818 = t820 * t956 - t821 * t951;
t874 = t894 * t956 - t895 * t951;
t832 = qJD(6) * t874 + t851 * t951 + t852 * t956;
t875 = t894 * t951 + t895 * t956;
t846 = -mrSges(7,1) * t874 + mrSges(7,2) * t875;
t898 = qJD(6) + t900;
t853 = -mrSges(7,2) * t898 + mrSges(7,3) * t874;
t866 = qJDD(6) + t871;
t813 = m(7) * t818 + mrSges(7,1) * t866 - mrSges(7,3) * t832 - t846 * t875 + t853 * t898;
t819 = t820 * t951 + t821 * t956;
t831 = -qJD(6) * t875 + t851 * t956 - t852 * t951;
t854 = mrSges(7,1) * t898 - mrSges(7,3) * t875;
t814 = m(7) * t819 - mrSges(7,2) * t866 + mrSges(7,3) * t831 + t846 * t874 - t854 * t898;
t805 = t956 * t813 + t951 * t814;
t876 = -mrSges(6,1) * t894 + mrSges(6,2) * t895;
t877 = -mrSges(6,2) * t900 + mrSges(6,3) * t894;
t803 = m(6) * t822 + mrSges(6,1) * t871 - mrSges(6,3) * t852 - t876 * t895 + t877 * t900 + t805;
t878 = mrSges(6,1) * t900 - mrSges(6,3) * t895;
t973 = -t813 * t951 + t956 * t814;
t804 = m(6) * t823 - mrSges(6,2) * t871 + mrSges(6,3) * t851 + t876 * t894 - t878 * t900 + t973;
t974 = -t803 * t952 + t957 * t804;
t795 = m(5) * t838 - mrSges(5,2) * t945 + mrSges(5,3) * t872 + t885 * t904 - t897 * t946 + t974;
t837 = t845 * t958 - t953 * t849;
t896 = -mrSges(5,2) * t946 + mrSges(5,3) * t904;
t826 = -pkin(4) * t945 - pkin(9) * t944 + t905 * t886 - t837;
t824 = -pkin(5) * t851 - pkin(10) * t892 + t879 * t895 + t826;
t968 = m(7) * t824 - t831 * mrSges(7,1) + mrSges(7,2) * t832 - t874 * t853 + t854 * t875;
t965 = -m(6) * t826 + t851 * mrSges(6,1) - mrSges(6,2) * t852 + t894 * t877 - t878 * t895 - t968;
t809 = m(5) * t837 + mrSges(5,1) * t945 - mrSges(5,3) * t873 - t885 * t905 + t896 * t946 + t965;
t785 = t953 * t795 + t958 * t809;
t908 = -mrSges(4,1) * t923 + mrSges(4,2) * t924;
t915 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t923;
t783 = m(4) * t861 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t912 + qJD(2) * t915 - t908 * t924 + t785;
t916 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t924;
t975 = t958 * t795 - t809 * t953;
t784 = m(4) * t862 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t911 - qJD(2) * t916 + t908 * t923 + t975;
t778 = t950 * t783 + t949 * t784;
t902 = Ifges(4,4) * t924 + Ifges(4,2) * t923 + Ifges(4,6) * qJD(2);
t903 = Ifges(4,1) * t924 + Ifges(4,4) * t923 + Ifges(4,5) * qJD(2);
t913 = -t959 * g(3) - t983;
t926 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t954 + Ifges(3,2) * t959) * qJD(1);
t927 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t954 + Ifges(3,4) * t959) * qJD(1);
t840 = Ifges(7,5) * t875 + Ifges(7,6) * t874 + Ifges(7,3) * t898;
t842 = Ifges(7,1) * t875 + Ifges(7,4) * t874 + Ifges(7,5) * t898;
t806 = -mrSges(7,1) * t824 + mrSges(7,3) * t819 + Ifges(7,4) * t832 + Ifges(7,2) * t831 + Ifges(7,6) * t866 - t840 * t875 + t842 * t898;
t841 = Ifges(7,4) * t875 + Ifges(7,2) * t874 + Ifges(7,6) * t898;
t807 = mrSges(7,2) * t824 - mrSges(7,3) * t818 + Ifges(7,1) * t832 + Ifges(7,4) * t831 + Ifges(7,5) * t866 + t840 * t874 - t841 * t898;
t855 = Ifges(6,5) * t895 + Ifges(6,6) * t894 + Ifges(6,3) * t900;
t857 = Ifges(6,1) * t895 + Ifges(6,4) * t894 + Ifges(6,5) * t900;
t787 = -mrSges(6,1) * t826 + mrSges(6,3) * t823 + Ifges(6,4) * t852 + Ifges(6,2) * t851 + Ifges(6,6) * t871 - pkin(5) * t968 + pkin(10) * t973 + t956 * t806 + t951 * t807 - t895 * t855 + t900 * t857;
t856 = Ifges(6,4) * t895 + Ifges(6,2) * t894 + Ifges(6,6) * t900;
t789 = mrSges(6,2) * t826 - mrSges(6,3) * t822 + Ifges(6,1) * t852 + Ifges(6,4) * t851 + Ifges(6,5) * t871 - pkin(10) * t805 - t806 * t951 + t807 * t956 + t855 * t894 - t856 * t900;
t881 = Ifges(5,4) * t905 + Ifges(5,2) * t904 + Ifges(5,6) * t946;
t882 = Ifges(5,1) * t905 + Ifges(5,4) * t904 + Ifges(5,5) * t946;
t966 = -mrSges(5,1) * t837 + mrSges(5,2) * t838 - Ifges(5,5) * t873 - Ifges(5,6) * t872 - Ifges(5,3) * t945 - pkin(4) * t965 - pkin(9) * t974 - t957 * t787 - t952 * t789 - t905 * t881 + t904 * t882;
t986 = mrSges(3,1) * t913 + mrSges(4,1) * t861 - mrSges(3,2) * t914 - mrSges(4,2) * t862 + Ifges(3,5) * t934 + Ifges(4,5) * t912 + Ifges(3,6) * t935 + Ifges(4,6) * t911 + pkin(2) * t778 + pkin(3) * t785 + (t926 * t954 - t927 * t959) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t924 * t902 - t923 * t903 - t966;
t933 = (-mrSges(3,1) * t959 + mrSges(3,2) * t954) * qJD(1);
t980 = qJD(1) * t959;
t939 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t980;
t776 = m(3) * t913 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t934 + qJD(2) * t939 - t933 * t981 + t778;
t938 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t981;
t976 = -t783 * t949 + t950 * t784;
t777 = m(3) * t914 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t935 - qJD(2) * t938 + t933 * t980 + t976;
t977 = -t776 * t954 + t959 * t777;
t769 = m(2) * t941 - mrSges(2,1) * t961 - qJDD(1) * mrSges(2,2) + t977;
t798 = t957 * t803 + t952 * t804;
t970 = m(5) * t860 - t872 * mrSges(5,1) + t873 * mrSges(5,2) - t904 * t896 + t905 * t897 + t798;
t796 = m(4) * t893 - t911 * mrSges(4,1) + mrSges(4,2) * t912 - t923 * t915 + t916 * t924 + t970;
t928 = -t961 * pkin(7) + t971;
t964 = -m(3) * t928 + t935 * mrSges(3,1) - mrSges(3,2) * t934 - t938 * t981 + t939 * t980 - t796;
t791 = m(2) * t940 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t961 + t964;
t982 = t955 * t769 + t960 * t791;
t771 = t959 * t776 + t954 * t777;
t978 = t960 * t769 - t791 * t955;
t880 = Ifges(5,5) * t905 + Ifges(5,6) * t904 + Ifges(5,3) * t946;
t772 = mrSges(5,2) * t860 - mrSges(5,3) * t837 + Ifges(5,1) * t873 + Ifges(5,4) * t872 + Ifges(5,5) * t945 - pkin(9) * t798 - t787 * t952 + t789 * t957 + t880 * t904 - t881 * t946;
t967 = -mrSges(7,1) * t818 + mrSges(7,2) * t819 - Ifges(7,5) * t832 - Ifges(7,6) * t831 - Ifges(7,3) * t866 - t875 * t841 + t874 * t842;
t963 = mrSges(6,1) * t822 - mrSges(6,2) * t823 + Ifges(6,5) * t852 + Ifges(6,6) * t851 + Ifges(6,3) * t871 + pkin(5) * t805 + t895 * t856 - t894 * t857 - t967;
t779 = -mrSges(5,1) * t860 + mrSges(5,3) * t838 + Ifges(5,4) * t873 + Ifges(5,2) * t872 + Ifges(5,6) * t945 - pkin(4) * t798 - t905 * t880 + t946 * t882 - t963;
t901 = Ifges(4,5) * t924 + Ifges(4,6) * t923 + Ifges(4,3) * qJD(2);
t765 = -mrSges(4,1) * t893 + mrSges(4,3) * t862 + Ifges(4,4) * t912 + Ifges(4,2) * t911 + Ifges(4,6) * qJDD(2) - pkin(3) * t970 + pkin(8) * t975 + qJD(2) * t903 + t953 * t772 + t958 * t779 - t924 * t901;
t766 = mrSges(4,2) * t893 - mrSges(4,3) * t861 + Ifges(4,1) * t912 + Ifges(4,4) * t911 + Ifges(4,5) * qJDD(2) - pkin(8) * t785 - qJD(2) * t902 + t772 * t958 - t779 * t953 + t901 * t923;
t925 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t954 + Ifges(3,6) * t959) * qJD(1);
t761 = -mrSges(3,1) * t928 + mrSges(3,3) * t914 + Ifges(3,4) * t934 + Ifges(3,2) * t935 + Ifges(3,6) * qJDD(2) - pkin(2) * t796 + qJ(3) * t976 + qJD(2) * t927 + t950 * t765 + t949 * t766 - t925 * t981;
t763 = mrSges(3,2) * t928 - mrSges(3,3) * t913 + Ifges(3,1) * t934 + Ifges(3,4) * t935 + Ifges(3,5) * qJDD(2) - qJ(3) * t778 - qJD(2) * t926 - t765 * t949 + t766 * t950 + t925 * t980;
t969 = mrSges(2,1) * t940 - mrSges(2,2) * t941 + Ifges(2,3) * qJDD(1) + pkin(1) * t964 + pkin(7) * t977 + t959 * t761 + t954 * t763;
t764 = mrSges(2,1) * g(3) + mrSges(2,3) * t941 + t961 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t771 - t986;
t759 = -mrSges(2,2) * g(3) - mrSges(2,3) * t940 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t961 - pkin(7) * t771 - t761 * t954 + t763 * t959;
t1 = [-m(1) * g(1) + t978; -m(1) * g(2) + t982; (-m(1) - m(2)) * g(3) + t771; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t982 + t960 * t759 - t955 * t764; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t978 + t955 * t759 + t960 * t764; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t969; t969; t986; t796; -t966; t963; -t967;];
tauJB  = t1;
