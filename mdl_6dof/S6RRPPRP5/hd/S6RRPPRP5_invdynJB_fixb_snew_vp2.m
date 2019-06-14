% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRP5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-05-06 09:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:25:41
% EndTime: 2019-05-06 09:25:52
% DurationCPUTime: 7.54s
% Computational Cost: add. (76639->344), mult. (166881->405), div. (0->0), fcn. (99855->8), ass. (0->137)
t1003 = Ifges(6,1) + Ifges(7,1);
t985 = Ifges(6,4) - Ifges(7,5);
t999 = Ifges(7,4) + Ifges(6,5);
t1002 = Ifges(6,2) + Ifges(7,3);
t997 = Ifges(6,6) - Ifges(7,6);
t1001 = -2 * qJD(3);
t1000 = Ifges(3,1) + Ifges(4,2);
t986 = Ifges(3,4) + Ifges(4,6);
t984 = Ifges(3,5) - Ifges(4,4);
t998 = Ifges(3,2) + Ifges(4,3);
t983 = Ifges(3,6) - Ifges(4,5);
t996 = Ifges(3,3) + Ifges(4,1);
t995 = Ifges(6,3) + Ifges(7,2);
t944 = sin(pkin(9));
t945 = cos(pkin(9));
t949 = cos(qJ(2));
t974 = qJD(1) * t949;
t911 = -t944 * qJD(2) - t945 * t974;
t912 = t945 * qJD(2) - t944 * t974;
t946 = sin(qJ(5));
t989 = cos(qJ(5));
t876 = -t911 * t989 + t946 * t912;
t877 = t946 * t911 + t912 * t989;
t947 = sin(qJ(2));
t973 = t947 * qJD(1);
t934 = qJD(5) + t973;
t994 = t1002 * t876 - t985 * t877 - t997 * t934;
t993 = t1003 * t877 - t985 * t876 + t999 * t934;
t948 = sin(qJ(1));
t950 = cos(qJ(1));
t931 = -t950 * g(1) - t948 * g(2);
t952 = qJD(1) ^ 2;
t903 = -t952 * pkin(1) + qJDD(1) * pkin(7) + t931;
t881 = -t947 * g(3) + t949 * t903;
t917 = (-pkin(2) * t949 - qJ(3) * t947) * qJD(1);
t951 = qJD(2) ^ 2;
t865 = t951 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t1001 - t917 * t974 - t881;
t972 = qJD(1) * qJD(2);
t969 = t947 * t972;
t921 = t949 * qJDD(1) - t969;
t927 = pkin(3) * t973 - qJD(2) * qJ(4);
t943 = t949 ^ 2;
t849 = -t943 * t952 * qJ(4) + t921 * pkin(3) + qJD(2) * t927 + qJDD(4) - t865;
t885 = -t944 * qJDD(2) - t945 * t921;
t887 = pkin(4) * t973 - t912 * pkin(8);
t910 = t911 ^ 2;
t833 = -t885 * pkin(4) - t910 * pkin(8) + t912 * t887 + t849;
t886 = t945 * qJDD(2) - t944 * t921;
t844 = t877 * qJD(5) - t885 * t989 + t946 * t886;
t845 = -t876 * qJD(5) + t946 * t885 + t886 * t989;
t825 = (t876 * t934 - t845) * qJ(6) - 0.2e1 * qJD(6) * t877 + (t877 * t934 + t844) * pkin(5) + t833;
t868 = -t876 * mrSges(7,2) + t934 * mrSges(7,3);
t870 = -t934 * mrSges(7,1) + t877 * mrSges(7,2);
t816 = m(7) * t825 + t844 * mrSges(7,1) - t845 * mrSges(7,3) + t876 * t868 - t877 * t870;
t968 = t949 * t972;
t920 = t947 * qJDD(1) + t968;
t930 = t948 * g(1) - t950 * g(2);
t963 = -qJDD(1) * pkin(1) - t930;
t958 = pkin(2) * t969 + t973 * t1001 + (-t920 - t968) * qJ(3) + t963;
t837 = -t927 * t973 + (-pkin(3) * t943 - pkin(7)) * t952 + (-pkin(2) - qJ(4)) * t921 + t958;
t880 = -t949 * g(3) - t947 * t903;
t866 = -qJDD(2) * pkin(2) - t951 * qJ(3) + t917 * t973 + qJDD(3) - t880;
t856 = (-t947 * t949 * t952 - qJDD(2)) * qJ(4) + (t920 - t968) * pkin(3) + t866;
t830 = -0.2e1 * qJD(4) * t912 - t944 * t837 + t945 * t856;
t827 = (t911 * t973 - t886) * pkin(8) + (t911 * t912 + t920) * pkin(4) + t830;
t831 = 0.2e1 * qJD(4) * t911 + t945 * t837 + t944 * t856;
t829 = -t910 * pkin(4) + t885 * pkin(8) - t887 * t973 + t831;
t823 = t946 * t827 + t829 * t989;
t859 = t876 * pkin(5) - t877 * qJ(6);
t916 = qJDD(5) + t920;
t932 = t934 ^ 2;
t819 = -t932 * pkin(5) + t916 * qJ(6) + 0.2e1 * qJD(6) * t934 - t876 * t859 + t823;
t979 = t876 * t997 - t877 * t999 - t934 * t995;
t801 = -mrSges(6,1) * t833 - mrSges(7,1) * t825 + mrSges(7,2) * t819 + mrSges(6,3) * t823 - pkin(5) * t816 - t1002 * t844 + t985 * t845 + t979 * t877 + t997 * t916 + t993 * t934;
t822 = t827 * t989 - t946 * t829;
t820 = -t916 * pkin(5) - t932 * qJ(6) + t877 * t859 + qJDD(6) - t822;
t802 = mrSges(6,2) * t833 + mrSges(7,2) * t820 - mrSges(6,3) * t822 - mrSges(7,3) * t825 - qJ(6) * t816 + t1003 * t845 - t985 * t844 + t979 * t876 + t999 * t916 + t994 * t934;
t871 = Ifges(5,5) * t912 + Ifges(5,6) * t911 + Ifges(5,3) * t973;
t873 = Ifges(5,1) * t912 + Ifges(5,4) * t911 + Ifges(5,5) * t973;
t867 = -t934 * mrSges(6,2) - t876 * mrSges(6,3);
t869 = t934 * mrSges(6,1) - t877 * mrSges(6,3);
t957 = m(6) * t833 + t844 * mrSges(6,1) + t845 * mrSges(6,2) + t876 * t867 + t877 * t869 + t816;
t970 = m(7) * t819 + t916 * mrSges(7,3) + t934 * t870;
t860 = t876 * mrSges(7,1) - t877 * mrSges(7,3);
t978 = -t876 * mrSges(6,1) - t877 * mrSges(6,2) - t860;
t987 = -mrSges(6,3) - mrSges(7,2);
t805 = m(6) * t823 - t916 * mrSges(6,2) + t844 * t987 - t934 * t869 + t876 * t978 + t970;
t964 = -m(7) * t820 + t916 * mrSges(7,1) + t934 * t868;
t807 = m(6) * t822 + t916 * mrSges(6,1) + t845 * t987 + t934 * t867 + t877 * t978 + t964;
t965 = t805 * t989 - t946 * t807;
t781 = -mrSges(5,1) * t849 + mrSges(5,3) * t831 + Ifges(5,4) * t886 + Ifges(5,2) * t885 + Ifges(5,6) * t920 - pkin(4) * t957 + pkin(8) * t965 + t801 * t989 + t946 * t802 - t912 * t871 + t873 * t973;
t800 = t946 * t805 + t807 * t989;
t872 = Ifges(5,4) * t912 + Ifges(5,2) * t911 + Ifges(5,6) * t973;
t782 = mrSges(5,2) * t849 - mrSges(5,3) * t830 + Ifges(5,1) * t886 + Ifges(5,4) * t885 + Ifges(5,5) * t920 - pkin(8) * t800 - t946 * t801 + t802 * t989 + t911 * t871 - t872 * t973;
t918 = (mrSges(4,2) * t949 - mrSges(4,3) * t947) * qJD(1);
t928 = -mrSges(4,1) * t974 - qJD(2) * mrSges(4,3);
t878 = -t911 * mrSges(5,1) + t912 * mrSges(5,2);
t883 = -mrSges(5,2) * t973 + t911 * mrSges(5,3);
t798 = m(5) * t830 + t920 * mrSges(5,1) - t886 * mrSges(5,3) - t912 * t878 + t883 * t973 + t800;
t884 = mrSges(5,1) * t973 - t912 * mrSges(5,3);
t799 = m(5) * t831 - t920 * mrSges(5,2) + t885 * mrSges(5,3) + t911 * t878 - t884 * t973 + t965;
t795 = t945 * t798 + t944 * t799;
t960 = -m(4) * t866 - t920 * mrSges(4,1) - t795;
t794 = qJDD(2) * mrSges(4,2) + qJD(2) * t928 + t918 * t973 - t960;
t812 = m(5) * t849 - t885 * mrSges(5,1) + t886 * mrSges(5,2) - t911 * t883 + t912 * t884 + t957;
t929 = mrSges(4,1) * t973 + qJD(2) * mrSges(4,2);
t953 = -m(4) * t865 + qJDD(2) * mrSges(4,3) + qJD(2) * t929 + t918 * t974 + t812;
t975 = t984 * qJD(2) + (t1000 * t947 + t949 * t986) * qJD(1);
t976 = t983 * qJD(2) + (t947 * t986 + t949 * t998) * qJD(1);
t992 = (t947 * t976 - t949 * t975) * qJD(1) + t996 * qJDD(2) + t984 * t920 + t983 * t921 + mrSges(3,1) * t880 - mrSges(3,2) * t881 + mrSges(4,2) * t866 - mrSges(4,3) * t865 - pkin(2) * t794 + qJ(3) * (t921 * mrSges(4,1) + t953) - qJ(4) * t795 - t944 * t781 + t945 * t782;
t988 = t952 * pkin(7);
t919 = (-mrSges(3,1) * t949 + mrSges(3,2) * t947) * qJD(1);
t926 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t974;
t792 = m(3) * t880 - t920 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t926 - t928) * qJD(2) + (-t918 - t919) * t973 + t960;
t925 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t973;
t810 = (mrSges(3,3) + mrSges(4,1)) * t921 + t953 - qJD(2) * t925 + m(3) * t881 - qJDD(2) * mrSges(3,2) + t919 * t974;
t966 = -t947 * t792 + t949 * t810;
t785 = m(2) * t931 - t952 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t966;
t902 = t963 - t988;
t862 = -t921 * pkin(2) + t958 - t988;
t980 = -t944 * t798 + t945 * t799;
t962 = -m(4) * t862 - t921 * mrSges(4,2) + t929 * t973 - t980;
t956 = -m(3) * t902 + t926 * t974 + t921 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t920 + (-t925 * t947 - t928 * t949) * qJD(1) + t962;
t789 = m(2) * t930 + qJDD(1) * mrSges(2,1) - t952 * mrSges(2,2) + t956;
t981 = t948 * t785 + t950 * t789;
t787 = t949 * t792 + t947 * t810;
t977 = t996 * qJD(2) + (t947 * t984 + t949 * t983) * qJD(1);
t967 = t950 * t785 - t948 * t789;
t793 = -t920 * mrSges(4,3) + t928 * t974 - t962;
t778 = -mrSges(3,1) * t902 - mrSges(4,1) * t865 + mrSges(4,2) * t862 + mrSges(3,3) * t881 - pkin(2) * t793 + pkin(3) * t812 - qJ(4) * t980 + t975 * qJD(2) + t983 * qJDD(2) - t945 * t781 - t944 * t782 + t986 * t920 + t921 * t998 - t977 * t973;
t815 = t845 * mrSges(7,2) + t877 * t860 - t964;
t955 = mrSges(6,1) * t822 - mrSges(7,1) * t820 - mrSges(6,2) * t823 + mrSges(7,3) * t819 - pkin(5) * t815 + qJ(6) * t970 + t995 * t916 - t994 * t877 + (-qJ(6) * t860 + t993) * t876 + t999 * t845 + (-mrSges(7,2) * qJ(6) - t997) * t844;
t780 = (Ifges(5,3) + t1000) * t920 + t977 * t974 + t986 * t921 + t955 + t984 * qJDD(2) + t912 * t872 + mrSges(3,2) * t902 - t911 * t873 - mrSges(3,3) * t880 + Ifges(5,6) * t885 + Ifges(5,5) * t886 + mrSges(4,1) * t866 - mrSges(4,3) * t862 + mrSges(5,1) * t830 - mrSges(5,2) * t831 + pkin(4) * t800 + pkin(3) * t795 - qJ(3) * t793 - t976 * qJD(2);
t959 = mrSges(2,1) * t930 - mrSges(2,2) * t931 + Ifges(2,3) * qJDD(1) + pkin(1) * t956 + pkin(7) * t966 + t949 * t778 + t947 * t780;
t776 = mrSges(2,1) * g(3) + mrSges(2,3) * t931 + t952 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t787 - t992;
t775 = -mrSges(2,2) * g(3) - mrSges(2,3) * t930 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t952 - pkin(7) * t787 - t778 * t947 + t780 * t949;
t1 = [-m(1) * g(1) + t967; -m(1) * g(2) + t981; (-m(1) - m(2)) * g(3) + t787; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t981 + t950 * t775 - t948 * t776; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t967 + t948 * t775 + t950 * t776; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t959; t959; t992; t794; t812; t955; t815;];
tauJB  = t1;
