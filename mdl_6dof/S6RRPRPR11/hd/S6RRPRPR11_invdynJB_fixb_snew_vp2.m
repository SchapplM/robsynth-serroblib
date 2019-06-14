% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 16:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR11_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:00:45
% EndTime: 2019-05-06 16:01:00
% DurationCPUTime: 12.99s
% Computational Cost: add. (208211->368), mult. (441954->447), div. (0->0), fcn. (281407->10), ass. (0->147)
t1014 = -2 * qJD(3);
t1013 = Ifges(3,1) + Ifges(4,2);
t1006 = Ifges(3,4) + Ifges(4,6);
t1005 = Ifges(3,5) - Ifges(4,4);
t1012 = Ifges(3,2) + Ifges(4,3);
t1004 = Ifges(3,6) - Ifges(4,5);
t1011 = Ifges(3,3) + Ifges(4,1);
t972 = cos(qJ(2));
t968 = sin(qJ(2));
t996 = qJD(1) * qJD(2);
t993 = t968 * t996;
t940 = qJDD(1) * t972 - t993;
t956 = t968 * qJD(1);
t948 = pkin(3) * t956 - qJD(2) * pkin(8);
t963 = t972 ^ 2;
t975 = qJD(1) ^ 2;
t994 = t972 * t996;
t939 = qJDD(1) * t968 + t994;
t969 = sin(qJ(1));
t973 = cos(qJ(1));
t949 = t969 * g(1) - t973 * g(2);
t988 = -qJDD(1) * pkin(1) - t949;
t981 = pkin(2) * t993 + t956 * t1014 + (-t939 - t994) * qJ(3) + t988;
t862 = -t948 * t956 + (-pkin(3) * t963 - pkin(7)) * t975 + (-pkin(2) - pkin(8)) * t940 + t981;
t950 = -g(1) * t973 - g(2) * t969;
t923 = -pkin(1) * t975 + qJDD(1) * pkin(7) + t950;
t908 = -t972 * g(3) - t968 * t923;
t936 = (-pkin(2) * t972 - qJ(3) * t968) * qJD(1);
t974 = qJD(2) ^ 2;
t885 = -qJDD(2) * pkin(2) - t974 * qJ(3) + t936 * t956 + qJDD(3) - t908;
t872 = (-t968 * t972 * t975 - qJDD(2)) * pkin(8) + (t939 - t994) * pkin(3) + t885;
t967 = sin(qJ(4));
t971 = cos(qJ(4));
t851 = -t967 * t862 + t971 * t872;
t997 = qJD(1) * t972;
t934 = -qJD(2) * t967 - t971 * t997;
t900 = qJD(4) * t934 + qJDD(2) * t971 - t940 * t967;
t933 = qJDD(4) + t939;
t935 = qJD(2) * t971 - t967 * t997;
t953 = t956 + qJD(4);
t841 = (t934 * t953 - t900) * qJ(5) + (t934 * t935 + t933) * pkin(4) + t851;
t852 = t971 * t862 + t967 * t872;
t899 = -qJD(4) * t935 - qJDD(2) * t967 - t940 * t971;
t906 = pkin(4) * t953 - qJ(5) * t935;
t932 = t934 ^ 2;
t843 = -pkin(4) * t932 + qJ(5) * t899 - t906 * t953 + t852;
t964 = sin(pkin(10));
t965 = cos(pkin(10));
t903 = t934 * t964 + t935 * t965;
t835 = -0.2e1 * qJD(5) * t903 + t965 * t841 - t964 * t843;
t877 = t899 * t964 + t900 * t965;
t902 = t934 * t965 - t935 * t964;
t832 = (t902 * t953 - t877) * pkin(9) + (t902 * t903 + t933) * pkin(5) + t835;
t836 = 0.2e1 * qJD(5) * t902 + t964 * t841 + t965 * t843;
t876 = t899 * t965 - t900 * t964;
t888 = pkin(5) * t953 - pkin(9) * t903;
t901 = t902 ^ 2;
t833 = -pkin(5) * t901 + pkin(9) * t876 - t888 * t953 + t836;
t966 = sin(qJ(6));
t970 = cos(qJ(6));
t831 = t832 * t966 + t833 * t970;
t909 = -g(3) * t968 + t972 * t923;
t884 = pkin(2) * t974 - qJDD(2) * qJ(3) + qJD(2) * t1014 - t936 * t997 - t909;
t871 = -pkin(8) * t963 * t975 + pkin(3) * t940 + qJD(2) * t948 - t884;
t854 = -pkin(4) * t899 - qJ(5) * t932 + t935 * t906 + qJDD(5) + t871;
t838 = -pkin(5) * t876 - pkin(9) * t901 + t888 * t903 + t854;
t880 = t902 * t966 + t903 * t970;
t847 = -qJD(6) * t880 + t876 * t970 - t877 * t966;
t879 = t902 * t970 - t903 * t966;
t848 = qJD(6) * t879 + t876 * t966 + t877 * t970;
t951 = qJD(6) + t953;
t855 = Ifges(7,5) * t880 + Ifges(7,6) * t879 + Ifges(7,3) * t951;
t857 = Ifges(7,1) * t880 + Ifges(7,4) * t879 + Ifges(7,5) * t951;
t926 = qJDD(6) + t933;
t817 = -mrSges(7,1) * t838 + mrSges(7,3) * t831 + Ifges(7,4) * t848 + Ifges(7,2) * t847 + Ifges(7,6) * t926 - t855 * t880 + t857 * t951;
t830 = t832 * t970 - t833 * t966;
t856 = Ifges(7,4) * t880 + Ifges(7,2) * t879 + Ifges(7,6) * t951;
t818 = mrSges(7,2) * t838 - mrSges(7,3) * t830 + Ifges(7,1) * t848 + Ifges(7,4) * t847 + Ifges(7,5) * t926 + t855 * t879 - t856 * t951;
t873 = Ifges(6,5) * t903 + Ifges(6,6) * t902 + Ifges(6,3) * t953;
t875 = Ifges(6,1) * t903 + Ifges(6,4) * t902 + Ifges(6,5) * t953;
t863 = -mrSges(7,2) * t951 + mrSges(7,3) * t879;
t864 = mrSges(7,1) * t951 - mrSges(7,3) * t880;
t985 = m(7) * t838 - t847 * mrSges(7,1) + t848 * mrSges(7,2) - t879 * t863 + t880 * t864;
t859 = -mrSges(7,1) * t879 + mrSges(7,2) * t880;
t822 = m(7) * t830 + mrSges(7,1) * t926 - mrSges(7,3) * t848 - t859 * t880 + t863 * t951;
t823 = m(7) * t831 - mrSges(7,2) * t926 + mrSges(7,3) * t847 + t859 * t879 - t864 * t951;
t989 = -t822 * t966 + t970 * t823;
t802 = -mrSges(6,1) * t854 + mrSges(6,3) * t836 + Ifges(6,4) * t877 + Ifges(6,2) * t876 + Ifges(6,6) * t933 - pkin(5) * t985 + pkin(9) * t989 + t970 * t817 + t966 * t818 - t903 * t873 + t953 * t875;
t816 = t970 * t822 + t966 * t823;
t874 = Ifges(6,4) * t903 + Ifges(6,2) * t902 + Ifges(6,6) * t953;
t803 = mrSges(6,2) * t854 - mrSges(6,3) * t835 + Ifges(6,1) * t877 + Ifges(6,4) * t876 + Ifges(6,5) * t933 - pkin(9) * t816 - t817 * t966 + t818 * t970 + t873 * t902 - t874 * t953;
t886 = -mrSges(6,2) * t953 + mrSges(6,3) * t902;
t887 = mrSges(6,1) * t953 - mrSges(6,3) * t903;
t828 = m(6) * t854 - t876 * mrSges(6,1) + t877 * mrSges(6,2) - t902 * t886 + t903 * t887 + t985;
t889 = Ifges(5,5) * t935 + Ifges(5,6) * t934 + Ifges(5,3) * t953;
t891 = Ifges(5,1) * t935 + Ifges(5,4) * t934 + Ifges(5,5) * t953;
t881 = -mrSges(6,1) * t902 + mrSges(6,2) * t903;
t813 = m(6) * t835 + mrSges(6,1) * t933 - mrSges(6,3) * t877 - t881 * t903 + t886 * t953 + t816;
t814 = m(6) * t836 - mrSges(6,2) * t933 + mrSges(6,3) * t876 + t881 * t902 - t887 * t953 + t989;
t990 = -t813 * t964 + t965 * t814;
t787 = -mrSges(5,1) * t871 + mrSges(5,3) * t852 + Ifges(5,4) * t900 + Ifges(5,2) * t899 + Ifges(5,6) * t933 - pkin(4) * t828 + qJ(5) * t990 + t965 * t802 + t964 * t803 - t935 * t889 + t953 * t891;
t809 = t965 * t813 + t964 * t814;
t890 = Ifges(5,4) * t935 + Ifges(5,2) * t934 + Ifges(5,6) * t953;
t788 = mrSges(5,2) * t871 - mrSges(5,3) * t851 + Ifges(5,1) * t900 + Ifges(5,4) * t899 + Ifges(5,5) * t933 - qJ(5) * t809 - t802 * t964 + t803 * t965 + t889 * t934 - t890 * t953;
t937 = (mrSges(4,2) * t972 - mrSges(4,3) * t968) * qJD(1);
t946 = -mrSges(4,1) * t997 - qJD(2) * mrSges(4,3);
t904 = -mrSges(5,1) * t934 + mrSges(5,2) * t935;
t905 = -mrSges(5,2) * t953 + mrSges(5,3) * t934;
t806 = m(5) * t851 + mrSges(5,1) * t933 - mrSges(5,3) * t900 - t904 * t935 + t905 * t953 + t809;
t907 = mrSges(5,1) * t953 - mrSges(5,3) * t935;
t807 = m(5) * t852 - mrSges(5,2) * t933 + mrSges(5,3) * t899 + t904 * t934 - t907 * t953 + t990;
t801 = t971 * t806 + t967 * t807;
t984 = -m(4) * t885 - t939 * mrSges(4,1) - t801;
t800 = qJDD(2) * mrSges(4,2) + qJD(2) * t946 + t937 * t956 - t984;
t947 = mrSges(4,1) * t956 + qJD(2) * mrSges(4,2);
t978 = -m(5) * t871 + t899 * mrSges(5,1) - t900 * mrSges(5,2) + t934 * t905 - t935 * t907 - t828;
t977 = -m(4) * t884 + qJDD(2) * mrSges(4,3) + qJD(2) * t947 + t937 * t997 - t978;
t998 = t1005 * qJD(2) + (t1006 * t972 + t1013 * t968) * qJD(1);
t999 = t1004 * qJD(2) + (t1006 * t968 + t1012 * t972) * qJD(1);
t1010 = (t999 * t968 - t998 * t972) * qJD(1) + t1011 * qJDD(2) + t1004 * t940 + t1005 * t939 + mrSges(3,1) * t908 - mrSges(3,2) * t909 + mrSges(4,2) * t885 - mrSges(4,3) * t884 - pkin(2) * t800 - pkin(8) * t801 + qJ(3) * (t940 * mrSges(4,1) + t977) - t967 * t787 + t971 * t788;
t1007 = t975 * pkin(7);
t938 = (-mrSges(3,1) * t972 + mrSges(3,2) * t968) * qJD(1);
t945 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t997;
t798 = m(3) * t908 - t939 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t945 - t946) * qJD(2) + (-t937 - t938) * t956 + t984;
t944 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t956;
t826 = t977 - qJDD(2) * mrSges(3,2) + t938 * t997 - qJD(2) * t944 + m(3) * t909 + (mrSges(3,3) + mrSges(4,1)) * t940;
t991 = -t798 * t968 + t972 * t826;
t791 = m(2) * t950 - mrSges(2,1) * t975 - qJDD(1) * mrSges(2,2) + t991;
t922 = t988 - t1007;
t1001 = -t967 * t806 + t971 * t807;
t882 = -t940 * pkin(2) - t1007 + t981;
t987 = -m(4) * t882 - t940 * mrSges(4,2) + t947 * t956 - t1001;
t980 = -m(3) * t922 + t945 * t997 + t940 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t939 + (-t944 * t968 - t946 * t972) * qJD(1) + t987;
t795 = m(2) * t949 + qJDD(1) * mrSges(2,1) - t975 * mrSges(2,2) + t980;
t1002 = t969 * t791 + t973 * t795;
t793 = t972 * t798 + t968 * t826;
t1000 = t1011 * qJD(2) + (t1004 * t972 + t1005 * t968) * qJD(1);
t992 = t973 * t791 - t795 * t969;
t799 = -t939 * mrSges(4,3) + t946 * t997 - t987;
t784 = -mrSges(3,1) * t922 - mrSges(4,1) * t884 + mrSges(4,2) * t882 + mrSges(3,3) * t909 - pkin(2) * t799 - pkin(3) * t978 - pkin(8) * t1001 + t998 * qJD(2) + t1004 * qJDD(2) - t1000 * t956 + t1006 * t939 + t1012 * t940 - t971 * t787 - t967 * t788;
t982 = mrSges(7,1) * t830 - mrSges(7,2) * t831 + Ifges(7,5) * t848 + Ifges(7,6) * t847 + Ifges(7,3) * t926 + t880 * t856 - t879 * t857;
t976 = mrSges(5,1) * t851 + mrSges(6,1) * t835 - mrSges(5,2) * t852 - mrSges(6,2) * t836 + Ifges(5,5) * t900 + Ifges(6,5) * t877 + Ifges(5,6) * t899 + Ifges(6,6) * t876 + pkin(4) * t809 + pkin(5) * t816 + t903 * t874 - t902 * t875 + t935 * t890 - t934 * t891 + t982 + (Ifges(6,3) + Ifges(5,3)) * t933;
t786 = mrSges(4,1) * t885 + mrSges(3,2) * t922 - mrSges(3,3) * t908 - mrSges(4,3) * t882 + pkin(3) * t801 - qJ(3) * t799 - t999 * qJD(2) + t1005 * qJDD(2) + t1000 * t997 + t1006 * t940 + t1013 * t939 + t976;
t983 = mrSges(2,1) * t949 - mrSges(2,2) * t950 + Ifges(2,3) * qJDD(1) + pkin(1) * t980 + pkin(7) * t991 + t972 * t784 + t968 * t786;
t782 = mrSges(2,1) * g(3) + mrSges(2,3) * t950 + t975 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t793 - t1010;
t781 = -mrSges(2,2) * g(3) - mrSges(2,3) * t949 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t975 - pkin(7) * t793 - t784 * t968 + t786 * t972;
t1 = [-m(1) * g(1) + t992; -m(1) * g(2) + t1002; (-m(1) - m(2)) * g(3) + t793; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1002 + t973 * t781 - t969 * t782; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t992 + t969 * t781 + t973 * t782; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t983; t983; t1010; t800; t976; t828; t982;];
tauJB  = t1;
