% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-05-05 04:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:01:05
% EndTime: 2019-05-05 04:01:13
% DurationCPUTime: 5.73s
% Computational Cost: add. (59221->306), mult. (116088->359), div. (0->0), fcn. (68999->10), ass. (0->137)
t973 = Ifges(6,4) + Ifges(7,4);
t992 = Ifges(6,2) + Ifges(7,2);
t986 = Ifges(6,6) + Ifges(7,6);
t991 = -2 * qJD(4);
t990 = Ifges(4,1) + Ifges(5,2);
t989 = Ifges(6,1) + Ifges(7,1);
t974 = Ifges(4,4) + Ifges(5,6);
t972 = Ifges(4,5) - Ifges(5,4);
t988 = Ifges(6,5) + Ifges(7,5);
t987 = Ifges(4,2) + Ifges(5,3);
t971 = Ifges(4,6) - Ifges(5,5);
t985 = Ifges(4,3) + Ifges(5,1);
t984 = Ifges(6,3) + Ifges(7,3);
t925 = sin(qJ(5));
t928 = cos(qJ(5));
t929 = cos(qJ(3));
t956 = qJD(2) * t929;
t894 = -qJD(3) * t925 - t928 * t956;
t895 = qJD(3) * t928 - t925 * t956;
t926 = sin(qJ(3));
t957 = qJD(2) * t926;
t912 = qJD(5) + t957;
t983 = t992 * t894 + t973 * t895 + t986 * t912;
t955 = qJD(2) * qJD(3);
t949 = t926 * t955;
t900 = qJDD(2) * t929 - t949;
t860 = -qJD(5) * t895 - qJDD(3) * t925 - t900 * t928;
t866 = -mrSges(7,2) * t912 + mrSges(7,3) * t894;
t921 = sin(pkin(10));
t923 = cos(pkin(10));
t903 = g(1) * t921 - g(2) * t923;
t904 = -g(1) * t923 - g(2) * t921;
t920 = -g(3) + qJDD(1);
t930 = cos(qJ(2));
t924 = cos(pkin(6));
t927 = sin(qJ(2));
t967 = t924 * t927;
t922 = sin(pkin(6));
t968 = t922 * t927;
t844 = t903 * t967 + t930 * t904 + t920 * t968;
t932 = qJD(2) ^ 2;
t841 = -pkin(2) * t932 + qJDD(2) * pkin(8) + t844;
t871 = -t903 * t922 + t920 * t924;
t835 = t929 * t841 + t926 * t871;
t896 = (-pkin(3) * t929 - qJ(4) * t926) * qJD(2);
t931 = qJD(3) ^ 2;
t831 = pkin(3) * t931 - qJDD(3) * qJ(4) + qJD(3) * t991 - t896 * t956 - t835;
t910 = pkin(4) * t957 - qJD(3) * pkin(9);
t919 = t929 ^ 2;
t977 = pkin(9) * t932;
t827 = pkin(4) * t900 + qJD(3) * t910 - t919 * t977 - t831;
t868 = pkin(5) * t912 - qJ(6) * t895;
t890 = t894 ^ 2;
t823 = -t860 * pkin(5) - qJ(6) * t890 + t868 * t895 + qJDD(6) + t827;
t861 = qJD(5) * t894 + qJDD(3) * t928 - t900 * t925;
t869 = mrSges(7,1) * t912 - mrSges(7,3) * t895;
t951 = m(7) * t823 + t861 * mrSges(7,2) + t895 * t869;
t816 = -t860 * mrSges(7,1) - t894 * t866 + t951;
t950 = t929 * t955;
t899 = qJDD(2) * t926 + t950;
t838 = t926 * t841;
t944 = -t931 * qJ(4) + t896 * t957 + qJDD(4) + t838;
t978 = -pkin(3) - pkin(9);
t828 = t899 * pkin(4) + t978 * qJDD(3) + (-pkin(4) * t955 - t926 * t977 - t871) * t929 + t944;
t843 = -t927 * t904 + (t903 * t924 + t920 * t922) * t930;
t938 = -qJDD(2) * pkin(2) - t843;
t936 = pkin(3) * t949 + t957 * t991 + (-t899 - t950) * qJ(4) + t938;
t830 = -t910 * t957 + (-pkin(4) * t919 - pkin(8)) * t932 + t978 * t900 + t936;
t821 = t925 * t828 + t928 * t830;
t818 = -pkin(5) * t890 + t860 * qJ(6) + 0.2e1 * qJD(6) * t894 - t868 * t912 + t821;
t891 = qJDD(5) + t899;
t863 = -mrSges(7,1) * t894 + mrSges(7,2) * t895;
t952 = m(7) * t818 + t860 * mrSges(7,3) + t894 * t863;
t962 = -t973 * t894 - t989 * t895 - t988 * t912;
t963 = -t986 * t894 - t988 * t895 - t984 * t912;
t793 = -mrSges(6,1) * t827 + mrSges(6,3) * t821 - mrSges(7,1) * t823 + mrSges(7,3) * t818 - pkin(5) * t816 + qJ(6) * t952 + (-qJ(6) * t869 - t962) * t912 + t963 * t895 + (-mrSges(7,2) * qJ(6) + t986) * t891 + t973 * t861 + t992 * t860;
t897 = (mrSges(5,2) * t929 - mrSges(5,3) * t926) * qJD(2);
t907 = -mrSges(5,1) * t956 - qJD(3) * mrSges(5,3);
t820 = t928 * t828 - t925 * t830;
t864 = -mrSges(6,1) * t894 + mrSges(6,2) * t895;
t867 = -mrSges(6,2) * t912 + mrSges(6,3) * t894;
t815 = -0.2e1 * qJD(6) * t895 + (t894 * t912 - t861) * qJ(6) + (t894 * t895 + t891) * pkin(5) + t820;
t953 = m(7) * t815 + t891 * mrSges(7,1) + t912 * t866;
t807 = m(6) * t820 + t891 * mrSges(6,1) + t912 * t867 + (-t863 - t864) * t895 + (-mrSges(6,3) - mrSges(7,3)) * t861 + t953;
t870 = mrSges(6,1) * t912 - mrSges(6,3) * t895;
t809 = m(6) * t821 + t860 * mrSges(6,3) + t894 * t864 + (-t869 - t870) * t912 + (-mrSges(6,2) - mrSges(7,2)) * t891 + t952;
t802 = t928 * t807 + t925 * t809;
t966 = t929 * t871;
t832 = -qJDD(3) * pkin(3) + t944 - t966;
t939 = -m(5) * t832 - t899 * mrSges(5,1) - t802;
t800 = qJDD(3) * mrSges(5,2) + qJD(3) * t907 + t897 * t957 - t939;
t812 = -t861 * mrSges(7,3) - t895 * t863 + t953;
t801 = mrSges(6,2) * t827 + mrSges(7,2) * t823 - mrSges(6,3) * t820 - mrSges(7,3) * t815 - qJ(6) * t812 + t973 * t860 + t989 * t861 + t988 * t891 - t963 * t894 - t983 * t912;
t834 = -t838 + t966;
t908 = mrSges(5,1) * t957 + qJD(3) * mrSges(5,2);
t980 = -m(6) * t827 - t861 * mrSges(6,2) + (mrSges(6,1) + mrSges(7,1)) * t860 - t895 * t870 + (t866 + t867) * t894 - t951;
t937 = -m(5) * t831 + qJDD(3) * mrSges(5,3) + qJD(3) * t908 + t897 * t956 - t980;
t958 = t972 * qJD(3) + (t990 * t926 + t974 * t929) * qJD(2);
t959 = t971 * qJD(3) + (t974 * t926 + t987 * t929) * qJD(2);
t982 = (t959 * t926 - t958 * t929) * qJD(2) + t985 * qJDD(3) + t972 * t899 + t971 * t900 + mrSges(4,1) * t834 - mrSges(4,2) * t835 + mrSges(5,2) * t832 - mrSges(5,3) * t831 - pkin(3) * t800 - pkin(9) * t802 + qJ(4) * (t900 * mrSges(5,1) + t937) - t925 * t793 + t928 * t801;
t976 = t932 * pkin(8);
t840 = t938 - t976;
t905 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t957;
t906 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t956;
t833 = -t900 * pkin(3) + t936 - t976;
t964 = -t925 * t807 + t928 * t809;
t942 = -m(5) * t833 - t900 * mrSges(5,2) + t908 * t957 - t964;
t934 = -m(4) * t840 + t906 * t956 + t900 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t899 + (-t905 * t926 - t907 * t929) * qJD(2) + t942;
t796 = m(3) * t843 + qJDD(2) * mrSges(3,1) - t932 * mrSges(3,2) + t934;
t969 = t796 * t930;
t898 = (-mrSges(4,1) * t929 + mrSges(4,2) * t926) * qJD(2);
t798 = m(4) * t834 - t899 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t906 - t907) * qJD(3) + (-t897 - t898) * t957 + t939;
t805 = t898 * t956 + t937 - qJDD(3) * mrSges(4,2) + (mrSges(4,3) + mrSges(5,1)) * t900 + m(4) * t835 - qJD(3) * t905;
t947 = -t798 * t926 + t929 * t805;
t789 = m(3) * t844 - mrSges(3,1) * t932 - qJDD(2) * mrSges(3,2) + t947;
t792 = t929 * t798 + t926 * t805;
t791 = m(3) * t871 + t792;
t780 = t789 * t967 - t791 * t922 + t924 * t969;
t778 = m(2) * t903 + t780;
t785 = t930 * t789 - t796 * t927;
t784 = m(2) * t904 + t785;
t965 = t923 * t778 + t921 * t784;
t960 = t985 * qJD(3) + (t972 * t926 + t971 * t929) * qJD(2);
t779 = t789 * t968 + t924 * t791 + t922 * t969;
t948 = -t778 * t921 + t923 * t784;
t946 = m(2) * t920 + t779;
t799 = -t899 * mrSges(5,3) + t907 * t956 - t942;
t776 = -mrSges(4,1) * t840 - mrSges(5,1) * t831 + mrSges(5,2) * t833 + mrSges(4,3) * t835 - pkin(3) * t799 - pkin(4) * t980 - pkin(9) * t964 + t958 * qJD(3) + t971 * qJDD(3) - t928 * t793 - t925 * t801 + t974 * t899 + t987 * t900 - t960 * t957;
t935 = mrSges(6,1) * t820 + mrSges(7,1) * t815 - mrSges(6,2) * t821 - mrSges(7,2) * t818 + pkin(5) * t812 + t986 * t860 + t988 * t861 + t984 * t891 + t962 * t894 + t983 * t895;
t781 = mrSges(5,1) * t832 + mrSges(4,2) * t840 - mrSges(4,3) * t834 - mrSges(5,3) * t833 + pkin(4) * t802 - qJ(4) * t799 - t959 * qJD(3) + t972 * qJDD(3) + t990 * t899 + t974 * t900 + t960 * t956 + t935;
t774 = mrSges(3,2) * t871 - mrSges(3,3) * t843 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t932 - pkin(8) * t792 - t776 * t926 + t781 * t929;
t775 = -mrSges(3,1) * t871 + mrSges(3,3) * t844 + t932 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t792 - t982;
t940 = pkin(7) * t785 + t774 * t927 + t775 * t930;
t773 = mrSges(3,1) * t843 - mrSges(3,2) * t844 + Ifges(3,3) * qJDD(2) + pkin(2) * t934 + pkin(8) * t947 + t929 * t776 + t926 * t781;
t772 = mrSges(2,2) * t920 - mrSges(2,3) * t903 + t930 * t774 - t927 * t775 + (-t779 * t922 - t780 * t924) * pkin(7);
t771 = -mrSges(2,1) * t920 + mrSges(2,3) * t904 - pkin(1) * t779 - t922 * t773 + t940 * t924;
t1 = [-m(1) * g(1) + t948; -m(1) * g(2) + t965; -m(1) * g(3) + t946; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t965 - t921 * t771 + t923 * t772; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t948 + t923 * t771 + t921 * t772; -mrSges(1,1) * g(2) + mrSges(2,1) * t903 + mrSges(1,2) * g(1) - mrSges(2,2) * t904 + pkin(1) * t780 + t924 * t773 + t940 * t922; t946; t773; t982; t800; t935; t816;];
tauJB  = t1;
