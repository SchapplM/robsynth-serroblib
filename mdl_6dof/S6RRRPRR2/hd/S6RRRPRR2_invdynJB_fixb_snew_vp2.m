% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 10:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:58:15
% EndTime: 2019-05-07 09:58:53
% DurationCPUTime: 35.90s
% Computational Cost: add. (569189->389), mult. (1279020->490), div. (0->0), fcn. (961351->12), ass. (0->154)
t971 = sin(qJ(1));
t976 = cos(qJ(1));
t956 = -g(1) * t976 - g(2) * t971;
t977 = qJD(1) ^ 2;
t944 = -pkin(1) * t977 + qJDD(1) * pkin(7) + t956;
t970 = sin(qJ(2));
t1001 = t970 * t944;
t1002 = pkin(2) * t977;
t975 = cos(qJ(2));
t996 = qJD(1) * qJD(2);
t949 = qJDD(1) * t970 + t975 * t996;
t904 = qJDD(2) * pkin(2) - t949 * pkin(8) - t1001 + (pkin(8) * t996 + t970 * t1002 - g(3)) * t975;
t931 = -g(3) * t970 + t975 * t944;
t950 = qJDD(1) * t975 - t970 * t996;
t999 = qJD(1) * t970;
t954 = qJD(2) * pkin(2) - pkin(8) * t999;
t964 = t975 ^ 2;
t905 = pkin(8) * t950 - qJD(2) * t954 - t964 * t1002 + t931;
t969 = sin(qJ(3));
t974 = cos(qJ(3));
t876 = t974 * t904 - t969 * t905;
t941 = (-t969 * t970 + t974 * t975) * qJD(1);
t914 = qJD(3) * t941 + t949 * t974 + t950 * t969;
t942 = (t969 * t975 + t970 * t974) * qJD(1);
t961 = qJDD(2) + qJDD(3);
t962 = qJD(2) + qJD(3);
t857 = (t941 * t962 - t914) * qJ(4) + (t941 * t942 + t961) * pkin(3) + t876;
t877 = t969 * t904 + t974 * t905;
t913 = -qJD(3) * t942 - t949 * t969 + t950 * t974;
t933 = pkin(3) * t962 - qJ(4) * t942;
t937 = t941 ^ 2;
t861 = -pkin(3) * t937 + qJ(4) * t913 - t933 * t962 + t877;
t965 = sin(pkin(11));
t966 = cos(pkin(11));
t928 = t941 * t965 + t942 * t966;
t841 = -0.2e1 * qJD(4) * t928 + t857 * t966 - t965 * t861;
t927 = t941 * t966 - t965 * t942;
t842 = 0.2e1 * qJD(4) * t927 + t965 * t857 + t966 * t861;
t887 = t913 * t966 - t965 * t914;
t898 = -mrSges(5,1) * t927 + mrSges(5,2) * t928;
t917 = mrSges(5,1) * t962 - mrSges(5,3) * t928;
t899 = -pkin(4) * t927 - pkin(9) * t928;
t960 = t962 ^ 2;
t839 = -pkin(4) * t960 + pkin(9) * t961 + t899 * t927 + t842;
t955 = t971 * g(1) - t976 * g(2);
t987 = -qJDD(1) * pkin(1) - t955;
t915 = -t950 * pkin(2) + t954 * t999 + (-pkin(8) * t964 - pkin(7)) * t977 + t987;
t866 = -t913 * pkin(3) - t937 * qJ(4) + t942 * t933 + qJDD(4) + t915;
t888 = t913 * t965 + t914 * t966;
t845 = (-t927 * t962 - t888) * pkin(9) + (t928 * t962 - t887) * pkin(4) + t866;
t968 = sin(qJ(5));
t973 = cos(qJ(5));
t834 = -t968 * t839 + t973 * t845;
t911 = -t928 * t968 + t962 * t973;
t864 = qJD(5) * t911 + t888 * t973 + t961 * t968;
t886 = qJDD(5) - t887;
t912 = t928 * t973 + t962 * t968;
t921 = qJD(5) - t927;
t832 = (t911 * t921 - t864) * pkin(10) + (t911 * t912 + t886) * pkin(5) + t834;
t835 = t973 * t839 + t968 * t845;
t863 = -qJD(5) * t912 - t888 * t968 + t961 * t973;
t892 = pkin(5) * t921 - pkin(10) * t912;
t907 = t911 ^ 2;
t833 = -pkin(5) * t907 + pkin(10) * t863 - t892 * t921 + t835;
t967 = sin(qJ(6));
t972 = cos(qJ(6));
t830 = t832 * t972 - t833 * t967;
t880 = t911 * t972 - t912 * t967;
t850 = qJD(6) * t880 + t863 * t967 + t864 * t972;
t881 = t911 * t967 + t912 * t972;
t858 = -mrSges(7,1) * t880 + mrSges(7,2) * t881;
t918 = qJD(6) + t921;
t867 = -mrSges(7,2) * t918 + mrSges(7,3) * t880;
t879 = qJDD(6) + t886;
t825 = m(7) * t830 + mrSges(7,1) * t879 - mrSges(7,3) * t850 - t858 * t881 + t867 * t918;
t831 = t832 * t967 + t833 * t972;
t849 = -qJD(6) * t881 + t863 * t972 - t864 * t967;
t868 = mrSges(7,1) * t918 - mrSges(7,3) * t881;
t826 = m(7) * t831 - mrSges(7,2) * t879 + mrSges(7,3) * t849 + t858 * t880 - t868 * t918;
t817 = t972 * t825 + t967 * t826;
t889 = -mrSges(6,1) * t911 + mrSges(6,2) * t912;
t890 = -mrSges(6,2) * t921 + mrSges(6,3) * t911;
t815 = m(6) * t834 + mrSges(6,1) * t886 - mrSges(6,3) * t864 - t889 * t912 + t890 * t921 + t817;
t891 = mrSges(6,1) * t921 - mrSges(6,3) * t912;
t990 = -t825 * t967 + t972 * t826;
t816 = m(6) * t835 - mrSges(6,2) * t886 + mrSges(6,3) * t863 + t889 * t911 - t891 * t921 + t990;
t991 = -t815 * t968 + t973 * t816;
t807 = m(5) * t842 - mrSges(5,2) * t961 + mrSges(5,3) * t887 + t898 * t927 - t917 * t962 + t991;
t916 = -mrSges(5,2) * t962 + mrSges(5,3) * t927;
t838 = -pkin(4) * t961 - pkin(9) * t960 + t928 * t899 - t841;
t836 = -pkin(5) * t863 - pkin(10) * t907 + t892 * t912 + t838;
t985 = m(7) * t836 - t849 * mrSges(7,1) + mrSges(7,2) * t850 - t880 * t867 + t868 * t881;
t982 = -m(6) * t838 + t863 * mrSges(6,1) - mrSges(6,2) * t864 + t911 * t890 - t891 * t912 - t985;
t821 = m(5) * t841 + mrSges(5,1) * t961 - mrSges(5,3) * t888 - t898 * t928 + t916 * t962 + t982;
t797 = t965 * t807 + t966 * t821;
t929 = -mrSges(4,1) * t941 + mrSges(4,2) * t942;
t932 = -mrSges(4,2) * t962 + mrSges(4,3) * t941;
t794 = m(4) * t876 + mrSges(4,1) * t961 - mrSges(4,3) * t914 - t929 * t942 + t932 * t962 + t797;
t934 = mrSges(4,1) * t962 - mrSges(4,3) * t942;
t992 = t966 * t807 - t821 * t965;
t795 = m(4) * t877 - mrSges(4,2) * t961 + mrSges(4,3) * t913 + t929 * t941 - t934 * t962 + t992;
t789 = t974 * t794 + t969 * t795;
t930 = -t975 * g(3) - t1001;
t939 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t970 + Ifges(3,2) * t975) * qJD(1);
t940 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t970 + Ifges(3,4) * t975) * qJD(1);
t852 = Ifges(7,5) * t881 + Ifges(7,6) * t880 + Ifges(7,3) * t918;
t854 = Ifges(7,1) * t881 + Ifges(7,4) * t880 + Ifges(7,5) * t918;
t818 = -mrSges(7,1) * t836 + mrSges(7,3) * t831 + Ifges(7,4) * t850 + Ifges(7,2) * t849 + Ifges(7,6) * t879 - t852 * t881 + t854 * t918;
t853 = Ifges(7,4) * t881 + Ifges(7,2) * t880 + Ifges(7,6) * t918;
t819 = mrSges(7,2) * t836 - mrSges(7,3) * t830 + Ifges(7,1) * t850 + Ifges(7,4) * t849 + Ifges(7,5) * t879 + t852 * t880 - t853 * t918;
t869 = Ifges(6,5) * t912 + Ifges(6,6) * t911 + Ifges(6,3) * t921;
t871 = Ifges(6,1) * t912 + Ifges(6,4) * t911 + Ifges(6,5) * t921;
t799 = -mrSges(6,1) * t838 + mrSges(6,3) * t835 + Ifges(6,4) * t864 + Ifges(6,2) * t863 + Ifges(6,6) * t886 - pkin(5) * t985 + pkin(10) * t990 + t972 * t818 + t967 * t819 - t912 * t869 + t921 * t871;
t870 = Ifges(6,4) * t912 + Ifges(6,2) * t911 + Ifges(6,6) * t921;
t801 = mrSges(6,2) * t838 - mrSges(6,3) * t834 + Ifges(6,1) * t864 + Ifges(6,4) * t863 + Ifges(6,5) * t886 - pkin(10) * t817 - t818 * t967 + t819 * t972 + t869 * t911 - t870 * t921;
t894 = Ifges(5,4) * t928 + Ifges(5,2) * t927 + Ifges(5,6) * t962;
t895 = Ifges(5,1) * t928 + Ifges(5,4) * t927 + Ifges(5,5) * t962;
t923 = Ifges(4,4) * t942 + Ifges(4,2) * t941 + Ifges(4,6) * t962;
t924 = Ifges(4,1) * t942 + Ifges(4,4) * t941 + Ifges(4,5) * t962;
t981 = -mrSges(4,1) * t876 - mrSges(5,1) * t841 + mrSges(4,2) * t877 + mrSges(5,2) * t842 - pkin(3) * t797 - pkin(4) * t982 - pkin(9) * t991 - t973 * t799 - t968 * t801 - t928 * t894 + t927 * t895 + t941 * t924 - Ifges(5,6) * t887 - Ifges(5,5) * t888 - t942 * t923 - Ifges(4,6) * t913 - Ifges(4,5) * t914 + (-Ifges(4,3) - Ifges(5,3)) * t961;
t1003 = mrSges(3,1) * t930 - mrSges(3,2) * t931 + Ifges(3,5) * t949 + Ifges(3,6) * t950 + Ifges(3,3) * qJDD(2) + pkin(2) * t789 + (t939 * t970 - t940 * t975) * qJD(1) - t981;
t948 = (-mrSges(3,1) * t975 + mrSges(3,2) * t970) * qJD(1);
t998 = qJD(1) * t975;
t953 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t998;
t787 = m(3) * t930 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t949 + qJD(2) * t953 - t948 * t999 + t789;
t952 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t999;
t993 = -t794 * t969 + t974 * t795;
t788 = m(3) * t931 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t950 - qJD(2) * t952 + t948 * t998 + t993;
t994 = -t787 * t970 + t975 * t788;
t780 = m(2) * t956 - mrSges(2,1) * t977 - qJDD(1) * mrSges(2,2) + t994;
t943 = -t977 * pkin(7) + t987;
t810 = t973 * t815 + t968 * t816;
t808 = m(5) * t866 - t887 * mrSges(5,1) + t888 * mrSges(5,2) - t927 * t916 + t928 * t917 + t810;
t983 = m(4) * t915 - t913 * mrSges(4,1) + mrSges(4,2) * t914 - t941 * t932 + t934 * t942 + t808;
t980 = -m(3) * t943 + t950 * mrSges(3,1) - mrSges(3,2) * t949 - t952 * t999 + t953 * t998 - t983;
t803 = m(2) * t955 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t977 + t980;
t1000 = t971 * t780 + t976 * t803;
t782 = t975 * t787 + t970 * t788;
t995 = t976 * t780 - t803 * t971;
t893 = Ifges(5,5) * t928 + Ifges(5,6) * t927 + Ifges(5,3) * t962;
t783 = mrSges(5,2) * t866 - mrSges(5,3) * t841 + Ifges(5,1) * t888 + Ifges(5,4) * t887 + Ifges(5,5) * t961 - pkin(9) * t810 - t799 * t968 + t801 * t973 + t893 * t927 - t894 * t962;
t984 = -mrSges(7,1) * t830 + mrSges(7,2) * t831 - Ifges(7,5) * t850 - Ifges(7,6) * t849 - Ifges(7,3) * t879 - t881 * t853 + t880 * t854;
t979 = mrSges(6,1) * t834 - mrSges(6,2) * t835 + Ifges(6,5) * t864 + Ifges(6,6) * t863 + Ifges(6,3) * t886 + pkin(5) * t817 + t912 * t870 - t911 * t871 - t984;
t790 = -mrSges(5,1) * t866 + mrSges(5,3) * t842 + Ifges(5,4) * t888 + Ifges(5,2) * t887 + Ifges(5,6) * t961 - pkin(4) * t810 - t928 * t893 + t962 * t895 - t979;
t922 = Ifges(4,5) * t942 + Ifges(4,6) * t941 + Ifges(4,3) * t962;
t776 = -mrSges(4,1) * t915 + mrSges(4,3) * t877 + Ifges(4,4) * t914 + Ifges(4,2) * t913 + Ifges(4,6) * t961 - pkin(3) * t808 + qJ(4) * t992 + t965 * t783 + t966 * t790 - t942 * t922 + t962 * t924;
t777 = mrSges(4,2) * t915 - mrSges(4,3) * t876 + Ifges(4,1) * t914 + Ifges(4,4) * t913 + Ifges(4,5) * t961 - qJ(4) * t797 + t783 * t966 - t790 * t965 + t922 * t941 - t923 * t962;
t938 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t970 + Ifges(3,6) * t975) * qJD(1);
t772 = -mrSges(3,1) * t943 + mrSges(3,3) * t931 + Ifges(3,4) * t949 + Ifges(3,2) * t950 + Ifges(3,6) * qJDD(2) - pkin(2) * t983 + pkin(8) * t993 + qJD(2) * t940 + t974 * t776 + t969 * t777 - t938 * t999;
t774 = mrSges(3,2) * t943 - mrSges(3,3) * t930 + Ifges(3,1) * t949 + Ifges(3,4) * t950 + Ifges(3,5) * qJDD(2) - pkin(8) * t789 - qJD(2) * t939 - t776 * t969 + t777 * t974 + t938 * t998;
t986 = mrSges(2,1) * t955 - mrSges(2,2) * t956 + Ifges(2,3) * qJDD(1) + pkin(1) * t980 + pkin(7) * t994 + t975 * t772 + t970 * t774;
t775 = mrSges(2,1) * g(3) + mrSges(2,3) * t956 + t977 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t782 - t1003;
t770 = -mrSges(2,2) * g(3) - mrSges(2,3) * t955 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t977 - pkin(7) * t782 - t772 * t970 + t774 * t975;
t1 = [-m(1) * g(1) + t995; -m(1) * g(2) + t1000; (-m(1) - m(2)) * g(3) + t782; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1000 + t976 * t770 - t971 * t775; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t995 + t971 * t770 + t976 * t775; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t986; t986; t1003; -t981; t808; t979; -t984;];
tauJB  = t1;
