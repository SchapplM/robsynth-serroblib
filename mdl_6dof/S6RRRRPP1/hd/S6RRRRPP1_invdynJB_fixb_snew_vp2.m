% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:00
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:55:08
% EndTime: 2019-05-07 17:55:24
% DurationCPUTime: 16.20s
% Computational Cost: add. (270461->365), mult. (546802->448), div. (0->0), fcn. (393607->10), ass. (0->146)
t976 = Ifges(6,1) + Ifges(7,1);
t968 = Ifges(6,4) - Ifges(7,5);
t967 = Ifges(6,5) + Ifges(7,4);
t975 = -Ifges(6,2) - Ifges(7,3);
t974 = -Ifges(7,2) - Ifges(6,3);
t966 = Ifges(6,6) - Ifges(7,6);
t929 = sin(qJ(3));
t930 = sin(qJ(2));
t933 = cos(qJ(3));
t934 = cos(qJ(2));
t904 = (t929 * t930 - t933 * t934) * qJD(1);
t956 = qJD(1) * qJD(2);
t912 = t930 * qJDD(1) + t934 * t956;
t931 = sin(qJ(1));
t935 = cos(qJ(1));
t919 = -t935 * g(1) - t931 * g(2);
t936 = qJD(1) ^ 2;
t907 = -t936 * pkin(1) + qJDD(1) * pkin(7) + t919;
t964 = t930 * t907;
t970 = pkin(2) * t936;
t869 = qJDD(2) * pkin(2) - t912 * pkin(8) - t964 + (pkin(8) * t956 + t930 * t970 - g(3)) * t934;
t895 = -t930 * g(3) + t934 * t907;
t913 = t934 * qJDD(1) - t930 * t956;
t958 = qJD(1) * t930;
t917 = qJD(2) * pkin(2) - pkin(8) * t958;
t926 = t934 ^ 2;
t870 = t913 * pkin(8) - qJD(2) * t917 - t926 * t970 + t895;
t845 = t929 * t869 + t933 * t870;
t905 = (t929 * t934 + t930 * t933) * qJD(1);
t878 = -t905 * qJD(3) - t929 * t912 + t933 * t913;
t888 = t904 * mrSges(4,1) + t905 * mrSges(4,2);
t924 = qJD(2) + qJD(3);
t897 = t924 * mrSges(4,1) - t905 * mrSges(4,3);
t923 = qJDD(2) + qJDD(3);
t879 = -t904 * qJD(3) + t933 * t912 + t929 * t913;
t918 = t931 * g(1) - t935 * g(2);
t945 = -qJDD(1) * pkin(1) - t918;
t880 = -t913 * pkin(2) + t917 * t958 + (-pkin(8) * t926 - pkin(7)) * t936 + t945;
t822 = (t904 * t924 - t879) * pkin(9) + (t905 * t924 - t878) * pkin(3) + t880;
t889 = t904 * pkin(3) - t905 * pkin(9);
t922 = t924 ^ 2;
t830 = -t922 * pkin(3) + t923 * pkin(9) - t904 * t889 + t845;
t928 = sin(qJ(4));
t932 = cos(qJ(4));
t817 = t932 * t822 - t928 * t830;
t892 = -t928 * t905 + t932 * t924;
t851 = t892 * qJD(4) + t932 * t879 + t928 * t923;
t877 = qJDD(4) - t878;
t893 = t932 * t905 + t928 * t924;
t900 = qJD(4) + t904;
t814 = (t892 * t900 - t851) * qJ(5) + (t892 * t893 + t877) * pkin(4) + t817;
t818 = t928 * t822 + t932 * t830;
t850 = -t893 * qJD(4) - t928 * t879 + t932 * t923;
t882 = t900 * pkin(4) - t893 * qJ(5);
t891 = t892 ^ 2;
t816 = -t891 * pkin(4) + t850 * qJ(5) - t900 * t882 + t818;
t927 = sin(pkin(10));
t965 = cos(pkin(10));
t862 = -t965 * t892 + t927 * t893;
t971 = -2 * qJD(5);
t810 = t927 * t814 + t965 * t816 + t862 * t971;
t826 = -t965 * t850 + t927 * t851;
t863 = t927 * t892 + t965 * t893;
t854 = t900 * mrSges(6,1) - t863 * mrSges(6,3);
t840 = t862 * pkin(5) - t863 * qJ(6);
t899 = t900 ^ 2;
t807 = -t899 * pkin(5) + t877 * qJ(6) + 0.2e1 * qJD(6) * t900 - t862 * t840 + t810;
t855 = -t900 * mrSges(7,1) + t863 * mrSges(7,2);
t954 = m(7) * t807 + t877 * mrSges(7,3) + t900 * t855;
t841 = t862 * mrSges(7,1) - t863 * mrSges(7,3);
t959 = -t862 * mrSges(6,1) - t863 * mrSges(6,2) - t841;
t969 = -mrSges(6,3) - mrSges(7,2);
t795 = m(6) * t810 - t877 * mrSges(6,2) + t969 * t826 - t900 * t854 + t959 * t862 + t954;
t944 = t965 * t814 - t927 * t816;
t809 = t863 * t971 + t944;
t827 = t927 * t850 + t965 * t851;
t853 = -t900 * mrSges(6,2) - t862 * mrSges(6,3);
t808 = -t877 * pkin(5) - t899 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t840) * t863 - t944;
t852 = -t862 * mrSges(7,2) + t900 * mrSges(7,3);
t948 = -m(7) * t808 + t877 * mrSges(7,1) + t900 * t852;
t797 = m(6) * t809 + t877 * mrSges(6,1) + t969 * t827 + t900 * t853 + t959 * t863 + t948;
t790 = t927 * t795 + t965 * t797;
t867 = -t892 * mrSges(5,1) + t893 * mrSges(5,2);
t881 = -t900 * mrSges(5,2) + t892 * mrSges(5,3);
t788 = m(5) * t817 + t877 * mrSges(5,1) - t851 * mrSges(5,3) - t893 * t867 + t900 * t881 + t790;
t883 = t900 * mrSges(5,1) - t893 * mrSges(5,3);
t949 = t965 * t795 - t927 * t797;
t789 = m(5) * t818 - t877 * mrSges(5,2) + t850 * mrSges(5,3) + t892 * t867 - t900 * t883 + t949;
t950 = -t928 * t788 + t932 * t789;
t781 = m(4) * t845 - t923 * mrSges(4,2) + t878 * mrSges(4,3) - t904 * t888 - t924 * t897 + t950;
t844 = t933 * t869 - t929 * t870;
t896 = -t924 * mrSges(4,2) - t904 * mrSges(4,3);
t829 = -t923 * pkin(3) - t922 * pkin(9) + t905 * t889 - t844;
t819 = -t850 * pkin(4) - t891 * qJ(5) + t893 * t882 + qJDD(5) + t829;
t812 = -0.2e1 * qJD(6) * t863 + (t862 * t900 - t827) * qJ(6) + (t863 * t900 + t826) * pkin(5) + t819;
t805 = m(7) * t812 + t826 * mrSges(7,1) - t827 * mrSges(7,3) + t862 * t852 - t863 * t855;
t802 = m(6) * t819 + t826 * mrSges(6,1) + t827 * mrSges(6,2) + t862 * t853 + t863 * t854 + t805;
t939 = -m(5) * t829 + t850 * mrSges(5,1) - t851 * mrSges(5,2) + t892 * t881 - t893 * t883 - t802;
t799 = m(4) * t844 + t923 * mrSges(4,1) - t879 * mrSges(4,3) - t905 * t888 + t924 * t896 + t939;
t775 = t929 * t781 + t933 * t799;
t894 = -t934 * g(3) - t964;
t902 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t930 + Ifges(3,2) * t934) * qJD(1);
t903 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t930 + Ifges(3,4) * t934) * qJD(1);
t960 = -t968 * t862 + t976 * t863 + t967 * t900;
t962 = t966 * t862 - t967 * t863 + t974 * t900;
t791 = -mrSges(6,1) * t819 - mrSges(7,1) * t812 + mrSges(7,2) * t807 + mrSges(6,3) * t810 - pkin(5) * t805 + t975 * t826 + t968 * t827 + t962 * t863 + t966 * t877 + t960 * t900;
t961 = t975 * t862 + t968 * t863 + t966 * t900;
t792 = mrSges(6,2) * t819 + mrSges(7,2) * t808 - mrSges(6,3) * t809 - mrSges(7,3) * t812 - qJ(6) * t805 - t968 * t826 + t976 * t827 + t962 * t862 + t967 * t877 - t961 * t900;
t856 = Ifges(5,5) * t893 + Ifges(5,6) * t892 + Ifges(5,3) * t900;
t858 = Ifges(5,1) * t893 + Ifges(5,4) * t892 + Ifges(5,5) * t900;
t767 = -mrSges(5,1) * t829 + mrSges(5,3) * t818 + Ifges(5,4) * t851 + Ifges(5,2) * t850 + Ifges(5,6) * t877 - pkin(4) * t802 + qJ(5) * t949 + t965 * t791 + t927 * t792 - t893 * t856 + t900 * t858;
t857 = Ifges(5,4) * t893 + Ifges(5,2) * t892 + Ifges(5,6) * t900;
t769 = mrSges(5,2) * t829 - mrSges(5,3) * t817 + Ifges(5,1) * t851 + Ifges(5,4) * t850 + Ifges(5,5) * t877 - qJ(5) * t790 - t927 * t791 + t965 * t792 + t892 * t856 - t900 * t857;
t885 = Ifges(4,4) * t905 - Ifges(4,2) * t904 + Ifges(4,6) * t924;
t886 = Ifges(4,1) * t905 - Ifges(4,4) * t904 + Ifges(4,5) * t924;
t941 = -mrSges(4,1) * t844 + mrSges(4,2) * t845 - Ifges(4,5) * t879 - Ifges(4,6) * t878 - Ifges(4,3) * t923 - pkin(3) * t939 - pkin(9) * t950 - t932 * t767 - t928 * t769 - t905 * t885 - t904 * t886;
t973 = mrSges(3,1) * t894 - mrSges(3,2) * t895 + Ifges(3,5) * t912 + Ifges(3,6) * t913 + Ifges(3,3) * qJDD(2) + pkin(2) * t775 + (t930 * t902 - t934 * t903) * qJD(1) - t941;
t804 = t827 * mrSges(7,2) + t863 * t841 - t948;
t972 = -t966 * t826 + t967 * t827 + t960 * t862 + t961 * t863 + (Ifges(5,3) - t974) * t877 + mrSges(5,1) * t817 + mrSges(6,1) * t809 - mrSges(7,1) * t808 - mrSges(5,2) * t818 - mrSges(6,2) * t810 + mrSges(7,3) * t807 + Ifges(5,5) * t851 + Ifges(5,6) * t850 + pkin(4) * t790 - pkin(5) * t804 + qJ(6) * (-t826 * mrSges(7,2) - t862 * t841 + t954) + t893 * t857 - t892 * t858;
t911 = (-mrSges(3,1) * t934 + mrSges(3,2) * t930) * qJD(1);
t957 = qJD(1) * t934;
t916 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t957;
t773 = m(3) * t894 + qJDD(2) * mrSges(3,1) - t912 * mrSges(3,3) + qJD(2) * t916 - t911 * t958 + t775;
t915 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t958;
t951 = t933 * t781 - t929 * t799;
t774 = m(3) * t895 - qJDD(2) * mrSges(3,2) + t913 * mrSges(3,3) - qJD(2) * t915 + t911 * t957 + t951;
t952 = -t930 * t773 + t934 * t774;
t762 = m(2) * t919 - t936 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t952;
t906 = -t936 * pkin(7) + t945;
t783 = t932 * t788 + t928 * t789;
t942 = m(4) * t880 - t878 * mrSges(4,1) + t879 * mrSges(4,2) + t904 * t896 + t905 * t897 + t783;
t940 = -m(3) * t906 + t913 * mrSges(3,1) - t912 * mrSges(3,2) - t915 * t958 + t916 * t957 - t942;
t777 = m(2) * t918 + qJDD(1) * mrSges(2,1) - t936 * mrSges(2,2) + t940;
t963 = t931 * t762 + t935 * t777;
t764 = t934 * t773 + t930 * t774;
t953 = t935 * t762 - t931 * t777;
t884 = Ifges(4,5) * t905 - Ifges(4,6) * t904 + Ifges(4,3) * t924;
t759 = mrSges(4,2) * t880 - mrSges(4,3) * t844 + Ifges(4,1) * t879 + Ifges(4,4) * t878 + Ifges(4,5) * t923 - pkin(9) * t783 - t928 * t767 + t932 * t769 - t904 * t884 - t924 * t885;
t765 = -mrSges(4,1) * t880 + mrSges(4,3) * t845 + Ifges(4,4) * t879 + Ifges(4,2) * t878 + Ifges(4,6) * t923 - pkin(3) * t783 - t905 * t884 + t924 * t886 - t972;
t901 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t930 + Ifges(3,6) * t934) * qJD(1);
t755 = -mrSges(3,1) * t906 + mrSges(3,3) * t895 + Ifges(3,4) * t912 + Ifges(3,2) * t913 + Ifges(3,6) * qJDD(2) - pkin(2) * t942 + pkin(8) * t951 + qJD(2) * t903 + t929 * t759 + t933 * t765 - t901 * t958;
t758 = mrSges(3,2) * t906 - mrSges(3,3) * t894 + Ifges(3,1) * t912 + Ifges(3,4) * t913 + Ifges(3,5) * qJDD(2) - pkin(8) * t775 - qJD(2) * t902 + t933 * t759 - t929 * t765 + t901 * t957;
t943 = mrSges(2,1) * t918 - mrSges(2,2) * t919 + Ifges(2,3) * qJDD(1) + pkin(1) * t940 + pkin(7) * t952 + t934 * t755 + t930 * t758;
t756 = mrSges(2,1) * g(3) + mrSges(2,3) * t919 + t936 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t764 - t973;
t753 = -mrSges(2,2) * g(3) - mrSges(2,3) * t918 + Ifges(2,5) * qJDD(1) - t936 * Ifges(2,6) - pkin(7) * t764 - t930 * t755 + t934 * t758;
t1 = [-m(1) * g(1) + t953; -m(1) * g(2) + t963; (-m(1) - m(2)) * g(3) + t764; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t963 + t935 * t753 - t931 * t756; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t953 + t931 * t753 + t935 * t756; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t943; t943; t973; -t941; t972; t802; t804;];
tauJB  = t1;
