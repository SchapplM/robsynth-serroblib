% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRR5
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 12:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:58:25
% EndTime: 2019-05-05 11:59:18
% DurationCPUTime: 50.39s
% Computational Cost: add. (920978->354), mult. (1894028->467), div. (0->0), fcn. (1534383->16), ass. (0->160)
t911 = sin(pkin(7));
t919 = sin(qJ(3));
t924 = cos(qJ(3));
t942 = qJD(2) * qJD(3);
t894 = (-qJDD(2) * t924 + t919 * t942) * t911;
t910 = sin(pkin(13));
t913 = cos(pkin(13));
t901 = g(1) * t910 - g(2) * t913;
t902 = -g(1) * t913 - g(2) * t910;
t909 = -g(3) + qJDD(1);
t920 = sin(qJ(2));
t915 = cos(pkin(6));
t925 = cos(qJ(2));
t945 = t915 * t925;
t912 = sin(pkin(6));
t948 = t912 * t925;
t868 = t901 * t945 - t902 * t920 + t909 * t948;
t926 = qJD(2) ^ 2;
t952 = pkin(9) * t911;
t864 = qJDD(2) * pkin(2) + t926 * t952 + t868;
t946 = t915 * t920;
t949 = t912 * t920;
t869 = t901 * t946 + t925 * t902 + t909 * t949;
t865 = -pkin(2) * t926 + qJDD(2) * t952 + t869;
t886 = -t901 * t912 + t909 * t915;
t914 = cos(pkin(7));
t832 = -t919 * t865 + (t864 * t914 + t886 * t911) * t924;
t907 = qJD(2) * t914 + qJD(3);
t943 = qJD(2) * t911;
t940 = t924 * t943;
t890 = -mrSges(4,2) * t907 + mrSges(4,3) * t940;
t891 = (-mrSges(4,1) * t924 + mrSges(4,2) * t919) * t943;
t893 = (qJDD(2) * t919 + t924 * t942) * t911;
t906 = qJDD(2) * t914 + qJDD(3);
t947 = t914 * t919;
t950 = t911 * t919;
t833 = t864 * t947 + t924 * t865 + t886 * t950;
t892 = (-pkin(3) * t924 - pkin(10) * t919) * t943;
t905 = t907 ^ 2;
t828 = -pkin(3) * t905 + pkin(10) * t906 + t892 * t940 + t833;
t880 = t914 * t886;
t831 = t894 * pkin(3) - t893 * pkin(10) + t880 + (-t864 + (pkin(3) * t919 - pkin(10) * t924) * t907 * qJD(2)) * t911;
t918 = sin(qJ(4));
t923 = cos(qJ(4));
t819 = t923 * t828 + t918 * t831;
t941 = t919 * t943;
t884 = t907 * t923 - t918 * t941;
t885 = t907 * t918 + t923 * t941;
t867 = -pkin(4) * t884 - pkin(11) * t885;
t888 = qJDD(4) + t894;
t900 = qJD(4) - t940;
t898 = t900 ^ 2;
t809 = -pkin(4) * t898 + pkin(11) * t888 + t867 * t884 + t819;
t827 = -t906 * pkin(3) - t905 * pkin(10) + t892 * t941 - t832;
t859 = -t885 * qJD(4) - t918 * t893 + t906 * t923;
t860 = qJD(4) * t884 + t893 * t923 + t906 * t918;
t812 = (-t884 * t900 - t860) * pkin(11) + (t885 * t900 - t859) * pkin(4) + t827;
t917 = sin(qJ(5));
t922 = cos(qJ(5));
t804 = -t917 * t809 + t922 * t812;
t871 = -t885 * t917 + t900 * t922;
t836 = qJD(5) * t871 + t860 * t922 + t888 * t917;
t857 = qJDD(5) - t859;
t872 = t885 * t922 + t900 * t917;
t883 = qJD(5) - t884;
t802 = (t871 * t883 - t836) * pkin(12) + (t871 * t872 + t857) * pkin(5) + t804;
t805 = t922 * t809 + t917 * t812;
t835 = -qJD(5) * t872 - t860 * t917 + t888 * t922;
t850 = pkin(5) * t883 - pkin(12) * t872;
t870 = t871 ^ 2;
t803 = -pkin(5) * t870 + pkin(12) * t835 - t850 * t883 + t805;
t916 = sin(qJ(6));
t921 = cos(qJ(6));
t800 = t802 * t921 - t803 * t916;
t843 = t871 * t921 - t872 * t916;
t817 = qJD(6) * t843 + t835 * t916 + t836 * t921;
t844 = t871 * t916 + t872 * t921;
t829 = -mrSges(7,1) * t843 + mrSges(7,2) * t844;
t881 = qJD(6) + t883;
t837 = -mrSges(7,2) * t881 + mrSges(7,3) * t843;
t852 = qJDD(6) + t857;
t796 = m(7) * t800 + mrSges(7,1) * t852 - mrSges(7,3) * t817 - t829 * t844 + t837 * t881;
t801 = t802 * t916 + t803 * t921;
t816 = -qJD(6) * t844 + t835 * t921 - t836 * t916;
t838 = mrSges(7,1) * t881 - mrSges(7,3) * t844;
t797 = m(7) * t801 - mrSges(7,2) * t852 + mrSges(7,3) * t816 + t829 * t843 - t838 * t881;
t788 = t921 * t796 + t916 * t797;
t845 = -mrSges(6,1) * t871 + mrSges(6,2) * t872;
t848 = -mrSges(6,2) * t883 + mrSges(6,3) * t871;
t786 = m(6) * t804 + mrSges(6,1) * t857 - mrSges(6,3) * t836 - t845 * t872 + t848 * t883 + t788;
t849 = mrSges(6,1) * t883 - mrSges(6,3) * t872;
t937 = -t796 * t916 + t921 * t797;
t787 = m(6) * t805 - mrSges(6,2) * t857 + mrSges(6,3) * t835 + t845 * t871 - t849 * t883 + t937;
t783 = t786 * t922 + t787 * t917;
t873 = -mrSges(5,2) * t900 + mrSges(5,3) * t884;
t874 = mrSges(5,1) * t900 - mrSges(5,3) * t885;
t929 = -m(5) * t827 + t859 * mrSges(5,1) - mrSges(5,2) * t860 + t884 * t873 - t874 * t885 - t783;
t779 = m(4) * t832 + mrSges(4,1) * t906 - mrSges(4,3) * t893 + t890 * t907 - t891 * t941 + t929;
t951 = t779 * t924;
t889 = mrSges(4,1) * t907 - mrSges(4,3) * t941;
t784 = -t786 * t917 + t922 * t787;
t866 = -mrSges(5,1) * t884 + mrSges(5,2) * t885;
t782 = m(5) * t819 - mrSges(5,2) * t888 + mrSges(5,3) * t859 + t866 * t884 - t874 * t900 + t784;
t818 = -t918 * t828 + t831 * t923;
t808 = -pkin(4) * t888 - pkin(11) * t898 + t885 * t867 - t818;
t806 = -pkin(5) * t835 - pkin(12) * t870 + t850 * t872 + t808;
t931 = m(7) * t806 - t816 * mrSges(7,1) + mrSges(7,2) * t817 - t843 * t837 + t838 * t844;
t798 = -m(6) * t808 + t835 * mrSges(6,1) - mrSges(6,2) * t836 + t871 * t848 - t849 * t872 - t931;
t792 = m(5) * t818 + mrSges(5,1) * t888 - mrSges(5,3) * t860 - t866 * t885 + t873 * t900 + t798;
t938 = t923 * t782 - t792 * t918;
t771 = m(4) * t833 - mrSges(4,2) * t906 - mrSges(4,3) * t894 - t889 * t907 + t891 * t940 + t938;
t774 = t918 * t782 + t923 * t792;
t846 = -t911 * t864 + t880;
t773 = m(4) * t846 + t894 * mrSges(4,1) + t893 * mrSges(4,2) + (t889 * t919 - t890 * t924) * t943 + t774;
t760 = t771 * t947 - t773 * t911 + t914 * t951;
t756 = m(3) * t868 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t926 + t760;
t759 = t771 * t950 + t914 * t773 + t911 * t951;
t758 = m(3) * t886 + t759;
t766 = t924 * t771 - t779 * t919;
t765 = m(3) * t869 - mrSges(3,1) * t926 - qJDD(2) * mrSges(3,2) + t766;
t746 = t756 * t945 - t758 * t912 + t765 * t946;
t744 = m(2) * t901 + t746;
t752 = -t756 * t920 + t925 * t765;
t751 = m(2) * t902 + t752;
t944 = t913 * t744 + t910 * t751;
t745 = t756 * t948 + t915 * t758 + t765 * t949;
t939 = -t744 * t910 + t913 * t751;
t936 = m(2) * t909 + t745;
t821 = Ifges(7,5) * t844 + Ifges(7,6) * t843 + Ifges(7,3) * t881;
t823 = Ifges(7,1) * t844 + Ifges(7,4) * t843 + Ifges(7,5) * t881;
t789 = -mrSges(7,1) * t806 + mrSges(7,3) * t801 + Ifges(7,4) * t817 + Ifges(7,2) * t816 + Ifges(7,6) * t852 - t821 * t844 + t823 * t881;
t822 = Ifges(7,4) * t844 + Ifges(7,2) * t843 + Ifges(7,6) * t881;
t790 = mrSges(7,2) * t806 - mrSges(7,3) * t800 + Ifges(7,1) * t817 + Ifges(7,4) * t816 + Ifges(7,5) * t852 + t821 * t843 - t822 * t881;
t839 = Ifges(6,5) * t872 + Ifges(6,6) * t871 + Ifges(6,3) * t883;
t841 = Ifges(6,1) * t872 + Ifges(6,4) * t871 + Ifges(6,5) * t883;
t775 = -mrSges(6,1) * t808 + mrSges(6,3) * t805 + Ifges(6,4) * t836 + Ifges(6,2) * t835 + Ifges(6,6) * t857 - pkin(5) * t931 + pkin(12) * t937 + t921 * t789 + t916 * t790 - t872 * t839 + t883 * t841;
t840 = Ifges(6,4) * t872 + Ifges(6,2) * t871 + Ifges(6,6) * t883;
t776 = mrSges(6,2) * t808 - mrSges(6,3) * t804 + Ifges(6,1) * t836 + Ifges(6,4) * t835 + Ifges(6,5) * t857 - pkin(12) * t788 - t789 * t916 + t790 * t921 + t839 * t871 - t840 * t883;
t853 = Ifges(5,5) * t885 + Ifges(5,6) * t884 + Ifges(5,3) * t900;
t854 = Ifges(5,4) * t885 + Ifges(5,2) * t884 + Ifges(5,6) * t900;
t761 = mrSges(5,2) * t827 - mrSges(5,3) * t818 + Ifges(5,1) * t860 + Ifges(5,4) * t859 + Ifges(5,5) * t888 - pkin(11) * t783 - t775 * t917 + t776 * t922 + t853 * t884 - t854 * t900;
t855 = Ifges(5,1) * t885 + Ifges(5,4) * t884 + Ifges(5,5) * t900;
t930 = -mrSges(7,1) * t800 + mrSges(7,2) * t801 - Ifges(7,5) * t817 - Ifges(7,6) * t816 - Ifges(7,3) * t852 - t844 * t822 + t843 * t823;
t927 = mrSges(6,1) * t804 - mrSges(6,2) * t805 + Ifges(6,5) * t836 + Ifges(6,6) * t835 + Ifges(6,3) * t857 + pkin(5) * t788 + t872 * t840 - t871 * t841 - t930;
t767 = -mrSges(5,1) * t827 + mrSges(5,3) * t819 + Ifges(5,4) * t860 + Ifges(5,2) * t859 + Ifges(5,6) * t888 - pkin(4) * t783 - t885 * t853 + t900 * t855 - t927;
t877 = Ifges(4,6) * t907 + (Ifges(4,4) * t919 + Ifges(4,2) * t924) * t943;
t878 = Ifges(4,5) * t907 + (Ifges(4,1) * t919 + Ifges(4,4) * t924) * t943;
t747 = Ifges(4,5) * t893 - Ifges(4,6) * t894 + Ifges(4,3) * t906 + mrSges(4,1) * t832 - mrSges(4,2) * t833 + t918 * t761 + t923 * t767 + pkin(3) * t929 + pkin(10) * t938 + (t877 * t919 - t878 * t924) * t943;
t876 = Ifges(4,3) * t907 + (Ifges(4,5) * t919 + Ifges(4,6) * t924) * t943;
t748 = mrSges(4,2) * t846 - mrSges(4,3) * t832 + Ifges(4,1) * t893 - Ifges(4,4) * t894 + Ifges(4,5) * t906 - pkin(10) * t774 + t761 * t923 - t767 * t918 + t876 * t940 - t877 * t907;
t928 = mrSges(5,1) * t818 - mrSges(5,2) * t819 + Ifges(5,5) * t860 + Ifges(5,6) * t859 + Ifges(5,3) * t888 + pkin(4) * t798 + pkin(11) * t784 + t922 * t775 + t917 * t776 + t885 * t854 - t884 * t855;
t753 = -mrSges(4,1) * t846 + mrSges(4,3) * t833 + Ifges(4,4) * t893 - Ifges(4,2) * t894 + Ifges(4,6) * t906 - pkin(3) * t774 - t876 * t941 + t907 * t878 - t928;
t932 = pkin(9) * t766 + t748 * t919 + t753 * t924;
t741 = -mrSges(3,1) * t886 + mrSges(3,3) * t869 + t926 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t759 - t911 * t747 + t914 * t932;
t742 = mrSges(3,2) * t886 - mrSges(3,3) * t868 + Ifges(3,5) * qJDD(2) - t926 * Ifges(3,6) + t924 * t748 - t919 * t753 + (-t759 * t911 - t760 * t914) * pkin(9);
t933 = pkin(8) * t752 + t741 * t925 + t742 * t920;
t740 = mrSges(3,1) * t868 - mrSges(3,2) * t869 + Ifges(3,3) * qJDD(2) + pkin(2) * t760 + t914 * t747 + t911 * t932;
t739 = mrSges(2,2) * t909 - mrSges(2,3) * t901 - t920 * t741 + t925 * t742 + (-t745 * t912 - t746 * t915) * pkin(8);
t738 = -mrSges(2,1) * t909 + mrSges(2,3) * t902 - pkin(1) * t745 - t912 * t740 + t915 * t933;
t1 = [-m(1) * g(1) + t939; -m(1) * g(2) + t944; -m(1) * g(3) + t936; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t944 - t910 * t738 + t913 * t739; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t939 + t913 * t738 + t910 * t739; -mrSges(1,1) * g(2) + mrSges(2,1) * t901 + mrSges(1,2) * g(1) - mrSges(2,2) * t902 + pkin(1) * t746 + t915 * t740 + t912 * t933; t936; t740; t747; t928; t927; -t930;];
tauJB  = t1;
