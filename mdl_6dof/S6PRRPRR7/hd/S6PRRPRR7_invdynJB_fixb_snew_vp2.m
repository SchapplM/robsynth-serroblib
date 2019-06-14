% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 06:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:01:17
% EndTime: 2019-05-05 06:01:26
% DurationCPUTime: 8.20s
% Computational Cost: add. (133395->326), mult. (264803->401), div. (0->0), fcn. (167716->12), ass. (0->147)
t969 = -2 * qJD(4);
t959 = Ifges(4,4) + Ifges(5,6);
t958 = Ifges(4,5) - Ifges(5,4);
t968 = Ifges(4,2) + Ifges(5,3);
t967 = Ifges(5,2) + Ifges(4,1);
t957 = Ifges(4,6) - Ifges(5,5);
t966 = Ifges(4,3) + Ifges(5,1);
t909 = sin(pkin(11));
t911 = cos(pkin(11));
t889 = t909 * g(1) - t911 * g(2);
t908 = -g(3) + qJDD(1);
t910 = sin(pkin(6));
t912 = cos(pkin(6));
t857 = -t910 * t889 + t912 * t908;
t915 = sin(qJ(3));
t919 = cos(qJ(3));
t945 = qJD(2) * qJD(3);
t942 = t919 * t945;
t885 = t915 * qJDD(2) + t942;
t890 = -t911 * g(1) - t909 * g(2);
t920 = cos(qJ(2));
t916 = sin(qJ(2));
t953 = t912 * t916;
t954 = t910 * t916;
t840 = t889 * t953 + t920 * t890 + t908 * t954;
t922 = qJD(2) ^ 2;
t836 = -t922 * pkin(2) + qJDD(2) * pkin(8) + t840;
t833 = t915 * t836;
t882 = (-pkin(3) * t919 - qJ(4) * t915) * qJD(2);
t902 = t915 * qJD(2);
t921 = qJD(3) ^ 2;
t936 = -t921 * qJ(4) + t882 * t902 + qJDD(4) + t833;
t961 = pkin(9) * t922;
t962 = -pkin(3) - pkin(9);
t813 = t885 * pkin(4) + t962 * qJDD(3) + (-pkin(4) * t945 - t915 * t961 - t857) * t919 + t936;
t943 = t915 * t945;
t886 = t919 * qJDD(2) - t943;
t896 = pkin(4) * t902 - qJD(3) * pkin(9);
t907 = t919 ^ 2;
t839 = -t916 * t890 + (t889 * t912 + t908 * t910) * t920;
t929 = -qJDD(2) * pkin(2) - t839;
t927 = pkin(3) * t943 + t902 * t969 + (-t885 - t942) * qJ(4) + t929;
t815 = -t896 * t902 + (-pkin(4) * t907 - pkin(8)) * t922 + t962 * t886 + t927;
t914 = sin(qJ(5));
t918 = cos(qJ(5));
t805 = t918 * t813 - t914 * t815;
t946 = qJD(2) * t919;
t880 = -t914 * qJD(3) - t918 * t946;
t849 = t880 * qJD(5) + t918 * qJDD(3) - t914 * t886;
t877 = qJDD(5) + t885;
t881 = t918 * qJD(3) - t914 * t946;
t899 = t902 + qJD(5);
t802 = (t880 * t899 - t849) * pkin(10) + (t880 * t881 + t877) * pkin(5) + t805;
t806 = t914 * t813 + t918 * t815;
t848 = -t881 * qJD(5) - t914 * qJDD(3) - t918 * t886;
t856 = t899 * pkin(5) - t881 * pkin(10);
t876 = t880 ^ 2;
t803 = -t876 * pkin(5) + t848 * pkin(10) - t899 * t856 + t806;
t913 = sin(qJ(6));
t917 = cos(qJ(6));
t801 = t913 * t802 + t917 * t803;
t830 = t919 * t836 + t915 * t857;
t822 = t921 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t969 - t882 * t946 - t830;
t812 = t886 * pkin(4) + qJD(3) * t896 - t907 * t961 - t822;
t808 = -t848 * pkin(5) - t876 * pkin(10) + t881 * t856 + t812;
t851 = t913 * t880 + t917 * t881;
t820 = -t851 * qJD(6) + t917 * t848 - t913 * t849;
t850 = t917 * t880 - t913 * t881;
t821 = t850 * qJD(6) + t913 * t848 + t917 * t849;
t897 = qJD(6) + t899;
t825 = Ifges(7,5) * t851 + Ifges(7,6) * t850 + Ifges(7,3) * t897;
t827 = Ifges(7,1) * t851 + Ifges(7,4) * t850 + Ifges(7,5) * t897;
t870 = qJDD(6) + t877;
t788 = -mrSges(7,1) * t808 + mrSges(7,3) * t801 + Ifges(7,4) * t821 + Ifges(7,2) * t820 + Ifges(7,6) * t870 - t851 * t825 + t897 * t827;
t800 = t917 * t802 - t913 * t803;
t826 = Ifges(7,4) * t851 + Ifges(7,2) * t850 + Ifges(7,6) * t897;
t789 = mrSges(7,2) * t808 - mrSges(7,3) * t800 + Ifges(7,1) * t821 + Ifges(7,4) * t820 + Ifges(7,5) * t870 + t850 * t825 - t897 * t826;
t841 = Ifges(6,5) * t881 + Ifges(6,6) * t880 + Ifges(6,3) * t899;
t843 = Ifges(6,1) * t881 + Ifges(6,4) * t880 + Ifges(6,5) * t899;
t837 = -t897 * mrSges(7,2) + t850 * mrSges(7,3);
t838 = t897 * mrSges(7,1) - t851 * mrSges(7,3);
t932 = m(7) * t808 - t820 * mrSges(7,1) + t821 * mrSges(7,2) - t850 * t837 + t851 * t838;
t831 = -t850 * mrSges(7,1) + t851 * mrSges(7,2);
t797 = m(7) * t800 + t870 * mrSges(7,1) - t821 * mrSges(7,3) - t851 * t831 + t897 * t837;
t798 = m(7) * t801 - t870 * mrSges(7,2) + t820 * mrSges(7,3) + t850 * t831 - t897 * t838;
t939 = -t913 * t797 + t917 * t798;
t772 = -mrSges(6,1) * t812 + mrSges(6,3) * t806 + Ifges(6,4) * t849 + Ifges(6,2) * t848 + Ifges(6,6) * t877 - pkin(5) * t932 + pkin(10) * t939 + t917 * t788 + t913 * t789 - t881 * t841 + t899 * t843;
t787 = t917 * t797 + t913 * t798;
t842 = Ifges(6,4) * t881 + Ifges(6,2) * t880 + Ifges(6,6) * t899;
t773 = mrSges(6,2) * t812 - mrSges(6,3) * t805 + Ifges(6,1) * t849 + Ifges(6,4) * t848 + Ifges(6,5) * t877 - pkin(10) * t787 - t913 * t788 + t917 * t789 + t880 * t841 - t899 * t842;
t883 = (mrSges(5,2) * t919 - mrSges(5,3) * t915) * qJD(2);
t893 = -mrSges(5,1) * t946 - qJD(3) * mrSges(5,3);
t852 = -t880 * mrSges(6,1) + t881 * mrSges(6,2);
t854 = -t899 * mrSges(6,2) + t880 * mrSges(6,3);
t784 = m(6) * t805 + t877 * mrSges(6,1) - t849 * mrSges(6,3) - t881 * t852 + t899 * t854 + t787;
t855 = t899 * mrSges(6,1) - t881 * mrSges(6,3);
t785 = m(6) * t806 - t877 * mrSges(6,2) + t848 * mrSges(6,3) + t880 * t852 - t899 * t855 + t939;
t781 = t918 * t784 + t914 * t785;
t952 = t919 * t857;
t823 = -qJDD(3) * pkin(3) + t936 - t952;
t931 = -m(5) * t823 - t885 * mrSges(5,1) - t781;
t780 = qJDD(3) * mrSges(5,2) + qJD(3) * t893 + t883 * t902 - t931;
t829 = -t833 + t952;
t894 = mrSges(5,1) * t902 + qJD(3) * mrSges(5,2);
t928 = -m(6) * t812 + t848 * mrSges(6,1) - t849 * mrSges(6,2) + t880 * t854 - t881 * t855 - t932;
t925 = -m(5) * t822 + qJDD(3) * mrSges(5,3) + qJD(3) * t894 + t883 * t946 - t928;
t947 = t958 * qJD(3) + (t967 * t915 + t959 * t919) * qJD(2);
t948 = t957 * qJD(3) + (t959 * t915 + t968 * t919) * qJD(2);
t965 = (t915 * t948 - t919 * t947) * qJD(2) + t966 * qJDD(3) + t958 * t885 + t957 * t886 + mrSges(4,1) * t829 - mrSges(4,2) * t830 + mrSges(5,2) * t823 - mrSges(5,3) * t822 - pkin(3) * t780 - pkin(9) * t781 + qJ(4) * (t886 * mrSges(5,1) + t925) - t914 * t772 + t918 * t773;
t960 = t922 * pkin(8);
t835 = t929 - t960;
t891 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t902;
t892 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t946;
t824 = -t886 * pkin(3) + t927 - t960;
t950 = -t914 * t784 + t918 * t785;
t935 = -m(5) * t824 - t886 * mrSges(5,2) + t894 * t902 - t950;
t924 = -m(4) * t835 + t892 * t946 + t886 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t885 + (-t891 * t915 - t893 * t919) * qJD(2) + t935;
t776 = m(3) * t839 + qJDD(2) * mrSges(3,1) - t922 * mrSges(3,2) + t924;
t955 = t776 * t920;
t884 = (-mrSges(4,1) * t919 + mrSges(4,2) * t915) * qJD(2);
t778 = m(4) * t829 - t885 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t892 - t893) * qJD(3) + (-t883 - t884) * t902 + t931;
t792 = t925 + t884 * t946 - qJD(3) * t891 + m(4) * t830 - qJDD(3) * mrSges(4,2) + (mrSges(4,3) + mrSges(5,1)) * t886;
t940 = -t915 * t778 + t919 * t792;
t768 = m(3) * t840 - t922 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t940;
t771 = t919 * t778 + t915 * t792;
t770 = m(3) * t857 + t771;
t759 = t768 * t953 - t910 * t770 + t912 * t955;
t757 = m(2) * t889 + t759;
t764 = t920 * t768 - t916 * t776;
t763 = m(2) * t890 + t764;
t951 = t911 * t757 + t909 * t763;
t949 = t966 * qJD(3) + (t958 * t915 + t957 * t919) * qJD(2);
t758 = t768 * t954 + t912 * t770 + t910 * t955;
t941 = -t909 * t757 + t911 * t763;
t938 = m(2) * t908 + t758;
t779 = -t885 * mrSges(5,3) + t893 * t946 - t935;
t755 = -mrSges(4,1) * t835 - mrSges(5,1) * t822 + mrSges(5,2) * t824 + mrSges(4,3) * t830 - pkin(3) * t779 - pkin(4) * t928 - pkin(9) * t950 + t947 * qJD(3) + t957 * qJDD(3) - t918 * t772 - t914 * t773 + t959 * t885 + t968 * t886 - t949 * t902;
t930 = mrSges(7,1) * t800 - mrSges(7,2) * t801 + Ifges(7,5) * t821 + Ifges(7,6) * t820 + Ifges(7,3) * t870 + t851 * t826 - t850 * t827;
t926 = mrSges(6,1) * t805 - mrSges(6,2) * t806 + Ifges(6,5) * t849 + Ifges(6,6) * t848 + Ifges(6,3) * t877 + pkin(5) * t787 + t881 * t842 - t880 * t843 + t930;
t760 = mrSges(5,1) * t823 + mrSges(4,2) * t835 - mrSges(4,3) * t829 - mrSges(5,3) * t824 + pkin(4) * t781 - qJ(4) * t779 - t948 * qJD(3) + t958 * qJDD(3) + t967 * t885 + t959 * t886 + t949 * t946 + t926;
t753 = mrSges(3,2) * t857 - mrSges(3,3) * t839 + Ifges(3,5) * qJDD(2) - t922 * Ifges(3,6) - pkin(8) * t771 - t915 * t755 + t919 * t760;
t754 = -mrSges(3,1) * t857 + mrSges(3,3) * t840 + t922 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t771 - t965;
t933 = pkin(7) * t764 + t753 * t916 + t754 * t920;
t752 = mrSges(3,1) * t839 - mrSges(3,2) * t840 + Ifges(3,3) * qJDD(2) + pkin(2) * t924 + pkin(8) * t940 + t919 * t755 + t915 * t760;
t751 = mrSges(2,2) * t908 - mrSges(2,3) * t889 + t920 * t753 - t916 * t754 + (-t758 * t910 - t759 * t912) * pkin(7);
t750 = -mrSges(2,1) * t908 + mrSges(2,3) * t890 - pkin(1) * t758 - t910 * t752 + t912 * t933;
t1 = [-m(1) * g(1) + t941; -m(1) * g(2) + t951; -m(1) * g(3) + t938; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t951 - t909 * t750 + t911 * t751; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t941 + t911 * t750 + t909 * t751; -mrSges(1,1) * g(2) + mrSges(2,1) * t889 + mrSges(1,2) * g(1) - mrSges(2,2) * t890 + pkin(1) * t759 + t912 * t752 + t910 * t933; t938; t752; t965; t780; t926; t930;];
tauJB  = t1;
