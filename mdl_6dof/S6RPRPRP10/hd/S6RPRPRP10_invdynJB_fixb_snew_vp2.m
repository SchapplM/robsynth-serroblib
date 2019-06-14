% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-05-05 18:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRP10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:09:53
% EndTime: 2019-05-05 18:09:59
% DurationCPUTime: 2.98s
% Computational Cost: add. (21643->305), mult. (42247->345), div. (0->0), fcn. (20963->6), ass. (0->128)
t937 = Ifges(6,1) + Ifges(7,1);
t916 = Ifges(6,4) - Ifges(7,5);
t934 = Ifges(7,4) + Ifges(6,5);
t936 = Ifges(6,2) + Ifges(7,3);
t931 = Ifges(6,6) - Ifges(7,6);
t935 = -2 * qJD(4);
t933 = Ifges(4,5) - Ifges(5,4);
t932 = Ifges(4,6) - Ifges(5,5);
t930 = Ifges(4,3) + Ifges(5,1);
t929 = Ifges(6,3) + Ifges(7,2);
t871 = sin(qJ(5));
t872 = sin(qJ(3));
t908 = qJD(1) * t872;
t922 = cos(qJ(5));
t833 = t871 * qJD(3) - t922 * t908;
t834 = t922 * qJD(3) + t871 * t908;
t874 = cos(qJ(3));
t907 = qJD(1) * t874;
t854 = qJD(5) + t907;
t928 = t936 * t833 - t916 * t834 - t931 * t854;
t927 = -t916 * t833 + t937 * t834 + t934 * t854;
t873 = sin(qJ(1));
t875 = cos(qJ(1));
t849 = -t875 * g(1) - t873 * g(2);
t891 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t849;
t877 = qJD(1) ^ 2;
t924 = -pkin(1) - pkin(7);
t903 = t924 * t877;
t804 = t903 + t891;
t906 = qJD(1) * qJD(3);
t901 = t874 * t906;
t838 = qJDD(1) * t872 + t901;
t856 = t872 * t906;
t839 = qJDD(1) * t874 - t856;
t843 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t908;
t844 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t907;
t885 = pkin(3) * t901 + t907 * t935 + t891 + (-t839 + t856) * qJ(4);
t768 = t838 * pkin(3) + t885 + t903;
t845 = mrSges(5,1) * t908 - qJD(3) * mrSges(5,3);
t846 = mrSges(5,1) * t907 + qJD(3) * mrSges(5,2);
t847 = pkin(4) * t907 - qJD(3) * pkin(8);
t868 = t872 ^ 2;
t923 = pkin(3) + pkin(8);
t762 = -t847 * t907 + t923 * t838 + (-pkin(4) * t868 + t924) * t877 + t885;
t835 = (pkin(3) * t872 - qJ(4) * t874) * qJD(1);
t876 = qJD(3) ^ 2;
t848 = t873 * g(1) - t875 * g(2);
t890 = -t877 * qJ(2) + qJDD(2) - t848;
t805 = t924 * qJDD(1) + t890;
t912 = t874 * t805;
t889 = -t876 * qJ(4) + t835 * t907 + qJDD(4) - t912;
t921 = pkin(8) * t877;
t766 = t839 * pkin(4) - t923 * qJDD(3) + (pkin(4) * t906 + t874 * t921 - g(3)) * t872 + t889;
t760 = t922 * t762 + t871 * t766;
t788 = t834 * qJD(5) + t871 * qJDD(3) - t922 * t838;
t799 = mrSges(6,1) * t854 - mrSges(6,3) * t834;
t832 = qJDD(5) + t839;
t793 = pkin(5) * t833 - qJ(6) * t834;
t850 = t854 ^ 2;
t754 = -pkin(5) * t850 + t832 * qJ(6) + 0.2e1 * qJD(6) * t854 - t793 * t833 + t760;
t800 = -mrSges(7,1) * t854 + mrSges(7,2) * t834;
t902 = m(7) * t754 + t832 * mrSges(7,3) + t854 * t800;
t794 = mrSges(7,1) * t833 - mrSges(7,3) * t834;
t909 = -mrSges(6,1) * t833 - mrSges(6,2) * t834 - t794;
t917 = -mrSges(6,3) - mrSges(7,2);
t745 = m(6) * t760 - t832 * mrSges(6,2) + t917 * t788 - t854 * t799 + t909 * t833 + t902;
t759 = -t871 * t762 + t922 * t766;
t789 = -t833 * qJD(5) + t922 * qJDD(3) + t871 * t838;
t798 = -mrSges(6,2) * t854 - mrSges(6,3) * t833;
t755 = -t832 * pkin(5) - t850 * qJ(6) + t834 * t793 + qJDD(6) - t759;
t801 = -mrSges(7,2) * t833 + mrSges(7,3) * t854;
t895 = -m(7) * t755 + t832 * mrSges(7,1) + t854 * t801;
t746 = m(6) * t759 + t832 * mrSges(6,1) + t917 * t789 + t854 * t798 + t909 * t834 + t895;
t898 = t922 * t745 - t871 * t746;
t884 = m(5) * t768 - t839 * mrSges(5,3) - (t845 * t872 + t846 * t874) * qJD(1) + t898;
t918 = mrSges(4,1) - mrSges(5,2);
t926 = -m(4) * t804 - t839 * mrSges(4,2) - t918 * t838 - t843 * t908 - t844 * t907 - t884;
t797 = -t874 * g(3) + t872 * t805;
t769 = t876 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t935 + t835 * t908 - t797;
t920 = t872 * g(3);
t919 = mrSges(2,1) - mrSges(3,2);
t915 = Ifges(2,5) - Ifges(3,4);
t914 = -Ifges(2,6) + Ifges(3,5);
t913 = -Ifges(5,6) - Ifges(4,4);
t796 = t912 + t920;
t741 = t871 * t745 + t922 * t746;
t771 = -qJDD(3) * pkin(3) + t889 - t920;
t887 = m(5) * t771 + t839 * mrSges(5,1) + t741;
t836 = (-mrSges(5,2) * t872 - mrSges(5,3) * t874) * qJD(1);
t896 = qJD(1) * (-t836 - (mrSges(4,1) * t872 + mrSges(4,2) * t874) * qJD(1));
t734 = m(4) * t796 - t839 * mrSges(4,3) + t918 * qJDD(3) + (t843 - t845) * qJD(3) + t874 * t896 - t887;
t765 = -t838 * pkin(4) + qJD(3) * t847 - t868 * t921 - t769;
t757 = -0.2e1 * qJD(6) * t834 + (t833 * t854 - t789) * qJ(6) + (t834 * t854 + t788) * pkin(5) + t765;
t751 = m(7) * t757 + t788 * mrSges(7,1) - t789 * mrSges(7,3) - t834 * t800 + t833 * t801;
t883 = m(6) * t765 + t788 * mrSges(6,1) + t789 * mrSges(6,2) + t833 * t798 + t834 * t799 + t751;
t881 = -m(5) * t769 + qJDD(3) * mrSges(5,3) + qJD(3) * t846 + t883;
t743 = t881 + (-mrSges(4,3) - mrSges(5,1)) * t838 - qJDD(3) * mrSges(4,2) + t872 * t896 - qJD(3) * t844 + m(4) * t797;
t729 = t874 * t734 + t872 * t743;
t810 = -qJDD(1) * pkin(1) + t890;
t888 = -m(3) * t810 + t877 * mrSges(3,3) - t729;
t725 = m(2) * t848 - t877 * mrSges(2,2) + t919 * qJDD(1) + t888;
t808 = t877 * pkin(1) - t891;
t879 = -m(3) * t808 + t877 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t926;
t732 = m(2) * t849 - t877 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t879;
t911 = t875 * t725 + t873 * t732;
t910 = t931 * t833 - t934 * t834 - t929 * t854;
t900 = -t725 * t873 + t875 * t732;
t899 = -t872 * t734 + t874 * t743;
t897 = qJD(1) * (-t930 * qJD(3) + (t932 * t872 - t933 * t874) * qJD(1));
t735 = -t838 * mrSges(5,2) + t884;
t738 = -mrSges(6,1) * t765 - mrSges(7,1) * t757 + mrSges(7,2) * t754 + mrSges(6,3) * t760 - pkin(5) * t751 - t936 * t788 + t916 * t789 + t931 * t832 + t910 * t834 + t927 * t854;
t740 = mrSges(6,2) * t765 + mrSges(7,2) * t755 - mrSges(6,3) * t759 - mrSges(7,3) * t757 - qJ(6) * t751 - t916 * t788 + t937 * t789 + t934 * t832 + t910 * t833 + t928 * t854;
t814 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t874 - Ifges(4,4) * t872) * qJD(1);
t816 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t874 + Ifges(5,6) * t872) * qJD(1);
t721 = -mrSges(4,1) * t804 + mrSges(4,3) * t797 - mrSges(5,1) * t769 + mrSges(5,2) * t768 - t871 * t740 - t922 * t738 + pkin(4) * t883 - pkin(8) * t898 - pkin(3) * t735 - t913 * t839 + (-Ifges(4,2) - Ifges(5,3)) * t838 + t932 * qJDD(3) + (t814 - t816) * qJD(3) + t874 * t897;
t813 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t874 - Ifges(4,2) * t872) * qJD(1);
t815 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t874 + Ifges(5,3) * t872) * qJD(1);
t750 = t789 * mrSges(7,2) + t834 * t794 - t895;
t880 = mrSges(6,1) * t759 - mrSges(7,1) * t755 - mrSges(6,2) * t760 + mrSges(7,3) * t754 - pkin(5) * t750 + qJ(6) * t902 - t928 * t834 + (-qJ(6) * t794 + t927) * t833 + t929 * t832 + t934 * t789 + (-mrSges(7,2) * qJ(6) - t931) * t788;
t723 = t880 + t933 * qJDD(3) + (Ifges(5,2) + Ifges(4,1)) * t839 + t913 * t838 + (-t813 + t815) * qJD(3) + t872 * t897 - mrSges(4,3) * t796 + mrSges(4,2) * t804 - mrSges(5,3) * t768 + mrSges(5,1) * t771 + pkin(4) * t741 - qJ(4) * t735;
t727 = qJDD(1) * mrSges(3,2) - t888;
t882 = mrSges(2,1) * t848 - mrSges(2,2) * t849 + mrSges(3,2) * t810 - mrSges(3,3) * t808 - pkin(1) * t727 - pkin(7) * t729 + qJ(2) * t879 - t721 * t872 + t874 * t723 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t737 = qJDD(3) * mrSges(5,2) + qJD(3) * t845 + t836 * t907 + t887;
t878 = -mrSges(4,2) * t797 - mrSges(5,3) * t769 - pkin(8) * t741 - t871 * t738 - pkin(3) * t737 + t922 * t740 + qJ(4) * (-t836 * t908 + t881) + mrSges(5,2) * t771 + mrSges(4,1) * t796 + t814 * t908 + t813 * t907 + (-t815 * t874 - t816 * t872) * qJD(1) + t933 * t839 + (-mrSges(5,1) * qJ(4) - t932) * t838 + t930 * qJDD(3);
t728 = -m(3) * g(3) + t899;
t720 = t878 + t915 * qJDD(1) + t914 * t877 - mrSges(2,3) * t848 + mrSges(3,1) * t810 + pkin(2) * t729 - qJ(2) * t728 + (mrSges(3,3) - mrSges(2,2)) * g(3);
t719 = -mrSges(3,1) * t808 + mrSges(2,3) * t849 - pkin(1) * t728 - pkin(2) * t926 - pkin(7) * t899 + t919 * g(3) - t914 * qJDD(1) - t874 * t721 - t872 * t723 + t915 * t877;
t1 = [-m(1) * g(1) + t900; -m(1) * g(2) + t911; (-m(1) - m(2) - m(3)) * g(3) + t899; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t911 - t873 * t719 + t875 * t720; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t900 + t875 * t719 + t873 * t720; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t882; t882; t727; t878; t737; t880; t750;];
tauJB  = t1;
