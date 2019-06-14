% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-05-05 17:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:15:31
% EndTime: 2019-05-05 17:15:36
% DurationCPUTime: 4.96s
% Computational Cost: add. (45193->321), mult. (99421->374), div. (0->0), fcn. (61173->8), ass. (0->133)
t915 = Ifges(5,4) + Ifges(6,6);
t928 = -Ifges(5,2) - Ifges(6,3);
t923 = Ifges(5,6) - Ifges(6,5);
t927 = -2 * qJD(4);
t926 = Ifges(5,1) + Ifges(6,2);
t925 = -Ifges(6,1) - Ifges(5,3);
t924 = Ifges(5,5) - Ifges(6,4);
t872 = sin(qJ(3));
t875 = cos(qJ(3));
t901 = qJD(1) * qJD(3);
t847 = -t872 * qJDD(1) - t875 * t901;
t904 = qJD(1) * t875;
t851 = qJD(3) * pkin(3) - qJ(4) * t904;
t867 = t872 ^ 2;
t878 = qJD(1) ^ 2;
t873 = sin(qJ(1));
t876 = cos(qJ(1));
t854 = -t876 * g(1) - t873 * g(2);
t892 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t854;
t918 = -pkin(1) - pkin(7);
t777 = -t847 * pkin(3) + qJDD(4) + t851 * t904 + (-qJ(4) * t867 + t918) * t878 + t892;
t899 = t872 * t901;
t848 = t875 * qJDD(1) - t899;
t870 = cos(pkin(9));
t911 = sin(pkin(9));
t805 = -t870 * t847 + t911 * t848;
t806 = t911 * t847 + t870 * t848;
t905 = qJD(1) * t872;
t835 = t870 * t904 - t911 * t905;
t818 = qJD(3) * mrSges(5,1) - t835 * mrSges(5,3);
t834 = (t870 * t872 + t911 * t875) * qJD(1);
t903 = qJD(3) * t834;
t919 = -2 * qJD(5);
t881 = (-t806 + t903) * qJ(5) + t777 + (qJD(3) * pkin(4) + t919) * t835;
t757 = t805 * pkin(4) + t881;
t820 = t835 * mrSges(6,1) + qJD(3) * mrSges(6,2);
t853 = t873 * g(1) - t876 * g(2);
t889 = -t878 * qJ(2) + qJDD(2) - t853;
t822 = t918 * qJDD(1) + t889;
t808 = t872 * g(3) + t875 * t822;
t774 = (-t848 - t899) * qJ(4) + (-t872 * t875 * t878 + qJDD(3)) * pkin(3) + t808;
t809 = -t875 * g(3) + t872 * t822;
t775 = -t867 * t878 * pkin(3) + t847 * qJ(4) - qJD(3) * t851 + t809;
t759 = t870 * t774 - t911 * t775 + t835 * t927;
t793 = t834 * pkin(4) - t835 * qJ(5);
t877 = qJD(3) ^ 2;
t755 = -qJDD(3) * pkin(4) - t877 * qJ(5) + t835 * t793 + qJDD(5) - t759;
t749 = (t834 * t835 - qJDD(3)) * pkin(8) + (t806 + t903) * pkin(5) + t755;
t821 = t835 * pkin(5) - qJD(3) * pkin(8);
t833 = t834 ^ 2;
t752 = t881 - t835 * t821 - t833 * pkin(5) + (pkin(4) + pkin(8)) * t805;
t871 = sin(qJ(6));
t874 = cos(qJ(6));
t747 = t874 * t749 - t871 * t752;
t810 = -t871 * qJD(3) + t874 * t834;
t770 = t810 * qJD(6) + t874 * qJDD(3) + t871 * t805;
t811 = t874 * qJD(3) + t871 * t834;
t780 = -t810 * mrSges(7,1) + t811 * mrSges(7,2);
t831 = qJD(6) + t835;
t783 = -t831 * mrSges(7,2) + t810 * mrSges(7,3);
t804 = qJDD(6) + t806;
t744 = m(7) * t747 + t804 * mrSges(7,1) - t770 * mrSges(7,3) - t811 * t780 + t831 * t783;
t748 = t871 * t749 + t874 * t752;
t769 = -t811 * qJD(6) - t871 * qJDD(3) + t874 * t805;
t784 = t831 * mrSges(7,1) - t811 * mrSges(7,3);
t745 = m(7) * t748 - t804 * mrSges(7,2) + t769 * mrSges(7,3) + t810 * t780 - t831 * t784;
t896 = -t871 * t744 + t874 * t745;
t891 = -m(6) * t757 + t806 * mrSges(6,3) + t835 * t820 - t896;
t819 = t834 * mrSges(6,1) - qJD(3) * mrSges(6,3);
t906 = -qJD(3) * mrSges(5,2) - t834 * mrSges(5,3) - t819;
t916 = mrSges(5,1) - mrSges(6,2);
t731 = m(5) * t777 + t806 * mrSges(5,2) + t916 * t805 + t835 * t818 + t906 * t834 - t891;
t816 = t918 * t878 + t892;
t850 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t905;
t852 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t904;
t922 = -m(4) * t816 + t847 * mrSges(4,1) - t848 * mrSges(4,2) - t850 * t905 - t852 * t904 - t731;
t921 = t923 * qJD(3) + t928 * t834 + t915 * t835;
t917 = mrSges(2,1) - mrSges(3,2);
t914 = Ifges(2,5) - Ifges(3,4);
t913 = -Ifges(2,6) + Ifges(3,5);
t794 = t834 * mrSges(5,1) + t835 * mrSges(5,2);
t735 = t874 * t744 + t871 * t745;
t795 = -t834 * mrSges(6,2) - t835 * mrSges(6,3);
t886 = -m(6) * t755 - t806 * mrSges(6,1) - t835 * t795 - t735;
t730 = m(5) * t759 - t806 * mrSges(5,3) + t906 * qJD(3) + t916 * qJDD(3) - t835 * t794 + t886;
t828 = t834 * t927;
t909 = t911 * t774 + t870 * t775;
t760 = t828 + t909;
t890 = t877 * pkin(4) - qJDD(3) * qJ(5) - t909;
t753 = qJD(3) * t919 + ((2 * qJD(4)) + t793) * t834 + t890;
t751 = -t805 * pkin(5) - t833 * pkin(8) - t834 * t793 + t828 + ((2 * qJD(5)) + t821) * qJD(3) - t890;
t888 = -m(7) * t751 + t769 * mrSges(7,1) - t770 * mrSges(7,2) + t810 * t783 - t811 * t784;
t883 = -m(6) * t753 + qJDD(3) * mrSges(6,3) + qJD(3) * t820 - t888;
t741 = m(5) * t760 - qJDD(3) * mrSges(5,2) - qJD(3) * t818 + (-t794 - t795) * t834 + (-mrSges(5,3) - mrSges(6,1)) * t805 + t883;
t724 = t870 * t730 + t911 * t741;
t846 = (mrSges(4,1) * t872 + mrSges(4,2) * t875) * qJD(1);
t721 = m(4) * t808 + qJDD(3) * mrSges(4,1) - t848 * mrSges(4,3) + qJD(3) * t850 - t846 * t904 + t724;
t894 = -t911 * t730 + t870 * t741;
t722 = m(4) * t809 - qJDD(3) * mrSges(4,2) + t847 * mrSges(4,3) - qJD(3) * t852 - t846 * t905 + t894;
t718 = t875 * t721 + t872 * t722;
t832 = -qJDD(1) * pkin(1) + t889;
t887 = -m(3) * t832 + t878 * mrSges(3,3) - t718;
t714 = m(2) * t853 - t878 * mrSges(2,2) + t917 * qJDD(1) + t887;
t825 = t878 * pkin(1) - t892;
t880 = -m(3) * t825 + t878 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t922;
t727 = m(2) * t854 - t878 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t880;
t910 = t876 * t714 + t873 * t727;
t908 = qJD(3) * t925 + t834 * t923 - t835 * t924;
t907 = qJD(3) * t924 - t834 * t915 + t835 * t926;
t898 = -t873 * t714 + t876 * t727;
t897 = -t872 * t721 + t875 * t722;
t763 = Ifges(7,4) * t811 + Ifges(7,2) * t810 + Ifges(7,6) * t831;
t764 = Ifges(7,1) * t811 + Ifges(7,4) * t810 + Ifges(7,5) * t831;
t884 = mrSges(7,1) * t747 - mrSges(7,2) * t748 + Ifges(7,5) * t770 + Ifges(7,6) * t769 + Ifges(7,3) * t804 + t811 * t763 - t810 * t764;
t732 = -t805 * mrSges(6,2) - t834 * t819 - t891;
t762 = Ifges(7,5) * t811 + Ifges(7,6) * t810 + Ifges(7,3) * t831;
t737 = -mrSges(7,1) * t751 + mrSges(7,3) * t748 + Ifges(7,4) * t770 + Ifges(7,2) * t769 + Ifges(7,6) * t804 - t811 * t762 + t831 * t764;
t738 = mrSges(7,2) * t751 - mrSges(7,3) * t747 + Ifges(7,1) * t770 + Ifges(7,4) * t769 + Ifges(7,5) * t804 + t810 * t762 - t831 * t763;
t712 = -mrSges(5,1) * t777 - mrSges(6,1) * t753 + mrSges(6,2) * t757 + mrSges(5,3) * t760 - pkin(4) * t732 - pkin(5) * t888 - pkin(8) * t896 + t907 * qJD(3) + t923 * qJDD(3) - t874 * t737 - t871 * t738 + t928 * t805 + t915 * t806 + t908 * t835;
t719 = mrSges(6,1) * t755 + mrSges(5,2) * t777 - mrSges(5,3) * t759 - mrSges(6,3) * t757 + pkin(5) * t735 - qJ(5) * t732 - qJD(3) * t921 + qJDD(3) * t924 - t915 * t805 + t806 * t926 + t908 * t834 + t884;
t836 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t875 - Ifges(4,6) * t872) * qJD(1);
t838 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t875 - Ifges(4,4) * t872) * qJD(1);
t709 = -mrSges(4,1) * t816 + mrSges(4,3) * t809 + Ifges(4,4) * t848 + Ifges(4,2) * t847 + Ifges(4,6) * qJDD(3) - pkin(3) * t731 + qJ(4) * t894 + qJD(3) * t838 + t870 * t712 + t911 * t719 - t836 * t904;
t837 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t875 - Ifges(4,2) * t872) * qJD(1);
t711 = mrSges(4,2) * t816 - mrSges(4,3) * t808 + Ifges(4,1) * t848 + Ifges(4,4) * t847 + Ifges(4,5) * qJDD(3) - qJ(4) * t724 - qJD(3) * t837 - t911 * t712 + t870 * t719 - t836 * t905;
t716 = qJDD(1) * mrSges(3,2) - t887;
t882 = mrSges(2,1) * t853 - mrSges(2,2) * t854 + mrSges(3,2) * t832 - mrSges(3,3) * t825 - pkin(1) * t716 - pkin(7) * t718 + qJ(2) * t880 - t872 * t709 + t875 * t711 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t734 = qJDD(3) * mrSges(6,2) + qJD(3) * t819 - t886;
t879 = -t871 * t737 + t874 * t738 + Ifges(4,6) * t847 + Ifges(4,5) * t848 + mrSges(4,1) * t808 - mrSges(4,2) * t809 + mrSges(5,1) * t759 - mrSges(5,2) * t760 - mrSges(6,3) * t753 + mrSges(6,2) * t755 - pkin(4) * t734 - pkin(8) * t735 + pkin(3) * t724 + t907 * t834 + qJ(5) * (-t834 * t795 + t883) + t837 * t904 + t838 * t905 + t921 * t835 + t924 * t806 + (-qJ(5) * mrSges(6,1) - t923) * t805 + (Ifges(4,3) - t925) * qJDD(3);
t717 = -m(3) * g(3) + t897;
t708 = t879 - mrSges(2,3) * t853 + mrSges(3,1) * t832 + pkin(2) * t718 - qJ(2) * t717 + t914 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t913 * t878;
t707 = -mrSges(3,1) * t825 + mrSges(2,3) * t854 - pkin(1) * t717 - pkin(2) * t922 - pkin(7) * t897 + t917 * g(3) - t913 * qJDD(1) - t875 * t709 - t872 * t711 + t914 * t878;
t1 = [-m(1) * g(1) + t898; -m(1) * g(2) + t910; (-m(1) - m(2) - m(3)) * g(3) + t897; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t910 - t873 * t707 + t876 * t708; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t898 + t876 * t707 + t873 * t708; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t882; t882; t716; t879; t731; t734; t884;];
tauJB  = t1;
