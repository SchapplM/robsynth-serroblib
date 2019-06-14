% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-05-04 23:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:42:15
% EndTime: 2019-05-04 23:42:25
% DurationCPUTime: 9.52s
% Computational Cost: add. (147679->298), mult. (325607->367), div. (0->0), fcn. (241562->12), ass. (0->140)
t915 = Ifges(6,1) + Ifges(7,1);
t909 = Ifges(6,4) + Ifges(7,4);
t908 = Ifges(6,5) + Ifges(7,5);
t914 = Ifges(6,2) + Ifges(7,2);
t907 = Ifges(6,6) + Ifges(7,6);
t913 = Ifges(6,3) + Ifges(7,3);
t869 = qJD(2) ^ 2;
t856 = sin(pkin(11));
t859 = cos(pkin(11));
t863 = sin(qJ(4));
t866 = cos(qJ(4));
t877 = t856 * t863 - t859 * t866;
t837 = t877 * qJD(2);
t857 = sin(pkin(10));
t860 = cos(pkin(10));
t844 = g(1) * t857 - g(2) * t860;
t845 = -g(1) * t860 - g(2) * t857;
t855 = -g(3) + qJDD(1);
t858 = sin(pkin(6));
t861 = cos(pkin(6));
t864 = sin(qJ(2));
t867 = cos(qJ(2));
t818 = -t864 * t845 + (t844 * t861 + t855 * t858) * t867;
t878 = t856 * t866 + t859 * t863;
t838 = t878 * qJD(2);
t892 = t838 * qJD(4);
t826 = -qJDD(2) * t877 - t892;
t893 = t837 * qJD(4);
t827 = qJDD(2) * t878 - t893;
t862 = sin(qJ(5));
t865 = cos(qJ(5));
t829 = qJD(4) * t865 - t838 * t862;
t801 = qJD(5) * t829 + qJDD(4) * t862 + t827 * t865;
t830 = qJD(4) * t862 + t838 * t865;
t804 = -mrSges(7,1) * t829 + mrSges(7,2) * t830;
t901 = t861 * t864;
t902 = t858 * t864;
t819 = t844 * t901 + t867 * t845 + t855 * t902;
t814 = -pkin(2) * t869 + qJDD(2) * qJ(3) + t819;
t835 = -t844 * t858 + t855 * t861;
t891 = qJD(2) * qJD(3);
t895 = t859 * t835 - 0.2e1 * t856 * t891;
t911 = pkin(3) * t859;
t784 = (-pkin(8) * qJDD(2) + t869 * t911 - t814) * t856 + t895;
t787 = t856 * t835 + (t814 + 0.2e1 * t891) * t859;
t890 = qJDD(2) * t859;
t854 = t859 ^ 2;
t903 = t854 * t869;
t785 = -pkin(3) * t903 + pkin(8) * t890 + t787;
t777 = t863 * t784 + t866 * t785;
t825 = pkin(4) * t837 - pkin(9) * t838;
t868 = qJD(4) ^ 2;
t775 = -pkin(4) * t868 + qJDD(4) * pkin(9) - t825 * t837 + t777;
t853 = t856 ^ 2;
t874 = qJDD(3) - t818;
t802 = (-pkin(2) - t911) * qJDD(2) + (-qJ(3) + (-t853 - t854) * pkin(8)) * t869 + t874;
t780 = (-t827 + t893) * pkin(9) + (-t826 + t892) * pkin(4) + t802;
t770 = -t862 * t775 + t865 * t780;
t824 = qJDD(5) - t826;
t836 = qJD(5) + t837;
t767 = -0.2e1 * qJD(6) * t830 + (t829 * t836 - t801) * qJ(6) + (t829 * t830 + t824) * pkin(5) + t770;
t809 = -mrSges(7,2) * t836 + mrSges(7,3) * t829;
t889 = m(7) * t767 + t824 * mrSges(7,1) + t836 * t809;
t764 = -t801 * mrSges(7,3) - t830 * t804 + t889;
t771 = t865 * t775 + t862 * t780;
t800 = -qJD(5) * t830 + qJDD(4) * t865 - t827 * t862;
t811 = pkin(5) * t836 - qJ(6) * t830;
t828 = t829 ^ 2;
t769 = -pkin(5) * t828 + qJ(6) * t800 + 0.2e1 * qJD(6) * t829 - t811 * t836 + t771;
t897 = t909 * t829 + t915 * t830 + t908 * t836;
t898 = -t914 * t829 - t909 * t830 - t907 * t836;
t912 = mrSges(6,1) * t770 + mrSges(7,1) * t767 - mrSges(6,2) * t771 - mrSges(7,2) * t769 + pkin(5) * t764 + t800 * t907 + t801 * t908 + t913 * t824 - t829 * t897 - t830 * t898;
t910 = -mrSges(6,2) - mrSges(7,2);
t905 = mrSges(4,2) * t856;
t808 = -qJDD(2) * pkin(2) - t869 * qJ(3) + t874;
t805 = -mrSges(6,1) * t829 + mrSges(6,2) * t830;
t810 = -mrSges(6,2) * t836 + mrSges(6,3) * t829;
t757 = m(6) * t770 + t824 * mrSges(6,1) + t836 * t810 + (-t804 - t805) * t830 + (-mrSges(6,3) - mrSges(7,3)) * t801 + t889;
t888 = m(7) * t769 + t800 * mrSges(7,3) + t829 * t804;
t812 = mrSges(7,1) * t836 - mrSges(7,3) * t830;
t896 = -mrSges(6,1) * t836 + mrSges(6,3) * t830 - t812;
t761 = m(6) * t771 + t800 * mrSges(6,3) + t829 * t805 + t824 * t910 + t836 * t896 + t888;
t754 = t865 * t757 + t862 * t761;
t833 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t837;
t834 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t838;
t873 = m(5) * t802 - t826 * mrSges(5,1) + t827 * mrSges(5,2) + t837 * t833 + t838 * t834 + t754;
t871 = -m(4) * t808 + mrSges(4,1) * t890 - t873 + (t853 * t869 + t903) * mrSges(4,3);
t748 = t871 + (mrSges(3,1) - t905) * qJDD(2) - t869 * mrSges(3,2) + m(3) * t818;
t904 = t748 * t867;
t755 = -t757 * t862 + t865 * t761;
t822 = mrSges(5,1) * t837 + mrSges(5,2) * t838;
t751 = m(5) * t777 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t826 - qJD(4) * t834 - t822 * t837 + t755;
t776 = t784 * t866 - t863 * t785;
t774 = -qJDD(4) * pkin(4) - pkin(9) * t868 + t838 * t825 - t776;
t772 = -pkin(5) * t800 - qJ(6) * t828 + t811 * t830 + qJDD(6) + t774;
t883 = -m(7) * t772 + t800 * mrSges(7,1) + t829 * t809;
t763 = -m(6) * t774 + t800 * mrSges(6,1) + t801 * t910 + t829 * t810 + t830 * t896 + t883;
t762 = m(5) * t776 + qJDD(4) * mrSges(5,1) - t827 * mrSges(5,3) + qJD(4) * t833 - t838 * t822 + t763;
t744 = t863 * t751 + t866 * t762;
t786 = -t814 * t856 + t895;
t876 = mrSges(4,3) * qJDD(2) + t869 * (-mrSges(4,1) * t859 + t905);
t742 = m(4) * t786 - t856 * t876 + t744;
t885 = t866 * t751 - t863 * t762;
t743 = m(4) * t787 + t859 * t876 + t885;
t886 = -t742 * t856 + t859 * t743;
t734 = m(3) * t819 - mrSges(3,1) * t869 - qJDD(2) * mrSges(3,2) + t886;
t737 = t859 * t742 + t856 * t743;
t736 = m(3) * t835 + t737;
t724 = t734 * t901 - t736 * t858 + t861 * t904;
t722 = m(2) * t844 + t724;
t729 = t867 * t734 - t748 * t864;
t728 = m(2) * t845 + t729;
t900 = t860 * t722 + t857 * t728;
t899 = -t907 * t829 - t908 * t830 - t913 * t836;
t880 = Ifges(4,5) * t856 + Ifges(4,6) * t859;
t894 = t869 * t880;
t723 = t734 * t902 + t861 * t736 + t858 * t904;
t887 = -t722 * t857 + t860 * t728;
t884 = m(2) * t855 + t723;
t882 = Ifges(4,1) * t856 + Ifges(4,4) * t859;
t881 = Ifges(4,4) * t856 + Ifges(4,2) * t859;
t765 = t801 * mrSges(7,2) + t830 * t812 - t883;
t745 = -mrSges(6,1) * t774 + mrSges(6,3) * t771 - mrSges(7,1) * t772 + mrSges(7,3) * t769 - pkin(5) * t765 + qJ(6) * t888 + (-qJ(6) * t812 + t897) * t836 + t899 * t830 + (-mrSges(7,2) * qJ(6) + t907) * t824 + t909 * t801 + t914 * t800;
t753 = mrSges(6,2) * t774 + mrSges(7,2) * t772 - mrSges(6,3) * t770 - mrSges(7,3) * t767 - qJ(6) * t764 + t909 * t800 + t915 * t801 + t908 * t824 - t899 * t829 + t898 * t836;
t815 = Ifges(5,5) * t838 - Ifges(5,6) * t837 + Ifges(5,3) * qJD(4);
t816 = Ifges(5,4) * t838 - Ifges(5,2) * t837 + Ifges(5,6) * qJD(4);
t730 = mrSges(5,2) * t802 - mrSges(5,3) * t776 + Ifges(5,1) * t827 + Ifges(5,4) * t826 + Ifges(5,5) * qJDD(4) - pkin(9) * t754 - qJD(4) * t816 - t745 * t862 + t753 * t865 - t815 * t837;
t817 = Ifges(5,1) * t838 - Ifges(5,4) * t837 + Ifges(5,5) * qJD(4);
t738 = -mrSges(5,1) * t802 + mrSges(5,3) * t777 + Ifges(5,4) * t827 + Ifges(5,2) * t826 + Ifges(5,6) * qJDD(4) - pkin(4) * t754 + qJD(4) * t817 - t838 * t815 - t912;
t720 = -mrSges(4,1) * t808 + mrSges(4,3) * t787 - pkin(3) * t873 + pkin(8) * t885 + qJDD(2) * t881 + t863 * t730 + t866 * t738 - t856 * t894;
t725 = mrSges(4,2) * t808 - mrSges(4,3) * t786 - pkin(8) * t744 + qJDD(2) * t882 + t866 * t730 - t863 * t738 + t859 * t894;
t718 = mrSges(3,2) * t835 - mrSges(3,3) * t818 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t869 - qJ(3) * t737 - t720 * t856 + t725 * t859;
t870 = mrSges(5,1) * t776 - mrSges(5,2) * t777 + Ifges(5,5) * t827 + Ifges(5,6) * t826 + Ifges(5,3) * qJDD(4) + pkin(4) * t763 + pkin(9) * t755 + t865 * t745 + t862 * t753 + t838 * t816 + t837 * t817;
t719 = -t870 + (Ifges(3,6) - t880) * qJDD(2) - pkin(2) * t737 - mrSges(3,1) * t835 + mrSges(3,3) * t819 - mrSges(4,1) * t786 + mrSges(4,2) * t787 - pkin(3) * t744 + (-t856 * t881 + t859 * t882 + Ifges(3,5)) * t869;
t875 = pkin(7) * t729 + t718 * t864 + t719 * t867;
t752 = qJDD(2) * t905 - t871;
t717 = mrSges(3,1) * t818 - mrSges(3,2) * t819 + Ifges(3,3) * qJDD(2) - pkin(2) * t752 + qJ(3) * t886 + t859 * t720 + t856 * t725;
t716 = mrSges(2,2) * t855 - mrSges(2,3) * t844 + t867 * t718 - t864 * t719 + (-t723 * t858 - t724 * t861) * pkin(7);
t715 = -mrSges(2,1) * t855 + mrSges(2,3) * t845 - pkin(1) * t723 - t858 * t717 + t861 * t875;
t1 = [-m(1) * g(1) + t887; -m(1) * g(2) + t900; -m(1) * g(3) + t884; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t900 - t857 * t715 + t860 * t716; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t887 + t860 * t715 + t857 * t716; -mrSges(1,1) * g(2) + mrSges(2,1) * t844 + mrSges(1,2) * g(1) - mrSges(2,2) * t845 + pkin(1) * t724 + t861 * t717 + t858 * t875; t884; t717; t752; t870; t912; t765;];
tauJB  = t1;
