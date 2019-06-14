% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-05 22:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:15:45
% EndTime: 2019-05-05 22:15:55
% DurationCPUTime: 7.31s
% Computational Cost: add. (89448->323), mult. (170335->386), div. (0->0), fcn. (106123->10), ass. (0->134)
t906 = Ifges(5,1) + Ifges(6,1);
t897 = Ifges(5,4) - Ifges(6,5);
t896 = Ifges(5,5) + Ifges(6,4);
t905 = -Ifges(5,2) - Ifges(6,3);
t895 = Ifges(5,6) - Ifges(6,6);
t904 = Ifges(5,3) + Ifges(6,2);
t860 = sin(qJ(4));
t861 = sin(qJ(3));
t887 = qJD(1) * t861;
t899 = cos(qJ(4));
t830 = -t899 * qJD(3) + t860 * t887;
t864 = cos(qJ(3));
t885 = qJD(1) * qJD(3);
t883 = t864 * t885;
t836 = t861 * qJDD(1) + t883;
t798 = -t830 * qJD(4) + t860 * qJDD(3) + t899 * t836;
t862 = sin(qJ(1));
t865 = cos(qJ(1));
t842 = t862 * g(1) - t865 * g(2);
t832 = qJDD(1) * pkin(1) + t842;
t843 = -t865 * g(1) - t862 * g(2);
t867 = qJD(1) ^ 2;
t834 = -t867 * pkin(1) + t843;
t857 = sin(pkin(10));
t858 = cos(pkin(10));
t803 = t857 * t832 + t858 * t834;
t786 = -t867 * pkin(2) + qJDD(1) * pkin(7) + t803;
t856 = -g(3) + qJDD(2);
t777 = -t861 * t786 + t864 * t856;
t835 = (-pkin(3) * t864 - pkin(8) * t861) * qJD(1);
t866 = qJD(3) ^ 2;
t874 = qJDD(3) * pkin(3) + t866 * pkin(8) - t835 * t887 + t777;
t886 = t864 * qJD(1);
t846 = qJD(4) - t886;
t893 = t830 * t846;
t903 = (-t798 + t893) * qJ(5) - t874;
t802 = t858 * t832 - t857 * t834;
t785 = -qJDD(1) * pkin(2) - t867 * pkin(7) - t802;
t884 = t861 * t885;
t837 = t864 * qJDD(1) - t884;
t769 = (-t836 - t883) * pkin(8) + (-t837 + t884) * pkin(3) + t785;
t778 = t864 * t786 + t861 * t856;
t775 = -t866 * pkin(3) + qJDD(3) * pkin(8) + t835 * t886 + t778;
t756 = t899 * t769 - t860 * t775;
t831 = t860 * qJD(3) + t899 * t887;
t806 = t830 * pkin(4) - t831 * qJ(5);
t829 = qJDD(4) - t837;
t845 = t846 ^ 2;
t754 = -t829 * pkin(4) - t845 * qJ(5) + t831 * t806 + qJDD(5) - t756;
t748 = (-t798 - t893) * pkin(9) + (t830 * t831 - t829) * pkin(5) + t754;
t757 = t860 * t769 + t899 * t775;
t900 = 2 * qJD(5);
t753 = -t845 * pkin(4) + t829 * qJ(5) - t830 * t806 + t846 * t900 + t757;
t797 = t831 * qJD(4) - t899 * qJDD(3) + t860 * t836;
t813 = -t846 * pkin(5) - t831 * pkin(9);
t828 = t830 ^ 2;
t749 = -t828 * pkin(5) + t797 * pkin(9) + t846 * t813 + t753;
t859 = sin(qJ(6));
t863 = cos(qJ(6));
t747 = t859 * t748 + t863 * t749;
t751 = -t828 * pkin(9) + (-pkin(4) - pkin(5)) * t797 + (-pkin(4) * t846 + t813 + t900) * t831 - t903;
t801 = t859 * t830 + t863 * t831;
t763 = -t801 * qJD(6) + t863 * t797 - t859 * t798;
t800 = t863 * t830 - t859 * t831;
t764 = t800 * qJD(6) + t859 * t797 + t863 * t798;
t844 = qJD(6) - t846;
t766 = Ifges(7,5) * t801 + Ifges(7,6) * t800 + Ifges(7,3) * t844;
t768 = Ifges(7,1) * t801 + Ifges(7,4) * t800 + Ifges(7,5) * t844;
t825 = qJDD(6) - t829;
t736 = -mrSges(7,1) * t751 + mrSges(7,3) * t747 + Ifges(7,4) * t764 + Ifges(7,2) * t763 + Ifges(7,6) * t825 - t801 * t766 + t844 * t768;
t746 = t863 * t748 - t859 * t749;
t767 = Ifges(7,4) * t801 + Ifges(7,2) * t800 + Ifges(7,6) * t844;
t737 = mrSges(7,2) * t751 - mrSges(7,3) * t746 + Ifges(7,1) * t764 + Ifges(7,4) * t763 + Ifges(7,5) * t825 + t800 * t766 - t844 * t767;
t755 = -0.2e1 * qJD(5) * t831 + (t831 * t846 + t797) * pkin(4) + t903;
t811 = -t846 * mrSges(6,1) + t831 * mrSges(6,2);
t812 = -t830 * mrSges(6,2) + t846 * mrSges(6,3);
t779 = -t844 * mrSges(7,2) + t800 * mrSges(7,3);
t780 = t844 * mrSges(7,1) - t801 * mrSges(7,3);
t877 = -m(7) * t751 + t763 * mrSges(7,1) - t764 * mrSges(7,2) + t800 * t779 - t801 * t780;
t741 = m(6) * t755 + t797 * mrSges(6,1) - t798 * mrSges(6,3) - t831 * t811 + t830 * t812 + t877;
t776 = -t800 * mrSges(7,1) + t801 * mrSges(7,2);
t743 = m(7) * t746 + t825 * mrSges(7,1) - t764 * mrSges(7,3) - t801 * t776 + t844 * t779;
t744 = m(7) * t747 - t825 * mrSges(7,2) + t763 * mrSges(7,3) + t800 * t776 - t844 * t780;
t879 = -t859 * t743 + t863 * t744;
t889 = -t897 * t830 + t906 * t831 + t896 * t846;
t891 = t895 * t830 - t896 * t831 - t904 * t846;
t714 = mrSges(5,1) * t874 - mrSges(6,1) * t755 + mrSges(6,2) * t753 + mrSges(5,3) * t757 - pkin(4) * t741 - pkin(5) * t877 - pkin(9) * t879 - t863 * t736 - t859 * t737 + t905 * t797 + t897 * t798 + t895 * t829 + t891 * t831 + t889 * t846;
t735 = t863 * t743 + t859 * t744;
t890 = t905 * t830 + t897 * t831 + t895 * t846;
t715 = -mrSges(5,2) * t874 + mrSges(6,2) * t754 - mrSges(5,3) * t756 - mrSges(6,3) * t755 - pkin(9) * t735 - qJ(5) * t741 - t859 * t736 + t863 * t737 - t897 * t797 + t906 * t798 + t896 * t829 + t891 * t830 - t890 * t846;
t810 = t846 * mrSges(5,1) - t831 * mrSges(5,3);
t875 = m(6) * t753 + t829 * mrSges(6,3) + t846 * t811 + t879;
t807 = t830 * mrSges(6,1) - t831 * mrSges(6,3);
t888 = -t830 * mrSges(5,1) - t831 * mrSges(5,2) - t807;
t898 = -mrSges(5,3) - mrSges(6,2);
t731 = m(5) * t757 - t829 * mrSges(5,2) + t898 * t797 - t846 * t810 + t888 * t830 + t875;
t809 = -t846 * mrSges(5,2) - t830 * mrSges(5,3);
t873 = -m(6) * t754 + t829 * mrSges(6,1) + t846 * t812 - t735;
t732 = m(5) * t756 + t829 * mrSges(5,1) + t898 * t798 + t846 * t809 + t888 * t831 + t873;
t729 = t899 * t731 - t860 * t732;
t740 = m(5) * t874 - t797 * mrSges(5,1) - t798 * mrSges(5,2) - t830 * t809 - t831 * t810 - t741;
t820 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t861 + Ifges(4,2) * t864) * qJD(1);
t821 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t861 + Ifges(4,4) * t864) * qJD(1);
t902 = mrSges(4,1) * t777 - mrSges(4,2) * t778 + Ifges(4,5) * t836 + Ifges(4,6) * t837 + Ifges(4,3) * qJDD(3) + pkin(3) * t740 + pkin(8) * t729 + (t820 * t861 - t821 * t864) * qJD(1) + t899 * t714 + t860 * t715;
t734 = t798 * mrSges(6,2) + t831 * t807 - t873;
t872 = -mrSges(7,1) * t746 + mrSges(7,2) * t747 - Ifges(7,5) * t764 - Ifges(7,6) * t763 - Ifges(7,3) * t825 - t801 * t767 + t800 * t768;
t901 = -t895 * t797 + t896 * t798 + t904 * t829 + t889 * t830 + t890 * t831 + mrSges(5,1) * t756 - mrSges(6,1) * t754 - mrSges(5,2) * t757 + mrSges(6,3) * t753 - pkin(4) * t734 - pkin(5) * t735 + qJ(5) * (-t797 * mrSges(6,2) - t830 * t807 + t875) + t872;
t833 = (-mrSges(4,1) * t864 + mrSges(4,2) * t861) * qJD(1);
t839 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t887;
t727 = m(4) * t778 - qJDD(3) * mrSges(4,2) + t837 * mrSges(4,3) - qJD(3) * t839 + t833 * t886 + t729;
t840 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t886;
t739 = m(4) * t777 + qJDD(3) * mrSges(4,1) - t836 * mrSges(4,3) + qJD(3) * t840 - t833 * t887 + t740;
t880 = t864 * t727 - t861 * t739;
t718 = m(3) * t803 - t867 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t880;
t728 = t860 * t731 + t899 * t732;
t870 = -m(4) * t785 + t837 * mrSges(4,1) - t836 * mrSges(4,2) - t839 * t887 + t840 * t886 - t728;
t723 = m(3) * t802 + qJDD(1) * mrSges(3,1) - t867 * mrSges(3,2) + t870;
t713 = t857 * t718 + t858 * t723;
t710 = m(2) * t842 + qJDD(1) * mrSges(2,1) - t867 * mrSges(2,2) + t713;
t881 = t858 * t718 - t857 * t723;
t711 = m(2) * t843 - t867 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t881;
t892 = t865 * t710 + t862 * t711;
t721 = t861 * t727 + t864 * t739;
t719 = m(3) * t856 + t721;
t882 = -t862 * t710 + t865 * t711;
t819 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t861 + Ifges(4,6) * t864) * qJD(1);
t704 = mrSges(4,2) * t785 - mrSges(4,3) * t777 + Ifges(4,1) * t836 + Ifges(4,4) * t837 + Ifges(4,5) * qJDD(3) - pkin(8) * t728 - qJD(3) * t820 - t860 * t714 + t899 * t715 + t819 * t886;
t706 = -mrSges(4,1) * t785 + mrSges(4,3) * t778 + Ifges(4,4) * t836 + Ifges(4,2) * t837 + Ifges(4,6) * qJDD(3) - pkin(3) * t728 + qJD(3) * t821 - t819 * t887 - t901;
t871 = mrSges(2,1) * t842 + mrSges(3,1) * t802 - mrSges(2,2) * t843 - mrSges(3,2) * t803 + pkin(1) * t713 + pkin(2) * t870 + pkin(7) * t880 + t861 * t704 + t864 * t706 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t702 = -mrSges(3,1) * t856 + mrSges(3,3) * t803 + t867 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t721 - t902;
t701 = mrSges(3,2) * t856 - mrSges(3,3) * t802 + Ifges(3,5) * qJDD(1) - t867 * Ifges(3,6) - pkin(7) * t721 + t864 * t704 - t861 * t706;
t700 = -mrSges(2,2) * g(3) - mrSges(2,3) * t842 + Ifges(2,5) * qJDD(1) - t867 * Ifges(2,6) - qJ(2) * t713 + t858 * t701 - t857 * t702;
t699 = mrSges(2,1) * g(3) + mrSges(2,3) * t843 + t867 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t719 + qJ(2) * t881 + t857 * t701 + t858 * t702;
t1 = [-m(1) * g(1) + t882; -m(1) * g(2) + t892; (-m(1) - m(2)) * g(3) + t719; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t892 - t862 * t699 + t865 * t700; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t882 + t865 * t699 + t862 * t700; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t871; t871; t719; t902; t901; t734; -t872;];
tauJB  = t1;
