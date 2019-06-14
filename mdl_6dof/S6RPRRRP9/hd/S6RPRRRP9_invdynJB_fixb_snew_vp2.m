% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-05-06 01:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRP9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:51:25
% EndTime: 2019-05-06 01:51:35
% DurationCPUTime: 6.75s
% Computational Cost: add. (77421->316), mult. (151018->369), div. (0->0), fcn. (97237->8), ass. (0->128)
t886 = Ifges(6,4) + Ifges(7,4);
t898 = Ifges(6,2) + Ifges(7,2);
t894 = Ifges(6,6) + Ifges(7,6);
t851 = sin(qJ(4));
t855 = cos(qJ(4));
t856 = cos(qJ(3));
t881 = qJD(1) * t856;
t826 = t851 * qJD(3) + t855 * t881;
t852 = sin(qJ(3));
t880 = qJD(1) * qJD(3);
t876 = t852 * t880;
t830 = t856 * qJDD(1) - t876;
t789 = -t826 * qJD(4) + t855 * qJDD(3) - t851 * t830;
t825 = t855 * qJD(3) - t851 * t881;
t790 = t825 * qJD(4) + t851 * qJDD(3) + t855 * t830;
t850 = sin(qJ(5));
t854 = cos(qJ(5));
t792 = t854 * t825 - t850 * t826;
t758 = t792 * qJD(5) + t850 * t789 + t854 * t790;
t793 = t850 * t825 + t854 * t826;
t769 = -t792 * mrSges(7,1) + t793 * mrSges(7,2);
t859 = qJD(1) ^ 2;
t889 = -pkin(1) - pkin(7);
t853 = sin(qJ(1));
t857 = cos(qJ(1));
t835 = -t857 * g(1) - t853 * g(2);
t890 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t835;
t803 = t889 * t859 - t890;
t838 = t856 * t880;
t829 = -t852 * qJDD(1) - t838;
t773 = (-t830 + t876) * pkin(8) + (-t829 + t838) * pkin(3) + t803;
t834 = t853 * g(1) - t857 * g(2);
t869 = -t859 * qJ(2) + qJDD(2) - t834;
t804 = t889 * qJDD(1) + t869;
t797 = -t856 * g(3) + t852 * t804;
t828 = (pkin(3) * t852 - pkin(8) * t856) * qJD(1);
t840 = t852 * qJD(1);
t858 = qJD(3) ^ 2;
t777 = -t858 * pkin(3) + qJDD(3) * pkin(8) - t828 * t840 + t797;
t745 = t855 * t773 - t851 * t777;
t824 = qJDD(4) - t829;
t837 = t840 + qJD(4);
t741 = (t825 * t837 - t790) * pkin(9) + (t825 * t826 + t824) * pkin(4) + t745;
t746 = t851 * t773 + t855 * t777;
t802 = t837 * pkin(4) - t826 * pkin(9);
t823 = t825 ^ 2;
t743 = -t823 * pkin(4) + t789 * pkin(9) - t837 * t802 + t746;
t735 = t854 * t741 - t850 * t743;
t817 = qJDD(5) + t824;
t836 = qJD(5) + t837;
t731 = -0.2e1 * qJD(6) * t793 + (t792 * t836 - t758) * qJ(6) + (t792 * t793 + t817) * pkin(5) + t735;
t778 = -t836 * mrSges(7,2) + t792 * mrSges(7,3);
t878 = m(7) * t731 + t817 * mrSges(7,1) + t836 * t778;
t727 = -t758 * mrSges(7,3) - t793 * t769 + t878;
t736 = t850 * t741 + t854 * t743;
t757 = -t793 * qJD(5) + t854 * t789 - t850 * t790;
t780 = t836 * pkin(5) - t793 * qJ(6);
t791 = t792 ^ 2;
t733 = -t791 * pkin(5) + t757 * qJ(6) + 0.2e1 * qJD(6) * t792 - t836 * t780 + t736;
t895 = Ifges(6,5) + Ifges(7,5);
t896 = Ifges(6,1) + Ifges(7,1);
t882 = -t886 * t792 - t896 * t793 - t895 * t836;
t892 = t898 * t792 + t886 * t793 + t894 * t836;
t893 = Ifges(6,3) + Ifges(7,3);
t897 = mrSges(6,1) * t735 + mrSges(7,1) * t731 - mrSges(6,2) * t736 - mrSges(7,2) * t733 + pkin(5) * t727 + t894 * t757 + t895 * t758 + t882 * t792 + t892 * t793 + t893 * t817;
t770 = -t792 * mrSges(6,1) + t793 * mrSges(6,2);
t779 = -t836 * mrSges(6,2) + t792 * mrSges(6,3);
t720 = m(6) * t735 + t817 * mrSges(6,1) + t836 * t779 + (-t769 - t770) * t793 + (-mrSges(6,3) - mrSges(7,3)) * t758 + t878;
t781 = t836 * mrSges(7,1) - t793 * mrSges(7,3);
t782 = t836 * mrSges(6,1) - t793 * mrSges(6,3);
t877 = m(7) * t733 + t757 * mrSges(7,3) + t792 * t769;
t723 = m(6) * t736 + t757 * mrSges(6,3) + t792 * t770 + (-t781 - t782) * t836 + (-mrSges(6,2) - mrSges(7,2)) * t817 + t877;
t718 = t854 * t720 + t850 * t723;
t784 = Ifges(5,4) * t826 + Ifges(5,2) * t825 + Ifges(5,6) * t837;
t785 = Ifges(5,1) * t826 + Ifges(5,4) * t825 + Ifges(5,5) * t837;
t891 = mrSges(5,1) * t745 - mrSges(5,2) * t746 + Ifges(5,5) * t790 + Ifges(5,6) * t789 + Ifges(5,3) * t824 + pkin(4) * t718 + t826 * t784 - t825 * t785 + t897;
t888 = mrSges(2,1) - mrSges(3,2);
t887 = -Ifges(3,4) + Ifges(2,5);
t885 = Ifges(3,5) - Ifges(2,6);
t827 = (mrSges(4,1) * t852 + mrSges(4,2) * t856) * qJD(1);
t833 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t881;
t795 = -t825 * mrSges(5,1) + t826 * mrSges(5,2);
t798 = -t837 * mrSges(5,2) + t825 * mrSges(5,3);
t715 = m(5) * t745 + t824 * mrSges(5,1) - t790 * mrSges(5,3) - t826 * t795 + t837 * t798 + t718;
t799 = t837 * mrSges(5,1) - t826 * mrSges(5,3);
t871 = -t850 * t720 + t854 * t723;
t716 = m(5) * t746 - t824 * mrSges(5,2) + t789 * mrSges(5,3) + t825 * t795 - t837 * t799 + t871;
t872 = -t851 * t715 + t855 * t716;
t708 = m(4) * t797 - qJDD(3) * mrSges(4,2) + t829 * mrSges(4,3) - qJD(3) * t833 - t827 * t840 + t872;
t796 = t852 * g(3) + t856 * t804;
t832 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t840;
t776 = -qJDD(3) * pkin(3) - t858 * pkin(8) + t828 * t881 - t796;
t744 = -t789 * pkin(4) - t823 * pkin(9) + t826 * t802 + t776;
t738 = -t757 * pkin(5) - t791 * qJ(6) + t793 * t780 + qJDD(6) + t744;
t728 = m(7) * t738 - t757 * mrSges(7,1) + t758 * mrSges(7,2) - t792 * t778 + t793 * t781;
t864 = m(6) * t744 - t757 * mrSges(6,1) + t758 * mrSges(6,2) - t792 * t779 + t793 * t782 + t728;
t861 = -m(5) * t776 + t789 * mrSges(5,1) - t790 * mrSges(5,2) + t825 * t798 - t826 * t799 - t864;
t724 = m(4) * t796 + qJDD(3) * mrSges(4,1) - t830 * mrSges(4,3) + qJD(3) * t832 - t827 * t881 + t861;
t702 = t852 * t708 + t856 * t724;
t809 = -qJDD(1) * pkin(1) + t869;
t868 = -m(3) * t809 + t859 * mrSges(3,3) - t702;
t698 = m(2) * t834 - t859 * mrSges(2,2) + t888 * qJDD(1) + t868;
t807 = t859 * pkin(1) + t890;
t710 = t855 * t715 + t851 * t716;
t867 = -m(4) * t803 + t829 * mrSges(4,1) - t830 * mrSges(4,2) - t832 * t840 - t833 * t881 - t710;
t865 = -m(3) * t807 + t859 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t867;
t705 = m(2) * t835 - t859 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t865;
t884 = t857 * t698 + t853 * t705;
t883 = -t894 * t792 - t895 * t793 - t893 * t836;
t874 = -t853 * t698 + t857 * t705;
t873 = t856 * t708 - t852 * t724;
t711 = -mrSges(6,1) * t744 + mrSges(6,3) * t736 - mrSges(7,1) * t738 + mrSges(7,3) * t733 - pkin(5) * t728 + qJ(6) * t877 + (-qJ(6) * t781 - t882) * t836 + (-qJ(6) * mrSges(7,2) + t894) * t817 + t883 * t793 + t886 * t758 + t898 * t757;
t717 = mrSges(6,2) * t744 + mrSges(7,2) * t738 - mrSges(6,3) * t735 - mrSges(7,3) * t731 - qJ(6) * t727 + t886 * t757 + t896 * t758 - t883 * t792 + t895 * t817 - t892 * t836;
t783 = Ifges(5,5) * t826 + Ifges(5,6) * t825 + Ifges(5,3) * t837;
t694 = -mrSges(5,1) * t776 + mrSges(5,3) * t746 + Ifges(5,4) * t790 + Ifges(5,2) * t789 + Ifges(5,6) * t824 - pkin(4) * t864 + pkin(9) * t871 + t854 * t711 + t850 * t717 - t826 * t783 + t837 * t785;
t696 = mrSges(5,2) * t776 - mrSges(5,3) * t745 + Ifges(5,1) * t790 + Ifges(5,4) * t789 + Ifges(5,5) * t824 - pkin(9) * t718 - t850 * t711 + t854 * t717 + t825 * t783 - t837 * t784;
t815 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t856 - Ifges(4,2) * t852) * qJD(1);
t816 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t856 - Ifges(4,4) * t852) * qJD(1);
t866 = mrSges(4,1) * t796 - mrSges(4,2) * t797 + Ifges(4,5) * t830 + Ifges(4,6) * t829 + Ifges(4,3) * qJDD(3) + pkin(3) * t861 + pkin(8) * t872 + t855 * t694 + t851 * t696 + t815 * t881 + t816 * t840;
t814 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t856 - Ifges(4,6) * t852) * qJD(1);
t691 = mrSges(4,2) * t803 - mrSges(4,3) * t796 + Ifges(4,1) * t830 + Ifges(4,4) * t829 + Ifges(4,5) * qJDD(3) - pkin(8) * t710 - qJD(3) * t815 - t851 * t694 + t855 * t696 - t814 * t840;
t692 = -mrSges(4,1) * t803 + mrSges(4,3) * t797 + Ifges(4,4) * t830 + Ifges(4,2) * t829 + Ifges(4,6) * qJDD(3) - pkin(3) * t710 + qJD(3) * t816 - t814 * t881 - t891;
t700 = qJDD(1) * mrSges(3,2) - t868;
t862 = mrSges(2,1) * t834 - mrSges(2,2) * t835 + mrSges(3,2) * t809 - mrSges(3,3) * t807 - pkin(1) * t700 - pkin(7) * t702 + qJ(2) * t865 + t856 * t691 - t852 * t692 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t701 = -m(3) * g(3) + t873;
t689 = t885 * t859 + t866 - mrSges(2,3) * t834 + mrSges(3,1) * t809 + t887 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - qJ(2) * t701 + pkin(2) * t702;
t688 = -mrSges(3,1) * t807 + mrSges(2,3) * t835 - pkin(1) * t701 - pkin(2) * t867 - pkin(7) * t873 + t888 * g(3) - t885 * qJDD(1) - t852 * t691 - t856 * t692 + t887 * t859;
t1 = [-m(1) * g(1) + t874; -m(1) * g(2) + t884; (-m(1) - m(2) - m(3)) * g(3) + t873; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t884 - t853 * t688 + t857 * t689; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t874 + t857 * t688 + t853 * t689; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t862; t862; t700; t866; t891; t897; t728;];
tauJB  = t1;
