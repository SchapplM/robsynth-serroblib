% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:14:41
% EndTime: 2019-05-06 09:14:45
% DurationCPUTime: 2.80s
% Computational Cost: add. (19201->332), mult. (39621->372), div. (0->0), fcn. (19499->6), ass. (0->132)
t881 = Ifges(6,1) + Ifges(7,1);
t862 = Ifges(6,4) + Ifges(7,4);
t861 = Ifges(6,5) + Ifges(7,5);
t880 = -Ifges(6,2) - Ifges(7,2);
t860 = Ifges(6,6) + Ifges(7,6);
t879 = Ifges(6,3) + Ifges(7,3);
t811 = sin(qJ(5));
t814 = cos(qJ(5));
t815 = cos(qJ(2));
t848 = qJD(1) * t815;
t771 = qJD(2) * t811 + t814 * t848;
t812 = sin(qJ(2));
t846 = qJD(1) * qJD(2);
t832 = t812 * t846;
t778 = qJDD(1) * t815 - t832;
t724 = qJD(5) * t771 - qJDD(2) * t814 + t778 * t811;
t872 = (mrSges(6,1) + mrSges(7,1)) * t724;
t770 = -qJD(2) * t814 + t811 * t848;
t849 = qJD(1) * t812;
t794 = qJD(5) + t849;
t729 = -mrSges(7,2) * t794 + mrSges(7,3) * t770;
t730 = -mrSges(6,2) * t794 + mrSges(6,3) * t770;
t878 = t872 + (t729 + t730) * t770;
t877 = Ifges(3,1) + Ifges(4,1) + Ifges(5,2);
t843 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t842 = Ifges(3,5) + Ifges(4,4) + Ifges(5,6);
t876 = -Ifges(4,2) - Ifges(3,3) - Ifges(5,3);
t841 = Ifges(3,6) - Ifges(4,6) + Ifges(5,5);
t875 = Ifges(4,3) + Ifges(5,1) + Ifges(3,2);
t874 = qJ(3) + pkin(4);
t783 = -qJD(2) * pkin(3) - qJ(4) * t849;
t845 = qJD(1) * qJD(4);
t818 = qJD(1) ^ 2;
t857 = t815 ^ 2 * t818;
t813 = sin(qJ(1));
t816 = cos(qJ(1));
t791 = -g(1) * t816 - g(2) * t813;
t752 = -pkin(1) * t818 + qJDD(1) * pkin(7) + t791;
t735 = -t812 * g(3) + t815 * t752;
t772 = (-pkin(2) * t815 - qJ(3) * t812) * qJD(1);
t870 = qJDD(2) * qJ(3) + t772 * t848 + t735;
t871 = pkin(3) * t857 + t778 * qJ(4) - qJD(2) * t783 + 0.2e1 * t815 * t845 - t870;
t834 = t815 * t846;
t777 = qJDD(1) * t812 + t834;
t869 = -0.2e1 * t812 * t845 + (-t777 + t834) * qJ(4);
t868 = 2 * qJD(6);
t867 = -pkin(2) - pkin(8);
t866 = pkin(3) + pkin(8);
t817 = qJD(2) ^ 2;
t865 = pkin(2) * t817;
t863 = mrSges(3,3) + mrSges(4,2);
t859 = t770 * t729;
t856 = t815 * t818;
t734 = -t815 * g(3) - t812 * t752;
t773 = (-mrSges(4,1) * t815 - mrSges(4,3) * t812) * qJD(1);
t774 = (-mrSges(3,1) * t815 + mrSges(3,2) * t812) * qJD(1);
t787 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t848;
t788 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t848;
t829 = t772 * t849 + qJDD(3) - t734;
t711 = -qJDD(2) * pkin(2) - t817 * qJ(3) + t829;
t789 = mrSges(4,2) * t848 + qJD(2) * mrSges(4,3);
t707 = (-t812 * t856 - qJDD(2)) * pkin(3) + t711 + t869;
t775 = (mrSges(5,1) * t812 - mrSges(5,2) * t815) * qJD(1);
t790 = t813 * g(1) - t816 * g(2);
t751 = -qJDD(1) * pkin(1) - t818 * pkin(7) - t790;
t825 = -t778 * pkin(2) + t751 + (-t777 - t834) * qJ(3);
t821 = -qJ(4) * t857 + qJDD(4) - t825 + ((2 * qJD(3)) + t783) * t849;
t699 = t821 + t866 * t778 + (pkin(4) * t815 + t867 * t812) * t846 + t777 * pkin(4);
t776 = (pkin(4) * t812 + pkin(8) * t815) * qJD(1);
t703 = -t874 * t817 + (-pkin(3) * t856 - qJD(1) * t776) * t812 + (-pkin(2) - t866) * qJDD(2) + t829 + t869;
t693 = t814 * t699 - t811 * t703;
t725 = qJD(5) * t770 - qJDD(2) * t811 - t778 * t814;
t727 = -mrSges(7,1) * t770 - mrSges(7,2) * t771;
t728 = -mrSges(6,1) * t770 - mrSges(6,2) * t771;
t768 = qJDD(5) + t777;
t690 = t771 * t868 + (t770 * t794 - t725) * qJ(6) + (-t770 * t771 + t768) * pkin(5) + t693;
t840 = m(7) * t690 + t768 * mrSges(7,1) + t794 * t729;
t684 = m(6) * t693 + t768 * mrSges(6,1) + t794 * t730 + (t727 + t728) * t771 + (-mrSges(6,3) - mrSges(7,3)) * t725 + t840;
t694 = t811 * t699 + t814 * t703;
t732 = mrSges(7,1) * t794 + mrSges(7,3) * t771;
t733 = mrSges(6,1) * t794 + mrSges(6,3) * t771;
t731 = pkin(5) * t794 + qJ(6) * t771;
t767 = t770 ^ 2;
t692 = -pkin(5) * t767 + qJ(6) * t724 - t731 * t794 + t770 * t868 + t694;
t839 = m(7) * t692 + t724 * mrSges(7,3) + t770 * t727;
t687 = m(6) * t694 + t724 * mrSges(6,3) + t770 * t728 + (-t732 - t733) * t794 + (-mrSges(6,2) - mrSges(7,2)) * t768 + t839;
t854 = -t811 * t684 + t814 * t687;
t828 = -m(5) * t707 + t775 * t849 - t854;
t824 = -m(4) * t711 + qJDD(2) * mrSges(4,1) + qJD(2) * t789 + t828;
t674 = m(3) * t734 + (mrSges(3,1) - mrSges(5,2)) * qJDD(2) + (-t787 + t788) * qJD(2) + (-t773 - t774) * t849 + (mrSges(5,3) - t863) * t777 + t824;
t785 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t849;
t844 = qJD(3) * qJD(2);
t799 = 0.2e1 * t844;
t710 = t799 - t865 + t870;
t786 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t849;
t706 = -0.2e1 * t844 + t865 + t871;
t784 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t849;
t702 = qJDD(2) * pkin(4) - t776 * t848 + t867 * t817 + t799 - t871;
t696 = -t724 * pkin(5) - t767 * qJ(6) - t771 * t731 + qJDD(6) + t702;
t838 = -m(7) * t696 - t725 * mrSges(7,2) + t771 * t732;
t827 = -m(6) * t702 - t725 * mrSges(6,2) + t771 * t733 + t838;
t822 = -m(5) * t706 + qJDD(2) * mrSges(5,1) - t778 * mrSges(5,3) + qJD(2) * t784 - t827;
t820 = m(4) * t710 + qJDD(2) * mrSges(4,3) + qJD(2) * t786 + t773 * t848 + t822;
t681 = t820 - qJDD(2) * mrSges(3,2) + t863 * t778 - qJD(2) * t785 + m(3) * t735 + (t774 - t775) * t848 - t878;
t830 = -t674 * t812 + t815 * t681;
t668 = m(2) * t791 - mrSges(2,1) * t818 - qJDD(1) * mrSges(2,2) + t830;
t708 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t849 + t825;
t678 = t814 * t684 + t811 * t687;
t705 = -pkin(2) * t832 + t778 * pkin(3) + t821;
t826 = -m(5) * t705 - t777 * mrSges(5,1) + t778 * mrSges(5,2) - t784 * t849 + t787 * t848 - t678;
t675 = m(4) * t708 - mrSges(4,1) * t778 - t777 * mrSges(4,3) - t786 * t849 - t789 * t848 + t826;
t819 = -m(3) * t751 + t778 * mrSges(3,1) - mrSges(3,2) * t777 - t785 * t849 + t788 * t848 - t675;
t672 = m(2) * t790 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t818 + t819;
t855 = t813 * t668 + t816 * t672;
t669 = t815 * t674 + t812 * t681;
t853 = t860 * t770 - t861 * t771 + t879 * t794;
t852 = t880 * t770 + t862 * t771 - t860 * t794;
t851 = -t862 * t770 + t881 * t771 - t861 * t794;
t847 = qJD(2) * t787;
t837 = t876 * qJD(2) + (-t842 * t812 - t841 * t815) * qJD(1);
t836 = -t841 * qJD(2) + (-t843 * t812 - t875 * t815) * qJD(1);
t835 = t842 * qJD(2) + (t877 * t812 + t843 * t815) * qJD(1);
t831 = t816 * t668 - t672 * t813;
t688 = -t725 * mrSges(7,3) + t771 * t727 + t840;
t677 = mrSges(6,2) * t702 + mrSges(7,2) * t696 - mrSges(6,3) * t693 - mrSges(7,3) * t690 - qJ(6) * t688 + t862 * t724 + t881 * t725 + t861 * t768 + t853 * t770 + t852 * t794;
t676 = qJDD(2) * mrSges(5,2) - t777 * mrSges(5,3) - t828 + t847;
t670 = -mrSges(6,1) * t702 + mrSges(6,3) * t694 - mrSges(7,1) * t696 + mrSges(7,3) * t692 - pkin(5) * (-t838 - t859) + qJ(6) * t839 + (-qJ(6) * t732 - t851) * t794 + t853 * t771 + (-mrSges(7,2) * qJ(6) + t860) * t768 + t862 * t725 + (mrSges(7,1) * pkin(5) - t880) * t724;
t665 = mrSges(5,1) * t705 + mrSges(6,1) * t693 + mrSges(7,1) * t690 + mrSges(3,2) * t751 + mrSges(4,2) * t711 - mrSges(6,2) * t694 - mrSges(7,2) * t692 - mrSges(3,3) * t734 - mrSges(4,3) * t708 - mrSges(5,3) * t707 + pkin(4) * t678 + pkin(5) * t688 - qJ(3) * t675 - qJ(4) * t676 + t852 * t771 + t851 * t770 + t879 * t768 + t861 * t725 + t860 * t724 + t843 * t778 + t877 * t777 + t842 * qJDD(2) + t836 * qJD(2) - t837 * t848;
t664 = -qJ(4) * (-t770 * t730 + t822 - t859 - t872) + t811 * t670 - t814 * t677 - pkin(3) * t826 + mrSges(3,3) * t735 - mrSges(3,1) * t751 + mrSges(4,2) * t710 - mrSges(5,2) * t705 + mrSges(5,3) * t706 - mrSges(4,1) * t708 + pkin(8) * t678 - pkin(2) * t675 + t875 * t778 + t843 * t777 + t841 * qJDD(2) + t835 * qJD(2) + (qJ(4) * t815 * t775 + t837 * t812) * qJD(1);
t663 = Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + t818 * Ifges(2,5) + t811 * t677 + t814 * t670 + mrSges(2,3) * t791 - qJ(3) * t820 - pkin(2) * (t824 - t847) + pkin(4) * t827 - mrSges(3,1) * t734 + mrSges(3,2) * t735 - mrSges(4,3) * t710 + mrSges(4,1) * t711 + mrSges(5,1) * t706 - mrSges(5,2) * t707 + pkin(8) * t854 + pkin(3) * t676 - pkin(1) * t669 + (-mrSges(4,2) * qJ(3) - t841) * t778 + (mrSges(5,2) * pkin(2) + t876) * qJDD(2) + (-pkin(2) * (-mrSges(4,2) + mrSges(5,3)) - t842) * t777 + ((qJ(3) * t775 + t835) * t815 + (pkin(2) * t773 + t836) * t812) * qJD(1) + t878 * t874;
t662 = -mrSges(2,2) * g(3) - mrSges(2,3) * t790 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t818 - pkin(7) * t669 - t664 * t812 + t665 * t815;
t1 = [-m(1) * g(1) + t831; -m(1) * g(2) + t855; (-m(1) - m(2)) * g(3) + t669; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t855 + t816 * t662 - t813 * t663; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t831 + t813 * t662 + t816 * t663; -mrSges(1,1) * g(2) + mrSges(2,1) * t790 + mrSges(1,2) * g(1) - mrSges(2,2) * t791 + Ifges(2,3) * qJDD(1) + pkin(1) * t819 + pkin(7) * t830 + t815 * t664 + t812 * t665;];
tauB  = t1;
