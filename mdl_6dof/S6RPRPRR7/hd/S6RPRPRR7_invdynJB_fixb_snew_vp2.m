% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 19:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:12:57
% EndTime: 2019-05-05 19:13:06
% DurationCPUTime: 9.33s
% Computational Cost: add. (144547->338), mult. (321463->417), div. (0->0), fcn. (223952->10), ass. (0->135)
t844 = sin(qJ(1));
t848 = cos(qJ(1));
t822 = -t848 * g(1) - t844 * g(2);
t860 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t822;
t876 = -pkin(1) - pkin(7);
t875 = mrSges(2,1) - mrSges(3,2);
t874 = Ifges(2,5) - Ifges(3,4);
t873 = (-Ifges(2,6) + Ifges(3,5));
t821 = t844 * g(1) - t848 * g(2);
t849 = qJD(1) ^ 2;
t859 = -t849 * qJ(2) + qJDD(2) - t821;
t794 = qJDD(1) * t876 + t859;
t843 = sin(qJ(3));
t847 = cos(qJ(3));
t785 = t843 * g(3) + t847 * t794;
t869 = qJD(1) * qJD(3);
t867 = t843 * t869;
t816 = qJDD(1) * t847 - t867;
t761 = (-t816 - t867) * qJ(4) + (-t843 * t847 * t849 + qJDD(3)) * pkin(3) + t785;
t786 = -g(3) * t847 + t843 * t794;
t815 = -qJDD(1) * t843 - t847 * t869;
t870 = qJD(1) * t847;
t819 = qJD(3) * pkin(3) - qJ(4) * t870;
t836 = t843 ^ 2;
t762 = -pkin(3) * t836 * t849 + qJ(4) * t815 - qJD(3) * t819 + t786;
t839 = sin(pkin(10));
t840 = cos(pkin(10));
t804 = (-t839 * t843 + t840 * t847) * qJD(1);
t738 = -0.2e1 * qJD(4) * t804 + t840 * t761 - t839 * t762;
t783 = t815 * t839 + t816 * t840;
t803 = (-t839 * t847 - t840 * t843) * qJD(1);
t727 = (qJD(3) * t803 - t783) * pkin(8) + (t803 * t804 + qJDD(3)) * pkin(4) + t738;
t739 = 0.2e1 * qJD(4) * t803 + t839 * t761 + t840 * t762;
t782 = t815 * t840 - t816 * t839;
t793 = qJD(3) * pkin(4) - pkin(8) * t804;
t802 = t803 ^ 2;
t729 = -pkin(4) * t802 + pkin(8) * t782 - qJD(3) * t793 + t739;
t842 = sin(qJ(5));
t846 = cos(qJ(5));
t724 = t842 * t727 + t846 * t729;
t775 = t803 * t842 + t804 * t846;
t747 = -qJD(5) * t775 + t782 * t846 - t783 * t842;
t774 = t803 * t846 - t804 * t842;
t756 = -mrSges(6,1) * t774 + mrSges(6,2) * t775;
t829 = qJD(3) + qJD(5);
t769 = mrSges(6,1) * t829 - mrSges(6,3) * t775;
t828 = qJDD(3) + qJDD(5);
t757 = -pkin(5) * t774 - pkin(9) * t775;
t827 = t829 ^ 2;
t721 = -pkin(5) * t827 + pkin(9) * t828 + t757 * t774 + t724;
t764 = -t815 * pkin(3) + qJDD(4) + t819 * t870 + (-qJ(4) * t836 + t876) * t849 + t860;
t741 = -t782 * pkin(4) - t802 * pkin(8) + t804 * t793 + t764;
t748 = qJD(5) * t774 + t782 * t842 + t783 * t846;
t725 = t741 + (-t774 * t829 - t748) * pkin(9) + (t775 * t829 - t747) * pkin(5);
t841 = sin(qJ(6));
t845 = cos(qJ(6));
t718 = -t721 * t841 + t725 * t845;
t766 = -t775 * t841 + t829 * t845;
t732 = qJD(6) * t766 + t748 * t845 + t828 * t841;
t746 = qJDD(6) - t747;
t767 = t775 * t845 + t829 * t841;
t749 = -mrSges(7,1) * t766 + mrSges(7,2) * t767;
t770 = qJD(6) - t774;
t750 = -mrSges(7,2) * t770 + mrSges(7,3) * t766;
t715 = m(7) * t718 + mrSges(7,1) * t746 - t732 * mrSges(7,3) - t749 * t767 + t750 * t770;
t719 = t721 * t845 + t725 * t841;
t731 = -qJD(6) * t767 - t748 * t841 + t828 * t845;
t751 = mrSges(7,1) * t770 - mrSges(7,3) * t767;
t716 = m(7) * t719 - mrSges(7,2) * t746 + t731 * mrSges(7,3) + t749 * t766 - t751 * t770;
t862 = -t715 * t841 + t845 * t716;
t702 = m(6) * t724 - mrSges(6,2) * t828 + mrSges(6,3) * t747 + t756 * t774 - t769 * t829 + t862;
t723 = t727 * t846 - t729 * t842;
t768 = -mrSges(6,2) * t829 + mrSges(6,3) * t774;
t720 = -pkin(5) * t828 - pkin(9) * t827 + t757 * t775 - t723;
t857 = -m(7) * t720 + t731 * mrSges(7,1) - t732 * mrSges(7,2) + t766 * t750 - t751 * t767;
t711 = m(6) * t723 + mrSges(6,1) * t828 - mrSges(6,3) * t748 - t756 * t775 + t768 * t829 + t857;
t695 = t842 * t702 + t846 * t711;
t778 = -mrSges(5,1) * t803 + mrSges(5,2) * t804;
t791 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t803;
t692 = m(5) * t738 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t783 + qJD(3) * t791 - t778 * t804 + t695;
t792 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t804;
t863 = t846 * t702 - t711 * t842;
t693 = m(5) * t739 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t782 - qJD(3) * t792 + t778 * t803 + t863;
t686 = t840 * t692 + t839 * t693;
t814 = (mrSges(4,1) * t843 + mrSges(4,2) * t847) * qJD(1);
t871 = qJD(1) * t843;
t818 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t871;
t683 = m(4) * t785 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t816 + qJD(3) * t818 - t814 * t870 + t686;
t820 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t870;
t864 = -t692 * t839 + t840 * t693;
t684 = m(4) * t786 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t815 - qJD(3) * t820 - t814 * t871 + t864;
t680 = t847 * t683 + t843 * t684;
t801 = -qJDD(1) * pkin(1) + t859;
t858 = -m(3) * t801 + (t849 * mrSges(3,3)) - t680;
t676 = m(2) * t821 - (t849 * mrSges(2,2)) + qJDD(1) * t875 + t858;
t797 = t849 * pkin(1) - t860;
t705 = t845 * t715 + t841 * t716;
t856 = m(6) * t741 - t747 * mrSges(6,1) + t748 * mrSges(6,2) - t774 * t768 + t775 * t769 + t705;
t703 = m(5) * t764 - t782 * mrSges(5,1) + t783 * mrSges(5,2) - t803 * t791 + t804 * t792 + t856;
t790 = t849 * t876 + t860;
t852 = -m(4) * t790 + t815 * mrSges(4,1) - t816 * mrSges(4,2) - t818 * t871 - t820 * t870 - t703;
t851 = -m(3) * t797 + (t849 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t852;
t698 = m(2) * t822 - (t849 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t851;
t872 = t848 * t676 + t844 * t698;
t866 = -t676 * t844 + t848 * t698;
t865 = -t843 * t683 + t847 * t684;
t734 = Ifges(7,5) * t767 + Ifges(7,6) * t766 + Ifges(7,3) * t770;
t736 = Ifges(7,1) * t767 + Ifges(7,4) * t766 + Ifges(7,5) * t770;
t708 = -mrSges(7,1) * t720 + mrSges(7,3) * t719 + Ifges(7,4) * t732 + Ifges(7,2) * t731 + Ifges(7,6) * t746 - t734 * t767 + t736 * t770;
t735 = Ifges(7,4) * t767 + Ifges(7,2) * t766 + Ifges(7,6) * t770;
t709 = mrSges(7,2) * t720 - mrSges(7,3) * t718 + Ifges(7,1) * t732 + Ifges(7,4) * t731 + Ifges(7,5) * t746 + t734 * t766 - t735 * t770;
t753 = Ifges(6,4) * t775 + Ifges(6,2) * t774 + Ifges(6,6) * t829;
t754 = Ifges(6,1) * t775 + Ifges(6,4) * t774 + Ifges(6,5) * t829;
t855 = mrSges(6,1) * t723 - mrSges(6,2) * t724 + Ifges(6,5) * t748 + Ifges(6,6) * t747 + Ifges(6,3) * t828 + pkin(5) * t857 + pkin(9) * t862 + t845 * t708 + t841 * t709 + t775 * t753 - t774 * t754;
t752 = Ifges(6,5) * t775 + Ifges(6,6) * t774 + Ifges(6,3) * t829;
t687 = mrSges(6,2) * t741 - mrSges(6,3) * t723 + Ifges(6,1) * t748 + Ifges(6,4) * t747 + Ifges(6,5) * t828 - pkin(9) * t705 - t708 * t841 + t709 * t845 + t752 * t774 - t753 * t829;
t853 = mrSges(7,1) * t718 - mrSges(7,2) * t719 + Ifges(7,5) * t732 + Ifges(7,6) * t731 + Ifges(7,3) * t746 + t735 * t767 - t736 * t766;
t688 = -mrSges(6,1) * t741 + mrSges(6,3) * t724 + Ifges(6,4) * t748 + Ifges(6,2) * t747 + Ifges(6,6) * t828 - pkin(5) * t705 - t752 * t775 + t754 * t829 - t853;
t771 = Ifges(5,5) * t804 + Ifges(5,6) * t803 + (Ifges(5,3) * qJD(3));
t773 = Ifges(5,1) * t804 + Ifges(5,4) * t803 + Ifges(5,5) * qJD(3);
t674 = -mrSges(5,1) * t764 + mrSges(5,3) * t739 + Ifges(5,4) * t783 + Ifges(5,2) * t782 + Ifges(5,6) * qJDD(3) - pkin(4) * t856 + pkin(8) * t863 + qJD(3) * t773 + t842 * t687 + t846 * t688 - t804 * t771;
t772 = Ifges(5,4) * t804 + Ifges(5,2) * t803 + Ifges(5,6) * qJD(3);
t681 = mrSges(5,2) * t764 - mrSges(5,3) * t738 + Ifges(5,1) * t783 + Ifges(5,4) * t782 + Ifges(5,5) * qJDD(3) - pkin(8) * t695 - qJD(3) * t772 + t687 * t846 - t688 * t842 + t771 * t803;
t805 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t847 - Ifges(4,6) * t843) * qJD(1);
t807 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t847 - Ifges(4,4) * t843) * qJD(1);
t671 = -mrSges(4,1) * t790 + mrSges(4,3) * t786 + Ifges(4,4) * t816 + Ifges(4,2) * t815 + Ifges(4,6) * qJDD(3) - pkin(3) * t703 + qJ(4) * t864 + qJD(3) * t807 + t840 * t674 + t839 * t681 - t805 * t870;
t806 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t847 - Ifges(4,2) * t843) * qJD(1);
t673 = mrSges(4,2) * t790 - mrSges(4,3) * t785 + Ifges(4,1) * t816 + Ifges(4,4) * t815 + Ifges(4,5) * qJDD(3) - qJ(4) * t686 - qJD(3) * t806 - t674 * t839 + t681 * t840 - t805 * t871;
t678 = qJDD(1) * mrSges(3,2) - t858;
t854 = mrSges(2,1) * t821 - mrSges(2,2) * t822 + mrSges(3,2) * t801 - mrSges(3,3) * t797 - pkin(1) * t678 - pkin(7) * t680 + qJ(2) * t851 - t671 * t843 + t847 * t673 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t850 = mrSges(4,1) * t785 + mrSges(5,1) * t738 - mrSges(4,2) * t786 - mrSges(5,2) * t739 + Ifges(4,5) * t816 + Ifges(5,5) * t783 + Ifges(4,6) * t815 + Ifges(5,6) * t782 + pkin(3) * t686 + pkin(4) * t695 + t804 * t772 - t803 * t773 + t806 * t870 + t807 * t871 + t855 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3);
t679 = -m(3) * g(3) + t865;
t670 = -qJ(2) * t679 + t850 + pkin(2) * t680 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t873 * t849) - mrSges(2,3) * t821 + mrSges(3,1) * t801 + t874 * qJDD(1);
t669 = -mrSges(3,1) * t797 + mrSges(2,3) * t822 - pkin(1) * t679 - pkin(2) * t852 - pkin(7) * t865 + g(3) * t875 - qJDD(1) * t873 - t847 * t671 - t843 * t673 + t849 * t874;
t1 = [-m(1) * g(1) + t866; -m(1) * g(2) + t872; (-m(1) - m(2) - m(3)) * g(3) + t865; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t872 - t844 * t669 + t848 * t670; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t866 + t848 * t669 + t844 * t670; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t854; t854; t678; t850; t703; t855; t853;];
tauJB  = t1;
