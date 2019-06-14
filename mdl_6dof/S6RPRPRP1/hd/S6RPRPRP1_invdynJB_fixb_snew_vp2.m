% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:28:58
% EndTime: 2019-05-05 17:29:05
% DurationCPUTime: 6.74s
% Computational Cost: add. (100838->321), mult. (216703->391), div. (0->0), fcn. (140775->10), ass. (0->130)
t886 = -2 * qJD(4);
t885 = Ifges(6,1) + Ifges(7,1);
t879 = Ifges(6,4) + Ifges(7,4);
t878 = Ifges(6,5) + Ifges(7,5);
t884 = Ifges(6,2) + Ifges(7,2);
t877 = Ifges(6,6) + Ifges(7,6);
t883 = Ifges(6,3) + Ifges(7,3);
t844 = sin(qJ(1));
t847 = cos(qJ(1));
t828 = t844 * g(1) - t847 * g(2);
t819 = qJDD(1) * pkin(1) + t828;
t829 = -t847 * g(1) - t844 * g(2);
t849 = qJD(1) ^ 2;
t821 = -t849 * pkin(1) + t829;
t839 = sin(pkin(9));
t841 = cos(pkin(9));
t797 = t839 * t819 + t841 * t821;
t788 = -t849 * pkin(2) + qJDD(1) * pkin(7) + t797;
t837 = -g(3) + qJDD(2);
t843 = sin(qJ(3));
t846 = cos(qJ(3));
t776 = -t843 * t788 + t846 * t837;
t866 = qJD(1) * qJD(3);
t863 = t846 * t866;
t822 = t843 * qJDD(1) + t863;
t755 = (-t822 + t863) * qJ(4) + (t843 * t846 * t849 + qJDD(3)) * pkin(3) + t776;
t777 = t846 * t788 + t843 * t837;
t823 = t846 * qJDD(1) - t843 * t866;
t869 = qJD(1) * t843;
t825 = qJD(3) * pkin(3) - qJ(4) * t869;
t836 = t846 ^ 2;
t756 = -t836 * t849 * pkin(3) + t823 * qJ(4) - qJD(3) * t825 + t777;
t838 = sin(pkin(10));
t840 = cos(pkin(10));
t809 = (t838 * t846 + t840 * t843) * qJD(1);
t747 = t840 * t755 - t838 * t756 + t809 * t886;
t808 = (t838 * t843 - t840 * t846) * qJD(1);
t748 = t838 * t755 + t840 * t756 + t808 * t886;
t791 = t808 * pkin(4) - t809 * pkin(8);
t848 = qJD(3) ^ 2;
t746 = -t848 * pkin(4) + qJDD(3) * pkin(8) - t808 * t791 + t748;
t796 = t841 * t819 - t839 * t821;
t854 = -qJDD(1) * pkin(2) - t796;
t758 = -t823 * pkin(3) + qJDD(4) + t825 * t869 + (-qJ(4) * t836 - pkin(7)) * t849 + t854;
t798 = -t838 * t822 + t840 * t823;
t799 = t840 * t822 + t838 * t823;
t751 = (qJD(3) * t808 - t799) * pkin(8) + (qJD(3) * t809 - t798) * pkin(4) + t758;
t842 = sin(qJ(5));
t845 = cos(qJ(5));
t741 = -t842 * t746 + t845 * t751;
t801 = t845 * qJD(3) - t842 * t809;
t771 = t801 * qJD(5) + t842 * qJDD(3) + t845 * t799;
t802 = t842 * qJD(3) + t845 * t809;
t773 = -t801 * mrSges(7,1) + t802 * mrSges(7,2);
t774 = -t801 * mrSges(6,1) + t802 * mrSges(6,2);
t807 = qJD(5) + t808;
t779 = -t807 * mrSges(6,2) + t801 * mrSges(6,3);
t795 = qJDD(5) - t798;
t738 = -0.2e1 * qJD(6) * t802 + (t801 * t807 - t771) * qJ(6) + (t801 * t802 + t795) * pkin(5) + t741;
t778 = -t807 * mrSges(7,2) + t801 * mrSges(7,3);
t865 = m(7) * t738 + t795 * mrSges(7,1) + t807 * t778;
t728 = m(6) * t741 + t795 * mrSges(6,1) + t807 * t779 + (-t773 - t774) * t802 + (-mrSges(6,3) - mrSges(7,3)) * t771 + t865;
t742 = t845 * t746 + t842 * t751;
t770 = -t802 * qJD(5) + t845 * qJDD(3) - t842 * t799;
t780 = t807 * pkin(5) - t802 * qJ(6);
t800 = t801 ^ 2;
t740 = -t800 * pkin(5) + t770 * qJ(6) + 0.2e1 * qJD(6) * t801 - t807 * t780 + t742;
t864 = m(7) * t740 + t770 * mrSges(7,3) + t801 * t773;
t781 = t807 * mrSges(7,1) - t802 * mrSges(7,3);
t870 = -t807 * mrSges(6,1) + t802 * mrSges(6,3) - t781;
t880 = -mrSges(6,2) - mrSges(7,2);
t731 = m(6) * t742 + t770 * mrSges(6,3) + t801 * t774 + t880 * t795 + t870 * t807 + t864;
t726 = -t842 * t728 + t845 * t731;
t790 = t808 * mrSges(5,1) + t809 * mrSges(5,2);
t804 = qJD(3) * mrSges(5,1) - t809 * mrSges(5,3);
t722 = m(5) * t748 - qJDD(3) * mrSges(5,2) + t798 * mrSges(5,3) - qJD(3) * t804 - t808 * t790 + t726;
t745 = -qJDD(3) * pkin(4) - t848 * pkin(8) + t809 * t791 - t747;
t743 = -t770 * pkin(5) - t800 * qJ(6) + t802 * t780 + qJDD(6) + t745;
t857 = -m(7) * t743 + t770 * mrSges(7,1) + t801 * t778;
t734 = -m(6) * t745 + t770 * mrSges(6,1) + t880 * t771 + t801 * t779 + t870 * t802 + t857;
t803 = -qJD(3) * mrSges(5,2) - t808 * mrSges(5,3);
t733 = m(5) * t747 + qJDD(3) * mrSges(5,1) - t799 * mrSges(5,3) + qJD(3) * t803 - t809 * t790 + t734;
t715 = t838 * t722 + t840 * t733;
t736 = t771 * mrSges(7,2) + t802 * t781 - t857;
t871 = t879 * t801 + t885 * t802 + t878 * t807;
t873 = -t877 * t801 - t878 * t802 - t883 * t807;
t716 = -mrSges(6,1) * t745 + mrSges(6,3) * t742 - mrSges(7,1) * t743 + mrSges(7,3) * t740 - pkin(5) * t736 + qJ(6) * t864 + (-qJ(6) * t781 + t871) * t807 + t873 * t802 + (-qJ(6) * mrSges(7,2) + t877) * t795 + t879 * t771 + t884 * t770;
t735 = -t771 * mrSges(7,3) - t802 * t773 + t865;
t872 = -t884 * t801 - t879 * t802 - t877 * t807;
t723 = mrSges(6,2) * t745 + mrSges(7,2) * t743 - mrSges(6,3) * t741 - mrSges(7,3) * t738 - qJ(6) * t735 + t879 * t770 + t885 * t771 + t878 * t795 - t873 * t801 + t872 * t807;
t785 = Ifges(5,4) * t809 - Ifges(5,2) * t808 + Ifges(5,6) * qJD(3);
t786 = Ifges(5,1) * t809 - Ifges(5,4) * t808 + Ifges(5,5) * qJD(3);
t814 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t843 + Ifges(4,2) * t846) * qJD(1);
t815 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t843 + Ifges(4,4) * t846) * qJD(1);
t882 = (t843 * t814 - t846 * t815) * qJD(1) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + mrSges(4,1) * t776 + mrSges(5,1) * t747 - mrSges(4,2) * t777 - mrSges(5,2) * t748 + Ifges(4,5) * t822 + Ifges(5,5) * t799 + Ifges(4,6) * t823 + Ifges(5,6) * t798 + pkin(3) * t715 + pkin(4) * t734 + pkin(8) * t726 + t845 * t716 + t842 * t723 + t809 * t785 + t808 * t786;
t881 = mrSges(6,1) * t741 + mrSges(7,1) * t738 - mrSges(6,2) * t742 - mrSges(7,2) * t740 + pkin(5) * t735 + t877 * t770 + t878 * t771 + t883 * t795 - t871 * t801 - t872 * t802;
t820 = (-mrSges(4,1) * t846 + mrSges(4,2) * t843) * qJD(1);
t868 = qJD(1) * t846;
t827 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t868;
t713 = m(4) * t776 + qJDD(3) * mrSges(4,1) - t822 * mrSges(4,3) + qJD(3) * t827 - t820 * t869 + t715;
t826 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t869;
t859 = t840 * t722 - t838 * t733;
t714 = m(4) * t777 - qJDD(3) * mrSges(4,2) + t823 * mrSges(4,3) - qJD(3) * t826 + t820 * t868 + t859;
t860 = -t843 * t713 + t846 * t714;
t705 = m(3) * t797 - t849 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t860;
t725 = t845 * t728 + t842 * t731;
t724 = m(5) * t758 - t798 * mrSges(5,1) + t799 * mrSges(5,2) + t808 * t803 + t809 * t804 + t725;
t787 = -t849 * pkin(7) + t854;
t851 = -m(4) * t787 + t823 * mrSges(4,1) - t822 * mrSges(4,2) - t826 * t869 + t827 * t868 - t724;
t718 = m(3) * t796 + qJDD(1) * mrSges(3,1) - t849 * mrSges(3,2) + t851;
t701 = t839 * t705 + t841 * t718;
t698 = m(2) * t828 + qJDD(1) * mrSges(2,1) - t849 * mrSges(2,2) + t701;
t861 = t841 * t705 - t839 * t718;
t699 = m(2) * t829 - t849 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t861;
t874 = t847 * t698 + t844 * t699;
t708 = t846 * t713 + t843 * t714;
t706 = m(3) * t837 + t708;
t862 = -t844 * t698 + t847 * t699;
t784 = Ifges(5,5) * t809 - Ifges(5,6) * t808 + Ifges(5,3) * qJD(3);
t702 = mrSges(5,2) * t758 - mrSges(5,3) * t747 + Ifges(5,1) * t799 + Ifges(5,4) * t798 + Ifges(5,5) * qJDD(3) - pkin(8) * t725 - qJD(3) * t785 - t842 * t716 + t845 * t723 - t808 * t784;
t709 = -mrSges(5,1) * t758 + mrSges(5,3) * t748 + Ifges(5,4) * t799 + Ifges(5,2) * t798 + Ifges(5,6) * qJDD(3) - pkin(4) * t725 + qJD(3) * t786 - t809 * t784 - t881;
t813 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t843 + Ifges(4,6) * t846) * qJD(1);
t692 = -mrSges(4,1) * t787 + mrSges(4,3) * t777 + Ifges(4,4) * t822 + Ifges(4,2) * t823 + Ifges(4,6) * qJDD(3) - pkin(3) * t724 + qJ(4) * t859 + qJD(3) * t815 + t838 * t702 + t840 * t709 - t813 * t869;
t694 = mrSges(4,2) * t787 - mrSges(4,3) * t776 + Ifges(4,1) * t822 + Ifges(4,4) * t823 + Ifges(4,5) * qJDD(3) - qJ(4) * t715 - qJD(3) * t814 + t840 * t702 - t838 * t709 + t813 * t868;
t852 = mrSges(2,1) * t828 + mrSges(3,1) * t796 - mrSges(2,2) * t829 - mrSges(3,2) * t797 + pkin(1) * t701 + pkin(2) * t851 + pkin(7) * t860 + t846 * t692 + t843 * t694 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t690 = -mrSges(3,1) * t837 + mrSges(3,3) * t797 + t849 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t708 - t882;
t689 = mrSges(3,2) * t837 - mrSges(3,3) * t796 + Ifges(3,5) * qJDD(1) - t849 * Ifges(3,6) - pkin(7) * t708 - t843 * t692 + t846 * t694;
t688 = -mrSges(2,2) * g(3) - mrSges(2,3) * t828 + Ifges(2,5) * qJDD(1) - t849 * Ifges(2,6) - qJ(2) * t701 + t841 * t689 - t839 * t690;
t687 = mrSges(2,1) * g(3) + mrSges(2,3) * t829 + t849 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t706 + qJ(2) * t861 + t839 * t689 + t841 * t690;
t1 = [-m(1) * g(1) + t862; -m(1) * g(2) + t874; (-m(1) - m(2)) * g(3) + t706; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t874 - t844 * t687 + t847 * t688; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t862 + t847 * t687 + t844 * t688; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t852; t852; t706; t882; t724; t881; t736;];
tauJB  = t1;
