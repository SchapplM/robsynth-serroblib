% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-05 14:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:31:47
% EndTime: 2019-05-05 14:31:55
% DurationCPUTime: 7.66s
% Computational Cost: add. (112524->313), mult. (264500->382), div. (0->0), fcn. (186242->10), ass. (0->135)
t841 = sin(qJ(1));
t844 = cos(qJ(1));
t816 = t841 * g(1) - t844 * g(2);
t846 = qJD(1) ^ 2;
t855 = -t846 * qJ(2) + qJDD(2) - t816;
t878 = -pkin(1) - qJ(3);
t886 = -(2 * qJD(1) * qJD(3)) + t878 * qJDD(1) + t855;
t817 = -t844 * g(1) - t841 * g(2);
t885 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t817;
t807 = t846 * pkin(1) - t885;
t884 = -m(3) * t807 + t846 * mrSges(3,2) + qJDD(1) * mrSges(3,3);
t882 = pkin(3) * t846;
t881 = mrSges(2,1) - mrSges(3,2);
t880 = -Ifges(3,4) + Ifges(2,5);
t879 = -Ifges(2,6) + Ifges(3,5);
t838 = cos(pkin(9));
t877 = mrSges(4,2) * t838;
t836 = sin(pkin(9));
t791 = t836 * g(3) + t886 * t838;
t771 = (-pkin(7) * qJDD(1) - t836 * t882) * t838 + t791;
t792 = -t838 * g(3) + t886 * t836;
t826 = t836 ^ 2;
t870 = qJDD(1) * t836;
t772 = -pkin(7) * t870 - t826 * t882 + t792;
t840 = sin(qJ(4));
t843 = cos(qJ(4));
t756 = t840 * t771 + t843 * t772;
t859 = t836 * t843 + t838 * t840;
t811 = t859 * qJD(1);
t858 = -t836 * t840 + t838 * t843;
t812 = t858 * qJD(1);
t786 = t811 * mrSges(5,1) + t812 * mrSges(5,2);
t872 = t812 * qJD(4);
t793 = t859 * qJDD(1) + t872;
t805 = qJD(4) * mrSges(5,1) - t812 * mrSges(5,3);
t785 = t811 * pkin(4) - t812 * qJ(5);
t845 = qJD(4) ^ 2;
t742 = -t845 * pkin(4) + qJDD(4) * qJ(5) - t811 * t785 + t756;
t854 = qJDD(3) + t885;
t875 = -t838 ^ 2 - t826;
t780 = pkin(3) * t870 + (t875 * pkin(7) + t878) * t846 + t854;
t873 = t811 * qJD(4);
t794 = t858 * qJDD(1) - t873;
t748 = (-t794 + t873) * qJ(5) + (t793 + t872) * pkin(4) + t780;
t835 = sin(pkin(10));
t837 = cos(pkin(10));
t799 = t835 * qJD(4) + t837 * t812;
t737 = -0.2e1 * qJD(5) * t799 - t835 * t742 + t837 * t748;
t778 = t835 * qJDD(4) + t837 * t794;
t798 = t837 * qJD(4) - t835 * t812;
t735 = (t798 * t811 - t778) * pkin(8) + (t798 * t799 + t793) * pkin(5) + t737;
t738 = 0.2e1 * qJD(5) * t798 + t837 * t742 + t835 * t748;
t776 = t811 * pkin(5) - t799 * pkin(8);
t777 = t837 * qJDD(4) - t835 * t794;
t797 = t798 ^ 2;
t736 = -t797 * pkin(5) + t777 * pkin(8) - t811 * t776 + t738;
t839 = sin(qJ(6));
t842 = cos(qJ(6));
t733 = t842 * t735 - t839 * t736;
t763 = t842 * t798 - t839 * t799;
t747 = t763 * qJD(6) + t839 * t777 + t842 * t778;
t764 = t839 * t798 + t842 * t799;
t753 = -t763 * mrSges(7,1) + t764 * mrSges(7,2);
t809 = qJD(6) + t811;
t757 = -t809 * mrSges(7,2) + t763 * mrSges(7,3);
t790 = qJDD(6) + t793;
t729 = m(7) * t733 + t790 * mrSges(7,1) - t747 * mrSges(7,3) - t764 * t753 + t809 * t757;
t734 = t839 * t735 + t842 * t736;
t746 = -t764 * qJD(6) + t842 * t777 - t839 * t778;
t758 = t809 * mrSges(7,1) - t764 * mrSges(7,3);
t730 = m(7) * t734 - t790 * mrSges(7,2) + t746 * mrSges(7,3) + t763 * t753 - t809 * t758;
t721 = t842 * t729 + t839 * t730;
t766 = -t798 * mrSges(6,1) + t799 * mrSges(6,2);
t774 = -t811 * mrSges(6,2) + t798 * mrSges(6,3);
t719 = m(6) * t737 + t793 * mrSges(6,1) - t778 * mrSges(6,3) - t799 * t766 + t811 * t774 + t721;
t775 = t811 * mrSges(6,1) - t799 * mrSges(6,3);
t863 = -t839 * t729 + t842 * t730;
t720 = m(6) * t738 - t793 * mrSges(6,2) + t777 * mrSges(6,3) + t798 * t766 - t811 * t775 + t863;
t864 = -t835 * t719 + t837 * t720;
t713 = m(5) * t756 - qJDD(4) * mrSges(5,2) - t793 * mrSges(5,3) - qJD(4) * t805 - t811 * t786 + t864;
t755 = t843 * t771 - t840 * t772;
t741 = -qJDD(4) * pkin(4) - t845 * qJ(5) + t812 * t785 + qJDD(5) - t755;
t739 = -t777 * pkin(5) - t797 * pkin(8) + t799 * t776 + t741;
t852 = m(7) * t739 - t746 * mrSges(7,1) + t747 * mrSges(7,2) - t763 * t757 + t764 * t758;
t732 = m(6) * t741 - t777 * mrSges(6,1) + t778 * mrSges(6,2) - t798 * t774 + t799 * t775 + t852;
t804 = -qJD(4) * mrSges(5,2) - t811 * mrSges(5,3);
t725 = m(5) * t755 + qJDD(4) * mrSges(5,1) - t794 * mrSges(5,3) + qJD(4) * t804 - t812 * t786 - t732;
t701 = t840 * t713 + t843 * t725;
t857 = -qJDD(1) * mrSges(4,3) - t846 * (mrSges(4,1) * t836 + t877);
t699 = m(4) * t791 + t857 * t838 + t701;
t865 = t843 * t713 - t840 * t725;
t700 = m(4) * t792 + t857 * t836 + t865;
t696 = t838 * t699 + t836 * t700;
t810 = -qJDD(1) * pkin(1) + t855;
t853 = -m(3) * t810 + t846 * mrSges(3,3) - t696;
t692 = m(2) * t816 - t846 * mrSges(2,2) + t881 * qJDD(1) + t853;
t803 = t878 * t846 + t854;
t715 = t837 * t719 + t835 * t720;
t851 = m(5) * t780 + t793 * mrSges(5,1) + t794 * mrSges(5,2) + t811 * t804 + t812 * t805 + t715;
t850 = m(4) * t803 + mrSges(4,1) * t870 + qJDD(1) * t877 + t851;
t868 = t875 * mrSges(4,3);
t708 = t850 + m(2) * t817 + (-mrSges(2,1) + t868) * t846 - qJDD(1) * mrSges(2,2) + t884;
t876 = t844 * t692 + t841 * t708;
t860 = Ifges(4,5) * t838 - Ifges(4,6) * t836;
t874 = t846 * t860;
t867 = -t841 * t692 + t844 * t708;
t866 = -t836 * t699 + t838 * t700;
t862 = Ifges(4,1) * t838 - Ifges(4,4) * t836;
t861 = Ifges(4,4) * t838 - Ifges(4,2) * t836;
t749 = Ifges(7,5) * t764 + Ifges(7,6) * t763 + Ifges(7,3) * t809;
t751 = Ifges(7,1) * t764 + Ifges(7,4) * t763 + Ifges(7,5) * t809;
t722 = -mrSges(7,1) * t739 + mrSges(7,3) * t734 + Ifges(7,4) * t747 + Ifges(7,2) * t746 + Ifges(7,6) * t790 - t764 * t749 + t809 * t751;
t750 = Ifges(7,4) * t764 + Ifges(7,2) * t763 + Ifges(7,6) * t809;
t723 = mrSges(7,2) * t739 - mrSges(7,3) * t733 + Ifges(7,1) * t747 + Ifges(7,4) * t746 + Ifges(7,5) * t790 + t763 * t749 - t809 * t750;
t759 = Ifges(6,5) * t799 + Ifges(6,6) * t798 + Ifges(6,3) * t811;
t761 = Ifges(6,1) * t799 + Ifges(6,4) * t798 + Ifges(6,5) * t811;
t703 = -mrSges(6,1) * t741 + mrSges(6,3) * t738 + Ifges(6,4) * t778 + Ifges(6,2) * t777 + Ifges(6,6) * t793 - pkin(5) * t852 + pkin(8) * t863 + t842 * t722 + t839 * t723 - t799 * t759 + t811 * t761;
t760 = Ifges(6,4) * t799 + Ifges(6,2) * t798 + Ifges(6,6) * t811;
t705 = mrSges(6,2) * t741 - mrSges(6,3) * t737 + Ifges(6,1) * t778 + Ifges(6,4) * t777 + Ifges(6,5) * t793 - pkin(8) * t721 - t839 * t722 + t842 * t723 + t798 * t759 - t811 * t760;
t782 = Ifges(5,4) * t812 - Ifges(5,2) * t811 + Ifges(5,6) * qJD(4);
t783 = Ifges(5,1) * t812 - Ifges(5,4) * t811 + Ifges(5,5) * qJD(4);
t849 = mrSges(5,1) * t755 - mrSges(5,2) * t756 + Ifges(5,5) * t794 - Ifges(5,6) * t793 + Ifges(5,3) * qJDD(4) - pkin(4) * t732 + qJ(5) * t864 + t837 * t703 + t835 * t705 + t812 * t782 + t811 * t783;
t781 = Ifges(5,5) * t812 - Ifges(5,6) * t811 + Ifges(5,3) * qJD(4);
t690 = mrSges(5,2) * t780 - mrSges(5,3) * t755 + Ifges(5,1) * t794 - Ifges(5,4) * t793 + Ifges(5,5) * qJDD(4) - qJ(5) * t715 - qJD(4) * t782 - t835 * t703 + t837 * t705 - t811 * t781;
t847 = mrSges(7,1) * t733 - mrSges(7,2) * t734 + Ifges(7,5) * t747 + Ifges(7,6) * t746 + Ifges(7,3) * t790 + t764 * t750 - t763 * t751;
t697 = -t847 + (-Ifges(5,2) - Ifges(6,3)) * t793 + Ifges(5,6) * qJDD(4) - t812 * t781 + t798 * t761 - t799 * t760 + Ifges(5,4) * t794 - Ifges(6,5) * t778 - mrSges(5,1) * t780 + qJD(4) * t783 - Ifges(6,6) * t777 + mrSges(5,3) * t756 - mrSges(6,1) * t737 + mrSges(6,2) * t738 - pkin(5) * t721 - pkin(4) * t715;
t687 = -mrSges(4,1) * t803 + mrSges(4,3) * t792 - pkin(3) * t851 + pkin(7) * t865 + t861 * qJDD(1) + t840 * t690 + t843 * t697 - t838 * t874;
t689 = mrSges(4,2) * t803 - mrSges(4,3) * t791 - pkin(7) * t701 + t862 * qJDD(1) + t843 * t690 - t840 * t697 - t836 * t874;
t694 = qJDD(1) * mrSges(3,2) - t853;
t710 = t846 * t868 + t850;
t848 = -mrSges(2,2) * t817 - mrSges(3,3) * t807 - pkin(1) * t694 - qJ(3) * t696 - t836 * t687 + t838 * t689 + qJ(2) * (t710 + t884) + mrSges(3,2) * t810 + mrSges(2,1) * t816 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t695 = -m(3) * g(3) + t866;
t686 = (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t860 + t880) * qJDD(1) + t849 - mrSges(2,3) * t816 + mrSges(3,1) * t810 + mrSges(4,1) * t791 - mrSges(4,2) * t792 + pkin(3) * t701 - qJ(2) * t695 + pkin(2) * t696 + (t836 * t862 + t838 * t861 + t879) * t846;
t685 = -mrSges(3,1) * t807 + mrSges(2,3) * t817 - pkin(1) * t695 + pkin(2) * t710 + t881 * g(3) - qJ(3) * t866 - t879 * qJDD(1) - t838 * t687 - t836 * t689 + t880 * t846;
t1 = [-m(1) * g(1) + t867; -m(1) * g(2) + t876; (-m(1) - m(2) - m(3)) * g(3) + t866; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t876 - t841 * t685 + t844 * t686; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t867 + t844 * t685 + t841 * t686; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t848; t848; t694; t710; t849; t732; t847;];
tauJB  = t1;
