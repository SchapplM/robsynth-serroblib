% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR10
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
% Datum: 2019-05-05 20:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:02:14
% EndTime: 2019-05-05 20:02:24
% DurationCPUTime: 9.94s
% Computational Cost: add. (166053->338), mult. (350921->412), div. (0->0), fcn. (234097->10), ass. (0->136)
t835 = sin(qJ(1));
t839 = cos(qJ(1));
t817 = -t839 * g(1) - t835 * g(2);
t868 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t817;
t867 = -pkin(1) - pkin(7);
t866 = mrSges(2,1) - mrSges(3,2);
t865 = -Ifges(3,4) + Ifges(2,5);
t864 = (Ifges(3,5) - Ifges(2,6));
t816 = t835 * g(1) - t839 * g(2);
t841 = qJD(1) ^ 2;
t850 = -t841 * qJ(2) + qJDD(2) - t816;
t786 = qJDD(1) * t867 + t850;
t834 = sin(qJ(3));
t838 = cos(qJ(3));
t776 = -t838 * g(3) + t834 * t786;
t810 = (mrSges(4,1) * t834 + mrSges(4,2) * t838) * qJD(1);
t861 = qJD(1) * qJD(3);
t858 = t838 * t861;
t811 = t834 * qJDD(1) + t858;
t862 = qJD(1) * t838;
t815 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t862;
t820 = t834 * qJD(1);
t780 = t841 * t867 - t868;
t859 = t834 * t861;
t812 = t838 * qJDD(1) - t859;
t758 = (-t812 + t859) * qJ(4) + (t811 + t858) * pkin(3) + t780;
t809 = (pkin(3) * t834 - qJ(4) * t838) * qJD(1);
t840 = qJD(3) ^ 2;
t761 = -t840 * pkin(3) + qJDD(3) * qJ(4) - t809 * t820 + t776;
t830 = sin(pkin(10));
t831 = cos(pkin(10));
t804 = t830 * qJD(3) + t831 * t862;
t741 = -0.2e1 * qJD(4) * t804 + t831 * t758 - t830 * t761;
t784 = t830 * qJDD(3) + t831 * t812;
t803 = t831 * qJD(3) - t830 * t862;
t732 = (t803 * t820 - t784) * pkin(8) + (t803 * t804 + t811) * pkin(4) + t741;
t742 = 0.2e1 * qJD(4) * t803 + t830 * t758 + t831 * t761;
t783 = t831 * qJDD(3) - t830 * t812;
t785 = pkin(4) * t820 - t804 * pkin(8);
t802 = t803 ^ 2;
t734 = -t802 * pkin(4) + t783 * pkin(8) - t785 * t820 + t742;
t833 = sin(qJ(5));
t837 = cos(qJ(5));
t719 = t837 * t732 - t833 * t734;
t771 = t837 * t803 - t833 * t804;
t748 = t771 * qJD(5) + t833 * t783 + t837 * t784;
t772 = t833 * t803 + t837 * t804;
t808 = qJDD(5) + t811;
t819 = t820 + qJD(5);
t717 = (t771 * t819 - t748) * pkin(9) + (t771 * t772 + t808) * pkin(5) + t719;
t720 = t833 * t732 + t837 * t734;
t747 = -t772 * qJD(5) + t837 * t783 - t833 * t784;
t764 = t819 * pkin(5) - t772 * pkin(9);
t770 = t771 ^ 2;
t718 = -t770 * pkin(5) + t747 * pkin(9) - t819 * t764 + t720;
t832 = sin(qJ(6));
t836 = cos(qJ(6));
t715 = t836 * t717 - t832 * t718;
t753 = t836 * t771 - t832 * t772;
t728 = t753 * qJD(6) + t832 * t747 + t836 * t748;
t754 = t832 * t771 + t836 * t772;
t740 = -t753 * mrSges(7,1) + t754 * mrSges(7,2);
t818 = qJD(6) + t819;
t744 = -t818 * mrSges(7,2) + t753 * mrSges(7,3);
t799 = qJDD(6) + t808;
t710 = m(7) * t715 + t799 * mrSges(7,1) - t728 * mrSges(7,3) - t754 * t740 + t818 * t744;
t716 = t832 * t717 + t836 * t718;
t727 = -t754 * qJD(6) + t836 * t747 - t832 * t748;
t745 = t818 * mrSges(7,1) - t754 * mrSges(7,3);
t711 = m(7) * t716 - t799 * mrSges(7,2) + t727 * mrSges(7,3) + t753 * t740 - t818 * t745;
t703 = t836 * t710 + t832 * t711;
t755 = -t771 * mrSges(6,1) + t772 * mrSges(6,2);
t762 = -t819 * mrSges(6,2) + t771 * mrSges(6,3);
t701 = m(6) * t719 + t808 * mrSges(6,1) - t748 * mrSges(6,3) - t772 * t755 + t819 * t762 + t703;
t763 = t819 * mrSges(6,1) - t772 * mrSges(6,3);
t853 = -t832 * t710 + t836 * t711;
t702 = m(6) * t720 - t808 * mrSges(6,2) + t747 * mrSges(6,3) + t771 * t755 - t819 * t763 + t853;
t697 = t837 * t701 + t833 * t702;
t773 = -t803 * mrSges(5,1) + t804 * mrSges(5,2);
t781 = -mrSges(5,2) * t820 + t803 * mrSges(5,3);
t695 = m(5) * t741 + t811 * mrSges(5,1) - t784 * mrSges(5,3) - t804 * t773 + t781 * t820 + t697;
t782 = mrSges(5,1) * t820 - t804 * mrSges(5,3);
t854 = -t833 * t701 + t837 * t702;
t696 = m(5) * t742 - t811 * mrSges(5,2) + t783 * mrSges(5,3) + t803 * t773 - t782 * t820 + t854;
t855 = -t830 * t695 + t831 * t696;
t687 = m(4) * t776 - qJDD(3) * mrSges(4,2) - t811 * mrSges(4,3) - qJD(3) * t815 - t810 * t820 + t855;
t775 = t834 * g(3) + t838 * t786;
t760 = -qJDD(3) * pkin(3) - t840 * qJ(4) + t809 * t862 + qJDD(4) - t775;
t743 = -t783 * pkin(4) - t802 * pkin(8) + t804 * t785 + t760;
t722 = -t747 * pkin(5) - t770 * pkin(9) + t772 * t764 + t743;
t852 = m(7) * t722 - t727 * mrSges(7,1) + t728 * mrSges(7,2) - t753 * t744 + t754 * t745;
t844 = m(6) * t743 - t747 * mrSges(6,1) + t748 * mrSges(6,2) - t771 * t762 + t772 * t763 + t852;
t713 = m(5) * t760 - t783 * mrSges(5,1) + t784 * mrSges(5,2) - t803 * t781 + t804 * t782 + t844;
t814 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t820;
t706 = m(4) * t775 + qJDD(3) * mrSges(4,1) - t812 * mrSges(4,3) + qJD(3) * t814 - t810 * t862 - t713;
t681 = t834 * t687 + t838 * t706;
t791 = -qJDD(1) * pkin(1) + t850;
t849 = -m(3) * t791 + (t841 * mrSges(3,3)) - t681;
t677 = m(2) * t816 - (t841 * mrSges(2,2)) + qJDD(1) * t866 + t849;
t789 = t841 * pkin(1) + t868;
t689 = t831 * t695 + t830 * t696;
t848 = -m(4) * t780 - t811 * mrSges(4,1) - t812 * mrSges(4,2) - t814 * t820 - t815 * t862 - t689;
t845 = -m(3) * t789 + (t841 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t848;
t684 = m(2) * t817 - (t841 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t845;
t863 = t839 * t677 + t835 * t684;
t857 = -t835 * t677 + t839 * t684;
t856 = t838 * t687 - t834 * t706;
t736 = Ifges(7,4) * t754 + Ifges(7,2) * t753 + Ifges(7,6) * t818;
t737 = Ifges(7,1) * t754 + Ifges(7,4) * t753 + Ifges(7,5) * t818;
t847 = -mrSges(7,1) * t715 + mrSges(7,2) * t716 - Ifges(7,5) * t728 - Ifges(7,6) * t727 - Ifges(7,3) * t799 - t754 * t736 + t753 * t737;
t735 = Ifges(7,5) * t754 + Ifges(7,6) * t753 + Ifges(7,3) * t818;
t704 = -mrSges(7,1) * t722 + mrSges(7,3) * t716 + Ifges(7,4) * t728 + Ifges(7,2) * t727 + Ifges(7,6) * t799 - t754 * t735 + t818 * t737;
t705 = mrSges(7,2) * t722 - mrSges(7,3) * t715 + Ifges(7,1) * t728 + Ifges(7,4) * t727 + Ifges(7,5) * t799 + t753 * t735 - t818 * t736;
t749 = Ifges(6,5) * t772 + Ifges(6,6) * t771 + Ifges(6,3) * t819;
t751 = Ifges(6,1) * t772 + Ifges(6,4) * t771 + Ifges(6,5) * t819;
t690 = -mrSges(6,1) * t743 + mrSges(6,3) * t720 + Ifges(6,4) * t748 + Ifges(6,2) * t747 + Ifges(6,6) * t808 - pkin(5) * t852 + pkin(9) * t853 + t836 * t704 + t832 * t705 - t772 * t749 + t819 * t751;
t750 = Ifges(6,4) * t772 + Ifges(6,2) * t771 + Ifges(6,6) * t819;
t691 = mrSges(6,2) * t743 - mrSges(6,3) * t719 + Ifges(6,1) * t748 + Ifges(6,4) * t747 + Ifges(6,5) * t808 - pkin(9) * t703 - t832 * t704 + t836 * t705 + t771 * t749 - t819 * t750;
t765 = Ifges(5,5) * t804 + Ifges(5,6) * t803 + Ifges(5,3) * t820;
t767 = Ifges(5,1) * t804 + Ifges(5,4) * t803 + Ifges(5,5) * t820;
t673 = -mrSges(5,1) * t760 + mrSges(5,3) * t742 + Ifges(5,4) * t784 + Ifges(5,2) * t783 + Ifges(5,6) * t811 - pkin(4) * t844 + pkin(8) * t854 + t837 * t690 + t833 * t691 - t804 * t765 + t767 * t820;
t766 = Ifges(5,4) * t804 + Ifges(5,2) * t803 + Ifges(5,6) * t820;
t675 = mrSges(5,2) * t760 - mrSges(5,3) * t741 + Ifges(5,1) * t784 + Ifges(5,4) * t783 + Ifges(5,5) * t811 - pkin(8) * t697 - t833 * t690 + t837 * t691 + t803 * t765 - t766 * t820;
t797 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t838 - Ifges(4,2) * t834) * qJD(1);
t798 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t838 - Ifges(4,4) * t834) * qJD(1);
t846 = mrSges(4,1) * t775 - mrSges(4,2) * t776 + Ifges(4,5) * t812 - Ifges(4,6) * t811 + Ifges(4,3) * qJDD(3) - pkin(3) * t713 + qJ(4) * t855 + t831 * t673 + t830 * t675 + t797 * t862 + t798 * t820;
t796 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t838 - Ifges(4,6) * t834) * qJD(1);
t670 = mrSges(4,2) * t780 - mrSges(4,3) * t775 + Ifges(4,1) * t812 - Ifges(4,4) * t811 + Ifges(4,5) * qJDD(3) - qJ(4) * t689 - qJD(3) * t797 - t830 * t673 + t831 * t675 - t796 * t820;
t842 = mrSges(6,1) * t719 - mrSges(6,2) * t720 + Ifges(6,5) * t748 + Ifges(6,6) * t747 + Ifges(6,3) * t808 + pkin(5) * t703 + t772 * t750 - t771 * t751 - t847;
t671 = -t842 - t796 * t862 + Ifges(4,6) * qJDD(3) + (-Ifges(5,3) - Ifges(4,2)) * t811 + Ifges(4,4) * t812 + t803 * t767 - t804 * t766 - Ifges(5,5) * t784 + qJD(3) * t798 + mrSges(4,3) * t776 - mrSges(4,1) * t780 - Ifges(5,6) * t783 - mrSges(5,1) * t741 + mrSges(5,2) * t742 - pkin(4) * t697 - pkin(3) * t689;
t679 = qJDD(1) * mrSges(3,2) - t849;
t843 = mrSges(2,1) * t816 - mrSges(2,2) * t817 + mrSges(3,2) * t791 - mrSges(3,3) * t789 - pkin(1) * t679 - pkin(7) * t681 + qJ(2) * t845 + t838 * t670 - t834 * t671 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t680 = -m(3) * g(3) + t856;
t668 = t846 + (t864 * t841) - mrSges(2,3) * t816 + mrSges(3,1) * t791 + pkin(2) * t681 + t865 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - qJ(2) * t680;
t667 = -mrSges(3,1) * t789 + mrSges(2,3) * t817 - pkin(1) * t680 - pkin(2) * t848 - pkin(7) * t856 + g(3) * t866 - qJDD(1) * t864 - t834 * t670 - t838 * t671 + t841 * t865;
t1 = [-m(1) * g(1) + t857; -m(1) * g(2) + t863; (-m(1) - m(2) - m(3)) * g(3) + t856; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t863 - t835 * t667 + t839 * t668; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t857 + t839 * t667 + t835 * t668; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t843; t843; t679; t846; t713; t842; -t847;];
tauJB  = t1;
