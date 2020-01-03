% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:51
% EndTime: 2019-12-31 19:23:58
% DurationCPUTime: 5.45s
% Computational Cost: add. (51743->292), mult. (130073->356), div. (0->0), fcn. (90418->8), ass. (0->124)
t858 = -2 * qJD(3);
t857 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t856 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t832 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t855 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t831 = Ifges(4,6) - Ifges(5,5) - Ifges(6,4);
t854 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t801 = cos(qJ(2));
t845 = cos(pkin(5));
t820 = qJD(1) * t845;
t798 = sin(pkin(5));
t836 = qJD(2) * t798;
t770 = (t801 * t820 + t836) * qJ(3);
t800 = sin(qJ(1));
t802 = cos(qJ(1));
t790 = g(1) * t800 - t802 * g(2);
t803 = qJD(1) ^ 2;
t774 = -qJDD(1) * pkin(1) - pkin(7) * t803 - t790;
t799 = sin(qJ(2));
t843 = qJ(3) * t799;
t782 = qJD(2) * pkin(2) - t820 * t843;
t834 = qJD(1) * qJD(2);
t784 = qJDD(1) * t799 + t801 * t834;
t785 = qJDD(1) * t801 - t799 * t834;
t722 = -qJ(3) * t784 * t798 - pkin(2) * t785 + (-t770 * t801 + t782 * t799) * qJD(1) + t774;
t791 = -g(1) * t802 - g(2) * t800;
t775 = -pkin(1) * t803 + qJDD(1) * pkin(7) + t791;
t776 = (-pkin(2) * t801 - t798 * t843) * qJD(1);
t824 = qJ(3) * t845;
t848 = t801 * g(3);
t723 = -t784 * t824 + qJDD(2) * pkin(2) - t848 + qJD(2) * t770 + (-qJD(1) * t776 - t775) * t799;
t763 = -t799 * g(3) + t801 * t775;
t812 = qJDD(2) * t798 + t845 * t785;
t837 = qJD(1) * t801;
t724 = t812 * qJ(3) - qJD(2) * t782 + t776 * t837 + t763;
t797 = sin(pkin(8));
t823 = t797 * t845;
t844 = cos(pkin(8));
t761 = t797 * t836 + (t844 * t799 + t801 * t823) * qJD(1);
t816 = t845 * t844;
t822 = t798 * t844;
t712 = t722 * t822 + t723 * t816 - t797 * t724 + t761 * t858;
t838 = qJD(1) * t799;
t760 = -qJD(2) * t822 + t797 * t838 - t816 * t837;
t740 = pkin(3) * t760 - qJ(4) * t761;
t767 = t845 * qJDD(2) - t798 * t785;
t779 = -t845 * qJD(2) + t798 * t837;
t778 = t779 ^ 2;
t709 = -t767 * pkin(3) - t778 * qJ(4) + t761 * t740 + qJDD(4) - t712;
t742 = -mrSges(5,2) * t760 - mrSges(5,3) * t761;
t754 = t844 * t784 + t812 * t797;
t853 = -m(5) * t709 - t754 * mrSges(5,1) - t761 * t742;
t739 = -mrSges(6,2) * t761 + mrSges(6,3) * t760;
t741 = mrSges(4,1) * t760 + mrSges(4,2) * t761;
t842 = t760 * t779;
t850 = 2 * qJD(5);
t703 = t779 * t850 + (t760 * t761 - t767) * qJ(5) + (t754 - t842) * pkin(4) + t709;
t748 = -mrSges(6,1) * t760 - mrSges(6,2) * t779;
t817 = -m(6) * t703 + t767 * mrSges(6,3) - t779 * t748;
t747 = mrSges(5,1) * t760 + mrSges(5,3) * t779;
t839 = mrSges(4,2) * t779 - mrSges(4,3) * t760 - t747;
t846 = -mrSges(6,1) - mrSges(4,3);
t847 = mrSges(4,1) - mrSges(5,2);
t694 = m(4) * t712 - t839 * t779 + t847 * t767 + (-t739 - t741) * t761 + t846 * t754 + t817 + t853;
t757 = t760 * t858;
t828 = t798 * t797 * t722 + t723 * t823 + t844 * t724;
t713 = t757 + t828;
t744 = -mrSges(4,1) * t779 - mrSges(4,3) * t761;
t753 = -qJDD(2) * t822 + t784 * t797 - t785 * t816;
t810 = pkin(3) * t778 - qJ(4) * t767 - t828;
t708 = 0.2e1 * qJD(4) * t779 + ((2 * qJD(3)) + t740) * t760 + t810;
t749 = mrSges(5,1) * t761 - mrSges(5,2) * t779;
t745 = pkin(4) * t761 + qJ(5) * t779;
t759 = t760 ^ 2;
t851 = -0.2e1 * qJD(4);
t705 = -pkin(4) * t753 - qJ(5) * t759 - t740 * t760 + qJDD(5) + t757 + (t851 - t745) * t779 - t810;
t746 = mrSges(6,1) * t761 + mrSges(6,3) * t779;
t830 = m(6) * t705 + t767 * mrSges(6,2) - t779 * t746;
t811 = -m(5) * t708 + t767 * mrSges(5,3) - t779 * t749 + t830;
t840 = -t739 - t742;
t697 = m(4) * t713 - mrSges(4,2) * t767 + t744 * t779 + (-t741 + t840) * t760 + (-mrSges(5,1) + t846) * t753 + t811;
t714 = t845 * t722 - t723 * t798 + qJDD(3);
t806 = (-t754 - t842) * qJ(4) + t714 + (-t779 * pkin(3) + t851) * t761;
t711 = pkin(3) * t753 + t806;
t707 = -pkin(4) * t759 + t760 * t850 - t745 * t761 + (pkin(3) + qJ(5)) * t753 + t806;
t829 = m(6) * t707 + t753 * mrSges(6,3) + t760 * t748;
t814 = m(5) * t711 - t754 * mrSges(5,3) - t761 * t749 + t829;
t698 = m(4) * t714 + (t744 - t746) * t761 + t839 * t760 + (mrSges(4,2) - mrSges(6,2)) * t754 + t847 * t753 + t814;
t686 = (t844 * t694 + t697 * t797) * t798 + t845 * t698;
t701 = mrSges(6,1) * t754 + t739 * t761 - t817;
t699 = mrSges(5,2) * t767 - t747 * t779 + t701 - t853;
t825 = t856 * t760 - t857 * t761 + t832 * t779;
t826 = t855 * t760 + t856 * t761 - t831 * t779;
t679 = mrSges(4,1) * t712 - mrSges(4,2) * t713 + mrSges(5,2) * t709 - mrSges(5,3) * t708 + mrSges(6,2) * t705 - mrSges(6,3) * t703 - qJ(5) * t701 - pkin(3) * t699 + qJ(4) * t811 + t854 * t767 + t826 * t761 + t832 * t754 + (qJ(4) * t840 - t825) * t760 + (qJ(4) * (-mrSges(5,1) - mrSges(6,1)) - t831) * t753;
t700 = -mrSges(5,2) * t753 - mrSges(6,2) * t754 - t746 * t761 - t747 * t760 + t814;
t702 = -mrSges(6,1) * t753 - t739 * t760 + t830;
t827 = t831 * t760 - t832 * t761 + t854 * t779;
t680 = -mrSges(4,1) * t714 + mrSges(4,3) * t713 - mrSges(5,1) * t708 + mrSges(5,2) * t711 + mrSges(6,1) * t705 - mrSges(6,3) * t707 + pkin(4) * t702 - qJ(5) * t829 - pkin(3) * t700 + t825 * t779 + t831 * t767 + (qJ(5) * t746 + t827) * t761 + (qJ(5) * mrSges(6,2) + t856) * t754 + t855 * t753;
t687 = t694 * t816 + t697 * t823 - t698 * t798;
t688 = mrSges(5,1) * t709 + mrSges(6,1) * t703 + mrSges(4,2) * t714 - mrSges(6,2) * t707 - mrSges(4,3) * t712 - mrSges(5,3) * t711 + pkin(4) * t701 - qJ(4) * t700 - t753 * t856 + t857 * t754 + t827 * t760 + t832 * t767 + t826 * t779;
t692 = -t694 * t797 + t844 * t697;
t762 = -t799 * t775 - t848;
t772 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t799 + Ifges(3,2) * t801) * qJD(1);
t773 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t799 + Ifges(3,4) * t801) * qJD(1);
t852 = mrSges(3,1) * t762 - mrSges(3,2) * t763 + Ifges(3,5) * t784 + Ifges(3,6) * t785 + Ifges(3,3) * qJDD(2) + pkin(2) * t687 + (t772 * t799 - t773 * t801) * qJD(1) + t845 * t679 + (qJ(3) * t692 + t844 * t680 + t688 * t797) * t798;
t783 = (-mrSges(3,1) * t801 + mrSges(3,2) * t799) * qJD(1);
t789 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t837;
t685 = m(3) * t762 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t784 + qJD(2) * t789 - t783 * t838 + t687;
t788 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t838;
t691 = m(3) * t763 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t785 - qJD(2) * t788 + t783 * t837 + t692;
t818 = -t685 * t799 + t801 * t691;
t676 = m(2) * t791 - mrSges(2,1) * t803 - qJDD(1) * mrSges(2,2) + t818;
t804 = -m(3) * t774 + t785 * mrSges(3,1) - t784 * mrSges(3,2) - t788 * t838 + t789 * t837 - t686;
t682 = m(2) * t790 + qJDD(1) * mrSges(2,1) - t803 * mrSges(2,2) + t804;
t841 = t800 * t676 + t802 * t682;
t678 = t801 * t685 + t799 * t691;
t819 = t802 * t676 - t682 * t800;
t771 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t799 + Ifges(3,6) * t801) * qJD(1);
t671 = -mrSges(3,1) * t774 + mrSges(3,3) * t763 + Ifges(3,4) * t784 + Ifges(3,2) * t785 + Ifges(3,6) * qJDD(2) - pkin(2) * t686 + qJD(2) * t773 - t798 * t679 + t680 * t816 + t688 * t823 + t692 * t824 - t771 * t838;
t673 = t771 * t837 + mrSges(3,2) * t774 - mrSges(3,3) * t762 + t844 * t688 + Ifges(3,1) * t784 + Ifges(3,4) * t785 + Ifges(3,5) * qJDD(2) - qJD(2) * t772 - t797 * t680 + (-t686 * t798 - t845 * t687) * qJ(3);
t809 = mrSges(2,1) * t790 - mrSges(2,2) * t791 + Ifges(2,3) * qJDD(1) + pkin(1) * t804 + pkin(7) * t818 + t801 * t671 + t799 * t673;
t669 = mrSges(2,1) * g(3) + mrSges(2,3) * t791 + t803 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t678 - t852;
t668 = -mrSges(2,2) * g(3) - mrSges(2,3) * t790 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t803 - pkin(7) * t678 - t671 * t799 + t673 * t801;
t1 = [-m(1) * g(1) + t819; -m(1) * g(2) + t841; (-m(1) - m(2)) * g(3) + t678; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t841 + t802 * t668 - t800 * t669; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t819 + t800 * t668 + t802 * t669; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t809; t809; t852; t698; t699; t702;];
tauJB = t1;
