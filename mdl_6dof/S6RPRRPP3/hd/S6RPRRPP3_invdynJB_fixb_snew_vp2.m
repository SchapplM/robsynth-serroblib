% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:25:46
% EndTime: 2019-05-05 21:25:51
% DurationCPUTime: 4.03s
% Computational Cost: add. (44714->298), mult. (84077->342), div. (0->0), fcn. (48977->8), ass. (0->122)
t858 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t841 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t840 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t857 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t839 = -Ifges(5,6) + Ifges(6,5) + Ifges(7,4);
t856 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t813 = sin(qJ(1));
t815 = cos(qJ(1));
t796 = t813 * g(1) - g(2) * t815;
t787 = qJDD(1) * pkin(1) + t796;
t797 = -g(1) * t815 - g(2) * t813;
t817 = qJD(1) ^ 2;
t789 = -pkin(1) * t817 + t797;
t809 = sin(pkin(9));
t810 = cos(pkin(9));
t753 = t787 * t810 - t809 * t789;
t730 = -qJDD(1) * pkin(2) - pkin(7) * t817 - t753;
t812 = sin(qJ(3));
t814 = cos(qJ(3));
t842 = qJD(1) * qJD(3);
t832 = t814 * t842;
t791 = qJDD(1) * t812 + t832;
t833 = t812 * t842;
t792 = qJDD(1) * t814 - t833;
t719 = (-t791 - t832) * pkin(8) + (-t792 + t833) * pkin(3) + t730;
t754 = t809 * t787 + t810 * t789;
t731 = -pkin(2) * t817 + qJDD(1) * pkin(7) + t754;
t808 = -g(3) + qJDD(2);
t725 = t814 * t731 + t812 * t808;
t790 = (-pkin(3) * t814 - pkin(8) * t812) * qJD(1);
t816 = qJD(3) ^ 2;
t843 = qJD(1) * t814;
t723 = -pkin(3) * t816 + qJDD(3) * pkin(8) + t790 * t843 + t725;
t811 = sin(qJ(4));
t850 = cos(qJ(4));
t716 = t719 * t850 - t811 * t723;
t844 = qJD(1) * t812;
t785 = -qJD(3) * t850 + t811 * t844;
t786 = t811 * qJD(3) + t844 * t850;
t759 = pkin(4) * t785 - qJ(5) * t786;
t784 = qJDD(4) - t792;
t799 = -qJD(4) + t843;
t798 = t799 ^ 2;
t714 = -t784 * pkin(4) - t798 * qJ(5) + t786 * t759 + qJDD(5) - t716;
t751 = -t785 * qJD(4) + t811 * qJDD(3) + t791 * t850;
t761 = -mrSges(6,2) * t785 - mrSges(6,3) * t786;
t855 = -m(6) * t714 - t751 * mrSges(6,1) - t786 * t761;
t768 = mrSges(6,1) * t786 - mrSges(6,2) * t799;
t750 = qJD(4) * t786 - qJDD(3) * t850 + t791 * t811;
t724 = -t812 * t731 + t808 * t814;
t722 = -qJDD(3) * pkin(3) - pkin(8) * t816 + t790 * t844 - t724;
t847 = t785 * t799;
t852 = -2 * qJD(5);
t820 = (-t751 - t847) * qJ(5) + t722 + (-t799 * pkin(4) + t852) * t786;
t715 = pkin(4) * t750 + t820;
t766 = mrSges(6,1) * t785 + mrSges(6,3) * t799;
t764 = pkin(5) * t786 + qJ(6) * t799;
t783 = t785 ^ 2;
t851 = 2 * qJD(6);
t712 = -pkin(5) * t783 + t785 * t851 - t764 * t786 + (pkin(4) + qJ(6)) * t750 + t820;
t765 = mrSges(7,1) * t786 + mrSges(7,3) * t799;
t767 = -mrSges(7,1) * t785 - mrSges(7,2) * t799;
t827 = m(7) * t712 - t751 * mrSges(7,2) + t750 * mrSges(7,3) - t786 * t765 + t785 * t767;
t824 = -m(6) * t715 + t750 * mrSges(6,2) + t785 * t766 - t827;
t704 = -mrSges(6,3) * t751 - t768 * t786 - t824;
t758 = -mrSges(7,2) * t786 + mrSges(7,3) * t785;
t717 = t811 * t719 + t850 * t723;
t823 = -pkin(4) * t798 + qJ(5) * t784 - t759 * t785 + t717;
t710 = -pkin(5) * t750 - qJ(6) * t783 + qJDD(6) + (t852 - t764) * t799 + t823;
t837 = m(7) * t710 + t784 * mrSges(7,2) - t799 * t765;
t707 = -mrSges(7,1) * t750 - t758 * t785 + t837;
t713 = 0.2e1 * qJD(5) * t799 - t823;
t834 = -t841 * t785 + t858 * t786 - t840 * t799;
t836 = -t839 * t785 - t840 * t786 + t856 * t799;
t681 = -mrSges(5,1) * t722 - mrSges(6,1) * t713 + mrSges(7,1) * t710 + mrSges(6,2) * t715 + mrSges(5,3) * t717 - mrSges(7,3) * t712 - pkin(4) * t704 + pkin(5) * t707 - qJ(6) * t827 + t857 * t750 + t841 * t751 - t839 * t784 + t836 * t786 - t834 * t799;
t708 = t799 * t851 + (t785 * t786 - t784) * qJ(6) + (t751 - t847) * pkin(5) + t714;
t828 = -m(7) * t708 + t784 * mrSges(7,3) - t799 * t767;
t706 = mrSges(7,1) * t751 + t758 * t786 - t828;
t835 = t857 * t785 + t841 * t786 + t839 * t799;
t688 = mrSges(6,1) * t714 + mrSges(7,1) * t708 + mrSges(5,2) * t722 - mrSges(7,2) * t712 - mrSges(5,3) * t716 - mrSges(6,3) * t715 + pkin(5) * t706 - qJ(5) * t704 - t841 * t750 + t858 * t751 + t840 * t784 + t836 * t785 + t835 * t799;
t760 = mrSges(5,1) * t785 + mrSges(5,2) * t786;
t762 = mrSges(5,2) * t799 - mrSges(5,3) * t785;
t848 = -mrSges(7,1) - mrSges(5,3);
t699 = m(5) * t716 + (-t762 + t766) * t799 + (-t758 - t760) * t786 + (mrSges(5,1) - mrSges(6,2)) * t784 + t848 * t751 + t828 + t855;
t763 = -mrSges(5,1) * t799 - mrSges(5,3) * t786;
t825 = -m(6) * t713 + t784 * mrSges(6,3) - t799 * t768 + t837;
t845 = -t758 - t761;
t701 = m(5) * t717 - mrSges(5,2) * t784 + t763 * t799 + (-t760 + t845) * t785 + (-mrSges(6,1) + t848) * t750 + t825;
t696 = -t699 * t811 + t850 * t701;
t702 = -m(5) * t722 - t750 * mrSges(5,1) - t785 * t762 + (-t763 + t768) * t786 + (-mrSges(5,2) + mrSges(6,3)) * t751 + t824;
t776 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t812 + Ifges(4,2) * t814) * qJD(1);
t777 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t812 + Ifges(4,4) * t814) * qJD(1);
t854 = mrSges(4,1) * t724 - mrSges(4,2) * t725 + Ifges(4,5) * t791 + Ifges(4,6) * t792 + Ifges(4,3) * qJDD(3) + pkin(3) * t702 + pkin(8) * t696 + (t776 * t812 - t777 * t814) * qJD(1) + t681 * t850 + t811 * t688;
t703 = mrSges(6,2) * t784 - t766 * t799 + t706 - t855;
t853 = t750 * t839 + t751 * t840 + t856 * t784 + t785 * t834 + t786 * t835 + mrSges(5,1) * t716 - mrSges(5,2) * t717 + mrSges(6,2) * t714 + mrSges(7,2) * t710 - mrSges(6,3) * t713 - mrSges(7,3) * t708 - pkin(4) * t703 + qJ(5) * (t845 * t785 + (-mrSges(6,1) - mrSges(7,1)) * t750 + t825) - qJ(6) * t706;
t788 = (-mrSges(4,1) * t814 + mrSges(4,2) * t812) * qJD(1);
t794 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t844;
t694 = m(4) * t725 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t792 - qJD(3) * t794 + t788 * t843 + t696;
t795 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t843;
t698 = m(4) * t724 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t791 + qJD(3) * t795 - t788 * t844 + t702;
t829 = t814 * t694 - t698 * t812;
t684 = m(3) * t754 - mrSges(3,1) * t817 - qJDD(1) * mrSges(3,2) + t829;
t695 = t699 * t850 + t811 * t701;
t821 = -m(4) * t730 + t792 * mrSges(4,1) - t791 * mrSges(4,2) - t794 * t844 + t795 * t843 - t695;
t690 = m(3) * t753 + qJDD(1) * mrSges(3,1) - t817 * mrSges(3,2) + t821;
t680 = t809 * t684 + t810 * t690;
t677 = m(2) * t796 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t817 + t680;
t830 = t810 * t684 - t690 * t809;
t678 = m(2) * t797 - mrSges(2,1) * t817 - qJDD(1) * mrSges(2,2) + t830;
t846 = t815 * t677 + t813 * t678;
t687 = t812 * t694 + t814 * t698;
t685 = m(3) * t808 + t687;
t831 = -t677 * t813 + t815 * t678;
t775 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t812 + Ifges(4,6) * t814) * qJD(1);
t671 = mrSges(4,2) * t730 - mrSges(4,3) * t724 + Ifges(4,1) * t791 + Ifges(4,4) * t792 + Ifges(4,5) * qJDD(3) - pkin(8) * t695 - qJD(3) * t776 - t811 * t681 + t688 * t850 + t775 * t843;
t673 = -mrSges(4,1) * t730 + mrSges(4,3) * t725 + Ifges(4,4) * t791 + Ifges(4,2) * t792 + Ifges(4,6) * qJDD(3) - pkin(3) * t695 + qJD(3) * t777 - t775 * t844 - t853;
t822 = mrSges(2,1) * t796 + mrSges(3,1) * t753 - mrSges(2,2) * t797 - mrSges(3,2) * t754 + pkin(1) * t680 + pkin(2) * t821 + pkin(7) * t829 + t812 * t671 + t814 * t673 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t669 = -mrSges(3,1) * t808 + mrSges(3,3) * t754 + t817 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t687 - t854;
t668 = mrSges(3,2) * t808 - mrSges(3,3) * t753 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t817 - pkin(7) * t687 + t671 * t814 - t673 * t812;
t667 = -mrSges(2,2) * g(3) - mrSges(2,3) * t796 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t817 - qJ(2) * t680 + t668 * t810 - t669 * t809;
t666 = mrSges(2,1) * g(3) + mrSges(2,3) * t797 + t817 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t685 + qJ(2) * t830 + t809 * t668 + t810 * t669;
t1 = [-m(1) * g(1) + t831; -m(1) * g(2) + t846; (-m(1) - m(2)) * g(3) + t685; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t846 - t813 * t666 + t815 * t667; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t831 + t815 * t666 + t813 * t667; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t822; t822; t685; t854; t853; t703; t707;];
tauJB  = t1;
