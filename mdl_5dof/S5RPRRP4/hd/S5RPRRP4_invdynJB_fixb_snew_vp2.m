% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:43
% EndTime: 2020-01-03 11:49:49
% DurationCPUTime: 4.33s
% Computational Cost: add. (35957->268), mult. (86495->338), div. (0->0), fcn. (56474->8), ass. (0->118)
t816 = Ifges(5,4) + Ifges(6,4);
t826 = Ifges(5,2) + Ifges(6,2);
t821 = Ifges(5,6) + Ifges(6,6);
t775 = sin(qJ(1));
t778 = cos(qJ(1));
t754 = -t775 * g(2) + t778 * g(3);
t779 = qJD(1) ^ 2;
t825 = -t779 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t754;
t773 = sin(qJ(4));
t774 = sin(qJ(3));
t776 = cos(qJ(4));
t777 = cos(qJ(3));
t771 = sin(pkin(8));
t806 = t771 * qJD(1);
t732 = (-t773 * t777 - t774 * t776) * t806;
t803 = qJD(1) * qJD(3);
t741 = (-qJDD(1) * t774 - t777 * t803) * t771;
t742 = (qJDD(1) * t777 - t774 * t803) * t771;
t701 = t732 * qJD(4) + t773 * t741 + t776 * t742;
t733 = (-t773 * t774 + t776 * t777) * t806;
t714 = -t732 * mrSges(6,1) + t733 * mrSges(6,2);
t772 = cos(pkin(8));
t725 = -t771 * g(1) + t825 * t772;
t790 = -pkin(2) * t772 - pkin(6) * t771;
t748 = t790 * qJD(1);
t805 = t772 * qJD(1);
t713 = t748 * t805 + t725;
t755 = -t778 * g(2) - t775 * g(3);
t786 = -t779 * qJ(2) + qJDD(2) - t755;
t726 = (-pkin(1) + t790) * qJDD(1) + t786;
t723 = t777 * t726;
t802 = t772 * qJDD(1);
t757 = qJDD(3) - t802;
t758 = qJD(3) - t805;
t812 = t771 ^ 2 * t779;
t684 = t757 * pkin(3) - t742 * pkin(7) + t723 + (-pkin(3) * t777 * t812 - pkin(7) * t758 * t806 - t713) * t774;
t688 = t777 * t713 + t774 * t726;
t797 = t777 * t806;
t740 = t758 * pkin(3) - pkin(7) * t797;
t801 = t774 ^ 2 * t812;
t685 = -pkin(3) * t801 + t741 * pkin(7) - t758 * t740 + t688;
t677 = t776 * t684 - t773 * t685;
t753 = qJDD(4) + t757;
t756 = qJD(4) + t758;
t673 = -0.2e1 * qJD(5) * t733 + (t732 * t756 - t701) * qJ(5) + (t732 * t733 + t753) * pkin(4) + t677;
t717 = -t756 * mrSges(6,2) + t732 * mrSges(6,3);
t800 = m(6) * t673 + t753 * mrSges(6,1) + t756 * t717;
t669 = -t701 * mrSges(6,3) - t733 * t714 + t800;
t678 = t773 * t684 + t776 * t685;
t700 = -t733 * qJD(4) + t776 * t741 - t773 * t742;
t719 = t756 * pkin(4) - t733 * qJ(5);
t731 = t732 ^ 2;
t675 = -t731 * pkin(4) + t700 * qJ(5) + 0.2e1 * qJD(5) * t732 - t756 * t719 + t678;
t822 = Ifges(5,5) + Ifges(6,5);
t823 = Ifges(5,1) + Ifges(6,1);
t809 = -t816 * t732 - t823 * t733 - t822 * t756;
t819 = t826 * t732 + t816 * t733 + t821 * t756;
t820 = Ifges(5,3) + Ifges(6,3);
t824 = mrSges(5,1) * t677 + mrSges(6,1) * t673 - mrSges(5,2) * t678 - mrSges(6,2) * t675 + pkin(4) * t669 + t821 * t700 + t822 * t701 + t809 * t732 + t819 * t733 + t820 * t753;
t715 = -t732 * mrSges(5,1) + t733 * mrSges(5,2);
t718 = -t756 * mrSges(5,2) + t732 * mrSges(5,3);
t662 = m(5) * t677 + t753 * mrSges(5,1) + t756 * t718 + (-t714 - t715) * t733 + (-mrSges(5,3) - mrSges(6,3)) * t701 + t800;
t720 = t756 * mrSges(6,1) - t733 * mrSges(6,3);
t721 = t756 * mrSges(5,1) - t733 * mrSges(5,3);
t799 = m(6) * t675 + t700 * mrSges(6,3) + t732 * t714;
t665 = m(5) * t678 + t700 * mrSges(5,3) + t732 * t715 + (-t720 - t721) * t756 + (-mrSges(5,2) - mrSges(6,2)) * t753 + t799;
t660 = t776 * t662 + t773 * t665;
t687 = -t774 * t713 + t723;
t818 = -mrSges(4,1) * t687 + mrSges(4,2) * t688 - Ifges(4,5) * t742 - Ifges(4,6) * t741 - Ifges(4,3) * t757 - pkin(3) * t660 - t824;
t724 = -t772 * g(1) - t825 * t771;
t798 = t774 * t806;
t737 = -t758 * mrSges(4,2) - mrSges(4,3) * t798;
t738 = t758 * mrSges(4,1) - mrSges(4,3) * t797;
t817 = -t737 * t774 - t738 * t777;
t815 = mrSges(3,2) * t771;
t746 = (-mrSges(3,1) * t772 + t815) * qJD(1);
t739 = (mrSges(4,1) * t774 + mrSges(4,2) * t777) * t806;
t657 = m(4) * t687 + t757 * mrSges(4,1) - t742 * mrSges(4,3) + t758 * t737 - t739 * t797 + t660;
t791 = -t773 * t662 + t776 * t665;
t658 = m(4) * t688 - t757 * mrSges(4,2) + t741 * mrSges(4,3) - t758 * t738 - t739 * t798 + t791;
t792 = -t774 * t657 + t777 * t658;
t807 = qJDD(1) * mrSges(3,3);
t651 = m(3) * t725 + (qJD(1) * t746 + t807) * t772 + t792;
t712 = t748 * t806 - t724;
t686 = -t741 * pkin(3) - pkin(7) * t801 + t740 * t797 + t712;
t680 = -t700 * pkin(4) - t731 * qJ(5) + t733 * t719 + qJDD(5) + t686;
t670 = m(6) * t680 - t700 * mrSges(6,1) + t701 * mrSges(6,2) - t732 * t717 + t733 * t720;
t783 = m(5) * t686 - t700 * mrSges(5,1) + t701 * mrSges(5,2) - t732 * t718 + t733 * t721 + t670;
t781 = -m(4) * t712 + t741 * mrSges(4,1) - t742 * mrSges(4,2) - t783;
t667 = t781 + m(3) * t724 + (-t807 + (-t746 + t817) * qJD(1)) * t771;
t793 = t772 * t651 - t771 * t667;
t643 = m(2) * t754 - t779 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t793;
t654 = t777 * t657 + t774 * t658;
t744 = -qJDD(1) * pkin(1) + t786;
t784 = -m(3) * t744 + mrSges(3,1) * t802 - t654 + (t772 ^ 2 * t779 + t812) * mrSges(3,3);
t648 = m(2) * t755 - t779 * mrSges(2,2) + (mrSges(2,1) - t815) * qJDD(1) + t784;
t811 = t775 * t643 + t778 * t648;
t645 = t771 * t651 + t772 * t667;
t810 = -t821 * t732 - t822 * t733 - t820 * t756;
t794 = -t778 * t643 + t775 * t648;
t789 = Ifges(3,1) * t771 + Ifges(3,4) * t772;
t788 = Ifges(3,5) * t771 + Ifges(3,6) * t772;
t728 = Ifges(4,6) * t758 + (Ifges(4,4) * t777 - Ifges(4,2) * t774) * t806;
t729 = Ifges(4,5) * t758 + (Ifges(4,1) * t777 - Ifges(4,4) * t774) * t806;
t787 = t728 * t777 + t729 * t774;
t655 = -mrSges(5,1) * t686 + mrSges(5,3) * t678 - mrSges(6,1) * t680 + mrSges(6,3) * t675 - pkin(4) * t670 + qJ(5) * t799 + (-qJ(5) * t720 - t809) * t756 + (-qJ(5) * mrSges(6,2) + t821) * t753 + t810 * t733 + t816 * t701 + t826 * t700;
t659 = mrSges(5,2) * t686 + mrSges(6,2) * t680 - mrSges(5,3) * t677 - mrSges(6,3) * t673 - qJ(5) * t669 + t816 * t700 + t823 * t701 - t810 * t732 + t822 * t753 - t819 * t756;
t727 = Ifges(4,3) * t758 + (Ifges(4,5) * t777 - Ifges(4,6) * t774) * t806;
t640 = -mrSges(4,1) * t712 + mrSges(4,3) * t688 + Ifges(4,4) * t742 + Ifges(4,2) * t741 + Ifges(4,6) * t757 - pkin(3) * t783 + pkin(7) * t791 + t776 * t655 + t773 * t659 - t727 * t797 + t758 * t729;
t641 = mrSges(4,2) * t712 - mrSges(4,3) * t687 + Ifges(4,1) * t742 + Ifges(4,4) * t741 + Ifges(4,5) * t757 - pkin(7) * t660 - t773 * t655 + t776 * t659 - t727 * t798 - t758 * t728;
t747 = t788 * qJD(1);
t637 = mrSges(3,2) * t744 - mrSges(3,3) * t724 - pkin(6) * t654 + t789 * qJDD(1) - t774 * t640 + t777 * t641 + t747 * t805;
t639 = mrSges(3,3) * t725 - mrSges(3,1) * t744 - pkin(2) * t654 + Ifges(3,2) * t802 + (Ifges(3,4) * qJDD(1) + (-t747 - t787) * qJD(1)) * t771 + t818;
t653 = qJDD(1) * t815 - t784;
t785 = mrSges(2,1) * t755 - mrSges(2,2) * t754 + Ifges(2,3) * qJDD(1) - pkin(1) * t653 + qJ(2) * t793 + t771 * t637 + t772 * t639;
t635 = t779 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t754 - mrSges(3,1) * t724 + mrSges(3,2) * t725 - t774 * t641 - t777 * t640 - pkin(2) * t781 - pkin(6) * t792 - pkin(1) * t645 + (Ifges(2,6) - t788) * qJDD(1) + (-pkin(2) * t817 * t771 + (-t771 * (Ifges(3,4) * t771 + Ifges(3,2) * t772) + t772 * t789) * qJD(1)) * qJD(1);
t634 = -mrSges(2,2) * g(1) - mrSges(2,3) * t755 + Ifges(2,5) * qJDD(1) - t779 * Ifges(2,6) - qJ(2) * t645 + t772 * t637 - t771 * t639;
t1 = [(-m(1) - m(2)) * g(1) + t645; -m(1) * g(2) + t811; -m(1) * g(3) + t794; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t785; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t794 + t775 * t634 + t778 * t635; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t811 - t778 * t634 + t775 * t635; t785; t653; t787 * t806 - t818; t824; t670;];
tauJB = t1;
