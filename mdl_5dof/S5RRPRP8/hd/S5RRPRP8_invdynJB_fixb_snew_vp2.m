% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:34
% EndTime: 2019-12-31 20:03:37
% DurationCPUTime: 2.47s
% Computational Cost: add. (17676->268), mult. (37154->314), div. (0->0), fcn. (20151->6), ass. (0->108)
t826 = Ifges(5,4) + Ifges(6,4);
t842 = Ifges(5,2) + Ifges(6,2);
t836 = Ifges(5,6) + Ifges(6,6);
t789 = sin(qJ(4));
t790 = sin(qJ(2));
t792 = cos(qJ(4));
t793 = cos(qJ(2));
t743 = (-t789 * t790 - t792 * t793) * qJD(1);
t813 = qJD(1) * qJD(2);
t809 = t793 * t813;
t758 = t790 * qJDD(1) + t809;
t810 = t790 * t813;
t759 = t793 * qJDD(1) - t810;
t710 = t743 * qJD(4) + t792 * t758 - t789 * t759;
t744 = (-t789 * t793 + t790 * t792) * qJD(1);
t721 = -t743 * mrSges(6,1) + t744 * mrSges(6,2);
t791 = sin(qJ(1));
t794 = cos(qJ(1));
t769 = -t794 * g(1) - t791 * g(2);
t796 = qJD(1) ^ 2;
t746 = -t796 * pkin(1) + qJDD(1) * pkin(6) + t769;
t725 = -t790 * g(3) + t793 * t746;
t755 = (-pkin(2) * t793 - qJ(3) * t790) * qJD(1);
t795 = qJD(2) ^ 2;
t814 = qJD(1) * t793;
t829 = 2 * qJD(3);
t698 = -t795 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t829 + t755 * t814 + t725;
t815 = qJD(1) * t790;
t767 = -qJD(2) * pkin(3) - pkin(7) * t815;
t822 = t793 ^ 2 * t796;
t691 = -pkin(3) * t822 - t759 * pkin(7) + qJD(2) * t767 + t698;
t724 = -t793 * g(3) - t790 * t746;
t708 = -qJDD(2) * pkin(2) - t795 * qJ(3) + t755 * t815 + qJDD(3) - t724;
t692 = (-t758 + t809) * pkin(7) + (-t790 * t793 * t796 - qJDD(2)) * pkin(3) + t708;
t684 = -t789 * t691 + t792 * t692;
t780 = -qJDD(2) + qJDD(4);
t781 = -qJD(2) + qJD(4);
t678 = -0.2e1 * qJD(5) * t744 + (t743 * t781 - t710) * qJ(5) + (t743 * t744 + t780) * pkin(4) + t684;
t726 = -t781 * mrSges(6,2) + t743 * mrSges(6,3);
t812 = m(6) * t678 + t780 * mrSges(6,1) + t781 * t726;
t674 = -t710 * mrSges(6,3) - t744 * t721 + t812;
t685 = t792 * t691 + t789 * t692;
t709 = -t744 * qJD(4) - t789 * t758 - t792 * t759;
t728 = t781 * pkin(4) - t744 * qJ(5);
t736 = t743 ^ 2;
t680 = -t736 * pkin(4) + t709 * qJ(5) + 0.2e1 * qJD(5) * t743 - t781 * t728 + t685;
t838 = Ifges(5,5) + Ifges(6,5);
t839 = Ifges(5,1) + Ifges(6,1);
t819 = t826 * t743 + t839 * t744 + t838 * t781;
t833 = t842 * t743 + t826 * t744 + t836 * t781;
t834 = Ifges(5,3) + Ifges(6,3);
t841 = mrSges(5,1) * t684 + mrSges(6,1) * t678 - mrSges(5,2) * t685 - mrSges(6,2) * t680 + pkin(4) * t674 + t836 * t709 + t838 * t710 - t819 * t743 + t833 * t744 + t834 * t780;
t840 = Ifges(3,1) + Ifges(4,1);
t827 = Ifges(3,4) - Ifges(4,5);
t825 = Ifges(3,5) + Ifges(4,4);
t837 = Ifges(3,2) + Ifges(4,3);
t824 = Ifges(3,6) - Ifges(4,6);
t835 = Ifges(3,3) + Ifges(4,2);
t756 = (-mrSges(4,1) * t793 - mrSges(4,3) * t790) * qJD(1);
t722 = -t743 * mrSges(5,1) + t744 * mrSges(5,2);
t727 = -t781 * mrSges(5,2) + t743 * mrSges(5,3);
t670 = m(5) * t684 + t780 * mrSges(5,1) + t781 * t727 + (-t721 - t722) * t744 + (-mrSges(5,3) - mrSges(6,3)) * t710 + t812;
t729 = t781 * mrSges(6,1) - t744 * mrSges(6,3);
t730 = t781 * mrSges(5,1) - t744 * mrSges(5,3);
t811 = m(6) * t680 + t709 * mrSges(6,3) + t743 * t721;
t672 = m(5) * t685 + t709 * mrSges(5,3) + t743 * t722 + (-t729 - t730) * t781 + (-mrSges(5,2) - mrSges(6,2)) * t780 + t811;
t665 = t792 * t670 + t789 * t672;
t766 = mrSges(4,2) * t814 + qJD(2) * mrSges(4,3);
t801 = -m(4) * t708 + qJDD(2) * mrSges(4,1) + qJD(2) * t766 - t665;
t664 = t758 * mrSges(4,2) + t756 * t815 - t801;
t764 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t815;
t805 = -t789 * t670 + t792 * t672;
t803 = m(4) * t698 + qJDD(2) * mrSges(4,3) + qJD(2) * t764 + t756 * t814 + t805;
t816 = t825 * qJD(2) + (t840 * t790 + t827 * t793) * qJD(1);
t817 = -t824 * qJD(2) + (-t827 * t790 - t837 * t793) * qJD(1);
t832 = -(t817 * t790 + t816 * t793) * qJD(1) + t835 * qJDD(2) + t825 * t758 + t824 * t759 + mrSges(3,1) * t724 - mrSges(4,1) * t708 - mrSges(3,2) * t725 + mrSges(4,3) * t698 - pkin(2) * t664 - pkin(3) * t665 + qJ(3) * (t759 * mrSges(4,2) + t803) - t841;
t828 = mrSges(3,3) + mrSges(4,2);
t757 = (-mrSges(3,1) * t793 + mrSges(3,2) * t790) * qJD(1);
t763 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t815;
t660 = m(3) * t725 - qJDD(2) * mrSges(3,2) - qJD(2) * t763 + t757 * t814 + t828 * t759 + t803;
t765 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t814;
t661 = m(3) * t724 + qJDD(2) * mrSges(3,1) + qJD(2) * t765 - t828 * t758 + (-t756 - t757) * t815 + t801;
t806 = t793 * t660 - t790 * t661;
t653 = m(2) * t769 - t796 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t806;
t768 = t791 * g(1) - t794 * g(2);
t745 = -qJDD(1) * pkin(1) - t796 * pkin(6) - t768;
t804 = -t759 * pkin(2) + t745 + (-t758 - t809) * qJ(3);
t693 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t815 + t804;
t687 = -pkin(2) * t810 + t759 * pkin(3) - pkin(7) * t822 - t804 + (t767 + t829) * t815;
t682 = -t709 * pkin(4) - t736 * qJ(5) + t744 * t728 + qJDD(5) + t687;
t675 = m(6) * t682 - t709 * mrSges(6,1) + t710 * mrSges(6,2) - t743 * t726 + t744 * t729;
t800 = -m(5) * t687 + t709 * mrSges(5,1) - t710 * mrSges(5,2) + t743 * t727 - t744 * t730 - t675;
t669 = m(4) * t693 - t759 * mrSges(4,1) - t758 * mrSges(4,3) - t764 * t815 - t766 * t814 + t800;
t798 = -m(3) * t745 + t759 * mrSges(3,1) - t758 * mrSges(3,2) - t763 * t815 + t765 * t814 - t669;
t667 = m(2) * t768 + qJDD(1) * mrSges(2,1) - t796 * mrSges(2,2) + t798;
t821 = t791 * t653 + t794 * t667;
t655 = t790 * t660 + t793 * t661;
t820 = -t836 * t743 - t838 * t744 - t834 * t781;
t818 = t835 * qJD(2) + (t825 * t790 + t824 * t793) * qJD(1);
t807 = t794 * t653 - t791 * t667;
t656 = -mrSges(5,1) * t687 + mrSges(5,3) * t685 - mrSges(6,1) * t682 + mrSges(6,3) * t680 - pkin(4) * t675 + qJ(5) * t811 + (-qJ(5) * t729 + t819) * t781 + (-qJ(5) * mrSges(6,2) + t836) * t780 + t820 * t744 + t826 * t710 + t842 * t709;
t662 = mrSges(5,2) * t687 + mrSges(6,2) * t682 - mrSges(5,3) * t684 - mrSges(6,3) * t678 - qJ(5) * t674 + t826 * t709 + t839 * t710 - t820 * t743 + t838 * t780 - t833 * t781;
t648 = -mrSges(3,1) * t745 - mrSges(4,1) * t693 + mrSges(4,2) * t698 + mrSges(3,3) * t725 - pkin(2) * t669 - pkin(3) * t800 - pkin(7) * t805 + t816 * qJD(2) + t824 * qJDD(2) - t792 * t656 - t789 * t662 + t827 * t758 + t837 * t759 - t818 * t815;
t650 = mrSges(3,2) * t745 + mrSges(4,2) * t708 - mrSges(3,3) * t724 - mrSges(4,3) * t693 - pkin(7) * t665 - qJ(3) * t669 + t817 * qJD(2) + t825 * qJDD(2) - t789 * t656 + t792 * t662 + t840 * t758 + t827 * t759 + t818 * t814;
t802 = mrSges(2,1) * t768 - mrSges(2,2) * t769 + Ifges(2,3) * qJDD(1) + pkin(1) * t798 + pkin(6) * t806 + t793 * t648 + t790 * t650;
t646 = mrSges(2,1) * g(3) + mrSges(2,3) * t769 + t796 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t655 - t832;
t645 = -mrSges(2,2) * g(3) - mrSges(2,3) * t768 + Ifges(2,5) * qJDD(1) - t796 * Ifges(2,6) - pkin(6) * t655 - t790 * t648 + t793 * t650;
t1 = [-m(1) * g(1) + t807; -m(1) * g(2) + t821; (-m(1) - m(2)) * g(3) + t655; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t821 + t794 * t645 - t791 * t646; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t807 + t791 * t645 + t794 * t646; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t802; t802; t832; t664; t841; t675;];
tauJB = t1;
