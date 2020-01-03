% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR11_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:49
% EndTime: 2019-12-31 21:32:55
% DurationCPUTime: 3.91s
% Computational Cost: add. (41510->292), mult. (81414->351), div. (0->0), fcn. (51320->8), ass. (0->120)
t800 = Ifges(4,1) + Ifges(5,1);
t791 = Ifges(4,4) - Ifges(5,5);
t790 = Ifges(4,5) + Ifges(5,4);
t799 = -Ifges(4,2) - Ifges(5,3);
t789 = Ifges(4,6) - Ifges(5,6);
t798 = Ifges(4,3) + Ifges(5,2);
t755 = sin(qJ(3));
t756 = sin(qJ(2));
t781 = qJD(1) * t756;
t793 = cos(qJ(3));
t733 = -t793 * qJD(2) + t755 * t781;
t759 = cos(qJ(2));
t779 = qJD(1) * qJD(2);
t777 = t759 * t779;
t737 = t756 * qJDD(1) + t777;
t701 = -t733 * qJD(3) + t755 * qJDD(2) + t793 * t737;
t757 = sin(qJ(1));
t760 = cos(qJ(1));
t744 = -t760 * g(1) - t757 * g(2);
t762 = qJD(1) ^ 2;
t724 = -t762 * pkin(1) + qJDD(1) * pkin(6) + t744;
t713 = -t759 * g(3) - t756 * t724;
t736 = (-pkin(2) * t759 - pkin(7) * t756) * qJD(1);
t761 = qJD(2) ^ 2;
t769 = qJDD(2) * pkin(2) + t761 * pkin(7) - t736 * t781 + t713;
t780 = t759 * qJD(1);
t747 = qJD(3) - t780;
t787 = t733 * t747;
t797 = (-t701 + t787) * qJ(4) - t769;
t743 = t757 * g(1) - t760 * g(2);
t723 = -qJDD(1) * pkin(1) - t762 * pkin(6) - t743;
t778 = t756 * t779;
t738 = t759 * qJDD(1) - t778;
t681 = (-t737 - t777) * pkin(7) + (-t738 + t778) * pkin(2) + t723;
t714 = -t756 * g(3) + t759 * t724;
t685 = -t761 * pkin(2) + qJDD(2) * pkin(7) + t736 * t780 + t714;
t672 = t793 * t681 - t755 * t685;
t734 = t755 * qJD(2) + t793 * t781;
t706 = t733 * pkin(3) - t734 * qJ(4);
t732 = qJDD(3) - t738;
t746 = t747 ^ 2;
t664 = -t732 * pkin(3) - t746 * qJ(4) + t734 * t706 + qJDD(4) - t672;
t657 = (-t701 - t787) * pkin(8) + (t733 * t734 - t732) * pkin(4) + t664;
t673 = t755 * t681 + t793 * t685;
t794 = 2 * qJD(4);
t662 = -t746 * pkin(3) + t732 * qJ(4) - t733 * t706 + t747 * t794 + t673;
t700 = t734 * qJD(3) - t793 * qJDD(2) + t755 * t737;
t715 = -t747 * pkin(4) - t734 * pkin(8);
t731 = t733 ^ 2;
t659 = -t731 * pkin(4) + t700 * pkin(8) + t747 * t715 + t662;
t754 = sin(qJ(5));
t758 = cos(qJ(5));
t656 = t754 * t657 + t758 * t659;
t660 = -t731 * pkin(8) + (-pkin(3) - pkin(4)) * t700 + (-pkin(3) * t747 + t715 + t794) * t734 - t797;
t703 = t754 * t733 + t758 * t734;
t670 = -t703 * qJD(5) + t758 * t700 - t754 * t701;
t702 = t758 * t733 - t754 * t734;
t671 = t702 * qJD(5) + t754 * t700 + t758 * t701;
t745 = qJD(5) - t747;
t674 = Ifges(6,5) * t703 + Ifges(6,6) * t702 + Ifges(6,3) * t745;
t676 = Ifges(6,1) * t703 + Ifges(6,4) * t702 + Ifges(6,5) * t745;
t728 = qJDD(5) - t732;
t645 = -mrSges(6,1) * t660 + mrSges(6,3) * t656 + Ifges(6,4) * t671 + Ifges(6,2) * t670 + Ifges(6,6) * t728 - t703 * t674 + t745 * t676;
t655 = t758 * t657 - t754 * t659;
t675 = Ifges(6,4) * t703 + Ifges(6,2) * t702 + Ifges(6,6) * t745;
t646 = mrSges(6,2) * t660 - mrSges(6,3) * t655 + Ifges(6,1) * t671 + Ifges(6,4) * t670 + Ifges(6,5) * t728 + t702 * t674 - t745 * t675;
t663 = -0.2e1 * qJD(4) * t734 + (t734 * t747 + t700) * pkin(3) + t797;
t711 = -t747 * mrSges(5,1) + t734 * mrSges(5,2);
t712 = -t733 * mrSges(5,2) + t747 * mrSges(5,3);
t686 = -t745 * mrSges(6,2) + t702 * mrSges(6,3);
t687 = t745 * mrSges(6,1) - t703 * mrSges(6,3);
t772 = -m(6) * t660 + t670 * mrSges(6,1) - t671 * mrSges(6,2) + t702 * t686 - t703 * t687;
t650 = m(5) * t663 + t700 * mrSges(5,1) - t701 * mrSges(5,3) - t734 * t711 + t733 * t712 + t772;
t679 = -t702 * mrSges(6,1) + t703 * mrSges(6,2);
t652 = m(6) * t655 + t728 * mrSges(6,1) - t671 * mrSges(6,3) - t703 * t679 + t745 * t686;
t653 = m(6) * t656 - t728 * mrSges(6,2) + t670 * mrSges(6,3) + t702 * t679 - t745 * t687;
t774 = -t754 * t652 + t758 * t653;
t783 = -t791 * t733 + t800 * t734 + t790 * t747;
t785 = t789 * t733 - t790 * t734 - t798 * t747;
t624 = mrSges(4,1) * t769 - mrSges(5,1) * t663 + mrSges(5,2) * t662 + mrSges(4,3) * t673 - pkin(3) * t650 - pkin(4) * t772 - pkin(8) * t774 - t758 * t645 - t754 * t646 + t799 * t700 + t791 * t701 + t789 * t732 + t785 * t734 + t783 * t747;
t644 = t758 * t652 + t754 * t653;
t784 = t799 * t733 + t791 * t734 + t789 * t747;
t625 = -mrSges(4,2) * t769 + mrSges(5,2) * t664 - mrSges(4,3) * t672 - mrSges(5,3) * t663 - pkin(8) * t644 - qJ(4) * t650 - t754 * t645 + t758 * t646 - t791 * t700 + t800 * t701 + t790 * t732 + t785 * t733 - t784 * t747;
t710 = t747 * mrSges(4,1) - t734 * mrSges(4,3);
t770 = m(5) * t662 + t732 * mrSges(5,3) + t747 * t711 + t774;
t707 = t733 * mrSges(5,1) - t734 * mrSges(5,3);
t782 = -t733 * mrSges(4,1) - t734 * mrSges(4,2) - t707;
t792 = -mrSges(4,3) - mrSges(5,2);
t640 = m(4) * t673 - t732 * mrSges(4,2) + t792 * t700 - t747 * t710 + t782 * t733 + t770;
t709 = -t747 * mrSges(4,2) - t733 * mrSges(4,3);
t767 = -m(5) * t664 + t732 * mrSges(5,1) + t747 * t712 - t644;
t641 = m(4) * t672 + t732 * mrSges(4,1) + t792 * t701 + t747 * t709 + t782 * t734 + t767;
t638 = t793 * t640 - t755 * t641;
t649 = m(4) * t769 - t700 * mrSges(4,1) - t701 * mrSges(4,2) - t733 * t709 - t734 * t710 - t650;
t721 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t756 + Ifges(3,2) * t759) * qJD(1);
t722 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t756 + Ifges(3,4) * t759) * qJD(1);
t796 = mrSges(3,1) * t713 - mrSges(3,2) * t714 + Ifges(3,5) * t737 + Ifges(3,6) * t738 + Ifges(3,3) * qJDD(2) + pkin(2) * t649 + pkin(7) * t638 + (t721 * t756 - t722 * t759) * qJD(1) + t793 * t624 + t755 * t625;
t643 = t701 * mrSges(5,2) + t734 * t707 - t767;
t766 = -mrSges(6,1) * t655 + mrSges(6,2) * t656 - Ifges(6,5) * t671 - Ifges(6,6) * t670 - Ifges(6,3) * t728 - t703 * t675 + t702 * t676;
t795 = -t789 * t700 + t790 * t701 + t798 * t732 + t783 * t733 + t784 * t734 + mrSges(4,1) * t672 - mrSges(5,1) * t664 - mrSges(4,2) * t673 + mrSges(5,3) * t662 - pkin(3) * t643 - pkin(4) * t644 + qJ(4) * (-t700 * mrSges(5,2) - t733 * t707 + t770) + t766;
t735 = (-mrSges(3,1) * t759 + mrSges(3,2) * t756) * qJD(1);
t740 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t781;
t636 = m(3) * t714 - qJDD(2) * mrSges(3,2) + t738 * mrSges(3,3) - qJD(2) * t740 + t735 * t780 + t638;
t741 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t780;
t648 = m(3) * t713 + qJDD(2) * mrSges(3,1) - t737 * mrSges(3,3) + qJD(2) * t741 - t735 * t781 + t649;
t775 = t759 * t636 - t756 * t648;
t628 = m(2) * t744 - t762 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t775;
t637 = t755 * t640 + t793 * t641;
t765 = -m(3) * t723 + t738 * mrSges(3,1) - t737 * mrSges(3,2) - t740 * t781 + t741 * t780 - t637;
t632 = m(2) * t743 + qJDD(1) * mrSges(2,1) - t762 * mrSges(2,2) + t765;
t786 = t757 * t628 + t760 * t632;
t630 = t756 * t636 + t759 * t648;
t776 = t760 * t628 - t757 * t632;
t720 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t756 + Ifges(3,6) * t759) * qJD(1);
t621 = mrSges(3,2) * t723 - mrSges(3,3) * t713 + Ifges(3,1) * t737 + Ifges(3,4) * t738 + Ifges(3,5) * qJDD(2) - pkin(7) * t637 - qJD(2) * t721 - t755 * t624 + t793 * t625 + t720 * t780;
t623 = -mrSges(3,1) * t723 + mrSges(3,3) * t714 + Ifges(3,4) * t737 + Ifges(3,2) * t738 + Ifges(3,6) * qJDD(2) - pkin(2) * t637 + qJD(2) * t722 - t720 * t781 - t795;
t768 = mrSges(2,1) * t743 - mrSges(2,2) * t744 + Ifges(2,3) * qJDD(1) + pkin(1) * t765 + pkin(6) * t775 + t756 * t621 + t759 * t623;
t619 = mrSges(2,1) * g(3) + mrSges(2,3) * t744 + t762 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t630 - t796;
t618 = -mrSges(2,2) * g(3) - mrSges(2,3) * t743 + Ifges(2,5) * qJDD(1) - t762 * Ifges(2,6) - pkin(6) * t630 + t759 * t621 - t756 * t623;
t1 = [-m(1) * g(1) + t776; -m(1) * g(2) + t786; (-m(1) - m(2)) * g(3) + t630; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t786 + t760 * t618 - t757 * t619; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t776 + t757 * t618 + t760 * t619; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t768; t768; t796; t795; t643; -t766;];
tauJB = t1;
