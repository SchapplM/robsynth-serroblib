% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:30
% EndTime: 2019-12-31 17:53:32
% DurationCPUTime: 1.88s
% Computational Cost: add. (11487->227), mult. (27090->271), div. (0->0), fcn. (15956->6), ass. (0->100)
t785 = Ifges(5,1) + Ifges(6,1);
t769 = Ifges(5,4) - Ifges(6,5);
t780 = Ifges(5,5) + Ifges(6,4);
t784 = Ifges(5,2) + Ifges(6,3);
t778 = Ifges(5,6) - Ifges(6,6);
t732 = sin(qJ(1));
t734 = cos(qJ(1));
t702 = -t734 * g(1) - t732 * g(2);
t736 = qJD(1) ^ 2;
t783 = -t736 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t702;
t729 = sin(pkin(7));
t730 = cos(pkin(7));
t782 = (Ifges(3,6) - Ifges(4,6)) * t730 + (Ifges(4,4) + Ifges(3,5)) * t729;
t777 = Ifges(5,3) + Ifges(6,2);
t731 = sin(qJ(4));
t733 = cos(qJ(4));
t744 = t729 * t731 + t730 * t733;
t692 = t744 * qJD(1);
t745 = t729 * t733 - t730 * t731;
t693 = t745 * qJD(1);
t776 = -t778 * qJD(4) + t784 * t692 - t769 * t693;
t775 = t780 * qJD(4) - t769 * t692 + t785 * t693;
t720 = t729 ^ 2;
t721 = t730 ^ 2;
t765 = t721 * t736;
t774 = t720 * t736 + t765;
t770 = Ifges(3,4) - Ifges(4,5);
t773 = t770 * t729;
t701 = t732 * g(1) - t734 * g(2);
t691 = -qJDD(1) * pkin(1) - t736 * qJ(2) + qJDD(2) - t701;
t753 = qJDD(1) * t730;
t754 = qJDD(1) * t729;
t757 = t729 * qJD(1);
t675 = -pkin(2) * t753 - qJ(3) * t754 - 0.2e1 * qJD(3) * t757 + t691;
t678 = -t730 * g(3) - t783 * t729;
t772 = -mrSges(5,3) - mrSges(6,2);
t771 = Ifges(3,1) + Ifges(4,1);
t768 = Ifges(3,2) + Ifges(4,3);
t767 = mrSges(3,2) * t729;
t696 = (-mrSges(4,1) * t730 - mrSges(4,3) * t729) * qJD(1);
t697 = (-mrSges(3,1) * t730 + t767) * qJD(1);
t695 = (-pkin(2) * t730 - qJ(3) * t729) * qJD(1);
t651 = t695 * t757 + qJDD(3) - t678;
t644 = (-pkin(3) * t730 * t736 - pkin(6) * qJDD(1)) * t729 + t651;
t679 = -t729 * g(3) + t783 * t730;
t756 = t730 * qJD(1);
t653 = t695 * t756 + t679;
t646 = -pkin(3) * t765 - pkin(6) * t753 + t653;
t642 = t731 * t644 + t733 * t646;
t758 = t693 * qJD(4);
t676 = t744 * qJDD(1) + t758;
t683 = qJD(4) * mrSges(5,1) - t693 * mrSges(5,3);
t664 = t692 * pkin(4) - t693 * qJ(5);
t735 = qJD(4) ^ 2;
t636 = -t735 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t692 * t664 + t642;
t684 = -qJD(4) * mrSges(6,1) + t693 * mrSges(6,2);
t752 = m(6) * t636 + qJDD(4) * mrSges(6,3) + qJD(4) * t684;
t665 = t692 * mrSges(6,1) - t693 * mrSges(6,3);
t762 = -t692 * mrSges(5,1) - t693 * mrSges(5,2) - t665;
t628 = m(5) * t642 - qJDD(4) * mrSges(5,2) - qJD(4) * t683 + t772 * t676 + t762 * t692 + t752;
t641 = t733 * t644 - t731 * t646;
t759 = t692 * qJD(4);
t677 = t745 * qJDD(1) - t759;
t682 = -qJD(4) * mrSges(5,2) - t692 * mrSges(5,3);
t637 = -qJDD(4) * pkin(4) - t735 * qJ(5) + t693 * t664 + qJDD(5) - t641;
t685 = -t692 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t746 = -m(6) * t637 + qJDD(4) * mrSges(6,1) + qJD(4) * t685;
t629 = m(5) * t641 + qJDD(4) * mrSges(5,1) + qJD(4) * t682 + t772 * t677 + t762 * t693 + t746;
t621 = t731 * t628 + t733 * t629;
t741 = m(4) * t651 + t621;
t616 = m(3) * t678 + ((-mrSges(4,2) - mrSges(3,3)) * qJDD(1) + (-t696 - t697) * qJD(1)) * t729 - t741;
t747 = t733 * t628 - t731 * t629;
t743 = m(4) * t653 + mrSges(4,2) * t753 + t696 * t756 + t747;
t617 = m(3) * t679 + (qJDD(1) * mrSges(3,3) + qJD(1) * t697) * t730 + t743;
t748 = -t729 * t616 + t730 * t617;
t610 = m(2) * t702 - t736 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t748;
t650 = pkin(3) * t753 + (-t720 - t721) * t736 * pkin(6) - t675;
t639 = -0.2e1 * qJD(5) * t693 + (-t677 + t759) * qJ(5) + (t676 + t758) * pkin(4) + t650;
t630 = m(6) * t639 + t676 * mrSges(6,1) - t677 * mrSges(6,3) - t693 * t684 + t692 * t685;
t739 = -m(5) * t650 - t676 * mrSges(5,1) - t677 * mrSges(5,2) - t692 * t682 - t693 * t683 - t630;
t626 = m(4) * t675 - mrSges(4,1) * t753 - t774 * mrSges(4,2) - mrSges(4,3) * t754 + t739;
t738 = -m(3) * t691 + mrSges(3,1) * t753 + t774 * mrSges(3,3) - t626;
t623 = t738 + (mrSges(2,1) - t767) * qJDD(1) - t736 * mrSges(2,2) + m(2) * t701;
t764 = t732 * t610 + t734 * t623;
t612 = t730 * t616 + t729 * t617;
t763 = -t777 * qJD(4) + t778 * t692 - t780 * t693;
t760 = t782 * qJD(1);
t749 = t734 * t610 - t732 * t623;
t618 = -mrSges(5,1) * t650 - mrSges(6,1) * t639 + mrSges(6,2) * t636 + mrSges(5,3) * t642 - pkin(4) * t630 + t775 * qJD(4) + t778 * qJDD(4) - t784 * t676 + t769 * t677 + t763 * t693;
t619 = mrSges(5,2) * t650 + mrSges(6,2) * t637 - mrSges(5,3) * t641 - mrSges(6,3) * t639 - qJ(5) * t630 + t776 * qJD(4) + t780 * qJDD(4) - t769 * t676 + t785 * t677 + t763 * t692;
t605 = -mrSges(3,1) * t691 + mrSges(3,3) * t679 - mrSges(4,1) * t675 + mrSges(4,2) * t653 - t731 * t619 - t733 * t618 - pkin(3) * t739 - pkin(6) * t747 - pkin(2) * t626 - t760 * t757 + (t768 * t730 + t773) * qJDD(1);
t607 = mrSges(3,2) * t691 + mrSges(4,2) * t651 - mrSges(3,3) * t678 - mrSges(4,3) * t675 - pkin(6) * t621 - qJ(3) * t626 - t731 * t618 + t733 * t619 + t760 * t756 + (t771 * t729 + t770 * t730) * qJDD(1);
t625 = mrSges(3,2) * t754 - t738;
t740 = mrSges(2,1) * t701 - mrSges(2,2) * t702 + Ifges(2,3) * qJDD(1) - pkin(1) * t625 + qJ(2) * t748 + t730 * t605 + t729 * t607;
t633 = t677 * mrSges(6,2) + t693 * t665 - t746;
t737 = mrSges(5,1) * t641 - mrSges(6,1) * t637 - mrSges(5,2) * t642 + mrSges(6,3) * t636 - pkin(4) * t633 + qJ(5) * t752 - t776 * t693 + (-qJ(5) * t665 + t775) * t692 + t780 * t677 + (-qJ(5) * mrSges(6,2) - t778) * t676 + t777 * qJDD(4);
t620 = (qJDD(1) * mrSges(4,2) + qJD(1) * t696) * t729 + t741;
t603 = t737 + mrSges(2,1) * g(3) + (Ifges(2,6) - t782) * qJDD(1) - qJ(3) * t743 + mrSges(2,3) * t702 - mrSges(3,1) * t678 + mrSges(3,2) * t679 + mrSges(4,1) * t651 - mrSges(4,3) * t653 + pkin(3) * t621 + pkin(2) * t620 - pkin(1) * t612 + (t770 * t721 + (-t773 + (-t768 + t771) * t730) * t729 + Ifges(2,5)) * t736;
t602 = -mrSges(2,2) * g(3) - mrSges(2,3) * t701 + Ifges(2,5) * qJDD(1) - t736 * Ifges(2,6) - qJ(2) * t612 - t729 * t605 + t730 * t607;
t1 = [-m(1) * g(1) + t749; -m(1) * g(2) + t764; (-m(1) - m(2)) * g(3) + t612; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t764 + t734 * t602 - t732 * t603; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t749 + t732 * t602 + t734 * t603; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t740; t740; t625; t620; t737; t633;];
tauJB = t1;
