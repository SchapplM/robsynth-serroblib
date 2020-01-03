% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:40
% EndTime: 2019-12-31 20:20:48
% DurationCPUTime: 7.36s
% Computational Cost: add. (105900->314), mult. (240647->399), div. (0->0), fcn. (166322->10), ass. (0->125)
t746 = sin(qJ(2));
t750 = cos(qJ(2));
t767 = qJD(1) * qJD(2);
t730 = t746 * qJDD(1) + t750 * t767;
t747 = sin(qJ(1));
t751 = cos(qJ(1));
t737 = -t751 * g(1) - t747 * g(2);
t753 = qJD(1) ^ 2;
t725 = -t753 * pkin(1) + qJDD(1) * pkin(6) + t737;
t772 = t746 * t725;
t774 = pkin(2) * t753;
t686 = qJDD(2) * pkin(2) - t730 * qJ(3) - t772 + (qJ(3) * t767 + t746 * t774 - g(3)) * t750;
t710 = -t746 * g(3) + t750 * t725;
t731 = t750 * qJDD(1) - t746 * t767;
t770 = qJD(1) * t746;
t733 = qJD(2) * pkin(2) - qJ(3) * t770;
t741 = t750 ^ 2;
t687 = t731 * qJ(3) - qJD(2) * t733 - t741 * t774 + t710;
t742 = sin(pkin(9));
t743 = cos(pkin(9));
t719 = (t742 * t750 + t743 * t746) * qJD(1);
t669 = -0.2e1 * qJD(3) * t719 + t743 * t686 - t742 * t687;
t769 = qJD(1) * t750;
t718 = -t742 * t770 + t743 * t769;
t670 = 0.2e1 * qJD(3) * t718 + t742 * t686 + t743 * t687;
t700 = -t718 * pkin(3) - t719 * pkin(7);
t752 = qJD(2) ^ 2;
t661 = -t752 * pkin(3) + qJDD(2) * pkin(7) + t718 * t700 + t670;
t736 = t747 * g(1) - t751 * g(2);
t760 = -qJDD(1) * pkin(1) - t736;
t690 = -t731 * pkin(2) + qJDD(3) + t733 * t770 + (-qJ(3) * t741 - pkin(6)) * t753 + t760;
t704 = -t742 * t730 + t743 * t731;
t705 = t743 * t730 + t742 * t731;
t664 = (-qJD(2) * t718 - t705) * pkin(7) + (qJD(2) * t719 - t704) * pkin(3) + t690;
t745 = sin(qJ(4));
t749 = cos(qJ(4));
t651 = -t745 * t661 + t749 * t664;
t707 = t749 * qJD(2) - t745 * t719;
t680 = t707 * qJD(4) + t745 * qJDD(2) + t749 * t705;
t703 = qJDD(4) - t704;
t708 = t745 * qJD(2) + t749 * t719;
t717 = qJD(4) - t718;
t648 = (t707 * t717 - t680) * pkin(8) + (t707 * t708 + t703) * pkin(4) + t651;
t652 = t749 * t661 + t745 * t664;
t679 = -t708 * qJD(4) + t749 * qJDD(2) - t745 * t705;
t693 = t717 * pkin(4) - t708 * pkin(8);
t706 = t707 ^ 2;
t649 = -t706 * pkin(4) + t679 * pkin(8) - t717 * t693 + t652;
t744 = sin(qJ(5));
t748 = cos(qJ(5));
t647 = t744 * t648 + t748 * t649;
t660 = -qJDD(2) * pkin(3) - t752 * pkin(7) + t719 * t700 - t669;
t650 = -t679 * pkin(4) - t706 * pkin(8) + t708 * t693 + t660;
t685 = t744 * t707 + t748 * t708;
t656 = -t685 * qJD(5) + t748 * t679 - t744 * t680;
t684 = t748 * t707 - t744 * t708;
t657 = t684 * qJD(5) + t744 * t679 + t748 * t680;
t713 = qJD(5) + t717;
t665 = Ifges(6,5) * t685 + Ifges(6,6) * t684 + Ifges(6,3) * t713;
t667 = Ifges(6,1) * t685 + Ifges(6,4) * t684 + Ifges(6,5) * t713;
t701 = qJDD(5) + t703;
t635 = -mrSges(6,1) * t650 + mrSges(6,3) * t647 + Ifges(6,4) * t657 + Ifges(6,2) * t656 + Ifges(6,6) * t701 - t685 * t665 + t713 * t667;
t646 = t748 * t648 - t744 * t649;
t666 = Ifges(6,4) * t685 + Ifges(6,2) * t684 + Ifges(6,6) * t713;
t636 = mrSges(6,2) * t650 - mrSges(6,3) * t646 + Ifges(6,1) * t657 + Ifges(6,4) * t656 + Ifges(6,5) * t701 + t684 * t665 - t713 * t666;
t674 = Ifges(5,5) * t708 + Ifges(5,6) * t707 + Ifges(5,3) * t717;
t676 = Ifges(5,1) * t708 + Ifges(5,4) * t707 + Ifges(5,5) * t717;
t672 = -t713 * mrSges(6,2) + t684 * mrSges(6,3);
t673 = t713 * mrSges(6,1) - t685 * mrSges(6,3);
t758 = m(6) * t650 - t656 * mrSges(6,1) + t657 * mrSges(6,2) - t684 * t672 + t685 * t673;
t671 = -t684 * mrSges(6,1) + t685 * mrSges(6,2);
t642 = m(6) * t646 + t701 * mrSges(6,1) - t657 * mrSges(6,3) - t685 * t671 + t713 * t672;
t643 = m(6) * t647 - t701 * mrSges(6,2) + t656 * mrSges(6,3) + t684 * t671 - t713 * t673;
t763 = -t744 * t642 + t748 * t643;
t617 = -mrSges(5,1) * t660 + mrSges(5,3) * t652 + Ifges(5,4) * t680 + Ifges(5,2) * t679 + Ifges(5,6) * t703 - pkin(4) * t758 + pkin(8) * t763 + t748 * t635 + t744 * t636 - t708 * t674 + t717 * t676;
t634 = t748 * t642 + t744 * t643;
t688 = -t707 * mrSges(5,1) + t708 * mrSges(5,2);
t691 = -t717 * mrSges(5,2) + t707 * mrSges(5,3);
t632 = m(5) * t651 + t703 * mrSges(5,1) - t680 * mrSges(5,3) - t708 * t688 + t717 * t691 + t634;
t692 = t717 * mrSges(5,1) - t708 * mrSges(5,3);
t633 = m(5) * t652 - t703 * mrSges(5,2) + t679 * mrSges(5,3) + t707 * t688 - t717 * t692 + t763;
t628 = -t745 * t632 + t749 * t633;
t698 = -t718 * mrSges(4,1) + t719 * mrSges(4,2);
t712 = qJD(2) * mrSges(4,1) - t719 * mrSges(4,3);
t625 = m(4) * t670 - qJDD(2) * mrSges(4,2) + t704 * mrSges(4,3) - qJD(2) * t712 + t718 * t698 + t628;
t644 = -m(5) * t660 + t679 * mrSges(5,1) - t680 * mrSges(5,2) + t707 * t691 - t708 * t692 - t758;
t711 = -qJD(2) * mrSges(4,2) + t718 * mrSges(4,3);
t638 = m(4) * t669 + qJDD(2) * mrSges(4,1) - t705 * mrSges(4,3) + qJD(2) * t711 - t719 * t698 + t644;
t618 = t742 * t625 + t743 * t638;
t675 = Ifges(5,4) * t708 + Ifges(5,2) * t707 + Ifges(5,6) * t717;
t619 = mrSges(5,2) * t660 - mrSges(5,3) * t651 + Ifges(5,1) * t680 + Ifges(5,4) * t679 + Ifges(5,5) * t703 - pkin(8) * t634 - t744 * t635 + t748 * t636 + t707 * t674 - t717 * t675;
t695 = Ifges(4,4) * t719 + Ifges(4,2) * t718 + Ifges(4,6) * qJD(2);
t696 = Ifges(4,1) * t719 + Ifges(4,4) * t718 + Ifges(4,5) * qJD(2);
t709 = -t750 * g(3) - t772;
t721 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t746 + Ifges(3,2) * t750) * qJD(1);
t722 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t746 + Ifges(3,4) * t750) * qJD(1);
t775 = mrSges(3,1) * t709 + mrSges(4,1) * t669 - mrSges(3,2) * t710 - mrSges(4,2) * t670 + Ifges(3,5) * t730 + Ifges(4,5) * t705 + Ifges(3,6) * t731 + Ifges(4,6) * t704 + pkin(2) * t618 + pkin(3) * t644 + pkin(7) * t628 + t749 * t617 + t745 * t619 + t719 * t695 - t718 * t696 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t746 * t721 - t750 * t722) * qJD(1);
t729 = (-mrSges(3,1) * t750 + mrSges(3,2) * t746) * qJD(1);
t735 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t769;
t615 = m(3) * t709 + qJDD(2) * mrSges(3,1) - t730 * mrSges(3,3) + qJD(2) * t735 - t729 * t770 + t618;
t734 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t770;
t764 = t743 * t625 - t742 * t638;
t616 = m(3) * t710 - qJDD(2) * mrSges(3,2) + t731 * mrSges(3,3) - qJD(2) * t734 + t729 * t769 + t764;
t765 = -t746 * t615 + t750 * t616;
t608 = m(2) * t737 - t753 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t765;
t627 = t749 * t632 + t745 * t633;
t626 = m(4) * t690 - t704 * mrSges(4,1) + t705 * mrSges(4,2) - t718 * t711 + t719 * t712 + t627;
t724 = -t753 * pkin(6) + t760;
t756 = -m(3) * t724 + t731 * mrSges(3,1) - t730 * mrSges(3,2) - t734 * t770 + t735 * t769 - t626;
t621 = m(2) * t736 + qJDD(1) * mrSges(2,1) - t753 * mrSges(2,2) + t756;
t771 = t747 * t608 + t751 * t621;
t610 = t750 * t615 + t746 * t616;
t766 = t751 * t608 - t747 * t621;
t694 = Ifges(4,5) * t719 + Ifges(4,6) * t718 + Ifges(4,3) * qJD(2);
t605 = mrSges(4,2) * t690 - mrSges(4,3) * t669 + Ifges(4,1) * t705 + Ifges(4,4) * t704 + Ifges(4,5) * qJDD(2) - pkin(7) * t627 - qJD(2) * t695 - t745 * t617 + t749 * t619 + t718 * t694;
t757 = -mrSges(6,1) * t646 + mrSges(6,2) * t647 - Ifges(6,5) * t657 - Ifges(6,6) * t656 - Ifges(6,3) * t701 - t685 * t666 + t684 * t667;
t755 = mrSges(5,1) * t651 - mrSges(5,2) * t652 + Ifges(5,5) * t680 + Ifges(5,6) * t679 + Ifges(5,3) * t703 + pkin(4) * t634 + t708 * t675 - t707 * t676 - t757;
t611 = -mrSges(4,1) * t690 + mrSges(4,3) * t670 + Ifges(4,4) * t705 + Ifges(4,2) * t704 + Ifges(4,6) * qJDD(2) - pkin(3) * t627 + qJD(2) * t696 - t719 * t694 - t755;
t720 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t746 + Ifges(3,6) * t750) * qJD(1);
t601 = -mrSges(3,1) * t724 + mrSges(3,3) * t710 + Ifges(3,4) * t730 + Ifges(3,2) * t731 + Ifges(3,6) * qJDD(2) - pkin(2) * t626 + qJ(3) * t764 + qJD(2) * t722 + t742 * t605 + t743 * t611 - t720 * t770;
t604 = mrSges(3,2) * t724 - mrSges(3,3) * t709 + Ifges(3,1) * t730 + Ifges(3,4) * t731 + Ifges(3,5) * qJDD(2) - qJ(3) * t618 - qJD(2) * t721 + t743 * t605 - t742 * t611 + t720 * t769;
t759 = mrSges(2,1) * t736 - mrSges(2,2) * t737 + Ifges(2,3) * qJDD(1) + pkin(1) * t756 + pkin(6) * t765 + t750 * t601 + t746 * t604;
t602 = mrSges(2,1) * g(3) + mrSges(2,3) * t737 + t753 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t610 - t775;
t599 = -mrSges(2,2) * g(3) - mrSges(2,3) * t736 + Ifges(2,5) * qJDD(1) - t753 * Ifges(2,6) - pkin(6) * t610 - t746 * t601 + t750 * t604;
t1 = [-m(1) * g(1) + t766; -m(1) * g(2) + t771; (-m(1) - m(2)) * g(3) + t610; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t771 + t751 * t599 - t747 * t602; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t766 + t747 * t599 + t751 * t602; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t759; t759; t775; t626; t755; -t757;];
tauJB = t1;
