% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:03
% EndTime: 2019-12-31 22:25:11
% DurationCPUTime: 8.31s
% Computational Cost: add. (131751->316), mult. (264899->398), div. (0->0), fcn. (186114->10), ass. (0->128)
t761 = sin(qJ(2));
t766 = cos(qJ(2));
t785 = qJD(1) * qJD(2);
t742 = t761 * qJDD(1) + t766 * t785;
t762 = sin(qJ(1));
t767 = cos(qJ(1));
t749 = -t767 * g(1) - t762 * g(2);
t768 = qJD(1) ^ 2;
t736 = -t768 * pkin(1) + qJDD(1) * pkin(6) + t749;
t789 = t761 * t736;
t790 = pkin(2) * t768;
t698 = qJDD(2) * pkin(2) - t742 * pkin(7) - t789 + (pkin(7) * t785 + t761 * t790 - g(3)) * t766;
t723 = -t761 * g(3) + t766 * t736;
t743 = t766 * qJDD(1) - t761 * t785;
t787 = qJD(1) * t761;
t747 = qJD(2) * pkin(2) - pkin(7) * t787;
t757 = t766 ^ 2;
t699 = t743 * pkin(7) - qJD(2) * t747 - t757 * t790 + t723;
t760 = sin(qJ(3));
t765 = cos(qJ(3));
t681 = t760 * t698 + t765 * t699;
t734 = (t760 * t766 + t761 * t765) * qJD(1);
t707 = -t734 * qJD(3) - t760 * t742 + t765 * t743;
t786 = qJD(1) * t766;
t733 = -t760 * t787 + t765 * t786;
t717 = -t733 * mrSges(4,1) + t734 * mrSges(4,2);
t755 = qJD(2) + qJD(3);
t725 = t755 * mrSges(4,1) - t734 * mrSges(4,3);
t754 = qJDD(2) + qJDD(3);
t708 = t733 * qJD(3) + t765 * t742 + t760 * t743;
t748 = t762 * g(1) - t767 * g(2);
t778 = -qJDD(1) * pkin(1) - t748;
t709 = -t743 * pkin(2) + t747 * t787 + (-pkin(7) * t757 - pkin(6)) * t768 + t778;
t670 = (-t733 * t755 - t708) * pkin(8) + (t734 * t755 - t707) * pkin(3) + t709;
t718 = -t733 * pkin(3) - t734 * pkin(8);
t753 = t755 ^ 2;
t673 = -t753 * pkin(3) + t754 * pkin(8) + t733 * t718 + t681;
t759 = sin(qJ(4));
t764 = cos(qJ(4));
t659 = t764 * t670 - t759 * t673;
t720 = -t759 * t734 + t764 * t755;
t684 = t720 * qJD(4) + t764 * t708 + t759 * t754;
t706 = qJDD(4) - t707;
t721 = t764 * t734 + t759 * t755;
t729 = qJD(4) - t733;
t657 = (t720 * t729 - t684) * pkin(9) + (t720 * t721 + t706) * pkin(4) + t659;
t660 = t759 * t670 + t764 * t673;
t683 = -t721 * qJD(4) - t759 * t708 + t764 * t754;
t712 = t729 * pkin(4) - t721 * pkin(9);
t719 = t720 ^ 2;
t658 = -t719 * pkin(4) + t683 * pkin(9) - t729 * t712 + t660;
t758 = sin(qJ(5));
t763 = cos(qJ(5));
t655 = t763 * t657 - t758 * t658;
t691 = t763 * t720 - t758 * t721;
t666 = t691 * qJD(5) + t758 * t683 + t763 * t684;
t692 = t758 * t720 + t763 * t721;
t678 = -t691 * mrSges(6,1) + t692 * mrSges(6,2);
t727 = qJD(5) + t729;
t685 = -t727 * mrSges(6,2) + t691 * mrSges(6,3);
t701 = qJDD(5) + t706;
t650 = m(6) * t655 + t701 * mrSges(6,1) - t666 * mrSges(6,3) - t692 * t678 + t727 * t685;
t656 = t758 * t657 + t763 * t658;
t665 = -t692 * qJD(5) + t763 * t683 - t758 * t684;
t686 = t727 * mrSges(6,1) - t692 * mrSges(6,3);
t651 = m(6) * t656 - t701 * mrSges(6,2) + t665 * mrSges(6,3) + t691 * t678 - t727 * t686;
t642 = t763 * t650 + t758 * t651;
t696 = -t720 * mrSges(5,1) + t721 * mrSges(5,2);
t710 = -t729 * mrSges(5,2) + t720 * mrSges(5,3);
t640 = m(5) * t659 + t706 * mrSges(5,1) - t684 * mrSges(5,3) - t721 * t696 + t729 * t710 + t642;
t711 = t729 * mrSges(5,1) - t721 * mrSges(5,3);
t780 = -t758 * t650 + t763 * t651;
t641 = m(5) * t660 - t706 * mrSges(5,2) + t683 * mrSges(5,3) + t720 * t696 - t729 * t711 + t780;
t781 = -t759 * t640 + t764 * t641;
t633 = m(4) * t681 - t754 * mrSges(4,2) + t707 * mrSges(4,3) + t733 * t717 - t755 * t725 + t781;
t680 = t765 * t698 - t760 * t699;
t724 = -t755 * mrSges(4,2) + t733 * mrSges(4,3);
t672 = -t754 * pkin(3) - t753 * pkin(8) + t734 * t718 - t680;
t661 = -t683 * pkin(4) - t719 * pkin(9) + t721 * t712 + t672;
t776 = m(6) * t661 - t665 * mrSges(6,1) + t666 * mrSges(6,2) - t691 * t685 + t692 * t686;
t772 = -m(5) * t672 + t683 * mrSges(5,1) - t684 * mrSges(5,2) + t720 * t710 - t721 * t711 - t776;
t646 = m(4) * t680 + t754 * mrSges(4,1) - t708 * mrSges(4,3) - t734 * t717 + t755 * t724 + t772;
t625 = t760 * t633 + t765 * t646;
t722 = -t766 * g(3) - t789;
t731 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t761 + Ifges(3,2) * t766) * qJD(1);
t732 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t761 + Ifges(3,4) * t766) * qJD(1);
t674 = Ifges(6,5) * t692 + Ifges(6,6) * t691 + Ifges(6,3) * t727;
t676 = Ifges(6,1) * t692 + Ifges(6,4) * t691 + Ifges(6,5) * t727;
t643 = -mrSges(6,1) * t661 + mrSges(6,3) * t656 + Ifges(6,4) * t666 + Ifges(6,2) * t665 + Ifges(6,6) * t701 - t692 * t674 + t727 * t676;
t675 = Ifges(6,4) * t692 + Ifges(6,2) * t691 + Ifges(6,6) * t727;
t644 = mrSges(6,2) * t661 - mrSges(6,3) * t655 + Ifges(6,1) * t666 + Ifges(6,4) * t665 + Ifges(6,5) * t701 + t691 * t674 - t727 * t675;
t687 = Ifges(5,5) * t721 + Ifges(5,6) * t720 + Ifges(5,3) * t729;
t689 = Ifges(5,1) * t721 + Ifges(5,4) * t720 + Ifges(5,5) * t729;
t624 = -mrSges(5,1) * t672 + mrSges(5,3) * t660 + Ifges(5,4) * t684 + Ifges(5,2) * t683 + Ifges(5,6) * t706 - pkin(4) * t776 + pkin(9) * t780 + t763 * t643 + t758 * t644 - t721 * t687 + t729 * t689;
t688 = Ifges(5,4) * t721 + Ifges(5,2) * t720 + Ifges(5,6) * t729;
t627 = mrSges(5,2) * t672 - mrSges(5,3) * t659 + Ifges(5,1) * t684 + Ifges(5,4) * t683 + Ifges(5,5) * t706 - pkin(9) * t642 - t758 * t643 + t763 * t644 + t720 * t687 - t729 * t688;
t714 = Ifges(4,4) * t734 + Ifges(4,2) * t733 + Ifges(4,6) * t755;
t715 = Ifges(4,1) * t734 + Ifges(4,4) * t733 + Ifges(4,5) * t755;
t773 = -mrSges(4,1) * t680 + mrSges(4,2) * t681 - Ifges(4,5) * t708 - Ifges(4,6) * t707 - Ifges(4,3) * t754 - pkin(3) * t772 - pkin(8) * t781 - t764 * t624 - t759 * t627 - t734 * t714 + t733 * t715;
t791 = mrSges(3,1) * t722 - mrSges(3,2) * t723 + Ifges(3,5) * t742 + Ifges(3,6) * t743 + Ifges(3,3) * qJDD(2) + pkin(2) * t625 + (t761 * t731 - t766 * t732) * qJD(1) - t773;
t741 = (-mrSges(3,1) * t766 + mrSges(3,2) * t761) * qJD(1);
t746 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t786;
t621 = m(3) * t722 + qJDD(2) * mrSges(3,1) - t742 * mrSges(3,3) + qJD(2) * t746 - t741 * t787 + t625;
t745 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t787;
t782 = t765 * t633 - t760 * t646;
t622 = m(3) * t723 - qJDD(2) * mrSges(3,2) + t743 * mrSges(3,3) - qJD(2) * t745 + t741 * t786 + t782;
t783 = -t761 * t621 + t766 * t622;
t614 = m(2) * t749 - t768 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t783;
t735 = -t768 * pkin(6) + t778;
t635 = t764 * t640 + t759 * t641;
t774 = m(4) * t709 - t707 * mrSges(4,1) + t708 * mrSges(4,2) - t733 * t724 + t734 * t725 + t635;
t771 = -m(3) * t735 + t743 * mrSges(3,1) - t742 * mrSges(3,2) - t745 * t787 + t746 * t786 - t774;
t629 = m(2) * t748 + qJDD(1) * mrSges(2,1) - t768 * mrSges(2,2) + t771;
t788 = t762 * t614 + t767 * t629;
t616 = t766 * t621 + t761 * t622;
t784 = t767 * t614 - t762 * t629;
t713 = Ifges(4,5) * t734 + Ifges(4,6) * t733 + Ifges(4,3) * t755;
t611 = mrSges(4,2) * t709 - mrSges(4,3) * t680 + Ifges(4,1) * t708 + Ifges(4,4) * t707 + Ifges(4,5) * t754 - pkin(8) * t635 - t759 * t624 + t764 * t627 + t733 * t713 - t755 * t714;
t775 = -mrSges(6,1) * t655 + mrSges(6,2) * t656 - Ifges(6,5) * t666 - Ifges(6,6) * t665 - Ifges(6,3) * t701 - t692 * t675 + t691 * t676;
t769 = mrSges(5,1) * t659 - mrSges(5,2) * t660 + Ifges(5,5) * t684 + Ifges(5,6) * t683 + Ifges(5,3) * t706 + pkin(4) * t642 + t721 * t688 - t720 * t689 - t775;
t617 = -mrSges(4,1) * t709 + mrSges(4,3) * t681 + Ifges(4,4) * t708 + Ifges(4,2) * t707 + Ifges(4,6) * t754 - pkin(3) * t635 - t734 * t713 + t755 * t715 - t769;
t730 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t761 + Ifges(3,6) * t766) * qJD(1);
t607 = -mrSges(3,1) * t735 + mrSges(3,3) * t723 + Ifges(3,4) * t742 + Ifges(3,2) * t743 + Ifges(3,6) * qJDD(2) - pkin(2) * t774 + pkin(7) * t782 + qJD(2) * t732 + t760 * t611 + t765 * t617 - t730 * t787;
t610 = mrSges(3,2) * t735 - mrSges(3,3) * t722 + Ifges(3,1) * t742 + Ifges(3,4) * t743 + Ifges(3,5) * qJDD(2) - pkin(7) * t625 - qJD(2) * t731 + t765 * t611 - t760 * t617 + t730 * t786;
t777 = mrSges(2,1) * t748 - mrSges(2,2) * t749 + Ifges(2,3) * qJDD(1) + pkin(1) * t771 + pkin(6) * t783 + t766 * t607 + t761 * t610;
t608 = mrSges(2,1) * g(3) + mrSges(2,3) * t749 + t768 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t616 - t791;
t605 = -mrSges(2,2) * g(3) - mrSges(2,3) * t748 + Ifges(2,5) * qJDD(1) - t768 * Ifges(2,6) - pkin(6) * t616 - t761 * t607 + t766 * t610;
t1 = [-m(1) * g(1) + t784; -m(1) * g(2) + t788; (-m(1) - m(2)) * g(3) + t616; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t788 + t767 * t605 - t762 * t608; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t784 + t762 * t605 + t767 * t608; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t777; t777; t791; -t773; t769; -t775;];
tauJB = t1;
