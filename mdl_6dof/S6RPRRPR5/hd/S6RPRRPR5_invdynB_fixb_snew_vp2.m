% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 22:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:39:07
% EndTime: 2019-05-05 22:39:18
% DurationCPUTime: 10.12s
% Computational Cost: add. (135224->345), mult. (329634->413), div. (0->0), fcn. (254132->10), ass. (0->143)
t780 = Ifges(5,1) + Ifges(6,2);
t779 = -Ifges(6,1) - Ifges(5,3);
t774 = Ifges(5,4) + Ifges(6,6);
t773 = Ifges(5,5) - Ifges(6,4);
t778 = Ifges(5,2) + Ifges(6,3);
t772 = Ifges(5,6) - Ifges(6,5);
t739 = qJD(1) ^ 2;
t777 = -2 * qJD(5);
t776 = cos(qJ(4));
t731 = cos(pkin(10));
t775 = pkin(2) * t731;
t730 = sin(pkin(10));
t771 = mrSges(3,2) * t730;
t734 = sin(qJ(3));
t737 = cos(qJ(3));
t749 = -t730 * t734 + t731 * t737;
t711 = t749 * qJD(1);
t750 = t730 * t737 + t731 * t734;
t712 = t750 * qJD(1);
t733 = sin(qJ(4));
t692 = -t776 * t711 + t733 * t712;
t729 = qJD(3) + qJD(4);
t770 = t692 * t729;
t728 = t731 ^ 2;
t769 = t728 * t739;
t735 = sin(qJ(1));
t738 = cos(qJ(1));
t717 = -t738 * g(1) - t735 * g(2);
t713 = -t739 * pkin(1) + qJDD(1) * qJ(2) + t717;
t761 = qJD(1) * qJD(2);
t759 = -t731 * g(3) - 0.2e1 * t730 * t761;
t683 = (-pkin(7) * qJDD(1) + t739 * t775 - t713) * t730 + t759;
t703 = -t730 * g(3) + (t713 + 0.2e1 * t761) * t731;
t760 = qJDD(1) * t731;
t685 = -pkin(2) * t769 + pkin(7) * t760 + t703;
t655 = t737 * t683 - t734 * t685;
t762 = t711 * qJD(3);
t701 = t750 * qJDD(1) + t762;
t631 = (-t701 + t762) * pkin(8) + (t711 * t712 + qJDD(3)) * pkin(3) + t655;
t656 = t734 * t683 + t737 * t685;
t700 = -t712 * qJD(3) + t749 * qJDD(1);
t706 = qJD(3) * pkin(3) - t712 * pkin(8);
t710 = t711 ^ 2;
t637 = -t710 * pkin(3) + t700 * pkin(8) - qJD(3) * t706 + t656;
t628 = t776 * t631 - t733 * t637;
t653 = -t692 * qJD(4) + t733 * t700 + t776 * t701;
t693 = t733 * t711 + t776 * t712;
t670 = t692 * mrSges(5,1) + t693 * mrSges(5,2);
t677 = -t729 * mrSges(5,2) - t692 * mrSges(5,3);
t679 = t692 * mrSges(6,1) - t729 * mrSges(6,3);
t726 = qJDD(3) + qJDD(4);
t669 = t692 * pkin(4) - t693 * qJ(5);
t725 = t729 ^ 2;
t625 = -t726 * pkin(4) - t725 * qJ(5) + t693 * t669 + qJDD(5) - t628;
t620 = (t692 * t693 - t726) * pkin(9) + (t653 + t770) * pkin(5) + t625;
t652 = t693 * qJD(4) - t776 * t700 + t733 * t701;
t681 = t693 * pkin(5) - t729 * pkin(9);
t688 = t692 ^ 2;
t727 = t730 ^ 2;
t716 = t735 * g(1) - t738 * g(2);
t754 = qJDD(2) - t716;
t699 = (-pkin(1) - t775) * qJDD(1) + (-qJ(2) + (-t727 - t728) * pkin(7)) * t739 + t754;
t643 = -t700 * pkin(3) - t710 * pkin(8) + t712 * t706 + t699;
t740 = (-t653 + t770) * qJ(5) + t643 + (pkin(4) * t729 + t777) * t693;
t621 = t740 + (pkin(4) + pkin(9)) * t652 - t693 * t681 - t688 * pkin(5);
t732 = sin(qJ(6));
t736 = cos(qJ(6));
t618 = t736 * t620 - t732 * t621;
t673 = t736 * t692 - t732 * t729;
t634 = t673 * qJD(6) + t732 * t652 + t736 * t726;
t651 = qJDD(6) + t653;
t674 = t732 * t692 + t736 * t729;
t654 = -t673 * mrSges(7,1) + t674 * mrSges(7,2);
t687 = qJD(6) + t693;
t657 = -t687 * mrSges(7,2) + t673 * mrSges(7,3);
t616 = m(7) * t618 + t651 * mrSges(7,1) - t634 * mrSges(7,3) - t674 * t654 + t687 * t657;
t619 = t732 * t620 + t736 * t621;
t633 = -t674 * qJD(6) + t736 * t652 - t732 * t726;
t658 = t687 * mrSges(7,1) - t674 * mrSges(7,3);
t617 = m(7) * t619 - t651 * mrSges(7,2) + t633 * mrSges(7,3) + t673 * t654 - t687 * t658;
t608 = t736 * t616 + t732 * t617;
t671 = -t692 * mrSges(6,2) - t693 * mrSges(6,3);
t746 = -m(6) * t625 - t653 * mrSges(6,1) - t693 * t671 - t608;
t606 = m(5) * t628 - t653 * mrSges(5,3) - t693 * t670 + (t677 - t679) * t729 + (mrSges(5,1) - mrSges(6,2)) * t726 + t746;
t629 = t733 * t631 + t776 * t637;
t678 = t729 * mrSges(5,1) - t693 * mrSges(5,3);
t745 = -t725 * pkin(4) + t726 * qJ(5) - t692 * t669 + t629;
t624 = t729 * t777 - t745;
t680 = t693 * mrSges(6,1) + t729 * mrSges(6,2);
t623 = -t652 * pkin(5) - t688 * pkin(9) + ((2 * qJD(5)) + t681) * t729 + t745;
t747 = -m(7) * t623 + t633 * mrSges(7,1) - t634 * mrSges(7,2) + t673 * t657 - t674 * t658;
t743 = -m(6) * t624 + t726 * mrSges(6,3) + t729 * t680 - t747;
t613 = m(5) * t629 - t726 * mrSges(5,2) - t729 * t678 + (-t670 - t671) * t692 + (-mrSges(5,3) - mrSges(6,1)) * t652 + t743;
t602 = t776 * t606 + t733 * t613;
t697 = -t711 * mrSges(4,1) + t712 * mrSges(4,2);
t704 = -qJD(3) * mrSges(4,2) + t711 * mrSges(4,3);
t600 = m(4) * t655 + qJDD(3) * mrSges(4,1) - t701 * mrSges(4,3) + qJD(3) * t704 - t712 * t697 + t602;
t705 = qJD(3) * mrSges(4,1) - t712 * mrSges(4,3);
t755 = -t733 * t606 + t776 * t613;
t601 = m(4) * t656 - qJDD(3) * mrSges(4,2) + t700 * mrSges(4,3) - qJD(3) * t705 + t711 * t697 + t755;
t595 = t737 * t600 + t734 * t601;
t702 = -t730 * t713 + t759;
t748 = mrSges(3,3) * qJDD(1) + t739 * (-mrSges(3,1) * t731 + t771);
t593 = m(3) * t702 - t748 * t730 + t595;
t756 = -t734 * t600 + t737 * t601;
t594 = m(3) * t703 + t748 * t731 + t756;
t757 = -t730 * t593 + t731 * t594;
t587 = m(2) * t717 - t739 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t757;
t709 = -qJDD(1) * pkin(1) - t739 * qJ(2) + t754;
t627 = t652 * pkin(4) + t740;
t767 = -t732 * t616 + t736 * t617;
t607 = m(6) * t627 - t652 * mrSges(6,2) - t653 * mrSges(6,3) - t692 * t679 - t693 * t680 + t767;
t744 = m(5) * t643 + t652 * mrSges(5,1) + t653 * mrSges(5,2) + t692 * t677 + t693 * t678 + t607;
t742 = m(4) * t699 - t700 * mrSges(4,1) + t701 * mrSges(4,2) - t711 * t704 + t712 * t705 + t744;
t741 = -m(3) * t709 + mrSges(3,1) * t760 - t742 + (t727 * t739 + t769) * mrSges(3,3);
t604 = t741 + (mrSges(2,1) - t771) * qJDD(1) - t739 * mrSges(2,2) + m(2) * t716;
t768 = t735 * t587 + t738 * t604;
t588 = t731 * t593 + t730 * t594;
t766 = t772 * t692 - t773 * t693 + t779 * t729;
t765 = t778 * t692 - t774 * t693 - t772 * t729;
t764 = -t774 * t692 + t780 * t693 + t773 * t729;
t751 = Ifges(3,5) * t730 + Ifges(3,6) * t731;
t763 = t739 * t751;
t758 = t738 * t587 - t735 * t604;
t753 = Ifges(3,1) * t730 + Ifges(3,4) * t731;
t752 = Ifges(3,4) * t730 + Ifges(3,2) * t731;
t691 = Ifges(4,1) * t712 + Ifges(4,4) * t711 + Ifges(4,5) * qJD(3);
t690 = Ifges(4,4) * t712 + Ifges(4,2) * t711 + Ifges(4,6) * qJD(3);
t689 = Ifges(4,5) * t712 + Ifges(4,6) * t711 + Ifges(4,3) * qJD(3);
t640 = Ifges(7,1) * t674 + Ifges(7,4) * t673 + Ifges(7,5) * t687;
t639 = Ifges(7,4) * t674 + Ifges(7,2) * t673 + Ifges(7,6) * t687;
t638 = Ifges(7,5) * t674 + Ifges(7,6) * t673 + Ifges(7,3) * t687;
t610 = mrSges(7,2) * t623 - mrSges(7,3) * t618 + Ifges(7,1) * t634 + Ifges(7,4) * t633 + Ifges(7,5) * t651 + t673 * t638 - t687 * t639;
t609 = -mrSges(7,1) * t623 + mrSges(7,3) * t619 + Ifges(7,4) * t634 + Ifges(7,2) * t633 + Ifges(7,6) * t651 - t674 * t638 + t687 * t640;
t596 = mrSges(6,1) * t625 + mrSges(7,1) * t618 + mrSges(5,2) * t643 - mrSges(7,2) * t619 - mrSges(5,3) * t628 - mrSges(6,3) * t627 + Ifges(7,5) * t634 + Ifges(7,6) * t633 + Ifges(7,3) * t651 + pkin(5) * t608 - qJ(5) * t607 + t674 * t639 - t673 * t640 + t765 * t729 + t773 * t726 + t766 * t692 + t780 * t653 - t774 * t652;
t589 = -mrSges(5,1) * t643 - mrSges(6,1) * t624 + mrSges(6,2) * t627 + mrSges(5,3) * t629 - pkin(4) * t607 - pkin(5) * t747 - pkin(9) * t767 - t736 * t609 - t732 * t610 - t778 * t652 + t774 * t653 + t766 * t693 + t772 * t726 + t764 * t729;
t584 = mrSges(4,2) * t699 - mrSges(4,3) * t655 + Ifges(4,1) * t701 + Ifges(4,4) * t700 + Ifges(4,5) * qJDD(3) - pkin(8) * t602 - qJD(3) * t690 - t733 * t589 + t776 * t596 + t711 * t689;
t583 = -mrSges(4,1) * t699 + mrSges(4,3) * t656 + Ifges(4,4) * t701 + Ifges(4,2) * t700 + Ifges(4,6) * qJDD(3) - pkin(3) * t744 + pkin(8) * t755 + qJD(3) * t691 + t776 * t589 + t733 * t596 - t712 * t689;
t582 = (-t730 * t752 + t731 * t753 + Ifges(2,5)) * t739 + (qJ(5) * mrSges(6,1) + t772) * t652 - t773 * t653 + (Ifges(2,6) - t751) * qJDD(1) - pkin(1) * t588 + t765 * t693 - Ifges(4,3) * qJDD(3) + mrSges(2,1) * g(3) - t736 * t610 + t732 * t609 + mrSges(2,3) * t717 + t711 * t691 - t712 * t690 - Ifges(4,6) * t700 - Ifges(4,5) * t701 - mrSges(3,1) * t702 + mrSges(3,2) * t703 - mrSges(4,1) * t655 + mrSges(4,2) * t656 - mrSges(5,1) * t628 + mrSges(5,2) * t629 + mrSges(6,3) * t624 - mrSges(6,2) * t625 + pkin(9) * t608 - pkin(4) * (-t729 * t679 + t746) - pkin(3) * t602 + (qJ(5) * t671 - t764) * t692 - pkin(2) * t595 + (pkin(4) * mrSges(6,2) + t779) * t726 - qJ(5) * t743;
t581 = mrSges(3,2) * t709 - mrSges(3,3) * t702 - pkin(7) * t595 + t753 * qJDD(1) - t734 * t583 + t737 * t584 + t731 * t763;
t580 = -mrSges(3,1) * t709 + mrSges(3,3) * t703 - pkin(2) * t742 + pkin(7) * t756 + t752 * qJDD(1) + t737 * t583 + t734 * t584 - t730 * t763;
t579 = -mrSges(2,2) * g(3) - mrSges(2,3) * t716 + Ifges(2,5) * qJDD(1) - t739 * Ifges(2,6) - qJ(2) * t588 - t730 * t580 + t731 * t581;
t1 = [-m(1) * g(1) + t758; -m(1) * g(2) + t768; (-m(1) - m(2)) * g(3) + t588; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t768 + t738 * t579 - t735 * t582; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t758 + t735 * t579 + t738 * t582; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t716 - mrSges(2,2) * t717 + t730 * t581 + t731 * t580 + pkin(1) * (-qJDD(1) * t771 + t741) + qJ(2) * t757;];
tauB  = t1;
