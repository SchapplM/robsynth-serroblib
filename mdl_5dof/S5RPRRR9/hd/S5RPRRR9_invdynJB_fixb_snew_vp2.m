% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:51
% EndTime: 2019-12-31 19:07:58
% DurationCPUTime: 6.96s
% Computational Cost: add. (102178->292), mult. (246148->366), div. (0->0), fcn. (185860->10), ass. (0->127)
t751 = qJD(1) ^ 2;
t742 = cos(pkin(9));
t780 = pkin(2) * t742;
t741 = sin(pkin(9));
t779 = mrSges(3,2) * t741;
t738 = t742 ^ 2;
t778 = t738 * t751;
t746 = sin(qJ(1));
t750 = cos(qJ(1));
t727 = -t750 * g(1) - t746 * g(2);
t722 = -t751 * pkin(1) + qJDD(1) * qJ(2) + t727;
t774 = qJD(1) * qJD(2);
t772 = -t742 * g(3) - 0.2e1 * t741 * t774;
t697 = (-pkin(6) * qJDD(1) + t751 * t780 - t722) * t741 + t772;
t713 = -t741 * g(3) + (t722 + 0.2e1 * t774) * t742;
t773 = qJDD(1) * t742;
t698 = -pkin(2) * t778 + pkin(6) * t773 + t713;
t745 = sin(qJ(3));
t749 = cos(qJ(3));
t679 = t749 * t697 - t745 * t698;
t762 = t741 * t749 + t742 * t745;
t761 = -t741 * t745 + t742 * t749;
t720 = t761 * qJD(1);
t775 = t720 * qJD(3);
t711 = t762 * qJDD(1) + t775;
t721 = t762 * qJD(1);
t658 = (-t711 + t775) * pkin(7) + (t720 * t721 + qJDD(3)) * pkin(3) + t679;
t680 = t745 * t697 + t749 * t698;
t710 = -t721 * qJD(3) + t761 * qJDD(1);
t716 = qJD(3) * pkin(3) - t721 * pkin(7);
t719 = t720 ^ 2;
t663 = -t719 * pkin(3) + t710 * pkin(7) - qJD(3) * t716 + t680;
t744 = sin(qJ(4));
t748 = cos(qJ(4));
t656 = t744 * t658 + t748 * t663;
t704 = t744 * t720 + t748 * t721;
t676 = -t704 * qJD(4) + t748 * t710 - t744 * t711;
t703 = t748 * t720 - t744 * t721;
t688 = -t703 * mrSges(5,1) + t704 * mrSges(5,2);
t739 = qJD(3) + qJD(4);
t695 = t739 * mrSges(5,1) - t704 * mrSges(5,3);
t736 = qJDD(3) + qJDD(4);
t689 = -t703 * pkin(4) - t704 * pkin(8);
t735 = t739 ^ 2;
t652 = -t735 * pkin(4) + t736 * pkin(8) + t703 * t689 + t656;
t737 = t741 ^ 2;
t726 = t746 * g(1) - t750 * g(2);
t766 = qJDD(2) - t726;
t709 = (-pkin(1) - t780) * qJDD(1) + (-qJ(2) + (-t737 - t738) * pkin(6)) * t751 + t766;
t670 = -t710 * pkin(3) - t719 * pkin(7) + t721 * t716 + t709;
t677 = t703 * qJD(4) + t744 * t710 + t748 * t711;
t653 = (-t703 * t739 - t677) * pkin(8) + (t704 * t739 - t676) * pkin(4) + t670;
t743 = sin(qJ(5));
t747 = cos(qJ(5));
t649 = -t743 * t652 + t747 * t653;
t690 = -t743 * t704 + t747 * t739;
t661 = t690 * qJD(5) + t747 * t677 + t743 * t736;
t675 = qJDD(5) - t676;
t691 = t747 * t704 + t743 * t739;
t678 = -t690 * mrSges(6,1) + t691 * mrSges(6,2);
t699 = qJD(5) - t703;
t681 = -t699 * mrSges(6,2) + t690 * mrSges(6,3);
t645 = m(6) * t649 + t675 * mrSges(6,1) - t661 * mrSges(6,3) - t691 * t678 + t699 * t681;
t650 = t747 * t652 + t743 * t653;
t660 = -t691 * qJD(5) - t743 * t677 + t747 * t736;
t682 = t699 * mrSges(6,1) - t691 * mrSges(6,3);
t646 = m(6) * t650 - t675 * mrSges(6,2) + t660 * mrSges(6,3) + t690 * t678 - t699 * t682;
t767 = -t743 * t645 + t747 * t646;
t632 = m(5) * t656 - t736 * mrSges(5,2) + t676 * mrSges(5,3) + t703 * t688 - t739 * t695 + t767;
t655 = t748 * t658 - t744 * t663;
t694 = -t739 * mrSges(5,2) + t703 * mrSges(5,3);
t651 = -t736 * pkin(4) - t735 * pkin(8) + t704 * t689 - t655;
t757 = -m(6) * t651 + t660 * mrSges(6,1) - t661 * mrSges(6,2) + t690 * t681 - t691 * t682;
t641 = m(5) * t655 + t736 * mrSges(5,1) - t677 * mrSges(5,3) - t704 * t688 + t739 * t694 + t757;
t625 = t744 * t632 + t748 * t641;
t707 = -t720 * mrSges(4,1) + t721 * mrSges(4,2);
t714 = -qJD(3) * mrSges(4,2) + t720 * mrSges(4,3);
t623 = m(4) * t679 + qJDD(3) * mrSges(4,1) - t711 * mrSges(4,3) + qJD(3) * t714 - t721 * t707 + t625;
t715 = qJD(3) * mrSges(4,1) - t721 * mrSges(4,3);
t768 = t748 * t632 - t744 * t641;
t624 = m(4) * t680 - qJDD(3) * mrSges(4,2) + t710 * mrSges(4,3) - qJD(3) * t715 + t720 * t707 + t768;
t617 = t749 * t623 + t745 * t624;
t712 = -t741 * t722 + t772;
t760 = mrSges(3,3) * qJDD(1) + t751 * (-mrSges(3,1) * t742 + t779);
t615 = m(3) * t712 - t760 * t741 + t617;
t769 = -t745 * t623 + t749 * t624;
t616 = m(3) * t713 + t760 * t742 + t769;
t770 = -t741 * t615 + t742 * t616;
t608 = m(2) * t727 - t751 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t770;
t718 = -qJDD(1) * pkin(1) - t751 * qJ(2) + t766;
t634 = t747 * t645 + t743 * t646;
t759 = m(5) * t670 - t676 * mrSges(5,1) + t677 * mrSges(5,2) - t703 * t694 + t704 * t695 + t634;
t755 = m(4) * t709 - t710 * mrSges(4,1) + t711 * mrSges(4,2) - t720 * t714 + t721 * t715 + t759;
t753 = -m(3) * t718 + mrSges(3,1) * t773 - t755 + (t737 * t751 + t778) * mrSges(3,3);
t627 = t753 + (mrSges(2,1) - t779) * qJDD(1) - t751 * mrSges(2,2) + m(2) * t726;
t777 = t746 * t608 + t750 * t627;
t610 = t742 * t615 + t741 * t616;
t763 = Ifges(3,5) * t741 + Ifges(3,6) * t742;
t776 = t751 * t763;
t771 = t750 * t608 - t746 * t627;
t765 = Ifges(3,1) * t741 + Ifges(3,4) * t742;
t764 = Ifges(3,4) * t741 + Ifges(3,2) * t742;
t664 = Ifges(6,5) * t691 + Ifges(6,6) * t690 + Ifges(6,3) * t699;
t666 = Ifges(6,1) * t691 + Ifges(6,4) * t690 + Ifges(6,5) * t699;
t638 = -mrSges(6,1) * t651 + mrSges(6,3) * t650 + Ifges(6,4) * t661 + Ifges(6,2) * t660 + Ifges(6,6) * t675 - t691 * t664 + t699 * t666;
t665 = Ifges(6,4) * t691 + Ifges(6,2) * t690 + Ifges(6,6) * t699;
t639 = mrSges(6,2) * t651 - mrSges(6,3) * t649 + Ifges(6,1) * t661 + Ifges(6,4) * t660 + Ifges(6,5) * t675 + t690 * t664 - t699 * t665;
t683 = Ifges(5,5) * t704 + Ifges(5,6) * t703 + Ifges(5,3) * t739;
t684 = Ifges(5,4) * t704 + Ifges(5,2) * t703 + Ifges(5,6) * t739;
t618 = mrSges(5,2) * t670 - mrSges(5,3) * t655 + Ifges(5,1) * t677 + Ifges(5,4) * t676 + Ifges(5,5) * t736 - pkin(8) * t634 - t743 * t638 + t747 * t639 + t703 * t683 - t739 * t684;
t685 = Ifges(5,1) * t704 + Ifges(5,4) * t703 + Ifges(5,5) * t739;
t754 = mrSges(6,1) * t649 - mrSges(6,2) * t650 + Ifges(6,5) * t661 + Ifges(6,6) * t660 + Ifges(6,3) * t675 + t691 * t665 - t690 * t666;
t619 = -mrSges(5,1) * t670 + mrSges(5,3) * t656 + Ifges(5,4) * t677 + Ifges(5,2) * t676 + Ifges(5,6) * t736 - pkin(4) * t634 - t704 * t683 + t739 * t685 - t754;
t700 = Ifges(4,5) * t721 + Ifges(4,6) * t720 + Ifges(4,3) * qJD(3);
t702 = Ifges(4,1) * t721 + Ifges(4,4) * t720 + Ifges(4,5) * qJD(3);
t605 = -mrSges(4,1) * t709 + mrSges(4,3) * t680 + Ifges(4,4) * t711 + Ifges(4,2) * t710 + Ifges(4,6) * qJDD(3) - pkin(3) * t759 + pkin(7) * t768 + qJD(3) * t702 + t744 * t618 + t748 * t619 - t721 * t700;
t701 = Ifges(4,4) * t721 + Ifges(4,2) * t720 + Ifges(4,6) * qJD(3);
t611 = mrSges(4,2) * t709 - mrSges(4,3) * t679 + Ifges(4,1) * t711 + Ifges(4,4) * t710 + Ifges(4,5) * qJDD(3) - pkin(7) * t625 - qJD(3) * t701 + t748 * t618 - t744 * t619 + t720 * t700;
t601 = -mrSges(3,1) * t718 + mrSges(3,3) * t713 - pkin(2) * t755 + pkin(6) * t769 + t764 * qJDD(1) + t749 * t605 + t745 * t611 - t741 * t776;
t603 = mrSges(3,2) * t718 - mrSges(3,3) * t712 - pkin(6) * t617 + t765 * qJDD(1) - t745 * t605 + t749 * t611 + t742 * t776;
t629 = qJDD(1) * t779 - t753;
t758 = mrSges(2,1) * t726 - mrSges(2,2) * t727 + Ifges(2,3) * qJDD(1) - pkin(1) * t629 + qJ(2) * t770 + t742 * t601 + t741 * t603;
t756 = -mrSges(5,1) * t655 + mrSges(5,2) * t656 - Ifges(5,5) * t677 - Ifges(5,6) * t676 - Ifges(5,3) * t736 - pkin(4) * t757 - pkin(8) * t767 - t747 * t638 - t743 * t639 - t704 * t684 + t703 * t685;
t752 = mrSges(4,1) * t679 - mrSges(4,2) * t680 + Ifges(4,5) * t711 + Ifges(4,6) * t710 + Ifges(4,3) * qJDD(3) + pkin(3) * t625 + t721 * t701 - t720 * t702 - t756;
t604 = -t752 - pkin(2) * t617 + mrSges(2,1) * g(3) - pkin(1) * t610 + (Ifges(2,6) - t763) * qJDD(1) + mrSges(2,3) * t727 - mrSges(3,1) * t712 + mrSges(3,2) * t713 + (-t741 * t764 + t742 * t765 + Ifges(2,5)) * t751;
t599 = -mrSges(2,2) * g(3) - mrSges(2,3) * t726 + Ifges(2,5) * qJDD(1) - t751 * Ifges(2,6) - qJ(2) * t610 - t741 * t601 + t742 * t603;
t1 = [-m(1) * g(1) + t771; -m(1) * g(2) + t777; (-m(1) - m(2)) * g(3) + t610; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t777 + t750 * t599 - t746 * t604; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t771 + t746 * t599 + t750 * t604; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t758; t758; t629; t752; -t756; t754;];
tauJB = t1;
