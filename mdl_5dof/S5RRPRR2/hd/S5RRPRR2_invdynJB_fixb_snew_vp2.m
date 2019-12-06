% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR2
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:59
% EndTime: 2019-12-05 18:28:12
% DurationCPUTime: 10.96s
% Computational Cost: add. (156687->315), mult. (373443->400), div. (0->0), fcn. (265933->10), ass. (0->125)
t762 = sin(qJ(2));
t766 = cos(qJ(2));
t783 = qJD(1) * qJD(2);
t741 = t762 * qJDD(1) + t766 * t783;
t763 = sin(qJ(1));
t767 = cos(qJ(1));
t748 = -t767 * g(1) - t763 * g(2);
t768 = qJD(1) ^ 2;
t736 = -t768 * pkin(1) + qJDD(1) * pkin(6) + t748;
t787 = t762 * t736;
t789 = pkin(2) * t768;
t703 = qJDD(2) * pkin(2) - t741 * qJ(3) - t787 + (qJ(3) * t783 + t762 * t789 - g(3)) * t766;
t722 = -t762 * g(3) + t766 * t736;
t742 = t766 * qJDD(1) - t762 * t783;
t785 = qJD(1) * t762;
t744 = qJD(2) * pkin(2) - qJ(3) * t785;
t757 = t766 ^ 2;
t704 = t742 * qJ(3) - qJD(2) * t744 - t757 * t789 + t722;
t758 = sin(pkin(9));
t759 = cos(pkin(9));
t731 = (t758 * t766 + t759 * t762) * qJD(1);
t681 = -0.2e1 * qJD(3) * t731 + t759 * t703 - t758 * t704;
t720 = t759 * t741 + t758 * t742;
t730 = (-t758 * t762 + t759 * t766) * qJD(1);
t670 = (qJD(2) * t730 - t720) * pkin(7) + (t730 * t731 + qJDD(2)) * pkin(3) + t681;
t682 = 0.2e1 * qJD(3) * t730 + t758 * t703 + t759 * t704;
t719 = -t758 * t741 + t759 * t742;
t725 = qJD(2) * pkin(3) - t731 * pkin(7);
t729 = t730 ^ 2;
t672 = -t729 * pkin(3) + t719 * pkin(7) - qJD(2) * t725 + t682;
t761 = sin(qJ(4));
t765 = cos(qJ(4));
t657 = t765 * t670 - t761 * t672;
t713 = t765 * t730 - t761 * t731;
t689 = t713 * qJD(4) + t761 * t719 + t765 * t720;
t714 = t761 * t730 + t765 * t731;
t754 = qJDD(2) + qJDD(4);
t755 = qJD(2) + qJD(4);
t654 = (t713 * t755 - t689) * pkin(8) + (t713 * t714 + t754) * pkin(4) + t657;
t658 = t761 * t670 + t765 * t672;
t688 = -t714 * qJD(4) + t765 * t719 - t761 * t720;
t708 = t755 * pkin(4) - t714 * pkin(8);
t709 = t713 ^ 2;
t655 = -t709 * pkin(4) + t688 * pkin(8) - t755 * t708 + t658;
t760 = sin(qJ(5));
t764 = cos(qJ(5));
t652 = t764 * t654 - t760 * t655;
t697 = t764 * t713 - t760 * t714;
t666 = t697 * qJD(5) + t760 * t688 + t764 * t689;
t698 = t760 * t713 + t764 * t714;
t678 = -t697 * mrSges(6,1) + t698 * mrSges(6,2);
t752 = qJD(5) + t755;
t690 = -t752 * mrSges(6,2) + t697 * mrSges(6,3);
t751 = qJDD(5) + t754;
t649 = m(6) * t652 + t751 * mrSges(6,1) - t666 * mrSges(6,3) - t698 * t678 + t752 * t690;
t653 = t760 * t654 + t764 * t655;
t665 = -t698 * qJD(5) + t764 * t688 - t760 * t689;
t691 = t752 * mrSges(6,1) - t698 * mrSges(6,3);
t650 = m(6) * t653 - t751 * mrSges(6,2) + t665 * mrSges(6,3) + t697 * t678 - t752 * t691;
t639 = t764 * t649 + t760 * t650;
t699 = -t713 * mrSges(5,1) + t714 * mrSges(5,2);
t706 = -t755 * mrSges(5,2) + t713 * mrSges(5,3);
t636 = m(5) * t657 + t754 * mrSges(5,1) - t689 * mrSges(5,3) - t714 * t699 + t755 * t706 + t639;
t707 = t755 * mrSges(5,1) - t714 * mrSges(5,3);
t778 = -t760 * t649 + t764 * t650;
t637 = m(5) * t658 - t754 * mrSges(5,2) + t688 * mrSges(5,3) + t713 * t699 - t755 * t707 + t778;
t632 = t765 * t636 + t761 * t637;
t717 = -t730 * mrSges(4,1) + t731 * mrSges(4,2);
t723 = -qJD(2) * mrSges(4,2) + t730 * mrSges(4,3);
t630 = m(4) * t681 + qJDD(2) * mrSges(4,1) - t720 * mrSges(4,3) + qJD(2) * t723 - t731 * t717 + t632;
t724 = qJD(2) * mrSges(4,1) - t731 * mrSges(4,3);
t779 = -t761 * t636 + t765 * t637;
t631 = m(4) * t682 - qJDD(2) * mrSges(4,2) + t719 * mrSges(4,3) - qJD(2) * t724 + t730 * t717 + t779;
t624 = t759 * t630 + t758 * t631;
t711 = Ifges(4,4) * t731 + Ifges(4,2) * t730 + Ifges(4,6) * qJD(2);
t712 = Ifges(4,1) * t731 + Ifges(4,4) * t730 + Ifges(4,5) * qJD(2);
t721 = -t766 * g(3) - t787;
t733 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t762 + Ifges(3,2) * t766) * qJD(1);
t734 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t762 + Ifges(3,4) * t766) * qJD(1);
t693 = Ifges(5,4) * t714 + Ifges(5,2) * t713 + Ifges(5,6) * t755;
t694 = Ifges(5,1) * t714 + Ifges(5,4) * t713 + Ifges(5,5) * t755;
t674 = Ifges(6,4) * t698 + Ifges(6,2) * t697 + Ifges(6,6) * t752;
t675 = Ifges(6,1) * t698 + Ifges(6,4) * t697 + Ifges(6,5) * t752;
t772 = -mrSges(6,1) * t652 + mrSges(6,2) * t653 - Ifges(6,5) * t666 - Ifges(6,6) * t665 - Ifges(6,3) * t751 - t698 * t674 + t697 * t675;
t771 = -mrSges(5,1) * t657 + mrSges(5,2) * t658 - Ifges(5,5) * t689 - Ifges(5,6) * t688 - Ifges(5,3) * t754 - pkin(4) * t639 - t714 * t693 + t713 * t694 + t772;
t790 = mrSges(3,1) * t721 + mrSges(4,1) * t681 - mrSges(3,2) * t722 - mrSges(4,2) * t682 + Ifges(3,5) * t741 + Ifges(4,5) * t720 + Ifges(3,6) * t742 + Ifges(4,6) * t719 + pkin(2) * t624 + pkin(3) * t632 + (t762 * t733 - t766 * t734) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t731 * t711 - t730 * t712 - t771;
t740 = (-mrSges(3,1) * t766 + mrSges(3,2) * t762) * qJD(1);
t784 = qJD(1) * t766;
t746 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t784;
t622 = m(3) * t721 + qJDD(2) * mrSges(3,1) - t741 * mrSges(3,3) + qJD(2) * t746 - t740 * t785 + t624;
t745 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t785;
t780 = -t758 * t630 + t759 * t631;
t623 = m(3) * t722 - qJDD(2) * mrSges(3,2) + t742 * mrSges(3,3) - qJD(2) * t745 + t740 * t784 + t780;
t781 = -t762 * t622 + t766 * t623;
t614 = m(2) * t748 - t768 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t781;
t747 = t763 * g(1) - t767 * g(2);
t775 = -qJDD(1) * pkin(1) - t747;
t705 = -t742 * pkin(2) + qJDD(3) + t744 * t785 + (-qJ(3) * t757 - pkin(6)) * t768 + t775;
t680 = -t719 * pkin(3) - t729 * pkin(7) + t731 * t725 + t705;
t660 = -t688 * pkin(4) - t709 * pkin(8) + t714 * t708 + t680;
t777 = m(6) * t660 - t665 * mrSges(6,1) + t666 * mrSges(6,2) - t697 * t690 + t698 * t691;
t773 = m(5) * t680 - t688 * mrSges(5,1) + t689 * mrSges(5,2) - t713 * t706 + t714 * t707 + t777;
t645 = m(4) * t705 - t719 * mrSges(4,1) + t720 * mrSges(4,2) - t730 * t723 + t731 * t724 + t773;
t735 = -t768 * pkin(6) + t775;
t770 = -m(3) * t735 + t742 * mrSges(3,1) - t741 * mrSges(3,2) - t745 * t785 + t746 * t784 - t645;
t643 = m(2) * t747 + qJDD(1) * mrSges(2,1) - t768 * mrSges(2,2) + t770;
t786 = t763 * t614 + t767 * t643;
t616 = t766 * t622 + t762 * t623;
t782 = t767 * t614 - t763 * t643;
t673 = Ifges(6,5) * t698 + Ifges(6,6) * t697 + Ifges(6,3) * t752;
t640 = -mrSges(6,1) * t660 + mrSges(6,3) * t653 + Ifges(6,4) * t666 + Ifges(6,2) * t665 + Ifges(6,6) * t751 - t698 * t673 + t752 * t675;
t641 = mrSges(6,2) * t660 - mrSges(6,3) * t652 + Ifges(6,1) * t666 + Ifges(6,4) * t665 + Ifges(6,5) * t751 + t697 * t673 - t752 * t674;
t692 = Ifges(5,5) * t714 + Ifges(5,6) * t713 + Ifges(5,3) * t755;
t625 = -mrSges(5,1) * t680 + mrSges(5,3) * t658 + Ifges(5,4) * t689 + Ifges(5,2) * t688 + Ifges(5,6) * t754 - pkin(4) * t777 + pkin(8) * t778 + t764 * t640 + t760 * t641 - t714 * t692 + t755 * t694;
t626 = mrSges(5,2) * t680 - mrSges(5,3) * t657 + Ifges(5,1) * t689 + Ifges(5,4) * t688 + Ifges(5,5) * t754 - pkin(8) * t639 - t760 * t640 + t764 * t641 + t713 * t692 - t755 * t693;
t710 = Ifges(4,5) * t731 + Ifges(4,6) * t730 + Ifges(4,3) * qJD(2);
t617 = -mrSges(4,1) * t705 + mrSges(4,3) * t682 + Ifges(4,4) * t720 + Ifges(4,2) * t719 + Ifges(4,6) * qJDD(2) - pkin(3) * t773 + pkin(7) * t779 + qJD(2) * t712 + t765 * t625 + t761 * t626 - t731 * t710;
t618 = mrSges(4,2) * t705 - mrSges(4,3) * t681 + Ifges(4,1) * t720 + Ifges(4,4) * t719 + Ifges(4,5) * qJDD(2) - pkin(7) * t632 - qJD(2) * t711 - t761 * t625 + t765 * t626 + t730 * t710;
t732 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t762 + Ifges(3,6) * t766) * qJD(1);
t608 = -mrSges(3,1) * t735 + mrSges(3,3) * t722 + Ifges(3,4) * t741 + Ifges(3,2) * t742 + Ifges(3,6) * qJDD(2) - pkin(2) * t645 + qJ(3) * t780 + qJD(2) * t734 + t759 * t617 + t758 * t618 - t732 * t785;
t610 = mrSges(3,2) * t735 - mrSges(3,3) * t721 + Ifges(3,1) * t741 + Ifges(3,4) * t742 + Ifges(3,5) * qJDD(2) - qJ(3) * t624 - qJD(2) * t733 - t758 * t617 + t759 * t618 + t732 * t784;
t774 = mrSges(2,1) * t747 - mrSges(2,2) * t748 + Ifges(2,3) * qJDD(1) + pkin(1) * t770 + pkin(6) * t781 + t766 * t608 + t762 * t610;
t611 = mrSges(2,1) * g(3) + mrSges(2,3) * t748 + t768 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t616 - t790;
t606 = -mrSges(2,2) * g(3) - mrSges(2,3) * t747 + Ifges(2,5) * qJDD(1) - t768 * Ifges(2,6) - pkin(6) * t616 - t762 * t608 + t766 * t610;
t1 = [-m(1) * g(1) + t782; -m(1) * g(2) + t786; (-m(1) - m(2)) * g(3) + t616; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t786 + t767 * t606 - t763 * t611; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t782 + t763 * t606 + t767 * t611; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t774; t774; t790; t645; -t771; -t772;];
tauJB = t1;
