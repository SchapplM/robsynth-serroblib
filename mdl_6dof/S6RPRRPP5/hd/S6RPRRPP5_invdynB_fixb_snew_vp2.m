% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPP5
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:36:38
% EndTime: 2019-05-05 21:36:46
% DurationCPUTime: 5.62s
% Computational Cost: add. (57036->321), mult. (133564->377), div. (0->0), fcn. (95722->8), ass. (0->132)
t778 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t756 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t755 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t777 = Ifges(5,2) + Ifges(6,3) + Ifges(7,2);
t754 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t776 = -Ifges(5,3) - Ifges(6,2) - Ifges(7,3);
t728 = qJD(1) ^ 2;
t721 = sin(pkin(9));
t722 = cos(pkin(9));
t724 = sin(qJ(3));
t772 = cos(qJ(3));
t734 = t772 * t721 + t722 * t724;
t749 = t722 * t772;
t760 = t721 * qJD(1);
t702 = qJD(1) * t749 - t724 * t760;
t762 = t702 * qJD(3);
t689 = qJDD(1) * t734 + t762;
t703 = t734 * qJD(1);
t723 = sin(qJ(4));
t771 = cos(qJ(4));
t693 = -t771 * qJD(3) + t723 * t703;
t655 = -t693 * qJD(4) + t723 * qJDD(3) + t771 * t689;
t725 = sin(qJ(1));
t726 = cos(qJ(1));
t708 = -t726 * g(1) - t725 * g(2);
t704 = -t728 * pkin(1) + qJDD(1) * qJ(2) + t708;
t759 = qJD(1) * qJD(2);
t748 = -t722 * g(3) - 0.2e1 * t721 * t759;
t770 = pkin(2) * t722;
t674 = (-pkin(7) * qJDD(1) + t728 * t770 - t704) * t721 + t748;
t691 = -t721 * g(3) + (t704 + 0.2e1 * t759) * t722;
t757 = qJDD(1) * t722;
t720 = t722 ^ 2;
t765 = t720 * t728;
t675 = -pkin(2) * t765 + pkin(7) * t757 + t691;
t629 = t772 * t674 - t724 * t675;
t686 = -t702 * pkin(3) - t703 * pkin(8);
t727 = qJD(3) ^ 2;
t733 = qJDD(3) * pkin(3) + t727 * pkin(8) - t703 * t686 + t629;
t700 = qJD(4) - t702;
t767 = t693 * t700;
t775 = (-t655 + t767) * qJ(5) - t733;
t694 = t723 * qJD(3) + t771 * t703;
t774 = -0.2e1 * t694;
t773 = 2 * qJD(5);
t769 = -mrSges(5,3) - mrSges(6,2);
t768 = mrSges(3,2) * t721;
t666 = t700 * mrSges(7,2) + t693 * mrSges(7,3);
t766 = t700 * t666;
t630 = t724 * t674 + t772 * t675;
t683 = -t702 * mrSges(4,1) + t703 * mrSges(4,2);
t758 = qJDD(1) * t721;
t761 = t703 * qJD(3);
t688 = qJDD(1) * t749 - t724 * t758 - t761;
t696 = qJD(3) * mrSges(4,1) - t703 * mrSges(4,3);
t626 = -t727 * pkin(3) + qJDD(3) * pkin(8) + t702 * t686 + t630;
t719 = t721 ^ 2;
t707 = t725 * g(1) - t726 * g(2);
t741 = qJDD(2) - t707;
t687 = (-pkin(1) - t770) * qJDD(1) + (-qJ(2) + (-t719 - t720) * pkin(7)) * t728 + t741;
t628 = (-t689 - t762) * pkin(8) + (-t688 + t761) * pkin(3) + t687;
t621 = -t723 * t626 + t771 * t628;
t667 = -t700 * mrSges(5,2) - t693 * mrSges(5,3);
t685 = qJDD(4) - t688;
t658 = t693 * pkin(4) - t694 * qJ(5);
t699 = t700 ^ 2;
t619 = -t685 * pkin(4) - t699 * qJ(5) + t694 * t658 + qJDD(5) - t621;
t665 = -t693 * mrSges(6,2) + t700 * mrSges(6,3);
t612 = qJD(6) * t774 + (-t655 - t767) * qJ(6) + (t693 * t694 - t685) * pkin(5) + t619;
t660 = -t693 * mrSges(7,1) + t694 * mrSges(7,2);
t742 = -m(7) * t612 + t655 * mrSges(7,3) + t694 * t660;
t732 = -m(6) * t619 + t685 * mrSges(6,1) + t700 * t665 + t742;
t659 = t693 * mrSges(6,1) - t694 * mrSges(6,3);
t763 = -t693 * mrSges(5,1) - t694 * mrSges(5,2) - t659;
t606 = m(5) * t621 + (t666 + t667) * t700 + t763 * t694 + (mrSges(5,1) + mrSges(7,1)) * t685 + t769 * t655 + t732;
t622 = t771 * t626 + t723 * t628;
t654 = t694 * qJD(4) - t771 * qJDD(3) + t723 * t689;
t669 = -t700 * mrSges(7,1) - t694 * mrSges(7,3);
t670 = t700 * mrSges(5,1) - t694 * mrSges(5,3);
t618 = -t699 * pkin(4) + t685 * qJ(5) - t693 * t658 + t700 * t773 + t622;
t671 = -t700 * mrSges(6,1) + t694 * mrSges(6,2);
t668 = -t700 * pkin(5) - t694 * qJ(6);
t692 = t693 ^ 2;
t614 = -t692 * pkin(5) + t654 * qJ(6) + 0.2e1 * qJD(6) * t693 + t700 * t668 + t618;
t753 = m(7) * t614 + t654 * mrSges(7,3) + t693 * t660;
t735 = m(6) * t618 + t685 * mrSges(6,3) + t700 * t671 + t753;
t609 = m(5) * t622 + (t669 - t670) * t700 + t763 * t693 + (-mrSges(5,2) + mrSges(7,2)) * t685 + t769 * t654 + t735;
t744 = -t723 * t606 + t771 * t609;
t601 = m(4) * t630 - qJDD(3) * mrSges(4,2) + t688 * mrSges(4,3) - qJD(3) * t696 + t702 * t683 + t744;
t695 = -qJD(3) * mrSges(4,2) + t702 * mrSges(4,3);
t620 = qJD(5) * t774 + (t694 * t700 + t654) * pkin(4) + t775;
t616 = -t692 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t654 + (-pkin(4) * t700 + t668 + t773) * t694 - t775;
t737 = -m(7) * t616 + t654 * mrSges(7,1) - t655 * mrSges(7,2) + t693 * t666 - t694 * t669;
t610 = m(6) * t620 + t654 * mrSges(6,1) - t655 * mrSges(6,3) + t693 * t665 - t694 * t671 + t737;
t729 = m(5) * t733 - t654 * mrSges(5,1) - t655 * mrSges(5,2) - t693 * t667 - t694 * t670 - t610;
t604 = m(4) * t629 + qJDD(3) * mrSges(4,1) - t689 * mrSges(4,3) + qJD(3) * t695 - t703 * t683 + t729;
t595 = t724 * t601 + t772 * t604;
t690 = -t721 * t704 + t748;
t736 = mrSges(3,3) * qJDD(1) + t728 * (-mrSges(3,1) * t722 + t768);
t593 = m(3) * t690 - t721 * t736 + t595;
t745 = t772 * t601 - t724 * t604;
t594 = m(3) * t691 + t722 * t736 + t745;
t746 = -t721 * t593 + t722 * t594;
t586 = m(2) * t708 - t728 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t746;
t701 = -qJDD(1) * pkin(1) - t728 * qJ(2) + t741;
t602 = t771 * t606 + t723 * t609;
t731 = m(4) * t687 - t688 * mrSges(4,1) + t689 * mrSges(4,2) - t702 * t695 + t703 * t696 + t602;
t730 = -m(3) * t701 + mrSges(3,1) * t757 - t731 + (t719 * t728 + t765) * mrSges(3,3);
t598 = (mrSges(2,1) - t768) * qJDD(1) + t730 - t728 * mrSges(2,2) + m(2) * t707;
t764 = t725 * t586 + t726 * t598;
t587 = t722 * t593 + t721 * t594;
t752 = t754 * t693 - t755 * t694 + t776 * t700;
t751 = t777 * t693 - t756 * t694 - t754 * t700;
t750 = t756 * t693 - t778 * t694 - t755 * t700;
t747 = t726 * t586 - t725 * t598;
t740 = Ifges(3,1) * t721 + Ifges(3,4) * t722;
t739 = Ifges(3,4) * t721 + Ifges(3,2) * t722;
t738 = Ifges(3,5) * t721 + Ifges(3,6) * t722;
t706 = t738 * qJD(1);
t678 = Ifges(4,1) * t703 + Ifges(4,4) * t702 + Ifges(4,5) * qJD(3);
t677 = Ifges(4,4) * t703 + Ifges(4,2) * t702 + Ifges(4,6) * qJD(3);
t676 = Ifges(4,5) * t703 + Ifges(4,6) * t702 + Ifges(4,3) * qJD(3);
t611 = -t685 * mrSges(7,1) - t742 - t766;
t596 = -mrSges(5,2) * t733 + mrSges(6,2) * t619 + mrSges(7,2) * t616 - mrSges(5,3) * t621 - mrSges(6,3) * t620 - mrSges(7,3) * t612 - qJ(5) * t610 - qJ(6) * t611 - t756 * t654 + t778 * t655 + t755 * t685 + t752 * t693 + t751 * t700;
t589 = mrSges(5,1) * t733 + mrSges(5,3) * t622 - mrSges(6,1) * t620 + mrSges(6,2) * t618 + mrSges(7,1) * t616 - mrSges(7,3) * t614 - pkin(5) * t737 - qJ(6) * t753 - pkin(4) * t610 + (-qJ(6) * t669 - t750) * t700 + t752 * t694 + (-qJ(6) * mrSges(7,2) + t754) * t685 + t756 * t655 - t777 * t654;
t588 = Ifges(4,6) * qJDD(3) - t703 * t676 - mrSges(4,1) * t687 + Ifges(4,2) * t688 + Ifges(4,4) * t689 + qJD(3) * t678 + mrSges(4,3) * t630 + mrSges(5,2) * t622 - mrSges(6,3) * t618 + mrSges(6,1) * t619 - mrSges(5,1) * t621 - mrSges(7,2) * t614 + pkin(5) * t611 + mrSges(7,1) * t612 - pkin(3) * t602 - pkin(4) * (t732 + t766) - qJ(5) * (t700 * t669 + t735) + (pkin(4) * t659 + t751) * t694 + (qJ(5) * t659 + t750) * t693 + (-pkin(4) * mrSges(7,1) - qJ(5) * mrSges(7,2) + t776) * t685 + (pkin(4) * mrSges(6,2) - t755) * t655 + (qJ(5) * mrSges(6,2) + t754) * t654;
t583 = mrSges(4,2) * t687 - mrSges(4,3) * t629 + Ifges(4,1) * t689 + Ifges(4,4) * t688 + Ifges(4,5) * qJDD(3) - pkin(8) * t602 - qJD(3) * t677 - t723 * t589 + t771 * t596 + t702 * t676;
t582 = t722 * qJD(1) * t706 + mrSges(3,2) * t701 - mrSges(3,3) * t690 - pkin(7) * t595 + t740 * qJDD(1) + t772 * t583 - t724 * t588;
t581 = mrSges(2,1) * g(3) - pkin(1) * t587 + mrSges(2,3) * t708 - pkin(2) * t595 - mrSges(3,1) * t690 + mrSges(3,2) * t691 - pkin(8) * t744 - t723 * t596 - t771 * t589 - pkin(3) * t729 - mrSges(4,1) * t629 + mrSges(4,2) * t630 - Ifges(4,5) * t689 - Ifges(4,6) * t688 - Ifges(4,3) * qJDD(3) - t703 * t677 + t702 * t678 + (Ifges(2,6) - t738) * qJDD(1) + (-t721 * t739 + t722 * t740 + Ifges(2,5)) * t728;
t580 = -mrSges(3,1) * t701 + mrSges(3,3) * t691 - pkin(2) * t731 + pkin(7) * t745 + t739 * qJDD(1) + t724 * t583 + t772 * t588 - t706 * t760;
t579 = -mrSges(2,2) * g(3) - mrSges(2,3) * t707 + Ifges(2,5) * qJDD(1) - t728 * Ifges(2,6) - qJ(2) * t587 - t721 * t580 + t722 * t582;
t1 = [-m(1) * g(1) + t747; -m(1) * g(2) + t764; (-m(1) - m(2)) * g(3) + t587; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t764 + t726 * t579 - t725 * t581; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t747 + t725 * t579 + t726 * t581; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t707 - mrSges(2,2) * t708 + t721 * t582 + t722 * t580 + pkin(1) * (-mrSges(3,2) * t758 + t730) + qJ(2) * t746;];
tauB  = t1;
