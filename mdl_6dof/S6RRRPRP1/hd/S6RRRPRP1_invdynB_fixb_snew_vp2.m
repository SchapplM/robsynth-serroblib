% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:25:24
% EndTime: 2019-05-07 07:25:40
% DurationCPUTime: 13.28s
% Computational Cost: add. (196761->365), mult. (441950->448), div. (0->0), fcn. (323832->10), ass. (0->140)
t766 = Ifges(6,1) + Ifges(7,1);
t760 = Ifges(6,4) + Ifges(7,4);
t759 = Ifges(6,5) + Ifges(7,5);
t765 = Ifges(6,2) + Ifges(7,2);
t764 = Ifges(6,6) + Ifges(7,6);
t763 = Ifges(6,3) + Ifges(7,3);
t727 = sin(qJ(2));
t731 = cos(qJ(2));
t748 = qJD(1) * qJD(2);
t711 = qJDD(1) * t727 + t731 * t748;
t728 = sin(qJ(1));
t732 = cos(qJ(1));
t717 = -g(1) * t732 - g(2) * t728;
t733 = qJD(1) ^ 2;
t706 = -pkin(1) * t733 + qJDD(1) * pkin(7) + t717;
t757 = t727 * t706;
t762 = pkin(2) * t733;
t671 = qJDD(2) * pkin(2) - t711 * pkin(8) - t757 + (pkin(8) * t748 + t727 * t762 - g(3)) * t731;
t694 = -g(3) * t727 + t706 * t731;
t712 = qJDD(1) * t731 - t727 * t748;
t751 = qJD(1) * t727;
t715 = qJD(2) * pkin(2) - pkin(8) * t751;
t722 = t731 ^ 2;
t672 = pkin(8) * t712 - qJD(2) * t715 - t722 * t762 + t694;
t726 = sin(qJ(3));
t730 = cos(qJ(3));
t646 = t671 * t730 - t726 * t672;
t703 = (-t726 * t727 + t730 * t731) * qJD(1);
t678 = qJD(3) * t703 + t711 * t730 + t712 * t726;
t704 = (t726 * t731 + t727 * t730) * qJD(1);
t720 = qJDD(2) + qJDD(3);
t721 = qJD(2) + qJD(3);
t624 = (t703 * t721 - t678) * qJ(4) + (t703 * t704 + t720) * pkin(3) + t646;
t647 = t671 * t726 + t672 * t730;
t677 = -qJD(3) * t704 - t711 * t726 + t712 * t730;
t696 = pkin(3) * t721 - qJ(4) * t704;
t699 = t703 ^ 2;
t627 = -pkin(3) * t699 + qJ(4) * t677 - t696 * t721 + t647;
t723 = sin(pkin(10));
t724 = cos(pkin(10));
t691 = t703 * t723 + t704 * t724;
t618 = -0.2e1 * qJD(4) * t691 + t624 * t724 - t627 * t723;
t761 = -mrSges(6,2) - mrSges(7,2);
t690 = t703 * t724 - t704 * t723;
t619 = 0.2e1 * qJD(4) * t690 + t624 * t723 + t627 * t724;
t652 = t677 * t724 - t678 * t723;
t666 = -mrSges(5,1) * t690 + mrSges(5,2) * t691;
t681 = mrSges(5,1) * t721 - mrSges(5,3) * t691;
t667 = -pkin(4) * t690 - pkin(9) * t691;
t719 = t721 ^ 2;
t617 = -pkin(4) * t719 + pkin(9) * t720 + t667 * t690 + t619;
t716 = t728 * g(1) - g(2) * t732;
t738 = -qJDD(1) * pkin(1) - t716;
t679 = -t712 * pkin(2) + t715 * t751 + (-pkin(8) * t722 - pkin(7)) * t733 + t738;
t634 = -t677 * pkin(3) - t699 * qJ(4) + t696 * t704 + qJDD(4) + t679;
t653 = t677 * t723 + t678 * t724;
t622 = (-t690 * t721 - t653) * pkin(9) + (t691 * t721 - t652) * pkin(4) + t634;
t725 = sin(qJ(5));
t729 = cos(qJ(5));
t612 = -t725 * t617 + t622 * t729;
t675 = -t691 * t725 + t721 * t729;
t632 = qJD(5) * t675 + t653 * t729 + t720 * t725;
t651 = qJDD(5) - t652;
t676 = t691 * t729 + t721 * t725;
t654 = -mrSges(7,1) * t675 + mrSges(7,2) * t676;
t655 = -mrSges(6,1) * t675 + mrSges(6,2) * t676;
t684 = qJD(5) - t690;
t657 = -mrSges(6,2) * t684 + mrSges(6,3) * t675;
t609 = -0.2e1 * qJD(6) * t676 + (t675 * t684 - t632) * qJ(6) + (t675 * t676 + t651) * pkin(5) + t612;
t656 = -mrSges(7,2) * t684 + mrSges(7,3) * t675;
t747 = m(7) * t609 + mrSges(7,1) * t651 + t656 * t684;
t601 = m(6) * t612 + t651 * mrSges(6,1) + t684 * t657 + (-t654 - t655) * t676 + (-mrSges(6,3) - mrSges(7,3)) * t632 + t747;
t613 = t617 * t729 + t622 * t725;
t631 = -qJD(5) * t676 - t653 * t725 + t720 * t729;
t658 = pkin(5) * t684 - qJ(6) * t676;
t673 = t675 ^ 2;
t611 = -pkin(5) * t673 + qJ(6) * t631 + 0.2e1 * qJD(6) * t675 - t658 * t684 + t613;
t746 = m(7) * t611 + mrSges(7,3) * t631 + t654 * t675;
t659 = mrSges(7,1) * t684 - mrSges(7,3) * t676;
t752 = -mrSges(6,1) * t684 + mrSges(6,3) * t676 - t659;
t604 = m(6) * t613 + t631 * mrSges(6,3) + t651 * t761 + t675 * t655 + t684 * t752 + t746;
t741 = -t601 * t725 + t604 * t729;
t597 = m(5) * t619 - mrSges(5,2) * t720 + mrSges(5,3) * t652 + t666 * t690 - t681 * t721 + t741;
t680 = -mrSges(5,2) * t721 + mrSges(5,3) * t690;
t616 = -pkin(4) * t720 - pkin(9) * t719 + t667 * t691 - t618;
t614 = -pkin(5) * t631 - qJ(6) * t673 + t658 * t676 + qJDD(6) + t616;
t739 = m(7) * t614 - mrSges(7,1) * t631 - t656 * t675;
t735 = -m(6) * t616 + mrSges(6,1) * t631 + t632 * t761 + t657 * t675 + t676 * t752 - t739;
t606 = m(5) * t618 + t720 * mrSges(5,1) - t653 * mrSges(5,3) - t691 * t666 + t721 * t680 + t735;
t591 = t597 * t723 + t606 * t724;
t692 = -mrSges(4,1) * t703 + mrSges(4,2) * t704;
t695 = -mrSges(4,2) * t721 + mrSges(4,3) * t703;
t589 = m(4) * t646 + mrSges(4,1) * t720 - mrSges(4,3) * t678 - t692 * t704 + t695 * t721 + t591;
t697 = mrSges(4,1) * t721 - mrSges(4,3) * t704;
t742 = t597 * t724 - t606 * t723;
t590 = m(4) * t647 - mrSges(4,2) * t720 + mrSges(4,3) * t677 + t692 * t703 - t697 * t721 + t742;
t584 = t589 * t730 + t590 * t726;
t693 = -t731 * g(3) - t757;
t710 = (-mrSges(3,1) * t731 + mrSges(3,2) * t727) * qJD(1);
t750 = qJD(1) * t731;
t714 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t750;
t582 = m(3) * t693 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t711 + qJD(2) * t714 - t710 * t751 + t584;
t713 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t751;
t743 = -t589 * t726 + t590 * t730;
t583 = m(3) * t694 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t712 - qJD(2) * t713 + t710 * t750 + t743;
t744 = -t582 * t727 + t583 * t731;
t576 = m(2) * t717 - mrSges(2,1) * t733 - qJDD(1) * mrSges(2,2) + t744;
t705 = -t733 * pkin(7) + t738;
t599 = t601 * t729 + t604 * t725;
t737 = m(5) * t634 - mrSges(5,1) * t652 + mrSges(5,2) * t653 - t680 * t690 + t681 * t691 + t599;
t736 = m(4) * t679 - mrSges(4,1) * t677 + mrSges(4,2) * t678 - t695 * t703 + t697 * t704 + t737;
t734 = -m(3) * t705 + mrSges(3,1) * t712 - mrSges(3,2) * t711 - t713 * t751 + t714 * t750 - t736;
t594 = m(2) * t716 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t733 + t734;
t756 = t576 * t728 + t594 * t732;
t577 = t582 * t731 + t583 * t727;
t755 = t675 * t764 + t676 * t759 + t684 * t763;
t754 = -t675 * t765 - t676 * t760 - t684 * t764;
t753 = t675 * t760 + t676 * t766 + t684 * t759;
t745 = t576 * t732 - t594 * t728;
t702 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t727 + Ifges(3,4) * t731) * qJD(1);
t701 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t727 + Ifges(3,2) * t731) * qJD(1);
t700 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t727 + Ifges(3,6) * t731) * qJD(1);
t687 = Ifges(4,1) * t704 + Ifges(4,4) * t703 + Ifges(4,5) * t721;
t686 = Ifges(4,4) * t704 + Ifges(4,2) * t703 + Ifges(4,6) * t721;
t685 = Ifges(4,5) * t704 + Ifges(4,6) * t703 + Ifges(4,3) * t721;
t663 = Ifges(5,1) * t691 + Ifges(5,4) * t690 + Ifges(5,5) * t721;
t662 = Ifges(5,4) * t691 + Ifges(5,2) * t690 + Ifges(5,6) * t721;
t661 = Ifges(5,5) * t691 + Ifges(5,6) * t690 + Ifges(5,3) * t721;
t607 = -t632 * mrSges(7,3) - t676 * t654 + t747;
t598 = mrSges(6,2) * t616 + mrSges(7,2) * t614 - mrSges(6,3) * t612 - mrSges(7,3) * t609 - qJ(6) * t607 + t631 * t760 + t632 * t766 + t651 * t759 + t675 * t755 + t684 * t754;
t592 = -mrSges(6,1) * t616 + mrSges(6,3) * t613 - mrSges(7,1) * t614 + mrSges(7,3) * t611 - pkin(5) * t739 + qJ(6) * t746 + (-qJ(6) * t659 + t753) * t684 + (-pkin(5) * t659 - t755) * t676 + (-mrSges(7,2) * qJ(6) + t764) * t651 + (-mrSges(7,2) * pkin(5) + t760) * t632 + t765 * t631;
t585 = -mrSges(5,1) * t634 - mrSges(6,1) * t612 - mrSges(7,1) * t609 + mrSges(6,2) * t613 + mrSges(7,2) * t611 + mrSges(5,3) * t619 + Ifges(5,4) * t653 + Ifges(5,2) * t652 + Ifges(5,6) * t720 - pkin(4) * t599 - pkin(5) * t607 - t691 * t661 + t721 * t663 + t754 * t676 + t753 * t675 - t763 * t651 - t759 * t632 - t764 * t631;
t578 = mrSges(5,2) * t634 - mrSges(5,3) * t618 + Ifges(5,1) * t653 + Ifges(5,4) * t652 + Ifges(5,5) * t720 - pkin(9) * t599 - t592 * t725 + t598 * t729 + t661 * t690 - t662 * t721;
t573 = mrSges(4,2) * t679 - mrSges(4,3) * t646 + Ifges(4,1) * t678 + Ifges(4,4) * t677 + Ifges(4,5) * t720 - qJ(4) * t591 + t578 * t724 - t585 * t723 + t685 * t703 - t686 * t721;
t572 = -mrSges(4,1) * t679 + mrSges(4,3) * t647 + Ifges(4,4) * t678 + Ifges(4,2) * t677 + Ifges(4,6) * t720 - pkin(3) * t737 + qJ(4) * t742 + t723 * t578 + t724 * t585 - t704 * t685 + t721 * t687;
t571 = -pkin(9) * t741 - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) - pkin(1) * t577 - pkin(4) * t735 - pkin(2) * t584 + (-Ifges(4,3) - Ifges(5,3)) * t720 + (-t701 * t727 + t702 * t731) * qJD(1) + t733 * Ifges(2,5) - t729 * t592 - t725 * t598 + mrSges(2,3) * t717 - t704 * t686 - Ifges(3,5) * t711 - Ifges(3,6) * t712 - mrSges(3,1) * t693 + mrSges(3,2) * t694 + t703 * t687 + Ifges(2,6) * qJDD(1) + t690 * t663 - t691 * t662 - Ifges(4,6) * t677 - Ifges(4,5) * t678 - Ifges(5,6) * t652 - Ifges(5,5) * t653 - mrSges(4,1) * t646 + mrSges(4,2) * t647 - mrSges(5,1) * t618 + mrSges(5,2) * t619 - pkin(3) * t591;
t570 = mrSges(3,2) * t705 - mrSges(3,3) * t693 + Ifges(3,1) * t711 + Ifges(3,4) * t712 + Ifges(3,5) * qJDD(2) - pkin(8) * t584 - qJD(2) * t701 - t572 * t726 + t573 * t730 + t700 * t750;
t569 = -mrSges(3,1) * t705 + mrSges(3,3) * t694 + Ifges(3,4) * t711 + Ifges(3,2) * t712 + Ifges(3,6) * qJDD(2) - pkin(2) * t736 + pkin(8) * t743 + qJD(2) * t702 + t730 * t572 + t726 * t573 - t700 * t751;
t568 = -mrSges(2,2) * g(3) - mrSges(2,3) * t716 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t733 - pkin(7) * t577 - t569 * t727 + t570 * t731;
t1 = [-m(1) * g(1) + t745; -m(1) * g(2) + t756; (-m(1) - m(2)) * g(3) + t577; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t756 + t568 * t732 - t571 * t728; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t745 + t728 * t568 + t732 * t571; -mrSges(1,1) * g(2) + mrSges(2,1) * t716 + mrSges(1,2) * g(1) - mrSges(2,2) * t717 + Ifges(2,3) * qJDD(1) + pkin(1) * t734 + pkin(7) * t744 + t731 * t569 + t727 * t570;];
tauB  = t1;
