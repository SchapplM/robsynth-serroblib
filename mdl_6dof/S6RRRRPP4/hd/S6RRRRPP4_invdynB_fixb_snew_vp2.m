% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:15:58
% EndTime: 2019-05-07 18:16:18
% DurationCPUTime: 15.05s
% Computational Cost: add. (233415->362), mult. (480891->442), div. (0->0), fcn. (346743->10), ass. (0->140)
t758 = Ifges(6,1) + Ifges(7,1);
t753 = Ifges(6,4) - Ifges(7,5);
t752 = Ifges(7,4) + Ifges(6,5);
t757 = Ifges(6,2) + Ifges(7,3);
t756 = -Ifges(7,2) - Ifges(6,3);
t751 = Ifges(6,6) - Ifges(7,6);
t755 = -2 * qJD(5);
t754 = -mrSges(6,3) - mrSges(7,2);
t750 = cos(pkin(10));
t722 = sin(qJ(1));
t726 = cos(qJ(1));
t711 = -g(1) * t726 - g(2) * t722;
t728 = qJD(1) ^ 2;
t696 = -pkin(1) * t728 + qJDD(1) * pkin(7) + t711;
t721 = sin(qJ(2));
t725 = cos(qJ(2));
t685 = -g(3) * t721 + t725 * t696;
t704 = (-mrSges(3,1) * t725 + mrSges(3,2) * t721) * qJD(1);
t742 = qJD(1) * qJD(2);
t715 = t721 * t742;
t707 = qJDD(1) * t725 - t715;
t744 = qJD(1) * t721;
t708 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t744;
t710 = t722 * g(1) - t726 * g(2);
t695 = -qJDD(1) * pkin(1) - t728 * pkin(7) - t710;
t740 = t725 * t742;
t706 = qJDD(1) * t721 + t740;
t662 = (-t706 - t740) * pkin(8) + (-t707 + t715) * pkin(2) + t695;
t705 = (-pkin(2) * t725 - pkin(8) * t721) * qJD(1);
t727 = qJD(2) ^ 2;
t743 = qJD(1) * t725;
t665 = -pkin(2) * t727 + qJDD(2) * pkin(8) + t705 * t743 + t685;
t720 = sin(qJ(3));
t724 = cos(qJ(3));
t643 = t724 * t662 - t720 * t665;
t702 = qJD(2) * t724 - t720 * t744;
t677 = qJD(3) * t702 + qJDD(2) * t720 + t706 * t724;
t701 = qJDD(3) - t707;
t703 = qJD(2) * t720 + t724 * t744;
t714 = qJD(3) - t743;
t622 = (t702 * t714 - t677) * pkin(9) + (t702 * t703 + t701) * pkin(3) + t643;
t644 = t720 * t662 + t724 * t665;
t676 = -qJD(3) * t703 + qJDD(2) * t724 - t706 * t720;
t686 = pkin(3) * t714 - pkin(9) * t703;
t700 = t702 ^ 2;
t624 = -pkin(3) * t700 + pkin(9) * t676 - t686 * t714 + t644;
t719 = sin(qJ(4));
t723 = cos(qJ(4));
t610 = t723 * t622 - t719 * t624;
t679 = t702 * t723 - t703 * t719;
t642 = qJD(4) * t679 + t676 * t719 + t677 * t723;
t680 = t702 * t719 + t703 * t723;
t697 = qJDD(4) + t701;
t713 = qJD(4) + t714;
t607 = (t679 * t713 - t642) * qJ(5) + (t679 * t680 + t697) * pkin(4) + t610;
t611 = t719 * t622 + t723 * t624;
t641 = -qJD(4) * t680 + t676 * t723 - t677 * t719;
t667 = pkin(4) * t713 - qJ(5) * t680;
t678 = t679 ^ 2;
t609 = -pkin(4) * t678 + qJ(5) * t641 - t667 * t713 + t611;
t718 = sin(pkin(10));
t657 = -t750 * t679 + t680 * t718;
t603 = t718 * t607 + t750 * t609 + t657 * t755;
t618 = -t750 * t641 + t642 * t718;
t658 = t718 * t679 + t750 * t680;
t649 = mrSges(6,1) * t713 - mrSges(6,3) * t658;
t635 = pkin(5) * t657 - qJ(6) * t658;
t712 = t713 ^ 2;
t600 = -pkin(5) * t712 + qJ(6) * t697 + 0.2e1 * qJD(6) * t713 - t635 * t657 + t603;
t650 = -mrSges(7,1) * t713 + mrSges(7,2) * t658;
t741 = m(7) * t600 + t697 * mrSges(7,3) + t713 * t650;
t636 = mrSges(7,1) * t657 - mrSges(7,3) * t658;
t745 = -mrSges(6,1) * t657 - mrSges(6,2) * t658 - t636;
t593 = m(6) * t603 - t697 * mrSges(6,2) + t754 * t618 - t713 * t649 + t745 * t657 + t741;
t733 = t750 * t607 - t718 * t609;
t602 = t658 * t755 + t733;
t619 = t718 * t641 + t750 * t642;
t648 = -mrSges(6,2) * t713 - mrSges(6,3) * t657;
t601 = -t697 * pkin(5) - t712 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t635) * t658 - t733;
t647 = -mrSges(7,2) * t657 + mrSges(7,3) * t713;
t734 = -m(7) * t601 + t697 * mrSges(7,1) + t713 * t647;
t595 = m(6) * t602 + t697 * mrSges(6,1) + t754 * t619 + t713 * t648 + t745 * t658 + t734;
t588 = t718 * t593 + t750 * t595;
t659 = -mrSges(5,1) * t679 + mrSges(5,2) * t680;
t666 = -mrSges(5,2) * t713 + mrSges(5,3) * t679;
t586 = m(5) * t610 + mrSges(5,1) * t697 - mrSges(5,3) * t642 - t659 * t680 + t666 * t713 + t588;
t668 = mrSges(5,1) * t713 - mrSges(5,3) * t680;
t735 = t750 * t593 - t595 * t718;
t587 = m(5) * t611 - mrSges(5,2) * t697 + mrSges(5,3) * t641 + t659 * t679 - t668 * t713 + t735;
t582 = t723 * t586 + t719 * t587;
t681 = -mrSges(4,1) * t702 + mrSges(4,2) * t703;
t682 = -mrSges(4,2) * t714 + mrSges(4,3) * t702;
t580 = m(4) * t643 + mrSges(4,1) * t701 - mrSges(4,3) * t677 - t681 * t703 + t682 * t714 + t582;
t683 = mrSges(4,1) * t714 - mrSges(4,3) * t703;
t736 = -t586 * t719 + t723 * t587;
t581 = m(4) * t644 - mrSges(4,2) * t701 + mrSges(4,3) * t676 + t681 * t702 - t683 * t714 + t736;
t737 = -t580 * t720 + t724 * t581;
t575 = m(3) * t685 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t707 - qJD(2) * t708 + t704 * t743 + t737;
t684 = -t725 * g(3) - t721 * t696;
t709 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t743;
t664 = -qJDD(2) * pkin(2) - pkin(8) * t727 + t705 * t744 - t684;
t638 = -pkin(3) * t676 - pkin(9) * t700 + t703 * t686 + t664;
t613 = -pkin(4) * t641 - qJ(5) * t678 + t680 * t667 + qJDD(5) + t638;
t605 = (t657 * t713 - t619) * qJ(6) + (t658 * t713 + t618) * pkin(5) - 0.2e1 * qJD(6) * t658 + t613;
t598 = m(7) * t605 + t618 * mrSges(7,1) - t619 * mrSges(7,3) + t657 * t647 - t658 * t650;
t732 = m(6) * t613 + t618 * mrSges(6,1) + t619 * mrSges(6,2) + t657 * t648 + t658 * t649 + t598;
t730 = m(5) * t638 - t641 * mrSges(5,1) + t642 * mrSges(5,2) - t679 * t666 + t680 * t668 + t732;
t729 = -m(4) * t664 + t676 * mrSges(4,1) - t677 * mrSges(4,2) + t702 * t682 - t703 * t683 - t730;
t597 = m(3) * t684 + qJDD(2) * mrSges(3,1) - t706 * mrSges(3,3) + qJD(2) * t709 - t704 * t744 + t729;
t738 = t725 * t575 - t597 * t721;
t569 = m(2) * t711 - mrSges(2,1) * t728 - qJDD(1) * mrSges(2,2) + t738;
t576 = t580 * t724 + t581 * t720;
t731 = -m(3) * t695 + t707 * mrSges(3,1) - mrSges(3,2) * t706 - t708 * t744 + t709 * t743 - t576;
t572 = m(2) * t710 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t728 + t731;
t749 = t722 * t569 + t726 * t572;
t570 = t721 * t575 + t725 * t597;
t748 = t757 * t657 - t753 * t658 - t751 * t713;
t747 = t751 * t657 - t752 * t658 + t756 * t713;
t746 = -t753 * t657 + t758 * t658 + t752 * t713;
t739 = t726 * t569 - t572 * t722;
t694 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t721 + Ifges(3,4) * t725) * qJD(1);
t693 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t721 + Ifges(3,2) * t725) * qJD(1);
t692 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t721 + Ifges(3,6) * t725) * qJD(1);
t671 = Ifges(4,1) * t703 + Ifges(4,4) * t702 + Ifges(4,5) * t714;
t670 = Ifges(4,4) * t703 + Ifges(4,2) * t702 + Ifges(4,6) * t714;
t669 = Ifges(4,5) * t703 + Ifges(4,6) * t702 + Ifges(4,3) * t714;
t653 = Ifges(5,1) * t680 + Ifges(5,4) * t679 + Ifges(5,5) * t713;
t652 = Ifges(5,4) * t680 + Ifges(5,2) * t679 + Ifges(5,6) * t713;
t651 = Ifges(5,5) * t680 + Ifges(5,6) * t679 + Ifges(5,3) * t713;
t590 = mrSges(6,2) * t613 + mrSges(7,2) * t601 - mrSges(6,3) * t602 - mrSges(7,3) * t605 - qJ(6) * t598 - t753 * t618 + t758 * t619 + t747 * t657 + t752 * t697 + t748 * t713;
t589 = -mrSges(6,1) * t613 - mrSges(7,1) * t605 + mrSges(7,2) * t600 + mrSges(6,3) * t603 - pkin(5) * t598 - t757 * t618 + t753 * t619 + t747 * t658 + t751 * t697 + t746 * t713;
t578 = mrSges(5,2) * t638 - mrSges(5,3) * t610 + Ifges(5,1) * t642 + Ifges(5,4) * t641 + Ifges(5,5) * t697 - qJ(5) * t588 - t718 * t589 + t750 * t590 + t679 * t651 - t713 * t652;
t577 = -mrSges(5,1) * t638 + mrSges(5,3) * t611 + Ifges(5,4) * t642 + Ifges(5,2) * t641 + Ifges(5,6) * t697 - pkin(4) * t732 + qJ(5) * t735 + t750 * t589 + t718 * t590 - t680 * t651 + t713 * t653;
t566 = mrSges(4,2) * t664 - mrSges(4,3) * t643 + Ifges(4,1) * t677 + Ifges(4,4) * t676 + Ifges(4,5) * t701 - pkin(9) * t582 - t577 * t719 + t578 * t723 + t669 * t702 - t670 * t714;
t565 = -mrSges(4,1) * t664 + mrSges(4,3) * t644 + Ifges(4,4) * t677 + Ifges(4,2) * t676 + Ifges(4,6) * t701 - pkin(3) * t730 + pkin(9) * t736 + t723 * t577 + t719 * t578 - t703 * t669 + t714 * t671;
t564 = Ifges(3,6) * qJDD(2) - pkin(2) * t576 + (-Ifges(5,3) + t756) * t697 + Ifges(3,4) * t706 + Ifges(3,2) * t707 + (mrSges(7,2) * qJ(6) + t751) * t618 + (mrSges(7,2) * pkin(5) - t752) * t619 - Ifges(4,6) * t676 - Ifges(4,5) * t677 + t679 * t653 - t680 * t652 - pkin(3) * t582 + mrSges(7,1) * t601 - mrSges(6,1) * t602 + (pkin(5) * t636 + t748) * t658 - mrSges(5,1) * t610 + mrSges(5,2) * t611 - mrSges(7,3) * t600 - qJ(6) * t741 - pkin(4) * t588 + mrSges(3,3) * t685 + qJD(2) * t694 - mrSges(3,1) * t695 - t692 * t744 + (qJ(6) * t636 - t746) * t657 + mrSges(6,2) * t603 - pkin(5) * t734 - Ifges(4,3) * t701 + t702 * t671 - t703 * t670 - Ifges(5,6) * t641 - Ifges(5,5) * t642 - mrSges(4,1) * t643 + mrSges(4,2) * t644;
t563 = mrSges(3,2) * t695 - mrSges(3,3) * t684 + Ifges(3,1) * t706 + Ifges(3,4) * t707 + Ifges(3,5) * qJDD(2) - pkin(8) * t576 - qJD(2) * t693 - t565 * t720 + t566 * t724 + t692 * t743;
t562 = Ifges(2,6) * qJDD(1) + t728 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t711 - Ifges(3,5) * t706 - Ifges(3,6) * t707 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t684 + mrSges(3,2) * t685 - t720 * t566 - t724 * t565 - pkin(2) * t729 - pkin(8) * t737 - pkin(1) * t570 + (-t693 * t721 + t694 * t725) * qJD(1);
t561 = -mrSges(2,2) * g(3) - mrSges(2,3) * t710 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t728 - pkin(7) * t570 + t563 * t725 - t564 * t721;
t1 = [-m(1) * g(1) + t739; -m(1) * g(2) + t749; (-m(1) - m(2)) * g(3) + t570; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t749 + t726 * t561 - t722 * t562; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t739 + t722 * t561 + t726 * t562; -mrSges(1,1) * g(2) + mrSges(2,1) * t710 + mrSges(1,2) * g(1) - mrSges(2,2) * t711 + Ifges(2,3) * qJDD(1) + pkin(1) * t731 + pkin(7) * t738 + t721 * t563 + t725 * t564;];
tauB  = t1;
