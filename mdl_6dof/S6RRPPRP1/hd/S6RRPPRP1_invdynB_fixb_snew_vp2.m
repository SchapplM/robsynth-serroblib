% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-05-06 09:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:03:05
% EndTime: 2019-05-06 09:03:16
% DurationCPUTime: 11.04s
% Computational Cost: add. (155797->362), mult. (367464->447), div. (0->0), fcn. (259272->10), ass. (0->138)
t756 = -2 * qJD(3);
t755 = Ifges(6,1) + Ifges(7,1);
t748 = Ifges(6,4) - Ifges(7,5);
t754 = -Ifges(6,5) - Ifges(7,4);
t753 = Ifges(6,2) + Ifges(7,3);
t746 = Ifges(6,6) - Ifges(7,6);
t752 = -Ifges(6,3) - Ifges(7,2);
t716 = sin(qJ(2));
t718 = cos(qJ(2));
t735 = qJD(1) * qJD(2);
t702 = t716 * qJDD(1) + t718 * t735;
t717 = sin(qJ(1));
t719 = cos(qJ(1));
t708 = -t719 * g(1) - t717 * g(2);
t721 = qJD(1) ^ 2;
t697 = -t721 * pkin(1) + qJDD(1) * pkin(7) + t708;
t744 = t716 * t697;
t750 = pkin(2) * t721;
t651 = qJDD(2) * pkin(2) - t702 * qJ(3) - t744 + (qJ(3) * t735 + t716 * t750 - g(3)) * t718;
t682 = -t716 * g(3) + t718 * t697;
t703 = t718 * qJDD(1) - t716 * t735;
t738 = qJD(1) * t716;
t704 = qJD(2) * pkin(2) - qJ(3) * t738;
t711 = t718 ^ 2;
t653 = t703 * qJ(3) - qJD(2) * t704 - t711 * t750 + t682;
t713 = sin(pkin(9));
t745 = cos(pkin(9));
t691 = (t713 * t718 + t745 * t716) * qJD(1);
t630 = t745 * t651 - t713 * t653 + t691 * t756;
t751 = cos(qJ(5));
t749 = -mrSges(6,3) - mrSges(7,2);
t737 = qJD(1) * t718;
t690 = t713 * t738 - t745 * t737;
t631 = t713 * t651 + t745 * t653 + t690 * t756;
t668 = t690 * mrSges(4,1) + t691 * mrSges(4,2);
t673 = t713 * t702 - t745 * t703;
t684 = qJD(2) * mrSges(4,1) - t691 * mrSges(4,3);
t667 = t690 * pkin(3) - t691 * qJ(4);
t720 = qJD(2) ^ 2;
t612 = -t720 * pkin(3) + qJDD(2) * qJ(4) - t690 * t667 + t631;
t707 = t717 * g(1) - t719 * g(2);
t726 = -qJDD(1) * pkin(1) - t707;
t655 = -t703 * pkin(2) + qJDD(3) + t704 * t738 + (-qJ(3) * t711 - pkin(7)) * t721 + t726;
t674 = t745 * t702 + t713 * t703;
t615 = (qJD(2) * t690 - t674) * qJ(4) + (qJD(2) * t691 + t673) * pkin(3) + t655;
t712 = sin(pkin(10));
t714 = cos(pkin(10));
t680 = t712 * qJD(2) + t714 * t691;
t607 = -0.2e1 * qJD(4) * t680 - t712 * t612 + t714 * t615;
t662 = t712 * qJDD(2) + t714 * t674;
t679 = t714 * qJD(2) - t712 * t691;
t604 = (t679 * t690 - t662) * pkin(8) + (t679 * t680 + t673) * pkin(4) + t607;
t608 = 0.2e1 * qJD(4) * t679 + t714 * t612 + t712 * t615;
t658 = t690 * pkin(4) - t680 * pkin(8);
t661 = t714 * qJDD(2) - t712 * t674;
t678 = t679 ^ 2;
t606 = -t678 * pkin(4) + t661 * pkin(8) - t690 * t658 + t608;
t715 = sin(qJ(5));
t600 = t715 * t604 + t751 * t606;
t648 = t715 * t679 + t751 * t680;
t619 = t648 * qJD(5) - t751 * t661 + t715 * t662;
t689 = qJD(5) + t690;
t639 = t689 * mrSges(6,1) - t648 * mrSges(6,3);
t647 = -t751 * t679 + t715 * t680;
t672 = qJDD(5) + t673;
t632 = t647 * pkin(5) - t648 * qJ(6);
t688 = t689 ^ 2;
t597 = -t688 * pkin(5) + t672 * qJ(6) + 0.2e1 * qJD(6) * t689 - t647 * t632 + t600;
t640 = -t689 * mrSges(7,1) + t648 * mrSges(7,2);
t734 = m(7) * t597 + t672 * mrSges(7,3) + t689 * t640;
t633 = t647 * mrSges(7,1) - t648 * mrSges(7,3);
t739 = -t647 * mrSges(6,1) - t648 * mrSges(6,2) - t633;
t590 = m(6) * t600 - t672 * mrSges(6,2) + t749 * t619 - t689 * t639 + t739 * t647 + t734;
t599 = t751 * t604 - t715 * t606;
t620 = -t647 * qJD(5) + t715 * t661 + t751 * t662;
t638 = -t689 * mrSges(6,2) - t647 * mrSges(6,3);
t598 = -t672 * pkin(5) - t688 * qJ(6) + t648 * t632 + qJDD(6) - t599;
t637 = -t647 * mrSges(7,2) + t689 * mrSges(7,3);
t728 = -m(7) * t598 + t672 * mrSges(7,1) + t689 * t637;
t592 = m(6) * t599 + t672 * mrSges(6,1) + t749 * t620 + t689 * t638 + t739 * t648 + t728;
t585 = t715 * t590 + t751 * t592;
t652 = -t679 * mrSges(5,1) + t680 * mrSges(5,2);
t656 = -t690 * mrSges(5,2) + t679 * mrSges(5,3);
t583 = m(5) * t607 + t673 * mrSges(5,1) - t662 * mrSges(5,3) - t680 * t652 + t690 * t656 + t585;
t657 = t690 * mrSges(5,1) - t680 * mrSges(5,3);
t729 = t751 * t590 - t715 * t592;
t584 = m(5) * t608 - t673 * mrSges(5,2) + t661 * mrSges(5,3) + t679 * t652 - t690 * t657 + t729;
t730 = -t712 * t583 + t714 * t584;
t578 = m(4) * t631 - qJDD(2) * mrSges(4,2) - t673 * mrSges(4,3) - qJD(2) * t684 - t690 * t668 + t730;
t683 = -qJD(2) * mrSges(4,2) - t690 * mrSges(4,3);
t611 = -qJDD(2) * pkin(3) - t720 * qJ(4) + t691 * t667 + qJDD(4) - t630;
t609 = -t661 * pkin(4) - t678 * pkin(8) + t680 * t658 + t611;
t602 = -0.2e1 * qJD(6) * t648 + (t647 * t689 - t620) * qJ(6) + (t648 * t689 + t619) * pkin(5) + t609;
t595 = m(7) * t602 + t619 * mrSges(7,1) - t620 * mrSges(7,3) + t647 * t637 - t648 * t640;
t724 = m(6) * t609 + t619 * mrSges(6,1) + t620 * mrSges(6,2) + t647 * t638 + t648 * t639 + t595;
t722 = -m(5) * t611 + t661 * mrSges(5,1) - t662 * mrSges(5,2) + t679 * t656 - t680 * t657 - t724;
t594 = m(4) * t630 + qJDD(2) * mrSges(4,1) - t674 * mrSges(4,3) + qJD(2) * t683 - t691 * t668 + t722;
t573 = t713 * t578 + t745 * t594;
t681 = -t718 * g(3) - t744;
t701 = (-mrSges(3,1) * t718 + mrSges(3,2) * t716) * qJD(1);
t706 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t737;
t571 = m(3) * t681 + qJDD(2) * mrSges(3,1) - t702 * mrSges(3,3) + qJD(2) * t706 - t701 * t738 + t573;
t705 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t738;
t731 = t745 * t578 - t713 * t594;
t572 = m(3) * t682 - qJDD(2) * mrSges(3,2) + t703 * mrSges(3,3) - qJD(2) * t705 + t701 * t737 + t731;
t732 = -t716 * t571 + t718 * t572;
t563 = m(2) * t708 - t721 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t732;
t696 = -t721 * pkin(7) + t726;
t579 = t714 * t583 + t712 * t584;
t725 = m(4) * t655 + t673 * mrSges(4,1) + t674 * mrSges(4,2) + t690 * t683 + t691 * t684 + t579;
t723 = -m(3) * t696 + t703 * mrSges(3,1) - t702 * mrSges(3,2) - t705 * t738 + t706 * t737 - t725;
t575 = m(2) * t707 + qJDD(1) * mrSges(2,1) - t721 * mrSges(2,2) + t723;
t743 = t717 * t563 + t719 * t575;
t564 = t718 * t571 + t716 * t572;
t742 = t753 * t647 - t748 * t648 - t746 * t689;
t741 = t746 * t647 + t754 * t648 + t752 * t689;
t740 = -t748 * t647 + t755 * t648 - t754 * t689;
t733 = t719 * t563 - t717 * t575;
t694 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t716 + Ifges(3,4) * t718) * qJD(1);
t693 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t716 + Ifges(3,2) * t718) * qJD(1);
t692 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t716 + Ifges(3,6) * t718) * qJD(1);
t665 = Ifges(4,1) * t691 - Ifges(4,4) * t690 + Ifges(4,5) * qJD(2);
t664 = Ifges(4,4) * t691 - Ifges(4,2) * t690 + Ifges(4,6) * qJD(2);
t663 = Ifges(4,5) * t691 - Ifges(4,6) * t690 + Ifges(4,3) * qJD(2);
t643 = Ifges(5,1) * t680 + Ifges(5,4) * t679 + Ifges(5,5) * t690;
t642 = Ifges(5,4) * t680 + Ifges(5,2) * t679 + Ifges(5,6) * t690;
t641 = Ifges(5,5) * t680 + Ifges(5,6) * t679 + Ifges(5,3) * t690;
t587 = mrSges(6,2) * t609 + mrSges(7,2) * t598 - mrSges(6,3) * t599 - mrSges(7,3) * t602 - qJ(6) * t595 - t748 * t619 + t755 * t620 + t741 * t647 - t672 * t754 + t742 * t689;
t586 = -mrSges(6,1) * t609 - mrSges(7,1) * t602 + mrSges(7,2) * t597 + mrSges(6,3) * t600 - pkin(5) * t595 - t753 * t619 + t748 * t620 + t741 * t648 + t746 * t672 + t740 * t689;
t567 = mrSges(5,2) * t611 - mrSges(5,3) * t607 + Ifges(5,1) * t662 + Ifges(5,4) * t661 + Ifges(5,5) * t673 - pkin(8) * t585 - t715 * t586 + t751 * t587 + t679 * t641 - t690 * t642;
t566 = -mrSges(5,1) * t611 + mrSges(5,3) * t608 + Ifges(5,4) * t662 + Ifges(5,2) * t661 + Ifges(5,6) * t673 - pkin(4) * t724 + pkin(8) * t729 + t751 * t586 + t715 * t587 - t680 * t641 + t690 * t643;
t565 = (qJ(6) * t633 - t740) * t647 + (pkin(5) * t633 + t742) * t648 - qJ(6) * t734 - pkin(5) * t728 + t752 * t672 + (-Ifges(5,3) - Ifges(4,2)) * t673 + Ifges(4,6) * qJDD(2) - t691 * t663 - t680 * t642 + (pkin(5) * mrSges(7,2) + t754) * t620 + Ifges(4,4) * t674 + t679 * t643 - Ifges(5,6) * t661 - Ifges(5,5) * t662 + qJD(2) * t665 - mrSges(4,1) * t655 + mrSges(4,3) * t631 - mrSges(5,1) * t607 + mrSges(5,2) * t608 - mrSges(6,1) * t599 + mrSges(6,2) * t600 - mrSges(7,3) * t597 + mrSges(7,1) * t598 + (qJ(6) * mrSges(7,2) + t746) * t619 - pkin(4) * t585 - pkin(3) * t579;
t560 = mrSges(4,2) * t655 - mrSges(4,3) * t630 + Ifges(4,1) * t674 - Ifges(4,4) * t673 + Ifges(4,5) * qJDD(2) - qJ(4) * t579 - qJD(2) * t664 - t712 * t566 + t714 * t567 - t690 * t663;
t559 = mrSges(3,2) * t696 - mrSges(3,3) * t681 + Ifges(3,1) * t702 + Ifges(3,4) * t703 + Ifges(3,5) * qJDD(2) - qJ(3) * t573 - qJD(2) * t693 + t745 * t560 - t713 * t565 + t692 * t737;
t558 = -pkin(1) * t564 + mrSges(2,3) * t708 - pkin(2) * t573 - Ifges(3,5) * t702 - Ifges(3,6) * t703 - mrSges(3,1) * t681 + mrSges(3,2) * t682 - pkin(3) * t722 - qJ(4) * t730 - t712 * t567 - t714 * t566 - Ifges(4,5) * t674 + Ifges(4,6) * t673 - mrSges(4,1) * t630 + mrSges(4,2) * t631 + mrSges(2,1) * g(3) + t721 * Ifges(2,5) - t691 * t664 - t690 * t665 + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t716 * t693 + t718 * t694) * qJD(1);
t557 = -mrSges(3,1) * t696 + mrSges(3,3) * t682 + Ifges(3,4) * t702 + Ifges(3,2) * t703 + Ifges(3,6) * qJDD(2) - pkin(2) * t725 + qJ(3) * t731 + qJD(2) * t694 + t713 * t560 + t745 * t565 - t692 * t738;
t556 = -mrSges(2,2) * g(3) - mrSges(2,3) * t707 + Ifges(2,5) * qJDD(1) - t721 * Ifges(2,6) - pkin(7) * t564 - t716 * t557 + t718 * t559;
t1 = [-m(1) * g(1) + t733; -m(1) * g(2) + t743; (-m(1) - m(2)) * g(3) + t564; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t743 + t719 * t556 - t717 * t558; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t733 + t717 * t556 + t719 * t558; -mrSges(1,1) * g(2) + mrSges(2,1) * t707 + mrSges(1,2) * g(1) - mrSges(2,2) * t708 + Ifges(2,3) * qJDD(1) + pkin(1) * t723 + pkin(7) * t732 + t718 * t557 + t716 * t559;];
tauB  = t1;
