% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR13_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:59
% EndTime: 2024-09-27 17:31:08
% DurationCPUTime: 3.83s
% Computational Cost: add. (177085->245), mult. (206851->324), div. (0->0), fcn. (131433->12), ass. (0->115)
t706 = -m(3) - m(4);
t670 = sin(pkin(5));
t705 = pkin(9) * t670;
t671 = cos(pkin(5));
t704 = t671 * g(3);
t667 = qJD(1) + qJD(2);
t662 = qJD(3) + t667;
t660 = t662 ^ 2;
t703 = t660 * t670 ^ 2;
t702 = t662 * t670;
t673 = sin(qJ(4));
t701 = t670 * t673;
t678 = cos(qJ(4));
t700 = t670 * t678;
t699 = t671 * t673;
t698 = t671 * t678;
t665 = qJDD(1) + qJDD(2);
t661 = qJDD(3) + t665;
t696 = qJD(4) * t662;
t641 = (t661 * t673 + t678 * t696) * t670;
t651 = t671 * t661 + qJDD(4);
t653 = t671 * t662 + qJD(4);
t676 = sin(qJ(1));
t681 = cos(qJ(1));
t654 = t676 * g(1) - t681 * g(2);
t648 = qJDD(1) * pkin(1) + t654;
t655 = -t681 * g(1) - t676 * g(2);
t682 = qJD(1) ^ 2;
t649 = -t682 * pkin(1) + t655;
t675 = sin(qJ(2));
t680 = cos(qJ(2));
t633 = t680 * t648 - t675 * t649;
t630 = t665 * pkin(2) + t633;
t634 = t675 * t648 + t680 * t649;
t664 = t667 ^ 2;
t631 = -t664 * pkin(2) + t634;
t674 = sin(qJ(3));
t679 = cos(qJ(3));
t615 = t679 * t630 - t674 * t631;
t607 = t661 * pkin(3) + t660 * t705 + t615;
t616 = t674 * t630 + t679 * t631;
t608 = -t660 * pkin(3) + t661 * t705 + t616;
t688 = t607 * t698 - t673 * t608;
t596 = t651 * pkin(4) - t641 * pkin(10) + (pkin(4) * t673 * t703 + (pkin(10) * t653 * t662 - g(3)) * t670) * t678 + t688;
t599 = -g(3) * t701 + t607 * t699 + t678 * t608;
t694 = t662 * t701;
t640 = t653 * pkin(4) - pkin(10) * t694;
t642 = (t661 * t678 - t673 * t696) * t670;
t695 = t678 ^ 2 * t703;
t597 = -pkin(4) * t695 + t642 * pkin(10) - t653 * t640 + t599;
t672 = sin(qJ(5));
t677 = cos(qJ(5));
t594 = t677 * t596 - t672 * t597;
t635 = (-t672 * t673 + t677 * t678) * t702;
t613 = t635 * qJD(5) + t677 * t641 + t672 * t642;
t636 = (t672 * t678 + t673 * t677) * t702;
t621 = -t635 * mrSges(6,1) + t636 * mrSges(6,2);
t650 = qJD(5) + t653;
t622 = -t650 * mrSges(6,2) + t635 * mrSges(6,3);
t647 = qJDD(5) + t651;
t589 = m(6) * t594 + t647 * mrSges(6,1) - t613 * mrSges(6,3) - t636 * t621 + t650 * t622;
t595 = t672 * t596 + t677 * t597;
t612 = -t636 * qJD(5) - t672 * t641 + t677 * t642;
t623 = t650 * mrSges(6,1) - t636 * mrSges(6,3);
t590 = m(6) * t595 - t647 * mrSges(6,2) + t612 * mrSges(6,3) + t635 * t621 - t650 * t623;
t583 = t677 * t589 + t672 * t590;
t598 = -g(3) * t700 + t688;
t693 = t662 * t700;
t638 = -t653 * mrSges(5,2) + mrSges(5,3) * t693;
t639 = (-mrSges(5,1) * t678 + mrSges(5,2) * t673) * t702;
t581 = m(5) * t598 + t651 * mrSges(5,1) - t641 * mrSges(5,3) + t653 * t638 - t639 * t694 + t583;
t637 = t653 * mrSges(5,1) - mrSges(5,3) * t694;
t689 = -t672 * t589 + t677 * t590;
t582 = m(5) * t599 - t651 * mrSges(5,2) + t642 * mrSges(5,3) - t653 * t637 + t639 * t693 + t689;
t602 = -t670 * t607 - t704;
t601 = -pkin(10) * t695 - t642 * pkin(4) - t704 + (t640 * t662 * t673 - t607) * t670;
t687 = m(6) * t601 - t612 * mrSges(6,1) + t613 * mrSges(6,2) - t635 * t622 + t636 * t623;
t592 = m(5) * t602 - t642 * mrSges(5,1) + t641 * mrSges(5,2) + (t637 * t673 - t638 * t678) * t702 + t687;
t566 = t581 * t698 + t582 * t699 - t670 * t592;
t563 = m(4) * t615 + t661 * mrSges(4,1) - t660 * mrSges(4,2) + t566;
t573 = -t673 * t581 + t678 * t582;
t571 = m(4) * t616 - t660 * mrSges(4,1) - t661 * mrSges(4,2) + t573;
t559 = t679 * t563 + t674 * t571;
t556 = m(3) * t633 + t665 * mrSges(3,1) - t664 * mrSges(3,2) + t559;
t690 = -t674 * t563 + t679 * t571;
t557 = m(3) * t634 - t664 * mrSges(3,1) - t665 * mrSges(3,2) + t690;
t552 = t680 * t556 + t675 * t557;
t549 = m(2) * t654 + qJDD(1) * mrSges(2,1) - t682 * mrSges(2,2) + t552;
t691 = -t675 * t556 + t680 * t557;
t550 = m(2) * t655 - t682 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t691;
t697 = t681 * t549 + t676 * t550;
t565 = t581 * t700 + t582 * t701 + t671 * t592;
t692 = -t676 * t549 + t681 * t550;
t617 = Ifges(6,5) * t636 + Ifges(6,6) * t635 + Ifges(6,3) * t650;
t619 = Ifges(6,1) * t636 + Ifges(6,4) * t635 + Ifges(6,5) * t650;
t584 = -mrSges(6,1) * t601 + mrSges(6,3) * t595 + Ifges(6,4) * t613 + Ifges(6,2) * t612 + Ifges(6,6) * t647 - t636 * t617 + t650 * t619;
t618 = Ifges(6,4) * t636 + Ifges(6,2) * t635 + Ifges(6,6) * t650;
t585 = mrSges(6,2) * t601 - mrSges(6,3) * t594 + Ifges(6,1) * t613 + Ifges(6,4) * t612 + Ifges(6,5) * t647 + t635 * t617 - t650 * t618;
t627 = Ifges(5,3) * t653 + (Ifges(5,5) * t673 + Ifges(5,6) * t678) * t702;
t629 = Ifges(5,5) * t653 + (Ifges(5,1) * t673 + Ifges(5,4) * t678) * t702;
t561 = -mrSges(5,1) * t602 + mrSges(5,3) * t599 + Ifges(5,4) * t641 + Ifges(5,2) * t642 + Ifges(5,6) * t651 - pkin(4) * t687 + pkin(10) * t689 + t677 * t584 + t672 * t585 - t627 * t694 + t653 * t629;
t628 = Ifges(5,6) * t653 + (Ifges(5,4) * t673 + Ifges(5,2) * t678) * t702;
t568 = mrSges(5,2) * t602 - mrSges(5,3) * t598 + Ifges(5,1) * t641 + Ifges(5,4) * t642 + Ifges(5,5) * t651 - pkin(10) * t583 - t672 * t584 + t677 * t585 + t627 * t693 - t653 * t628;
t685 = mrSges(6,1) * t594 - mrSges(6,2) * t595 + Ifges(6,5) * t613 + Ifges(6,6) * t612 + Ifges(6,3) * t647 + t636 * t618 - t635 * t619;
t575 = mrSges(5,1) * t598 - mrSges(5,2) * t599 + Ifges(5,5) * t641 + Ifges(5,6) * t642 + Ifges(5,3) * t651 + pkin(4) * t583 + (t628 * t673 - t629 * t678) * t702 + t685;
t686 = mrSges(4,1) * t615 - mrSges(4,2) * t616 + Ifges(4,3) * t661 + pkin(3) * t566 + t561 * t700 + t568 * t701 + t573 * t705 + t671 * t575;
t684 = mrSges(3,1) * t633 - mrSges(3,2) * t634 + Ifges(3,3) * t665 + pkin(2) * t559 + t686;
t683 = mrSges(2,1) * t654 - mrSges(2,2) * t655 + Ifges(2,3) * qJDD(1) + pkin(1) * t552 + t684;
t545 = -mrSges(4,2) * g(3) - mrSges(4,3) * t615 + Ifges(4,5) * t661 - t660 * Ifges(4,6) - t673 * t561 + t678 * t568 + (-t565 * t670 - t566 * t671) * pkin(9);
t544 = mrSges(4,1) * g(3) + mrSges(4,3) * t616 + t660 * Ifges(4,5) + Ifges(4,6) * t661 - pkin(3) * t565 - t670 * t575 + (pkin(9) * t573 + t561 * t678 + t568 * t673) * t671;
t543 = -mrSges(3,2) * g(3) - mrSges(3,3) * t633 + Ifges(3,5) * t665 - t664 * Ifges(3,6) - pkin(8) * t559 - t674 * t544 + t679 * t545;
t542 = Ifges(3,6) * t665 + t664 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t634 + t674 * t545 + t679 * t544 - pkin(2) * (-m(4) * g(3) + t565) + pkin(8) * t690;
t541 = -mrSges(2,2) * g(3) - mrSges(2,3) * t654 + Ifges(2,5) * qJDD(1) - t682 * Ifges(2,6) - pkin(7) * t552 - t675 * t542 + t680 * t543;
t540 = Ifges(2,6) * qJDD(1) + t682 * Ifges(2,5) + mrSges(2,3) * t655 + t675 * t543 + t680 * t542 - pkin(1) * t565 + pkin(7) * t691 + (-pkin(1) * t706 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t692; -m(1) * g(2) + t697; (-m(1) - m(2) + t706) * g(3) + t565; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t697 - t676 * t540 + t681 * t541; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t692 + t681 * t540 + t676 * t541; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t683; t683; t684; t686; t575; t685;];
tauJB = t1;
