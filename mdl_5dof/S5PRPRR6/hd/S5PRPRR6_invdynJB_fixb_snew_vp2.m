% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRR6
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:58
% EndTime: 2019-12-05 15:57:05
% DurationCPUTime: 5.99s
% Computational Cost: add. (69901->246), mult. (151891->319), div. (0->0), fcn. (110347->12), ass. (0->120)
t683 = qJD(2) ^ 2;
t670 = sin(pkin(10));
t673 = cos(pkin(10));
t677 = sin(qJ(4));
t680 = cos(qJ(4));
t691 = t670 * t677 - t673 * t680;
t651 = t691 * qJD(2);
t671 = sin(pkin(9));
t674 = cos(pkin(9));
t658 = t671 * g(1) - t674 * g(2);
t659 = -t674 * g(1) - t671 * g(2);
t669 = -g(3) + qJDD(1);
t672 = sin(pkin(5));
t675 = cos(pkin(5));
t678 = sin(qJ(2));
t681 = cos(qJ(2));
t634 = -t678 * t659 + (t658 * t675 + t669 * t672) * t681;
t692 = t670 * t680 + t673 * t677;
t652 = t692 * qJD(2);
t703 = t652 * qJD(4);
t641 = -t691 * qJDD(2) - t703;
t713 = pkin(3) * t673;
t712 = mrSges(4,2) * t670;
t688 = qJDD(3) - t634;
t627 = -qJDD(2) * pkin(2) - t683 * qJ(3) + t688;
t667 = t670 ^ 2;
t708 = t675 * t678;
t709 = t672 * t678;
t635 = t658 * t708 + t681 * t659 + t669 * t709;
t630 = -t683 * pkin(2) + qJDD(2) * qJ(3) + t635;
t649 = -t672 * t658 + t675 * t669;
t702 = qJD(2) * qJD(3);
t706 = t673 * t649 - 0.2e1 * t670 * t702;
t613 = (-pkin(7) * qJDD(2) + t683 * t713 - t630) * t670 + t706;
t616 = t670 * t649 + (t630 + 0.2e1 * t702) * t673;
t701 = qJDD(2) * t673;
t668 = t673 ^ 2;
t710 = t668 * t683;
t614 = -pkin(3) * t710 + pkin(7) * t701 + t616;
t609 = t677 * t613 + t680 * t614;
t640 = t651 * pkin(4) - t652 * pkin(8);
t682 = qJD(4) ^ 2;
t607 = -t682 * pkin(4) + qJDD(4) * pkin(8) - t651 * t640 + t609;
t624 = (-pkin(2) - t713) * qJDD(2) + (-qJ(3) + (-t667 - t668) * pkin(7)) * t683 + t688;
t704 = t651 * qJD(4);
t642 = t692 * qJDD(2) - t704;
t610 = (-t642 + t704) * pkin(8) + (-t641 + t703) * pkin(4) + t624;
t676 = sin(qJ(5));
t679 = cos(qJ(5));
t604 = -t676 * t607 + t679 * t610;
t643 = t679 * qJD(4) - t676 * t652;
t623 = t643 * qJD(5) + t676 * qJDD(4) + t679 * t642;
t644 = t676 * qJD(4) + t679 * t652;
t625 = -t643 * mrSges(6,1) + t644 * mrSges(6,2);
t650 = qJD(5) + t651;
t628 = -t650 * mrSges(6,2) + t643 * mrSges(6,3);
t639 = qJDD(5) - t641;
t601 = m(6) * t604 + t639 * mrSges(6,1) - t623 * mrSges(6,3) - t644 * t625 + t650 * t628;
t605 = t679 * t607 + t676 * t610;
t622 = -t644 * qJD(5) + t679 * qJDD(4) - t676 * t642;
t629 = t650 * mrSges(6,1) - t644 * mrSges(6,3);
t602 = m(6) * t605 - t639 * mrSges(6,2) + t622 * mrSges(6,3) + t643 * t625 - t650 * t629;
t592 = t679 * t601 + t676 * t602;
t647 = -qJD(4) * mrSges(5,2) - t651 * mrSges(5,3);
t648 = qJD(4) * mrSges(5,1) - t652 * mrSges(5,3);
t687 = m(5) * t624 - t641 * mrSges(5,1) + t642 * mrSges(5,2) + t651 * t647 + t652 * t648 + t592;
t686 = -m(4) * t627 + mrSges(4,1) * t701 - t687 + (t667 * t683 + t710) * mrSges(4,3);
t587 = t686 + (mrSges(3,1) - t712) * qJDD(2) - t683 * mrSges(3,2) + m(3) * t634;
t711 = t587 * t681;
t593 = -t676 * t601 + t679 * t602;
t637 = t651 * mrSges(5,1) + t652 * mrSges(5,2);
t590 = m(5) * t609 - qJDD(4) * mrSges(5,2) + t641 * mrSges(5,3) - qJD(4) * t648 - t651 * t637 + t593;
t608 = t680 * t613 - t677 * t614;
t606 = -qJDD(4) * pkin(4) - t682 * pkin(8) + t652 * t640 - t608;
t603 = -m(6) * t606 + t622 * mrSges(6,1) - t623 * mrSges(6,2) + t643 * t628 - t644 * t629;
t597 = m(5) * t608 + qJDD(4) * mrSges(5,1) - t642 * mrSges(5,3) + qJD(4) * t647 - t652 * t637 + t603;
t584 = t677 * t590 + t680 * t597;
t615 = -t670 * t630 + t706;
t690 = mrSges(4,3) * qJDD(2) + t683 * (-mrSges(4,1) * t673 + t712);
t582 = m(4) * t615 - t690 * t670 + t584;
t698 = t680 * t590 - t677 * t597;
t583 = m(4) * t616 + t690 * t673 + t698;
t699 = -t670 * t582 + t673 * t583;
t573 = m(3) * t635 - t683 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t699;
t576 = t673 * t582 + t670 * t583;
t575 = m(3) * t649 + t576;
t563 = t573 * t708 - t672 * t575 + t675 * t711;
t561 = m(2) * t658 + t563;
t569 = t681 * t573 - t678 * t587;
t568 = m(2) * t659 + t569;
t707 = t674 * t561 + t671 * t568;
t694 = Ifges(4,5) * t670 + Ifges(4,6) * t673;
t705 = t683 * t694;
t562 = t573 * t709 + t675 * t575 + t672 * t711;
t700 = -t671 * t561 + t674 * t568;
t697 = m(2) * t669 + t562;
t696 = Ifges(4,1) * t670 + Ifges(4,4) * t673;
t695 = Ifges(4,4) * t670 + Ifges(4,2) * t673;
t617 = Ifges(6,5) * t644 + Ifges(6,6) * t643 + Ifges(6,3) * t650;
t619 = Ifges(6,1) * t644 + Ifges(6,4) * t643 + Ifges(6,5) * t650;
t594 = -mrSges(6,1) * t606 + mrSges(6,3) * t605 + Ifges(6,4) * t623 + Ifges(6,2) * t622 + Ifges(6,6) * t639 - t644 * t617 + t650 * t619;
t618 = Ifges(6,4) * t644 + Ifges(6,2) * t643 + Ifges(6,6) * t650;
t595 = mrSges(6,2) * t606 - mrSges(6,3) * t604 + Ifges(6,1) * t623 + Ifges(6,4) * t622 + Ifges(6,5) * t639 + t643 * t617 - t650 * t618;
t631 = Ifges(5,5) * t652 - Ifges(5,6) * t651 + Ifges(5,3) * qJD(4);
t632 = Ifges(5,4) * t652 - Ifges(5,2) * t651 + Ifges(5,6) * qJD(4);
t577 = mrSges(5,2) * t624 - mrSges(5,3) * t608 + Ifges(5,1) * t642 + Ifges(5,4) * t641 + Ifges(5,5) * qJDD(4) - pkin(8) * t592 - qJD(4) * t632 - t676 * t594 + t679 * t595 - t651 * t631;
t633 = Ifges(5,1) * t652 - Ifges(5,4) * t651 + Ifges(5,5) * qJD(4);
t685 = mrSges(6,1) * t604 - mrSges(6,2) * t605 + Ifges(6,5) * t623 + Ifges(6,6) * t622 + Ifges(6,3) * t639 + t644 * t618 - t643 * t619;
t578 = -mrSges(5,1) * t624 + mrSges(5,3) * t609 + Ifges(5,4) * t642 + Ifges(5,2) * t641 + Ifges(5,6) * qJDD(4) - pkin(4) * t592 + qJD(4) * t633 - t652 * t631 - t685;
t564 = -mrSges(4,1) * t627 + mrSges(4,3) * t616 - pkin(3) * t687 + pkin(7) * t698 + t695 * qJDD(2) + t677 * t577 + t680 * t578 - t670 * t705;
t565 = mrSges(4,2) * t627 - mrSges(4,3) * t615 - pkin(7) * t584 + t696 * qJDD(2) + t680 * t577 - t677 * t578 + t673 * t705;
t558 = mrSges(3,2) * t649 - mrSges(3,3) * t634 + Ifges(3,5) * qJDD(2) - t683 * Ifges(3,6) - qJ(3) * t576 - t670 * t564 + t673 * t565;
t684 = mrSges(5,1) * t608 - mrSges(5,2) * t609 + Ifges(5,5) * t642 + Ifges(5,6) * t641 + Ifges(5,3) * qJDD(4) + pkin(4) * t603 + pkin(8) * t593 + t679 * t594 + t676 * t595 + t652 * t632 + t651 * t633;
t559 = -t684 + (Ifges(3,6) - t694) * qJDD(2) - mrSges(3,1) * t649 + mrSges(3,3) * t635 - mrSges(4,1) * t615 + mrSges(4,2) * t616 - pkin(3) * t584 - pkin(2) * t576 + (-t670 * t695 + t673 * t696 + Ifges(3,5)) * t683;
t689 = pkin(6) * t569 + t558 * t678 + t559 * t681;
t591 = qJDD(2) * t712 - t686;
t557 = mrSges(3,1) * t634 - mrSges(3,2) * t635 + Ifges(3,3) * qJDD(2) - pkin(2) * t591 + qJ(3) * t699 + t673 * t564 + t670 * t565;
t556 = mrSges(2,2) * t669 - mrSges(2,3) * t658 + t681 * t558 - t678 * t559 + (-t562 * t672 - t563 * t675) * pkin(6);
t555 = -mrSges(2,1) * t669 + mrSges(2,3) * t659 - pkin(1) * t562 - t672 * t557 + t689 * t675;
t1 = [-m(1) * g(1) + t700; -m(1) * g(2) + t707; -m(1) * g(3) + t697; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t707 - t671 * t555 + t674 * t556; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t700 + t674 * t555 + t671 * t556; -mrSges(1,1) * g(2) + mrSges(2,1) * t658 + mrSges(1,2) * g(1) - mrSges(2,2) * t659 + pkin(1) * t563 + t675 * t557 + t689 * t672; t697; t557; t591; t684; t685;];
tauJB = t1;
