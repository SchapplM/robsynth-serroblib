% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRP7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:33
% EndTime: 2019-12-05 16:54:39
% DurationCPUTime: 3.55s
% Computational Cost: add. (39689->245), mult. (75204->305), div. (0->0), fcn. (48922->10), ass. (0->110)
t699 = Ifges(5,1) + Ifges(6,1);
t693 = Ifges(5,4) + Ifges(6,4);
t692 = Ifges(5,5) + Ifges(6,5);
t698 = Ifges(5,2) + Ifges(6,2);
t691 = Ifges(5,6) + Ifges(6,6);
t697 = Ifges(5,3) + Ifges(6,3);
t653 = sin(pkin(9));
t655 = cos(pkin(9));
t645 = t653 * g(1) - t655 * g(2);
t646 = -t655 * g(1) - t653 * g(2);
t652 = -g(3) + qJDD(1);
t654 = sin(pkin(5));
t656 = cos(pkin(5));
t659 = sin(qJ(2));
t662 = cos(qJ(2));
t602 = -t659 * t646 + (t645 * t656 + t652 * t654) * t662;
t657 = sin(qJ(4));
t660 = cos(qJ(4));
t658 = sin(qJ(3));
t681 = qJD(2) * t658;
t639 = t660 * qJD(3) - t657 * t681;
t661 = cos(qJ(3));
t679 = qJD(2) * qJD(3);
t675 = t661 * t679;
t643 = t658 * qJDD(2) + t675;
t616 = t639 * qJD(4) + t657 * qJDD(3) + t660 * t643;
t640 = t657 * qJD(3) + t660 * t681;
t680 = t661 * qJD(2);
t650 = qJD(4) - t680;
t624 = t650 * mrSges(6,1) - t640 * mrSges(6,3);
t687 = t656 * t659;
t688 = t654 * t659;
t603 = t645 * t687 + t662 * t646 + t652 * t688;
t664 = qJD(2) ^ 2;
t600 = -t664 * pkin(2) + qJDD(2) * pkin(7) + t603;
t626 = -t654 * t645 + t656 * t652;
t595 = -t658 * t600 + t661 * t626;
t642 = (-pkin(3) * t661 - pkin(8) * t658) * qJD(2);
t663 = qJD(3) ^ 2;
t590 = -qJDD(3) * pkin(3) - t663 * pkin(8) + t642 * t681 - t595;
t615 = -t640 * qJD(4) + t660 * qJDD(3) - t657 * t643;
t623 = t650 * pkin(4) - t640 * qJ(5);
t635 = t639 ^ 2;
t588 = -t615 * pkin(4) - t635 * qJ(5) + t640 * t623 + qJDD(5) + t590;
t621 = -t650 * mrSges(6,2) + t639 * mrSges(6,3);
t671 = -m(6) * t588 + t615 * mrSges(6,1) + t639 * t621;
t581 = t616 * mrSges(6,2) + t640 * t624 - t671;
t596 = t661 * t600 + t658 * t626;
t591 = -t663 * pkin(3) + qJDD(3) * pkin(8) + t642 * t680 + t596;
t599 = -qJDD(2) * pkin(2) - t664 * pkin(7) - t602;
t676 = t658 * t679;
t644 = t661 * qJDD(2) - t676;
t594 = (-t643 - t675) * pkin(8) + (-t644 + t676) * pkin(3) + t599;
t587 = t660 * t591 + t657 * t594;
t585 = -t635 * pkin(4) + t615 * qJ(5) + 0.2e1 * qJD(5) * t639 - t650 * t623 + t587;
t636 = qJDD(4) - t644;
t618 = -t639 * mrSges(6,1) + t640 * mrSges(6,2);
t677 = m(6) * t585 + t615 * mrSges(6,3) + t639 * t618;
t683 = t693 * t639 + t699 * t640 + t692 * t650;
t685 = -t691 * t639 - t692 * t640 - t697 * t650;
t564 = -mrSges(5,1) * t590 + mrSges(5,3) * t587 - mrSges(6,1) * t588 + mrSges(6,3) * t585 - pkin(4) * t581 + qJ(5) * t677 + (-qJ(5) * t624 + t683) * t650 + t685 * t640 + (-qJ(5) * mrSges(6,2) + t691) * t636 + t693 * t616 + t698 * t615;
t586 = -t657 * t591 + t660 * t594;
t583 = -0.2e1 * qJD(5) * t640 + (t639 * t650 - t616) * qJ(5) + (t639 * t640 + t636) * pkin(4) + t586;
t678 = m(6) * t583 + t636 * mrSges(6,1) + t650 * t621;
t580 = -t616 * mrSges(6,3) - t640 * t618 + t678;
t684 = -t698 * t639 - t693 * t640 - t691 * t650;
t571 = mrSges(5,2) * t590 + mrSges(6,2) * t588 - mrSges(5,3) * t586 - mrSges(6,3) * t583 - qJ(5) * t580 + t693 * t615 + t699 * t616 + t692 * t636 - t685 * t639 + t684 * t650;
t619 = -t639 * mrSges(5,1) + t640 * mrSges(5,2);
t622 = -t650 * mrSges(5,2) + t639 * mrSges(5,3);
t574 = m(5) * t586 + t636 * mrSges(5,1) + t650 * t622 + (-t618 - t619) * t640 + (-mrSges(5,3) - mrSges(6,3)) * t616 + t678;
t682 = -t650 * mrSges(5,1) + t640 * mrSges(5,3) - t624;
t694 = -mrSges(5,2) - mrSges(6,2);
t576 = m(5) * t587 + t615 * mrSges(5,3) + t639 * t619 + t694 * t636 + t682 * t650 + t677;
t573 = -t657 * t574 + t660 * t576;
t579 = -m(5) * t590 + t615 * mrSges(5,1) + t694 * t616 + t639 * t622 + t682 * t640 + t671;
t630 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t658 + Ifges(4,2) * t661) * qJD(2);
t631 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t658 + Ifges(4,4) * t661) * qJD(2);
t696 = mrSges(4,1) * t595 - mrSges(4,2) * t596 + Ifges(4,5) * t643 + Ifges(4,6) * t644 + Ifges(4,3) * qJDD(3) + pkin(3) * t579 + pkin(8) * t573 + t660 * t564 + t657 * t571 + (t658 * t630 - t661 * t631) * qJD(2);
t695 = mrSges(5,1) * t586 + mrSges(6,1) * t583 - mrSges(5,2) * t587 - mrSges(6,2) * t585 + pkin(4) * t580 + t691 * t615 + t692 * t616 + t697 * t636 - t683 * t639 - t684 * t640;
t572 = t660 * t574 + t657 * t576;
t647 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t681;
t648 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t680;
t666 = -m(4) * t599 + t644 * mrSges(4,1) - t643 * mrSges(4,2) - t647 * t681 + t648 * t680 - t572;
t567 = m(3) * t602 + qJDD(2) * mrSges(3,1) - t664 * mrSges(3,2) + t666;
t689 = t567 * t662;
t641 = (-mrSges(4,1) * t661 + mrSges(4,2) * t658) * qJD(2);
t570 = m(4) * t596 - qJDD(3) * mrSges(4,2) + t644 * mrSges(4,3) - qJD(3) * t647 + t641 * t680 + t573;
t578 = m(4) * t595 + qJDD(3) * mrSges(4,1) - t643 * mrSges(4,3) + qJD(3) * t648 - t641 * t681 + t579;
t673 = t661 * t570 - t658 * t578;
t560 = m(3) * t603 - t664 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t673;
t563 = t658 * t570 + t661 * t578;
t562 = m(3) * t626 + t563;
t550 = t560 * t687 - t654 * t562 + t656 * t689;
t548 = m(2) * t645 + t550;
t555 = t662 * t560 - t659 * t567;
t554 = m(2) * t646 + t555;
t686 = t655 * t548 + t653 * t554;
t549 = t560 * t688 + t656 * t562 + t654 * t689;
t674 = -t653 * t548 + t655 * t554;
t672 = m(2) * t652 + t549;
t629 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t658 + Ifges(4,6) * t661) * qJD(2);
t551 = mrSges(4,2) * t599 - mrSges(4,3) * t595 + Ifges(4,1) * t643 + Ifges(4,4) * t644 + Ifges(4,5) * qJDD(3) - pkin(8) * t572 - qJD(3) * t630 - t657 * t564 + t660 * t571 + t629 * t680;
t556 = -mrSges(4,1) * t599 + mrSges(4,3) * t596 + Ifges(4,4) * t643 + Ifges(4,2) * t644 + Ifges(4,6) * qJDD(3) - pkin(3) * t572 + qJD(3) * t631 - t629 * t681 - t695;
t545 = mrSges(3,2) * t626 - mrSges(3,3) * t602 + Ifges(3,5) * qJDD(2) - t664 * Ifges(3,6) - pkin(7) * t563 + t661 * t551 - t658 * t556;
t546 = -mrSges(3,1) * t626 + mrSges(3,3) * t603 + t664 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t563 - t696;
t668 = pkin(6) * t555 + t545 * t659 + t546 * t662;
t544 = mrSges(3,1) * t602 - mrSges(3,2) * t603 + Ifges(3,3) * qJDD(2) + pkin(2) * t666 + pkin(7) * t673 + t658 * t551 + t661 * t556;
t543 = mrSges(2,2) * t652 - mrSges(2,3) * t645 + t662 * t545 - t659 * t546 + (-t549 * t654 - t550 * t656) * pkin(6);
t542 = -mrSges(2,1) * t652 + mrSges(2,3) * t646 - pkin(1) * t549 - t654 * t544 + t668 * t656;
t1 = [-m(1) * g(1) + t674; -m(1) * g(2) + t686; -m(1) * g(3) + t672; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t686 - t653 * t542 + t655 * t543; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t674 + t655 * t542 + t653 * t543; -mrSges(1,1) * g(2) + mrSges(2,1) * t645 + mrSges(1,2) * g(1) - mrSges(2,2) * t646 + pkin(1) * t550 + t656 * t544 + t668 * t654; t672; t544; t696; t695; t581;];
tauJB = t1;
