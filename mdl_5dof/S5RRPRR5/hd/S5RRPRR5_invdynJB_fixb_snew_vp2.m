% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR5
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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:41
% EndTime: 2022-01-20 11:02:48
% DurationCPUTime: 6.60s
% Computational Cost: add. (94999->251), mult. (131436->313), div. (0->0), fcn. (88562->10), ass. (0->114)
t667 = qJD(1) + qJD(2);
t661 = t667 ^ 2;
t670 = cos(pkin(9));
t707 = pkin(3) * t670;
t669 = sin(pkin(9));
t706 = mrSges(4,2) * t669;
t691 = Ifges(4,5) * t669 + Ifges(4,6) * t670;
t705 = t691 * t661;
t665 = t670 ^ 2;
t704 = t661 * t665;
t663 = qJDD(1) + qJDD(2);
t703 = t663 * t670;
t674 = sin(qJ(1));
t678 = cos(qJ(1));
t653 = t674 * g(1) - g(2) * t678;
t648 = qJDD(1) * pkin(1) + t653;
t654 = -g(1) * t678 - g(2) * t674;
t679 = qJD(1) ^ 2;
t649 = -pkin(1) * t679 + t654;
t673 = sin(qJ(2));
t677 = cos(qJ(2));
t636 = t673 * t648 + t677 * t649;
t633 = -pkin(2) * t661 + qJ(3) * t663 + t636;
t701 = qJD(3) * t667;
t699 = -g(3) * t670 - 0.2e1 * t669 * t701;
t614 = (-pkin(7) * t663 + t661 * t707 - t633) * t669 + t699;
t618 = -g(3) * t669 + (t633 + 0.2e1 * t701) * t670;
t615 = -pkin(3) * t704 + pkin(7) * t703 + t618;
t672 = sin(qJ(4));
t676 = cos(qJ(4));
t596 = t676 * t614 - t615 * t672;
t687 = t669 * t676 + t670 * t672;
t686 = -t669 * t672 + t670 * t676;
t641 = t686 * t667;
t700 = qJD(4) * t641;
t632 = t687 * t663 + t700;
t642 = t687 * t667;
t592 = (-t632 + t700) * pkin(8) + (t641 * t642 + qJDD(4)) * pkin(4) + t596;
t597 = t672 * t614 + t676 * t615;
t631 = -qJD(4) * t642 + t686 * t663;
t639 = qJD(4) * pkin(4) - pkin(8) * t642;
t640 = t641 ^ 2;
t593 = -pkin(4) * t640 + pkin(8) * t631 - qJD(4) * t639 + t597;
t671 = sin(qJ(5));
t675 = cos(qJ(5));
t590 = t592 * t675 - t593 * t671;
t624 = t641 * t675 - t642 * t671;
t604 = qJD(5) * t624 + t631 * t671 + t632 * t675;
t625 = t641 * t671 + t642 * t675;
t610 = -mrSges(6,1) * t624 + mrSges(6,2) * t625;
t666 = qJD(4) + qJD(5);
t619 = -mrSges(6,2) * t666 + mrSges(6,3) * t624;
t662 = qJDD(4) + qJDD(5);
t587 = m(6) * t590 + mrSges(6,1) * t662 - t604 * mrSges(6,3) - t610 * t625 + t619 * t666;
t591 = t592 * t671 + t593 * t675;
t603 = -qJD(5) * t625 + t631 * t675 - t632 * t671;
t620 = mrSges(6,1) * t666 - mrSges(6,3) * t625;
t588 = m(6) * t591 - mrSges(6,2) * t662 + t603 * mrSges(6,3) + t610 * t624 - t620 * t666;
t577 = t675 * t587 + t671 * t588;
t629 = -mrSges(5,1) * t641 + mrSges(5,2) * t642;
t637 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t641;
t575 = m(5) * t596 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t632 + qJD(4) * t637 - t629 * t642 + t577;
t638 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t642;
t694 = -t587 * t671 + t675 * t588;
t576 = m(5) * t597 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t631 - qJD(4) * t638 + t629 * t641 + t694;
t571 = t676 * t575 + t672 * t576;
t617 = -t633 * t669 + t699;
t688 = mrSges(4,3) * t663 + (-mrSges(4,1) * t670 + t706) * t661;
t569 = m(4) * t617 - t688 * t669 + t571;
t695 = -t575 * t672 + t676 * t576;
t570 = m(4) * t618 + t688 * t670 + t695;
t696 = -t569 * t669 + t670 * t570;
t561 = m(3) * t636 - mrSges(3,1) * t661 - mrSges(3,2) * t663 + t696;
t635 = t648 * t677 - t673 * t649;
t690 = qJDD(3) - t635;
t630 = -pkin(2) * t663 - qJ(3) * t661 + t690;
t664 = t669 ^ 2;
t616 = (-pkin(2) - t707) * t663 + (-qJ(3) + (-t664 - t665) * pkin(7)) * t661 + t690;
t595 = -pkin(4) * t631 - pkin(8) * t640 + t639 * t642 + t616;
t689 = m(6) * t595 - t603 * mrSges(6,1) + t604 * mrSges(6,2) - t624 * t619 + t625 * t620;
t682 = m(5) * t616 - t631 * mrSges(5,1) + mrSges(5,2) * t632 - t641 * t637 + t638 * t642 + t689;
t681 = -m(4) * t630 + mrSges(4,1) * t703 - t682 + (t661 * t664 + t704) * mrSges(4,3);
t581 = m(3) * t635 - mrSges(3,2) * t661 + t681 + (mrSges(3,1) - t706) * t663;
t558 = t673 * t561 + t677 * t581;
t555 = m(2) * t653 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t679 + t558;
t697 = t677 * t561 - t581 * t673;
t556 = m(2) * t654 - mrSges(2,1) * t679 - qJDD(1) * mrSges(2,2) + t697;
t702 = t678 * t555 + t674 * t556;
t563 = t670 * t569 + t669 * t570;
t698 = -t555 * t674 + t678 * t556;
t693 = Ifges(4,1) * t669 + Ifges(4,4) * t670;
t692 = Ifges(4,4) * t669 + Ifges(4,2) * t670;
t605 = Ifges(6,5) * t625 + Ifges(6,6) * t624 + Ifges(6,3) * t666;
t607 = Ifges(6,1) * t625 + Ifges(6,4) * t624 + Ifges(6,5) * t666;
t578 = -mrSges(6,1) * t595 + mrSges(6,3) * t591 + Ifges(6,4) * t604 + Ifges(6,2) * t603 + Ifges(6,6) * t662 - t605 * t625 + t607 * t666;
t606 = Ifges(6,4) * t625 + Ifges(6,2) * t624 + Ifges(6,6) * t666;
t579 = mrSges(6,2) * t595 - mrSges(6,3) * t590 + Ifges(6,1) * t604 + Ifges(6,4) * t603 + Ifges(6,5) * t662 + t605 * t624 - t606 * t666;
t621 = Ifges(5,5) * t642 + Ifges(5,6) * t641 + Ifges(5,3) * qJD(4);
t623 = Ifges(5,1) * t642 + Ifges(5,4) * t641 + Ifges(5,5) * qJD(4);
t564 = -mrSges(5,1) * t616 + mrSges(5,3) * t597 + Ifges(5,4) * t632 + Ifges(5,2) * t631 + Ifges(5,6) * qJDD(4) - pkin(4) * t689 + pkin(8) * t694 + qJD(4) * t623 + t675 * t578 + t671 * t579 - t642 * t621;
t622 = Ifges(5,4) * t642 + Ifges(5,2) * t641 + Ifges(5,6) * qJD(4);
t565 = mrSges(5,2) * t616 - mrSges(5,3) * t596 + Ifges(5,1) * t632 + Ifges(5,4) * t631 + Ifges(5,5) * qJDD(4) - pkin(8) * t577 - qJD(4) * t622 - t578 * t671 + t579 * t675 + t621 * t641;
t549 = -mrSges(4,1) * t630 + mrSges(4,3) * t618 - pkin(3) * t682 + pkin(7) * t695 + t676 * t564 + t672 * t565 + t692 * t663 - t669 * t705;
t551 = mrSges(4,2) * t630 - mrSges(4,3) * t617 - pkin(7) * t571 - t564 * t672 + t565 * t676 + t693 * t663 + t670 * t705;
t583 = t663 * t706 - t681;
t685 = mrSges(3,1) * t635 - mrSges(3,2) * t636 + Ifges(3,3) * t663 - pkin(2) * t583 + qJ(3) * t696 + t670 * t549 + t669 * t551;
t684 = -mrSges(6,1) * t590 + mrSges(6,2) * t591 - Ifges(6,5) * t604 - Ifges(6,6) * t603 - Ifges(6,3) * t662 - t625 * t606 + t624 * t607;
t683 = mrSges(2,1) * t653 - mrSges(2,2) * t654 + Ifges(2,3) * qJDD(1) + pkin(1) * t558 + t685;
t680 = mrSges(5,1) * t596 - mrSges(5,2) * t597 + Ifges(5,5) * t632 + Ifges(5,6) * t631 + Ifges(5,3) * qJDD(4) + pkin(4) * t577 + t642 * t622 - t641 * t623 - t684;
t547 = -t680 + (Ifges(3,6) - t691) * t663 + mrSges(3,3) * t636 - mrSges(4,1) * t617 + mrSges(4,2) * t618 - pkin(3) * t571 + mrSges(3,1) * g(3) - pkin(2) * t563 + (-t669 * t692 + t670 * t693 + Ifges(3,5)) * t661;
t546 = -mrSges(3,2) * g(3) - mrSges(3,3) * t635 + Ifges(3,5) * t663 - Ifges(3,6) * t661 - qJ(3) * t563 - t549 * t669 + t551 * t670;
t545 = -mrSges(2,2) * g(3) - mrSges(2,3) * t653 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t679 - pkin(6) * t558 + t546 * t677 - t547 * t673;
t544 = Ifges(2,6) * qJDD(1) + t679 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t654 + t673 * t546 + t677 * t547 - pkin(1) * (-m(3) * g(3) + t563) + pkin(6) * t697;
t1 = [-m(1) * g(1) + t698; -m(1) * g(2) + t702; (-m(1) - m(2) - m(3)) * g(3) + t563; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t702 - t674 * t544 + t678 * t545; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t698 + t678 * t544 + t674 * t545; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t683; t683; t685; t583; t680; -t684;];
tauJB = t1;
