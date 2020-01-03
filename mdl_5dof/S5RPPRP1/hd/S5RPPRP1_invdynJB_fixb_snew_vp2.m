% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:28
% EndTime: 2020-01-03 11:25:31
% DurationCPUTime: 2.28s
% Computational Cost: add. (17596->229), mult. (36515->286), div. (0->0), fcn. (20211->8), ass. (0->107)
t652 = sin(qJ(1));
t654 = cos(qJ(1));
t631 = -t654 * g(2) - t652 * g(3);
t626 = qJDD(1) * pkin(1) + t631;
t630 = -t652 * g(2) + t654 * g(3);
t655 = qJD(1) ^ 2;
t627 = -t655 * pkin(1) + t630;
t648 = sin(pkin(7));
t650 = cos(pkin(7));
t599 = t648 * t626 + t650 * t627;
t703 = -t655 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t599;
t702 = Ifges(5,1) + Ifges(6,1);
t695 = Ifges(5,4) + Ifges(6,4);
t694 = Ifges(5,5) + Ifges(6,5);
t701 = Ifges(5,2) + Ifges(6,2);
t693 = Ifges(5,6) + Ifges(6,6);
t700 = Ifges(5,3) + Ifges(6,3);
t646 = -g(1) + qJDD(2);
t647 = sin(pkin(8));
t649 = cos(pkin(8));
t588 = t649 * t646 - t703 * t647;
t651 = sin(qJ(4));
t653 = cos(qJ(4));
t680 = t649 * qJD(1);
t633 = qJD(4) - t680;
t681 = t647 * qJD(1);
t684 = (-t695 * t651 + t702 * t653) * t681 + t694 * t633;
t685 = (-t701 * t651 + t695 * t653) * t681 + t693 * t633;
t699 = t684 * t651 + t685 * t653;
t615 = (mrSges(6,1) * t651 + mrSges(6,2) * t653) * t681;
t678 = qJD(1) * qJD(4);
t618 = (qJDD(1) * t653 - t651 * t678) * t647;
t673 = t653 * t681;
t589 = t647 * t646 + t703 * t649;
t664 = -pkin(3) * t649 - pkin(6) * t647;
t625 = t664 * qJD(1);
t587 = t625 * t680 + t589;
t598 = t650 * t626 - t648 * t627;
t660 = -t655 * qJ(3) + qJDD(3) - t598;
t592 = (-pkin(2) + t664) * qJDD(1) + t660;
t591 = t653 * t592;
t677 = t649 * qJDD(1);
t632 = qJDD(4) - t677;
t665 = -0.2e1 * qJD(5) * t681;
t689 = t647 ^ 2 * t655;
t579 = t653 * t665 + t632 * pkin(4) - t618 * qJ(5) + t591 + (-pkin(4) * t653 * t689 - qJ(5) * t633 * t681 - t587) * t651;
t674 = t651 * t681;
t610 = -t633 * mrSges(6,2) - mrSges(6,3) * t674;
t675 = m(6) * t579 + t632 * mrSges(6,1) + t633 * t610;
t576 = -t618 * mrSges(6,3) - t615 * t673 + t675;
t583 = t653 * t587 + t651 * t592;
t612 = t633 * pkin(4) - qJ(5) * t673;
t617 = (-qJDD(1) * t651 - t653 * t678) * t647;
t676 = t651 ^ 2 * t689;
t581 = -pkin(4) * t676 + t617 * qJ(5) - t633 * t612 + t651 * t665 + t583;
t582 = -t651 * t587 + t591;
t698 = mrSges(5,1) * t582 + mrSges(6,1) * t579 - mrSges(5,2) * t583 - mrSges(6,2) * t581 + pkin(4) * t576 + t693 * t617 + t694 * t618 + t700 * t632;
t697 = t649 ^ 2;
t696 = -mrSges(5,2) - mrSges(6,2);
t691 = mrSges(4,2) * t647;
t690 = t618 * mrSges(6,2);
t623 = (-mrSges(4,1) * t649 + t691) * qJD(1);
t611 = -t633 * mrSges(5,2) - mrSges(5,3) * t674;
t663 = (-t615 - (mrSges(5,1) * t651 + mrSges(5,2) * t653) * t681) * t681;
t572 = m(5) * t582 + t632 * mrSges(5,1) + t633 * t611 + (-mrSges(5,3) - mrSges(6,3)) * t618 + t653 * t663 + t675;
t613 = t633 * mrSges(6,1) - mrSges(6,3) * t673;
t683 = -t633 * mrSges(5,1) + mrSges(5,3) * t673 - t613;
t687 = m(6) * t581 + t617 * mrSges(6,3);
t573 = m(5) * t583 + t617 * mrSges(5,3) + t696 * t632 + t683 * t633 + t651 * t663 + t687;
t667 = -t651 * t572 + t653 * t573;
t682 = qJDD(1) * mrSges(4,3);
t566 = m(4) * t589 + (qJD(1) * t623 + t682) * t649 + t667;
t659 = t683 * t653 + (-t610 - t611) * t651;
t586 = t625 * t681 - t588;
t584 = -t617 * pkin(4) - qJ(5) * t676 + t612 * t673 + qJDD(5) + t586;
t671 = m(6) * t584 - t617 * mrSges(6,1);
t661 = -m(5) * t586 + t617 * mrSges(5,1) - t671;
t575 = m(4) * t588 + t696 * t618 + (-t682 + (-t623 + t659) * qJD(1)) * t647 + t661;
t668 = t649 * t566 - t647 * t575;
t557 = m(3) * t599 - t655 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t668;
t570 = t653 * t572 + t651 * t573;
t595 = -qJDD(1) * pkin(2) + t660;
t658 = -m(4) * t595 + mrSges(4,1) * t677 - t570 + (t655 * t697 + t689) * mrSges(4,3);
t563 = m(3) * t598 - t655 * mrSges(3,2) + (mrSges(3,1) - t691) * qJDD(1) + t658;
t669 = t650 * t557 - t648 * t563;
t549 = m(2) * t630 - t655 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t669;
t552 = t648 * t557 + t650 * t563;
t550 = m(2) * t631 + qJDD(1) * mrSges(2,1) - t655 * mrSges(2,2) + t552;
t688 = t652 * t549 + t654 * t550;
t560 = t647 * t566 + t649 * t575;
t686 = (t693 * t651 - t694 * t653) * t681 - t700 * t633;
t558 = m(3) * t646 + t560;
t670 = -t654 * t549 + t652 * t550;
t662 = Ifges(4,5) * t647 + Ifges(4,6) * t649;
t577 = t690 + (t610 * t651 + t613 * t653) * t681 + t671;
t561 = -mrSges(5,1) * t586 + mrSges(5,3) * t583 - mrSges(6,1) * t584 + mrSges(6,3) * t581 - pkin(4) * t577 + qJ(5) * t687 + (-qJ(5) * t613 + t684) * t633 + (-qJ(5) * mrSges(6,2) + t693) * t632 + t695 * t618 + t701 * t617 + (-qJ(5) * t615 * t651 + t686 * t653) * t681;
t569 = mrSges(5,2) * t586 + mrSges(6,2) * t584 - mrSges(5,3) * t582 - mrSges(6,3) * t579 - qJ(5) * t576 + t695 * t617 + t702 * t618 + t694 * t632 - t685 * t633 + t686 * t674;
t624 = t662 * qJD(1);
t545 = t624 * t680 + mrSges(4,2) * t595 - mrSges(4,3) * t588 - pkin(6) * t570 - t651 * t561 + t653 * t569 + (Ifges(4,1) * t647 + Ifges(4,4) * t649) * qJDD(1);
t554 = Ifges(4,2) * t677 - mrSges(4,1) * t595 + mrSges(4,3) * t589 - pkin(3) * t570 + (Ifges(4,4) * qJDD(1) + (-t624 - t699) * qJD(1)) * t647 - t698;
t568 = qJDD(1) * t691 - t658;
t656 = mrSges(2,1) * t631 + mrSges(3,1) * t598 - mrSges(2,2) * t630 - mrSges(3,2) * t599 + pkin(1) * t552 - pkin(2) * t568 + qJ(3) * t668 + t647 * t545 + t649 * t554 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t543 = t655 * Ifges(3,5) - mrSges(3,1) * t646 + mrSges(3,3) * t599 - mrSges(4,1) * t588 + mrSges(4,2) * t589 - t651 * t569 - t653 * t561 - pkin(3) * (-t618 * mrSges(5,2) + t661 - t690) - pkin(6) * t667 - pkin(2) * t560 + (Ifges(3,6) - t662) * qJDD(1) + (Ifges(4,4) * t697 * qJD(1) + (-pkin(3) * t659 + (-Ifges(4,4) * t647 + (Ifges(4,1) - Ifges(4,2)) * t649) * qJD(1)) * t647) * qJD(1);
t542 = mrSges(3,2) * t646 - mrSges(3,3) * t598 + Ifges(3,5) * qJDD(1) - t655 * Ifges(3,6) - qJ(3) * t560 + t649 * t545 - t647 * t554;
t541 = -mrSges(2,2) * g(1) - mrSges(2,3) * t631 + Ifges(2,5) * qJDD(1) - t655 * Ifges(2,6) - qJ(2) * t552 + t650 * t542 - t648 * t543;
t540 = mrSges(2,1) * g(1) + mrSges(2,3) * t630 + t655 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t558 + qJ(2) * t669 + t648 * t542 + t650 * t543;
t1 = [(-m(1) - m(2)) * g(1) + t558; -m(1) * g(2) + t688; -m(1) * g(3) + t670; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t656; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t670 + t654 * t540 + t652 * t541; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t688 + t652 * t540 - t654 * t541; t656; t558; t568; t699 * t681 + t698; t577;];
tauJB = t1;
