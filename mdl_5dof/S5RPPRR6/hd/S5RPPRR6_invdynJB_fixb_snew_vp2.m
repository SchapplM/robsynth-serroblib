% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:57
% EndTime: 2019-12-31 17:58:01
% DurationCPUTime: 3.89s
% Computational Cost: add. (42717->248), mult. (93230->310), div. (0->0), fcn. (61188->10), ass. (0->113)
t656 = qJD(1) ^ 2;
t645 = sin(pkin(9));
t647 = cos(pkin(9));
t650 = sin(qJ(4));
t653 = cos(qJ(4));
t663 = t645 * t650 - t647 * t653;
t619 = t663 * qJD(1);
t664 = t645 * t653 + t647 * t650;
t620 = t664 * qJD(1);
t675 = t620 * qJD(4);
t609 = -t663 * qJDD(1) - t675;
t682 = pkin(3) * t647;
t681 = mrSges(4,2) * t645;
t641 = t647 ^ 2;
t680 = t641 * t656;
t651 = sin(qJ(1));
t654 = cos(qJ(1));
t629 = t651 * g(1) - t654 * g(2);
t626 = qJDD(1) * pkin(1) + t629;
t630 = -t654 * g(1) - t651 * g(2);
t627 = -t656 * pkin(1) + t630;
t646 = sin(pkin(8));
t648 = cos(pkin(8));
t613 = t646 * t626 + t648 * t627;
t603 = -t656 * pkin(2) + qJDD(1) * qJ(3) + t613;
t644 = -g(3) + qJDD(2);
t674 = qJD(1) * qJD(3);
t678 = t647 * t644 - 0.2e1 * t645 * t674;
t588 = (-pkin(6) * qJDD(1) + t656 * t682 - t603) * t645 + t678;
t594 = t645 * t644 + (t603 + 0.2e1 * t674) * t647;
t673 = qJDD(1) * t647;
t591 = -pkin(3) * t680 + pkin(6) * t673 + t594;
t580 = t650 * t588 + t653 * t591;
t608 = t619 * pkin(4) - t620 * pkin(7);
t655 = qJD(4) ^ 2;
t577 = -t655 * pkin(4) + qJDD(4) * pkin(7) - t619 * t608 + t580;
t640 = t645 ^ 2;
t612 = t648 * t626 - t646 * t627;
t665 = qJDD(3) - t612;
t592 = (-pkin(2) - t682) * qJDD(1) + (-qJ(3) + (-t640 - t641) * pkin(6)) * t656 + t665;
t676 = t619 * qJD(4);
t610 = t664 * qJDD(1) - t676;
t578 = (-t610 + t676) * pkin(7) + (-t609 + t675) * pkin(4) + t592;
t649 = sin(qJ(5));
t652 = cos(qJ(5));
t574 = -t649 * t577 + t652 * t578;
t614 = t652 * qJD(4) - t649 * t620;
t590 = t614 * qJD(5) + t649 * qJDD(4) + t652 * t610;
t615 = t649 * qJD(4) + t652 * t620;
t595 = -t614 * mrSges(6,1) + t615 * mrSges(6,2);
t618 = qJD(5) + t619;
t596 = -t618 * mrSges(6,2) + t614 * mrSges(6,3);
t607 = qJDD(5) - t609;
t571 = m(6) * t574 + t607 * mrSges(6,1) - t590 * mrSges(6,3) - t615 * t595 + t618 * t596;
t575 = t652 * t577 + t649 * t578;
t589 = -t615 * qJD(5) + t652 * qJDD(4) - t649 * t610;
t597 = t618 * mrSges(6,1) - t615 * mrSges(6,3);
t572 = m(6) * t575 - t607 * mrSges(6,2) + t589 * mrSges(6,3) + t614 * t595 - t618 * t597;
t563 = -t649 * t571 + t652 * t572;
t605 = t619 * mrSges(5,1) + t620 * mrSges(5,2);
t617 = qJD(4) * mrSges(5,1) - t620 * mrSges(5,3);
t561 = m(5) * t580 - qJDD(4) * mrSges(5,2) + t609 * mrSges(5,3) - qJD(4) * t617 - t619 * t605 + t563;
t579 = t653 * t588 - t650 * t591;
t576 = -qJDD(4) * pkin(4) - t655 * pkin(7) + t620 * t608 - t579;
t573 = -m(6) * t576 + t589 * mrSges(6,1) - t590 * mrSges(6,2) + t614 * t596 - t615 * t597;
t616 = -qJD(4) * mrSges(5,2) - t619 * mrSges(5,3);
t567 = m(5) * t579 + qJDD(4) * mrSges(5,1) - t610 * mrSges(5,3) + qJD(4) * t616 - t620 * t605 + t573;
t554 = t650 * t561 + t653 * t567;
t593 = -t645 * t603 + t678;
t662 = mrSges(4,3) * qJDD(1) + t656 * (-mrSges(4,1) * t647 + t681);
t552 = m(4) * t593 - t662 * t645 + t554;
t669 = t653 * t561 - t650 * t567;
t553 = m(4) * t594 + t662 * t647 + t669;
t670 = -t645 * t552 + t647 * t553;
t543 = m(3) * t613 - t656 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t670;
t599 = -qJDD(1) * pkin(2) - t656 * qJ(3) + t665;
t562 = t652 * t571 + t649 * t572;
t661 = m(5) * t592 - t609 * mrSges(5,1) + t610 * mrSges(5,2) + t619 * t616 + t620 * t617 + t562;
t659 = -m(4) * t599 + mrSges(4,1) * t673 - t661 + (t640 * t656 + t680) * mrSges(4,3);
t556 = -t656 * mrSges(3,2) + m(3) * t612 + t659 + (mrSges(3,1) - t681) * qJDD(1);
t540 = t646 * t543 + t648 * t556;
t537 = m(2) * t629 + qJDD(1) * mrSges(2,1) - t656 * mrSges(2,2) + t540;
t671 = t648 * t543 - t646 * t556;
t538 = m(2) * t630 - t656 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t671;
t679 = t654 * t537 + t651 * t538;
t546 = t647 * t552 + t645 * t553;
t666 = Ifges(4,5) * t645 + Ifges(4,6) * t647;
t677 = t656 * t666;
t544 = m(3) * t644 + t546;
t672 = -t651 * t537 + t654 * t538;
t668 = Ifges(4,1) * t645 + Ifges(4,4) * t647;
t667 = Ifges(4,4) * t645 + Ifges(4,2) * t647;
t582 = Ifges(6,5) * t615 + Ifges(6,6) * t614 + Ifges(6,3) * t618;
t584 = Ifges(6,1) * t615 + Ifges(6,4) * t614 + Ifges(6,5) * t618;
t564 = -mrSges(6,1) * t576 + mrSges(6,3) * t575 + Ifges(6,4) * t590 + Ifges(6,2) * t589 + Ifges(6,6) * t607 - t615 * t582 + t618 * t584;
t583 = Ifges(6,4) * t615 + Ifges(6,2) * t614 + Ifges(6,6) * t618;
t565 = mrSges(6,2) * t576 - mrSges(6,3) * t574 + Ifges(6,1) * t590 + Ifges(6,4) * t589 + Ifges(6,5) * t607 + t614 * t582 - t618 * t583;
t600 = Ifges(5,5) * t620 - Ifges(5,6) * t619 + Ifges(5,3) * qJD(4);
t601 = Ifges(5,4) * t620 - Ifges(5,2) * t619 + Ifges(5,6) * qJD(4);
t547 = mrSges(5,2) * t592 - mrSges(5,3) * t579 + Ifges(5,1) * t610 + Ifges(5,4) * t609 + Ifges(5,5) * qJDD(4) - pkin(7) * t562 - qJD(4) * t601 - t649 * t564 + t652 * t565 - t619 * t600;
t602 = Ifges(5,1) * t620 - Ifges(5,4) * t619 + Ifges(5,5) * qJD(4);
t658 = mrSges(6,1) * t574 - mrSges(6,2) * t575 + Ifges(6,5) * t590 + Ifges(6,6) * t589 + Ifges(6,3) * t607 + t615 * t583 - t614 * t584;
t548 = -mrSges(5,1) * t592 + mrSges(5,3) * t580 + Ifges(5,4) * t610 + Ifges(5,2) * t609 + Ifges(5,6) * qJDD(4) - pkin(4) * t562 + qJD(4) * t602 - t620 * t600 - t658;
t531 = -mrSges(4,1) * t599 + mrSges(4,3) * t594 - pkin(3) * t661 + pkin(6) * t669 + t667 * qJDD(1) + t650 * t547 + t653 * t548 - t645 * t677;
t533 = mrSges(4,2) * t599 - mrSges(4,3) * t593 - pkin(6) * t554 + t668 * qJDD(1) + t653 * t547 - t650 * t548 + t647 * t677;
t558 = qJDD(1) * t681 - t659;
t660 = mrSges(2,1) * t629 + mrSges(3,1) * t612 - mrSges(2,2) * t630 - mrSges(3,2) * t613 + pkin(1) * t540 - pkin(2) * t558 + qJ(3) * t670 + t647 * t531 + t645 * t533 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t657 = mrSges(5,1) * t579 - mrSges(5,2) * t580 + Ifges(5,5) * t610 + Ifges(5,6) * t609 + Ifges(5,3) * qJDD(4) + pkin(4) * t573 + pkin(7) * t563 + t652 * t564 + t649 * t565 + t620 * t601 + t619 * t602;
t529 = -mrSges(3,1) * t644 - pkin(2) * t546 - pkin(3) * t554 - mrSges(4,1) * t593 + mrSges(4,2) * t594 - t657 + (Ifges(3,6) - t666) * qJDD(1) + mrSges(3,3) * t613 + (-t645 * t667 + t647 * t668 + Ifges(3,5)) * t656;
t528 = mrSges(3,2) * t644 - mrSges(3,3) * t612 + Ifges(3,5) * qJDD(1) - t656 * Ifges(3,6) - qJ(3) * t546 - t645 * t531 + t647 * t533;
t527 = -mrSges(2,2) * g(3) - mrSges(2,3) * t629 + Ifges(2,5) * qJDD(1) - t656 * Ifges(2,6) - qJ(2) * t540 + t648 * t528 - t646 * t529;
t526 = mrSges(2,1) * g(3) + mrSges(2,3) * t630 + t656 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t544 + qJ(2) * t671 + t646 * t528 + t648 * t529;
t1 = [-m(1) * g(1) + t672; -m(1) * g(2) + t679; (-m(1) - m(2)) * g(3) + t544; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t679 - t651 * t526 + t654 * t527; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t672 + t654 * t526 + t651 * t527; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t660; t660; t544; t558; t657; t658;];
tauJB = t1;
