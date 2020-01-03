% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:06
% EndTime: 2019-12-31 17:55:08
% DurationCPUTime: 1.64s
% Computational Cost: add. (11861->216), mult. (25916->249), div. (0->0), fcn. (15145->6), ass. (0->99)
t702 = Ifges(5,1) + Ifges(6,1);
t687 = Ifges(5,4) - Ifges(6,5);
t700 = Ifges(5,5) + Ifges(6,4);
t701 = Ifges(5,2) + Ifges(6,3);
t699 = Ifges(5,6) - Ifges(6,6);
t698 = Ifges(5,3) + Ifges(6,2);
t646 = sin(pkin(7));
t647 = cos(pkin(7));
t648 = sin(qJ(4));
t650 = cos(qJ(4));
t664 = t646 * t650 + t647 * t648;
t619 = t664 * qJD(1);
t663 = -t646 * t648 + t647 * t650;
t620 = t663 * qJD(1);
t697 = -t699 * qJD(4) + t701 * t619 - t687 * t620;
t696 = t700 * qJD(4) - t687 * t619 + t702 * t620;
t649 = sin(qJ(1));
t651 = cos(qJ(1));
t624 = t649 * g(1) - t651 * g(2);
t653 = qJD(1) ^ 2;
t660 = -t653 * qJ(2) + qJDD(2) - t624;
t685 = -pkin(1) - qJ(3);
t695 = -0.2e1 * qJD(1) * qJD(3) + t685 * qJDD(1) + t660;
t625 = -t651 * g(1) - t649 * g(2);
t694 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t625;
t616 = t653 * pkin(1) - t694;
t693 = -m(3) * t616 + t653 * mrSges(3,2) + qJDD(1) * mrSges(3,3);
t691 = pkin(3) * t653;
t690 = mrSges(2,1) - mrSges(3,2);
t689 = -mrSges(5,3) - mrSges(6,2);
t688 = -Ifges(3,4) + Ifges(2,5);
t686 = -Ifges(2,6) + Ifges(3,5);
t684 = mrSges(4,2) * t647;
t601 = t646 * g(3) + t695 * t647;
t580 = (-pkin(6) * qJDD(1) - t646 * t691) * t647 + t601;
t602 = -t647 * g(3) + t695 * t646;
t634 = t646 ^ 2;
t675 = qJDD(1) * t646;
t581 = -pkin(6) * t675 - t634 * t691 + t602;
t575 = t648 * t580 + t650 * t581;
t677 = t620 * qJD(4);
t603 = t664 * qJDD(1) + t677;
t612 = qJD(4) * mrSges(5,1) - t620 * mrSges(5,3);
t592 = t619 * pkin(4) - t620 * qJ(5);
t652 = qJD(4) ^ 2;
t571 = -t652 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t619 * t592 + t575;
t613 = -qJD(4) * mrSges(6,1) + t620 * mrSges(6,2);
t673 = m(6) * t571 + qJDD(4) * mrSges(6,3) + qJD(4) * t613;
t593 = t619 * mrSges(6,1) - t620 * mrSges(6,3);
t681 = -t619 * mrSges(5,1) - t620 * mrSges(5,2) - t593;
t561 = m(5) * t575 - qJDD(4) * mrSges(5,2) - qJD(4) * t612 + t689 * t603 + t681 * t619 + t673;
t574 = t650 * t580 - t648 * t581;
t678 = t619 * qJD(4);
t604 = t663 * qJDD(1) - t678;
t611 = -qJD(4) * mrSges(5,2) - t619 * mrSges(5,3);
t572 = -qJDD(4) * pkin(4) - t652 * qJ(5) + t620 * t592 + qJDD(5) - t574;
t614 = -t619 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t668 = -m(6) * t572 + qJDD(4) * mrSges(6,1) + qJD(4) * t614;
t562 = m(5) * t574 + qJDD(4) * mrSges(5,1) + qJD(4) * t611 + t689 * t604 + t681 * t620 + t668;
t552 = t648 * t561 + t650 * t562;
t662 = -qJDD(1) * mrSges(4,3) - t653 * (mrSges(4,1) * t646 + t684);
t550 = m(4) * t601 + t662 * t647 + t552;
t669 = t650 * t561 - t648 * t562;
t551 = m(4) * t602 + t662 * t646 + t669;
t546 = t647 * t550 + t646 * t551;
t618 = -qJDD(1) * pkin(1) + t660;
t658 = -m(3) * t618 + t653 * mrSges(3,3) - t546;
t542 = m(2) * t624 - t653 * mrSges(2,2) + t690 * qJDD(1) + t658;
t659 = qJDD(3) + t694;
t608 = t685 * t653 + t659;
t680 = -t647 ^ 2 - t634;
t583 = pkin(3) * t675 + (t680 * pkin(6) + t685) * t653 + t659;
t568 = -0.2e1 * qJD(5) * t620 + (-t604 + t678) * qJ(5) + (t603 + t677) * pkin(4) + t583;
t563 = m(6) * t568 + t603 * mrSges(6,1) - t604 * mrSges(6,3) - t620 * t613 + t619 * t614;
t657 = m(5) * t583 + t603 * mrSges(5,1) + t604 * mrSges(5,2) + t619 * t611 + t620 * t612 + t563;
t655 = m(4) * t608 + mrSges(4,1) * t675 + qJDD(1) * t684 + t657;
t672 = t680 * mrSges(4,3);
t555 = (-mrSges(2,1) + t672) * t653 + t655 - qJDD(1) * mrSges(2,2) + m(2) * t625 + t693;
t683 = t651 * t542 + t649 * t555;
t682 = -t698 * qJD(4) + t699 * t619 - t700 * t620;
t665 = Ifges(4,5) * t647 - Ifges(4,6) * t646;
t679 = t653 * t665;
t671 = -t649 * t542 + t651 * t555;
t670 = -t646 * t550 + t647 * t551;
t667 = Ifges(4,1) * t647 - Ifges(4,4) * t646;
t666 = Ifges(4,4) * t647 - Ifges(4,2) * t646;
t547 = -mrSges(5,1) * t583 - mrSges(6,1) * t568 + mrSges(6,2) * t571 + mrSges(5,3) * t575 - pkin(4) * t563 + t696 * qJD(4) + t699 * qJDD(4) - t701 * t603 + t687 * t604 + t682 * t620;
t548 = mrSges(5,2) * t583 + mrSges(6,2) * t572 - mrSges(5,3) * t574 - mrSges(6,3) * t568 - qJ(5) * t563 + t697 * qJD(4) + t700 * qJDD(4) - t687 * t603 + t702 * t604 + t682 * t619;
t538 = -mrSges(4,1) * t608 + mrSges(4,3) * t602 - pkin(3) * t657 + pkin(6) * t669 + t666 * qJDD(1) + t650 * t547 + t648 * t548 - t647 * t679;
t540 = mrSges(4,2) * t608 - mrSges(4,3) * t601 - pkin(6) * t552 + t667 * qJDD(1) - t648 * t547 + t650 * t548 - t646 * t679;
t544 = qJDD(1) * mrSges(3,2) - t658;
t557 = t653 * t672 + t655;
t656 = -mrSges(2,2) * t625 - mrSges(3,3) * t616 - pkin(1) * t544 - qJ(3) * t546 - t646 * t538 + t647 * t540 + qJ(2) * (t557 + t693) + mrSges(3,2) * t618 + mrSges(2,1) * t624 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t566 = t604 * mrSges(6,2) + t620 * t593 - t668;
t654 = mrSges(5,1) * t574 - mrSges(6,1) * t572 - mrSges(5,2) * t575 + mrSges(6,3) * t571 - pkin(4) * t566 + qJ(5) * t673 - t697 * t620 + (-qJ(5) * t593 + t696) * t619 + t700 * t604 + (-mrSges(6,2) * qJ(5) - t699) * t603 + t698 * qJDD(4);
t545 = -m(3) * g(3) + t670;
t537 = t654 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t665 + t688) * qJDD(1) - mrSges(2,3) * t624 + mrSges(3,1) * t618 + mrSges(4,1) * t601 - mrSges(4,2) * t602 + pkin(3) * t552 + pkin(2) * t546 - qJ(2) * t545 + (t646 * t667 + t647 * t666 + t686) * t653;
t536 = -mrSges(3,1) * t616 + mrSges(2,3) * t625 - pkin(1) * t545 + pkin(2) * t557 + t690 * g(3) - qJ(3) * t670 - t686 * qJDD(1) - t647 * t538 - t646 * t540 + t688 * t653;
t1 = [-m(1) * g(1) + t671; -m(1) * g(2) + t683; (-m(1) - m(2) - m(3)) * g(3) + t670; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t683 - t649 * t536 + t651 * t537; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t671 + t651 * t536 + t649 * t537; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t656; t656; t544; t557; t654; t566;];
tauJB = t1;
