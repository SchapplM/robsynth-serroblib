% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:54
% EndTime: 2019-12-31 17:40:55
% DurationCPUTime: 1.24s
% Computational Cost: add. (8118->212), mult. (15967->259), div. (0->0), fcn. (7437->6), ass. (0->92)
t694 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t677 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t676 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t693 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t675 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t692 = Ifges(4,3) + Ifges(5,2) + Ifges(6,3);
t650 = sin(qJ(3));
t652 = cos(qJ(3));
t617 = (-mrSges(5,1) * t652 - mrSges(5,3) * t650) * qJD(2);
t678 = qJD(2) * qJD(3);
t620 = t650 * qJDD(2) + t652 * t678;
t647 = -g(3) + qJDD(1);
t648 = sin(pkin(7));
t649 = cos(pkin(7));
t626 = t648 * g(1) - t649 * g(2);
t627 = -t649 * g(1) - t648 * g(2);
t651 = sin(qJ(2));
t653 = cos(qJ(2));
t587 = t651 * t626 + t653 * t627;
t655 = qJD(2) ^ 2;
t584 = -t655 * pkin(2) + qJDD(2) * pkin(6) + t587;
t581 = t650 * t584;
t616 = (-pkin(3) * t652 - qJ(4) * t650) * qJD(2);
t654 = qJD(3) ^ 2;
t680 = qJD(2) * t650;
t662 = -t654 * qJ(4) + t616 * t680 + qJDD(4) + t581;
t673 = -0.2e1 * qJD(2) * qJD(5);
t686 = pkin(4) * t655;
t687 = pkin(3) + pkin(4);
t574 = t650 * t673 - t620 * qJ(5) - t687 * qJDD(3) + (qJ(5) * t678 - t650 * t686 - t647) * t652 + t662;
t618 = (mrSges(6,1) * t652 + mrSges(6,2) * t650) * qJD(2);
t679 = qJD(2) * t652;
t632 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t679;
t569 = m(6) * t574 - qJDD(3) * mrSges(6,1) - t620 * mrSges(6,3) - qJD(3) * t632 - t618 * t680;
t682 = t652 * t647;
t578 = -qJDD(3) * pkin(3) + t662 - t682;
t634 = mrSges(5,2) * t679 + qJD(3) * mrSges(5,3);
t658 = -m(5) * t578 + qJDD(3) * mrSges(5,1) + qJD(3) * t634 - t569;
t567 = t620 * mrSges(5,2) + t617 * t680 - t658;
t580 = t652 * t584 + t650 * t647;
t688 = 2 * qJD(4);
t577 = -t654 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t688 + t616 * t679 + t580;
t621 = t652 * qJDD(2) - t650 * t678;
t628 = -qJD(3) * pkin(4) - qJ(5) * t680;
t646 = t652 ^ 2;
t573 = -t621 * qJ(5) + qJD(3) * t628 - t646 * t686 + t652 * t673 + t577;
t579 = -t581 + t682;
t631 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t680;
t629 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t680;
t665 = m(6) * t573 + qJDD(3) * mrSges(6,2) - t621 * mrSges(6,3) + qJD(3) * t629;
t660 = m(5) * t577 + qJDD(3) * mrSges(5,3) + qJD(3) * t631 + t617 * t679 + t665;
t670 = t676 * qJD(3) + (t694 * t650 + t677 * t652) * qJD(2);
t671 = -t675 * qJD(3) + (-t677 * t650 - t693 * t652) * qJD(2);
t691 = -(t671 * t650 + t670 * t652) * qJD(2) + t692 * qJDD(3) + t676 * t620 + t675 * t621 + mrSges(4,1) * t579 - mrSges(5,1) * t578 - mrSges(6,1) * t574 - mrSges(4,2) * t580 + mrSges(6,2) * t573 + mrSges(5,3) * t577 - pkin(3) * t567 - pkin(4) * t569 + qJ(4) * (t621 * mrSges(5,2) - t618 * t679 + t660);
t685 = t655 * pkin(6);
t684 = mrSges(4,3) + mrSges(5,2);
t683 = qJ(4) * t652;
t619 = (-mrSges(4,1) * t652 + mrSges(4,2) * t650) * qJD(2);
t630 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t680;
t563 = m(4) * t580 - qJDD(3) * mrSges(4,2) - qJD(3) * t630 + t684 * t621 + (-t618 + t619) * t679 + t660;
t633 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t679;
t564 = m(4) * t579 + qJDD(3) * mrSges(4,1) + qJD(3) * t633 - t684 * t620 + (-t617 - t619) * t680 + t658;
t666 = t652 * t563 - t650 * t564;
t554 = m(3) * t587 - t655 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t666;
t586 = t653 * t626 - t651 * t627;
t663 = qJDD(2) * pkin(2) + t586;
t661 = -t620 * qJ(4) - t663;
t571 = qJDD(5) + (-qJ(5) * t646 + pkin(6)) * t655 + t687 * t621 + (qJD(3) * t683 + (-pkin(3) * qJD(3) + t628 + t688) * t650) * qJD(2) - t661;
t568 = m(6) * t571 + t621 * mrSges(6,1) + t620 * mrSges(6,2) + t629 * t680 + t632 * t679;
t575 = -t621 * pkin(3) - t685 + (-0.2e1 * qJD(4) * t650 + (pkin(3) * t650 - t683) * qJD(3)) * qJD(2) + t661;
t565 = m(5) * t575 - t621 * mrSges(5,1) - t620 * mrSges(5,3) - t631 * t680 - t634 * t679 - t568;
t583 = -t663 - t685;
t656 = -m(4) * t583 + t621 * mrSges(4,1) - t620 * mrSges(4,2) - t630 * t680 + t633 * t679 - t565;
t558 = m(3) * t586 + qJDD(2) * mrSges(3,1) - t655 * mrSges(3,2) + t656;
t549 = t651 * t554 + t653 * t558;
t547 = m(2) * t626 + t549;
t667 = t653 * t554 - t651 * t558;
t548 = m(2) * t627 + t667;
t681 = t649 * t547 + t648 * t548;
t556 = t650 * t563 + t652 * t564;
t672 = -t692 * qJD(3) + (-t676 * t650 - t675 * t652) * qJD(2);
t669 = m(3) * t647 + t556;
t668 = -t648 * t547 + t649 * t548;
t664 = m(2) * t647 + t669;
t543 = -mrSges(4,1) * t583 + mrSges(4,3) * t580 - mrSges(5,1) * t575 + mrSges(5,2) * t577 + mrSges(6,1) * t571 - mrSges(6,3) * t573 + pkin(4) * t568 - qJ(5) * t665 - pkin(3) * t565 + t693 * t621 + t677 * t620 + t675 * qJDD(3) + t670 * qJD(3) + (qJ(5) * t618 * t652 + t672 * t650) * qJD(2);
t551 = mrSges(4,2) * t583 + mrSges(5,2) * t578 + mrSges(6,2) * t571 - mrSges(4,3) * t579 - mrSges(5,3) * t575 - mrSges(6,3) * t574 - qJ(4) * t565 - qJ(5) * t569 + t671 * qJD(3) + t676 * qJDD(3) + t694 * t620 + t677 * t621 - t672 * t679;
t659 = mrSges(3,1) * t586 - mrSges(3,2) * t587 + Ifges(3,3) * qJDD(2) + pkin(2) * t656 + pkin(6) * t666 + t652 * t543 + t650 * t551;
t541 = -mrSges(3,1) * t647 + mrSges(3,3) * t587 + t655 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t556 - t691;
t540 = mrSges(3,2) * t647 - mrSges(3,3) * t586 + Ifges(3,5) * qJDD(2) - t655 * Ifges(3,6) - pkin(6) * t556 - t650 * t543 + t652 * t551;
t539 = mrSges(2,2) * t647 - mrSges(2,3) * t626 - pkin(5) * t549 + t653 * t540 - t651 * t541;
t538 = -mrSges(2,1) * t647 + mrSges(2,3) * t627 - pkin(1) * t669 + pkin(5) * t667 + t651 * t540 + t653 * t541;
t1 = [-m(1) * g(1) + t668; -m(1) * g(2) + t681; -m(1) * g(3) + t664; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t681 - t648 * t538 + t649 * t539; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t668 + t649 * t538 + t648 * t539; -mrSges(1,1) * g(2) + mrSges(2,1) * t626 + mrSges(1,2) * g(1) - mrSges(2,2) * t627 + pkin(1) * t549 + t659; t664; t659; t691; t567; t568;];
tauJB = t1;
