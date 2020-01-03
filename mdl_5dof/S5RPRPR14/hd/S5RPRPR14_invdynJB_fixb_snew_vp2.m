% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR14_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR14_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:46
% EndTime: 2019-12-31 18:34:49
% DurationCPUTime: 2.66s
% Computational Cost: add. (27448->265), mult. (58747->326), div. (0->0), fcn. (35739->8), ass. (0->108)
t653 = sin(qJ(1));
t656 = cos(qJ(1));
t635 = -t656 * g(1) - t653 * g(2);
t667 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t635;
t649 = sin(pkin(8));
t650 = cos(pkin(8));
t652 = sin(qJ(3));
t655 = cos(qJ(3));
t616 = (t649 * t655 + t650 * t652) * qJD(1);
t685 = 2 * qJD(4);
t684 = -pkin(1) - pkin(6);
t683 = mrSges(2,1) - mrSges(3,2);
t682 = Ifges(2,5) - Ifges(3,4);
t681 = (-Ifges(2,6) + Ifges(3,5));
t634 = t653 * g(1) - t656 * g(2);
t658 = qJD(1) ^ 2;
t666 = -t658 * qJ(2) + qJDD(2) - t634;
t608 = qJDD(1) * t684 + t666;
t598 = t652 * g(3) + t655 * t608;
t677 = qJD(1) * qJD(3);
t675 = t652 * t677;
t629 = t655 * qJDD(1) - t675;
t577 = (-t629 - t675) * qJ(4) + (-t652 * t655 * t658 + qJDD(3)) * pkin(3) + t598;
t599 = -t655 * g(3) + t652 * t608;
t628 = -t652 * qJDD(1) - t655 * t677;
t678 = qJD(1) * t655;
t632 = qJD(3) * pkin(3) - qJ(4) * t678;
t646 = t652 ^ 2;
t578 = -t646 * t658 * pkin(3) + t628 * qJ(4) - qJD(3) * t632 + t599;
t567 = t649 * t577 + t650 * t578 - t616 * t685;
t679 = qJD(1) * t652;
t617 = -t649 * t679 + t650 * t678;
t589 = t616 * mrSges(5,1) + t617 * mrSges(5,2);
t595 = t650 * t628 - t649 * t629;
t607 = qJD(3) * mrSges(5,1) - t617 * mrSges(5,3);
t590 = t616 * pkin(4) - t617 * pkin(7);
t657 = qJD(3) ^ 2;
t563 = -t657 * pkin(4) + qJDD(3) * pkin(7) - t616 * t590 + t567;
t580 = -t628 * pkin(3) + qJDD(4) + t632 * t678 + (-qJ(4) * t646 + t684) * t658 + t667;
t596 = t649 * t628 + t650 * t629;
t564 = (qJD(3) * t616 - t596) * pkin(7) + (qJD(3) * t617 - t595) * pkin(4) + t580;
t651 = sin(qJ(5));
t654 = cos(qJ(5));
t560 = -t651 * t563 + t654 * t564;
t600 = t654 * qJD(3) - t651 * t617;
t574 = t600 * qJD(5) + t651 * qJDD(3) + t654 * t596;
t601 = t651 * qJD(3) + t654 * t617;
t582 = -t600 * mrSges(6,1) + t601 * mrSges(6,2);
t614 = qJD(5) + t616;
t583 = -t614 * mrSges(6,2) + t600 * mrSges(6,3);
t594 = qJDD(5) - t595;
t557 = m(6) * t560 + t594 * mrSges(6,1) - t574 * mrSges(6,3) - t601 * t582 + t614 * t583;
t561 = t654 * t563 + t651 * t564;
t573 = -t601 * qJD(5) + t654 * qJDD(3) - t651 * t596;
t584 = t614 * mrSges(6,1) - t601 * mrSges(6,3);
t558 = m(6) * t561 - t594 * mrSges(6,2) + t573 * mrSges(6,3) + t600 * t582 - t614 * t584;
t671 = -t651 * t557 + t654 * t558;
t544 = m(5) * t567 - qJDD(3) * mrSges(5,2) + t595 * mrSges(5,3) - qJD(3) * t607 - t616 * t589 + t671;
t670 = -t650 * t577 + t649 * t578;
t566 = -0.2e1 * qJD(4) * t617 - t670;
t606 = -qJD(3) * mrSges(5,2) - t616 * mrSges(5,3);
t562 = -qJDD(3) * pkin(4) - t657 * pkin(7) + (t685 + t590) * t617 + t670;
t664 = -m(6) * t562 + t573 * mrSges(6,1) - t574 * mrSges(6,2) + t600 * t583 - t601 * t584;
t553 = m(5) * t566 + qJDD(3) * mrSges(5,1) - t596 * mrSges(5,3) + qJD(3) * t606 - t617 * t589 + t664;
t537 = t649 * t544 + t650 * t553;
t627 = (mrSges(4,1) * t652 + mrSges(4,2) * t655) * qJD(1);
t631 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t679;
t534 = m(4) * t598 + qJDD(3) * mrSges(4,1) - t629 * mrSges(4,3) + qJD(3) * t631 - t627 * t678 + t537;
t633 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t678;
t672 = t650 * t544 - t649 * t553;
t535 = m(4) * t599 - qJDD(3) * mrSges(4,2) + t628 * mrSges(4,3) - qJD(3) * t633 - t627 * t679 + t672;
t530 = t655 * t534 + t652 * t535;
t615 = -qJDD(1) * pkin(1) + t666;
t665 = -m(3) * t615 + (t658 * mrSges(3,3)) - t530;
t526 = m(2) * t634 - (t658 * mrSges(2,2)) + qJDD(1) * t683 + t665;
t611 = t658 * pkin(1) - t667;
t547 = t654 * t557 + t651 * t558;
t545 = m(5) * t580 - t595 * mrSges(5,1) + t596 * mrSges(5,2) + t616 * t606 + t617 * t607 + t547;
t605 = t658 * t684 + t667;
t663 = -m(4) * t605 + t628 * mrSges(4,1) - t629 * mrSges(4,2) - t631 * t679 - t633 * t678 - t545;
t660 = -m(3) * t611 + (t658 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t663;
t540 = m(2) * t635 - (t658 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t660;
t680 = t656 * t526 + t653 * t540;
t674 = -t653 * t526 + t656 * t540;
t673 = -t652 * t534 + t655 * t535;
t568 = Ifges(6,5) * t601 + Ifges(6,6) * t600 + Ifges(6,3) * t614;
t570 = Ifges(6,1) * t601 + Ifges(6,4) * t600 + Ifges(6,5) * t614;
t550 = -mrSges(6,1) * t562 + mrSges(6,3) * t561 + Ifges(6,4) * t574 + Ifges(6,2) * t573 + Ifges(6,6) * t594 - t601 * t568 + t614 * t570;
t569 = Ifges(6,4) * t601 + Ifges(6,2) * t600 + Ifges(6,6) * t614;
t551 = mrSges(6,2) * t562 - mrSges(6,3) * t560 + Ifges(6,1) * t574 + Ifges(6,4) * t573 + Ifges(6,5) * t594 + t600 * t568 - t614 * t569;
t585 = Ifges(5,5) * t617 - Ifges(5,6) * t616 + (Ifges(5,3) * qJD(3));
t586 = Ifges(5,4) * t617 - Ifges(5,2) * t616 + Ifges(5,6) * qJD(3);
t531 = mrSges(5,2) * t580 - mrSges(5,3) * t566 + Ifges(5,1) * t596 + Ifges(5,4) * t595 + Ifges(5,5) * qJDD(3) - pkin(7) * t547 - qJD(3) * t586 - t651 * t550 + t654 * t551 - t616 * t585;
t587 = Ifges(5,1) * t617 - Ifges(5,4) * t616 + Ifges(5,5) * qJD(3);
t661 = mrSges(6,1) * t560 - mrSges(6,2) * t561 + Ifges(6,5) * t574 + Ifges(6,6) * t573 + Ifges(6,3) * t594 + t601 * t569 - t600 * t570;
t532 = -mrSges(5,1) * t580 + mrSges(5,3) * t567 + Ifges(5,4) * t596 + Ifges(5,2) * t595 + Ifges(5,6) * qJDD(3) - pkin(4) * t547 + qJD(3) * t587 - t617 * t585 - t661;
t618 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t655 - Ifges(4,6) * t652) * qJD(1);
t620 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t655 - Ifges(4,4) * t652) * qJD(1);
t522 = -mrSges(4,1) * t605 + mrSges(4,3) * t599 + Ifges(4,4) * t629 + Ifges(4,2) * t628 + Ifges(4,6) * qJDD(3) - pkin(3) * t545 + qJ(4) * t672 + qJD(3) * t620 + t649 * t531 + t650 * t532 - t618 * t678;
t619 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t655 - Ifges(4,2) * t652) * qJD(1);
t524 = mrSges(4,2) * t605 - mrSges(4,3) * t598 + Ifges(4,1) * t629 + Ifges(4,4) * t628 + Ifges(4,5) * qJDD(3) - qJ(4) * t537 - qJD(3) * t619 + t650 * t531 - t649 * t532 - t618 * t679;
t528 = qJDD(1) * mrSges(3,2) - t665;
t662 = mrSges(2,1) * t634 - mrSges(2,2) * t635 + mrSges(3,2) * t615 - mrSges(3,3) * t611 - pkin(1) * t528 - pkin(6) * t530 + qJ(2) * t660 - t652 * t522 + t655 * t524 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t659 = mrSges(5,1) * t566 - mrSges(4,2) * t599 - mrSges(5,2) * t567 + Ifges(5,5) * t596 + Ifges(5,6) * t595 + pkin(3) * t537 + pkin(4) * t664 + pkin(7) * t671 + t654 * t550 + t651 * t551 + t617 * t586 + t616 * t587 + mrSges(4,1) * t598 + t620 * t679 + t619 * t678 + Ifges(4,6) * t628 + Ifges(4,5) * t629 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3);
t529 = -m(3) * g(3) + t673;
t521 = t659 + (t681 * t658) - mrSges(2,3) * t634 + mrSges(3,1) * t615 + pkin(2) * t530 - qJ(2) * t529 + t682 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t520 = -mrSges(3,1) * t611 + mrSges(2,3) * t635 - pkin(1) * t529 - pkin(2) * t663 - pkin(6) * t673 + g(3) * t683 - qJDD(1) * t681 - t655 * t522 - t652 * t524 + t658 * t682;
t1 = [-m(1) * g(1) + t674; -m(1) * g(2) + t680; (-m(1) - m(2) - m(3)) * g(3) + t673; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t680 - t653 * t520 + t656 * t521; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t674 + t656 * t520 + t653 * t521; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t662; t662; t528; t659; t545; t661;];
tauJB = t1;
