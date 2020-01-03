% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:26
% EndTime: 2019-12-31 18:14:30
% DurationCPUTime: 1.97s
% Computational Cost: add. (13783->241), mult. (29167->284), div. (0->0), fcn. (16131->6), ass. (0->98)
t687 = Ifges(5,1) + Ifges(6,1);
t676 = Ifges(5,4) - Ifges(6,5);
t685 = Ifges(5,5) + Ifges(6,4);
t686 = Ifges(5,2) + Ifges(6,3);
t683 = Ifges(5,6) - Ifges(6,6);
t684 = -Ifges(6,2) - Ifges(5,3);
t644 = sin(pkin(7));
t645 = cos(pkin(7));
t646 = sin(qJ(3));
t648 = cos(qJ(3));
t607 = (t644 * t648 + t645 * t646) * qJD(1);
t669 = qJD(1) * t648;
t670 = qJD(1) * t646;
t608 = -t644 * t670 + t645 * t669;
t682 = -t683 * qJD(3) + t686 * t607 - t676 * t608;
t681 = t685 * qJD(3) - t676 * t607 + t687 * t608;
t647 = sin(qJ(1));
t649 = cos(qJ(1));
t627 = -t649 * g(1) - t647 * g(2);
t658 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t627;
t680 = -2 * qJD(4);
t679 = -pkin(1) - pkin(6);
t678 = mrSges(2,1) - mrSges(3,2);
t677 = -mrSges(5,3) - mrSges(6,2);
t675 = Ifges(2,5) - Ifges(3,4);
t674 = -Ifges(2,6) + Ifges(3,5);
t626 = t647 * g(1) - t649 * g(2);
t651 = qJD(1) ^ 2;
t657 = -t651 * qJ(2) + qJDD(2) - t626;
t600 = t679 * qJDD(1) + t657;
t588 = t646 * g(3) + t648 * t600;
t668 = qJD(1) * qJD(3);
t665 = t646 * t668;
t621 = t648 * qJDD(1) - t665;
t562 = (-t621 - t665) * qJ(4) + (-t646 * t648 * t651 + qJDD(3)) * pkin(3) + t588;
t589 = -t648 * g(3) + t646 * t600;
t620 = -t646 * qJDD(1) - t648 * t668;
t624 = qJD(3) * pkin(3) - qJ(4) * t669;
t641 = t646 ^ 2;
t563 = -t641 * t651 * pkin(3) + t620 * qJ(4) - qJD(3) * t624 + t589;
t559 = t644 * t562 + t645 * t563 + t607 * t680;
t585 = -t645 * t620 + t644 * t621;
t597 = qJD(3) * mrSges(5,1) - t608 * mrSges(5,3);
t576 = t607 * pkin(4) - t608 * qJ(5);
t650 = qJD(3) ^ 2;
t553 = -t650 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t607 * t576 + t559;
t598 = -qJD(3) * mrSges(6,1) + t608 * mrSges(6,2);
t666 = m(6) * t553 + qJDD(3) * mrSges(6,3) + qJD(3) * t598;
t577 = t607 * mrSges(6,1) - t608 * mrSges(6,3);
t671 = -t607 * mrSges(5,1) - t608 * mrSges(5,2) - t577;
t544 = m(5) * t559 - qJDD(3) * mrSges(5,2) - qJD(3) * t597 + t677 * t585 + t671 * t607 + t666;
t660 = -t645 * t562 + t644 * t563;
t558 = t608 * t680 - t660;
t586 = t644 * t620 + t645 * t621;
t596 = -qJD(3) * mrSges(5,2) - t607 * mrSges(5,3);
t554 = -qJDD(3) * pkin(4) - t650 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t576) * t608 + t660;
t599 = -t607 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t661 = -m(6) * t554 + qJDD(3) * mrSges(6,1) + qJD(3) * t599;
t545 = m(5) * t558 + qJDD(3) * mrSges(5,1) + qJD(3) * t596 + t677 * t586 + t671 * t608 + t661;
t536 = t644 * t544 + t645 * t545;
t619 = (mrSges(4,1) * t646 + mrSges(4,2) * t648) * qJD(1);
t623 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t670;
t533 = m(4) * t588 + qJDD(3) * mrSges(4,1) - t621 * mrSges(4,3) + qJD(3) * t623 - t619 * t669 + t536;
t625 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t669;
t662 = t645 * t544 - t644 * t545;
t534 = m(4) * t589 - qJDD(3) * mrSges(4,2) + t620 * mrSges(4,3) - qJD(3) * t625 - t619 * t670 + t662;
t529 = t648 * t533 + t646 * t534;
t606 = -qJDD(1) * pkin(1) + t657;
t656 = -m(3) * t606 + t651 * mrSges(3,3) - t529;
t525 = m(2) * t626 - t651 * mrSges(2,2) + t678 * qJDD(1) + t656;
t603 = t651 * pkin(1) - t658;
t565 = -t620 * pkin(3) + qJDD(4) + t624 * t669 + (-qJ(4) * t641 + t679) * t651 + t658;
t556 = -0.2e1 * qJD(5) * t608 + (qJD(3) * t607 - t586) * qJ(5) + (qJD(3) * t608 + t585) * pkin(4) + t565;
t548 = m(6) * t556 + t585 * mrSges(6,1) - t586 * mrSges(6,3) - t608 * t598 + t607 * t599;
t546 = m(5) * t565 + t585 * mrSges(5,1) + t586 * mrSges(5,2) + t607 * t596 + t608 * t597 + t548;
t595 = t679 * t651 + t658;
t654 = -m(4) * t595 + t620 * mrSges(4,1) - t621 * mrSges(4,2) - t623 * t670 - t625 * t669 - t546;
t653 = -m(3) * t603 + t651 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t654;
t539 = m(2) * t627 - t651 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t653;
t673 = t649 * t525 + t647 * t539;
t672 = t684 * qJD(3) + t683 * t607 - t685 * t608;
t664 = -t647 * t525 + t649 * t539;
t663 = -t646 * t533 + t648 * t534;
t530 = -mrSges(5,1) * t565 - mrSges(6,1) * t556 + mrSges(6,2) * t553 + mrSges(5,3) * t559 - pkin(4) * t548 + t681 * qJD(3) + t683 * qJDD(3) - t686 * t585 + t676 * t586 + t672 * t608;
t531 = mrSges(5,2) * t565 + mrSges(6,2) * t554 - mrSges(5,3) * t558 - mrSges(6,3) * t556 - qJ(5) * t548 + t682 * qJD(3) + t685 * qJDD(3) - t676 * t585 + t687 * t586 + t672 * t607;
t609 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t648 - Ifges(4,6) * t646) * qJD(1);
t611 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t648 - Ifges(4,4) * t646) * qJD(1);
t521 = -mrSges(4,1) * t595 + mrSges(4,3) * t589 + Ifges(4,4) * t621 + Ifges(4,2) * t620 + Ifges(4,6) * qJDD(3) - pkin(3) * t546 + qJ(4) * t662 + qJD(3) * t611 + t645 * t530 + t644 * t531 - t609 * t669;
t610 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t648 - Ifges(4,2) * t646) * qJD(1);
t523 = mrSges(4,2) * t595 - mrSges(4,3) * t588 + Ifges(4,1) * t621 + Ifges(4,4) * t620 + Ifges(4,5) * qJDD(3) - qJ(4) * t536 - qJD(3) * t610 - t644 * t530 + t645 * t531 - t609 * t670;
t527 = qJDD(1) * mrSges(3,2) - t656;
t655 = mrSges(2,1) * t626 - mrSges(2,2) * t627 + mrSges(3,2) * t606 - mrSges(3,3) * t603 - pkin(1) * t527 - pkin(6) * t529 + qJ(2) * t653 - t646 * t521 + t648 * t523 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t550 = t586 * mrSges(6,2) + t608 * t577 - t661;
t652 = mrSges(4,1) * t588 + mrSges(5,1) * t558 - mrSges(6,1) * t554 - mrSges(4,2) * t589 - mrSges(5,2) * t559 + mrSges(6,3) * t553 + Ifges(4,5) * t621 + Ifges(4,6) * t620 + pkin(3) * t536 - pkin(4) * t550 + qJ(5) * t666 + t610 * t669 + t611 * t670 - t682 * t608 + (-qJ(5) * t577 + t681) * t607 + t685 * t586 + (-qJ(5) * mrSges(6,2) - t683) * t585 + (Ifges(4,3) - t684) * qJDD(3);
t528 = -m(3) * g(3) + t663;
t520 = -qJ(2) * t528 + t674 * t651 - mrSges(2,3) * t626 + pkin(2) * t529 + t675 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t652 + mrSges(3,1) * t606;
t519 = -mrSges(3,1) * t603 + mrSges(2,3) * t627 - pkin(1) * t528 - pkin(2) * t654 - pkin(6) * t663 + t678 * g(3) - t674 * qJDD(1) - t648 * t521 - t646 * t523 + t675 * t651;
t1 = [-m(1) * g(1) + t664; -m(1) * g(2) + t673; (-m(1) - m(2) - m(3)) * g(3) + t663; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t673 - t647 * t519 + t649 * t520; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t664 + t649 * t519 + t647 * t520; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t655; t655; t527; t652; t546; t550;];
tauJB = t1;
