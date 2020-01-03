% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP13_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP13_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:43
% EndTime: 2019-12-31 18:58:46
% DurationCPUTime: 1.74s
% Computational Cost: add. (13922->240), mult. (26072->279), div. (0->0), fcn. (14063->6), ass. (0->100)
t675 = Ifges(5,1) + Ifges(6,1);
t665 = Ifges(5,4) - Ifges(6,5);
t663 = -Ifges(5,5) - Ifges(6,4);
t674 = Ifges(5,2) + Ifges(6,3);
t662 = Ifges(5,6) - Ifges(6,6);
t673 = -Ifges(5,3) - Ifges(6,2);
t631 = sin(qJ(1));
t633 = cos(qJ(1));
t615 = -t633 * g(1) - t631 * g(2);
t672 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t615;
t629 = sin(qJ(4));
t632 = cos(qJ(3));
t654 = qJD(1) * t632;
t669 = cos(qJ(4));
t605 = -t669 * qJD(3) + t629 * t654;
t630 = sin(qJ(3));
t653 = qJD(1) * qJD(3);
t650 = t630 * t653;
t610 = qJDD(1) * t632 - t650;
t571 = -t605 * qJD(4) + t629 * qJDD(3) + t669 * t610;
t606 = t629 * qJD(3) + t669 * t654;
t576 = mrSges(6,1) * t605 - mrSges(6,3) * t606;
t635 = qJD(1) ^ 2;
t670 = (-pkin(1) - pkin(6));
t586 = (t670 * t635) - t672;
t649 = t632 * t653;
t609 = -qJDD(1) * t630 - t649;
t556 = (-t610 + t650) * pkin(7) + (-t609 + t649) * pkin(3) + t586;
t614 = t631 * g(1) - t633 * g(2);
t643 = -t635 * qJ(2) + qJDD(2) - t614;
t587 = t670 * qJDD(1) + t643;
t579 = -g(3) * t632 + t630 * t587;
t608 = (pkin(3) * t630 - pkin(7) * t632) * qJD(1);
t634 = qJD(3) ^ 2;
t655 = qJD(1) * t630;
t559 = -pkin(3) * t634 + qJDD(3) * pkin(7) - t608 * t655 + t579;
t553 = t669 * t556 - t629 * t559;
t575 = pkin(4) * t605 - qJ(5) * t606;
t604 = qJDD(4) - t609;
t617 = qJD(4) + t655;
t616 = t617 ^ 2;
t551 = -t604 * pkin(4) - t616 * qJ(5) + t606 * t575 + qJDD(5) - t553;
t583 = -mrSges(6,2) * t605 + mrSges(6,3) * t617;
t645 = -m(6) * t551 + t604 * mrSges(6,1) + t617 * t583;
t548 = t571 * mrSges(6,2) + t606 * t576 - t645;
t554 = t629 * t556 + t669 * t559;
t550 = -pkin(4) * t616 + qJ(5) * t604 + 0.2e1 * qJD(5) * t617 - t575 * t605 + t554;
t570 = qJD(4) * t606 - t669 * qJDD(3) + t610 * t629;
t582 = -mrSges(6,1) * t617 + mrSges(6,2) * t606;
t651 = m(6) * t550 + t604 * mrSges(6,3) + t617 * t582;
t657 = t665 * t605 - t675 * t606 + t663 * t617;
t658 = t674 * t605 - t665 * t606 - t662 * t617;
t671 = -t662 * t570 - t663 * t571 - t673 * t604 - t657 * t605 - t658 * t606 + mrSges(5,1) * t553 - mrSges(6,1) * t551 - mrSges(5,2) * t554 + mrSges(6,3) * t550 - pkin(4) * t548 + qJ(5) * (-t570 * mrSges(6,2) - t605 * t576 + t651);
t668 = mrSges(2,1) - mrSges(3,2);
t667 = -mrSges(5,3) - mrSges(6,2);
t666 = -Ifges(3,4) + Ifges(2,5);
t664 = (Ifges(3,5) - Ifges(2,6));
t607 = (mrSges(4,1) * t630 + mrSges(4,2) * t632) * qJD(1);
t613 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t654;
t581 = mrSges(5,1) * t617 - mrSges(5,3) * t606;
t656 = -mrSges(5,1) * t605 - mrSges(5,2) * t606 - t576;
t542 = m(5) * t554 - t604 * mrSges(5,2) + t667 * t570 - t617 * t581 + t656 * t605 + t651;
t580 = -mrSges(5,2) * t617 - mrSges(5,3) * t605;
t544 = m(5) * t553 + t604 * mrSges(5,1) + t667 * t571 + t617 * t580 + t656 * t606 + t645;
t646 = t669 * t542 - t544 * t629;
t536 = m(4) * t579 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t609 - qJD(3) * t613 - t607 * t655 + t646;
t578 = t630 * g(3) + t632 * t587;
t612 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t655;
t558 = -qJDD(3) * pkin(3) - t634 * pkin(7) + t608 * t654 - t578;
t552 = -0.2e1 * qJD(5) * t606 + (t605 * t617 - t571) * qJ(5) + (t606 * t617 + t570) * pkin(4) + t558;
t546 = m(6) * t552 + t570 * mrSges(6,1) - t571 * mrSges(6,3) - t606 * t582 + t583 * t605;
t636 = -m(5) * t558 - t570 * mrSges(5,1) - t571 * mrSges(5,2) - t605 * t580 - t581 * t606 - t546;
t539 = m(4) * t578 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t610 + qJD(3) * t612 - t607 * t654 + t636;
t526 = t630 * t536 + t632 * t539;
t592 = -qJDD(1) * pkin(1) + t643;
t642 = -m(3) * t592 + (t635 * mrSges(3,3)) - t526;
t522 = m(2) * t614 - (t635 * mrSges(2,2)) + t668 * qJDD(1) + t642;
t590 = t635 * pkin(1) + t672;
t538 = t629 * t542 + t669 * t544;
t641 = -m(4) * t586 + mrSges(4,1) * t609 - t610 * mrSges(4,2) - t612 * t655 - t613 * t654 - t538;
t639 = -m(3) * t590 + (t635 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t641;
t529 = m(2) * t615 - (mrSges(2,1) * t635) - qJDD(1) * mrSges(2,2) + t639;
t660 = t633 * t522 + t631 * t529;
t659 = t662 * t605 + t663 * t606 + t673 * t617;
t648 = -t522 * t631 + t633 * t529;
t647 = t632 * t536 - t630 * t539;
t531 = -mrSges(5,1) * t558 - mrSges(6,1) * t552 + mrSges(6,2) * t550 + mrSges(5,3) * t554 - pkin(4) * t546 - t674 * t570 + t665 * t571 + t662 * t604 + t659 * t606 - t657 * t617;
t533 = mrSges(5,2) * t558 + mrSges(6,2) * t551 - mrSges(5,3) * t553 - mrSges(6,3) * t552 - qJ(5) * t546 - t665 * t570 + t675 * t571 - t663 * t604 + t659 * t605 + t658 * t617;
t595 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t632 - Ifges(4,2) * t630) * qJD(1);
t596 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t632 - Ifges(4,4) * t630) * qJD(1);
t640 = mrSges(4,1) * t578 - mrSges(4,2) * t579 + Ifges(4,5) * t610 + Ifges(4,6) * t609 + Ifges(4,3) * qJDD(3) + pkin(3) * t636 + pkin(7) * t646 + t669 * t531 + t629 * t533 + t595 * t654 + t596 * t655;
t594 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t632 - Ifges(4,6) * t630) * qJD(1);
t519 = mrSges(4,2) * t586 - mrSges(4,3) * t578 + Ifges(4,1) * t610 + Ifges(4,4) * t609 + Ifges(4,5) * qJDD(3) - pkin(7) * t538 - qJD(3) * t595 - t629 * t531 + t669 * t533 - t594 * t655;
t520 = -mrSges(4,1) * t586 + mrSges(4,3) * t579 + Ifges(4,4) * t610 + Ifges(4,2) * t609 + Ifges(4,6) * qJDD(3) - pkin(3) * t538 + qJD(3) * t596 - t594 * t654 - t671;
t524 = qJDD(1) * mrSges(3,2) - t642;
t637 = mrSges(2,1) * t614 - mrSges(2,2) * t615 + mrSges(3,2) * t592 - mrSges(3,3) * t590 - pkin(1) * t524 - pkin(6) * t526 + qJ(2) * t639 + t632 * t519 - t520 * t630 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t525 = -m(3) * g(3) + t647;
t517 = t640 + (t664 * t635) + mrSges(3,1) * t592 - mrSges(2,3) * t614 + t666 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - qJ(2) * t525 + pkin(2) * t526;
t516 = -mrSges(3,1) * t590 + mrSges(2,3) * t615 - pkin(1) * t525 - pkin(2) * t641 - pkin(6) * t647 + t668 * g(3) - t664 * qJDD(1) - t630 * t519 - t632 * t520 + t666 * t635;
t1 = [-m(1) * g(1) + t648; -m(1) * g(2) + t660; (-m(1) - m(2) - m(3)) * g(3) + t647; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t660 - t631 * t516 + t633 * t517; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t648 + t633 * t516 + t631 * t517; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t637; t637; t524; t640; t671; t548;];
tauJB = t1;
