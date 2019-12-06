% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRR5
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:56
% EndTime: 2019-12-05 15:54:04
% DurationCPUTime: 4.41s
% Computational Cost: add. (47276->239), mult. (107414->302), div. (0->0), fcn. (76474->10), ass. (0->110)
t641 = qJD(2) ^ 2;
t634 = cos(pkin(9));
t669 = pkin(3) * t634;
t632 = sin(pkin(9));
t668 = mrSges(4,2) * t632;
t667 = cos(pkin(8));
t629 = t634 ^ 2;
t666 = t629 * t641;
t633 = sin(pkin(8));
t617 = -t667 * g(1) - t633 * g(2);
t631 = -g(3) + qJDD(1);
t637 = sin(qJ(2));
t640 = cos(qJ(2));
t606 = t640 * t617 + t637 * t631;
t601 = -t641 * pkin(2) + qJDD(2) * qJ(3) + t606;
t616 = t633 * g(1) - t667 * g(2);
t661 = qJD(2) * qJD(3);
t664 = -t634 * t616 - 0.2e1 * t632 * t661;
t582 = (-pkin(6) * qJDD(2) + t641 * t669 - t601) * t632 + t664;
t585 = -t632 * t616 + (t601 + 0.2e1 * t661) * t634;
t660 = qJDD(2) * t634;
t583 = -pkin(3) * t666 + pkin(6) * t660 + t585;
t636 = sin(qJ(4));
t639 = cos(qJ(4));
t566 = t639 * t582 - t636 * t583;
t650 = t632 * t639 + t634 * t636;
t649 = -t632 * t636 + t634 * t639;
t608 = t649 * qJD(2);
t662 = t608 * qJD(4);
t598 = t650 * qJDD(2) + t662;
t609 = t650 * qJD(2);
t563 = (-t598 + t662) * pkin(7) + (t608 * t609 + qJDD(4)) * pkin(4) + t566;
t567 = t636 * t582 + t639 * t583;
t597 = -t609 * qJD(4) + t649 * qJDD(2);
t604 = qJD(4) * pkin(4) - t609 * pkin(7);
t607 = t608 ^ 2;
t564 = -t607 * pkin(4) + t597 * pkin(7) - qJD(4) * t604 + t567;
t635 = sin(qJ(5));
t638 = cos(qJ(5));
t561 = t638 * t563 - t635 * t564;
t592 = t638 * t608 - t635 * t609;
t573 = t592 * qJD(5) + t635 * t597 + t638 * t598;
t593 = t635 * t608 + t638 * t609;
t578 = -t592 * mrSges(6,1) + t593 * mrSges(6,2);
t630 = qJD(4) + qJD(5);
t586 = -t630 * mrSges(6,2) + t592 * mrSges(6,3);
t627 = qJDD(4) + qJDD(5);
t558 = m(6) * t561 + t627 * mrSges(6,1) - t573 * mrSges(6,3) - t593 * t578 + t630 * t586;
t562 = t635 * t563 + t638 * t564;
t572 = -t593 * qJD(5) + t638 * t597 - t635 * t598;
t587 = t630 * mrSges(6,1) - t593 * mrSges(6,3);
t559 = m(6) * t562 - t627 * mrSges(6,2) + t572 * mrSges(6,3) + t592 * t578 - t630 * t587;
t549 = t638 * t558 + t635 * t559;
t595 = -t608 * mrSges(5,1) + t609 * mrSges(5,2);
t602 = -qJD(4) * mrSges(5,2) + t608 * mrSges(5,3);
t547 = m(5) * t566 + qJDD(4) * mrSges(5,1) - t598 * mrSges(5,3) + qJD(4) * t602 - t609 * t595 + t549;
t603 = qJD(4) * mrSges(5,1) - t609 * mrSges(5,3);
t655 = -t635 * t558 + t638 * t559;
t548 = m(5) * t567 - qJDD(4) * mrSges(5,2) + t597 * mrSges(5,3) - qJD(4) * t603 + t608 * t595 + t655;
t543 = t639 * t547 + t636 * t548;
t584 = -t632 * t601 + t664;
t648 = mrSges(4,3) * qJDD(2) + t641 * (-mrSges(4,1) * t634 + t668);
t541 = m(4) * t584 - t648 * t632 + t543;
t656 = -t636 * t547 + t639 * t548;
t542 = m(4) * t585 + t648 * t634 + t656;
t537 = -t632 * t541 + t634 * t542;
t533 = m(3) * t606 - t641 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t537;
t605 = -t637 * t617 + t640 * t631;
t651 = qJDD(3) - t605;
t600 = -qJDD(2) * pkin(2) - t641 * qJ(3) + t651;
t628 = t632 ^ 2;
t588 = (-pkin(2) - t669) * qJDD(2) + (-qJ(3) + (-t628 - t629) * pkin(6)) * t641 + t651;
t565 = -t597 * pkin(4) - t607 * pkin(7) + t609 * t604 + t588;
t647 = m(6) * t565 - t572 * mrSges(6,1) + t573 * mrSges(6,2) - t592 * t586 + t593 * t587;
t644 = m(5) * t588 - t597 * mrSges(5,1) + t598 * mrSges(5,2) - t608 * t602 + t609 * t603 + t647;
t643 = -m(4) * t600 + mrSges(4,1) * t660 - t644 + (t628 * t641 + t666) * mrSges(4,3);
t553 = (mrSges(3,1) - t668) * qJDD(2) + t643 - t641 * mrSges(3,2) + m(3) * t605;
t657 = t640 * t533 - t637 * t553;
t529 = m(2) * t617 + t657;
t536 = t634 * t541 + t632 * t542;
t535 = (m(2) + m(3)) * t616 - t536;
t665 = t633 * t529 + t667 * t535;
t530 = t637 * t533 + t640 * t553;
t652 = Ifges(4,5) * t632 + Ifges(4,6) * t634;
t663 = t641 * t652;
t659 = m(2) * t631 + t530;
t658 = t667 * t529 - t633 * t535;
t654 = Ifges(4,1) * t632 + Ifges(4,4) * t634;
t653 = Ifges(4,4) * t632 + Ifges(4,2) * t634;
t575 = Ifges(6,4) * t593 + Ifges(6,2) * t592 + Ifges(6,6) * t630;
t576 = Ifges(6,1) * t593 + Ifges(6,4) * t592 + Ifges(6,5) * t630;
t646 = -mrSges(6,1) * t561 + mrSges(6,2) * t562 - Ifges(6,5) * t573 - Ifges(6,6) * t572 - Ifges(6,3) * t627 - t593 * t575 + t592 * t576;
t574 = Ifges(6,5) * t593 + Ifges(6,6) * t592 + Ifges(6,3) * t630;
t550 = -mrSges(6,1) * t565 + mrSges(6,3) * t562 + Ifges(6,4) * t573 + Ifges(6,2) * t572 + Ifges(6,6) * t627 - t593 * t574 + t630 * t576;
t551 = mrSges(6,2) * t565 - mrSges(6,3) * t561 + Ifges(6,1) * t573 + Ifges(6,4) * t572 + Ifges(6,5) * t627 + t592 * t574 - t630 * t575;
t589 = Ifges(5,5) * t609 + Ifges(5,6) * t608 + Ifges(5,3) * qJD(4);
t591 = Ifges(5,1) * t609 + Ifges(5,4) * t608 + Ifges(5,5) * qJD(4);
t538 = -mrSges(5,1) * t588 + mrSges(5,3) * t567 + Ifges(5,4) * t598 + Ifges(5,2) * t597 + Ifges(5,6) * qJDD(4) - pkin(4) * t647 + pkin(7) * t655 + qJD(4) * t591 + t638 * t550 + t635 * t551 - t609 * t589;
t590 = Ifges(5,4) * t609 + Ifges(5,2) * t608 + Ifges(5,6) * qJD(4);
t539 = mrSges(5,2) * t588 - mrSges(5,3) * t566 + Ifges(5,1) * t598 + Ifges(5,4) * t597 + Ifges(5,5) * qJDD(4) - pkin(7) * t549 - qJD(4) * t590 - t635 * t550 + t638 * t551 + t608 * t589;
t525 = -mrSges(4,1) * t600 + mrSges(4,3) * t585 - pkin(3) * t644 + pkin(6) * t656 + t653 * qJDD(2) + t639 * t538 + t636 * t539 - t632 * t663;
t526 = mrSges(4,2) * t600 - mrSges(4,3) * t584 - pkin(6) * t543 + t654 * qJDD(2) - t636 * t538 + t639 * t539 + t634 * t663;
t554 = qJDD(2) * t668 - t643;
t645 = mrSges(3,1) * t605 - mrSges(3,2) * t606 + Ifges(3,3) * qJDD(2) - pkin(2) * t554 + qJ(3) * t537 + t634 * t525 + t632 * t526;
t642 = mrSges(5,1) * t566 - mrSges(5,2) * t567 + Ifges(5,5) * t598 + Ifges(5,6) * t597 + Ifges(5,3) * qJDD(4) + pkin(4) * t549 + t609 * t590 - t608 * t591 - t646;
t524 = -t642 + (Ifges(3,6) - t652) * qJDD(2) + mrSges(3,1) * t616 + mrSges(3,3) * t606 - mrSges(4,1) * t584 + mrSges(4,2) * t585 - pkin(3) * t543 - pkin(2) * t536 + (-t632 * t653 + t634 * t654 + Ifges(3,5)) * t641;
t523 = -mrSges(3,2) * t616 - mrSges(3,3) * t605 + Ifges(3,5) * qJDD(2) - t641 * Ifges(3,6) - qJ(3) * t536 - t632 * t525 + t634 * t526;
t522 = -mrSges(2,1) * t631 + mrSges(2,3) * t617 - pkin(1) * t530 - t645;
t521 = mrSges(2,2) * t631 - mrSges(2,3) * t616 - pkin(5) * t530 + t640 * t523 - t637 * t524;
t1 = [-m(1) * g(1) + t658; -m(1) * g(2) + t665; -m(1) * g(3) + t659; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t665 + t667 * t521 - t633 * t522; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t658 + t633 * t521 + t667 * t522; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t616 - mrSges(2,2) * t617 + t637 * t523 + t640 * t524 + pkin(1) * (m(3) * t616 - t536) + pkin(5) * t657; t659; t645; t554; t642; -t646;];
tauJB = t1;
