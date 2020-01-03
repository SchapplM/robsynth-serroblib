% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR15
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR15_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR15_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR15_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:54
% EndTime: 2019-12-31 18:36:58
% DurationCPUTime: 2.90s
% Computational Cost: add. (30202->263), mult. (62389->322), div. (0->0), fcn. (37103->8), ass. (0->107)
t630 = sin(qJ(1));
t633 = cos(qJ(1));
t615 = -t633 * g(1) - t630 * g(2);
t660 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t615;
t659 = -pkin(1) - pkin(6);
t658 = mrSges(2,1) - mrSges(3,2);
t657 = -Ifges(3,4) + Ifges(2,5);
t656 = (Ifges(3,5) - Ifges(2,6));
t614 = t630 * g(1) - t633 * g(2);
t635 = qJD(1) ^ 2;
t643 = -t635 * qJ(2) + qJDD(2) - t614;
t586 = t659 * qJDD(1) + t643;
t629 = sin(qJ(3));
t632 = cos(qJ(3));
t576 = -t632 * g(3) + t629 * t586;
t608 = (mrSges(4,1) * t629 + mrSges(4,2) * t632) * qJD(1);
t652 = qJD(1) * qJD(3);
t649 = t632 * t652;
t609 = t629 * qJDD(1) + t649;
t654 = qJD(1) * t632;
t613 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t654;
t580 = t659 * t635 - t660;
t650 = t629 * t652;
t610 = t632 * qJDD(1) - t650;
t561 = (-t610 + t650) * qJ(4) + (t609 + t649) * pkin(3) + t580;
t607 = (pkin(3) * t629 - qJ(4) * t632) * qJD(1);
t634 = qJD(3) ^ 2;
t653 = t629 * qJD(1);
t564 = -t634 * pkin(3) + qJDD(3) * qJ(4) - t607 * t653 + t576;
t626 = sin(pkin(8));
t627 = cos(pkin(8));
t602 = t626 * qJD(3) + t627 * t654;
t548 = -0.2e1 * qJD(4) * t602 + t627 * t561 - t626 * t564;
t584 = t626 * qJDD(3) + t627 * t610;
t601 = t627 * qJD(3) - t626 * t654;
t546 = (t601 * t653 - t584) * pkin(7) + (t601 * t602 + t609) * pkin(4) + t548;
t549 = 0.2e1 * qJD(4) * t601 + t626 * t561 + t627 * t564;
t583 = t627 * qJDD(3) - t626 * t610;
t585 = pkin(4) * t653 - t602 * pkin(7);
t600 = t601 ^ 2;
t547 = -t600 * pkin(4) + t583 * pkin(7) - t585 * t653 + t549;
t628 = sin(qJ(5));
t631 = cos(qJ(5));
t544 = t631 * t546 - t628 * t547;
t571 = t631 * t601 - t628 * t602;
t553 = t571 * qJD(5) + t628 * t583 + t631 * t584;
t572 = t628 * t601 + t631 * t602;
t558 = -t571 * mrSges(6,1) + t572 * mrSges(6,2);
t616 = qJD(5) + t653;
t565 = -t616 * mrSges(6,2) + t571 * mrSges(6,3);
t606 = qJDD(5) + t609;
t540 = m(6) * t544 + t606 * mrSges(6,1) - t553 * mrSges(6,3) - t572 * t558 + t616 * t565;
t545 = t628 * t546 + t631 * t547;
t552 = -t572 * qJD(5) + t631 * t583 - t628 * t584;
t566 = t616 * mrSges(6,1) - t572 * mrSges(6,3);
t541 = m(6) * t545 - t606 * mrSges(6,2) + t552 * mrSges(6,3) + t571 * t558 - t616 * t566;
t533 = t631 * t540 + t628 * t541;
t573 = -t601 * mrSges(5,1) + t602 * mrSges(5,2);
t581 = -mrSges(5,2) * t653 + t601 * mrSges(5,3);
t531 = m(5) * t548 + t609 * mrSges(5,1) - t584 * mrSges(5,3) - t602 * t573 + t581 * t653 + t533;
t582 = mrSges(5,1) * t653 - t602 * mrSges(5,3);
t645 = -t628 * t540 + t631 * t541;
t532 = m(5) * t549 - t609 * mrSges(5,2) + t583 * mrSges(5,3) + t601 * t573 - t582 * t653 + t645;
t646 = -t626 * t531 + t627 * t532;
t525 = m(4) * t576 - qJDD(3) * mrSges(4,2) - t609 * mrSges(4,3) - qJD(3) * t613 - t608 * t653 + t646;
t575 = t629 * g(3) + t632 * t586;
t563 = -qJDD(3) * pkin(3) - t634 * qJ(4) + t607 * t654 + qJDD(4) - t575;
t550 = -t583 * pkin(4) - t600 * pkin(7) + t602 * t585 + t563;
t640 = m(6) * t550 - t552 * mrSges(6,1) + t553 * mrSges(6,2) - t571 * t565 + t572 * t566;
t543 = m(5) * t563 - t583 * mrSges(5,1) + t584 * mrSges(5,2) - t601 * t581 + t602 * t582 + t640;
t612 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t653;
t536 = m(4) * t575 + qJDD(3) * mrSges(4,1) - t610 * mrSges(4,3) + qJD(3) * t612 - t608 * t654 - t543;
t517 = t629 * t525 + t632 * t536;
t591 = -qJDD(1) * pkin(1) + t643;
t642 = -m(3) * t591 + (t635 * mrSges(3,3)) - t517;
t511 = m(2) * t614 - (t635 * mrSges(2,2)) + t658 * qJDD(1) + t642;
t589 = t635 * pkin(1) + t660;
t527 = t627 * t531 + t626 * t532;
t641 = -m(4) * t580 - t609 * mrSges(4,1) - t610 * mrSges(4,2) - t612 * t653 - t613 * t654 - t527;
t638 = -m(3) * t589 + (t635 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t641;
t522 = m(2) * t615 - (t635 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t638;
t655 = t633 * t511 + t630 * t522;
t648 = -t630 * t511 + t633 * t522;
t647 = t632 * t525 - t629 * t536;
t554 = Ifges(6,5) * t572 + Ifges(6,6) * t571 + Ifges(6,3) * t616;
t556 = Ifges(6,1) * t572 + Ifges(6,4) * t571 + Ifges(6,5) * t616;
t534 = -mrSges(6,1) * t550 + mrSges(6,3) * t545 + Ifges(6,4) * t553 + Ifges(6,2) * t552 + Ifges(6,6) * t606 - t572 * t554 + t616 * t556;
t555 = Ifges(6,4) * t572 + Ifges(6,2) * t571 + Ifges(6,6) * t616;
t535 = mrSges(6,2) * t550 - mrSges(6,3) * t544 + Ifges(6,1) * t553 + Ifges(6,4) * t552 + Ifges(6,5) * t606 + t571 * t554 - t616 * t555;
t567 = Ifges(5,5) * t602 + Ifges(5,6) * t601 + Ifges(5,3) * t653;
t569 = Ifges(5,1) * t602 + Ifges(5,4) * t601 + Ifges(5,5) * t653;
t513 = -mrSges(5,1) * t563 + mrSges(5,3) * t549 + Ifges(5,4) * t584 + Ifges(5,2) * t583 + Ifges(5,6) * t609 - pkin(4) * t640 + pkin(7) * t645 + t631 * t534 + t628 * t535 - t602 * t567 + t569 * t653;
t568 = Ifges(5,4) * t602 + Ifges(5,2) * t601 + Ifges(5,6) * t653;
t519 = mrSges(5,2) * t563 - mrSges(5,3) * t548 + Ifges(5,1) * t584 + Ifges(5,4) * t583 + Ifges(5,5) * t609 - pkin(7) * t533 - t628 * t534 + t631 * t535 + t601 * t567 - t568 * t653;
t596 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t632 - Ifges(4,2) * t629) * qJD(1);
t597 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t632 - Ifges(4,4) * t629) * qJD(1);
t639 = mrSges(4,1) * t575 - mrSges(4,2) * t576 + Ifges(4,5) * t610 - Ifges(4,6) * t609 + Ifges(4,3) * qJDD(3) - pkin(3) * t543 + qJ(4) * t646 + t627 * t513 + t626 * t519 + t596 * t654 + t597 * t653;
t595 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t632 - Ifges(4,6) * t629) * qJD(1);
t508 = mrSges(4,2) * t580 - mrSges(4,3) * t575 + Ifges(4,1) * t610 - Ifges(4,4) * t609 + Ifges(4,5) * qJDD(3) - qJ(4) * t527 - qJD(3) * t596 - t626 * t513 + t627 * t519 - t595 * t653;
t636 = mrSges(6,1) * t544 - mrSges(6,2) * t545 + Ifges(6,5) * t553 + Ifges(6,6) * t552 + Ifges(6,3) * t606 + t572 * t555 - t571 * t556;
t509 = qJD(3) * t597 + t601 * t569 + mrSges(4,3) * t576 - mrSges(4,1) * t580 - Ifges(5,6) * t583 - Ifges(5,5) * t584 - mrSges(5,1) * t548 + mrSges(5,2) * t549 - t636 + (-Ifges(4,2) - Ifges(5,3)) * t609 - t602 * t568 + Ifges(4,4) * t610 - pkin(3) * t527 + Ifges(4,6) * qJDD(3) - pkin(4) * t533 - t595 * t654;
t515 = qJDD(1) * mrSges(3,2) - t642;
t637 = mrSges(2,1) * t614 - mrSges(2,2) * t615 + mrSges(3,2) * t591 - mrSges(3,3) * t589 - pkin(1) * t515 - pkin(6) * t517 + qJ(2) * t638 + t632 * t508 - t629 * t509 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t516 = -m(3) * g(3) + t647;
t506 = mrSges(3,1) * t591 + t639 - mrSges(2,3) * t614 + (t656 * t635) + pkin(2) * t517 - qJ(2) * t516 + t657 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t505 = -mrSges(3,1) * t589 + mrSges(2,3) * t615 - pkin(1) * t516 - pkin(2) * t641 - pkin(6) * t647 + t658 * g(3) - t656 * qJDD(1) - t629 * t508 - t632 * t509 + t657 * t635;
t1 = [-m(1) * g(1) + t648; -m(1) * g(2) + t655; (-m(1) - m(2) - m(3)) * g(3) + t647; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t655 - t630 * t505 + t633 * t506; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t648 + t633 * t505 + t630 * t506; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t637; t637; t515; t639; t543; t636;];
tauJB = t1;
