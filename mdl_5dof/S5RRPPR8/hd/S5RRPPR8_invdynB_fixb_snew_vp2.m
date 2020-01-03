% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:25
% EndTime: 2019-12-31 19:38:30
% DurationCPUTime: 3.36s
% Computational Cost: add. (29330->289), mult. (66265->357), div. (0->0), fcn. (38983->8), ass. (0->112)
t646 = Ifges(3,1) + Ifges(4,1);
t641 = Ifges(3,4) - Ifges(4,5);
t640 = Ifges(3,5) + Ifges(4,4);
t645 = Ifges(3,2) + Ifges(4,3);
t639 = Ifges(3,6) - Ifges(4,6);
t644 = Ifges(3,3) + Ifges(4,2);
t643 = 2 * qJD(3);
t642 = mrSges(3,3) + mrSges(4,2);
t615 = cos(qJ(2));
t618 = qJD(1) ^ 2;
t638 = t615 ^ 2 * t618;
t613 = sin(qJ(1));
t616 = cos(qJ(1));
t593 = -t616 * g(1) - t613 * g(2);
t571 = -t618 * pkin(1) + qJDD(1) * pkin(6) + t593;
t612 = sin(qJ(2));
t550 = -t612 * g(3) + t615 * t571;
t582 = (-mrSges(3,1) * t615 + mrSges(3,2) * t612) * qJD(1);
t631 = qJD(1) * qJD(2);
t630 = t612 * t631;
t584 = t615 * qJDD(1) - t630;
t633 = qJD(1) * t612;
t588 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t633;
t580 = (-pkin(2) * t615 - qJ(3) * t612) * qJD(1);
t617 = qJD(2) ^ 2;
t632 = qJD(1) * t615;
t533 = -t617 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t643 + t580 * t632 + t550;
t581 = (-mrSges(4,1) * t615 - mrSges(4,3) * t612) * qJD(1);
t589 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t633;
t587 = -qJD(2) * pkin(3) - qJ(4) * t633;
t529 = -pkin(3) * t638 - t584 * qJ(4) + qJD(2) * t587 + t533;
t549 = -t615 * g(3) - t612 * t571;
t534 = -qJDD(2) * pkin(2) - t617 * qJ(3) + t580 * t633 + qJDD(3) - t549;
t629 = t615 * t631;
t583 = t612 * qJDD(1) + t629;
t530 = (-t583 + t629) * qJ(4) + (-t612 * t615 * t618 - qJDD(2)) * pkin(3) + t534;
t608 = sin(pkin(8));
t609 = cos(pkin(8));
t563 = (-t608 * t615 + t609 * t612) * qJD(1);
t510 = -0.2e1 * qJD(4) * t563 - t608 * t529 + t609 * t530;
t548 = t609 * t583 - t608 * t584;
t562 = (-t608 * t612 - t609 * t615) * qJD(1);
t508 = (-qJD(2) * t562 - t548) * pkin(7) + (t562 * t563 - qJDD(2)) * pkin(4) + t510;
t511 = 0.2e1 * qJD(4) * t562 + t609 * t529 + t608 * t530;
t547 = -t608 * t583 - t609 * t584;
t553 = -qJD(2) * pkin(4) - t563 * pkin(7);
t561 = t562 ^ 2;
t509 = -t561 * pkin(4) + t547 * pkin(7) + qJD(2) * t553 + t511;
t611 = sin(qJ(5));
t614 = cos(qJ(5));
t506 = t614 * t508 - t611 * t509;
t540 = t614 * t562 - t611 * t563;
t517 = t540 * qJD(5) + t611 * t547 + t614 * t548;
t541 = t611 * t562 + t614 * t563;
t528 = -t540 * mrSges(6,1) + t541 * mrSges(6,2);
t602 = -qJD(2) + qJD(5);
t535 = -t602 * mrSges(6,2) + t540 * mrSges(6,3);
t601 = -qJDD(2) + qJDD(5);
t504 = m(6) * t506 + t601 * mrSges(6,1) - t517 * mrSges(6,3) - t541 * t528 + t602 * t535;
t507 = t611 * t508 + t614 * t509;
t516 = -t541 * qJD(5) + t614 * t547 - t611 * t548;
t536 = t602 * mrSges(6,1) - t541 * mrSges(6,3);
t505 = m(6) * t507 - t601 * mrSges(6,2) + t516 * mrSges(6,3) + t540 * t528 - t602 * t536;
t495 = t614 * t504 + t611 * t505;
t544 = -t562 * mrSges(5,1) + t563 * mrSges(5,2);
t551 = qJD(2) * mrSges(5,2) + t562 * mrSges(5,3);
t493 = m(5) * t510 - qJDD(2) * mrSges(5,1) - t548 * mrSges(5,3) - qJD(2) * t551 - t563 * t544 + t495;
t552 = -qJD(2) * mrSges(5,1) - t563 * mrSges(5,3);
t625 = -t611 * t504 + t614 * t505;
t494 = m(5) * t511 + qJDD(2) * mrSges(5,2) + t547 * mrSges(5,3) + qJD(2) * t552 + t562 * t544 + t625;
t626 = -t608 * t493 + t609 * t494;
t622 = m(4) * t533 + qJDD(2) * mrSges(4,3) + qJD(2) * t589 + t581 * t632 + t626;
t489 = m(3) * t550 - qJDD(2) * mrSges(3,2) - qJD(2) * t588 + t582 * t632 + t642 * t584 + t622;
t590 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t632;
t491 = t609 * t493 + t608 * t494;
t591 = mrSges(4,2) * t632 + qJD(2) * mrSges(4,3);
t621 = -m(4) * t534 + qJDD(2) * mrSges(4,1) + qJD(2) * t591 - t491;
t490 = m(3) * t549 + qJDD(2) * mrSges(3,1) + qJD(2) * t590 - t642 * t583 + (-t581 - t582) * t633 + t621;
t627 = t615 * t489 - t612 * t490;
t482 = m(2) * t593 - t618 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t627;
t592 = t613 * g(1) - t616 * g(2);
t570 = -qJDD(1) * pkin(1) - t618 * pkin(6) - t592;
t623 = -t584 * pkin(2) + t570 + (-t583 - t629) * qJ(3);
t531 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t633 + t623;
t519 = -pkin(2) * t630 + t584 * pkin(3) - qJ(4) * t638 + qJDD(4) - t623 + (t587 + t643) * t633;
t513 = -t547 * pkin(4) - t561 * pkin(7) + t563 * t553 + t519;
t624 = m(6) * t513 - t516 * mrSges(6,1) + t517 * mrSges(6,2) - t540 * t535 + t541 * t536;
t620 = -m(5) * t519 + t547 * mrSges(5,1) - t548 * mrSges(5,2) + t562 * t551 - t563 * t552 - t624;
t500 = m(4) * t531 - t584 * mrSges(4,1) - t583 * mrSges(4,3) - t589 * t633 - t591 * t632 + t620;
t619 = -m(3) * t570 + t584 * mrSges(3,1) - t583 * mrSges(3,2) - t588 * t633 + t590 * t632 - t500;
t499 = m(2) * t592 + qJDD(1) * mrSges(2,1) - t618 * mrSges(2,2) + t619;
t637 = t613 * t482 + t616 * t499;
t483 = t612 * t489 + t615 * t490;
t636 = t644 * qJD(2) + (t640 * t612 + t639 * t615) * qJD(1);
t635 = -t639 * qJD(2) + (-t641 * t612 - t645 * t615) * qJD(1);
t634 = t640 * qJD(2) + (t646 * t612 + t641 * t615) * qJD(1);
t628 = t616 * t482 - t613 * t499;
t539 = Ifges(5,1) * t563 + Ifges(5,4) * t562 - Ifges(5,5) * qJD(2);
t538 = Ifges(5,4) * t563 + Ifges(5,2) * t562 - Ifges(5,6) * qJD(2);
t537 = Ifges(5,5) * t563 + Ifges(5,6) * t562 - Ifges(5,3) * qJD(2);
t522 = Ifges(6,1) * t541 + Ifges(6,4) * t540 + Ifges(6,5) * t602;
t521 = Ifges(6,4) * t541 + Ifges(6,2) * t540 + Ifges(6,6) * t602;
t520 = Ifges(6,5) * t541 + Ifges(6,6) * t540 + Ifges(6,3) * t602;
t497 = mrSges(6,2) * t513 - mrSges(6,3) * t506 + Ifges(6,1) * t517 + Ifges(6,4) * t516 + Ifges(6,5) * t601 + t540 * t520 - t602 * t521;
t496 = -mrSges(6,1) * t513 + mrSges(6,3) * t507 + Ifges(6,4) * t517 + Ifges(6,2) * t516 + Ifges(6,6) * t601 - t541 * t520 + t602 * t522;
t485 = mrSges(5,2) * t519 - mrSges(5,3) * t510 + Ifges(5,1) * t548 + Ifges(5,4) * t547 - Ifges(5,5) * qJDD(2) - pkin(7) * t495 + qJD(2) * t538 - t611 * t496 + t614 * t497 + t562 * t537;
t484 = -mrSges(5,1) * t519 + mrSges(5,3) * t511 + Ifges(5,4) * t548 + Ifges(5,2) * t547 - Ifges(5,6) * qJDD(2) - pkin(4) * t624 + pkin(7) * t625 - qJD(2) * t539 + t614 * t496 + t611 * t497 - t563 * t537;
t479 = mrSges(3,2) * t570 + mrSges(4,2) * t534 - mrSges(3,3) * t549 - mrSges(4,3) * t531 - qJ(3) * t500 - qJ(4) * t491 + t635 * qJD(2) + t640 * qJDD(2) - t608 * t484 + t609 * t485 + t646 * t583 + t641 * t584 + t636 * t632;
t478 = -mrSges(3,1) * t570 - mrSges(4,1) * t531 + mrSges(4,2) * t533 + mrSges(3,3) * t550 - pkin(2) * t500 - pkin(3) * t620 - qJ(4) * t626 + t634 * qJD(2) + t639 * qJDD(2) - t609 * t484 - t608 * t485 + t641 * t583 + t645 * t584 - t636 * t633;
t477 = Ifges(6,6) * t516 + Ifges(6,5) * t517 + mrSges(5,1) * t510 - mrSges(5,2) * t511 + mrSges(6,1) * t506 - mrSges(6,2) * t507 + mrSges(2,3) * t593 + Ifges(6,3) * t601 + pkin(4) * t495 + (-qJ(3) * mrSges(4,2) - t639) * t584 + (pkin(2) * mrSges(4,2) - t640) * t583 + t618 * Ifges(2,5) + pkin(3) * t491 - t562 * t539 + t563 * t538 + Ifges(5,6) * t547 + Ifges(5,5) * t548 - mrSges(3,1) * t549 + mrSges(3,2) * t550 + mrSges(4,1) * t534 - t540 * t522 + t541 * t521 - mrSges(4,3) * t533 - pkin(2) * t621 - qJ(3) * t622 + (t634 * t615 + (pkin(2) * t581 + t635) * t612) * qJD(1) + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + (-Ifges(5,3) - t644) * qJDD(2) - pkin(1) * t483;
t476 = -mrSges(2,2) * g(3) - mrSges(2,3) * t592 + Ifges(2,5) * qJDD(1) - t618 * Ifges(2,6) - pkin(6) * t483 - t612 * t478 + t615 * t479;
t1 = [-m(1) * g(1) + t628; -m(1) * g(2) + t637; (-m(1) - m(2)) * g(3) + t483; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t637 + t616 * t476 - t613 * t477; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t628 + t613 * t476 + t616 * t477; -mrSges(1,1) * g(2) + mrSges(2,1) * t592 + mrSges(1,2) * g(1) - mrSges(2,2) * t593 + Ifges(2,3) * qJDD(1) + pkin(1) * t619 + pkin(6) * t627 + t615 * t478 + t612 * t479;];
tauB = t1;
