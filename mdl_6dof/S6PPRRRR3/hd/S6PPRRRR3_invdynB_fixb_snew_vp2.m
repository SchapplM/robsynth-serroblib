% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 21:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PPRRRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_invdynB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:13:47
% EndTime: 2019-05-04 21:14:32
% DurationCPUTime: 44.56s
% Computational Cost: add. (855540->306), mult. (1627082->416), div. (0->0), fcn. (1378843->18), ass. (0->147)
t611 = sin(pkin(13));
t616 = cos(pkin(13));
t605 = -g(1) * t616 - g(2) * t611;
t610 = sin(pkin(14));
t615 = cos(pkin(14));
t604 = g(1) * t611 - g(2) * t616;
t609 = -g(3) + qJDD(1);
t614 = sin(pkin(6));
t619 = cos(pkin(6));
t635 = t604 * t619 + t609 * t614;
t579 = -t610 * t605 + t615 * t635;
t580 = t615 * t605 + t610 * t635;
t592 = -t604 * t614 + t609 * t619 + qJDD(2);
t623 = sin(qJ(3));
t618 = cos(pkin(7));
t627 = cos(qJ(3));
t645 = t618 * t627;
t613 = sin(pkin(7));
t648 = t613 * t627;
t554 = t579 * t645 - t580 * t623 + t592 * t648;
t628 = qJD(3) ^ 2;
t612 = sin(pkin(8));
t654 = pkin(10) * t612;
t549 = qJDD(3) * pkin(3) + t628 * t654 + t554;
t646 = t618 * t623;
t649 = t613 * t623;
t555 = t579 * t646 + t627 * t580 + t592 * t649;
t550 = -pkin(3) * t628 + qJDD(3) * t654 + t555;
t563 = -t579 * t613 + t592 * t618;
t626 = cos(qJ(4));
t617 = cos(pkin(8));
t622 = sin(qJ(4));
t647 = t617 * t622;
t650 = t612 * t622;
t542 = t549 * t647 + t626 * t550 + t563 * t650;
t608 = qJD(3) * t617 + qJD(4);
t643 = qJD(3) * t612;
t641 = t622 * t643;
t596 = mrSges(5,1) * t608 - mrSges(5,3) * t641;
t598 = (-mrSges(5,1) * t626 + mrSges(5,2) * t622) * t643;
t642 = qJD(3) * qJD(4);
t601 = (-qJDD(3) * t626 + t622 * t642) * t612;
t607 = qJDD(3) * t617 + qJDD(4);
t599 = (-pkin(4) * t626 - pkin(11) * t622) * t643;
t606 = t608 ^ 2;
t640 = t626 * t643;
t540 = -pkin(4) * t606 + pkin(11) * t607 + t599 * t640 + t542;
t562 = t617 * t563;
t600 = (qJDD(3) * t622 + t626 * t642) * t612;
t544 = t601 * pkin(4) - t600 * pkin(11) + t562 + (-t549 + (pkin(4) * t622 - pkin(11) * t626) * t608 * qJD(3)) * t612;
t621 = sin(qJ(5));
t625 = cos(qJ(5));
t536 = t625 * t540 + t621 * t544;
t594 = t608 * t621 + t625 * t641;
t572 = -qJD(5) * t594 - t600 * t621 + t607 * t625;
t593 = t608 * t625 - t621 * t641;
t574 = -mrSges(6,1) * t593 + mrSges(6,2) * t594;
t603 = qJD(5) - t640;
t584 = mrSges(6,1) * t603 - mrSges(6,3) * t594;
t595 = qJDD(5) + t601;
t575 = -pkin(5) * t593 - pkin(12) * t594;
t602 = t603 ^ 2;
t534 = -pkin(5) * t602 + pkin(12) * t595 + t575 * t593 + t536;
t541 = -t622 * t550 + (t549 * t617 + t563 * t612) * t626;
t539 = -t607 * pkin(4) - t606 * pkin(11) + t599 * t641 - t541;
t573 = qJD(5) * t593 + t600 * t625 + t607 * t621;
t537 = (-t593 * t603 - t573) * pkin(12) + (t594 * t603 - t572) * pkin(5) + t539;
t620 = sin(qJ(6));
t624 = cos(qJ(6));
t531 = -t534 * t620 + t537 * t624;
t581 = -t594 * t620 + t603 * t624;
t553 = qJD(6) * t581 + t573 * t624 + t595 * t620;
t582 = t594 * t624 + t603 * t620;
t560 = -mrSges(7,1) * t581 + mrSges(7,2) * t582;
t591 = qJD(6) - t593;
t564 = -mrSges(7,2) * t591 + mrSges(7,3) * t581;
t570 = qJDD(6) - t572;
t529 = m(7) * t531 + mrSges(7,1) * t570 - mrSges(7,3) * t553 - t560 * t582 + t564 * t591;
t532 = t534 * t624 + t537 * t620;
t552 = -qJD(6) * t582 - t573 * t620 + t595 * t624;
t565 = mrSges(7,1) * t591 - mrSges(7,3) * t582;
t530 = m(7) * t532 - mrSges(7,2) * t570 + mrSges(7,3) * t552 + t560 * t581 - t565 * t591;
t637 = -t529 * t620 + t624 * t530;
t522 = m(6) * t536 - mrSges(6,2) * t595 + mrSges(6,3) * t572 + t574 * t593 - t584 * t603 + t637;
t535 = -t540 * t621 + t544 * t625;
t583 = -mrSges(6,2) * t603 + mrSges(6,3) * t593;
t533 = -pkin(5) * t595 - pkin(12) * t602 + t575 * t594 - t535;
t630 = -m(7) * t533 + t552 * mrSges(7,1) - mrSges(7,2) * t553 + t581 * t564 - t565 * t582;
t527 = m(6) * t535 + mrSges(6,1) * t595 - mrSges(6,3) * t573 - t574 * t594 + t583 * t603 + t630;
t638 = t625 * t522 - t527 * t621;
t513 = m(5) * t542 - mrSges(5,2) * t607 - mrSges(5,3) * t601 - t596 * t608 + t598 * t640 + t638;
t516 = t621 * t522 + t625 * t527;
t545 = -t612 * t549 + t562;
t597 = -mrSges(5,2) * t608 + mrSges(5,3) * t640;
t515 = m(5) * t545 + t601 * mrSges(5,1) + t600 * mrSges(5,2) + (t596 * t622 - t597 * t626) * t643 + t516;
t523 = t529 * t624 + t530 * t620;
t629 = -m(6) * t539 + t572 * mrSges(6,1) - mrSges(6,2) * t573 + t593 * t583 - t584 * t594 - t523;
t519 = m(5) * t541 + mrSges(5,1) * t607 - mrSges(5,3) * t600 + t597 * t608 - t598 * t641 + t629;
t651 = t519 * t626;
t502 = t513 * t647 - t515 * t612 + t617 * t651;
t498 = m(4) * t554 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t628 + t502;
t501 = t513 * t650 + t617 * t515 + t612 * t651;
t500 = m(4) * t563 + t501;
t507 = t626 * t513 - t519 * t622;
t506 = m(4) * t555 - mrSges(4,1) * t628 - qJDD(3) * mrSges(4,2) + t507;
t487 = t498 * t645 - t500 * t613 + t506 * t646;
t483 = m(3) * t579 + t487;
t492 = -t498 * t623 + t627 * t506;
t491 = m(3) * t580 + t492;
t655 = t483 * t615 + t491 * t610;
t486 = t498 * t648 + t618 * t500 + t506 * t649;
t485 = m(3) * t592 + t486;
t473 = -t485 * t614 + t619 * t655;
t471 = m(2) * t604 + t473;
t480 = -t483 * t610 + t615 * t491;
t479 = m(2) * t605 + t480;
t644 = t616 * t471 + t611 * t479;
t472 = t619 * t485 + t614 * t655;
t639 = -t471 * t611 + t616 * t479;
t556 = Ifges(7,5) * t582 + Ifges(7,6) * t581 + Ifges(7,3) * t591;
t558 = Ifges(7,1) * t582 + Ifges(7,4) * t581 + Ifges(7,5) * t591;
t524 = -mrSges(7,1) * t533 + mrSges(7,3) * t532 + Ifges(7,4) * t553 + Ifges(7,2) * t552 + Ifges(7,6) * t570 - t556 * t582 + t558 * t591;
t557 = Ifges(7,4) * t582 + Ifges(7,2) * t581 + Ifges(7,6) * t591;
t525 = mrSges(7,2) * t533 - mrSges(7,3) * t531 + Ifges(7,1) * t553 + Ifges(7,4) * t552 + Ifges(7,5) * t570 + t556 * t581 - t557 * t591;
t566 = Ifges(6,5) * t594 + Ifges(6,6) * t593 + Ifges(6,3) * t603;
t567 = Ifges(6,4) * t594 + Ifges(6,2) * t593 + Ifges(6,6) * t603;
t508 = mrSges(6,2) * t539 - mrSges(6,3) * t535 + Ifges(6,1) * t573 + Ifges(6,4) * t572 + Ifges(6,5) * t595 - pkin(12) * t523 - t524 * t620 + t525 * t624 + t566 * t593 - t567 * t603;
t568 = Ifges(6,1) * t594 + Ifges(6,4) * t593 + Ifges(6,5) * t603;
t509 = -mrSges(6,1) * t539 - mrSges(7,1) * t531 + mrSges(7,2) * t532 + mrSges(6,3) * t536 + Ifges(6,4) * t573 - Ifges(7,5) * t553 + Ifges(6,2) * t572 + Ifges(6,6) * t595 - Ifges(7,6) * t552 - Ifges(7,3) * t570 - pkin(5) * t523 - t557 * t582 + t558 * t581 - t566 * t594 + t568 * t603;
t588 = Ifges(5,6) * t608 + (Ifges(5,4) * t622 + Ifges(5,2) * t626) * t643;
t589 = Ifges(5,5) * t608 + (Ifges(5,1) * t622 + Ifges(5,4) * t626) * t643;
t493 = Ifges(5,5) * t600 - Ifges(5,6) * t601 + Ifges(5,3) * t607 + mrSges(5,1) * t541 - mrSges(5,2) * t542 + t621 * t508 + t625 * t509 + pkin(4) * t629 + pkin(11) * t638 + (t588 * t622 - t589 * t626) * t643;
t587 = Ifges(5,3) * t608 + (Ifges(5,5) * t622 + Ifges(5,6) * t626) * t643;
t494 = mrSges(5,2) * t545 - mrSges(5,3) * t541 + Ifges(5,1) * t600 - Ifges(5,4) * t601 + Ifges(5,5) * t607 - pkin(11) * t516 + t508 * t625 - t509 * t621 + t587 * t640 - t588 * t608;
t495 = Ifges(5,4) * t600 - Ifges(5,2) * t601 + Ifges(5,6) * t607 - t587 * t641 + t608 * t589 - mrSges(5,1) * t545 + mrSges(5,3) * t542 - Ifges(6,5) * t573 - Ifges(6,6) * t572 - Ifges(6,3) * t595 - t594 * t567 + t593 * t568 - mrSges(6,1) * t535 + mrSges(6,2) * t536 - t620 * t525 - t624 * t524 - pkin(5) * t630 - pkin(12) * t637 - pkin(4) * t516;
t632 = pkin(10) * t507 + t494 * t622 + t495 * t626;
t475 = -mrSges(4,1) * t563 + mrSges(4,3) * t555 + t628 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t501 - t612 * t493 + t617 * t632;
t476 = mrSges(4,2) * t563 - mrSges(4,3) * t554 + Ifges(4,5) * qJDD(3) - t628 * Ifges(4,6) + t626 * t494 - t622 * t495 + (-t501 * t612 - t502 * t617) * pkin(10);
t633 = pkin(9) * t492 + t475 * t627 + t476 * t623;
t474 = mrSges(4,1) * t554 - mrSges(4,2) * t555 + Ifges(4,3) * qJDD(3) + pkin(3) * t502 + t617 * t493 + t612 * t632;
t468 = -mrSges(3,1) * t592 + mrSges(3,3) * t580 - pkin(2) * t486 - t613 * t474 + t618 * t633;
t469 = mrSges(3,2) * t592 - mrSges(3,3) * t579 - t623 * t475 + t627 * t476 + (-t486 * t613 - t487 * t618) * pkin(9);
t631 = qJ(2) * t480 + t468 * t615 + t469 * t610;
t467 = mrSges(3,1) * t579 - mrSges(3,2) * t580 + pkin(2) * t487 + t618 * t474 + t613 * t633;
t466 = mrSges(2,2) * t609 - mrSges(2,3) * t604 - t610 * t468 + t615 * t469 + (-t472 * t614 - t473 * t619) * qJ(2);
t465 = -mrSges(2,1) * t609 + mrSges(2,3) * t605 - pkin(1) * t472 - t614 * t467 + t619 * t631;
t1 = [-m(1) * g(1) + t639; -m(1) * g(2) + t644; -m(1) * g(3) + m(2) * t609 + t472; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t644 - t611 * t465 + t616 * t466; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t639 + t616 * t465 + t611 * t466; -mrSges(1,1) * g(2) + mrSges(2,1) * t604 + mrSges(1,2) * g(1) - mrSges(2,2) * t605 + pkin(1) * t473 + t619 * t467 + t614 * t631;];
tauB  = t1;
