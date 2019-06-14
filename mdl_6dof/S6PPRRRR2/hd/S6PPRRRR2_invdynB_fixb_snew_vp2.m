% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PPRRRR2
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PPRRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:55:23
% EndTime: 2019-05-04 20:55:46
% DurationCPUTime: 21.96s
% Computational Cost: add. (398700->294), mult. (726743->386), div. (0->0), fcn. (569540->16), ass. (0->134)
t579 = sin(pkin(12));
t583 = cos(pkin(12));
t571 = -t583 * g(1) - t579 * g(2);
t578 = sin(pkin(13));
t582 = cos(pkin(13));
t570 = t579 * g(1) - t583 * g(2);
t577 = -g(3) + qJDD(1);
t581 = sin(pkin(6));
t585 = cos(pkin(6));
t601 = t570 * t585 + t577 * t581;
t538 = -t578 * t571 + t582 * t601;
t539 = t582 * t571 + t578 * t601;
t554 = -t581 * t570 + t585 * t577 + qJDD(2);
t593 = cos(qJ(3));
t584 = cos(pkin(7));
t589 = sin(qJ(3));
t612 = t584 * t589;
t580 = sin(pkin(7));
t613 = t580 * t589;
t522 = t538 * t612 + t593 * t539 + t554 * t613;
t595 = qJD(3) ^ 2;
t520 = -t595 * pkin(3) + qJDD(3) * pkin(9) + t522;
t531 = -t580 * t538 + t584 * t554;
t588 = sin(qJ(4));
t592 = cos(qJ(4));
t513 = t592 * t520 + t588 * t531;
t566 = (-mrSges(5,1) * t592 + mrSges(5,2) * t588) * qJD(3);
t608 = qJD(3) * qJD(4);
t576 = t588 * t608;
t569 = t592 * qJDD(3) - t576;
t610 = qJD(3) * t588;
t572 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t610;
t567 = (-pkin(4) * t592 - pkin(10) * t588) * qJD(3);
t594 = qJD(4) ^ 2;
t609 = t592 * qJD(3);
t511 = -t594 * pkin(4) + qJDD(4) * pkin(10) + t567 * t609 + t513;
t521 = -t589 * t539 + (t538 * t584 + t554 * t580) * t593;
t519 = -qJDD(3) * pkin(3) - t595 * pkin(9) - t521;
t607 = t592 * t608;
t568 = t588 * qJDD(3) + t607;
t516 = (-t568 - t607) * pkin(10) + (-t569 + t576) * pkin(4) + t519;
t587 = sin(qJ(5));
t591 = cos(qJ(5));
t506 = -t587 * t511 + t591 * t516;
t564 = t591 * qJD(4) - t587 * t610;
t546 = t564 * qJD(5) + t587 * qJDD(4) + t591 * t568;
t563 = qJDD(5) - t569;
t565 = t587 * qJD(4) + t591 * t610;
t575 = qJD(5) - t609;
t504 = (t564 * t575 - t546) * pkin(11) + (t564 * t565 + t563) * pkin(5) + t506;
t507 = t591 * t511 + t587 * t516;
t545 = -t565 * qJD(5) + t591 * qJDD(4) - t587 * t568;
t553 = t575 * pkin(5) - t565 * pkin(11);
t562 = t564 ^ 2;
t505 = -t562 * pkin(5) + t545 * pkin(11) - t575 * t553 + t507;
t586 = sin(qJ(6));
t590 = cos(qJ(6));
t502 = t590 * t504 - t586 * t505;
t547 = t590 * t564 - t586 * t565;
t525 = t547 * qJD(6) + t586 * t545 + t590 * t546;
t548 = t586 * t564 + t590 * t565;
t532 = -mrSges(7,1) * t547 + mrSges(7,2) * t548;
t574 = qJD(6) + t575;
t534 = -t574 * mrSges(7,2) + t547 * mrSges(7,3);
t559 = qJDD(6) + t563;
t500 = m(7) * t502 + t559 * mrSges(7,1) - t525 * mrSges(7,3) - t548 * t532 + t574 * t534;
t503 = t586 * t504 + t590 * t505;
t524 = -t548 * qJD(6) + t590 * t545 - t586 * t546;
t535 = t574 * mrSges(7,1) - t548 * mrSges(7,3);
t501 = m(7) * t503 - t559 * mrSges(7,2) + t524 * mrSges(7,3) + t547 * t532 - t574 * t535;
t492 = t590 * t500 + t586 * t501;
t549 = -t564 * mrSges(6,1) + t565 * mrSges(6,2);
t551 = -t575 * mrSges(6,2) + t564 * mrSges(6,3);
t490 = m(6) * t506 + t563 * mrSges(6,1) - t546 * mrSges(6,3) - t565 * t549 + t575 * t551 + t492;
t552 = t575 * mrSges(6,1) - t565 * mrSges(6,3);
t603 = -t586 * t500 + t590 * t501;
t491 = m(6) * t507 - t563 * mrSges(6,2) + t545 * mrSges(6,3) + t564 * t549 - t575 * t552 + t603;
t604 = -t587 * t490 + t591 * t491;
t487 = m(5) * t513 - qJDD(4) * mrSges(5,2) + t569 * mrSges(5,3) - qJD(4) * t572 + t566 * t609 + t604;
t512 = -t588 * t520 + t592 * t531;
t573 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t609;
t510 = -qJDD(4) * pkin(4) - t594 * pkin(10) + t567 * t610 - t512;
t508 = -t545 * pkin(5) - t562 * pkin(11) + t565 * t553 + t510;
t598 = m(7) * t508 - t524 * mrSges(7,1) + t525 * mrSges(7,2) - t547 * t534 + t548 * t535;
t596 = -m(6) * t510 + t545 * mrSges(6,1) - t546 * mrSges(6,2) + t564 * t551 - t565 * t552 - t598;
t496 = m(5) * t512 + qJDD(4) * mrSges(5,1) - t568 * mrSges(5,3) + qJD(4) * t573 - t566 * t610 + t596;
t605 = t592 * t487 - t588 * t496;
t476 = m(4) * t522 - t595 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t605;
t479 = t588 * t487 + t592 * t496;
t478 = m(4) * t531 + t479;
t488 = t591 * t490 + t587 * t491;
t597 = -m(5) * t519 + t569 * mrSges(5,1) - t568 * mrSges(5,2) - t572 * t610 + t573 * t609 - t488;
t484 = m(4) * t521 + qJDD(3) * mrSges(4,1) - t595 * mrSges(4,2) + t597;
t614 = t484 * t593;
t465 = t476 * t612 - t580 * t478 + t584 * t614;
t461 = m(3) * t538 + t465;
t471 = t593 * t476 - t589 * t484;
t470 = m(3) * t539 + t471;
t617 = t461 * t582 + t470 * t578;
t464 = t476 * t613 + t584 * t478 + t580 * t614;
t463 = m(3) * t554 + t464;
t451 = -t581 * t463 + t585 * t617;
t449 = m(2) * t570 + t451;
t457 = -t578 * t461 + t582 * t470;
t456 = m(2) * t571 + t457;
t611 = t583 * t449 + t579 * t456;
t450 = t585 * t463 + t581 * t617;
t606 = -t579 * t449 + t583 * t456;
t526 = Ifges(7,5) * t548 + Ifges(7,6) * t547 + Ifges(7,3) * t574;
t528 = Ifges(7,1) * t548 + Ifges(7,4) * t547 + Ifges(7,5) * t574;
t493 = -mrSges(7,1) * t508 + mrSges(7,3) * t503 + Ifges(7,4) * t525 + Ifges(7,2) * t524 + Ifges(7,6) * t559 - t548 * t526 + t574 * t528;
t527 = Ifges(7,4) * t548 + Ifges(7,2) * t547 + Ifges(7,6) * t574;
t494 = mrSges(7,2) * t508 - mrSges(7,3) * t502 + Ifges(7,1) * t525 + Ifges(7,4) * t524 + Ifges(7,5) * t559 + t547 * t526 - t574 * t527;
t540 = Ifges(6,5) * t565 + Ifges(6,6) * t564 + Ifges(6,3) * t575;
t542 = Ifges(6,1) * t565 + Ifges(6,4) * t564 + Ifges(6,5) * t575;
t480 = -mrSges(6,1) * t510 + mrSges(6,3) * t507 + Ifges(6,4) * t546 + Ifges(6,2) * t545 + Ifges(6,6) * t563 - pkin(5) * t598 + pkin(11) * t603 + t590 * t493 + t586 * t494 - t565 * t540 + t575 * t542;
t541 = Ifges(6,4) * t565 + Ifges(6,2) * t564 + Ifges(6,6) * t575;
t481 = mrSges(6,2) * t510 - mrSges(6,3) * t506 + Ifges(6,1) * t546 + Ifges(6,4) * t545 + Ifges(6,5) * t563 - pkin(11) * t492 - t586 * t493 + t590 * t494 + t564 * t540 - t575 * t541;
t556 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t588 + Ifges(5,6) * t592) * qJD(3);
t557 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t588 + Ifges(5,2) * t592) * qJD(3);
t466 = mrSges(5,2) * t519 - mrSges(5,3) * t512 + Ifges(5,1) * t568 + Ifges(5,4) * t569 + Ifges(5,5) * qJDD(4) - pkin(10) * t488 - qJD(4) * t557 - t587 * t480 + t591 * t481 + t556 * t609;
t558 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t588 + Ifges(5,4) * t592) * qJD(3);
t472 = Ifges(5,4) * t568 + Ifges(5,2) * t569 + Ifges(5,6) * qJDD(4) - t556 * t610 + qJD(4) * t558 - mrSges(5,1) * t519 + mrSges(5,3) * t513 - Ifges(6,5) * t546 - Ifges(6,6) * t545 - Ifges(6,3) * t563 - t565 * t541 + t564 * t542 - mrSges(6,1) * t506 + mrSges(6,2) * t507 - Ifges(7,5) * t525 - Ifges(7,6) * t524 - Ifges(7,3) * t559 - t548 * t527 + t547 * t528 - mrSges(7,1) * t502 + mrSges(7,2) * t503 - pkin(5) * t492 - pkin(4) * t488;
t453 = mrSges(4,2) * t531 - mrSges(4,3) * t521 + Ifges(4,5) * qJDD(3) - t595 * Ifges(4,6) - pkin(9) * t479 + t592 * t466 - t588 * t472;
t458 = Ifges(4,6) * qJDD(3) + t595 * Ifges(4,5) - mrSges(4,1) * t531 + mrSges(4,3) * t522 - Ifges(5,5) * t568 - Ifges(5,6) * t569 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t512 + mrSges(5,2) * t513 - t587 * t481 - t591 * t480 - pkin(4) * t596 - pkin(10) * t604 - pkin(3) * t479 + (-t588 * t557 + t592 * t558) * qJD(3);
t600 = pkin(8) * t471 + t453 * t589 + t458 * t593;
t452 = mrSges(4,1) * t521 - mrSges(4,2) * t522 + Ifges(4,3) * qJDD(3) + pkin(3) * t597 + pkin(9) * t605 + t588 * t466 + t592 * t472;
t446 = -mrSges(3,1) * t554 + mrSges(3,3) * t539 - pkin(2) * t464 - t580 * t452 + t584 * t600;
t447 = mrSges(3,2) * t554 - mrSges(3,3) * t538 + t593 * t453 - t589 * t458 + (-t464 * t580 - t465 * t584) * pkin(8);
t599 = qJ(2) * t457 + t446 * t582 + t447 * t578;
t445 = mrSges(3,1) * t538 - mrSges(3,2) * t539 + pkin(2) * t465 + t584 * t452 + t580 * t600;
t444 = mrSges(2,2) * t577 - mrSges(2,3) * t570 - t578 * t446 + t582 * t447 + (-t450 * t581 - t451 * t585) * qJ(2);
t443 = -mrSges(2,1) * t577 + mrSges(2,3) * t571 - pkin(1) * t450 - t581 * t445 + t585 * t599;
t1 = [-m(1) * g(1) + t606; -m(1) * g(2) + t611; -m(1) * g(3) + m(2) * t577 + t450; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t611 - t579 * t443 + t583 * t444; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t606 + t583 * t443 + t579 * t444; -mrSges(1,1) * g(2) + mrSges(2,1) * t570 + mrSges(1,2) * g(1) - mrSges(2,2) * t571 + pkin(1) * t451 + t585 * t445 + t581 * t599;];
tauB  = t1;
