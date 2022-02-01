% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:21
% EndTime: 2022-01-23 08:59:26
% DurationCPUTime: 4.57s
% Computational Cost: add. (36210->282), mult. (101411->386), div. (0->0), fcn. (68373->10), ass. (0->128)
t563 = sin(pkin(7));
t566 = cos(pkin(7));
t568 = sin(qJ(1));
t570 = cos(qJ(1));
t548 = -t570 * g(1) - t568 * g(2);
t571 = qJD(1) ^ 2;
t616 = -t571 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t548;
t519 = -t563 * g(3) + t616 * t566;
t582 = -pkin(2) * t566 - qJ(3) * t563;
t542 = t582 * qJD(1);
t601 = t566 * qJD(1);
t510 = t542 * t601 + t519;
t562 = sin(pkin(8));
t565 = cos(pkin(8));
t547 = t568 * g(1) - t570 * g(2);
t577 = -t571 * qJ(2) + qJDD(2) - t547;
t602 = t563 * qJD(1);
t615 = (-pkin(1) + t582) * qJDD(1) + t577 - 0.2e1 * qJD(3) * t602;
t490 = -t562 * t510 + t615 * t565;
t518 = -t566 * g(3) - t616 * t563;
t561 = sin(pkin(9));
t564 = cos(pkin(9));
t606 = t563 * t565;
t578 = t561 * t606 + t564 * t566;
t531 = t578 * qJD(1);
t529 = t578 * qJDD(1);
t614 = 2 * qJD(4);
t613 = mrSges(3,2) * t563;
t612 = Ifges(4,4) * t565;
t611 = Ifges(3,6) * t566;
t610 = Ifges(4,6) * t566;
t609 = t563 ^ 2 * t571;
t608 = t566 ^ 2 * t571;
t607 = t562 * t563;
t605 = t566 * t571;
t543 = (-mrSges(3,1) * t566 + t613) * qJD(1);
t491 = t565 * t510 + t615 * t562;
t585 = mrSges(4,1) * t562 + mrSges(4,2) * t565;
t536 = t585 * t602;
t580 = -mrSges(4,1) * t566 - mrSges(4,3) * t606;
t540 = t580 * qJD(1);
t579 = mrSges(4,2) * t566 - mrSges(4,3) * t607;
t535 = (pkin(3) * t562 - qJ(4) * t565) * t602;
t596 = t562 * t602;
t598 = qJDD(1) * t566;
t488 = -pkin(3) * t608 - qJ(4) * t598 - t535 * t596 + t491;
t509 = t542 * t602 + qJDD(3) - t518;
t496 = ((-qJDD(1) * t565 - t562 * t605) * qJ(4) + (qJDD(1) * t562 - t565 * t605) * pkin(3)) * t563 + t509;
t484 = t564 * t488 + t561 * t496 - t531 * t614;
t595 = t565 * t602;
t532 = -t561 * t601 + t564 * t595;
t511 = t531 * mrSges(5,1) + t532 * mrSges(5,2);
t517 = mrSges(5,1) * t596 - t532 * mrSges(5,3);
t512 = t531 * pkin(4) - t532 * pkin(6);
t599 = qJDD(1) * t563;
t592 = t562 * t599;
t597 = t562 ^ 2 * t609;
t482 = -pkin(4) * t597 + pkin(6) * t592 - t531 * t512 + t484;
t487 = pkin(3) * t598 - qJ(4) * t608 + t535 * t595 + qJDD(4) - t490;
t530 = (-t561 * t566 + t564 * t606) * qJDD(1);
t485 = (t531 * t596 - t530) * pkin(6) + (t532 * t596 + t529) * pkin(4) + t487;
t567 = sin(qJ(5));
t569 = cos(qJ(5));
t479 = -t567 * t482 + t569 * t485;
t513 = -t567 * t532 + t569 * t596;
t514 = t569 * t532 + t567 * t596;
t498 = -t513 * mrSges(6,1) + t514 * mrSges(6,2);
t500 = t513 * qJD(5) + t569 * t530 + t567 * t592;
t528 = qJD(5) + t531;
t501 = -t528 * mrSges(6,2) + t513 * mrSges(6,3);
t527 = qJDD(5) + t529;
t477 = m(6) * t479 + t527 * mrSges(6,1) - t500 * mrSges(6,3) - t514 * t498 + t528 * t501;
t480 = t569 * t482 + t567 * t485;
t499 = -t514 * qJD(5) - t567 * t530 + t569 * t592;
t502 = t528 * mrSges(6,1) - t514 * mrSges(6,3);
t478 = m(6) * t480 - t527 * mrSges(6,2) + t499 * mrSges(6,3) + t513 * t498 - t528 * t502;
t587 = -t567 * t477 + t569 * t478;
t471 = m(5) * t484 - t529 * mrSges(5,3) - t531 * t511 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t517) * t607 + t587;
t581 = t561 * t488 - t564 * t496;
t483 = -0.2e1 * qJD(4) * t532 - t581;
t516 = -mrSges(5,2) * t596 - t531 * mrSges(5,3);
t481 = -pkin(4) * t592 - pkin(6) * t597 + (t614 + t512) * t532 + t581;
t575 = -m(6) * t481 + t499 * mrSges(6,1) - t500 * mrSges(6,2) + t513 * t501 - t514 * t502;
t475 = m(5) * t483 - t530 * mrSges(5,3) - t532 * t511 + (mrSges(5,1) * qJDD(1) + qJD(1) * t516) * t607 + t575;
t588 = t564 * t471 - t561 * t475;
t467 = m(4) * t491 + t579 * qJDD(1) + (-t536 * t607 + t540 * t566) * qJD(1) + t588;
t539 = t579 * qJD(1);
t472 = t569 * t477 + t567 * t478;
t572 = -m(5) * t487 - t529 * mrSges(5,1) - t530 * mrSges(5,2) - t531 * t516 - t532 * t517 - t472;
t469 = m(4) * t490 + t580 * qJDD(1) + (-t536 * t606 - t539 * t566) * qJD(1) + t572;
t589 = t565 * t467 - t562 * t469;
t460 = m(3) * t519 + (qJDD(1) * mrSges(3,3) + qJD(1) * t543) * t566 + t589;
t468 = t561 * t471 + t564 * t475;
t576 = -m(4) * t509 - t468;
t465 = m(3) * t518 + ((-mrSges(3,3) - t585) * qJDD(1) + (-t539 * t562 - t540 * t565 - t543) * qJD(1)) * t563 + t576;
t590 = t566 * t460 - t563 * t465;
t454 = m(2) * t548 - t571 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t590;
t461 = t562 * t467 + t565 * t469;
t538 = -qJDD(1) * pkin(1) + t577;
t573 = -m(3) * t538 + mrSges(3,1) * t598 - t461 + (t608 + t609) * mrSges(3,3);
t457 = m(2) * t547 - t571 * mrSges(2,2) + (mrSges(2,1) - t613) * qJDD(1) + t573;
t604 = t568 * t454 + t570 * t457;
t455 = t563 * t460 + t566 * t465;
t591 = t570 * t454 - t568 * t457;
t584 = Ifges(3,1) * t563 + Ifges(3,4) * t566;
t583 = Ifges(4,5) * t565 - Ifges(4,6) * t562;
t574 = -Ifges(4,5) * t566 + (Ifges(4,1) * t565 - Ifges(4,4) * t562) * t563;
t544 = (Ifges(3,5) * t563 + t611) * qJD(1);
t524 = t574 * qJD(1);
t523 = (-t610 + (-Ifges(4,2) * t562 + t612) * t563) * qJD(1);
t522 = (-Ifges(4,3) * t566 + t583 * t563) * qJD(1);
t505 = Ifges(5,1) * t532 - Ifges(5,4) * t531 + Ifges(5,5) * t596;
t504 = Ifges(5,4) * t532 - Ifges(5,2) * t531 + Ifges(5,6) * t596;
t503 = Ifges(5,5) * t532 - Ifges(5,6) * t531 + Ifges(5,3) * t596;
t494 = Ifges(6,1) * t514 + Ifges(6,4) * t513 + Ifges(6,5) * t528;
t493 = Ifges(6,4) * t514 + Ifges(6,2) * t513 + Ifges(6,6) * t528;
t492 = Ifges(6,5) * t514 + Ifges(6,6) * t513 + Ifges(6,3) * t528;
t474 = mrSges(6,2) * t481 - mrSges(6,3) * t479 + Ifges(6,1) * t500 + Ifges(6,4) * t499 + Ifges(6,5) * t527 + t513 * t492 - t528 * t493;
t473 = -mrSges(6,1) * t481 + mrSges(6,3) * t480 + Ifges(6,4) * t500 + Ifges(6,2) * t499 + Ifges(6,6) * t527 - t514 * t492 + t528 * t494;
t463 = -mrSges(5,1) * t487 - mrSges(6,1) * t479 + mrSges(6,2) * t480 + mrSges(5,3) * t484 + Ifges(5,4) * t530 - Ifges(6,5) * t500 - Ifges(5,2) * t529 - Ifges(6,6) * t499 - Ifges(6,3) * t527 - pkin(4) * t472 - t514 * t493 + t513 * t494 - t532 * t503 + (Ifges(5,6) * qJDD(1) + qJD(1) * t505) * t607;
t462 = mrSges(5,2) * t487 - mrSges(5,3) * t483 + Ifges(5,1) * t530 - Ifges(5,4) * t529 - pkin(6) * t472 - t567 * t473 + t569 * t474 - t531 * t503 + (Ifges(5,5) * qJDD(1) - qJD(1) * t504) * t607;
t451 = -mrSges(4,1) * t509 + mrSges(4,3) * t491 - Ifges(5,5) * t530 + Ifges(5,6) * t529 - t532 * t504 - t531 * t505 - mrSges(5,1) * t483 + mrSges(5,2) * t484 - t567 * t474 - t569 * t473 - pkin(4) * t575 - pkin(6) * t587 - pkin(3) * t468 + (-t522 * t606 - t566 * t524) * qJD(1) + (-t610 + (t612 + (-Ifges(4,2) - Ifges(5,3)) * t562) * t563) * qJDD(1);
t450 = mrSges(4,2) * t509 - mrSges(4,3) * t490 - qJ(4) * t468 + t564 * t462 - t561 * t463 + (-t522 * t607 + t523 * t566) * qJD(1) + t574 * qJDD(1);
t449 = -mrSges(3,1) * t538 + mrSges(3,3) * t519 - mrSges(4,1) * t490 + mrSges(4,2) * t491 - t561 * t462 - t564 * t463 - pkin(3) * t572 - qJ(4) * t588 - pkin(2) * t461 + (Ifges(3,2) + Ifges(4,3)) * t598 + ((Ifges(3,4) - t583) * qJDD(1) + (-t523 * t565 - t524 * t562 - t544) * qJD(1)) * t563;
t448 = mrSges(3,2) * t538 - mrSges(3,3) * t518 - qJ(3) * t461 + t584 * qJDD(1) + t565 * t450 - t562 * t451 + t544 * t601;
t447 = t571 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t548 - mrSges(3,1) * t518 + mrSges(3,2) * t519 - t562 * t450 - t565 * t451 - pkin(2) * t576 - qJ(3) * t589 - pkin(1) * t455 + (-t611 + Ifges(2,6) + (pkin(2) * t585 - Ifges(3,5)) * t563) * qJDD(1) + (-pkin(2) * (-t539 * t607 - t540 * t606) + (-t563 * (Ifges(3,4) * t563 + Ifges(3,2) * t566) + t566 * t584) * qJD(1)) * qJD(1);
t446 = -mrSges(2,2) * g(3) - mrSges(2,3) * t547 + Ifges(2,5) * qJDD(1) - t571 * Ifges(2,6) - qJ(2) * t455 + t566 * t448 - t563 * t449;
t1 = [-m(1) * g(1) + t591; -m(1) * g(2) + t604; (-m(1) - m(2)) * g(3) + t455; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t604 + t570 * t446 - t568 * t447; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t591 + t568 * t446 + t570 * t447; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t547 - mrSges(2,2) * t548 + t563 * t448 + t566 * t449 + pkin(1) * (-mrSges(3,2) * t599 + t573) + qJ(2) * t590;];
tauB = t1;
