% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR13
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR13_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:43
% EndTime: 2019-12-31 20:32:51
% DurationCPUTime: 6.92s
% Computational Cost: add. (101333->311), mult. (220446->394), div. (0->0), fcn. (152839->10), ass. (0->120)
t565 = sin(qJ(1));
t569 = cos(qJ(1));
t554 = -t569 * g(1) - t565 * g(2);
t571 = qJD(1) ^ 2;
t539 = -t571 * pkin(1) + qJDD(1) * pkin(6) + t554;
t564 = sin(qJ(2));
t568 = cos(qJ(2));
t523 = -t564 * g(3) + t568 * t539;
t548 = (-mrSges(3,1) * t568 + mrSges(3,2) * t564) * qJD(1);
t582 = qJD(1) * qJD(2);
t557 = t564 * t582;
t550 = t568 * qJDD(1) - t557;
t584 = qJD(1) * t564;
t551 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t584;
t553 = t565 * g(1) - t569 * g(2);
t538 = -qJDD(1) * pkin(1) - t571 * pkin(6) - t553;
t581 = t568 * t582;
t549 = t564 * qJDD(1) + t581;
t506 = (-t549 - t581) * qJ(3) + (-t550 + t557) * pkin(2) + t538;
t547 = (-pkin(2) * t568 - qJ(3) * t564) * qJD(1);
t570 = qJD(2) ^ 2;
t583 = t568 * qJD(1);
t509 = -t570 * pkin(2) + qJDD(2) * qJ(3) + t547 * t583 + t523;
t560 = sin(pkin(9));
t561 = cos(pkin(9));
t544 = t560 * qJD(2) + t561 * t584;
t489 = -0.2e1 * qJD(3) * t544 + t561 * t506 - t560 * t509;
t528 = t560 * qJDD(2) + t561 * t549;
t543 = t561 * qJD(2) - t560 * t584;
t480 = (-t543 * t583 - t528) * pkin(7) + (t543 * t544 - t550) * pkin(3) + t489;
t490 = 0.2e1 * qJD(3) * t543 + t560 * t506 + t561 * t509;
t527 = t561 * qJDD(2) - t560 * t549;
t529 = -pkin(3) * t583 - t544 * pkin(7);
t542 = t543 ^ 2;
t482 = -t542 * pkin(3) + t527 * pkin(7) + t529 * t583 + t490;
t563 = sin(qJ(4));
t567 = cos(qJ(4));
t470 = t567 * t480 - t563 * t482;
t519 = t567 * t543 - t563 * t544;
t496 = t519 * qJD(4) + t563 * t527 + t567 * t528;
t520 = t563 * t543 + t567 * t544;
t546 = qJDD(4) - t550;
t556 = qJD(4) - t583;
t468 = (t519 * t556 - t496) * pkin(8) + (t519 * t520 + t546) * pkin(4) + t470;
t471 = t563 * t480 + t567 * t482;
t495 = -t520 * qJD(4) + t567 * t527 - t563 * t528;
t512 = t556 * pkin(4) - t520 * pkin(8);
t518 = t519 ^ 2;
t469 = -t518 * pkin(4) + t495 * pkin(8) - t556 * t512 + t471;
t562 = sin(qJ(5));
t566 = cos(qJ(5));
t466 = t566 * t468 - t562 * t469;
t501 = t566 * t519 - t562 * t520;
t477 = t501 * qJD(5) + t562 * t495 + t566 * t496;
t502 = t562 * t519 + t566 * t520;
t488 = -t501 * mrSges(6,1) + t502 * mrSges(6,2);
t555 = qJD(5) + t556;
t492 = -t555 * mrSges(6,2) + t501 * mrSges(6,3);
t540 = qJDD(5) + t546;
t464 = m(6) * t466 + t540 * mrSges(6,1) - t477 * mrSges(6,3) - t502 * t488 + t555 * t492;
t467 = t562 * t468 + t566 * t469;
t476 = -t502 * qJD(5) + t566 * t495 - t562 * t496;
t493 = t555 * mrSges(6,1) - t502 * mrSges(6,3);
t465 = m(6) * t467 - t540 * mrSges(6,2) + t476 * mrSges(6,3) + t501 * t488 - t555 * t493;
t456 = t566 * t464 + t562 * t465;
t503 = -t519 * mrSges(5,1) + t520 * mrSges(5,2);
t510 = -t556 * mrSges(5,2) + t519 * mrSges(5,3);
t454 = m(5) * t470 + t546 * mrSges(5,1) - t496 * mrSges(5,3) - t520 * t503 + t556 * t510 + t456;
t511 = t556 * mrSges(5,1) - t520 * mrSges(5,3);
t576 = -t562 * t464 + t566 * t465;
t455 = m(5) * t471 - t546 * mrSges(5,2) + t495 * mrSges(5,3) + t519 * t503 - t556 * t511 + t576;
t450 = t567 * t454 + t563 * t455;
t521 = -t543 * mrSges(4,1) + t544 * mrSges(4,2);
t525 = mrSges(4,2) * t583 + t543 * mrSges(4,3);
t448 = m(4) * t489 - t550 * mrSges(4,1) - t528 * mrSges(4,3) - t544 * t521 - t525 * t583 + t450;
t526 = -mrSges(4,1) * t583 - t544 * mrSges(4,3);
t577 = -t563 * t454 + t567 * t455;
t449 = m(4) * t490 + t550 * mrSges(4,2) + t527 * mrSges(4,3) + t543 * t521 + t526 * t583 + t577;
t578 = -t560 * t448 + t561 * t449;
t443 = m(3) * t523 - qJDD(2) * mrSges(3,2) + t550 * mrSges(3,3) - qJD(2) * t551 + t548 * t583 + t578;
t522 = -t568 * g(3) - t564 * t539;
t552 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t583;
t508 = -qJDD(2) * pkin(2) - t570 * qJ(3) + t547 * t584 + qJDD(3) - t522;
t491 = -t527 * pkin(3) - t542 * pkin(7) + t544 * t529 + t508;
t473 = -t495 * pkin(4) - t518 * pkin(8) + t520 * t512 + t491;
t575 = m(6) * t473 - t476 * mrSges(6,1) + t477 * mrSges(6,2) - t501 * t492 + t502 * t493;
t574 = m(5) * t491 - t495 * mrSges(5,1) + t496 * mrSges(5,2) - t519 * t510 + t520 * t511 + t575;
t572 = -m(4) * t508 + t527 * mrSges(4,1) - t528 * mrSges(4,2) + t543 * t525 - t544 * t526 - t574;
t460 = m(3) * t522 + qJDD(2) * mrSges(3,1) - t549 * mrSges(3,3) + qJD(2) * t552 - t548 * t584 + t572;
t579 = t568 * t443 - t564 * t460;
t437 = m(2) * t554 - t571 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t579;
t444 = t561 * t448 + t560 * t449;
t573 = -m(3) * t538 + t550 * mrSges(3,1) - t549 * mrSges(3,2) - t551 * t584 + t552 * t583 - t444;
t440 = m(2) * t553 + qJDD(1) * mrSges(2,1) - t571 * mrSges(2,2) + t573;
t585 = t565 * t437 + t569 * t440;
t438 = t564 * t443 + t568 * t460;
t580 = t569 * t437 - t565 * t440;
t537 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t564 + Ifges(3,4) * t568) * qJD(1);
t536 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t564 + Ifges(3,2) * t568) * qJD(1);
t535 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t564 + Ifges(3,6) * t568) * qJD(1);
t515 = Ifges(4,1) * t544 + Ifges(4,4) * t543 - Ifges(4,5) * t583;
t514 = Ifges(4,4) * t544 + Ifges(4,2) * t543 - Ifges(4,6) * t583;
t513 = Ifges(4,5) * t544 + Ifges(4,6) * t543 - Ifges(4,3) * t583;
t499 = Ifges(5,1) * t520 + Ifges(5,4) * t519 + Ifges(5,5) * t556;
t498 = Ifges(5,4) * t520 + Ifges(5,2) * t519 + Ifges(5,6) * t556;
t497 = Ifges(5,5) * t520 + Ifges(5,6) * t519 + Ifges(5,3) * t556;
t485 = Ifges(6,1) * t502 + Ifges(6,4) * t501 + Ifges(6,5) * t555;
t484 = Ifges(6,4) * t502 + Ifges(6,2) * t501 + Ifges(6,6) * t555;
t483 = Ifges(6,5) * t502 + Ifges(6,6) * t501 + Ifges(6,3) * t555;
t458 = mrSges(6,2) * t473 - mrSges(6,3) * t466 + Ifges(6,1) * t477 + Ifges(6,4) * t476 + Ifges(6,5) * t540 + t501 * t483 - t555 * t484;
t457 = -mrSges(6,1) * t473 + mrSges(6,3) * t467 + Ifges(6,4) * t477 + Ifges(6,2) * t476 + Ifges(6,6) * t540 - t502 * t483 + t555 * t485;
t446 = mrSges(5,2) * t491 - mrSges(5,3) * t470 + Ifges(5,1) * t496 + Ifges(5,4) * t495 + Ifges(5,5) * t546 - pkin(8) * t456 - t562 * t457 + t566 * t458 + t519 * t497 - t556 * t498;
t445 = -mrSges(5,1) * t491 + mrSges(5,3) * t471 + Ifges(5,4) * t496 + Ifges(5,2) * t495 + Ifges(5,6) * t546 - pkin(4) * t575 + pkin(8) * t576 + t566 * t457 + t562 * t458 - t520 * t497 + t556 * t499;
t434 = mrSges(4,2) * t508 - mrSges(4,3) * t489 + Ifges(4,1) * t528 + Ifges(4,4) * t527 - Ifges(4,5) * t550 - pkin(7) * t450 - t563 * t445 + t567 * t446 + t543 * t513 + t514 * t583;
t433 = -mrSges(4,1) * t508 + mrSges(4,3) * t490 + Ifges(4,4) * t528 + Ifges(4,2) * t527 - Ifges(4,6) * t550 - pkin(3) * t574 + pkin(7) * t577 + t567 * t445 + t563 * t446 - t544 * t513 - t515 * t583;
t432 = -pkin(3) * t450 + Ifges(3,6) * qJDD(2) - pkin(2) * t444 - pkin(4) * t456 + (Ifges(3,2) + Ifges(4,3)) * t550 - mrSges(6,1) * t466 + mrSges(6,2) * t467 - mrSges(5,1) * t470 + mrSges(5,2) * t471 - Ifges(6,6) * t476 - Ifges(6,5) * t477 - mrSges(4,1) * t489 + mrSges(4,2) * t490 - Ifges(5,6) * t495 - Ifges(5,5) * t496 - t535 * t584 + t501 * t485 - t502 * t484 + t519 * t499 - t520 * t498 + mrSges(3,3) * t523 - Ifges(4,6) * t527 - Ifges(4,5) * t528 + qJD(2) * t537 - mrSges(3,1) * t538 - Ifges(6,3) * t540 + t543 * t515 - t544 * t514 - Ifges(5,3) * t546 + Ifges(3,4) * t549;
t431 = mrSges(3,2) * t538 - mrSges(3,3) * t522 + Ifges(3,1) * t549 + Ifges(3,4) * t550 + Ifges(3,5) * qJDD(2) - qJ(3) * t444 - qJD(2) * t536 - t560 * t433 + t561 * t434 + t535 * t583;
t430 = Ifges(2,6) * qJDD(1) + t571 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t554 - Ifges(3,5) * t549 - Ifges(3,6) * t550 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t522 + mrSges(3,2) * t523 - t560 * t434 - t561 * t433 - pkin(2) * t572 - qJ(3) * t578 - pkin(1) * t438 + (-t564 * t536 + t568 * t537) * qJD(1);
t429 = -mrSges(2,2) * g(3) - mrSges(2,3) * t553 + Ifges(2,5) * qJDD(1) - t571 * Ifges(2,6) - pkin(6) * t438 + t568 * t431 - t564 * t432;
t1 = [-m(1) * g(1) + t580; -m(1) * g(2) + t585; (-m(1) - m(2)) * g(3) + t438; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t585 + t569 * t429 - t565 * t430; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t580 + t565 * t429 + t569 * t430; -mrSges(1,1) * g(2) + mrSges(2,1) * t553 + mrSges(1,2) * g(1) - mrSges(2,2) * t554 + Ifges(2,3) * qJDD(1) + pkin(1) * t573 + pkin(6) * t579 + t564 * t431 + t568 * t432;];
tauB = t1;
