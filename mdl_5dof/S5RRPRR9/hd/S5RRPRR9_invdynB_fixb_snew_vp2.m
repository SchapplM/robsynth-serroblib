% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR9
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:32
% EndTime: 2019-12-31 20:20:40
% DurationCPUTime: 7.46s
% Computational Cost: add. (81341->313), mult. (184830->399), div. (0->0), fcn. (127734->10), ass. (0->121)
t573 = sin(qJ(2));
t577 = cos(qJ(2));
t592 = qJD(1) * qJD(2);
t559 = t573 * qJDD(1) + t577 * t592;
t574 = sin(qJ(1));
t578 = cos(qJ(1));
t565 = -t578 * g(1) - t574 * g(2);
t580 = qJD(1) ^ 2;
t554 = -t580 * pkin(1) + qJDD(1) * pkin(6) + t565;
t597 = t573 * t554;
t598 = pkin(2) * t580;
t516 = qJDD(2) * pkin(2) - t559 * qJ(3) - t597 + (qJ(3) * t592 + t573 * t598 - g(3)) * t577;
t539 = -t573 * g(3) + t577 * t554;
t560 = t577 * qJDD(1) - t573 * t592;
t595 = qJD(1) * t573;
t561 = qJD(2) * pkin(2) - qJ(3) * t595;
t568 = t577 ^ 2;
t517 = t560 * qJ(3) - qJD(2) * t561 - t568 * t598 + t539;
t569 = sin(pkin(9));
t570 = cos(pkin(9));
t548 = (t569 * t577 + t570 * t573) * qJD(1);
t499 = -0.2e1 * qJD(3) * t548 + t570 * t516 - t569 * t517;
t594 = qJD(1) * t577;
t547 = -t569 * t595 + t570 * t594;
t500 = 0.2e1 * qJD(3) * t547 + t569 * t516 + t570 * t517;
t528 = -t547 * mrSges(4,1) + t548 * mrSges(4,2);
t533 = -t569 * t559 + t570 * t560;
t541 = qJD(2) * mrSges(4,1) - t548 * mrSges(4,3);
t529 = -t547 * pkin(3) - t548 * pkin(7);
t579 = qJD(2) ^ 2;
t491 = -t579 * pkin(3) + qJDD(2) * pkin(7) + t547 * t529 + t500;
t564 = t574 * g(1) - t578 * g(2);
t585 = -qJDD(1) * pkin(1) - t564;
t520 = -t560 * pkin(2) + qJDD(3) + t561 * t595 + (-qJ(3) * t568 - pkin(6)) * t580 + t585;
t534 = t570 * t559 + t569 * t560;
t494 = (-qJD(2) * t547 - t534) * pkin(7) + (qJD(2) * t548 - t533) * pkin(3) + t520;
t572 = sin(qJ(4));
t576 = cos(qJ(4));
t484 = -t572 * t491 + t576 * t494;
t536 = t576 * qJD(2) - t572 * t548;
t510 = t536 * qJD(4) + t572 * qJDD(2) + t576 * t534;
t532 = qJDD(4) - t533;
t537 = t572 * qJD(2) + t576 * t548;
t546 = qJD(4) - t547;
t481 = (t536 * t546 - t510) * pkin(8) + (t536 * t537 + t532) * pkin(4) + t484;
t485 = t576 * t491 + t572 * t494;
t509 = -t537 * qJD(4) + t576 * qJDD(2) - t572 * t534;
t523 = t546 * pkin(4) - t537 * pkin(8);
t535 = t536 ^ 2;
t482 = -t535 * pkin(4) + t509 * pkin(8) - t546 * t523 + t485;
t571 = sin(qJ(5));
t575 = cos(qJ(5));
t479 = t575 * t481 - t571 * t482;
t514 = t575 * t536 - t571 * t537;
t488 = t514 * qJD(5) + t571 * t509 + t575 * t510;
t515 = t571 * t536 + t575 * t537;
t501 = -t514 * mrSges(6,1) + t515 * mrSges(6,2);
t542 = qJD(5) + t546;
t502 = -t542 * mrSges(6,2) + t514 * mrSges(6,3);
t530 = qJDD(5) + t532;
t477 = m(6) * t479 + t530 * mrSges(6,1) - t488 * mrSges(6,3) - t515 * t501 + t542 * t502;
t480 = t571 * t481 + t575 * t482;
t487 = -t515 * qJD(5) + t575 * t509 - t571 * t510;
t503 = t542 * mrSges(6,1) - t515 * mrSges(6,3);
t478 = m(6) * t480 - t530 * mrSges(6,2) + t487 * mrSges(6,3) + t514 * t501 - t542 * t503;
t469 = t575 * t477 + t571 * t478;
t518 = -t536 * mrSges(5,1) + t537 * mrSges(5,2);
t521 = -t546 * mrSges(5,2) + t536 * mrSges(5,3);
t467 = m(5) * t484 + t532 * mrSges(5,1) - t510 * mrSges(5,3) - t537 * t518 + t546 * t521 + t469;
t522 = t546 * mrSges(5,1) - t537 * mrSges(5,3);
t587 = -t571 * t477 + t575 * t478;
t468 = m(5) * t485 - t532 * mrSges(5,2) + t509 * mrSges(5,3) + t536 * t518 - t546 * t522 + t587;
t588 = -t572 * t467 + t576 * t468;
t462 = m(4) * t500 - qJDD(2) * mrSges(4,2) + t533 * mrSges(4,3) - qJD(2) * t541 + t547 * t528 + t588;
t540 = -qJD(2) * mrSges(4,2) + t547 * mrSges(4,3);
t490 = -qJDD(2) * pkin(3) - t579 * pkin(7) + t548 * t529 - t499;
t483 = -t509 * pkin(4) - t535 * pkin(8) + t537 * t523 + t490;
t584 = m(6) * t483 - t487 * mrSges(6,1) + t488 * mrSges(6,2) - t514 * t502 + t515 * t503;
t582 = -m(5) * t490 + t509 * mrSges(5,1) - t510 * mrSges(5,2) + t536 * t521 - t537 * t522 - t584;
t473 = m(4) * t499 + qJDD(2) * mrSges(4,1) - t534 * mrSges(4,3) + qJD(2) * t540 - t548 * t528 + t582;
t456 = t569 * t462 + t570 * t473;
t538 = -t577 * g(3) - t597;
t558 = (-mrSges(3,1) * t577 + mrSges(3,2) * t573) * qJD(1);
t563 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t594;
t453 = m(3) * t538 + qJDD(2) * mrSges(3,1) - t559 * mrSges(3,3) + qJD(2) * t563 - t558 * t595 + t456;
t562 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t595;
t589 = t570 * t462 - t569 * t473;
t454 = m(3) * t539 - qJDD(2) * mrSges(3,2) + t560 * mrSges(3,3) - qJD(2) * t562 + t558 * t594 + t589;
t590 = -t573 * t453 + t577 * t454;
t447 = m(2) * t565 - t580 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t590;
t553 = -t580 * pkin(6) + t585;
t463 = t576 * t467 + t572 * t468;
t583 = m(4) * t520 - t533 * mrSges(4,1) + t534 * mrSges(4,2) - t547 * t540 + t548 * t541 + t463;
t581 = -m(3) * t553 + t560 * mrSges(3,1) - t559 * mrSges(3,2) - t562 * t595 + t563 * t594 - t583;
t459 = m(2) * t564 + qJDD(1) * mrSges(2,1) - t580 * mrSges(2,2) + t581;
t596 = t574 * t447 + t578 * t459;
t448 = t577 * t453 + t573 * t454;
t591 = t578 * t447 - t574 * t459;
t551 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t573 + Ifges(3,4) * t577) * qJD(1);
t550 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t573 + Ifges(3,2) * t577) * qJD(1);
t549 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t573 + Ifges(3,6) * t577) * qJD(1);
t526 = Ifges(4,1) * t548 + Ifges(4,4) * t547 + Ifges(4,5) * qJD(2);
t525 = Ifges(4,4) * t548 + Ifges(4,2) * t547 + Ifges(4,6) * qJD(2);
t524 = Ifges(4,5) * t548 + Ifges(4,6) * t547 + Ifges(4,3) * qJD(2);
t506 = Ifges(5,1) * t537 + Ifges(5,4) * t536 + Ifges(5,5) * t546;
t505 = Ifges(5,4) * t537 + Ifges(5,2) * t536 + Ifges(5,6) * t546;
t504 = Ifges(5,5) * t537 + Ifges(5,6) * t536 + Ifges(5,3) * t546;
t497 = Ifges(6,1) * t515 + Ifges(6,4) * t514 + Ifges(6,5) * t542;
t496 = Ifges(6,4) * t515 + Ifges(6,2) * t514 + Ifges(6,6) * t542;
t495 = Ifges(6,5) * t515 + Ifges(6,6) * t514 + Ifges(6,3) * t542;
t471 = mrSges(6,2) * t483 - mrSges(6,3) * t479 + Ifges(6,1) * t488 + Ifges(6,4) * t487 + Ifges(6,5) * t530 + t514 * t495 - t542 * t496;
t470 = -mrSges(6,1) * t483 + mrSges(6,3) * t480 + Ifges(6,4) * t488 + Ifges(6,2) * t487 + Ifges(6,6) * t530 - t515 * t495 + t542 * t497;
t457 = mrSges(5,2) * t490 - mrSges(5,3) * t484 + Ifges(5,1) * t510 + Ifges(5,4) * t509 + Ifges(5,5) * t532 - pkin(8) * t469 - t571 * t470 + t575 * t471 + t536 * t504 - t546 * t505;
t455 = -mrSges(5,1) * t490 + mrSges(5,3) * t485 + Ifges(5,4) * t510 + Ifges(5,2) * t509 + Ifges(5,6) * t532 - pkin(4) * t584 + pkin(8) * t587 + t575 * t470 + t571 * t471 - t537 * t504 + t546 * t506;
t449 = Ifges(4,4) * t534 + Ifges(4,2) * t533 + Ifges(4,6) * qJDD(2) - t548 * t524 + qJD(2) * t526 - mrSges(4,1) * t520 + mrSges(4,3) * t500 - Ifges(5,5) * t510 - Ifges(5,6) * t509 - Ifges(5,3) * t532 - t537 * t505 + t536 * t506 - mrSges(5,1) * t484 + mrSges(5,2) * t485 - Ifges(6,5) * t488 - Ifges(6,6) * t487 - Ifges(6,3) * t530 - t515 * t496 + t514 * t497 - mrSges(6,1) * t479 + mrSges(6,2) * t480 - pkin(4) * t469 - pkin(3) * t463;
t444 = mrSges(4,2) * t520 - mrSges(4,3) * t499 + Ifges(4,1) * t534 + Ifges(4,4) * t533 + Ifges(4,5) * qJDD(2) - pkin(7) * t463 - qJD(2) * t525 - t572 * t455 + t576 * t457 + t547 * t524;
t443 = mrSges(3,2) * t553 - mrSges(3,3) * t538 + Ifges(3,1) * t559 + Ifges(3,4) * t560 + Ifges(3,5) * qJDD(2) - qJ(3) * t456 - qJD(2) * t550 + t570 * t444 - t569 * t449 + t549 * t594;
t442 = -pkin(1) * t448 + mrSges(2,3) * t565 - pkin(2) * t456 - Ifges(3,5) * t559 - Ifges(3,6) * t560 - mrSges(3,1) * t538 + mrSges(3,2) * t539 - t572 * t457 - t576 * t455 - pkin(3) * t582 - pkin(7) * t588 - Ifges(4,5) * t534 - Ifges(4,6) * t533 - mrSges(4,1) * t499 + mrSges(4,2) * t500 + mrSges(2,1) * g(3) - t548 * t525 + t547 * t526 + t580 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t573 * t550 + t577 * t551) * qJD(1);
t441 = -mrSges(3,1) * t553 + mrSges(3,3) * t539 + Ifges(3,4) * t559 + Ifges(3,2) * t560 + Ifges(3,6) * qJDD(2) - pkin(2) * t583 + qJ(3) * t589 + qJD(2) * t551 + t569 * t444 + t570 * t449 - t549 * t595;
t440 = -mrSges(2,2) * g(3) - mrSges(2,3) * t564 + Ifges(2,5) * qJDD(1) - t580 * Ifges(2,6) - pkin(6) * t448 - t573 * t441 + t577 * t443;
t1 = [-m(1) * g(1) + t591; -m(1) * g(2) + t596; (-m(1) - m(2)) * g(3) + t448; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t596 + t578 * t440 - t574 * t442; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t591 + t574 * t440 + t578 * t442; -mrSges(1,1) * g(2) + mrSges(2,1) * t564 + mrSges(1,2) * g(1) - mrSges(2,2) * t565 + Ifges(2,3) * qJDD(1) + pkin(1) * t581 + pkin(6) * t590 + t577 * t441 + t573 * t443;];
tauB = t1;
