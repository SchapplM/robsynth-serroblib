% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR2
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:49
% EndTime: 2019-12-05 18:27:59
% DurationCPUTime: 8.90s
% Computational Cost: add. (123573->313), mult. (294410->400), div. (0->0), fcn. (209494->10), ass. (0->121)
t589 = qJD(1) ^ 2;
t605 = pkin(2) * t589;
t584 = sin(qJ(1));
t588 = cos(qJ(1));
t572 = -g(1) * t588 - g(2) * t584;
t561 = -pkin(1) * t589 + qJDD(1) * pkin(6) + t572;
t583 = sin(qJ(2));
t604 = t561 * t583;
t587 = cos(qJ(2));
t600 = qJD(1) * qJD(2);
t566 = qJDD(1) * t583 + t587 * t600;
t528 = qJDD(2) * pkin(2) - qJ(3) * t566 - t604 + (qJ(3) * t600 + t583 * t605 - g(3)) * t587;
t547 = -g(3) * t583 + t587 * t561;
t567 = qJDD(1) * t587 - t583 * t600;
t602 = qJD(1) * t583;
t568 = qJD(2) * pkin(2) - qJ(3) * t602;
t578 = t587 ^ 2;
t529 = qJ(3) * t567 - qJD(2) * t568 - t578 * t605 + t547;
t579 = sin(pkin(9));
t580 = cos(pkin(9));
t556 = (t579 * t587 + t580 * t583) * qJD(1);
t509 = -0.2e1 * qJD(3) * t556 + t580 * t528 - t529 * t579;
t545 = t566 * t580 + t567 * t579;
t555 = (-t579 * t583 + t580 * t587) * qJD(1);
t498 = (qJD(2) * t555 - t545) * pkin(7) + (t555 * t556 + qJDD(2)) * pkin(3) + t509;
t510 = 0.2e1 * qJD(3) * t555 + t579 * t528 + t580 * t529;
t544 = -t566 * t579 + t567 * t580;
t550 = qJD(2) * pkin(3) - pkin(7) * t556;
t554 = t555 ^ 2;
t500 = -pkin(3) * t554 + pkin(7) * t544 - qJD(2) * t550 + t510;
t582 = sin(qJ(4));
t586 = cos(qJ(4));
t488 = t586 * t498 - t500 * t582;
t538 = t555 * t586 - t556 * t582;
t514 = qJD(4) * t538 + t544 * t582 + t545 * t586;
t539 = t555 * t582 + t556 * t586;
t576 = qJDD(2) + qJDD(4);
t577 = qJD(2) + qJD(4);
t486 = (t538 * t577 - t514) * pkin(8) + (t538 * t539 + t576) * pkin(4) + t488;
t489 = t582 * t498 + t586 * t500;
t513 = -qJD(4) * t539 + t544 * t586 - t545 * t582;
t533 = pkin(4) * t577 - pkin(8) * t539;
t534 = t538 ^ 2;
t487 = -pkin(4) * t534 + pkin(8) * t513 - t533 * t577 + t489;
t581 = sin(qJ(5));
t585 = cos(qJ(5));
t484 = t486 * t585 - t487 * t581;
t522 = t538 * t585 - t539 * t581;
t495 = qJD(5) * t522 + t513 * t581 + t514 * t585;
t523 = t538 * t581 + t539 * t585;
t506 = -mrSges(6,1) * t522 + mrSges(6,2) * t523;
t574 = qJD(5) + t577;
t515 = -mrSges(6,2) * t574 + mrSges(6,3) * t522;
t573 = qJDD(5) + t576;
t482 = m(6) * t484 + mrSges(6,1) * t573 - mrSges(6,3) * t495 - t506 * t523 + t515 * t574;
t485 = t486 * t581 + t487 * t585;
t494 = -qJD(5) * t523 + t513 * t585 - t514 * t581;
t516 = mrSges(6,1) * t574 - mrSges(6,3) * t523;
t483 = m(6) * t485 - mrSges(6,2) * t573 + mrSges(6,3) * t494 + t506 * t522 - t516 * t574;
t474 = t585 * t482 + t581 * t483;
t524 = -mrSges(5,1) * t538 + mrSges(5,2) * t539;
t531 = -mrSges(5,2) * t577 + mrSges(5,3) * t538;
t472 = m(5) * t488 + mrSges(5,1) * t576 - mrSges(5,3) * t514 - t524 * t539 + t531 * t577 + t474;
t532 = mrSges(5,1) * t577 - mrSges(5,3) * t539;
t595 = -t482 * t581 + t585 * t483;
t473 = m(5) * t489 - mrSges(5,2) * t576 + mrSges(5,3) * t513 + t524 * t538 - t532 * t577 + t595;
t468 = t586 * t472 + t582 * t473;
t542 = -mrSges(4,1) * t555 + mrSges(4,2) * t556;
t548 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t555;
t466 = m(4) * t509 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t545 + qJD(2) * t548 - t542 * t556 + t468;
t549 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t556;
t596 = -t472 * t582 + t586 * t473;
t467 = m(4) * t510 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t544 - qJD(2) * t549 + t542 * t555 + t596;
t460 = t580 * t466 + t579 * t467;
t546 = -g(3) * t587 - t604;
t565 = (-mrSges(3,1) * t587 + mrSges(3,2) * t583) * qJD(1);
t601 = qJD(1) * t587;
t570 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t601;
t458 = m(3) * t546 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t566 + qJD(2) * t570 - t565 * t602 + t460;
t569 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t602;
t597 = -t579 * t466 + t580 * t467;
t459 = m(3) * t547 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t567 - qJD(2) * t569 + t565 * t601 + t597;
t598 = -t458 * t583 + t587 * t459;
t451 = m(2) * t572 - mrSges(2,1) * t589 - qJDD(1) * mrSges(2,2) + t598;
t571 = g(1) * t584 - t588 * g(2);
t593 = -qJDD(1) * pkin(1) - t571;
t560 = -pkin(6) * t589 + t593;
t530 = -pkin(2) * t567 + qJDD(3) + t568 * t602 + (-qJ(3) * t578 - pkin(6)) * t589 + t593;
t508 = -pkin(3) * t544 - pkin(7) * t554 + t556 * t550 + t530;
t491 = -pkin(4) * t513 - pkin(8) * t534 + t533 * t539 + t508;
t594 = m(6) * t491 - t494 * mrSges(6,1) + t495 * mrSges(6,2) - t522 * t515 + t523 * t516;
t592 = m(5) * t508 - t513 * mrSges(5,1) + t514 * mrSges(5,2) - t538 * t531 + t539 * t532 + t594;
t591 = m(4) * t530 - t544 * mrSges(4,1) + t545 * mrSges(4,2) - t555 * t548 + t556 * t549 + t592;
t590 = -m(3) * t560 + t567 * mrSges(3,1) - t566 * mrSges(3,2) - t569 * t602 + t570 * t601 - t591;
t478 = m(2) * t571 + qJDD(1) * mrSges(2,1) - t589 * mrSges(2,2) + t590;
t603 = t584 * t451 + t588 * t478;
t452 = t587 * t458 + t583 * t459;
t599 = t588 * t451 - t478 * t584;
t559 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t583 + Ifges(3,4) * t587) * qJD(1);
t558 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t583 + Ifges(3,2) * t587) * qJD(1);
t557 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t583 + Ifges(3,6) * t587) * qJD(1);
t537 = Ifges(4,1) * t556 + Ifges(4,4) * t555 + Ifges(4,5) * qJD(2);
t536 = Ifges(4,4) * t556 + Ifges(4,2) * t555 + Ifges(4,6) * qJD(2);
t535 = Ifges(4,5) * t556 + Ifges(4,6) * t555 + Ifges(4,3) * qJD(2);
t519 = Ifges(5,1) * t539 + Ifges(5,4) * t538 + Ifges(5,5) * t577;
t518 = Ifges(5,4) * t539 + Ifges(5,2) * t538 + Ifges(5,6) * t577;
t517 = Ifges(5,5) * t539 + Ifges(5,6) * t538 + Ifges(5,3) * t577;
t503 = Ifges(6,1) * t523 + Ifges(6,4) * t522 + Ifges(6,5) * t574;
t502 = Ifges(6,4) * t523 + Ifges(6,2) * t522 + Ifges(6,6) * t574;
t501 = Ifges(6,5) * t523 + Ifges(6,6) * t522 + Ifges(6,3) * t574;
t476 = mrSges(6,2) * t491 - mrSges(6,3) * t484 + Ifges(6,1) * t495 + Ifges(6,4) * t494 + Ifges(6,5) * t573 + t501 * t522 - t502 * t574;
t475 = -mrSges(6,1) * t491 + mrSges(6,3) * t485 + Ifges(6,4) * t495 + Ifges(6,2) * t494 + Ifges(6,6) * t573 - t501 * t523 + t503 * t574;
t462 = mrSges(5,2) * t508 - mrSges(5,3) * t488 + Ifges(5,1) * t514 + Ifges(5,4) * t513 + Ifges(5,5) * t576 - pkin(8) * t474 - t475 * t581 + t476 * t585 + t517 * t538 - t518 * t577;
t461 = -mrSges(5,1) * t508 + mrSges(5,3) * t489 + Ifges(5,4) * t514 + Ifges(5,2) * t513 + Ifges(5,6) * t576 - pkin(4) * t594 + pkin(8) * t595 + t585 * t475 + t581 * t476 - t539 * t517 + t577 * t519;
t454 = mrSges(4,2) * t530 - mrSges(4,3) * t509 + Ifges(4,1) * t545 + Ifges(4,4) * t544 + Ifges(4,5) * qJDD(2) - pkin(7) * t468 - qJD(2) * t536 - t461 * t582 + t462 * t586 + t535 * t555;
t453 = -mrSges(4,1) * t530 + mrSges(4,3) * t510 + Ifges(4,4) * t545 + Ifges(4,2) * t544 + Ifges(4,6) * qJDD(2) - pkin(3) * t592 + pkin(7) * t596 + qJD(2) * t537 + t586 * t461 + t582 * t462 - t556 * t535;
t448 = Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) - pkin(1) * t452 + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t558 * t583 + t559 * t587) * qJD(1) + t589 * Ifges(2,5) - Ifges(6,3) * t573 - Ifges(5,3) * t576 - Ifges(3,5) * t566 - Ifges(3,6) * t567 + mrSges(2,3) * t572 + t555 * t537 - t556 * t536 - Ifges(4,6) * t544 - Ifges(4,5) * t545 - mrSges(3,1) * t546 + mrSges(3,2) * t547 + t538 * t519 - t539 * t518 + t522 * t503 - t523 * t502 - Ifges(5,6) * t513 - Ifges(5,5) * t514 - mrSges(4,1) * t509 + mrSges(4,2) * t510 - Ifges(6,6) * t494 - Ifges(6,5) * t495 - mrSges(5,1) * t488 + mrSges(5,2) * t489 + mrSges(6,2) * t485 - mrSges(6,1) * t484 - pkin(4) * t474 - pkin(3) * t468 - pkin(2) * t460;
t447 = mrSges(3,2) * t560 - mrSges(3,3) * t546 + Ifges(3,1) * t566 + Ifges(3,4) * t567 + Ifges(3,5) * qJDD(2) - qJ(3) * t460 - qJD(2) * t558 - t453 * t579 + t454 * t580 + t557 * t601;
t446 = -mrSges(3,1) * t560 + mrSges(3,3) * t547 + Ifges(3,4) * t566 + Ifges(3,2) * t567 + Ifges(3,6) * qJDD(2) - pkin(2) * t591 + qJ(3) * t597 + qJD(2) * t559 + t580 * t453 + t579 * t454 - t557 * t602;
t445 = -mrSges(2,2) * g(3) - mrSges(2,3) * t571 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t589 - pkin(6) * t452 - t446 * t583 + t447 * t587;
t1 = [-m(1) * g(1) + t599; -m(1) * g(2) + t603; (-m(1) - m(2)) * g(3) + t452; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t603 + t588 * t445 - t584 * t448; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t599 + t584 * t445 + t588 * t448; -mrSges(1,1) * g(2) + mrSges(2,1) * t571 + mrSges(1,2) * g(1) - mrSges(2,2) * t572 + Ifges(2,3) * qJDD(1) + pkin(1) * t590 + pkin(6) * t598 + t587 * t446 + t583 * t447;];
tauB = t1;
