% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP11_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP11_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP11_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:42
% EndTime: 2019-12-31 18:53:48
% DurationCPUTime: 3.61s
% Computational Cost: add. (31660->267), mult. (74120->323), div. (0->0), fcn. (51660->8), ass. (0->113)
t595 = Ifges(5,1) + Ifges(6,1);
t588 = Ifges(5,4) - Ifges(6,5);
t594 = -Ifges(5,5) - Ifges(6,4);
t593 = Ifges(5,2) + Ifges(6,3);
t586 = Ifges(5,6) - Ifges(6,6);
t592 = -Ifges(5,3) - Ifges(6,2);
t556 = qJD(1) ^ 2;
t548 = sin(pkin(8));
t549 = cos(pkin(8));
t551 = sin(qJ(3));
t553 = cos(qJ(3));
t561 = t548 * t551 - t549 * t553;
t531 = t561 * qJD(1);
t562 = t548 * t553 + t549 * t551;
t532 = t562 * qJD(1);
t576 = t532 * qJD(3);
t519 = -t561 * qJDD(1) - t576;
t591 = cos(qJ(4));
t590 = pkin(2) * t549;
t589 = -mrSges(5,3) - mrSges(6,2);
t585 = mrSges(3,2) * t548;
t547 = t549 ^ 2;
t584 = t547 * t556;
t552 = sin(qJ(1));
t554 = cos(qJ(1));
t537 = -t554 * g(1) - t552 * g(2);
t533 = -t556 * pkin(1) + qJDD(1) * qJ(2) + t537;
t575 = qJD(1) * qJD(2);
t572 = -t549 * g(3) - 0.2e1 * t548 * t575;
t506 = (-pkin(6) * qJDD(1) + t556 * t590 - t533) * t548 + t572;
t522 = -t548 * g(3) + (t533 + 0.2e1 * t575) * t549;
t574 = qJDD(1) * t549;
t507 = -pkin(2) * t584 + pkin(6) * t574 + t522;
t480 = t551 * t506 + t553 * t507;
t514 = t531 * mrSges(4,1) + t532 * mrSges(4,2);
t526 = qJD(3) * mrSges(4,1) - t532 * mrSges(4,3);
t517 = t531 * pkin(3) - t532 * pkin(7);
t555 = qJD(3) ^ 2;
t476 = -t555 * pkin(3) + qJDD(3) * pkin(7) - t531 * t517 + t480;
t546 = t548 ^ 2;
t536 = t552 * g(1) - t554 * g(2);
t566 = qJDD(2) - t536;
t518 = (-pkin(1) - t590) * qJDD(1) + (-qJ(2) + (-t546 - t547) * pkin(6)) * t556 + t566;
t577 = t531 * qJD(3);
t520 = t562 * qJDD(1) - t577;
t478 = (-t520 + t577) * pkin(7) + (-t519 + t576) * pkin(3) + t518;
t550 = sin(qJ(4));
t473 = t591 * t476 + t550 * t478;
t524 = t550 * qJD(3) + t591 * t532;
t491 = t524 * qJD(4) - t591 * qJDD(3) + t550 * t520;
t529 = qJD(4) + t531;
t502 = t529 * mrSges(5,1) - t524 * mrSges(5,3);
t516 = qJDD(4) - t519;
t523 = -t591 * qJD(3) + t550 * t532;
t495 = t523 * pkin(4) - t524 * qJ(5);
t528 = t529 ^ 2;
t469 = -t528 * pkin(4) + t516 * qJ(5) + 0.2e1 * qJD(5) * t529 - t523 * t495 + t473;
t503 = -t529 * mrSges(6,1) + t524 * mrSges(6,2);
t573 = m(6) * t469 + t516 * mrSges(6,3) + t529 * t503;
t496 = t523 * mrSges(6,1) - t524 * mrSges(6,3);
t579 = -t523 * mrSges(5,1) - t524 * mrSges(5,2) - t496;
t464 = m(5) * t473 - t516 * mrSges(5,2) + t589 * t491 - t529 * t502 + t579 * t523 + t573;
t472 = -t550 * t476 + t591 * t478;
t492 = -t523 * qJD(4) + t550 * qJDD(3) + t591 * t520;
t501 = -t529 * mrSges(5,2) - t523 * mrSges(5,3);
t470 = -t516 * pkin(4) - t528 * qJ(5) + t524 * t495 + qJDD(5) - t472;
t500 = -t523 * mrSges(6,2) + t529 * mrSges(6,3);
t567 = -m(6) * t470 + t516 * mrSges(6,1) + t529 * t500;
t466 = m(5) * t472 + t516 * mrSges(5,1) + t589 * t492 + t529 * t501 + t579 * t524 + t567;
t568 = t591 * t464 - t550 * t466;
t458 = m(4) * t480 - qJDD(3) * mrSges(4,2) + t519 * mrSges(4,3) - qJD(3) * t526 - t531 * t514 + t568;
t479 = t553 * t506 - t551 * t507;
t525 = -qJD(3) * mrSges(4,2) - t531 * mrSges(4,3);
t475 = -qJDD(3) * pkin(3) - t555 * pkin(7) + t532 * t517 - t479;
t471 = -0.2e1 * qJD(5) * t524 + (t523 * t529 - t492) * qJ(5) + (t524 * t529 + t491) * pkin(4) + t475;
t467 = m(6) * t471 + t491 * mrSges(6,1) - t492 * mrSges(6,3) + t523 * t500 - t524 * t503;
t557 = -m(5) * t475 - t491 * mrSges(5,1) - t492 * mrSges(5,2) - t523 * t501 - t524 * t502 - t467;
t461 = m(4) * t479 + qJDD(3) * mrSges(4,1) - t520 * mrSges(4,3) + qJD(3) * t525 - t532 * t514 + t557;
t451 = t551 * t458 + t553 * t461;
t521 = -t548 * t533 + t572;
t560 = mrSges(3,3) * qJDD(1) + t556 * (-mrSges(3,1) * t549 + t585);
t449 = m(3) * t521 - t560 * t548 + t451;
t569 = t553 * t458 - t551 * t461;
t450 = m(3) * t522 + t560 * t549 + t569;
t570 = -t548 * t449 + t549 * t450;
t443 = m(2) * t537 - t556 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t570;
t530 = -qJDD(1) * pkin(1) - t556 * qJ(2) + t566;
t459 = t550 * t464 + t591 * t466;
t559 = m(4) * t518 - t519 * mrSges(4,1) + t520 * mrSges(4,2) + t531 * t525 + t532 * t526 + t459;
t558 = -m(3) * t530 + mrSges(3,1) * t574 - t559 + (t546 * t556 + t584) * mrSges(3,3);
t453 = -t556 * mrSges(2,2) + m(2) * t536 + t558 + (mrSges(2,1) - t585) * qJDD(1);
t583 = t552 * t443 + t554 * t453;
t444 = t549 * t449 + t548 * t450;
t582 = t593 * t523 - t588 * t524 - t586 * t529;
t581 = t586 * t523 + t594 * t524 + t592 * t529;
t580 = -t588 * t523 + t595 * t524 - t594 * t529;
t563 = Ifges(3,5) * t548 + Ifges(3,6) * t549;
t578 = t556 * t563;
t571 = t554 * t443 - t552 * t453;
t565 = Ifges(3,1) * t548 + Ifges(3,4) * t549;
t564 = Ifges(3,4) * t548 + Ifges(3,2) * t549;
t510 = Ifges(4,1) * t532 - Ifges(4,4) * t531 + Ifges(4,5) * qJD(3);
t509 = Ifges(4,4) * t532 - Ifges(4,2) * t531 + Ifges(4,6) * qJD(3);
t508 = Ifges(4,5) * t532 - Ifges(4,6) * t531 + Ifges(4,3) * qJD(3);
t455 = mrSges(5,2) * t475 + mrSges(6,2) * t470 - mrSges(5,3) * t472 - mrSges(6,3) * t471 - qJ(5) * t467 - t588 * t491 + t595 * t492 - t516 * t594 + t581 * t523 + t582 * t529;
t454 = -mrSges(5,1) * t475 - mrSges(6,1) * t471 + mrSges(6,2) * t469 + mrSges(5,3) * t473 - pkin(4) * t467 - t593 * t491 + t588 * t492 + t586 * t516 + t581 * t524 + t580 * t529;
t445 = Ifges(4,4) * t520 + Ifges(4,2) * t519 + Ifges(4,6) * qJDD(3) - t532 * t508 + qJD(3) * t510 - mrSges(4,1) * t518 + mrSges(4,3) * t480 - mrSges(5,1) * t472 + mrSges(5,2) * t473 + mrSges(6,1) * t470 - mrSges(6,3) * t469 - pkin(4) * t567 - qJ(5) * t573 - pkin(3) * t459 + (pkin(4) * t496 + t582) * t524 + (qJ(5) * t496 - t580) * t523 + t592 * t516 + (pkin(4) * mrSges(6,2) + t594) * t492 + (qJ(5) * mrSges(6,2) + t586) * t491;
t440 = mrSges(4,2) * t518 - mrSges(4,3) * t479 + Ifges(4,1) * t520 + Ifges(4,4) * t519 + Ifges(4,5) * qJDD(3) - pkin(7) * t459 - qJD(3) * t509 - t550 * t454 + t591 * t455 - t531 * t508;
t439 = mrSges(3,2) * t530 - mrSges(3,3) * t521 - pkin(6) * t451 + t565 * qJDD(1) + t553 * t440 - t551 * t445 + t549 * t578;
t438 = -mrSges(3,1) * t530 + mrSges(3,3) * t522 - pkin(2) * t559 + pkin(6) * t569 + t564 * qJDD(1) + t551 * t440 + t553 * t445 - t548 * t578;
t437 = -pkin(1) * t444 + mrSges(2,1) * g(3) + mrSges(2,3) * t537 - pkin(2) * t451 - mrSges(3,1) * t521 + mrSges(3,2) * t522 - pkin(7) * t568 - t550 * t455 - t591 * t454 - pkin(3) * t557 - Ifges(4,5) * t520 - Ifges(4,6) * t519 - Ifges(4,3) * qJDD(3) - t532 * t509 - t531 * t510 - mrSges(4,1) * t479 + mrSges(4,2) * t480 + (Ifges(2,6) - t563) * qJDD(1) + (-t548 * t564 + t549 * t565 + Ifges(2,5)) * t556;
t436 = -mrSges(2,2) * g(3) - mrSges(2,3) * t536 + Ifges(2,5) * qJDD(1) - t556 * Ifges(2,6) - qJ(2) * t444 - t548 * t438 + t549 * t439;
t1 = [-m(1) * g(1) + t571; -m(1) * g(2) + t583; (-m(1) - m(2)) * g(3) + t444; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t583 + t554 * t436 - t552 * t437; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t571 + t552 * t436 + t554 * t437; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t536 - mrSges(2,2) * t537 + t548 * t439 + t549 * t438 + pkin(1) * (-qJDD(1) * t585 + t558) + qJ(2) * t570;];
tauB = t1;
