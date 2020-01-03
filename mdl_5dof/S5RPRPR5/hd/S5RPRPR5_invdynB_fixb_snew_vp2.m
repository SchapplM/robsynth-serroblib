% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:42:15
% EndTime: 2020-01-03 11:42:27
% DurationCPUTime: 6.96s
% Computational Cost: add. (65060->288), mult. (167253->380), div. (0->0), fcn. (114048->10), ass. (0->121)
t583 = sin(qJ(1));
t586 = cos(qJ(1));
t561 = -t583 * g(2) + t586 * g(3);
t587 = qJD(1) ^ 2;
t618 = -t587 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t561;
t578 = sin(pkin(8));
t580 = cos(pkin(8));
t533 = -t580 * g(1) - t618 * t578;
t608 = t580 * qJD(1);
t565 = qJD(3) - t608;
t582 = sin(qJ(3));
t609 = t578 * qJD(1);
t603 = t582 * t609;
t548 = -t565 * mrSges(4,2) - mrSges(4,3) * t603;
t585 = cos(qJ(3));
t602 = t585 * t609;
t550 = t565 * mrSges(4,1) - mrSges(4,3) * t602;
t617 = -t548 * t582 - t550 * t585;
t616 = mrSges(3,2) * t578;
t613 = t578 ^ 2 * t587;
t534 = -t578 * g(1) + t618 * t580;
t557 = (-mrSges(3,1) * t580 + t616) * qJD(1);
t595 = -pkin(2) * t580 - pkin(6) * t578;
t559 = t595 * qJD(1);
t523 = t559 * t608 + t534;
t562 = -t586 * g(2) - t583 * g(3);
t591 = -t587 * qJ(2) + qJDD(2) - t562;
t535 = (-pkin(1) + t595) * qJDD(1) + t591;
t532 = t585 * t535;
t606 = qJD(1) * qJD(3);
t553 = (qJDD(1) * t585 - t582 * t606) * t578;
t605 = t580 * qJDD(1);
t564 = qJDD(3) - t605;
t502 = t564 * pkin(3) - t553 * qJ(4) + t532 + (-pkin(3) * t585 * t613 - qJ(4) * t565 * t609 - t523) * t582;
t512 = t585 * t523 + t582 * t535;
t549 = t565 * pkin(3) - qJ(4) * t602;
t552 = (-qJDD(1) * t582 - t585 * t606) * t578;
t604 = t582 ^ 2 * t613;
t503 = -pkin(3) * t604 + t552 * qJ(4) - t565 * t549 + t512;
t577 = sin(pkin(9));
t579 = cos(pkin(9));
t544 = (-t577 * t582 + t579 * t585) * t609;
t491 = -0.2e1 * qJD(4) * t544 + t579 * t502 - t577 * t503;
t527 = t577 * t552 + t579 * t553;
t543 = (-t577 * t585 - t579 * t582) * t609;
t489 = (t543 * t565 - t527) * pkin(7) + (t543 * t544 + t564) * pkin(4) + t491;
t492 = 0.2e1 * qJD(4) * t543 + t577 * t502 + t579 * t503;
t526 = t579 * t552 - t577 * t553;
t530 = t565 * pkin(4) - t544 * pkin(7);
t542 = t543 ^ 2;
t490 = -t542 * pkin(4) + t526 * pkin(7) - t565 * t530 + t492;
t581 = sin(qJ(5));
t584 = cos(qJ(5));
t487 = t584 * t489 - t581 * t490;
t520 = t584 * t543 - t581 * t544;
t498 = t520 * qJD(5) + t581 * t526 + t584 * t527;
t521 = t581 * t543 + t584 * t544;
t509 = -t520 * mrSges(6,1) + t521 * mrSges(6,2);
t563 = qJD(5) + t565;
t513 = -t563 * mrSges(6,2) + t520 * mrSges(6,3);
t560 = qJDD(5) + t564;
t485 = m(6) * t487 + t560 * mrSges(6,1) - t498 * mrSges(6,3) - t521 * t509 + t563 * t513;
t488 = t581 * t489 + t584 * t490;
t497 = -t521 * qJD(5) + t584 * t526 - t581 * t527;
t514 = t563 * mrSges(6,1) - t521 * mrSges(6,3);
t486 = m(6) * t488 - t560 * mrSges(6,2) + t497 * mrSges(6,3) + t520 * t509 - t563 * t514;
t477 = t584 * t485 + t581 * t486;
t524 = -t543 * mrSges(5,1) + t544 * mrSges(5,2);
t528 = -t565 * mrSges(5,2) + t543 * mrSges(5,3);
t475 = m(5) * t491 + t564 * mrSges(5,1) - t527 * mrSges(5,3) - t544 * t524 + t565 * t528 + t477;
t529 = t565 * mrSges(5,1) - t544 * mrSges(5,3);
t596 = -t581 * t485 + t584 * t486;
t476 = m(5) * t492 - t564 * mrSges(5,2) + t526 * mrSges(5,3) + t543 * t524 - t565 * t529 + t596;
t471 = t579 * t475 + t577 * t476;
t511 = -t582 * t523 + t532;
t551 = (mrSges(4,1) * t582 + mrSges(4,2) * t585) * t609;
t469 = m(4) * t511 + t564 * mrSges(4,1) - t553 * mrSges(4,3) + t565 * t548 - t551 * t602 + t471;
t597 = -t577 * t475 + t579 * t476;
t470 = m(4) * t512 - t564 * mrSges(4,2) + t552 * mrSges(4,3) - t565 * t550 - t551 * t603 + t597;
t598 = -t582 * t469 + t585 * t470;
t610 = qJDD(1) * mrSges(3,3);
t464 = m(3) * t534 + (qJD(1) * t557 + t610) * t580 + t598;
t522 = t559 * t609 - t533;
t510 = -t552 * pkin(3) - qJ(4) * t604 + t549 * t602 + qJDD(4) + t522;
t494 = -t526 * pkin(4) - t542 * pkin(7) + t544 * t530 + t510;
t592 = m(6) * t494 - t497 * mrSges(6,1) + t498 * mrSges(6,2) - t520 * t513 + t521 * t514;
t589 = m(5) * t510 - t526 * mrSges(5,1) + t527 * mrSges(5,2) - t543 * t528 + t544 * t529 + t592;
t588 = -m(4) * t522 + t552 * mrSges(4,1) - t553 * mrSges(4,2) - t589;
t484 = t588 + (-t610 + (-t557 + t617) * qJD(1)) * t578 + m(3) * t533;
t599 = t580 * t464 - t578 * t484;
t457 = m(2) * t561 - t587 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t599;
t465 = t585 * t469 + t582 * t470;
t555 = -qJDD(1) * pkin(1) + t591;
t590 = -m(3) * t555 + mrSges(3,1) * t605 - t465 + (t580 ^ 2 * t587 + t613) * mrSges(3,3);
t461 = m(2) * t562 - t587 * mrSges(2,2) + (mrSges(2,1) - t616) * qJDD(1) + t590;
t612 = t583 * t457 + t586 * t461;
t458 = t578 * t464 + t580 * t484;
t600 = -t586 * t457 + t583 * t461;
t594 = Ifges(3,1) * t578 + Ifges(3,4) * t580;
t593 = Ifges(3,5) * t578 + Ifges(3,6) * t580;
t558 = t593 * qJD(1);
t538 = Ifges(4,5) * t565 + (Ifges(4,1) * t585 - Ifges(4,4) * t582) * t609;
t537 = Ifges(4,6) * t565 + (Ifges(4,4) * t585 - Ifges(4,2) * t582) * t609;
t536 = Ifges(4,3) * t565 + (Ifges(4,5) * t585 - Ifges(4,6) * t582) * t609;
t517 = Ifges(5,1) * t544 + Ifges(5,4) * t543 + Ifges(5,5) * t565;
t516 = Ifges(5,4) * t544 + Ifges(5,2) * t543 + Ifges(5,6) * t565;
t515 = Ifges(5,5) * t544 + Ifges(5,6) * t543 + Ifges(5,3) * t565;
t506 = Ifges(6,1) * t521 + Ifges(6,4) * t520 + Ifges(6,5) * t563;
t505 = Ifges(6,4) * t521 + Ifges(6,2) * t520 + Ifges(6,6) * t563;
t504 = Ifges(6,5) * t521 + Ifges(6,6) * t520 + Ifges(6,3) * t563;
t479 = mrSges(6,2) * t494 - mrSges(6,3) * t487 + Ifges(6,1) * t498 + Ifges(6,4) * t497 + Ifges(6,5) * t560 + t520 * t504 - t563 * t505;
t478 = -mrSges(6,1) * t494 + mrSges(6,3) * t488 + Ifges(6,4) * t498 + Ifges(6,2) * t497 + Ifges(6,6) * t560 - t521 * t504 + t563 * t506;
t467 = mrSges(5,2) * t510 - mrSges(5,3) * t491 + Ifges(5,1) * t527 + Ifges(5,4) * t526 + Ifges(5,5) * t564 - pkin(7) * t477 - t581 * t478 + t584 * t479 + t543 * t515 - t565 * t516;
t466 = -mrSges(5,1) * t510 + mrSges(5,3) * t492 + Ifges(5,4) * t527 + Ifges(5,2) * t526 + Ifges(5,6) * t564 - pkin(4) * t592 + pkin(7) * t596 + t584 * t478 + t581 * t479 - t544 * t515 + t565 * t517;
t455 = mrSges(4,2) * t522 - mrSges(4,3) * t511 + Ifges(4,1) * t553 + Ifges(4,4) * t552 + Ifges(4,5) * t564 - qJ(4) * t471 - t577 * t466 + t579 * t467 - t536 * t603 - t565 * t537;
t454 = -mrSges(4,1) * t522 + mrSges(4,3) * t512 + Ifges(4,4) * t553 + Ifges(4,2) * t552 + Ifges(4,6) * t564 - pkin(3) * t589 + qJ(4) * t597 + t579 * t466 + t577 * t467 - t536 * t602 + t565 * t538;
t453 = (Ifges(3,4) * qJDD(1) + (-t537 * t585 - t538 * t582 - t558) * qJD(1)) * t578 + (-Ifges(4,3) - Ifges(5,3)) * t564 - Ifges(6,3) * t560 - Ifges(4,6) * t552 - Ifges(4,5) * t553 - mrSges(3,1) * t555 + t543 * t517 - t544 * t516 - Ifges(5,6) * t526 - Ifges(5,5) * t527 + mrSges(3,3) * t534 + t520 * t506 - t521 * t505 - mrSges(4,1) * t511 + mrSges(4,2) * t512 - Ifges(6,6) * t497 - Ifges(6,5) * t498 - mrSges(5,1) * t491 + mrSges(5,2) * t492 + mrSges(6,2) * t488 - mrSges(6,1) * t487 - pkin(4) * t477 - pkin(3) * t471 - pkin(2) * t465 + Ifges(3,2) * t605;
t452 = mrSges(3,2) * t555 - mrSges(3,3) * t533 - pkin(6) * t465 + qJDD(1) * t594 - t582 * t454 + t585 * t455 + t558 * t608;
t451 = t587 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t561 - mrSges(3,1) * t533 + mrSges(3,2) * t534 - t582 * t455 - t585 * t454 - pkin(2) * t588 - pkin(6) * t598 - pkin(1) * t458 + (Ifges(2,6) - t593) * qJDD(1) + (-pkin(2) * t617 * t578 + (-t578 * (Ifges(3,4) * t578 + Ifges(3,2) * t580) + t580 * t594) * qJD(1)) * qJD(1);
t450 = -mrSges(2,2) * g(1) - mrSges(2,3) * t562 + Ifges(2,5) * qJDD(1) - t587 * Ifges(2,6) - qJ(2) * t458 + t580 * t452 - t578 * t453;
t1 = [(-m(1) - m(2)) * g(1) + t458; -m(1) * g(2) + t612; -m(1) * g(3) + t600; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t562 - mrSges(2,2) * t561 + t578 * t452 + t580 * t453 + pkin(1) * (-qJDD(1) * t616 + t590) + qJ(2) * t599; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t600 + t583 * t450 + t586 * t451; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t612 - t586 * t450 + t583 * t451;];
tauB = t1;
