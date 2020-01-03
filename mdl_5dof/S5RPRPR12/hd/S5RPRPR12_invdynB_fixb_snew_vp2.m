% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR12
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR12_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR12_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:49
% EndTime: 2019-12-31 18:29:56
% DurationCPUTime: 6.62s
% Computational Cost: add. (69114->291), mult. (169768->370), div. (0->0), fcn. (122932->10), ass. (0->123)
t572 = qJD(1) ^ 2;
t600 = cos(qJ(3));
t565 = cos(pkin(8));
t599 = pkin(2) * t565;
t563 = sin(pkin(8));
t598 = mrSges(3,2) * t563;
t561 = t565 ^ 2;
t597 = t561 * t572;
t568 = sin(qJ(1));
t570 = cos(qJ(1));
t551 = -t570 * g(1) - t568 * g(2);
t547 = -t572 * pkin(1) + qJDD(1) * qJ(2) + t551;
t592 = qJD(1) * qJD(2);
t588 = -t565 * g(3) - 0.2e1 * t563 * t592;
t514 = (-pkin(6) * qJDD(1) + t572 * t599 - t547) * t563 + t588;
t534 = -t563 * g(3) + (t547 + 0.2e1 * t592) * t565;
t590 = qJDD(1) * t565;
t521 = -pkin(2) * t597 + pkin(6) * t590 + t534;
t567 = sin(qJ(3));
t500 = t567 * t514 + t600 * t521;
t589 = t565 * t600;
t593 = t563 * qJD(1);
t545 = -qJD(1) * t589 + t567 * t593;
t577 = t600 * t563 + t565 * t567;
t546 = t577 * qJD(1);
t527 = t545 * mrSges(4,1) + t546 * mrSges(4,2);
t591 = qJDD(1) * t563;
t594 = t546 * qJD(3);
t531 = -qJDD(1) * t589 + t567 * t591 + t594;
t541 = qJD(3) * mrSges(4,1) - t546 * mrSges(4,3);
t526 = t545 * pkin(3) - t546 * qJ(4);
t571 = qJD(3) ^ 2;
t490 = -t571 * pkin(3) + qJDD(3) * qJ(4) - t545 * t526 + t500;
t560 = t563 ^ 2;
t550 = t568 * g(1) - t570 * g(2);
t582 = qJDD(2) - t550;
t530 = (-pkin(1) - t599) * qJDD(1) + (-qJ(2) + (-t560 - t561) * pkin(6)) * t572 + t582;
t595 = t545 * qJD(3);
t532 = t577 * qJDD(1) - t595;
t493 = (-t532 + t595) * qJ(4) + (t531 + t594) * pkin(3) + t530;
t562 = sin(pkin(9));
t564 = cos(pkin(9));
t539 = t562 * qJD(3) + t564 * t546;
t482 = -0.2e1 * qJD(4) * t539 - t562 * t490 + t564 * t493;
t520 = t562 * qJDD(3) + t564 * t532;
t538 = t564 * qJD(3) - t562 * t546;
t480 = (t538 * t545 - t520) * pkin(7) + (t538 * t539 + t531) * pkin(4) + t482;
t483 = 0.2e1 * qJD(4) * t538 + t564 * t490 + t562 * t493;
t518 = t545 * pkin(4) - t539 * pkin(7);
t519 = t564 * qJDD(3) - t562 * t532;
t537 = t538 ^ 2;
t481 = -t537 * pkin(4) + t519 * pkin(7) - t545 * t518 + t483;
t566 = sin(qJ(5));
t569 = cos(qJ(5));
t478 = t569 * t480 - t566 * t481;
t507 = t569 * t538 - t566 * t539;
t489 = t507 * qJD(5) + t566 * t519 + t569 * t520;
t508 = t566 * t538 + t569 * t539;
t498 = -t507 * mrSges(6,1) + t508 * mrSges(6,2);
t543 = qJD(5) + t545;
t501 = -t543 * mrSges(6,2) + t507 * mrSges(6,3);
t529 = qJDD(5) + t531;
t476 = m(6) * t478 + t529 * mrSges(6,1) - t489 * mrSges(6,3) - t508 * t498 + t543 * t501;
t479 = t566 * t480 + t569 * t481;
t488 = -t508 * qJD(5) + t569 * t519 - t566 * t520;
t502 = t543 * mrSges(6,1) - t508 * mrSges(6,3);
t477 = m(6) * t479 - t529 * mrSges(6,2) + t488 * mrSges(6,3) + t507 * t498 - t543 * t502;
t468 = t569 * t476 + t566 * t477;
t509 = -t538 * mrSges(5,1) + t539 * mrSges(5,2);
t516 = -t545 * mrSges(5,2) + t538 * mrSges(5,3);
t466 = m(5) * t482 + t531 * mrSges(5,1) - t520 * mrSges(5,3) - t539 * t509 + t545 * t516 + t468;
t517 = t545 * mrSges(5,1) - t539 * mrSges(5,3);
t583 = -t566 * t476 + t569 * t477;
t467 = m(5) * t483 - t531 * mrSges(5,2) + t519 * mrSges(5,3) + t538 * t509 - t545 * t517 + t583;
t584 = -t562 * t466 + t564 * t467;
t461 = m(4) * t500 - qJDD(3) * mrSges(4,2) - t531 * mrSges(4,3) - qJD(3) * t541 - t545 * t527 + t584;
t499 = t600 * t514 - t567 * t521;
t540 = -qJD(3) * mrSges(4,2) - t545 * mrSges(4,3);
t487 = -qJDD(3) * pkin(3) - t571 * qJ(4) + t546 * t526 + qJDD(4) - t499;
t484 = -t519 * pkin(4) - t537 * pkin(7) + t539 * t518 + t487;
t576 = m(6) * t484 - t488 * mrSges(6,1) + t489 * mrSges(6,2) - t507 * t501 + t508 * t502;
t573 = -m(5) * t487 + t519 * mrSges(5,1) - t520 * mrSges(5,2) + t538 * t516 - t539 * t517 - t576;
t472 = m(4) * t499 + qJDD(3) * mrSges(4,1) - t532 * mrSges(4,3) + qJD(3) * t540 - t546 * t527 + t573;
t454 = t567 * t461 + t600 * t472;
t533 = -t563 * t547 + t588;
t578 = mrSges(3,3) * qJDD(1) + t572 * (-mrSges(3,1) * t565 + t598);
t452 = m(3) * t533 - t578 * t563 + t454;
t585 = t600 * t461 - t567 * t472;
t453 = m(3) * t534 + t578 * t565 + t585;
t586 = -t563 * t452 + t565 * t453;
t446 = m(2) * t551 - t572 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t586;
t544 = -qJDD(1) * pkin(1) - t572 * qJ(2) + t582;
t462 = t564 * t466 + t562 * t467;
t575 = m(4) * t530 + t531 * mrSges(4,1) + t532 * mrSges(4,2) + t545 * t540 + t546 * t541 + t462;
t574 = -m(3) * t544 + mrSges(3,1) * t590 - t575 + (t560 * t572 + t597) * mrSges(3,3);
t458 = t574 + (mrSges(2,1) - t598) * qJDD(1) + m(2) * t550 - t572 * mrSges(2,2);
t596 = t568 * t446 + t570 * t458;
t447 = t565 * t452 + t563 * t453;
t587 = t570 * t446 - t568 * t458;
t581 = Ifges(3,1) * t563 + Ifges(3,4) * t565;
t580 = Ifges(3,4) * t563 + Ifges(3,2) * t565;
t579 = Ifges(3,5) * t563 + Ifges(3,6) * t565;
t549 = t579 * qJD(1);
t524 = Ifges(4,1) * t546 - Ifges(4,4) * t545 + Ifges(4,5) * qJD(3);
t523 = Ifges(4,4) * t546 - Ifges(4,2) * t545 + Ifges(4,6) * qJD(3);
t522 = Ifges(4,5) * t546 - Ifges(4,6) * t545 + Ifges(4,3) * qJD(3);
t505 = Ifges(5,1) * t539 + Ifges(5,4) * t538 + Ifges(5,5) * t545;
t504 = Ifges(5,4) * t539 + Ifges(5,2) * t538 + Ifges(5,6) * t545;
t503 = Ifges(5,5) * t539 + Ifges(5,6) * t538 + Ifges(5,3) * t545;
t496 = Ifges(6,1) * t508 + Ifges(6,4) * t507 + Ifges(6,5) * t543;
t495 = Ifges(6,4) * t508 + Ifges(6,2) * t507 + Ifges(6,6) * t543;
t494 = Ifges(6,5) * t508 + Ifges(6,6) * t507 + Ifges(6,3) * t543;
t470 = mrSges(6,2) * t484 - mrSges(6,3) * t478 + Ifges(6,1) * t489 + Ifges(6,4) * t488 + Ifges(6,5) * t529 + t507 * t494 - t543 * t495;
t469 = -mrSges(6,1) * t484 + mrSges(6,3) * t479 + Ifges(6,4) * t489 + Ifges(6,2) * t488 + Ifges(6,6) * t529 - t508 * t494 + t543 * t496;
t456 = mrSges(5,2) * t487 - mrSges(5,3) * t482 + Ifges(5,1) * t520 + Ifges(5,4) * t519 + Ifges(5,5) * t531 - pkin(7) * t468 - t566 * t469 + t569 * t470 + t538 * t503 - t545 * t504;
t455 = -mrSges(5,1) * t487 + mrSges(5,3) * t483 + Ifges(5,4) * t520 + Ifges(5,2) * t519 + Ifges(5,6) * t531 - pkin(4) * t576 + pkin(7) * t583 + t569 * t469 + t566 * t470 - t539 * t503 + t545 * t505;
t448 = Ifges(4,4) * t532 + Ifges(4,6) * qJDD(3) - t546 * t522 + qJD(3) * t524 - mrSges(4,1) * t530 + mrSges(4,3) * t500 - Ifges(5,5) * t520 - Ifges(5,6) * t519 - t539 * t504 + t538 * t505 - mrSges(5,1) * t482 + mrSges(5,2) * t483 - Ifges(6,5) * t489 - Ifges(6,6) * t488 - Ifges(6,3) * t529 - t508 * t495 + t507 * t496 - mrSges(6,1) * t478 + mrSges(6,2) * t479 - pkin(4) * t468 - pkin(3) * t462 + (-Ifges(4,2) - Ifges(5,3)) * t531;
t443 = mrSges(4,2) * t530 - mrSges(4,3) * t499 + Ifges(4,1) * t532 - Ifges(4,4) * t531 + Ifges(4,5) * qJDD(3) - qJ(4) * t462 - qJD(3) * t523 - t562 * t455 + t564 * t456 - t545 * t522;
t442 = t565 * qJD(1) * t549 + mrSges(3,2) * t544 - mrSges(3,3) * t533 - pkin(6) * t454 + t581 * qJDD(1) + t600 * t443 - t567 * t448;
t441 = mrSges(2,1) * g(3) - pkin(1) * t447 + mrSges(2,3) * t551 - pkin(2) * t454 - mrSges(3,1) * t533 + mrSges(3,2) * t534 - t562 * t456 - t564 * t455 - pkin(3) * t573 - qJ(4) * t584 - Ifges(4,5) * t532 + Ifges(4,6) * t531 - Ifges(4,3) * qJDD(3) - t546 * t523 - t545 * t524 - mrSges(4,1) * t499 + mrSges(4,2) * t500 + (Ifges(2,6) - t579) * qJDD(1) + (-t563 * t580 + t565 * t581 + Ifges(2,5)) * t572;
t440 = -mrSges(3,1) * t544 + mrSges(3,3) * t534 - pkin(2) * t575 + pkin(6) * t585 + t580 * qJDD(1) + t567 * t443 + t600 * t448 - t549 * t593;
t439 = -mrSges(2,2) * g(3) - mrSges(2,3) * t550 + Ifges(2,5) * qJDD(1) - t572 * Ifges(2,6) - qJ(2) * t447 - t563 * t440 + t565 * t442;
t1 = [-m(1) * g(1) + t587; -m(1) * g(2) + t596; (-m(1) - m(2)) * g(3) + t447; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t596 + t570 * t439 - t568 * t441; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t587 + t568 * t439 + t570 * t441; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t550 - mrSges(2,2) * t551 + t563 * t442 + t565 * t440 + pkin(1) * (-mrSges(3,2) * t591 + t574) + qJ(2) * t586;];
tauB = t1;
