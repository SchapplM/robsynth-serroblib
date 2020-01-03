% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:53
% EndTime: 2019-12-31 21:23:04
% DurationCPUTime: 7.22s
% Computational Cost: add. (108302->311), mult. (225571->394), div. (0->0), fcn. (156123->10), ass. (0->120)
t567 = sin(qJ(1));
t571 = cos(qJ(1));
t556 = -t571 * g(1) - t567 * g(2);
t573 = qJD(1) ^ 2;
t541 = -t573 * pkin(1) + qJDD(1) * pkin(6) + t556;
t566 = sin(qJ(2));
t570 = cos(qJ(2));
t533 = -t566 * g(3) + t570 * t541;
t549 = (-mrSges(3,1) * t570 + mrSges(3,2) * t566) * qJD(1);
t584 = qJD(1) * qJD(2);
t559 = t566 * t584;
t552 = t570 * qJDD(1) - t559;
t586 = qJD(1) * t566;
t553 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t586;
t555 = t567 * g(1) - t571 * g(2);
t540 = -qJDD(1) * pkin(1) - t573 * pkin(6) - t555;
t583 = t570 * t584;
t551 = t566 * qJDD(1) + t583;
t508 = (-t551 - t583) * pkin(7) + (-t552 + t559) * pkin(2) + t540;
t550 = (-pkin(2) * t570 - pkin(7) * t566) * qJD(1);
t572 = qJD(2) ^ 2;
t585 = t570 * qJD(1);
t511 = -t572 * pkin(2) + qJDD(2) * pkin(7) + t550 * t585 + t533;
t565 = sin(qJ(3));
t569 = cos(qJ(3));
t492 = t569 * t508 - t565 * t511;
t547 = t569 * qJD(2) - t565 * t586;
t524 = t547 * qJD(3) + t565 * qJDD(2) + t569 * t551;
t546 = qJDD(3) - t552;
t548 = t565 * qJD(2) + t569 * t586;
t558 = qJD(3) - t585;
t482 = (t547 * t558 - t524) * qJ(4) + (t547 * t548 + t546) * pkin(3) + t492;
t493 = t565 * t508 + t569 * t511;
t523 = -t548 * qJD(3) + t569 * qJDD(2) - t565 * t551;
t530 = t558 * pkin(3) - t548 * qJ(4);
t545 = t547 ^ 2;
t484 = -t545 * pkin(3) + t523 * qJ(4) - t558 * t530 + t493;
t562 = sin(pkin(9));
t563 = cos(pkin(9));
t527 = t562 * t547 + t563 * t548;
t472 = -0.2e1 * qJD(4) * t527 + t563 * t482 - t562 * t484;
t501 = t562 * t523 + t563 * t524;
t526 = t563 * t547 - t562 * t548;
t470 = (t526 * t558 - t501) * pkin(8) + (t526 * t527 + t546) * pkin(4) + t472;
t473 = 0.2e1 * qJD(4) * t526 + t562 * t482 + t563 * t484;
t500 = t563 * t523 - t562 * t524;
t514 = t558 * pkin(4) - t527 * pkin(8);
t525 = t526 ^ 2;
t471 = -t525 * pkin(4) + t500 * pkin(8) - t558 * t514 + t473;
t564 = sin(qJ(5));
t568 = cos(qJ(5));
t468 = t568 * t470 - t564 * t471;
t503 = t568 * t526 - t564 * t527;
t479 = t503 * qJD(5) + t564 * t500 + t568 * t501;
t504 = t564 * t526 + t568 * t527;
t490 = -t503 * mrSges(6,1) + t504 * mrSges(6,2);
t557 = qJD(5) + t558;
t494 = -t557 * mrSges(6,2) + t503 * mrSges(6,3);
t542 = qJDD(5) + t546;
t466 = m(6) * t468 + t542 * mrSges(6,1) - t479 * mrSges(6,3) - t504 * t490 + t557 * t494;
t469 = t564 * t470 + t568 * t471;
t478 = -t504 * qJD(5) + t568 * t500 - t564 * t501;
t495 = t557 * mrSges(6,1) - t504 * mrSges(6,3);
t467 = m(6) * t469 - t542 * mrSges(6,2) + t478 * mrSges(6,3) + t503 * t490 - t557 * t495;
t458 = t568 * t466 + t564 * t467;
t505 = -t526 * mrSges(5,1) + t527 * mrSges(5,2);
t512 = -t558 * mrSges(5,2) + t526 * mrSges(5,3);
t456 = m(5) * t472 + t546 * mrSges(5,1) - t501 * mrSges(5,3) - t527 * t505 + t558 * t512 + t458;
t513 = t558 * mrSges(5,1) - t527 * mrSges(5,3);
t578 = -t564 * t466 + t568 * t467;
t457 = m(5) * t473 - t546 * mrSges(5,2) + t500 * mrSges(5,3) + t526 * t505 - t558 * t513 + t578;
t452 = t563 * t456 + t562 * t457;
t528 = -t547 * mrSges(4,1) + t548 * mrSges(4,2);
t529 = -t558 * mrSges(4,2) + t547 * mrSges(4,3);
t450 = m(4) * t492 + t546 * mrSges(4,1) - t524 * mrSges(4,3) - t548 * t528 + t558 * t529 + t452;
t531 = t558 * mrSges(4,1) - t548 * mrSges(4,3);
t579 = -t562 * t456 + t563 * t457;
t451 = m(4) * t493 - t546 * mrSges(4,2) + t523 * mrSges(4,3) + t547 * t528 - t558 * t531 + t579;
t580 = -t565 * t450 + t569 * t451;
t445 = m(3) * t533 - qJDD(2) * mrSges(3,2) + t552 * mrSges(3,3) - qJD(2) * t553 + t549 * t585 + t580;
t532 = -t570 * g(3) - t566 * t541;
t554 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t585;
t510 = -qJDD(2) * pkin(2) - t572 * pkin(7) + t550 * t586 - t532;
t491 = -t523 * pkin(3) - t545 * qJ(4) + t548 * t530 + qJDD(4) + t510;
t475 = -t500 * pkin(4) - t525 * pkin(8) + t527 * t514 + t491;
t577 = m(6) * t475 - t478 * mrSges(6,1) + t479 * mrSges(6,2) - t503 * t494 + t504 * t495;
t576 = m(5) * t491 - t500 * mrSges(5,1) + t501 * mrSges(5,2) - t526 * t512 + t527 * t513 + t577;
t574 = -m(4) * t510 + t523 * mrSges(4,1) - t524 * mrSges(4,2) + t547 * t529 - t548 * t531 - t576;
t462 = m(3) * t532 + qJDD(2) * mrSges(3,1) - t551 * mrSges(3,3) + qJD(2) * t554 - t549 * t586 + t574;
t581 = t570 * t445 - t566 * t462;
t439 = m(2) * t556 - t573 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t581;
t446 = t569 * t450 + t565 * t451;
t575 = -m(3) * t540 + t552 * mrSges(3,1) - t551 * mrSges(3,2) - t553 * t586 + t554 * t585 - t446;
t442 = m(2) * t555 + qJDD(1) * mrSges(2,1) - t573 * mrSges(2,2) + t575;
t587 = t567 * t439 + t571 * t442;
t440 = t566 * t445 + t570 * t462;
t582 = t571 * t439 - t567 * t442;
t539 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t566 + Ifges(3,4) * t570) * qJD(1);
t538 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t566 + Ifges(3,2) * t570) * qJD(1);
t537 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t566 + Ifges(3,6) * t570) * qJD(1);
t517 = Ifges(4,1) * t548 + Ifges(4,4) * t547 + Ifges(4,5) * t558;
t516 = Ifges(4,4) * t548 + Ifges(4,2) * t547 + Ifges(4,6) * t558;
t515 = Ifges(4,5) * t548 + Ifges(4,6) * t547 + Ifges(4,3) * t558;
t499 = Ifges(5,1) * t527 + Ifges(5,4) * t526 + Ifges(5,5) * t558;
t498 = Ifges(5,4) * t527 + Ifges(5,2) * t526 + Ifges(5,6) * t558;
t497 = Ifges(5,5) * t527 + Ifges(5,6) * t526 + Ifges(5,3) * t558;
t487 = Ifges(6,1) * t504 + Ifges(6,4) * t503 + Ifges(6,5) * t557;
t486 = Ifges(6,4) * t504 + Ifges(6,2) * t503 + Ifges(6,6) * t557;
t485 = Ifges(6,5) * t504 + Ifges(6,6) * t503 + Ifges(6,3) * t557;
t460 = mrSges(6,2) * t475 - mrSges(6,3) * t468 + Ifges(6,1) * t479 + Ifges(6,4) * t478 + Ifges(6,5) * t542 + t503 * t485 - t557 * t486;
t459 = -mrSges(6,1) * t475 + mrSges(6,3) * t469 + Ifges(6,4) * t479 + Ifges(6,2) * t478 + Ifges(6,6) * t542 - t504 * t485 + t557 * t487;
t448 = mrSges(5,2) * t491 - mrSges(5,3) * t472 + Ifges(5,1) * t501 + Ifges(5,4) * t500 + Ifges(5,5) * t546 - pkin(8) * t458 - t564 * t459 + t568 * t460 + t526 * t497 - t558 * t498;
t447 = -mrSges(5,1) * t491 + mrSges(5,3) * t473 + Ifges(5,4) * t501 + Ifges(5,2) * t500 + Ifges(5,6) * t546 - pkin(4) * t577 + pkin(8) * t578 + t568 * t459 + t564 * t460 - t527 * t497 + t558 * t499;
t436 = mrSges(4,2) * t510 - mrSges(4,3) * t492 + Ifges(4,1) * t524 + Ifges(4,4) * t523 + Ifges(4,5) * t546 - qJ(4) * t452 - t562 * t447 + t563 * t448 + t547 * t515 - t558 * t516;
t435 = -mrSges(4,1) * t510 + mrSges(4,3) * t493 + Ifges(4,4) * t524 + Ifges(4,2) * t523 + Ifges(4,6) * t546 - pkin(3) * t576 + qJ(4) * t579 + t563 * t447 + t562 * t448 - t548 * t515 + t558 * t517;
t434 = Ifges(3,6) * qJDD(2) + (-Ifges(4,3) - Ifges(5,3)) * t546 - t548 * t516 + Ifges(3,4) * t551 + Ifges(3,2) * t552 - mrSges(3,1) * t540 - Ifges(6,3) * t542 + t547 * t517 + t526 * t499 - t527 * t498 + mrSges(3,3) * t533 + qJD(2) * t539 - Ifges(4,6) * t523 - Ifges(4,5) * t524 + t503 * t487 - t504 * t486 - Ifges(5,6) * t500 - Ifges(5,5) * t501 - mrSges(4,1) * t492 + mrSges(4,2) * t493 - Ifges(6,5) * t479 - mrSges(5,1) * t472 + mrSges(5,2) * t473 - Ifges(6,6) * t478 + mrSges(6,2) * t469 - mrSges(6,1) * t468 - pkin(4) * t458 - pkin(3) * t452 - pkin(2) * t446 - t537 * t586;
t433 = mrSges(3,2) * t540 - mrSges(3,3) * t532 + Ifges(3,1) * t551 + Ifges(3,4) * t552 + Ifges(3,5) * qJDD(2) - pkin(7) * t446 - qJD(2) * t538 - t565 * t435 + t569 * t436 + t537 * t585;
t432 = Ifges(2,6) * qJDD(1) + t573 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t556 - Ifges(3,5) * t551 - Ifges(3,6) * t552 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t532 + mrSges(3,2) * t533 - t565 * t436 - t569 * t435 - pkin(2) * t574 - pkin(7) * t580 - pkin(1) * t440 + (-t566 * t538 + t570 * t539) * qJD(1);
t431 = -mrSges(2,2) * g(3) - mrSges(2,3) * t555 + Ifges(2,5) * qJDD(1) - t573 * Ifges(2,6) - pkin(6) * t440 + t570 * t433 - t566 * t434;
t1 = [-m(1) * g(1) + t582; -m(1) * g(2) + t587; (-m(1) - m(2)) * g(3) + t440; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t587 + t571 * t431 - t567 * t432; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t582 + t567 * t431 + t571 * t432; -mrSges(1,1) * g(2) + mrSges(2,1) * t555 + mrSges(1,2) * g(1) - mrSges(2,2) * t556 + Ifges(2,3) * qJDD(1) + pkin(1) * t575 + pkin(6) * t581 + t566 * t433 + t570 * t434;];
tauB = t1;
