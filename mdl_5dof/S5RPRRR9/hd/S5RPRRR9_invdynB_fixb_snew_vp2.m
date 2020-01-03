% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:44
% EndTime: 2019-12-31 19:07:51
% DurationCPUTime: 7.51s
% Computational Cost: add. (80518->291), mult. (193932->366), div. (0->0), fcn. (146354->10), ass. (0->122)
t573 = qJD(1) ^ 2;
t564 = cos(pkin(9));
t598 = pkin(2) * t564;
t563 = sin(pkin(9));
t597 = mrSges(3,2) * t563;
t561 = t564 ^ 2;
t596 = t561 * t573;
t568 = sin(qJ(1));
t572 = cos(qJ(1));
t551 = -t572 * g(1) - t568 * g(2);
t547 = -t573 * pkin(1) + qJDD(1) * qJ(2) + t551;
t592 = qJD(1) * qJD(2);
t590 = -t564 * g(3) - 0.2e1 * t563 * t592;
t522 = (-pkin(6) * qJDD(1) + t573 * t598 - t547) * t563 + t590;
t538 = -t563 * g(3) + (t547 + 0.2e1 * t592) * t564;
t591 = qJDD(1) * t564;
t523 = -pkin(2) * t596 + pkin(6) * t591 + t538;
t567 = sin(qJ(3));
t571 = cos(qJ(3));
t504 = t571 * t522 - t567 * t523;
t580 = t563 * t571 + t564 * t567;
t579 = -t563 * t567 + t564 * t571;
t545 = t579 * qJD(1);
t593 = t545 * qJD(3);
t536 = t580 * qJDD(1) + t593;
t546 = t580 * qJD(1);
t486 = (-t536 + t593) * pkin(7) + (t545 * t546 + qJDD(3)) * pkin(3) + t504;
t505 = t567 * t522 + t571 * t523;
t535 = -t546 * qJD(3) + t579 * qJDD(1);
t541 = qJD(3) * pkin(3) - t546 * pkin(7);
t544 = t545 ^ 2;
t491 = -t544 * pkin(3) + t535 * pkin(7) - qJD(3) * t541 + t505;
t566 = sin(qJ(4));
t570 = cos(qJ(4));
t484 = t566 * t486 + t570 * t491;
t529 = t566 * t545 + t570 * t546;
t501 = -t529 * qJD(4) + t570 * t535 - t566 * t536;
t528 = t570 * t545 - t566 * t546;
t513 = -t528 * mrSges(5,1) + t529 * mrSges(5,2);
t562 = qJD(3) + qJD(4);
t520 = t562 * mrSges(5,1) - t529 * mrSges(5,3);
t559 = qJDD(3) + qJDD(4);
t514 = -t528 * pkin(4) - t529 * pkin(8);
t558 = t562 ^ 2;
t481 = -t558 * pkin(4) + t559 * pkin(8) + t528 * t514 + t484;
t560 = t563 ^ 2;
t550 = t568 * g(1) - t572 * g(2);
t584 = qJDD(2) - t550;
t534 = (-pkin(1) - t598) * qJDD(1) + (-qJ(2) + (-t560 - t561) * pkin(6)) * t573 + t584;
t497 = -t535 * pkin(3) - t544 * pkin(7) + t546 * t541 + t534;
t502 = t528 * qJD(4) + t566 * t535 + t570 * t536;
t482 = (-t528 * t562 - t502) * pkin(8) + (t529 * t562 - t501) * pkin(4) + t497;
t565 = sin(qJ(5));
t569 = cos(qJ(5));
t478 = -t565 * t481 + t569 * t482;
t515 = -t565 * t529 + t569 * t562;
t489 = t515 * qJD(5) + t569 * t502 + t565 * t559;
t500 = qJDD(5) - t501;
t516 = t569 * t529 + t565 * t562;
t503 = -t515 * mrSges(6,1) + t516 * mrSges(6,2);
t524 = qJD(5) - t528;
t506 = -t524 * mrSges(6,2) + t515 * mrSges(6,3);
t476 = m(6) * t478 + t500 * mrSges(6,1) - t489 * mrSges(6,3) - t516 * t503 + t524 * t506;
t479 = t569 * t481 + t565 * t482;
t488 = -t516 * qJD(5) - t565 * t502 + t569 * t559;
t507 = t524 * mrSges(6,1) - t516 * mrSges(6,3);
t477 = m(6) * t479 - t500 * mrSges(6,2) + t488 * mrSges(6,3) + t515 * t503 - t524 * t507;
t585 = -t565 * t476 + t569 * t477;
t467 = m(5) * t484 - t559 * mrSges(5,2) + t501 * mrSges(5,3) + t528 * t513 - t562 * t520 + t585;
t483 = t570 * t486 - t566 * t491;
t519 = -t562 * mrSges(5,2) + t528 * mrSges(5,3);
t480 = -t559 * pkin(4) - t558 * pkin(8) + t529 * t514 - t483;
t576 = -m(6) * t480 + t488 * mrSges(6,1) - t489 * mrSges(6,2) + t515 * t506 - t516 * t507;
t472 = m(5) * t483 + t559 * mrSges(5,1) - t502 * mrSges(5,3) - t529 * t513 + t562 * t519 + t576;
t462 = t566 * t467 + t570 * t472;
t532 = -t545 * mrSges(4,1) + t546 * mrSges(4,2);
t539 = -qJD(3) * mrSges(4,2) + t545 * mrSges(4,3);
t460 = m(4) * t504 + qJDD(3) * mrSges(4,1) - t536 * mrSges(4,3) + qJD(3) * t539 - t546 * t532 + t462;
t540 = qJD(3) * mrSges(4,1) - t546 * mrSges(4,3);
t586 = t570 * t467 - t566 * t472;
t461 = m(4) * t505 - qJDD(3) * mrSges(4,2) + t535 * mrSges(4,3) - qJD(3) * t540 + t545 * t532 + t586;
t454 = t571 * t460 + t567 * t461;
t537 = -t563 * t547 + t590;
t578 = mrSges(3,3) * qJDD(1) + t573 * (-mrSges(3,1) * t564 + t597);
t452 = m(3) * t537 - t578 * t563 + t454;
t587 = -t567 * t460 + t571 * t461;
t453 = m(3) * t538 + t578 * t564 + t587;
t588 = -t563 * t452 + t564 * t453;
t446 = m(2) * t551 - t573 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t588;
t543 = -qJDD(1) * pkin(1) - t573 * qJ(2) + t584;
t468 = t569 * t476 + t565 * t477;
t577 = m(5) * t497 - t501 * mrSges(5,1) + t502 * mrSges(5,2) - t528 * t519 + t529 * t520 + t468;
t575 = m(4) * t534 - t535 * mrSges(4,1) + t536 * mrSges(4,2) - t545 * t539 + t546 * t540 + t577;
t574 = -m(3) * t543 + mrSges(3,1) * t591 - t575 + (t560 * t573 + t596) * mrSges(3,3);
t464 = t574 - t573 * mrSges(2,2) + m(2) * t550 + (mrSges(2,1) - t597) * qJDD(1);
t595 = t568 * t446 + t572 * t464;
t447 = t564 * t452 + t563 * t453;
t581 = Ifges(3,5) * t563 + Ifges(3,6) * t564;
t594 = t573 * t581;
t589 = t572 * t446 - t568 * t464;
t583 = Ifges(3,1) * t563 + Ifges(3,4) * t564;
t582 = Ifges(3,4) * t563 + Ifges(3,2) * t564;
t527 = Ifges(4,1) * t546 + Ifges(4,4) * t545 + Ifges(4,5) * qJD(3);
t526 = Ifges(4,4) * t546 + Ifges(4,2) * t545 + Ifges(4,6) * qJD(3);
t525 = Ifges(4,5) * t546 + Ifges(4,6) * t545 + Ifges(4,3) * qJD(3);
t510 = Ifges(5,1) * t529 + Ifges(5,4) * t528 + Ifges(5,5) * t562;
t509 = Ifges(5,4) * t529 + Ifges(5,2) * t528 + Ifges(5,6) * t562;
t508 = Ifges(5,5) * t529 + Ifges(5,6) * t528 + Ifges(5,3) * t562;
t494 = Ifges(6,1) * t516 + Ifges(6,4) * t515 + Ifges(6,5) * t524;
t493 = Ifges(6,4) * t516 + Ifges(6,2) * t515 + Ifges(6,6) * t524;
t492 = Ifges(6,5) * t516 + Ifges(6,6) * t515 + Ifges(6,3) * t524;
t470 = mrSges(6,2) * t480 - mrSges(6,3) * t478 + Ifges(6,1) * t489 + Ifges(6,4) * t488 + Ifges(6,5) * t500 + t515 * t492 - t524 * t493;
t469 = -mrSges(6,1) * t480 + mrSges(6,3) * t479 + Ifges(6,4) * t489 + Ifges(6,2) * t488 + Ifges(6,6) * t500 - t516 * t492 + t524 * t494;
t456 = -mrSges(5,1) * t497 - mrSges(6,1) * t478 + mrSges(6,2) * t479 + mrSges(5,3) * t484 + Ifges(5,4) * t502 - Ifges(6,5) * t489 + Ifges(5,2) * t501 + Ifges(5,6) * t559 - Ifges(6,6) * t488 - Ifges(6,3) * t500 - pkin(4) * t468 - t516 * t493 + t515 * t494 - t529 * t508 + t562 * t510;
t455 = mrSges(5,2) * t497 - mrSges(5,3) * t483 + Ifges(5,1) * t502 + Ifges(5,4) * t501 + Ifges(5,5) * t559 - pkin(8) * t468 - t565 * t469 + t569 * t470 + t528 * t508 - t562 * t509;
t448 = mrSges(4,2) * t534 - mrSges(4,3) * t504 + Ifges(4,1) * t536 + Ifges(4,4) * t535 + Ifges(4,5) * qJDD(3) - pkin(7) * t462 - qJD(3) * t526 + t570 * t455 - t566 * t456 + t545 * t525;
t443 = -mrSges(4,1) * t534 + mrSges(4,3) * t505 + Ifges(4,4) * t536 + Ifges(4,2) * t535 + Ifges(4,6) * qJDD(3) - pkin(3) * t577 + pkin(7) * t586 + qJD(3) * t527 + t566 * t455 + t570 * t456 - t546 * t525;
t442 = -Ifges(4,3) * qJDD(3) + mrSges(2,1) * g(3) - t565 * t470 - t569 * t469 - Ifges(5,3) * t559 + t545 * t527 - t546 * t526 + mrSges(2,3) * t551 - Ifges(4,6) * t535 - Ifges(4,5) * t536 - mrSges(3,1) * t537 + mrSges(3,2) * t538 + t528 * t510 - t529 * t509 - Ifges(5,6) * t501 - Ifges(5,5) * t502 - mrSges(4,1) * t504 + mrSges(4,2) * t505 + mrSges(5,2) * t484 - mrSges(5,1) * t483 - pkin(3) * t462 - pkin(2) * t454 - pkin(1) * t447 + (-t563 * t582 + t564 * t583 + Ifges(2,5)) * t573 - pkin(8) * t585 + (Ifges(2,6) - t581) * qJDD(1) - pkin(4) * t576;
t441 = mrSges(3,2) * t543 - mrSges(3,3) * t537 - pkin(6) * t454 + t583 * qJDD(1) - t567 * t443 + t571 * t448 + t564 * t594;
t440 = -mrSges(3,1) * t543 + mrSges(3,3) * t538 - pkin(2) * t575 + pkin(6) * t587 + t582 * qJDD(1) + t571 * t443 + t567 * t448 - t563 * t594;
t439 = -mrSges(2,2) * g(3) - mrSges(2,3) * t550 + Ifges(2,5) * qJDD(1) - t573 * Ifges(2,6) - qJ(2) * t447 - t563 * t440 + t564 * t441;
t1 = [-m(1) * g(1) + t589; -m(1) * g(2) + t595; (-m(1) - m(2)) * g(3) + t447; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t595 + t572 * t439 - t568 * t442; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t589 + t568 * t439 + t572 * t442; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t550 - mrSges(2,2) * t551 + t563 * t441 + t564 * t440 + pkin(1) * (-qJDD(1) * t597 + t574) + qJ(2) * t588;];
tauB = t1;
