% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:53:11
% EndTime: 2019-12-31 21:53:18
% DurationCPUTime: 4.34s
% Computational Cost: add. (42755->292), mult. (86342->357), div. (0->0), fcn. (57986->8), ass. (0->114)
t601 = Ifges(5,1) + Ifges(6,1);
t595 = Ifges(5,4) + Ifges(6,4);
t594 = Ifges(5,5) + Ifges(6,5);
t600 = Ifges(5,2) + Ifges(6,2);
t599 = Ifges(5,6) + Ifges(6,6);
t598 = Ifges(5,3) + Ifges(6,3);
t564 = sin(qJ(3));
t565 = sin(qJ(2));
t568 = cos(qJ(3));
t569 = cos(qJ(2));
t543 = (t564 * t565 - t568 * t569) * qJD(1);
t571 = qJD(1) ^ 2;
t597 = pkin(2) * t571;
t596 = -mrSges(5,2) - mrSges(6,2);
t566 = sin(qJ(1));
t570 = cos(qJ(1));
t557 = -t570 * g(1) - t566 * g(2);
t546 = -t571 * pkin(1) + qJDD(1) * pkin(6) + t557;
t592 = t565 * t546;
t584 = qJD(1) * qJD(2);
t551 = t565 * qJDD(1) + t569 * t584;
t511 = qJDD(2) * pkin(2) - t551 * pkin(7) - t592 + (pkin(7) * t584 + t565 * t597 - g(3)) * t569;
t535 = -t565 * g(3) + t569 * t546;
t552 = t569 * qJDD(1) - t565 * t584;
t586 = qJD(1) * t565;
t555 = qJD(2) * pkin(2) - pkin(7) * t586;
t562 = t569 ^ 2;
t512 = t552 * pkin(7) - qJD(2) * t555 - t562 * t597 + t535;
t490 = t564 * t511 + t568 * t512;
t544 = (t564 * t569 + t565 * t568) * qJD(1);
t517 = -t544 * qJD(3) - t564 * t551 + t568 * t552;
t529 = t543 * mrSges(4,1) + t544 * mrSges(4,2);
t561 = qJD(2) + qJD(3);
t537 = t561 * mrSges(4,1) - t544 * mrSges(4,3);
t560 = qJDD(2) + qJDD(3);
t518 = -t543 * qJD(3) + t568 * t551 + t564 * t552;
t556 = t566 * g(1) - t570 * g(2);
t575 = -qJDD(1) * pkin(1) - t556;
t519 = -t552 * pkin(2) + t555 * t586 + (-pkin(7) * t562 - pkin(6)) * t571 + t575;
t485 = (t543 * t561 - t518) * pkin(8) + (t544 * t561 - t517) * pkin(3) + t519;
t530 = t543 * pkin(3) - t544 * pkin(8);
t559 = t561 ^ 2;
t488 = -t559 * pkin(3) + t560 * pkin(8) - t543 * t530 + t490;
t563 = sin(qJ(4));
t567 = cos(qJ(4));
t480 = t567 * t485 - t563 * t488;
t532 = -t563 * t544 + t567 * t561;
t495 = t532 * qJD(4) + t567 * t518 + t563 * t560;
t533 = t567 * t544 + t563 * t561;
t509 = -t532 * mrSges(6,1) + t533 * mrSges(6,2);
t510 = -t532 * mrSges(5,1) + t533 * mrSges(5,2);
t516 = qJDD(4) - t517;
t539 = qJD(4) + t543;
t521 = -t539 * mrSges(5,2) + t532 * mrSges(5,3);
t477 = -0.2e1 * qJD(5) * t533 + (t532 * t539 - t495) * qJ(5) + (t532 * t533 + t516) * pkin(4) + t480;
t520 = -t539 * mrSges(6,2) + t532 * mrSges(6,3);
t583 = m(6) * t477 + t516 * mrSges(6,1) + t539 * t520;
t469 = m(5) * t480 + t516 * mrSges(5,1) + t539 * t521 + (-t509 - t510) * t533 + (-mrSges(5,3) - mrSges(6,3)) * t495 + t583;
t481 = t563 * t485 + t567 * t488;
t494 = -t533 * qJD(4) - t563 * t518 + t567 * t560;
t522 = t539 * pkin(4) - t533 * qJ(5);
t531 = t532 ^ 2;
t479 = -t531 * pkin(4) + t494 * qJ(5) + 0.2e1 * qJD(5) * t532 - t539 * t522 + t481;
t582 = m(6) * t479 + t494 * mrSges(6,3) + t532 * t509;
t523 = t539 * mrSges(6,1) - t533 * mrSges(6,3);
t587 = -t539 * mrSges(5,1) + t533 * mrSges(5,3) - t523;
t472 = m(5) * t481 + t494 * mrSges(5,3) + t532 * t510 + t596 * t516 + t587 * t539 + t582;
t578 = -t563 * t469 + t567 * t472;
t465 = m(4) * t490 - t560 * mrSges(4,2) + t517 * mrSges(4,3) - t543 * t529 - t561 * t537 + t578;
t489 = t568 * t511 - t564 * t512;
t536 = -t561 * mrSges(4,2) - t543 * mrSges(4,3);
t487 = -t560 * pkin(3) - t559 * pkin(8) + t544 * t530 - t489;
t482 = -t494 * pkin(4) - t531 * qJ(5) + t533 * t522 + qJDD(5) + t487;
t577 = m(6) * t482 - t494 * mrSges(6,1) - t532 * t520;
t573 = -m(5) * t487 + t494 * mrSges(5,1) + t596 * t495 + t532 * t521 + t587 * t533 - t577;
t474 = m(4) * t489 + t560 * mrSges(4,1) - t518 * mrSges(4,3) - t544 * t529 + t561 * t536 + t573;
t459 = t564 * t465 + t568 * t474;
t534 = -t569 * g(3) - t592;
t550 = (-mrSges(3,1) * t569 + mrSges(3,2) * t565) * qJD(1);
t585 = qJD(1) * t569;
t554 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t585;
t457 = m(3) * t534 + qJDD(2) * mrSges(3,1) - t551 * mrSges(3,3) + qJD(2) * t554 - t550 * t586 + t459;
t553 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t586;
t579 = t568 * t465 - t564 * t474;
t458 = m(3) * t535 - qJDD(2) * mrSges(3,2) + t552 * mrSges(3,3) - qJD(2) * t553 + t550 * t585 + t579;
t580 = -t565 * t457 + t569 * t458;
t451 = m(2) * t557 - t571 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t580;
t545 = -t571 * pkin(6) + t575;
t467 = t567 * t469 + t563 * t472;
t574 = m(4) * t519 - t517 * mrSges(4,1) + t518 * mrSges(4,2) + t543 * t536 + t544 * t537 + t467;
t572 = -m(3) * t545 + t552 * mrSges(3,1) - t551 * mrSges(3,2) - t553 * t586 + t554 * t585 - t574;
t462 = m(2) * t556 + qJDD(1) * mrSges(2,1) - t571 * mrSges(2,2) + t572;
t591 = t566 * t451 + t570 * t462;
t452 = t569 * t457 + t565 * t458;
t590 = t599 * t532 + t594 * t533 + t598 * t539;
t589 = -t600 * t532 - t595 * t533 - t599 * t539;
t588 = t595 * t532 + t601 * t533 + t594 * t539;
t581 = t570 * t451 - t566 * t462;
t542 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t565 + Ifges(3,4) * t569) * qJD(1);
t541 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t565 + Ifges(3,2) * t569) * qJD(1);
t540 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t565 + Ifges(3,6) * t569) * qJD(1);
t527 = Ifges(4,1) * t544 - Ifges(4,4) * t543 + Ifges(4,5) * t561;
t526 = Ifges(4,4) * t544 - Ifges(4,2) * t543 + Ifges(4,6) * t561;
t525 = Ifges(4,5) * t544 - Ifges(4,6) * t543 + Ifges(4,3) * t561;
t475 = -t495 * mrSges(6,3) - t533 * t509 + t583;
t466 = mrSges(5,2) * t487 + mrSges(6,2) * t482 - mrSges(5,3) * t480 - mrSges(6,3) * t477 - qJ(5) * t475 + t595 * t494 + t601 * t495 + t594 * t516 + t590 * t532 + t589 * t539;
t460 = -mrSges(5,1) * t487 + mrSges(5,3) * t481 - mrSges(6,1) * t482 + mrSges(6,3) * t479 - pkin(4) * t577 + qJ(5) * t582 + (-qJ(5) * t523 + t588) * t539 + (-pkin(4) * t523 - t590) * t533 + (-qJ(5) * mrSges(6,2) + t599) * t516 + (-pkin(4) * mrSges(6,2) + t595) * t495 + t600 * t494;
t453 = -mrSges(4,1) * t519 - mrSges(5,1) * t480 - mrSges(6,1) * t477 + mrSges(5,2) * t481 + mrSges(6,2) * t479 + mrSges(4,3) * t490 + Ifges(4,4) * t518 + Ifges(4,2) * t517 + Ifges(4,6) * t560 - pkin(3) * t467 - pkin(4) * t475 - t544 * t525 + t561 * t527 + t589 * t533 + t588 * t532 - t598 * t516 - t594 * t495 - t599 * t494;
t448 = mrSges(4,2) * t519 - mrSges(4,3) * t489 + Ifges(4,1) * t518 + Ifges(4,4) * t517 + Ifges(4,5) * t560 - pkin(8) * t467 - t563 * t460 + t567 * t466 - t543 * t525 - t561 * t526;
t447 = mrSges(3,2) * t545 - mrSges(3,3) * t534 + Ifges(3,1) * t551 + Ifges(3,4) * t552 + Ifges(3,5) * qJDD(2) - pkin(7) * t459 - qJD(2) * t541 + t568 * t448 - t564 * t453 + t540 * t585;
t446 = -mrSges(3,1) * t545 + mrSges(3,3) * t535 + Ifges(3,4) * t551 + Ifges(3,2) * t552 + Ifges(3,6) * qJDD(2) - pkin(2) * t574 + pkin(7) * t579 + qJD(2) * t542 + t564 * t448 + t568 * t453 - t540 * t586;
t445 = -pkin(1) * t452 + mrSges(2,3) * t557 - pkin(2) * t459 - Ifges(3,5) * t551 - Ifges(3,6) * t552 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t534 + mrSges(3,2) * t535 - t567 * t460 - pkin(3) * t573 - pkin(8) * t578 - Ifges(4,5) * t518 - Ifges(4,6) * t517 - Ifges(4,3) * t560 - mrSges(4,1) * t489 + mrSges(4,2) * t490 - t563 * t466 + mrSges(2,1) * g(3) + t571 * Ifges(2,5) - t544 * t526 - t543 * t527 + Ifges(2,6) * qJDD(1) + (-t565 * t541 + t569 * t542) * qJD(1);
t444 = -mrSges(2,2) * g(3) - mrSges(2,3) * t556 + Ifges(2,5) * qJDD(1) - t571 * Ifges(2,6) - pkin(6) * t452 - t565 * t446 + t569 * t447;
t1 = [-m(1) * g(1) + t581; -m(1) * g(2) + t591; (-m(1) - m(2)) * g(3) + t452; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t591 + t570 * t444 - t566 * t445; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t581 + t566 * t444 + t570 * t445; -mrSges(1,1) * g(2) + mrSges(2,1) * t556 + mrSges(1,2) * g(1) - mrSges(2,2) * t557 + Ifges(2,3) * qJDD(1) + pkin(1) * t572 + pkin(6) * t580 + t569 * t446 + t565 * t447;];
tauB = t1;
