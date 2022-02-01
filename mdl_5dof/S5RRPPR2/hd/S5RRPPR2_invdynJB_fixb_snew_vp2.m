% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:43
% EndTime: 2022-01-20 10:05:45
% DurationCPUTime: 2.92s
% Computational Cost: add. (39285->208), mult. (51725->271), div. (0->0), fcn. (28744->10), ass. (0->103)
t598 = 2 * qJD(4);
t563 = sin(qJ(1));
t566 = cos(qJ(1));
t542 = t563 * g(1) - t566 * g(2);
t536 = qJDD(1) * pkin(1) + t542;
t543 = -t566 * g(1) - t563 * g(2);
t567 = qJD(1) ^ 2;
t537 = -t567 * pkin(1) + t543;
t562 = sin(qJ(2));
t565 = cos(qJ(2));
t522 = t565 * t536 - t562 * t537;
t553 = qJDD(1) + qJDD(2);
t516 = t553 * pkin(2) + t522;
t523 = t562 * t536 + t565 * t537;
t554 = (qJD(1) + qJD(2));
t552 = t554 ^ 2;
t517 = -t552 * pkin(2) + t523;
t558 = sin(pkin(8));
t560 = cos(pkin(8));
t512 = t558 * t516 + t560 * t517;
t509 = -t552 * pkin(3) + t553 * qJ(4) + t512;
t597 = (t554 * t598) + t509;
t557 = sin(pkin(9));
t596 = mrSges(5,2) * t557;
t594 = mrSges(5,3) * t553;
t593 = t557 * t554;
t561 = sin(qJ(5));
t592 = t557 * t561;
t564 = cos(qJ(5));
t591 = t557 * t564;
t559 = cos(pkin(9));
t590 = t559 * t553;
t589 = t559 * t554;
t556 = -g(3) + qJDD(3);
t588 = t559 * t556;
t505 = t557 * t556 + t597 * t559;
t530 = (-mrSges(5,1) * t559 + t596) * t554;
t577 = -pkin(4) * t559 - pkin(7) * t557;
t532 = t577 * t554;
t503 = t532 * t589 + t505;
t511 = t560 * t516 - t558 * t517;
t572 = -t552 * qJ(4) + qJDD(4) - t511;
t506 = (-pkin(3) + t577) * t553 + t572;
t500 = -t561 * t503 + t564 * t506;
t540 = qJD(5) - t589;
t584 = t554 * t592;
t525 = -t540 * mrSges(6,2) - mrSges(6,3) * t584;
t527 = (mrSges(6,1) * t561 + mrSges(6,2) * t564) * t593;
t585 = qJD(5) * t554;
t529 = (t553 * t564 - t561 * t585) * t557;
t539 = qJDD(5) - t590;
t583 = t554 * t591;
t498 = m(6) * t500 + t539 * mrSges(6,1) - t529 * mrSges(6,3) + t540 * t525 - t527 * t583;
t501 = t564 * t503 + t561 * t506;
t526 = t540 * mrSges(6,1) - mrSges(6,3) * t583;
t528 = (-t553 * t561 - t564 * t585) * t557;
t499 = m(6) * t501 - t539 * mrSges(6,2) + t528 * mrSges(6,3) - t540 * t526 - t527 * t584;
t578 = -t561 * t498 + t564 * t499;
t489 = m(5) * t505 + (t530 * t554 + t594) * t559 + t578;
t504 = -t597 * t557 + t588;
t502 = -t588 + (t509 + (t598 + t532) * t554) * t557;
t573 = -m(6) * t502 + t528 * mrSges(6,1) - t529 * mrSges(6,2);
t496 = m(5) * t504 + (-t594 + (-t525 * t561 - t526 * t564 - t530) * t554) * t557 + t573;
t579 = t559 * t489 - t557 * t496;
t481 = m(4) * t512 - t552 * mrSges(4,1) - t553 * mrSges(4,2) + t579;
t492 = t564 * t498 + t561 * t499;
t508 = -t553 * pkin(3) + t572;
t571 = -m(5) * t508 + mrSges(5,1) * t590 - t492 + (t557 ^ 2 + t559 ^ 2) * mrSges(5,3) * t552;
t486 = m(4) * t511 - t552 * mrSges(4,2) + (mrSges(4,1) - t596) * t553 + t571;
t474 = t558 * t481 + t560 * t486;
t471 = m(3) * t522 + t553 * mrSges(3,1) - t552 * mrSges(3,2) + t474;
t580 = t560 * t481 - t558 * t486;
t472 = m(3) * t523 - t552 * mrSges(3,1) - t553 * mrSges(3,2) + t580;
t466 = t565 * t471 + t562 * t472;
t463 = m(2) * t542 + qJDD(1) * mrSges(2,1) - t567 * mrSges(2,2) + t466;
t581 = -t562 * t471 + t565 * t472;
t464 = m(2) * t543 - t567 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t581;
t587 = t566 * t463 + t563 * t464;
t484 = t557 * t489 + t559 * t496;
t482 = m(4) * t556 + t484;
t582 = -t563 * t463 + t566 * t464;
t576 = Ifges(5,1) * t557 + Ifges(5,4) * t559;
t575 = Ifges(5,5) * t557 + Ifges(5,6) * t559;
t519 = Ifges(6,6) * t540 + (Ifges(6,4) * t564 - Ifges(6,2) * t561) * t593;
t520 = Ifges(6,5) * t540 + (Ifges(6,1) * t564 - Ifges(6,4) * t561) * t593;
t574 = t519 * t564 + t520 * t561;
t570 = mrSges(6,1) * t500 - mrSges(6,2) * t501 + Ifges(6,5) * t529 + Ifges(6,6) * t528 + Ifges(6,3) * t539;
t518 = Ifges(6,3) * t540 + (Ifges(6,5) * t564 - Ifges(6,6) * t561) * t593;
t493 = -mrSges(6,1) * t502 + mrSges(6,3) * t501 + Ifges(6,4) * t529 + Ifges(6,2) * t528 + Ifges(6,6) * t539 - t518 * t583 + t540 * t520;
t494 = mrSges(6,2) * t502 - mrSges(6,3) * t500 + Ifges(6,1) * t529 + Ifges(6,4) * t528 + Ifges(6,5) * t539 - t518 * t584 - t540 * t519;
t531 = t575 * t554;
t476 = mrSges(5,2) * t508 - mrSges(5,3) * t504 - pkin(7) * t492 - t561 * t493 + t564 * t494 + t531 * t589 + t576 * t553;
t478 = Ifges(5,2) * t590 - mrSges(5,1) * t508 + mrSges(5,3) * t505 - pkin(4) * t492 + (Ifges(5,4) * t553 + (-t531 - t574) * t554) * t557 - t570;
t491 = t553 * t596 - t571;
t569 = mrSges(3,1) * t522 + mrSges(4,1) * t511 - mrSges(3,2) * t523 - mrSges(4,2) * t512 + pkin(2) * t474 - pkin(3) * t491 + qJ(4) * t579 + t557 * t476 + t559 * t478 + (Ifges(3,3) + Ifges(4,3)) * t553;
t568 = mrSges(2,1) * t542 - mrSges(2,2) * t543 + Ifges(2,3) * qJDD(1) + pkin(1) * t466 + t569;
t467 = t552 * Ifges(4,5) - mrSges(4,1) * t556 + mrSges(4,3) * t512 - mrSges(5,1) * t504 + mrSges(5,2) * t505 - t561 * t494 - t564 * t493 - pkin(4) * t573 - pkin(7) * t578 - pkin(3) * t484 + (Ifges(4,6) - t575) * t553 + (-pkin(4) * (-t525 * t592 - t526 * t591) + (-t557 * (Ifges(5,4) * t557 + Ifges(5,2) * t559) + t559 * t576) * t554) * t554;
t459 = mrSges(4,2) * t556 - mrSges(4,3) * t511 + Ifges(4,5) * t553 - t552 * Ifges(4,6) - qJ(4) * t484 + t559 * t476 - t557 * t478;
t458 = -mrSges(3,2) * g(3) - mrSges(3,3) * t522 + Ifges(3,5) * t553 - t552 * Ifges(3,6) - qJ(3) * t474 + t560 * t459 - t558 * t467;
t457 = mrSges(3,1) * g(3) + mrSges(3,3) * t523 + t552 * Ifges(3,5) + Ifges(3,6) * t553 - pkin(2) * t482 + qJ(3) * t580 + t558 * t459 + t560 * t467;
t456 = -mrSges(2,2) * g(3) - mrSges(2,3) * t542 + Ifges(2,5) * qJDD(1) - t567 * Ifges(2,6) - pkin(6) * t466 - t562 * t457 + t565 * t458;
t455 = Ifges(2,6) * qJDD(1) + t567 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t543 + t562 * t458 + t565 * t457 - pkin(1) * (-m(3) * g(3) + t482) + pkin(6) * t581;
t1 = [-m(1) * g(1) + t582; -m(1) * g(2) + t587; (-m(1) - m(2) - m(3)) * g(3) + t482; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t587 - t563 * t455 + t566 * t456; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t582 + t566 * t455 + t563 * t456; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t568; t568; t569; t482; t491; t574 * t593 + t570;];
tauJB = t1;
