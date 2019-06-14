% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:50:55
% EndTime: 2019-05-05 14:51:00
% DurationCPUTime: 2.69s
% Computational Cost: add. (26888->270), mult. (48766->315), div. (0->0), fcn. (26457->8), ass. (0->109)
t602 = Ifges(6,1) + Ifges(7,1);
t592 = Ifges(6,4) - Ifges(7,5);
t601 = -Ifges(6,5) - Ifges(7,4);
t600 = Ifges(6,2) + Ifges(7,3);
t589 = Ifges(6,6) - Ifges(7,6);
t599 = -Ifges(6,3) - Ifges(7,2);
t560 = sin(qJ(1));
t562 = cos(qJ(1));
t539 = t560 * g(1) - t562 * g(2);
t531 = qJDD(1) * pkin(1) + t539;
t540 = -t562 * g(1) - t560 * g(2);
t564 = qJD(1) ^ 2;
t533 = -t564 * pkin(1) + t540;
t556 = sin(pkin(9));
t557 = cos(pkin(9));
t506 = t556 * t531 + t557 * t533;
t598 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t506;
t597 = -pkin(2) - pkin(7);
t596 = cos(qJ(5));
t595 = mrSges(3,1) - mrSges(4,2);
t594 = -mrSges(6,3) - mrSges(7,2);
t593 = -Ifges(4,4) + Ifges(3,5);
t590 = Ifges(4,5) - Ifges(3,6);
t505 = t557 * t531 - t556 * t533;
t569 = -t564 * qJ(3) + qJDD(3) - t505;
t490 = t597 * qJDD(1) + t569;
t553 = -g(3) + qJDD(2);
t559 = sin(qJ(4));
t561 = cos(qJ(4));
t486 = t559 * t490 + t561 * t553;
t532 = (mrSges(5,1) * t559 + mrSges(5,2) * t561) * qJD(1);
t581 = qJD(1) * qJD(4);
t576 = t561 * t581;
t535 = -t559 * qJDD(1) - t576;
t583 = qJD(1) * t561;
t538 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t583;
t489 = t597 * t564 - t598;
t577 = t559 * t581;
t536 = t561 * qJDD(1) - t577;
t481 = (-t536 + t577) * pkin(8) + (-t535 + t576) * pkin(4) + t489;
t534 = (pkin(4) * t559 - pkin(8) * t561) * qJD(1);
t563 = qJD(4) ^ 2;
t582 = t559 * qJD(1);
t484 = -pkin(4) * t563 + qJDD(4) * pkin(8) - t534 * t582 + t486;
t558 = sin(qJ(5));
t479 = t558 * t481 + t596 * t484;
t530 = t558 * qJD(4) + t596 * t583;
t503 = t530 * qJD(5) - t596 * qJDD(4) + t558 * t536;
t542 = qJD(5) + t582;
t513 = t542 * mrSges(6,1) - t530 * mrSges(6,3);
t528 = qJDD(5) - t535;
t529 = -t596 * qJD(4) + t558 * t583;
t509 = t529 * pkin(5) - t530 * qJ(6);
t541 = t542 ^ 2;
t475 = -pkin(5) * t541 + qJ(6) * t528 + 0.2e1 * qJD(6) * t542 - t509 * t529 + t479;
t514 = -t542 * mrSges(7,1) + t530 * mrSges(7,2);
t579 = m(7) * t475 + t528 * mrSges(7,3) + t542 * t514;
t510 = t529 * mrSges(7,1) - t530 * mrSges(7,3);
t584 = -t529 * mrSges(6,1) - t530 * mrSges(6,2) - t510;
t470 = m(6) * t479 - t528 * mrSges(6,2) + t594 * t503 - t542 * t513 + t584 * t529 + t579;
t478 = t596 * t481 - t558 * t484;
t504 = -t529 * qJD(5) + t558 * qJDD(4) + t596 * t536;
t512 = -t542 * mrSges(6,2) - t529 * mrSges(6,3);
t476 = -t528 * pkin(5) - t541 * qJ(6) + t530 * t509 + qJDD(6) - t478;
t515 = -t529 * mrSges(7,2) + t542 * mrSges(7,3);
t571 = -m(7) * t476 + t528 * mrSges(7,1) + t542 * t515;
t472 = m(6) * t478 + t528 * mrSges(6,1) + t594 * t504 + t542 * t512 + t584 * t530 + t571;
t572 = t596 * t470 - t472 * t558;
t463 = m(5) * t486 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t535 - qJD(4) * t538 - t532 * t582 + t572;
t485 = t561 * t490 - t559 * t553;
t537 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t582;
t483 = -qJDD(4) * pkin(4) - t563 * pkin(8) + t534 * t583 - t485;
t477 = -0.2e1 * qJD(6) * t530 + (t529 * t542 - t504) * qJ(6) + (t530 * t542 + t503) * pkin(5) + t483;
t473 = m(7) * t477 + mrSges(7,1) * t503 - t504 * mrSges(7,3) - t530 * t514 + t515 * t529;
t565 = -m(6) * t483 - t503 * mrSges(6,1) - mrSges(6,2) * t504 - t529 * t512 - t513 * t530 - t473;
t467 = m(5) * t485 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t536 + qJD(4) * t537 - t532 * t583 + t565;
t458 = t559 * t463 + t561 * t467;
t492 = -qJDD(1) * pkin(2) + t569;
t568 = -m(4) * t492 + t564 * mrSges(4,3) - t458;
t456 = m(3) * t505 - t564 * mrSges(3,2) + t595 * qJDD(1) + t568;
t491 = t564 * pkin(2) + t598;
t466 = t558 * t470 + t596 * t472;
t567 = -m(5) * t489 + mrSges(5,1) * t535 - t536 * mrSges(5,2) - t537 * t582 - t538 * t583 - t466;
t566 = -m(4) * t491 + t564 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t567;
t461 = m(3) * t506 - mrSges(3,1) * t564 - qJDD(1) * mrSges(3,2) + t566;
t453 = t557 * t456 + t556 * t461;
t451 = m(2) * t539 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t564 + t453;
t574 = -t456 * t556 + t557 * t461;
t452 = m(2) * t540 - mrSges(2,1) * t564 - qJDD(1) * mrSges(2,2) + t574;
t588 = t562 * t451 + t560 * t452;
t587 = t600 * t529 - t592 * t530 - t589 * t542;
t586 = t589 * t529 + t601 * t530 + t599 * t542;
t585 = -t592 * t529 + t602 * t530 - t601 * t542;
t575 = -t451 * t560 + t562 * t452;
t573 = t561 * t463 - t559 * t467;
t457 = m(4) * t553 + t573;
t570 = m(3) * t553 + t457;
t522 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t561 - Ifges(5,4) * t559) * qJD(1);
t521 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t561 - Ifges(5,2) * t559) * qJD(1);
t520 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t561 - Ifges(5,6) * t559) * qJD(1);
t465 = mrSges(6,2) * t483 + mrSges(7,2) * t476 - mrSges(6,3) * t478 - mrSges(7,3) * t477 - qJ(6) * t473 - t592 * t503 + t602 * t504 - t528 * t601 + t586 * t529 + t587 * t542;
t464 = -mrSges(6,1) * t483 - mrSges(7,1) * t477 + mrSges(7,2) * t475 + mrSges(6,3) * t479 - pkin(5) * t473 - t600 * t503 + t592 * t504 + t589 * t528 + t586 * t530 + t585 * t542;
t454 = Ifges(5,4) * t536 + Ifges(5,2) * t535 + Ifges(5,6) * qJDD(4) - t520 * t583 + qJD(4) * t522 - mrSges(5,1) * t489 + mrSges(5,3) * t486 - mrSges(6,1) * t478 + mrSges(6,2) * t479 + mrSges(7,1) * t476 - mrSges(7,3) * t475 - pkin(5) * t571 - qJ(6) * t579 - pkin(4) * t466 + (pkin(5) * t510 + t587) * t530 + (qJ(6) * t510 - t585) * t529 + t599 * t528 + (mrSges(7,2) * pkin(5) + t601) * t504 + (mrSges(7,2) * qJ(6) + t589) * t503;
t447 = mrSges(5,2) * t489 - mrSges(5,3) * t485 + Ifges(5,1) * t536 + Ifges(5,4) * t535 + Ifges(5,5) * qJDD(4) - pkin(8) * t466 - qJD(4) * t521 - t558 * t464 + t596 * t465 - t520 * t582;
t446 = -qJ(3) * t457 - mrSges(3,3) * t505 + pkin(3) * t458 + mrSges(4,1) * t492 + t558 * t465 + t596 * t464 + pkin(4) * t565 + pkin(8) * t572 + mrSges(5,1) * t485 - mrSges(5,2) * t486 + Ifges(5,5) * t536 + Ifges(5,6) * t535 + Ifges(5,3) * qJDD(4) + t590 * t564 + (mrSges(3,2) - mrSges(4,3)) * t553 + t593 * qJDD(1) + (t521 * t561 + t522 * t559) * qJD(1);
t445 = -mrSges(4,1) * t491 + mrSges(3,3) * t506 - pkin(2) * t457 - pkin(3) * t567 - pkin(7) * t573 - t590 * qJDD(1) - t559 * t447 - t561 * t454 - t595 * t553 + t593 * t564;
t444 = -mrSges(2,2) * g(3) - mrSges(2,3) * t539 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t564 - qJ(2) * t453 - t445 * t556 + t446 * t557;
t443 = mrSges(2,1) * g(3) + mrSges(2,3) * t540 + t564 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t570 + qJ(2) * t574 + t557 * t445 + t556 * t446;
t1 = [-m(1) * g(1) + t575; -m(1) * g(2) + t588; (-m(1) - m(2)) * g(3) + t570; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t588 - t560 * t443 + t562 * t444; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t575 + t562 * t443 + t560 * t444; pkin(1) * t453 + mrSges(2,1) * t539 - mrSges(2,2) * t540 + pkin(2) * t568 + qJ(3) * t566 - pkin(7) * t458 + mrSges(3,1) * t505 - mrSges(3,2) * t506 + t561 * t447 - t559 * t454 + mrSges(4,2) * t492 - mrSges(4,3) * t491 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(4,2) * pkin(2) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
