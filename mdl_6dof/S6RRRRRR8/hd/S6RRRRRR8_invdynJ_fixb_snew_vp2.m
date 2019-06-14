% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRR8
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 14:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 13:39:50
% EndTime: 2019-05-08 13:40:47
% DurationCPUTime: 40.26s
% Computational Cost: add. (643391->380), mult. (1577062->506), div. (0->0), fcn. (1346246->16), ass. (0->170)
t580 = cos(pkin(6));
t575 = qJD(1) * t580 + qJD(2);
t577 = sin(pkin(7));
t579 = cos(pkin(7));
t578 = sin(pkin(6));
t591 = cos(qJ(2));
t613 = qJD(1) * t591;
t610 = t578 * t613;
t601 = t575 * t577 + t579 * t610;
t560 = t601 * pkin(10);
t585 = sin(qJ(2));
t615 = qJD(1) * t578;
t626 = pkin(10) * t577;
t565 = (-pkin(2) * t591 - t585 * t626) * t615;
t612 = qJD(1) * qJD(2);
t571 = (qJDD(1) * t585 + t591 * t612) * t578;
t574 = qJDD(1) * t580 + qJDD(2);
t593 = qJD(1) ^ 2;
t586 = sin(qJ(1));
t592 = cos(qJ(1));
t604 = -g(1) * t592 - g(2) * t586;
t627 = pkin(9) * t578;
t569 = -pkin(1) * t593 + qJDD(1) * t627 + t604;
t609 = t586 * g(1) - g(2) * t592;
t568 = qJDD(1) * pkin(1) + t593 * t627 + t609;
t622 = t568 * t580;
t605 = -t585 * t569 + t591 * t622;
t614 = qJD(1) * t585;
t625 = pkin(10) * t579;
t517 = -t571 * t625 + t574 * pkin(2) + t575 * t560 + (-g(3) * t591 - t565 * t614) * t578 + t605;
t611 = t578 * t614;
t563 = pkin(2) * t575 - t611 * t625;
t572 = (qJDD(1) * t591 - t585 * t612) * t578;
t602 = t572 * t579 + t574 * t577;
t616 = t591 * t569 + t585 * t622;
t518 = -t575 * t563 + (-g(3) * t585 + t565 * t613) * t578 + t602 * pkin(10) + t616;
t624 = t580 * g(3);
t524 = -t571 * t626 - t572 * pkin(2) - t624 + (-t568 + (-t560 * t591 + t563 * t585) * qJD(1)) * t578;
t584 = sin(qJ(3));
t590 = cos(qJ(3));
t489 = -t584 * t518 + (t517 * t579 + t524 * t577) * t590;
t550 = -t584 * t611 + t601 * t590;
t618 = t579 * t584;
t621 = t577 * t584;
t551 = t575 * t621 + (t585 * t590 + t591 * t618) * t615;
t535 = -t551 * qJD(3) - t584 * t571 + t602 * t590;
t536 = t550 * qJD(3) + t590 * t571 + t602 * t584;
t537 = -mrSges(4,1) * t550 + mrSges(4,2) * t551;
t561 = t575 * t579 - t577 * t610 + qJD(3);
t542 = -mrSges(4,2) * t561 + mrSges(4,3) * t550;
t552 = -t572 * t577 + t574 * t579 + qJDD(3);
t538 = -pkin(3) * t550 - pkin(11) * t551;
t559 = t561 ^ 2;
t473 = -t552 * pkin(3) - t559 * pkin(11) + t551 * t538 - t489;
t583 = sin(qJ(4));
t589 = cos(qJ(4));
t541 = t551 * t589 + t561 * t583;
t502 = -qJD(4) * t541 - t536 * t583 + t552 * t589;
t540 = -t551 * t583 + t561 * t589;
t503 = qJD(4) * t540 + t536 * t589 + t552 * t583;
t549 = qJD(4) - t550;
t526 = -mrSges(5,2) * t549 + mrSges(5,3) * t540;
t527 = mrSges(5,1) * t549 - mrSges(5,3) * t541;
t490 = t517 * t618 + t590 * t518 + t524 * t621;
t474 = -pkin(3) * t559 + pkin(11) * t552 + t538 * t550 + t490;
t500 = -t577 * t517 + t579 * t524;
t477 = (-t550 * t561 - t536) * pkin(11) + (t551 * t561 - t535) * pkin(3) + t500;
t466 = -t583 * t474 + t589 * t477;
t534 = qJDD(4) - t535;
t463 = (t540 * t549 - t503) * pkin(12) + (t540 * t541 + t534) * pkin(4) + t466;
t467 = t589 * t474 + t583 * t477;
t528 = pkin(4) * t549 - pkin(12) * t541;
t539 = t540 ^ 2;
t465 = -pkin(4) * t539 + pkin(12) * t502 - t528 * t549 + t467;
t582 = sin(qJ(5));
t588 = cos(qJ(5));
t460 = t582 * t463 + t588 * t465;
t519 = t540 * t588 - t541 * t582;
t520 = t540 * t582 + t541 * t588;
t499 = -pkin(5) * t519 - pkin(13) * t520;
t533 = qJDD(5) + t534;
t547 = qJD(5) + t549;
t546 = t547 ^ 2;
t457 = -pkin(5) * t546 + pkin(13) * t533 + t499 * t519 + t460;
t468 = -t502 * pkin(4) - t539 * pkin(12) + t541 * t528 + t473;
t482 = -qJD(5) * t520 + t502 * t588 - t503 * t582;
t483 = qJD(5) * t519 + t502 * t582 + t503 * t588;
t461 = (-t519 * t547 - t483) * pkin(13) + (t520 * t547 - t482) * pkin(5) + t468;
t581 = sin(qJ(6));
t587 = cos(qJ(6));
t454 = -t457 * t581 + t461 * t587;
t504 = -t520 * t581 + t547 * t587;
t471 = qJD(6) * t504 + t483 * t587 + t533 * t581;
t481 = qJDD(6) - t482;
t505 = t520 * t587 + t547 * t581;
t491 = -mrSges(7,1) * t504 + mrSges(7,2) * t505;
t516 = qJD(6) - t519;
t492 = -mrSges(7,2) * t516 + mrSges(7,3) * t504;
t450 = m(7) * t454 + mrSges(7,1) * t481 - mrSges(7,3) * t471 - t491 * t505 + t492 * t516;
t455 = t457 * t587 + t461 * t581;
t470 = -qJD(6) * t505 - t483 * t581 + t533 * t587;
t493 = mrSges(7,1) * t516 - mrSges(7,3) * t505;
t451 = m(7) * t455 - mrSges(7,2) * t481 + mrSges(7,3) * t470 + t491 * t504 - t493 * t516;
t439 = t587 * t450 + t581 * t451;
t506 = -mrSges(6,2) * t547 + mrSges(6,3) * t519;
t507 = mrSges(6,1) * t547 - mrSges(6,3) * t520;
t598 = m(6) * t468 - t482 * mrSges(6,1) + mrSges(6,2) * t483 - t519 * t506 + t507 * t520 + t439;
t595 = -m(5) * t473 + t502 * mrSges(5,1) - mrSges(5,2) * t503 + t540 * t526 - t527 * t541 - t598;
t434 = m(4) * t489 + mrSges(4,1) * t552 - mrSges(4,3) * t536 - t537 * t551 + t542 * t561 + t595;
t623 = t434 * t590;
t620 = t578 * t585;
t619 = t578 * t591;
t498 = -mrSges(6,1) * t519 + mrSges(6,2) * t520;
t606 = -t450 * t581 + t587 * t451;
t437 = m(6) * t460 - mrSges(6,2) * t533 + mrSges(6,3) * t482 + t498 * t519 - t507 * t547 + t606;
t459 = t463 * t588 - t465 * t582;
t456 = -pkin(5) * t533 - pkin(13) * t546 + t499 * t520 - t459;
t599 = -m(7) * t456 + t470 * mrSges(7,1) - mrSges(7,2) * t471 + t504 * t492 - t493 * t505;
t446 = m(6) * t459 + mrSges(6,1) * t533 - mrSges(6,3) * t483 - t498 * t520 + t506 * t547 + t599;
t431 = t582 * t437 + t588 * t446;
t521 = -mrSges(5,1) * t540 + mrSges(5,2) * t541;
t429 = m(5) * t466 + mrSges(5,1) * t534 - mrSges(5,3) * t503 - t521 * t541 + t526 * t549 + t431;
t607 = t588 * t437 - t446 * t582;
t430 = m(5) * t467 - mrSges(5,2) * t534 + mrSges(5,3) * t502 + t521 * t540 - t527 * t549 + t607;
t423 = t589 * t429 + t583 * t430;
t543 = mrSges(4,1) * t561 - mrSges(4,3) * t551;
t608 = -t429 * t583 + t589 * t430;
t420 = m(4) * t490 - mrSges(4,2) * t552 + mrSges(4,3) * t535 + t537 * t550 - t543 * t561 + t608;
t422 = m(4) * t500 - mrSges(4,1) * t535 + mrSges(4,2) * t536 - t542 * t550 + t543 * t551 + t423;
t411 = t420 * t621 + t579 * t422 + t577 * t623;
t416 = t590 * t420 - t434 * t584;
t412 = t420 * t618 - t422 * t577 + t579 * t623;
t484 = Ifges(7,5) * t505 + Ifges(7,6) * t504 + Ifges(7,3) * t516;
t486 = Ifges(7,1) * t505 + Ifges(7,4) * t504 + Ifges(7,5) * t516;
t443 = -mrSges(7,1) * t456 + mrSges(7,3) * t455 + Ifges(7,4) * t471 + Ifges(7,2) * t470 + Ifges(7,6) * t481 - t484 * t505 + t486 * t516;
t485 = Ifges(7,4) * t505 + Ifges(7,2) * t504 + Ifges(7,6) * t516;
t444 = mrSges(7,2) * t456 - mrSges(7,3) * t454 + Ifges(7,1) * t471 + Ifges(7,4) * t470 + Ifges(7,5) * t481 + t484 * t504 - t485 * t516;
t494 = Ifges(6,5) * t520 + Ifges(6,6) * t519 + Ifges(6,3) * t547;
t495 = Ifges(6,4) * t520 + Ifges(6,2) * t519 + Ifges(6,6) * t547;
t424 = mrSges(6,2) * t468 - mrSges(6,3) * t459 + Ifges(6,1) * t483 + Ifges(6,4) * t482 + Ifges(6,5) * t533 - pkin(13) * t439 - t443 * t581 + t444 * t587 + t494 * t519 - t495 * t547;
t496 = Ifges(6,1) * t520 + Ifges(6,4) * t519 + Ifges(6,5) * t547;
t596 = mrSges(7,1) * t454 - mrSges(7,2) * t455 + Ifges(7,5) * t471 + Ifges(7,6) * t470 + Ifges(7,3) * t481 + t485 * t505 - t486 * t504;
t425 = -mrSges(6,1) * t468 + mrSges(6,3) * t460 + Ifges(6,4) * t483 + Ifges(6,2) * t482 + Ifges(6,6) * t533 - pkin(5) * t439 - t494 * t520 + t496 * t547 - t596;
t508 = Ifges(5,5) * t541 + Ifges(5,6) * t540 + Ifges(5,3) * t549;
t510 = Ifges(5,1) * t541 + Ifges(5,4) * t540 + Ifges(5,5) * t549;
t413 = -mrSges(5,1) * t473 + mrSges(5,3) * t467 + Ifges(5,4) * t503 + Ifges(5,2) * t502 + Ifges(5,6) * t534 - pkin(4) * t598 + pkin(12) * t607 + t582 * t424 + t588 * t425 - t541 * t508 + t549 * t510;
t509 = Ifges(5,4) * t541 + Ifges(5,2) * t540 + Ifges(5,6) * t549;
t414 = mrSges(5,2) * t473 - mrSges(5,3) * t466 + Ifges(5,1) * t503 + Ifges(5,4) * t502 + Ifges(5,5) * t534 - pkin(12) * t431 + t424 * t588 - t425 * t582 + t508 * t540 - t509 * t549;
t529 = Ifges(4,5) * t551 + Ifges(4,6) * t550 + Ifges(4,3) * t561;
t530 = Ifges(4,4) * t551 + Ifges(4,2) * t550 + Ifges(4,6) * t561;
t408 = mrSges(4,2) * t500 - mrSges(4,3) * t489 + Ifges(4,1) * t536 + Ifges(4,4) * t535 + Ifges(4,5) * t552 - pkin(11) * t423 - t413 * t583 + t414 * t589 + t529 * t550 - t530 * t561;
t531 = Ifges(4,1) * t551 + Ifges(4,4) * t550 + Ifges(4,5) * t561;
t597 = -mrSges(6,1) * t459 + mrSges(6,2) * t460 - Ifges(6,5) * t483 - Ifges(6,6) * t482 - Ifges(6,3) * t533 - pkin(5) * t599 - pkin(13) * t606 - t587 * t443 - t581 * t444 - t520 * t495 + t519 * t496;
t594 = mrSges(5,1) * t466 - mrSges(5,2) * t467 + Ifges(5,5) * t503 + Ifges(5,6) * t502 + Ifges(5,3) * t534 + pkin(4) * t431 + t541 * t509 - t540 * t510 - t597;
t409 = -mrSges(4,1) * t500 + mrSges(4,3) * t490 + Ifges(4,4) * t536 + Ifges(4,2) * t535 + Ifges(4,6) * t552 - pkin(3) * t423 - t551 * t529 + t561 * t531 - t594;
t600 = pkin(10) * t416 + t408 * t584 + t409 * t590;
t570 = (-mrSges(3,1) * t591 + mrSges(3,2) * t585) * t615;
t567 = -mrSges(3,2) * t575 + mrSges(3,3) * t610;
t566 = mrSges(3,1) * t575 - mrSges(3,3) * t611;
t556 = -t578 * t568 - t624;
t555 = Ifges(3,5) * t575 + (Ifges(3,1) * t585 + Ifges(3,4) * t591) * t615;
t554 = Ifges(3,6) * t575 + (Ifges(3,4) * t585 + Ifges(3,2) * t591) * t615;
t553 = Ifges(3,3) * t575 + (Ifges(3,5) * t585 + Ifges(3,6) * t591) * t615;
t545 = -g(3) * t620 + t616;
t544 = -g(3) * t619 + t605;
t415 = m(3) * t545 - mrSges(3,2) * t574 + mrSges(3,3) * t572 - t566 * t575 + t570 * t610 + t416;
t410 = m(3) * t544 + mrSges(3,1) * t574 - mrSges(3,3) * t571 + t567 * t575 - t570 * t611 + t412;
t407 = mrSges(4,1) * t489 - mrSges(4,2) * t490 + Ifges(4,5) * t536 + Ifges(4,6) * t535 + Ifges(4,3) * t552 + pkin(3) * t595 + pkin(11) * t608 + t589 * t413 + t583 * t414 + t551 * t530 - t550 * t531;
t406 = mrSges(3,1) * t544 - mrSges(3,2) * t545 + Ifges(3,5) * t571 + Ifges(3,6) * t572 + Ifges(3,3) * t574 + pkin(2) * t412 + t579 * t407 + (t554 * t585 - t555 * t591) * t615 + t600 * t577;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t609 - mrSges(2,2) * t604 + (t553 * t610 + mrSges(3,2) * t556 - mrSges(3,3) * t544 + Ifges(3,1) * t571 + Ifges(3,4) * t572 + Ifges(3,5) * t574 + t590 * t408 - t584 * t409 - t575 * t554 + (-t411 * t577 - t412 * t579) * pkin(10)) * t620 + (-mrSges(3,1) * t556 + mrSges(3,3) * t545 + Ifges(3,4) * t571 + Ifges(3,2) * t572 + Ifges(3,6) * t574 - pkin(2) * t411 - t577 * t407 - t553 * t611 + t575 * t555 + t579 * t600) * t619 + t580 * t406 + pkin(1) * ((t410 * t591 + t415 * t585) * t580 + (-m(3) * t556 + t572 * mrSges(3,1) - t571 * mrSges(3,2) + (-t566 * t585 + t567 * t591) * t615 - t411) * t578) + (-t410 * t585 + t415 * t591) * t627; t406; t407; t594; -t597; t596;];
tauJ  = t1;
