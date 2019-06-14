% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 23:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:22:27
% EndTime: 2019-05-05 23:22:52
% DurationCPUTime: 25.20s
% Computational Cost: add. (368846->356), mult. (1157619->492), div. (0->0), fcn. (996501->16), ass. (0->167)
t571 = sin(pkin(12));
t573 = sin(pkin(6));
t575 = cos(pkin(12));
t577 = cos(pkin(6));
t580 = sin(qJ(3));
t576 = cos(pkin(7));
t584 = cos(qJ(3));
t614 = t576 * t584;
t572 = sin(pkin(7));
t618 = t572 * t584;
t590 = t573 * (-t571 * t580 + t575 * t614) + t577 * t618;
t547 = t590 * qJD(1);
t615 = t576 * t580;
t619 = t572 * t580;
t592 = t577 * t619 + (t571 * t584 + t575 * t615) * t573;
t548 = t592 * qJD(1);
t536 = -t548 * qJD(3) + qJDD(1) * t590;
t616 = t573 * t576;
t559 = (t572 * t577 + t575 * t616) * qJD(1) * pkin(9);
t586 = qJD(1) ^ 2;
t581 = sin(qJ(1));
t585 = cos(qJ(1));
t604 = -g(1) * t585 - g(2) * t581;
t622 = qJ(2) * t573;
t563 = -pkin(1) * t586 + qJDD(1) * t622 + t604;
t625 = pkin(9) * t572;
t599 = -pkin(2) * t575 - t571 * t625;
t613 = qJD(1) * t573;
t623 = pkin(9) * qJDD(1);
t594 = qJD(1) * t599 * t613 + t576 * t623;
t609 = qJD(2) * t613;
t617 = t573 * t575;
t608 = t581 * g(1) - g(2) * t585;
t562 = qJDD(1) * pkin(1) + t586 * t622 + t608;
t621 = t562 * t577;
t600 = -g(3) * t617 - 0.2e1 * t571 * t609 + t575 * t621;
t516 = (pkin(2) * qJDD(1) + qJD(1) * t559) * t577 + (-t573 * t594 - t563) * t571 + t600;
t626 = pkin(9) * t571;
t564 = (pkin(2) * t577 - t616 * t626) * qJD(1);
t628 = 0.2e1 * t575;
t610 = t575 * t563 + t571 * t621 + t609 * t628;
t517 = (-qJD(1) * t564 + t572 * t623) * t577 + (-g(3) * t571 + t575 * t594) * t573 + t610;
t607 = -g(3) * t577 + qJDD(2);
t526 = (-t562 + t599 * qJDD(1) + (-t559 * t575 + t564 * t571) * qJD(1)) * t573 + t607;
t485 = -t580 * t517 + (t516 * t576 + t526 * t572) * t584;
t486 = t516 * t615 + t584 * t517 + t526 * t619;
t535 = -pkin(3) * t547 - pkin(10) * t548;
t595 = -t572 * t617 + t576 * t577;
t560 = qJD(1) * t595 + qJD(3);
t556 = t560 ^ 2;
t557 = qJDD(1) * t595 + qJDD(3);
t477 = -pkin(3) * t556 + pkin(10) * t557 + t535 * t547 + t486;
t500 = -t572 * t516 + t576 * t526;
t537 = t547 * qJD(3) + qJDD(1) * t592;
t480 = (-t547 * t560 - t537) * pkin(10) + (t548 * t560 - t536) * pkin(3) + t500;
t579 = sin(qJ(4));
t583 = cos(qJ(4));
t469 = -t579 * t477 + t583 * t480;
t541 = -t548 * t579 + t560 * t583;
t512 = qJD(4) * t541 + t537 * t583 + t557 * t579;
t533 = qJDD(4) - t536;
t542 = t548 * t583 + t560 * t579;
t546 = qJD(4) - t547;
t466 = (t541 * t546 - t512) * qJ(5) + (t541 * t542 + t533) * pkin(4) + t469;
t470 = t583 * t477 + t579 * t480;
t511 = -qJD(4) * t542 - t537 * t579 + t557 * t583;
t528 = pkin(4) * t546 - qJ(5) * t542;
t540 = t541 ^ 2;
t468 = -pkin(4) * t540 + qJ(5) * t511 - t528 * t546 + t470;
t570 = sin(pkin(13));
t574 = cos(pkin(13));
t520 = t541 * t574 - t542 * t570;
t627 = 2 * qJD(5);
t463 = t570 * t466 + t574 * t468 + t520 * t627;
t521 = t541 * t570 + t542 * t574;
t499 = -pkin(5) * t520 - pkin(11) * t521;
t545 = t546 ^ 2;
t461 = -pkin(5) * t545 + pkin(11) * t533 + t499 * t520 + t463;
t476 = -t557 * pkin(3) - t556 * pkin(10) + t548 * t535 - t485;
t471 = -pkin(4) * t511 - qJ(5) * t540 + t542 * t528 + qJDD(5) + t476;
t490 = t511 * t574 - t512 * t570;
t491 = t511 * t570 + t512 * t574;
t464 = (-t520 * t546 - t491) * pkin(11) + (t521 * t546 - t490) * pkin(5) + t471;
t578 = sin(qJ(6));
t582 = cos(qJ(6));
t458 = -t461 * t578 + t464 * t582;
t501 = -t521 * t578 + t546 * t582;
t474 = qJD(6) * t501 + t491 * t582 + t533 * t578;
t502 = t521 * t582 + t546 * t578;
t487 = -mrSges(7,1) * t501 + mrSges(7,2) * t502;
t489 = qJDD(6) - t490;
t519 = qJD(6) - t520;
t492 = -mrSges(7,2) * t519 + mrSges(7,3) * t501;
t455 = m(7) * t458 + mrSges(7,1) * t489 - mrSges(7,3) * t474 - t487 * t502 + t492 * t519;
t459 = t461 * t582 + t464 * t578;
t473 = -qJD(6) * t502 - t491 * t578 + t533 * t582;
t493 = mrSges(7,1) * t519 - mrSges(7,3) * t502;
t456 = m(7) * t459 - mrSges(7,2) * t489 + mrSges(7,3) * t473 + t487 * t501 - t493 * t519;
t447 = -t455 * t578 + t582 * t456;
t498 = -mrSges(6,1) * t520 + mrSges(6,2) * t521;
t504 = mrSges(6,1) * t546 - mrSges(6,3) * t521;
t443 = m(6) * t463 - mrSges(6,2) * t533 + mrSges(6,3) * t490 + t498 * t520 - t504 * t546 + t447;
t602 = -t574 * t466 + t570 * t468;
t460 = -t533 * pkin(5) - t545 * pkin(11) + (t627 + t499) * t521 + t602;
t457 = -m(7) * t460 + t473 * mrSges(7,1) - mrSges(7,2) * t474 + t501 * t492 - t493 * t502;
t462 = -0.2e1 * qJD(5) * t521 - t602;
t503 = -mrSges(6,2) * t546 + mrSges(6,3) * t520;
t451 = m(6) * t462 + mrSges(6,1) * t533 - mrSges(6,3) * t491 - t498 * t521 + t503 * t546 + t457;
t438 = t570 * t443 + t574 * t451;
t481 = Ifges(7,5) * t502 + Ifges(7,6) * t501 + Ifges(7,3) * t519;
t483 = Ifges(7,1) * t502 + Ifges(7,4) * t501 + Ifges(7,5) * t519;
t448 = -mrSges(7,1) * t460 + mrSges(7,3) * t459 + Ifges(7,4) * t474 + Ifges(7,2) * t473 + Ifges(7,6) * t489 - t481 * t502 + t483 * t519;
t482 = Ifges(7,4) * t502 + Ifges(7,2) * t501 + Ifges(7,6) * t519;
t449 = mrSges(7,2) * t460 - mrSges(7,3) * t458 + Ifges(7,1) * t474 + Ifges(7,4) * t473 + Ifges(7,5) * t489 + t481 * t501 - t482 * t519;
t495 = Ifges(6,4) * t521 + Ifges(6,2) * t520 + Ifges(6,6) * t546;
t496 = Ifges(6,1) * t521 + Ifges(6,4) * t520 + Ifges(6,5) * t546;
t506 = Ifges(5,4) * t542 + Ifges(5,2) * t541 + Ifges(5,6) * t546;
t507 = Ifges(5,1) * t542 + Ifges(5,4) * t541 + Ifges(5,5) * t546;
t629 = Ifges(5,5) * t512 + Ifges(5,6) * t511 + t542 * t506 - t541 * t507 + mrSges(5,1) * t469 - mrSges(5,2) * t470 + Ifges(6,5) * t491 + Ifges(6,6) * t490 + t521 * t495 - t520 * t496 + mrSges(6,1) * t462 - mrSges(6,2) * t463 + t578 * t449 + t582 * t448 + pkin(5) * t457 + pkin(11) * t447 + pkin(4) * t438 + (Ifges(5,3) + Ifges(6,3)) * t533;
t620 = t571 * t573;
t522 = -mrSges(5,1) * t541 + mrSges(5,2) * t542;
t527 = -mrSges(5,2) * t546 + mrSges(5,3) * t541;
t436 = m(5) * t469 + mrSges(5,1) * t533 - mrSges(5,3) * t512 - t522 * t542 + t527 * t546 + t438;
t529 = mrSges(5,1) * t546 - mrSges(5,3) * t542;
t605 = t574 * t443 - t451 * t570;
t437 = m(5) * t470 - mrSges(5,2) * t533 + mrSges(5,3) * t511 + t522 * t541 - t529 * t546 + t605;
t430 = t583 * t436 + t579 * t437;
t446 = t582 * t455 + t578 * t456;
t534 = -mrSges(4,1) * t547 + mrSges(4,2) * t548;
t544 = mrSges(4,1) * t560 - mrSges(4,3) * t548;
t606 = -t436 * t579 + t583 * t437;
t427 = m(4) * t486 - mrSges(4,2) * t557 + mrSges(4,3) * t536 + t534 * t547 - t544 * t560 + t606;
t543 = -mrSges(4,2) * t560 + mrSges(4,3) * t547;
t429 = m(4) * t500 - mrSges(4,1) * t536 + mrSges(4,2) * t537 - t543 * t547 + t544 * t548 + t430;
t445 = m(6) * t471 - t490 * mrSges(6,1) + mrSges(6,2) * t491 - t520 * t503 + t504 * t521 + t446;
t588 = -m(5) * t476 + t511 * mrSges(5,1) - mrSges(5,2) * t512 + t541 * t527 - t529 * t542 - t445;
t444 = m(4) * t485 + mrSges(4,1) * t557 - mrSges(4,3) * t537 - t534 * t548 + t543 * t560 + t588;
t418 = t427 * t619 + t576 * t429 + t444 * t618;
t423 = t584 * t427 - t444 * t580;
t419 = t427 * t615 - t429 * t572 + t444 * t614;
t603 = -mrSges(3,1) * t575 + mrSges(3,2) * t571;
t598 = mrSges(3,1) * t577 - mrSges(3,3) * t620;
t597 = -mrSges(3,2) * t577 + mrSges(3,3) * t617;
t494 = Ifges(6,5) * t521 + Ifges(6,6) * t520 + Ifges(6,3) * t546;
t431 = mrSges(6,2) * t471 - mrSges(6,3) * t462 + Ifges(6,1) * t491 + Ifges(6,4) * t490 + Ifges(6,5) * t533 - pkin(11) * t446 - t448 * t578 + t449 * t582 + t494 * t520 - t495 * t546;
t589 = mrSges(7,1) * t458 - mrSges(7,2) * t459 + Ifges(7,5) * t474 + Ifges(7,6) * t473 + Ifges(7,3) * t489 + t482 * t502 - t483 * t501;
t432 = -mrSges(6,1) * t471 + mrSges(6,3) * t463 + Ifges(6,4) * t491 + Ifges(6,2) * t490 + Ifges(6,6) * t533 - pkin(5) * t446 - t494 * t521 + t496 * t546 - t589;
t505 = Ifges(5,5) * t542 + Ifges(5,6) * t541 + Ifges(5,3) * t546;
t420 = -mrSges(5,1) * t476 + mrSges(5,3) * t470 + Ifges(5,4) * t512 + Ifges(5,2) * t511 + Ifges(5,6) * t533 - pkin(4) * t445 + qJ(5) * t605 + t570 * t431 + t574 * t432 - t542 * t505 + t546 * t507;
t421 = mrSges(5,2) * t476 - mrSges(5,3) * t469 + Ifges(5,1) * t512 + Ifges(5,4) * t511 + Ifges(5,5) * t533 - qJ(5) * t438 + t431 * t574 - t432 * t570 + t505 * t541 - t506 * t546;
t530 = Ifges(4,5) * t548 + Ifges(4,6) * t547 + Ifges(4,3) * t560;
t531 = Ifges(4,4) * t548 + Ifges(4,2) * t547 + Ifges(4,6) * t560;
t414 = mrSges(4,2) * t500 - mrSges(4,3) * t485 + Ifges(4,1) * t537 + Ifges(4,4) * t536 + Ifges(4,5) * t557 - pkin(10) * t430 - t420 * t579 + t421 * t583 + t530 * t547 - t531 * t560;
t532 = Ifges(4,1) * t548 + Ifges(4,4) * t547 + Ifges(4,5) * t560;
t415 = -mrSges(4,1) * t500 + mrSges(4,3) * t486 + Ifges(4,4) * t537 + Ifges(4,2) * t536 + Ifges(4,6) * t557 - pkin(3) * t430 - t548 * t530 + t560 * t532 - t629;
t593 = pkin(9) * t423 + t414 * t580 + t415 * t584;
t566 = t597 * qJD(1);
t565 = t598 * qJD(1);
t561 = t603 * t613;
t549 = -t562 * t573 + t607;
t539 = -g(3) * t620 + t610;
t538 = -t563 * t571 + t600;
t422 = m(3) * t539 + t597 * qJDD(1) + (t561 * t617 - t565 * t577) * qJD(1) + t423;
t417 = m(3) * t549 + (t603 * qJDD(1) + (t565 * t571 - t566 * t575) * qJD(1)) * t573 + t418;
t416 = m(3) * t538 + t598 * qJDD(1) + (-t561 * t620 + t566 * t577) * qJD(1) + t419;
t413 = mrSges(4,1) * t485 - mrSges(4,2) * t486 + Ifges(4,5) * t537 + Ifges(4,6) * t536 + Ifges(4,3) * t557 + pkin(3) * t588 + pkin(10) * t606 + t583 * t420 + t579 * t421 + t548 * t531 - t547 * t532;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t608 - mrSges(2,2) * t604 + (mrSges(3,1) * t538 - mrSges(3,2) * t539 + pkin(2) * t419 + t576 * t413 + pkin(1) * (t416 * t575 + t422 * t571) + Ifges(3,3) * t577 * qJDD(1) + t593 * t572) * t577 + (t571 * (mrSges(3,2) * t549 - mrSges(3,3) * t538 + t584 * t414 - t580 * t415 - t418 * t625) + t575 * (-mrSges(3,1) * t549 + mrSges(3,3) * t539 - pkin(2) * t418 - t572 * t413) - pkin(1) * t417 + qJ(2) * (-t416 * t571 + t422 * t575) + (-t419 * t626 + t575 * t593) * t576 + ((Ifges(3,2) * t575 ^ 2 + (Ifges(3,1) * t571 + Ifges(3,4) * t628) * t571) * t573 + 0.2e1 * t577 * (Ifges(3,5) * t571 + Ifges(3,6) * t575)) * qJDD(1)) * t573; t417; t413; t629; t445; t589;];
tauJ  = t1;
