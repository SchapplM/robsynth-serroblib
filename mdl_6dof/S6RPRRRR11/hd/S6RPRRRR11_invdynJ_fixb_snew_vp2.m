% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 05:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR11_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR11_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 05:34:22
% EndTime: 2019-05-06 05:34:50
% DurationCPUTime: 27.92s
% Computational Cost: add. (404353->357), mult. (1253819->491), div. (0->0), fcn. (1081514->16), ass. (0->168)
t573 = sin(pkin(13));
t575 = sin(pkin(6));
t576 = cos(pkin(13));
t578 = cos(pkin(6));
t582 = sin(qJ(3));
t577 = cos(pkin(7));
t587 = cos(qJ(3));
t618 = t577 * t587;
t574 = sin(pkin(7));
t622 = t574 * t587;
t593 = t575 * (-t573 * t582 + t576 * t618) + t578 * t622;
t550 = t593 * qJD(1);
t619 = t577 * t582;
t623 = t574 * t582;
t595 = t578 * t623 + (t573 * t587 + t576 * t619) * t575;
t551 = t595 * qJD(1);
t536 = -qJD(3) * t551 + qJDD(1) * t593;
t620 = t575 * t577;
t562 = (t574 * t578 + t576 * t620) * qJD(1) * pkin(9);
t589 = qJD(1) ^ 2;
t583 = sin(qJ(1));
t588 = cos(qJ(1));
t608 = -g(1) * t588 - g(2) * t583;
t626 = qJ(2) * t575;
t566 = -pkin(1) * t589 + qJDD(1) * t626 + t608;
t628 = pkin(9) * t574;
t604 = -pkin(2) * t576 - t573 * t628;
t617 = qJD(1) * t575;
t627 = pkin(9) * qJDD(1);
t599 = qJD(1) * t604 * t617 + t577 * t627;
t613 = qJD(2) * t617;
t621 = t575 * t576;
t612 = t583 * g(1) - g(2) * t588;
t565 = qJDD(1) * pkin(1) + t589 * t626 + t612;
t625 = t565 * t578;
t605 = -g(3) * t621 - 0.2e1 * t573 * t613 + t576 * t625;
t516 = (pkin(2) * qJDD(1) + qJD(1) * t562) * t578 + (-t575 * t599 - t566) * t573 + t605;
t629 = pkin(9) * t573;
t567 = (pkin(2) * t578 - t620 * t629) * qJD(1);
t630 = 0.2e1 * t576;
t614 = t576 * t566 + t573 * t625 + t613 * t630;
t517 = (-qJD(1) * t567 + t574 * t627) * t578 + (-g(3) * t573 + t576 * t599) * t575 + t614;
t611 = -g(3) * t578 + qJDD(2);
t526 = (-t565 + t604 * qJDD(1) + (-t562 * t576 + t567 * t573) * qJD(1)) * t575 + t611;
t484 = -t582 * t517 + (t516 * t577 + t526 * t574) * t587;
t624 = t573 * t575;
t485 = t516 * t619 + t587 * t517 + t526 * t623;
t535 = -pkin(3) * t550 - pkin(10) * t551;
t600 = -t574 * t621 + t577 * t578;
t563 = qJD(1) * t600 + qJD(3);
t559 = t563 ^ 2;
t560 = qJDD(1) * t600 + qJDD(3);
t476 = -pkin(3) * t559 + pkin(10) * t560 + t535 * t550 + t485;
t495 = -t516 * t574 + t577 * t526;
t537 = qJD(3) * t550 + qJDD(1) * t595;
t481 = (-t550 * t563 - t537) * pkin(10) + (t551 * t563 - t536) * pkin(3) + t495;
t581 = sin(qJ(4));
t586 = cos(qJ(4));
t466 = t586 * t476 + t581 * t481;
t543 = -t581 * t551 + t563 * t586;
t544 = t551 * t586 + t563 * t581;
t519 = -pkin(4) * t543 - pkin(11) * t544;
t533 = qJDD(4) - t536;
t549 = qJD(4) - t550;
t548 = t549 ^ 2;
t461 = -pkin(4) * t548 + pkin(11) * t533 + t519 * t543 + t466;
t475 = -pkin(3) * t560 - pkin(10) * t559 + t551 * t535 - t484;
t511 = -t544 * qJD(4) - t581 * t537 + t560 * t586;
t512 = qJD(4) * t543 + t537 * t586 + t560 * t581;
t464 = (-t543 * t549 - t512) * pkin(11) + (t544 * t549 - t511) * pkin(4) + t475;
t580 = sin(qJ(5));
t585 = cos(qJ(5));
t456 = -t461 * t580 + t585 * t464;
t524 = -t544 * t580 + t549 * t585;
t488 = qJD(5) * t524 + t512 * t585 + t533 * t580;
t509 = qJDD(5) - t511;
t525 = t544 * t585 + t549 * t580;
t540 = qJD(5) - t543;
t454 = (t524 * t540 - t488) * pkin(12) + (t524 * t525 + t509) * pkin(5) + t456;
t457 = t585 * t461 + t580 * t464;
t487 = -qJD(5) * t525 - t512 * t580 + t533 * t585;
t502 = pkin(5) * t540 - pkin(12) * t525;
t523 = t524 ^ 2;
t455 = -pkin(5) * t523 + pkin(12) * t487 - t502 * t540 + t457;
t579 = sin(qJ(6));
t584 = cos(qJ(6));
t452 = t454 * t584 - t455 * t579;
t496 = t524 * t584 - t525 * t579;
t471 = qJD(6) * t496 + t487 * t579 + t488 * t584;
t497 = t524 * t579 + t525 * t584;
t483 = -mrSges(7,1) * t496 + mrSges(7,2) * t497;
t538 = qJD(6) + t540;
t489 = -mrSges(7,2) * t538 + mrSges(7,3) * t496;
t504 = qJDD(6) + t509;
t448 = m(7) * t452 + mrSges(7,1) * t504 - mrSges(7,3) * t471 - t483 * t497 + t489 * t538;
t453 = t454 * t579 + t455 * t584;
t470 = -qJD(6) * t497 + t487 * t584 - t488 * t579;
t490 = mrSges(7,1) * t538 - mrSges(7,3) * t497;
t449 = m(7) * t453 - mrSges(7,2) * t504 + mrSges(7,3) * t470 + t483 * t496 - t490 * t538;
t440 = t584 * t448 + t579 * t449;
t498 = -mrSges(6,1) * t524 + mrSges(6,2) * t525;
t500 = -mrSges(6,2) * t540 + mrSges(6,3) * t524;
t438 = m(6) * t456 + mrSges(6,1) * t509 - mrSges(6,3) * t488 - t498 * t525 + t500 * t540 + t440;
t501 = mrSges(6,1) * t540 - mrSges(6,3) * t525;
t609 = -t448 * t579 + t584 * t449;
t439 = m(6) * t457 - mrSges(6,2) * t509 + mrSges(6,3) * t487 + t498 * t524 - t501 * t540 + t609;
t436 = -t438 * t580 + t585 * t439;
t518 = -mrSges(5,1) * t543 + mrSges(5,2) * t544;
t528 = mrSges(5,1) * t549 - mrSges(5,3) * t544;
t434 = m(5) * t466 - mrSges(5,2) * t533 + mrSges(5,3) * t511 + t518 * t543 - t528 * t549 + t436;
t465 = -t581 * t476 + t481 * t586;
t460 = -pkin(4) * t533 - pkin(11) * t548 + t544 * t519 - t465;
t458 = -pkin(5) * t487 - pkin(12) * t523 + t502 * t525 + t460;
t597 = m(7) * t458 - t470 * mrSges(7,1) + mrSges(7,2) * t471 - t496 * t489 + t490 * t497;
t450 = -m(6) * t460 + t487 * mrSges(6,1) - mrSges(6,2) * t488 + t524 * t500 - t501 * t525 - t597;
t527 = -mrSges(5,2) * t549 + mrSges(5,3) * t543;
t444 = m(5) * t465 + mrSges(5,1) * t533 - mrSges(5,3) * t512 - t518 * t544 + t527 * t549 + t450;
t426 = t581 * t434 + t586 * t444;
t534 = -mrSges(4,1) * t550 + mrSges(4,2) * t551;
t546 = mrSges(4,1) * t563 - mrSges(4,3) * t551;
t610 = t586 * t434 - t444 * t581;
t423 = m(4) * t485 - mrSges(4,2) * t560 + mrSges(4,3) * t536 + t534 * t550 - t546 * t563 + t610;
t545 = -mrSges(4,2) * t563 + mrSges(4,3) * t550;
t425 = m(4) * t495 - mrSges(4,1) * t536 + mrSges(4,2) * t537 - t545 * t550 + t546 * t551 + t426;
t435 = t438 * t585 + t439 * t580;
t592 = -m(5) * t475 + t511 * mrSges(5,1) - mrSges(5,2) * t512 + t543 * t527 - t528 * t544 - t435;
t431 = m(4) * t484 + mrSges(4,1) * t560 - mrSges(4,3) * t537 - t534 * t551 + t545 * t563 + t592;
t414 = t423 * t623 + t577 * t425 + t431 * t622;
t418 = t587 * t423 - t431 * t582;
t415 = t423 * t619 - t425 * t574 + t431 * t618;
t607 = -mrSges(3,1) * t576 + mrSges(3,2) * t573;
t603 = mrSges(3,1) * t578 - mrSges(3,3) * t624;
t602 = -mrSges(3,2) * t578 + mrSges(3,3) * t621;
t478 = Ifges(7,5) * t497 + Ifges(7,6) * t496 + Ifges(7,3) * t538;
t480 = Ifges(7,1) * t497 + Ifges(7,4) * t496 + Ifges(7,5) * t538;
t441 = -mrSges(7,1) * t458 + mrSges(7,3) * t453 + Ifges(7,4) * t471 + Ifges(7,2) * t470 + Ifges(7,6) * t504 - t478 * t497 + t480 * t538;
t479 = Ifges(7,4) * t497 + Ifges(7,2) * t496 + Ifges(7,6) * t538;
t442 = mrSges(7,2) * t458 - mrSges(7,3) * t452 + Ifges(7,1) * t471 + Ifges(7,4) * t470 + Ifges(7,5) * t504 + t478 * t496 - t479 * t538;
t491 = Ifges(6,5) * t525 + Ifges(6,6) * t524 + Ifges(6,3) * t540;
t493 = Ifges(6,1) * t525 + Ifges(6,4) * t524 + Ifges(6,5) * t540;
t427 = -mrSges(6,1) * t460 + mrSges(6,3) * t457 + Ifges(6,4) * t488 + Ifges(6,2) * t487 + Ifges(6,6) * t509 - pkin(5) * t597 + pkin(12) * t609 + t584 * t441 + t579 * t442 - t525 * t491 + t540 * t493;
t492 = Ifges(6,4) * t525 + Ifges(6,2) * t524 + Ifges(6,6) * t540;
t428 = mrSges(6,2) * t460 - mrSges(6,3) * t456 + Ifges(6,1) * t488 + Ifges(6,4) * t487 + Ifges(6,5) * t509 - pkin(12) * t440 - t441 * t579 + t442 * t584 + t491 * t524 - t492 * t540;
t505 = Ifges(5,5) * t544 + Ifges(5,6) * t543 + Ifges(5,3) * t549;
t506 = Ifges(5,4) * t544 + Ifges(5,2) * t543 + Ifges(5,6) * t549;
t416 = mrSges(5,2) * t475 - mrSges(5,3) * t465 + Ifges(5,1) * t512 + Ifges(5,4) * t511 + Ifges(5,5) * t533 - pkin(11) * t435 - t427 * t580 + t428 * t585 + t505 * t543 - t506 * t549;
t507 = Ifges(5,1) * t544 + Ifges(5,4) * t543 + Ifges(5,5) * t549;
t596 = -mrSges(7,1) * t452 + mrSges(7,2) * t453 - Ifges(7,5) * t471 - Ifges(7,6) * t470 - Ifges(7,3) * t504 - t497 * t479 + t496 * t480;
t590 = mrSges(6,1) * t456 - mrSges(6,2) * t457 + Ifges(6,5) * t488 + Ifges(6,6) * t487 + Ifges(6,3) * t509 + pkin(5) * t440 + t525 * t492 - t524 * t493 - t596;
t419 = -mrSges(5,1) * t475 + mrSges(5,3) * t466 + Ifges(5,4) * t512 + Ifges(5,2) * t511 + Ifges(5,6) * t533 - pkin(4) * t435 - t544 * t505 + t549 * t507 - t590;
t529 = Ifges(4,5) * t551 + Ifges(4,6) * t550 + Ifges(4,3) * t563;
t530 = Ifges(4,4) * t551 + Ifges(4,2) * t550 + Ifges(4,6) * t563;
t410 = mrSges(4,2) * t495 - mrSges(4,3) * t484 + Ifges(4,1) * t537 + Ifges(4,4) * t536 + Ifges(4,5) * t560 - pkin(10) * t426 + t416 * t586 - t419 * t581 + t529 * t550 - t530 * t563;
t531 = Ifges(4,1) * t551 + Ifges(4,4) * t550 + Ifges(4,5) * t563;
t591 = mrSges(5,1) * t465 - mrSges(5,2) * t466 + Ifges(5,5) * t512 + Ifges(5,6) * t511 + Ifges(5,3) * t533 + pkin(4) * t450 + pkin(11) * t436 + t585 * t427 + t580 * t428 + t544 * t506 - t543 * t507;
t411 = -mrSges(4,1) * t495 + mrSges(4,3) * t485 + Ifges(4,4) * t537 + Ifges(4,2) * t536 + Ifges(4,6) * t560 - pkin(3) * t426 - t551 * t529 + t563 * t531 - t591;
t598 = pkin(9) * t418 + t410 * t582 + t411 * t587;
t569 = t602 * qJD(1);
t568 = t603 * qJD(1);
t564 = t607 * t617;
t552 = -t565 * t575 + t611;
t542 = -g(3) * t624 + t614;
t541 = -t566 * t573 + t605;
t417 = m(3) * t542 + t602 * qJDD(1) + (t564 * t621 - t568 * t578) * qJD(1) + t418;
t413 = m(3) * t552 + (t607 * qJDD(1) + (t568 * t573 - t569 * t576) * qJD(1)) * t575 + t414;
t412 = m(3) * t541 + t603 * qJDD(1) + (-t564 * t624 + t569 * t578) * qJD(1) + t415;
t409 = mrSges(4,1) * t484 - mrSges(4,2) * t485 + Ifges(4,5) * t537 + Ifges(4,6) * t536 + Ifges(4,3) * t560 + pkin(3) * t592 + pkin(10) * t610 + t581 * t416 + t586 * t419 + t551 * t530 - t550 * t531;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t612 - mrSges(2,2) * t608 + (mrSges(3,1) * t541 - mrSges(3,2) * t542 + pkin(2) * t415 + t577 * t409 + pkin(1) * (t412 * t576 + t417 * t573) + Ifges(3,3) * t578 * qJDD(1) + t598 * t574) * t578 + (t573 * (mrSges(3,2) * t552 - mrSges(3,3) * t541 + t587 * t410 - t582 * t411 - t414 * t628) + t576 * (-mrSges(3,1) * t552 + mrSges(3,3) * t542 - pkin(2) * t414 - t574 * t409) - pkin(1) * t413 + qJ(2) * (-t412 * t573 + t417 * t576) + (-t415 * t629 + t576 * t598) * t577 + ((Ifges(3,2) * t576 ^ 2 + (Ifges(3,1) * t573 + Ifges(3,4) * t630) * t573) * t575 + 0.2e1 * t578 * (Ifges(3,5) * t573 + Ifges(3,6) * t576)) * qJDD(1)) * t575; t413; t409; t591; t590; -t596;];
tauJ  = t1;
