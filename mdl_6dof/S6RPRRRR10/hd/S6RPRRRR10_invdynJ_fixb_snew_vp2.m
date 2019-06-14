% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRR10
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
% Datum: 2019-05-06 04:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:44:55
% EndTime: 2019-05-06 04:45:22
% DurationCPUTime: 26.13s
% Computational Cost: add. (389142->358), mult. (1206665->490), div. (0->0), fcn. (1043579->16), ass. (0->171)
t586 = sin(pkin(13));
t588 = sin(pkin(6));
t589 = cos(pkin(13));
t591 = cos(pkin(6));
t600 = cos(qJ(3));
t590 = cos(pkin(7));
t595 = sin(qJ(3));
t633 = t590 * t595;
t587 = sin(pkin(7));
t636 = t587 * t595;
t607 = t591 * t636 + (t586 * t600 + t589 * t633) * t588;
t562 = t607 * qJD(1);
t634 = t588 * t590;
t615 = t587 * t591 + t589 * t634;
t611 = t615 * t600;
t637 = t586 * t588;
t630 = t595 * t637;
t548 = -t562 * qJD(3) + (t611 - t630) * qJDD(1);
t610 = t615 * qJD(1);
t573 = pkin(9) * t610;
t602 = qJD(1) ^ 2;
t596 = sin(qJ(1));
t601 = cos(qJ(1));
t622 = -g(1) * t601 - g(2) * t596;
t640 = qJ(2) * t588;
t577 = -pkin(1) * t602 + qJDD(1) * t640 + t622;
t642 = pkin(9) * t587;
t618 = -pkin(2) * t589 - t586 * t642;
t631 = qJD(1) * t588;
t641 = pkin(9) * qJDD(1);
t613 = qJD(1) * t618 * t631 + t590 * t641;
t628 = qJD(2) * t631;
t635 = t588 * t589;
t627 = t596 * g(1) - g(2) * t601;
t576 = qJDD(1) * pkin(1) + t602 * t640 + t627;
t638 = t576 * t591;
t619 = -g(3) * t635 - 0.2e1 * t586 * t628 + t589 * t638;
t527 = (pkin(2) * qJDD(1) + qJD(1) * t573) * t591 + (-t588 * t613 - t577) * t586 + t619;
t643 = pkin(9) * t586;
t578 = (pkin(2) * t591 - t634 * t643) * qJD(1);
t644 = 0.2e1 * t589;
t629 = t589 * t577 + t586 * t638 + t628 * t644;
t528 = (-qJD(1) * t578 + t587 * t641) * t591 + (-g(3) * t586 + t589 * t613) * t588 + t629;
t626 = -t591 * g(3) + qJDD(2);
t536 = (-t576 + t618 * qJDD(1) + (-t573 * t589 + t578 * t586) * qJD(1)) * t588 + t626;
t500 = -t595 * t528 + (t527 * t590 + t536 * t587) * t600;
t582 = qJD(1) * t630;
t561 = t600 * t610 - t582;
t546 = -mrSges(4,1) * t561 + mrSges(4,2) * t562;
t549 = t561 * qJD(3) + qJDD(1) * t607;
t614 = -t587 * t635 + t590 * t591;
t574 = qJD(1) * t614 + qJD(3);
t555 = -mrSges(4,2) * t574 + mrSges(4,3) * t561;
t571 = qJDD(1) * t614 + qJDD(3);
t547 = -pkin(3) * t561 - pkin(10) * t562;
t570 = t574 ^ 2;
t484 = -t571 * pkin(3) - t570 * pkin(10) + t562 * t547 - t500;
t594 = sin(qJ(4));
t599 = cos(qJ(4));
t554 = t562 * t599 + t574 * t594;
t522 = -qJD(4) * t554 - t549 * t594 + t571 * t599;
t553 = -t562 * t594 + t574 * t599;
t523 = qJD(4) * t553 + t549 * t599 + t571 * t594;
t560 = -qJD(1) * t611 + qJD(4) + t582;
t537 = -mrSges(5,2) * t560 + mrSges(5,3) * t553;
t538 = mrSges(5,1) * t560 - mrSges(5,3) * t554;
t501 = t527 * t633 + t600 * t528 + t536 * t636;
t485 = -pkin(3) * t570 + pkin(10) * t571 + t547 * t561 + t501;
t511 = -t587 * t527 + t590 * t536;
t488 = (-t561 * t574 - t549) * pkin(10) + (t562 * t574 - t548) * pkin(3) + t511;
t477 = -t594 * t485 + t599 * t488;
t545 = qJDD(4) - t548;
t474 = (t553 * t560 - t523) * pkin(11) + (t553 * t554 + t545) * pkin(4) + t477;
t478 = t599 * t485 + t594 * t488;
t539 = pkin(4) * t560 - pkin(11) * t554;
t552 = t553 ^ 2;
t476 = -pkin(4) * t552 + pkin(11) * t522 - t539 * t560 + t478;
t593 = sin(qJ(5));
t598 = cos(qJ(5));
t471 = t593 * t474 + t598 * t476;
t530 = t553 * t598 - t554 * t593;
t531 = t553 * t593 + t554 * t598;
t510 = -pkin(5) * t530 - pkin(12) * t531;
t544 = qJDD(5) + t545;
t558 = qJD(5) + t560;
t557 = t558 ^ 2;
t468 = -pkin(5) * t557 + pkin(12) * t544 + t510 * t530 + t471;
t479 = -t522 * pkin(4) - t552 * pkin(11) + t554 * t539 + t484;
t497 = -qJD(5) * t531 + t522 * t598 - t523 * t593;
t498 = qJD(5) * t530 + t522 * t593 + t523 * t598;
t472 = (-t530 * t558 - t498) * pkin(12) + (t531 * t558 - t497) * pkin(5) + t479;
t592 = sin(qJ(6));
t597 = cos(qJ(6));
t465 = -t468 * t592 + t472 * t597;
t512 = -t531 * t592 + t558 * t597;
t482 = qJD(6) * t512 + t498 * t597 + t544 * t592;
t496 = qJDD(6) - t497;
t513 = t531 * t597 + t558 * t592;
t502 = -mrSges(7,1) * t512 + mrSges(7,2) * t513;
t529 = qJD(6) - t530;
t503 = -mrSges(7,2) * t529 + mrSges(7,3) * t512;
t461 = m(7) * t465 + mrSges(7,1) * t496 - mrSges(7,3) * t482 - t502 * t513 + t503 * t529;
t466 = t468 * t597 + t472 * t592;
t481 = -qJD(6) * t513 - t498 * t592 + t544 * t597;
t504 = mrSges(7,1) * t529 - mrSges(7,3) * t513;
t462 = m(7) * t466 - mrSges(7,2) * t496 + mrSges(7,3) * t481 + t502 * t512 - t504 * t529;
t450 = t597 * t461 + t592 * t462;
t514 = -mrSges(6,2) * t558 + mrSges(6,3) * t530;
t515 = mrSges(6,1) * t558 - mrSges(6,3) * t531;
t608 = m(6) * t479 - t497 * mrSges(6,1) + mrSges(6,2) * t498 - t530 * t514 + t515 * t531 + t450;
t604 = -m(5) * t484 + t522 * mrSges(5,1) - mrSges(5,2) * t523 + t553 * t537 - t538 * t554 - t608;
t445 = m(4) * t500 + mrSges(4,1) * t571 - mrSges(4,3) * t549 - t546 * t562 + t555 * t574 + t604;
t639 = t445 * t600;
t509 = -mrSges(6,1) * t530 + mrSges(6,2) * t531;
t623 = -t461 * t592 + t597 * t462;
t448 = m(6) * t471 - mrSges(6,2) * t544 + mrSges(6,3) * t497 + t509 * t530 - t515 * t558 + t623;
t470 = t474 * t598 - t476 * t593;
t467 = -pkin(5) * t544 - pkin(12) * t557 + t510 * t531 - t470;
t609 = -m(7) * t467 + t481 * mrSges(7,1) - mrSges(7,2) * t482 + t512 * t503 - t504 * t513;
t457 = m(6) * t470 + mrSges(6,1) * t544 - mrSges(6,3) * t498 - t509 * t531 + t514 * t558 + t609;
t442 = t593 * t448 + t598 * t457;
t532 = -mrSges(5,1) * t553 + mrSges(5,2) * t554;
t440 = m(5) * t477 + mrSges(5,1) * t545 - mrSges(5,3) * t523 - t532 * t554 + t537 * t560 + t442;
t624 = t598 * t448 - t457 * t593;
t441 = m(5) * t478 - mrSges(5,2) * t545 + mrSges(5,3) * t522 + t532 * t553 - t538 * t560 + t624;
t434 = t599 * t440 + t594 * t441;
t556 = mrSges(4,1) * t574 - mrSges(4,3) * t562;
t625 = -t440 * t594 + t599 * t441;
t431 = m(4) * t501 - mrSges(4,2) * t571 + mrSges(4,3) * t548 + t546 * t561 - t556 * t574 + t625;
t433 = m(4) * t511 - mrSges(4,1) * t548 + mrSges(4,2) * t549 - t555 * t561 + t556 * t562 + t434;
t422 = t431 * t636 + t590 * t433 + t587 * t639;
t427 = t600 * t431 - t595 * t445;
t423 = t431 * t633 - t587 * t433 + t590 * t639;
t621 = -mrSges(3,1) * t589 + mrSges(3,2) * t586;
t617 = mrSges(3,1) * t591 - mrSges(3,3) * t637;
t616 = -mrSges(3,2) * t591 + mrSges(3,3) * t635;
t489 = Ifges(7,5) * t513 + Ifges(7,6) * t512 + Ifges(7,3) * t529;
t491 = Ifges(7,1) * t513 + Ifges(7,4) * t512 + Ifges(7,5) * t529;
t454 = -mrSges(7,1) * t467 + mrSges(7,3) * t466 + Ifges(7,4) * t482 + Ifges(7,2) * t481 + Ifges(7,6) * t496 - t489 * t513 + t491 * t529;
t490 = Ifges(7,4) * t513 + Ifges(7,2) * t512 + Ifges(7,6) * t529;
t455 = mrSges(7,2) * t467 - mrSges(7,3) * t465 + Ifges(7,1) * t482 + Ifges(7,4) * t481 + Ifges(7,5) * t496 + t489 * t512 - t490 * t529;
t505 = Ifges(6,5) * t531 + Ifges(6,6) * t530 + Ifges(6,3) * t558;
t506 = Ifges(6,4) * t531 + Ifges(6,2) * t530 + Ifges(6,6) * t558;
t435 = mrSges(6,2) * t479 - mrSges(6,3) * t470 + Ifges(6,1) * t498 + Ifges(6,4) * t497 + Ifges(6,5) * t544 - pkin(12) * t450 - t454 * t592 + t455 * t597 + t505 * t530 - t506 * t558;
t507 = Ifges(6,1) * t531 + Ifges(6,4) * t530 + Ifges(6,5) * t558;
t605 = mrSges(7,1) * t465 - mrSges(7,2) * t466 + Ifges(7,5) * t482 + Ifges(7,6) * t481 + Ifges(7,3) * t496 + t490 * t513 - t491 * t512;
t436 = -mrSges(6,1) * t479 + mrSges(6,3) * t471 + Ifges(6,4) * t498 + Ifges(6,2) * t497 + Ifges(6,6) * t544 - pkin(5) * t450 - t505 * t531 + t507 * t558 - t605;
t516 = Ifges(5,5) * t554 + Ifges(5,6) * t553 + Ifges(5,3) * t560;
t518 = Ifges(5,1) * t554 + Ifges(5,4) * t553 + Ifges(5,5) * t560;
t424 = -mrSges(5,1) * t484 + mrSges(5,3) * t478 + Ifges(5,4) * t523 + Ifges(5,2) * t522 + Ifges(5,6) * t545 - pkin(4) * t608 + pkin(11) * t624 + t593 * t435 + t598 * t436 - t554 * t516 + t560 * t518;
t517 = Ifges(5,4) * t554 + Ifges(5,2) * t553 + Ifges(5,6) * t560;
t425 = mrSges(5,2) * t484 - mrSges(5,3) * t477 + Ifges(5,1) * t523 + Ifges(5,4) * t522 + Ifges(5,5) * t545 - pkin(11) * t442 + t435 * t598 - t436 * t593 + t516 * t553 - t517 * t560;
t540 = Ifges(4,5) * t562 + Ifges(4,6) * t561 + Ifges(4,3) * t574;
t541 = Ifges(4,4) * t562 + Ifges(4,2) * t561 + Ifges(4,6) * t574;
t418 = mrSges(4,2) * t511 - mrSges(4,3) * t500 + Ifges(4,1) * t549 + Ifges(4,4) * t548 + Ifges(4,5) * t571 - pkin(10) * t434 - t424 * t594 + t425 * t599 + t540 * t561 - t541 * t574;
t542 = Ifges(4,1) * t562 + Ifges(4,4) * t561 + Ifges(4,5) * t574;
t606 = -mrSges(6,1) * t470 + mrSges(6,2) * t471 - Ifges(6,5) * t498 - Ifges(6,6) * t497 - Ifges(6,3) * t544 - pkin(5) * t609 - pkin(12) * t623 - t597 * t454 - t592 * t455 - t531 * t506 + t530 * t507;
t603 = mrSges(5,1) * t477 - mrSges(5,2) * t478 + Ifges(5,5) * t523 + Ifges(5,6) * t522 + Ifges(5,3) * t545 + pkin(4) * t442 + t554 * t517 - t553 * t518 - t606;
t419 = -mrSges(4,1) * t511 + mrSges(4,3) * t501 + Ifges(4,4) * t549 + Ifges(4,2) * t548 + Ifges(4,6) * t571 - pkin(3) * t434 - t562 * t540 + t574 * t542 - t603;
t612 = pkin(9) * t427 + t418 * t595 + t419 * t600;
t580 = t616 * qJD(1);
t579 = t617 * qJD(1);
t575 = t621 * t631;
t563 = -t588 * t576 + t626;
t551 = -g(3) * t637 + t629;
t550 = -t586 * t577 + t619;
t426 = m(3) * t551 + t616 * qJDD(1) + (t575 * t635 - t579 * t591) * qJD(1) + t427;
t421 = m(3) * t563 + (t621 * qJDD(1) + (t579 * t586 - t580 * t589) * qJD(1)) * t588 + t422;
t420 = m(3) * t550 + t617 * qJDD(1) + (-t575 * t637 + t580 * t591) * qJD(1) + t423;
t417 = mrSges(4,1) * t500 - mrSges(4,2) * t501 + Ifges(4,5) * t549 + Ifges(4,6) * t548 + Ifges(4,3) * t571 + pkin(3) * t604 + pkin(10) * t625 + t599 * t424 + t594 * t425 + t562 * t541 - t561 * t542;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t627 - mrSges(2,2) * t622 + (mrSges(3,1) * t550 - mrSges(3,2) * t551 + pkin(2) * t423 + t590 * t417 + pkin(1) * (t420 * t589 + t426 * t586) + Ifges(3,3) * t591 * qJDD(1) + t612 * t587) * t591 + (t586 * (mrSges(3,2) * t563 - mrSges(3,3) * t550 + t600 * t418 - t595 * t419 - t422 * t642) + t589 * (-mrSges(3,1) * t563 + mrSges(3,3) * t551 - pkin(2) * t422 - t587 * t417) - pkin(1) * t421 + qJ(2) * (-t420 * t586 + t426 * t589) + (-t423 * t643 + t589 * t612) * t590 + ((Ifges(3,2) * t589 ^ 2 + (Ifges(3,1) * t586 + Ifges(3,4) * t644) * t586) * t588 + 0.2e1 * t591 * (Ifges(3,5) * t586 + Ifges(3,6) * t589)) * qJDD(1)) * t588; t421; t417; t603; -t606; t605;];
tauJ  = t1;
