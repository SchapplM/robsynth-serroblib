% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRR12
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 07:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR12_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_invdynJ_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR12_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 06:31:38
% EndTime: 2019-05-06 06:32:42
% DurationCPUTime: 63.64s
% Computational Cost: add. (970507->371), mult. (3165314->521), div. (0->0), fcn. (2783205->18), ass. (0->181)
t587 = sin(pkin(7));
t589 = cos(pkin(14));
t592 = cos(pkin(6));
t588 = sin(pkin(6));
t591 = cos(pkin(7));
t632 = t588 * t591;
t574 = (t587 * t592 + t589 * t632) * qJD(1) * pkin(10);
t603 = qJD(1) ^ 2;
t597 = sin(qJ(1));
t602 = cos(qJ(1));
t622 = -g(1) * t602 - g(2) * t597;
t641 = qJ(2) * t588;
t578 = -pkin(1) * t603 + qJDD(1) * t641 + t622;
t585 = sin(pkin(14));
t645 = pkin(10) * t587;
t615 = -pkin(2) * t589 - t585 * t645;
t628 = qJD(1) * t588;
t642 = pkin(10) * qJDD(1);
t611 = qJD(1) * t615 * t628 + t591 * t642;
t626 = qJD(2) * t628;
t633 = t588 * t589;
t625 = t597 * g(1) - g(2) * t602;
t577 = qJDD(1) * pkin(1) + t603 * t641 + t625;
t638 = t577 * t592;
t616 = -g(3) * t633 - 0.2e1 * t585 * t626 + t589 * t638;
t537 = (pkin(2) * qJDD(1) + qJD(1) * t574) * t592 + (-t588 * t611 - t578) * t585 + t616;
t646 = pkin(10) * t585;
t579 = (pkin(2) * t592 - t632 * t646) * qJD(1);
t647 = 0.2e1 * t589;
t627 = t589 * t578 + t585 * t638 + t626 * t647;
t538 = (-qJD(1) * t579 + t587 * t642) * t592 + (-g(3) * t585 + t589 * t611) * t588 + t627;
t624 = -g(3) * t592 + qJDD(2);
t546 = (-t577 + t615 * qJDD(1) + (-t574 * t589 + t579 * t585) * qJD(1)) * t588 + t624;
t596 = sin(qJ(3));
t601 = cos(qJ(3));
t629 = t591 * t601;
t634 = t587 * t601;
t511 = t537 * t629 - t538 * t596 + t546 * t634;
t607 = t592 * t634 + (-t585 * t596 + t589 * t629) * t588;
t563 = t607 * qJD(1);
t630 = t591 * t596;
t635 = t587 * t596;
t608 = t592 * t635 + (t585 * t601 + t589 * t630) * t588;
t564 = t608 * qJD(1);
t586 = sin(pkin(8));
t644 = pkin(11) * t586;
t551 = -pkin(3) * t563 - t564 * t644;
t554 = qJD(3) * t563 + qJDD(1) * t608;
t612 = -t587 * t633 + t591 * t592;
t575 = qJD(1) * t612 + qJD(3);
t590 = cos(pkin(8));
t618 = t563 * t590 + t575 * t586;
t556 = t618 * pkin(11);
t572 = qJDD(1) * t612 + qJDD(3);
t643 = pkin(11) * t590;
t495 = pkin(3) * t572 - t551 * t564 - t554 * t643 + t556 * t575 + t511;
t512 = t537 * t630 + t601 * t538 + t546 * t635;
t560 = pkin(3) * t575 - t564 * t643;
t553 = -qJD(3) * t564 + qJDD(1) * t607;
t619 = t553 * t590 + t572 * t586;
t496 = pkin(11) * t619 + t551 * t563 - t560 * t575 + t512;
t526 = -t537 * t587 + t591 * t546;
t499 = -pkin(3) * t553 - t554 * t644 - t556 * t563 + t560 * t564 + t526;
t595 = sin(qJ(4));
t600 = cos(qJ(4));
t482 = -t595 * t496 + (t495 * t590 + t499 * t586) * t600;
t545 = t564 * t600 + t595 * t618;
t519 = -qJD(4) * t545 - t554 * t595 + t600 * t619;
t544 = -t564 * t595 + t600 * t618;
t520 = qJD(4) * t544 + t554 * t600 + t595 * t619;
t527 = -mrSges(5,1) * t544 + mrSges(5,2) * t545;
t557 = -t563 * t586 + t575 * t590 + qJD(4);
t532 = -mrSges(5,2) * t557 + mrSges(5,3) * t544;
t547 = -t553 * t586 + t572 * t590 + qJDD(4);
t631 = t590 * t595;
t636 = t586 * t595;
t483 = t495 * t631 + t600 * t496 + t499 * t636;
t528 = -pkin(4) * t544 - pkin(12) * t545;
t555 = t557 ^ 2;
t479 = -pkin(4) * t555 + pkin(12) * t547 + t528 * t544 + t483;
t484 = -t495 * t586 + t590 * t499;
t481 = (-t544 * t557 - t520) * pkin(12) + (t545 * t557 - t519) * pkin(4) + t484;
t594 = sin(qJ(5));
t599 = cos(qJ(5));
t475 = t599 * t479 + t594 * t481;
t530 = -t545 * t594 + t557 * t599;
t531 = t545 * t599 + t557 * t594;
t514 = -pkin(5) * t530 - pkin(13) * t531;
t518 = qJDD(5) - t519;
t542 = qJD(5) - t544;
t541 = t542 ^ 2;
t473 = -pkin(5) * t541 + pkin(13) * t518 + t514 * t530 + t475;
t478 = -pkin(4) * t547 - pkin(12) * t555 + t545 * t528 - t482;
t503 = -qJD(5) * t531 - t520 * t594 + t547 * t599;
t504 = qJD(5) * t530 + t520 * t599 + t547 * t594;
t476 = (-t530 * t542 - t504) * pkin(13) + (t531 * t542 - t503) * pkin(5) + t478;
t593 = sin(qJ(6));
t598 = cos(qJ(6));
t469 = -t473 * t593 + t476 * t598;
t516 = -t531 * t593 + t542 * t598;
t487 = qJD(6) * t516 + t504 * t598 + t518 * t593;
t517 = t531 * t598 + t542 * t593;
t500 = -mrSges(7,1) * t516 + mrSges(7,2) * t517;
t502 = qJDD(6) - t503;
t529 = qJD(6) - t530;
t505 = -mrSges(7,2) * t529 + mrSges(7,3) * t516;
t467 = m(7) * t469 + mrSges(7,1) * t502 - mrSges(7,3) * t487 - t500 * t517 + t505 * t529;
t470 = t473 * t598 + t476 * t593;
t486 = -qJD(6) * t517 - t504 * t593 + t518 * t598;
t506 = mrSges(7,1) * t529 - mrSges(7,3) * t517;
t468 = m(7) * t470 - mrSges(7,2) * t502 + mrSges(7,3) * t486 + t500 * t516 - t506 * t529;
t460 = t467 * t598 + t468 * t593;
t521 = -mrSges(6,2) * t542 + mrSges(6,3) * t530;
t522 = mrSges(6,1) * t542 - mrSges(6,3) * t531;
t606 = -m(6) * t478 + t503 * mrSges(6,1) - mrSges(6,2) * t504 + t530 * t521 - t522 * t531 - t460;
t456 = m(5) * t482 + mrSges(5,1) * t547 - mrSges(5,3) * t520 - t527 * t545 + t532 * t557 + t606;
t640 = t456 * t600;
t637 = t585 * t588;
t461 = -t467 * t593 + t598 * t468;
t513 = -mrSges(6,1) * t530 + mrSges(6,2) * t531;
t459 = m(6) * t475 - mrSges(6,2) * t518 + mrSges(6,3) * t503 + t513 * t530 - t522 * t542 + t461;
t474 = -t479 * t594 + t481 * t599;
t472 = -pkin(5) * t518 - pkin(13) * t541 + t514 * t531 - t474;
t471 = -m(7) * t472 + t486 * mrSges(7,1) - mrSges(7,2) * t487 + t516 * t505 - t506 * t517;
t465 = m(6) * t474 + mrSges(6,1) * t518 - mrSges(6,3) * t504 - t513 * t531 + t521 * t542 + t471;
t453 = t594 * t459 + t599 * t465;
t533 = mrSges(5,1) * t557 - mrSges(5,3) * t545;
t623 = t599 * t459 - t465 * t594;
t450 = m(5) * t483 - mrSges(5,2) * t547 + mrSges(5,3) * t519 + t527 * t544 - t533 * t557 + t623;
t452 = m(5) * t484 - mrSges(5,1) * t519 + mrSges(5,2) * t520 - t532 * t544 + t533 * t545 + t453;
t439 = t450 * t631 - t452 * t586 + t590 * t640;
t552 = -mrSges(4,1) * t563 + mrSges(4,2) * t564;
t561 = -mrSges(4,2) * t575 + mrSges(4,3) * t563;
t435 = m(4) * t511 + mrSges(4,1) * t572 - mrSges(4,3) * t554 - t552 * t564 + t561 * t575 + t439;
t438 = t450 * t636 + t590 * t452 + t586 * t640;
t562 = mrSges(4,1) * t575 - mrSges(4,3) * t564;
t437 = m(4) * t526 - mrSges(4,1) * t553 + mrSges(4,2) * t554 - t561 * t563 + t562 * t564 + t438;
t444 = t600 * t450 - t456 * t595;
t443 = m(4) * t512 - mrSges(4,2) * t572 + mrSges(4,3) * t553 + t552 * t563 - t562 * t575 + t444;
t426 = t435 * t634 + t591 * t437 + t443 * t635;
t429 = -t435 * t596 + t601 * t443;
t427 = t435 * t629 - t437 * t587 + t443 * t630;
t621 = -mrSges(3,1) * t589 + mrSges(3,2) * t585;
t614 = mrSges(3,1) * t592 - mrSges(3,3) * t637;
t613 = -mrSges(3,2) * t592 + mrSges(3,3) * t633;
t488 = Ifges(7,5) * t517 + Ifges(7,6) * t516 + Ifges(7,3) * t529;
t490 = Ifges(7,1) * t517 + Ifges(7,4) * t516 + Ifges(7,5) * t529;
t462 = -mrSges(7,1) * t472 + mrSges(7,3) * t470 + Ifges(7,4) * t487 + Ifges(7,2) * t486 + Ifges(7,6) * t502 - t488 * t517 + t490 * t529;
t489 = Ifges(7,4) * t517 + Ifges(7,2) * t516 + Ifges(7,6) * t529;
t463 = mrSges(7,2) * t472 - mrSges(7,3) * t469 + Ifges(7,1) * t487 + Ifges(7,4) * t486 + Ifges(7,5) * t502 + t488 * t516 - t489 * t529;
t507 = Ifges(6,5) * t531 + Ifges(6,6) * t530 + Ifges(6,3) * t542;
t508 = Ifges(6,4) * t531 + Ifges(6,2) * t530 + Ifges(6,6) * t542;
t445 = mrSges(6,2) * t478 - mrSges(6,3) * t474 + Ifges(6,1) * t504 + Ifges(6,4) * t503 + Ifges(6,5) * t518 - pkin(13) * t460 - t462 * t593 + t463 * t598 + t507 * t530 - t508 * t542;
t509 = Ifges(6,1) * t531 + Ifges(6,4) * t530 + Ifges(6,5) * t542;
t605 = mrSges(7,1) * t469 - mrSges(7,2) * t470 + Ifges(7,5) * t487 + Ifges(7,6) * t486 + Ifges(7,3) * t502 + t489 * t517 - t490 * t516;
t446 = -mrSges(6,1) * t478 + mrSges(6,3) * t475 + Ifges(6,4) * t504 + Ifges(6,2) * t503 + Ifges(6,6) * t518 - pkin(5) * t460 - t507 * t531 + t509 * t542 - t605;
t524 = Ifges(5,4) * t545 + Ifges(5,2) * t544 + Ifges(5,6) * t557;
t525 = Ifges(5,1) * t545 + Ifges(5,4) * t544 + Ifges(5,5) * t557;
t430 = mrSges(5,1) * t482 - mrSges(5,2) * t483 + Ifges(5,5) * t520 + Ifges(5,6) * t519 + Ifges(5,3) * t547 + pkin(4) * t606 + pkin(12) * t623 + t594 * t445 + t599 * t446 + t545 * t524 - t544 * t525;
t548 = Ifges(4,5) * t564 + Ifges(4,6) * t563 + Ifges(4,3) * t575;
t550 = Ifges(4,1) * t564 + Ifges(4,4) * t563 + Ifges(4,5) * t575;
t523 = Ifges(5,5) * t545 + Ifges(5,6) * t544 + Ifges(5,3) * t557;
t431 = mrSges(5,2) * t484 - mrSges(5,3) * t482 + Ifges(5,1) * t520 + Ifges(5,4) * t519 + Ifges(5,5) * t547 - pkin(12) * t453 + t445 * t599 - t446 * t594 + t523 * t544 - t524 * t557;
t604 = mrSges(6,1) * t474 - mrSges(6,2) * t475 + Ifges(6,5) * t504 + Ifges(6,6) * t503 + Ifges(6,3) * t518 + pkin(5) * t471 + pkin(13) * t461 + t598 * t462 + t593 * t463 + t531 * t508 - t530 * t509;
t432 = -mrSges(5,1) * t484 + mrSges(5,3) * t483 + Ifges(5,4) * t520 + Ifges(5,2) * t519 + Ifges(5,6) * t547 - pkin(4) * t453 - t545 * t523 + t557 * t525 - t604;
t609 = pkin(11) * t444 + t431 * t595 + t432 * t600;
t422 = -mrSges(4,1) * t526 + mrSges(4,3) * t512 + Ifges(4,4) * t554 + Ifges(4,2) * t553 + Ifges(4,6) * t572 - pkin(3) * t438 - t430 * t586 - t548 * t564 + t550 * t575 + t590 * t609;
t549 = Ifges(4,4) * t564 + Ifges(4,2) * t563 + Ifges(4,6) * t575;
t423 = mrSges(4,2) * t526 - mrSges(4,3) * t511 + Ifges(4,1) * t554 + Ifges(4,4) * t553 + Ifges(4,5) * t572 + t431 * t600 - t432 * t595 + t548 * t563 - t549 * t575 + (-t438 * t586 - t439 * t590) * pkin(11);
t610 = pkin(10) * t429 + t422 * t601 + t423 * t596;
t581 = t613 * qJD(1);
t580 = t614 * qJD(1);
t576 = t621 * t628;
t565 = -t577 * t588 + t624;
t559 = -g(3) * t637 + t627;
t558 = -t578 * t585 + t616;
t428 = m(3) * t559 + t613 * qJDD(1) + (t576 * t633 - t580 * t592) * qJD(1) + t429;
t425 = m(3) * t565 + (t621 * qJDD(1) + (t580 * t585 - t581 * t589) * qJD(1)) * t588 + t426;
t424 = m(3) * t558 + t614 * qJDD(1) + (-t576 * t637 + t581 * t592) * qJD(1) + t427;
t421 = mrSges(4,1) * t511 - mrSges(4,2) * t512 + Ifges(4,5) * t554 + Ifges(4,6) * t553 + Ifges(4,3) * t572 + pkin(3) * t439 + t430 * t590 + t549 * t564 - t550 * t563 + t586 * t609;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t625 - mrSges(2,2) * t622 + (mrSges(3,1) * t558 - mrSges(3,2) * t559 + pkin(2) * t427 + t421 * t591 + pkin(1) * (t424 * t589 + t428 * t585) + Ifges(3,3) * qJDD(1) * t592 + t610 * t587) * t592 + (t585 * (mrSges(3,2) * t565 - mrSges(3,3) * t558 - t422 * t596 + t601 * t423 - t426 * t645) + t589 * (-mrSges(3,1) * t565 + mrSges(3,3) * t559 - pkin(2) * t426 - t421 * t587) - pkin(1) * t425 + qJ(2) * (-t424 * t585 + t428 * t589) + (-t427 * t646 + t589 * t610) * t591 + ((Ifges(3,2) * t589 ^ 2 + (Ifges(3,1) * t585 + Ifges(3,4) * t647) * t585) * t588 + 0.2e1 * t592 * (Ifges(3,5) * t585 + Ifges(3,6) * t589)) * qJDD(1)) * t588; t425; t421; t430; t604; t605;];
tauJ  = t1;
