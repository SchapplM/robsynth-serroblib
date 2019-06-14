% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-09 12:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRR14_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_invdynJ_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-09 11:27:42
% EndTime: 2019-05-09 11:30:14
% DurationCPUTime: 88.98s
% Computational Cost: add. (1356389->394), mult. (3678748->542), div. (0->0), fcn. (3171749->18), ass. (0->184)
t594 = cos(pkin(6));
t585 = qJD(1) * t594 + qJD(2);
t589 = sin(pkin(7));
t593 = cos(pkin(7));
t590 = sin(pkin(6));
t603 = cos(qJ(2));
t623 = qJD(1) * t603;
t620 = t590 * t623;
t570 = (t585 * t589 + t593 * t620) * qJ(3);
t598 = sin(qJ(2));
t625 = qJD(1) * t590;
t640 = qJ(3) * t589;
t576 = (-pkin(2) * t603 - t598 * t640) * t625;
t622 = qJD(1) * qJD(2);
t582 = (qJDD(1) * t598 + t603 * t622) * t590;
t584 = qJDD(1) * t594 + qJDD(2);
t605 = qJD(1) ^ 2;
t599 = sin(qJ(1));
t604 = cos(qJ(1));
t616 = -g(1) * t604 - g(2) * t599;
t644 = pkin(10) * t590;
t580 = -pkin(1) * t605 + qJDD(1) * t644 + t616;
t619 = t599 * g(1) - g(2) * t604;
t579 = qJDD(1) * pkin(1) + t605 * t644 + t619;
t637 = t579 * t594;
t617 = -t598 * t580 + t603 * t637;
t624 = qJD(1) * t598;
t639 = qJ(3) * t593;
t535 = -t582 * t639 + t584 * pkin(2) + t585 * t570 + (-g(3) * t603 - t576 * t624) * t590 + t617;
t621 = t590 * t624;
t575 = pkin(2) * t585 - t621 * t639;
t583 = (qJDD(1) * t603 - t598 * t622) * t590;
t612 = t583 * t593 + t584 * t589;
t626 = t603 * t580 + t598 * t637;
t536 = -t585 * t575 + (-g(3) * t598 + t576 * t623) * t590 + t612 * qJ(3) + t626;
t641 = t594 * g(3);
t540 = -t582 * t640 - t583 * pkin(2) - t641 + (-t579 + (-t570 * t603 + t575 * t598) * qJD(1)) * t590;
t587 = sin(pkin(14));
t591 = cos(pkin(14));
t628 = t593 * t603;
t636 = t587 * t589;
t564 = t585 * t636 + (t587 * t628 + t591 * t598) * t625;
t630 = t591 * t593;
t633 = t589 * t591;
t509 = -0.2e1 * qJD(3) * t564 + t535 * t630 - t536 * t587 + t540 * t633;
t563 = t585 * t633 + (-t587 * t598 + t591 * t628) * t625;
t588 = sin(pkin(8));
t643 = pkin(11) * t588;
t549 = -pkin(3) * t563 - t564 * t643;
t573 = t585 * t593 - t589 * t620;
t592 = cos(pkin(8));
t613 = t563 * t592 + t573 * t588;
t552 = t613 * pkin(11);
t558 = t591 * t582 + t587 * t612;
t565 = -t583 * t589 + t584 * t593;
t642 = pkin(11) * t592;
t493 = pkin(3) * t565 - t549 * t564 + t552 * t573 - t558 * t642 + t509;
t635 = t587 * t593;
t510 = 0.2e1 * qJD(3) * t563 + t535 * t635 + t591 * t536 + t540 * t636;
t554 = pkin(3) * t573 - t564 * t642;
t557 = -t587 * t582 + t591 * t612;
t614 = t557 * t592 + t565 * t588;
t494 = pkin(11) * t614 + t563 * t549 - t573 * t554 + t510;
t524 = -t535 * t589 + t593 * t540 + qJDD(3);
t498 = -pkin(3) * t557 - t552 * t563 + t554 * t564 - t558 * t643 + t524;
t597 = sin(qJ(4));
t602 = cos(qJ(4));
t480 = -t597 * t494 + (t493 * t592 + t498 * t588) * t602;
t544 = t602 * t564 + t597 * t613;
t522 = -t544 * qJD(4) - t597 * t558 + t602 * t614;
t543 = -t597 * t564 + t602 * t613;
t523 = t543 * qJD(4) + t602 * t558 + t597 * t614;
t525 = -mrSges(5,1) * t543 + mrSges(5,2) * t544;
t553 = -t563 * t588 + t573 * t592 + qJD(4);
t530 = -mrSges(5,2) * t553 + mrSges(5,3) * t543;
t548 = -t557 * t588 + t565 * t592 + qJDD(4);
t629 = t592 * t597;
t634 = t588 * t597;
t481 = t493 * t629 + t602 * t494 + t498 * t634;
t526 = -pkin(4) * t543 - pkin(12) * t544;
t551 = t553 ^ 2;
t477 = -pkin(4) * t551 + pkin(12) * t548 + t526 * t543 + t481;
t482 = -t588 * t493 + t592 * t498;
t479 = (-t543 * t553 - t523) * pkin(12) + (t544 * t553 - t522) * pkin(4) + t482;
t596 = sin(qJ(5));
t601 = cos(qJ(5));
t473 = t601 * t477 + t596 * t479;
t528 = -t544 * t596 + t553 * t601;
t529 = t544 * t601 + t553 * t596;
t512 = -pkin(5) * t528 - pkin(13) * t529;
t521 = qJDD(5) - t522;
t542 = qJD(5) - t543;
t541 = t542 ^ 2;
t471 = -pkin(5) * t541 + pkin(13) * t521 + t512 * t528 + t473;
t476 = -t548 * pkin(4) - t551 * pkin(12) + t544 * t526 - t480;
t501 = -qJD(5) * t529 - t523 * t596 + t548 * t601;
t502 = qJD(5) * t528 + t523 * t601 + t548 * t596;
t474 = (-t528 * t542 - t502) * pkin(13) + (t529 * t542 - t501) * pkin(5) + t476;
t595 = sin(qJ(6));
t600 = cos(qJ(6));
t467 = -t471 * t595 + t474 * t600;
t514 = -t529 * t595 + t542 * t600;
t485 = qJD(6) * t514 + t502 * t600 + t521 * t595;
t515 = t529 * t600 + t542 * t595;
t495 = -mrSges(7,1) * t514 + mrSges(7,2) * t515;
t500 = qJDD(6) - t501;
t527 = qJD(6) - t528;
t503 = -mrSges(7,2) * t527 + mrSges(7,3) * t514;
t465 = m(7) * t467 + mrSges(7,1) * t500 - mrSges(7,3) * t485 - t495 * t515 + t503 * t527;
t468 = t471 * t600 + t474 * t595;
t484 = -qJD(6) * t515 - t502 * t595 + t521 * t600;
t504 = mrSges(7,1) * t527 - mrSges(7,3) * t515;
t466 = m(7) * t468 - mrSges(7,2) * t500 + mrSges(7,3) * t484 + t495 * t514 - t504 * t527;
t458 = t465 * t600 + t466 * t595;
t516 = -mrSges(6,2) * t542 + mrSges(6,3) * t528;
t517 = mrSges(6,1) * t542 - mrSges(6,3) * t529;
t608 = -m(6) * t476 + t501 * mrSges(6,1) - mrSges(6,2) * t502 + t528 * t516 - t517 * t529 - t458;
t454 = m(5) * t480 + mrSges(5,1) * t548 - mrSges(5,3) * t523 - t525 * t544 + t530 * t553 + t608;
t638 = t454 * t602;
t632 = t590 * t598;
t631 = t590 * t603;
t459 = -t465 * t595 + t600 * t466;
t511 = -mrSges(6,1) * t528 + mrSges(6,2) * t529;
t457 = m(6) * t473 - mrSges(6,2) * t521 + mrSges(6,3) * t501 + t511 * t528 - t517 * t542 + t459;
t472 = -t477 * t596 + t479 * t601;
t470 = -pkin(5) * t521 - pkin(13) * t541 + t512 * t529 - t472;
t469 = -m(7) * t470 + t484 * mrSges(7,1) - mrSges(7,2) * t485 + t514 * t503 - t504 * t515;
t463 = m(6) * t472 + mrSges(6,1) * t521 - mrSges(6,3) * t502 - t511 * t529 + t516 * t542 + t469;
t451 = t596 * t457 + t601 * t463;
t531 = mrSges(5,1) * t553 - mrSges(5,3) * t544;
t618 = t601 * t457 - t463 * t596;
t448 = m(5) * t481 - mrSges(5,2) * t548 + mrSges(5,3) * t522 + t525 * t543 - t531 * t553 + t618;
t450 = m(5) * t482 - mrSges(5,1) * t522 + mrSges(5,2) * t523 - t530 * t543 + t531 * t544 + t451;
t437 = t448 * t629 - t450 * t588 + t592 * t638;
t550 = -mrSges(4,1) * t563 + mrSges(4,2) * t564;
t555 = -mrSges(4,2) * t573 + mrSges(4,3) * t563;
t433 = m(4) * t509 + mrSges(4,1) * t565 - mrSges(4,3) * t558 - t550 * t564 + t555 * t573 + t437;
t436 = t448 * t634 + t592 * t450 + t588 * t638;
t556 = mrSges(4,1) * t573 - mrSges(4,3) * t564;
t435 = m(4) * t524 - mrSges(4,1) * t557 + mrSges(4,2) * t558 - t555 * t563 + t556 * t564 + t436;
t442 = t602 * t448 - t454 * t597;
t441 = m(4) * t510 - mrSges(4,2) * t565 + mrSges(4,3) * t557 + t550 * t563 - t556 * t573 + t442;
t424 = t433 * t633 + t593 * t435 + t441 * t636;
t427 = -t433 * t587 + t591 * t441;
t425 = t433 * t630 - t435 * t589 + t441 * t635;
t486 = Ifges(7,5) * t515 + Ifges(7,6) * t514 + Ifges(7,3) * t527;
t488 = Ifges(7,1) * t515 + Ifges(7,4) * t514 + Ifges(7,5) * t527;
t460 = -mrSges(7,1) * t470 + mrSges(7,3) * t468 + Ifges(7,4) * t485 + Ifges(7,2) * t484 + Ifges(7,6) * t500 - t486 * t515 + t488 * t527;
t487 = Ifges(7,4) * t515 + Ifges(7,2) * t514 + Ifges(7,6) * t527;
t461 = mrSges(7,2) * t470 - mrSges(7,3) * t467 + Ifges(7,1) * t485 + Ifges(7,4) * t484 + Ifges(7,5) * t500 + t486 * t514 - t487 * t527;
t505 = Ifges(6,5) * t529 + Ifges(6,6) * t528 + Ifges(6,3) * t542;
t506 = Ifges(6,4) * t529 + Ifges(6,2) * t528 + Ifges(6,6) * t542;
t443 = mrSges(6,2) * t476 - mrSges(6,3) * t472 + Ifges(6,1) * t502 + Ifges(6,4) * t501 + Ifges(6,5) * t521 - pkin(13) * t458 - t460 * t595 + t461 * t600 + t505 * t528 - t506 * t542;
t507 = Ifges(6,1) * t529 + Ifges(6,4) * t528 + Ifges(6,5) * t542;
t607 = mrSges(7,1) * t467 - mrSges(7,2) * t468 + Ifges(7,5) * t485 + Ifges(7,6) * t484 + Ifges(7,3) * t500 + t487 * t515 - t488 * t514;
t444 = -mrSges(6,1) * t476 + mrSges(6,3) * t473 + Ifges(6,4) * t502 + Ifges(6,2) * t501 + Ifges(6,6) * t521 - pkin(5) * t458 - t505 * t529 + t507 * t542 - t607;
t518 = Ifges(5,5) * t544 + Ifges(5,6) * t543 + Ifges(5,3) * t553;
t519 = Ifges(5,4) * t544 + Ifges(5,2) * t543 + Ifges(5,6) * t553;
t429 = mrSges(5,2) * t482 - mrSges(5,3) * t480 + Ifges(5,1) * t523 + Ifges(5,4) * t522 + Ifges(5,5) * t548 - pkin(12) * t451 + t443 * t601 - t444 * t596 + t518 * t543 - t519 * t553;
t520 = Ifges(5,1) * t544 + Ifges(5,4) * t543 + Ifges(5,5) * t553;
t606 = mrSges(6,1) * t472 - mrSges(6,2) * t473 + Ifges(6,5) * t502 + Ifges(6,6) * t501 + Ifges(6,3) * t521 + pkin(5) * t469 + pkin(13) * t459 + t600 * t460 + t595 * t461 + t529 * t506 - t528 * t507;
t430 = -mrSges(5,1) * t482 + mrSges(5,3) * t481 + Ifges(5,4) * t523 + Ifges(5,2) * t522 + Ifges(5,6) * t548 - pkin(4) * t451 - t544 * t518 + t553 * t520 - t606;
t610 = pkin(11) * t442 + t429 * t597 + t430 * t602;
t428 = mrSges(5,1) * t480 - mrSges(5,2) * t481 + Ifges(5,5) * t523 + Ifges(5,6) * t522 + Ifges(5,3) * t548 + pkin(4) * t608 + pkin(12) * t618 + t596 * t443 + t601 * t444 + t544 * t519 - t543 * t520;
t545 = Ifges(4,5) * t564 + Ifges(4,6) * t563 + Ifges(4,3) * t573;
t547 = Ifges(4,1) * t564 + Ifges(4,4) * t563 + Ifges(4,5) * t573;
t421 = -mrSges(4,1) * t524 + mrSges(4,3) * t510 + Ifges(4,4) * t558 + Ifges(4,2) * t557 + Ifges(4,6) * t565 - pkin(3) * t436 - t588 * t428 - t564 * t545 + t573 * t547 + t592 * t610;
t546 = Ifges(4,4) * t564 + Ifges(4,2) * t563 + Ifges(4,6) * t573;
t422 = mrSges(4,2) * t524 - mrSges(4,3) * t509 + Ifges(4,1) * t558 + Ifges(4,4) * t557 + Ifges(4,5) * t565 + t602 * t429 - t597 * t430 + t563 * t545 - t573 * t546 + (-t436 * t588 - t437 * t592) * pkin(11);
t609 = qJ(3) * t427 + t421 * t591 + t422 * t587;
t581 = (-mrSges(3,1) * t603 + mrSges(3,2) * t598) * t625;
t578 = -mrSges(3,2) * t585 + mrSges(3,3) * t620;
t577 = mrSges(3,1) * t585 - mrSges(3,3) * t621;
t569 = -t590 * t579 - t641;
t568 = Ifges(3,5) * t585 + (Ifges(3,1) * t598 + Ifges(3,4) * t603) * t625;
t567 = Ifges(3,6) * t585 + (Ifges(3,4) * t598 + Ifges(3,2) * t603) * t625;
t566 = Ifges(3,3) * t585 + (Ifges(3,5) * t598 + Ifges(3,6) * t603) * t625;
t560 = -g(3) * t632 + t626;
t559 = -g(3) * t631 + t617;
t426 = m(3) * t560 - mrSges(3,2) * t584 + mrSges(3,3) * t583 - t577 * t585 + t581 * t620 + t427;
t423 = m(3) * t559 + mrSges(3,1) * t584 - mrSges(3,3) * t582 + t578 * t585 - t581 * t621 + t425;
t420 = mrSges(4,1) * t509 - mrSges(4,2) * t510 + Ifges(4,5) * t558 + Ifges(4,6) * t557 + Ifges(4,3) * t565 + pkin(3) * t437 + t592 * t428 + t564 * t546 - t563 * t547 + t588 * t610;
t419 = mrSges(3,1) * t559 - mrSges(3,2) * t560 + Ifges(3,5) * t582 + Ifges(3,6) * t583 + Ifges(3,3) * t584 + pkin(2) * t425 + t593 * t420 + (t567 * t598 - t568 * t603) * t625 + t609 * t589;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t619 - mrSges(2,2) * t616 + (t566 * t620 + mrSges(3,2) * t569 - mrSges(3,3) * t559 + Ifges(3,1) * t582 + Ifges(3,4) * t583 + Ifges(3,5) * t584 - t587 * t421 + t591 * t422 - t585 * t567 + (-t424 * t589 - t425 * t593) * qJ(3)) * t632 + (-mrSges(3,1) * t569 + mrSges(3,3) * t560 + Ifges(3,4) * t582 + Ifges(3,2) * t583 + Ifges(3,6) * t584 - pkin(2) * t424 - t589 * t420 - t566 * t621 + t585 * t568 + t593 * t609) * t631 + t594 * t419 + pkin(1) * ((t423 * t603 + t426 * t598) * t594 + (-m(3) * t569 + t583 * mrSges(3,1) - t582 * mrSges(3,2) + (-t577 * t598 + t578 * t603) * t625 - t424) * t590) + (-t423 * t598 + t426 * t603) * t644; t419; t435; t428; t606; t607;];
tauJ  = t1;
