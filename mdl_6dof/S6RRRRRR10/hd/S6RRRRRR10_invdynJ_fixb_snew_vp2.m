% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_invdynJ_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 18:20:55
% EndTime: 2019-05-08 18:22:37
% DurationCPUTime: 99.70s
% Computational Cost: add. (1625731->395), mult. (4168229->540), div. (0->0), fcn. (3612559->18), ass. (0->184)
t585 = cos(pkin(6));
t578 = qJD(1) * t585 + qJD(2);
t581 = sin(pkin(7));
t584 = cos(pkin(7));
t582 = sin(pkin(6));
t596 = cos(qJ(2));
t616 = qJD(1) * t596;
t613 = t582 * t616;
t565 = (t578 * t581 + t584 * t613) * pkin(11);
t590 = sin(qJ(2));
t618 = qJD(1) * t582;
t636 = pkin(11) * t581;
t569 = (-pkin(2) * t596 - t590 * t636) * t618;
t615 = qJD(1) * qJD(2);
t575 = (qJDD(1) * t590 + t596 * t615) * t582;
t577 = qJDD(1) * t585 + qJDD(2);
t598 = qJD(1) ^ 2;
t591 = sin(qJ(1));
t597 = cos(qJ(1));
t609 = -g(1) * t597 - g(2) * t591;
t637 = pkin(10) * t582;
t573 = -pkin(1) * t598 + qJDD(1) * t637 + t609;
t612 = t591 * g(1) - g(2) * t597;
t572 = qJDD(1) * pkin(1) + t598 * t637 + t612;
t630 = t572 * t585;
t610 = -t590 * t573 + t596 * t630;
t617 = qJD(1) * t590;
t635 = pkin(11) * t584;
t530 = -t575 * t635 + t577 * pkin(2) + t578 * t565 + (-g(3) * t596 - t569 * t617) * t582 + t610;
t614 = t582 * t617;
t568 = pkin(2) * t578 - t614 * t635;
t576 = (qJDD(1) * t596 - t590 * t615) * t582;
t605 = t576 * t584 + t577 * t581;
t619 = t596 * t573 + t590 * t630;
t531 = -t578 * t568 + (-g(3) * t590 + t569 * t616) * t582 + t605 * pkin(11) + t619;
t632 = t585 * g(3);
t536 = -t575 * t636 - t576 * pkin(2) - t632 + (-t572 + (-t565 * t596 + t568 * t590) * qJD(1)) * t582;
t589 = sin(qJ(3));
t595 = cos(qJ(3));
t622 = t584 * t595;
t627 = t581 * t595;
t504 = t530 * t622 - t531 * t589 + t536 * t627;
t621 = t584 * t596;
t556 = t578 * t627 + (-t589 * t590 + t595 * t621) * t618;
t545 = t556 * qJD(3) + t595 * t575 + t589 * t605;
t628 = t581 * t589;
t557 = t578 * t628 + (t589 * t621 + t590 * t595) * t618;
t580 = sin(pkin(8));
t634 = pkin(12) * t580;
t546 = -pkin(3) * t556 - t557 * t634;
t566 = t578 * t584 - t581 * t613 + qJD(3);
t583 = cos(pkin(8));
t606 = t556 * t583 + t566 * t580;
t549 = t606 * pkin(12);
t558 = -t576 * t581 + t577 * t584 + qJDD(3);
t633 = pkin(12) * t583;
t488 = pkin(3) * t558 - t545 * t633 - t546 * t557 + t549 * t566 + t504;
t623 = t584 * t589;
t505 = t530 * t623 + t595 * t531 + t536 * t628;
t551 = pkin(3) * t566 - t557 * t633;
t544 = -t557 * qJD(3) - t589 * t575 + t595 * t605;
t607 = t544 * t583 + t558 * t580;
t489 = pkin(12) * t607 + t556 * t546 - t566 * t551 + t505;
t519 = -t530 * t581 + t584 * t536;
t496 = -pkin(3) * t544 - t545 * t634 - t549 * t556 + t551 * t557 + t519;
t588 = sin(qJ(4));
t594 = cos(qJ(4));
t475 = -t588 * t489 + (t488 * t583 + t496 * t580) * t594;
t540 = t594 * t557 + t588 * t606;
t509 = -t540 * qJD(4) - t588 * t545 + t594 * t607;
t539 = -t588 * t557 + t594 * t606;
t510 = t539 * qJD(4) + t594 * t545 + t588 * t607;
t520 = -mrSges(5,1) * t539 + mrSges(5,2) * t540;
t550 = -t556 * t580 + t566 * t583 + qJD(4);
t525 = -mrSges(5,2) * t550 + mrSges(5,3) * t539;
t532 = -t544 * t580 + t558 * t583 + qJDD(4);
t624 = t583 * t588;
t629 = t580 * t588;
t476 = t488 * t624 + t594 * t489 + t496 * t629;
t521 = -pkin(4) * t539 - pkin(13) * t540;
t548 = t550 ^ 2;
t472 = -pkin(4) * t548 + pkin(13) * t532 + t521 * t539 + t476;
t477 = -t580 * t488 + t583 * t496;
t474 = (-t539 * t550 - t510) * pkin(13) + (t540 * t550 - t509) * pkin(4) + t477;
t587 = sin(qJ(5));
t593 = cos(qJ(5));
t468 = t593 * t472 + t587 * t474;
t523 = -t540 * t587 + t550 * t593;
t524 = t540 * t593 + t550 * t587;
t507 = -pkin(5) * t523 - pkin(14) * t524;
t508 = qJDD(5) - t509;
t538 = qJD(5) - t539;
t537 = t538 ^ 2;
t466 = -pkin(5) * t537 + pkin(14) * t508 + t507 * t523 + t468;
t471 = -t532 * pkin(4) - t548 * pkin(13) + t540 * t521 - t475;
t492 = -qJD(5) * t524 - t510 * t587 + t532 * t593;
t493 = qJD(5) * t523 + t510 * t593 + t532 * t587;
t469 = (-t523 * t538 - t493) * pkin(14) + (t524 * t538 - t492) * pkin(5) + t471;
t586 = sin(qJ(6));
t592 = cos(qJ(6));
t462 = -t466 * t586 + t469 * t592;
t512 = -t524 * t586 + t538 * t592;
t480 = qJD(6) * t512 + t493 * t592 + t508 * t586;
t491 = qJDD(6) - t492;
t513 = t524 * t592 + t538 * t586;
t497 = -mrSges(7,1) * t512 + mrSges(7,2) * t513;
t522 = qJD(6) - t523;
t498 = -mrSges(7,2) * t522 + mrSges(7,3) * t512;
t460 = m(7) * t462 + mrSges(7,1) * t491 - mrSges(7,3) * t480 - t497 * t513 + t498 * t522;
t463 = t466 * t592 + t469 * t586;
t479 = -qJD(6) * t513 - t493 * t586 + t508 * t592;
t499 = mrSges(7,1) * t522 - mrSges(7,3) * t513;
t461 = m(7) * t463 - mrSges(7,2) * t491 + mrSges(7,3) * t479 + t497 * t512 - t499 * t522;
t453 = t460 * t592 + t461 * t586;
t514 = -mrSges(6,2) * t538 + mrSges(6,3) * t523;
t515 = mrSges(6,1) * t538 - mrSges(6,3) * t524;
t601 = -m(6) * t471 + t492 * mrSges(6,1) - mrSges(6,2) * t493 + t523 * t514 - t515 * t524 - t453;
t449 = m(5) * t475 + mrSges(5,1) * t532 - mrSges(5,3) * t510 - t520 * t540 + t525 * t550 + t601;
t631 = t449 * t594;
t626 = t582 * t590;
t625 = t582 * t596;
t454 = -t460 * t586 + t592 * t461;
t506 = -mrSges(6,1) * t523 + mrSges(6,2) * t524;
t452 = m(6) * t468 - mrSges(6,2) * t508 + mrSges(6,3) * t492 + t506 * t523 - t515 * t538 + t454;
t467 = -t472 * t587 + t474 * t593;
t465 = -pkin(5) * t508 - pkin(14) * t537 + t507 * t524 - t467;
t464 = -m(7) * t465 + t479 * mrSges(7,1) - mrSges(7,2) * t480 + t512 * t498 - t499 * t513;
t458 = m(6) * t467 + mrSges(6,1) * t508 - mrSges(6,3) * t493 - t506 * t524 + t514 * t538 + t464;
t446 = t587 * t452 + t593 * t458;
t526 = mrSges(5,1) * t550 - mrSges(5,3) * t540;
t611 = t593 * t452 - t458 * t587;
t443 = m(5) * t476 - mrSges(5,2) * t532 + mrSges(5,3) * t509 + t520 * t539 - t526 * t550 + t611;
t445 = m(5) * t477 - mrSges(5,1) * t509 + mrSges(5,2) * t510 - t525 * t539 + t526 * t540 + t446;
t432 = t443 * t624 - t445 * t580 + t583 * t631;
t547 = -mrSges(4,1) * t556 + mrSges(4,2) * t557;
t552 = -mrSges(4,2) * t566 + mrSges(4,3) * t556;
t428 = m(4) * t504 + mrSges(4,1) * t558 - mrSges(4,3) * t545 - t547 * t557 + t552 * t566 + t432;
t431 = t443 * t629 + t583 * t445 + t580 * t631;
t553 = mrSges(4,1) * t566 - mrSges(4,3) * t557;
t430 = m(4) * t519 - mrSges(4,1) * t544 + mrSges(4,2) * t545 - t552 * t556 + t553 * t557 + t431;
t437 = t594 * t443 - t449 * t588;
t436 = m(4) * t505 - mrSges(4,2) * t558 + mrSges(4,3) * t544 + t547 * t556 - t553 * t566 + t437;
t419 = t428 * t627 + t584 * t430 + t436 * t628;
t422 = -t428 * t589 + t595 * t436;
t420 = t428 * t622 - t430 * t581 + t436 * t623;
t481 = Ifges(7,5) * t513 + Ifges(7,6) * t512 + Ifges(7,3) * t522;
t483 = Ifges(7,1) * t513 + Ifges(7,4) * t512 + Ifges(7,5) * t522;
t455 = -mrSges(7,1) * t465 + mrSges(7,3) * t463 + Ifges(7,4) * t480 + Ifges(7,2) * t479 + Ifges(7,6) * t491 - t481 * t513 + t483 * t522;
t482 = Ifges(7,4) * t513 + Ifges(7,2) * t512 + Ifges(7,6) * t522;
t456 = mrSges(7,2) * t465 - mrSges(7,3) * t462 + Ifges(7,1) * t480 + Ifges(7,4) * t479 + Ifges(7,5) * t491 + t481 * t512 - t482 * t522;
t500 = Ifges(6,5) * t524 + Ifges(6,6) * t523 + Ifges(6,3) * t538;
t501 = Ifges(6,4) * t524 + Ifges(6,2) * t523 + Ifges(6,6) * t538;
t438 = mrSges(6,2) * t471 - mrSges(6,3) * t467 + Ifges(6,1) * t493 + Ifges(6,4) * t492 + Ifges(6,5) * t508 - pkin(14) * t453 - t455 * t586 + t456 * t592 + t500 * t523 - t501 * t538;
t502 = Ifges(6,1) * t524 + Ifges(6,4) * t523 + Ifges(6,5) * t538;
t600 = mrSges(7,1) * t462 - mrSges(7,2) * t463 + Ifges(7,5) * t480 + Ifges(7,6) * t479 + Ifges(7,3) * t491 + t482 * t513 - t483 * t512;
t439 = -mrSges(6,1) * t471 + mrSges(6,3) * t468 + Ifges(6,4) * t493 + Ifges(6,2) * t492 + Ifges(6,6) * t508 - pkin(5) * t453 - t500 * t524 + t502 * t538 - t600;
t517 = Ifges(5,4) * t540 + Ifges(5,2) * t539 + Ifges(5,6) * t550;
t518 = Ifges(5,1) * t540 + Ifges(5,4) * t539 + Ifges(5,5) * t550;
t423 = mrSges(5,1) * t475 - mrSges(5,2) * t476 + Ifges(5,5) * t510 + Ifges(5,6) * t509 + Ifges(5,3) * t532 + pkin(4) * t601 + pkin(13) * t611 + t587 * t438 + t593 * t439 + t540 * t517 - t539 * t518;
t541 = Ifges(4,5) * t557 + Ifges(4,6) * t556 + Ifges(4,3) * t566;
t543 = Ifges(4,1) * t557 + Ifges(4,4) * t556 + Ifges(4,5) * t566;
t516 = Ifges(5,5) * t540 + Ifges(5,6) * t539 + Ifges(5,3) * t550;
t424 = mrSges(5,2) * t477 - mrSges(5,3) * t475 + Ifges(5,1) * t510 + Ifges(5,4) * t509 + Ifges(5,5) * t532 - pkin(13) * t446 + t438 * t593 - t439 * t587 + t516 * t539 - t517 * t550;
t599 = mrSges(6,1) * t467 - mrSges(6,2) * t468 + Ifges(6,5) * t493 + Ifges(6,6) * t492 + Ifges(6,3) * t508 + pkin(5) * t464 + pkin(14) * t454 + t592 * t455 + t586 * t456 + t524 * t501 - t523 * t502;
t425 = -mrSges(5,1) * t477 + mrSges(5,3) * t476 + Ifges(5,4) * t510 + Ifges(5,2) * t509 + Ifges(5,6) * t532 - pkin(4) * t446 - t540 * t516 + t550 * t518 - t599;
t602 = pkin(12) * t437 + t424 * t588 + t425 * t594;
t416 = -mrSges(4,1) * t519 + mrSges(4,3) * t505 + Ifges(4,4) * t545 + Ifges(4,2) * t544 + Ifges(4,6) * t558 - pkin(3) * t431 - t580 * t423 - t557 * t541 + t566 * t543 + t583 * t602;
t542 = Ifges(4,4) * t557 + Ifges(4,2) * t556 + Ifges(4,6) * t566;
t417 = mrSges(4,2) * t519 - mrSges(4,3) * t504 + Ifges(4,1) * t545 + Ifges(4,4) * t544 + Ifges(4,5) * t558 + t594 * t424 - t588 * t425 + t556 * t541 - t566 * t542 + (-t431 * t580 - t432 * t583) * pkin(12);
t603 = pkin(11) * t422 + t416 * t595 + t417 * t589;
t574 = (-mrSges(3,1) * t596 + mrSges(3,2) * t590) * t618;
t571 = -mrSges(3,2) * t578 + mrSges(3,3) * t613;
t570 = mrSges(3,1) * t578 - mrSges(3,3) * t614;
t562 = -t582 * t572 - t632;
t561 = Ifges(3,5) * t578 + (Ifges(3,1) * t590 + Ifges(3,4) * t596) * t618;
t560 = Ifges(3,6) * t578 + (Ifges(3,4) * t590 + Ifges(3,2) * t596) * t618;
t559 = Ifges(3,3) * t578 + (Ifges(3,5) * t590 + Ifges(3,6) * t596) * t618;
t555 = -g(3) * t626 + t619;
t554 = -g(3) * t625 + t610;
t421 = m(3) * t555 - mrSges(3,2) * t577 + mrSges(3,3) * t576 - t570 * t578 + t574 * t613 + t422;
t418 = m(3) * t554 + mrSges(3,1) * t577 - mrSges(3,3) * t575 + t571 * t578 - t574 * t614 + t420;
t415 = mrSges(4,1) * t504 - mrSges(4,2) * t505 + Ifges(4,5) * t545 + Ifges(4,6) * t544 + Ifges(4,3) * t558 + pkin(3) * t432 + t583 * t423 + t557 * t542 - t556 * t543 + t580 * t602;
t414 = mrSges(3,1) * t554 - mrSges(3,2) * t555 + Ifges(3,5) * t575 + Ifges(3,6) * t576 + Ifges(3,3) * t577 + pkin(2) * t420 + t584 * t415 + (t560 * t590 - t561 * t596) * t618 + t603 * t581;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t612 - mrSges(2,2) * t609 + (t559 * t613 + mrSges(3,2) * t562 - mrSges(3,3) * t554 + Ifges(3,1) * t575 + Ifges(3,4) * t576 + Ifges(3,5) * t577 - t589 * t416 + t595 * t417 - t578 * t560 + (-t419 * t581 - t420 * t584) * pkin(11)) * t626 + (-mrSges(3,1) * t562 + mrSges(3,3) * t555 + Ifges(3,4) * t575 + Ifges(3,2) * t576 + Ifges(3,6) * t577 - pkin(2) * t419 - t581 * t415 - t559 * t614 + t578 * t561 + t584 * t603) * t625 + t585 * t414 + pkin(1) * ((t418 * t596 + t421 * t590) * t585 + (-m(3) * t562 + t576 * mrSges(3,1) - t575 * mrSges(3,2) + (-t570 * t590 + t571 * t596) * t618 - t419) * t582) + (-t418 * t590 + t421 * t596) * t637; t414; t415; t423; t599; t600;];
tauJ  = t1;
