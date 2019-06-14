% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR15
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-08 03:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR15_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR15_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR15_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 03:22:36
% EndTime: 2019-05-08 03:22:57
% DurationCPUTime: 14.75s
% Computational Cost: add. (225617->362), mult. (557501->465), div. (0->0), fcn. (467247->14), ass. (0->164)
t629 = Ifges(5,1) + Ifges(6,2);
t619 = Ifges(5,4) + Ifges(6,6);
t618 = Ifges(5,5) - Ifges(6,4);
t628 = -Ifges(5,2) - Ifges(6,3);
t617 = Ifges(5,6) - Ifges(6,5);
t627 = Ifges(5,3) + Ifges(6,1);
t569 = cos(pkin(6));
t564 = qJD(1) * t569 + qJD(2);
t566 = sin(pkin(7));
t568 = cos(pkin(7));
t567 = sin(pkin(6));
t577 = cos(qJ(2));
t601 = qJD(1) * t577;
t598 = t567 * t601;
t594 = t568 * t598;
t549 = (t564 * t566 + t594) * pkin(10);
t573 = sin(qJ(2));
t603 = qJD(1) * t567;
t622 = pkin(10) * t566;
t553 = (-pkin(2) * t577 - t573 * t622) * t603;
t600 = qJD(1) * qJD(2);
t559 = (qJDD(1) * t573 + t577 * t600) * t567;
t563 = qJDD(1) * t569 + qJDD(2);
t579 = qJD(1) ^ 2;
t574 = sin(qJ(1));
t578 = cos(qJ(1));
t593 = -g(1) * t578 - g(2) * t574;
t623 = pkin(9) * t567;
t557 = -pkin(1) * t579 + qJDD(1) * t623 + t593;
t597 = t574 * g(1) - g(2) * t578;
t556 = qJDD(1) * pkin(1) + t579 * t623 + t597;
t614 = t556 * t569;
t595 = -t573 * t557 + t577 * t614;
t602 = qJD(1) * t573;
t621 = pkin(10) * t568;
t500 = -t559 * t621 + t563 * pkin(2) + t564 * t549 + (-g(3) * t577 - t553 * t602) * t567 + t595;
t599 = t567 * t602;
t552 = pkin(2) * t564 - t599 * t621;
t560 = (qJDD(1) * t577 - t573 * t600) * t567;
t591 = t560 * t568 + t563 * t566;
t604 = t577 * t557 + t573 * t614;
t501 = -t564 * t552 + (-g(3) * t573 + t553 * t601) * t567 + t591 * pkin(10) + t604;
t620 = t569 * g(3);
t507 = -t559 * t622 - t560 * pkin(2) - t620 + (-t556 + (-t549 * t577 + t552 * t573) * qJD(1)) * t567;
t572 = sin(qJ(3));
t576 = cos(qJ(3));
t468 = -t572 * t501 + (t500 * t568 + t507 * t566) * t576;
t609 = t568 * t572;
t613 = t566 * t572;
t539 = t564 * t613 + (t573 * t576 + t577 * t609) * t603;
t522 = -t539 * qJD(3) - t572 * t559 + t576 * t591;
t550 = t564 * t568 - t566 * t598 + qJD(3);
t571 = sin(qJ(4));
t624 = cos(qJ(4));
t529 = t571 * t539 - t550 * t624;
t612 = t566 * t576;
t538 = t564 * t612 - t572 * t599 + t576 * t594;
t536 = qJD(4) - t538;
t511 = mrSges(6,1) * t529 - mrSges(6,3) * t536;
t521 = qJDD(4) - t522;
t469 = t500 * t609 + t576 * t501 + t507 * t613;
t525 = -pkin(3) * t538 - pkin(11) * t539;
t540 = -t560 * t566 + t563 * t568 + qJDD(3);
t548 = t550 ^ 2;
t460 = -pkin(3) * t548 + pkin(11) * t540 + t525 * t538 + t469;
t475 = -t566 * t500 + t568 * t507;
t523 = t538 * qJD(3) + t576 * t559 + t572 * t591;
t462 = (-t538 * t550 - t523) * pkin(11) + (t539 * t550 - t522) * pkin(3) + t475;
t455 = -t571 * t460 + t462 * t624;
t530 = t539 * t624 + t571 * t550;
t502 = pkin(4) * t529 - qJ(5) * t530;
t535 = t536 ^ 2;
t453 = -t521 * pkin(4) - t535 * qJ(5) + t530 * t502 + qJDD(5) - t455;
t483 = -t529 * qJD(4) + t523 * t624 + t571 * t540;
t615 = t529 * t536;
t448 = (t529 * t530 - t521) * pkin(12) + (t483 + t615) * pkin(5) + t453;
t482 = t530 * qJD(4) + t571 * t523 - t540 * t624;
t515 = pkin(5) * t530 - pkin(12) * t536;
t528 = t529 ^ 2;
t459 = -t540 * pkin(3) - t548 * pkin(11) + t539 * t525 - t468;
t625 = -2 * qJD(5);
t581 = (-t483 + t615) * qJ(5) + t459 + (t536 * pkin(4) + t625) * t530;
t451 = -t528 * pkin(5) - t530 * t515 + (pkin(4) + pkin(12)) * t482 + t581;
t570 = sin(qJ(6));
t575 = cos(qJ(6));
t446 = t448 * t575 - t451 * t570;
t509 = t529 * t575 - t536 * t570;
t467 = qJD(6) * t509 + t482 * t570 + t521 * t575;
t510 = t529 * t570 + t536 * t575;
t476 = -mrSges(7,1) * t509 + mrSges(7,2) * t510;
t481 = qJDD(6) + t483;
t527 = qJD(6) + t530;
t486 = -mrSges(7,2) * t527 + mrSges(7,3) * t509;
t443 = m(7) * t446 + mrSges(7,1) * t481 - mrSges(7,3) * t467 - t476 * t510 + t486 * t527;
t447 = t448 * t570 + t451 * t575;
t466 = -qJD(6) * t510 + t482 * t575 - t521 * t570;
t487 = mrSges(7,1) * t527 - mrSges(7,3) * t510;
t444 = m(7) * t447 - mrSges(7,2) * t481 + mrSges(7,3) * t466 + t476 * t509 - t487 * t527;
t434 = t575 * t443 + t570 * t444;
t504 = -mrSges(6,2) * t529 - mrSges(6,3) * t530;
t586 = -m(6) * t453 - t483 * mrSges(6,1) - t530 * t504 - t434;
t432 = t521 * mrSges(6,2) + t536 * t511 - t586;
t456 = t624 * t460 + t571 * t462;
t585 = -t535 * pkin(4) + t521 * qJ(5) - t529 * t502 + t456;
t450 = -t482 * pkin(5) - t528 * pkin(12) + ((2 * qJD(5)) + t515) * t536 + t585;
t471 = Ifges(7,5) * t510 + Ifges(7,6) * t509 + Ifges(7,3) * t527;
t473 = Ifges(7,1) * t510 + Ifges(7,4) * t509 + Ifges(7,5) * t527;
t435 = -mrSges(7,1) * t450 + mrSges(7,3) * t447 + Ifges(7,4) * t467 + Ifges(7,2) * t466 + Ifges(7,6) * t481 - t471 * t510 + t473 * t527;
t472 = Ifges(7,4) * t510 + Ifges(7,2) * t509 + Ifges(7,6) * t527;
t436 = mrSges(7,2) * t450 - mrSges(7,3) * t446 + Ifges(7,1) * t467 + Ifges(7,4) * t466 + Ifges(7,5) * t481 + t471 * t509 - t472 * t527;
t452 = t536 * t625 - t585;
t512 = mrSges(6,1) * t530 + mrSges(6,2) * t536;
t587 = -m(7) * t450 + t466 * mrSges(7,1) - t467 * mrSges(7,2) + t509 * t486 - t510 * t487;
t583 = -m(6) * t452 + t521 * mrSges(6,3) + t536 * t512 - t587;
t605 = -t619 * t529 + t629 * t530 + t618 * t536;
t606 = t628 * t529 + t619 * t530 + t617 * t536;
t626 = -t482 * t617 + t483 * t618 + t627 * t521 + t529 * t605 + t530 * t606 + mrSges(5,1) * t455 - mrSges(5,2) * t456 + mrSges(6,2) * t453 - mrSges(6,3) * t452 - pkin(4) * t432 - pkin(12) * t434 + qJ(5) * (-t482 * mrSges(6,1) - t529 * t504 + t583) - t570 * t435 + t575 * t436;
t611 = t567 * t573;
t610 = t567 * t577;
t503 = mrSges(5,1) * t529 + mrSges(5,2) * t530;
t513 = -mrSges(5,2) * t536 - mrSges(5,3) * t529;
t431 = m(5) * t455 - t483 * mrSges(5,3) - t530 * t503 + (-t511 + t513) * t536 + (mrSges(5,1) - mrSges(6,2)) * t521 + t586;
t514 = mrSges(5,1) * t536 - mrSges(5,3) * t530;
t439 = m(5) * t456 - t521 * mrSges(5,2) - t536 * t514 + (-t503 - t504) * t529 + (-mrSges(5,3) - mrSges(6,1)) * t482 + t583;
t426 = t624 * t431 + t571 * t439;
t608 = -t570 * t443 + t575 * t444;
t607 = t617 * t529 - t618 * t530 - t627 * t536;
t524 = -mrSges(4,1) * t538 + mrSges(4,2) * t539;
t532 = mrSges(4,1) * t550 - mrSges(4,3) * t539;
t596 = -t431 * t571 + t624 * t439;
t423 = m(4) * t469 - mrSges(4,2) * t540 + mrSges(4,3) * t522 + t524 * t538 - t532 * t550 + t596;
t531 = -mrSges(4,2) * t550 + mrSges(4,3) * t538;
t425 = m(4) * t475 - mrSges(4,1) * t522 + mrSges(4,2) * t523 - t531 * t538 + t532 * t539 + t426;
t454 = t482 * pkin(4) + t581;
t589 = -m(6) * t454 + t482 * mrSges(6,2) + t529 * t511 - t608;
t582 = -m(5) * t459 - t482 * mrSges(5,1) - t529 * t513 + (t512 - t514) * t530 + (-mrSges(5,2) + mrSges(6,3)) * t483 + t589;
t429 = m(4) * t468 + t540 * mrSges(4,1) - t523 * mrSges(4,3) - t539 * t524 + t550 * t531 + t582;
t414 = t423 * t613 + t568 * t425 + t429 * t612;
t419 = t576 * t423 - t429 * t572;
t415 = t568 * t576 * t429 + t423 * t609 - t425 * t566;
t433 = -t483 * mrSges(6,3) - t530 * t512 - t589;
t416 = -mrSges(5,1) * t459 - mrSges(6,1) * t452 + mrSges(6,2) * t454 + mrSges(5,3) * t456 - pkin(4) * t433 - pkin(5) * t587 - pkin(12) * t608 - t575 * t435 - t570 * t436 + t628 * t482 + t619 * t483 + t617 * t521 + t607 * t530 + t605 * t536;
t584 = mrSges(7,1) * t446 - mrSges(7,2) * t447 + Ifges(7,5) * t467 + Ifges(7,6) * t466 + Ifges(7,3) * t481 + t510 * t472 - t509 * t473;
t417 = mrSges(6,1) * t453 + mrSges(5,2) * t459 - mrSges(5,3) * t455 - mrSges(6,3) * t454 + pkin(5) * t434 - qJ(5) * t433 - t619 * t482 + t629 * t483 + t618 * t521 + t607 * t529 - t606 * t536 + t584;
t517 = Ifges(4,5) * t539 + Ifges(4,6) * t538 + Ifges(4,3) * t550;
t518 = Ifges(4,4) * t539 + Ifges(4,2) * t538 + Ifges(4,6) * t550;
t411 = mrSges(4,2) * t475 - mrSges(4,3) * t468 + Ifges(4,1) * t523 + Ifges(4,4) * t522 + Ifges(4,5) * t540 - pkin(11) * t426 - t571 * t416 + t417 * t624 + t538 * t517 - t550 * t518;
t519 = Ifges(4,1) * t539 + Ifges(4,4) * t538 + Ifges(4,5) * t550;
t412 = -mrSges(4,1) * t475 + mrSges(4,3) * t469 + Ifges(4,4) * t523 + Ifges(4,2) * t522 + Ifges(4,6) * t540 - pkin(3) * t426 - t539 * t517 + t550 * t519 - t626;
t588 = pkin(10) * t419 + t411 * t572 + t412 * t576;
t558 = (-mrSges(3,1) * t577 + mrSges(3,2) * t573) * t603;
t555 = -mrSges(3,2) * t564 + mrSges(3,3) * t598;
t554 = mrSges(3,1) * t564 - mrSges(3,3) * t599;
t544 = -t567 * t556 - t620;
t543 = Ifges(3,5) * t564 + (Ifges(3,1) * t573 + Ifges(3,4) * t577) * t603;
t542 = Ifges(3,6) * t564 + (Ifges(3,4) * t573 + Ifges(3,2) * t577) * t603;
t541 = Ifges(3,3) * t564 + (Ifges(3,5) * t573 + Ifges(3,6) * t577) * t603;
t534 = -g(3) * t611 + t604;
t533 = -g(3) * t610 + t595;
t418 = m(3) * t534 - mrSges(3,2) * t563 + mrSges(3,3) * t560 - t554 * t564 + t558 * t598 + t419;
t413 = m(3) * t533 + mrSges(3,1) * t563 - mrSges(3,3) * t559 + t555 * t564 - t558 * t599 + t415;
t410 = mrSges(4,1) * t468 - mrSges(4,2) * t469 + Ifges(4,5) * t523 + Ifges(4,6) * t522 + Ifges(4,3) * t540 + pkin(3) * t582 + pkin(11) * t596 + t416 * t624 + t571 * t417 + t539 * t518 - t538 * t519;
t409 = mrSges(3,1) * t533 - mrSges(3,2) * t534 + Ifges(3,5) * t559 + Ifges(3,6) * t560 + Ifges(3,3) * t563 + pkin(2) * t415 + t568 * t410 + (t542 * t573 - t543 * t577) * t603 + t588 * t566;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t597 - mrSges(2,2) * t593 + (t541 * t598 + mrSges(3,2) * t544 - mrSges(3,3) * t533 + Ifges(3,1) * t559 + Ifges(3,4) * t560 + Ifges(3,5) * t563 + t576 * t411 - t572 * t412 - t564 * t542 + (-t414 * t566 - t415 * t568) * pkin(10)) * t611 + (-mrSges(3,1) * t544 + mrSges(3,3) * t534 + Ifges(3,4) * t559 + Ifges(3,2) * t560 + Ifges(3,6) * t563 - pkin(2) * t414 - t566 * t410 - t541 * t599 + t564 * t543 + t568 * t588) * t610 + t569 * t409 + pkin(1) * ((t413 * t577 + t418 * t573) * t569 + (-m(3) * t544 + t560 * mrSges(3,1) - t559 * mrSges(3,2) + (-t554 * t573 + t555 * t577) * t603 - t414) * t567) + (-t413 * t573 + t418 * t577) * t623; t409; t410; t626; t432; t584;];
tauJ  = t1;
