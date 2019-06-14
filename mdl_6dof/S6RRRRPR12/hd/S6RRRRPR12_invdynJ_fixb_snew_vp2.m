% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-08 00:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR12_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 00:16:03
% EndTime: 2019-05-08 00:16:56
% DurationCPUTime: 38.50s
% Computational Cost: add. (612958->380), mult. (1518606->511), div. (0->0), fcn. (1290104->16), ass. (0->169)
t569 = cos(pkin(6));
t562 = qJD(1) * t569 + qJD(2);
t565 = sin(pkin(7));
t568 = cos(pkin(7));
t566 = sin(pkin(6));
t578 = cos(qJ(2));
t599 = qJD(1) * t578;
t595 = t566 * t599;
t549 = (t562 * t565 + t568 * t595) * pkin(10);
t573 = sin(qJ(2));
t601 = qJD(1) * t566;
t613 = pkin(10) * t565;
t553 = (-pkin(2) * t578 - t573 * t613) * t601;
t598 = qJD(1) * qJD(2);
t559 = (qJDD(1) * t573 + t578 * t598) * t566;
t561 = qJDD(1) * t569 + qJDD(2);
t580 = qJD(1) ^ 2;
t574 = sin(qJ(1));
t579 = cos(qJ(1));
t590 = -g(1) * t579 - g(2) * t574;
t614 = pkin(9) * t566;
t557 = -pkin(1) * t580 + qJDD(1) * t614 + t590;
t594 = t574 * g(1) - g(2) * t579;
t556 = qJDD(1) * pkin(1) + t580 * t614 + t594;
t609 = t556 * t569;
t591 = -t573 * t557 + t578 * t609;
t600 = qJD(1) * t573;
t612 = pkin(10) * t568;
t510 = -t559 * t612 + t561 * pkin(2) + t562 * t549 + (-g(3) * t578 - t553 * t600) * t566 + t591;
t596 = t566 * t600;
t552 = pkin(2) * t562 - t596 * t612;
t560 = (qJDD(1) * t578 - t573 * t598) * t566;
t587 = t560 * t568 + t561 * t565;
t602 = t578 * t557 + t573 * t609;
t511 = -t562 * t552 + (-g(3) * t573 + t553 * t599) * t566 + t587 * pkin(10) + t602;
t611 = t569 * g(3);
t517 = -t559 * t613 - t560 * pkin(2) - t611 + (-t556 + (-t549 * t578 + t552 * t573) * qJD(1)) * t566;
t572 = sin(qJ(3));
t577 = cos(qJ(3));
t477 = -t572 * t511 + (t510 * t568 + t517 * t565) * t577;
t604 = t568 * t572;
t608 = t565 * t572;
t478 = t510 * t604 + t577 * t511 + t517 * t608;
t603 = t568 * t578;
t607 = t565 * t577;
t539 = (-t572 * t573 + t577 * t603) * t601 + t562 * t607;
t540 = t562 * t608 + (t572 * t603 + t573 * t577) * t601;
t529 = -pkin(3) * t539 - pkin(11) * t540;
t541 = -t560 * t565 + t561 * t568 + qJDD(3);
t550 = t562 * t568 - t565 * t595 + qJD(3);
t548 = t550 ^ 2;
t469 = -pkin(3) * t548 + pkin(11) * t541 + t529 * t539 + t478;
t492 = -t565 * t510 + t568 * t517;
t526 = -t540 * qJD(3) - t572 * t559 + t577 * t587;
t527 = t539 * qJD(3) + t577 * t559 + t572 * t587;
t472 = (-t539 * t550 - t527) * pkin(11) + (t540 * t550 - t526) * pkin(3) + t492;
t571 = sin(qJ(4));
t576 = cos(qJ(4));
t461 = -t571 * t469 + t576 * t472;
t531 = -t540 * t571 + t550 * t576;
t495 = qJD(4) * t531 + t527 * t576 + t541 * t571;
t525 = qJDD(4) - t526;
t532 = t540 * t576 + t550 * t571;
t538 = qJD(4) - t539;
t458 = (t531 * t538 - t495) * qJ(5) + (t531 * t532 + t525) * pkin(4) + t461;
t462 = t576 * t469 + t571 * t472;
t494 = -qJD(4) * t532 - t527 * t571 + t541 * t576;
t520 = pkin(4) * t538 - qJ(5) * t532;
t530 = t531 ^ 2;
t460 = -pkin(4) * t530 + qJ(5) * t494 - t520 * t538 + t462;
t564 = sin(pkin(13));
t567 = cos(pkin(13));
t512 = t531 * t567 - t532 * t564;
t615 = 2 * qJD(5);
t455 = t564 * t458 + t567 * t460 + t512 * t615;
t513 = t531 * t564 + t532 * t567;
t491 = -pkin(5) * t512 - pkin(12) * t513;
t537 = t538 ^ 2;
t453 = -pkin(5) * t537 + pkin(12) * t525 + t491 * t512 + t455;
t468 = -t541 * pkin(3) - t548 * pkin(11) + t540 * t529 - t477;
t463 = -t494 * pkin(4) - t530 * qJ(5) + t532 * t520 + qJDD(5) + t468;
t481 = t494 * t567 - t495 * t564;
t482 = t494 * t564 + t495 * t567;
t456 = (-t512 * t538 - t482) * pkin(12) + (t513 * t538 - t481) * pkin(5) + t463;
t570 = sin(qJ(6));
t575 = cos(qJ(6));
t450 = -t453 * t570 + t456 * t575;
t496 = -t513 * t570 + t538 * t575;
t466 = qJD(6) * t496 + t482 * t575 + t525 * t570;
t480 = qJDD(6) - t481;
t497 = t513 * t575 + t538 * t570;
t483 = -mrSges(7,1) * t496 + mrSges(7,2) * t497;
t509 = qJD(6) - t512;
t484 = -mrSges(7,2) * t509 + mrSges(7,3) * t496;
t447 = m(7) * t450 + mrSges(7,1) * t480 - mrSges(7,3) * t466 - t483 * t497 + t484 * t509;
t451 = t453 * t575 + t456 * t570;
t465 = -qJD(6) * t497 - t482 * t570 + t525 * t575;
t485 = mrSges(7,1) * t509 - mrSges(7,3) * t497;
t448 = m(7) * t451 - mrSges(7,2) * t480 + mrSges(7,3) * t465 + t483 * t496 - t485 * t509;
t439 = -t447 * t570 + t575 * t448;
t490 = -mrSges(6,1) * t512 + mrSges(6,2) * t513;
t499 = mrSges(6,1) * t538 - mrSges(6,3) * t513;
t433 = m(6) * t455 - mrSges(6,2) * t525 + mrSges(6,3) * t481 + t490 * t512 - t499 * t538 + t439;
t589 = -t567 * t458 + t564 * t460;
t452 = -t525 * pkin(5) - t537 * pkin(12) + (t615 + t491) * t513 + t589;
t449 = -m(7) * t452 + t465 * mrSges(7,1) - mrSges(7,2) * t466 + t496 * t484 - t485 * t497;
t454 = -0.2e1 * qJD(5) * t513 - t589;
t498 = -mrSges(6,2) * t538 + mrSges(6,3) * t512;
t443 = m(6) * t454 + mrSges(6,1) * t525 - mrSges(6,3) * t482 - t490 * t513 + t498 * t538 + t449;
t430 = t564 * t433 + t567 * t443;
t473 = Ifges(7,5) * t497 + Ifges(7,6) * t496 + Ifges(7,3) * t509;
t475 = Ifges(7,1) * t497 + Ifges(7,4) * t496 + Ifges(7,5) * t509;
t440 = -mrSges(7,1) * t452 + mrSges(7,3) * t451 + Ifges(7,4) * t466 + Ifges(7,2) * t465 + Ifges(7,6) * t480 - t473 * t497 + t475 * t509;
t474 = Ifges(7,4) * t497 + Ifges(7,2) * t496 + Ifges(7,6) * t509;
t441 = mrSges(7,2) * t452 - mrSges(7,3) * t450 + Ifges(7,1) * t466 + Ifges(7,4) * t465 + Ifges(7,5) * t480 + t473 * t496 - t474 * t509;
t487 = Ifges(6,4) * t513 + Ifges(6,2) * t512 + Ifges(6,6) * t538;
t488 = Ifges(6,1) * t513 + Ifges(6,4) * t512 + Ifges(6,5) * t538;
t501 = Ifges(5,4) * t532 + Ifges(5,2) * t531 + Ifges(5,6) * t538;
t502 = Ifges(5,1) * t532 + Ifges(5,4) * t531 + Ifges(5,5) * t538;
t616 = Ifges(5,5) * t495 + Ifges(5,6) * t494 + t532 * t501 - t531 * t502 + mrSges(5,1) * t461 - mrSges(5,2) * t462 + Ifges(6,5) * t482 + Ifges(6,6) * t481 + t513 * t487 - t512 * t488 + mrSges(6,1) * t454 - mrSges(6,2) * t455 + t570 * t441 + t575 * t440 + pkin(5) * t449 + pkin(12) * t439 + pkin(4) * t430 + (Ifges(5,3) + Ifges(6,3)) * t525;
t606 = t566 * t573;
t605 = t566 * t578;
t514 = -mrSges(5,1) * t531 + mrSges(5,2) * t532;
t519 = -mrSges(5,2) * t538 + mrSges(5,3) * t531;
t428 = m(5) * t461 + mrSges(5,1) * t525 - mrSges(5,3) * t495 - t514 * t532 + t519 * t538 + t430;
t521 = mrSges(5,1) * t538 - mrSges(5,3) * t532;
t592 = t567 * t433 - t443 * t564;
t429 = m(5) * t462 - mrSges(5,2) * t525 + mrSges(5,3) * t494 + t514 * t531 - t521 * t538 + t592;
t422 = t576 * t428 + t571 * t429;
t438 = t575 * t447 + t570 * t448;
t528 = -mrSges(4,1) * t539 + mrSges(4,2) * t540;
t534 = mrSges(4,1) * t550 - mrSges(4,3) * t540;
t593 = -t428 * t571 + t576 * t429;
t419 = m(4) * t478 - mrSges(4,2) * t541 + mrSges(4,3) * t526 + t528 * t539 - t534 * t550 + t593;
t533 = -mrSges(4,2) * t550 + mrSges(4,3) * t539;
t421 = m(4) * t492 - mrSges(4,1) * t526 + mrSges(4,2) * t527 - t533 * t539 + t534 * t540 + t422;
t437 = m(6) * t463 - t481 * mrSges(6,1) + mrSges(6,2) * t482 - t512 * t498 + t499 * t513 + t438;
t582 = -m(5) * t468 + t494 * mrSges(5,1) - mrSges(5,2) * t495 + t531 * t519 - t521 * t532 - t437;
t436 = m(4) * t477 + mrSges(4,1) * t541 - mrSges(4,3) * t527 - t528 * t540 + t533 * t550 + t582;
t410 = t419 * t608 + t568 * t421 + t436 * t607;
t415 = t577 * t419 - t436 * t572;
t411 = t568 * t577 * t436 + t419 * t604 - t421 * t565;
t486 = Ifges(6,5) * t513 + Ifges(6,6) * t512 + Ifges(6,3) * t538;
t423 = mrSges(6,2) * t463 - mrSges(6,3) * t454 + Ifges(6,1) * t482 + Ifges(6,4) * t481 + Ifges(6,5) * t525 - pkin(12) * t438 - t440 * t570 + t441 * t575 + t486 * t512 - t487 * t538;
t583 = mrSges(7,1) * t450 - mrSges(7,2) * t451 + Ifges(7,5) * t466 + Ifges(7,6) * t465 + Ifges(7,3) * t480 + t474 * t497 - t475 * t496;
t424 = -mrSges(6,1) * t463 + mrSges(6,3) * t455 + Ifges(6,4) * t482 + Ifges(6,2) * t481 + Ifges(6,6) * t525 - pkin(5) * t438 - t486 * t513 + t488 * t538 - t583;
t500 = Ifges(5,5) * t532 + Ifges(5,6) * t531 + Ifges(5,3) * t538;
t412 = -mrSges(5,1) * t468 + mrSges(5,3) * t462 + Ifges(5,4) * t495 + Ifges(5,2) * t494 + Ifges(5,6) * t525 - pkin(4) * t437 + qJ(5) * t592 + t564 * t423 + t567 * t424 - t532 * t500 + t538 * t502;
t413 = mrSges(5,2) * t468 - mrSges(5,3) * t461 + Ifges(5,1) * t495 + Ifges(5,4) * t494 + Ifges(5,5) * t525 - qJ(5) * t430 + t423 * t567 - t424 * t564 + t500 * t531 - t501 * t538;
t522 = Ifges(4,5) * t540 + Ifges(4,6) * t539 + Ifges(4,3) * t550;
t523 = Ifges(4,4) * t540 + Ifges(4,2) * t539 + Ifges(4,6) * t550;
t407 = mrSges(4,2) * t492 - mrSges(4,3) * t477 + Ifges(4,1) * t527 + Ifges(4,4) * t526 + Ifges(4,5) * t541 - pkin(11) * t422 - t412 * t571 + t413 * t576 + t522 * t539 - t523 * t550;
t524 = Ifges(4,1) * t540 + Ifges(4,4) * t539 + Ifges(4,5) * t550;
t408 = -mrSges(4,1) * t492 + mrSges(4,3) * t478 + Ifges(4,4) * t527 + Ifges(4,2) * t526 + Ifges(4,6) * t541 - pkin(3) * t422 - t540 * t522 + t550 * t524 - t616;
t584 = pkin(10) * t415 + t407 * t572 + t408 * t577;
t558 = (-mrSges(3,1) * t578 + mrSges(3,2) * t573) * t601;
t555 = -mrSges(3,2) * t562 + mrSges(3,3) * t595;
t554 = mrSges(3,1) * t562 - mrSges(3,3) * t596;
t545 = -t566 * t556 - t611;
t544 = Ifges(3,5) * t562 + (Ifges(3,1) * t573 + Ifges(3,4) * t578) * t601;
t543 = Ifges(3,6) * t562 + (Ifges(3,4) * t573 + Ifges(3,2) * t578) * t601;
t542 = Ifges(3,3) * t562 + (Ifges(3,5) * t573 + Ifges(3,6) * t578) * t601;
t536 = -g(3) * t606 + t602;
t535 = -g(3) * t605 + t591;
t414 = m(3) * t536 - mrSges(3,2) * t561 + mrSges(3,3) * t560 - t554 * t562 + t558 * t595 + t415;
t409 = m(3) * t535 + mrSges(3,1) * t561 - mrSges(3,3) * t559 + t555 * t562 - t558 * t596 + t411;
t406 = mrSges(4,1) * t477 - mrSges(4,2) * t478 + Ifges(4,5) * t527 + Ifges(4,6) * t526 + Ifges(4,3) * t541 + pkin(3) * t582 + pkin(11) * t593 + t576 * t412 + t571 * t413 + t540 * t523 - t539 * t524;
t405 = mrSges(3,1) * t535 - mrSges(3,2) * t536 + Ifges(3,5) * t559 + Ifges(3,6) * t560 + Ifges(3,3) * t561 + pkin(2) * t411 + t568 * t406 + (t543 * t573 - t544 * t578) * t601 + t584 * t565;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t594 - mrSges(2,2) * t590 + (t542 * t595 + mrSges(3,2) * t545 - mrSges(3,3) * t535 + Ifges(3,1) * t559 + Ifges(3,4) * t560 + Ifges(3,5) * t561 + t577 * t407 - t572 * t408 - t562 * t543 + (-t410 * t565 - t411 * t568) * pkin(10)) * t606 + (-mrSges(3,1) * t545 + mrSges(3,3) * t536 + Ifges(3,4) * t559 + Ifges(3,2) * t560 + Ifges(3,6) * t561 - pkin(2) * t410 - t565 * t406 - t542 * t596 + t562 * t544 + t568 * t584) * t605 + t569 * t405 + pkin(1) * ((t409 * t578 + t414 * t573) * t569 + (-m(3) * t545 + t560 * mrSges(3,1) - t559 * mrSges(3,2) + (-t554 * t573 + t555 * t578) * t601 - t410) * t566) + (-t409 * t573 + t414 * t578) * t614; t405; t406; t616; t437; t583;];
tauJ  = t1;
