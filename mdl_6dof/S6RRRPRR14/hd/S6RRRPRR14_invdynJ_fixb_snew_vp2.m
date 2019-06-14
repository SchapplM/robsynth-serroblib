% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR14
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 16:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR14_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR14_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 16:34:44
% EndTime: 2019-05-07 16:34:56
% DurationCPUTime: 6.80s
% Computational Cost: add. (94487->343), mult. (202244->428), div. (0->0), fcn. (157879->12), ass. (0->147)
t605 = Ifges(4,1) + Ifges(5,2);
t598 = Ifges(4,5) - Ifges(5,4);
t604 = -Ifges(4,2) - Ifges(5,3);
t597 = Ifges(4,6) - Ifges(5,5);
t596 = -Ifges(5,6) - Ifges(4,4);
t603 = Ifges(4,3) + Ifges(5,1);
t552 = sin(pkin(6));
t557 = sin(qJ(2));
t561 = cos(qJ(2));
t582 = qJD(1) * qJD(2);
t540 = (-qJDD(1) * t561 + t557 * t582) * t552;
t585 = qJD(1) * t552;
t538 = (-pkin(2) * t561 - pkin(9) * t557) * t585;
t553 = cos(pkin(6));
t549 = qJD(1) * t553 + qJD(2);
t547 = t549 ^ 2;
t548 = qJDD(1) * t553 + qJDD(2);
t584 = qJD(1) * t561;
t563 = qJD(1) ^ 2;
t558 = sin(qJ(1));
t562 = cos(qJ(1));
t576 = -g(1) * t562 - g(2) * t558;
t600 = pkin(8) * t552;
t536 = -pkin(1) * t563 + qJDD(1) * t600 + t576;
t579 = t558 * g(1) - g(2) * t562;
t535 = qJDD(1) * pkin(1) + t563 * t600 + t579;
t593 = t535 * t553;
t586 = t561 * t536 + t557 * t593;
t472 = -t547 * pkin(2) + t548 * pkin(9) + (-g(3) * t557 + t538 * t584) * t552 + t586;
t539 = (qJDD(1) * t557 + t561 * t582) * t552;
t599 = t553 * g(3);
t473 = t540 * pkin(2) - t539 * pkin(9) - t599 + (-t535 + (pkin(2) * t557 - pkin(9) * t561) * t549 * qJD(1)) * t552;
t556 = sin(qJ(3));
t601 = cos(qJ(3));
t450 = t601 * t472 + t556 * t473;
t581 = t557 * t585;
t526 = -t601 * t549 + t556 * t581;
t527 = t556 * t549 + t601 * t581;
t502 = pkin(3) * t526 - qJ(4) * t527;
t532 = qJDD(3) + t540;
t580 = t552 * t584;
t545 = -qJD(3) + t580;
t544 = t545 ^ 2;
t445 = pkin(3) * t544 - t532 * qJ(4) + 0.2e1 * qJD(4) * t545 + t526 * t502 - t450;
t449 = -t556 * t472 + t601 * t473;
t446 = -t532 * pkin(3) - t544 * qJ(4) + t527 * t502 + qJDD(4) - t449;
t499 = -t526 * qJD(3) + t601 * t539 + t556 * t548;
t594 = t526 * t545;
t435 = (t526 * t527 - t532) * pkin(10) + (t499 - t594) * pkin(4) + t446;
t498 = t527 * qJD(3) + t556 * t539 - t601 * t548;
t514 = pkin(4) * t527 + pkin(10) * t545;
t525 = t526 ^ 2;
t591 = t552 * t561;
t500 = -g(3) * t591 - t557 * t536 + t561 * t593;
t471 = -t548 * pkin(2) - t547 * pkin(9) + t538 * t581 - t500;
t567 = (-t499 - t594) * qJ(4) + t471 + (-pkin(3) * t545 - 0.2e1 * qJD(4)) * t527;
t439 = -t525 * pkin(4) - t527 * t514 + (pkin(3) + pkin(10)) * t498 + t567;
t555 = sin(qJ(5));
t560 = cos(qJ(5));
t429 = t560 * t435 - t555 * t439;
t508 = t526 * t560 + t545 * t555;
t461 = qJD(5) * t508 + t498 * t555 + t532 * t560;
t495 = qJDD(5) + t499;
t509 = t526 * t555 - t545 * t560;
t524 = qJD(5) + t527;
t426 = (t508 * t524 - t461) * pkin(11) + (t508 * t509 + t495) * pkin(5) + t429;
t430 = t555 * t435 + t560 * t439;
t460 = -qJD(5) * t509 + t498 * t560 - t532 * t555;
t482 = pkin(5) * t524 - pkin(11) * t509;
t507 = t508 ^ 2;
t427 = -pkin(5) * t507 + pkin(11) * t460 - t482 * t524 + t430;
t554 = sin(qJ(6));
t559 = cos(qJ(6));
t425 = t426 * t554 + t427 * t559;
t438 = -pkin(4) * t498 - pkin(10) * t525 - t545 * t514 - t445;
t432 = -pkin(5) * t460 - pkin(11) * t507 + t482 * t509 + t438;
t475 = t508 * t554 + t509 * t559;
t443 = -qJD(6) * t475 + t460 * t559 - t461 * t554;
t474 = t508 * t559 - t509 * t554;
t444 = qJD(6) * t474 + t460 * t554 + t461 * t559;
t522 = qJD(6) + t524;
t451 = Ifges(7,5) * t475 + Ifges(7,6) * t474 + Ifges(7,3) * t522;
t453 = Ifges(7,1) * t475 + Ifges(7,4) * t474 + Ifges(7,5) * t522;
t485 = qJDD(6) + t495;
t412 = -mrSges(7,1) * t432 + mrSges(7,3) * t425 + Ifges(7,4) * t444 + Ifges(7,2) * t443 + Ifges(7,6) * t485 - t451 * t475 + t453 * t522;
t424 = t426 * t559 - t427 * t554;
t452 = Ifges(7,4) * t475 + Ifges(7,2) * t474 + Ifges(7,6) * t522;
t413 = mrSges(7,2) * t432 - mrSges(7,3) * t424 + Ifges(7,1) * t444 + Ifges(7,4) * t443 + Ifges(7,5) * t485 + t451 * t474 - t452 * t522;
t464 = Ifges(6,5) * t509 + Ifges(6,6) * t508 + Ifges(6,3) * t524;
t466 = Ifges(6,1) * t509 + Ifges(6,4) * t508 + Ifges(6,5) * t524;
t462 = -mrSges(7,2) * t522 + mrSges(7,3) * t474;
t463 = mrSges(7,1) * t522 - mrSges(7,3) * t475;
t572 = m(7) * t432 - t443 * mrSges(7,1) + t444 * mrSges(7,2) - t474 * t462 + t475 * t463;
t455 = -mrSges(7,1) * t474 + mrSges(7,2) * t475;
t421 = m(7) * t424 + mrSges(7,1) * t485 - mrSges(7,3) * t444 - t455 * t475 + t462 * t522;
t422 = m(7) * t425 - mrSges(7,2) * t485 + mrSges(7,3) * t443 + t455 * t474 - t463 * t522;
t577 = -t421 * t554 + t559 * t422;
t398 = -mrSges(6,1) * t438 + mrSges(6,3) * t430 + Ifges(6,4) * t461 + Ifges(6,2) * t460 + Ifges(6,6) * t495 - pkin(5) * t572 + pkin(11) * t577 + t559 * t412 + t554 * t413 - t509 * t464 + t524 * t466;
t411 = t559 * t421 + t554 * t422;
t465 = Ifges(6,4) * t509 + Ifges(6,2) * t508 + Ifges(6,6) * t524;
t399 = mrSges(6,2) * t438 - mrSges(6,3) * t429 + Ifges(6,1) * t461 + Ifges(6,4) * t460 + Ifges(6,5) * t495 - pkin(11) * t411 - t412 * t554 + t413 * t559 + t464 * t508 - t465 * t524;
t512 = mrSges(5,1) * t526 + mrSges(5,3) * t545;
t476 = -mrSges(6,1) * t508 + mrSges(6,2) * t509;
t480 = -mrSges(6,2) * t524 + mrSges(6,3) * t508;
t408 = m(6) * t429 + mrSges(6,1) * t495 - mrSges(6,3) * t461 - t476 * t509 + t480 * t524 + t411;
t481 = mrSges(6,1) * t524 - mrSges(6,3) * t509;
t409 = m(6) * t430 - mrSges(6,2) * t495 + mrSges(6,3) * t460 + t476 * t508 - t481 * t524 + t577;
t405 = t560 * t408 + t555 * t409;
t504 = -mrSges(5,2) * t526 - mrSges(5,3) * t527;
t571 = -m(5) * t446 - t499 * mrSges(5,1) - t527 * t504 - t405;
t404 = t532 * mrSges(5,2) - t545 * t512 - t571;
t513 = mrSges(5,1) * t527 - mrSges(5,2) * t545;
t569 = -m(6) * t438 + t460 * mrSges(6,1) - t461 * mrSges(6,2) + t508 * t480 - t509 * t481 - t572;
t566 = -m(5) * t445 + t532 * mrSges(5,3) - t545 * t513 - t569;
t587 = t596 * t526 + t605 * t527 - t598 * t545;
t588 = t604 * t526 - t596 * t527 - t597 * t545;
t602 = -t597 * t498 + t598 * t499 + t587 * t526 + t588 * t527 + t603 * t532 + mrSges(4,1) * t449 - mrSges(4,2) * t450 + mrSges(5,2) * t446 - mrSges(5,3) * t445 - pkin(3) * t404 - pkin(10) * t405 + qJ(4) * (-t498 * mrSges(5,1) - t526 * t504 + t566) - t555 * t398 + t560 * t399;
t592 = t552 * t557;
t590 = -t555 * t408 + t560 * t409;
t503 = mrSges(4,1) * t526 + mrSges(4,2) * t527;
t510 = mrSges(4,2) * t545 - mrSges(4,3) * t526;
t402 = m(4) * t449 - t499 * mrSges(4,3) - t527 * t503 + (-t510 + t512) * t545 + (mrSges(4,1) - mrSges(5,2)) * t532 + t571;
t511 = -mrSges(4,1) * t545 - mrSges(4,3) * t527;
t416 = (-t503 - t504) * t526 + t566 + (-mrSges(4,3) - mrSges(5,1)) * t498 + t545 * t511 - t532 * mrSges(4,2) + m(4) * t450;
t397 = t601 * t402 + t556 * t416;
t589 = t597 * t526 - t598 * t527 + t603 * t545;
t578 = -t402 * t556 + t601 * t416;
t447 = t498 * pkin(3) + t567;
t575 = -m(5) * t447 + t498 * mrSges(5,2) + t526 * t512 - t590;
t570 = mrSges(7,1) * t424 - mrSges(7,2) * t425 + Ifges(7,5) * t444 + Ifges(7,6) * t443 + Ifges(7,3) * t485 + t475 * t452 - t474 * t453;
t568 = -m(4) * t471 - t498 * mrSges(4,1) - t526 * t510 + (-t511 + t513) * t527 + (-mrSges(4,2) + mrSges(5,3)) * t499 + t575;
t565 = mrSges(6,1) * t429 - mrSges(6,2) * t430 + Ifges(6,5) * t461 + Ifges(6,6) * t460 + Ifges(6,3) * t495 + pkin(5) * t411 + t509 * t465 - t508 * t466 + t570;
t537 = (-mrSges(3,1) * t561 + mrSges(3,2) * t557) * t585;
t534 = -mrSges(3,2) * t549 + mrSges(3,3) * t580;
t533 = mrSges(3,1) * t549 - mrSges(3,3) * t581;
t519 = -t552 * t535 - t599;
t518 = Ifges(3,5) * t549 + (Ifges(3,1) * t557 + Ifges(3,4) * t561) * t585;
t517 = Ifges(3,6) * t549 + (Ifges(3,4) * t557 + Ifges(3,2) * t561) * t585;
t516 = Ifges(3,3) * t549 + (Ifges(3,5) * t557 + Ifges(3,6) * t561) * t585;
t501 = -g(3) * t592 + t586;
t403 = -t499 * mrSges(5,3) - t527 * t513 - t575;
t400 = m(3) * t500 + t548 * mrSges(3,1) - t539 * mrSges(3,3) + t549 * t534 - t537 * t581 + t568;
t396 = m(3) * t501 - mrSges(3,2) * t548 - mrSges(3,3) * t540 - t533 * t549 + t537 * t580 + t578;
t395 = mrSges(5,1) * t446 + mrSges(4,2) * t471 - mrSges(4,3) * t449 - mrSges(5,3) * t447 + pkin(4) * t405 - qJ(4) * t403 + t596 * t498 + t605 * t499 + t589 * t526 + t598 * t532 + t588 * t545 + t565;
t394 = -mrSges(4,1) * t471 - mrSges(5,1) * t445 + mrSges(5,2) * t447 + mrSges(4,3) * t450 - pkin(3) * t403 - pkin(4) * t569 - pkin(10) * t590 - t560 * t398 - t555 * t399 + t604 * t498 - t596 * t499 + t589 * t527 + t597 * t532 - t587 * t545;
t393 = Ifges(3,5) * t539 - Ifges(3,6) * t540 + Ifges(3,3) * t548 + mrSges(3,1) * t500 - mrSges(3,2) * t501 + t556 * t395 + t601 * t394 + pkin(2) * t568 + pkin(9) * t578 + (t517 * t557 - t518 * t561) * t585;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t579 - mrSges(2,2) * t576 + (mrSges(3,2) * t519 - mrSges(3,3) * t500 + Ifges(3,1) * t539 - Ifges(3,4) * t540 + Ifges(3,5) * t548 - pkin(9) * t397 - t556 * t394 + t601 * t395 + t516 * t580 - t549 * t517) * t592 + (-mrSges(3,1) * t519 + mrSges(3,3) * t501 + Ifges(3,4) * t539 - Ifges(3,2) * t540 + Ifges(3,6) * t548 - pkin(2) * t397 - t516 * t581 + t549 * t518 - t602) * t591 + t553 * t393 + pkin(1) * ((t396 * t557 + t400 * t561) * t553 + (-m(3) * t519 - t540 * mrSges(3,1) - t539 * mrSges(3,2) + (-t533 * t557 + t534 * t561) * t585 - t397) * t552) + (t396 * t561 - t400 * t557) * t600; t393; t602; t404; t565; t570;];
tauJ  = t1;
