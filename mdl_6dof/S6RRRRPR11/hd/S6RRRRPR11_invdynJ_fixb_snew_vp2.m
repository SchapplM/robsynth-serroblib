% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 23:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR11_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR11_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR11_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 23:25:58
% EndTime: 2019-05-07 23:26:23
% DurationCPUTime: 19.35s
% Computational Cost: add. (315539->362), mult. (678906->472), div. (0->0), fcn. (550224->14), ass. (0->152)
t543 = sin(pkin(6));
t549 = sin(qJ(2));
t554 = cos(qJ(2));
t570 = qJD(1) * qJD(2);
t533 = (-qJDD(1) * t554 + t549 * t570) * t543;
t572 = qJD(1) * t543;
t531 = (-pkin(2) * t554 - pkin(9) * t549) * t572;
t545 = cos(pkin(6));
t539 = qJD(1) * t545 + qJD(2);
t537 = t539 ^ 2;
t538 = qJDD(1) * t545 + qJDD(2);
t571 = qJD(1) * t554;
t556 = qJD(1) ^ 2;
t550 = sin(qJ(1));
t555 = cos(qJ(1));
t563 = -g(1) * t555 - g(2) * t550;
t579 = pkin(8) * t543;
t529 = -pkin(1) * t556 + qJDD(1) * t579 + t563;
t567 = t550 * g(1) - g(2) * t555;
t528 = qJDD(1) * pkin(1) + t556 * t579 + t567;
t576 = t528 * t545;
t573 = t554 * t529 + t549 * t576;
t480 = -t537 * pkin(2) + t538 * pkin(9) + (-g(3) * t549 + t531 * t571) * t543 + t573;
t532 = (qJDD(1) * t549 + t554 * t570) * t543;
t578 = t545 * g(3);
t481 = t533 * pkin(2) - t532 * pkin(9) - t578 + (-t528 + (pkin(2) * t549 - pkin(9) * t554) * t539 * qJD(1)) * t543;
t548 = sin(qJ(3));
t553 = cos(qJ(3));
t457 = t553 * t480 + t548 * t481;
t569 = t549 * t572;
t520 = t539 * t553 - t548 * t569;
t521 = t539 * t548 + t553 * t569;
t505 = -pkin(3) * t520 - pkin(10) * t521;
t525 = qJDD(3) + t533;
t568 = t543 * t571;
t536 = qJD(3) - t568;
t534 = t536 ^ 2;
t447 = -pkin(3) * t534 + pkin(10) * t525 + t505 * t520 + t457;
t574 = t543 * t554;
t502 = -g(3) * t574 - t549 * t529 + t554 * t576;
t479 = -t538 * pkin(2) - t537 * pkin(9) + t531 * t569 - t502;
t500 = -t521 * qJD(3) - t548 * t532 + t538 * t553;
t501 = qJD(3) * t520 + t532 * t553 + t538 * t548;
t450 = (-t520 * t536 - t501) * pkin(10) + (t521 * t536 - t500) * pkin(3) + t479;
t547 = sin(qJ(4));
t552 = cos(qJ(4));
t435 = -t547 * t447 + t552 * t450;
t507 = -t521 * t547 + t536 * t552;
t467 = qJD(4) * t507 + t501 * t552 + t525 * t547;
t498 = qJDD(4) - t500;
t508 = t521 * t552 + t536 * t547;
t519 = qJD(4) - t520;
t426 = (t507 * t519 - t467) * qJ(5) + (t507 * t508 + t498) * pkin(4) + t435;
t436 = t552 * t447 + t547 * t450;
t466 = -qJD(4) * t508 - t501 * t547 + t525 * t552;
t490 = pkin(4) * t519 - qJ(5) * t508;
t506 = t507 ^ 2;
t428 = -pkin(4) * t506 + qJ(5) * t466 - t490 * t519 + t436;
t542 = sin(pkin(12));
t544 = cos(pkin(12));
t486 = t507 * t542 + t508 * t544;
t420 = -0.2e1 * qJD(5) * t486 + t544 * t426 - t542 * t428;
t453 = t466 * t542 + t467 * t544;
t485 = t507 * t544 - t508 * t542;
t418 = (t485 * t519 - t453) * pkin(11) + (t485 * t486 + t498) * pkin(5) + t420;
t421 = 0.2e1 * qJD(5) * t485 + t542 * t426 + t544 * t428;
t452 = t466 * t544 - t467 * t542;
t470 = pkin(5) * t519 - pkin(11) * t486;
t484 = t485 ^ 2;
t419 = -pkin(5) * t484 + pkin(11) * t452 - t470 * t519 + t421;
t546 = sin(qJ(6));
t551 = cos(qJ(6));
t415 = t418 * t551 - t419 * t546;
t462 = t485 * t551 - t486 * t546;
t434 = qJD(6) * t462 + t452 * t546 + t453 * t551;
t463 = t485 * t546 + t486 * t551;
t444 = -mrSges(7,1) * t462 + mrSges(7,2) * t463;
t517 = qJD(6) + t519;
t454 = -mrSges(7,2) * t517 + mrSges(7,3) * t462;
t493 = qJDD(6) + t498;
t409 = m(7) * t415 + mrSges(7,1) * t493 - mrSges(7,3) * t434 - t444 * t463 + t454 * t517;
t416 = t418 * t546 + t419 * t551;
t433 = -qJD(6) * t463 + t452 * t551 - t453 * t546;
t455 = mrSges(7,1) * t517 - mrSges(7,3) * t463;
t410 = m(7) * t416 - mrSges(7,2) * t493 + mrSges(7,3) * t433 + t444 * t462 - t455 * t517;
t403 = t551 * t409 + t546 * t410;
t464 = -mrSges(6,1) * t485 + mrSges(6,2) * t486;
t468 = -mrSges(6,2) * t519 + mrSges(6,3) * t485;
t401 = m(6) * t420 + mrSges(6,1) * t498 - mrSges(6,3) * t453 - t464 * t486 + t468 * t519 + t403;
t469 = mrSges(6,1) * t519 - mrSges(6,3) * t486;
t564 = -t409 * t546 + t551 * t410;
t402 = m(6) * t421 - mrSges(6,2) * t498 + mrSges(6,3) * t452 + t464 * t485 - t469 * t519 + t564;
t397 = t544 * t401 + t542 * t402;
t459 = Ifges(6,4) * t486 + Ifges(6,2) * t485 + Ifges(6,6) * t519;
t460 = Ifges(6,1) * t486 + Ifges(6,4) * t485 + Ifges(6,5) * t519;
t472 = Ifges(5,4) * t508 + Ifges(5,2) * t507 + Ifges(5,6) * t519;
t473 = Ifges(5,1) * t508 + Ifges(5,4) * t507 + Ifges(5,5) * t519;
t440 = Ifges(7,4) * t463 + Ifges(7,2) * t462 + Ifges(7,6) * t517;
t441 = Ifges(7,1) * t463 + Ifges(7,4) * t462 + Ifges(7,5) * t517;
t560 = -mrSges(7,1) * t415 + mrSges(7,2) * t416 - Ifges(7,5) * t434 - Ifges(7,6) * t433 - Ifges(7,3) * t493 - t463 * t440 + t462 * t441;
t580 = mrSges(5,1) * t435 + mrSges(6,1) * t420 - mrSges(5,2) * t436 - mrSges(6,2) * t421 + Ifges(5,5) * t467 + Ifges(6,5) * t453 + Ifges(5,6) * t466 + Ifges(6,6) * t452 + pkin(4) * t397 + pkin(5) * t403 + t486 * t459 - t485 * t460 + t508 * t472 - t507 * t473 + (Ifges(5,3) + Ifges(6,3)) * t498 - t560;
t575 = t543 * t549;
t487 = -mrSges(5,1) * t507 + mrSges(5,2) * t508;
t489 = -mrSges(5,2) * t519 + mrSges(5,3) * t507;
t395 = m(5) * t435 + mrSges(5,1) * t498 - mrSges(5,3) * t467 - t487 * t508 + t489 * t519 + t397;
t491 = mrSges(5,1) * t519 - mrSges(5,3) * t508;
t565 = -t401 * t542 + t544 * t402;
t396 = m(5) * t436 - mrSges(5,2) * t498 + mrSges(5,3) * t466 + t487 * t507 - t491 * t519 + t565;
t391 = -t395 * t547 + t552 * t396;
t504 = -mrSges(4,1) * t520 + mrSges(4,2) * t521;
t510 = mrSges(4,1) * t536 - mrSges(4,3) * t521;
t389 = m(4) * t457 - mrSges(4,2) * t525 + mrSges(4,3) * t500 + t504 * t520 - t510 * t536 + t391;
t456 = -t548 * t480 + t481 * t553;
t446 = -pkin(3) * t525 - pkin(10) * t534 + t521 * t505 - t456;
t437 = -pkin(4) * t466 - qJ(5) * t506 + t508 * t490 + qJDD(5) + t446;
t423 = -pkin(5) * t452 - pkin(11) * t484 + t470 * t486 + t437;
t562 = m(7) * t423 - t433 * mrSges(7,1) + t434 * mrSges(7,2) - t462 * t454 + t463 * t455;
t417 = m(6) * t437 - t452 * mrSges(6,1) + mrSges(6,2) * t453 - t485 * t468 + t469 * t486 + t562;
t413 = -m(5) * t446 + t466 * mrSges(5,1) - mrSges(5,2) * t467 + t507 * t489 - t491 * t508 - t417;
t509 = -mrSges(4,2) * t536 + mrSges(4,3) * t520;
t412 = m(4) * t456 + mrSges(4,1) * t525 - mrSges(4,3) * t501 - t504 * t521 + t509 * t536 + t413;
t385 = t548 * t389 + t553 * t412;
t566 = t553 * t389 - t412 * t548;
t390 = t395 * t552 + t396 * t547;
t559 = -m(4) * t479 + t500 * mrSges(4,1) - mrSges(4,2) * t501 + t520 * t509 - t510 * t521 - t390;
t439 = Ifges(7,5) * t463 + Ifges(7,6) * t462 + Ifges(7,3) * t517;
t404 = -mrSges(7,1) * t423 + mrSges(7,3) * t416 + Ifges(7,4) * t434 + Ifges(7,2) * t433 + Ifges(7,6) * t493 - t439 * t463 + t441 * t517;
t405 = mrSges(7,2) * t423 - mrSges(7,3) * t415 + Ifges(7,1) * t434 + Ifges(7,4) * t433 + Ifges(7,5) * t493 + t439 * t462 - t440 * t517;
t458 = Ifges(6,5) * t486 + Ifges(6,6) * t485 + Ifges(6,3) * t519;
t392 = -mrSges(6,1) * t437 + mrSges(6,3) * t421 + Ifges(6,4) * t453 + Ifges(6,2) * t452 + Ifges(6,6) * t498 - pkin(5) * t562 + pkin(11) * t564 + t551 * t404 + t546 * t405 - t486 * t458 + t519 * t460;
t393 = mrSges(6,2) * t437 - mrSges(6,3) * t420 + Ifges(6,1) * t453 + Ifges(6,4) * t452 + Ifges(6,5) * t498 - pkin(11) * t403 - t404 * t546 + t405 * t551 + t458 * t485 - t459 * t519;
t471 = Ifges(5,5) * t508 + Ifges(5,6) * t507 + Ifges(5,3) * t519;
t382 = -mrSges(5,1) * t446 + mrSges(5,3) * t436 + Ifges(5,4) * t467 + Ifges(5,2) * t466 + Ifges(5,6) * t498 - pkin(4) * t417 + qJ(5) * t565 + t544 * t392 + t542 * t393 - t508 * t471 + t519 * t473;
t383 = mrSges(5,2) * t446 - mrSges(5,3) * t435 + Ifges(5,1) * t467 + Ifges(5,4) * t466 + Ifges(5,5) * t498 - qJ(5) * t397 - t392 * t542 + t393 * t544 + t471 * t507 - t472 * t519;
t495 = Ifges(4,4) * t521 + Ifges(4,2) * t520 + Ifges(4,6) * t536;
t496 = Ifges(4,1) * t521 + Ifges(4,4) * t520 + Ifges(4,5) * t536;
t558 = mrSges(4,1) * t456 - mrSges(4,2) * t457 + Ifges(4,5) * t501 + Ifges(4,6) * t500 + Ifges(4,3) * t525 + pkin(3) * t413 + pkin(10) * t391 + t552 * t382 + t547 * t383 + t521 * t495 - t520 * t496;
t530 = (-mrSges(3,1) * t554 + mrSges(3,2) * t549) * t572;
t527 = -mrSges(3,2) * t539 + mrSges(3,3) * t568;
t526 = mrSges(3,1) * t539 - mrSges(3,3) * t569;
t514 = -t543 * t528 - t578;
t513 = Ifges(3,5) * t539 + (Ifges(3,1) * t549 + Ifges(3,4) * t554) * t572;
t512 = Ifges(3,6) * t539 + (Ifges(3,4) * t549 + Ifges(3,2) * t554) * t572;
t511 = Ifges(3,3) * t539 + (Ifges(3,5) * t549 + Ifges(3,6) * t554) * t572;
t503 = -g(3) * t575 + t573;
t494 = Ifges(4,5) * t521 + Ifges(4,6) * t520 + Ifges(4,3) * t536;
t386 = m(3) * t502 + mrSges(3,1) * t538 - mrSges(3,3) * t532 + t527 * t539 - t530 * t569 + t559;
t384 = m(3) * t503 - mrSges(3,2) * t538 - mrSges(3,3) * t533 - t526 * t539 + t530 * t568 + t566;
t381 = -mrSges(4,1) * t479 + mrSges(4,3) * t457 + Ifges(4,4) * t501 + Ifges(4,2) * t500 + Ifges(4,6) * t525 - pkin(3) * t390 - t521 * t494 + t536 * t496 - t580;
t380 = mrSges(4,2) * t479 - mrSges(4,3) * t456 + Ifges(4,1) * t501 + Ifges(4,4) * t500 + Ifges(4,5) * t525 - pkin(10) * t390 - t382 * t547 + t383 * t552 + t494 * t520 - t495 * t536;
t379 = Ifges(3,5) * t532 - Ifges(3,6) * t533 + Ifges(3,3) * t538 + mrSges(3,1) * t502 - mrSges(3,2) * t503 + t548 * t380 + t553 * t381 + pkin(2) * t559 + pkin(9) * t566 + (t512 * t549 - t513 * t554) * t572;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t567 - mrSges(2,2) * t563 + (mrSges(3,2) * t514 - mrSges(3,3) * t502 + Ifges(3,1) * t532 - Ifges(3,4) * t533 + Ifges(3,5) * t538 - pkin(9) * t385 + t380 * t553 - t381 * t548 + t511 * t568 - t512 * t539) * t575 + (-mrSges(3,1) * t514 + mrSges(3,3) * t503 + Ifges(3,4) * t532 - Ifges(3,2) * t533 + Ifges(3,6) * t538 - pkin(2) * t385 - t511 * t569 + t539 * t513 - t558) * t574 + t545 * t379 + pkin(1) * ((t384 * t549 + t386 * t554) * t545 + (-m(3) * t514 - t533 * mrSges(3,1) - t532 * mrSges(3,2) + (-t526 * t549 + t527 * t554) * t572 - t385) * t543) + (t384 * t554 - t386 * t549) * t579; t379; t558; t580; t417; -t560;];
tauJ  = t1;
