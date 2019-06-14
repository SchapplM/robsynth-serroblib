% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:40:38
% EndTime: 2019-05-06 18:40:47
% DurationCPUTime: 7.68s
% Computational Cost: add. (103494->337), mult. (234344->429), div. (0->0), fcn. (187937->12), ass. (0->142)
t574 = Ifges(6,1) + Ifges(7,1);
t566 = Ifges(6,4) - Ifges(7,5);
t565 = -Ifges(6,5) - Ifges(7,4);
t573 = Ifges(6,2) + Ifges(7,3);
t564 = Ifges(6,6) - Ifges(7,6);
t572 = -Ifges(6,3) - Ifges(7,2);
t532 = cos(pkin(6));
t526 = qJD(1) * t532 + qJD(2);
t529 = sin(pkin(11));
t531 = cos(pkin(11));
t535 = sin(qJ(2));
t530 = sin(pkin(6));
t554 = qJD(1) * t530;
t550 = t535 * t554;
t507 = t526 * t531 - t529 * t550;
t508 = t526 * t529 + t531 * t550;
t534 = sin(qJ(4));
t537 = cos(qJ(4));
t488 = t507 * t537 - t508 * t534;
t538 = cos(qJ(2));
t553 = qJD(1) * t538;
t518 = (qJD(2) * t553 + qJDD(1) * t535) * t530;
t525 = qJDD(1) * t532 + qJDD(2);
t494 = -t518 * t529 + t525 * t531;
t495 = t518 * t531 + t525 * t529;
t458 = qJD(4) * t488 + t494 * t534 + t495 * t537;
t489 = t507 * t534 + t508 * t537;
t549 = t530 * t553;
t522 = qJD(4) - t549;
t533 = sin(qJ(5));
t570 = cos(qJ(5));
t475 = t489 * t533 - t570 * t522;
t552 = qJDD(1) * t530;
t519 = -qJD(2) * t550 + t538 * t552;
t511 = qJDD(4) - t519;
t437 = -t475 * qJD(5) + t570 * t458 + t533 * t511;
t476 = t570 * t489 + t533 * t522;
t451 = mrSges(7,1) * t475 - mrSges(7,3) * t476;
t516 = (-pkin(2) * t538 - qJ(3) * t535) * t554;
t524 = t526 ^ 2;
t540 = qJD(1) ^ 2;
t536 = sin(qJ(1));
t539 = cos(qJ(1));
t544 = -g(1) * t539 - g(2) * t536;
t515 = -pkin(1) * t540 + pkin(8) * t552 + t544;
t548 = t536 * g(1) - g(2) * t539;
t569 = pkin(8) * t530;
t514 = qJDD(1) * pkin(1) + t540 * t569 + t548;
t562 = t514 * t532;
t555 = t538 * t515 + t535 * t562;
t473 = -t524 * pkin(2) + t525 * qJ(3) + (-g(3) * t535 + t516 * t553) * t530 + t555;
t568 = t532 * g(3);
t474 = -t519 * pkin(2) - t568 - t518 * qJ(3) + (-t514 + (pkin(2) * t535 - qJ(3) * t538) * t526 * qJD(1)) * t530;
t432 = -0.2e1 * qJD(3) * t508 - t529 * t473 + t531 * t474;
t428 = (-t507 * t549 - t495) * pkin(9) + (t507 * t508 - t519) * pkin(3) + t432;
t433 = 0.2e1 * qJD(3) * t507 + t531 * t473 + t529 * t474;
t496 = -pkin(3) * t549 - pkin(9) * t508;
t505 = t507 ^ 2;
t431 = -pkin(3) * t505 + pkin(9) * t494 + t496 * t549 + t433;
t424 = t534 * t428 + t537 * t431;
t468 = -pkin(4) * t488 - pkin(10) * t489;
t521 = t522 ^ 2;
t422 = -pkin(4) * t521 + pkin(10) * t511 + t468 * t488 + t424;
t560 = t530 * t538;
t484 = -g(3) * t560 - t535 * t515 + t538 * t562;
t472 = -t525 * pkin(2) - t524 * qJ(3) + t516 * t550 + qJDD(3) - t484;
t438 = -t494 * pkin(3) - t505 * pkin(9) + t508 * t496 + t472;
t457 = -qJD(4) * t489 + t494 * t537 - t495 * t534;
t426 = (-t488 * t522 - t458) * pkin(10) + (t489 * t522 - t457) * pkin(4) + t438;
t418 = -t533 * t422 + t570 * t426;
t450 = pkin(5) * t475 - qJ(6) * t476;
t456 = qJDD(5) - t457;
t487 = qJD(5) - t488;
t486 = t487 ^ 2;
t416 = -t456 * pkin(5) - t486 * qJ(6) + t476 * t450 + qJDD(6) - t418;
t459 = -mrSges(7,2) * t475 + mrSges(7,3) * t487;
t545 = -m(7) * t416 + t456 * mrSges(7,1) + t487 * t459;
t412 = t437 * mrSges(7,2) + t476 * t451 - t545;
t419 = t570 * t422 + t533 * t426;
t415 = -pkin(5) * t486 + qJ(6) * t456 + 0.2e1 * qJD(6) * t487 - t450 * t475 + t419;
t436 = qJD(5) * t476 + t458 * t533 - t570 * t511;
t462 = -mrSges(7,1) * t487 + mrSges(7,2) * t476;
t551 = m(7) * t415 + t456 * mrSges(7,3) + t487 * t462;
t557 = t566 * t475 - t574 * t476 + t565 * t487;
t558 = t573 * t475 - t566 * t476 - t564 * t487;
t571 = -t564 * t436 - t565 * t437 - t572 * t456 - t557 * t475 - t558 * t476 + mrSges(6,1) * t418 - mrSges(7,1) * t416 - mrSges(6,2) * t419 + mrSges(7,3) * t415 - pkin(5) * t412 + qJ(6) * (-t436 * mrSges(7,2) - t475 * t451 + t551);
t567 = -mrSges(6,3) - mrSges(7,2);
t561 = t530 * t535;
t461 = mrSges(6,1) * t487 - mrSges(6,3) * t476;
t556 = -mrSges(6,1) * t475 - mrSges(6,2) * t476 - t451;
t407 = m(6) * t419 - t456 * mrSges(6,2) + t567 * t436 - t487 * t461 + t556 * t475 + t551;
t460 = -mrSges(6,2) * t487 - mrSges(6,3) * t475;
t409 = m(6) * t418 + t456 * mrSges(6,1) + t567 * t437 + t487 * t460 + t556 * t476 + t545;
t402 = t570 * t407 - t409 * t533;
t467 = -mrSges(5,1) * t488 + mrSges(5,2) * t489;
t478 = mrSges(5,1) * t522 - mrSges(5,3) * t489;
t397 = m(5) * t424 - mrSges(5,2) * t511 + mrSges(5,3) * t457 + t467 * t488 - t478 * t522 + t402;
t423 = t537 * t428 - t534 * t431;
t421 = -t511 * pkin(4) - t521 * pkin(10) + t489 * t468 - t423;
t417 = -0.2e1 * qJD(6) * t476 + (t475 * t487 - t437) * qJ(6) + (t476 * t487 + t436) * pkin(5) + t421;
t413 = m(7) * t417 + mrSges(7,1) * t436 - t437 * mrSges(7,3) + t459 * t475 - t476 * t462;
t410 = -m(6) * t421 - t436 * mrSges(6,1) - mrSges(6,2) * t437 - t475 * t460 - t461 * t476 - t413;
t477 = -mrSges(5,2) * t522 + mrSges(5,3) * t488;
t404 = m(5) * t423 + mrSges(5,1) * t511 - mrSges(5,3) * t458 - t467 * t489 + t477 * t522 + t410;
t393 = t534 * t397 + t537 * t404;
t490 = -mrSges(4,1) * t507 + mrSges(4,2) * t508;
t492 = mrSges(4,2) * t549 + mrSges(4,3) * t507;
t391 = m(4) * t432 - mrSges(4,1) * t519 - mrSges(4,3) * t495 - t490 * t508 - t492 * t549 + t393;
t493 = -mrSges(4,1) * t549 - mrSges(4,3) * t508;
t546 = t537 * t397 - t404 * t534;
t392 = m(4) * t433 + mrSges(4,2) * t519 + mrSges(4,3) * t494 + t490 * t507 + t493 * t549 + t546;
t385 = t531 * t391 + t529 * t392;
t401 = t533 * t407 + t570 * t409;
t559 = t564 * t475 + t565 * t476 + t572 * t487;
t547 = -t391 * t529 + t531 * t392;
t543 = m(5) * t438 - t457 * mrSges(5,1) + mrSges(5,2) * t458 - t488 * t477 + t478 * t489 + t401;
t398 = m(4) * t472 - t494 * mrSges(4,1) + mrSges(4,2) * t495 - t507 * t492 + t493 * t508 + t543;
t399 = -mrSges(6,1) * t421 - mrSges(7,1) * t417 + mrSges(7,2) * t415 + mrSges(6,3) * t419 - pkin(5) * t413 - t573 * t436 + t566 * t437 + t564 * t456 + t559 * t476 - t557 * t487;
t400 = mrSges(6,2) * t421 + mrSges(7,2) * t416 - mrSges(6,3) * t418 - mrSges(7,3) * t417 - qJ(6) * t413 - t566 * t436 + t574 * t437 - t565 * t456 + t559 * t475 + t558 * t487;
t464 = Ifges(5,4) * t489 + Ifges(5,2) * t488 + Ifges(5,6) * t522;
t465 = Ifges(5,1) * t489 + Ifges(5,4) * t488 + Ifges(5,5) * t522;
t541 = mrSges(5,1) * t423 - mrSges(5,2) * t424 + Ifges(5,5) * t458 + Ifges(5,6) * t457 + Ifges(5,3) * t511 + pkin(4) * t410 + pkin(10) * t402 + t570 * t399 + t533 * t400 + t489 * t464 - t488 * t465;
t517 = (-mrSges(3,1) * t538 + mrSges(3,2) * t535) * t554;
t513 = -mrSges(3,2) * t526 + mrSges(3,3) * t549;
t512 = mrSges(3,1) * t526 - mrSges(3,3) * t550;
t500 = -t530 * t514 - t568;
t499 = Ifges(3,5) * t526 + (Ifges(3,1) * t535 + Ifges(3,4) * t538) * t554;
t498 = Ifges(3,6) * t526 + (t535 * Ifges(3,4) + Ifges(3,2) * t538) * t554;
t497 = Ifges(3,3) * t526 + (Ifges(3,5) * t535 + Ifges(3,6) * t538) * t554;
t485 = -g(3) * t561 + t555;
t481 = Ifges(4,1) * t508 + Ifges(4,4) * t507 - Ifges(4,5) * t549;
t480 = Ifges(4,4) * t508 + Ifges(4,2) * t507 - Ifges(4,6) * t549;
t479 = Ifges(4,5) * t508 + Ifges(4,6) * t507 - Ifges(4,3) * t549;
t463 = Ifges(5,5) * t489 + Ifges(5,6) * t488 + Ifges(5,3) * t522;
t394 = m(3) * t484 + mrSges(3,1) * t525 - mrSges(3,3) * t518 + t513 * t526 - t517 * t550 - t398;
t387 = -mrSges(5,1) * t438 + mrSges(5,3) * t424 + Ifges(5,4) * t458 + Ifges(5,2) * t457 + Ifges(5,6) * t511 - pkin(4) * t401 - t489 * t463 + t522 * t465 - t571;
t386 = mrSges(5,2) * t438 - mrSges(5,3) * t423 + Ifges(5,1) * t458 + Ifges(5,4) * t457 + Ifges(5,5) * t511 - pkin(10) * t401 - t533 * t399 + t570 * t400 + t488 * t463 - t522 * t464;
t384 = m(3) * t485 - mrSges(3,2) * t525 + mrSges(3,3) * t519 - t512 * t526 + t517 * t549 + t547;
t383 = mrSges(4,2) * t472 - mrSges(4,3) * t432 + Ifges(4,1) * t495 + Ifges(4,4) * t494 - Ifges(4,5) * t519 - pkin(9) * t393 + t386 * t537 - t387 * t534 + t479 * t507 + t480 * t549;
t382 = -mrSges(4,1) * t472 + mrSges(4,3) * t433 + Ifges(4,4) * t495 + Ifges(4,2) * t494 - Ifges(4,6) * t519 - pkin(3) * t543 + pkin(9) * t546 + t534 * t386 + t537 * t387 - t508 * t479 - t481 * t549;
t381 = Ifges(3,5) * t518 + Ifges(3,6) * t519 + Ifges(3,3) * t525 + mrSges(3,1) * t484 - mrSges(3,2) * t485 + t529 * t383 + t531 * t382 - pkin(2) * t398 + qJ(3) * t547 + (t498 * t535 - t499 * t538) * t554;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t548 - mrSges(2,2) * t544 + (mrSges(3,2) * t500 - mrSges(3,3) * t484 + Ifges(3,1) * t518 + Ifges(3,4) * t519 + Ifges(3,5) * t525 - qJ(3) * t385 - t382 * t529 + t383 * t531 + t497 * t549 - t498 * t526) * t561 + (Ifges(3,6) * t525 + t526 * t499 + Ifges(3,4) * t518 + t507 * t481 - t508 * t480 - Ifges(4,6) * t494 - Ifges(4,5) * t495 - mrSges(3,1) * t500 + mrSges(3,3) * t485 - mrSges(4,1) * t432 + mrSges(4,2) * t433 - pkin(3) * t393 + (Ifges(3,2) + Ifges(4,3)) * t519 - t497 * t550 - t541 - pkin(2) * t385) * t560 + t532 * t381 + pkin(1) * ((t384 * t535 + t394 * t538) * t532 + (-m(3) * t500 + t519 * mrSges(3,1) - t518 * mrSges(3,2) + (-t512 * t535 + t513 * t538) * t554 - t385) * t530) + (t384 * t538 - t394 * t535) * t569; t381; t398; t541; t571; t412;];
tauJ  = t1;
