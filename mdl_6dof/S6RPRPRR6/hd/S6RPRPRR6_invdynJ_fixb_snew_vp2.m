% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 19:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:02:22
% EndTime: 2019-05-05 19:02:30
% DurationCPUTime: 7.39s
% Computational Cost: add. (106295->323), mult. (262602->413), div. (0->0), fcn. (202285->12), ass. (0->136)
t532 = qJD(1) ^ 2;
t558 = cos(qJ(3));
t557 = pkin(2) * t532;
t556 = pkin(7) * qJDD(1);
t527 = sin(qJ(1));
t530 = cos(qJ(1));
t542 = -g(1) * t530 - g(2) * t527;
t510 = -pkin(1) * t532 + qJDD(1) * qJ(2) + t542;
t521 = sin(pkin(10));
t523 = cos(pkin(10));
t550 = qJD(1) * qJD(2);
t547 = -t523 * g(3) - 0.2e1 * t521 * t550;
t476 = (t523 * t557 - t510 - t556) * t521 + t547;
t496 = -g(3) * t521 + (t510 + 0.2e1 * t550) * t523;
t519 = t523 ^ 2;
t483 = -t519 * t557 + t523 * t556 + t496;
t526 = sin(qJ(3));
t459 = t526 * t476 + t558 * t483;
t549 = t523 * t558;
t553 = qJD(1) * t521;
t508 = -qJD(1) * t549 + t526 * t553;
t537 = t558 * t521 + t523 * t526;
t509 = t537 * qJD(1);
t490 = mrSges(4,1) * t508 + mrSges(4,2) * t509;
t552 = qJD(3) * t509;
t493 = t552 + (t521 * t526 - t549) * qJDD(1);
t503 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t509;
t489 = pkin(3) * t508 - qJ(4) * t509;
t531 = qJD(3) ^ 2;
t445 = -pkin(3) * t531 + qJDD(3) * qJ(4) - t489 * t508 + t459;
t548 = t527 * g(1) - t530 * g(2);
t541 = qJDD(2) - t548;
t554 = -t521 ^ 2 - t519;
t492 = (-pkin(2) * t523 - pkin(1)) * qJDD(1) + (t554 * pkin(7) - qJ(2)) * t532 + t541;
t551 = t508 * qJD(3);
t494 = t537 * qJDD(1) - t551;
t450 = (-t494 + t551) * qJ(4) + (t493 + t552) * pkin(3) + t492;
t520 = sin(pkin(11));
t522 = cos(pkin(11));
t501 = qJD(3) * t520 + t509 * t522;
t431 = -0.2e1 * qJD(4) * t501 - t520 * t445 + t522 * t450;
t482 = qJDD(3) * t520 + t494 * t522;
t500 = qJD(3) * t522 - t509 * t520;
t421 = (t500 * t508 - t482) * pkin(8) + (t500 * t501 + t493) * pkin(4) + t431;
t432 = 0.2e1 * qJD(4) * t500 + t522 * t445 + t520 * t450;
t480 = pkin(4) * t508 - pkin(8) * t501;
t481 = qJDD(3) * t522 - t494 * t520;
t499 = t500 ^ 2;
t423 = -pkin(4) * t499 + pkin(8) * t481 - t480 * t508 + t432;
t525 = sin(qJ(5));
t529 = cos(qJ(5));
t415 = t529 * t421 - t525 * t423;
t469 = t500 * t529 - t501 * t525;
t444 = qJD(5) * t469 + t481 * t525 + t482 * t529;
t470 = t500 * t525 + t501 * t529;
t491 = qJDD(5) + t493;
t506 = qJD(5) + t508;
t413 = (t469 * t506 - t444) * pkin(9) + (t469 * t470 + t491) * pkin(5) + t415;
t416 = t525 * t421 + t529 * t423;
t443 = -qJD(5) * t470 + t481 * t529 - t482 * t525;
t462 = pkin(5) * t506 - pkin(9) * t470;
t468 = t469 ^ 2;
t414 = -pkin(5) * t468 + pkin(9) * t443 - t462 * t506 + t416;
t524 = sin(qJ(6));
t528 = cos(qJ(6));
t411 = t413 * t528 - t414 * t524;
t455 = t469 * t528 - t470 * t524;
t429 = qJD(6) * t455 + t443 * t524 + t444 * t528;
t456 = t469 * t524 + t470 * t528;
t439 = -mrSges(7,1) * t455 + mrSges(7,2) * t456;
t505 = qJD(6) + t506;
t448 = -mrSges(7,2) * t505 + mrSges(7,3) * t455;
t488 = qJDD(6) + t491;
t405 = m(7) * t411 + mrSges(7,1) * t488 - mrSges(7,3) * t429 - t439 * t456 + t448 * t505;
t412 = t413 * t524 + t414 * t528;
t428 = -qJD(6) * t456 + t443 * t528 - t444 * t524;
t449 = mrSges(7,1) * t505 - mrSges(7,3) * t456;
t406 = m(7) * t412 - mrSges(7,2) * t488 + mrSges(7,3) * t428 + t439 * t455 - t449 * t505;
t399 = t528 * t405 + t524 * t406;
t457 = -mrSges(6,1) * t469 + mrSges(6,2) * t470;
t460 = -mrSges(6,2) * t506 + mrSges(6,3) * t469;
t397 = m(6) * t415 + mrSges(6,1) * t491 - mrSges(6,3) * t444 - t457 * t470 + t460 * t506 + t399;
t461 = mrSges(6,1) * t506 - mrSges(6,3) * t470;
t543 = -t405 * t524 + t528 * t406;
t398 = m(6) * t416 - mrSges(6,2) * t491 + mrSges(6,3) * t443 + t457 * t469 - t461 * t506 + t543;
t393 = t529 * t397 + t525 * t398;
t471 = -mrSges(5,1) * t500 + mrSges(5,2) * t501;
t478 = -mrSges(5,2) * t508 + mrSges(5,3) * t500;
t391 = m(5) * t431 + mrSges(5,1) * t493 - mrSges(5,3) * t482 - t471 * t501 + t478 * t508 + t393;
t479 = mrSges(5,1) * t508 - mrSges(5,3) * t501;
t544 = -t397 * t525 + t529 * t398;
t392 = m(5) * t432 - mrSges(5,2) * t493 + mrSges(5,3) * t481 + t471 * t500 - t479 * t508 + t544;
t545 = -t391 * t520 + t522 * t392;
t384 = m(4) * t459 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t493 - qJD(3) * t503 - t490 * t508 + t545;
t458 = t558 * t476 - t526 * t483;
t442 = -qJDD(3) * pkin(3) - t531 * qJ(4) + t509 * t489 + qJDD(4) - t458;
t433 = -t481 * pkin(4) - t499 * pkin(8) + t501 * t480 + t442;
t418 = -t443 * pkin(5) - t468 * pkin(9) + t470 * t462 + t433;
t539 = m(7) * t418 - t428 * mrSges(7,1) + t429 * mrSges(7,2) - t455 * t448 + t456 * t449;
t534 = m(6) * t433 - t443 * mrSges(6,1) + mrSges(6,2) * t444 - t469 * t460 + t461 * t470 + t539;
t409 = m(5) * t442 - t481 * mrSges(5,1) + mrSges(5,2) * t482 - t500 * t478 + t479 * t501 + t534;
t502 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t508;
t408 = m(4) * t458 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t494 + qJD(3) * t502 - t490 * t509 - t409;
t555 = t526 * t384 + t558 * t408;
t385 = t522 * t391 + t520 * t392;
t546 = t558 * t384 - t526 * t408;
t540 = -mrSges(3,1) * t523 + mrSges(3,2) * t521;
t538 = mrSges(3,3) * qJDD(1) + t532 * t540;
t435 = Ifges(7,4) * t456 + Ifges(7,2) * t455 + Ifges(7,6) * t505;
t436 = Ifges(7,1) * t456 + Ifges(7,4) * t455 + Ifges(7,5) * t505;
t536 = -mrSges(7,1) * t411 + mrSges(7,2) * t412 - Ifges(7,5) * t429 - Ifges(7,6) * t428 - Ifges(7,3) * t488 - t456 * t435 + t455 * t436;
t535 = m(4) * t492 + t493 * mrSges(4,1) + t494 * mrSges(4,2) + t508 * t502 + t509 * t503 + t385;
t452 = Ifges(6,4) * t470 + Ifges(6,2) * t469 + Ifges(6,6) * t506;
t453 = Ifges(6,1) * t470 + Ifges(6,4) * t469 + Ifges(6,5) * t506;
t533 = mrSges(6,1) * t415 - mrSges(6,2) * t416 + Ifges(6,5) * t444 + Ifges(6,6) * t443 + Ifges(6,3) * t491 + pkin(5) * t399 + t470 * t452 - t469 * t453 - t536;
t512 = (Ifges(3,5) * t521 + Ifges(3,6) * t523) * qJD(1);
t507 = -qJDD(1) * pkin(1) - t532 * qJ(2) + t541;
t495 = -t521 * t510 + t547;
t486 = Ifges(4,1) * t509 - Ifges(4,4) * t508 + Ifges(4,5) * qJD(3);
t485 = Ifges(4,4) * t509 - Ifges(4,2) * t508 + Ifges(4,6) * qJD(3);
t484 = Ifges(4,5) * t509 - Ifges(4,6) * t508 + Ifges(4,3) * qJD(3);
t465 = Ifges(5,1) * t501 + Ifges(5,4) * t500 + Ifges(5,5) * t508;
t464 = Ifges(5,4) * t501 + Ifges(5,2) * t500 + Ifges(5,6) * t508;
t463 = Ifges(5,5) * t501 + Ifges(5,6) * t500 + Ifges(5,3) * t508;
t451 = Ifges(6,5) * t470 + Ifges(6,6) * t469 + Ifges(6,3) * t506;
t434 = Ifges(7,5) * t456 + Ifges(7,6) * t455 + Ifges(7,3) * t505;
t401 = mrSges(7,2) * t418 - mrSges(7,3) * t411 + Ifges(7,1) * t429 + Ifges(7,4) * t428 + Ifges(7,5) * t488 + t434 * t455 - t435 * t505;
t400 = -mrSges(7,1) * t418 + mrSges(7,3) * t412 + Ifges(7,4) * t429 + Ifges(7,2) * t428 + Ifges(7,6) * t488 - t434 * t456 + t436 * t505;
t387 = mrSges(6,2) * t433 - mrSges(6,3) * t415 + Ifges(6,1) * t444 + Ifges(6,4) * t443 + Ifges(6,5) * t491 - pkin(9) * t399 - t400 * t524 + t401 * t528 + t451 * t469 - t452 * t506;
t386 = -mrSges(6,1) * t433 + mrSges(6,3) * t416 + Ifges(6,4) * t444 + Ifges(6,2) * t443 + Ifges(6,6) * t491 - pkin(5) * t539 + pkin(9) * t543 + t528 * t400 + t524 * t401 - t470 * t451 + t506 * t453;
t381 = t554 * t532 * mrSges(3,3) + m(3) * t507 + t540 * qJDD(1) + t535;
t380 = mrSges(5,2) * t442 - mrSges(5,3) * t431 + Ifges(5,1) * t482 + Ifges(5,4) * t481 + Ifges(5,5) * t493 - pkin(8) * t393 - t386 * t525 + t387 * t529 + t463 * t500 - t464 * t508;
t379 = -mrSges(5,1) * t442 + mrSges(5,3) * t432 + Ifges(5,4) * t482 + Ifges(5,2) * t481 + Ifges(5,6) * t493 - pkin(4) * t534 + pkin(8) * t544 + t529 * t386 + t525 * t387 - t501 * t463 + t508 * t465;
t378 = (-Ifges(5,3) - Ifges(4,2)) * t493 - t533 - t509 * t484 + t500 * t465 - t501 * t464 - mrSges(4,1) * t492 + Ifges(4,4) * t494 - Ifges(5,6) * t481 - Ifges(5,5) * t482 + qJD(3) * t486 + mrSges(4,3) * t459 - mrSges(5,1) * t431 + mrSges(5,2) * t432 + Ifges(4,6) * qJDD(3) - pkin(4) * t393 - pkin(3) * t385;
t377 = mrSges(4,2) * t492 - mrSges(4,3) * t458 + Ifges(4,1) * t494 - Ifges(4,4) * t493 + Ifges(4,5) * qJDD(3) - qJ(4) * t385 - qJD(3) * t485 - t379 * t520 + t380 * t522 - t484 * t508;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t548 - mrSges(2,2) * t542 + t521 * (t523 * qJD(1) * t512 + mrSges(3,2) * t507 - mrSges(3,3) * t495 + t558 * t377 - t526 * t378 - pkin(7) * t555 + (Ifges(3,1) * t521 + Ifges(3,4) * t523) * qJDD(1)) + t523 * (-t512 * t553 - mrSges(3,1) * t507 + mrSges(3,3) * t496 + t526 * t377 + t558 * t378 - pkin(2) * t535 + pkin(7) * t546 + (Ifges(3,4) * t521 + Ifges(3,2) * t523) * qJDD(1)) - pkin(1) * t381 + qJ(2) * ((m(3) * t496 + t538 * t523 + t546) * t523 + (-m(3) * t495 + t538 * t521 - t555) * t521); t381; mrSges(4,1) * t458 - mrSges(4,2) * t459 + Ifges(4,5) * t494 - Ifges(4,6) * t493 + Ifges(4,3) * qJDD(3) - pkin(3) * t409 + qJ(4) * t545 + t522 * t379 + t520 * t380 + t509 * t485 + t508 * t486; t409; t533; -t536;];
tauJ  = t1;
