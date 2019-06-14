% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 22:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:49:15
% EndTime: 2019-05-05 22:49:23
% DurationCPUTime: 7.88s
% Computational Cost: add. (113531->323), mult. (273919->414), div. (0->0), fcn. (211276->12), ass. (0->136)
t533 = qJD(1) ^ 2;
t527 = sin(qJ(1));
t531 = cos(qJ(1));
t543 = -t531 * g(1) - t527 * g(2);
t510 = -t533 * pkin(1) + qJDD(1) * qJ(2) + t543;
t521 = sin(pkin(10));
t523 = cos(pkin(10));
t551 = qJD(1) * qJD(2);
t548 = -t523 * g(3) - 0.2e1 * t521 * t551;
t558 = pkin(2) * t533;
t483 = (-pkin(7) * qJDD(1) + t523 * t558 - t510) * t521 + t548;
t497 = -t521 * g(3) + (t510 + 0.2e1 * t551) * t523;
t519 = t523 ^ 2;
t550 = qJDD(1) * t523;
t484 = pkin(7) * t550 - t519 * t558 + t497;
t526 = sin(qJ(3));
t530 = cos(qJ(3));
t458 = t526 * t483 + t530 * t484;
t553 = qJD(1) * t523;
t554 = qJD(1) * t521;
t508 = -t526 * t554 + t530 * t553;
t539 = t521 * t530 + t523 * t526;
t509 = t539 * qJD(1);
t492 = -pkin(3) * t508 - pkin(8) * t509;
t532 = qJD(3) ^ 2;
t441 = -t532 * pkin(3) + qJDD(3) * pkin(8) + t492 * t508 + t458;
t549 = t527 * g(1) - t531 * g(2);
t542 = qJDD(2) - t549;
t555 = -t521 ^ 2 - t519;
t493 = (-pkin(2) * t523 - pkin(1)) * qJDD(1) + (pkin(7) * t555 - qJ(2)) * t533 + t542;
t505 = t509 * qJD(3);
t494 = -t526 * t521 * qJDD(1) + t530 * t550 - t505;
t552 = t508 * qJD(3);
t495 = t539 * qJDD(1) + t552;
t449 = (-t495 - t552) * pkin(8) + (-t494 + t505) * pkin(3) + t493;
t525 = sin(qJ(4));
t529 = cos(qJ(4));
t431 = -t525 * t441 + t529 * t449;
t499 = t529 * qJD(3) - t525 * t509;
t469 = t499 * qJD(4) + t525 * qJDD(3) + t529 * t495;
t491 = qJDD(4) - t494;
t500 = t525 * qJD(3) + t529 * t509;
t506 = qJD(4) - t508;
t420 = (t499 * t506 - t469) * qJ(5) + (t499 * t500 + t491) * pkin(4) + t431;
t432 = t529 * t441 + t525 * t449;
t468 = -t500 * qJD(4) + t529 * qJDD(3) - t525 * t495;
t479 = pkin(4) * t506 - qJ(5) * t500;
t498 = t499 ^ 2;
t422 = -pkin(4) * t498 + qJ(5) * t468 - t479 * t506 + t432;
t520 = sin(pkin(11));
t522 = cos(pkin(11));
t474 = t520 * t499 + t522 * t500;
t414 = -0.2e1 * qJD(5) * t474 + t522 * t420 - t520 * t422;
t444 = t520 * t468 + t522 * t469;
t473 = t522 * t499 - t520 * t500;
t412 = (t473 * t506 - t444) * pkin(9) + (t473 * t474 + t491) * pkin(5) + t414;
t415 = 0.2e1 * qJD(5) * t473 + t520 * t420 + t522 * t422;
t443 = t522 * t468 - t520 * t469;
t461 = pkin(5) * t506 - pkin(9) * t474;
t472 = t473 ^ 2;
t413 = -pkin(5) * t472 + pkin(9) * t443 - t461 * t506 + t415;
t524 = sin(qJ(6));
t528 = cos(qJ(6));
t410 = t528 * t412 - t524 * t413;
t454 = t528 * t473 - t524 * t474;
t428 = t454 * qJD(6) + t524 * t443 + t528 * t444;
t455 = t524 * t473 + t528 * t474;
t438 = -mrSges(7,1) * t454 + mrSges(7,2) * t455;
t504 = qJD(6) + t506;
t447 = -mrSges(7,2) * t504 + mrSges(7,3) * t454;
t489 = qJDD(6) + t491;
t404 = m(7) * t410 + mrSges(7,1) * t489 - t428 * mrSges(7,3) - t438 * t455 + t447 * t504;
t411 = t524 * t412 + t528 * t413;
t427 = -t455 * qJD(6) + t528 * t443 - t524 * t444;
t448 = mrSges(7,1) * t504 - mrSges(7,3) * t455;
t405 = m(7) * t411 - mrSges(7,2) * t489 + t427 * mrSges(7,3) + t438 * t454 - t448 * t504;
t398 = t528 * t404 + t524 * t405;
t456 = -mrSges(6,1) * t473 + mrSges(6,2) * t474;
t459 = -mrSges(6,2) * t506 + mrSges(6,3) * t473;
t396 = m(6) * t414 + mrSges(6,1) * t491 - mrSges(6,3) * t444 - t456 * t474 + t459 * t506 + t398;
t460 = mrSges(6,1) * t506 - mrSges(6,3) * t474;
t544 = -t524 * t404 + t528 * t405;
t397 = m(6) * t415 - t491 * mrSges(6,2) + t443 * mrSges(6,3) + t473 * t456 - t506 * t460 + t544;
t392 = t522 * t396 + t520 * t397;
t451 = Ifges(6,4) * t474 + Ifges(6,2) * t473 + Ifges(6,6) * t506;
t452 = Ifges(6,1) * t474 + Ifges(6,4) * t473 + Ifges(6,5) * t506;
t463 = Ifges(5,4) * t500 + Ifges(5,2) * t499 + Ifges(5,6) * t506;
t464 = Ifges(5,1) * t500 + Ifges(5,4) * t499 + Ifges(5,5) * t506;
t434 = Ifges(7,4) * t455 + Ifges(7,2) * t454 + Ifges(7,6) * t504;
t435 = Ifges(7,1) * t455 + Ifges(7,4) * t454 + Ifges(7,5) * t504;
t537 = -mrSges(7,1) * t410 + mrSges(7,2) * t411 - Ifges(7,5) * t428 - Ifges(7,6) * t427 - Ifges(7,3) * t489 - t455 * t434 + t454 * t435;
t559 = mrSges(5,1) * t431 + mrSges(6,1) * t414 - mrSges(5,2) * t432 - mrSges(6,2) * t415 + Ifges(5,5) * t469 + Ifges(6,5) * t444 + Ifges(5,6) * t468 + Ifges(6,6) * t443 + pkin(4) * t392 + pkin(5) * t398 + t474 * t451 - t473 * t452 + t500 * t463 - t499 * t464 + (Ifges(5,3) + Ifges(6,3)) * t491 - t537;
t490 = -mrSges(4,1) * t508 + mrSges(4,2) * t509;
t502 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t509;
t475 = -mrSges(5,1) * t499 + mrSges(5,2) * t500;
t478 = -mrSges(5,2) * t506 + mrSges(5,3) * t499;
t390 = m(5) * t431 + mrSges(5,1) * t491 - mrSges(5,3) * t469 - t475 * t500 + t478 * t506 + t392;
t480 = mrSges(5,1) * t506 - mrSges(5,3) * t500;
t545 = -t520 * t396 + t522 * t397;
t391 = m(5) * t432 - t491 * mrSges(5,2) + t468 * mrSges(5,3) + t499 * t475 - t506 * t480 + t545;
t546 = -t525 * t390 + t529 * t391;
t383 = m(4) * t458 - qJDD(3) * mrSges(4,2) + t494 * mrSges(4,3) - qJD(3) * t502 + t508 * t490 + t546;
t457 = t530 * t483 - t526 * t484;
t501 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t508;
t440 = -qJDD(3) * pkin(3) - t532 * pkin(8) + t509 * t492 - t457;
t430 = -t468 * pkin(4) - t498 * qJ(5) + t500 * t479 + qJDD(5) + t440;
t417 = -t443 * pkin(5) - t472 * pkin(9) + t474 * t461 + t430;
t540 = m(7) * t417 - t427 * mrSges(7,1) + t428 * mrSges(7,2) - t454 * t447 + t455 * t448;
t408 = m(6) * t430 - t443 * mrSges(6,1) + mrSges(6,2) * t444 - t473 * t459 + t460 * t474 + t540;
t535 = -m(5) * t440 + t468 * mrSges(5,1) - mrSges(5,2) * t469 + t499 * t478 - t480 * t500 - t408;
t407 = m(4) * t457 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t495 + qJD(3) * t501 - t490 * t509 + t535;
t556 = t526 * t383 + t530 * t407;
t384 = t529 * t390 + t525 * t391;
t547 = t530 * t383 - t526 * t407;
t541 = -t523 * mrSges(3,1) + t521 * mrSges(3,2);
t538 = mrSges(3,3) * qJDD(1) + t533 * t541;
t536 = m(4) * t493 - t494 * mrSges(4,1) + t495 * mrSges(4,2) - t508 * t501 + t509 * t502 + t384;
t512 = (Ifges(3,5) * t521 + Ifges(3,6) * t523) * qJD(1);
t507 = -qJDD(1) * pkin(1) - t533 * qJ(2) + t542;
t496 = -t521 * t510 + t548;
t487 = Ifges(4,1) * t509 + Ifges(4,4) * t508 + Ifges(4,5) * qJD(3);
t486 = Ifges(4,4) * t509 + Ifges(4,2) * t508 + Ifges(4,6) * qJD(3);
t485 = Ifges(4,5) * t509 + Ifges(4,6) * t508 + Ifges(4,3) * qJD(3);
t462 = Ifges(5,5) * t500 + Ifges(5,6) * t499 + Ifges(5,3) * t506;
t450 = Ifges(6,5) * t474 + Ifges(6,6) * t473 + Ifges(6,3) * t506;
t433 = Ifges(7,5) * t455 + Ifges(7,6) * t454 + Ifges(7,3) * t504;
t400 = mrSges(7,2) * t417 - mrSges(7,3) * t410 + Ifges(7,1) * t428 + Ifges(7,4) * t427 + Ifges(7,5) * t489 + t433 * t454 - t434 * t504;
t399 = -mrSges(7,1) * t417 + mrSges(7,3) * t411 + Ifges(7,4) * t428 + Ifges(7,2) * t427 + Ifges(7,6) * t489 - t433 * t455 + t435 * t504;
t386 = mrSges(6,2) * t430 - mrSges(6,3) * t414 + Ifges(6,1) * t444 + Ifges(6,4) * t443 + Ifges(6,5) * t491 - pkin(9) * t398 - t524 * t399 + t528 * t400 + t473 * t450 - t506 * t451;
t385 = -mrSges(6,1) * t430 + mrSges(6,3) * t415 + Ifges(6,4) * t444 + Ifges(6,2) * t443 + Ifges(6,6) * t491 - pkin(5) * t540 + pkin(9) * t544 + t528 * t399 + t524 * t400 - t474 * t450 + t506 * t452;
t380 = mrSges(3,3) * t533 * t555 + m(3) * t507 + qJDD(1) * t541 + t536;
t379 = mrSges(5,2) * t440 - mrSges(5,3) * t431 + Ifges(5,1) * t469 + Ifges(5,4) * t468 + Ifges(5,5) * t491 - qJ(5) * t392 - t520 * t385 + t522 * t386 + t499 * t462 - t506 * t463;
t378 = -mrSges(5,1) * t440 + mrSges(5,3) * t432 + Ifges(5,4) * t469 + Ifges(5,2) * t468 + Ifges(5,6) * t491 - pkin(4) * t408 + qJ(5) * t545 + t522 * t385 + t520 * t386 - t500 * t462 + t506 * t464;
t377 = -mrSges(4,1) * t493 + mrSges(4,3) * t458 + Ifges(4,4) * t495 + Ifges(4,2) * t494 + Ifges(4,6) * qJDD(3) - pkin(3) * t384 + qJD(3) * t487 - t509 * t485 - t559;
t376 = mrSges(4,2) * t493 - mrSges(4,3) * t457 + Ifges(4,1) * t495 + Ifges(4,4) * t494 + Ifges(4,5) * qJDD(3) - pkin(8) * t384 - qJD(3) * t486 - t525 * t378 + t529 * t379 + t508 * t485;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t549 - mrSges(2,2) * t543 + t521 * (t512 * t553 + mrSges(3,2) * t507 - mrSges(3,3) * t496 + t530 * t376 - t526 * t377 - pkin(7) * t556 + (Ifges(3,1) * t521 + Ifges(3,4) * t523) * qJDD(1)) + t523 * (-t512 * t554 - mrSges(3,1) * t507 + mrSges(3,3) * t497 + t526 * t376 + t530 * t377 - pkin(2) * t536 + pkin(7) * t547 + (Ifges(3,4) * t521 + Ifges(3,2) * t523) * qJDD(1)) - pkin(1) * t380 + qJ(2) * ((m(3) * t497 + t523 * t538 + t547) * t523 + (-m(3) * t496 + t521 * t538 - t556) * t521); t380; mrSges(4,1) * t457 - mrSges(4,2) * t458 + Ifges(4,5) * t495 + Ifges(4,6) * t494 + Ifges(4,3) * qJDD(3) + pkin(3) * t535 + pkin(8) * t546 + t529 * t378 + t525 * t379 + t509 * t486 - t508 * t487; t559; t408; -t537;];
tauJ  = t1;
