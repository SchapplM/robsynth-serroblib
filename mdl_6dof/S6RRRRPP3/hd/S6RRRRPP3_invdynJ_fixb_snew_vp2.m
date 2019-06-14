% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:08:13
% EndTime: 2019-05-07 18:08:19
% DurationCPUTime: 2.88s
% Computational Cost: add. (25803->305), mult. (51657->359), div. (0->0), fcn. (35667->8), ass. (0->121)
t557 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t541 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t540 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t556 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t539 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t555 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t512 = sin(qJ(3));
t513 = sin(qJ(2));
t515 = cos(qJ(3));
t516 = cos(qJ(2));
t492 = (t512 * t516 + t513 * t515) * qJD(1);
t542 = qJD(1) * qJD(2);
t497 = qJDD(1) * t513 + t516 * t542;
t498 = qJDD(1) * t516 - t513 * t542;
t461 = -qJD(3) * t492 - t497 * t512 + t498 * t515;
t543 = qJD(1) * t516;
t544 = qJD(1) * t513;
t491 = -t512 * t544 + t515 * t543;
t462 = qJD(3) * t491 + t497 * t515 + t498 * t512;
t501 = qJD(2) * pkin(2) - pkin(8) * t544;
t510 = t516 ^ 2;
t518 = qJD(1) ^ 2;
t514 = sin(qJ(1));
t517 = cos(qJ(1));
t533 = t514 * g(1) - t517 * g(2);
t527 = -qJDD(1) * pkin(1) - t533;
t463 = -t498 * pkin(2) + t501 * t544 + (-pkin(8) * t510 - pkin(7)) * t518 + t527;
t509 = qJD(2) + qJD(3);
t411 = (-t491 * t509 - t462) * pkin(9) + (t492 * t509 - t461) * pkin(3) + t463;
t529 = -g(1) * t517 - g(2) * t514;
t494 = -pkin(1) * t518 + qJDD(1) * pkin(7) + t529;
t546 = t513 * t494;
t549 = pkin(2) * t518;
t452 = qJDD(2) * pkin(2) - t497 * pkin(8) - t546 + (pkin(8) * t542 + t513 * t549 - g(3)) * t516;
t482 = -g(3) * t513 + t516 * t494;
t453 = pkin(8) * t498 - qJD(2) * t501 - t510 * t549 + t482;
t418 = t512 * t452 + t515 * t453;
t476 = -pkin(3) * t491 - pkin(9) * t492;
t507 = t509 ^ 2;
t508 = qJDD(2) + qJDD(3);
t415 = -pkin(3) * t507 + pkin(9) * t508 + t476 * t491 + t418;
t511 = sin(qJ(4));
t551 = cos(qJ(4));
t408 = t551 * t411 - t511 * t415;
t479 = t511 * t492 - t551 * t509;
t480 = t551 * t492 + t511 * t509;
t448 = pkin(4) * t479 - qJ(5) * t480;
t460 = qJDD(4) - t461;
t487 = qJD(4) - t491;
t486 = t487 ^ 2;
t406 = -t460 * pkin(4) - t486 * qJ(5) + t480 * t448 + qJDD(5) - t408;
t425 = -t479 * qJD(4) + t551 * t462 + t511 * t508;
t450 = -mrSges(6,2) * t479 - mrSges(6,3) * t480;
t554 = -m(6) * t406 - t425 * mrSges(6,1) - t480 * t450;
t447 = -mrSges(7,2) * t480 + mrSges(7,3) * t479;
t547 = t479 * t487;
t400 = -0.2e1 * qJD(6) * t487 + (t479 * t480 - t460) * qJ(6) + (t425 + t547) * pkin(5) + t406;
t467 = -mrSges(7,1) * t479 + mrSges(7,2) * t487;
t530 = -m(7) * t400 + t460 * mrSges(7,3) + t487 * t467;
t398 = t425 * mrSges(7,1) + t480 * t447 - t530;
t466 = mrSges(6,1) * t479 - mrSges(6,3) * t487;
t395 = t460 * mrSges(6,2) + t487 * t466 + t398 - t554;
t424 = t480 * qJD(4) + t511 * t462 - t551 * t508;
t464 = pkin(5) * t480 - qJ(6) * t487;
t478 = t479 ^ 2;
t409 = t511 * t411 + t551 * t415;
t524 = -t486 * pkin(4) + t460 * qJ(5) - t479 * t448 + t409;
t402 = -t424 * pkin(5) - t478 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t464) * t487 + t524;
t552 = -2 * qJD(5);
t405 = t487 * t552 - t524;
t468 = mrSges(6,1) * t480 + mrSges(6,2) * t487;
t465 = mrSges(7,1) * t480 - mrSges(7,3) * t487;
t537 = m(7) * t402 + t460 * mrSges(7,2) + t487 * t465;
t526 = -m(6) * t405 + t460 * mrSges(6,3) + t487 * t468 + t537;
t534 = -t541 * t479 + t480 * t557 + t540 * t487;
t535 = t479 * t556 + t480 * t541 + t487 * t539;
t545 = -t447 - t450;
t553 = -t539 * t424 + t540 * t425 + t555 * t460 + t534 * t479 + t535 * t480 + mrSges(5,1) * t408 - mrSges(5,2) * t409 + mrSges(6,2) * t406 + mrSges(7,2) * t402 - mrSges(6,3) * t405 - mrSges(7,3) * t400 - pkin(4) * t395 + qJ(5) * (t545 * t479 + (-mrSges(6,1) - mrSges(7,1)) * t424 + t526) - qJ(6) * t398;
t548 = -mrSges(7,1) - mrSges(5,3);
t475 = -mrSges(4,1) * t491 + mrSges(4,2) * t492;
t484 = mrSges(4,1) * t509 - mrSges(4,3) * t492;
t449 = mrSges(5,1) * t479 + mrSges(5,2) * t480;
t469 = -mrSges(5,2) * t487 - mrSges(5,3) * t479;
t390 = m(5) * t408 + (-t466 + t469) * t487 + (-t447 - t449) * t480 + (mrSges(5,1) - mrSges(6,2)) * t460 + t548 * t425 + t530 + t554;
t470 = mrSges(5,1) * t487 - mrSges(5,3) * t480;
t393 = m(5) * t409 - t460 * mrSges(5,2) - t487 * t470 + (-t449 + t545) * t479 + (-mrSges(6,1) + t548) * t424 + t526;
t531 = -t390 * t511 + t551 * t393;
t384 = m(4) * t418 - mrSges(4,2) * t508 + mrSges(4,3) * t461 + t475 * t491 - t484 * t509 + t531;
t417 = t515 * t452 - t512 * t453;
t483 = -mrSges(4,2) * t509 + mrSges(4,3) * t491;
t414 = -t508 * pkin(3) - t507 * pkin(9) + t492 * t476 - t417;
t521 = (-t425 + t547) * qJ(5) + t414 + (t487 * pkin(4) + t552) * t480;
t407 = t424 * pkin(4) + t521;
t404 = -t478 * pkin(5) + 0.2e1 * qJD(6) * t479 - t480 * t464 + (pkin(4) + qJ(6)) * t424 + t521;
t528 = m(7) * t404 - t425 * mrSges(7,2) + t424 * mrSges(7,3) - t480 * t465 + t479 * t467;
t525 = -m(6) * t407 + t424 * mrSges(6,2) + t479 * t466 - t528;
t520 = -m(5) * t414 - t424 * mrSges(5,1) - t479 * t469 + (t468 - t470) * t480 + (-mrSges(5,2) + mrSges(6,3)) * t425 + t525;
t388 = m(4) * t417 + t508 * mrSges(4,1) - t462 * mrSges(4,3) - t492 * t475 + t509 * t483 + t520;
t379 = t512 * t384 + t515 * t388;
t386 = t551 * t390 + t511 * t393;
t536 = t479 * t539 - t480 * t540 - t487 * t555;
t532 = t515 * t384 - t388 * t512;
t397 = -t425 * mrSges(6,3) - t480 * t468 - t525;
t399 = -t424 * mrSges(7,1) - t479 * t447 + t537;
t378 = -mrSges(5,1) * t414 - mrSges(6,1) * t405 + mrSges(7,1) * t402 + mrSges(6,2) * t407 + mrSges(5,3) * t409 - mrSges(7,3) * t404 - pkin(4) * t397 + pkin(5) * t399 - qJ(6) * t528 + t424 * t556 + t541 * t425 + t539 * t460 + t536 * t480 + t534 * t487;
t381 = mrSges(6,1) * t406 + mrSges(7,1) * t400 + mrSges(5,2) * t414 - mrSges(7,2) * t404 - mrSges(5,3) * t408 - mrSges(6,3) * t407 + pkin(5) * t398 - qJ(5) * t397 - t541 * t424 + t425 * t557 + t540 * t460 + t536 * t479 - t535 * t487;
t473 = Ifges(4,4) * t492 + Ifges(4,2) * t491 + Ifges(4,6) * t509;
t474 = Ifges(4,1) * t492 + Ifges(4,4) * t491 + Ifges(4,5) * t509;
t523 = mrSges(4,1) * t417 - mrSges(4,2) * t418 + Ifges(4,5) * t462 + Ifges(4,6) * t461 + Ifges(4,3) * t508 + pkin(3) * t520 + pkin(9) * t531 + t551 * t378 + t511 * t381 + t492 * t473 - t491 * t474;
t522 = m(4) * t463 - t461 * mrSges(4,1) + t462 * mrSges(4,2) - t491 * t483 + t492 * t484 + t386;
t500 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t543;
t499 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t544;
t496 = (-mrSges(3,1) * t516 + mrSges(3,2) * t513) * qJD(1);
t493 = -t518 * pkin(7) + t527;
t490 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t513 + Ifges(3,4) * t516) * qJD(1);
t489 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t513 + Ifges(3,2) * t516) * qJD(1);
t481 = -t516 * g(3) - t546;
t472 = Ifges(4,5) * t492 + Ifges(4,6) * t491 + Ifges(4,3) * t509;
t376 = -mrSges(4,1) * t463 + mrSges(4,3) * t418 + Ifges(4,4) * t462 + Ifges(4,2) * t461 + Ifges(4,6) * t508 - pkin(3) * t386 - t492 * t472 + t509 * t474 - t553;
t375 = mrSges(4,2) * t463 - mrSges(4,3) * t417 + Ifges(4,1) * t462 + Ifges(4,4) * t461 + Ifges(4,5) * t508 - pkin(9) * t386 - t511 * t378 + t551 * t381 + t491 * t472 - t509 * t473;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t533 - mrSges(2,2) * t529 + t513 * (mrSges(3,2) * t493 - mrSges(3,3) * t481 + Ifges(3,1) * t497 + Ifges(3,4) * t498 + Ifges(3,5) * qJDD(2) - pkin(8) * t379 - qJD(2) * t489 + t515 * t375 - t512 * t376) + t516 * (-mrSges(3,1) * t493 + mrSges(3,3) * t482 + Ifges(3,4) * t497 + Ifges(3,2) * t498 + Ifges(3,6) * qJDD(2) - pkin(2) * t522 + pkin(8) * t532 + qJD(2) * t490 + t512 * t375 + t515 * t376) + pkin(1) * (-m(3) * t493 + t498 * mrSges(3,1) - t497 * mrSges(3,2) + (-t499 * t513 + t500 * t516) * qJD(1) - t522) + pkin(7) * (t516 * (m(3) * t482 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t498 - qJD(2) * t499 + t496 * t543 + t532) - t513 * (m(3) * t481 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t497 + qJD(2) * t500 - t496 * t544 + t379)); t523 + (t513 * t489 - t516 * t490) * qJD(1) + Ifges(3,3) * qJDD(2) + Ifges(3,5) * t497 + Ifges(3,6) * t498 + mrSges(3,1) * t481 - mrSges(3,2) * t482 + pkin(2) * t379; t523; t553; t395; t399;];
tauJ  = t1;
