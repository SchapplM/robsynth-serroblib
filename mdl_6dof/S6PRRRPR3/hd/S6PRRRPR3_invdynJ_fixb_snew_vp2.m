% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 07:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:28:02
% EndTime: 2019-05-05 07:28:06
% DurationCPUTime: 2.39s
% Computational Cost: add. (15320->271), mult. (30665->331), div. (0->0), fcn. (21263->12), ass. (0->117)
t528 = Ifges(5,4) + Ifges(6,6);
t537 = -Ifges(5,2) - Ifges(6,3);
t534 = Ifges(5,6) - Ifges(6,5);
t536 = Ifges(5,1) + Ifges(6,2);
t535 = Ifges(5,5) - Ifges(6,4);
t533 = Ifges(5,3) + Ifges(6,1);
t498 = sin(qJ(4));
t502 = cos(qJ(3));
t520 = qJD(2) * t502;
t499 = sin(qJ(3));
t521 = qJD(2) * t499;
t529 = cos(qJ(4));
t467 = t498 * t521 - t520 * t529;
t468 = (t498 * t502 + t499 * t529) * qJD(2);
t490 = qJD(3) + qJD(4);
t532 = t537 * t467 + t528 * t468 + t534 * t490;
t493 = sin(pkin(11));
t495 = cos(pkin(11));
t478 = g(1) * t493 - g(2) * t495;
t492 = -g(3) + qJDD(1);
t494 = sin(pkin(6));
t496 = cos(pkin(6));
t531 = t478 * t496 + t492 * t494;
t479 = -g(1) * t495 - g(2) * t493;
t500 = sin(qJ(2));
t503 = cos(qJ(2));
t440 = -t500 * t479 + t503 * t531;
t530 = -2 * qJD(5);
t527 = t467 * t490;
t441 = t503 * t479 + t500 * t531;
t504 = qJD(2) ^ 2;
t433 = -pkin(2) * t504 + qJDD(2) * pkin(8) + t441;
t458 = -t478 * t494 + t492 * t496;
t410 = -t499 * t433 + t502 * t458;
t519 = qJD(2) * qJD(3);
t518 = t502 * t519;
t476 = qJDD(2) * t499 + t518;
t397 = (-t476 + t518) * pkin(9) + (t499 * t502 * t504 + qJDD(3)) * pkin(3) + t410;
t411 = t502 * t433 + t499 * t458;
t477 = qJDD(2) * t502 - t499 * t519;
t483 = qJD(3) * pkin(3) - pkin(9) * t521;
t491 = t502 ^ 2;
t398 = -pkin(3) * t491 * t504 + pkin(9) * t477 - qJD(3) * t483 + t411;
t392 = t397 * t529 - t498 * t398;
t429 = -t467 * qJD(4) + t476 * t529 + t498 * t477;
t446 = mrSges(5,1) * t467 + mrSges(5,2) * t468;
t453 = -mrSges(5,2) * t490 - mrSges(5,3) * t467;
t455 = mrSges(6,1) * t467 - mrSges(6,3) * t490;
t489 = qJDD(3) + qJDD(4);
t445 = pkin(4) * t467 - qJ(5) * t468;
t488 = t490 ^ 2;
t388 = -t489 * pkin(4) - t488 * qJ(5) + t468 * t445 + qJDD(5) - t392;
t382 = (t467 * t468 - t489) * pkin(10) + (t429 + t527) * pkin(5) + t388;
t428 = t468 * qJD(4) + t498 * t476 - t477 * t529;
t457 = pkin(5) * t468 - pkin(10) * t490;
t463 = t467 ^ 2;
t510 = -qJDD(2) * pkin(2) - t440;
t405 = -t477 * pkin(3) + t483 * t521 + (-pkin(9) * t491 - pkin(8)) * t504 + t510;
t505 = (-t429 + t527) * qJ(5) + t405 + (pkin(4) * t490 + t530) * t468;
t385 = -t463 * pkin(5) - t468 * t457 + (pkin(4) + pkin(10)) * t428 + t505;
t497 = sin(qJ(6));
t501 = cos(qJ(6));
t380 = t382 * t501 - t385 * t497;
t449 = t467 * t501 - t490 * t497;
t404 = qJD(6) * t449 + t428 * t497 + t489 * t501;
t450 = t467 * t497 + t490 * t501;
t412 = -mrSges(7,1) * t449 + mrSges(7,2) * t450;
t426 = qJDD(6) + t429;
t461 = qJD(6) + t468;
t430 = -mrSges(7,2) * t461 + mrSges(7,3) * t449;
t377 = m(7) * t380 + mrSges(7,1) * t426 - mrSges(7,3) * t404 - t412 * t450 + t430 * t461;
t381 = t382 * t497 + t385 * t501;
t403 = -qJD(6) * t450 + t428 * t501 - t489 * t497;
t431 = mrSges(7,1) * t461 - mrSges(7,3) * t450;
t378 = m(7) * t381 - mrSges(7,2) * t426 + mrSges(7,3) * t403 + t412 * t449 - t431 * t461;
t367 = t501 * t377 + t497 * t378;
t447 = -mrSges(6,2) * t467 - mrSges(6,3) * t468;
t513 = -m(6) * t388 - t429 * mrSges(6,1) - t468 * t447 - t367;
t363 = m(5) * t392 - t429 * mrSges(5,3) - t468 * t446 + (t453 - t455) * t490 + (mrSges(5,1) - mrSges(6,2)) * t489 + t513;
t393 = t498 * t397 + t529 * t398;
t454 = mrSges(5,1) * t490 - mrSges(5,3) * t468;
t512 = -t488 * pkin(4) + t489 * qJ(5) - t467 * t445 + t393;
t386 = t490 * t530 - t512;
t456 = mrSges(6,1) * t468 + mrSges(6,2) * t490;
t384 = -t428 * pkin(5) - t463 * pkin(10) + ((2 * qJD(5)) + t457) * t490 + t512;
t514 = -m(7) * t384 + t403 * mrSges(7,1) - t404 * mrSges(7,2) + t449 * t430 - t450 * t431;
t509 = -m(6) * t386 + t489 * mrSges(6,3) + t490 * t456 - t514;
t373 = m(5) * t393 - t489 * mrSges(5,2) - t490 * t454 + (-t446 - t447) * t467 + (-mrSges(5,3) - mrSges(6,1)) * t428 + t509;
t361 = t529 * t363 + t498 * t373;
t524 = -t497 * t377 + t501 * t378;
t523 = t467 * t534 - t468 * t535 - t490 * t533;
t522 = -t528 * t467 + t468 * t536 + t535 * t490;
t475 = (-mrSges(4,1) * t502 + mrSges(4,2) * t499) * qJD(2);
t481 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t520;
t359 = m(4) * t410 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t476 + qJD(3) * t481 - t475 * t521 + t361;
t480 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t521;
t516 = -t363 * t498 + t529 * t373;
t360 = m(4) * t411 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t477 - qJD(3) * t480 + t475 * t520 + t516;
t517 = -t359 * t499 + t502 * t360;
t390 = t428 * pkin(4) + t505;
t365 = m(6) * t390 - t428 * mrSges(6,2) - t429 * mrSges(6,3) - t467 * t455 - t468 * t456 + t524;
t407 = Ifges(7,4) * t450 + Ifges(7,2) * t449 + Ifges(7,6) * t461;
t408 = Ifges(7,1) * t450 + Ifges(7,4) * t449 + Ifges(7,5) * t461;
t511 = mrSges(7,1) * t380 - mrSges(7,2) * t381 + Ifges(7,5) * t404 + Ifges(7,6) * t403 + Ifges(7,3) * t426 + t450 * t407 - t449 * t408;
t508 = m(5) * t405 + t428 * mrSges(5,1) + mrSges(5,2) * t429 + t467 * t453 + t454 * t468 + t365;
t432 = -t504 * pkin(8) + t510;
t507 = -m(4) * t432 + t477 * mrSges(4,1) - mrSges(4,2) * t476 - t480 * t521 + t481 * t520 - t508;
t366 = t489 * mrSges(6,2) + t490 * t455 - t513;
t406 = Ifges(7,5) * t450 + Ifges(7,6) * t449 + Ifges(7,3) * t461;
t369 = -mrSges(7,1) * t384 + mrSges(7,3) * t381 + Ifges(7,4) * t404 + Ifges(7,2) * t403 + Ifges(7,6) * t426 - t406 * t450 + t408 * t461;
t370 = mrSges(7,2) * t384 - mrSges(7,3) * t380 + Ifges(7,1) * t404 + Ifges(7,4) * t403 + Ifges(7,5) * t426 + t406 * t449 - t407 * t461;
t506 = -mrSges(5,2) * t393 - mrSges(6,3) * t386 - pkin(4) * t366 - pkin(10) * t367 - t497 * t369 + t501 * t370 + t522 * t467 + qJ(5) * (-t467 * t447 + t509) + mrSges(6,2) * t388 + mrSges(5,1) * t392 + t533 * t489 + t532 * t468 + t535 * t429 + (-mrSges(6,1) * qJ(5) - t534) * t428;
t466 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t499 + Ifges(4,4) * t502) * qJD(2);
t465 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t499 + Ifges(4,2) * t502) * qJD(2);
t357 = mrSges(6,1) * t388 + mrSges(5,2) * t405 - mrSges(5,3) * t392 - mrSges(6,3) * t390 + pkin(5) * t367 - qJ(5) * t365 - t528 * t428 + t429 * t536 + t523 * t467 + t535 * t489 - t532 * t490 + t511;
t356 = -mrSges(5,1) * t405 - mrSges(6,1) * t386 + mrSges(6,2) * t390 + mrSges(5,3) * t393 - pkin(4) * t365 - pkin(5) * t514 - pkin(10) * t524 - t501 * t369 - t497 * t370 + t537 * t428 + t528 * t429 + t523 * t468 + t534 * t489 + t522 * t490;
t1 = [m(2) * t492 + t496 * (m(3) * t458 + t359 * t502 + t360 * t499) + (t500 * (m(3) * t441 - mrSges(3,1) * t504 - qJDD(2) * mrSges(3,2) + t517) + t503 * (m(3) * t440 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t504 + t507)) * t494; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t440 - mrSges(3,2) * t441 + t499 * (mrSges(4,2) * t432 - mrSges(4,3) * t410 + Ifges(4,1) * t476 + Ifges(4,4) * t477 + Ifges(4,5) * qJDD(3) - pkin(9) * t361 - qJD(3) * t465 - t498 * t356 + t357 * t529) + t502 * (-mrSges(4,1) * t432 + mrSges(4,3) * t411 + Ifges(4,4) * t476 + Ifges(4,2) * t477 + Ifges(4,6) * qJDD(3) - pkin(3) * t508 + pkin(9) * t516 + qJD(3) * t466 + t356 * t529 + t498 * t357) + pkin(2) * t507 + pkin(8) * t517; Ifges(4,3) * qJDD(3) + t506 + Ifges(4,6) * t477 + Ifges(4,5) * t476 - mrSges(4,2) * t411 + mrSges(4,1) * t410 + pkin(3) * t361 + (t465 * t499 - t466 * t502) * qJD(2); t506; t366; t511;];
tauJ  = t1;
