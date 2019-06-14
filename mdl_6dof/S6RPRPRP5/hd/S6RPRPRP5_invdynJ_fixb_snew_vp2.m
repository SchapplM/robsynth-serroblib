% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRP5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:47:33
% EndTime: 2019-05-05 17:47:39
% DurationCPUTime: 4.46s
% Computational Cost: add. (41559->299), mult. (102037->371), div. (0->0), fcn. (76359->10), ass. (0->126)
t557 = Ifges(6,1) + Ifges(7,1);
t549 = Ifges(6,4) - Ifges(7,5);
t548 = Ifges(7,4) + Ifges(6,5);
t556 = Ifges(6,2) + Ifges(7,3);
t555 = Ifges(7,2) + Ifges(6,3);
t546 = Ifges(6,6) - Ifges(7,6);
t518 = qJD(1) ^ 2;
t510 = sin(pkin(9));
t512 = cos(pkin(9));
t514 = sin(qJ(3));
t553 = cos(qJ(3));
t522 = t510 * t553 + t512 * t514;
t498 = t522 * qJD(1);
t509 = sin(pkin(10));
t511 = cos(pkin(10));
t489 = t511 * qJD(3) - t509 * t498;
t490 = t509 * qJD(3) + t511 * t498;
t513 = sin(qJ(5));
t552 = cos(qJ(5));
t456 = -t489 * t552 + t513 * t490;
t533 = t512 * t553;
t538 = qJD(1) * t510;
t497 = -qJD(1) * t533 + t514 * t538;
t537 = t497 * qJD(3);
t482 = qJDD(1) * t522 - t537;
t469 = t511 * qJDD(3) - t509 * t482;
t470 = t509 * qJDD(3) + t511 * t482;
t426 = -t456 * qJD(5) + t513 * t469 + t470 * t552;
t457 = t513 * t489 + t490 * t552;
t441 = t456 * mrSges(7,1) - t457 * mrSges(7,3);
t515 = sin(qJ(1));
t516 = cos(qJ(1));
t526 = -t516 * g(1) - t515 * g(2);
t499 = -t518 * pkin(1) + qJDD(1) * qJ(2) + t526;
t535 = qJD(1) * qJD(2);
t531 = -t512 * g(3) - 0.2e1 * t510 * t535;
t545 = pkin(7) * qJDD(1);
t551 = pkin(2) * t518;
t464 = (t512 * t551 - t499 - t545) * t510 + t531;
t484 = -t510 * g(3) + (t499 + 0.2e1 * t535) * t512;
t508 = t512 ^ 2;
t471 = -t508 * t551 + t512 * t545 + t484;
t446 = t514 * t464 + t553 * t471;
t477 = t497 * pkin(3) - t498 * qJ(4);
t517 = qJD(3) ^ 2;
t427 = -t517 * pkin(3) + qJDD(3) * qJ(4) - t497 * t477 + t446;
t532 = t515 * g(1) - t516 * g(2);
t525 = qJDD(2) - t532;
t539 = -t510 ^ 2 - t508;
t480 = (-pkin(2) * t512 - pkin(1)) * qJDD(1) + (pkin(7) * t539 - qJ(2)) * t518 + t525;
t536 = t498 * qJD(3);
t481 = t536 + (t510 * t514 - t533) * qJDD(1);
t430 = (-t482 + t537) * qJ(4) + (t481 + t536) * pkin(3) + t480;
t417 = -0.2e1 * qJD(4) * t490 - t509 * t427 + t511 * t430;
t414 = (t489 * t497 - t470) * pkin(8) + (t489 * t490 + t481) * pkin(4) + t417;
t418 = 0.2e1 * qJD(4) * t489 + t511 * t427 + t509 * t430;
t468 = t497 * pkin(4) - t490 * pkin(8);
t488 = t489 ^ 2;
t416 = -t488 * pkin(4) + t469 * pkin(8) - t497 * t468 + t418;
t409 = t414 * t552 - t513 * t416;
t440 = t456 * pkin(5) - t457 * qJ(6);
t479 = qJDD(5) + t481;
t495 = qJD(5) + t497;
t494 = t495 ^ 2;
t408 = -t479 * pkin(5) - t494 * qJ(6) + t457 * t440 + qJDD(6) - t409;
t447 = -t456 * mrSges(7,2) + t495 * mrSges(7,3);
t527 = -m(7) * t408 + t479 * mrSges(7,1) + t495 * t447;
t404 = t426 * mrSges(7,2) + t457 * t441 - t527;
t410 = t513 * t414 + t552 * t416;
t407 = -t494 * pkin(5) + t479 * qJ(6) + 0.2e1 * qJD(6) * t495 - t456 * t440 + t410;
t425 = t457 * qJD(5) - t469 * t552 + t513 * t470;
t450 = -t495 * mrSges(7,1) + t457 * mrSges(7,2);
t534 = m(7) * t407 + t479 * mrSges(7,3) + t495 * t450;
t541 = -t549 * t456 + t557 * t457 + t548 * t495;
t542 = t556 * t456 - t549 * t457 - t546 * t495;
t554 = t546 * t425 - t548 * t426 - t541 * t456 + t542 * t457 - t555 * t479 - mrSges(6,1) * t409 + mrSges(7,1) * t408 + mrSges(6,2) * t410 - mrSges(7,3) * t407 + pkin(5) * t404 - qJ(6) * (-t425 * mrSges(7,2) - t456 * t441 + t534);
t550 = -mrSges(6,3) - mrSges(7,2);
t478 = t497 * mrSges(4,1) + t498 * mrSges(4,2);
t492 = qJD(3) * mrSges(4,1) - t498 * mrSges(4,3);
t449 = t495 * mrSges(6,1) - t457 * mrSges(6,3);
t540 = -t456 * mrSges(6,1) - t457 * mrSges(6,2) - t441;
t397 = m(6) * t410 - t479 * mrSges(6,2) + t425 * t550 - t495 * t449 + t456 * t540 + t534;
t448 = -t495 * mrSges(6,2) - t456 * mrSges(6,3);
t399 = m(6) * t409 + t479 * mrSges(6,1) + t426 * t550 + t495 * t448 + t457 * t540 + t527;
t392 = t513 * t397 + t552 * t399;
t458 = -t489 * mrSges(5,1) + t490 * mrSges(5,2);
t466 = -t497 * mrSges(5,2) + t489 * mrSges(5,3);
t390 = m(5) * t417 + t481 * mrSges(5,1) - t470 * mrSges(5,3) - t490 * t458 + t497 * t466 + t392;
t467 = t497 * mrSges(5,1) - t490 * mrSges(5,3);
t528 = t552 * t397 - t513 * t399;
t391 = m(5) * t418 - t481 * mrSges(5,2) + t469 * mrSges(5,3) + t489 * t458 - t497 * t467 + t528;
t529 = -t509 * t390 + t511 * t391;
t385 = m(4) * t446 - qJDD(3) * mrSges(4,2) - t481 * mrSges(4,3) - qJD(3) * t492 - t497 * t478 + t529;
t445 = t464 * t553 - t514 * t471;
t424 = -qJDD(3) * pkin(3) - t517 * qJ(4) + t498 * t477 + qJDD(4) - t445;
t419 = -t469 * pkin(4) - t488 * pkin(8) + t490 * t468 + t424;
t412 = -0.2e1 * qJD(6) * t457 + (t456 * t495 - t426) * qJ(6) + (t457 * t495 + t425) * pkin(5) + t419;
t405 = m(7) * t412 + t425 * mrSges(7,1) - t426 * mrSges(7,3) + t456 * t447 - t457 * t450;
t520 = m(6) * t419 + t425 * mrSges(6,1) + t426 * mrSges(6,2) + t456 * t448 + t457 * t449 + t405;
t402 = m(5) * t424 - t469 * mrSges(5,1) + t470 * mrSges(5,2) - t489 * t466 + t490 * t467 + t520;
t491 = -qJD(3) * mrSges(4,2) - t497 * mrSges(4,3);
t401 = m(4) * t445 + qJDD(3) * mrSges(4,1) - t482 * mrSges(4,3) + qJD(3) * t491 - t498 * t478 - t402;
t544 = t514 * t385 + t553 * t401;
t386 = t511 * t390 + t509 * t391;
t543 = t546 * t456 - t457 * t548 - t555 * t495;
t530 = t553 * t385 - t514 * t401;
t524 = -t512 * mrSges(3,1) + t510 * mrSges(3,2);
t523 = mrSges(3,3) * qJDD(1) + t518 * t524;
t521 = m(4) * t480 + t481 * mrSges(4,1) + t482 * mrSges(4,2) + t497 * t491 + t498 * t492 + t386;
t501 = (Ifges(3,5) * t510 + Ifges(3,6) * t512) * qJD(1);
t496 = -qJDD(1) * pkin(1) - t518 * qJ(2) + t525;
t483 = -t510 * t499 + t531;
t474 = Ifges(4,1) * t498 - Ifges(4,4) * t497 + Ifges(4,5) * qJD(3);
t473 = Ifges(4,4) * t498 - Ifges(4,2) * t497 + Ifges(4,6) * qJD(3);
t472 = Ifges(4,5) * t498 - Ifges(4,6) * t497 + Ifges(4,3) * qJD(3);
t453 = Ifges(5,1) * t490 + Ifges(5,4) * t489 + Ifges(5,5) * t497;
t452 = Ifges(5,4) * t490 + Ifges(5,2) * t489 + Ifges(5,6) * t497;
t451 = Ifges(5,5) * t490 + Ifges(5,6) * t489 + Ifges(5,3) * t497;
t394 = mrSges(6,2) * t419 + mrSges(7,2) * t408 - mrSges(6,3) * t409 - mrSges(7,3) * t412 - qJ(6) * t405 - t549 * t425 + t557 * t426 + t543 * t456 + t548 * t479 + t542 * t495;
t393 = -mrSges(6,1) * t419 - mrSges(7,1) * t412 + mrSges(7,2) * t407 + mrSges(6,3) * t410 - pkin(5) * t405 - t425 * t556 + t549 * t426 + t543 * t457 + t546 * t479 + t541 * t495;
t382 = mrSges(3,3) * t518 * t539 + m(3) * t496 + qJDD(1) * t524 + t521;
t381 = mrSges(5,2) * t424 - mrSges(5,3) * t417 + Ifges(5,1) * t470 + Ifges(5,4) * t469 + Ifges(5,5) * t481 - pkin(8) * t392 - t513 * t393 + t394 * t552 + t489 * t451 - t497 * t452;
t380 = -mrSges(5,1) * t424 + mrSges(5,3) * t418 + Ifges(5,4) * t470 + Ifges(5,2) * t469 + Ifges(5,6) * t481 - pkin(4) * t520 + pkin(8) * t528 + t393 * t552 + t513 * t394 - t490 * t451 + t497 * t453;
t379 = -pkin(3) * t386 + Ifges(4,6) * qJDD(3) + t554 - t498 * t472 + t489 * t453 - t490 * t452 + Ifges(4,4) * t482 - mrSges(4,1) * t480 - Ifges(5,6) * t469 - Ifges(5,5) * t470 + qJD(3) * t474 + mrSges(4,3) * t446 - mrSges(5,1) * t417 + mrSges(5,2) * t418 + (-Ifges(5,3) - Ifges(4,2)) * t481 - pkin(4) * t392;
t378 = mrSges(4,2) * t480 - mrSges(4,3) * t445 + Ifges(4,1) * t482 - Ifges(4,4) * t481 + Ifges(4,5) * qJDD(3) - qJ(4) * t386 - qJD(3) * t473 - t509 * t380 + t511 * t381 - t497 * t472;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t532 - mrSges(2,2) * t526 + t510 * (t512 * qJD(1) * t501 + mrSges(3,2) * t496 - mrSges(3,3) * t483 + t553 * t378 - t514 * t379 - pkin(7) * t544 + (Ifges(3,1) * t510 + Ifges(3,4) * t512) * qJDD(1)) + t512 * (-t501 * t538 - mrSges(3,1) * t496 + mrSges(3,3) * t484 + t514 * t378 + t553 * t379 - pkin(2) * t521 + pkin(7) * t530 + (Ifges(3,4) * t510 + Ifges(3,2) * t512) * qJDD(1)) - pkin(1) * t382 + qJ(2) * ((m(3) * t484 + t512 * t523 + t530) * t512 + (-m(3) * t483 + t510 * t523 - t544) * t510); t382; mrSges(4,1) * t445 - mrSges(4,2) * t446 + Ifges(4,5) * t482 - Ifges(4,6) * t481 + Ifges(4,3) * qJDD(3) - pkin(3) * t402 + qJ(4) * t529 + t511 * t380 + t509 * t381 + t498 * t473 + t497 * t474; t402; -t554; t404;];
tauJ  = t1;
