% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 09:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRPR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:06:16
% EndTime: 2019-05-05 09:06:20
% DurationCPUTime: 3.31s
% Computational Cost: add. (29054->283), mult. (60090->353), div. (0->0), fcn. (46797->14), ass. (0->129)
t547 = Ifges(5,1) + Ifges(6,2);
t537 = Ifges(5,4) + Ifges(6,6);
t536 = Ifges(5,5) - Ifges(6,4);
t546 = -Ifges(5,2) - Ifges(6,3);
t535 = Ifges(5,6) - Ifges(6,5);
t543 = Ifges(5,3) + Ifges(6,1);
t493 = sin(pkin(12));
t496 = cos(pkin(12));
t485 = g(1) * t493 - g(2) * t496;
t492 = -g(3) + qJDD(1);
t495 = sin(pkin(6));
t498 = cos(pkin(6));
t545 = t485 * t498 + t492 * t495;
t497 = cos(pkin(7));
t491 = qJD(2) * t497 + qJD(3);
t500 = sin(qJ(4));
t501 = sin(qJ(3));
t494 = sin(pkin(7));
t523 = qJD(2) * t494;
t521 = t501 * t523;
t539 = cos(qJ(4));
t466 = -t491 * t539 + t500 * t521;
t504 = cos(qJ(3));
t520 = t504 * t523;
t483 = -qJD(4) + t520;
t453 = mrSges(6,1) * t466 + mrSges(6,3) * t483;
t522 = qJD(2) * qJD(3);
t476 = (-qJDD(2) * t504 + t501 * t522) * t494;
t470 = qJDD(4) + t476;
t486 = -g(1) * t496 - g(2) * t493;
t502 = sin(qJ(2));
t505 = cos(qJ(2));
t447 = -t486 * t502 + t505 * t545;
t506 = qJD(2) ^ 2;
t538 = pkin(9) * t494;
t441 = qJDD(2) * pkin(2) + t506 * t538 + t447;
t448 = t505 * t486 + t502 * t545;
t442 = -pkin(2) * t506 + qJDD(2) * t538 + t448;
t529 = t497 * t501;
t468 = -t485 * t495 + t492 * t498;
t532 = t468 * t494;
t407 = t441 * t529 + t504 * t442 + t501 * t532;
t474 = (-pkin(3) * t504 - pkin(10) * t501) * t523;
t489 = t491 ^ 2;
t490 = qJDD(2) * t497 + qJDD(3);
t402 = -pkin(3) * t489 + pkin(10) * t490 + t474 * t520 + t407;
t462 = t497 * t468;
t475 = (qJDD(2) * t501 + t504 * t522) * t494;
t404 = t476 * pkin(3) - t475 * pkin(10) + t462 + (-t441 + (pkin(3) * t501 - pkin(10) * t504) * t491 * qJD(2)) * t494;
t397 = -t500 * t402 + t404 * t539;
t467 = t500 * t491 + t521 * t539;
t443 = pkin(4) * t466 - qJ(5) * t467;
t482 = t483 ^ 2;
t395 = -t470 * pkin(4) - t482 * qJ(5) + t467 * t443 + qJDD(5) - t397;
t437 = -t466 * qJD(4) + t475 * t539 + t500 * t490;
t533 = t466 * t483;
t390 = (t466 * t467 - t470) * pkin(11) + (t437 - t533) * pkin(5) + t395;
t436 = qJD(4) * t467 + t475 * t500 - t490 * t539;
t455 = pkin(5) * t467 + pkin(11) * t483;
t465 = t466 ^ 2;
t406 = -t501 * t442 + t504 * (t441 * t497 + t532);
t401 = -t490 * pkin(3) - t489 * pkin(10) + t474 * t521 - t406;
t540 = -2 * qJD(5);
t508 = (-t437 - t533) * qJ(5) + t401 + (-t483 * pkin(4) + t540) * t467;
t393 = -pkin(5) * t465 - t455 * t467 + (pkin(4) + pkin(11)) * t436 + t508;
t499 = sin(qJ(6));
t503 = cos(qJ(6));
t388 = t390 * t503 - t393 * t499;
t449 = t466 * t503 + t483 * t499;
t412 = qJD(6) * t449 + t436 * t499 + t470 * t503;
t450 = t466 * t499 - t483 * t503;
t417 = -mrSges(7,1) * t449 + mrSges(7,2) * t450;
t464 = qJD(6) + t467;
t421 = -mrSges(7,2) * t464 + mrSges(7,3) * t449;
t433 = qJDD(6) + t437;
t385 = m(7) * t388 + mrSges(7,1) * t433 - mrSges(7,3) * t412 - t417 * t450 + t421 * t464;
t389 = t390 * t499 + t393 * t503;
t411 = -qJD(6) * t450 + t436 * t503 - t470 * t499;
t422 = mrSges(7,1) * t464 - mrSges(7,3) * t450;
t386 = m(7) * t389 - mrSges(7,2) * t433 + mrSges(7,3) * t411 + t417 * t449 - t422 * t464;
t376 = t385 * t503 + t386 * t499;
t445 = -mrSges(6,2) * t466 - mrSges(6,3) * t467;
t513 = -m(6) * t395 - t437 * mrSges(6,1) - t467 * t445 - t376;
t374 = mrSges(6,2) * t470 - t453 * t483 - t513;
t398 = t539 * t402 + t500 * t404;
t512 = -t482 * pkin(4) + t470 * qJ(5) - t466 * t443 + t398;
t392 = -t436 * pkin(5) - t465 * pkin(11) + (t540 - t455) * t483 + t512;
t413 = Ifges(7,5) * t450 + Ifges(7,6) * t449 + Ifges(7,3) * t464;
t415 = Ifges(7,1) * t450 + Ifges(7,4) * t449 + Ifges(7,5) * t464;
t377 = -mrSges(7,1) * t392 + mrSges(7,3) * t389 + Ifges(7,4) * t412 + Ifges(7,2) * t411 + Ifges(7,6) * t433 - t413 * t450 + t415 * t464;
t414 = Ifges(7,4) * t450 + Ifges(7,2) * t449 + Ifges(7,6) * t464;
t378 = mrSges(7,2) * t392 - mrSges(7,3) * t388 + Ifges(7,1) * t412 + Ifges(7,4) * t411 + Ifges(7,5) * t433 + t413 * t449 - t414 * t464;
t394 = 0.2e1 * qJD(5) * t483 - t512;
t454 = mrSges(6,1) * t467 - mrSges(6,2) * t483;
t514 = -m(7) * t392 + mrSges(7,1) * t411 - t412 * mrSges(7,2) + t421 * t449 - t450 * t422;
t510 = -m(6) * t394 + t470 * mrSges(6,3) - t483 * t454 - t514;
t524 = -t537 * t466 + t547 * t467 - t536 * t483;
t525 = t546 * t466 + t537 * t467 - t535 * t483;
t544 = -t535 * t436 + t536 * t437 + t543 * t470 + mrSges(5,1) * t397 - mrSges(5,2) * t398 + mrSges(6,2) * t395 - mrSges(6,3) * t394 - pkin(4) * t374 - pkin(11) * t376 + qJ(5) * (-mrSges(6,1) * t436 - t445 * t466 + t510) - t499 * t377 + t503 * t378 + t525 * t467 + t524 * t466;
t472 = -mrSges(4,2) * t491 + mrSges(4,3) * t520;
t473 = (-mrSges(4,1) * t504 + mrSges(4,2) * t501) * t523;
t451 = mrSges(5,2) * t483 - mrSges(5,3) * t466;
t452 = -mrSges(5,1) * t483 - mrSges(5,3) * t467;
t396 = pkin(4) * t436 + t508;
t527 = -t499 * t385 + t503 * t386;
t516 = -m(6) * t396 + t436 * mrSges(6,2) + t466 * t453 - t527;
t509 = -m(5) * t401 - t436 * mrSges(5,1) - t466 * t451 + (-t452 + t454) * t467 + (-mrSges(5,2) + mrSges(6,3)) * t437 + t516;
t371 = m(4) * t406 + mrSges(4,1) * t490 - mrSges(4,3) * t475 + t472 * t491 - t473 * t521 + t509;
t534 = t371 * t504;
t471 = mrSges(4,1) * t491 - mrSges(4,3) * t521;
t444 = mrSges(5,1) * t466 + mrSges(5,2) * t467;
t373 = m(5) * t397 - mrSges(5,3) * t437 - t444 * t467 + (-t451 + t453) * t483 + (mrSges(5,1) - mrSges(6,2)) * t470 + t513;
t381 = m(5) * t398 - mrSges(5,2) * t470 + t452 * t483 + (-t444 - t445) * t466 + (-mrSges(5,3) - mrSges(6,1)) * t436 + t510;
t518 = -t373 * t500 + t539 * t381;
t367 = m(4) * t407 - mrSges(4,2) * t490 - mrSges(4,3) * t476 - t471 * t491 + t473 * t520 + t518;
t528 = t367 * t529 + t497 * t534;
t369 = t539 * t373 + t500 * t381;
t526 = t535 * t466 - t467 * t536 + t543 * t483;
t519 = t504 * t367 - t371 * t501;
t511 = mrSges(7,1) * t388 - mrSges(7,2) * t389 + Ifges(7,5) * t412 + Ifges(7,6) * t411 + Ifges(7,3) * t433 + t450 * t414 - t449 * t415;
t460 = Ifges(4,5) * t491 + (Ifges(4,1) * t501 + Ifges(4,4) * t504) * t523;
t459 = Ifges(4,6) * t491 + (Ifges(4,4) * t501 + Ifges(4,2) * t504) * t523;
t418 = -t494 * t441 + t462;
t375 = -mrSges(6,3) * t437 - t454 * t467 - t516;
t368 = m(4) * t418 + t476 * mrSges(4,1) + t475 * mrSges(4,2) + (t471 * t501 - t472 * t504) * t523 + t369;
t364 = mrSges(6,1) * t395 + mrSges(5,2) * t401 - mrSges(5,3) * t397 - mrSges(6,3) * t396 + pkin(5) * t376 - qJ(5) * t375 - t537 * t436 + t547 * t437 + t526 * t466 + t536 * t470 + t525 * t483 + t511;
t363 = -mrSges(5,1) * t401 - mrSges(6,1) * t394 + mrSges(6,2) * t396 + mrSges(5,3) * t398 - pkin(4) * t375 - pkin(5) * t514 - pkin(11) * t527 - t503 * t377 - t499 * t378 + t546 * t436 + t537 * t437 + t526 * t467 + t535 * t470 - t524 * t483;
t362 = Ifges(4,5) * t475 - Ifges(4,6) * t476 + Ifges(4,3) * t490 + mrSges(4,1) * t406 - mrSges(4,2) * t407 + t500 * t364 + t539 * t363 + pkin(3) * t509 + pkin(10) * t518 + (t459 * t501 - t460 * t504) * t523;
t1 = [m(2) * t492 + t498 * (m(3) * t468 + t497 * t368 + (t367 * t501 + t534) * t494) + (t502 * (m(3) * t448 - mrSges(3,1) * t506 - qJDD(2) * mrSges(3,2) + t519) + t505 * (m(3) * t447 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t506 - t368 * t494 + t528)) * t495; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t447 - mrSges(3,2) * t448 + t497 * t362 + pkin(2) * t528 + (t501 * (mrSges(4,2) * t418 - mrSges(4,3) * t406 + Ifges(4,1) * t475 - Ifges(4,4) * t476 + Ifges(4,5) * t490 - pkin(10) * t369 - t500 * t363 + t364 * t539 - t491 * t459) + t504 * (-mrSges(4,1) * t418 + mrSges(4,3) * t407 + Ifges(4,4) * t475 - Ifges(4,2) * t476 + Ifges(4,6) * t490 - pkin(3) * t369 + t491 * t460 - t544) - pkin(2) * t368 + pkin(9) * t519) * t494; t362; t544; t374; t511;];
tauJ  = t1;
