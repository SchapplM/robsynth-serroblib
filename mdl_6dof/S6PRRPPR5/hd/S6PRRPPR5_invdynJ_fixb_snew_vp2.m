% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-05-05 03:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:24:57
% EndTime: 2019-05-05 03:25:01
% DurationCPUTime: 2.15s
% Computational Cost: add. (13261->272), mult. (27708->334), div. (0->0), fcn. (17353->12), ass. (0->119)
t500 = cos(qJ(3));
t534 = Ifges(4,4) + Ifges(5,6);
t544 = t500 * t534;
t497 = sin(qJ(3));
t543 = t497 * t534 + t500 * (Ifges(4,2) + Ifges(5,3));
t542 = -2 * qJD(4);
t541 = Ifges(4,1) + Ifges(5,2);
t533 = (Ifges(4,5) - Ifges(5,4));
t532 = (Ifges(4,6) - Ifges(5,5));
t491 = sin(pkin(10));
t494 = cos(pkin(10));
t473 = g(1) * t491 - g(2) * t494;
t489 = -g(3) + qJDD(1);
t492 = sin(pkin(6));
t495 = cos(pkin(6));
t539 = t473 * t495 + t489 * t492;
t474 = -g(1) * t494 - g(2) * t491;
t498 = sin(qJ(2));
t501 = cos(qJ(2));
t426 = t501 * t474 + t539 * t498;
t503 = qJD(2) ^ 2;
t422 = -pkin(2) * t503 + qJDD(2) * pkin(8) + t426;
t441 = -t473 * t492 + t489 * t495;
t417 = t500 * t422 + t497 * t441;
t466 = (-t500 * pkin(3) - t497 * qJ(4)) * qJD(2);
t502 = qJD(3) ^ 2;
t523 = qJD(2) * t500;
t404 = pkin(3) * t502 - qJDD(3) * qJ(4) + (qJD(3) * t542) - t466 * t523 - t417;
t425 = -t498 * t474 + t501 * t539;
t506 = -qJDD(2) * pkin(2) - t425;
t536 = pkin(8) * t503;
t421 = t506 - t536;
t521 = qJD(2) * qJD(3);
t519 = t500 * t521;
t469 = qJDD(2) * t497 + t519;
t518 = t497 * t521;
t470 = qJDD(2) * t500 - t518;
t522 = t497 * qJD(2);
t475 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t522;
t476 = -(qJD(3) * mrSges(4,2)) + mrSges(4,3) * t523;
t478 = -mrSges(5,1) * t523 - (qJD(3) * mrSges(5,3));
t505 = pkin(3) * t518 + t522 * t542 + (-t469 - t519) * qJ(4) + t506;
t406 = -pkin(3) * t470 + t505 - t536;
t479 = mrSges(5,1) * t522 + (qJD(3) * mrSges(5,2));
t419 = t497 * t422;
t514 = -qJ(4) * t502 + t466 * t522 + qJDD(4) + t419;
t530 = qJ(5) * t503;
t531 = -pkin(3) - qJ(5);
t400 = pkin(4) * t469 + t531 * qJDD(3) + (-pkin(4) * t521 - t497 * t530 - t441) * t500 + t514;
t477 = pkin(4) * t522 - qJD(3) * qJ(5);
t488 = t500 ^ 2;
t402 = -t477 * t522 + (-pkin(4) * t488 - pkin(8)) * t503 + t531 * t470 + t505;
t490 = sin(pkin(11));
t493 = cos(pkin(11));
t459 = qJD(3) * t493 - t490 * t523;
t392 = -0.2e1 * qJD(5) * t459 + t493 * t400 - t402 * t490;
t439 = qJDD(3) * t493 - t470 * t490;
t458 = -qJD(3) * t490 - t493 * t523;
t390 = (t458 * t522 - t439) * pkin(9) + (t458 * t459 + t469) * pkin(5) + t392;
t393 = 0.2e1 * qJD(5) * t458 + t490 * t400 + t493 * t402;
t438 = -qJDD(3) * t490 - t470 * t493;
t440 = pkin(5) * t522 - pkin(9) * t459;
t457 = t458 ^ 2;
t391 = -pkin(5) * t457 + pkin(9) * t438 - t440 * t522 + t393;
t496 = sin(qJ(6));
t499 = cos(qJ(6));
t388 = t390 * t499 - t391 * t496;
t431 = t458 * t499 - t459 * t496;
t411 = qJD(6) * t431 + t438 * t496 + t439 * t499;
t432 = t458 * t496 + t459 * t499;
t418 = -mrSges(7,1) * t431 + mrSges(7,2) * t432;
t482 = qJD(6) + t522;
t423 = -mrSges(7,2) * t482 + mrSges(7,3) * t431;
t463 = qJDD(6) + t469;
t384 = m(7) * t388 + mrSges(7,1) * t463 - mrSges(7,3) * t411 - t418 * t432 + t423 * t482;
t389 = t390 * t496 + t391 * t499;
t410 = -qJD(6) * t432 + t438 * t499 - t439 * t496;
t424 = mrSges(7,1) * t482 - mrSges(7,3) * t432;
t385 = m(7) * t389 - mrSges(7,2) * t463 + mrSges(7,3) * t410 + t418 * t431 - t424 * t482;
t376 = t499 * t384 + t496 * t385;
t433 = -mrSges(6,1) * t458 + mrSges(6,2) * t459;
t436 = -mrSges(6,2) * t522 + mrSges(6,3) * t458;
t374 = m(6) * t392 + mrSges(6,1) * t469 - mrSges(6,3) * t439 - t433 * t459 + t436 * t522 + t376;
t437 = mrSges(6,1) * t522 - mrSges(6,3) * t459;
t516 = -t384 * t496 + t499 * t385;
t375 = m(6) * t393 - mrSges(6,2) * t469 + mrSges(6,3) * t438 + t433 * t458 - t437 * t522 + t516;
t526 = -t490 * t374 + t493 * t375;
t513 = m(5) * t406 + t470 * mrSges(5,2) - t479 * t522 + t526;
t538 = -m(4) * t421 + t470 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t469 + t476 * t523 + (-t475 * t497 - t478 * t500) * qJD(2) - t513;
t537 = t497 * (-t543 * qJD(2) - (t532 * qJD(3))) + t500 * ((t533 * qJD(3)) + (t541 * t497 + t544) * qJD(2));
t529 = t441 * t500;
t416 = -t419 + t529;
t467 = (t500 * mrSges(5,2) - t497 * mrSges(5,3)) * qJD(2);
t468 = (-t500 * mrSges(4,1) + t497 * mrSges(4,2)) * qJD(2);
t371 = t374 * t493 + t375 * t490;
t405 = -qJDD(3) * pkin(3) + t514 - t529;
t509 = -m(5) * t405 - t469 * mrSges(5,1) - t371;
t368 = m(4) * t416 - mrSges(4,3) * t469 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t476 - t478) * qJD(3) + (-t467 - t468) * t522 + t509;
t399 = pkin(4) * t470 + qJD(3) * t477 - t488 * t530 + qJDD(5) - t404;
t395 = -pkin(5) * t438 - pkin(9) * t457 + t440 * t459 + t399;
t510 = m(7) * t395 - mrSges(7,1) * t410 + t411 * mrSges(7,2) - t423 * t431 + t432 * t424;
t386 = m(6) * t399 - mrSges(6,1) * t438 + t439 * mrSges(6,2) - t436 * t458 + t459 * t437 + t510;
t504 = -m(5) * t404 + qJDD(3) * mrSges(5,3) + qJD(3) * t479 + t467 * t523 + t386;
t380 = t468 * t523 + t504 + m(4) * t417 - qJD(3) * t475 + (mrSges(4,3) + mrSges(5,1)) * t470 - qJDD(3) * mrSges(4,2);
t517 = -t368 * t497 + t500 * t380;
t413 = Ifges(7,4) * t432 + Ifges(7,2) * t431 + Ifges(7,6) * t482;
t414 = Ifges(7,1) * t432 + Ifges(7,4) * t431 + Ifges(7,5) * t482;
t508 = mrSges(7,1) * t388 - mrSges(7,2) * t389 + Ifges(7,5) * t411 + Ifges(7,6) * t410 + Ifges(7,3) * t463 + t432 * t413 - t431 * t414;
t429 = Ifges(6,1) * t459 + Ifges(6,4) * t458 + Ifges(6,5) * t522;
t428 = Ifges(6,4) * t459 + Ifges(6,2) * t458 + Ifges(6,6) * t522;
t427 = Ifges(6,5) * t459 + Ifges(6,6) * t458 + Ifges(6,3) * t522;
t412 = Ifges(7,5) * t432 + Ifges(7,6) * t431 + Ifges(7,3) * t482;
t378 = mrSges(7,2) * t395 - mrSges(7,3) * t388 + Ifges(7,1) * t411 + Ifges(7,4) * t410 + Ifges(7,5) * t463 + t412 * t431 - t413 * t482;
t377 = -mrSges(7,1) * t395 + mrSges(7,3) * t389 + Ifges(7,4) * t411 + Ifges(7,2) * t410 + Ifges(7,6) * t463 - t412 * t432 + t414 * t482;
t370 = qJDD(3) * mrSges(5,2) + qJD(3) * t478 + t467 * t522 - t509;
t369 = -mrSges(5,3) * t469 + t478 * t523 + t513;
t367 = mrSges(6,2) * t399 - mrSges(6,3) * t392 + Ifges(6,1) * t439 + Ifges(6,4) * t438 + Ifges(6,5) * t469 - pkin(9) * t376 - t377 * t496 + t378 * t499 + t427 * t458 - t428 * t522;
t366 = -mrSges(6,1) * t399 + mrSges(6,3) * t393 + Ifges(6,4) * t439 + Ifges(6,2) * t438 + Ifges(6,6) * t469 - pkin(5) * t510 + pkin(9) * t516 + t499 * t377 + t496 * t378 - t459 * t427 + t429 * t522;
t1 = [m(2) * t489 + t495 * (m(3) * t441 + t368 * t500 + t380 * t497) + (t498 * (m(3) * t426 - mrSges(3,1) * t503 - qJDD(2) * mrSges(3,2) + t517) + t501 * (m(3) * t425 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t503 + t538)) * t492; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t425 - mrSges(3,2) * t426 + t497 * (mrSges(5,1) * t405 + mrSges(6,1) * t392 + mrSges(4,2) * t421 - mrSges(6,2) * t393 - mrSges(4,3) * t416 - mrSges(5,3) * t406 + Ifges(6,5) * t439 + Ifges(6,6) * t438 + pkin(4) * t371 + pkin(5) * t376 - qJ(4) * t369 + t459 * t428 - t458 * t429 + t508) + t500 * (-mrSges(4,1) * t421 - mrSges(5,1) * t404 + mrSges(5,2) * t406 + mrSges(4,3) * t417 - pkin(3) * t369 + pkin(4) * t386 - qJ(5) * t526 - t493 * t366 - t490 * t367) + pkin(8) * t517 + t543 * t470 + (t497 * t533 + t500 * t532) * qJDD(3) + t537 * qJD(3) + (t497 * (Ifges(6,3) + t541) + t544) * t469 + t538 * pkin(2); mrSges(4,1) * t416 - mrSges(4,2) * t417 + mrSges(5,2) * t405 - mrSges(5,3) * t404 + t493 * t367 - t490 * t366 - qJ(5) * t371 - pkin(3) * t370 + qJ(4) * t504 + (qJ(4) * mrSges(5,1) + t532) * t470 + t533 * t469 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) - t537 * qJD(2); t370; t386; t508;];
tauJ  = t1;
