% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-05-04 23:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:57:21
% EndTime: 2019-05-04 22:57:23
% DurationCPUTime: 1.68s
% Computational Cost: add. (10827->244), mult. (24311->297), div. (0->0), fcn. (17522->12), ass. (0->116)
t533 = Ifges(5,1) + Ifges(6,2);
t525 = Ifges(5,4) + Ifges(6,6);
t524 = Ifges(5,5) - Ifges(6,4);
t532 = -Ifges(5,2) - Ifges(6,3);
t523 = Ifges(5,6) - Ifges(6,5);
t531 = Ifges(5,3) + Ifges(6,1);
t481 = sin(pkin(10));
t484 = cos(pkin(10));
t466 = g(1) * t481 - g(2) * t484;
t479 = -g(3) + qJDD(1);
t482 = sin(pkin(6));
t485 = cos(pkin(6));
t530 = t466 * t485 + t479 * t482;
t492 = qJD(2) ^ 2;
t467 = -g(1) * t484 - g(2) * t481;
t488 = sin(qJ(2));
t490 = cos(qJ(2));
t428 = -t488 * t467 + t530 * t490;
t483 = cos(pkin(11));
t477 = t483 ^ 2;
t529 = 0.2e1 * t483;
t528 = -2 * qJD(5);
t527 = cos(qJ(4));
t526 = pkin(3) * t483;
t480 = sin(pkin(11));
t522 = mrSges(4,2) * t480;
t520 = t477 * t492;
t429 = t490 * t467 + t530 * t488;
t421 = -pkin(2) * t492 + qJDD(2) * qJ(3) + t429;
t454 = -t466 * t482 + t479 * t485;
t510 = qJD(2) * qJD(3);
t513 = t483 * t454 - 0.2e1 * t480 * t510;
t400 = (-pkin(8) * qJDD(2) + t492 * t526 - t421) * t480 + t513;
t403 = t483 * t421 + t480 * t454 + t510 * t529;
t508 = qJDD(2) * t483;
t401 = -pkin(3) * t520 + pkin(8) * t508 + t403;
t487 = sin(qJ(4));
t392 = t527 * t400 - t487 * t401;
t507 = t483 * t527;
t459 = (t480 * t487 - t507) * qJD(2);
t502 = t527 * t480 + t483 * t487;
t460 = t502 * qJD(2);
t435 = mrSges(5,1) * t459 + mrSges(5,2) * t460;
t511 = t459 * qJD(4);
t443 = t502 * qJDD(2) - t511;
t449 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t459;
t451 = mrSges(6,1) * t459 - qJD(4) * mrSges(6,3);
t434 = pkin(4) * t459 - qJ(5) * t460;
t491 = qJD(4) ^ 2;
t391 = -qJDD(4) * pkin(4) - t491 * qJ(5) + t460 * t434 + qJDD(5) - t392;
t386 = (t459 * t460 - qJDD(4)) * pkin(9) + (t443 + t511) * pkin(5) + t391;
t509 = qJDD(2) * t480;
t512 = qJD(4) * t460;
t442 = -qJDD(2) * t507 + t487 * t509 + t512;
t453 = pkin(5) * t460 - qJD(4) * pkin(9);
t458 = t459 ^ 2;
t476 = t480 ^ 2;
t499 = qJDD(3) - t428;
t413 = (-pkin(2) - t526) * qJDD(2) + (-qJ(3) + (-t476 - t477) * pkin(8)) * t492 + t499;
t493 = pkin(4) * t512 + t460 * t528 + (-t443 + t511) * qJ(5) + t413;
t389 = -pkin(5) * t458 - t453 * t460 + (pkin(4) + pkin(9)) * t442 + t493;
t486 = sin(qJ(6));
t489 = cos(qJ(6));
t384 = t386 * t489 - t389 * t486;
t444 = -qJD(4) * t486 + t459 * t489;
t412 = qJD(6) * t444 + qJDD(4) * t489 + t442 * t486;
t445 = qJD(4) * t489 + t459 * t486;
t414 = -mrSges(7,1) * t444 + mrSges(7,2) * t445;
t457 = qJD(6) + t460;
t419 = -mrSges(7,2) * t457 + mrSges(7,3) * t444;
t441 = qJDD(6) + t443;
t381 = m(7) * t384 + mrSges(7,1) * t441 - mrSges(7,3) * t412 - t414 * t445 + t419 * t457;
t385 = t386 * t486 + t389 * t489;
t411 = -qJD(6) * t445 - qJDD(4) * t486 + t442 * t489;
t420 = mrSges(7,1) * t457 - mrSges(7,3) * t445;
t382 = m(7) * t385 - mrSges(7,2) * t441 + mrSges(7,3) * t411 + t414 * t444 - t420 * t457;
t373 = t381 * t489 + t382 * t486;
t436 = -mrSges(6,2) * t459 - mrSges(6,3) * t460;
t500 = -m(6) * t391 - t443 * mrSges(6,1) - t460 * t436 - t373;
t369 = m(5) * t392 - mrSges(5,3) * t443 - t435 * t460 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t449 - t451) * qJD(4) + t500;
t393 = t487 * t400 + t527 * t401;
t450 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t460;
t498 = -pkin(4) * t491 + qJDD(4) * qJ(5) - t434 * t459 + t393;
t390 = qJD(4) * t528 - t498;
t452 = mrSges(6,1) * t460 + qJD(4) * mrSges(6,2);
t388 = -pkin(5) * t442 - pkin(9) * t458 + ((2 * qJD(5)) + t453) * qJD(4) + t498;
t501 = -m(7) * t388 + mrSges(7,1) * t411 - t412 * mrSges(7,2) + t419 * t444 - t445 * t420;
t496 = -m(6) * t390 + qJDD(4) * mrSges(6,3) + qJD(4) * t452 - t501;
t378 = m(5) * t393 - qJDD(4) * mrSges(5,2) - qJD(4) * t450 + (-t435 - t436) * t459 + (-mrSges(5,3) - mrSges(6,1)) * t442 + t496;
t518 = t527 * t369 + t487 * t378;
t517 = -t486 * t381 + t489 * t382;
t516 = -t531 * qJD(4) + t523 * t459 - t524 * t460;
t515 = t523 * qJD(4) + t532 * t459 + t525 * t460;
t514 = t524 * qJD(4) - t525 * t459 + t533 * t460;
t402 = -t421 * t480 + t513;
t503 = mrSges(4,3) * qJDD(2) + t492 * (-mrSges(4,1) * t483 + t522);
t366 = m(4) * t402 - t503 * t480 + t518;
t505 = -t369 * t487 + t527 * t378;
t367 = m(4) * t403 + t503 * t483 + t505;
t506 = -t366 * t480 + t483 * t367;
t395 = pkin(4) * t442 + t493;
t372 = m(6) * t395 - t442 * mrSges(6,2) - t443 * mrSges(6,3) - t459 * t451 - t460 * t452 + t517;
t405 = Ifges(7,4) * t445 + Ifges(7,2) * t444 + Ifges(7,6) * t457;
t406 = Ifges(7,1) * t445 + Ifges(7,4) * t444 + Ifges(7,5) * t457;
t497 = mrSges(7,1) * t384 - mrSges(7,2) * t385 + Ifges(7,5) * t412 + Ifges(7,6) * t411 + Ifges(7,3) * t441 + t445 * t405 - t444 * t406;
t495 = m(5) * t413 + t442 * mrSges(5,1) + mrSges(5,2) * t443 + t459 * t449 + t450 * t460 + t372;
t418 = -qJDD(2) * pkin(2) - qJ(3) * t492 + t499;
t494 = -m(4) * t418 + mrSges(4,1) * t508 - t495 + (t476 * t492 + t520) * mrSges(4,3);
t404 = Ifges(7,5) * t445 + Ifges(7,6) * t444 + Ifges(7,3) * t457;
t375 = mrSges(7,2) * t388 - mrSges(7,3) * t384 + Ifges(7,1) * t412 + Ifges(7,4) * t411 + Ifges(7,5) * t441 + t404 * t444 - t405 * t457;
t374 = -mrSges(7,1) * t388 + mrSges(7,3) * t385 + Ifges(7,4) * t412 + Ifges(7,2) * t411 + Ifges(7,6) * t441 - t404 * t445 + t406 * t457;
t371 = qJDD(4) * mrSges(6,2) + qJD(4) * t451 - t500;
t370 = mrSges(4,2) * t509 - t494;
t364 = mrSges(6,1) * t391 + mrSges(5,2) * t413 - mrSges(5,3) * t392 - mrSges(6,3) * t395 + pkin(5) * t373 - qJ(5) * t372 - t515 * qJD(4) + t524 * qJDD(4) - t525 * t442 + t533 * t443 + t516 * t459 + t497;
t363 = -mrSges(5,1) * t413 - mrSges(6,1) * t390 + mrSges(6,2) * t395 + mrSges(5,3) * t393 - pkin(4) * t372 - pkin(5) * t501 - pkin(9) * t517 + t514 * qJD(4) + t523 * qJDD(4) - t489 * t374 - t486 * t375 + t532 * t442 + t525 * t443 + t516 * t460;
t1 = [m(2) * t479 + t485 * (m(3) * t454 + t366 * t483 + t367 * t480) + (t488 * (m(3) * t429 - mrSges(3,1) * t492 - qJDD(2) * mrSges(3,2) + t506) + t490 * (t494 + (mrSges(3,1) - t522) * qJDD(2) + m(3) * t428 - mrSges(3,2) * t492)) * t482; mrSges(3,1) * t428 - mrSges(3,2) * t429 + t480 * (mrSges(4,2) * t418 - mrSges(4,3) * t402 - pkin(8) * t518 - t487 * t363 + t527 * t364) + t483 * (-mrSges(4,1) * t418 + mrSges(4,3) * t403 - pkin(3) * t495 + pkin(8) * t505 + t527 * t363 + t487 * t364) - pkin(2) * t370 + qJ(3) * t506 + (Ifges(4,2) * t477 + Ifges(3,3) + (Ifges(4,1) * t480 + Ifges(4,4) * t529) * t480) * qJDD(2); t370; mrSges(5,1) * t392 - mrSges(5,2) * t393 + mrSges(6,2) * t391 - mrSges(6,3) * t390 + t489 * t375 - t486 * t374 - pkin(9) * t373 - pkin(4) * t371 + qJ(5) * t496 + t515 * t460 + (-qJ(5) * t436 + t514) * t459 + t524 * t443 + (-qJ(5) * mrSges(6,1) - t523) * t442 + t531 * qJDD(4); t371; t497;];
tauJ  = t1;
