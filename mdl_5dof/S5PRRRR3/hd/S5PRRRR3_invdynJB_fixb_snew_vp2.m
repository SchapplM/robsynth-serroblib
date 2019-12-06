% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:15
% EndTime: 2019-12-05 17:06:17
% DurationCPUTime: 2.35s
% Computational Cost: add. (33321->174), mult. (42036->224), div. (0->0), fcn. (27482->10), ass. (0->85)
t471 = sin(pkin(9));
t472 = cos(pkin(9));
t455 = g(1) * t471 - g(2) * t472;
t456 = -g(1) * t472 - g(2) * t471;
t476 = sin(qJ(2));
t480 = cos(qJ(2));
t440 = t480 * t455 - t456 * t476;
t437 = qJDD(2) * pkin(2) + t440;
t441 = t476 * t455 + t480 * t456;
t481 = qJD(2) ^ 2;
t438 = -pkin(2) * t481 + t441;
t475 = sin(qJ(3));
t479 = cos(qJ(3));
t432 = t479 * t437 - t438 * t475;
t467 = qJDD(2) + qJDD(3);
t429 = pkin(3) * t467 + t432;
t433 = t475 * t437 + t479 * t438;
t468 = qJD(2) + qJD(3);
t466 = t468 ^ 2;
t430 = -pkin(3) * t466 + t433;
t474 = sin(qJ(4));
t478 = cos(qJ(4));
t426 = t474 * t429 + t478 * t430;
t461 = qJD(4) + t468;
t459 = t461 ^ 2;
t460 = qJDD(4) + t467;
t423 = -pkin(4) * t459 + pkin(8) * t460 + t426;
t470 = -g(3) + qJDD(1);
t473 = sin(qJ(5));
t477 = cos(qJ(5));
t420 = -t423 * t473 + t470 * t477;
t421 = t423 * t477 + t470 * t473;
t443 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t473 + Ifges(6,2) * t477) * t461;
t444 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t473 + Ifges(6,4) * t477) * t461;
t497 = qJD(5) * t461;
t448 = t460 * t473 + t477 * t497;
t449 = t460 * t477 - t473 * t497;
t501 = mrSges(6,1) * t420 - mrSges(6,2) * t421 + Ifges(6,5) * t448 + Ifges(6,6) * t449 + Ifges(6,3) * qJDD(5) + (t443 * t473 - t444 * t477) * t461;
t500 = t461 * t473;
t499 = t461 * t477;
t447 = (-mrSges(6,1) * t477 + mrSges(6,2) * t473) * t461;
t451 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t499;
t418 = m(6) * t420 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t448 + qJD(5) * t451 - t447 * t500;
t450 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t500;
t419 = m(6) * t421 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t449 - qJD(5) * t450 + t447 * t499;
t491 = -t418 * t473 + t477 * t419;
t405 = m(5) * t426 - mrSges(5,1) * t459 - mrSges(5,2) * t460 + t491;
t425 = t429 * t478 - t430 * t474;
t422 = -pkin(4) * t460 - pkin(8) * t459 - t425;
t485 = -m(6) * t422 + t449 * mrSges(6,1) - mrSges(6,2) * t448 - t450 * t500 + t451 * t499;
t413 = m(5) * t425 + mrSges(5,1) * t460 - mrSges(5,2) * t459 + t485;
t402 = t405 * t474 + t478 * t413;
t398 = m(4) * t432 + mrSges(4,1) * t467 - mrSges(4,2) * t466 + t402;
t492 = t405 * t478 - t413 * t474;
t399 = m(4) * t433 - mrSges(4,1) * t466 - mrSges(4,2) * t467 + t492;
t393 = t398 * t479 + t399 * t475;
t390 = m(3) * t440 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t481 + t393;
t493 = -t398 * t475 + t399 * t479;
t391 = m(3) * t441 - mrSges(3,1) * t481 - qJDD(2) * mrSges(3,2) + t493;
t384 = t390 * t480 + t391 * t476;
t382 = m(2) * t455 + t384;
t494 = -t390 * t476 + t391 * t480;
t383 = m(2) * t456 + t494;
t498 = t382 * t472 + t383 * t471;
t407 = t477 * t418 + t473 * t419;
t496 = m(5) * t470 + t407;
t495 = -t382 * t471 + t383 * t472;
t490 = m(4) * t470 + t496;
t489 = m(3) * t470 + t490;
t487 = m(2) * t470 + t489;
t442 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t473 + Ifges(6,6) * t477) * t461;
t410 = -mrSges(6,1) * t422 + mrSges(6,3) * t421 + Ifges(6,4) * t448 + Ifges(6,2) * t449 + Ifges(6,6) * qJDD(5) + qJD(5) * t444 - t442 * t500;
t411 = mrSges(6,2) * t422 - mrSges(6,3) * t420 + Ifges(6,1) * t448 + Ifges(6,4) * t449 + Ifges(6,5) * qJDD(5) - qJD(5) * t443 + t442 * t499;
t486 = mrSges(5,1) * t425 - mrSges(5,2) * t426 + Ifges(5,3) * t460 + pkin(4) * t485 + pkin(8) * t491 + t410 * t477 + t411 * t473;
t483 = mrSges(4,1) * t432 - mrSges(4,2) * t433 + Ifges(4,3) * t467 + pkin(3) * t402 + t486;
t482 = mrSges(3,1) * t440 - mrSges(3,2) * t441 + Ifges(3,3) * qJDD(2) + pkin(2) * t393 + t483;
t400 = -mrSges(5,1) * t470 + mrSges(5,3) * t426 + t459 * Ifges(5,5) + Ifges(5,6) * t460 - pkin(4) * t407 - t501;
t394 = mrSges(5,2) * t470 - mrSges(5,3) * t425 + Ifges(5,5) * t460 - Ifges(5,6) * t459 - pkin(8) * t407 - t410 * t473 + t411 * t477;
t386 = mrSges(4,2) * t470 - mrSges(4,3) * t432 + Ifges(4,5) * t467 - Ifges(4,6) * t466 - pkin(7) * t402 + t394 * t478 - t400 * t474;
t385 = -mrSges(4,1) * t470 + mrSges(4,3) * t433 + t466 * Ifges(4,5) + Ifges(4,6) * t467 - pkin(3) * t496 + pkin(7) * t492 + t474 * t394 + t478 * t400;
t378 = mrSges(3,2) * t470 - mrSges(3,3) * t440 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t481 - pkin(6) * t393 - t385 * t475 + t386 * t479;
t377 = -mrSges(3,1) * t470 + mrSges(3,3) * t441 + t481 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t490 + pkin(6) * t493 + t479 * t385 + t475 * t386;
t376 = mrSges(2,2) * t470 - mrSges(2,3) * t455 - pkin(5) * t384 - t377 * t476 + t378 * t480;
t375 = -mrSges(2,1) * t470 + mrSges(2,3) * t456 - pkin(1) * t489 + pkin(5) * t494 + t480 * t377 + t476 * t378;
t1 = [-m(1) * g(1) + t495; -m(1) * g(2) + t498; -m(1) * g(3) + t487; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t498 - t375 * t471 + t376 * t472; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t495 + t472 * t375 + t471 * t376; -mrSges(1,1) * g(2) + mrSges(2,1) * t455 + mrSges(1,2) * g(1) - mrSges(2,2) * t456 + pkin(1) * t384 + t482; t487; t482; t483; t486; t501;];
tauJB = t1;
