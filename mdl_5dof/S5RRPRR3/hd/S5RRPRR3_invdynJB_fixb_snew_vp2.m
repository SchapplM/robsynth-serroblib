% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:32
% EndTime: 2020-01-03 12:00:34
% DurationCPUTime: 2.44s
% Computational Cost: add. (44455->186), mult. (53565->233), div. (0->0), fcn. (28670->10), ass. (0->85)
t497 = sin(qJ(1));
t501 = cos(qJ(1));
t477 = -t501 * g(2) - t497 * g(3);
t473 = qJDD(1) * pkin(1) + t477;
t476 = -t497 * g(2) + t501 * g(3);
t502 = qJD(1) ^ 2;
t474 = -t502 * pkin(1) + t476;
t496 = sin(qJ(2));
t500 = cos(qJ(2));
t458 = t500 * t473 - t496 * t474;
t488 = qJDD(1) + qJDD(2);
t455 = t488 * pkin(2) + t458;
t459 = t496 * t473 + t500 * t474;
t489 = qJD(1) + qJD(2);
t487 = t489 ^ 2;
t456 = -t487 * pkin(2) + t459;
t492 = sin(pkin(9));
t493 = cos(pkin(9));
t450 = t493 * t455 - t492 * t456;
t447 = t488 * pkin(3) + t450;
t451 = t492 * t455 + t493 * t456;
t448 = -t487 * pkin(3) + t451;
t495 = sin(qJ(4));
t499 = cos(qJ(4));
t444 = t495 * t447 + t499 * t448;
t483 = qJD(4) + t489;
t481 = t483 ^ 2;
t482 = qJDD(4) + t488;
t441 = -t481 * pkin(4) + t482 * pkin(8) + t444;
t491 = -g(1) + qJDD(3);
t494 = sin(qJ(5));
t498 = cos(qJ(5));
t438 = -t494 * t441 + t498 * t491;
t439 = t498 * t441 + t494 * t491;
t461 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t494 + Ifges(6,2) * t498) * t483;
t462 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t494 + Ifges(6,4) * t498) * t483;
t515 = qJD(5) * t483;
t466 = t494 * t482 + t498 * t515;
t467 = t498 * t482 - t494 * t515;
t519 = mrSges(6,1) * t438 - mrSges(6,2) * t439 + Ifges(6,5) * t466 + Ifges(6,6) * t467 + Ifges(6,3) * qJDD(5) + (t461 * t494 - t462 * t498) * t483;
t518 = t483 * t494;
t517 = t483 * t498;
t465 = (-mrSges(6,1) * t498 + mrSges(6,2) * t494) * t483;
t472 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t517;
t436 = m(6) * t438 + qJDD(5) * mrSges(6,1) - t466 * mrSges(6,3) + qJD(5) * t472 - t465 * t518;
t471 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t518;
t437 = m(6) * t439 - qJDD(5) * mrSges(6,2) + t467 * mrSges(6,3) - qJD(5) * t471 + t465 * t517;
t509 = -t494 * t436 + t498 * t437;
t422 = m(5) * t444 - t481 * mrSges(5,1) - t482 * mrSges(5,2) + t509;
t443 = t499 * t447 - t495 * t448;
t440 = -t482 * pkin(4) - t481 * pkin(8) - t443;
t506 = -m(6) * t440 + t467 * mrSges(6,1) - t466 * mrSges(6,2) - t471 * t518 + t472 * t517;
t431 = m(5) * t443 + t482 * mrSges(5,1) - t481 * mrSges(5,2) + t506;
t419 = t495 * t422 + t499 * t431;
t415 = m(4) * t450 + t488 * mrSges(4,1) - t487 * mrSges(4,2) + t419;
t510 = t499 * t422 - t495 * t431;
t416 = m(4) * t451 - t487 * mrSges(4,1) - t488 * mrSges(4,2) + t510;
t410 = t493 * t415 + t492 * t416;
t407 = m(3) * t458 + t488 * mrSges(3,1) - t487 * mrSges(3,2) + t410;
t511 = -t492 * t415 + t493 * t416;
t408 = m(3) * t459 - t487 * mrSges(3,1) - t488 * mrSges(3,2) + t511;
t512 = -t407 * t496 + t500 * t408;
t398 = m(2) * t476 - mrSges(2,1) * t502 - qJDD(1) * mrSges(2,2) + t512;
t401 = t500 * t407 + t496 * t408;
t399 = m(2) * t477 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t502 + t401;
t516 = t497 * t398 + t501 * t399;
t425 = t498 * t436 + t494 * t437;
t514 = m(5) * t491 + t425;
t513 = -t398 * t501 + t497 * t399;
t423 = m(4) * t491 + t514;
t460 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t494 + Ifges(6,6) * t498) * t483;
t428 = -mrSges(6,1) * t440 + mrSges(6,3) * t439 + Ifges(6,4) * t466 + Ifges(6,2) * t467 + Ifges(6,6) * qJDD(5) + qJD(5) * t462 - t460 * t518;
t429 = mrSges(6,2) * t440 - mrSges(6,3) * t438 + Ifges(6,1) * t466 + Ifges(6,4) * t467 + Ifges(6,5) * qJDD(5) - qJD(5) * t461 + t460 * t517;
t507 = mrSges(5,1) * t443 - mrSges(5,2) * t444 + Ifges(5,3) * t482 + pkin(4) * t506 + pkin(8) * t509 + t498 * t428 + t494 * t429;
t504 = mrSges(3,1) * t458 + mrSges(4,1) * t450 - mrSges(3,2) * t459 - mrSges(4,2) * t451 + pkin(2) * t410 + pkin(3) * t419 + t507 + (Ifges(3,3) + Ifges(4,3)) * t488;
t503 = mrSges(2,1) * t477 - mrSges(2,2) * t476 + Ifges(2,3) * qJDD(1) + pkin(1) * t401 + t504;
t417 = -mrSges(5,1) * t491 + mrSges(5,3) * t444 + t481 * Ifges(5,5) + Ifges(5,6) * t482 - pkin(4) * t425 - t519;
t411 = mrSges(5,2) * t491 - mrSges(5,3) * t443 + Ifges(5,5) * t482 - t481 * Ifges(5,6) - pkin(8) * t425 - t494 * t428 + t498 * t429;
t403 = mrSges(4,2) * t491 - mrSges(4,3) * t450 + Ifges(4,5) * t488 - t487 * Ifges(4,6) - pkin(7) * t419 + t499 * t411 - t495 * t417;
t402 = -mrSges(4,1) * t491 + mrSges(4,3) * t451 + t487 * Ifges(4,5) + Ifges(4,6) * t488 - pkin(3) * t514 + pkin(7) * t510 + t495 * t411 + t499 * t417;
t394 = -mrSges(3,2) * g(1) - mrSges(3,3) * t458 + Ifges(3,5) * t488 - Ifges(3,6) * t487 - qJ(3) * t410 - t402 * t492 + t403 * t493;
t393 = mrSges(3,1) * g(1) + mrSges(3,3) * t459 + t487 * Ifges(3,5) + Ifges(3,6) * t488 - pkin(2) * t423 + qJ(3) * t511 + t493 * t402 + t492 * t403;
t392 = -mrSges(2,2) * g(1) - mrSges(2,3) * t477 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t502 - pkin(6) * t401 - t393 * t496 + t394 * t500;
t391 = Ifges(2,6) * qJDD(1) + t502 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t476 + t496 * t394 + t500 * t393 - pkin(1) * (-m(3) * g(1) + t423) + pkin(6) * t512;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t423; -m(1) * g(2) + t516; -m(1) * g(3) + t513; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t503; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t513 + t501 * t391 + t497 * t392; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t516 + t497 * t391 - t501 * t392; t503; t504; t423; t507; t519;];
tauJB = t1;
