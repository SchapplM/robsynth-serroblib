% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:47
% EndTime: 2019-12-31 18:17:48
% DurationCPUTime: 1.11s
% Computational Cost: add. (14737->178), mult. (19585->214), div. (0->0), fcn. (9500->8), ass. (0->81)
t497 = -pkin(3) - pkin(7);
t496 = mrSges(4,1) - mrSges(5,2);
t495 = -Ifges(5,4) + Ifges(4,5);
t494 = Ifges(5,5) - Ifges(4,6);
t463 = qJD(1) + qJD(3);
t470 = sin(qJ(5));
t493 = t463 * t470;
t473 = cos(qJ(5));
t492 = t463 * t473;
t472 = sin(qJ(1));
t475 = cos(qJ(1));
t450 = t472 * g(1) - t475 * g(2);
t445 = qJDD(1) * pkin(1) + t450;
t451 = -t475 * g(1) - t472 * g(2);
t476 = qJD(1) ^ 2;
t446 = -t476 * pkin(1) + t451;
t468 = sin(pkin(8));
t469 = cos(pkin(8));
t427 = t469 * t445 - t468 * t446;
t424 = qJDD(1) * pkin(2) + t427;
t428 = t468 * t445 + t469 * t446;
t425 = -t476 * pkin(2) + t428;
t471 = sin(qJ(3));
t474 = cos(qJ(3));
t419 = t474 * t424 - t471 * t425;
t461 = t463 ^ 2;
t462 = qJDD(1) + qJDD(3);
t482 = -t461 * qJ(4) + qJDD(4) - t419;
t414 = t462 * t497 + t482;
t467 = -g(3) + qJDD(2);
t410 = t473 * t414 - t470 * t467;
t439 = (mrSges(6,1) * t470 + mrSges(6,2) * t473) * t463;
t490 = qJD(5) * t463;
t441 = t473 * t462 - t470 * t490;
t447 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t493;
t407 = m(6) * t410 + qJDD(5) * mrSges(6,1) - t441 * mrSges(6,3) + qJD(5) * t447 - t439 * t492;
t411 = t470 * t414 + t473 * t467;
t440 = -t470 * t462 - t473 * t490;
t448 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t492;
t408 = m(6) * t411 - qJDD(5) * mrSges(6,2) + t440 * mrSges(6,3) - qJD(5) * t448 - t439 * t493;
t401 = t473 * t407 + t470 * t408;
t417 = -t462 * pkin(3) + t482;
t481 = -m(5) * t417 + t461 * mrSges(5,3) - t401;
t393 = m(4) * t419 - t461 * mrSges(4,2) + t462 * t496 + t481;
t420 = t471 * t424 + t474 * t425;
t483 = t462 * qJ(4) + 0.2e1 * qJD(4) * t463 + t420;
t415 = t461 * pkin(3) - t483;
t413 = t461 * t497 + t483;
t484 = -m(6) * t413 + t440 * mrSges(6,1) - t441 * mrSges(6,2) - t447 * t493 - t448 * t492;
t479 = -m(5) * t415 + t461 * mrSges(5,2) + t462 * mrSges(5,3) - t484;
t398 = m(4) * t420 - t461 * mrSges(4,1) - t462 * mrSges(4,2) + t479;
t391 = t474 * t393 + t471 * t398;
t388 = m(3) * t427 + qJDD(1) * mrSges(3,1) - t476 * mrSges(3,2) + t391;
t487 = -t471 * t393 + t474 * t398;
t389 = m(3) * t428 - t476 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t487;
t382 = t469 * t388 + t468 * t389;
t379 = m(2) * t450 + qJDD(1) * mrSges(2,1) - t476 * mrSges(2,2) + t382;
t488 = -t468 * t388 + t469 * t389;
t380 = m(2) * t451 - t476 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t488;
t491 = t475 * t379 + t472 * t380;
t489 = -t472 * t379 + t475 * t380;
t486 = -t470 * t407 + t473 * t408;
t400 = m(5) * t467 + t486;
t485 = m(4) * t467 + t400;
t399 = m(3) * t467 + t485;
t432 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t473 - Ifges(6,2) * t470) * t463;
t433 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t473 - Ifges(6,4) * t470) * t463;
t480 = mrSges(6,1) * t410 - mrSges(6,2) * t411 + Ifges(6,5) * t441 + Ifges(6,6) * t440 + Ifges(6,3) * qJDD(5) + t432 * t492 + t433 * t493;
t395 = t462 * mrSges(5,2) - t481;
t431 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t473 - Ifges(6,6) * t470) * t463;
t403 = -mrSges(6,1) * t413 + mrSges(6,3) * t411 + Ifges(6,4) * t441 + Ifges(6,2) * t440 + Ifges(6,6) * qJDD(5) + qJD(5) * t433 - t431 * t492;
t404 = mrSges(6,2) * t413 - mrSges(6,3) * t410 + Ifges(6,1) * t441 + Ifges(6,4) * t440 + Ifges(6,5) * qJDD(5) - qJD(5) * t432 - t431 * t493;
t478 = mrSges(4,1) * t419 - mrSges(4,2) * t420 + mrSges(5,2) * t417 - mrSges(5,3) * t415 - pkin(3) * t395 - pkin(7) * t401 + qJ(4) * t479 - t470 * t403 + t473 * t404 + (Ifges(4,3) + Ifges(5,1)) * t462;
t477 = mrSges(2,1) * t450 + mrSges(3,1) * t427 - mrSges(2,2) * t451 - mrSges(3,2) * t428 + pkin(1) * t382 + pkin(2) * t391 + t478 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t384 = mrSges(5,1) * t417 - mrSges(4,3) * t419 + pkin(4) * t401 - qJ(4) * t400 + (mrSges(4,2) - mrSges(5,3)) * t467 + t495 * t462 + t494 * t461 + t480;
t383 = -mrSges(5,1) * t415 + mrSges(4,3) * t420 - pkin(3) * t400 - pkin(4) * t484 - pkin(7) * t486 - t473 * t403 - t470 * t404 + t495 * t461 - t494 * t462 - t496 * t467;
t375 = mrSges(3,2) * t467 - mrSges(3,3) * t427 + Ifges(3,5) * qJDD(1) - t476 * Ifges(3,6) - pkin(6) * t391 - t471 * t383 + t474 * t384;
t374 = -mrSges(3,1) * t467 + mrSges(3,3) * t428 + t476 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t485 + pkin(6) * t487 + t474 * t383 + t471 * t384;
t373 = -mrSges(2,2) * g(3) - mrSges(2,3) * t450 + Ifges(2,5) * qJDD(1) - t476 * Ifges(2,6) - qJ(2) * t382 - t468 * t374 + t469 * t375;
t372 = mrSges(2,1) * g(3) + mrSges(2,3) * t451 + t476 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t399 + qJ(2) * t488 + t469 * t374 + t468 * t375;
t1 = [-m(1) * g(1) + t489; -m(1) * g(2) + t491; (-m(1) - m(2)) * g(3) + t399; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t491 - t472 * t372 + t475 * t373; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t489 + t475 * t372 + t472 * t373; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t477; t477; t399; t478; t395; t480;];
tauJB = t1;
