% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR3_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR3_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:54
% EndTime: 2019-12-31 16:37:55
% DurationCPUTime: 1.01s
% Computational Cost: add. (9822->174), mult. (20493->222), div. (0->0), fcn. (11927->8), ass. (0->82)
t474 = qJD(1) ^ 2;
t468 = cos(pkin(7));
t497 = pkin(3) * t468;
t466 = sin(pkin(7));
t496 = mrSges(4,2) * t466;
t462 = t468 ^ 2;
t495 = t462 * t474;
t471 = sin(qJ(1));
t473 = cos(qJ(1));
t450 = t471 * g(1) - t473 * g(2);
t447 = qJDD(1) * pkin(1) + t450;
t451 = -t473 * g(1) - t471 * g(2);
t448 = -t474 * pkin(1) + t451;
t467 = sin(pkin(6));
t469 = cos(pkin(6));
t437 = t467 * t447 + t469 * t448;
t429 = -t474 * pkin(2) + qJDD(1) * qJ(3) + t437;
t465 = -g(3) + qJDD(2);
t491 = qJD(1) * qJD(3);
t493 = t468 * t465 - 0.2e1 * t466 * t491;
t419 = (-pkin(5) * qJDD(1) + t474 * t497 - t429) * t466 + t493;
t423 = t466 * t465 + (t429 + 0.2e1 * t491) * t468;
t490 = qJDD(1) * t468;
t420 = -pkin(3) * t495 + pkin(5) * t490 + t423;
t470 = sin(qJ(4));
t472 = cos(qJ(4));
t417 = t472 * t419 - t470 * t420;
t480 = -t466 * t470 + t468 * t472;
t440 = t480 * qJD(1);
t481 = t466 * t472 + t468 * t470;
t441 = t481 * qJD(1);
t431 = -t440 * mrSges(5,1) + t441 * mrSges(5,2);
t434 = t440 * qJD(4) + t481 * qJDD(1);
t438 = -qJD(4) * mrSges(5,2) + t440 * mrSges(5,3);
t415 = m(5) * t417 + qJDD(4) * mrSges(5,1) - t434 * mrSges(5,3) + qJD(4) * t438 - t441 * t431;
t418 = t470 * t419 + t472 * t420;
t433 = -t441 * qJD(4) + t480 * qJDD(1);
t439 = qJD(4) * mrSges(5,1) - t441 * mrSges(5,3);
t416 = m(5) * t418 - qJDD(4) * mrSges(5,2) + t433 * mrSges(5,3) - qJD(4) * t439 + t440 * t431;
t405 = t472 * t415 + t470 * t416;
t422 = -t466 * t429 + t493;
t479 = mrSges(4,3) * qJDD(1) + t474 * (-mrSges(4,1) * t468 + t496);
t403 = m(4) * t422 - t479 * t466 + t405;
t486 = -t470 * t415 + t472 * t416;
t404 = m(4) * t423 + t479 * t468 + t486;
t487 = -t466 * t403 + t468 * t404;
t396 = m(3) * t437 - t474 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t487;
t436 = t469 * t447 - t467 * t448;
t482 = qJDD(3) - t436;
t425 = -qJDD(1) * pkin(2) - t474 * qJ(3) + t482;
t461 = t466 ^ 2;
t421 = (-pkin(2) - t497) * qJDD(1) + (-qJ(3) + (-t461 - t462) * pkin(5)) * t474 + t482;
t478 = m(5) * t421 - t433 * mrSges(5,1) + t434 * mrSges(5,2) - t440 * t438 + t441 * t439;
t476 = -m(4) * t425 + mrSges(4,1) * t490 - t478 + (t461 * t474 + t495) * mrSges(4,3);
t409 = m(3) * t436 - t474 * mrSges(3,2) + (mrSges(3,1) - t496) * qJDD(1) + t476;
t391 = t467 * t396 + t469 * t409;
t388 = m(2) * t450 + qJDD(1) * mrSges(2,1) - t474 * mrSges(2,2) + t391;
t488 = t469 * t396 - t467 * t409;
t389 = m(2) * t451 - t474 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t488;
t494 = t473 * t388 + t471 * t389;
t399 = t468 * t403 + t466 * t404;
t483 = Ifges(4,5) * t466 + Ifges(4,6) * t468;
t492 = t474 * t483;
t397 = m(3) * t465 + t399;
t489 = -t471 * t388 + t473 * t389;
t485 = Ifges(4,1) * t466 + Ifges(4,4) * t468;
t484 = Ifges(4,4) * t466 + Ifges(4,2) * t468;
t426 = Ifges(5,5) * t441 + Ifges(5,6) * t440 + Ifges(5,3) * qJD(4);
t428 = Ifges(5,1) * t441 + Ifges(5,4) * t440 + Ifges(5,5) * qJD(4);
t406 = -mrSges(5,1) * t421 + mrSges(5,3) * t418 + Ifges(5,4) * t434 + Ifges(5,2) * t433 + Ifges(5,6) * qJDD(4) + qJD(4) * t428 - t441 * t426;
t427 = Ifges(5,4) * t441 + Ifges(5,2) * t440 + Ifges(5,6) * qJD(4);
t407 = mrSges(5,2) * t421 - mrSges(5,3) * t417 + Ifges(5,1) * t434 + Ifges(5,4) * t433 + Ifges(5,5) * qJDD(4) - qJD(4) * t427 + t440 * t426;
t384 = -mrSges(4,1) * t425 + mrSges(4,3) * t423 - pkin(3) * t478 + pkin(5) * t486 + t484 * qJDD(1) + t472 * t406 + t470 * t407 - t466 * t492;
t393 = mrSges(4,2) * t425 - mrSges(4,3) * t422 - pkin(5) * t405 + t485 * qJDD(1) - t470 * t406 + t472 * t407 + t468 * t492;
t411 = qJDD(1) * t496 - t476;
t477 = mrSges(2,1) * t450 + mrSges(3,1) * t436 - mrSges(2,2) * t451 - mrSges(3,2) * t437 + pkin(1) * t391 - pkin(2) * t411 + qJ(3) * t487 + t468 * t384 + t466 * t393 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t475 = mrSges(5,1) * t417 - mrSges(5,2) * t418 + Ifges(5,5) * t434 + Ifges(5,6) * t433 + Ifges(5,3) * qJDD(4) + t441 * t427 - t440 * t428;
t382 = -mrSges(3,1) * t465 - mrSges(4,1) * t422 + mrSges(4,2) * t423 + mrSges(3,3) * t437 - pkin(2) * t399 - pkin(3) * t405 + (Ifges(3,6) - t483) * qJDD(1) - t475 + (-t466 * t484 + t468 * t485 + Ifges(3,5)) * t474;
t381 = mrSges(3,2) * t465 - mrSges(3,3) * t436 + Ifges(3,5) * qJDD(1) - t474 * Ifges(3,6) - qJ(3) * t399 - t466 * t384 + t468 * t393;
t380 = -mrSges(2,2) * g(3) - mrSges(2,3) * t450 + Ifges(2,5) * qJDD(1) - t474 * Ifges(2,6) - qJ(2) * t391 + t469 * t381 - t467 * t382;
t379 = mrSges(2,1) * g(3) + mrSges(2,3) * t451 + t474 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t397 + qJ(2) * t488 + t467 * t381 + t469 * t382;
t1 = [-m(1) * g(1) + t489; -m(1) * g(2) + t494; (-m(1) - m(2)) * g(3) + t397; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t494 - t471 * t379 + t473 * t380; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t489 + t473 * t379 + t471 * t380; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t477; t477; t397; t411; t475;];
tauJB = t1;
