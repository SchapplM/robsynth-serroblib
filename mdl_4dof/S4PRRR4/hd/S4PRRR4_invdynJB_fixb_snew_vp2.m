% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRRR4
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR4_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:35
% EndTime: 2019-12-31 16:32:37
% DurationCPUTime: 1.18s
% Computational Cost: add. (10721->186), mult. (21197->244), div. (0->0), fcn. (12879->8), ass. (0->82)
t446 = sin(pkin(7));
t447 = cos(pkin(7));
t431 = t446 * g(1) - t447 * g(2);
t432 = -t447 * g(1) - t446 * g(2);
t450 = sin(qJ(2));
t453 = cos(qJ(2));
t415 = t450 * t431 + t453 * t432;
t454 = qJD(2) ^ 2;
t412 = -t454 * pkin(2) + qJDD(2) * pkin(5) + t415;
t445 = -g(3) + qJDD(1);
t449 = sin(qJ(3));
t452 = cos(qJ(3));
t403 = -t449 * t412 + t452 * t445;
t469 = qJD(2) * qJD(3);
t467 = t452 * t469;
t429 = t449 * qJDD(2) + t467;
t394 = (-t429 + t467) * pkin(6) + (t449 * t452 * t454 + qJDD(3)) * pkin(3) + t403;
t404 = t452 * t412 + t449 * t445;
t430 = t452 * qJDD(2) - t449 * t469;
t471 = qJD(2) * t449;
t435 = qJD(3) * pkin(3) - pkin(6) * t471;
t444 = t452 ^ 2;
t395 = -t444 * t454 * pkin(3) + t430 * pkin(6) - qJD(3) * t435 + t404;
t448 = sin(qJ(4));
t451 = cos(qJ(4));
t392 = t451 * t394 - t448 * t395;
t421 = (-t448 * t449 + t451 * t452) * qJD(2);
t402 = t421 * qJD(4) + t451 * t429 + t448 * t430;
t422 = (t448 * t452 + t449 * t451) * qJD(2);
t410 = -t421 * mrSges(5,1) + t422 * mrSges(5,2);
t442 = qJD(3) + qJD(4);
t416 = -t442 * mrSges(5,2) + t421 * mrSges(5,3);
t441 = qJDD(3) + qJDD(4);
t389 = m(5) * t392 + t441 * mrSges(5,1) - t402 * mrSges(5,3) - t422 * t410 + t442 * t416;
t393 = t448 * t394 + t451 * t395;
t401 = -t422 * qJD(4) - t448 * t429 + t451 * t430;
t417 = t442 * mrSges(5,1) - t422 * mrSges(5,3);
t390 = m(5) * t393 - t441 * mrSges(5,2) + t401 * mrSges(5,3) + t421 * t410 - t442 * t417;
t380 = t451 * t389 + t448 * t390;
t419 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t449 + Ifges(4,2) * t452) * qJD(2);
t420 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t449 + Ifges(4,4) * t452) * qJD(2);
t406 = Ifges(5,4) * t422 + Ifges(5,2) * t421 + Ifges(5,6) * t442;
t407 = Ifges(5,1) * t422 + Ifges(5,4) * t421 + Ifges(5,5) * t442;
t457 = -mrSges(5,1) * t392 + mrSges(5,2) * t393 - Ifges(5,5) * t402 - Ifges(5,6) * t401 - Ifges(5,3) * t441 - t422 * t406 + t421 * t407;
t473 = mrSges(4,1) * t403 - mrSges(4,2) * t404 + Ifges(4,5) * t429 + Ifges(4,6) * t430 + Ifges(4,3) * qJDD(3) + pkin(3) * t380 + (t449 * t419 - t452 * t420) * qJD(2) - t457;
t428 = (-mrSges(4,1) * t452 + mrSges(4,2) * t449) * qJD(2);
t470 = qJD(2) * t452;
t434 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t470;
t378 = m(4) * t403 + qJDD(3) * mrSges(4,1) - t429 * mrSges(4,3) + qJD(3) * t434 - t428 * t471 + t380;
t433 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t471;
t463 = -t448 * t389 + t451 * t390;
t379 = m(4) * t404 - qJDD(3) * mrSges(4,2) + t430 * mrSges(4,3) - qJD(3) * t433 + t428 * t470 + t463;
t464 = -t449 * t378 + t452 * t379;
t372 = m(3) * t415 - t454 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t464;
t414 = t453 * t431 - t450 * t432;
t460 = -qJDD(2) * pkin(2) - t414;
t411 = -t454 * pkin(5) + t460;
t396 = t435 * t471 - t430 * pkin(3) + (-pkin(6) * t444 - pkin(5)) * t454 + t460;
t458 = m(5) * t396 - t401 * mrSges(5,1) + t402 * mrSges(5,2) - t421 * t416 + t422 * t417;
t456 = -m(4) * t411 + t430 * mrSges(4,1) - t429 * mrSges(4,2) - t433 * t471 + t434 * t470 - t458;
t384 = m(3) * t414 + qJDD(2) * mrSges(3,1) - t454 * mrSges(3,2) + t456;
t367 = t450 * t372 + t453 * t384;
t365 = m(2) * t431 + t367;
t465 = t453 * t372 - t450 * t384;
t366 = m(2) * t432 + t465;
t472 = t447 * t365 + t446 * t366;
t374 = t452 * t378 + t449 * t379;
t468 = m(3) * t445 + t374;
t466 = -t446 * t365 + t447 * t366;
t462 = m(2) * t445 + t468;
t405 = Ifges(5,5) * t422 + Ifges(5,6) * t421 + Ifges(5,3) * t442;
t381 = -mrSges(5,1) * t396 + mrSges(5,3) * t393 + Ifges(5,4) * t402 + Ifges(5,2) * t401 + Ifges(5,6) * t441 - t422 * t405 + t442 * t407;
t382 = mrSges(5,2) * t396 - mrSges(5,3) * t392 + Ifges(5,1) * t402 + Ifges(5,4) * t401 + Ifges(5,5) * t441 + t421 * t405 - t442 * t406;
t418 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t449 + Ifges(4,6) * t452) * qJD(2);
t361 = -mrSges(4,1) * t411 + mrSges(4,3) * t404 + Ifges(4,4) * t429 + Ifges(4,2) * t430 + Ifges(4,6) * qJDD(3) - pkin(3) * t458 + pkin(6) * t463 + qJD(3) * t420 + t451 * t381 + t448 * t382 - t418 * t471;
t369 = mrSges(4,2) * t411 - mrSges(4,3) * t403 + Ifges(4,1) * t429 + Ifges(4,4) * t430 + Ifges(4,5) * qJDD(3) - pkin(6) * t380 - qJD(3) * t419 - t448 * t381 + t451 * t382 + t418 * t470;
t459 = mrSges(3,1) * t414 - mrSges(3,2) * t415 + Ifges(3,3) * qJDD(2) + pkin(2) * t456 + pkin(5) * t464 + t452 * t361 + t449 * t369;
t359 = -mrSges(3,1) * t445 + mrSges(3,3) * t415 + t454 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t374 - t473;
t358 = mrSges(3,2) * t445 - mrSges(3,3) * t414 + Ifges(3,5) * qJDD(2) - t454 * Ifges(3,6) - pkin(5) * t374 - t449 * t361 + t452 * t369;
t357 = mrSges(2,2) * t445 - mrSges(2,3) * t431 - pkin(4) * t367 + t453 * t358 - t450 * t359;
t356 = -mrSges(2,1) * t445 + mrSges(2,3) * t432 - pkin(1) * t468 + pkin(4) * t465 + t450 * t358 + t453 * t359;
t1 = [-m(1) * g(1) + t466; -m(1) * g(2) + t472; -m(1) * g(3) + t462; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t472 - t446 * t356 + t447 * t357; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t466 + t447 * t356 + t446 * t357; -mrSges(1,1) * g(2) + mrSges(2,1) * t431 + mrSges(1,2) * g(1) - mrSges(2,2) * t432 + pkin(1) * t367 + t459; t462; t459; t473; -t457;];
tauJB = t1;
