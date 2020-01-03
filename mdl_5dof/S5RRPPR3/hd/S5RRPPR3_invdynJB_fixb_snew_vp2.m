% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:33
% EndTime: 2019-12-31 19:26:35
% DurationCPUTime: 1.16s
% Computational Cost: add. (15987->179), mult. (19865->214), div. (0->0), fcn. (9638->8), ass. (0->80)
t494 = -pkin(3) - pkin(7);
t493 = mrSges(4,1) - mrSges(5,2);
t492 = -Ifges(5,4) + Ifges(4,5);
t491 = Ifges(5,5) - Ifges(4,6);
t462 = qJD(1) + qJD(2);
t468 = sin(qJ(5));
t490 = t462 * t468;
t471 = cos(qJ(5));
t489 = t462 * t471;
t470 = sin(qJ(1));
t473 = cos(qJ(1));
t449 = t470 * g(1) - t473 * g(2);
t444 = qJDD(1) * pkin(1) + t449;
t450 = -t473 * g(1) - t470 * g(2);
t474 = qJD(1) ^ 2;
t445 = -t474 * pkin(1) + t450;
t469 = sin(qJ(2));
t472 = cos(qJ(2));
t426 = t472 * t444 - t469 * t445;
t461 = qJDD(1) + qJDD(2);
t423 = t461 * pkin(2) + t426;
t427 = t469 * t444 + t472 * t445;
t460 = t462 ^ 2;
t424 = -t460 * pkin(2) + t427;
t466 = sin(pkin(8));
t467 = cos(pkin(8));
t418 = t467 * t423 - t466 * t424;
t480 = -t460 * qJ(4) + qJDD(4) - t418;
t413 = t494 * t461 + t480;
t465 = -g(3) + qJDD(3);
t409 = t471 * t413 - t468 * t465;
t438 = (mrSges(6,1) * t468 + mrSges(6,2) * t471) * t462;
t487 = qJD(5) * t462;
t440 = t471 * t461 - t468 * t487;
t446 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t490;
t406 = m(6) * t409 + qJDD(5) * mrSges(6,1) - t440 * mrSges(6,3) + qJD(5) * t446 - t438 * t489;
t410 = t468 * t413 + t471 * t465;
t439 = -t468 * t461 - t471 * t487;
t447 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t489;
t407 = m(6) * t410 - qJDD(5) * mrSges(6,2) + t439 * mrSges(6,3) - qJD(5) * t447 - t438 * t490;
t400 = t471 * t406 + t468 * t407;
t416 = -t461 * pkin(3) + t480;
t479 = -m(5) * t416 + t460 * mrSges(5,3) - t400;
t392 = m(4) * t418 - t460 * mrSges(4,2) + t493 * t461 + t479;
t419 = t466 * t423 + t467 * t424;
t481 = t461 * qJ(4) + 0.2e1 * qJD(4) * t462 + t419;
t414 = t460 * pkin(3) - t481;
t412 = t494 * t460 + t481;
t482 = -m(6) * t412 + t439 * mrSges(6,1) - t440 * mrSges(6,2) - t446 * t490 - t447 * t489;
t477 = -m(5) * t414 + t460 * mrSges(5,2) + t461 * mrSges(5,3) - t482;
t397 = m(4) * t419 - t460 * mrSges(4,1) - t461 * mrSges(4,2) + t477;
t390 = t467 * t392 + t466 * t397;
t387 = m(3) * t426 + t461 * mrSges(3,1) - t460 * mrSges(3,2) + t390;
t484 = -t466 * t392 + t467 * t397;
t388 = m(3) * t427 - t460 * mrSges(3,1) - t461 * mrSges(3,2) + t484;
t381 = t472 * t387 + t469 * t388;
t378 = m(2) * t449 + qJDD(1) * mrSges(2,1) - t474 * mrSges(2,2) + t381;
t485 = -t469 * t387 + t472 * t388;
t379 = m(2) * t450 - t474 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t485;
t488 = t473 * t378 + t470 * t379;
t486 = -t470 * t378 + t473 * t379;
t483 = -t468 * t406 + t471 * t407;
t399 = m(5) * t465 + t483;
t398 = m(4) * t465 + t399;
t431 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t471 - Ifges(6,2) * t468) * t462;
t432 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t471 - Ifges(6,4) * t468) * t462;
t478 = mrSges(6,1) * t409 - mrSges(6,2) * t410 + Ifges(6,5) * t440 + Ifges(6,6) * t439 + Ifges(6,3) * qJDD(5) + t431 * t489 + t432 * t490;
t394 = t461 * mrSges(5,2) - t479;
t430 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t471 - Ifges(6,6) * t468) * t462;
t402 = -mrSges(6,1) * t412 + mrSges(6,3) * t410 + Ifges(6,4) * t440 + Ifges(6,2) * t439 + Ifges(6,6) * qJDD(5) + qJD(5) * t432 - t430 * t489;
t403 = mrSges(6,2) * t412 - mrSges(6,3) * t409 + Ifges(6,1) * t440 + Ifges(6,4) * t439 + Ifges(6,5) * qJDD(5) - qJD(5) * t431 - t430 * t490;
t476 = mrSges(3,1) * t426 + mrSges(4,1) * t418 - mrSges(3,2) * t427 - mrSges(4,2) * t419 + mrSges(5,2) * t416 - mrSges(5,3) * t414 + pkin(2) * t390 - pkin(3) * t394 - pkin(7) * t400 + qJ(4) * t477 - t468 * t402 + t471 * t403 + (Ifges(4,3) + Ifges(3,3) + Ifges(5,1)) * t461;
t475 = mrSges(2,1) * t449 - mrSges(2,2) * t450 + Ifges(2,3) * qJDD(1) + pkin(1) * t381 + t476;
t383 = mrSges(5,1) * t416 - mrSges(4,3) * t418 + pkin(4) * t400 - qJ(4) * t399 + (mrSges(4,2) - mrSges(5,3)) * t465 + t492 * t461 + t491 * t460 + t478;
t382 = -mrSges(5,1) * t414 + mrSges(4,3) * t419 - pkin(3) * t399 - pkin(4) * t482 - pkin(7) * t483 - t471 * t402 - t468 * t403 + t492 * t460 - t491 * t461 - t493 * t465;
t374 = -mrSges(3,2) * g(3) - mrSges(3,3) * t426 + Ifges(3,5) * t461 - t460 * Ifges(3,6) - qJ(3) * t390 - t466 * t382 + t467 * t383;
t373 = mrSges(3,1) * g(3) + mrSges(3,3) * t427 + t460 * Ifges(3,5) + Ifges(3,6) * t461 - pkin(2) * t398 + qJ(3) * t484 + t467 * t382 + t466 * t383;
t372 = -mrSges(2,2) * g(3) - mrSges(2,3) * t449 + Ifges(2,5) * qJDD(1) - t474 * Ifges(2,6) - pkin(6) * t381 - t469 * t373 + t472 * t374;
t371 = Ifges(2,6) * qJDD(1) + t474 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t450 + t469 * t374 + t472 * t373 - pkin(1) * (-m(3) * g(3) + t398) + pkin(6) * t485;
t1 = [-m(1) * g(1) + t486; -m(1) * g(2) + t488; (-m(1) - m(2) - m(3)) * g(3) + t398; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t488 - t470 * t371 + t473 * t372; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t486 + t473 * t371 + t470 * t372; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t475; t475; t476; t398; t394; t478;];
tauJB = t1;
