% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRR4
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:21
% EndTime: 2019-12-31 16:50:22
% DurationCPUTime: 1.19s
% Computational Cost: add. (10944->196), mult. (20435->249), div. (0->0), fcn. (11231->8), ass. (0->84)
t447 = sin(qJ(1));
t450 = cos(qJ(1));
t434 = t447 * g(1) - t450 * g(2);
t425 = qJDD(1) * pkin(1) + t434;
t435 = -t450 * g(1) - t447 * g(2);
t452 = qJD(1) ^ 2;
t427 = -t452 * pkin(1) + t435;
t443 = sin(pkin(7));
t444 = cos(pkin(7));
t409 = t444 * t425 - t443 * t427;
t399 = -qJDD(1) * pkin(2) - t452 * pkin(5) - t409;
t446 = sin(qJ(3));
t449 = cos(qJ(3));
t463 = qJD(1) * qJD(3);
t461 = t449 * t463;
t429 = t446 * qJDD(1) + t461;
t462 = t446 * t463;
t430 = t449 * qJDD(1) - t462;
t393 = (-t429 - t461) * pkin(6) + (-t430 + t462) * pkin(3) + t399;
t410 = t443 * t425 + t444 * t427;
t400 = -t452 * pkin(2) + qJDD(1) * pkin(5) + t410;
t442 = -g(3) + qJDD(2);
t397 = t449 * t400 + t446 * t442;
t428 = (-pkin(3) * t449 - pkin(6) * t446) * qJD(1);
t451 = qJD(3) ^ 2;
t464 = t449 * qJD(1);
t395 = -t451 * pkin(3) + qJDD(3) * pkin(6) + t428 * t464 + t397;
t445 = sin(qJ(4));
t448 = cos(qJ(4));
t391 = t448 * t393 - t445 * t395;
t465 = qJD(1) * t446;
t423 = t448 * qJD(3) - t445 * t465;
t407 = t423 * qJD(4) + t445 * qJDD(3) + t448 * t429;
t424 = t445 * qJD(3) + t448 * t465;
t411 = -t423 * mrSges(5,1) + t424 * mrSges(5,2);
t436 = qJD(4) - t464;
t412 = -t436 * mrSges(5,2) + t423 * mrSges(5,3);
t422 = qJDD(4) - t430;
t388 = m(5) * t391 + t422 * mrSges(5,1) - t407 * mrSges(5,3) - t424 * t411 + t436 * t412;
t392 = t445 * t393 + t448 * t395;
t406 = -t424 * qJD(4) + t448 * qJDD(3) - t445 * t429;
t413 = t436 * mrSges(5,1) - t424 * mrSges(5,3);
t389 = m(5) * t392 - t422 * mrSges(5,2) + t406 * mrSges(5,3) + t423 * t411 - t436 * t413;
t382 = -t445 * t388 + t448 * t389;
t467 = t449 * t442;
t394 = -qJDD(3) * pkin(3) - t451 * pkin(6) - t467 + (qJD(1) * t428 + t400) * t446;
t401 = Ifges(5,5) * t424 + Ifges(5,6) * t423 + Ifges(5,3) * t436;
t403 = Ifges(5,1) * t424 + Ifges(5,4) * t423 + Ifges(5,5) * t436;
t383 = -mrSges(5,1) * t394 + mrSges(5,3) * t392 + Ifges(5,4) * t407 + Ifges(5,2) * t406 + Ifges(5,6) * t422 - t424 * t401 + t436 * t403;
t402 = Ifges(5,4) * t424 + Ifges(5,2) * t423 + Ifges(5,6) * t436;
t384 = mrSges(5,2) * t394 - mrSges(5,3) * t391 + Ifges(5,1) * t407 + Ifges(5,4) * t406 + Ifges(5,5) * t422 + t423 * t401 - t436 * t402;
t390 = -m(5) * t394 + t406 * mrSges(5,1) - t407 * mrSges(5,2) + t423 * t412 - t424 * t413;
t396 = -t446 * t400 + t467;
t418 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t446 + Ifges(4,2) * t449) * qJD(1);
t419 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t446 + Ifges(4,4) * t449) * qJD(1);
t468 = mrSges(4,1) * t396 - mrSges(4,2) * t397 + Ifges(4,5) * t429 + Ifges(4,6) * t430 + Ifges(4,3) * qJDD(3) + pkin(3) * t390 + pkin(6) * t382 + t448 * t383 + t445 * t384 + (t446 * t418 - t449 * t419) * qJD(1);
t426 = (-mrSges(4,1) * t449 + mrSges(4,2) * t446) * qJD(1);
t432 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t465;
t380 = m(4) * t397 - qJDD(3) * mrSges(4,2) + t430 * mrSges(4,3) - qJD(3) * t432 + t426 * t464 + t382;
t433 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t464;
t386 = m(4) * t396 + qJDD(3) * mrSges(4,1) - t429 * mrSges(4,3) + qJD(3) * t433 - t426 * t465 + t390;
t458 = t449 * t380 - t446 * t386;
t371 = m(3) * t410 - t452 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t458;
t381 = t448 * t388 + t445 * t389;
t455 = -m(4) * t399 + t430 * mrSges(4,1) - t429 * mrSges(4,2) - t432 * t465 + t433 * t464 - t381;
t376 = m(3) * t409 + qJDD(1) * mrSges(3,1) - t452 * mrSges(3,2) + t455;
t364 = t443 * t371 + t444 * t376;
t361 = m(2) * t434 + qJDD(1) * mrSges(2,1) - t452 * mrSges(2,2) + t364;
t459 = t444 * t371 - t443 * t376;
t362 = m(2) * t435 - t452 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t459;
t466 = t450 * t361 + t447 * t362;
t374 = t446 * t380 + t449 * t386;
t372 = m(3) * t442 + t374;
t460 = -t447 * t361 + t450 * t362;
t417 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t446 + Ifges(4,6) * t449) * qJD(1);
t366 = mrSges(4,2) * t399 - mrSges(4,3) * t396 + Ifges(4,1) * t429 + Ifges(4,4) * t430 + Ifges(4,5) * qJDD(3) - pkin(6) * t381 - qJD(3) * t418 - t445 * t383 + t448 * t384 + t417 * t464;
t454 = mrSges(5,1) * t391 - mrSges(5,2) * t392 + Ifges(5,5) * t407 + Ifges(5,6) * t406 + Ifges(5,3) * t422 + t424 * t402 - t423 * t403;
t368 = -mrSges(4,1) * t399 + mrSges(4,3) * t397 + Ifges(4,4) * t429 + Ifges(4,2) * t430 + Ifges(4,6) * qJDD(3) - pkin(3) * t381 + qJD(3) * t419 - t417 * t465 - t454;
t456 = mrSges(2,1) * t434 + mrSges(3,1) * t409 - mrSges(2,2) * t435 - mrSges(3,2) * t410 + pkin(1) * t364 + pkin(2) * t455 + pkin(5) * t458 + t446 * t366 + t449 * t368 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t357 = -mrSges(3,1) * t442 + mrSges(3,3) * t410 + t452 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t374 - t468;
t356 = mrSges(3,2) * t442 - mrSges(3,3) * t409 + Ifges(3,5) * qJDD(1) - t452 * Ifges(3,6) - pkin(5) * t374 + t449 * t366 - t446 * t368;
t355 = -mrSges(2,2) * g(3) - mrSges(2,3) * t434 + Ifges(2,5) * qJDD(1) - t452 * Ifges(2,6) - qJ(2) * t364 + t444 * t356 - t443 * t357;
t354 = mrSges(2,1) * g(3) + mrSges(2,3) * t435 + t452 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t372 + qJ(2) * t459 + t443 * t356 + t444 * t357;
t1 = [-m(1) * g(1) + t460; -m(1) * g(2) + t466; (-m(1) - m(2)) * g(3) + t372; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t466 - t447 * t354 + t450 * t355; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t460 + t450 * t354 + t447 * t355; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t456; t456; t372; t468; t454;];
tauJB = t1;
