% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:49
% EndTime: 2019-12-05 17:04:50
% DurationCPUTime: 0.88s
% Computational Cost: add. (12278->161), mult. (13805->202), div. (0->0), fcn. (7236->8), ass. (0->76)
t422 = sin(qJ(2));
t426 = cos(qJ(2));
t402 = t422 * g(1) - t426 * g(2);
t399 = qJDD(2) * pkin(2) + t402;
t403 = -t426 * g(1) - t422 * g(2);
t427 = qJD(2) ^ 2;
t400 = -t427 * pkin(2) + t403;
t421 = sin(qJ(3));
t425 = cos(qJ(3));
t386 = t425 * t399 - t421 * t400;
t415 = qJDD(2) + qJDD(3);
t383 = t415 * pkin(3) + t386;
t387 = t421 * t399 + t425 * t400;
t416 = qJD(2) + qJD(3);
t414 = t416 ^ 2;
t384 = -t414 * pkin(3) + t387;
t420 = sin(qJ(4));
t424 = cos(qJ(4));
t380 = t420 * t383 + t424 * t384;
t407 = qJDD(4) + t415;
t376 = t407 * pkin(6) + t380;
t418 = -g(3) + qJDD(1);
t419 = sin(qJ(5));
t423 = cos(qJ(5));
t374 = -t419 * t376 + t423 * t418;
t375 = t423 * t376 + t419 * t418;
t408 = qJD(4) + t416;
t389 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t419 + Ifges(6,2) * t423) * t408;
t390 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t419 + Ifges(6,4) * t423) * t408;
t442 = qJD(5) * t408;
t392 = t419 * t407 + t423 * t442;
t393 = t423 * t407 - t419 * t442;
t448 = mrSges(6,1) * t374 - mrSges(6,2) * t375 + Ifges(6,5) * t392 + Ifges(6,6) * t393 + Ifges(6,3) * qJDD(5) + (t389 * t419 - t390 * t423) * t408;
t447 = -m(1) - m(2);
t446 = t408 * t419;
t445 = t408 * t423;
t406 = t408 ^ 2;
t391 = (-mrSges(6,1) * t423 + mrSges(6,2) * t419) * t408;
t398 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t445;
t372 = m(6) * t374 + qJDD(5) * mrSges(6,1) - t392 * mrSges(6,3) + qJD(5) * t398 - t391 * t446;
t397 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t446;
t373 = m(6) * t375 - qJDD(5) * mrSges(6,2) + t393 * mrSges(6,3) - qJD(5) * t397 + t391 * t445;
t437 = -t419 * t372 + t423 * t373;
t360 = m(5) * t380 - t406 * mrSges(5,1) - t407 * mrSges(5,2) + t437;
t379 = t424 * t383 - t420 * t384;
t377 = -t406 * pkin(6) - t379;
t368 = m(5) * t379 - m(6) * t377 + t407 * mrSges(5,1) + t393 * mrSges(6,1) - t406 * mrSges(5,2) - t392 * mrSges(6,2) + (-t397 * t419 + t398 * t423) * t408;
t357 = t420 * t360 + t424 * t368;
t354 = m(4) * t386 + t415 * mrSges(4,1) - t414 * mrSges(4,2) + t357;
t438 = t424 * t360 - t420 * t368;
t355 = m(4) * t387 - t414 * mrSges(4,1) - t415 * mrSges(4,2) + t438;
t349 = t425 * t354 + t421 * t355;
t346 = m(3) * t402 + qJDD(2) * mrSges(3,1) - t427 * mrSges(3,2) + t349;
t439 = -t421 * t354 + t425 * t355;
t347 = m(3) * t403 - t427 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t439;
t444 = t426 * t346 + t422 * t347;
t443 = t423 * t372 + t419 * t373;
t441 = m(5) * t418 + t443;
t440 = -t422 * t346 + t426 * t347;
t436 = m(4) * t418 + t441;
t435 = qJ(1) * m(2) + mrSges(1,3) + mrSges(2,3);
t434 = m(3) * t418 + t436;
t432 = m(2) * t418 + t434;
t388 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t419 + Ifges(6,6) * t423) * t408;
t365 = -mrSges(6,1) * t377 + mrSges(6,3) * t375 + Ifges(6,4) * t392 + Ifges(6,2) * t393 + Ifges(6,6) * qJDD(5) + qJD(5) * t390 - t388 * t446;
t366 = mrSges(6,2) * t377 - mrSges(6,3) * t374 + Ifges(6,1) * t392 + Ifges(6,4) * t393 + Ifges(6,5) * qJDD(5) - qJD(5) * t389 + t388 * t445;
t431 = mrSges(5,1) * t379 - mrSges(5,2) * t380 + Ifges(5,3) * t407 + pkin(6) * t437 + t423 * t365 + t419 * t366;
t429 = mrSges(4,1) * t386 - mrSges(4,2) * t387 + Ifges(4,3) * t415 + pkin(3) * t357 + t431;
t428 = mrSges(3,1) * t402 - mrSges(3,2) * t403 + Ifges(3,3) * qJDD(2) + pkin(2) * t349 + t429;
t361 = -mrSges(5,1) * t418 + mrSges(5,3) * t380 + t406 * Ifges(5,5) + Ifges(5,6) * t407 - t448;
t350 = mrSges(5,2) * t418 - mrSges(5,3) * t379 + Ifges(5,5) * t407 - t406 * Ifges(5,6) - pkin(6) * t443 - t419 * t365 + t423 * t366;
t342 = mrSges(4,2) * t418 - mrSges(4,3) * t386 + Ifges(4,5) * t415 - t414 * Ifges(4,6) - pkin(5) * t357 + t424 * t350 - t420 * t361;
t341 = -mrSges(4,1) * t418 + mrSges(4,3) * t387 + t414 * Ifges(4,5) + Ifges(4,6) * t415 - pkin(3) * t441 + pkin(5) * t438 + t420 * t350 + t424 * t361;
t340 = mrSges(3,2) * t418 - mrSges(3,3) * t402 + Ifges(3,5) * qJDD(2) - t427 * Ifges(3,6) - pkin(4) * t349 - t421 * t341 + t425 * t342;
t339 = -mrSges(3,1) * t418 + mrSges(3,3) * t403 + t427 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t436 + pkin(4) * t439 + t425 * t341 + t421 * t342;
t1 = [t447 * g(1) + t440; t447 * g(2) + t444; -m(1) * g(3) + t432; -mrSges(1,2) * g(3) + mrSges(2,2) * t418 + t435 * g(2) - qJ(1) * t444 - t422 * t339 + t426 * t340; mrSges(1,1) * g(3) - mrSges(2,1) * t418 - pkin(1) * t434 - t435 * g(1) + qJ(1) * t440 + t426 * t339 + t422 * t340; t428 + (-mrSges(1,1) - mrSges(2,1)) * g(2) + pkin(1) * t444 + (mrSges(2,2) + mrSges(1,2)) * g(1); t432; t428; t429; t431; t448;];
tauJB = t1;
