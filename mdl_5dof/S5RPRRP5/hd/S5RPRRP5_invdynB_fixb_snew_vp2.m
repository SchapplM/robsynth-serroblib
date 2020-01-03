% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:49
% EndTime: 2019-12-31 18:40:50
% DurationCPUTime: 1.60s
% Computational Cost: add. (18784->204), mult. (25449->249), div. (0->0), fcn. (12352->8), ass. (0->85)
t454 = Ifges(5,1) + Ifges(6,1);
t450 = Ifges(5,4) - Ifges(6,5);
t449 = Ifges(5,5) + Ifges(6,4);
t453 = Ifges(5,2) + Ifges(6,3);
t448 = Ifges(5,6) - Ifges(6,6);
t452 = Ifges(5,3) + Ifges(6,2);
t451 = mrSges(5,3) + mrSges(6,2);
t417 = qJD(1) + qJD(3);
t423 = sin(qJ(4));
t447 = t417 * t423;
t426 = cos(qJ(4));
t446 = t417 * t426;
t420 = -g(3) + qJDD(2);
t445 = t426 * t420;
t425 = sin(qJ(1));
t428 = cos(qJ(1));
t408 = t425 * g(1) - t428 * g(2);
t402 = qJDD(1) * pkin(1) + t408;
t409 = -t428 * g(1) - t425 * g(2);
t430 = qJD(1) ^ 2;
t403 = -t430 * pkin(1) + t409;
t421 = sin(pkin(8));
t422 = cos(pkin(8));
t378 = t422 * t402 - t421 * t403;
t376 = qJDD(1) * pkin(2) + t378;
t379 = t421 * t402 + t422 * t403;
t377 = -t430 * pkin(2) + t379;
t424 = sin(qJ(3));
t427 = cos(qJ(3));
t372 = t424 * t376 + t427 * t377;
t415 = t417 ^ 2;
t416 = qJDD(1) + qJDD(3);
t370 = -t415 * pkin(3) + t416 * pkin(7) + t372;
t367 = t426 * t370 + t423 * t420;
t394 = (-mrSges(5,1) * t426 + mrSges(5,2) * t423) * t417;
t440 = qJD(4) * t417;
t396 = t426 * t416 - t423 * t440;
t404 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t447;
t392 = (-pkin(4) * t426 - qJ(5) * t423) * t417;
t429 = qJD(4) ^ 2;
t364 = -t429 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t392 * t446 + t367;
t393 = (-mrSges(6,1) * t426 - mrSges(6,3) * t423) * t417;
t405 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t447;
t434 = m(6) * t364 + qJDD(4) * mrSges(6,3) + qJD(4) * t405 + t393 * t446;
t359 = m(5) * t367 - qJDD(4) * mrSges(5,2) - qJD(4) * t404 + t394 * t446 + t451 * t396 + t434;
t366 = -t423 * t370 + t445;
t395 = t423 * t416 + t426 * t440;
t406 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t446;
t365 = -qJDD(4) * pkin(4) - t429 * qJ(5) - t445 + qJDD(5) + (t392 * t417 + t370) * t423;
t407 = mrSges(6,2) * t446 + qJD(4) * mrSges(6,3);
t432 = -m(6) * t365 + qJDD(4) * mrSges(6,1) + qJD(4) * t407;
t360 = m(5) * t366 + qJDD(4) * mrSges(5,1) + qJD(4) * t406 + (-t393 - t394) * t447 - t451 * t395 + t432;
t435 = t426 * t359 - t423 * t360;
t352 = m(4) * t372 - t415 * mrSges(4,1) - t416 * mrSges(4,2) + t435;
t371 = t427 * t376 - t424 * t377;
t369 = -t416 * pkin(3) - t415 * pkin(7) - t371;
t362 = -t396 * pkin(4) - t395 * qJ(5) + (-0.2e1 * qJD(5) * t423 + (pkin(4) * t423 - qJ(5) * t426) * qJD(4)) * t417 + t369;
t361 = m(6) * t362 - t396 * mrSges(6,1) - t395 * mrSges(6,3) - t405 * t447 - t407 * t446;
t431 = -m(5) * t369 + t396 * mrSges(5,1) - t395 * mrSges(5,2) - t404 * t447 + t406 * t446 - t361;
t355 = m(4) * t371 + t416 * mrSges(4,1) - t415 * mrSges(4,2) + t431;
t347 = t424 * t352 + t427 * t355;
t345 = m(3) * t378 + qJDD(1) * mrSges(3,1) - t430 * mrSges(3,2) + t347;
t436 = t427 * t352 - t424 * t355;
t346 = m(3) * t379 - t430 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t436;
t339 = t422 * t345 + t421 * t346;
t337 = m(2) * t408 + qJDD(1) * mrSges(2,1) - t430 * mrSges(2,2) + t339;
t437 = -t421 * t345 + t422 * t346;
t338 = m(2) * t409 - t430 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t437;
t444 = t428 * t337 + t425 * t338;
t353 = t423 * t359 + t426 * t360;
t443 = (-t450 * t423 - t453 * t426) * t417 - t448 * qJD(4);
t442 = (t449 * t423 + t448 * t426) * t417 + t452 * qJD(4);
t441 = (t454 * t423 + t450 * t426) * t417 + t449 * qJD(4);
t439 = m(4) * t420 + t353;
t438 = -t425 * t337 + t428 * t338;
t433 = m(3) * t420 + t439;
t349 = mrSges(5,2) * t369 + mrSges(6,2) * t365 - mrSges(5,3) * t366 - mrSges(6,3) * t362 - qJ(5) * t361 + t443 * qJD(4) + t449 * qJDD(4) + t454 * t395 + t450 * t396 + t442 * t446;
t348 = -mrSges(5,1) * t369 - mrSges(6,1) * t362 + mrSges(6,2) * t364 + mrSges(5,3) * t367 - pkin(4) * t361 + t441 * qJD(4) + t448 * qJDD(4) + t450 * t395 + t453 * t396 - t442 * t447;
t341 = Ifges(4,6) * t416 + t415 * Ifges(4,5) - mrSges(4,1) * t420 + mrSges(4,3) * t372 - mrSges(5,1) * t366 + mrSges(5,2) * t367 + mrSges(6,1) * t365 - mrSges(6,3) * t364 - pkin(4) * t432 - qJ(5) * t434 - pkin(3) * t353 + (-qJ(5) * mrSges(6,2) - t448) * t396 + (pkin(4) * mrSges(6,2) - t449) * t395 - t452 * qJDD(4) + (t441 * t426 + (pkin(4) * t393 + t443) * t423) * t417;
t340 = mrSges(4,2) * t420 - mrSges(4,3) * t371 + Ifges(4,5) * t416 - t415 * Ifges(4,6) - pkin(7) * t353 - t423 * t348 + t426 * t349;
t333 = mrSges(3,2) * t420 - mrSges(3,3) * t378 + Ifges(3,5) * qJDD(1) - t430 * Ifges(3,6) - pkin(6) * t347 + t427 * t340 - t424 * t341;
t332 = -mrSges(3,1) * t420 + mrSges(3,3) * t379 + t430 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t439 + pkin(6) * t436 + t424 * t340 + t427 * t341;
t331 = -mrSges(2,2) * g(3) - mrSges(2,3) * t408 + Ifges(2,5) * qJDD(1) - t430 * Ifges(2,6) - qJ(2) * t339 - t421 * t332 + t422 * t333;
t330 = mrSges(2,1) * g(3) + mrSges(2,3) * t409 + t430 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t433 + qJ(2) * t437 + t422 * t332 + t421 * t333;
t1 = [-m(1) * g(1) + t438; -m(1) * g(2) + t444; (-m(1) - m(2)) * g(3) + t433; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t444 - t425 * t330 + t428 * t331; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t438 + t428 * t330 + t425 * t331; pkin(1) * t339 + mrSges(2,1) * t408 - mrSges(2,2) * t409 + pkin(2) * t347 + mrSges(3,1) * t378 - mrSges(3,2) * t379 + pkin(7) * t435 + t423 * t349 + t426 * t348 + pkin(3) * t431 + mrSges(4,1) * t371 - mrSges(4,2) * t372 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t416 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
