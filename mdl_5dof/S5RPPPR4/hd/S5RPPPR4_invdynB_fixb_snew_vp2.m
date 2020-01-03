% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:12
% EndTime: 2019-12-31 17:45:13
% DurationCPUTime: 1.32s
% Computational Cost: add. (11933->197), mult. (23143->240), div. (0->0), fcn. (12887->8), ass. (0->89)
t425 = sin(qJ(1));
t427 = cos(qJ(1));
t402 = t425 * g(1) - t427 * g(2);
t400 = qJDD(1) * pkin(1) + t402;
t403 = -t427 * g(1) - t425 * g(2);
t428 = qJD(1) ^ 2;
t401 = -t428 * pkin(1) + t403;
t421 = sin(pkin(7));
t423 = cos(pkin(7));
t389 = t423 * t400 - t421 * t401;
t433 = -t428 * qJ(3) + qJDD(3) - t389;
t456 = -pkin(2) - qJ(4);
t462 = -(2 * qJD(1) * qJD(4)) + t456 * qJDD(1) + t433;
t420 = sin(pkin(8));
t413 = t420 ^ 2;
t422 = cos(pkin(8));
t453 = t422 ^ 2 + t413;
t446 = t453 * mrSges(5,3);
t390 = t421 * t400 + t423 * t401;
t461 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t390;
t460 = pkin(4) * t428;
t459 = mrSges(3,1) - mrSges(4,2);
t458 = -Ifges(4,4) + Ifges(3,5);
t457 = Ifges(3,6) - Ifges(4,5);
t417 = -g(3) + qJDD(2);
t449 = qJDD(1) * t422;
t454 = t462 * t422;
t368 = -pkin(6) * t449 + (-t422 * t460 - t417) * t420 + t454;
t373 = t422 * t417 + t462 * t420;
t450 = qJDD(1) * t420;
t369 = -pkin(6) * t450 - t413 * t460 + t373;
t424 = sin(qJ(5));
t426 = cos(qJ(5));
t366 = t426 * t368 - t424 * t369;
t438 = -t420 * t426 - t422 * t424;
t393 = t438 * qJD(1);
t437 = -t420 * t424 + t422 * t426;
t394 = t437 * qJD(1);
t385 = -t393 * mrSges(6,1) + t394 * mrSges(6,2);
t388 = t393 * qJD(5) + qJDD(1) * t437;
t391 = -qJD(5) * mrSges(6,2) + t393 * mrSges(6,3);
t364 = m(6) * t366 + qJDD(5) * mrSges(6,1) - t388 * mrSges(6,3) + qJD(5) * t391 - t394 * t385;
t367 = t424 * t368 + t426 * t369;
t387 = -t394 * qJD(5) + qJDD(1) * t438;
t392 = qJD(5) * mrSges(6,1) - t394 * mrSges(6,3);
t365 = m(6) * t367 - qJDD(5) * mrSges(6,2) + t387 * mrSges(6,3) - qJD(5) * t392 + t393 * t385;
t355 = t426 * t364 + t424 * t365;
t372 = -t420 * t417 + t454;
t435 = -qJDD(1) * mrSges(5,3) - t428 * (mrSges(5,1) * t420 + mrSges(5,2) * t422);
t353 = m(5) * t372 + t422 * t435 + t355;
t442 = -t424 * t364 + t426 * t365;
t354 = m(5) * t373 + t420 * t435 + t442;
t351 = t422 * t353 + t420 * t354;
t380 = -qJDD(1) * pkin(2) + t433;
t431 = -m(4) * t380 + t428 * mrSges(4,3) - t351;
t349 = m(3) * t389 - t428 * mrSges(3,2) + t459 * qJDD(1) + t431;
t379 = t428 * pkin(2) - t461;
t436 = qJDD(4) + t461;
t377 = t456 * t428 + t436;
t371 = pkin(4) * t450 + (-t453 * pkin(6) + t456) * t428 + t436;
t432 = m(6) * t371 - t387 * mrSges(6,1) + t388 * mrSges(6,2) - t393 * t391 + t394 * t392;
t430 = -m(5) * t377 - mrSges(5,1) * t450 - mrSges(5,2) * t449 - t432;
t429 = -m(4) * t379 + t428 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t430;
t360 = -qJDD(1) * mrSges(3,2) + (-mrSges(3,1) - t446) * t428 + t429 + m(3) * t390;
t345 = t423 * t349 + t421 * t360;
t343 = m(2) * t402 + qJDD(1) * mrSges(2,1) - t428 * mrSges(2,2) + t345;
t444 = -t421 * t349 + t423 * t360;
t344 = m(2) * t403 - t428 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t444;
t455 = t427 * t343 + t425 * t344;
t439 = Ifges(5,5) * t422 - Ifges(5,6) * t420;
t452 = t428 * t439;
t445 = -t425 * t343 + t427 * t344;
t443 = -t420 * t353 + t422 * t354;
t350 = m(4) * t417 + t443;
t441 = Ifges(5,1) * t422 - Ifges(5,4) * t420;
t440 = Ifges(5,4) * t422 - Ifges(5,2) * t420;
t434 = m(3) * t417 + t350;
t383 = Ifges(6,1) * t394 + Ifges(6,4) * t393 + Ifges(6,5) * qJD(5);
t382 = Ifges(6,4) * t394 + Ifges(6,2) * t393 + Ifges(6,6) * qJD(5);
t381 = Ifges(6,5) * t394 + Ifges(6,6) * t393 + Ifges(6,3) * qJD(5);
t357 = mrSges(6,2) * t371 - mrSges(6,3) * t366 + Ifges(6,1) * t388 + Ifges(6,4) * t387 + Ifges(6,5) * qJDD(5) - qJD(5) * t382 + t393 * t381;
t356 = -mrSges(6,1) * t371 + mrSges(6,3) * t367 + Ifges(6,4) * t388 + Ifges(6,2) * t387 + Ifges(6,6) * qJDD(5) + qJD(5) * t383 - t394 * t381;
t347 = mrSges(5,2) * t377 - mrSges(5,3) * t372 - pkin(6) * t355 + qJDD(1) * t441 - t424 * t356 + t426 * t357 - t420 * t452;
t346 = -mrSges(5,1) * t377 + mrSges(5,3) * t373 - pkin(4) * t432 + pkin(6) * t442 + t440 * qJDD(1) + t426 * t356 + t424 * t357 - t422 * t452;
t339 = mrSges(4,1) * t380 + mrSges(5,1) * t372 + mrSges(6,1) * t366 - mrSges(5,2) * t373 - mrSges(6,2) * t367 - mrSges(3,3) * t389 + Ifges(6,5) * t388 + Ifges(6,6) * t387 + Ifges(6,3) * qJDD(5) + pkin(3) * t351 + pkin(4) * t355 - qJ(3) * t350 + t394 * t382 - t393 * t383 + (mrSges(3,2) - mrSges(4,3)) * t417 + (t439 + t458) * qJDD(1) + (t420 * t441 + t422 * t440 - t457) * t428;
t338 = mrSges(3,3) * t390 - mrSges(4,1) * t379 - t420 * t347 - t422 * t346 - pkin(3) * t430 - qJ(4) * t443 - pkin(2) * t350 - t459 * t417 + t457 * qJDD(1) + (-pkin(3) * t446 + t458) * t428;
t337 = -mrSges(2,2) * g(3) - mrSges(2,3) * t402 + Ifges(2,5) * qJDD(1) - t428 * Ifges(2,6) - qJ(2) * t345 - t421 * t338 + t423 * t339;
t336 = mrSges(2,1) * g(3) + mrSges(2,3) * t403 + t428 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t434 + qJ(2) * t444 + t423 * t338 + t421 * t339;
t1 = [-m(1) * g(1) + t445; -m(1) * g(2) + t455; (-m(1) - m(2)) * g(3) + t434; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t455 - t425 * t336 + t427 * t337; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t445 + t427 * t336 + t425 * t337; pkin(1) * t345 + mrSges(2,1) * t402 - mrSges(2,2) * t403 + pkin(2) * t431 + qJ(3) * (-t428 * t446 + t429) + t422 * t347 - t420 * t346 - qJ(4) * t351 + mrSges(3,1) * t389 - mrSges(3,2) * t390 + mrSges(4,2) * t380 - mrSges(4,3) * t379 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * mrSges(4,2) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
