% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP5_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP5_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP5_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:44
% EndTime: 2019-12-31 17:16:46
% DurationCPUTime: 1.19s
% Computational Cost: add. (8877->216), mult. (18384->266), div. (0->0), fcn. (10782->6), ass. (0->85)
t455 = Ifges(4,1) + Ifges(5,1);
t449 = Ifges(4,4) - Ifges(5,5);
t448 = Ifges(5,4) + Ifges(4,5);
t454 = Ifges(4,2) + Ifges(5,3);
t453 = -Ifges(5,2) - Ifges(4,3);
t447 = Ifges(4,6) - Ifges(5,6);
t452 = cos(qJ(3));
t429 = qJD(1) ^ 2;
t451 = pkin(2) * t429;
t450 = -mrSges(4,3) - mrSges(5,2);
t426 = sin(qJ(1));
t428 = cos(qJ(1));
t415 = -t428 * g(1) - t426 * g(2);
t403 = -t429 * pkin(1) + qJDD(1) * pkin(5) + t415;
t425 = sin(qJ(2));
t446 = t425 * t403;
t427 = cos(qJ(2));
t438 = qJD(1) * qJD(2);
t409 = t425 * qJDD(1) + t427 * t438;
t368 = qJDD(2) * pkin(2) - t409 * pkin(6) - t446 + (pkin(6) * t438 + t425 * t451 - g(3)) * t427;
t391 = -t425 * g(3) + t427 * t403;
t410 = t427 * qJDD(1) - t425 * t438;
t440 = qJD(1) * t425;
t413 = qJD(2) * pkin(2) - pkin(6) * t440;
t423 = t427 ^ 2;
t369 = t410 * pkin(6) - qJD(2) * t413 - t423 * t451 + t391;
t424 = sin(qJ(3));
t365 = t424 * t368 + t452 * t369;
t401 = (t424 * t427 + t452 * t425) * qJD(1);
t373 = t401 * qJD(3) + t424 * t409 - t452 * t410;
t422 = qJD(2) + qJD(3);
t393 = t422 * mrSges(4,1) - t401 * mrSges(4,3);
t439 = qJD(1) * t427;
t400 = t424 * t440 - t452 * t439;
t421 = qJDD(2) + qJDD(3);
t385 = t400 * pkin(3) - t401 * qJ(4);
t420 = t422 ^ 2;
t362 = -t420 * pkin(3) + t421 * qJ(4) + 0.2e1 * qJD(4) * t422 - t400 * t385 + t365;
t394 = -t422 * mrSges(5,1) + t401 * mrSges(5,2);
t437 = m(5) * t362 + t421 * mrSges(5,3) + t422 * t394;
t386 = t400 * mrSges(5,1) - t401 * mrSges(5,3);
t441 = -t400 * mrSges(4,1) - t401 * mrSges(4,2) - t386;
t355 = m(4) * t365 - t421 * mrSges(4,2) + t450 * t373 - t422 * t393 + t441 * t400 + t437;
t364 = t452 * t368 - t424 * t369;
t374 = -t400 * qJD(3) + t452 * t409 + t424 * t410;
t392 = -t422 * mrSges(4,2) - t400 * mrSges(4,3);
t363 = -t421 * pkin(3) - t420 * qJ(4) + t401 * t385 + qJDD(4) - t364;
t395 = -t400 * mrSges(5,2) + t422 * mrSges(5,3);
t433 = -m(5) * t363 + t421 * mrSges(5,1) + t422 * t395;
t357 = m(4) * t364 + t421 * mrSges(4,1) + t450 * t374 + t422 * t392 + t441 * t401 + t433;
t350 = t424 * t355 + t452 * t357;
t390 = -t427 * g(3) - t446;
t408 = (-mrSges(3,1) * t427 + mrSges(3,2) * t425) * qJD(1);
t412 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t439;
t348 = m(3) * t390 + qJDD(2) * mrSges(3,1) - t409 * mrSges(3,3) + qJD(2) * t412 - t408 * t440 + t350;
t411 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t440;
t434 = t452 * t355 - t424 * t357;
t349 = m(3) * t391 - qJDD(2) * mrSges(3,2) + t410 * mrSges(3,3) - qJD(2) * t411 + t408 * t439 + t434;
t435 = -t425 * t348 + t427 * t349;
t341 = m(2) * t415 - t429 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t435;
t414 = t426 * g(1) - t428 * g(2);
t432 = -qJDD(1) * pkin(1) - t414;
t402 = -t429 * pkin(5) + t432;
t375 = -t410 * pkin(2) + t413 * t440 + (-pkin(6) * t423 - pkin(5)) * t429 + t432;
t360 = -0.2e1 * qJD(4) * t401 + (t400 * t422 - t374) * qJ(4) + (t401 * t422 + t373) * pkin(3) + t375;
t358 = m(5) * t360 + t373 * mrSges(5,1) - t374 * mrSges(5,3) - t401 * t394 + t400 * t395;
t431 = m(4) * t375 + t373 * mrSges(4,1) + t374 * mrSges(4,2) + t400 * t392 + t401 * t393 + t358;
t430 = -m(3) * t402 + t410 * mrSges(3,1) - t409 * mrSges(3,2) - t411 * t440 + t412 * t439 - t431;
t352 = m(2) * t414 + qJDD(1) * mrSges(2,1) - t429 * mrSges(2,2) + t430;
t445 = t426 * t341 + t428 * t352;
t342 = t427 * t348 + t425 * t349;
t444 = t454 * t400 - t449 * t401 - t447 * t422;
t443 = t447 * t400 - t448 * t401 + t453 * t422;
t442 = -t449 * t400 + t455 * t401 + t448 * t422;
t436 = t428 * t341 - t426 * t352;
t399 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t425 + Ifges(3,4) * t427) * qJD(1);
t398 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t425 + Ifges(3,2) * t427) * qJD(1);
t397 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t425 + Ifges(3,6) * t427) * qJD(1);
t344 = mrSges(4,2) * t375 + mrSges(5,2) * t363 - mrSges(4,3) * t364 - mrSges(5,3) * t360 - qJ(4) * t358 - t449 * t373 + t455 * t374 + t443 * t400 + t448 * t421 + t444 * t422;
t343 = -mrSges(4,1) * t375 - mrSges(5,1) * t360 + mrSges(5,2) * t362 + mrSges(4,3) * t365 - pkin(3) * t358 - t454 * t373 + t449 * t374 + t443 * t401 + t447 * t421 + t442 * t422;
t338 = mrSges(3,2) * t402 - mrSges(3,3) * t390 + Ifges(3,1) * t409 + Ifges(3,4) * t410 + Ifges(3,5) * qJDD(2) - pkin(6) * t350 - qJD(2) * t398 - t424 * t343 + t452 * t344 + t397 * t439;
t337 = -mrSges(3,1) * t402 + mrSges(3,3) * t391 + Ifges(3,4) * t409 + Ifges(3,2) * t410 + Ifges(3,6) * qJDD(2) - pkin(2) * t431 + pkin(6) * t434 + qJD(2) * t399 + t452 * t343 + t424 * t344 - t397 * t440;
t336 = Ifges(2,6) * qJDD(1) - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) - qJ(4) * t437 + t429 * Ifges(2,5) - pkin(3) * t433 - Ifges(3,5) * t409 - Ifges(3,6) * t410 + mrSges(2,3) * t415 - mrSges(3,1) * t390 + mrSges(3,2) * t391 - mrSges(4,1) * t364 + mrSges(4,2) * t365 - mrSges(5,3) * t362 + mrSges(5,1) * t363 - pkin(2) * t350 - pkin(1) * t342 + t453 * t421 + (pkin(3) * t386 + t444) * t401 + (qJ(4) * t386 - t442) * t400 + (pkin(3) * mrSges(5,2) - t448) * t374 + (qJ(4) * mrSges(5,2) + t447) * t373 + (-t425 * t398 + t427 * t399) * qJD(1);
t335 = -mrSges(2,2) * g(3) - mrSges(2,3) * t414 + Ifges(2,5) * qJDD(1) - t429 * Ifges(2,6) - pkin(5) * t342 - t425 * t337 + t427 * t338;
t1 = [-m(1) * g(1) + t436; -m(1) * g(2) + t445; (-m(1) - m(2)) * g(3) + t342; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t445 + t428 * t335 - t426 * t336; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t436 + t426 * t335 + t428 * t336; -mrSges(1,1) * g(2) + mrSges(2,1) * t414 + mrSges(1,2) * g(1) - mrSges(2,2) * t415 + Ifges(2,3) * qJDD(1) + pkin(1) * t430 + pkin(5) * t435 + t427 * t337 + t425 * t338;];
tauB = t1;
