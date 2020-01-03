% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRP5
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP5_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:51
% EndTime: 2019-12-31 16:44:52
% DurationCPUTime: 1.04s
% Computational Cost: add. (6772->194), mult. (16006->237), div. (0->0), fcn. (9740->6), ass. (0->87)
t455 = Ifges(4,1) + Ifges(5,1);
t449 = Ifges(4,4) - Ifges(5,5);
t448 = Ifges(4,5) + Ifges(5,4);
t454 = Ifges(4,2) + Ifges(5,3);
t447 = Ifges(4,6) - Ifges(5,6);
t453 = -Ifges(4,3) - Ifges(5,2);
t418 = qJD(1) ^ 2;
t452 = cos(qJ(3));
t413 = cos(pkin(6));
t451 = pkin(2) * t413;
t450 = -mrSges(4,3) - mrSges(5,2);
t412 = sin(pkin(6));
t446 = mrSges(3,2) * t412;
t409 = t413 ^ 2;
t445 = t409 * t418;
t415 = sin(qJ(1));
t416 = cos(qJ(1));
t399 = -g(1) * t416 - g(2) * t415;
t395 = -pkin(1) * t418 + qJDD(1) * qJ(2) + t399;
t436 = qJD(1) * qJD(2);
t431 = -g(3) * t413 - 0.2e1 * t412 * t436;
t363 = (-pkin(5) * qJDD(1) + t418 * t451 - t395) * t412 + t431;
t384 = -g(3) * t412 + (t395 + 0.2e1 * t436) * t413;
t434 = qJDD(1) * t413;
t364 = -pkin(2) * t445 + pkin(5) * t434 + t384;
t414 = sin(qJ(3));
t360 = t414 * t363 + t452 * t364;
t432 = t413 * t452;
t435 = qJDD(1) * t412;
t421 = t452 * t412 + t413 * t414;
t394 = t421 * qJD(1);
t438 = qJD(3) * t394;
t381 = -qJDD(1) * t432 + t414 * t435 + t438;
t388 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t394;
t439 = qJD(1) * t412;
t393 = -qJD(1) * t432 + t414 * t439;
t374 = pkin(3) * t393 - qJ(4) * t394;
t417 = qJD(3) ^ 2;
t357 = -pkin(3) * t417 + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) - t374 * t393 + t360;
t389 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t394;
t433 = m(5) * t357 + qJDD(3) * mrSges(5,3) + qJD(3) * t389;
t375 = mrSges(5,1) * t393 - mrSges(5,3) * t394;
t440 = -mrSges(4,1) * t393 - mrSges(4,2) * t394 - t375;
t351 = m(4) * t360 - qJDD(3) * mrSges(4,2) - qJD(3) * t388 + t450 * t381 + t440 * t393 + t433;
t359 = t452 * t363 - t414 * t364;
t437 = t393 * qJD(3);
t382 = t421 * qJDD(1) - t437;
t387 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t393;
t358 = -qJDD(3) * pkin(3) - t417 * qJ(4) + t394 * t374 + qJDD(4) - t359;
t390 = -mrSges(5,2) * t393 + qJD(3) * mrSges(5,3);
t427 = -m(5) * t358 + qJDD(3) * mrSges(5,1) + qJD(3) * t390;
t352 = m(4) * t359 + qJDD(3) * mrSges(4,1) + qJD(3) * t387 + t450 * t382 + t440 * t394 + t427;
t345 = t414 * t351 + t452 * t352;
t383 = -t395 * t412 + t431;
t422 = mrSges(3,3) * qJDD(1) + t418 * (-mrSges(3,1) * t413 + t446);
t343 = m(3) * t383 - t422 * t412 + t345;
t428 = t452 * t351 - t352 * t414;
t344 = m(3) * t384 + t422 * t413 + t428;
t429 = -t343 * t412 + t413 * t344;
t336 = m(2) * t399 - mrSges(2,1) * t418 - qJDD(1) * mrSges(2,2) + t429;
t398 = g(1) * t415 - t416 * g(2);
t426 = qJDD(2) - t398;
t392 = -qJDD(1) * pkin(1) - qJ(2) * t418 + t426;
t408 = t412 ^ 2;
t380 = (-pkin(1) - t451) * qJDD(1) + (-qJ(2) + (-t408 - t409) * pkin(5)) * t418 + t426;
t355 = -0.2e1 * qJD(4) * t394 + (-t382 + t437) * qJ(4) + (t381 + t438) * pkin(3) + t380;
t353 = m(5) * t355 + t381 * mrSges(5,1) - t382 * mrSges(5,3) - t394 * t389 + t393 * t390;
t420 = m(4) * t380 + t381 * mrSges(4,1) + mrSges(4,2) * t382 + t393 * t387 + t388 * t394 + t353;
t419 = -m(3) * t392 + mrSges(3,1) * t434 - t420 + (t408 * t418 + t445) * mrSges(3,3);
t347 = t419 + m(2) * t398 - mrSges(2,2) * t418 + (mrSges(2,1) - t446) * qJDD(1);
t444 = t415 * t336 + t416 * t347;
t337 = t413 * t343 + t412 * t344;
t443 = -t447 * qJD(3) + t454 * t393 - t449 * t394;
t442 = t453 * qJD(3) + t447 * t393 - t448 * t394;
t441 = t448 * qJD(3) - t449 * t393 + t455 * t394;
t430 = t416 * t336 - t347 * t415;
t425 = Ifges(3,1) * t412 + Ifges(3,4) * t413;
t424 = Ifges(3,4) * t412 + Ifges(3,2) * t413;
t423 = Ifges(3,5) * t412 + Ifges(3,6) * t413;
t397 = t423 * qJD(1);
t339 = mrSges(4,2) * t380 + mrSges(5,2) * t358 - mrSges(4,3) * t359 - mrSges(5,3) * t355 - qJ(4) * t353 + t443 * qJD(3) + t448 * qJDD(3) - t449 * t381 + t455 * t382 + t442 * t393;
t338 = -mrSges(4,1) * t380 - mrSges(5,1) * t355 + mrSges(5,2) * t357 + mrSges(4,3) * t360 - pkin(3) * t353 + t441 * qJD(3) + t447 * qJDD(3) - t454 * t381 + t449 * t382 + t442 * t394;
t333 = t413 * qJD(1) * t397 + mrSges(3,2) * t392 - mrSges(3,3) * t383 - pkin(5) * t345 + t425 * qJDD(1) - t414 * t338 + t452 * t339;
t332 = -mrSges(3,1) * t392 + mrSges(3,3) * t384 - pkin(2) * t420 + pkin(5) * t428 + t424 * qJDD(1) + t452 * t338 + t414 * t339 - t397 * t439;
t331 = mrSges(2,1) * g(3) - mrSges(3,1) * t383 + mrSges(3,2) * t384 + mrSges(2,3) * t399 - qJ(4) * t433 - pkin(3) * t427 - pkin(1) * t337 - pkin(2) * t345 - mrSges(5,3) * t357 + mrSges(5,1) * t358 - mrSges(4,1) * t359 + mrSges(4,2) * t360 + (pkin(3) * t375 + t443) * t394 + (qJ(4) * t375 - t441) * t393 + (mrSges(5,2) * pkin(3) - t448) * t382 + (mrSges(5,2) * qJ(4) + t447) * t381 + t453 * qJDD(3) + (Ifges(2,6) - t423) * qJDD(1) + (-t412 * t424 + t413 * t425 + Ifges(2,5)) * t418;
t330 = -mrSges(2,2) * g(3) - mrSges(2,3) * t398 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t418 - qJ(2) * t337 - t332 * t412 + t333 * t413;
t1 = [-m(1) * g(1) + t430; -m(1) * g(2) + t444; (-m(1) - m(2)) * g(3) + t337; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t444 + t416 * t330 - t415 * t331; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t430 + t415 * t330 + t416 * t331; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t398 - mrSges(2,2) * t399 + t412 * t333 + t413 * t332 + pkin(1) * (-mrSges(3,2) * t435 + t419) + qJ(2) * t429;];
tauB = t1;
