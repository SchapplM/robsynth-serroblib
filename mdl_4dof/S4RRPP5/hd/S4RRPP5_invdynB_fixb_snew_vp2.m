% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPP5
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP5_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:15
% EndTime: 2019-12-31 17:00:16
% DurationCPUTime: 0.80s
% Computational Cost: add. (2830->196), mult. (5765->229), div. (0->0), fcn. (2301->4), ass. (0->75)
t450 = Ifges(3,1) + Ifges(4,2) + Ifges(5,3);
t438 = Ifges(3,4) + Ifges(4,6) - Ifges(5,6);
t437 = Ifges(3,5) - Ifges(4,4) + Ifges(5,5);
t449 = -Ifges(3,2) - Ifges(5,2) - Ifges(4,3);
t436 = Ifges(3,6) - Ifges(4,5) - Ifges(5,4);
t448 = Ifges(3,3) + Ifges(4,1) + Ifges(5,1);
t414 = sin(qJ(2));
t439 = qJD(1) * qJD(2);
t431 = t414 * t439;
t441 = qJD(1) * t414;
t446 = -2 * qJD(3);
t447 = pkin(2) * t431 + t441 * t446;
t445 = -2 * qJD(4);
t444 = -mrSges(5,1) - mrSges(3,3);
t415 = sin(qJ(1));
t417 = cos(qJ(1));
t404 = -g(1) * t417 - g(2) * t415;
t419 = qJD(1) ^ 2;
t377 = -pkin(1) * t419 + qJDD(1) * pkin(5) + t404;
t416 = cos(qJ(2));
t361 = -t414 * g(3) + t416 * t377;
t389 = (-mrSges(3,1) * t416 + mrSges(3,2) * t414) * qJD(1);
t392 = qJDD(1) * t416 - t431;
t396 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t441;
t387 = (-pkin(2) * t416 - qJ(3) * t414) * qJD(1);
t418 = qJD(2) ^ 2;
t440 = qJD(1) * t416;
t421 = -t418 * pkin(2) + qJDD(2) * qJ(3) + t387 * t440 + t361;
t358 = qJD(2) * t446 - t421;
t388 = (mrSges(4,2) * t416 - mrSges(4,3) * t414) * qJD(1);
t402 = mrSges(4,1) * t441 + qJD(2) * mrSges(4,2);
t398 = pkin(3) * t441 - qJD(2) * qJ(4);
t413 = t416 ^ 2;
t356 = -t413 * t419 * qJ(4) + t392 * pkin(3) + qJDD(4) + ((2 * qJD(3)) + t398) * qJD(2) + t421;
t390 = (-mrSges(5,2) * t414 - mrSges(5,3) * t416) * qJD(1);
t399 = mrSges(5,1) * t441 - qJD(2) * mrSges(5,3);
t428 = -m(5) * t356 - qJDD(2) * mrSges(5,2) - qJD(2) * t399 - t390 * t440;
t422 = -m(4) * t358 + qJDD(2) * mrSges(4,3) + qJD(2) * t402 + t388 * t440 - t428;
t348 = t389 * t440 + m(3) * t361 - qJDD(2) * mrSges(3,2) - qJD(2) * t396 + (mrSges(4,1) - t444) * t392 + t422;
t360 = -t416 * g(3) - t414 * t377;
t432 = t416 * t439;
t391 = qJDD(1) * t414 + t432;
t397 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t440;
t400 = -mrSges(4,1) * t440 - qJD(2) * mrSges(4,3);
t359 = -qJDD(2) * pkin(2) - t418 * qJ(3) + t387 * t441 + qJDD(3) - t360;
t355 = qJD(2) * t445 + (-t414 * t416 * t419 - qJDD(2)) * qJ(4) + (t391 - t432) * pkin(3) + t359;
t401 = mrSges(5,1) * t440 + qJD(2) * mrSges(5,2);
t427 = m(5) * t355 - qJDD(2) * mrSges(5,3) - qJD(2) * t401;
t424 = -m(4) * t359 - t391 * mrSges(4,1) - t427;
t442 = -t388 - t390;
t349 = m(3) * t360 + t444 * t391 + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t397 - t400) * qJD(2) + (-t389 + t442) * t441 + t424;
t429 = t416 * t348 - t349 * t414;
t341 = m(2) * t404 - mrSges(2,1) * t419 - qJDD(1) * mrSges(2,2) + t429;
t403 = t415 * g(1) - t417 * g(2);
t425 = -qJDD(1) * pkin(1) - t403;
t376 = -t419 * pkin(5) + t425;
t357 = -t392 * pkin(2) + (-t391 - t432) * qJ(3) + t376 + t447;
t353 = -t391 * qJ(3) + (-pkin(3) * t413 - pkin(5)) * t419 + (-pkin(2) - qJ(4)) * t392 + (-t398 * t414 + (-qJ(3) * qJD(2) + t445) * t416) * qJD(1) + t425 + t447;
t426 = m(5) * t353 - t391 * mrSges(5,2) - t392 * mrSges(5,3) - t399 * t441 - t401 * t440;
t423 = -m(4) * t357 - t392 * mrSges(4,2) + t402 * t441 - t426;
t420 = -m(3) * t376 + t397 * t440 + t392 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t391 + (-t396 * t414 - t400 * t416) * qJD(1) + t423;
t344 = m(2) * t403 + qJDD(1) * mrSges(2,1) - t419 * mrSges(2,2) + t420;
t443 = t415 * t341 + t417 * t344;
t342 = t414 * t348 + t416 * t349;
t435 = t448 * qJD(2) + (t437 * t414 + t436 * t416) * qJD(1);
t434 = -t436 * qJD(2) + (-t438 * t414 + t449 * t416) * qJD(1);
t433 = t437 * qJD(2) + (t450 * t414 + t438 * t416) * qJD(1);
t430 = t417 * t341 - t344 * t415;
t351 = t391 * mrSges(5,1) + t390 * t441 + t427;
t350 = -t391 * mrSges(4,3) + t400 * t440 - t423;
t338 = mrSges(4,1) * t359 + mrSges(5,1) * t355 + mrSges(3,2) * t376 - mrSges(5,2) * t353 - mrSges(3,3) * t360 - mrSges(4,3) * t357 + pkin(3) * t351 - qJ(3) * t350 + t434 * qJD(2) + t437 * qJDD(2) + t450 * t391 + t438 * t392 + t435 * t440;
t337 = -mrSges(3,1) * t376 + mrSges(3,3) * t361 - mrSges(4,1) * t358 + mrSges(4,2) * t357 + mrSges(5,1) * t356 - mrSges(5,3) * t353 - pkin(3) * t428 - qJ(4) * t426 - pkin(2) * t350 + (mrSges(5,1) * pkin(3) - t449) * t392 + t438 * t391 + t436 * qJDD(2) + t433 * qJD(2) - t435 * t441;
t336 = qJ(4) * t351 + mrSges(5,3) * t355 - mrSges(5,2) * t356 - pkin(1) * t342 + mrSges(4,3) * t358 - mrSges(4,2) * t359 - mrSges(3,1) * t360 + mrSges(3,2) * t361 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + mrSges(2,3) * t404 - qJ(3) * t422 - pkin(2) * (-qJD(2) * t400 + t424) + t419 * Ifges(2,5) + (mrSges(5,1) * pkin(2) - t437) * t391 + (mrSges(4,2) * pkin(2) - t448) * qJDD(2) + (-qJ(3) * (mrSges(4,1) + mrSges(5,1)) - t436) * t392 + (t433 * t416 + (-pkin(2) * t442 + t434) * t414) * qJD(1);
t335 = -mrSges(2,2) * g(3) - mrSges(2,3) * t403 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t419 - pkin(5) * t342 - t337 * t414 + t338 * t416;
t1 = [-m(1) * g(1) + t430; -m(1) * g(2) + t443; (-m(1) - m(2)) * g(3) + t342; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t443 + t417 * t335 - t415 * t336; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t430 + t415 * t335 + t417 * t336; -mrSges(1,1) * g(2) + mrSges(2,1) * t403 + mrSges(1,2) * g(1) - mrSges(2,2) * t404 + Ifges(2,3) * qJDD(1) + pkin(1) * t420 + pkin(5) * t429 + t416 * t337 + t414 * t338;];
tauB = t1;
