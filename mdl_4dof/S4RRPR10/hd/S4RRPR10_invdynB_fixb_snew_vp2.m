% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPR10
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR10_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:16
% EndTime: 2019-12-31 17:11:18
% DurationCPUTime: 1.05s
% Computational Cost: add. (5954->222), mult. (12132->267), div. (0->0), fcn. (5964->6), ass. (0->89)
t450 = Ifges(3,1) + Ifges(4,2);
t445 = Ifges(3,4) + Ifges(4,6);
t444 = Ifges(3,5) - Ifges(4,4);
t449 = Ifges(3,2) + Ifges(4,3);
t443 = Ifges(3,6) - Ifges(4,5);
t448 = Ifges(3,3) + Ifges(4,1);
t447 = -2 * qJD(3);
t422 = qJD(1) ^ 2;
t446 = t422 * pkin(5);
t417 = sin(qJ(1));
t420 = cos(qJ(1));
t406 = -t420 * g(1) - t417 * g(2);
t385 = -t422 * pkin(1) + qJDD(1) * pkin(5) + t406;
t416 = sin(qJ(2));
t419 = cos(qJ(2));
t372 = -t419 * g(3) - t416 * t385;
t395 = (mrSges(4,2) * t419 - mrSges(4,3) * t416) * qJD(1);
t396 = (-mrSges(3,1) * t419 + mrSges(3,2) * t416) * qJD(1);
t435 = qJD(1) * qJD(2);
t433 = t419 * t435;
t397 = t416 * qJDD(1) + t433;
t437 = qJD(1) * t419;
t401 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t437;
t402 = -mrSges(4,1) * t437 - qJD(2) * mrSges(4,3);
t434 = t416 * t435;
t398 = t419 * qJDD(1) - t434;
t436 = t416 * qJD(1);
t404 = pkin(3) * t436 - qJD(2) * pkin(6);
t414 = t419 ^ 2;
t405 = t417 * g(1) - t420 * g(2);
t430 = -qJDD(1) * pkin(1) - t405;
t424 = pkin(2) * t434 + t436 * t447 + (-t397 - t433) * qJ(3) + t430;
t355 = -t404 * t436 + (-pkin(3) * t414 - pkin(5)) * t422 + (-pkin(2) - pkin(6)) * t398 + t424;
t394 = (-pkin(2) * t419 - qJ(3) * t416) * qJD(1);
t421 = qJD(2) ^ 2;
t361 = -qJDD(2) * pkin(2) - t421 * qJ(3) + t394 * t436 + qJDD(3) - t372;
t358 = (-t416 * t419 * t422 - qJDD(2)) * pkin(6) + (t397 - t433) * pkin(3) + t361;
t415 = sin(qJ(4));
t418 = cos(qJ(4));
t353 = -t415 * t355 + t418 * t358;
t392 = -t415 * qJD(2) - t418 * t437;
t368 = t392 * qJD(4) + t418 * qJDD(2) - t415 * t398;
t393 = t418 * qJD(2) - t415 * t437;
t369 = -t392 * mrSges(5,1) + t393 * mrSges(5,2);
t408 = qJD(4) + t436;
t370 = -t408 * mrSges(5,2) + t392 * mrSges(5,3);
t391 = qJDD(4) + t397;
t351 = m(5) * t353 + t391 * mrSges(5,1) - t368 * mrSges(5,3) - t393 * t369 + t408 * t370;
t354 = t418 * t355 + t415 * t358;
t367 = -t393 * qJD(4) - t415 * qJDD(2) - t418 * t398;
t371 = t408 * mrSges(5,1) - t393 * mrSges(5,3);
t352 = m(5) * t354 - t391 * mrSges(5,2) + t367 * mrSges(5,3) + t392 * t369 - t408 * t371;
t343 = t418 * t351 + t415 * t352;
t427 = -m(4) * t361 - t397 * mrSges(4,1) - t343;
t341 = m(3) * t372 - t397 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t401 - t402) * qJD(2) + (-t395 - t396) * t436 + t427;
t373 = -t416 * g(3) + t419 * t385;
t400 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t436;
t426 = -t421 * pkin(2) + qJDD(2) * qJ(3) + t394 * t437 + t373;
t360 = qJD(2) * t447 - t426;
t403 = mrSges(4,1) * t436 + qJD(2) * mrSges(4,2);
t357 = -t414 * t422 * pkin(6) + t398 * pkin(3) + ((2 * qJD(3)) + t404) * qJD(2) + t426;
t428 = -m(5) * t357 + t367 * mrSges(5,1) - t368 * mrSges(5,2) + t392 * t370 - t393 * t371;
t425 = -m(4) * t360 + qJDD(2) * mrSges(4,3) + qJD(2) * t403 + t395 * t437 - t428;
t348 = t396 * t437 + m(3) * t373 - qJDD(2) * mrSges(3,2) - qJD(2) * t400 + (mrSges(3,3) + mrSges(4,1)) * t398 + t425;
t431 = -t416 * t341 + t419 * t348;
t336 = m(2) * t406 - t422 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t431;
t384 = t430 - t446;
t359 = -t398 * pkin(2) + t424 - t446;
t441 = -t415 * t351 + t418 * t352;
t429 = -m(4) * t359 - t398 * mrSges(4,2) + t403 * t436 - t441;
t423 = -m(3) * t384 + t401 * t437 + t398 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t397 + (-t400 * t416 - t402 * t419) * qJD(1) + t429;
t339 = m(2) * t405 + qJDD(1) * mrSges(2,1) - t422 * mrSges(2,2) + t423;
t442 = t417 * t336 + t420 * t339;
t337 = t419 * t341 + t416 * t348;
t440 = t448 * qJD(2) + (t444 * t416 + t443 * t419) * qJD(1);
t439 = -t443 * qJD(2) + (-t445 * t416 - t449 * t419) * qJD(1);
t438 = t444 * qJD(2) + (t450 * t416 + t445 * t419) * qJD(1);
t432 = t420 * t336 - t417 * t339;
t364 = Ifges(5,1) * t393 + Ifges(5,4) * t392 + Ifges(5,5) * t408;
t363 = Ifges(5,4) * t393 + Ifges(5,2) * t392 + Ifges(5,6) * t408;
t362 = Ifges(5,5) * t393 + Ifges(5,6) * t392 + Ifges(5,3) * t408;
t345 = mrSges(5,2) * t357 - mrSges(5,3) * t353 + Ifges(5,1) * t368 + Ifges(5,4) * t367 + Ifges(5,5) * t391 + t392 * t362 - t408 * t363;
t344 = -mrSges(5,1) * t357 + mrSges(5,3) * t354 + Ifges(5,4) * t368 + Ifges(5,2) * t367 + Ifges(5,6) * t391 - t393 * t362 + t408 * t364;
t342 = -t397 * mrSges(4,3) + t402 * t437 - t429;
t333 = mrSges(4,1) * t361 + mrSges(5,1) * t353 + mrSges(3,2) * t384 - mrSges(5,2) * t354 - mrSges(3,3) * t372 - mrSges(4,3) * t359 + Ifges(5,5) * t368 + Ifges(5,6) * t367 + Ifges(5,3) * t391 + pkin(3) * t343 - qJ(3) * t342 + t393 * t363 - t392 * t364 + t445 * t398 + t450 * t397 + t444 * qJDD(2) + t439 * qJD(2) + t440 * t437;
t332 = -mrSges(3,1) * t384 - mrSges(4,1) * t360 + mrSges(4,2) * t359 + mrSges(3,3) * t373 - pkin(2) * t342 - pkin(3) * t428 - pkin(6) * t441 + t438 * qJD(2) + t443 * qJDD(2) - t418 * t344 - t415 * t345 + t445 * t397 + t449 * t398 - t440 * t436;
t331 = -pkin(1) * t337 + mrSges(2,3) * t406 + t415 * t344 + pkin(6) * t343 - mrSges(3,1) * t372 + mrSges(3,2) * t373 - pkin(2) * (-qJD(2) * t402 + t427) - qJ(3) * t425 - mrSges(4,2) * t361 + mrSges(4,3) * t360 - t418 * t345 + mrSges(2,1) * g(3) + t422 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-qJ(3) * mrSges(4,1) - t443) * t398 - t444 * t397 + (pkin(2) * mrSges(4,2) - t448) * qJDD(2) + (t438 * t419 + (pkin(2) * t395 + t439) * t416) * qJD(1);
t330 = -mrSges(2,2) * g(3) - mrSges(2,3) * t405 + Ifges(2,5) * qJDD(1) - t422 * Ifges(2,6) - pkin(5) * t337 - t416 * t332 + t419 * t333;
t1 = [-m(1) * g(1) + t432; -m(1) * g(2) + t442; (-m(1) - m(2)) * g(3) + t337; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t442 + t420 * t330 - t417 * t331; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t432 + t417 * t330 + t420 * t331; -mrSges(1,1) * g(2) + mrSges(2,1) * t405 + mrSges(1,2) * g(1) - mrSges(2,2) * t406 + Ifges(2,3) * qJDD(1) + pkin(1) * t423 + pkin(5) * t431 + t419 * t332 + t416 * t333;];
tauB = t1;
