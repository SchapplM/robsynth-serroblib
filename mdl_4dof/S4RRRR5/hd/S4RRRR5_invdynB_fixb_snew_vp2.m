% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRR5
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR5_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR5_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR5_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:33
% EndTime: 2019-12-31 17:27:36
% DurationCPUTime: 1.82s
% Computational Cost: add. (19971->238), mult. (39749->303), div. (0->0), fcn. (25288->8), ass. (0->94)
t411 = sin(qJ(1));
t415 = cos(qJ(1));
t402 = -g(1) * t415 - g(2) * t411;
t417 = qJD(1) ^ 2;
t387 = -pkin(1) * t417 + qJDD(1) * pkin(5) + t402;
t410 = sin(qJ(2));
t414 = cos(qJ(2));
t378 = -g(3) * t410 + t414 * t387;
t395 = (-mrSges(3,1) * t414 + mrSges(3,2) * t410) * qJD(1);
t426 = qJD(1) * qJD(2);
t405 = t410 * t426;
t398 = qJDD(1) * t414 - t405;
t428 = qJD(1) * t410;
t399 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t428;
t401 = g(1) * t411 - t415 * g(2);
t386 = -qJDD(1) * pkin(1) - pkin(5) * t417 - t401;
t425 = t414 * t426;
t397 = qJDD(1) * t410 + t425;
t359 = (-t397 - t425) * pkin(6) + (-t398 + t405) * pkin(2) + t386;
t396 = (-pkin(2) * t414 - pkin(6) * t410) * qJD(1);
t416 = qJD(2) ^ 2;
t427 = qJD(1) * t414;
t362 = -pkin(2) * t416 + qJDD(2) * pkin(6) + t396 * t427 + t378;
t409 = sin(qJ(3));
t413 = cos(qJ(3));
t350 = t413 * t359 - t362 * t409;
t393 = qJD(2) * t413 - t409 * t428;
t371 = qJD(3) * t393 + qJDD(2) * t409 + t397 * t413;
t392 = qJDD(3) - t398;
t394 = qJD(2) * t409 + t413 * t428;
t404 = qJD(3) - t427;
t344 = (t393 * t404 - t371) * pkin(7) + (t393 * t394 + t392) * pkin(3) + t350;
t351 = t409 * t359 + t413 * t362;
t370 = -qJD(3) * t394 + qJDD(2) * t413 - t397 * t409;
t379 = pkin(3) * t404 - pkin(7) * t394;
t391 = t393 ^ 2;
t345 = -pkin(3) * t391 + pkin(7) * t370 - t379 * t404 + t351;
t408 = sin(qJ(4));
t412 = cos(qJ(4));
t342 = t344 * t412 - t345 * t408;
t372 = t393 * t412 - t394 * t408;
t349 = qJD(4) * t372 + t370 * t408 + t371 * t412;
t373 = t393 * t408 + t394 * t412;
t356 = -mrSges(5,1) * t372 + mrSges(5,2) * t373;
t403 = qJD(4) + t404;
t363 = -mrSges(5,2) * t403 + mrSges(5,3) * t372;
t388 = qJDD(4) + t392;
t340 = m(5) * t342 + mrSges(5,1) * t388 - mrSges(5,3) * t349 - t356 * t373 + t363 * t403;
t343 = t344 * t408 + t345 * t412;
t348 = -qJD(4) * t373 + t370 * t412 - t371 * t408;
t364 = mrSges(5,1) * t403 - mrSges(5,3) * t373;
t341 = m(5) * t343 - mrSges(5,2) * t388 + mrSges(5,3) * t348 + t356 * t372 - t364 * t403;
t332 = t412 * t340 + t408 * t341;
t374 = -mrSges(4,1) * t393 + mrSges(4,2) * t394;
t375 = -mrSges(4,2) * t404 + mrSges(4,3) * t393;
t330 = m(4) * t350 + mrSges(4,1) * t392 - mrSges(4,3) * t371 - t374 * t394 + t375 * t404 + t332;
t376 = mrSges(4,1) * t404 - mrSges(4,3) * t394;
t421 = -t340 * t408 + t412 * t341;
t331 = m(4) * t351 - mrSges(4,2) * t392 + mrSges(4,3) * t370 + t374 * t393 - t376 * t404 + t421;
t422 = -t409 * t330 + t413 * t331;
t327 = m(3) * t378 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t398 - qJD(2) * t399 + t395 * t427 + t422;
t377 = -t414 * g(3) - t410 * t387;
t400 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t427;
t361 = -qJDD(2) * pkin(2) - pkin(6) * t416 + t396 * t428 - t377;
t346 = -pkin(3) * t370 - pkin(7) * t391 + t379 * t394 + t361;
t420 = m(5) * t346 - t348 * mrSges(5,1) + mrSges(5,2) * t349 - t372 * t363 + t364 * t373;
t418 = -m(4) * t361 + t370 * mrSges(4,1) - mrSges(4,2) * t371 + t393 * t375 - t376 * t394 - t420;
t336 = m(3) * t377 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t397 + qJD(2) * t400 - t395 * t428 + t418;
t423 = t414 * t327 - t336 * t410;
t320 = m(2) * t402 - mrSges(2,1) * t417 - qJDD(1) * mrSges(2,2) + t423;
t328 = t330 * t413 + t331 * t409;
t419 = -m(3) * t386 + t398 * mrSges(3,1) - mrSges(3,2) * t397 - t399 * t428 + t400 * t427 - t328;
t324 = m(2) * t401 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t417 + t419;
t429 = t411 * t320 + t415 * t324;
t321 = t410 * t327 + t414 * t336;
t424 = t415 * t320 - t324 * t411;
t385 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t410 + Ifges(3,4) * t414) * qJD(1);
t384 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t410 + Ifges(3,2) * t414) * qJD(1);
t383 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t410 + Ifges(3,6) * t414) * qJD(1);
t367 = Ifges(4,1) * t394 + Ifges(4,4) * t393 + Ifges(4,5) * t404;
t366 = Ifges(4,4) * t394 + Ifges(4,2) * t393 + Ifges(4,6) * t404;
t365 = Ifges(4,5) * t394 + Ifges(4,6) * t393 + Ifges(4,3) * t404;
t354 = Ifges(5,1) * t373 + Ifges(5,4) * t372 + Ifges(5,5) * t403;
t353 = Ifges(5,4) * t373 + Ifges(5,2) * t372 + Ifges(5,6) * t403;
t352 = Ifges(5,5) * t373 + Ifges(5,6) * t372 + Ifges(5,3) * t403;
t334 = mrSges(5,2) * t346 - mrSges(5,3) * t342 + Ifges(5,1) * t349 + Ifges(5,4) * t348 + Ifges(5,5) * t388 + t352 * t372 - t353 * t403;
t333 = -mrSges(5,1) * t346 + mrSges(5,3) * t343 + Ifges(5,4) * t349 + Ifges(5,2) * t348 + Ifges(5,6) * t388 - t352 * t373 + t354 * t403;
t322 = mrSges(4,2) * t361 - mrSges(4,3) * t350 + Ifges(4,1) * t371 + Ifges(4,4) * t370 + Ifges(4,5) * t392 - pkin(7) * t332 - t333 * t408 + t334 * t412 + t365 * t393 - t366 * t404;
t317 = -mrSges(4,1) * t361 + mrSges(4,3) * t351 + Ifges(4,4) * t371 + Ifges(4,2) * t370 + Ifges(4,6) * t392 - pkin(3) * t420 + pkin(7) * t421 + t412 * t333 + t408 * t334 - t394 * t365 + t404 * t367;
t316 = Ifges(3,4) * t397 + Ifges(3,2) * t398 + Ifges(3,6) * qJDD(2) - t383 * t428 + qJD(2) * t385 - mrSges(3,1) * t386 + mrSges(3,3) * t378 - Ifges(4,5) * t371 - Ifges(4,6) * t370 - Ifges(4,3) * t392 - t394 * t366 + t393 * t367 - mrSges(4,1) * t350 + mrSges(4,2) * t351 - Ifges(5,5) * t349 - Ifges(5,6) * t348 - Ifges(5,3) * t388 - t373 * t353 + t372 * t354 - mrSges(5,1) * t342 + mrSges(5,2) * t343 - pkin(3) * t332 - pkin(2) * t328;
t315 = mrSges(3,2) * t386 - mrSges(3,3) * t377 + Ifges(3,1) * t397 + Ifges(3,4) * t398 + Ifges(3,5) * qJDD(2) - pkin(6) * t328 - qJD(2) * t384 - t317 * t409 + t322 * t413 + t383 * t427;
t314 = Ifges(2,6) * qJDD(1) + t417 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t402 - Ifges(3,5) * t397 - Ifges(3,6) * t398 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t377 + mrSges(3,2) * t378 - t409 * t322 - t413 * t317 - pkin(2) * t418 - pkin(6) * t422 - pkin(1) * t321 + (-t384 * t410 + t385 * t414) * qJD(1);
t313 = -mrSges(2,2) * g(3) - mrSges(2,3) * t401 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t417 - pkin(5) * t321 + t315 * t414 - t316 * t410;
t1 = [-m(1) * g(1) + t424; -m(1) * g(2) + t429; (-m(1) - m(2)) * g(3) + t321; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t429 + t415 * t313 - t411 * t314; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t424 + t411 * t313 + t415 * t314; -mrSges(1,1) * g(2) + mrSges(2,1) * t401 + mrSges(1,2) * g(1) - mrSges(2,2) * t402 + Ifges(2,3) * qJDD(1) + pkin(1) * t419 + pkin(5) * t423 + t410 * t315 + t414 * t316;];
tauB = t1;
