% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:17
% EndTime: 2019-12-31 17:18:19
% DurationCPUTime: 1.25s
% Computational Cost: add. (8218->216), mult. (16148->261), div. (0->0), fcn. (9258->6), ass. (0->85)
t436 = Ifges(4,1) + Ifges(5,1);
t431 = Ifges(4,4) + Ifges(5,4);
t430 = Ifges(4,5) + Ifges(5,5);
t435 = Ifges(4,2) + Ifges(5,2);
t434 = Ifges(4,6) + Ifges(5,6);
t433 = Ifges(4,3) + Ifges(5,3);
t432 = -mrSges(4,2) - mrSges(5,2);
t405 = sin(qJ(1));
t408 = cos(qJ(1));
t399 = -t408 * g(1) - t405 * g(2);
t410 = qJD(1) ^ 2;
t384 = -t410 * pkin(1) + qJDD(1) * pkin(5) + t399;
t404 = sin(qJ(2));
t407 = cos(qJ(2));
t376 = -t404 * g(3) + t407 * t384;
t392 = (-mrSges(3,1) * t407 + mrSges(3,2) * t404) * qJD(1);
t421 = qJD(1) * qJD(2);
t418 = t404 * t421;
t395 = t407 * qJDD(1) - t418;
t423 = qJD(1) * t404;
t396 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t423;
t398 = t405 * g(1) - t408 * g(2);
t383 = -qJDD(1) * pkin(1) - t410 * pkin(5) - t398;
t417 = t407 * t421;
t394 = t404 * qJDD(1) + t417;
t349 = (-t394 - t417) * pkin(6) + (-t395 + t418) * pkin(2) + t383;
t393 = (-pkin(2) * t407 - pkin(6) * t404) * qJD(1);
t409 = qJD(2) ^ 2;
t422 = t407 * qJD(1);
t352 = -t409 * pkin(2) + qJDD(2) * pkin(6) + t393 * t422 + t376;
t403 = sin(qJ(3));
t406 = cos(qJ(3));
t345 = t406 * t349 - t403 * t352;
t390 = t406 * qJD(2) - t403 * t423;
t366 = t390 * qJD(3) + t403 * qJDD(2) + t406 * t394;
t391 = t403 * qJD(2) + t406 * t423;
t368 = -t390 * mrSges(5,1) + t391 * mrSges(5,2);
t369 = -t390 * mrSges(4,1) + t391 * mrSges(4,2);
t400 = qJD(3) - t422;
t371 = -t400 * mrSges(4,2) + t390 * mrSges(4,3);
t389 = qJDD(3) - t395;
t341 = -0.2e1 * qJD(4) * t391 + (t390 * t400 - t366) * qJ(4) + (t390 * t391 + t389) * pkin(3) + t345;
t370 = -t400 * mrSges(5,2) + t390 * mrSges(5,3);
t420 = m(5) * t341 + t389 * mrSges(5,1) + t400 * t370;
t334 = m(4) * t345 + t389 * mrSges(4,1) + t400 * t371 + (-t368 - t369) * t391 + (-mrSges(4,3) - mrSges(5,3)) * t366 + t420;
t346 = t403 * t349 + t406 * t352;
t365 = -t391 * qJD(3) + t406 * qJDD(2) - t403 * t394;
t372 = t400 * pkin(3) - t391 * qJ(4);
t388 = t390 ^ 2;
t343 = -t388 * pkin(3) + t365 * qJ(4) + 0.2e1 * qJD(4) * t390 - t400 * t372 + t346;
t419 = m(5) * t343 + t365 * mrSges(5,3) + t390 * t368;
t373 = t400 * mrSges(5,1) - t391 * mrSges(5,3);
t424 = -t400 * mrSges(4,1) + t391 * mrSges(4,3) - t373;
t336 = m(4) * t346 + t365 * mrSges(4,3) + t390 * t369 + t432 * t389 + t424 * t400 + t419;
t414 = -t403 * t334 + t406 * t336;
t332 = m(3) * t376 - qJDD(2) * mrSges(3,2) + t395 * mrSges(3,3) - qJD(2) * t396 + t392 * t422 + t414;
t375 = -t407 * g(3) - t404 * t384;
t397 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t422;
t351 = -qJDD(2) * pkin(2) - t409 * pkin(6) + t393 * t423 - t375;
t344 = -t365 * pkin(3) - t388 * qJ(4) + t391 * t372 + qJDD(4) + t351;
t413 = m(5) * t344 - t365 * mrSges(5,1) - t390 * t370;
t411 = -m(4) * t351 + t365 * mrSges(4,1) + t432 * t366 + t390 * t371 + t424 * t391 - t413;
t338 = m(3) * t375 + qJDD(2) * mrSges(3,1) - t394 * mrSges(3,3) + qJD(2) * t397 - t392 * t423 + t411;
t415 = t407 * t332 - t404 * t338;
t324 = m(2) * t399 - t410 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t415;
t333 = t406 * t334 + t403 * t336;
t412 = -m(3) * t383 + t395 * mrSges(3,1) - t394 * mrSges(3,2) - t396 * t423 + t397 * t422 - t333;
t328 = m(2) * t398 + qJDD(1) * mrSges(2,1) - t410 * mrSges(2,2) + t412;
t428 = t405 * t324 + t408 * t328;
t325 = t404 * t332 + t407 * t338;
t427 = t434 * t390 + t430 * t391 + t433 * t400;
t426 = -t435 * t390 - t431 * t391 - t434 * t400;
t425 = t431 * t390 + t436 * t391 + t430 * t400;
t416 = t408 * t324 - t405 * t328;
t382 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t404 + Ifges(3,4) * t407) * qJD(1);
t381 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t404 + Ifges(3,2) * t407) * qJD(1);
t380 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t404 + Ifges(3,6) * t407) * qJD(1);
t339 = -t366 * mrSges(5,3) - t391 * t368 + t420;
t329 = mrSges(4,2) * t351 + mrSges(5,2) * t344 - mrSges(4,3) * t345 - mrSges(5,3) * t341 - qJ(4) * t339 + t431 * t365 + t436 * t366 + t430 * t389 + t427 * t390 + t426 * t400;
t326 = -mrSges(4,1) * t351 + mrSges(4,3) * t346 - mrSges(5,1) * t344 + mrSges(5,3) * t343 - pkin(3) * t413 + qJ(4) * t419 + (-qJ(4) * t373 + t425) * t400 + (-pkin(3) * t373 - t427) * t391 + (-qJ(4) * mrSges(5,2) + t434) * t389 + (-pkin(3) * mrSges(5,2) + t431) * t366 + t435 * t365;
t321 = -t380 * t423 - mrSges(3,1) * t383 - mrSges(4,1) * t345 - mrSges(5,1) * t341 + mrSges(4,2) * t346 + mrSges(5,2) * t343 + mrSges(3,3) * t376 + Ifges(3,4) * t394 + Ifges(3,2) * t395 + Ifges(3,6) * qJDD(2) - pkin(2) * t333 - pkin(3) * t339 + qJD(2) * t382 + t426 * t391 + t425 * t390 - t433 * t389 - t430 * t366 - t434 * t365;
t320 = mrSges(3,2) * t383 - mrSges(3,3) * t375 + Ifges(3,1) * t394 + Ifges(3,4) * t395 + Ifges(3,5) * qJDD(2) - pkin(6) * t333 - qJD(2) * t381 - t403 * t326 + t406 * t329 + t380 * t422;
t319 = Ifges(2,6) * qJDD(1) + t410 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t399 - Ifges(3,5) * t394 - Ifges(3,6) * t395 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t375 + mrSges(3,2) * t376 - t403 * t329 - t406 * t326 - pkin(2) * t411 - pkin(6) * t414 - pkin(1) * t325 + (-t404 * t381 + t407 * t382) * qJD(1);
t318 = -mrSges(2,2) * g(3) - mrSges(2,3) * t398 + Ifges(2,5) * qJDD(1) - t410 * Ifges(2,6) - pkin(5) * t325 + t407 * t320 - t404 * t321;
t1 = [-m(1) * g(1) + t416; -m(1) * g(2) + t428; (-m(1) - m(2)) * g(3) + t325; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t428 + t408 * t318 - t405 * t319; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t416 + t405 * t318 + t408 * t319; -mrSges(1,1) * g(2) + mrSges(2,1) * t398 + mrSges(1,2) * g(1) - mrSges(2,2) * t399 + Ifges(2,3) * qJDD(1) + pkin(1) * t412 + pkin(5) * t415 + t404 * t320 + t407 * t321;];
tauB = t1;
