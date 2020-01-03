% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:14
% EndTime: 2019-12-31 18:47:15
% DurationCPUTime: 0.98s
% Computational Cost: add. (8840->202), mult. (11290->235), div. (0->0), fcn. (4042->6), ass. (0->82)
t426 = Ifges(5,1) + Ifges(6,1);
t417 = Ifges(5,4) - Ifges(6,5);
t416 = Ifges(5,5) + Ifges(6,4);
t425 = Ifges(5,2) + Ifges(6,3);
t414 = Ifges(5,6) - Ifges(6,6);
t424 = Ifges(5,3) + Ifges(6,2);
t423 = -m(3) - m(4);
t422 = -pkin(1) - pkin(2);
t392 = cos(qJ(4));
t421 = t392 * g(3);
t420 = -mrSges(2,1) - mrSges(3,1);
t419 = mrSges(5,3) + mrSges(6,2);
t418 = Ifges(3,4) + Ifges(2,5);
t415 = Ifges(2,6) - Ifges(3,6);
t382 = -qJD(1) + qJD(3);
t389 = sin(qJ(4));
t413 = t382 * t389;
t412 = t382 * t392;
t391 = sin(qJ(1));
t394 = cos(qJ(1));
t375 = -t394 * g(1) - t391 * g(2);
t396 = qJD(1) ^ 2;
t400 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t375;
t356 = -t396 * pkin(1) + t400;
t347 = t422 * t396 + t400;
t374 = t391 * g(1) - t394 * g(2);
t399 = -t396 * qJ(2) + qJDD(2) - t374;
t349 = t422 * qJDD(1) + t399;
t390 = sin(qJ(3));
t393 = cos(qJ(3));
t342 = t393 * t347 + t390 * t349;
t380 = t382 ^ 2;
t381 = -qJDD(1) + qJDD(3);
t340 = -t380 * pkin(3) + t381 * pkin(7) + t342;
t337 = t389 * g(3) + t392 * t340;
t365 = (-mrSges(5,1) * t392 + mrSges(5,2) * t389) * t382;
t407 = qJD(4) * t382;
t367 = t392 * t381 - t389 * t407;
t370 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t413;
t363 = (-pkin(4) * t392 - qJ(5) * t389) * t382;
t395 = qJD(4) ^ 2;
t334 = -t395 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t363 * t412 + t337;
t364 = (-mrSges(6,1) * t392 - mrSges(6,3) * t389) * t382;
t371 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t413;
t403 = m(6) * t334 + qJDD(4) * mrSges(6,3) + qJD(4) * t371 + t364 * t412;
t329 = m(5) * t337 - qJDD(4) * mrSges(5,2) - qJD(4) * t370 + t365 * t412 + t419 * t367 + t403;
t336 = -t389 * t340 + t421;
t366 = t389 * t381 + t392 * t407;
t372 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t412;
t335 = -qJDD(4) * pkin(4) - t421 - t395 * qJ(5) + qJDD(5) + (t363 * t382 + t340) * t389;
t373 = mrSges(6,2) * t412 + qJD(4) * mrSges(6,3);
t402 = -m(6) * t335 + qJDD(4) * mrSges(6,1) + qJD(4) * t373;
t330 = m(5) * t336 + qJDD(4) * mrSges(5,1) + qJD(4) * t372 + (-t364 - t365) * t413 - t419 * t366 + t402;
t404 = t392 * t329 - t389 * t330;
t324 = m(4) * t342 - t380 * mrSges(4,1) - t381 * mrSges(4,2) + t404;
t341 = -t390 * t347 + t393 * t349;
t339 = -t381 * pkin(3) - t380 * pkin(7) - t341;
t332 = -t367 * pkin(4) - t366 * qJ(5) + (-0.2e1 * qJD(5) * t389 + (pkin(4) * t389 - qJ(5) * t392) * qJD(4)) * t382 + t339;
t331 = m(6) * t332 - t367 * mrSges(6,1) - t366 * mrSges(6,3) - t371 * t413 - t373 * t412;
t397 = -m(5) * t339 + t367 * mrSges(5,1) - t366 * mrSges(5,2) - t370 * t413 + t372 * t412 - t331;
t327 = m(4) * t341 + t381 * mrSges(4,1) - t380 * mrSges(4,2) + t397;
t405 = t393 * t324 - t390 * t327;
t401 = m(3) * t356 + qJDD(1) * mrSges(3,3) + t405;
t318 = m(2) * t375 - qJDD(1) * mrSges(2,2) + t420 * t396 + t401;
t320 = t390 * t324 + t393 * t327;
t362 = -qJDD(1) * pkin(1) + t399;
t398 = -m(3) * t362 + qJDD(1) * mrSges(3,1) + t396 * mrSges(3,3) - t320;
t319 = m(2) * t374 + qJDD(1) * mrSges(2,1) - t396 * mrSges(2,2) + t398;
t411 = t391 * t318 + t394 * t319;
t410 = (-t417 * t389 - t425 * t392) * t382 - t414 * qJD(4);
t409 = (t416 * t389 + t414 * t392) * t382 + t424 * qJD(4);
t408 = (t426 * t389 + t417 * t392) * t382 + t416 * qJD(4);
t406 = t394 * t318 - t391 * t319;
t326 = t389 * t329 + t392 * t330;
t325 = t423 * g(3) - t326;
t322 = mrSges(5,2) * t339 + mrSges(6,2) * t335 - mrSges(5,3) * t336 - mrSges(6,3) * t332 - qJ(5) * t331 + t410 * qJD(4) + t416 * qJDD(4) + t426 * t366 + t417 * t367 + t409 * t412;
t321 = -mrSges(5,1) * t339 - mrSges(6,1) * t332 + mrSges(6,2) * t334 + mrSges(5,3) * t337 - pkin(4) * t331 + t408 * qJD(4) + t414 * qJDD(4) + t417 * t366 + t425 * t367 - t409 * t413;
t314 = Ifges(4,6) * t381 + t380 * Ifges(4,5) - mrSges(4,1) * g(3) + mrSges(4,3) * t342 - mrSges(5,1) * t336 + mrSges(5,2) * t337 + mrSges(6,1) * t335 - mrSges(6,3) * t334 - pkin(4) * t402 - qJ(5) * t403 - pkin(3) * t326 + (-qJ(5) * mrSges(6,2) - t414) * t367 + (pkin(4) * mrSges(6,2) - t416) * t366 - t424 * qJDD(4) + (t408 * t392 + (pkin(4) * t364 + t410) * t389) * t382;
t313 = mrSges(4,2) * g(3) - mrSges(4,3) * t341 + Ifges(4,5) * t381 - t380 * Ifges(4,6) - pkin(7) * t326 - t389 * t321 + t392 * t322;
t312 = mrSges(3,2) * t362 - mrSges(2,3) * t374 - pkin(6) * t320 - qJ(2) * t325 + t393 * t313 - t390 * t314 - t415 * t396 + t418 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t311 = mrSges(2,3) * t375 + mrSges(3,2) * t356 - t390 * t313 - t393 * t314 + pkin(2) * t326 - pkin(6) * t405 - pkin(1) * t325 + t418 * t396 + t415 * qJDD(1) + (pkin(2) * m(4) - t420) * g(3);
t1 = [-m(1) * g(1) + t406; -m(1) * g(2) + t411; (-m(1) - m(2) + t423) * g(3) - t326; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t411 - t391 * t311 + t394 * t312; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t406 + t394 * t311 + t391 * t312; pkin(1) * t398 + qJ(2) * (-t396 * mrSges(3,1) + t401) + mrSges(2,1) * t374 - mrSges(2,2) * t375 - pkin(2) * t320 - mrSges(3,1) * t362 + mrSges(3,3) * t356 - t392 * t321 - pkin(3) * t397 - pkin(7) * t404 - mrSges(4,1) * t341 + mrSges(4,2) * t342 - t389 * t322 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - Ifges(4,3) * t381 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
