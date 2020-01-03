% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPR9
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR9_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:31
% EndTime: 2019-12-31 17:09:34
% DurationCPUTime: 1.82s
% Computational Cost: add. (17795->237), mult. (38231->304), div. (0->0), fcn. (23887->8), ass. (0->92)
t406 = sin(qJ(1));
t409 = cos(qJ(1));
t397 = -t409 * g(1) - t406 * g(2);
t411 = qJD(1) ^ 2;
t383 = -t411 * pkin(1) + qJDD(1) * pkin(5) + t397;
t405 = sin(qJ(2));
t408 = cos(qJ(2));
t367 = -t405 * g(3) + t408 * t383;
t391 = (-mrSges(3,1) * t408 + mrSges(3,2) * t405) * qJD(1);
t420 = qJD(1) * qJD(2);
t399 = t405 * t420;
t393 = t408 * qJDD(1) - t399;
t422 = qJD(1) * t405;
t394 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t422;
t396 = t406 * g(1) - t409 * g(2);
t382 = -qJDD(1) * pkin(1) - t411 * pkin(5) - t396;
t419 = t408 * t420;
t392 = t405 * qJDD(1) + t419;
t353 = (-t392 - t419) * qJ(3) + (-t393 + t399) * pkin(2) + t382;
t390 = (-pkin(2) * t408 - qJ(3) * t405) * qJD(1);
t410 = qJD(2) ^ 2;
t421 = t408 * qJD(1);
t356 = -t410 * pkin(2) + qJDD(2) * qJ(3) + t390 * t421 + t367;
t402 = sin(pkin(7));
t403 = cos(pkin(7));
t387 = t402 * qJD(2) + t403 * t422;
t340 = -0.2e1 * qJD(3) * t387 + t403 * t353 - t402 * t356;
t372 = t402 * qJDD(2) + t403 * t392;
t386 = t403 * qJD(2) - t402 * t422;
t338 = (-t386 * t421 - t372) * pkin(6) + (t386 * t387 - t393) * pkin(3) + t340;
t341 = 0.2e1 * qJD(3) * t386 + t402 * t353 + t403 * t356;
t371 = t403 * qJDD(2) - t402 * t392;
t373 = -pkin(3) * t421 - t387 * pkin(6);
t385 = t386 ^ 2;
t339 = -t385 * pkin(3) + t371 * pkin(6) + t373 * t421 + t341;
t404 = sin(qJ(4));
t407 = cos(qJ(4));
t336 = t407 * t338 - t404 * t339;
t363 = t407 * t386 - t404 * t387;
t345 = t363 * qJD(4) + t404 * t371 + t407 * t372;
t364 = t404 * t386 + t407 * t387;
t350 = -t363 * mrSges(5,1) + t364 * mrSges(5,2);
t398 = qJD(4) - t421;
t357 = -t398 * mrSges(5,2) + t363 * mrSges(5,3);
t389 = qJDD(4) - t393;
t334 = m(5) * t336 + t389 * mrSges(5,1) - t345 * mrSges(5,3) - t364 * t350 + t398 * t357;
t337 = t404 * t338 + t407 * t339;
t344 = -t364 * qJD(4) + t407 * t371 - t404 * t372;
t358 = t398 * mrSges(5,1) - t364 * mrSges(5,3);
t335 = m(5) * t337 - t389 * mrSges(5,2) + t344 * mrSges(5,3) + t363 * t350 - t398 * t358;
t326 = t407 * t334 + t404 * t335;
t365 = -t386 * mrSges(4,1) + t387 * mrSges(4,2);
t369 = mrSges(4,2) * t421 + t386 * mrSges(4,3);
t324 = m(4) * t340 - t393 * mrSges(4,1) - t372 * mrSges(4,3) - t387 * t365 - t369 * t421 + t326;
t370 = -mrSges(4,1) * t421 - t387 * mrSges(4,3);
t415 = -t404 * t334 + t407 * t335;
t325 = m(4) * t341 + t393 * mrSges(4,2) + t371 * mrSges(4,3) + t386 * t365 + t370 * t421 + t415;
t416 = -t402 * t324 + t403 * t325;
t321 = m(3) * t367 - qJDD(2) * mrSges(3,2) + t393 * mrSges(3,3) - qJD(2) * t394 + t391 * t421 + t416;
t366 = -t408 * g(3) - t405 * t383;
t395 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t421;
t355 = -qJDD(2) * pkin(2) - t410 * qJ(3) + t390 * t422 + qJDD(3) - t366;
t342 = -t371 * pkin(3) - t385 * pkin(6) + t387 * t373 + t355;
t414 = m(5) * t342 - t344 * mrSges(5,1) + t345 * mrSges(5,2) - t363 * t357 + t364 * t358;
t412 = -m(4) * t355 + t371 * mrSges(4,1) - t372 * mrSges(4,2) + t386 * t369 - t387 * t370 - t414;
t330 = m(3) * t366 + qJDD(2) * mrSges(3,1) - t392 * mrSges(3,3) + qJD(2) * t395 - t391 * t422 + t412;
t417 = t408 * t321 - t405 * t330;
t314 = m(2) * t397 - t411 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t417;
t322 = t403 * t324 + t402 * t325;
t413 = -m(3) * t382 + t393 * mrSges(3,1) - t392 * mrSges(3,2) - t394 * t422 + t395 * t421 - t322;
t318 = m(2) * t396 + qJDD(1) * mrSges(2,1) - t411 * mrSges(2,2) + t413;
t423 = t406 * t314 + t409 * t318;
t315 = t405 * t321 + t408 * t330;
t418 = t409 * t314 - t406 * t318;
t381 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t405 + Ifges(3,4) * t408) * qJD(1);
t380 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t405 + Ifges(3,2) * t408) * qJD(1);
t379 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t405 + Ifges(3,6) * t408) * qJD(1);
t361 = Ifges(4,1) * t387 + Ifges(4,4) * t386 - Ifges(4,5) * t421;
t360 = Ifges(4,4) * t387 + Ifges(4,2) * t386 - Ifges(4,6) * t421;
t359 = Ifges(4,5) * t387 + Ifges(4,6) * t386 - Ifges(4,3) * t421;
t348 = Ifges(5,1) * t364 + Ifges(5,4) * t363 + Ifges(5,5) * t398;
t347 = Ifges(5,4) * t364 + Ifges(5,2) * t363 + Ifges(5,6) * t398;
t346 = Ifges(5,5) * t364 + Ifges(5,6) * t363 + Ifges(5,3) * t398;
t328 = mrSges(5,2) * t342 - mrSges(5,3) * t336 + Ifges(5,1) * t345 + Ifges(5,4) * t344 + Ifges(5,5) * t389 + t363 * t346 - t398 * t347;
t327 = -mrSges(5,1) * t342 + mrSges(5,3) * t337 + Ifges(5,4) * t345 + Ifges(5,2) * t344 + Ifges(5,6) * t389 - t364 * t346 + t398 * t348;
t316 = mrSges(4,2) * t355 - mrSges(4,3) * t340 + Ifges(4,1) * t372 + Ifges(4,4) * t371 - Ifges(4,5) * t393 - pkin(6) * t326 - t404 * t327 + t407 * t328 + t386 * t359 + t360 * t421;
t311 = -mrSges(4,1) * t355 + mrSges(4,3) * t341 + Ifges(4,4) * t372 + Ifges(4,2) * t371 - Ifges(4,6) * t393 - pkin(3) * t414 + pkin(6) * t415 + t407 * t327 + t404 * t328 - t387 * t359 - t361 * t421;
t310 = Ifges(3,4) * t392 + Ifges(3,6) * qJDD(2) - t379 * t422 + qJD(2) * t381 - mrSges(3,1) * t382 + mrSges(3,3) * t367 - Ifges(4,5) * t372 - Ifges(4,6) * t371 - t387 * t360 + t386 * t361 - mrSges(4,1) * t340 + mrSges(4,2) * t341 - Ifges(5,5) * t345 - Ifges(5,6) * t344 - Ifges(5,3) * t389 - t364 * t347 + t363 * t348 - mrSges(5,1) * t336 + mrSges(5,2) * t337 - pkin(3) * t326 - pkin(2) * t322 + (Ifges(3,2) + Ifges(4,3)) * t393;
t309 = mrSges(3,2) * t382 - mrSges(3,3) * t366 + Ifges(3,1) * t392 + Ifges(3,4) * t393 + Ifges(3,5) * qJDD(2) - qJ(3) * t322 - qJD(2) * t380 - t402 * t311 + t403 * t316 + t379 * t421;
t308 = Ifges(2,6) * qJDD(1) + t411 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t397 - Ifges(3,5) * t392 - Ifges(3,6) * t393 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t366 + mrSges(3,2) * t367 - t402 * t316 - t403 * t311 - pkin(2) * t412 - qJ(3) * t416 - pkin(1) * t315 + (-t405 * t380 + t408 * t381) * qJD(1);
t307 = -mrSges(2,2) * g(3) - mrSges(2,3) * t396 + Ifges(2,5) * qJDD(1) - t411 * Ifges(2,6) - pkin(5) * t315 + t408 * t309 - t405 * t310;
t1 = [-m(1) * g(1) + t418; -m(1) * g(2) + t423; (-m(1) - m(2)) * g(3) + t315; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t423 + t409 * t307 - t406 * t308; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t418 + t406 * t307 + t409 * t308; -mrSges(1,1) * g(2) + mrSges(2,1) * t396 + mrSges(1,2) * g(1) - mrSges(2,2) * t397 + Ifges(2,3) * qJDD(1) + pkin(1) * t413 + pkin(5) * t417 + t405 * t309 + t408 * t310;];
tauB = t1;
