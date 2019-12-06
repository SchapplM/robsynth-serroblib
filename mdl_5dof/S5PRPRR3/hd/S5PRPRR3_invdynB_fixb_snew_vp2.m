% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:51
% EndTime: 2019-12-05 15:46:53
% DurationCPUTime: 2.06s
% Computational Cost: add. (23811->217), mult. (42670->279), div. (0->0), fcn. (26262->10), ass. (0->91)
t412 = sin(pkin(8));
t414 = cos(pkin(8));
t400 = -t414 * g(1) - t412 * g(2);
t410 = -g(3) + qJDD(1);
t417 = sin(qJ(2));
t420 = cos(qJ(2));
t383 = -t417 * t400 + t420 * t410;
t379 = qJDD(2) * pkin(2) + t383;
t384 = t420 * t400 + t417 * t410;
t421 = qJD(2) ^ 2;
t380 = -t421 * pkin(2) + t384;
t411 = sin(pkin(9));
t413 = cos(pkin(9));
t367 = t411 * t379 + t413 * t380;
t365 = -t421 * pkin(3) + qJDD(2) * pkin(6) + t367;
t399 = t412 * g(1) - t414 * g(2);
t398 = qJDD(3) - t399;
t416 = sin(qJ(4));
t419 = cos(qJ(4));
t361 = -t416 * t365 + t419 * t398;
t432 = qJD(2) * qJD(4);
t431 = t419 * t432;
t396 = t416 * qJDD(2) + t431;
t358 = (-t396 + t431) * pkin(7) + (t416 * t419 * t421 + qJDD(4)) * pkin(4) + t361;
t362 = t419 * t365 + t416 * t398;
t397 = t419 * qJDD(2) - t416 * t432;
t434 = qJD(2) * t416;
t403 = qJD(4) * pkin(4) - pkin(7) * t434;
t409 = t419 ^ 2;
t359 = -t409 * t421 * pkin(4) + t397 * pkin(7) - qJD(4) * t403 + t362;
t415 = sin(qJ(5));
t418 = cos(qJ(5));
t356 = t418 * t358 - t415 * t359;
t388 = (-t415 * t416 + t418 * t419) * qJD(2);
t370 = t388 * qJD(5) + t418 * t396 + t415 * t397;
t389 = (t415 * t419 + t416 * t418) * qJD(2);
t375 = -t388 * mrSges(6,1) + t389 * mrSges(6,2);
t408 = qJD(4) + qJD(5);
t381 = -t408 * mrSges(6,2) + t388 * mrSges(6,3);
t407 = qJDD(4) + qJDD(5);
t354 = m(6) * t356 + t407 * mrSges(6,1) - t370 * mrSges(6,3) - t389 * t375 + t408 * t381;
t357 = t415 * t358 + t418 * t359;
t369 = -t389 * qJD(5) - t415 * t396 + t418 * t397;
t382 = t408 * mrSges(6,1) - t389 * mrSges(6,3);
t355 = m(6) * t357 - t407 * mrSges(6,2) + t369 * mrSges(6,3) + t388 * t375 - t408 * t382;
t346 = t418 * t354 + t415 * t355;
t395 = (-mrSges(5,1) * t419 + mrSges(5,2) * t416) * qJD(2);
t433 = qJD(2) * t419;
t402 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t433;
t344 = m(5) * t361 + qJDD(4) * mrSges(5,1) - t396 * mrSges(5,3) + qJD(4) * t402 - t395 * t434 + t346;
t401 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t434;
t426 = -t415 * t354 + t418 * t355;
t345 = m(5) * t362 - qJDD(4) * mrSges(5,2) + t397 * mrSges(5,3) - qJD(4) * t401 + t395 * t433 + t426;
t427 = -t416 * t344 + t419 * t345;
t337 = m(4) * t367 - t421 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t427;
t366 = t413 * t379 - t411 * t380;
t424 = -qJDD(2) * pkin(3) - t366;
t364 = -t421 * pkin(6) + t424;
t360 = t403 * t434 - t397 * pkin(4) + (-pkin(7) * t409 - pkin(6)) * t421 + t424;
t423 = m(6) * t360 - t369 * mrSges(6,1) + t370 * mrSges(6,2) - t388 * t381 + t389 * t382;
t422 = -m(5) * t364 + t397 * mrSges(5,1) - t396 * mrSges(5,2) - t401 * t434 + t402 * t433 - t423;
t350 = m(4) * t366 + qJDD(2) * mrSges(4,1) - t421 * mrSges(4,2) + t422;
t333 = t411 * t337 + t413 * t350;
t331 = m(3) * t383 + qJDD(2) * mrSges(3,1) - t421 * mrSges(3,2) + t333;
t428 = t413 * t337 - t411 * t350;
t332 = m(3) * t384 - t421 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t428;
t429 = -t417 * t331 + t420 * t332;
t324 = m(2) * t400 + t429;
t340 = t419 * t344 + t416 * t345;
t425 = m(4) * t398 + t340;
t339 = (m(2) + m(3)) * t399 - t425;
t435 = t412 * t324 + t414 * t339;
t325 = t420 * t331 + t417 * t332;
t430 = t414 * t324 - t412 * t339;
t387 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t416 + Ifges(5,4) * t419) * qJD(2);
t386 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t416 + Ifges(5,2) * t419) * qJD(2);
t385 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t416 + Ifges(5,6) * t419) * qJD(2);
t373 = Ifges(6,1) * t389 + Ifges(6,4) * t388 + Ifges(6,5) * t408;
t372 = Ifges(6,4) * t389 + Ifges(6,2) * t388 + Ifges(6,6) * t408;
t371 = Ifges(6,5) * t389 + Ifges(6,6) * t388 + Ifges(6,3) * t408;
t348 = mrSges(6,2) * t360 - mrSges(6,3) * t356 + Ifges(6,1) * t370 + Ifges(6,4) * t369 + Ifges(6,5) * t407 + t388 * t371 - t408 * t372;
t347 = -mrSges(6,1) * t360 + mrSges(6,3) * t357 + Ifges(6,4) * t370 + Ifges(6,2) * t369 + Ifges(6,6) * t407 - t389 * t371 + t408 * t373;
t334 = mrSges(5,2) * t364 - mrSges(5,3) * t361 + Ifges(5,1) * t396 + Ifges(5,4) * t397 + Ifges(5,5) * qJDD(4) - pkin(7) * t346 - qJD(4) * t386 - t415 * t347 + t418 * t348 + t385 * t433;
t327 = -mrSges(5,1) * t364 + mrSges(5,3) * t362 + Ifges(5,4) * t396 + Ifges(5,2) * t397 + Ifges(5,6) * qJDD(4) - pkin(4) * t423 + pkin(7) * t426 + qJD(4) * t387 + t418 * t347 + t415 * t348 - t385 * t434;
t326 = Ifges(4,6) * qJDD(2) + t421 * Ifges(4,5) - mrSges(4,1) * t398 + mrSges(4,3) * t367 - Ifges(5,5) * t396 - Ifges(5,6) * t397 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t361 + mrSges(5,2) * t362 - Ifges(6,5) * t370 - Ifges(6,6) * t369 - Ifges(6,3) * t407 - t389 * t372 + t388 * t373 - mrSges(6,1) * t356 + mrSges(6,2) * t357 - pkin(4) * t346 - pkin(3) * t340 + (-t416 * t386 + t419 * t387) * qJD(2);
t321 = mrSges(4,2) * t398 - mrSges(4,3) * t366 + Ifges(4,5) * qJDD(2) - t421 * Ifges(4,6) - pkin(6) * t340 - t416 * t327 + t419 * t334;
t320 = -mrSges(3,2) * t399 - mrSges(3,3) * t383 + Ifges(3,5) * qJDD(2) - t421 * Ifges(3,6) - qJ(3) * t333 + t413 * t321 - t411 * t326;
t319 = -pkin(1) * t325 + mrSges(2,3) * t400 - pkin(2) * t333 - mrSges(3,1) * t383 + mrSges(3,2) * t384 - t416 * t334 - t419 * t327 - pkin(3) * t422 - pkin(6) * t427 - mrSges(4,1) * t366 + mrSges(4,2) * t367 - mrSges(2,1) * t410 + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2);
t318 = mrSges(3,1) * t399 + mrSges(3,3) * t384 + t421 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t425 + qJ(3) * t428 + t411 * t321 + t413 * t326;
t317 = mrSges(2,2) * t410 - mrSges(2,3) * t399 - pkin(5) * t325 - t417 * t318 + t420 * t320;
t1 = [-m(1) * g(1) + t430; -m(1) * g(2) + t435; -m(1) * g(3) + m(2) * t410 + t325; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t435 + t414 * t317 - t412 * t319; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t430 + t412 * t317 + t414 * t319; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t399 - mrSges(2,2) * t400 + t417 * t320 + t420 * t318 + pkin(1) * (m(3) * t399 - t425) + pkin(5) * t429;];
tauB = t1;
