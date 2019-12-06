% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRP4
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:14
% EndTime: 2019-12-05 15:35:17
% DurationCPUTime: 1.39s
% Computational Cost: add. (10181->193), mult. (17324->238), div. (0->0), fcn. (9126->8), ass. (0->80)
t431 = Ifges(5,1) + Ifges(6,1);
t426 = Ifges(5,4) - Ifges(6,5);
t425 = Ifges(6,4) + Ifges(5,5);
t430 = Ifges(5,2) + Ifges(6,3);
t429 = Ifges(6,6) - Ifges(5,6);
t428 = Ifges(5,3) + Ifges(6,2);
t427 = mrSges(5,3) + mrSges(6,2);
t399 = sin(pkin(7));
t401 = cos(pkin(7));
t385 = t399 * g(1) - t401 * g(2);
t384 = qJDD(3) - t385;
t404 = cos(qJ(4));
t423 = t404 * t384;
t386 = -t401 * g(1) - t399 * g(2);
t397 = -g(3) + qJDD(1);
t403 = sin(qJ(2));
t405 = cos(qJ(2));
t361 = -t403 * t386 + t405 * t397;
t359 = qJDD(2) * pkin(2) + t361;
t362 = t405 * t386 + t403 * t397;
t407 = qJD(2) ^ 2;
t360 = -t407 * pkin(2) + t362;
t398 = sin(pkin(8));
t400 = cos(pkin(8));
t355 = t398 * t359 + t400 * t360;
t353 = -t407 * pkin(3) + qJDD(2) * pkin(6) + t355;
t402 = sin(qJ(4));
t350 = t404 * t353 + t402 * t384;
t379 = (-mrSges(5,1) * t404 + mrSges(5,2) * t402) * qJD(2);
t416 = qJD(2) * qJD(4);
t381 = t404 * qJDD(2) - t402 * t416;
t418 = qJD(2) * t402;
t387 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t418;
t377 = (-pkin(4) * t404 - qJ(5) * t402) * qJD(2);
t406 = qJD(4) ^ 2;
t417 = qJD(2) * t404;
t347 = -t406 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t377 * t417 + t350;
t378 = (-mrSges(6,1) * t404 - mrSges(6,3) * t402) * qJD(2);
t388 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t418;
t411 = m(6) * t347 + qJDD(4) * mrSges(6,3) + qJD(4) * t388 + t378 * t417;
t342 = m(5) * t350 - qJDD(4) * mrSges(5,2) - qJD(4) * t387 + t379 * t417 + t427 * t381 + t411;
t349 = -t402 * t353 + t423;
t380 = t402 * qJDD(2) + t404 * t416;
t389 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t417;
t348 = -qJDD(4) * pkin(4) - t406 * qJ(5) - t423 + qJDD(5) + (qJD(2) * t377 + t353) * t402;
t390 = mrSges(6,2) * t417 + qJD(4) * mrSges(6,3);
t409 = -m(6) * t348 + qJDD(4) * mrSges(6,1) + qJD(4) * t390;
t343 = m(5) * t349 + qJDD(4) * mrSges(5,1) + qJD(4) * t389 - t427 * t380 + (-t378 - t379) * t418 + t409;
t412 = t404 * t342 - t402 * t343;
t333 = m(4) * t355 - t407 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t412;
t354 = t400 * t359 - t398 * t360;
t352 = -qJDD(2) * pkin(3) - t407 * pkin(6) - t354;
t345 = -t381 * pkin(4) - t380 * qJ(5) + (-0.2e1 * qJD(5) * t402 + (pkin(4) * t402 - qJ(5) * t404) * qJD(4)) * qJD(2) + t352;
t344 = m(6) * t345 - t381 * mrSges(6,1) - t380 * mrSges(6,3) - t388 * t418 - t390 * t417;
t408 = -m(5) * t352 + t381 * mrSges(5,1) - t380 * mrSges(5,2) - t387 * t418 + t389 * t417 - t344;
t338 = m(4) * t354 + qJDD(2) * mrSges(4,1) - t407 * mrSges(4,2) + t408;
t328 = t398 * t333 + t400 * t338;
t326 = m(3) * t361 + qJDD(2) * mrSges(3,1) - t407 * mrSges(3,2) + t328;
t413 = t400 * t333 - t398 * t338;
t327 = m(3) * t362 - t407 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t413;
t414 = -t403 * t326 + t405 * t327;
t319 = m(2) * t386 + t414;
t336 = t402 * t342 + t404 * t343;
t410 = m(4) * t384 + t336;
t335 = (m(2) + m(3)) * t385 - t410;
t422 = t399 * t319 + t401 * t335;
t320 = t405 * t326 + t403 * t327;
t421 = t429 * qJD(4) + (-t426 * t402 - t430 * t404) * qJD(2);
t420 = t428 * qJD(4) + (t425 * t402 - t429 * t404) * qJD(2);
t419 = t425 * qJD(4) + (t431 * t402 + t426 * t404) * qJD(2);
t415 = t401 * t319 - t399 * t335;
t330 = mrSges(5,2) * t352 + mrSges(6,2) * t348 - mrSges(5,3) * t349 - mrSges(6,3) * t345 - qJ(5) * t344 + t421 * qJD(4) + t425 * qJDD(4) + t431 * t380 + t426 * t381 + t420 * t417;
t329 = -mrSges(5,1) * t352 - mrSges(6,1) * t345 + mrSges(6,2) * t347 + mrSges(5,3) * t350 - pkin(4) * t344 + t419 * qJD(4) - qJDD(4) * t429 + t426 * t380 + t430 * t381 - t420 * t418;
t322 = Ifges(4,6) * qJDD(2) + t407 * Ifges(4,5) - mrSges(4,1) * t384 + mrSges(4,3) * t355 - mrSges(5,1) * t349 + mrSges(5,2) * t350 + mrSges(6,1) * t348 - mrSges(6,3) * t347 - pkin(4) * t409 - qJ(5) * t411 - pkin(3) * t336 + (-qJ(5) * mrSges(6,2) + t429) * t381 + (pkin(4) * mrSges(6,2) - t425) * t380 - t428 * qJDD(4) + (t419 * t404 + (pkin(4) * t378 + t421) * t402) * qJD(2);
t321 = mrSges(4,2) * t384 - mrSges(4,3) * t354 + Ifges(4,5) * qJDD(2) - t407 * Ifges(4,6) - pkin(6) * t336 - t402 * t329 + t404 * t330;
t316 = -mrSges(3,2) * t385 - mrSges(3,3) * t361 + Ifges(3,5) * qJDD(2) - t407 * Ifges(3,6) - qJ(3) * t328 + t400 * t321 - t398 * t322;
t315 = mrSges(3,1) * t385 + mrSges(3,3) * t362 + t407 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t410 + qJ(3) * t413 + t398 * t321 + t400 * t322;
t314 = -pkin(1) * t320 + mrSges(2,3) * t386 - pkin(2) * t328 - mrSges(3,1) * t361 + mrSges(3,2) * t362 - t404 * t329 - pkin(3) * t408 - pkin(6) * t412 - mrSges(4,1) * t354 + mrSges(4,2) * t355 - t402 * t330 - mrSges(2,1) * t397 + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2);
t313 = mrSges(2,2) * t397 - mrSges(2,3) * t385 - pkin(5) * t320 - t403 * t315 + t405 * t316;
t1 = [-m(1) * g(1) + t415; -m(1) * g(2) + t422; -m(1) * g(3) + m(2) * t397 + t320; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t422 + t401 * t313 - t399 * t314; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t415 + t399 * t313 + t401 * t314; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t385 - mrSges(2,2) * t386 + t403 * t316 + t405 * t315 + pkin(1) * (m(3) * t385 - t410) + pkin(5) * t414;];
tauB = t1;
