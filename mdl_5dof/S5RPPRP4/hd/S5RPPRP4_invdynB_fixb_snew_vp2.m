% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:08
% EndTime: 2019-12-31 17:52:10
% DurationCPUTime: 1.10s
% Computational Cost: add. (6668->205), mult. (11645->235), div. (0->0), fcn. (4214->6), ass. (0->86)
t439 = Ifges(5,1) + Ifges(6,1);
t429 = Ifges(5,4) + Ifges(6,4);
t428 = Ifges(5,5) + Ifges(6,5);
t438 = Ifges(5,2) + Ifges(6,2);
t437 = Ifges(5,6) + Ifges(6,6);
t436 = (Ifges(5,3) + Ifges(6,3));
t402 = qJD(1) ^ 2;
t399 = sin(qJ(1));
t401 = cos(qJ(1));
t383 = -t401 * g(1) - t399 * g(2);
t407 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t383;
t434 = -pkin(1) - pkin(2);
t355 = t434 * t402 + t407;
t382 = t399 * g(1) - t401 * g(2);
t406 = -t402 * qJ(2) + qJDD(2) - t382;
t357 = t434 * qJDD(1) + t406;
t396 = sin(pkin(7));
t397 = cos(pkin(7));
t350 = -t396 * t355 + t397 * t357;
t409 = qJDD(1) * pkin(3) - t350;
t348 = -t402 * pkin(6) + t409;
t398 = sin(qJ(4));
t400 = cos(qJ(4));
t418 = qJD(1) * qJD(4);
t416 = t400 * t418;
t374 = -t398 * qJDD(1) - t416;
t375 = -t400 * qJDD(1) + t398 * t418;
t420 = qJD(1) * t398;
t379 = (qJD(4) * mrSges(5,1)) + mrSges(5,3) * t420;
t377 = (qJD(4) * pkin(4)) + qJ(5) * t420;
t392 = t400 ^ 2;
t344 = -t377 * t420 - t375 * pkin(4) + qJDD(5) + (-qJ(5) * t392 - pkin(6)) * t402 + t409;
t378 = (qJD(4) * mrSges(6,1)) + mrSges(6,3) * t420;
t410 = m(6) * t344 - t375 * mrSges(6,1) - t378 * t420;
t431 = -mrSges(5,2) - mrSges(6,2);
t435 = -m(5) * t348 + t375 * mrSges(5,1) + t431 * t374 + t379 * t420 - t410;
t433 = pkin(4) * t402;
t432 = mrSges(2,1) + mrSges(3,1);
t430 = Ifges(3,4) + Ifges(2,5);
t427 = Ifges(2,6) - Ifges(3,6);
t358 = -t402 * pkin(1) + t407;
t351 = t397 * t355 + t396 * t357;
t349 = -t402 * pkin(3) - qJDD(1) * pkin(6) + t351;
t393 = g(3) + qJDD(3);
t385 = t400 * t393;
t345 = -t398 * t349 + t385;
t373 = (mrSges(5,1) * t400 - mrSges(5,2) * t398) * qJD(1);
t419 = qJD(1) * t400;
t381 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t419;
t417 = qJD(1) * qJD(5);
t342 = (qJDD(4) * pkin(4)) + t385 + (-t374 - t416) * qJ(5) + (t400 * t433 - t349 + (2 * t417)) * t398;
t372 = (mrSges(6,1) * t400 - mrSges(6,2) * t398) * qJD(1);
t380 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t419;
t411 = m(6) * t342 + (qJDD(4) * mrSges(6,1)) + qJD(4) * t380 + t372 * t420;
t337 = t373 * t420 + m(5) * t345 + (qJDD(4) * mrSges(5,1)) + qJD(4) * t381 + (-mrSges(5,3) - mrSges(6,3)) * t374 + t411;
t346 = t400 * t349 + t398 * t393;
t343 = t375 * qJ(5) - qJD(4) * t377 - t392 * t433 - 0.2e1 * t400 * t417 + t346;
t424 = m(6) * t343 + t375 * mrSges(6,3);
t338 = m(5) * t346 + t375 * mrSges(5,3) + t431 * qJDD(4) + (-t378 - t379) * qJD(4) + (-t372 - t373) * t419 + t424;
t413 = -t398 * t337 + t400 * t338;
t332 = m(4) * t351 - t402 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t413;
t412 = qJD(1) * (-t380 - t381);
t335 = m(4) * t350 - qJDD(1) * mrSges(4,1) - t402 * mrSges(4,2) + t400 * t412 + t435;
t414 = t397 * t332 - t396 * t335;
t408 = m(3) * t358 + qJDD(1) * mrSges(3,3) + t414;
t326 = m(2) * t383 - qJDD(1) * mrSges(2,2) - t432 * t402 + t408;
t328 = t396 * t332 + t397 * t335;
t359 = -qJDD(1) * pkin(1) + t406;
t403 = -m(3) * t359 + qJDD(1) * mrSges(3,1) + t402 * mrSges(3,3) - t328;
t327 = m(2) * t382 + qJDD(1) * mrSges(2,1) - t402 * mrSges(2,2) + t403;
t425 = t399 * t326 + t401 * t327;
t423 = (t436 * qJD(4)) + (-t428 * t398 - t437 * t400) * qJD(1);
t422 = t437 * qJD(4) + (-t429 * t398 - t438 * t400) * qJD(1);
t421 = -t428 * qJD(4) + (t439 * t398 + t429 * t400) * qJD(1);
t415 = t401 * t326 - t399 * t327;
t334 = t400 * t337 + t398 * t338;
t405 = -m(4) * t393 - t334;
t339 = -t374 * mrSges(6,3) + t411;
t333 = -m(3) * g(3) + t405;
t330 = mrSges(5,2) * t348 + mrSges(6,2) * t344 - mrSges(5,3) * t345 - mrSges(6,3) * t342 - qJ(5) * t339 - t422 * qJD(4) + t428 * qJDD(4) + t439 * t374 + t429 * t375 - t423 * t419;
t329 = -mrSges(5,1) * t348 + mrSges(5,3) * t346 - mrSges(6,1) * t344 + mrSges(6,3) * t343 - pkin(4) * t410 + qJ(5) * t424 + t438 * t375 + (-pkin(4) * mrSges(6,2) + t429) * t374 + (-qJ(5) * mrSges(6,2) + t437) * qJDD(4) + (-qJ(5) * t378 - t421) * qJD(4) + ((-pkin(4) * t380 - qJ(5) * t372) * t400 + t423 * t398) * qJD(1);
t322 = -mrSges(4,1) * t393 - mrSges(5,1) * t345 - mrSges(6,1) * t342 + mrSges(5,2) * t346 + mrSges(6,2) * t343 + mrSges(4,3) * t351 + t402 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t334 - pkin(4) * t339 - t437 * t375 - t428 * t374 - (t436 * qJDD(4)) + (t422 * t398 + t421 * t400) * qJD(1);
t321 = mrSges(4,2) * t393 - mrSges(4,3) * t350 - Ifges(4,5) * qJDD(1) - t402 * Ifges(4,6) - pkin(6) * t334 - t398 * t329 + t400 * t330;
t320 = mrSges(3,2) * t359 - mrSges(2,3) * t382 - qJ(2) * t333 - qJ(3) * t328 + t397 * t321 - t396 * t322 - t427 * t402 + t430 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t319 = mrSges(3,2) * t358 + mrSges(2,3) * t383 - pkin(1) * t333 - pkin(2) * t405 + t432 * g(3) - qJ(3) * t414 + t427 * qJDD(1) - t396 * t321 - t397 * t322 + t430 * t402;
t1 = [-m(1) * g(1) + t415; -m(1) * g(2) + t425; (-m(1) - m(2) - m(3)) * g(3) + t405; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t425 - t399 * t319 + t401 * t320; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t415 + t401 * t319 + t399 * t320; qJ(2) * (-t402 * mrSges(3,1) + t408) + pkin(1) * t403 + mrSges(2,1) * t382 - mrSges(2,2) * t383 - pkin(2) * t328 - mrSges(3,1) * t359 + mrSges(3,3) * t358 - pkin(3) * t435 - pkin(6) * t413 + mrSges(4,2) * t351 - t398 * t330 - mrSges(4,1) * t350 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(3) * t412 - t329) * t400 + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);];
tauB = t1;
