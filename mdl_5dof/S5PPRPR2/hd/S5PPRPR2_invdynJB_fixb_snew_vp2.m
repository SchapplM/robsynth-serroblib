% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:06
% EndTime: 2019-12-05 15:03:07
% DurationCPUTime: 0.81s
% Computational Cost: add. (7062->156), mult. (10987->193), div. (0->0), fcn. (6550->8), ass. (0->74)
t394 = sin(pkin(7));
t396 = cos(pkin(7));
t383 = t394 * g(1) - t396 * g(2);
t382 = qJDD(2) - t383;
t384 = -t396 * g(1) - t394 * g(2);
t390 = -g(3) + qJDD(1);
t393 = sin(pkin(8));
t395 = cos(pkin(8));
t366 = -t393 * t384 + t395 * t390;
t367 = t395 * t384 + t393 * t390;
t398 = sin(qJ(3));
t400 = cos(qJ(3));
t361 = t400 * t366 - t398 * t367;
t401 = qJD(3) ^ 2;
t406 = -t401 * qJ(4) + qJDD(4) - t361;
t421 = -pkin(3) - pkin(6);
t358 = t421 * qJDD(3) + t406;
t397 = sin(qJ(5));
t399 = cos(qJ(5));
t354 = t399 * t358 - t397 * t382;
t379 = (mrSges(6,1) * t397 + mrSges(6,2) * t399) * qJD(3);
t412 = qJD(3) * qJD(5);
t381 = t399 * qJDD(3) - t397 * t412;
t414 = qJD(3) * t397;
t385 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t414;
t413 = qJD(3) * t399;
t351 = m(6) * t354 + qJDD(5) * mrSges(6,1) - t381 * mrSges(6,3) + qJD(5) * t385 - t379 * t413;
t355 = t397 * t358 + t399 * t382;
t380 = -t397 * qJDD(3) - t399 * t412;
t386 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t413;
t352 = m(6) * t355 - qJDD(5) * mrSges(6,2) + t380 * mrSges(6,3) - qJD(5) * t386 - t379 * t414;
t415 = -t397 * t351 + t399 * t352;
t341 = m(5) * t382 + t415;
t340 = (m(3) + m(4)) * t382 + t341;
t342 = t399 * t351 + t397 * t352;
t360 = -qJDD(3) * pkin(3) + t406;
t404 = -m(5) * t360 + t401 * mrSges(5,3) - t342;
t337 = qJDD(3) * mrSges(5,2) - t404;
t362 = t398 * t366 + t400 * t367;
t405 = qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t362;
t357 = t421 * t401 + t405;
t370 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t399 - Ifges(6,6) * t397) * qJD(3);
t372 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t399 - Ifges(6,4) * t397) * qJD(3);
t346 = -mrSges(6,1) * t357 + mrSges(6,3) * t355 + Ifges(6,4) * t381 + Ifges(6,2) * t380 + Ifges(6,6) * qJDD(5) + qJD(5) * t372 - t370 * t413;
t371 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t399 - Ifges(6,2) * t397) * qJD(3);
t347 = mrSges(6,2) * t357 - mrSges(6,3) * t354 + Ifges(6,1) * t381 + Ifges(6,4) * t380 + Ifges(6,5) * qJDD(5) - qJD(5) * t371 - t370 * t414;
t359 = t401 * pkin(3) - t405;
t407 = -m(6) * t357 + t380 * mrSges(6,1) - t381 * mrSges(6,2) - t385 * t414 - t386 * t413;
t348 = -m(5) * t359 + t401 * mrSges(5,2) + qJDD(3) * mrSges(5,3) - t407;
t423 = mrSges(4,1) * t361 - mrSges(4,2) * t362 + mrSges(5,2) * t360 - mrSges(5,3) * t359 - pkin(3) * t337 - pkin(6) * t342 + qJ(4) * t348 - t397 * t346 + t399 * t347 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3);
t420 = mrSges(4,1) - mrSges(5,2);
t419 = -Ifges(5,4) + Ifges(4,5);
t418 = Ifges(5,5) - Ifges(4,6);
t336 = m(4) * t361 - t401 * mrSges(4,2) + t420 * qJDD(3) + t404;
t345 = m(4) * t362 - t401 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t348;
t334 = t400 * t336 + t398 * t345;
t332 = m(3) * t366 + t334;
t408 = -t398 * t336 + t400 * t345;
t333 = m(3) * t367 + t408;
t409 = -t393 * t332 + t395 * t333;
t325 = m(2) * t384 + t409;
t339 = m(2) * t383 - t340;
t416 = t394 * t325 + t396 * t339;
t326 = t395 * t332 + t393 * t333;
t411 = m(2) * t390 + t326;
t410 = t396 * t325 - t394 * t339;
t403 = mrSges(6,1) * t354 - mrSges(6,2) * t355 + Ifges(6,5) * t381 + Ifges(6,6) * t380 + Ifges(6,3) * qJDD(5) + t371 * t413 + t372 * t414;
t328 = mrSges(5,1) * t360 - mrSges(4,3) * t361 + pkin(4) * t342 - qJ(4) * t341 + t418 * t401 + (mrSges(4,2) - mrSges(5,3)) * t382 + t419 * qJDD(3) + t403;
t327 = -mrSges(5,1) * t359 + mrSges(4,3) * t362 - pkin(3) * t341 - pkin(4) * t407 - pkin(6) * t415 - t418 * qJDD(3) - t399 * t346 - t397 * t347 - t420 * t382 + t419 * t401;
t322 = mrSges(3,2) * t382 - mrSges(3,3) * t366 - pkin(5) * t334 - t398 * t327 + t400 * t328;
t321 = -mrSges(3,1) * t382 + mrSges(3,3) * t367 + t398 * t328 + t400 * t327 - pkin(2) * (m(4) * t382 + t341) + pkin(5) * t408;
t320 = -mrSges(2,1) * t390 - mrSges(3,1) * t366 + mrSges(3,2) * t367 + mrSges(2,3) * t384 - pkin(1) * t326 - pkin(2) * t334 - t423;
t319 = mrSges(2,2) * t390 - mrSges(2,3) * t383 - qJ(2) * t326 - t393 * t321 + t395 * t322;
t1 = [-m(1) * g(1) + t410; -m(1) * g(2) + t416; -m(1) * g(3) + t411; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t416 + t396 * t319 - t394 * t320; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t410 + t394 * t319 + t396 * t320; -mrSges(1,1) * g(2) + mrSges(2,1) * t383 + mrSges(1,2) * g(1) - mrSges(2,2) * t384 - pkin(1) * t340 + qJ(2) * t409 + t395 * t321 + t393 * t322; t411; t340; t423; t337; t403;];
tauJB = t1;
