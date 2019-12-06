% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_invdynJB_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:47
% EndTime: 2019-12-05 18:08:50
% DurationCPUTime: 0.95s
% Computational Cost: add. (6983->222), mult. (12365->275), div. (0->0), fcn. (8636->8), ass. (0->82)
t380 = sin(qJ(4));
t384 = cos(qJ(4));
t381 = sin(qJ(3));
t398 = qJD(1) * t381;
t365 = qJD(3) * t384 - t380 * t398;
t385 = cos(qJ(3));
t396 = qJD(1) * qJD(3);
t368 = qJDD(1) * t381 + t385 * t396;
t347 = t365 * qJD(4) + qJDD(3) * t380 + t368 * t384;
t366 = qJD(3) * t380 + t384 * t398;
t397 = qJD(1) * t385;
t375 = qJD(4) - t397;
t379 = sin(qJ(5));
t383 = cos(qJ(5));
t350 = t366 * t383 + t379 * t375;
t369 = qJDD(1) * t385 - t381 * t396;
t364 = qJDD(4) - t369;
t330 = -t350 * qJD(5) - t379 * t347 + t364 * t383;
t349 = -t379 * t366 + t375 * t383;
t331 = t349 * qJD(5) + t347 * t383 + t379 * t364;
t382 = sin(qJ(1));
t386 = cos(qJ(1));
t374 = -g(1) * t386 - g(2) * t382;
t358 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t374;
t352 = -g(3) * t381 + t358 * t385;
t373 = g(1) * t382 - g(2) * t386;
t387 = qJD(1) ^ 2;
t363 = -qJ(2) * t387 + qJDD(2) - t373;
t339 = t352 * t384 + t363 * t380;
t351 = g(3) * t385 + t358 * t381;
t333 = t339 * t383 + t379 * t351;
t362 = qJD(5) - t365;
t334 = Ifges(6,5) * t350 + Ifges(6,6) * t349 + Ifges(6,3) * t362;
t336 = Ifges(6,1) * t350 + Ifges(6,4) * t349 + Ifges(6,5) * t362;
t338 = t352 * t380 - t384 * t363;
t346 = -t366 * qJD(4) + qJDD(3) * t384 - t368 * t380;
t345 = qJDD(5) - t346;
t326 = -mrSges(6,1) * t338 + mrSges(6,3) * t333 + Ifges(6,4) * t331 + Ifges(6,2) * t330 + Ifges(6,6) * t345 - t334 * t350 + t336 * t362;
t332 = -t379 * t339 + t351 * t383;
t335 = Ifges(6,4) * t350 + Ifges(6,2) * t349 + Ifges(6,6) * t362;
t327 = mrSges(6,2) * t338 - mrSges(6,3) * t332 + Ifges(6,1) * t331 + Ifges(6,4) * t330 + Ifges(6,5) * t345 + t334 * t349 - t335 * t362;
t342 = Ifges(5,5) * t366 + Ifges(5,6) * t365 + Ifges(5,3) * t375;
t343 = Ifges(5,4) * t366 + Ifges(5,2) * t365 + Ifges(5,6) * t375;
t320 = mrSges(5,2) * t351 + mrSges(5,3) * t338 + Ifges(5,1) * t347 + Ifges(5,4) * t346 + Ifges(5,5) * t364 - t379 * t326 + t327 * t383 + t365 * t342 - t375 * t343;
t344 = Ifges(5,1) * t366 + Ifges(5,4) * t365 + Ifges(5,5) * t375;
t390 = mrSges(6,1) * t332 - mrSges(6,2) * t333 + Ifges(6,5) * t331 + Ifges(6,6) * t330 + Ifges(6,3) * t345 + t335 * t350 - t336 * t349;
t324 = -mrSges(5,1) * t351 + mrSges(5,3) * t339 + Ifges(5,4) * t347 + Ifges(5,2) * t346 + Ifges(5,6) * t364 - t342 * t366 + t344 * t375 - t390;
t360 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t381 + Ifges(4,2) * t385) * qJD(1);
t361 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t381 + Ifges(4,4) * t385) * qJD(1);
t406 = -mrSges(4,1) * t351 - mrSges(4,2) * t352 + Ifges(4,5) * t368 + Ifges(4,6) * t369 + Ifges(4,3) * qJDD(3) + t380 * t320 + t384 * t324 + (t360 * t381 - t361 * t385) * qJD(1);
t337 = -mrSges(6,1) * t349 + mrSges(6,2) * t350;
t340 = -mrSges(6,2) * t362 + mrSges(6,3) * t349;
t328 = m(6) * t332 + mrSges(6,1) * t345 - mrSges(6,3) * t331 - t337 * t350 + t340 * t362;
t341 = mrSges(6,1) * t362 - mrSges(6,3) * t350;
t329 = m(6) * t333 - mrSges(6,2) * t345 + mrSges(6,3) * t330 + t337 * t349 - t341 * t362;
t348 = -mrSges(5,1) * t365 + mrSges(5,2) * t366;
t354 = mrSges(5,1) * t375 - mrSges(5,3) * t366;
t323 = m(5) * t339 - t364 * mrSges(5,2) + t346 * mrSges(5,3) - t379 * t328 + t329 * t383 + t365 * t348 - t375 * t354;
t353 = -mrSges(5,2) * t375 + mrSges(5,3) * t365;
t325 = mrSges(5,1) * t364 + mrSges(6,1) * t330 - mrSges(6,2) * t331 - mrSges(5,3) * t347 + t340 * t349 - t341 * t350 - t348 * t366 + t353 * t375 + (-m(5) - m(6)) * t338;
t371 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t398;
t372 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t397;
t405 = -t369 * mrSges(4,1) + t368 * mrSges(4,2) + (t371 * t381 - t372 * t385) * qJD(1) + t323 * t380 + t325 * t384 + (m(3) + m(4)) * t363;
t403 = mrSges(2,1) + mrSges(3,1);
t402 = -mrSges(2,2) + mrSges(3,3);
t401 = Ifges(3,4) + Ifges(2,5);
t400 = Ifges(2,6) - Ifges(3,6);
t367 = (-mrSges(4,1) * t385 + mrSges(4,2) * t381) * qJD(1);
t317 = m(4) * t352 - qJDD(3) * mrSges(4,2) + t369 * mrSges(4,3) - qJD(3) * t371 + t323 * t384 - t325 * t380 + t367 * t397;
t322 = -t367 * t398 + qJDD(3) * mrSges(4,1) + t346 * mrSges(5,1) - t347 * mrSges(5,2) - t368 * mrSges(4,3) + qJD(3) * t372 - t328 * t383 - t379 * t329 + t365 * t353 - t366 * t354 + (-m(4) - m(5)) * t351;
t399 = t381 * t317 + t385 * t322;
t393 = m(3) * t358 + qJDD(1) * mrSges(3,3) + t385 * t317 - t322 * t381;
t359 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t381 + Ifges(4,6) * t385) * qJD(1);
t313 = mrSges(4,2) * t363 + mrSges(4,3) * t351 + Ifges(4,1) * t368 + Ifges(4,4) * t369 + Ifges(4,5) * qJDD(3) - qJD(3) * t360 + t320 * t384 - t324 * t380 + t359 * t397;
t388 = mrSges(5,1) * t338 + mrSges(5,2) * t339 - Ifges(5,5) * t347 - Ifges(5,6) * t346 - Ifges(5,3) * t364 - t383 * t326 - t379 * t327 - t366 * t343 + t365 * t344;
t319 = -mrSges(4,1) * t363 + mrSges(4,3) * t352 + Ifges(4,4) * t368 + Ifges(4,2) * t369 + Ifges(4,6) * qJDD(3) + qJD(3) * t361 - t359 * t398 + t388;
t391 = -mrSges(3,1) * t363 - mrSges(2,2) * t374 + qJ(2) * (-mrSges(3,1) * t387 + t393) + t381 * t313 + t385 * t319 + mrSges(3,3) * t358 + mrSges(2,1) * t373 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);
t314 = m(2) * t373 + t403 * qJDD(1) + t402 * t387 - t405;
t311 = mrSges(3,2) * t358 + mrSges(2,3) * t374 + t403 * g(3) + t400 * qJDD(1) + t401 * t387 - t406;
t309 = m(2) * t374 - qJDD(1) * mrSges(2,2) - t403 * t387 + t393;
t308 = -mrSges(2,3) * t373 + mrSges(3,2) * t363 + t385 * t313 - t381 * t319 - qJ(2) * t399 - t400 * t387 + t401 * qJDD(1) + (m(3) * qJ(2) + t402) * g(3);
t1 = [-m(1) * g(1) + t309 * t386 - t314 * t382; -m(1) * g(2) + t309 * t382 + t314 * t386; (-m(1) - m(2) - m(3)) * g(3) + t399; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t308 * t386 - t311 * t382; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t308 * t382 + t311 * t386; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t391; t391; -qJDD(1) * mrSges(3,1) - mrSges(3,3) * t387 + t405; t406; -t388; t390;];
tauJB = t1;
