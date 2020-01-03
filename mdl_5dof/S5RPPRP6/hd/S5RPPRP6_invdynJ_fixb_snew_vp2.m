% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRP6
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:03
% EndTime: 2019-12-31 17:55:05
% DurationCPUTime: 0.75s
% Computational Cost: add. (2449->160), mult. (5383->193), div. (0->0), fcn. (3157->6), ass. (0->76)
t393 = Ifges(5,1) + Ifges(6,1);
t385 = Ifges(5,4) - Ifges(6,5);
t384 = Ifges(5,5) + Ifges(6,4);
t392 = -Ifges(5,2) - Ifges(6,3);
t383 = Ifges(5,6) - Ifges(6,6);
t391 = Ifges(5,3) + Ifges(6,2);
t353 = qJD(1) ^ 2;
t349 = sin(qJ(1));
t351 = cos(qJ(1));
t367 = g(1) * t349 - t351 * g(2);
t357 = -qJ(2) * t353 + qJDD(2) - t367;
t382 = -pkin(1) - qJ(3);
t390 = -(2 * qJD(1) * qJD(3)) + t382 * qJDD(1) + t357;
t347 = cos(pkin(7));
t389 = t347 ^ 2;
t363 = -g(1) * t351 - g(2) * t349;
t388 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t363;
t387 = pkin(3) * t353;
t386 = -mrSges(5,3) - mrSges(6,2);
t381 = t347 * mrSges(4,2);
t346 = sin(pkin(7));
t312 = t346 * g(3) + t390 * t347;
t295 = (-pkin(6) * qJDD(1) - t346 * t387) * t347 + t312;
t313 = -g(3) * t347 + t390 * t346;
t341 = t346 ^ 2;
t370 = t346 * qJDD(1);
t296 = -pkin(6) * t370 - t341 * t387 + t313;
t348 = sin(qJ(4));
t350 = cos(qJ(4));
t292 = t348 * t295 + t350 * t296;
t361 = t346 * t350 + t347 * t348;
t360 = -t346 * t348 + t347 * t350;
t330 = t360 * qJD(1);
t372 = qJD(4) * t330;
t314 = t361 * qJDD(1) + t372;
t323 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t330;
t329 = t361 * qJD(1);
t307 = pkin(4) * t329 - qJ(5) * t330;
t352 = qJD(4) ^ 2;
t289 = -pkin(4) * t352 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t307 * t329 + t292;
t324 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t330;
t368 = m(6) * t289 + qJDD(4) * mrSges(6,3) + qJD(4) * t324;
t308 = mrSges(6,1) * t329 - mrSges(6,3) * t330;
t376 = -mrSges(5,1) * t329 - mrSges(5,2) * t330 - t308;
t282 = m(5) * t292 - qJDD(4) * mrSges(5,2) - qJD(4) * t323 + t386 * t314 + t376 * t329 + t368;
t291 = t295 * t350 - t296 * t348;
t373 = qJD(4) * t329;
t315 = t360 * qJDD(1) - t373;
t322 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t329;
t290 = -qJDD(4) * pkin(4) - qJ(5) * t352 + t307 * t330 + qJDD(5) - t291;
t325 = -mrSges(6,2) * t329 + qJD(4) * mrSges(6,3);
t364 = -m(6) * t290 + qJDD(4) * mrSges(6,1) + qJD(4) * t325;
t283 = m(5) * t291 + qJDD(4) * mrSges(5,1) + qJD(4) * t322 + t386 * t315 + t376 * t330 + t364;
t380 = t348 * t282 + t350 * t283;
t379 = -t391 * qJD(4) + t383 * t329 - t384 * t330;
t378 = t383 * qJD(4) + t392 * t329 + t385 * t330;
t377 = t384 * qJD(4) - t385 * t329 + t393 * t330;
t375 = -t341 - t389;
t366 = t375 * mrSges(4,3);
t365 = t350 * t282 - t283 * t348;
t359 = -qJDD(1) * mrSges(4,3) - t353 * (t346 * mrSges(4,1) + t381);
t362 = (m(4) * t312 + t359 * t347 + t380) * t347 + (m(4) * t313 + t359 * t346 + t365) * t346;
t356 = qJDD(3) + t388;
t298 = pkin(3) * t370 + (t375 * pkin(6) + t382) * t353 + t356;
t287 = -0.2e1 * qJD(5) * t330 + (-t315 + t373) * qJ(5) + (t314 + t372) * pkin(4) + t298;
t284 = m(6) * t287 + t314 * mrSges(6,1) - mrSges(6,3) * t315 - t324 * t330 + t329 * t325;
t355 = m(5) * t298 + mrSges(5,1) * t314 + t315 * mrSges(5,2) + t322 * t329 + t330 * t323 + t284;
t319 = t382 * t353 + t356;
t354 = m(4) * t319 + mrSges(4,1) * t370 + qJDD(1) * t381 + t355;
t328 = -qJDD(1) * pkin(1) + t357;
t327 = pkin(1) * t353 - t388;
t285 = mrSges(6,2) * t315 + t308 * t330 - t364;
t276 = mrSges(5,2) * t298 + mrSges(6,2) * t290 - mrSges(5,3) * t291 - mrSges(6,3) * t287 - qJ(5) * t284 - t378 * qJD(4) + t384 * qJDD(4) - t385 * t314 + t393 * t315 + t379 * t329;
t275 = -mrSges(5,1) * t298 - mrSges(6,1) * t287 + mrSges(6,2) * t289 + mrSges(5,3) * t292 - pkin(4) * t284 + t377 * qJD(4) + t383 * qJDD(4) + t392 * t314 + t385 * t315 + t379 * t330;
t274 = m(3) * t328 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t353 + t362;
t1 = [mrSges(2,1) * t367 - mrSges(2,2) * t363 + mrSges(3,2) * t328 - mrSges(3,3) * t327 + t347 * (mrSges(4,2) * t319 - mrSges(4,3) * t312 - pkin(6) * t380 - t348 * t275 + t350 * t276) - t346 * (-mrSges(4,1) * t319 + mrSges(4,3) * t313 - pkin(3) * t355 + pkin(6) * t365 + t350 * t275 + t348 * t276) - qJ(3) * t362 - pkin(1) * t274 + qJ(2) * (t354 + (mrSges(3,2) + t366) * t353 - m(3) * t327) + (Ifges(4,1) * t389 + qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t347 + Ifges(4,2) * t346) * t346) * qJDD(1); t274; t353 * t366 + t354; mrSges(5,1) * t291 - mrSges(5,2) * t292 - mrSges(6,1) * t290 + mrSges(6,3) * t289 - pkin(4) * t285 + qJ(5) * t368 + t378 * t330 + (-qJ(5) * t308 + t377) * t329 + t384 * t315 + (-qJ(5) * mrSges(6,2) - t383) * t314 + t391 * qJDD(4); t285;];
tauJ = t1;
