% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPR4
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR4_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR4_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:30
% EndTime: 2019-12-31 17:02:31
% DurationCPUTime: 1.06s
% Computational Cost: add. (12881->176), mult. (17882->224), div. (0->0), fcn. (10403->8), ass. (0->80)
t351 = qJD(1) + qJD(2);
t347 = t351 ^ 2;
t353 = cos(pkin(7));
t381 = pkin(3) * t353;
t352 = sin(pkin(7));
t380 = mrSges(4,2) * t352;
t350 = t353 ^ 2;
t379 = t347 * t350;
t348 = qJDD(1) + qJDD(2);
t378 = t348 * t353;
t367 = Ifges(4,5) * t352 + Ifges(4,6) * t353;
t377 = t347 * t367;
t356 = sin(qJ(1));
t359 = cos(qJ(1));
t341 = t356 * g(1) - t359 * g(2);
t337 = qJDD(1) * pkin(1) + t341;
t342 = -t359 * g(1) - t356 * g(2);
t360 = qJD(1) ^ 2;
t338 = -t360 * pkin(1) + t342;
t355 = sin(qJ(2));
t358 = cos(qJ(2));
t327 = t355 * t337 + t358 * t338;
t325 = -t347 * pkin(2) + t348 * qJ(3) + t327;
t375 = qJD(3) * t351;
t374 = -t353 * g(3) - 0.2e1 * t352 * t375;
t310 = (-pkin(6) * t348 + t347 * t381 - t325) * t352 + t374;
t314 = -t352 * g(3) + (t325 + 0.2e1 * t375) * t353;
t311 = -pkin(3) * t379 + pkin(6) * t378 + t314;
t354 = sin(qJ(4));
t357 = cos(qJ(4));
t308 = t357 * t310 - t354 * t311;
t363 = -t352 * t354 + t353 * t357;
t330 = t363 * t351;
t364 = t352 * t357 + t353 * t354;
t331 = t364 * t351;
t321 = -t330 * mrSges(5,1) + t331 * mrSges(5,2);
t324 = t330 * qJD(4) + t364 * t348;
t328 = -qJD(4) * mrSges(5,2) + t330 * mrSges(5,3);
t306 = m(5) * t308 + qJDD(4) * mrSges(5,1) - t324 * mrSges(5,3) + qJD(4) * t328 - t331 * t321;
t309 = t354 * t310 + t357 * t311;
t323 = -t331 * qJD(4) + t363 * t348;
t329 = qJD(4) * mrSges(5,1) - t331 * mrSges(5,3);
t307 = m(5) * t309 - qJDD(4) * mrSges(5,2) + t323 * mrSges(5,3) - qJD(4) * t329 + t330 * t321;
t298 = t357 * t306 + t354 * t307;
t313 = -t352 * t325 + t374;
t365 = mrSges(4,3) * t348 + (-mrSges(4,1) * t353 + t380) * t347;
t296 = m(4) * t313 - t365 * t352 + t298;
t370 = -t354 * t306 + t357 * t307;
t297 = m(4) * t314 + t365 * t353 + t370;
t371 = -t352 * t296 + t353 * t297;
t291 = m(3) * t327 - t347 * mrSges(3,1) - t348 * mrSges(3,2) + t371;
t326 = t358 * t337 - t355 * t338;
t366 = qJDD(3) - t326;
t322 = -t348 * pkin(2) - t347 * qJ(3) + t366;
t349 = t352 ^ 2;
t312 = (-pkin(2) - t381) * t348 + (-qJ(3) + (-t349 - t350) * pkin(6)) * t347 + t366;
t362 = m(5) * t312 - t323 * mrSges(5,1) + t324 * mrSges(5,2) - t330 * t328 + t331 * t329;
t361 = -m(4) * t322 + mrSges(4,1) * t378 - t362 + (t347 * t349 + t379) * mrSges(4,3);
t302 = m(3) * t326 - t347 * mrSges(3,2) + (mrSges(3,1) - t380) * t348 + t361;
t287 = t355 * t291 + t358 * t302;
t285 = m(2) * t341 + qJDD(1) * mrSges(2,1) - t360 * mrSges(2,2) + t287;
t372 = t358 * t291 - t355 * t302;
t286 = m(2) * t342 - t360 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t372;
t376 = t359 * t285 + t356 * t286;
t292 = t353 * t296 + t352 * t297;
t373 = -t356 * t285 + t359 * t286;
t369 = Ifges(4,1) * t352 + Ifges(4,4) * t353;
t368 = Ifges(4,4) * t352 + Ifges(4,2) * t353;
t317 = Ifges(5,1) * t331 + Ifges(5,4) * t330 + Ifges(5,5) * qJD(4);
t316 = Ifges(5,4) * t331 + Ifges(5,2) * t330 + Ifges(5,6) * qJD(4);
t315 = Ifges(5,5) * t331 + Ifges(5,6) * t330 + Ifges(5,3) * qJD(4);
t300 = mrSges(5,2) * t312 - mrSges(5,3) * t308 + Ifges(5,1) * t324 + Ifges(5,4) * t323 + Ifges(5,5) * qJDD(4) - qJD(4) * t316 + t330 * t315;
t299 = -mrSges(5,1) * t312 + mrSges(5,3) * t309 + Ifges(5,4) * t324 + Ifges(5,2) * t323 + Ifges(5,6) * qJDD(4) + qJD(4) * t317 - t331 * t315;
t288 = mrSges(4,2) * t322 - mrSges(4,3) * t313 - pkin(6) * t298 - t354 * t299 + t357 * t300 + t369 * t348 + t353 * t377;
t281 = -mrSges(4,1) * t322 + mrSges(4,3) * t314 - pkin(3) * t362 + pkin(6) * t370 + t357 * t299 + t354 * t300 + t368 * t348 - t352 * t377;
t280 = mrSges(3,1) * g(3) - mrSges(4,1) * t313 - mrSges(5,1) * t308 + mrSges(4,2) * t314 + mrSges(5,2) * t309 + mrSges(3,3) * t327 - Ifges(5,5) * t324 - Ifges(5,6) * t323 - Ifges(5,3) * qJDD(4) - pkin(2) * t292 - pkin(3) * t298 - t331 * t316 + t330 * t317 + (Ifges(3,6) - t367) * t348 + (-t352 * t368 + t353 * t369 + Ifges(3,5)) * t347;
t279 = -mrSges(3,2) * g(3) - mrSges(3,3) * t326 + Ifges(3,5) * t348 - t347 * Ifges(3,6) - qJ(3) * t292 - t352 * t281 + t353 * t288;
t278 = -mrSges(2,2) * g(3) - mrSges(2,3) * t341 + Ifges(2,5) * qJDD(1) - t360 * Ifges(2,6) - pkin(5) * t287 + t358 * t279 - t355 * t280;
t277 = Ifges(2,6) * qJDD(1) + t360 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t342 + t355 * t279 + t358 * t280 - pkin(1) * (-m(3) * g(3) + t292) + pkin(5) * t372;
t1 = [-m(1) * g(1) + t373; -m(1) * g(2) + t376; (-m(1) - m(2) - m(3)) * g(3) + t292; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t376 - t356 * t277 + t359 * t278; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t373 + t359 * t277 + t356 * t278; pkin(1) * t287 + mrSges(2,1) * t341 - mrSges(2,2) * t342 + t352 * t288 + t353 * t281 + pkin(2) * (-t348 * t380 + t361) + qJ(3) * t371 + mrSges(3,1) * t326 - mrSges(3,2) * t327 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(3,3) * t348 + Ifges(2,3) * qJDD(1);];
tauB = t1;
