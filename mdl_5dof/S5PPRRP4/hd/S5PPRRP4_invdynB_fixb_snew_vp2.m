% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRRP4
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:30
% EndTime: 2019-12-31 17:34:31
% DurationCPUTime: 0.87s
% Computational Cost: add. (4691->181), mult. (8520->214), div. (0->0), fcn. (4214->6), ass. (0->78)
t383 = Ifges(5,1) + Ifges(6,1);
t375 = Ifges(5,4) + Ifges(6,4);
t374 = Ifges(5,5) + Ifges(6,5);
t382 = Ifges(5,2) + Ifges(6,2);
t381 = Ifges(5,6) + Ifges(6,6);
t380 = Ifges(5,3) + Ifges(6,3);
t351 = qJD(3) ^ 2;
t345 = sin(pkin(7));
t346 = cos(pkin(7));
t332 = t345 * g(1) - t346 * g(2);
t330 = qJDD(2) - t332;
t333 = -t346 * g(1) - t345 * g(2);
t348 = sin(qJ(3));
t350 = cos(qJ(3));
t308 = t350 * t330 - t348 * t333;
t355 = -qJDD(3) * pkin(3) - t308;
t306 = -t351 * pkin(6) + t355;
t347 = sin(qJ(4));
t349 = cos(qJ(4));
t365 = qJD(3) * qJD(4);
t361 = t349 * t365;
t327 = t347 * qJDD(3) + t361;
t328 = t349 * qJDD(3) - t347 * t365;
t366 = qJD(3) * t349;
t338 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t366;
t367 = qJD(3) * t347;
t334 = qJD(4) * pkin(4) - qJ(5) * t367;
t343 = t349 ^ 2;
t302 = t334 * t367 - t328 * pkin(4) + qJDD(5) + (-qJ(5) * t343 - pkin(6)) * t351 + t355;
t337 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t366;
t357 = m(6) * t302 - t328 * mrSges(6,1) - t337 * t366;
t376 = -mrSges(5,2) - mrSges(6,2);
t379 = -m(5) * t306 + t328 * mrSges(5,1) + t376 * t327 + t338 * t366 - t357;
t378 = pkin(4) * t351;
t377 = -mrSges(2,2) + mrSges(3,3);
t309 = t348 * t330 + t350 * t333;
t307 = -t351 * pkin(3) + qJDD(3) * pkin(6) + t309;
t344 = g(3) - qJDD(1);
t304 = t349 * t307 + t347 * t344;
t326 = (-mrSges(5,1) * t349 + mrSges(5,2) * t347) * qJD(3);
t364 = qJD(3) * qJD(5);
t301 = t328 * qJ(5) - qJD(4) * t334 - t343 * t378 + 0.2e1 * t349 * t364 + t304;
t325 = (-mrSges(6,1) * t349 + mrSges(6,2) * t347) * qJD(3);
t362 = m(6) * t301 + t328 * mrSges(6,3) + t325 * t366;
t335 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t367;
t368 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t367 - t335;
t296 = m(5) * t304 + t328 * mrSges(5,3) + t368 * qJD(4) + t376 * qJDD(4) + t326 * t366 + t362;
t294 = t349 * t296;
t340 = t349 * t344;
t303 = -t347 * t307 + t340;
t300 = qJDD(4) * pkin(4) + t340 + (-t327 + t361) * qJ(5) + (t349 * t378 - t307 - 0.2e1 * t364) * t347;
t363 = m(6) * t300 + qJDD(4) * mrSges(6,1) + qJD(4) * t337;
t295 = m(5) * t303 + qJDD(4) * mrSges(5,1) + qJD(4) * t338 + (-mrSges(5,3) - mrSges(6,3)) * t327 + (-t325 - t326) * t367 + t363;
t290 = m(4) * t309 - t351 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t347 * t295 + t294;
t358 = qJD(3) * t368;
t293 = m(4) * t308 + qJDD(3) * mrSges(4,1) - t351 * mrSges(4,2) + t347 * t358 + t379;
t286 = t348 * t290 + t350 * t293;
t353 = -m(3) * t330 - t286;
t284 = m(2) * t332 + t353;
t359 = t350 * t290 - t348 * t293;
t356 = m(3) * t333 + t359;
t285 = m(2) * t333 + t356;
t372 = t346 * t284 + t345 * t285;
t371 = t380 * qJD(4) + (t374 * t347 + t381 * t349) * qJD(3);
t370 = -t381 * qJD(4) + (-t375 * t347 - t382 * t349) * qJD(3);
t369 = t374 * qJD(4) + (t383 * t347 + t375 * t349) * qJD(3);
t360 = -t345 * t284 + t346 * t285;
t292 = t349 * t295 + t347 * t296;
t354 = -m(3) * t344 - t292;
t297 = -t327 * mrSges(6,3) - t325 * t367 + t363;
t291 = -m(4) * t344 + t354;
t288 = mrSges(5,2) * t306 + mrSges(6,2) * t302 - mrSges(5,3) * t303 - mrSges(6,3) * t300 - qJ(5) * t297 + t370 * qJD(4) + t374 * qJDD(4) + t383 * t327 + t375 * t328 + t371 * t366;
t287 = -mrSges(5,1) * t306 + mrSges(5,3) * t304 - mrSges(6,1) * t302 + mrSges(6,3) * t301 - pkin(4) * t357 + qJ(5) * t362 + t382 * t328 + (-pkin(4) * mrSges(6,2) + t375) * t327 + (-qJ(5) * mrSges(6,2) + t381) * qJDD(4) + (-qJ(5) * t335 + t369) * qJD(4) + (-pkin(4) * t335 - t371) * t367;
t280 = -mrSges(4,1) * t344 - mrSges(5,1) * t303 - mrSges(6,1) * t300 + mrSges(5,2) * t304 + mrSges(6,2) * t301 + mrSges(4,3) * t309 + t351 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t292 - pkin(4) * t297 - t381 * t328 - t374 * t327 - t380 * qJDD(4) + (t370 * t347 + t369 * t349) * qJD(3);
t279 = mrSges(4,2) * t344 - mrSges(4,3) * t308 + Ifges(4,5) * qJDD(3) - t351 * Ifges(4,6) - pkin(6) * t292 - t347 * t287 + t349 * t288;
t278 = mrSges(3,2) * t330 - mrSges(2,3) * t332 - pkin(5) * t286 - qJ(2) * t291 + t350 * t279 - t348 * t280 + t377 * t344;
t277 = -t348 * t279 - t350 * t280 + pkin(2) * t292 - pkin(5) * t359 - pkin(1) * t291 + (pkin(2) * m(4) + mrSges(2,1) + mrSges(3,1)) * t344 + (mrSges(3,2) + mrSges(2,3)) * t333;
t1 = [-m(1) * g(1) + t360; -m(1) * g(2) + t372; -m(1) * g(3) + (-m(2) - m(4)) * t344 + t354; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t372 - t345 * t277 + t346 * t278; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t360 + t346 * t277 + t345 * t278; pkin(1) * t353 + qJ(2) * t356 + mrSges(2,1) * t332 - pkin(2) * t286 - mrSges(3,1) * t330 - t349 * t287 - pkin(3) * t379 - pkin(6) * t294 - mrSges(4,1) * t308 + mrSges(4,2) * t309 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - Ifges(4,3) * qJDD(3) + (-pkin(3) * t358 + pkin(6) * t295 - t288) * t347 + t377 * t333;];
tauB = t1;
