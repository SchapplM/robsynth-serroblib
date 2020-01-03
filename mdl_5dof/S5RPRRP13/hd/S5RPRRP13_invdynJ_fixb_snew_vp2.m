% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP13
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP13_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:39
% EndTime: 2019-12-31 18:58:41
% DurationCPUTime: 0.97s
% Computational Cost: add. (3234->191), mult. (6060->230), div. (0->0), fcn. (3279->6), ass. (0->79)
t379 = Ifges(5,1) + Ifges(6,1);
t367 = Ifges(5,4) - Ifges(6,5);
t375 = -Ifges(5,5) - Ifges(6,4);
t378 = Ifges(5,2) + Ifges(6,3);
t365 = Ifges(5,6) - Ifges(6,6);
t338 = sin(qJ(4));
t341 = cos(qJ(3));
t360 = qJD(1) * t341;
t369 = cos(qJ(4));
t323 = -qJD(3) * t369 + t338 * t360;
t339 = sin(qJ(3));
t358 = qJD(1) * qJD(3);
t355 = t339 * t358;
t328 = t341 * qJDD(1) - t355;
t299 = -t323 * qJD(4) + t338 * qJDD(3) + t328 * t369;
t324 = t338 * qJD(3) + t360 * t369;
t303 = t323 * mrSges(6,1) - t324 * mrSges(6,3);
t344 = qJD(1) ^ 2;
t370 = -pkin(1) - pkin(6);
t340 = sin(qJ(1));
t342 = cos(qJ(1));
t350 = -t342 * g(1) - t340 * g(2);
t371 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t350;
t312 = t344 * t370 - t371;
t354 = t341 * t358;
t327 = -t339 * qJDD(1) - t354;
t284 = (-t328 + t355) * pkin(7) + (-t327 + t354) * pkin(3) + t312;
t353 = t340 * g(1) - t342 * g(2);
t347 = -t344 * qJ(2) + qJDD(2) - t353;
t313 = qJDD(1) * t370 + t347;
t306 = -t341 * g(3) + t339 * t313;
t326 = (t339 * pkin(3) - t341 * pkin(7)) * qJD(1);
t343 = qJD(3) ^ 2;
t359 = t339 * qJD(1);
t287 = -t343 * pkin(3) + qJDD(3) * pkin(7) - t326 * t359 + t306;
t281 = t284 * t369 - t338 * t287;
t302 = t323 * pkin(4) - t324 * qJ(5);
t322 = qJDD(4) - t327;
t332 = qJD(4) + t359;
t331 = t332 ^ 2;
t279 = -t322 * pkin(4) - t331 * qJ(5) + t324 * t302 + qJDD(5) - t281;
t310 = -t323 * mrSges(6,2) + t332 * mrSges(6,3);
t351 = -m(6) * t279 + t322 * mrSges(6,1) + t332 * t310;
t276 = t299 * mrSges(6,2) + t324 * t303 - t351;
t282 = t338 * t284 + t369 * t287;
t278 = -t331 * pkin(4) + t322 * qJ(5) + 0.2e1 * qJD(5) * t332 - t323 * t302 + t282;
t298 = t324 * qJD(4) - qJDD(3) * t369 + t338 * t328;
t309 = -t332 * mrSges(6,1) + t324 * mrSges(6,2);
t356 = m(6) * t278 + t322 * mrSges(6,3) + t332 * t309;
t363 = t378 * t323 - t367 * t324 - t365 * t332;
t372 = t367 * t323 - t379 * t324 + t375 * t332;
t374 = -Ifges(5,3) - Ifges(6,2);
t377 = -t375 * t299 - t372 * t323 - t365 * t298 - t374 * t322 + mrSges(5,1) * t281 - mrSges(6,1) * t279 - mrSges(5,2) * t282 + mrSges(6,3) * t278 - pkin(4) * t276 + qJ(5) * (-t298 * mrSges(6,2) - t323 * t303 + t356) - t363 * t324;
t368 = -mrSges(5,3) - mrSges(6,2);
t308 = t332 * mrSges(5,1) - t324 * mrSges(5,3);
t361 = -t323 * mrSges(5,1) - t324 * mrSges(5,2) - t303;
t271 = m(5) * t282 - t322 * mrSges(5,2) + t298 * t368 - t332 * t308 + t323 * t361 + t356;
t307 = -t332 * mrSges(5,2) - t323 * mrSges(5,3);
t273 = m(5) * t281 + t322 * mrSges(5,1) + t299 * t368 + t332 * t307 + t324 * t361 + t351;
t267 = t338 * t271 + t369 * t273;
t364 = t365 * t323 + t375 * t324 + t374 * t332;
t352 = t369 * t271 - t338 * t273;
t305 = t339 * g(3) + t341 * t313;
t325 = (t339 * mrSges(4,1) + t341 * mrSges(4,2)) * qJD(1);
t329 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t359;
t330 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t360;
t286 = -qJDD(3) * pkin(3) - t343 * pkin(7) + t326 * t360 - t305;
t280 = -0.2e1 * qJD(5) * t324 + (t323 * t332 - t299) * qJ(5) + (t324 * t332 + t298) * pkin(4) + t286;
t274 = m(6) * t280 + t298 * mrSges(6,1) - t299 * mrSges(6,3) - t324 * t309 + t323 * t310;
t345 = -m(5) * t286 - t298 * mrSges(5,1) - t299 * mrSges(5,2) - t323 * t307 - t324 * t308 - t274;
t349 = t339 * (m(4) * t306 - qJDD(3) * mrSges(4,2) + t327 * mrSges(4,3) - qJD(3) * t330 - t325 * t359 + t352) + t341 * (m(4) * t305 + qJDD(3) * mrSges(4,1) - t328 * mrSges(4,3) + qJD(3) * t329 - t325 * t360 + t345);
t319 = Ifges(4,5) * qJD(3) + (t341 * Ifges(4,1) - t339 * Ifges(4,4)) * qJD(1);
t318 = Ifges(4,6) * qJD(3) + (t341 * Ifges(4,4) - t339 * Ifges(4,2)) * qJD(1);
t315 = -qJDD(1) * pkin(1) + t347;
t314 = t344 * pkin(1) + t371;
t265 = mrSges(5,2) * t286 + mrSges(6,2) * t279 - mrSges(5,3) * t281 - mrSges(6,3) * t280 - qJ(5) * t274 - t367 * t298 + t379 * t299 - t375 * t322 + t364 * t323 + t363 * t332;
t264 = -mrSges(5,1) * t286 - mrSges(6,1) * t280 + mrSges(6,2) * t278 + mrSges(5,3) * t282 - pkin(4) * t274 - t378 * t298 + t367 * t299 + t365 * t322 + t364 * t324 - t372 * t332;
t263 = m(3) * t315 + qJDD(1) * mrSges(3,2) - t344 * mrSges(3,3) + t349;
t1 = [mrSges(2,1) * t353 - mrSges(2,2) * t350 + mrSges(3,2) * t315 - mrSges(3,3) * t314 + t341 * (mrSges(4,2) * t312 - mrSges(4,3) * t305 + Ifges(4,1) * t328 + Ifges(4,4) * t327 + Ifges(4,5) * qJDD(3) - pkin(7) * t267 - qJD(3) * t318 - t338 * t264 + t265 * t369) - t339 * (-mrSges(4,1) * t312 + mrSges(4,3) * t306 + Ifges(4,4) * t328 + Ifges(4,2) * t327 + Ifges(4,6) * qJDD(3) - pkin(3) * t267 + qJD(3) * t319 - t377) - pkin(6) * t349 - pkin(1) * t263 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t314 + m(4) * t312 - t327 * mrSges(4,1) + t344 * mrSges(3,2) + t328 * mrSges(4,2) + t267 + qJDD(1) * mrSges(3,3) + (t329 * t339 + t330 * t341) * qJD(1)) * qJ(2); t263; Ifges(4,5) * t328 + Ifges(4,6) * t327 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t305 - mrSges(4,2) * t306 + t338 * t265 + t369 * t264 + pkin(3) * t345 + pkin(7) * t352 + (t341 * t318 + t339 * t319) * qJD(1); t377; t276;];
tauJ = t1;
