% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP3_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP3_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:05
% EndTime: 2019-12-31 17:14:06
% DurationCPUTime: 0.76s
% Computational Cost: add. (6254->174), mult. (7890->214), div. (0->0), fcn. (3528->6), ass. (0->70)
t371 = Ifges(4,1) + Ifges(5,1);
t366 = Ifges(4,4) - Ifges(5,5);
t365 = Ifges(4,5) + Ifges(5,4);
t370 = Ifges(4,2) + Ifges(5,3);
t364 = Ifges(4,6) - Ifges(5,6);
t369 = Ifges(4,3) + Ifges(5,2);
t346 = cos(qJ(3));
t368 = t346 * g(3);
t367 = mrSges(4,3) + mrSges(5,2);
t340 = qJD(1) + qJD(2);
t343 = sin(qJ(3));
t363 = t340 * t343;
t362 = t340 * t346;
t345 = sin(qJ(1));
t348 = cos(qJ(1));
t334 = t345 * g(1) - t348 * g(2);
t328 = qJDD(1) * pkin(1) + t334;
t335 = -t348 * g(1) - t345 * g(2);
t350 = qJD(1) ^ 2;
t329 = -t350 * pkin(1) + t335;
t344 = sin(qJ(2));
t347 = cos(qJ(2));
t305 = t344 * t328 + t347 * t329;
t338 = t340 ^ 2;
t339 = qJDD(1) + qJDD(2);
t303 = -t338 * pkin(2) + t339 * pkin(6) + t305;
t300 = -t343 * g(3) + t346 * t303;
t320 = (-mrSges(4,1) * t346 + mrSges(4,2) * t343) * t340;
t357 = qJD(3) * t340;
t322 = t346 * t339 - t343 * t357;
t330 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t363;
t318 = (-pkin(3) * t346 - qJ(4) * t343) * t340;
t349 = qJD(3) ^ 2;
t297 = -t349 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t318 * t362 + t300;
t319 = (-mrSges(5,1) * t346 - mrSges(5,3) * t343) * t340;
t331 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t363;
t353 = m(5) * t297 + qJDD(3) * mrSges(5,3) + qJD(3) * t331 + t319 * t362;
t292 = m(4) * t300 - qJDD(3) * mrSges(4,2) - qJD(3) * t330 + t320 * t362 + t367 * t322 + t353;
t299 = -t343 * t303 - t368;
t321 = t343 * t339 + t346 * t357;
t332 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t362;
t298 = -qJDD(3) * pkin(3) + t368 - t349 * qJ(4) + qJDD(4) + (t318 * t340 + t303) * t343;
t333 = mrSges(5,2) * t362 + qJD(3) * mrSges(5,3);
t352 = -m(5) * t298 + qJDD(3) * mrSges(5,1) + qJD(3) * t333;
t293 = m(4) * t299 + qJDD(3) * mrSges(4,1) + qJD(3) * t332 + (-t319 - t320) * t363 - t367 * t321 + t352;
t354 = t346 * t292 - t343 * t293;
t285 = m(3) * t305 - t338 * mrSges(3,1) - t339 * mrSges(3,2) + t354;
t304 = t347 * t328 - t344 * t329;
t302 = -t339 * pkin(2) - t338 * pkin(6) - t304;
t295 = -t322 * pkin(3) - t321 * qJ(4) + (-0.2e1 * qJD(4) * t343 + (pkin(3) * t343 - qJ(4) * t346) * qJD(3)) * t340 + t302;
t294 = m(5) * t295 - t322 * mrSges(5,1) - t321 * mrSges(5,3) - t331 * t363 - t333 * t362;
t351 = -m(4) * t302 + t322 * mrSges(4,1) - t321 * mrSges(4,2) - t330 * t363 + t332 * t362 - t294;
t288 = m(3) * t304 + t339 * mrSges(3,1) - t338 * mrSges(3,2) + t351;
t280 = t344 * t285 + t347 * t288;
t278 = m(2) * t334 + qJDD(1) * mrSges(2,1) - t350 * mrSges(2,2) + t280;
t355 = t347 * t285 - t344 * t288;
t279 = m(2) * t335 - t350 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t355;
t361 = t348 * t278 + t345 * t279;
t286 = t343 * t292 + t346 * t293;
t360 = (-t366 * t343 - t370 * t346) * t340 - t364 * qJD(3);
t359 = (t365 * t343 + t364 * t346) * t340 + t369 * qJD(3);
t358 = (t371 * t343 + t366 * t346) * t340 + t365 * qJD(3);
t356 = -t345 * t278 + t348 * t279;
t282 = mrSges(4,2) * t302 + mrSges(5,2) * t298 - mrSges(4,3) * t299 - mrSges(5,3) * t295 - qJ(4) * t294 + t360 * qJD(3) + t365 * qJDD(3) + t371 * t321 + t366 * t322 + t359 * t362;
t281 = -mrSges(4,1) * t302 - mrSges(5,1) * t295 + mrSges(5,2) * t297 + mrSges(4,3) * t300 - pkin(3) * t294 + t358 * qJD(3) + t364 * qJDD(3) + t366 * t321 + t370 * t322 - t359 * t363;
t274 = Ifges(3,6) * t339 + t338 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t305 - mrSges(4,1) * t299 + mrSges(4,2) * t300 + mrSges(5,1) * t298 - mrSges(5,3) * t297 - pkin(3) * t352 - qJ(4) * t353 - pkin(2) * t286 + (-qJ(4) * mrSges(5,2) - t364) * t322 + (pkin(3) * mrSges(5,2) - t365) * t321 - t369 * qJDD(3) + (t358 * t346 + (pkin(3) * t319 + t360) * t343) * t340;
t273 = -mrSges(3,2) * g(3) - mrSges(3,3) * t304 + Ifges(3,5) * t339 - t338 * Ifges(3,6) - pkin(6) * t286 - t343 * t281 + t346 * t282;
t272 = -mrSges(2,2) * g(3) - mrSges(2,3) * t334 + Ifges(2,5) * qJDD(1) - t350 * Ifges(2,6) - pkin(5) * t280 + t347 * t273 - t344 * t274;
t271 = Ifges(2,6) * qJDD(1) + t350 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t335 + t344 * t273 + t347 * t274 - pkin(1) * (-m(3) * g(3) + t286) + pkin(5) * t355;
t1 = [-m(1) * g(1) + t356; -m(1) * g(2) + t361; (-m(1) - m(2) - m(3)) * g(3) + t286; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t361 - t345 * t271 + t348 * t272; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t356 + t348 * t271 + t345 * t272; -mrSges(1,1) * g(2) + mrSges(2,1) * t334 + mrSges(3,1) * t304 + mrSges(1,2) * g(1) - mrSges(2,2) * t335 - mrSges(3,2) * t305 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t339 + pkin(1) * t280 + pkin(2) * t351 + pkin(6) * t354 + t346 * t281 + t343 * t282;];
tauB = t1;
