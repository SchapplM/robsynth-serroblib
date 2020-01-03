% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRR8
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR8_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR8_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:10
% EndTime: 2019-12-31 16:55:10
% DurationCPUTime: 0.76s
% Computational Cost: add. (5321->191), mult. (10243->237), div. (0->0), fcn. (5498->6), ass. (0->76)
t347 = sin(qJ(1));
t350 = cos(qJ(1));
t334 = -t350 * g(1) - t347 * g(2);
t357 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t334;
t371 = -pkin(1) - pkin(5);
t370 = mrSges(2,1) - mrSges(3,2);
t369 = Ifges(2,5) - Ifges(3,4);
t368 = (-Ifges(2,6) + Ifges(3,5));
t333 = t347 * g(1) - t350 * g(2);
t351 = qJD(1) ^ 2;
t356 = -t351 * qJ(2) + qJDD(2) - t333;
t316 = t371 * qJDD(1) + t356;
t346 = sin(qJ(3));
t349 = cos(qJ(3));
t308 = t346 * g(3) + t349 * t316;
t364 = qJD(1) * qJD(3);
t362 = t346 * t364;
t329 = t349 * qJDD(1) - t362;
t296 = (-t329 - t362) * pkin(6) + (-t346 * t349 * t351 + qJDD(3)) * pkin(3) + t308;
t309 = -t349 * g(3) + t346 * t316;
t328 = -t346 * qJDD(1) - t349 * t364;
t365 = qJD(1) * t349;
t332 = qJD(3) * pkin(3) - pkin(6) * t365;
t342 = t346 ^ 2;
t297 = -t342 * t351 * pkin(3) + t328 * pkin(6) - qJD(3) * t332 + t309;
t345 = sin(qJ(4));
t348 = cos(qJ(4));
t294 = t348 * t296 - t345 * t297;
t322 = (-t345 * t349 - t346 * t348) * qJD(1);
t302 = t322 * qJD(4) + t345 * t328 + t348 * t329;
t323 = (-t345 * t346 + t348 * t349) * qJD(1);
t307 = -t322 * mrSges(5,1) + t323 * mrSges(5,2);
t339 = qJD(3) + qJD(4);
t313 = -t339 * mrSges(5,2) + t322 * mrSges(5,3);
t338 = qJDD(3) + qJDD(4);
t292 = m(5) * t294 + t338 * mrSges(5,1) - t302 * mrSges(5,3) - t323 * t307 + t339 * t313;
t295 = t345 * t296 + t348 * t297;
t301 = -t323 * qJD(4) + t348 * t328 - t345 * t329;
t314 = t339 * mrSges(5,1) - t323 * mrSges(5,3);
t293 = m(5) * t295 - t338 * mrSges(5,2) + t301 * mrSges(5,3) + t322 * t307 - t339 * t314;
t283 = t348 * t292 + t345 * t293;
t327 = (mrSges(4,1) * t346 + mrSges(4,2) * t349) * qJD(1);
t366 = qJD(1) * t346;
t330 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t366;
t281 = m(4) * t308 + qJDD(3) * mrSges(4,1) - t329 * mrSges(4,3) + qJD(3) * t330 - t327 * t365 + t283;
t331 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t365;
t359 = -t345 * t292 + t348 * t293;
t282 = m(4) * t309 - qJDD(3) * mrSges(4,2) + t328 * mrSges(4,3) - qJD(3) * t331 - t327 * t366 + t359;
t279 = t349 * t281 + t346 * t282;
t318 = -qJDD(1) * pkin(1) + t356;
t354 = -m(3) * t318 + (t351 * mrSges(3,3)) - t279;
t277 = m(2) * t333 - (t351 * mrSges(2,2)) + t370 * qJDD(1) + t354;
t317 = t351 * pkin(1) - t357;
t315 = t371 * t351 + t357;
t299 = t332 * t365 - t328 * pkin(3) + (-pkin(6) * t342 + t371) * t351 + t357;
t355 = m(5) * t299 - t301 * mrSges(5,1) + t302 * mrSges(5,2) - t322 * t313 + t323 * t314;
t353 = -m(4) * t315 + t328 * mrSges(4,1) - t329 * mrSges(4,2) - t330 * t366 - t331 * t365 - t355;
t352 = -m(3) * t317 + (t351 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t353;
t288 = m(2) * t334 - (t351 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t352;
t367 = t350 * t277 + t347 * t288;
t361 = -t347 * t277 + t350 * t288;
t360 = -t346 * t281 + t349 * t282;
t321 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t349 - Ifges(4,4) * t346) * qJD(1);
t320 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t349 - Ifges(4,2) * t346) * qJD(1);
t319 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t349 - Ifges(4,6) * t346) * qJD(1);
t305 = Ifges(5,1) * t323 + Ifges(5,4) * t322 + Ifges(5,5) * t339;
t304 = Ifges(5,4) * t323 + Ifges(5,2) * t322 + Ifges(5,6) * t339;
t303 = Ifges(5,5) * t323 + Ifges(5,6) * t322 + Ifges(5,3) * t339;
t285 = mrSges(5,2) * t299 - mrSges(5,3) * t294 + Ifges(5,1) * t302 + Ifges(5,4) * t301 + Ifges(5,5) * t338 + t322 * t303 - t339 * t304;
t284 = -mrSges(5,1) * t299 + mrSges(5,3) * t295 + Ifges(5,4) * t302 + Ifges(5,2) * t301 + Ifges(5,6) * t338 - t323 * t303 + t339 * t305;
t278 = -m(3) * g(3) + t360;
t275 = mrSges(4,2) * t315 - mrSges(4,3) * t308 + Ifges(4,1) * t329 + Ifges(4,4) * t328 + Ifges(4,5) * qJDD(3) - pkin(6) * t283 - qJD(3) * t320 - t345 * t284 + t348 * t285 - t319 * t366;
t274 = -mrSges(4,1) * t315 + mrSges(4,3) * t309 + Ifges(4,4) * t329 + Ifges(4,2) * t328 + Ifges(4,6) * qJDD(3) - pkin(3) * t355 + pkin(6) * t359 + qJD(3) * t321 + t348 * t284 + t345 * t285 - t319 * t365;
t273 = Ifges(4,3) * qJDD(3) + Ifges(5,3) * t338 + t323 * t304 + Ifges(4,6) * t328 + Ifges(4,5) * t329 - mrSges(2,3) * t333 - t322 * t305 - mrSges(4,2) * t309 + mrSges(3,1) * t318 + Ifges(5,6) * t301 + Ifges(5,5) * t302 + mrSges(4,1) * t308 + mrSges(5,1) * t294 - mrSges(5,2) * t295 + pkin(3) * t283 + pkin(2) * t279 - qJ(2) * t278 + (t368 * t351) + t369 * qJDD(1) + (t349 * t320 + t346 * t321) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t272 = -mrSges(3,1) * t317 + mrSges(2,3) * t334 - pkin(1) * t278 - pkin(2) * t353 - pkin(5) * t360 + t370 * g(3) - t368 * qJDD(1) - t349 * t274 - t346 * t275 + t369 * t351;
t1 = [-m(1) * g(1) + t361; -m(1) * g(2) + t367; (-m(1) - m(2) - m(3)) * g(3) + t360; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t367 - t347 * t272 + t350 * t273; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t361 + t350 * t272 + t347 * t273; pkin(1) * t354 + qJ(2) * t352 + t349 * t275 - t346 * t274 - pkin(5) * t279 + mrSges(2,1) * t333 - mrSges(2,2) * t334 + mrSges(3,2) * t318 - mrSges(3,3) * t317 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
