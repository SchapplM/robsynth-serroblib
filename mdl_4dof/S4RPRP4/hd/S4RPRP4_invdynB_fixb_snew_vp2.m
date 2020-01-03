% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRP4
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP4_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP4_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:41
% EndTime: 2019-12-31 16:43:42
% DurationCPUTime: 0.74s
% Computational Cost: add. (4258->172), mult. (7890->212), div. (0->0), fcn. (3528->6), ass. (0->69)
t365 = Ifges(4,1) + Ifges(5,1);
t361 = Ifges(4,4) - Ifges(5,5);
t360 = Ifges(5,4) + Ifges(4,5);
t364 = Ifges(4,2) + Ifges(5,3);
t359 = Ifges(4,6) - Ifges(5,6);
t363 = Ifges(4,3) + Ifges(5,2);
t362 = mrSges(4,3) + mrSges(5,2);
t335 = -g(3) + qJDD(2);
t340 = cos(qJ(3));
t358 = t340 * t335;
t339 = sin(qJ(1));
t341 = cos(qJ(1));
t327 = t339 * g(1) - t341 * g(2);
t314 = qJDD(1) * pkin(1) + t327;
t328 = -t341 * g(1) - t339 * g(2);
t343 = qJD(1) ^ 2;
t318 = -t343 * pkin(1) + t328;
t336 = sin(pkin(6));
t337 = cos(pkin(6));
t298 = t336 * t314 + t337 * t318;
t296 = -t343 * pkin(2) + qJDD(1) * pkin(5) + t298;
t338 = sin(qJ(3));
t293 = t340 * t296 + t338 * t335;
t317 = (-mrSges(4,1) * t340 + mrSges(4,2) * t338) * qJD(1);
t351 = qJD(1) * qJD(3);
t320 = t340 * qJDD(1) - t338 * t351;
t353 = qJD(1) * t338;
t323 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t353;
t315 = (-pkin(3) * t340 - qJ(4) * t338) * qJD(1);
t342 = qJD(3) ^ 2;
t352 = qJD(1) * t340;
t290 = -t342 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t315 * t352 + t293;
t316 = (-mrSges(5,1) * t340 - mrSges(5,3) * t338) * qJD(1);
t324 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t353;
t346 = m(5) * t290 + qJDD(3) * mrSges(5,3) + qJD(3) * t324 + t316 * t352;
t285 = m(4) * t293 - qJDD(3) * mrSges(4,2) - qJD(3) * t323 + t317 * t352 + t362 * t320 + t346;
t292 = -t338 * t296 + t358;
t319 = t338 * qJDD(1) + t340 * t351;
t325 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t352;
t291 = -qJDD(3) * pkin(3) - t342 * qJ(4) - t358 + qJDD(4) + (qJD(1) * t315 + t296) * t338;
t326 = mrSges(5,2) * t352 + qJD(3) * mrSges(5,3);
t345 = -m(5) * t291 + qJDD(3) * mrSges(5,1) + qJD(3) * t326;
t286 = m(4) * t292 + qJDD(3) * mrSges(4,1) + qJD(3) * t325 - t362 * t319 + (-t316 - t317) * t353 + t345;
t347 = t340 * t285 - t338 * t286;
t278 = m(3) * t298 - t343 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t347;
t297 = t337 * t314 - t336 * t318;
t295 = -qJDD(1) * pkin(2) - t343 * pkin(5) - t297;
t288 = -t320 * pkin(3) - t319 * qJ(4) + (-0.2e1 * qJD(4) * t338 + (pkin(3) * t338 - qJ(4) * t340) * qJD(3)) * qJD(1) + t295;
t287 = m(5) * t288 - t320 * mrSges(5,1) - t319 * mrSges(5,3) - t324 * t353 - t326 * t352;
t344 = -m(4) * t295 + t320 * mrSges(4,1) - t319 * mrSges(4,2) - t323 * t353 + t325 * t352 - t287;
t281 = m(3) * t297 + qJDD(1) * mrSges(3,1) - t343 * mrSges(3,2) + t344;
t273 = t336 * t278 + t337 * t281;
t271 = m(2) * t327 + qJDD(1) * mrSges(2,1) - t343 * mrSges(2,2) + t273;
t348 = t337 * t278 - t336 * t281;
t272 = m(2) * t328 - t343 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t348;
t357 = t341 * t271 + t339 * t272;
t279 = t338 * t285 + t340 * t286;
t356 = -t359 * qJD(3) + (-t361 * t338 - t364 * t340) * qJD(1);
t355 = t363 * qJD(3) + (t360 * t338 + t359 * t340) * qJD(1);
t354 = t360 * qJD(3) + (t365 * t338 + t361 * t340) * qJD(1);
t350 = m(3) * t335 + t279;
t349 = -t339 * t271 + t341 * t272;
t275 = mrSges(4,2) * t295 + mrSges(5,2) * t291 - mrSges(4,3) * t292 - mrSges(5,3) * t288 - qJ(4) * t287 + t356 * qJD(3) + t360 * qJDD(3) + t365 * t319 + t361 * t320 + t355 * t352;
t274 = -mrSges(4,1) * t295 - mrSges(5,1) * t288 + mrSges(5,2) * t290 + mrSges(4,3) * t293 - pkin(3) * t287 + t354 * qJD(3) + t359 * qJDD(3) + t361 * t319 + t364 * t320 - t355 * t353;
t267 = Ifges(3,6) * qJDD(1) + t343 * Ifges(3,5) - mrSges(3,1) * t335 + mrSges(3,3) * t298 - mrSges(4,1) * t292 + mrSges(4,2) * t293 + mrSges(5,1) * t291 - mrSges(5,3) * t290 - pkin(3) * t345 - qJ(4) * t346 - pkin(2) * t279 + (-qJ(4) * mrSges(5,2) - t359) * t320 + (pkin(3) * mrSges(5,2) - t360) * t319 - t363 * qJDD(3) + (t354 * t340 + (pkin(3) * t316 + t356) * t338) * qJD(1);
t266 = mrSges(3,2) * t335 - mrSges(3,3) * t297 + Ifges(3,5) * qJDD(1) - t343 * Ifges(3,6) - pkin(5) * t279 - t338 * t274 + t340 * t275;
t265 = -mrSges(2,2) * g(3) - mrSges(2,3) * t327 + Ifges(2,5) * qJDD(1) - t343 * Ifges(2,6) - qJ(2) * t273 + t337 * t266 - t336 * t267;
t264 = mrSges(2,1) * g(3) + mrSges(2,3) * t328 + t343 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t350 + qJ(2) * t348 + t336 * t266 + t337 * t267;
t1 = [-m(1) * g(1) + t349; -m(1) * g(2) + t357; (-m(1) - m(2)) * g(3) + t350; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t357 - t339 * t264 + t341 * t265; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t349 + t341 * t264 + t339 * t265; pkin(1) * t273 + mrSges(2,1) * t327 - mrSges(2,2) * t328 + mrSges(3,1) * t297 - mrSges(3,2) * t298 + t338 * t275 + t340 * t274 + pkin(2) * t344 + pkin(5) * t347 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
