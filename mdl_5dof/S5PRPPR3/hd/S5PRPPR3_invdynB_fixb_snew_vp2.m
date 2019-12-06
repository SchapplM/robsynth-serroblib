% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:25
% EndTime: 2019-12-05 15:26:27
% DurationCPUTime: 0.74s
% Computational Cost: add. (7191->167), mult. (11162->204), div. (0->0), fcn. (5968->8), ass. (0->72)
t343 = -pkin(3) - pkin(6);
t342 = mrSges(4,1) - mrSges(5,2);
t341 = -Ifges(5,4) + Ifges(4,5);
t340 = Ifges(5,5) - Ifges(4,6);
t318 = sin(pkin(7));
t320 = cos(pkin(7));
t307 = -t320 * g(1) - t318 * g(2);
t314 = -g(3) + qJDD(1);
t322 = sin(qJ(2));
t324 = cos(qJ(2));
t292 = -t322 * t307 + t324 * t314;
t290 = qJDD(2) * pkin(2) + t292;
t293 = t324 * t307 + t322 * t314;
t325 = qJD(2) ^ 2;
t291 = -t325 * pkin(2) + t293;
t317 = sin(pkin(8));
t319 = cos(pkin(8));
t285 = t319 * t290 - t317 * t291;
t329 = -t325 * qJ(4) + qJDD(4) - t285;
t282 = t343 * qJDD(2) + t329;
t306 = t318 * g(1) - t320 * g(2);
t305 = qJDD(3) - t306;
t321 = sin(qJ(5));
t323 = cos(qJ(5));
t278 = t323 * t282 - t321 * t305;
t302 = (mrSges(6,1) * t321 + mrSges(6,2) * t323) * qJD(2);
t335 = qJD(2) * qJD(5);
t304 = t323 * qJDD(2) - t321 * t335;
t337 = qJD(2) * t321;
t308 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t337;
t336 = qJD(2) * t323;
t276 = m(6) * t278 + qJDD(5) * mrSges(6,1) - t304 * mrSges(6,3) + qJD(5) * t308 - t302 * t336;
t279 = t321 * t282 + t323 * t305;
t303 = -t321 * qJDD(2) - t323 * t335;
t309 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t336;
t277 = m(6) * t279 - qJDD(5) * mrSges(6,2) + t303 * mrSges(6,3) - qJD(5) * t309 - t302 * t337;
t268 = t323 * t276 + t321 * t277;
t284 = -qJDD(2) * pkin(3) + t329;
t327 = -m(5) * t284 + t325 * mrSges(5,3) - t268;
t264 = m(4) * t285 - t325 * mrSges(4,2) + t342 * qJDD(2) + t327;
t286 = t317 * t290 + t319 * t291;
t328 = qJDD(2) * qJ(4) + 0.2e1 * qJD(4) * qJD(2) + t286;
t283 = t325 * pkin(3) - t328;
t281 = t343 * t325 + t328;
t330 = -m(6) * t281 + t303 * mrSges(6,1) - t304 * mrSges(6,2) - t308 * t337 - t309 * t336;
t326 = -m(5) * t283 + t325 * mrSges(5,2) + qJDD(2) * mrSges(5,3) - t330;
t271 = m(4) * t286 - t325 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t326;
t262 = t319 * t264 + t317 * t271;
t260 = m(3) * t292 + qJDD(2) * mrSges(3,1) - t325 * mrSges(3,2) + t262;
t332 = -t317 * t264 + t319 * t271;
t261 = m(3) * t293 - t325 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t332;
t333 = -t322 * t260 + t324 * t261;
t253 = m(2) * t307 + t333;
t338 = -t321 * t276 + t323 * t277;
t267 = m(5) * t305 + t338;
t331 = m(4) * t305 + t267;
t266 = (m(2) + m(3)) * t306 - t331;
t339 = t318 * t253 + t320 * t266;
t254 = t324 * t260 + t322 * t261;
t334 = t320 * t253 - t318 * t266;
t296 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t323 - Ifges(6,4) * t321) * qJD(2);
t295 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t323 - Ifges(6,2) * t321) * qJD(2);
t294 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t323 - Ifges(6,6) * t321) * qJD(2);
t273 = mrSges(6,2) * t281 - mrSges(6,3) * t278 + Ifges(6,1) * t304 + Ifges(6,4) * t303 + Ifges(6,5) * qJDD(5) - qJD(5) * t295 - t294 * t337;
t272 = -mrSges(6,1) * t281 + mrSges(6,3) * t279 + Ifges(6,4) * t304 + Ifges(6,2) * t303 + Ifges(6,6) * qJDD(5) + qJD(5) * t296 - t294 * t336;
t256 = mrSges(5,1) * t284 + mrSges(6,1) * t278 - mrSges(6,2) * t279 - mrSges(4,3) * t285 + Ifges(6,5) * t304 + Ifges(6,6) * t303 + Ifges(6,3) * qJDD(5) + pkin(4) * t268 - qJ(4) * t267 + t340 * t325 + (mrSges(4,2) - mrSges(5,3)) * t305 + t341 * qJDD(2) + (t323 * t295 + t321 * t296) * qJD(2);
t255 = -mrSges(5,1) * t283 + mrSges(4,3) * t286 - pkin(3) * t267 - pkin(4) * t330 - pkin(6) * t338 - t340 * qJDD(2) - t323 * t272 - t321 * t273 - t342 * t305 + t341 * t325;
t250 = -mrSges(3,2) * t306 - mrSges(3,3) * t292 + Ifges(3,5) * qJDD(2) - t325 * Ifges(3,6) - qJ(3) * t262 - t317 * t255 + t319 * t256;
t249 = mrSges(3,1) * t306 + mrSges(3,3) * t293 + t325 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t331 + qJ(3) * t332 + t319 * t255 + t317 * t256;
t248 = -pkin(1) * t254 - pkin(2) * t262 - mrSges(3,1) * t292 + mrSges(3,2) * t293 - pkin(3) * t327 - qJ(4) * t326 - mrSges(4,1) * t285 + mrSges(4,2) * t286 - mrSges(5,2) * t284 + mrSges(5,3) * t283 - t323 * t273 + t321 * t272 + pkin(6) * t268 + mrSges(2,3) * t307 - mrSges(2,1) * t314 + (pkin(3) * mrSges(5,2) - Ifges(5,1) - Ifges(3,3) - Ifges(4,3)) * qJDD(2);
t247 = mrSges(2,2) * t314 - mrSges(2,3) * t306 - pkin(5) * t254 - t322 * t249 + t324 * t250;
t1 = [-m(1) * g(1) + t334; -m(1) * g(2) + t339; -m(1) * g(3) + m(2) * t314 + t254; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t339 + t320 * t247 - t318 * t248; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t334 + t318 * t247 + t320 * t248; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t306 - mrSges(2,2) * t307 + t322 * t250 + t324 * t249 + pkin(1) * (m(3) * t306 - t331) + pkin(5) * t333;];
tauB = t1;
