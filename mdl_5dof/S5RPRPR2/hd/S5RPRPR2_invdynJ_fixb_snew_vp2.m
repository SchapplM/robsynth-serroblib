% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:56
% EndTime: 2022-01-23 09:18:57
% DurationCPUTime: 0.54s
% Computational Cost: add. (4896->137), mult. (7045->180), div. (0->0), fcn. (4082->10), ass. (0->73)
t306 = qJD(1) + qJD(3);
t302 = t306 ^ 2;
t310 = cos(pkin(9));
t338 = pkin(4) * t310;
t308 = sin(pkin(9));
t337 = mrSges(5,2) * t308;
t305 = t310 ^ 2;
t336 = t302 * t305;
t303 = qJDD(1) + qJDD(3);
t335 = t303 * t310;
t314 = sin(qJ(1));
t317 = cos(qJ(1));
t329 = t314 * g(1) - t317 * g(2);
t291 = qJDD(1) * pkin(1) + t329;
t318 = qJD(1) ^ 2;
t326 = -t317 * g(1) - t314 * g(2);
t292 = -t318 * pkin(1) + t326;
t309 = sin(pkin(8));
t311 = cos(pkin(8));
t280 = t311 * t291 - t309 * t292;
t278 = qJDD(1) * pkin(2) + t280;
t281 = t309 * t291 + t311 * t292;
t279 = -t318 * pkin(2) + t281;
t313 = sin(qJ(3));
t316 = cos(qJ(3));
t266 = t313 * t278 + t316 * t279;
t263 = -t302 * pkin(3) + t303 * qJ(4) + t266;
t307 = -g(3) + qJDD(2);
t330 = qJD(4) * t306;
t331 = t310 * t307 - 0.2e1 * t308 * t330;
t259 = -t308 * t263 + t331;
t324 = mrSges(5,3) * t303 + (-mrSges(5,1) * t310 + t337) * t302;
t256 = (-pkin(7) * t303 + t302 * t338 - t263) * t308 + t331;
t260 = t308 * t307 + (t263 + 0.2e1 * t330) * t310;
t257 = -pkin(4) * t336 + pkin(7) * t335 + t260;
t312 = sin(qJ(5));
t315 = cos(qJ(5));
t254 = t315 * t256 - t312 * t257;
t322 = -t308 * t312 + t310 * t315;
t284 = t322 * t306;
t323 = t308 * t315 + t310 * t312;
t285 = t323 * t306;
t272 = -t284 * mrSges(6,1) + t285 * mrSges(6,2);
t274 = t284 * qJD(5) + t323 * t303;
t282 = -qJD(5) * mrSges(6,2) + t284 * mrSges(6,3);
t252 = m(6) * t254 + qJDD(5) * mrSges(6,1) - t274 * mrSges(6,3) + qJD(5) * t282 - t285 * t272;
t255 = t312 * t256 + t315 * t257;
t273 = -t285 * qJD(5) + t322 * t303;
t283 = qJD(5) * mrSges(6,1) - t285 * mrSges(6,3);
t253 = m(6) * t255 - qJDD(5) * mrSges(6,2) + t273 * mrSges(6,3) - qJD(5) * t283 + t284 * t272;
t332 = t315 * t252 + t312 * t253;
t241 = m(5) * t259 - t324 * t308 + t332;
t327 = -t312 * t252 + t315 * t253;
t242 = m(5) * t260 + t324 * t310 + t327;
t328 = -t308 * t241 + t310 * t242;
t238 = m(4) * t266 - t302 * mrSges(4,1) - t303 * mrSges(4,2) + t328;
t265 = t316 * t278 - t313 * t279;
t325 = qJDD(4) - t265;
t262 = -t303 * pkin(3) - t302 * qJ(4) + t325;
t304 = t308 ^ 2;
t258 = (-pkin(3) - t338) * t303 + (-qJ(4) + (-t304 - t305) * pkin(7)) * t302 + t325;
t320 = m(6) * t258 - t273 * mrSges(6,1) + t274 * mrSges(6,2) - t284 * t282 + t285 * t283;
t319 = -m(5) * t262 + mrSges(5,1) * t335 - t320 + (t302 * t304 + t336) * mrSges(5,3);
t246 = m(4) * t265 - t302 * mrSges(4,2) + (mrSges(4,1) - t337) * t303 + t319;
t333 = t313 * t238 + t316 * t246;
t267 = Ifges(6,5) * t285 + Ifges(6,6) * t284 + Ifges(6,3) * qJD(5);
t269 = Ifges(6,1) * t285 + Ifges(6,4) * t284 + Ifges(6,5) * qJD(5);
t243 = -mrSges(6,1) * t258 + mrSges(6,3) * t255 + Ifges(6,4) * t274 + Ifges(6,2) * t273 + Ifges(6,6) * qJDD(5) + qJD(5) * t269 - t285 * t267;
t268 = Ifges(6,4) * t285 + Ifges(6,2) * t284 + Ifges(6,6) * qJD(5);
t244 = mrSges(6,2) * t258 - mrSges(6,3) * t254 + Ifges(6,1) * t274 + Ifges(6,4) * t273 + Ifges(6,5) * qJDD(5) - qJD(5) * t268 + t284 * t267;
t248 = t303 * t337 - t319;
t321 = -mrSges(4,2) * t266 + t310 * (-mrSges(5,1) * t262 + mrSges(5,3) * t260 + t312 * t244 + t315 * t243 - pkin(4) * t320 + pkin(7) * t327 + (Ifges(5,4) * t308 + Ifges(5,2) * t310) * t303) + t308 * (mrSges(5,2) * t262 - mrSges(5,3) * t259 + t315 * t244 - t312 * t243 - pkin(7) * t332 + (Ifges(5,1) * t308 + Ifges(5,4) * t310) * t303) + qJ(4) * t328 - pkin(3) * t248 + mrSges(4,1) * t265 + Ifges(4,3) * t303;
t1 = [pkin(1) * (t309 * (m(3) * t281 - t318 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t316 * t238 - t313 * t246) + t311 * (m(3) * t280 + qJDD(1) * mrSges(3,1) - t318 * mrSges(3,2) + t333)) + mrSges(2,1) * t329 - mrSges(2,2) * t326 + pkin(2) * t333 + mrSges(3,1) * t280 - mrSges(3,2) * t281 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + t321; t310 * t241 + t308 * t242 + (m(3) + m(4)) * t307; t321; t248; mrSges(6,1) * t254 - mrSges(6,2) * t255 + Ifges(6,5) * t274 + Ifges(6,6) * t273 + Ifges(6,3) * qJDD(5) + t285 * t268 - t284 * t269;];
tauJ = t1;
