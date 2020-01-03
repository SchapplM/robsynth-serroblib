% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRR3
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:25
% EndTime: 2019-12-31 17:24:27
% DurationCPUTime: 0.88s
% Computational Cost: add. (6391->202), mult. (13865->265), div. (0->0), fcn. (9136->8), ass. (0->83)
t316 = qJD(1) ^ 2;
t330 = pkin(2) * t316;
t311 = sin(qJ(1));
t315 = cos(qJ(1));
t322 = -g(1) * t315 - g(2) * t311;
t292 = -pkin(1) * t316 + qJDD(1) * pkin(5) + t322;
t310 = sin(qJ(2));
t329 = t292 * t310;
t314 = cos(qJ(2));
t326 = qJD(1) * qJD(2);
t295 = qJDD(1) * t310 + t314 * t326;
t264 = qJDD(2) * pkin(2) - pkin(6) * t295 - t329 + (pkin(6) * t326 + t310 * t330 - g(3)) * t314;
t280 = -g(3) * t310 + t314 * t292;
t296 = qJDD(1) * t314 - t310 * t326;
t328 = qJD(1) * t310;
t299 = qJD(2) * pkin(2) - pkin(6) * t328;
t307 = t314 ^ 2;
t265 = pkin(6) * t296 - qJD(2) * t299 - t307 * t330 + t280;
t309 = sin(qJ(3));
t313 = cos(qJ(3));
t253 = t313 * t264 - t265 * t309;
t289 = (-t309 * t310 + t313 * t314) * qJD(1);
t269 = qJD(3) * t289 + t295 * t313 + t296 * t309;
t290 = (t309 * t314 + t310 * t313) * qJD(1);
t305 = qJDD(2) + qJDD(3);
t306 = qJD(2) + qJD(3);
t242 = (t289 * t306 - t269) * pkin(7) + (t289 * t290 + t305) * pkin(3) + t253;
t254 = t309 * t264 + t313 * t265;
t268 = -qJD(3) * t290 - t295 * t309 + t296 * t313;
t283 = pkin(3) * t306 - pkin(7) * t290;
t285 = t289 ^ 2;
t243 = -pkin(3) * t285 + pkin(7) * t268 - t283 * t306 + t254;
t308 = sin(qJ(4));
t312 = cos(qJ(4));
t240 = t242 * t312 - t243 * t308;
t276 = t289 * t312 - t290 * t308;
t250 = qJD(4) * t276 + t268 * t308 + t269 * t312;
t277 = t289 * t308 + t290 * t312;
t259 = -mrSges(5,1) * t276 + mrSges(5,2) * t277;
t303 = qJD(4) + t306;
t271 = -mrSges(5,2) * t303 + mrSges(5,3) * t276;
t302 = qJDD(4) + t305;
t237 = m(5) * t240 + mrSges(5,1) * t302 - mrSges(5,3) * t250 - t259 * t277 + t271 * t303;
t241 = t242 * t308 + t243 * t312;
t249 = -qJD(4) * t277 + t268 * t312 - t269 * t308;
t272 = mrSges(5,1) * t303 - mrSges(5,3) * t277;
t238 = m(5) * t241 - mrSges(5,2) * t302 + mrSges(5,3) * t249 + t259 * t276 - t272 * t303;
t231 = t312 * t237 + t308 * t238;
t278 = -mrSges(4,1) * t289 + mrSges(4,2) * t290;
t281 = -mrSges(4,2) * t306 + mrSges(4,3) * t289;
t228 = m(4) * t253 + mrSges(4,1) * t305 - mrSges(4,3) * t269 - t278 * t290 + t281 * t306 + t231;
t282 = mrSges(4,1) * t306 - mrSges(4,3) * t290;
t323 = -t237 * t308 + t312 * t238;
t229 = m(4) * t254 - mrSges(4,2) * t305 + mrSges(4,3) * t268 + t278 * t289 - t282 * t306 + t323;
t224 = t313 * t228 + t309 * t229;
t327 = qJD(1) * t314;
t325 = g(1) * t311 - t315 * g(2);
t324 = -t228 * t309 + t313 * t229;
t321 = -qJDD(1) * pkin(1) - t325;
t270 = -pkin(2) * t296 + t299 * t328 + (-pkin(6) * t307 - pkin(5)) * t316 + t321;
t245 = -pkin(3) * t268 - pkin(7) * t285 + t283 * t290 + t270;
t320 = -m(5) * t245 + mrSges(5,1) * t249 - t250 * mrSges(5,2) + t276 * t271 - t277 * t272;
t256 = Ifges(5,4) * t277 + Ifges(5,2) * t276 + Ifges(5,6) * t303;
t257 = Ifges(5,1) * t277 + Ifges(5,4) * t276 + Ifges(5,5) * t303;
t319 = mrSges(5,1) * t240 - mrSges(5,2) * t241 + Ifges(5,5) * t250 + Ifges(5,6) * t249 + Ifges(5,3) * t302 + t277 * t256 - t276 * t257;
t318 = m(4) * t270 - mrSges(4,1) * t268 + mrSges(4,2) * t269 - t281 * t289 + t282 * t290 - t320;
t274 = Ifges(4,4) * t290 + Ifges(4,2) * t289 + Ifges(4,6) * t306;
t275 = Ifges(4,1) * t290 + Ifges(4,4) * t289 + Ifges(4,5) * t306;
t317 = mrSges(4,1) * t253 - mrSges(4,2) * t254 + Ifges(4,5) * t269 + Ifges(4,6) * t268 + Ifges(4,3) * t305 + pkin(3) * t231 + t290 * t274 - t289 * t275 + t319;
t298 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t327;
t297 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t328;
t294 = (-mrSges(3,1) * t314 + mrSges(3,2) * t310) * qJD(1);
t291 = -pkin(5) * t316 + t321;
t288 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t310 + Ifges(3,4) * t314) * qJD(1);
t287 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t310 + Ifges(3,2) * t314) * qJD(1);
t279 = -g(3) * t314 - t329;
t273 = Ifges(4,5) * t290 + Ifges(4,6) * t289 + Ifges(4,3) * t306;
t255 = Ifges(5,5) * t277 + Ifges(5,6) * t276 + Ifges(5,3) * t303;
t233 = mrSges(5,2) * t245 - mrSges(5,3) * t240 + Ifges(5,1) * t250 + Ifges(5,4) * t249 + Ifges(5,5) * t302 + t255 * t276 - t256 * t303;
t232 = -mrSges(5,1) * t245 + mrSges(5,3) * t241 + Ifges(5,4) * t250 + Ifges(5,2) * t249 + Ifges(5,6) * t302 - t255 * t277 + t257 * t303;
t223 = mrSges(4,2) * t270 - mrSges(4,3) * t253 + Ifges(4,1) * t269 + Ifges(4,4) * t268 + Ifges(4,5) * t305 - pkin(7) * t231 - t232 * t308 + t233 * t312 + t273 * t289 - t274 * t306;
t222 = -mrSges(4,1) * t270 + mrSges(4,3) * t254 + Ifges(4,4) * t269 + Ifges(4,2) * t268 + Ifges(4,6) * t305 + pkin(3) * t320 + pkin(7) * t323 + t312 * t232 + t308 * t233 - t290 * t273 + t306 * t275;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t325 - mrSges(2,2) * t322 + t310 * (mrSges(3,2) * t291 - mrSges(3,3) * t279 + Ifges(3,1) * t295 + Ifges(3,4) * t296 + Ifges(3,5) * qJDD(2) - pkin(6) * t224 - qJD(2) * t287 - t309 * t222 + t313 * t223) + t314 * (-mrSges(3,1) * t291 + mrSges(3,3) * t280 + Ifges(3,4) * t295 + Ifges(3,2) * t296 + Ifges(3,6) * qJDD(2) - pkin(2) * t318 + pkin(6) * t324 + qJD(2) * t288 + t313 * t222 + t309 * t223) + pkin(1) * (-t318 - m(3) * t291 + mrSges(3,1) * t296 - mrSges(3,2) * t295 + (-t297 * t310 + t298 * t314) * qJD(1)) + pkin(5) * (t314 * (m(3) * t280 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t296 - qJD(2) * t297 + t294 * t327 + t324) - t310 * (m(3) * t279 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t295 + qJD(2) * t298 - t294 * t328 + t224)); t317 + pkin(2) * t224 + Ifges(3,6) * t296 + Ifges(3,5) * t295 + mrSges(3,1) * t279 - mrSges(3,2) * t280 + Ifges(3,3) * qJDD(2) + (t310 * t287 - t314 * t288) * qJD(1); t317; t319;];
tauJ = t1;
