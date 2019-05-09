% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 12:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:16:11
% EndTime: 2019-05-06 12:16:47
% DurationCPUTime: 18.35s
% Computational Cost: add. (299714->382), mult. (692046->471), div. (0->0), fcn. (491923->10), ass. (0->143)
t372 = -2 * qJD(3);
t330 = sin(qJ(2));
t333 = cos(qJ(2));
t360 = qJD(1) * qJD(2);
t314 = qJDD(1) * t330 + t333 * t360;
t331 = sin(qJ(1));
t334 = cos(qJ(1));
t321 = -g(1) * t334 - g(2) * t331;
t336 = qJD(1) ^ 2;
t309 = -pkin(1) * t336 + qJDD(1) * pkin(7) + t321;
t366 = t330 * t309;
t369 = pkin(2) * t336;
t265 = qJDD(2) * pkin(2) - t314 * qJ(3) - t366 + (qJ(3) * t360 + t330 * t369 - g(3)) * t333;
t295 = -g(3) * t330 + t309 * t333;
t315 = qJDD(1) * t333 - t330 * t360;
t363 = qJD(1) * t330;
t317 = qJD(2) * pkin(2) - qJ(3) * t363;
t325 = t333 ^ 2;
t266 = qJ(3) * t315 - qJD(2) * t317 - t325 * t369 + t295;
t327 = sin(pkin(9));
t328 = cos(pkin(9));
t304 = (t327 * t333 + t328 * t330) * qJD(1);
t238 = t328 * t265 - t266 * t327 + t304 * t372;
t303 = (t327 * t330 - t328 * t333) * qJD(1);
t239 = t265 * t327 + t266 * t328 + t303 * t372;
t277 = mrSges(4,1) * t303 + mrSges(4,2) * t304;
t288 = -t314 * t327 + t315 * t328;
t297 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t304;
t278 = pkin(3) * t303 - pkin(8) * t304;
t335 = qJD(2) ^ 2;
t212 = -pkin(3) * t335 + qJDD(2) * pkin(8) - t278 * t303 + t239;
t320 = t331 * g(1) - g(2) * t334;
t349 = -qJDD(1) * pkin(1) - t320;
t269 = -t315 * pkin(2) + qJDD(3) + t317 * t363 + (-qJ(3) * t325 - pkin(7)) * t336 + t349;
t289 = t314 * t328 + t315 * t327;
t216 = (qJD(2) * t303 - t289) * pkin(8) + (qJD(2) * t304 - t288) * pkin(3) + t269;
t329 = sin(qJ(4));
t332 = cos(qJ(4));
t206 = -t329 * t212 + t216 * t332;
t292 = qJD(2) * t332 - t304 * t329;
t257 = qJD(4) * t292 + qJDD(2) * t329 + t289 * t332;
t287 = qJDD(4) - t288;
t293 = qJD(2) * t329 + t304 * t332;
t302 = qJD(4) + t303;
t202 = (t292 * t302 - t257) * qJ(5) + (t292 * t293 + t287) * pkin(4) + t206;
t207 = t212 * t332 + t216 * t329;
t256 = -qJD(4) * t293 + qJDD(2) * t332 - t289 * t329;
t271 = pkin(4) * t302 - qJ(5) * t293;
t291 = t292 ^ 2;
t204 = -pkin(4) * t291 + qJ(5) * t256 - t271 * t302 + t207;
t326 = sin(pkin(10));
t367 = cos(pkin(10));
t262 = -t292 * t367 + t293 * t326;
t370 = -2 * qJD(5);
t198 = t202 * t326 + t204 * t367 + t262 * t370;
t226 = -t256 * t367 + t257 * t326;
t263 = t292 * t326 + t293 * t367;
t247 = mrSges(6,1) * t302 - mrSges(6,3) * t263;
t240 = pkin(5) * t262 - qJ(6) * t263;
t301 = t302 ^ 2;
t193 = -pkin(5) * t301 + qJ(6) * t287 + 0.2e1 * qJD(6) * t302 - t240 * t262 + t198;
t248 = -mrSges(7,1) * t302 + mrSges(7,2) * t263;
t359 = m(7) * t193 + mrSges(7,3) * t287 + t248 * t302;
t241 = mrSges(7,1) * t262 - mrSges(7,3) * t263;
t364 = -mrSges(6,1) * t262 - mrSges(6,2) * t263 - t241;
t368 = -mrSges(6,3) - mrSges(7,2);
t179 = m(6) * t198 - t287 * mrSges(6,2) + t226 * t368 - t302 * t247 + t262 * t364 + t359;
t348 = t202 * t367 - t326 * t204;
t197 = t263 * t370 + t348;
t227 = t256 * t326 + t257 * t367;
t246 = -mrSges(6,2) * t302 - mrSges(6,3) * t262;
t195 = -t287 * pkin(5) - t301 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t240) * t263 - t348;
t245 = -mrSges(7,2) * t262 + mrSges(7,3) * t302;
t353 = -m(7) * t195 + mrSges(7,1) * t287 + t245 * t302;
t181 = m(6) * t197 + t287 * mrSges(6,1) + t227 * t368 + t302 * t246 + t263 * t364 + t353;
t174 = t179 * t326 + t181 * t367;
t267 = -mrSges(5,1) * t292 + mrSges(5,2) * t293;
t270 = -mrSges(5,2) * t302 + mrSges(5,3) * t292;
t172 = m(5) * t206 + mrSges(5,1) * t287 - mrSges(5,3) * t257 - t267 * t293 + t270 * t302 + t174;
t272 = mrSges(5,1) * t302 - mrSges(5,3) * t293;
t355 = t179 * t367 - t181 * t326;
t173 = m(5) * t207 - mrSges(5,2) * t287 + mrSges(5,3) * t256 + t267 * t292 - t272 * t302 + t355;
t356 = -t172 * t329 + t173 * t332;
t165 = m(4) * t239 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t288 - qJD(2) * t297 - t277 * t303 + t356;
t296 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t303;
t211 = -qJDD(2) * pkin(3) - t335 * pkin(8) + t278 * t304 - t238;
t205 = -t256 * pkin(4) - t291 * qJ(5) + t271 * t293 + qJDD(5) + t211;
t200 = -0.2e1 * qJD(6) * t263 + (t262 * t302 - t227) * qJ(6) + (t263 * t302 + t226) * pkin(5) + t205;
t190 = m(7) * t200 + mrSges(7,1) * t226 - mrSges(7,3) * t227 + t245 * t262 - t248 * t263;
t342 = m(6) * t205 + mrSges(6,1) * t226 + mrSges(6,2) * t227 + t246 * t262 + t247 * t263 + t190;
t339 = -m(5) * t211 + mrSges(5,1) * t256 - mrSges(5,2) * t257 + t270 * t292 - t272 * t293 - t342;
t183 = m(4) * t238 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t289 + qJD(2) * t296 - t277 * t304 + t339;
t160 = t165 * t327 + t183 * t328;
t294 = -t333 * g(3) - t366;
t306 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t330 + Ifges(3,2) * t333) * qJD(1);
t307 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t330 + Ifges(3,4) * t333) * qJD(1);
t232 = Ifges(7,1) * t263 + Ifges(7,4) * t302 + Ifges(7,5) * t262;
t233 = Ifges(6,1) * t263 - Ifges(6,4) * t262 + Ifges(6,5) * t302;
t352 = -mrSges(7,1) * t200 + mrSges(7,2) * t193;
t230 = Ifges(7,4) * t263 + Ifges(7,2) * t302 + Ifges(7,6) * t262;
t365 = -Ifges(6,5) * t263 + Ifges(6,6) * t262 - Ifges(6,3) * t302 - t230;
t175 = -mrSges(6,1) * t205 + mrSges(6,3) * t198 - pkin(5) * t190 + (t232 + t233) * t302 + (Ifges(6,6) - Ifges(7,6)) * t287 + t365 * t263 + (Ifges(6,4) - Ifges(7,5)) * t227 + (-Ifges(6,2) - Ifges(7,3)) * t226 + t352;
t231 = Ifges(6,4) * t263 - Ifges(6,2) * t262 + Ifges(6,6) * t302;
t228 = Ifges(7,5) * t263 + Ifges(7,6) * t302 + Ifges(7,3) * t262;
t347 = mrSges(7,2) * t195 - mrSges(7,3) * t200 + Ifges(7,1) * t227 + Ifges(7,4) * t287 + Ifges(7,5) * t226 + t228 * t302;
t176 = mrSges(6,2) * t205 - mrSges(6,3) * t197 + Ifges(6,1) * t227 - Ifges(6,4) * t226 + Ifges(6,5) * t287 - qJ(6) * t190 - t302 * t231 + t262 * t365 + t347;
t249 = Ifges(5,5) * t293 + Ifges(5,6) * t292 + Ifges(5,3) * t302;
t251 = Ifges(5,1) * t293 + Ifges(5,4) * t292 + Ifges(5,5) * t302;
t154 = -mrSges(5,1) * t211 + mrSges(5,3) * t207 + Ifges(5,4) * t257 + Ifges(5,2) * t256 + Ifges(5,6) * t287 - pkin(4) * t342 + qJ(5) * t355 + t175 * t367 + t326 * t176 - t293 * t249 + t302 * t251;
t250 = Ifges(5,4) * t293 + Ifges(5,2) * t292 + Ifges(5,6) * t302;
t156 = mrSges(5,2) * t211 - mrSges(5,3) * t206 + Ifges(5,1) * t257 + Ifges(5,4) * t256 + Ifges(5,5) * t287 - qJ(5) * t174 - t175 * t326 + t176 * t367 + t249 * t292 - t250 * t302;
t274 = Ifges(4,4) * t304 - Ifges(4,2) * t303 + Ifges(4,6) * qJD(2);
t275 = Ifges(4,1) * t304 - Ifges(4,4) * t303 + Ifges(4,5) * qJD(2);
t343 = -mrSges(4,1) * t238 + mrSges(4,2) * t239 - Ifges(4,5) * t289 - Ifges(4,6) * t288 - Ifges(4,3) * qJDD(2) - pkin(3) * t339 - pkin(8) * t356 - t154 * t332 - t156 * t329 - t274 * t304 - t303 * t275;
t371 = mrSges(3,1) * t294 - mrSges(3,2) * t295 + Ifges(3,5) * t314 + Ifges(3,6) * t315 + Ifges(3,3) * qJDD(2) + pkin(2) * t160 + (t306 * t330 - t307 * t333) * qJD(1) - t343;
t167 = t172 * t332 + t173 * t329;
t362 = qJD(1) * t333;
t313 = (-mrSges(3,1) * t333 + mrSges(3,2) * t330) * qJD(1);
t319 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t362;
t158 = m(3) * t294 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t314 + qJD(2) * t319 - t313 * t363 + t160;
t318 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t363;
t357 = t165 * t328 - t183 * t327;
t159 = m(3) * t295 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t315 - qJD(2) * t318 + t313 * t362 + t357;
t358 = -t158 * t330 + t159 * t333;
t273 = Ifges(4,5) * t304 - Ifges(4,6) * t303 + Ifges(4,3) * qJD(2);
t148 = mrSges(4,2) * t269 - mrSges(4,3) * t238 + Ifges(4,1) * t289 + Ifges(4,4) * t288 + Ifges(4,5) * qJDD(2) - pkin(8) * t167 - qJD(2) * t274 - t154 * t329 + t156 * t332 - t273 * t303;
t345 = mrSges(7,1) * t195 - mrSges(7,3) * t193 - Ifges(7,4) * t227 - Ifges(7,2) * t287 - Ifges(7,6) * t226 + t228 * t263 - t232 * t262;
t340 = mrSges(6,2) * t198 - t262 * t233 - qJ(6) * (-t226 * mrSges(7,2) - t262 * t241 + t359) - pkin(5) * (-t227 * mrSges(7,2) - t263 * t241 + t353) - mrSges(6,1) * t197 - t263 * t231 + Ifges(6,6) * t226 - Ifges(6,5) * t227 - Ifges(6,3) * t287 + t345;
t337 = mrSges(5,1) * t206 - mrSges(5,2) * t207 + Ifges(5,5) * t257 + Ifges(5,6) * t256 + Ifges(5,3) * t287 + pkin(4) * t174 + t250 * t293 - t251 * t292 - t340;
t152 = -mrSges(4,1) * t269 + mrSges(4,3) * t239 + Ifges(4,4) * t289 + Ifges(4,2) * t288 + Ifges(4,6) * qJDD(2) - pkin(3) * t167 + qJD(2) * t275 - t273 * t304 - t337;
t305 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t330 + Ifges(3,6) * t333) * qJD(1);
t308 = -t336 * pkin(7) + t349;
t344 = m(4) * t269 - mrSges(4,1) * t288 + mrSges(4,2) * t289 + t296 * t303 + t297 * t304 + t167;
t144 = -mrSges(3,1) * t308 + mrSges(3,3) * t295 + Ifges(3,4) * t314 + Ifges(3,2) * t315 + Ifges(3,6) * qJDD(2) - pkin(2) * t344 + qJ(3) * t357 + qJD(2) * t307 + t327 * t148 + t328 * t152 - t305 * t363;
t147 = mrSges(3,2) * t308 - mrSges(3,3) * t294 + Ifges(3,1) * t314 + Ifges(3,4) * t315 + Ifges(3,5) * qJDD(2) - qJ(3) * t160 - qJD(2) * t306 + t148 * t328 - t152 * t327 + t305 * t362;
t341 = -m(3) * t308 + mrSges(3,1) * t315 - mrSges(3,2) * t314 - t318 * t363 + t319 * t362 - t344;
t346 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t341 + pkin(7) * t358 + t144 * t333 + t147 * t330;
t161 = m(2) * t320 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t336 + t341;
t151 = t158 * t333 + t159 * t330;
t149 = m(2) * t321 - mrSges(2,1) * t336 - qJDD(1) * mrSges(2,2) + t358;
t145 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t336 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t151 - t371;
t142 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t336 - pkin(7) * t151 - t144 * t330 + t147 * t333;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t334 * t142 - t331 * t145 - pkin(6) * (t149 * t331 + t161 * t334), t142, t147, t148, t156, t176, -t230 * t262 + t347; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t331 * t142 + t334 * t145 + pkin(6) * (t149 * t334 - t161 * t331), t145, t144, t152, t154, t175, -t345; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t346, t346, t371, -t343, t337, -t340, Ifges(7,5) * t227 + Ifges(7,6) * t287 + Ifges(7,3) * t226 + t263 * t230 - t302 * t232 - t352;];
m_new  = t1;
