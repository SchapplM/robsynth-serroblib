% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:36:39
% EndTime: 2019-05-06 17:37:31
% DurationCPUTime: 18.58s
% Computational Cost: add. (310802->382), mult. (705822->469), div. (0->0), fcn. (505783->10), ass. (0->142)
t333 = sin(qJ(2));
t336 = cos(qJ(2));
t361 = qJD(1) * qJD(2);
t316 = qJDD(1) * t333 + t336 * t361;
t334 = sin(qJ(1));
t337 = cos(qJ(1));
t323 = -g(1) * t337 - g(2) * t334;
t339 = qJD(1) ^ 2;
t311 = -pkin(1) * t339 + qJDD(1) * pkin(7) + t323;
t367 = t311 * t333;
t369 = pkin(2) * t339;
t264 = qJDD(2) * pkin(2) - qJ(3) * t316 - t367 + (qJ(3) * t361 + t333 * t369 - g(3)) * t336;
t295 = -g(3) * t333 + t336 * t311;
t317 = qJDD(1) * t336 - t333 * t361;
t364 = qJD(1) * t333;
t319 = qJD(2) * pkin(2) - qJ(3) * t364;
t328 = t336 ^ 2;
t265 = qJ(3) * t317 - qJD(2) * t319 - t328 * t369 + t295;
t329 = sin(pkin(10));
t330 = cos(pkin(10));
t305 = (t329 * t336 + t330 * t333) * qJD(1);
t238 = -0.2e1 * qJD(3) * t305 + t264 * t330 - t329 * t265;
t363 = qJD(1) * t336;
t304 = -t329 * t364 + t330 * t363;
t239 = 0.2e1 * qJD(3) * t304 + t329 * t264 + t330 * t265;
t276 = -mrSges(4,1) * t304 + mrSges(4,2) * t305;
t288 = -t329 * t316 + t317 * t330;
t297 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t305;
t282 = -pkin(3) * t304 - pkin(8) * t305;
t338 = qJD(2) ^ 2;
t223 = -pkin(3) * t338 + qJDD(2) * pkin(8) + t282 * t304 + t239;
t322 = g(1) * t334 - t337 * g(2);
t351 = -qJDD(1) * pkin(1) - t322;
t268 = -pkin(2) * t317 + qJDD(3) + t319 * t364 + (-qJ(3) * t328 - pkin(7)) * t339 + t351;
t289 = t316 * t330 + t317 * t329;
t227 = (-qJD(2) * t304 - t289) * pkin(8) + (qJD(2) * t305 - t288) * pkin(3) + t268;
t332 = sin(qJ(4));
t335 = cos(qJ(4));
t206 = -t223 * t332 + t335 * t227;
t292 = qJD(2) * t335 - t305 * t332;
t257 = qJD(4) * t292 + qJDD(2) * t332 + t289 * t335;
t287 = qJDD(4) - t288;
t293 = qJD(2) * t332 + t305 * t335;
t303 = qJD(4) - t304;
t202 = (t292 * t303 - t257) * pkin(9) + (t292 * t293 + t287) * pkin(4) + t206;
t207 = t335 * t223 + t332 * t227;
t256 = -qJD(4) * t293 + qJDD(2) * t335 - t289 * t332;
t271 = pkin(4) * t303 - pkin(9) * t293;
t291 = t292 ^ 2;
t204 = -pkin(4) * t291 + pkin(9) * t256 - t271 * t303 + t207;
t331 = sin(qJ(5));
t370 = cos(qJ(5));
t198 = t331 * t202 + t370 * t204;
t262 = t331 * t292 + t370 * t293;
t217 = qJD(5) * t262 - t256 * t370 + t257 * t331;
t299 = qJD(5) + t303;
t247 = mrSges(6,1) * t299 - mrSges(6,3) * t262;
t261 = -t370 * t292 + t293 * t331;
t283 = qJDD(5) + t287;
t240 = pkin(5) * t261 - t262 * qJ(6);
t298 = t299 ^ 2;
t193 = -pkin(5) * t298 + qJ(6) * t283 + 0.2e1 * qJD(6) * t299 - t240 * t261 + t198;
t248 = -mrSges(7,1) * t299 + mrSges(7,2) * t262;
t360 = m(7) * t193 + t283 * mrSges(7,3) + t299 * t248;
t241 = mrSges(7,1) * t261 - mrSges(7,3) * t262;
t365 = -mrSges(6,1) * t261 - mrSges(6,2) * t262 - t241;
t368 = -mrSges(6,3) - mrSges(7,2);
t179 = m(6) * t198 - mrSges(6,2) * t283 + t217 * t368 - t247 * t299 + t261 * t365 + t360;
t197 = t370 * t202 - t331 * t204;
t218 = -t261 * qJD(5) + t331 * t256 + t257 * t370;
t246 = -mrSges(6,2) * t299 - mrSges(6,3) * t261;
t195 = -t283 * pkin(5) - t298 * qJ(6) + t262 * t240 + qJDD(6) - t197;
t245 = -mrSges(7,2) * t261 + mrSges(7,3) * t299;
t354 = -m(7) * t195 + t283 * mrSges(7,1) + t299 * t245;
t183 = m(6) * t197 + mrSges(6,1) * t283 + t218 * t368 + t246 * t299 + t262 * t365 + t354;
t176 = t331 * t179 + t370 * t183;
t266 = -mrSges(5,1) * t292 + mrSges(5,2) * t293;
t269 = -mrSges(5,2) * t303 + mrSges(5,3) * t292;
t172 = m(5) * t206 + mrSges(5,1) * t287 - mrSges(5,3) * t257 - t266 * t293 + t269 * t303 + t176;
t270 = mrSges(5,1) * t303 - mrSges(5,3) * t293;
t356 = t370 * t179 - t183 * t331;
t173 = m(5) * t207 - mrSges(5,2) * t287 + mrSges(5,3) * t256 + t266 * t292 - t270 * t303 + t356;
t357 = -t172 * t332 + t335 * t173;
t165 = m(4) * t239 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t288 - qJD(2) * t297 + t276 * t304 + t357;
t296 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t304;
t222 = -qJDD(2) * pkin(3) - pkin(8) * t338 + t305 * t282 - t238;
t205 = -pkin(4) * t256 - pkin(9) * t291 + t293 * t271 + t222;
t200 = -0.2e1 * qJD(6) * t262 + (t261 * t299 - t218) * qJ(6) + (t262 * t299 + t217) * pkin(5) + t205;
t190 = m(7) * t200 + t217 * mrSges(7,1) - t218 * mrSges(7,3) + t261 * t245 - t262 * t248;
t345 = m(6) * t205 + t217 * mrSges(6,1) + mrSges(6,2) * t218 + t261 * t246 + t247 * t262 + t190;
t342 = -m(5) * t222 + t256 * mrSges(5,1) - mrSges(5,2) * t257 + t292 * t269 - t270 * t293 - t345;
t181 = m(4) * t238 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t289 + qJD(2) * t296 - t276 * t305 + t342;
t160 = t329 * t165 + t330 * t181;
t294 = -g(3) * t336 - t367;
t307 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t333 + Ifges(3,2) * t336) * qJD(1);
t308 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t333 + Ifges(3,4) * t336) * qJD(1);
t232 = Ifges(7,1) * t262 + Ifges(7,4) * t299 + Ifges(7,5) * t261;
t233 = Ifges(6,1) * t262 - Ifges(6,4) * t261 + Ifges(6,5) * t299;
t353 = -mrSges(7,1) * t200 + mrSges(7,2) * t193;
t230 = Ifges(7,4) * t262 + Ifges(7,2) * t299 + Ifges(7,6) * t261;
t366 = -Ifges(6,5) * t262 + Ifges(6,6) * t261 - Ifges(6,3) * t299 - t230;
t174 = -mrSges(6,1) * t205 + mrSges(6,3) * t198 - pkin(5) * t190 + (t232 + t233) * t299 + (Ifges(6,6) - Ifges(7,6)) * t283 + t366 * t262 + (Ifges(6,4) - Ifges(7,5)) * t218 + (-Ifges(6,2) - Ifges(7,3)) * t217 + t353;
t231 = Ifges(6,4) * t262 - Ifges(6,2) * t261 + Ifges(6,6) * t299;
t228 = Ifges(7,5) * t262 + Ifges(7,6) * t299 + Ifges(7,3) * t261;
t350 = mrSges(7,2) * t195 - mrSges(7,3) * t200 + Ifges(7,1) * t218 + Ifges(7,4) * t283 + Ifges(7,5) * t217 + t299 * t228;
t175 = mrSges(6,2) * t205 - mrSges(6,3) * t197 + Ifges(6,1) * t218 - Ifges(6,4) * t217 + Ifges(6,5) * t283 - qJ(6) * t190 - t231 * t299 + t261 * t366 + t350;
t249 = Ifges(5,5) * t293 + Ifges(5,6) * t292 + Ifges(5,3) * t303;
t251 = Ifges(5,1) * t293 + Ifges(5,4) * t292 + Ifges(5,5) * t303;
t154 = -mrSges(5,1) * t222 + mrSges(5,3) * t207 + Ifges(5,4) * t257 + Ifges(5,2) * t256 + Ifges(5,6) * t287 - pkin(4) * t345 + pkin(9) * t356 + t174 * t370 + t331 * t175 - t293 * t249 + t303 * t251;
t250 = Ifges(5,4) * t293 + Ifges(5,2) * t292 + Ifges(5,6) * t303;
t156 = mrSges(5,2) * t222 - mrSges(5,3) * t206 + Ifges(5,1) * t257 + Ifges(5,4) * t256 + Ifges(5,5) * t287 - pkin(9) * t176 - t331 * t174 + t175 * t370 + t292 * t249 - t303 * t250;
t273 = Ifges(4,4) * t305 + Ifges(4,2) * t304 + Ifges(4,6) * qJD(2);
t274 = Ifges(4,1) * t305 + Ifges(4,4) * t304 + Ifges(4,5) * qJD(2);
t346 = -mrSges(4,1) * t238 + mrSges(4,2) * t239 - Ifges(4,5) * t289 - Ifges(4,6) * t288 - Ifges(4,3) * qJDD(2) - pkin(3) * t342 - pkin(8) * t357 - t335 * t154 - t332 * t156 - t305 * t273 + t304 * t274;
t371 = mrSges(3,1) * t294 - mrSges(3,2) * t295 + Ifges(3,5) * t316 + Ifges(3,6) * t317 + Ifges(3,3) * qJDD(2) + pkin(2) * t160 + (t307 * t333 - t308 * t336) * qJD(1) - t346;
t167 = t335 * t172 + t332 * t173;
t315 = (-mrSges(3,1) * t336 + mrSges(3,2) * t333) * qJD(1);
t321 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t363;
t158 = m(3) * t294 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t316 + qJD(2) * t321 - t315 * t364 + t160;
t320 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t364;
t358 = t330 * t165 - t181 * t329;
t159 = m(3) * t295 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t317 - qJD(2) * t320 + t315 * t363 + t358;
t359 = -t158 * t333 + t336 * t159;
t272 = Ifges(4,5) * t305 + Ifges(4,6) * t304 + Ifges(4,3) * qJD(2);
t148 = mrSges(4,2) * t268 - mrSges(4,3) * t238 + Ifges(4,1) * t289 + Ifges(4,4) * t288 + Ifges(4,5) * qJDD(2) - pkin(8) * t167 - qJD(2) * t273 - t154 * t332 + t156 * t335 + t272 * t304;
t348 = mrSges(7,1) * t195 - mrSges(7,3) * t193 - Ifges(7,4) * t218 - Ifges(7,2) * t283 - Ifges(7,6) * t217 + t262 * t228 - t261 * t232;
t343 = mrSges(6,2) * t198 - t261 * t233 - qJ(6) * (-mrSges(7,2) * t217 - t241 * t261 + t360) - pkin(5) * (-mrSges(7,2) * t218 - t241 * t262 + t354) - mrSges(6,1) * t197 + Ifges(6,6) * t217 - Ifges(6,5) * t218 - t262 * t231 - Ifges(6,3) * t283 + t348;
t340 = mrSges(5,1) * t206 - mrSges(5,2) * t207 + Ifges(5,5) * t257 + Ifges(5,6) * t256 + Ifges(5,3) * t287 + pkin(4) * t176 + t293 * t250 - t292 * t251 - t343;
t152 = -mrSges(4,1) * t268 + mrSges(4,3) * t239 + Ifges(4,4) * t289 + Ifges(4,2) * t288 + Ifges(4,6) * qJDD(2) - pkin(3) * t167 + qJD(2) * t274 - t305 * t272 - t340;
t306 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t333 + Ifges(3,6) * t336) * qJD(1);
t310 = -pkin(7) * t339 + t351;
t347 = m(4) * t268 - t288 * mrSges(4,1) + mrSges(4,2) * t289 - t304 * t296 + t297 * t305 + t167;
t144 = -mrSges(3,1) * t310 + mrSges(3,3) * t295 + Ifges(3,4) * t316 + Ifges(3,2) * t317 + Ifges(3,6) * qJDD(2) - pkin(2) * t347 + qJ(3) * t358 + qJD(2) * t308 + t329 * t148 + t330 * t152 - t306 * t364;
t147 = mrSges(3,2) * t310 - mrSges(3,3) * t294 + Ifges(3,1) * t316 + Ifges(3,4) * t317 + Ifges(3,5) * qJDD(2) - qJ(3) * t160 - qJD(2) * t307 + t148 * t330 - t152 * t329 + t306 * t363;
t344 = -m(3) * t310 + t317 * mrSges(3,1) - mrSges(3,2) * t316 - t320 * t364 + t321 * t363 - t347;
t349 = mrSges(2,1) * t322 - mrSges(2,2) * t323 + Ifges(2,3) * qJDD(1) + pkin(1) * t344 + pkin(7) * t359 + t336 * t144 + t333 * t147;
t161 = m(2) * t322 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t339 + t344;
t151 = t158 * t336 + t159 * t333;
t149 = m(2) * t323 - mrSges(2,1) * t339 - qJDD(1) * mrSges(2,2) + t359;
t145 = mrSges(2,1) * g(3) + mrSges(2,3) * t323 + t339 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t151 - t371;
t142 = -mrSges(2,2) * g(3) - mrSges(2,3) * t322 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t339 - pkin(7) * t151 - t144 * t333 + t147 * t336;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t337 * t142 - t334 * t145 - pkin(6) * (t149 * t334 + t161 * t337), t142, t147, t148, t156, t175, -t230 * t261 + t350; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t334 * t142 + t337 * t145 + pkin(6) * (t149 * t337 - t161 * t334), t145, t144, t152, t154, t174, -t348; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t349, t349, t371, -t346, t340, -t343, Ifges(7,5) * t218 + Ifges(7,6) * t283 + Ifges(7,3) * t217 + t230 * t262 - t232 * t299 - t353;];
m_new  = t1;
