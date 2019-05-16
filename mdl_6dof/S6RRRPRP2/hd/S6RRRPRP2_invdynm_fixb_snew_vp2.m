% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:30:36
% EndTime: 2019-05-07 07:31:13
% DurationCPUTime: 20.64s
% Computational Cost: add. (352873->382), mult. (792702->470), div. (0->0), fcn. (580361->10), ass. (0->142)
t332 = sin(qJ(2));
t335 = cos(qJ(2));
t359 = qJD(1) * qJD(2);
t310 = qJDD(1) * t332 + t335 * t359;
t333 = sin(qJ(1));
t336 = cos(qJ(1));
t317 = -g(1) * t336 - g(2) * t333;
t337 = qJD(1) ^ 2;
t305 = -pkin(1) * t337 + qJDD(1) * pkin(7) + t317;
t366 = t332 * t305;
t368 = pkin(2) * t337;
t268 = qJDD(2) * pkin(2) - t310 * pkin(8) - t366 + (pkin(8) * t359 + t332 * t368 - g(3)) * t335;
t293 = -g(3) * t332 + t335 * t305;
t311 = qJDD(1) * t335 - t332 * t359;
t362 = qJD(1) * t332;
t315 = qJD(2) * pkin(2) - pkin(8) * t362;
t327 = t335 ^ 2;
t269 = pkin(8) * t311 - qJD(2) * t315 - t327 * t368 + t293;
t331 = sin(qJ(3));
t334 = cos(qJ(3));
t237 = t334 * t268 - t331 * t269;
t302 = (-t331 * t332 + t334 * t335) * qJD(1);
t276 = qJD(3) * t302 + t310 * t334 + t311 * t331;
t303 = (t331 * t335 + t332 * t334) * qJD(1);
t324 = qJDD(2) + qJDD(3);
t325 = qJD(2) + qJD(3);
t208 = (t302 * t325 - t276) * qJ(4) + (t302 * t303 + t324) * pkin(3) + t237;
t238 = t331 * t268 + t334 * t269;
t275 = -qJD(3) * t303 - t310 * t331 + t311 * t334;
t295 = pkin(3) * t325 - qJ(4) * t303;
t298 = t302 ^ 2;
t211 = -pkin(3) * t298 + qJ(4) * t275 - t295 * t325 + t238;
t328 = sin(pkin(10));
t329 = cos(pkin(10));
t290 = t302 * t328 + t303 * t329;
t203 = -0.2e1 * qJD(4) * t290 + t329 * t208 - t328 * t211;
t289 = t302 * t329 - t303 * t328;
t204 = 0.2e1 * qJD(4) * t289 + t328 * t208 + t329 * t211;
t263 = -pkin(4) * t289 - pkin(9) * t290;
t323 = t325 ^ 2;
t201 = -pkin(4) * t323 + pkin(9) * t324 + t263 * t289 + t204;
t316 = t333 * g(1) - t336 * g(2);
t349 = -qJDD(1) * pkin(1) - t316;
t277 = -t311 * pkin(2) + t315 * t362 + (-pkin(8) * t327 - pkin(7)) * t337 + t349;
t223 = -t275 * pkin(3) - t298 * qJ(4) + t303 * t295 + qJDD(4) + t277;
t248 = t275 * t329 - t276 * t328;
t249 = t275 * t328 + t276 * t329;
t206 = (-t289 * t325 - t249) * pkin(9) + (t290 * t325 - t248) * pkin(4) + t223;
t330 = sin(qJ(5));
t369 = cos(qJ(5));
t197 = -t330 * t201 + t369 * t206;
t198 = t369 * t201 + t330 * t206;
t274 = t369 * t290 + t330 * t325;
t220 = t274 * qJD(5) + t330 * t249 - t369 * t324;
t273 = t330 * t290 - t369 * t325;
t221 = -t273 * qJD(5) + t369 * t249 + t330 * t324;
t283 = qJD(5) - t289;
t224 = Ifges(7,5) * t274 + Ifges(7,6) * t283 + Ifges(7,3) * t273;
t227 = Ifges(6,4) * t274 - Ifges(6,2) * t273 + Ifges(6,6) * t283;
t229 = Ifges(6,1) * t274 - Ifges(6,4) * t273 + Ifges(6,5) * t283;
t247 = qJDD(5) - t248;
t251 = mrSges(7,1) * t273 - mrSges(7,3) * t274;
t250 = pkin(5) * t273 - qJ(6) * t274;
t282 = t283 ^ 2;
t193 = -pkin(5) * t282 + qJ(6) * t247 + 0.2e1 * qJD(6) * t283 - t250 * t273 + t198;
t195 = -t247 * pkin(5) - t282 * qJ(6) + t274 * t250 + qJDD(6) - t197;
t228 = Ifges(7,1) * t274 + Ifges(7,4) * t283 + Ifges(7,5) * t273;
t347 = mrSges(7,1) * t195 - mrSges(7,3) * t193 - Ifges(7,4) * t221 - Ifges(7,2) * t247 - Ifges(7,6) * t220 - t273 * t228;
t253 = -mrSges(7,2) * t273 + mrSges(7,3) * t283;
t352 = -m(7) * t195 + t247 * mrSges(7,1) + t283 * t253;
t256 = -mrSges(7,1) * t283 + mrSges(7,2) * t274;
t358 = m(7) * t193 + t247 * mrSges(7,3) + t283 * t256;
t371 = -(-t227 + t224) * t274 + mrSges(6,1) * t197 - mrSges(6,2) * t198 + Ifges(6,5) * t221 - Ifges(6,6) * t220 + Ifges(6,3) * t247 + pkin(5) * (-t221 * mrSges(7,2) - t274 * t251 + t352) + qJ(6) * (-t220 * mrSges(7,2) - t273 * t251 + t358) + t273 * t229 - t347;
t262 = -mrSges(5,1) * t289 + mrSges(5,2) * t290;
t279 = mrSges(5,1) * t325 - mrSges(5,3) * t290;
t255 = mrSges(6,1) * t283 - mrSges(6,3) * t274;
t363 = -mrSges(6,1) * t273 - mrSges(6,2) * t274 - t251;
t367 = -mrSges(6,3) - mrSges(7,2);
t183 = m(6) * t198 - t247 * mrSges(6,2) + t367 * t220 - t283 * t255 + t363 * t273 + t358;
t254 = -mrSges(6,2) * t283 - mrSges(6,3) * t273;
t185 = m(6) * t197 + t247 * mrSges(6,1) + t367 * t221 + t283 * t254 + t363 * t274 + t352;
t354 = t369 * t183 - t185 * t330;
t171 = m(5) * t204 - mrSges(5,2) * t324 + mrSges(5,3) * t248 + t262 * t289 - t279 * t325 + t354;
t278 = -mrSges(5,2) * t325 + mrSges(5,3) * t289;
t200 = -t324 * pkin(4) - t323 * pkin(9) + t290 * t263 - t203;
t196 = -0.2e1 * qJD(6) * t274 + (t273 * t283 - t221) * qJ(6) + (t274 * t283 + t220) * pkin(5) + t200;
t190 = m(7) * t196 + mrSges(7,1) * t220 - t221 * mrSges(7,3) + t253 * t273 - t274 * t256;
t342 = -m(6) * t200 - t220 * mrSges(6,1) - mrSges(6,2) * t221 - t273 * t254 - t255 * t274 - t190;
t180 = m(5) * t203 + mrSges(5,1) * t324 - mrSges(5,3) * t249 - t262 * t290 + t278 * t325 + t342;
t166 = t328 * t171 + t329 * t180;
t291 = -mrSges(4,1) * t302 + mrSges(4,2) * t303;
t294 = -mrSges(4,2) * t325 + mrSges(4,3) * t302;
t163 = m(4) * t237 + mrSges(4,1) * t324 - mrSges(4,3) * t276 - t291 * t303 + t294 * t325 + t166;
t296 = mrSges(4,1) * t325 - mrSges(4,3) * t303;
t355 = t329 * t171 - t180 * t328;
t164 = m(4) * t238 - mrSges(4,2) * t324 + mrSges(4,3) * t275 + t291 * t302 - t296 * t325 + t355;
t157 = t334 * t163 + t331 * t164;
t292 = -t335 * g(3) - t366;
t300 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t332 + Ifges(3,2) * t335) * qJD(1);
t301 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t332 + Ifges(3,4) * t335) * qJD(1);
t285 = Ifges(4,4) * t303 + Ifges(4,2) * t302 + Ifges(4,6) * t325;
t286 = Ifges(4,1) * t303 + Ifges(4,4) * t302 + Ifges(4,5) * t325;
t351 = -mrSges(7,1) * t196 + mrSges(7,2) * t193;
t226 = Ifges(7,4) * t274 + Ifges(7,2) * t283 + Ifges(7,6) * t273;
t365 = -Ifges(6,5) * t274 + Ifges(6,6) * t273 - Ifges(6,3) * t283 - t226;
t173 = -mrSges(6,1) * t200 + mrSges(6,3) * t198 - pkin(5) * t190 + (t228 + t229) * t283 + t365 * t274 + (Ifges(6,6) - Ifges(7,6)) * t247 + (Ifges(6,4) - Ifges(7,5)) * t221 + (-Ifges(6,2) - Ifges(7,3)) * t220 + t351;
t346 = mrSges(7,2) * t195 - mrSges(7,3) * t196 + Ifges(7,1) * t221 + Ifges(7,4) * t247 + Ifges(7,5) * t220 + t283 * t224;
t175 = mrSges(6,2) * t200 - mrSges(6,3) * t197 + Ifges(6,1) * t221 - Ifges(6,4) * t220 + Ifges(6,5) * t247 - qJ(6) * t190 - t283 * t227 + t365 * t273 + t346;
t258 = Ifges(5,4) * t290 + Ifges(5,2) * t289 + Ifges(5,6) * t325;
t259 = Ifges(5,1) * t290 + Ifges(5,4) * t289 + Ifges(5,5) * t325;
t344 = -mrSges(5,1) * t203 + mrSges(5,2) * t204 - Ifges(5,5) * t249 - Ifges(5,6) * t248 - Ifges(5,3) * t324 - pkin(4) * t342 - pkin(9) * t354 - t369 * t173 - t330 * t175 - t290 * t258 + t289 * t259;
t341 = -mrSges(4,1) * t237 + mrSges(4,2) * t238 - Ifges(4,5) * t276 - Ifges(4,6) * t275 - Ifges(4,3) * t324 - pkin(3) * t166 - t303 * t285 + t302 * t286 + t344;
t370 = mrSges(3,1) * t292 - mrSges(3,2) * t293 + Ifges(3,5) * t310 + Ifges(3,6) * t311 + Ifges(3,3) * qJDD(2) + pkin(2) * t157 + (t300 * t332 - t301 * t335) * qJD(1) - t341;
t177 = t330 * t183 + t369 * t185;
t361 = qJD(1) * t335;
t309 = (-mrSges(3,1) * t335 + mrSges(3,2) * t332) * qJD(1);
t314 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t361;
t155 = m(3) * t292 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t310 + qJD(2) * t314 - t309 * t362 + t157;
t313 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t362;
t356 = -t163 * t331 + t334 * t164;
t156 = m(3) * t293 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t311 - qJD(2) * t313 + t309 * t361 + t356;
t357 = -t155 * t332 + t335 * t156;
t348 = m(5) * t223 - t248 * mrSges(5,1) + t249 * mrSges(5,2) - t289 * t278 + t290 * t279 + t177;
t257 = Ifges(5,5) * t290 + Ifges(5,6) * t289 + Ifges(5,3) * t325;
t158 = mrSges(5,2) * t223 - mrSges(5,3) * t203 + Ifges(5,1) * t249 + Ifges(5,4) * t248 + Ifges(5,5) * t324 - pkin(9) * t177 - t330 * t173 + t369 * t175 + t289 * t257 - t325 * t258;
t159 = -mrSges(5,1) * t223 + mrSges(5,3) * t204 + Ifges(5,4) * t249 + Ifges(5,2) * t248 + Ifges(5,6) * t324 - pkin(4) * t177 - t290 * t257 + t325 * t259 - t371;
t284 = Ifges(4,5) * t303 + Ifges(4,6) * t302 + Ifges(4,3) * t325;
t149 = -mrSges(4,1) * t277 + mrSges(4,3) * t238 + Ifges(4,4) * t276 + Ifges(4,2) * t275 + Ifges(4,6) * t324 - pkin(3) * t348 + qJ(4) * t355 + t328 * t158 + t329 * t159 - t303 * t284 + t325 * t286;
t150 = mrSges(4,2) * t277 - mrSges(4,3) * t237 + Ifges(4,1) * t276 + Ifges(4,4) * t275 + Ifges(4,5) * t324 - qJ(4) * t166 + t158 * t329 - t159 * t328 + t284 * t302 - t285 * t325;
t299 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t332 + Ifges(3,6) * t335) * qJD(1);
t304 = -t337 * pkin(7) + t349;
t343 = m(4) * t277 - t275 * mrSges(4,1) + mrSges(4,2) * t276 - t302 * t294 + t296 * t303 + t348;
t145 = -mrSges(3,1) * t304 + mrSges(3,3) * t293 + Ifges(3,4) * t310 + Ifges(3,2) * t311 + Ifges(3,6) * qJDD(2) - pkin(2) * t343 + pkin(8) * t356 + qJD(2) * t301 + t334 * t149 + t331 * t150 - t299 * t362;
t147 = mrSges(3,2) * t304 - mrSges(3,3) * t292 + Ifges(3,1) * t310 + Ifges(3,4) * t311 + Ifges(3,5) * qJDD(2) - pkin(8) * t157 - qJD(2) * t300 - t149 * t331 + t150 * t334 + t299 * t361;
t340 = -m(3) * t304 + t311 * mrSges(3,1) - mrSges(3,2) * t310 - t313 * t362 + t314 * t361 - t343;
t345 = mrSges(2,1) * t316 - mrSges(2,2) * t317 + Ifges(2,3) * qJDD(1) + pkin(1) * t340 + pkin(7) * t357 + t335 * t145 + t332 * t147;
t167 = m(2) * t316 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t337 + t340;
t153 = t155 * t335 + t156 * t332;
t151 = m(2) * t317 - mrSges(2,1) * t337 - qJDD(1) * mrSges(2,2) + t357;
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t317 + t337 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t153 - t370;
t143 = -mrSges(2,2) * g(3) - mrSges(2,3) * t316 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t337 - pkin(7) * t153 - t145 * t332 + t147 * t335;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t336 * t143 - t333 * t148 - pkin(6) * (t151 * t333 + t167 * t336), t143, t147, t150, t158, t175, -t226 * t273 + t346; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t333 * t143 + t336 * t148 + pkin(6) * (t151 * t336 - t167 * t333), t148, t145, t149, t159, t173, -t274 * t224 - t347; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t345, t345, t370, -t341, -t344, t371, Ifges(7,5) * t221 + Ifges(7,6) * t247 + Ifges(7,3) * t220 + t274 * t226 - t283 * t228 - t351;];
m_new  = t1;
