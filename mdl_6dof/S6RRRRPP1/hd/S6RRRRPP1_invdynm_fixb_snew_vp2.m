% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:00
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:54:05
% EndTime: 2019-05-07 17:54:44
% DurationCPUTime: 19.75s
% Computational Cost: add. (370247->383), mult. (748275->470), div. (0->0), fcn. (539235->10), ass. (0->144)
t326 = sin(qJ(3));
t327 = sin(qJ(2));
t330 = cos(qJ(3));
t331 = cos(qJ(2));
t301 = (t326 * t327 - t330 * t331) * qJD(1);
t356 = qJD(1) * qJD(2);
t309 = qJDD(1) * t327 + t331 * t356;
t328 = sin(qJ(1));
t332 = cos(qJ(1));
t316 = -g(1) * t332 - g(2) * t328;
t333 = qJD(1) ^ 2;
t304 = -pkin(1) * t333 + qJDD(1) * pkin(7) + t316;
t361 = t327 * t304;
t364 = pkin(2) * t333;
t263 = qJDD(2) * pkin(2) - t309 * pkin(8) - t361 + (pkin(8) * t356 + t327 * t364 - g(3)) * t331;
t292 = -g(3) * t327 + t331 * t304;
t310 = qJDD(1) * t331 - t327 * t356;
t358 = qJD(1) * t327;
t314 = qJD(2) * pkin(2) - pkin(8) * t358;
t323 = t331 ^ 2;
t264 = pkin(8) * t310 - qJD(2) * t314 - t323 * t364 + t292;
t239 = t326 * t263 + t330 * t264;
t302 = (t326 * t331 + t327 * t330) * qJD(1);
t275 = -qJD(3) * t302 - t309 * t326 + t310 * t330;
t285 = mrSges(4,1) * t301 + mrSges(4,2) * t302;
t321 = qJD(2) + qJD(3);
t294 = mrSges(4,1) * t321 - mrSges(4,3) * t302;
t320 = qJDD(2) + qJDD(3);
t276 = -qJD(3) * t301 + t309 * t330 + t310 * t326;
t315 = t328 * g(1) - t332 * g(2);
t346 = -qJDD(1) * pkin(1) - t315;
t277 = -t310 * pkin(2) + t314 * t358 + (-pkin(8) * t323 - pkin(7)) * t333 + t346;
t209 = (t301 * t321 - t276) * pkin(9) + (t302 * t321 - t275) * pkin(3) + t277;
t286 = pkin(3) * t301 - pkin(9) * t302;
t319 = t321 ^ 2;
t223 = -pkin(3) * t319 + pkin(9) * t320 - t286 * t301 + t239;
t325 = sin(qJ(4));
t329 = cos(qJ(4));
t202 = t329 * t209 - t325 * t223;
t289 = -t302 * t325 + t321 * t329;
t245 = qJD(4) * t289 + t276 * t329 + t320 * t325;
t274 = qJDD(4) - t275;
t290 = t302 * t329 + t321 * t325;
t297 = qJD(4) + t301;
t199 = (t289 * t297 - t245) * qJ(5) + (t289 * t290 + t274) * pkin(4) + t202;
t203 = t325 * t209 + t329 * t223;
t244 = -qJD(4) * t290 - t276 * t325 + t320 * t329;
t279 = pkin(4) * t297 - qJ(5) * t290;
t288 = t289 ^ 2;
t201 = -pkin(4) * t288 + qJ(5) * t244 - t279 * t297 + t203;
t324 = sin(pkin(10));
t362 = cos(pkin(10));
t256 = -t362 * t289 + t290 * t324;
t365 = -2 * qJD(5);
t195 = t324 * t199 + t362 * t201 + t256 * t365;
t219 = -t362 * t244 + t245 * t324;
t257 = t324 * t289 + t362 * t290;
t248 = mrSges(6,1) * t297 - mrSges(6,3) * t257;
t234 = pkin(5) * t256 - qJ(6) * t257;
t296 = t297 ^ 2;
t190 = -pkin(5) * t296 + qJ(6) * t274 + 0.2e1 * qJD(6) * t297 - t234 * t256 + t195;
t249 = -mrSges(7,1) * t297 + mrSges(7,2) * t257;
t355 = m(7) * t190 + t274 * mrSges(7,3) + t297 * t249;
t235 = mrSges(7,1) * t256 - mrSges(7,3) * t257;
t359 = -mrSges(6,1) * t256 - mrSges(6,2) * t257 - t235;
t363 = -mrSges(6,3) - mrSges(7,2);
t176 = m(6) * t195 - t274 * mrSges(6,2) + t363 * t219 - t297 * t248 + t359 * t256 + t355;
t345 = t362 * t199 - t324 * t201;
t194 = t257 * t365 + t345;
t220 = t324 * t244 + t362 * t245;
t247 = -mrSges(6,2) * t297 - mrSges(6,3) * t256;
t192 = -t274 * pkin(5) - t296 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t234) * t257 - t345;
t246 = -mrSges(7,2) * t256 + mrSges(7,3) * t297;
t350 = -m(7) * t192 + t274 * mrSges(7,1) + t297 * t246;
t178 = m(6) * t194 + t274 * mrSges(6,1) + t363 * t220 + t297 * t247 + t359 * t257 + t350;
t171 = t324 * t176 + t362 * t178;
t261 = -mrSges(5,1) * t289 + mrSges(5,2) * t290;
t278 = -mrSges(5,2) * t297 + mrSges(5,3) * t289;
t169 = m(5) * t202 + mrSges(5,1) * t274 - mrSges(5,3) * t245 - t261 * t290 + t278 * t297 + t171;
t280 = mrSges(5,1) * t297 - mrSges(5,3) * t290;
t351 = t362 * t176 - t178 * t324;
t170 = m(5) * t203 - mrSges(5,2) * t274 + mrSges(5,3) * t244 + t261 * t289 - t280 * t297 + t351;
t352 = -t169 * t325 + t329 * t170;
t162 = m(4) * t239 - mrSges(4,2) * t320 + mrSges(4,3) * t275 - t285 * t301 - t294 * t321 + t352;
t238 = t330 * t263 - t326 * t264;
t293 = -mrSges(4,2) * t321 - mrSges(4,3) * t301;
t222 = -t320 * pkin(3) - t319 * pkin(9) + t302 * t286 - t238;
t204 = -t244 * pkin(4) - t288 * qJ(5) + t290 * t279 + qJDD(5) + t222;
t197 = -0.2e1 * qJD(6) * t257 + (t256 * t297 - t220) * qJ(6) + (t257 * t297 + t219) * pkin(5) + t204;
t187 = m(7) * t197 + t219 * mrSges(7,1) - t220 * mrSges(7,3) + t256 * t246 - t257 * t249;
t339 = m(6) * t204 + t219 * mrSges(6,1) + mrSges(6,2) * t220 + t256 * t247 + t248 * t257 + t187;
t336 = -m(5) * t222 + t244 * mrSges(5,1) - mrSges(5,2) * t245 + t289 * t278 - t280 * t290 - t339;
t180 = m(4) * t238 + mrSges(4,1) * t320 - mrSges(4,3) * t276 - t285 * t302 + t293 * t321 + t336;
t157 = t326 * t162 + t330 * t180;
t291 = -t331 * g(3) - t361;
t299 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t327 + Ifges(3,2) * t331) * qJD(1);
t300 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t327 + Ifges(3,4) * t331) * qJD(1);
t229 = Ifges(7,1) * t257 + Ifges(7,4) * t297 + Ifges(7,5) * t256;
t230 = Ifges(6,1) * t257 - Ifges(6,4) * t256 + Ifges(6,5) * t297;
t349 = -mrSges(7,1) * t197 + mrSges(7,2) * t190;
t227 = Ifges(7,4) * t257 + Ifges(7,2) * t297 + Ifges(7,6) * t256;
t360 = -Ifges(6,5) * t257 + Ifges(6,6) * t256 - Ifges(6,3) * t297 - t227;
t172 = -mrSges(6,1) * t204 + mrSges(6,3) * t195 - pkin(5) * t187 + (t229 + t230) * t297 + (Ifges(6,6) - Ifges(7,6)) * t274 + t360 * t257 + (Ifges(6,4) - Ifges(7,5)) * t220 + (-Ifges(6,2) - Ifges(7,3)) * t219 + t349;
t228 = Ifges(6,4) * t257 - Ifges(6,2) * t256 + Ifges(6,6) * t297;
t225 = Ifges(7,5) * t257 + Ifges(7,6) * t297 + Ifges(7,3) * t256;
t344 = mrSges(7,2) * t192 - mrSges(7,3) * t197 + Ifges(7,1) * t220 + Ifges(7,4) * t274 + Ifges(7,5) * t219 + t297 * t225;
t173 = mrSges(6,2) * t204 - mrSges(6,3) * t194 + Ifges(6,1) * t220 - Ifges(6,4) * t219 + Ifges(6,5) * t274 - qJ(6) * t187 - t297 * t228 + t360 * t256 + t344;
t250 = Ifges(5,5) * t290 + Ifges(5,6) * t289 + Ifges(5,3) * t297;
t252 = Ifges(5,1) * t290 + Ifges(5,4) * t289 + Ifges(5,5) * t297;
t151 = -mrSges(5,1) * t222 + mrSges(5,3) * t203 + Ifges(5,4) * t245 + Ifges(5,2) * t244 + Ifges(5,6) * t274 - pkin(4) * t339 + qJ(5) * t351 + t362 * t172 + t324 * t173 - t290 * t250 + t297 * t252;
t251 = Ifges(5,4) * t290 + Ifges(5,2) * t289 + Ifges(5,6) * t297;
t153 = mrSges(5,2) * t222 - mrSges(5,3) * t202 + Ifges(5,1) * t245 + Ifges(5,4) * t244 + Ifges(5,5) * t274 - qJ(5) * t171 - t324 * t172 + t362 * t173 + t289 * t250 - t297 * t251;
t282 = Ifges(4,4) * t302 - Ifges(4,2) * t301 + Ifges(4,6) * t321;
t283 = Ifges(4,1) * t302 - Ifges(4,4) * t301 + Ifges(4,5) * t321;
t340 = -mrSges(4,1) * t238 + mrSges(4,2) * t239 - Ifges(4,5) * t276 - Ifges(4,6) * t275 - Ifges(4,3) * t320 - pkin(3) * t336 - pkin(9) * t352 - t329 * t151 - t325 * t153 - t302 * t282 - t301 * t283;
t366 = mrSges(3,1) * t291 - mrSges(3,2) * t292 + Ifges(3,5) * t309 + Ifges(3,6) * t310 + Ifges(3,3) * qJDD(2) + pkin(2) * t157 + (t299 * t327 - t300 * t331) * qJD(1) - t340;
t164 = t329 * t169 + t325 * t170;
t357 = qJD(1) * t331;
t308 = (-mrSges(3,1) * t331 + mrSges(3,2) * t327) * qJD(1);
t313 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t357;
t155 = m(3) * t291 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t309 + qJD(2) * t313 - t308 * t358 + t157;
t312 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t358;
t353 = t330 * t162 - t180 * t326;
t156 = m(3) * t292 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t310 - qJD(2) * t312 + t308 * t357 + t353;
t354 = -t155 * t327 + t331 * t156;
t281 = Ifges(4,5) * t302 - Ifges(4,6) * t301 + Ifges(4,3) * t321;
t145 = mrSges(4,2) * t277 - mrSges(4,3) * t238 + Ifges(4,1) * t276 + Ifges(4,4) * t275 + Ifges(4,5) * t320 - pkin(9) * t164 - t151 * t325 + t153 * t329 - t281 * t301 - t282 * t321;
t342 = mrSges(7,1) * t192 - mrSges(7,3) * t190 - Ifges(7,4) * t220 - Ifges(7,2) * t274 - Ifges(7,6) * t219 + t257 * t225 - t256 * t229;
t337 = mrSges(6,2) * t195 - t256 * t230 - qJ(6) * (-t219 * mrSges(7,2) - t256 * t235 + t355) - pkin(5) * (-t220 * mrSges(7,2) - t257 * t235 + t350) - mrSges(6,1) * t194 - t257 * t228 + Ifges(6,6) * t219 - Ifges(6,5) * t220 - Ifges(6,3) * t274 + t342;
t334 = mrSges(5,1) * t202 - mrSges(5,2) * t203 + Ifges(5,5) * t245 + Ifges(5,6) * t244 + Ifges(5,3) * t274 + pkin(4) * t171 + t290 * t251 - t289 * t252 - t337;
t149 = -mrSges(4,1) * t277 + mrSges(4,3) * t239 + Ifges(4,4) * t276 + Ifges(4,2) * t275 + Ifges(4,6) * t320 - pkin(3) * t164 - t302 * t281 + t321 * t283 - t334;
t298 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t327 + Ifges(3,6) * t331) * qJD(1);
t303 = -t333 * pkin(7) + t346;
t341 = m(4) * t277 - t275 * mrSges(4,1) + mrSges(4,2) * t276 + t301 * t293 + t294 * t302 + t164;
t141 = -mrSges(3,1) * t303 + mrSges(3,3) * t292 + Ifges(3,4) * t309 + Ifges(3,2) * t310 + Ifges(3,6) * qJDD(2) - pkin(2) * t341 + pkin(8) * t353 + qJD(2) * t300 + t326 * t145 + t330 * t149 - t298 * t358;
t144 = mrSges(3,2) * t303 - mrSges(3,3) * t291 + Ifges(3,1) * t309 + Ifges(3,4) * t310 + Ifges(3,5) * qJDD(2) - pkin(8) * t157 - qJD(2) * t299 + t145 * t330 - t149 * t326 + t298 * t357;
t338 = -m(3) * t303 + t310 * mrSges(3,1) - mrSges(3,2) * t309 - t312 * t358 + t313 * t357 - t341;
t343 = mrSges(2,1) * t315 - mrSges(2,2) * t316 + Ifges(2,3) * qJDD(1) + pkin(1) * t338 + pkin(7) * t354 + t331 * t141 + t327 * t144;
t158 = m(2) * t315 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t333 + t338;
t148 = t155 * t331 + t156 * t327;
t146 = m(2) * t316 - mrSges(2,1) * t333 - qJDD(1) * mrSges(2,2) + t354;
t142 = mrSges(2,1) * g(3) + mrSges(2,3) * t316 + t333 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t148 - t366;
t139 = -mrSges(2,2) * g(3) - mrSges(2,3) * t315 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t333 - pkin(7) * t148 - t141 * t327 + t144 * t331;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t332 * t139 - t328 * t142 - pkin(6) * (t146 * t328 + t158 * t332), t139, t144, t145, t153, t173, -t227 * t256 + t344; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t328 * t139 + t332 * t142 + pkin(6) * (t146 * t332 - t158 * t328), t142, t141, t149, t151, t172, -t342; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t343, t343, t366, -t340, t334, -t337, Ifges(7,5) * t220 + Ifges(7,6) * t274 + Ifges(7,3) * t219 + t257 * t227 - t297 * t229 - t349;];
m_new  = t1;
