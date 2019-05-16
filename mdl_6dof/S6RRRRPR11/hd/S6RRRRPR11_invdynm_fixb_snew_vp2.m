% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR11
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 23:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 23:22:55
% EndTime: 2019-05-07 23:25:58
% DurationCPUTime: 99.44s
% Computational Cost: add. (1862346->398), mult. (4007021->518), div. (0->0), fcn. (3246752->14), ass. (0->165)
t341 = sin(pkin(6));
t347 = sin(qJ(2));
t352 = cos(qJ(2));
t369 = qJD(1) * qJD(2);
t327 = (-qJDD(1) * t352 + t347 * t369) * t341;
t378 = pkin(8) * t341;
t343 = cos(pkin(6));
t377 = g(3) * t343;
t376 = t341 * t347;
t375 = t341 * t352;
t374 = t343 * t347;
t373 = t343 * t352;
t371 = qJD(1) * t341;
t325 = (-pkin(2) * t352 - pkin(9) * t347) * t371;
t336 = qJD(1) * t343 + qJD(2);
t334 = t336 ^ 2;
t335 = qJDD(1) * t343 + qJDD(2);
t370 = qJD(1) * t352;
t348 = sin(qJ(1));
t353 = cos(qJ(1));
t332 = t348 * g(1) - g(2) * t353;
t354 = qJD(1) ^ 2;
t322 = qJDD(1) * pkin(1) + t354 * t378 + t332;
t333 = -g(1) * t353 - g(2) * t348;
t323 = -pkin(1) * t354 + qJDD(1) * t378 + t333;
t372 = t322 * t374 + t352 * t323;
t273 = -pkin(2) * t334 + pkin(9) * t335 + (-g(3) * t347 + t325 * t370) * t341 + t372;
t326 = (qJDD(1) * t347 + t352 * t369) * t341;
t274 = pkin(2) * t327 - pkin(9) * t326 - t377 + (-t322 + (pkin(2) * t347 - pkin(9) * t352) * t336 * qJD(1)) * t341;
t346 = sin(qJ(3));
t351 = cos(qJ(3));
t250 = t351 * t273 + t346 * t274;
t368 = t347 * t371;
t314 = t336 * t351 - t346 * t368;
t315 = t336 * t346 + t351 * t368;
t299 = -pkin(3) * t314 - pkin(10) * t315;
t319 = qJDD(3) + t327;
t367 = t341 * t370;
t331 = qJD(3) - t367;
t329 = t331 ^ 2;
t237 = -pkin(3) * t329 + pkin(10) * t319 + t299 * t314 + t250;
t296 = -g(3) * t375 + t322 * t373 - t347 * t323;
t272 = -pkin(2) * t335 - pkin(9) * t334 + t325 * t368 - t296;
t294 = -t315 * qJD(3) - t346 * t326 + t335 * t351;
t295 = qJD(3) * t314 + t326 * t351 + t335 * t346;
t241 = (-t314 * t331 - t295) * pkin(10) + (t315 * t331 - t294) * pkin(3) + t272;
t345 = sin(qJ(4));
t350 = cos(qJ(4));
t225 = -t237 * t345 + t350 * t241;
t301 = -t315 * t345 + t331 * t350;
t260 = qJD(4) * t301 + t295 * t350 + t319 * t345;
t292 = qJDD(4) - t294;
t302 = t315 * t350 + t331 * t345;
t313 = qJD(4) - t314;
t216 = (t301 * t313 - t260) * qJ(5) + (t301 * t302 + t292) * pkin(4) + t225;
t226 = t350 * t237 + t345 * t241;
t259 = -qJD(4) * t302 - t295 * t345 + t319 * t350;
t283 = pkin(4) * t313 - qJ(5) * t302;
t300 = t301 ^ 2;
t218 = -pkin(4) * t300 + qJ(5) * t259 - t283 * t313 + t226;
t340 = sin(pkin(12));
t342 = cos(pkin(12));
t279 = t301 * t340 + t302 * t342;
t210 = -0.2e1 * qJD(5) * t279 + t342 * t216 - t218 * t340;
t246 = t259 * t340 + t260 * t342;
t278 = t301 * t342 - t302 * t340;
t207 = (t278 * t313 - t246) * pkin(11) + (t278 * t279 + t292) * pkin(5) + t210;
t211 = 0.2e1 * qJD(5) * t278 + t340 * t216 + t342 * t218;
t245 = t259 * t342 - t260 * t340;
t263 = pkin(5) * t313 - pkin(11) * t279;
t277 = t278 ^ 2;
t208 = -pkin(5) * t277 + pkin(11) * t245 - t263 * t313 + t211;
t344 = sin(qJ(6));
t349 = cos(qJ(6));
t205 = t207 * t349 - t208 * t344;
t255 = t278 * t349 - t279 * t344;
t224 = qJD(6) * t255 + t245 * t344 + t246 * t349;
t256 = t278 * t344 + t279 * t349;
t234 = -mrSges(7,1) * t255 + mrSges(7,2) * t256;
t311 = qJD(6) + t313;
t247 = -mrSges(7,2) * t311 + mrSges(7,3) * t255;
t287 = qJDD(6) + t292;
t199 = m(7) * t205 + mrSges(7,1) * t287 - mrSges(7,3) * t224 - t234 * t256 + t247 * t311;
t206 = t207 * t344 + t208 * t349;
t223 = -qJD(6) * t256 + t245 * t349 - t246 * t344;
t248 = mrSges(7,1) * t311 - mrSges(7,3) * t256;
t200 = m(7) * t206 - mrSges(7,2) * t287 + mrSges(7,3) * t223 + t234 * t255 - t248 * t311;
t193 = t349 * t199 + t344 * t200;
t257 = -mrSges(6,1) * t278 + mrSges(6,2) * t279;
t261 = -mrSges(6,2) * t313 + mrSges(6,3) * t278;
t190 = m(6) * t210 + mrSges(6,1) * t292 - mrSges(6,3) * t246 - t257 * t279 + t261 * t313 + t193;
t262 = mrSges(6,1) * t313 - mrSges(6,3) * t279;
t364 = -t199 * t344 + t349 * t200;
t191 = m(6) * t211 - mrSges(6,2) * t292 + mrSges(6,3) * t245 + t257 * t278 - t262 * t313 + t364;
t186 = t342 * t190 + t340 * t191;
t280 = -mrSges(5,1) * t301 + mrSges(5,2) * t302;
t282 = -mrSges(5,2) * t313 + mrSges(5,3) * t301;
t184 = m(5) * t225 + mrSges(5,1) * t292 - mrSges(5,3) * t260 - t280 * t302 + t282 * t313 + t186;
t284 = mrSges(5,1) * t313 - mrSges(5,3) * t302;
t365 = -t190 * t340 + t342 * t191;
t185 = m(5) * t226 - mrSges(5,2) * t292 + mrSges(5,3) * t259 + t280 * t301 - t284 * t313 + t365;
t180 = -t184 * t345 + t350 * t185;
t298 = -mrSges(4,1) * t314 + mrSges(4,2) * t315;
t304 = mrSges(4,1) * t331 - mrSges(4,3) * t315;
t178 = m(4) * t250 - mrSges(4,2) * t319 + mrSges(4,3) * t294 + t298 * t314 - t304 * t331 + t180;
t249 = -t346 * t273 + t274 * t351;
t236 = -pkin(3) * t319 - pkin(10) * t329 + t315 * t299 - t249;
t227 = -pkin(4) * t259 - qJ(5) * t300 + t302 * t283 + qJDD(5) + t236;
t213 = -pkin(5) * t245 - pkin(11) * t277 + t263 * t279 + t227;
t363 = m(7) * t213 - t223 * mrSges(7,1) + t224 * mrSges(7,2) - t255 * t247 + t256 * t248;
t359 = m(6) * t227 - t245 * mrSges(6,1) + mrSges(6,2) * t246 - t278 * t261 + t262 * t279 + t363;
t203 = -m(5) * t236 + t259 * mrSges(5,1) - mrSges(5,2) * t260 + t301 * t282 - t284 * t302 - t359;
t303 = -mrSges(4,2) * t331 + mrSges(4,3) * t314;
t202 = m(4) * t249 + mrSges(4,1) * t319 - mrSges(4,3) * t295 - t298 * t315 + t303 * t331 + t203;
t173 = t346 * t178 + t351 * t202;
t297 = -g(3) * t376 + t372;
t320 = mrSges(3,1) * t336 - mrSges(3,3) * t368;
t324 = (-mrSges(3,1) * t352 + mrSges(3,2) * t347) * t371;
t366 = t351 * t178 - t202 * t346;
t171 = m(3) * t297 - mrSges(3,2) * t335 - mrSges(3,3) * t327 - t320 * t336 + t324 * t367 + t366;
t321 = -mrSges(3,2) * t336 + mrSges(3,3) * t367;
t179 = t184 * t350 + t185 * t345;
t358 = -m(4) * t272 + t294 * mrSges(4,1) - mrSges(4,2) * t295 + t314 * t303 - t304 * t315 - t179;
t175 = m(3) * t296 + mrSges(3,1) * t335 - mrSges(3,3) * t326 + t321 * t336 - t324 * t368 + t358;
t165 = t352 * t171 - t175 * t347;
t308 = -t322 * t341 - t377;
t172 = m(3) * t308 + mrSges(3,1) * t327 + mrSges(3,2) * t326 + (t320 * t347 - t321 * t352) * t371 + t173;
t162 = t171 * t374 - t172 * t341 + t175 * t373;
t229 = Ifges(7,5) * t256 + Ifges(7,6) * t255 + Ifges(7,3) * t311;
t231 = Ifges(7,1) * t256 + Ifges(7,4) * t255 + Ifges(7,5) * t311;
t194 = -mrSges(7,1) * t213 + mrSges(7,3) * t206 + Ifges(7,4) * t224 + Ifges(7,2) * t223 + Ifges(7,6) * t287 - t229 * t256 + t231 * t311;
t230 = Ifges(7,4) * t256 + Ifges(7,2) * t255 + Ifges(7,6) * t311;
t195 = mrSges(7,2) * t213 - mrSges(7,3) * t205 + Ifges(7,1) * t224 + Ifges(7,4) * t223 + Ifges(7,5) * t287 + t229 * t255 - t230 * t311;
t251 = Ifges(6,5) * t279 + Ifges(6,6) * t278 + Ifges(6,3) * t313;
t253 = Ifges(6,1) * t279 + Ifges(6,4) * t278 + Ifges(6,5) * t313;
t181 = -mrSges(6,1) * t227 + mrSges(6,3) * t211 + Ifges(6,4) * t246 + Ifges(6,2) * t245 + Ifges(6,6) * t292 - pkin(5) * t363 + pkin(11) * t364 + t349 * t194 + t344 * t195 - t279 * t251 + t313 * t253;
t252 = Ifges(6,4) * t279 + Ifges(6,2) * t278 + Ifges(6,6) * t313;
t182 = mrSges(6,2) * t227 - mrSges(6,3) * t210 + Ifges(6,1) * t246 + Ifges(6,4) * t245 + Ifges(6,5) * t292 - pkin(11) * t193 - t194 * t344 + t195 * t349 + t251 * t278 - t252 * t313;
t264 = Ifges(5,5) * t302 + Ifges(5,6) * t301 + Ifges(5,3) * t313;
t266 = Ifges(5,1) * t302 + Ifges(5,4) * t301 + Ifges(5,5) * t313;
t167 = -mrSges(5,1) * t236 + mrSges(5,3) * t226 + Ifges(5,4) * t260 + Ifges(5,2) * t259 + Ifges(5,6) * t292 - pkin(4) * t359 + qJ(5) * t365 + t342 * t181 + t340 * t182 - t302 * t264 + t313 * t266;
t265 = Ifges(5,4) * t302 + Ifges(5,2) * t301 + Ifges(5,6) * t313;
t168 = mrSges(5,2) * t236 - mrSges(5,3) * t225 + Ifges(5,1) * t260 + Ifges(5,4) * t259 + Ifges(5,5) * t292 - qJ(5) * t186 - t181 * t340 + t182 * t342 + t264 * t301 - t265 * t313;
t288 = Ifges(4,5) * t315 + Ifges(4,6) * t314 + Ifges(4,3) * t331;
t289 = Ifges(4,4) * t315 + Ifges(4,2) * t314 + Ifges(4,6) * t331;
t158 = mrSges(4,2) * t272 - mrSges(4,3) * t249 + Ifges(4,1) * t295 + Ifges(4,4) * t294 + Ifges(4,5) * t319 - pkin(10) * t179 - t167 * t345 + t168 * t350 + t288 * t314 - t289 * t331;
t290 = Ifges(4,1) * t315 + Ifges(4,4) * t314 + Ifges(4,5) * t331;
t360 = -mrSges(7,1) * t205 + mrSges(7,2) * t206 - Ifges(7,5) * t224 - Ifges(7,6) * t223 - Ifges(7,3) * t287 - t256 * t230 + t255 * t231;
t357 = -mrSges(6,1) * t210 + mrSges(6,2) * t211 - Ifges(6,5) * t246 - Ifges(6,6) * t245 - Ifges(6,3) * t292 - pkin(5) * t193 - t279 * t252 + t278 * t253 + t360;
t355 = mrSges(5,1) * t225 - mrSges(5,2) * t226 + Ifges(5,5) * t260 + Ifges(5,6) * t259 + Ifges(5,3) * t292 + pkin(4) * t186 + t302 * t265 - t301 * t266 - t357;
t166 = -mrSges(4,1) * t272 + mrSges(4,3) * t250 + Ifges(4,4) * t295 + Ifges(4,2) * t294 + Ifges(4,6) * t319 - pkin(3) * t179 - t315 * t288 + t331 * t290 - t355;
t306 = Ifges(3,6) * t336 + (Ifges(3,4) * t347 + Ifges(3,2) * t352) * t371;
t307 = Ifges(3,5) * t336 + (Ifges(3,1) * t347 + Ifges(3,4) * t352) * t371;
t153 = Ifges(3,5) * t326 - Ifges(3,6) * t327 + Ifges(3,3) * t335 + mrSges(3,1) * t296 - mrSges(3,2) * t297 + t346 * t158 + t351 * t166 + pkin(2) * t358 + pkin(9) * t366 + (t306 * t347 - t307 * t352) * t371;
t305 = Ifges(3,3) * t336 + (Ifges(3,5) * t347 + Ifges(3,6) * t352) * t371;
t155 = mrSges(3,2) * t308 - mrSges(3,3) * t296 + Ifges(3,1) * t326 - Ifges(3,4) * t327 + Ifges(3,5) * t335 - pkin(9) * t173 + t158 * t351 - t166 * t346 + t305 * t367 - t306 * t336;
t356 = mrSges(4,1) * t249 - mrSges(4,2) * t250 + Ifges(4,5) * t295 + Ifges(4,6) * t294 + Ifges(4,3) * t319 + pkin(3) * t203 + pkin(10) * t180 + t350 * t167 + t345 * t168 + t315 * t289 - t314 * t290;
t157 = -mrSges(3,1) * t308 + mrSges(3,3) * t297 + Ifges(3,4) * t326 - Ifges(3,2) * t327 + Ifges(3,6) * t335 - pkin(2) * t173 - t305 * t368 + t336 * t307 - t356;
t361 = mrSges(2,1) * t332 - mrSges(2,2) * t333 + Ifges(2,3) * qJDD(1) + pkin(1) * t162 + t343 * t153 + t155 * t376 + t157 * t375 + t165 * t378;
t163 = m(2) * t333 - mrSges(2,1) * t354 - qJDD(1) * mrSges(2,2) + t165;
t161 = t172 * t343 + (t171 * t347 + t175 * t352) * t341;
t159 = m(2) * t332 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t354 + t162;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t332 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t354 + t155 * t352 - t157 * t347 + (-t161 * t341 - t162 * t343) * pkin(8);
t150 = mrSges(2,1) * g(3) + mrSges(2,3) * t333 + Ifges(2,5) * t354 + Ifges(2,6) * qJDD(1) - pkin(1) * t161 - t153 * t341 + (pkin(8) * t165 + t155 * t347 + t157 * t352) * t343;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t353 * t151 - t348 * t150 - pkin(7) * (t159 * t353 + t163 * t348), t151, t155, t158, t168, t182, t195; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t348 * t151 + t353 * t150 + pkin(7) * (-t159 * t348 + t163 * t353), t150, t157, t166, t167, t181, t194; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t361, t361, t153, t356, t355, -t357, -t360;];
m_new  = t1;
