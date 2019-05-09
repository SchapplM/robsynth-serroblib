% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 19:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:52:08
% EndTime: 2019-05-06 19:53:38
% DurationCPUTime: 45.18s
% Computational Cost: add. (766816->386), mult. (1771822->488), div. (0->0), fcn. (1341623->12), ass. (0->153)
t339 = sin(qJ(2));
t344 = cos(qJ(2));
t365 = qJD(1) * qJD(2);
t318 = t339 * qJDD(1) + t344 * t365;
t340 = sin(qJ(1));
t345 = cos(qJ(1));
t325 = -t345 * g(1) - t340 * g(2);
t346 = qJD(1) ^ 2;
t313 = -t346 * pkin(1) + qJDD(1) * pkin(7) + t325;
t368 = t339 * t313;
t369 = pkin(2) * t346;
t272 = qJDD(2) * pkin(2) - t318 * qJ(3) - t368 + (qJ(3) * t365 + t339 * t369 - g(3)) * t344;
t298 = -t339 * g(3) + t344 * t313;
t319 = t344 * qJDD(1) - t339 * t365;
t367 = qJD(1) * t339;
t321 = qJD(2) * pkin(2) - qJ(3) * t367;
t333 = t344 ^ 2;
t273 = t319 * qJ(3) - qJD(2) * t321 - t333 * t369 + t298;
t334 = sin(pkin(11));
t335 = cos(pkin(11));
t308 = (t334 * t344 + t335 * t339) * qJD(1);
t242 = -0.2e1 * qJD(3) * t308 + t335 * t272 - t334 * t273;
t296 = t335 * t318 + t334 * t319;
t307 = (-t334 * t339 + t335 * t344) * qJD(1);
t225 = (qJD(2) * t307 - t296) * pkin(8) + (t307 * t308 + qJDD(2)) * pkin(3) + t242;
t243 = 0.2e1 * qJD(3) * t307 + t334 * t272 + t335 * t273;
t295 = -t334 * t318 + t335 * t319;
t301 = qJD(2) * pkin(3) - t308 * pkin(8);
t306 = t307 ^ 2;
t229 = -t306 * pkin(3) + t295 * pkin(8) - qJD(2) * t301 + t243;
t338 = sin(qJ(4));
t343 = cos(qJ(4));
t218 = t338 * t225 + t343 * t229;
t287 = t338 * t307 + t343 * t308;
t253 = -t287 * qJD(4) + t343 * t295 - t338 * t296;
t286 = t343 * t307 - t338 * t308;
t266 = -t286 * mrSges(5,1) + t287 * mrSges(5,2);
t330 = qJD(2) + qJD(4);
t279 = t330 * mrSges(5,1) - t287 * mrSges(5,3);
t329 = qJDD(2) + qJDD(4);
t267 = -t286 * pkin(4) - t287 * pkin(9);
t328 = t330 ^ 2;
t207 = -t328 * pkin(4) + t329 * pkin(9) + t286 * t267 + t218;
t324 = t340 * g(1) - t345 * g(2);
t358 = -qJDD(1) * pkin(1) - t324;
t275 = -t319 * pkin(2) + qJDD(3) + t321 * t367 + (-qJ(3) * t333 - pkin(7)) * t346 + t358;
t240 = -t295 * pkin(3) - t306 * pkin(8) + t308 * t301 + t275;
t254 = t286 * qJD(4) + t338 * t295 + t343 * t296;
t215 = (-t286 * t330 - t254) * pkin(9) + (t287 * t330 - t253) * pkin(4) + t240;
t337 = sin(qJ(5));
t342 = cos(qJ(5));
t202 = -t337 * t207 + t342 * t215;
t276 = -t337 * t287 + t342 * t330;
t232 = t276 * qJD(5) + t342 * t254 + t337 * t329;
t252 = qJDD(5) - t253;
t277 = t342 * t287 + t337 * t330;
t282 = qJD(5) - t286;
t200 = (t276 * t282 - t232) * pkin(10) + (t276 * t277 + t252) * pkin(5) + t202;
t203 = t342 * t207 + t337 * t215;
t231 = -t277 * qJD(5) - t337 * t254 + t342 * t329;
t260 = t282 * pkin(5) - t277 * pkin(10);
t274 = t276 ^ 2;
t201 = -t274 * pkin(5) + t231 * pkin(10) - t282 * t260 + t203;
t336 = sin(qJ(6));
t341 = cos(qJ(6));
t198 = t341 * t200 - t336 * t201;
t255 = t341 * t276 - t336 * t277;
t212 = t255 * qJD(6) + t336 * t231 + t341 * t232;
t256 = t336 * t276 + t341 * t277;
t226 = -t255 * mrSges(7,1) + t256 * mrSges(7,2);
t280 = qJD(6) + t282;
t233 = -t280 * mrSges(7,2) + t255 * mrSges(7,3);
t247 = qJDD(6) + t252;
t193 = m(7) * t198 + t247 * mrSges(7,1) - t212 * mrSges(7,3) - t256 * t226 + t280 * t233;
t199 = t336 * t200 + t341 * t201;
t211 = -t256 * qJD(6) + t341 * t231 - t336 * t232;
t234 = t280 * mrSges(7,1) - t256 * mrSges(7,3);
t194 = m(7) * t199 - t247 * mrSges(7,2) + t211 * mrSges(7,3) + t255 * t226 - t280 * t234;
t185 = t341 * t193 + t336 * t194;
t257 = -t276 * mrSges(6,1) + t277 * mrSges(6,2);
t258 = -t282 * mrSges(6,2) + t276 * mrSges(6,3);
t183 = m(6) * t202 + t252 * mrSges(6,1) - t232 * mrSges(6,3) - t277 * t257 + t282 * t258 + t185;
t259 = t282 * mrSges(6,1) - t277 * mrSges(6,3);
t360 = -t336 * t193 + t341 * t194;
t184 = m(6) * t203 - t252 * mrSges(6,2) + t231 * mrSges(6,3) + t276 * t257 - t282 * t259 + t360;
t361 = -t337 * t183 + t342 * t184;
t176 = m(5) * t218 - t329 * mrSges(5,2) + t253 * mrSges(5,3) + t286 * t266 - t330 * t279 + t361;
t217 = t343 * t225 - t338 * t229;
t278 = -t330 * mrSges(5,2) + t286 * mrSges(5,3);
t206 = -t329 * pkin(4) - t328 * pkin(9) + t287 * t267 - t217;
t204 = -t231 * pkin(5) - t274 * pkin(10) + t277 * t260 + t206;
t355 = m(7) * t204 - t211 * mrSges(7,1) + t212 * mrSges(7,2) - t255 * t233 + t256 * t234;
t351 = -m(6) * t206 + t231 * mrSges(6,1) - t232 * mrSges(6,2) + t276 * t258 - t277 * t259 - t355;
t189 = m(5) * t217 + t329 * mrSges(5,1) - t254 * mrSges(5,3) - t287 * t266 + t330 * t278 + t351;
t167 = t338 * t176 + t343 * t189;
t290 = -t307 * mrSges(4,1) + t308 * mrSges(4,2);
t299 = -qJD(2) * mrSges(4,2) + t307 * mrSges(4,3);
t164 = m(4) * t242 + qJDD(2) * mrSges(4,1) - t296 * mrSges(4,3) + qJD(2) * t299 - t308 * t290 + t167;
t300 = qJD(2) * mrSges(4,1) - t308 * mrSges(4,3);
t362 = t343 * t176 - t338 * t189;
t165 = m(4) * t243 - qJDD(2) * mrSges(4,2) + t295 * mrSges(4,3) - qJD(2) * t300 + t307 * t290 + t362;
t159 = t335 * t164 + t334 * t165;
t297 = -t344 * g(3) - t368;
t310 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t339 + Ifges(3,2) * t344) * qJD(1);
t311 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t339 + Ifges(3,4) * t344) * qJD(1);
t284 = Ifges(4,4) * t308 + Ifges(4,2) * t307 + Ifges(4,6) * qJD(2);
t285 = Ifges(4,1) * t308 + Ifges(4,4) * t307 + Ifges(4,5) * qJD(2);
t220 = Ifges(7,5) * t256 + Ifges(7,6) * t255 + Ifges(7,3) * t280;
t222 = Ifges(7,1) * t256 + Ifges(7,4) * t255 + Ifges(7,5) * t280;
t186 = -mrSges(7,1) * t204 + mrSges(7,3) * t199 + Ifges(7,4) * t212 + Ifges(7,2) * t211 + Ifges(7,6) * t247 - t256 * t220 + t280 * t222;
t221 = Ifges(7,4) * t256 + Ifges(7,2) * t255 + Ifges(7,6) * t280;
t187 = mrSges(7,2) * t204 - mrSges(7,3) * t198 + Ifges(7,1) * t212 + Ifges(7,4) * t211 + Ifges(7,5) * t247 + t255 * t220 - t280 * t221;
t235 = Ifges(6,5) * t277 + Ifges(6,6) * t276 + Ifges(6,3) * t282;
t237 = Ifges(6,1) * t277 + Ifges(6,4) * t276 + Ifges(6,5) * t282;
t169 = -mrSges(6,1) * t206 + mrSges(6,3) * t203 + Ifges(6,4) * t232 + Ifges(6,2) * t231 + Ifges(6,6) * t252 - pkin(5) * t355 + pkin(10) * t360 + t341 * t186 + t336 * t187 - t277 * t235 + t282 * t237;
t236 = Ifges(6,4) * t277 + Ifges(6,2) * t276 + Ifges(6,6) * t282;
t171 = mrSges(6,2) * t206 - mrSges(6,3) * t202 + Ifges(6,1) * t232 + Ifges(6,4) * t231 + Ifges(6,5) * t252 - pkin(10) * t185 - t336 * t186 + t341 * t187 + t276 * t235 - t282 * t236;
t262 = Ifges(5,4) * t287 + Ifges(5,2) * t286 + Ifges(5,6) * t330;
t263 = Ifges(5,1) * t287 + Ifges(5,4) * t286 + Ifges(5,5) * t330;
t353 = -mrSges(5,1) * t217 + mrSges(5,2) * t218 - Ifges(5,5) * t254 - Ifges(5,6) * t253 - Ifges(5,3) * t329 - pkin(4) * t351 - pkin(9) * t361 - t342 * t169 - t337 * t171 - t287 * t262 + t286 * t263;
t350 = -mrSges(4,1) * t242 + mrSges(4,2) * t243 - Ifges(4,5) * t296 - Ifges(4,6) * t295 - Ifges(4,3) * qJDD(2) - pkin(3) * t167 - t308 * t284 + t307 * t285 + t353;
t370 = mrSges(3,1) * t297 - mrSges(3,2) * t298 + Ifges(3,5) * t318 + Ifges(3,6) * t319 + Ifges(3,3) * qJDD(2) + pkin(2) * t159 + (t339 * t310 - t344 * t311) * qJD(1) - t350;
t178 = t342 * t183 + t337 * t184;
t366 = qJD(1) * t344;
t317 = (-mrSges(3,1) * t344 + mrSges(3,2) * t339) * qJD(1);
t323 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t366;
t157 = m(3) * t297 + qJDD(2) * mrSges(3,1) - t318 * mrSges(3,3) + qJD(2) * t323 - t317 * t367 + t159;
t322 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t367;
t363 = -t334 * t164 + t335 * t165;
t158 = m(3) * t298 - qJDD(2) * mrSges(3,2) + t319 * mrSges(3,3) - qJD(2) * t322 + t317 * t366 + t363;
t364 = -t339 * t157 + t344 * t158;
t357 = m(5) * t240 - t253 * mrSges(5,1) + t254 * mrSges(5,2) - t286 * t278 + t287 * t279 + t178;
t261 = Ifges(5,5) * t287 + Ifges(5,6) * t286 + Ifges(5,3) * t330;
t155 = mrSges(5,2) * t240 - mrSges(5,3) * t217 + Ifges(5,1) * t254 + Ifges(5,4) * t253 + Ifges(5,5) * t329 - pkin(9) * t178 - t337 * t169 + t342 * t171 + t286 * t261 - t330 * t262;
t354 = -mrSges(7,1) * t198 + mrSges(7,2) * t199 - Ifges(7,5) * t212 - Ifges(7,6) * t211 - Ifges(7,3) * t247 - t256 * t221 + t255 * t222;
t348 = mrSges(6,1) * t202 - mrSges(6,2) * t203 + Ifges(6,5) * t232 + Ifges(6,6) * t231 + Ifges(6,3) * t252 + pkin(5) * t185 + t277 * t236 - t276 * t237 - t354;
t160 = -mrSges(5,1) * t240 + mrSges(5,3) * t218 + Ifges(5,4) * t254 + Ifges(5,2) * t253 + Ifges(5,6) * t329 - pkin(4) * t178 - t287 * t261 + t330 * t263 - t348;
t283 = Ifges(4,5) * t308 + Ifges(4,6) * t307 + Ifges(4,3) * qJD(2);
t150 = -mrSges(4,1) * t275 + mrSges(4,3) * t243 + Ifges(4,4) * t296 + Ifges(4,2) * t295 + Ifges(4,6) * qJDD(2) - pkin(3) * t357 + pkin(8) * t362 + qJD(2) * t285 + t338 * t155 + t343 * t160 - t308 * t283;
t151 = mrSges(4,2) * t275 - mrSges(4,3) * t242 + Ifges(4,1) * t296 + Ifges(4,4) * t295 + Ifges(4,5) * qJDD(2) - pkin(8) * t167 - qJD(2) * t284 + t343 * t155 - t338 * t160 + t307 * t283;
t309 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t339 + Ifges(3,6) * t344) * qJD(1);
t312 = -t346 * pkin(7) + t358;
t352 = m(4) * t275 - t295 * mrSges(4,1) + t296 * mrSges(4,2) - t307 * t299 + t308 * t300 + t357;
t146 = -mrSges(3,1) * t312 + mrSges(3,3) * t298 + Ifges(3,4) * t318 + Ifges(3,2) * t319 + Ifges(3,6) * qJDD(2) - pkin(2) * t352 + qJ(3) * t363 + qJD(2) * t311 + t335 * t150 + t334 * t151 - t309 * t367;
t148 = mrSges(3,2) * t312 - mrSges(3,3) * t297 + Ifges(3,1) * t318 + Ifges(3,4) * t319 + Ifges(3,5) * qJDD(2) - qJ(3) * t159 - qJD(2) * t310 - t334 * t150 + t335 * t151 + t309 * t366;
t349 = -m(3) * t312 + t319 * mrSges(3,1) - t318 * mrSges(3,2) - t322 * t367 + t323 * t366 - t352;
t356 = mrSges(2,1) * t324 - mrSges(2,2) * t325 + Ifges(2,3) * qJDD(1) + pkin(1) * t349 + pkin(7) * t364 + t344 * t146 + t339 * t148;
t172 = m(2) * t324 + qJDD(1) * mrSges(2,1) - t346 * mrSges(2,2) + t349;
t154 = t344 * t157 + t339 * t158;
t152 = m(2) * t325 - t346 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t364;
t149 = mrSges(2,1) * g(3) + mrSges(2,3) * t325 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t154 - t370;
t144 = -mrSges(2,2) * g(3) - mrSges(2,3) * t324 + Ifges(2,5) * qJDD(1) - t346 * Ifges(2,6) - pkin(7) * t154 - t339 * t146 + t344 * t148;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t144 - t340 * t149 - pkin(6) * (t340 * t152 + t345 * t172), t144, t148, t151, t155, t171, t187; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t340 * t144 + t345 * t149 + pkin(6) * (t345 * t152 - t340 * t172), t149, t146, t150, t160, t169, t186; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, t370, -t350, -t353, t348, -t354;];
m_new  = t1;
