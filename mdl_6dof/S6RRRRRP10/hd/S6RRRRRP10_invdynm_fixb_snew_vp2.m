% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 06:39
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:18:33
% EndTime: 2019-05-08 06:20:20
% DurationCPUTime: 40.58s
% Computational Cost: add. (752782->395), mult. (1597166->498), div. (0->0), fcn. (1277567->12), ass. (0->156)
t336 = sin(pkin(6));
t341 = sin(qJ(2));
t345 = cos(qJ(2));
t364 = qJD(1) * qJD(2);
t323 = (-qJDD(1) * t345 + t341 * t364) * t336;
t377 = cos(qJ(5));
t376 = pkin(8) * t336;
t337 = cos(pkin(6));
t375 = t337 * g(3);
t374 = -mrSges(6,3) - mrSges(7,2);
t373 = t336 * t341;
t372 = t336 * t345;
t371 = t337 * t341;
t370 = t337 * t345;
t366 = qJD(1) * t336;
t321 = (-pkin(2) * t345 - pkin(9) * t341) * t366;
t332 = t337 * qJD(1) + qJD(2);
t330 = t332 ^ 2;
t331 = t337 * qJDD(1) + qJDD(2);
t365 = qJD(1) * t345;
t342 = sin(qJ(1));
t346 = cos(qJ(1));
t328 = t342 * g(1) - t346 * g(2);
t347 = qJD(1) ^ 2;
t318 = qJDD(1) * pkin(1) + t347 * t376 + t328;
t329 = -t346 * g(1) - t342 * g(2);
t319 = -t347 * pkin(1) + qJDD(1) * t376 + t329;
t367 = t318 * t371 + t345 * t319;
t267 = -t330 * pkin(2) + t331 * pkin(9) + (-g(3) * t341 + t321 * t365) * t336 + t367;
t322 = (qJDD(1) * t341 + t345 * t364) * t336;
t268 = t323 * pkin(2) - t322 * pkin(9) - t375 + (-t318 + (pkin(2) * t341 - pkin(9) * t345) * t332 * qJD(1)) * t336;
t340 = sin(qJ(3));
t344 = cos(qJ(3));
t235 = t344 * t267 + t340 * t268;
t362 = t341 * t366;
t310 = t344 * t332 - t340 * t362;
t311 = t340 * t332 + t344 * t362;
t293 = -t310 * pkin(3) - t311 * pkin(10);
t315 = qJDD(3) + t323;
t361 = t336 * t365;
t327 = qJD(3) - t361;
t325 = t327 ^ 2;
t227 = -t325 * pkin(3) + t315 * pkin(10) + t310 * t293 + t235;
t290 = -g(3) * t372 + t318 * t370 - t341 * t319;
t266 = -t331 * pkin(2) - t330 * pkin(9) + t321 * t362 - t290;
t288 = -t311 * qJD(3) - t340 * t322 + t344 * t331;
t289 = t310 * qJD(3) + t344 * t322 + t340 * t331;
t232 = (-t310 * t327 - t289) * pkin(10) + (t311 * t327 - t288) * pkin(3) + t266;
t339 = sin(qJ(4));
t343 = cos(qJ(4));
t211 = -t339 * t227 + t343 * t232;
t296 = -t339 * t311 + t343 * t327;
t253 = t296 * qJD(4) + t343 * t289 + t339 * t315;
t286 = qJDD(4) - t288;
t297 = t343 * t311 + t339 * t327;
t309 = qJD(4) - t310;
t208 = (t296 * t309 - t253) * pkin(11) + (t296 * t297 + t286) * pkin(4) + t211;
t212 = t343 * t227 + t339 * t232;
t252 = -t297 * qJD(4) - t339 * t289 + t343 * t315;
t275 = t309 * pkin(4) - t297 * pkin(11);
t295 = t296 ^ 2;
t210 = -t295 * pkin(4) + t252 * pkin(11) - t309 * t275 + t212;
t338 = sin(qJ(5));
t204 = t338 * t208 + t377 * t210;
t270 = t338 * t296 + t377 * t297;
t223 = t270 * qJD(5) - t377 * t252 + t338 * t253;
t307 = qJD(5) + t309;
t256 = t307 * mrSges(6,1) - t270 * mrSges(6,3);
t269 = -t377 * t296 + t338 * t297;
t281 = qJDD(5) + t286;
t245 = t269 * pkin(5) - t270 * qJ(6);
t305 = t307 ^ 2;
t199 = -t305 * pkin(5) + t281 * qJ(6) + 0.2e1 * qJD(6) * t307 - t269 * t245 + t204;
t257 = -t307 * mrSges(7,1) + t270 * mrSges(7,2);
t363 = m(7) * t199 + t281 * mrSges(7,3) + t307 * t257;
t246 = t269 * mrSges(7,1) - t270 * mrSges(7,3);
t368 = -t269 * mrSges(6,1) - t270 * mrSges(6,2) - t246;
t186 = m(6) * t204 - t281 * mrSges(6,2) + t374 * t223 - t307 * t256 + t368 * t269 + t363;
t203 = t377 * t208 - t338 * t210;
t224 = -t269 * qJD(5) + t338 * t252 + t377 * t253;
t255 = -t307 * mrSges(6,2) - t269 * mrSges(6,3);
t201 = -t281 * pkin(5) - t305 * qJ(6) + t270 * t245 + qJDD(6) - t203;
t254 = -t269 * mrSges(7,2) + t307 * mrSges(7,3);
t358 = -m(7) * t201 + t281 * mrSges(7,1) + t307 * t254;
t188 = m(6) * t203 + t281 * mrSges(6,1) + t374 * t224 + t307 * t255 + t368 * t270 + t358;
t183 = t338 * t186 + t377 * t188;
t271 = -t296 * mrSges(5,1) + t297 * mrSges(5,2);
t273 = -t309 * mrSges(5,2) + t296 * mrSges(5,3);
t179 = m(5) * t211 + t286 * mrSges(5,1) - t253 * mrSges(5,3) - t297 * t271 + t309 * t273 + t183;
t274 = t309 * mrSges(5,1) - t297 * mrSges(5,3);
t359 = t377 * t186 - t338 * t188;
t180 = m(5) * t212 - t286 * mrSges(5,2) + t252 * mrSges(5,3) + t296 * t271 - t309 * t274 + t359;
t177 = -t339 * t179 + t343 * t180;
t292 = -t310 * mrSges(4,1) + t311 * mrSges(4,2);
t299 = t327 * mrSges(4,1) - t311 * mrSges(4,3);
t175 = m(4) * t235 - t315 * mrSges(4,2) + t288 * mrSges(4,3) + t310 * t292 - t327 * t299 + t177;
t234 = -t340 * t267 + t344 * t268;
t226 = -t315 * pkin(3) - t325 * pkin(10) + t311 * t293 - t234;
t213 = -t252 * pkin(4) - t295 * pkin(11) + t297 * t275 + t226;
t206 = -0.2e1 * qJD(6) * t270 + (t269 * t307 - t224) * qJ(6) + (t270 * t307 + t223) * pkin(5) + t213;
t196 = m(7) * t206 + t223 * mrSges(7,1) - t224 * mrSges(7,3) + t269 * t254 - t270 * t257;
t352 = m(6) * t213 + t223 * mrSges(6,1) + t224 * mrSges(6,2) + t269 * t255 + t270 * t256 + t196;
t191 = -m(5) * t226 + t252 * mrSges(5,1) - t253 * mrSges(5,2) + t296 * t273 - t297 * t274 - t352;
t298 = -t327 * mrSges(4,2) + t310 * mrSges(4,3);
t190 = m(4) * t234 + t315 * mrSges(4,1) - t289 * mrSges(4,3) - t311 * t292 + t327 * t298 + t191;
t170 = t340 * t175 + t344 * t190;
t238 = Ifges(7,4) * t270 + Ifges(7,2) * t307 + Ifges(7,6) * t269;
t369 = -Ifges(6,5) * t270 + Ifges(6,6) * t269 - Ifges(6,3) * t307 - t238;
t291 = -g(3) * t373 + t367;
t316 = t332 * mrSges(3,1) - mrSges(3,3) * t362;
t320 = (-mrSges(3,1) * t345 + mrSges(3,2) * t341) * t366;
t360 = t344 * t175 - t340 * t190;
t168 = m(3) * t291 - t331 * mrSges(3,2) - t323 * mrSges(3,3) - t332 * t316 + t320 * t361 + t360;
t317 = -t332 * mrSges(3,2) + mrSges(3,3) * t361;
t176 = t343 * t179 + t339 * t180;
t351 = -m(4) * t266 + t288 * mrSges(4,1) - t289 * mrSges(4,2) + t310 * t298 - t311 * t299 - t176;
t172 = m(3) * t290 + t331 * mrSges(3,1) - t322 * mrSges(3,3) + t332 * t317 - t320 * t362 + t351;
t162 = t345 * t168 - t341 * t172;
t303 = -t336 * t318 - t375;
t169 = m(3) * t303 + t323 * mrSges(3,1) + t322 * mrSges(3,2) + (t316 * t341 - t317 * t345) * t366 + t170;
t159 = t168 * t371 - t336 * t169 + t172 * t370;
t357 = -mrSges(7,1) * t206 + mrSges(7,2) * t199;
t236 = Ifges(7,5) * t270 + Ifges(7,6) * t307 + Ifges(7,3) * t269;
t355 = mrSges(7,2) * t201 - mrSges(7,3) * t206 + Ifges(7,1) * t224 + Ifges(7,4) * t281 + Ifges(7,5) * t223 + t307 * t236;
t240 = Ifges(7,1) * t270 + Ifges(7,4) * t307 + Ifges(7,5) * t269;
t241 = Ifges(6,1) * t270 - Ifges(6,4) * t269 + Ifges(6,5) * t307;
t181 = -mrSges(6,1) * t213 + mrSges(6,3) * t204 - pkin(5) * t196 + (t240 + t241) * t307 + (Ifges(6,6) - Ifges(7,6)) * t281 + t369 * t270 + (Ifges(6,4) - Ifges(7,5)) * t224 + (-Ifges(6,2) - Ifges(7,3)) * t223 + t357;
t239 = Ifges(6,4) * t270 - Ifges(6,2) * t269 + Ifges(6,6) * t307;
t182 = mrSges(6,2) * t213 - mrSges(6,3) * t203 + Ifges(6,1) * t224 - Ifges(6,4) * t223 + Ifges(6,5) * t281 - qJ(6) * t196 - t307 * t239 + t369 * t269 + t355;
t258 = Ifges(5,5) * t297 + Ifges(5,6) * t296 + Ifges(5,3) * t309;
t260 = Ifges(5,1) * t297 + Ifges(5,4) * t296 + Ifges(5,5) * t309;
t164 = -mrSges(5,1) * t226 + mrSges(5,3) * t212 + Ifges(5,4) * t253 + Ifges(5,2) * t252 + Ifges(5,6) * t286 - pkin(4) * t352 + pkin(11) * t359 + t377 * t181 + t338 * t182 - t297 * t258 + t309 * t260;
t259 = Ifges(5,4) * t297 + Ifges(5,2) * t296 + Ifges(5,6) * t309;
t165 = mrSges(5,2) * t226 - mrSges(5,3) * t211 + Ifges(5,1) * t253 + Ifges(5,4) * t252 + Ifges(5,5) * t286 - pkin(11) * t183 - t338 * t181 + t377 * t182 + t296 * t258 - t309 * t259;
t282 = Ifges(4,5) * t311 + Ifges(4,6) * t310 + Ifges(4,3) * t327;
t283 = Ifges(4,4) * t311 + Ifges(4,2) * t310 + Ifges(4,6) * t327;
t155 = mrSges(4,2) * t266 - mrSges(4,3) * t234 + Ifges(4,1) * t289 + Ifges(4,4) * t288 + Ifges(4,5) * t315 - pkin(10) * t176 - t339 * t164 + t343 * t165 + t310 * t282 - t327 * t283;
t284 = Ifges(4,1) * t311 + Ifges(4,4) * t310 + Ifges(4,5) * t327;
t353 = mrSges(7,1) * t201 - mrSges(7,3) * t199 - Ifges(7,4) * t224 - Ifges(7,2) * t281 - Ifges(7,6) * t223 + t270 * t236 - t269 * t240;
t350 = mrSges(6,2) * t204 - t269 * t241 - qJ(6) * (-t223 * mrSges(7,2) - t269 * t246 + t363) - pkin(5) * (-t224 * mrSges(7,2) - t270 * t246 + t358) - mrSges(6,1) * t203 + Ifges(6,6) * t223 - Ifges(6,5) * t224 - t270 * t239 - Ifges(6,3) * t281 + t353;
t348 = mrSges(5,1) * t211 - mrSges(5,2) * t212 + Ifges(5,5) * t253 + Ifges(5,6) * t252 + Ifges(5,3) * t286 + pkin(4) * t183 + t297 * t259 - t296 * t260 - t350;
t163 = -mrSges(4,1) * t266 + mrSges(4,3) * t235 + Ifges(4,4) * t289 + Ifges(4,2) * t288 + Ifges(4,6) * t315 - pkin(3) * t176 - t311 * t282 + t327 * t284 - t348;
t301 = Ifges(3,6) * t332 + (Ifges(3,4) * t341 + Ifges(3,2) * t345) * t366;
t302 = Ifges(3,5) * t332 + (Ifges(3,1) * t341 + Ifges(3,4) * t345) * t366;
t150 = Ifges(3,5) * t322 - Ifges(3,6) * t323 + Ifges(3,3) * t331 + mrSges(3,1) * t290 - mrSges(3,2) * t291 + t340 * t155 + t344 * t163 + pkin(2) * t351 + pkin(9) * t360 + (t301 * t341 - t302 * t345) * t366;
t300 = Ifges(3,3) * t332 + (Ifges(3,5) * t341 + Ifges(3,6) * t345) * t366;
t152 = mrSges(3,2) * t303 - mrSges(3,3) * t290 + Ifges(3,1) * t322 - Ifges(3,4) * t323 + Ifges(3,5) * t331 - pkin(9) * t170 + t344 * t155 - t340 * t163 + t300 * t361 - t332 * t301;
t349 = mrSges(4,1) * t234 - mrSges(4,2) * t235 + Ifges(4,5) * t289 + Ifges(4,6) * t288 + Ifges(4,3) * t315 + pkin(3) * t191 + pkin(10) * t177 + t343 * t164 + t339 * t165 + t311 * t283 - t310 * t284;
t154 = -mrSges(3,1) * t303 + mrSges(3,3) * t291 + Ifges(3,4) * t322 - Ifges(3,2) * t323 + Ifges(3,6) * t331 - pkin(2) * t170 - t300 * t362 + t332 * t302 - t349;
t354 = mrSges(2,1) * t328 - mrSges(2,2) * t329 + Ifges(2,3) * qJDD(1) + pkin(1) * t159 + t337 * t150 + t152 * t373 + t154 * t372 + t162 * t376;
t160 = m(2) * t329 - t347 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t162;
t158 = t337 * t169 + (t168 * t341 + t172 * t345) * t336;
t156 = m(2) * t328 + qJDD(1) * mrSges(2,1) - t347 * mrSges(2,2) + t159;
t148 = -mrSges(2,2) * g(3) - mrSges(2,3) * t328 + Ifges(2,5) * qJDD(1) - t347 * Ifges(2,6) + t345 * t152 - t341 * t154 + (-t158 * t336 - t159 * t337) * pkin(8);
t147 = mrSges(2,1) * g(3) + mrSges(2,3) * t329 + t347 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t158 - t336 * t150 + (pkin(8) * t162 + t152 * t341 + t154 * t345) * t337;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t346 * t148 - t342 * t147 - pkin(7) * (t346 * t156 + t342 * t160), t148, t152, t155, t165, t182, -t269 * t238 + t355; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t342 * t148 + t346 * t147 + pkin(7) * (-t342 * t156 + t346 * t160), t147, t154, t163, t164, t181, -t353; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t354, t354, t150, t349, t348, -t350, Ifges(7,5) * t224 + Ifges(7,6) * t281 + Ifges(7,3) * t223 + t270 * t238 - t307 * t240 - t357;];
m_new  = t1;
