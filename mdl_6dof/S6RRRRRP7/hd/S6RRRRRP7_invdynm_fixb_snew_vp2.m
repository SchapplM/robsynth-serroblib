% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP7
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
% Datum: 2019-05-08 05:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:20:41
% EndTime: 2019-05-08 05:22:30
% DurationCPUTime: 40.38s
% Computational Cost: add. (752426->394), mult. (1598498->497), div. (0->0), fcn. (1290945->12), ass. (0->157)
t339 = cos(pkin(6));
t334 = qJD(1) * t339 + qJD(2);
t342 = sin(qJ(3));
t347 = cos(qJ(3));
t343 = sin(qJ(2));
t338 = sin(pkin(6));
t371 = qJD(1) * t338;
t366 = t343 * t371;
t312 = t334 * t342 + t347 * t366;
t348 = cos(qJ(2));
t370 = qJD(1) * t348;
t323 = (qJD(2) * t370 + qJDD(1) * t343) * t338;
t333 = qJDD(1) * t339 + qJDD(2);
t292 = -qJD(3) * t312 - t323 * t342 + t333 * t347;
t311 = t334 * t347 - t342 * t366;
t293 = qJD(3) * t311 + t323 * t347 + t333 * t342;
t341 = sin(qJ(4));
t346 = cos(qJ(4));
t297 = t311 * t346 - t312 * t341;
t259 = qJD(4) * t297 + t292 * t341 + t293 * t346;
t298 = t311 * t341 + t312 * t346;
t365 = t338 * t370;
t328 = qJD(3) - t365;
t326 = qJD(4) + t328;
t340 = sin(qJ(5));
t345 = cos(qJ(5));
t282 = -t298 * t340 + t326 * t345;
t369 = qJDD(1) * t338;
t324 = -qJD(2) * t366 + t348 * t369;
t316 = qJDD(3) - t324;
t315 = qJDD(4) + t316;
t233 = qJD(5) * t282 + t259 * t345 + t315 * t340;
t283 = t298 * t345 + t326 * t340;
t262 = -mrSges(7,1) * t282 + mrSges(7,2) * t283;
t322 = (-pkin(2) * t348 - pkin(9) * t343) * t371;
t332 = t334 ^ 2;
t344 = sin(qJ(1));
t349 = cos(qJ(1));
t329 = g(1) * t344 - g(2) * t349;
t350 = qJD(1) ^ 2;
t381 = pkin(8) * t338;
t319 = qJDD(1) * pkin(1) + t350 * t381 + t329;
t330 = -g(1) * t349 - g(2) * t344;
t320 = -pkin(1) * t350 + pkin(8) * t369 + t330;
t376 = t339 * t343;
t372 = t319 * t376 + t320 * t348;
t279 = -t332 * pkin(2) + t333 * pkin(9) + (-g(3) * t343 + t322 * t370) * t338 + t372;
t380 = t339 * g(3);
t280 = -t324 * pkin(2) - t323 * pkin(9) - t380 + (-t319 + (pkin(2) * t343 - pkin(9) * t348) * t334 * qJD(1)) * t338;
t239 = -t342 * t279 + t280 * t347;
t219 = (t311 * t328 - t293) * pkin(10) + (t311 * t312 + t316) * pkin(3) + t239;
t240 = t279 * t347 + t280 * t342;
t302 = pkin(3) * t328 - pkin(10) * t312;
t310 = t311 ^ 2;
t222 = -pkin(3) * t310 + pkin(10) * t292 - t302 * t328 + t240;
t217 = t219 * t341 + t222 * t346;
t274 = -pkin(4) * t297 - pkin(11) * t298;
t325 = t326 ^ 2;
t211 = -pkin(4) * t325 + pkin(11) * t315 + t274 * t297 + t217;
t375 = t339 * t348;
t377 = t338 * t348;
t294 = -g(3) * t377 + t319 * t375 - t320 * t343;
t278 = -t333 * pkin(2) - t332 * pkin(9) + t322 * t366 - t294;
t237 = -t292 * pkin(3) - t310 * pkin(10) + t302 * t312 + t278;
t258 = -qJD(4) * t298 + t292 * t346 - t293 * t341;
t214 = (-t297 * t326 - t259) * pkin(11) + (t298 * t326 - t258) * pkin(4) + t237;
t205 = -t340 * t211 + t214 * t345;
t257 = qJDD(5) - t258;
t296 = qJD(5) - t297;
t201 = -0.2e1 * qJD(6) * t283 + (t282 * t296 - t233) * qJ(6) + (t282 * t283 + t257) * pkin(5) + t205;
t264 = -mrSges(7,2) * t296 + mrSges(7,3) * t282;
t368 = m(7) * t201 + mrSges(7,1) * t257 + t264 * t296;
t198 = -t233 * mrSges(7,3) - t283 * t262 + t368;
t206 = t211 * t345 + t214 * t340;
t232 = -qJD(5) * t283 - t259 * t340 + t315 * t345;
t244 = Ifges(6,4) * t283 + Ifges(6,2) * t282 + Ifges(6,6) * t296;
t245 = Ifges(7,1) * t283 + Ifges(7,4) * t282 + Ifges(7,5) * t296;
t246 = Ifges(6,1) * t283 + Ifges(6,4) * t282 + Ifges(6,5) * t296;
t266 = pkin(5) * t296 - qJ(6) * t283;
t281 = t282 ^ 2;
t204 = -pkin(5) * t281 + qJ(6) * t232 + 0.2e1 * qJD(6) * t282 - t266 * t296 + t206;
t243 = Ifges(7,4) * t283 + Ifges(7,2) * t282 + Ifges(7,6) * t296;
t359 = -mrSges(7,1) * t201 + mrSges(7,2) * t204 - Ifges(7,5) * t233 - Ifges(7,6) * t232 - Ifges(7,3) * t257 - t243 * t283;
t382 = mrSges(6,1) * t205 - mrSges(6,2) * t206 + Ifges(6,5) * t233 + Ifges(6,6) * t232 + Ifges(6,3) * t257 + pkin(5) * t198 + t283 * t244 - (t246 + t245) * t282 - t359;
t379 = -mrSges(6,2) - mrSges(7,2);
t378 = t338 * t343;
t273 = -mrSges(5,1) * t297 + mrSges(5,2) * t298;
t285 = mrSges(5,1) * t326 - mrSges(5,3) * t298;
t263 = -mrSges(6,1) * t282 + mrSges(6,2) * t283;
t265 = -mrSges(6,2) * t296 + mrSges(6,3) * t282;
t190 = m(6) * t205 + t257 * mrSges(6,1) + t296 * t265 + (-t262 - t263) * t283 + (-mrSges(6,3) - mrSges(7,3)) * t233 + t368;
t367 = m(7) * t204 + mrSges(7,3) * t232 + t262 * t282;
t267 = mrSges(7,1) * t296 - mrSges(7,3) * t283;
t373 = -mrSges(6,1) * t296 + mrSges(6,3) * t283 - t267;
t193 = m(6) * t206 + t232 * mrSges(6,3) + t257 * t379 + t282 * t263 + t296 * t373 + t367;
t362 = -t190 * t340 + t193 * t345;
t183 = m(5) * t217 - mrSges(5,2) * t315 + mrSges(5,3) * t258 + t273 * t297 - t285 * t326 + t362;
t216 = t219 * t346 - t222 * t341;
t284 = -mrSges(5,2) * t326 + mrSges(5,3) * t297;
t210 = -pkin(4) * t315 - pkin(11) * t325 + t274 * t298 - t216;
t208 = -pkin(5) * t232 - qJ(6) * t281 + t266 * t283 + qJDD(6) + t210;
t361 = -m(7) * t208 + mrSges(7,1) * t232 + t264 * t282;
t354 = -m(6) * t210 + mrSges(6,1) * t232 + t233 * t379 + t265 * t282 + t283 * t373 + t361;
t195 = m(5) * t216 + t315 * mrSges(5,1) - t259 * mrSges(5,3) - t298 * t273 + t326 * t284 + t354;
t176 = t183 * t341 + t195 * t346;
t299 = -mrSges(4,1) * t311 + mrSges(4,2) * t312;
t300 = -mrSges(4,2) * t328 + mrSges(4,3) * t311;
t174 = m(4) * t239 + mrSges(4,1) * t316 - mrSges(4,3) * t293 - t299 * t312 + t300 * t328 + t176;
t301 = mrSges(4,1) * t328 - mrSges(4,3) * t312;
t363 = t183 * t346 - t195 * t341;
t175 = m(4) * t240 - mrSges(4,2) * t316 + mrSges(4,3) * t292 + t299 * t311 - t301 * t328 + t363;
t169 = t174 * t347 + t175 * t342;
t187 = t190 * t345 + t193 * t340;
t295 = -g(3) * t378 + t372;
t317 = mrSges(3,1) * t334 - mrSges(3,3) * t366;
t321 = (-mrSges(3,1) * t348 + mrSges(3,2) * t343) * t371;
t364 = -t174 * t342 + t175 * t347;
t167 = m(3) * t295 - mrSges(3,2) * t333 + mrSges(3,3) * t324 - t317 * t334 + t321 * t365 + t364;
t318 = -mrSges(3,2) * t334 + mrSges(3,3) * t365;
t356 = m(5) * t237 - mrSges(5,1) * t258 + mrSges(5,2) * t259 - t284 * t297 + t285 * t298 + t187;
t353 = -m(4) * t278 + mrSges(4,1) * t292 - mrSges(4,2) * t293 + t300 * t311 - t301 * t312 - t356;
t180 = m(3) * t294 + mrSges(3,1) * t333 - mrSges(3,3) * t323 + t318 * t334 - t321 * t366 + t353;
t163 = t167 * t348 - t180 * t343;
t306 = -t338 * t319 - t380;
t168 = m(3) * t306 - t324 * mrSges(3,1) + t323 * mrSges(3,2) + (t317 * t343 - t318 * t348) * t371 + t169;
t159 = t167 * t376 - t168 * t338 + t180 * t375;
t360 = -mrSges(7,1) * t208 + mrSges(7,3) * t204 + Ifges(7,4) * t233 + Ifges(7,2) * t232 + Ifges(7,6) * t257 + t245 * t296;
t241 = Ifges(7,5) * t283 + Ifges(7,6) * t282 + Ifges(7,3) * t296;
t358 = mrSges(7,2) * t208 - mrSges(7,3) * t201 + Ifges(7,1) * t233 + Ifges(7,4) * t232 + Ifges(7,5) * t257 + t241 * t282;
t242 = Ifges(6,5) * t283 + Ifges(6,6) * t282 + Ifges(6,3) * t296;
t178 = Ifges(6,4) * t233 + Ifges(6,2) * t232 + Ifges(6,6) * t257 + t296 * t246 - mrSges(6,1) * t210 + mrSges(6,3) * t206 - pkin(5) * (t233 * mrSges(7,2) - t361) + qJ(6) * (-t257 * mrSges(7,2) - t296 * t267 + t367) + (-pkin(5) * t267 - t241 - t242) * t283 + t360;
t185 = mrSges(6,2) * t210 - mrSges(6,3) * t205 + Ifges(6,1) * t233 + Ifges(6,4) * t232 + Ifges(6,5) * t257 - qJ(6) * t198 + t282 * t242 + (-t243 - t244) * t296 + t358;
t269 = Ifges(5,5) * t298 + Ifges(5,6) * t297 + Ifges(5,3) * t326;
t270 = Ifges(5,4) * t298 + Ifges(5,2) * t297 + Ifges(5,6) * t326;
t164 = mrSges(5,2) * t237 - mrSges(5,3) * t216 + Ifges(5,1) * t259 + Ifges(5,4) * t258 + Ifges(5,5) * t315 - pkin(11) * t187 - t178 * t340 + t185 * t345 + t269 * t297 - t270 * t326;
t271 = Ifges(5,1) * t298 + Ifges(5,4) * t297 + Ifges(5,5) * t326;
t170 = -mrSges(5,1) * t237 + mrSges(5,3) * t217 + Ifges(5,4) * t259 + Ifges(5,2) * t258 + Ifges(5,6) * t315 - pkin(4) * t187 - t298 * t269 + t326 * t271 - t382;
t286 = Ifges(4,5) * t312 + Ifges(4,6) * t311 + Ifges(4,3) * t328;
t288 = Ifges(4,1) * t312 + Ifges(4,4) * t311 + Ifges(4,5) * t328;
t155 = -mrSges(4,1) * t278 + mrSges(4,3) * t240 + Ifges(4,4) * t293 + Ifges(4,2) * t292 + Ifges(4,6) * t316 - pkin(3) * t356 + pkin(10) * t363 + t341 * t164 + t346 * t170 - t312 * t286 + t328 * t288;
t287 = Ifges(4,4) * t312 + Ifges(4,2) * t311 + Ifges(4,6) * t328;
t160 = mrSges(4,2) * t278 - mrSges(4,3) * t239 + Ifges(4,1) * t293 + Ifges(4,4) * t292 + Ifges(4,5) * t316 - pkin(10) * t176 + t164 * t346 - t170 * t341 + t286 * t311 - t287 * t328;
t304 = Ifges(3,6) * t334 + (Ifges(3,4) * t343 + Ifges(3,2) * t348) * t371;
t305 = Ifges(3,5) * t334 + (Ifges(3,1) * t343 + Ifges(3,4) * t348) * t371;
t150 = Ifges(3,5) * t323 + Ifges(3,6) * t324 + Ifges(3,3) * t333 + mrSges(3,1) * t294 - mrSges(3,2) * t295 + t342 * t160 + t347 * t155 + pkin(2) * t353 + pkin(9) * t364 + (t304 * t343 - t305 * t348) * t371;
t303 = Ifges(3,3) * t334 + (Ifges(3,5) * t343 + Ifges(3,6) * t348) * t371;
t152 = mrSges(3,2) * t306 - mrSges(3,3) * t294 + Ifges(3,1) * t323 + Ifges(3,4) * t324 + Ifges(3,5) * t333 - pkin(9) * t169 - t155 * t342 + t160 * t347 + t303 * t365 - t304 * t334;
t355 = -mrSges(5,1) * t216 + mrSges(5,2) * t217 - Ifges(5,5) * t259 - Ifges(5,6) * t258 - Ifges(5,3) * t315 - pkin(4) * t354 - pkin(11) * t362 - t178 * t345 - t185 * t340 - t270 * t298 + t271 * t297;
t351 = mrSges(4,1) * t239 - mrSges(4,2) * t240 + Ifges(4,5) * t293 + Ifges(4,6) * t292 + Ifges(4,3) * t316 + pkin(3) * t176 + t287 * t312 - t288 * t311 - t355;
t154 = -mrSges(3,1) * t306 + mrSges(3,3) * t295 + Ifges(3,4) * t323 + Ifges(3,2) * t324 + Ifges(3,6) * t333 - pkin(2) * t169 - t303 * t366 + t305 * t334 - t351;
t357 = mrSges(2,1) * t329 - mrSges(2,2) * t330 + Ifges(2,3) * qJDD(1) + pkin(1) * t159 + t150 * t339 + t152 * t378 + t154 * t377 + t163 * t381;
t161 = m(2) * t330 - mrSges(2,1) * t350 - qJDD(1) * mrSges(2,2) + t163;
t158 = t339 * t168 + (t167 * t343 + t180 * t348) * t338;
t156 = m(2) * t329 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t350 + t159;
t148 = -mrSges(2,2) * g(3) - mrSges(2,3) * t329 + Ifges(2,5) * qJDD(1) - t350 * Ifges(2,6) + t348 * t152 - t343 * t154 + (-t158 * t338 - t159 * t339) * pkin(8);
t147 = mrSges(2,1) * g(3) + mrSges(2,3) * t330 + t350 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t158 - t338 * t150 + (pkin(8) * t163 + t152 * t343 + t154 * t348) * t339;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t349 * t148 - t344 * t147 - pkin(7) * (t156 * t349 + t161 * t344), t148, t152, t160, t164, t185, -t243 * t296 + t358; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t344 * t148 + t349 * t147 + pkin(7) * (-t156 * t344 + t161 * t349), t147, t154, t155, t170, t178, -t283 * t241 + t360; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t357, t357, t150, t351, -t355, t382, -t282 * t245 - t359;];
m_new  = t1;
