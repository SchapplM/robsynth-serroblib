% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 20:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 20:39:36
% EndTime: 2019-05-06 20:41:35
% DurationCPUTime: 70.41s
% Computational Cost: add. (1181879->396), mult. (3084446->517), div. (0->0), fcn. (2481821->14), ass. (0->164)
t350 = sin(pkin(6));
t356 = sin(qJ(2));
t361 = cos(qJ(2));
t382 = qJD(1) * qJD(2);
t335 = (qJDD(1) * t356 + t361 * t382) * t350;
t352 = cos(pkin(6));
t343 = qJDD(1) * t352 + qJDD(2);
t344 = qJD(1) * t352 + qJD(2);
t357 = sin(qJ(1));
t362 = cos(qJ(1));
t340 = t357 * g(1) - g(2) * t362;
t363 = qJD(1) ^ 2;
t390 = pkin(8) * t350;
t332 = qJDD(1) * pkin(1) + t363 * t390 + t340;
t341 = -g(1) * t362 - g(2) * t357;
t333 = -pkin(1) * t363 + qJDD(1) * t390 + t341;
t385 = t352 * t361;
t373 = t332 * t385 - t356 * t333;
t389 = t350 ^ 2 * t363;
t269 = t343 * pkin(2) - t335 * qJ(3) + (pkin(2) * t356 * t389 + (qJ(3) * qJD(1) * t344 - g(3)) * t350) * t361 + t373;
t386 = t352 * t356;
t388 = t350 * t356;
t297 = -g(3) * t388 + t332 * t386 + t361 * t333;
t384 = qJD(1) * t350;
t380 = t356 * t384;
t329 = pkin(2) * t344 - qJ(3) * t380;
t336 = (qJDD(1) * t361 - t356 * t382) * t350;
t381 = t361 ^ 2 * t389;
t272 = -pkin(2) * t381 + qJ(3) * t336 - t329 * t344 + t297;
t349 = sin(pkin(12));
t351 = cos(pkin(12));
t326 = (t349 * t361 + t351 * t356) * t384;
t250 = -0.2e1 * qJD(3) * t326 + t351 * t269 - t349 * t272;
t387 = t350 * t361;
t379 = t361 * t384;
t325 = -t349 * t380 + t351 * t379;
t251 = 0.2e1 * qJD(3) * t325 + t349 * t269 + t351 * t272;
t298 = -mrSges(4,1) * t325 + mrSges(4,2) * t326;
t306 = -t349 * t335 + t336 * t351;
t312 = mrSges(4,1) * t344 - mrSges(4,3) * t326;
t300 = -pkin(3) * t325 - pkin(9) * t326;
t342 = t344 ^ 2;
t237 = -pkin(3) * t342 + pkin(9) * t343 + t300 * t325 + t251;
t316 = -t352 * g(3) - t350 * t332;
t282 = -t336 * pkin(2) - qJ(3) * t381 + t329 * t380 + qJDD(3) + t316;
t307 = t335 * t351 + t336 * t349;
t254 = (-t325 * t344 - t307) * pkin(9) + (t326 * t344 - t306) * pkin(3) + t282;
t355 = sin(qJ(4));
t360 = cos(qJ(4));
t229 = -t355 * t237 + t360 * t254;
t309 = -t326 * t355 + t344 * t360;
t280 = qJD(4) * t309 + t307 * t360 + t343 * t355;
t305 = qJDD(4) - t306;
t310 = t326 * t360 + t344 * t355;
t324 = qJD(4) - t325;
t226 = (t309 * t324 - t280) * pkin(10) + (t309 * t310 + t305) * pkin(4) + t229;
t230 = t360 * t237 + t355 * t254;
t279 = -qJD(4) * t310 - t307 * t355 + t343 * t360;
t291 = pkin(4) * t324 - pkin(10) * t310;
t308 = t309 ^ 2;
t228 = -pkin(4) * t308 + pkin(10) * t279 - t291 * t324 + t230;
t354 = sin(qJ(5));
t359 = cos(qJ(5));
t223 = t354 * t226 + t359 * t228;
t285 = t309 * t354 + t310 * t359;
t246 = -qJD(5) * t285 + t279 * t359 - t280 * t354;
t284 = t309 * t359 - t310 * t354;
t262 = -mrSges(6,1) * t284 + mrSges(6,2) * t285;
t319 = qJD(5) + t324;
t271 = mrSges(6,1) * t319 - mrSges(6,3) * t285;
t301 = qJDD(5) + t305;
t263 = -pkin(5) * t284 - pkin(11) * t285;
t318 = t319 ^ 2;
t220 = -pkin(5) * t318 + pkin(11) * t301 + t263 * t284 + t223;
t236 = -t343 * pkin(3) - t342 * pkin(9) + t326 * t300 - t250;
t231 = -t279 * pkin(4) - t308 * pkin(10) + t310 * t291 + t236;
t247 = qJD(5) * t284 + t279 * t354 + t280 * t359;
t224 = (-t284 * t319 - t247) * pkin(11) + (t285 * t319 - t246) * pkin(5) + t231;
t353 = sin(qJ(6));
t358 = cos(qJ(6));
t217 = -t220 * t353 + t224 * t358;
t265 = -t285 * t353 + t319 * t358;
t234 = qJD(6) * t265 + t247 * t358 + t301 * t353;
t245 = qJDD(6) - t246;
t266 = t285 * t358 + t319 * t353;
t255 = -mrSges(7,1) * t265 + mrSges(7,2) * t266;
t283 = qJD(6) - t284;
t256 = -mrSges(7,2) * t283 + mrSges(7,3) * t265;
t213 = m(7) * t217 + mrSges(7,1) * t245 - mrSges(7,3) * t234 - t255 * t266 + t256 * t283;
t218 = t220 * t358 + t224 * t353;
t233 = -qJD(6) * t266 - t247 * t353 + t301 * t358;
t257 = mrSges(7,1) * t283 - mrSges(7,3) * t266;
t214 = m(7) * t218 - mrSges(7,2) * t245 + mrSges(7,3) * t233 + t255 * t265 - t257 * t283;
t375 = -t213 * t353 + t358 * t214;
t199 = m(6) * t223 - mrSges(6,2) * t301 + mrSges(6,3) * t246 + t262 * t284 - t271 * t319 + t375;
t222 = t226 * t359 - t228 * t354;
t270 = -mrSges(6,2) * t319 + mrSges(6,3) * t284;
t219 = -pkin(5) * t301 - pkin(11) * t318 + t263 * t285 - t222;
t372 = -m(7) * t219 + t233 * mrSges(7,1) - mrSges(7,2) * t234 + t265 * t256 - t257 * t266;
t209 = m(6) * t222 + mrSges(6,1) * t301 - mrSges(6,3) * t247 - t262 * t285 + t270 * t319 + t372;
t194 = t354 * t199 + t359 * t209;
t287 = -mrSges(5,1) * t309 + mrSges(5,2) * t310;
t289 = -mrSges(5,2) * t324 + mrSges(5,3) * t309;
t192 = m(5) * t229 + mrSges(5,1) * t305 - mrSges(5,3) * t280 - t287 * t310 + t289 * t324 + t194;
t290 = mrSges(5,1) * t324 - mrSges(5,3) * t310;
t376 = t359 * t199 - t209 * t354;
t193 = m(5) * t230 - mrSges(5,2) * t305 + mrSges(5,3) * t279 + t287 * t309 - t290 * t324 + t376;
t377 = -t192 * t355 + t360 * t193;
t183 = m(4) * t251 - mrSges(4,2) * t343 + mrSges(4,3) * t306 + t298 * t325 - t312 * t344 + t377;
t311 = -mrSges(4,2) * t344 + mrSges(4,3) * t325;
t202 = t358 * t213 + t353 * t214;
t369 = m(6) * t231 - t246 * mrSges(6,1) + mrSges(6,2) * t247 - t284 * t270 + t271 * t285 + t202;
t365 = -m(5) * t236 + t279 * mrSges(5,1) - mrSges(5,2) * t280 + t309 * t289 - t290 * t310 - t369;
t196 = m(4) * t250 + mrSges(4,1) * t343 - mrSges(4,3) * t307 - t298 * t326 + t311 * t344 + t365;
t180 = t349 * t183 + t351 * t196;
t186 = t360 * t192 + t355 * t193;
t296 = -g(3) * t387 + t373;
t331 = -mrSges(3,2) * t344 + mrSges(3,3) * t379;
t334 = (-mrSges(3,1) * t361 + mrSges(3,2) * t356) * t384;
t178 = m(3) * t296 + mrSges(3,1) * t343 - mrSges(3,3) * t335 + t331 * t344 - t334 * t380 + t180;
t330 = mrSges(3,1) * t344 - mrSges(3,3) * t380;
t378 = t351 * t183 - t196 * t349;
t179 = m(3) * t297 - mrSges(3,2) * t343 + mrSges(3,3) * t336 - t330 * t344 + t334 * t379 + t378;
t169 = -t178 * t356 + t361 * t179;
t370 = m(4) * t282 - t306 * mrSges(4,1) + t307 * mrSges(4,2) - t325 * t311 + t326 * t312 + t186;
t184 = m(3) * t316 - t336 * mrSges(3,1) + t335 * mrSges(3,2) + (t330 * t356 - t331 * t361) * t384 + t370;
t166 = t178 * t385 + t179 * t386 - t184 * t350;
t238 = Ifges(7,5) * t266 + Ifges(7,6) * t265 + Ifges(7,3) * t283;
t240 = Ifges(7,1) * t266 + Ifges(7,4) * t265 + Ifges(7,5) * t283;
t206 = -mrSges(7,1) * t219 + mrSges(7,3) * t218 + Ifges(7,4) * t234 + Ifges(7,2) * t233 + Ifges(7,6) * t245 - t238 * t266 + t240 * t283;
t239 = Ifges(7,4) * t266 + Ifges(7,2) * t265 + Ifges(7,6) * t283;
t207 = mrSges(7,2) * t219 - mrSges(7,3) * t217 + Ifges(7,1) * t234 + Ifges(7,4) * t233 + Ifges(7,5) * t245 + t238 * t265 - t239 * t283;
t258 = Ifges(6,5) * t285 + Ifges(6,6) * t284 + Ifges(6,3) * t319;
t259 = Ifges(6,4) * t285 + Ifges(6,2) * t284 + Ifges(6,6) * t319;
t187 = mrSges(6,2) * t231 - mrSges(6,3) * t222 + Ifges(6,1) * t247 + Ifges(6,4) * t246 + Ifges(6,5) * t301 - pkin(11) * t202 - t206 * t353 + t207 * t358 + t258 * t284 - t259 * t319;
t260 = Ifges(6,1) * t285 + Ifges(6,4) * t284 + Ifges(6,5) * t319;
t366 = mrSges(7,1) * t217 - mrSges(7,2) * t218 + Ifges(7,5) * t234 + Ifges(7,6) * t233 + Ifges(7,3) * t245 + t239 * t266 - t240 * t265;
t188 = -mrSges(6,1) * t231 + mrSges(6,3) * t223 + Ifges(6,4) * t247 + Ifges(6,2) * t246 + Ifges(6,6) * t301 - pkin(5) * t202 - t258 * t285 + t260 * t319 - t366;
t273 = Ifges(5,5) * t310 + Ifges(5,6) * t309 + Ifges(5,3) * t324;
t275 = Ifges(5,1) * t310 + Ifges(5,4) * t309 + Ifges(5,5) * t324;
t172 = -mrSges(5,1) * t236 + mrSges(5,3) * t230 + Ifges(5,4) * t280 + Ifges(5,2) * t279 + Ifges(5,6) * t305 - pkin(4) * t369 + pkin(10) * t376 + t354 * t187 + t359 * t188 - t310 * t273 + t324 * t275;
t274 = Ifges(5,4) * t310 + Ifges(5,2) * t309 + Ifges(5,6) * t324;
t174 = mrSges(5,2) * t236 - mrSges(5,3) * t229 + Ifges(5,1) * t280 + Ifges(5,4) * t279 + Ifges(5,5) * t305 - pkin(10) * t194 + t187 * t359 - t188 * t354 + t273 * t309 - t274 * t324;
t292 = Ifges(4,5) * t326 + Ifges(4,6) * t325 + Ifges(4,3) * t344;
t293 = Ifges(4,4) * t326 + Ifges(4,2) * t325 + Ifges(4,6) * t344;
t162 = mrSges(4,2) * t282 - mrSges(4,3) * t250 + Ifges(4,1) * t307 + Ifges(4,4) * t306 + Ifges(4,5) * t343 - pkin(9) * t186 - t172 * t355 + t174 * t360 + t292 * t325 - t293 * t344;
t294 = Ifges(4,1) * t326 + Ifges(4,4) * t325 + Ifges(4,5) * t344;
t367 = -mrSges(6,1) * t222 + mrSges(6,2) * t223 - Ifges(6,5) * t247 - Ifges(6,6) * t246 - Ifges(6,3) * t301 - pkin(5) * t372 - pkin(11) * t375 - t358 * t206 - t353 * t207 - t285 * t259 + t284 * t260;
t364 = mrSges(5,1) * t229 - mrSges(5,2) * t230 + Ifges(5,5) * t280 + Ifges(5,6) * t279 + Ifges(5,3) * t305 + pkin(4) * t194 + t310 * t274 - t309 * t275 - t367;
t170 = -mrSges(4,1) * t282 + mrSges(4,3) * t251 + Ifges(4,4) * t307 + Ifges(4,2) * t306 + Ifges(4,6) * t343 - pkin(3) * t186 - t326 * t292 + t344 * t294 - t364;
t313 = Ifges(3,3) * t344 + (Ifges(3,5) * t356 + Ifges(3,6) * t361) * t384;
t315 = Ifges(3,5) * t344 + (Ifges(3,1) * t356 + Ifges(3,4) * t361) * t384;
t157 = -mrSges(3,1) * t316 + mrSges(3,3) * t297 + Ifges(3,4) * t335 + Ifges(3,2) * t336 + Ifges(3,6) * t343 - pkin(2) * t370 + qJ(3) * t378 + t349 * t162 + t351 * t170 - t313 * t380 + t344 * t315;
t314 = Ifges(3,6) * t344 + (Ifges(3,4) * t356 + Ifges(3,2) * t361) * t384;
t159 = mrSges(3,2) * t316 - mrSges(3,3) * t296 + Ifges(3,1) * t335 + Ifges(3,4) * t336 + Ifges(3,5) * t343 - qJ(3) * t180 + t162 * t351 - t170 * t349 + t313 * t379 - t314 * t344;
t368 = mrSges(4,1) * t250 - mrSges(4,2) * t251 + Ifges(4,5) * t307 + Ifges(4,6) * t306 + Ifges(4,3) * t343 + pkin(3) * t365 + pkin(9) * t377 + t360 * t172 + t355 * t174 + t326 * t293 - t325 * t294;
t161 = t368 + (t314 * t356 - t315 * t361) * t384 + pkin(2) * t180 + mrSges(3,1) * t296 - mrSges(3,2) * t297 + Ifges(3,5) * t335 + Ifges(3,6) * t336 + Ifges(3,3) * t343;
t371 = mrSges(2,1) * t340 - mrSges(2,2) * t341 + Ifges(2,3) * qJDD(1) + pkin(1) * t166 + t157 * t387 + t159 * t388 + t352 * t161 + t169 * t390;
t167 = m(2) * t341 - mrSges(2,1) * t363 - qJDD(1) * mrSges(2,2) + t169;
t165 = t352 * t184 + (t178 * t361 + t179 * t356) * t350;
t163 = m(2) * t340 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t363 + t166;
t155 = -mrSges(2,2) * g(3) - mrSges(2,3) * t340 + Ifges(2,5) * qJDD(1) - t363 * Ifges(2,6) - t356 * t157 + t361 * t159 + (-t165 * t350 - t166 * t352) * pkin(8);
t154 = mrSges(2,1) * g(3) + mrSges(2,3) * t341 + t363 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t165 - t350 * t161 + (pkin(8) * t169 + t157 * t361 + t159 * t356) * t352;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t362 * t155 - t357 * t154 - pkin(7) * (t163 * t362 + t167 * t357), t155, t159, t162, t174, t187, t207; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t357 * t155 + t362 * t154 + pkin(7) * (-t163 * t357 + t167 * t362), t154, t157, t170, t172, t188, t206; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t371, t371, t161, t368, t364, -t367, t366;];
m_new  = t1;
