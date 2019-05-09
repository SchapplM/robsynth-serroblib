% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:34:13
% EndTime: 2019-05-06 01:34:48
% DurationCPUTime: 18.74s
% Computational Cost: add. (303498->360), mult. (717094->439), div. (0->0), fcn. (544602->10), ass. (0->146)
t332 = sin(pkin(10));
t333 = cos(pkin(10));
t336 = sin(qJ(3));
t340 = cos(qJ(3));
t357 = t332 * t340 + t333 * t336;
t376 = qJD(1) * t333;
t377 = qJD(1) * t332;
t310 = -t336 * t377 + t340 * t376;
t375 = t310 * qJD(3);
t297 = qJDD(1) * t357 + t375;
t311 = t357 * qJD(1);
t335 = sin(qJ(4));
t339 = cos(qJ(4));
t302 = t335 * qJD(3) + t339 * t311;
t263 = -t302 * qJD(4) + t339 * qJDD(3) - t335 * t297;
t301 = t339 * qJD(3) - t335 * t311;
t264 = t301 * qJD(4) + t335 * qJDD(3) + t339 * t297;
t334 = sin(qJ(5));
t338 = cos(qJ(5));
t266 = t338 * t301 - t334 * t302;
t225 = t266 * qJD(5) + t334 * t263 + t338 * t264;
t267 = t334 * t301 + t338 * t302;
t246 = -t266 * mrSges(7,1) + t267 * mrSges(7,2);
t337 = sin(qJ(1));
t341 = cos(qJ(1));
t319 = -t341 * g(1) - t337 * g(2);
t343 = qJD(1) ^ 2;
t312 = -t343 * pkin(1) + qJDD(1) * qJ(2) + t319;
t374 = qJD(1) * qJD(2);
t369 = -t333 * g(3) - 0.2e1 * t332 * t374;
t381 = pkin(2) * t333;
t277 = (-pkin(7) * qJDD(1) + t343 * t381 - t312) * t332 + t369;
t299 = -t332 * g(3) + (t312 + 0.2e1 * t374) * t333;
t372 = qJDD(1) * t333;
t329 = t333 ^ 2;
t379 = t329 * t343;
t278 = -pkin(2) * t379 + pkin(7) * t372 + t299;
t251 = t336 * t277 + t340 * t278;
t294 = -t310 * pkin(3) - t311 * pkin(8);
t342 = qJD(3) ^ 2;
t233 = -t342 * pkin(3) + qJDD(3) * pkin(8) + t310 * t294 + t251;
t328 = t332 ^ 2;
t318 = t337 * g(1) - t341 * g(2);
t363 = qJDD(2) - t318;
t295 = (-pkin(1) - t381) * qJDD(1) + (-qJ(2) + (-t328 - t329) * pkin(7)) * t343 + t363;
t307 = t311 * qJD(3);
t373 = qJDD(1) * t332;
t296 = -t336 * t373 + t340 * t372 - t307;
t236 = (-t297 - t375) * pkin(8) + (-t296 + t307) * pkin(3) + t295;
t210 = -t335 * t233 + t339 * t236;
t293 = qJDD(4) - t296;
t308 = qJD(4) - t310;
t206 = (t301 * t308 - t264) * pkin(9) + (t301 * t302 + t293) * pkin(4) + t210;
t211 = t339 * t233 + t335 * t236;
t276 = t308 * pkin(4) - t302 * pkin(9);
t300 = t301 ^ 2;
t208 = -t300 * pkin(4) + t263 * pkin(9) - t308 * t276 + t211;
t199 = t338 * t206 - t334 * t208;
t288 = qJDD(5) + t293;
t306 = qJD(5) + t308;
t194 = -0.2e1 * qJD(6) * t267 + (t266 * t306 - t225) * qJ(6) + (t266 * t267 + t288) * pkin(5) + t199;
t252 = -t306 * mrSges(7,2) + t266 * mrSges(7,3);
t371 = m(7) * t194 + t288 * mrSges(7,1) + t306 * t252;
t191 = -t225 * mrSges(7,3) - t267 * t246 + t371;
t200 = t334 * t206 + t338 * t208;
t224 = -t267 * qJD(5) + t338 * t263 - t334 * t264;
t240 = Ifges(6,4) * t267 + Ifges(6,2) * t266 + Ifges(6,6) * t306;
t241 = Ifges(7,1) * t267 + Ifges(7,4) * t266 + Ifges(7,5) * t306;
t242 = Ifges(6,1) * t267 + Ifges(6,4) * t266 + Ifges(6,5) * t306;
t254 = t306 * pkin(5) - t267 * qJ(6);
t265 = t266 ^ 2;
t197 = -t265 * pkin(5) + t224 * qJ(6) + 0.2e1 * qJD(6) * t266 - t306 * t254 + t200;
t239 = Ifges(7,4) * t267 + Ifges(7,2) * t266 + Ifges(7,6) * t306;
t354 = -mrSges(7,1) * t194 + mrSges(7,2) * t197 - Ifges(7,5) * t225 - Ifges(7,6) * t224 - Ifges(7,3) * t288 - t267 * t239;
t384 = mrSges(6,1) * t199 - mrSges(6,2) * t200 + Ifges(6,5) * t225 + Ifges(6,6) * t224 + Ifges(6,3) * t288 + pkin(5) * t191 + t267 * t240 - t354 + (-t242 - t241) * t266;
t247 = -t266 * mrSges(6,1) + t267 * mrSges(6,2);
t253 = -t306 * mrSges(6,2) + t266 * mrSges(6,3);
t182 = m(6) * t199 + t288 * mrSges(6,1) + t306 * t253 + (-t246 - t247) * t267 + (-mrSges(6,3) - mrSges(7,3)) * t225 + t371;
t255 = t306 * mrSges(7,1) - t267 * mrSges(7,3);
t256 = t306 * mrSges(6,1) - t267 * mrSges(6,3);
t370 = m(7) * t197 + t224 * mrSges(7,3) + t266 * t246;
t185 = m(6) * t200 + t224 * mrSges(6,3) + t266 * t247 + (-t255 - t256) * t306 + (-mrSges(6,2) - mrSges(7,2)) * t288 + t370;
t180 = t338 * t182 + t334 * t185;
t258 = Ifges(5,4) * t302 + Ifges(5,2) * t301 + Ifges(5,6) * t308;
t259 = Ifges(5,1) * t302 + Ifges(5,4) * t301 + Ifges(5,5) * t308;
t383 = mrSges(5,1) * t210 - mrSges(5,2) * t211 + Ifges(5,5) * t264 + Ifges(5,6) * t263 + Ifges(5,3) * t293 + pkin(4) * t180 + t302 * t258 - t301 * t259 + t384;
t289 = -t310 * mrSges(4,1) + t311 * mrSges(4,2);
t304 = qJD(3) * mrSges(4,1) - t311 * mrSges(4,3);
t269 = -t301 * mrSges(5,1) + t302 * mrSges(5,2);
t272 = -t308 * mrSges(5,2) + t301 * mrSges(5,3);
t177 = m(5) * t210 + t293 * mrSges(5,1) - t264 * mrSges(5,3) - t302 * t269 + t308 * t272 + t180;
t273 = t308 * mrSges(5,1) - t302 * mrSges(5,3);
t364 = -t334 * t182 + t338 * t185;
t178 = m(5) * t211 - t293 * mrSges(5,2) + t263 * mrSges(5,3) + t301 * t269 - t308 * t273 + t364;
t365 = -t335 * t177 + t339 * t178;
t169 = m(4) * t251 - qJDD(3) * mrSges(4,2) + t296 * mrSges(4,3) - qJD(3) * t304 + t310 * t289 + t365;
t250 = t340 * t277 - t336 * t278;
t303 = -qJD(3) * mrSges(4,2) + t310 * mrSges(4,3);
t232 = -qJDD(3) * pkin(3) - t342 * pkin(8) + t311 * t294 - t250;
t209 = -t263 * pkin(4) - t300 * pkin(9) + t302 * t276 + t232;
t203 = -t224 * pkin(5) - t265 * qJ(6) + t267 * t254 + qJDD(6) + t209;
t359 = m(7) * t203 - t224 * mrSges(7,1) + t225 * mrSges(7,2) - t266 * t252 + t267 * t255;
t349 = m(6) * t209 - t224 * mrSges(6,1) + t225 * mrSges(6,2) - t266 * t253 + t267 * t256 + t359;
t345 = -m(5) * t232 + t263 * mrSges(5,1) - t264 * mrSges(5,2) + t301 * t272 - t302 * t273 - t349;
t187 = m(4) * t250 + qJDD(3) * mrSges(4,1) - t297 * mrSges(4,3) + qJD(3) * t303 - t311 * t289 + t345;
t164 = t336 * t169 + t340 * t187;
t298 = -t332 * t312 + t369;
t237 = Ifges(7,5) * t267 + Ifges(7,6) * t266 + Ifges(7,3) * t306;
t238 = Ifges(6,5) * t267 + Ifges(6,6) * t266 + Ifges(6,3) * t306;
t355 = -mrSges(7,1) * t203 + mrSges(7,3) * t197 + Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t288 + t306 * t241;
t173 = Ifges(6,4) * t225 + Ifges(6,2) * t224 + Ifges(6,6) * t288 + t306 * t242 - mrSges(6,1) * t209 + mrSges(6,3) * t200 - pkin(5) * t359 + qJ(6) * (-t288 * mrSges(7,2) - t306 * t255 + t370) + (-t238 - t237) * t267 + t355;
t353 = mrSges(7,2) * t203 - mrSges(7,3) * t194 + Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t288 + t266 * t237;
t179 = mrSges(6,2) * t209 - mrSges(6,3) * t199 + Ifges(6,1) * t225 + Ifges(6,4) * t224 + Ifges(6,5) * t288 - qJ(6) * t191 + t266 * t238 + (-t239 - t240) * t306 + t353;
t257 = Ifges(5,5) * t302 + Ifges(5,6) * t301 + Ifges(5,3) * t308;
t158 = -mrSges(5,1) * t232 + mrSges(5,3) * t211 + Ifges(5,4) * t264 + Ifges(5,2) * t263 + Ifges(5,6) * t293 - pkin(4) * t349 + pkin(9) * t364 + t338 * t173 + t334 * t179 - t302 * t257 + t308 * t259;
t160 = mrSges(5,2) * t232 - mrSges(5,3) * t210 + Ifges(5,1) * t264 + Ifges(5,4) * t263 + Ifges(5,5) * t293 - pkin(9) * t180 - t334 * t173 + t338 * t179 + t301 * t257 - t308 * t258;
t280 = Ifges(4,4) * t311 + Ifges(4,2) * t310 + Ifges(4,6) * qJD(3);
t281 = Ifges(4,1) * t311 + Ifges(4,4) * t310 + Ifges(4,5) * qJD(3);
t350 = -mrSges(4,1) * t250 + mrSges(4,2) * t251 - Ifges(4,5) * t297 - Ifges(4,6) * t296 - Ifges(4,3) * qJDD(3) - pkin(3) * t345 - pkin(8) * t365 - t339 * t158 - t335 * t160 - t311 * t280 + t310 * t281;
t361 = Ifges(3,4) * t332 + Ifges(3,2) * t333;
t362 = Ifges(3,1) * t332 + Ifges(3,4) * t333;
t382 = -mrSges(3,1) * t298 + mrSges(3,2) * t299 - pkin(2) * t164 - (t361 * t377 - t362 * t376) * qJD(1) + t350;
t380 = mrSges(3,2) * t332;
t171 = t339 * t177 + t335 * t178;
t356 = mrSges(3,3) * qJDD(1) + t343 * (-mrSges(3,1) * t333 + t380);
t162 = m(3) * t298 - t332 * t356 + t164;
t366 = t340 * t169 - t336 * t187;
t163 = m(3) * t299 + t333 * t356 + t366;
t367 = -t332 * t162 + t333 * t163;
t360 = Ifges(3,5) * t332 + Ifges(3,6) * t333;
t279 = Ifges(4,5) * t311 + Ifges(4,6) * t310 + Ifges(4,3) * qJD(3);
t152 = mrSges(4,2) * t295 - mrSges(4,3) * t250 + Ifges(4,1) * t297 + Ifges(4,4) * t296 + Ifges(4,5) * qJDD(3) - pkin(8) * t171 - qJD(3) * t280 - t335 * t158 + t339 * t160 + t310 * t279;
t156 = -mrSges(4,1) * t295 + mrSges(4,3) * t251 + Ifges(4,4) * t297 + Ifges(4,2) * t296 + Ifges(4,6) * qJDD(3) - pkin(3) * t171 + qJD(3) * t281 - t311 * t279 - t383;
t309 = -qJDD(1) * pkin(1) - t343 * qJ(2) + t363;
t314 = t360 * qJD(1);
t351 = m(4) * t295 - t296 * mrSges(4,1) + t297 * mrSges(4,2) - t310 * t303 + t311 * t304 + t171;
t148 = -mrSges(3,1) * t309 + mrSges(3,3) * t299 - pkin(2) * t351 + pkin(7) * t366 + qJDD(1) * t361 + t336 * t152 + t340 * t156 - t314 * t377;
t151 = mrSges(3,2) * t309 - mrSges(3,3) * t298 - pkin(7) * t164 + qJDD(1) * t362 + t340 * t152 - t336 * t156 + t314 * t376;
t347 = -m(3) * t309 + mrSges(3,1) * t372 - t351 + (t328 * t343 + t379) * mrSges(3,3);
t352 = -mrSges(2,2) * t319 + qJ(2) * t367 + t333 * t148 + t332 * t151 + pkin(1) * (-mrSges(3,2) * t373 + t347) + mrSges(2,1) * t318 + Ifges(2,3) * qJDD(1);
t165 = -t343 * mrSges(2,2) + m(2) * t318 + t347 + (mrSges(2,1) - t380) * qJDD(1);
t155 = t333 * t162 + t332 * t163;
t153 = m(2) * t319 - t343 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t367;
t149 = (Ifges(2,6) - t360) * qJDD(1) + mrSges(2,1) * g(3) + t343 * Ifges(2,5) + mrSges(2,3) * t319 - pkin(1) * t155 + t382;
t146 = -mrSges(2,2) * g(3) - mrSges(2,3) * t318 + Ifges(2,5) * qJDD(1) - t343 * Ifges(2,6) - qJ(2) * t155 - t332 * t148 + t333 * t151;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t341 * t146 - t337 * t149 - pkin(6) * (t337 * t153 + t341 * t165), t146, t151, t152, t160, t179, -t306 * t239 + t353; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t337 * t146 + t341 * t149 + pkin(6) * (t341 * t153 - t337 * t165), t149, t148, t156, t158, t173, -t267 * t237 + t355; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t352, t352, qJDD(1) * t360 - t382, -t350, t383, t384, -t266 * t241 - t354;];
m_new  = t1;
