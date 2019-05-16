% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 05:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:26:53
% EndTime: 2019-05-07 05:27:52
% DurationCPUTime: 31.03s
% Computational Cost: add. (540263->402), mult. (1197481->500), div. (0->0), fcn. (949103->12), ass. (0->161)
t353 = sin(qJ(2));
t357 = cos(qJ(2));
t349 = sin(pkin(6));
t383 = qJD(1) * t349;
t333 = (-pkin(2) * t357 - pkin(9) * t353) * t383;
t350 = cos(pkin(6));
t344 = qJD(1) * t350 + qJD(2);
t342 = t344 ^ 2;
t343 = qJDD(1) * t350 + qJDD(2);
t382 = qJD(1) * t357;
t354 = sin(qJ(1));
t358 = cos(qJ(1));
t340 = t354 * g(1) - g(2) * t358;
t359 = qJD(1) ^ 2;
t396 = pkin(8) * t349;
t330 = qJDD(1) * pkin(1) + t359 * t396 + t340;
t341 = -g(1) * t358 - g(2) * t354;
t331 = -pkin(1) * t359 + qJDD(1) * t396 + t341;
t389 = t350 * t353;
t384 = t330 * t389 + t357 * t331;
t279 = -t342 * pkin(2) + t343 * pkin(9) + (-g(3) * t353 + t333 * t382) * t349 + t384;
t380 = qJD(1) * qJD(2);
t334 = (qJDD(1) * t353 + t357 * t380) * t349;
t335 = (-qJDD(1) * t357 + t353 * t380) * t349;
t395 = t350 * g(3);
t280 = t335 * pkin(2) - t334 * pkin(9) - t395 + (-t330 + (pkin(2) * t353 - pkin(9) * t357) * t344 * qJD(1)) * t349;
t352 = sin(qJ(3));
t356 = cos(qJ(3));
t238 = -t352 * t279 + t356 * t280;
t379 = t353 * t383;
t320 = t344 * t356 - t352 * t379;
t298 = qJD(3) * t320 + t334 * t356 + t343 * t352;
t321 = t344 * t352 + t356 * t379;
t327 = qJDD(3) + t335;
t378 = t349 * t382;
t338 = -qJD(3) + t378;
t227 = (-t320 * t338 - t298) * qJ(4) + (t320 * t321 + t327) * pkin(3) + t238;
t239 = t356 * t279 + t352 * t280;
t297 = -qJD(3) * t321 - t334 * t352 + t343 * t356;
t310 = -pkin(3) * t338 - qJ(4) * t321;
t319 = t320 ^ 2;
t230 = -pkin(3) * t319 + qJ(4) * t297 + t310 * t338 + t239;
t348 = sin(pkin(11));
t393 = cos(pkin(11));
t307 = t348 * t320 + t321 * t393;
t399 = -2 * qJD(4);
t222 = t227 * t393 - t348 * t230 + t307 * t399;
t306 = -t320 * t393 + t321 * t348;
t302 = t306 * t399;
t387 = t348 * t227 + t393 * t230;
t223 = t302 + t387;
t261 = -t297 * t393 + t298 * t348;
t262 = t348 * t297 + t298 * t393;
t266 = Ifges(5,4) * t307 - Ifges(5,2) * t306 - Ifges(5,6) * t338;
t274 = -mrSges(6,2) * t306 - mrSges(6,3) * t307;
t284 = mrSges(6,1) * t306 + mrSges(6,3) * t338;
t272 = pkin(4) * t306 - qJ(5) * t307;
t337 = t338 ^ 2;
t220 = -t327 * pkin(4) - t337 * qJ(5) + t307 * t272 + qJDD(5) - t222;
t392 = t306 * t338;
t214 = (t306 * t307 - t327) * pkin(10) + (t262 - t392) * pkin(5) + t220;
t288 = pkin(5) * t307 + pkin(10) * t338;
t305 = t306 ^ 2;
t388 = t350 * t357;
t390 = t349 * t357;
t299 = -g(3) * t390 + t330 * t388 - t353 * t331;
t278 = -t343 * pkin(2) - t342 * pkin(9) + t333 * t379 - t299;
t232 = -t297 * pkin(3) - t319 * qJ(4) + t321 * t310 + qJDD(4) + t278;
t397 = -2 * qJD(5);
t361 = (-t262 - t392) * qJ(5) + t232 + (-t338 * pkin(4) + t397) * t307;
t217 = -t307 * t288 - t305 * pkin(5) + (pkin(4) + pkin(10)) * t261 + t361;
t351 = sin(qJ(6));
t355 = cos(qJ(6));
t211 = t214 * t355 - t217 * t351;
t282 = t306 * t355 + t338 * t351;
t237 = qJD(6) * t282 + t261 * t351 + t327 * t355;
t283 = t306 * t351 - t338 * t355;
t248 = -mrSges(7,1) * t282 + mrSges(7,2) * t283;
t304 = qJD(6) + t307;
t249 = -mrSges(7,2) * t304 + mrSges(7,3) * t282;
t260 = qJDD(6) + t262;
t208 = m(7) * t211 + mrSges(7,1) * t260 - mrSges(7,3) * t237 - t248 * t283 + t249 * t304;
t212 = t214 * t351 + t217 * t355;
t236 = -qJD(6) * t283 + t261 * t355 - t327 * t351;
t250 = mrSges(7,1) * t304 - mrSges(7,3) * t283;
t209 = m(7) * t212 - t260 * mrSges(7,2) + mrSges(7,3) * t236 + t248 * t282 - t250 * t304;
t196 = t355 * t208 + t351 * t209;
t372 = t337 * pkin(4) - t327 * qJ(5) - t387;
t216 = -t261 * pkin(5) - t305 * pkin(10) - t306 * t272 + t302 + (t397 - t288) * t338 - t372;
t240 = Ifges(7,5) * t283 + Ifges(7,6) * t282 + Ifges(7,3) * t304;
t242 = Ifges(7,1) * t283 + Ifges(7,4) * t282 + Ifges(7,5) * t304;
t201 = -mrSges(7,1) * t216 + mrSges(7,3) * t212 + Ifges(7,4) * t237 + Ifges(7,2) * t236 + Ifges(7,6) * t260 - t240 * t283 + t242 * t304;
t241 = Ifges(7,4) * t283 + Ifges(7,2) * t282 + Ifges(7,6) * t304;
t202 = mrSges(7,2) * t216 - mrSges(7,3) * t211 + Ifges(7,1) * t237 + Ifges(7,4) * t236 + Ifges(7,5) * t260 + t240 * t282 - t241 * t304;
t218 = 0.2e1 * qJD(5) * t338 + ((2 * qJD(4)) + t272) * t306 + t372;
t263 = -Ifges(6,5) * t338 - Ifges(6,6) * t307 + Ifges(6,3) * t306;
t367 = -mrSges(6,2) * t220 + mrSges(6,3) * t218 - Ifges(6,1) * t327 + Ifges(6,4) * t262 - Ifges(6,5) * t261 + pkin(10) * t196 + t351 * t201 - t355 * t202 + t307 * t263;
t213 = -m(7) * t216 + t236 * mrSges(7,1) - t237 * mrSges(7,2) + t282 * t249 - t283 * t250;
t285 = mrSges(6,1) * t307 - mrSges(6,2) * t338;
t368 = -m(6) * t218 + t327 * mrSges(6,3) - t338 * t285 - t213;
t371 = -m(6) * t220 - t262 * mrSges(6,1) - t307 * t274 - t196;
t265 = -Ifges(6,4) * t338 - Ifges(6,2) * t307 + Ifges(6,6) * t306;
t385 = -Ifges(5,1) * t307 + Ifges(5,4) * t306 + Ifges(5,5) * t338 + t265;
t400 = -t385 * t306 - mrSges(5,2) * t223 + pkin(4) * (-t327 * mrSges(6,2) + t338 * t284 + t371) + qJ(5) * (-t261 * mrSges(6,1) - t306 * t274 + t368) + mrSges(5,1) * t222 + t307 * t266 - Ifges(5,6) * t261 + Ifges(5,5) * t262 + Ifges(5,3) * t327 - t367;
t273 = mrSges(5,1) * t306 + mrSges(5,2) * t307;
t286 = mrSges(5,2) * t338 - mrSges(5,3) * t306;
t190 = m(5) * t222 - t262 * mrSges(5,3) - t307 * t273 + (t284 - t286) * t338 + (mrSges(5,1) - mrSges(6,2)) * t327 + t371;
t287 = -mrSges(5,1) * t338 - mrSges(5,3) * t307;
t203 = m(5) * t223 - t327 * mrSges(5,2) + t338 * t287 + (-t273 - t274) * t306 + (-mrSges(5,3) - mrSges(6,1)) * t261 + t368;
t188 = t393 * t190 + t348 * t203;
t292 = Ifges(4,4) * t321 + Ifges(4,2) * t320 - Ifges(4,6) * t338;
t293 = Ifges(4,1) * t321 + Ifges(4,4) * t320 - Ifges(4,5) * t338;
t398 = mrSges(4,1) * t238 - mrSges(4,2) * t239 + Ifges(4,5) * t298 + Ifges(4,6) * t297 + Ifges(4,3) * t327 + pkin(3) * t188 + t321 * t292 - t320 * t293 + t400;
t394 = Ifges(5,4) + Ifges(6,6);
t391 = t349 * t353;
t308 = -mrSges(4,1) * t320 + mrSges(4,2) * t321;
t309 = mrSges(4,2) * t338 + mrSges(4,3) * t320;
t186 = m(4) * t238 + mrSges(4,1) * t327 - mrSges(4,3) * t298 - t308 * t321 - t309 * t338 + t188;
t311 = -mrSges(4,1) * t338 - mrSges(4,3) * t321;
t375 = -t190 * t348 + t393 * t203;
t187 = m(4) * t239 - mrSges(4,2) * t327 + mrSges(4,3) * t297 + t308 * t320 + t311 * t338 + t375;
t181 = t356 * t186 + t352 * t187;
t197 = -t351 * t208 + t355 * t209;
t267 = -Ifges(6,1) * t338 - Ifges(6,4) * t307 + Ifges(6,5) * t306;
t386 = -Ifges(5,5) * t307 + Ifges(5,6) * t306 + Ifges(5,3) * t338 - t267;
t300 = -g(3) * t391 + t384;
t328 = mrSges(3,1) * t344 - mrSges(3,3) * t379;
t332 = (-mrSges(3,1) * t357 + mrSges(3,2) * t353) * t383;
t376 = -t186 * t352 + t356 * t187;
t179 = m(3) * t300 - mrSges(3,2) * t343 - mrSges(3,3) * t335 - t328 * t344 + t332 * t378 + t376;
t329 = -mrSges(3,2) * t344 + mrSges(3,3) * t378;
t225 = t261 * pkin(4) + t361;
t195 = m(6) * t225 - t261 * mrSges(6,2) - t262 * mrSges(6,3) - t306 * t284 - t307 * t285 + t197;
t365 = m(5) * t232 + t261 * mrSges(5,1) + t262 * mrSges(5,2) + t306 * t286 + t307 * t287 + t195;
t362 = -m(4) * t278 + t297 * mrSges(4,1) - t298 * mrSges(4,2) + t320 * t309 - t321 * t311 - t365;
t192 = m(3) * t299 + t343 * mrSges(3,1) - t334 * mrSges(3,3) + t344 * t329 - t332 * t379 + t362;
t175 = t357 * t179 - t192 * t353;
t315 = -t349 * t330 - t395;
t180 = m(3) * t315 + t335 * mrSges(3,1) + t334 * mrSges(3,2) + (t328 * t353 - t329 * t357) * t383 + t181;
t172 = t179 * t389 - t180 * t349 + t192 * t388;
t366 = -mrSges(6,1) * t218 + mrSges(6,2) * t225 - pkin(5) * t213 - pkin(10) * t197 - t355 * t201 - t351 * t202;
t176 = -mrSges(5,1) * t232 + mrSges(5,3) * t223 - pkin(4) * t195 + t385 * t338 + (Ifges(5,6) - Ifges(6,5)) * t327 + t386 * t307 + t394 * t262 + (-Ifges(5,2) - Ifges(6,3)) * t261 + t366;
t369 = mrSges(7,1) * t211 - mrSges(7,2) * t212 + Ifges(7,5) * t237 + Ifges(7,6) * t236 + Ifges(7,3) * t260 + t283 * t241 - t282 * t242;
t364 = mrSges(6,1) * t220 - mrSges(6,3) * t225 + pkin(5) * t196 + t369;
t182 = (t266 - t263) * t338 + (Ifges(5,5) - Ifges(6,4)) * t327 + t386 * t306 + (Ifges(5,1) + Ifges(6,2)) * t262 - t394 * t261 + t364 + mrSges(5,2) * t232 - mrSges(5,3) * t222 - qJ(5) * t195;
t291 = Ifges(4,5) * t321 + Ifges(4,6) * t320 - Ifges(4,3) * t338;
t165 = -mrSges(4,1) * t278 + mrSges(4,3) * t239 + Ifges(4,4) * t298 + Ifges(4,2) * t297 + Ifges(4,6) * t327 - pkin(3) * t365 + qJ(4) * t375 + t176 * t393 + t348 * t182 - t321 * t291 - t338 * t293;
t168 = mrSges(4,2) * t278 - mrSges(4,3) * t238 + Ifges(4,1) * t298 + Ifges(4,4) * t297 + Ifges(4,5) * t327 - qJ(4) * t188 - t348 * t176 + t182 * t393 + t320 * t291 + t338 * t292;
t313 = Ifges(3,6) * t344 + (Ifges(3,4) * t353 + Ifges(3,2) * t357) * t383;
t314 = Ifges(3,5) * t344 + (Ifges(3,1) * t353 + Ifges(3,4) * t357) * t383;
t162 = Ifges(3,5) * t334 - Ifges(3,6) * t335 + Ifges(3,3) * t343 + mrSges(3,1) * t299 - mrSges(3,2) * t300 + t352 * t168 + t356 * t165 + pkin(2) * t362 + pkin(9) * t376 + (t313 * t353 - t314 * t357) * t383;
t312 = Ifges(3,3) * t344 + (Ifges(3,5) * t353 + Ifges(3,6) * t357) * t383;
t164 = mrSges(3,2) * t315 - mrSges(3,3) * t299 + Ifges(3,1) * t334 - Ifges(3,4) * t335 + Ifges(3,5) * t343 - pkin(9) * t181 - t165 * t352 + t168 * t356 + t312 * t378 - t313 * t344;
t167 = -mrSges(3,1) * t315 + mrSges(3,3) * t300 + Ifges(3,4) * t334 - Ifges(3,2) * t335 + Ifges(3,6) * t343 - pkin(2) * t181 - t312 * t379 + t344 * t314 - t398;
t370 = mrSges(2,1) * t340 - mrSges(2,2) * t341 + Ifges(2,3) * qJDD(1) + pkin(1) * t172 + t350 * t162 + t164 * t391 + t167 * t390 + t175 * t396;
t173 = m(2) * t341 - mrSges(2,1) * t359 - qJDD(1) * mrSges(2,2) + t175;
t171 = t350 * t180 + (t179 * t353 + t192 * t357) * t349;
t169 = m(2) * t340 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t359 + t172;
t160 = -mrSges(2,2) * g(3) - mrSges(2,3) * t340 + Ifges(2,5) * qJDD(1) - t359 * Ifges(2,6) + t357 * t164 - t353 * t167 + (-t171 * t349 - t172 * t350) * pkin(8);
t159 = mrSges(2,1) * g(3) + mrSges(2,3) * t341 + t359 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t171 - t349 * t162 + (pkin(8) * t175 + t164 * t353 + t167 * t357) * t350;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t358 * t160 - t354 * t159 - pkin(7) * (t169 * t358 + t173 * t354), t160, t164, t168, t182, -t306 * t265 - t367, t202; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t354 * t160 + t358 * t159 + pkin(7) * (-t169 * t354 + t173 * t358), t159, t167, t165, t176, Ifges(6,4) * t327 - Ifges(6,2) * t262 + Ifges(6,6) * t261 + t338 * t263 + t306 * t267 - t364, t201; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t370, t370, t162, t398, t400, Ifges(6,5) * t327 - Ifges(6,6) * t262 + Ifges(6,3) * t261 - t338 * t265 + t307 * t267 - t366, t369;];
m_new  = t1;
