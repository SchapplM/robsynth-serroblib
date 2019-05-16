% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPP5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:23:54
% EndTime: 2019-05-07 18:24:13
% DurationCPUTime: 8.20s
% Computational Cost: add. (145524->382), mult. (289796->444), div. (0->0), fcn. (200309->8), ass. (0->137)
t337 = sin(qJ(1));
t340 = cos(qJ(1));
t324 = t337 * g(1) - t340 * g(2);
t342 = qJD(1) ^ 2;
t307 = -qJDD(1) * pkin(1) - t342 * pkin(7) - t324;
t336 = sin(qJ(2));
t339 = cos(qJ(2));
t364 = qJD(1) * qJD(2);
t362 = t339 * t364;
t318 = t336 * qJDD(1) + t362;
t330 = t336 * t364;
t319 = t339 * qJDD(1) - t330;
t260 = (-t318 - t362) * pkin(8) + (-t319 + t330) * pkin(2) + t307;
t325 = -t340 * g(1) - t337 * g(2);
t308 = -t342 * pkin(1) + qJDD(1) * pkin(7) + t325;
t289 = -t336 * g(3) + t339 * t308;
t317 = (-pkin(2) * t339 - pkin(8) * t336) * qJD(1);
t341 = qJD(2) ^ 2;
t365 = t339 * qJD(1);
t265 = -t341 * pkin(2) + qJDD(2) * pkin(8) + t317 * t365 + t289;
t335 = sin(qJ(3));
t338 = cos(qJ(3));
t232 = t338 * t260 - t335 * t265;
t366 = qJD(1) * t336;
t314 = t338 * qJD(2) - t335 * t366;
t281 = t314 * qJD(3) + t335 * qJDD(2) + t338 * t318;
t313 = qJDD(3) - t319;
t315 = t335 * qJD(2) + t338 * t366;
t329 = qJD(3) - t365;
t201 = (t314 * t329 - t281) * pkin(9) + (t314 * t315 + t313) * pkin(3) + t232;
t233 = t335 * t260 + t338 * t265;
t280 = -t315 * qJD(3) + t338 * qJDD(2) - t335 * t318;
t290 = t329 * pkin(3) - t315 * pkin(9);
t312 = t314 ^ 2;
t204 = -t312 * pkin(3) + t280 * pkin(9) - t329 * t290 + t233;
t334 = sin(qJ(4));
t373 = cos(qJ(4));
t199 = t334 * t201 + t373 * t204;
t283 = -t373 * t314 + t334 * t315;
t284 = t334 * t314 + t373 * t315;
t254 = t283 * pkin(4) - t284 * qJ(5);
t309 = -qJDD(4) - t313;
t327 = -qJD(4) - t329;
t326 = t327 ^ 2;
t374 = -2 * qJD(5);
t192 = -t326 * pkin(4) - t309 * qJ(5) - t283 * t254 + t327 * t374 + t199;
t272 = t327 * mrSges(6,1) + t284 * mrSges(6,2);
t378 = m(6) * t192 - t309 * mrSges(6,3) - t327 * t272;
t198 = t373 * t201 - t334 * t204;
t194 = t309 * pkin(4) - t326 * qJ(5) + t284 * t254 + qJDD(5) - t198;
t266 = -t283 * mrSges(6,2) - t327 * mrSges(6,3);
t377 = -m(6) * t194 - t309 * mrSges(6,1) - t327 * t266;
t231 = -t283 * qJD(4) + t334 * t280 + t373 * t281;
t288 = -t339 * g(3) - t336 * t308;
t353 = qJDD(2) * pkin(2) + t341 * pkin(8) - t317 * t366 + t288;
t351 = t280 * pkin(3) + t312 * pkin(9) - t315 * t290 + t353;
t370 = t283 * t327;
t376 = -(t231 + t370) * qJ(5) - t351;
t230 = t284 * qJD(4) - t373 * t280 + t334 * t281;
t269 = t327 * pkin(5) - t284 * qJ(6);
t282 = t283 ^ 2;
t189 = -t282 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t230 + (pkin(4) * t327 + (2 * qJD(5)) + t269) * t284 - t376;
t267 = -t327 * mrSges(7,2) + t283 * mrSges(7,3);
t270 = t327 * mrSges(7,1) - t284 * mrSges(7,3);
t179 = -m(7) * t189 + t230 * mrSges(7,1) - t231 * mrSges(7,2) + t283 * t267 - t284 * t270;
t196 = t284 * t374 + (-t284 * t327 + t230) * pkin(4) + t376;
t173 = m(6) * t196 + t230 * mrSges(6,1) - t231 * mrSges(6,3) + t283 * t266 - t284 * t272 + t179;
t240 = Ifges(7,5) * t284 + Ifges(7,6) * t283 + Ifges(7,3) * t327;
t247 = Ifges(6,1) * t284 - Ifges(6,4) * t327 + Ifges(6,5) * t283;
t248 = Ifges(5,1) * t284 - Ifges(5,4) * t283 - Ifges(5,5) * t327;
t186 = -t282 * pkin(5) + t230 * qJ(6) + 0.2e1 * qJD(6) * t283 - t327 * t269 + t192;
t256 = -t283 * mrSges(7,1) + t284 * mrSges(7,2);
t363 = m(7) * t186 + t230 * mrSges(7,3) + t283 * t256;
t180 = -t309 * mrSges(7,2) - t327 * t270 + t363;
t246 = Ifges(7,1) * t284 + Ifges(7,4) * t283 + Ifges(7,5) * t327;
t355 = mrSges(7,1) * t189 - mrSges(7,3) * t186 - Ifges(7,4) * t231 - Ifges(7,2) * t230 - Ifges(7,6) * t309 - t327 * t246;
t349 = mrSges(6,1) * t196 - mrSges(6,2) * t192 + pkin(5) * t179 + qJ(6) * t180 - t355;
t244 = Ifges(6,4) * t284 - Ifges(6,2) * t327 + Ifges(6,6) * t283;
t368 = -Ifges(5,5) * t284 + Ifges(5,6) * t283 + Ifges(5,3) * t327 - t244;
t160 = -t349 + (-t248 - t247) * t327 + (-Ifges(5,6) + Ifges(6,6)) * t309 + (t240 + t368) * t284 + (Ifges(5,4) - Ifges(6,5)) * t231 + (-Ifges(5,2) - Ifges(6,3)) * t230 + mrSges(5,1) * t351 + mrSges(5,3) * t199 - pkin(4) * t173;
t243 = Ifges(7,4) * t284 + Ifges(7,2) * t283 + Ifges(7,6) * t327;
t245 = Ifges(5,4) * t284 - Ifges(5,2) * t283 - Ifges(5,6) * t327;
t182 = -0.2e1 * qJD(6) * t284 + (-t231 + t370) * qJ(6) + (t283 * t284 + t309) * pkin(5) + t194;
t358 = -m(7) * t182 + t231 * mrSges(7,3) + t284 * t256;
t178 = t309 * mrSges(7,1) + t327 * t267 - t358;
t241 = Ifges(6,5) * t284 - Ifges(6,6) * t327 + Ifges(6,3) * t283;
t354 = mrSges(7,2) * t189 - mrSges(7,3) * t182 + Ifges(7,1) * t231 + Ifges(7,4) * t230 + Ifges(7,5) * t309 + t283 * t240;
t348 = mrSges(6,2) * t194 - mrSges(6,3) * t196 + Ifges(6,1) * t231 - Ifges(6,4) * t309 + Ifges(6,5) * t230 - qJ(6) * t178 - t327 * t241 + t354;
t161 = t368 * t283 + (t245 - t243) * t327 - Ifges(5,5) * t309 - Ifges(5,4) * t230 + Ifges(5,1) * t231 - mrSges(5,2) * t351 - mrSges(5,3) * t198 - qJ(5) * t173 + t348;
t274 = Ifges(4,5) * t315 + Ifges(4,6) * t314 + Ifges(4,3) * t329;
t276 = Ifges(4,1) * t315 + Ifges(4,4) * t314 + Ifges(4,5) * t329;
t268 = t327 * mrSges(5,2) - t283 * mrSges(5,3);
t271 = -t327 * mrSges(5,1) - t284 * mrSges(5,3);
t347 = -m(5) * t351 + t230 * mrSges(5,1) + t231 * mrSges(5,2) + t283 * t268 + t284 * t271 + t173;
t255 = t283 * mrSges(6,1) - t284 * mrSges(6,3);
t367 = -t283 * mrSges(5,1) - t284 * mrSges(5,2) - t255;
t371 = -mrSges(5,3) - mrSges(6,2);
t169 = m(5) * t198 + (-t267 - t268) * t327 + (-mrSges(5,1) - mrSges(7,1)) * t309 + t367 * t284 + t371 * t231 + t358 + t377;
t172 = m(5) * t199 + (-t270 + t271) * t327 + (mrSges(5,2) - mrSges(7,2)) * t309 + t367 * t283 + t371 * t230 + t363 + t378;
t360 = -t334 * t169 + t373 * t172;
t149 = mrSges(4,1) * t353 + mrSges(4,3) * t233 + Ifges(4,4) * t281 + Ifges(4,2) * t280 + Ifges(4,6) * t313 - pkin(3) * t347 + pkin(9) * t360 + t373 * t160 + t334 * t161 - t315 * t274 + t329 * t276;
t165 = t373 * t169 + t334 * t172;
t275 = Ifges(4,4) * t315 + Ifges(4,2) * t314 + Ifges(4,6) * t329;
t150 = -mrSges(4,2) * t353 - mrSges(4,3) * t232 + Ifges(4,1) * t281 + Ifges(4,4) * t280 + Ifges(4,5) * t313 - pkin(9) * t165 - t334 * t160 + t373 * t161 + t314 * t274 - t329 * t275;
t285 = -t314 * mrSges(4,1) + t315 * mrSges(4,2);
t286 = -t329 * mrSges(4,2) + t314 * mrSges(4,3);
t163 = m(4) * t232 + t313 * mrSges(4,1) - t281 * mrSges(4,3) - t315 * t285 + t329 * t286 + t165;
t287 = t329 * mrSges(4,1) - t315 * mrSges(4,3);
t164 = m(4) * t233 - t313 * mrSges(4,2) + t280 * mrSges(4,3) + t314 * t285 - t329 * t287 + t360;
t159 = -t335 * t163 + t338 * t164;
t167 = m(4) * t353 + t280 * mrSges(4,1) - t281 * mrSges(4,2) + t314 * t286 - t315 * t287 - t347;
t305 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t336 + Ifges(3,2) * t339) * qJD(1);
t306 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t336 + Ifges(3,4) * t339) * qJD(1);
t375 = mrSges(3,1) * t288 - mrSges(3,2) * t289 + Ifges(3,5) * t318 + Ifges(3,6) * t319 + Ifges(3,3) * qJDD(2) + pkin(2) * t167 + pkin(8) * t159 + t338 * t149 + t335 * t150 + (t336 * t305 - t339 * t306) * qJD(1);
t369 = t327 * t243;
t316 = (-mrSges(3,1) * t339 + mrSges(3,2) * t336) * qJD(1);
t322 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t366;
t157 = m(3) * t289 - qJDD(2) * mrSges(3,2) + t319 * mrSges(3,3) - qJD(2) * t322 + t316 * t365 + t159;
t323 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t365;
t166 = m(3) * t288 + qJDD(2) * mrSges(3,1) - t318 * mrSges(3,3) + qJD(2) * t323 - t316 * t366 + t167;
t361 = t339 * t157 - t336 * t166;
t158 = t338 * t163 + t335 * t164;
t356 = mrSges(7,1) * t182 - mrSges(7,2) * t186 + Ifges(7,5) * t231 + Ifges(7,6) * t230 + Ifges(7,3) * t309 + t284 * t243 - t283 * t246;
t304 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t336 + Ifges(3,6) * t339) * qJD(1);
t146 = mrSges(3,2) * t307 - mrSges(3,3) * t288 + Ifges(3,1) * t318 + Ifges(3,4) * t319 + Ifges(3,5) * qJDD(2) - pkin(8) * t158 - qJD(2) * t305 - t335 * t149 + t338 * t150 + t304 * t365;
t346 = mrSges(6,1) * t194 - mrSges(6,3) * t192 - Ifges(6,4) * t231 + Ifges(6,2) * t309 - Ifges(6,6) * t230 + pkin(5) * t178 + t284 * t241 - t283 * t247 + t356;
t344 = mrSges(5,2) * t199 - t283 * t248 - pkin(4) * (-t231 * mrSges(6,2) - t284 * t255 - t178 + t377) - qJ(5) * (-t230 * mrSges(6,2) - t283 * t255 + t180 + t378) - mrSges(5,1) * t198 + Ifges(5,6) * t230 - Ifges(5,5) * t231 - t284 * t245 + Ifges(5,3) * t309 + t346;
t343 = mrSges(4,1) * t232 - mrSges(4,2) * t233 + Ifges(4,5) * t281 + Ifges(4,6) * t280 + Ifges(4,3) * t313 + pkin(3) * t165 + t315 * t275 - t314 * t276 - t344;
t148 = -mrSges(3,1) * t307 + mrSges(3,3) * t289 + Ifges(3,4) * t318 + Ifges(3,2) * t319 + Ifges(3,6) * qJDD(2) - pkin(2) * t158 + qJD(2) * t306 - t304 * t366 - t343;
t350 = -m(3) * t307 + t319 * mrSges(3,1) - t318 * mrSges(3,2) - t322 * t366 + t323 * t365 - t158;
t352 = mrSges(2,1) * t324 - mrSges(2,2) * t325 + Ifges(2,3) * qJDD(1) + pkin(1) * t350 + pkin(7) * t361 + t336 * t146 + t339 * t148;
t154 = m(2) * t324 + qJDD(1) * mrSges(2,1) - t342 * mrSges(2,2) + t350;
t153 = t336 * t157 + t339 * t166;
t151 = m(2) * t325 - t342 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t361;
t144 = mrSges(2,1) * g(3) + mrSges(2,3) * t325 + t342 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t153 - t375;
t143 = -mrSges(2,2) * g(3) - mrSges(2,3) * t324 + Ifges(2,5) * qJDD(1) - t342 * Ifges(2,6) - pkin(7) * t153 + t339 * t146 - t336 * t148;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t340 * t143 - t337 * t144 - pkin(6) * (t337 * t151 + t340 * t154), t143, t146, t150, t161, -t283 * t244 + t348 - t369, t354 - t369; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t337 * t143 + t340 * t144 + pkin(6) * (t340 * t151 - t337 * t154), t144, t148, t149, t160, -t346, -t284 * t240 - t355; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t352, t352, t375, t343, -t344, t349 + t327 * t247 - Ifges(6,6) * t309 + Ifges(6,3) * t230 + Ifges(6,5) * t231 + (-t240 + t244) * t284, t356;];
m_new  = t1;
