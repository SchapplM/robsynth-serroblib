% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-05-05 17:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:22:16
% EndTime: 2019-05-05 17:22:22
% DurationCPUTime: 3.17s
% Computational Cost: add. (27520->359), mult. (55017->399), div. (0->0), fcn. (24462->6), ass. (0->144)
t362 = (qJD(2) * qJD(1));
t301 = 2 * t362;
t331 = qJD(1) ^ 2;
t326 = sin(qJ(1));
t329 = cos(qJ(1));
t294 = -t329 * g(1) - t326 * g(2);
t360 = t331 * pkin(1) - qJDD(1) * qJ(2) - t294;
t358 = -t331 * pkin(7) - t360;
t225 = t301 + t358;
t325 = sin(qJ(3));
t328 = cos(qJ(3));
t363 = qJD(1) * qJD(3);
t280 = t325 * qJDD(1) + t328 * t363;
t299 = t325 * t363;
t281 = t328 * qJDD(1) - t299;
t366 = qJD(1) * t325;
t287 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t366;
t365 = t328 * qJD(1);
t290 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t365;
t355 = t280 * pkin(3) - t281 * qJ(4) + t358;
t359 = pkin(3) * t328 + qJ(4) * t325;
t391 = -2 * qJD(4);
t201 = t301 + (qJD(3) * t359 + t328 * t391) * qJD(1) + t355;
t289 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t365;
t291 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t365;
t292 = -mrSges(5,2) * t366 + qJD(3) * mrSges(5,3);
t288 = -qJD(3) * pkin(4) - qJ(5) * t365;
t302 = -2 * t362;
t373 = t325 ^ 2 * t331;
t390 = 2 * qJD(4);
t342 = -qJ(5) * t373 + qJDD(5) + t302 - t355 + (t288 + t390) * t365;
t377 = -pkin(5) - qJ(4);
t387 = -pkin(4) - pkin(8);
t388 = -pkin(3) - pkin(8);
t189 = (t377 * t325 + t388 * t328) * t363 + t342 + t281 * pkin(5) + t387 * t280;
t293 = t326 * g(1) - t329 * g(2);
t354 = -t331 * qJ(2) + qJDD(2) - t293;
t226 = (-pkin(1) - pkin(7)) * qJDD(1) + t354;
t279 = (pkin(5) * t328 - pkin(8) * t325) * qJD(1);
t330 = qJD(3) ^ 2;
t385 = t325 * g(3);
t389 = -2 * qJD(5);
t276 = (pkin(3) * t325 - qJ(4) * t328) * qJD(1);
t395 = t276 * t365 + qJDD(4);
t343 = t328 * t331 * t325 * pkin(4) + t365 * t389 + (-t281 - t299) * qJ(5) - t385 + t395;
t192 = t377 * t330 + (-qJD(1) * t279 - t226) * t328 + (-pkin(3) + t387) * qJDD(3) + t343;
t324 = sin(qJ(6));
t327 = cos(qJ(6));
t186 = t327 * t189 - t324 * t192;
t273 = -t327 * qJD(3) - t324 * t366;
t216 = t273 * qJD(6) - t324 * qJDD(3) + t327 * t280;
t274 = -t324 * qJD(3) + t327 * t366;
t218 = -t273 * mrSges(7,1) + t274 * mrSges(7,2);
t296 = qJD(6) + t365;
t221 = -t296 * mrSges(7,2) + t273 * mrSges(7,3);
t271 = qJDD(6) + t281;
t182 = m(7) * t186 + t271 * mrSges(7,1) - t216 * mrSges(7,3) - t274 * t218 + t296 * t221;
t187 = t324 * t189 + t327 * t192;
t215 = -t274 * qJD(6) - t327 * qJDD(3) - t324 * t280;
t222 = t296 * mrSges(7,1) - t274 * mrSges(7,3);
t183 = m(7) * t187 - t271 * mrSges(7,2) + t215 * mrSges(7,3) + t273 * t218 - t296 * t222;
t166 = t327 * t182 + t324 * t183;
t195 = -t280 * pkin(4) - t359 * t363 + t342;
t356 = -m(6) * t195 - t281 * mrSges(6,1) - t166;
t286 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t366;
t372 = t325 * t286;
t338 = m(5) * t201 + t280 * mrSges(5,1) - t281 * mrSges(5,3) - ((t289 + t291) * t328 + t372) * qJD(1) + t292 * t366 + t356;
t382 = mrSges(4,1) - mrSges(6,2);
t398 = -m(4) * t225 - t281 * mrSges(4,2) - t382 * t280 - t287 * t366 - t290 * t365 - t338;
t220 = -t328 * g(3) + t325 * t226;
t397 = -qJDD(3) * qJ(4) - t220;
t396 = Ifges(6,4) * t280 - Ifges(6,2) * t281;
t371 = t328 * t226;
t219 = t371 + t385;
t370 = t330 * qJ(4);
t206 = -qJDD(3) * pkin(3) - t219 - t370 + t395;
t394 = -m(5) * t206 + qJDD(3) * mrSges(5,1) + qJD(3) * t292;
t384 = t330 * pkin(3);
t383 = mrSges(2,1) - mrSges(3,2);
t381 = -mrSges(4,3) - mrSges(5,2);
t380 = Ifges(2,5) - Ifges(3,4);
t379 = -Ifges(2,6) + Ifges(3,5);
t378 = Ifges(5,6) - Ifges(6,5);
t374 = t280 * mrSges(6,2);
t369 = -t324 * t182 + t327 * t183;
t246 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t328 + Ifges(5,6) * t325) * qJD(1);
t368 = -Ifges(4,3) * qJD(3) - (Ifges(4,5) * t328 - Ifges(4,6) * t325) * qJD(1) - t246;
t364 = t389 + t276;
t199 = -t370 - t371 + (-pkin(3) - pkin(4)) * qJDD(3) + t343;
t275 = (mrSges(6,1) * t328 + mrSges(6,2) * t325) * qJD(1);
t357 = -m(6) * t199 + t275 * t365 - t369;
t277 = (mrSges(5,1) * t325 - mrSges(5,3) * t328) * qJD(1);
t361 = qJD(1) * (-t277 - (mrSges(4,1) * t325 + mrSges(4,2) * t328) * qJD(1));
t158 = m(4) * t219 + t382 * qJDD(3) + (-t286 + t287) * qJD(3) + t328 * t361 + (mrSges(6,3) + t381) * t281 + t357 + t394;
t350 = pkin(4) * t373 - t280 * qJ(5) + t397;
t198 = t384 + (t391 - t288) * qJD(3) + t364 * t366 + t350;
t300 = qJD(3) * t390;
t191 = qJDD(3) * pkin(5) + qJD(3) * t288 + t300 + t388 * t330 + (t279 - t364) * t366 - t350;
t353 = -m(7) * t191 + t215 * mrSges(7,1) - t216 * mrSges(7,2) + t273 * t221 - t274 * t222;
t177 = -m(6) * t198 + qJDD(3) * mrSges(6,1) + t280 * mrSges(6,3) + qJD(3) * t289 + t275 * t366 - t353;
t204 = -t276 * t366 + t300 - t384 - t397;
t340 = m(5) * t204 + qJDD(3) * mrSges(5,3) + qJD(3) * t291 + t177;
t169 = m(4) * t220 - qJDD(3) * mrSges(4,2) - qJD(3) * t290 + t381 * t280 + t325 * t361 + t340;
t154 = -t325 * t158 + t328 * t169;
t153 = t328 * t158 + t325 * t169;
t238 = -qJDD(1) * pkin(1) + t354;
t352 = -m(3) * t238 + t331 * mrSges(3,3) - t153;
t209 = Ifges(7,4) * t274 + Ifges(7,2) * t273 + Ifges(7,6) * t296;
t210 = Ifges(7,1) * t274 + Ifges(7,4) * t273 + Ifges(7,5) * t296;
t351 = mrSges(7,1) * t186 - mrSges(7,2) * t187 + Ifges(7,5) * t216 + Ifges(7,6) * t215 + Ifges(7,3) * t271 + t274 * t209 - t273 * t210;
t159 = t338 - t374;
t242 = -Ifges(6,3) * qJD(3) + (Ifges(6,5) * t325 - Ifges(6,6) * t328) * qJD(1);
t249 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t328 + Ifges(5,5) * t325) * qJD(1);
t250 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t328 - Ifges(4,4) * t325) * qJD(1);
t208 = Ifges(7,5) * t274 + Ifges(7,6) * t273 + Ifges(7,3) * t296;
t173 = -mrSges(7,1) * t191 + mrSges(7,3) * t187 + Ifges(7,4) * t216 + Ifges(7,2) * t215 + Ifges(7,6) * t271 - t274 * t208 + t296 * t210;
t174 = mrSges(7,2) * t191 - mrSges(7,3) * t186 + Ifges(7,1) * t216 + Ifges(7,4) * t215 + Ifges(7,5) * t271 + t273 * t208 - t296 * t209;
t245 = -Ifges(6,6) * qJD(3) + (Ifges(6,4) * t325 - Ifges(6,2) * t328) * qJD(1);
t344 = -mrSges(6,2) * t195 + mrSges(6,3) * t198 - Ifges(6,1) * t280 + Ifges(6,4) * t281 + pkin(8) * t166 - qJD(3) * t245 + t324 * t173 - t327 * t174;
t337 = mrSges(5,1) * t201 - mrSges(5,2) * t204 + pkin(4) * (-t374 + (-t328 * t289 - t372) * qJD(1) + t356) + qJ(5) * t177 - t344;
t147 = (-Ifges(5,3) - Ifges(4,2)) * t280 + (Ifges(4,6) - t378) * qJDD(3) - t337 + (t249 + t250) * qJD(3) + mrSges(4,3) * t220 - mrSges(4,1) * t225 - pkin(3) * t159 + (t242 + t368) * t365 + (Ifges(4,4) - Ifges(5,5)) * t281;
t247 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t328 - Ifges(4,2) * t325) * qJD(1);
t348 = -qJDD(3) * mrSges(6,2) - qJD(3) * t286 + t357;
t163 = -t281 * mrSges(6,3) - t348;
t243 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t328 + Ifges(5,3) * t325) * qJD(1);
t248 = -Ifges(6,5) * qJD(3) + (Ifges(6,1) * t325 - Ifges(6,4) * t328) * qJD(1);
t339 = mrSges(6,1) * t195 - mrSges(6,3) * t199 + Ifges(6,6) * qJDD(3) + pkin(5) * t166 + qJD(3) * t248 + t242 * t366 + t351;
t335 = mrSges(5,2) * t206 - mrSges(5,3) * t201 + Ifges(5,1) * t281 + Ifges(5,4) * qJDD(3) + Ifges(5,5) * t280 - qJ(5) * t163 + qJD(3) * t243 + t339;
t149 = t335 + Ifges(4,5) * qJDD(3) - qJD(3) * t247 - mrSges(4,3) * t219 + mrSges(4,2) * t225 - qJ(4) * t159 + t368 * t366 + (-Ifges(6,4) - Ifges(4,4)) * t280 + (Ifges(4,1) + Ifges(6,2)) * t281;
t233 = t302 + t360;
t349 = mrSges(3,2) * t238 - mrSges(3,3) * t233 + Ifges(3,1) * qJDD(1) - pkin(7) * t153 - t325 * t147 + t328 * t149;
t346 = mrSges(6,1) * t198 - mrSges(6,2) * t199 + Ifges(6,5) * t280 - Ifges(6,6) * t281 - Ifges(6,3) * qJDD(3) + pkin(5) * t353 + pkin(8) * t369 + t327 * t173 + t324 * t174 + t245 * t366 + t248 * t365;
t345 = -mrSges(3,1) * t233 - pkin(2) * t398 - pkin(7) * t154 - t328 * t147 - t325 * t149;
t334 = -m(3) * t233 + t331 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t398;
t341 = -mrSges(2,2) * t294 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t352) + qJ(2) * t334 + mrSges(2,1) * t293 + Ifges(2,3) * qJDD(1) + t349;
t336 = mrSges(5,1) * t206 - mrSges(5,3) * t204 - Ifges(5,4) * t281 - Ifges(5,2) * qJDD(3) - Ifges(5,6) * t280 + pkin(4) * t163 + t243 * t365 - t249 * t366 + t346;
t333 = -mrSges(4,2) * t220 + pkin(3) * (-t277 * t365 + (-mrSges(5,2) + mrSges(6,3)) * t281 + t348 + t394) + qJ(4) * (-t280 * mrSges(5,2) - t277 * t366 + t340) + mrSges(4,1) * t219 + t250 * t366 + t247 * t365 - Ifges(4,6) * t280 + Ifges(4,5) * t281 + Ifges(4,3) * qJDD(3) - t336;
t332 = -mrSges(3,1) * t238 - pkin(2) * t153 - t333;
t155 = m(2) * t294 - t331 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t334;
t152 = -m(3) * g(3) + t154;
t150 = m(2) * t293 - t331 * mrSges(2,2) + t383 * qJDD(1) + t352;
t146 = -t332 - qJ(2) * t152 - mrSges(2,3) * t293 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t380 * qJDD(1) + t379 * t331;
t145 = mrSges(2,3) * t294 - pkin(1) * t152 + t383 * g(3) - t379 * qJDD(1) + t380 * t331 + t345;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t329 * t146 - t326 * t145 - pkin(6) * (t329 * t150 + t326 * t155), t146, t349, t149, -t246 * t366 + t335 - t396, -Ifges(6,5) * qJDD(3) - t242 * t365 - t344, t174; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t326 * t146 + t329 * t145 + pkin(6) * (-t326 * t150 + t329 * t155), t145, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t331 * Ifges(3,5) + t332, t147, -t336, -t339 + t396, t173; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t341, t341, mrSges(3,2) * g(3) + t331 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t345, t333, t378 * qJDD(3) + t337 + Ifges(5,3) * t280 + Ifges(5,5) * t281 - qJD(3) * t249 + (-t242 + t246) * t365, t346, t351;];
m_new  = t1;
