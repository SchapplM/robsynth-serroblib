% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 08:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:38:58
% EndTime: 2019-05-07 08:39:18
% DurationCPUTime: 7.00s
% Computational Cost: add. (116865->377), mult. (229420->443), div. (0->0), fcn. (151376->8), ass. (0->136)
t351 = sin(qJ(1));
t353 = cos(qJ(1));
t334 = t351 * g(1) - t353 * g(2);
t355 = qJD(1) ^ 2;
t312 = -qJDD(1) * pkin(1) - t355 * pkin(7) - t334;
t350 = sin(qJ(2));
t352 = cos(qJ(2));
t376 = qJD(1) * qJD(2);
t373 = t352 * t376;
t328 = t350 * qJDD(1) + t373;
t374 = t350 * t376;
t329 = t352 * qJDD(1) - t374;
t249 = (-t328 - t373) * pkin(8) + (-t329 + t374) * pkin(2) + t312;
t335 = -t353 * g(1) - t351 * g(2);
t313 = -t355 * pkin(1) + qJDD(1) * pkin(7) + t335;
t296 = -t350 * g(3) + t352 * t313;
t327 = (-pkin(2) * t352 - pkin(8) * t350) * qJD(1);
t354 = qJD(2) ^ 2;
t377 = t352 * qJD(1);
t256 = -t354 * pkin(2) + qJDD(2) * pkin(8) + t327 * t377 + t296;
t349 = sin(qJ(3));
t389 = cos(qJ(3));
t232 = t389 * t249 - t349 * t256;
t233 = t349 * t249 + t389 * t256;
t378 = qJD(1) * t350;
t324 = -t389 * qJD(2) + t349 * t378;
t325 = t349 * qJD(2) + t389 * t378;
t339 = qJD(3) - t377;
t265 = Ifges(5,5) * t325 + Ifges(5,6) * t339 + Ifges(5,3) * t324;
t268 = Ifges(4,4) * t325 - Ifges(4,2) * t324 + Ifges(4,6) * t339;
t270 = Ifges(4,1) * t325 - Ifges(4,4) * t324 + Ifges(4,5) * t339;
t282 = t325 * qJD(3) - t389 * qJDD(2) + t349 * t328;
t283 = -t324 * qJD(3) + t349 * qJDD(2) + t389 * t328;
t289 = t324 * mrSges(5,1) - t325 * mrSges(5,3);
t323 = qJDD(3) - t329;
t288 = t324 * pkin(3) - t325 * qJ(4);
t338 = t339 ^ 2;
t215 = -t323 * pkin(3) - t338 * qJ(4) + t325 * t288 + qJDD(4) - t232;
t384 = t324 * t339;
t205 = (-t283 - t384) * pkin(9) + (t324 * t325 - t323) * pkin(4) + t215;
t390 = 2 * qJD(4);
t212 = -t338 * pkin(3) + t323 * qJ(4) - t324 * t288 + t339 * t390 + t233;
t297 = -t339 * pkin(4) - t325 * pkin(9);
t322 = t324 ^ 2;
t208 = -t322 * pkin(4) + t282 * pkin(9) + t339 * t297 + t212;
t348 = sin(qJ(5));
t388 = cos(qJ(5));
t203 = t348 * t205 + t388 * t208;
t285 = t348 * t324 + t388 * t325;
t230 = t285 * qJD(5) - t388 * t282 + t348 * t283;
t337 = qJD(5) - t339;
t260 = t337 * mrSges(6,1) - t285 * mrSges(6,3);
t284 = -t388 * t324 + t348 * t325;
t319 = qJDD(5) - t323;
t245 = t284 * pkin(5) - t285 * qJ(6);
t336 = t337 ^ 2;
t194 = -t336 * pkin(5) + t319 * qJ(6) + 0.2e1 * qJD(6) * t337 - t284 * t245 + t203;
t261 = -t337 * mrSges(7,1) + t285 * mrSges(7,2);
t375 = m(7) * t194 + t319 * mrSges(7,3) + t337 * t261;
t246 = t284 * mrSges(7,1) - t285 * mrSges(7,3);
t382 = -t284 * mrSges(6,1) - t285 * mrSges(6,2) - t246;
t385 = -mrSges(6,3) - mrSges(7,2);
t184 = m(6) * t203 - t319 * mrSges(6,2) + t385 * t230 - t337 * t260 + t382 * t284 + t375;
t202 = t388 * t205 - t348 * t208;
t231 = -t284 * qJD(5) + t348 * t282 + t388 * t283;
t259 = -t337 * mrSges(6,2) - t284 * mrSges(6,3);
t197 = -t319 * pkin(5) - t336 * qJ(6) + t285 * t245 + qJDD(6) - t202;
t258 = -t284 * mrSges(7,2) + t337 * mrSges(7,3);
t371 = -m(7) * t197 + t319 * mrSges(7,1) + t337 * t258;
t185 = m(6) * t202 + t319 * mrSges(6,1) + t385 * t231 + t337 * t259 + t382 * t285 + t371;
t178 = t348 * t184 + t388 * t185;
t269 = Ifges(5,1) * t325 + Ifges(5,4) * t339 + Ifges(5,5) * t324;
t238 = Ifges(6,4) * t285 - Ifges(6,2) * t284 + Ifges(6,6) * t337;
t240 = Ifges(6,1) * t285 - Ifges(6,4) * t284 + Ifges(6,5) * t337;
t235 = Ifges(7,5) * t285 + Ifges(7,6) * t337 + Ifges(7,3) * t284;
t239 = Ifges(7,1) * t285 + Ifges(7,4) * t337 + Ifges(7,5) * t284;
t367 = -mrSges(7,1) * t197 + mrSges(7,3) * t194 + Ifges(7,4) * t231 + Ifges(7,2) * t319 + Ifges(7,6) * t230 - t285 * t235 + t284 * t239;
t360 = Ifges(6,5) * t231 - Ifges(6,6) * t230 + Ifges(6,3) * t319 + t285 * t238 + t284 * t240 + mrSges(6,1) * t202 - mrSges(6,2) * t203 + pkin(5) * (-t231 * mrSges(7,2) - t285 * t246 + t371) + qJ(6) * (-t230 * mrSges(7,2) - t284 * t246 + t375) + t367;
t357 = mrSges(5,1) * t215 - mrSges(5,3) * t212 - Ifges(5,4) * t283 - Ifges(5,2) * t323 - Ifges(5,6) * t282 + pkin(4) * t178 - t324 * t269 + t360;
t294 = -t324 * mrSges(5,2) + t339 * mrSges(5,3);
t363 = -m(5) * t215 + t323 * mrSges(5,1) + t339 * t294 - t178;
t179 = t388 * t184 - t348 * t185;
t293 = -t339 * mrSges(5,1) + t325 * mrSges(5,2);
t368 = m(5) * t212 + t323 * mrSges(5,3) + t339 * t293 + t179;
t392 = -(t265 - t268) * t325 + mrSges(4,1) * t232 - mrSges(4,2) * t233 + Ifges(4,5) * t283 - Ifges(4,6) * t282 + Ifges(4,3) * t323 + pkin(3) * (-t283 * mrSges(5,2) - t325 * t289 + t363) + qJ(4) * (-t282 * mrSges(5,2) - t324 * t289 + t368) + t324 * t270 - t357;
t295 = -t352 * g(3) - t350 * t313;
t255 = -qJDD(2) * pkin(2) - t354 * pkin(8) + t327 * t378 - t295;
t365 = t282 * pkin(3) + t255 + (-t283 + t384) * qJ(4);
t387 = pkin(3) * t339;
t209 = -t282 * pkin(4) - t322 * pkin(9) - t365 + (t297 - t387 + t390) * t325;
t199 = -0.2e1 * qJD(6) * t285 + (t285 * t337 + t230) * pkin(5) + (t284 * t337 - t231) * qJ(6) + t209;
t187 = m(7) * t199 + t230 * mrSges(7,1) - t231 * mrSges(7,3) + t284 * t258 - t285 * t261;
t186 = -m(6) * t209 - t230 * mrSges(6,1) - t231 * mrSges(6,2) - t284 * t259 - t285 * t260 - t187;
t214 = (-(2 * qJD(4)) + t387) * t325 + t365;
t182 = m(5) * t214 + t282 * mrSges(5,1) - t283 * mrSges(5,3) - t325 * t293 + t324 * t294 + t186;
t370 = -mrSges(7,1) * t199 + mrSges(7,2) * t194;
t237 = Ifges(7,4) * t285 + Ifges(7,2) * t337 + Ifges(7,6) * t284;
t383 = -Ifges(6,5) * t285 + Ifges(6,6) * t284 - Ifges(6,3) * t337 - t237;
t173 = -mrSges(6,1) * t209 + mrSges(6,3) * t203 - pkin(5) * t187 + (t239 + t240) * t337 + (Ifges(6,6) - Ifges(7,6)) * t319 + t383 * t285 + (Ifges(6,4) - Ifges(7,5)) * t231 + (-Ifges(6,2) - Ifges(7,3)) * t230 + t370;
t366 = mrSges(7,2) * t197 - mrSges(7,3) * t199 + Ifges(7,1) * t231 + Ifges(7,4) * t319 + Ifges(7,5) * t230 + t337 * t235;
t175 = mrSges(6,2) * t209 - mrSges(6,3) * t202 + Ifges(6,1) * t231 - Ifges(6,4) * t230 + Ifges(6,5) * t319 - qJ(6) * t187 - t337 * t238 + t383 * t284 + t366;
t361 = mrSges(5,1) * t214 - mrSges(5,2) * t212 + pkin(4) * t186 + pkin(9) * t179 + t388 * t173 + t348 * t175;
t267 = Ifges(5,4) * t325 + Ifges(5,2) * t339 + Ifges(5,6) * t324;
t380 = -Ifges(4,5) * t325 + Ifges(4,6) * t324 - Ifges(4,3) * t339 - t267;
t159 = -mrSges(4,1) * t255 + mrSges(4,3) * t233 - pkin(3) * t182 + (t270 + t269) * t339 + t380 * t325 + (Ifges(4,6) - Ifges(5,6)) * t323 + (Ifges(4,4) - Ifges(5,5)) * t283 + (-Ifges(4,2) - Ifges(5,3)) * t282 - t361;
t362 = mrSges(5,2) * t215 - mrSges(5,3) * t214 + Ifges(5,1) * t283 + Ifges(5,4) * t323 + Ifges(5,5) * t282 - pkin(9) * t178 - t348 * t173 + t388 * t175 + t339 * t265;
t160 = mrSges(4,2) * t255 - mrSges(4,3) * t232 + Ifges(4,1) * t283 - Ifges(4,4) * t282 + Ifges(4,5) * t323 - qJ(4) * t182 - t339 * t268 + t380 * t324 + t362;
t292 = t339 * mrSges(4,1) - t325 * mrSges(4,3);
t379 = -t324 * mrSges(4,1) - t325 * mrSges(4,2) - t289;
t386 = -mrSges(4,3) - mrSges(5,2);
t171 = m(4) * t233 - t323 * mrSges(4,2) + t386 * t282 - t339 * t292 + t379 * t324 + t368;
t291 = -t339 * mrSges(4,2) - t324 * mrSges(4,3);
t172 = m(4) * t232 + t323 * mrSges(4,1) + t386 * t283 + t339 * t291 + t379 * t325 + t363;
t169 = t389 * t171 - t349 * t172;
t181 = -m(4) * t255 - t282 * mrSges(4,1) - t283 * mrSges(4,2) - t324 * t291 - t325 * t292 - t182;
t310 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t350 + Ifges(3,2) * t352) * qJD(1);
t311 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t350 + Ifges(3,4) * t352) * qJD(1);
t391 = mrSges(3,1) * t295 - mrSges(3,2) * t296 + Ifges(3,5) * t328 + Ifges(3,6) * t329 + Ifges(3,3) * qJDD(2) + pkin(2) * t181 + pkin(8) * t169 + (t310 * t350 - t311 * t352) * qJD(1) + t389 * t159 + t349 * t160;
t326 = (-mrSges(3,1) * t352 + mrSges(3,2) * t350) * qJD(1);
t331 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t378;
t167 = m(3) * t296 - qJDD(2) * mrSges(3,2) + t329 * mrSges(3,3) - qJD(2) * t331 + t326 * t377 + t169;
t332 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t377;
t180 = m(3) * t295 + qJDD(2) * mrSges(3,1) - t328 * mrSges(3,3) + qJD(2) * t332 - t326 * t378 + t181;
t372 = t352 * t167 - t350 * t180;
t168 = t349 * t171 + t389 * t172;
t309 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t350 + Ifges(3,6) * t352) * qJD(1);
t156 = mrSges(3,2) * t312 - mrSges(3,3) * t295 + Ifges(3,1) * t328 + Ifges(3,4) * t329 + Ifges(3,5) * qJDD(2) - pkin(8) * t168 - qJD(2) * t310 - t349 * t159 + t389 * t160 + t309 * t377;
t158 = -mrSges(3,1) * t312 + mrSges(3,3) * t296 + Ifges(3,4) * t328 + Ifges(3,2) * t329 + Ifges(3,6) * qJDD(2) - pkin(2) * t168 + qJD(2) * t311 - t309 * t378 - t392;
t359 = -m(3) * t312 + t329 * mrSges(3,1) - t328 * mrSges(3,2) - t331 * t378 + t332 * t377 - t168;
t364 = mrSges(2,1) * t334 - mrSges(2,2) * t335 + Ifges(2,3) * qJDD(1) + pkin(1) * t359 + pkin(7) * t372 + t350 * t156 + t352 * t158;
t164 = m(2) * t334 + qJDD(1) * mrSges(2,1) - t355 * mrSges(2,2) + t359;
t163 = t350 * t167 + t352 * t180;
t161 = m(2) * t335 - t355 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t372;
t154 = mrSges(2,1) * g(3) + mrSges(2,3) * t335 + t355 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t163 - t391;
t153 = -mrSges(2,2) * g(3) - mrSges(2,3) * t334 + Ifges(2,5) * qJDD(1) - t355 * Ifges(2,6) - pkin(7) * t163 + t352 * t156 - t350 * t158;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t353 * t153 - t351 * t154 - pkin(6) * (t351 * t161 + t353 * t164), t153, t156, t160, -t324 * t267 + t362, t175, -t284 * t237 + t366; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t351 * t153 + t353 * t154 + pkin(6) * (t353 * t161 - t351 * t164), t154, t158, t159, -t325 * t265 - t357, t173, t367; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t364, t364, t391, t392, Ifges(5,5) * t283 + Ifges(5,6) * t323 + Ifges(5,3) * t282 + t325 * t267 - t339 * t269 + t361, t360, Ifges(7,5) * t231 + Ifges(7,6) * t319 + Ifges(7,3) * t230 + t285 * t237 - t337 * t239 - t370;];
m_new  = t1;
