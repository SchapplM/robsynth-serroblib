% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 19:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:01:31
% EndTime: 2019-05-05 19:02:22
% DurationCPUTime: 39.83s
% Computational Cost: add. (675332->366), mult. (1668212->460), div. (0->0), fcn. (1285310->12), ass. (0->156)
t339 = qJD(1) ^ 2;
t328 = sin(pkin(10));
t370 = qJD(1) * t328;
t330 = cos(pkin(10));
t376 = qJD(1) * t330;
t334 = sin(qJ(1));
t337 = cos(qJ(1));
t314 = -t337 * g(1) - t334 * g(2);
t307 = -t339 * pkin(1) + qJDD(1) * qJ(2) + t314;
t367 = qJD(1) * qJD(2);
t363 = -t330 * g(3) - 0.2e1 * t328 * t367;
t373 = pkin(2) * t330;
t268 = (-pkin(7) * qJDD(1) + t339 * t373 - t307) * t328 + t363;
t293 = -t328 * g(3) + (t307 + 0.2e1 * t367) * t330;
t365 = qJDD(1) * t330;
t324 = t330 ^ 2;
t371 = t324 * t339;
t275 = -pkin(2) * t371 + pkin(7) * t365 + t293;
t333 = sin(qJ(3));
t374 = cos(qJ(3));
t250 = t333 * t268 + t374 * t275;
t364 = t330 * t374;
t305 = -qJD(1) * t364 + t333 * t370;
t350 = t374 * t328 + t330 * t333;
t306 = t350 * qJD(1);
t284 = t305 * mrSges(4,1) + t306 * mrSges(4,2);
t366 = qJDD(1) * t328;
t368 = t306 * qJD(3);
t290 = -qJDD(1) * t364 + t333 * t366 + t368;
t300 = qJD(3) * mrSges(4,1) - t306 * mrSges(4,3);
t283 = t305 * pkin(3) - t306 * qJ(4);
t338 = qJD(3) ^ 2;
t235 = -t338 * pkin(3) + qJDD(3) * qJ(4) - t305 * t283 + t250;
t323 = t328 ^ 2;
t313 = t334 * g(1) - t337 * g(2);
t357 = qJDD(2) - t313;
t289 = (-pkin(1) - t373) * qJDD(1) + (-qJ(2) + (-t323 - t324) * pkin(7)) * t339 + t357;
t369 = t305 * qJD(3);
t291 = t350 * qJDD(1) - t369;
t240 = (-t291 + t369) * qJ(4) + (t290 + t368) * pkin(3) + t289;
t327 = sin(pkin(11));
t329 = cos(pkin(11));
t298 = t327 * qJD(3) + t329 * t306;
t218 = -0.2e1 * qJD(4) * t298 - t327 * t235 + t329 * t240;
t274 = t327 * qJDD(3) + t329 * t291;
t297 = t329 * qJD(3) - t327 * t306;
t208 = (t297 * t305 - t274) * pkin(8) + (t297 * t298 + t290) * pkin(4) + t218;
t219 = 0.2e1 * qJD(4) * t297 + t329 * t235 + t327 * t240;
t272 = t305 * pkin(4) - t298 * pkin(8);
t273 = t329 * qJDD(3) - t327 * t291;
t296 = t297 ^ 2;
t210 = -t296 * pkin(4) + t273 * pkin(8) - t305 * t272 + t219;
t332 = sin(qJ(5));
t336 = cos(qJ(5));
t202 = t336 * t208 - t332 * t210;
t260 = t336 * t297 - t332 * t298;
t234 = t260 * qJD(5) + t332 * t273 + t336 * t274;
t261 = t332 * t297 + t336 * t298;
t288 = qJDD(5) + t290;
t303 = qJD(5) + t305;
t199 = (t260 * t303 - t234) * pkin(9) + (t260 * t261 + t288) * pkin(5) + t202;
t203 = t332 * t208 + t336 * t210;
t233 = -t261 * qJD(5) + t336 * t273 - t332 * t274;
t253 = t303 * pkin(5) - t261 * pkin(9);
t259 = t260 ^ 2;
t200 = -t259 * pkin(5) + t233 * pkin(9) - t303 * t253 + t203;
t331 = sin(qJ(6));
t335 = cos(qJ(6));
t197 = t335 * t199 - t331 * t200;
t245 = t335 * t260 - t331 * t261;
t216 = t245 * qJD(6) + t331 * t233 + t335 * t234;
t246 = t331 * t260 + t335 * t261;
t226 = -t245 * mrSges(7,1) + t246 * mrSges(7,2);
t302 = qJD(6) + t303;
t238 = -t302 * mrSges(7,2) + t245 * mrSges(7,3);
t282 = qJDD(6) + t288;
t190 = m(7) * t197 + t282 * mrSges(7,1) - t216 * mrSges(7,3) - t246 * t226 + t302 * t238;
t198 = t331 * t199 + t335 * t200;
t215 = -t246 * qJD(6) + t335 * t233 - t331 * t234;
t239 = t302 * mrSges(7,1) - t246 * mrSges(7,3);
t191 = m(7) * t198 - t282 * mrSges(7,2) + t215 * mrSges(7,3) + t245 * t226 - t302 * t239;
t184 = t335 * t190 + t331 * t191;
t247 = -t260 * mrSges(6,1) + t261 * mrSges(6,2);
t251 = -t303 * mrSges(6,2) + t260 * mrSges(6,3);
t181 = m(6) * t202 + t288 * mrSges(6,1) - t234 * mrSges(6,3) - t261 * t247 + t303 * t251 + t184;
t252 = t303 * mrSges(6,1) - t261 * mrSges(6,3);
t358 = -t331 * t190 + t335 * t191;
t182 = m(6) * t203 - t288 * mrSges(6,2) + t233 * mrSges(6,3) + t260 * t247 - t303 * t252 + t358;
t177 = t336 * t181 + t332 * t182;
t263 = -t297 * mrSges(5,1) + t298 * mrSges(5,2);
t270 = -t305 * mrSges(5,2) + t297 * mrSges(5,3);
t175 = m(5) * t218 + t290 * mrSges(5,1) - t274 * mrSges(5,3) - t298 * t263 + t305 * t270 + t177;
t271 = t305 * mrSges(5,1) - t298 * mrSges(5,3);
t359 = -t332 * t181 + t336 * t182;
t176 = m(5) * t219 - t290 * mrSges(5,2) + t273 * mrSges(5,3) + t297 * t263 - t305 * t271 + t359;
t360 = -t327 * t175 + t329 * t176;
t166 = m(4) * t250 - qJDD(3) * mrSges(4,2) - t290 * mrSges(4,3) - qJD(3) * t300 - t305 * t284 + t360;
t249 = t374 * t268 - t333 * t275;
t299 = -qJD(3) * mrSges(4,2) - t305 * mrSges(4,3);
t232 = -qJDD(3) * pkin(3) - t338 * qJ(4) + t306 * t283 + qJDD(4) - t249;
t220 = -t273 * pkin(4) - t296 * pkin(8) + t298 * t272 + t232;
t205 = -t233 * pkin(5) - t259 * pkin(9) + t261 * t253 + t220;
t353 = m(7) * t205 - t215 * mrSges(7,1) + t216 * mrSges(7,2) - t245 * t238 + t246 * t239;
t345 = m(6) * t220 - t233 * mrSges(6,1) + t234 * mrSges(6,2) - t260 * t251 + t261 * t252 + t353;
t341 = -m(5) * t232 + t273 * mrSges(5,1) - t274 * mrSges(5,2) + t297 * t270 - t298 * t271 - t345;
t193 = m(4) * t249 + qJDD(3) * mrSges(4,1) - t291 * mrSges(4,3) + qJD(3) * t299 - t306 * t284 + t341;
t161 = t333 * t166 + t374 * t193;
t292 = -t328 * t307 + t363;
t221 = Ifges(7,5) * t246 + Ifges(7,6) * t245 + Ifges(7,3) * t302;
t223 = Ifges(7,1) * t246 + Ifges(7,4) * t245 + Ifges(7,5) * t302;
t185 = -mrSges(7,1) * t205 + mrSges(7,3) * t198 + Ifges(7,4) * t216 + Ifges(7,2) * t215 + Ifges(7,6) * t282 - t246 * t221 + t302 * t223;
t222 = Ifges(7,4) * t246 + Ifges(7,2) * t245 + Ifges(7,6) * t302;
t186 = mrSges(7,2) * t205 - mrSges(7,3) * t197 + Ifges(7,1) * t216 + Ifges(7,4) * t215 + Ifges(7,5) * t282 + t245 * t221 - t302 * t222;
t241 = Ifges(6,5) * t261 + Ifges(6,6) * t260 + Ifges(6,3) * t303;
t243 = Ifges(6,1) * t261 + Ifges(6,4) * t260 + Ifges(6,5) * t303;
t170 = -mrSges(6,1) * t220 + mrSges(6,3) * t203 + Ifges(6,4) * t234 + Ifges(6,2) * t233 + Ifges(6,6) * t288 - pkin(5) * t353 + pkin(9) * t358 + t335 * t185 + t331 * t186 - t261 * t241 + t303 * t243;
t242 = Ifges(6,4) * t261 + Ifges(6,2) * t260 + Ifges(6,6) * t303;
t171 = mrSges(6,2) * t220 - mrSges(6,3) * t202 + Ifges(6,1) * t234 + Ifges(6,4) * t233 + Ifges(6,5) * t288 - pkin(9) * t184 - t331 * t185 + t335 * t186 + t260 * t241 - t303 * t242;
t254 = Ifges(5,5) * t298 + Ifges(5,6) * t297 + Ifges(5,3) * t305;
t256 = Ifges(5,1) * t298 + Ifges(5,4) * t297 + Ifges(5,5) * t305;
t155 = -mrSges(5,1) * t232 + mrSges(5,3) * t219 + Ifges(5,4) * t274 + Ifges(5,2) * t273 + Ifges(5,6) * t290 - pkin(4) * t345 + pkin(8) * t359 + t336 * t170 + t332 * t171 - t298 * t254 + t305 * t256;
t255 = Ifges(5,4) * t298 + Ifges(5,2) * t297 + Ifges(5,6) * t305;
t157 = mrSges(5,2) * t232 - mrSges(5,3) * t218 + Ifges(5,1) * t274 + Ifges(5,4) * t273 + Ifges(5,5) * t290 - pkin(8) * t177 - t332 * t170 + t336 * t171 + t297 * t254 - t305 * t255;
t277 = Ifges(4,4) * t306 - Ifges(4,2) * t305 + Ifges(4,6) * qJD(3);
t278 = Ifges(4,1) * t306 - Ifges(4,4) * t305 + Ifges(4,5) * qJD(3);
t346 = -mrSges(4,1) * t249 + mrSges(4,2) * t250 - Ifges(4,5) * t291 + Ifges(4,6) * t290 - Ifges(4,3) * qJDD(3) - pkin(3) * t341 - qJ(4) * t360 - t329 * t155 - t327 * t157 - t306 * t277 - t305 * t278;
t355 = Ifges(3,4) * t328 + Ifges(3,2) * t330;
t356 = Ifges(3,1) * t328 + Ifges(3,4) * t330;
t375 = -mrSges(3,1) * t292 + mrSges(3,2) * t293 - pkin(2) * t161 - (t355 * t370 - t356 * t376) * qJD(1) + t346;
t372 = mrSges(3,2) * t328;
t168 = t329 * t175 + t327 * t176;
t351 = mrSges(3,3) * qJDD(1) + t339 * (-mrSges(3,1) * t330 + t372);
t159 = m(3) * t292 - t351 * t328 + t161;
t361 = t374 * t166 - t333 * t193;
t160 = m(3) * t293 + t351 * t330 + t361;
t362 = -t328 * t159 + t330 * t160;
t354 = Ifges(3,5) * t328 + Ifges(3,6) * t330;
t276 = Ifges(4,5) * t306 - Ifges(4,6) * t305 + Ifges(4,3) * qJD(3);
t149 = mrSges(4,2) * t289 - mrSges(4,3) * t249 + Ifges(4,1) * t291 - Ifges(4,4) * t290 + Ifges(4,5) * qJDD(3) - qJ(4) * t168 - qJD(3) * t277 - t327 * t155 + t329 * t157 - t305 * t276;
t348 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t216 - Ifges(7,6) * t215 - Ifges(7,3) * t282 - t246 * t222 + t245 * t223;
t343 = -mrSges(6,1) * t202 + mrSges(6,2) * t203 - Ifges(6,5) * t234 - Ifges(6,6) * t233 - Ifges(6,3) * t288 - pkin(5) * t184 - t261 * t242 + t260 * t243 + t348;
t340 = mrSges(5,1) * t218 - mrSges(5,2) * t219 + Ifges(5,5) * t274 + Ifges(5,6) * t273 + pkin(4) * t177 + t298 * t255 - t297 * t256 - t343;
t153 = (-Ifges(5,3) - Ifges(4,2)) * t290 - t306 * t276 - mrSges(4,1) * t289 + Ifges(4,4) * t291 + qJD(3) * t278 + mrSges(4,3) * t250 - t340 + Ifges(4,6) * qJDD(3) - pkin(3) * t168;
t304 = -qJDD(1) * pkin(1) - t339 * qJ(2) + t357;
t309 = t354 * qJD(1);
t347 = m(4) * t289 + t290 * mrSges(4,1) + t291 * mrSges(4,2) + t305 * t299 + t306 * t300 + t168;
t145 = -mrSges(3,1) * t304 + mrSges(3,3) * t293 - pkin(2) * t347 + pkin(7) * t361 + t355 * qJDD(1) + t333 * t149 + t374 * t153 - t309 * t370;
t148 = mrSges(3,2) * t304 - mrSges(3,3) * t292 - pkin(7) * t161 + t356 * qJDD(1) + t374 * t149 - t333 * t153 + t309 * t376;
t344 = -m(3) * t304 + mrSges(3,1) * t365 - t347 + (t323 * t339 + t371) * mrSges(3,3);
t349 = -mrSges(2,2) * t314 + qJ(2) * t362 + t330 * t145 + t328 * t148 + pkin(1) * (-mrSges(3,2) * t366 + t344) + mrSges(2,1) * t313 + Ifges(2,3) * qJDD(1);
t162 = (mrSges(2,1) - t372) * qJDD(1) + t344 - t339 * mrSges(2,2) + m(2) * t313;
t152 = t330 * t159 + t328 * t160;
t150 = m(2) * t314 - t339 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t362;
t146 = mrSges(2,1) * g(3) + (Ifges(2,6) - t354) * qJDD(1) + t339 * Ifges(2,5) + mrSges(2,3) * t314 - pkin(1) * t152 + t375;
t143 = -mrSges(2,2) * g(3) - mrSges(2,3) * t313 + Ifges(2,5) * qJDD(1) - t339 * Ifges(2,6) - qJ(2) * t152 - t328 * t145 + t330 * t148;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t337 * t143 - t334 * t146 - pkin(6) * (t334 * t150 + t337 * t162), t143, t148, t149, t157, t171, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t334 * t143 + t337 * t146 + pkin(6) * (t337 * t150 - t334 * t162), t146, t145, t153, t155, t170, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t349, t349, t354 * qJDD(1) - t375, -t346, Ifges(5,3) * t290 + t340, -t343, -t348;];
m_new  = t1;
