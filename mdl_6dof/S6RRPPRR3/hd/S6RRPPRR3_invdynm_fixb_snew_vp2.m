% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-06 10:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:59:40
% EndTime: 2019-05-06 10:01:20
% DurationCPUTime: 66.06s
% Computational Cost: add. (1081565->396), mult. (2921744->518), div. (0->0), fcn. (2330153->14), ass. (0->163)
t392 = -2 * qJD(3);
t352 = sin(pkin(6));
t357 = sin(qJ(2));
t361 = cos(qJ(2));
t382 = qJD(1) * qJD(2);
t336 = (qJDD(1) * t357 + t361 * t382) * t352;
t354 = cos(pkin(6));
t344 = t354 * qJDD(1) + qJDD(2);
t345 = t354 * qJD(1) + qJD(2);
t358 = sin(qJ(1));
t362 = cos(qJ(1));
t341 = t358 * g(1) - t362 * g(2);
t363 = qJD(1) ^ 2;
t391 = pkin(8) * t352;
t333 = qJDD(1) * pkin(1) + t363 * t391 + t341;
t342 = -t362 * g(1) - t358 * g(2);
t334 = -t363 * pkin(1) + qJDD(1) * t391 + t342;
t385 = t354 * t361;
t374 = t333 * t385 - t357 * t334;
t389 = t352 ^ 2 * t363;
t270 = t344 * pkin(2) - t336 * qJ(3) + (pkin(2) * t357 * t389 + (qJ(3) * qJD(1) * t345 - g(3)) * t352) * t361 + t374;
t386 = t354 * t357;
t388 = t352 * t357;
t298 = -g(3) * t388 + t333 * t386 + t361 * t334;
t384 = qJD(1) * t352;
t380 = t357 * t384;
t330 = t345 * pkin(2) - qJ(3) * t380;
t337 = (qJDD(1) * t361 - t357 * t382) * t352;
t381 = t361 ^ 2 * t389;
t273 = -pkin(2) * t381 + t337 * qJ(3) - t345 * t330 + t298;
t351 = sin(pkin(11));
t390 = cos(pkin(11));
t327 = (t351 * t361 + t357 * t390) * t384;
t245 = t270 * t390 - t351 * t273 + t327 * t392;
t387 = t352 * t361;
t379 = t361 * t384;
t326 = t351 * t380 - t390 * t379;
t246 = t351 * t270 + t390 * t273 + t326 * t392;
t300 = t326 * mrSges(4,1) + t327 * mrSges(4,2);
t306 = t351 * t336 - t390 * t337;
t314 = t345 * mrSges(4,1) - t327 * mrSges(4,3);
t299 = t326 * pkin(3) - t327 * qJ(4);
t343 = t345 ^ 2;
t238 = -t343 * pkin(3) + t344 * qJ(4) - t326 * t299 + t246;
t318 = -t354 * g(3) - t352 * t333;
t281 = -t337 * pkin(2) - qJ(3) * t381 + t330 * t380 + qJDD(3) + t318;
t307 = t336 * t390 + t351 * t337;
t249 = (t326 * t345 - t307) * qJ(4) + (t327 * t345 + t306) * pkin(3) + t281;
t350 = sin(pkin(12));
t353 = cos(pkin(12));
t312 = t353 * t327 + t350 * t345;
t230 = -0.2e1 * qJD(4) * t312 - t350 * t238 + t353 * t249;
t292 = t353 * t307 + t350 * t344;
t311 = -t350 * t327 + t353 * t345;
t227 = (t311 * t326 - t292) * pkin(9) + (t311 * t312 + t306) * pkin(4) + t230;
t231 = 0.2e1 * qJD(4) * t311 + t353 * t238 + t350 * t249;
t289 = t326 * pkin(4) - t312 * pkin(9);
t291 = -t350 * t307 + t353 * t344;
t310 = t311 ^ 2;
t229 = -t310 * pkin(4) + t291 * pkin(9) - t326 * t289 + t231;
t356 = sin(qJ(5));
t360 = cos(qJ(5));
t224 = t356 * t227 + t360 * t229;
t283 = t356 * t311 + t360 * t312;
t255 = -t283 * qJD(5) + t360 * t291 - t356 * t292;
t282 = t360 * t311 - t356 * t312;
t263 = -t282 * mrSges(6,1) + t283 * mrSges(6,2);
t325 = qJD(5) + t326;
t272 = t325 * mrSges(6,1) - t283 * mrSges(6,3);
t305 = qJDD(5) + t306;
t264 = -t282 * pkin(5) - t283 * pkin(10);
t324 = t325 ^ 2;
t221 = -t324 * pkin(5) + t305 * pkin(10) + t282 * t264 + t224;
t237 = -t344 * pkin(3) - t343 * qJ(4) + t327 * t299 + qJDD(4) - t245;
t232 = -t291 * pkin(4) - t310 * pkin(9) + t312 * t289 + t237;
t256 = t282 * qJD(5) + t356 * t291 + t360 * t292;
t225 = (-t282 * t325 - t256) * pkin(10) + (t283 * t325 - t255) * pkin(5) + t232;
t355 = sin(qJ(6));
t359 = cos(qJ(6));
t218 = -t355 * t221 + t359 * t225;
t266 = -t355 * t283 + t359 * t325;
t235 = t266 * qJD(6) + t359 * t256 + t355 * t305;
t267 = t359 * t283 + t355 * t325;
t250 = -t266 * mrSges(7,1) + t267 * mrSges(7,2);
t254 = qJDD(6) - t255;
t280 = qJD(6) - t282;
t257 = -t280 * mrSges(7,2) + t266 * mrSges(7,3);
t214 = m(7) * t218 + t254 * mrSges(7,1) - t235 * mrSges(7,3) - t267 * t250 + t280 * t257;
t219 = t359 * t221 + t355 * t225;
t234 = -t267 * qJD(6) - t355 * t256 + t359 * t305;
t258 = t280 * mrSges(7,1) - t267 * mrSges(7,3);
t215 = m(7) * t219 - t254 * mrSges(7,2) + t234 * mrSges(7,3) + t266 * t250 - t280 * t258;
t375 = -t355 * t214 + t359 * t215;
t200 = m(6) * t224 - t305 * mrSges(6,2) + t255 * mrSges(6,3) + t282 * t263 - t325 * t272 + t375;
t223 = t360 * t227 - t356 * t229;
t271 = -t325 * mrSges(6,2) + t282 * mrSges(6,3);
t220 = -t305 * pkin(5) - t324 * pkin(10) + t283 * t264 - t223;
t372 = -m(7) * t220 + t234 * mrSges(7,1) - t235 * mrSges(7,2) + t266 * t257 - t267 * t258;
t210 = m(6) * t223 + t305 * mrSges(6,1) - t256 * mrSges(6,3) - t283 * t263 + t325 * t271 + t372;
t195 = t356 * t200 + t360 * t210;
t285 = -t311 * mrSges(5,1) + t312 * mrSges(5,2);
t287 = -t326 * mrSges(5,2) + t311 * mrSges(5,3);
t193 = m(5) * t230 + t306 * mrSges(5,1) - t292 * mrSges(5,3) - t312 * t285 + t326 * t287 + t195;
t288 = t326 * mrSges(5,1) - t312 * mrSges(5,3);
t376 = t360 * t200 - t356 * t210;
t194 = m(5) * t231 - t306 * mrSges(5,2) + t291 * mrSges(5,3) + t311 * t285 - t326 * t288 + t376;
t377 = -t350 * t193 + t353 * t194;
t184 = m(4) * t246 - t344 * mrSges(4,2) - t306 * mrSges(4,3) - t326 * t300 - t345 * t314 + t377;
t313 = -t345 * mrSges(4,2) - t326 * mrSges(4,3);
t203 = t359 * t214 + t355 * t215;
t369 = m(6) * t232 - t255 * mrSges(6,1) + t256 * mrSges(6,2) - t282 * t271 + t283 * t272 + t203;
t365 = -m(5) * t237 + t291 * mrSges(5,1) - t292 * mrSges(5,2) + t311 * t287 - t312 * t288 - t369;
t197 = m(4) * t245 + t344 * mrSges(4,1) - t307 * mrSges(4,3) - t327 * t300 + t345 * t313 + t365;
t181 = t351 * t184 + t390 * t197;
t187 = t353 * t193 + t350 * t194;
t297 = -g(3) * t387 + t374;
t332 = -t345 * mrSges(3,2) + mrSges(3,3) * t379;
t335 = (-mrSges(3,1) * t361 + mrSges(3,2) * t357) * t384;
t179 = m(3) * t297 + t344 * mrSges(3,1) - t336 * mrSges(3,3) + t345 * t332 - t335 * t380 + t181;
t331 = t345 * mrSges(3,1) - mrSges(3,3) * t380;
t378 = t390 * t184 - t351 * t197;
t180 = m(3) * t298 - t344 * mrSges(3,2) + t337 * mrSges(3,3) - t345 * t331 + t335 * t379 + t378;
t170 = -t357 * t179 + t361 * t180;
t370 = m(4) * t281 + t306 * mrSges(4,1) + t307 * mrSges(4,2) + t326 * t313 + t327 * t314 + t187;
t185 = m(3) * t318 - t337 * mrSges(3,1) + t336 * mrSges(3,2) + (t331 * t357 - t332 * t361) * t384 + t370;
t167 = t179 * t385 + t180 * t386 - t352 * t185;
t239 = Ifges(7,5) * t267 + Ifges(7,6) * t266 + Ifges(7,3) * t280;
t241 = Ifges(7,1) * t267 + Ifges(7,4) * t266 + Ifges(7,5) * t280;
t207 = -mrSges(7,1) * t220 + mrSges(7,3) * t219 + Ifges(7,4) * t235 + Ifges(7,2) * t234 + Ifges(7,6) * t254 - t267 * t239 + t280 * t241;
t240 = Ifges(7,4) * t267 + Ifges(7,2) * t266 + Ifges(7,6) * t280;
t208 = mrSges(7,2) * t220 - mrSges(7,3) * t218 + Ifges(7,1) * t235 + Ifges(7,4) * t234 + Ifges(7,5) * t254 + t266 * t239 - t280 * t240;
t259 = Ifges(6,5) * t283 + Ifges(6,6) * t282 + Ifges(6,3) * t325;
t260 = Ifges(6,4) * t283 + Ifges(6,2) * t282 + Ifges(6,6) * t325;
t188 = mrSges(6,2) * t232 - mrSges(6,3) * t223 + Ifges(6,1) * t256 + Ifges(6,4) * t255 + Ifges(6,5) * t305 - pkin(10) * t203 - t355 * t207 + t359 * t208 + t282 * t259 - t325 * t260;
t261 = Ifges(6,1) * t283 + Ifges(6,4) * t282 + Ifges(6,5) * t325;
t366 = mrSges(7,1) * t218 - mrSges(7,2) * t219 + Ifges(7,5) * t235 + Ifges(7,6) * t234 + Ifges(7,3) * t254 + t267 * t240 - t266 * t241;
t189 = -mrSges(6,1) * t232 + mrSges(6,3) * t224 + Ifges(6,4) * t256 + Ifges(6,2) * t255 + Ifges(6,6) * t305 - pkin(5) * t203 - t283 * t259 + t325 * t261 - t366;
t274 = Ifges(5,5) * t312 + Ifges(5,6) * t311 + Ifges(5,3) * t326;
t276 = Ifges(5,1) * t312 + Ifges(5,4) * t311 + Ifges(5,5) * t326;
t173 = -mrSges(5,1) * t237 + mrSges(5,3) * t231 + Ifges(5,4) * t292 + Ifges(5,2) * t291 + Ifges(5,6) * t306 - pkin(4) * t369 + pkin(9) * t376 + t356 * t188 + t360 * t189 - t312 * t274 + t326 * t276;
t275 = Ifges(5,4) * t312 + Ifges(5,2) * t311 + Ifges(5,6) * t326;
t175 = mrSges(5,2) * t237 - mrSges(5,3) * t230 + Ifges(5,1) * t292 + Ifges(5,4) * t291 + Ifges(5,5) * t306 - pkin(9) * t195 + t360 * t188 - t356 * t189 + t311 * t274 - t326 * t275;
t293 = Ifges(4,5) * t327 - Ifges(4,6) * t326 + Ifges(4,3) * t345;
t294 = Ifges(4,4) * t327 - Ifges(4,2) * t326 + Ifges(4,6) * t345;
t163 = mrSges(4,2) * t281 - mrSges(4,3) * t245 + Ifges(4,1) * t307 - Ifges(4,4) * t306 + Ifges(4,5) * t344 - qJ(4) * t187 - t350 * t173 + t353 * t175 - t326 * t293 - t345 * t294;
t295 = Ifges(4,1) * t327 - Ifges(4,4) * t326 + Ifges(4,5) * t345;
t367 = -mrSges(6,1) * t223 + mrSges(6,2) * t224 - Ifges(6,5) * t256 - Ifges(6,6) * t255 - Ifges(6,3) * t305 - pkin(5) * t372 - pkin(10) * t375 - t359 * t207 - t355 * t208 - t283 * t260 + t282 * t261;
t364 = -mrSges(5,1) * t230 + mrSges(5,2) * t231 - Ifges(5,5) * t292 - Ifges(5,6) * t291 - pkin(4) * t195 - t312 * t275 + t311 * t276 + t367;
t171 = -pkin(3) * t187 + t364 + (-Ifges(5,3) - Ifges(4,2)) * t306 + Ifges(4,6) * t344 + t345 * t295 - t327 * t293 + mrSges(4,3) * t246 - mrSges(4,1) * t281 + Ifges(4,4) * t307;
t315 = Ifges(3,3) * t345 + (Ifges(3,5) * t357 + Ifges(3,6) * t361) * t384;
t317 = Ifges(3,5) * t345 + (Ifges(3,1) * t357 + Ifges(3,4) * t361) * t384;
t158 = -mrSges(3,1) * t318 + mrSges(3,3) * t298 + Ifges(3,4) * t336 + Ifges(3,2) * t337 + Ifges(3,6) * t344 - pkin(2) * t370 + qJ(3) * t378 + t351 * t163 + t171 * t390 - t315 * t380 + t345 * t317;
t316 = Ifges(3,6) * t345 + (Ifges(3,4) * t357 + Ifges(3,2) * t361) * t384;
t160 = mrSges(3,2) * t318 - mrSges(3,3) * t297 + Ifges(3,1) * t336 + Ifges(3,4) * t337 + Ifges(3,5) * t344 - qJ(3) * t181 + t163 * t390 - t351 * t171 + t315 * t379 - t345 * t316;
t368 = mrSges(4,1) * t245 - mrSges(4,2) * t246 + Ifges(4,5) * t307 - Ifges(4,6) * t306 + Ifges(4,3) * t344 + pkin(3) * t365 + qJ(4) * t377 + t353 * t173 + t350 * t175 + t327 * t294 + t326 * t295;
t162 = t368 + pkin(2) * t181 + Ifges(3,3) * t344 + Ifges(3,6) * t337 + Ifges(3,5) * t336 + mrSges(3,1) * t297 - mrSges(3,2) * t298 + (t316 * t357 - t317 * t361) * t384;
t371 = mrSges(2,1) * t341 - mrSges(2,2) * t342 + Ifges(2,3) * qJDD(1) + pkin(1) * t167 + t158 * t387 + t160 * t388 + t354 * t162 + t170 * t391;
t168 = m(2) * t342 - t363 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t170;
t166 = t354 * t185 + (t179 * t361 + t180 * t357) * t352;
t164 = m(2) * t341 + qJDD(1) * mrSges(2,1) - t363 * mrSges(2,2) + t167;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t341 + Ifges(2,5) * qJDD(1) - t363 * Ifges(2,6) - t357 * t158 + t361 * t160 + (-t166 * t352 - t167 * t354) * pkin(8);
t155 = mrSges(2,1) * g(3) + mrSges(2,3) * t342 + t363 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t166 - t352 * t162 + (pkin(8) * t170 + t158 * t361 + t160 * t357) * t354;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t362 * t156 - t358 * t155 - pkin(7) * (t362 * t164 + t358 * t168), t156, t160, t163, t175, t188, t208; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t358 * t156 + t362 * t155 + pkin(7) * (-t358 * t164 + t362 * t168), t155, t158, t171, t173, t189, t207; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t371, t371, t162, t368, Ifges(5,3) * t306 - t364, -t367, t366;];
m_new  = t1;
