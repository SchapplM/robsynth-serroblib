% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 10:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:10:50
% EndTime: 2019-05-07 10:11:34
% DurationCPUTime: 15.06s
% Computational Cost: add. (255362->382), mult. (547941->465), div. (0->0), fcn. (389272->10), ass. (0->146)
t356 = sin(qJ(2));
t360 = cos(qJ(2));
t381 = qJD(1) * qJD(2);
t323 = qJDD(1) * t356 + t360 * t381;
t362 = qJD(1) ^ 2;
t357 = sin(qJ(1));
t361 = cos(qJ(1));
t330 = -g(1) * t361 - g(2) * t357;
t316 = -pkin(1) * t362 + qJDD(1) * pkin(7) + t330;
t386 = t356 * t316;
t258 = qJDD(2) * pkin(2) - t323 * pkin(8) - t386 + (pkin(2) * t356 * t362 + pkin(8) * t381 - g(3)) * t360;
t300 = -g(3) * t356 + t316 * t360;
t324 = qJDD(1) * t360 - t356 * t381;
t383 = qJD(1) * t356;
t328 = qJD(2) * pkin(2) - pkin(8) * t383;
t387 = t360 ^ 2 * t362;
t259 = -pkin(2) * t387 + pkin(8) * t324 - qJD(2) * t328 + t300;
t355 = sin(qJ(3));
t391 = cos(qJ(3));
t242 = t258 * t355 + t259 * t391;
t314 = (t355 * t360 + t356 * t391) * qJD(1);
t275 = qJD(3) * t314 + t323 * t355 - t324 * t391;
t348 = qJD(2) + qJD(3);
t302 = mrSges(4,1) * t348 - mrSges(4,3) * t314;
t382 = qJD(1) * t360;
t313 = t355 * t383 - t382 * t391;
t347 = qJDD(2) + qJDD(3);
t241 = t258 * t391 - t259 * t355;
t294 = pkin(3) * t313 - qJ(4) * t314;
t346 = t348 ^ 2;
t226 = -t347 * pkin(3) - t346 * qJ(4) + t294 * t314 + qJDD(4) - t241;
t276 = -qJD(3) * t313 + t323 * t391 + t324 * t355;
t388 = t313 * t348;
t214 = (-t276 - t388) * pkin(9) + (t313 * t314 - t347) * pkin(4) + t226;
t392 = 2 * qJD(4);
t224 = -pkin(3) * t346 + qJ(4) * t347 - t294 * t313 + t348 * t392 + t242;
t305 = -pkin(4) * t348 - pkin(9) * t314;
t309 = t313 ^ 2;
t216 = -pkin(4) * t309 + pkin(9) * t275 + t305 * t348 + t224;
t354 = sin(qJ(5));
t359 = cos(qJ(5));
t210 = t214 * t354 + t216 * t359;
t293 = t313 * t354 + t314 * t359;
t236 = -qJD(5) * t293 + t275 * t359 - t276 * t354;
t292 = t313 * t359 - t314 * t354;
t251 = -mrSges(6,1) * t292 + mrSges(6,2) * t293;
t343 = qJD(5) - t348;
t279 = mrSges(6,1) * t343 - mrSges(6,3) * t293;
t342 = qJDD(5) - t347;
t252 = -pkin(5) * t292 - pkin(10) * t293;
t341 = t343 ^ 2;
t205 = -pkin(5) * t341 + pkin(10) * t342 + t252 * t292 + t210;
t329 = g(1) * t357 - g(2) * t361;
t315 = -qJDD(1) * pkin(1) - pkin(7) * t362 - t329;
t277 = -pkin(2) * t324 - pkin(8) * t387 + t328 * t383 + t315;
t372 = t275 * pkin(3) + t277 + (-t276 + t388) * qJ(4);
t390 = pkin(3) * t348;
t212 = -t275 * pkin(4) - t309 * pkin(9) - t372 + (t305 - t390 + t392) * t314;
t237 = qJD(5) * t292 + t275 * t354 + t276 * t359;
t206 = t212 + (-t292 * t343 - t237) * pkin(10) + (t293 * t343 - t236) * pkin(5);
t353 = sin(qJ(6));
t358 = cos(qJ(6));
t202 = -t205 * t353 + t206 * t358;
t262 = -t293 * t353 + t343 * t358;
t219 = qJD(6) * t262 + t237 * t358 + t342 * t353;
t234 = qJDD(6) - t236;
t263 = t293 * t358 + t343 * t353;
t243 = -mrSges(7,1) * t262 + mrSges(7,2) * t263;
t287 = qJD(6) - t292;
t244 = -mrSges(7,2) * t287 + mrSges(7,3) * t262;
t198 = m(7) * t202 + mrSges(7,1) * t234 - mrSges(7,3) * t219 - t243 * t263 + t244 * t287;
t203 = t205 * t358 + t206 * t353;
t218 = -qJD(6) * t263 - t237 * t353 + t342 * t358;
t245 = mrSges(7,1) * t287 - mrSges(7,3) * t263;
t199 = m(7) * t203 - mrSges(7,2) * t234 + mrSges(7,3) * t218 + t243 * t262 - t245 * t287;
t378 = -t198 * t353 + t199 * t358;
t185 = m(6) * t210 - mrSges(6,2) * t342 + mrSges(6,3) * t236 + t251 * t292 - t279 * t343 + t378;
t209 = t214 * t359 - t216 * t354;
t278 = -mrSges(6,2) * t343 + mrSges(6,3) * t292;
t204 = -pkin(5) * t342 - pkin(10) * t341 + t252 * t293 - t209;
t374 = -m(7) * t204 + mrSges(7,1) * t218 - mrSges(7,2) * t219 + t244 * t262 - t245 * t263;
t194 = m(6) * t209 + mrSges(6,1) * t342 - mrSges(6,3) * t237 - t251 * t293 + t278 * t343 + t374;
t180 = t185 * t359 - t354 * t194;
t303 = -mrSges(5,1) * t348 + mrSges(5,2) * t314;
t376 = m(5) * t224 + mrSges(5,3) * t347 + t303 * t348 + t180;
t295 = mrSges(5,1) * t313 - mrSges(5,3) * t314;
t384 = -mrSges(4,1) * t313 - mrSges(4,2) * t314 - t295;
t389 = -mrSges(4,3) - mrSges(5,2);
t172 = m(4) * t242 - t347 * mrSges(4,2) + t275 * t389 - t348 * t302 + t313 * t384 + t376;
t301 = -mrSges(4,2) * t348 - mrSges(4,3) * t313;
t179 = t354 * t185 + t359 * t194;
t304 = -mrSges(5,2) * t313 + mrSges(5,3) * t348;
t373 = -m(5) * t226 + mrSges(5,1) * t347 + t304 * t348 - t179;
t174 = m(4) * t241 + t347 * mrSges(4,1) + t276 * t389 + t348 * t301 + t314 * t384 + t373;
t166 = t172 * t355 + t174 * t391;
t299 = -t360 * g(3) - t386;
t311 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t356 + Ifges(3,2) * t360) * qJD(1);
t312 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t356 + Ifges(3,4) * t360) * qJD(1);
t284 = Ifges(4,4) * t314 - Ifges(4,2) * t313 + Ifges(4,6) * t348;
t286 = Ifges(4,1) * t314 - Ifges(4,4) * t313 + Ifges(4,5) * t348;
t281 = Ifges(5,5) * t314 + Ifges(5,6) * t348 + Ifges(5,3) * t313;
t285 = Ifges(5,1) * t314 + Ifges(5,4) * t348 + Ifges(5,5) * t313;
t227 = Ifges(7,5) * t263 + Ifges(7,6) * t262 + Ifges(7,3) * t287;
t229 = Ifges(7,1) * t263 + Ifges(7,4) * t262 + Ifges(7,5) * t287;
t192 = -mrSges(7,1) * t204 + mrSges(7,3) * t203 + Ifges(7,4) * t219 + Ifges(7,2) * t218 + Ifges(7,6) * t234 - t227 * t263 + t229 * t287;
t228 = Ifges(7,4) * t263 + Ifges(7,2) * t262 + Ifges(7,6) * t287;
t193 = mrSges(7,2) * t204 - mrSges(7,3) * t202 + Ifges(7,1) * t219 + Ifges(7,4) * t218 + Ifges(7,5) * t234 + t227 * t262 - t228 * t287;
t247 = Ifges(6,4) * t293 + Ifges(6,2) * t292 + Ifges(6,6) * t343;
t248 = Ifges(6,1) * t293 + Ifges(6,4) * t292 + Ifges(6,5) * t343;
t371 = mrSges(6,1) * t209 - mrSges(6,2) * t210 + Ifges(6,5) * t237 + Ifges(6,6) * t236 + Ifges(6,3) * t342 + pkin(5) * t374 + pkin(10) * t378 + t192 * t358 + t193 * t353 + t247 * t293 - t248 * t292;
t366 = mrSges(5,1) * t226 - mrSges(5,3) * t224 - Ifges(5,4) * t276 - Ifges(5,2) * t347 - Ifges(5,6) * t275 + pkin(4) * t179 + t314 * t281 - t285 * t313 + t371;
t364 = -mrSges(4,2) * t242 + t313 * t286 + qJ(4) * (-t275 * mrSges(5,2) - t313 * t295 + t376) + pkin(3) * (-t276 * mrSges(5,2) - t314 * t295 + t373) + mrSges(4,1) * t241 + t314 * t284 - Ifges(4,6) * t275 + Ifges(4,5) * t276 + Ifges(4,3) * t347 - t366;
t393 = mrSges(3,1) * t299 - mrSges(3,2) * t300 + Ifges(3,5) * t323 + Ifges(3,6) * t324 + Ifges(3,3) * qJDD(2) + pkin(2) * t166 + (t311 * t356 - t312 * t360) * qJD(1) + t364;
t188 = t198 * t358 + t199 * t353;
t283 = Ifges(5,4) * t314 + Ifges(5,2) * t348 + Ifges(5,6) * t313;
t385 = -Ifges(4,5) * t314 + Ifges(4,6) * t313 - Ifges(4,3) * t348 - t283;
t322 = (-mrSges(3,1) * t360 + mrSges(3,2) * t356) * qJD(1);
t327 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t382;
t164 = m(3) * t299 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t323 + qJD(2) * t327 - t322 * t383 + t166;
t326 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t383;
t379 = t172 * t391 - t174 * t355;
t165 = m(3) * t300 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t324 - qJD(2) * t326 + t322 * t382 + t379;
t380 = -t164 * t356 + t165 * t360;
t186 = -m(6) * t212 + mrSges(6,1) * t236 - mrSges(6,2) * t237 + t278 * t292 - t279 * t293 - t188;
t221 = (-(2 * qJD(4)) + t390) * t314 + t372;
t183 = m(5) * t221 + mrSges(5,1) * t275 - mrSges(5,3) * t276 - t303 * t314 + t304 * t313 + t186;
t246 = Ifges(6,5) * t293 + Ifges(6,6) * t292 + Ifges(6,3) * t343;
t168 = mrSges(6,2) * t212 - mrSges(6,3) * t209 + Ifges(6,1) * t237 + Ifges(6,4) * t236 + Ifges(6,5) * t342 - pkin(10) * t188 - t192 * t353 + t193 * t358 + t246 * t292 - t247 * t343;
t368 = mrSges(7,1) * t202 - mrSges(7,2) * t203 + Ifges(7,5) * t219 + Ifges(7,6) * t218 + Ifges(7,3) * t234 + t228 * t263 - t229 * t262;
t169 = -mrSges(6,1) * t212 + mrSges(6,3) * t210 + Ifges(6,4) * t237 + Ifges(6,2) * t236 + Ifges(6,6) * t342 - pkin(5) * t188 - t246 * t293 + t248 * t343 - t368;
t369 = -mrSges(5,1) * t221 + mrSges(5,2) * t224 - pkin(4) * t186 - pkin(9) * t180 - t354 * t168 - t359 * t169;
t158 = -mrSges(4,1) * t277 + mrSges(4,3) * t242 - pkin(3) * t183 + (t286 + t285) * t348 + (Ifges(4,6) - Ifges(5,6)) * t347 + t385 * t314 + (Ifges(4,4) - Ifges(5,5)) * t276 + (-Ifges(4,2) - Ifges(5,3)) * t275 + t369;
t370 = mrSges(5,2) * t226 - mrSges(5,3) * t221 + Ifges(5,1) * t276 + Ifges(5,4) * t347 + Ifges(5,5) * t275 - pkin(9) * t179 + t168 * t359 - t354 * t169 + t281 * t348;
t159 = mrSges(4,2) * t277 - mrSges(4,3) * t241 + Ifges(4,1) * t276 - Ifges(4,4) * t275 + Ifges(4,5) * t347 - qJ(4) * t183 - t348 * t284 + t313 * t385 + t370;
t310 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t356 + Ifges(3,6) * t360) * qJD(1);
t367 = m(4) * t277 + mrSges(4,1) * t275 + mrSges(4,2) * t276 + t301 * t313 + t302 * t314 + t183;
t154 = -mrSges(3,1) * t315 + mrSges(3,3) * t300 + Ifges(3,4) * t323 + Ifges(3,2) * t324 + Ifges(3,6) * qJDD(2) - pkin(2) * t367 + pkin(8) * t379 + qJD(2) * t312 + t158 * t391 + t355 * t159 - t310 * t383;
t156 = mrSges(3,2) * t315 - mrSges(3,3) * t299 + Ifges(3,1) * t323 + Ifges(3,4) * t324 + Ifges(3,5) * qJDD(2) - pkin(8) * t166 - qJD(2) * t311 - t158 * t355 + t159 * t391 + t310 * t382;
t365 = -m(3) * t315 + mrSges(3,1) * t324 - mrSges(3,2) * t323 - t326 * t383 + t327 * t382 - t367;
t375 = mrSges(2,1) * t329 - mrSges(2,2) * t330 + Ifges(2,3) * qJDD(1) + pkin(1) * t365 + pkin(7) * t380 + t154 * t360 + t156 * t356;
t181 = m(2) * t329 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t362 + t365;
t162 = t164 * t360 + t165 * t356;
t160 = m(2) * t330 - mrSges(2,1) * t362 - qJDD(1) * mrSges(2,2) + t380;
t157 = mrSges(2,1) * g(3) + mrSges(2,3) * t330 + t362 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t393;
t152 = -mrSges(2,2) * g(3) - mrSges(2,3) * t329 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t362 - pkin(7) * t162 - t154 * t356 + t156 * t360;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t361 * t152 - t357 * t157 - pkin(6) * (t160 * t357 + t181 * t361), t152, t156, t159, -t283 * t313 + t370, t168, t193; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t357 * t152 + t361 * t157 + pkin(6) * (t160 * t361 - t181 * t357), t157, t154, t158, -t366, t169, t192; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t375, t375, t393, t364, Ifges(5,5) * t276 + Ifges(5,6) * t347 + Ifges(5,3) * t275 + t314 * t283 - t348 * t285 - t369, t371, t368;];
m_new  = t1;
