% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR10
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
% Datum: 2019-05-07 14:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 13:58:18
% EndTime: 2019-05-07 13:59:00
% DurationCPUTime: 15.58s
% Computational Cost: add. (285874->381), mult. (571352->461), div. (0->0), fcn. (394711->10), ass. (0->147)
t357 = sin(qJ(1));
t361 = cos(qJ(1));
t339 = t357 * g(1) - t361 * g(2);
t363 = qJD(1) ^ 2;
t316 = -qJDD(1) * pkin(1) - t363 * pkin(7) - t339;
t356 = sin(qJ(2));
t360 = cos(qJ(2));
t382 = qJD(1) * qJD(2);
t380 = t360 * t382;
t332 = t356 * qJDD(1) + t380;
t381 = t356 * t382;
t333 = t360 * qJDD(1) - t381;
t259 = (-t332 - t380) * pkin(8) + (-t333 + t381) * pkin(2) + t316;
t340 = -t361 * g(1) - t357 * g(2);
t317 = -t363 * pkin(1) + qJDD(1) * pkin(7) + t340;
t303 = -t356 * g(3) + t360 * t317;
t331 = (-pkin(2) * t360 - pkin(8) * t356) * qJD(1);
t347 = t360 * qJD(1);
t362 = qJD(2) ^ 2;
t264 = -t362 * pkin(2) + qJDD(2) * pkin(8) + t331 * t347 + t303;
t355 = sin(qJ(3));
t390 = cos(qJ(3));
t246 = t390 * t259 - t355 * t264;
t247 = t355 * t259 + t390 * t264;
t383 = qJD(1) * t356;
t328 = -t390 * qJD(2) + t355 * t383;
t329 = t355 * qJD(2) + t390 * t383;
t343 = -t347 + qJD(3);
t271 = Ifges(5,5) * t329 + Ifges(5,6) * t343 + Ifges(5,3) * t328;
t274 = Ifges(4,4) * t329 - Ifges(4,2) * t328 + Ifges(4,6) * t343;
t276 = Ifges(4,1) * t329 - Ifges(4,4) * t328 + Ifges(4,5) * t343;
t288 = t329 * qJD(3) - t390 * qJDD(2) + t355 * t332;
t289 = -t328 * qJD(3) + t355 * qJDD(2) + t390 * t332;
t296 = t328 * mrSges(5,1) - t329 * mrSges(5,3);
t327 = qJDD(3) - t333;
t295 = t328 * pkin(3) - t329 * qJ(4);
t342 = t343 ^ 2;
t231 = -t327 * pkin(3) - t342 * qJ(4) + t329 * t295 + qJDD(4) - t246;
t387 = t328 * t343;
t219 = (-t289 - t387) * pkin(9) + (t328 * t329 - t327) * pkin(4) + t231;
t391 = 2 * qJD(4);
t228 = -t342 * pkin(3) + t327 * qJ(4) - t328 * t295 + t343 * t391 + t247;
t304 = -t343 * pkin(4) - t329 * pkin(9);
t326 = t328 ^ 2;
t222 = -t326 * pkin(4) + t288 * pkin(9) + t343 * t304 + t228;
t354 = sin(qJ(5));
t359 = cos(qJ(5));
t207 = t359 * t219 - t354 * t222;
t291 = t359 * t328 - t354 * t329;
t245 = t291 * qJD(5) + t354 * t288 + t359 * t289;
t292 = t354 * t328 + t359 * t329;
t323 = qJDD(5) - t327;
t341 = qJD(5) - t343;
t203 = (t291 * t341 - t245) * pkin(10) + (t291 * t292 + t323) * pkin(5) + t207;
t208 = t354 * t219 + t359 * t222;
t244 = -t292 * qJD(5) + t359 * t288 - t354 * t289;
t268 = t341 * pkin(5) - t292 * pkin(10);
t290 = t291 ^ 2;
t204 = -t290 * pkin(5) + t244 * pkin(10) - t341 * t268 + t208;
t353 = sin(qJ(6));
t358 = cos(qJ(6));
t201 = t358 * t203 - t353 * t204;
t255 = t358 * t291 - t353 * t292;
t216 = t255 * qJD(6) + t353 * t244 + t358 * t245;
t256 = t353 * t291 + t358 * t292;
t237 = -t255 * mrSges(7,1) + t256 * mrSges(7,2);
t335 = qJD(6) + t341;
t248 = -t335 * mrSges(7,2) + t255 * mrSges(7,3);
t311 = qJDD(6) + t323;
t196 = m(7) * t201 + t311 * mrSges(7,1) - t216 * mrSges(7,3) - t256 * t237 + t335 * t248;
t202 = t353 * t203 + t358 * t204;
t215 = -t256 * qJD(6) + t358 * t244 - t353 * t245;
t249 = t335 * mrSges(7,1) - t256 * mrSges(7,3);
t197 = m(7) * t202 - t311 * mrSges(7,2) + t215 * mrSges(7,3) + t255 * t237 - t335 * t249;
t187 = t358 * t196 + t353 * t197;
t257 = -t291 * mrSges(6,1) + t292 * mrSges(6,2);
t266 = -t341 * mrSges(6,2) + t291 * mrSges(6,3);
t184 = m(6) * t207 + t323 * mrSges(6,1) - t245 * mrSges(6,3) - t292 * t257 + t341 * t266 + t187;
t267 = t341 * mrSges(6,1) - t292 * mrSges(6,3);
t378 = -t353 * t196 + t358 * t197;
t185 = m(6) * t208 - t323 * mrSges(6,2) + t244 * mrSges(6,3) + t291 * t257 - t341 * t267 + t378;
t181 = t359 * t184 + t354 * t185;
t275 = Ifges(5,1) * t329 + Ifges(5,4) * t343 + Ifges(5,5) * t328;
t251 = Ifges(6,4) * t292 + Ifges(6,2) * t291 + Ifges(6,6) * t341;
t252 = Ifges(6,1) * t292 + Ifges(6,4) * t291 + Ifges(6,5) * t341;
t233 = Ifges(7,4) * t256 + Ifges(7,2) * t255 + Ifges(7,6) * t335;
t234 = Ifges(7,1) * t256 + Ifges(7,4) * t255 + Ifges(7,5) * t335;
t374 = mrSges(7,1) * t201 - mrSges(7,2) * t202 + Ifges(7,5) * t216 + Ifges(7,6) * t215 + Ifges(7,3) * t311 + t256 * t233 - t255 * t234;
t368 = mrSges(6,1) * t207 - mrSges(6,2) * t208 + Ifges(6,5) * t245 + Ifges(6,6) * t244 + Ifges(6,3) * t323 + pkin(5) * t187 + t292 * t251 - t291 * t252 + t374;
t365 = mrSges(5,1) * t231 - mrSges(5,3) * t228 - Ifges(5,4) * t289 - Ifges(5,2) * t327 - Ifges(5,6) * t288 + pkin(4) * t181 - t328 * t275 + t368;
t301 = -t328 * mrSges(5,2) + t343 * mrSges(5,3);
t371 = -m(5) * t231 + t327 * mrSges(5,1) + t343 * t301 - t181;
t182 = -t354 * t184 + t359 * t185;
t300 = -t343 * mrSges(5,1) + t329 * mrSges(5,2);
t375 = m(5) * t228 + t327 * mrSges(5,3) + t343 * t300 + t182;
t393 = -(t271 - t274) * t329 + mrSges(4,1) * t246 - mrSges(4,2) * t247 + Ifges(4,5) * t289 - Ifges(4,6) * t288 + Ifges(4,3) * t327 + pkin(3) * (-t289 * mrSges(5,2) - t329 * t296 + t371) + qJ(4) * (-t288 * mrSges(5,2) - t328 * t296 + t375) + t328 * t276 - t365;
t302 = -t360 * g(3) - t356 * t317;
t263 = -qJDD(2) * pkin(2) - t362 * pkin(8) + t331 * t383 - t302;
t373 = t288 * pkin(3) + t263 + (-t289 + t387) * qJ(4);
t389 = pkin(3) * t343;
t223 = -t288 * pkin(4) - t326 * pkin(9) - t373 + (t304 - t389 + t391) * t329;
t210 = -t244 * pkin(5) - t290 * pkin(10) + t292 * t268 + t223;
t377 = m(7) * t210 - t215 * mrSges(7,1) + t216 * mrSges(7,2) - t255 * t248 + t256 * t249;
t198 = -m(6) * t223 + t244 * mrSges(6,1) - t245 * mrSges(6,2) + t291 * t266 - t292 * t267 - t377;
t230 = (-(2 * qJD(4)) + t389) * t329 + t373;
t192 = m(5) * t230 + t288 * mrSges(5,1) - t289 * mrSges(5,3) - t329 * t300 + t328 * t301 + t198;
t232 = Ifges(7,5) * t256 + Ifges(7,6) * t255 + Ifges(7,3) * t335;
t188 = -mrSges(7,1) * t210 + mrSges(7,3) * t202 + Ifges(7,4) * t216 + Ifges(7,2) * t215 + Ifges(7,6) * t311 - t256 * t232 + t335 * t234;
t189 = mrSges(7,2) * t210 - mrSges(7,3) * t201 + Ifges(7,1) * t216 + Ifges(7,4) * t215 + Ifges(7,5) * t311 + t255 * t232 - t335 * t233;
t250 = Ifges(6,5) * t292 + Ifges(6,6) * t291 + Ifges(6,3) * t341;
t173 = -mrSges(6,1) * t223 + mrSges(6,3) * t208 + Ifges(6,4) * t245 + Ifges(6,2) * t244 + Ifges(6,6) * t323 - pkin(5) * t377 + pkin(10) * t378 + t358 * t188 + t353 * t189 - t292 * t250 + t341 * t252;
t175 = mrSges(6,2) * t223 - mrSges(6,3) * t207 + Ifges(6,1) * t245 + Ifges(6,4) * t244 + Ifges(6,5) * t323 - pkin(10) * t187 - t353 * t188 + t358 * t189 + t291 * t250 - t341 * t251;
t369 = -mrSges(5,1) * t230 + mrSges(5,2) * t228 - pkin(4) * t198 - pkin(9) * t182 - t359 * t173 - t354 * t175;
t273 = Ifges(5,4) * t329 + Ifges(5,2) * t343 + Ifges(5,6) * t328;
t385 = -Ifges(4,5) * t329 + Ifges(4,6) * t328 - Ifges(4,3) * t343 - t273;
t162 = -mrSges(4,1) * t263 + mrSges(4,3) * t247 - pkin(3) * t192 + (t276 + t275) * t343 + t385 * t329 + (Ifges(4,6) - Ifges(5,6)) * t327 + (Ifges(4,4) - Ifges(5,5)) * t289 + (-Ifges(4,2) - Ifges(5,3)) * t288 + t369;
t370 = mrSges(5,2) * t231 - mrSges(5,3) * t230 + Ifges(5,1) * t289 + Ifges(5,4) * t327 + Ifges(5,5) * t288 - pkin(9) * t181 - t354 * t173 + t359 * t175 + t343 * t271;
t163 = mrSges(4,2) * t263 - mrSges(4,3) * t246 + Ifges(4,1) * t289 - Ifges(4,4) * t288 + Ifges(4,5) * t327 - qJ(4) * t192 - t343 * t274 + t385 * t328 + t370;
t299 = t343 * mrSges(4,1) - t329 * mrSges(4,3);
t384 = -t328 * mrSges(4,1) - t329 * mrSges(4,2) - t296;
t388 = -mrSges(4,3) - mrSges(5,2);
t177 = m(4) * t247 - t327 * mrSges(4,2) + t388 * t288 - t343 * t299 + t384 * t328 + t375;
t298 = -t343 * mrSges(4,2) - t328 * mrSges(4,3);
t178 = m(4) * t246 + t327 * mrSges(4,1) + t388 * t289 + t343 * t298 + t384 * t329 + t371;
t172 = t390 * t177 - t355 * t178;
t191 = -m(4) * t263 - t288 * mrSges(4,1) - t289 * mrSges(4,2) - t328 * t298 - t329 * t299 - t192;
t314 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t356 + Ifges(3,2) * t360) * qJD(1);
t315 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t356 + Ifges(3,4) * t360) * qJD(1);
t392 = mrSges(3,1) * t302 - mrSges(3,2) * t303 + Ifges(3,5) * t332 + Ifges(3,6) * t333 + Ifges(3,3) * qJDD(2) + pkin(2) * t191 + pkin(8) * t172 + (t314 * t356 - t315 * t360) * qJD(1) + t390 * t162 + t355 * t163;
t330 = (-mrSges(3,1) * t360 + mrSges(3,2) * t356) * qJD(1);
t336 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t383;
t170 = m(3) * t303 - qJDD(2) * mrSges(3,2) + t333 * mrSges(3,3) - qJD(2) * t336 + t330 * t347 + t172;
t337 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t347;
t190 = m(3) * t302 + qJDD(2) * mrSges(3,1) - t332 * mrSges(3,3) + qJD(2) * t337 - t330 * t383 + t191;
t379 = t360 * t170 - t356 * t190;
t171 = t355 * t177 + t390 * t178;
t313 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t356 + Ifges(3,6) * t360) * qJD(1);
t159 = mrSges(3,2) * t316 - mrSges(3,3) * t302 + Ifges(3,1) * t332 + Ifges(3,4) * t333 + Ifges(3,5) * qJDD(2) - pkin(8) * t171 - qJD(2) * t314 - t355 * t162 + t390 * t163 + t313 * t347;
t161 = -mrSges(3,1) * t316 + mrSges(3,3) * t303 + Ifges(3,4) * t332 + Ifges(3,2) * t333 + Ifges(3,6) * qJDD(2) - pkin(2) * t171 + qJD(2) * t315 - t313 * t383 - t393;
t367 = -m(3) * t316 + t333 * mrSges(3,1) - t332 * mrSges(3,2) - t336 * t383 + t337 * t347 - t171;
t372 = mrSges(2,1) * t339 - mrSges(2,2) * t340 + Ifges(2,3) * qJDD(1) + pkin(1) * t367 + pkin(7) * t379 + t356 * t159 + t360 * t161;
t167 = m(2) * t339 + qJDD(1) * mrSges(2,1) - t363 * mrSges(2,2) + t367;
t166 = t356 * t170 + t360 * t190;
t164 = m(2) * t340 - t363 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t379;
t157 = mrSges(2,1) * g(3) + mrSges(2,3) * t340 + t363 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t166 - t392;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t339 + Ifges(2,5) * qJDD(1) - t363 * Ifges(2,6) - pkin(7) * t166 + t360 * t159 - t356 * t161;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t361 * t156 - t357 * t157 - pkin(6) * (t357 * t164 + t361 * t167), t156, t159, t163, -t328 * t273 + t370, t175, t189; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t357 * t156 + t361 * t157 + pkin(6) * (t361 * t164 - t357 * t167), t157, t161, t162, -t329 * t271 - t365, t173, t188; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t372, t372, t392, t393, Ifges(5,5) * t289 + Ifges(5,6) * t327 + Ifges(5,3) * t288 + t329 * t273 - t343 * t275 - t369, t368, t374;];
m_new  = t1;
