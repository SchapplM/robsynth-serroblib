% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-05-07 05:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:52:00
% EndTime: 2019-05-07 05:52:33
% DurationCPUTime: 15.16s
% Computational Cost: add. (271567->380), mult. (553655->463), div. (0->0), fcn. (377551->10), ass. (0->145)
t358 = sin(qJ(1));
t361 = cos(qJ(1));
t339 = t358 * g(1) - t361 * g(2);
t363 = qJD(1) ^ 2;
t316 = -qJDD(1) * pkin(1) - t363 * pkin(7) - t339;
t357 = sin(qJ(2));
t360 = cos(qJ(2));
t381 = qJD(1) * qJD(2);
t380 = t360 * t381;
t333 = t357 * qJDD(1) + t380;
t345 = t357 * t381;
t334 = t360 * qJDD(1) - t345;
t259 = (-t333 - t380) * pkin(8) + (-t334 + t345) * pkin(2) + t316;
t340 = -t361 * g(1) - t358 * g(2);
t317 = -t363 * pkin(1) + qJDD(1) * pkin(7) + t340;
t306 = -t357 * g(3) + t360 * t317;
t332 = (-pkin(2) * t360 - pkin(8) * t357) * qJD(1);
t362 = qJD(2) ^ 2;
t382 = t360 * qJD(1);
t264 = -t362 * pkin(2) + qJDD(2) * pkin(8) + t332 * t382 + t306;
t356 = sin(qJ(3));
t390 = cos(qJ(3));
t240 = t390 * t259 - t356 * t264;
t241 = t356 * t259 + t390 * t264;
t383 = qJD(1) * t357;
t329 = -t390 * qJD(2) + t356 * t383;
t330 = t356 * qJD(2) + t390 * t383;
t343 = -qJD(3) + t382;
t271 = Ifges(5,5) * t330 - Ifges(5,6) * t343 + Ifges(5,3) * t329;
t274 = Ifges(4,4) * t330 - Ifges(4,2) * t329 - Ifges(4,6) * t343;
t276 = Ifges(4,1) * t330 - Ifges(4,4) * t329 - Ifges(4,5) * t343;
t290 = t330 * qJD(3) - t390 * qJDD(2) + t356 * t333;
t291 = -t329 * qJD(3) + t356 * qJDD(2) + t390 * t333;
t298 = t329 * mrSges(5,1) - t330 * mrSges(5,3);
t328 = -qJDD(3) + t334;
t297 = t329 * pkin(3) - t330 * qJ(4);
t342 = t343 ^ 2;
t231 = t328 * pkin(3) - t342 * qJ(4) + t330 * t297 + qJDD(4) - t240;
t387 = t329 * t343;
t213 = (-t291 + t387) * qJ(5) + (t329 * t330 + t328) * pkin(4) + t231;
t391 = -2 * qJD(4);
t228 = -t342 * pkin(3) - t328 * qJ(4) - t329 * t297 + t343 * t391 + t241;
t301 = t343 * pkin(4) - t330 * qJ(5);
t327 = t329 ^ 2;
t217 = -t327 * pkin(4) + t290 * qJ(5) - t343 * t301 + t228;
t352 = sin(pkin(10));
t353 = cos(pkin(10));
t294 = t352 * t329 + t353 * t330;
t207 = -0.2e1 * qJD(5) * t294 + t353 * t213 - t352 * t217;
t252 = t352 * t290 + t353 * t291;
t293 = t353 * t329 - t352 * t330;
t203 = (t293 * t343 - t252) * pkin(9) + (t293 * t294 + t328) * pkin(5) + t207;
t208 = 0.2e1 * qJD(5) * t293 + t352 * t213 + t353 * t217;
t251 = t353 * t290 - t352 * t291;
t268 = t343 * pkin(5) - t294 * pkin(9);
t292 = t293 ^ 2;
t204 = -t292 * pkin(5) + t251 * pkin(9) - t343 * t268 + t208;
t355 = sin(qJ(6));
t359 = cos(qJ(6));
t201 = t359 * t203 - t355 * t204;
t255 = t359 * t293 - t355 * t294;
t223 = t255 * qJD(6) + t355 * t251 + t359 * t252;
t256 = t355 * t293 + t359 * t294;
t237 = -t255 * mrSges(7,1) + t256 * mrSges(7,2);
t341 = qJD(6) + t343;
t242 = -t341 * mrSges(7,2) + t255 * mrSges(7,3);
t324 = qJDD(6) + t328;
t196 = m(7) * t201 + t324 * mrSges(7,1) - t223 * mrSges(7,3) - t256 * t237 + t341 * t242;
t202 = t355 * t203 + t359 * t204;
t222 = -t256 * qJD(6) + t359 * t251 - t355 * t252;
t243 = t341 * mrSges(7,1) - t256 * mrSges(7,3);
t197 = m(7) * t202 - t324 * mrSges(7,2) + t222 * mrSges(7,3) + t255 * t237 - t341 * t243;
t187 = t359 * t196 + t355 * t197;
t257 = -t293 * mrSges(6,1) + t294 * mrSges(6,2);
t266 = -t343 * mrSges(6,2) + t293 * mrSges(6,3);
t184 = m(6) * t207 + t328 * mrSges(6,1) - t252 * mrSges(6,3) - t294 * t257 + t343 * t266 + t187;
t267 = t343 * mrSges(6,1) - t294 * mrSges(6,3);
t378 = -t355 * t196 + t359 * t197;
t185 = m(6) * t208 - t328 * mrSges(6,2) + t251 * mrSges(6,3) + t293 * t257 - t343 * t267 + t378;
t181 = t353 * t184 + t352 * t185;
t275 = Ifges(5,1) * t330 - Ifges(5,4) * t343 + Ifges(5,5) * t329;
t249 = Ifges(6,4) * t294 + Ifges(6,2) * t293 + Ifges(6,6) * t343;
t250 = Ifges(6,1) * t294 + Ifges(6,4) * t293 + Ifges(6,5) * t343;
t233 = Ifges(7,4) * t256 + Ifges(7,2) * t255 + Ifges(7,6) * t341;
t234 = Ifges(7,1) * t256 + Ifges(7,4) * t255 + Ifges(7,5) * t341;
t374 = mrSges(7,1) * t201 - mrSges(7,2) * t202 + Ifges(7,5) * t223 + Ifges(7,6) * t222 + Ifges(7,3) * t324 + t256 * t233 - t255 * t234;
t368 = mrSges(6,1) * t207 - mrSges(6,2) * t208 + Ifges(6,5) * t252 + Ifges(6,6) * t251 + Ifges(6,3) * t328 + pkin(5) * t187 + t294 * t249 - t293 * t250 + t374;
t365 = mrSges(5,1) * t231 - mrSges(5,3) * t228 - Ifges(5,4) * t291 + Ifges(5,2) * t328 - Ifges(5,6) * t290 + pkin(4) * t181 - t329 * t275 + t368;
t304 = -t329 * mrSges(5,2) - t343 * mrSges(5,3);
t371 = -m(5) * t231 - t328 * mrSges(5,1) - t343 * t304 - t181;
t182 = -t352 * t184 + t353 * t185;
t303 = t343 * mrSges(5,1) + t330 * mrSges(5,2);
t375 = m(5) * t228 - t328 * mrSges(5,3) - t343 * t303 + t182;
t393 = -(t271 - t274) * t330 + mrSges(4,1) * t240 - mrSges(4,2) * t241 + Ifges(4,5) * t291 - Ifges(4,6) * t290 - Ifges(4,3) * t328 + pkin(3) * (-t291 * mrSges(5,2) - t330 * t298 + t371) + qJ(4) * (-t290 * mrSges(5,2) - t329 * t298 + t375) + t329 * t276 - t365;
t305 = -t360 * g(3) - t357 * t317;
t263 = -qJDD(2) * pkin(2) - t362 * pkin(8) + t332 * t383 - t305;
t373 = t290 * pkin(3) + t263 + (-t291 - t387) * qJ(4);
t389 = pkin(3) * t343;
t216 = -t290 * pkin(4) - t327 * qJ(5) + qJDD(5) - t373 + ((2 * qJD(4)) + t301 + t389) * t330;
t210 = -t251 * pkin(5) - t292 * pkin(9) + t294 * t268 + t216;
t377 = m(7) * t210 - t222 * mrSges(7,1) + t223 * mrSges(7,2) - t255 * t242 + t256 * t243;
t198 = -m(6) * t216 + t251 * mrSges(6,1) - t252 * mrSges(6,2) + t293 * t266 - t294 * t267 - t377;
t230 = (t391 - t389) * t330 + t373;
t192 = m(5) * t230 + t290 * mrSges(5,1) - t291 * mrSges(5,3) - t330 * t303 + t329 * t304 + t198;
t232 = Ifges(7,5) * t256 + Ifges(7,6) * t255 + Ifges(7,3) * t341;
t188 = -mrSges(7,1) * t210 + mrSges(7,3) * t202 + Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t324 - t256 * t232 + t341 * t234;
t189 = mrSges(7,2) * t210 - mrSges(7,3) * t201 + Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t324 + t255 * t232 - t341 * t233;
t248 = Ifges(6,5) * t294 + Ifges(6,6) * t293 + Ifges(6,3) * t343;
t173 = -mrSges(6,1) * t216 + mrSges(6,3) * t208 + Ifges(6,4) * t252 + Ifges(6,2) * t251 + Ifges(6,6) * t328 - pkin(5) * t377 + pkin(9) * t378 + t359 * t188 + t355 * t189 - t294 * t248 + t343 * t250;
t175 = mrSges(6,2) * t216 - mrSges(6,3) * t207 + Ifges(6,1) * t252 + Ifges(6,4) * t251 + Ifges(6,5) * t328 - pkin(9) * t187 - t355 * t188 + t359 * t189 + t293 * t248 - t343 * t249;
t369 = -mrSges(5,1) * t230 + mrSges(5,2) * t228 - pkin(4) * t198 - qJ(5) * t182 - t353 * t173 - t352 * t175;
t273 = Ifges(5,4) * t330 - Ifges(5,2) * t343 + Ifges(5,6) * t329;
t385 = -Ifges(4,5) * t330 + Ifges(4,6) * t329 + Ifges(4,3) * t343 - t273;
t162 = -mrSges(4,1) * t263 + mrSges(4,3) * t241 - pkin(3) * t192 + (-t276 - t275) * t343 + t385 * t330 + (-Ifges(4,6) + Ifges(5,6)) * t328 + (Ifges(4,4) - Ifges(5,5)) * t291 + (-Ifges(4,2) - Ifges(5,3)) * t290 + t369;
t370 = mrSges(5,2) * t231 - mrSges(5,3) * t230 + Ifges(5,1) * t291 - Ifges(5,4) * t328 + Ifges(5,5) * t290 - qJ(5) * t181 - t352 * t173 + t353 * t175 - t343 * t271;
t163 = mrSges(4,2) * t263 - mrSges(4,3) * t240 + Ifges(4,1) * t291 - Ifges(4,4) * t290 - Ifges(4,5) * t328 - qJ(4) * t192 + t343 * t274 + t385 * t329 + t370;
t302 = -t343 * mrSges(4,1) - t330 * mrSges(4,3);
t384 = -t329 * mrSges(4,1) - t330 * mrSges(4,2) - t298;
t388 = -mrSges(4,3) - mrSges(5,2);
t177 = m(4) * t241 + t328 * mrSges(4,2) + t388 * t290 + t343 * t302 + t384 * t329 + t375;
t300 = t343 * mrSges(4,2) - t329 * mrSges(4,3);
t178 = m(4) * t240 - t328 * mrSges(4,1) + t388 * t291 - t343 * t300 + t384 * t330 + t371;
t172 = t390 * t177 - t356 * t178;
t191 = -m(4) * t263 - t290 * mrSges(4,1) - t291 * mrSges(4,2) - t329 * t300 - t330 * t302 - t192;
t314 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t357 + Ifges(3,2) * t360) * qJD(1);
t315 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t357 + Ifges(3,4) * t360) * qJD(1);
t392 = mrSges(3,1) * t305 - mrSges(3,2) * t306 + Ifges(3,5) * t333 + Ifges(3,6) * t334 + Ifges(3,3) * qJDD(2) + pkin(2) * t191 + pkin(8) * t172 + (t314 * t357 - t315 * t360) * qJD(1) + t390 * t162 + t356 * t163;
t331 = (-mrSges(3,1) * t360 + mrSges(3,2) * t357) * qJD(1);
t336 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t383;
t170 = m(3) * t306 - qJDD(2) * mrSges(3,2) + t334 * mrSges(3,3) - qJD(2) * t336 + t331 * t382 + t172;
t337 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t382;
t190 = m(3) * t305 + qJDD(2) * mrSges(3,1) - t333 * mrSges(3,3) + qJD(2) * t337 - t331 * t383 + t191;
t379 = t360 * t170 - t357 * t190;
t171 = t356 * t177 + t390 * t178;
t313 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t357 + Ifges(3,6) * t360) * qJD(1);
t159 = mrSges(3,2) * t316 - mrSges(3,3) * t305 + Ifges(3,1) * t333 + Ifges(3,4) * t334 + Ifges(3,5) * qJDD(2) - pkin(8) * t171 - qJD(2) * t314 - t356 * t162 + t390 * t163 + t313 * t382;
t161 = -mrSges(3,1) * t316 + mrSges(3,3) * t306 + Ifges(3,4) * t333 + Ifges(3,2) * t334 + Ifges(3,6) * qJDD(2) - pkin(2) * t171 + qJD(2) * t315 - t313 * t383 - t393;
t367 = -m(3) * t316 + t334 * mrSges(3,1) - t333 * mrSges(3,2) - t336 * t383 + t337 * t382 - t171;
t372 = mrSges(2,1) * t339 - mrSges(2,2) * t340 + Ifges(2,3) * qJDD(1) + pkin(1) * t367 + pkin(7) * t379 + t357 * t159 + t360 * t161;
t167 = m(2) * t339 + qJDD(1) * mrSges(2,1) - t363 * mrSges(2,2) + t367;
t166 = t357 * t170 + t360 * t190;
t164 = m(2) * t340 - t363 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t379;
t157 = mrSges(2,1) * g(3) + mrSges(2,3) * t340 + t363 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t166 - t392;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t339 + Ifges(2,5) * qJDD(1) - t363 * Ifges(2,6) - pkin(7) * t166 + t360 * t159 - t357 * t161;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t361 * t156 - t358 * t157 - pkin(6) * (t358 * t164 + t361 * t167), t156, t159, t163, -t329 * t273 + t370, t175, t189; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t358 * t156 + t361 * t157 + pkin(6) * (t361 * t164 - t358 * t167), t157, t161, t162, -t330 * t271 - t365, t173, t188; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t372, t372, t392, t393, Ifges(5,5) * t291 - Ifges(5,6) * t328 + Ifges(5,3) * t290 + t330 * t273 + t343 * t275 - t369, t368, t374;];
m_new  = t1;
