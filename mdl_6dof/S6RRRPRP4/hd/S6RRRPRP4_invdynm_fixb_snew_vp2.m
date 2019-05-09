% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP4
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
% Datum: 2019-05-07 07:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:44:24
% EndTime: 2019-05-07 07:44:49
% DurationCPUTime: 8.02s
% Computational Cost: add. (123457->384), mult. (254093->447), div. (0->0), fcn. (170644->8), ass. (0->136)
t336 = sin(qJ(2));
t338 = cos(qJ(2));
t365 = qJD(1) * qJD(2);
t314 = qJDD(1) * t336 + t338 * t365;
t337 = sin(qJ(1));
t339 = cos(qJ(1));
t321 = -g(1) * t339 - g(2) * t337;
t340 = qJD(1) ^ 2;
t308 = -pkin(1) * t340 + qJDD(1) * pkin(7) + t321;
t373 = t336 * t308;
t377 = pkin(2) * t340;
t245 = qJDD(2) * pkin(2) - t314 * pkin(8) - t373 + (pkin(8) * t365 + t336 * t377 - g(3)) * t338;
t289 = -g(3) * t336 + t308 * t338;
t315 = qJDD(1) * t338 - t336 * t365;
t368 = qJD(1) * t336;
t319 = qJD(2) * pkin(2) - pkin(8) * t368;
t333 = t338 ^ 2;
t246 = pkin(8) * t315 - qJD(2) * t319 - t333 * t377 + t289;
t335 = sin(qJ(3));
t379 = cos(qJ(3));
t212 = t245 * t379 - t246 * t335;
t213 = t245 * t335 + t246 * t379;
t306 = (t335 * t338 + t336 * t379) * qJD(1);
t265 = qJD(3) * t306 + t314 * t335 - t315 * t379;
t367 = qJD(1) * t338;
t305 = t335 * t368 - t367 * t379;
t266 = -qJD(3) * t305 + t314 * t379 + t315 * t335;
t331 = qJD(2) + qJD(3);
t275 = Ifges(4,4) * t306 - Ifges(4,2) * t305 + Ifges(4,6) * t331;
t283 = -mrSges(5,2) * t305 - mrSges(5,3) * t306;
t292 = mrSges(5,1) * t305 - mrSges(5,3) * t331;
t330 = qJDD(2) + qJDD(3);
t281 = pkin(3) * t305 - qJ(4) * t306;
t329 = t331 ^ 2;
t382 = -2 * qJD(4);
t208 = t329 * pkin(3) - qJ(4) * t330 + t281 * t305 + t331 * t382 - t213;
t293 = mrSges(5,1) * t306 + mrSges(5,2) * t331;
t294 = pkin(4) * t306 - pkin(9) * t331;
t301 = t305 ^ 2;
t205 = -t265 * pkin(4) - t301 * pkin(9) + t294 * t331 - t208;
t334 = sin(qJ(5));
t378 = cos(qJ(5));
t287 = t305 * t334 + t331 * t378;
t225 = qJD(5) * t287 - t265 * t378 + t330 * t334;
t286 = -t305 * t378 + t331 * t334;
t226 = -qJD(5) * t286 + t265 * t334 + t330 * t378;
t300 = qJD(5) + t306;
t199 = -0.2e1 * qJD(6) * t287 + (t286 * t300 - t226) * qJ(6) + (t287 * t300 + t225) * pkin(5) + t205;
t268 = -mrSges(7,2) * t286 + mrSges(7,3) * t300;
t271 = -mrSges(7,1) * t300 + mrSges(7,2) * t287;
t189 = m(7) * t199 + mrSges(7,1) * t225 - t226 * mrSges(7,3) + t268 * t286 - t287 * t271;
t269 = -mrSges(6,2) * t300 - mrSges(6,3) * t286;
t270 = mrSges(6,1) * t300 - mrSges(6,3) * t287;
t349 = m(6) * t205 + t225 * mrSges(6,1) + mrSges(6,2) * t226 + t286 * t269 + t270 * t287 + t189;
t347 = -m(5) * t208 + mrSges(5,3) * t330 + t293 * t331 + t349;
t320 = t337 * g(1) - g(2) * t339;
t357 = -qJDD(1) * pkin(1) - t320;
t267 = -t315 * pkin(2) + t319 * t368 + (-pkin(8) * t333 - pkin(7)) * t340 + t357;
t374 = t305 * t331;
t346 = (-t266 + t374) * qJ(4) + t267 + (pkin(3) * t331 + t382) * t306;
t201 = -t301 * pkin(4) - t306 * t294 + (pkin(3) + pkin(9)) * t265 + t346;
t210 = -t330 * pkin(3) - t329 * qJ(4) + t281 * t306 + qJDD(4) - t212;
t203 = (t305 * t306 - t330) * pkin(9) + (t266 + t374) * pkin(4) + t210;
t197 = t201 * t378 + t203 * t334;
t231 = Ifges(7,1) * t287 + Ifges(7,4) * t300 + Ifges(7,5) * t286;
t232 = Ifges(6,1) * t287 - Ifges(6,4) * t286 + Ifges(6,5) * t300;
t264 = qJDD(5) + t266;
t240 = pkin(5) * t286 - qJ(6) * t287;
t299 = t300 ^ 2;
t192 = -pkin(5) * t299 + qJ(6) * t264 + 0.2e1 * qJD(6) * t300 - t240 * t286 + t197;
t359 = -mrSges(7,1) * t199 + mrSges(7,2) * t192;
t229 = Ifges(7,4) * t287 + Ifges(7,2) * t300 + Ifges(7,6) * t286;
t372 = -Ifges(6,5) * t287 + Ifges(6,6) * t286 - Ifges(6,3) * t300 - t229;
t170 = -mrSges(6,1) * t205 + mrSges(6,3) * t197 - pkin(5) * t189 + (t231 + t232) * t300 + t372 * t287 + (Ifges(6,6) - Ifges(7,6)) * t264 + (Ifges(6,4) - Ifges(7,5)) * t226 + (-Ifges(6,2) - Ifges(7,3)) * t225 + t359;
t196 = -t201 * t334 + t203 * t378;
t230 = Ifges(6,4) * t287 - Ifges(6,2) * t286 + Ifges(6,6) * t300;
t194 = -pkin(5) * t264 - qJ(6) * t299 + t240 * t287 + qJDD(6) - t196;
t227 = Ifges(7,5) * t287 + Ifges(7,6) * t300 + Ifges(7,3) * t286;
t355 = mrSges(7,2) * t194 - mrSges(7,3) * t199 + Ifges(7,1) * t226 + Ifges(7,4) * t264 + Ifges(7,5) * t225 + t227 * t300;
t172 = mrSges(6,2) * t205 - mrSges(6,3) * t196 + Ifges(6,1) * t226 - Ifges(6,4) * t225 + Ifges(6,5) * t264 - qJ(6) * t189 - t300 * t230 + t286 * t372 + t355;
t364 = m(7) * t192 + mrSges(7,3) * t264 + t271 * t300;
t241 = mrSges(7,1) * t286 - mrSges(7,3) * t287;
t371 = -mrSges(6,1) * t286 - mrSges(6,2) * t287 - t241;
t376 = -mrSges(6,3) - mrSges(7,2);
t179 = m(6) * t197 - t264 * mrSges(6,2) + t225 * t376 - t300 * t270 + t286 * t371 + t364;
t360 = -m(7) * t194 + mrSges(7,1) * t264 + t268 * t300;
t181 = m(6) * t196 + t264 * mrSges(6,1) + t226 * t376 + t300 * t269 + t287 * t371 + t360;
t173 = t334 * t179 + t181 * t378;
t272 = Ifges(5,5) * t331 - Ifges(5,6) * t306 + Ifges(5,3) * t305;
t351 = -mrSges(5,2) * t210 + mrSges(5,3) * t208 - Ifges(5,1) * t330 + Ifges(5,4) * t266 - Ifges(5,5) * t265 + pkin(9) * t173 + t334 * t170 - t172 * t378 + t272 * t306;
t352 = -m(5) * t210 - mrSges(5,1) * t266 - t283 * t306 - t173;
t274 = Ifges(5,4) * t331 - Ifges(5,2) * t306 + Ifges(5,6) * t305;
t369 = Ifges(4,1) * t306 - Ifges(4,4) * t305 + Ifges(4,5) * t331 - t274;
t383 = t369 * t305 - mrSges(4,2) * t213 + pkin(3) * (-t330 * mrSges(5,2) - t331 * t292 + t352) + qJ(4) * (-t265 * mrSges(5,1) - t305 * t283 + t347) + mrSges(4,1) * t212 + t306 * t275 - Ifges(4,6) * t265 + Ifges(4,5) * t266 + Ifges(4,3) * t330 - t351;
t282 = mrSges(4,1) * t305 + mrSges(4,2) * t306;
t290 = -mrSges(4,2) * t331 - mrSges(4,3) * t305;
t166 = m(4) * t212 - t266 * mrSges(4,3) - t306 * t282 + (t290 - t292) * t331 + (mrSges(4,1) - mrSges(5,2)) * t330 + t352;
t291 = mrSges(4,1) * t331 - mrSges(4,3) * t306;
t177 = -t330 * mrSges(4,2) - t331 * t291 + m(4) * t213 + (-t282 - t283) * t305 + (-mrSges(4,3) - mrSges(5,1)) * t265 + t347;
t162 = t166 * t379 + t177 * t335;
t288 = -t338 * g(3) - t373;
t303 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t336 + Ifges(3,2) * t338) * qJD(1);
t304 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t336 + Ifges(3,4) * t338) * qJD(1);
t381 = mrSges(3,1) * t288 - mrSges(3,2) * t289 + Ifges(3,5) * t314 + Ifges(3,6) * t315 + Ifges(3,3) * qJDD(2) + pkin(2) * t162 + (t303 * t336 - t304 * t338) * qJD(1) + t383;
t375 = -Ifges(5,6) - Ifges(4,4);
t174 = t179 * t378 - t181 * t334;
t276 = Ifges(5,1) * t331 - Ifges(5,4) * t306 + Ifges(5,5) * t305;
t370 = -Ifges(4,5) * t306 + Ifges(4,6) * t305 - Ifges(4,3) * t331 - t276;
t313 = (-mrSges(3,1) * t338 + mrSges(3,2) * t336) * qJD(1);
t318 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t367;
t160 = m(3) * t288 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t314 + qJD(2) * t318 - t313 * t368 + t162;
t317 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t368;
t361 = -t166 * t335 + t177 * t379;
t161 = m(3) * t289 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t315 - qJD(2) * t317 + t313 * t367 + t361;
t362 = -t160 * t336 + t161 * t338;
t207 = t265 * pkin(3) + t346;
t167 = m(5) * t207 - mrSges(5,2) * t265 - mrSges(5,3) * t266 - t292 * t305 - t293 * t306 + t174;
t350 = mrSges(5,1) * t208 - mrSges(5,2) * t207 - pkin(4) * t349 + pkin(9) * t174 + t170 * t378 + t334 * t172;
t154 = -mrSges(4,1) * t267 + mrSges(4,3) * t213 - pkin(3) * t167 + t369 * t331 + (Ifges(4,6) - Ifges(5,5)) * t330 + t370 * t306 - t375 * t266 + (-Ifges(4,2) - Ifges(5,3)) * t265 - t350;
t353 = mrSges(7,1) * t194 - mrSges(7,3) * t192 - Ifges(7,4) * t226 - Ifges(7,2) * t264 - Ifges(7,6) * t225 + t287 * t227 - t231 * t286;
t345 = mrSges(6,2) * t197 - t286 * t232 - qJ(6) * (-t225 * mrSges(7,2) - t286 * t241 + t364) - pkin(5) * (-t226 * mrSges(7,2) - t287 * t241 + t360) - mrSges(6,1) * t196 - t287 * t230 + Ifges(6,6) * t225 - Ifges(6,5) * t226 - Ifges(6,3) * t264 + t353;
t342 = -mrSges(5,1) * t210 + mrSges(5,3) * t207 - pkin(4) * t173 + t345;
t155 = mrSges(4,2) * t267 - mrSges(4,3) * t212 - qJ(4) * t167 + (t272 - t275) * t331 + (Ifges(4,5) - Ifges(5,4)) * t330 + t370 * t305 + (Ifges(4,1) + Ifges(5,2)) * t266 + t375 * t265 - t342;
t302 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t336 + Ifges(3,6) * t338) * qJD(1);
t307 = -t340 * pkin(7) + t357;
t348 = m(4) * t267 + mrSges(4,1) * t265 + mrSges(4,2) * t266 + t290 * t305 + t291 * t306 + t167;
t150 = -mrSges(3,1) * t307 + mrSges(3,3) * t289 + Ifges(3,4) * t314 + Ifges(3,2) * t315 + Ifges(3,6) * qJDD(2) - pkin(2) * t348 + pkin(8) * t361 + qJD(2) * t304 + t154 * t379 + t335 * t155 - t302 * t368;
t152 = mrSges(3,2) * t307 - mrSges(3,3) * t288 + Ifges(3,1) * t314 + Ifges(3,4) * t315 + Ifges(3,5) * qJDD(2) - pkin(8) * t162 - qJD(2) * t303 - t154 * t335 + t155 * t379 + t302 * t367;
t343 = -m(3) * t307 + mrSges(3,1) * t315 - mrSges(3,2) * t314 - t317 * t368 + t318 * t367 - t348;
t354 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t343 + pkin(7) * t362 + t150 * t338 + t152 * t336;
t163 = m(2) * t320 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t340 + t343;
t158 = t160 * t338 + t161 * t336;
t156 = m(2) * t321 - mrSges(2,1) * t340 - qJDD(1) * mrSges(2,2) + t362;
t153 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t340 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t158 - t381;
t148 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t340 - pkin(7) * t158 - t150 * t336 + t152 * t338;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t339 * t148 - t337 * t153 - pkin(6) * (t156 * t337 + t163 * t339), t148, t152, t155, -t305 * t274 - t351, t172, -t229 * t286 + t355; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t337 * t148 + t339 * t153 + pkin(6) * (t156 * t339 - t163 * t337), t153, t150, t154, Ifges(5,4) * t330 - Ifges(5,2) * t266 + Ifges(5,6) * t265 - t331 * t272 + t305 * t276 + t342, t170, -t353; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t354, t354, t381, t383, Ifges(5,5) * t330 - Ifges(5,6) * t266 + Ifges(5,3) * t265 + t331 * t274 + t306 * t276 + t350, -t345, Ifges(7,5) * t226 + Ifges(7,6) * t264 + Ifges(7,3) * t225 + t287 * t229 - t300 * t231 - t359;];
m_new  = t1;
