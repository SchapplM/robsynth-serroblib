% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-05-05 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:57:39
% EndTime: 2019-05-05 02:57:49
% DurationCPUTime: 5.51s
% Computational Cost: add. (66496->350), mult. (132050->417), div. (0->0), fcn. (72122->10), ass. (0->143)
t340 = sin(pkin(10));
t342 = cos(pkin(10));
t311 = g(1) * t340 - g(2) * t342;
t312 = -g(1) * t342 - g(2) * t340;
t338 = -g(3) + qJDD(1);
t350 = cos(qJ(2));
t343 = cos(pkin(6));
t347 = sin(qJ(2));
t387 = t343 * t347;
t341 = sin(pkin(6));
t389 = t341 * t347;
t237 = t311 * t387 + t350 * t312 + t338 * t389;
t352 = qJD(2) ^ 2;
t235 = -pkin(2) * t352 + qJDD(2) * pkin(8) + t237;
t251 = -t311 * t341 + t338 * t343;
t346 = sin(qJ(3));
t349 = cos(qJ(3));
t229 = -t346 * t235 + t349 * t251;
t230 = t349 * t235 + t346 * t251;
t265 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t346 - Ifges(5,3) * t349) * qJD(2);
t269 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t346 + Ifges(4,2) * t349) * qJD(2);
t301 = (-mrSges(5,1) * t349 - mrSges(5,3) * t346) * qJD(2);
t378 = qJD(2) * qJD(3);
t374 = t349 * t378;
t305 = qJDD(2) * t346 + t374;
t373 = t346 * t378;
t306 = qJDD(2) * t349 - t373;
t379 = qJD(2) * t349;
t317 = -qJD(3) * mrSges(6,1) + mrSges(6,3) * t379;
t351 = qJD(3) ^ 2;
t300 = (-pkin(3) * t349 - qJ(4) * t346) * qJD(2);
t380 = qJD(2) * t346;
t370 = t300 * t380 + qJDD(4) - t229;
t227 = -qJDD(3) * pkin(3) - t351 * qJ(4) + t370;
t385 = t349 * t352;
t377 = qJD(2) * qJD(5);
t405 = -0.2e1 * t346 * t377 + (-t305 + t374) * qJ(5);
t218 = (-t346 * t385 - qJDD(3)) * pkin(4) + t227 + t405;
t303 = (mrSges(6,1) * t346 - mrSges(6,2) * t349) * qJD(2);
t313 = -qJD(3) * pkin(4) - qJ(5) * t380;
t386 = t343 * t350;
t388 = t341 * t350;
t236 = t311 * t386 - t347 * t312 + t338 * t388;
t234 = -qJDD(2) * pkin(2) - t352 * pkin(8) - t236;
t369 = -t306 * pkin(3) + t234 + (-t305 - t374) * qJ(4);
t390 = t349 ^ 2 * t352;
t401 = 2 * qJD(4);
t360 = -qJ(5) * t390 + qJDD(5) - t369 + (t313 + t401) * t380;
t399 = pkin(4) + pkin(9);
t400 = -pkin(3) - pkin(9);
t211 = t360 + t305 * pkin(5) + t399 * t306 + (pkin(5) * t349 + t346 * t400) * t378;
t304 = (pkin(5) * t346 + pkin(9) * t349) * qJD(2);
t214 = (-pkin(5) - qJ(4)) * t351 + (-pkin(4) * t385 - qJD(2) * t304) * t346 + (-pkin(3) - t399) * qJDD(3) + t370 + t405;
t345 = sin(qJ(6));
t348 = cos(qJ(6));
t207 = t211 * t348 - t214 * t345;
t298 = -qJD(3) * t348 + t345 * t379;
t246 = qJD(6) * t298 - qJDD(3) * t345 - t306 * t348;
t299 = -qJD(3) * t345 - t348 * t379;
t247 = -mrSges(7,1) * t298 + mrSges(7,2) * t299;
t324 = qJD(6) + t380;
t249 = -mrSges(7,2) * t324 + mrSges(7,3) * t298;
t294 = qJDD(6) + t305;
t204 = m(7) * t207 + mrSges(7,1) * t294 - mrSges(7,3) * t246 - t247 * t299 + t249 * t324;
t208 = t211 * t345 + t214 * t348;
t245 = -qJD(6) * t299 - qJDD(3) * t348 + t306 * t345;
t250 = mrSges(7,1) * t324 - mrSges(7,3) * t299;
t205 = m(7) * t208 - mrSges(7,2) * t294 + mrSges(7,3) * t245 + t247 * t298 - t250 * t324;
t384 = -t345 * t204 + t348 * t205;
t371 = -m(6) * t218 + t303 * t380 - t384;
t364 = -qJDD(3) * mrSges(6,2) - qJD(3) * t317 + t371;
t185 = -t305 * mrSges(6,3) - t364;
t328 = qJD(3) * t401;
t396 = t351 * pkin(3);
t407 = qJDD(3) * qJ(4) + t300 * t379 + t230;
t225 = t328 - t396 + t407;
t366 = pkin(4) * t390 + t306 * qJ(5) - t407;
t213 = qJDD(3) * pkin(5) + qJD(3) * t313 + t328 + t400 * t351 + (-0.2e1 * qJD(5) - t304) * t379 - t366;
t238 = Ifges(7,5) * t299 + Ifges(7,6) * t298 + Ifges(7,3) * t324;
t240 = Ifges(7,1) * t299 + Ifges(7,4) * t298 + Ifges(7,5) * t324;
t196 = -mrSges(7,1) * t213 + mrSges(7,3) * t208 + Ifges(7,4) * t246 + Ifges(7,2) * t245 + Ifges(7,6) * t294 - t238 * t299 + t240 * t324;
t239 = Ifges(7,4) * t299 + Ifges(7,2) * t298 + Ifges(7,6) * t324;
t197 = mrSges(7,2) * t213 - mrSges(7,3) * t207 + Ifges(7,1) * t246 + Ifges(7,4) * t245 + Ifges(7,5) * t294 + t238 * t298 - t239 * t324;
t402 = -2 * qJD(4);
t217 = 0.2e1 * t349 * t377 + t396 + (t402 - t313) * qJD(3) + t366;
t267 = -Ifges(6,6) * qJD(3) + (-Ifges(6,4) * t349 - Ifges(6,2) * t346) * qJD(2);
t270 = -Ifges(6,5) * qJD(3) + (-Ifges(6,1) * t349 - Ifges(6,4) * t346) * qJD(2);
t368 = -m(7) * t213 + t245 * mrSges(7,1) - t246 * mrSges(7,2) + t298 * t249 - t299 * t250;
t363 = mrSges(6,1) * t217 - mrSges(6,2) * t218 - Ifges(6,5) * t306 - Ifges(6,6) * t305 - Ifges(6,3) * qJDD(3) + pkin(5) * t368 + pkin(9) * t384 + t348 * t196 + t345 * t197 - t267 * t379 + t270 * t380;
t357 = -mrSges(5,1) * t227 + mrSges(5,3) * t225 + Ifges(5,4) * t305 + Ifges(5,2) * qJDD(3) - Ifges(5,6) * t306 - pkin(4) * t185 - t363;
t316 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t380;
t314 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t380;
t361 = -m(6) * t217 + qJDD(3) * mrSges(6,1) - t306 * mrSges(6,3) + qJD(3) * t314 - t368;
t358 = m(5) * t225 + qJDD(3) * mrSges(5,3) + qJD(3) * t316 + t301 * t379 + t361;
t375 = t303 * t379;
t271 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t346 - Ifges(5,5) * t349) * qJD(2);
t381 = t271 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t346 + Ifges(4,4) * t349) * qJD(2);
t319 = mrSges(5,2) * t379 + qJD(3) * mrSges(5,3);
t406 = -m(5) * t227 + qJDD(3) * mrSges(5,1) + qJD(3) * t319;
t408 = -((t265 - t269) * t346 + t349 * t381) * qJD(2) + mrSges(4,1) * t229 - mrSges(4,2) * t230 + Ifges(4,5) * t305 + Ifges(4,6) * t306 + Ifges(4,3) * qJDD(3) + pkin(3) * (-t301 * t380 + (-mrSges(5,2) + mrSges(6,3)) * t305 + t364 + t406) + qJ(4) * (t306 * mrSges(5,2) + t358 - t375) + t357;
t264 = -Ifges(6,3) * qJD(3) + (-Ifges(6,5) * t349 - Ifges(6,6) * t346) * qJD(2);
t404 = Ifges(6,4) * t306 + Ifges(6,2) * t305 - t264 * t379;
t302 = (-mrSges(4,1) * t349 + mrSges(4,2) * t346) * qJD(2);
t318 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t379;
t395 = mrSges(4,3) + mrSges(5,2);
t181 = m(4) * t229 + (mrSges(4,1) - mrSges(6,2)) * qJDD(3) + (-t317 + t318) * qJD(3) + (-t301 - t302) * t380 + (mrSges(6,3) - t395) * t305 + t371 + t406;
t315 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t380;
t192 = t358 - qJDD(3) * mrSges(4,2) - qJD(3) * t315 + m(4) * t230 + t395 * t306 + (t302 - t303) * t379;
t372 = -t181 * t346 + t349 * t192;
t175 = m(3) * t237 - mrSges(3,1) * t352 - qJDD(2) * mrSges(3,2) + t372;
t188 = t348 * t204 + t345 * t205;
t222 = -pkin(3) * t373 + t306 * pkin(4) + t360;
t184 = -m(6) * t222 - t305 * mrSges(6,1) + t306 * mrSges(6,2) - t314 * t380 + t317 * t379 - t188;
t228 = (pkin(3) * qJD(3) + t402) * t380 + t369;
t182 = m(5) * t228 - mrSges(5,1) * t306 - t305 * mrSges(5,3) - t316 * t380 - t319 * t379 + t184;
t355 = -m(4) * t234 + t306 * mrSges(4,1) - mrSges(4,2) * t305 - t315 * t380 + t318 * t379 - t182;
t179 = m(3) * t236 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t352 + t355;
t172 = t350 * t175 - t179 * t347;
t397 = pkin(7) * t172;
t394 = Ifges(5,6) - Ifges(6,5);
t177 = t349 * t181 + t346 * t192;
t268 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t346 - Ifges(5,6) * t349) * qJD(2);
t383 = -t264 + t268;
t176 = m(3) * t251 + t177;
t168 = t175 * t387 - t176 * t341 + t179 * t386;
t266 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t346 + Ifges(4,6) * t349) * qJD(2);
t362 = -mrSges(6,2) * t222 + mrSges(6,3) * t217 + Ifges(6,1) * t306 + Ifges(6,4) * t305 + pkin(9) * t188 - qJD(3) * t267 + t345 * t196 - t348 * t197;
t356 = mrSges(5,1) * t228 - mrSges(5,2) * t225 + pkin(4) * t184 + qJ(5) * (t361 - t375) - t362;
t164 = -t356 + (Ifges(5,3) + Ifges(4,2)) * t306 + (Ifges(4,4) - Ifges(5,5)) * t305 + (Ifges(4,6) - t394) * qJDD(3) + t381 * qJD(3) + mrSges(4,3) * t230 - mrSges(4,1) * t234 - pkin(3) * t182 + (-t266 - t383) * t380;
t365 = mrSges(7,1) * t207 - mrSges(7,2) * t208 + Ifges(7,5) * t246 + Ifges(7,6) * t245 + Ifges(7,3) * t294 + t299 * t239 - t298 * t240;
t359 = mrSges(6,1) * t222 - mrSges(6,3) * t218 + Ifges(6,6) * qJDD(3) + pkin(5) * t188 + qJD(3) * t270 + t365;
t354 = mrSges(5,2) * t227 - mrSges(5,3) * t228 + Ifges(5,1) * t305 + Ifges(5,4) * qJDD(3) - Ifges(5,5) * t306 - qJ(5) * t185 + qJD(3) * t265 + t268 * t379 + t359;
t169 = t354 + (-t264 + t266) * t379 + Ifges(4,5) * qJDD(3) + (Ifges(6,4) + Ifges(4,4)) * t306 + (Ifges(4,1) + Ifges(6,2)) * t305 - qJD(3) * t269 - mrSges(4,3) * t229 + mrSges(4,2) * t234 - qJ(4) * t182;
t159 = mrSges(3,1) * t236 - mrSges(3,2) * t237 + Ifges(3,3) * qJDD(2) + pkin(2) * t355 + pkin(8) * t372 + t349 * t164 + t346 * t169;
t161 = mrSges(3,2) * t251 - mrSges(3,3) * t236 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t352 - pkin(8) * t177 - t164 * t346 + t169 * t349;
t163 = -mrSges(3,1) * t251 + mrSges(3,3) * t237 + t352 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t177 - t408;
t367 = mrSges(2,1) * t311 - mrSges(2,2) * t312 + pkin(1) * t168 + t343 * t159 + t161 * t389 + t163 * t388 + t341 * t397;
t170 = m(2) * t312 + t172;
t167 = t343 * t176 + (t175 * t347 + t179 * t350) * t341;
t165 = m(2) * t311 + t168;
t157 = mrSges(2,2) * t338 - mrSges(2,3) * t311 + t350 * t161 - t347 * t163 + (-t167 * t341 - t168 * t343) * pkin(7);
t156 = -mrSges(2,1) * t338 + mrSges(2,3) * t312 - pkin(1) * t167 - t341 * t159 + (t161 * t347 + t163 * t350 + t397) * t343;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t342 * t157 - t340 * t156 - qJ(1) * (t165 * t342 + t170 * t340), t157, t161, t169, t354 + t404, -Ifges(6,5) * qJDD(3) - t264 * t380 - t362, t197; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t340 * t157 + t342 * t156 + qJ(1) * (-t165 * t340 + t170 * t342), t156, t163, t164, t357 + (-t346 * t265 - t349 * t271) * qJD(2), -t359 - t404, t196; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t367, t367, t159, t408, Ifges(5,5) * t305 - Ifges(5,3) * t306 - qJD(3) * t271 + qJDD(3) * t394 + t380 * t383 + t356, t363, t365;];
m_new  = t1;
