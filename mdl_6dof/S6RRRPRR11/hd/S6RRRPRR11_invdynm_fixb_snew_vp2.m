% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 14:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 14:22:08
% EndTime: 2019-05-07 14:23:07
% DurationCPUTime: 27.04s
% Computational Cost: add. (490869->396), mult. (1051718->496), div. (0->0), fcn. (822660->12), ass. (0->160)
t367 = sin(qJ(2));
t371 = cos(qJ(2));
t362 = sin(pkin(6));
t392 = qJD(1) * t362;
t340 = (-pkin(2) * t371 - pkin(9) * t367) * t392;
t363 = cos(pkin(6));
t358 = t363 * qJD(1) + qJD(2);
t356 = t358 ^ 2;
t357 = t363 * qJDD(1) + qJDD(2);
t391 = qJD(1) * t371;
t368 = sin(qJ(1));
t372 = cos(qJ(1));
t351 = t368 * g(1) - t372 * g(2);
t373 = qJD(1) ^ 2;
t404 = pkin(8) * t362;
t337 = qJDD(1) * pkin(1) + t373 * t404 + t351;
t352 = -t372 * g(1) - t368 * g(2);
t390 = qJDD(1) * t362;
t338 = -t373 * pkin(1) + pkin(8) * t390 + t352;
t398 = t363 * t367;
t393 = t337 * t398 + t371 * t338;
t269 = -t356 * pkin(2) + t357 * pkin(9) + (-g(3) * t367 + t340 * t391) * t362 + t393;
t341 = (qJD(2) * t391 + qJDD(1) * t367) * t362;
t389 = t367 * t392;
t342 = -qJD(2) * t389 + t371 * t390;
t403 = t363 * g(3);
t270 = -t342 * pkin(2) - t341 * pkin(9) - t403 + (-t337 + (pkin(2) * t367 - pkin(9) * t371) * t358 * qJD(1)) * t362;
t366 = sin(qJ(3));
t406 = cos(qJ(3));
t240 = -t366 * t269 + t406 * t270;
t241 = t406 * t269 + t366 * t270;
t324 = -t406 * t358 + t366 * t389;
t325 = t366 * t358 + t406 * t389;
t388 = t362 * t391;
t350 = qJD(3) - t388;
t279 = Ifges(5,5) * t325 + Ifges(5,6) * t350 + Ifges(5,3) * t324;
t282 = Ifges(4,4) * t325 - Ifges(4,2) * t324 + Ifges(4,6) * t350;
t284 = Ifges(4,1) * t325 - Ifges(4,4) * t324 + Ifges(4,5) * t350;
t296 = t325 * qJD(3) + t366 * t341 - t406 * t357;
t297 = -t324 * qJD(3) + t406 * t341 + t366 * t357;
t304 = t324 * mrSges(5,1) - t325 * mrSges(5,3);
t334 = qJDD(3) - t342;
t303 = t324 * pkin(3) - t325 * qJ(4);
t349 = t350 ^ 2;
t235 = -t334 * pkin(3) - t349 * qJ(4) + t325 * t303 + qJDD(4) - t240;
t401 = t324 * t350;
t226 = (-t297 - t401) * pkin(10) + (t324 * t325 - t334) * pkin(4) + t235;
t407 = 2 * qJD(4);
t233 = -t349 * pkin(3) + t334 * qJ(4) - t324 * t303 + t350 * t407 + t241;
t312 = -t350 * pkin(4) - t325 * pkin(10);
t323 = t324 ^ 2;
t229 = -t323 * pkin(4) + t296 * pkin(10) + t350 * t312 + t233;
t365 = sin(qJ(5));
t370 = cos(qJ(5));
t224 = t365 * t226 + t370 * t229;
t302 = t365 * t324 + t370 * t325;
t251 = -t302 * qJD(5) + t370 * t296 - t365 * t297;
t301 = t370 * t324 - t365 * t325;
t263 = -t301 * mrSges(6,1) + t302 * mrSges(6,2);
t345 = qJD(5) - t350;
t278 = t345 * mrSges(6,1) - t302 * mrSges(6,3);
t333 = qJDD(5) - t334;
t264 = -t301 * pkin(5) - t302 * pkin(11);
t344 = t345 ^ 2;
t219 = -t344 * pkin(5) + t333 * pkin(11) + t301 * t264 + t224;
t397 = t363 * t371;
t399 = t362 * t371;
t298 = -g(3) * t399 + t337 * t397 - t367 * t338;
t268 = -t357 * pkin(2) - t356 * pkin(9) + t340 * t389 - t298;
t384 = t296 * pkin(3) + t268 + (-t297 + t401) * qJ(4);
t405 = pkin(3) * t350;
t230 = -t296 * pkin(4) - t323 * pkin(10) - t384 + (t312 - t405 + t407) * t325;
t252 = t301 * qJD(5) + t365 * t296 + t370 * t297;
t220 = t230 + (-t301 * t345 - t252) * pkin(11) + (t302 * t345 - t251) * pkin(5);
t364 = sin(qJ(6));
t369 = cos(qJ(6));
t216 = -t364 * t219 + t369 * t220;
t275 = -t364 * t302 + t369 * t345;
t239 = t275 * qJD(6) + t369 * t252 + t364 * t333;
t250 = qJDD(6) - t251;
t276 = t369 * t302 + t364 * t345;
t255 = -t275 * mrSges(7,1) + t276 * mrSges(7,2);
t300 = qJD(6) - t301;
t256 = -t300 * mrSges(7,2) + t275 * mrSges(7,3);
t212 = m(7) * t216 + t250 * mrSges(7,1) - t239 * mrSges(7,3) - t276 * t255 + t300 * t256;
t217 = t369 * t219 + t364 * t220;
t238 = -t276 * qJD(6) - t364 * t252 + t369 * t333;
t257 = t300 * mrSges(7,1) - t276 * mrSges(7,3);
t213 = m(7) * t217 - t250 * mrSges(7,2) + t238 * mrSges(7,3) + t275 * t255 - t300 * t257;
t386 = -t364 * t212 + t369 * t213;
t199 = m(6) * t224 - t333 * mrSges(6,2) + t251 * mrSges(6,3) + t301 * t263 - t345 * t278 + t386;
t223 = t370 * t226 - t365 * t229;
t277 = -t345 * mrSges(6,2) + t301 * mrSges(6,3);
t218 = -t333 * pkin(5) - t344 * pkin(11) + t302 * t264 - t223;
t383 = -m(7) * t218 + t238 * mrSges(7,1) - t239 * mrSges(7,2) + t275 * t256 - t276 * t257;
t208 = m(6) * t223 + t333 * mrSges(6,1) - t252 * mrSges(6,3) - t302 * t263 + t345 * t277 + t383;
t193 = t365 * t199 + t370 * t208;
t283 = Ifges(5,1) * t325 + Ifges(5,4) * t350 + Ifges(5,5) * t324;
t242 = Ifges(7,5) * t276 + Ifges(7,6) * t275 + Ifges(7,3) * t300;
t244 = Ifges(7,1) * t276 + Ifges(7,4) * t275 + Ifges(7,5) * t300;
t206 = -mrSges(7,1) * t218 + mrSges(7,3) * t217 + Ifges(7,4) * t239 + Ifges(7,2) * t238 + Ifges(7,6) * t250 - t276 * t242 + t300 * t244;
t243 = Ifges(7,4) * t276 + Ifges(7,2) * t275 + Ifges(7,6) * t300;
t207 = mrSges(7,2) * t218 - mrSges(7,3) * t216 + Ifges(7,1) * t239 + Ifges(7,4) * t238 + Ifges(7,5) * t250 + t275 * t242 - t300 * t243;
t259 = Ifges(6,4) * t302 + Ifges(6,2) * t301 + Ifges(6,6) * t345;
t260 = Ifges(6,1) * t302 + Ifges(6,4) * t301 + Ifges(6,5) * t345;
t380 = mrSges(6,1) * t223 - mrSges(6,2) * t224 + Ifges(6,5) * t252 + Ifges(6,6) * t251 + Ifges(6,3) * t333 + pkin(5) * t383 + pkin(11) * t386 + t369 * t206 + t364 * t207 + t302 * t259 - t301 * t260;
t376 = mrSges(5,1) * t235 - mrSges(5,3) * t233 - Ifges(5,4) * t297 - Ifges(5,2) * t334 - Ifges(5,6) * t296 + pkin(4) * t193 - t324 * t283 + t380;
t311 = -t324 * mrSges(5,2) + t350 * mrSges(5,3);
t382 = -m(5) * t235 + t334 * mrSges(5,1) + t350 * t311 - t193;
t194 = t370 * t199 - t365 * t208;
t310 = -t350 * mrSges(5,1) + t325 * mrSges(5,2);
t385 = m(5) * t233 + t334 * mrSges(5,3) + t350 * t310 + t194;
t408 = (t282 - t279) * t325 + mrSges(4,1) * t240 - mrSges(4,2) * t241 + Ifges(4,5) * t297 - Ifges(4,6) * t296 + Ifges(4,3) * t334 + pkin(3) * (-t297 * mrSges(5,2) - t325 * t304 + t382) + qJ(4) * (-t296 * mrSges(5,2) - t324 * t304 + t385) + t324 * t284 - t376;
t402 = -mrSges(4,3) - mrSges(5,2);
t400 = t362 * t367;
t309 = t350 * mrSges(4,1) - t325 * mrSges(4,3);
t394 = -t324 * mrSges(4,1) - t325 * mrSges(4,2) - t304;
t188 = m(4) * t241 - t334 * mrSges(4,2) + t402 * t296 - t350 * t309 + t394 * t324 + t385;
t308 = -t350 * mrSges(4,2) - t324 * mrSges(4,3);
t190 = m(4) * t240 + t334 * mrSges(4,1) + t402 * t297 + t350 * t308 + t394 * t325 + t382;
t182 = t366 * t188 + t406 * t190;
t202 = t369 * t212 + t364 * t213;
t281 = Ifges(5,4) * t325 + Ifges(5,2) * t350 + Ifges(5,6) * t324;
t396 = -Ifges(4,5) * t325 + Ifges(4,6) * t324 - Ifges(4,3) * t350 - t281;
t299 = -g(3) * t400 + t393;
t335 = t358 * mrSges(3,1) - mrSges(3,3) * t389;
t339 = (-mrSges(3,1) * t371 + mrSges(3,2) * t367) * t392;
t387 = t406 * t188 - t366 * t190;
t180 = m(3) * t299 - t357 * mrSges(3,2) + t342 * mrSges(3,3) - t358 * t335 + t339 * t388 + t387;
t336 = -t358 * mrSges(3,2) + mrSges(3,3) * t388;
t200 = -m(6) * t230 + t251 * mrSges(6,1) - t252 * mrSges(6,2) + t301 * t277 - t302 * t278 - t202;
t236 = (-(2 * qJD(4)) + t405) * t325 + t384;
t197 = m(5) * t236 + t296 * mrSges(5,1) - t297 * mrSges(5,3) - t325 * t310 + t324 * t311 + t200;
t375 = -m(4) * t268 - t296 * mrSges(4,1) - t297 * mrSges(4,2) - t324 * t308 - t325 * t309 - t197;
t196 = m(3) * t298 + t357 * mrSges(3,1) - t341 * mrSges(3,3) + t358 * t336 - t339 * t389 + t375;
t177 = t371 * t180 - t367 * t196;
t317 = -t362 * t337 - t403;
t181 = m(3) * t317 - t342 * mrSges(3,1) + t341 * mrSges(3,2) + (t335 * t367 - t336 * t371) * t392 + t182;
t173 = t180 * t398 - t362 * t181 + t196 * t397;
t258 = Ifges(6,5) * t302 + Ifges(6,6) * t301 + Ifges(6,3) * t345;
t184 = mrSges(6,2) * t230 - mrSges(6,3) * t223 + Ifges(6,1) * t252 + Ifges(6,4) * t251 + Ifges(6,5) * t333 - pkin(11) * t202 - t364 * t206 + t369 * t207 + t301 * t258 - t345 * t259;
t377 = mrSges(7,1) * t216 - mrSges(7,2) * t217 + Ifges(7,5) * t239 + Ifges(7,6) * t238 + Ifges(7,3) * t250 + t276 * t243 - t275 * t244;
t185 = -mrSges(6,1) * t230 + mrSges(6,3) * t224 + Ifges(6,4) * t252 + Ifges(6,2) * t251 + Ifges(6,6) * t333 - pkin(5) * t202 - t302 * t258 + t345 * t260 - t377;
t378 = -mrSges(5,1) * t236 + mrSges(5,2) * t233 - pkin(4) * t200 - pkin(10) * t194 - t365 * t184 - t370 * t185;
t169 = -mrSges(4,1) * t268 + mrSges(4,3) * t241 - pkin(3) * t197 + (t284 + t283) * t350 + (Ifges(4,6) - Ifges(5,6)) * t334 + t396 * t325 + (Ifges(4,4) - Ifges(5,5)) * t297 + (-Ifges(4,2) - Ifges(5,3)) * t296 + t378;
t379 = mrSges(5,2) * t235 - mrSges(5,3) * t236 + Ifges(5,1) * t297 + Ifges(5,4) * t334 + Ifges(5,5) * t296 - pkin(10) * t193 + t370 * t184 - t365 * t185 + t350 * t279;
t174 = mrSges(4,2) * t268 - mrSges(4,3) * t240 + Ifges(4,1) * t297 - Ifges(4,4) * t296 + Ifges(4,5) * t334 - qJ(4) * t197 - t350 * t282 + t396 * t324 + t379;
t315 = Ifges(3,6) * t358 + (Ifges(3,4) * t367 + Ifges(3,2) * t371) * t392;
t316 = Ifges(3,5) * t358 + (Ifges(3,1) * t367 + Ifges(3,4) * t371) * t392;
t164 = Ifges(3,5) * t341 + Ifges(3,6) * t342 + Ifges(3,3) * t357 + mrSges(3,1) * t298 - mrSges(3,2) * t299 + t366 * t174 + t406 * t169 + pkin(2) * t375 + pkin(9) * t387 + (t315 * t367 - t316 * t371) * t392;
t314 = Ifges(3,3) * t358 + (Ifges(3,5) * t367 + Ifges(3,6) * t371) * t392;
t166 = mrSges(3,2) * t317 - mrSges(3,3) * t298 + Ifges(3,1) * t341 + Ifges(3,4) * t342 + Ifges(3,5) * t357 - pkin(9) * t182 - t366 * t169 + t406 * t174 + t314 * t388 - t358 * t315;
t168 = -mrSges(3,1) * t317 + mrSges(3,3) * t299 + Ifges(3,4) * t341 + Ifges(3,2) * t342 + Ifges(3,6) * t357 - pkin(2) * t182 - t314 * t389 + t358 * t316 - t408;
t381 = mrSges(2,1) * t351 - mrSges(2,2) * t352 + Ifges(2,3) * qJDD(1) + pkin(1) * t173 + t363 * t164 + t166 * t400 + t168 * t399 + t177 * t404;
t175 = m(2) * t352 - t373 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t177;
t172 = t363 * t181 + (t180 * t367 + t196 * t371) * t362;
t170 = m(2) * t351 + qJDD(1) * mrSges(2,1) - t373 * mrSges(2,2) + t173;
t162 = -mrSges(2,2) * g(3) - mrSges(2,3) * t351 + Ifges(2,5) * qJDD(1) - t373 * Ifges(2,6) + t371 * t166 - t367 * t168 + (-t172 * t362 - t173 * t363) * pkin(8);
t161 = mrSges(2,1) * g(3) + mrSges(2,3) * t352 + t373 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t172 - t362 * t164 + (pkin(8) * t177 + t166 * t367 + t168 * t371) * t363;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t372 * t162 - t368 * t161 - pkin(7) * (t372 * t170 + t368 * t175), t162, t166, t174, -t324 * t281 + t379, t184, t207; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t368 * t162 + t372 * t161 + pkin(7) * (-t368 * t170 + t372 * t175), t161, t168, t169, -t325 * t279 - t376, t185, t206; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t381, t381, t164, t408, Ifges(5,5) * t297 + Ifges(5,6) * t334 + Ifges(5,3) * t296 + t325 * t281 - t350 * t283 - t378, t380, t377;];
m_new  = t1;
