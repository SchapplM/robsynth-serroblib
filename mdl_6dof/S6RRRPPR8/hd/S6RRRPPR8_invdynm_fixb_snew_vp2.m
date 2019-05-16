% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-05-07 06:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:09:28
% EndTime: 2019-05-07 06:09:53
% DurationCPUTime: 11.56s
% Computational Cost: add. (188756->405), mult. (407343->480), div. (0->0), fcn. (303932->10), ass. (0->160)
t380 = cos(pkin(6));
t375 = qJD(1) * t380 + qJD(2);
t382 = sin(qJ(3));
t383 = sin(qJ(2));
t379 = sin(pkin(6));
t414 = qJD(1) * t379;
t411 = t383 * t414;
t434 = cos(qJ(3));
t340 = t382 * t375 + t434 * t411;
t386 = cos(qJ(2));
t413 = qJD(1) * t386;
t410 = t379 * t413;
t366 = -qJD(3) + t410;
t322 = pkin(4) * t366 - qJ(5) * t340;
t443 = (2 * qJD(4)) + t322;
t438 = -2 * qJD(4);
t360 = t366 * t438;
t357 = (-pkin(2) * t386 - pkin(9) * t383) * t414;
t373 = t375 ^ 2;
t374 = qJDD(1) * t380 + qJDD(2);
t384 = sin(qJ(1));
t387 = cos(qJ(1));
t368 = t384 * g(1) - g(2) * t387;
t388 = qJD(1) ^ 2;
t432 = pkin(8) * t379;
t354 = qJDD(1) * pkin(1) + t388 * t432 + t368;
t369 = -g(1) * t387 - g(2) * t384;
t412 = qJDD(1) * t379;
t355 = -pkin(1) * t388 + pkin(8) * t412 + t369;
t422 = t380 * t383;
t415 = t354 * t422 + t386 * t355;
t260 = -t373 * pkin(2) + t374 * pkin(9) + (-g(3) * t383 + t357 * t413) * t379 + t415;
t358 = (qJD(2) * t413 + qJDD(1) * t383) * t379;
t359 = -qJD(2) * t411 + t386 * t412;
t431 = t380 * g(3);
t261 = -t359 * pkin(2) - t358 * pkin(9) - t431 + (-t354 + (pkin(2) * t383 - pkin(9) * t386) * t375 * qJD(1)) * t379;
t245 = t434 * t260 + t382 * t261;
t339 = -t434 * t375 + t382 * t411;
t310 = pkin(3) * t339 - qJ(4) * t340;
t351 = -qJDD(3) + t359;
t439 = t366 ^ 2;
t407 = pkin(3) * t439 + t351 * qJ(4) + t339 * t310 - t245;
t240 = t360 - t407;
t324 = mrSges(5,1) * t366 + mrSges(5,2) * t340;
t442 = m(5) * t240 - t351 * mrSges(5,3) - t366 * t324;
t325 = -mrSges(6,2) * t366 - mrSges(6,3) * t340;
t305 = qJD(3) * t340 + t358 * t382 - t434 * t374;
t338 = t339 ^ 2;
t399 = t338 * pkin(4) - t305 * qJ(5) + t407;
t436 = -2 * qJD(5);
t233 = t339 * t436 + t366 * t443 + t399;
t309 = mrSges(6,1) * t340 + mrSges(6,2) * t339;
t313 = pkin(5) * t340 - pkin(10) * t339;
t229 = -t351 * pkin(5) - t439 * pkin(10) - t366 * t322 + t360 + ((2 * qJD(5)) + t313) * t339 - t399;
t381 = sin(qJ(6));
t385 = cos(qJ(6));
t319 = t339 * t385 + t366 * t381;
t250 = -qJD(6) * t319 - t305 * t381 + t351 * t385;
t318 = -t339 * t381 + t366 * t385;
t251 = qJD(6) * t318 + t305 * t385 + t351 * t381;
t337 = qJD(6) + t340;
t272 = -mrSges(7,2) * t337 + mrSges(7,3) * t318;
t273 = mrSges(7,1) * t337 - mrSges(7,3) * t319;
t406 = -m(7) * t229 + t250 * mrSges(7,1) - t251 * mrSges(7,2) + t318 * t272 - t319 * t273;
t398 = -m(6) * t233 + t305 * mrSges(6,3) + t339 * t309 - t406;
t214 = -t351 * mrSges(6,1) - t366 * t325 + t398;
t244 = -t382 * t260 + t434 * t261;
t276 = Ifges(5,5) * t340 - Ifges(5,6) * t366 + Ifges(5,3) * t339;
t280 = Ifges(4,4) * t340 - Ifges(4,2) * t339 - Ifges(4,6) * t366;
t283 = Ifges(4,1) * t340 - Ifges(4,4) * t339 - Ifges(4,5) * t366;
t306 = -t339 * qJD(3) + t434 * t358 + t382 * t374;
t311 = mrSges(5,1) * t339 - mrSges(5,3) * t340;
t403 = -qJ(4) * t439 + t340 * t310 + qJDD(4) - t244;
t426 = t339 * t366;
t394 = (-t306 + t426) * qJ(5) + t403 + (t339 * pkin(4) + t436) * t340;
t231 = (pkin(3) + pkin(4)) * t351 + t394;
t320 = mrSges(6,1) * t366 - mrSges(6,3) * t339;
t421 = t380 * t386;
t423 = t379 * t386;
t307 = -g(3) * t423 + t354 * t421 - t383 * t355;
t259 = -t374 * pkin(2) - t373 * pkin(9) + t357 * t411 - t307;
t405 = t305 * pkin(3) + t259 + (-t306 - t426) * qJ(4);
t396 = -t338 * qJ(5) + t340 * t443 + qJDD(5) - t405;
t435 = -pkin(4) - pkin(10);
t226 = t396 + t306 * pkin(5) + (pkin(5) * t339 + (pkin(3) + pkin(10)) * t340) * t366 + t435 * t305;
t227 = -t439 * pkin(5) - t340 * t313 + (pkin(3) - t435) * t351 + t394;
t222 = t226 * t385 - t227 * t381;
t262 = -mrSges(7,1) * t318 + mrSges(7,2) * t319;
t298 = qJDD(6) + t306;
t219 = m(7) * t222 + mrSges(7,1) * t298 - mrSges(7,3) * t251 - t262 * t319 + t272 * t337;
t223 = t226 * t381 + t227 * t385;
t220 = m(7) * t223 - mrSges(7,2) * t298 + mrSges(7,3) * t250 + t262 * t318 - t273 * t337;
t420 = -t381 * t219 + t385 * t220;
t408 = -m(6) * t231 + t351 * mrSges(6,2) + t366 * t320 - t420;
t199 = -t306 * mrSges(6,3) - t340 * t309 - t408;
t242 = t351 * pkin(3) + t403;
t282 = Ifges(5,1) * t340 - Ifges(5,4) * t366 + Ifges(5,5) * t339;
t252 = Ifges(7,5) * t319 + Ifges(7,6) * t318 + Ifges(7,3) * t337;
t254 = Ifges(7,1) * t319 + Ifges(7,4) * t318 + Ifges(7,5) * t337;
t211 = -mrSges(7,1) * t229 + mrSges(7,3) * t223 + Ifges(7,4) * t251 + Ifges(7,2) * t250 + Ifges(7,6) * t298 - t252 * t319 + t254 * t337;
t253 = Ifges(7,4) * t319 + Ifges(7,2) * t318 + Ifges(7,6) * t337;
t212 = mrSges(7,2) * t229 - mrSges(7,3) * t222 + Ifges(7,1) * t251 + Ifges(7,4) * t250 + Ifges(7,5) * t298 + t252 * t318 - t253 * t337;
t278 = Ifges(6,4) * t339 - Ifges(6,2) * t340 + Ifges(6,6) * t366;
t281 = Ifges(6,1) * t339 - Ifges(6,4) * t340 + Ifges(6,5) * t366;
t400 = mrSges(6,1) * t233 - mrSges(6,2) * t231 + Ifges(6,5) * t305 - Ifges(6,6) * t306 + Ifges(6,3) * t351 + pkin(5) * t406 + pkin(10) * t420 + t385 * t211 + t381 * t212 + t339 * t278 + t340 * t281;
t393 = mrSges(5,1) * t242 - mrSges(5,3) * t240 - Ifges(5,4) * t306 + Ifges(5,2) * t351 - Ifges(5,6) * t305 + pkin(4) * t199 - t339 * t282 + t400;
t326 = -mrSges(5,2) * t339 - mrSges(5,3) * t366;
t401 = -m(5) * t242 - t351 * mrSges(5,1) - t366 * t326 + t408;
t441 = -(t276 - t280) * t340 + mrSges(4,1) * t244 - mrSges(4,2) * t245 + Ifges(4,5) * t306 - Ifges(4,6) * t305 - Ifges(4,3) * t351 + pkin(3) * ((t309 - t311) * t340 + (-mrSges(5,2) + mrSges(6,3)) * t306 + t401) + qJ(4) * (-t305 * mrSges(5,2) - t339 * t311 + t214 + t442) + t339 * t283 - t393;
t440 = Ifges(6,4) * t305 - Ifges(6,2) * t306 + Ifges(6,6) * t351 + t366 * t281;
t433 = pkin(3) * t366;
t430 = -mrSges(4,3) - mrSges(5,2);
t424 = t379 * t383;
t321 = mrSges(4,2) * t366 - mrSges(4,3) * t339;
t416 = -mrSges(4,1) * t339 - mrSges(4,2) * t340 - t311;
t196 = m(4) * t244 - t351 * mrSges(4,1) - t366 * t321 + (t309 + t416) * t340 + (mrSges(6,3) + t430) * t306 + t401;
t323 = -mrSges(4,1) * t366 - mrSges(4,3) * t340;
t207 = t398 + m(4) * t245 + (t323 - t325) * t366 + (mrSges(4,2) - mrSges(6,1)) * t351 + t416 * t339 + t430 * t305 + t442;
t192 = t434 * t196 + t382 * t207;
t203 = t385 * t219 + t381 * t220;
t279 = Ifges(5,4) * t340 - Ifges(5,2) * t366 + Ifges(5,6) * t339;
t418 = -Ifges(4,5) * t340 + Ifges(4,6) * t339 + Ifges(4,3) * t366 - t279;
t417 = -t278 + t282;
t308 = -g(3) * t424 + t415;
t352 = mrSges(3,1) * t375 - mrSges(3,3) * t411;
t356 = (-mrSges(3,1) * t386 + mrSges(3,2) * t383) * t414;
t409 = -t196 * t382 + t434 * t207;
t190 = m(3) * t308 - mrSges(3,2) * t374 + mrSges(3,3) * t359 - t352 * t375 + t356 * t410 + t409;
t353 = -mrSges(3,2) * t375 + mrSges(3,3) * t410;
t237 = -t305 * pkin(4) + t340 * t433 + t396;
t200 = -m(6) * t237 - t306 * mrSges(6,1) - t305 * mrSges(6,2) - t339 * t320 - t340 * t325 - t203;
t243 = (t438 - t433) * t340 + t405;
t197 = m(5) * t243 + t305 * mrSges(5,1) - t306 * mrSges(5,3) - t340 * t324 + t339 * t326 + t200;
t390 = -m(4) * t259 - t305 * mrSges(4,1) - t306 * mrSges(4,2) - t339 * t321 - t340 * t323 - t197;
t194 = m(3) * t307 + t374 * mrSges(3,1) - t358 * mrSges(3,3) + t375 * t353 - t356 * t411 + t390;
t187 = t386 * t190 - t194 * t383;
t331 = -t379 * t354 - t431;
t191 = m(3) * t331 - t359 * mrSges(3,1) + t358 * mrSges(3,2) + (t352 * t383 - t353 * t386) * t414 + t192;
t183 = t190 * t422 - t191 * t379 + t194 * t421;
t275 = Ifges(6,5) * t339 - Ifges(6,6) * t340 + Ifges(6,3) * t366;
t397 = -mrSges(6,2) * t237 + mrSges(6,3) * t233 - Ifges(6,1) * t305 + Ifges(6,4) * t306 - Ifges(6,5) * t351 + pkin(10) * t203 + t381 * t211 - t385 * t212 + t340 * t275;
t392 = mrSges(5,1) * t243 - mrSges(5,2) * t240 + pkin(4) * t200 + qJ(5) * t214 - t397;
t179 = -t392 + (-t283 - t417) * t366 + (-Ifges(4,6) + Ifges(5,6)) * t351 + t418 * t340 + (Ifges(4,4) - Ifges(5,5)) * t306 + (-Ifges(5,3) - Ifges(4,2)) * t305 - mrSges(4,1) * t259 + mrSges(4,3) * t245 - pkin(3) * t197;
t402 = mrSges(7,1) * t222 - mrSges(7,2) * t223 + Ifges(7,5) * t251 + Ifges(7,6) * t250 + Ifges(7,3) * t298 + t319 * t253 - t318 * t254;
t395 = mrSges(6,1) * t237 - mrSges(6,3) * t231 + pkin(5) * t203 + t339 * t275 + t402;
t391 = mrSges(5,2) * t242 - mrSges(5,3) * t243 + Ifges(5,1) * t306 - Ifges(5,4) * t351 + Ifges(5,5) * t305 - qJ(5) * t199 - t366 * t276 + t395;
t184 = (t280 - t281) * t366 + (-Ifges(6,6) - Ifges(4,5)) * t351 + t418 * t339 + (Ifges(4,1) + Ifges(6,2)) * t306 + (-Ifges(6,4) - Ifges(4,4)) * t305 + mrSges(4,2) * t259 - mrSges(4,3) * t244 + t391 - qJ(4) * t197;
t329 = Ifges(3,6) * t375 + (Ifges(3,4) * t383 + Ifges(3,2) * t386) * t414;
t330 = Ifges(3,5) * t375 + (Ifges(3,1) * t383 + Ifges(3,4) * t386) * t414;
t174 = Ifges(3,5) * t358 + Ifges(3,6) * t359 + Ifges(3,3) * t374 + mrSges(3,1) * t307 - mrSges(3,2) * t308 + t382 * t184 + t434 * t179 + pkin(2) * t390 + pkin(9) * t409 + (t329 * t383 - t330 * t386) * t414;
t328 = Ifges(3,3) * t375 + (Ifges(3,5) * t383 + Ifges(3,6) * t386) * t414;
t176 = mrSges(3,2) * t331 - mrSges(3,3) * t307 + Ifges(3,1) * t358 + Ifges(3,4) * t359 + Ifges(3,5) * t374 - pkin(9) * t192 - t382 * t179 + t434 * t184 + t328 * t410 - t375 * t329;
t178 = -mrSges(3,1) * t331 + mrSges(3,3) * t308 + Ifges(3,4) * t358 + Ifges(3,2) * t359 + Ifges(3,6) * t374 - pkin(2) * t192 - t328 * t411 + t375 * t330 - t441;
t404 = mrSges(2,1) * t368 - mrSges(2,2) * t369 + Ifges(2,3) * qJDD(1) + pkin(1) * t183 + t380 * t174 + t176 * t424 + t178 * t423 + t187 * t432;
t185 = m(2) * t369 - mrSges(2,1) * t388 - qJDD(1) * mrSges(2,2) + t187;
t182 = t380 * t191 + (t190 * t383 + t194 * t386) * t379;
t180 = m(2) * t368 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t388 + t183;
t172 = -mrSges(2,2) * g(3) - mrSges(2,3) * t368 + Ifges(2,5) * qJDD(1) - t388 * Ifges(2,6) + t386 * t176 - t383 * t178 + (-t182 * t379 - t183 * t380) * pkin(8);
t171 = mrSges(2,1) * g(3) + mrSges(2,3) * t369 + t388 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t182 - t379 * t174 + (pkin(8) * t187 + t176 * t383 + t178 * t386) * t380;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t387 * t172 - t384 * t171 - pkin(7) * (t180 * t387 + t185 * t384), t172, t176, t184, -t339 * t279 + t391 - t440, -t366 * t278 - t397, t212; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t384 * t172 + t387 * t171 + pkin(7) * (-t180 * t384 + t185 * t387), t171, t178, t179, -t340 * t276 - t393, -t395 + t440, t211; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t404, t404, t174, t441, Ifges(5,5) * t306 - Ifges(5,6) * t351 + Ifges(5,3) * t305 + t340 * t279 + t417 * t366 + t392, t400, t402;];
m_new  = t1;
