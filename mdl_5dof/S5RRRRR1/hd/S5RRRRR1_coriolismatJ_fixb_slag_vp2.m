% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:59
% EndTime: 2019-03-08 18:37:08
% DurationCPUTime: 5.02s
% Computational Cost: add. (9863->428), mult. (21102->564), div. (0->0), fcn. (24558->8), ass. (0->257)
t230 = sin(qJ(5));
t232 = cos(qJ(5));
t379 = Ifges(6,4) * t232;
t380 = Ifges(6,4) * t230;
t393 = t232 / 0.2e1;
t395 = -t230 / 0.2e1;
t447 = (Ifges(6,2) * t232 + t380) * t395 + (Ifges(6,1) * t230 + t379) * t393;
t388 = sin(qJ(3));
t389 = sin(qJ(2));
t391 = cos(qJ(3));
t392 = cos(qJ(2));
t212 = t388 * t392 + t389 * t391;
t231 = sin(qJ(4));
t254 = t388 * t389 - t391 * t392;
t390 = cos(qJ(4));
t171 = t212 * t231 + t390 * t254;
t279 = -Ifges(6,2) * t230 + t379;
t241 = -t212 * t390 + t231 * t254;
t376 = Ifges(6,6) * t241;
t47 = t171 * t279 + t376;
t449 = t232 * t47;
t281 = Ifges(6,1) * t232 - t380;
t378 = Ifges(6,5) * t241;
t49 = t171 * t281 + t378;
t450 = t230 * t49;
t452 = t449 / 0.2e1 + t450 / 0.2e1;
t451 = t449 / 0.4e1 + t450 / 0.4e1;
t422 = pkin(4) / 0.2e1;
t377 = Ifges(6,5) * t171;
t327 = -t377 / 0.2e1;
t365 = t230 * mrSges(6,1);
t384 = mrSges(6,2) * t232;
t282 = t365 + t384;
t448 = t282 * t171;
t318 = t388 * t231;
t204 = (t390 * t391 - t318) * pkin(2);
t228 = t230 ^ 2;
t229 = t232 ^ 2;
t340 = t228 + t229;
t297 = t340 * t204;
t344 = t230 * t171;
t437 = mrSges(6,2) * t241;
t108 = -mrSges(6,3) * t344 - t437;
t342 = t232 * t108;
t355 = t171 * t232;
t438 = mrSges(6,1) * t241;
t111 = -mrSges(6,3) * t355 + t438;
t346 = t230 * t111;
t446 = -t346 / 0.2e1 + t342 / 0.2e1;
t226 = Ifges(6,5) * t232;
t374 = Ifges(6,6) * t230;
t436 = t226 - t374;
t441 = t241 / 0.2e1;
t442 = -t171 / 0.2e1;
t445 = t436 * t441 - Ifges(5,4) * t241 + (Ifges(5,2) + Ifges(6,3)) * t442;
t118 = pkin(4) * t241 - t171 * pkin(6);
t215 = Ifges(6,5) * t230 + Ifges(6,6) * t232;
t357 = t241 * t215;
t224 = pkin(2) * t391 + pkin(3);
t292 = t388 * t390;
t198 = pkin(2) * t292 + t231 * t224;
t197 = -pkin(2) * t318 + t224 * t390;
t192 = -pkin(4) - t197;
t406 = t192 / 0.2e1;
t106 = t282 * t241;
t411 = t106 / 0.2e1;
t434 = t198 * t411 + t406 * t448;
t382 = mrSges(6,3) * t232;
t110 = -t171 * t382 + t438;
t112 = -mrSges(6,1) * t171 - t241 * t382;
t403 = -t197 / 0.2e1;
t193 = pkin(6) + t198;
t405 = -t193 / 0.2e1;
t433 = t110 * t405 + t112 * t403;
t383 = mrSges(6,3) * t230;
t107 = -t171 * t383 - t437;
t109 = mrSges(6,2) * t171 - t241 * t383;
t402 = t197 / 0.2e1;
t404 = t193 / 0.2e1;
t432 = t107 * t404 + t109 * t402;
t361 = t232 * mrSges(6,1);
t364 = t230 * mrSges(6,2);
t431 = -t364 / 0.2e1 + t361 / 0.2e1;
t375 = Ifges(6,6) * t171;
t48 = t241 * t279 - t375;
t363 = t230 * t48;
t50 = t241 * t281 - t377;
t429 = Ifges(5,4) * t171 + Ifges(5,1) * t441 - t363 / 0.2e1 + t50 * t393;
t428 = (-t388 * mrSges(4,1) - t391 * mrSges(4,2)) * pkin(2);
t427 = mrSges(6,3) * t297;
t426 = (-t226 / 0.2e1 + t374 / 0.2e1) * t171;
t394 = t230 / 0.2e1;
t425 = t279 * t393 + t281 * t394 + t447;
t263 = t380 + (-Ifges(6,1) + Ifges(6,2)) * t232;
t166 = Ifges(5,6) * t241;
t167 = Ifges(5,5) * t171;
t424 = t166 / 0.2e1 - t167 / 0.2e1 - t357 / 0.4e1 + t448 * t422;
t423 = -pkin(4) / 0.2e1;
t421 = -pkin(6) / 0.2e1;
t420 = m(6) * pkin(3);
t419 = -mrSges(6,1) / 0.2e1;
t418 = mrSges(6,2) / 0.2e1;
t417 = -Ifges(5,5) / 0.2e1;
t416 = Ifges(6,3) / 0.2e1;
t415 = t47 / 0.4e1;
t414 = -t48 / 0.4e1;
t413 = t49 / 0.4e1;
t412 = t50 / 0.4e1;
t410 = t109 / 0.2e1;
t409 = -t171 / 0.4e1;
t408 = t171 / 0.2e1;
t176 = t192 * t282;
t407 = -t176 / 0.2e1;
t338 = t390 * pkin(3);
t223 = -t338 - pkin(4);
t201 = t223 * t282;
t401 = -t201 / 0.2e1;
t203 = (t231 * t391 + t292) * pkin(2);
t400 = t203 / 0.2e1;
t399 = -t204 / 0.2e1;
t398 = t215 / 0.2e1;
t387 = pkin(3) * t231;
t222 = pkin(6) + t387;
t397 = -t222 / 0.2e1;
t396 = t223 / 0.2e1;
t385 = t212 * pkin(3);
t225 = t392 * pkin(2) + pkin(1);
t373 = Ifges(6,3) * t241;
t114 = mrSges(5,1) * t241 + mrSges(5,2) * t171;
t115 = -mrSges(5,1) * t171 + mrSges(5,2) * t241;
t186 = -pkin(3) * t254 + t225;
t337 = t389 * pkin(2);
t191 = -t337 - t385;
t332 = t416 + Ifges(5,2) / 0.2e1;
t233 = t225 * (-t212 * mrSges(4,1) + mrSges(4,2) * t254) + t254 ^ 2 * Ifges(4,4) + (t436 * t442 + t429) * t171 + (-t332 * t171 + (Ifges(6,1) * t355 + t378) * t393 + (-Ifges(6,2) * t344 + t376) * t395 + Ifges(5,1) * t408 - t344 * t379 + t445) * t241;
t238 = -Ifges(4,4) * t212 + (-Ifges(4,1) + Ifges(4,2)) * t254;
t348 = t230 * t109;
t270 = t232 * t112 + t348;
t271 = t230 * t108 + t232 * t111;
t300 = m(6) * t340;
t76 = -pkin(4) * t171 - pkin(6) * t241 + t186;
t78 = t118 - t385;
t77 = -t337 + t78;
t1 = t191 * t115 + t270 * t77 + (m(5) * t191 + t114) * t186 + (t300 * t77 + t271) * t76 + (mrSges(4,2) * t337 + t238) * t212 + t233 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t392) * t392 + (-pkin(1) * mrSges(3,1) + (-m(4) * t225 + mrSges(4,1) * t254) * pkin(2) + (-Ifges(3,2) + Ifges(3,1)) * t392 - Ifges(3,4) * t389) * t389;
t372 = t1 * qJD(1);
t371 = t197 * mrSges(5,2);
t370 = t198 * mrSges(5,1);
t2 = t186 * t114 + t270 * t78 + (t300 * t78 + t271) * t76 + ((-m(5) * t186 - t115) * pkin(3) + t238) * t212 + t233;
t369 = t2 * qJD(1);
t368 = t203 * mrSges(5,1);
t367 = t204 * mrSges(5,2);
t3 = (t107 * t230 + t110 * t232) * t76 + (t186 * mrSges(5,1) + t49 * t393 + t47 * t395 + t445) * t241 - (-t186 * mrSges(5,2) - t426 + (-Ifges(5,1) / 0.2e1 + t332) * t241 - t429) * t171 + (t300 * t76 + t270) * t118;
t359 = t3 * qJD(1);
t341 = t232 * t109;
t345 = t230 * t112;
t16 = (-t171 * t398 + t241 * t447 + t48 * t393 + t50 * t394) * t241 + (t345 - t341) * t76;
t358 = t16 * qJD(1);
t356 = t241 * t231;
t354 = t192 * t448;
t353 = t197 * t171;
t352 = t198 * t241;
t283 = t361 - t364;
t351 = t198 * t283;
t350 = t203 * t283;
t349 = t223 * t448;
t347 = t230 * t110;
t343 = t232 * t107;
t339 = t283 * t387;
t336 = -t387 / 0.2e1;
t103 = t283 * t241;
t335 = t103 * t423;
t334 = t112 * t421;
t333 = -Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.4e1;
t331 = t403 + t422;
t328 = -t384 / 0.2e1;
t319 = t390 * t171;
t317 = t226 * t409;
t316 = t103 * t406;
t309 = t106 * t400;
t308 = t103 * t396;
t307 = t448 * t396;
t306 = -t348 / 0.2e1;
t304 = t231 * t411;
t302 = t408 + t409;
t301 = t404 + t397;
t298 = t340 * t197;
t296 = t340 * t222;
t295 = mrSges(6,3) * t338;
t294 = -t338 / 0.2e1;
t291 = t390 * t410;
t290 = t390 * t395;
t289 = t333 * t230;
t288 = t333 * t228;
t287 = (Ifges(6,1) / 0.4e1 - Ifges(6,2) / 0.2e1) * t232;
t286 = mrSges(6,3) * (-t229 / 0.2e1 - t228 / 0.2e1);
t285 = (t215 / 0.4e1 - Ifges(5,6) / 0.2e1) * t241;
t276 = t282 * t423 + t425;
t247 = mrSges(6,3) * t298 - t351 - t370 - t371;
t19 = m(6) * (t192 * t198 + t193 * t298) + t247;
t237 = t413 - t378 / 0.4e1 + pkin(6) * t111 / 0.2e1 + t263 * t302;
t240 = -t171 * t417 + t285 + t424;
t244 = -t302 * t379 + t415 - t376 / 0.4e1 + t108 * t421;
t4 = (t244 + t432) * t232 + (t237 + t433) * t230 + t240 + t434;
t275 = t4 * qJD(1) + t19 * qJD(2);
t249 = (-0.3e1 / 0.2e1 * t380 + t287) * t241 + t412;
t243 = t112 * t405 + t249;
t246 = (t171 / 0.4e1 + t408) * Ifges(6,6) + t241 * t289 + t414;
t264 = -t373 / 0.2e1 + t317;
t267 = t193 * t286;
t10 = t316 + t241 * t267 + (t109 * t405 + t418 * t77 + t246) * t230 + (t419 * t77 + t243 + t327) * t232 + t264;
t252 = -Ifges(6,4) * t229 + t263 * t230;
t88 = -t176 + t252;
t274 = t10 * qJD(1) - t88 * qJD(2);
t245 = -t352 / 0.2e1 + t241 * t400 - t171 * t399 - t353 / 0.2e1;
t15 = t309 + (-t223 / 0.2e1 + t406) * t448 + (t108 * t301 + t204 * t410) * t232 + (-t111 * t301 + t112 * t399) * t230 + ((t319 / 0.2e1 + t356 / 0.2e1) * pkin(3) + t245) * mrSges(5,3);
t20 = (-mrSges(5,1) - t283) * t203 + t428 + (mrSges(6,3) * t340 - mrSges(5,2)) * t204 + m(5) * (-t197 * t203 + t198 * t204) + m(6) * (t192 * t203 + t193 * t297);
t273 = t15 * qJD(1) + t20 * qJD(2);
t272 = t342 - t346;
t269 = pkin(6) * t286;
t268 = t294 + t422;
t266 = t222 * t286;
t265 = t328 - t365 / 0.2e1;
t260 = t340 * t390;
t234 = m(6) * (t223 * t198 + t197 * t296 + (t192 * t231 + t193 * t260) * pkin(3)) / 0.2e1 - t371 / 0.2e1 - t370 / 0.2e1 - t351 / 0.2e1 + mrSges(5,1) * t336 - t339 / 0.2e1 + mrSges(5,2) * t294 + t340 * (mrSges(6,3) * t402 + t295 / 0.2e1);
t236 = -m(6) * (-pkin(4) * t203 + pkin(6) * t297) / 0.2e1 + t368 / 0.2e1 + t350 / 0.2e1 + t367 / 0.2e1 - t427 / 0.2e1;
t18 = t234 + t236;
t6 = t307 + pkin(3) * t304 + (pkin(3) * t291 + t222 * t107 / 0.2e1 + t244) * t232 + (t110 * t397 + t112 * t294 + t237) * t230 + t240;
t239 = -mrSges(5,1) * t387 - mrSges(5,2) * t338 + t295 * t340 - t339;
t89 = (t222 * t260 + t223 * t231) * t420 + t239;
t258 = t6 * qJD(1) + t18 * qJD(2) + t89 * qJD(3);
t113 = -t201 + t252;
t242 = t112 * t397 + t249;
t12 = t308 + t241 * t266 + (t109 * t397 + t418 * t78 + t246) * t230 + (t419 * t78 + t242 + t327) * t232 + t264;
t31 = t401 + t407 + (mrSges(6,2) * t399 - t379) * t232 + (mrSges(6,1) * t399 + t263) * t230;
t257 = t12 * qJD(1) - t31 * qJD(2) - t113 * qJD(3);
t256 = t357 / 0.2e1 + t167 + Ifges(4,6) * t212 + Ifges(4,5) * t254 - t166 + t447 * t171 + t452;
t255 = t436 * t409 - t363 / 0.4e1;
t119 = (pkin(4) * mrSges(6,2) - t379) * t232 + (pkin(4) * mrSges(6,1) + t263) * t230;
t14 = t317 + t335 + (-Ifges(6,3) / 0.2e1 + t269) * t241 + (t118 * t419 + t241 * t287 + t327 + t334 + t412) * t232 + (0.3e1 / 0.4e1 * t375 + t414 + t109 * t421 + t118 * t418 + (-0.3e1 / 0.2e1 * t379 + t289) * t241) * t230;
t33 = t407 + (mrSges(6,2) * t331 - t379) * t232 + (mrSges(6,1) * t331 + t263) * t230;
t55 = t401 + (mrSges(6,2) * t268 - t379) * t232 + (mrSges(6,1) * t268 + t263) * t230;
t253 = t14 * qJD(1) - t33 * qJD(2) - t55 * qJD(3) - t119 * qJD(4);
t248 = -Ifges(6,6) * t344 / 0.2e1 + Ifges(6,5) * t355 / 0.2e1 + t373 / 0.2e1 + t255;
t235 = t285 - t424 + t446 * pkin(6) + (-t417 + t447) * t171 + t451;
t196 = t201 / 0.2e1;
t174 = t176 / 0.2e1;
t56 = t196 + (mrSges(6,1) * t290 + t328 * t390) * pkin(3) + t276;
t34 = t197 * t265 + t174 + t276;
t32 = t204 * t265 + t174 + t196 + t425;
t17 = t234 - t236;
t13 = pkin(6) * t306 + t335 + (t334 + t249) * t232 + t255 + (t416 + t269 + t288) * t241 - t426 + t431 * t118;
t11 = t222 * t306 + t308 + (t266 + t288) * t241 + t242 * t232 + t248 + t431 * t78;
t9 = t193 * t306 + t316 + (t267 + t288) * t241 + t243 * t232 + t248 + t431 * t77;
t8 = t349 / 0.2e1 + t309 + t354 / 0.2e1 + (-t345 / 0.2e1 + t341 / 0.2e1) * t204 + t256 + (t171 * t294 + t241 * t336 + t245) * mrSges(5,3) + (t222 + t193) * t446;
t7 = (t343 / 0.2e1 - t347 / 0.2e1) * t222 + t235 + t307 + (t112 * t290 + t232 * t291 + t304) * pkin(3) + t451;
t5 = t235 + (t413 + t433) * t230 + (t415 + t432) * t232 + t434;
t21 = [qJD(2) * t1 + qJD(3) * t2 + qJD(4) * t3 - qJD(5) * t16, t372 + (-Ifges(3,5) * t392 + Ifges(3,6) * t389 + t354 + t272 * t193 + (-t352 - t353) * mrSges(5,3) + t256 + (t212 * t388 - t254 * t391) * pkin(2) * mrSges(4,3)) * qJD(2) + t8 * qJD(3) + t5 * qJD(4) + t9 * qJD(5), t369 + t8 * qJD(2) + (t349 + t272 * t222 + (-t319 - t356) * pkin(3) * mrSges(5,3) + t256) * qJD(3) + t7 * qJD(4) + t11 * qJD(5), t5 * qJD(2) + t7 * qJD(3) + t13 * qJD(5) + t359 + (-pkin(4) * t448 - (-Ifges(5,5) - t447) * t171 + (t398 - Ifges(5,6)) * t241 + (-t347 + t343) * pkin(6) + t452) * qJD(4), -t358 + t9 * qJD(2) + t11 * qJD(3) + t13 * qJD(4) + (-t282 * t76 - t357) * qJD(5); qJD(3) * t15 + qJD(4) * t4 + qJD(5) * t10 - t372, qJD(3) * t20 + qJD(4) * t19 - qJD(5) * t88 (m(5) * (-t203 * t390 + t204 * t231) * pkin(3) + m(6) * (t223 * t203 + t204 * t296) - t350 - t367 - t368 + t427 + t428) * qJD(3) + t17 * qJD(4) + t32 * qJD(5) + t273, t17 * qJD(3) + (m(6) * (-pkin(4) * t198 + pkin(6) * t298) + t247) * qJD(4) + t34 * qJD(5) + t275, t32 * qJD(3) + t34 * qJD(4) + (-t193 * t283 + t436) * qJD(5) + t274; -qJD(2) * t15 + qJD(4) * t6 + qJD(5) * t12 - t369, qJD(4) * t18 - qJD(5) * t31 - t273, qJD(4) * t89 - qJD(5) * t113 ((-pkin(4) * t231 + pkin(6) * t260) * t420 + t239) * qJD(4) + t56 * qJD(5) + t258, t56 * qJD(4) + (-t222 * t283 + t436) * qJD(5) + t257; -qJD(2) * t4 - qJD(3) * t6 + qJD(5) * t14 - t359, -qJD(3) * t18 - qJD(5) * t33 - t275, -qJD(5) * t55 - t258, -t119 * qJD(5) (-pkin(6) * t283 + t436) * qJD(5) + t253; -qJD(2) * t10 - qJD(3) * t12 - qJD(4) * t14 + t358, qJD(3) * t31 + qJD(4) * t33 - t274, qJD(4) * t55 - t257, -t253, 0;];
Cq  = t21;
