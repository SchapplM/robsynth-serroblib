% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:17
% EndTime: 2019-03-09 02:28:24
% DurationCPUTime: 3.78s
% Computational Cost: add. (9096->376), mult. (16752->496), div. (0->0), fcn. (16567->6), ass. (0->247)
t234 = sin(qJ(5));
t235 = sin(qJ(4));
t237 = cos(qJ(4));
t399 = cos(qJ(5));
t206 = t234 * t237 + t235 * t399;
t233 = sin(qJ(6));
t227 = t233 ^ 2;
t236 = cos(qJ(6));
t228 = t236 ^ 2;
t329 = t227 + t228;
t300 = t329 * t206;
t444 = mrSges(7,3) * t300;
t205 = t234 * t235 - t237 * t399;
t385 = Ifges(7,4) * t236;
t290 = Ifges(7,1) * t233 + t385;
t270 = t236 * t290;
t386 = Ifges(7,4) * t233;
t288 = Ifges(7,2) * t236 + t386;
t271 = t233 * t288;
t287 = Ifges(7,5) * t233 + Ifges(7,6) * t236;
t272 = t205 * t287;
t400 = t236 / 0.2e1;
t402 = t233 / 0.2e1;
t405 = t206 / 0.2e1;
t289 = -Ifges(7,2) * t233 + t385;
t97 = -Ifges(7,6) * t205 - t206 * t289;
t291 = Ifges(7,1) * t236 - t386;
t99 = -Ifges(7,5) * t205 - t206 * t291;
t268 = t271 * t405 + t97 * t400 + t99 * t402 - t272 / 0.2e1 + Ifges(6,6) * t205 + (-t270 / 0.2e1 - Ifges(6,5)) * t206;
t370 = t236 * mrSges(7,1);
t374 = t233 * mrSges(7,2);
t212 = -t370 + t374;
t229 = qJ(2) - pkin(7);
t391 = -pkin(8) + t229;
t211 = t391 * t235;
t303 = t391 * t237;
t257 = t211 * t399 + t234 * t303;
t430 = t257 * t212;
t432 = t257 * mrSges(6,1);
t151 = t211 * t234 - t399 * t303;
t440 = t151 * mrSges(6,2);
t443 = t268 + t430 - t432 + t440;
t442 = -t430 / 0.2e1 + t432 / 0.2e1 - t440 / 0.2e1;
t369 = t236 * mrSges(7,2);
t375 = t233 * mrSges(7,1);
t213 = t369 + t375;
t404 = -t213 / 0.2e1;
t230 = pkin(1) + qJ(3);
t217 = t235 * pkin(4) + t230;
t441 = m(6) * t217;
t439 = t151 * t233;
t438 = t151 * t234;
t437 = t151 * t236;
t358 = t151 * t257;
t136 = t206 * t213;
t347 = t206 * t233;
t145 = mrSges(7,2) * t205 + mrSges(7,3) * t347;
t346 = t206 * t236;
t147 = -mrSges(7,1) * t205 + mrSges(7,3) * t346;
t155 = -t205 * mrSges(6,1) - t206 * mrSges(6,2);
t100 = Ifges(7,5) * t206 - t205 * t291;
t339 = t236 * t100;
t98 = Ifges(7,6) * t206 - t205 * t289;
t372 = t233 * t98;
t225 = Ifges(7,5) * t236;
t384 = Ifges(7,6) * t233;
t279 = t225 / 0.2e1 - t384 / 0.2e1;
t429 = t279 * t206;
t394 = pkin(5) * t206;
t134 = pkin(9) * t205 + t217 + t394;
t56 = t134 * t236 - t233 * t257;
t57 = t134 * t233 + t236 * t257;
t434 = (Ifges(6,4) * t206 - t339 / 0.2e1 + t372 / 0.2e1 - t429) * t206 - t151 * t136 + t57 * t145 + t56 * t147 + t217 * t155;
t433 = pkin(5) * t257;
t324 = t399 * pkin(4);
t221 = -t324 - pkin(5);
t431 = t221 * t257;
t277 = t369 / 0.2e1 + t375 / 0.2e1;
t79 = (t404 + t277) * t205;
t332 = t79 * qJD(3);
t276 = t386 + (-Ifges(7,1) + Ifges(7,2)) * t236;
t189 = t221 * t213;
t393 = pkin(5) * t213;
t293 = -t189 / 0.2e1 + t393 / 0.2e1;
t409 = -mrSges(7,2) / 0.2e1;
t296 = t399 * t409;
t410 = -mrSges(7,1) / 0.2e1;
t45 = (pkin(4) * t296 - t385) * t236 + (t324 * t410 + t276) * t233 + t293;
t254 = -Ifges(7,4) * t228 + t233 * t276;
t72 = t254 + t393;
t305 = t228 / 0.2e1 + t227 / 0.2e1;
t387 = mrSges(7,3) * t205;
t349 = t205 * t236;
t148 = t206 * mrSges(7,1) + mrSges(7,3) * t349;
t336 = t236 * t148;
t351 = t205 * t233;
t146 = -t206 * mrSges(7,2) + mrSges(7,3) * t351;
t342 = t233 * t146;
t425 = -t336 / 0.2e1 - t342 / 0.2e1;
t252 = t305 * t387 + t425;
t137 = t288 * t205;
t138 = t290 * t205;
t373 = t233 * mrSges(7,3);
t316 = t373 / 0.2e1;
t317 = -t373 / 0.2e1;
t214 = t225 - t384;
t421 = -t206 * t214 / 0.4e1 + t151 * t404;
t253 = t138 * t402 + t137 * t400 + t289 * t351 / 0.4e1 - t291 * t349 / 0.4e1 - t372 / 0.4e1 + t339 / 0.4e1 - t421 + (t316 + t317) * t57;
t135 = t205 * t212;
t407 = -t135 / 0.2e1;
t242 = pkin(5) * t407 + pkin(9) * t252 + t253;
t395 = pkin(5) * t205;
t156 = pkin(9) * t206 - t395;
t70 = t156 * t236 + t439;
t71 = t156 * t233 - t437;
t280 = t70 * t410 + t71 * mrSges(7,2) / 0.2e1;
t383 = Ifges(7,3) * t205;
t9 = t242 + t429 + t383 / 0.2e1 + t280;
t428 = -t9 * qJD(1) + t45 * qJD(4) + t72 * qJD(5) + t332;
t412 = m(6) * pkin(4);
t327 = -t412 / 0.2e1;
t413 = m(7) / 0.2e1;
t392 = t237 * pkin(4);
t139 = t156 + t392;
t61 = t139 * t236 + t439;
t62 = t139 * t233 - t437;
t427 = (-t233 * t62 - t236 * t61) * t413 + t237 * t327;
t315 = Ifges(7,6) * t347 / 0.2e1 - Ifges(7,5) * t346 / 0.2e1 - t383 / 0.2e1;
t256 = t61 * mrSges(7,1) / 0.2e1 + t62 * t409 + t315;
t292 = t305 * mrSges(7,3);
t322 = Ifges(7,1) / 0.4e1 - Ifges(7,2) / 0.4e1;
t344 = t221 * t135;
t397 = pkin(4) * t234;
t220 = pkin(9) + t397;
t352 = t205 * t220;
t403 = t220 / 0.2e1;
t5 = -t344 / 0.2e1 - t292 * t352 + (-t100 / 0.4e1 - t137 / 0.4e1 + t148 * t403 + t322 * t349) * t236 + (-t138 / 0.4e1 + t98 / 0.4e1 + t146 * t403 + (-t233 * t322 - t385) * t205) * t233 + t256 + t421;
t60 = -t189 + t254;
t426 = -t5 * qJD(1) - t60 * qJD(4) - t332;
t338 = t236 * t145;
t341 = t233 * t147;
t273 = -t338 / 0.2e1 + t341 / 0.2e1;
t424 = -t235 * mrSges(5,1) - t237 * mrSges(5,2);
t417 = t206 ^ 2;
t418 = t205 ^ 2;
t423 = t417 + t418;
t172 = t206 * t212;
t422 = -t329 * t387 + t172;
t401 = -t236 / 0.2e1;
t420 = -t213 * t257 + (-Ifges(7,3) - Ifges(6,2) + Ifges(6,1)) * t206 + (-Ifges(6,4) + t279) * t205 + t97 * t402 + t99 * t401;
t337 = t236 * t146;
t340 = t233 * t148;
t274 = t340 / 0.2e1 - t337 / 0.2e1;
t286 = t233 * t56 - t236 * t57;
t414 = -m(7) / 0.2e1;
t419 = ((t257 + t286) * t414 + t136 - t274) * t205;
t411 = m(7) * pkin(4);
t408 = mrSges(7,3) / 0.2e1;
t406 = t205 / 0.2e1;
t301 = t329 * t205;
t398 = m(7) * (-pkin(9) * t301 - t394);
t396 = pkin(5) * t136;
t366 = t62 * t236;
t367 = t61 * t233;
t285 = t366 - t367;
t11 = ((t151 + t285) * t414 + t273) * t206 + t419;
t364 = t71 * t236;
t284 = -t70 * t233 + t364;
t13 = ((t151 + t284) * t414 + t273) * t206 + t419;
t390 = -t11 * qJD(4) - t13 * qJD(5);
t40 = m(7) * (0.1e1 - t329) * t206 * t205;
t334 = t40 * qJD(3);
t80 = t205 * t277 + t213 * t406;
t389 = t80 * qJD(6) + t334;
t388 = -t79 * qJD(6) - t334;
t7 = t151 * t135 + t56 * t146 - t57 * t148 + (-t286 * mrSges(7,3) + t138 * t401 + t287 * t405 + t98 * t400 + (t100 + t137) * t402) * t205;
t365 = t7 * qJD(1);
t304 = t235 ^ 2 + t237 ^ 2;
t294 = m(5) * t304;
t14 = t304 * mrSges(5,3) + t418 * t213 - t229 * t294 - mrSges(4,2) - mrSges(3,3) + (-m(7) - m(6)) * t151 * t205 + (-m(6) * t257 + m(7) * t286 - t337 + t340) * t206 + (-m(4) - m(3)) * qJ(2) + t423 * mrSges(6,3);
t363 = qJD(1) * t14;
t197 = t206 * mrSges(6,1);
t302 = t205 * mrSges(6,2) - t197;
t259 = t302 + t424;
t25 = t342 + t336 + mrSges(4,3) + m(7) * (t233 * t57 + t236 * t56) + t441 + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t230 - t259;
t362 = qJD(1) * t25;
t361 = t11 * qJD(1);
t360 = t13 * qJD(1);
t224 = t237 * mrSges(5,1);
t246 = (t221 * t205 + t220 * t300) * t414 + (-t205 * t399 + t206 * t234) * t327;
t275 = -t233 * t145 / 0.2e1 + t147 * t401;
t15 = t235 * mrSges(5,2) - t155 - t224 + t246 + t275 + t407 - t444 / 0.2e1 + t427;
t359 = t15 * qJD(1);
t248 = t135 * t406 + t206 * t252;
t278 = -t374 / 0.2e1 + t370 / 0.2e1;
t19 = t248 - t278;
t353 = t19 * qJD(1);
t350 = t205 * t234;
t258 = m(7) * (pkin(9) * t300 - t395);
t265 = (-t233 * t71 - t236 * t70) * t413;
t21 = (mrSges(6,1) - t212 / 0.2e1) * t205 + (mrSges(6,2) - t292) * t206 + t265 - t258 / 0.2e1 + t275;
t345 = t21 * qJD(1);
t343 = t221 * t136;
t263 = t277 * t206;
t26 = -t263 - t274;
t335 = t26 * qJD(1);
t247 = t294 / 0.2e1 + t423 * m(6) / 0.2e1 + (t329 * t417 + t418) * t413;
t267 = -m(5) / 0.2e1 - m(6) / 0.2e1 + t329 * t414;
t43 = -m(4) - t247 + t267;
t333 = t43 * qJD(1);
t328 = qJD(4) + qJD(5);
t326 = pkin(4) * t350;
t325 = mrSges(7,3) * t366;
t323 = t398 / 0.2e1;
t321 = t220 * t341;
t320 = t220 * t338;
t314 = t233 * t399;
t313 = t236 * t399;
t312 = t399 * t206;
t295 = -t314 / 0.2e1;
t1 = t62 * t146 + t61 * t148 + t230 * t224 + (-Ifges(5,4) * t237 + pkin(4) * t197) * t237 + t392 * t441 + m(7) * (t56 * t61 + t57 * t62 + t358) + (-t230 * mrSges(5,2) + Ifges(5,4) * t235 + (-Ifges(5,1) + Ifges(5,2)) * t237) * t235 + (-mrSges(6,2) * t392 + t420) * t205 + t434;
t283 = t1 * qJD(1) - t11 * qJD(3);
t3 = m(7) * (t56 * t70 + t57 * t71 + t358) + t71 * t146 + t70 * t148 + t420 * t205 + t434;
t282 = t3 * qJD(1) - t13 * qJD(3);
t269 = t302 + t422;
t264 = t329 * t399;
t260 = t221 * t206 - t220 * t301;
t243 = m(7) * ((t206 * t264 + t350) * pkin(4) + t260);
t29 = t323 - t243 / 0.2e1;
t240 = (t284 * t220 + t431) * t413 - t343 / 0.2e1 - t321 / 0.2e1 + t320 / 0.2e1 + t70 * t317 + t364 * t408 + t326 * t404 + ((t313 * t57 - t314 * t56 + t438) * t413 + t148 * t295 + t146 * t313 / 0.2e1) * pkin(4) - t442;
t241 = -t414 * t433 - t396 / 0.2e1 + t61 * t316 - t325 / 0.2e1 + (t285 * t414 + t273) * pkin(9) + t442;
t4 = t240 + t241;
t244 = (-mrSges(6,1) + t212) * t397 + (t329 * mrSges(7,3) - mrSges(6,2)) * t324;
t48 = (t220 * t264 + t221 * t234) * t411 + t244;
t262 = t4 * qJD(1) - t29 * qJD(3) + t48 * qJD(4);
t249 = t135 / 0.2e1 + t275 + t444 / 0.2e1;
t46 = -t271 / 0.2e1 + t289 * t400 + t270 / 0.2e1 + t291 * t402 + (mrSges(7,1) * t295 + t236 * t296) * pkin(4) - t293;
t44 = t247 + t267;
t27 = -t263 + t274;
t24 = t243 / 0.2e1 + t323 + t269;
t22 = t265 + t258 / 0.2e1 + t249;
t20 = t248 + t278;
t16 = -t246 + t249 + t427;
t8 = t242 - t280 + t315;
t6 = t253 + t344 / 0.2e1 + t256 + t329 * t352 * t408 + t425 * t220;
t2 = -t241 + t240 + t268;
t10 = [-qJD(2) * t14 + qJD(3) * t25 + qJD(4) * t1 + qJD(5) * t3 + qJD(6) * t7, qJD(3) * t44 + qJD(4) * t16 + qJD(5) * t22 + qJD(6) * t27 - t363, qJD(2) * t44 + qJD(6) * t20 + t362 + t390, t16 * qJD(2) + ((-t257 * t399 - t438) * t412 + t320 + t325 - mrSges(7,3) * t367 - t321 + m(7) * (t220 * t285 + t431) - t343 - Ifges(5,5) * t235 - Ifges(5,6) * t237 + t424 * t229 + (pkin(4) * t312 + t326) * mrSges(6,3) + t443) * qJD(4) + t2 * qJD(5) + t6 * qJD(6) + t283, t22 * qJD(2) + t2 * qJD(4) + (-m(7) * t433 + t284 * mrSges(7,3) + t396 + (m(7) * t284 + t338 - t341) * pkin(9) + t443) * qJD(5) + t8 * qJD(6) + t282, t365 + t27 * qJD(2) + t20 * qJD(3) + t6 * qJD(4) + t8 * qJD(5) + (-mrSges(7,1) * t57 - mrSges(7,2) * t56 + t272) * qJD(6); qJD(3) * t43 + qJD(4) * t15 + qJD(5) * t21 - qJD(6) * t26 + t363, 0, t333, t359, t345, qJD(6) * t213 - t335; -qJD(2) * t43 + qJD(6) * t19 - t362 + t390, -t333, t328 * t40, -t361 + (m(7) * t260 + (-t312 - t350) * t412 + t259 + t422) * qJD(4) + t24 * qJD(5) + t389, -t360 + t24 * qJD(4) + (t269 + t398) * qJD(5) + t389, qJD(6) * t172 + t328 * t80 + t353; -qJD(2) * t15 + qJD(5) * t4 - qJD(6) * t5 - t283, -t359, -qJD(5) * t29 + t361 + t388, qJD(5) * t48 - qJD(6) * t60 ((-pkin(5) * t234 + pkin(9) * t264) * t411 + t244) * qJD(5) + t46 * qJD(6) + t262, t46 * qJD(5) + (t212 * t220 + t214) * qJD(6) + t426; -qJD(2) * t21 - qJD(4) * t4 + qJD(6) * t9 - t282, -t345, qJD(4) * t29 + t360 + t388, -qJD(6) * t45 - t262, -t72 * qJD(6) (pkin(9) * t212 + t214) * qJD(6) - t428; qJD(2) * t26 - qJD(3) * t19 + qJD(4) * t5 - qJD(5) * t9 - t365, t335, t328 * t79 - t353, qJD(5) * t45 - t426, t428, 0;];
Cq  = t10;
