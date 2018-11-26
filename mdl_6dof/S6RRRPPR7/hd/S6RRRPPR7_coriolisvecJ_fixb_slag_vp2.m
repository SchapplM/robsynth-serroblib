% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:37:42
% EndTime: 2018-11-23 17:37:58
% DurationCPUTime: 16.20s
% Computational Cost: add. (10086->762), mult. (24220->1012), div. (0->0), fcn. (16190->8), ass. (0->351)
t305 = sin(qJ(6));
t308 = cos(qJ(6));
t306 = sin(qJ(3));
t307 = sin(qJ(2));
t379 = qJD(1) * t307;
t363 = t306 * t379;
t309 = cos(qJ(3));
t370 = t309 * qJD(2);
t258 = t363 - t370;
t361 = t309 * t379;
t259 = qJD(2) * t306 + t361;
t302 = sin(pkin(10));
t303 = cos(pkin(10));
t354 = t258 * t302 + t303 * t259;
t355 = t303 * t258 - t259 * t302;
t495 = t305 * t355 + t308 * t354;
t496 = -t305 * t354 + t308 * t355;
t81 = Ifges(7,4) * t496;
t515 = Ifges(7,2) * t495 - t81;
t417 = pkin(8) - qJ(5);
t272 = t417 * t309;
t423 = pkin(7) * t306;
t364 = -pkin(3) - t423;
t353 = -pkin(4) + t364;
t310 = cos(qJ(2));
t389 = t309 * t310;
t346 = pkin(2) * t307 - pkin(8) * t310;
t263 = t346 * qJD(1);
t392 = t263 * t309;
t514 = t392 - (-qJ(5) * t389 + t307 * t353) * qJD(1) + qJD(3) * t272 - qJD(5) * t306;
t235 = t306 * t263;
t290 = qJ(4) * t379;
t371 = qJD(5) * t309;
t375 = qJD(3) * t306;
t390 = t307 * t309;
t391 = t306 * t310;
t513 = t235 + t290 + (-pkin(7) * t390 + qJ(5) * t391) * qJD(1) + t375 * t417 + t371;
t378 = qJD(1) * t310;
t286 = -qJD(3) + t378;
t279 = qJD(6) + t286;
t432 = -t279 / 0.2e1;
t273 = -qJD(2) * pkin(2) + pkin(7) * t379;
t152 = t258 * pkin(3) - t259 * qJ(4) + t273;
t114 = -pkin(4) * t258 + qJD(5) - t152;
t66 = -pkin(5) * t355 + t114;
t268 = -pkin(2) * t310 - pkin(8) * t307 - pkin(1);
t237 = t268 * qJD(1);
t295 = pkin(7) * t378;
t274 = qJD(2) * pkin(8) + t295;
t177 = t306 * t237 + t309 * t274;
t131 = qJ(5) * t258 + t177;
t276 = t286 * qJ(4);
t112 = t131 - t276;
t176 = t309 * t237 - t306 * t274;
t130 = qJ(5) * t259 + t176;
t448 = pkin(3) + pkin(4);
t97 = t286 * t448 + qJD(4) - t130;
t44 = -t112 * t302 + t303 * t97;
t475 = pkin(9) * t354;
t33 = pkin(5) * t286 + t44 - t475;
t45 = t303 * t112 + t302 * t97;
t494 = pkin(9) * t355;
t37 = t45 + t494;
t8 = -t305 * t37 + t308 * t33;
t9 = t305 * t33 + t308 * t37;
t512 = (t495 * t9 + t496 * t8) * mrSges(7,3) + (Ifges(7,5) * t496 - Ifges(7,6) * t495) * t432 - t66 * (mrSges(7,1) * t495 + mrSges(7,2) * t496);
t323 = t302 * t309 - t303 * t306;
t318 = t310 * t323;
t209 = qJD(1) * t318;
t229 = t323 * qJD(3);
t511 = t209 - t229;
t322 = t302 * t306 + t303 * t309;
t317 = t322 * t310;
t210 = qJD(1) * t317;
t230 = t322 * qJD(3);
t385 = t210 - t230;
t424 = Ifges(7,4) * t495;
t508 = Ifges(7,1) * t496 - t424;
t36 = Ifges(7,1) * t495 + Ifges(7,5) * t279 + t81;
t506 = t36 / 0.2e1;
t450 = -t495 / 0.2e1;
t486 = t302 * t513 + t303 * t514;
t485 = t302 * t514 - t303 * t513;
t501 = qJD(2) / 0.2e1;
t493 = Ifges(4,1) + Ifges(5,1);
t489 = Ifges(5,4) + Ifges(4,5);
t500 = Ifges(4,6) - Ifges(5,6);
t499 = pkin(5) * t379 + t385 * pkin(9) + t486;
t498 = t511 * pkin(9) + t485;
t373 = qJD(3) * t309;
t497 = -t306 * qJD(4) - t295 + (t309 * t378 - t373) * qJ(4);
t358 = Ifges(3,5) * t501;
t492 = Ifges(6,4) * t355;
t362 = t306 * t378;
t365 = t448 * t306;
t482 = -qJD(3) * t365 + t362 * t448 - t497;
t76 = Ifges(6,1) * t354 + t286 * Ifges(6,5) + t492;
t491 = t76 / 0.2e1;
t490 = -qJD(2) / 0.2e1;
t271 = t417 * t306;
t187 = t303 * t271 - t272 * t302;
t144 = pkin(9) * t323 + t187;
t188 = t302 * t271 + t303 * t272;
t145 = -pkin(9) * t322 + t188;
t63 = t144 * t305 + t145 * t308;
t488 = -qJD(6) * t63 - t305 * t498 + t308 * t499;
t62 = t144 * t308 - t145 * t305;
t487 = qJD(6) * t62 + t305 * t499 + t308 * t498;
t374 = qJD(3) * t307;
t359 = t306 * t374;
t360 = t310 * t370;
t369 = qJD(2) * qJD(3);
t199 = t309 * t369 + (-t359 + t360) * qJD(1);
t376 = qJD(2) * t310;
t316 = t306 * t376 + t307 * t373;
t200 = qJD(1) * t316 + t306 * t369;
t377 = qJD(2) * t307;
t356 = qJD(1) * t377;
t484 = t489 * t356 + (-Ifges(4,4) + Ifges(5,5)) * t200 + t493 * t199;
t483 = -pkin(5) * t511 + t482;
t481 = (-t362 + t375) * pkin(3) + t497;
t480 = t306 * t489 + t309 * t500;
t408 = Ifges(5,5) * t309;
t412 = Ifges(4,4) * t309;
t479 = t306 * t493 - t408 + t412;
t111 = t199 * t303 + t200 * t302;
t264 = t346 * qJD(2);
t239 = qJD(1) * t264;
t352 = pkin(7) * t356;
t99 = t237 * t373 + t306 * t239 - t274 * t375 - t309 * t352;
t78 = qJ(4) * t356 - t286 * qJD(4) + t99;
t51 = qJ(5) * t200 + qJD(5) * t258 + t78;
t321 = t353 * qJD(2);
t345 = t237 * t375 - t239 * t309 + t274 * t373;
t54 = -qJ(5) * t199 - qJD(5) * t259 + t321 * t379 + t345;
t14 = -t302 * t51 + t303 * t54;
t10 = -pkin(5) * t356 - pkin(9) * t111 + t14;
t110 = -t199 * t302 + t200 * t303;
t15 = t302 * t54 + t303 * t51;
t11 = pkin(9) * t110 + t15;
t1 = qJD(6) * t8 + t10 * t305 + t11 * t308;
t2 = -qJD(6) * t9 + t10 * t308 - t11 * t305;
t29 = qJD(6) * t496 + t110 * t305 + t111 * t308;
t30 = -qJD(6) * t495 + t110 * t308 - t111 * t305;
t478 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t29 + Ifges(7,6) * t30;
t293 = Ifges(3,4) * t378;
t245 = Ifges(5,5) * t259;
t153 = -t286 * Ifges(5,6) + t258 * Ifges(5,3) + t245;
t402 = t259 * Ifges(4,4);
t156 = -t258 * Ifges(4,2) - t286 * Ifges(4,6) + t402;
t325 = t176 * t309 + t177 * t306;
t461 = qJD(4) - t176;
t146 = pkin(3) * t286 + t461;
t148 = -t276 + t177;
t326 = t146 * t309 - t148 * t306;
t331 = Ifges(5,3) * t306 + t408;
t335 = -Ifges(4,2) * t306 + t412;
t340 = mrSges(5,1) * t306 - mrSges(5,3) * t309;
t342 = mrSges(4,1) * t306 + mrSges(4,2) * t309;
t406 = Ifges(5,6) * t306;
t407 = Ifges(4,6) * t306;
t410 = Ifges(4,5) * t309;
t411 = Ifges(5,4) * t309;
t425 = t309 / 0.2e1;
t427 = t306 / 0.2e1;
t428 = -t306 / 0.2e1;
t430 = -t286 / 0.2e1;
t433 = t259 / 0.2e1;
t435 = t258 / 0.2e1;
t436 = -t258 / 0.2e1;
t409 = Ifges(5,5) * t306;
t413 = Ifges(4,4) * t306;
t462 = t309 * t493 + t409 - t413;
t246 = Ifges(4,4) * t258;
t403 = t258 * Ifges(5,5);
t466 = t259 * t493 - t286 * t489 - t246 + t403;
t312 = t326 * mrSges(5,2) - t325 * mrSges(4,3) + t152 * t340 + t153 * t427 + t156 * t428 + t273 * t342 + t331 * t435 + t335 * t436 + t462 * t433 + (-t407 + t410 + t406 + t411) * t430 + t466 * t425;
t477 = t312 + Ifges(3,1) * t379 / 0.2e1 + t293 / 0.2e1 + t358;
t456 = t29 / 0.2e1;
t455 = t30 / 0.2e1;
t445 = t110 / 0.2e1;
t444 = t111 / 0.2e1;
t441 = -t354 / 0.2e1;
t476 = -t356 / 0.2e1;
t357 = Ifges(3,6) * t490;
t473 = Ifges(6,4) * t354;
t265 = -qJ(4) * t302 - t303 * t448;
t260 = -pkin(5) + t265;
t266 = t303 * qJ(4) - t302 * t448;
t175 = t260 * t305 + t266 * t308;
t253 = t302 * t308 + t303 * t305;
t59 = -t130 * t302 + t303 * t131;
t42 = t59 + t494;
t60 = t303 * t130 + t302 * t131;
t43 = t60 + t475;
t472 = -qJD(4) * t253 - qJD(6) * t175 + t305 * t43 - t308 * t42;
t174 = t260 * t308 - t266 * t305;
t324 = t302 * t305 - t303 * t308;
t471 = -qJD(4) * t324 + qJD(6) * t174 - t305 * t42 - t308 * t43;
t465 = t279 * t324;
t464 = t279 * t253;
t463 = qJ(4) * t377 - qJD(4) * t310;
t100 = t306 * t352 - t345;
t347 = t307 * t364;
t320 = qJD(2) * t347;
t92 = qJD(1) * t320 + t345;
t460 = -t100 * mrSges(4,1) + t92 * mrSges(5,1) + t14 * mrSges(6,1) + t99 * mrSges(4,2) - t15 * mrSges(6,2) - t78 * mrSges(5,3) + Ifges(6,5) * t111 + Ifges(6,6) * t110 + t478;
t366 = -Ifges(5,6) / 0.2e1 + Ifges(4,6) / 0.2e1;
t367 = Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1;
t368 = Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t414 = Ifges(3,4) * t307;
t459 = -t366 * t258 + t368 * t259 - (Ifges(6,3) / 0.2e1 + t367) * t286 + t148 * mrSges(5,3) + t176 * mrSges(4,1) + t45 * mrSges(6,2) + t9 * mrSges(7,2) + Ifges(4,6) * t436 + Ifges(5,6) * t435 + t357 - (t310 * Ifges(3,2) + t414) * qJD(1) / 0.2e1 - t279 * Ifges(7,3) - t495 * Ifges(7,5) - t496 * Ifges(7,6) - t354 * Ifges(6,5) - t355 * Ifges(6,6) - t146 * mrSges(5,1) - t177 * mrSges(4,2) - t44 * mrSges(6,1) - t8 * mrSges(7,1) + t489 * t433 + (Ifges(4,3) + Ifges(5,2) + Ifges(6,3)) * t430;
t458 = Ifges(7,4) * t456 + Ifges(7,2) * t455 + Ifges(7,6) * t476;
t457 = Ifges(7,1) * t456 + Ifges(7,4) * t455 + Ifges(7,5) * t476;
t454 = Ifges(6,4) * t444 + Ifges(6,2) * t445 + Ifges(6,6) * t476;
t453 = Ifges(6,1) * t444 + Ifges(6,4) * t445 + Ifges(6,5) * t476;
t452 = -t496 / 0.2e1;
t451 = t496 / 0.2e1;
t449 = t495 / 0.2e1;
t447 = pkin(1) * mrSges(3,1);
t446 = pkin(1) * mrSges(3,2);
t443 = -t355 / 0.2e1;
t442 = t355 / 0.2e1;
t440 = t354 / 0.2e1;
t439 = t199 / 0.2e1;
t438 = -t200 / 0.2e1;
t437 = t200 / 0.2e1;
t434 = -t259 / 0.2e1;
t431 = t279 / 0.2e1;
t429 = t286 / 0.2e1;
t426 = -t309 / 0.2e1;
t418 = Ifges(6,3) + Ifges(7,3);
t288 = pkin(7) * t389;
t383 = qJD(3) * t288 + t268 * t375;
t72 = (-qJ(5) * t376 - t264) * t309 + (qJ(5) * t375 + t321 - t371) * t307 + t383;
t384 = t306 * t264 + t268 * t373;
t73 = (-pkin(7) * qJD(2) + qJ(5) * qJD(3)) * t390 + (qJD(5) * t307 + (-pkin(7) * qJD(3) + qJ(5) * qJD(2)) * t310) * t306 + t384 + t463;
t32 = t302 * t72 + t303 * t73;
t416 = mrSges(4,3) * t258;
t415 = mrSges(4,3) * t259;
t118 = -t209 * t308 - t210 * t305;
t163 = -t305 * t322 - t308 * t323;
t83 = -qJD(6) * t163 - t229 * t308 - t230 * t305;
t400 = t118 - t83;
t119 = -t209 * t305 + t210 * t308;
t162 = t305 * t323 - t308 * t322;
t82 = qJD(6) * t162 - t229 * t305 + t230 * t308;
t399 = t119 - t82;
t396 = qJ(5) * t307;
t395 = qJD(2) * mrSges(3,2);
t298 = t306 * qJ(4);
t287 = pkin(7) * t391;
t301 = t310 * pkin(3);
t159 = pkin(4) * t310 + t287 + t301 + (-t268 - t396) * t309;
t212 = t306 * t268 + t288;
t197 = -qJ(4) * t310 + t212;
t173 = t306 * t396 + t197;
t80 = t302 * t159 + t303 * t173;
t203 = mrSges(4,2) * t286 - t416;
t206 = -mrSges(5,2) * t258 - mrSges(5,3) * t286;
t388 = -t203 - t206;
t204 = -mrSges(4,1) * t286 - t415;
t205 = mrSges(5,1) * t286 + mrSges(5,2) * t259;
t387 = -t204 + t205;
t181 = t259 * pkin(3) + t258 * qJ(4);
t381 = qJ(4) * t360 + qJD(4) * t390;
t267 = -t309 * pkin(3) - pkin(2) - t298;
t7 = -t30 * mrSges(7,1) + t29 * mrSges(7,2);
t31 = -t302 * t73 + t303 * t72;
t53 = -t110 * mrSges(6,1) + t111 * mrSges(6,2);
t79 = t303 * t159 - t173 * t302;
t211 = t268 * t309 - t287;
t244 = t309 * pkin(4) - t267;
t351 = -pkin(7) - t365;
t348 = m(4) * t273 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t258 + mrSges(4,2) * t259 + mrSges(3,3) * t379;
t137 = -pkin(4) * t259 - t181;
t344 = -t264 * t309 + t383;
t343 = mrSges(4,1) * t309 - mrSges(4,2) * t306;
t341 = mrSges(5,1) * t309 + mrSges(5,3) * t306;
t334 = Ifges(4,2) * t309 + t413;
t330 = -Ifges(5,3) * t309 + t409;
t329 = t302 * t44 - t303 * t45;
t220 = t322 * t307;
t57 = pkin(5) * t310 - pkin(9) * t220 + t79;
t219 = t323 * t307;
t61 = -pkin(9) * t219 + t80;
t24 = -t305 * t61 + t308 * t57;
t25 = t305 * t57 + t308 * t61;
t328 = t306 * t92 + t309 * t78;
t327 = -t100 * t306 + t309 * t99;
t135 = -t219 * t308 - t220 * t305;
t136 = -t219 * t305 + t220 * t308;
t202 = -pkin(7) * t361 + t235;
t165 = -mrSges(5,1) * t356 + t199 * mrSges(5,2);
t133 = -mrSges(6,2) * t286 + mrSges(6,3) * t355;
t134 = mrSges(6,1) * t286 - mrSges(6,3) * t354;
t319 = t133 * t303 - t134 * t302 + t206;
t285 = qJ(4) * t390;
t196 = t307 * t351 + t285;
t84 = t200 * pkin(3) - t199 * qJ(4) + qJD(2) * t295 - t259 * qJD(4);
t58 = -pkin(4) * t200 - t84;
t128 = (-t307 * t370 - t310 * t375) * pkin(7) + t384;
t98 = (-t309 * t448 - t298) * t374 + t351 * t376 + t381;
t284 = Ifges(5,2) * t356;
t283 = Ifges(4,3) * t356;
t270 = mrSges(3,3) * t378 - t395;
t218 = -t285 + (pkin(3) * t306 + pkin(7)) * t307;
t201 = pkin(7) * t363 + t392;
t198 = -t211 + t301;
t194 = Ifges(5,4) * t199;
t193 = Ifges(4,5) * t199;
t192 = Ifges(4,6) * t200;
t191 = Ifges(5,6) * t200;
t184 = pkin(5) * t322 + t244;
t182 = mrSges(5,1) * t258 - mrSges(5,3) * t259;
t180 = qJD(1) * t347 - t392;
t179 = t202 + t290;
t166 = -mrSges(4,2) * t356 - mrSges(4,3) * t200;
t164 = mrSges(4,1) * t356 - mrSges(4,3) * t199;
t161 = -mrSges(5,2) * t200 + mrSges(5,3) * t356;
t141 = qJD(2) * t317 + t323 * t374;
t140 = -qJD(2) * t318 + t230 * t307;
t132 = pkin(5) * t219 + t196;
t129 = t377 * t423 - t344;
t127 = pkin(3) * t316 + pkin(7) * t376 + qJ(4) * t359 - t381;
t117 = t320 + t344;
t116 = mrSges(4,1) * t200 + mrSges(4,2) * t199;
t115 = mrSges(5,1) * t200 - mrSges(5,3) * t199;
t109 = t128 + t463;
t102 = t199 * Ifges(4,4) - t200 * Ifges(4,2) + Ifges(4,6) * t356;
t101 = t199 * Ifges(5,5) + Ifges(5,6) * t356 + t200 * Ifges(5,3);
t96 = -mrSges(6,1) * t356 - mrSges(6,3) * t111;
t95 = mrSges(6,2) * t356 + mrSges(6,3) * t110;
t91 = -mrSges(6,1) * t355 + mrSges(6,2) * t354;
t77 = -pkin(5) * t354 + t137;
t75 = Ifges(6,2) * t355 + t286 * Ifges(6,6) + t473;
t68 = mrSges(7,1) * t279 - mrSges(7,3) * t495;
t67 = -mrSges(7,2) * t279 + mrSges(7,3) * t496;
t56 = -pkin(5) * t140 + t98;
t41 = -qJD(6) * t136 + t140 * t308 - t141 * t305;
t40 = qJD(6) * t135 + t140 * t305 + t141 * t308;
t39 = -mrSges(7,1) * t496 + mrSges(7,2) * t495;
t38 = -pkin(5) * t110 + t58;
t35 = Ifges(7,2) * t496 + Ifges(7,6) * t279 + t424;
t23 = mrSges(7,2) * t356 + mrSges(7,3) * t30;
t22 = -mrSges(7,1) * t356 - mrSges(7,3) * t29;
t21 = pkin(9) * t140 + t32;
t18 = -pkin(5) * t377 - pkin(9) * t141 + t31;
t4 = -qJD(6) * t25 + t18 * t308 - t21 * t305;
t3 = qJD(6) * t24 + t18 * t305 + t21 * t308;
t5 = [t40 * t506 + m(7) * (t1 * t25 + t132 * t38 + t2 * t24 + t3 * t9 + t4 * t8 + t56 * t66) + m(6) * (t114 * t98 + t14 * t79 + t15 * t80 + t196 * t58 + t31 * t44 + t32 * t45) + m(5) * (t109 * t148 + t117 * t146 + t127 * t152 + t197 * t78 + t198 * t92 + t218 * t84) + t141 * t491 + (-t14 * t220 + t140 * t45 - t141 * t44 - t15 * t219) * mrSges(6,3) + (Ifges(6,1) * t220 - Ifges(6,4) * t219) * t444 + (Ifges(6,4) * t220 - Ifges(6,2) * t219) * t445 + t58 * (mrSges(6,1) * t219 + mrSges(6,2) * t220) + (pkin(7) * t116 + t101 * t427 + t84 * t340 + t331 * t437 + t335 * t438 + (-t100 * t309 - t306 * t99) * mrSges(4,3) + (-t306 * t78 + t309 * t92) * mrSges(5,2) + (t273 * t343 + t334 * t435 + t330 * t436 + t152 * t341 + t156 * t426 + (t176 * t306 - t177 * t309) * mrSges(4,3) + (-t146 * t306 - t148 * t309) * mrSges(5,2) + t479 * t434 + t480 * t429) * qJD(3) + (-pkin(7) * t270 + t357 + ((-0.3e1 / 0.2e1 * Ifges(3,4) + t411 / 0.2e1 + t406 / 0.2e1 + t410 / 0.2e1 - t407 / 0.2e1) * t307 - Ifges(6,5) * t220 / 0.2e1 + Ifges(6,6) * t219 / 0.2e1 - Ifges(7,5) * t136 / 0.2e1 - Ifges(7,6) * t135 / 0.2e1 - 0.2e1 * t447 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + (m(4) * pkin(7) + t342) * pkin(7) - t367 - t418) * t310) * qJD(1) + t459) * qJD(2) + t462 * t439 + (qJD(3) * t466 + t102) * t428 + (qJD(3) * t153 + t484) * t425) * t307 + (-t283 / 0.2e1 - t284 / 0.2e1 + ((0.3e1 / 0.2e1 * Ifges(3,4) * t310 - 0.2e1 * t446) * qJD(1) + t358 + t348 * pkin(7) + t477) * qJD(2) - t193 / 0.2e1 - t194 / 0.2e1 - t191 / 0.2e1 + t192 / 0.2e1 - t368 * t199 + t366 * t200 + t460) * t310 + (t1 * t135 - t136 * t2 - t40 * t8 + t41 * t9) * mrSges(7,3) + (Ifges(7,4) * t40 + Ifges(7,2) * t41) * t451 + t220 * t453 - t219 * t454 + (Ifges(7,4) * t136 + Ifges(7,2) * t135) * t455 + (Ifges(7,1) * t136 + Ifges(7,4) * t135) * t456 + t136 * t457 + t135 * t458 + (Ifges(6,1) * t141 + Ifges(6,4) * t140) * t440 + (Ifges(6,4) * t141 + Ifges(6,2) * t140) * t442 + (Ifges(7,1) * t40 + Ifges(7,4) * t41) * t449 + (Ifges(7,5) * t40 + Ifges(7,6) * t41) * t431 + (Ifges(6,5) * t141 + Ifges(6,6) * t140) * t429 + m(4) * (t100 * t211 + t128 * t177 + t129 * t176 + t212 * t99) + t24 * t22 + t25 * t23 + t41 * t35 / 0.2e1 + t56 * t39 + t66 * (-mrSges(7,1) * t41 + mrSges(7,2) * t40) + t3 * t67 + t4 * t68 + t80 * t95 + t79 * t96 + t98 * t91 + t132 * t7 + t32 * t133 + t31 * t134 + t38 * (-mrSges(7,1) * t135 + mrSges(7,2) * t136) + t140 * t75 / 0.2e1 + t114 * (-mrSges(6,1) * t140 + mrSges(6,2) * t141) + t127 * t182 + t196 * t53 + t197 * t161 + t198 * t165 + t128 * t203 + t129 * t204 + t117 * t205 + t109 * t206 + t211 * t164 + t212 * t166 + t218 * t115; (t14 * t323 - t15 * t322 + t385 * t44 + t45 * t511) * mrSges(6,3) + (-mrSges(6,1) * t511 - mrSges(6,2) * t385) * t114 + (Ifges(6,4) * t210 - Ifges(6,2) * t209) * t443 + (Ifges(6,5) * t210 - Ifges(6,6) * t209) * t430 + (-Ifges(6,1) * t323 - Ifges(6,4) * t322) * t444 + (-Ifges(6,4) * t323 - Ifges(6,2) * t322) * t445 + t58 * (mrSges(6,1) * t322 - mrSges(6,2) * t323) - t323 * t453 + (t209 / 0.2e1 - t229 / 0.2e1) * t75 + (Ifges(6,1) * t210 - Ifges(6,4) * t209) * t441 - t322 * t454 - m(4) * (t176 * t201 + t177 * t202) + (-t210 / 0.2e1 + t230 / 0.2e1) * t76 + t487 * t67 + (t1 * t63 + t184 * t38 + t2 * t62 + t483 * t66 + t487 * t9 + t488 * t8) * m(7) + t488 * t68 + (Ifges(6,1) * t230 - Ifges(6,4) * t229) * t440 + (Ifges(6,4) * t230 - Ifges(6,2) * t229) * t442 + (Ifges(6,5) * t230 - Ifges(6,6) * t229) * t429 + t479 * t439 + t481 * t182 + (-t146 * t180 - t148 * t179 + t481 * t152 + t267 * t84) * m(5) + t482 * t91 + t483 * t39 + t484 * t427 + t485 * t133 + t486 * t134 + (t114 * t482 + t14 * t187 + t15 * t188 + t244 * t58 + t44 * t486 + t45 * t485) * m(6) + (t82 / 0.2e1 - t119 / 0.2e1) * t36 + (Ifges(7,4) * t82 + Ifges(7,2) * t83) * t451 + (Ifges(7,4) * t119 + Ifges(7,2) * t118) * t452 + (Ifges(7,4) * t163 + Ifges(7,2) * t162) * t455 + (Ifges(7,1) * t163 + Ifges(7,4) * t162) * t456 + t163 * t457 + t162 * t458 + (Ifges(7,1) * t82 + Ifges(7,4) * t83) * t449 + (Ifges(7,1) * t119 + Ifges(7,4) * t118) * t450 + (Ifges(7,5) * t82 + Ifges(7,6) * t83) * t431 + (Ifges(7,5) * t119 + Ifges(7,6) * t118) * t432 + t330 * t437 + t334 * t438 + t102 * t425 + t101 * t426 + (mrSges(7,1) * t400 - mrSges(7,2) * t399) * t66 + (t1 * t162 - t163 * t2 + t399 * t8 - t400 * t9) * mrSges(7,3) + (t83 / 0.2e1 - t118 / 0.2e1) * t35 + t267 * t115 + t244 * t53 + ((t357 + (t447 + t414 / 0.2e1) * qJD(1) + (t270 + t395) * pkin(7) + (-Ifges(6,5) * t323 + Ifges(7,5) * t163 - Ifges(6,6) * t322 + Ifges(7,6) * t162) * t490 + t480 * t501 - t459) * t307 + (-t293 / 0.2e1 + (t446 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t307) * qJD(1) + t358 + ((-m(4) * pkin(2) - mrSges(3,1) - t343) * qJD(2) - t348) * pkin(7) - t477) * t310) * qJD(1) + t327 * mrSges(4,3) + t328 * mrSges(5,2) + t62 * t22 + t63 * t23 - t84 * t341 - pkin(2) * t116 + t38 * (-mrSges(7,1) * t162 + mrSges(7,2) * t163) + t184 * t7 + t187 * t96 + t188 * t95 - t202 * t203 - t201 * t204 - t180 * t205 - t179 * t206 + t312 * qJD(3) + ((t161 + t166) * t309 + (-t164 + t165) * t306 + (-m(4) * t325 + m(5) * t326 + t306 * t388 + t309 * t387) * qJD(3) + m(5) * t328 + m(4) * t327) * pkin(8); (t35 - t508) * t450 - t460 + (-Ifges(6,5) * t355 + Ifges(6,6) * t354) * t430 + (-t354 * t45 - t355 * t44) * mrSges(6,3) - t114 * (-mrSges(6,1) * t354 - mrSges(6,2) * t355) + (-Ifges(6,1) * t355 + t473 + t75) * t441 + t355 * t491 + (Ifges(6,2) * t354 - t492) * t443 + t515 * t452 + t283 + t284 + (t146 * t258 + t148 * t259) * mrSges(5,2) + t193 + t194 + t191 - t192 + t496 * t506 + (-Ifges(4,2) * t259 - t246 + t466) * t435 + (-pkin(3) * t92 + qJ(4) * t78 - t146 * t177 + t148 * t461 - t152 * t181) * m(5) + (-qJD(4) * t329 - t114 * t137 + t14 * t265 + t15 * t266 - t44 * t59 - t45 * t60) * m(6) + t471 * t67 + t472 * t68 + (t1 * t175 + t174 * t2 + t471 * t9 + t472 * t8 - t66 * t77) * m(7) + t418 * t356 + t156 * t433 + (Ifges(5,3) * t259 - t403) * t436 - t512 + (-t387 + t415) * t177 + (t388 - t416) * t176 - t273 * (mrSges(4,1) * t259 - mrSges(4,2) * t258) + t265 * t96 + t266 * t95 - t152 * (mrSges(5,1) * t259 + mrSges(5,3) * t258) + (-t489 * t258 - t259 * t500) * t429 + (-t258 * t493 + t153 + t245 - t402) * t434 + t319 * qJD(4) - t77 * t39 - t60 * t133 - t59 * t134 - t137 * t91 + qJ(4) * t161 - pkin(3) * t165 + t174 * t22 + t175 * t23 - t181 * t182; -t324 * t22 + t253 * t23 + t302 * t95 + t303 * t96 - t464 * t68 - t465 * t67 + t319 * t286 + (t182 - t39 - t91) * t259 + t165 + (t1 * t253 - t2 * t324 - t259 * t66 - t464 * t8 - t465 * t9) * m(7) + (-t114 * t259 + t14 * t303 + t15 * t302 - t286 * t329) * m(6) + (t148 * t286 + t152 * t259 + t92) * m(5); -t355 * t133 + t354 * t134 - t496 * t67 + t495 * t68 + t53 + t7 + (t495 * t8 - t496 * t9 + t38) * m(7) + (t354 * t44 - t355 * t45 + t58) * m(6); -Ifges(7,3) * t356 + t508 * t450 + t35 * t449 - t8 * t67 + t9 * t68 + (t36 - t515) * t452 + t478 + t512;];
tauc  = t5(:);
