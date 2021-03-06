% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 06:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRPPR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:39:05
% EndTime: 2019-05-07 06:39:48
% DurationCPUTime: 21.34s
% Computational Cost: add. (84203->645), mult. (185131->896), div. (0->0), fcn. (147640->12), ass. (0->412)
t361 = sin(pkin(11));
t363 = cos(pkin(11));
t364 = cos(pkin(6));
t449 = qJD(1) * t364;
t356 = qJD(2) + t449;
t366 = sin(qJ(3));
t369 = cos(qJ(3));
t362 = sin(pkin(6));
t367 = sin(qJ(2));
t448 = qJD(1) * t367;
t427 = t362 * t448;
t333 = t356 * t366 + t369 * t427;
t370 = cos(qJ(2));
t447 = qJD(1) * t370;
t426 = t362 * t447;
t349 = -qJD(3) + t426;
t308 = t361 * t333 + t349 * t363;
t310 = t363 * t333 - t361 * t349;
t273 = t310 * t308;
t439 = qJDD(1) * t362;
t339 = qJD(2) * t426 + t367 * t439;
t355 = qJDD(1) * t364 + qJDD(2);
t417 = t339 * t366 - t369 * t355;
t292 = qJD(3) * t333 + t417;
t525 = t273 + t292;
t478 = t525 * t363;
t307 = t310 ^ 2;
t331 = -t369 * t356 + t366 * t427;
t508 = t331 ^ 2;
t524 = -t307 - t508;
t166 = t361 * t524 + t478;
t405 = -t339 * t369 - t355 * t366;
t293 = -qJD(3) * t331 - t405;
t438 = qJDD(1) * t370;
t389 = (qJD(2) * t448 - t438) * t362;
t380 = qJDD(3) + t389;
t268 = t363 * t293 + t361 * t380;
t466 = t308 * t331;
t526 = t268 - t466;
t117 = t166 * t369 - t366 * t526;
t479 = t525 * t361;
t164 = -t363 * t524 + t479;
t601 = pkin(1) * (t117 * t367 - t164 * t370);
t115 = t166 * t366 + t369 * t526;
t600 = -pkin(8) * (t117 * t370 + t367 * t164) + pkin(1) * t115;
t598 = pkin(2) * t115;
t597 = pkin(9) * t115;
t596 = pkin(2) * t164 - pkin(9) * t117;
t419 = t293 * t361 - t363 * t380;
t463 = t331 * t310;
t398 = t419 + t463;
t482 = t526 * t361;
t154 = t398 * t363 + t482;
t509 = t308 ^ 2;
t270 = t307 - t509;
t119 = t154 * t366 + t270 * t369;
t151 = -t398 * t361 + t363 * t526;
t593 = t364 * t119 + (t370 * t151 + t367 * (t154 * t369 - t270 * t366)) * t362;
t281 = t509 - t508;
t188 = t281 * t363 - t479;
t399 = -t419 + t463;
t140 = t188 * t366 - t399 * t369;
t184 = t281 * t361 + t478;
t592 = t364 * t140 + (t367 * (t188 * t369 + t399 * t366) - t370 * t184) * t362;
t212 = t466 + t268;
t282 = -t307 + t508;
t420 = -t292 + t273;
t476 = t420 * t363;
t570 = -t282 * t361 - t476;
t477 = t420 * t361;
t571 = t282 * t363 - t477;
t579 = -t212 * t369 + t366 * t570;
t591 = t364 * t579 + (t367 * (t212 * t366 + t369 * t570) - t370 * t571) * t362;
t587 = pkin(3) * t164;
t586 = qJ(4) * t164;
t585 = qJ(4) * t166;
t521 = -t508 - t509;
t536 = t363 * t521 + t477;
t551 = t366 * t536 - t369 * t398;
t577 = pkin(2) * t551;
t230 = -t509 - t307;
t515 = t212 * t361 + t363 * t399;
t552 = -t230 * t369 + t366 * t515;
t576 = pkin(2) * t552;
t575 = pkin(9) * t551;
t574 = pkin(9) * t552;
t553 = t230 * t366 + t369 * t515;
t573 = pkin(9) * t553;
t537 = t361 * t521 - t476;
t550 = t366 * t398 + t369 * t536;
t569 = -pkin(2) * t537 + pkin(9) * t550;
t516 = -t212 * t363 + t361 * t399;
t568 = pkin(1) * (t367 * t553 - t370 * t516);
t567 = pkin(1) * (t367 * t550 - t370 * t537);
t566 = pkin(8) * (t367 * t516 + t370 * t553) - pkin(1) * t552;
t565 = pkin(8) * (t367 * t537 + t370 * t550) - pkin(1) * t551;
t563 = pkin(3) * t537;
t562 = qJ(4) * t536;
t561 = qJ(4) * t537;
t554 = -pkin(3) * t230 + qJ(4) * t515;
t549 = 2 * qJD(5);
t547 = qJ(4) * t516;
t546 = qJ(5) * t526;
t462 = t331 * t361;
t203 = t308 * t462 - t363 * t419;
t461 = t331 * t363;
t432 = t308 * t461;
t387 = t361 * t419 + t432;
t434 = t366 * t273;
t433 = t369 * t273;
t513 = t366 * t387 + t433;
t535 = t364 * t513 + (-t370 * t203 + t367 * (t369 * t387 - t434)) * t362;
t276 = t310 * t462;
t411 = t276 - t432;
t467 = t292 * t366;
t288 = t369 * t292;
t517 = t366 * t411 - t288;
t527 = (t308 * t361 + t310 * t363) * t331;
t534 = t364 * t517 + (t370 * t527 + t367 * (t369 * t411 + t467)) * t362;
t365 = sin(qJ(6));
t290 = qJDD(6) - t292;
t368 = cos(qJ(6));
t261 = -t368 * t308 + t310 * t365;
t263 = t308 * t365 + t310 * t368;
t472 = t263 * t261;
t379 = t290 - t472;
t533 = t365 * t379;
t460 = t333 * t331;
t376 = t380 - t460;
t532 = t366 * t376;
t530 = t368 * t379;
t529 = t369 * t376;
t316 = t331 * t349;
t254 = t293 + t316;
t269 = pkin(4) * t308 - qJ(5) * t310;
t502 = sin(qJ(1));
t503 = cos(qJ(1));
t396 = g(1) * t503 + g(2) * t502;
t506 = qJD(1) ^ 2;
t335 = -pkin(1) * t506 + pkin(8) * t439 - t396;
t501 = pkin(2) * t370;
t414 = -pkin(9) * t367 - t501;
t450 = qJD(1) * t362;
t338 = t414 * t450;
t354 = t356 ^ 2;
t395 = g(1) * t502 - g(2) * t503;
t443 = t362 * t506;
t375 = qJDD(1) * pkin(1) + pkin(8) * t443 + t395;
t373 = t364 * t375;
t457 = t362 * t367;
t372 = -g(3) * t457 + t367 * t373;
t246 = t355 * pkin(9) - t354 * pkin(2) + (t338 * t450 + t335) * t370 + t372;
t496 = t364 * g(3);
t371 = -t339 * pkin(9) - t496 + ((-pkin(1) - t501) * qJDD(1) + (-pkin(8) * t450 - t370 * t356 * pkin(9) + (qJD(2) + t356) * t367 * pkin(2)) * qJD(1) - t395) * t362;
t191 = t369 * t246 + t366 * t371;
t297 = pkin(3) * t331 - qJ(4) * t333;
t507 = t349 ^ 2;
t149 = -pkin(3) * t507 + qJ(4) * t380 - t331 * t297 + t191;
t418 = t367 * t335 - t370 * t373;
t245 = -t355 * pkin(2) - t354 * pkin(9) + (g(3) * t370 + t338 * t448) * t362 + t418;
t160 = -t254 * qJ(4) + (-t333 * t349 + t292) * pkin(3) + t245;
t451 = t363 * t149 + t361 * t160;
t523 = t292 * qJ(5) - t308 * t269 + t331 * t549 + t451;
t456 = t362 * t370;
t295 = g(3) * t456 + t418;
t296 = t370 * t335 + t372;
t522 = t367 * t295 + t370 * t296;
t422 = t149 * t361 - t363 * t160;
t394 = -t292 * pkin(4) - qJ(5) * t508 + qJDD(5) + t422;
t378 = -pkin(10) * t212 + t394;
t505 = 2 * qJD(4);
t440 = t505 + t269;
t520 = t378 - pkin(5) * t292 + (pkin(5) * t308 + t440) * t310;
t423 = qJ(5) * t361 + pkin(3);
t446 = qJD(4) * t308;
t303 = -0.2e1 * t446;
t400 = t303 + t523;
t77 = -pkin(4) * t508 + t400;
t78 = t310 * t440 + t394;
t43 = t361 * t78 + t363 * t77;
t499 = pkin(4) * t363;
t190 = t366 * t246 - t369 * t371;
t148 = -t380 * pkin(3) - t507 * qJ(4) + t333 * t297 + qJDD(4) + t190;
t377 = t419 * pkin(4) + t148 - t546;
t96 = (pkin(4) * t331 - (2 * qJD(5))) * t310 + t377;
t519 = -(t423 + t499) * t96 + qJ(4) * t43;
t374 = t310 * t549 - t377;
t79 = -pkin(4) * t463 + t374 + t546;
t518 = t585 + t526 * (pkin(3) + t499) + t361 * t79;
t249 = (qJD(3) + t349) * t333 + t417;
t206 = t268 * t361 + t310 * t461;
t207 = t268 * t363 - t276;
t412 = t366 * t207 - t433;
t514 = t364 * t412 + (t369 * t207 + t434) * t457 - t206 * t456;
t80 = (-t398 - t463) * pkin(4) + t374;
t512 = t363 * t80 - t398 * t423 + t562;
t259 = t261 ^ 2;
t260 = t263 ^ 2;
t325 = -qJD(6) + t331;
t322 = t325 ^ 2;
t330 = t333 ^ 2;
t358 = t362 ^ 2;
t504 = pkin(4) + pkin(5);
t500 = pkin(3) * t366;
t497 = pkin(8) * t362;
t413 = -pkin(5) * t331 - pkin(10) * t310;
t75 = pkin(5) * t419 + pkin(10) * t509 - t310 * t413 + t96;
t494 = t365 * t75;
t67 = -pkin(5) * t509 + pkin(10) * t419 + t331 * t413 + t77;
t493 = t368 * t67;
t492 = t368 * t75;
t491 = qJ(5) * t363;
t490 = t148 * t361;
t489 = t148 * t363;
t172 = -t290 - t472;
t488 = t172 * t365;
t487 = t172 * t368;
t475 = t245 * t366;
t474 = t245 * t369;
t473 = t261 * t325;
t277 = -t380 - t460;
t471 = t277 * t366;
t470 = t277 * t369;
t465 = t325 * t365;
t464 = t325 * t368;
t459 = t349 * t366;
t458 = t349 * t369;
t444 = t358 * t506;
t348 = t370 * t367 * t444;
t336 = t348 + t355;
t454 = t367 * t336;
t337 = -t348 + t355;
t452 = t370 * t337;
t442 = qJD(3) - t349;
t437 = t310 * t505;
t436 = t366 * t472;
t435 = t369 * t472;
t431 = t370 * t460;
t99 = t303 + t451;
t162 = -t261 * qJD(6) + t368 * t268 + t365 * t419;
t343 = t356 * t426;
t430 = t343 + t339;
t429 = -pkin(3) * t369 - pkin(2);
t359 = t367 ^ 2;
t425 = t359 * t444;
t360 = t370 ^ 2;
t424 = t360 * t444;
t98 = t422 + t437;
t57 = t361 * t98 + t363 * t99;
t34 = t365 * t67 - t368 * t520;
t114 = t190 * t366 + t369 * t191;
t421 = t268 * t365 - t368 * t419;
t416 = -t356 + t449;
t410 = -pkin(3) * t148 + qJ(4) * t57;
t35 = t365 * t520 + t493;
t19 = -t34 * t368 + t35 * t365;
t20 = t34 * t365 + t35 * t368;
t56 = t361 * t99 - t363 * t98;
t407 = t190 * t369 - t191 * t366;
t404 = -pkin(1) + t414;
t402 = t162 + t473;
t393 = -pkin(3) * t398 - t489 + t562;
t392 = -pkin(3) * t526 + t490 - t585;
t391 = (-qJD(6) - t325) * t263 - t421;
t390 = t554 + t57;
t12 = t19 * t361 + t20 * t363;
t13 = -pkin(10) * t20 - t504 * t75;
t17 = -pkin(10) * t19 - qJ(5) * t75;
t386 = -pkin(3) * t75 + qJ(4) * t12 + t13 * t363 + t17 * t361;
t129 = (qJD(6) - t325) * t263 + t421;
t192 = -t322 - t259;
t108 = t192 * t368 - t533;
t37 = -pkin(10) * t108 + t129 * t504 - t492;
t107 = t192 * t365 + t530;
t45 = -pkin(10) * t107 + qJ(5) * t129 - t494;
t74 = t107 * t361 + t108 * t363;
t385 = pkin(3) * t129 + qJ(4) * t74 + t361 * t45 + t363 * t37;
t222 = -t260 - t322;
t136 = -t222 * t365 + t487;
t39 = -pkin(10) * t136 + t402 * t504 + t494;
t135 = t222 * t368 + t488;
t46 = -pkin(10) * t135 + qJ(5) * t402 - t492;
t86 = t135 * t361 + t136 * t363;
t384 = pkin(3) * t402 + qJ(4) * t86 + t361 * t46 + t363 * t39;
t163 = -t259 - t260;
t133 = -t162 + t473;
t90 = -t133 * t365 + t368 * t391;
t14 = -pkin(10) * t90 + t163 * t504 - t20;
t88 = t133 * t368 + t365 * t391;
t15 = -pkin(10) * t88 + qJ(5) * t163 - t19;
t53 = t361 * t88 + t363 * t90;
t383 = pkin(3) * t163 + qJ(4) * t53 + t14 * t363 + t15 * t361;
t68 = (-t230 - t508) * pkin(4) + t400;
t69 = -qJ(5) * t230 + t78;
t382 = t361 * t69 + t363 * t68 + t554;
t342 = t356 * t427;
t341 = (t359 - t360) * t444;
t340 = -t354 - t424;
t323 = -t425 - t354;
t317 = t362 * t375 + t496;
t315 = -t342 - t389;
t314 = t342 - t389;
t313 = -t343 + t339;
t312 = -t330 + t507;
t311 = t508 - t507;
t300 = -t330 - t507;
t298 = t330 - t508;
t294 = -t508 - t507;
t275 = t508 + t330;
t274 = (t331 * t366 + t333 * t369) * t349;
t253 = t331 * t442 + t405;
t252 = t293 - t316;
t250 = -t333 * t442 - t417;
t248 = t293 * t366 - t333 * t458;
t247 = -t331 * t459 - t288;
t239 = t311 * t366 - t470;
t238 = t312 * t369 + t532;
t235 = -t260 + t322;
t234 = t259 - t322;
t233 = -t300 * t366 + t470;
t232 = t300 * t369 + t471;
t228 = t294 * t369 - t532;
t227 = t294 * t366 + t529;
t198 = t260 - t259;
t197 = -t249 * t369 + t252 * t366;
t195 = t250 * t366 + t254 * t369;
t180 = (t261 * t368 - t263 * t365) * t325;
t179 = (-t261 * t365 - t263 * t368) * t325;
t161 = -qJD(6) * t263 - t421;
t147 = pkin(2) * t253 + pkin(9) * t233 + t475;
t145 = pkin(2) * t250 + pkin(9) * t228 - t474;
t144 = t234 * t368 + t488;
t143 = -t235 * t365 + t530;
t142 = -t234 * t365 + t487;
t141 = -t235 * t368 - t533;
t128 = t162 * t368 + t263 * t465;
t127 = -t162 * t365 + t263 * t464;
t126 = -t161 * t365 - t261 * t464;
t125 = -t161 * t368 + t261 * t465;
t106 = -t179 * t361 + t180 * t363;
t105 = t179 * t363 + t180 * t361;
t104 = -pkin(2) * t245 + pkin(9) * t114;
t103 = t106 * t366 + t290 * t369;
t102 = t489 + t586;
t101 = t490 - t561;
t100 = -pkin(3) * t516 + pkin(4) * t212 - qJ(5) * t399;
t95 = pkin(2) * t275 + pkin(9) * t197 + t114;
t94 = -t142 * t361 + t144 * t363;
t93 = -t141 * t361 + t143 * t363;
t92 = t142 * t363 + t144 * t361;
t91 = t141 * t363 + t143 * t361;
t89 = -t129 * t368 - t365 * t402;
t87 = t129 * t365 - t368 * t402;
t85 = -t135 * t363 + t136 * t361;
t84 = -t127 * t361 + t128 * t363;
t83 = -t125 * t361 + t126 * t363;
t82 = t127 * t363 + t128 * t361;
t81 = t125 * t363 + t126 * t361;
t76 = t99 + t587;
t73 = -t107 * t363 + t108 * t361;
t72 = t98 - t563;
t71 = t366 * t84 + t435;
t70 = t366 * t83 - t435;
t66 = t366 * t94 + t369 * t391;
t65 = -t133 * t369 + t366 * t93;
t63 = -t366 * t402 + t369 * t86;
t62 = t366 * t86 + t369 * t402;
t61 = -t361 * t80 - t398 * t491 - t561;
t60 = -pkin(4) * t482 + t363 * t79 - t586;
t59 = -t129 * t366 + t369 * t74;
t58 = t129 * t369 + t366 * t74;
t55 = pkin(4) * t420 - qJ(5) * t521 - t563 + t78;
t54 = -t587 - qJ(5) * t525 + 0.2e1 * t446 + (t524 + t508) * pkin(4) - t523;
t52 = -t361 * t87 + t363 * t89;
t51 = t361 * t90 - t363 * t88;
t50 = t361 * t89 + t363 * t87;
t49 = -t56 - t547;
t48 = t148 * t366 + t369 * t57;
t47 = -t148 * t369 + t366 * t57;
t44 = t198 * t369 + t366 * t52;
t42 = t361 * t77 - t363 * t78;
t41 = -t163 * t366 + t369 * t53;
t40 = t163 * t369 + t366 * t53;
t38 = t102 * t366 + t369 * t76 + t596;
t36 = t101 * t366 + t369 * t72 + t569;
t33 = -t361 * t68 + t363 * t69 - t547;
t32 = t366 * t96 + t369 * t43;
t31 = t366 * t43 - t369 * t96;
t30 = t366 * t49 + t429 * t516 + t573;
t29 = -qJ(4) * t42 + (pkin(4) * t361 - t491) * t96;
t28 = t366 * t61 + t369 * t55 + t569;
t27 = t366 * t60 + t369 * t54 - t596;
t26 = -pkin(3) * t42 + pkin(4) * t78 - qJ(5) * t77;
t25 = -pkin(2) * t516 + t100 * t369 + t33 * t366 + t573;
t24 = -pkin(3) * t51 - qJ(5) * t90 + t504 * t88;
t23 = -pkin(3) * t85 + pkin(4) * t135 - qJ(5) * t136 - t493 - t365 * (t269 * t310 + t378 + t437) + (-t365 * t420 + t135) * pkin(5);
t22 = -qJ(4) * t85 - t361 * t39 + t363 * t46;
t21 = -pkin(3) * t73 - qJ(5) * t108 + t107 * t504 - t34;
t18 = -qJ(4) * t73 - t361 * t37 + t363 * t45;
t16 = pkin(9) * t48 + (-qJ(4) * t366 + t429) * t56;
t11 = -t19 * t363 + t20 * t361;
t10 = t12 * t369 + t366 * t75;
t9 = t12 * t366 - t369 * t75;
t8 = -pkin(2) * t85 + pkin(9) * t63 + t22 * t366 + t23 * t369;
t7 = -pkin(2) * t42 + pkin(9) * t32 + t26 * t369 + t29 * t366;
t6 = -pkin(2) * t73 + pkin(9) * t59 + t18 * t366 + t21 * t369;
t5 = -qJ(4) * t51 - t14 * t361 + t15 * t363;
t4 = -pkin(2) * t51 + pkin(9) * t41 + t24 * t369 + t366 * t5;
t3 = -qJ(4) * t11 - t13 * t361 + t17 * t363;
t2 = -pkin(3) * t11 - qJ(5) * t20 + t19 * t504;
t1 = -pkin(2) * t11 + pkin(9) * t10 + t2 * t369 + t3 * t366;
t64 = [0, 0, 0, 0, 0, qJDD(1), t395, t396, 0, 0, (-t358 * t416 * t447 + t339 * t362) * t367, t364 * t341 + (t367 * t315 + t370 * t430) * t362, t364 * t313 + (t454 + t370 * (t354 - t425)) * t362, (t438 + (-qJD(2) + t416) * t448) * t358 * t370, t364 * t314 + (t367 * (-t354 + t424) + t452) * t362, t364 * t355, (-t295 + pkin(1) * (t336 * t370 + t340 * t367)) * t364 + (t370 * t317 + pkin(1) * t315 + pkin(8) * (t340 * t370 - t454)) * t362, -t317 * t457 - t364 * t296 + pkin(1) * (-t362 * t430 + (t323 * t370 - t337 * t367) * t364) + (-t367 * t323 - t452) * t497, pkin(1) * ((-t313 * t370 + t314 * t367) * t364 - (-t359 - t360) * t358 * t443) + (t367 * t313 + t314 * t370) * t497 + t522 * t362, pkin(1) * (t362 * t317 + (-t295 * t370 + t296 * t367) * t364) + t522 * t497, t364 * t248 + (t367 * (t293 * t369 + t333 * t459) - t431) * t362, t364 * t195 + (t367 * (t250 * t369 - t254 * t366) - t370 * t298) * t362, t364 * t238 + (t367 * (-t312 * t366 + t529) - t370 * t252) * t362, t364 * t247 + (t367 * (-t331 * t458 + t467) + t431) * t362, t364 * t239 + (t367 * (t311 * t369 + t471) + t370 * t249) * t362, -t380 * t456 + t364 * t274 + (t331 * t369 - t333 * t366) * t349 * t457, (t145 + pkin(1) * (t228 * t367 + t250 * t370)) * t364 + (t367 * (-pkin(9) * t227 + t475) + t370 * (-pkin(2) * t227 + t190) - pkin(1) * t227 + pkin(8) * (t228 * t370 - t367 * t250)) * t362, (t147 + pkin(1) * (t233 * t367 + t253 * t370)) * t364 + (t367 * (-pkin(9) * t232 + t474) + t370 * (-pkin(2) * t232 + t191) - pkin(1) * t232 + pkin(8) * (t233 * t370 - t367 * t253)) * t362, (t95 + pkin(1) * (t197 * t367 + t275 * t370)) * t364 + (t367 * t407 + pkin(8) * (t197 * t370 - t367 * t275) + t404 * (-t249 * t366 - t252 * t369)) * t362, (t104 + pkin(1) * (t114 * t367 - t245 * t370)) * t364 + (pkin(8) * (t114 * t370 + t367 * t245) - t404 * t407) * t362, t514, -t593, t591, t535, t592, t534, (t36 + t567) * t364 + (t367 * (t101 * t369 - t366 * t72 - t575) + t370 * (-t393 - t577) + t565) * t362, (t38 - t601) * t364 + (t367 * (t102 * t369 - t366 * t76 + t597) + t370 * (-t392 + t598) + t600) * t362, (t30 + t568) * t364 + (t367 * (t369 * t49 + t500 * t516 - t574) + t370 * (-t390 - t576) + t566) * t362, (t16 + pkin(1) * (t367 * t48 - t370 * t56)) * t364 + (t367 * (-pkin(9) * t47 + (-qJ(4) * t369 + t500) * t56) + t370 * (-pkin(2) * t47 - t410) - pkin(1) * t47 + pkin(8) * (t367 * t56 + t370 * t48)) * t362, t514, t591, t593, t534, -t592, t535, (t28 + t567) * t364 + (t367 * (-t366 * t55 + t369 * t61 - t575) + t370 * (-t512 - t577) + t565) * t362, (t25 + t568) * t364 + (t367 * (-t100 * t366 + t33 * t369 - t574) + t370 * (-t382 - t576) + t566) * t362, (t27 + t601) * t364 + (t367 * (-t366 * t54 + t369 * t60 - t597) + t370 * (-t518 - t598) - t600) * t362, (t7 + pkin(1) * (t32 * t367 - t370 * t42)) * t364 + (-pkin(1) * t31 + (pkin(8) * t42 - pkin(9) * t31 - t26 * t366 + t29 * t369) * t367 + (-pkin(2) * t31 + pkin(8) * t32 - t519) * t370) * t362, t364 * t71 + (t367 * (t369 * t84 - t436) - t370 * t82) * t362, t364 * t44 + (t367 * (-t198 * t366 + t369 * t52) - t370 * t50) * t362, t364 * t65 + (t367 * (t133 * t366 + t369 * t93) - t370 * t91) * t362, t364 * t70 + (t367 * (t369 * t83 + t436) - t370 * t81) * t362, t364 * t66 + (t367 * (-t366 * t391 + t369 * t94) - t370 * t92) * t362, t364 * t103 + (t367 * (t106 * t369 - t290 * t366) - t370 * t105) * t362, (t6 + pkin(1) * (t367 * t59 - t370 * t73)) * t364 + (t367 * (-pkin(9) * t58 + t18 * t369 - t21 * t366) + t370 * (-pkin(2) * t58 - t385) - pkin(1) * t58 + pkin(8) * (t367 * t73 + t370 * t59)) * t362, (t8 + pkin(1) * (t367 * t63 - t370 * t85)) * t364 + (t367 * (-pkin(9) * t62 + t22 * t369 - t23 * t366) + t370 * (-pkin(2) * t62 - t384) - pkin(1) * t62 + pkin(8) * (t367 * t85 + t370 * t63)) * t362, (t4 + pkin(1) * (t367 * t41 - t370 * t51)) * t364 + (t367 * (-pkin(9) * t40 - t24 * t366 + t369 * t5) + t370 * (-pkin(2) * t40 - t383) - pkin(1) * t40 + pkin(8) * (t367 * t51 + t370 * t41)) * t362, (t1 + pkin(1) * (t10 * t367 - t11 * t370)) * t364 + (t367 * (-pkin(9) * t9 - t2 * t366 + t3 * t369) + t370 * (-pkin(2) * t9 - t386) - pkin(1) * t9 + pkin(8) * (t10 * t370 + t367 * t11)) * t362; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t348, t341, t313, t348, t314, t355, -t295, -t296, 0, 0, t248, t195, t238, t247, t239, t274, t145, t147, t95, t104, t412, -t119, t579, t513, t140, t517, t36, t38, t30, t16, t412, t579, t119, t517, -t140, t513, t28, t25, t27, t7, t71, t44, t65, t70, t66, t103, t6, t8, t4, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t460, t298, t252, -t460, -t249, t380, -t190, -t191, 0, 0, t206, t151, t571, t203, t184, -t527, t393, t392, t390, t410, t206, t571, -t151, -t527, -t184, t203, t512, t382, t518, t519, t82, t50, t91, t81, t92, t105, t385, t384, t383, t386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t398, t526, t230, t148, 0, 0, 0, 0, 0, 0, t398, t230, -t526, t96, 0, 0, 0, 0, 0, 0, -t129, -t402, -t163, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t420, t212, t524, t78, 0, 0, 0, 0, 0, 0, t107, t135, t88, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t472, t198, -t133, -t472, t391, t290, -t34, -t35, 0, 0;];
tauJ_reg  = t64;
