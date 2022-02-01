% Calculate vector of inverse dynamics joint torques for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:17
% EndTime: 2022-01-23 09:16:42
% DurationCPUTime: 18.18s
% Computational Cost: add. (21430->797), mult. (23803->1062), div. (0->0), fcn. (22780->10), ass. (0->367)
t374 = sin(qJ(1));
t371 = cos(pkin(9));
t346 = t371 * pkin(3) + pkin(2);
t367 = pkin(9) + qJ(4);
t354 = cos(t367);
t581 = pkin(4) * t354;
t303 = t346 + t581;
t372 = cos(pkin(8));
t520 = t303 * t372;
t441 = pkin(1) + t520;
t353 = sin(t367);
t549 = pkin(4) * t353;
t369 = sin(pkin(9));
t550 = pkin(3) * t369;
t307 = t549 + t550;
t373 = qJ(3) + pkin(6);
t365 = -pkin(7) - t373;
t375 = cos(qJ(1));
t370 = sin(pkin(8));
t518 = t370 * t374;
t575 = t375 * t307 + t365 * t518;
t361 = t375 * qJ(2);
t288 = t346 * t372 + t370 * t373 + pkin(1);
t343 = qJ(2) + t550;
t572 = -t288 * t374 + t343 * t375;
t588 = -t361 + t572;
t110 = t374 * t441 - t575 + t588;
t356 = qJ(5) + t367;
t344 = sin(t356);
t345 = cos(t356);
t516 = t372 * t374;
t246 = t344 * t516 + t345 * t375;
t514 = t375 * t344;
t247 = t345 * t516 - t514;
t501 = -t247 * rSges(6,1) + t246 * rSges(6,2);
t146 = rSges(6,3) * t518 - t501;
t201 = (t365 + t373) * t372 + (t303 - t346) * t370;
t223 = -rSges(6,3) * t372 + (rSges(6,1) * t345 - rSges(6,2) * t344) * t370;
t368 = qJD(4) + qJD(5);
t432 = t368 * t370;
t295 = t374 * t432;
t304 = -t368 * t372 + qJD(1);
t478 = qJD(4) * t372;
t342 = qJD(1) - t478;
t477 = qJD(4) * t374;
t449 = t370 * t477;
t608 = -t110 * t342 - t146 * t304 + t201 * t449 + t223 * t295;
t314 = pkin(2) * t372 + qJ(3) * t370 + pkin(1);
t300 = t314 * t375;
t359 = t374 * qJ(2);
t497 = t300 + t359;
t590 = t288 * t375 + t343 * t374;
t172 = -t497 + t590;
t364 = t375 * pkin(1);
t279 = -t364 + t300;
t321 = t364 + t359;
t310 = qJD(1) * t321;
t480 = qJD(3) * t374;
t334 = t370 * t480;
t358 = qJD(2) * t375;
t487 = t358 - t334;
t428 = -qJD(1) * t279 - t310 + t487;
t607 = qJD(1) * t172 - t428;
t270 = t353 * t516 + t354 * t375;
t513 = t375 * t353;
t271 = t354 * t516 - t513;
t159 = Icges(5,5) * t271 - Icges(5,6) * t270 + Icges(5,3) * t518;
t252 = Icges(5,4) * t271;
t162 = -Icges(5,2) * t270 + Icges(5,6) * t518 + t252;
t251 = Icges(5,4) * t270;
t166 = -Icges(5,1) * t271 - Icges(5,5) * t518 + t251;
t66 = -t159 * t372 - t370 * (t162 * t353 + t166 * t354);
t137 = Icges(6,5) * t247 - Icges(6,6) * t246 + Icges(6,3) * t518;
t229 = Icges(6,4) * t247;
t140 = -Icges(6,2) * t246 + Icges(6,6) * t518 + t229;
t228 = Icges(6,4) * t246;
t144 = -Icges(6,1) * t247 - Icges(6,5) * t518 + t228;
t59 = -t137 * t372 - t370 * (t140 * t344 + t144 * t345);
t299 = t314 * t374;
t604 = t299 + t588;
t248 = t345 * t374 - t372 * t514;
t515 = t372 * t375;
t249 = t344 * t374 + t345 * t515;
t517 = t370 * t375;
t148 = t249 * rSges(6,1) + t248 * rSges(6,2) + rSges(6,3) * t517;
t296 = t375 * t432;
t297 = t374 * t307;
t392 = t303 * t515 - t365 * t517 + t297 + t321;
t574 = t392 - t590;
t602 = t148 * t304 - t223 * t296 + t574 * t342 + t607;
t272 = t354 * t374 - t372 * t513;
t273 = t353 * t374 + t354 * t515;
t170 = t273 * rSges(5,1) + t272 * rSges(5,2) + rSges(5,3) * t517;
t601 = t170 * t342 + t607;
t49 = t137 * t517 + t248 * t140 - t144 * t249;
t139 = Icges(6,5) * t249 + Icges(6,6) * t248 + Icges(6,3) * t517;
t523 = Icges(6,4) * t249;
t142 = Icges(6,2) * t248 + Icges(6,6) * t517 + t523;
t230 = Icges(6,4) * t248;
t145 = Icges(6,1) * t249 + Icges(6,5) * t517 + t230;
t50 = t139 * t517 + t248 * t142 + t249 * t145;
t220 = -Icges(6,3) * t372 + (Icges(6,5) * t345 - Icges(6,6) * t344) * t370;
t521 = Icges(6,4) * t345;
t221 = -Icges(6,6) * t372 + (-Icges(6,2) * t344 + t521) * t370;
t522 = Icges(6,4) * t344;
t222 = -Icges(6,5) * t372 + (Icges(6,1) * t345 - t522) * t370;
t72 = t220 * t517 + t221 * t248 + t222 * t249;
t22 = t295 * t49 + t296 * t50 + t72 * t304;
t168 = t271 * rSges(5,1) - t270 * rSges(5,2) + rSges(5,3) * t518;
t227 = -rSges(5,3) * t372 + (rSges(5,1) * t354 - rSges(5,2) * t353) * t370;
t600 = -t168 * t342 + t227 * t449;
t598 = -t162 * t270 - t166 * t271;
t597 = t272 * t162 - t166 * t273;
t47 = t137 * t518 - t140 * t246 - t144 * t247;
t591 = t159 * t517;
t589 = (t369 * t516 + t371 * t375) * rSges(4,2) - (-t369 * t375 + t371 * t516) * rSges(4,1);
t551 = pkin(1) * t374;
t278 = -t299 + t551;
t335 = qJD(3) * t517;
t319 = -t361 + t551;
t357 = qJD(2) * t374;
t488 = -qJD(1) * t319 + t357;
t429 = qJD(1) * t278 + t335 + t488;
t492 = t335 + t357;
t587 = -qJD(1) * t604 - t429 + t492;
t586 = -(t369 * t374 + t371 * t515) * rSges(4,1) - (-t369 * t515 + t371 * t374) * rSges(4,2);
t56 = t591 + t597;
t578 = t56 - t591;
t181 = t248 * rSges(6,1) - t249 * rSges(6,2);
t415 = -rSges(6,1) * t344 - rSges(6,2) * t345;
t250 = t415 * t370;
t483 = qJD(1) * t374;
t451 = t370 * t483;
t576 = -t304 * t181 + t223 * t451 + t250 * t296;
t180 = -t246 * rSges(6,1) - t247 * rSges(6,2);
t217 = t368 * t250;
t482 = qJD(1) * t375;
t450 = t370 * t482;
t154 = qJD(1) * t248 - t247 * t368;
t155 = qJD(1) * t249 - t246 * t368;
t417 = rSges(6,1) * t155 + rSges(6,2) * t154;
t86 = rSges(6,3) * t450 + t417;
t571 = t180 * t304 + t217 * t518 + t223 * t450 - t295 * t250 + t372 * t86;
t276 = rSges(3,1) * t515 - rSges(3,2) * t517 + t374 * rSges(3,3);
t569 = qJD(1) * t276;
t48 = t139 * t518 - t246 * t142 + t247 * t145;
t161 = Icges(5,5) * t273 + Icges(5,6) * t272 + Icges(5,3) * t517;
t526 = Icges(5,4) * t273;
t164 = Icges(5,2) * t272 + Icges(5,6) * t517 + t526;
t253 = Icges(5,4) * t272;
t167 = Icges(5,1) * t273 + Icges(5,5) * t517 + t253;
t55 = t161 * t518 - t270 * t164 + t271 * t167;
t197 = rSges(4,3) * t518 - t589;
t499 = t278 - t319;
t455 = t499 + t604;
t393 = qJD(1) * t455 + t492;
t64 = t393 + t600;
t533 = t375 * t64;
t476 = qJD(4) * t375;
t448 = t370 * t476;
t65 = -t227 * t448 + t601;
t568 = (t374 * t65 + t533) * qJD(1);
t566 = t374 * (-Icges(5,2) * t271 - t166 - t251) + t375 * (-Icges(5,2) * t273 + t167 + t253);
t244 = (-Icges(6,2) * t345 - t522) * t370;
t565 = t295 * (-Icges(6,2) * t247 - t144 - t228) + t296 * (-Icges(6,2) * t249 + t145 + t230) + t304 * (t222 + t244);
t564 = t372 ^ 2;
t479 = qJD(4) * t370;
t444 = qJD(1) * t479;
t470 = qJDD(4) * t370;
t289 = t374 * t470 + t375 * t444;
t208 = (qJD(5) * t482 + qJDD(5) * t374) * t370 + t289;
t563 = t208 / 0.2e1;
t331 = t375 * t470;
t209 = t331 + (qJDD(5) * t375 - t368 * t483) * t370;
t562 = t209 / 0.2e1;
t561 = t289 / 0.2e1;
t290 = -t374 * t444 + t331;
t560 = t290 / 0.2e1;
t559 = -t295 / 0.2e1;
t558 = t295 / 0.2e1;
t556 = t296 / 0.2e1;
t554 = -t372 / 0.2e1;
t553 = t374 / 0.2e1;
t552 = -t375 / 0.2e1;
t548 = g(1) * t374;
t302 = qJDD(1) + (-qJDD(4) - qJDD(5)) * t372;
t243 = (-Icges(6,5) * t344 - Icges(6,6) * t345) * t370;
t214 = t368 * t243;
t215 = t368 * t244;
t245 = (-Icges(6,1) * t344 - t521) * t370;
t216 = t368 * t245;
t63 = -t214 * t372 + ((-t221 * t368 + t216) * t345 + (-t222 * t368 - t215) * t344) * t370;
t91 = -t220 * t372 + (-t221 * t344 + t222 * t345) * t370;
t547 = t91 * t302 + t63 * t304;
t224 = -Icges(5,3) * t372 + (Icges(5,5) * t354 - Icges(5,6) * t353) * t370;
t524 = Icges(5,4) * t354;
t225 = -Icges(5,6) * t372 + (-Icges(5,2) * t353 + t524) * t370;
t525 = Icges(5,4) * t353;
t226 = -Icges(5,5) * t372 + (Icges(5,1) * t354 - t525) * t370;
t102 = -t224 * t372 + (-t225 * t353 + t226 * t354) * t370;
t336 = -qJDD(4) * t372 + qJDD(1);
t267 = (-Icges(5,5) * t353 - Icges(5,6) * t354) * t370;
t239 = qJD(4) * t267;
t268 = (-Icges(5,2) * t354 - t525) * t370;
t240 = qJD(4) * t268;
t269 = (-Icges(5,1) * t353 - t524) * t370;
t241 = qJD(4) * t269;
t69 = -t239 * t372 + (-t240 * t353 + t241 * t354 + (-t225 * t354 - t226 * t353) * qJD(4)) * t370;
t546 = t102 * t336 + t69 * t342;
t152 = qJD(1) * t246 - t249 * t368;
t153 = -qJD(1) * t247 + t248 * t368;
t509 = t153 * rSges(6,1) + t152 * rSges(6,2);
t85 = -rSges(6,3) * t451 + t509;
t312 = t343 * t482;
t350 = qJ(2) * t482;
t463 = pkin(4) * t476;
t382 = -t353 * t372 * t463 + t307 * t482 + t365 * t451 + t477 * t581 + t350;
t434 = t288 - t520;
t88 = -t312 + (-pkin(1) + t434) * t483 + t382;
t545 = -t85 - t88;
t544 = rSges(3,1) * t372;
t80 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t450;
t82 = Icges(6,4) * t155 + Icges(6,2) * t154 + Icges(6,6) * t450;
t84 = Icges(6,1) * t155 + Icges(6,4) * t154 + Icges(6,5) * t450;
t29 = -t372 * t80 + ((-t140 * t368 + t84) * t345 + (t144 * t368 - t82) * t344) * t370;
t543 = t29 * t295;
t79 = Icges(6,5) * t153 + Icges(6,6) * t152 - Icges(6,3) * t451;
t81 = Icges(6,4) * t153 + Icges(6,2) * t152 - Icges(6,6) * t451;
t83 = Icges(6,1) * t153 + Icges(6,4) * t152 - Icges(6,5) * t451;
t30 = -t372 * t79 + ((-t142 * t368 + t83) * t345 + (-t145 * t368 - t81) * t344) * t370;
t540 = t30 * t296;
t187 = qJD(1) * t272 - qJD(4) * t271;
t390 = t270 * qJD(4);
t188 = qJD(1) * t273 - t390;
t418 = rSges(5,1) * t188 + rSges(5,2) * t187;
t100 = rSges(5,3) * t450 + t418;
t274 = (-rSges(5,1) * t353 - rSges(5,2) * t354) * t370;
t242 = qJD(4) * t274;
t472 = qJDD(3) * t370;
t474 = qJD(1) * qJD(2);
t490 = qJDD(2) * t374 + t375 * t474;
t453 = t375 * t472 + t490;
t313 = t343 * t483;
t349 = qJ(2) * t483;
t496 = t313 - t349;
t498 = t288 - t314;
t352 = pkin(1) * t482;
t486 = t358 - t349;
t291 = t352 - t486;
t505 = -t314 * t482 - t291 - t334 + t352;
t378 = t455 * qJDD(1) + (-t482 * t498 - t334 - t496 + t505) * qJD(1) + t453;
t32 = -t100 * t342 - t168 * t336 + t227 * t289 + t242 * t449 + t378;
t538 = t32 * t374;
t489 = t350 + t357;
t458 = qJD(1) * (-pkin(1) * t483 + t489) + qJDD(1) * t321 + t374 * t474;
t394 = qJDD(1) * t279 + t374 * t472 + t458 + (0.2e1 * t335 + (pkin(1) - t314) * t483) * qJD(1);
t385 = qJD(1) * (-t483 * t498 + t312 - t350) + qJDD(1) * t172 + t394;
t185 = qJD(1) * t270 - qJD(4) * t273;
t186 = -qJD(1) * t271 + qJD(4) * t272;
t506 = t186 * rSges(5,1) + t185 * rSges(5,2);
t99 = -rSges(5,3) * t451 + t506;
t33 = t170 * t336 - t227 * t290 + t342 * t99 + (-t242 * t479 - qJDD(2)) * t375 + t385;
t537 = t33 * t375;
t77 = t224 * t518 - t225 * t270 + t226 * t271;
t536 = t342 * t77;
t54 = t159 * t518 + t598;
t535 = t374 * t54;
t40 = t393 + t608;
t532 = t40 * t201;
t41 = -t201 * t448 + t602;
t531 = t41 * t217;
t530 = t59 * t208;
t60 = -t139 * t372 + (-t142 * t344 + t145 * t345) * t370;
t529 = t60 * t209;
t528 = t66 * t289;
t67 = -t161 * t372 + (-t164 * t353 + t167 * t354) * t370;
t527 = t67 * t290;
t512 = -t574 - t148;
t503 = -t225 + t269;
t502 = t226 + t268;
t500 = t589 * qJD(1);
t469 = rSges(3,1) * t516;
t491 = rSges(3,2) * t518 + t375 * rSges(3,3);
t275 = t469 - t491;
t494 = -t319 - t275;
t493 = rSges(3,2) * t451 + rSges(3,3) * t482;
t481 = qJD(3) * t372;
t475 = -m(4) - m(5) - m(6);
t473 = qJDD(2) * t375;
t471 = qJDD(3) * t372;
t57 = t161 * t517 + t272 * t164 + t273 * t167;
t198 = rSges(4,3) * t517 - t586;
t454 = -t197 + t499;
t452 = -t334 + t486;
t447 = t518 / 0.2e1;
t446 = t517 / 0.2e1;
t445 = -pkin(1) - t544;
t443 = -t479 / 0.2e1;
t442 = t479 / 0.2e1;
t438 = t296 * t180 - t181 * t295;
t433 = t370 ^ 2 * qJD(4) ^ 2 * t549;
t427 = -t451 / 0.2e1;
t426 = qJD(1) * t446;
t425 = t374 * t443;
t424 = t374 * t442;
t423 = t375 * t443;
t422 = t375 * t442;
t322 = rSges(2,1) * t375 - rSges(2,2) * t374;
t320 = rSges(2,1) * t374 + rSges(2,2) * t375;
t419 = t586 * qJD(1);
t93 = Icges(5,5) * t186 + Icges(5,6) * t185 - Icges(5,3) * t451;
t94 = Icges(5,5) * t188 + Icges(5,6) * t187 + Icges(5,3) * t450;
t95 = Icges(5,4) * t186 + Icges(5,2) * t185 - Icges(5,6) * t451;
t96 = Icges(5,4) * t188 + Icges(5,2) * t187 + Icges(5,6) * t450;
t97 = Icges(5,1) * t186 + Icges(5,4) * t185 - Icges(5,5) * t451;
t98 = Icges(5,1) * t188 + Icges(5,4) * t187 + Icges(5,5) * t450;
t414 = (t162 * t185 - t166 * t186 + t272 * t96 + t273 * t98 + (-t159 * t483 + t375 * t94) * t370) * t374 + (t164 * t185 + t167 * t186 + t272 * t95 + t273 * t97 + (-t161 * t483 + t375 * t93) * t370) * t375;
t413 = (t162 * t187 - t166 * t188 - t270 * t96 + t271 * t98 + (t159 * t482 + t374 * t94) * t370) * t374 + (t164 * t187 + t167 * t188 - t270 * t95 + t271 * t97 + (t161 * t482 + t374 * t93) * t370) * t375;
t34 = -t372 * t94 + (-t353 * t96 + t354 * t98 + (-t162 * t354 + t166 * t353) * qJD(4)) * t370;
t35 = -t372 * t93 + (-t353 * t95 + t354 * t97 + (-t164 * t354 - t167 * t353) * qJD(4)) * t370;
t412 = t34 * t374 + t35 * t375;
t411 = t374 * t64 - t375 * t65;
t409 = t100 * t375 - t374 * t99;
t408 = t168 * t375 - t170 * t374;
t407 = (-Icges(5,5) * t270 - Icges(5,6) * t271) * t374 + (Icges(5,5) * t272 - Icges(5,6) * t273) * t375;
t255 = t272 * pkin(4);
t398 = (t375 * t55 + t535) * t370;
t397 = (t374 * t56 + t375 * t57) * t370;
t395 = t370 * t407;
t254 = t270 * pkin(4);
t391 = (-Icges(6,5) * t246 - Icges(6,6) * t247) * t295 + (Icges(6,5) * t248 - Icges(6,6) * t249) * t296 + t243 * t304;
t389 = (Icges(5,1) * t272 - t164 - t526) * t375 + (-Icges(5,1) * t270 - t162 - t252) * t374;
t89 = t352 - pkin(4) * t390 + (t297 + (-t365 * t370 - t434) * t375) * qJD(1) - t496;
t11 = -t336 * t110 - t302 * t146 + t289 * t201 + t208 * t223 + t295 * t217 - t304 * t86 - t342 * t89 - t374 * t433 + t378;
t44 = -t481 + t146 * t296 - t148 * t295 + (t110 * t375 - t374 * t574) * t479;
t9 = -t471 + t110 * t290 - t574 * t289 + t146 * t209 - t148 * t208 - t295 * t85 + t296 * t86 + (-t374 * t88 + t375 * t89) * t479;
t387 = t11 * (t372 * t146 + t223 * t518) + (t146 * t9 + t44 * t86) * t517;
t71 = t220 * t518 - t221 * t246 + t222 * t247;
t383 = t370 * t391;
t13 = t140 * t152 - t144 * t153 + t248 * t82 + t249 * t84 + (-t137 * t483 + t375 * t80) * t370;
t14 = t142 * t152 + t145 * t153 + t248 * t81 + t249 * t83 + (-t139 * t483 + t375 * t79) * t370;
t15 = t140 * t154 - t144 * t155 - t246 * t82 + t247 * t84 + (t137 * t482 + t374 * t80) * t370;
t16 = t142 * t154 + t145 * t155 - t246 * t81 + t247 * t83 + (t139 * t482 + t374 * t79) * t370;
t380 = (Icges(6,1) * t248 - t142 - t523) * t296 + (-Icges(6,1) * t246 - t140 - t229) * t295 + (-t221 + t245) * t304;
t42 = t152 * t221 + t153 * t222 + t215 * t248 + t216 * t249 + (t214 * t375 - t220 * t483) * t370;
t43 = t154 * t221 + t155 * t222 - t215 * t246 + t216 * t247 + (t214 * t374 + t220 * t482) * t370;
t381 = (t15 * t295 + t16 * t296 + t208 * t47 + t209 * t48 + t302 * t71 + t304 * t43) * t447 - (-t391 * t372 + (-t344 * t565 + t345 * t380) * t370) * t304 / 0.2e1 + t22 * t427 + (t295 * t47 + t296 * t48 + t304 * t71) * t426 + (-t372 * t71 + (t374 * t47 + t375 * t48) * t370) * t563 + (t13 * t295 + t14 * t296 + t208 * t49 + t209 * t50 + t302 * t72 + t304 * t42) * t446 + (-t372 * t72 + (t374 * t49 + t375 * t50) * t370) * t562 + (-t372 * t43 + (t15 * t374 + t16 * t375 + (-t374 * t48 + t375 * t47) * qJD(1)) * t370) * t558 + t302 * (-t372 * t91 + (t374 * t59 + t375 * t60) * t370) / 0.2e1 + (-t372 * t42 + (t13 * t374 + t14 * t375 + (-t374 * t50 + t375 * t49) * qJD(1)) * t370) * t556 + (t529 + t530 + t540 + t543 + t547) * t554 + t304 * (-t372 * t63 + (t29 * t374 + t30 * t375 + (-t374 * t60 + t375 * t59) * qJD(1)) * t370) / 0.2e1 + (-t246 * t565 + t247 * t380 + t374 * t383) * t559 - (t248 * t565 + t380 * t249 + t375 * t383) * t296 / 0.2e1;
t204 = t310 - t358 + t569;
t203 = qJD(1) * t494 + t357;
t196 = rSges(5,1) * t272 - rSges(5,2) * t273;
t195 = -rSges(5,1) * t270 - rSges(5,2) * t271;
t108 = -t473 + qJDD(1) * t276 + qJD(1) * (-qJD(1) * t469 + t493) + t458;
t107 = t494 * qJDD(1) + (-t291 - t569) * qJD(1) + t490;
t106 = qJD(1) * t198 - t428;
t105 = qJD(1) * t454 + t492;
t78 = t224 * t517 + t225 * t272 + t226 * t273;
t76 = t408 * t479 - t481;
t73 = t78 * t342;
t52 = -t473 + qJDD(1) * t198 + qJD(1) * (-rSges(4,3) * t451 + t500) + t394;
t51 = t454 * qJDD(1) + ((-rSges(4,3) * t482 - t480) * t370 + t419 + t505) * qJD(1) + t453;
t46 = t187 * t225 + t188 * t226 - t240 * t270 + t241 * t271 + (t224 * t482 + t239 * t374) * t370;
t45 = t185 * t225 + t186 * t226 + t240 * t272 + t241 * t273 + (-t224 * t483 + t239 * t375) * t370;
t36 = t168 * t290 - t170 * t289 + t409 * t479 - t471;
t28 = qJD(4) * t397 + t73;
t27 = qJD(4) * t398 + t536;
t12 = t385 + (-qJDD(2) + t433) * t375 + t574 * t336 + t148 * t302 - t201 * t290 - t209 * t223 - t217 * t296 + t304 * t85 + t342 * t88;
t1 = [(-t536 + ((-t55 + t578 - t597) * t375 - t535) * t479 + t27) * t423 + (-t203 * t291 + t204 * (t489 + t493) + (t203 * (rSges(3,2) * t370 - t544) * t375 + (-t203 * rSges(3,3) + t204 * t445) * t374) * qJD(1) - (-qJD(1) * t275 - t203 + t488) * t204 + (t108 - g(2)) * (t276 + t321) + (t107 - g(1)) * (t374 * t445 + t361 + t491)) * m(3) + (m(2) * (t320 ^ 2 + t322 ^ 2) + Icges(2,3) + (Icges(3,2) + Icges(4,3)) * t564 + ((Icges(3,1) + Icges(4,1) * t371 ^ 2 + (-0.2e1 * Icges(4,4) * t371 + Icges(4,2) * t369) * t369) * t370 + 0.2e1 * (-Icges(4,5) * t371 + Icges(4,6) * t369 + Icges(3,4)) * t372) * t370) * qJDD(1) + (-(-qJD(1) * t197 - t105 + t429) * t106 + t105 * (t419 + t452) + t106 * (t335 + t489 + t500) + (t105 * t375 + t106 * t374) * qJD(1) * (-rSges(4,3) * t370 - t314) + (t52 - g(2)) * (t198 + t497) + (t51 - g(1)) * (-t197 + t361 - t299)) * m(4) + (-t227 * t533 * t479 + (-rSges(5,3) * t370 - t288) * t568 + (t312 + t506 + t587 - t600) * t65 + (-t313 - t418 + t487 + t601) * t64 + (t33 - g(2)) * (t170 + t590) + (t32 - g(1)) * (-t168 + t572)) * m(5) + t543 / 0.2e1 + t540 / 0.2e1 + (t35 + t45) * t422 + (t43 + t22) * t558 + (t34 + t46 + t28) * t424 + (-t532 * t375 * t479 + (-g(2) + t12) * (t392 + t148) + (t382 + t509 + t587 - t608) * t41 + (-t548 + (qJD(1) * t41 + t11) * t374) * (-rSges(6,3) * t370 - t441) + (t354 * t463 - t352 - t417 + t452 + t478 * t549 * t374 + ((-t520 + (-rSges(6,3) + t365) * t370) * t375 - t297) * qJD(1) + t602) * t40 + (-g(1) + t11) * (t361 + t575 + t501)) * m(6) + t22 * t559 - m(2) * (-g(1) * t320 + g(2) * t322) + (t73 + ((t54 + t57 - t598) * t375 + t578 * t374) * t479) * t425 + t529 / 0.2e1 + t530 / 0.2e1 + t527 / 0.2e1 + t528 / 0.2e1 + t42 * t556 + t78 * t560 + t77 * t561 + t72 * t562 + t71 * t563 + t546 + t547; (-m(3) + t475) * (-g(2) * t375 + t548) + 0.2e1 * (t11 * t553 + t12 * t552) * m(6) + 0.2e1 * (t538 / 0.2e1 - t537 / 0.2e1) * m(5) + 0.2e1 * (t51 * t553 + t52 * t552) * m(4) + 0.2e1 * (t107 * t553 + t108 * t552) * m(3); t475 * (-g(3) * t372 + (g(1) * t375 + g(2) * t374) * t370) + m(4) * (qJDD(3) * t564 + t51 * t517 + t518 * t52) + m(5) * (t32 * t517 + t33 * t518 - t36 * t372) + m(6) * (t11 * t517 + t12 * t518 - t372 * t9); (t289 * t54 + t290 * t55 + t336 * t77 + t342 * t46 + t413 * t479) * t447 + (-t372 * t77 + t398) * t561 + (t412 * t479 + t527 + t528 + t546) * t554 + t381 + t27 * t426 + t28 * t427 - t342 * (-t372 * t267 * t342 + ((-t353 * t502 + t354 * t503) * t342 + ((-t353 * t566 + t354 * t389) * t370 - t407 * t372) * qJD(4)) * t370) / 0.2e1 + t336 * (-t102 * t372 + (t374 * t66 + t375 * t67) * t370) / 0.2e1 + ((t267 * t517 + t272 * t502 + t273 * t503) * t342 + (t272 * t566 + t389 * t273 + t375 * t395) * t479) * t423 + ((t267 * t518 - t270 * t502 + t271 * t503) * t342 + (-t270 * t566 + t271 * t389 + t374 * t395) * t479) * t425 + (-t372 * t46 + ((-t55 * t374 + t54 * t375) * qJD(1) + t413) * t370) * t424 + (-t372 * t45 + ((-t57 * t374 + t56 * t375) * qJD(1) + t414) * t370) * t422 + (t289 * t56 + t290 * t57 + t336 * t78 + t342 * t45 + t414 * t479) * t446 + t342 * (-t372 * t69 + ((-t374 * t67 + t66 * t375) * qJD(1) + t412) * t370) / 0.2e1 + (-t372 * t78 + t397) * t560 + (-t44 * ((-t254 * t375 - t255 * t374) * t479 + t438) + (t11 * t110 + t12 * t512) * t372 + t387 - g(1) * (t255 + t181) - g(2) * (-t254 + t180) + (-t255 * t342 + t372 * t545 + t576) * t41 + ((t12 * (-t201 - t223) - t531 + t9 * t110 + t44 * t89 + (t44 * t512 + t532) * qJD(1)) * t375 + (t11 * t201 + t9 * t512 + t44 * t545 + (t41 * t201 + t44 * (-t110 - t146)) * qJD(1)) * t374 - g(3) * (t415 - t549)) * t370 + (-t254 * t342 + t372 * t89 + t571) * t40) * m(6) + (-(-t195 * t64 + t196 * t65) * t342 - (t76 * (t195 * t375 - t196 * t374) + t411 * t274) * t479 + (t100 * t64 + t168 * t32 - t170 * t33 - t65 * t99) * t372 + (t36 * t408 + t76 * (-t168 * t483 - t170 * t482 + t409) + t411 * t242 + (t538 - t537 + t568) * t227) * t370 - g(1) * t196 - g(2) * t195 - g(3) * t274) * m(5); t381 + (-t12 * t148 * t372 + ((-qJD(1) * t148 * t44 - t12 * t223 - t531) * t375 + (-t9 * t148 + t44 * (-qJD(1) * t146 - t85)) * t374) * t370 + t387 - t44 * t438 - g(1) * t181 - g(2) * t180 - g(3) * t250 + (-t372 * t85 + t576) * t41 + t571 * t40) * m(6);];
tau = t1;
