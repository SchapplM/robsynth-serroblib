% Calculate vector of inverse dynamics joint torques for
% S5PRPPR3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:18
% EndTime: 2019-12-05 15:26:53
% DurationCPUTime: 28.04s
% Computational Cost: add. (9248->627), mult. (14022->911), div. (0->0), fcn. (12830->8), ass. (0->324)
t572 = Icges(5,4) - Icges(4,5);
t571 = Icges(5,5) - Icges(4,6);
t288 = qJ(2) + pkin(8);
t283 = sin(t288);
t284 = cos(t288);
t293 = sin(qJ(2));
t295 = cos(qJ(2));
t378 = -Icges(3,5) * t293 - Icges(3,6) * t295;
t566 = t572 * t283 + t571 * t284 + t378;
t570 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t569 = t566 * qJD(2);
t568 = Icges(3,5) * t295 - Icges(3,6) * t293 + t571 * t283 - t572 * t284;
t292 = sin(qJ(5));
t294 = cos(qJ(5));
t405 = rSges(6,1) * t292 + rSges(6,2) * t294;
t549 = t284 * t405;
t289 = sin(pkin(7));
t290 = cos(pkin(7));
t567 = t289 * t290;
t286 = t289 ^ 2;
t287 = t290 ^ 2;
t528 = t286 + t287;
t269 = rSges(3,1) * t293 + rSges(3,2) * t295;
t565 = t269 * t528;
t133 = t284 * rSges(6,3) + t283 * t405;
t564 = t568 * t289 - t570 * t290;
t563 = t570 * t289 + t568 * t290;
t562 = t569 * t289;
t561 = t569 * t290;
t483 = Icges(5,6) * t284;
t373 = Icges(5,3) * t283 - t483;
t144 = -Icges(5,5) * t290 + t289 * t373;
t479 = t283 * t289;
t263 = Icges(5,6) * t479;
t477 = t284 * t289;
t146 = -Icges(5,4) * t290 - Icges(5,2) * t477 + t263;
t491 = Icges(3,4) * t295;
t386 = -Icges(3,2) * t293 + t491;
t180 = -Icges(3,6) * t290 + t289 * t386;
t492 = Icges(3,4) * t293;
t391 = Icges(3,1) * t295 - t492;
t182 = -Icges(3,5) * t290 + t289 * t391;
t385 = -Icges(3,2) * t295 - t492;
t390 = -Icges(3,1) * t293 - t491;
t457 = qJD(2) * t284;
t458 = qJD(2) * t283;
t489 = Icges(4,4) * t284;
t490 = Icges(4,4) * t283;
t384 = -Icges(4,2) * t283 + t489;
t139 = -Icges(4,6) * t290 + t289 * t384;
t389 = Icges(4,1) * t284 - t490;
t141 = -Icges(4,5) * t290 + t289 * t389;
t543 = t139 * t284 + t141 * t283;
t559 = (-(-t293 * t385 + t295 * t390) * qJD(2) + (-Icges(5,6) * t283 - t490 + (-Icges(4,2) - Icges(5,3)) * t284) * t458 + (t483 + t489 + (Icges(4,1) + Icges(5,2)) * t283) * t457) * t289 + (-t144 * t284 - t146 * t283 + t180 * t295 + t182 * t293 + t543) * qJD(2);
t557 = t180 * t293 - t182 * t295 + (-t141 + t146) * t284 + (t139 - t144) * t283;
t556 = 0.2e1 * qJD(2);
t555 = 2 * qJDD(2);
t553 = t528 * t458;
t476 = t284 * t290;
t478 = t283 * t290;
t552 = ((Icges(5,5) * t289 + t290 * t373) * t289 - t144 * t290) * t284 + ((Icges(5,4) * t289 + 0.2e1 * Icges(5,6) * t478 + (-Icges(5,2) + Icges(5,3)) * t476) * t289 - (Icges(5,3) * t477 + t146 + t263) * t290) * t283 + t289 * (-(Icges(4,6) * t289 + t290 * t384) * t284 - (Icges(4,5) * t289 + t290 * t389) * t283) + t290 * t543;
t451 = qJD(4) * t289;
t256 = t283 * t451;
t512 = pkin(3) * t283;
t240 = -qJ(4) * t284 + t512;
t347 = qJD(2) * t240;
t116 = -t289 * t347 + t256;
t450 = qJD(4) * t290;
t258 = t283 * t450;
t117 = -t290 * t347 + t258;
t261 = qJ(4) * t477;
t262 = qJ(4) * t476;
t455 = qJD(2) * t290;
t456 = qJD(2) * t289;
t497 = pkin(2) * qJD(2);
t436 = t293 * t497;
t453 = qJD(3) * t290;
t238 = -t289 * t436 - t453;
t282 = qJD(3) * t289;
t239 = -t290 * t436 + t282;
t465 = t289 * t238 + t290 * t239;
t551 = t289 * t116 + t290 * t117 - (-pkin(3) * t479 + t261) * t456 - (-pkin(3) * t478 + t262) * t455 - qJD(4) * t283 + t465;
t445 = qJD(2) * qJD(4);
t550 = qJDD(4) * t283 + t284 * t445;
t545 = t566 * t289;
t544 = t566 * t290;
t403 = rSges(5,2) * t283 + rSges(5,3) * t284;
t541 = t528 * qJD(2) * t403;
t285 = t295 * pkin(2);
t530 = t284 * rSges(4,1) - rSges(4,2) * t283;
t532 = t530 + t285;
t531 = -rSges(5,2) * t284 + t283 * rSges(5,3);
t529 = t284 * pkin(3) + t283 * qJ(4);
t431 = t285 + t529;
t527 = g(1) * t290 + g(2) * t289;
t375 = Icges(6,5) * t292 + Icges(6,6) * t294;
t314 = -Icges(6,3) * t283 + t284 * t375;
t486 = Icges(6,4) * t292;
t380 = Icges(6,2) * t294 + t486;
t315 = -Icges(6,6) * t283 + t284 * t380;
t485 = Icges(6,4) * t294;
t387 = Icges(6,1) * t292 + t485;
t316 = -Icges(6,5) * t283 + t284 * t387;
t526 = -t293 * (t385 * t289 + t182) - t295 * (-t390 * t289 + t180);
t472 = t290 * t294;
t475 = t289 * t292;
t217 = t283 * t472 - t475;
t199 = Icges(6,4) * t217;
t473 = t290 * t292;
t474 = t289 * t294;
t219 = t283 * t474 + t473;
t200 = Icges(6,4) * t219;
t204 = (Icges(6,2) * t292 - t485) * t284;
t218 = t283 * t473 + t474;
t220 = -t283 * t475 + t472;
t448 = qJD(5) * t284;
t232 = t290 * t448 + t456;
t233 = t289 * t448 - t455;
t449 = qJD(5) * t283;
t86 = Icges(6,1) * t218 + Icges(6,5) * t476 + t199;
t87 = -Icges(6,1) * t220 + Icges(6,5) * t477 + t200;
t303 = t232 * (-Icges(6,2) * t218 + t199 + t86) + t233 * (Icges(6,2) * t220 + t200 + t87) + t449 * (-t316 + t204);
t205 = (-Icges(6,1) * t294 + t486) * t284;
t487 = Icges(6,4) * t220;
t488 = Icges(6,4) * t218;
t84 = Icges(6,2) * t217 + Icges(6,6) * t476 + t488;
t85 = Icges(6,2) * t219 + Icges(6,6) * t477 - t487;
t304 = t232 * (-Icges(6,1) * t217 + t488 + t84) + t233 * (-Icges(6,1) * t219 - t487 + t85) + t449 * (-t315 - t205);
t296 = qJD(2) ^ 2;
t525 = -m(5) - m(6);
t444 = qJD(2) * qJD(5);
t333 = qJDD(5) * t284 - t283 * t444;
t443 = qJDD(2) * t289;
t124 = t290 * t333 + t443;
t524 = t124 / 0.2e1;
t442 = qJDD(2) * t290;
t125 = t289 * t333 - t442;
t523 = t125 / 0.2e1;
t522 = -t553 / 0.2e1;
t225 = qJDD(5) * t283 + t284 * t444;
t521 = t225 / 0.2e1;
t520 = -t232 / 0.2e1;
t519 = t232 / 0.2e1;
t518 = -t233 / 0.2e1;
t517 = t233 / 0.2e1;
t516 = -t284 / 0.2e1;
t513 = pkin(2) * t293;
t511 = pkin(6) * t283;
t510 = pkin(6) * t284;
t83 = -Icges(6,5) * t220 + Icges(6,6) * t219 + Icges(6,3) * t477;
t22 = t217 * t85 + t218 * t87 + t476 * t83;
t496 = t22 * t289;
t82 = Icges(6,5) * t218 + Icges(6,6) * t217 + Icges(6,3) * t476;
t23 = t219 * t84 - t220 * t86 + t477 * t82;
t495 = t23 * t290;
t46 = -t217 * t315 - t218 * t316 - t314 * t476;
t494 = t46 * t284;
t47 = -t219 * t315 + t220 * t316 - t314 * t477;
t493 = t47 * t284;
t480 = t314 * t283;
t471 = t295 * t296;
t135 = -qJ(3) * t290 + t285 * t289;
t136 = qJ(3) * t289 + t285 * t290;
t470 = t289 * t135 + t290 * t136;
t464 = t549 * t289;
t463 = t549 * t290;
t462 = t258 + t282;
t461 = rSges(5,2) * t479 + rSges(5,3) * t477;
t460 = rSges(5,2) * t478 + rSges(5,3) * t476;
t452 = qJD(4) * t284;
t446 = -m(4) + t525;
t441 = qJDD(3) * t290;
t435 = t295 * t497;
t281 = qJDD(3) * t289;
t432 = t290 * t550 + t281;
t430 = t283 * t456;
t429 = t283 * t455;
t428 = t292 * t457;
t427 = t294 * t457;
t426 = t135 * t456 + t136 * t455 + qJD(1);
t423 = -t456 / 0.2e1;
t421 = -t449 / 0.2e1;
t420 = t449 / 0.2e1;
t419 = -t240 - t513;
t242 = rSges(4,1) * t283 + rSges(4,2) * t284;
t323 = -t242 - t513;
t417 = t256 - t453;
t201 = t529 * t289;
t202 = t529 * t290;
t415 = t289 * t201 + t290 * t202 + t470;
t414 = t528 * t511;
t413 = t528 * t513;
t412 = t403 + t419;
t159 = qJD(2) * t529 - t452;
t410 = -t159 - t435;
t224 = t530 * qJD(2);
t409 = -t224 - t435;
t408 = -t512 - t513;
t407 = -t285 - t510;
t270 = rSges(3,1) * t295 - rSges(3,2) * t293;
t134 = rSges(6,3) * t283 - t549;
t358 = t419 - t511;
t300 = -qJD(2) * t159 + qJDD(2) * t358 + t296 * t407;
t393 = t289 * t550 - t441;
t110 = qJD(5) * t217 + t290 * t428;
t111 = -qJD(5) * t218 + t290 * t427;
t65 = rSges(6,1) * t110 + rSges(6,2) * t111 - rSges(6,3) * t429;
t208 = (-rSges(6,1) * t294 + rSges(6,2) * t292) * t284;
t81 = qJD(2) * t133 + qJD(5) * t208;
t88 = rSges(6,1) * t218 + rSges(6,2) * t217 + rSges(6,3) * t476;
t13 = -t124 * t134 + t225 * t88 - t232 * t81 + t289 * t300 + t449 * t65 + t393;
t112 = qJD(5) * t219 + t289 * t428;
t113 = qJD(5) * t220 + t289 * t427;
t66 = rSges(6,1) * t112 + rSges(6,2) * t113 - rSges(6,3) * t430;
t89 = -rSges(6,1) * t220 + rSges(6,2) * t219 + rSges(6,3) * t477;
t14 = t125 * t134 - t225 * t89 + t233 * t81 + t290 * t300 - t449 * t66 + t432;
t401 = -t13 * t290 + t14 * t289;
t21 = t217 * t84 + t218 * t86 + t476 * t82;
t400 = t21 * t290 + t496;
t24 = t219 * t85 - t220 * t87 + t477 * t83;
t399 = t24 * t289 + t495;
t398 = -t82 * t232 - t83 * t233;
t395 = t292 * t86 + t294 * t84;
t26 = t283 * t82 - t284 * t395;
t394 = t292 * t87 + t294 * t85;
t27 = t283 * t83 - t284 * t394;
t397 = t26 * t290 + t27 * t289;
t396 = t289 * t88 - t290 * t89;
t392 = qJD(2) * t323;
t372 = -t292 * t316 - t294 * t315;
t155 = -rSges(4,3) * t290 + t289 * t530;
t156 = rSges(4,3) * t289 + t290 * t530;
t367 = t155 * t289 + t156 * t290;
t157 = rSges(5,1) * t289 + t290 * t531;
t158 = -rSges(5,1) * t290 + t289 * t531;
t366 = t157 * t290 + t158 * t289;
t363 = t528 * t270;
t362 = qJD(2) * t565;
t236 = pkin(4) * t289 + pkin(6) * t476;
t237 = -pkin(4) * t290 + pkin(6) * t477;
t361 = t236 * t290 + t237 * t289;
t360 = t135 * t443 + t136 * t442 + t238 * t456 + t239 * t455 + qJDD(1);
t223 = t531 * qJD(2);
t359 = -t223 + t410;
t356 = qJD(2) * t412;
t59 = Icges(6,5) * t110 + Icges(6,6) * t111 - Icges(6,3) * t429;
t355 = t284 * t59 - t458 * t82;
t60 = Icges(6,5) * t112 + Icges(6,6) * t113 - Icges(6,3) * t430;
t354 = t284 * t60 - t458 * t83;
t126 = Icges(6,3) * t284 + t283 * t375;
t203 = (-Icges(6,5) * t294 + Icges(6,6) * t292) * t284;
t78 = qJD(2) * t126 + qJD(5) * t203;
t353 = t284 * t78 + t314 * t458;
t320 = t201 * t456 + t202 * t455 + t426 - t452;
t19 = qJD(2) * t361 + t232 * t89 - t233 * t88 + t320;
t352 = t19 * t396;
t48 = qJD(2) * t367 + t426;
t351 = t48 * t242;
t349 = qJD(2) * t242;
t340 = -t134 + t358;
t322 = qJD(2) * t358;
t321 = (t126 + t372) * t283;
t319 = t378 * t289 + ((-t386 + t390) * t295 + (-t385 - t391) * t293) * t290;
t318 = t203 * t449 + t232 * (Icges(6,5) * t217 - Icges(6,6) * t218) + t233 * (Icges(6,5) * t219 + Icges(6,6) * t220);
t317 = -pkin(6) * t457 + t410 - t81;
t130 = Icges(6,5) * t284 + t283 * t387;
t128 = Icges(6,6) * t284 + t283 * t380;
t313 = t284 * t318;
t306 = -qJD(2) * t224 - qJDD(2) * t242 + (-qJDD(2) * t293 - t471) * pkin(2);
t305 = -qJDD(4) * t284 + t116 * t456 + t117 * t455 + t201 * t443 + t202 * t442 + t283 * t445 + t360;
t302 = -pkin(2) * t471 + (-t159 - t223) * qJD(2) + t412 * qJDD(2);
t301 = (t314 * t290 + t395) * t232 + (t314 * t289 + t394) * t233;
t297 = (qJD(5) * t321 + t301) * t284;
t259 = t284 * t450;
t257 = t284 * t451;
t235 = t269 * t290;
t234 = t269 * t289;
t177 = t290 * t349;
t175 = t289 * t349;
t115 = t290 * t392 + t282;
t114 = t289 * t392 - t453;
t109 = -rSges(6,3) * t478 + t463;
t108 = -rSges(6,3) * t479 + t464;
t107 = t316 * t290;
t106 = t316 * t289;
t105 = t315 * t290;
t104 = t315 * t289;
t97 = rSges(6,1) * t219 + rSges(6,2) * t220;
t96 = rSges(6,1) * t217 - rSges(6,2) * t218;
t80 = qJD(2) * t130 + qJD(5) * t205;
t79 = qJD(2) * t128 + qJD(5) * t204;
t76 = t290 * t356 + t462;
t75 = t289 * t356 + t417;
t74 = t290 * t306 + t281;
t73 = t289 * t306 - t441;
t64 = Icges(6,1) * t112 + Icges(6,4) * t113 - Icges(6,5) * t430;
t63 = Icges(6,1) * t110 + Icges(6,4) * t111 - Icges(6,5) * t429;
t62 = Icges(6,4) * t112 + Icges(6,2) * t113 - Icges(6,6) * t430;
t61 = Icges(6,4) * t110 + Icges(6,2) * t111 - Icges(6,6) * t429;
t50 = -qJD(2) * t362 + qJDD(2) * t363 + qJDD(1);
t49 = -t284 * t372 - t480;
t41 = t290 * t302 + t432;
t40 = t289 * t302 + t393;
t39 = t134 * t233 + t290 * t322 - t449 * t89 + t462;
t38 = -t134 * t232 + t289 * t322 + t449 * t88 + t417;
t25 = qJD(2) * t366 + t320;
t20 = t367 * qJDD(2) + (-t175 * t289 - t177 * t290) * qJD(2) + t360;
t18 = (qJD(2) * t372 + t78) * t283 + (-qJD(2) * t314 - t292 * t80 - t294 * t79 + (-t292 * t315 + t294 * t316) * qJD(5)) * t284;
t17 = qJD(2) * t541 + t366 * qJDD(2) + t305;
t16 = -t112 * t316 - t113 * t315 + t219 * t79 - t220 * t80 + t289 * t353;
t15 = -t110 * t316 - t111 * t315 + t217 * t79 + t218 * t80 + t290 * t353;
t12 = t232 * t26 + t233 * t27 + t449 * t49;
t11 = t112 * t87 + t113 * t85 + t219 * t62 - t220 * t64 + t289 * t354;
t10 = t112 * t86 + t113 * t84 + t219 * t61 - t220 * t63 + t289 * t355;
t9 = t110 * t87 + t111 * t85 + t217 * t62 + t218 * t64 + t290 * t354;
t8 = t110 * t86 + t111 * t84 + t217 * t61 + t218 * t63 + t290 * t355;
t7 = t23 * t232 + t233 * t24 + t449 * t47;
t6 = t21 * t232 + t22 * t233 + t449 * t46;
t5 = (qJD(2) * t394 + t60) * t283 + (qJD(2) * t83 - t292 * t64 - t294 * t62 + (t292 * t85 - t294 * t87) * qJD(5)) * t284;
t4 = (qJD(2) * t395 + t59) * t283 + (qJD(2) * t82 - t292 * t63 - t294 * t61 + (t292 * t84 - t294 * t86) * qJD(5)) * t284;
t3 = qJDD(2) * t361 + t124 * t89 - t125 * t88 + t232 * t66 - t233 * t65 - t296 * t414 + t305;
t2 = t10 * t232 + t233 * t11 + t124 * t23 + t125 * t24 + t16 * t449 + t225 * t47;
t1 = t124 * t21 + t125 * t22 + t15 * t449 + t225 * t46 + t232 * t8 + t233 * t9;
t28 = [m(2) * qJDD(1) + (-m(2) - m(3) + t446) * g(3) + m(3) * t50 + m(4) * t20 + m(5) * t17 + m(6) * t3; -t12 * t448 / 0.2e1 + (t26 * t289 - t27 * t290) * t521 + (t23 * t289 - t24 * t290) * t523 + (t21 * t289 - t22 * t290) * t524 + (t10 * t289 - t11 * t290) * t517 + ((t105 * t219 - t107 * t220) * t232 + (t104 * t219 - t106 * t220) * t233 + (t493 + (t128 * t219 - t130 * t220 - t495) * t283) * qJD(5) + (((-t24 + t480) * qJD(5) + t398) * t283 + t297) * t289) * t518 + (t289 * t8 - t290 * t9) * t519 + ((t105 * t217 + t107 * t218) * t232 + (t104 * t217 + t106 * t218) * t233 + (t494 + (t128 * t217 + t130 * t218 - t496) * t283) * qJD(5) + (((-t21 + t480) * qJD(5) + t398) * t283 + t297) * t290) * t520 + (((-t105 * t294 - t107 * t292 + t82) * t232 + (-t104 * t294 - t106 * t292 + t83) * t233 + t49 * qJD(5)) * t284 + ((t321 + (-t128 * t294 - t130 * t292 - t314) * t284 - t397) * qJD(5) + t301) * t283) * t421 + ((-t526 * t290 + (t319 - t545) * t289 + t552) * t455 + t544 * t286 * qJD(2)) * t423 + ((t319 * t289 + (-t526 - t544) * t290 + t552) * t456 + t545 * qJD(2) * t287) * t455 / 0.2e1 + ((-t5 + t6) * t290 + (t4 + t7) * t289) * t420 + (-g(1) * (-t290 * t513 + t262 + t463) - g(2) * (-t289 * t513 + t261 + t464) - g(3) * (t133 + t431 + t510) - t527 * t283 * (-rSges(6,3) - pkin(3) - pkin(6)) + t3 * t415 + (t14 * t340 + t39 * t317 + t3 * (t236 + t88)) * t290 + (t13 * t340 + t38 * t317 + t3 * (t237 + t89)) * t289 - t39 * (t133 * t233 + t259) - t38 * (-t133 * t232 + t257) - ((t38 * t88 - t39 * t89) * t284 + (t39 * (-t134 * t289 - t108) + t38 * (t134 * t290 + t109) + t352) * t283) * qJD(5) - (t289 * t38 + t290 * t39) * (-t529 + t407) * qJD(2) + (-pkin(6) * t553 + t65 * t290 + t66 * t289 - t108 * t232 + t109 * t233 - (-t414 - t413) * qJD(2) + t551) * t19) * m(6) + (t17 * t415 + (t17 * t157 + t359 * t76 + t41 * t412) * t290 + (t17 * t158 + t359 * t75 + t40 * t412) * t289 - g(1) * (t290 * t408 + t262 + t460) - g(2) * (t289 * t408 + t261 + t461) - t76 * t259 - t75 * t257 + (-(t289 * t461 + t290 * t460 - t413) * qJD(2) + t541 + t551) * t25 + (-g(3) + (t289 * t75 + t290 * t76) * qJD(2)) * (t531 + t431)) * m(5) + (-g(3) * t532 - t527 * t323 - (-t48 * t413 + (-t115 * t532 - t290 * t351) * t290 + (-t114 * t532 - t289 * t351) * t289) * qJD(2) + t20 * t470 + t48 * t465 + (t115 * t409 + t20 * t156 - t48 * t177 + t323 * t74) * t290 + (t114 * t409 + t20 * t155 - t48 * t175 + t323 * t73) * t289) * m(4) + (g(1) * t235 + g(2) * t234 - g(3) * t270 + t50 * t363 + (-t362 - (-t234 * t289 - t235 * t290) * qJD(2)) * (qJD(2) * t363 + qJD(1)) + (qJDD(2) * t269 + t270 * t296) * t565) * m(3) + (t1 + (t559 * t287 + (t561 * t289 - t290 * t562) * t289) * t556 + (t557 * t287 + (t563 * t289 - t290 * t564) * t289) * t555) * t289 / 0.2e1 - (t2 + (t562 * t287 + (t559 - t561) * t567) * t556 + (t564 * t287 + (t557 - t563) * t567) * t555) * t290 / 0.2e1; t446 * (g(1) * t289 - g(2) * t290) + m(4) * (t289 * t74 - t290 * t73) + m(5) * (t289 * t41 - t290 * t40) + m(6) * t401; -t525 * g(3) * t284 + 0.2e1 * (t19 * t522 + t3 * t516) * m(6) + 0.2e1 * (t17 * t516 + t25 * t522) * m(5) + (t525 * t527 + m(5) * (qJD(2) * t25 + t289 * t40 + t290 * t41) + m(6) * (qJD(2) * t19 + t13 * t289 + t14 * t290)) * t283; -t6 * t429 / 0.2e1 + t1 * t476 / 0.2e1 + (t283 * t46 + t284 * t400) * t524 + (t15 * t283 + (t289 * t9 + t290 * t8) * t284 + (-t283 * t400 + t494) * qJD(2)) * t519 + t283 * t7 * t423 + t2 * t477 / 0.2e1 + (t283 * t47 + t284 * t399) * t523 + (t16 * t283 + (t10 * t290 + t11 * t289) * t284 + (-t283 * t399 + t493) * qJD(2)) * t517 + t12 * t457 / 0.2e1 + t283 * (t124 * t26 + t125 * t27 + t18 * t449 + t225 * t49 + t232 * t4 + t233 * t5) / 0.2e1 + (t283 * t49 + t284 * t397) * t521 + (t18 * t283 + (t289 * t5 + t290 * t4) * t284 + (-t283 * t397 + t49 * t284) * qJD(2)) * t420 + (t303 * t217 - t218 * t304 + t290 * t313) * t520 + (t219 * t303 + t220 * t304 + t289 * t313) * t518 + (t318 * t283 + (t304 * t292 - t294 * t303) * t284) * t421 + ((t13 * t88 - t14 * t89 + t38 * t65 - t39 * t66 + (t352 + (-t289 * t39 + t290 * t38) * t134) * qJD(2)) * t283 + (t39 * (-qJD(2) * t89 + t289 * t81) + t38 * (qJD(2) * t88 - t290 * t81) - t3 * t396 + t19 * (-t289 * t65 + t290 * t66) + t401 * t134) * t284 - t39 * (t208 * t233 - t449 * t97) - t38 * (-t208 * t232 + t449 * t96) - t19 * (t232 * t97 - t233 * t96) - g(1) * t96 - g(2) * t97 - g(3) * t208) * m(6);];
tau = t28;
