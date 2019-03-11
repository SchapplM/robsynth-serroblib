% Calculate time derivative of joint inertia matrix for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:33:53
% EndTime: 2019-03-09 03:34:13
% DurationCPUTime: 14.59s
% Computational Cost: add. (35354->728), mult. (26072->999), div. (0->0), fcn. (24063->12), ass. (0->395)
t558 = Icges(4,3) + Icges(5,3);
t307 = qJ(3) + pkin(11);
t298 = sin(t307);
t300 = cos(t307);
t311 = sin(qJ(3));
t314 = cos(qJ(3));
t557 = Icges(4,5) * t314 + Icges(5,5) * t300 - Icges(4,6) * t311 - Icges(5,6) * t298;
t308 = qJ(1) + pkin(10);
t299 = sin(t308);
t301 = cos(t308);
t553 = t299 * t558 + t557 * t301;
t497 = Icges(5,4) * t300;
t374 = -Icges(5,2) * t298 + t497;
t209 = Icges(5,6) * t299 + t301 * t374;
t498 = Icges(5,4) * t298;
t380 = Icges(5,1) * t300 - t498;
t211 = Icges(5,5) * t299 + t301 * t380;
t499 = Icges(4,4) * t314;
t376 = -Icges(4,2) * t311 + t499;
t223 = Icges(4,6) * t299 + t301 * t376;
t500 = Icges(4,4) * t311;
t382 = Icges(4,1) * t314 - t500;
t225 = Icges(4,5) * t299 + t301 * t382;
t530 = t209 * t298 - t211 * t300 + t223 * t311 - t225 * t314;
t556 = t299 * t530;
t555 = t530 * t301;
t554 = t557 * t299 - t301 * t558;
t552 = (-Icges(4,5) * t311 - Icges(5,5) * t298 - Icges(4,6) * t314 - Icges(5,6) * t300) * qJD(3);
t208 = -Icges(5,6) * t301 + t299 * t374;
t210 = -Icges(5,5) * t301 + t299 * t380;
t222 = -Icges(4,6) * t301 + t299 * t376;
t224 = -Icges(4,5) * t301 + t299 * t382;
t551 = t208 * t298 - t210 * t300 + t222 * t311 - t224 * t314;
t302 = qJ(5) + t307;
t293 = sin(t302);
t306 = qJD(3) + qJD(5);
t313 = cos(qJ(6));
t310 = sin(qJ(6));
t493 = Icges(7,4) * t313;
t371 = -Icges(7,2) * t310 + t493;
t294 = cos(t302);
t464 = t294 * t306;
t494 = Icges(7,4) * t310;
t144 = t371 * t464 + (Icges(7,6) * t306 + (-Icges(7,2) * t313 - t494) * qJD(6)) * t293;
t215 = -Icges(7,6) * t294 + t293 * t371;
t377 = Icges(7,1) * t313 - t494;
t216 = -Icges(7,5) * t294 + t293 * t377;
t550 = -t144 * t310 + (-t215 * t313 - t216 * t310) * qJD(6);
t367 = Icges(7,5) * t313 - Icges(7,6) * t310;
t143 = t367 * t464 + (Icges(7,3) * t306 + (-Icges(7,5) * t310 - Icges(7,6) * t313) * qJD(6)) * t293;
t456 = t306 * t310;
t549 = -t215 * t456 - t143;
t548 = t553 * qJD(1);
t547 = t554 * t299 - t556 + (-t551 - t553) * t301;
t296 = t299 ^ 2;
t297 = t301 ^ 2;
t248 = Icges(6,5) * t293 + Icges(6,6) * t294;
t338 = t306 * t248;
t368 = Icges(6,5) * t294 - Icges(6,6) * t293;
t195 = -Icges(6,3) * t301 + t299 * t368;
t527 = qJD(1) * t195;
t126 = -t301 * t338 - t527;
t196 = Icges(6,3) * t299 + t301 * t368;
t438 = qJD(1) * t196;
t127 = -t299 * t338 + t438;
t495 = Icges(6,4) * t294;
t372 = -Icges(6,2) * t293 + t495;
t197 = -Icges(6,6) * t301 + t299 * t372;
t496 = Icges(6,4) * t293;
t249 = Icges(6,2) * t294 + t496;
t470 = t249 * t306;
t128 = -qJD(1) * t197 - t301 * t470;
t378 = Icges(6,1) * t294 - t496;
t199 = -Icges(6,5) * t301 + t299 * t378;
t250 = Icges(6,1) * t293 + t495;
t469 = t250 * t306;
t130 = -qJD(1) * t199 - t301 * t469;
t198 = Icges(6,6) * t299 + t301 * t372;
t200 = Icges(6,5) * t299 + t301 * t378;
t357 = t198 * t293 - t200 * t294;
t131 = qJD(1) * t200 - t299 * t469;
t403 = t197 * t306 - t131;
t129 = qJD(1) * t198 - t299 * t470;
t405 = t199 * t306 + t129;
t466 = t293 * t306;
t358 = t197 * t293 - t199 * t294;
t536 = t301 * t358;
t88 = -t195 * t301 - t299 * t358;
t538 = t299 * t357;
t89 = -t196 * t301 - t538;
t401 = -qJD(6) * t294 + qJD(1);
t322 = t293 * t456 + t313 * t401;
t400 = qJD(1) * t294 - qJD(6);
t458 = t301 * t310;
t135 = t299 * t322 - t400 * t458;
t455 = t306 * t313;
t321 = -t293 * t455 + t310 * t401;
t457 = t301 * t313;
t136 = t299 * t321 + t400 * t457;
t462 = t299 * t310;
t236 = -t294 * t462 - t457;
t461 = t299 * t313;
t237 = t294 * t461 - t458;
t468 = t293 * t299;
t148 = Icges(7,5) * t237 + Icges(7,6) * t236 + Icges(7,3) * t468;
t150 = Icges(7,4) * t237 + Icges(7,2) * t236 + Icges(7,6) * t468;
t152 = Icges(7,1) * t237 + Icges(7,4) * t236 + Icges(7,5) * t468;
t434 = qJD(1) * t301;
t463 = t299 * t306;
t327 = t293 * t434 + t294 * t463;
t72 = Icges(7,5) * t136 + Icges(7,6) * t135 + Icges(7,3) * t327;
t74 = Icges(7,4) * t136 + Icges(7,2) * t135 + Icges(7,6) * t327;
t76 = Icges(7,1) * t136 + Icges(7,4) * t135 + Icges(7,5) * t327;
t16 = t135 * t150 + t136 * t152 + t148 * t327 + t236 * t74 + t237 * t76 + t468 * t72;
t238 = -t294 * t458 + t461;
t239 = t294 * t457 + t462;
t467 = t293 * t301;
t149 = Icges(7,5) * t239 + Icges(7,6) * t238 + Icges(7,3) * t467;
t151 = Icges(7,4) * t239 + Icges(7,2) * t238 + Icges(7,6) * t467;
t153 = Icges(7,1) * t239 + Icges(7,4) * t238 + Icges(7,5) * t467;
t133 = t301 * t322 + t400 * t462;
t134 = t301 * t321 - t400 * t461;
t435 = qJD(1) * t299;
t413 = t293 * t435;
t460 = t301 * t306;
t419 = t294 * t460;
t326 = -t413 + t419;
t71 = Icges(7,5) * t134 + Icges(7,6) * t133 + Icges(7,3) * t326;
t73 = Icges(7,4) * t134 + Icges(7,2) * t133 + Icges(7,6) * t326;
t75 = Icges(7,1) * t134 + Icges(7,4) * t133 + Icges(7,5) * t326;
t17 = t135 * t151 + t136 * t153 + t149 * t327 + t236 * t73 + t237 * t75 + t468 * t71;
t54 = t148 * t468 + t150 * t236 + t152 * t237;
t55 = t149 * t468 + t151 * t236 + t153 * t237;
t387 = t299 * t54 + t301 * t55;
t9 = qJD(1) * t387 - t16 * t301 + t17 * t299;
t545 = -(t127 * t301 + (t89 + t536) * qJD(1)) * t301 - (t88 * qJD(1) + (-t128 * t293 + t130 * t294 - t198 * t464 - t200 * t466 + t438) * t299 + (-t126 + t403 * t294 + t405 * t293 + (-t195 - t357) * qJD(1)) * t301) * t299 - t9;
t303 = t314 * pkin(3);
t295 = t303 + pkin(2);
t273 = t301 * t295;
t309 = -qJ(4) - pkin(7);
t509 = -pkin(7) - t309;
t194 = -pkin(2) * t301 + t299 * t509 + t273;
t504 = rSges(5,2) * t298;
t507 = rSges(5,1) * t300;
t393 = -t504 + t507;
t212 = -rSges(5,3) * t301 + t299 * t393;
t290 = t299 * rSges(5,3);
t532 = -t301 * t504 + t290;
t213 = t301 * t507 + t532;
t260 = rSges(5,1) * t298 + rSges(5,2) * t300;
t341 = qJD(3) * t260;
t292 = t301 * pkin(7);
t459 = t301 * t309;
t510 = pkin(2) - t295;
t537 = t299 * t510;
t193 = t292 + t459 - t537;
t287 = qJD(4) * t299;
t429 = qJD(3) * t311;
t424 = pkin(3) * t429;
t402 = t301 * t424;
t288 = qJD(4) * t301;
t442 = t299 * t424 + t309 * t435;
t414 = -t288 - t442;
t511 = pkin(7) * t299;
t415 = t299 * ((-t301 * t510 - t511) * qJD(1) + t414) + t301 * (-t402 + t287 + (t301 * t509 + t537) * qJD(1)) + t193 * t434;
t412 = t298 * t435;
t446 = rSges(5,2) * t412 + rSges(5,3) * t434;
t45 = (qJD(1) * t212 - t301 * t341 + t446) * t301 + (-t299 * t341 + (-t194 - t213 + t532) * qJD(1)) * t299 + t415;
t524 = 2 * m(5);
t544 = t45 * t524;
t525 = 2 * m(4);
t505 = rSges(4,2) * t311;
t508 = rSges(4,1) * t314;
t394 = -t505 + t508;
t503 = rSges(4,3) * t301;
t226 = t299 * t394 - t503;
t426 = t301 * t505;
t291 = t299 * rSges(4,3);
t441 = t301 * t508 + t291;
t227 = -t426 + t441;
t282 = rSges(4,1) * t311 + rSges(4,2) * t314;
t342 = qJD(3) * t282;
t433 = qJD(1) * t311;
t411 = t299 * t433;
t319 = rSges(4,2) * t411 + rSges(4,3) * t434 - t301 * t342;
t68 = (qJD(1) * t226 + t319) * t301 + (-t299 * t342 + (-t227 - t426 + t291) * qJD(1)) * t299;
t543 = t525 * t68;
t541 = t299 * t553 - t555;
t540 = t299 * t551 + t301 * t554;
t251 = rSges(6,1) * t293 + rSges(6,2) * t294;
t343 = t251 * t306;
t390 = -rSges(7,1) * t237 - rSges(7,2) * t236;
t155 = rSges(7,3) * t468 - t390;
t156 = t239 * rSges(7,1) + t238 * rSges(7,2) + rSges(7,3) * t467;
t535 = -t299 * t155 - t301 * t156;
t534 = -qJD(1) * t554 + t301 * t552;
t533 = -t299 * t552 - t548;
t516 = sin(qJ(1)) * pkin(1);
t531 = t292 - t516;
t515 = pkin(3) * t311;
t398 = -pkin(4) * t298 - t515;
t350 = -t251 + t398;
t180 = t350 * t299;
t181 = t350 * t301;
t528 = t180 * t299 + t181 * t301;
t352 = t249 * t293 - t250 * t294;
t526 = qJD(1) * t352 + t368 * t306;
t305 = -pkin(8) + t309;
t439 = t305 - t309;
t514 = pkin(4) * t300;
t270 = t295 + t514;
t444 = t270 - t295;
t178 = t299 * t444 + t301 * t439;
t523 = 2 * m(6);
t522 = 2 * m(7);
t521 = t299 / 0.2e1;
t520 = -t301 / 0.2e1;
t519 = -rSges(7,3) - pkin(9);
t518 = m(4) * t282;
t517 = m(6) * t251;
t513 = pkin(5) * t293;
t512 = pkin(5) * t294;
t304 = cos(qJ(1)) * pkin(1);
t506 = rSges(6,1) * t294;
t365 = -t150 * t310 + t152 * t313;
t20 = (t306 * t365 - t72) * t294 + (t148 * t306 - t310 * t74 + t313 * t76 + (-t150 * t313 - t152 * t310) * qJD(6)) * t293;
t502 = t20 * t301;
t364 = -t151 * t310 + t153 * t313;
t21 = (t306 * t364 - t71) * t294 + (t149 * t306 - t310 * t73 + t313 * t75 + (-t151 * t313 - t153 * t310) * qJD(6)) * t293;
t501 = t21 * t299;
t289 = t299 * rSges(6,3);
t480 = t208 * t300;
t479 = t209 * t300;
t478 = t210 * t298;
t477 = t211 * t298;
t474 = t222 * t314;
t473 = t223 * t314;
t472 = t224 * t311;
t471 = t225 * t311;
t465 = t294 * t301;
t389 = rSges(7,1) * t313 - rSges(7,2) * t310;
t154 = t389 * t464 + (rSges(7,3) * t306 + (-rSges(7,1) * t310 - rSges(7,2) * t313) * qJD(6)) * t293;
t397 = pkin(9) * t293 + t512;
t454 = -t397 * t306 - t154;
t229 = pkin(5) * t465 + pkin(9) * t467;
t453 = -t156 - t229;
t243 = t301 * t270;
t179 = -t299 * t439 + t243 - t273;
t452 = -t179 - t194;
t451 = t299 * t193 + t301 * t194;
t392 = -rSges(6,2) * t293 + t506;
t201 = -rSges(6,3) * t301 + t299 * t392;
t202 = rSges(6,1) * t465 - rSges(6,2) * t467 + t289;
t118 = t299 * t201 + t301 * t202;
t217 = -rSges(7,3) * t294 + t293 * t389;
t192 = t217 * t435;
t253 = -pkin(9) * t294 + t513;
t450 = t253 * t435 + t192;
t449 = -t217 - t253;
t258 = t398 * qJD(3);
t242 = t301 * t258;
t448 = t242 + t287;
t447 = rSges(6,2) * t413 + rSges(6,3) * t434;
t277 = pkin(3) * t411;
t445 = pkin(4) * t412 + t277;
t275 = t305 * t435;
t443 = t275 + t288;
t440 = t296 + t297;
t432 = qJD(3) * t298;
t431 = qJD(3) * t299;
t430 = qJD(3) * t300;
t428 = qJD(3) * t314;
t63 = -t148 * t294 + t293 * t365;
t214 = -Icges(7,3) * t294 + t293 * t367;
t84 = t214 * t468 + t215 * t236 + t216 * t237;
t423 = t63 / 0.2e1 + t84 / 0.2e1;
t64 = -t149 * t294 + t293 * t364;
t85 = t214 * t467 + t215 * t238 + t216 * t239;
t422 = t85 / 0.2e1 + t64 / 0.2e1;
t420 = t293 * t460;
t328 = -t294 * t435 - t420;
t418 = t299 * (-t299 * t343 + (t301 * t392 + t289) * qJD(1)) + t301 * (rSges(6,1) * t328 - rSges(6,2) * t419 + t447) + t201 * t434;
t145 = t377 * t464 + (Icges(7,5) * t306 + (-Icges(7,1) * t310 - t493) * qJD(6)) * t293;
t417 = t293 * t313 * t145 + t294 * t216 * t455 + t214 * t466;
t416 = t134 * rSges(7,1) + t133 * rSges(7,2) + rSges(7,3) * t419;
t410 = t464 / 0.2e1;
t409 = t435 / 0.2e1;
t408 = t434 / 0.2e1;
t407 = -t260 - t515;
t166 = t449 * t301;
t406 = t200 * t306 + t128;
t404 = -t198 * t306 + t130;
t399 = t299 * t178 + t301 * t179 + t451;
t228 = t397 * t299;
t67 = t299 * t228 + t301 * t229 - t535;
t396 = -t299 * t305 + t243 + t304;
t219 = t407 * t301;
t391 = t136 * rSges(7,1) + t135 * rSges(7,2);
t388 = -t301 * t305 - t516;
t39 = t299 * t55 - t301 * t54;
t56 = t148 * t467 + t150 * t238 + t152 * t239;
t57 = t149 * t467 + t151 * t238 + t153 * t239;
t40 = t299 * t57 - t301 * t56;
t386 = t299 * t56 + t301 * t57;
t385 = t299 * t64 - t301 * t63;
t384 = t299 * t63 + t301 * t64;
t14 = t133 * t150 + t134 * t152 + t148 * t326 + t238 * t74 + t239 * t76 + t467 * t72;
t15 = t133 * t151 + t134 * t153 + t149 * t326 + t238 * t73 + t239 * t75 + t467 * t71;
t8 = qJD(1) * t386 - t14 * t301 + t15 * t299;
t90 = t195 * t299 - t536;
t91 = t196 * t299 - t301 * t357;
t383 = t39 * t435 + t40 * t434 + (-t90 * t434 - t88 * t435) * t301 + (t8 + (t91 * qJD(1) + (t129 * t293 - t131 * t294 + t197 * t464 + t199 * t466 - t527) * t301) * t301 + t89 * t435 + t91 * t434 + ((t90 + t538) * qJD(1) + (-t127 + t404 * t294 - t406 * t293 + (t196 - t358) * qJD(1)) * t301 + t126 * t299) * t299) * t299;
t381 = Icges(4,1) * t311 + t499;
t379 = Icges(5,1) * t298 + t497;
t375 = Icges(4,2) * t314 + t500;
t373 = Icges(5,2) * t300 + t498;
t347 = -t270 - t392;
t146 = -t516 + (rSges(6,3) - t305) * t301 + t347 * t299;
t147 = t202 + t396;
t366 = t146 * t301 + t147 * t299;
t363 = t155 * t301 - t156 * t299;
t351 = -pkin(2) - t394;
t254 = pkin(9) * t419;
t77 = -rSges(7,3) * t413 + t416;
t78 = rSges(7,3) * t327 + t391;
t349 = (t155 + t228) * t434 + (pkin(5) * t328 - pkin(9) * t413 + t254 + t77) * t301 + (t78 + t327 * pkin(9) + (-t293 * t463 + t294 * t434) * pkin(5)) * t299;
t348 = -t295 - t393;
t345 = t299 * (t299 * t258 + t434 * t444 - t275 + t442) + t301 * (-qJD(1) * t178 + t242 + t402) + t178 * t434 + t415;
t344 = (-t303 - t514) * qJD(3);
t340 = t398 + t449;
t336 = qJD(3) * t381;
t335 = qJD(3) * t379;
t334 = qJD(3) * t375;
t333 = qJD(3) * t373;
t330 = t293 * t519 - t270 - t512;
t235 = t392 * t306;
t325 = -t235 + t344;
t117 = t340 * t301;
t324 = t301 * t545 + t383;
t323 = t344 + t454;
t26 = t293 * t387 - t294 * t84;
t27 = t293 * t386 - t294 * t85;
t31 = t133 * t215 + t134 * t216 + t143 * t467 + t144 * t238 + t145 * t239 + t214 * t326;
t3 = (t306 * t386 - t31) * t294 + (-qJD(1) * t40 + t14 * t299 + t15 * t301 + t306 * t85) * t293;
t32 = t135 * t215 + t136 * t216 + t143 * t468 + t144 * t236 + t145 * t237 + t214 * t327;
t4 = (t306 * t387 - t32) * t294 + (-qJD(1) * t39 + t16 * t299 + t17 * t301 + t306 * t84) * t293;
t320 = t3 * t521 + t4 * t520 - t294 * (qJD(1) * t384 + t501 - t502) / 0.2e1 + t26 * t409 + t27 * t408 + t385 * t466 / 0.2e1 + t9 * t468 / 0.2e1 + t8 * t467 / 0.2e1 + (t301 * t410 - t413 / 0.2e1) * t40 + (t293 * t408 + t299 * t410) * t39;
t318 = t299 * t330 + t388;
t233 = t372 * t306;
t234 = t378 * t306;
t317 = qJD(1) * t248 + (t234 - t470) * t294 + (-t233 - t469) * t293;
t316 = -t502 / 0.2e1 + t501 / 0.2e1 + (t293 * t404 + t294 * t406 + t299 * t526 + t317 * t301 + t31) * t521 + (-t293 * t403 + t294 * t405 + t317 * t299 - t301 * t526 + t32) * t520 + (t197 * t294 + t199 * t293 - t248 * t301 - t299 * t352 + t63 + t84) * t409 + (t198 * t294 + t200 * t293 + t248 * t299 - t301 * t352 + t64 + t85) * t408;
t271 = t394 * qJD(3);
t247 = t393 * qJD(3);
t218 = t407 * t299;
t183 = t511 + t304 + (pkin(2) - t505) * t301 + t441;
t182 = t299 * t351 + t503 + t531;
t169 = -t299 * t309 + t213 + t273 + t304;
t168 = -t516 + (rSges(5,3) - t309) * t301 + t348 * t299;
t165 = t449 * t299;
t142 = -t260 * t434 - t247 * t299 + (-t299 * t428 - t301 * t433) * pkin(3);
t141 = t260 * t435 + t277 + (-pkin(3) * t428 - t247) * t301;
t123 = t282 * t431 + (-t304 + (-rSges(4,3) - pkin(7)) * t299 + t351 * t301) * qJD(1);
t122 = ((-pkin(2) - t508) * t299 + t531) * qJD(1) + t319;
t116 = t340 * t299;
t109 = qJD(1) * t181 + t299 * t325;
t108 = t251 * t435 + t301 * t325 + t445;
t103 = t260 * t431 + (t301 * t348 - t290 - t304) * qJD(1) - t414;
t102 = t287 + qJD(3) * t219 + (-t516 - t459 + (-t295 - t507) * t299) * qJD(1) + t446;
t101 = -t156 * t294 - t217 * t467;
t100 = t155 * t294 + t217 * t468;
t99 = -t214 * t294 + (-t215 * t310 + t216 * t313) * t293;
t98 = t396 - t453;
t97 = t318 + t390;
t92 = t99 * t466;
t87 = (-t258 + t343) * t299 + (t301 * t347 - t289 - t304) * qJD(1) + t443;
t86 = -t301 * t343 + ((-t270 - t506) * t299 + t388) * qJD(1) + t447 + t448;
t83 = t363 * t293;
t82 = qJD(1) * t166 + t299 * t454;
t81 = t301 * t454 + t450;
t66 = qJD(1) * t117 + t299 * t323;
t65 = t301 * t323 + t445 + t450;
t60 = -t202 * t435 + t418;
t51 = t399 + t118;
t48 = (-t258 + (t294 * t519 + t513) * t306) * t299 + (t301 * t330 - t304) * qJD(1) - t391 + t443;
t47 = -pkin(5) * t420 + qJD(1) * t318 + t254 + t416 + t448;
t46 = t67 + t399;
t44 = (t217 * t463 + t78) * t294 + (t154 * t299 - t155 * t306 + t217 * t434) * t293;
t43 = (-t217 * t460 - t77) * t294 + (-t154 * t301 + t156 * t306 + t192) * t293;
t42 = t293 * t550 + t549 * t294 + t417;
t28 = t363 * t464 + (qJD(1) * t535 - t299 * t77 + t301 * t78) * t293;
t23 = t435 * t453 + t349;
t22 = (-t202 + t452) * t435 + t345 + t418;
t13 = (t452 + t453) * t435 + t345 + t349;
t1 = [t250 * t464 - t249 * t466 + (t122 * t183 + t123 * t182) * t525 + (t102 * t169 + t103 * t168) * t524 + (t146 * t87 + t147 * t86) * t523 + (t47 * t98 + t48 * t97) * t522 + t417 + (-t373 + t380) * t432 + (t374 + t379) * t430 + (t382 - t375) * t429 + (t381 + t376) * t428 + (t233 + t549) * t294 + (t234 + t550) * t293; 0; 0; m(5) * (t102 * t218 + t103 * t219 + t141 * t168 + t142 * t169) + m(6) * (t108 * t146 + t109 * t147 + t180 * t86 + t181 * t87) + m(7) * (t116 * t47 + t117 * t48 + t65 * t97 + t66 * t98) + m(4) * ((-t122 * t299 - t123 * t301) * t282 + (-t182 * t301 - t183 * t299) * t271) + t316 + ((t473 / 0.2e1 + t471 / 0.2e1 + t479 / 0.2e1 + t477 / 0.2e1 - t183 * t518) * t301 + (t474 / 0.2e1 + t472 / 0.2e1 + t182 * t518 + t480 / 0.2e1 + t478 / 0.2e1) * t299) * qJD(1) + (-qJD(3) * t530 + (-qJD(1) * t208 - t301 * t333) * t300 + (-qJD(1) * t210 - t301 * t335) * t298 + (-qJD(1) * t222 - t301 * t334) * t314 + (-qJD(1) * t224 - t301 * t336) * t311) * t521 + (-qJD(3) * t551 + (qJD(1) * t209 - t299 * t333) * t300 + (qJD(1) * t211 - t299 * t335) * t298 + (qJD(1) * t223 - t299 * t334) * t314 + (qJD(1) * t225 - t299 * t336) * t311) * t520 + t557 * qJD(3) * (t296 / 0.2e1 + t297 / 0.2e1); m(4) * t68 + m(5) * t45 + m(6) * t22 + m(7) * t13; (t116 * t66 + t117 * t65 + t13 * t46) * t522 + (t108 * t181 + t109 * t180 + t22 * t51) * t523 + (t219 * t141 + t218 * t142 + t45 * t451) * t524 + t440 * t282 * t271 * t525 + t383 + (t212 * t544 + t226 * t543 + t534 * t296 + t541 * t434 + (t547 + t556) * t435) * t299 + (t213 * t544 + t227 * t543 + t533 * t297 + t540 * t435 + ((t208 * t430 + t210 * t432 + t222 * t428 + t224 * t429 + t534) * t301 + (t209 * t430 + t211 * t432 + t223 * t428 + t225 * t429 + t533 - t548) * t299 + ((-t472 - t474 - t478 - t480) * t301 + (-t471 - t473 - t477 - t479) * t299) * qJD(3) + (t555 + (-t551 + t553) * t299 + t540 + t541) * qJD(1)) * t299 + t545 + (-t301 * t551 - t547) * t434) * t301; m(7) * (t299 * t48 - t301 * t47 + (t299 * t98 + t301 * t97) * qJD(1)) + m(6) * (qJD(1) * t366 + t299 * t87 - t301 * t86) + m(5) * (-t102 * t301 + t103 * t299 + (t168 * t301 + t169 * t299) * qJD(1)); 0; m(7) * (t299 * t65 - t301 * t66 + (t116 * t299 + t117 * t301) * qJD(1)) + m(6) * (qJD(1) * t528 + t108 * t299 - t109 * t301) + m(5) * (t141 * t299 - t142 * t301 + (t218 * t299 + t219 * t301) * qJD(1)); 0; t316 + m(7) * (t165 * t47 + t166 * t48 + t81 * t97 + t82 * t98) - m(6) * t366 * t235 + (-t299 * t86 - t301 * t87 + (t146 * t299 - t147 * t301) * qJD(1)) * t517; m(6) * t60 + m(7) * t23; m(7) * (t116 * t82 + t117 * t81 + t13 * t67 + t165 * t66 + t166 * t65 + t23 * t46) + m(6) * (t118 * t22 - t235 * t528 + t51 * t60) + (-t108 * t301 - t109 * t299 + (-t180 * t301 + t181 * t299) * qJD(1)) * t517 + t324; m(7) * (t299 * t81 - t301 * t82 + (t165 * t299 + t166 * t301) * qJD(1)); (t235 * t251 * t440 + t118 * t60) * t523 + (t165 * t82 + t166 * t81 + t23 * t67) * t522 + t324; m(7) * (t100 * t48 + t101 * t47 + t43 * t98 + t44 * t97) + t92 + (-t42 + (t299 * t423 + t301 * t422) * t306) * t294 + ((t31 / 0.2e1 + t21 / 0.2e1) * t301 + (t20 / 0.2e1 + t32 / 0.2e1) * t299 + (-t299 * t422 + t301 * t423) * qJD(1)) * t293; m(7) * t28; t320 + m(7) * (t100 * t65 + t101 * t66 + t116 * t43 + t117 * t44 + t13 * t83 + t28 * t46); m(7) * (t299 * t44 - t301 * t43 + (t100 * t301 + t101 * t299) * qJD(1)); t320 + m(7) * (t100 * t81 + t101 * t82 + t165 * t43 + t166 * t44 + t23 * t83 + t28 * t67); (t100 * t44 + t101 * t43 + t28 * t83) * t522 + (t42 * t294 - t92 + (t299 * t26 + t301 * t27 - t294 * t384) * t306) * t294 + (t301 * t3 + t299 * t4 + t384 * t466 + (-t20 * t299 - t21 * t301 - t306 * t99) * t294 + (t301 * t26 - t299 * t27 + t294 * t385) * qJD(1)) * t293;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
