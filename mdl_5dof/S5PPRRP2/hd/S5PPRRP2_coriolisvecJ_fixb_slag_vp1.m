% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:32
% EndTime: 2019-12-05 15:09:04
% DurationCPUTime: 24.55s
% Computational Cost: add. (14033->574), mult. (20831->835), div. (0->0), fcn. (20296->6), ass. (0->301)
t502 = Icges(6,4) + Icges(5,5);
t501 = Icges(5,6) - Icges(6,6);
t239 = pkin(8) + qJ(3);
t238 = cos(t239);
t241 = cos(pkin(7));
t242 = sin(qJ(4));
t384 = t241 * t242;
t240 = sin(pkin(7));
t243 = cos(qJ(4));
t385 = t240 * t243;
t218 = t238 * t385 - t384;
t200 = Icges(6,5) * t218;
t383 = t241 * t243;
t386 = t240 * t242;
t217 = t238 * t386 + t383;
t237 = sin(t239);
t394 = t237 * t240;
t92 = Icges(6,6) * t394 + Icges(6,3) * t217 + t200;
t408 = Icges(5,4) * t218;
t98 = -Icges(5,2) * t217 + Icges(5,6) * t394 + t408;
t493 = t92 - t98;
t220 = t238 * t383 + t386;
t201 = Icges(6,5) * t220;
t219 = t238 * t384 - t385;
t393 = t237 * t241;
t93 = Icges(6,6) * t393 + Icges(6,3) * t219 + t201;
t407 = Icges(5,4) * t220;
t99 = -Icges(5,2) * t219 + Icges(5,6) * t393 + t407;
t492 = t93 - t99;
t94 = Icges(5,5) * t218 - Icges(5,6) * t217 + Icges(5,3) * t394;
t96 = Icges(6,4) * t218 + Icges(6,2) * t394 + Icges(6,6) * t217;
t491 = t96 + t94;
t95 = Icges(5,5) * t220 - Icges(5,6) * t219 + Icges(5,3) * t393;
t97 = Icges(6,4) * t220 + Icges(6,2) * t393 + Icges(6,6) * t219;
t490 = t97 + t95;
t402 = Icges(6,5) * t217;
t100 = Icges(6,1) * t218 + Icges(6,4) * t394 + t402;
t202 = Icges(5,4) * t217;
t102 = Icges(5,1) * t218 + Icges(5,5) * t394 - t202;
t512 = t100 + t102;
t401 = Icges(6,5) * t219;
t101 = Icges(6,1) * t220 + Icges(6,4) * t393 + t401;
t203 = Icges(5,4) * t219;
t103 = Icges(5,1) * t220 + Icges(5,5) * t393 - t203;
t511 = t101 + t103;
t480 = Icges(5,2) + Icges(6,3);
t532 = Icges(6,2) + Icges(5,3);
t391 = t237 * t243;
t232 = Icges(6,5) * t391;
t392 = t237 * t242;
t396 = Icges(6,6) * t238;
t150 = Icges(6,3) * t392 + t232 - t396;
t405 = Icges(5,4) * t243;
t306 = -Icges(5,2) * t242 + t405;
t156 = -Icges(5,6) * t238 + t237 * t306;
t510 = t150 - t156;
t303 = Icges(5,5) * t243 - Icges(5,6) * t242;
t152 = -Icges(5,3) * t238 + t237 * t303;
t305 = Icges(6,4) * t243 + Icges(6,6) * t242;
t154 = -Icges(6,2) * t238 + t237 * t305;
t531 = t152 + t154;
t400 = Icges(6,5) * t242;
t309 = Icges(6,1) * t243 + t400;
t158 = -Icges(6,4) * t238 + t237 * t309;
t406 = Icges(5,4) * t242;
t310 = Icges(5,1) * t243 - t406;
t160 = -Icges(5,5) * t238 + t237 * t310;
t478 = t158 + t160;
t530 = (-t502 * t242 - t501 * t243) * t237;
t499 = t492 * t217 + t511 * t218 + t490 * t394;
t498 = t493 * t219 + t512 * t220 + t491 * t393;
t500 = t493 * t217 + t512 * t218 + t491 * t394;
t464 = t492 * t219 + t511 * t220 + t490 * t393;
t529 = t217 * t510 + t218 * t478 + t394 * t531;
t528 = t219 * t510 + t220 * t478 + t393 * t531;
t357 = qJD(4) * t243;
t345 = t238 * t357;
t359 = qJD(4) * t241;
t365 = qJD(3) * t240;
t144 = -t240 * t345 + (t237 * t365 + t359) * t242;
t362 = qJD(3) * t243;
t350 = t237 * t362;
t145 = -qJD(4) * t217 - t240 * t350;
t349 = t238 * t365;
t527 = t501 * t144 + t502 * t145 + t532 * t349;
t358 = qJD(4) * t242;
t363 = qJD(3) * t242;
t146 = -t240 * t358 - t241 * t345 + t363 * t393;
t147 = -qJD(4) * t219 - t241 * t350;
t364 = qJD(3) * t241;
t348 = t238 * t364;
t526 = t501 * t146 + t502 * t147 + t532 * t348;
t153 = Icges(5,3) * t237 + t238 * t303;
t155 = Icges(6,2) * t237 + t238 * t305;
t525 = t530 * qJD(4) + (t153 + t155) * qJD(3);
t524 = Icges(5,1) + Icges(6,1);
t523 = Icges(5,4) - Icges(6,5);
t399 = Icges(6,5) * t243;
t302 = Icges(6,3) * t242 + t399;
t151 = Icges(6,6) * t237 + t238 * t302;
t157 = Icges(5,6) * t237 + t238 * t306;
t457 = -t151 + t157;
t159 = Icges(6,4) * t237 + t238 * t309;
t161 = Icges(5,5) * t237 + t238 * t310;
t456 = -t159 - t161;
t521 = (-t480 * t243 + t400 - t406) * t237;
t520 = t499 * t241;
t519 = t498 * t240;
t518 = -t480 * t144 - t523 * t145 - t501 * t349;
t517 = -t480 * t146 - t523 * t147 - t501 * t348;
t516 = t523 * t144 + t524 * t145 + t502 * t349;
t515 = t523 * t146 + t524 * t147 + t502 * t348;
t514 = t457 * qJD(3) + t521 * qJD(4);
t211 = (-Icges(5,1) * t242 - t405) * t237;
t361 = qJD(4) * t237;
t513 = -(-Icges(6,1) * t242 + t399) * t361 - qJD(4) * t211 + t456 * qJD(3);
t366 = qJD(3) * t238;
t509 = -t525 * t237 - t366 * t531;
t508 = t526 * t237 + t490 * t366;
t507 = t527 * t237 + t491 * t366;
t506 = t464 * t241 + t519;
t505 = t500 * t240 + t520;
t504 = t528 * t237;
t503 = t529 * t237;
t223 = t237 * t359 + t365;
t430 = t223 / 0.2e1;
t224 = t240 * t361 - t364;
t428 = t224 / 0.2e1;
t470 = -t493 * t144 + t512 * t145 + t518 * t217 + t516 * t218 + t507 * t240;
t469 = -t492 * t144 + t511 * t145 + t517 * t217 + t515 * t218 + t508 * t240;
t468 = -t493 * t146 + t512 * t147 + t518 * t219 + t516 * t220 + t507 * t241;
t467 = -t492 * t146 + t511 * t147 + t517 * t219 + t515 * t220 + t508 * t241;
t316 = t100 * t243 + t242 * t92;
t47 = t237 * t316 - t238 * t96;
t314 = t102 * t243 - t242 * t98;
t49 = t237 * t314 - t238 * t94;
t497 = t47 + t49;
t315 = t101 * t243 + t242 * t93;
t48 = t237 * t315 - t238 * t97;
t313 = t103 * t243 - t242 * t99;
t50 = t237 * t313 - t238 * t95;
t496 = t48 + t50;
t300 = t150 * t242 + t158 * t243;
t389 = t238 * t154;
t61 = t237 * t300 - t389;
t299 = -t156 * t242 + t160 * t243;
t390 = t238 * t152;
t62 = t237 * t299 - t390;
t482 = t61 + t62;
t328 = rSges(5,1) * t243 - rSges(5,2) * t242;
t165 = -t238 * rSges(5,3) + t237 * t328;
t489 = t165 * t223;
t488 = t165 * t224;
t485 = (t510 * t146 - t478 * t147 + t514 * t219 + t513 * t220 + t509 * t241) * t238 + (t506 * t238 + t504) * qJD(3);
t484 = (t510 * t144 - t478 * t145 + t514 * t217 + t513 * t218 + t509 * t240) * t238 + (t505 * t238 + t503) * qJD(3);
t483 = t468 * t428 + t467 * t430 + t485 * qJD(4) / 0.2e1;
t360 = qJD(4) * t238;
t477 = t530 * t360 + (t217 * t502 + t501 * t218) * t224 + (t219 * t502 + t501 * t220) * t223;
t317 = t49 * t240 + t50 * t241;
t318 = t47 * t240 + t48 * t241;
t476 = t317 + t318;
t432 = t241 ^ 2;
t433 = t240 ^ 2;
t437 = t432 + t433;
t448 = qJD(3) * t237;
t475 = qJD(4) * t484 + t223 * t469 + t224 * t470;
t445 = ((t314 + t316) * qJD(3) - t527) * t238 + (t516 * t243 + t518 * t242 + (-t242 * t512 + t243 * t493) * qJD(4) + t491 * qJD(3)) * t237;
t444 = ((t313 + t315) * qJD(3) - t526) * t238 + (t515 * t243 + t517 * t242 + (-t242 * t511 + t243 * t492) * qJD(4) + t490 * qJD(3)) * t237;
t473 = rSges(6,1) + pkin(4);
t472 = t223 * t499 + t224 * t500 - t360 * t529;
t471 = t223 * t464 + t224 * t498 - t360 * t528;
t466 = t223 * t496 + t224 * t497 - t360 * t482;
t463 = rSges(6,3) + qJ(5);
t225 = rSges(4,1) * t237 + rSges(4,2) * t238;
t462 = qJD(3) * t225 * t437;
t263 = -t237 * t302 + t396;
t128 = t263 * t240;
t134 = t156 * t240;
t461 = t128 + t134;
t129 = t263 * t241;
t135 = t156 * t241;
t460 = t129 + t135;
t136 = t158 * t240;
t138 = t160 * t240;
t459 = -t136 - t138;
t137 = t158 * t241;
t139 = t160 * t241;
t458 = -t137 - t139;
t283 = t153 - t299;
t284 = -t155 + t300;
t434 = (t154 * t241 + t315) * t223 + (t154 * t240 + t316) * t224;
t435 = (t152 * t241 + t313) * t223 + (t152 * t240 + t314) * t224;
t455 = (-t434 - t435 + (-t283 + t284) * t360) * t237;
t326 = pkin(4) * t243 + qJ(5) * t242;
t327 = rSges(6,1) * t243 + rSges(6,3) * t242;
t454 = t238 * rSges(6,2) + (-t326 - t327) * t237;
t453 = t223 * t490 + t224 * t491;
t346 = t237 * t357;
t452 = t238 * t363 + t346;
t227 = pkin(3) * t237 - pkin(6) * t238;
t282 = qJD(3) * t227;
t192 = t240 * t282;
t193 = t241 * t282;
t451 = -t240 * t192 - t241 * t193 + t437 * t282;
t450 = -t389 - t390;
t409 = Icges(4,4) * t238;
t410 = Icges(4,4) * t237;
t449 = -(-Icges(4,2) * t238 - t410) * t448 + (-Icges(4,1) * t237 - t409) * t366;
t308 = -Icges(4,2) * t237 + t409;
t312 = Icges(4,1) * t238 - t410;
t447 = (-Icges(4,6) * t241 + t240 * t308) * t238 + (-Icges(4,5) * t241 + t240 * t312) * t237;
t446 = -(Icges(4,6) * t240 + t241 * t308) * t238 - (Icges(4,5) * t240 + t241 * t312) * t237;
t441 = t482 * t448 + (t525 * t238 + (t513 * t243 + t514 * t242 + (t242 * t478 - t243 * t510) * qJD(4)) * t237 + ((-t299 - t300) * t238 - t531 * t237 + t476) * qJD(3)) * t238;
t440 = (t478 + t521) * t360 + (t218 * t480 + t202 - t402 - t512) * t224 + (t220 * t480 + t203 - t401 - t511) * t223;
t439 = (Icges(6,1) * t392 - t211 - t232 - t510) * t360 + (-t217 * t524 + t200 - t408 + t493) * t224 + (-t219 * t524 + t201 - t407 + t492) * t223;
t438 = t477 * t237;
t354 = qJD(5) * t242;
t231 = t237 * t354;
t228 = pkin(3) * t238 + pkin(6) * t237;
t204 = t228 * t240;
t205 = t228 * t241;
t342 = t204 * t365 + t205 * t364 + qJD(1);
t381 = rSges(6,2) * t393 + t219 * t463 + t220 * t473;
t382 = rSges(6,2) * t394 + t217 * t463 + t218 * t473;
t33 = t223 * t382 - t224 * t381 + t231 + t342;
t355 = qJD(5) * t219;
t423 = rSges(6,2) * t348 - t146 * t463 + t147 * t473 + t355;
t376 = -t192 * t365 - t193 * t364;
t356 = qJD(5) * t217;
t424 = rSges(6,2) * t349 - t144 * t463 + t145 * t473 + t356;
t5 = qJD(5) * t346 - t423 * t224 + t424 * t223 + (t354 + (-t240 * t381 + t241 * t382) * qJD(4)) * t366 + t376;
t436 = t33 * t423 + t381 * t5;
t431 = -t223 / 0.2e1;
t429 = -t224 / 0.2e1;
t166 = t237 * rSges(6,2) + t238 * t327;
t213 = (-rSges(6,1) * t242 + rSges(6,3) * t243) * t237;
t422 = t231 + t452 * qJ(5) + (-t237 * t358 + t238 * t362) * pkin(4) + qJD(3) * t166 + qJD(4) * t213;
t222 = qJD(3) * t228;
t167 = t237 * rSges(5,3) + t238 * t328;
t214 = (-rSges(5,1) * t242 - rSges(5,2) * t243) * t237;
t91 = qJD(3) * t167 + qJD(4) * t214;
t411 = -t222 - t91;
t388 = t238 * t240;
t387 = t238 * t241;
t380 = -t217 * t473 + t218 * t463;
t379 = t219 * t473 - t220 * t463;
t378 = t454 * t240;
t377 = t454 * t241;
t373 = -t165 - t227;
t370 = t240 * t204 + t241 * t205;
t369 = (-pkin(4) * t242 + qJ(5) * t243) * t237 + t213;
t368 = qJD(2) * t241;
t353 = qJD(3) * qJD(4);
t352 = -t222 - t422;
t351 = -t227 + t454;
t340 = t364 / 0.2e1;
t338 = -t360 / 0.2e1;
t337 = t360 / 0.2e1;
t336 = t353 / 0.2e1;
t331 = t237 * t336;
t330 = t238 * t336;
t236 = qJD(2) * t240;
t329 = -t227 * t364 + t236;
t107 = rSges(5,1) * t220 - rSges(5,2) * t219 + rSges(5,3) * t393;
t80 = rSges(5,1) * t147 + rSges(5,2) * t146 + rSges(5,3) * t348;
t34 = -t80 * t360 - t223 * t91 + (-t222 * t240 + (t107 * t237 - t165 * t387) * qJD(4)) * qJD(3);
t105 = rSges(5,1) * t218 - rSges(5,2) * t217 + rSges(5,3) * t394;
t78 = rSges(5,1) * t145 + rSges(5,2) * t144 + rSges(5,3) * t349;
t35 = t78 * t360 + t224 * t91 + (-t222 * t241 + (-t105 * t237 + t165 * t388) * qJD(4)) * qJD(3);
t323 = t240 * t35 - t241 * t34;
t304 = -Icges(4,5) * t237 - Icges(4,6) * t238;
t301 = t105 * t241 - t107 * t240;
t296 = t240 * t330;
t295 = t241 * t330;
t280 = -t227 * t365 - t368;
t46 = t105 * t223 - t107 * t224 + t342;
t277 = t46 * t301;
t274 = qJD(3) * t304;
t273 = t33 * t424 + t382 * t5;
t36 = -t224 * t454 + t360 * t382 + t329 + t355;
t37 = t223 * t454 - t360 * t381 + t280 + t356;
t272 = -t36 * t382 + t37 * t381;
t259 = -qJD(3) * t447 + t240 * t449;
t258 = qJD(3) * t446 + t241 * t449;
t247 = (t33 * t382 + t37 * t454) * t241 + (-t33 * t381 - t36 * t454) * t240;
t246 = t240 * t446 + t241 * t447;
t195 = t304 * t241;
t194 = t304 * t240;
t185 = t241 * t274;
t184 = t240 * t274;
t175 = -t225 * t365 - t368;
t174 = -t225 * t364 + t236;
t126 = -rSges(5,1) * t219 - rSges(5,2) * t220;
t122 = -rSges(5,1) * t217 - rSges(5,2) * t218;
t82 = t462 * qJD(3);
t60 = -t107 * t360 + t280 - t489;
t59 = t105 * t360 + t329 + t488;
t26 = t238 * t301 * t353 + t223 * t78 - t224 * t80 + t376;
t15 = -t222 * t364 - qJD(5) * t146 + t422 * t224 + (t424 * t238 + (-t237 * t382 - t388 * t454) * qJD(3)) * qJD(4);
t14 = -t222 * t365 - qJD(5) * t144 - t422 * t223 + (-t423 * t238 + (t237 * t381 + t387 * t454) * qJD(3)) * qJD(4);
t1 = [-m(4) * t82 + m(5) * t26 + m(6) * t5; m(5) * t323 + m(6) * (-t14 * t241 + t15 * t240); -(t240 * (-t185 * t241 + t240 * t258) - t241 * (-t184 * t241 + t240 * t259)) * t364 + (t194 * qJD(3) * t432 + (-t241 * t195 + t246) * t365) * t340 - (t195 * qJD(3) * t433 + (-t240 * t194 + t246) * t364) * t365 / 0.2e1 + (t240 * (t185 * t240 + t241 * t258) - t241 * (t184 * t240 + t241 * t259)) * t365 + (((t219 * t457 + t220 * t456 + t519) * t238 + t504) * qJD(4) + (((t450 + t464) * qJD(4) + t453) * t238 + t455) * t241 + (t219 * t461 + t220 * t459) * t224 + (t219 * t460 + t220 * t458) * t223) * t431 + (t240 * t467 - t241 * t468) * t430 + (((t217 * t457 + t218 * t456 + t520) * t238 + t503) * qJD(4) + (((t450 + t500) * qJD(4) + t453) * t238 + t455) * t240 + (t217 * t461 + t218 * t459) * t224 + (t217 * t460 + t218 * t458) * t223) * t429 + (t240 * t469 - t241 * t470) * t428 + t240 * t483 - t475 * t241 / 0.2e1 - t466 * t361 / 0.2e1 + (((t129 * t242 - t137 * t243 + t97) * t223 + (t128 * t242 - t136 * t243 + t96) * t224 + t61 * qJD(4)) * t237 + ((-t284 * t238 + (-t242 * t151 - t159 * t243 - t154) * t237 + t318) * qJD(4) + t434) * t238 + ((t135 * t242 - t139 * t243 + t95) * t223 + (t134 * t242 - t138 * t243 + t94) * t224 + t62 * qJD(4)) * t237 + ((t283 * t238 + (t242 * t157 - t161 * t243 - t152) * t237 + t317) * qJD(4) + t435) * t238) * t337 + (t240 * t496 - t241 * t497) * t331 + (t240 * t499 - t241 * t500) * t296 + (t464 * t240 - t241 * t498) * t295 + (t5 * t370 + (t15 * t351 + t352 * t36 + t436) * t241 + (t14 * t351 + t352 * t37 + t273) * t240 - (t272 * t237 + (t36 * t378 - t37 * t377 + t247) * t238) * qJD(4) - (t240 * t37 + t241 * t36) * (-t231 - t222) + (t37 * t223 - t36 * t224) * (t326 * t238 + t166) + (-t378 * t223 + t377 * t224 - t238 * t354 + t451) * t33) * m(6) + (-t59 * (t167 * t224 - t228 * t364) - t60 * (-t167 * t223 - t228 * t365) - ((-t105 * t59 + t107 * t60) * t237 + t277 * t238) * qJD(4) + t26 * t370 + t451 * t46 + (t26 * t107 + t35 * t373 + t411 * t59 + (-t488 + t80) * t46) * t241 + (t26 * t105 + t34 * t373 + t411 * t60 + (t489 + t78) * t46) * t240) * m(5) + (-t82 * t437 + (-t174 * t241 - t175 * t240 + t462) * qJD(3) + t174 * t364 + t175 * t365) * (rSges(4,1) * t238 - rSges(4,2) * t237) * m(4) + ((-t445 + t471) * t241 + (t444 + t472) * t240) * t338; (t219 * t440 + t220 * t439 - t241 * t438) * t431 + ((t240 * t468 + t241 * t467) * t237 + t485) * t430 + (t217 * t440 + t218 * t439 - t240 * t438) * t429 + ((t240 * t470 + t241 * t469) * t237 + t484) * t428 - (qJD(4) * t441 + t223 * t444 + t224 * t445) * t238 / 0.2e1 + t475 * t394 / 0.2e1 + t393 * t483 + t466 * t448 / 0.2e1 + ((t240 * t445 + t241 * t444) * t237 + t441) * t338 + (t477 * t238 + (t242 * t440 + t243 * t439) * t237) * t337 + t472 * t349 / 0.2e1 + t471 * t238 * t340 + (t237 * t476 - t238 * t482) * t331 + (t505 * t237 - t238 * t529) * t296 + (t506 * t237 - t238 * t528) * t295 + ((qJD(3) * t247 - t14 * t381 + t15 * t382 + t36 * t424 - t37 * t423) * t238 + (t272 * qJD(3) + (t14 * t454 - t37 * t422 + t273) * t241 + (-t15 * t454 + t36 * t422 - t436) * t240) * t237 - (t218 * t37 + t220 * t36 + t33 * t391) * qJD(5) - (t33 * t379 + t36 * t369) * t224 - (t33 * t380 - t369 * t37) * t223 - (t36 * t380 + t37 * t379) * t360) * m(6) + ((t35 * t105 - t34 * t107 + t59 * t78 - t60 * t80 + (t277 + (t240 * t59 - t241 * t60) * t165) * qJD(3)) * t238 + (t59 * (-qJD(3) * t105 + t240 * t91) + t60 * (qJD(3) * t107 - t241 * t91) + t26 * t301 + t46 * (-t240 * t80 + t241 * t78) + t323 * t165) * t237 - t59 * (t122 * t360 + t214 * t224) - t60 * (-t126 * t360 - t214 * t223) - t46 * (t122 * t223 - t126 * t224)) * m(5); (t5 * t392 + t14 * t217 + t15 * t219 + (t219 * t360 + t223 * t392 - t144) * t37 + (-t217 * t360 - t224 * t392 - t146) * t36 + (-t217 * t223 + t219 * t224 + t452) * t33) * m(6);];
tauc = t1(:);
