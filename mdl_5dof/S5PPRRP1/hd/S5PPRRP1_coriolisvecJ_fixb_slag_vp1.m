% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRP1
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:39
% EndTime: 2019-12-05 15:07:08
% DurationCPUTime: 24.40s
% Computational Cost: add. (14471->565), mult. (20801->830), div. (0->0), fcn. (20013->6), ass. (0->298)
t533 = Icges(5,5) + Icges(6,5);
t506 = Icges(5,6) + Icges(6,6);
t235 = pkin(8) + qJ(3);
t234 = cos(t235);
t237 = cos(pkin(7));
t240 = cos(qJ(4));
t382 = t237 * t240;
t236 = sin(pkin(7));
t239 = sin(qJ(4));
t385 = t236 * t239;
t216 = -t234 * t385 - t382;
t383 = t237 * t239;
t384 = t236 * t240;
t217 = t234 * t384 - t383;
t233 = sin(t235);
t391 = t233 * t236;
t94 = Icges(6,5) * t217 + Icges(6,6) * t216 + Icges(6,3) * t391;
t96 = Icges(5,5) * t217 + Icges(5,6) * t216 + Icges(5,3) * t391;
t497 = t94 + t96;
t218 = -t234 * t383 + t384;
t219 = t234 * t382 + t385;
t390 = t233 * t237;
t95 = Icges(6,5) * t219 + Icges(6,6) * t218 + Icges(6,3) * t390;
t97 = Icges(5,5) * t219 + Icges(5,6) * t218 + Icges(5,3) * t390;
t496 = t95 + t97;
t405 = Icges(5,4) * t217;
t100 = Icges(5,2) * t216 + Icges(5,6) * t391 + t405;
t401 = Icges(6,4) * t217;
t98 = Icges(6,2) * t216 + Icges(6,6) * t391 + t401;
t519 = t100 + t98;
t404 = Icges(5,4) * t219;
t101 = Icges(5,2) * t218 + Icges(5,6) * t390 + t404;
t400 = Icges(6,4) * t219;
t99 = Icges(6,2) * t218 + Icges(6,6) * t390 + t400;
t518 = t101 + t99;
t200 = Icges(6,4) * t216;
t102 = Icges(6,1) * t217 + Icges(6,5) * t391 + t200;
t202 = Icges(5,4) * t216;
t104 = Icges(5,1) * t217 + Icges(5,5) * t391 + t202;
t517 = t102 + t104;
t201 = Icges(6,4) * t218;
t103 = Icges(6,1) * t219 + Icges(6,5) * t390 + t201;
t203 = Icges(5,4) * t218;
t105 = Icges(5,1) * t219 + Icges(5,5) * t390 + t203;
t516 = t103 + t105;
t484 = Icges(5,1) + Icges(6,1);
t532 = Icges(5,2) + Icges(6,2);
t541 = Icges(5,3) + Icges(6,3);
t302 = Icges(6,5) * t240 - Icges(6,6) * t239;
t152 = -Icges(6,3) * t234 + t233 * t302;
t303 = Icges(5,5) * t240 - Icges(5,6) * t239;
t154 = -Icges(5,3) * t234 + t233 * t303;
t531 = -t154 - t152;
t398 = Icges(6,4) * t240;
t305 = -Icges(6,2) * t239 + t398;
t156 = -Icges(6,6) * t234 + t233 * t305;
t402 = Icges(5,4) * t240;
t306 = -Icges(5,2) * t239 + t402;
t158 = -Icges(5,6) * t234 + t233 * t306;
t483 = t156 + t158;
t399 = Icges(6,4) * t239;
t309 = Icges(6,1) * t240 - t399;
t160 = -Icges(6,5) * t234 + t233 * t309;
t403 = Icges(5,4) * t239;
t310 = Icges(5,1) * t240 - t403;
t162 = -Icges(5,5) * t234 + t233 * t310;
t489 = t160 + t162;
t540 = (-t533 * t239 - t506 * t240) * t233;
t503 = t518 * t216 + t516 * t217 + t496 * t391;
t502 = t519 * t218 + t517 * t219 + t497 * t390;
t504 = t519 * t216 + t517 * t217 + t391 * t497;
t466 = t218 * t518 + t219 * t516 + t390 * t496;
t539 = t216 * t483 + t217 * t489 - t391 * t531;
t538 = t218 * t483 + t219 * t489 - t390 * t531;
t362 = qJD(3) * t233;
t349 = t239 * t362;
t144 = -qJD(4) * t217 + t236 * t349;
t348 = t240 * t362;
t145 = qJD(4) * t216 - t236 * t348;
t360 = qJD(3) * t236;
t347 = t234 * t360;
t537 = t144 * t506 + t145 * t533 + t347 * t541;
t146 = -qJD(4) * t219 + t237 * t349;
t147 = qJD(4) * t218 - t237 * t348;
t359 = qJD(3) * t237;
t346 = t234 * t359;
t536 = t146 * t506 + t147 * t533 + t346 * t541;
t153 = Icges(6,3) * t233 + t234 * t302;
t155 = Icges(5,3) * t233 + t234 * t303;
t535 = t540 * qJD(4) + (t153 + t155) * qJD(3);
t534 = Icges(5,4) + Icges(6,4);
t157 = Icges(6,6) * t233 + t234 * t305;
t159 = Icges(5,6) * t233 + t234 * t306;
t530 = t157 + t159;
t161 = Icges(6,5) * t233 + t234 * t309;
t163 = Icges(5,5) * t233 + t234 * t310;
t458 = -t161 - t163;
t529 = (-t240 * t532 - t399 - t403) * t233;
t528 = (t239 * t484 + t398 + t402) * t233;
t527 = t503 * t237;
t526 = t502 * t236;
t525 = t144 * t532 + t145 * t534 + t347 * t506;
t524 = t146 * t532 + t147 * t534 + t346 * t506;
t523 = t534 * t144 + t145 * t484 + t533 * t347;
t522 = t534 * t146 + t147 * t484 + t533 * t346;
t521 = qJD(3) * t530 + qJD(4) * t529;
t520 = qJD(3) * t458 + qJD(4) * t528;
t361 = qJD(3) * t234;
t515 = -t535 * t233 + t531 * t361;
t514 = t536 * t233 + t496 * t361;
t513 = t537 * t233 + t497 * t361;
t512 = t466 * t237 + t526;
t511 = t504 * t236 + t527;
t510 = t538 * t233;
t509 = t539 * t233;
t358 = qJD(4) * t233;
t344 = t237 * t358;
t222 = t344 + t360;
t345 = t236 * t358;
t223 = t345 - t359;
t505 = t222 * t236 - t223 * t237;
t432 = t222 / 0.2e1;
t430 = t223 / 0.2e1;
t472 = t519 * t144 + t517 * t145 + t216 * t525 + t523 * t217 + t513 * t236;
t471 = t144 * t518 + t145 * t516 + t216 * t524 + t217 * t522 + t236 * t514;
t470 = t519 * t146 + t517 * t147 + t218 * t525 + t523 * t219 + t513 * t237;
t469 = t146 * t518 + t147 * t516 + t218 * t524 + t219 * t522 + t237 * t514;
t314 = t102 * t240 - t239 * t98;
t47 = t233 * t314 - t234 * t94;
t301 = -t100 * t239 + t104 * t240;
t49 = t233 * t301 - t234 * t96;
t501 = t47 + t49;
t313 = t103 * t240 - t239 * t99;
t48 = t233 * t313 - t234 * t95;
t300 = -t101 * t239 + t105 * t240;
t50 = t233 * t300 - t234 * t97;
t500 = t48 + t50;
t298 = -t156 * t239 + t160 * t240;
t389 = t234 * t152;
t61 = t233 * t298 - t389;
t297 = -t158 * t239 + t162 * t240;
t388 = t234 * t154;
t62 = t233 * t297 - t388;
t485 = t61 + t62;
t488 = (-t146 * t483 - t147 * t489 - t218 * t521 + t219 * t520 + t237 * t515) * t234 + (t234 * t512 + t510) * qJD(3);
t487 = (-t144 * t483 - t145 * t489 - t216 * t521 + t217 * t520 + t236 * t515) * t234 + (t234 * t511 + t509) * qJD(3);
t357 = qJD(4) * t234;
t444 = t222 * (-t219 * t532 + t201 + t203 + t516) + t223 * (-t217 * t532 + t200 + t202 + t517) - t357 * (t489 + t529);
t486 = t470 * t430 + t469 * t432 + t488 * qJD(4) / 0.2e1;
t482 = t540 * t357 + (-t216 * t533 + t506 * t217) * t223 + (-t218 * t533 + t506 * t219) * t222;
t315 = t49 * t236 + t50 * t237;
t316 = t47 * t236 + t48 * t237;
t481 = t315 + t316;
t436 = t237 ^ 2;
t437 = t236 ^ 2;
t441 = t436 + t437;
t477 = qJD(4) * t487 + t222 * t471 + t223 * t472;
t449 = ((t301 + t314) * qJD(3) - t537) * t234 + (t523 * t240 - t525 * t239 + (-t239 * t517 - t240 * t519) * qJD(4) + t497 * qJD(3)) * t233;
t448 = ((t300 + t313) * qJD(3) - t536) * t234 + (t522 * t240 - t524 * t239 + (-t239 * t516 - t240 * t518) * qJD(4) + t496 * qJD(3)) * t233;
t475 = rSges(6,1) + pkin(4);
t474 = t222 * t503 + t223 * t504 - t357 * t539;
t473 = t222 * t466 + t223 * t502 - t357 * t538;
t468 = t222 * t500 + t223 * t501 - t357 * t485;
t224 = rSges(4,1) * t233 + rSges(4,2) * t234;
t465 = qJD(3) * t224 * t441;
t132 = t156 * t236;
t134 = t158 * t236;
t464 = -t132 - t134;
t133 = t156 * t237;
t135 = t158 * t237;
t463 = -t133 - t135;
t136 = t160 * t236;
t138 = t162 * t236;
t462 = -t136 - t138;
t137 = t160 * t237;
t139 = t162 * t237;
t461 = -t137 - t139;
t426 = pkin(4) * t240;
t149 = qJ(5) * t233 + t234 * t426;
t324 = rSges(6,1) * t240 - rSges(6,2) * t239;
t460 = t233 * rSges(6,3) + t234 * t324 + t149;
t280 = t155 - t297;
t281 = t153 - t298;
t438 = (t152 * t237 + t313) * t222 + (t152 * t236 + t314) * t223;
t439 = -(-t154 * t237 - t300) * t222 - (-t154 * t236 - t301) * t223;
t457 = (-t438 - t439 + (-t280 - t281) * t357) * t233;
t270 = qJ(5) * t234 - t233 * t426;
t456 = t234 * rSges(6,3) - t233 * t324 + t270;
t455 = t222 * t496 + t223 * t497;
t226 = pkin(3) * t233 - pkin(6) * t234;
t279 = qJD(3) * t226;
t192 = t236 * t279;
t193 = t237 * t279;
t454 = -t236 * t192 - t237 * t193 + t441 * t279;
t453 = -t388 - t389;
t406 = Icges(4,4) * t234;
t407 = Icges(4,4) * t233;
t452 = -(-Icges(4,2) * t234 - t407) * t362 + (-Icges(4,1) * t233 - t406) * t361;
t308 = -Icges(4,2) * t233 + t406;
t312 = Icges(4,1) * t234 - t407;
t451 = (-Icges(4,6) * t237 + t236 * t308) * t234 + (-Icges(4,5) * t237 + t236 * t312) * t233;
t450 = -(Icges(4,6) * t236 + t237 * t308) * t234 - (Icges(4,5) * t236 + t237 * t312) * t233;
t227 = pkin(3) * t234 + pkin(6) * t233;
t221 = qJD(3) * t227;
t355 = qJD(5) * t234;
t334 = -t221 + t355;
t445 = t485 * t362 + (t535 * t234 + (t520 * t240 + t521 * t239 + (t239 * t489 + t483 * t240) * qJD(4)) * t233 + ((-t297 - t298) * t234 + t531 * t233 + t481) * qJD(3)) * t234;
t443 = (t483 + t528) * t357 + (t216 * t484 - t401 - t405 - t519) * t223 + (t218 * t484 - t400 - t404 - t518) * t222;
t442 = t482 * t233;
t325 = rSges(5,1) * t240 - rSges(5,2) * t239;
t167 = -t234 * rSges(5,3) + t233 * t325;
t206 = t227 * t236;
t207 = t227 * t237;
t341 = t206 * t360 + t207 * t359 + qJD(1);
t409 = rSges(6,1) * t219 + rSges(6,2) * t218 + rSges(6,3) * t390 + pkin(4) * t385 + t149 * t237;
t410 = rSges(6,1) * t217 + rSges(6,2) * t216 + rSges(6,3) * t391 - pkin(4) * t383 + t149 * t236;
t33 = t222 * t410 - t223 * t409 + t341 - t355;
t356 = qJD(5) * t233;
t229 = t237 * t356;
t421 = pkin(4) * qJD(4);
t352 = t239 * t421;
t245 = qJD(3) * t270 - t234 * t352;
t351 = t240 * t421;
t423 = rSges(6,1) * t147 + rSges(6,2) * t146 + rSges(6,3) * t346 + t236 * t351 + t237 * t245 + t229;
t371 = -t192 * t360 - t193 * t359;
t228 = t236 * t356;
t424 = rSges(6,1) * t145 + rSges(6,2) * t144 + rSges(6,3) * t347 + t236 * t245 - t237 * t351 + t228;
t5 = -t423 * t223 + t424 * t222 + (t356 + (-t236 * t409 + t237 * t410) * t357) * qJD(3) + t371;
t440 = t33 * t423 + t409 * t5;
t435 = -m(6) / 0.2e1;
t434 = m(6) / 0.2e1;
t433 = -t222 / 0.2e1;
t431 = -t223 / 0.2e1;
t214 = (-rSges(6,1) * t239 - rSges(6,2) * t240) * t233;
t422 = qJD(3) * t460 + qJD(4) * t214 - t233 * t352 - t355;
t169 = t233 * rSges(5,3) + t234 * t325;
t215 = (-rSges(5,1) * t239 - rSges(5,2) * t240) * t233;
t93 = qJD(3) * t169 + qJD(4) * t215;
t408 = -t221 - t93;
t387 = t234 * t236;
t386 = t234 * t237;
t377 = t456 * t236;
t376 = t456 * t237;
t375 = -rSges(6,2) * t217 + t216 * t475;
t374 = rSges(6,2) * t219 - t218 * t475;
t367 = -t167 - t226;
t365 = t236 * t206 + t237 * t207;
t364 = qJD(2) * t237;
t354 = qJD(3) * qJD(4);
t353 = -t221 - t422;
t350 = -t226 + t456;
t339 = t359 / 0.2e1;
t337 = -t357 / 0.2e1;
t336 = t357 / 0.2e1;
t335 = t354 / 0.2e1;
t329 = pkin(4) * t233 * t239 - t214;
t328 = t233 * t335;
t327 = t234 * t335;
t232 = qJD(2) * t236;
t326 = -t226 * t359 + t232;
t109 = rSges(5,1) * t219 + rSges(5,2) * t218 + rSges(5,3) * t390;
t80 = rSges(5,1) * t147 + rSges(5,2) * t146 + rSges(5,3) * t346;
t34 = -t80 * t357 - t222 * t93 + (-t221 * t236 + (t109 * t233 - t167 * t386) * qJD(4)) * qJD(3);
t107 = rSges(5,1) * t217 + rSges(5,2) * t216 + rSges(5,3) * t391;
t78 = rSges(5,1) * t145 + rSges(5,2) * t144 + rSges(5,3) * t347;
t35 = t78 * t357 + t223 * t93 + (-t221 * t237 + (-t107 * t233 + t167 * t387) * qJD(4)) * qJD(3);
t321 = t236 * t35 - t237 * t34;
t304 = -Icges(4,5) * t233 - Icges(4,6) * t234;
t299 = t107 * t237 - t109 * t236;
t294 = t236 * t327;
t293 = t237 * t327;
t278 = -t226 * t360 - t364;
t46 = t107 * t222 - t109 * t223 + t341;
t276 = t46 * t299;
t273 = qJD(3) * t304;
t272 = t33 * t424 + t410 * t5;
t36 = -t223 * t456 + t357 * t410 + t229 + t326;
t37 = t222 * t456 - t357 * t409 + t228 + t278;
t271 = -t36 * t410 + t37 * t409;
t257 = -qJD(3) * t451 + t236 * t452;
t256 = qJD(3) * t450 + t237 * t452;
t244 = (t33 * t410 + t37 * t456) * t237 + (-t33 * t409 - t36 * t456) * t236;
t243 = t236 * t450 + t237 * t451;
t195 = t304 * t237;
t194 = t304 * t236;
t185 = t237 * t273;
t184 = t236 * t273;
t175 = -t224 * t360 - t364;
t174 = -t224 * t359 + t232;
t127 = rSges(5,1) * t218 - rSges(5,2) * t219;
t125 = rSges(5,1) * t216 - rSges(5,2) * t217;
t83 = t465 * qJD(3);
t60 = -t109 * t357 - t167 * t222 + t278;
t59 = t107 * t357 + t167 * t223 + t326;
t26 = t234 * t299 * t354 + t222 * t78 - t223 * t80 + t371;
t15 = t422 * t223 + t424 * t357 + (t334 * t237 + (-t233 * t410 - t387 * t456) * qJD(4)) * qJD(3);
t14 = -t422 * t222 - t423 * t357 + (t334 * t236 + (t233 * t409 + t386 * t456) * qJD(4)) * qJD(3);
t1 = [-m(4) * t83 + m(5) * t26 + m(6) * t5; m(5) * t321 + m(6) * (-t14 * t237 + t15 * t236); -(t236 * (-t185 * t237 + t236 * t256) - t237 * (-t184 * t237 + t236 * t257)) * t359 + (t194 * qJD(3) * t436 + (-t237 * t195 + t243) * t360) * t339 - (t195 * qJD(3) * t437 + (-t236 * t194 + t243) * t359) * t360 / 0.2e1 + (t236 * (t185 * t236 + t237 * t256) - t237 * (t184 * t236 + t237 * t257)) * t360 + (((-t218 * t530 + t219 * t458 + t526) * t234 + t510) * qJD(4) + (((t453 + t466) * qJD(4) + t455) * t234 + t457) * t237 + (t218 * t464 + t219 * t462) * t223 + (t218 * t463 + t219 * t461) * t222) * t433 + (t236 * t469 - t237 * t470) * t432 + (((-t216 * t530 + t217 * t458 + t527) * t234 + t509) * qJD(4) + (((t453 + t504) * qJD(4) + t455) * t234 + t457) * t236 + (t216 * t464 + t217 * t462) * t223 + (t216 * t463 + t217 * t461) * t222) * t431 + (t236 * t471 - t237 * t472) * t430 + t236 * t486 - t477 * t237 / 0.2e1 - t468 * t358 / 0.2e1 + (((t133 * t239 - t137 * t240 + t95) * t222 + (t132 * t239 - t136 * t240 + t94) * t223 + t61 * qJD(4)) * t233 + ((t281 * t234 + (t157 * t239 - t161 * t240 - t152) * t233 + t316) * qJD(4) + t438) * t234 + ((t135 * t239 - t139 * t240 + t97) * t222 + (t134 * t239 - t138 * t240 + t96) * t223 + t62 * qJD(4)) * t233 + ((t280 * t234 + (t159 * t239 - t163 * t240 - t154) * t233 + t315) * qJD(4) + t439) * t234) * t336 + (t236 * t500 - t237 * t501) * t328 + (t236 * t503 - t237 * t504) * t294 + (t466 * t236 - t237 * t502) * t293 + (t5 * t365 + (t15 * t350 + t353 * t36 + t440) * t237 + (t14 * t350 + t353 * t37 + t272) * t236 - (t271 * t233 + (t36 * t377 - t37 * t376 + t244) * t234) * qJD(4) - (t236 * t37 + t237 * t36) * t334 + (t37 * t222 - t36 * t223) * t460 + (-t377 * t222 + t376 * t223 - t356 + t454) * t33) * m(6) + (-t59 * (t169 * t223 - t227 * t359) - t60 * (-t169 * t222 - t227 * t360) - ((-t107 * t59 + t109 * t60) * t233 + t276 * t234) * qJD(4) + t26 * t365 + (t26 * t109 + t35 * t367 + t408 * t59) * t237 + (t26 * t107 + t34 * t367 + t408 * t60) * t236 + (t167 * t505 + t78 * t236 + t80 * t237 + t454) * t46) * m(5) + (-t83 * t441 + (-t174 * t237 - t175 * t236 + t465) * qJD(3) + t174 * t359 + t175 * t360) * (rSges(4,1) * t234 - rSges(4,2) * t233) * m(4) + ((-t449 + t473) * t237 + (t448 + t474) * t236) * t337; (t218 * t444 + t219 * t443 - t237 * t442) * t433 + ((t236 * t470 + t237 * t469) * t233 + t488) * t432 + (t216 * t444 + t217 * t443 - t236 * t442) * t431 + ((t236 * t472 + t237 * t471) * t233 + t487) * t430 - (qJD(4) * t445 + t222 * t448 + t223 * t449) * t234 / 0.2e1 + t477 * t391 / 0.2e1 + t390 * t486 + t468 * t362 / 0.2e1 + ((t236 * t449 + t237 * t448) * t233 + t445) * t337 + (t482 * t234 + (-t239 * t444 + t443 * t240) * t233) * t336 + t474 * t347 / 0.2e1 + t473 * t234 * t339 + (t233 * t481 - t234 * t485) * t328 + (t233 * t511 - t234 * t539) * t294 + (t233 * t512 - t234 * t538) * t293 + ((qJD(3) * t244 - t14 * t409 + t15 * t410 + t36 * t424 - t37 * t423) * t234 + (t271 * qJD(3) + (t14 * t456 - t37 * t422 + t272) * t237 + (-t15 * t456 + t36 * t422 - t440) * t236) * t233 - (-t329 * t36 + t33 * t374) * t223 - (t329 * t37 + t33 * t375) * t222 - (t36 * t375 + t37 * t374) * t357) * m(6) + (-t59 * (t125 * t357 + t215 * t223) - t60 * (-t127 * t357 - t215 * t222) - t46 * (t125 * t222 - t127 * t223) + (t35 * t107 - t34 * t109 + t59 * t78 - t60 * t80 + (t276 + (t236 * t59 - t237 * t60) * t167) * qJD(3)) * t234 + (t59 * (-qJD(3) * t107 + t236 * t93) + t60 * (qJD(3) * t109 - t237 * t93) + t26 * t299 + t46 * (-t236 * t80 + t237 * t78) + t321 * t167) * t233) * m(5); 0.2e1 * ((qJD(3) * t33 + t14 * t236 + t15 * t237) * t434 + t33 * t505 * t435) * t233 + 0.2e1 * ((t359 * t36 + t360 * t37 - t5) * t434 + (t36 * (-t223 + t345) + t37 * (t222 - t344)) * t435) * t234;];
tauc = t1(:);
