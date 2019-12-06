% Calculate vector of inverse dynamics joint torques for
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:46
% EndTime: 2019-12-05 15:45:12
% DurationCPUTime: 18.80s
% Computational Cost: add. (18214->672), mult. (18527->1009), div. (0->0), fcn. (17292->10), ass. (0->363)
t332 = qJ(2) + pkin(9);
t324 = sin(t332);
t325 = cos(t332);
t337 = sin(qJ(2));
t339 = cos(qJ(2));
t598 = -Icges(3,5) * t337 - Icges(4,5) * t324 - Icges(3,6) * t339 - Icges(4,6) * t325;
t602 = Icges(3,3) + Icges(4,3);
t601 = t598 * qJD(2);
t600 = Icges(3,5) * t339 + Icges(4,5) * t325 - Icges(3,6) * t337 - Icges(4,6) * t324;
t536 = Icges(4,4) * t325;
t537 = Icges(4,4) * t324;
t538 = Icges(3,4) * t339;
t539 = Icges(3,4) * t337;
t599 = (-t324 * (-Icges(4,2) * t325 - t537) + t325 * (-Icges(4,1) * t324 - t536) - t337 * (-Icges(3,2) * t339 - t539) + t339 * (-Icges(3,1) * t337 - t538)) * qJD(2);
t334 = cos(pkin(8));
t330 = t334 ^ 2;
t306 = rSges(3,1) * t337 + rSges(3,2) * t339;
t333 = sin(pkin(8));
t329 = t333 ^ 2;
t565 = t329 + t330;
t597 = t306 * t565;
t596 = t600 * t333 - t602 * t334;
t595 = t602 * t333 + t600 * t334;
t594 = t601 * t333;
t593 = t601 * t334;
t436 = -Icges(4,2) * t324 + t536;
t206 = Icges(4,6) * t333 + t334 * t436;
t444 = Icges(4,1) * t325 - t537;
t208 = Icges(4,5) * t333 + t334 * t444;
t438 = -Icges(3,2) * t337 + t538;
t244 = Icges(3,6) * t333 + t334 * t438;
t446 = Icges(3,1) * t339 - t539;
t246 = Icges(3,5) * t333 + t334 * t446;
t592 = t599 * t334 + (-t206 * t325 - t208 * t324 - t244 * t339 - t246 * t337) * qJD(2);
t205 = -Icges(4,6) * t334 + t333 * t436;
t207 = -Icges(4,5) * t334 + t333 * t444;
t243 = -Icges(3,6) * t334 + t333 * t438;
t245 = -Icges(3,5) * t334 + t333 * t446;
t591 = -t599 * t333 + (t205 * t325 + t207 * t324 + t243 * t339 + t245 * t337) * qJD(2);
t590 = -t206 * t324 + t208 * t325 - t244 * t337 + t246 * t339;
t589 = t205 * t324 - t207 * t325 + t243 * t337 - t245 * t339;
t326 = qJ(4) + t332;
t315 = sin(t326);
t316 = cos(t326);
t534 = Icges(5,4) * t316;
t434 = -Icges(5,2) * t315 + t534;
t190 = -Icges(5,6) * t334 + t333 * t434;
t441 = -Icges(5,1) * t315 - t534;
t588 = t441 * t333 - t190;
t191 = Icges(5,6) * t333 + t334 * t434;
t587 = t441 * t334 - t191;
t535 = Icges(5,4) * t315;
t442 = Icges(5,1) * t316 - t535;
t192 = -Icges(5,5) * t334 + t333 * t442;
t433 = -Icges(5,2) * t316 - t535;
t586 = -t433 * t333 - t192;
t193 = Icges(5,5) * t333 + t334 * t442;
t585 = -t433 * t334 - t193;
t584 = t598 * t333;
t583 = t598 * t334;
t338 = cos(qJ(5));
t516 = t334 * t338;
t336 = sin(qJ(5));
t519 = t333 * t336;
t255 = -t316 * t519 - t516;
t517 = t334 * t336;
t518 = t333 * t338;
t256 = t316 * t518 - t517;
t525 = t315 * t333;
t110 = rSges(6,1) * t256 + rSges(6,2) * t255 + rSges(6,3) * t525;
t257 = -t316 * t517 + t518;
t258 = t316 * t516 + t519;
t524 = t315 * t334;
t111 = rSges(6,1) * t258 + rSges(6,2) * t257 + rSges(6,3) * t524;
t484 = t315 * t518;
t485 = t315 * t519;
t522 = t316 * t333;
t504 = rSges(6,2) * t485 + rSges(6,3) * t522;
t140 = -rSges(6,1) * t484 + t504;
t482 = t315 * t516;
t483 = t315 * t517;
t521 = t316 * t334;
t503 = rSges(6,2) * t483 + rSges(6,3) * t521;
t141 = -rSges(6,1) * t482 + t503;
t547 = rSges(6,1) * t338;
t455 = -rSges(6,2) * t336 + t547;
t178 = -rSges(6,3) * t316 + t315 * t455;
t179 = t315 * rSges(6,3) + t316 * t455;
t323 = qJD(2) * t333;
t298 = qJD(4) * t333 + t323;
t495 = qJD(5) * t315;
t229 = t334 * t495 + t298;
t492 = qJD(5) * t333;
t331 = qJD(2) + qJD(4);
t520 = t331 * t334;
t230 = t315 * t492 - t520;
t476 = t316 * t492;
t494 = qJD(5) * t316;
t274 = pkin(4) * t315 - pkin(7) * t316;
t552 = pkin(2) * t337;
t295 = -pkin(3) * t324 - t552;
t288 = t295 * qJD(2);
t496 = qJD(3) * t334;
t360 = t288 * t333 - t496;
t54 = -t111 * t494 - t178 * t229 - t274 * t298 + t360;
t322 = qJD(3) * t333;
t367 = t288 * t334 + t322;
t477 = t110 * t494;
t55 = t178 * t230 - t274 * t520 + t367 + t477;
t564 = g(1) * t334 + g(2) * t333;
t567 = t316 * pkin(4) + t315 * pkin(7);
t582 = -t54 * (-t179 * t229 - t567 * t298 + t111 * t495 + (-t178 * t334 - t141) * t494) - t55 * (-t110 * t495 + t140 * t494 + t178 * t476 + t230 * t179 - t520 * t567) - t564 * (-pkin(4) - t547) * t315;
t104 = Icges(6,5) * t256 + Icges(6,6) * t255 + Icges(6,3) * t525;
t533 = Icges(6,4) * t256;
t106 = Icges(6,2) * t255 + Icges(6,6) * t525 + t533;
t238 = Icges(6,4) * t255;
t108 = Icges(6,1) * t256 + Icges(6,5) * t525 + t238;
t44 = t104 * t525 + t106 * t255 + t108 * t256;
t105 = Icges(6,5) * t258 + Icges(6,6) * t257 + Icges(6,3) * t524;
t532 = Icges(6,4) * t258;
t107 = Icges(6,2) * t257 + Icges(6,6) * t524 + t532;
t239 = Icges(6,4) * t257;
t109 = Icges(6,1) * t258 + Icges(6,5) * t524 + t239;
t45 = t105 * t525 + t107 * t255 + t109 * t256;
t46 = t104 * t524 + t106 * t257 + t108 * t258;
t47 = t105 * t524 + t107 * t257 + t109 * t258;
t424 = Icges(6,5) * t338 - Icges(6,6) * t336;
t166 = -Icges(6,3) * t316 + t315 * t424;
t530 = Icges(6,4) * t338;
t432 = -Icges(6,2) * t336 + t530;
t168 = -Icges(6,6) * t316 + t315 * t432;
t531 = Icges(6,4) * t336;
t440 = Icges(6,1) * t338 - t531;
t170 = -Icges(6,5) * t316 + t315 * t440;
t60 = t166 * t525 + t168 * t255 + t170 * t256;
t61 = t166 * t524 + t168 * t257 + t170 * t258;
t579 = (t229 * t47 + t230 * t46 - t494 * t61) * t334 + (t229 * t45 + t230 * t44 - t494 * t60) * t333;
t327 = t339 * pkin(2);
t566 = t325 * rSges(4,1) - rSges(4,2) * t324;
t569 = t566 + t327;
t568 = t316 * rSges(5,1) - rSges(5,2) * t315;
t420 = -t107 * t336 + t109 * t338;
t421 = -t106 * t336 + t108 * t338;
t563 = -(-t166 * t334 - t420) * t229 - (-t166 * t333 - t421) * t230;
t431 = -Icges(6,2) * t338 - t531;
t349 = t229 * (-Icges(6,2) * t258 + t109 + t239) + t230 * (-Icges(6,2) * t256 + t108 + t238) - t494 * (t431 * t315 + t170);
t340 = qJD(2) ^ 2;
t320 = qJDD(2) * t333;
t296 = qJDD(4) * t333 + t320;
t493 = qJD(5) * t331;
t390 = qJDD(5) * t315 + t316 * t493;
t146 = t334 * t390 + t296;
t561 = t146 / 0.2e1;
t297 = (-qJDD(2) - qJDD(4)) * t334;
t147 = t333 * t390 + t297;
t560 = t147 / 0.2e1;
t559 = -t229 / 0.2e1;
t558 = t229 / 0.2e1;
t557 = -t230 / 0.2e1;
t556 = t230 / 0.2e1;
t240 = -qJDD(5) * t316 + t315 * t493;
t555 = t240 / 0.2e1;
t554 = t333 / 0.2e1;
t553 = -t334 / 0.2e1;
t313 = pkin(3) * t325;
t543 = pkin(2) * qJD(2);
t542 = t333 * t46;
t541 = t334 * t45;
t218 = t567 * t331;
t454 = -rSges(6,1) * t336 - rSges(6,2) * t338;
t523 = t316 * t331;
t99 = t455 * t523 + (rSges(6,3) * t331 + qJD(5) * t454) * t315;
t540 = -t218 - t99;
t526 = t166 * t316;
t272 = rSges(5,1) * t315 + rSges(5,2) * t316;
t227 = t272 * t333;
t180 = t331 * t227;
t228 = t272 * t334;
t181 = t331 * t228;
t513 = -t333 * t180 - t334 * t181;
t195 = -rSges(5,3) * t334 + t333 * t568;
t196 = rSges(5,3) * t333 + t334 * t568;
t512 = t333 * t195 + t334 * t196;
t510 = -t178 - t274;
t199 = -qJ(3) * t334 + t327 * t333;
t200 = qJ(3) * t333 + t327 * t334;
t509 = t333 * t199 + t334 * t200;
t284 = t333 * t295;
t285 = t334 * t295;
t497 = qJD(2) * t334;
t508 = (t333 * t552 + t284) * t323 + (t334 * t552 + t285) * t497;
t487 = t337 * t543;
t286 = -t333 * t487 - t496;
t287 = -t334 * t487 + t322;
t505 = t333 * t286 + t334 * t287;
t501 = t313 + t327;
t502 = t501 - t327;
t490 = -m(4) - m(5) - m(6);
t489 = qJDD(2) * t334;
t488 = qJDD(3) * t334;
t486 = t339 * t543;
t481 = t331 * t522;
t480 = t316 * t520;
t302 = pkin(7) * t522;
t479 = t302 + t504;
t303 = pkin(7) * t521;
t478 = t303 + t503;
t474 = t199 * t323 + t200 * t497 + qJD(1);
t469 = -t494 / 0.2e1;
t468 = t494 / 0.2e1;
t289 = rSges(4,1) * t324 + rSges(4,2) * t325;
t372 = -t289 - t552;
t466 = t337 * t565;
t465 = -t298 * t227 - t228 * t520;
t401 = t331 * t274;
t184 = t333 * t401;
t185 = t334 * t401;
t142 = -qJD(5) * t256 + t331 * t485;
t143 = qJD(5) * t255 - t331 * t484;
t87 = rSges(6,1) * t143 + rSges(6,2) * t142 + rSges(6,3) * t481;
t144 = -qJD(5) * t258 + t331 * t483;
t145 = qJD(5) * t257 - t331 * t482;
t88 = rSges(6,1) * t145 + rSges(6,2) * t144 + rSges(6,3) * t480;
t464 = (-t185 + t88) * t334 + (-t184 + t87) * t333;
t132 = -pkin(6) * t334 + t333 * t502;
t133 = pkin(6) * t333 + t334 * t502;
t463 = t333 * t132 + t334 * t133 + t509;
t459 = t288 + t487;
t197 = t459 * t333;
t198 = t459 * t334;
t462 = t333 * t197 + t334 * t198 + t505;
t231 = t567 * t333;
t233 = t567 * t334;
t461 = (t111 + t233) * t334 + (t110 + t231) * t333;
t271 = t566 * qJD(2);
t460 = -t271 - t486;
t307 = rSges(3,1) * t339 - rSges(3,2) * t337;
t385 = (-qJDD(2) * t337 - t339 * t340) * pkin(2);
t348 = (-qJDD(2) * t324 - t325 * t340) * pkin(3) + t385;
t345 = t333 * t348 - t488;
t28 = t111 * t240 - t146 * t178 - t218 * t298 - t229 * t99 - t274 * t296 - t494 * t88 + t345;
t319 = qJDD(3) * t333;
t347 = t334 * t348 + t319;
t29 = -t110 * t240 + t147 * t178 - t218 * t520 + t230 * t99 + t274 * t297 + t494 * t87 + t347;
t453 = -t28 * t334 + t29 * t333;
t452 = t333 * t44 + t541;
t451 = t334 * t47 + t542;
t52 = -t104 * t316 + t315 * t421;
t53 = -t105 * t316 + t315 * t420;
t450 = t52 * t333 + t53 * t334;
t449 = -t333 * t54 - t334 * t55;
t448 = qJD(2) * t372;
t447 = t132 * t323 + t133 * t497 + t474;
t439 = -Icges(6,1) * t336 - t530;
t426 = Icges(5,5) * t316 - Icges(5,6) * t315;
t425 = -Icges(5,5) * t315 - Icges(5,6) * t316;
t423 = -Icges(6,5) * t336 - Icges(6,6) * t338;
t422 = t104 * t230 + t105 * t229;
t419 = t110 * t334 - t111 * t333;
t112 = -t272 * t298 + t360;
t113 = -t272 * t520 + t367;
t418 = -t112 * t333 - t113 * t334;
t417 = -t168 * t336 + t170 * t338;
t416 = -t190 * t315 + t192 * t316;
t415 = -t191 * t315 + t193 * t316;
t209 = -rSges(4,3) * t334 + t333 * t566;
t210 = rSges(4,3) * t333 + t334 * t566;
t412 = t209 * t333 + t210 * t334;
t221 = t425 * t333;
t222 = t425 * t334;
t411 = -t221 * t520 + t222 * t298;
t408 = t565 * t307;
t407 = qJD(2) * t597;
t406 = t199 * t320 + t200 * t489 + t286 * t323 + t287 * t497 + qJDD(1);
t405 = -t272 + t295;
t81 = Icges(6,5) * t143 + Icges(6,6) * t142 + Icges(6,3) * t481;
t404 = t104 * t523 + t315 * t81;
t82 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t480;
t403 = t105 * t523 + t315 * t82;
t386 = t424 * t316;
t94 = t331 * t386 + (Icges(6,3) * t331 + qJD(5) * t423) * t315;
t402 = t166 * t523 + t315 * t94;
t398 = Icges(6,3) * t315 + t386 - t417;
t66 = qJD(2) * t412 + t474;
t397 = t66 * t289;
t395 = qJD(2) * t289;
t389 = t295 + t510;
t388 = t440 * t316;
t387 = t432 * t316;
t384 = (-t112 * t298 - t113 * t520) * t568;
t373 = -qJD(2) * t313 - t486;
t371 = t179 + t567;
t211 = t568 * t331;
t370 = -t211 + t373;
t368 = t373 + t540;
t366 = -(Icges(6,5) * t255 - Icges(6,6) * t256) * t230 - (Icges(6,5) * t257 - Icges(6,6) * t258) * t229 + t423 * t315 * t494;
t365 = t132 * t320 + t133 * t489 + t197 * t323 + t198 * t497 + t406;
t361 = t315 * t366;
t359 = (t315 * t586 + t316 * t588) * t331;
t358 = (t315 * t585 + t316 * t587) * t331;
t353 = -t111 * t476 + t229 * t140 - t141 * t230 + (-pkin(4) * t524 + t303) * t520 + t334 * t477 + t298 * (-pkin(4) * t525 + t302);
t352 = -qJD(2) * t271 - qJDD(2) * t289 + t385;
t350 = (Icges(6,1) * t257 - t107 - t532) * t229 + (Icges(6,1) * t255 - t106 - t533) * t230 - (t439 * t315 - t168) * t494;
t136 = t168 * t333;
t137 = t168 * t334;
t138 = t170 * t333;
t139 = t170 * t334;
t83 = Icges(6,4) * t143 + Icges(6,2) * t142 + Icges(6,6) * t481;
t85 = Icges(6,1) * t143 + Icges(6,4) * t142 + Icges(6,5) * t481;
t14 = (t331 * t421 - t81) * t316 + (t104 * t331 - t336 * t83 + t338 * t85 + (-t106 * t338 - t108 * t336) * qJD(5)) * t315;
t84 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t480;
t86 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t480;
t15 = (t331 * t420 - t82) * t316 + (t105 * t331 - t336 * t84 + t338 * t86 + (-t107 * t338 - t109 * t336) * qJD(5)) * t315;
t169 = Icges(6,6) * t315 + t387;
t171 = Icges(6,5) * t315 + t388;
t20 = t106 * t142 + t108 * t143 + t255 * t83 + t256 * t85 + t333 * t404;
t21 = t107 * t142 + t109 * t143 + t255 * t84 + t256 * t86 + t333 * t403;
t22 = t106 * t144 + t108 * t145 + t257 * t83 + t258 * t85 + t334 * t404;
t23 = t107 * t144 + t109 * t145 + t257 * t84 + t258 * t86 + t334 * t403;
t67 = t315 * t417 - t526;
t25 = t229 * t53 + t230 * t52 - t494 * t67;
t95 = t331 * t387 + (Icges(6,6) * t331 + qJD(5) * t431) * t315;
t96 = t331 * t388 + (Icges(6,5) * t331 + qJD(5) * t439) * t315;
t31 = t142 * t168 + t143 * t170 + t255 * t95 + t256 * t96 + t333 * t402;
t3 = t146 * t45 + t147 * t44 + t20 * t230 + t21 * t229 + t240 * t60 - t31 * t494;
t341 = (-t398 * t494 - t563) * t315;
t344 = (t298 * t587 - t520 * t588) * t316 + (t298 * t585 - t520 * t586) * t315;
t32 = t144 * t168 + t145 * t170 + t257 * t95 + t258 * t96 + t334 * t402;
t4 = t146 * t47 + t147 * t46 + t22 * t230 + t229 * t23 + t240 * t61 - t32 * t494;
t172 = t331 * t221;
t48 = -t172 * t334 + t333 * t359;
t173 = t331 * t222;
t49 = -t173 * t334 + t333 * t358;
t50 = t172 * t333 + t334 * t359;
t51 = t173 * t333 + t334 * t358;
t188 = -Icges(5,3) * t334 + t333 * t426;
t70 = -t188 * t334 + t333 * t416;
t189 = Icges(5,3) * t333 + t334 * t426;
t71 = -t189 * t334 + t333 * t415;
t72 = t188 * t333 + t334 * t416;
t73 = t189 * t333 + t334 * t415;
t346 = (-t20 * t334 + t21 * t333) * t556 - t25 * t495 / 0.2e1 + (t333 * t47 - t334 * t46) * t561 + (t333 * t45 - t334 * t44) * t560 + t298 * (t333 * t51 - t334 * t50) / 0.2e1 - t520 * (t333 * t49 - t334 * t48) / 0.2e1 + (t333 * t53 - t334 * t52) * t555 - t298 * (t333 * t411 + t334 * t344) / 0.2e1 + t520 * (t333 * t344 - t334 * t411) / 0.2e1 + t296 * (t333 * t73 - t334 * t72) / 0.2e1 + t297 * (t333 * t71 - t334 * t70) / 0.2e1 + ((-t137 * t257 - t139 * t258) * t229 + (-t136 * t257 - t138 * t258) * t230 + (t61 * t315 + (-t169 * t257 - t171 * t258 + t542) * t316) * qJD(5) + (((t47 - t526) * qJD(5) + t422) * t316 + t341) * t334) * t559 + ((-t137 * t255 - t139 * t256) * t229 + (-t136 * t255 - t138 * t256) * t230 + (t60 * t315 + (-t169 * t255 - t171 * t256 + t541) * t316) * qJD(5) + (((t44 - t526) * qJD(5) + t422) * t316 + t341) * t333) * t557 + (((t137 * t336 - t139 * t338 + t105) * t229 + (t136 * t336 - t138 * t338 + t104) * t230 + t67 * qJD(5)) * t315 + ((t398 * t316 + (t169 * t336 - t171 * t338 - t166) * t315 + t450) * qJD(5) + t563) * t316) * t468 + (-t22 * t334 + t23 * t333) * t558 + (t296 * t73 + t297 * t72 + t298 * t51 - t50 * t520 + t4) * t554 + (t296 * t71 + t297 * t70 + t298 * t49 - t48 * t520 + t3) * t553 + (-t14 * t334 + t15 * t333 + t579) * t469;
t283 = t306 * t334;
t282 = t306 * t333;
t237 = t454 * t315;
t220 = t334 * t395;
t219 = t333 * t395;
t151 = t334 * t448 + t322;
t150 = t333 * t448 - t496;
t130 = rSges(6,1) * t257 - rSges(6,2) * t258;
t129 = rSges(6,1) * t255 - rSges(6,2) * t256;
t98 = t334 * t352 + t319;
t97 = t333 * t352 - t488;
t74 = -qJD(2) * t407 + qJDD(2) * t408 + qJDD(1);
t69 = -t211 * t520 + t272 * t297 + t347;
t68 = -t211 * t298 - t272 * t296 + t345;
t43 = t412 * qJDD(2) + (-t219 * t333 - t220 * t334) * qJD(2) + t406;
t42 = t195 * t298 + t196 * t520 + t447;
t37 = t110 * t229 - t111 * t230 + t231 * t298 + t233 * t520 + t447;
t33 = -t180 * t298 - t181 * t520 + t195 * t296 - t196 * t297 + t365;
t30 = (t331 * t417 - t94) * t316 + (t166 * t331 - t336 * t95 + t338 * t96 + (-t168 * t338 - t170 * t336) * qJD(5)) * t315;
t11 = t110 * t146 - t111 * t147 - t184 * t298 - t185 * t520 + t229 * t87 - t230 * t88 + t231 * t296 - t233 * t297 + t365;
t1 = [m(2) * qJDD(1) + (-m(2) - m(3) + t490) * g(3) + m(3) * t74 + m(4) * t43 + m(5) * t33 + m(6) * t11; t346 - (t583 * qJD(2) * t329 - t333 * t584 * t497) * t323 / 0.2e1 + (t584 * t330 * qJD(2) - t583 * t334 * t323) * t497 / 0.2e1 + (-g(1) * (t285 + t478) - g(2) * (t284 + t479) - g(3) * (t371 + t501) - t37 * (t353 + t508) - (t449 * t313 + (t339 * t449 - t37 * t466) * pkin(2)) * qJD(2) + t11 * (t461 + t463) + t37 * (t462 + t464) + (t29 * t389 + t368 * t55) * t334 + (t28 * t389 + t368 * t54) * t333 + t582) * m(6) + (-g(1) * (t285 - t228) - g(2) * (t284 - t227) - g(3) * (t568 + t501) - t42 * (t465 + t508) - t384 - (t418 * t313 + (t339 * t418 - t42 * t466) * pkin(2)) * qJD(2) + t33 * (t463 + t512) + t42 * (t462 + t513) + (t113 * t370 + t405 * t69) * t334 + (t112 * t370 + t405 * t68) * t333) * m(5) + (-g(3) * t569 - t564 * t372 - (-t66 * pkin(2) * t466 + (-t151 * t569 - t334 * t397) * t334 + (-t150 * t569 - t333 * t397) * t333) * qJD(2) + t43 * t509 + t66 * t505 + (t151 * t460 + t43 * t210 - t66 * t220 + t372 * t98) * t334 + (t150 * t460 + t43 * t209 - t66 * t219 + t372 * t97) * t333) * m(4) + (g(1) * t283 + g(2) * t282 - g(3) * t307 + t74 * t408 + (qJDD(2) * t306 + t307 * t340) * t597 + (-(-t282 * t333 - t283 * t334) * qJD(2) - t407) * (qJD(2) * t408 + qJD(1))) * m(3) + 0.2e1 * ((t589 * t330 + (t595 * t333 + (t590 - t596) * t334) * t333) * qJDD(2) + (t591 * t330 + (t593 * t333 + (t592 - t594) * t334) * t333) * qJD(2)) * t554 + 0.2e1 * ((t596 * t330 + (t590 * t333 + (t589 - t595) * t334) * t333) * qJDD(2) + (t594 * t330 + (t592 * t333 + (t591 - t593) * t334) * t333) * qJD(2)) * t553; t490 * (g(1) * t333 - g(2) * t334) + m(4) * (t333 * t98 - t334 * t97) + m(5) * (t333 * t69 - t334 * t68) + m(6) * t453; t346 + (t11 * t461 + (t29 * t510 + t540 * t55) * t334 + (t28 * t510 + t54 * t540) * t333 - g(1) * t478 - g(2) * t479 - g(3) * t371 + (-t353 + t464) * t37 + t582) * m(6) + (-t384 + t33 * t512 + (-t333 * t68 - t334 * t69) * t272 + t418 * t211 + g(1) * t228 + g(2) * t227 - g(3) * t568 + (-t465 + t513) * t42) * m(5); t4 * t524 / 0.2e1 + (t315 * t451 - t316 * t61) * t561 + ((t331 * t451 - t32) * t316 + (t22 * t333 + t23 * t334 + t331 * t61) * t315) * t558 + t3 * t525 / 0.2e1 + (t315 * t452 - t316 * t60) * t560 + ((t331 * t452 - t31) * t316 + (t20 * t333 + t21 * t334 + t331 * t60) * t315) * t556 + t331 * t315 * t25 / 0.2e1 - t316 * (t14 * t230 + t146 * t53 + t147 * t52 + t15 * t229 + t240 * t67 - t30 * t494) / 0.2e1 + (t315 * t450 - t316 * t67) * t555 + ((t331 * t450 - t30) * t316 + (t14 * t333 + t15 * t334 + t331 * t67) * t315) * t469 + (t257 * t349 + t258 * t350 - t334 * t361) * t559 + (t255 * t349 + t256 * t350 - t333 * t361) * t557 + (t366 * t316 + (-t336 * t349 + t350 * t338) * t315) * t468 + t579 * t523 / 0.2e1 + ((t29 * t110 - t28 * t111 - t54 * t88 + t55 * t87 + (t37 * t419 + (t333 * t55 - t334 * t54) * t178) * t331) * t316 + (t55 * (-t110 * t331 + t333 * t99) + t54 * (t111 * t331 - t334 * t99) + t11 * t419 + t37 * (-t333 * t88 + t334 * t87) + t453 * t178) * t315 - t55 * (t129 * t494 + t230 * t237) - t54 * (-t130 * t494 - t229 * t237) - t37 * (t129 * t229 - t130 * t230) - g(1) * t130 - g(2) * t129 - g(3) * t237) * m(6);];
tau = t1;
