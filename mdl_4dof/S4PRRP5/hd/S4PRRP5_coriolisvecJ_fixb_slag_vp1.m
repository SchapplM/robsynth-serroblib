% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP5_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:51
% EndTime: 2019-12-31 16:29:17
% DurationCPUTime: 23.65s
% Computational Cost: add. (7529->550), mult. (20333->817), div. (0->0), fcn. (19657->6), ass. (0->285)
t534 = Icges(4,5) + Icges(5,5);
t505 = Icges(4,6) + Icges(5,6);
t232 = sin(pkin(6));
t233 = cos(pkin(6));
t237 = cos(qJ(3));
t235 = sin(qJ(3));
t238 = cos(qJ(2));
t380 = t235 * t238;
t198 = -t232 * t380 - t233 * t237;
t379 = t237 * t238;
t383 = t233 * t235;
t199 = t232 * t379 - t383;
t236 = sin(qJ(2));
t385 = t232 * t236;
t83 = Icges(5,5) * t199 + Icges(5,6) * t198 + Icges(5,3) * t385;
t85 = Icges(4,5) * t199 + Icges(4,6) * t198 + Icges(4,3) * t385;
t496 = t83 + t85;
t200 = t232 * t237 - t233 * t380;
t386 = t232 * t235;
t201 = t233 * t379 + t386;
t382 = t233 * t236;
t84 = Icges(5,5) * t201 + Icges(5,6) * t200 + Icges(5,3) * t382;
t86 = Icges(4,5) * t201 + Icges(4,6) * t200 + Icges(4,3) * t382;
t495 = t84 + t86;
t398 = Icges(5,4) * t199;
t87 = Icges(5,2) * t198 + Icges(5,6) * t385 + t398;
t402 = Icges(4,4) * t199;
t89 = Icges(4,2) * t198 + Icges(4,6) * t385 + t402;
t520 = t87 + t89;
t397 = Icges(5,4) * t201;
t88 = Icges(5,2) * t200 + Icges(5,6) * t382 + t397;
t401 = Icges(4,4) * t201;
t90 = Icges(4,2) * t200 + Icges(4,6) * t382 + t401;
t519 = t88 + t90;
t192 = Icges(5,4) * t198;
t91 = Icges(5,1) * t199 + Icges(5,5) * t385 + t192;
t194 = Icges(4,4) * t198;
t93 = Icges(4,1) * t199 + Icges(4,5) * t385 + t194;
t518 = t91 + t93;
t193 = Icges(5,4) * t200;
t92 = Icges(5,1) * t201 + Icges(5,5) * t382 + t193;
t195 = Icges(4,4) * t200;
t94 = Icges(4,1) * t201 + Icges(4,5) * t382 + t195;
t517 = t92 + t94;
t482 = Icges(4,1) + Icges(5,1);
t533 = Icges(4,2) + Icges(5,2);
t541 = Icges(4,3) + Icges(5,3);
t297 = Icges(5,5) * t237 - Icges(5,6) * t235;
t164 = -Icges(5,3) * t238 + t236 * t297;
t298 = Icges(4,5) * t237 - Icges(4,6) * t235;
t166 = -Icges(4,3) * t238 + t236 * t298;
t531 = -t166 - t164;
t395 = Icges(5,4) * t237;
t300 = -Icges(5,2) * t235 + t395;
t168 = -Icges(5,6) * t238 + t236 * t300;
t399 = Icges(4,4) * t237;
t301 = -Icges(4,2) * t235 + t399;
t170 = -Icges(4,6) * t238 + t236 * t301;
t481 = t168 + t170;
t396 = Icges(5,4) * t235;
t304 = Icges(5,1) * t237 - t396;
t172 = -Icges(5,5) * t238 + t236 * t304;
t400 = Icges(4,4) * t235;
t305 = Icges(4,1) * t237 - t400;
t174 = -Icges(4,5) * t238 + t236 * t305;
t488 = t172 + t174;
t540 = (-t534 * t235 - t505 * t237) * t236;
t502 = t519 * t198 + t517 * t199 + t495 * t385;
t501 = t520 * t200 + t518 * t201 + t496 * t382;
t503 = t520 * t198 + t518 * t199 + t385 * t496;
t465 = t200 * t519 + t201 * t517 + t382 * t495;
t539 = t198 * t481 + t199 * t488 - t385 * t531;
t538 = t200 * t481 + t201 * t488 - t382 * t531;
t359 = qJD(2) * t236;
t343 = t235 * t359;
t128 = -qJD(3) * t199 + t232 * t343;
t342 = t237 * t359;
t129 = qJD(3) * t198 - t232 * t342;
t358 = qJD(2) * t238;
t345 = t232 * t358;
t537 = t128 * t505 + t129 * t534 + t345 * t541;
t130 = -qJD(3) * t201 + t233 * t343;
t131 = qJD(3) * t200 - t233 * t342;
t344 = t233 * t358;
t536 = t130 * t505 + t131 * t534 + t344 * t541;
t535 = Icges(4,4) + Icges(5,4);
t165 = Icges(5,3) * t236 + t238 * t297;
t167 = Icges(4,3) * t236 + t238 * t298;
t532 = t540 * qJD(3) + (t165 + t167) * qJD(2);
t169 = Icges(5,6) * t236 + t238 * t300;
t171 = Icges(4,6) * t236 + t238 * t301;
t530 = t169 + t171;
t173 = Icges(5,5) * t236 + t238 * t304;
t175 = Icges(4,5) * t236 + t238 * t305;
t457 = -t173 - t175;
t529 = (-t237 * t533 - t396 - t400) * t236;
t528 = (t482 * t235 + t395 + t399) * t236;
t527 = t502 * t233;
t526 = t501 * t232;
t230 = t232 ^ 2;
t231 = t233 ^ 2;
t444 = t230 + t231;
t525 = m(3) * qJD(2) ^ 2 * (rSges(3,1) * t236 + rSges(3,2) * t238) * t444;
t524 = t533 * t128 + t535 * t129 + t505 * t345;
t523 = t533 * t130 + t535 * t131 + t505 * t344;
t522 = t535 * t128 + t482 * t129 + t534 * t345;
t521 = t535 * t130 + t482 * t131 + t534 * t344;
t516 = t530 * qJD(2) + t529 * qJD(3);
t515 = t457 * qJD(2) + t528 * qJD(3);
t514 = -t532 * t236 + t531 * t358;
t513 = t536 * t236 + t495 * t358;
t512 = t537 * t236 + t496 * t358;
t511 = t465 * t233 + t526;
t510 = t503 * t232 + t527;
t509 = t538 * t236;
t508 = t539 * t236;
t357 = qJD(3) * t236;
t340 = t233 * t357;
t361 = qJD(2) * t232;
t218 = t340 + t361;
t341 = t232 * t357;
t360 = qJD(2) * t233;
t219 = t341 - t360;
t504 = t218 * t232 - t219 * t233;
t431 = t218 / 0.2e1;
t429 = t219 / 0.2e1;
t471 = t128 * t520 + t129 * t518 + t198 * t524 + t199 * t522 + t232 * t512;
t470 = t128 * t519 + t129 * t517 + t198 * t523 + t199 * t521 + t232 * t513;
t469 = t130 * t520 + t131 * t518 + t200 * t524 + t201 * t522 + t233 * t512;
t468 = t130 * t519 + t131 * t517 + t200 * t523 + t201 * t521 + t233 * t513;
t312 = -t235 * t87 + t237 * t91;
t47 = t236 * t312 - t238 * t83;
t310 = -t235 * t89 + t237 * t93;
t49 = t236 * t310 - t238 * t85;
t500 = t47 + t49;
t311 = -t235 * t88 + t237 * t92;
t48 = t236 * t311 - t238 * t84;
t309 = -t235 * t90 + t237 * t94;
t50 = t236 * t309 - t238 * t86;
t499 = t48 + t50;
t295 = -t168 * t235 + t172 * t237;
t388 = t164 * t238;
t61 = t236 * t295 - t388;
t294 = -t170 * t235 + t174 * t237;
t387 = t166 * t238;
t62 = t236 * t294 - t387;
t484 = t61 + t62;
t487 = (-t130 * t481 - t131 * t488 - t200 * t516 + t201 * t515 + t233 * t514) * t238 + (t238 * t511 + t509) * qJD(2);
t486 = (-t128 * t481 - t129 * t488 - t516 * t198 + t199 * t515 + t232 * t514) * t238 + (t238 * t510 + t508) * qJD(2);
t356 = qJD(3) * t238;
t441 = t218 * (-t201 * t533 + t193 + t195 + t517) + t219 * (-t199 * t533 + t192 + t194 + t518) - t356 * (t488 + t529);
t485 = t469 * t429 + t468 * t431 + t487 * qJD(3) / 0.2e1;
t480 = t540 * t356 + (-t198 * t534 + t505 * t199) * t219 + (-t200 * t534 + t505 * t201) * t218;
t313 = t49 * t232 + t50 * t233;
t314 = t47 * t232 + t48 * t233;
t479 = t313 + t314;
t475 = (t218 * t502 + t219 * t503 - t356 * t539) * t232 + (t218 * t465 + t219 * t501 - t356 * t538) * t233;
t474 = qJD(3) * t486 + t218 * t470 + t219 * t471;
t448 = ((t310 + t312) * qJD(2) - t537) * t238 + (t522 * t237 - t524 * t235 + (-t235 * t518 - t237 * t520) * qJD(3) + t496 * qJD(2)) * t236;
t447 = ((t309 + t311) * qJD(2) - t536) * t238 + (t521 * t237 - t523 * t235 + (-t235 * t517 - t237 * t519) * qJD(3) + t495 * qJD(2)) * t236;
t472 = rSges(5,1) + pkin(3);
t467 = t218 * t499 + t219 * t500 - t356 * t484;
t136 = t168 * t232;
t138 = t170 * t232;
t463 = -t136 - t138;
t137 = t168 * t233;
t139 = t170 * t233;
t462 = -t137 - t139;
t140 = t172 * t232;
t142 = t174 * t232;
t461 = -t140 - t142;
t141 = t172 * t233;
t143 = t174 * t233;
t460 = -t141 - t143;
t425 = pkin(3) * t237;
t149 = qJ(4) * t236 + t238 * t425;
t321 = rSges(5,1) * t237 - rSges(5,2) * t235;
t459 = rSges(5,3) * t236 + t238 * t321 + t149;
t280 = t167 - t294;
t281 = t165 - t295;
t437 = -(-t164 * t233 - t311) * t218 - (-t164 * t232 - t312) * t219;
t438 = -(-t166 * t233 - t309) * t218 - (-t166 * t232 - t310) * t219;
t456 = (-t437 - t438 + (-t280 - t281) * t356) * t236;
t269 = qJ(4) * t238 - t236 * t425;
t455 = rSges(5,3) * t238 - t236 * t321 + t269;
t454 = t218 * t495 + t219 * t496;
t224 = pkin(2) * t236 - pkin(5) * t238;
t279 = qJD(2) * t224;
t190 = t232 * t279;
t191 = t233 * t279;
t453 = -t232 * t190 - t233 * t191 + t444 * t279;
t452 = -t387 - t388;
t403 = Icges(3,4) * t238;
t404 = Icges(3,4) * t236;
t451 = -(-Icges(3,2) * t238 - t404) * t359 + (-Icges(3,1) * t236 - t403) * t358;
t225 = pkin(2) * t238 + pkin(5) * t236;
t221 = qJD(2) * t225;
t354 = qJD(4) * t238;
t330 = -t221 + t354;
t443 = t484 * t359 + (t532 * t238 + (t515 * t237 + t516 * t235 + (t235 * t488 + t481 * t237) * qJD(3)) * t236 + ((-t294 - t295) * t238 + t531 * t236 + t479) * qJD(2)) * t238;
t442 = (t481 + t528) * t356 + (t198 * t482 - t398 - t402 - t520) * t219 + (t200 * t482 - t397 - t401 - t519) * t218;
t440 = t480 * t236;
t322 = rSges(4,1) * t237 - rSges(4,2) * t235;
t177 = -rSges(4,3) * t238 + t236 * t322;
t208 = t225 * t232;
t209 = t225 * t233;
t337 = t208 * t361 + t209 * t360 + qJD(1);
t420 = rSges(5,1) * t201 + rSges(5,2) * t200 + rSges(5,3) * t382 + pkin(3) * t386 + t149 * t233;
t421 = rSges(5,1) * t199 + rSges(5,2) * t198 + rSges(5,3) * t385 - pkin(3) * t383 + t149 * t232;
t31 = t218 * t421 - t219 * t420 + t337 - t354;
t355 = qJD(4) * t236;
t227 = t233 * t355;
t417 = pkin(3) * qJD(3);
t352 = t235 * t417;
t243 = qJD(2) * t269 - t238 * t352;
t351 = t237 * t417;
t422 = rSges(5,1) * t131 + rSges(5,2) * t130 + rSges(5,3) * t344 + t232 * t351 + t233 * t243 + t227;
t370 = -t190 * t361 - t191 * t360;
t226 = t232 * t355;
t423 = rSges(5,1) * t129 + rSges(5,2) * t128 + rSges(5,3) * t345 + t232 * t243 - t233 * t351 + t226;
t5 = -t422 * t219 + t423 * t218 + (t355 + (-t232 * t420 + t233 * t421) * t356) * qJD(2) + t370;
t439 = t31 * t422 + t420 * t5;
t434 = -m(5) / 0.2e1;
t433 = m(5) / 0.2e1;
t432 = -t218 / 0.2e1;
t430 = -t219 / 0.2e1;
t384 = t232 * t238;
t381 = t233 * t238;
t216 = (-rSges(5,1) * t235 - rSges(5,2) * t237) * t236;
t378 = qJD(2) * t459 + qJD(3) * t216 - t236 * t352 - t354;
t179 = rSges(4,3) * t236 + t238 * t322;
t217 = (-rSges(4,1) * t235 - rSges(4,2) * t237) * t236;
t109 = qJD(2) * t179 + qJD(3) * t217;
t377 = -t109 - t221;
t376 = t455 * t232;
t375 = t455 * t233;
t374 = -rSges(5,2) * t199 + t198 * t472;
t373 = rSges(5,2) * t201 - t200 * t472;
t363 = -t177 - t224;
t362 = t232 * t208 + t233 * t209;
t353 = qJD(2) * qJD(3);
t349 = -t221 - t378;
t348 = -t224 + t455;
t347 = t224 * t361;
t346 = t224 * t360;
t333 = -t356 / 0.2e1;
t332 = t356 / 0.2e1;
t331 = t353 / 0.2e1;
t325 = pkin(3) * t235 * t236 - t216;
t324 = t236 * t331;
t323 = t238 * t331;
t100 = rSges(4,1) * t201 + rSges(4,2) * t200 + rSges(4,3) * t382;
t98 = rSges(4,1) * t199 + rSges(4,2) * t198 + rSges(4,3) * t385;
t308 = -t100 * t232 + t233 * t98;
t307 = Icges(3,1) * t238 - t404;
t303 = -Icges(3,2) * t236 + t403;
t299 = -Icges(3,5) * t236 - Icges(3,6) * t238;
t292 = t232 * t323;
t291 = t233 * t323;
t46 = -t100 * t219 + t218 * t98 + t337;
t276 = t46 * t308;
t272 = qJD(2) * t299;
t271 = t31 * t423 + t421 * t5;
t36 = t218 * t455 - t356 * t420 + t226 - t347;
t37 = -t219 * t455 + t356 * t421 + t227 - t346;
t270 = t36 * t420 - t37 * t421;
t255 = (-(-Icges(3,6) * t233 + t232 * t303) * t238 - (-Icges(3,5) * t233 + t232 * t307) * t236) * qJD(2) + t451 * t232;
t254 = (-(Icges(3,6) * t232 + t233 * t303) * t238 - (Icges(3,5) * t232 + t233 * t307) * t236) * qJD(2) + t451 * t233;
t242 = (t31 * t421 + t36 * t455) * t233 + (-t31 * t420 - t37 * t455) * t232;
t203 = t299 * t233;
t202 = t299 * t232;
t183 = t233 * t272;
t182 = t232 * t272;
t127 = rSges(4,1) * t200 - rSges(4,2) * t201;
t125 = rSges(4,1) * t198 - rSges(4,2) * t199;
t78 = rSges(4,1) * t131 + rSges(4,2) * t130 + rSges(4,3) * t344;
t76 = rSges(4,1) * t129 + rSges(4,2) * t128 + rSges(4,3) * t345;
t60 = t177 * t219 + t356 * t98 - t346;
t59 = -t100 * t356 - t177 * t218 - t347;
t35 = t76 * t356 + t109 * t219 + (-t221 * t233 + (t177 * t384 - t236 * t98) * qJD(3)) * qJD(2);
t34 = -t78 * t356 - t109 * t218 + (-t221 * t232 + (t100 * t236 - t177 * t381) * qJD(3)) * qJD(2);
t26 = t238 * t308 * t353 + t218 * t76 - t219 * t78 + t370;
t23 = t378 * t219 + t423 * t356 + (t330 * t233 + (-t236 * t421 - t384 * t455) * qJD(3)) * qJD(2);
t22 = -t378 * t218 - t422 * t356 + (t330 * t232 + (t236 * t420 + t381 * t455) * qJD(3)) * qJD(2);
t1 = [m(4) * t26 + m(5) * t5 - t525; (t202 * qJD(2) * t231 - t203 * t233 * t361) * t360 / 0.2e1 - (t232 * (-t183 * t233 + t232 * t254) - t233 * (-t182 * t233 + t232 * t255)) * t360 - (t203 * qJD(2) * t230 - t202 * t232 * t360) * t361 / 0.2e1 + (t232 * (t183 * t232 + t233 * t254) - t233 * (t182 * t232 + t233 * t255)) * t361 + (((-t200 * t530 + t201 * t457 + t526) * t238 + t509) * qJD(3) + (((t452 + t465) * qJD(3) + t454) * t238 + t456) * t233 + (t200 * t463 + t201 * t461) * t219 + (t200 * t462 + t201 * t460) * t218) * t432 + (t232 * t468 - t233 * t469) * t431 + (((-t198 * t530 + t199 * t457 + t527) * t238 + t508) * qJD(3) + (((t452 + t503) * qJD(3) + t454) * t238 + t456) * t232 + (t198 * t463 + t199 * t461) * t219 + (t198 * t462 + t199 * t460) * t218) * t430 + (t232 * t470 - t233 * t471) * t429 + t232 * t485 - t474 * t233 / 0.2e1 - t467 * t357 / 0.2e1 + (((t137 * t235 - t141 * t237 + t84) * t218 + (t136 * t235 - t140 * t237 + t83) * t219 + t61 * qJD(3)) * t236 + ((t281 * t238 + (t169 * t235 - t173 * t237 - t164) * t236 + t314) * qJD(3) + t437) * t238 + ((t139 * t235 - t143 * t237 + t86) * t218 + (t138 * t235 - t142 * t237 + t85) * t219 + t62 * qJD(3)) * t236 + ((t280 * t238 + (t171 * t235 - t175 * t237 - t166) * t236 + t313) * qJD(3) + t438) * t238) * t332 + (t232 * t499 - t233 * t500) * t324 + (t232 * t502 - t233 * t503) * t292 + (t465 * t232 - t233 * t501) * t291 + (-(t270 * t236 + (-t36 * t375 + t37 * t376 + t242) * t238) * qJD(3) - (t232 * t36 + t233 * t37) * t330 + t5 * t362 + (t23 * t348 + t349 * t37 + t439) * t233 + (t22 * t348 + t349 * t36 + t271) * t232 + (t218 * t36 - t219 * t37) * t459 + (-t218 * t376 + t219 * t375 - t355 + t453) * t31) * m(5) + (-t60 * (t179 * t219 - t225 * t360) - t59 * (-t179 * t218 - t225 * t361) - ((t100 * t59 - t60 * t98) * t236 + t276 * t238) * qJD(3) + t26 * t362 + (t26 * t100 + t35 * t363 + t377 * t60) * t233 + (t26 * t98 + t34 * t363 + t377 * t59) * t232 + (t177 * t504 + t232 * t76 + t233 * t78 + t453) * t46) * m(4) + (t232 * t447 - t233 * t448 + t475) * t333 + (-t444 + 0.1e1) * (rSges(3,1) * t238 - rSges(3,2) * t236) * t525; (t200 * t441 + t201 * t442 - t233 * t440) * t432 + ((t232 * t469 + t233 * t468) * t236 + t487) * t431 + (t198 * t441 + t199 * t442 - t232 * t440) * t430 + ((t232 * t471 + t233 * t470) * t236 + t486) * t429 - (qJD(3) * t443 + t218 * t447 + t219 * t448) * t238 / 0.2e1 + t474 * t385 / 0.2e1 + t382 * t485 + t467 * t359 / 0.2e1 + ((t232 * t448 + t233 * t447) * t236 + t443) * t333 + (t480 * t238 + (-t235 * t441 + t442 * t237) * t236) * t332 + (t236 * t479 - t238 * t484) * t324 + (t236 * t510 - t238 * t539) * t292 + (t236 * t511 - t238 * t538) * t291 + (-(t31 * t373 - t325 * t37) * t219 - (t31 * t374 + t325 * t36) * t218 - (t36 * t373 + t37 * t374) * t356 + (qJD(2) * t242 - t22 * t420 + t23 * t421 - t36 * t422 + t37 * t423) * t238 + (t270 * qJD(2) + (t22 * t455 - t36 * t378 + t271) * t233 + (-t23 * t455 + t37 * t378 - t439) * t232) * t236) * m(5) + (-t60 * (t125 * t356 + t217 * t219) - t59 * (-t127 * t356 - t217 * t218) - t46 * (t125 * t218 - t127 * t219) + (-t34 * t100 + t35 * t98 - t59 * t78 + t60 * t76 + (t276 + (t232 * t60 - t233 * t59) * t177) * qJD(2)) * t238 + (t60 * (-qJD(2) * t98 + t109 * t232) + t59 * (qJD(2) * t100 - t109 * t233) + t26 * t308 + t46 * (-t232 * t78 + t233 * t76) + (t232 * t35 - t233 * t34) * t177) * t236) * m(4) + t475 * t358 / 0.2e1; 0.2e1 * ((qJD(2) * t31 + t22 * t232 + t23 * t233) * t433 + t31 * t504 * t434) * t236 + 0.2e1 * ((t36 * t361 + t360 * t37 - t5) * t433 + (t37 * (-t219 + t341) + t36 * (t218 - t340)) * t434) * t238;];
tauc = t1(:);
