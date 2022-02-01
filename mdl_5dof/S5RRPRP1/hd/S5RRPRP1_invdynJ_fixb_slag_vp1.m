% Calculate vector of inverse dynamics joint torques for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:34
% EndTime: 2022-01-20 10:19:55
% DurationCPUTime: 14.50s
% Computational Cost: add. (13377->505), mult. (9752->591), div. (0->0), fcn. (7457->8), ass. (0->295)
t284 = sin(qJ(4));
t286 = cos(qJ(4));
t228 = Icges(6,5) * t286 - Icges(6,6) * t284;
t230 = Icges(5,5) * t286 - Icges(5,6) * t284;
t539 = t228 + t230;
t560 = Icges(5,5) + Icges(6,5);
t559 = Icges(5,6) + Icges(6,6);
t558 = Icges(5,3) + Icges(6,3);
t282 = qJ(1) + qJ(2);
t271 = pkin(8) + t282;
t265 = cos(t271);
t264 = sin(t271);
t424 = t264 * t286;
t425 = t264 * t284;
t136 = Icges(6,4) * t424 - Icges(6,2) * t425 - Icges(6,6) * t265;
t138 = Icges(5,4) * t424 - Icges(5,2) * t425 - Icges(5,6) * t265;
t551 = t136 + t138;
t220 = Icges(6,4) * t425;
t140 = Icges(6,1) * t424 - Icges(6,5) * t265 - t220;
t221 = Icges(5,4) * t425;
t142 = Icges(5,1) * t424 - Icges(5,5) * t265 - t221;
t556 = t140 + t142;
t445 = Icges(6,4) * t284;
t236 = Icges(6,1) * t286 - t445;
t323 = t236 * t265;
t141 = Icges(6,5) * t264 + t323;
t446 = Icges(5,4) * t284;
t238 = Icges(5,1) * t286 - t446;
t324 = t238 * t265;
t143 = Icges(5,5) * t264 + t324;
t547 = t141 + t143;
t557 = t539 * t265;
t531 = t558 * t265 - t424 * t560 + t559 * t425;
t530 = t264 * t558 + t557;
t231 = Icges(6,2) * t286 + t445;
t233 = Icges(5,2) * t286 + t446;
t555 = t231 + t233;
t274 = Icges(6,4) * t286;
t235 = Icges(6,1) * t284 + t274;
t275 = Icges(5,4) * t286;
t237 = Icges(5,1) * t284 + t275;
t554 = -t235 - t237;
t553 = t551 * t284;
t552 = t547 * t424;
t342 = -Icges(6,2) * t284 + t274;
t321 = t342 * t265;
t137 = Icges(6,6) * t264 + t321;
t343 = -Icges(5,2) * t284 + t275;
t322 = t343 * t265;
t139 = Icges(5,6) * t264 + t322;
t550 = t137 + t139;
t549 = t236 + t238;
t527 = t556 * t286 - t553;
t548 = t342 + t343;
t227 = Icges(6,5) * t284 + Icges(6,6) * t286;
t229 = Icges(5,5) * t284 + Icges(5,6) * t286;
t546 = t227 + t229;
t281 = qJD(1) + qJD(2);
t545 = -t555 * qJD(4) + t281 * t559;
t544 = t554 * qJD(4) + t281 * t560;
t538 = t555 * t284 + t554 * t286;
t543 = -t530 * t265 + t552;
t421 = t265 * t286;
t502 = t530 * t264 + t547 * t421;
t542 = t531 * t264 - t556 * t421;
t489 = t527 * t264 + t531 * t265;
t488 = -t425 * t550 + t543;
t422 = t265 * t284;
t487 = -t422 * t551 - t542;
t486 = -t422 * t550 + t502;
t541 = t548 * qJD(4);
t540 = t549 * qJD(4);
t537 = t554 * t284 - t555 * t286;
t536 = t550 * t284;
t485 = t556 * t284 + t551 * t286;
t484 = t547 * t284 + t286 * t550;
t427 = t264 * t281;
t535 = t545 * t265 - t427 * t548;
t534 = (t321 + t322) * t281 + t545 * t264;
t533 = t544 * t265 - t427 * t549;
t532 = (t323 + t324) * t281 + t544 * t264;
t429 = t229 * t265;
t432 = t227 * t265;
t508 = -t538 * t264 - t429 - t432;
t430 = t229 * t264;
t433 = t227 * t264;
t507 = -t538 * t265 + t430 + t433;
t529 = -t546 * qJD(4) + t281 * t558;
t528 = -t547 * t286 + t536;
t526 = t537 * qJD(4) + t546 * t281 - t541 * t284 + t540 * t286;
t525 = t486 * t264 - t487 * t265;
t524 = t488 * t264 - t489 * t265;
t283 = -qJ(5) - pkin(7);
t240 = t265 * t283;
t277 = t286 * pkin(4);
t267 = t277 + pkin(3);
t523 = -rSges(6,1) * t424 - t264 * t267 - t240;
t522 = t539 * qJD(4) + t538 * t281;
t521 = t507 * t281;
t261 = t265 * pkin(7);
t185 = pkin(3) * t264 - t261;
t480 = -rSges(6,2) * t425 - t265 * rSges(6,3);
t412 = -t185 + t480 - t523;
t520 = -t484 * qJD(4) + t530 * t281 - t535 * t284 + t533 * t286;
t519 = t485 * qJD(4) + t531 * t281 + t534 * t284 - t532 * t286;
t518 = t508 * t281;
t517 = (-t527 + t557) * t281 + t529 * t264;
t516 = t529 * t265 + t528 * t281 - t539 * t427;
t515 = t524 * qJD(4) + t518;
t514 = t525 * qJD(4) + t521;
t285 = sin(qJ(1));
t287 = cos(qJ(1));
t289 = qJD(1) ^ 2;
t318 = (-qJDD(1) * t285 - t287 * t289) * pkin(1);
t513 = t318 - g(1);
t512 = t527 * qJD(4) + t532 * t284 + t534 * t286;
t511 = -t528 * qJD(4) + t533 * t284 + t535 * t286;
t510 = t522 * t264 + t526 * t265;
t509 = t526 * t264 - t522 * t265;
t328 = rSges(6,1) * t421 + t264 * rSges(6,3);
t506 = -t265 * t267 - t328;
t272 = sin(t282);
t457 = pkin(2) * t272;
t505 = -t457 + t523;
t504 = t531 + t536;
t390 = qJD(4) * t286;
t418 = t281 * t284;
t503 = t264 * t390 + t265 * t418;
t451 = pkin(1) * qJD(1);
t385 = t285 * t451;
t273 = cos(t282);
t192 = rSges(3,1) * t272 + rSges(3,2) * t273;
t434 = t192 * t281;
t151 = -t385 - t434;
t501 = -t264 * t517 + t265 * t519;
t500 = t264 * t516 + t265 * t520;
t260 = t264 * pkin(7);
t262 = t265 * pkin(3);
t186 = t262 + t260;
t153 = t186 * t281;
t392 = qJD(4) * t281;
t173 = -qJDD(4) * t265 + t264 * t392;
t244 = t286 * rSges(6,1) - rSges(6,2) * t284;
t203 = t244 * qJD(4);
t450 = t286 * rSges(6,2);
t241 = rSges(6,1) * t284 + t450;
t248 = qJD(5) * t265;
t280 = qJDD(1) + qJDD(2);
t279 = t281 ^ 2;
t419 = t273 * t279;
t303 = -pkin(2) * t419 + t318;
t365 = -t185 - t457;
t332 = t365 - t412;
t393 = qJD(4) * t265;
t417 = t286 * qJD(4) ^ 2;
t391 = qJD(4) * t284;
t377 = t264 * t391;
t400 = pkin(4) * t377 + t248;
t426 = t264 * t283;
t327 = rSges(6,1) * t377 + rSges(6,2) * t503 + t281 * t426 + t400;
t456 = pkin(3) - t267;
t454 = -t327 + (-t265 * t456 - t260 + t328) * t281;
t14 = -t203 * t393 + qJDD(5) * t264 + t173 * t241 + (t173 * t284 - t265 * t417) * pkin(4) + (-t153 + t248 - t454) * t281 + t332 * t280 + t303;
t499 = t14 - g(1);
t172 = qJDD(4) * t264 + t265 * t392;
t423 = t265 * t281;
t215 = pkin(7) * t423;
t266 = pkin(2) * t273;
t278 = t287 * pkin(1);
t458 = pkin(1) * t285;
t356 = qJDD(1) * t278 - t289 * t458;
t312 = t280 * t266 - t279 * t457 + t356;
t302 = t281 * (-pkin(3) * t427 + t215) + t280 * t186 + t312;
t363 = -pkin(4) * t284 - t241;
t483 = -rSges(6,2) * t422 - t426 - t506;
t411 = -t186 + t483;
t375 = t265 * t391;
t316 = -t281 * t424 - t375;
t247 = qJD(5) * t264;
t330 = -pkin(4) * t375 + t247;
t374 = t265 * t390;
t384 = t264 * t418;
t402 = rSges(6,2) * t384 + rSges(6,3) * t423;
t455 = -t215 + (t264 * t456 - t240) * t281 + t330 + rSges(6,1) * t316 - rSges(6,2) * t374 + t402;
t15 = -qJDD(5) * t265 + t455 * t281 + t411 * t280 + t363 * t172 + (-pkin(4) * t417 - qJD(4) * t203 + qJD(5) * t281) * t264 + t302;
t498 = t15 - g(2);
t453 = rSges(5,1) * t286;
t245 = -rSges(5,2) * t284 + t453;
t204 = t245 * qJD(4);
t242 = rSges(5,1) * t284 + rSges(5,2) * t286;
t399 = rSges(5,2) * t425 + t265 * rSges(5,3);
t145 = rSges(5,1) * t424 - t399;
t354 = -t145 + t365;
t380 = rSges(5,1) * t377 + rSges(5,2) * t503;
t225 = rSges(5,1) * t421;
t479 = t264 * rSges(5,3) + t225;
t96 = t281 * t479 - t380;
t26 = -t204 * t393 + t173 * t242 + (-t153 - t96) * t281 + t354 * t280 + t303;
t497 = t26 - g(1);
t147 = -rSges(5,2) * t422 + t479;
t394 = qJD(4) * t264;
t325 = rSges(5,3) * t423 + (-t374 + t384) * rSges(5,2);
t94 = rSges(5,1) * t316 + t325;
t27 = t147 * t280 - t172 * t242 - t204 * t394 + t281 * t94 + t302;
t496 = t27 - g(2);
t183 = rSges(4,1) * t264 + rSges(4,2) * t265;
t212 = rSges(4,2) * t427;
t495 = -t280 * t183 - t281 * (rSges(4,1) * t423 - t212) + (-t272 * t280 - t419) * pkin(2) + t513;
t259 = t265 * rSges(4,1);
t184 = -rSges(4,2) * t264 + t259;
t494 = -t183 * t279 + t280 * t184 - g(2) + t312;
t493 = t264 * t519 + t265 * t517;
t263 = t273 * rSges(3,1);
t420 = t272 * t281;
t177 = -rSges(3,2) * t420 + t263 * t281;
t492 = -t177 * t281 - t192 * t280 + t513;
t491 = t264 * t520 - t265 * t516;
t193 = -rSges(3,2) * t272 + t263;
t490 = t193 * t280 - t281 * t434 - g(2) + t356;
t364 = t186 + t266;
t482 = t147 + t364;
t481 = t184 + t266;
t478 = t244 + t277;
t129 = t281 * t145;
t178 = t281 * t185;
t477 = -rSges(5,1) * t375 + t129 + t178 + t215 + t325;
t472 = -t242 * t394 + t281 * t482;
t469 = t281 * (t364 + t411) - t241 * t394 - t400;
t468 = t412 * t281 + t178 + t402;
t404 = -Icges(5,2) * t424 + t142 - t221;
t408 = t237 * t264 + t138;
t467 = -t284 * t404 - t286 * t408;
t406 = -Icges(6,2) * t424 + t140 - t220;
t410 = t235 * t264 + t136;
t466 = -t284 * t406 - t286 * t410;
t465 = t172 / 0.2e1;
t464 = t173 / 0.2e1;
t431 = t228 * t281;
t428 = t230 * t281;
t409 = -t235 * t265 - t137;
t407 = -t237 * t265 - t139;
t405 = -t231 * t265 + t141;
t403 = -t233 * t265 + t143;
t398 = -t231 + t236;
t397 = t235 + t342;
t396 = -t233 + t238;
t395 = t237 + t343;
t386 = t287 * t451;
t65 = t386 + t472;
t389 = t65 * t457;
t388 = pkin(2) * t420;
t378 = t242 * t393;
t371 = -pkin(3) - t453;
t370 = -t394 / 0.2e1;
t369 = t394 / 0.2e1;
t368 = -t393 / 0.2e1;
t367 = t393 / 0.2e1;
t154 = -t183 - t457;
t355 = -pkin(4) * t422 - t241 * t265;
t352 = t247 - t385;
t351 = -pkin(4) * t390 - t203;
t246 = rSges(2,1) * t287 - rSges(2,2) * t285;
t243 = rSges(2,1) * t285 + rSges(2,2) * t287;
t317 = -t378 - t385;
t64 = t281 * t354 + t317;
t346 = -t264 * t65 - t265 * t64;
t337 = t145 * t264 + t147 * t265;
t326 = -t450 + (-rSges(6,1) - pkin(4)) * t284;
t315 = t265 * t326;
t314 = -t284 * t405 + t286 * t409;
t313 = -t284 * t403 + t286 * t407;
t101 = -t480 + t505;
t311 = (-t284 * t397 + t286 * t398) * t281;
t310 = (-t284 * t395 + t286 * t396) * t281;
t102 = t266 + t483;
t111 = t264 * t371 + t261 + t399 - t457;
t301 = t363 * t393 + t352;
t126 = t154 * t281 - t385;
t127 = t281 * t481 + t386;
t293 = (t126 * (-t259 - t266) + t127 * t154) * t281;
t292 = (t64 * (-t225 - t262 - t266) - t389 + (t64 * (-rSges(5,3) - pkin(7)) + t65 * t371) * t264) * t281;
t291 = ((t502 * t264 + ((t530 + t553) * t265 + t488 + t542 - t552) * t265) * qJD(4) + t521) * t367 + (-t538 * qJD(4) + t540 * t284 + t541 * t286) * t281 + (t507 + t484) * t465 + (t508 + t485) * t464 + (((t265 * t504 + t486 - t502) * t265 + (t264 * t504 + t487 - t543) * t264) * qJD(4) + t515 - t518) * t370 + (t510 + t511) * t369 + (Icges(3,3) + Icges(4,3) - t537) * t280 + (t509 + t512 + t514) * t368;
t50 = t281 * t332 + t301;
t51 = t386 + t469;
t290 = (t50 * (-t266 + t506) + t51 * t505) * t281 + t51 * qJD(4) * t315;
t175 = t281 * t183;
t171 = t242 * t265;
t169 = t242 * t264;
t168 = t241 * t264;
t152 = t193 * t281 + t386;
t72 = qJD(4) * t337 + qJD(3);
t44 = qJD(3) + (t264 * t412 + t265 * t411) * qJD(4);
t25 = t145 * t172 - t147 * t173 + qJDD(3) + (t264 * t96 + t265 * t94) * qJD(4);
t5 = qJDD(3) - t411 * t173 + t412 * t172 + (t264 * t454 + t265 * t455) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t291 + (t490 * (t193 + t278) + t492 * (-t192 - t458) + (-t177 - t386 + t152) * t151) * m(3) + ((t243 ^ 2 + t246 ^ 2) * qJDD(1) + g(1) * t243 - g(2) * t246) * m(2) + (t50 * (t327 - t386) + t290 + t498 * (t102 + t278) + t499 * (t101 - t458) + (t352 - t301 + t50 + t388 + t468) * t51) * m(6) + (t64 * (t380 - t386) + t292 + (-t385 - t317 + t64 + t388 + t477) * t65 + t496 * (t482 + t278) + t497 * (t111 - t458)) * m(5) + (t126 * (t212 - t386) + t293 + t494 * (t481 + t278) + t495 * (t154 - t458) + (t126 + t175 + t388) * t127) * m(4); t291 + (t290 + (t241 * t393 + t281 * t457 + t247 - t330 + t468) * t51 + (t327 + t469) * t50 + t498 * t102 + t499 * t101) * m(6) + (t389 * t281 + t292 + (t378 + t477) * t65 + (t380 + t472) * t64 + t496 * t482 + t497 * t111) * m(5) + (t127 * t175 - (-t126 * t481 - t127 * t457) * t281 + t126 * t212 + t293 + t494 * t481 + t495 * t154) * m(4) + (-t151 * t177 - t152 * t434 + (t151 * t281 + t490) * t193 + (t152 * t281 - t492) * t192) * m(3); m(4) * qJDD(3) + m(5) * t25 + m(6) * t5 + (-m(4) - m(5) - m(6)) * g(3); t525 * t465 + t524 * t464 + (t510 * t281 + t507 * t280 + t487 * t173 + t486 * t172 + (t500 * t264 + t501 * t265) * qJD(4)) * t264 / 0.2e1 - (t509 * t281 + t508 * t280 + t489 * t173 + t488 * t172 + (t491 * t264 + t493 * t265) * qJD(4)) * t265 / 0.2e1 + (t484 * t264 - t485 * t265) * t280 / 0.2e1 - (((t395 + t397) * t286 + (t396 + t398) * t284) * t281 + (((-t404 - t406) * t265 + (t403 + t405) * t264) * t286 + ((t408 + t410) * t265 + (t407 + t409) * t264) * t284) * qJD(4)) * t281 / 0.2e1 + ((t281 * t484 - t512) * t265 + (t281 * t485 + t511) * t264) * t281 / 0.2e1 + t515 * t427 / 0.2e1 + t514 * t423 / 0.2e1 + ((-t394 * t429 + t428) * t264 + (t310 + (-t467 * t265 + (t430 + t313) * t264) * qJD(4)) * t265 + (-t394 * t432 + t431) * t264 + (t311 + (-t466 * t265 + (t433 + t314) * t264) * qJD(4)) * t265) * t370 + ((t281 * t486 + t501) * t265 + (t281 * t487 + t500) * t264) * t369 + ((t281 * t488 + t493) * t265 + (t281 * t489 + t491) * t264) * t368 + ((-t393 * t430 - t428) * t265 + (t310 + (t313 * t264 + (t429 - t467) * t265) * qJD(4)) * t264 + (-t393 * t433 - t431) * t265 + (t311 + (t314 * t264 + (t432 - t466) * t265) * qJD(4)) * t264) * t367 + (-(t169 * t64 - t171 * t65) * t281 - (t72 * (-t169 * t264 - t171 * t265) + t346 * t245) * qJD(4) + t25 * t337 + t72 * ((t94 + t129) * t265 + (-t147 * t281 + t96) * t264) + t346 * t204 + ((-t281 * t65 - t26) * t265 + (t281 * t64 - t27) * t264) * t242 + g(1) * t171 + g(2) * t169 - g(3) * t245) * m(5) + (-(t50 * t168 + t355 * t51) * t281 - g(3) * t478 - g(1) * t315 + (t14 * t363 + t50 * t351 + t5 * t411 + t44 * t455 + (t363 * t51 + t412 * t44) * t281 - (t355 * t44 - t478 * t50) * qJD(4)) * t265 + (t15 * t363 + t51 * t351 + t5 * t412 + t44 * t454 + (t50 * t241 - t411 * t44) * t281 - (-t51 * t478 + (-pkin(4) * t425 - t168) * t44) * qJD(4) - g(2) * t326) * t264) * m(6); (t264 * t499 - t265 * t498) * m(6);];
tau = t1;
