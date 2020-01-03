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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:58:54
% EndTime: 2020-01-03 11:59:11
% DurationCPUTime: 13.67s
% Computational Cost: add. (13377->517), mult. (9752->607), div. (0->0), fcn. (7457->8), ass. (0->299)
t279 = qJ(1) + qJ(2);
t266 = pkin(8) + t279;
t255 = sin(t266);
t256 = cos(t266);
t283 = cos(qJ(4));
t410 = t256 * t283;
t281 = sin(qJ(4));
t411 = t256 * t281;
t132 = Icges(6,4) * t410 - Icges(6,2) * t411 + Icges(6,6) * t255;
t134 = Icges(5,4) * t410 - Icges(5,2) * t411 + Icges(5,6) * t255;
t561 = t132 + t134;
t568 = Icges(5,5) + Icges(6,5);
t567 = -Icges(5,6) - Icges(6,6);
t566 = Icges(5,3) + Icges(6,3);
t214 = Icges(6,4) * t411;
t136 = Icges(6,1) * t410 + Icges(6,5) * t255 - t214;
t215 = Icges(5,4) * t411;
t138 = Icges(5,1) * t410 + Icges(5,5) * t255 - t215;
t564 = t136 + t138;
t226 = Icges(6,5) * t283 - Icges(6,6) * t281;
t228 = Icges(5,5) * t283 - Icges(5,6) * t281;
t565 = t226 + t228;
t269 = Icges(6,4) * t283;
t334 = -Icges(6,2) * t281 + t269;
t316 = t334 * t255;
t131 = -Icges(6,6) * t256 + t316;
t270 = Icges(5,4) * t283;
t335 = -Icges(5,2) * t281 + t270;
t133 = -Icges(5,6) * t256 + t255 * t335;
t562 = t131 + t133;
t435 = Icges(6,4) * t281;
t234 = Icges(6,1) * t283 - t435;
t318 = t234 * t255;
t135 = -Icges(6,5) * t256 + t318;
t436 = Icges(5,4) * t281;
t236 = Icges(5,1) * t283 - t436;
t319 = t236 * t255;
t137 = -Icges(5,5) * t256 + t319;
t560 = t135 + t137;
t229 = Icges(6,2) * t283 + t435;
t231 = Icges(5,2) * t283 + t436;
t559 = t229 + t231;
t563 = t561 * t281;
t475 = Icges(5,1) * t281 + t270;
t476 = Icges(6,1) * t281 + t269;
t557 = t475 + t476;
t509 = t255 * t565 - t256 * t566;
t508 = t255 * t566 + t410 * t568 + t411 * t567;
t558 = t234 + t236;
t538 = -t283 * t564 + t563;
t413 = t255 * t283;
t511 = t560 * t413;
t556 = t564 * t413;
t510 = t562 * t411;
t225 = Icges(6,5) * t281 + Icges(6,6) * t283;
t227 = Icges(5,5) * t281 + Icges(5,6) * t283;
t555 = t225 + t227;
t550 = -t281 * t559 + t283 * t557;
t278 = qJD(1) + qJD(2);
t554 = t557 * qJD(4) - t278 * t568;
t553 = t559 * qJD(4) + t278 * t567;
t414 = t255 * t281;
t494 = -t256 * t509 - t414 * t562 + t511;
t493 = t256 * t508 + t414 * t561 - t556;
t492 = -t255 * t509 - t410 * t560 + t510;
t491 = t255 * t508 - t256 * t538;
t552 = (t334 + t335) * qJD(4);
t551 = t558 * qJD(4);
t549 = t281 * t557 + t283 * t559;
t504 = t560 * t283;
t503 = t562 * t281;
t490 = t281 * t560 + t283 * t562;
t317 = t335 * t278;
t548 = t255 * t317 + t256 * t553 + t278 * t316;
t412 = t256 * t278;
t547 = -t255 * t553 + t256 * t317 + t334 * t412;
t546 = (t318 + t319) * t278 + t554 * t256;
t545 = -t255 * t554 + t412 * t558;
t280 = -qJ(5) - pkin(7);
t544 = rSges(6,3) - t280;
t152 = t225 * t256;
t154 = t227 * t256;
t513 = t255 * t550 - t152 - t154;
t416 = t227 * t255;
t417 = t225 * t255;
t512 = t410 * t557 - t411 * t559 + t416 + t417;
t314 = t226 * t278;
t315 = t228 * t278;
t543 = t314 + t315;
t542 = -t503 + t504;
t489 = t281 * t564 + t283 * t561;
t541 = t555 * qJD(4) - t278 * t566;
t218 = rSges(6,2) * t411;
t274 = t283 * pkin(4);
t259 = t274 + pkin(3);
t540 = rSges(6,1) * t410 + t256 * t259 - t218;
t539 = -qJD(4) * t565 + t550 * t278;
t537 = -t491 * t255 - t256 * t492;
t536 = -t255 * t493 - t494 * t256;
t535 = qJD(4) * t549 - t278 * t555 + t281 * t552 - t283 * t551;
t534 = t513 * t278;
t448 = pkin(3) * t256;
t404 = t448 - t540 + (pkin(7) - t544) * t255;
t533 = -t255 * t541 + t256 * t543 - t278 * t542;
t532 = t255 * t543 + t256 * t541 - t278 * t538;
t267 = sin(t279);
t257 = pkin(2) * t267;
t344 = -rSges(4,1) * t255 - rSges(4,2) * t256;
t531 = t344 - t257;
t268 = cos(t279);
t478 = rSges(3,1) * t267 + rSges(3,2) * t268;
t171 = t478 * t278;
t282 = sin(qJ(1));
t443 = pkin(1) * qJD(1);
t376 = t282 * t443;
t146 = t376 + t171;
t530 = t512 * t278;
t529 = qJD(4) * t489 - t278 * t508 - t281 * t548 + t283 * t546;
t528 = qJD(4) * t490 - t278 * t509 + t281 * t547 - t283 * t545;
t382 = qJD(4) * t255;
t208 = pkin(7) * t412;
t415 = t255 * t278;
t148 = pkin(3) * t415 - t208;
t380 = qJD(4) * t278;
t167 = -qJDD(4) * t255 - t256 * t380;
t445 = rSges(5,1) * t283;
t245 = -rSges(5,2) * t281 + t445;
t198 = t245 * qJD(4);
t242 = rSges(5,1) * t281 + rSges(5,2) * t283;
t277 = qJDD(1) + qJDD(2);
t273 = t282 * pkin(1);
t284 = cos(qJ(1));
t286 = qJD(1) ^ 2;
t426 = pkin(1) * qJDD(1);
t355 = -t273 * t286 + t284 * t426;
t449 = pkin(2) * t277;
t276 = t278 ^ 2;
t450 = pkin(2) * t276;
t292 = -t267 * t450 + t268 * t449 + t355;
t219 = rSges(5,2) * t411;
t142 = rSges(5,1) * t410 + rSges(5,3) * t255 - t219;
t181 = pkin(7) * t255 + t448;
t395 = t142 + t181;
t379 = qJD(4) * t281;
t369 = t256 * t379;
t407 = t278 * t283;
t311 = t255 * t407 + t369;
t378 = qJD(4) * t283;
t368 = t256 * t378;
t408 = t278 * t281;
t375 = t255 * t408;
t393 = rSges(5,2) * t375 + rSges(5,3) * t412;
t93 = rSges(5,1) * t311 + rSges(5,2) * t368 - t393;
t27 = -t198 * t382 + t167 * t242 + (-t148 - t93) * t278 + t395 * t277 + t292;
t527 = -g(2) + t27;
t217 = rSges(5,1) * t413;
t350 = -rSges(5,2) * t414 + t217;
t140 = -rSges(5,3) * t256 + t350;
t168 = -qJDD(4) * t256 + t255 * t380;
t251 = t255 * pkin(3);
t180 = -pkin(7) * t256 + t251;
t275 = t284 * pkin(1);
t383 = t275 * t286 + t282 * t426;
t357 = t267 * t449 + t268 * t450 + t383;
t389 = pkin(3) * t412 + pkin(7) * t415;
t321 = t180 * t277 + t278 * t389 + t357;
t381 = qJD(4) * t256;
t310 = -t255 * t378 - t256 * t408;
t370 = t255 * t379;
t374 = t256 * t407;
t391 = rSges(5,1) * t374 + rSges(5,3) * t415;
t95 = -rSges(5,1) * t370 + rSges(5,2) * t310 + t391;
t26 = t140 * t277 - t168 * t242 + t198 * t381 + t278 * t95 + t321;
t526 = -g(3) + t26;
t244 = rSges(6,1) * t283 - rSges(6,2) * t281;
t197 = t244 * qJD(4);
t440 = t283 * rSges(6,2);
t241 = rSges(6,1) * t281 + t440;
t359 = pkin(4) * t281 + t241;
t240 = t256 * t280;
t507 = -rSges(6,1) * t413 - t255 * t259 - t240;
t485 = -rSges(6,2) * t414 - rSges(6,3) * t256 - t507;
t405 = -t180 + t485;
t406 = t283 * qJD(4) ^ 2;
t377 = qJD(5) * t256;
t484 = rSges(6,1) * t374 + rSges(6,3) * t415 + t259 * t412;
t446 = -t377 + (-pkin(4) * t379 - t278 * t280) * t255 - t389 - rSges(6,1) * t370 + rSges(6,2) * t310 + t484;
t14 = -qJDD(5) * t255 + t446 * t278 + t405 * t277 - t359 * t168 + (pkin(4) * t406 + qJD(4) * t197 - qJD(5) * t278) * t256 + t321;
t498 = t14 - g(3);
t247 = qJD(5) * t255;
t373 = t181 - t404;
t394 = rSges(6,2) * t375 + rSges(6,3) * t412;
t480 = t247 + t394;
t447 = -pkin(4) * t369 - t208 - (t240 + (-pkin(3) + t259) * t255) * t278 - rSges(6,1) * t311 - rSges(6,2) * t368 + t480;
t15 = -t197 * t382 - qJDD(5) * t256 + t167 * t241 + (t167 * t281 - t255 * t406) * pkin(4) + t373 * t277 + (-t148 + t247 + t447) * t278 + t292;
t497 = t15 - g(2);
t206 = rSges(4,1) * t412;
t525 = -t277 * t344 + t278 * (-rSges(4,2) * t415 + t206) + t357 - g(3);
t248 = t255 * rSges(4,2);
t179 = rSges(4,1) * t256 - t248;
t524 = t179 * t277 + t276 * t344 - g(2) + t292;
t252 = t267 * rSges(3,2);
t409 = t268 * t278;
t172 = rSges(3,1) * t409 - t252 * t278;
t523 = t172 * t278 + t277 * t478 - g(3) + t383;
t188 = rSges(3,1) * t268 - t252;
t522 = -t171 * t278 + t188 * t277 - g(2) + t355;
t521 = qJD(4) * t536 + t534;
t520 = qJD(4) * t537 - t530;
t519 = qJD(4) * t542 + t281 * t545 + t283 * t547;
t518 = qJD(4) * t538 + t281 * t546 + t283 * t548;
t517 = t255 * t539 + t256 * t535;
t516 = -t255 * t535 + t256 * t539;
t387 = t476 + t334;
t388 = t229 - t234;
t515 = (t281 * t387 + t283 * t388) * t278;
t385 = t475 + t335;
t386 = t231 - t236;
t514 = (t281 * t385 + t283 * t386) * t278;
t173 = t278 * t181;
t239 = pkin(2) * t409;
t358 = t239 - t377;
t506 = t278 * t404 - t173 + t358 + t484;
t124 = t278 * t142;
t347 = -t242 * t382 + t239;
t505 = -t124 - t173 - t347 + t239 + t389 + t391;
t164 = t242 * t255;
t166 = t242 * t256;
t360 = t180 + t257;
t352 = t140 + t360;
t371 = t242 * t381;
t63 = t278 * t352 + t371 + t376;
t265 = t284 * t443;
t64 = t278 * t395 + t265 + t347;
t502 = -t164 * t63 - t166 * t64;
t501 = -t359 * t382 + t358;
t500 = t255 * t533 - t256 * t528;
t499 = t255 * t532 + t256 * t529;
t496 = t255 * t528 + t256 * t533;
t495 = -t255 * t529 + t256 * t532;
t487 = t278 * (t360 + t405);
t486 = t278 * t531;
t222 = pkin(4) * t411;
t483 = rSges(6,1) * t411 + rSges(6,2) * t410 + t222;
t170 = t278 * t179;
t482 = t239 + t170;
t481 = t244 + t274;
t479 = -t251 - t257;
t458 = t167 / 0.2e1;
t457 = t168 / 0.2e1;
t451 = rSges(5,3) + pkin(7);
t441 = t278 * t64;
t403 = t255 * t476 + t131;
t402 = t256 * t476 + t132;
t401 = t255 * t475 + t133;
t400 = t256 * t475 + t134;
t399 = -t229 * t255 + t135;
t398 = -Icges(6,2) * t410 + t136 - t214;
t397 = -t231 * t255 + t137;
t396 = -Icges(5,2) * t410 + t138 - t215;
t372 = t208 + t393;
t365 = -t382 / 0.2e1;
t364 = t382 / 0.2e1;
t363 = -t381 / 0.2e1;
t362 = t381 / 0.2e1;
t147 = t188 * t278 + t265;
t354 = -pkin(4) * t414 - t241 * t255;
t349 = -t247 + t376;
t258 = pkin(2) * t268;
t150 = t179 + t258;
t291 = t265 + t501;
t50 = t278 * t373 + t291;
t348 = t50 * t359;
t246 = rSges(2,1) * t284 - rSges(2,2) * t282;
t243 = rSges(2,1) * t282 + rSges(2,2) * t284;
t339 = -t255 * t64 + t256 * t63;
t327 = t140 * t255 + t142 * t256;
t320 = -t440 + (-rSges(6,1) - pkin(4)) * t281;
t305 = t255 * t320;
t296 = t281 * t399 + t283 * t403;
t295 = t281 * t398 + t283 * t402;
t294 = t281 * t397 + t283 * t401;
t293 = t281 * t396 + t283 * t400;
t100 = t257 + t485;
t110 = -t256 * t451 + t350 - t479;
t101 = t255 * t544 + t258 + t540;
t111 = -t219 + t258 + (pkin(3) + t445) * t256 + t451 * t255;
t121 = t376 - t486;
t122 = t265 + t482;
t290 = (-t121 * t248 + t122 * t531) * t278;
t289 = -t512 * t167 / 0.2e1 - t489 * t458 + ((((t509 - t538) * t256 - t511 - t491) * t256 + ((t509 - t563) * t255 + (t503 + t504) * t256 - t510 + t492 + t556) * t255) * qJD(4) + t534) * t364 + (qJD(4) * t550 + t281 * t551 + t283 * t552) * t278 + (t513 + t490) * t457 + (t516 + t519) * t363 + ((((-t504 + t508) * t256 + t510 - t493) * t256 + ((t503 + t508) * t255 - t511 + t494) * t255) * qJD(4) + t520 + t530) * t362 + (Icges(3,3) + Icges(4,3) + t549) * t277 + (t517 + t518 + t521) * t365;
t288 = (t64 * (-t217 + t479) - t63 * t219) * t278 + t502 * qJD(4);
t49 = t359 * t381 + t349 + t487;
t287 = (t50 * (-t257 + t507) + t49 * (-t255 * t280 - t218)) * t278 + (t256 * t320 * t50 + t305 * t49) * qJD(4);
t71 = qJD(4) * t327 + qJD(3);
t44 = qJD(3) + (t255 * t405 - t256 * t404) * qJD(4);
t25 = -t140 * t167 - t142 * t168 + qJDD(3) + (t255 * t95 - t256 * t93) * qJD(4);
t5 = qJDD(3) + t404 * t168 - t405 * t167 + (t255 * t446 + t256 * t447) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t289 + (t522 * (t188 + t275) + t523 * (t273 + t478) + (-t147 + t172 + t265) * t146) * m(3) + ((t243 ^ 2 + t246 ^ 2) * qJDD(1) - g(2) * t246 - g(3) * t243) * m(2) + (t50 * (-t349 + t394) + t287 + t497 * (t275 + t101) + t498 * (t100 + t273) + (t50 - t291 + t265 + t506) * t49) * m(6) + (t64 * (t372 - t376) + t288 + t527 * (t275 + t111) + t526 * (t273 + t110) + (t64 + t505) * t63) * m(5) + (-t122 * t376 + t290 + t524 * (t150 + t275) + t525 * (t273 - t531) + (t122 - t170 + t206) * t121) * m(4); t289 + (t348 * t381 + t287 + (-t247 + t487 + t480) * t50 + (-t501 + t506) * t49 + t497 * t101 + t498 * t100) * m(6) + (t352 * t441 + t288 + (t371 + t372) * t64 + t505 * t63 + t527 * t111 + t526 * t110) * m(5) + (-t122 * t486 + t290 + t524 * t150 - t525 * t531 + (-t482 + t206 + t239) * t121) * m(4) + (t146 * t172 - t147 * t171 + (-t146 * t278 + t522) * t188 + (t147 * t278 + t523) * t478) * m(3); m(4) * qJDD(3) + m(5) * t25 + m(6) * t5 + (-m(4) - m(5) - m(6)) * g(1); t537 * t458 + t536 * t457 - (t517 * t278 - t512 * t277 + t492 * t168 + t491 * t167 + (t499 * t255 + t500 * t256) * qJD(4)) * t255 / 0.2e1 - (t516 * t278 + t513 * t277 + t494 * t168 + t493 * t167 + (t495 * t255 + t496 * t256) * qJD(4)) * t256 / 0.2e1 + (t255 * t489 - t256 * t490) * t277 / 0.2e1 - (((t385 + t387) * t283 + (-t386 - t388) * t281) * t278 + (((-t397 - t399) * t256 + (t396 + t398) * t255) * t283 + ((t401 + t403) * t256 + (-t400 - t402) * t255) * t281) * qJD(4)) * t278 / 0.2e1 + ((t278 * t489 - t519) * t256 + (t278 * t490 - t518) * t255) * t278 / 0.2e1 + t521 * t415 / 0.2e1 - t520 * t412 / 0.2e1 + ((-t278 * t491 + t500) * t256 + (t278 * t492 + t499) * t255) * t365 + ((t154 * t382 - t315) * t255 + (t514 + (-t294 * t256 + (-t416 + t293) * t255) * qJD(4)) * t256 + (t152 * t382 - t314) * t255 + (t515 + (-t296 * t256 + (-t417 + t295) * t255) * qJD(4)) * t256) * t364 + ((-t278 * t493 + t496) * t256 + (t278 * t494 + t495) * t255) * t363 + ((-t381 * t416 - t315) * t256 + (-t514 + (-t293 * t255 + (t154 + t294) * t256) * qJD(4)) * t255 + (-t381 * t417 - t314) * t256 + (-t515 + (-t295 * t255 + (t152 + t296) * t256) * qJD(4)) * t255) * t362 + (-t502 * t278 - (t71 * (-t164 * t255 - t166 * t256) + t339 * t245) * qJD(4) + t25 * t327 + t71 * ((t140 * t278 - t93) * t256 + (t95 - t124) * t255) + t339 * t198 + ((t26 - t441) * t256 + (-t278 * t63 - t27) * t255) * t242 - g(1) * t245 + g(2) * t164 - g(3) * t166) * m(5) + (-g(1) * t481 - g(3) * t483 - g(2) * t305 + t14 * t222 + (-t5 * t404 + t44 * t447 + t14 * t241 + t49 * t197 + (t405 * t44 - t348) * t278) * t256 + (t5 * t405 + t44 * t446 - t15 * t359 + t50 * (-pkin(4) * t378 - t197) + (-t359 * t49 + t404 * t44) * t278) * t255 - (t354 * t49 - t483 * t50) * t278 - ((t244 * t49 - t44 * t483) * t256 + (t354 * t44 - t481 * t50) * t255) * qJD(4)) * m(6); (-t255 * t498 - t256 * t497) * m(6);];
tau = t1;
