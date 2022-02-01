% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:12
% EndTime: 2022-01-20 11:02:34
% DurationCPUTime: 9.04s
% Computational Cost: add. (18296->573), mult. (12244->729), div. (0->0), fcn. (9476->10), ass. (0->350)
t303 = pkin(9) + qJ(4);
t297 = qJ(5) + t303;
t291 = sin(t297);
t306 = qJ(1) + qJ(2);
t299 = cos(t306);
t473 = t291 * t299;
t425 = rSges(6,2) * t473;
t298 = sin(t306);
t292 = cos(t297);
t471 = t292 * t299;
t443 = rSges(6,1) * t471 + t298 * rSges(6,3);
t155 = -t425 + t443;
t308 = cos(pkin(9));
t293 = t308 * pkin(3) + pkin(2);
t255 = t299 * t293;
t309 = -pkin(7) - qJ(3);
t270 = t298 * t309;
t438 = t270 - t255;
t296 = cos(t303);
t503 = pkin(4) * t296;
t238 = t293 + t503;
t302 = pkin(8) - t309;
t540 = t299 * t238 + t302 * t298;
t530 = t438 + t540;
t543 = -t530 - t155;
t271 = t299 * t309;
t441 = t298 * t293 + t271;
t125 = t238 * t298 - t299 * t302 - t441;
t305 = qJD(1) + qJD(2);
t542 = t125 * t305;
t467 = t296 * t299;
t364 = rSges(5,1) * t467 + t298 * rSges(5,3);
t295 = sin(t303);
t469 = t295 * t299;
t170 = -rSges(5,2) * t469 + t364;
t541 = t170 - t438;
t272 = Icges(6,4) * t292;
t211 = Icges(6,1) * t291 + t272;
t372 = -Icges(6,2) * t291 + t272;
t539 = t211 + t372;
t310 = sin(qJ(1));
t495 = pkin(1) * qJD(1);
t423 = t310 * t495;
t233 = rSges(3,1) * t298 + rSges(3,2) * t299;
t475 = t233 * t305;
t182 = -t423 - t475;
t304 = qJD(4) + qJD(5);
t496 = rSges(6,2) * t292;
t424 = t304 * t496;
t500 = rSges(6,1) * t291;
t538 = -t304 * t500 - t424;
t537 = -t438 - t543;
t536 = 2 * qJD(4);
t535 = rSges(5,2) * t295;
t498 = rSges(4,2) * sin(pkin(9));
t426 = t299 * t498;
t276 = t298 * qJ(3);
t529 = t299 * pkin(2) + t276;
t501 = rSges(4,1) * t308;
t531 = -t298 * rSges(4,3) - t299 * t501;
t348 = -t426 + t529 - t531;
t534 = t348 * t305;
t213 = t496 + t500;
t230 = t298 * t304;
t435 = qJD(4) * t298;
t413 = t295 * t435;
t244 = pkin(4) * t413;
t533 = -t230 * t213 - t244;
t277 = t299 * qJ(3);
t232 = pkin(2) * t298 - t277;
t146 = t232 - t441;
t202 = t305 * t232;
t532 = -t305 * t146 + t202;
t273 = Icges(5,4) * t296;
t373 = -Icges(5,2) * t295 + t273;
t226 = Icges(5,1) * t295 + t273;
t231 = t299 * t304;
t474 = t291 * t298;
t245 = rSges(6,2) * t474;
t463 = t299 * t305;
t449 = rSges(6,3) * t463 + t305 * t245;
t466 = t298 * t305;
t100 = -t299 * t424 + (-t291 * t231 - t292 * t466) * rSges(6,1) + t449;
t462 = t302 * t305;
t247 = t299 * t462;
t434 = qJD(4) * t299;
t412 = t295 * t434;
t389 = pkin(4) * t412;
t445 = t238 - t293;
t103 = -t389 + t247 + (-t445 * t298 + t271) * t305;
t497 = rSges(6,2) * t291;
t499 = rSges(6,1) * t292;
t214 = -t497 + t499;
t187 = t214 * t304;
t197 = t305 * t231;
t260 = qJ(3) * t463;
t313 = qJD(1) ^ 2;
t505 = pkin(1) * t310;
t432 = t313 * t505;
t437 = qJD(3) * t305;
t274 = qJD(3) * t298;
t440 = t260 + t274;
t365 = t305 * (-pkin(2) * t466 + t440) + t298 * t437 - t432;
t502 = pkin(2) - t293;
t362 = t305 * (-t260 + (t502 * t298 - t271) * t305) + t365;
t433 = (qJD(4) ^ 2) * t503;
t30 = -t298 * t433 - t187 * t230 - t197 * t213 + (t100 + t103 - t389) * t305 + t362;
t415 = t538 * t298 - t305 * t425;
t101 = t443 * t305 + t415;
t250 = t305 * t270;
t442 = t298 * t462 - t244;
t104 = t445 * t463 + t250 + t442;
t196 = t305 * t230;
t311 = cos(qJ(1));
t301 = t311 * pkin(1);
t431 = t313 * t301;
t385 = t299 * t437 - t431;
t275 = qJD(3) * t299;
t171 = t529 * t305 - t275;
t461 = t250 - (-t502 * t299 - t276) * t305 - t171;
t31 = -t299 * t433 - t187 * t231 + t196 * t213 + (-t101 - t104 + t244 + t461) * t305 + t385;
t401 = -t238 - t499;
t444 = t299 * rSges(6,3) + t245;
t472 = t292 * t298;
t154 = rSges(6,1) * t472 - t444;
t344 = -t213 * t231 + t274 - t389;
t456 = t146 - t232;
t59 = -t423 + (-t125 - t154 + t456) * t305 + t344;
t494 = t305 * t59;
t504 = pkin(4) * t295;
t422 = t311 * t495;
t383 = -t275 + t422;
t361 = t383 + t533;
t60 = t537 * t305 + t361;
t528 = (t31 * t401 + (-t59 * rSges(6,3) + t60 * t401) * t305) * t298 + (t31 * t302 - t30 * t497 + t60 * (-qJD(4) * t504 + t538) + t401 * t494) * t299 + t537 * t494;
t429 = t298 * t501;
t267 = t298 * t498;
t439 = t299 * rSges(4,3) + t267;
t178 = t429 - t439;
t446 = rSges(4,3) * t463 + t305 * t267;
t527 = t305 * t178 + t202 + t446;
t142 = t305 * t154;
t526 = t142 + t247 - t344 + t449 + t532 + t542;
t488 = Icges(6,4) * t291;
t212 = Icges(6,1) * t292 - t488;
t355 = t212 * t299;
t153 = Icges(6,5) * t298 + t355;
t208 = Icges(6,5) * t292 - Icges(6,6) * t291;
t207 = Icges(6,5) * t291 + Icges(6,6) * t292;
t336 = Icges(6,3) * t305 - t207 * t304;
t353 = t372 * t299;
t151 = Icges(6,6) * t298 + t353;
t484 = t151 * t291;
t525 = -t208 * t466 + t336 * t299 + (-t153 * t292 + t484) * t305;
t351 = t208 * t299;
t150 = Icges(6,4) * t472 - Icges(6,2) * t474 - Icges(6,6) * t299;
t243 = Icges(6,4) * t474;
t152 = Icges(6,1) * t472 - Icges(6,5) * t299 - t243;
t371 = t150 * t291 - t152 * t292;
t524 = t336 * t298 + (t351 + t371) * t305;
t223 = Icges(5,5) * t296 - Icges(5,6) * t295;
t222 = Icges(5,5) * t295 + Icges(5,6) * t296;
t333 = Icges(5,3) * t305 - t222 * qJD(4);
t489 = Icges(5,4) * t295;
t227 = Icges(5,1) * t296 - t489;
t356 = t227 * t299;
t166 = Icges(5,5) * t298 + t356;
t354 = t373 * t299;
t164 = Icges(5,6) * t298 + t354;
t482 = t164 * t295;
t368 = -t166 * t296 + t482;
t523 = -t223 * t466 + t333 * t299 + t368 * t305;
t352 = t223 * t299;
t470 = t295 * t298;
t253 = Icges(5,4) * t470;
t468 = t296 * t298;
t165 = Icges(5,1) * t468 - Icges(5,5) * t299 - t253;
t163 = Icges(5,4) * t468 - Icges(5,2) * t470 - Icges(5,6) * t299;
t483 = t163 * t295;
t369 = -t165 * t296 + t483;
t522 = t333 * t298 + (t352 + t369) * t305;
t209 = Icges(6,2) * t292 + t488;
t367 = t209 * t291 - t211 * t292;
t521 = t208 * t304 + t367 * t305;
t224 = Icges(5,2) * t296 + t489;
t366 = t224 * t295 - t226 * t296;
t520 = t223 * qJD(4) + t366 * t305;
t161 = Icges(5,5) * t468 - Icges(5,6) * t470 - Icges(5,3) * t299;
t72 = -t299 * t161 - t369 * t298;
t428 = rSges(5,1) * t468;
t169 = -rSges(5,2) * t470 - t299 * rSges(5,3) + t428;
t145 = t305 * t169;
t436 = qJD(4) * t296;
t421 = rSges(5,2) * t436;
t357 = rSges(5,3) * t463 - t299 * t421 + t466 * t535;
t519 = -rSges(5,1) * t412 + t145 + t274 + t357 + t532;
t12 = t100 * t231 + t101 * t230 + t154 * t197 - t155 * t196 + ((t103 + t542) * t299 + (-t305 * t530 + t104) * t298) * qJD(4);
t180 = t213 * t298;
t181 = t213 * t299;
t58 = t154 * t230 + t155 * t231 + (t125 * t298 + t299 * t530) * qJD(4);
t518 = -t59 * (t180 * t305 - t231 * t214) - t58 * (-t230 * t180 - t181 * t231) - t60 * (-t305 * t181 - t214 * t230) + t12 * (t298 * t154 + t299 * t155);
t517 = (-t224 * t299 + t166) * t298 - (-Icges(5,2) * t468 + t165 - t253) * t299;
t516 = (-t209 * t299 + t153) * t230 - (-Icges(6,2) * t472 + t152 - t243) * t231 + t539 * t305;
t515 = t196 / 0.2e1;
t514 = t197 / 0.2e1;
t513 = -t230 / 0.2e1;
t512 = t230 / 0.2e1;
t511 = -t231 / 0.2e1;
t510 = t231 / 0.2e1;
t509 = t298 / 0.2e1;
t508 = -t299 / 0.2e1;
t507 = -t305 / 0.2e1;
t506 = t305 / 0.2e1;
t228 = rSges(5,1) * t295 + rSges(5,2) * t296;
t381 = -t228 * t434 + t274;
t343 = t381 - t423;
t69 = (-t169 + t456) * t305 + t343;
t493 = t305 * t69;
t480 = t207 * t299;
t90 = -t367 * t298 - t480;
t492 = t90 * t305;
t477 = t222 * t299;
t108 = -t366 * t298 - t477;
t485 = t108 * t305;
t481 = t207 * t298;
t479 = t209 * t304;
t478 = t222 * t298;
t476 = t223 * t305;
t148 = Icges(6,5) * t472 - Icges(6,6) * t474 - Icges(6,3) * t299;
t465 = t299 * t148;
t460 = -t298 * t148 - t152 * t471;
t149 = Icges(6,3) * t298 + t351;
t459 = t298 * t149 + t153 * t471;
t458 = -t298 * t161 - t165 * t467;
t162 = Icges(5,3) * t298 + t352;
t457 = t298 * t162 + t166 * t467;
t448 = -t224 + t227;
t447 = t226 + t373;
t430 = t298 * t101 + (t100 + t142) * t299;
t419 = t295 * t463;
t414 = -rSges(5,1) * t413 - rSges(5,2) * t419 - t298 * t421;
t198 = t228 * t435;
t411 = t466 / 0.2e1;
t410 = t463 / 0.2e1;
t409 = -pkin(2) - t501;
t408 = -t435 / 0.2e1;
t405 = t434 / 0.2e1;
t403 = -t213 - t504;
t235 = t299 * rSges(3,1) - rSges(3,2) * t298;
t338 = Icges(6,5) * t305 - t211 * t304;
t400 = -t150 * t304 + t338 * t298 + t305 * t355;
t399 = -t151 * t304 - t212 * t466 + t338 * t299;
t337 = Icges(6,6) * t305 - t479;
t398 = t152 * t304 + t337 * t298 + t305 * t353;
t397 = t153 * t304 + t337 * t299 - t372 * t466;
t129 = t153 * t472;
t396 = t299 * t149 - t129;
t134 = t166 * t468;
t395 = t299 * t162 - t134;
t394 = -t148 + t484;
t392 = -t161 + t482;
t391 = t539 * t304;
t390 = t212 * t304 - t479;
t388 = t443 + t540;
t200 = rSges(3,1) * t463 - rSges(3,2) * t466;
t384 = t274 - t423;
t382 = -pkin(4) * t436 - t187;
t379 = rSges(5,1) * t296 - t535;
t378 = -t298 * t60 - t299 * t59;
t70 = t305 * t541 - t198 + t383;
t377 = -t298 * t70 - t299 * t69;
t376 = t250 + t275 - t414;
t92 = t150 * t292 + t152 * t291;
t106 = t163 * t296 + t165 * t295;
t107 = t164 * t296 + t166 * t295;
t363 = t275 - t415 - t442;
t73 = -t164 * t470 - t395;
t359 = (t298 * t73 - t299 * t72) * qJD(4);
t74 = -t163 * t469 - t458;
t75 = -t164 * t469 + t457;
t358 = (t298 * t75 - t299 * t74) * qJD(4);
t350 = t371 * t298;
t105 = (t169 * t298 + t170 * t299) * qJD(4);
t349 = -t169 - t441;
t346 = -t208 * t305 + t230 * t480 - t231 * t481;
t342 = t163 * t299 - t164 * t298;
t341 = t409 * t298 + t277 + t439;
t325 = t148 * t305 - t398 * t291 + t400 * t292;
t14 = t524 * t298 + t325 * t299;
t324 = t149 * t305 - t397 * t291 + t399 * t292;
t15 = t525 * t298 + t324 * t299;
t16 = t325 * t298 - t524 * t299;
t17 = t324 * t298 - t525 * t299;
t65 = -t350 - t465;
t66 = -t151 * t474 - t396;
t28 = t230 * t66 - t231 * t65 + t492;
t67 = -t150 * t473 - t460;
t68 = -t151 * t473 + t459;
t91 = -t367 * t299 + t481;
t85 = t91 * t305;
t29 = t230 * t68 - t231 * t67 + t85;
t328 = (-t211 * t299 - t151) * t230 - (-t211 * t298 - t150) * t231 + (-t209 + t212) * t305;
t315 = -t516 * t291 + t328 * t292;
t323 = t207 * t305 - t391 * t291 + t390 * t292;
t42 = t521 * t298 + t323 * t299;
t43 = t323 * t298 - t521 * t299;
t44 = t400 * t291 + t398 * t292;
t45 = t399 * t291 + t397 * t292;
t93 = t151 * t292 + t153 * t291;
t340 = (-t14 * t231 + t15 * t230 + t196 * t67 + t197 * t68 + t305 * t42) * t509 + (-t346 * t298 + t315 * t299) * t513 + (t315 * t298 + t346 * t299) * t510 + (-t16 * t231 + t17 * t230 + t196 * t65 + t197 * t66 + t305 * t43) * t508 + (t328 * t291 + t516 * t292) * t507 + t28 * t411 + t29 * t410 + ((t305 * t68 - t14) * t299 + (t305 * t67 + t15) * t298) * t512 + (t298 * t66 - t299 * t65) * t515 + (t298 * t68 - t299 * t67) * t514 + ((t305 * t66 - t16) * t299 + (t305 * t65 + t17) * t298) * t511 + ((t305 * t93 - t44) * t299 + (t305 * t92 + t45) * t298) * t506;
t339 = (-t447 * t295 + t448 * t296) * t305;
t335 = Icges(5,5) * t305 - qJD(4) * t226;
t334 = Icges(5,6) * t305 - t224 * qJD(4);
t114 = t334 * t299 - t373 * t466;
t116 = -t227 * t466 + t335 * t299;
t322 = -t107 * qJD(4) - t114 * t295 + t116 * t296 + t162 * t305;
t115 = t334 * t298 + t305 * t354;
t117 = t335 * t298 + t305 * t356;
t321 = -t106 * qJD(4) - t115 * t295 + t117 * t296 + t161 * t305;
t204 = t373 * qJD(4);
t205 = t227 * qJD(4);
t320 = -t204 * t295 + t205 * t296 + t222 * t305 + (-t224 * t296 - t226 * t295) * qJD(4);
t120 = (-t178 - t232) * t305 + t384;
t121 = t383 + t534;
t319 = (t120 * t409 * t299 + (t120 * (-rSges(4,3) - qJ(3)) + t121 * t409) * t298) * t305;
t318 = -t517 * t295 + t342 * t296;
t317 = (t69 * (-t364 - t255) + t70 * (-t428 - t441)) * t305;
t109 = -t366 * t299 + t478;
t102 = t109 * t305;
t34 = t359 + t485;
t35 = t102 + t358;
t48 = -t369 * qJD(4) + t115 * t296 + t117 * t295;
t49 = -t368 * qJD(4) + t114 * t296 + t116 * t295;
t54 = t520 * t298 + t320 * t299;
t55 = t320 * t298 - t520 * t299;
t316 = (t102 + ((t73 - t134 + (t162 + t483) * t299 + t458) * t299 + t457 * t298) * qJD(4)) * t405 + (t85 + (t66 + (t150 * t299 + t151 * t298) * t291 + t396 + t460) * t231 + (-t152 * t472 + t465 + t65 + (t150 * t298 - t151 * t299) * t291 + t459) * t230) * t510 + (t92 + t90) * t515 + (t91 + t93) * t514 + (-t492 + (t68 - t350 - t459) * t231 + (t394 * t298 - t129 + t67) * t230 + ((t149 + t371) * t230 + t394 * t231) * t299 + t28) * t513 + (t45 + t42) * t512 + (-t485 + ((t392 * t299 - t457 + t75) * t299 + (t392 * t298 + t395 + t74) * t298) * qJD(4) + t34) * t408 + (t49 + t54) * t435 / 0.2e1 + (-t366 * qJD(4) + t204 * t296 + t205 * t295 + t390 * t291 + t391 * t292) * t305 + (t44 + t43 + t29) * t511 - (t48 + t55 + t35) * t434 / 0.2e1 + ((t106 + t108) * t298 + (t107 + t109) * t299) * qJD(4) * t506;
t237 = t305 * t426;
t206 = t379 * qJD(4);
t195 = t228 * t299;
t194 = t228 * t298;
t183 = t235 * t305 + t422;
t160 = -t200 * t305 - t431;
t159 = -t305 * t475 - t432;
t119 = t364 * t305 + t414;
t118 = (-t296 * t466 - t412) * rSges(5,1) + t357;
t82 = (t531 * t305 - t171 + t237) * t305 + t385;
t81 = t305 * (-t305 * t429 + t446) + t365;
t53 = -t206 * t434 + (-t119 + t198 + t461) * t305 + t385;
t52 = t118 * t305 + (-t206 * t298 - t228 * t463) * qJD(4) + t362;
t1 = [t316 + m(3) * (t160 * (-t233 - t505) + t159 * (t235 + t301) + (-t200 - t422 + t183) * t182) + (t31 * (t444 - t505) + t30 * (t301 + t388) + (t384 + t423 + t526) * t60 + (t363 - t422 + t361) * t59 + t528) * m(6) + (t53 * (t349 - t505) + t69 * (t376 - t422) + t52 * (t301 + t541) + t317 + (-t423 - t343 + t69 + t519) * t70) * m(5) + (t82 * (t341 - t505) + t120 * (t237 - t383) + t81 * (t301 + t348) + t319 + (t260 + t120 + t527) * t121) * m(4); t316 + (t30 * t388 + t31 * t444 + (t274 + t526) * t60 + (t363 - t275 + t533) * t59 + t528) * m(6) + (t53 * t349 + t317 + (-t381 + t519) * t70 + (-t198 - t275 + t376) * t69 + (t52 + t493) * t541) * m(5) + (t82 * t341 + t81 * t348 + t319 + (-t274 + t440 + t527) * t121 + (t237 + t534) * t120) * m(4) + (t159 * t235 - t160 * t233 - t182 * t200 - t183 * t475 - (-t182 * t235 - t183 * t233) * t305) * m(3); 0.2e1 * (t30 * t508 + t31 * t509) * m(6) + 0.2e1 * (t52 * t508 + t53 * t509) * m(5) + 0.2e1 * (t81 * t508 + t82 * t509) * m(4); ((-t434 * t478 - t476) * t299 + (t339 + (t299 * t477 + t318) * qJD(4)) * t298) * t405 + ((-t435 * t477 + t476) * t298 + (t339 + (t298 * t478 + t318) * qJD(4)) * t299) * t408 + t340 + ((t448 * t295 + t447 * t296) * t305 + (t342 * t295 + t517 * t296) * qJD(4)) * t507 + ((t107 * t305 - t48) * t299 + (t106 * t305 + t49) * t298) * t506 + (t305 * t54 + ((-t522 * t298 - t321 * t299 + t305 * t75) * t299 + (t523 * t298 + t322 * t299 + t305 * t74) * t298) * t536) * t509 + (t305 * t55 + ((-t321 * t298 + t522 * t299 + t305 * t73) * t299 + (t322 * t298 - t523 * t299 + t305 * t72) * t298) * t536) * t508 + (t34 + t359) * t411 + (t358 + t35) * t410 + (-(-t60 * t419 + (t378 * t296 + t58 * (-t298 ^ 2 - t299 ^ 2) * t295) * qJD(4)) * pkin(4) + t58 * t430 + (t31 * t403 + t59 * t382 + t12 * t530 + t58 * t103 + (t58 * t125 + t60 * t403) * t305) * t299 + (t30 * t403 + t60 * t382 + t12 * t125 + t58 * t104 + (t59 * t213 + t543 * t58) * t305) * t298 + t518) * m(6) + (0.2e1 * t105 * ((t118 + t145) * t299 + (-t170 * t305 + t119) * t298) + t377 * t206 + ((-t305 * t70 - t53) * t299 + (-t52 + t493) * t298) * t228 - (t194 * t69 - t195 * t70) * t305 - (t105 * (-t194 * t298 - t195 * t299) + t377 * t379) * qJD(4)) * m(5); t340 + (t58 * (-t155 * t466 + t430) + t378 * t187 + ((-t305 * t60 - t31) * t299 + (-t30 + t494) * t298) * t213 + t518) * m(6);];
tauc = t1(:);
