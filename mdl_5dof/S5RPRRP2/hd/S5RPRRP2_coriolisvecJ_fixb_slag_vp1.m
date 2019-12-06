% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:23
% EndTime: 2019-12-05 18:01:41
% DurationCPUTime: 11.82s
% Computational Cost: add. (12066->408), mult. (8983->496), div. (0->0), fcn. (6890->8), ass. (0->260)
t545 = Icges(5,3) + Icges(6,3);
t248 = qJ(1) + pkin(8);
t243 = qJ(3) + t248;
t237 = sin(t243);
t238 = cos(t243);
t252 = cos(qJ(4));
t391 = t238 * t252;
t250 = sin(qJ(4));
t392 = t238 * t250;
t128 = Icges(6,4) * t391 - Icges(6,2) * t392 + Icges(6,6) * t237;
t130 = Icges(5,4) * t391 - Icges(5,2) * t392 + Icges(5,6) * t237;
t533 = t128 + t130;
t204 = Icges(6,4) * t392;
t132 = Icges(6,1) * t391 + Icges(6,5) * t237 - t204;
t206 = Icges(5,4) * t392;
t134 = Icges(5,1) * t391 + Icges(5,5) * t237 - t206;
t531 = t132 + t134;
t544 = -Icges(5,5) - Icges(6,5);
t543 = Icges(5,6) + Icges(6,6);
t415 = Icges(6,4) * t250;
t215 = Icges(6,2) * t252 + t415;
t416 = Icges(5,4) * t250;
t217 = Icges(5,2) * t252 + t416;
t539 = t215 + t217;
t245 = Icges(5,4) * t252;
t469 = Icges(5,1) * t250 + t245;
t244 = Icges(6,4) * t252;
t470 = Icges(6,1) * t250 + t244;
t535 = t469 + t470;
t212 = Icges(6,5) * t252 - Icges(6,6) * t250;
t214 = Icges(5,5) * t252 - Icges(5,6) * t250;
t286 = t214 * t237;
t541 = -t212 * t237 + t545 * t238 - t286;
t315 = -Icges(6,2) * t250 + t244;
t127 = Icges(6,6) * t238 - t237 * t315;
t316 = -Icges(5,2) * t250 + t245;
t129 = Icges(5,6) * t238 - t237 * t316;
t540 = t127 + t129;
t395 = t237 * t250;
t203 = Icges(6,4) * t395;
t394 = t237 * t252;
t131 = -Icges(6,1) * t394 + Icges(6,5) * t238 + t203;
t205 = Icges(5,4) * t395;
t133 = -Icges(5,1) * t394 + Icges(5,5) * t238 + t205;
t532 = -t131 - t133;
t220 = Icges(6,1) * t252 - t415;
t222 = Icges(5,1) * t252 - t416;
t538 = t220 + t222;
t537 = t533 * t250 - t531 * t252;
t536 = t315 + t316;
t504 = -t237 * t545 + t544 * t391 + t543 * t392;
t211 = Icges(6,5) * t250 + Icges(6,6) * t252;
t213 = Icges(5,5) * t250 + Icges(5,6) * t252;
t530 = t211 + t213;
t526 = t539 * t250 - t535 * t252;
t247 = qJD(1) + qJD(3);
t529 = (t535 * qJD(4) + t544 * t247) * t250 + (t539 * qJD(4) - t543 * t247) * t252;
t528 = t536 * qJD(4);
t527 = t538 * qJD(4);
t479 = -t541 * t237 + t532 * t391;
t478 = t541 * t238 + t540 * t395;
t525 = t532 * t252;
t477 = t540 * t250;
t522 = t537 * t238;
t521 = (-t538 * t250 - t536 * t252) * t247;
t255 = qJD(1) ^ 2;
t520 = t394 * t532 + t478;
t488 = -t238 * t504 - t531 * t394 + t533 * t395;
t519 = -t392 * t540 - t479;
t507 = -t237 * t504 - t522;
t399 = t213 * t238;
t401 = t211 * t238;
t518 = t237 * t526 + t399 + t401;
t143 = t211 * t237;
t145 = t213 * t237;
t517 = -t238 * t526 + t143 + t145;
t515 = t477 + t525;
t514 = t530 * qJD(4) - t247 * t545;
t393 = t238 * t247;
t400 = t212 * t247;
t513 = (-t214 * t393 + t514 * t237 - t400 * t238 + t515 * t247) * t238;
t512 = t526 * t247 + (t212 + t214) * qJD(4);
t433 = pkin(4) * t252;
t239 = pkin(3) + t433;
t183 = t238 * t239;
t360 = rSges(6,1) * t391;
t302 = rSges(6,3) * t237 + t360;
t511 = -t302 - t183;
t510 = -t527 * t252 + t528 * t250 - t530 * t247 + (t535 * t250 + t539 * t252) * qJD(4);
t506 = t518 * t247;
t209 = rSges(6,2) * t392;
t249 = -qJ(5) - pkin(7);
t429 = pkin(7) + t249;
t434 = pkin(3) * t238;
t384 = t237 * t429 + t209 + t434 + t511;
t505 = t247 * t384;
t502 = -t400 * t237 + (-t286 + t537) * t247 - t514 * t238;
t501 = (t507 * t237 + t519 * t238) * qJD(4);
t500 = (t488 * t237 + t520 * t238) * qJD(4);
t472 = rSges(6,2) * t395 + t238 * rSges(6,3);
t499 = t429 * t238 - t472;
t498 = t517 * t247;
t495 = 0.2e1 * qJD(4);
t494 = t500 + t506;
t493 = t498 + t501;
t492 = -t515 * qJD(4) + t237 * t529 + t521 * t238;
t491 = -qJD(4) * t537 + t521 * t237 - t238 * t529;
t490 = t512 * t237 - t510 * t238;
t489 = t510 * t237 + t512 * t238;
t487 = t250 * t532 - t252 * t540;
t486 = t250 * t531 + t252 * t533;
t370 = t470 + t315;
t371 = t215 - t220;
t485 = (t250 * t370 + t252 * t371) * t247;
t368 = t469 + t316;
t369 = t217 - t222;
t484 = (t250 * t368 + t252 * t369) * t247;
t483 = t504 + t525;
t423 = rSges(6,1) * t252;
t334 = -t239 - t423;
t482 = (-pkin(3) - t334) * t237 + t499;
t364 = qJD(4) * t250;
t346 = t238 * t364;
t481 = -t247 * t394 - t346;
t363 = qJD(4) * t252;
t480 = t237 * t363 + t247 * t392;
t242 = cos(t248);
t253 = cos(qJ(1));
t437 = pkin(1) * t253;
t327 = pkin(2) * t242 + t437;
t463 = t327 * qJD(1);
t241 = sin(t248);
t435 = pkin(2) * t241;
t251 = sin(qJ(1));
t438 = pkin(1) * t251;
t328 = t435 + t438;
t297 = t328 * qJD(1);
t425 = rSges(4,1) * t238;
t164 = -t237 * rSges(4,2) + t425;
t390 = t247 * t164;
t118 = -t463 - t390;
t321 = rSges(4,1) * t237 + rSges(4,2) * t238;
t476 = t118 * t321;
t138 = t321 * t247;
t210 = rSges(5,2) * t392;
t361 = rSges(5,1) * t391;
t303 = -rSges(5,3) * t237 - t361;
t136 = -t210 - t303;
t122 = t247 * t136;
t224 = rSges(5,1) * t250 + rSges(5,2) * t252;
t366 = qJD(4) * t237;
t163 = t224 * t366;
t473 = -t122 + t163;
t372 = rSges(5,2) * t395 + t238 * rSges(5,3);
t236 = t241 * rSges(3,2);
t426 = rSges(3,1) * t242;
t325 = -t426 - t437;
t471 = t236 + t325;
t229 = qJD(5) * t237;
t223 = rSges(6,1) * t250 + rSges(6,2) * t252;
t335 = pkin(4) * t250 + t223;
t319 = qJD(4) * t335;
t468 = t238 * t319 - t229;
t345 = t238 * t363;
t389 = t247 * t249;
t396 = t237 * t247;
t466 = rSges(6,1) * t481 - rSges(6,2) * t345 - t238 * t389 - t239 * t396 + t229;
t348 = t237 * t364;
t230 = qJD(5) * t238;
t374 = pkin(4) * t348 + t230;
t465 = rSges(6,1) * t348 + rSges(6,2) * t480 + t237 * t389 + t374;
t422 = rSges(5,2) * t250;
t424 = rSges(5,1) * t252;
t320 = -t422 + t424;
t185 = t320 * qJD(4);
t464 = -t185 * t238 + t224 * t396;
t356 = t223 * t366 + t374;
t462 = t356 + t505;
t376 = -Icges(5,2) * t391 + t134 - t206;
t380 = t238 * t469 + t130;
t447 = t250 * t376 + t252 * t380;
t377 = Icges(5,2) * t394 + t133 + t205;
t381 = -t237 * t469 + t129;
t446 = -t250 * t377 - t252 * t381;
t378 = -Icges(6,2) * t391 + t132 - t204;
t382 = t238 * t470 + t128;
t445 = t250 * t378 + t252 * t382;
t379 = Icges(6,2) * t394 + t131 + t203;
t383 = -t237 * t470 + t127;
t444 = -t250 * t379 - t252 * t383;
t439 = -rSges(5,3) - pkin(7);
t436 = pkin(1) * t255;
t432 = pkin(7) * t237;
t431 = t238 * pkin(7);
t430 = pkin(3) - t239;
t198 = pkin(3) * t396;
t428 = t198 + (-pkin(4) * t364 - pkin(7) * t247) * t238 + t247 * t472 + t466;
t427 = -t465 + (-t238 * t430 + t302 - t432) * t247;
t165 = t432 + t434;
t60 = t163 + (-t136 - t165) * t247 - t463;
t421 = t237 * t60;
t398 = t214 * t247;
t156 = t224 * t237;
t397 = t224 * t238;
t240 = t251 * t436;
t367 = t255 * t435 + t240;
t365 = qJD(4) * t238;
t354 = rSges(5,1) * t481 - rSges(5,2) * t345;
t352 = rSges(5,1) * t348 + rSges(5,2) * t480;
t349 = t224 * t365;
t342 = -pkin(3) - t424;
t340 = -t366 / 0.2e1;
t338 = -t365 / 0.2e1;
t337 = t365 / 0.2e1;
t329 = t372 + t431;
t326 = -t238 * t249 + t472;
t322 = rSges(3,1) * t241 + rSges(3,2) * t242;
t226 = -rSges(6,2) * t250 + t423;
t306 = rSges(5,1) * t394 - t372;
t305 = -t183 + t209 - t360;
t120 = t247 * t306;
t160 = t247 * (-pkin(3) * t237 + t431);
t304 = t120 - t160 + t349;
t184 = t226 * qJD(4);
t299 = (t433 * qJD(4) + t184) * qJD(4);
t298 = t327 * t255;
t284 = t210 - t361 - t434;
t278 = -t247 ^ 2 * t165 - t298;
t161 = t247 * t165;
t275 = t161 + t463;
t267 = t238 * t136 + t237 * t306;
t266 = -t160 + t468 + (rSges(6,1) * t394 - t237 * t430 + t499) * t247;
t91 = t247 * t372 + t354;
t93 = t247 * t303 + t352;
t261 = (-t93 - t122) * t237 + (t91 + t120) * t238;
t260 = t237 * t482 - t384 * t238;
t259 = (((t478 + t507 + t522) * t238 + ((-t477 + t483) * t238 - t479 + t488 - t519) * t237) * qJD(4) + t506) * t340 + ((((t477 + t504) * t238 + t479 + t488) * t238 + (t237 * t483 + t478 - t520) * t237) * qJD(4) + t493 - t498) * t338 + (t489 + t492) * t337 + (t490 + t491 + t494) * t366 / 0.2e1 + (t528 * t252 + t527 * t250 - t526 * qJD(4) + (-t487 + t518) * t340 + (t486 + t517) * t337) * t247;
t258 = (t427 + t505) * t237 + (t247 * t482 + t428) * t238;
t140 = pkin(7) * t393 - t198;
t47 = t185 * t366 + (-t140 - t91 + t349) * t247 + t367;
t48 = qJD(4) * t464 + t247 * t93 + t278;
t59 = t297 + t304;
t257 = t60 * (t198 - t354) + (t48 * t342 + t47 * t439) * t237 + ((-t422 * t60 - t439 * t59) * t237 + (-t342 * t59 + t439 * t60) * t238) * t247 - t59 * t352;
t23 = t299 * t237 + (-t140 - t428 + t468) * t247 + t367;
t24 = -t299 * t238 + (t237 * t319 + t230 - t427) * t247 + t278;
t49 = t297 + t266;
t50 = -t463 + (-t165 + t384) * t247 + t356;
t256 = t50 * (pkin(4) * t346 - t466) + (t23 * (-rSges(6,3) + t249) + t24 * t334) * t237 + (-t50 * t472 - t49 * t511) * t247 - t49 * t465;
t197 = rSges(4,2) * t396;
t157 = t223 * t238;
t155 = t223 * t237;
t139 = -rSges(4,1) * t393 + t197;
t117 = t297 + t138;
t103 = t139 * t247 - t298;
t102 = t138 * t247 + t367;
t65 = qJD(4) * t267 + qJD(2);
t42 = qJD(4) * t260 + qJD(2);
t25 = t261 * qJD(4);
t5 = t258 * qJD(4);
t1 = [t259 + m(3) * ((t255 * t322 + t240) * t471 + (t253 * t436 + (-0.2e1 * t236 - t325 + t426 + t471) * t255) * (t322 + t438)) + (t23 * (t305 - t327) + t24 * (t326 - t328) + (t327 * t49 + t328 * t50) * qJD(1) + t256 - (t50 + t275 - t462) * t49) * m(6) + (t47 * (t284 - t327) + t48 * (-t328 + t329) + (t327 * t59 + t328 * t60) * qJD(1) + t257 - (t60 + t275 - t473) * t59) * m(5) + (t102 * (-t164 - t327) + t103 * (-t321 - t328) + t476 * t247 + t118 * t297 + (t247 * t425 - t118 - t197 - t390) * t117) * m(4); m(5) * t25 + m(6) * t5; t259 + (t23 * t305 + t24 * t326 + t256 - t50 * t266 + t49 * (-t161 + t462)) * m(6) + (t47 * t284 + t48 * t329 + t257 - t60 * t304 + t59 * (-t161 + t473)) * m(5) + (-(t117 * t164 + t476) * t247 - t102 * t164 - t103 * t321 - t117 * t139 + t118 * t138) * m(4); -(((t368 + t370) * t252 + (-t369 - t371) * t250) * t247 + (((t377 + t379) * t238 + (t376 + t378) * t237) * t252 + ((-t381 - t383) * t238 + (-t380 - t382) * t237) * t250) * qJD(4)) * t247 / 0.2e1 + ((t247 * t486 + t492) * t238 + (t247 * t487 + t491) * t237) * t247 / 0.2e1 + ((-t366 * t401 + t400) * t237 + (-t485 + (t444 * t238 + (t143 - t445) * t237) * qJD(4)) * t238 + (-t366 * t399 + t398) * t237 + (-t484 + (t446 * t238 + (t145 - t447) * t237) * qJD(4)) * t238) * t340 + ((t143 * t365 + t400) * t238 + (t485 + (t445 * t237 + (-t401 - t444) * t238) * qJD(4)) * t237 + (t145 * t365 + t398) * t238 + (t484 + (t447 * t237 + (-t399 - t446) * t238) * qJD(4)) * t237) * t338 + (-((t50 * t226 + (-pkin(4) * t395 - t155) * t42) * t237 + (-t49 * (-t226 - t433) + (-pkin(4) * t392 - t157) * t42) * t238) * qJD(4) + t5 * t260 + t42 * t258 + t50 * (t184 * t237 + t223 * t393) - t49 * (t223 * t396 + (-pkin(4) * t363 - t184) * t238) + (t23 * t237 - t24 * t238) * t335 + (t49 * t155 - t50 * t157) * t247) * m(6) + (-(-t156 * t59 + t397 * t60) * t247 - (t65 * (-t156 * t237 - t238 * t397) + (t238 * t59 + t421) * t320) * qJD(4) + t60 * t224 * t393 + t47 * t156 + t185 * t421 + t25 * t267 + t65 * t261 - t48 * t397 - t464 * t59) * m(5) + (t490 * t247 + (t507 * t393 + (t502 * t237 - t247 * t519 + t513) * t237) * t495) * t237 / 0.2e1 + (t489 * t247 + ((t488 * t247 + t513) * t238 + (t502 * t238 - t247 * t520) * t237) * t495) * t238 / 0.2e1 - (t494 + t500) * t396 / 0.2e1 + (t493 + t501) * t393 / 0.2e1; m(6) * (t23 * t238 + t237 * t24);];
tauc = t1(:);
