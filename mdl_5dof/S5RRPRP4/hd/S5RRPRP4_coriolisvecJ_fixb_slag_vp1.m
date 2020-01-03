% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:37
% EndTime: 2019-12-31 19:52:52
% DurationCPUTime: 12.07s
% Computational Cost: add. (9168->472), mult. (10011->566), div. (0->0), fcn. (7603->6), ass. (0->268)
t276 = qJ(1) + qJ(2);
t270 = sin(t276);
t271 = cos(t276);
t279 = cos(qJ(4));
t423 = t271 * t279;
t277 = sin(qJ(4));
t424 = t271 * t277;
t128 = Icges(6,5) * t424 - Icges(6,6) * t270 - Icges(6,3) * t423;
t447 = Icges(5,4) * t277;
t348 = Icges(5,2) * t279 + t447;
t134 = -Icges(5,6) * t270 + t271 * t348;
t559 = t128 - t134;
t242 = Icges(6,5) * t423;
t136 = Icges(6,1) * t424 - Icges(6,4) * t270 - t242;
t446 = Icges(5,4) * t279;
t350 = Icges(5,1) * t277 + t446;
t138 = -Icges(5,5) * t270 + t271 * t350;
t564 = t136 + t138;
t346 = Icges(5,5) * t277 + Icges(5,6) * t279;
t129 = Icges(5,3) * t271 + t270 * t346;
t347 = Icges(6,4) * t277 - Icges(6,6) * t279;
t315 = t347 * t270;
t131 = Icges(6,2) * t271 + t315;
t563 = t129 + t131;
t221 = -Icges(5,2) * t277 + t446;
t272 = Icges(6,5) * t279;
t492 = Icges(6,3) * t277 + t272;
t562 = -t221 + t492;
t443 = Icges(6,5) * t277;
t223 = Icges(6,1) * t279 + t443;
t225 = Icges(5,1) * t279 - t447;
t561 = t223 + t225;
t537 = t564 * t277 - t559 * t279;
t344 = -Icges(6,3) * t279 + t443;
t560 = t344 - t348;
t130 = -Icges(5,3) * t270 + t271 * t346;
t132 = Icges(6,4) * t424 - Icges(6,2) * t270 - Icges(6,6) * t423;
t539 = t130 + t132;
t335 = t221 * t279 + t277 * t225;
t553 = t223 * t277 - t279 * t492 + t335;
t427 = t270 * t277;
t241 = Icges(6,5) * t427;
t426 = t270 * t279;
t127 = Icges(6,6) * t271 - Icges(6,3) * t426 + t241;
t512 = -t127 * t423 - t563 * t270;
t557 = t560 * qJD(4);
t349 = Icges(6,1) * t277 - t272;
t556 = (-t349 - t350) * qJD(4);
t217 = Icges(5,5) * t279 - Icges(5,6) * t277;
t219 = Icges(6,4) * t279 + Icges(6,6) * t277;
t555 = -t219 - t217;
t275 = qJD(1) + qJD(2);
t554 = t560 * t275;
t552 = (Icges(5,6) - Icges(6,6)) * t275 + t562 * qJD(4);
t551 = (-Icges(6,4) - Icges(5,5)) * t275 + t561 * qJD(4);
t550 = t537 * t271;
t533 = rSges(6,1) + pkin(4);
t317 = t349 * t270;
t135 = Icges(6,4) * t271 + t317;
t421 = t271 * t131 + t135 * t427;
t51 = -t127 * t426 + t421;
t133 = Icges(5,6) * t271 + t270 * t348;
t243 = Icges(5,4) * t426;
t137 = Icges(5,1) * t427 + Icges(5,5) * t271 + t243;
t53 = t271 * t129 + t133 * t426 + t137 * t427;
t549 = t51 + t53;
t519 = -t539 * t271 + t426 * t559 - t427 * t564;
t340 = t133 * t279 + t137 * t277;
t548 = -t135 * t424 - t271 * t340 - t512;
t517 = -t270 * t539 + t550;
t516 = rSges(6,3) + qJ(5);
t547 = t516 * t279;
t429 = t219 * t271;
t431 = t217 * t271;
t546 = t270 * t553 + t429 + t431;
t150 = t270 * t217;
t152 = t270 * t219;
t545 = t223 * t424 + t271 * t335 - t423 * t492 - t150 - t152;
t509 = t135 * t277 + t340;
t544 = t270 * t554 - t271 * t552;
t543 = t270 * t552 + t271 * t554;
t318 = t350 * t275;
t542 = t270 * t318 - t271 * t551 + t275 * t317;
t425 = t271 * t275;
t541 = t270 * t551 + t271 * t318 + t349 * t425;
t540 = -t277 * t559 - t279 * t564;
t515 = (-t135 - t137) * t279 + (-t127 + t133) * t277;
t538 = t553 * t275 + (-t346 - t347) * qJD(4);
t435 = t127 * t279;
t536 = -t435 + t509;
t535 = t557 * t279 + t556 * t277 + t555 * t275 + (t562 * t277 + t561 * t279) * qJD(4);
t534 = (Icges(6,2) + Icges(5,3)) * t275 + t555 * qJD(4);
t495 = -t270 * rSges(6,2) - rSges(6,3) * t423;
t410 = -qJ(5) * t423 + t424 * t533 + t495;
t532 = t410 * t275;
t531 = t546 * t275;
t399 = t516 * t277 + t279 * t533;
t530 = (t270 * t517 + t271 * t548) * qJD(4);
t529 = (t270 * t519 + t271 * t549) * qJD(4);
t528 = t545 * t275;
t527 = t533 * t277;
t526 = -t536 * qJD(4) + t540 * t275 + t543 * t277 + t541 * t279;
t463 = t270 / 0.2e1;
t462 = -t271 / 0.2e1;
t525 = t529 + t531;
t524 = -t528 + t530;
t522 = t537 * qJD(4) + t544 * t277 + t542 * t279;
t521 = t538 * t270 - t535 * t271;
t520 = t535 * t270 + t538 * t271;
t404 = t492 - t349;
t405 = -t344 - t223;
t514 = (t277 * t404 - t279 * t405) * t275;
t402 = t221 + t350;
t403 = -t348 + t225;
t513 = (t277 * t402 - t279 * t403) * t275;
t263 = t271 * rSges(5,3);
t140 = rSges(5,1) * t427 + rSges(5,2) * t426 + t263;
t268 = t271 * pkin(2);
t183 = t270 * qJ(3) + t268;
t267 = t271 * pkin(7);
t493 = t267 + t183;
t325 = t140 + t493;
t400 = -t527 + t547;
t494 = t271 * rSges(6,2) + t533 * t427;
t497 = t346 * t275;
t511 = -t534 * t270 + t271 * t497 + t536 * t275 + t347 * t425;
t510 = t270 * t497 + (t315 - t537) * t275 + t534 * t271;
t392 = qJD(4) * t279;
t507 = t270 * t392 + t275 * t424;
t278 = sin(qJ(1));
t452 = pkin(1) * qJD(1);
t385 = t278 * t452;
t182 = rSges(3,1) * t270 + rSges(3,2) * t271;
t433 = t182 * t275;
t143 = -t385 - t433;
t506 = t515 * qJD(4) + t563 * t275 - t541 * t277 + t543 * t279;
t505 = t540 * qJD(4) + t539 * t275 + t542 * t277 - t544 * t279;
t503 = 0.2e1 * qJD(4);
t269 = qJD(5) * t277;
t411 = t516 * t426 - t494;
t50 = t269 + (t411 * t270 - t410 * t271) * qJD(4);
t502 = qJD(4) * t50;
t328 = -rSges(4,2) * t271 + t270 * rSges(4,3) + t183;
t500 = t275 * t328;
t456 = pkin(7) * t270;
t499 = t275 * t456;
t395 = qJD(4) * t270;
t496 = t399 * t395;
t352 = rSges(5,1) * t277 + rSges(5,2) * t279;
t142 = -t270 * rSges(5,3) + t271 * t352;
t125 = t275 * t142;
t393 = qJD(4) * t277;
t377 = t270 * t393;
t422 = t275 * t279;
t380 = t271 * rSges(5,2) * t422 + t507 * rSges(5,1);
t213 = qJ(3) * t425;
t249 = qJD(3) * t270;
t406 = t213 + t249;
t491 = -rSges(5,2) * t377 - t125 + t380 + t406;
t391 = qJD(5) * t279;
t372 = t271 * t391;
t374 = t271 * t392;
t375 = t271 * t393;
t490 = -t372 + t516 * (t270 * t422 + t375) + t533 * t374;
t489 = t516 * t377 + t533 * t507;
t181 = t270 * rSges(4,2) + t271 * rSges(4,3);
t428 = t270 * t275;
t401 = rSges(4,2) * t428 + rSges(4,3) * t425;
t488 = -t275 * t181 + t401;
t487 = t406 + t489 - t532;
t274 = t275 ^ 2;
t281 = qJD(1) ^ 2;
t458 = pkin(1) * t278;
t388 = t281 * t458;
t396 = qJD(3) * t275;
t333 = t275 * (-pkin(2) * t428 + t406) + t270 * t396 - t388;
t295 = -t274 * t456 + t333;
t407 = t400 * qJD(4) + t269;
t354 = t269 + t407;
t373 = t270 * t391;
t454 = -t495 * t275 - (-qJ(5) * t425 - qJD(5) * t270) * t279 - t489;
t14 = (-t373 - t454) * t275 + (-t271 * t354 + t399 * t428) * qJD(4) + t295;
t250 = qJD(3) * t271;
t126 = t183 * t275 - t250;
t280 = cos(qJ(1));
t273 = t280 * pkin(1);
t387 = t281 * t273;
t358 = t271 * t396 - t387;
t305 = -t267 * t274 + t358;
t455 = t494 * t275 - t490;
t15 = t354 * t395 + (-t126 + (qJD(4) * t399 - t391) * t271 - t455) * t275 + t305;
t486 = t14 * t462 + t15 * t463;
t236 = rSges(5,1) * t279 - rSges(5,2) * t277;
t394 = qJD(4) * t271;
t177 = t236 * t394;
t483 = t275 * t325 - t177;
t482 = t275 * (t493 - t411) - t250 + t372 - t399 * t394;
t412 = t221 * t271 + t138;
t416 = -t225 * t271 + t134;
t469 = t277 * t416 - t279 * t412;
t413 = -Icges(5,2) * t427 + t137 + t243;
t417 = -t225 * t270 + t133;
t468 = t277 * t417 - t279 * t413;
t414 = -Icges(6,3) * t424 + t136 - t242;
t418 = t223 * t271 + t128;
t467 = -t277 * t418 - t279 * t414;
t415 = -t270 * t492 + t135;
t419 = Icges(6,1) * t426 + t127 + t241;
t466 = -t277 * t419 - t279 * t415;
t465 = t270 ^ 2;
t464 = -pkin(2) - pkin(7);
t460 = -t275 / 0.2e1;
t430 = t347 * t275;
t409 = t399 * t270;
t408 = t399 * t271;
t252 = t271 * qJ(3);
t180 = pkin(2) * t270 - t252;
t332 = -t180 + t181;
t169 = t275 * t180;
t397 = t249 - t169;
t390 = -rSges(6,2) + t464;
t389 = -rSges(5,3) + t464;
t386 = t280 * t452;
t369 = -t395 / 0.2e1;
t367 = -t394 / 0.2e1;
t364 = -t180 - t456;
t185 = t271 * rSges(3,1) - rSges(3,2) * t270;
t362 = t132 - t435;
t147 = rSges(3,1) * t425 - rSges(3,2) * t428;
t356 = t249 - t385;
t355 = -t250 + t386;
t174 = t236 * t395;
t329 = t174 + t356;
t61 = (t142 + t364) * t275 + t329;
t62 = t355 + t483;
t351 = t270 * t61 - t271 * t62;
t330 = -t169 + t356;
t327 = -rSges(5,1) * t374 + rSges(5,2) * t375;
t326 = t493 + t494;
t319 = t250 - t327;
t312 = t249 - t373 + t496;
t71 = (-t140 * t270 - t142 * t271) * qJD(4);
t293 = t252 + t410;
t292 = t250 + t490;
t286 = t270 * t464 + t142 + t252;
t106 = t275 * t332 + t356;
t107 = t355 + t500;
t285 = (-t106 * t268 + (t106 * (-rSges(4,3) - qJ(3)) - t107 * pkin(2)) * t270) * t275;
t284 = (t61 * t389 * t271 + (t61 * (-qJ(3) - t352) + t62 * t389) * t270) * t275;
t283 = (((t421 + t53 + t517 - t550) * t271 + ((t362 + t130 - t509) * t271 - t512 - t548 + t519) * t270) * qJD(4) + t531) * t369 + (-qJD(4) * t553 - t277 * t557 + t279 * t556) * t275 + ((t130 * t465 + (t270 * t362 + t421 - t51) * t270 + ((t509 + t539) * t271 + t512 + t519) * t271) * qJD(4) + t524 + t528) * t367 + (t521 + t522 + t525) * t395 / 0.2e1 + (t520 + t526) * t394 / 0.2e1 + (t545 * t271 + (-t515 + t546) * t270) * qJD(4) * t460;
t47 = -t385 + (t364 + t410) * t275 + t312;
t48 = t386 + t482;
t282 = (t15 * t464 + (-t48 * qJD(5) - t14 * t516) * t279) * t270 + ((t390 * t47 - t48 * t547) * t271 + (t47 * (-qJ(3) - t527) + t48 * t390) * t270) * t275;
t229 = rSges(4,2) * t425;
t198 = t352 * qJD(4);
t167 = t236 * t271;
t163 = t236 * t270;
t144 = t185 * t275 + t386;
t117 = -t147 * t275 - t387;
t116 = -t275 * t433 - t388;
t97 = (-rSges(5,2) * t393 - rSges(5,3) * t275) * t270 + t380;
t95 = (t270 * t352 + t263) * t275 + t327;
t64 = (-rSges(4,3) * t428 - t126 + t229) * t275 + t358;
t63 = t275 * t401 + t333;
t42 = -t198 * t395 + (-t126 - t95 + t177) * t275 + t305;
t41 = t275 * t97 + (t198 * t271 + t236 * t428) * qJD(4) + t295;
t5 = (t391 + (t411 * t275 + t455) * t271 + (t454 + t532) * t270) * qJD(4);
t1 = [m(3) * (t117 * (-t182 - t458) + t116 * (t185 + t273) + (-t147 - t386 + t144) * t143) + t283 + (t15 * (t293 - t458) + t47 * (t292 - t386) + t14 * (t273 + t326) + t282 + (t47 - (-pkin(7) * t275 - t391) * t270 - t330 - t385 + t487 - t496) * t48) * m(6) + (t42 * (t286 - t458) + t61 * (t319 - t386) + t41 * (t273 + t325) + t284 + (pkin(7) * t428 + t169 - t329 - t385 + t491 + t61) * t62) * m(5) + (t64 * (t332 - t458) + t106 * (t229 - t355) + t63 * (t273 + t328) + t285 + (t213 + t356 + t106 - t330 + t488) * t107) * m(4); t283 + (t14 * t326 + t15 * t293 + t282 + (t169 - t312 + t487 + t499) * t48 + (t292 + t482) * t47) * m(6) + (t42 * t286 + t41 * t325 + t284 + (-t174 - t397 + t491 + t499) * t62 + (-t250 + t319 + t483) * t61) * m(5) + (t63 * t328 + t64 * t332 + t285 + (-t397 + t406 + t488) * t107 + (t229 + t500) * t106) * m(4) + (-(-t143 * t185 - t144 * t182) * t275 + t116 * t185 - t117 * t182 - t143 * t147 - t144 * t433) * m(3); 0.2e1 * t486 * m(6) + 0.2e1 * (t41 * t462 + t42 * t463) * m(5) + 0.2e1 * (t462 * t63 + t463 * t64) * m(4); (((-t402 + t404) * t279 + (-t403 + t405) * t277) * t275 + (((-t417 + t419) * t271 + (t416 - t418) * t270) * t279 + ((-t413 - t415) * t271 + (t412 + t414) * t270) * t277) * qJD(4)) * t460 + (t526 * t271 + (t515 * t275 + t522) * t270) * t275 / 0.2e1 + ((-t395 * t429 - t430) * t270 + (-t514 + (t466 * t271 + (t152 - t467) * t270) * qJD(4)) * t271 + (-t395 * t431 - t497) * t270 + (t513 + (t468 * t271 + (t150 - t469) * t270) * qJD(4)) * t271) * t369 + ((t152 * t394 - t430) * t271 + (t514 + (t467 * t270 + (-t429 - t466) * t271) * qJD(4)) * t270 + (t150 * t394 - t497) * t271 + (-t513 + (t469 * t270 + (-t431 - t468) * t271) * qJD(4)) * t270) * t367 + (-(t279 * t50 + (t270 * t47 - t271 * t48) * t277) * qJD(5) - (t408 * t47 + t409 * t48) * t275 - ((-t400 * t48 - t408 * t50) * t271 + (t400 * t47 - t409 * t50) * t270) * qJD(4) + (-t14 * t399 - t48 * t407 - t5 * t410 + t50 * t455 + (t399 * t47 + t411 * t50) * t275) * t271 + (t15 * t399 + t47 * t407 + t5 * t411 + t50 * t454 + (t399 * t48 + t410 * t50) * t275) * t270) * m(6) + (0.2e1 * t71 * ((-t140 * t275 + t95) * t271 + (-t97 + t125) * t270) - t351 * t198 + ((t275 * t61 - t41) * t271 + (t275 * t62 + t42) * t270) * t236 - (t163 * t62 + t167 * t61) * t275 - (t71 * (-t163 * t270 - t167 * t271) - t351 * t352) * qJD(4)) * m(5) + (t521 * t275 + ((t506 * t271 + t517 * t275) * t271 + (t510 * t270 - t548 * t275 + (-t505 + t511) * t271) * t270) * t503) * t463 + (t520 * t275 + ((t511 * t271 + t519 * t275) * t271 + (t505 * t270 - t549 * t275 + (-t506 + t510) * t271) * t270) * t503) * t271 / 0.2e1 - (t525 + t529) * t428 / 0.2e1 + (t524 + t530) * t425 / 0.2e1; (t277 * t5 + 0.2e1 * (t502 / 0.2e1 - (t271 ^ 2 + t465) * t502 / 0.2e1 - t486) * t279) * m(6);];
tauc = t1(:);
