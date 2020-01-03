% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:55
% EndTime: 2019-12-31 19:51:10
% DurationCPUTime: 11.28s
% Computational Cost: add. (13097->447), mult. (10421->539), div. (0->0), fcn. (7957->8), ass. (0->270)
t559 = Icges(6,4) + Icges(5,5);
t558 = Icges(5,6) - Icges(6,6);
t274 = pkin(8) + qJ(4);
t270 = cos(t274);
t269 = sin(t274);
t248 = Icges(6,5) * t269;
t346 = Icges(6,3) * t270 - t248;
t445 = Icges(5,4) * t269;
t553 = Icges(5,2) * t270 + t346 + t445;
t444 = Icges(6,5) * t270;
t199 = Icges(6,1) * t269 - t444;
t249 = Icges(5,4) * t270;
t557 = Icges(5,1) * t269 + t199 + t249;
t276 = qJ(1) + qJ(2);
t271 = sin(t276);
t194 = Icges(5,5) * t270 - Icges(5,6) * t269;
t272 = cos(t276);
t321 = t194 * t272;
t134 = Icges(5,3) * t271 + t321;
t196 = Icges(6,4) * t270 + Icges(6,6) * t269;
t322 = t196 * t272;
t136 = Icges(6,2) * t271 + t322;
t556 = t134 + t136;
t348 = Icges(6,1) * t270 + t248;
t139 = -Icges(6,4) * t272 + t271 * t348;
t430 = t269 * t271;
t228 = Icges(5,4) * t430;
t428 = t270 * t271;
t141 = Icges(5,1) * t428 - Icges(5,5) * t272 - t228;
t547 = t139 + t141;
t324 = t348 * t272;
t140 = Icges(6,4) * t271 + t324;
t202 = Icges(5,1) * t270 - t445;
t325 = t202 * t272;
t142 = Icges(5,5) * t271 + t325;
t555 = t140 + t142;
t192 = Icges(6,3) * t269 + t444;
t347 = -Icges(5,2) * t269 + t249;
t554 = t192 - t347;
t542 = t559 * t269 + t558 * t270;
t552 = t202 + t348;
t537 = -t553 * t269 + t557 * t270;
t427 = t270 * t272;
t227 = Icges(6,5) * t427;
t429 = t269 * t272;
t132 = Icges(6,6) * t271 + Icges(6,3) * t429 + t227;
t551 = t132 * t429 + t556 * t271 + t555 * t427;
t135 = -Icges(6,2) * t272 + t196 * t271;
t122 = t271 * t135;
t131 = -Icges(6,6) * t272 + t192 * t271;
t133 = Icges(5,5) * t428 - Icges(5,6) * t430 - Icges(5,3) * t272;
t550 = -t131 * t429 - t271 * t133 - t547 * t427 - t122;
t137 = Icges(5,4) * t428 - Icges(5,2) * t430 - Icges(5,6) * t272;
t549 = -t131 + t137;
t323 = t347 * t272;
t138 = Icges(5,6) * t271 + t323;
t548 = t132 - t138;
t545 = t554 * qJD(4);
t544 = t552 * qJD(4);
t525 = -t194 - t196;
t275 = qJD(1) + qJD(2);
t539 = t553 * qJD(4) - t558 * t275;
t538 = -qJD(4) * t557 + t559 * t275;
t482 = t542 * t272;
t483 = t542 * t271;
t439 = t137 * t269;
t341 = -t141 * t270 + t439;
t440 = t135 * t272;
t344 = t131 * t269 + t139 * t270;
t486 = t271 * t344;
t55 = -t440 + t486;
t536 = -t133 * t272 - t271 * t341 + t55;
t501 = -t137 * t429 - t550;
t500 = -t138 * t429 + t551;
t535 = t537 * t271 - t482;
t534 = t537 * t272 + t483;
t356 = -t132 * t430 + t136 * t272 - t140 * t428;
t112 = t142 * t428;
t365 = t134 * t272 - t112;
t58 = -t138 * t430 - t365;
t533 = -t356 + t58;
t508 = rSges(6,1) + pkin(4);
t532 = t508 * t269;
t426 = t271 * t275;
t530 = t539 * t272 - t426 * t554;
t424 = t272 * t275;
t529 = t192 * t424 + t539 * t271 - t275 * t323;
t499 = rSges(6,3) + qJ(5);
t528 = t538 * t272 - t426 * t552;
t527 = (-t324 - t325) * t275 - t538 * t271;
t526 = -t269 * t555 + t548 * t270;
t497 = t547 * t269 + t549 * t270;
t207 = rSges(6,1) * t270 + rSges(6,3) * t269;
t524 = pkin(4) * t270 + qJ(5) * t269 + t207;
t523 = t542 * t275 + t544 * t270 + t545 * t269 + (-t269 * t557 - t270 * t553) * qJD(4);
t522 = (Icges(6,2) + Icges(5,3)) * t275 - t542 * qJD(4);
t438 = t138 * t269;
t521 = -t132 * t269 - t270 * t555 + t438;
t520 = t341 - t344;
t519 = t525 * qJD(4) + t537 * t275;
t518 = (-t272 * t553 + t555) * t271 - (-Icges(5,2) * t428 - t346 * t271 - t228 + t547) * t272;
t517 = t534 * t275;
t516 = t549 * t272 + (-Icges(6,1) * t429 + t199 * t272 + t227 + t548) * t271;
t515 = (t500 * t271 - t501 * t272) * qJD(4);
t514 = (t533 * t271 - t536 * t272) * qJD(4);
t513 = t557 - t554;
t512 = -t553 + t552;
t509 = t535 * t275;
t507 = t509 + t514;
t506 = t515 + t517;
t505 = t520 * qJD(4) + t527 * t269 + t529 * t270;
t504 = -t521 * qJD(4) + t528 * t269 - t530 * t270;
t503 = -t519 * t271 + t523 * t272;
t502 = t523 * t271 + t519 * t272;
t405 = -t499 * t270 + t532;
t366 = t272 * t405;
t265 = t272 * rSges(6,2);
t414 = t524 * t271 - t265;
t496 = t414 * t275;
t495 = t271 * rSges(6,2) + pkin(4) * t427;
t494 = (-t133 - t135) * t275 + t527 * t270 - t529 * t269 + t497 * qJD(4);
t493 = t526 * qJD(4) + t530 * t269 + t528 * t270 + t556 * t275;
t492 = t440 + t551;
t395 = qJD(4) * t271;
t377 = t270 * t395;
t491 = t269 * t424 + t377;
t280 = sin(qJ(1));
t454 = pkin(1) * qJD(1);
t385 = t280 * t454;
t210 = rSges(3,1) * t271 + rSges(3,2) * t272;
t431 = t210 * t275;
t153 = -t385 - t431;
t490 = (t321 + t322 + t520) * t275 + t522 * t271;
t489 = t522 * t272 + t521 * t275 + t525 * t426;
t488 = 0.2e1 * qJD(4);
t463 = t271 / 0.2e1;
t392 = qJD(5) * t270;
t413 = rSges(6,1) * t427 + t499 * t429 + t495;
t50 = -t392 + (t414 * t271 + t413 * t272) * qJD(4);
t487 = qJD(4) * t50;
t455 = rSges(4,2) * sin(pkin(8));
t387 = t272 * t455;
t252 = t271 * qJ(3);
t477 = t272 * pkin(2) + t252;
t278 = cos(pkin(8));
t456 = rSges(4,1) * t278;
t478 = -t271 * rSges(4,3) - t272 * t456;
t316 = -t387 + t477 - t478;
t485 = t275 * t316;
t334 = rSges(5,1) * t427 + t271 * rSges(5,3);
t147 = -rSges(5,2) * t429 + t334;
t267 = pkin(3) * t278 + pkin(2);
t230 = t272 * t267;
t279 = -pkin(7) - qJ(3);
t425 = t271 * t279;
t363 = t230 - t425;
t484 = t147 + t363;
t253 = t272 * qJ(3);
t209 = pkin(2) * t271 - t253;
t247 = t272 * t279;
t399 = -t271 * t267 - t247;
t126 = t209 + t399;
t180 = t275 * t209;
t481 = -t275 * t126 + t180;
t394 = qJD(4) * t272;
t376 = t270 * t394;
t480 = rSges(6,2) * t424 + t499 * t376;
t379 = t269 * t395;
t479 = t508 * t379;
t476 = -t518 * t269 + t516 * t270;
t475 = (-t513 * t269 + t512 * t270) * t275;
t393 = qJD(5) * t269;
t223 = t272 * t393;
t250 = qJD(3) * t271;
t400 = t223 + t250;
t319 = -qJD(4) * t366 + t400;
t474 = -t319 + t400 + t480 + t481 + t496;
t473 = t525 * t275;
t389 = t271 * t456;
t245 = t271 * t455;
t397 = t272 * rSges(4,3) + t245;
t150 = t389 - t397;
t403 = rSges(4,3) * t424 + t275 * t245;
t472 = t275 * t150 + t180 + t403;
t251 = qJD(3) * t272;
t375 = t271 * t393;
t465 = t275 * (t413 + t363) - t251 + t375 - t405 * t395;
t388 = rSges(5,1) * t428;
t145 = -rSges(5,2) * t430 - t272 * rSges(5,3) + t388;
t125 = t275 * t145;
t384 = t269 * t426;
t326 = rSges(5,3) * t424 + (-t376 + t384) * rSges(5,2);
t378 = t269 * t394;
t464 = -rSges(5,1) * t378 + t125 + t250 + t326 + t481;
t462 = -t272 / 0.2e1;
t460 = t275 / 0.2e1;
t458 = pkin(1) * t280;
t281 = cos(qJ(1));
t273 = t281 * pkin(1);
t457 = pkin(2) - t267;
t205 = rSges(5,1) * t269 + rSges(5,2) * t270;
t355 = -t205 * t394 + t250;
t315 = t355 - t385;
t420 = t126 - t209;
t53 = (-t145 + t420) * t275 + t315;
t453 = t275 * t53;
t318 = -t270 * t426 - t378;
t449 = t508 * t318 - t499 * t384 + t223 + t480;
t448 = rSges(6,3) * t377 + t375 + t491 * qJ(5) - t479 + (t207 * t272 + t495) * t275;
t149 = t275 * t477 - t251;
t222 = t275 * t425;
t423 = t222 - (-t272 * t457 - t252) * t275 - t149;
t412 = -qJD(4) * t524 + t392;
t411 = t405 * t271;
t401 = t222 + t251;
t238 = qJ(3) * t424;
t398 = t238 + t250;
t396 = qJD(3) * t275;
t282 = qJD(1) ^ 2;
t391 = t282 * t458;
t390 = t282 * t273;
t386 = t281 * t454;
t381 = -rSges(5,1) * t379 - t491 * rSges(5,2);
t380 = t265 + t399;
t177 = t205 * t395;
t372 = -pkin(2) - t456;
t371 = -t395 / 0.2e1;
t368 = t394 / 0.2e1;
t212 = t272 * rSges(3,1) - rSges(3,2) * t271;
t364 = -t133 + t438;
t359 = t272 * t396 - t390;
t179 = rSges(3,1) * t424 - rSges(3,2) * t426;
t357 = -t251 + t386;
t354 = t392 + t412;
t352 = rSges(5,1) * t270 - rSges(5,2) * t269;
t54 = t275 * t484 - t177 + t357;
t351 = -t271 * t54 - t272 * t53;
t337 = t275 * (-pkin(2) * t426 + t398) + t271 * t396 - t391;
t336 = t222 - t357;
t333 = t230 + t413;
t332 = t275 * (-t238 + (t271 * t457 - t247) * t275) + t337;
t77 = (t145 * t271 + t147 * t272) * qJD(4);
t317 = -t145 + t399;
t312 = t271 * t372 + t253 + t397;
t302 = -t267 - t524;
t106 = (-t150 - t209) * t275 + t250 - t385;
t107 = t357 + t485;
t288 = (t106 * t372 * t272 + (t106 * (-rSges(4,3) - qJ(3)) + t107 * t372) * t271) * t275;
t285 = (t53 * (-t334 - t230) + t54 * (-t388 + t399)) * t275;
t284 = (((t58 - t112 + (t134 + t439) * t272 + t550) * t272 + (t55 - t486 + t492) * t271) * qJD(4) + t517) * t368 + (t537 * qJD(4) + t544 * t269 - t545 * t270) * t275 + (t503 + t504) * t395 / 0.2e1 + (((t272 * t364 - t492 + t500) * t272 + (t271 * t364 - t122 + t356 + t365 + t501) * t271) * qJD(4) + t507 - t509) * t371 - (t502 - t505 + t506) * t394 / 0.2e1 + ((t497 + t535) * t271 + (-t526 + t534) * t272) * qJD(4) * t460;
t14 = (t223 + t449) * t275 + (t271 * t354 - t275 * t366) * qJD(4) + t332;
t15 = t354 * t394 + ((qJD(4) * t405 - t393) * t271 + t423 - t448) * t275 + t359;
t45 = -t385 + (-t414 + t420) * t275 + t319;
t46 = t386 + t465;
t283 = (-t46 * qJD(4) * t532 + (-t46 * t279 + t302 * t45) * t275) * t272 + (-t14 * t279 + (-t45 * qJD(5) - t15 * t499) * t269 + (-qJD(4) * t45 * t499 - t15 * t508) * t270 + (-t45 * rSges(6,2) + t302 * t46) * t275) * t271;
t215 = t275 * t387;
t188 = t352 * qJD(4);
t173 = t205 * t272;
t169 = t205 * t271;
t154 = t212 * t275 + t386;
t130 = -t179 * t275 - t390;
t129 = -t275 * t431 - t391;
t105 = t275 * t334 + t381;
t103 = rSges(5,1) * t318 + t326;
t66 = (t275 * t478 - t149 + t215) * t275 + t359;
t65 = t275 * (-t275 * t389 + t403) + t337;
t38 = -t188 * t394 + (-t105 + t177 + t423) * t275 + t359;
t37 = t103 * t275 + (-t188 * t271 - t205 * t424) * qJD(4) + t332;
t5 = (t393 + (t449 + t496) * t272 + (-t413 * t275 + t448) * t271) * qJD(4);
t1 = [t284 + m(3) * (t130 * (-t210 - t458) + t129 * (t212 + t273) + (-t179 - t386 + t154) * t153) + (t15 * (t380 - t458) + t45 * (t336 + t479) + t14 * (t273 + t333) + t283 + (t45 + t474) * t46) * m(6) + (t38 * (t317 - t458) + t53 * (t336 - t381) + t37 * (t273 + t484) + t285 + (-t385 - t315 + t53 + t464) * t54) * m(5) + (t66 * (t312 - t458) + t106 * (t215 - t357) + t65 * (t273 + t316) + t288 + (t238 + t106 + t472) * t107) * m(4); t284 + (t14 * t333 + t15 * t380 + t283 + t474 * t46 + (t401 + t465 + t479) * t45) * m(6) + (t38 * t317 + t285 + (-t355 + t464) * t54 + (-t177 - t251 - t381 + t401) * t53 + (t37 + t453) * t484) * m(5) + (t66 * t312 + t65 * t316 + t288 + (-t250 + t398 + t472) * t107 + (t215 + t485) * t106) * m(4) + (t129 * t212 - t130 * t210 - t153 * t179 - t154 * t431 - (-t153 * t212 - t154 * t210) * t275) * m(3); 0.2e1 * (t14 * t462 + t15 * t463) * m(6) + 0.2e1 * (t37 * t462 + t38 * t463) * m(5) + 0.2e1 * (t65 * t462 + t463 * t66) * m(4); -((t512 * t269 + t513 * t270) * t275 + (t516 * t269 + t518 * t270) * qJD(4)) * t275 / 0.2e1 + ((-t275 * t526 + t505) * t272 + (t497 * t275 + t504) * t271) * t460 + ((-t395 * t482 - t473) * t271 + ((t271 * t483 + t476) * qJD(4) + t475) * t272) * t371 + ((-t394 * t483 + t473) * t272 + ((t272 * t482 + t476) * qJD(4) + t475) * t271) * t368 + ((-t15 * t405 + t45 * t412 + t5 * t413 + t50 * t449 + (-t405 * t46 + t414 * t50) * t275) * t272 + (-t14 * t405 + t46 * t412 + t5 * t414 + t50 * t448 + (t405 * t45 - t413 * t50) * t275) * t271 - (t50 * t269 + (t271 * t46 + t272 * t45) * t270) * qJD(5) - (-t366 * t46 + t411 * t45) * t275 - ((-t366 * t50 - t45 * t524) * t272 + (-t411 * t50 - t46 * t524) * t271) * qJD(4)) * m(6) + (0.2e1 * t77 * ((t103 + t125) * t272 + (-t147 * t275 + t105) * t271) + t351 * t188 + ((-t275 * t54 - t38) * t272 + (-t37 + t453) * t271) * t205 - (t169 * t53 - t173 * t54) * t275 - (t77 * (-t169 * t271 - t173 * t272) + t351 * t352) * qJD(4)) * m(5) + (t503 * t275 + ((t494 * t272 + t500 * t275) * t272 + (t489 * t271 + t501 * t275 + (-t490 + t493) * t272) * t271) * t488) * t463 + (t502 * t275 + ((t490 * t272 + t533 * t275) * t272 + (t493 * t271 + t536 * t275 + (-t489 + t494) * t272) * t271) * t488) * t462 + (t507 + t514) * t426 / 0.2e1 + (t506 + t515) * t424 / 0.2e1; (-t270 * t5 + 0.2e1 * (t487 / 0.2e1 + t14 * t463 + t15 * t272 / 0.2e1 - (t271 ^ 2 + t272 ^ 2) * t487 / 0.2e1) * t269) * m(6);];
tauc = t1(:);
