% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:09
% EndTime: 2019-12-05 18:58:23
% DurationCPUTime: 9.07s
% Computational Cost: add. (22097->567), mult. (13212->707), div. (0->0), fcn. (10256->10), ass. (0->358)
t276 = qJ(4) + qJ(5);
t266 = sin(t276);
t268 = cos(t276);
t213 = rSges(6,1) * t266 + rSges(6,2) * t268;
t277 = qJ(1) + qJ(2);
t271 = qJ(3) + t277;
t262 = cos(t271);
t280 = cos(qJ(4));
t501 = pkin(4) * t280;
t263 = pkin(3) + t501;
t229 = t262 * t263;
t261 = sin(t271);
t282 = -pkin(9) - pkin(8);
t497 = pkin(8) + t282;
t502 = pkin(3) * t262;
t130 = t261 * t497 - t229 + t502;
t450 = t262 * t266;
t231 = rSges(6,2) * t450;
t449 = t262 * t268;
t403 = rSges(6,1) * t449;
t335 = rSges(6,3) * t261 + t403;
t142 = -t231 + t335;
t438 = t130 - t142;
t456 = t261 * t266;
t227 = Icges(6,4) * t456;
t455 = t261 * t268;
t140 = -Icges(6,1) * t455 + Icges(6,5) * t262 + t227;
t228 = Icges(6,4) * t450;
t141 = Icges(6,1) * t449 + Icges(6,5) * t261 - t228;
t274 = qJD(4) + qJD(5);
t193 = t261 * t274;
t194 = t262 * t274;
t275 = qJD(1) + qJD(2);
t265 = qJD(3) + t275;
t259 = Icges(6,4) * t268;
t349 = -Icges(6,2) * t266 + t259;
t542 = Icges(6,1) * t266 + t259;
t554 = t542 + t349;
t294 = t193 * (-Icges(6,2) * t449 + t141 - t228) + t194 * (Icges(6,2) * t455 + t140 + t227) + t265 * t554;
t138 = Icges(6,6) * t262 - t261 * t349;
t139 = Icges(6,4) * t449 - Icges(6,2) * t450 + Icges(6,6) * t261;
t483 = Icges(6,4) * t266;
t209 = Icges(6,2) * t268 + t483;
t212 = Icges(6,1) * t268 - t483;
t520 = t193 * (t262 * t542 + t139) + t194 * (-t261 * t542 + t138) + t265 * (t209 - t212);
t561 = t294 * t266 + t268 * t520;
t560 = 2 * qJD(4);
t166 = t213 * t261;
t462 = t213 * t262;
t409 = rSges(6,1) * t455;
t556 = t213 * t274;
t397 = -t556 * t262 - t265 * t409;
t543 = rSges(6,2) * t456 + t262 * rSges(6,3);
t85 = t265 * t543 + t397;
t559 = t166 * t193 + t194 * t462 + t262 * t85;
t278 = sin(qJ(4));
t448 = t262 * t278;
t270 = Icges(5,4) * t280;
t350 = -Icges(5,2) * t278 + t270;
t541 = Icges(5,1) * t278 + t270;
t421 = t541 + t350;
t484 = Icges(5,4) * t278;
t246 = Icges(5,2) * t280 + t484;
t249 = Icges(5,1) * t280 - t484;
t422 = t246 - t249;
t558 = (t278 * t421 + t422 * t280) * t265;
t267 = sin(t277);
t269 = cos(t277);
t356 = rSges(3,1) * t267 + rSges(3,2) * t269;
t186 = t356 * t275;
t279 = sin(qJ(1));
t490 = pkin(1) * qJD(1);
t406 = t279 * t490;
t168 = t186 + t406;
t493 = rSges(6,1) * t268;
t214 = -rSges(6,2) * t266 + t493;
t457 = t261 * t265;
t557 = t265 * t166 - t194 * t214 - t213 * t457;
t500 = pkin(8) * t261;
t200 = t500 + t502;
t181 = t265 * t200;
t396 = t265 * t231 + t556 * t261;
t417 = qJD(4) * t278;
t392 = t261 * t417;
t236 = pkin(4) * t392;
t445 = t265 * t282;
t427 = t261 * t445 + t236;
t555 = t438 * t265 - t181 - t396 - t427;
t243 = rSges(5,2) * t448;
t447 = t262 * t280;
t410 = rSges(5,1) * t447;
t336 = -rSges(5,3) * t261 - t410;
t154 = -t243 - t336;
t135 = t265 * t154;
t416 = qJD(4) * t280;
t394 = rSges(5,1) * t392 + (t261 * t416 + t265 * t448) * rSges(5,2);
t553 = -t135 - t181 - t394;
t498 = pkin(3) - t263;
t545 = t497 * t262;
t552 = -t261 * t498 + t545;
t185 = t214 * t274;
t451 = t262 * t265;
t366 = t193 * t213 + t236;
t443 = t269 * t275;
t413 = pkin(2) * t443;
t320 = t366 - t413;
t281 = cos(qJ(1));
t405 = t281 * t490;
t296 = t320 - t405;
t58 = (-t200 + t438) * t265 + t296;
t551 = (t185 * t261 - t193 * t214 + t213 * t451 - t265 * t462) * t58;
t137 = Icges(6,5) * t449 - Icges(6,6) * t450 + Icges(6,3) * t261;
t62 = t262 * t137 + t139 * t456 - t141 * t455;
t345 = t209 * t266 - t268 * t542;
t207 = Icges(6,5) * t266 + Icges(6,6) * t268;
t466 = t207 * t262;
t91 = t261 * t345 + t466;
t550 = t193 * t62 + t91 * t265;
t151 = Icges(5,4) * t447 - Icges(5,2) * t448 + Icges(5,6) * t261;
t241 = Icges(5,4) * t448;
t153 = Icges(5,1) * t447 + Icges(5,5) * t261 - t241;
t346 = t151 * t278 - t153 * t280;
t547 = t346 * t262;
t348 = t139 * t266 - t141 * t268;
t546 = t348 * t262;
t454 = t261 * t278;
t423 = rSges(5,2) * t454 + t262 * rSges(5,3);
t195 = rSges(4,1) * t262 - t261 * rSges(4,2);
t177 = t265 * t195;
t322 = t405 + t413;
t129 = -t322 - t177;
t492 = rSges(5,2) * t278;
t495 = rSges(5,1) * t280;
t354 = -t492 + t495;
t222 = t354 * qJD(4);
t252 = rSges(5,1) * t278 + rSges(5,2) * t280;
t540 = -t222 * t262 + t252 * t457;
t463 = t209 * t274;
t538 = -Icges(6,6) * t265 + t463;
t536 = -Icges(6,3) * t265 + t207 * t274;
t535 = -Icges(6,5) * t265 + t274 * t542;
t531 = -Icges(5,6) * t265 + qJD(4) * t246;
t100 = t261 * t531 - t350 * t451;
t329 = t265 * t249;
t529 = -Icges(5,5) * t265 + qJD(4) * t541;
t102 = t261 * t529 - t262 * t329;
t327 = t350 * t261;
t150 = Icges(5,6) * t262 - t327;
t240 = Icges(5,4) * t454;
t453 = t261 * t280;
t152 = -Icges(5,1) * t453 + Icges(5,5) * t262 + t240;
t106 = t150 * t280 + t152 * t278;
t245 = Icges(5,5) * t280 - Icges(5,6) * t278;
t325 = t245 * t261;
t148 = Icges(5,3) * t262 - t325;
t534 = qJD(4) * t106 + t100 * t278 - t102 * t280 - t148 * t265;
t101 = -t261 * t329 - t262 * t529;
t107 = t151 * t280 + t153 * t278;
t149 = Icges(5,5) * t447 - Icges(5,6) * t448 + Icges(5,3) * t261;
t99 = -t262 * t531 - t265 * t327;
t533 = qJD(4) * t107 - t101 * t280 - t149 * t265 + t278 * t99;
t244 = Icges(5,5) * t278 + Icges(5,6) * t280;
t532 = -Icges(5,3) * t265 + qJD(4) * t244;
t217 = t350 * qJD(4);
t218 = t249 * qJD(4);
t530 = t217 * t278 - t218 * t280 - t244 * t265 + (t246 * t280 + t278 * t541) * qJD(4);
t528 = t413 + t555;
t419 = qJD(4) * t261;
t189 = t252 * t419;
t386 = -pkin(3) - t495;
t509 = -rSges(5,3) - pkin(8);
t389 = t262 * t416;
t390 = t262 * t417;
t411 = rSges(5,1) * t453;
t395 = -rSges(5,1) * t390 - rSges(5,2) * t389 - t265 * t411;
t104 = t265 * t423 + t395;
t224 = pkin(3) * t457;
t159 = pkin(8) * t451 - t224;
t418 = qJD(4) * t262;
t393 = t252 * t418;
t506 = pkin(1) * qJD(1) ^ 2;
t264 = t279 * t506;
t503 = pkin(2) * t275 ^ 2;
t420 = t267 * t503 + t264;
t55 = t222 * t419 + (-t104 - t159 + t393) * t265 + t420;
t105 = t265 * t336 + t394;
t415 = t281 * t506;
t331 = -t269 * t503 - t415;
t315 = -t265 ^ 2 * t200 + t331;
t56 = qJD(4) * t540 + t105 * t265 + t315;
t342 = t411 - t423;
t134 = t265 * t342;
t499 = t262 * pkin(8);
t180 = t265 * (-pkin(3) * t261 + t499);
t338 = t134 - t180 + t393;
t444 = t267 * t275;
t414 = pkin(2) * t444;
t300 = t338 + t414;
t74 = t300 + t406;
t310 = t189 - t322;
t75 = (-t154 - t200) * t265 + t310;
t285 = (t56 * t386 + t55 * t509) * t261 + ((-t492 * t75 - t509 * t74) * t261 + (-t386 * t74 + t509 * t75) * t262) * t265;
t527 = t285 + (t189 + t553) * t74;
t430 = -Icges(5,2) * t447 + t153 - t241;
t432 = t262 * t541 + t151;
t525 = t430 * t278 + t280 * t432;
t431 = Icges(5,2) * t453 + t152 + t240;
t433 = -t261 * t541 + t150;
t524 = -t278 * t431 - t280 * t433;
t367 = -t212 * t274 + t463;
t368 = t554 * t274;
t523 = -t207 * t265 + t266 * t368 + t268 * t367;
t326 = t349 * t265;
t373 = t141 * t274 - t261 * t326 - t262 * t538;
t328 = t212 * t265;
t375 = t139 * t274 + t261 * t328 + t262 * t535;
t522 = -t137 * t265 + t266 * t373 + t268 * t375;
t208 = Icges(6,5) * t268 - Icges(6,6) * t266;
t136 = Icges(6,3) * t262 - t208 * t261;
t374 = t140 * t274 + t261 * t538 - t262 * t326;
t376 = t138 * t274 - t261 * t535 + t262 * t328;
t521 = -t136 * t265 + t266 * t374 + t268 * t376;
t155 = t265 * t193;
t519 = -t155 / 0.2e1;
t156 = t265 * t194;
t518 = t156 / 0.2e1;
t517 = -t193 / 0.2e1;
t516 = t193 / 0.2e1;
t515 = -t194 / 0.2e1;
t514 = t194 / 0.2e1;
t513 = t261 / 0.2e1;
t512 = t262 / 0.2e1;
t511 = -t265 / 0.2e1;
t510 = t265 / 0.2e1;
t508 = pkin(1) * t279;
t507 = pkin(1) * t281;
t505 = pkin(2) * t267;
t504 = pkin(2) * t269;
t86 = -t265 * t335 + t396;
t385 = t262 * t498;
t96 = (t385 + t500) * t265 + t427;
t496 = -t96 - t86;
t488 = t261 * t75;
t160 = t207 * t261;
t92 = -t262 * t345 + t160;
t487 = t92 * t265;
t170 = t244 * t261;
t343 = t278 * t246 - t280 * t541;
t110 = -t262 * t343 + t170;
t474 = t110 * t265;
t471 = t138 * t266;
t470 = t140 * t268;
t469 = t150 * t278;
t468 = t152 * t280;
t464 = t208 * t265;
t460 = t244 * t262;
t459 = t245 * t265;
t178 = t252 * t261;
t458 = t252 * t262;
t452 = t262 * t130;
t442 = t262 * t136 + t138 * t456;
t441 = -t261 * t136 - t140 * t449;
t440 = t262 * t148 + t150 * t454;
t439 = t261 * t148 + t152 * t447;
t429 = -t262 * t445 - t263 * t457;
t157 = rSges(4,1) * t457 + rSges(4,2) * t451;
t412 = (qJD(4) ^ 2) * t501;
t404 = pkin(4) * t417;
t388 = -t457 / 0.2e1;
t387 = t451 / 0.2e1;
t382 = -t419 / 0.2e1;
t380 = -t418 / 0.2e1;
t379 = t418 / 0.2e1;
t215 = rSges(3,1) * t269 - t267 * rSges(3,2);
t377 = -t263 - t493;
t369 = -t149 - t468;
t365 = pkin(4) * t390;
t364 = t224 - t395;
t363 = -pkin(3) - t377;
t360 = t423 + t499;
t187 = -rSges(3,1) * t443 + rSges(3,2) * t444;
t158 = -rSges(4,1) * t451 + rSges(4,2) * t457;
t359 = -t262 * t282 + t543;
t355 = -rSges(4,1) * t261 - rSges(4,2) * t262;
t88 = t139 * t268 + t141 * t266;
t347 = -t468 + t469;
t341 = -t229 + t231 - t403;
t340 = t409 - t543;
t63 = -t138 * t450 - t441;
t66 = t262 * t149 + t151 * t454 - t153 * t453;
t339 = -t195 - t504;
t337 = t545 - t543;
t169 = -t215 * t275 - t405;
t65 = -t152 * t453 + t440;
t333 = (t261 * t66 + t262 * t65) * qJD(4);
t67 = -t150 * t448 + t439;
t68 = t149 * t261 - t547;
t332 = (t261 * t68 + t262 * t67) * qJD(4);
t330 = t360 - t505;
t323 = t406 + t414;
t321 = t243 - t410 - t502;
t319 = t359 - t505;
t318 = t355 - t505;
t317 = t364 + t414;
t313 = t160 * t194 - t193 * t466 + t464;
t311 = t158 - t413;
t309 = t341 - t504;
t308 = -t261 * t464 - t262 * t536 + t265 * t348;
t307 = -t262 * t464 + t536 * t261 + (-t470 + t471) * t265;
t306 = -t262 * t532 + (-t325 + t346) * t265;
t305 = -t245 * t451 + t261 * t532 + t265 * t347;
t304 = t208 * t274 + t265 * t345;
t303 = t245 * qJD(4) + t265 * t343;
t301 = t321 - t504;
t299 = t365 - t397 - t429;
t13 = t307 * t261 - t262 * t521;
t14 = t308 * t261 - t262 * t522;
t15 = t261 * t521 + t307 * t262;
t16 = t261 * t522 + t308 * t262;
t61 = -t140 * t455 + t442;
t28 = t194 * t61 + t550;
t64 = t137 * t261 - t546;
t29 = t193 * t64 + t194 * t63 + t487;
t40 = -t266 * t376 + t268 * t374;
t41 = -t266 * t375 + t268 * t373;
t44 = t304 * t261 - t262 * t523;
t45 = t261 * t523 + t304 * t262;
t87 = t138 * t268 + t140 * t266;
t298 = (t13 * t194 + t14 * t193 - t155 * t63 + t156 * t64 + t265 * t44) * t513 + (t313 * t261 - t561 * t262) * t517 + (t561 * t261 + t313 * t262) * t515 + (t15 * t194 - t155 * t61 + t156 * t62 + t16 * t193 + t265 * t45) * t512 + (-t266 * t520 + t268 * t294) * t511 + t28 * t388 + t29 * t387 + ((t265 * t64 + t13) * t262 + (-t265 * t63 + t14) * t261) * t516 + (t261 * t62 + t262 * t61) * t519 + (t261 * t64 + t262 * t63) * t518 + ((t265 * t62 + t15) * t262 + (-t265 * t61 + t16) * t261) * t514 + ((t265 * t88 + t40) * t262 + (-t265 * t87 + t41) * t261) * t510;
t297 = t194 * t213 - t180 + t365 + (t340 + t552) * t265;
t90 = (t262 * t154 + t261 * t342) * qJD(4);
t291 = t299 + t414;
t290 = t297 + t414;
t289 = (-t137 - t470) * t261 + t546 + t442;
t109 = t261 * t343 + t460;
t108 = t109 * t265;
t34 = t108 + t333;
t35 = t332 + t474;
t48 = -qJD(4) * t347 + t100 * t280 + t102 * t278;
t49 = -qJD(4) * t346 + t101 * t278 + t280 * t99;
t53 = t303 * t261 - t262 * t530;
t54 = t261 * t530 + t303 * t262;
t287 = (t108 + ((t440 + t68 + t547) * t262 + (-t67 + (t369 - t469) * t262 + t66 + t439) * t261) * qJD(4)) * t382 + ((t64 + t289) * t194 + t550) * t517 + (t87 + t91) * t519 + (t88 + t92) * t518 + (-t487 + (t62 + (-t137 + t471) * t262 - t348 * t261 + t441) * t194 + (-t61 + t289) * t193 + t29) * t515 + (t40 + t45) * t514 + (t35 - t474 + ((t66 + (-t149 + t469) * t262 - t439) * t262 + (t261 * t369 + t440 - t65) * t261) * qJD(4)) * t380 + (t48 + t54) * t379 + (t41 + t44 + t28) * t516 + (t49 + t53 + t34) * t419 / 0.2e1 + ((t106 + t109) * t382 + (t107 + t110) * t379 - qJD(4) * t343 + t217 * t280 + t218 * t278 - t266 * t367 + t268 * t368) * t265;
t95 = t224 + (-pkin(8) * t265 - t404) * t262 + t429;
t32 = t261 * t412 + t156 * t213 + t185 * t193 + (-t159 - t85 - t95 + t365) * t265 + t420;
t33 = -t262 * t412 + t155 * t213 - t185 * t194 + (t236 - t496) * t265 + t315;
t57 = t290 + t406;
t286 = (t32 * (-rSges(6,3) + t282) + t33 * t377) * t261 + (-t58 * t543 - t57 * (-t335 - t229)) * t265;
t176 = t265 * t355;
t144 = t187 * t275 - t415;
t143 = t186 * t275 + t264;
t128 = -t176 + t323;
t122 = t262 * t142;
t112 = t158 * t265 + t331;
t111 = t157 * t265 + t420;
t50 = t194 * t142 + t193 * t340 + (t552 * t261 - t452) * qJD(4);
t12 = -t155 * t142 + t194 * t85 + t156 * t340 - t193 * t86 + ((t265 * t545 + t95) * t262 + (-t96 + (t130 - t385) * t265) * t261) * qJD(4);
t1 = [m(3) * (t143 * (-t215 - t507) + t144 * (-t356 - t508) + (t169 - t187 + t405) * t168) + t287 + (t32 * (t309 - t507) + t58 * (t291 + t406) + t33 * (t319 - t508) + t286 + (t296 - t58 + t405 + t528) * t57) * m(6) + (t55 * (t301 - t507) + t75 * (t317 + t406) + t56 * (t330 - t508) + t285 + (t310 - t75 + t322 + t553) * t74) * m(5) + (t111 * (t339 - t507) + t129 * (t323 + t157) + t112 * (t318 - t508) + (-t311 + t405) * t128) * m(4); t287 + (t32 * t309 + t33 * t319 + t286 + (-t290 + t291) * t58 + (t320 + t528) * t57) * m(6) + (t55 * t301 + t56 * t330 + (-t300 + t317) * t75 + t527) * m(5) + (t111 * t339 + t112 * t318 + (-t177 - t413 - t311) * t128 + (t176 + t157) * t129) * m(4) + (-(t168 * t215 + t169 * t356) * t275 - t143 * t215 - t144 * t356 - t168 * t187 + t169 * t186) * m(3); t287 + (t32 * t341 + t33 * t359 + t286 + (-t297 + t299) * t58 + (t366 + t555) * t57) * m(6) + (t55 * t321 + t56 * t360 + (-t338 + t364) * t75 + t527) * m(5) + (-(t128 * t195 - t129 * t355) * t265 - t111 * t195 + t112 * t355 - t128 * t158 + t129 * t157) * m(4); ((-t419 * t460 + t459) * t261 + (-t558 + (t524 * t262 + (t170 - t525) * t261) * qJD(4)) * t262) * t382 + t298 + ((t107 * t265 + t48) * t262 + (-t106 * t265 + t49) * t261) * t510 + ((t170 * t418 + t459) * t262 + (t558 + (t525 * t261 + (-t460 - t524) * t262) * qJD(4)) * t261) * t380 + ((-t278 * t422 + t280 * t421) * t265 + ((t261 * t430 + t262 * t431) * t280 + (-t261 * t432 - t262 * t433) * t278) * qJD(4)) * t511 + (t265 * t53 + ((t305 * t261 - t262 * t534 + t265 * t68) * t262 + (t306 * t261 - t262 * t533 - t265 * t67) * t261) * t560) * t513 + (t265 * t54 + ((t261 * t534 + t305 * t262 + t265 * t66) * t262 + (t261 * t533 + t306 * t262 - t265 * t65) * t261) * t560) * t512 + (t333 + t34) * t388 + (t332 + t35) * t387 + (t12 * (-t452 + t122 + (t261 * t363 + t337) * t261) + (t261 * t32 - t262 * t33) * (pkin(4) * t278 + t213) + (-(-t261 ^ 2 - t262 ^ 2) * t404 + t262 * t95 + t496 * t261 + (t337 * t262 + (t262 * t363 + t438) * t261) * t265 + t559) * t50 + t551 + (-pkin(4) * t389 - (-pkin(4) * t416 - t185) * t262 + t557) * t57) * m(6) + (t75 * t252 * t451 + t55 * t178 + t222 * t488 + 0.2e1 * t90 * ((-t105 - t135) * t261 + (t104 + t134) * t262) - t56 * t458 - t540 * t74 - (-t178 * t74 + t458 * t75) * t265 - (t90 * (-t178 * t261 - t262 * t458) + (t262 * t74 + t488) * t354) * qJD(4)) * m(5); t298 + (t12 * (t261 * t340 + t122) + t32 * t166 - t33 * t462 + (-t261 * t86 + (-t261 * t142 + t262 * t340) * t265 + t559) * t50 + t551 + (t185 * t262 + t557) * t57) * m(6);];
tauc = t1(:);
