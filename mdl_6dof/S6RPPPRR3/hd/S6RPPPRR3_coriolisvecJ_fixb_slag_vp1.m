% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:28
% EndTime: 2019-03-09 01:33:51
% DurationCPUTime: 19.30s
% Computational Cost: add. (20971->786), mult. (35811->1033), div. (0->0), fcn. (40059->10), ass. (0->372)
t498 = sin(pkin(9));
t499 = cos(pkin(9));
t522 = sin(qJ(1));
t523 = cos(qJ(1));
t268 = t523 * t498 - t522 * t499;
t300 = pkin(10) + qJ(5);
t289 = sin(t300);
t441 = qJD(6) * t289;
t267 = -t498 * t522 - t499 * t523;
t445 = qJD(5) * t267;
t192 = t268 * t441 + t445;
t444 = qJD(5) * t268;
t193 = -t267 * t441 + t444;
t479 = t267 * t289;
t305 = cos(qJ(6));
t290 = cos(t300);
t304 = sin(qJ(6));
t473 = t290 * t304;
t188 = -t267 * t305 + t268 * t473;
t173 = Icges(7,4) * t188;
t472 = t290 * t305;
t187 = t267 * t304 + t268 * t472;
t476 = t268 * t289;
t100 = Icges(7,1) * t187 + Icges(7,5) * t476 - t173;
t190 = t267 * t473 + t268 * t305;
t191 = -t267 * t472 + t268 * t304;
t495 = Icges(7,4) * t187;
t97 = -Icges(7,2) * t188 + Icges(7,6) * t476 + t495;
t569 = t191 * t100 + t190 * t97;
t94 = Icges(7,5) * t187 - Icges(7,6) * t188 + Icges(7,3) * t476;
t25 = t479 * t94 - t569;
t174 = Icges(7,4) * t190;
t102 = Icges(7,1) * t191 - Icges(7,5) * t479 + t174;
t96 = Icges(7,5) * t191 + Icges(7,6) * t190 - Icges(7,3) * t479;
t494 = Icges(7,4) * t191;
t99 = Icges(7,2) * t190 - Icges(7,6) * t479 + t494;
t26 = t191 * t102 + t190 * t99 - t96 * t479;
t440 = qJD(6) * t290;
t278 = qJD(1) + t440;
t493 = Icges(7,4) * t304;
t277 = t289 * t493;
t474 = t289 * t305;
t490 = Icges(7,5) * t290;
t207 = -Icges(7,1) * t474 + t277 + t490;
t372 = Icges(7,5) * t305 - Icges(7,6) * t304;
t320 = -Icges(7,3) * t290 + t289 * t372;
t492 = Icges(7,4) * t305;
t375 = -Icges(7,2) * t304 + t492;
t321 = -Icges(7,6) * t290 + t289 * t375;
t61 = -t190 * t321 + t191 * t207 + t320 * t479;
t11 = -t192 * t25 + t193 * t26 + t61 * t278;
t374 = Icges(6,5) * t290 - Icges(6,6) * t289;
t146 = Icges(6,3) * t267 + t268 * t374;
t496 = Icges(6,4) * t290;
t377 = -Icges(6,2) * t289 + t496;
t149 = Icges(6,6) * t267 + t268 * t377;
t497 = Icges(6,4) * t289;
t380 = Icges(6,1) * t290 - t497;
t152 = Icges(6,5) * t267 + t268 * t380;
t559 = t149 * t289 - t152 * t290;
t55 = t267 * t146 - t268 * t559;
t382 = t100 * t305 - t304 * t97;
t39 = t289 * t382 - t290 * t94;
t23 = t100 * t187 - t188 * t97 + t476 * t94;
t292 = qJD(2) * t523;
t338 = t523 * pkin(1) + t522 * qJ(2);
t242 = qJD(1) * t338 - t292;
t390 = rSges(7,1) * t187 - t188 * rSges(7,2);
t105 = -rSges(7,3) * t476 - t390;
t389 = rSges(7,1) * t305 - rSges(7,2) * t304;
t325 = -rSges(7,3) * t290 + t289 * t389;
t520 = t289 * pkin(5);
t397 = -pkin(8) * t290 + t520;
t249 = qJD(4) * t268;
t291 = qJD(2) * t522;
t454 = t249 + t291;
t566 = -t278 * t105 + t192 * t325 + t397 * t445 + t454;
t552 = t187 * t207 + t188 * t321 - t320 * t476;
t565 = t192 * t23 + t552 * t278;
t148 = Icges(6,3) * t268 - t267 * t374;
t562 = t268 * t148;
t419 = qJD(1) * t523;
t558 = -pkin(2) * t419 - t242;
t373 = Icges(6,5) * t289 + Icges(6,6) * t290;
t168 = t373 * t267;
t376 = Icges(6,2) * t290 + t497;
t379 = Icges(6,1) * t289 + t496;
t364 = -t289 * t376 + t290 * t379;
t108 = t268 * t364 + t168;
t561 = qJD(1) * t108;
t75 = -t149 * t290 - t152 * t289;
t531 = -t192 / 0.2e1;
t167 = t373 * t268;
t442 = qJD(5) * t290;
t247 = t268 * qJD(1);
t483 = t247 * t289;
t349 = t267 * t442 - t483;
t443 = qJD(5) * t289;
t423 = t267 * t443;
t482 = t247 * t290;
t348 = t423 + t482;
t302 = cos(pkin(10));
t557 = rSges(5,2) * sin(pkin(10)) - rSges(5,1) * t302;
t553 = t187 * t102 - t188 * t99;
t551 = 0.2e1 * qJD(5);
t519 = t290 * pkin(5);
t398 = pkin(8) * t289 + t519;
t550 = t398 * t268;
t283 = pkin(4) * t302 + pkin(3);
t516 = pkin(3) - t283;
t549 = t516 * t268;
t294 = t523 * qJ(2);
t436 = t522 * pkin(1);
t272 = t436 - t294;
t251 = t267 * qJ(4);
t547 = t268 * pkin(3) + t251;
t548 = t547 - t272;
t435 = t522 * pkin(2);
t331 = t268 * rSges(4,1) - t267 * rSges(4,2) - t435;
t546 = -t272 + t331;
t339 = t523 * rSges(3,1) + t522 * rSges(3,3);
t545 = t338 + t339;
t246 = t267 * qJD(1);
t407 = t247 * pkin(3) + t246 * qJ(4);
t543 = t407 + t249;
t270 = qJD(1) * t272;
t448 = t291 - t270;
t449 = qJ(2) * t419 + t291;
t542 = t449 - t448;
t513 = rSges(6,2) * t289;
t392 = rSges(6,1) * t290 - t513;
t155 = t267 * rSges(6,3) + t268 * t392;
t303 = -pkin(7) - qJ(4);
t481 = t247 * t303;
t540 = t481 + t558;
t459 = -t246 * t303 + t247 * t283;
t539 = -qJD(1) * (-t549 + (-qJ(4) - t303) * t267) + t449 + t459 + t249;
t465 = t376 * t268 - t152;
t467 = -t379 * t268 - t149;
t538 = t289 * t465 + t290 * t467;
t204 = -Icges(7,3) * t289 - t290 * t372;
t365 = -t207 * t305 - t304 * t321;
t381 = -t102 * t305 + t304 * t99;
t537 = t193 * (-t320 * t267 - t381) - t192 * (-t320 * t268 - t382) - t278 * (t204 + t365);
t226 = (Icges(7,1) * t304 + t492) * t289;
t315 = t193 * (Icges(7,1) * t190 - t494 - t99) - t192 * (Icges(7,1) * t188 + t495 + t97) - t278 * (-t321 - t226);
t306 = qJD(1) ^ 2;
t421 = t267 * t440;
t142 = t247 * t441 + (t246 - t421) * qJD(5);
t536 = t142 / 0.2e1;
t420 = t268 * t440;
t143 = -t246 * t441 + (t247 - t420) * qJD(5);
t535 = t143 / 0.2e1;
t534 = -t193 / 0.2e1;
t533 = t193 / 0.2e1;
t532 = t192 / 0.2e1;
t530 = t246 / 0.2e1;
t529 = t247 / 0.2e1;
t526 = -t278 / 0.2e1;
t525 = t278 / 0.2e1;
t524 = rSges(7,3) + pkin(8);
t484 = t246 * t289;
t351 = t268 * t442 + t484;
t350 = -t246 * t290 + t268 * t443;
t328 = -qJD(6) * t267 + t350;
t395 = t247 + t420;
t82 = -t304 * t328 + t305 * t395;
t83 = t304 * t395 + t305 * t328;
t43 = Icges(7,5) * t83 + Icges(7,6) * t82 - Icges(7,3) * t351;
t45 = Icges(7,4) * t83 + Icges(7,2) * t82 - Icges(7,6) * t351;
t47 = Icges(7,1) * t83 + Icges(7,4) * t82 - Icges(7,5) * t351;
t7 = (qJD(5) * t382 + t43) * t290 + (qJD(5) * t94 + t304 * t45 - t305 * t47 + (-t100 * t304 - t305 * t97) * qJD(6)) * t289;
t518 = t7 * t192;
t327 = qJD(6) * t268 + t348;
t396 = t246 + t421;
t84 = -t304 * t327 + t305 * t396;
t85 = t304 * t396 + t305 * t327;
t44 = Icges(7,5) * t85 + Icges(7,6) * t84 - Icges(7,3) * t349;
t46 = Icges(7,4) * t85 + Icges(7,2) * t84 - Icges(7,6) * t349;
t48 = Icges(7,1) * t85 + Icges(7,4) * t84 - Icges(7,5) * t349;
t8 = (qJD(5) * t381 + t44) * t290 + (-qJD(5) * t96 + t304 * t46 - t305 * t48 + (t102 * t304 + t305 * t99) * qJD(6)) * t289;
t517 = t8 * t193;
t512 = rSges(7,3) * t289;
t477 = t267 * t303;
t345 = -t251 - t435 - t477 + t548 - t549;
t31 = (t550 + t345) * qJD(1) + t566;
t510 = t246 * t31;
t391 = rSges(6,1) * t289 + rSges(6,2) * t290;
t53 = t391 * t445 + (t155 + t345) * qJD(1) + t454;
t509 = t246 * t53;
t508 = t247 * rSges(5,3);
t507 = t247 * rSges(6,3);
t506 = t267 * rSges(5,3);
t504 = t39 * t143;
t40 = t289 * t381 + t290 * t96;
t503 = t40 * t142;
t111 = pkin(5) * t348 - pkin(8) * t349;
t50 = t85 * rSges(7,1) + t84 * rSges(7,2) - rSges(7,3) * t349;
t502 = t111 + t50;
t151 = Icges(6,6) * t268 - t267 * t377;
t485 = t151 * t289;
t478 = t267 * t290;
t475 = t268 * t290;
t471 = t105 - t550;
t106 = t191 * rSges(7,1) + t190 * rSges(7,2) - rSges(7,3) * t479;
t185 = -pkin(5) * t478 - pkin(8) * t479;
t470 = t106 + t185;
t154 = Icges(6,5) * t268 - t267 * t380;
t469 = -t267 * t148 - t154 * t475;
t468 = -t154 * t478 + t562;
t466 = -t379 * t267 + t151;
t464 = t376 * t267 + t154;
t233 = t247 * qJ(4);
t439 = t267 * qJD(4);
t463 = t246 * pkin(3) - t233 - t242 + t439;
t210 = -t290 * t389 - t512;
t228 = (rSges(7,1) * t304 + rSges(7,2) * t305) * t289;
t166 = qJD(5) * t210 + qJD(6) * t228;
t245 = t398 * qJD(5);
t462 = -t166 + t245;
t460 = -t325 - t397;
t458 = rSges(6,1) * t482 + t246 * rSges(6,3);
t457 = -t267 * t283 - t268 * t303;
t418 = qJD(1) * t522;
t456 = (-pkin(1) * t418 + t291 + t449) * qJD(1);
t455 = t247 * rSges(4,1) - t246 * rSges(4,2);
t453 = -rSges(6,1) * t478 + t268 * rSges(6,3);
t452 = -t267 * rSges(4,1) - t268 * rSges(4,2);
t451 = -t376 + t380;
t450 = -t377 - t379;
t446 = qJD(1) * t374;
t24 = -t96 * t476 - t553;
t438 = -t523 / 0.2e1;
t437 = t522 / 0.2e1;
t298 = t523 * pkin(2);
t432 = t522 * rSges(3,1);
t427 = -t246 * t516 + t233 + t463 + t481;
t424 = t298 + t338;
t417 = qJD(5) * t441;
t413 = t445 / 0.2e1;
t412 = -t444 / 0.2e1;
t410 = -t443 / 0.2e1;
t409 = -t442 / 0.2e1;
t256 = t267 * pkin(3);
t406 = qJ(4) * t268 - t256;
t197 = qJD(1) * t547;
t405 = t197 + t249 + t448;
t404 = qJD(1) * t406 - t558;
t402 = qJD(6) * t410;
t399 = t83 * rSges(7,1) + t82 * rSges(7,2);
t394 = t246 * rSges(4,1) + t247 * rSges(4,2);
t388 = -t23 * t268 - t24 * t267;
t387 = -t25 * t268 - t26 * t267;
t386 = -t267 * t40 - t268 * t39;
t157 = rSges(6,2) * t479 + t453;
t384 = -qJD(1) * (-t406 + t457) - t404;
t334 = -qJD(1) * t157 + t384;
t54 = t391 * t444 - t334 - t439;
t385 = t267 * t53 + t268 * t54;
t383 = t424 + t457;
t378 = Icges(7,1) * t305 - t493;
t371 = -t105 * t267 + t106 * t268;
t367 = -t154 * t290 + t485;
t366 = -t155 * t268 + t157 * t267;
t285 = qJD(1) * t292;
t363 = -t298 * t306 + t285;
t362 = t246 * rSges(5,3) - t247 * t557;
t164 = t268 * rSges(5,3) + t267 * t557;
t361 = pkin(3) - t557;
t57 = -t146 * t268 - t149 * t479 + t152 * t478;
t103 = t289 * t365 - t290 * t320;
t224 = (Icges(7,5) * t304 + Icges(7,6) * t305) * t289;
t161 = qJD(5) * t204 + qJD(6) * t224;
t206 = -Icges(7,6) * t289 - t290 * t375;
t162 = (Icges(7,2) * t305 + t493) * t441 + t206 * qJD(5);
t208 = -Icges(7,5) * t289 - t290 * t378;
t163 = qJD(5) * t208 + qJD(6) * t226;
t36 = (qJD(5) * t365 + t161) * t290 + (qJD(5) * t320 + t162 * t304 - t163 * t305 + (t207 * t304 - t305 * t321) * qJD(6)) * t289;
t359 = -t103 * t417 + t36 * t278;
t357 = t283 + t392;
t356 = -t289 * t43 + t442 * t94;
t355 = -t289 * t44 - t442 * t96;
t347 = -t306 * t435 + t456;
t346 = qJD(4) * t246 + t363;
t56 = t151 * t476 + t469;
t344 = (-t267 * t55 + t268 * t56) * qJD(5);
t58 = t151 * t479 + t468;
t343 = (-t267 * t57 + t268 * t58) * qJD(5);
t342 = -qJD(1) * t164 - t404;
t341 = -t436 - t435;
t340 = -t432 - t436;
t336 = -t192 * t94 - t193 * t96 + t278 * t320;
t335 = -(Icges(7,5) * t188 + Icges(7,6) * t187) * t192 + (Icges(7,5) * t190 - Icges(7,6) * t191) * t193 + t224 * t278;
t333 = t294 + t341;
t332 = t289 * t464 + t290 * t466;
t330 = qJD(1) * t543 + qJD(4) * t247 + t347;
t326 = t289 * t335;
t324 = qJD(1) * (-t407 + t459) + t330;
t12 = qJD(1) * t111 + t142 * t325 - t193 * t166 + t278 * t50 + (-t106 * t441 + t245 * t268 + t246 * t397) * qJD(5) + t324;
t110 = pkin(5) * t350 - pkin(8) * t351;
t49 = -rSges(7,3) * t351 + t399;
t13 = -t143 * t325 - t192 * t166 - t278 * t49 + (t105 * t441 + t245 * t267 - t247 * t397) * qJD(5) + (-t110 + t427) * qJD(1) + t346;
t308 = -qJD(1) * t185 - t106 * t278 - t193 * t325 - t397 * t444 + t384;
t32 = -t439 - t308;
t323 = t12 * t267 - t13 * t268 - t247 * t32 - t510;
t322 = t289 * t378 - t490;
t319 = (t289 * t450 + t290 * t451) * qJD(1);
t92 = rSges(6,1) * t350 + rSges(6,2) * t351 + t507;
t93 = rSges(6,1) * t423 + rSges(6,2) * t349 + t458;
t318 = -t155 * t246 - t157 * t247 + t267 * t93 + t268 * t92;
t312 = (-Icges(7,2) * t191 + t102 + t174) * t193 - (Icges(7,2) * t187 - t100 + t173) * t192 + (Icges(7,2) * t474 + t207 + t277) * t278;
t240 = t377 * qJD(5);
t241 = t380 * qJD(5);
t311 = -t240 * t289 + t241 * t290 + (-t289 * t379 - t290 * t376) * qJD(5);
t310 = -t268 * t557 - t435 + t506;
t30 = t105 * t193 + t106 * t192 - qJD(3) + (t185 * t267 - t268 * t550) * qJD(5);
t309 = t30 * t371 - (t267 * t32 - t268 * t31) * t325;
t307 = t537 * t289;
t296 = t523 * rSges(3,3);
t287 = rSges(3,3) * t419;
t273 = t432 - t296;
t243 = t392 * qJD(5);
t239 = t374 * qJD(5);
t201 = t291 + (-t272 - t273) * qJD(1);
t184 = t397 * t267;
t182 = t397 * t268;
t180 = -qJD(1) * t242 - t306 * t339 + t285;
t179 = qJD(1) * (-rSges(3,1) * t418 + t287) + t456;
t176 = t391 * t267;
t175 = t391 * t268;
t158 = qJD(1) * t546 + t291;
t141 = t325 * t267;
t140 = t325 * t268;
t139 = t322 * t267;
t138 = t322 * t268;
t137 = t321 * t267;
t136 = t321 * t268;
t129 = (-t242 + t394) * qJD(1) + t363;
t128 = qJD(1) * t455 + t347;
t119 = rSges(7,1) * t190 - rSges(7,2) * t191;
t118 = rSges(7,1) * t188 + rSges(7,2) * t187;
t109 = t267 * t364 - t167;
t107 = t109 * qJD(1);
t87 = Icges(6,5) * t348 + Icges(6,6) * t349 + Icges(6,3) * t246;
t86 = Icges(6,5) * t350 + Icges(6,6) * t351 + Icges(6,3) * t247;
t76 = t151 * t290 + t154 * t289;
t74 = -t439 - t342;
t73 = (t310 + t548) * qJD(1) + t454;
t64 = qJD(5) * t366 - qJD(3);
t63 = (-t246 * t557 + t463 - t508) * qJD(1) + t346;
t62 = qJD(1) * t362 + t330;
t42 = -t239 * t268 - t246 * t373 - t247 * t364 + t267 * t311;
t41 = t239 * t267 + t246 * t364 - t247 * t373 + t268 * t311;
t38 = qJD(5) * t367 - t289 * (Icges(6,1) * t348 + Icges(6,4) * t349 + Icges(6,5) * t246) - t290 * (Icges(6,4) * t348 + Icges(6,2) * t349 + Icges(6,6) * t246);
t37 = -qJD(5) * t559 - t289 * (Icges(6,1) * t350 + Icges(6,4) * t351 + Icges(6,5) * t247) - t290 * (Icges(6,4) * t350 + Icges(6,2) * t351 + Icges(6,6) * t247);
t35 = (t243 * t267 - t247 * t391) * qJD(5) + (-t92 + t427) * qJD(1) + t346;
t34 = qJD(1) * t93 + (t243 * t268 + t246 * t391) * qJD(5) + t324;
t27 = t318 * qJD(5);
t22 = t107 + t343;
t21 = t344 + t561;
t20 = -t161 * t479 + t162 * t190 + t163 * t191 + t207 * t85 + t320 * t349 - t321 * t84;
t19 = -t161 * t476 + t162 * t188 - t163 * t187 + t207 * t83 + t320 * t351 - t321 * t82;
t14 = t103 * t278 - t192 * t39 + t193 * t40;
t10 = t193 * t24 - t565;
t9 = t105 * t142 - t106 * t143 + t193 * t49 + t192 * t50 + (t110 * t268 + t111 * t267 - t185 * t247 - t246 * t550) * qJD(5);
t6 = t102 * t85 + t190 * t46 + t191 * t48 + t267 * t355 + t483 * t96 + t84 * t99;
t5 = -t100 * t85 + t190 * t45 + t191 * t47 + t267 * t356 - t483 * t94 - t84 * t97;
t4 = t102 * t83 - t187 * t48 + t188 * t46 + t268 * t355 - t484 * t96 + t82 * t99;
t3 = -t100 * t83 - t187 * t47 + t188 * t45 + t268 * t356 + t484 * t94 - t82 * t97;
t2 = t142 * t26 + t143 * t25 - t192 * t5 + t193 * t6 + t20 * t278 - t417 * t61;
t1 = t142 * t24 + t143 * t23 + t19 * t278 - t192 * t3 + t193 * t4 + t417 * t552;
t15 = [((t25 + (-t267 * t94 + t268 * t96) * t289 + t553 + t569) * t193 + t10 + t565) * t534 + t107 * t413 + t11 * t532 + (t12 * (t383 + t470) - (-t289 * t524 - t283 - t519) * t510 + (t333 + t390 - t477 + (t283 + t398 + t512) * t268) * t13 + (-t197 + t270 + t502 + (t341 + t435 - t550) * qJD(1) + t539 - t566) * t32 + (-t399 + (t290 * t524 - t520) * t444 - t308 + t540) * t31) * m(7) + (t240 * t290 + t241 * t289) * qJD(1) + t503 / 0.2e1 + t504 / 0.2e1 + t359 + t517 / 0.2e1 - t518 / 0.2e1 + (t129 * t546 + t128 * (t424 + t452) + (t394 + t558) * t158 + (t158 + t455 + (-t331 + t341) * qJD(1) + t542) * (qJD(1) * t452 - t558)) * m(4) + (t38 + t42) * t444 / 0.2e1 - (t37 + t41 + t22) * t445 / 0.2e1 + (t21 - t561) * t412 + (t34 * (t383 + t453) + t357 * t509 + (t357 * t268 + t333) * t35 + (-rSges(6,2) * t483 + pkin(2) * t418 - t405 + t458 + (-t155 + t341) * qJD(1) + t539) * t54 + (t35 * (rSges(6,3) - t303) + t34 * t513) * t267 + (-t334 - t507 + t540) * t53) * m(6) + t20 * t533 - t552 * t535 + t61 * t536 + (t63 * (t251 + t333 + t506) + t62 * (t164 - t256 + t424) + (t62 * qJ(4) + t361 * t63) * t268 + (t362 - t405 + t449 + (-t310 + t341) * qJD(1) + t543) * t74 + (t246 * t361 - t233 - t342 - t508 + t558) * t73) * m(5) + (t180 * (t294 + t296 + t340) + t179 * t545 + (-qJD(1) * t545 + t292) * t201 + (t201 + t287 + (t273 + t340) * qJD(1) + t542) * (qJD(1) * t339 + t242)) * m(3) + (t11 + t19) * t531 + (((t56 - t57 - t469) * t267 + t468 * t268) * t413 + (t109 - t76) * t530 + t364 * qJD(1) + (t108 - t75) * t529 + ((t57 + (t146 - t367) * t268) * t268 + (t58 + (t146 - t485) * t267 + t562 - t468) * t267) * t412) * qJD(5); 0.2e1 * (t12 * t438 + t13 * t437) * m(7) + 0.2e1 * (t34 * t438 + t35 * t437) * m(6) + 0.2e1 * (t437 * t63 + t438 * t62) * m(5) + 0.2e1 * (t128 * t438 + t129 * t437) * m(4) + 0.2e1 * (t179 * t438 + t180 * t437) * m(3); -m(6) * t27 - m(7) * t9; m(5) * (t246 * t73 + t247 * t74 - t267 * t62 + t268 * t63) - m(7) * t323 + (t247 * t54 - t267 * t34 + t268 * t35 + t509) * m(6) + (-m(5) * (t267 * t73 + t268 * t74) - m(7) * (t267 * t31 + t268 * t32) - m(6) * t385) * qJD(1); t14 * t441 / 0.2e1 + qJD(1) * (-t246 * t76 - t247 * t75 - t267 * t37 + t268 * t38) / 0.2e1 + ((t168 * t444 - t446) * t268 + (t319 + (-t538 * t267 + (-t167 + t332) * t268) * qJD(5)) * t267) * t412 - qJD(1) * ((t289 * t451 - t290 * t450) * qJD(1) + ((t267 * t465 - t268 * t464) * t290 + (-t267 * t467 + t268 * t466) * t289) * qJD(5)) / 0.2e1 + ((t167 * t445 + t446) * t267 + (t319 + (t332 * t268 + (-t168 - t538) * t267) * qJD(5)) * t268) * t413 + (-t267 * t39 + t268 * t40) * t402 + (t246 * t40 + t247 * t39 - t267 * t7 + t268 * t8) * t525 + (((t137 * t304 - t139 * t305 - t96) * t193 - (t136 * t304 - t138 * t305 + t94) * t192 + (t206 * t304 - t208 * t305 + t320) * t278 - t103 * qJD(6)) * t289 + (qJD(6) * t386 - t537) * t290) * t526 + (t23 * t247 + t24 * t246 - t267 * t3 + t268 * t4) * t531 + ((t137 * t188 - t139 * t187) * t193 - (t136 * t188 - t138 * t187) * t192 + (-t187 * t208 + t188 * t206) * t278 + (-t24 * t478 + t289 * t552) * qJD(6) + ((-qJD(6) * t23 + t336) * t290 + t307) * t268) * t532 + (t246 * t26 + t247 * t25 - t267 * t5 + t268 * t6) * t533 + ((t137 * t190 + t139 * t191) * t193 - (t136 * t190 + t138 * t191) * t192 + (t190 * t206 + t191 * t208) * t278 + (-t25 * t475 - t289 * t61) * qJD(6) + ((-qJD(6) * t26 + t336) * t290 + t307) * t267) * t534 + (-t23 * t267 + t24 * t268) * t535 + (-t25 * t267 + t26 * t268) * t536 + (t268 * t10 + t267 * t11) * t440 / 0.2e1 + ((-t30 * t470 + t31 * t460) * t247 - (-t30 * t471 + t32 * t460) * t246 + (-t12 * t460 + t32 * t462 + t9 * t471 + t30 * (t110 + t49)) * t268 + (-t13 * t460 + t30 * t502 + t31 * t462 + t470 * t9) * t267 - t31 * (-qJD(1) * t182 - t140 * t278 - t192 * t210 + t398 * t445) - t32 * (qJD(1) * t184 + t141 * t278 - t193 * t210 + t398 * t444) - t30 * (t140 * t193 + t141 * t192 + t182 * t444 + t184 * t445) - ((t105 * t31 - t106 * t32) * t289 + t309 * t290) * qJD(6)) * m(7) + (-(-t175 * t53 + t176 * t54) * qJD(1) - (t64 * (t175 * t268 + t176 * t267) + t385 * t392) * qJD(5) + t27 * t366 + t64 * t318 + t385 * t243 - (-t246 * t54 + t247 * t53 - t267 * t35 - t268 * t34) * t391) * m(6) - (qJD(1) * t41 + t1 + (-(-t146 * t247 - t246 * t559 - t267 * t86) * t267 + (t148 * t247 + t246 * t367 - t267 * t87) * t268 + t246 * t56 + t247 * t55) * t551) * t267 / 0.2e1 + (qJD(1) * t42 + t2 + (-(-t146 * t246 + t247 * t559 + t268 * t86) * t267 + (t148 * t246 - t247 * t367 + t268 * t87) * t268 + t246 * t58 + t247 * t57) * t551) * t268 / 0.2e1 + (t343 + t22 + t11) * t530 + (t344 + t21 + t10) * t529; -t2 * t479 / 0.2e1 + (t289 * t387 + t290 * t61) * t536 + ((qJD(5) * t387 + t20) * t290 + (-qJD(5) * t61 - t246 * t25 + t247 * t26 - t267 * t6 - t268 * t5) * t289) * t533 - t1 * t476 / 0.2e1 + (t289 * t388 - t290 * t552) * t535 + ((qJD(5) * t388 + t19) * t290 + (qJD(5) * t552 - t23 * t246 + t24 * t247 - t267 * t4 - t268 * t3) * t289) * t531 + t14 * t410 + t290 * (t359 + t503 + t504 + t517 - t518) / 0.2e1 + (t103 * t290 + t289 * t386) * t402 + ((qJD(5) * t386 + t36) * t290 + (-qJD(5) * t103 - t246 * t39 + t247 * t40 - t267 * t8 - t268 * t7) * t289) * t525 + (t190 * t312 + t191 * t315 - t267 * t326) * t534 + (-t187 * t315 + t188 * t312 - t268 * t326) * t532 + (t335 * t290 + (t312 * t304 - t305 * t315) * t289) * t526 + (t483 / 0.2e1 + t267 * t409) * t11 + (-t484 / 0.2e1 + t268 * t409) * t10 + ((qJD(5) * t309 - t13 * t105 + t12 * t106 - t31 * t49 + t32 * t50) * t290 + (t31 * (qJD(5) * t105 - t166 * t268) + t32 * (-qJD(5) * t106 + t166 * t267) + t9 * t371 + t30 * (t105 * t247 + t106 * t246 - t267 * t49 + t268 * t50) - t323 * t325) * t289 - t31 * (-t118 * t278 - t192 * t228) - t32 * (t119 * t278 - t193 * t228) - t30 * (t118 * t193 + t119 * t192)) * m(7);];
tauc  = t15(:);
