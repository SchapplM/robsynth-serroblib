% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:58
% EndTime: 2019-12-31 17:55:15
% DurationCPUTime: 14.18s
% Computational Cost: add. (6141->458), mult. (8746->569), div. (0->0), fcn. (6540->6), ass. (0->261)
t529 = Icges(6,4) + Icges(5,5);
t528 = Icges(5,6) - Icges(6,6);
t245 = sin(qJ(1));
t246 = cos(qJ(1));
t241 = pkin(7) + qJ(4);
t222 = cos(t241);
t221 = sin(t241);
t405 = Icges(5,4) * t221;
t297 = Icges(5,2) * t222 + t405;
t527 = t528 * t245 - t246 * t297;
t404 = Icges(5,4) * t222;
t299 = Icges(5,1) * t221 + t404;
t526 = -t529 * t245 + t246 * t299;
t385 = t222 * t246;
t387 = t221 * t246;
t519 = Icges(6,5) * t387 - Icges(6,3) * t385 + t527;
t197 = Icges(6,5) * t385;
t525 = Icges(6,1) * t387 - t197 + t526;
t401 = Icges(6,5) * t221;
t168 = Icges(6,1) * t222 + t401;
t170 = Icges(5,1) * t222 - t405;
t511 = t168 + t170;
t295 = Icges(5,5) * t221 + Icges(5,6) * t222;
t103 = Icges(5,3) * t246 + t245 * t295;
t296 = Icges(6,4) * t221 - Icges(6,6) * t222;
t105 = Icges(6,2) * t246 + t245 * t296;
t524 = t103 + t105;
t523 = -t528 * t221 + t529 * t222;
t495 = t525 * t221 - t519 * t222;
t293 = -Icges(6,3) * t222 + t401;
t522 = t293 - t297;
t215 = Icges(6,5) * t222;
t298 = Icges(6,1) * t221 - t215;
t521 = -t298 - t299;
t388 = t221 * t245;
t196 = Icges(6,5) * t388;
t386 = t222 * t245;
t398 = Icges(6,6) * t246;
t101 = -Icges(6,3) * t386 + t196 + t398;
t107 = Icges(5,6) * t246 + t245 * t297;
t520 = -t101 + t107;
t104 = -Icges(5,3) * t245 + t246 * t295;
t400 = Icges(6,2) * t245;
t106 = Icges(6,4) * t387 - Icges(6,6) * t385 - t400;
t518 = -t104 - t106;
t517 = t511 * t246;
t166 = -Icges(5,2) * t221 + t404;
t284 = t222 * t166 + t221 * t170;
t447 = Icges(6,3) * t221 + t215;
t510 = t221 * t168 - t222 * t447 + t284;
t472 = -t101 * t385 - t524 * t245;
t109 = Icges(6,4) * t246 + t245 * t298;
t198 = Icges(5,4) * t386;
t402 = Icges(5,5) * t246;
t111 = Icges(5,1) * t388 + t198 + t402;
t516 = -t109 - t111;
t470 = t523 * t245;
t514 = t522 * qJD(4);
t513 = t521 * qJD(4);
t512 = -t166 + t447;
t509 = -t295 - t296;
t464 = t523 * t246;
t508 = t495 * t246;
t421 = t246 * t105 + t109 * t388;
t35 = -t101 * t386 + t421;
t37 = t246 * t103 + t107 * t386 + t111 * t388;
t507 = t35 + t37;
t476 = t518 * t246 + t519 * t386 - t388 * t525;
t289 = t107 * t222 + t111 * t221;
t506 = -t109 * t387 - t289 * t246 - t472;
t474 = t518 * t245 + t508;
t505 = t510 * t245 + t464;
t504 = t168 * t387 + t246 * t284 - t385 * t447 - t470;
t132 = t170 * t245;
t503 = (Icges(6,1) * t386 + t132 + t196 - t520) * t246 + (-t517 - t519) * t245;
t467 = -t109 * t221 - t289;
t129 = t166 * t246;
t343 = qJD(4) * t246;
t502 = qJD(4) * t129 - t447 * t343 + (t245 * t293 - t107 + t398) * qJD(1);
t122 = t447 * t245;
t344 = qJD(4) * t245;
t501 = qJD(4) * t122 - t166 * t344 + (t246 * t293 + t527) * qJD(1);
t500 = -t517 * qJD(4) + (t245 * t299 + t109 + t402) * qJD(1);
t499 = qJD(4) * t132 + t168 * t344 + (t246 * t298 + t526) * qJD(1);
t473 = -t519 * t221 - t222 * t525;
t471 = t520 * t221 + t516 * t222;
t407 = rSges(6,3) + qJ(5);
t498 = t407 * t221;
t497 = t510 * qJD(1) + qJD(4) * t509;
t496 = (Icges(5,2) * t388 + t122 - t198 + t516) * t246 + (-Icges(6,3) * t387 + t129 - t197 + t525) * t245;
t392 = t101 * t222;
t494 = t392 + t467;
t493 = -t512 - t521;
t492 = -t511 - t522;
t491 = t514 * t222 + t513 * t221 + (t221 * t512 + t222 * t511) * qJD(4) - t523 * qJD(1);
t483 = rSges(6,1) + pkin(4);
t490 = t505 * qJD(1);
t489 = (t245 * t474 + t246 * t506) * qJD(4);
t488 = (t476 * t245 + t246 * t507) * qJD(4);
t487 = t524 * qJD(1);
t486 = t504 * qJD(1);
t485 = (t221 * t493 + t222 * t492) * qJD(1);
t484 = -t221 * t503 + t222 * t496;
t433 = t245 / 0.2e1;
t432 = -t246 / 0.2e1;
t482 = t488 + t490;
t481 = -t486 + t489;
t480 = t245 * t497 - t246 * t491;
t479 = t245 * t491 + t246 * t497;
t478 = qJD(4) * t494 + t221 * t501 + t222 * t499;
t477 = t495 * qJD(4) + t221 * t502 + t500 * t222;
t412 = rSges(6,3) * t222;
t418 = rSges(6,1) * t221;
t302 = -t412 + t418;
t393 = qJ(5) * t222;
t427 = pkin(4) * t221;
t360 = t393 - t427 - t302;
t174 = pkin(4) * t222 + qJ(5) * t221;
t175 = rSges(6,1) * t222 + rSges(6,3) * t221;
t359 = t174 + t175;
t348 = qJD(1) * t104;
t469 = t348 + t470 * qJD(4) + (t246 * t296 - t400 - t494) * qJD(1);
t468 = -qJD(1) * t495 - qJD(4) * t464 + t487;
t345 = qJD(1) * t246;
t466 = t221 * t345 + t222 * t344;
t465 = t509 * qJD(1);
t462 = qJD(4) * t471 - t221 * t499 + t222 * t501 + t487;
t461 = qJD(1) * t106 + t473 * qJD(4) + t500 * t221 - t222 * t502 + t348;
t459 = 0.2e1 * qJD(4);
t212 = qJD(5) * t221;
t357 = -t245 * rSges(6,2) - rSges(6,3) * t385;
t449 = -qJ(5) * t385 + t357;
t369 = t387 * t483 + t449;
t239 = t246 * rSges(6,2);
t448 = t388 * t483 + t239;
t370 = t386 * t407 - t448;
t31 = t212 + (t245 * t370 - t246 * t369) * qJD(4);
t458 = qJD(4) * t31;
t342 = qJD(5) * t222;
t273 = qJD(4) * t359 - t342;
t457 = t246 * t273;
t191 = t246 * pkin(1) + t245 * qJ(2);
t395 = qJ(3) * t246;
t242 = sin(pkin(7));
t384 = t242 * t245;
t416 = rSges(4,2) * cos(pkin(7));
t445 = -rSges(4,1) * t384 - t246 * rSges(4,3) - t245 * t416;
t451 = t191 + t395 - t445;
t315 = -rSges(3,2) * t246 + t245 * rSges(3,3);
t450 = t191 + t315;
t211 = pkin(3) * t384;
t330 = t211 + t191;
t446 = t344 * t498 + t466 * t483;
t337 = qJD(1) * qJD(3);
t338 = qJD(1) * qJD(2);
t346 = qJD(1) * t245;
t218 = qJ(2) * t345;
t224 = qJD(2) * t245;
t353 = t218 + t224;
t365 = qJD(1) * (-pkin(1) * t346 + t353) + t245 * t338;
t394 = qJ(3) * qJD(1) ^ 2;
t274 = -t245 * t394 + 0.2e1 * t246 * t337 + t365;
t244 = -pkin(6) - qJ(3);
t326 = t242 * t345;
t355 = pkin(3) * t326 + t244 * t346;
t266 = qJD(1) * (qJ(3) * t346 + t355) + t274;
t371 = qJD(4) * t360 + t212;
t306 = t212 + t371;
t341 = qJD(5) * t245;
t422 = -(-qJ(5) * t345 - t341) * t222 - t357 * qJD(1) - t446;
t10 = -t306 * t343 + (t245 * t273 - t422) * qJD(1) + t266;
t217 = t246 * t338;
t282 = -0.2e1 * t245 * t337 - t246 * t394 + t217;
t225 = qJD(2) * t246;
t156 = qJD(1) * t191 - t225;
t208 = t244 * t345;
t368 = t208 - (t211 - t395) * qJD(1) - t156;
t140 = t175 * t246;
t423 = (pkin(4) * t346 - qJ(5) * t343) * t221 + (-qJ(5) * t346 + (-pkin(4) * qJD(4) + qJD(5)) * t246) * t222 - qJD(4) * t140 + (t245 * t302 + t239) * qJD(1);
t11 = t306 * t344 + (t368 - t423 + t457) * qJD(1) + t282;
t444 = t10 * t432 + t11 * t433;
t434 = t245 ^ 2;
t430 = rSges(3,2) - pkin(1);
t429 = -rSges(6,2) - pkin(1);
t428 = pkin(3) * t242;
t426 = -qJD(1) / 0.2e1;
t424 = -pkin(1) + t244;
t419 = rSges(5,1) * t222;
t415 = rSges(5,2) * t221;
t414 = rSges(5,2) * t222;
t413 = rSges(3,3) * t246;
t176 = -t415 + t419;
t145 = t176 * t344;
t303 = rSges(5,1) * t221 + t414;
t158 = t303 * qJD(4);
t334 = qJD(4) * t415;
t280 = -rSges(5,3) * qJD(1) - t334;
t331 = rSges(5,1) * t466 + t345 * t414;
t80 = t245 * t280 + t331;
t24 = t158 * t343 + (t80 + t145) * qJD(1) + t266;
t409 = t24 * t246;
t237 = t246 * rSges(5,3);
t325 = t176 * t343;
t141 = t176 * t246;
t78 = -qJD(4) * t141 + (t245 * t303 + t237) * qJD(1);
t25 = -t158 * t344 + (-t78 + t325 + t368) * qJD(1) + t282;
t408 = t25 * t245;
t396 = qJ(3) * t245;
t383 = t242 * t246;
t380 = qJ(3) + t244;
t367 = t359 * t245;
t366 = -t174 * t246 - t140;
t335 = t246 * t416;
t356 = rSges(4,1) * t326 + qJD(1) * t335;
t354 = t208 + t225;
t352 = rSges(3,2) * t346 + rSges(3,3) * t345;
t223 = qJD(3) * t246;
t351 = t223 + t224;
t227 = t246 * qJ(2);
t188 = pkin(1) * t245 - t227;
t178 = qJD(1) * t188;
t350 = t224 - t178;
t336 = -rSges(4,3) - pkin(1) - qJ(3);
t115 = rSges(5,1) * t388 + rSges(5,2) * t386 + t237;
t329 = t218 + t351;
t328 = t223 + t350;
t320 = -t344 / 0.2e1;
t318 = -t343 / 0.2e1;
t314 = -t188 - t396;
t313 = t106 - t392;
t312 = qJD(3) * t245 - t225;
t146 = pkin(3) * t383 + t245 * t380;
t308 = t146 + t314;
t307 = -t246 * t380 + t330 + t395;
t304 = rSges(4,1) * t242 + t416;
t234 = t245 * rSges(5,3);
t117 = t246 * t303 - t234;
t43 = t145 + (t117 + t308) * qJD(1) + t351;
t44 = -t325 + (t115 + t307) * qJD(1) + t312;
t301 = t245 * t43 - t246 * t44;
t300 = t329 + t355;
t57 = (-t115 * t245 - t117 * t246) * qJD(4);
t272 = t418 + t427 + t428;
t271 = t303 + t428;
t265 = -t222 * t341 + t344 * t359 + t351;
t189 = rSges(3,2) * t245 + t413;
t137 = t176 * t245;
t121 = -rSges(4,3) * t245 + t246 * t304;
t120 = qJD(1) * t146;
t93 = qJD(1) * t450 - t225;
t92 = t224 + (-t188 + t189) * qJD(1);
t82 = t217 + (-qJD(1) * t315 - t156) * qJD(1);
t81 = qJD(1) * t352 + t365;
t76 = qJD(1) * t451 + t312;
t75 = (t121 + t314) * qJD(1) + t351;
t48 = (qJD(1) * t445 - t156) * qJD(1) + t282;
t47 = qJD(1) * (-rSges(4,3) * t346 + t356) + t274;
t30 = -t457 + (t307 - t370) * qJD(1) + t312;
t29 = (t308 + t369) * qJD(1) + t265;
t1 = (t342 + t423 * t246 + t422 * t245 + (t245 * t369 + t246 * t370) * qJD(1)) * qJD(4);
t2 = [(((t37 + t421 + t474 - t508) * t246 + ((t104 + t313 + t467) * t246 - t472 - t506 + t476) * t245) * qJD(4) + t490) * t320 + (-t510 * qJD(4) - t514 * t221 + t513 * t222) * qJD(1) + (-(t120 - t178 - t29 + (t369 - t396) * qJD(1) + t265) * t30 + t11 * (t227 + t449) + t29 * t354 + t10 * (t330 + t448) + t30 * (t300 + t446) + (t11 * t424 - t29 * qJD(3) + (-t30 * qJD(5) - t10 * t407) * t222) * t245 + (-t10 * t244 + t11 * t272 + (-t342 + (t222 * t483 + t498) * qJD(4)) * t29) * t246 + ((-t222 * t30 * t407 + t29 * t429) * t246 + (t29 * (-qJ(2) - t272 + t393 + t412) + t30 * t429) * t245) * qJD(1)) * m(6) + ((-t244 * t246 + t115 + t330) * t24 + (t424 * t245 + t271 * t246 + t227 - t234) * t25 + (t354 + (-pkin(1) * qJD(1) + qJD(4) * t419 + t280) * t246 + (-qJD(3) + (-qJ(2) - t271) * qJD(1)) * t245) * t43 + (t300 + t331 + (-t334 + (-rSges(5,3) - pkin(1)) * qJD(1)) * t245 - t120 - t145 + t43 - (t117 - t396) * qJD(1) - t328) * t44) * m(5) + (-(-t75 + (t121 - t396) * qJD(1) + t328) * t76 + t48 * (rSges(4,1) * t383 + t227 + t335) + t75 * t225 + t47 * t451 + t76 * (t329 + t356) + (-t75 * qJD(3) + t48 * t336) * t245 + (t75 * t336 * t246 + (t75 * (-qJ(2) - t304) + t76 * t336) * t245) * qJD(1)) * m(4) + (-(qJD(1) * t189 + t350 - t92) * t93 + t82 * (t245 * t430 + t227 + t413) + t92 * t225 + t81 * t450 + t93 * (t352 + t353) + (t92 * t430 * t246 + (t92 * (-rSges(3,3) - qJ(2)) - t93 * pkin(1)) * t245) * qJD(1)) * m(3) + ((t434 * t104 + (t245 * t313 - t35 + t421) * t245 + ((-t467 - t518) * t246 + t472 + t476) * t246) * qJD(4) + t481 + t486) * t318 + (t477 + t480 + t482) * t344 / 0.2e1 + (qJD(1) * t473 + t478 + t479) * t343 / 0.2e1 + (t504 * t246 + (-t471 + t505) * t245) * qJD(4) * t426; 0.2e1 * t444 * m(6) + 0.2e1 * (-t409 / 0.2e1 + t408 / 0.2e1) * m(5) + 0.2e1 * (t432 * t47 + t433 * t48) * m(4) + 0.2e1 * (t432 * t81 + t433 * t82) * m(3); m(4) * (t245 * t47 + t246 * t48) + m(5) * (t24 * t245 + t246 * t25) + m(6) * (t10 * t245 + t11 * t246); ((t496 * t221 + t222 * t503) * qJD(4) + (t221 * t492 - t222 * t493) * qJD(1)) * t426 + (t478 * t246 + t477 * t245 + (t471 * t245 + t473 * t246) * qJD(1)) * qJD(1) / 0.2e1 + ((-t344 * t464 + t465) * t245 + ((t245 * t470 + t484) * qJD(4) + t485) * t246) * t320 + ((t343 * t470 + t465) * t246 + ((-t246 * t464 - t484) * qJD(4) - t485) * t245) * t318 + ((-t10 * t359 - t30 * t371 - t1 * t369 + t31 * t423 + (t29 * t359 + t31 * t370) * qJD(1)) * t246 + (t11 * t359 + t29 * t371 + t1 * t370 + t31 * t422 + (t30 * t359 + t31 * t369) * qJD(1)) * t245 - (t222 * t31 + (t245 * t29 - t246 * t30) * t221) * qJD(5) - (-t29 * t366 + t30 * t367) * qJD(1) - ((-t30 * t360 + t31 * t366) * t246 + (t29 * t360 - t31 * t367) * t245) * qJD(4)) * m(6) + (-(t137 * t44 + t141 * t43) * qJD(1) - (t57 * (-t137 * t245 - t141 * t246) - t301 * t303) * qJD(4) + 0.2e1 * t57 * (-t245 * t80 + t246 * t78 + (-t115 * t246 + t117 * t245) * qJD(1)) - t301 * t158 + (-t409 + t408 + (t245 * t44 + t246 * t43) * qJD(1)) * t176) * m(5) + (t480 * qJD(1) + ((t474 * qJD(1) + t462 * t246) * t246 + (t468 * t245 - t506 * qJD(1) + (-t461 + t469) * t246) * t245) * t459) * t433 + (t479 * qJD(1) + ((t476 * qJD(1) + t469 * t246) * t246 + (t461 * t245 - t507 * qJD(1) + (-t462 + t468) * t246) * t245) * t459) * t246 / 0.2e1 - (t482 + t488) * t346 / 0.2e1 + (t481 + t489) * t345 / 0.2e1; (t1 * t221 + 0.2e1 * (t458 / 0.2e1 - (t246 ^ 2 + t434) * t458 / 0.2e1 - t444) * t222) * m(6);];
tauc = t2(:);
