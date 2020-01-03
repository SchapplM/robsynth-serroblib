% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:29
% EndTime: 2020-01-03 11:31:00
% DurationCPUTime: 11.56s
% Computational Cost: add. (40177->479), mult. (39793->696), div. (0->0), fcn. (43098->10), ass. (0->288)
t367 = pkin(9) + qJ(4);
t357 = sin(t367);
t369 = sin(pkin(8));
t371 = cos(pkin(8));
t358 = cos(t367);
t468 = Icges(5,4) * t358;
t292 = -Icges(5,6) * t371 + (-Icges(5,2) * t357 + t468) * t369;
t469 = Icges(5,4) * t357;
t293 = -Icges(5,5) * t371 + (Icges(5,1) * t358 - t469) * t369;
t318 = (-Icges(5,2) * t358 - t469) * t369;
t319 = (-Icges(5,1) * t357 - t468) * t369;
t317 = (-Icges(5,5) * t357 - Icges(5,6) * t358) * t369;
t451 = t371 * t317;
t563 = t451 / 0.2e1 - (-(t293 / 0.2e1 + t318 / 0.2e1) * t357 + (t319 / 0.2e1 - t292 / 0.2e1) * t358) * t369;
t370 = cos(pkin(9));
t356 = t370 * pkin(3) + pkin(2);
t337 = pkin(4) * t358 + t356;
t373 = cos(qJ(1));
t362 = t373 * qJ(2);
t372 = sin(qJ(1));
t359 = qJ(5) + t367;
t355 = cos(t359);
t354 = sin(t359);
t449 = t372 * t354;
t305 = -t355 * t373 - t371 * t449;
t444 = t373 * t354;
t448 = t372 * t355;
t306 = t371 * t448 - t444;
t390 = rSges(6,1) * t306 + rSges(6,2) * t305;
t479 = pkin(6) + qJ(3);
t365 = -pkin(7) - t479;
t480 = pkin(4) * t357;
t368 = sin(pkin(9));
t481 = pkin(3) * t368;
t392 = t480 + t481;
t454 = t369 * t372;
t416 = -t365 * t454 - t373 * t392;
t167 = -t362 + (rSges(6,3) * t369 + t337 * t371 + pkin(1)) * t372 + t390 + t416;
t264 = t305 * rSges(6,1) - t306 * rSges(6,2);
t307 = t371 * t444 - t448;
t450 = t371 * t373;
t308 = t355 * t450 + t449;
t265 = t307 * rSges(6,1) + t308 * rSges(6,2);
t453 = t369 * t373;
t241 = t308 * rSges(6,1) - t307 * rSges(6,2) + rSges(6,3) * t453;
t413 = t373 * pkin(1) + t372 * qJ(2);
t417 = -t337 * t450 - t372 * t392;
t531 = -t365 * t453 + t241 + t413 - t417;
t562 = m(6) * (t167 * t264 - t265 * t531);
t396 = t369 * t479;
t445 = t372 * t368;
t394 = pkin(3) * t445 + t373 * t396;
t455 = t356 * t371;
t208 = (t365 * t369 + t455) * t373 + t394 + t417;
t228 = t371 * t241;
t290 = -t371 * rSges(6,3) + (rSges(6,1) * t355 - rSges(6,2) * t354) * t369;
t414 = -t356 + t337;
t422 = (t365 + t479) * t371 + t414 * t369 + t290;
t393 = t422 * t453;
t106 = t371 * t208 - t228 - t393;
t107 = (t208 - t241) * t371 - t393;
t442 = t106 - t107;
t561 = m(6) * t442;
t232 = Icges(6,5) * t308 - Icges(6,6) * t307 + Icges(6,3) * t453;
t560 = t232 * t453;
t231 = Icges(6,5) * t306 + Icges(6,6) * t305 + Icges(6,3) * t454;
t559 = t231 * t453;
t558 = -t290 * t453 - t228;
t557 = -t369 / 0.2e1;
t466 = Icges(6,4) * t354;
t289 = -Icges(6,5) * t371 + (Icges(6,1) * t355 - t466) * t369;
t420 = t289 + (-Icges(6,2) * t355 - t466) * t369;
t556 = t354 * t420;
t443 = t373 * t357;
t446 = t372 * t358;
t322 = t371 * t443 - t446;
t447 = t372 * t357;
t323 = t358 * t450 + t447;
t247 = Icges(5,5) * t323 - Icges(5,6) * t322 + Icges(5,3) * t453;
t555 = t247 * t453;
t256 = t323 * rSges(5,1) - t322 * rSges(5,2) + rSges(5,3) * t453;
t294 = -t371 * rSges(5,3) + (rSges(5,1) * t358 - rSges(5,2) * t357) * t369;
t554 = -t371 * t256 - t294 * t453;
t297 = Icges(6,4) * t308;
t236 = Icges(6,2) * t307 - Icges(6,6) * t453 - t297;
t296 = Icges(6,4) * t307;
t238 = Icges(6,1) * t308 + Icges(6,5) * t453 - t296;
t553 = -t305 * t236 + t306 * t238;
t312 = Icges(5,4) * t323;
t251 = Icges(5,2) * t322 - Icges(5,6) * t453 - t312;
t311 = Icges(5,4) * t322;
t253 = Icges(5,1) * t323 + Icges(5,5) * t453 - t311;
t320 = -t358 * t373 - t371 * t447;
t321 = t371 * t446 - t443;
t551 = -t320 * t251 + t321 * t253;
t520 = m(5) / 0.2e1;
t518 = m(6) / 0.2e1;
t538 = t369 / 0.2e1;
t460 = t232 * t372;
t461 = t231 * t373;
t528 = t373 ^ 2;
t546 = t369 * t528;
t97 = -t232 * t454 - t553;
t549 = (t232 * t546 - t560 * t373 + (t97 + t553 + (t460 + t461) * t369 - t559) * t372) * t369;
t548 = rSges(4,3) + qJ(3);
t246 = Icges(5,5) * t321 + Icges(5,6) * t320 + Icges(5,3) * t454;
t543 = t246 * t453;
t104 = t106 * t372;
t415 = -t372 * t396 + t373 * t481;
t207 = t371 * t372 * t414 + t415 + t416;
t240 = rSges(6,3) * t454 + t390;
t105 = (t207 + t240) * t371 + t422 * t454;
t391 = rSges(5,1) * t321 + rSges(5,2) * t320;
t255 = rSges(5,3) * t454 + t391;
t192 = t255 * t371 + t294 * t454;
t471 = (t192 * t373 + t554 * t372) * t520 + (t105 * t373 + t104) * t518;
t477 = (-t107 * t372 + t104) * t518;
t541 = t471 - t477;
t126 = -t192 * t454 + t453 * t554;
t64 = -t105 * t454 + t106 * t453;
t472 = t126 * t520 + t64 * t518;
t478 = ((t105 * t372 - t107 * t373) * t369 + t64) * t518 + ((t192 * t372 - t373 * t554) * t369 + t126) * t520;
t540 = t472 - t478;
t302 = (-Icges(6,5) * t354 - Icges(6,6) * t355) * t369;
t452 = t371 * t302;
t539 = -t452 / 0.2e1 + t556 * t557;
t497 = -t371 / 0.2e1;
t314 = t322 * pkin(4);
t163 = t264 * t453 + t265 * t454;
t309 = (-rSges(6,1) * t354 - rSges(6,2) * t355) * t369;
t195 = -t371 * t264 - t309 * t454;
t245 = t371 * t265;
t196 = -t309 * t453 + t245;
t465 = Icges(6,4) * t355;
t288 = -Icges(6,6) * t371 + (-Icges(6,2) * t354 + t465) * t369;
t304 = (-Icges(6,1) * t354 - t465) * t369;
t421 = t288 - t304;
t102 = t302 * t454 + t305 * t420 - t306 * t421;
t103 = -t302 * t453 + t307 * t420 + t308 * t421;
t258 = Icges(6,5) * t305 - Icges(6,6) * t306;
t259 = Icges(6,5) * t307 + Icges(6,6) * t308;
t400 = -t453 / 0.2e1;
t401 = t454 / 0.2e1;
t428 = Icges(6,2) * t308 - t238 + t296;
t295 = Icges(6,4) * t305;
t237 = Icges(6,1) * t306 + Icges(6,5) * t454 + t295;
t429 = -Icges(6,2) * t306 + t237 + t295;
t430 = -Icges(6,1) * t307 + t236 - t297;
t467 = Icges(6,4) * t306;
t234 = Icges(6,2) * t305 + Icges(6,6) * t454 + t467;
t431 = -Icges(6,1) * t305 + t234 + t467;
t137 = -t452 + (-t355 * t421 - t556) * t369;
t464 = t137 * t371;
t80 = -t258 * t371 + (-t354 * t429 - t355 * t431) * t369;
t81 = -t259 * t371 + (-t354 * t428 - t355 * t430) * t369;
t411 = (-t464 + (t372 * t80 - t373 * t81) * t369) * t497 + (-t102 * t371 + (t258 * t454 + t305 * t429 - t306 * t431) * t454 - (t259 * t454 + t305 * t428 - t306 * t430) * t453) * t401 + (-t103 * t371 + (-t258 * t453 + t307 * t429 + t308 * t431) * t454 - (-t259 * t453 + t307 * t428 + t308 * t430) * t453) * t400;
t150 = t240 * t453 - t241 * t454;
t86 = (t207 * t373 + t208 * t372) * t369 + t150;
t6 = t411 + m(6) * (-t105 * t195 + t106 * t196 + t86 * t163);
t536 = t6 * qJD(5);
t272 = rSges(5,1) * t320 - rSges(5,2) * t321;
t273 = rSges(5,1) * t322 + rSges(5,2) * t323;
t173 = (t272 * t373 + t273 * t372) * t369;
t530 = t356 * t450 + t256 + t394 + t413;
t217 = -t265 - t314;
t313 = t320 * pkin(4);
t423 = -t264 - t313;
t439 = (-t217 * t373 + t372 * t423) * t518 + (-t372 * t272 + t273 * t373) * t520;
t440 = m(6) * (t217 * t372 + t373 * t423) * t538 - m(5) * t173 / 0.2e1;
t336 = (t372 ^ 2 + t528) * t369;
t527 = 0.2e1 * t336;
t526 = 2 * qJD(1);
t525 = 4 * qJD(1);
t524 = 2 * qJD(4);
t523 = 4 * qJD(4);
t522 = m(4) / 0.2e1;
t519 = m(5) / 0.4e1;
t517 = m(6) / 0.4e1;
t190 = -t362 + (rSges(5,3) * t369 + pkin(1) + t455) * t372 + t391 - t415;
t513 = m(5) * (t190 * t272 - t273 * t530);
t182 = t371 * t240 + t290 * t454;
t512 = t182 * t561;
t511 = t105 * t561;
t441 = t195 * t167 + t196 * t531;
t508 = m(6) * (-t105 * t264 - t106 * t265 + t441);
t506 = m(6) * (t182 * t423 + t217 * t558 + t441);
t115 = -t182 * t454 + t453 * t558;
t505 = m(6) * ((t182 * t372 - t373 * t558) * t369 + t115);
t501 = m(6) * (-t167 * t423 + t217 * t531);
t92 = t167 * t454 + t453 * t531;
t500 = m(6) * t92;
t499 = m(6) * (-t167 * t373 + t531 * t372);
t496 = m(3) * (-(-rSges(3,2) * t454 - t373 * rSges(3,3) + t372 * pkin(1) - t362) * t373 + (-rSges(3,2) * t453 + t372 * rSges(3,3) + t413) * t372);
t219 = -t362 + (-rSges(4,1) * t368 - rSges(4,2) * t370) * t373 + (pkin(1) + t548 * t369 + (rSges(4,1) * t370 - rSges(4,2) * t368 + pkin(2)) * t371) * t372;
t376 = rSges(4,1) * (t370 * t450 + t445) - rSges(4,2) * (t368 * t450 - t372 * t370) + pkin(2) * t450 + t413 + t548 * t453;
t142 = t219 * t454 + t376 * t453;
t495 = m(4) * t142;
t494 = m(4) * (-t219 * t373 + t376 * t372);
t125 = t190 * t454 + t453 * t530;
t493 = m(5) * t125;
t491 = m(5) * (-t190 * t373 + t530 * t372);
t487 = m(6) * t115;
t486 = m(6) * (t182 * t373 + t558 * t372);
t483 = m(6) * t163;
t482 = m(6) * (-t372 * t264 + t265 * t373);
t476 = m(6) * qJD(4);
t475 = m(6) * qJD(5);
t470 = Icges(5,4) * t321;
t463 = ((-Icges(6,3) * t371 + (Icges(6,5) * t355 - Icges(6,6) * t354) * t369) * t454 + t288 * t305 + t289 * t306) * t371;
t459 = t246 * t373;
t458 = t247 * t372;
t456 = t355 * t369;
t249 = Icges(5,2) * t320 + Icges(5,6) * t454 + t470;
t427 = -Icges(5,1) * t320 + t249 + t470;
t426 = -Icges(5,1) * t322 + t251 - t312;
t310 = Icges(5,4) * t320;
t252 = Icges(5,1) * t321 + Icges(5,5) * t454 + t310;
t425 = -Icges(5,2) * t321 + t252 + t310;
t424 = Icges(5,2) * t323 - t253 + t311;
t419 = t292 - t319;
t418 = t293 + t318;
t171 = (t522 + t520 + t518) * t527;
t412 = t171 * qJD(1);
t110 = -t247 * t454 - t551;
t410 = (t247 * t546 - t555 * t373 + (t110 + (t458 + t459) * t369 - t543 + t551) * t372) * t538 + (t497 + t371 / 0.2e1) * ((-Icges(5,3) * t371 + (Icges(5,5) * t358 - Icges(5,6) * t357) * t369) * t453 - t322 * t292 + t323 * t293);
t109 = t246 * t454 + t320 * t249 + t321 * t252;
t409 = (t109 * t372 - t110 * t373) * t538 + ((t109 + t555) * t372 + ((-t458 + t459) * t369 - t110 - t543) * t373) * t557;
t96 = t231 * t454 + t305 * t234 + t306 * t237;
t404 = -t456 / 0.2e1;
t403 = t456 / 0.2e1;
t399 = t453 / 0.2e1;
t21 = -t463 + ((t96 + t560) * t372 + ((-t460 + t461) * t369 - t97 - t559) * t373) * t369;
t51 = -t463 + (t372 * t96 - t373 * t97) * t369;
t395 = t21 * t400 + t51 * t399 + t549 * t401;
t389 = t288 * t404 + t304 * t403 + t539;
t388 = t512 / 0.2e1 + t395;
t384 = (t480 * t369 - t309) * t369;
t379 = t21 * t399 - t549 * t454 / 0.2e1 + (t80 + t102) * t401 + (t51 + t81 + t103) * t400;
t377 = t288 * t403 + t304 * t404 - t539;
t375 = t379 - t464;
t324 = (-rSges(5,1) * t357 - rSges(5,2) * t358) * t369;
t267 = Icges(5,5) * t322 + Icges(5,6) * t323;
t266 = Icges(5,5) * t320 - Icges(5,6) * t321;
t204 = t371 * t273 - t324 * t453;
t203 = -t272 * t371 - t324 * t454;
t172 = (m(4) / 0.4e1 + t519 + t517) * t527 - (m(4) + m(5) + m(6)) * t336 / 0.2e1;
t165 = t482 / 0.2e1;
t159 = -t483 / 0.2e1;
t156 = t371 * t314 + t373 * t384 + t245;
t155 = t371 * t423 + t372 * t384;
t141 = -t451 + (-t357 * t418 - t358 * t419) * t369;
t136 = (t313 * t373 + t314 * t372) * t369 + t163;
t135 = -t195 * t372 - t196 * t373;
t132 = t135 * t475;
t118 = t486 / 0.2e1;
t117 = -t317 * t453 + t322 * t418 + t323 * t419;
t116 = t317 * t454 + t320 * t418 - t321 * t419;
t108 = t487 / 0.2e1;
t89 = -t267 * t371 + (-t357 * t424 - t358 * t426) * t369;
t88 = -t266 * t371 + (-t357 * t425 - t358 * t427) * t369;
t85 = -t163 * t371 + (-t195 * t373 + t196 * t372) * t369;
t84 = t85 * t475;
t62 = t389 + t562;
t58 = t505 / 0.2e1;
t54 = t493 + t495 + t500;
t45 = t491 + t494 + t496 + t499;
t44 = t118 - t482 / 0.2e1;
t43 = t165 + t118;
t42 = t165 - t486 / 0.2e1;
t40 = t506 / 0.2e1;
t38 = t108 + t58 + t483 / 0.2e1;
t37 = t159 + t108 - t505 / 0.2e1;
t36 = t159 + t58 - t487 / 0.2e1;
t32 = t508 / 0.2e1;
t31 = t513 + t501 + t389 - t563;
t14 = t471 + t477 - t439;
t13 = t439 + t541;
t12 = t439 - t541;
t11 = t472 + t478 - t440;
t10 = t440 + t540;
t9 = t440 - t540;
t8 = m(6) * (t150 * t163 - t182 * t195 + t196 * t558) + t411;
t7 = t8 * qJD(5);
t4 = t40 - t508 / 0.2e1 + t388;
t3 = t32 - t506 / 0.2e1 + t388;
t2 = t32 + t40 - t512 / 0.2e1 + t375;
t1 = t511 + (t372 * t410 + t373 * t409) * t369 + t395;
t5 = [t45 * qJD(2) + t54 * qJD(3) + t31 * qJD(4) + t62 * qJD(5), qJD(1) * t45 + qJD(3) * t172 + qJD(4) * t13 + qJD(5) * t43, qJD(1) * t54 + qJD(2) * t172 + qJD(4) * t10 + qJD(5) * t37, t31 * qJD(1) + t13 * qJD(2) + t10 * qJD(3) + t2 * qJD(5) - t511 * t523 / 0.4e1 + ((t190 * t203 - t192 * t272 + t204 * t530 - t273 * t554) * t520 + (t105 * t423 + t106 * t217 + t155 * t167 + t156 * t531) * t518) * t524 + (t379 + (-t141 - t137) * t371 + ((-t89 / 0.2e1 - t117 / 0.2e1 - t409) * t373 + (t116 / 0.2e1 + t88 / 0.2e1 - t410) * t372) * t369) * qJD(4), t62 * qJD(1) + t43 * qJD(2) + t37 * qJD(3) + t2 * qJD(4) + ((-t182 * t264 - t265 * t558 + t441) * m(6) + t375) * qJD(5); -t171 * qJD(3) + t12 * qJD(4) + t42 * qJD(5) + (-t496 / 0.4e1 - t494 / 0.4e1 - t491 / 0.4e1 - t499 / 0.4e1) * t525, 0, -t412, t12 * qJD(1) + ((-t203 * t372 - t204 * t373) * t520 + (-t155 * t372 - t156 * t373) * t518) * t524 + t132, t42 * qJD(1) + t135 * t476 + t132; t171 * qJD(2) + t9 * qJD(4) + t36 * qJD(5) + (-t495 / 0.4e1 - t493 / 0.4e1 - t500 / 0.4e1) * t525 + (((-t219 * t372 - t373 * t376) * t369 + t142) * t522 + ((-t190 * t372 - t373 * t530) * t369 + t125) * t520 + ((-t167 * t372 - t373 * t531) * t369 + t92) * t518) * t526, t412, 0, t9 * qJD(1) + ((-t173 * t371 + (-t203 * t373 + t204 * t372) * t369) * t520 + (-t136 * t371 + (-t155 * t373 + t156 * t372) * t369) * t518) * t524 + t84, t36 * qJD(1) + t85 * t476 + t84; t14 * qJD(2) + t11 * qJD(3) + t1 * qJD(4) + t3 * qJD(5) + (-t513 / 0.4e1 - t501 / 0.4e1) * t525 - t442 * t167 * t518 * t526 + (t377 + t563) * qJD(1), t14 * qJD(1), t11 * qJD(1), t1 * qJD(1) + ((-t141 * t371 + (t372 * t88 - t373 * t89) * t369) * t497 + (-t116 * t371 + (t266 * t454 + t320 * t425 - t321 * t427) * t454 - (t267 * t454 + t320 * t424 - t321 * t426) * t453) * t401 + (-t117 * t371 + (-t266 * t453 + t322 * t425 + t323 * t427) * t454 - (-t267 * t453 + t322 * t424 + t323 * t426) * t453) * t400 + t411) * qJD(4) + t536 + (((t255 * t373 - t256 * t372) * t369 * t173 - t192 * t203 + t204 * t554) * t519 + (-t105 * t155 + t106 * t156 + t136 * t86) * t517) * t523, t3 * qJD(1) + t6 * qJD(4) + t536; (t377 - t562) * qJD(1) + t44 * qJD(2) + t38 * qJD(3) + t4 * qJD(4) + t395 * qJD(5), t44 * qJD(1), t38 * qJD(1), t4 * qJD(1) + ((t136 * t150 - t155 * t182 + t156 * t558) * m(6) + t411) * qJD(4) + t7, qJD(1) * t395 + qJD(4) * t8 + t7;];
Cq = t5;
