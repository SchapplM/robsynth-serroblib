% Calculate time derivative of joint inertia matrix for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:18:04
% EndTime: 2019-03-09 03:18:32
% DurationCPUTime: 20.24s
% Computational Cost: add. (21842->890), mult. (31387->1235), div. (0->0), fcn. (29230->8), ass. (0->411)
t296 = pkin(9) + qJ(3);
t289 = cos(t296);
t288 = sin(t296);
t478 = Icges(5,6) * t288;
t488 = Icges(4,4) * t288;
t556 = -t478 - t488 + (-Icges(4,2) - Icges(5,3)) * t289;
t477 = Icges(5,6) * t289;
t487 = Icges(4,4) * t289;
t555 = -t477 - t487 + (-Icges(4,1) - Icges(5,2)) * t288;
t303 = sin(qJ(5));
t305 = cos(qJ(5));
t363 = Icges(7,5) * t303 + Icges(7,6) * t305;
t193 = Icges(7,3) * t288 - t289 * t363;
t364 = Icges(6,5) * t303 + Icges(6,6) * t305;
t194 = Icges(6,3) * t288 - t289 * t364;
t554 = t194 + t193;
t483 = Icges(7,4) * t303;
t366 = Icges(7,2) * t305 + t483;
t195 = Icges(7,6) * t288 - t289 * t366;
t485 = Icges(6,4) * t303;
t367 = Icges(6,2) * t305 + t485;
t196 = Icges(6,6) * t288 - t289 * t367;
t553 = t195 + t196;
t482 = Icges(7,4) * t305;
t371 = Icges(7,1) * t303 + t482;
t197 = Icges(7,5) * t288 - t289 * t371;
t484 = Icges(6,4) * t305;
t372 = Icges(6,1) * t303 + t484;
t198 = Icges(6,5) * t288 - t289 * t372;
t544 = -t198 - t197;
t552 = t556 * qJD(3);
t551 = t555 * qJD(3);
t301 = -qJ(6) - pkin(8);
t490 = rSges(7,3) - t301;
t304 = sin(qJ(1));
t306 = cos(qJ(1));
t370 = -Icges(4,2) * t288 + t487;
t204 = -Icges(4,6) * t306 + t304 * t370;
t360 = -Icges(5,3) * t288 + t477;
t523 = Icges(5,5) * t306 + t304 * t360;
t550 = -t204 - t523;
t205 = Icges(4,6) * t304 + t306 * t370;
t208 = Icges(5,5) * t304 - t306 * t360;
t549 = t205 - t208;
t374 = Icges(4,1) * t289 - t488;
t206 = -Icges(4,5) * t306 + t304 * t374;
t362 = Icges(5,2) * t289 - t478;
t522 = Icges(5,4) * t306 + t304 * t362;
t548 = t206 + t522;
t207 = Icges(4,5) * t304 + t306 * t374;
t210 = Icges(5,4) * t304 - t306 * t362;
t547 = t207 - t210;
t546 = (t544 * t303 - t553 * t305) * t289 + t554 * t288;
t436 = qJD(5) * t289;
t133 = (Icges(7,2) * t303 - t482) * t436 + (Icges(7,6) * t289 + t288 * t366) * qJD(3);
t134 = (Icges(6,2) * t303 - t484) * t436 + (Icges(6,6) * t289 + t288 * t367) * qJD(3);
t545 = -t134 - t133;
t131 = (-Icges(7,5) * t305 + Icges(7,6) * t303) * t436 + (Icges(7,3) * t289 + t288 * t363) * qJD(3);
t132 = (-Icges(6,5) * t305 + Icges(6,6) * t303) * t436 + (Icges(6,3) * t289 + t288 * t364) * qJD(3);
t435 = qJD(5) * t303;
t439 = qJD(3) * t305;
t441 = qJD(3) * t303;
t442 = qJD(3) * t289;
t543 = t553 * t289 * t435 + t554 * t442 + (t553 * t439 - t544 * t441 + t131 + t132) * t288;
t282 = pkin(5) * t305 + pkin(4);
t295 = t306 * pkin(4);
t502 = -pkin(8) - t301;
t505 = pkin(5) * t303;
t317 = t288 * t505 + t289 * t502;
t463 = t304 * t305;
t465 = t303 * t306;
t232 = t288 * t463 + t465;
t462 = t305 * t306;
t464 = t304 * t303;
t233 = t288 * t464 - t462;
t389 = -t233 * rSges(7,1) - t232 * rSges(7,2);
t471 = t289 * t304;
t457 = rSges(7,3) * t471 - t282 * t306 + t304 * t317 + t295 - t389;
t404 = t457 * t306;
t294 = t304 * pkin(4);
t469 = t289 * t306;
t253 = pkin(8) * t469 + t294;
t230 = t288 * t462 - t464;
t425 = t288 * t465;
t231 = t425 + t463;
t527 = t231 * rSges(7,1) + t230 * rSges(7,2) + pkin(5) * t425 + t304 * t282 + t469 * t490;
t458 = -t253 + t527;
t542 = -t304 * t458 + t404;
t438 = qJD(3) * t306;
t414 = t288 * t438;
t445 = qJD(1) * t304;
t540 = t289 * t445 + t414;
t135 = (-Icges(7,1) * t305 + t483) * t436 + (Icges(7,5) * t289 + t288 * t371) * qJD(3);
t136 = (-Icges(6,1) * t305 + t485) * t436 + (Icges(6,5) * t289 + t288 * t372) * qJD(3);
t539 = (-t135 - t136) * t303;
t156 = Icges(7,5) * t233 + Icges(7,6) * t232 + Icges(7,3) * t471;
t160 = Icges(7,4) * t233 + Icges(7,2) * t232 + Icges(7,6) * t471;
t164 = Icges(7,1) * t233 + Icges(7,4) * t232 + Icges(7,5) * t471;
t350 = t160 * t305 + t164 * t303;
t64 = t156 * t288 - t289 * t350;
t158 = Icges(6,5) * t233 + Icges(6,6) * t232 + Icges(6,3) * t471;
t162 = Icges(6,4) * t233 + Icges(6,2) * t232 + Icges(6,6) * t471;
t166 = Icges(6,1) * t233 + Icges(6,4) * t232 + Icges(6,5) * t471;
t348 = t162 * t305 + t166 * t303;
t66 = t158 * t288 - t289 * t348;
t498 = t64 + t66;
t155 = Icges(7,5) * t231 + Icges(7,6) * t230 + Icges(7,3) * t469;
t159 = Icges(7,4) * t231 + Icges(7,2) * t230 + Icges(7,6) * t469;
t163 = Icges(7,1) * t231 + Icges(7,4) * t230 + Icges(7,5) * t469;
t351 = t159 * t305 + t163 * t303;
t63 = t155 * t288 - t289 * t351;
t157 = Icges(6,5) * t231 + Icges(6,6) * t230 + Icges(6,3) * t469;
t161 = Icges(6,4) * t231 + Icges(6,2) * t230 + Icges(6,6) * t469;
t165 = Icges(6,1) * t231 + Icges(6,4) * t230 + Icges(6,5) * t469;
t349 = t161 * t305 + t165 * t303;
t65 = t157 * t288 - t289 * t349;
t499 = t63 + t65;
t538 = t304 * t499 - t306 * t498;
t537 = t304 * t498 + t306 * t499;
t536 = qJD(3) / 0.2e1;
t400 = qJD(1) * t288 + qJD(5);
t412 = t289 * t438;
t310 = -t304 * t400 + t412;
t401 = qJD(5) * t288 + qJD(1);
t339 = t401 * t303;
t139 = t305 * t310 - t306 * t339;
t340 = t305 * t401;
t140 = t303 * t310 + t306 * t340;
t70 = Icges(7,5) * t140 + Icges(7,6) * t139 - Icges(7,3) * t540;
t74 = Icges(7,4) * t140 + Icges(7,2) * t139 - Icges(7,6) * t540;
t78 = Icges(7,1) * t140 + Icges(7,4) * t139 - Icges(7,5) * t540;
t19 = (qJD(3) * t351 + t70) * t288 + (qJD(3) * t155 - t303 * t78 - t305 * t74 + (t159 * t303 - t163 * t305) * qJD(5)) * t289;
t72 = Icges(6,5) * t140 + Icges(6,6) * t139 - Icges(6,3) * t540;
t76 = Icges(6,4) * t140 + Icges(6,2) * t139 - Icges(6,6) * t540;
t80 = Icges(6,1) * t140 + Icges(6,4) * t139 - Icges(6,5) * t540;
t21 = (qJD(3) * t349 + t72) * t288 + (qJD(3) * t157 - t303 * t80 - t305 * t76 + (t161 * t303 - t165 * t305) * qJD(5)) * t289;
t535 = t19 + t21;
t338 = t400 * t306;
t137 = t305 * t338 + (t289 * t439 - t339) * t304;
t440 = qJD(3) * t304;
t413 = t289 * t440;
t138 = t304 * t340 + (t338 + t413) * t303;
t416 = t288 * t440;
t444 = qJD(1) * t306;
t419 = t289 * t444;
t316 = -t416 + t419;
t69 = Icges(7,5) * t138 + Icges(7,6) * t137 + Icges(7,3) * t316;
t73 = Icges(7,4) * t138 + Icges(7,2) * t137 + Icges(7,6) * t316;
t77 = Icges(7,1) * t138 + Icges(7,4) * t137 + Icges(7,5) * t316;
t20 = (qJD(3) * t350 + t69) * t288 + (qJD(3) * t156 - t303 * t77 - t305 * t73 + (t160 * t303 - t164 * t305) * qJD(5)) * t289;
t71 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t316;
t75 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t316;
t79 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t316;
t22 = (qJD(3) * t348 + t71) * t288 + (qJD(3) * t158 - t303 * t79 - t305 * t75 + (t162 * t303 - t166 * t305) * qJD(5)) * t289;
t534 = t20 + t22;
t54 = t157 * t469 + t230 * t161 + t231 * t165;
t55 = t158 * t469 + t230 * t162 + t231 * t166;
t382 = t304 * t55 + t306 * t54;
t52 = t155 * t469 + t230 * t159 + t231 * t163;
t53 = t156 * t469 + t230 * t160 + t231 * t164;
t383 = t304 * t53 + t306 * t52;
t87 = t193 * t469 + t230 * t195 + t231 * t197;
t88 = t194 * t469 + t230 * t196 + t231 * t198;
t533 = (t382 + t383) * t289 + (t87 + t88) * t288;
t58 = t157 * t471 + t161 * t232 + t165 * t233;
t59 = t158 * t471 + t162 * t232 + t166 * t233;
t380 = t304 * t59 + t306 * t58;
t56 = t155 * t471 + t159 * t232 + t163 * t233;
t57 = t156 * t471 + t160 * t232 + t164 * t233;
t381 = t304 * t57 + t306 * t56;
t89 = t193 * t471 + t195 * t232 + t197 * t233;
t90 = t194 * t471 + t196 * t232 + t198 * t233;
t500 = (t380 + t381) * t289 + (t89 + t90) * t288;
t343 = t208 * t288 - t210 * t289;
t532 = t304 * t343;
t344 = t205 * t288 - t207 * t289;
t531 = t304 * t344;
t342 = -t288 * t523 + t289 * t522;
t529 = t306 * t342;
t345 = t204 * t288 - t206 * t289;
t528 = t306 * t345;
t293 = t304 * rSges(5,1);
t526 = -rSges(5,2) * t469 + t293;
t434 = qJD(5) * t305;
t411 = t288 * t434;
t433 = qJD(6) * t289;
t525 = t282 * t444 + t412 * t505 + (pkin(5) * t411 + t433) * t306 + t540 * t301 + t140 * rSges(7,1) + t139 * rSges(7,2);
t292 = t304 * rSges(4,3);
t472 = t288 * t306;
t524 = -rSges(4,2) * t472 + t292;
t365 = Icges(4,5) * t289 - Icges(4,6) * t288;
t202 = -Icges(4,3) * t306 + t304 * t365;
t368 = Icges(5,4) * t289 - Icges(5,5) * t288;
t521 = Icges(5,1) * t306 + t304 * t368;
t300 = cos(pkin(9));
t336 = rSges(3,1) * t300 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t492 = rSges(3,3) + qJ(2);
t192 = t304 * t492 + t306 * t336;
t520 = -t304 * t457 - t306 * t458;
t519 = 2 * m(4);
t518 = 2 * m(5);
t517 = 2 * m(6);
t516 = 2 * m(7);
t515 = m(5) / 0.2e1;
t514 = m(6) / 0.2e1;
t513 = m(7) / 0.2e1;
t511 = t304 / 0.2e1;
t510 = -t306 / 0.2e1;
t508 = -rSges(7,3) - pkin(3);
t252 = rSges(4,1) * t288 + rSges(4,2) * t289;
t507 = m(4) * t252;
t506 = pkin(3) * t289;
t265 = pkin(8) * t416;
t390 = t138 * rSges(7,1) + t137 * rSges(7,2);
t428 = pkin(5) * t435;
t443 = qJD(3) * t288;
t497 = rSges(7,3) * t316 + t390 + t265 + (qJD(1) * t317 + t428) * t306 + (t301 * t443 + t433 + (-pkin(4) + t282) * qJD(1) + (t289 * t441 + t411) * pkin(5)) * t304;
t285 = pkin(4) * t444;
t496 = -rSges(7,3) * t540 + pkin(8) * t414 - t285 + (pkin(8) * qJD(1) * t289 - t400 * t505) * t304 + t525;
t495 = rSges(4,1) * t289;
t494 = rSges(5,2) * t288;
t493 = rSges(6,3) * t288;
t491 = -rSges(5,3) - qJ(4);
t475 = qJ(4) * t288;
t474 = qJ(4) * t289;
t392 = -t233 * rSges(6,1) - t232 * rSges(6,2);
t170 = rSges(6,3) * t471 - t392;
t473 = t170 * t306;
t302 = -pkin(7) - qJ(2);
t468 = t302 * t306;
t460 = t140 * rSges(6,1) + t139 * rSges(6,2);
t388 = rSges(7,1) * t303 + rSges(7,2) * t305;
t409 = t289 * t434;
t459 = (-rSges(7,1) * t305 + rSges(7,2) * t303) * t436 - pkin(5) * t409 + qJD(6) * t288 + (rSges(7,3) * t289 + t288 * t388 + t317) * qJD(3);
t456 = (-t388 - t505) * t289 + (rSges(7,3) + t502) * t288;
t385 = t475 + t506;
t215 = qJD(3) * t385 - qJD(4) * t289;
t387 = -rSges(5,2) * t289 + rSges(5,3) * t288;
t455 = -t387 * qJD(3) - t215;
t222 = t385 * t304;
t223 = pkin(3) * t469 + qJ(4) * t472;
t454 = t304 * t222 + t306 * t223;
t453 = -t223 - t253;
t250 = pkin(3) * t288 - t474;
t228 = t250 * t445;
t421 = t288 * t445;
t452 = pkin(8) * t421 + t228;
t386 = rSges(5,3) * t289 + t494;
t451 = -t250 + t386;
t437 = qJD(4) * t288;
t450 = qJ(4) * t412 + t306 * t437;
t291 = qJD(2) * t306;
t449 = t302 * t445 + t291;
t448 = t304 ^ 2 + t306 ^ 2;
t203 = Icges(4,3) * t304 + t306 * t365;
t447 = qJD(1) * t203;
t212 = Icges(5,1) * t304 - t306 * t368;
t446 = qJD(1) * t212;
t432 = -rSges(6,3) - pkin(3) - pkin(8);
t431 = -0.1e1 + t448;
t35 = t52 * t304 - t306 * t53;
t36 = t54 * t304 - t306 * t55;
t427 = -t35 / 0.2e1 - t36 / 0.2e1;
t37 = t56 * t304 - t306 * t57;
t38 = t58 * t304 - t306 * t59;
t426 = t37 / 0.2e1 + t38 / 0.2e1;
t266 = pkin(3) * t416;
t424 = t304 * (pkin(3) * t419 + t304 * t437 - t266 + (t288 * t444 + t413) * qJ(4)) + t306 * (-pkin(3) * t540 - qJ(4) * t421 + t450) + t222 * t444;
t168 = t231 * rSges(6,1) + t230 * rSges(6,2) + rSges(6,3) * t469;
t290 = qJD(2) * t304;
t423 = t290 + t450;
t422 = t266 + t449;
t391 = rSges(6,1) * t303 + rSges(6,2) * t305;
t200 = -t289 * t391 + t493;
t418 = t200 * t445;
t408 = -t368 * qJD(3) / 0.2e1 + t365 * t536;
t407 = -qJ(4) - t505;
t406 = -pkin(8) * t288 - t250;
t405 = t546 * t442 + (((qJD(5) * t544 + t545) * t305 + t539) * t289 + t543) * t288;
t186 = t451 * t306;
t281 = pkin(2) * t300 + pkin(1);
t403 = t306 * t281 - t304 * t302;
t402 = qJD(1) * t456;
t254 = pkin(8) * t471 - t295;
t399 = t306 * t253 + t304 * t254 + t454;
t398 = rSges(5,1) * t444 + rSges(5,2) * t540 + rSges(5,3) * t412;
t397 = -t200 + t406;
t396 = -pkin(8) * t442 - t215;
t395 = t304 * t402;
t394 = -rSges(4,2) * t288 + t495;
t393 = t138 * rSges(6,1) + t137 * rSges(6,2);
t326 = t288 * t407 - t281;
t44 = -t304 * t428 + t508 * t414 + (-t468 + (t289 * t508 + t326) * t304) * qJD(1) + t423 + t525;
t308 = (-pkin(3) - t490) * t289 + t326;
t45 = (qJD(1) * t308 - t428) * t306 + (-qJD(1) * t282 - t433 + (-pkin(5) * t434 - qJD(4)) * t288 + (t288 * t490 + t289 * t407) * qJD(3)) * t304 - t390 + t422;
t384 = t304 * t44 + t306 * t45;
t61 = -t288 * t457 + t456 * t471;
t62 = t288 * t458 - t456 * t469;
t379 = t304 * t62 + t306 * t61;
t93 = (t282 - t302) * t306 + t308 * t304 + t389;
t334 = t403 + t223;
t94 = t334 + t527;
t378 = t304 * t94 + t306 * t93;
t377 = t288 * t491 - t281;
t312 = t289 * t432 - t281 - t475;
t307 = t304 * t312 - t468;
t101 = t295 + t307 + t392;
t102 = t334 + t168 + t253;
t358 = t101 * t306 + t102 * t304;
t341 = t406 - t456;
t115 = t341 * t304;
t116 = t341 * t306;
t357 = t115 * t304 + t116 * t306;
t311 = (rSges(5,2) - pkin(3)) * t289 + t377;
t121 = (rSges(5,1) - t302) * t306 + t311 * t304;
t218 = rSges(5,3) * t472 + t526;
t122 = t218 + t334;
t356 = t121 * t306 + t122 * t304;
t347 = -t168 * t306 - t304 * t170;
t346 = t168 * t304 - t473;
t217 = rSges(4,1) * t469 + t524;
t142 = (-rSges(6,1) * t305 + rSges(6,2) * t303) * t436 + (rSges(6,3) * t289 + t288 * t391) * qJD(3);
t337 = -t142 + t396;
t31 = t131 * t471 + t232 * t133 + t233 * t135 + t137 * t195 + t138 * t197 + t193 * t316;
t32 = t132 * t471 + t232 * t134 + t233 * t136 + t137 * t196 + t138 * t198 + t194 * t316;
t333 = t22 / 0.2e1 + t20 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1;
t33 = t131 * t469 + t230 * t133 + t231 * t135 + t139 * t195 + t140 * t197 - t193 * t540;
t34 = t132 * t469 + t230 * t134 + t231 * t136 + t139 * t196 + t140 * t198 - t194 * t540;
t332 = t33 / 0.2e1 + t34 / 0.2e1 + t21 / 0.2e1 + t19 / 0.2e1;
t331 = -t87 / 0.2e1 - t88 / 0.2e1 - t65 / 0.2e1 - t63 / 0.2e1;
t330 = t89 / 0.2e1 + t90 / 0.2e1 + t66 / 0.2e1 + t64 / 0.2e1;
t124 = t397 * t306;
t329 = -t281 - t394;
t328 = t304 * (qJD(1) * t253 - t265) + t306 * (-pkin(8) * t540 + t285) + t254 * t444 + t424;
t325 = t396 - t459;
t324 = qJD(3) * t252;
t321 = qJD(3) * (Icges(5,4) * t288 + Icges(5,5) * t289);
t320 = qJD(3) * (-Icges(4,5) * t288 - Icges(4,6) * t289);
t24 = (t438 * t456 + t496) * t288 + (qJD(3) * t458 - t306 * t459 + t395) * t289;
t25 = (-t440 * t456 - t497) * t288 + (-qJD(3) * t457 + t304 * t459 + t306 * t402) * t289;
t51 = t542 * t289;
t314 = qJD(3) * t51 + t24 * t304 + t25 * t306;
t46 = t399 - t520;
t47 = t306 * t325 + t395 + t452;
t48 = qJD(1) * t116 + t304 * t325;
t313 = qJD(3) * t46 + t304 * t48 + t306 * t47;
t309 = rSges(4,2) * t421 + rSges(4,3) * t444 - t306 * t324;
t191 = -t304 * t336 + t306 * t492;
t243 = t394 * qJD(3);
t219 = -rSges(5,1) * t306 + t304 * t387;
t216 = -rSges(4,3) * t306 + t304 * t394;
t185 = t451 * t304;
t178 = t217 + t403;
t177 = (rSges(4,3) - t302) * t306 + t329 * t304;
t176 = -qJD(1) * t192 + t291;
t175 = qJD(1) * t191 + t290;
t174 = t431 * t288 * t442;
t154 = qJD(1) * t521 + t306 * t321;
t153 = t304 * t321 + t446;
t144 = t304 * t320 + t447;
t143 = -qJD(1) * t202 + t306 * t320;
t123 = t397 * t304;
t118 = t252 * t440 + (t306 * t329 - t292) * qJD(1) + t449;
t117 = t290 + (-t468 + (-t281 - t495) * t304) * qJD(1) + t309;
t114 = qJD(1) * t186 + t304 * t455;
t113 = t306 * t455 - t386 * t445 + t228;
t112 = t288 * t168 - t200 * t469;
t111 = -t170 * t288 + t200 * t471;
t110 = t304 * t342 + t306 * t521;
t109 = -t212 * t306 + t532;
t108 = -t304 * t521 + t529;
t107 = t304 * t212 + t343 * t306;
t106 = t304 * t203 - t344 * t306;
t105 = t304 * t202 - t528;
t104 = -t203 * t306 - t531;
t103 = -t202 * t306 - t304 * t345;
t100 = t218 * t306 + t304 * t219 + t454;
t95 = t346 * t289;
t86 = (-t437 + (t289 * t491 - t494) * qJD(3)) * t304 + (t306 * t311 - t293) * qJD(1) + t422;
t85 = -pkin(3) * t414 + (-t468 + (t377 - t506) * t304) * qJD(1) + t398 + t423;
t84 = -rSges(6,3) * t540 + t460;
t82 = rSges(6,3) * t316 + t393;
t68 = qJD(1) * t124 + t304 * t337;
t67 = t306 * t337 + t418 + t452;
t60 = -t347 + t399;
t50 = t265 + (-t437 + (-t474 + t493) * qJD(3)) * t304 + (t306 * t312 - t294) * qJD(1) - t393 + t422;
t49 = qJD(1) * t307 + t414 * t432 + t285 + t423 + t460;
t43 = (-t200 * t440 - t82) * t288 + (-qJD(3) * t170 + t304 * t142 + t200 * t444) * t289;
t42 = (t200 * t438 + t84) * t288 + (qJD(3) * t168 - t142 * t306 + t418) * t289;
t41 = (qJD(1) * t219 + t398) * t306 + (t386 * t440 + (-t218 - t223 + t526) * qJD(1)) * t304 + t424;
t30 = t346 * t443 + (qJD(1) * t347 - t304 * t84 + t306 * t82) * t289;
t23 = t304 * t82 + t306 * t84 + (t473 + (-t168 + t453) * t304) * qJD(1) + t328;
t18 = t139 * t162 + t140 * t166 - t158 * t540 + t230 * t75 + t231 * t79 + t469 * t71;
t17 = t139 * t161 + t140 * t165 - t157 * t540 + t230 * t76 + t231 * t80 + t469 * t72;
t16 = t139 * t160 + t140 * t164 - t156 * t540 + t230 * t73 + t231 * t77 + t469 * t69;
t15 = t139 * t159 + t140 * t163 - t155 * t540 + t230 * t74 + t231 * t78 + t469 * t70;
t14 = t137 * t162 + t138 * t166 + t158 * t316 + t232 * t75 + t233 * t79 + t471 * t71;
t13 = t137 * t161 + t138 * t165 + t157 * t316 + t232 * t76 + t233 * t80 + t471 * t72;
t12 = t137 * t160 + t138 * t164 + t156 * t316 + t232 * t73 + t233 * t77 + t471 * t69;
t11 = t137 * t159 + t138 * t163 + t155 * t316 + t232 * t74 + t233 * t78 + t471 * t70;
t10 = t496 * t306 + t497 * t304 + (t404 + (t453 - t458) * t304) * qJD(1) + t328;
t9 = -t542 * t443 + (qJD(1) * t520 - t496 * t304 + t497 * t306) * t289;
t8 = qJD(1) * t382 + t17 * t304 - t18 * t306;
t7 = qJD(1) * t383 + t15 * t304 - t16 * t306;
t6 = qJD(1) * t380 + t13 * t304 - t14 * t306;
t5 = qJD(1) * t381 + t11 * t304 - t12 * t306;
t4 = (-qJD(3) * t382 + t34) * t288 + (-qJD(1) * t36 + qJD(3) * t88 + t17 * t306 + t18 * t304) * t289;
t3 = (-qJD(3) * t383 + t33) * t288 + (-qJD(1) * t35 + qJD(3) * t87 + t15 * t306 + t16 * t304) * t289;
t2 = (-qJD(3) * t380 + t32) * t288 + (-qJD(1) * t38 + qJD(3) * t90 + t13 * t306 + t14 * t304) * t289;
t1 = (-qJD(3) * t381 + t31) * t288 + (-qJD(1) * t37 + qJD(3) * t89 + t11 * t306 + t12 * t304) * t289;
t26 = [0.2e1 * m(3) * (t175 * t192 + t176 * t191) + (t117 * t178 + t118 * t177) * t519 + (t121 * t86 + t122 * t85) * t518 + (t44 * t94 + t45 * t93) * t516 + (t101 * t50 + t102 * t49) * t517 + t544 * t409 + (t374 + t362 + t556) * t443 + (t370 + t360 - t555) * t442 + (t305 * t545 + t539) * t289 + t543; m(7) * (qJD(1) * t378 + t304 * t45 - t306 * t44) + m(6) * (qJD(1) * t358 + t304 * t50 - t306 * t49) + m(4) * (-t117 * t306 + t304 * t118 + (t177 * t306 + t178 * t304) * qJD(1)) + m(5) * (qJD(1) * t356 + t304 * t86 - t306 * t85) + m(3) * (-t175 * t306 + t304 * t176 + (t191 * t306 + t192 * t304) * qJD(1)); 0; m(5) * (t113 * t121 + t114 * t122 + t185 * t85 + t186 * t86) + m(6) * (t101 * t67 + t102 * t68 + t123 * t49 + t124 * t50) + m(7) * (t115 * t44 + t116 * t45 + t47 * t93 + t48 * t94) + (m(4) * (-t118 * t252 - t177 * t243) + t408 * t306 - t333) * t306 + (m(4) * (-t117 * t252 - t178 * t243) + t408 * t304 + t332) * t304 + ((-t549 * qJD(3) + t551 * t306) * t511 + (t550 * qJD(3) + t551 * t304) * t510 + (t510 * t547 - t511 * t548) * qJD(1)) * t288 + ((t547 * qJD(3) + t552 * t306) * t511 + (t548 * qJD(3) + t552 * t304) * t510 + (t510 * t549 + t511 * t550) * qJD(1)) * t289 + ((-t178 * t507 + (t205 / 0.2e1 - t208 / 0.2e1) * t289 + (t207 / 0.2e1 - t210 / 0.2e1) * t288 - t331) * t306 + (t177 * t507 + (t204 / 0.2e1 + t523 / 0.2e1) * t289 + (t206 / 0.2e1 + t522 / 0.2e1) * t288 + t330) * t304) * qJD(1); m(5) * (t113 * t304 - t114 * t306 + (t185 * t304 + t186 * t306) * qJD(1)) + m(6) * (t67 * t304 - t306 * t68 + (t123 * t304 + t124 * t306) * qJD(1)) + m(7) * (qJD(1) * t357 + t47 * t304 - t306 * t48); (t10 * t46 + t115 * t48 + t116 * t47) * t516 + (t123 * t68 + t124 * t67 + t23 * t60) * t517 - t306 * t5 - t306 * t6 + t304 * t8 + t304 * t7 + (t100 * t41 + t113 * t186 + t114 * t185) * t518 + t304 * ((t304 * t143 + (t105 + t531) * qJD(1)) * t304 + (t106 * qJD(1) + (t204 * t442 + t206 * t443) * t306 + (-t144 + (-t205 * t289 - t207 * t288) * qJD(3) + (t203 - t345) * qJD(1)) * t304) * t306) - t306 * ((t144 * t306 + (t104 + t528) * qJD(1)) * t306 + (t103 * qJD(1) + (-t205 * t442 - t207 * t443 + t447) * t304 + (-t143 + (t204 * t289 + t206 * t288) * qJD(3) - t344 * qJD(1)) * t306) * t304) + ((t304 * t216 + t217 * t306) * ((qJD(1) * t216 + t309) * t306 + (-t304 * t324 + (-t217 + t524) * qJD(1)) * t304) + t448 * t252 * t243) * t519 + t304 * ((t304 * t154 + (t108 - t532) * qJD(1)) * t304 + (t107 * qJD(1) + (t442 * t523 + t443 * t522) * t306 + (-t153 + (t208 * t289 + t210 * t288) * qJD(3) + (t212 + t342) * qJD(1)) * t304) * t306) - t306 * ((t153 * t306 + (t109 - t529) * qJD(1)) * t306 + (t110 * qJD(1) + (t208 * t442 + t210 * t443 + t446) * t304 + (-t154 + (t288 * t522 + t289 * t523) * qJD(3) + t343 * qJD(1)) * t306) * t304) + (t37 + t38 + (-t103 - t110) * t306 + (t104 + t109) * t304) * t445 + (t35 + t36 + (-t105 - t108) * t306 + (t106 + t107) * t304) * t444; 0.2e1 * (t356 * t515 + t358 * t514 + t378 * t513) * t442 + 0.2e1 * ((t444 * t94 - t445 * t93 + t384) * t513 + (-t101 * t445 + t102 * t444 + t304 * t49 + t306 * t50) * t514 + (-t121 * t445 + t122 * t444 + t304 * t85 + t306 * t86) * t515) * t288; 0; 0.2e1 * ((t115 * t440 + t116 * t438 - t10) * t513 + (t123 * t440 + t124 * t438 - t23) * t514 + (t185 * t440 + t186 * t438 - t41) * t515) * t289 + 0.2e1 * ((t115 * t444 - t116 * t445 + t313) * t513 + (qJD(3) * t60 + t123 * t444 - t124 * t445 + t304 * t68 + t306 * t67) * t514 + (qJD(3) * t100 + t113 * t306 + t114 * t304 + t185 * t444 - t186 * t445) * t515) * t288; 0.4e1 * (t515 + t514 + t513) * t174; m(6) * (t101 * t43 + t102 * t42 + t111 * t50 + t112 * t49) + m(7) * (t24 * t94 + t25 * t93 + t44 * t62 + t45 * t61) + (-t304 * t330 + t306 * t331) * t443 + (t332 * t306 + t333 * t304 + (t304 * t331 + t306 * t330) * qJD(1)) * t289 + t405; m(6) * (t43 * t304 - t306 * t42 + (t111 * t306 + t112 * t304) * qJD(1)) + m(7) * (qJD(1) * t379 - t24 * t306 + t25 * t304); m(6) * (t111 * t67 + t112 * t68 + t123 * t42 + t124 * t43 - t23 * t95 + t30 * t60) + m(7) * (t10 * t51 + t115 * t24 + t116 * t25 + t46 * t9 + t47 * t61 + t48 * t62) + (-t2 / 0.2e1 - t1 / 0.2e1 + t427 * t443) * t306 + (t3 / 0.2e1 + t4 / 0.2e1 - t426 * t443) * t304 + ((t304 * t427 + t306 * t426) * qJD(1) + (t5 + t6) * t511 + (t7 + t8) * t306 / 0.2e1 + t538 * t536) * t289 + (qJD(1) * t537 + t535 * t304 - t534 * t306) * t288 / 0.2e1 + (t500 * t304 + t533 * t306) * qJD(1) / 0.2e1; 0.2e1 * ((t111 * t438 + t112 * t440 - t30) * t514 + (t438 * t61 + t440 * t62 - t9) * t513) * t289 + 0.2e1 * ((-qJD(3) * t95 - t111 * t445 + t112 * t444 + t304 * t42 + t306 * t43) * t514 + (t444 * t62 - t445 * t61 + t314) * t513) * t288; (t24 * t62 + t25 * t61 + t51 * t9) * t516 + (t111 * t43 + t112 * t42 - t30 * t95) * t517 + (((-t288 * t499 - t533) * t306 + (-t288 * t498 - t500) * t304) * qJD(3) + t405) * t288 + ((t3 + t4) * t306 + (t1 + t2) * t304 + (t534 * t304 + t535 * t306) * t288 + (t288 * t546 + t537 * t289) * qJD(3) + (-t288 * t538 - t304 * t533 + t500 * t306) * qJD(1)) * t289; m(7) * (-t378 * t443 + ((-t304 * t93 + t306 * t94) * qJD(1) + t384) * t289); 0; m(7) * ((-qJD(3) * t357 + t10) * t288 + ((t115 * t306 - t116 * t304) * qJD(1) + t313) * t289); m(7) * (-t288 ^ 2 + t289 ^ 2) * t431 * qJD(3); m(7) * ((-qJD(3) * t379 + t9) * t288 + ((-t304 * t61 + t306 * t62) * qJD(1) + t314) * t289); -0.2e1 * m(7) * t174;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t26(1) t26(2) t26(4) t26(7) t26(11) t26(16); t26(2) t26(3) t26(5) t26(8) t26(12) t26(17); t26(4) t26(5) t26(6) t26(9) t26(13) t26(18); t26(7) t26(8) t26(9) t26(10) t26(14) t26(19); t26(11) t26(12) t26(13) t26(14) t26(15) t26(20); t26(16) t26(17) t26(18) t26(19) t26(20) t26(21);];
Mq  = res;
