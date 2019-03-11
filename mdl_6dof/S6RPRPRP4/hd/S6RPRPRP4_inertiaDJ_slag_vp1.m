% Calculate time derivative of joint inertia matrix for
% S6RPRPRP4
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:14
% EndTime: 2019-03-09 03:11:36
% DurationCPUTime: 17.16s
% Computational Cost: add. (22609->803), mult. (30974->1095), div. (0->0), fcn. (29386->8), ass. (0->383)
t562 = Icges(5,1) + Icges(4,3);
t296 = cos(qJ(3));
t293 = sin(qJ(3));
t454 = Icges(5,6) * t293;
t464 = Icges(4,4) * t293;
t561 = -t454 - t464 + (-Icges(4,2) - Icges(5,3)) * t296;
t453 = Icges(5,6) * t296;
t463 = Icges(4,4) * t296;
t560 = -t453 - t463 + (-Icges(4,1) - Icges(5,2)) * t293;
t348 = Icges(4,5) * t296 - Icges(4,6) * t293;
t351 = Icges(5,4) * t296 - Icges(5,5) * t293;
t559 = t348 - t351;
t290 = qJ(1) + pkin(9);
t288 = cos(t290);
t287 = sin(t290);
t357 = Icges(4,1) * t296 - t464;
t193 = Icges(4,5) * t287 + t288 * t357;
t345 = Icges(5,2) * t296 - t454;
t196 = Icges(5,4) * t287 - t288 * t345;
t507 = t193 - t196;
t353 = -Icges(4,2) * t293 + t463;
t191 = Icges(4,6) * t287 + t288 * t353;
t343 = -Icges(5,3) * t293 + t453;
t194 = Icges(5,5) * t287 - t288 * t343;
t509 = t191 - t194;
t545 = -t509 * t293 + t507 * t296;
t558 = t288 * t545;
t292 = sin(qJ(5));
t295 = cos(qJ(5));
t456 = Icges(7,5) * t295;
t354 = Icges(7,1) * t292 - t456;
t225 = Icges(7,4) * t293 - t296 * t354;
t460 = Icges(6,4) * t295;
t355 = Icges(6,1) * t292 + t460;
t226 = Icges(6,5) * t293 - t296 * t355;
t555 = t225 + t226;
t557 = t555 * t292;
t190 = -Icges(4,6) * t288 + t287 * t353;
t501 = Icges(5,5) * t288 + t287 * t343;
t510 = t190 + t501;
t192 = -Icges(4,5) * t288 + t287 * t357;
t500 = Icges(5,4) * t288 + t287 * t345;
t508 = t192 + t500;
t547 = t287 * t562 + t559 * t288;
t347 = Icges(6,5) * t292 + Icges(6,6) * t295;
t222 = Icges(6,3) * t293 - t296 * t347;
t349 = Icges(7,4) * t292 - Icges(7,6) * t295;
t223 = Icges(7,2) * t293 - t296 * t349;
t556 = t222 + t223;
t554 = t561 * qJD(3);
t553 = t560 * qJD(3);
t381 = qJD(5) * t293 + qJD(1);
t415 = qJD(3) * t296;
t392 = t295 * t415;
t420 = qJD(1) * t293;
t399 = t288 * t420;
t444 = t288 * t295;
t448 = t287 * t292;
t129 = -qJD(5) * t444 - t287 * t392 - t295 * t399 + t381 * t448;
t301 = t292 * t415 + t295 * t381;
t380 = qJD(5) + t420;
t446 = t288 * t292;
t130 = t287 * t301 + t380 * t446;
t440 = t293 * t295;
t219 = t287 * t440 + t446;
t466 = rSges(7,3) + qJ(6);
t541 = rSges(7,1) + pkin(5);
t552 = t219 * qJD(6) - t466 * t129 - t541 * t130;
t442 = t292 * t293;
t220 = t287 * t442 - t444;
t551 = -t219 * t466 + t541 * t220;
t503 = -t510 * t293 + t508 * t296;
t549 = t288 * t503;
t548 = t559 * t287 - t288 * t562;
t546 = ((Icges(5,5) - Icges(4,6)) * t296 + (Icges(5,4) - Icges(4,5)) * t293) * qJD(3);
t217 = -t288 * t440 + t448;
t218 = t287 * t295 + t288 * t442;
t443 = t288 * t296;
t135 = Icges(6,5) * t218 - Icges(6,6) * t217 + Icges(6,3) * t443;
t139 = Icges(6,4) * t218 - Icges(6,2) * t217 + Icges(6,6) * t443;
t143 = Icges(6,1) * t218 - Icges(6,4) * t217 + Icges(6,5) * t443;
t447 = t287 * t296;
t57 = t135 * t447 + t139 * t219 + t143 * t220;
t136 = Icges(6,5) * t220 + Icges(6,6) * t219 + Icges(6,3) * t447;
t140 = Icges(6,4) * t220 + Icges(6,2) * t219 + Icges(6,6) * t447;
t144 = Icges(6,1) * t220 + Icges(6,4) * t219 + Icges(6,5) * t447;
t58 = t136 * t447 + t140 * t219 + t144 * t220;
t360 = t287 * t58 + t288 * t57;
t133 = Icges(7,5) * t218 + Icges(7,6) * t443 + Icges(7,3) * t217;
t137 = Icges(7,4) * t218 + Icges(7,2) * t443 + Icges(7,6) * t217;
t141 = Icges(7,1) * t218 + Icges(7,4) * t443 + Icges(7,5) * t217;
t55 = -t133 * t219 + t137 * t447 + t141 * t220;
t134 = Icges(7,5) * t220 + Icges(7,6) * t447 - Icges(7,3) * t219;
t138 = Icges(7,4) * t220 + Icges(7,2) * t447 - Icges(7,6) * t219;
t142 = Icges(7,1) * t220 + Icges(7,4) * t447 - Icges(7,5) * t219;
t56 = -t134 * t219 + t138 * t447 + t142 * t220;
t361 = t287 * t56 + t288 * t55;
t544 = t360 + t361;
t53 = t135 * t443 - t139 * t217 + t143 * t218;
t54 = t136 * t443 - t140 * t217 + t144 * t218;
t362 = t287 * t54 + t288 * t53;
t51 = t133 * t217 + t137 * t443 + t141 * t218;
t52 = t134 * t217 + t138 * t443 + t142 * t218;
t363 = t287 * t52 + t288 * t51;
t543 = t362 + t363;
t457 = Icges(7,5) * t292;
t346 = -Icges(7,3) * t295 + t457;
t221 = Icges(7,6) * t293 - t296 * t346;
t461 = Icges(6,4) * t292;
t350 = Icges(6,2) * t295 + t461;
t224 = Icges(6,6) * t293 - t296 * t350;
t540 = ((t221 - t224) * t295 - t557) * t296 + t556 * t293;
t413 = qJD(5) * t296;
t175 = (-Icges(7,1) * t295 - t457) * t413 + (Icges(7,4) * t296 + t293 * t354) * qJD(3);
t176 = (-Icges(6,1) * t295 + t461) * t413 + (Icges(6,5) * t296 + t293 * t355) * qJD(3);
t539 = -t176 - t175;
t171 = (-Icges(7,3) * t292 - t456) * t413 + (Icges(7,6) * t296 + t293 * t346) * qJD(3);
t172 = (-Icges(6,5) * t295 + Icges(6,6) * t292) * t413 + (Icges(6,3) * t296 + t293 * t347) * qJD(3);
t173 = (-Icges(7,4) * t295 - Icges(7,6) * t292) * t413 + (Icges(7,2) * t296 + t293 * t349) * qJD(3);
t416 = qJD(3) * t293;
t393 = t295 * t416;
t439 = t295 * t296;
t390 = t292 * t413;
t535 = t390 + t393;
t537 = t171 * t439 - t221 * t393 + t535 * t224 + t556 * t415 + t416 * t557 + (t172 + t173) * t293;
t435 = rSges(7,2) * t447 + t551;
t384 = t435 * t288;
t436 = rSges(7,2) * t443 + t466 * t217 + t541 * t218;
t536 = -t287 * t436 + t384;
t396 = t288 * t416;
t419 = qJD(1) * t296;
t534 = t287 * t419 + t396;
t533 = t547 * qJD(1);
t532 = -t547 * t288 + t549 + (t545 + t548) * t287;
t340 = t134 * t295 - t142 * t292;
t63 = t138 * t293 + t296 * t340;
t338 = t140 * t295 + t144 * t292;
t65 = t136 * t293 - t296 * t338;
t476 = t63 + t65;
t341 = t133 * t295 - t141 * t292;
t62 = t137 * t293 + t296 * t341;
t339 = t139 * t295 + t143 * t292;
t64 = t135 * t293 - t296 * t339;
t477 = t62 + t64;
t531 = t287 * t476 + t288 * t477;
t530 = t287 * t507 + t288 * t508;
t529 = t287 * t509 + t288 * t510;
t528 = t287 ^ 2;
t527 = t288 ^ 2;
t397 = t287 * t416;
t398 = t288 * t419;
t307 = -t397 + t398;
t131 = qJD(1) * t219 + qJD(5) * t218 - t288 * t392;
t132 = t288 * t301 - t380 * t448;
t71 = Icges(7,5) * t132 - Icges(7,6) * t534 + Icges(7,3) * t131;
t75 = Icges(7,4) * t132 - Icges(7,2) * t534 + Icges(7,6) * t131;
t79 = Icges(7,1) * t132 - Icges(7,4) * t534 + Icges(7,5) * t131;
t11 = t129 * t133 + t130 * t141 + t137 * t307 - t219 * t71 + t220 * t79 + t447 * t75;
t70 = Icges(7,5) * t130 + Icges(7,6) * t307 + Icges(7,3) * t129;
t74 = Icges(7,4) * t130 + Icges(7,2) * t307 + Icges(7,6) * t129;
t78 = Icges(7,1) * t130 + Icges(7,4) * t307 + Icges(7,5) * t129;
t12 = t129 * t134 + t130 * t142 + t138 * t307 - t219 * t70 + t220 * t78 + t447 * t74;
t73 = Icges(6,5) * t132 - Icges(6,6) * t131 - Icges(6,3) * t534;
t77 = Icges(6,4) * t132 - Icges(6,2) * t131 - Icges(6,6) * t534;
t81 = Icges(6,1) * t132 - Icges(6,4) * t131 - Icges(6,5) * t534;
t13 = -t129 * t139 + t130 * t143 + t135 * t307 + t219 * t77 + t220 * t81 + t447 * t73;
t72 = Icges(6,5) * t130 - Icges(6,6) * t129 + Icges(6,3) * t307;
t76 = Icges(6,4) * t130 - Icges(6,2) * t129 + Icges(6,6) * t307;
t80 = Icges(6,1) * t130 - Icges(6,4) * t129 + Icges(6,5) * t307;
t14 = -t129 * t140 + t130 * t144 + t136 * t307 + t219 * t76 + t220 * t80 + t447 * t72;
t526 = (-t12 - t14) * t288 + (t11 + t13) * t287 + t544 * qJD(1);
t15 = t131 * t133 + t132 * t141 - t137 * t534 + t217 * t71 + t218 * t79 + t443 * t75;
t16 = t131 * t134 + t132 * t142 - t138 * t534 + t217 * t70 + t218 * t78 + t443 * t74;
t17 = -t131 * t139 + t132 * t143 - t135 * t534 - t217 * t77 + t218 * t81 + t443 * t73;
t18 = -t131 * t140 + t132 * t144 - t136 * t534 - t217 * t76 + t218 * t80 + t443 * t72;
t525 = (-t16 - t18) * t288 + (t15 + t17) * t287 + t543 * qJD(1);
t524 = qJD(3) / 0.2e1;
t19 = (-qJD(3) * t341 + t75) * t293 + (qJD(3) * t137 - t292 * t79 + t295 * t71 + (-t133 * t292 - t141 * t295) * qJD(5)) * t296;
t21 = (qJD(3) * t339 + t73) * t293 + (qJD(3) * t135 - t292 * t81 - t295 * t77 + (t139 * t292 - t143 * t295) * qJD(5)) * t296;
t523 = t19 + t21;
t20 = (-qJD(3) * t340 + t74) * t293 + (qJD(3) * t138 - t292 * t78 + t295 * t70 + (-t134 * t292 - t142 * t295) * qJD(5)) * t296;
t22 = (qJD(3) * t338 + t72) * t293 + (qJD(3) * t136 - t292 * t80 - t295 * t76 + (t140 * t292 - t144 * t295) * qJD(5)) * t296;
t522 = t20 + t22;
t94 = t217 * t221 + t218 * t225 + t223 * t443;
t95 = -t217 * t224 + t218 * t226 + t222 * t443;
t521 = t543 * t296 + (t94 + t95) * t293;
t96 = -t219 * t221 + t220 * t225 + t223 * t447;
t97 = t219 * t224 + t220 * t226 + t222 * t447;
t520 = t544 * t296 + (t96 + t97) * t293;
t496 = 2 * m(4);
t473 = rSges(4,1) * t296;
t371 = -rSges(4,2) * t293 + t473;
t469 = rSges(4,3) * t288;
t201 = t287 * t371 - t469;
t445 = t288 * t293;
t506 = -rSges(4,2) * t445 + t287 * rSges(4,3);
t202 = rSges(4,1) * t443 + t506;
t266 = rSges(4,1) * t293 + rSges(4,2) * t296;
t316 = qJD(3) * t266;
t401 = t287 * t420;
t421 = qJD(1) * t288;
t300 = rSges(4,2) * t401 + rSges(4,3) * t421 - t288 * t316;
t61 = (qJD(1) * t201 + t300) * t288 + (-t287 * t316 + (-t202 + t506) * qJD(1)) * t287;
t519 = t496 * t61;
t518 = qJD(6) * t217 + t466 * t131 + t132 * t541;
t516 = -t287 * t503 + t288 * t548;
t513 = t287 * t547 + t558;
t512 = -qJD(1) * t548 + t288 * t546;
t511 = -t287 * t546 - t533;
t505 = t287 * rSges(5,1) - rSges(5,2) * t443;
t233 = t287 * pkin(4) + pkin(8) * t443;
t282 = t288 * pkin(7);
t283 = t288 * pkin(4);
t504 = t282 + t283;
t422 = qJD(1) * t287;
t498 = t293 * t476 + t520;
t497 = -t287 * t435 - t288 * t436;
t495 = 2 * m(5);
t494 = 2 * m(6);
t493 = 2 * m(7);
t492 = m(5) / 0.2e1;
t491 = m(6) / 0.2e1;
t490 = m(7) / 0.2e1;
t489 = -pkin(3) - pkin(8);
t488 = t287 / 0.2e1;
t487 = -t288 / 0.2e1;
t483 = m(4) * t266;
t482 = sin(qJ(1)) * pkin(1);
t481 = pkin(3) * t296;
t289 = cos(qJ(1)) * pkin(1);
t475 = rSges(7,2) * t307 - t552;
t474 = -rSges(7,2) * t534 + t518;
t472 = rSges(5,1) * t288;
t471 = rSges(5,2) * t293;
t470 = rSges(7,2) * t293;
t468 = rSges(6,3) * t293;
t467 = -rSges(5,3) - qJ(4);
t451 = qJ(4) * t293;
t450 = qJ(4) * t296;
t369 = -rSges(6,1) * t220 - rSges(6,2) * t219;
t148 = rSges(6,3) * t447 - t369;
t449 = t148 * t288;
t437 = t132 * rSges(6,1) - t131 * rSges(6,2);
t367 = rSges(7,1) * t292 - rSges(7,3) * t295;
t434 = (pkin(5) * t416 - qJ(6) * t413) * t292 + (-qJ(6) * t416 + (-pkin(5) * qJD(5) + qJD(6)) * t296) * t295 + (-rSges(7,1) * t295 - rSges(7,3) * t292) * t413 + (rSges(7,2) * t296 + t293 * t367) * qJD(3);
t364 = t451 + t481;
t227 = t364 * t287;
t228 = pkin(3) * t443 + qJ(4) * t445;
t433 = t287 * t227 + t288 * t228;
t432 = -t228 - t233;
t431 = t470 + (-pkin(5) * t292 + qJ(6) * t295 - t367) * t296;
t231 = qJD(3) * t364 - qJD(4) * t296;
t366 = -rSges(5,2) * t296 + rSges(5,3) * t293;
t430 = -t366 * qJD(3) - t231;
t264 = pkin(3) * t293 - t450;
t232 = t264 * t422;
t429 = pkin(8) * t401 + t232;
t395 = t288 * t415;
t414 = qJD(4) * t293;
t428 = qJ(4) * t395 + t288 * t414;
t253 = pkin(8) * t397;
t254 = pkin(3) * t397;
t427 = t253 + t254;
t365 = rSges(5,3) * t296 + t471;
t426 = -t264 + t365;
t425 = t527 + t528;
t418 = qJD(3) * t287;
t417 = qJD(3) * t288;
t411 = -rSges(7,2) + t489;
t410 = -rSges(6,3) + t489;
t31 = t287 * t51 - t288 * t52;
t32 = t287 * t53 - t288 * t54;
t407 = -t31 / 0.2e1 - t32 / 0.2e1;
t33 = t287 * t55 - t288 * t56;
t34 = t287 * t57 - t288 * t58;
t406 = t33 / 0.2e1 + t34 / 0.2e1;
t405 = t287 * (pkin(3) * t398 + t287 * t414 - t254 + (t287 * t415 + t399) * qJ(4)) + t288 * (-pkin(3) * t534 - qJ(4) * t401 + t428) + t227 * t421;
t146 = t218 * rSges(6,1) - t217 * rSges(6,2) + rSges(6,3) * t443;
t276 = pkin(7) * t421;
t404 = t276 + t428;
t403 = t288 * pkin(2) + t287 * pkin(7) + t289;
t368 = rSges(6,1) * t292 + rSges(6,2) * t295;
t230 = -t296 * t368 + t468;
t402 = t230 * t422;
t391 = qJD(5) * t296 ^ 2 * t292;
t388 = -t351 * qJD(3) / 0.2e1 + t348 * t524;
t387 = t282 - t482;
t386 = -pkin(2) - t451;
t385 = -pkin(8) * t293 - t264;
t183 = t426 * t288;
t383 = qJD(1) * t431;
t174 = (Icges(6,2) * t292 - t460) * t413 + (Icges(6,6) * t296 + t293 * t350) * qJD(3);
t382 = t540 * t415 + (((-qJD(5) * t555 - t174) * t295 + (-t221 * qJD(5) + t539) * t292) * t296 + t537) * t293;
t234 = pkin(8) * t447 - t283;
t379 = t288 * t233 + t287 * t234 + t433;
t277 = pkin(4) * t421;
t378 = t277 + t404;
t377 = rSges(5,1) * t421 + t534 * rSges(5,2) + rSges(5,3) * t395;
t376 = -t230 + t385;
t375 = -pkin(8) * t415 - t231;
t374 = t287 * t383;
t373 = t293 * t467 - pkin(2);
t370 = t130 * rSges(6,1) - t129 * rSges(6,2);
t358 = t403 + t228;
t337 = -t146 * t288 - t148 * t287;
t336 = t146 * t287 - t449;
t327 = t385 - t431;
t203 = rSges(5,3) * t445 + t505;
t178 = (-rSges(6,1) * t295 + rSges(6,2) * t292) * t413 + (rSges(6,3) * t296 + t293 * t368) * qJD(3);
t326 = -t178 + t375;
t325 = -t293 * t477 - t521;
t324 = -pkin(2) - t371;
t37 = t131 * t221 + t132 * t225 + t171 * t217 + t173 * t443 + t175 * t218 - t223 * t534;
t38 = -t131 * t224 + t132 * t226 + t172 * t443 - t174 * t217 + t176 * t218 - t222 * t534;
t323 = t21 / 0.2e1 + t19 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1;
t35 = t129 * t221 + t130 * t225 - t171 * t219 + t173 * t447 + t175 * t220 + t223 * t307;
t36 = -t129 * t224 + t130 * t226 + t172 * t447 + t174 * t219 + t176 * t220 + t222 * t307;
t322 = t22 / 0.2e1 + t20 / 0.2e1 + t35 / 0.2e1 + t36 / 0.2e1;
t321 = -t94 / 0.2e1 - t95 / 0.2e1 - t64 / 0.2e1 - t62 / 0.2e1;
t320 = -t96 / 0.2e1 - t97 / 0.2e1 - t65 / 0.2e1 - t63 / 0.2e1;
t162 = t376 * t288;
t319 = -t289 + (-pkin(4) - pkin(7)) * t287;
t318 = t287 * (t233 * qJD(1) - t253) + t288 * (-pkin(8) * t534 + t277) + t234 * t421 + t405;
t317 = t375 - t434;
t119 = t327 * t288;
t309 = t358 + t233;
t305 = t296 * t411 + t386;
t304 = t296 * t410 + t386;
t303 = (rSges(5,2) - pkin(3)) * t296 + t373;
t299 = t287 * t305 - t482;
t298 = t287 * t304 - t482;
t243 = t371 * qJD(3);
t204 = t287 * t366 - t472;
t182 = t426 * t287;
t169 = t202 + t403;
t168 = t287 * t324 + t387 + t469;
t161 = t376 * t287;
t121 = t203 + t358;
t120 = t287 * t303 + t387 + t472;
t118 = t327 * t287;
t117 = t266 * t418 + (-t289 + (-rSges(4,3) - pkin(7)) * t287 + t324 * t288) * qJD(1);
t116 = t276 + (-t482 + (-pkin(2) - t473) * t287) * qJD(1) + t300;
t115 = qJD(1) * t183 + t287 * t430;
t114 = t288 * t430 - t365 * t422 + t232;
t109 = t146 * t293 - t230 * t443;
t108 = -t148 * t293 + t230 * t447;
t99 = t309 + t146;
t98 = t298 + t369 + t504;
t93 = t203 * t288 + t204 * t287 + t433;
t92 = t336 * t296;
t91 = t254 + (-t414 + (t296 * t467 - t471) * qJD(3)) * t287 + (-t289 + (-rSges(5,1) - pkin(7)) * t287 + t303 * t288) * qJD(1);
t90 = -pkin(3) * t396 + (-t482 + (t373 - t481) * t287) * qJD(1) + t377 + t404;
t89 = qJD(1) * t162 + t287 * t326;
t88 = t288 * t326 + t402 + t429;
t87 = -rSges(6,3) * t534 + t437;
t85 = rSges(6,3) * t307 + t370;
t69 = t309 + t436;
t68 = t299 + t504 - t551;
t67 = t293 * t436 - t431 * t443;
t66 = -t293 * t435 + t431 * t447;
t60 = qJD(1) * t119 + t287 * t317;
t59 = t288 * t317 + t374 + t429;
t50 = -t337 + t379;
t49 = t536 * t296;
t48 = (-t414 + (-t450 + t468) * qJD(3)) * t287 + (t288 * t304 + t319) * qJD(1) - t370 + t427;
t47 = qJD(1) * t298 + t396 * t410 + t378 + t437;
t44 = (-t230 * t418 - t85) * t293 + (-qJD(3) * t148 + t178 * t287 + t230 * t421) * t296;
t43 = (t230 * t417 + t87) * t293 + (qJD(3) * t146 - t178 * t288 + t402) * t296;
t42 = t379 - t497;
t41 = (qJD(1) * t204 + t377) * t288 + (t365 * t418 + (-t203 - t228 + t505) * qJD(1)) * t287 + t405;
t40 = (-t414 + (-t450 + t470) * qJD(3)) * t287 + (t288 * t305 + t319) * qJD(1) + t427 + t552;
t39 = qJD(1) * t299 + t396 * t411 + t378 + t518;
t30 = t336 * t416 + (qJD(1) * t337 - t287 * t87 + t288 * t85) * t296;
t25 = (-t418 * t431 - t475) * t293 + (-qJD(3) * t435 + t287 * t434 + t288 * t383) * t296;
t24 = (t417 * t431 + t474) * t293 + (qJD(3) * t436 - t288 * t434 + t374) * t296;
t23 = t287 * t85 + t288 * t87 + (t449 + (-t146 + t432) * t287) * qJD(1) + t318;
t10 = -t536 * t416 + (qJD(1) * t497 - t474 * t287 + t475 * t288) * t296;
t9 = t474 * t288 + t475 * t287 + (t384 + (t432 - t436) * t287) * qJD(1) + t318;
t4 = (-qJD(3) * t362 + t38) * t293 + (-qJD(1) * t32 + qJD(3) * t95 + t17 * t288 + t18 * t287) * t296;
t3 = (-qJD(3) * t363 + t37) * t293 + (-qJD(1) * t31 + qJD(3) * t94 + t15 * t288 + t16 * t287) * t296;
t2 = (-qJD(3) * t360 + t36) * t293 + (-qJD(1) * t34 + qJD(3) * t97 + t13 * t288 + t14 * t287) * t296;
t1 = (-qJD(3) * t361 + t35) * t293 + (-qJD(1) * t33 + qJD(3) * t96 + t11 * t288 + t12 * t287) * t296;
t5 = [-t174 * t439 + (t116 * t169 + t117 * t168) * t496 + (t120 * t91 + t121 * t90) * t495 + (t47 * t99 + t48 * t98) * t494 + (t39 * t69 + t40 * t68) * t493 - t221 * t390 + t539 * t292 * t296 - t555 * t295 * t413 + (t357 + t345 + t561) * t416 + (t353 + t343 - t560) * t415 + t537; 0; 0; m(5) * (t114 * t120 + t115 * t121 + t182 * t90 + t183 * t91) + m(6) * (t161 * t47 + t162 * t48 + t88 * t98 + t89 * t99) + m(7) * (t118 * t39 + t119 * t40 + t59 * t68 + t60 * t69) + (m(4) * (-t117 * t266 - t168 * t243) + t388 * t288 - t322) * t288 + (m(4) * (-t116 * t266 - t169 * t243) + t388 * t287 + t323) * t287 + ((-t509 * qJD(3) + t553 * t288) * t488 + (-t510 * qJD(3) + t553 * t287) * t487 + (t487 * t507 - t488 * t508) * qJD(1)) * t293 + ((t507 * qJD(3) + t554 * t288) * t488 + (t508 * qJD(3) + t554 * t287) * t487 + (t487 * t509 - t488 * t510) * qJD(1)) * t296 + ((-t169 * t483 + (-t194 / 0.2e1 + t191 / 0.2e1) * t296 + (-t196 / 0.2e1 + t193 / 0.2e1) * t293 - t321) * t288 + (t168 * t483 + (t501 / 0.2e1 + t190 / 0.2e1) * t296 + (t500 / 0.2e1 + t192 / 0.2e1) * t293 - t320) * t287) * qJD(1); m(4) * t61 + m(5) * t41 + m(6) * t23 + m(7) * t9; (t118 * t60 + t119 * t59 + t42 * t9) * t493 + (t161 * t89 + t162 * t88 + t23 * t50) * t494 + (t114 * t183 + t115 * t182 + t41 * t93) * t495 + t425 * t266 * t243 * t496 + (t34 + t33) * t422 + (t32 + t31) * t421 + (t202 * t519 + t516 * t422 + t511 * t527 - t526 + (-t532 + t549) * t421) * t288 + (t201 * t519 + t512 * t528 + t513 * t421 + ((t511 - t533) * t287 + t512 * t288 + t530 * t416 + t529 * t415 + (-t530 * t293 - t529 * t296) * qJD(3) + ((t503 + t547) * t287 - t558 + t513 + t516) * qJD(1)) * t288 + t525 + (-t287 * t545 + t532) * t422) * t287; 0.2e1 * ((t287 * t69 + t288 * t68) * t490 + (t287 * t99 + t288 * t98) * t491 + (t120 * t288 + t121 * t287) * t492) * t415 + 0.2e1 * ((t287 * t39 + t288 * t40 + t421 * t69 - t422 * t68) * t490 + (t287 * t47 + t288 * t48 + t421 * t99 - t422 * t98) * t491 + (-t120 * t422 + t121 * t421 + t287 * t90 + t288 * t91) * t492) * t293; (m(5) + m(6) + m(7)) * t416; 0.2e1 * ((t118 * t418 + t119 * t417 - t9) * t490 + (t161 * t418 + t162 * t417 - t23) * t491 + (t182 * t418 + t183 * t417 - t41) * t492) * t296 + 0.2e1 * ((qJD(3) * t42 + t118 * t421 - t119 * t422 + t287 * t60 + t288 * t59) * t490 + (qJD(3) * t50 + t161 * t421 - t162 * t422 + t287 * t89 + t288 * t88) * t491 + (qJD(3) * t93 + t114 * t288 + t115 * t287 + t182 * t421 - t183 * t422) * t492) * t293; 0.4e1 * (t492 + t491 + t490) * (-0.1e1 + t425) * t293 * t415; m(6) * (t108 * t48 + t109 * t47 + t43 * t99 + t44 * t98) + m(7) * (t24 * t69 + t25 * t68 + t39 * t67 + t40 * t66) + (t287 * t320 + t288 * t321) * t416 + (t323 * t288 + t322 * t287 + (t287 * t321 - t288 * t320) * qJD(1)) * t296 + t382; m(6) * t30 + m(7) * t10; m(6) * (t108 * t88 + t109 * t89 + t161 * t43 + t162 * t44 - t23 * t92 + t30 * t50) + m(7) * (t10 * t42 + t118 * t24 + t119 * t25 + t49 * t9 + t59 * t66 + t60 * t67) + (-t2 / 0.2e1 - t1 / 0.2e1 + t407 * t416) * t288 + (t3 / 0.2e1 + t4 / 0.2e1 - t406 * t416) * t287 + ((t287 * t407 + t288 * t406) * qJD(1) + t526 * t488 + t525 * t288 / 0.2e1 + (t287 * t477 - t288 * t476) * t524) * t296 + (t531 * qJD(1) + t523 * t287 - t522 * t288) * t293 / 0.2e1 + (t287 * t520 + t288 * t521) * qJD(1) / 0.2e1; 0.2e1 * ((t108 * t417 + t109 * t418 - t30) * t491 + (t417 * t66 + t418 * t67 - t10) * t490) * t296 + 0.2e1 * ((-qJD(3) * t92 - t108 * t422 + t109 * t421 + t287 * t43 + t288 * t44) * t491 + (qJD(3) * t49 + t24 * t287 + t25 * t288 + t421 * t67 - t422 * t66) * t490) * t293; (t10 * t49 + t24 * t67 + t25 * t66) * t493 + (t108 * t44 + t109 * t43 - t30 * t92) * t494 + ((-t287 * t498 + t325 * t288) * qJD(3) + t382) * t293 + ((t293 * t523 + t3 + t4) * t288 + (t293 * t522 + t1 + t2) * t287 + (t540 * t293 + t531 * t296) * qJD(3) + (t325 * t287 + t288 * t498) * qJD(1)) * t296; m(7) * (t129 * t69 + t131 * t68 + t217 * t40 - t219 * t39); -t535 * m(7); m(7) * (-t42 * t390 + t118 * t129 + t119 * t131 + t217 * t59 - t219 * t60 + (t296 * t9 - t416 * t42) * t295); m(7) * (t391 + (t217 * t288 - t219 * t287) * t415 + (0.2e1 * t392 + t129 * t287 + t131 * t288 + (-t217 * t287 - t219 * t288) * qJD(1)) * t293); m(7) * (-t49 * t390 + t129 * t67 + t131 * t66 + t217 * t25 - t219 * t24 + (t10 * t296 - t416 * t49) * t295); (-t129 * t219 + t131 * t217 + (-t293 * t392 - t391) * t295) * t493;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
