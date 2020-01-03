% Calculate time derivative of joint inertia matrix for
% S5RRPRP10
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP10_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP10_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP10_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:26
% EndTime: 2019-12-31 20:09:55
% DurationCPUTime: 18.74s
% Computational Cost: add. (11259->843), mult. (29705->1169), div. (0->0), fcn. (27874->6), ass. (0->391)
t296 = cos(qJ(2));
t293 = sin(qJ(2));
t458 = Icges(4,6) * t293;
t468 = Icges(3,4) * t293;
t538 = -t458 - t468 + (-Icges(3,2) - Icges(4,3)) * t296;
t457 = Icges(4,6) * t296;
t467 = Icges(3,4) * t296;
t537 = -t457 - t467 + (-Icges(3,1) - Icges(4,2)) * t293;
t292 = sin(qJ(4));
t295 = cos(qJ(4));
t462 = Icges(6,4) * t295;
t356 = Icges(6,1) * t292 + t462;
t197 = Icges(6,5) * t293 - t296 * t356;
t464 = Icges(5,4) * t295;
t357 = Icges(5,1) * t292 + t464;
t198 = Icges(5,5) * t293 - t296 * t357;
t524 = -t198 - t197;
t536 = t524 * t292;
t348 = Icges(6,5) * t292 + Icges(6,6) * t295;
t189 = Icges(6,3) * t293 - t296 * t348;
t349 = Icges(5,5) * t292 + Icges(5,6) * t295;
t190 = Icges(5,3) * t293 - t296 * t349;
t535 = t190 + t189;
t463 = Icges(6,4) * t292;
t351 = Icges(6,2) * t295 + t463;
t193 = Icges(6,6) * t293 - t296 * t351;
t465 = Icges(5,4) * t292;
t352 = Icges(5,2) * t295 + t465;
t194 = Icges(5,6) * t293 - t296 * t352;
t534 = t193 + t194;
t533 = t538 * qJD(2);
t532 = t537 * qJD(2);
t291 = -qJ(5) - pkin(7);
t470 = rSges(6,3) - t291;
t294 = sin(qJ(1));
t297 = cos(qJ(1));
t355 = -Icges(3,2) * t293 + t467;
t195 = -Icges(3,6) * t297 + t294 * t355;
t345 = -Icges(4,3) * t293 + t457;
t503 = Icges(4,5) * t297 + t294 * t345;
t531 = -t195 - t503;
t196 = Icges(3,6) * t294 + t297 * t355;
t201 = Icges(4,5) * t294 - t297 * t345;
t530 = t196 - t201;
t359 = Icges(3,1) * t296 - t468;
t199 = -Icges(3,5) * t297 + t294 * t359;
t347 = Icges(4,2) * t296 - t458;
t502 = Icges(4,4) * t297 + t294 * t347;
t529 = t199 + t502;
t200 = Icges(3,5) * t294 + t297 * t359;
t203 = Icges(4,4) * t294 - t297 * t347;
t528 = t200 - t203;
t527 = (-t534 * t295 + t536) * t296 + t535 * t293;
t417 = qJD(4) * t296;
t155 = (Icges(6,2) * t292 - t462) * t417 + (Icges(6,6) * t296 + t293 * t351) * qJD(2);
t156 = (Icges(5,2) * t292 - t464) * t417 + (Icges(5,6) * t296 + t293 * t352) * qJD(2);
t526 = -t156 - t155;
t159 = (-Icges(6,1) * t295 + t463) * t417 + (Icges(6,5) * t296 + t293 * t356) * qJD(2);
t160 = (-Icges(5,1) * t295 + t465) * t417 + (Icges(5,5) * t296 + t293 * t357) * qJD(2);
t525 = -t160 - t159;
t151 = (-Icges(6,5) * t295 + Icges(6,6) * t292) * t417 + (Icges(6,3) * t296 + t293 * t348) * qJD(2);
t152 = (-Icges(5,5) * t295 + Icges(5,6) * t292) * t417 + (Icges(5,3) * t296 + t293 * t349) * qJD(2);
t421 = qJD(2) * t296;
t423 = qJD(2) * t293;
t523 = t534 * (t292 * t417 + t295 * t423) + t535 * t421 - t423 * t536 + (t151 + t152) * t293;
t285 = t297 * pkin(3);
t482 = -pkin(7) - t291;
t485 = pkin(4) * t292;
t307 = t293 * t485 + t296 * t482;
t447 = t294 * t296;
t448 = t294 * t295;
t225 = t292 * t297 + t293 * t448;
t445 = t295 * t297;
t449 = t294 * t292;
t226 = t293 * t449 - t445;
t275 = pkin(4) * t295 + pkin(3);
t505 = -t226 * rSges(6,1) - t225 * rSges(6,2) + t275 * t297;
t440 = rSges(6,3) * t447 + t294 * t307 + t285 - t505;
t388 = t440 * t297;
t444 = t296 * t297;
t240 = t294 * pkin(3) + pkin(7) * t444;
t223 = t293 * t445 - t449;
t450 = t293 * t297;
t408 = t292 * t450;
t224 = t408 + t448;
t508 = t224 * rSges(6,1) + t223 * rSges(6,2) + pkin(4) * t408 + t294 * t275 + t470 * t444;
t441 = -t240 + t508;
t522 = -t294 * t441 + t388;
t420 = qJD(2) * t297;
t398 = t293 * t420;
t425 = qJD(1) * t296;
t520 = t294 * t425 + t398;
t134 = Icges(6,5) * t226 + Icges(6,6) * t225 + Icges(6,3) * t447;
t138 = Icges(6,4) * t226 + Icges(6,2) * t225 + Icges(6,6) * t447;
t142 = Icges(6,1) * t226 + Icges(6,4) * t225 + Icges(6,5) * t447;
t342 = t138 * t295 + t142 * t292;
t62 = t134 * t293 - t296 * t342;
t136 = Icges(5,5) * t226 + Icges(5,6) * t225 + Icges(5,3) * t447;
t140 = Icges(5,4) * t226 + Icges(5,2) * t225 + Icges(5,6) * t447;
t144 = Icges(5,1) * t226 + Icges(5,4) * t225 + Icges(5,5) * t447;
t340 = t140 * t295 + t144 * t292;
t64 = t136 * t293 - t296 * t340;
t478 = t62 + t64;
t133 = Icges(6,5) * t224 + Icges(6,6) * t223 + Icges(6,3) * t444;
t137 = Icges(6,4) * t224 + Icges(6,2) * t223 + Icges(6,6) * t444;
t141 = Icges(6,1) * t224 + Icges(6,4) * t223 + Icges(6,5) * t444;
t343 = t137 * t295 + t141 * t292;
t61 = t133 * t293 - t296 * t343;
t135 = Icges(5,5) * t224 + Icges(5,6) * t223 + Icges(5,3) * t444;
t139 = Icges(5,4) * t224 + Icges(5,2) * t223 + Icges(5,6) * t444;
t143 = Icges(5,1) * t224 + Icges(5,4) * t223 + Icges(5,5) * t444;
t341 = t139 * t295 + t143 * t292;
t63 = t135 * t293 - t296 * t341;
t479 = t61 + t63;
t519 = t294 * t479 - t297 * t478;
t518 = t294 * t478 + t297 * t479;
t517 = qJD(2) / 0.2e1;
t385 = qJD(1) * t293 + qJD(4);
t396 = t296 * t420;
t299 = -t294 * t385 + t396;
t386 = qJD(4) * t293 + qJD(1);
t328 = t386 * t292;
t125 = t295 * t299 - t297 * t328;
t327 = t386 * t295;
t126 = t292 * t299 + t297 * t327;
t70 = Icges(6,5) * t126 + Icges(6,6) * t125 - Icges(6,3) * t520;
t74 = Icges(6,4) * t126 + Icges(6,2) * t125 - Icges(6,6) * t520;
t78 = Icges(6,1) * t126 + Icges(6,4) * t125 - Icges(6,5) * t520;
t19 = (qJD(2) * t343 + t70) * t293 + (qJD(2) * t133 - t292 * t78 - t295 * t74 + (t137 * t292 - t141 * t295) * qJD(4)) * t296;
t72 = Icges(5,5) * t126 + Icges(5,6) * t125 - Icges(5,3) * t520;
t76 = Icges(5,4) * t126 + Icges(5,2) * t125 - Icges(5,6) * t520;
t80 = Icges(5,1) * t126 + Icges(5,4) * t125 - Icges(5,5) * t520;
t21 = (qJD(2) * t341 + t72) * t293 + (qJD(2) * t135 - t292 * t80 - t295 * t76 + (t139 * t292 - t143 * t295) * qJD(4)) * t296;
t516 = t19 + t21;
t326 = t385 * t297;
t123 = t295 * t326 + (t295 * t421 - t328) * t294;
t397 = t294 * t421;
t124 = t294 * t327 + (t326 + t397) * t292;
t422 = qJD(2) * t294;
t400 = t293 * t422;
t424 = qJD(1) * t297;
t402 = t296 * t424;
t306 = -t400 + t402;
t69 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t306;
t73 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t306;
t77 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t306;
t20 = (qJD(2) * t342 + t69) * t293 + (qJD(2) * t134 - t292 * t77 - t295 * t73 + (t138 * t292 - t142 * t295) * qJD(4)) * t296;
t71 = Icges(5,5) * t124 + Icges(5,6) * t123 + Icges(5,3) * t306;
t75 = Icges(5,4) * t124 + Icges(5,2) * t123 + Icges(5,6) * t306;
t79 = Icges(5,1) * t124 + Icges(5,4) * t123 + Icges(5,5) * t306;
t22 = (qJD(2) * t340 + t71) * t293 + (qJD(2) * t136 - t292 * t79 - t295 * t75 + (t140 * t292 - t144 * t295) * qJD(4)) * t296;
t515 = t20 + t22;
t54 = t135 * t444 + t223 * t139 + t224 * t143;
t55 = t136 * t444 + t223 * t140 + t224 * t144;
t365 = t294 * t55 + t297 * t54;
t52 = t133 * t444 + t223 * t137 + t224 * t141;
t53 = t134 * t444 + t223 * t138 + t224 * t142;
t366 = t294 * t53 + t297 * t52;
t87 = t189 * t444 + t223 * t193 + t224 * t197;
t88 = t190 * t444 + t223 * t194 + t224 * t198;
t514 = (t365 + t366) * t296 + (t87 + t88) * t293;
t58 = t135 * t447 + t139 * t225 + t143 * t226;
t59 = t136 * t447 + t140 * t225 + t144 * t226;
t363 = t294 * t59 + t297 * t58;
t56 = t133 * t447 + t137 * t225 + t141 * t226;
t57 = t134 * t447 + t138 * t225 + t142 * t226;
t364 = t294 * t57 + t297 * t56;
t89 = t189 * t447 + t193 * t225 + t197 * t226;
t90 = t190 * t447 + t194 * t225 + t198 * t226;
t480 = (t363 + t364) * t296 + (t89 + t90) * t293;
t331 = t201 * t293 - t203 * t296;
t513 = t294 * t331;
t332 = t196 * t293 - t200 * t296;
t512 = t294 * t332;
t330 = -t293 * t503 + t296 * t502;
t510 = t297 * t330;
t333 = t195 * t293 - t199 * t296;
t509 = t297 * t333;
t507 = t294 * rSges(4,1) - rSges(4,2) * t444;
t418 = qJD(4) * t295;
t394 = t293 * t418;
t416 = qJD(5) * t296;
t506 = t275 * t424 + t396 * t485 + (pkin(4) * t394 + t416) * t297 + t520 * t291 + t126 * rSges(6,1) + t125 * rSges(6,2);
t504 = -rSges(3,2) * t450 + t294 * rSges(3,3);
t350 = Icges(3,5) * t296 - Icges(3,6) * t293;
t191 = -Icges(3,3) * t297 + t294 * t350;
t353 = Icges(4,4) * t296 - Icges(4,5) * t293;
t501 = Icges(4,1) * t297 + t294 * t353;
t500 = -t294 * t440 - t297 * t441;
t499 = 2 * m(3);
t498 = 2 * m(4);
t497 = 2 * m(5);
t496 = 2 * m(6);
t495 = m(4) / 0.2e1;
t494 = m(5) / 0.2e1;
t493 = m(6) / 0.2e1;
t491 = t294 / 0.2e1;
t490 = -t297 / 0.2e1;
t488 = -rSges(6,3) - pkin(2);
t255 = rSges(3,1) * t293 + rSges(3,2) * t296;
t487 = m(3) * t255;
t486 = pkin(2) * t296;
t263 = pkin(7) * t400;
t373 = t124 * rSges(6,1) + t123 * rSges(6,2);
t411 = qJD(4) * t485;
t477 = rSges(6,3) * t306 + t373 + t263 + (qJD(1) * t307 + t411) * t297 + (t291 * t423 + t416 + (-pkin(3) + t275) * qJD(1) + (t292 * t421 + t394) * pkin(4)) * t294;
t279 = pkin(3) * t424;
t476 = -rSges(6,3) * t520 + pkin(7) * t398 - t279 + (pkin(7) * t425 - t385 * t485) * t294 + t506;
t475 = rSges(4,1) * t297;
t474 = rSges(4,2) * t293;
t473 = rSges(3,3) * t297;
t472 = rSges(5,3) * t293;
t471 = -rSges(4,3) - qJ(3);
t455 = qJ(3) * t293;
t454 = qJ(3) * t296;
t375 = -rSges(5,1) * t226 - rSges(5,2) * t225;
t148 = rSges(5,3) * t447 - t375;
t453 = t148 * t297;
t451 = t292 * t296;
t442 = t126 * rSges(5,1) + t125 * rSges(5,2);
t371 = rSges(6,1) * t292 + rSges(6,2) * t295;
t393 = t295 * t417;
t439 = (-rSges(6,1) * t295 + rSges(6,2) * t292) * t417 - pkin(4) * t393 + qJD(5) * t293 + (rSges(6,3) * t296 + t293 * t371 + t307) * qJD(2);
t438 = -pkin(4) * t451 - t296 * t371 + (rSges(6,3) + t482) * t293;
t368 = t455 + t486;
t229 = t368 * t294;
t230 = pkin(2) * t444 + qJ(3) * t450;
t437 = t294 * t229 + t297 * t230;
t217 = qJD(2) * t368 - qJD(3) * t296;
t370 = -rSges(4,2) * t296 + rSges(4,3) * t293;
t436 = -t370 * qJD(2) - t217;
t435 = -t230 - t240;
t253 = pkin(2) * t293 - t454;
t426 = qJD(1) * t294;
t231 = t253 * t426;
t405 = t293 * t426;
t434 = pkin(7) * t405 + t231;
t369 = rSges(4,3) * t296 + t474;
t433 = -t253 + t369;
t419 = qJD(3) * t293;
t432 = qJ(3) * t396 + t297 * t419;
t431 = rSges(3,2) * t405 + rSges(3,3) * t424;
t430 = t297 * pkin(1) + t294 * pkin(6);
t429 = t294 ^ 2 + t297 ^ 2;
t192 = Icges(3,3) * t294 + t297 * t350;
t428 = qJD(1) * t192;
t205 = Icges(4,1) * t294 - t297 * t353;
t427 = qJD(1) * t205;
t415 = -rSges(5,3) - pkin(2) - pkin(7);
t414 = -0.1e1 + t429;
t35 = t52 * t294 - t297 * t53;
t36 = t54 * t294 - t297 * t55;
t410 = -t35 / 0.2e1 - t36 / 0.2e1;
t37 = t56 * t294 - t297 * t57;
t38 = t58 * t294 - t297 * t59;
t409 = t37 / 0.2e1 + t38 / 0.2e1;
t264 = pkin(2) * t400;
t407 = t294 * (pkin(2) * t402 + t294 * t419 - t264 + (t293 * t424 + t397) * qJ(3)) + t297 * (-pkin(2) * t520 - qJ(3) * t405 + t432) + t229 * t424;
t146 = t224 * rSges(5,1) + t223 * rSges(5,2) + rSges(5,3) * t444;
t278 = pkin(6) * t424;
t406 = t278 + t432;
t374 = rSges(5,1) * t292 + rSges(5,2) * t295;
t210 = -t296 * t374 + t472;
t404 = t210 * t426;
t392 = -t353 * qJD(2) / 0.2e1 + t350 * t517;
t391 = -qJ(3) - t485;
t390 = -pkin(7) * t293 - t253;
t389 = t527 * t421 + (((qJD(4) * t524 + t526) * t295 + t525 * t292) * t296 + t523) * t293;
t184 = t433 * t297;
t387 = qJD(1) * t438;
t241 = pkin(7) * t447 - t285;
t384 = t297 * t240 + t294 * t241 + t437;
t383 = rSges(4,1) * t424 + rSges(4,2) * t520 + rSges(4,3) * t396;
t382 = t430 + t230;
t381 = -t210 + t390;
t380 = -pkin(7) * t421 - t217;
t379 = t294 * t387;
t378 = t293 * t471 - pkin(1);
t377 = rSges(3,1) * t296 - rSges(3,2) * t293;
t376 = t124 * rSges(5,1) + t123 * rSges(5,2);
t316 = t293 * t391 - pkin(1);
t44 = t488 * t398 + (-t411 + (t296 * t488 + t316) * qJD(1)) * t294 + t406 + t506;
t298 = (-pkin(2) - t470) * t296 + t316;
t45 = t264 + (qJD(1) * t298 - t411) * t297 + (-t416 + (-pkin(4) * t418 - qJD(3)) * t293 + (-pkin(6) - t275) * qJD(1) + (t293 * t470 + t296 * t391) * qJD(2)) * t294 - t373;
t367 = t294 * t44 + t297 * t45;
t284 = t297 * pkin(6);
t94 = t294 * t298 + t284 + t505;
t95 = t382 + t508;
t362 = t294 * t95 + t297 * t94;
t339 = -t146 * t297 - t294 * t148;
t338 = t146 * t294 - t453;
t329 = t390 - t438;
t211 = rSges(3,1) * t444 + t504;
t212 = rSges(4,3) * t450 + t507;
t170 = (-rSges(5,1) * t295 + rSges(5,2) * t292) * t417 + (rSges(5,3) * t296 + t293 * t374) * qJD(2);
t325 = -t170 + t380;
t324 = -pkin(1) - t377;
t33 = t125 * t193 + t126 * t197 + t151 * t444 + t223 * t155 + t224 * t159 - t189 * t520;
t34 = t125 * t194 + t126 * t198 + t152 * t444 + t223 * t156 + t224 * t160 - t190 * t520;
t322 = t21 / 0.2e1 + t19 / 0.2e1 + t33 / 0.2e1 + t34 / 0.2e1;
t31 = t123 * t193 + t124 * t197 + t151 * t447 + t225 * t155 + t226 * t159 + t189 * t306;
t32 = t123 * t194 + t124 * t198 + t152 * t447 + t225 * t156 + t226 * t160 + t190 * t306;
t321 = t22 / 0.2e1 + t20 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1;
t320 = -t63 / 0.2e1 - t61 / 0.2e1 - t87 / 0.2e1 - t88 / 0.2e1;
t319 = t64 / 0.2e1 + t62 / 0.2e1 + t89 / 0.2e1 + t90 / 0.2e1;
t132 = t381 * t297;
t318 = t294 * (t240 * qJD(1) - t263) + t297 * (-pkin(7) * t520 + t279) + t241 * t424 + t407;
t315 = t380 - t439;
t314 = qJD(2) * t255;
t311 = qJD(2) * (Icges(4,4) * t293 + Icges(4,5) * t296);
t310 = qJD(2) * (-Icges(3,5) * t293 - Icges(3,6) * t296);
t116 = t329 * t297;
t24 = (t420 * t438 + t476) * t293 + (qJD(2) * t441 - t297 * t439 + t379) * t296;
t25 = (-t422 * t438 - t477) * t293 + (-qJD(2) * t440 + t294 * t439 + t297 * t387) * t296;
t49 = t522 * t296;
t304 = qJD(2) * t49 + t24 * t294 + t25 * t297;
t46 = t384 - t500;
t50 = t297 * t315 + t379 + t434;
t51 = qJD(1) * t116 + t294 * t315;
t303 = qJD(2) * t46 + t294 * t51 + t297 * t50;
t302 = t296 * t415 - pkin(1) - t455;
t301 = (rSges(4,2) - pkin(2)) * t296 + t378;
t300 = t302 * t294;
t239 = t377 * qJD(2);
t213 = t294 * t370 - t475;
t208 = t294 * t377 - t473;
t183 = t433 * t294;
t176 = t211 + t430;
t175 = t294 * t324 + t284 + t473;
t174 = t414 * t293 * t421;
t168 = t501 * qJD(1) + t297 * t311;
t167 = t294 * t311 + t427;
t154 = t294 * t310 + t428;
t153 = -qJD(1) * t191 + t297 * t310;
t150 = t212 + t382;
t149 = t294 * t301 + t284 + t475;
t131 = t381 * t294;
t118 = t255 * t422 + ((-rSges(3,3) - pkin(6)) * t294 + t324 * t297) * qJD(1);
t117 = -rSges(3,1) * t520 - rSges(3,2) * t396 - pkin(1) * t426 + t278 + t431;
t115 = t329 * t294;
t114 = qJD(1) * t184 + t294 * t436;
t113 = t297 * t436 - t369 * t426 + t231;
t112 = t293 * t146 - t210 * t444;
t111 = -t148 * t293 + t210 * t447;
t110 = t294 * t330 + t297 * t501;
t109 = -t205 * t297 + t513;
t108 = -t294 * t501 + t510;
t107 = t294 * t205 + t331 * t297;
t106 = t294 * t192 - t332 * t297;
t105 = t294 * t191 - t509;
t102 = -t192 * t297 - t512;
t101 = -t191 * t297 - t294 * t333;
t100 = t382 + t146 + t240;
t99 = t284 + t285 + t300 + t375;
t96 = t212 * t297 + t294 * t213 + t437;
t93 = t264 + (-t419 + (t296 * t471 - t474) * qJD(2)) * t294 + ((-rSges(4,1) - pkin(6)) * t294 + t301 * t297) * qJD(1);
t92 = -pkin(2) * t398 + (t378 - t486) * t426 + t383 + t406;
t91 = t338 * t296;
t84 = -rSges(5,3) * t520 + t442;
t82 = rSges(5,3) * t306 + t376;
t68 = qJD(1) * t132 + t294 * t325;
t67 = t297 * t325 + t404 + t434;
t66 = t293 * t441 - t438 * t444;
t65 = -t293 * t440 + t438 * t447;
t60 = -t339 + t384;
t48 = t263 + t264 + (-t419 + (-t454 + t472) * qJD(2)) * t294 + ((-pkin(3) - pkin(6)) * t294 + t302 * t297) * qJD(1) - t376;
t47 = qJD(1) * t300 + t398 * t415 + t279 + t406 + t442;
t43 = (-t210 * t422 - t82) * t293 + (-qJD(2) * t148 + t294 * t170 + t210 * t424) * t296;
t42 = (t210 * t420 + t84) * t293 + (qJD(2) * t146 - t170 * t297 + t404) * t296;
t41 = (qJD(1) * t213 + t383) * t297 + (t369 * t422 + (-t212 - t230 + t507) * qJD(1)) * t294 + t407;
t30 = t338 * t423 + (qJD(1) * t339 - t294 * t84 + t297 * t82) * t296;
t23 = t294 * t82 + t297 * t84 + (t453 + (-t146 + t435) * t294) * qJD(1) + t318;
t18 = t125 * t140 + t126 * t144 - t136 * t520 + t223 * t75 + t224 * t79 + t444 * t71;
t17 = t125 * t139 + t126 * t143 - t135 * t520 + t223 * t76 + t224 * t80 + t444 * t72;
t16 = t125 * t138 + t126 * t142 - t134 * t520 + t223 * t73 + t224 * t77 + t444 * t69;
t15 = t125 * t137 + t126 * t141 - t133 * t520 + t223 * t74 + t224 * t78 + t444 * t70;
t14 = t123 * t140 + t124 * t144 + t136 * t306 + t225 * t75 + t226 * t79 + t447 * t71;
t13 = t123 * t139 + t124 * t143 + t135 * t306 + t225 * t76 + t226 * t80 + t447 * t72;
t12 = t123 * t138 + t124 * t142 + t134 * t306 + t225 * t73 + t226 * t77 + t447 * t69;
t11 = t123 * t137 + t124 * t141 + t133 * t306 + t225 * t74 + t226 * t78 + t447 * t70;
t10 = -t522 * t423 + (t500 * qJD(1) - t476 * t294 + t477 * t297) * t296;
t9 = t476 * t297 + t477 * t294 + (t388 + (t435 - t441) * t294) * qJD(1) + t318;
t8 = qJD(1) * t365 + t17 * t294 - t18 * t297;
t7 = qJD(1) * t366 + t15 * t294 - t16 * t297;
t6 = qJD(1) * t363 + t13 * t294 - t14 * t297;
t5 = qJD(1) * t364 + t11 * t294 - t12 * t297;
t4 = (-qJD(2) * t365 + t34) * t293 + (-qJD(1) * t36 + qJD(2) * t88 + t17 * t297 + t18 * t294) * t296;
t3 = (-qJD(2) * t366 + t33) * t293 + (-qJD(1) * t35 + qJD(2) * t87 + t15 * t297 + t16 * t294) * t296;
t2 = (-qJD(2) * t363 + t32) * t293 + (-qJD(1) * t38 + qJD(2) * t90 + t13 * t297 + t14 * t294) * t296;
t1 = (-qJD(2) * t364 + t31) * t293 + (-qJD(1) * t37 + qJD(2) * t89 + t11 * t297 + t12 * t294) * t296;
t26 = [(t117 * t176 + t118 * t175) * t499 + (t149 * t93 + t150 * t92) * t498 + (t44 * t95 + t45 * t94) * t496 + (t100 * t47 + t48 * t99) * t497 + t525 * t451 + t526 * t295 * t296 + t524 * t393 + (t359 + t347 + t538) * t423 + (t345 + t355 - t537) * t421 + t523; m(4) * (t113 * t149 + t114 * t150 + t183 * t92 + t184 * t93) + m(5) * (t100 * t68 + t131 * t47 + t132 * t48 + t67 * t99) + m(6) * (t115 * t44 + t116 * t45 + t50 * t94 + t51 * t95) + (m(3) * (-t118 * t255 - t175 * t239) + t392 * t297 - t321) * t297 + (m(3) * (-t117 * t255 - t176 * t239) + t392 * t294 + t322) * t294 + ((-t530 * qJD(2) + t532 * t297) * t491 + (t531 * qJD(2) + t532 * t294) * t490 + (t528 * t490 - t529 * t491) * qJD(1)) * t293 + ((t528 * qJD(2) + t533 * t297) * t491 + (t529 * qJD(2) + t533 * t294) * t490 + (t530 * t490 + t531 * t491) * qJD(1)) * t296 + ((-t176 * t487 + (t196 / 0.2e1 - t201 / 0.2e1) * t296 + (t200 / 0.2e1 - t203 / 0.2e1) * t293 - t320) * t297 + (t175 * t487 + (t503 / 0.2e1 + t195 / 0.2e1) * t296 + (t502 / 0.2e1 + t199 / 0.2e1) * t293 + t319) * t294) * qJD(1); t294 * t7 + (t115 * t51 + t116 * t50 + t46 * t9) * t496 + (t131 * t68 + t132 * t67 + t60 * t23) * t497 - t297 * t6 + t294 * t8 - t297 * t5 + (t113 * t184 + t114 * t183 + t41 * t96) * t498 + ((t294 * t208 + t211 * t297) * ((qJD(1) * t208 - t297 * t314 + t431) * t297 + (-t294 * t314 + (-t211 + t504) * qJD(1)) * t294) + t429 * t255 * t239) * t499 + t294 * ((t294 * t153 + (t105 + t512) * qJD(1)) * t294 + (t106 * qJD(1) + (t195 * t421 + t199 * t423) * t297 + (-t154 + (-t196 * t296 - t200 * t293) * qJD(2) + (t192 - t333) * qJD(1)) * t294) * t297) - t297 * ((t154 * t297 + (t102 + t509) * qJD(1)) * t297 + (t101 * qJD(1) + (-t196 * t421 - t200 * t423 + t428) * t294 + (-t153 + (t195 * t296 + t199 * t293) * qJD(2) - t332 * qJD(1)) * t297) * t294) + t294 * ((t294 * t168 + (t108 - t513) * qJD(1)) * t294 + (t107 * qJD(1) + (t421 * t503 + t423 * t502) * t297 + (-t167 + (t201 * t296 + t203 * t293) * qJD(2) + (t205 + t330) * qJD(1)) * t294) * t297) - t297 * ((t167 * t297 + (t109 - t510) * qJD(1)) * t297 + (t110 * qJD(1) + (t201 * t421 + t203 * t423 + t427) * t294 + (-t168 + (t293 * t502 + t296 * t503) * qJD(2) + t331 * qJD(1)) * t297) * t294) + (t37 + t38 + (-t101 - t110) * t297 + (t102 + t109) * t294) * t426 + (t35 + t36 + (-t105 - t108) * t297 + (t106 + t107) * t294) * t424; 0.2e1 * (t362 * t493 + (t100 * t294 + t297 * t99) * t494 + (t149 * t297 + t150 * t294) * t495) * t421 + 0.2e1 * ((t424 * t95 - t426 * t94 + t367) * t493 + (t100 * t424 + t294 * t47 + t297 * t48 - t426 * t99) * t494 + (-t149 * t426 + t150 * t424 + t294 * t92 + t297 * t93) * t495) * t293; 0.2e1 * ((t115 * t422 + t116 * t420 - t9) * t493 + (t131 * t422 + t132 * t420 - t23) * t494 + (t183 * t422 + t184 * t420 - t41) * t495) * t296 + 0.2e1 * ((t115 * t424 - t116 * t426 + t303) * t493 + (qJD(2) * t60 + t131 * t424 - t132 * t426 + t294 * t68 + t297 * t67) * t494 + (qJD(2) * t96 + t113 * t297 + t114 * t294 + t183 * t424 - t184 * t426) * t495) * t293; 0.4e1 * (t495 + t494 + t493) * t174; m(5) * (t100 * t42 + t111 * t48 + t112 * t47 + t43 * t99) + m(6) * (t24 * t95 + t25 * t94 + t44 * t66 + t45 * t65) + (-t294 * t319 + t297 * t320) * t423 + (t322 * t297 + t321 * t294 + (t294 * t320 + t297 * t319) * qJD(1)) * t296 + t389; m(5) * (t111 * t67 + t112 * t68 + t131 * t42 + t132 * t43 - t23 * t91 + t30 * t60) + m(6) * (t10 * t46 + t115 * t24 + t116 * t25 + t49 * t9 + t50 * t65 + t51 * t66) + (-t2 / 0.2e1 - t1 / 0.2e1 + t410 * t423) * t297 + (t4 / 0.2e1 + t3 / 0.2e1 - t409 * t423) * t294 + ((t294 * t410 + t297 * t409) * qJD(1) + (t5 + t6) * t491 + (t7 + t8) * t297 / 0.2e1 + t519 * t517) * t296 + (qJD(1) * t518 + t516 * t294 - t515 * t297) * t293 / 0.2e1 + (t480 * t294 + t514 * t297) * qJD(1) / 0.2e1; 0.2e1 * ((t111 * t420 + t112 * t422 - t30) * t494 + (t420 * t65 + t422 * t66 - t10) * t493) * t296 + 0.2e1 * ((-qJD(2) * t91 - t111 * t426 + t112 * t424 + t294 * t42 + t297 * t43) * t494 + (t424 * t66 - t426 * t65 + t304) * t493) * t293; (t49 * t10 + t24 * t66 + t25 * t65) * t496 + (t111 * t43 + t112 * t42 - t30 * t91) * t497 + (((-t293 * t479 - t514) * t297 + (-t293 * t478 - t480) * t294) * qJD(2) + t389) * t293 + ((t3 + t4) * t297 + (t1 + t2) * t294 + (t515 * t294 + t516 * t297) * t293 + (t293 * t527 + t518 * t296) * qJD(2) + (-t293 * t519 - t294 * t514 + t480 * t297) * qJD(1)) * t296; m(6) * (-t362 * t423 + ((-t294 * t94 + t297 * t95) * qJD(1) + t367) * t296); m(6) * ((t9 + (-t115 * t294 - t116 * t297) * qJD(2)) * t293 + ((t115 * t297 - t116 * t294) * qJD(1) + t303) * t296); m(6) * (-t293 ^ 2 + t296 ^ 2) * t414 * qJD(2); m(6) * ((t10 + (-t294 * t66 - t297 * t65) * qJD(2)) * t293 + ((-t294 * t65 + t297 * t66) * qJD(1) + t304) * t296); -0.2e1 * m(6) * t174;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t26(1), t26(2), t26(4), t26(7), t26(11); t26(2), t26(3), t26(5), t26(8), t26(12); t26(4), t26(5), t26(6), t26(9), t26(13); t26(7), t26(8), t26(9), t26(10), t26(14); t26(11), t26(12), t26(13), t26(14), t26(15);];
Mq = res;
