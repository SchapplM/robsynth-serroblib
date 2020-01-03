% Calculate time derivative of joint inertia matrix for
% S5RRPRP9
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP9_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP9_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP9_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:43
% EndTime: 2019-12-31 20:06:12
% DurationCPUTime: 16.60s
% Computational Cost: add. (20371->899), mult. (31612->1273), div. (0->0), fcn. (30295->8), ass. (0->413)
t311 = sin(qJ(2));
t313 = cos(qJ(2));
t304 = pkin(8) + qJ(4);
t298 = sin(t304);
t299 = cos(t304);
t358 = Icges(5,5) * t299 - Icges(5,6) * t298;
t215 = -Icges(5,3) * t313 + t311 * t358;
t361 = Icges(6,4) * t299 + Icges(6,6) * t298;
t216 = -Icges(6,2) * t313 + t311 * t361;
t532 = t215 + t216;
t314 = cos(qJ(1));
t312 = sin(qJ(1));
t458 = t312 * t313;
t232 = t298 * t458 + t299 * t314;
t467 = t298 * t314;
t233 = t299 * t458 - t467;
t527 = rSges(6,3) + qJ(5);
t529 = rSges(6,1) + pkin(4);
t531 = t527 * t232 + t529 * t233;
t431 = qJD(2) * t312;
t403 = t311 * t431;
t426 = qJD(4) * t312;
t433 = qJD(1) * t314;
t435 = qJD(1) * t312;
t139 = -t298 * t403 - qJD(4) * t467 - t299 * t435 + (t298 * t433 + t299 * t426) * t313;
t425 = qJD(4) * t313;
t346 = t298 * (qJD(1) - t425);
t434 = qJD(1) * t313;
t389 = -qJD(4) + t434;
t140 = t312 * t346 + (t314 * t389 - t403) * t299;
t530 = -t232 * qJD(5) - t527 * t139 - t529 * t140;
t500 = t312 / 0.2e1;
t498 = -t314 / 0.2e1;
t528 = -qJD(1) / 0.2e1;
t473 = Icges(6,5) * t299;
t357 = Icges(6,3) * t298 + t473;
t214 = -Icges(6,6) * t313 + t311 * t357;
t476 = Icges(5,4) * t299;
t362 = -Icges(5,2) * t298 + t476;
t217 = -Icges(5,6) * t313 + t311 * t362;
t474 = Icges(6,5) * t298;
t366 = Icges(6,1) * t299 + t474;
t218 = -Icges(6,4) * t313 + t311 * t366;
t477 = Icges(5,4) * t298;
t367 = Icges(5,1) * t299 - t477;
t219 = -Icges(5,5) * t313 + t311 * t367;
t526 = t532 * t313 + ((-t218 - t219) * t299 + (-t214 + t217) * t298) * t311;
t309 = cos(pkin(8));
t295 = pkin(3) * t309 + pkin(2);
t491 = pkin(2) - t295;
t395 = t491 * t313;
t470 = qJ(3) * t311;
t525 = t395 + t470;
t430 = qJD(2) * t313;
t400 = t312 * t430;
t322 = t311 * t433 + t400;
t457 = t313 * t314;
t234 = t298 * t457 - t312 * t299;
t235 = t312 * t298 + t299 * t457;
t460 = t311 * t314;
t149 = Icges(6,4) * t235 + Icges(6,2) * t460 + Icges(6,6) * t234;
t145 = Icges(6,5) * t235 + Icges(6,6) * t460 + Icges(6,3) * t234;
t153 = Icges(6,1) * t235 + Icges(6,4) * t460 + Icges(6,5) * t234;
t355 = t145 * t298 + t153 * t299;
t64 = -t149 * t313 + t311 * t355;
t147 = Icges(5,5) * t235 - Icges(5,6) * t234 + Icges(5,3) * t460;
t151 = Icges(5,4) * t235 - Icges(5,2) * t234 + Icges(5,6) * t460;
t155 = Icges(5,1) * t235 - Icges(5,4) * t234 + Icges(5,5) * t460;
t353 = -t151 * t298 + t155 * t299;
t66 = -t147 * t313 + t311 * t353;
t487 = t64 + t66;
t461 = t311 * t312;
t148 = Icges(6,4) * t233 + Icges(6,2) * t461 + Icges(6,6) * t232;
t144 = Icges(6,5) * t233 + Icges(6,6) * t461 + Icges(6,3) * t232;
t152 = Icges(6,1) * t233 + Icges(6,4) * t461 + Icges(6,5) * t232;
t356 = t144 * t298 + t152 * t299;
t63 = -t148 * t313 + t311 * t356;
t146 = Icges(5,5) * t233 - Icges(5,6) * t232 + Icges(5,3) * t461;
t150 = Icges(5,4) * t233 - Icges(5,2) * t232 + Icges(5,6) * t461;
t154 = Icges(5,1) * t233 - Icges(5,4) * t232 + Icges(5,5) * t461;
t354 = -t150 * t298 + t154 * t299;
t65 = -t146 * t313 + t311 * t354;
t488 = t63 + t65;
t524 = t312 * t488 + t314 * t487;
t76 = Icges(6,5) * t140 + Icges(6,6) * t322 + Icges(6,3) * t139;
t80 = Icges(6,4) * t140 + Icges(6,2) * t322 + Icges(6,6) * t139;
t84 = Icges(6,1) * t140 + Icges(6,4) * t322 + Icges(6,5) * t139;
t20 = (qJD(2) * t356 - t80) * t313 + (qJD(2) * t148 + t298 * t76 + t299 * t84 + (t144 * t299 - t152 * t298) * qJD(4)) * t311;
t78 = Icges(5,5) * t140 - Icges(5,6) * t139 + Icges(5,3) * t322;
t82 = Icges(5,4) * t140 - Icges(5,2) * t139 + Icges(5,6) * t322;
t86 = Icges(5,1) * t140 - Icges(5,4) * t139 + Icges(5,5) * t322;
t22 = (qJD(2) * t354 - t78) * t313 + (qJD(2) * t146 - t298 * t82 + t299 * t86 + (-t150 * t299 - t154 * t298) * qJD(4)) * t311;
t523 = -t20 - t22;
t397 = t299 * t425;
t429 = qJD(2) * t314;
t401 = t311 * t429;
t137 = qJD(1) * t232 - t314 * t397 + (t401 - t426) * t298;
t138 = t314 * t346 + (-t312 * t389 - t401) * t299;
t396 = t313 * t429;
t406 = t311 * t435;
t321 = t396 - t406;
t75 = Icges(6,5) * t138 + Icges(6,6) * t321 - Icges(6,3) * t137;
t79 = Icges(6,4) * t138 + Icges(6,2) * t321 - Icges(6,6) * t137;
t83 = Icges(6,1) * t138 + Icges(6,4) * t321 - Icges(6,5) * t137;
t21 = (qJD(2) * t355 - t79) * t313 + (qJD(2) * t149 + t298 * t75 + t299 * t83 + (t145 * t299 - t153 * t298) * qJD(4)) * t311;
t77 = Icges(5,5) * t138 + Icges(5,6) * t137 + Icges(5,3) * t321;
t81 = Icges(5,4) * t138 + Icges(5,2) * t137 + Icges(5,6) * t321;
t85 = Icges(5,1) * t138 + Icges(5,4) * t137 + Icges(5,5) * t321;
t23 = (qJD(2) * t353 - t77) * t313 + (qJD(2) * t147 - t298 * t81 + t299 * t85 + (-t151 * t299 - t155 * t298) * qJD(4)) * t311;
t522 = t21 + t23;
t55 = t146 * t461 - t150 * t232 + t154 * t233;
t56 = t147 * t461 - t151 * t232 + t155 * t233;
t373 = t312 * t55 + t314 * t56;
t53 = t144 * t232 + t148 * t461 + t152 * t233;
t54 = t145 * t232 + t149 * t461 + t153 * t233;
t374 = t312 * t53 + t314 * t54;
t96 = t214 * t232 + t216 * t461 + t218 * t233;
t97 = t215 * t461 - t217 * t232 + t219 * t233;
t521 = (-t96 - t97) * t313 + (t373 + t374) * t311;
t59 = t146 * t460 - t234 * t150 + t235 * t154;
t60 = t147 * t460 - t234 * t151 + t235 * t155;
t371 = t312 * t59 + t314 * t60;
t57 = t234 * t144 + t148 * t460 + t235 * t152;
t58 = t234 * t145 + t149 * t460 + t235 * t153;
t372 = t312 * t57 + t314 * t58;
t98 = t234 * t214 + t216 * t460 + t235 * t218;
t99 = t215 * t460 - t234 * t217 + t235 * t219;
t520 = (-t98 - t99) * t313 + (t371 + t372) * t311;
t427 = qJD(4) * t311;
t160 = (Icges(6,3) * t299 - t474) * t427 + (Icges(6,6) * t311 + t313 * t357) * qJD(2);
t161 = (-Icges(5,5) * t298 - Icges(5,6) * t299) * t427 + (Icges(5,3) * t311 + t313 * t358) * qJD(2);
t162 = (-Icges(6,4) * t298 + Icges(6,6) * t299) * t427 + (Icges(6,2) * t311 + t313 * t361) * qJD(2);
t164 = (-Icges(6,1) * t298 + t473) * t427 + (Icges(6,4) * t311 + t313 * t366) * qJD(2);
t165 = (-Icges(5,1) * t298 - t476) * t427 + (Icges(5,5) * t311 + t313 * t367) * qJD(2);
t398 = t299 * t427;
t399 = t298 * t427;
t404 = t299 * t430;
t405 = t298 * t430;
t432 = qJD(2) * t311;
t468 = t298 * t311;
t519 = -t217 * t398 + t219 * t404 + t160 * t468 + (-t399 + t404) * t218 + (t398 + t405) * t214 + (t165 + t164) * t299 * t311 + t532 * t432 + (-t161 - t162) * t313;
t518 = rSges(6,2) * t396 + t234 * qJD(5) - t527 * t137 + t529 * t138;
t517 = t311 * t491;
t478 = Icges(3,4) * t313;
t365 = -Icges(3,2) * t311 + t478;
t239 = Icges(3,6) * t312 + t314 * t365;
t479 = Icges(3,4) * t311;
t370 = Icges(3,1) * t313 - t479;
t241 = Icges(3,5) * t312 + t314 * t370;
t347 = t239 * t311 - t241 * t313;
t330 = t347 * t312;
t238 = -Icges(3,6) * t314 + t312 * t365;
t240 = -Icges(3,5) * t314 + t312 * t370;
t348 = t238 * t311 - t240 * t313;
t331 = t348 * t314;
t516 = -rSges(3,2) * t460 + t312 * rSges(3,3);
t360 = Icges(3,5) * t313 - Icges(3,6) * t311;
t236 = -Icges(3,3) * t314 + t312 * t360;
t515 = t313 * t487 - t520;
t452 = rSges(6,2) * t460 + t527 * t234 + t529 * t235;
t453 = rSges(6,2) * t461 + t531;
t514 = -t312 * t453 - t314 * t452;
t513 = 2 * m(3);
t512 = 2 * m(4);
t511 = 2 * m(5);
t510 = 2 * m(6);
t509 = m(4) / 0.2e1;
t508 = m(5) / 0.2e1;
t507 = m(6) / 0.2e1;
t308 = sin(pkin(8));
t363 = Icges(4,4) * t309 - Icges(4,2) * t308;
t229 = -Icges(4,6) * t313 + t311 * t363;
t506 = t229 / 0.2e1;
t368 = Icges(4,1) * t309 - Icges(4,4) * t308;
t230 = -Icges(4,5) * t313 + t311 * t368;
t505 = t230 / 0.2e1;
t253 = -t308 * t457 + t312 * t309;
t504 = t253 / 0.2e1;
t459 = t312 * t308;
t254 = t309 * t457 + t459;
t503 = t254 / 0.2e1;
t502 = -t308 / 0.2e1;
t501 = t309 / 0.2e1;
t497 = t314 / 0.2e1;
t276 = rSges(3,1) * t311 + rSges(3,2) * t313;
t495 = m(3) * t276;
t494 = pkin(2) * t313;
t301 = t312 * pkin(6);
t163 = (-Icges(5,2) * t299 - t477) * t427 + (Icges(5,6) * t311 + t313 * t362) * qJD(2);
t489 = (-t217 * t430 + (-qJD(4) * t219 - t163) * t311) * t298 + t519;
t486 = -rSges(6,2) * t406 + t518;
t485 = rSges(6,2) * t322 - t530;
t484 = rSges(6,2) * t311;
t483 = rSges(3,3) * t314;
t482 = rSges(5,3) * t311;
t481 = -rSges(4,3) - qJ(3);
t378 = -rSges(5,1) * t233 + rSges(5,2) * t232;
t157 = rSges(5,3) * t461 - t378;
t469 = t157 * t314;
t465 = t308 * t313;
t464 = t308 * t314;
t463 = t309 * t313;
t310 = -pkin(7) - qJ(3);
t462 = t310 * t313;
t456 = qJ(3) + t310;
t455 = t526 * t432;
t376 = rSges(6,1) * t299 + rSges(6,3) * t298;
t454 = (pkin(4) * t430 + qJ(5) * t427) * t299 + (qJ(5) * t430 + (-pkin(4) * qJD(4) + qJD(5)) * t311) * t298 + (-rSges(6,1) * t298 + rSges(6,3) * t299) * t427 + (t313 * t376 + t484) * qJD(2);
t182 = t254 * rSges(4,1) + t253 * rSges(4,2) + rSges(4,3) * t460;
t294 = pkin(2) * t457;
t256 = qJ(3) * t460 + t294;
t451 = -t182 - t256;
t287 = pkin(3) * t459;
t343 = t295 * t457 - t310 * t460 + t287;
t184 = t343 - t256;
t450 = -t184 - t256;
t375 = t470 + t494;
t250 = qJD(2) * t375 - qJD(3) * t313;
t449 = -(-t311 * t456 - t395) * qJD(2) - t250;
t210 = t313 * t456 - t517;
t275 = pkin(2) * t311 - qJ(3) * t313;
t257 = t275 * t435;
t448 = t210 * t435 + t257;
t447 = -t210 - t275;
t380 = rSges(4,1) * t309 - rSges(4,2) * t308;
t446 = -(rSges(4,3) * t311 + t313 * t380) * qJD(2) - t250;
t445 = -rSges(6,2) * t313 + (pkin(4) * t299 + qJ(5) * t298 + t376) * t311;
t231 = -rSges(4,3) * t313 + t311 * t380;
t444 = -t231 - t275;
t255 = t375 * t312;
t443 = t312 * t255 + t314 * t256;
t288 = pkin(3) * t464;
t442 = qJD(1) * t288 + t310 * t406;
t441 = -t310 * t461 - t288;
t440 = rSges(3,2) * t406 + rSges(3,3) * t433;
t428 = qJD(3) * t311;
t286 = t314 * t428;
t297 = pkin(6) * t433;
t439 = t286 + t297;
t438 = t314 * pkin(1) + t301;
t437 = t312 ^ 2 + t314 ^ 2;
t237 = Icges(3,3) * t312 + t314 * t360;
t436 = qJD(1) * t237;
t35 = t54 * t312 - t314 * t53;
t36 = t56 * t312 - t314 * t55;
t422 = t35 / 0.2e1 + t36 / 0.2e1;
t37 = t58 * t312 - t314 * t57;
t38 = t60 * t312 - t314 * t59;
t421 = t37 / 0.2e1 + t38 / 0.2e1;
t252 = t309 * t458 - t464;
t335 = t308 * t458 + t309 * t314;
t175 = Icges(4,5) * t252 - Icges(4,6) * t335 + Icges(4,3) * t461;
t420 = t175 * t461;
t419 = t175 * t460;
t176 = Icges(4,5) * t254 + Icges(4,6) * t253 + Icges(4,3) * t460;
t418 = t176 * t461;
t417 = t176 * t460;
t416 = t138 * rSges(5,1) + t137 * rSges(5,2) + rSges(5,3) * t396;
t377 = rSges(5,1) * t299 - rSges(5,2) * t298;
t167 = (-rSges(5,1) * t298 - rSges(5,2) * t299) * t427 + (t313 * t377 + t482) * qJD(2);
t414 = -t167 + t449;
t278 = qJ(3) * t396;
t284 = pkin(2) * t403;
t323 = -t312 * t434 - t401;
t413 = t312 * (qJ(3) * t322 + qJD(1) * t294 + t312 * t428 - t284) + t314 * (pkin(2) * t323 - qJ(3) * t406 + t278 + t286) + t255 * t433;
t201 = qJD(1) * t335 + t308 * t401;
t202 = -qJD(1) * t252 - t309 * t401;
t412 = t202 * rSges(4,1) + t201 * rSges(4,2) + rSges(4,3) * t396;
t222 = -rSges(5,3) * t313 + t311 * t377;
t411 = -t222 + t447;
t159 = t235 * rSges(5,1) - t234 * rSges(5,2) + rSges(5,3) * t460;
t410 = t295 * t403 + t322 * t310;
t302 = t314 * pkin(6);
t409 = t302 - t441;
t407 = t222 * t435;
t402 = t311 * t430;
t394 = -t295 * t313 - pkin(1);
t393 = t312 * t445;
t392 = t314 * t445;
t391 = t453 * t314;
t192 = t444 * t314;
t390 = qJD(1) * t445;
t388 = t449 - t454;
t183 = -t312 * t525 + t441;
t387 = t312 * t183 + t314 * t184 + t443;
t386 = -t445 + t447;
t126 = t411 * t314;
t385 = (-pkin(3) * t308 - pkin(6)) * t312;
t384 = t312 * t390;
t383 = rSges(3,1) * t313 - rSges(3,2) * t311;
t203 = qJD(1) * t253 + t308 * t403;
t204 = qJD(1) * t254 - t309 * t403;
t382 = -t204 * rSges(4,1) - t203 * rSges(4,2);
t381 = -rSges(4,1) * t252 + rSges(4,2) * t335;
t379 = t140 * rSges(5,1) - t139 * rSges(5,2);
t369 = Icges(3,1) * t311 + t478;
t364 = Icges(3,2) * t313 + t479;
t359 = Icges(4,5) * t309 - Icges(4,6) * t308;
t352 = -t159 * t312 + t469;
t351 = -t312 * t157 - t159 * t314;
t245 = rSges(3,1) * t457 + t516;
t108 = t386 * t314;
t345 = -t313 * t488 + t521;
t344 = -pkin(1) - t383;
t33 = t139 * t214 + t140 * t218 + t232 * t160 + t162 * t461 + t233 * t164 + t216 * t322;
t34 = -t139 * t217 + t140 * t219 + t161 * t461 - t232 * t163 + t233 * t165 + t215 * t322;
t342 = t20 / 0.2e1 + t22 / 0.2e1 + t33 / 0.2e1 + t34 / 0.2e1;
t31 = -t137 * t214 + t138 * t218 + t234 * t160 + t162 * t460 + t235 * t164 + t216 * t321;
t32 = t137 * t217 + t138 * t219 + t161 * t460 - t234 * t163 + t235 * t165 + t215 * t321;
t341 = t21 / 0.2e1 + t23 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1;
t340 = t63 / 0.2e1 + t65 / 0.2e1 + t97 / 0.2e1 + t96 / 0.2e1;
t339 = t98 / 0.2e1 + t99 / 0.2e1 + t64 / 0.2e1 + t66 / 0.2e1;
t338 = t394 - t484;
t337 = t394 - t482;
t336 = t312 * (-qJ(3) * t400 + t284 + (-t314 * t525 + t287) * qJD(1) - t410) + t314 * (-t278 + (-t462 + t517) * t429 + t525 * t435 + t442) + t183 * t433 + t413;
t332 = qJD(2) * t276;
t329 = qJD(2) * t369;
t328 = qJD(2) * t364;
t327 = qJD(2) * (-Icges(3,5) * t311 - Icges(3,6) * t313);
t326 = t311 * t481 - pkin(1) - t494;
t325 = t338 * t312;
t324 = t337 * t312;
t320 = t343 + t438;
t318 = -t312 * t452 + t391;
t317 = t326 * t312;
t315 = (-t295 * t311 - t462) * t429 + t439 + t442;
t305 = t311 ^ 2;
t262 = t383 * qJD(2);
t259 = t360 * qJD(2);
t244 = t312 * t383 - t483;
t213 = (Icges(4,5) * t311 + t313 * t368) * qJD(2);
t212 = (Icges(4,6) * t311 + t313 * t363) * qJD(2);
t209 = t245 + t438;
t208 = t312 * t344 + t302 + t483;
t191 = t444 * t312;
t186 = t312 * t327 + t436;
t185 = -qJD(1) * t236 + t314 * t327;
t181 = rSges(4,3) * t461 - t381;
t180 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t460;
t179 = Icges(4,1) * t252 - Icges(4,4) * t335 + Icges(4,5) * t461;
t178 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t460;
t177 = Icges(4,4) * t252 - Icges(4,2) * t335 + Icges(4,6) * t461;
t143 = t276 * t431 + ((-rSges(3,3) - pkin(6)) * t312 + t344 * t314) * qJD(1);
t142 = rSges(3,1) * t323 - rSges(3,2) * t396 - pkin(1) * t435 + t297 + t440;
t128 = t438 - t451;
t127 = t302 + t317 + t381;
t125 = t411 * t312;
t124 = t312 * t237 - t347 * t314;
t123 = t312 * t236 - t331;
t122 = -t237 * t314 - t330;
t121 = -t236 * t314 - t312 * t348;
t120 = Icges(4,1) * t204 + Icges(4,4) * t203 + Icges(4,5) * t322;
t119 = Icges(4,1) * t202 + Icges(4,4) * t201 + Icges(4,5) * t321;
t118 = Icges(4,4) * t204 + Icges(4,2) * t203 + Icges(4,6) * t322;
t117 = Icges(4,4) * t202 + Icges(4,2) * t201 + Icges(4,6) * t321;
t114 = -t313 * t159 - t222 * t460;
t113 = t157 * t313 + t222 * t461;
t112 = qJD(1) * t192 + t312 * t446;
t111 = t231 * t435 + t314 * t446 + t257;
t110 = t320 + t159;
t109 = t324 + t378 + t409;
t107 = t386 * t312;
t100 = t352 * t311;
t95 = t312 * t181 + t182 * t314 + t443;
t94 = t320 + t452;
t93 = t325 + t409 - t531;
t92 = t284 + (t430 * t481 - t428) * t312 + (t314 * t326 - t301) * qJD(1) + t382;
t91 = -pkin(2) * t401 + qJD(1) * t317 + t278 + t412 + t439;
t90 = rSges(5,3) * t322 + t379;
t88 = -rSges(5,3) * t406 + t416;
t72 = -t311 * t392 - t313 * t452;
t71 = t311 * t393 + t313 * t453;
t70 = t253 * t178 + t254 * t180 + t417;
t69 = t253 * t177 + t254 * t179 + t419;
t68 = -t178 * t335 + t180 * t252 + t418;
t67 = -t177 * t335 + t179 * t252 + t420;
t62 = qJD(1) * t126 + t312 * t414;
t61 = t314 * t414 + t407 + t448;
t52 = (-rSges(5,3) * t430 - t428) * t312 + (t314 * t337 + t385) * qJD(1) - t379 + t410;
t51 = qJD(1) * t324 + t315 + t416;
t50 = t318 * t311;
t49 = -t351 + t387;
t48 = qJD(1) * t108 + t312 * t388;
t47 = t314 * t388 + t384 + t448;
t46 = (t222 * t431 + t90) * t313 + (-qJD(2) * t157 + t312 * t167 + t222 * t433) * t311;
t45 = (-t222 * t429 - t88) * t313 + (qJD(2) * t159 - t167 * t314 + t407) * t311;
t44 = t387 - t514;
t41 = (-rSges(6,2) * t430 - t428) * t312 + (t314 * t338 + t385) * qJD(1) + t410 + t530;
t40 = qJD(1) * t325 + t315 + t518;
t39 = t312 * (rSges(4,3) * t400 - t382) + t314 * t412 + (t314 * t181 + t312 * t451) * qJD(1) + t413;
t30 = t352 * t430 + (qJD(1) * t351 - t312 * t88 + t314 * t90) * t311;
t25 = (qJD(2) * t393 + t485) * t313 + (-qJD(2) * t453 + t312 * t454 + t314 * t390) * t311;
t24 = (-qJD(2) * t392 - t486) * t313 + (qJD(2) * t452 - t314 * t454 + t384) * t311;
t19 = -t139 * t151 + t140 * t155 + t147 * t322 - t232 * t81 + t233 * t85 + t461 * t77;
t18 = -t139 * t150 + t140 * t154 + t146 * t322 - t232 * t82 + t233 * t86 + t461 * t78;
t17 = t139 * t145 + t140 * t153 + t149 * t322 + t232 * t75 + t233 * t83 + t461 * t79;
t16 = t139 * t144 + t140 * t152 + t148 * t322 + t232 * t76 + t233 * t84 + t461 * t80;
t15 = t137 * t151 + t138 * t155 + t147 * t321 - t234 * t81 + t235 * t85 + t460 * t77;
t14 = t137 * t150 + t138 * t154 + t146 * t321 - t234 * t82 + t235 * t86 + t460 * t78;
t13 = -t137 * t145 + t138 * t153 + t149 * t321 + t234 * t75 + t235 * t83 + t460 * t79;
t12 = -t137 * t144 + t138 * t152 + t148 * t321 + t234 * t76 + t235 * t84 + t460 * t80;
t11 = t312 * t90 + t314 * t88 + (t469 + (-t159 + t450) * t312) * qJD(1) + t336;
t10 = t318 * t430 + (qJD(1) * t514 - t486 * t312 + t485 * t314) * t311;
t9 = t486 * t314 + t485 * t312 + (t391 + (t450 - t452) * t312) * qJD(1) + t336;
t8 = qJD(1) * t373 - t18 * t314 + t19 * t312;
t7 = qJD(1) * t374 - t16 * t314 + t17 * t312;
t6 = qJD(1) * t371 - t14 * t314 + t15 * t312;
t5 = qJD(1) * t372 - t12 * t314 + t13 * t312;
t4 = (qJD(2) * t373 - t34) * t313 + (-qJD(1) * t36 + qJD(2) * t97 + t18 * t312 + t19 * t314) * t311;
t3 = (qJD(2) * t374 - t33) * t313 + (-qJD(1) * t35 + qJD(2) * t96 + t16 * t312 + t17 * t314) * t311;
t2 = (qJD(2) * t371 - t32) * t313 + (-qJD(1) * t38 + qJD(2) * t99 + t14 * t312 + t15 * t314) * t311;
t1 = (qJD(2) * t372 - t31) * t313 + (-qJD(1) * t37 + qJD(2) * t98 + t12 * t312 + t13 * t314) * t311;
t26 = [-t217 * t405 - t163 * t468 - t219 * t399 + (t142 * t209 + t143 * t208) * t513 + (t127 * t92 + t128 * t91) * t512 + (t109 * t52 + t110 * t51) * t511 + (t40 * t94 + t41 * t93) * t510 + (-t212 * t308 + t213 * t309) * t311 + (-Icges(4,3) * t313 + t311 * t359 - t364 + t370) * t432 + (-Icges(4,3) * t311 - t229 * t308 + t230 * t309 - t313 * t359 + t365 + t369) * t430 + t519; (-t203 * t229 / 0.2e1 - t204 * t230 / 0.2e1 + t335 * t212 / 0.2e1 - t252 * t213 / 0.2e1 + t259 * t497 - t342) * t314 + (t201 * t506 + t202 * t505 + t212 * t504 + t213 * t503 + t259 * t500 + t341) * t312 + m(3) * ((-t142 * t312 - t143 * t314) * t276 + (-t208 * t314 - t209 * t312) * t262) + m(4) * (t111 * t127 + t112 * t128 + t191 * t91 + t192 * t92) + m(5) * (t109 * t61 + t110 * t62 + t125 * t51 + t126 * t52) + m(6) * (t107 * t40 + t108 * t41 + t47 * t93 + t48 * t94) + ((Icges(4,5) * t204 / 0.2e1 + Icges(4,6) * t203 / 0.2e1 + Icges(4,3) * t322 / 0.2e1 + t239 * t528 + t328 * t500) * t314 + (-Icges(4,5) * t202 / 0.2e1 - Icges(4,6) * t201 / 0.2e1 - Icges(4,3) * t321 / 0.2e1 + t238 * t528 + t328 * t498) * t312) * t313 + ((-qJD(1) * t240 - t117 * t308 + t119 * t309 - t314 * t329) * t500 + (qJD(1) * t241 - t118 * t308 + t120 * t309 - t312 * t329) * t498) * t311 + ((t175 * t311 - t177 * t465 + t179 * t463) * t498 + t331 / 0.2e1 + (t176 * t311 - t178 * t465 + t180 * t463) * t500 - t330 / 0.2e1) * qJD(2) + ((-t209 * t495 + t229 * t504 + t230 * t503 + (-t176 / 0.2e1 + t239 / 0.2e1) * t313 + (t178 * t502 + t180 * t501 + t241 / 0.2e1) * t311 + t339) * t314 + (t208 * t495 - t335 * t506 + t252 * t505 + (-t175 / 0.2e1 + t238 / 0.2e1) * t313 + (t177 * t502 + t179 * t501 + t240 / 0.2e1) * t311 + t340) * t312) * qJD(1); (t107 * t48 + t108 * t47 + t44 * t9) * t510 + (t11 * t49 + t125 * t62 + t126 * t61) * t511 - t314 * t7 + t312 * t6 - t314 * t8 + t312 * t5 + (t111 * t192 + t112 * t191 + t39 * t95) * t512 + ((t312 * t244 + t245 * t314) * ((qJD(1) * t244 - t314 * t332 + t440) * t314 + (-t312 * t332 + (-t245 + t516) * qJD(1)) * t312) + t437 * t276 * t262) * t513 - t314 * ((t335 * t118 - t252 * t120 - t203 * t177 - t204 * t179 + (t68 - t419) * qJD(1)) * t314 + (-t335 * t117 + t252 * t119 + t203 * t178 + t204 * t180 + (t67 + t417) * qJD(1)) * t312) + t312 * ((t312 * t185 + (t123 + t330) * qJD(1)) * t312 + (t124 * qJD(1) + (t238 * t430 + t240 * t432) * t314 + (-t186 + (-t239 * t313 - t241 * t311) * qJD(2) + (t237 - t348) * qJD(1)) * t312) * t314) - t314 * ((t314 * t186 + (t122 + t331) * qJD(1)) * t314 + (t121 * qJD(1) + (-t239 * t430 - t241 * t432 + t436) * t312 + (-t185 + (t238 * t313 + t240 * t311) * qJD(2) - t347 * qJD(1)) * t314) * t312) + t312 * ((t253 * t117 + t254 * t119 + t201 * t178 + t202 * t180 + (t69 - t418) * qJD(1)) * t312 + (-t253 * t118 - t254 * t120 - t201 * t177 - t202 * t179 + (t70 + t420) * qJD(1)) * t314) + (t35 + t36 + (-t121 - t67) * t314 + (t122 + t68) * t312) * t435 + (t37 + t38 + (-t123 - t69) * t314 + (t124 + t70) * t312) * t433; 0.2e1 * ((t312 * t94 + t314 * t93) * t507 + (t109 * t314 + t110 * t312) * t508 + (t127 * t314 + t128 * t312) * t509) * t430 + 0.2e1 * ((t312 * t40 + t314 * t41 + t433 * t94 - t435 * t93) * t507 + (-t109 * t435 + t110 * t433 + t312 * t51 + t314 * t52) * t508 + (-t127 * t435 + t128 * t433 + t312 * t91 + t314 * t92) * t509) * t311; 0.2e1 * ((t107 * t431 + t108 * t429 - t9) * t507 + (t125 * t431 + t126 * t429 - t11) * t508 + (t191 * t431 + t192 * t429 - t39) * t509) * t313 + 0.2e1 * ((qJD(2) * t44 + t107 * t433 - t108 * t435 + t312 * t48 + t314 * t47) * t507 + (qJD(2) * t49 + t125 * t433 - t126 * t435 + t312 * t62 + t314 * t61) * t508 + (qJD(2) * t95 + t111 * t314 + t112 * t312 + t191 * t433 - t192 * t435) * t509) * t311; 0.4e1 * (t509 + t508 + t507) * (-0.1e1 + t437) * t402; m(5) * (t109 * t46 + t110 * t45 + t113 * t52 + t114 * t51) + m(6) * (t24 * t94 + t25 * t93 + t40 * t72 + t41 * t71) + ((t312 * t340 + t314 * t339) * qJD(2) - t489) * t313 + (t341 * t314 + t342 * t312 + (-t312 * t339 + t314 * t340) * qJD(1)) * t311 - t455; m(5) * (t100 * t11 + t113 * t61 + t114 * t62 + t125 * t45 + t126 * t46 + t30 * t49) + m(6) * (t10 * t44 + t107 * t24 + t108 * t25 + t47 * t71 + t48 * t72 + t50 * t9) + (-t3 / 0.2e1 - t4 / 0.2e1 + t421 * t430) * t314 + (t1 / 0.2e1 + t2 / 0.2e1 + t422 * t430) * t312 + ((-t312 * t421 + t314 * t422) * qJD(1) + (t7 + t8) * t500 + (t5 + t6) * t497 + (t312 * t487 - t314 * t488) * qJD(2) / 0.2e1) * t311 - (t524 * qJD(1) + t522 * t312 + t523 * t314) * t313 / 0.2e1 + (t521 * t312 + t520 * t314) * qJD(1) / 0.2e1; 0.2e1 * ((t113 * t429 + t114 * t431 - t30) * t508 + (t429 * t71 + t431 * t72 - t10) * t507) * t313 + 0.2e1 * ((qJD(2) * t100 - t113 * t435 + t114 * t433 + t312 * t45 + t314 * t46) * t508 + (qJD(2) * t50 + t24 * t312 + t25 * t314 + t433 * t72 - t435 * t71) * t507) * t311; (t10 * t50 + t24 * t72 + t25 * t71) * t510 + (t100 * t30 + t113 * t46 + t114 * t45) * t511 + (t489 * t313 + (t345 * t312 - t314 * t515) * qJD(2) + t455) * t313 + ((-t313 * t522 + t1 + t2) * t314 + (t313 * t523 + t3 + t4) * t312 + (t524 * t311 + t526 * t313) * qJD(2) + (t312 * t515 + t345 * t314) * qJD(1)) * t311; m(6) * (-t137 * t93 + t139 * t94 + t232 * t40 + t234 * t41); m(6) * (t44 * t398 + t107 * t139 - t108 * t137 + t232 * t48 + t234 * t47 + (t311 * t9 + t430 * t44) * t298); m(6) * ((t298 * t305 + (t232 * t312 + t234 * t314 - t298 * t313) * t313) * qJD(2) + (-t397 - t137 * t314 + t139 * t312 + (t232 * t314 - t234 * t312) * qJD(1)) * t311); m(6) * (t50 * t398 - t137 * t71 + t139 * t72 + t232 * t24 + t234 * t25 + (t10 * t311 + t430 * t50) * t298); (-t137 * t234 + t139 * t232 + (qJD(4) * t299 * t305 + t298 * t402) * t298) * t510;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t26(1), t26(2), t26(4), t26(7), t26(11); t26(2), t26(3), t26(5), t26(8), t26(12); t26(4), t26(5), t26(6), t26(9), t26(13); t26(7), t26(8), t26(9), t26(10), t26(14); t26(11), t26(12), t26(13), t26(14), t26(15);];
Mq = res;
