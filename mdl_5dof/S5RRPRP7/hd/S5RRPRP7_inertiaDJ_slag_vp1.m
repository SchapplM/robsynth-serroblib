% Calculate time derivative of joint inertia matrix for
% S5RRPRP7
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP7_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP7_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:56
% EndTime: 2019-12-31 20:00:23
% DurationCPUTime: 14.87s
% Computational Cost: add. (18545->807), mult. (28180->1132), div. (0->0), fcn. (26932->8), ass. (0->392)
t281 = qJ(2) + pkin(8);
t274 = sin(t281);
t275 = cos(t281);
t283 = sin(qJ(4));
t286 = cos(qJ(4));
t332 = Icges(5,5) * t286 - Icges(5,6) * t283;
t188 = -Icges(5,3) * t275 + t274 * t332;
t335 = Icges(6,4) * t286 + Icges(6,6) * t283;
t189 = -Icges(6,2) * t275 + t274 * t335;
t509 = t188 + t189;
t288 = cos(qJ(1));
t421 = t286 * t288;
t285 = sin(qJ(1));
t423 = t285 * t283;
t227 = t275 * t423 + t421;
t422 = t285 * t286;
t424 = t283 * t288;
t228 = t275 * t422 - t424;
t504 = rSges(6,3) + qJ(5);
t506 = rSges(6,1) + pkin(4);
t508 = -t504 * t227 - t506 * t228;
t403 = qJD(2) * t285;
t380 = t274 * t403;
t397 = qJD(4) * t286;
t398 = qJD(4) * t283;
t408 = qJD(1) * t288;
t409 = qJD(1) * t285;
t134 = -t283 * t380 - t288 * t398 - t286 * t409 + (t283 * t408 + t285 * t397) * t275;
t315 = (-qJD(4) * t275 + qJD(1)) * t283;
t364 = qJD(1) * t275 - qJD(4);
t402 = qJD(2) * t286;
t135 = t364 * t421 + (-t274 * t402 + t315) * t285;
t507 = t227 * qJD(5) + t504 * t134 + t506 * t135;
t444 = Icges(6,5) * t286;
t331 = Icges(6,3) * t283 + t444;
t187 = -Icges(6,6) * t275 + t274 * t331;
t448 = Icges(5,4) * t286;
t336 = -Icges(5,2) * t283 + t448;
t190 = -Icges(5,6) * t275 + t274 * t336;
t445 = Icges(6,5) * t283;
t341 = Icges(6,1) * t286 + t445;
t191 = -Icges(6,4) * t275 + t274 * t341;
t449 = Icges(5,4) * t283;
t342 = Icges(5,1) * t286 - t449;
t192 = -Icges(5,5) * t275 + t274 * t342;
t505 = t509 * t275 + ((-t191 - t192) * t286 + (-t187 + t190) * t283) * t274;
t405 = qJD(2) * t283;
t378 = t275 * t405;
t490 = -t274 * t397 - t378;
t229 = t275 * t424 - t422;
t389 = t275 * t421;
t230 = t389 + t423;
t427 = t274 * t288;
t151 = Icges(6,4) * t230 + Icges(6,2) * t427 + Icges(6,6) * t229;
t147 = Icges(6,5) * t230 + Icges(6,6) * t427 + Icges(6,3) * t229;
t155 = Icges(6,1) * t230 + Icges(6,4) * t427 + Icges(6,5) * t229;
t327 = t147 * t283 + t155 * t286;
t61 = -t151 * t275 + t274 * t327;
t149 = Icges(5,5) * t230 - Icges(5,6) * t229 + Icges(5,3) * t427;
t153 = Icges(5,4) * t230 - Icges(5,2) * t229 + Icges(5,6) * t427;
t157 = Icges(5,1) * t230 - Icges(5,4) * t229 + Icges(5,5) * t427;
t325 = -t153 * t283 + t157 * t286;
t63 = -t149 * t275 + t274 * t325;
t462 = t61 + t63;
t429 = t274 * t285;
t150 = Icges(6,4) * t228 + Icges(6,2) * t429 + Icges(6,6) * t227;
t146 = Icges(6,5) * t228 + Icges(6,6) * t429 + Icges(6,3) * t227;
t154 = Icges(6,1) * t228 + Icges(6,4) * t429 + Icges(6,5) * t227;
t328 = t146 * t283 + t154 * t286;
t60 = -t150 * t275 + t274 * t328;
t148 = Icges(5,5) * t228 - Icges(5,6) * t227 + Icges(5,3) * t429;
t152 = Icges(5,4) * t228 - Icges(5,2) * t227 + Icges(5,6) * t429;
t156 = Icges(5,1) * t228 - Icges(5,4) * t227 + Icges(5,5) * t429;
t326 = -t152 * t283 + t156 * t286;
t62 = -t148 * t275 + t274 * t326;
t463 = t60 + t62;
t503 = t285 * t463 + t288 * t462;
t502 = t285 * t462 - t288 * t463;
t373 = t275 * t403;
t295 = t274 * t408 + t373;
t69 = Icges(6,5) * t135 + Icges(6,6) * t295 + Icges(6,3) * t134;
t73 = Icges(6,4) * t135 + Icges(6,2) * t295 + Icges(6,6) * t134;
t77 = Icges(6,1) * t135 + Icges(6,4) * t295 + Icges(6,5) * t134;
t19 = (qJD(2) * t328 - t73) * t275 + (qJD(2) * t150 + t283 * t69 + t286 * t77 + (t146 * t286 - t154 * t283) * qJD(4)) * t274;
t71 = Icges(5,5) * t135 - Icges(5,6) * t134 + Icges(5,3) * t295;
t75 = Icges(5,4) * t135 - Icges(5,2) * t134 + Icges(5,6) * t295;
t79 = Icges(5,1) * t135 - Icges(5,4) * t134 + Icges(5,5) * t295;
t21 = (qJD(2) * t326 - t71) * t275 + (qJD(2) * t148 - t283 * t75 + t286 * t79 + (-t152 * t286 - t156 * t283) * qJD(4)) * t274;
t501 = -t19 - t21;
t400 = qJD(2) * t288;
t379 = t274 * t400;
t132 = qJD(1) * t227 - qJD(4) * t389 + t283 * t379 - t285 * t398;
t133 = t288 * t315 + (-t285 * t364 - t379) * t286;
t376 = t275 * t400;
t383 = t274 * t409;
t294 = t376 - t383;
t68 = Icges(6,5) * t133 + Icges(6,6) * t294 - Icges(6,3) * t132;
t72 = Icges(6,4) * t133 + Icges(6,2) * t294 - Icges(6,6) * t132;
t76 = Icges(6,1) * t133 + Icges(6,4) * t294 - Icges(6,5) * t132;
t20 = (qJD(2) * t327 - t72) * t275 + (qJD(2) * t151 + t283 * t68 + t286 * t76 + (t147 * t286 - t155 * t283) * qJD(4)) * t274;
t70 = Icges(5,5) * t133 + Icges(5,6) * t132 + Icges(5,3) * t294;
t74 = Icges(5,4) * t133 + Icges(5,2) * t132 + Icges(5,6) * t294;
t78 = Icges(5,1) * t133 + Icges(5,4) * t132 + Icges(5,5) * t294;
t22 = (qJD(2) * t325 - t70) * t275 + (qJD(2) * t149 - t283 * t74 + t286 * t78 + (-t153 * t286 - t157 * t283) * qJD(4)) * t274;
t500 = t20 + t22;
t54 = t148 * t429 - t152 * t227 + t156 * t228;
t55 = t149 * t429 - t153 * t227 + t157 * t228;
t349 = t285 * t54 + t288 * t55;
t52 = t146 * t227 + t150 * t429 + t154 * t228;
t53 = t147 * t227 + t151 * t429 + t155 * t228;
t350 = t285 * t52 + t288 * t53;
t86 = t187 * t227 + t189 * t429 + t191 * t228;
t87 = t188 * t429 - t190 * t227 + t192 * t228;
t466 = (-t86 - t87) * t275 + (t349 + t350) * t274;
t58 = t148 * t427 - t229 * t152 + t230 * t156;
t59 = t149 * t427 - t229 * t153 + t230 * t157;
t347 = t285 * t58 + t288 * t59;
t56 = t229 * t146 + t150 * t427 + t230 * t154;
t57 = t229 * t147 + t151 * t427 + t230 * t155;
t348 = t285 * t56 + t288 * t57;
t88 = t229 * t187 + t189 * t427 + t230 * t191;
t89 = t188 * t427 - t229 * t190 + t230 * t192;
t499 = (-t88 - t89) * t275 + (t347 + t348) * t274;
t498 = rSges(6,2) * t376 + t229 * qJD(5) - t504 * t132 + t506 * t133;
t399 = qJD(4) * t274;
t126 = (Icges(6,3) * t286 - t445) * t399 + (Icges(6,6) * t274 + t275 * t331) * qJD(2);
t128 = (-Icges(6,4) * t283 + Icges(6,6) * t286) * t399 + (Icges(6,2) * t274 + t275 * t335) * qJD(2);
t130 = (-Icges(6,1) * t283 + t444) * t399 + (Icges(6,4) * t274 + t275 * t341) * qJD(2);
t131 = (-Icges(5,1) * t283 - t448) * t399 + (Icges(5,5) * t274 + t275 * t342) * qJD(2);
t375 = t274 * t398;
t377 = t275 * t402;
t407 = qJD(2) * t274;
t430 = t274 * t283;
t497 = t192 * t377 + t126 * t430 - t275 * t128 + (-t375 + t377) * t191 - t490 * t187 + (t131 + t130) * t274 * t286 + t509 * t407;
t284 = sin(qJ(2));
t287 = cos(qJ(2));
t452 = Icges(3,4) * t287;
t340 = -Icges(3,2) * t284 + t452;
t213 = Icges(3,6) * t285 + t288 * t340;
t453 = Icges(3,4) * t284;
t346 = Icges(3,1) * t287 - t453;
t215 = Icges(3,5) * t285 + t288 * t346;
t317 = t213 * t284 - t215 * t287;
t496 = t285 * t317;
t450 = Icges(4,4) * t275;
t338 = -Icges(4,2) * t274 + t450;
t202 = Icges(4,6) * t285 + t288 * t338;
t451 = Icges(4,4) * t274;
t344 = Icges(4,1) * t275 - t451;
t204 = Icges(4,5) * t285 + t288 * t344;
t319 = t202 * t274 - t204 * t275;
t495 = t285 * t319;
t271 = pkin(2) * t287 + pkin(1);
t468 = pkin(1) - t271;
t494 = t285 * t468;
t212 = -Icges(3,6) * t288 + t285 * t340;
t214 = -Icges(3,5) * t288 + t285 * t346;
t318 = t212 * t284 - t214 * t287;
t493 = t288 * t318;
t201 = -Icges(4,6) * t288 + t285 * t338;
t203 = -Icges(4,5) * t288 + t285 * t344;
t320 = t201 * t274 - t203 * t275;
t492 = t288 * t320;
t278 = t285 * rSges(4,3);
t491 = -rSges(4,2) * t427 + t278;
t333 = Icges(4,5) * t275 - Icges(4,6) * t274;
t199 = -Icges(4,3) * t288 + t285 * t333;
t334 = Icges(3,5) * t287 - Icges(3,6) * t284;
t210 = -Icges(3,3) * t288 + t285 * t334;
t418 = rSges(6,2) * t427 + t504 * t229 + t506 * t230;
t419 = rSges(6,2) * t429 - t508;
t489 = -t285 * t419 - t288 * t418;
t488 = 2 * m(3);
t487 = 2 * m(4);
t486 = 2 * m(5);
t485 = 2 * m(6);
t484 = t285 ^ 2;
t483 = t288 ^ 2;
t481 = t285 / 0.2e1;
t477 = -rSges(6,2) - pkin(7);
t476 = -rSges(5,3) - pkin(7);
t255 = rSges(3,1) * t284 + rSges(3,2) * t287;
t475 = m(3) * t255;
t474 = pkin(2) * t284;
t473 = pkin(3) * t274;
t472 = pkin(3) * t275;
t471 = t285 * pkin(6);
t280 = t288 * pkin(6);
t282 = -qJ(3) - pkin(6);
t467 = -pkin(6) - t282;
t127 = (-Icges(5,5) * t283 - Icges(5,6) * t286) * t399 + (Icges(5,3) * t274 + t275 * t332) * qJD(2);
t129 = (-Icges(5,2) * t286 - t449) * t399 + (Icges(5,6) * t274 + t275 * t336) * qJD(2);
t464 = (-t190 * t405 - t127) * t275 + (-t129 * t283 + (-t190 * t286 - t192 * t283) * qJD(4)) * t274 + t497;
t461 = -rSges(6,2) * t383 + t498;
t460 = rSges(6,2) * t295 + t507;
t459 = t505 * t407;
t458 = rSges(3,1) * t287;
t457 = rSges(4,1) * t275;
t456 = rSges(3,2) * t284;
t455 = rSges(3,3) * t288;
t279 = t285 * rSges(3,3);
t354 = -t228 * rSges(5,1) + t227 * rSges(5,2);
t159 = rSges(5,3) * t429 - t354;
t439 = t159 * t288;
t438 = t201 * t275;
t437 = t202 * t275;
t436 = t203 * t274;
t435 = t204 * t274;
t434 = t212 * t287;
t433 = t213 * t287;
t432 = t214 * t284;
t431 = t215 * t284;
t426 = t275 * t288;
t425 = t282 * t288;
t351 = pkin(4) * t286 + qJ(5) * t283;
t352 = rSges(6,1) * t286 + rSges(6,3) * t283;
t406 = qJD(2) * t275;
t420 = t351 * t406 + (qJD(5) * t283 + (-pkin(4) * t283 + qJ(5) * t286) * qJD(4)) * t274 + (-rSges(6,1) * t283 + rSges(6,3) * t286) * t399 + (rSges(6,2) * t274 + t275 * t352) * qJD(2);
t197 = t280 + t425 - t494;
t261 = t288 * t271;
t198 = -pkin(1) * t288 + t285 * t467 + t261;
t417 = t285 * t197 + t288 * t198;
t416 = -rSges(6,2) * t275 + (t351 + t352) * t274;
t260 = pkin(3) * t426;
t225 = pkin(7) * t427 + t260;
t415 = -t198 - t225;
t240 = -pkin(7) * t275 + t473;
t382 = t284 * t409;
t264 = pkin(2) * t382;
t414 = t240 * t409 + t264;
t413 = rSges(4,2) * t383 + rSges(4,3) * t408;
t412 = t288 * t458 + t279;
t200 = Icges(4,3) * t285 + t288 * t333;
t411 = qJD(1) * t200;
t211 = Icges(3,3) * t285 + t288 * t334;
t410 = qJD(1) * t211;
t404 = qJD(2) * t284;
t401 = qJD(2) * t287;
t395 = t288 * t456;
t393 = pkin(2) * t404;
t392 = pkin(2) * t401;
t35 = t53 * t285 - t288 * t52;
t36 = t55 * t285 - t288 * t54;
t391 = t35 / 0.2e1 + t36 / 0.2e1;
t37 = t57 * t285 - t288 * t56;
t38 = t59 * t285 - t288 * t58;
t390 = t37 / 0.2e1 + t38 / 0.2e1;
t387 = t133 * rSges(5,1) + t132 * rSges(5,2) + rSges(5,3) * t376;
t276 = qJD(3) * t285;
t384 = qJD(3) * t288 + t282 * t409 + t285 * t393;
t385 = t285 * ((-t288 * t468 - t471) * qJD(1) - t384) + t288 * (-t288 * t393 + t276 + (t288 * t467 + t494) * qJD(1)) + t197 * t408;
t161 = t230 * rSges(5,1) - t229 * rSges(5,2) + rSges(5,3) * t427;
t353 = rSges(5,1) * t286 - rSges(5,2) * t283;
t194 = -rSges(5,3) * t275 + t274 * t353;
t381 = t194 * t409;
t239 = rSges(4,1) * t274 + rSges(4,2) * t275;
t372 = -t239 - t474;
t371 = -t240 - t474;
t370 = -t271 - t472;
t369 = t285 * t416;
t368 = t288 * t416;
t367 = t419 * t288;
t366 = -t285 * t282 + t261;
t365 = qJD(1) * t416;
t359 = pkin(7) * t274 + t472;
t224 = t359 * t285;
t363 = t285 * t224 + t288 * t225 + t417;
t250 = pkin(3) * t380;
t362 = t250 + t384;
t361 = -t194 + t371;
t360 = -t359 * qJD(2) - t392;
t358 = t285 * t365;
t357 = -t456 + t458;
t356 = -rSges(4,2) * t274 + t457;
t355 = t135 * rSges(5,1) - t134 * rSges(5,2);
t345 = Icges(3,1) * t284 + t452;
t343 = Icges(4,1) * t274 + t450;
t339 = Icges(3,2) * t287 + t453;
t337 = Icges(4,2) * t275 + t451;
t324 = -t161 * t285 + t439;
t323 = -t285 * t159 - t161 * t288;
t316 = t371 - t416;
t206 = rSges(4,1) * t426 + t491;
t137 = (-rSges(5,1) * t283 - rSges(5,2) * t286) * t399 + (rSges(5,3) * t274 + t275 * t353) * qJD(2);
t314 = -t137 + t360;
t313 = -pkin(1) - t357;
t312 = t366 + t225;
t311 = t86 / 0.2e1 + t87 / 0.2e1 + t62 / 0.2e1 + t60 / 0.2e1;
t310 = t89 / 0.2e1 + t63 / 0.2e1 + t61 / 0.2e1 + t88 / 0.2e1;
t123 = t361 * t288;
t309 = -t271 - t356;
t251 = pkin(7) * t376;
t307 = t285 * (pkin(7) * t295 + qJD(1) * t260 - t250) + t288 * (-pkin(7) * t383 + t251 + (-t275 * t409 - t379) * pkin(3)) + t224 * t408 + t385;
t306 = t360 - t420;
t305 = qJD(2) * t255;
t304 = qJD(2) * t239;
t303 = qJD(2) * t345;
t302 = qJD(2) * t343;
t301 = qJD(2) * t339;
t300 = qJD(2) * t337;
t299 = qJD(2) * (-Icges(3,5) * t284 - Icges(3,6) * t287);
t298 = qJD(2) * (-Icges(4,5) * t274 - Icges(4,6) * t275);
t108 = t316 * t288;
t297 = t274 * t477 + t370;
t296 = t274 * t476 + t370;
t293 = -t285 * t418 + t367;
t292 = t251 + t276 + (-t473 - t474) * t400;
t291 = rSges(3,2) * t382 + rSges(3,3) * t408 - t288 * t305;
t290 = t285 * t297 - t425;
t289 = t285 * t296 - t425;
t244 = t357 * qJD(2);
t234 = t356 * qJD(2);
t222 = -t395 + t412;
t221 = t285 * t357 - t455;
t205 = -rSges(4,3) * t288 + t285 * t356;
t196 = t372 * t288;
t195 = t372 * t285;
t186 = t471 + (pkin(1) - t456) * t288 + t412;
t185 = t285 * t313 + t280 + t455;
t175 = t206 + t366;
t174 = (rSges(4,3) - t282) * t288 + t309 * t285;
t167 = t285 * t299 + t410;
t166 = -qJD(1) * t210 + t288 * t299;
t141 = t285 * t298 + t411;
t140 = -qJD(1) * t199 + t288 * t298;
t139 = t255 * t403 + ((-rSges(3,3) - pkin(6)) * t285 + t313 * t288) * qJD(1);
t138 = (t280 + (-pkin(1) - t458) * t285) * qJD(1) + t291;
t122 = t361 * t285;
t117 = -t239 * t408 - t285 * t234 + (-t284 * t408 - t285 * t401) * pkin(2);
t116 = t239 * t409 + t264 + (-t234 - t392) * t288;
t112 = t285 * t211 - t317 * t288;
t111 = t285 * t210 - t493;
t110 = -t211 * t288 - t496;
t109 = -t210 * t288 - t285 * t318;
t107 = t316 * t285;
t106 = t239 * t403 + (t288 * t309 - t278) * qJD(1) + t384;
t105 = t276 + (-t271 - t457) * t409 + (-qJD(1) * t282 + qJD(2) * t372) * t288 + t413;
t104 = t312 + t161;
t103 = t289 + t354;
t102 = -t275 * t161 - t194 * t427;
t101 = t159 * t275 + t194 * t429;
t100 = t285 * t200 - t319 * t288;
t99 = t285 * t199 - t492;
t98 = -t200 * t288 - t495;
t97 = -t199 * t288 - t285 * t320;
t92 = t324 * t274;
t91 = t312 + t418;
t90 = t290 + t508;
t83 = rSges(5,3) * t295 + t355;
t81 = -rSges(5,3) * t383 + t387;
t67 = qJD(1) * t123 + t285 * t314;
t66 = t288 * t314 + t381 + t414;
t65 = -t274 * t368 - t275 * t418;
t64 = t274 * t369 + t275 * t419;
t51 = t296 * t408 + t373 * t476 - t355 + t362;
t50 = qJD(1) * t289 + t292 + t387;
t49 = -t323 + t363;
t48 = t293 * t274;
t47 = qJD(1) * t108 + t285 * t306;
t46 = t288 * t306 + t358 + t414;
t45 = t363 - t489;
t44 = (t194 * t403 + t83) * t275 + (-qJD(2) * t159 + t285 * t137 + t194 * t408) * t274;
t43 = (-t194 * t400 - t81) * t275 + (qJD(2) * t161 - t137 * t288 + t381) * t274;
t42 = t297 * t408 + t373 * t477 + t362 - t507;
t41 = qJD(1) * t290 + t292 + t498;
t34 = t127 * t429 - t227 * t129 + t228 * t131 - t134 * t190 + t135 * t192 + t188 * t295;
t33 = t227 * t126 + t128 * t429 + t228 * t130 + t134 * t187 + t135 * t191 + t189 * t295;
t32 = t127 * t427 - t229 * t129 + t230 * t131 + t132 * t190 + t133 * t192 + t188 * t294;
t31 = t229 * t126 + t128 * t427 + t230 * t130 - t132 * t187 + t133 * t191 + t189 * t294;
t30 = t324 * t406 + (qJD(1) * t323 - t285 * t81 + t288 * t83) * t274;
t25 = (qJD(2) * t369 + t460) * t275 + (-qJD(2) * t419 + t285 * t420 + t288 * t365) * t274;
t24 = (-qJD(2) * t368 - t461) * t275 + (qJD(2) * t418 - t288 * t420 + t358) * t274;
t23 = t285 * t83 + t288 * t81 + (t439 + (-t161 + t415) * t285) * qJD(1) + t307;
t18 = -t134 * t153 + t135 * t157 + t149 * t295 - t227 * t74 + t228 * t78 + t429 * t70;
t17 = -t134 * t152 + t135 * t156 + t148 * t295 - t227 * t75 + t228 * t79 + t429 * t71;
t16 = t134 * t147 + t135 * t155 + t151 * t295 + t227 * t68 + t228 * t76 + t429 * t72;
t15 = t134 * t146 + t135 * t154 + t150 * t295 + t227 * t69 + t228 * t77 + t429 * t73;
t14 = t132 * t153 + t133 * t157 + t149 * t294 - t229 * t74 + t230 * t78 + t427 * t70;
t13 = t132 * t152 + t133 * t156 + t148 * t294 - t229 * t75 + t230 * t79 + t427 * t71;
t12 = -t132 * t147 + t133 * t155 + t151 * t294 + t229 * t68 + t230 * t76 + t427 * t72;
t11 = -t132 * t146 + t133 * t154 + t150 * t294 + t229 * t69 + t230 * t77 + t427 * t73;
t10 = t293 * t406 + (qJD(1) * t489 - t461 * t285 + t460 * t288) * t274;
t9 = t461 * t288 + t460 * t285 + (t367 + (t415 - t418) * t285) * qJD(1) + t307;
t8 = qJD(1) * t349 - t17 * t288 + t18 * t285;
t7 = qJD(1) * t350 - t15 * t288 + t16 * t285;
t6 = qJD(1) * t347 - t13 * t288 + t14 * t285;
t5 = qJD(1) * t348 - t11 * t288 + t12 * t285;
t4 = (qJD(2) * t349 - t34) * t275 + (-qJD(1) * t36 + qJD(2) * t87 + t17 * t285 + t18 * t288) * t274;
t3 = (qJD(2) * t350 - t33) * t275 + (-qJD(1) * t35 + qJD(2) * t86 + t15 * t285 + t16 * t288) * t274;
t2 = (qJD(2) * t347 - t32) * t275 + (-qJD(1) * t38 + qJD(2) * t89 + t13 * t285 + t14 * t288) * t274;
t1 = (qJD(2) * t348 - t31) * t275 + (-qJD(1) * t37 + qJD(2) * t88 + t11 * t285 + t12 * t288) * t274;
t26 = [-t192 * t375 - t275 * t127 + (t138 * t186 + t139 * t185) * t488 + (t105 * t175 + t106 * t174) * t487 + (t41 * t91 + t42 * t90) * t485 + (t103 * t51 + t104 * t50) * t486 - t129 * t430 + (t344 - t337) * t407 + (t343 + t338) * t406 + (t346 - t339) * t404 + (t345 + t340) * t401 + t490 * t190 + t497; m(3) * ((-t138 * t285 - t139 * t288) * t255 + (-t185 * t288 - t186 * t285) * t244) + m(4) * (t105 * t195 + t106 * t196 + t116 * t174 + t117 * t175) + m(6) * (t107 * t41 + t108 * t42 + t46 * t90 + t47 * t91) + m(5) * (t103 * t66 + t104 * t67 + t122 * t50 + t123 * t51) + ((t433 / 0.2e1 + t431 / 0.2e1 - t186 * t475 + t437 / 0.2e1 + t435 / 0.2e1 + t310) * t288 + (t185 * t475 + t434 / 0.2e1 + t432 / 0.2e1 + t438 / 0.2e1 + t436 / 0.2e1 + t311) * t285) * qJD(1) + (t334 + t333) * qJD(2) * (t484 / 0.2e1 + t483 / 0.2e1) + ((-qJD(1) * t201 - t288 * t300) * t275 + (-qJD(1) * t203 - t288 * t302) * t274 + (-qJD(1) * t212 - t288 * t301) * t287 + (-qJD(1) * t214 - t288 * t303) * t284 + t31 + t32 + (-t317 - t319) * qJD(2) + t500) * t481 - ((qJD(1) * t202 - t285 * t300) * t275 + (qJD(1) * t204 - t285 * t302) * t274 + (qJD(1) * t213 - t285 * t301) * t287 + (qJD(1) * t215 - t285 * t303) * t284 + t33 + t34 + (-t318 - t320) * qJD(2) - t501) * t288 / 0.2e1; -t288 * t7 - t288 * t8 + t285 * t5 + t285 * t6 + (t107 * t47 + t108 * t46 + t45 * t9) * t485 + (t122 * t67 + t123 * t66 + t23 * t49) * t486 + (t196 * t116 + t195 * t117 + (t285 * t205 + t206 * t288 + t417) * ((qJD(1) * t205 - t288 * t304 + t413) * t288 + (-t285 * t304 + (-t198 - t206 + t491) * qJD(1)) * t285 + t385)) * t487 + t285 * ((t285 * t166 + (t111 + t496) * qJD(1)) * t285 + (t112 * qJD(1) + (t212 * t401 + t214 * t404) * t288 + (-t167 + (-t431 - t433) * qJD(2) + (t211 - t318) * qJD(1)) * t285) * t288) - t288 * ((t288 * t167 + (t110 + t493) * qJD(1)) * t288 + (t109 * qJD(1) + (-t213 * t401 - t215 * t404 + t410) * t285 + (-t166 + (t432 + t434) * qJD(2) - t317 * qJD(1)) * t288) * t285) + t285 * ((t285 * t140 + (t99 + t495) * qJD(1)) * t285 + (t100 * qJD(1) + (t201 * t406 + t203 * t407) * t288 + (-t141 + (-t435 - t437) * qJD(2) + (t200 - t320) * qJD(1)) * t285) * t288) - t288 * ((t288 * t141 + (t98 + t492) * qJD(1)) * t288 + (t97 * qJD(1) + (-t202 * t406 - t204 * t407 + t411) * t285 + (-t140 + (t436 + t438) * qJD(2) - t319 * qJD(1)) * t288) * t285) + ((t285 * t221 + t222 * t288) * ((qJD(1) * t221 + t291) * t288 + (-t285 * t305 + (-t222 - t395 + t279) * qJD(1)) * t285) + (t483 + t484) * t255 * t244) * t488 + (t35 + t36 + (-t109 - t97) * t288 + (t110 + t98) * t285) * t409 + (t37 + t38 + (-t111 - t99) * t288 + (t100 + t112) * t285) * t408; m(6) * (t285 * t42 - t288 * t41 + (t285 * t91 + t288 * t90) * qJD(1)) + m(5) * (t285 * t51 - t288 * t50 + (t103 * t288 + t104 * t285) * qJD(1)) + m(4) * (-t105 * t288 + t285 * t106 + (t174 * t288 + t175 * t285) * qJD(1)); m(6) * (t285 * t46 - t288 * t47 + (t107 * t285 + t108 * t288) * qJD(1)) + m(5) * (t285 * t66 - t288 * t67 + (t122 * t285 + t123 * t288) * qJD(1)) + m(4) * (t285 * t116 - t117 * t288 + (t195 * t285 + t196 * t288) * qJD(1)); 0; m(6) * (t24 * t91 + t25 * t90 + t41 * t65 + t42 * t64) + m(5) * (t101 * t51 + t102 * t50 + t103 * t44 + t104 * t43) + ((t285 * t311 + t288 * t310) * qJD(2) - t464) * t275 + ((t20 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1 + t22 / 0.2e1) * t288 + (t21 / 0.2e1 + t19 / 0.2e1 + t33 / 0.2e1 + t34 / 0.2e1) * t285 + (-t285 * t310 + t288 * t311) * qJD(1)) * t274 - t459; m(6) * (t10 * t45 + t107 * t24 + t108 * t25 + t46 * t64 + t47 * t65 + t48 * t9) + m(5) * (t101 * t66 + t102 * t67 + t122 * t43 + t123 * t44 + t23 * t92 + t30 * t49) + (-t4 / 0.2e1 - t3 / 0.2e1 + t390 * t406) * t288 + (t2 / 0.2e1 + t1 / 0.2e1 + t391 * t406) * t285 + ((-t285 * t390 + t288 * t391) * qJD(1) + (t7 + t8) * t481 + (t5 + t6) * t288 / 0.2e1 + t502 * qJD(2) / 0.2e1) * t274 - (t503 * qJD(1) + t500 * t285 + t501 * t288) * t275 / 0.2e1 + (t466 * t285 + t499 * t288) * qJD(1) / 0.2e1; m(5) * (t44 * t285 - t288 * t43 + (t101 * t288 + t102 * t285) * qJD(1)) + m(6) * (-t24 * t288 + t25 * t285 + (t285 * t65 + t288 * t64) * qJD(1)); (t10 * t48 + t24 * t65 + t25 * t64) * t485 + (t101 * t44 + t102 * t43 + t30 * t92) * t486 + (t464 * t275 + ((-t275 * t462 + t499) * t288 + (-t275 * t463 + t466) * t285) * qJD(2) + t459) * t275 + ((t1 + t2) * t288 + (t3 + t4) * t285 + (t501 * t285 - t500 * t288) * t275 + (t503 * t274 + t505 * t275) * qJD(2) + (t502 * t275 - t285 * t499 + t466 * t288) * qJD(1)) * t274; m(6) * (-t132 * t90 + t134 * t91 + t227 * t41 + t229 * t42); m(6) * (t45 * t378 + t107 * t134 - t108 * t132 + t227 * t47 + t229 * t46 + (t283 * t9 + t397 * t45) * t274); m(6) * (-t132 * t285 - t134 * t288 + (t227 * t285 + t229 * t288) * qJD(1)); m(6) * (t48 * t378 - t132 * t64 + t134 * t65 + t227 * t24 + t229 * t25 + (t10 * t283 + t397 * t48) * t274); (-t132 * t229 + t134 * t227 - t430 * t490) * t485;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t26(1), t26(2), t26(4), t26(7), t26(11); t26(2), t26(3), t26(5), t26(8), t26(12); t26(4), t26(5), t26(6), t26(9), t26(13); t26(7), t26(8), t26(9), t26(10), t26(14); t26(11), t26(12), t26(13), t26(14), t26(15);];
Mq = res;
