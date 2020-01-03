% Calculate time derivative of joint inertia matrix for
% S5RRPRP11
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP11_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP11_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP11_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:24
% EndTime: 2019-12-31 20:12:54
% DurationCPUTime: 19.38s
% Computational Cost: add. (10958->829), mult. (29872->1156), div. (0->0), fcn. (28487->6), ass. (0->380)
t291 = cos(qJ(2));
t288 = sin(qJ(2));
t443 = Icges(4,6) * t288;
t453 = Icges(3,4) * t288;
t523 = -t443 - t453 + (-Icges(3,2) - Icges(4,3)) * t291;
t442 = Icges(4,6) * t291;
t452 = Icges(3,4) * t291;
t522 = -t442 - t452 + (-Icges(3,1) - Icges(4,2)) * t288;
t287 = sin(qJ(4));
t290 = cos(qJ(4));
t445 = Icges(6,5) * t290;
t345 = Icges(6,1) * t287 - t445;
t197 = Icges(6,4) * t288 - t291 * t345;
t449 = Icges(5,4) * t290;
t346 = Icges(5,1) * t287 + t449;
t198 = Icges(5,5) * t288 - t291 * t346;
t506 = -t198 - t197;
t521 = t506 * t287;
t338 = Icges(5,5) * t287 + Icges(5,6) * t290;
t190 = Icges(5,3) * t288 - t291 * t338;
t340 = Icges(6,4) * t287 - Icges(6,6) * t290;
t193 = Icges(6,2) * t288 - t291 * t340;
t520 = t190 + t193;
t519 = t523 * qJD(2);
t518 = t522 * qJD(2);
t289 = sin(qJ(1));
t371 = qJD(4) * t288 + qJD(1);
t406 = qJD(2) * t291;
t382 = t290 * t406;
t292 = cos(qJ(1));
t409 = qJD(1) * t292;
t391 = t288 * t409;
t431 = t290 * t292;
t437 = t287 * t289;
t125 = -qJD(4) * t431 - t289 * t382 - t290 * t391 + t371 * t437;
t317 = t290 * t371;
t370 = qJD(1) * t288 + qJD(4);
t383 = t289 * t406;
t126 = t289 * t317 + (t292 * t370 + t383) * t287;
t434 = t289 * t290;
t225 = t287 * t292 + t288 * t434;
t455 = rSges(6,3) + qJ(5);
t509 = rSges(6,1) + pkin(4);
t517 = t225 * qJD(5) - t125 * t455 - t126 * t509;
t226 = t288 * t437 - t431;
t516 = -t225 * t455 + t226 * t509;
t344 = -Icges(3,2) * t288 + t452;
t195 = -Icges(3,6) * t292 + t289 * t344;
t334 = -Icges(4,3) * t288 + t442;
t488 = Icges(4,5) * t292 + t289 * t334;
t514 = -t195 - t488;
t196 = Icges(3,6) * t289 + t292 * t344;
t201 = Icges(4,5) * t289 - t292 * t334;
t513 = t196 - t201;
t348 = Icges(3,1) * t291 - t453;
t199 = -Icges(3,5) * t292 + t289 * t348;
t336 = Icges(4,2) * t291 - t443;
t487 = Icges(4,4) * t292 + t289 * t336;
t512 = t199 + t487;
t200 = Icges(3,5) * t289 + t292 * t348;
t203 = Icges(4,4) * t289 - t292 * t336;
t511 = t200 - t203;
t446 = Icges(6,5) * t287;
t337 = -Icges(6,3) * t290 + t446;
t189 = Icges(6,6) * t288 - t291 * t337;
t450 = Icges(5,4) * t287;
t341 = Icges(5,2) * t290 + t450;
t194 = Icges(5,6) * t288 - t291 * t341;
t508 = ((t189 - t194) * t290 + t521) * t291 + t520 * t288;
t403 = qJD(4) * t291;
t162 = (-Icges(6,1) * t290 - t446) * t403 + (Icges(6,4) * t291 + t288 * t345) * qJD(2);
t163 = (-Icges(5,1) * t290 + t450) * t403 + (Icges(5,5) * t291 + t288 * t346) * qJD(2);
t507 = -t163 - t162;
t154 = (-Icges(6,3) * t287 - t445) * t403 + (Icges(6,6) * t291 + t288 * t337) * qJD(2);
t155 = (-Icges(5,5) * t290 + Icges(5,6) * t287) * t403 + (Icges(5,3) * t291 + t288 * t338) * qJD(2);
t158 = (-Icges(6,4) * t290 - Icges(6,6) * t287) * t403 + (Icges(6,2) * t291 + t288 * t340) * qJD(2);
t380 = t287 * t403;
t408 = qJD(2) * t288;
t385 = t290 * t408;
t432 = t290 * t291;
t505 = t154 * t432 - t189 * t385 + (t380 + t385) * t194 + t520 * t406 - t408 * t521 + (t158 + t155) * t288;
t433 = t289 * t291;
t425 = rSges(6,2) * t433 + t516;
t373 = t425 * t292;
t223 = -t288 * t431 + t437;
t435 = t288 * t292;
t224 = t287 * t435 + t434;
t430 = t291 * t292;
t426 = rSges(6,2) * t430 + t223 * t455 + t224 * t509;
t504 = -t289 * t426 + t373;
t405 = qJD(2) * t292;
t384 = t288 * t405;
t410 = qJD(1) * t289;
t503 = t291 * t410 + t384;
t140 = Icges(6,4) * t226 + Icges(6,2) * t433 - Icges(6,6) * t225;
t136 = Icges(6,5) * t226 + Icges(6,6) * t433 - Icges(6,3) * t225;
t144 = Icges(6,1) * t226 + Icges(6,4) * t433 - Icges(6,5) * t225;
t331 = t136 * t290 - t144 * t287;
t62 = t140 * t288 + t291 * t331;
t138 = Icges(5,5) * t226 + Icges(5,6) * t225 + Icges(5,3) * t433;
t142 = Icges(5,4) * t226 + Icges(5,2) * t225 + Icges(5,6) * t433;
t146 = Icges(5,1) * t226 + Icges(5,4) * t225 + Icges(5,5) * t433;
t329 = t142 * t290 + t146 * t287;
t64 = t138 * t288 - t291 * t329;
t464 = t62 + t64;
t139 = Icges(6,4) * t224 + Icges(6,2) * t430 + Icges(6,6) * t223;
t135 = Icges(6,5) * t224 + Icges(6,6) * t430 + Icges(6,3) * t223;
t143 = Icges(6,1) * t224 + Icges(6,4) * t430 + Icges(6,5) * t223;
t332 = t135 * t290 - t143 * t287;
t61 = t139 * t288 + t291 * t332;
t137 = Icges(5,5) * t224 - Icges(5,6) * t223 + Icges(5,3) * t430;
t141 = Icges(5,4) * t224 - Icges(5,2) * t223 + Icges(5,6) * t430;
t145 = Icges(5,1) * t224 - Icges(5,4) * t223 + Icges(5,5) * t430;
t330 = t141 * t290 + t145 * t287;
t63 = t137 * t288 - t291 * t330;
t465 = t61 + t63;
t502 = t289 * t465 - t292 * t464;
t501 = t289 * t464 + t292 * t465;
t500 = qJD(2) / 0.2e1;
t381 = t291 * t405;
t127 = qJD(1) * t225 + qJD(4) * t224 - t290 * t381;
t128 = t292 * t317 + (-t289 * t370 + t381) * t287;
t70 = Icges(6,5) * t128 - Icges(6,6) * t503 + Icges(6,3) * t127;
t74 = Icges(6,4) * t128 - Icges(6,2) * t503 + Icges(6,6) * t127;
t78 = Icges(6,1) * t128 - Icges(6,4) * t503 + Icges(6,5) * t127;
t19 = (-qJD(2) * t332 + t74) * t288 + (qJD(2) * t139 - t287 * t78 + t290 * t70 + (-t135 * t287 - t143 * t290) * qJD(4)) * t291;
t72 = Icges(5,5) * t128 - Icges(5,6) * t127 - Icges(5,3) * t503;
t76 = Icges(5,4) * t128 - Icges(5,2) * t127 - Icges(5,6) * t503;
t80 = Icges(5,1) * t128 - Icges(5,4) * t127 - Icges(5,5) * t503;
t21 = (qJD(2) * t330 + t72) * t288 + (qJD(2) * t137 - t287 * t80 - t290 * t76 + (t141 * t287 - t145 * t290) * qJD(4)) * t291;
t499 = t19 + t21;
t407 = qJD(2) * t289;
t386 = t288 * t407;
t388 = t291 * t409;
t300 = -t386 + t388;
t69 = Icges(6,5) * t126 + Icges(6,6) * t300 + Icges(6,3) * t125;
t73 = Icges(6,4) * t126 + Icges(6,2) * t300 + Icges(6,6) * t125;
t77 = Icges(6,1) * t126 + Icges(6,4) * t300 + Icges(6,5) * t125;
t20 = (-qJD(2) * t331 + t73) * t288 + (qJD(2) * t140 - t287 * t77 + t290 * t69 + (-t136 * t287 - t144 * t290) * qJD(4)) * t291;
t71 = Icges(5,5) * t126 - Icges(5,6) * t125 + Icges(5,3) * t300;
t75 = Icges(5,4) * t126 - Icges(5,2) * t125 + Icges(5,6) * t300;
t79 = Icges(5,1) * t126 - Icges(5,4) * t125 + Icges(5,5) * t300;
t22 = (qJD(2) * t329 + t71) * t288 + (qJD(2) * t138 - t287 * t79 - t290 * t75 + (t142 * t287 - t146 * t290) * qJD(4)) * t291;
t498 = t20 + t22;
t54 = t137 * t430 - t141 * t223 + t145 * t224;
t55 = t138 * t430 - t142 * t223 + t146 * t224;
t352 = t289 * t55 + t292 * t54;
t52 = t135 * t223 + t139 * t430 + t143 * t224;
t53 = t136 * t223 + t140 * t430 + t144 * t224;
t353 = t289 * t53 + t292 * t52;
t89 = t189 * t223 + t193 * t430 + t197 * t224;
t90 = t190 * t430 - t194 * t223 + t198 * t224;
t497 = (t352 + t353) * t291 + (t89 + t90) * t288;
t58 = t137 * t433 + t141 * t225 + t145 * t226;
t59 = t138 * t433 + t142 * t225 + t146 * t226;
t350 = t289 * t59 + t292 * t58;
t56 = -t135 * t225 + t139 * t433 + t143 * t226;
t57 = -t136 * t225 + t140 * t433 + t144 * t226;
t351 = t289 * t57 + t292 * t56;
t91 = -t189 * t225 + t193 * t433 + t197 * t226;
t92 = t190 * t433 + t194 * t225 + t198 * t226;
t466 = (t350 + t351) * t291 + (t91 + t92) * t288;
t496 = t223 * qJD(5) + t127 * t455 + t128 * t509;
t320 = t201 * t288 - t203 * t291;
t495 = t289 * t320;
t321 = t196 * t288 - t200 * t291;
t494 = t289 * t321;
t319 = -t288 * t488 + t291 * t487;
t492 = t292 * t319;
t322 = t195 * t288 - t199 * t291;
t491 = t292 * t322;
t490 = rSges(4,1) * t289 - rSges(4,2) * t430;
t241 = pkin(3) * t289 + pkin(7) * t430;
t489 = -rSges(3,2) * t435 + rSges(3,3) * t289;
t339 = Icges(3,5) * t291 - Icges(3,6) * t288;
t191 = -Icges(3,3) * t292 + t289 * t339;
t342 = Icges(4,4) * t291 - Icges(4,5) * t288;
t486 = Icges(4,1) * t292 + t289 * t342;
t485 = -t289 * t425 - t292 * t426;
t484 = 2 * m(3);
t483 = 2 * m(4);
t482 = 2 * m(5);
t481 = 2 * m(6);
t480 = m(4) / 0.2e1;
t479 = m(5) / 0.2e1;
t478 = m(6) / 0.2e1;
t477 = -pkin(2) - pkin(7);
t475 = t289 / 0.2e1;
t474 = -t292 / 0.2e1;
t254 = rSges(3,1) * t288 + rSges(3,2) * t291;
t471 = m(3) * t254;
t470 = pkin(2) * t291;
t463 = rSges(6,2) * t300 - t517;
t462 = -rSges(6,2) * t503 + t496;
t461 = rSges(4,1) * t292;
t460 = rSges(4,2) * t288;
t459 = rSges(6,2) * t288;
t458 = rSges(3,3) * t292;
t457 = rSges(5,3) * t288;
t456 = -rSges(4,3) - qJ(3);
t440 = qJ(3) * t288;
t439 = qJ(3) * t291;
t359 = -rSges(5,1) * t226 - rSges(5,2) * t225;
t151 = rSges(5,3) * t433 - t359;
t438 = t151 * t292;
t428 = rSges(5,1) * t128 - rSges(5,2) * t127;
t357 = rSges(6,1) * t287 - rSges(6,3) * t290;
t427 = (pkin(4) * t408 - qJ(5) * t403) * t287 + (-qJ(5) * t408 + (-pkin(4) * qJD(4) + qJD(5)) * t291) * t290 + (-rSges(6,1) * t290 - rSges(6,3) * t287) * t403 + (rSges(6,2) * t291 + t288 * t357) * qJD(2);
t424 = t459 + (-pkin(4) * t287 + qJ(5) * t290 - t357) * t291;
t354 = t440 + t470;
t229 = t354 * t289;
t231 = pkin(2) * t430 + qJ(3) * t435;
t423 = t229 * t289 + t231 * t292;
t216 = qJD(2) * t354 - qJD(3) * t291;
t356 = -rSges(4,2) * t291 + rSges(4,3) * t288;
t422 = -qJD(2) * t356 - t216;
t421 = -t231 - t241;
t252 = pkin(2) * t288 - t439;
t232 = t252 * t410;
t392 = t288 * t410;
t420 = pkin(7) * t392 + t232;
t355 = rSges(4,3) * t291 + t460;
t419 = -t252 + t355;
t404 = qJD(3) * t288;
t418 = qJ(3) * t381 + t292 * t404;
t417 = rSges(3,2) * t392 + rSges(3,3) * t409;
t261 = pkin(7) * t386;
t262 = pkin(2) * t386;
t416 = t261 + t262;
t415 = pkin(1) * t292 + pkin(6) * t289;
t281 = t292 * pkin(6);
t282 = t292 * pkin(3);
t414 = t281 + t282;
t413 = t289 ^ 2 + t292 ^ 2;
t192 = Icges(3,3) * t289 + t292 * t339;
t412 = qJD(1) * t192;
t205 = Icges(4,1) * t289 - t292 * t342;
t411 = qJD(1) * t205;
t401 = -rSges(6,2) + t477;
t400 = -rSges(5,3) + t477;
t35 = t289 * t52 - t292 * t53;
t36 = t289 * t54 - t292 * t55;
t397 = -t36 / 0.2e1 - t35 / 0.2e1;
t37 = t289 * t56 - t292 * t57;
t38 = t289 * t58 - t292 * t59;
t396 = t37 / 0.2e1 + t38 / 0.2e1;
t395 = (-pkin(3) - pkin(6)) * t289;
t394 = t289 * (pkin(2) * t388 + t289 * t404 - t262 + (t383 + t391) * qJ(3)) + t292 * (-pkin(2) * t503 - qJ(3) * t392 + t418) + t229 * t409;
t149 = rSges(5,1) * t224 - rSges(5,2) * t223 + rSges(5,3) * t430;
t275 = pkin(6) * t409;
t393 = t275 + t418;
t358 = rSges(5,1) * t287 + rSges(5,2) * t290;
t210 = -t291 * t358 + t457;
t390 = t210 * t410;
t378 = t291 ^ 2 * qJD(4) * t287;
t377 = -t342 * qJD(2) / 0.2e1 + t339 * t500;
t376 = -pkin(1) - t440;
t375 = -pkin(7) * t288 - t252;
t159 = (Icges(5,2) * t287 - t449) * t403 + (Icges(5,6) * t291 + t288 * t341) * qJD(2);
t374 = t508 * t406 + (((qJD(4) * t506 - t159) * t290 + (-qJD(4) * t189 + t507) * t287) * t291 + t505) * t288;
t184 = t419 * t292;
t372 = qJD(1) * t424;
t242 = pkin(7) * t433 - t282;
t369 = t241 * t292 + t242 * t289 + t423;
t276 = pkin(3) * t409;
t368 = t276 + t393;
t367 = rSges(4,1) * t409 + rSges(4,2) * t503 + rSges(4,3) * t381;
t366 = t415 + t231;
t365 = -t210 + t375;
t364 = -pkin(7) * t406 - t216;
t363 = t289 * t372;
t362 = t288 * t456 - pkin(1);
t361 = rSges(3,1) * t291 - rSges(3,2) * t288;
t360 = rSges(5,1) * t126 - rSges(5,2) * t125;
t328 = -t149 * t292 - t151 * t289;
t327 = t149 * t289 - t438;
t318 = t375 - t424;
t211 = rSges(3,1) * t430 + t489;
t212 = rSges(4,3) * t435 + t490;
t173 = (-rSges(5,1) * t290 + rSges(5,2) * t287) * t403 + (rSges(5,3) * t291 + t288 * t358) * qJD(2);
t316 = -t173 + t364;
t315 = -pkin(1) - t361;
t33 = t127 * t189 + t128 * t197 + t154 * t223 + t158 * t430 + t162 * t224 - t193 * t503;
t34 = -t127 * t194 + t128 * t198 + t155 * t430 - t159 * t223 + t163 * t224 - t190 * t503;
t314 = t21 / 0.2e1 + t19 / 0.2e1 + t33 / 0.2e1 + t34 / 0.2e1;
t31 = t125 * t189 + t126 * t197 - t154 * t225 + t158 * t433 + t162 * t226 + t193 * t300;
t32 = -t125 * t194 + t126 * t198 + t155 * t433 + t159 * t225 + t163 * t226 + t190 * t300;
t313 = t22 / 0.2e1 + t20 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1;
t312 = -t63 / 0.2e1 - t61 / 0.2e1 - t89 / 0.2e1 - t90 / 0.2e1;
t311 = -t92 / 0.2e1 - t64 / 0.2e1 - t62 / 0.2e1 - t91 / 0.2e1;
t134 = t365 * t292;
t310 = t289 * (qJD(1) * t241 - t261) + t292 * (-pkin(7) * t503 + t276) + t242 * t409 + t394;
t309 = t366 + t241;
t308 = t364 - t427;
t307 = qJD(2) * t254;
t304 = qJD(2) * (Icges(4,4) * t288 + Icges(4,5) * t291);
t303 = qJD(2) * (-Icges(3,5) * t288 - Icges(3,6) * t291);
t116 = t318 * t292;
t298 = t291 * t401 + t376;
t297 = t291 * t400 + t376;
t296 = (rSges(4,2) - pkin(2)) * t291 + t362;
t294 = t298 * t289;
t293 = t297 * t289;
t240 = t361 * qJD(2);
t213 = t289 * t356 - t461;
t208 = t289 * t361 - t458;
t183 = t419 * t289;
t178 = t211 + t415;
t177 = t289 * t315 + t281 + t458;
t171 = qJD(1) * t486 + t292 * t304;
t170 = t289 * t304 + t411;
t157 = t289 * t303 + t412;
t156 = -qJD(1) * t191 + t292 * t303;
t153 = t212 + t366;
t152 = t289 * t296 + t281 + t461;
t133 = t365 * t289;
t118 = t254 * t407 + ((-rSges(3,3) - pkin(6)) * t289 + t315 * t292) * qJD(1);
t117 = -rSges(3,1) * t503 - rSges(3,2) * t381 - pkin(1) * t410 + t275 + t417;
t115 = t318 * t289;
t114 = qJD(1) * t184 + t289 * t422;
t113 = t292 * t422 - t355 * t410 + t232;
t112 = t149 * t288 - t210 * t430;
t111 = -t151 * t288 + t210 * t433;
t110 = t289 * t319 + t292 * t486;
t109 = -t205 * t292 + t495;
t108 = -t289 * t486 + t492;
t107 = t205 * t289 + t292 * t320;
t106 = t192 * t289 - t292 * t321;
t105 = t191 * t289 - t491;
t102 = -t192 * t292 - t494;
t101 = -t191 * t292 - t289 * t322;
t100 = t309 + t149;
t99 = t293 + t359 + t414;
t96 = t212 * t292 + t213 * t289 + t423;
t95 = t262 + (-t404 + (t291 * t456 - t460) * qJD(2)) * t289 + ((-rSges(4,1) - pkin(6)) * t289 + t296 * t292) * qJD(1);
t94 = -pkin(2) * t384 + (t362 - t470) * t410 + t367 + t393;
t93 = t327 * t291;
t88 = t309 + t426;
t87 = t294 + t414 - t516;
t86 = -rSges(5,3) * t503 + t428;
t84 = rSges(5,3) * t300 + t360;
t68 = qJD(1) * t134 + t289 * t316;
t67 = t292 * t316 + t390 + t420;
t66 = t288 * t426 - t424 * t430;
t65 = -t288 * t425 + t424 * t433;
t60 = -t328 + t369;
t51 = t504 * t291;
t50 = qJD(1) * t116 + t289 * t308;
t49 = t292 * t308 + t363 + t420;
t48 = (-t404 + (-t439 + t457) * qJD(2)) * t289 + (t292 * t297 + t395) * qJD(1) - t360 + t416;
t47 = qJD(1) * t293 + t384 * t400 + t368 + t428;
t46 = t369 - t485;
t45 = (-t210 * t407 - t84) * t288 + (-qJD(2) * t151 + t173 * t289 + t210 * t409) * t291;
t44 = (t210 * t405 + t86) * t288 + (qJD(2) * t149 - t173 * t292 + t390) * t291;
t43 = (qJD(1) * t213 + t367) * t292 + (t355 * t407 + (-t212 - t231 + t490) * qJD(1)) * t289 + t394;
t40 = (-t404 + (-t439 + t459) * qJD(2)) * t289 + (t292 * t298 + t395) * qJD(1) + t416 + t517;
t39 = qJD(1) * t294 + t384 * t401 + t368 + t496;
t30 = t327 * t408 + (qJD(1) * t328 - t289 * t86 + t292 * t84) * t291;
t25 = (-t407 * t424 - t463) * t288 + (-qJD(2) * t425 + t289 * t427 + t292 * t372) * t291;
t24 = (t405 * t424 + t462) * t288 + (qJD(2) * t426 - t292 * t427 + t363) * t291;
t23 = t289 * t84 + t292 * t86 + (t438 + (-t149 + t421) * t289) * qJD(1) + t310;
t18 = -t127 * t142 + t128 * t146 - t138 * t503 - t223 * t75 + t224 * t79 + t430 * t71;
t17 = -t127 * t141 + t128 * t145 - t137 * t503 - t223 * t76 + t224 * t80 + t430 * t72;
t16 = t127 * t136 + t128 * t144 - t140 * t503 + t223 * t69 + t224 * t77 + t430 * t73;
t15 = t127 * t135 + t128 * t143 - t139 * t503 + t223 * t70 + t224 * t78 + t430 * t74;
t14 = -t125 * t142 + t126 * t146 + t138 * t300 + t225 * t75 + t226 * t79 + t433 * t71;
t13 = -t125 * t141 + t126 * t145 + t137 * t300 + t225 * t76 + t226 * t80 + t433 * t72;
t12 = t125 * t136 + t126 * t144 + t140 * t300 - t225 * t69 + t226 * t77 + t433 * t73;
t11 = t125 * t135 + t126 * t143 + t139 * t300 - t225 * t70 + t226 * t78 + t433 * t74;
t10 = -t504 * t408 + (qJD(1) * t485 - t462 * t289 + t463 * t292) * t291;
t9 = t462 * t292 + t463 * t289 + (t373 + (t421 - t426) * t289) * qJD(1) + t310;
t8 = qJD(1) * t352 + t17 * t289 - t18 * t292;
t7 = qJD(1) * t353 + t15 * t289 - t16 * t292;
t6 = qJD(1) * t350 + t13 * t289 - t14 * t292;
t5 = qJD(1) * t351 + t11 * t289 - t12 * t292;
t4 = (-qJD(2) * t352 + t34) * t288 + (-qJD(1) * t36 + qJD(2) * t90 + t17 * t292 + t18 * t289) * t291;
t3 = (-qJD(2) * t353 + t33) * t288 + (-qJD(1) * t35 + qJD(2) * t89 + t15 * t292 + t16 * t289) * t291;
t2 = (-qJD(2) * t350 + t32) * t288 + (-qJD(1) * t38 + qJD(2) * t92 + t13 * t292 + t14 * t289) * t291;
t1 = (-qJD(2) * t351 + t31) * t288 + (-qJD(1) * t37 + qJD(2) * t91 + t11 * t292 + t12 * t289) * t291;
t26 = [-t189 * t380 - t159 * t432 + (t117 * t178 + t118 * t177) * t484 + (t152 * t95 + t153 * t94) * t483 + (t100 * t47 + t48 * t99) * t482 + (t39 * t88 + t40 * t87) * t481 + t507 * t287 * t291 + t506 * t290 * t403 + (t348 + t336 + t523) * t408 + (t344 + t334 - t522) * t406 + t505; m(4) * (t113 * t152 + t114 * t153 + t183 * t94 + t184 * t95) + m(5) * (t100 * t68 + t133 * t47 + t134 * t48 + t67 * t99) + m(6) * (t115 * t39 + t116 * t40 + t49 * t87 + t50 * t88) + (m(3) * (-t118 * t254 - t177 * t240) + t377 * t292 - t313) * t292 + (m(3) * (-t117 * t254 - t178 * t240) + t377 * t289 + t314) * t289 + ((-t513 * qJD(2) + t292 * t518) * t475 + (t514 * qJD(2) + t289 * t518) * t474 + (t474 * t511 - t475 * t512) * qJD(1)) * t288 + ((t511 * qJD(2) + t292 * t519) * t475 + (t512 * qJD(2) + t289 * t519) * t474 + (t474 * t513 + t475 * t514) * qJD(1)) * t291 + ((-t178 * t471 + (t196 / 0.2e1 - t201 / 0.2e1) * t291 + (t200 / 0.2e1 - t203 / 0.2e1) * t288 - t312) * t292 + (t177 * t471 + (t195 / 0.2e1 + t488 / 0.2e1) * t291 + (t199 / 0.2e1 + t487 / 0.2e1) * t288 - t311) * t289) * qJD(1); (t115 * t50 + t116 * t49 + t46 * t9) * t481 + t289 * t8 - t292 * t5 + (t133 * t68 + t134 * t67 + t23 * t60) * t482 + t289 * t7 - t292 * t6 + (t113 * t184 + t114 * t183 + t43 * t96) * t483 - t292 * ((t170 * t292 + (t109 - t492) * qJD(1)) * t292 + (t110 * qJD(1) + (t201 * t406 + t203 * t408 + t411) * t289 + (-t171 + (t288 * t487 + t291 * t488) * qJD(2) + t320 * qJD(1)) * t292) * t289) - t292 * ((t157 * t292 + (t102 + t491) * qJD(1)) * t292 + (t101 * qJD(1) + (-t196 * t406 - t200 * t408 + t412) * t289 + (-t156 + (t195 * t291 + t199 * t288) * qJD(2) - t321 * qJD(1)) * t292) * t289) + t289 * ((t289 * t171 + (t108 - t495) * qJD(1)) * t289 + (t107 * qJD(1) + (t406 * t488 + t408 * t487) * t292 + (-t170 + (t201 * t291 + t203 * t288) * qJD(2) + (t205 + t319) * qJD(1)) * t289) * t292) + ((t208 * t289 + t211 * t292) * ((qJD(1) * t208 - t292 * t307 + t417) * t292 + (-t289 * t307 + (-t211 + t489) * qJD(1)) * t289) + t413 * t254 * t240) * t484 + t289 * ((t289 * t156 + (t105 + t494) * qJD(1)) * t289 + (t106 * qJD(1) + (t195 * t406 + t199 * t408) * t292 + (-t157 + (-t196 * t291 - t200 * t288) * qJD(2) + (t192 - t322) * qJD(1)) * t289) * t292) + (t37 + t38 + (-t101 - t110) * t292 + (t102 + t109) * t289) * t410 + (t35 + t36 + (-t105 - t108) * t292 + (t106 + t107) * t289) * t409; 0.2e1 * ((t100 * t289 + t292 * t99) * t479 + (t289 * t88 + t292 * t87) * t478 + (t152 * t292 + t153 * t289) * t480) * t406 + 0.2e1 * ((t100 * t409 + t289 * t47 + t292 * t48 - t410 * t99) * t479 + (t289 * t39 + t292 * t40 + t409 * t88 - t410 * t87) * t478 + (-t152 * t410 + t153 * t409 + t289 * t94 + t292 * t95) * t480) * t288; 0.2e1 * ((t115 * t407 + t116 * t405 - t9) * t478 + (t133 * t407 + t134 * t405 - t23) * t479 + (t183 * t407 + t184 * t405 - t43) * t480) * t291 + 0.2e1 * ((qJD(2) * t46 + t115 * t409 - t116 * t410 + t289 * t50 + t292 * t49) * t478 + (qJD(2) * t60 + t133 * t409 - t134 * t410 + t289 * t68 + t292 * t67) * t479 + (qJD(2) * t96 + t113 * t292 + t114 * t289 + t183 * t409 - t184 * t410) * t480) * t288; 0.4e1 * (t480 + t479 + t478) * (-0.1e1 + t413) * t288 * t406; m(5) * (t100 * t44 + t111 * t48 + t112 * t47 + t45 * t99) + m(6) * (t24 * t88 + t25 * t87 + t39 * t66 + t40 * t65) + (t289 * t311 + t292 * t312) * t408 + (t314 * t292 + t313 * t289 + (t289 * t312 - t292 * t311) * qJD(1)) * t291 + t374; m(5) * (t111 * t67 + t112 * t68 + t133 * t44 + t134 * t45 - t23 * t93 + t30 * t60) + m(6) * (t10 * t46 + t115 * t24 + t116 * t25 + t49 * t65 + t50 * t66 + t51 * t9) + (-t2 / 0.2e1 - t1 / 0.2e1 + t397 * t408) * t292 + (t4 / 0.2e1 + t3 / 0.2e1 - t396 * t408) * t289 + ((t289 * t397 + t292 * t396) * qJD(1) + (t5 + t6) * t475 + (t7 + t8) * t292 / 0.2e1 + t502 * t500) * t291 + (qJD(1) * t501 + t499 * t289 - t498 * t292) * t288 / 0.2e1 + (t466 * t289 + t497 * t292) * qJD(1) / 0.2e1; 0.2e1 * ((t111 * t405 + t112 * t407 - t30) * t479 + (t405 * t65 + t407 * t66 - t10) * t478) * t291 + 0.2e1 * ((-qJD(2) * t93 - t111 * t410 + t112 * t409 + t289 * t44 + t292 * t45) * t479 + (qJD(2) * t51 + t24 * t289 + t25 * t292 + t409 * t66 - t410 * t65) * t478) * t288; (t10 * t51 + t24 * t66 + t25 * t65) * t481 + (t111 * t45 + t112 * t44 - t30 * t93) * t482 + (((-t288 * t465 - t497) * t292 + (-t288 * t464 - t466) * t289) * qJD(2) + t374) * t288 + ((t3 + t4) * t292 + (t1 + t2) * t289 + (t498 * t289 + t499 * t292) * t288 + (t288 * t508 + t291 * t501) * qJD(2) + (-t288 * t502 - t289 * t497 + t466 * t292) * qJD(1)) * t291; m(6) * (t125 * t88 + t127 * t87 + t223 * t40 - t225 * t39); m(6) * (-t46 * t380 + t115 * t125 + t116 * t127 + t223 * t49 - t225 * t50 + (t291 * t9 - t408 * t46) * t290); m(6) * (t378 + (t223 * t292 - t225 * t289) * t406 + (0.2e1 * t382 + t125 * t289 + t127 * t292 + (-t223 * t289 - t225 * t292) * qJD(1)) * t288); m(6) * (-t51 * t380 + t125 * t66 + t127 * t65 + t223 * t25 - t225 * t24 + (t10 * t291 - t408 * t51) * t290); (-t125 * t225 + t127 * t223 + (-t288 * t382 - t378) * t290) * t481;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t26(1), t26(2), t26(4), t26(7), t26(11); t26(2), t26(3), t26(5), t26(8), t26(12); t26(4), t26(5), t26(6), t26(9), t26(13); t26(7), t26(8), t26(9), t26(10), t26(14); t26(11), t26(12), t26(13), t26(14), t26(15);];
Mq = res;
