% Calculate time derivative of joint inertia matrix for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR11_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR11_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR11_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:27
% EndTime: 2019-12-31 19:46:56
% DurationCPUTime: 16.44s
% Computational Cost: add. (11126->755), mult. (19495->1076), div. (0->0), fcn. (17967->8), ass. (0->381)
t270 = cos(qJ(2));
t268 = sin(qJ(2));
t435 = Icges(4,6) * t268;
t443 = Icges(3,4) * t268;
t502 = t435 + t443 + (Icges(3,2) + Icges(4,3)) * t270;
t434 = Icges(4,6) * t270;
t442 = Icges(3,4) * t270;
t501 = t434 + t442 + (Icges(3,1) + Icges(4,2)) * t268;
t500 = -Icges(5,3) / 0.2e1;
t499 = t502 * qJD(2);
t498 = t501 * qJD(2);
t266 = cos(pkin(8));
t265 = sin(pkin(8));
t269 = sin(qJ(1));
t422 = t269 * t265;
t271 = cos(qJ(1));
t425 = t268 * t271;
t194 = t266 * t425 - t422;
t399 = qJD(2) * t270;
t371 = t269 * t399;
t141 = qJD(1) * t194 + t266 * t371;
t497 = -t141 / 0.2e1;
t380 = t265 * t425;
t421 = t269 * t266;
t195 = t380 + t421;
t142 = qJD(1) * t195 + t265 * t371;
t496 = -t142 / 0.2e1;
t495 = t269 / 0.2e1;
t494 = -t271 / 0.2e1;
t493 = -qJD(1) / 0.2e1;
t454 = qJD(1) / 0.2e1;
t398 = qJD(2) * t271;
t372 = t268 * t398;
t403 = qJD(1) * t269;
t492 = t270 * t403 + t372;
t260 = pkin(8) + qJ(5);
t252 = cos(t260);
t251 = sin(t260);
t440 = Icges(6,4) * t251;
t323 = Icges(6,2) * t252 + t440;
t155 = Icges(6,6) * t268 - t270 * t323;
t439 = Icges(6,4) * t252;
t328 = Icges(6,1) * t251 + t439;
t156 = Icges(6,5) * t268 - t270 * t328;
t491 = t155 * t252 + t156 * t251;
t406 = t269 ^ 2 + t271 ^ 2;
t490 = -0.1e1 + t406;
t489 = qJD(2) / 0.2e1;
t246 = pkin(4) * t266 + pkin(3);
t267 = -pkin(7) - qJ(4);
t419 = t270 * t271;
t379 = t267 * t419;
t424 = t269 * t251;
t166 = t252 * t425 - t424;
t423 = t269 * t252;
t167 = t251 * t425 + t423;
t99 = t167 * rSges(6,1) + t166 * rSges(6,2) + rSges(6,3) * t419;
t488 = pkin(4) * t380 + t269 * t246 - t379 + t99;
t317 = -Icges(4,3) * t268 + t434;
t176 = Icges(4,5) * t269 - t271 * t317;
t319 = Icges(4,2) * t270 - t435;
t178 = Icges(4,4) * t269 - t271 * t319;
t309 = t176 * t268 - t178 * t270;
t487 = t269 * t309;
t327 = -Icges(3,2) * t268 + t442;
t173 = Icges(3,6) * t269 + t271 * t327;
t331 = Icges(3,1) * t270 - t443;
t175 = Icges(3,5) * t269 + t271 * t331;
t310 = t173 * t268 - t175 * t270;
t486 = t269 * t310;
t479 = Icges(4,4) * t271 + t269 * t319;
t480 = Icges(4,5) * t271 + t269 * t317;
t308 = -t268 * t480 + t270 * t479;
t485 = t271 * t308;
t172 = -Icges(3,6) * t271 + t269 * t327;
t174 = -Icges(3,5) * t271 + t269 * t331;
t311 = t172 * t268 - t174 * t270;
t484 = t271 * t311;
t483 = t269 * rSges(4,1) - rSges(4,2) * t419;
t211 = t269 * pkin(3) + qJ(4) * t419;
t482 = -rSges(3,2) * t425 + t269 * rSges(3,3);
t168 = t251 * t271 + t268 * t423;
t169 = -t252 * t271 + t268 * t424;
t346 = -t169 * rSges(6,1) - t168 * rSges(6,2);
t420 = t269 * t270;
t100 = rSges(6,3) * t420 - t346;
t332 = t271 * t100 - t269 * t99;
t401 = qJD(2) * t268;
t400 = qJD(2) * t269;
t373 = t268 * t400;
t402 = qJD(1) * t271;
t374 = t270 * t402;
t290 = -t373 + t374;
t363 = qJD(1) * t268 + qJD(5);
t277 = t271 * t363 + t371;
t364 = qJD(5) * t268 + qJD(1);
t307 = t251 * t364;
t86 = t252 * t277 - t269 * t307;
t306 = t252 * t364;
t87 = t251 * t277 + t269 * t306;
t356 = t87 * rSges(6,1) + t86 * rSges(6,2);
t48 = rSges(6,3) * t290 + t356;
t370 = t270 * t398;
t276 = -t269 * t363 + t370;
t88 = t252 * t276 - t271 * t307;
t89 = t251 * t276 + t271 * t306;
t452 = t89 * rSges(6,1) + t88 * rSges(6,2);
t49 = -rSges(6,3) * t492 + t452;
t14 = -t332 * t401 + (-t269 * t49 + t271 * t48 + (-t269 * t100 - t271 * t99) * qJD(1)) * t270;
t345 = rSges(6,1) * t251 + rSges(6,2) * t252;
t158 = rSges(6,3) * t268 - t270 * t345;
t61 = -t100 * t268 + t158 * t420;
t62 = -t158 * t419 + t268 * t99;
t481 = qJD(2) * (t269 * t62 + t271 * t61) - t14;
t322 = Icges(3,5) * t270 - Icges(3,6) * t268;
t170 = -Icges(3,3) * t271 + t269 * t322;
t325 = Icges(4,4) * t270 - Icges(4,5) * t268;
t478 = Icges(4,1) * t271 + t269 * t325;
t477 = 2 * m(3);
t476 = 2 * m(4);
t475 = 2 * m(5);
t474 = 2 * m(6);
t473 = m(4) / 0.2e1;
t472 = -m(5) / 0.2e1;
t471 = m(5) / 0.2e1;
t470 = -m(6) / 0.2e1;
t469 = m(6) / 0.2e1;
t324 = Icges(5,4) * t265 + Icges(5,2) * t266;
t163 = Icges(5,6) * t268 - t270 * t324;
t468 = t163 / 0.2e1;
t329 = Icges(5,1) * t265 + Icges(5,4) * t266;
t164 = Icges(5,5) * t268 - t270 * t329;
t467 = t164 / 0.2e1;
t466 = t194 / 0.2e1;
t465 = t195 / 0.2e1;
t464 = -t265 / 0.2e1;
t463 = t265 / 0.2e1;
t462 = -t266 / 0.2e1;
t461 = t266 / 0.2e1;
t460 = t268 / 0.2e1;
t459 = t271 / 0.2e1;
t458 = -rSges(6,3) - pkin(2);
t226 = rSges(3,1) * t268 + rSges(3,2) * t270;
t457 = m(3) * t226;
t456 = pkin(2) * t270;
t455 = pkin(4) * t265;
t395 = qJD(5) * t270;
t102 = (Icges(6,2) * t251 - t439) * t395 + (Icges(6,6) * t270 + t268 * t323) * qJD(2);
t320 = Icges(6,5) * t251 + Icges(6,6) * t252;
t101 = (-Icges(6,5) * t252 + Icges(6,6) * t251) * t395 + (Icges(6,3) * t270 + t268 * t320) * qJD(2);
t154 = Icges(6,3) * t268 - t270 * t320;
t351 = t251 * t155 * t395 + t268 * t101 + t154 * t399 + t401 * t491;
t103 = (-Icges(6,1) * t252 + t440) * t395 + (Icges(6,5) * t270 + t268 * t328) * qJD(2);
t427 = t251 * t103;
t58 = t154 * t268 - t270 * t491;
t453 = ((-t427 + (-qJD(5) * t156 - t102) * t252) * t270 + t351) * t268 + t58 * t399;
t451 = rSges(4,1) * t271;
t450 = rSges(4,2) * t268;
t449 = rSges(3,3) * t271;
t448 = rSges(5,3) * t268;
t447 = -rSges(4,3) - qJ(3);
t446 = rSges(6,3) - t267;
t445 = -t211 + t488;
t432 = qJ(3) * t268;
t431 = qJ(3) * t270;
t428 = t246 * t271;
t426 = t265 * t270;
t418 = -qJ(4) - t267;
t258 = t271 * pkin(3);
t291 = t268 * t455 + t270 * t418;
t417 = t269 * t291 + t100 + t258 - t428;
t196 = t265 * t271 + t268 * t421;
t143 = -qJD(1) * t196 + t266 * t370;
t197 = -t266 * t271 + t268 * t422;
t357 = t265 * t370;
t144 = -qJD(1) * t197 + t357;
t416 = t144 * rSges(5,1) + t143 * rSges(5,2);
t391 = pkin(4) * t426;
t415 = t268 * t418 + t158 - t391;
t342 = t432 + t456;
t200 = t342 * t269;
t201 = pkin(2) * t419 + qJ(3) * t425;
t414 = t269 * t200 + t271 * t201;
t192 = qJD(2) * t342 - qJD(3) * t270;
t344 = -rSges(4,2) * t270 + rSges(4,3) * t268;
t413 = -t344 * qJD(2) - t192;
t412 = -t201 - t211;
t224 = pkin(2) * t268 - t431;
t202 = t224 * t403;
t376 = t268 * t403;
t411 = qJ(4) * t376 + t202;
t343 = rSges(4,3) * t270 + t450;
t410 = -t224 + t343;
t397 = qJD(3) * t268;
t409 = qJ(3) * t370 + t271 * t397;
t408 = rSges(3,2) * t376 + rSges(3,3) * t402;
t407 = t271 * pkin(1) + t269 * pkin(6);
t171 = Icges(3,3) * t269 + t271 * t322;
t405 = qJD(1) * t171;
t180 = Icges(4,1) * t269 - t271 * t325;
t404 = qJD(1) * t180;
t396 = qJD(4) * t270;
t394 = t471 + t469;
t393 = -rSges(5,3) - pkin(2) - qJ(4);
t96 = Icges(6,4) * t169 + Icges(6,2) * t168 + Icges(6,6) * t420;
t98 = Icges(6,1) * t169 + Icges(6,4) * t168 + Icges(6,5) * t420;
t340 = t251 * t98 + t252 * t96;
t42 = Icges(6,5) * t87 + Icges(6,6) * t86 + Icges(6,3) * t290;
t44 = Icges(6,4) * t87 + Icges(6,2) * t86 + Icges(6,6) * t290;
t46 = Icges(6,1) * t87 + Icges(6,4) * t86 + Icges(6,5) * t290;
t94 = Icges(6,5) * t169 + Icges(6,6) * t168 + Icges(6,3) * t420;
t11 = (qJD(2) * t340 + t42) * t268 + (qJD(2) * t94 - t251 * t46 - t252 * t44 + (t251 * t96 - t252 * t98) * qJD(5)) * t270;
t16 = t101 * t420 + t168 * t102 + t169 * t103 + t154 * t290 + t86 * t155 + t87 * t156;
t388 = t11 / 0.2e1 + t16 / 0.2e1;
t95 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t419;
t97 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t419;
t341 = t251 * t97 + t252 * t95;
t43 = Icges(6,5) * t89 + Icges(6,6) * t88 - Icges(6,3) * t492;
t45 = Icges(6,4) * t89 + Icges(6,2) * t88 - Icges(6,6) * t492;
t47 = Icges(6,1) * t89 + Icges(6,4) * t88 - Icges(6,5) * t492;
t93 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t419;
t10 = (qJD(2) * t341 + t43) * t268 + (qJD(2) * t93 - t251 * t47 - t252 * t45 + (t251 * t95 - t252 * t97) * qJD(5)) * t270;
t17 = t101 * t419 + t166 * t102 + t167 * t103 - t154 * t492 + t88 * t155 + t89 * t156;
t387 = t17 / 0.2e1 + t10 / 0.2e1;
t34 = t268 * t94 - t270 * t340;
t51 = t154 * t420 + t155 * t168 + t156 * t169;
t386 = t34 / 0.2e1 + t51 / 0.2e1;
t33 = t268 * t93 - t270 * t341;
t50 = t154 * t419 + t166 * t155 + t167 * t156;
t385 = -t50 / 0.2e1 - t33 / 0.2e1;
t109 = Icges(5,5) * t195 + Icges(5,6) * t194 + Icges(5,3) * t419;
t384 = t109 * t420;
t383 = t109 * t419;
t110 = Icges(5,5) * t197 + Icges(5,6) * t196 + Icges(5,3) * t420;
t382 = t110 * t420;
t381 = t110 * t419;
t235 = pkin(2) * t373;
t378 = t269 * (pkin(2) * t374 + t269 * t397 - t235 + (t268 * t402 + t371) * qJ(3)) + t271 * (-pkin(2) * t492 - qJ(3) * t376 + t409) + t200 * t402;
t115 = t195 * rSges(5,1) + t194 * rSges(5,2) + rSges(5,3) * t419;
t249 = pkin(6) * t402;
t377 = t249 + t409;
t368 = -t325 * qJD(2) / 0.2e1 + t322 * t489;
t367 = -t401 / 0.2e1;
t366 = -qJ(3) - t455;
t150 = t410 * t271;
t365 = -qJ(4) * t268 - t224;
t212 = qJ(4) * t420 - t258;
t362 = t271 * t211 + t269 * t212 + t414;
t361 = pkin(4) * t357 + t246 * t402 + t267 * t492;
t237 = t271 * t396;
t360 = t237 + t377;
t359 = rSges(4,1) * t402 + rSges(4,2) * t492 + rSges(4,3) * t370;
t358 = t407 + t201;
t347 = rSges(5,1) * t265 + rSges(5,2) * t266;
t165 = -t270 * t347 + t448;
t355 = -t165 + t365;
t354 = t268 * t447 - pkin(1);
t353 = t109 / 0.2e1 + t175 / 0.2e1 - t178 / 0.2e1;
t352 = t110 / 0.2e1 + t174 / 0.2e1 + t479 / 0.2e1;
t350 = rSges(3,1) * t270 - rSges(3,2) * t268;
t349 = -t142 * rSges(5,1) - t141 * rSges(5,2);
t348 = -rSges(5,1) * t197 - rSges(5,2) * t196;
t27 = t166 * t95 + t167 * t97 + t419 * t93;
t28 = t166 * t96 + t167 * t98 + t419 * t94;
t339 = t269 * t28 + t27 * t271;
t18 = t27 * t269 - t271 * t28;
t29 = t168 * t95 + t169 * t97 + t420 * t93;
t30 = t168 * t96 + t169 * t98 + t420 * t94;
t338 = t269 * t30 + t271 * t29;
t19 = t29 * t269 - t271 * t30;
t337 = t269 * t34 + t271 * t33;
t336 = t33 * t269 - t271 * t34;
t257 = t271 * pkin(6);
t300 = t268 * t366 - pkin(1);
t275 = (-pkin(2) - t446) * t270 + t300;
t55 = t269 * t275 + t257 + t346 + t428;
t56 = t358 + t488;
t335 = t269 * t56 + t271 * t55;
t285 = t270 * t393 - pkin(1) - t432;
t278 = t285 * t269;
t70 = t257 + t258 + t278 + t348;
t71 = t358 + t115 + t211;
t333 = t269 * t71 + t271 * t70;
t321 = Icges(5,5) * t265 + Icges(5,6) * t266;
t305 = -t396 - t397;
t184 = rSges(3,1) * t419 + t482;
t185 = rSges(4,3) * t425 + t483;
t304 = -pkin(1) - t350;
t302 = t365 - t415;
t228 = qJ(4) * t373;
t250 = pkin(3) * t402;
t301 = t269 * (qJD(1) * t211 + t269 * t396 - t228) + t271 * (-qJ(4) * t492 + t237 + t250) + t212 * t402 + t378;
t108 = t355 * t271;
t298 = qJD(2) * t226;
t295 = qJD(2) * (Icges(4,4) * t268 + Icges(4,5) * t270);
t294 = qJD(2) * (-Icges(3,5) * t268 - Icges(3,6) * t270);
t81 = t302 * t271;
t288 = -qJ(4) * t399 - qJD(4) * t268 - t192;
t5 = (qJ(4) * t372 - t250 + t361 + t49) * t271 + (t228 + t48 + (t267 * t268 + t391) * t400) * t269 + (t417 * t271 + (-t379 + (-pkin(3) + t246) * t269 + t412 - t445) * t269) * qJD(1) + t301;
t80 = t302 * t269;
t286 = t398 * t81 + t400 * t80 - t5;
t107 = t355 * t269;
t116 = rSges(5,3) * t420 - t348;
t15 = t269 * (-rSges(5,3) * t373 - t349) + t271 * (-rSges(5,3) * t372 + t416) + (t271 * t116 + (-t115 + t412) * t269) * qJD(1) + t301;
t284 = t107 * t400 + t108 * t398 - t15;
t283 = (rSges(4,2) - pkin(2)) * t270 + t354;
t282 = -(rSges(5,3) * t270 + t268 * t347) * qJD(2) + t288;
t111 = Icges(5,4) * t195 + Icges(5,2) * t194 + Icges(5,6) * t419;
t113 = Icges(5,1) * t195 + Icges(5,4) * t194 + Icges(5,5) * t419;
t281 = t111 * t461 + t113 * t463 - t173 / 0.2e1 + t176 / 0.2e1;
t112 = Icges(5,4) * t197 + Icges(5,2) * t196 + Icges(5,6) * t420;
t114 = Icges(5,1) * t197 + Icges(5,4) * t196 + Icges(5,5) * t420;
t280 = t112 * t462 + t114 * t464 + t172 / 0.2e1 + t480 / 0.2e1;
t104 = (-rSges(6,1) * t252 + rSges(6,2) * t251) * t395 + (rSges(6,3) * t270 + t268 * t345) * qJD(2);
t279 = -t291 * qJD(2) - t104 + t288;
t21 = (t158 * t398 + t49) * t268 + (qJD(2) * t99 - t104 * t271 + t158 * t403) * t270;
t22 = (-t158 * t400 - t48) * t268 + (-qJD(2) * t100 + t269 * t104 + t158 * t402) * t270;
t52 = t332 * t270;
t274 = qJD(2) * t52 + t21 * t269 + t22 * t271 + (-t269 * t61 + t271 * t62) * qJD(1);
t25 = t458 * t372 + (t270 * t458 + t300) * t403 + t360 + t361 + t452;
t26 = t235 + ((t268 * t446 + t270 * t366) * qJD(2) + t305) * t269 + ((-pkin(6) - t246) * t269 + t275 * t271) * qJD(1) - t356;
t39 = qJD(1) * t278 + t372 * t393 + t250 + t360 + t416;
t40 = t228 + t235 + ((-t431 + t448) * qJD(2) + t305) * t269 + ((-pkin(3) - pkin(6)) * t269 + t285 * t271) * qJD(1) + t349;
t273 = (t25 * t269 + t26 * t271 + t402 * t56 - t403 * t55) * t469 + (t269 * t39 + t271 * t40 + t402 * t71 - t403 * t70) * t471;
t24 = t269 * t417 + t271 * t445 + t362;
t31 = t271 * t279 + t403 * t415 + t411;
t32 = qJD(1) * t81 + t269 * t279;
t41 = t115 * t271 + t269 * t116 + t362;
t53 = t165 * t403 + t271 * t282 + t411;
t54 = qJD(1) * t108 + t269 * t282;
t272 = (qJD(2) * t24 + t269 * t32 + t271 * t31 + t402 * t80 - t403 * t81) * t469 + (qJD(2) * t41 + t107 * t402 - t108 * t403 + t269 * t54 + t271 * t53) * t471;
t210 = t350 * qJD(2);
t186 = t269 * t344 - t451;
t183 = t269 * t350 - t449;
t153 = (Icges(5,5) * t270 + t268 * t329) * qJD(2);
t152 = (Icges(5,6) * t270 + t268 * t324) * qJD(2);
t149 = t410 * t269;
t148 = t184 + t407;
t147 = t269 * t304 + t257 + t449;
t146 = t490 * t268 * t399;
t130 = qJD(1) * t478 + t271 * t295;
t129 = t269 * t295 + t404;
t120 = t269 * t294 + t405;
t119 = -qJD(1) * t170 + t271 * t294;
t118 = t185 + t358;
t117 = t269 * t283 + t257 + t451;
t92 = t226 * t400 + ((-rSges(3,3) - pkin(6)) * t269 + t304 * t271) * qJD(1);
t91 = -rSges(3,1) * t492 - rSges(3,2) * t370 - pkin(1) * t403 + t249 + t408;
t83 = qJD(1) * t150 + t269 * t413;
t82 = t271 * t413 - t343 * t403 + t202;
t79 = t269 * t308 + t271 * t478;
t78 = -t180 * t271 + t487;
t77 = -t269 * t478 + t485;
t76 = t269 * t180 + t271 * t309;
t75 = t269 * t171 - t271 * t310;
t74 = t269 * t170 - t484;
t73 = -t171 * t271 - t486;
t72 = -t170 * t271 - t269 * t311;
t69 = Icges(5,1) * t144 + Icges(5,4) * t143 - Icges(5,5) * t492;
t68 = Icges(5,1) * t142 + Icges(5,4) * t141 + Icges(5,5) * t290;
t67 = Icges(5,4) * t144 + Icges(5,2) * t143 - Icges(5,6) * t492;
t66 = Icges(5,4) * t142 + Icges(5,2) * t141 + Icges(5,6) * t290;
t63 = t185 * t271 + t269 * t186 + t414;
t60 = t235 + (-t397 + (t270 * t447 - t450) * qJD(2)) * t269 + ((-rSges(4,1) - pkin(6)) * t269 + t283 * t271) * qJD(1);
t59 = -pkin(2) * t372 + (t354 - t456) * t403 + t359 + t377;
t38 = t112 * t196 + t114 * t197 + t382;
t37 = t111 * t196 + t113 * t197 + t384;
t36 = t194 * t112 + t195 * t114 + t381;
t35 = t194 * t111 + t195 * t113 + t383;
t23 = (qJD(1) * t186 + t359) * t271 + (t343 * t400 + (-t185 - t201 + t483) * qJD(1)) * t269 + t378;
t13 = t51 * t268 + t270 * t338;
t12 = t50 * t268 + t270 * t339;
t9 = -t94 * t372 + t166 * t44 + t167 * t46 + t88 * t96 + t89 * t98 + (t271 * t42 - t403 * t94) * t270;
t8 = -t93 * t372 + t166 * t45 + t167 * t47 + t88 * t95 + t89 * t97 + (t271 * t43 - t403 * t93) * t270;
t7 = t94 * t374 + t168 * t44 + t169 * t46 + t86 * t96 + t87 * t98 + (t270 * t42 - t401 * t94) * t269;
t6 = t93 * t374 + t168 * t45 + t169 * t47 + t86 * t95 + t87 * t97 + (t270 * t43 - t401 * t93) * t269;
t4 = qJD(1) * t339 + t8 * t269 - t271 * t9;
t3 = qJD(1) * t338 + t6 * t269 - t271 * t7;
t2 = (-qJD(2) * t339 + t17) * t268 + (-qJD(1) * t18 + qJD(2) * t50 + t269 * t9 + t271 * t8) * t270;
t1 = (-qJD(2) * t338 + t16) * t268 + (-qJD(1) * t19 + qJD(2) * t51 + t269 * t7 + t271 * t6) * t270;
t20 = [(t147 * t92 + t148 * t91) * t477 + (t117 * t60 + t118 * t59) * t476 + (t39 * t71 + t40 * t70) * t475 + (t25 * t56 + t26 * t55) * t474 + t351 - t252 * t156 * t395 - t153 * t426 + (-t102 * t252 - t152 * t266 - t427) * t270 + (Icges(5,3) * t268 - t270 * t321 + t317 + t327 + t501) * t399 + (Icges(5,3) * t270 + t163 * t266 + t164 * t265 + t268 * t321 + t319 + t331 - t502) * t401; (t163 * t497 + t164 * t496 - t196 * t152 / 0.2e1 - t197 * t153 / 0.2e1 + t368 * t271 - t388) * t271 + (t143 * t468 + t144 * t467 + t152 * t466 + t153 * t465 + t269 * t368 + t387) * t269 + m(3) * ((-t269 * t91 - t271 * t92) * t226 + (-t147 * t271 - t148 * t269) * t210) + m(4) * (t117 * t82 + t118 * t83 + t149 * t59 + t150 * t60) + m(5) * (t107 * t39 + t108 * t40 + t53 * t70 + t54 * t71) + m(6) * (t25 * t80 + t26 * t81 + t31 * t55 + t32 * t56) + ((Icges(5,5) * t496 + Icges(5,6) * t497 + t175 * t493 + t178 * t454 + t290 * t500 + t498 * t495) * t271 + (Icges(5,5) * t144 / 0.2e1 + Icges(5,6) * t143 / 0.2e1 + t492 * t500 + t498 * t494 + (t174 + t479) * t493) * t269) * t268 + ((t173 * t493 + t176 * t454 + t66 * t461 + t68 * t463 + t499 * t495) * t271 + (t67 * t462 + t69 * t464 + t499 * t494 + (t172 + t480) * t493) * t269) * t270 + ((t269 * t353 - t271 * t352) * t270 + (t269 * t281 + t271 * t280) * t268) * qJD(2) + ((-t148 * t457 + t163 * t466 + t164 * t465 + t268 * t353 - t270 * t281 - t385) * t271 + (t147 * t457 + t196 * t468 + t197 * t467 + t268 * t352 + t270 * t280 + t386) * t269) * qJD(1); (t24 * t5 + t31 * t81 + t32 * t80) * t474 - t271 * t3 + t269 * t4 + (t107 * t54 + t108 * t53 + t15 * t41) * t475 + (t149 * t83 + t150 * t82 + t23 * t63) * t476 + ((t269 * t183 + t184 * t271) * ((qJD(1) * t183 - t271 * t298 + t408) * t271 + (-t269 * t298 + (-t184 + t482) * qJD(1)) * t269) + t406 * t226 * t210) * t477 - t271 * ((-t141 * t112 - t142 * t114 - t196 * t66 - t197 * t68 + (t37 - t381) * qJD(1)) * t271 + (t141 * t111 + t142 * t113 + t196 * t67 + t197 * t69 + (t38 + t383) * qJD(1)) * t269) + t269 * ((t269 * t119 + (t74 + t486) * qJD(1)) * t269 + (t75 * qJD(1) + (t172 * t399 + t174 * t401) * t271 + (-t120 + (-t173 * t270 - t175 * t268) * qJD(2) + (t171 - t311) * qJD(1)) * t269) * t271) + t269 * ((t269 * t130 + (t77 - t487) * qJD(1)) * t269 + (t76 * qJD(1) + (t399 * t480 + t401 * t479) * t271 + (-t129 + (t176 * t270 + t178 * t268) * qJD(2) + (t180 + t308) * qJD(1)) * t269) * t271) + t269 * ((t143 * t111 + t144 * t113 + t194 * t67 + t195 * t69 + (t36 - t384) * qJD(1)) * t269 + (-t143 * t112 - t144 * t114 - t194 * t66 - t195 * t68 + (t35 + t382) * qJD(1)) * t271) - t271 * ((t129 * t271 + (t78 - t485) * qJD(1)) * t271 + (t79 * qJD(1) + (t176 * t399 + t178 * t401 + t404) * t269 + (-t130 + (t268 * t479 + t270 * t480) * qJD(2) + t309 * qJD(1)) * t271) * t269) - t271 * ((t120 * t271 + (t73 + t484) * qJD(1)) * t271 + (t72 * qJD(1) + (-t173 * t399 - t175 * t401 + t405) * t269 + (-t119 + (t172 * t270 + t174 * t268) * qJD(2) - t310 * qJD(1)) * t271) * t269) + (t19 + (-t38 - t72 - t79) * t271 + (t37 + t73 + t78) * t269) * t403 + (t18 + (-t36 - t74 - t77) * t271 + (t35 + t75 + t76) * t269) * t402; 0.2e1 * (t335 * t469 + (t117 * t271 + t118 * t269) * t473 + t333 * t471) * t399 + 0.2e1 * ((-t117 * t403 + t118 * t402 + t269 * t59 + t271 * t60) * t473 + t273) * t268; 0.2e1 * (t286 * t469 + t284 * t471 + (t149 * t400 + t150 * t398 - t23) * t473) * t270 + 0.2e1 * ((qJD(2) * t63 + t149 * t402 - t150 * t403 + t269 * t83 + t271 * t82) * t473 + t272) * t268; 0.4e1 * (t473 + t394) * t146; 0.2e1 * (t333 * t472 + t335 * t470) * t401 + 0.2e1 * t273 * t270; 0.2e1 * (t284 * t472 + t286 * t470) * t268 + 0.2e1 * t272 * t270; 0.2e1 * t394 * (-t268 ^ 2 + t270 ^ 2) * t490 * qJD(2); -0.4e1 * t394 * t146; m(6) * (t21 * t56 + t22 * t55 + t25 * t62 + t26 * t61) + (-t269 * t386 + t271 * t385) * t401 + (t387 * t271 + t388 * t269 + (t269 * t385 + t271 * t386) * qJD(1)) * t270 + t453; m(6) * (t14 * t24 + t21 * t80 + t22 * t81 + t31 * t61 + t32 * t62 + t5 * t52) + (-t1 / 0.2e1 + t12 * t454 + t18 * t367 + (t33 * qJD(1) - t11) * t460) * t271 + (t19 * t367 + t13 * t454 + t2 / 0.2e1 + (t34 * qJD(1) + t10) * t460) * t269 + (t3 * t495 + t4 * t459 + t336 * t489 + (t19 * t459 - t269 * t18 / 0.2e1) * qJD(1)) * t270; m(6) * (t274 * t268 + t270 * t481); m(6) * (-t268 * t481 + t274 * t270); (t14 * t52 + t21 * t62 + t22 * t61) * t474 + ((-t271 * t12 - t269 * t13 - t268 * t337) * qJD(2) + t453) * t268 + (t271 * t2 + t269 * t1 + t268 * (t10 * t271 + t11 * t269) + (t58 * t268 + t270 * t337) * qJD(2) + (-t269 * t12 + t271 * t13 - t268 * t336) * qJD(1)) * t270;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t20(1), t20(2), t20(4), t20(7), t20(11); t20(2), t20(3), t20(5), t20(8), t20(12); t20(4), t20(5), t20(6), t20(9), t20(13); t20(7), t20(8), t20(9), t20(10), t20(14); t20(11), t20(12), t20(13), t20(14), t20(15);];
Mq = res;
