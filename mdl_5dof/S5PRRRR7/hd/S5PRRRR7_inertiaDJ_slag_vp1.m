% Calculate time derivative of joint inertia matrix for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR7_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR7_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:49
% EndTime: 2019-12-05 17:12:11
% DurationCPUTime: 10.32s
% Computational Cost: add. (38512->747), mult. (43914->1119), div. (0->0), fcn. (41906->10), ass. (0->397)
t323 = cos(pkin(9));
t492 = t323 ^ 2;
t322 = sin(pkin(9));
t493 = t322 ^ 2;
t500 = t492 + t493;
t325 = sin(qJ(2));
t327 = cos(qJ(2));
t499 = qJD(2) * (rSges(3,1) * t325 + rSges(3,2) * t327);
t321 = qJ(3) + qJ(4);
t316 = cos(t321);
t431 = pkin(4) * t316;
t251 = -pkin(8) * t327 + t325 * t431;
t326 = cos(qJ(3));
t485 = t326 * pkin(3);
t268 = -pkin(7) * t327 + t325 * t485;
t496 = 2 * m(4);
t495 = 2 * m(5);
t494 = 2 * m(6);
t491 = t322 / 0.2e1;
t490 = -t323 / 0.2e1;
t489 = t323 / 0.2e1;
t488 = -t327 / 0.2e1;
t317 = qJ(5) + t321;
t311 = sin(t317);
t312 = cos(t317);
t472 = Icges(6,4) * t312;
t372 = -Icges(6,2) * t311 + t472;
t265 = -Icges(6,6) * t327 + t325 * t372;
t473 = Icges(6,4) * t311;
t376 = Icges(6,1) * t312 - t473;
t266 = -Icges(6,5) * t327 + t325 * t376;
t359 = t265 * t311 - t266 * t312;
t369 = Icges(6,5) * t312 - Icges(6,6) * t311;
t264 = -Icges(6,3) * t327 + t325 * t369;
t469 = t264 * t327;
t137 = -t325 * t359 - t469;
t319 = qJD(3) + qJD(4);
t314 = qJD(5) + t319;
t465 = t314 * t325;
t187 = (-Icges(6,2) * t312 - t473) * t465 + (Icges(6,6) * t325 + t327 * t372) * qJD(2);
t188 = (-Icges(6,1) * t311 - t472) * t465 + (Icges(6,5) * t325 + t327 * t376) * qJD(2);
t429 = qJD(2) * t325;
t409 = t322 * t429;
t344 = t314 * t323 + t409;
t459 = t322 * t327;
t422 = t314 * t459;
t224 = t311 * t344 - t312 * t422;
t225 = -t311 * t422 - t312 * t344;
t428 = qJD(2) * t327;
t408 = t322 * t428;
t128 = Icges(6,5) * t225 + Icges(6,6) * t224 + Icges(6,3) * t408;
t130 = Icges(6,4) * t225 + Icges(6,2) * t224 + Icges(6,6) * t408;
t132 = Icges(6,1) * t225 + Icges(6,4) * t224 + Icges(6,5) * t408;
t273 = -t311 * t459 - t312 * t323;
t274 = -t311 * t323 + t312 * t459;
t460 = t322 * t325;
t196 = Icges(6,5) * t274 + Icges(6,6) * t273 + Icges(6,3) * t460;
t198 = Icges(6,4) * t274 + Icges(6,2) * t273 + Icges(6,6) * t460;
t200 = Icges(6,1) * t274 + Icges(6,4) * t273 + Icges(6,5) * t460;
t366 = -t198 * t311 + t200 * t312;
t33 = (qJD(2) * t366 - t128) * t327 + (qJD(2) * t196 + (-t198 * t314 + t132) * t312 + (-t200 * t314 - t130) * t311) * t325;
t407 = t323 * t429;
t345 = -t314 * t322 + t407;
t456 = t323 * t327;
t421 = t314 * t456;
t226 = t311 * t345 - t312 * t421;
t227 = -t311 * t421 - t312 * t345;
t406 = t323 * t428;
t129 = Icges(6,5) * t227 + Icges(6,6) * t226 + Icges(6,3) * t406;
t131 = Icges(6,4) * t227 + Icges(6,2) * t226 + Icges(6,6) * t406;
t133 = Icges(6,1) * t227 + Icges(6,4) * t226 + Icges(6,5) * t406;
t275 = -t311 * t456 + t312 * t322;
t276 = t311 * t322 + t312 * t456;
t457 = t323 * t325;
t197 = Icges(6,5) * t276 + Icges(6,6) * t275 + Icges(6,3) * t457;
t199 = Icges(6,4) * t276 + Icges(6,2) * t275 + Icges(6,6) * t457;
t201 = Icges(6,1) * t276 + Icges(6,4) * t275 + Icges(6,5) * t457;
t365 = -t199 * t311 + t201 * t312;
t34 = (qJD(2) * t365 - t129) * t327 + (qJD(2) * t197 + (-t199 * t314 + t133) * t312 + (-t201 * t314 - t131) * t311) * t325;
t100 = -t197 * t327 + t325 * t365;
t99 = -t196 * t327 + t325 * t366;
t380 = t100 * t323 + t99 * t322;
t453 = t327 * ((-Icges(6,5) * t311 - Icges(6,6) * t312) * t465 + (Icges(6,3) * t325 + t327 * t369) * qJD(2));
t7 = (t453 + (t327 * t359 + t380) * qJD(2)) * t327 + (t34 * t323 + t33 * t322 - (qJD(2) * t264 + (-t265 * t314 + t188) * t312 + (-t266 * t314 - t187) * t311) * t327 + t137 * qJD(2)) * t325;
t486 = t327 * t7;
t85 = t196 * t457 + t198 * t275 + t200 * t276;
t483 = t322 * t85;
t315 = sin(t321);
t283 = -t315 * t459 - t316 * t323;
t284 = -t315 * t323 + t316 * t459;
t210 = Icges(5,5) * t284 + Icges(5,6) * t283 + Icges(5,3) * t460;
t212 = Icges(5,4) * t284 + Icges(5,2) * t283 + Icges(5,6) * t460;
t214 = Icges(5,1) * t284 + Icges(5,4) * t283 + Icges(5,5) * t460;
t285 = -t315 * t456 + t316 * t322;
t286 = t315 * t322 + t316 * t456;
t92 = t210 * t457 + t212 * t285 + t214 * t286;
t482 = t322 * t92;
t84 = t197 * t460 + t199 * t273 + t201 * t274;
t481 = t323 * t84;
t211 = Icges(5,5) * t286 + Icges(5,6) * t285 + Icges(5,3) * t457;
t213 = Icges(5,4) * t286 + Icges(5,2) * t285 + Icges(5,6) * t457;
t215 = Icges(5,1) * t286 + Icges(5,4) * t285 + Icges(5,5) * t457;
t91 = t211 * t460 + t213 * t283 + t215 * t284;
t480 = t323 * t91;
t324 = sin(qJ(3));
t477 = Icges(4,4) * t324;
t476 = Icges(4,4) * t326;
t475 = Icges(5,4) * t315;
t474 = Icges(5,4) * t316;
t455 = t324 * t327;
t301 = t322 * t326 - t323 * t455;
t454 = t326 * t327;
t461 = t322 * t324;
t302 = t323 * t454 + t461;
t232 = Icges(4,5) * t302 + Icges(4,6) * t301 + Icges(4,3) * t457;
t234 = Icges(4,4) * t302 + Icges(4,2) * t301 + Icges(4,6) * t457;
t236 = Icges(4,1) * t302 + Icges(4,4) * t301 + Icges(4,5) * t457;
t299 = -t322 * t455 - t323 * t326;
t458 = t323 * t324;
t300 = t322 * t454 - t458;
t104 = t232 * t460 + t234 * t299 + t236 * t300;
t471 = t104 * t323;
t231 = Icges(4,5) * t300 + Icges(4,6) * t299 + Icges(4,3) * t460;
t233 = Icges(4,4) * t300 + Icges(4,2) * t299 + Icges(4,6) * t460;
t235 = Icges(4,1) * t300 + Icges(4,4) * t299 + Icges(4,5) * t460;
t105 = t231 * t457 + t233 * t301 + t235 * t302;
t470 = t105 * t322;
t370 = Icges(5,5) * t316 - Icges(5,6) * t315;
t269 = -Icges(5,3) * t327 + t325 * t370;
t468 = t269 * t327;
t371 = Icges(4,5) * t326 - Icges(4,6) * t324;
t287 = -Icges(4,3) * t327 + t325 * t371;
t467 = t287 * t327;
t464 = t315 * t319;
t463 = t316 * t319;
t462 = t319 * t325;
t452 = t327 * ((-Icges(5,5) * t315 - Icges(5,6) * t316) * t462 + (Icges(5,3) * t325 + t327 * t370) * qJD(2));
t426 = qJD(3) * t325;
t451 = t327 * ((-Icges(4,5) * t324 - Icges(4,6) * t326) * t426 + (Icges(4,3) * t325 + t327 * t371) * qJD(2));
t135 = rSges(6,1) * t227 + rSges(6,2) * t226 + rSges(6,3) * t406;
t388 = pkin(4) * t464;
t329 = -t251 * qJD(2) - t388 * t327;
t387 = pkin(4) * t463;
t450 = t322 * t387 + t323 * t329 + t135;
t134 = rSges(6,1) * t225 + rSges(6,2) * t224 + rSges(6,3) * t408;
t202 = rSges(6,1) * t274 + rSges(6,2) * t273 + rSges(6,3) * t460;
t449 = t134 * t457 + t202 * t406;
t342 = t319 * t323 + t409;
t420 = t319 * t459;
t246 = t315 * t342 - t316 * t420;
t247 = -t315 * t420 - t316 * t342;
t150 = rSges(5,1) * t247 + rSges(5,2) * t246 + rSges(5,3) * t408;
t216 = rSges(5,1) * t284 + rSges(5,2) * t283 + rSges(5,3) * t460;
t448 = t150 * t457 + t216 * t406;
t343 = -t319 * t322 + t407;
t419 = t319 * t456;
t248 = t315 * t343 - t316 * t419;
t249 = -t315 * t419 - t316 * t343;
t151 = rSges(5,1) * t249 + rSges(5,2) * t248 + rSges(5,3) * t406;
t427 = qJD(3) * t324;
t424 = pkin(3) * t427;
t330 = -t268 * qJD(2) - t327 * t424;
t425 = qJD(3) * t326;
t423 = pkin(3) * t425;
t184 = t322 * t423 + t323 * t330;
t447 = -t151 - t184;
t333 = pkin(8) * t325 + t327 * t431;
t399 = pkin(4) * t315;
t169 = t322 * t333 - t323 * t399;
t182 = t202 * t457;
t446 = t169 * t457 + t182;
t170 = t322 * t399 + t323 * t333;
t203 = rSges(6,1) * t276 + rSges(6,2) * t275 + rSges(6,3) * t457;
t185 = t203 * t429;
t445 = t170 * t429 + t185;
t444 = t169 + t202;
t443 = t170 + t203;
t173 = qJD(2) * t333 - t325 * t388;
t381 = rSges(6,1) * t312 - rSges(6,2) * t311;
t191 = (-rSges(6,1) * t311 - rSges(6,2) * t312) * t465 + (rSges(6,3) * t325 + t327 * t381) * qJD(2);
t442 = -t173 - t191;
t183 = t322 * t330 - t323 * t423;
t335 = pkin(7) * t325 + t327 * t485;
t240 = -pkin(3) * t458 + t322 * t335;
t441 = t183 * t457 + t240 * t406;
t382 = rSges(5,1) * t316 - rSges(5,2) * t315;
t209 = (-rSges(5,1) * t315 - rSges(5,2) * t316) * t462 + (rSges(5,3) * t325 + t327 * t382) * qJD(2);
t250 = qJD(2) * t335 - t325 * t424;
t440 = -t209 - t250;
t217 = rSges(5,1) * t286 + rSges(5,2) * t285 + rSges(5,3) * t457;
t241 = pkin(3) * t461 + t323 * t335;
t439 = -t217 - t241;
t438 = t327 * t240 + t268 * t460;
t383 = rSges(4,1) * t326 - rSges(4,2) * t324;
t245 = (-rSges(4,1) * t324 - rSges(4,2) * t326) * t426 + (rSges(4,3) * t325 + t327 * t383) * qJD(2);
t385 = pkin(2) * t327 + pkin(6) * t325;
t307 = t385 * qJD(2);
t437 = -t245 - t307;
t267 = -rSges(6,3) * t327 + t325 * t381;
t436 = -t251 - t267;
t142 = t327 * t202 + t267 * t460;
t272 = -rSges(5,3) * t327 + t325 * t382;
t153 = t327 * t216 + t272 * t460;
t435 = -t268 - t272;
t310 = t325 * pkin(2) - t327 * pkin(6);
t434 = t500 * qJD(2) * t310;
t290 = -rSges(4,3) * t327 + t325 * t383;
t433 = -t290 - t310;
t432 = t500 * t385;
t418 = -t184 - t450;
t417 = t327 * t134 + t191 * t460 + t267 * t408;
t416 = t327 * t150 + t209 * t460 + t272 * t408;
t415 = -t241 - t443;
t414 = -t250 + t442;
t413 = t327 * t183 + t250 * t460 + t268 * t408;
t412 = -t307 + t440;
t411 = -t268 + t436;
t410 = -t310 + t435;
t405 = t324 * t429;
t404 = t326 * t429;
t403 = t460 / 0.2e1;
t402 = t457 / 0.2e1;
t401 = t429 / 0.2e1;
t400 = t428 / 0.2e1;
t398 = t443 * t327;
t397 = t439 * t327;
t119 = t322 * t329 - t323 * t387;
t396 = t119 * t457 + t169 * t406 + t449;
t395 = -t307 + t414;
t394 = t322 * t183 + t323 * t184 - t434;
t393 = t322 * t240 + t323 * t241 + t432;
t88 = t327 * t169 + t251 * t460 + t142;
t392 = -t310 + t411;
t112 = t264 * t460 + t265 * t273 + t266 * t274;
t113 = t264 * t457 + t265 * t275 + t266 * t276;
t351 = t128 * t325 + t196 * t428;
t35 = t130 * t273 + t132 * t274 + t198 * t224 + t200 * t225 + t322 * t351;
t350 = t129 * t325 + t197 * t428;
t36 = t131 * t273 + t133 * t274 + t199 * t224 + t201 * t225 + t322 * t350;
t83 = t196 * t460 + t198 * t273 + t200 * t274;
t5 = -(t187 * t273 + t188 * t274 + t224 * t265 + t225 * t266) * t327 + (t36 * t323 + (t35 - t453) * t322) * t325 + (t112 * t325 + (t481 + (t83 - t469) * t322) * t327) * qJD(2);
t37 = t130 * t275 + t132 * t276 + t198 * t226 + t200 * t227 + t323 * t351;
t38 = t131 * t275 + t133 * t276 + t199 * t226 + t201 * t227 + t323 * t350;
t86 = t197 * t457 + t199 * t275 + t201 * t276;
t6 = -(t187 * t275 + t188 * t276 + t226 * t265 + t227 * t266) * t327 + (t37 * t322 + (t38 - t453) * t323) * t325 + (t113 * t325 + (t483 + (t86 - t469) * t323) * t327) * qJD(2);
t391 = t6 * t457 + t5 * t460 + (-t112 * t327 + (t322 * t83 + t481) * t325) * t408 + (-t113 * t327 + (t323 * t86 + t483) * t325) * t406 + (-t137 * t327 + t325 * t380) * t429;
t390 = t322 * t400;
t389 = t323 * t400;
t386 = t415 * t327;
t378 = Icges(4,1) * t326 - t477;
t377 = Icges(5,1) * t316 - t475;
t374 = -Icges(4,2) * t324 + t476;
t373 = -Icges(5,2) * t315 + t474;
t364 = -t212 * t315 + t214 * t316;
t101 = -t210 * t327 + t325 * t364;
t363 = -t213 * t315 + t215 * t316;
t102 = -t211 * t327 + t325 * t363;
t368 = t101 * t322 + t102 * t323;
t362 = -t233 * t324 + t235 * t326;
t109 = -t231 * t327 + t325 * t362;
t361 = -t234 * t324 + t236 * t326;
t110 = -t232 * t327 + t325 * t361;
t367 = t109 * t322 + t110 * t323;
t238 = rSges(4,1) * t300 + rSges(4,2) * t299 + rSges(4,3) * t460;
t239 = rSges(4,1) * t302 + rSges(4,2) * t301 + rSges(4,3) * t457;
t360 = t238 * t323 - t239 * t322;
t270 = -Icges(5,6) * t327 + t325 * t373;
t271 = -Icges(5,5) * t327 + t325 * t377;
t358 = t270 * t315 - t271 * t316;
t288 = -Icges(4,6) * t327 + t325 * t374;
t289 = -Icges(4,5) * t327 + t325 * t378;
t356 = t288 * t324 - t289 * t326;
t353 = t327 * t119 + t173 * t460 + t251 * t408 + t417;
t144 = Icges(5,5) * t247 + Icges(5,6) * t246 + Icges(5,3) * t408;
t349 = t144 * t325 + t210 * t428;
t145 = Icges(5,5) * t249 + Icges(5,6) * t248 + Icges(5,3) * t406;
t348 = t145 * t325 + t211 * t428;
t260 = -qJD(3) * t300 + t322 * t405;
t261 = qJD(3) * t299 - t322 * t404;
t163 = Icges(4,5) * t261 + Icges(4,6) * t260 + Icges(4,3) * t408;
t347 = t163 * t325 + t231 * t428;
t262 = -qJD(3) * t302 + t323 * t405;
t263 = qJD(3) * t301 - t323 * t404;
t164 = Icges(4,5) * t263 + Icges(4,6) * t262 + Icges(4,3) * t406;
t346 = t164 * t325 + t232 * t428;
t341 = t391 - t486;
t20 = t322 * t36 - t323 * t35;
t21 = t322 * t38 - t323 * t37;
t340 = t5 * t490 + t6 * t491 + (t322 * t34 - t323 * t33) * t488 + t20 * t403 + t21 * t402 + (t322 * t84 - t323 * t83) * t390 + (t322 * t86 - t323 * t85) * t389 + (t100 * t322 - t323 * t99) * t401;
t336 = qJD(2) * (-Icges(3,5) * t325 - Icges(3,6) * t327);
t117 = t269 * t460 + t270 * t283 + t271 * t284;
t118 = t269 * t457 + t270 * t285 + t271 * t286;
t207 = (-Icges(5,2) * t316 - t475) * t462 + (Icges(5,6) * t325 + t327 * t373) * qJD(2);
t208 = (-Icges(5,1) * t315 - t474) * t462 + (Icges(5,5) * t325 + t327 * t377) * qJD(2);
t146 = Icges(5,4) * t247 + Icges(5,2) * t246 + Icges(5,6) * t408;
t148 = Icges(5,1) * t247 + Icges(5,4) * t246 + Icges(5,5) * t408;
t46 = t146 * t283 + t148 * t284 + t212 * t246 + t214 * t247 + t322 * t349;
t147 = Icges(5,4) * t249 + Icges(5,2) * t248 + Icges(5,6) * t406;
t149 = Icges(5,1) * t249 + Icges(5,4) * t248 + Icges(5,5) * t406;
t47 = t147 * t283 + t149 * t284 + t213 * t246 + t215 * t247 + t322 * t348;
t90 = t210 * t460 + t212 * t283 + t214 * t284;
t12 = -(t207 * t283 + t208 * t284 + t246 * t270 + t247 * t271) * t327 + (t47 * t323 + (t46 - t452) * t322) * t325 + (t117 * t325 + (t480 + (t90 - t468) * t322) * t327) * qJD(2);
t48 = t146 * t285 + t148 * t286 + t212 * t248 + t214 * t249 + t323 * t349;
t49 = t147 * t285 + t149 * t286 + t213 * t248 + t215 * t249 + t323 * t348;
t93 = t211 * t457 + t213 * t285 + t215 * t286;
t13 = -(t207 * t285 + t208 * t286 + t248 * t270 + t249 * t271) * t327 + (t48 * t322 + (t49 - t452) * t323) * t325 + (t118 * t325 + (t482 + (t93 - t468) * t323) * t327) * qJD(2);
t152 = -t325 * t358 - t468;
t334 = t12 * t460 + t13 * t457 + t391 + (-t117 * t327 + (t322 * t90 + t480) * t325) * t408 + (-t118 * t327 + (t323 * t93 + t482) * t325) * t406 + (-t152 * t327 + t325 * t368) * t429;
t42 = (qJD(2) * t364 - t144) * t327 + (qJD(2) * t210 + (-t212 * t319 + t148) * t316 + (-t214 * t319 - t146) * t315) * t325;
t43 = (qJD(2) * t363 - t145) * t327 + (qJD(2) * t211 + (-t213 * t319 + t149) * t316 + (-t215 * t319 - t147) * t315) * t325;
t14 = (t452 + (t327 * t358 + t368) * qJD(2)) * t327 + (t43 * t323 + t42 * t322 - (qJD(2) * t269 - t207 * t315 + t208 * t316 - t270 * t463 - t271 * t464) * t327 + t152 * qJD(2)) * t325;
t332 = (-t7 - t14) * t327 + t334;
t25 = t322 * t47 - t323 * t46;
t26 = t322 * t49 - t323 * t48;
t331 = t12 * t490 + t13 * t491 + (t322 * t43 - t323 * t42) * t488 + t25 * t403 + t26 * t402 + t340 + (t322 * t91 - t323 * t90) * t390 + (t322 * t93 - t323 * t92) * t389 + (-t101 * t323 + t102 * t322) * t401;
t294 = t323 * t336;
t293 = t322 * t336;
t253 = t433 * t323;
t252 = t433 * t322;
t244 = (-Icges(4,1) * t324 - t476) * t426 + (Icges(4,5) * t325 + t327 * t378) * qJD(2);
t243 = (-Icges(4,2) * t326 - t477) * t426 + (Icges(4,6) * t325 + t327 * t374) * qJD(2);
t229 = t500 * t499;
t220 = t241 * t429;
t219 = t240 * t457;
t204 = t217 * t429;
t195 = t216 * t457;
t190 = t437 * t323;
t189 = t437 * t322;
t175 = t410 * t323;
t174 = t410 * t322;
t172 = rSges(4,1) * t263 + rSges(4,2) * t262 + rSges(4,3) * t406;
t171 = rSges(4,1) * t261 + rSges(4,2) * t260 + rSges(4,3) * t408;
t168 = Icges(4,1) * t263 + Icges(4,4) * t262 + Icges(4,5) * t406;
t167 = Icges(4,1) * t261 + Icges(4,4) * t260 + Icges(4,5) * t408;
t166 = Icges(4,4) * t263 + Icges(4,2) * t262 + Icges(4,6) * t406;
t165 = Icges(4,4) * t261 + Icges(4,2) * t260 + Icges(4,6) * t408;
t162 = -t239 * t327 - t290 * t457;
t161 = t238 * t327 + t290 * t460;
t159 = -t325 * t356 - t467;
t154 = -t217 * t327 - t272 * t457;
t143 = -t203 * t327 - t267 * t457;
t139 = t287 * t457 + t288 * t301 + t289 * t302;
t138 = t287 * t460 + t288 * t299 + t289 * t300;
t136 = t360 * t325;
t125 = t392 * t323;
t124 = t392 * t322;
t123 = t412 * t323;
t122 = t412 * t322;
t121 = -t217 * t460 + t195;
t114 = -t203 * t460 + t182;
t111 = t238 * t322 + t239 * t323 + t432;
t108 = t435 * t457 + t397;
t107 = t153 + t438;
t106 = t232 * t457 + t234 * t301 + t236 * t302;
t103 = t231 * t460 + t233 * t299 + t235 * t300;
t98 = t171 * t322 + t172 * t323 - t434;
t97 = t395 * t323;
t96 = t395 * t322;
t95 = -t245 * t457 - t172 * t327 + (t239 * t325 - t290 * t456) * qJD(2);
t94 = t245 * t460 + t171 * t327 + (-t238 * t325 + t290 * t459) * qJD(2);
t89 = t436 * t457 - t398;
t87 = t439 * t460 + t195 + t219;
t82 = -t151 * t327 + t204 + (-t209 * t325 - t272 * t428) * t323;
t81 = -t216 * t429 + t416;
t80 = t216 * t322 + t217 * t323 + t393;
t79 = (t171 * t323 - t172 * t322) * t325 + t360 * t428;
t78 = -t135 * t327 + t185 + (-t191 * t325 - t267 * t428) * t323;
t77 = -t202 * t429 + t417;
t76 = -t443 * t460 + t446;
t75 = t411 * t457 + t386;
t74 = t88 + t438;
t73 = (-t151 * t325 - t217 * t428) * t322 + t448;
t72 = (-t135 * t325 - t203 * t428) * t322 + t449;
t71 = t150 * t322 + t151 * t323 + t394;
t69 = t415 * t460 + t219 + t446;
t67 = t322 * t444 + t323 * t443 + t393;
t64 = t204 + t220 + t447 * t327 + (t325 * t440 + t428 * t435) * t323;
t63 = (-t216 - t240) * t429 + t413 + t416;
t62 = t166 * t301 + t168 * t302 + t234 * t262 + t236 * t263 + t323 * t346;
t61 = t165 * t301 + t167 * t302 + t233 * t262 + t235 * t263 + t323 * t347;
t60 = t166 * t299 + t168 * t300 + t234 * t260 + t236 * t261 + t322 * t346;
t59 = t165 * t299 + t167 * t300 + t233 * t260 + t235 * t261 + t322 * t347;
t55 = (qJD(2) * t361 - t164) * t327 + (qJD(2) * t232 - t166 * t324 + t168 * t326 + (-t234 * t326 - t236 * t324) * qJD(3)) * t325;
t54 = (qJD(2) * t362 - t163) * t327 + (qJD(2) * t231 - t165 * t324 + t167 * t326 + (-t233 * t326 - t235 * t324) * qJD(3)) * t325;
t50 = (qJD(2) * t397 + t325 * t447) * t322 + t441 + t448;
t45 = -t450 * t327 + (t325 * t442 + t428 * t436) * t323 + t445;
t44 = -t429 * t444 + t353;
t39 = t450 * t323 + (t119 + t134) * t322 + t394;
t32 = (-qJD(2) * t398 - t325 * t450) * t322 + t396;
t31 = t220 + t418 * t327 + (t325 * t414 + t411 * t428) * t323 + t445;
t30 = (-t240 - t444) * t429 + t353 + t413;
t29 = t322 * t62 - t323 * t61;
t28 = t322 * t60 - t323 * t59;
t27 = (qJD(2) * t386 + t325 * t418) * t322 + t396 + t441;
t16 = -(t243 * t301 + t244 * t302 + t262 * t288 + t263 * t289) * t327 + (t61 * t322 + (t62 - t451) * t323) * t325 + (t139 * t325 + (t470 + (t106 - t467) * t323) * t327) * qJD(2);
t15 = -(t243 * t299 + t244 * t300 + t260 * t288 + t261 * t289) * t327 + (t60 * t323 + (t59 - t451) * t322) * t325 + (t138 * t325 + (t471 + (t103 - t467) * t322) * t327) * qJD(2);
t1 = [0; -m(3) * t229 + m(4) * t98 + m(5) * t71 + m(6) * t39; (t124 * t96 + t125 * t97 + t39 * t67) * t494 + (t122 * t174 + t123 * t175 + t71 * t80) * t495 + (t111 * t98 + t189 * t252 + t190 * t253) * t496 + 0.2e1 * m(3) * (-t229 + t499) * t500 * (rSges(3,1) * t327 - rSges(3,2) * t325) + (-t492 * t293 - t20 - t25 - t28) * t323 + (t493 * t294 + t21 + t26 + t29 + (-t322 * t293 + t323 * t294) * t323) * t322; m(4) * t79 + m(5) * t50 + m(6) * t27; t331 + (t322 * t55 - t323 * t54) * t488 + t15 * t490 + t16 * t491 + (t325 * (-t109 * t323 + t110 * t322) / 0.2e1 + ((-t103 * t323 + t322 * t104) * t491 + (-t323 * t105 + t106 * t322) * t489) * t327) * qJD(2) + (t28 * t491 + t29 * t489) * t325 + m(4) * (t111 * t79 + t136 * t98 + t161 * t190 + t162 * t189 + t252 * t95 + t253 * t94) + m(5) * (t107 * t123 + t108 * t122 + t174 * t64 + t175 * t63 + t50 * t80 + t71 * t87) + m(6) * (t124 * t31 + t125 * t30 + t27 * t67 + t39 * t69 + t74 * t97 + t75 * t96); t334 + t16 * t457 + t15 * t460 - t486 - t327 * t14 + (-t138 * t327 + (t103 * t322 + t471) * t325) * t408 + (-t139 * t327 + (t106 * t323 + t470) * t325) * t406 + (-t159 * t327 + t325 * t367) * t429 + (t27 * t69 + t30 * t74 + t31 * t75) * t494 + (t107 * t63 + t108 * t64 + t50 * t87) * t495 + (t136 * t79 + t161 * t94 + t162 * t95) * t496 - t327 * ((t451 + (t327 * t356 + t367) * qJD(2)) * t327 + (t55 * t323 + t54 * t322 - (qJD(2) * t287 - t243 * t324 + t244 * t326 - t288 * t425 - t289 * t427) * t327 + t159 * qJD(2)) * t325); m(5) * t73 + m(6) * t32; t331 + m(6) * (t124 * t45 + t125 * t44 + t32 * t67 + t39 * t76 + t88 * t97 + t89 * t96) + m(5) * (t121 * t71 + t122 * t154 + t123 * t153 + t174 * t82 + t175 * t81 + t73 * t80); t332 + m(6) * (t27 * t76 + t30 * t88 + t31 * t89 + t32 * t69 + t44 * t74 + t45 * t75) + m(5) * (t107 * t81 + t108 * t82 + t121 * t50 + t153 * t63 + t154 * t64 + t73 * t87); t332 + (t32 * t76 + t44 * t88 + t45 * t89) * t494 + (t121 * t73 + t153 * t81 + t154 * t82) * t495; m(6) * t72; m(6) * (t114 * t39 + t124 * t78 + t125 * t77 + t142 * t97 + t143 * t96 + t67 * t72) + t340; m(6) * (t114 * t27 + t142 * t30 + t143 * t31 + t69 * t72 + t74 * t77 + t75 * t78) + t341; m(6) * (t114 * t32 + t142 * t44 + t143 * t45 + t72 * t76 + t77 * t88 + t78 * t89) + t341; (t114 * t72 + t142 * t77 + t143 * t78) * t494 + t341;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
