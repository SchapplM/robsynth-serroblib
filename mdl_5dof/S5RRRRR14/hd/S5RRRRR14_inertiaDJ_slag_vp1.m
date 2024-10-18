% Calculate time derivative of joint inertia matrix for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR14_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR14_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR14_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:24
% EndTime: 2024-09-27 18:42:33
% DurationCPUTime: 6.74s
% Computational Cost: add. (68487->758), mult. (40209->1022), div. (0->0), fcn. (32113->22), ass. (0->389)
t380 = qJD(3) + qJD(4);
t369 = qJD(5) + t380;
t485 = t369 / 0.2e1;
t484 = t380 / 0.2e1;
t383 = qJ(3) + qJ(4);
t370 = pkin(5) + t383;
t358 = sin(t370);
t494 = t358 * t484;
t363 = qJ(5) + t370;
t353 = sin(t363);
t493 = t353 * t485;
t351 = t358 / 0.2e1;
t371 = pkin(5) - t383;
t359 = sin(t371);
t322 = t351 + t359 / 0.2e1;
t361 = cos(t371);
t352 = t361 / 0.2e1;
t360 = cos(t370);
t325 = t352 - t360 / 0.2e1;
t386 = cos(pkin(5));
t263 = Icges(5,5) * t325 + Icges(5,6) * t322 + Icges(5,3) * t386;
t264 = Icges(5,4) * t325 + Icges(5,2) * t322 + Icges(5,6) * t386;
t265 = Icges(5,1) * t325 + Icges(5,4) * t322 + Icges(5,5) * t386;
t324 = t352 + t360 / 0.2e1;
t372 = sin(t383);
t384 = qJ(1) + qJ(2);
t373 = sin(t384);
t375 = cos(t384);
t276 = t324 * t375 - t372 * t373;
t323 = t351 - t359 / 0.2e1;
t374 = cos(t383);
t277 = t323 * t375 + t373 * t374;
t385 = sin(pkin(5));
t468 = t375 * t385;
t113 = -t263 * t468 + t264 * t276 + t265 * t277;
t278 = -t324 * t373 - t372 * t375;
t279 = -t323 * t373 + t374 * t375;
t471 = t373 * t385;
t114 = t263 * t471 + t264 * t278 + t265 * t279;
t381 = qJD(1) + qJD(2);
t408 = -t324 * t381 - t374 * t380;
t337 = t359 * t484;
t422 = t372 * t381 - t337 + t494;
t197 = t373 * t422 + t375 * t408;
t473 = t372 * t380;
t409 = -t323 * t381 - t473;
t336 = t360 * t484;
t466 = t380 * t361;
t423 = t374 * t381 + t336 + t466 / 0.2e1;
t198 = -t373 * t423 + t375 * t409;
t465 = t381 * t385;
t438 = t375 * t465;
t123 = Icges(5,5) * t198 + Icges(5,6) * t197 + Icges(5,3) * t438;
t199 = t373 * t408 - t375 * t422;
t200 = t373 * t409 + t375 * t423;
t439 = t373 * t465;
t124 = Icges(5,5) * t200 + Icges(5,6) * t199 + Icges(5,3) * t439;
t125 = Icges(5,4) * t198 + Icges(5,2) * t197 + Icges(5,6) * t438;
t126 = Icges(5,4) * t200 + Icges(5,2) * t199 + Icges(5,6) * t439;
t127 = Icges(5,1) * t198 + Icges(5,4) * t197 + Icges(5,5) * t438;
t128 = Icges(5,1) * t200 + Icges(5,4) * t199 + Icges(5,5) * t439;
t204 = Icges(5,5) * t277 + Icges(5,6) * t276 - Icges(5,3) * t468;
t205 = Icges(5,5) * t279 + Icges(5,6) * t278 + Icges(5,3) * t471;
t206 = Icges(5,4) * t277 + Icges(5,2) * t276 - Icges(5,6) * t468;
t207 = Icges(5,4) * t279 + Icges(5,2) * t278 + Icges(5,6) * t471;
t208 = Icges(5,1) * t277 + Icges(5,4) * t276 - Icges(5,5) * t468;
t209 = Icges(5,1) * t279 + Icges(5,4) * t278 + Icges(5,5) * t471;
t290 = t336 - t466 / 0.2e1;
t292 = t337 + t494;
t22 = t124 * t386 + t126 * t322 + t128 * t325 + t206 * t290 + t208 * t292;
t23 = t123 * t386 + t125 * t322 + t127 * t325 + t207 * t290 + t209 * t292;
t242 = Icges(5,5) * t292 + Icges(5,6) * t290;
t243 = Icges(5,4) * t292 + Icges(5,2) * t290;
t244 = Icges(5,1) * t292 + Icges(5,4) * t290;
t469 = t375 * t381;
t34 = t197 * t264 + t198 * t265 + t243 * t278 + t244 * t279 + (t242 * t373 + t263 * t469) * t385;
t472 = t373 * t381;
t412 = t386 * t242 + t322 * t243 + t325 * t244 + t290 * t264 + t292 * t265;
t57 = t412 * t386;
t60 = -t204 * t468 + t206 * t276 + t208 * t277;
t61 = -t205 * t468 + t207 * t276 + t209 * t277;
t62 = t204 * t471 + t206 * t278 + t208 * t279;
t63 = t205 * t471 + t207 * t278 + t209 * t279;
t77 = t204 * t386 + t206 * t322 + t208 * t325;
t78 = t205 * t386 + t207 * t322 + t209 * t325;
t492 = (t113 * t386 + (t373 * t61 - t375 * t60) * t385) * t439 + (t114 * t386 + (t373 * t63 - t375 * t62) * t385) * t438 + (t34 * t386 + ((t278 * t125 + t279 * t127 + t197 * t207 + t198 * t209) * t373 + t63 * t469 - (t278 * t126 + t279 * t128 + t197 * t206 + t198 * t208) * t375 + t62 * t472 + ((t123 * t373 + t205 * t469) * t373 - (t124 * t373 + t204 * t469) * t375) * t385) * t385) * t471 + t386 * (t57 + ((t381 * t78 - t22) * t375 + (t381 * t77 + t23) * t373) * t385);
t390 = cos(qJ(4));
t367 = pkin(4) * t390 + pkin(3);
t387 = sin(qJ(4));
t388 = sin(qJ(3));
t391 = cos(qJ(3));
t447 = qJD(3) * t388;
t397 = (-t387 * t447 + (-t387 * t388 + t390 * t391) * qJD(4)) * pkin(4);
t446 = qJD(3) * t391;
t254 = (t367 * t446 + t397) * t386;
t441 = pkin(3) * t446;
t419 = t386 * t441;
t368 = t391 * pkin(3) + pkin(2);
t335 = pkin(4) * t374 + t368;
t450 = t335 - t368;
t491 = t450 * t381 + t254 - t419;
t482 = pkin(2) - t368;
t490 = t482 * t381 - t419;
t489 = 2 * m(3);
t488 = 2 * m(4);
t487 = 2 * m(5);
t486 = 2 * m(6);
t379 = t385 ^ 2;
t393 = pkin(8) + pkin(9);
t389 = sin(qJ(1));
t483 = pkin(1) * t389;
t362 = t375 * pkin(2);
t481 = -pkin(3) + t367;
t480 = pkin(1) * qJD(1);
t442 = pkin(3) * t447;
t321 = -pkin(4) * t473 - t442;
t296 = t375 * t321;
t382 = pkin(10) + t393;
t445 = pkin(4) * t387 * t391;
t280 = -t385 * t382 + (t367 * t388 + t445) * t386;
t463 = t386 * t388;
t328 = pkin(3) * t463 - t385 * t393;
t452 = t280 - t328;
t109 = t296 + (-t452 * t381 + t442) * t375 - t491 * t373;
t364 = -qJ(5) + t371;
t356 = cos(t364);
t350 = t356 / 0.2e1;
t355 = cos(t363);
t319 = t350 + t355 / 0.2e1;
t376 = qJ(5) + t383;
t366 = cos(t376);
t410 = -t319 * t381 - t366 * t369;
t354 = sin(t364);
t327 = t354 * t485;
t365 = sin(t376);
t424 = t365 * t381 - t327 + t493;
t169 = t373 * t424 + t375 * t410;
t349 = t353 / 0.2e1;
t318 = t349 - t354 / 0.2e1;
t411 = -t318 * t381 - t365 * t369;
t326 = t355 * t485;
t474 = t369 * t356;
t425 = t366 * t381 + t326 + t474 / 0.2e1;
t170 = -t373 * t425 + t375 * t411;
t107 = t170 * rSges(6,1) + t169 * rSges(6,2) + rSges(6,3) * t438;
t97 = t386 * t107;
t479 = t386 * t109 + t97;
t478 = Icges(4,4) * t388;
t477 = Icges(4,4) * t391;
t246 = rSges(5,1) * t292 + rSges(5,2) * t290;
t476 = t246 * t373;
t470 = t375 * t328;
t464 = t385 * t388;
t462 = t386 * t391;
t171 = t373 * t410 - t375 * t424;
t172 = t373 * t411 + t375 * t425;
t414 = -rSges(6,1) * t172 - rSges(6,2) * t171;
t108 = rSges(6,3) * t439 - t414;
t262 = t280 * t472;
t451 = t328 * t472 + t373 * t442;
t110 = t321 * t373 + t491 * t375 - t262 + t451;
t461 = -t108 - t110;
t269 = t319 * t375 - t365 * t373;
t270 = t318 * t375 + t366 * t373;
t193 = rSges(6,1) * t270 + rSges(6,2) * t269 - rSges(6,3) * t468;
t271 = -t319 * t373 - t365 * t375;
t272 = -t318 * t373 + t366 * t375;
t194 = t272 * rSges(6,1) + t271 * rSges(6,2) + rSges(6,3) * t471;
t117 = t193 * t471 + t194 * t468;
t180 = t386 * t194;
t316 = t375 * t335;
t338 = t375 * t368;
t184 = -t373 * t452 + t316 - t338;
t460 = t386 * t184 + t180;
t183 = t373 * t450 + t375 * t452;
t459 = -t183 - t193;
t458 = -t184 - t194;
t210 = rSges(5,1) * t277 + rSges(5,2) * t276 - rSges(5,3) * t468;
t211 = t279 * rSges(5,1) + t278 * rSges(5,2) + rSges(5,3) * t471;
t131 = t210 * t471 + t211 * t468;
t282 = t326 - t474 / 0.2e1;
t284 = t327 + t493;
t226 = rSges(6,1) * t284 + rSges(6,2) * t282;
t457 = -t226 - (t481 * t446 + t397) * t385;
t317 = t349 + t354 / 0.2e1;
t320 = t350 - t355 / 0.2e1;
t261 = rSges(6,1) * t320 + rSges(6,2) * t317 + rSges(6,3) * t386;
t232 = t261 * t439;
t260 = (t382 - t393) * t386 + (t481 * t388 + t445) * t385;
t456 = t260 * t439 + t232;
t348 = pkin(8) * t468;
t255 = -t482 * t373 + t348 + t470;
t421 = -t328 * t373 + t338;
t449 = pkin(8) * t471 + t362;
t256 = t421 - t449;
t455 = t255 * t471 + t256 * t468;
t454 = -t260 - t261;
t266 = rSges(5,1) * t325 + rSges(5,2) * t322 + rSges(5,3) * t386;
t309 = pkin(3) * t464 + (-pkin(8) + t393) * t386;
t453 = -t266 - t309;
t448 = qJD(3) * t385;
t444 = t107 * t468 + t108 * t471 + t193 * t438;
t392 = cos(qJ(1));
t443 = t392 * t480;
t440 = t389 * t480;
t129 = t198 * rSges(5,1) + t197 * rSges(5,2) + rSges(5,3) * t438;
t130 = rSges(5,1) * t200 + rSges(5,2) * t199 + rSges(5,3) * t439;
t437 = t129 * t468 + t130 * t471 + t210 * t438;
t334 = pkin(8) * t438;
t181 = -t334 + (-t328 * t381 - t442) * t375 + t490 * t373;
t182 = -pkin(8) * t439 - t490 * t375 - t451;
t436 = t181 * t468 + t182 * t471 + t255 * t438;
t305 = -t373 * t388 + t375 * t462;
t404 = t373 * t463 - t375 * t391;
t250 = qJD(3) * t404 - t305 * t381;
t306 = t373 * t391 + t375 * t463;
t307 = -t373 * t462 - t375 * t388;
t400 = t307 * qJD(3);
t251 = -t306 * t381 + t400;
t153 = t251 * rSges(4,1) + t250 * rSges(4,2) + rSges(4,3) * t438;
t435 = -t309 + t454;
t241 = -rSges(4,1) * t404 + t307 * rSges(4,2) + rSges(4,3) * t471;
t101 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t438;
t102 = Icges(6,5) * t172 + Icges(6,6) * t171 + Icges(6,3) * t439;
t103 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t438;
t104 = Icges(6,4) * t172 + Icges(6,2) * t171 + Icges(6,6) * t439;
t105 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t438;
t106 = Icges(6,1) * t172 + Icges(6,4) * t171 + Icges(6,5) * t439;
t189 = Icges(6,4) * t270 + Icges(6,2) * t269 - Icges(6,6) * t468;
t191 = Icges(6,1) * t270 + Icges(6,4) * t269 - Icges(6,5) * t468;
t13 = t102 * t386 + t104 * t317 + t106 * t320 + t189 * t282 + t191 * t284;
t190 = Icges(6,4) * t272 + Icges(6,2) * t271 + Icges(6,6) * t471;
t192 = Icges(6,1) * t272 + Icges(6,4) * t271 + Icges(6,5) * t471;
t14 = t101 * t386 + t103 * t317 + t105 * t320 + t190 * t282 + t192 * t284;
t187 = Icges(6,5) * t270 + Icges(6,6) * t269 - Icges(6,3) * t468;
t188 = Icges(6,5) * t272 + Icges(6,6) * t271 + Icges(6,3) * t471;
t223 = Icges(6,5) * t284 + Icges(6,6) * t282;
t224 = Icges(6,4) * t284 + Icges(6,2) * t282;
t225 = Icges(6,1) * t284 + Icges(6,4) * t282;
t257 = Icges(6,5) * t320 + Icges(6,6) * t317 + Icges(6,3) * t386;
t258 = Icges(6,4) * t320 + Icges(6,2) * t317 + Icges(6,6) * t386;
t259 = Icges(6,1) * t320 + Icges(6,4) * t317 + Icges(6,5) * t386;
t29 = t169 * t258 + t170 * t259 + t224 * t271 + t225 * t272 + (t223 * t373 + t257 * t469) * t385;
t413 = t386 * t223 + t317 * t224 + t320 * t225 + t282 * t258 + t284 * t259;
t43 = t413 * t386;
t53 = -t187 * t468 + t189 * t269 + t191 * t270;
t54 = -t188 * t468 + t190 * t269 + t192 * t270;
t55 = t187 * t471 + t189 * t271 + t191 * t272;
t56 = t188 * t471 + t190 * t271 + t192 * t272;
t68 = t187 * t386 + t189 * t317 + t191 * t320;
t69 = t188 * t386 + t190 * t317 + t192 * t320;
t89 = -t257 * t468 + t258 * t269 + t259 * t270;
t90 = t257 * t471 + t258 * t271 + t259 * t272;
t434 = (t29 * t386 + ((t271 * t103 + t272 * t105 + t169 * t190 + t170 * t192) * t373 + t56 * t469 - (t271 * t104 + t272 * t106 + t169 * t189 + t170 * t191) * t375 + t55 * t472 + ((t101 * t373 + t188 * t469) * t373 - (t102 * t373 + t187 * t469) * t375) * t385) * t385) * t471 + (t90 * t386 + (t373 * t56 - t375 * t55) * t385) * t438 + t386 * (t43 + ((t381 * t69 - t13) * t375 + (t381 * t68 + t14) * t373) * t385) + (t89 * t386 + (t373 * t54 - t375 * t53) * t385) * t439;
t433 = t471 / 0.2e1;
t432 = -t468 / 0.2e1;
t431 = t465 / 0.2e1;
t330 = t375 * rSges(3,1) - rSges(3,2) * t373;
t430 = t457 * t373;
t429 = t457 * t385;
t428 = t454 * t385;
t427 = t453 * t385;
t426 = -t335 * t381 - t254;
t420 = t379 * t441;
t44 = t183 * t471 + t184 * t468 + t117;
t418 = t373 * t431;
t417 = t375 * t431;
t303 = -rSges(3,1) * t469 + rSges(3,2) * t472;
t416 = t435 * t385;
t329 = -rSges(3,1) * t373 - rSges(3,2) * t375;
t252 = -qJD(3) * t306 + t307 * t381;
t253 = qJD(3) * t305 - t404 * t381;
t415 = -rSges(4,1) * t253 - rSges(4,2) * t252;
t217 = t241 + t449;
t407 = -t373 * t368 - t470;
t406 = t109 * t468 + t110 * t471 + t183 * t438 + t444;
t240 = rSges(4,1) * t306 + rSges(4,2) * t305 - rSges(4,3) * t468;
t30 = t171 * t258 + t172 * t259 + t224 * t269 + t225 * t270 + (-t223 * t375 + t257 * t472) * t385;
t2 = t30 * t386 + ((t269 * t103 + t270 * t105 + t171 * t190 + t172 * t192) * t373 + t54 * t469 - (t104 * t269 + t106 * t270 + t171 * t189 + t172 * t191) * t375 + t53 * t472 + ((-t101 * t375 + t188 * t472) * t373 - (-t102 * t375 + t187 * t472) * t375) * t385) * t385;
t403 = -t2 * t468 + t434;
t302 = t329 * t381;
t402 = t386 * t181 - t373 * t420;
t141 = -t280 * t373 + t194 + t316;
t160 = t421 + t211;
t138 = -pkin(2) * t472 + t153 + t334;
t401 = t43 + (t14 + t29) * t433 + (t13 + t30) * t432 + (t68 + t89) * t418 + (t69 + t90) * t417;
t300 = Icges(4,6) * t386 + (Icges(4,2) * t391 + t478) * t385;
t301 = Icges(4,5) * t386 + (Icges(4,1) * t388 + t477) * t385;
t310 = (Icges(4,5) * t391 - Icges(4,6) * t388) * t448;
t311 = (-Icges(4,2) * t388 + t477) * t448;
t312 = (Icges(4,1) * t391 - t478) * t448;
t399 = t386 * t310 + t312 * t464 + (-t300 * t447 + t301 * t446 + t311 * t391) * t385;
t216 = -pkin(2) * t373 - t240 + t348;
t35 = t199 * t264 + t200 * t265 + t243 * t276 + t244 * t277 + (-t242 * t375 + t263 * t472) * t385;
t4 = t35 * t386 + ((t125 * t276 + t127 * t277 + t199 * t207 + t200 * t209) * t373 + t61 * t469 - (t276 * t126 + t277 * t128 + t199 * t206 + t200 * t208) * t375 + t60 * t472 + ((-t123 * t375 + t205 * t472) * t373 - (-t124 * t375 + t204 * t472) * t375) * t385) * t385;
t398 = (-t2 - t4) * t468 + t434 + t492;
t159 = -t210 + t407;
t140 = -t375 * t280 - t373 * t335 - t193;
t49 = -t280 * t469 + t373 * t426 + t107 + t296;
t139 = (-t362 + (-rSges(4,3) - pkin(8)) * t471) * t381 + t415;
t396 = t401 + t57 + (t23 + t34) * t433 + (t22 + t35) * t432 + (t113 + t77) * t418 + (t114 + t78) * t417;
t50 = t262 + t426 * t375 + (-rSges(6,3) * t465 - t321) * t373 + t414;
t395 = t399 + t412 + t413;
t81 = pkin(3) * t400 + t407 * t381 + t129;
t82 = (-t368 * t381 - t419) * t375 - t130 + t451;
t234 = Icges(4,5) * t306 + Icges(4,6) * t305 - Icges(4,3) * t468;
t236 = Icges(4,4) * t306 + Icges(4,2) * t305 - Icges(4,6) * t468;
t238 = Icges(4,1) * t306 + Icges(4,4) * t305 - Icges(4,5) * t468;
t121 = t386 * t234 + (t236 * t391 + t238 * t388) * t385;
t235 = -Icges(4,5) * t404 + Icges(4,6) * t307 + Icges(4,3) * t471;
t237 = -Icges(4,4) * t404 + Icges(4,2) * t307 + Icges(4,6) * t471;
t239 = -Icges(4,1) * t404 + Icges(4,4) * t307 + Icges(4,5) * t471;
t122 = t386 * t235 + (t237 * t391 + t239 * t388) * t385;
t142 = t399 * t386;
t299 = Icges(4,3) * t386 + (Icges(4,5) * t388 + Icges(4,6) * t391) * t385;
t155 = -t299 * t468 + t300 * t305 + t301 * t306;
t156 = t299 * t471 + t300 * t307 - t301 * t404;
t148 = Icges(4,5) * t253 + Icges(4,6) * t252 + Icges(4,3) * t439;
t150 = Icges(4,4) * t253 + Icges(4,2) * t252 + Icges(4,6) * t439;
t152 = Icges(4,1) * t253 + Icges(4,4) * t252 + Icges(4,5) * t439;
t38 = t148 * t386 + (t150 * t391 + t152 * t388 + (-t236 * t388 + t238 * t391) * qJD(3)) * t385;
t147 = Icges(4,5) * t251 + Icges(4,6) * t250 + Icges(4,3) * t438;
t149 = Icges(4,4) * t251 + Icges(4,2) * t250 + Icges(4,6) * t438;
t151 = Icges(4,1) * t251 + Icges(4,4) * t250 + Icges(4,5) * t438;
t39 = t147 * t386 + (t149 * t391 + t151 * t388 + (-t237 * t388 + t239 * t391) * qJD(3)) * t385;
t70 = t250 * t300 + t251 * t301 + t307 * t311 - t404 * t312 + (t299 * t469 + t310 * t373) * t385;
t71 = t252 * t300 + t253 * t301 + t305 * t311 + t306 * t312 + (t299 * t472 - t310 * t375) * t385;
t394 = t142 + t396 + (t39 + t70) * t433 + (t38 + t71) * t432 + (t121 + t155) * t418 + (t122 + t156) * t417;
t378 = t392 * pkin(1);
t315 = t330 + t378;
t314 = t329 - t483;
t313 = (rSges(4,1) * t391 - rSges(4,2) * t388) * t448;
t304 = t386 * rSges(4,3) + (rSges(4,1) * t388 + rSges(4,2) * t391) * t385;
t287 = t303 - t443;
t286 = t302 - t440;
t275 = t309 * t439;
t249 = t266 * t439;
t245 = t386 * t256;
t213 = t378 + t217;
t212 = t216 - t483;
t203 = t386 * t211;
t179 = -t240 * t386 - t304 * t468;
t178 = t241 * t386 - t304 * t471;
t158 = t160 + t378;
t157 = t159 - t483;
t154 = rSges(4,3) * t439 - t415;
t146 = -t210 * t386 - t266 * t468;
t145 = -t266 * t471 + t203;
t137 = -t193 * t386 - t261 * t468;
t136 = -t261 * t471 + t180;
t135 = t141 + t378;
t134 = t140 - t483;
t133 = t139 - t443;
t132 = t138 - t440;
t120 = t386 * t153 + (-t304 * t469 - t313 * t373) * t385;
t119 = -t386 * t154 + (t304 * t472 - t313 * t375) * t385;
t118 = t386 * t129;
t96 = (-t210 - t255) * t386 + t375 * t427;
t95 = t373 * t427 + t203 + t245;
t88 = t235 * t471 + t237 * t307 - t239 * t404;
t87 = t234 * t471 + t236 * t307 - t238 * t404;
t86 = -t235 * t468 + t237 * t305 + t239 * t306;
t85 = -t234 * t468 + t236 * t305 + t238 * t306;
t80 = t82 - t443;
t79 = t81 - t440;
t74 = t455 + t131;
t73 = t375 * t428 + t386 * t459;
t72 = t373 * t428 + t460;
t65 = t118 + (-t266 * t469 - t476) * t385;
t64 = -t130 * t386 - t246 * t468 + t249;
t52 = t97 + (-t226 * t373 - t261 * t469) * t385;
t51 = -t108 * t386 - t226 * t468 + t232;
t48 = t50 - t443;
t47 = t49 - t440;
t46 = (-t255 + t459) * t386 + t375 * t416;
t45 = t373 * t416 + t245 + t460;
t42 = t118 + (t453 * t469 - t476) * t385 + t402;
t41 = t249 + t275 + (-t130 - t182) * t386 + (-t246 * t385 - t420) * t375;
t40 = t44 + t455;
t31 = -t211 * t439 + t437;
t26 = -t194 * t439 + t444;
t25 = (t454 * t469 + t430) * t385 + t479;
t24 = t375 * t429 + t386 * t461 + t456;
t19 = (t435 * t469 + t430) * t385 + t402 + t479;
t18 = t275 + (-t182 + t461) * t386 + (-t420 + t429) * t375 + t456;
t15 = (-t211 - t256) * t439 + t436 + t437;
t8 = t439 * t458 + t406;
t7 = (-t256 + t458) * t439 + t406 + t436;
t1 = [t395 + (t134 * t48 + t135 * t47) * t486 + (t157 * t80 + t158 * t79) * t487 + (t132 * t213 + t133 * t212) * t488 + (t286 * t315 + t287 * t314) * t489; t395 + m(6) * (t134 * t50 + t135 * t49 + t140 * t48 + t141 * t47) + m(5) * (t157 * t82 + t158 * t81 + t159 * t80 + t160 * t79) + m(4) * (t132 * t217 + t133 * t216 + t138 * t213 + t139 * t212) + m(3) * (t286 * t330 + t287 * t329 + t302 * t315 + t303 * t314); t395 + (t140 * t50 + t141 * t49) * t486 + (t159 * t82 + t160 * t81) * t487 + (t138 * t217 + t139 * t216) * t488 + (t302 * t330 + t303 * t329) * t489; t394 + m(4) * (t119 * t212 + t120 * t213 + t132 * t178 + t133 * t179) + m(5) * (t157 * t41 + t158 * t42 + t79 * t95 + t80 * t96) + m(6) * (t134 * t18 + t135 * t19 + t45 * t47 + t46 * t48); t394 + m(4) * (t119 * t216 + t120 * t217 + t138 * t178 + t139 * t179) + m(5) * (t159 * t41 + t160 * t42 + t81 * t95 + t82 * t96) + m(6) * (t140 * t18 + t141 * t19 + t45 * t49 + t46 * t50); (t156 * t386 + (t373 * t88 - t375 * t87) * t385) * t438 + (t70 * t386 + ((t307 * t149 - t151 * t404 + t250 * t237 + t251 * t239) * t373 + t88 * t469 - (t307 * t150 - t152 * t404 + t250 * t236 + t251 * t238) * t375 + t87 * t472 + ((t147 * t373 + t235 * t469) * t373 - (t148 * t373 + t234 * t469) * t375) * t385) * t385) * t471 + (t155 * t386 + (t373 * t86 - t375 * t85) * t385) * t439 + (t18 * t46 + t19 * t45 + t40 * t7) * t486 + (t15 * t74 + t41 * t96 + t42 * t95) * t487 + t386 * (t142 + ((t122 * t381 - t38) * t375 + (t121 * t381 + t39) * t373) * t385) + (t179 * t119 + t178 * t120 + (t240 * t373 + t241 * t375) * ((t240 * t381 + t153) * t375 + (-t241 * t381 + t154) * t373) * t379) * t488 + t403 + (-t4 - t71 * t386 - ((t305 * t149 + t306 * t151 + t252 * t237 + t253 * t239) * t373 + t86 * t469 - (t305 * t150 + t306 * t152 + t252 * t236 + t253 * t238) * t375 + t85 * t472 + ((-t147 * t375 + t235 * t472) * t373 - (-t148 * t375 + t234 * t472) * t375) * t385) * t385) * t468 + t492; m(6) * (t134 * t24 + t135 * t25 + t47 * t72 + t48 * t73) + m(5) * (t145 * t79 + t146 * t80 + t157 * t64 + t158 * t65) + t396; t396 + m(6) * (t140 * t24 + t141 * t25 + t49 * t72 + t50 * t73) + m(5) * (t145 * t81 + t146 * t82 + t159 * t64 + t160 * t65); m(6) * (t18 * t73 + t19 * t72 + t24 * t46 + t25 * t45 + t40 * t8 + t44 * t7) + m(5) * (t131 * t15 + t145 * t42 + t146 * t41 + t31 * t74 + t64 * t96 + t65 * t95) + t398; (t131 * t31 + t145 * t65 + t146 * t64) * t487 + (t24 * t73 + t25 * t72 + t44 * t8) * t486 + t398; m(6) * (t134 * t51 + t135 * t52 + t136 * t47 + t137 * t48) + t401; m(6) * (t136 * t49 + t137 * t50 + t140 * t51 + t141 * t52) + t401; m(6) * (t117 * t7 + t136 * t19 + t137 * t18 + t26 * t40 + t45 * t52 + t46 * t51) + t403; m(6) * (t117 * t8 + t136 * t25 + t137 * t24 + t26 * t44 + t51 * t73 + t52 * t72) + t403; (t117 * t26 + t136 * t52 + t137 * t51) * t486 + t403;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
