% Calculate vector of inverse dynamics joint torques for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:08
% EndTime: 2022-01-23 09:14:26
% DurationCPUTime: 14.51s
% Computational Cost: add. (14957->631), mult. (10691->798), div. (0->0), fcn. (8283->10), ass. (0->346)
t278 = qJ(1) + pkin(8);
t269 = sin(t278);
t276 = pkin(9) + qJ(4);
t270 = cos(t276);
t417 = t269 * t270;
t374 = rSges(5,1) * t417;
t271 = cos(t278);
t281 = -pkin(6) - qJ(3);
t235 = t271 * t281;
t280 = cos(pkin(9));
t266 = t280 * pkin(3) + pkin(2);
t392 = t269 * t266 + t235;
t282 = sin(qJ(1));
t464 = pkin(1) * t282;
t497 = -t374 - t464 - t392;
t283 = cos(qJ(1));
t274 = t283 * pkin(1);
t198 = rSges(3,1) * t269 + rSges(3,2) * t271;
t168 = -t198 - t464;
t285 = qJD(1) ^ 2;
t496 = t285 * t274;
t255 = t269 * rSges(6,3);
t272 = qJ(5) + t276;
t265 = cos(t272);
t420 = t265 * t271;
t394 = rSges(6,1) * t420 + t255;
t264 = sin(t272);
t422 = t264 * t271;
t129 = -rSges(6,2) * t422 + t394;
t389 = -t271 * t266 + t269 * t281;
t263 = pkin(4) * t270;
t211 = t263 + t266;
t275 = pkin(7) - t281;
t228 = t275 * t269;
t486 = t271 * t211 + t228;
t489 = t486 + t389;
t446 = t129 + t489;
t247 = t271 * qJ(3);
t197 = pkin(2) * t269 - t247;
t120 = t197 - t392;
t186 = qJD(1) * t197;
t244 = qJD(3) * t269;
t495 = -qJD(1) * t120 + t186 + t244;
t423 = t264 * t269;
t209 = rSges(6,2) * t423;
t395 = t271 * rSges(6,3) + t209;
t421 = t265 * t269;
t128 = rSges(6,1) * t421 - t395;
t229 = t275 * t271;
t97 = t211 * t269 - t229 - t392;
t494 = -t128 - t97;
t246 = t269 * qJ(3);
t200 = t271 * pkin(2) + t246;
t121 = -t200 - t389;
t359 = t200 + t274;
t341 = t121 + t359;
t493 = -t341 - t446;
t381 = qJD(4) * t271;
t268 = sin(t276);
t492 = rSges(5,2) * t268;
t240 = Icges(6,4) * t265;
t179 = Icges(6,1) * t264 + t240;
t322 = -Icges(6,2) * t264 + t240;
t491 = t179 + t322;
t490 = t274 + t486;
t454 = rSges(6,2) * t265;
t459 = rSges(6,1) * t264;
t181 = t454 + t459;
t150 = t181 * t269;
t151 = t181 * t271;
t277 = qJD(4) + qJD(5);
t187 = t269 * t277;
t188 = t271 * t277;
t37 = t128 * t187 + t129 * t188 + qJD(2) + (t269 * t97 + t271 * t489) * qJD(4);
t360 = -t197 - t464;
t342 = t120 + t360;
t309 = t342 + t494;
t369 = t268 * t381;
t344 = pkin(4) * t369;
t346 = -t181 * t188 + t244;
t42 = qJD(1) * t309 - t344 + t346;
t453 = pkin(4) * qJD(4);
t371 = t268 * t453;
t207 = t269 * t371;
t245 = qJD(3) * t271;
t396 = t207 + t245;
t337 = t181 * t187 + t396;
t43 = -t493 * qJD(1) - t337;
t241 = t265 * rSges(6,1);
t455 = rSges(6,2) * t264;
t484 = t241 - t455;
t378 = qJD(1) * qJD(4);
t183 = qJDD(4) * t269 + t271 * t378;
t377 = qJD(1) * qJD(5);
t130 = qJDD(5) * t269 + t271 * t377 + t183;
t230 = t269 * t378;
t131 = t269 * t377 + t230 + (-qJDD(4) - qJDD(5)) * t271;
t184 = -qJDD(4) * t271 + t230;
t372 = t277 * t454;
t384 = qJD(1) * t269;
t383 = qJD(1) * t271;
t398 = rSges(6,3) * t383 + qJD(1) * t209;
t78 = -t271 * t372 + (-t188 * t264 - t265 * t384) * rSges(6,1) + t398;
t79 = -t277 * t150 + (t271 * t484 + t255) * qJD(1);
t222 = t275 * t383;
t393 = t211 - t266;
t80 = -t344 + t222 + (-t269 * t393 + t235) * qJD(1);
t223 = t281 * t384;
t81 = -t207 + t223 + (t271 * t393 + t228) * qJD(1);
t8 = t128 * t130 - t129 * t131 + t183 * t97 - t184 * t489 + t187 * t79 + t188 * t78 + qJDD(2) + (t269 * t81 + t271 * t80) * qJD(4);
t488 = -t42 * (qJD(1) * t150 - t188 * t484) - t37 * (-t187 * t150 - t151 * t188) - t43 * (-qJD(1) * t151 - t187 * t484) + t8 * (t269 * t128 + t271 * t129);
t279 = sin(pkin(9));
t457 = rSges(4,2) * t279;
t460 = rSges(4,1) * t280;
t143 = t269 * rSges(4,3) + (-t457 + t460) * t271;
t104 = t143 + t359;
t201 = t271 * rSges(3,1) - rSges(3,2) * t269;
t169 = t201 + t274;
t419 = t268 * t269;
t485 = -rSges(5,2) * t419 - t271 * rSges(5,3);
t254 = Icges(5,4) * t270;
t323 = -Icges(5,2) * t268 + t254;
t193 = Icges(5,1) * t268 + t254;
t190 = Icges(5,5) * t270 - Icges(5,6) * t268;
t189 = Icges(5,5) * t268 + Icges(5,6) * t270;
t303 = qJD(4) * t189;
t443 = Icges(5,4) * t268;
t194 = Icges(5,1) * t270 - t443;
t138 = Icges(5,5) * t269 + t194 * t271;
t136 = Icges(5,6) * t269 + t271 * t323;
t429 = t136 * t268;
t318 = -t138 * t270 + t429;
t435 = Icges(5,3) * t271;
t482 = -t271 * t303 + (-t190 * t269 + t318 + t435) * qJD(1);
t216 = Icges(5,4) * t419;
t441 = Icges(5,5) * t271;
t137 = Icges(5,1) * t417 - t216 - t441;
t437 = Icges(5,6) * t271;
t135 = Icges(5,4) * t417 - Icges(5,2) * t419 - t437;
t430 = t135 * t268;
t319 = -t137 * t270 + t430;
t134 = Icges(5,3) * t269 + t190 * t271;
t386 = qJD(1) * t134;
t481 = qJD(1) * t319 - t269 * t303 + t386;
t442 = Icges(6,4) * t264;
t180 = Icges(6,1) * t265 - t442;
t127 = Icges(6,5) * t269 + t180 * t271;
t176 = Icges(6,5) * t265 - Icges(6,6) * t264;
t175 = Icges(6,5) * t264 + Icges(6,6) * t265;
t427 = t175 * t271;
t125 = Icges(6,6) * t269 + t271 * t322;
t432 = t125 * t264;
t434 = Icges(6,3) * t271;
t480 = -t277 * t427 + (-t127 * t265 - t176 * t269 + t432 + t434) * qJD(1);
t436 = Icges(6,6) * t271;
t124 = Icges(6,4) * t421 - Icges(6,2) * t423 - t436;
t206 = Icges(6,4) * t423;
t440 = Icges(6,5) * t271;
t126 = Icges(6,1) * t421 - t206 - t440;
t321 = t124 * t264 - t126 * t265;
t123 = Icges(6,3) * t269 + t176 * t271;
t387 = qJD(1) * t123;
t428 = t175 * t269;
t479 = qJD(1) * t321 - t277 * t428 + t387;
t133 = Icges(5,5) * t417 - Icges(5,6) * t419 - t435;
t53 = -t133 * t271 - t269 * t319;
t177 = Icges(6,2) * t265 + t442;
t316 = t177 * t264 - t179 * t265;
t478 = qJD(1) * t316 + t176 * t277;
t191 = Icges(5,2) * t270 + t443;
t314 = t191 * t268 - t193 * t270;
t477 = t314 * qJD(1) + t190 * qJD(4);
t476 = t269 * (-t191 * t271 + t138) - t271 * (-Icges(5,2) * t417 + t137 - t216);
t475 = qJD(1) * t491 + t187 * (-t177 * t271 + t127) - t188 * (-Icges(6,2) * t421 + t126 - t206);
t474 = t130 / 0.2e1;
t473 = t131 / 0.2e1;
t472 = t183 / 0.2e1;
t471 = t184 / 0.2e1;
t470 = -t187 / 0.2e1;
t469 = t187 / 0.2e1;
t468 = -t188 / 0.2e1;
t467 = t188 / 0.2e1;
t466 = t269 / 0.2e1;
t465 = -t271 / 0.2e1;
t463 = -qJD(1) / 0.2e1;
t462 = qJD(1) / 0.2e1;
t461 = pkin(2) - t266;
t456 = rSges(5,2) * t270;
t196 = rSges(5,1) * t268 + t456;
t161 = t196 * t271;
t256 = t269 * rSges(5,3);
t416 = t270 * t271;
t418 = t268 * t271;
t140 = rSges(5,1) * t416 - rSges(5,2) * t418 + t256;
t382 = qJD(4) * t269;
t58 = -t196 * t382 - t245 + (t140 + t341) * qJD(1);
t452 = t161 * t58;
t451 = t268 * t43;
t139 = t374 + t485;
t312 = -t139 + t342;
t335 = -t196 * t381 + t244;
t57 = qJD(1) * t312 + t335;
t450 = t269 * t57;
t449 = t271 * t42;
t448 = t42 * t181;
t447 = qJDD(1) / 0.2e1;
t122 = Icges(6,5) * t421 - Icges(6,6) * t423 - t434;
t433 = t122 * t271;
t426 = t177 * t277;
t425 = t189 * t269;
t424 = t189 * t271;
t415 = t270 * qJD(4) ^ 2;
t62 = -t269 * t316 - t427;
t413 = t62 * qJD(1);
t76 = -t269 * t314 - t424;
t412 = t76 * qJD(1);
t411 = -t269 * t122 - t126 * t420;
t410 = t269 * t123 + t127 * t420;
t409 = -t269 * t133 - t137 * t416;
t408 = t269 * t134 + t138 * t416;
t159 = qJD(1) * t200 - t245;
t407 = t223 - (-t271 * t461 - t246) * qJD(1) - t159;
t400 = -t191 + t194;
t399 = t193 + t323;
t397 = rSges(5,3) * t383 + t384 * t492;
t224 = t269 * t457;
t391 = rSges(4,3) * t383 + qJD(1) * t224;
t390 = t271 * rSges(4,3) + t224;
t236 = qJ(3) * t383;
t388 = t236 + t244;
t385 = qJD(1) * t190;
t380 = -m(4) - m(5) - m(6);
t379 = qJD(1) * qJD(3);
t376 = t128 * t383 + t269 * t79 + t271 * t78;
t375 = t269 * t460;
t368 = -pkin(2) - t460;
t366 = t384 / 0.2e1;
t365 = t383 / 0.2e1;
t364 = -t382 / 0.2e1;
t363 = t382 / 0.2e1;
t362 = -t381 / 0.2e1;
t361 = t381 / 0.2e1;
t302 = -pkin(4) * t268 - t181;
t358 = -t211 - t241;
t308 = t179 * t277;
t357 = qJD(1) * t127 - t124 * t277 - t269 * t308;
t356 = -t125 * t277 - t271 * t308 + (-t180 * t269 + t440) * qJD(1);
t355 = qJD(1) * t125 + t126 * t277 - t269 * t426;
t354 = t127 * t277 - t271 * t426 + (-t269 * t322 + t436) * qJD(1);
t100 = t127 * t421;
t353 = t123 * t271 - t100;
t107 = t138 * t417;
t352 = t134 * t271 - t107;
t351 = -t122 + t432;
t350 = -t133 + t429;
t349 = t491 * t277;
t348 = t180 * t277 - t426;
t343 = qJDD(1) * t274 - t285 * t464;
t142 = t375 - t390;
t340 = -t142 + t360;
t339 = t395 - t464;
t165 = t484 * t277;
t338 = -t270 * t453 - t165;
t227 = rSges(2,1) * t283 - rSges(2,2) * t282;
t226 = rSges(2,1) * t282 + rSges(2,2) * t283;
t199 = rSges(5,1) * t270 - t492;
t68 = t136 * t270 + t138 * t268;
t304 = qJD(4) * t191;
t86 = -t271 * t304 + (-t269 * t323 + t437) * qJD(1);
t305 = qJD(4) * t193;
t88 = -t271 * t305 + (-t194 * t269 + t441) * qJD(1);
t289 = -qJD(4) * t68 - t268 * t86 + t270 * t88 + t386;
t67 = t135 * t270 + t137 * t268;
t87 = qJD(1) * t136 - t269 * t304;
t89 = qJD(1) * t138 - t269 * t305;
t290 = qJD(1) * t133 - qJD(4) * t67 - t268 * t87 + t270 * t89;
t332 = -(t269 * t481 + t290 * t271) * t271 + (t269 * t482 + t289 * t271) * t269;
t331 = -(t290 * t269 - t271 * t481) * t271 + (t289 * t269 - t271 * t482) * t269;
t330 = -t269 * t43 - t449;
t54 = -t136 * t419 - t352;
t329 = t269 * t54 - t271 * t53;
t55 = -t135 * t418 - t409;
t56 = -t136 * t418 + t408;
t328 = t269 * t56 - t271 * t55;
t327 = -t269 * t58 - t271 * t57;
t90 = -t381 * t456 + (-t270 * t384 - t369) * rSges(5,1) + t397;
t160 = t196 * t269;
t91 = -qJD(4) * t160 + (t199 * t271 + t256) * qJD(1);
t326 = t269 * t91 + t271 * t90;
t60 = t124 * t265 + t126 * t264;
t317 = t139 * t269 + t140 * t271;
t315 = t191 * t270 + t193 * t268;
t313 = qJDD(3) * t269 + t271 * t379 - t496;
t306 = t321 * t269;
t301 = qJD(1) * t176 - t187 * t427 + t188 * t428;
t300 = t135 * t271 - t136 * t269;
t299 = (-t268 * t399 + t270 * t400) * qJD(1);
t293 = qJD(1) * t122 - t264 * t355 + t265 * t357;
t11 = t269 * t479 + t293 * t271;
t292 = -t264 * t354 + t265 * t356 + t387;
t12 = t269 * t480 + t292 * t271;
t13 = t293 * t269 - t271 * t479;
t14 = t292 * t269 - t271 * t480;
t47 = -t306 - t433;
t48 = -t125 * t423 - t353;
t22 = t187 * t48 - t188 * t47 + t413;
t49 = -t124 * t422 - t411;
t50 = -t125 * t422 + t410;
t63 = -t271 * t316 + t428;
t59 = t63 * qJD(1);
t23 = t187 * t50 - t188 * t49 + t59;
t294 = (-t179 * t271 - t125) * t187 - (-t179 * t269 - t124) * t188 + (-t177 + t180) * qJD(1);
t286 = -t264 * t475 + t294 * t265;
t291 = qJD(1) * t175 - t264 * t349 + t265 * t348;
t30 = t269 * t478 + t291 * t271;
t31 = t291 * t269 - t271 * t478;
t32 = t264 * t357 + t265 * t355;
t33 = t264 * t356 + t265 * t354;
t61 = t125 * t265 + t127 * t264;
t298 = (qJD(1) * t30 + qJDD(1) * t63 - t11 * t188 + t12 * t187 + t130 * t50 + t131 * t49) * t466 + (t294 * t264 + t265 * t475) * t463 + (qJD(1) * t31 + qJDD(1) * t62 - t13 * t188 + t130 * t48 + t131 * t47 + t14 * t187) * t465 + t22 * t366 + t23 * t365 + (t269 * t50 - t271 * t49) * t474 + (t269 * t48 - t271 * t47) * t473 + (-t11 * t271 + t12 * t269 + (t269 * t49 + t271 * t50) * qJD(1)) * t469 + (-t13 * t271 + t14 * t269 + (t269 * t47 + t271 * t48) * qJD(1)) * t468 + (t269 * t61 - t271 * t60) * t447 + (t269 * t33 - t271 * t32 + (t269 * t60 + t271 * t61) * qJD(1)) * t462 + (t269 * t301 + t271 * t286) * t470 + (t269 * t286 - t271 * t301) * t467;
t297 = -qJDD(3) * t271 + qJD(1) * (-pkin(2) * t384 + t388) + qJDD(1) * t200 + t269 * t379 + t343;
t295 = qJD(1) * (-t236 + (t269 * t461 - t235) * qJD(1)) + qJDD(1) * t121 + t297;
t172 = t323 * qJD(4);
t173 = t194 * qJD(4);
t288 = qJD(1) * t189 - qJD(4) * t315 - t172 * t268 + t173 * t270;
t287 = -t268 * t476 + t300 * t270;
t174 = t199 * qJD(4);
t93 = qJD(1) * t104 - t245;
t92 = qJD(1) * t340 + t244;
t77 = -t271 * t314 + t425;
t69 = t77 * qJD(1);
t64 = qJD(4) * t317 + qJD(2);
t46 = qJDD(1) * t143 + qJD(1) * (-qJD(1) * t375 + t391) + t297;
t45 = t340 * qJDD(1) + (-qJD(1) * t143 - t159) * qJD(1) + t313;
t39 = t288 * t269 - t271 * t477;
t38 = t269 * t477 + t288 * t271;
t36 = -qJD(4) * t318 + t268 * t88 + t270 * t86;
t35 = -t319 * qJD(4) + t268 * t89 + t270 * t87;
t34 = qJD(4) * t326 + t139 * t183 - t140 * t184 + qJDD(2);
t29 = qJD(1) * t90 + qJDD(1) * t140 - t174 * t382 - t183 * t196 + t295;
t28 = -t174 * t381 + t184 * t196 + (-t91 + t407) * qJD(1) + t312 * qJDD(1) + t313;
t25 = qJD(4) * t328 + t69;
t24 = qJD(4) * t329 + t412;
t10 = -t130 * t181 - t165 * t187 + t446 * qJDD(1) + (t78 + t80) * qJD(1) + (-t183 * t268 - t269 * t415) * pkin(4) + t295;
t9 = t131 * t181 - t165 * t188 + (t184 * t268 - t271 * t415) * pkin(4) + (-t79 - t81 + t407) * qJD(1) + t309 * qJDD(1) + t313;
t1 = [(t59 + (t48 + (t124 * t271 + t125 * t269) * t264 + t353 + t411) * t188 + (-t126 * t421 + t433 + t47 + (t124 * t269 - t125 * t271) * t264 + t410) * t187) * t467 - m(2) * (-g(1) * t226 + g(2) * t227) + (t69 + ((t54 - t107 + (t134 + t430) * t271 + t409) * t271 + t408 * t269) * qJD(4)) * t361 + (t61 + t63) * t474 + (t60 + t62) * t473 + (t68 + t77) * t472 + (t67 + t76) * t471 + (-t413 + (t50 - t306 - t410) * t188 + (t351 * t269 - t100 + t49) * t187 + ((t123 + t321) * t187 + t351 * t188) * t271 + t22) * t470 + (t33 + t30) * t469 + (t24 - t412 + ((t271 * t350 - t408 + t56) * t271 + (t269 * t350 + t352 + t55) * t269) * qJD(4)) * t364 + (t36 + t38) * t363 + (-qJD(4) * t314 + t172 * t270 + t173 * t268 + t264 * t348 + t265 * t349) * qJD(1) + ((-t198 * t285 - g(2) + t343) * t169 + (-t496 + (-0.2e1 * t201 - t274 + t169) * t285 - g(1)) * t168) * m(3) + (t32 + t31 + t23) * t468 + (t35 + t39 + t25) * t362 + (t9 * t339 + t10 * (t394 + t490) + (-t10 * t455 + t9 * t275) * t271 + (t277 * t448 + t9 * t358) * t269 - g(1) * (t269 * t358 + t229 + t339) - g(2) * (t129 + t490) + pkin(4) * t451 * t381 + (t222 + t398 + (-t277 * t459 - t371 - t372) * t271 - t346 + t495) * t43 + (t396 - t337) * t42 + ((-t282 * t43 - t283 * t42) * pkin(1) + (-t211 - t484) * t449 + (t42 * (-rSges(6,3) - t275) + t43 * t358) * t269 - t42 * t493 - t43 * (-t464 + t494)) * qJD(1)) * m(6) + ((t196 * t450 - t452) * qJD(4) + (-g(2) + t29) * (t140 + t274 - t389) + (-g(1) + t28) * (-t485 + t497) + (t223 + t245 + (-t256 - t274 + (-t199 - t266) * t271) * qJD(1)) * t57 + (t397 + t57 - t335 + (t139 + t464 + t497) * qJD(1) + t495) * t58) * m(5) + (-(-t186 + t244 - t92 + (-t142 - t464) * qJD(1)) * t93 + t92 * t245 + t93 * (t388 + t391) + ((-t282 * t93 - t283 * t92) * pkin(1) + t92 * (t368 + t457) * t271 + (t92 * (-rSges(4,3) - qJ(3)) + t93 * t368) * t269) * qJD(1) + (-g(2) + t46) * t104 + (-g(1) + t45) * (t368 * t269 + t247 + t390 - t464)) * m(4) + (Icges(4,2) * t280 ^ 2 + (Icges(4,1) * t279 + 0.2e1 * Icges(4,4) * t280) * t279 + t315 + t177 * t265 + t179 * t264 + m(2) * (t226 ^ 2 + t227 ^ 2) + m(3) * (t168 ^ 2 + t201 * t169) + Icges(2,3) + Icges(3,3)) * qJDD(1); (m(3) + m(4)) * qJDD(2) + m(5) * t34 + m(6) * t8 + (-m(3) + t380) * g(3); t380 * (g(1) * t269 - g(2) * t271) + 0.2e1 * (t10 * t465 + t466 * t9) * m(6) + 0.2e1 * (t28 * t466 + t29 * t465) * m(5) + 0.2e1 * (t45 * t466 + t46 * t465) * m(4); t328 * t472 + t329 * t471 + (t269 * t68 - t271 * t67) * t447 + ((t268 * t400 + t270 * t399) * qJD(1) + (t300 * t268 + t270 * t476) * qJD(4)) * t463 + (t269 * t36 - t271 * t35 + (t67 * t269 + t271 * t68) * qJD(1)) * t462 + ((t55 * t269 + t56 * t271) * qJD(1) + t332) * t363 + ((t53 * t269 + t54 * t271) * qJD(1) + t331) * t362 + (qJD(1) * t39 + qJD(4) * t331 + qJDD(1) * t76 + t183 * t54 + t184 * t53) * t465 + t298 + ((-t381 * t425 - t385) * t271 + (t299 + (t271 * t424 + t287) * qJD(4)) * t269) * t361 + ((-t382 * t424 + t385) * t269 + (t299 + (t269 * t425 + t287) * qJD(4)) * t271) * t364 + (qJD(1) * t38 + qJD(4) * t332 + qJDD(1) * t77 + t183 * t56 + t184 * t55) * t466 + t24 * t366 + t25 * t365 + (-g(3) * (t484 + t263) - (g(1) * t271 + g(2) * t269) * t302 - (-t383 * t451 + (t330 * t270 + t37 * (-t269 ^ 2 - t271 ^ 2) * t268) * qJD(4)) * pkin(4) + t37 * t376 + (t9 * t302 + t42 * t338 + t8 * t489 + t37 * t80 + (t302 * t43 + t37 * t97) * qJD(1)) * t271 + (t10 * t302 + t43 * t338 + t8 * t97 + t37 * t81 + (-t37 * t446 + t448) * qJD(1)) * t269 + t488) * m(6) + (g(1) * t161 + g(2) * t160 - g(3) * t199 + t34 * t317 + t64 * ((t139 * t271 - t140 * t269) * qJD(1) + t326) + t327 * t174 + (-t29 * t269 - t28 * t271 + (-t271 * t58 + t450) * qJD(1)) * t196 - (t160 * t57 - t452) * qJD(1) - (t64 * (-t160 * t269 - t161 * t271) + t327 * t199) * qJD(4)) * m(5); t298 + (g(1) * t151 + g(2) * t150 - g(3) * t484 + t37 * (-t129 * t384 + t376) + t330 * t165 + (-t10 * t269 - t9 * t271 + (t269 * t42 - t271 * t43) * qJD(1)) * t181 + t488) * m(6);];
tau = t1;
