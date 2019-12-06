% Calculate vector of inverse dynamics joint torques for
% S5PRPRR1
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:38
% EndTime: 2019-12-05 15:42:53
% DurationCPUTime: 8.42s
% Computational Cost: add. (14784->597), mult. (10369->775), div. (0->0), fcn. (8119->8), ass. (0->329)
t262 = pkin(9) + qJ(4);
t255 = sin(t262);
t464 = rSges(5,2) * t255;
t259 = qJ(5) + t262;
t253 = cos(t259);
t229 = Icges(6,4) * t253;
t252 = sin(t259);
t178 = Icges(6,1) * t252 + t229;
t306 = -Icges(6,2) * t252 + t229;
t463 = t178 + t306;
t263 = pkin(8) + qJ(2);
t256 = sin(t263);
t401 = t252 * t256;
t208 = rSges(6,2) * t401;
t258 = cos(t263);
t399 = t253 * t256;
t128 = rSges(6,1) * t399 - rSges(6,3) * t258 - t208;
t244 = t256 * rSges(6,3);
t398 = t253 * t258;
t400 = t252 * t258;
t129 = rSges(6,1) * t398 - rSges(6,2) * t400 + t244;
t431 = rSges(6,2) * t253;
t435 = rSges(6,1) * t252;
t180 = t431 + t435;
t150 = t180 * t256;
t151 = t180 * t258;
t264 = qJD(4) + qJD(5);
t186 = t256 * t264;
t187 = t258 * t264;
t266 = cos(pkin(9));
t254 = pkin(3) * t266 + pkin(2);
t267 = -pkin(6) - qJ(3);
t261 = -pkin(7) + t267;
t257 = cos(t262);
t251 = pkin(4) * t257;
t351 = t251 + t254;
t377 = -t256 * t351 - t258 * t261;
t391 = t258 * t267;
t98 = t254 * t256 + t377 + t391;
t167 = t258 * t351;
t213 = t258 * t254;
t364 = t267 - t261;
t99 = t256 * t364 + t167 - t213;
t37 = t128 * t186 + t129 * t187 + qJD(1) + (-t256 * t98 + t258 * t99) * qJD(4);
t233 = qJD(3) * t256;
t357 = qJD(4) * t258;
t344 = t255 * t357;
t321 = pkin(4) * t344;
t283 = -t180 * t187 + t233 - t321;
t236 = t258 * qJ(3);
t438 = pkin(2) - t254;
t294 = t256 * t438 - t391;
t120 = -t236 + t294;
t196 = pkin(2) * t256 - t236;
t383 = t120 - t196;
t425 = -t128 + t98;
t320 = t383 + t425;
t42 = qJD(2) * t320 + t283;
t430 = pkin(4) * qJD(4);
t347 = t255 * t430;
t206 = t256 * t347;
t234 = qJD(3) * t258;
t370 = -t206 - t234;
t235 = t256 * qJ(3);
t199 = pkin(2) * t258 + t235;
t324 = -t256 * t267 + t213;
t121 = t324 - t199;
t382 = t121 + t199;
t424 = t129 + t99;
t43 = -t180 * t186 + (t382 + t424) * qJD(2) + t370;
t230 = t253 * rSges(6,1);
t432 = rSges(6,2) * t252;
t459 = t230 - t432;
t354 = qJD(2) * qJD(4);
t182 = qJDD(4) * t256 + t258 * t354;
t353 = qJD(2) * qJD(5);
t130 = qJDD(5) * t256 + t258 * t353 + t182;
t221 = t256 * t354;
t131 = t256 * t353 + t221 + (-qJDD(4) - qJDD(5)) * t258;
t183 = -qJDD(4) * t258 + t221;
t348 = t264 * t431;
t360 = qJD(2) * t256;
t359 = qJD(2) * t258;
t373 = rSges(6,3) * t359 + qJD(2) * t208;
t78 = -t258 * t348 + (-t187 * t252 - t253 * t360) * rSges(6,1) + t373;
t295 = t180 * t264;
t79 = -t256 * t295 + (t258 * t459 + t244) * qJD(2);
t323 = t254 - t351;
t80 = -t321 + (t256 * t323 + t258 * t364) * qJD(2);
t217 = t267 * t360;
t394 = t256 * t261;
t81 = t217 - t206 + (-t258 * t323 - t394) * qJD(2);
t8 = t128 * t130 - t129 * t131 - t182 * t98 - t183 * t99 + t186 * t79 + t187 * t78 + qJDD(1) + (t256 * t81 + t258 * t80) * qJD(4);
t462 = -t42 * (qJD(2) * t150 - t187 * t459) - t37 * (-t150 * t186 - t151 * t187) - t43 * (-qJD(2) * t151 - t186 * t459) + t8 * (t128 * t256 + t129 * t258);
t461 = -m(2) - m(3);
t185 = qJD(2) * t196;
t460 = qJD(2) * t120 - t185;
t243 = Icges(5,4) * t257;
t307 = -Icges(5,2) * t255 + t243;
t192 = Icges(5,1) * t255 + t243;
t265 = sin(pkin(9));
t434 = rSges(4,2) * t265;
t437 = rSges(4,1) * t266;
t143 = t256 * rSges(4,3) + (-t434 + t437) * t258;
t189 = Icges(5,5) * t257 - Icges(5,6) * t255;
t188 = Icges(5,5) * t255 + Icges(5,6) * t257;
t287 = qJD(4) * t188;
t421 = Icges(5,4) * t255;
t193 = Icges(5,1) * t257 - t421;
t138 = Icges(5,5) * t256 + t193 * t258;
t136 = Icges(5,6) * t256 + t258 * t307;
t407 = t136 * t255;
t302 = -t138 * t257 + t407;
t413 = Icges(5,3) * t258;
t458 = -t258 * t287 + (-t189 * t256 + t302 + t413) * qJD(2);
t397 = t255 * t256;
t212 = Icges(5,4) * t397;
t395 = t256 * t257;
t419 = Icges(5,5) * t258;
t137 = Icges(5,1) * t395 - t212 - t419;
t415 = Icges(5,6) * t258;
t135 = Icges(5,4) * t395 - Icges(5,2) * t397 - t415;
t408 = t135 * t255;
t303 = -t137 * t257 + t408;
t134 = Icges(5,3) * t256 + t189 * t258;
t362 = qJD(2) * t134;
t457 = qJD(2) * t303 - t256 * t287 + t362;
t420 = Icges(6,4) * t252;
t179 = Icges(6,1) * t253 - t420;
t127 = Icges(6,5) * t256 + t179 * t258;
t175 = Icges(6,5) * t253 - Icges(6,6) * t252;
t174 = Icges(6,5) * t252 + Icges(6,6) * t253;
t291 = t174 * t264;
t125 = Icges(6,6) * t256 + t258 * t306;
t410 = t125 * t252;
t412 = Icges(6,3) * t258;
t456 = -t258 * t291 + (-t127 * t253 - t175 * t256 + t410 + t412) * qJD(2);
t414 = Icges(6,6) * t258;
t124 = Icges(6,4) * t399 - Icges(6,2) * t401 - t414;
t205 = Icges(6,4) * t401;
t418 = Icges(6,5) * t258;
t126 = Icges(6,1) * t399 - t205 - t418;
t305 = t124 * t252 - t126 * t253;
t123 = Icges(6,3) * t256 + t175 * t258;
t363 = qJD(2) * t123;
t455 = qJD(2) * t305 - t256 * t291 + t363;
t133 = Icges(5,5) * t395 - Icges(5,6) * t397 - t413;
t53 = -t133 * t258 - t256 * t303;
t176 = Icges(6,2) * t253 + t420;
t300 = t176 * t252 - t178 * t253;
t454 = qJD(2) * t300 + t175 * t264;
t190 = Icges(5,2) * t257 + t421;
t298 = t190 * t255 - t192 * t257;
t453 = t298 * qJD(2) + qJD(4) * t189;
t452 = t256 * (-t190 * t258 + t138) - t258 * (-Icges(5,2) * t395 + t137 - t212);
t451 = qJD(2) * t463 + t186 * (-t176 * t258 + t127) - t187 * (-Icges(6,2) * t399 + t126 - t205);
t450 = t130 / 0.2e1;
t449 = t131 / 0.2e1;
t448 = t182 / 0.2e1;
t447 = t183 / 0.2e1;
t446 = -t186 / 0.2e1;
t445 = t186 / 0.2e1;
t444 = -t187 / 0.2e1;
t443 = t187 / 0.2e1;
t442 = t256 / 0.2e1;
t441 = -t258 / 0.2e1;
t440 = -qJD(2) / 0.2e1;
t439 = qJD(2) / 0.2e1;
t436 = rSges(5,1) * t257;
t433 = rSges(5,2) * t257;
t195 = rSges(5,1) * t255 + t433;
t161 = t195 * t258;
t245 = t256 * rSges(5,3);
t393 = t257 * t258;
t396 = t255 * t258;
t140 = rSges(5,1) * t393 - rSges(5,2) * t396 + t245;
t358 = qJD(4) * t256;
t58 = -t195 * t358 - t234 + (t140 + t382) * qJD(2);
t429 = t161 * t58;
t318 = -t195 * t357 + t233;
t369 = rSges(5,2) * t397 + rSges(5,3) * t258;
t139 = rSges(5,1) * t395 - t369;
t346 = -t139 + t383;
t57 = qJD(2) * t346 + t318;
t428 = t256 * t57;
t427 = t258 * t43;
t426 = qJDD(2) / 0.2e1;
t122 = Icges(6,5) * t399 - Icges(6,6) * t401 - t412;
t411 = t122 * t258;
t406 = t174 * t256;
t405 = t174 * t258;
t404 = t176 * t264;
t403 = t188 * t256;
t402 = t188 * t258;
t392 = t257 * qJD(4) ^ 2;
t62 = -t256 * t300 - t405;
t390 = t62 * qJD(2);
t76 = -t256 * t298 - t402;
t389 = t76 * qJD(2);
t388 = -t122 * t256 - t126 * t398;
t387 = t123 * t256 + t127 * t398;
t386 = -t133 * t256 - t137 * t393;
t385 = t134 * t256 + t138 * t393;
t159 = qJD(2) * t199 - t234;
t384 = t217 - (-t258 * t438 - t235) * qJD(2) - t159;
t109 = t143 + t199;
t375 = -t190 + t193;
t374 = t192 + t307;
t350 = t256 * t437;
t218 = t256 * t434;
t367 = rSges(4,3) * t258 + t218;
t142 = t350 - t367;
t372 = -t196 - t142;
t371 = rSges(5,3) * t359 + t360 * t464;
t368 = rSges(4,3) * t359 + qJD(2) * t218;
t355 = qJD(2) * qJD(3);
t366 = qJDD(3) * t256 + t258 * t355;
t225 = qJ(3) * t359;
t365 = t225 + t233;
t361 = qJD(2) * t189;
t356 = -m(4) - m(5) - m(6);
t352 = t128 * t359 + t256 * t79 + t258 * t78;
t343 = -pkin(2) - t437;
t341 = t360 / 0.2e1;
t340 = t359 / 0.2e1;
t339 = -t358 / 0.2e1;
t338 = t358 / 0.2e1;
t337 = -t357 / 0.2e1;
t336 = t357 / 0.2e1;
t286 = -pkin(4) * t255 - t180;
t335 = -t254 - t436;
t293 = t178 * t264;
t334 = qJD(2) * t127 - t124 * t264 - t256 * t293;
t333 = -t125 * t264 - t258 * t293 + (-t179 * t256 + t418) * qJD(2);
t332 = qJD(2) * t125 + t126 * t264 - t256 * t404;
t331 = t127 * t264 - t258 * t404 + (-t256 * t306 + t414) * qJD(2);
t100 = t127 * t399;
t330 = t123 * t258 - t100;
t105 = t138 * t395;
t329 = t134 * t258 - t105;
t328 = -t122 + t410;
t327 = -t133 + t407;
t326 = t463 * t264;
t325 = t179 * t264 - t404;
t165 = t459 * t264;
t319 = -t257 * t430 - t165;
t200 = rSges(3,1) * t258 - rSges(3,2) * t256;
t197 = rSges(3,1) * t256 + rSges(3,2) * t258;
t198 = t436 - t464;
t68 = t136 * t257 + t138 * t255;
t288 = qJD(4) * t190;
t84 = -t258 * t288 + (-t256 * t307 + t415) * qJD(2);
t289 = qJD(4) * t192;
t86 = -t258 * t289 + (-t193 * t256 + t419) * qJD(2);
t273 = -qJD(4) * t68 - t255 * t84 + t257 * t86 + t362;
t67 = t135 * t257 + t137 * t255;
t85 = qJD(2) * t136 - t256 * t288;
t87 = qJD(2) * t138 - t256 * t289;
t274 = qJD(2) * t133 - qJD(4) * t67 - t255 * t85 + t257 * t87;
t316 = -(t256 * t457 + t258 * t274) * t258 + (t256 * t458 + t258 * t273) * t256;
t315 = -(t256 * t274 - t258 * t457) * t258 + (t256 * t273 - t258 * t458) * t256;
t314 = -t256 * t43 - t258 * t42;
t54 = -t136 * t397 - t329;
t313 = t256 * t54 - t258 * t53;
t55 = -t135 * t396 - t386;
t56 = -t136 * t396 + t385;
t312 = t256 * t56 - t258 * t55;
t311 = -t256 * t58 - t258 * t57;
t88 = -t357 * t433 + (-t257 * t360 - t344) * rSges(5,1) + t371;
t160 = t195 * t256;
t89 = -qJD(4) * t160 + (t198 * t258 + t245) * qJD(2);
t310 = t256 * t89 + t258 * t88;
t60 = t124 * t253 + t126 * t252;
t301 = t139 * t256 + t140 * t258;
t299 = t190 * t257 + t192 * t255;
t297 = -t351 - t230;
t296 = -qJDD(3) * t258 + qJD(2) * (-pkin(2) * t360 + t365) + qJDD(2) * t199 + t256 * t355;
t290 = t305 * t256;
t285 = qJD(2) * t175 - t186 * t405 + t187 * t406;
t284 = qJD(2) * (qJD(2) * t294 - t225) + qJDD(2) * t121 + t296;
t282 = t135 * t258 - t136 * t256;
t281 = (-t255 * t374 + t257 * t375) * qJD(2);
t277 = qJD(2) * t122 - t252 * t332 + t253 * t334;
t11 = t256 * t455 + t258 * t277;
t276 = -t252 * t331 + t253 * t333 + t363;
t12 = t256 * t456 + t258 * t276;
t13 = t256 * t277 - t258 * t455;
t14 = t256 * t276 - t258 * t456;
t45 = -t290 - t411;
t46 = -t125 * t401 - t330;
t22 = t186 * t46 - t187 * t45 + t390;
t47 = -t124 * t400 - t388;
t48 = -t125 * t400 + t387;
t63 = -t258 * t300 + t406;
t59 = t63 * qJD(2);
t23 = t186 * t48 - t187 * t47 + t59;
t278 = (-t178 * t258 - t125) * t186 - (-t178 * t256 - t124) * t187 + (-t176 + t179) * qJD(2);
t270 = -t252 * t451 + t253 * t278;
t275 = qJD(2) * t174 - t252 * t326 + t253 * t325;
t30 = t256 * t454 + t258 * t275;
t31 = t256 * t275 - t258 * t454;
t32 = t252 * t334 + t253 * t332;
t33 = t252 * t333 + t253 * t331;
t61 = t125 * t253 + t127 * t252;
t280 = (qJD(2) * t30 + qJDD(2) * t63 - t11 * t187 + t12 * t186 + t130 * t48 + t131 * t47) * t442 + (t252 * t278 + t253 * t451) * t440 + (qJD(2) * t31 + qJDD(2) * t62 - t13 * t187 + t130 * t46 + t131 * t45 + t14 * t186) * t441 + t22 * t341 + t23 * t340 + (t256 * t48 - t258 * t47) * t450 + (t256 * t46 - t258 * t45) * t449 + (-t11 * t258 + t12 * t256 + (t256 * t47 + t258 * t48) * qJD(2)) * t445 + (-t13 * t258 + t14 * t256 + (t256 * t45 + t258 * t46) * qJD(2)) * t444 + (t256 * t61 - t258 * t60) * t426 + (t256 * t33 - t258 * t32 + (t256 * t60 + t258 * t61) * qJD(2)) * t439 + (t256 * t285 + t258 * t270) * t446 + (t256 * t270 - t258 * t285) * t443;
t169 = t307 * qJD(4);
t170 = t193 * qJD(4);
t272 = qJD(2) * t188 - qJD(4) * t299 - t169 * t255 + t170 * t257;
t271 = -t255 * t452 + t257 * t282;
t171 = t198 * qJD(4);
t95 = qJD(2) * t109 - t234;
t94 = qJD(2) * t372 + t233;
t77 = -t258 * t298 + t403;
t69 = t77 * qJD(2);
t64 = qJD(4) * t301 + qJD(1);
t52 = qJDD(2) * t143 + qJD(2) * (-qJD(2) * t350 + t368) + t296;
t51 = t372 * qJDD(2) + (-qJD(2) * t143 - t159) * qJD(2) + t366;
t39 = t256 * t272 - t258 * t453;
t38 = t256 * t453 + t258 * t272;
t36 = -qJD(4) * t302 + t255 * t86 + t257 * t84;
t35 = -qJD(4) * t303 + t255 * t87 + t257 * t85;
t34 = qJD(4) * t310 + t139 * t182 - t140 * t183 + qJDD(1);
t29 = qJD(2) * t88 + qJDD(2) * t140 - t171 * t358 - t182 * t195 + t284;
t28 = -t171 * t357 + t183 * t195 + t346 * qJDD(2) + (-t89 + t384) * qJD(2) + t366;
t25 = qJD(4) * t312 + t69;
t24 = qJD(4) * t313 + t389;
t10 = -t130 * t180 - t165 * t186 + t424 * qJDD(2) + (t78 + t80) * qJD(2) + (-t182 * t255 - t256 * t392) * pkin(4) + t284;
t9 = t131 * t180 - t165 * t187 + (t183 * t255 - t258 * t392) * pkin(4) + t320 * qJDD(2) + (-t79 - t81 + t384) * qJD(2) + t366;
t1 = [m(5) * t34 + m(6) * t8 + (m(4) - t461) * qJDD(1) + (t356 + t461) * g(3); (t69 + ((t54 - t105 + (t134 + t408) * t258 + t386) * t258 + t385 * t256) * qJD(4)) * t336 + (t59 + (t46 + (t124 * t258 + t125 * t256) * t252 + t330 + t388) * t187 + (-t126 * t399 + t411 + t45 + (t124 * t256 - t125 * t258) * t252 + t387) * t186) * t443 - m(3) * (-g(1) * t197 + g(2) * t200) + (t61 + t63) * t450 + (t60 + t62) * t449 + (t68 + t77) * t448 + (t67 + t76) * t447 + (-t390 + (t48 - t290 - t387) * t187 + (t328 * t256 - t100 + t47) * t186 + ((t123 + t305) * t186 + t328 * t187) * t258 + t22) * t446 + (t33 + t30) * t445 + (-t389 + ((t258 * t327 - t385 + t56) * t258 + (t256 * t327 + t329 + t55) * t256) * qJD(4) + t24) * t339 + (t36 + t38) * t338 + (-qJD(4) * t298 + t169 * t257 + t170 * t255 + t252 * t325 + t253 * t326) * qJD(2) + (t32 + t31 + t23) * t444 + (t35 + t39 + t25) * t337 + (t42 * (t186 * t435 + t256 * t348 - t370) + t43 * (t233 + t373) + (-t295 - t347) * t427 + ((t42 * (t297 + t432) - t43 * t261) * t258 + (t42 * (t261 - rSges(6,3)) + t43 * t297) * t256) * qJD(2) - (qJD(2) * t425 + t283 - t42 + t460) * t43 + (t10 - g(2)) * (t129 + t167 - t394) + (t9 - g(1)) * (-t128 + t377)) * m(6) + (-(-qJD(2) * t139 + t318 + t460 - t57) * t58 + t57 * (t217 + t234) + t58 * (t233 + t371) + (t195 * t428 - t429) * qJD(4) + ((-rSges(5,3) * t57 + t335 * t58) * t256 + (t57 * (-t198 - t254) - t58 * t267) * t258) * qJD(2) + (t29 - g(2)) * (t140 + t324) + (t28 - g(1)) * (t256 * t335 + t369 - t391)) * m(5) + (-(-qJD(2) * t142 - t185 + t233 - t94) * t95 + t94 * t234 + t95 * (t365 + t368) + (t94 * (t343 + t434) * t258 + (t94 * (-rSges(4,3) - qJ(3)) + t95 * t343) * t256) * qJD(2) + (t52 - g(2)) * t109 + (t51 - g(1)) * (t256 * t343 + t236 + t367)) * m(4) + (t299 + t176 * t253 + t178 * t252 + Icges(4,2) * t266 ^ 2 + (Icges(4,1) * t265 + 0.2e1 * Icges(4,4) * t266) * t265 + m(3) * (t197 ^ 2 + t200 ^ 2) + Icges(3,3)) * qJDD(2); t356 * (g(1) * t256 - g(2) * t258) + 0.2e1 * (t10 * t441 + t442 * t9) * m(6) + 0.2e1 * (t28 * t442 + t29 * t441) * m(5) + 0.2e1 * (t441 * t52 + t442 * t51) * m(4); ((t255 * t375 + t257 * t374) * qJD(2) + (t282 * t255 + t257 * t452) * qJD(4)) * t440 + (t256 * t36 - t258 * t35 + (t67 * t256 + t258 * t68) * qJD(2)) * t439 + ((t256 * t55 + t258 * t56) * qJD(2) + t316) * t338 + (t256 * t68 - t258 * t67) * t426 + t312 * t448 + t313 * t447 + ((t256 * t53 + t258 * t54) * qJD(2) + t315) * t337 + t280 + (qJD(2) * t39 + qJD(4) * t315 + qJDD(2) * t76 + t182 * t54 + t183 * t53) * t441 + (qJD(2) * t38 + qJD(4) * t316 + qJDD(2) * t77 + t182 * t56 + t183 * t55) * t442 + ((-t357 * t403 - t361) * t258 + (t281 + (t258 * t402 + t271) * qJD(4)) * t256) * t336 + ((-t358 * t402 + t361) * t256 + (t281 + (t256 * t403 + t271) * qJD(4)) * t258) * t339 + t24 * t341 + t25 * t340 + (t37 * t352 + (t9 * t286 + t42 * t319 + t8 * t99 + t37 * t80 + (t286 * t43 - t37 * t98) * qJD(2)) * t258 + (t10 * t286 + t43 * t319 - t8 * t98 + t37 * t81 + (t180 * t42 - t37 * t424) * qJD(2)) * t256 - g(3) * (t459 + t251) - (g(1) * t258 + g(2) * t256) * t286 - (-t43 * t255 * t359 + (t314 * t257 + t37 * (-t256 ^ 2 - t258 ^ 2) * t255) * qJD(4)) * pkin(4) + t462) * m(6) + (-(t160 * t57 - t429) * qJD(2) - (t64 * (-t160 * t256 - t161 * t258) + t311 * t198) * qJD(4) + t34 * t301 + t64 * ((t139 * t258 - t140 * t256) * qJD(2) + t310) + t311 * t171 + (-t29 * t256 - t28 * t258 + (-t258 * t58 + t428) * qJD(2)) * t195 + g(1) * t161 + g(2) * t160 - g(3) * t198) * m(5); t280 + (t37 * (-t129 * t360 + t352) + t314 * t165 + (-t10 * t256 - t9 * t258 + (t256 * t42 - t427) * qJD(2)) * t180 + g(1) * t151 + g(2) * t150 - g(3) * t459 + t462) * m(6);];
tau = t1;
