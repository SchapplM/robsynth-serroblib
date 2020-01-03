% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR10_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:59
% EndTime: 2019-12-31 18:04:14
% DurationCPUTime: 10.89s
% Computational Cost: add. (11275->639), mult. (19396->832), div. (0->0), fcn. (19066->8), ass. (0->304)
t299 = cos(qJ(1));
t464 = rSges(6,3) * t299;
t408 = qJ(4) + qJ(5);
t285 = sin(t408);
t295 = cos(pkin(8));
t294 = sin(pkin(8));
t361 = cos(t408);
t341 = t294 * t361;
t230 = -t295 * t285 + t341;
t313 = t294 * t285 + t295 * t361;
t183 = rSges(6,1) * t230 - rSges(6,2) * t313;
t296 = sin(qJ(4));
t411 = t295 * t296;
t372 = pkin(4) * t411;
t298 = cos(qJ(4));
t278 = pkin(4) * t298 + pkin(3);
t428 = pkin(3) - t278;
t211 = -t294 * t428 - t372;
t376 = qJD(4) * t299;
t198 = t211 * t376;
t293 = qJD(4) + qJD(5);
t254 = t293 * t299;
t463 = t183 * t254 + t198;
t297 = sin(qJ(1));
t410 = t295 * t297;
t413 = t294 * t298;
t223 = t296 * t410 - t297 * t413;
t414 = t294 * t296;
t324 = t295 * t298 + t414;
t224 = t324 * t297;
t144 = Icges(5,5) * t224 - Icges(5,6) * t223 + Icges(5,3) * t299;
t213 = Icges(5,4) * t224;
t147 = -Icges(5,2) * t223 + Icges(5,6) * t299 + t213;
t212 = Icges(5,4) * t223;
t151 = -Icges(5,1) * t224 - Icges(5,5) * t299 + t212;
t242 = -t411 + t413;
t225 = t242 * t299;
t226 = t324 * t299;
t407 = t225 * t147 - t151 * t226;
t51 = -t144 * t297 + t407;
t415 = t293 * t297;
t207 = t285 * t410 - t297 * t341;
t208 = t313 * t297;
t115 = Icges(6,5) * t208 - Icges(6,6) * t207 + Icges(6,3) * t299;
t201 = Icges(6,4) * t208;
t118 = -Icges(6,2) * t207 + Icges(6,6) * t299 + t201;
t200 = Icges(6,4) * t207;
t122 = -Icges(6,1) * t208 - Icges(6,5) * t299 + t200;
t209 = t230 * t299;
t210 = t313 * t299;
t47 = -t115 * t297 + t209 * t118 - t122 * t210;
t117 = Icges(6,5) * t210 + Icges(6,6) * t209 - Icges(6,3) * t297;
t419 = Icges(6,4) * t210;
t120 = Icges(6,2) * t209 - Icges(6,6) * t297 + t419;
t202 = Icges(6,4) * t209;
t123 = Icges(6,1) * t210 - Icges(6,5) * t297 + t202;
t48 = -t117 * t297 + t209 * t120 + t210 * t123;
t177 = Icges(6,5) * t230 - Icges(6,6) * t313;
t418 = Icges(6,4) * t230;
t179 = -Icges(6,2) * t313 + t418;
t222 = Icges(6,4) * t313;
t181 = Icges(6,1) * t230 - t222;
t55 = -t177 * t297 + t179 * t209 + t181 * t210;
t17 = t55 * qJD(1) + t254 * t47 - t415 * t48;
t59 = -t118 * t313 - t122 * t230;
t462 = -t147 * t223 - t151 * t224;
t461 = -t147 * t324 - t151 * t242;
t45 = t115 * t299 - t118 * t207 - t122 * t208;
t146 = Icges(5,5) * t226 + Icges(5,6) * t225 - Icges(5,3) * t297;
t421 = Icges(5,4) * t226;
t149 = Icges(5,2) * t225 - Icges(5,6) * t297 + t421;
t214 = Icges(5,4) * t225;
t152 = Icges(5,1) * t226 - Icges(5,5) * t297 + t214;
t406 = t225 * t149 + t226 * t152;
t342 = t146 * t297 - t406;
t49 = t144 * t299 + t462;
t460 = t342 - t49;
t124 = rSges(6,1) * t208 - rSges(6,2) * t207 + t464;
t300 = -pkin(7) - pkin(6);
t374 = pkin(3) * t410;
t373 = pkin(4) * t414;
t387 = -t278 * t410 - t297 * t373;
t172 = t374 + (pkin(6) + t300) * t299 + t387;
t402 = -t124 + t172;
t455 = (pkin(4) * t413 - t372) * qJD(4);
t454 = 0.2e1 * qJD(4);
t453 = -t223 * t149 + t224 * t152;
t379 = qJD(3) * t294;
t451 = (qJD(4) * t211 + t379) * t297;
t196 = rSges(5,1) * t242 - rSges(5,2) * t324;
t450 = (qJD(4) * t196 + t379) * t297;
t391 = t210 * rSges(6,1) + t209 * rSges(6,2);
t126 = -rSges(6,3) * t297 + t391;
t142 = -rSges(6,1) * t207 - rSges(6,2) * t208;
t143 = rSges(6,1) * t209 - rSges(6,2) * t210;
t381 = qJD(1) * t297;
t449 = t126 * t381 + t142 * t415 + t143 * t254;
t259 = t299 * pkin(1) + t297 * qJ(2);
t409 = t295 * t299;
t412 = t294 * t299;
t321 = rSges(3,1) * t409 - rSges(3,2) * t412 + t297 * rSges(3,3);
t448 = t259 + t321;
t346 = pkin(2) * t409 + qJ(3) * t412 + t259;
t364 = rSges(4,1) * t409 + t297 * rSges(4,2) + rSges(4,3) * t412;
t447 = t346 + t364;
t445 = qJD(4) * t242;
t46 = t299 * t117 - t207 * t120 + t208 * t123;
t444 = t295 * t428 - t373;
t443 = t297 * (-Icges(5,2) * t226 + t152 + t214) - t299 * (-Icges(5,2) * t224 - t151 - t212);
t304 = qJD(1) * (-Icges(6,2) * t230 + t181 - t222) - t415 * (-Icges(6,2) * t210 + t123 + t202) + t254 * (-Icges(6,2) * t208 - t122 - t200);
t349 = qJD(1) * t293;
t243 = t297 * t349;
t442 = -t243 / 0.2e1;
t244 = t299 * t349;
t441 = -t244 / 0.2e1;
t439 = -t415 / 0.2e1;
t438 = -t254 / 0.2e1;
t437 = t254 / 0.2e1;
t436 = -t297 / 0.2e1;
t435 = t297 / 0.2e1;
t434 = -t299 / 0.2e1;
t433 = t299 / 0.2e1;
t432 = -rSges(5,3) - pkin(6);
t431 = pkin(6) * t299;
t430 = -qJD(1) / 0.2e1;
t429 = qJD(1) / 0.2e1;
t380 = qJD(1) * t299;
t203 = t313 * t293;
t111 = -t203 * t299 - t230 * t381;
t112 = t209 * t293 - t313 * t381;
t405 = t112 * rSges(6,1) + t111 * rSges(6,2);
t74 = -rSges(6,3) * t380 + t405;
t315 = t299 * t455 + t300 * t380;
t94 = (t297 * t444 + t431) * qJD(1) + t315;
t427 = -t74 - t94;
t113 = qJD(1) * t209 - t203 * t297;
t114 = qJD(1) * t210 + t230 * t415;
t337 = -rSges(6,1) * t114 - rSges(6,2) * t113;
t75 = -rSges(6,3) * t381 - t337;
t281 = t297 * t300;
t282 = pkin(6) * t381;
t314 = t297 * t445;
t95 = t282 + pkin(4) * t314 + (-t299 * t444 + t281) * qJD(1);
t426 = -t75 - t95;
t425 = rSges(3,2) * t294;
t422 = -rSges(4,3) - qJ(3);
t420 = Icges(5,4) * t242;
t190 = Icges(5,5) * t242 - Icges(5,6) * t324;
t192 = -Icges(5,2) * t324 + t420;
t237 = Icges(5,4) * t324;
t194 = Icges(5,1) * t242 - t237;
t64 = t190 * t299 - t192 * t223 + t194 * t224;
t416 = qJD(1) * t64;
t271 = pkin(3) * t409;
t247 = -pkin(6) * t297 + t271;
t366 = t278 * t409 + t299 * t373 + t281;
t173 = -t247 + t366;
t401 = -t126 - t173;
t204 = t230 * t293;
t130 = -rSges(6,1) * t203 - rSges(6,2) * t204;
t233 = t324 * qJD(4);
t227 = pkin(4) * t233;
t400 = t130 - t227;
t160 = -t233 * t299 - t242 * t381;
t161 = t299 * t445 - t324 * t381;
t397 = t161 * rSges(5,1) + t160 * rSges(5,2);
t394 = -Icges(5,2) * t242 + t194 - t237;
t393 = -Icges(5,1) * t324 - t192 - t420;
t284 = qJD(2) * t299;
t231 = qJD(1) * t259 - t284;
t336 = pkin(2) * t295 + qJ(3) * t294;
t363 = t297 * t379;
t392 = -t336 * t380 - t231 - t363;
t390 = t226 * rSges(5,1) + t225 * rSges(5,2);
t375 = qJD(1) * qJD(2);
t283 = qJD(2) * t297;
t383 = qJ(2) * t380 + t283;
t389 = qJD(1) * (-pkin(1) * t381 + t383) + t297 * t375;
t235 = t336 * t297;
t287 = t299 * qJ(2);
t257 = pkin(1) * t297 - t287;
t388 = -t235 - t257;
t267 = t297 * t425;
t386 = rSges(3,3) * t380 + qJD(1) * t267;
t265 = t299 * t379;
t385 = t265 + t283;
t384 = t299 * rSges(3,3) + t267;
t382 = -qJD(1) * t257 + t283;
t378 = qJD(3) * t295;
t377 = qJD(4) * t297;
t50 = t299 * t146 + t453;
t371 = rSges(3,1) * t410;
t369 = -qJD(1) * t271 + t282 + t392;
t246 = t374 + t431;
t368 = -t246 + t388;
t367 = t247 + t346;
t365 = t265 + t383;
t188 = t196 * t376;
t362 = -rSges(3,1) * t295 - pkin(1);
t360 = -t381 / 0.2e1;
t359 = -t380 / 0.2e1;
t357 = t377 / 0.2e1;
t356 = -t376 / 0.2e1;
t182 = -rSges(6,1) * t313 - rSges(6,2) * t230;
t352 = qJD(1) * t143 + t182 * t415;
t348 = t389 + (-t336 * t381 + 0.2e1 * t265) * qJD(1);
t347 = -qJD(1) * t235 + t265 + t382;
t343 = -pkin(1) + (-rSges(4,1) - pkin(2)) * t295;
t339 = rSges(4,1) * t295 + rSges(4,3) * t294;
t162 = qJD(1) * t225 - t233 * t297;
t163 = qJD(1) * t226 + t314;
t338 = rSges(5,1) * t163 + rSges(5,2) * t162;
t332 = -qJD(1) ^ 2 * t246 + t348;
t28 = -t227 * t377 + t130 * t415 + t183 * t244 + (t198 - t427) * qJD(1) + t332;
t276 = t299 * t375;
t29 = -t227 * t376 + t130 * t254 - t183 * t243 + t276 + (t369 + t426 - t451) * qJD(1);
t335 = t28 * t297 + t29 * t299;
t187 = -rSges(5,1) * t233 - rSges(5,2) * t445;
t83 = -rSges(5,3) * t380 + t397;
t37 = t187 * t377 + (t83 + t188) * qJD(1) + t332;
t84 = -rSges(5,3) * t381 + t338;
t38 = t187 * t376 + t276 + (t369 - t84 - t450) * qJD(1);
t334 = t37 * t297 + t38 * t299;
t153 = rSges(5,1) * t224 - rSges(5,2) * t223 + rSges(5,3) * t299;
t61 = t188 + (-t153 + t368) * qJD(1) + t385;
t155 = -rSges(5,3) * t297 + t390;
t62 = -t284 + t450 + (t155 + t367) * qJD(1);
t333 = t297 * t62 + t299 * t61;
t331 = -qJD(1) * t246 + t347;
t326 = -t153 * t297 - t155 * t299;
t325 = (-Icges(5,5) * t223 - Icges(5,6) * t224) * t299 - (Icges(5,5) * t225 - Icges(5,6) * t226) * t297;
t320 = -pkin(1) - t336;
t318 = (-t297 * t50 + t299 * t49) * qJD(4);
t317 = (t297 * t342 + t299 * t51) * qJD(4);
t54 = t177 * t299 - t179 * t207 + t181 * t208;
t312 = qJD(1) * (-Icges(6,5) * t313 - Icges(6,6) * t230) + (-Icges(6,5) * t207 - Icges(6,6) * t208) * t254 - (Icges(6,5) * t209 - Icges(6,6) * t210) * t415;
t311 = -pkin(3) * t295 + t320;
t310 = -(Icges(5,1) * t225 - t149 - t421) * t297 + (-Icges(5,1) * t223 - t147 - t213) * t299;
t68 = Icges(6,5) * t112 + Icges(6,6) * t111 - Icges(6,3) * t380;
t70 = Icges(6,4) * t112 + Icges(6,2) * t111 - Icges(6,6) * t380;
t72 = Icges(6,1) * t112 + Icges(6,4) * t111 - Icges(6,5) * t380;
t10 = t111 * t120 + t112 * t123 - t117 * t380 + t209 * t70 + t210 * t72 - t297 * t68;
t69 = Icges(6,5) * t114 + Icges(6,6) * t113 - Icges(6,3) * t381;
t71 = Icges(6,4) * t114 + Icges(6,2) * t113 - Icges(6,6) * t381;
t73 = Icges(6,1) * t114 + Icges(6,4) * t113 - Icges(6,5) * t381;
t11 = t113 * t118 - t114 * t122 - t115 * t381 - t207 * t71 + t208 * t73 + t299 * t69;
t12 = t113 * t120 + t114 * t123 - t117 * t381 - t207 * t70 + t208 * t72 + t299 * t68;
t127 = -Icges(6,5) * t203 - Icges(6,6) * t204;
t128 = -Icges(6,4) * t203 - Icges(6,2) * t204;
t129 = -Icges(6,1) * t203 - Icges(6,4) * t204;
t26 = t111 * t179 + t112 * t181 - t127 * t297 + t128 * t209 + t129 * t210 - t177 * t380;
t27 = t113 * t179 + t114 * t181 + t127 * t299 - t128 * t207 + t129 * t208 - t177 * t381;
t30 = -t118 * t204 + t122 * t203 + t230 * t73 - t313 * t71;
t305 = -(Icges(6,1) * t209 - t120 - t419) * t415 + (-Icges(6,1) * t207 - t118 - t201) * t254 + (-Icges(6,1) * t313 - t179 - t418) * qJD(1);
t31 = -t120 * t204 - t123 * t203 + t230 * t72 - t313 * t70;
t60 = -t120 * t313 + t123 * t230;
t9 = t111 * t118 - t112 * t122 - t115 * t380 + t209 * t71 + t210 * t73 - t297 * t69;
t308 = (qJD(1) * t26 - t10 * t415 - t243 * t47 - t244 * t48 + t254 * t9) * t436 + (qJD(1) * t54 + t254 * t45 - t415 * t46) * t360 + t17 * t359 + (qJD(1) * t27 + t11 * t254 - t12 * t415 - t243 * t45 - t244 * t46) * t433 + (-t297 * t46 + t299 * t45) * t442 + (-t297 * t48 + t299 * t47) * t441 + (-t10 * t297 + t299 * t9 + (-t297 * t47 - t299 * t48) * qJD(1)) * t439 + (t11 * t299 - t12 * t297 + (-t297 * t45 - t299 * t46) * qJD(1)) * t437 + (-t297 * t31 + t299 * t30 + (-t297 * t59 - t299 * t60) * qJD(1)) * t429 + (t209 * t304 + t210 * t305 - t297 * t312) * t415 / 0.2e1 + (-t207 * t304 + t305 * t208 + t312 * t299) * t438 + (t230 * t305 - t304 * t313) * t430;
t307 = -t278 * t295 + t320 - t373;
t306 = -t297 * t84 - t299 * t83 + (-t153 * t299 + t155 * t297) * qJD(1);
t291 = t299 * rSges(4,2);
t280 = rSges(4,2) * t380;
t238 = t324 * pkin(4);
t220 = t371 - t384;
t219 = t297 * t339 - t291;
t216 = pkin(4) * t225;
t215 = t242 * t297 * pkin(4);
t189 = -Icges(5,5) * t324 - Icges(5,6) * t242;
t186 = -Icges(5,1) * t233 - Icges(5,4) * t445;
t185 = -Icges(5,4) * t233 - Icges(5,2) * t445;
t184 = -Icges(5,5) * t233 - Icges(5,6) * t445;
t175 = qJD(1) * t448 - t284;
t174 = t283 + (-t220 - t257) * qJD(1);
t171 = rSges(5,1) * t225 - rSges(5,2) * t226;
t170 = -rSges(5,1) * t223 - rSges(5,2) * t224;
t141 = t254 * t182;
t132 = t276 + (-qJD(1) * t321 - t231) * qJD(1);
t131 = qJD(1) * (-qJD(1) * t371 + t386) + t389;
t105 = qJD(1) * t447 - t284 + t363;
t104 = (-t219 + t388) * qJD(1) + t385;
t86 = t276 + (-qJD(1) * t364 - t363 + t392) * qJD(1);
t85 = (-t339 * t381 + t280) * qJD(1) + t348;
t82 = Icges(5,1) * t163 + Icges(5,4) * t162 - Icges(5,5) * t381;
t81 = Icges(5,1) * t161 + Icges(5,4) * t160 - Icges(5,5) * t380;
t80 = Icges(5,4) * t163 + Icges(5,2) * t162 - Icges(5,6) * t381;
t79 = Icges(5,4) * t161 + Icges(5,2) * t160 - Icges(5,6) * t380;
t78 = Icges(5,5) * t163 + Icges(5,6) * t162 - Icges(5,3) * t381;
t77 = Icges(5,5) * t161 + Icges(5,6) * t160 - Icges(5,3) * t380;
t76 = qJD(4) * t326 - t378;
t67 = -t149 * t324 + t152 * t242;
t65 = -t190 * t297 + t192 * t225 + t194 * t226;
t63 = t65 * qJD(1);
t44 = -t378 - t124 * t415 - t126 * t254 + (t172 * t297 - t173 * t299) * qJD(4);
t43 = t183 * t415 - t284 + t451 + (t367 - t401) * qJD(1);
t42 = (t368 + t402) * qJD(1) + t385 + t463;
t36 = t306 * qJD(4);
t35 = t162 * t192 + t163 * t194 + t184 * t299 - t185 * t223 + t186 * t224 - t190 * t381;
t34 = t160 * t192 + t161 * t194 - t184 * t297 + t185 * t225 + t186 * t226 - t190 * t380;
t33 = -t149 * t445 - t152 * t233 + t242 * t81 - t324 * t79;
t32 = -t147 * t445 + t151 * t233 + t242 * t82 - t324 * t80;
t23 = t63 + t317;
t22 = t318 + t416;
t13 = -t124 * t244 + t126 * t243 - t415 * t75 - t254 * t74 + (-t297 * t95 - t299 * t94 + (t172 * t299 + t173 * t297) * qJD(1)) * qJD(4);
t1 = [(t63 + (t407 * t299 + (t460 + t462) * t297) * qJD(4)) * t356 + t17 * t438 + (t59 + t54) * t442 + (t60 + t55) * t441 + (t31 + t26) * t439 - (t33 + t34) * t377 / 0.2e1 + (t22 - t416 + ((t406 + t460) * t299 + t453 * t297) * qJD(4)) * t357 + (-t128 * t313 + t129 * t230 - t179 * t204 - t181 * t203 - t185 * t324 + t186 * t242 - t192 * t445 - t194 * t233) * qJD(1) + (t29 * (t299 * t300 - t124 + t287 + t387) + t28 * (t346 + t366 + t391) + (-t28 * rSges(6,3) + t29 * t320) * t297 + (t284 + t337 + (-t379 - t455) * t297 + (t307 * t299 + (rSges(6,3) - qJ(2) - t300) * t297) * qJD(1)) * t42 + (t315 + t365 + t405 - t331 + t42 + (t307 * t297 - t402 - t464) * qJD(1) - t463) * t43) * m(6) + (t38 * (-t153 + t287 - t431) + t61 * (t282 + t284 - t338) + t37 * (t271 + t346 + t390) + t62 * (t365 + t397) + (t38 * t311 + t37 * t432 - t61 * t379) * t297 + ((t311 * t61 + t432 * t62) * t299 + (t61 * (rSges(5,3) - qJ(2)) + t62 * t311) * t297) * qJD(1) - (-qJD(1) * t153 + t188 + t331 - t61) * t62) * m(5) + (t86 * (t287 + t291) + t104 * t284 + t85 * t447 + t105 * (t280 + t365) + (t86 * t343 + (-t104 * qJD(3) + t422 * t86) * t294) * t297 + (t104 * (t294 * t422 + t343) * t299 + (t104 * (-rSges(4,2) - qJ(2)) + t105 * (t320 - t339)) * t297) * qJD(1) - (-qJD(1) * t219 - t104 + t347) * t105) * m(4) + (t132 * (t297 * t362 + t287 + t384) + t174 * t284 + t131 * t448 + t175 * (t383 + t386) + (t174 * (t362 + t425) * t299 + (t174 * (-rSges(3,3) - qJ(2)) + t175 * t362) * t297) * qJD(1) - (-qJD(1) * t220 - t174 + t382) * t175) * m(3) + (t27 + t17 + t30) * t437 + (t32 + t35 + t23) * t376 / 0.2e1 + ((t67 + t65) * t299 + (t461 + t64) * t297) * qJD(4) * t430; 0.2e1 * (t28 * t434 + t29 * t435) * m(6) + 0.2e1 * (t37 * t434 + t38 * t435) * m(5) + 0.2e1 * (t434 * t85 + t435 * t86) * m(4) + 0.2e1 * (t131 * t434 + t132 * t435) * m(3); 0.2e1 * (-m(5) * t36 / 0.2e1 - m(6) * t13 / 0.2e1) * t295 + 0.2e1 * (m(4) * (t297 * t85 + t299 * t86) / 0.2e1 + m(5) * t334 / 0.2e1 + m(6) * t335 / 0.2e1) * t294; (-t297 * t33 + t299 * t32 + (-t297 * t461 - t67 * t299) * qJD(1)) * t429 + t308 + ((t242 * t393 - t324 * t394) * qJD(1) + (t242 * t310 + t324 * t443) * qJD(4)) * t430 + ((-t297 * t189 + t225 * t394 + t226 * t393) * qJD(1) + (-t225 * t443 + t226 * t310 - t297 * t325) * qJD(4)) * t357 + ((t299 * t189 - t223 * t394 + t224 * t393) * qJD(1) + (t223 * t443 + t310 * t224 + t325 * t299) * qJD(4)) * t356 + (qJD(1) * t34 + ((-t144 * t380 + t147 * t160 - t151 * t161 + t225 * t80 + t226 * t82 - t297 * t78) * t299 - (-t146 * t380 + t149 * t160 + t152 * t161 + t225 * t79 + t226 * t81 - t297 * t77) * t297 + (-t51 * t297 + t299 * t342) * qJD(1)) * t454) * t436 + (qJD(1) * t35 + ((-t144 * t381 + t147 * t162 - t151 * t163 - t223 * t80 + t224 * t82 + t299 * t78) * t299 - (-t146 * t381 + t149 * t162 + t152 * t163 - t223 * t79 + t224 * t81 + t299 * t77) * t297 + (-t49 * t297 - t50 * t299) * qJD(1)) * t454) * t433 + (t318 + t22) * t360 + (t317 + t23) * t359 + (-t42 * (-t238 * t376 + t141 + (-t142 - t215) * qJD(1)) - t43 * (qJD(1) * t216 - t238 * t377 + t352) + (t13 * t402 + t43 * t400) * t297 + (t13 * t401 + t42 * t400) * t299 + ((-qJD(1) * t42 + t28) * t297 + (qJD(1) * t43 + t29) * t299) * (t183 + t211) + (-(-t215 * t297 - t216 * t299) * qJD(4) + (qJD(1) * t173 + t426) * t297 + (qJD(1) * t402 + t427) * t299 + t449) * t44) * m(6) + (-(-t170 * t61 + t171 * t62) * qJD(1) - (t76 * (-t170 * t297 - t171 * t299) + t333 * (-rSges(5,1) * t324 - rSges(5,2) * t242)) * qJD(4) + t36 * t326 + t76 * t306 + t333 * t187 + ((-t297 * t61 + t299 * t62) * qJD(1) + t334) * t196) * m(5); t308 + (-t42 * (-qJD(1) * t142 + t141) - t43 * t352 + t13 * (-t124 * t297 - t126 * t299) + (t297 * t43 + t299 * t42) * t130 + ((-t297 * t42 + t299 * t43) * qJD(1) + t335) * t183 + (-t297 * t75 + (-qJD(1) * t124 - t74) * t299 + t449) * t44) * m(6);];
tauc = t1(:);
