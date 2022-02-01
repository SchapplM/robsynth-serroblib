% Calculate vector of inverse dynamics joint torques for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:04
% EndTime: 2022-01-23 08:59:20
% DurationCPUTime: 10.90s
% Computational Cost: add. (11526->619), mult. (27090->785), div. (0->0), fcn. (30800->10), ass. (0->300)
t288 = sin(pkin(8));
t291 = cos(pkin(8));
t294 = sin(qJ(1));
t292 = cos(pkin(7));
t296 = cos(qJ(1));
t385 = t292 * t296;
t231 = t288 * t294 + t291 * t385;
t287 = sin(pkin(9));
t290 = cos(pkin(9));
t289 = sin(pkin(7));
t389 = t289 * t296;
t176 = t231 * t287 - t290 * t389;
t388 = t290 * t292;
t225 = t287 * t289 + t291 * t388;
t383 = t296 * t288;
t169 = t225 * t294 - t290 * t383;
t386 = t292 * t294;
t228 = t288 * t386 + t291 * t296;
t293 = sin(qJ(5));
t295 = cos(qJ(5));
t128 = -t169 * t293 + t228 * t295;
t229 = t291 * t386 - t383;
t387 = t290 * t294;
t173 = t229 * t287 - t289 * t387;
t393 = t288 * t293;
t170 = t225 * t295 + t292 * t393;
t392 = t288 * t295;
t227 = t290 * t392 - t291 * t293;
t435 = -t170 * t294 + t296 * t227;
t456 = t173 * (Icges(6,5) * t128 + Icges(6,6) * t435);
t131 = t170 * t296 + t227 * t294;
t171 = -t225 * t293 + t292 * t392;
t226 = t290 * t393 + t291 * t295;
t429 = t171 * t296 - t294 * t226;
t77 = Icges(6,5) * t429 - Icges(6,6) * t131;
t461 = t176 * t77 + t456;
t391 = t289 * t291;
t224 = -t292 * t287 + t290 * t391;
t324 = -t224 * t293 + t289 * t392;
t166 = t224 * t295 + t289 * t393;
t446 = Icges(6,4) * t166;
t110 = Icges(6,1) * t324 - t446;
t223 = t287 * t391 + t388;
t97 = Icges(6,2) * t324 + Icges(6,6) * t223 + t446;
t460 = t110 - t97;
t129 = t169 * t295 + t228 * t293;
t133 = t171 * t294 + t226 * t296;
t437 = Icges(6,2) * t133 + Icges(6,6) * t173;
t57 = Icges(6,4) * t129 + t437;
t381 = Icges(6,2) * t429 + Icges(6,6) * t176;
t458 = Icges(6,4) * t131;
t59 = t381 + t458;
t81 = Icges(6,1) * t429 - t458;
t454 = Icges(6,4) * t435;
t82 = Icges(6,1) * t128 + t454;
t459 = (-t59 + t81) * t176 + (-t57 + t82) * t173;
t247 = pkin(4) * t290 + pkin(6) * t287 + pkin(3);
t417 = t288 * qJ(4) + pkin(2);
t141 = (t247 * t291 + t417) * t292 + pkin(1) + (pkin(4) * t287 - pkin(6) * t290 + qJ(3)) * t289;
t137 = t141 * t296;
t327 = qJ(4) * t291 - qJ(2);
t190 = -t247 * t288 + t327;
t444 = t190 * t294;
t308 = t444 - t137;
t199 = rSges(3,1) * t385 - rSges(3,2) * t389 + t294 * rSges(3,3);
t282 = t294 * qJ(2);
t286 = t296 * pkin(1);
t252 = t286 + t282;
t244 = qJD(1) * t252;
t280 = qJD(2) * t296;
t222 = -t280 + t244;
t452 = -qJD(1) * t199 - t222;
t379 = rSges(6,2) * t429 + t176 * rSges(6,3);
t451 = t379 - t308;
t418 = t289 * qJ(3) + pkin(1);
t189 = (pkin(3) * t291 + t417) * t292 + t418;
t181 = t189 * t294;
t245 = -t288 * pkin(3) + t327;
t325 = t245 * t296 + t181;
t246 = pkin(2) * t292 + t418;
t236 = t246 * t294;
t283 = t296 * qJ(2);
t361 = t283 - t236;
t450 = -t361 - t325;
t436 = Icges(6,6) * t133 + Icges(6,3) * t173;
t54 = Icges(6,5) * t435 - t436;
t346 = t292 * t383;
t230 = -t291 * t294 + t346;
t300 = (t225 * t296 + t288 * t387) * t295;
t132 = t230 * t293 + t300;
t382 = Icges(6,6) * t429 + Icges(6,3) * t176;
t56 = Icges(6,5) * t132 + t382;
t449 = t173 * t56 + t176 * t54;
t433 = Icges(6,4) * t324;
t109 = -Icges(6,2) * t166 + t433;
t192 = qJD(5) * t223 + qJD(1);
t438 = Icges(6,4) * t133 + Icges(6,5) * t173;
t61 = Icges(6,1) * t129 + t438;
t445 = Icges(6,4) * t429;
t380 = Icges(6,5) * t176 + t445;
t63 = Icges(6,1) * t131 + t380;
t79 = -Icges(6,2) * t131 + t445;
t80 = Icges(6,4) * t128 + Icges(6,2) * t435;
t98 = Icges(6,1) * t166 + Icges(6,5) * t223 + t433;
t448 = (t109 + t98) * t192 + ((t63 + t79) * t176 + (t61 + t80) * t173) * qJD(5);
t354 = qJD(5) * t173;
t99 = rSges(6,1) * t166 + rSges(6,2) * t324 + rSges(6,3) * t223;
t447 = t354 * t99;
t182 = t189 * t296;
t234 = t245 * t294;
t443 = t182 - t234;
t390 = t289 * t294;
t343 = rSges(5,2) * t173 - t228 * rSges(5,3) - rSges(5,1) * (t229 * t290 + t287 * t390);
t442 = t133 * rSges(6,2) + rSges(6,3) * t173;
t96 = Icges(6,5) * t166 + Icges(6,6) * t324 + Icges(6,3) * t223;
t441 = t133 * t97 + t173 * t96;
t55 = Icges(6,5) * t131 + t382;
t440 = t133 * t59 + t173 * t55;
t53 = Icges(6,5) * t129 + t436;
t439 = t133 * t57 + t53 * t173;
t402 = t141 * t294;
t326 = -t190 * t296 - t402;
t237 = t246 * t296;
t370 = t237 + t282;
t430 = -t370 + t443;
t203 = -t286 + t237;
t357 = qJD(3) * t294;
t263 = t289 * t357;
t362 = t280 - t263;
t319 = -qJD(1) * t203 - t244 + t362;
t309 = -qJD(4) * t228 + t319;
t138 = t229 * rSges(4,1) - t228 * rSges(4,2) + rSges(4,3) * t390;
t427 = t292 ^ 2;
t426 = -m(5) - m(6);
t208 = t231 * qJD(1);
t359 = qJD(1) * t296;
t334 = t289 * t359;
t150 = t208 * t287 - t290 * t334;
t106 = qJD(5) * t150 + qJDD(5) * t173;
t425 = t106 / 0.2e1;
t206 = t229 * qJD(1);
t360 = qJD(1) * t294;
t335 = t289 * t360;
t149 = -t206 * t287 + t290 * t335;
t107 = qJD(5) * t149 + qJDD(5) * t176;
t424 = t107 / 0.2e1;
t423 = t294 / 0.2e1;
t422 = -t296 / 0.2e1;
t421 = pkin(1) * t294;
t419 = pkin(1) - t246;
t185 = qJDD(5) * t223 + qJDD(1);
t152 = t324 * qJD(5);
t153 = t166 * qJD(5);
t102 = Icges(6,5) * t152 - Icges(6,6) * t153;
t103 = Icges(6,4) * t152 - Icges(6,2) * t153;
t104 = Icges(6,1) * t152 - Icges(6,4) * t153;
t19 = t102 * t223 + t103 * t324 + t104 * t166 + t152 * t98 - t153 * t97;
t33 = t166 * t98 + t223 * t96 + t324 * t97;
t416 = t33 * t185 + t19 * t192;
t415 = t176 * t53 + t429 * t57;
t414 = t176 * t55 + t429 * t59;
t413 = t176 * t96 + t429 * t97;
t25 = t166 * t61 + t223 * t53 + t324 * t57;
t407 = t25 * t106;
t26 = t166 * t63 + t223 * t55 + t324 * t59;
t406 = t26 * t107;
t355 = qJD(4) * t288;
t358 = qJD(3) * t289;
t239 = t292 * t355 + t358;
t395 = t239 * t294;
t394 = t288 * t289;
t378 = t137 - t182;
t377 = t359 * t419 - t222 - t263;
t375 = t189 - t246;
t374 = t190 - t245;
t205 = t228 * qJD(1);
t373 = -t206 * rSges(4,1) + t205 * rSges(4,2);
t202 = -t236 + t421;
t250 = -t283 + t421;
t372 = t202 - t250;
t269 = qJD(4) * t291 - qJD(2);
t371 = t239 * t296 - t269 * t294;
t349 = rSges(3,1) * t386;
t366 = rSges(3,2) * t390 + t296 * rSges(3,3);
t198 = t349 - t366;
t369 = -t250 - t198;
t368 = rSges(3,2) * t335 + rSges(3,3) * t359;
t264 = t296 * t358;
t279 = qJD(2) * t294;
t367 = t264 + t279;
t351 = qJD(1) * qJD(2);
t365 = qJDD(2) * t294 + t296 * t351;
t364 = qJ(2) * t359 + t279;
t363 = -qJD(1) * t250 + t279;
t353 = qJD(5) * t176;
t352 = -m(4) + t426;
t350 = qJDD(3) * t289;
t72 = qJD(1) * t435 + qJD(5) * t429;
t154 = t170 * qJD(5);
t204 = t227 * qJD(5);
t75 = -qJD(1) * t133 - t154 * t296 - t204 * t294;
t40 = t72 * rSges(6,1) + t75 * rSges(6,2) + t149 * rSges(6,3);
t66 = rSges(6,1) * t435 - t442;
t344 = -t149 * rSges(5,2) + (-t206 * t290 - t287 * t335) * rSges(5,1) - t205 * rSges(5,3);
t101 = -t176 * rSges(5,2) + (t231 * t290 + t287 * t389) * rSges(5,1) + t230 * rSges(5,3);
t341 = t372 + t450;
t214 = qJD(4) * t230;
t340 = t214 + t367;
t139 = t231 * rSges(4,1) - t230 * rSges(4,2) + rSges(4,3) * t389;
t233 = t245 * t359;
t339 = -t233 + t371;
t338 = -t138 + t372;
t337 = t296 * t350 + t365;
t336 = t264 + t364;
t333 = -rSges(3,1) * t292 - pkin(1);
t332 = -t354 / 0.2e1;
t331 = t354 / 0.2e1;
t330 = -t353 / 0.2e1;
t329 = t353 / 0.2e1;
t328 = -rSges(4,3) * t289 - t246;
t323 = -t296 * t374 + t181 + t341 - t402;
t238 = -qJDD(3) * t292 + qJDD(4) * t394;
t321 = qJD(1) * t202 + t264 + t363;
t320 = t343 + t341;
t34 = Icges(6,5) * t72 + Icges(6,6) * t75 + Icges(6,3) * t149;
t207 = qJD(1) * t346 - t291 * t360;
t73 = qJD(1) * t300 + qJD(5) * t128 + t207 * t293;
t74 = qJD(1) * t429 - t154 * t294 + t204 * t296;
t35 = Icges(6,5) * t73 + Icges(6,6) * t74 + Icges(6,3) * t150;
t36 = Icges(6,4) * t72 + Icges(6,2) * t75 + Icges(6,6) * t149;
t37 = Icges(6,4) * t73 + Icges(6,2) * t74 + Icges(6,6) * t150;
t38 = Icges(6,1) * t72 + Icges(6,4) * t75 + Icges(6,5) * t149;
t39 = Icges(6,1) * t73 + Icges(6,4) * t74 + Icges(6,5) * t150;
t318 = (t131 * t39 + t149 * t53 + t176 * t35 + t37 * t429 + t57 * t75 + t61 * t72) * t173 + t176 * (t131 * t38 + t149 * t55 + t176 * t34 + t36 * t429 + t59 * t75 + t63 * t72);
t317 = t173 * (t129 * t39 + t133 * t37 + t150 * t53 + t173 * t35 + t57 * t74 + t61 * t73) + t176 * (t129 * t38 + t133 * t36 + t150 * t55 + t173 * t34 + t59 * t74 + t63 * t73);
t10 = t152 * t63 - t153 * t59 + t166 * t38 + t223 * t34 + t324 * t36;
t9 = t152 * t61 - t153 * t57 + t166 * t39 + t223 * t35 + t324 * t37;
t316 = t10 * t176 + t173 * t9;
t253 = rSges(2,1) * t296 - rSges(2,2) * t294;
t251 = rSges(2,1) * t294 + rSges(2,2) * t296;
t315 = -t208 * rSges(4,1) + t207 * rSges(4,2);
t14 = t129 * t61 + t439;
t15 = t129 * t63 + t440;
t314 = t14 * t173 + t15 * t176;
t16 = t131 * t61 + t415;
t17 = t131 * t63 + t414;
t313 = t16 * t173 + t17 * t176;
t41 = rSges(6,1) * t73 + rSges(6,2) * t74 + rSges(6,3) * t150;
t312 = -t173 * t40 + t176 * t41;
t65 = rSges(6,1) * t129 + t442;
t123 = t131 * rSges(6,1);
t67 = t123 + t379;
t311 = -t173 * t67 + t176 * t65;
t310 = -qJD(4) * t205 + qJDD(4) * t230 + t337;
t307 = -t269 * t296 - t395;
t305 = -qJDD(2) * t296 + qJD(1) * (-pkin(1) * t360 + t364) + qJDD(1) * t252 + t294 * t351;
t304 = qJD(1) * t450 + t214 + t321;
t303 = qJD(1) * t430 - t309;
t232 = t245 * t360;
t302 = t232 - (t296 * t375 - t282) * qJD(1) + t307 - t362 - t263 + t377;
t301 = -rSges(5,1) * (t208 * t290 + t287 * t334) + rSges(5,2) * t150 - rSges(5,3) * t207;
t299 = qJDD(1) * t203 + t294 * t350 + t305 + (t360 * t419 + 0.2e1 * t264) * qJD(1);
t298 = qJDD(1) * t430 + qJD(4) * t207 + qJDD(4) * t228 + t299 + qJD(1) * (-t360 * t375 - t336 + t339);
t180 = t190 * t359;
t144 = qJD(1) * t369 + t279;
t124 = t132 * rSges(6,1);
t111 = rSges(6,1) * t324 - rSges(6,2) * t166;
t108 = Icges(6,5) * t324 - Icges(6,6) * t166;
t105 = rSges(6,1) * t152 - rSges(6,2) * t153;
t95 = qJDD(1) * t199 + qJD(1) * (-qJD(1) * t349 + t368) + t305;
t94 = qJD(1) * t452 + t369 * qJDD(1) + t365;
t93 = qJD(1) * t139 - t319;
t92 = qJD(1) * t338 + t367;
t86 = t234 - t444 + t378;
t84 = rSges(6,1) * t128 + rSges(6,2) * t435;
t83 = rSges(6,1) * t429 - rSges(6,2) * t131;
t68 = t124 + t379;
t64 = Icges(6,1) * t132 + t380;
t62 = Icges(6,1) * t435 - t438;
t60 = Icges(6,4) * t132 + t381;
t58 = -t437 + t454;
t48 = qJD(1) * t101 + t303;
t47 = qJD(1) * t320 + t340;
t46 = qJDD(1) * t139 + qJD(1) * (-rSges(4,3) * t335 + t373) + t299;
t45 = t338 * qJDD(1) + ((-rSges(4,3) * t359 - t357) * t289 + t315 + t377) * qJD(1) + t337;
t31 = -qJD(3) * t292 + qJD(5) * t311 + t289 * t355;
t30 = t131 * t98 + t413;
t29 = t129 * t98 + t441;
t28 = qJD(1) * t344 + qJDD(1) * t101 + t298;
t27 = t320 * qJDD(1) + (t301 + t302) * qJD(1) + t310;
t24 = qJD(1) * t86 + t192 * t67 - t353 * t99 + t303;
t23 = qJD(1) * t323 - t192 * t65 + t340 + t447;
t13 = t102 * t173 + t103 * t133 + t104 * t129 + t150 * t96 + t73 * t98 + t74 * t97;
t12 = t102 * t176 + t103 * t429 + t104 * t131 + t149 * t96 + t72 * t98 + t75 * t97;
t11 = qJD(5) * t312 - t106 * t67 + t107 * t65 + t238;
t8 = t298 + (-t180 + t233 + (-t141 + t189) * t360) * qJD(1) - t105 * t353 + qJDD(1) * t86 - t107 * t99 + t185 * t67 + t192 * t40;
t7 = t105 * t354 + t106 * t99 - t185 * t65 - t192 * t41 + t323 * qJDD(1) + (-t232 + (t308 + t182) * qJD(1) + t302) * qJD(1) + t310;
t6 = qJD(5) * t313 + t192 * t30;
t5 = qJD(5) * t314 + t192 * t29;
t1 = [t407 / 0.2e1 + t406 / 0.2e1 - m(2) * (-g(1) * t251 + g(2) * t253) + ((t132 * t98 + t413) * t192 + ((t129 * t62 + t132 * t63 + t133 * t58 + t14 + t414) * t176 + (t129 * t64 + t132 * t61 + t133 * t60 - t15 + t415 + t449) * t173) * qJD(5)) * t332 + t416 + t30 * t424 + t29 * t425 - (t192 * (t166 * t62 + t223 * t54 + t324 * t58 + t25) * t176 + (t192 * (t166 * t64 + t223 * t56 + t324 * t60 - t26) - t6) * t173) * qJD(5) / 0.2e1 + (t13 + t9) * t331 + (t5 + (t435 * t98 - t441) * t192 + ((t131 * t62 + t429 * t58 + t435 * t63 + t16 - t440 + t449) * t176 + (t131 * t64 + t429 * t60 + t435 * t61 - t17 - t439) * t173) * qJD(5)) * t330 + (t10 + t12) * t329 + (t7 * (t326 - t65) + t23 * (qJD(1) * t308 + t307 - t41) + t8 * (t123 + t451) - g(1) * (t326 + t66) - g(2) * (t124 + t451) - t23 * (-t192 * t68 + t309 + (t294 * t374 - t378 - t430) * qJD(1)) - ((t23 * t99 + t31 * (-t67 + t68)) * t176 + t31 * (-t65 - t66) * t173) * qJD(5) + (-t141 * t360 - t180 + t371 + t40 - qJD(1) * (t326 + t325) - t192 * t66 - t304 - t447) * t24) * m(6) + ((-qJD(1) * t343 - t189 * t360 - t304 + t339 + t344) * t48 + (-t395 + t232 + (-qJD(1) * t189 - t269) * t296 + t301 - (-t101 - t430) * qJD(1) - t309) * t47 + (-g(2) + t28) * (t101 + t443) + (-g(1) + t27) * (t343 - t325)) * m(5) + (-(-qJD(1) * t138 + t321 - t92) * t93 + t92 * (t315 + t362) + t93 * (t336 + t373) + (t92 * t328 * t296 + (-t92 * qJ(2) + t328 * t93) * t294) * qJD(1) + (t46 - g(2)) * (t139 + t370) + (t45 - g(1)) * (-t138 + t361)) * m(4) + (t144 * t280 - t452 * (t364 + t368) + (t144 * (rSges(3,2) * t289 + t333) * t296 + (t144 * (-rSges(3,3) - qJ(2)) - t452 * t333) * t294) * qJD(1) + (-qJD(1) * t198 - t144 + t363) * t452 + (t95 - g(2)) * (t199 + t252) + (t94 - g(1)) * (t333 * t294 + t283 + t366)) * m(3) + ((Icges(5,1) * t224 + 0.2e1 * Icges(5,5) * t394) * t224 + (-0.2e1 * Icges(5,4) * t224 + Icges(5,2) * t223 - 0.2e1 * Icges(5,6) * t394) * t223 + m(2) * (t251 ^ 2 + t253 ^ 2) + Icges(2,3) + (Icges(4,3) + Icges(3,2)) * t427 + ((Icges(4,1) * t291 ^ 2 + Icges(3,1) + (-0.2e1 * Icges(4,4) * t291 + (Icges(4,2) + Icges(5,3)) * t288) * t288) * t289 + 0.2e1 * (-Icges(4,5) * t291 + Icges(4,6) * t288 + Icges(3,4)) * t292) * t289) * qJDD(1); (-m(3) + t352) * (g(1) * t294 - g(2) * t296) + 0.2e1 * (t422 * t8 + t423 * t7) * m(6) + 0.2e1 * (t27 * t423 + t28 * t422) * m(5) + 0.2e1 * (t422 * t46 + t423 * t45) * m(4) + 0.2e1 * (t422 * t95 + t423 * t94) * m(3); t352 * (-g(3) * t292 + (g(1) * t296 + g(2) * t294) * t289) + m(4) * (qJDD(3) * t427 + t389 * t45 + t390 * t46) + m(5) * (-t238 * t292 + t27 * t389 + t28 * t390) + m(6) * (-t11 * t292 + t389 * t7 + t390 * t8); t426 * (g(1) * t230 + g(2) * t228 + g(3) * t394) + m(5) * (-t205 * t47 + t207 * t48 + t228 * t28 + t230 * t27 + t238 * t394) + m(6) * (t11 * t394 - t205 * t23 + t207 * t24 + t228 * t8 + t230 * t7) + 0.2e1 * (-m(5) * (-t228 * t47 + t230 * t48) / 0.2e1 - m(6) * (-t228 * t23 + t230 * t24) / 0.2e1) * qJD(1); t149 * t6 / 0.2e1 + t176 * (qJD(5) * t318 + t106 * t16 + t107 * t17 + t12 * t192 + t185 * t30) / 0.2e1 + (t223 * t30 + t313) * t424 + (t12 * t223 + t149 * t17 + t150 * t16 + t318) * t329 + t150 * t5 / 0.2e1 + t173 * (qJD(5) * t317 + t106 * t14 + t107 * t15 + t13 * t192 + t185 * t29) / 0.2e1 + (t223 * t29 + t314) * t425 + (t13 * t223 + t14 * t150 + t149 * t15 + t317) * t331 + t223 * (qJD(5) * t316 + t406 + t407 + t416) / 0.2e1 + t185 * (t173 * t25 + t176 * t26 + t223 * t33) / 0.2e1 + t192 * (t149 * t26 + t150 * t25 + t19 * t223 + t316) / 0.2e1 + ((t108 * t176 + t460 * t131) * t192 + (t459 * t131 + t461 * t176) * qJD(5) + t448 * t429) * t330 + ((t108 * t173 + t109 * t133 + t110 * t129 + t128 * t98 + t435 * t97) * t192 + ((t128 * t63 + t129 * t81 + t133 * t79 + t173 * t77 + t435 * t59) * t176 + (t128 * t61 + t129 * t82 + t133 * t80 + t435 * t57 + t456) * t173) * qJD(5)) * t332 - t192 * ((t108 * t223 + t460 * t166) * t192 + (t459 * t166 + t461 * t223) * qJD(5) + t448 * t324) / 0.2e1 + (t7 * (t173 * t99 - t223 * t65) + t23 * (t105 * t173 + t150 * t99 - t223 * t41) + t8 * (-t176 * t99 + t223 * t67) + t24 * (-t105 * t176 - t149 * t99 + t223 * t40) + t11 * t311 + t31 * (t149 * t65 - t150 * t67 + t312) - (-t23 * t84 + t24 * t83) * t192 - (t31 * (-t173 * t83 + t176 * t84) + (t173 * t23 - t176 * t24) * t111) * qJD(5) - g(1) * t83 - g(2) * t84 - g(3) * t111) * m(6);];
tau = t1;
