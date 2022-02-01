% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:04
% EndTime: 2022-01-23 08:59:18
% DurationCPUTime: 8.49s
% Computational Cost: add. (10032->526), mult. (24383->687), div. (0->0), fcn. (27758->10), ass. (0->275)
t248 = sin(pkin(8));
t251 = cos(pkin(8));
t254 = sin(qJ(1));
t252 = cos(pkin(7));
t256 = cos(qJ(1));
t343 = t252 * t256;
t200 = t248 * t254 + t251 * t343;
t247 = sin(pkin(9));
t250 = cos(pkin(9));
t249 = sin(pkin(7));
t347 = t249 * t256;
t157 = t200 * t247 - t250 * t347;
t346 = t250 * t252;
t194 = t247 * t249 + t251 * t346;
t341 = t256 * t248;
t150 = t194 * t254 - t250 * t341;
t344 = t252 * t254;
t197 = t248 * t344 + t251 * t256;
t253 = sin(qJ(5));
t255 = cos(qJ(5));
t113 = -t150 * t253 + t197 * t255;
t198 = t251 * t344 - t341;
t345 = t250 * t254;
t154 = t198 * t247 - t249 * t345;
t351 = t248 * t253;
t151 = t194 * t255 + t252 * t351;
t350 = t248 * t255;
t196 = t250 * t350 - t251 * t253;
t384 = -t151 * t254 + t256 * t196;
t402 = t154 * (Icges(6,5) * t113 + Icges(6,6) * t384);
t116 = t151 * t256 + t196 * t254;
t152 = -t194 * t253 + t252 * t350;
t195 = t250 * t351 + t251 * t255;
t382 = t152 * t256 - t254 * t195;
t75 = Icges(6,5) * t382 - Icges(6,6) * t116;
t406 = t157 * t75 + t402;
t349 = t249 * t251;
t193 = -t252 * t247 + t250 * t349;
t285 = -t193 * t253 + t249 * t350;
t147 = t193 * t255 + t249 * t351;
t394 = Icges(6,4) * t147;
t101 = Icges(6,1) * t285 - t394;
t192 = t247 * t349 + t346;
t91 = Icges(6,2) * t285 + Icges(6,6) * t192 + t394;
t405 = t101 - t91;
t114 = t150 * t255 + t197 * t253;
t118 = t152 * t254 + t195 * t256;
t386 = Icges(6,2) * t118 + Icges(6,6) * t154;
t53 = Icges(6,4) * t114 + t386;
t339 = Icges(6,2) * t382 + Icges(6,6) * t157;
t403 = Icges(6,4) * t116;
t55 = t339 + t403;
t79 = Icges(6,1) * t382 - t403;
t400 = Icges(6,4) * t384;
t80 = Icges(6,1) * t113 + t400;
t404 = (-t55 + t79) * t157 + (-t53 + t80) * t154;
t242 = t254 * qJ(2);
t246 = t256 * pkin(1);
t318 = t246 + t242;
t209 = qJD(1) * t318;
t240 = qJD(2) * t256;
t191 = -t240 + t209;
t271 = rSges(3,1) * t343 - rSges(3,2) * t347 + t254 * rSges(3,3);
t399 = -qJD(1) * t271 - t191;
t372 = t248 * qJ(4) + pkin(2);
t373 = t249 * qJ(3) + pkin(1);
t165 = (pkin(3) * t251 + t372) * t252 + t373;
t160 = t165 * t254;
t212 = pkin(4) * t250 + pkin(6) * t247 + pkin(3);
t288 = qJ(4) * t251 - qJ(2);
t166 = -t212 * t248 + t288;
t210 = -t248 * pkin(3) + t288;
t333 = -t166 + t210;
t125 = (t212 * t251 + t372) * t252 + pkin(1) + (pkin(4) * t247 - pkin(6) * t250 + qJ(3)) * t249;
t358 = t125 * t254;
t398 = t256 * t333 + t160 - t358;
t385 = Icges(6,6) * t118 + Icges(6,3) * t154;
t50 = Icges(6,5) * t384 - t385;
t306 = t252 * t341;
t199 = -t251 * t254 + t306;
t263 = (t194 * t256 + t248 * t345) * t255;
t117 = t199 * t253 + t263;
t340 = Icges(6,6) * t382 + Icges(6,3) * t157;
t52 = Icges(6,5) * t117 + t340;
t397 = t154 * t52 + t157 * t50;
t383 = Icges(6,4) * t285;
t100 = -Icges(6,2) * t147 + t383;
t168 = qJD(5) * t192 + qJD(1);
t387 = Icges(6,4) * t118 + Icges(6,5) * t154;
t57 = Icges(6,1) * t114 + t387;
t393 = Icges(6,4) * t382;
t338 = Icges(6,5) * t157 + t393;
t59 = Icges(6,1) * t116 + t338;
t77 = -Icges(6,2) * t116 + t393;
t78 = Icges(6,4) * t113 + Icges(6,2) * t384;
t92 = Icges(6,1) * t147 + Icges(6,5) * t192 + t383;
t396 = (t100 + t92) * t168 + ((t59 + t77) * t157 + (t57 + t78) * t154) * qJD(5);
t311 = qJD(5) * t154;
t93 = rSges(6,1) * t147 + rSges(6,2) * t285 + rSges(6,3) * t192;
t395 = t93 * t311;
t161 = t165 * t256;
t328 = t210 * t254 - t161;
t348 = t249 * t254;
t264 = rSges(5,2) * t154 - t197 * rSges(5,3) - rSges(5,1) * (t198 * t250 + t247 * t348);
t391 = rSges(6,2) * t118 + rSges(6,3) * t154;
t90 = Icges(6,5) * t147 + Icges(6,6) * t285 + Icges(6,3) * t192;
t390 = t118 * t91 + t154 * t90;
t51 = Icges(6,5) * t116 + t340;
t389 = t118 * t55 + t154 * t51;
t49 = Icges(6,5) * t114 + t385;
t388 = t118 * t53 + t49 * t154;
t124 = t125 * t256;
t273 = t166 * t254 - t124;
t211 = pkin(2) * t252 + t373;
t205 = t211 * t256;
t327 = t205 + t242;
t106 = -t327 - t328;
t314 = qJD(3) * t254;
t224 = t249 * t314;
t320 = t240 - t224;
t282 = -qJD(1) * (-t246 + t205) - t209 + t320;
t274 = -qJD(4) * t197 + t282;
t379 = t254 / 0.2e1;
t378 = -t256 / 0.2e1;
t377 = pkin(1) * t254;
t375 = qJD(5) / 0.2e1;
t374 = pkin(1) - t211;
t371 = t157 * t49 + t382 * t53;
t370 = t157 * t51 + t382 * t55;
t369 = t157 * t90 + t382 * t91;
t312 = qJD(4) * t248;
t315 = qJD(3) * t249;
t206 = t252 * t312 + t315;
t352 = t206 * t254;
t337 = rSges(6,2) * t382 + t157 * rSges(6,3);
t316 = qJD(1) * t256;
t335 = t316 * t374 - t191 - t224;
t334 = t165 - t211;
t179 = t197 * qJD(1);
t180 = t198 * qJD(1);
t332 = -t180 * rSges(4,1) + t179 * rSges(4,2);
t204 = t211 * t254;
t176 = -t204 + t377;
t243 = t256 * qJ(2);
t215 = -t243 + t377;
t331 = t176 - t215;
t309 = qJD(1) * qJD(2);
t317 = qJD(1) * t254;
t239 = qJD(2) * t254;
t322 = qJ(2) * t316 + t239;
t330 = qJD(1) * (-pkin(1) * t317 + t322) + t254 * t309;
t230 = qJD(4) * t251 - qJD(2);
t329 = t206 * t256 - t230 * t254;
t298 = t249 * t317;
t326 = rSges(3,2) * t298 + rSges(3,3) * t316;
t225 = t256 * t315;
t325 = t225 + t239;
t324 = rSges(3,2) * t348 + t256 * rSges(3,3);
t234 = t256 * t309;
t323 = -qJD(4) * t179 + t234;
t321 = -qJD(1) * t215 + t239;
t319 = t243 - t204;
t310 = qJD(5) * t157;
t133 = -t180 * t247 + t250 * t298;
t70 = qJD(1) * t384 + t382 * qJD(5);
t138 = t151 * qJD(5);
t178 = t196 * qJD(5);
t73 = -qJD(1) * t118 - t138 * t256 - t178 * t254;
t38 = t70 * rSges(6,1) + t73 * rSges(6,2) + t133 * rSges(6,3);
t308 = rSges(3,1) * t344;
t63 = t116 * rSges(6,1) + t337;
t305 = -t133 * rSges(5,2) + (-t180 * t250 - t247 * t298) * rSges(5,1) - t179 * rSges(5,3);
t94 = -t157 * rSges(5,2) + (t200 * t250 + t247 * t347) * rSges(5,1) + t199 * rSges(5,3);
t286 = t210 * t256 + t160;
t105 = t286 + t319;
t303 = -t105 + t331;
t186 = qJD(4) * t199;
t302 = t186 + t325;
t301 = t200 * rSges(4,1) - t199 * rSges(4,2) + rSges(4,3) * t347;
t202 = t210 * t316;
t300 = -t202 + t329;
t299 = t225 + t322;
t297 = t249 * t316;
t296 = -rSges(3,1) * t252 - pkin(1);
t295 = t133 * t375;
t182 = t200 * qJD(1);
t134 = t182 * t247 - t250 * t297;
t294 = t134 * t375;
t293 = -t311 / 0.2e1;
t292 = t311 / 0.2e1;
t291 = -t310 / 0.2e1;
t290 = t310 / 0.2e1;
t289 = -rSges(4,3) * t249 - t211;
t284 = t330 + (t317 * t374 + 0.2e1 * t225) * qJD(1);
t283 = qJD(1) * t176 + t225 + t321;
t181 = qJD(1) * t306 - t251 * t317;
t280 = -t182 * rSges(4,1) + t181 * rSges(4,2);
t136 = t285 * qJD(5);
t137 = t147 * qJD(5);
t98 = rSges(6,1) * t136 - rSges(6,2) * t137;
t279 = -t133 * t93 - t157 * t98;
t278 = t134 * t93 + t154 * t98;
t14 = t114 * t57 + t388;
t15 = t114 * t59 + t389;
t277 = t14 * t154 + t15 * t157;
t16 = t116 * t57 + t371;
t17 = t116 * t59 + t370;
t276 = t154 * t16 + t157 * t17;
t61 = rSges(6,1) * t114 + t391;
t275 = -t154 * t63 + t157 * t61;
t272 = -t230 * t256 - t352;
t270 = qJD(4) * t181 + qJD(1) * (-t317 * t334 - t299 + t300) + t284;
t269 = -qJD(1) * t105 + t186 + t283;
t268 = qJD(1) * t106 - t274;
t201 = t210 * t317;
t266 = t201 - (t256 * t334 - t242) * qJD(1) + t272 - t320 - t224 + t335;
t71 = qJD(1) * t263 + qJD(5) * t113 + t181 * t253;
t72 = qJD(1) * t382 - t138 * t254 + t178 * t256;
t39 = rSges(6,1) * t71 + rSges(6,2) * t72 + rSges(6,3) * t134;
t265 = -rSges(5,1) * (t182 * t250 + t247 * t297) + rSges(5,2) * t134 - rSges(5,3) * t181;
t262 = -rSges(4,1) * t198 + rSges(4,2) * t197 - rSges(4,3) * t348;
t32 = Icges(6,5) * t70 + Icges(6,6) * t73 + Icges(6,3) * t133;
t33 = Icges(6,5) * t71 + Icges(6,6) * t72 + Icges(6,3) * t134;
t34 = Icges(6,4) * t70 + Icges(6,2) * t73 + Icges(6,6) * t133;
t35 = Icges(6,4) * t71 + Icges(6,2) * t72 + Icges(6,6) * t134;
t36 = Icges(6,1) * t70 + Icges(6,4) * t73 + Icges(6,5) * t133;
t37 = Icges(6,1) * t71 + Icges(6,4) * t72 + Icges(6,5) * t134;
t261 = (t116 * t37 + t133 * t49 + t157 * t33 + t35 * t382 + t53 * t73 + t57 * t70) * t154 + t133 * t17 + t134 * t16 + t157 * (t116 * t36 + t133 * t51 + t157 * t32 + t34 * t382 + t55 * t73 + t59 * t70);
t260 = t133 * t15 + t134 * t14 + t154 * (t114 * t37 + t118 * t35 + t134 * t49 + t154 * t33 + t53 * t72 + t57 * t71) + t157 * (t114 * t36 + t118 * t34 + t134 * t51 + t154 * t32 + t55 * t72 + t59 * t71);
t25 = t147 * t57 + t192 * t49 + t285 * t53;
t26 = t147 * t59 + t192 * t51 + t285 * t55;
t7 = t136 * t57 - t137 * t53 + t147 * t37 + t192 * t33 + t285 * t35;
t8 = t136 * t59 - t137 * t55 + t147 * t36 + t192 * t32 + t285 * t34;
t259 = t26 * t133 + t25 * t134 + t7 * t154 + t8 * t157;
t258 = t133 * t61 - t134 * t63 - t154 * t38 + t157 * t39;
t174 = t308 - t324;
t159 = t166 * t316;
t128 = t239 + (-t174 - t215) * qJD(1);
t123 = t399 * qJD(1) + t234;
t122 = qJD(1) * (-qJD(1) * t308 + t326) + t330;
t99 = Icges(6,5) * t285 - Icges(6,6) * t147;
t97 = Icges(6,1) * t136 - Icges(6,4) * t137;
t96 = Icges(6,4) * t136 - Icges(6,2) * t137;
t95 = Icges(6,5) * t136 - Icges(6,6) * t137;
t89 = qJD(1) * t301 - t282;
t88 = (t262 + t331) * qJD(1) + t325;
t82 = rSges(6,1) * t113 + rSges(6,2) * t384;
t81 = rSges(6,1) * t382 - rSges(6,2) * t116;
t67 = t234 + ((-rSges(4,3) * t316 - t314) * t249 + t280 + t335) * qJD(1);
t66 = qJD(1) * (-rSges(4,3) * t298 + t332) + t284;
t64 = rSges(6,1) * t117 + t337;
t62 = rSges(6,1) * t384 - t391;
t60 = Icges(6,1) * t117 + t338;
t58 = Icges(6,1) * t384 - t387;
t56 = Icges(6,4) * t117 + t339;
t54 = -t386 + t400;
t44 = qJD(1) * t94 + t268;
t43 = (t264 + t303) * qJD(1) + t302;
t31 = (t265 + t266) * qJD(1) + t323;
t30 = qJD(1) * t305 + t270;
t29 = -qJD(3) * t252 + qJD(5) * t275 + t249 * t312;
t28 = t116 * t92 + t369;
t27 = t114 * t92 + t390;
t24 = qJD(1) * (t328 - t273) - t93 * t310 + t168 * t63 + t268;
t23 = t395 - t168 * t61 + (t303 + t398) * qJD(1) + t302;
t19 = t136 * t92 - t137 * t91 + t147 * t97 + t192 * t95 + t285 * t96;
t18 = t19 * t168;
t13 = t168 * t38 + t279 * qJD(5) + (-t159 + t202 + (-t125 + t165) * t317) * qJD(1) + t270;
t12 = -t168 * t39 + t278 * qJD(5) + (-t201 + (t273 + t161) * qJD(1) + t266) * qJD(1) + t323;
t11 = t114 * t97 + t118 * t96 + t134 * t90 + t154 * t95 + t71 * t92 + t72 * t91;
t10 = t116 * t97 + t133 * t90 + t157 * t95 + t382 * t96 + t70 * t92 + t73 * t91;
t9 = t258 * qJD(5);
t6 = qJD(5) * t276 + t168 * t28;
t5 = qJD(5) * t277 + t168 * t27;
t1 = [t18 + ((t117 * t92 + t369) * t168 + ((t114 * t58 + t117 * t59 + t118 * t54 + t14 + t370) * t157 + (t114 * t60 + t117 * t57 + t118 * t56 - t15 + t371 + t397) * t154) * qJD(5)) * t293 - (t168 * (t147 * t58 + t192 * t50 + t285 * t54 + t25) * t157 + (t168 * (t147 * t60 + t192 * t52 + t285 * t56 - t26) - t6) * t154) * qJD(5) / 0.2e1 + (t26 + t28) * t295 + (t25 + t27) * t294 + (t7 + t11) * t292 + (t5 + (t384 * t92 - t390) * t168 + ((t116 * t58 + t382 * t54 + t384 * t59 + t16 - t389 + t397) * t157 + (t116 * t60 + t382 * t56 + t384 * t57 - t17 - t388) * t154) * qJD(5)) * t291 + (t8 + t10) * t290 + (-t23 * (-t168 * t64 + t274 + (-t254 * t333 - t106 - t124 + t161) * qJD(1)) - ((t23 * t93 + t29 * (-t63 + t64)) * t157 + t29 * (-t61 - t62) * t154) * qJD(5) + t12 * (-t166 * t256 - t358 - t61) + t23 * (qJD(1) * t273 + t272 - t39) + t13 * (t63 - t273) + (-qJD(1) * t398 - t125 * t317 - t168 * t62 - t159 - t269 + t329 + t38 - t395) * t24) * m(6) + (t31 * (t264 - t286) + t30 * (t94 - t328) + (-qJD(1) * t264 - t165 * t317 - t269 + t300 + t305) * t44 + (-(-t106 - t94) * qJD(1) - t274 - t352 + t201 + (-qJD(1) * t165 - t230) * t256 + t265) * t43) * m(5) + (-(qJD(1) * t262 + t283 - t88) * t89 + t67 * (t262 + t319) + t88 * (t280 + t320) + t66 * (t301 + t327) + t89 * (t299 + t332) + (t88 * t289 * t256 + (-t88 * qJ(2) + t289 * t89) * t254) * qJD(1)) * m(4) + ((-qJD(1) * t174 - t128 + t321) * t399 + t123 * (t254 * t296 + t243 + t324) + t128 * t240 + t122 * (t271 + t318) - t399 * (t322 + t326) + (t128 * (rSges(3,2) * t249 + t296) * t256 + (t128 * (-rSges(3,3) - qJ(2)) - t399 * t296) * t254) * qJD(1)) * m(3); 0.2e1 * (t12 * t379 + t13 * t378) * m(6) + 0.2e1 * (t30 * t378 + t31 * t379) * m(5) + 0.2e1 * (t66 * t378 + t379 * t67) * m(4) + 0.2e1 * (t122 * t378 + t123 * t379) * m(3); -m(6) * t9 * t252 + 0.2e1 * (m(4) * (t254 * t66 + t256 * t67) / 0.2e1 + m(5) * (t254 * t30 + t256 * t31) / 0.2e1 + m(6) * (t12 * t256 + t13 * t254) / 0.2e1) * t249; m(5) * (-t179 * t43 + t181 * t44 + t197 * t30 + t199 * t31) + m(6) * (t248 * t249 * t9 + t12 * t199 + t13 * t197 - t179 * t23 + t181 * t24) + 0.2e1 * (-m(5) * (-t197 * t43 + t199 * t44) / 0.2e1 - m(6) * (-t197 * t23 + t199 * t24) / 0.2e1) * qJD(1); t133 * t6 / 0.2e1 + t157 * (t261 * qJD(5) + t10 * t168) / 0.2e1 + (t192 * t28 + t276) * t295 + (t10 * t192 + t261) * t290 + t134 * t5 / 0.2e1 + t154 * (qJD(5) * t260 + t11 * t168) / 0.2e1 + (t192 * t27 + t277) * t294 + (t11 * t192 + t260) * t292 + t192 * (t259 * qJD(5) + t18) / 0.2e1 + t168 * (t19 * t192 + t259) / 0.2e1 + ((t405 * t116 + t157 * t99) * t168 + (t404 * t116 + t406 * t157) * qJD(5) + t396 * t382) * t291 + ((t100 * t118 + t101 * t114 + t113 * t92 + t154 * t99 + t384 * t91) * t168 + ((t113 * t59 + t114 * t79 + t118 * t77 + t154 * t75 + t384 * t55) * t157 + (t113 * t57 + t114 * t80 + t118 * t78 + t384 * t53 + t402) * t154) * qJD(5)) * t293 - t168 * ((t405 * t147 + t192 * t99) * t168 + (t404 * t147 + t406 * t192) * qJD(5) + t396 * t285) / 0.2e1 + (t12 * (t154 * t93 - t192 * t61) + t23 * (-t192 * t39 + t278) + t13 * (-t157 * t93 + t192 * t63) + t24 * (t192 * t38 + t279) + t9 * t275 + t29 * t258 - (-t23 * t82 + t24 * t81) * t168 - (t29 * (-t154 * t81 + t157 * t82) + (t154 * t23 - t157 * t24) * (rSges(6,1) * t285 - rSges(6,2) * t147)) * qJD(5)) * m(6);];
tauc = t1(:);
