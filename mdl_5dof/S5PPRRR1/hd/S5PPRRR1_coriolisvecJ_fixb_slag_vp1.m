% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:37
% EndTime: 2019-12-05 15:12:49
% DurationCPUTime: 9.01s
% Computational Cost: add. (16029->466), mult. (15481->739), div. (0->0), fcn. (14767->8), ass. (0->254)
t412 = -Icges(5,1) + Icges(5,2);
t227 = pkin(9) + qJ(3);
t222 = qJ(4) + t227;
t215 = cos(t222);
t230 = cos(pkin(8));
t233 = cos(qJ(5));
t352 = t230 * t233;
t229 = sin(pkin(8));
t232 = sin(qJ(5));
t355 = t229 * t232;
t202 = t215 * t352 + t355;
t214 = sin(t222);
t228 = qJD(3) + qJD(4);
t361 = t214 * t228;
t331 = t232 * t361;
t120 = -qJD(5) * t202 + t230 * t331;
t353 = t230 * t232;
t354 = t229 * t233;
t201 = -t215 * t353 + t354;
t330 = t233 * t361;
t121 = qJD(5) * t201 - t230 * t330;
t358 = t215 * t228;
t356 = t228 * t230;
t328 = t215 * t356;
t63 = Icges(6,5) * t121 + Icges(6,6) * t120 + Icges(6,3) * t328;
t359 = t214 * t230;
t82 = Icges(6,5) * t202 + Icges(6,6) * t201 + Icges(6,3) * t359;
t272 = t214 * t63 + t358 * t82;
t65 = Icges(6,4) * t121 + Icges(6,2) * t120 + Icges(6,6) * t328;
t67 = Icges(6,1) * t121 + Icges(6,4) * t120 + Icges(6,5) * t328;
t368 = Icges(6,4) * t202;
t84 = Icges(6,2) * t201 + Icges(6,6) * t359 + t368;
t192 = Icges(6,4) * t201;
t86 = Icges(6,1) * t202 + Icges(6,5) * t359 + t192;
t20 = t120 * t84 + t121 * t86 + t201 * t65 + t202 * t67 + t230 * t272;
t218 = qJD(3) * t229;
t211 = qJD(4) * t229 + t218;
t341 = qJD(5) * t214;
t182 = t230 * t341 + t211;
t387 = t182 / 0.2e1;
t411 = t20 * t387;
t309 = rSges(6,1) * t233 - rSges(6,2) * t232;
t142 = -rSges(6,3) * t215 + t214 * t309;
t338 = qJD(5) * t229;
t183 = t214 * t338 - t356;
t340 = qJD(5) * t215;
t199 = -t215 * t355 - t352;
t200 = t215 * t354 - t353;
t360 = t214 * t229;
t87 = rSges(6,1) * t200 + rSges(6,2) * t199 + rSges(6,3) * t360;
t410 = -t142 * t183 - t87 * t340;
t385 = t183 / 0.2e1;
t407 = 0.2e1 * Icges(5,4) * t214 + t412 * t215;
t406 = -0.2e1 * Icges(5,4) * t215 + t412 * t214;
t118 = -qJD(5) * t200 + t229 * t331;
t119 = qJD(5) * t199 - t229 * t330;
t357 = t228 * t229;
t329 = t215 * t357;
t62 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t329;
t81 = Icges(6,5) * t200 + Icges(6,6) * t199 + Icges(6,3) * t360;
t273 = t214 * t62 + t358 * t81;
t64 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t329;
t66 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t329;
t369 = Icges(6,4) * t200;
t83 = Icges(6,2) * t199 + Icges(6,6) * t360 + t369;
t191 = Icges(6,4) * t199;
t85 = Icges(6,1) * t200 + Icges(6,5) * t360 + t191;
t19 = t120 * t83 + t121 * t85 + t201 * t64 + t202 * t66 + t230 * t273;
t366 = Icges(6,4) * t233;
t289 = -Icges(6,2) * t232 + t366;
t132 = -Icges(6,6) * t215 + t214 * t289;
t367 = Icges(6,4) * t232;
t295 = Icges(6,1) * t233 - t367;
t134 = -Icges(6,5) * t215 + t214 * t295;
t285 = Icges(6,5) * t233 - Icges(6,6) * t232;
t130 = -Icges(6,3) * t215 + t214 * t285;
t261 = t285 * t215;
t284 = -Icges(6,5) * t232 - Icges(6,6) * t233;
t71 = t228 * t261 + (Icges(6,3) * t228 + qJD(5) * t284) * t214;
t269 = t130 * t358 + t214 * t71;
t41 = t201 * t83 + t202 * t85 + t359 * t81;
t378 = t229 * t41;
t42 = t201 * t84 + t202 * t86 + t359 * t82;
t304 = t230 * t42 + t378;
t262 = t289 * t215;
t288 = -Icges(6,2) * t233 - t367;
t72 = t228 * t262 + (Icges(6,6) * t228 + qJD(5) * t288) * t214;
t263 = t295 * t215;
t294 = -Icges(6,1) * t232 - t366;
t73 = t228 * t263 + (Icges(6,5) * t228 + qJD(5) * t294) * t214;
t252 = (-t120 * t132 - t121 * t134 - t201 * t72 - t202 * t73 + t228 * t304 - t230 * t269) * t215;
t57 = t130 * t359 + t132 * t201 + t134 * t202;
t405 = t411 + t19 * t385 + (t361 * t57 + t252) * qJD(5) / 0.2e1;
t286 = -Icges(5,5) * t214 - Icges(5,6) * t215;
t175 = t286 * t229;
t176 = t286 * t230;
t404 = t175 * t356 / 0.2e1 - t176 * t211 / 0.2e1;
t225 = t229 ^ 2;
t226 = t230 ^ 2;
t345 = t225 + t226;
t136 = t228 * t175;
t137 = t228 * t176;
t399 = -Icges(5,5) * t229 + t407 * t230;
t401 = -Icges(5,6) * t229 + t406 * t230;
t245 = (t399 * t214 + t401 * t215) * t228;
t400 = Icges(5,5) * t230 + t407 * t229;
t402 = Icges(5,6) * t230 + t406 * t229;
t246 = (t400 * t214 + t402 * t215) * t228;
t403 = -(t245 + t136) * t229 / 0.2e1 + (t137 / 0.2e1 - t246 / 0.2e1) * t230;
t220 = sin(t227);
t221 = cos(t227);
t208 = rSges(4,1) * t220 + rSges(4,2) * t221;
t390 = qJD(3) * t208 * t345;
t372 = Icges(4,4) * t221;
t373 = Icges(4,4) * t220;
t398 = (-t220 * (-Icges(4,2) * t221 - t373) + t221 * (-Icges(4,1) * t220 - t372)) * qJD(3);
t39 = t199 * t83 + t200 * t85 + t360 * t81;
t40 = t199 * t84 + t200 * t86 + t360 * t82;
t56 = t130 * t360 + t132 * t199 + t134 * t200;
t397 = t230 * (t182 * t42 + t183 * t41 - t340 * t57) + t229 * (t182 * t40 + t183 * t39 - t340 * t56);
t116 = t229 * t142;
t206 = pkin(4) * t214 - pkin(7) * t215;
t268 = t206 * t229;
t146 = t228 * t268;
t185 = t206 * t230;
t147 = t228 * t185;
t327 = t215 * t338;
t207 = pkin(4) * t215 + pkin(7) * t214;
t184 = t207 * t229;
t186 = t207 * t230;
t381 = pkin(3) * t221;
t108 = -pkin(6) * t230 + t229 * t381;
t109 = pkin(6) * t229 + t230 * t381;
t342 = qJD(3) * t230;
t333 = t108 * t218 + t109 * t342 + qJD(1);
t88 = rSges(6,1) * t202 + rSges(6,2) * t201 + rSges(6,3) * t359;
t38 = t182 * t87 - t183 * t88 + t184 * t211 + t186 * t356 + t333;
t68 = rSges(6,1) * t119 + rSges(6,2) * t118 + rSges(6,3) * t329;
t69 = rSges(6,1) * t121 + rSges(6,2) * t120 + rSges(6,3) * t328;
t396 = (t182 * t116 + t185 * t356 + t211 * t268 + t327 * t88 + (-t146 + t68) * t229 + (-t147 + t69 + t410) * t230) * t38;
t204 = rSges(5,1) * t214 + rSges(5,2) * t215;
t266 = t229 * t204;
t144 = t228 * t266;
t181 = t204 * t230;
t145 = t228 * t181;
t205 = rSges(5,1) * t215 - rSges(5,2) * t214;
t153 = -rSges(5,3) * t230 + t205 * t229;
t154 = rSges(5,3) * t229 + t205 * t230;
t395 = (-t229 * t144 - t230 * t145 + t181 * t356 + t211 * t266) * (t153 * t211 + t154 * t356 + t333);
t300 = -t232 * t84 + t233 * t86;
t301 = -t232 * t83 + t233 * t85;
t389 = -(-t130 * t230 - t300) * t182 - (-t130 * t229 - t301) * t183;
t240 = t182 * (-Icges(6,2) * t202 + t192 + t86) + t183 * (-Icges(6,2) * t200 + t191 + t85) - t340 * (t288 * t214 + t134);
t388 = -t182 / 0.2e1;
t386 = -t183 / 0.2e1;
t382 = pkin(3) * t220;
t380 = t229 * t108 + t230 * t109;
t377 = t230 * t40;
t172 = t207 * t228;
t265 = t309 * t215;
t308 = -rSges(6,1) * t232 - rSges(6,2) * t233;
t74 = t228 * t265 + (rSges(6,3) * t228 + qJD(5) * t308) * t214;
t374 = -t172 - t74;
t362 = t130 * t215;
t350 = t229 * t153 + t230 * t154;
t348 = -t142 - t206;
t344 = qJD(2) * t230;
t339 = qJD(5) * t228;
t234 = qJD(3) ^ 2;
t337 = t234 * t381;
t336 = qJD(3) * t381;
t335 = pkin(3) * t218;
t334 = pkin(3) * t342;
t325 = -t340 / 0.2e1;
t324 = t340 / 0.2e1;
t323 = t339 / 0.2e1;
t322 = -t204 - t382;
t320 = t229 * t337;
t319 = t230 * t337;
t317 = (t186 + t88) * t230 + (t184 + t87) * t229;
t316 = t221 * t335;
t315 = t221 * t334;
t313 = t214 * t323;
t312 = t215 * t323;
t311 = t348 - t382;
t165 = t205 * t228;
t310 = -t165 - t336;
t307 = t82 * t182 + t81 * t183;
t36 = -t320 - t172 * t211 - t182 * t74 + (t88 * t361 + (-t142 * t356 - t69) * t215) * qJD(5);
t37 = -t319 - t172 * t356 + t183 * t74 + (-t87 * t361 + (t142 * t357 + t68) * t215) * qJD(5);
t306 = t229 * t37 - t230 * t36;
t305 = t229 * t39 + t377;
t47 = t214 * t301 - t215 * t81;
t48 = t214 * t300 - t215 * t82;
t303 = t47 * t229 + t48 * t230;
t302 = -t229 * t88 + t230 * t87;
t299 = Icges(4,1) * t221 - t373;
t293 = -Icges(4,2) * t220 + t372;
t287 = -Icges(4,5) * t220 - Icges(4,6) * t221;
t283 = -t132 * t232 + t134 * t233;
t279 = t234 * t345 * t382;
t278 = t229 * t312;
t277 = t230 * t312;
t276 = -t336 + t374;
t219 = qJD(2) * t229;
t274 = -t220 * t334 + t219;
t267 = Icges(6,3) * t214 + t261 - t283;
t258 = qJD(3) * t287;
t257 = -t220 * t335 - t344;
t143 = rSges(6,3) * t214 + t265;
t256 = -t116 * t340 + t142 * t327 + t183 * t143 - t207 * t356 - t341 * t87;
t255 = -(Icges(6,5) * t201 - Icges(6,6) * t202) * t182 - t183 * (Icges(6,5) * t199 - Icges(6,6) * t200) + t284 * t214 * t340;
t254 = (t228 * t303 - (t228 * t283 - t71) * t215 - (t130 * t228 - t232 * t72 + t233 * t73 + (-t132 * t233 - t134 * t232) * qJD(5)) * t214) * t215;
t253 = (-t118 * t132 - t119 * t134 - t199 * t72 - t200 * t73 + t228 * t305 - t229 * t269) * t215;
t247 = t214 * t255;
t243 = (-(-Icges(4,6) * t230 + t229 * t293) * t221 - (-Icges(4,5) * t230 + t229 * t299) * t220) * qJD(3) + t398 * t229;
t242 = (-(Icges(4,6) * t229 + t230 * t293) * t221 - (Icges(4,5) * t229 + t230 * t299) * t220) * qJD(3) + t398 * t230;
t241 = -t143 * t182 - t207 * t211 + t88 * t341;
t239 = (Icges(6,1) * t201 - t368 - t84) * t182 + (Icges(6,1) * t199 - t369 - t83) * t183 - (t294 * t214 - t132) * t340;
t11 = (t228 * t301 - t62) * t215 + (t228 * t81 - t232 * t64 + t233 * t66 + (-t232 * t85 - t233 * t83) * qJD(5)) * t214;
t112 = t132 * t229;
t113 = t132 * t230;
t114 = t134 * t229;
t115 = t134 * t230;
t12 = (t228 * t300 - t63) * t215 + (t228 * t82 - t232 * t65 + t233 * t67 + (-t232 * t86 - t233 * t84) * qJD(5)) * t214;
t133 = Icges(6,6) * t214 + t262;
t135 = Icges(6,5) * t214 + t263;
t17 = t118 * t83 + t119 * t85 + t199 * t64 + t200 * t66 + t229 * t273;
t18 = t118 * t84 + t119 * t86 + t199 * t65 + t200 * t67 + t229 * t272;
t58 = t214 * t283 - t362;
t22 = t182 * t48 + t183 * t47 - t340 * t58;
t235 = (-t267 * t340 - t389) * t214;
t237 = (t401 * t211 - t402 * t356) * t215 + (t399 * t211 - t400 * t356) * t214;
t3 = t183 * t17 + t18 * t182 + (t361 * t56 + t253) * qJD(5);
t238 = -t22 * t341 / 0.2e1 + ((-t113 * t201 - t115 * t202) * t182 + (-t112 * t201 - t114 * t202) * t183 + (t57 * t214 + (-t133 * t201 - t135 * t202 + t378) * t215) * qJD(5)) * t388 + ((-t113 * t199 - t115 * t200) * t182 + (-t112 * t199 - t114 * t200) * t183 + (t56 * t214 + (-t133 * t199 - t135 * t200 + t377) * t215) * qJD(5)) * t386 + (((t113 * t232 - t115 * t233 + t82) * t182 + (t112 * t232 - t114 * t233 + t81) * t183 + t58 * qJD(5)) * t214 + ((t267 * t215 + (t133 * t232 - t135 * t233 - t130) * t214 + t303) * qJD(5) + t389) * t215) * t324 + t397 * t325 + (t18 * t385 + t40 * t278 + t42 * t277 + t48 * t313 + (((t39 - t362) * qJD(5) + t307) * t215 + t235) * t386 + t411 + t405 + t12 * t325 + (t237 / 0.2e1 + t403) * t356) * t229 + (-t17 * t385 - t39 * t278 - t41 * t277 - t47 * t313 + (((t42 - t362) * qJD(5) + t307) * t215 + t235) * t388 - t19 * t387 - t3 / 0.2e1 - t11 * t325 + (-t136 * t230 + t229 * t246 + t404) * t356) * t230 + ((t137 * t229 + t230 * t245 + t404) * t229 + (-t237 / 0.2e1 + t403) * t230) * t211;
t194 = t287 * t230;
t193 = t287 * t229;
t190 = t308 * t214;
t167 = t230 * t258;
t166 = t229 * t258;
t162 = -t208 * t218 - t344;
t161 = -t208 * t342 + t219;
t106 = rSges(6,1) * t201 - rSges(6,2) * t202;
t105 = rSges(6,1) * t199 - rSges(6,2) * t200;
t96 = -t204 * t211 + t257;
t95 = -t204 * t356 + t274;
t92 = -t165 * t356 - t319;
t91 = -t165 * t211 - t320;
t80 = t390 * qJD(3);
t59 = -t144 * t211 - t145 * t356 - t279;
t55 = -t142 * t182 - t206 * t211 - t340 * t88 + t257;
t54 = -t206 * t356 + t274 - t410;
t23 = t215 * t302 * t339 - t146 * t211 - t147 * t356 + t182 * t68 - t183 * t69 - t279;
t1 = [-m(4) * t80 + m(5) * t59 + m(6) * t23; m(5) * (t229 * t92 - t230 * t91) + m(6) * t306; -(t229 * (-t167 * t230 + t229 * t242) - t230 * (-t166 * t230 + t229 * t243)) * t342 + (t229 * (t167 * t229 + t230 * t242) - t230 * (t166 * t229 + t230 * t243)) * t218 + t238 + (t193 * qJD(3) * t226 - t194 * t218 * t230) * t342 / 0.2e1 - (t194 * qJD(3) * t225 - t229 * t193 * t342) * t218 / 0.2e1 + (t23 * (t317 + t380) + (t276 * t54 + t311 * t37) * t230 + (t276 * t55 + t311 * t36) * t229 - t54 * (t256 - t315) - t55 * (t241 - t316) + t396) * m(6) + (-t95 * (-t205 * t356 - t315) - t96 * (-t205 * t211 - t316) + t59 * (t350 + t380) + (t310 * t95 + t322 * t92) * t230 + (t310 * t96 + t322 * t91) * t229 + t395) * m(5) + (-t80 * t345 + (-t161 * t230 - t162 * t229 + t390) * qJD(3) + t161 * t342 + t162 * t218) * (rSges(4,1) * t221 - rSges(4,2) * t220) * m(4); t238 + (-t241 * t55 - t256 * t54 + t23 * t317 + (t348 * t37 + t374 * t54) * t230 + (t348 * t36 + t374 * t55) * t229 + t396) * m(6) + (-(-t211 * t96 - t356 * t95) * t205 + t59 * t350 + (-t229 * t91 - t230 * t92) * t204 + (-t229 * t96 - t230 * t95) * t165 + t395) * m(5); t359 * t405 + (t214 * t304 - t215 * t57) * t277 + (t252 + (t19 * t229 + t20 * t230 + t228 * t57) * t214) * t387 + t3 * t360 / 0.2e1 + (t214 * t305 - t215 * t56) * t278 + (t253 + (t17 * t229 + t18 * t230 + t228 * t56) * t214) * t385 + t22 * t361 / 0.2e1 - t215 * (t11 * t183 + t12 * t182 + (t361 * t58 + t254) * qJD(5)) / 0.2e1 + (t214 * t303 - t215 * t58) * t313 + (t254 + (t11 * t229 + t12 * t230 + t228 * t58) * t214) * t325 + (t201 * t240 + t202 * t239 - t230 * t247) * t388 + (t199 * t240 + t200 * t239 - t229 * t247) * t386 + (t255 * t215 + (-t232 * t240 + t239 * t233) * t214) * t324 + t397 * t358 / 0.2e1 + ((-t36 * t88 + t37 * t87 + t54 * t68 - t55 * t69 + (t38 * t302 + (t229 * t54 - t230 * t55) * t142) * t228) * t215 + (t54 * (-t228 * t87 + t229 * t74) + t55 * (t228 * t88 - t230 * t74) + t23 * t302 + t38 * (-t229 * t69 + t230 * t68) + t306 * t142) * t214 - t54 * (t105 * t340 + t183 * t190) - t55 * (-t106 * t340 - t182 * t190) - t38 * (t105 * t182 - t106 * t183)) * m(6);];
tauc = t1(:);
