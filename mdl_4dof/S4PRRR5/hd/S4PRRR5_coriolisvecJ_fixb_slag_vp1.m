% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:35
% EndTime: 2019-12-31 16:33:45
% DurationCPUTime: 8.14s
% Computational Cost: add. (10925->458), mult. (15203->733), div. (0->0), fcn. (14549->8), ass. (0->252)
t406 = -Icges(4,1) + Icges(4,2);
t219 = qJ(2) + qJ(3);
t215 = cos(t219);
t221 = cos(pkin(7));
t224 = cos(qJ(4));
t343 = t221 * t224;
t220 = sin(pkin(7));
t222 = sin(qJ(4));
t346 = t220 * t222;
t194 = t215 * t343 + t346;
t214 = sin(t219);
t218 = qJD(2) + qJD(3);
t352 = t214 * t218;
t324 = t222 * t352;
t103 = -qJD(4) * t194 + t221 * t324;
t344 = t221 * t222;
t345 = t220 * t224;
t193 = -t215 * t344 + t345;
t323 = t224 * t352;
t104 = qJD(4) * t193 - t221 * t323;
t349 = t215 * t218;
t347 = t218 * t221;
t321 = t215 * t347;
t63 = Icges(5,5) * t104 + Icges(5,6) * t103 + Icges(5,3) * t321;
t350 = t214 * t221;
t81 = Icges(5,5) * t194 + Icges(5,6) * t193 + Icges(5,3) * t350;
t264 = t214 * t63 + t349 * t81;
t65 = Icges(5,4) * t104 + Icges(5,2) * t103 + Icges(5,6) * t321;
t67 = Icges(5,1) * t104 + Icges(5,4) * t103 + Icges(5,5) * t321;
t359 = Icges(5,4) * t194;
t83 = Icges(5,2) * t193 + Icges(5,6) * t350 + t359;
t168 = Icges(5,4) * t193;
t85 = Icges(5,1) * t194 + Icges(5,5) * t350 + t168;
t20 = t103 * t83 + t104 * t85 + t193 * t65 + t194 * t67 + t221 * t264;
t213 = qJD(2) * t220;
t206 = qJD(3) * t220 + t213;
t333 = qJD(4) * t214;
t161 = t221 * t333 + t206;
t380 = t161 / 0.2e1;
t405 = t20 * t380;
t299 = rSges(5,1) * t224 - rSges(5,2) * t222;
t146 = -rSges(5,3) * t215 + t214 * t299;
t330 = qJD(4) * t220;
t162 = t214 * t330 - t347;
t332 = qJD(4) * t215;
t191 = -t215 * t346 - t343;
t192 = t215 * t345 - t344;
t351 = t214 * t220;
t86 = rSges(5,1) * t192 + rSges(5,2) * t191 + rSges(5,3) * t351;
t404 = -t146 * t162 - t86 * t332;
t378 = t162 / 0.2e1;
t401 = 0.2e1 * Icges(4,4) * t214 + t406 * t215;
t400 = -0.2e1 * Icges(4,4) * t215 + t406 * t214;
t101 = -qJD(4) * t192 + t220 * t324;
t102 = qJD(4) * t191 - t220 * t323;
t348 = t218 * t220;
t322 = t215 * t348;
t62 = Icges(5,5) * t102 + Icges(5,6) * t101 + Icges(5,3) * t322;
t80 = Icges(5,5) * t192 + Icges(5,6) * t191 + Icges(5,3) * t351;
t265 = t214 * t62 + t349 * t80;
t64 = Icges(5,4) * t102 + Icges(5,2) * t101 + Icges(5,6) * t322;
t66 = Icges(5,1) * t102 + Icges(5,4) * t101 + Icges(5,5) * t322;
t360 = Icges(5,4) * t192;
t82 = Icges(5,2) * t191 + Icges(5,6) * t351 + t360;
t167 = Icges(5,4) * t191;
t84 = Icges(5,1) * t192 + Icges(5,5) * t351 + t167;
t19 = t103 * t82 + t104 * t84 + t193 * t64 + t194 * t66 + t221 * t265;
t357 = Icges(5,4) * t224;
t280 = -Icges(5,2) * t222 + t357;
t142 = -Icges(5,6) * t215 + t214 * t280;
t358 = Icges(5,4) * t222;
t286 = Icges(5,1) * t224 - t358;
t144 = -Icges(5,5) * t215 + t214 * t286;
t276 = Icges(5,5) * t224 - Icges(5,6) * t222;
t140 = -Icges(5,3) * t215 + t214 * t276;
t253 = t276 * t215;
t275 = -Icges(5,5) * t222 - Icges(5,6) * t224;
t71 = t218 * t253 + (Icges(5,3) * t218 + qJD(4) * t275) * t214;
t261 = t140 * t349 + t214 * t71;
t41 = t193 * t82 + t194 * t84 + t350 * t80;
t367 = t220 * t41;
t42 = t193 * t83 + t194 * t85 + t350 * t81;
t295 = t221 * t42 + t367;
t254 = t280 * t215;
t279 = -Icges(5,2) * t224 - t358;
t72 = t218 * t254 + (Icges(5,6) * t218 + qJD(4) * t279) * t214;
t255 = t286 * t215;
t285 = -Icges(5,1) * t222 - t357;
t73 = t218 * t255 + (Icges(5,5) * t218 + qJD(4) * t285) * t214;
t245 = (-t103 * t142 - t104 * t144 - t193 * t72 - t194 * t73 + t218 * t295 - t221 * t261) * t215;
t57 = t140 * t350 + t142 * t193 + t144 * t194;
t399 = t405 + t19 * t378 + (t352 * t57 + t245) * qJD(4) / 0.2e1;
t277 = -Icges(4,5) * t214 - Icges(4,6) * t215;
t169 = t277 * t220;
t170 = t277 * t221;
t398 = t169 * t347 / 0.2e1 - t170 * t206 / 0.2e1;
t216 = t220 ^ 2;
t217 = t221 ^ 2;
t335 = t216 + t217;
t128 = t218 * t169;
t129 = t218 * t170;
t393 = -Icges(4,5) * t220 + t401 * t221;
t395 = -Icges(4,6) * t220 + t400 * t221;
t237 = (t214 * t393 + t215 * t395) * t218;
t394 = Icges(4,5) * t221 + t401 * t220;
t396 = Icges(4,6) * t221 + t400 * t220;
t238 = (t214 * t394 + t215 * t396) * t218;
t397 = -(t237 + t128) * t220 / 0.2e1 + (t129 / 0.2e1 - t238 / 0.2e1) * t221;
t223 = sin(qJ(2));
t225 = cos(qJ(2));
t363 = Icges(3,4) * t225;
t364 = Icges(3,4) * t223;
t392 = (-t223 * (-Icges(3,2) * t225 - t364) + t225 * (-Icges(3,1) * t223 - t363)) * qJD(2);
t39 = t191 * t82 + t192 * t84 + t351 * t80;
t40 = t191 * t83 + t192 * t85 + t351 * t81;
t56 = t140 * t351 + t142 * t191 + t144 * t192;
t391 = (t161 * t42 + t162 * t41 - t332 * t57) * t221 + (t161 * t40 + t162 * t39 - t332 * t56) * t220;
t114 = t146 * t220;
t203 = pkin(3) * t214 - pkin(6) * t215;
t260 = t203 * t220;
t136 = t218 * t260;
t177 = t203 * t221;
t137 = t218 * t177;
t320 = t215 * t330;
t204 = pkin(3) * t215 + pkin(6) * t214;
t176 = t204 * t220;
t178 = t204 * t221;
t374 = pkin(2) * t225;
t138 = -pkin(5) * t221 + t220 * t374;
t139 = pkin(5) * t220 + t221 * t374;
t334 = qJD(2) * t221;
t318 = t138 * t213 + t139 * t334 + qJD(1);
t87 = rSges(5,1) * t194 + rSges(5,2) * t193 + rSges(5,3) * t350;
t38 = t161 * t86 - t162 * t87 + t176 * t206 + t178 * t347 + t318;
t68 = rSges(5,1) * t102 + rSges(5,2) * t101 + rSges(5,3) * t322;
t69 = rSges(5,1) * t104 + rSges(5,2) * t103 + rSges(5,3) * t321;
t390 = (t161 * t114 + t177 * t347 + t206 * t260 + t320 * t87 + (-t136 + t68) * t220 + (-t137 + t69 + t404) * t221) * t38;
t201 = rSges(4,1) * t214 + rSges(4,2) * t215;
t258 = t220 * t201;
t134 = t218 * t258;
t175 = t201 * t221;
t135 = t218 * t175;
t202 = rSges(4,1) * t215 - rSges(4,2) * t214;
t153 = -rSges(4,3) * t221 + t202 * t220;
t154 = rSges(4,3) * t220 + t202 * t221;
t389 = (-t220 * t134 - t221 * t135 + t175 * t347 + t206 * t258) * (t153 * t206 + t154 * t347 + t318);
t208 = rSges(3,1) * t223 + rSges(3,2) * t225;
t256 = qJD(2) * t208;
t291 = -t222 * t83 + t224 * t85;
t292 = -t222 * t82 + t224 * t84;
t383 = -(-t140 * t221 - t291) * t161 - (-t140 * t220 - t292) * t162;
t231 = t161 * (-Icges(5,2) * t194 + t168 + t85) + t162 * (-Icges(5,2) * t192 + t167 + t84) - t332 * (t279 * t214 + t144);
t227 = qJD(2) ^ 2;
t381 = -t161 / 0.2e1;
t379 = -t162 / 0.2e1;
t375 = pkin(2) * t223;
t368 = pkin(2) * qJD(2);
t366 = t221 * t40;
t166 = t204 * t218;
t257 = t299 * t215;
t298 = -rSges(5,1) * t222 - rSges(5,2) * t224;
t74 = t218 * t257 + (rSges(5,3) * t218 + qJD(4) * t298) * t214;
t365 = -t166 - t74;
t353 = t140 * t215;
t341 = t220 * t138 + t221 * t139;
t340 = t220 * t153 + t221 * t154;
t338 = -t146 - t203;
t331 = qJD(4) * t218;
t329 = t227 * t374;
t328 = t223 * t368;
t327 = t225 * t368;
t317 = -t332 / 0.2e1;
t316 = t332 / 0.2e1;
t315 = t331 / 0.2e1;
t314 = -t201 - t375;
t312 = t220 * t329;
t311 = t221 * t329;
t309 = (t178 + t87) * t221 + (t176 + t86) * t220;
t308 = t220 * t328;
t307 = t220 * t327;
t306 = t221 * t328;
t305 = t221 * t327;
t303 = t214 * t315;
t302 = t215 * t315;
t301 = t338 - t375;
t163 = t202 * t218;
t300 = -t163 - t327;
t209 = rSges(3,1) * t225 - rSges(3,2) * t223;
t297 = t81 * t161 + t80 * t162;
t296 = t220 * t39 + t366;
t47 = t214 * t292 - t215 * t80;
t48 = t214 * t291 - t215 * t81;
t294 = t47 * t220 + t48 * t221;
t293 = -t220 * t87 + t221 * t86;
t290 = Icges(3,1) * t225 - t364;
t284 = -Icges(3,2) * t223 + t363;
t278 = -Icges(3,5) * t223 - Icges(3,6) * t225;
t274 = -t142 * t222 + t144 * t224;
t273 = t335 * t209;
t271 = t335 * t256;
t270 = t227 * t335 * t375;
t269 = t220 * t302;
t268 = t221 * t302;
t267 = -t327 + t365;
t259 = Icges(5,3) * t214 + t253 - t274;
t250 = qJD(2) * t278;
t147 = rSges(5,3) * t214 + t257;
t249 = -t114 * t332 + t146 * t320 + t162 * t147 - t204 * t347 - t333 * t86;
t247 = -t161 * (Icges(5,5) * t193 - Icges(5,6) * t194) - t162 * (Icges(5,5) * t191 - Icges(5,6) * t192) + t275 * t214 * t332;
t246 = (-t101 * t142 - t102 * t144 - t191 * t72 - t192 * t73 + t218 * t296 - t220 * t261) * t215;
t244 = (t218 * t294 - (t218 * t274 - t71) * t215 - (t140 * t218 - t222 * t72 + t224 * t73 + (-t142 * t224 - t144 * t222) * qJD(4)) * t214) * t215;
t239 = t214 * t247;
t235 = (-(-Icges(3,6) * t221 + t220 * t284) * t225 - (-Icges(3,5) * t221 + t220 * t290) * t223) * qJD(2) + t392 * t220;
t234 = (-(Icges(3,6) * t220 + t221 * t284) * t225 - (Icges(3,5) * t220 + t221 * t290) * t223) * qJD(2) + t392 * t221;
t233 = -t147 * t161 - t204 * t206 + t87 * t333;
t232 = (Icges(5,1) * t193 - t359 - t83) * t161 + (Icges(5,1) * t191 - t360 - t82) * t162 - (t285 * t214 - t142) * t332;
t108 = t142 * t220;
t109 = t142 * t221;
t11 = (t218 * t292 - t62) * t215 + (t218 * t80 - t222 * t64 + t224 * t66 + (-t222 * t84 - t224 * t82) * qJD(4)) * t214;
t110 = t144 * t220;
t111 = t144 * t221;
t12 = (t218 * t291 - t63) * t215 + (t218 * t81 - t222 * t65 + t224 * t67 + (-t222 * t85 - t224 * t83) * qJD(4)) * t214;
t143 = Icges(5,6) * t214 + t254;
t145 = Icges(5,5) * t214 + t255;
t17 = t101 * t82 + t102 * t84 + t191 * t64 + t192 * t66 + t220 * t265;
t18 = t101 * t83 + t102 * t85 + t191 * t65 + t192 * t67 + t220 * t264;
t58 = t214 * t274 - t353;
t22 = t161 * t48 + t162 * t47 - t332 * t58;
t228 = (-t259 * t332 - t383) * t214;
t229 = (t206 * t395 - t347 * t396) * t215 + (t206 * t393 - t347 * t394) * t214;
t3 = t161 * t18 + t162 * t17 + (t352 * t56 + t246) * qJD(4);
t230 = -t22 * t333 / 0.2e1 + ((-t109 * t193 - t111 * t194) * t161 + (-t108 * t193 - t110 * t194) * t162 + (t57 * t214 + (-t143 * t193 - t145 * t194 + t367) * t215) * qJD(4)) * t381 + ((-t109 * t191 - t111 * t192) * t161 + (-t108 * t191 - t110 * t192) * t162 + (t56 * t214 + (-t143 * t191 - t145 * t192 + t366) * t215) * qJD(4)) * t379 + (((t109 * t222 - t111 * t224 + t81) * t161 + (t108 * t222 - t110 * t224 + t80) * t162 + t58 * qJD(4)) * t214 + ((t259 * t215 + (t143 * t222 - t145 * t224 - t140) * t214 + t294) * qJD(4) + t383) * t215) * t316 + t391 * t317 + (t18 * t378 + t40 * t269 + t42 * t268 + t48 * t303 + (((t39 - t353) * qJD(4) + t297) * t215 + t228) * t379 + t405 + t399 + t12 * t317 + (t229 / 0.2e1 + t397) * t347) * t220 + (-t17 * t378 - t39 * t269 - t41 * t268 - t47 * t303 + (((t42 - t353) * qJD(4) + t297) * t215 + t228) * t381 - t19 * t380 - t3 / 0.2e1 - t11 * t317 + (-t128 * t221 + t220 * t238 + t398) * t347) * t221 + ((t129 * t220 + t221 * t237 + t398) * t220 + (-t229 / 0.2e1 + t397) * t221) * t206;
t196 = t278 * t221;
t195 = t278 * t220;
t188 = t298 * t214;
t183 = t221 * t250;
t182 = t220 * t250;
t113 = -t201 * t347 - t306;
t112 = -t201 * t206 - t308;
t100 = rSges(5,1) * t193 - rSges(5,2) * t194;
t99 = rSges(5,1) * t191 - rSges(5,2) * t192;
t91 = -t163 * t347 - t311;
t90 = -t163 * t206 - t312;
t88 = t271 * qJD(2);
t79 = qJD(2) * t273 + qJD(1);
t59 = -t134 * t206 - t135 * t347 - t270;
t50 = -t203 * t347 - t306 - t404;
t49 = -t146 * t161 - t203 * t206 - t332 * t87 - t308;
t37 = -t311 + t162 * t74 - t166 * t347 + (-t86 * t352 + (t146 * t348 + t68) * t215) * qJD(4);
t36 = -t312 - t161 * t74 - t166 * t206 + (t87 * t352 + (-t146 * t347 - t69) * t215) * qJD(4);
t23 = t215 * t293 * t331 - t136 * t206 - t137 * t347 + t161 * t68 - t162 * t69 - t270;
t1 = [-m(3) * t88 + m(4) * t59 + m(5) * t23; -(t220 * (-t183 * t221 + t220 * t234) - t221 * (-t182 * t221 + t220 * t235)) * t334 + (t220 * (t183 * t220 + t221 * t234) - t221 * (t182 * t220 + t221 * t235)) * t213 + (t195 * qJD(2) * t217 - t221 * t196 * t213) * t334 / 0.2e1 - (t196 * qJD(2) * t216 - t220 * t195 * t334) * t213 / 0.2e1 + t230 + (t23 * (t309 + t341) + (t267 * t50 + t301 * t37) * t221 + (t267 * t49 + t301 * t36) * t220 - t50 * (t249 - t305) - t49 * (t233 - t307) + t390) * m(5) + (-t113 * (-t202 * t347 - t305) - t112 * (-t202 * t206 - t307) + t59 * (t340 + t341) + (t113 * t300 + t314 * t91) * t221 + (t112 * t300 + t314 * t90) * t220 + t389) * m(4) + (-t271 * t79 - t273 * t88 + (t208 * t209 * t227 + t256 * t79) * t335) * m(3); t230 + (-t233 * t49 - t249 * t50 + t23 * t309 + (t338 * t37 + t365 * t50) * t221 + (t338 * t36 + t365 * t49) * t220 + t390) * m(5) + (-(-t112 * t206 - t113 * t347) * t202 + t59 * t340 + (-t220 * t90 - t221 * t91) * t201 + (-t112 * t220 - t113 * t221) * t163 + t389) * m(4); t350 * t399 + (t214 * t295 - t215 * t57) * t268 + (t245 + (t19 * t220 + t20 * t221 + t218 * t57) * t214) * t380 + t3 * t351 / 0.2e1 + (t214 * t296 - t215 * t56) * t269 + (t246 + (t17 * t220 + t18 * t221 + t218 * t56) * t214) * t378 + t22 * t352 / 0.2e1 - t215 * (t11 * t162 + t12 * t161 + (t352 * t58 + t244) * qJD(4)) / 0.2e1 + (t214 * t294 - t215 * t58) * t303 + (t244 + (t11 * t220 + t12 * t221 + t218 * t58) * t214) * t317 + (t193 * t231 + t194 * t232 - t221 * t239) * t381 + (t191 * t231 + t192 * t232 - t220 * t239) * t379 + (t247 * t215 + (-t222 * t231 + t232 * t224) * t214) * t316 + t391 * t349 / 0.2e1 + ((-t36 * t87 + t37 * t86 - t49 * t69 + t50 * t68 + (t38 * t293 + (t220 * t50 - t221 * t49) * t146) * t218) * t215 + (t50 * (-t218 * t86 + t220 * t74) + t49 * (t218 * t87 - t221 * t74) + t23 * t293 + t38 * (-t220 * t69 + t221 * t68) + (t220 * t37 - t221 * t36) * t146) * t214 - t50 * (t162 * t188 + t332 * t99) - t49 * (-t100 * t332 - t161 * t188) - t38 * (-t100 * t162 + t161 * t99)) * m(5);];
tauc = t1(:);
