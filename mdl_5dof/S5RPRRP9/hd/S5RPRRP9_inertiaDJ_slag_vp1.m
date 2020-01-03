% Calculate time derivative of joint inertia matrix for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP9_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP9_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP9_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:37
% EndTime: 2019-12-31 18:48:53
% DurationCPUTime: 7.99s
% Computational Cost: add. (9270->471), mult. (8481->649), div. (0->0), fcn. (6594->8), ass. (0->266)
t200 = pkin(8) + qJ(3);
t192 = qJ(4) + t200;
t185 = cos(t192);
t184 = sin(t192);
t329 = Icges(6,5) * t184;
t154 = -Icges(6,3) * t185 + t329;
t334 = Icges(5,4) * t184;
t157 = Icges(5,2) * t185 + t334;
t384 = -t157 + t154;
t328 = Icges(6,5) * t185;
t158 = Icges(6,1) * t184 - t328;
t333 = Icges(5,4) * t185;
t159 = Icges(5,1) * t184 + t333;
t383 = t159 + t158;
t348 = rSges(6,1) + pkin(4);
t382 = rSges(6,3) + qJ(5);
t207 = sin(qJ(1));
t208 = cos(qJ(1));
t241 = Icges(6,3) * t184 + t328;
t112 = Icges(6,6) * t207 + t241 * t208;
t245 = -Icges(5,2) * t184 + t333;
t118 = Icges(5,6) * t207 + t245 * t208;
t381 = t112 - t118;
t242 = Icges(5,5) * t185 - Icges(5,6) * t184;
t114 = Icges(5,3) * t207 + t242 * t208;
t244 = Icges(6,4) * t185 + Icges(6,6) * t184;
t116 = Icges(6,2) * t207 + t244 * t208;
t380 = t114 + t116;
t248 = Icges(6,1) * t185 + t329;
t120 = Icges(6,4) * t207 + t248 * t208;
t249 = Icges(5,1) * t185 - t334;
t122 = Icges(5,5) * t207 + t249 * t208;
t379 = t120 + t122;
t201 = qJD(3) + qJD(4);
t378 = (-t241 + t245) * t201;
t377 = (t248 + t249) * t201;
t155 = Icges(5,5) * t184 + Icges(5,6) * t185;
t156 = Icges(6,4) * t184 - Icges(6,6) * t185;
t371 = t155 + t156;
t368 = t384 * t184 + t383 * t185;
t111 = -Icges(6,6) * t208 + t241 * t207;
t117 = -Icges(5,6) * t208 + t245 * t207;
t376 = -t111 + t117;
t119 = -Icges(6,4) * t208 + t248 * t207;
t121 = -Icges(5,5) * t208 + t249 * t207;
t375 = t119 + t121;
t113 = -Icges(5,3) * t208 + t242 * t207;
t115 = -Icges(6,2) * t208 + t244 * t207;
t240 = t111 * t184 + t119 * t185;
t359 = t240 * t208;
t238 = t117 * t184 - t121 * t185;
t361 = t238 * t208;
t374 = -t359 + t361 + (-t113 - t115) * t207;
t237 = t118 * t184 - t122 * t185;
t239 = t112 * t184 + t120 * t185;
t373 = (-t237 + t239) * t208 + t380 * t207;
t198 = t207 * rSges(6,2);
t309 = t185 * t208;
t311 = t184 * t208;
t372 = t348 * t309 + t382 * t311 + t198;
t161 = rSges(6,1) * t184 - rSges(6,3) * t185;
t370 = pkin(4) * t184 - qJ(5) * t185 + t161;
t313 = t159 * t201;
t314 = t158 * t201;
t315 = t157 * t201;
t316 = t154 * t201;
t369 = (t316 - t315 + t377) * t185 + (-t314 - t313 - t378) * t184 + t371 * qJD(1);
t65 = -qJD(1) * t111 - t208 * t316;
t71 = -qJD(1) * t117 - t208 * t315;
t367 = -t379 * t201 + t65 - t71;
t73 = -qJD(1) * t119 - t208 * t314;
t75 = -qJD(1) * t121 - t208 * t313;
t366 = t381 * t201 + t73 + t75;
t365 = (-t242 - t244) * t201 + t368 * qJD(1);
t190 = sin(t200);
t191 = cos(t200);
t335 = Icges(4,4) * t191;
t247 = -Icges(4,2) * t190 + t335;
t133 = Icges(4,6) * t207 + t247 * t208;
t336 = Icges(4,4) * t190;
t251 = Icges(4,1) * t191 - t336;
t135 = Icges(4,5) * t207 + t251 * t208;
t235 = t133 * t190 - t135 * t191;
t364 = t235 * t207;
t132 = -Icges(4,6) * t208 + t247 * t207;
t134 = -Icges(4,5) * t208 + t251 * t207;
t236 = t132 * t190 - t134 * t191;
t363 = t236 * t208;
t362 = t237 * t207;
t360 = t239 * t207;
t206 = -pkin(6) - qJ(2);
t205 = cos(pkin(8));
t186 = t205 * pkin(2) + pkin(1);
t341 = rSges(4,2) * t190;
t343 = rSges(4,1) * t191;
t260 = -t341 + t343;
t229 = -t186 - t260;
t96 = (rSges(4,3) - t206) * t208 + t229 * t207;
t197 = t207 * rSges(4,3);
t296 = t208 * t343 + t197;
t97 = -t207 * t206 + (t186 - t341) * t208 + t296;
t358 = t207 * t97 + t208 * t96;
t357 = qJD(1) * t113;
t356 = qJD(1) * t115;
t243 = Icges(4,5) * t191 - Icges(4,6) * t190;
t130 = -Icges(4,3) * t208 + t243 * t207;
t230 = rSges(3,1) * t205 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t338 = rSges(3,3) + qJ(2);
t110 = t338 * t207 + t230 * t208;
t199 = -pkin(7) + t206;
t295 = t199 - t206;
t170 = pkin(3) * t191 + t186;
t299 = t170 - t186;
t98 = t299 * t207 + t295 * t208;
t353 = 2 * m(4);
t352 = 2 * m(5);
t351 = 2 * m(6);
t202 = t207 ^ 2;
t203 = t208 ^ 2;
t350 = t207 / 0.2e1;
t349 = -t208 / 0.2e1;
t168 = rSges(4,1) * t190 + rSges(4,2) * t191;
t347 = m(4) * t168;
t162 = rSges(5,1) * t184 + rSges(5,2) * t185;
t346 = m(5) * t162;
t345 = pkin(3) * t190;
t163 = t208 * t170;
t99 = -t186 * t208 - t295 * t207 + t163;
t344 = t207 * t98 + t208 * t99;
t342 = rSges(5,1) * t185;
t196 = t207 * rSges(5,3);
t273 = -t162 - t345;
t104 = t273 * t208;
t321 = t104 * t208;
t320 = t132 * t191;
t319 = t133 * t191;
t318 = t134 * t190;
t317 = t135 * t190;
t312 = t184 * t201;
t310 = t185 * t201;
t308 = t201 * t207;
t307 = t201 * t208;
t259 = -rSges(5,2) * t184 + t342;
t145 = t259 * t201;
t306 = t207 * t145;
t258 = rSges(6,1) * t185 + rSges(6,3) * t184;
t272 = pkin(4) * t201 - qJD(5);
t305 = -qJ(5) * t312 - t272 * t185 - t258 * t201;
t125 = -rSges(5,3) * t208 + t259 * t207;
t127 = rSges(5,1) * t309 - rSges(5,2) * t311 + t196;
t64 = t207 * t125 + t208 * t127;
t290 = qJD(1) * t207;
t303 = t370 * t290;
t279 = t185 * t307;
t285 = qJD(5) * t184;
t301 = qJ(5) * t279 + t208 * t285;
t289 = qJD(1) * t208;
t300 = rSges(6,2) * t289 + rSges(6,3) * t279;
t277 = t184 * t290;
t298 = rSges(5,2) * t277 + rSges(5,3) * t289;
t288 = qJD(3) * t190;
t281 = pkin(3) * t288;
t297 = -t199 * t290 - t207 * t281;
t294 = t202 + t203;
t293 = qJD(1) * t114;
t292 = qJD(1) * t116;
t131 = Icges(4,3) * t207 + t243 * t208;
t291 = qJD(1) * t131;
t287 = qJD(3) * t191;
t286 = qJD(3) * t207;
t183 = t206 * t290;
t284 = t207 * (t299 * t289 + t183 + t297) + t208 * (-qJD(1) * t98 - t208 * t281) + t98 * t289;
t216 = -t184 * t307 - t185 * t290;
t227 = t162 * t201;
t283 = t125 * t289 + t207 * (-t207 * t227 + (t259 * t208 + t196) * qJD(1)) + t208 * (t216 * rSges(5,1) - rSges(5,2) * t279 + t298);
t282 = t208 * t341;
t280 = pkin(3) * t287;
t194 = qJD(2) * t208;
t278 = t194 - t297;
t276 = t190 * t290;
t74 = t120 * qJD(1) - t207 * t314;
t271 = t111 * t201 + t74;
t76 = t122 * qJD(1) - t207 * t313;
t269 = t117 * t201 - t76;
t66 = t112 * qJD(1) - t207 * t316;
t267 = t119 * t201 - t66;
t72 = t118 * qJD(1) - t207 * t315;
t265 = t121 * t201 + t72;
t95 = t370 * t208;
t263 = -t207 * t199 + t163;
t124 = -rSges(6,2) * t208 + t258 * t207;
t146 = (pkin(4) * t185 + qJ(5) * t184) * t207;
t27 = t372 * t208 + (t124 + t146) * t207;
t262 = -t370 - t345;
t30 = -t115 * t208 + t240 * t207;
t31 = -t116 * t208 + t360;
t32 = -t113 * t208 - t238 * t207;
t33 = -t114 * t208 - t362;
t221 = t201 * t155;
t67 = -t208 * t221 - t357;
t68 = -t207 * t221 + t293;
t222 = t201 * t156;
t69 = -t208 * t222 - t356;
t70 = -t207 * t222 + t292;
t261 = ((-t30 - t32) * t290 + t374 * t289) * t208 + (((t69 + t67) * t207 + (-t360 + t362 - t374) * qJD(1)) * t207 + (t31 + t33) * t290 + t373 * t289 + ((-t68 - t70) * t207 + (t376 * t310 + t375 * t312 - t356 - t357) * t208 + (t366 * t207 + (-t74 - t76) * t208) * t185 + (t367 * t207 + (-t66 + t72) * t208) * t184 + ((t240 - t238 + t380) * t207 + t373) * qJD(1)) * t208) * t207;
t214 = -t184 * t382 - t348 * t185 - t170;
t212 = t214 * t207;
t54 = (rSges(6,2) - t199) * t208 + t212;
t55 = t263 + t372;
t255 = t207 * t55 + t208 * t54;
t85 = t262 * t207;
t86 = t262 * t208;
t254 = t207 * t85 + t208 * t86;
t228 = -t170 - t259;
t87 = (rSges(5,3) - t199) * t208 + t228 * t207;
t88 = t127 + t263;
t253 = t207 * t88 + t208 * t87;
t94 = t370 * t207;
t252 = -t207 * t94 - t208 * t95;
t250 = Icges(4,1) * t190 + t335;
t246 = Icges(4,2) * t191 + t336;
t232 = t124 * t289 + t146 * t289 + t207 * ((pkin(4) * t289 + qJ(5) * t308) * t185 + (qJ(5) * t289 - t272 * t207) * t184) + t208 * (t216 * pkin(4) - qJ(5) * t277 + t301) + t207 * (-t161 * t308 + (t258 * t208 + t198) * qJD(1)) + t208 * (t216 * rSges(6,1) - rSges(6,3) * t277 + t300);
t231 = -t280 + t305;
t226 = qJD(3) * t168;
t219 = qJD(3) * t250;
t218 = qJD(3) * t246;
t217 = qJD(3) * (-Icges(4,5) * t190 - Icges(4,6) * t191);
t3 = (t208 * t70 + (t31 - t359) * qJD(1)) * t208 + (t30 * qJD(1) + (t112 * t310 - t120 * t312 + t184 * t65 + t185 * t73 + t292) * t207 + (-t69 - t271 * t185 + t267 * t184 + (-t115 + t239) * qJD(1)) * t208) * t207;
t4 = (t208 * t68 + (t33 + t361) * qJD(1)) * t208 + (t32 * qJD(1) + (-t118 * t310 - t122 * t312 - t184 * t71 + t185 * t75 + t293) * t207 + (-t67 + t269 * t185 + t265 * t184 + (-t113 - t237) * qJD(1)) * t208) * t207;
t215 = (-t4 - t3) * t208 + t261;
t213 = rSges(4,2) * t276 + rSges(4,3) * t289 - t208 * t226;
t109 = -t230 * t207 + t338 * t208;
t211 = (t366 * t184 - t367 * t185 - t365 * t207 + t369 * t208) * t350 + (t365 * t208 + t369 * t207 + (t265 + t267) * t185 + (-t269 + t271) * t184) * t349 + (t375 * t184 + t376 * t185 + t368 * t207 - t371 * t208) * t290 / 0.2e1 + (t379 * t184 - t381 * t185 + t371 * t207 + t368 * t208) * t289 / 0.2e1;
t193 = qJD(2) * t207;
t176 = pkin(3) * t276;
t153 = t260 * qJD(3);
t143 = -t282 + t296;
t142 = -rSges(4,3) * t208 + t260 * t207;
t103 = t273 * t207;
t93 = -qJD(1) * t110 + t194;
t92 = t109 * qJD(1) + t193;
t80 = t207 * t217 + t291;
t79 = -qJD(1) * t130 + t208 * t217;
t57 = -t162 * t289 - t306 + (-t190 * t289 - t191 * t286) * pkin(3);
t56 = t162 * t290 + t176 + (-t145 - t280) * t208;
t53 = t183 + t194 + t168 * t286 + (t229 * t208 - t197) * qJD(1);
t52 = t193 + (-t206 * t208 + (-t186 - t343) * t207) * qJD(1) + t213;
t43 = t207 * t131 - t235 * t208;
t42 = t207 * t130 - t363;
t41 = -t131 * t208 - t364;
t40 = -t130 * t208 - t236 * t207;
t39 = t162 * t308 + (t228 * t208 - t196) * qJD(1) + t278;
t38 = t193 + (-t170 - t342) * t290 + (-qJD(1) * t199 - t227 - t281) * t208 + t298;
t29 = -qJD(1) * t95 + t305 * t207;
t28 = t305 * t208 + t303;
t26 = qJD(1) * t86 + t231 * t207;
t25 = t231 * t208 + t176 + t303;
t24 = t64 + t344;
t23 = (-t285 + (t348 * t184 - t185 * t382) * t201) * t207 + (t214 * t208 - t198) * qJD(1) + t278;
t22 = t193 + (-t348 * t312 - t281) * t208 + (-t199 * t208 + t212) * qJD(1) + t300 + t301;
t21 = -t127 * t290 + t283;
t20 = t27 + t344;
t7 = -t290 * t372 + t232;
t6 = (-t127 - t99) * t290 + t283 + t284;
t5 = (-t99 - t372) * t290 + t232 + t284;
t1 = [(t38 * t88 + t39 * t87) * t352 + (t22 * t55 + t23 * t54) * t351 + (t52 * t97 + t53 * t96) * t353 + 0.2e1 * m(3) * (t109 * t93 + t110 * t92) + t384 * t312 + t383 * t310 + (t251 - t246) * t288 + (t250 + t247) * t287 + t378 * t185 + t377 * t184; m(5) * (t253 * qJD(1) + t207 * t39 - t208 * t38) + m(6) * (t255 * qJD(1) + t207 * t23 - t208 * t22) + m(4) * (qJD(1) * t358 + t207 * t53 - t208 * t52) + m(3) * (t207 * t93 - t208 * t92 + (t109 * t208 + t110 * t207) * qJD(1)); 0; (-t235 * qJD(3) + t190 * (-qJD(1) * t134 - t208 * t219) + t191 * (-qJD(1) * t132 - t208 * t218)) * t350 + (-t236 * qJD(3) + t190 * (t135 * qJD(1) - t207 * t219) + t191 * (t133 * qJD(1) - t207 * t218)) * t349 + (t202 / 0.2e1 + t203 / 0.2e1) * t243 * qJD(3) + m(6) * (t22 * t85 + t23 * t86 + t25 * t54 + t26 * t55) + m(5) * (t103 * t38 + t104 * t39 + t56 * t87 + t57 * t88) + ((t319 / 0.2e1 + t317 / 0.2e1 - t97 * t347) * t208 + (t320 / 0.2e1 + t318 / 0.2e1 + t96 * t347) * t207) * qJD(1) + t211 + m(4) * ((-t207 * t52 - t208 * t53) * t168 - t358 * t153); m(5) * (t56 * t207 - t208 * t57 + (t103 * t207 + t321) * qJD(1)) + m(6) * (t254 * qJD(1) + t25 * t207 - t208 * t26); (t20 * t5 + t25 * t86 + t26 * t85) * t351 - t208 * t4 + (t103 * t57 + t104 * t56 + t24 * t6) * t352 - t208 * t3 + ((t207 * t142 + t143 * t208) * ((qJD(1) * t142 + t213) * t208 + (-t207 * t226 + (-t143 - t282 + t197) * qJD(1)) * t207) + t294 * t168 * t153) * t353 + (t41 * t207 - t208 * t40) * t290 - t208 * ((t208 * t80 + (t41 + t363) * qJD(1)) * t208 + (t40 * qJD(1) + (-t133 * t287 - t135 * t288 + t291) * t207 + (-t79 + (t318 + t320) * qJD(3) - t235 * qJD(1)) * t208) * t207) + (t43 * t207 - t208 * t42) * t289 + t207 * ((t207 * t79 + (t42 + t364) * qJD(1)) * t207 + (t43 * qJD(1) + (t132 * t287 + t134 * t288) * t208 + (-t80 + (-t317 - t319) * qJD(3) + (t131 - t236) * qJD(1)) * t207) * t208) + t261; m(6) * (-t22 * t94 - t23 * t95 + t28 * t54 + t29 * t55) + (-t207 * t38 - t208 * t39 + (t207 * t87 - t208 * t88) * qJD(1)) * t346 - m(5) * t253 * t145 + t211; m(6) * (t252 * qJD(1) + t28 * t207 - t208 * t29); m(6) * (t20 * t7 - t25 * t95 - t26 * t94 + t27 * t5 + t28 * t86 + t29 * t85) + m(5) * (-t103 * t306 - t145 * t321 + t21 * t24 + t64 * t6) + (-t207 * t57 - t208 * t56 + (-t103 * t208 + t104 * t207) * qJD(1)) * t346 + t215; (t294 * t162 * t145 + t21 * t64) * t352 + (t27 * t7 - t28 * t95 - t29 * t94) * t351 + t215; m(6) * (t255 * t310 + (t207 * t22 + t208 * t23 + (-t207 * t54 + t208 * t55) * qJD(1)) * t184); 0; m(6) * ((t254 * t201 - t5) * t185 + (t20 * t201 + t207 * t26 + t208 * t25 + (-t207 * t86 + t208 * t85) * qJD(1)) * t184); m(6) * ((t252 * t201 - t7) * t185 + (t201 * t27 + t207 * t29 + t208 * t28 + (t207 * t95 - t208 * t94) * qJD(1)) * t184); (-0.1e1 + t294) * t184 * t310 * t351;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
