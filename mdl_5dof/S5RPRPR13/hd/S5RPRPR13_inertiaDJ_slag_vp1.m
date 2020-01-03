% Calculate time derivative of joint inertia matrix for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR13_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR13_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR13_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:58
% EndTime: 2019-12-31 18:32:17
% DurationCPUTime: 11.07s
% Computational Cost: add. (10449->578), mult. (15248->842), div. (0->0), fcn. (14137->8), ass. (0->287)
t194 = pkin(8) + qJ(3);
t188 = cos(t194);
t187 = sin(t194);
t324 = Icges(5,6) * t187;
t332 = Icges(4,4) * t187;
t373 = -t324 - t332 + (-Icges(4,2) - Icges(5,3)) * t188;
t323 = Icges(5,6) * t188;
t331 = Icges(4,4) * t188;
t372 = -t323 - t331 + (-Icges(4,1) - Icges(5,2)) * t187;
t371 = t373 * qJD(3);
t370 = t372 * qJD(3);
t201 = sin(qJ(1));
t203 = cos(qJ(1));
t239 = -Icges(4,2) * t187 + t331;
t122 = -Icges(4,6) * t203 + t201 * t239;
t231 = -Icges(5,3) * t187 + t323;
t356 = Icges(5,5) * t203 + t201 * t231;
t369 = -t122 - t356;
t123 = Icges(4,6) * t201 + t203 * t239;
t126 = Icges(5,5) * t201 - t203 * t231;
t368 = t123 - t126;
t242 = Icges(4,1) * t188 - t332;
t124 = -Icges(4,5) * t203 + t201 * t242;
t233 = Icges(5,2) * t188 - t324;
t355 = Icges(5,4) * t203 + t201 * t233;
t367 = t124 + t355;
t125 = Icges(4,5) * t201 + t203 * t242;
t128 = Icges(5,4) * t201 - t203 * t233;
t366 = t125 - t128;
t292 = qJD(3) * t203;
t275 = t187 * t292;
t298 = qJD(1) * t201;
t365 = t188 * t298 + t275;
t364 = qJD(3) / 0.2e1;
t226 = t126 * t187 - t128 * t188;
t363 = t201 * t226;
t227 = t123 * t187 - t125 * t188;
t362 = t201 * t227;
t225 = -t187 * t356 + t188 * t355;
t361 = t203 * t225;
t228 = t122 * t187 - t124 * t188;
t360 = t203 * t228;
t192 = t201 * rSges(5,1);
t314 = t188 * t203;
t359 = -rSges(5,2) * t314 + t192;
t191 = t201 * rSges(4,3);
t316 = t187 * t203;
t358 = -rSges(4,2) * t316 + t191;
t199 = -pkin(6) - qJ(2);
t198 = cos(pkin(8));
t183 = pkin(2) * t198 + pkin(1);
t338 = rSges(4,1) * t188;
t262 = -rSges(4,2) * t187 + t338;
t218 = -t183 - t262;
t106 = (rSges(4,3) - t199) * t203 + t218 * t201;
t135 = rSges(4,1) * t314 + t358;
t269 = t183 * t203 - t199 * t201;
t107 = t135 + t269;
t357 = t106 * t203 + t107 * t201;
t235 = Icges(4,5) * t188 - Icges(4,6) * t187;
t120 = -Icges(4,3) * t203 + t201 * t235;
t237 = Icges(5,4) * t188 - Icges(5,5) * t187;
t354 = Icges(5,1) * t203 + t201 * t237;
t220 = rSges(3,1) * t198 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t335 = rSges(3,3) + qJ(2);
t115 = t201 * t335 + t203 * t220;
t353 = 2 * m(4);
t352 = 2 * m(5);
t351 = 2 * m(6);
t350 = m(5) / 0.2e1;
t349 = m(6) / 0.2e1;
t348 = t187 / 0.2e1;
t347 = t201 / 0.2e1;
t346 = -t203 / 0.2e1;
t345 = t203 / 0.2e1;
t344 = rSges(6,3) + pkin(7);
t165 = rSges(4,1) * t187 + rSges(4,2) * t188;
t343 = m(4) * t165;
t342 = pkin(3) * t188;
t193 = t201 * pkin(4);
t341 = qJD(1) / 0.2e1;
t200 = sin(qJ(5));
t202 = cos(qJ(5));
t328 = Icges(6,4) * t202;
t240 = Icges(6,1) * t200 + t328;
t118 = Icges(6,5) * t187 - t188 * t240;
t234 = Icges(6,5) * t200 + Icges(6,6) * t202;
t116 = Icges(6,3) * t187 - t188 * t234;
t329 = Icges(6,4) * t200;
t236 = Icges(6,2) * t202 + t329;
t117 = Icges(6,6) * t187 - t188 * t236;
t290 = qJD(5) * t188;
t293 = qJD(3) * t202;
t295 = qJD(3) * t188;
t296 = qJD(3) * t187;
t317 = t118 * t200;
t75 = (-Icges(6,5) * t202 + Icges(6,6) * t200) * t290 + (Icges(6,3) * t188 + t187 * t234) * qJD(3);
t263 = t116 * t295 + t187 * t75 + t296 * t317 + (t187 * t293 + t200 * t290) * t117;
t77 = (-Icges(6,1) * t202 + t329) * t290 + (Icges(6,5) * t188 + t187 * t240) * qJD(3);
t336 = t200 * t77;
t48 = t116 * t187 + (-t117 * t202 - t317) * t188;
t76 = (Icges(6,2) * t200 - t328) * t290 + (Icges(6,6) * t188 + t187 * t236) * qJD(3);
t340 = ((-t336 + (-qJD(5) * t118 - t76) * t202) * t188 + t263) * t187 + t48 * t295;
t267 = qJD(1) * t187 + qJD(5);
t273 = t188 * t292;
t205 = -t201 * t267 + t273;
t268 = qJD(5) * t187 + qJD(1);
t222 = t268 * t200;
t80 = t202 * t205 - t203 * t222;
t223 = t202 * t268;
t81 = t200 * t205 + t203 * t223;
t339 = rSges(6,1) * t81 + rSges(6,2) * t80;
t337 = rSges(5,2) * t187;
t334 = -rSges(5,3) - qJ(4);
t321 = qJ(4) * t187;
t320 = qJ(4) * t188;
t315 = t188 * t201;
t313 = t199 * t203;
t312 = t200 * t203;
t311 = t201 * t200;
t310 = t201 * t202;
t309 = t202 * t203;
t145 = t187 * t309 - t311;
t146 = t187 * t312 + t310;
t101 = rSges(6,1) * t146 + rSges(6,2) * t145 + rSges(6,3) * t314;
t308 = pkin(7) * t314 + t101 + t193;
t147 = t187 * t310 + t312;
t148 = t187 * t311 - t309;
t261 = -t148 * rSges(6,1) - t147 * rSges(6,2);
t102 = rSges(6,3) * t315 - t261;
t307 = -pkin(4) * t203 + pkin(7) * t315 + t102;
t257 = t321 + t342;
t133 = qJD(3) * t257 - qJD(4) * t188;
t259 = -rSges(5,2) * t188 + rSges(5,3) * t187;
t306 = -qJD(3) * t259 - t133;
t140 = t257 * t201;
t141 = pkin(3) * t314 + qJ(4) * t316;
t305 = t140 * t201 + t141 * t203;
t163 = pkin(3) * t187 - t320;
t258 = rSges(5,3) * t188 + t337;
t304 = -t163 + t258;
t291 = qJD(4) * t187;
t303 = qJ(4) * t273 + t203 * t291;
t190 = qJD(2) * t203;
t302 = t199 * t298 + t190;
t301 = t201 ^ 2 + t203 ^ 2;
t121 = Icges(4,3) * t201 + t203 * t235;
t300 = qJD(1) * t121;
t130 = Icges(5,1) * t201 - t203 * t237;
t299 = qJD(1) * t130;
t297 = qJD(1) * t203;
t294 = qJD(3) * t201;
t289 = -pkin(3) - t344;
t16 = -t116 * t365 + t117 * t80 + t118 * t81 + t145 * t76 + t146 * t77 + t314 * t75;
t97 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t314;
t99 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t314;
t252 = t200 * t99 + t202 * t97;
t35 = Icges(6,5) * t81 + Icges(6,6) * t80 - Icges(6,3) * t365;
t37 = Icges(6,4) * t81 + Icges(6,2) * t80 - Icges(6,6) * t365;
t39 = Icges(6,1) * t81 + Icges(6,4) * t80 - Icges(6,5) * t365;
t95 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t314;
t9 = (qJD(3) * t252 + t35) * t187 + (qJD(3) * t95 - t200 * t39 - t202 * t37 + (t200 * t97 - t202 * t99) * qJD(5)) * t188;
t288 = t9 / 0.2e1 + t16 / 0.2e1;
t276 = t187 * t294;
t173 = pkin(3) * t276;
t274 = t188 * t294;
t277 = t188 * t297;
t279 = t187 * t298;
t287 = t140 * t297 + t201 * (pkin(3) * t277 + t201 * t291 - t173 + (t187 * t297 + t274) * qJ(4)) + t203 * (-pkin(3) * t365 - qJ(4) * t279 + t303);
t100 = Icges(6,1) * t148 + Icges(6,4) * t147 + Icges(6,5) * t315;
t98 = Icges(6,4) * t148 + Icges(6,2) * t147 + Icges(6,6) * t315;
t243 = t100 * t200 + t202 * t98;
t210 = -t276 + t277;
t221 = t267 * t203;
t78 = t202 * t221 + (t188 * t293 - t222) * t201;
t79 = t201 * t223 + (t221 + t274) * t200;
t34 = Icges(6,5) * t79 + Icges(6,6) * t78 + Icges(6,3) * t210;
t36 = Icges(6,4) * t79 + Icges(6,2) * t78 + Icges(6,6) * t210;
t38 = Icges(6,1) * t79 + Icges(6,4) * t78 + Icges(6,5) * t210;
t96 = Icges(6,5) * t148 + Icges(6,6) * t147 + Icges(6,3) * t315;
t10 = (qJD(3) * t243 + t34) * t187 + (qJD(3) * t96 - t200 * t38 - t202 * t36 + (-t100 * t202 + t200 * t98) * qJD(5)) * t188;
t15 = t116 * t210 + t117 * t78 + t118 * t79 + t147 * t76 + t148 * t77 + t315 * t75;
t284 = t10 / 0.2e1 + t15 / 0.2e1;
t30 = t187 * t95 - t188 * t252;
t44 = t116 * t314 + t117 * t145 + t118 * t146;
t283 = -t30 / 0.2e1 - t44 / 0.2e1;
t31 = t187 * t96 - t188 * t243;
t45 = t116 * t315 + t117 * t147 + t118 * t148;
t282 = t31 / 0.2e1 + t45 / 0.2e1;
t189 = qJD(2) * t201;
t281 = t189 + t303;
t280 = t173 + t302;
t272 = -t237 * qJD(3) / 0.2e1 + t235 * t364;
t271 = -t296 / 0.2e1;
t260 = rSges(6,1) * t200 + rSges(6,2) * t202;
t119 = rSges(6,3) * t187 - t188 * t260;
t270 = pkin(7) * t187 + t119;
t112 = t304 * t203;
t266 = rSges(5,1) * t297 + rSges(5,2) * t365 + rSges(5,3) * t273;
t265 = -t163 - t270;
t264 = t79 * rSges(6,1) + t78 * rSges(6,2);
t25 = t145 * t97 + t146 * t99 + t314 * t95;
t26 = t100 * t146 + t145 * t98 + t314 * t96;
t251 = t201 * t26 + t203 * t25;
t17 = t201 * t25 - t203 * t26;
t27 = t147 * t97 + t148 * t99 + t315 * t95;
t28 = t100 * t148 + t147 * t98 + t315 * t96;
t250 = t201 * t28 + t203 * t27;
t18 = t201 * t27 - t203 * t28;
t249 = t201 * t31 + t203 * t30;
t248 = t201 * t30 - t203 * t31;
t208 = t188 * t289 - t183 - t321;
t206 = t208 * t201;
t50 = (pkin(4) - t199) * t203 + t206 + t261;
t219 = t269 + t141;
t51 = t219 + t308;
t247 = t201 * t51 + t203 * t50;
t60 = -t102 * t187 + t119 * t315;
t61 = t101 * t187 - t119 * t314;
t246 = t201 * t61 + t203 * t60;
t244 = t187 * t334 - t183;
t207 = (rSges(5,2) - pkin(3)) * t188 + t244;
t67 = (rSges(5,1) - t199) * t203 + t207 * t201;
t136 = rSges(5,3) * t316 + t359;
t68 = t136 + t219;
t245 = t201 * t68 + t203 * t67;
t229 = t101 * t201 - t102 * t203;
t82 = (-rSges(6,1) * t202 + rSges(6,2) * t200) * t290 + (rSges(6,3) * t188 + t187 * t260) * qJD(3);
t224 = -pkin(7) * t295 - t133 - t82;
t70 = t265 * t203;
t217 = qJD(3) * t165;
t214 = qJD(3) * (Icges(5,4) * t187 + Icges(5,5) * t188);
t213 = qJD(3) * (-Icges(4,5) * t187 - Icges(4,6) * t188);
t204 = rSges(4,2) * t279 + rSges(4,3) * t297 - t203 * t217;
t114 = -t201 * t220 + t203 * t335;
t186 = pkin(4) * t297;
t156 = t262 * qJD(3);
t144 = t163 * t298;
t137 = -rSges(5,1) * t203 + t201 * t259;
t134 = -rSges(4,3) * t203 + t201 * t262;
t111 = t304 * t201;
t105 = -qJD(1) * t115 + t190;
t104 = qJD(1) * t114 + t189;
t94 = qJD(1) * t354 + t203 * t214;
t93 = t201 * t214 + t299;
t84 = t201 * t213 + t300;
t83 = -qJD(1) * t120 + t203 * t213;
t69 = t265 * t201;
t65 = t165 * t294 + (t203 * t218 - t191) * qJD(1) + t302;
t64 = t189 + (-t313 + (-t183 - t338) * t201) * qJD(1) + t204;
t63 = qJD(1) * t112 + t201 * t306;
t62 = t203 * t306 - t258 * t298 + t144;
t59 = t201 * t225 + t203 * t354;
t58 = -t130 * t203 + t363;
t57 = -t201 * t354 + t361;
t56 = t130 * t201 + t203 * t226;
t55 = t121 * t201 - t203 * t227;
t54 = t120 * t201 - t360;
t53 = -t121 * t203 - t362;
t52 = -t120 * t203 - t201 * t228;
t49 = t136 * t203 + t137 * t201 + t305;
t46 = t229 * t188;
t43 = (-t291 + (t188 * t334 - t337) * qJD(3)) * t201 + (t203 * t207 - t192) * qJD(1) + t280;
t42 = -pkin(3) * t275 + (-t313 + (t244 - t342) * t201) * qJD(1) + t266 + t281;
t41 = -rSges(6,3) * t365 + t339;
t40 = rSges(6,3) * t210 + t264;
t33 = qJD(1) * t70 + t201 * t224;
t32 = t203 * t224 + t270 * t298 + t144;
t29 = t201 * t307 + t203 * t308 + t305;
t24 = (-t291 + (t187 * t344 - t320) * qJD(3)) * t201 + (t203 * t208 - t193) * qJD(1) - t264 + t280;
t23 = t186 + t289 * t275 + (t206 - t313) * qJD(1) + t281 + t339;
t22 = (-t119 * t294 - t40) * t187 + (-qJD(3) * t102 + t119 * t297 + t201 * t82) * t188;
t21 = (t119 * t292 + t41) * t187 + (qJD(3) * t101 + t119 * t298 - t203 * t82) * t188;
t20 = (qJD(1) * t137 + t266) * t203 + (t258 * t294 + (-t136 - t141 + t359) * qJD(1)) * t201 + t287;
t14 = t229 * t296 + (-t201 * t41 + t203 * t40 + (-t101 * t203 - t102 * t201) * qJD(1)) * t188;
t13 = t187 * t45 + t188 * t250;
t12 = t187 * t44 + t188 * t251;
t11 = (-pkin(7) * t275 + qJD(1) * t307 + t186 + t41) * t203 + (-pkin(7) * t276 + t40 + (-t141 - t308 + t193) * qJD(1)) * t201 + t287;
t8 = -t96 * t275 + t81 * t100 + t145 * t36 + t146 * t38 + t80 * t98 + (t203 * t34 - t298 * t96) * t188;
t7 = -t95 * t275 + t145 * t37 + t146 * t39 + t80 * t97 + t81 * t99 + (t203 * t35 - t298 * t95) * t188;
t6 = -t96 * t276 + t79 * t100 + t147 * t36 + t148 * t38 + t78 * t98 + (t201 * t34 + t297 * t96) * t188;
t5 = -t95 * t276 + t147 * t37 + t148 * t39 + t78 * t97 + t79 * t99 + (t201 * t35 + t297 * t95) * t188;
t4 = qJD(1) * t251 + t201 * t7 - t203 * t8;
t3 = qJD(1) * t250 + t201 * t5 - t203 * t6;
t2 = (-qJD(3) * t251 + t16) * t187 + (-qJD(1) * t17 + qJD(3) * t44 + t201 * t8 + t203 * t7) * t188;
t1 = (-qJD(3) * t250 + t15) * t187 + (-qJD(1) * t18 + qJD(3) * t45 + t201 * t6 + t203 * t5) * t188;
t19 = [-t202 * t118 * t290 + t263 + (t23 * t51 + t24 * t50) * t351 + (t42 * t68 + t43 * t67) * t352 + (t106 * t65 + t107 * t64) * t353 + 0.2e1 * m(3) * (t104 * t115 + t105 * t114) + (-t202 * t76 - t336) * t188 + (t233 + t242 + t373) * t296 + (t231 + t239 - t372) * t295; m(6) * (qJD(1) * t247 + t201 * t24 - t203 * t23) + m(5) * (qJD(1) * t245 + t201 * t43 - t203 * t42) + m(4) * (qJD(1) * t357 + t201 * t65 - t203 * t64) + m(3) * (-t104 * t203 + t201 * t105 + (t114 * t203 + t115 * t201) * qJD(1)); 0; (t272 * t203 - t284) * t203 + (t272 * t201 + t288) * t201 + m(4) * ((-t201 * t64 - t203 * t65) * t165 - t357 * t156) + m(5) * (t111 * t42 + t112 * t43 + t62 * t67 + t63 * t68) + m(6) * (t23 * t69 + t24 * t70 + t32 * t50 + t33 * t51) + ((-t368 * qJD(3) + t203 * t370) * t347 + (t369 * qJD(3) + t201 * t370) * t346 + (t346 * t366 - t347 * t367) * qJD(1)) * t187 + ((t366 * qJD(3) + t203 * t371) * t347 + (t367 * qJD(3) + t201 * t371) * t346 + (t346 * t368 + t347 * t369) * qJD(1)) * t188 + ((-t107 * t343 + (-t126 / 0.2e1 + t123 / 0.2e1) * t188 + (-t128 / 0.2e1 + t125 / 0.2e1) * t187 - t283) * t203 + (t106 * t343 + (t356 / 0.2e1 + t122 / 0.2e1) * t188 + (t355 / 0.2e1 + t124 / 0.2e1) * t187 + t282) * t201) * qJD(1); m(5) * (t62 * t201 - t203 * t63 + (t111 * t201 + t112 * t203) * qJD(1)) + m(6) * (t32 * t201 - t203 * t33 + (t201 * t69 + t203 * t70) * qJD(1)); t201 * t4 + (t11 * t29 + t32 * t70 + t33 * t69) * t351 - t203 * t3 + (t111 * t63 + t112 * t62 + t20 * t49) * t352 + ((t134 * t201 + t135 * t203) * ((qJD(1) * t134 + t204) * t203 + (-t201 * t217 + (-t135 + t358) * qJD(1)) * t201) + t301 * t165 * t156) * t353 + t201 * ((t201 * t83 + (t54 + t362) * qJD(1)) * t201 + (t55 * qJD(1) + (t122 * t295 + t124 * t296) * t203 + (-t84 + (-t123 * t188 - t125 * t187) * qJD(3) + (t121 - t228) * qJD(1)) * t201) * t203) + t201 * ((t201 * t94 + (t57 - t363) * qJD(1)) * t201 + (t56 * qJD(1) + (t295 * t356 + t296 * t355) * t203 + (-t93 + (t126 * t188 + t128 * t187) * qJD(3) + (t130 + t225) * qJD(1)) * t201) * t203) - t203 * ((t203 * t93 + (t58 - t361) * qJD(1)) * t203 + (t59 * qJD(1) + (t126 * t295 + t128 * t296 + t299) * t201 + (-t94 + (t187 * t355 + t188 * t356) * qJD(3) + t226 * qJD(1)) * t203) * t201) - t203 * ((t203 * t84 + (t53 + t360) * qJD(1)) * t203 + (t52 * qJD(1) + (-t123 * t295 - t125 * t296 + t300) * t201 + (-t83 + (t122 * t188 + t124 * t187) * qJD(3) - t227 * qJD(1)) * t203) * t201) + (t18 + (-t52 - t59) * t203 + (t53 + t58) * t201) * t298 + (t17 + (-t54 - t57) * t203 + (t55 + t56) * t201) * t297; 0.2e1 * (t245 * t350 + t247 * t349) * t295 + 0.2e1 * ((t201 * t23 + t203 * t24 + t297 * t51 - t298 * t50) * t349 + (t201 * t42 + t203 * t43 + t297 * t68 - t298 * t67) * t350) * t187; 0; 0.2e1 * ((t292 * t70 + t294 * t69 - t11) * t349 + (t111 * t294 + t112 * t292 - t20) * t350) * t188 + 0.2e1 * ((qJD(3) * t29 + t201 * t33 + t203 * t32 + t297 * t69 - t298 * t70) * t349 + (qJD(3) * t49 + t111 * t297 - t112 * t298 + t201 * t63 + t203 * t62) * t350) * t187; 0.4e1 * (t350 + t349) * (-0.1e1 + t301) * t187 * t295; m(6) * (t21 * t51 + t22 * t50 + t23 * t61 + t24 * t60) + (-t201 * t282 + t203 * t283) * t296 + (t288 * t203 + t284 * t201 + (t201 * t283 + t203 * t282) * qJD(1)) * t188 + t340; m(6) * (qJD(1) * t246 + t201 * t22 - t203 * t21); m(6) * (-t11 * t46 + t14 * t29 + t21 * t69 + t22 * t70 + t32 * t60 + t33 * t61) + ((qJD(1) * t30 - t10) * t348 + t12 * t341 - t1 / 0.2e1 + t17 * t271) * t203 + ((qJD(1) * t31 + t9) * t348 + t2 / 0.2e1 + t13 * t341 + t18 * t271) * t201 + (t248 * t364 + t4 * t345 + t3 * t347 + (-t201 * t17 / 0.2e1 + t18 * t345) * qJD(1)) * t188; m(6) * ((qJD(3) * t246 - t14) * t188 + (-qJD(3) * t46 + t201 * t21 + t203 * t22 + (-t201 * t60 + t203 * t61) * qJD(1)) * t187); (-t14 * t46 + t21 * t61 + t22 * t60) * t351 + ((-t12 * t203 - t13 * t201 - t187 * t249) * qJD(3) + t340) * t187 + (t203 * t2 + t201 * t1 + t187 * (t10 * t201 + t203 * t9) + (t187 * t48 + t188 * t249) * qJD(3) + (-t12 * t201 + t13 * t203 - t187 * t248) * qJD(1)) * t188;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
