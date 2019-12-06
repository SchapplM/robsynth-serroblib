% Calculate time derivative of joint inertia matrix for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:17
% EndTime: 2019-12-05 18:11:37
% DurationCPUTime: 6.80s
% Computational Cost: add. (13713->499), mult. (10466->705), div. (0->0), fcn. (8208->10), ass. (0->288)
t212 = sin(qJ(1));
t213 = cos(qJ(1));
t205 = pkin(9) + qJ(3);
t196 = qJ(4) + t205;
t189 = qJ(5) + t196;
t183 = sin(t189);
t349 = rSges(6,2) * t183;
t184 = cos(t189);
t353 = rSges(6,1) * t184;
t264 = -t349 + t353;
t128 = t212 * rSges(6,3) + t213 * t264;
t202 = t212 * rSges(5,3);
t186 = sin(t196);
t187 = cos(t196);
t354 = rSges(5,1) * t187;
t265 = -rSges(5,2) * t186 + t354;
t136 = t213 * t265 + t202;
t348 = rSges(6,2) * t184;
t157 = rSges(6,1) * t183 + t348;
t206 = qJD(3) + qJD(4);
t195 = qJD(5) + t206;
t232 = t157 * t195;
t193 = sin(t205);
t194 = cos(t205);
t343 = Icges(4,4) * t194;
t254 = -Icges(4,2) * t193 + t343;
t140 = Icges(4,6) * t212 + t213 * t254;
t344 = Icges(4,4) * t193;
t258 = Icges(4,1) * t194 - t344;
t142 = Icges(4,5) * t212 + t213 * t258;
t241 = t140 * t193 - t142 * t194;
t376 = t212 * t241;
t341 = Icges(5,4) * t187;
t252 = -Icges(5,2) * t186 + t341;
t132 = Icges(5,6) * t212 + t213 * t252;
t342 = Icges(5,4) * t186;
t256 = Icges(5,1) * t187 - t342;
t134 = Icges(5,5) * t212 + t213 * t256;
t243 = t132 * t186 - t134 * t187;
t375 = t212 * t243;
t339 = Icges(6,4) * t184;
t251 = -Icges(6,2) * t183 + t339;
t118 = Icges(6,6) * t212 + t213 * t251;
t340 = Icges(6,4) * t183;
t255 = Icges(6,1) * t184 - t340;
t120 = Icges(6,5) * t212 + t213 * t255;
t245 = t118 * t183 - t120 * t184;
t374 = t212 * t245;
t139 = -Icges(4,6) * t213 + t212 * t254;
t141 = -Icges(4,5) * t213 + t212 * t258;
t242 = t139 * t193 - t141 * t194;
t373 = t213 * t242;
t131 = -Icges(5,6) * t213 + t212 * t252;
t133 = -Icges(5,5) * t213 + t212 * t256;
t244 = t131 * t186 - t133 * t187;
t372 = t213 * t244;
t117 = -Icges(6,6) * t213 + t212 * t251;
t119 = -Icges(6,5) * t213 + t212 * t255;
t246 = t117 * t183 - t119 * t184;
t371 = t213 * t246;
t211 = -pkin(6) - qJ(2);
t210 = cos(pkin(9));
t188 = t210 * pkin(2) + pkin(1);
t352 = rSges(4,2) * t193;
t355 = rSges(4,1) * t194;
t266 = -t352 + t355;
t235 = -t188 - t266;
t101 = (rSges(4,3) - t211) * t213 + t235 * t212;
t203 = t212 * rSges(4,3);
t303 = t213 * t355 + t203;
t102 = -t212 * t211 + (t188 - t352) * t213 + t303;
t370 = t101 * t213 + t102 * t212;
t248 = Icges(6,5) * t184 - Icges(6,6) * t183;
t115 = -Icges(6,3) * t213 + t212 * t248;
t369 = qJD(1) * t115;
t249 = Icges(5,5) * t187 - Icges(5,6) * t186;
t129 = -Icges(5,3) * t213 + t212 * t249;
t368 = qJD(1) * t129;
t250 = Icges(4,5) * t194 - Icges(4,6) * t193;
t137 = -Icges(4,3) * t213 + t212 * t250;
t155 = Icges(6,2) * t184 + t340;
t156 = Icges(6,1) * t183 + t339;
t240 = t155 * t183 - t156 * t184;
t367 = qJD(1) * t240 + t248 * t195;
t163 = Icges(5,2) * t187 + t342;
t164 = Icges(5,1) * t186 + t341;
t239 = t163 * t186 - t164 * t187;
t366 = qJD(1) * t239 + t249 * t206;
t237 = rSges(3,1) * t210 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t347 = rSges(3,3) + qJ(2);
t126 = t212 * t347 + t213 * t237;
t204 = -pkin(7) + t211;
t301 = t204 - t211;
t172 = pkin(3) * t194 + t188;
t306 = t172 - t188;
t103 = t212 * t306 + t213 * t301;
t199 = -pkin(8) + t204;
t302 = t199 - t204;
t152 = pkin(4) * t187 + t172;
t308 = t152 - t172;
t84 = t212 * t308 + t213 * t302;
t365 = 2 * m(4);
t364 = 2 * m(5);
t363 = 2 * m(6);
t207 = t212 ^ 2;
t208 = t213 ^ 2;
t362 = t212 / 0.2e1;
t361 = -t213 / 0.2e1;
t170 = rSges(4,1) * t193 + rSges(4,2) * t194;
t360 = m(4) * t170;
t350 = rSges(5,2) * t187;
t165 = rSges(5,1) * t186 + t350;
t359 = m(5) * t165;
t358 = pkin(3) * t193;
t27 = -t115 * t213 - t212 * t246;
t322 = t155 * t195;
t64 = qJD(1) * t118 - t212 * t322;
t276 = t119 * t195 + t64;
t321 = t156 * t195;
t66 = qJD(1) * t120 - t212 * t321;
t278 = t117 * t195 - t66;
t116 = Icges(6,3) * t212 + t213 * t248;
t28 = -t116 * t213 - t374;
t299 = qJD(1) * t116;
t317 = t184 * t195;
t318 = t183 * t195;
t154 = Icges(6,5) * t183 + Icges(6,6) * t184;
t226 = t195 * t154;
t61 = -t213 * t226 - t369;
t62 = -t212 * t226 + t299;
t63 = -qJD(1) * t117 - t213 * t322;
t65 = -qJD(1) * t119 - t213 * t321;
t2 = (t213 * t62 + (t28 + t371) * qJD(1)) * t213 + (t27 * qJD(1) + (-t118 * t317 - t120 * t318 - t183 * t63 + t184 * t65 + t299) * t212 + (-t61 + t278 * t184 + t276 * t183 + (-t115 - t245) * qJD(1)) * t213) * t212;
t357 = t213 * t2;
t166 = t213 * t172;
t104 = -t188 * t213 - t212 * t301 + t166;
t356 = t212 * t103 + t213 * t104;
t346 = rSges(6,3) - t199;
t149 = t213 * t152;
t85 = -t212 * t302 + t149 - t166;
t345 = -t128 - t85;
t280 = -t165 - t358;
t112 = t280 * t213;
t327 = t112 * t213;
t326 = t139 * t194;
t325 = t140 * t194;
t324 = t141 * t193;
t323 = t142 * t193;
t320 = t163 * t206;
t319 = t164 * t206;
t316 = t186 * t206;
t315 = t187 * t206;
t314 = t195 * t213;
t313 = t206 * t212;
t312 = t206 * t213;
t150 = t265 * t206;
t311 = t212 * t150;
t310 = t212 * t199;
t127 = -rSges(6,3) * t213 + t212 * t264;
t67 = t212 * t127 + t213 * t128;
t135 = -rSges(5,3) * t213 + t212 * t265;
t73 = t212 * t135 + t213 * t136;
t296 = qJD(1) * t212;
t284 = t186 * t296;
t309 = pkin(4) * t284 + t157 * t296;
t295 = qJD(1) * t213;
t307 = rSges(6,3) * t295 + t296 * t349;
t305 = rSges(5,2) * t284 + rSges(5,3) * t295;
t294 = qJD(3) * t193;
t286 = pkin(3) * t294;
t304 = t204 * t296 + t212 * t286;
t300 = t207 + t208;
t130 = Icges(5,3) * t212 + t213 * t249;
t298 = qJD(1) * t130;
t138 = Icges(4,3) * t212 + t213 * t250;
t297 = qJD(1) * t138;
t293 = qJD(3) * t194;
t292 = qJD(3) * t212;
t275 = t120 * t195 + t63;
t277 = -t118 * t195 + t65;
t29 = t212 * t115 - t371;
t30 = t212 * t116 - t213 * t245;
t291 = t212 * ((t212 * t61 + (t29 + t374) * qJD(1)) * t212 + (t30 * qJD(1) + (t117 * t317 + t119 * t318 + t183 * t64 - t184 * t66 - t369) * t213 + (-t62 + t277 * t184 - t275 * t183 + (t116 - t246) * qJD(1)) * t212) * t213) + (t28 * t212 - t213 * t27) * t296 + (t30 * t212 - t213 * t29) * t295;
t182 = t211 * t296;
t270 = t213 * t286;
t290 = t212 * (t295 * t306 + t182 - t304) + t213 * (-qJD(1) * t103 - t270) + t103 * t295;
t289 = t127 * t295 + t212 * (t128 * qJD(1) - t212 * t232) + t213 * (-t314 * t348 + (-t183 * t314 - t184 * t296) * rSges(6,1) + t307);
t231 = t165 * t206;
t288 = t135 * t295 + t212 * (t136 * qJD(1) - t212 * t231) + t213 * (-t312 * t350 + (-t186 * t312 - t187 * t296) * rSges(5,1) + t305);
t287 = t213 * t352;
t285 = pkin(3) * t293;
t283 = t193 * t296;
t282 = t296 / 0.2e1;
t281 = t295 / 0.2e1;
t279 = -pkin(4) * t186 - t157;
t79 = qJD(1) * t134 - t212 * t319;
t274 = t131 * t206 - t79;
t78 = -qJD(1) * t133 - t213 * t319;
t273 = -t132 * t206 + t78;
t77 = qJD(1) * t132 - t212 * t320;
t272 = t133 * t206 + t77;
t76 = -qJD(1) * t131 - t213 * t320;
t271 = t134 * t206 + t76;
t23 = t212 * t84 + t213 * t85 + t67;
t124 = t264 * t195;
t269 = -pkin(4) * t315 - t124;
t268 = t291 - t357;
t33 = -t129 * t213 - t212 * t244;
t34 = -t130 * t213 - t375;
t35 = t212 * t129 - t372;
t36 = t212 * t130 - t213 * t243;
t162 = Icges(5,5) * t186 + Icges(5,6) * t187;
t223 = t206 * t162;
t74 = -t213 * t223 - t368;
t75 = -t212 * t223 + t298;
t267 = (t34 * t212 - t213 * t33) * t296 + (t36 * t212 - t213 * t35) * t295 + t212 * ((t212 * t74 + (t35 + t375) * qJD(1)) * t212 + (t36 * qJD(1) + (t131 * t315 + t133 * t316 + t186 * t77 - t187 * t79 - t368) * t213 + (-t75 + t273 * t187 - t271 * t186 + (t130 - t244) * qJD(1)) * t212) * t213) + t291;
t233 = -t152 - t264;
t70 = t212 * t233 + t213 * t346;
t71 = t128 + t149 - t310;
t261 = t212 * t71 + t213 * t70;
t234 = -t172 - t265;
t92 = (rSges(5,3) - t204) * t213 + t234 * t212;
t93 = -t212 * t204 + t136 + t166;
t260 = t212 * t93 + t213 * t92;
t236 = t279 - t358;
t97 = t236 * t212;
t98 = t236 * t213;
t259 = t212 * t97 + t213 * t98;
t257 = Icges(4,1) * t193 + t343;
t253 = Icges(4,2) * t194 + t344;
t105 = t279 * t212;
t106 = t279 * t213;
t247 = t105 * t212 + t106 * t213;
t153 = -pkin(4) * t316 - t286;
t151 = t213 * t153;
t238 = t212 * (t212 * t153 + (t213 * t308 - t310) * qJD(1) + t304) + t213 * (-qJD(1) * t84 + t151 + t270) + t84 * t295 + t289;
t230 = qJD(3) * t170;
t122 = t251 * t195;
t123 = t255 * t195;
t215 = qJD(1) * t154 + (t123 - t322) * t184 + (-t122 - t321) * t183;
t229 = (t183 * t277 + t184 * t275 + t212 * t367 + t215 * t213) * t362 + (-t183 * t278 + t184 * t276 + t215 * t212 - t367 * t213) * t361 + (t117 * t184 + t119 * t183 - t154 * t213 - t212 * t240) * t282 + (t118 * t184 + t120 * t183 + t212 * t154 - t213 * t240) * t281;
t222 = qJD(3) * t257;
t221 = qJD(3) * t253;
t220 = qJD(3) * (-Icges(4,5) * t193 - Icges(4,6) * t194);
t219 = t269 - t285;
t4 = (t213 * t75 + (t34 + t372) * qJD(1)) * t213 + (t33 * qJD(1) + (-t132 * t315 - t134 * t316 - t186 * t76 + t187 * t78 + t298) * t212 + (-t74 + t274 * t187 + t272 * t186 + (-t129 - t243) * qJD(1)) * t213) * t212;
t218 = (-t4 - t2) * t213 + t267;
t217 = rSges(4,2) * t283 + rSges(4,3) * t295 - t213 * t230;
t125 = -t212 * t237 + t213 * t347;
t144 = t252 * t206;
t145 = t256 * t206;
t214 = qJD(1) * t162 + (t145 - t320) * t187 + (-t144 - t319) * t186;
t216 = t229 + (t186 * t273 + t187 * t271 + t212 * t366 + t214 * t213) * t362 + (-t186 * t274 + t187 * t272 + t214 * t212 - t366 * t213) * t361 + (t131 * t187 + t133 * t186 - t162 * t213 - t212 * t239) * t282 + (t132 * t187 + t134 * t186 + t212 * t162 - t213 * t239) * t281;
t198 = qJD(2) * t213;
t197 = qJD(2) * t212;
t178 = pkin(3) * t283;
t161 = t266 * qJD(3);
t147 = -t287 + t303;
t146 = -rSges(4,3) * t213 + t212 * t266;
t111 = t280 * t212;
t100 = -qJD(1) * t126 + t198;
t99 = qJD(1) * t125 + t197;
t87 = t212 * t220 + t297;
t86 = -qJD(1) * t137 + t213 * t220;
t58 = -t165 * t295 - t311 + (-t193 * t295 - t194 * t292) * pkin(3);
t57 = t165 * t296 + t178 + (-t150 - t285) * t213;
t56 = t182 + t198 + t170 * t292 + (t213 * t235 - t203) * qJD(1);
t55 = t197 + (-t211 * t213 + (-t188 - t355) * t212) * qJD(1) + t217;
t48 = -t157 * t295 - t212 * t124 + (-t186 * t295 - t187 * t313) * pkin(4);
t47 = t213 * t269 + t309;
t44 = t212 * t138 - t213 * t241;
t43 = t212 * t137 - t373;
t42 = -t138 * t213 - t376;
t41 = -t137 * t213 - t212 * t242;
t38 = t198 + t165 * t313 + (t213 * t234 - t202) * qJD(1) + t304;
t37 = t197 + (-t172 - t354) * t296 + (-qJD(1) * t204 - t231 - t286) * t213 + t305;
t32 = qJD(1) * t98 + t212 * t219;
t31 = t213 * t219 + t178 + t309;
t26 = t198 + (-t153 + t232) * t212 + (-t212 * t346 + t213 * t233) * qJD(1);
t25 = t151 + t197 - t213 * t232 + (-t199 * t213 + (-t152 - t353) * t212) * qJD(1) + t307;
t24 = t73 + t356;
t22 = -t136 * t296 + t288;
t21 = -t128 * t296 + t289;
t10 = t23 + t356;
t7 = (-t104 - t136) * t296 + t288 + t290;
t6 = t296 * t345 + t238;
t5 = (-t104 + t345) * t296 + t238 + t290;
t1 = [(t25 * t71 + t26 * t70) * t363 + t156 * t317 + t183 * t123 - t155 * t318 + t184 * t122 + (t37 * t93 + t38 * t92) * t364 + t164 * t315 + t186 * t145 - t163 * t316 + t187 * t144 + (t101 * t56 + t102 * t55) * t365 + 0.2e1 * m(3) * (t100 * t125 + t126 * t99) + (t258 - t253) * t294 + (t257 + t254) * t293; m(6) * (qJD(1) * t261 + t212 * t26 - t213 * t25) + m(5) * (qJD(1) * t260 + t212 * t38 - t213 * t37) + m(4) * (qJD(1) * t370 + t212 * t56 - t213 * t55) + m(3) * (t212 * t100 - t213 * t99 + (t125 * t213 + t126 * t212) * qJD(1)); 0; t216 + ((t325 / 0.2e1 + t323 / 0.2e1 - t102 * t360) * t213 + (t326 / 0.2e1 + t324 / 0.2e1 + t101 * t360) * t212) * qJD(1) + m(6) * (t25 * t97 + t26 * t98 + t31 * t70 + t32 * t71) + m(5) * (t111 * t37 + t112 * t38 + t57 * t92 + t58 * t93) + m(4) * ((-t212 * t55 - t213 * t56) * t170 - t370 * t161) + (-qJD(3) * t241 + t193 * (-qJD(1) * t141 - t213 * t222) + t194 * (-qJD(1) * t139 - t213 * t221)) * t362 + (-qJD(3) * t242 + t193 * (qJD(1) * t142 - t212 * t222) + t194 * (qJD(1) * t140 - t212 * t221)) * t361 + (t208 / 0.2e1 + t207 / 0.2e1) * t250 * qJD(3); m(5) * (t57 * t212 - t213 * t58 + (t111 * t212 + t327) * qJD(1)) + m(6) * (qJD(1) * t259 + t31 * t212 - t213 * t32); (t10 * t5 + t31 * t98 + t32 * t97) * t363 - t357 + (t111 * t58 + t112 * t57 + t24 * t7) * t364 - t213 * t4 + (t42 * t212 - t213 * t41) * t296 - t213 * ((t213 * t87 + (t42 + t373) * qJD(1)) * t213 + (t41 * qJD(1) + (-t140 * t293 - t142 * t294 + t297) * t212 + (-t86 + (t324 + t326) * qJD(3) - t241 * qJD(1)) * t213) * t212) + (t44 * t212 - t213 * t43) * t295 + t212 * ((t212 * t86 + (t43 + t376) * qJD(1)) * t212 + (t44 * qJD(1) + (t139 * t293 + t141 * t294) * t213 + (-t87 + (-t323 - t325) * qJD(3) + (t138 - t242) * qJD(1)) * t212) * t213) + ((t212 * t146 + t147 * t213) * ((qJD(1) * t146 + t217) * t213 + (-t212 * t230 + (-t147 - t287 + t203) * qJD(1)) * t212) + t300 * t170 * t161) * t365 + t267; t216 + m(6) * (t105 * t25 + t106 * t26 + t47 * t70 + t48 * t71) + (-t212 * t37 - t213 * t38 + (t212 * t92 - t213 * t93) * qJD(1)) * t359 - m(5) * t260 * t150; m(6) * (qJD(1) * t247 + t47 * t212 - t213 * t48); m(6) * (t6 * t10 + t105 * t32 + t106 * t31 + t23 * t5 + t47 * t98 + t48 * t97) + m(5) * (-t111 * t311 - t150 * t327 + t22 * t24 + t73 * t7) + (-t212 * t58 - t213 * t57 + (-t111 * t213 + t112 * t212) * qJD(1)) * t359 + t218; (t150 * t165 * t300 + t73 * t22) * t364 + (t105 * t48 + t106 * t47 + t23 * t6) * t363 + t218; m(6) * (-t261 * t124 + (-t212 * t25 - t213 * t26 + (t212 * t70 - t213 * t71) * qJD(1)) * t157) + t229; 0; m(6) * (t21 * t10 + t67 * t5 - t259 * t124 + (-t212 * t32 - t213 * t31 + (t212 * t98 - t213 * t97) * qJD(1)) * t157) + t268; m(6) * (t21 * t23 + t67 * t6 - t247 * t124 + (-t212 * t48 - t213 * t47 + (-t105 * t213 + t106 * t212) * qJD(1)) * t157) + t268; (t124 * t157 * t300 + t67 * t21) * t363 + t268;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
