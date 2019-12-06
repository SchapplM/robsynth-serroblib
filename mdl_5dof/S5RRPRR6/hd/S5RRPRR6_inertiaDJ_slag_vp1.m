% Calculate time derivative of joint inertia matrix for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:42
% EndTime: 2019-12-05 18:35:53
% DurationCPUTime: 5.71s
% Computational Cost: add. (16061->523), mult. (15018->737), div. (0->0), fcn. (14002->10), ass. (0->262)
t215 = qJD(1) + qJD(2);
t218 = sin(pkin(9));
t224 = -pkin(8) - pkin(7);
t287 = t218 * t224;
t262 = t215 * t287;
t222 = cos(qJ(4));
t271 = qJD(4) * t222;
t317 = pkin(4) * t271 + t262;
t217 = qJ(1) + qJ(2);
t212 = cos(t217);
t219 = cos(pkin(9));
t207 = pkin(4) * t222 + pkin(3);
t308 = pkin(3) - t207;
t238 = pkin(7) * t218 + t219 * t308;
t210 = sin(t217);
t220 = sin(qJ(4));
t295 = t210 * t220;
t316 = -pkin(4) * t295 + t212 * t238;
t315 = 2 * m(3);
t314 = 2 * m(4);
t313 = 2 * m(5);
t312 = 2 * m(6);
t213 = t218 ^ 2;
t221 = sin(qJ(1));
t311 = pkin(1) * t221;
t223 = cos(qJ(1));
t310 = pkin(1) * t223;
t272 = qJD(4) * t220;
t268 = pkin(4) * t272;
t248 = t219 * t268;
t261 = t210 * t248 + t317 * t212;
t289 = t215 * t218;
t264 = t212 * t289;
t216 = qJ(4) + qJ(5);
t211 = cos(t216);
t214 = qJD(4) + qJD(5);
t251 = t214 * t219 - t215;
t240 = t251 * t210;
t209 = sin(t216);
t288 = t215 * t219;
t252 = -t214 + t288;
t243 = t209 * t252;
t117 = t211 * t240 + t212 * t243;
t242 = t252 * t211;
t118 = t209 * t240 - t212 * t242;
t283 = t118 * rSges(6,1) + t117 * rSges(6,2);
t69 = -rSges(6,3) * t264 + t283;
t307 = -t316 * t215 - t261 - t69;
t306 = pkin(1) * qJD(1);
t290 = t214 * t218;
t300 = Icges(6,4) * t209;
t146 = (-Icges(6,2) * t211 - t300) * t290;
t157 = -Icges(6,5) * t219 + (Icges(6,1) * t211 - t300) * t218;
t299 = Icges(6,4) * t211;
t156 = -Icges(6,6) * t219 + (-Icges(6,2) * t209 + t299) * t218;
t263 = t214 * t211 * t156;
t231 = -t263 + (-t157 * t214 - t146) * t209;
t147 = (-Icges(6,1) * t209 - t299) * t290;
t140 = t218 * t211 * t147;
t145 = (-Icges(6,5) * t209 - Icges(6,6) * t211) * t290;
t254 = -t219 * t145 + t140;
t53 = t218 * t231 + t254;
t305 = t53 * t219;
t304 = -rSges(4,3) - qJ(3);
t148 = (-rSges(6,1) * t209 - rSges(6,2) * t211) * t290;
t293 = t212 * t218;
t241 = t251 * t212;
t115 = t210 * t243 - t211 * t241;
t116 = -t209 * t241 - t210 * t242;
t266 = t210 * t289;
t68 = t116 * rSges(6,1) + t115 * rSges(6,2) - rSges(6,3) * t266;
t303 = t148 * t293 + t219 * t68;
t302 = Icges(5,4) * t220;
t301 = Icges(5,4) * t222;
t298 = t210 * t215;
t297 = t210 * t218;
t296 = t210 * t219;
t294 = t212 * t215;
t292 = t212 * t219;
t291 = t212 * t220;
t286 = t219 * t220;
t285 = t219 * t222;
t273 = qJD(4) * t218;
t176 = (-Icges(5,2) * t222 - t302) * t273;
t284 = t220 * t176;
t161 = -t209 * t292 + t210 * t211;
t162 = t209 * t210 + t211 * t292;
t244 = -rSges(6,1) * t162 - rSges(6,2) * t161;
t114 = rSges(6,3) * t293 - t244;
t158 = -rSges(6,3) * t219 + (rSges(6,1) * t211 - rSges(6,2) * t209) * t218;
t144 = t158 * t293;
t90 = t219 * t114 + t144;
t159 = t209 * t296 + t211 * t212;
t160 = t209 * t212 - t211 * t296;
t278 = t160 * rSges(6,1) + t159 * rSges(6,2);
t113 = -rSges(6,3) * t297 + t278;
t274 = pkin(4) * t291 + t210 * t287;
t125 = t210 * t238 + t274;
t282 = -t113 - t125;
t194 = t212 * t287;
t126 = -t194 - t316;
t281 = -t114 - t126;
t236 = t210 * t285 - t291;
t237 = -t210 * t222 + t212 * t286;
t133 = qJD(4) * t236 + t215 * t237;
t171 = t210 * t286 + t212 * t222;
t174 = t212 * t285 + t295;
t134 = qJD(4) * t171 - t174 * t215;
t280 = t134 * rSges(5,1) + t133 * rSges(5,2);
t279 = t215 * t144 + t148 * t297;
t277 = -rSges(5,1) * t236 + t171 * rSges(5,2);
t265 = t210 * t288;
t276 = -t207 * t265 - t212 * t248;
t275 = pkin(3) * t265 + pkin(7) * t266;
t168 = rSges(3,1) * t298 + rSges(3,2) * t294;
t269 = t223 * t306;
t200 = rSges(4,2) * t297;
t143 = t158 * t297;
t110 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t293;
t112 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t293;
t61 = Icges(6,5) * t116 + Icges(6,6) * t115 - Icges(6,3) * t266;
t63 = Icges(6,4) * t116 + Icges(6,2) * t115 - Icges(6,6) * t266;
t65 = Icges(6,1) * t116 + Icges(6,4) * t115 - Icges(6,5) * t266;
t10 = -t219 * t61 + ((-t110 * t214 + t65) * t211 + (-t112 * t214 - t63) * t209) * t218;
t107 = Icges(6,5) * t160 + Icges(6,6) * t159 - Icges(6,3) * t297;
t108 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t293;
t109 = Icges(6,4) * t160 + Icges(6,2) * t159 - Icges(6,6) * t297;
t111 = Icges(6,1) * t160 + Icges(6,4) * t159 - Icges(6,5) * t297;
t155 = -Icges(6,3) * t219 + (Icges(6,5) * t211 - Icges(6,6) * t209) * t218;
t20 = t115 * t156 + t116 * t157 + t146 * t161 + t147 * t162 + (t145 * t212 - t155 * t298) * t218;
t29 = t107 * t293 + t109 * t161 + t111 * t162;
t30 = t108 * t293 + t110 * t161 + t112 * t162;
t43 = -t107 * t219 + (-t109 * t209 + t111 * t211) * t218;
t44 = -t108 * t219 + (-t110 * t209 + t112 * t211) * t218;
t62 = Icges(6,5) * t118 + Icges(6,6) * t117 - Icges(6,3) * t264;
t64 = Icges(6,4) * t118 + Icges(6,2) * t117 - Icges(6,6) * t264;
t66 = Icges(6,1) * t118 + Icges(6,4) * t117 - Icges(6,5) * t264;
t9 = -t219 * t62 + ((-t109 * t214 + t66) * t211 + (-t111 * t214 - t64) * t209) * t218;
t260 = -t219 * (-t305 + ((-t215 * t43 + t10) * t212 + (-t215 * t44 - t9) * t210) * t218) + (-t20 * t219 + (-(t115 * t109 + t116 * t111 + t161 * t64 + t162 * t66) * t210 - t29 * t294 + (t115 * t110 + t116 * t112 + t161 * t63 + t162 * t65) * t212 - t30 * t298 + (-(-t107 * t298 + t212 * t62) * t210 + (-t108 * t298 + t212 * t61) * t212) * t218) * t218) * t293;
t259 = -t297 / 0.2e1;
t258 = t293 / 0.2e1;
t257 = -t289 / 0.2e1;
t256 = -rSges(4,1) * t219 - pkin(2);
t255 = -pkin(4) * t220 - qJ(3);
t184 = -rSges(3,1) * t212 + t210 * rSges(3,2);
t177 = (-Icges(5,1) * t220 - t301) * t273;
t154 = t218 * t222 * t177;
t175 = (-Icges(5,5) * t220 - Icges(5,6) * t222) * t273;
t253 = -t219 * t175 + t154;
t199 = pkin(2) * t298;
t250 = -qJD(3) * t210 + t199;
t249 = t213 * t268;
t247 = t210 * t257;
t246 = t212 * t257;
t169 = -rSges(3,1) * t294 + rSges(3,2) * t298;
t183 = -rSges(3,1) * t210 - rSges(3,2) * t212;
t245 = -rSges(5,1) * t174 + rSges(5,2) * t237;
t239 = -rSges(6,3) * t218 - t207 * t219 - pkin(2);
t131 = -qJD(4) * t174 + t171 * t215;
t132 = -qJD(4) * t237 - t215 * t236;
t81 = rSges(5,1) * t132 + rSges(5,2) * t131 - rSges(5,3) * t266;
t21 = t117 * t156 + t118 * t157 + t146 * t159 + t147 * t160 + (-t145 * t210 - t155 * t294) * t218;
t70 = -t155 * t297 + t156 * t159 + t157 * t160;
t71 = t155 * t293 + t156 * t161 + t157 * t162;
t235 = (t21 + t9) * t259 + (t10 + t20) * t258 + (t44 + t71) * t247 + (t43 + t70) * t246;
t234 = -pkin(3) * t219 - pkin(2) + (-rSges(5,3) - pkin(7)) * t218;
t204 = t212 * qJ(3);
t137 = t212 * rSges(4,3) + t210 * t256 + t200 + t204;
t233 = t235 - t305;
t232 = t210 * t304 + t212 * t256;
t166 = -Icges(5,6) * t219 + (-Icges(5,2) * t220 + t301) * t218;
t167 = -Icges(5,5) * t219 + (Icges(5,1) * t222 - t302) * t218;
t230 = -t284 + (-t166 * t222 - t167 * t220) * qJD(4);
t138 = rSges(4,2) * t293 + t232;
t229 = -qJ(3) * t210 + t212 * t234;
t97 = t210 * t234 + t204 + t277;
t87 = t210 * t239 + t204 + t274 + t278;
t228 = t210 * t255 + t212 * t239;
t203 = qJD(3) * t212;
t104 = rSges(4,2) * t264 + t215 * t232 + t203;
t27 = -t107 * t297 + t109 * t159 + t111 * t160;
t28 = -t108 * t297 + t110 * t159 + t112 * t160;
t2 = -t21 * t219 + (-(t117 * t109 + t118 * t111 + t159 * t64 + t160 * t66) * t210 - t27 * t294 + (t117 * t110 + t118 * t112 + t159 * t63 + t160 * t65) * t212 - t28 * t298 + (-(-t107 * t294 - t210 * t62) * t210 + (-t108 * t294 - t210 * t61) * t212) * t218) * t218;
t5 = -t219 * t70 + (-t210 * t27 + t212 * t28) * t218;
t6 = -t219 * t71 + (-t210 * t29 + t212 * t30) * t218;
t227 = (-t5 * t294 + (-t215 * t6 - t2) * t210) * t218 + t260;
t103 = rSges(4,1) * t265 + (t212 * t304 - t200) * t215 + t250;
t56 = -qJ(3) * t294 + t250 + t275 - t81;
t98 = t229 + t245;
t57 = t215 * t229 + t203 + t280;
t121 = -Icges(5,4) * t236 + Icges(5,2) * t171 - Icges(5,6) * t297;
t123 = -Icges(5,1) * t236 + Icges(5,4) * t171 - Icges(5,5) * t297;
t76 = Icges(5,5) * t134 + Icges(5,6) * t133 - Icges(5,3) * t264;
t78 = Icges(5,4) * t134 + Icges(5,2) * t133 - Icges(5,6) * t264;
t80 = Icges(5,1) * t134 + Icges(5,4) * t133 - Icges(5,5) * t264;
t13 = -t219 * t76 + (-t220 * t78 + t222 * t80 + (-t121 * t222 - t123 * t220) * qJD(4)) * t218;
t122 = Icges(5,4) * t174 - Icges(5,2) * t237 + Icges(5,6) * t293;
t124 = Icges(5,1) * t174 - Icges(5,4) * t237 + Icges(5,5) * t293;
t75 = Icges(5,5) * t132 + Icges(5,6) * t131 - Icges(5,3) * t266;
t77 = Icges(5,4) * t132 + Icges(5,2) * t131 - Icges(5,6) * t266;
t79 = Icges(5,1) * t132 + Icges(5,4) * t131 - Icges(5,5) * t266;
t14 = -t219 * t75 + (-t220 * t77 + t222 * t79 + (-t122 * t222 - t124 * t220) * qJD(4)) * t218;
t165 = -Icges(5,3) * t219 + (Icges(5,5) * t222 - Icges(5,6) * t220) * t218;
t24 = t131 * t166 + t132 * t167 - t237 * t176 + t174 * t177 + (-t165 * t298 + t175 * t212) * t218;
t25 = t133 * t166 + t134 * t167 + t171 * t176 - t236 * t177 + (-t165 * t294 - t175 * t210) * t218;
t119 = -Icges(5,5) * t236 + Icges(5,6) * t171 - Icges(5,3) * t297;
t49 = -t119 * t219 + (-t121 * t220 + t123 * t222) * t218;
t120 = Icges(5,5) * t174 - Icges(5,6) * t237 + Icges(5,3) * t293;
t50 = -t120 * t219 + (-t122 * t220 + t124 * t222) * t218;
t72 = t218 * t230 + t253;
t83 = -t165 * t297 + t166 * t171 - t167 * t236;
t84 = t165 * t293 - t166 * t237 + t167 * t174;
t226 = (-t72 - t53) * t219 + t235 + (t25 + t13) * t259 + (t24 + t14) * t258 + (t50 + t84) * t247 + (t49 + t83) * t246;
t88 = t194 + t228 + t244;
t34 = t215 * t228 + t203 + t261 + t283;
t33 = t199 + t255 * t294 + (-qJD(3) - t317) * t210 - t68 - t276;
t225 = t140 + t154 + (-t175 - t145) * t219 + (t230 + t231) * t218;
t208 = t221 * t306;
t180 = t184 - t310;
t179 = t183 - t311;
t178 = (-rSges(5,1) * t220 - rSges(5,2) * t222) * t273;
t170 = -rSges(5,3) * t219 + (rSges(5,1) * t222 - rSges(5,2) * t220) * t218;
t151 = (pkin(7) + t224) * t219 - t308 * t218;
t150 = t169 - t269;
t149 = t208 + t168;
t136 = t138 - t310;
t135 = t137 - t311;
t128 = rSges(5,3) * t293 - t245;
t127 = -rSges(5,3) * t297 + t277;
t101 = t104 - t269;
t100 = t208 + t103;
t99 = t113 * t266;
t96 = t128 * t219 + t170 * t293;
t95 = -t127 * t219 + t170 * t297;
t94 = t98 - t310;
t93 = t97 - t311;
t91 = t210 * t262 + (t210 * t271 + t215 * t291) * pkin(4) + t275 + t276;
t89 = -t113 * t219 + t143;
t86 = t88 - t310;
t85 = t87 - t311;
t82 = -rSges(5,3) * t264 + t280;
t67 = (-t113 * t212 - t114 * t210) * t218;
t55 = t57 - t269;
t54 = t208 + t56;
t52 = -t219 * t82 + (t170 * t294 + t178 * t210) * t218;
t51 = t219 * t81 + (-t170 * t298 + t178 * t212) * t218;
t48 = t126 * t219 + t151 * t293 + t90;
t47 = t151 * t297 + t219 * t282 + t143;
t42 = t120 * t293 - t122 * t237 + t124 * t174;
t41 = t119 * t293 - t121 * t237 + t123 * t174;
t40 = -t120 * t297 + t122 * t171 - t124 * t236;
t39 = -t119 * t297 + t121 * t171 - t123 * t236;
t38 = -t219 * t69 + t279;
t37 = -t143 * t215 + t303;
t32 = t34 - t269;
t31 = t208 + t33;
t26 = (t210 * t281 + t212 * t282) * t218;
t17 = t151 * t264 - t210 * t249 + t219 * t307 + t279;
t16 = -t212 * t249 + t219 * t91 + (-t151 - t158) * t266 + t303;
t15 = t99 + (-t210 * t68 + (-t114 * t215 - t69) * t212) * t218;
t4 = t99 + ((t125 * t215 - t68 - t91) * t210 + (t215 * t281 + t307) * t212) * t218;
t1 = [(t149 * t180 + t150 * t179) * t315 + (t100 * t136 + t101 * t135) * t314 + (t54 * t94 + t55 * t93) * t313 - t209 * t157 * t290 + (t31 * t86 + t32 * t85) * t312 + t253 + t254 + (-t209 * t146 - t166 * t271 - t167 * t272 - t263 - t284) * t218; m(3) * (t149 * t184 + t150 * t183 + t168 * t180 + t169 * t179) + m(4) * (t100 * t138 + t101 * t137 + t103 * t136 + t104 * t135) + m(5) * (t54 * t98 + t55 * t97 + t56 * t94 + t57 * t93) + m(6) * (t31 * t88 + t32 * t87 + t33 * t86 + t34 * t85) + t225; (t33 * t88 + t34 * t87) * t312 + (t56 * t98 + t57 * t97) * t313 + (t103 * t138 + t104 * t137) * t314 + (t168 * t184 + t169 * t183) * t315 + t225; m(4) * ((t135 * t215 + t100) * t212 + (-t136 * t215 + t101) * t210) + m(5) * ((t215 * t93 + t54) * t212 + (-t215 * t94 + t55) * t210) + m(6) * ((t215 * t85 + t31) * t212 + (-t215 * t86 + t32) * t210); m(6) * ((t215 * t87 + t33) * t212 + (-t215 * t88 + t34) * t210) + m(5) * ((t215 * t97 + t56) * t212 + (-t215 * t98 + t57) * t210) + m(4) * ((t137 * t215 + t103) * t212 + (-t138 * t215 + t104) * t210); 0; t226 + m(5) * (t51 * t94 + t52 * t93 + t54 * t96 + t55 * t95) + m(6) * (t16 * t86 + t17 * t85 + t31 * t48 + t32 * t47); m(6) * (t16 * t88 + t17 * t87 + t33 * t48 + t34 * t47) + m(5) * (t51 * t98 + t52 * t97 + t56 * t96 + t57 * t95) + t226; m(5) * ((t215 * t95 + t51) * t212 + (-t215 * t96 + t52) * t210) + m(6) * ((t215 * t47 + t16) * t212 + (-t215 * t48 + t17) * t210); (t96 * t51 + t95 * t52 + (-t127 * t212 - t128 * t210) * ((-t128 * t215 - t82) * t212 + (t127 * t215 - t81) * t210) * t213) * t313 - t219 * (-t72 * t219 + ((-t215 * t49 + t14) * t212 + (-t215 * t50 - t13) * t210) * t218) + (-t24 * t219 + (-(t131 * t121 + t132 * t123 + t174 * t80 - t237 * t78) * t210 - t41 * t294 + (t131 * t122 + t132 * t124 + t174 * t79 - t237 * t77) * t212 - t42 * t298 + (-(-t119 * t298 + t212 * t76) * t210 + (-t120 * t298 + t212 * t75) * t212) * t218) * t218) * t293 + (t16 * t48 + t17 * t47 + t26 * t4) * t312 + t260 + (t25 * t219 - (-(t133 * t121 + t134 * t123 + t171 * t78 - t236 * t80) * t210 - t39 * t294 + (t133 * t122 + t134 * t124 + t171 * t77 - t236 * t79) * t212 - t40 * t298 + (-(-t119 * t294 - t210 * t76) * t210 + (-t120 * t294 - t210 * t75) * t212) * t218) * t218 - t2) * t297 + (t84 * t219 - (-t210 * t41 + t212 * t42) * t218 - t6) * t266 + (t83 * t219 - (-t210 * t39 + t212 * t40) * t218 - t5) * t264; m(6) * (t31 * t90 + t32 * t89 + t37 * t86 + t38 * t85) + t233; m(6) * (t33 * t90 + t34 * t89 + t37 * t88 + t38 * t87) + t233; m(6) * ((t215 * t89 + t37) * t212 + (-t215 * t90 + t38) * t210); m(6) * (t15 * t26 + t16 * t90 + t17 * t89 + t37 * t48 + t38 * t47 + t4 * t67) + t227; (t15 * t67 + t37 * t90 + t38 * t89) * t312 + t227;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
