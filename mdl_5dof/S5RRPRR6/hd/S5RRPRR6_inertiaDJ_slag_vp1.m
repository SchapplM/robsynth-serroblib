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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:16:51
% EndTime: 2022-01-20 11:17:02
% DurationCPUTime: 4.93s
% Computational Cost: add. (16061->516), mult. (15018->732), div. (0->0), fcn. (14002->10), ass. (0->252)
t220 = qJD(1) + qJD(2);
t227 = cos(qJ(4));
t272 = qJD(4) * t227;
t223 = sin(pkin(9));
t229 = -pkin(8) - pkin(7);
t285 = t223 * t229;
t316 = pkin(4) * t272 + t220 * t285;
t221 = qJ(4) + qJ(5);
t213 = sin(t221);
t222 = qJ(1) + qJ(2);
t216 = cos(t222);
t224 = cos(pkin(9));
t289 = t216 * t224;
t214 = sin(t222);
t215 = cos(t221);
t297 = t214 * t215;
t161 = -t213 * t289 + t297;
t299 = t213 * t214;
t162 = t215 * t289 + t299;
t290 = t216 * t223;
t115 = t162 * rSges(6,1) + t161 * rSges(6,2) + rSges(6,3) * t290;
t225 = sin(qJ(4));
t293 = t214 * t225;
t202 = pkin(4) * t293;
t212 = pkin(4) * t227 + pkin(3);
t315 = -t212 * t289 + t216 * t285 - t115 - t202;
t307 = pkin(3) - t212;
t314 = pkin(7) * t223 + t224 * t307;
t313 = 2 * m(3);
t312 = 2 * m(4);
t311 = 2 * m(5);
t310 = 2 * m(6);
t218 = t223 ^ 2;
t226 = sin(qJ(1));
t309 = pkin(1) * t226;
t288 = t216 * t225;
t203 = pkin(4) * t288;
t273 = qJD(4) * t225;
t269 = pkin(4) * t273;
t254 = t224 * t269;
t236 = t220 * t203 + t316 * t214 - t216 * t254;
t296 = t214 * t220;
t286 = t220 * t223;
t267 = t214 * t286;
t219 = qJD(4) + qJD(5);
t256 = -t219 * t224 + t220;
t245 = t256 * t216;
t257 = t220 * t224 - t219;
t116 = t215 * t245 + t257 * t299;
t117 = t213 * t245 - t257 * t297;
t281 = t117 * rSges(6,1) + t116 * rSges(6,2);
t69 = -rSges(6,3) * t267 + t281;
t306 = -t296 * t314 - t236 - t69;
t305 = pkin(1) * qJD(1);
t287 = t219 * t223;
t301 = Icges(6,4) * t213;
t146 = (-Icges(6,2) * t215 - t301) * t287;
t157 = -Icges(6,5) * t224 + (Icges(6,1) * t215 - t301) * t223;
t300 = Icges(6,4) * t215;
t156 = -Icges(6,6) * t224 + (-Icges(6,2) * t213 + t300) * t223;
t265 = t219 * t215 * t156;
t233 = -t265 + (-t157 * t219 - t146) * t213;
t147 = (-Icges(6,1) * t213 - t300) * t287;
t142 = t223 * t215 * t147;
t145 = (-Icges(6,5) * t213 - Icges(6,6) * t215) * t287;
t259 = -t224 * t145 + t142;
t53 = t223 * t233 + t259;
t304 = t53 * t224;
t303 = Icges(5,4) * t225;
t302 = Icges(5,4) * t227;
t298 = t213 * t216;
t206 = t214 * qJ(3);
t295 = t214 * t223;
t294 = t214 * t224;
t292 = t215 * t216;
t291 = t216 * t220;
t284 = t224 * t225;
t283 = t224 * t227;
t274 = qJD(4) * t223;
t176 = (-Icges(5,2) * t227 - t303) * t274;
t282 = t225 * t176;
t276 = pkin(3) * t289 + pkin(7) * t290;
t280 = t276 + t315;
t174 = t216 * t283 + t293;
t241 = t214 * t284 + t216 * t227;
t132 = -qJD(4) * t174 + t220 * t241;
t172 = t214 * t283 - t288;
t173 = t214 * t227 - t216 * t284;
t133 = qJD(4) * t173 - t172 * t220;
t279 = t133 * rSges(5,1) + t132 * rSges(5,2);
t159 = -t213 * t294 - t292;
t160 = t215 * t294 - t298;
t247 = -rSges(6,1) * t160 - rSges(6,2) * t159;
t114 = rSges(6,3) * t295 - t247;
t158 = -rSges(6,3) * t224 + (rSges(6,1) * t215 - rSges(6,2) * t213) * t223;
t144 = t158 * t295;
t90 = t224 * t114 + t144;
t278 = -t214 * t285 - t203;
t277 = qJ(3) * t291 + qJD(3) * t214;
t275 = t216 * pkin(2) + t206;
t271 = t226 * t305;
t228 = cos(qJ(1));
t270 = t228 * t305;
t148 = (-rSges(6,1) * t213 - rSges(6,2) * t215) * t287;
t266 = t216 * t286;
t246 = t256 * t214;
t118 = t215 * t246 - t257 * t298;
t119 = t213 * t246 + t257 * t292;
t248 = -rSges(6,1) * t119 - rSges(6,2) * t118;
t70 = rSges(6,3) * t266 - t248;
t38 = t148 * t295 + t158 * t266 + t224 * t70;
t129 = t174 * rSges(5,1) + t173 * rSges(5,2) + rSges(5,3) * t290;
t263 = -t214 * t254 - t316 * t216;
t262 = t295 / 0.2e1;
t261 = t290 / 0.2e1;
t260 = -rSges(4,1) * t224 - pkin(2);
t185 = t216 * rSges(3,1) - rSges(3,2) * t214;
t177 = (-Icges(5,1) * t225 - t302) * t274;
t154 = t223 * t227 * t177;
t175 = (-Icges(5,5) * t225 - Icges(5,6) * t227) * t274;
t258 = -t224 * t175 + t154;
t255 = t218 * t269;
t253 = -t267 / 0.2e1;
t252 = t220 * t261;
t169 = -rSges(3,1) * t291 + rSges(3,2) * t296;
t251 = t260 * t214;
t184 = -rSges(3,1) * t214 - rSges(3,2) * t216;
t134 = -qJD(4) * t172 + t173 * t220;
t135 = -qJD(4) * t241 + t174 * t220;
t250 = -rSges(5,1) * t135 - rSges(5,2) * t134;
t249 = -rSges(5,1) * t172 + rSges(5,2) * t241;
t243 = -rSges(6,3) * t223 - t212 * t224 - pkin(2);
t168 = t184 * t220;
t111 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t290;
t113 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t290;
t62 = Icges(6,5) * t117 + Icges(6,6) * t116 - Icges(6,3) * t267;
t64 = Icges(6,4) * t117 + Icges(6,2) * t116 - Icges(6,6) * t267;
t66 = Icges(6,1) * t117 + Icges(6,4) * t116 - Icges(6,5) * t267;
t10 = -t224 * t62 + ((-t111 * t219 + t66) * t215 + (-t113 * t219 - t64) * t213) * t223;
t155 = -Icges(6,3) * t224 + (Icges(6,5) * t215 - Icges(6,6) * t213) * t223;
t20 = t116 * t156 + t117 * t157 + t146 * t161 + t147 * t162 + (t145 * t216 - t155 * t296) * t223;
t21 = t118 * t156 + t119 * t157 + t146 * t159 + t147 * t160 + (t145 * t214 + t155 * t291) * t223;
t108 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t295;
t110 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t295;
t112 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t295;
t43 = -t108 * t224 + (-t110 * t213 + t112 * t215) * t223;
t109 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t290;
t44 = -t109 * t224 + (-t111 * t213 + t113 * t215) * t223;
t71 = t155 * t295 + t156 * t159 + t157 * t160;
t72 = t155 * t290 + t156 * t161 + t157 * t162;
t63 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t266;
t65 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t266;
t67 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t266;
t9 = -t224 * t63 + ((-t110 * t219 + t67) * t215 + (-t112 * t219 - t65) * t213) * t223;
t240 = (t21 + t9) * t262 + (t10 + t20) * t261 + (t44 + t72) * t253 + (t43 + t71) * t252;
t239 = -pkin(3) * t224 - pkin(2) + (-rSges(5,3) - pkin(7)) * t223;
t99 = t129 + t275 + t276;
t238 = t243 * t214;
t139 = rSges(4,1) * t289 - rSges(4,2) * t290 + t214 * rSges(4,3) + t275;
t237 = t239 * t214;
t207 = t216 * qJ(3);
t138 = rSges(4,2) * t295 + t216 * rSges(4,3) + t207 + t251;
t27 = t108 * t295 + t110 * t159 + t112 * t160;
t28 = t109 * t295 + t111 * t159 + t113 * t160;
t29 = t108 * t290 + t110 * t161 + t112 * t162;
t30 = t109 * t290 + t111 * t161 + t113 * t162;
t235 = -t224 * (-t304 + ((t220 * t43 + t10) * t216 + (-t220 * t44 + t9) * t214) * t223) - t267 * (-t224 * t72 + (t214 * t29 + t216 * t30) * t223) + (-t21 * t224 + ((t118 * t111 + t119 * t113 + t159 * t64 + t160 * t66) * t216 - t28 * t296 + (t118 * t110 + t119 * t112 + t159 * t65 + t160 * t67) * t214 + t27 * t291 + ((t109 * t291 + t214 * t62) * t216 + (t108 * t291 + t214 * t63) * t214) * t223) * t223) * t295 + (-t20 * t224 + ((t116 * t111 + t117 * t113 + t161 * t64 + t162 * t66) * t216 - t30 * t296 + (t110 * t116 + t112 * t117 + t161 * t65 + t162 * t67) * t214 + t29 * t291 + ((-t109 * t296 + t216 * t62) * t216 + (-t108 * t296 + t216 * t63) * t214) * t223) * t223) * t290 + (-t224 * t71 + (t214 * t27 + t216 * t28) * t223) * t266;
t234 = t240 - t304;
t166 = -Icges(5,6) * t224 + (-Icges(5,2) * t225 + t302) * t223;
t167 = -Icges(5,5) * t224 + (Icges(5,1) * t227 - t303) * t223;
t232 = -t282 + (-t166 * t227 - t167 * t225) * qJD(4);
t104 = rSges(4,2) * t267 + rSges(4,3) * t291 + t220 * t251 + t277;
t89 = t275 - t315;
t205 = qJD(3) * t216;
t105 = rSges(4,2) * t266 + t205 + (t260 * t216 + (-rSges(4,3) - qJ(3)) * t214) * t220;
t56 = t220 * t237 + t277 + t279;
t98 = t207 + t237 + t249;
t88 = t207 + t238 + t247 - t278;
t122 = Icges(5,4) * t172 - Icges(5,2) * t241 + Icges(5,6) * t295;
t124 = Icges(5,1) * t172 - Icges(5,4) * t241 + Icges(5,5) * t295;
t77 = Icges(5,5) * t135 + Icges(5,6) * t134 + Icges(5,3) * t266;
t79 = Icges(5,4) * t135 + Icges(5,2) * t134 + Icges(5,6) * t266;
t81 = Icges(5,1) * t135 + Icges(5,4) * t134 + Icges(5,5) * t266;
t13 = -t224 * t77 + (-t225 * t79 + t227 * t81 + (-t122 * t227 - t124 * t225) * qJD(4)) * t223;
t123 = Icges(5,4) * t174 + Icges(5,2) * t173 + Icges(5,6) * t290;
t125 = Icges(5,1) * t174 + Icges(5,4) * t173 + Icges(5,5) * t290;
t76 = Icges(5,5) * t133 + Icges(5,6) * t132 - Icges(5,3) * t267;
t78 = Icges(5,4) * t133 + Icges(5,2) * t132 - Icges(5,6) * t267;
t80 = Icges(5,1) * t133 + Icges(5,4) * t132 - Icges(5,5) * t267;
t14 = -t224 * t76 + (-t225 * t78 + t227 * t80 + (-t123 * t227 - t125 * t225) * qJD(4)) * t223;
t165 = -Icges(5,3) * t224 + (Icges(5,5) * t227 - Icges(5,6) * t225) * t223;
t24 = t132 * t166 + t133 * t167 + t173 * t176 + t174 * t177 + (-t165 * t296 + t175 * t216) * t223;
t25 = t134 * t166 + t135 * t167 - t241 * t176 + t172 * t177 + (t165 * t291 + t175 * t214) * t223;
t120 = Icges(5,5) * t172 - Icges(5,6) * t241 + Icges(5,3) * t295;
t51 = -t120 * t224 + (-t122 * t225 + t124 * t227) * t223;
t121 = Icges(5,5) * t174 + Icges(5,6) * t173 + Icges(5,3) * t290;
t52 = -t121 * t224 + (-t123 * t225 + t125 * t227) * t223;
t73 = t223 * t232 + t258;
t84 = t165 * t295 - t166 * t241 + t167 * t172;
t85 = t165 * t290 + t166 * t173 + t167 * t174;
t231 = (-t53 - t73) * t224 + t240 + (t25 + t13) * t262 + (t24 + t14) * t261 + (t52 + t85) * t253 + (t51 + t84) * t252;
t57 = t205 + (t216 * t239 - t206) * t220 + t250;
t33 = t220 * t238 + t236 + t277 + t281;
t34 = t205 + ((-pkin(4) * t225 - qJ(3)) * t214 + t243 * t216) * t220 + t248 - t263;
t230 = t142 + t154 + (-t145 - t175) * t224 + (t232 + t233) * t223;
t217 = t228 * pkin(1);
t180 = t185 + t217;
t179 = t184 - t309;
t178 = (-rSges(5,1) * t225 - rSges(5,2) * t227) * t274;
t170 = -rSges(5,3) * t224 + (rSges(5,1) * t227 - rSges(5,2) * t225) * t223;
t151 = (pkin(7) + t229) * t224 - t307 * t223;
t150 = t169 - t270;
t149 = t168 - t271;
t140 = t220 * t144;
t137 = t139 + t217;
t136 = t138 - t309;
t128 = rSges(5,3) * t295 - t249;
t126 = -t214 * t314 + t278;
t102 = t105 - t270;
t101 = t104 - t271;
t100 = t114 * t290;
t97 = -t129 * t224 - t170 * t290;
t96 = t128 * t224 + t170 * t295;
t95 = t217 + t99;
t94 = t98 - t309;
t93 = (-t216 * t314 + t202) * t220 + t263;
t91 = -t115 * t224 - t158 * t290;
t87 = t217 + t89;
t86 = t88 - t309;
t83 = rSges(5,3) * t266 - t250;
t82 = -rSges(5,3) * t267 + t279;
t68 = -t115 * t295 + t100;
t58 = t70 * t290;
t55 = t57 - t270;
t54 = t56 - t271;
t50 = t224 * t83 + (t170 * t291 + t178 * t214) * t223;
t49 = -t224 * t82 + (t170 * t296 - t178 * t216) * t223;
t48 = t280 * t224 + (-t151 - t158) * t290;
t47 = t126 * t224 + t151 * t295 + t90;
t42 = t121 * t290 + t123 * t173 + t125 * t174;
t41 = t120 * t290 + t122 * t173 + t124 * t174;
t40 = t121 * t295 - t123 * t241 + t125 * t172;
t39 = t120 * t295 - t122 * t241 + t124 * t172;
t37 = -t148 * t290 - t224 * t69 + t140;
t32 = t34 - t270;
t31 = t33 - t271;
t26 = t100 + (t126 * t216 + t214 * t280) * t223;
t17 = t151 * t266 - t214 * t255 + t224 * t93 + t38;
t16 = t151 * t267 + t140 + t306 * t224 + (-t148 * t223 + t255) * t216;
t15 = t58 + (-t115 * t291 + (-t114 * t220 - t69) * t214) * t223;
t4 = t58 + ((t220 * t280 + t93) * t216 + ((-t114 - t126) * t220 + t306) * t214) * t223;
t1 = [(t31 * t87 + t32 * t86) * t310 - t213 * t157 * t287 + (t54 * t95 + t55 * t94) * t311 + (t101 * t137 + t102 * t136) * t312 + (t149 * t180 + t150 * t179) * t313 + t258 + t259 + (-t213 * t146 - t166 * t272 - t167 * t273 - t265 - t282) * t223; m(6) * (t31 * t89 + t32 * t88 + t33 * t87 + t34 * t86) + m(5) * (t54 * t99 + t55 * t98 + t56 * t95 + t57 * t94) + m(4) * (t101 * t139 + t102 * t138 + t104 * t137 + t105 * t136) + m(3) * (t149 * t185 + t150 * t184 + t168 * t180 + t169 * t179) + t230; (t33 * t89 + t34 * t88) * t310 + (t56 * t99 + t57 * t98) * t311 + (t104 * t139 + t105 * t138) * t312 + (t168 * t185 + t169 * t184) * t313 + t230; m(6) * ((t220 * t86 - t31) * t216 + (t220 * t87 + t32) * t214) + m(5) * ((t220 * t94 - t54) * t216 + (t220 * t95 + t55) * t214) + m(4) * ((t136 * t220 - t101) * t216 + (t137 * t220 + t102) * t214); m(6) * ((t220 * t88 - t33) * t216 + (t220 * t89 + t34) * t214) + m(5) * ((t220 * t98 - t56) * t216 + (t220 * t99 + t57) * t214) + m(4) * ((t138 * t220 - t104) * t216 + (t139 * t220 + t105) * t214); 0; m(6) * (t16 * t87 + t17 * t86 + t31 * t48 + t32 * t47) + m(5) * (t49 * t95 + t50 * t94 + t54 * t97 + t55 * t96) + t231; m(6) * (t16 * t89 + t17 * t88 + t33 * t48 + t34 * t47) + m(5) * (t49 * t99 + t50 * t98 + t56 * t97 + t57 * t96) + t231; m(5) * ((t220 * t96 - t49) * t216 + (t220 * t97 + t50) * t214) + m(6) * ((t220 * t47 - t16) * t216 + (t220 * t48 + t17) * t214); (t97 * t49 + t96 * t50 + (t128 * t216 - t129 * t214) * ((-t129 * t220 + t83) * t216 + (-t128 * t220 - t82) * t214) * t218) * t311 - (-t85 * t224 + (t214 * t41 + t216 * t42) * t223) * t267 + (-t24 * t224 + ((t132 * t123 + t133 * t125 + t173 * t78 + t174 * t80) * t216 - t42 * t296 + (t132 * t122 + t133 * t124 + t173 * t79 + t174 * t81) * t214 + t41 * t291 + ((-t121 * t296 + t216 * t76) * t216 + (-t120 * t296 + t216 * t77) * t214) * t223) * t223) * t290 + (-t84 * t224 + (t214 * t39 + t216 * t40) * t223) * t266 + (-t25 * t224 + ((t134 * t123 + t135 * t125 + t172 * t80 - t241 * t78) * t216 - t40 * t296 + (t134 * t122 + t135 * t124 + t172 * t81 - t241 * t79) * t214 + t39 * t291 + ((t121 * t291 + t214 * t76) * t216 + (t120 * t291 + t214 * t77) * t214) * t223) * t223) * t295 - t224 * (-t73 * t224 + ((t220 * t51 + t14) * t216 + (-t220 * t52 + t13) * t214) * t223) + (t16 * t48 + t17 * t47 + t26 * t4) * t310 + t235; m(6) * (t31 * t91 + t32 * t90 + t37 * t87 + t38 * t86) + t234; m(6) * (t33 * t91 + t34 * t90 + t37 * t89 + t38 * t88) + t234; m(6) * ((t220 * t90 - t37) * t216 + (t220 * t91 + t38) * t214); m(6) * (t15 * t26 + t16 * t91 + t17 * t90 + t37 * t48 + t38 * t47 + t4 * t68) + t235; (t15 * t68 + t37 * t91 + t38 * t90) * t310 + t235;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
