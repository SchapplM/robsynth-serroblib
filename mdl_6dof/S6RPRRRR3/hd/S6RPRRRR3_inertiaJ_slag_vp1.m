% Calculate joint inertia matrix for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:23
% EndTime: 2019-03-09 07:00:32
% DurationCPUTime: 3.64s
% Computational Cost: add. (16426->462), mult. (13082->638), div. (0->0), fcn. (14198->12), ass. (0->241)
t227 = sin(qJ(3));
t312 = Icges(4,5) * t227;
t311 = t312 / 0.2e1;
t225 = qJ(4) + qJ(5);
t220 = qJ(6) + t225;
t213 = sin(t220);
t214 = cos(t220);
t223 = qJ(1) + pkin(11);
t216 = sin(t223);
t217 = cos(t223);
t230 = cos(qJ(3));
t280 = t217 * t230;
t155 = -t213 * t280 + t214 * t216;
t156 = t213 * t216 + t214 * t280;
t281 = t217 * t227;
t109 = t156 * rSges(7,1) + t155 * rSges(7,2) + rSges(7,3) * t281;
t229 = cos(qJ(4));
t215 = t229 * pkin(4) + pkin(3);
t219 = cos(t225);
t190 = pkin(5) * t219 + t215;
t218 = sin(t225);
t226 = sin(qJ(4));
t191 = pkin(4) * t226 + pkin(5) * t218;
t310 = t190 * t280 + t216 * t191 + t109;
t309 = t216 ^ 2;
t308 = t217 ^ 2;
t232 = -pkin(9) - pkin(8);
t285 = t216 * t227;
t284 = t216 * t230;
t153 = -t213 * t284 - t214 * t217;
t154 = -t213 * t217 + t214 * t284;
t102 = Icges(7,5) * t154 + Icges(7,6) * t153 + Icges(7,3) * t285;
t104 = Icges(7,4) * t154 + Icges(7,2) * t153 + Icges(7,6) * t285;
t106 = Icges(7,1) * t154 + Icges(7,4) * t153 + Icges(7,5) * t285;
t35 = t102 * t285 + t104 * t153 + t106 * t154;
t103 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t281;
t105 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t281;
t107 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t281;
t36 = t103 * t285 + t105 * t153 + t107 * t154;
t157 = -Icges(7,3) * t230 + (Icges(7,5) * t214 - Icges(7,6) * t213) * t227;
t158 = -Icges(7,6) * t230 + (Icges(7,4) * t214 - Icges(7,2) * t213) * t227;
t159 = -Icges(7,5) * t230 + (Icges(7,1) * t214 - Icges(7,4) * t213) * t227;
t66 = t153 * t158 + t154 * t159 + t157 * t285;
t5 = -t66 * t230 + (t216 * t35 + t217 * t36) * t227;
t37 = t102 * t281 + t104 * t155 + t106 * t156;
t38 = t103 * t281 + t105 * t155 + t107 * t156;
t67 = t155 * t158 + t156 * t159 + t157 * t281;
t6 = -t67 * t230 + (t216 * t37 + t217 * t38) * t227;
t307 = t6 * t281 + t5 * t285;
t306 = t216 / 0.2e1;
t305 = -t217 / 0.2e1;
t304 = t217 / 0.2e1;
t303 = -t230 / 0.2e1;
t197 = rSges(4,1) * t227 + rSges(4,2) * t230;
t302 = m(4) * t197;
t301 = pkin(3) * t230;
t300 = pkin(8) * t227;
t228 = sin(qJ(1));
t299 = t228 * pkin(1);
t298 = -pkin(3) + t215;
t140 = t227 * t214 * t159;
t289 = t158 * t213;
t84 = -t230 * t157 - t227 * t289 + t140;
t293 = t84 * t230;
t48 = -t230 * t102 + (-t104 * t213 + t106 * t214) * t227;
t49 = -t230 * t103 + (-t105 * t213 + t107 * t214) * t227;
t13 = -t293 + (t216 * t48 + t217 * t49) * t227;
t279 = t218 * t230;
t167 = -t216 * t279 - t217 * t219;
t278 = t219 * t230;
t168 = t216 * t278 - t217 * t218;
t112 = Icges(6,5) * t168 + Icges(6,6) * t167 + Icges(6,3) * t285;
t114 = Icges(6,4) * t168 + Icges(6,2) * t167 + Icges(6,6) * t285;
t116 = Icges(6,1) * t168 + Icges(6,4) * t167 + Icges(6,5) * t285;
t56 = -t230 * t112 + (-t114 * t218 + t116 * t219) * t227;
t169 = t216 * t219 - t217 * t279;
t170 = t216 * t218 + t217 * t278;
t113 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t281;
t115 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t281;
t117 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t281;
t57 = -t230 * t113 + (-t115 * t218 + t117 * t219) * t227;
t173 = -Icges(6,5) * t230 + (Icges(6,1) * t219 - Icges(6,4) * t218) * t227;
t143 = t227 * t219 * t173;
t171 = -Icges(6,3) * t230 + (Icges(6,5) * t219 - Icges(6,6) * t218) * t227;
t172 = -Icges(6,6) * t230 + (Icges(6,4) * t219 - Icges(6,2) * t218) * t227;
t288 = t172 * t218;
t88 = -t230 * t171 - t227 * t288 + t143;
t297 = -t13 + t88 * t230 - (t216 * t56 + t217 * t57) * t227;
t296 = -t84 - t88;
t275 = t227 * t232;
t282 = t217 * t226;
t264 = pkin(4) * t282 + t216 * t275;
t265 = t190 - t215;
t224 = -pkin(10) + t232;
t277 = t224 * t227;
t283 = t217 * t191;
t94 = -t283 + (t230 * t265 - t277) * t216 + t264;
t243 = -t154 * rSges(7,1) - t153 * rSges(7,2);
t108 = rSges(7,3) * t285 - t243;
t97 = t108 * t281;
t295 = t94 * t281 + t97;
t294 = t217 * rSges(4,3);
t262 = t224 - t232;
t286 = t216 * t226;
t266 = -pkin(4) * t286 - t215 * t280;
t292 = -t262 * t281 + t266 + t310;
t290 = Icges(4,4) * t230;
t183 = -Icges(5,6) * t230 + (Icges(5,4) * t229 - Icges(5,2) * t226) * t227;
t287 = t183 * t226;
t276 = t226 * t230;
t274 = t229 * t230;
t120 = t170 * rSges(6,1) + t169 * rSges(6,2) + rSges(6,3) * t281;
t236 = -t217 * t275 - t266;
t263 = pkin(3) * t280 + pkin(8) * t281;
t133 = t236 - t263;
t273 = -t120 - t133;
t132 = (t230 * t298 - t300) * t216 - t264;
t166 = (pkin(8) + t232) * t230 + t298 * t227;
t272 = t230 * t132 + t166 * t285;
t139 = t227 * t265 + t230 * t262;
t161 = -t230 * rSges(7,3) + (rSges(7,1) * t214 - rSges(7,2) * t213) * t227;
t271 = -t139 - t161;
t82 = t230 * t108 + t161 * t285;
t244 = -t168 * rSges(6,1) - t167 * rSges(6,2);
t119 = rSges(6,3) * t285 - t244;
t174 = -t230 * rSges(6,3) + (rSges(6,1) * t219 - rSges(6,2) * t218) * t227;
t85 = t230 * t119 + t174 * t285;
t270 = t309 * (t300 + t301) + t217 * t263;
t269 = -t166 - t174;
t185 = -t230 * rSges(5,3) + (rSges(5,1) * t229 - rSges(5,2) * t226) * t227;
t204 = t227 * pkin(3) - t230 * pkin(8);
t267 = -t185 - t204;
t180 = t216 * t229 - t217 * t276;
t181 = t217 * t274 + t286;
t125 = Icges(5,5) * t181 + Icges(5,6) * t180 + Icges(5,3) * t281;
t127 = Icges(5,4) * t181 + Icges(5,2) * t180 + Icges(5,6) * t281;
t129 = Icges(5,1) * t181 + Icges(5,4) * t180 + Icges(5,5) * t281;
t61 = -t230 * t125 + (-t127 * t226 + t129 * t229) * t227;
t182 = -Icges(5,3) * t230 + (Icges(5,5) * t229 - Icges(5,6) * t226) * t227;
t184 = -Icges(5,5) * t230 + (Icges(5,1) * t229 - Icges(5,4) * t226) * t227;
t81 = t180 * t183 + t181 * t184 + t182 * t281;
t261 = t61 / 0.2e1 + t81 / 0.2e1;
t178 = -t216 * t276 - t217 * t229;
t179 = t216 * t274 - t282;
t124 = Icges(5,5) * t179 + Icges(5,6) * t178 + Icges(5,3) * t285;
t126 = Icges(5,4) * t179 + Icges(5,2) * t178 + Icges(5,6) * t285;
t128 = Icges(5,1) * t179 + Icges(5,4) * t178 + Icges(5,5) * t285;
t60 = -t230 * t124 + (-t126 * t226 + t128 * t229) * t227;
t80 = t178 * t183 + t179 * t184 + t182 * t285;
t260 = t80 / 0.2e1 + t60 / 0.2e1;
t259 = -t133 - t292;
t40 = t112 * t285 + t114 * t167 + t116 * t168;
t41 = t113 * t285 + t115 * t167 + t117 * t168;
t73 = t167 * t172 + t168 * t173 + t171 * t285;
t11 = -t73 * t230 + (t216 * t40 + t217 * t41) * t227;
t42 = t112 * t281 + t114 * t169 + t116 * t170;
t43 = t113 * t281 + t115 * t169 + t117 * t170;
t74 = t169 * t172 + t170 * t173 + t171 * t281;
t12 = -t74 * t230 + (t216 * t42 + t217 * t43) * t227;
t258 = t11 * t285 + t12 * t281 + t307;
t257 = -t166 + t271;
t256 = -t204 + t269;
t131 = t181 * rSges(5,1) + t180 * rSges(5,2) + rSges(5,3) * t281;
t231 = cos(qJ(1));
t222 = t231 * pkin(1);
t255 = t217 * pkin(2) + t216 * pkin(7) + t222;
t254 = t285 / 0.2e1;
t253 = t281 / 0.2e1;
t252 = t217 * pkin(7) - t299;
t251 = (t48 + t66) * t254 + (t49 + t67) * t253;
t44 = t139 * t285 + t230 * t94 + t82;
t250 = -t230 * t13 + t307;
t249 = t216 * t132 + t217 * t133 + t270;
t248 = -t204 + t257;
t18 = t216 * t36 - t217 * t35;
t19 = t216 * t38 - t217 * t37;
t247 = t18 * t254 + t19 * t253 + t5 * t305 + t6 * t306 + (t49 * t216 - t48 * t217) * t303;
t246 = rSges(4,1) * t230 - rSges(4,2) * t227;
t245 = -t179 * rSges(5,1) - t178 * rSges(5,2);
t241 = -Icges(4,2) * t227 + t290;
t240 = Icges(4,5) * t230 - Icges(4,6) * t227;
t237 = rSges(4,1) * t280 - rSges(4,2) * t281 + t216 * rSges(4,3);
t235 = t230 * t297 + t258;
t234 = t251 + (t56 + t73) * t254 + (t57 + t74) * t253;
t23 = t216 * t41 - t217 * t40;
t24 = t216 * t43 - t217 * t42;
t233 = t11 * t305 + t12 * t306 + t23 * t254 + t24 * t253 + t247 + (t57 * t216 - t56 * t217) * t303;
t199 = rSges(2,1) * t231 - rSges(2,2) * t228;
t198 = -rSges(2,1) * t228 - rSges(2,2) * t231;
t194 = Icges(4,6) * t230 + t312;
t187 = rSges(3,1) * t217 - rSges(3,2) * t216 + t222;
t186 = -rSges(3,1) * t216 - rSges(3,2) * t217 - t299;
t160 = t227 * t229 * t184;
t148 = Icges(4,3) * t216 + t217 * t240;
t147 = -Icges(4,3) * t217 + t216 * t240;
t138 = t267 * t217;
t137 = t267 * t216;
t136 = t237 + t255;
t135 = t294 + (-pkin(2) - t246) * t216 + t252;
t130 = rSges(5,3) * t285 - t245;
t118 = t132 * t281;
t111 = t217 * t237 + (t216 * t246 - t294) * t216;
t100 = t119 * t281;
t99 = t256 * t217;
t98 = t256 * t216;
t96 = -t230 * t182 - t227 * t287 + t160;
t93 = -t131 * t230 - t185 * t281;
t92 = t130 * t230 + t185 * t285;
t90 = t255 + t131 + t263;
t89 = (-t301 - pkin(2) + (-rSges(5,3) - pkin(8)) * t227) * t216 + t245 + t252;
t86 = -t230 * t120 - t174 * t281;
t83 = -t230 * t109 - t161 * t281;
t79 = t236 + t255 + t120;
t78 = (-rSges(6,3) * t227 - t215 * t230 - pkin(2)) * t216 + t244 + t252 + t264;
t77 = t248 * t217;
t76 = t248 * t216;
t75 = (t130 * t217 - t131 * t216) * t227;
t72 = -t217 * t277 + t255 + t310;
t71 = t283 + (-t190 * t230 - pkin(2) + (-rSges(7,3) + t224) * t227) * t216 + t243 + t252;
t68 = -t120 * t285 + t100;
t65 = -t109 * t285 + t97;
t62 = t130 * t216 + t131 * t217 + t270;
t59 = t230 * t273 + t269 * t281;
t58 = t85 + t272;
t55 = t125 * t281 + t127 * t180 + t129 * t181;
t54 = t124 * t281 + t126 * t180 + t128 * t181;
t53 = t125 * t285 + t127 * t178 + t129 * t179;
t52 = t124 * t285 + t126 * t178 + t128 * t179;
t45 = -t230 * t292 + t271 * t281;
t39 = t273 * t285 + t100 + t118;
t34 = t119 * t216 + t120 * t217 + t249;
t33 = t230 * t259 + t257 * t281;
t32 = t44 + t272;
t31 = -t285 * t292 + t295;
t29 = t216 * t55 - t217 * t54;
t28 = t216 * t53 - t217 * t52;
t27 = t259 * t285 + t118 + t295;
t25 = t292 * t217 + (t108 + t94) * t216 + t249;
t15 = -t81 * t230 + (t216 * t54 + t217 * t55) * t227;
t14 = -t80 * t230 + (t216 * t52 + t217 * t53) * t227;
t1 = [Icges(2,3) + Icges(3,3) + t140 + t143 + t160 + (Icges(4,4) * t227 + Icges(4,2) * t230 - t157 - t171 - t182) * t230 + (Icges(4,1) * t227 - t287 - t288 - t289 + t290) * t227 + m(7) * (t71 ^ 2 + t72 ^ 2) + m(6) * (t78 ^ 2 + t79 ^ 2) + m(5) * (t89 ^ 2 + t90 ^ 2) + m(4) * (t135 ^ 2 + t136 ^ 2) + m(2) * (t198 ^ 2 + t199 ^ 2) + m(3) * (t186 ^ 2 + t187 ^ 2); 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t71 * t77 + t72 * t76) + m(6) * (t78 * t99 + t79 * t98) + m(5) * (t137 * t90 + t138 * t89) + (-t48 / 0.2e1 - t66 / 0.2e1 - t73 / 0.2e1 - t56 / 0.2e1 + (-Icges(4,6) * t217 + t216 * t241) * t303 + t217 * t311 - t135 * t302 + t194 * t304 - t260) * t217 + (t67 / 0.2e1 + t49 / 0.2e1 + t57 / 0.2e1 + t74 / 0.2e1 - t136 * t302 + t230 * (Icges(4,6) * t216 + t217 * t241) / 0.2e1 + t216 * t311 + t194 * t306 + t261) * t216; m(4) * t111 + m(5) * t62 + m(6) * t34 + m(7) * t25; m(7) * (t25 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(6) * (t34 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(5) * (t137 ^ 2 + t138 ^ 2 + t62 ^ 2) + m(4) * (t111 ^ 2 + (t308 + t309) * t197 ^ 2) + (-t308 * t147 - t18 - t23 - t28) * t217 + (t309 * t148 + t19 + t24 + t29 + (-t216 * t147 + t217 * t148) * t217) * t216; (-t96 + t296) * t230 + m(7) * (t32 * t71 + t33 * t72) + m(6) * (t58 * t78 + t59 * t79) + m(5) * (t89 * t92 + t90 * t93) + (t216 * t260 + t217 * t261) * t227 + t234; m(5) * t75 + m(6) * t39 + m(7) * t27; m(7) * (t25 * t27 + t32 * t77 + t33 * t76) + m(6) * (t34 * t39 + t58 * t99 + t59 * t98) + m(5) * (t137 * t93 + t138 * t92 + t62 * t75) + (t28 * t306 + t29 * t304) * t227 + t233 + t15 * t306 + t14 * t305 + (t61 * t216 - t60 * t217) * t303; (t96 * t230 + t297) * t230 + (t216 * t14 + t217 * t15 - t230 * (t216 * t60 + t217 * t61)) * t227 + m(7) * (t27 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t39 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(5) * (t75 ^ 2 + t92 ^ 2 + t93 ^ 2) + t258; t296 * t230 + m(7) * (t44 * t71 + t45 * t72) + m(6) * (t78 * t85 + t79 * t86) + t234; m(6) * t68 + m(7) * t31; m(7) * (t25 * t31 + t44 * t77 + t45 * t76) + m(6) * (t34 * t68 + t85 * t99 + t86 * t98) + t233; m(7) * (t27 * t31 + t32 * t44 + t33 * t45) + m(6) * (t39 * t68 + t58 * t85 + t59 * t86) + t235; m(7) * (t31 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t68 ^ 2 + t85 ^ 2 + t86 ^ 2) + t235; m(7) * (t71 * t82 + t72 * t83) - t293 + t251; m(7) * t65; m(7) * (t25 * t65 + t76 * t83 + t77 * t82) + t247; m(7) * (t27 * t65 + t32 * t82 + t33 * t83) + t250; m(7) * (t31 * t65 + t44 * t82 + t45 * t83) + t250; m(7) * (t65 ^ 2 + t82 ^ 2 + t83 ^ 2) + t250;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
