% Calculate joint inertia matrix for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:56:01
% EndTime: 2019-03-08 21:56:13
% DurationCPUTime: 4.86s
% Computational Cost: add. (30003->592), mult. (46153->856), div. (0->0), fcn. (58560->14), ass. (0->278)
t247 = sin(pkin(6));
t315 = t247 ^ 2;
t246 = sin(pkin(11));
t248 = cos(pkin(11));
t256 = cos(qJ(2));
t249 = cos(pkin(6));
t253 = sin(qJ(2));
t295 = t249 * t253;
t231 = t246 * t256 + t248 * t295;
t280 = qJ(3) + pkin(12);
t242 = sin(t280);
t270 = cos(t280);
t261 = t247 * t270;
t214 = t231 * t242 + t248 * t261;
t314 = t214 / 0.2e1;
t233 = -t246 * t295 + t248 * t256;
t216 = t233 * t242 - t246 * t261;
t313 = t216 / 0.2e1;
t299 = t247 * t253;
t227 = t242 * t299 - t249 * t270;
t312 = t227 / 0.2e1;
t294 = t249 * t256;
t230 = t246 * t253 - t248 * t294;
t311 = t230 / 0.2e1;
t232 = t246 * t294 + t248 * t253;
t310 = t232 / 0.2e1;
t309 = t249 / 0.2e1;
t255 = cos(qJ(3));
t308 = t255 * pkin(3);
t254 = cos(qJ(5));
t307 = t254 * pkin(5);
t251 = sin(qJ(5));
t304 = t230 * t251;
t303 = t232 * t251;
t302 = t246 * t247;
t301 = t247 * t248;
t252 = sin(qJ(3));
t300 = t247 * t252;
t298 = t247 * t255;
t297 = t247 * t256;
t296 = t249 * t252;
t215 = t231 * t270 - t242 * t301;
t115 = pkin(5) * t304 + pkin(10) * t214 + t215 * t307;
t245 = qJ(5) + qJ(6);
t243 = sin(t245);
t244 = cos(t245);
t179 = -t215 * t243 + t230 * t244;
t180 = t215 * t244 + t230 * t243;
t123 = rSges(7,1) * t180 + rSges(7,2) * t179 + rSges(7,3) * t214;
t293 = t115 + t123;
t217 = t233 * t270 + t242 * t302;
t116 = pkin(5) * t303 + pkin(10) * t216 + t217 * t307;
t181 = -t217 * t243 + t232 * t244;
t182 = t217 * t244 + t232 * t243;
t124 = rSges(7,1) * t182 + rSges(7,2) * t181 + rSges(7,3) * t216;
t292 = t116 + t124;
t276 = t248 * t300;
t159 = -pkin(3) * t276 + qJ(4) * t230 + t231 * t308;
t136 = t232 * t159;
t176 = t215 * pkin(4) + t214 * pkin(9);
t291 = t232 * t176 + t136;
t228 = t249 * t242 + t253 * t261;
t277 = t251 * t297;
t137 = -pkin(5) * t277 + pkin(10) * t227 + t228 * t307;
t210 = -t228 * t243 - t244 * t297;
t211 = t228 * t244 - t243 * t297;
t142 = rSges(7,1) * t211 + rSges(7,2) * t210 + rSges(7,3) * t227;
t290 = t137 + t142;
t204 = pkin(3) * t296 + (-qJ(4) * t256 + t253 * t308) * t247;
t289 = t159 * t297 + t230 * t204;
t278 = t246 * t300;
t160 = pkin(3) * t278 + qJ(4) * t232 + t233 * t308;
t213 = t233 * pkin(2) + t232 * pkin(8);
t209 = t249 * t213;
t288 = t249 * t160 + t209;
t158 = rSges(5,1) * t217 - rSges(5,2) * t216 + rSges(5,3) * t232;
t287 = -t158 - t160;
t212 = t231 * pkin(2) + t230 * pkin(8);
t286 = -t159 - t212;
t177 = t217 * pkin(4) + t216 * pkin(9);
t285 = -t160 - t177;
t191 = t228 * rSges(5,1) - t227 * rSges(5,2) - rSges(5,3) * t297;
t284 = -t191 - t204;
t198 = t228 * pkin(4) + t227 * pkin(9);
t283 = -t198 - t204;
t282 = t212 * t302 + t213 * t301;
t281 = -m(5) - m(6) - m(7);
t117 = Icges(7,5) * t180 + Icges(7,6) * t179 + Icges(7,3) * t214;
t119 = Icges(7,4) * t180 + Icges(7,2) * t179 + Icges(7,6) * t214;
t121 = Icges(7,1) * t180 + Icges(7,4) * t179 + Icges(7,5) * t214;
t67 = t117 * t227 + t119 * t210 + t121 * t211;
t118 = Icges(7,5) * t182 + Icges(7,6) * t181 + Icges(7,3) * t216;
t120 = Icges(7,4) * t182 + Icges(7,2) * t181 + Icges(7,6) * t216;
t122 = Icges(7,1) * t182 + Icges(7,4) * t181 + Icges(7,5) * t216;
t68 = t118 * t227 + t120 * t210 + t122 * t211;
t139 = Icges(7,5) * t211 + Icges(7,6) * t210 + Icges(7,3) * t227;
t140 = Icges(7,4) * t211 + Icges(7,2) * t210 + Icges(7,6) * t227;
t141 = Icges(7,1) * t211 + Icges(7,4) * t210 + Icges(7,5) * t227;
t79 = t139 * t227 + t140 * t210 + t141 * t211;
t26 = t214 * t67 + t216 * t68 + t227 * t79;
t57 = t117 * t214 + t119 * t179 + t121 * t180;
t58 = t118 * t214 + t120 * t179 + t122 * t180;
t74 = t139 * t214 + t140 * t179 + t141 * t180;
t7 = t214 * t57 + t216 * t58 + t227 * t74;
t59 = t117 * t216 + t119 * t181 + t121 * t182;
t60 = t118 * t216 + t120 * t181 + t122 * t182;
t75 = t139 * t216 + t140 * t181 + t141 * t182;
t8 = t214 * t59 + t216 * t60 + t227 * t75;
t279 = t214 * t7 + t216 * t8 + t227 * t26;
t186 = -t217 * t251 + t232 * t254;
t187 = t217 * t254 + t303;
t132 = rSges(6,1) * t187 + rSges(6,2) * t186 + rSges(6,3) * t216;
t275 = -t132 + t285;
t218 = -t228 * t251 - t254 * t297;
t219 = t228 * t254 - t277;
t149 = rSges(6,1) * t219 + rSges(6,2) * t218 + rSges(6,3) * t227;
t274 = -t149 + t283;
t273 = t249 * t177 + t288;
t272 = -t176 + t286;
t271 = -t297 / 0.2e1;
t234 = t249 * t255 - t252 * t299;
t235 = t253 * t298 + t296;
t205 = t235 * rSges(4,1) + t234 * rSges(4,2) - rSges(4,3) * t297;
t236 = (pkin(2) * t253 - pkin(8) * t256) * t247;
t269 = (-t205 - t236) * t247;
t268 = t285 - t292;
t267 = t283 - t290;
t266 = t159 * t302 + t160 * t301 + t282;
t265 = t176 * t297 + t230 * t198 + t289;
t264 = (-t236 + t284) * t247;
t13 = t230 * t57 + t232 * t58 - t297 * t74;
t14 = t230 * t59 + t232 * t60 - t297 * t75;
t28 = t230 * t67 + t232 * t68 - t297 * t79;
t263 = t13 * t314 + t14 * t313 + t26 * t271 + t28 * t312 + t8 * t310 + t7 * t311;
t15 = t74 * t249 + (t246 * t58 - t248 * t57) * t247;
t16 = t75 * t249 + (t246 * t60 - t248 * t59) * t247;
t30 = t79 * t249 + (t246 * t68 - t248 * t67) * t247;
t262 = t15 * t314 + t16 * t313 + t26 * t309 + t30 * t312 + t8 * t302 / 0.2e1 - t7 * t301 / 0.2e1;
t260 = (-t236 + t274) * t247;
t259 = t176 * t302 + t177 * t301 + t266;
t258 = (-t236 + t267) * t247;
t229 = t249 * rSges(3,3) + (rSges(3,1) * t253 + rSges(3,2) * t256) * t247;
t226 = Icges(3,5) * t249 + (Icges(3,1) * t253 + Icges(3,4) * t256) * t247;
t225 = Icges(3,6) * t249 + (Icges(3,4) * t253 + Icges(3,2) * t256) * t247;
t224 = Icges(3,3) * t249 + (Icges(3,5) * t253 + Icges(3,6) * t256) * t247;
t223 = t233 * t255 + t278;
t222 = -t233 * t252 + t246 * t298;
t221 = t231 * t255 - t276;
t220 = -t231 * t252 - t248 * t298;
t203 = Icges(4,1) * t235 + Icges(4,4) * t234 - Icges(4,5) * t297;
t202 = Icges(4,4) * t235 + Icges(4,2) * t234 - Icges(4,6) * t297;
t201 = Icges(4,5) * t235 + Icges(4,6) * t234 - Icges(4,3) * t297;
t200 = t233 * rSges(3,1) - t232 * rSges(3,2) + rSges(3,3) * t302;
t199 = t231 * rSges(3,1) - t230 * rSges(3,2) - rSges(3,3) * t301;
t197 = Icges(3,1) * t233 - Icges(3,4) * t232 + Icges(3,5) * t302;
t196 = Icges(3,1) * t231 - Icges(3,4) * t230 - Icges(3,5) * t301;
t195 = Icges(3,4) * t233 - Icges(3,2) * t232 + Icges(3,6) * t302;
t194 = Icges(3,4) * t231 - Icges(3,2) * t230 - Icges(3,6) * t301;
t193 = Icges(3,5) * t233 - Icges(3,6) * t232 + Icges(3,3) * t302;
t192 = Icges(3,5) * t231 - Icges(3,6) * t230 - Icges(3,3) * t301;
t190 = Icges(5,1) * t228 - Icges(5,4) * t227 - Icges(5,5) * t297;
t189 = Icges(5,4) * t228 - Icges(5,2) * t227 - Icges(5,6) * t297;
t188 = Icges(5,5) * t228 - Icges(5,6) * t227 - Icges(5,3) * t297;
t185 = t215 * t254 + t304;
t184 = -t215 * t251 + t230 * t254;
t175 = -t249 * t199 - t229 * t301;
t174 = t249 * t200 - t229 * t302;
t169 = rSges(4,1) * t223 + rSges(4,2) * t222 + rSges(4,3) * t232;
t168 = rSges(4,1) * t221 + rSges(4,2) * t220 + rSges(4,3) * t230;
t167 = Icges(4,1) * t223 + Icges(4,4) * t222 + Icges(4,5) * t232;
t166 = Icges(4,1) * t221 + Icges(4,4) * t220 + Icges(4,5) * t230;
t165 = Icges(4,4) * t223 + Icges(4,2) * t222 + Icges(4,6) * t232;
t164 = Icges(4,4) * t221 + Icges(4,2) * t220 + Icges(4,6) * t230;
t163 = Icges(4,5) * t223 + Icges(4,6) * t222 + Icges(4,3) * t232;
t162 = Icges(4,5) * t221 + Icges(4,6) * t220 + Icges(4,3) * t230;
t157 = rSges(5,1) * t215 - rSges(5,2) * t214 + rSges(5,3) * t230;
t156 = Icges(5,1) * t217 - Icges(5,4) * t216 + Icges(5,5) * t232;
t155 = Icges(5,1) * t215 - Icges(5,4) * t214 + Icges(5,5) * t230;
t154 = Icges(5,4) * t217 - Icges(5,2) * t216 + Icges(5,6) * t232;
t153 = Icges(5,4) * t215 - Icges(5,2) * t214 + Icges(5,6) * t230;
t152 = Icges(5,5) * t217 - Icges(5,6) * t216 + Icges(5,3) * t232;
t151 = Icges(5,5) * t215 - Icges(5,6) * t214 + Icges(5,3) * t230;
t148 = Icges(6,1) * t219 + Icges(6,4) * t218 + Icges(6,5) * t227;
t147 = Icges(6,4) * t219 + Icges(6,2) * t218 + Icges(6,6) * t227;
t146 = Icges(6,5) * t219 + Icges(6,6) * t218 + Icges(6,3) * t227;
t138 = (t199 * t246 + t200 * t248) * t247;
t135 = t214 * t142;
t134 = -t169 * t297 - t232 * t205;
t133 = t168 * t297 + t230 * t205;
t131 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t214;
t130 = Icges(6,1) * t187 + Icges(6,4) * t186 + Icges(6,5) * t216;
t129 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t214;
t128 = Icges(6,4) * t187 + Icges(6,2) * t186 + Icges(6,6) * t216;
t127 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t214;
t126 = Icges(6,5) * t187 + Icges(6,6) * t186 + Icges(6,3) * t216;
t125 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t214;
t114 = -t201 * t297 + t202 * t234 + t203 * t235;
t113 = t227 * t124;
t112 = t168 * t232 - t169 * t230;
t111 = (-t168 - t212) * t249 + t248 * t269;
t110 = t249 * t169 + t246 * t269 + t209;
t109 = t216 * t123;
t108 = -t188 * t297 - t189 * t227 + t190 * t228;
t107 = t201 * t232 + t202 * t222 + t203 * t223;
t106 = t201 * t230 + t202 * t220 + t203 * t221;
t105 = t188 * t232 - t189 * t216 + t190 * t217;
t104 = t188 * t230 - t189 * t214 + t190 * t215;
t103 = (t168 * t246 + t169 * t248) * t247 + t282;
t102 = -t163 * t297 + t165 * t234 + t167 * t235;
t101 = -t162 * t297 + t164 * t234 + t166 * t235;
t100 = t132 * t227 - t149 * t216;
t99 = -t131 * t227 + t149 * t214;
t98 = -t142 * t216 + t113;
t97 = -t123 * t227 + t135;
t96 = -t152 * t297 - t154 * t227 + t156 * t228;
t95 = -t151 * t297 - t153 * t227 + t155 * t228;
t94 = t232 * t284 + t287 * t297;
t93 = t157 * t297 + t191 * t230 + t289;
t92 = t163 * t232 + t165 * t222 + t167 * t223;
t91 = t162 * t232 + t164 * t222 + t166 * t223;
t90 = t163 * t230 + t165 * t220 + t167 * t221;
t89 = t162 * t230 + t164 * t220 + t166 * t221;
t88 = (-t157 + t286) * t249 + t248 * t264;
t87 = t249 * t158 + t246 * t264 + t288;
t86 = t152 * t232 - t154 * t216 + t156 * t217;
t85 = t151 * t232 - t153 * t216 + t155 * t217;
t84 = t152 * t230 - t154 * t214 + t156 * t215;
t83 = t151 * t230 - t153 * t214 + t155 * t215;
t82 = t146 * t227 + t147 * t218 + t148 * t219;
t81 = t131 * t216 - t132 * t214;
t80 = -t124 * t214 + t109;
t78 = t157 * t232 + t230 * t287 + t136;
t77 = t146 * t216 + t147 * t186 + t148 * t187;
t76 = t146 * t214 + t147 * t184 + t148 * t185;
t73 = (t157 * t246 + t158 * t248) * t247 + t266;
t72 = t126 * t227 + t128 * t218 + t130 * t219;
t71 = t125 * t227 + t127 * t218 + t129 * t219;
t70 = t232 * t274 + t275 * t297;
t69 = t131 * t297 + t149 * t230 + t265;
t66 = (-t131 + t272) * t249 + t248 * t260;
t65 = t249 * t132 + t246 * t260 + t273;
t64 = t126 * t216 + t128 * t186 + t130 * t187;
t63 = t125 * t216 + t127 * t186 + t129 * t187;
t62 = t126 * t214 + t128 * t184 + t130 * t185;
t61 = t125 * t214 + t127 * t184 + t129 * t185;
t56 = t116 * t227 - t216 * t290 + t113;
t55 = t137 * t214 - t227 * t293 + t135;
t54 = t131 * t232 + t230 * t275 + t291;
t53 = (t131 * t246 + t132 * t248) * t247 + t259;
t52 = t114 * t249 + (-t101 * t248 + t102 * t246) * t247;
t51 = t101 * t230 + t102 * t232 - t114 * t297;
t50 = t115 * t216 - t214 * t292 + t109;
t49 = t108 * t249 + (t246 * t96 - t248 * t95) * t247;
t48 = t107 * t249 + (t246 * t92 - t248 * t91) * t247;
t47 = t106 * t249 + (t246 * t90 - t248 * t89) * t247;
t46 = t232 * t267 + t268 * t297;
t45 = t230 * t290 + t293 * t297 + t265;
t44 = (t272 - t293) * t249 + t248 * t258;
t43 = t246 * t258 + t249 * t292 + t273;
t42 = -t108 * t297 + t230 * t95 + t232 * t96;
t41 = -t107 * t297 + t230 * t91 + t232 * t92;
t40 = -t106 * t297 + t230 * t89 + t232 * t90;
t39 = t105 * t249 + (t246 * t86 - t248 * t85) * t247;
t38 = t104 * t249 + (t246 * t84 - t248 * t83) * t247;
t37 = -t105 * t297 + t230 * t85 + t232 * t86;
t36 = -t104 * t297 + t230 * t83 + t232 * t84;
t35 = t230 * t268 + t232 * t293 + t291;
t34 = (t246 * t293 + t248 * t292) * t247 + t259;
t33 = t82 * t249 + (t246 * t72 - t248 * t71) * t247;
t32 = t230 * t71 + t232 * t72 - t297 * t82;
t31 = t214 * t71 + t216 * t72 + t227 * t82;
t23 = t77 * t249 + (t246 * t64 - t248 * t63) * t247;
t22 = t76 * t249 + (t246 * t62 - t248 * t61) * t247;
t20 = t230 * t63 + t232 * t64 - t297 * t77;
t19 = t230 * t61 + t232 * t62 - t297 * t76;
t18 = t214 * t63 + t216 * t64 + t227 * t77;
t17 = t214 * t61 + t216 * t62 + t227 * t76;
t1 = [m(2) + m(3) + m(4) - t281; m(3) * t138 + m(4) * t103 + m(5) * t73 + m(6) * t53 + m(7) * t34; m(7) * (t34 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t53 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t73 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t103 ^ 2 + t110 ^ 2 + t111 ^ 2) + m(3) * (t138 ^ 2 + t174 ^ 2 + t175 ^ 2) + (t16 + t23 + t39 + t48 + (t193 * t302 - t195 * t232 + t197 * t233) * t302) * t302 + (-t15 - t22 - t47 - t38 + (-t192 * t301 - t194 * t230 + t196 * t231) * t301 + (-t192 * t302 + t193 * t301 + t194 * t232 + t195 * t230 - t196 * t233 - t197 * t231) * t302) * t301 + ((t224 * t302 - t232 * t225 + t233 * t226) * t302 - (-t224 * t301 - t230 * t225 + t231 * t226) * t301 + t30 + t33 + t52 + t49 + ((t195 * t256 + t197 * t253) * t246 - (t194 * t256 + t196 * t253) * t248) * t315 + ((-t192 * t248 + t193 * t246 + t225 * t256 + t226 * t253) * t247 + t249 * t224) * t249) * t249; m(4) * t112 + m(5) * t78 + m(6) * t54 + m(7) * t35; (t28 / 0.2e1 + t32 / 0.2e1 + t42 / 0.2e1 + t51 / 0.2e1) * t249 + (t16 / 0.2e1 + t23 / 0.2e1 + t39 / 0.2e1 + t48 / 0.2e1) * t232 + (t15 / 0.2e1 + t22 / 0.2e1 + t38 / 0.2e1 + t47 / 0.2e1) * t230 + m(7) * (t34 * t35 + t43 * t46 + t44 * t45) + m(6) * (t53 * t54 + t65 * t70 + t66 * t69) + m(5) * (t73 * t78 + t87 * t94 + t88 * t93) + m(4) * (t103 * t112 + t110 * t134 + t111 * t133) + ((-t49 / 0.2e1 - t52 / 0.2e1 - t30 / 0.2e1 - t33 / 0.2e1) * t256 + (-t36 / 0.2e1 - t40 / 0.2e1 - t13 / 0.2e1 - t19 / 0.2e1) * t248 + (t37 / 0.2e1 + t41 / 0.2e1 + t14 / 0.2e1 + t20 / 0.2e1) * t246) * t247; (-t28 - t32 - t42 - t51) * t297 + (t14 + t20 + t41 + t37) * t232 + (t13 + t19 + t36 + t40) * t230 + m(7) * (t35 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t54 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t78 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(4) * (t112 ^ 2 + t133 ^ 2 + t134 ^ 2); t281 * t297; m(7) * (t230 * t43 + t232 * t44 - t297 * t34) + m(6) * (t230 * t65 + t232 * t66 - t297 * t53) + m(5) * (t230 * t87 + t232 * t88 - t297 * t73); m(7) * (t230 * t46 + t232 * t45 - t297 * t35) + m(6) * (t230 * t70 + t232 * t69 - t297 * t54) + m(5) * (t230 * t94 + t232 * t93 - t297 * t78); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t315 * t256 ^ 2 + t230 ^ 2 + t232 ^ 2); m(6) * t81 + m(7) * t50; t33 * t312 + t22 * t314 + t23 * t313 + t31 * t309 + (t246 * t18 / 0.2e1 - t248 * t17 / 0.2e1) * t247 + m(7) * (t34 * t50 + t43 * t56 + t44 * t55) + m(6) * (t100 * t65 + t53 * t81 + t66 * t99) + t262; t31 * t271 + t18 * t310 + t20 * t313 + t32 * t312 + t19 * t314 + t17 * t311 + m(7) * (t35 * t50 + t45 * t55 + t46 * t56) + m(6) * (t100 * t70 + t54 * t81 + t69 * t99) + t263; m(6) * (t100 * t230 + t232 * t99 - t297 * t81) + m(7) * (t230 * t56 + t232 * t55 - t297 * t50); t214 * t17 + t216 * t18 + t227 * t31 + m(7) * (t50 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(6) * (t100 ^ 2 + t81 ^ 2 + t99 ^ 2) + t279; m(7) * t80; m(7) * (t34 * t80 + t43 * t98 + t44 * t97) + t262; m(7) * (t35 * t80 + t45 * t97 + t46 * t98) + t263; m(7) * (t230 * t98 + t232 * t97 - t297 * t80); m(7) * (t50 * t80 + t55 * t97 + t56 * t98) + t279; m(7) * (t80 ^ 2 + t97 ^ 2 + t98 ^ 2) + t279;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
