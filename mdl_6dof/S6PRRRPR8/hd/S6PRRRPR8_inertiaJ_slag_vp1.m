% Calculate joint inertia matrix for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:18
% EndTime: 2019-03-08 23:48:34
% DurationCPUTime: 7.26s
% Computational Cost: add. (44041->629), mult. (122182->892), div. (0->0), fcn. (160649->14), ass. (0->280)
t303 = m(6) + m(7);
t236 = sin(pkin(12));
t238 = cos(pkin(12));
t243 = sin(qJ(2));
t239 = cos(pkin(6));
t245 = cos(qJ(2));
t288 = t239 * t245;
t229 = -t236 * t243 + t238 * t288;
t289 = t239 * t243;
t230 = t236 * t245 + t238 * t289;
t242 = sin(qJ(3));
t237 = sin(pkin(6));
t293 = sin(pkin(7));
t266 = t237 * t293;
t294 = cos(pkin(7));
t296 = cos(qJ(3));
t198 = t230 * t296 + (t294 * t229 - t238 * t266) * t242;
t267 = t237 * t294;
t218 = -t229 * t293 - t238 * t267;
t241 = sin(qJ(4));
t295 = cos(qJ(4));
t184 = t198 * t295 + t218 * t241;
t231 = -t236 * t288 - t238 * t243;
t232 = -t236 * t289 + t238 * t245;
t200 = t232 * t296 + (t294 * t231 + t236 * t266) * t242;
t219 = -t231 * t293 + t236 * t267;
t186 = t200 * t295 + t219 * t241;
t217 = t239 * t293 * t242 + (t294 * t242 * t245 + t296 * t243) * t237;
t228 = t239 * t294 - t245 * t266;
t202 = t217 * t295 + t228 * t241;
t183 = t198 * t241 - t218 * t295;
t250 = t296 * t293;
t248 = t237 * t250;
t251 = t294 * t296;
t197 = -t229 * t251 + t230 * t242 + t238 * t248;
t240 = sin(qJ(6));
t244 = cos(qJ(6));
t148 = t183 * t244 - t197 * t240;
t149 = t183 * t240 + t197 * t244;
t101 = Icges(7,5) * t149 + Icges(7,6) * t148 + Icges(7,3) * t184;
t103 = Icges(7,4) * t149 + Icges(7,2) * t148 + Icges(7,6) * t184;
t105 = Icges(7,1) * t149 + Icges(7,4) * t148 + Icges(7,5) * t184;
t36 = t101 * t184 + t103 * t148 + t105 * t149;
t185 = t200 * t241 - t219 * t295;
t199 = -t231 * t251 + t232 * t242 - t236 * t248;
t150 = t185 * t244 - t199 * t240;
t151 = t185 * t240 + t199 * t244;
t102 = Icges(7,5) * t151 + Icges(7,6) * t150 + Icges(7,3) * t186;
t104 = Icges(7,4) * t151 + Icges(7,2) * t150 + Icges(7,6) * t186;
t106 = Icges(7,1) * t151 + Icges(7,4) * t150 + Icges(7,5) * t186;
t37 = t102 * t184 + t104 * t148 + t106 * t149;
t201 = t217 * t241 - t228 * t295;
t290 = t237 * t243;
t216 = -t237 * t245 * t251 - t239 * t250 + t242 * t290;
t181 = t201 * t244 - t216 * t240;
t182 = t201 * t240 + t216 * t244;
t117 = Icges(7,5) * t182 + Icges(7,6) * t181 + Icges(7,3) * t202;
t122 = Icges(7,4) * t182 + Icges(7,2) * t181 + Icges(7,6) * t202;
t127 = Icges(7,1) * t182 + Icges(7,4) * t181 + Icges(7,5) * t202;
t54 = t117 * t184 + t122 * t148 + t127 * t149;
t1 = t184 * t36 + t186 * t37 + t202 * t54;
t302 = t1 / 0.2e1;
t38 = t101 * t186 + t103 * t150 + t105 * t151;
t39 = t102 * t186 + t104 * t150 + t106 * t151;
t55 = t117 * t186 + t122 * t150 + t127 * t151;
t2 = t184 * t38 + t186 * t39 + t202 * t55;
t301 = t2 / 0.2e1;
t46 = t101 * t202 + t103 * t181 + t105 * t182;
t47 = t102 * t202 + t104 * t181 + t106 * t182;
t66 = t117 * t202 + t122 * t181 + t127 * t182;
t9 = t184 * t46 + t186 * t47 + t202 * t66;
t300 = t9 / 0.2e1;
t299 = t184 / 0.2e1;
t298 = t186 / 0.2e1;
t297 = t202 / 0.2e1;
t292 = t236 * t237;
t291 = t237 * t238;
t107 = rSges(7,1) * t149 + rSges(7,2) * t148 + rSges(7,3) * t184;
t287 = pkin(5) * t197 + pkin(11) * t184 + t107;
t108 = rSges(7,1) * t151 + rSges(7,2) * t150 + rSges(7,3) * t186;
t286 = pkin(5) * t199 + pkin(11) * t186 + t108;
t132 = rSges(6,1) * t197 - rSges(6,2) * t184 + rSges(6,3) * t183;
t145 = pkin(4) * t184 + qJ(5) * t183;
t285 = -t132 - t145;
t133 = rSges(6,1) * t199 - rSges(6,2) * t186 + rSges(6,3) * t185;
t146 = pkin(4) * t186 + qJ(5) * t185;
t284 = -t133 - t146;
t134 = rSges(7,1) * t182 + rSges(7,2) * t181 + rSges(7,3) * t202;
t283 = pkin(5) * t216 + pkin(11) * t202 + t134;
t135 = rSges(5,1) * t184 - rSges(5,2) * t183 + rSges(5,3) * t197;
t176 = pkin(3) * t198 + pkin(10) * t197;
t282 = -t135 - t176;
t170 = t219 * t176;
t281 = t219 * t145 + t170;
t177 = pkin(3) * t200 + pkin(10) * t199;
t172 = t228 * t177;
t280 = t228 * t146 + t172;
t168 = rSges(6,1) * t216 - rSges(6,2) * t202 + rSges(6,3) * t201;
t178 = pkin(4) * t202 + qJ(5) * t201;
t279 = -t168 - t178;
t169 = rSges(5,1) * t202 - rSges(5,2) * t201 + rSges(5,3) * t216;
t194 = pkin(3) * t217 + pkin(10) * t216;
t278 = -t169 - t194;
t180 = t218 * t194;
t277 = t218 * t178 + t180;
t205 = t232 * pkin(2) + pkin(9) * t219;
t203 = t239 * t205;
t276 = t239 * t177 + t203;
t204 = t230 * pkin(2) + pkin(9) * t218;
t275 = t204 * t292 + t205 * t291;
t273 = -t145 - t287;
t272 = -t146 - t286;
t271 = -t176 + t285;
t270 = -t178 - t283;
t269 = t239 * t146 + t276;
t268 = -t194 + t279;
t191 = rSges(4,1) * t217 - rSges(4,2) * t216 + rSges(4,3) * t228;
t220 = pkin(2) * t290 + pkin(9) * t228;
t265 = (-t191 - t220) * t237;
t264 = -t176 + t273;
t263 = -t194 + t270;
t262 = t176 * t292 + t177 * t291 + t275;
t118 = Icges(6,5) * t197 - Icges(6,6) * t184 + Icges(6,3) * t183;
t123 = Icges(6,4) * t197 - Icges(6,2) * t184 + Icges(6,6) * t183;
t128 = Icges(6,1) * t197 - Icges(6,4) * t184 + Icges(6,5) * t183;
t62 = t118 * t183 - t123 * t184 + t128 * t197;
t119 = Icges(6,5) * t199 - Icges(6,6) * t186 + Icges(6,3) * t185;
t124 = Icges(6,4) * t199 - Icges(6,2) * t186 + Icges(6,6) * t185;
t129 = Icges(6,1) * t199 - Icges(6,4) * t186 + Icges(6,5) * t185;
t63 = t119 * t183 - t124 * t184 + t129 * t197;
t160 = Icges(6,5) * t216 - Icges(6,6) * t202 + Icges(6,3) * t201;
t162 = Icges(6,4) * t216 - Icges(6,2) * t202 + Icges(6,6) * t201;
t164 = Icges(6,1) * t216 - Icges(6,4) * t202 + Icges(6,5) * t201;
t81 = t160 * t183 - t162 * t184 + t164 * t197;
t13 = t197 * t62 + t199 * t63 + t216 * t81;
t120 = Icges(5,5) * t184 - Icges(5,6) * t183 + Icges(5,3) * t197;
t125 = Icges(5,4) * t184 - Icges(5,2) * t183 + Icges(5,6) * t197;
t130 = Icges(5,1) * t184 - Icges(5,4) * t183 + Icges(5,5) * t197;
t67 = t120 * t197 - t125 * t183 + t130 * t184;
t121 = Icges(5,5) * t186 - Icges(5,6) * t185 + Icges(5,3) * t199;
t126 = Icges(5,4) * t186 - Icges(5,2) * t185 + Icges(5,6) * t199;
t131 = Icges(5,1) * t186 - Icges(5,4) * t185 + Icges(5,5) * t199;
t68 = t121 * t197 - t126 * t183 + t131 * t184;
t161 = Icges(5,5) * t202 - Icges(5,6) * t201 + Icges(5,3) * t216;
t163 = Icges(5,4) * t202 - Icges(5,2) * t201 + Icges(5,6) * t216;
t165 = Icges(5,1) * t202 - Icges(5,4) * t201 + Icges(5,5) * t216;
t83 = t161 * t197 - t163 * t183 + t165 * t184;
t15 = t197 * t67 + t199 * t68 + t216 * t83;
t3 = t197 * t36 + t199 * t37 + t216 * t54;
t261 = t3 / 0.2e1 + t15 / 0.2e1 + t13 / 0.2e1;
t64 = t118 * t185 - t123 * t186 + t128 * t199;
t65 = t119 * t185 - t124 * t186 + t129 * t199;
t82 = t160 * t185 - t162 * t186 + t164 * t199;
t14 = t197 * t64 + t199 * t65 + t216 * t82;
t69 = t120 * t199 - t125 * t185 + t130 * t186;
t70 = t121 * t199 - t126 * t185 + t131 * t186;
t84 = t161 * t199 - t163 * t185 + t165 * t186;
t16 = t197 * t69 + t199 * t70 + t216 * t84;
t4 = t197 * t38 + t199 * t39 + t216 * t55;
t260 = t4 / 0.2e1 + t14 / 0.2e1 + t16 / 0.2e1;
t17 = t218 * t62 + t219 * t63 + t228 * t81;
t19 = t218 * t67 + t219 * t68 + t228 * t83;
t5 = t218 * t36 + t219 * t37 + t228 * t54;
t259 = t5 / 0.2e1 + t17 / 0.2e1 + t19 / 0.2e1;
t18 = t218 * t64 + t219 * t65 + t228 * t82;
t20 = t218 * t69 + t219 * t70 + t228 * t84;
t6 = t218 * t38 + t219 * t39 + t228 * t55;
t258 = t6 / 0.2e1 + t18 / 0.2e1 + t20 / 0.2e1;
t21 = t239 * t81 + (t236 * t63 - t238 * t62) * t237;
t23 = t239 * t83 + (t236 * t68 - t238 * t67) * t237;
t7 = t239 * t54 + (t236 * t37 - t238 * t36) * t237;
t257 = t7 / 0.2e1 + t21 / 0.2e1 + t23 / 0.2e1;
t22 = t239 * t82 + (t236 * t65 - t238 * t64) * t237;
t24 = t239 * t84 + (t236 * t70 - t238 * t69) * t237;
t8 = t239 * t55 + (t236 * t39 - t238 * t38) * t237;
t256 = t8 / 0.2e1 + t22 / 0.2e1 + t24 / 0.2e1;
t10 = t197 * t46 + t199 * t47 + t216 * t66;
t71 = t118 * t201 - t123 * t202 + t128 * t216;
t72 = t119 * t201 - t124 * t202 + t129 * t216;
t94 = t160 * t201 - t162 * t202 + t164 * t216;
t25 = t197 * t71 + t199 * t72 + t216 * t94;
t73 = t120 * t216 - t125 * t201 + t130 * t202;
t74 = t121 * t216 - t126 * t201 + t131 * t202;
t95 = t161 * t216 - t163 * t201 + t165 * t202;
t26 = t197 * t73 + t199 * t74 + t216 * t95;
t255 = t10 / 0.2e1 + t26 / 0.2e1 + t25 / 0.2e1;
t11 = t218 * t46 + t219 * t47 + t228 * t66;
t27 = t218 * t71 + t219 * t72 + t228 * t94;
t28 = t218 * t73 + t219 * t74 + t228 * t95;
t254 = t11 / 0.2e1 + t28 / 0.2e1 + t27 / 0.2e1;
t12 = t239 * t66 + (t236 * t47 - t238 * t46) * t237;
t29 = t239 * t94 + (t236 * t72 - t238 * t71) * t237;
t30 = t239 * t95 + (t236 * t74 - t238 * t73) * t237;
t253 = t12 / 0.2e1 + t29 / 0.2e1 + t30 / 0.2e1;
t252 = (-t220 + t278) * t237;
t249 = (-t220 + t268) * t237;
t247 = t145 * t292 + t146 * t291 + t262;
t246 = (-t220 + t263) * t237;
t226 = t239 * rSges(3,3) + (rSges(3,1) * t243 + rSges(3,2) * t245) * t237;
t225 = Icges(3,5) * t239 + (Icges(3,1) * t243 + Icges(3,4) * t245) * t237;
t224 = Icges(3,6) * t239 + (Icges(3,4) * t243 + Icges(3,2) * t245) * t237;
t223 = Icges(3,3) * t239 + (Icges(3,5) * t243 + Icges(3,6) * t245) * t237;
t213 = rSges(3,1) * t232 + rSges(3,2) * t231 + rSges(3,3) * t292;
t212 = rSges(3,1) * t230 + rSges(3,2) * t229 - rSges(3,3) * t291;
t211 = Icges(3,1) * t232 + Icges(3,4) * t231 + Icges(3,5) * t292;
t210 = Icges(3,1) * t230 + Icges(3,4) * t229 - Icges(3,5) * t291;
t209 = Icges(3,4) * t232 + Icges(3,2) * t231 + Icges(3,6) * t292;
t208 = Icges(3,4) * t230 + Icges(3,2) * t229 - Icges(3,6) * t291;
t207 = Icges(3,5) * t232 + Icges(3,6) * t231 + Icges(3,3) * t292;
t206 = Icges(3,5) * t230 + Icges(3,6) * t229 - Icges(3,3) * t291;
t193 = -t212 * t239 - t226 * t291;
t192 = t213 * t239 - t226 * t292;
t190 = Icges(4,1) * t217 - Icges(4,4) * t216 + Icges(4,5) * t228;
t189 = Icges(4,4) * t217 - Icges(4,2) * t216 + Icges(4,6) * t228;
t188 = Icges(4,5) * t217 - Icges(4,6) * t216 + Icges(4,3) * t228;
t179 = (t212 * t236 + t213 * t238) * t237;
t167 = rSges(4,1) * t200 - rSges(4,2) * t199 + rSges(4,3) * t219;
t166 = rSges(4,1) * t198 - rSges(4,2) * t197 + rSges(4,3) * t218;
t159 = Icges(4,1) * t200 - Icges(4,4) * t199 + Icges(4,5) * t219;
t158 = Icges(4,1) * t198 - Icges(4,4) * t197 + Icges(4,5) * t218;
t157 = Icges(4,4) * t200 - Icges(4,2) * t199 + Icges(4,6) * t219;
t156 = Icges(4,4) * t198 - Icges(4,2) * t197 + Icges(4,6) * t218;
t155 = Icges(4,5) * t200 - Icges(4,6) * t199 + Icges(4,3) * t219;
t154 = Icges(4,5) * t198 - Icges(4,6) * t197 + Icges(4,3) * t218;
t147 = t197 * t178;
t139 = t216 * t146;
t137 = t199 * t145;
t136 = rSges(5,1) * t186 - rSges(5,2) * t185 + rSges(5,3) * t199;
t116 = t167 * t228 - t191 * t219;
t115 = -t166 * t228 + t191 * t218;
t114 = (-t166 - t204) * t239 + t238 * t265;
t113 = t167 * t239 + t236 * t265 + t203;
t112 = t188 * t228 - t189 * t216 + t190 * t217;
t111 = t166 * t219 - t167 * t218;
t110 = t188 * t219 - t189 * t199 + t190 * t200;
t109 = t188 * t218 - t189 * t197 + t190 * t198;
t100 = (t166 * t236 + t167 * t238) * t237 + t275;
t99 = t136 * t216 - t169 * t199;
t98 = -t135 * t216 + t169 * t197;
t97 = t155 * t228 - t157 * t216 + t159 * t217;
t96 = t154 * t228 - t156 * t216 + t158 * t217;
t93 = t155 * t219 - t157 * t199 + t159 * t200;
t92 = t154 * t219 - t156 * t199 + t158 * t200;
t91 = t155 * t218 - t157 * t197 + t159 * t198;
t90 = t154 * t218 - t156 * t197 + t158 * t198;
t89 = t135 * t199 - t136 * t197;
t88 = t136 * t228 + t278 * t219 + t172;
t87 = t169 * t218 + t282 * t228 + t180;
t86 = (-t204 + t282) * t239 + t238 * t252;
t85 = t136 * t239 + t236 * t252 + t276;
t80 = t108 * t202 - t134 * t186;
t79 = -t107 * t202 + t134 * t184;
t78 = t135 * t219 + t170 + (-t136 - t177) * t218;
t77 = (t135 * t236 + t136 * t238) * t237 + t262;
t76 = t133 * t216 + t279 * t199 + t139;
t75 = t168 * t197 + t285 * t216 + t147;
t61 = (-t204 + t271) * t239 + t238 * t249;
t60 = t133 * t239 + t236 * t249 + t269;
t59 = t107 * t186 - t108 * t184;
t58 = t133 * t228 + t268 * t219 + t280;
t57 = t168 * t218 + t271 * t228 + t277;
t56 = t132 * t199 + t284 * t197 + t137;
t53 = (t132 * t236 + t133 * t238) * t237 + t247;
t52 = t132 * t219 + (-t177 + t284) * t218 + t281;
t51 = t112 * t239 + (t236 * t97 - t238 * t96) * t237;
t50 = t112 * t228 + t218 * t96 + t219 * t97;
t49 = t270 * t199 + t286 * t216 + t139;
t48 = t283 * t197 + t273 * t216 + t147;
t45 = (-t204 + t264) * t239 + t238 * t246;
t44 = t236 * t246 + t286 * t239 + t269;
t43 = t263 * t219 + t286 * t228 + t280;
t42 = t283 * t218 + t264 * t228 + t277;
t41 = t110 * t239 + (t236 * t93 - t238 * t92) * t237;
t40 = t109 * t239 + (t236 * t91 - t238 * t90) * t237;
t35 = t110 * t228 + t218 * t92 + t219 * t93;
t34 = t109 * t228 + t218 * t90 + t219 * t91;
t33 = t272 * t197 + t287 * t199 + t137;
t32 = (t287 * t236 + t286 * t238) * t237 + t247;
t31 = t287 * t219 + (-t177 + t272) * t218 + t281;
t138 = [m(4) + m(5) + m(2) + m(3) + t303; m(3) * t179 + m(4) * t100 + m(5) * t77 + m(6) * t53 + m(7) * t32; m(7) * (t32 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t53 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t77 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(4) * (t100 ^ 2 + t113 ^ 2 + t114 ^ 2) + m(3) * (t179 ^ 2 + t192 ^ 2 + t193 ^ 2) + (t8 + t22 + t24 + t41 + (t207 * t292 + t209 * t231 + t211 * t232) * t292) * t292 + (-t7 - t21 - t23 - t40 + (-t206 * t291 + t208 * t229 + t210 * t230) * t291 + (-t206 * t292 + t207 * t291 - t208 * t231 - t209 * t229 - t210 * t232 - t211 * t230) * t292) * t291 + ((t223 * t292 + t231 * t224 + t232 * t225) * t292 - (-t223 * t291 + t229 * t224 + t230 * t225) * t291 + t12 + t30 + t29 + t51 + ((t209 * t245 + t211 * t243) * t236 - (t208 * t245 + t210 * t243) * t238) * t237 ^ 2 + ((-t206 * t238 + t207 * t236 + t224 * t245 + t225 * t243) * t237 + t239 * t223) * t239) * t239; m(4) * t111 + m(5) * t78 + m(6) * t52 + m(7) * t31; (t50 / 0.2e1 + t254) * t239 + (t51 / 0.2e1 + t253) * t228 + (t41 / 0.2e1 + t256) * t219 + (t40 / 0.2e1 + t257) * t218 + m(7) * (t31 * t32 + t42 * t45 + t43 * t44) + m(6) * (t52 * t53 + t57 * t61 + t58 * t60) + m(5) * (t77 * t78 + t85 * t88 + t86 * t87) + m(4) * (t100 * t111 + t113 * t116 + t114 * t115) + ((-t34 / 0.2e1 - t259) * t238 + (t35 / 0.2e1 + t258) * t236) * t237; (t11 + t28 + t27 + t50) * t228 + (t6 + t18 + t20 + t35) * t219 + (t5 + t19 + t17 + t34) * t218 + m(7) * (t31 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t52 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(5) * (t78 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t111 ^ 2 + t115 ^ 2 + t116 ^ 2); m(5) * t89 + m(6) * t56 + m(7) * t33; t255 * t239 + t253 * t216 + t256 * t199 + t257 * t197 + m(7) * (t32 * t33 + t44 * t49 + t45 * t48) + m(6) * (t53 * t56 + t60 * t76 + t61 * t75) + m(5) * (t77 * t89 + t85 * t99 + t86 * t98) + (t260 * t236 - t261 * t238) * t237; t255 * t228 + t260 * t219 + t261 * t218 + t254 * t216 + t258 * t199 + t259 * t197 + m(7) * (t31 * t33 + t42 * t48 + t43 * t49) + m(6) * (t52 * t56 + t57 * t75 + t58 * t76) + m(5) * (t78 * t89 + t87 * t98 + t88 * t99); (t10 + t25 + t26) * t216 + (t4 + t14 + t16) * t199 + (t3 + t15 + t13) * t197 + m(7) * (t33 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t56 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(5) * (t89 ^ 2 + t98 ^ 2 + t99 ^ 2); t201 * t303; m(7) * (t183 * t44 + t185 * t45 + t201 * t32) + m(6) * (t183 * t60 + t185 * t61 + t201 * t53); m(7) * (t183 * t43 + t185 * t42 + t201 * t31) + m(6) * (t183 * t58 + t185 * t57 + t201 * t52); m(7) * (t183 * t49 + t185 * t48 + t201 * t33) + m(6) * (t183 * t76 + t185 * t75 + t201 * t56); (t183 ^ 2 + t185 ^ 2 + t201 ^ 2) * t303; m(7) * t59; m(7) * (t32 * t59 + t44 * t80 + t45 * t79) + t7 * t299 + t8 * t298 + t239 * t300 + t12 * t297 + (t236 * t301 - t238 * t1 / 0.2e1) * t237; t218 * t302 + t6 * t298 + t11 * t297 + t5 * t299 + m(7) * (t31 * t59 + t42 * t79 + t43 * t80) + t228 * t300 + t219 * t301; t199 * t301 + t4 * t298 + t216 * t300 + t3 * t299 + t10 * t297 + t197 * t302 + m(7) * (t33 * t59 + t48 * t79 + t49 * t80); m(7) * (t183 * t80 + t185 * t79 + t201 * t59); t186 * t2 + t184 * t1 + t202 * t9 + m(7) * (t59 ^ 2 + t79 ^ 2 + t80 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t138(1) t138(2) t138(4) t138(7) t138(11) t138(16); t138(2) t138(3) t138(5) t138(8) t138(12) t138(17); t138(4) t138(5) t138(6) t138(9) t138(13) t138(18); t138(7) t138(8) t138(9) t138(10) t138(14) t138(19); t138(11) t138(12) t138(13) t138(14) t138(15) t138(20); t138(16) t138(17) t138(18) t138(19) t138(20) t138(21);];
Mq  = res;
