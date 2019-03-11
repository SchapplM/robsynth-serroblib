% Calculate joint inertia matrix for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:17
% EndTime: 2019-03-08 22:30:27
% DurationCPUTime: 4.82s
% Computational Cost: add. (20274->580), mult. (44899->832), div. (0->0), fcn. (57223->12), ass. (0->262)
t298 = m(5) + m(6) + m(7);
t239 = sin(pkin(11));
t241 = cos(pkin(11));
t246 = cos(qJ(2));
t242 = cos(pkin(6));
t244 = sin(qJ(2));
t281 = t242 * t244;
t225 = t239 * t246 + t241 * t281;
t240 = sin(pkin(6));
t290 = sin(qJ(3));
t261 = t240 * t290;
t291 = cos(qJ(3));
t214 = t225 * t291 - t241 * t261;
t297 = t214 / 0.2e1;
t227 = -t239 * t281 + t241 * t246;
t216 = t227 * t291 + t239 * t261;
t296 = t216 / 0.2e1;
t280 = t242 * t246;
t224 = t239 * t244 - t241 * t280;
t295 = t224 / 0.2e1;
t226 = t239 * t280 + t241 * t244;
t294 = t226 / 0.2e1;
t262 = t240 * t291;
t229 = t242 * t290 + t244 * t262;
t293 = t229 / 0.2e1;
t292 = t242 / 0.2e1;
t245 = cos(qJ(5));
t289 = pkin(5) * t245;
t213 = t225 * t290 + t241 * t262;
t243 = sin(qJ(5));
t287 = t213 * t243;
t215 = t227 * t290 - t239 * t262;
t286 = t215 * t243;
t228 = -t242 * t291 + t244 * t261;
t285 = t228 * t243;
t284 = t239 * t240;
t283 = t240 * t241;
t282 = t240 * t246;
t238 = qJ(5) + qJ(6);
t236 = sin(t238);
t237 = cos(t238);
t173 = t213 * t237 - t224 * t236;
t174 = t213 * t236 + t224 * t237;
t121 = rSges(7,1) * t174 + rSges(7,2) * t173 + rSges(7,3) * t214;
t133 = pkin(5) * t287 + pkin(10) * t214 + t224 * t289;
t279 = t121 + t133;
t175 = t215 * t237 - t226 * t236;
t176 = t215 * t236 + t226 * t237;
t122 = rSges(7,1) * t176 + rSges(7,2) * t175 + rSges(7,3) * t216;
t134 = pkin(5) * t286 + pkin(10) * t216 + t226 * t289;
t278 = t122 + t134;
t211 = t228 * t237 + t236 * t282;
t212 = t228 * t236 - t237 * t282;
t140 = rSges(7,1) * t212 + rSges(7,2) * t211 + rSges(7,3) * t229;
t163 = pkin(5) * t285 + pkin(10) * t229 - t282 * t289;
t277 = t140 + t163;
t154 = rSges(5,1) * t226 - rSges(5,2) * t216 + rSges(5,3) * t215;
t172 = pkin(3) * t216 + qJ(4) * t215;
t276 = -t154 - t172;
t171 = pkin(3) * t214 + qJ(4) * t213;
t157 = t226 * t171;
t186 = t224 * pkin(4) + t214 * pkin(9);
t275 = t226 * t186 + t157;
t210 = pkin(3) * t229 + qJ(4) * t228;
t274 = t171 * t282 + t224 * t210;
t209 = pkin(2) * t227 + pkin(8) * t226;
t207 = t242 * t209;
t273 = t242 * t172 + t207;
t208 = pkin(2) * t225 + pkin(8) * t224;
t272 = -t171 - t208;
t187 = t226 * pkin(4) + t216 * pkin(9);
t271 = -t172 - t187;
t202 = -rSges(5,1) * t282 - t229 * rSges(5,2) + t228 * rSges(5,3);
t270 = -t202 - t210;
t269 = t208 * t284 + t209 * t283;
t219 = -pkin(4) * t282 + t229 * pkin(9);
t268 = -t210 - t219;
t115 = Icges(7,5) * t174 + Icges(7,6) * t173 + Icges(7,3) * t214;
t117 = Icges(7,4) * t174 + Icges(7,2) * t173 + Icges(7,6) * t214;
t119 = Icges(7,1) * t174 + Icges(7,4) * t173 + Icges(7,5) * t214;
t65 = t115 * t229 + t117 * t211 + t119 * t212;
t116 = Icges(7,5) * t176 + Icges(7,6) * t175 + Icges(7,3) * t216;
t118 = Icges(7,4) * t176 + Icges(7,2) * t175 + Icges(7,6) * t216;
t120 = Icges(7,1) * t176 + Icges(7,4) * t175 + Icges(7,5) * t216;
t66 = t116 * t229 + t118 * t211 + t120 * t212;
t137 = Icges(7,5) * t212 + Icges(7,6) * t211 + Icges(7,3) * t229;
t138 = Icges(7,4) * t212 + Icges(7,2) * t211 + Icges(7,6) * t229;
t139 = Icges(7,1) * t212 + Icges(7,4) * t211 + Icges(7,5) * t229;
t80 = t137 * t229 + t138 * t211 + t139 * t212;
t26 = t214 * t65 + t216 * t66 + t229 * t80;
t55 = t115 * t214 + t117 * t173 + t119 * t174;
t56 = t116 * t214 + t118 * t173 + t120 * t174;
t73 = t137 * t214 + t138 * t173 + t139 * t174;
t7 = t214 * t55 + t216 * t56 + t229 * t73;
t57 = t115 * t216 + t117 * t175 + t119 * t176;
t58 = t116 * t216 + t118 * t175 + t120 * t176;
t74 = t137 * t216 + t138 * t175 + t139 * t176;
t8 = t214 * t57 + t216 * t58 + t229 * t74;
t267 = t214 * t7 + t216 * t8 + t229 * t26;
t183 = t215 * t245 - t226 * t243;
t184 = t226 * t245 + t286;
t132 = rSges(6,1) * t184 + rSges(6,2) * t183 + rSges(6,3) * t216;
t266 = -t132 + t271;
t217 = t228 * t245 + t243 * t282;
t218 = -t245 * t282 + t285;
t161 = rSges(6,1) * t218 + rSges(6,2) * t217 + rSges(6,3) * t229;
t265 = -t161 + t268;
t264 = t242 * t187 + t273;
t263 = -t186 + t272;
t260 = -t282 / 0.2e1;
t203 = t229 * rSges(4,1) - t228 * rSges(4,2) - rSges(4,3) * t282;
t230 = (pkin(2) * t244 - pkin(8) * t246) * t240;
t259 = (-t203 - t230) * t240;
t257 = t271 - t278;
t256 = t268 - t277;
t255 = t171 * t284 + t172 * t283 + t269;
t254 = t186 * t282 + t224 * t219 + t274;
t253 = (-t230 + t270) * t240;
t13 = t55 * t224 + t56 * t226 - t282 * t73;
t14 = t57 * t224 + t58 * t226 - t282 * t74;
t28 = t65 * t224 + t66 * t226 - t282 * t80;
t252 = t13 * t297 + t14 * t296 + t26 * t260 + t28 * t293 + t8 * t294 + t7 * t295;
t15 = t73 * t242 + (t239 * t56 - t241 * t55) * t240;
t16 = t74 * t242 + (t239 * t58 - t241 * t57) * t240;
t30 = t80 * t242 + (t239 * t66 - t241 * t65) * t240;
t251 = t15 * t297 + t16 * t296 + t26 * t292 + t30 * t293 + t8 * t284 / 0.2e1 - t7 * t283 / 0.2e1;
t250 = (-t230 + t265) * t240;
t249 = t186 * t284 + t187 * t283 + t255;
t248 = (-t230 + t256) * t240;
t223 = t242 * rSges(3,3) + (rSges(3,1) * t244 + rSges(3,2) * t246) * t240;
t222 = Icges(3,5) * t242 + (Icges(3,1) * t244 + Icges(3,4) * t246) * t240;
t221 = Icges(3,6) * t242 + (Icges(3,4) * t244 + Icges(3,2) * t246) * t240;
t220 = Icges(3,3) * t242 + (Icges(3,5) * t244 + Icges(3,6) * t246) * t240;
t201 = Icges(4,1) * t229 - Icges(4,4) * t228 - Icges(4,5) * t282;
t200 = Icges(4,4) * t229 - Icges(4,2) * t228 - Icges(4,6) * t282;
t199 = Icges(4,5) * t229 - Icges(4,6) * t228 - Icges(4,3) * t282;
t198 = -Icges(5,1) * t282 - Icges(5,4) * t229 + Icges(5,5) * t228;
t197 = -Icges(5,4) * t282 - Icges(5,2) * t229 + Icges(5,6) * t228;
t196 = -Icges(5,5) * t282 - Icges(5,6) * t229 + Icges(5,3) * t228;
t195 = t227 * rSges(3,1) - t226 * rSges(3,2) + rSges(3,3) * t284;
t194 = t225 * rSges(3,1) - t224 * rSges(3,2) - rSges(3,3) * t283;
t193 = Icges(3,1) * t227 - Icges(3,4) * t226 + Icges(3,5) * t284;
t192 = Icges(3,1) * t225 - Icges(3,4) * t224 - Icges(3,5) * t283;
t191 = Icges(3,4) * t227 - Icges(3,2) * t226 + Icges(3,6) * t284;
t190 = Icges(3,4) * t225 - Icges(3,2) * t224 - Icges(3,6) * t283;
t189 = Icges(3,5) * t227 - Icges(3,6) * t226 + Icges(3,3) * t284;
t188 = Icges(3,5) * t225 - Icges(3,6) * t224 - Icges(3,3) * t283;
t182 = t224 * t245 + t287;
t181 = t213 * t245 - t224 * t243;
t165 = -t194 * t242 - t223 * t283;
t164 = t195 * t242 - t223 * t284;
t160 = Icges(6,1) * t218 + Icges(6,4) * t217 + Icges(6,5) * t229;
t159 = Icges(6,4) * t218 + Icges(6,2) * t217 + Icges(6,6) * t229;
t158 = Icges(6,5) * t218 + Icges(6,6) * t217 + Icges(6,3) * t229;
t156 = rSges(4,1) * t216 - rSges(4,2) * t215 + rSges(4,3) * t226;
t155 = rSges(4,1) * t214 - rSges(4,2) * t213 + rSges(4,3) * t224;
t153 = rSges(5,1) * t224 - rSges(5,2) * t214 + rSges(5,3) * t213;
t152 = Icges(4,1) * t216 - Icges(4,4) * t215 + Icges(4,5) * t226;
t151 = Icges(4,1) * t214 - Icges(4,4) * t213 + Icges(4,5) * t224;
t150 = Icges(5,1) * t226 - Icges(5,4) * t216 + Icges(5,5) * t215;
t149 = Icges(5,1) * t224 - Icges(5,4) * t214 + Icges(5,5) * t213;
t148 = Icges(4,4) * t216 - Icges(4,2) * t215 + Icges(4,6) * t226;
t147 = Icges(4,4) * t214 - Icges(4,2) * t213 + Icges(4,6) * t224;
t146 = Icges(5,4) * t226 - Icges(5,2) * t216 + Icges(5,6) * t215;
t145 = Icges(5,4) * t224 - Icges(5,2) * t214 + Icges(5,6) * t213;
t144 = Icges(4,5) * t216 - Icges(4,6) * t215 + Icges(4,3) * t226;
t143 = Icges(4,5) * t214 - Icges(4,6) * t213 + Icges(4,3) * t224;
t142 = Icges(5,5) * t226 - Icges(5,6) * t216 + Icges(5,3) * t215;
t141 = Icges(5,5) * t224 - Icges(5,6) * t214 + Icges(5,3) * t213;
t136 = (t194 * t239 + t195 * t241) * t240;
t135 = t214 * t140;
t131 = rSges(6,1) * t182 + rSges(6,2) * t181 + rSges(6,3) * t214;
t130 = Icges(6,1) * t184 + Icges(6,4) * t183 + Icges(6,5) * t216;
t129 = Icges(6,1) * t182 + Icges(6,4) * t181 + Icges(6,5) * t214;
t128 = Icges(6,4) * t184 + Icges(6,2) * t183 + Icges(6,6) * t216;
t127 = Icges(6,4) * t182 + Icges(6,2) * t181 + Icges(6,6) * t214;
t126 = Icges(6,5) * t184 + Icges(6,6) * t183 + Icges(6,3) * t216;
t125 = Icges(6,5) * t182 + Icges(6,6) * t181 + Icges(6,3) * t214;
t124 = -t156 * t282 - t226 * t203;
t123 = t155 * t282 + t224 * t203;
t114 = t229 * t122;
t113 = -t199 * t282 - t228 * t200 + t229 * t201;
t112 = t228 * t196 - t229 * t197 - t198 * t282;
t111 = t216 * t121;
t110 = t155 * t226 - t156 * t224;
t109 = (-t155 - t208) * t242 + t241 * t259;
t108 = t242 * t156 + t239 * t259 + t207;
t107 = t199 * t226 - t200 * t215 + t201 * t216;
t106 = t199 * t224 - t200 * t213 + t201 * t214;
t105 = t196 * t215 - t197 * t216 + t198 * t226;
t104 = t196 * t213 - t197 * t214 + t198 * t224;
t103 = (t155 * t239 + t156 * t241) * t240 + t269;
t102 = t132 * t229 - t161 * t216;
t101 = -t131 * t229 + t161 * t214;
t100 = t226 * t270 + t276 * t282;
t99 = t153 * t282 + t224 * t202 + t274;
t98 = -t140 * t216 + t114;
t97 = -t121 * t229 + t135;
t96 = -t144 * t282 - t228 * t148 + t229 * t152;
t95 = -t143 * t282 - t228 * t147 + t229 * t151;
t94 = t228 * t142 - t229 * t146 - t150 * t282;
t93 = t228 * t141 - t229 * t145 - t149 * t282;
t92 = (-t153 + t272) * t242 + t241 * t253;
t91 = t242 * t154 + t239 * t253 + t273;
t90 = t158 * t229 + t159 * t217 + t160 * t218;
t89 = t144 * t226 - t148 * t215 + t152 * t216;
t88 = t143 * t226 - t147 * t215 + t151 * t216;
t87 = t144 * t224 - t148 * t213 + t152 * t214;
t86 = t143 * t224 - t147 * t213 + t151 * t214;
t85 = t142 * t215 - t146 * t216 + t150 * t226;
t84 = t141 * t215 - t145 * t216 + t149 * t226;
t83 = t142 * t213 - t146 * t214 + t150 * t224;
t82 = t141 * t213 - t145 * t214 + t149 * t224;
t81 = t131 * t216 - t132 * t214;
t79 = t226 * t153 + t224 * t276 + t157;
t78 = -t122 * t214 + t111;
t77 = t158 * t216 + t159 * t183 + t160 * t184;
t76 = t158 * t214 + t159 * t181 + t160 * t182;
t75 = (t153 * t239 + t154 * t241) * t240 + t255;
t72 = t226 * t265 + t266 * t282;
t71 = t131 * t282 + t224 * t161 + t254;
t70 = (-t131 + t263) * t242 + t241 * t250;
t69 = t242 * t132 + t239 * t250 + t264;
t68 = t126 * t229 + t128 * t217 + t130 * t218;
t67 = t125 * t229 + t127 * t217 + t129 * t218;
t64 = t229 * t134 - t216 * t277 + t114;
t63 = t214 * t163 - t229 * t279 + t135;
t62 = t126 * t216 + t128 * t183 + t130 * t184;
t61 = t125 * t216 + t127 * t183 + t129 * t184;
t60 = t126 * t214 + t128 * t181 + t130 * t182;
t59 = t125 * t214 + t127 * t181 + t129 * t182;
t54 = t226 * t131 + t224 * t266 + t275;
t53 = (t131 * t239 + t132 * t241) * t240 + t249;
t52 = t216 * t133 - t214 * t278 + t111;
t51 = t226 * t256 + t257 * t282;
t50 = t224 * t277 + t279 * t282 + t254;
t49 = (t263 - t279) * t242 + t241 * t248;
t48 = t239 * t248 + t242 * t278 + t264;
t47 = t113 * t242 + (t239 * t96 - t241 * t95) * t240;
t46 = t112 * t242 + (t239 * t94 - t241 * t93) * t240;
t45 = -t113 * t282 + t95 * t224 + t96 * t226;
t44 = -t112 * t282 + t93 * t224 + t94 * t226;
t43 = t107 * t242 + (t239 * t89 - t241 * t88) * t240;
t42 = t106 * t242 + (t239 * t87 - t241 * t86) * t240;
t41 = t105 * t242 + (t239 * t85 - t241 * t84) * t240;
t40 = t104 * t242 + (t239 * t83 - t241 * t82) * t240;
t39 = -t107 * t282 + t88 * t224 + t89 * t226;
t38 = -t106 * t282 + t86 * t224 + t87 * t226;
t37 = -t105 * t282 + t84 * t224 + t85 * t226;
t36 = -t104 * t282 + t82 * t224 + t83 * t226;
t35 = t224 * t257 + t226 * t279 + t275;
t34 = (t239 * t279 + t241 * t278) * t240 + t249;
t33 = t90 * t242 + (t239 * t68 - t241 * t67) * t240;
t32 = t67 * t224 + t68 * t226 - t282 * t90;
t31 = t214 * t67 + t216 * t68 + t229 * t90;
t22 = t77 * t242 + (t239 * t62 - t241 * t61) * t240;
t21 = t76 * t242 + (t239 * t60 - t241 * t59) * t240;
t20 = t61 * t224 + t62 * t226 - t282 * t77;
t19 = t59 * t224 + t60 * t226 - t282 * t76;
t18 = t214 * t61 + t216 * t62 + t229 * t77;
t17 = t214 * t59 + t216 * t60 + t229 * t76;
t1 = [m(2) + m(3) + m(4) + t298; m(3) * t136 + m(4) * t103 + m(5) * t75 + m(6) * t53 + m(7) * t34; m(7) * (t34 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t53 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t75 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(4) * (t103 ^ 2 + t108 ^ 2 + t109 ^ 2) + m(3) * (t136 ^ 2 + t164 ^ 2 + t165 ^ 2) + (t16 + t22 + t43 + t41 + (t189 * t284 - t226 * t191 + t227 * t193) * t284) * t284 + (-t15 - t21 - t40 - t42 + (-t188 * t283 - t224 * t190 + t225 * t192) * t283 + (-t188 * t284 + t189 * t283 + t226 * t190 + t224 * t191 - t227 * t192 - t225 * t193) * t284) * t283 + (-(-t220 * t283 - t224 * t221 + t225 * t222) * t283 + (t220 * t284 - t226 * t221 + t227 * t222) * t284 + t30 + t33 + t47 + t46 + ((t191 * t246 + t193 * t244) * t239 - (t190 * t246 + t192 * t244) * t241) * t240 ^ 2 + ((-t188 * t241 + t189 * t239 + t221 * t246 + t222 * t244) * t240 + t242 * t220) * t242) * t242; m(4) * t110 + m(5) * t79 + m(6) * t54 + m(7) * t35; (t28 / 0.2e1 + t32 / 0.2e1 + t44 / 0.2e1 + t45 / 0.2e1) * t242 + (t16 / 0.2e1 + t22 / 0.2e1 + t41 / 0.2e1 + t43 / 0.2e1) * t226 + (t15 / 0.2e1 + t21 / 0.2e1 + t40 / 0.2e1 + t42 / 0.2e1) * t224 + m(7) * (t34 * t35 + t48 * t51 + t49 * t50) + m(6) * (t53 * t54 + t69 * t72 + t70 * t71) + m(5) * (t100 * t91 + t75 * t79 + t92 * t99) + m(4) * (t103 * t110 + t108 * t124 + t109 * t123) + ((-t30 / 0.2e1 - t33 / 0.2e1 - t46 / 0.2e1 - t47 / 0.2e1) * t246 + (-t13 / 0.2e1 - t19 / 0.2e1 - t36 / 0.2e1 - t38 / 0.2e1) * t241 + (t14 / 0.2e1 + t20 / 0.2e1 + t37 / 0.2e1 + t39 / 0.2e1) * t239) * t240; (-t28 - t32 - t44 - t45) * t282 + (t14 + t20 + t39 + t37) * t226 + (t13 + t19 + t38 + t36) * t224 + m(7) * (t35 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t54 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(5) * (t100 ^ 2 + t79 ^ 2 + t99 ^ 2) + m(4) * (t110 ^ 2 + t123 ^ 2 + t124 ^ 2); t228 * t298; m(7) * (t213 * t48 + t215 * t49 + t228 * t34) + m(6) * (t213 * t69 + t215 * t70 + t228 * t53) + m(5) * (t213 * t91 + t215 * t92 + t228 * t75); m(7) * (t213 * t51 + t215 * t50 + t228 * t35) + m(6) * (t213 * t72 + t215 * t71 + t228 * t54) + m(5) * (t100 * t213 + t215 * t99 + t228 * t79); (t213 ^ 2 + t215 ^ 2 + t228 ^ 2) * t298; m(6) * t81 + m(7) * t52; t21 * t297 + t22 * t296 + t33 * t293 + t31 * t292 + (t239 * t18 / 0.2e1 - t241 * t17 / 0.2e1) * t240 + m(7) * (t34 * t52 + t48 * t64 + t49 * t63) + m(6) * (t101 * t70 + t102 * t69 + t53 * t81) + t251; t31 * t260 + t18 * t294 + t19 * t297 + t20 * t296 + t17 * t295 + t32 * t293 + m(7) * (t35 * t52 + t50 * t63 + t51 * t64) + m(6) * (t101 * t71 + t102 * t72 + t54 * t81) + t252; m(6) * (t101 * t215 + t102 * t213 + t228 * t81) + m(7) * (t213 * t64 + t215 * t63 + t228 * t52); t214 * t17 + t216 * t18 + t229 * t31 + m(7) * (t52 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(6) * (t101 ^ 2 + t102 ^ 2 + t81 ^ 2) + t267; m(7) * t78; m(7) * (t34 * t78 + t48 * t98 + t49 * t97) + t251; m(7) * (t35 * t78 + t50 * t97 + t51 * t98) + t252; m(7) * (t213 * t98 + t215 * t97 + t228 * t78); m(7) * (t52 * t78 + t63 * t97 + t64 * t98) + t267; m(7) * (t78 ^ 2 + t97 ^ 2 + t98 ^ 2) + t267;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
