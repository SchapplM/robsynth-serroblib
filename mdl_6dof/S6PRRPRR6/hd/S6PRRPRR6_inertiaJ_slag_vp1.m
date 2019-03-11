% Calculate joint inertia matrix for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:22:57
% EndTime: 2019-03-08 22:23:13
% DurationCPUTime: 7.03s
% Computational Cost: add. (46816->610), mult. (111544->882), div. (0->0), fcn. (146424->16), ass. (0->281)
t307 = m(5) + m(6) + m(7);
t236 = sin(pkin(12));
t239 = cos(pkin(12));
t244 = sin(qJ(2));
t240 = cos(pkin(6));
t246 = cos(qJ(2));
t289 = t240 * t246;
t226 = -t236 * t244 + t239 * t289;
t290 = t240 * t244;
t227 = t236 * t246 + t239 * t290;
t243 = sin(qJ(3));
t237 = sin(pkin(6));
t297 = sin(pkin(7));
t260 = t237 * t297;
t298 = cos(pkin(7));
t300 = cos(qJ(3));
t197 = t227 * t300 + (t226 * t298 - t239 * t260) * t243;
t261 = t237 * t298;
t215 = -t226 * t297 - t239 * t261;
t274 = pkin(13) + qJ(5);
t234 = sin(t274);
t259 = cos(t274);
t174 = t197 * t234 - t215 * t259;
t228 = -t236 * t289 - t239 * t244;
t229 = -t236 * t290 + t239 * t246;
t199 = t229 * t300 + (t228 * t298 + t236 * t260) * t243;
t216 = -t228 * t297 + t236 * t261;
t176 = t199 * t234 - t216 * t259;
t214 = t240 * t297 * t243 + (t243 * t246 * t298 + t244 * t300) * t237;
t225 = t240 * t298 - t246 * t260;
t190 = t214 * t234 - t225 * t259;
t175 = t197 * t259 + t215 * t234;
t251 = t300 * t297;
t249 = t237 * t251;
t252 = t298 * t300;
t196 = -t226 * t252 + t227 * t243 + t239 * t249;
t242 = sin(qJ(6));
t245 = cos(qJ(6));
t141 = -t175 * t242 + t196 * t245;
t142 = t175 * t245 + t196 * t242;
t94 = Icges(7,5) * t142 + Icges(7,6) * t141 + Icges(7,3) * t174;
t96 = Icges(7,4) * t142 + Icges(7,2) * t141 + Icges(7,6) * t174;
t98 = Icges(7,1) * t142 + Icges(7,4) * t141 + Icges(7,5) * t174;
t34 = t141 * t96 + t142 * t98 + t174 * t94;
t177 = t199 * t259 + t216 * t234;
t198 = -t228 * t252 + t229 * t243 - t236 * t249;
t143 = -t177 * t242 + t198 * t245;
t144 = t177 * t245 + t198 * t242;
t95 = Icges(7,5) * t144 + Icges(7,6) * t143 + Icges(7,3) * t176;
t97 = Icges(7,4) * t144 + Icges(7,2) * t143 + Icges(7,6) * t176;
t99 = Icges(7,1) * t144 + Icges(7,4) * t143 + Icges(7,5) * t176;
t35 = t141 * t97 + t142 * t99 + t174 * t95;
t191 = t214 * t259 + t225 * t234;
t291 = t237 * t244;
t213 = -t237 * t246 * t252 - t240 * t251 + t243 * t291;
t172 = -t191 * t242 + t213 * t245;
t173 = t191 * t245 + t213 * t242;
t118 = Icges(7,5) * t173 + Icges(7,6) * t172 + Icges(7,3) * t190;
t119 = Icges(7,4) * t173 + Icges(7,2) * t172 + Icges(7,6) * t190;
t120 = Icges(7,1) * t173 + Icges(7,4) * t172 + Icges(7,5) * t190;
t49 = t118 * t174 + t119 * t141 + t120 * t142;
t1 = t174 * t34 + t176 * t35 + t190 * t49;
t306 = t1 / 0.2e1;
t36 = t143 * t96 + t144 * t98 + t176 * t94;
t37 = t143 * t97 + t144 * t99 + t176 * t95;
t50 = t118 * t176 + t119 * t143 + t120 * t144;
t2 = t174 * t36 + t176 * t37 + t190 * t50;
t305 = t2 / 0.2e1;
t42 = t172 * t96 + t173 * t98 + t190 * t94;
t43 = t172 * t97 + t173 * t99 + t190 * t95;
t58 = t118 * t190 + t119 * t172 + t120 * t173;
t9 = t174 * t42 + t176 * t43 + t190 * t58;
t304 = t9 / 0.2e1;
t303 = t174 / 0.2e1;
t302 = t176 / 0.2e1;
t301 = t190 / 0.2e1;
t238 = cos(pkin(13));
t299 = pkin(4) * t238;
t235 = sin(pkin(13));
t296 = t215 * t235;
t295 = t216 * t235;
t294 = t225 * t235;
t293 = t236 * t237;
t292 = t237 * t239;
t100 = rSges(7,1) * t142 + rSges(7,2) * t141 + rSges(7,3) * t174;
t287 = pkin(5) * t175 + pkin(11) * t174 + t100;
t101 = rSges(7,1) * t144 + rSges(7,2) * t143 + rSges(7,3) * t176;
t286 = pkin(5) * t177 + pkin(11) * t176 + t101;
t116 = pkin(4) * t296 + pkin(10) * t196 + t197 * t299;
t169 = pkin(3) * t197 + qJ(4) * t196;
t162 = t216 * t169;
t285 = t116 * t216 + t162;
t117 = pkin(4) * t295 + pkin(10) * t198 + t199 * t299;
t170 = pkin(3) * t199 + qJ(4) * t198;
t164 = t225 * t170;
t284 = t117 * t225 + t164;
t283 = -t116 - t169;
t282 = -t117 - t170;
t121 = rSges(7,1) * t173 + rSges(7,2) * t172 + rSges(7,3) * t190;
t281 = pkin(5) * t191 + pkin(11) * t190 + t121;
t179 = -t197 * t235 + t215 * t238;
t180 = t197 * t238 + t296;
t136 = rSges(5,1) * t180 + rSges(5,2) * t179 + rSges(5,3) * t196;
t280 = -t136 - t169;
t145 = pkin(4) * t294 + pkin(10) * t213 + t214 * t299;
t189 = pkin(3) * t214 + qJ(4) * t213;
t178 = t215 * t189;
t279 = t145 * t215 + t178;
t278 = -t145 - t189;
t194 = -t214 * t235 + t225 * t238;
t195 = t214 * t238 + t294;
t153 = rSges(5,1) * t195 + rSges(5,2) * t194 + rSges(5,3) * t213;
t277 = -t153 - t189;
t202 = pkin(2) * t229 + pkin(9) * t216;
t200 = t240 * t202;
t276 = t170 * t240 + t200;
t201 = t227 * pkin(2) + pkin(9) * t215;
t275 = t201 * t293 + t202 * t292;
t122 = Icges(6,5) * t175 - Icges(6,6) * t174 + Icges(6,3) * t196;
t124 = Icges(6,4) * t175 - Icges(6,2) * t174 + Icges(6,6) * t196;
t126 = Icges(6,1) * t175 - Icges(6,4) * t174 + Icges(6,5) * t196;
t59 = t122 * t196 - t124 * t174 + t126 * t175;
t123 = Icges(6,5) * t177 - Icges(6,6) * t176 + Icges(6,3) * t198;
t125 = Icges(6,4) * t177 - Icges(6,2) * t176 + Icges(6,6) * t198;
t127 = Icges(6,1) * t177 - Icges(6,4) * t176 + Icges(6,5) * t198;
t60 = t123 * t196 - t125 * t174 + t127 * t175;
t146 = Icges(6,5) * t191 - Icges(6,6) * t190 + Icges(6,3) * t213;
t147 = Icges(6,4) * t191 - Icges(6,2) * t190 + Icges(6,6) * t213;
t148 = Icges(6,1) * t191 - Icges(6,4) * t190 + Icges(6,5) * t213;
t75 = t146 * t196 - t147 * t174 + t148 * t175;
t13 = t196 * t59 + t198 * t60 + t213 * t75;
t3 = t196 * t34 + t198 * t35 + t213 * t49;
t273 = t3 / 0.2e1 + t13 / 0.2e1;
t61 = t122 * t198 - t124 * t176 + t126 * t177;
t62 = t123 * t198 - t125 * t176 + t127 * t177;
t76 = t146 * t198 - t147 * t176 + t148 * t177;
t14 = t196 * t61 + t198 * t62 + t213 * t76;
t4 = t196 * t36 + t198 * t37 + t213 * t50;
t272 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t215 * t59 + t216 * t60 + t225 * t75;
t5 = t215 * t34 + t216 * t35 + t225 * t49;
t271 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = t215 * t61 + t216 * t62 + t225 * t76;
t6 = t215 * t36 + t216 * t37 + t225 * t50;
t270 = t6 / 0.2e1 + t16 / 0.2e1;
t17 = t240 * t75 + (t236 * t60 - t239 * t59) * t237;
t7 = t240 * t49 + (t236 * t35 - t239 * t34) * t237;
t269 = t7 / 0.2e1 + t17 / 0.2e1;
t18 = t240 * t76 + (t236 * t62 - t239 * t61) * t237;
t8 = t240 * t50 + (t236 * t37 - t239 * t36) * t237;
t268 = t8 / 0.2e1 + t18 / 0.2e1;
t10 = t196 * t42 + t198 * t43 + t213 * t58;
t67 = t122 * t213 - t124 * t190 + t126 * t191;
t68 = t123 * t213 - t125 * t190 + t127 * t191;
t84 = t146 * t213 - t147 * t190 + t148 * t191;
t23 = t196 * t67 + t198 * t68 + t213 * t84;
t267 = t10 / 0.2e1 + t23 / 0.2e1;
t11 = t215 * t42 + t216 * t43 + t225 * t58;
t24 = t215 * t67 + t216 * t68 + t225 * t84;
t266 = t11 / 0.2e1 + t24 / 0.2e1;
t12 = t240 * t58 + (t236 * t43 - t239 * t42) * t237;
t25 = t240 * t84 + (t236 * t68 - t239 * t67) * t237;
t265 = t12 / 0.2e1 + t25 / 0.2e1;
t264 = t117 * t240 + t276;
t128 = rSges(6,1) * t175 - rSges(6,2) * t174 + rSges(6,3) * t196;
t263 = -t128 + t283;
t149 = rSges(6,1) * t191 - rSges(6,2) * t190 + rSges(6,3) * t213;
t262 = -t149 + t278;
t186 = rSges(4,1) * t214 - rSges(4,2) * t213 + rSges(4,3) * t225;
t217 = pkin(2) * t291 + pkin(9) * t225;
t258 = (-t186 - t217) * t237;
t256 = t283 - t287;
t255 = t278 - t281;
t254 = t169 * t293 + t170 * t292 + t275;
t253 = (-t217 + t277) * t237;
t250 = (-t217 + t262) * t237;
t248 = t116 * t293 + t117 * t292 + t254;
t247 = (-t217 + t255) * t237;
t224 = t240 * rSges(3,3) + (rSges(3,1) * t244 + rSges(3,2) * t246) * t237;
t223 = Icges(3,5) * t240 + (Icges(3,1) * t244 + Icges(3,4) * t246) * t237;
t222 = Icges(3,6) * t240 + (Icges(3,4) * t244 + Icges(3,2) * t246) * t237;
t221 = Icges(3,3) * t240 + (Icges(3,5) * t244 + Icges(3,6) * t246) * t237;
t210 = rSges(3,1) * t229 + rSges(3,2) * t228 + rSges(3,3) * t293;
t209 = rSges(3,1) * t227 + rSges(3,2) * t226 - rSges(3,3) * t292;
t208 = Icges(3,1) * t229 + Icges(3,4) * t228 + Icges(3,5) * t293;
t207 = Icges(3,1) * t227 + Icges(3,4) * t226 - Icges(3,5) * t292;
t206 = Icges(3,4) * t229 + Icges(3,2) * t228 + Icges(3,6) * t293;
t205 = Icges(3,4) * t227 + Icges(3,2) * t226 - Icges(3,6) * t292;
t204 = Icges(3,5) * t229 + Icges(3,6) * t228 + Icges(3,3) * t293;
t203 = Icges(3,5) * t227 + Icges(3,6) * t226 - Icges(3,3) * t292;
t188 = -t209 * t240 - t224 * t292;
t187 = t210 * t240 - t224 * t293;
t185 = Icges(4,1) * t214 - Icges(4,4) * t213 + Icges(4,5) * t225;
t184 = Icges(4,4) * t214 - Icges(4,2) * t213 + Icges(4,6) * t225;
t183 = Icges(4,5) * t214 - Icges(4,6) * t213 + Icges(4,3) * t225;
t182 = t199 * t238 + t295;
t181 = -t199 * t235 + t216 * t238;
t171 = (t209 * t236 + t210 * t239) * t237;
t161 = rSges(4,1) * t199 - rSges(4,2) * t198 + rSges(4,3) * t216;
t160 = rSges(4,1) * t197 - rSges(4,2) * t196 + rSges(4,3) * t215;
t159 = Icges(4,1) * t199 - Icges(4,4) * t198 + Icges(4,5) * t216;
t158 = Icges(4,1) * t197 - Icges(4,4) * t196 + Icges(4,5) * t215;
t157 = Icges(4,4) * t199 - Icges(4,2) * t198 + Icges(4,6) * t216;
t156 = Icges(4,4) * t197 - Icges(4,2) * t196 + Icges(4,6) * t215;
t155 = Icges(4,5) * t199 - Icges(4,6) * t198 + Icges(4,3) * t216;
t154 = Icges(4,5) * t197 - Icges(4,6) * t196 + Icges(4,3) * t215;
t152 = Icges(5,1) * t195 + Icges(5,4) * t194 + Icges(5,5) * t213;
t151 = Icges(5,4) * t195 + Icges(5,2) * t194 + Icges(5,6) * t213;
t150 = Icges(5,5) * t195 + Icges(5,6) * t194 + Icges(5,3) * t213;
t137 = rSges(5,1) * t182 + rSges(5,2) * t181 + rSges(5,3) * t198;
t135 = Icges(5,1) * t182 + Icges(5,4) * t181 + Icges(5,5) * t198;
t134 = Icges(5,1) * t180 + Icges(5,4) * t179 + Icges(5,5) * t196;
t133 = Icges(5,4) * t182 + Icges(5,2) * t181 + Icges(5,6) * t198;
t132 = Icges(5,4) * t180 + Icges(5,2) * t179 + Icges(5,6) * t196;
t131 = Icges(5,5) * t182 + Icges(5,6) * t181 + Icges(5,3) * t198;
t130 = Icges(5,5) * t180 + Icges(5,6) * t179 + Icges(5,3) * t196;
t129 = rSges(6,1) * t177 - rSges(6,2) * t176 + rSges(6,3) * t198;
t112 = t161 * t225 - t186 * t216;
t111 = -t160 * t225 + t186 * t215;
t108 = (-t160 - t201) * t240 + t239 * t258;
t107 = t161 * t240 + t236 * t258 + t200;
t106 = t183 * t225 - t184 * t213 + t185 * t214;
t105 = t160 * t216 - t161 * t215;
t104 = t183 * t216 - t184 * t198 + t185 * t199;
t103 = t183 * t215 - t184 * t196 + t185 * t197;
t102 = (t160 * t236 + t161 * t239) * t237 + t275;
t93 = t129 * t213 - t149 * t198;
t92 = -t128 * t213 + t149 * t196;
t91 = t155 * t225 - t157 * t213 + t159 * t214;
t90 = t154 * t225 - t156 * t213 + t158 * t214;
t89 = t155 * t216 - t157 * t198 + t159 * t199;
t88 = t154 * t216 - t156 * t198 + t158 * t199;
t87 = t155 * t215 - t157 * t196 + t159 * t197;
t86 = t154 * t215 - t156 * t196 + t158 * t197;
t85 = t150 * t213 + t151 * t194 + t152 * t195;
t83 = t128 * t198 - t129 * t196;
t82 = t137 * t225 + t216 * t277 + t164;
t81 = t153 * t215 + t225 * t280 + t178;
t80 = (-t201 + t280) * t240 + t239 * t253;
t79 = t137 * t240 + t236 * t253 + t276;
t78 = t150 * t198 + t151 * t181 + t152 * t182;
t77 = t150 * t196 + t151 * t179 + t152 * t180;
t74 = t101 * t190 - t121 * t176;
t73 = -t100 * t190 + t121 * t174;
t72 = t136 * t216 + t162 + (-t137 - t170) * t215;
t71 = (t136 * t236 + t137 * t239) * t237 + t254;
t70 = t131 * t213 + t133 * t194 + t135 * t195;
t69 = t130 * t213 + t132 * t194 + t134 * t195;
t66 = t131 * t198 + t133 * t181 + t135 * t182;
t65 = t130 * t198 + t132 * t181 + t134 * t182;
t64 = t131 * t196 + t133 * t179 + t135 * t180;
t63 = t130 * t196 + t132 * t179 + t134 * t180;
t57 = t100 * t176 - t101 * t174;
t56 = (-t201 + t263) * t240 + t239 * t250;
t55 = t129 * t240 + t236 * t250 + t264;
t54 = t129 * t225 + t216 * t262 + t284;
t53 = t149 * t215 + t225 * t263 + t279;
t52 = -t198 * t281 + t213 * t286;
t51 = t196 * t281 - t213 * t287;
t48 = t106 * t240 + (t236 * t91 - t239 * t90) * t237;
t47 = (t128 * t236 + t129 * t239) * t237 + t248;
t46 = t106 * t225 + t215 * t90 + t216 * t91;
t45 = -t196 * t286 + t198 * t287;
t44 = t128 * t216 + (-t129 + t282) * t215 + t285;
t41 = t104 * t240 + (t236 * t89 - t239 * t88) * t237;
t40 = t103 * t240 + (t236 * t87 - t239 * t86) * t237;
t39 = t104 * t225 + t215 * t88 + t216 * t89;
t38 = t103 * t225 + t215 * t86 + t216 * t87;
t33 = (-t201 + t256) * t240 + t239 * t247;
t32 = t236 * t247 + t240 * t286 + t264;
t31 = t216 * t255 + t225 * t286 + t284;
t30 = t215 * t281 + t225 * t256 + t279;
t29 = (t236 * t287 + t239 * t286) * t237 + t248;
t28 = t287 * t216 + (t282 - t286) * t215 + t285;
t27 = t240 * t85 + (t236 * t70 - t239 * t69) * t237;
t26 = t215 * t69 + t216 * t70 + t225 * t85;
t22 = t240 * t78 + (t236 * t66 - t239 * t65) * t237;
t21 = t240 * t77 + (t236 * t64 - t239 * t63) * t237;
t20 = t215 * t65 + t216 * t66 + t225 * t78;
t19 = t215 * t63 + t216 * t64 + t225 * t77;
t109 = [m(4) + m(2) + m(3) + t307; m(3) * t171 + m(4) * t102 + m(5) * t71 + m(6) * t47 + m(7) * t29; m(7) * (t29 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t47 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t71 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(4) * (t102 ^ 2 + t107 ^ 2 + t108 ^ 2) + m(3) * (t171 ^ 2 + t187 ^ 2 + t188 ^ 2) + (t41 + t22 + t8 + t18 + (t204 * t293 + t206 * t228 + t208 * t229) * t293) * t293 + (-t17 - t40 - t21 - t7 + (-t203 * t292 + t205 * t226 + t207 * t227) * t292 + (-t203 * t293 + t204 * t292 - t205 * t228 - t206 * t226 - t207 * t229 - t208 * t227) * t293) * t292 + ((t221 * t293 + t222 * t228 + t223 * t229) * t293 - (-t221 * t292 + t222 * t226 + t223 * t227) * t292 + t12 + t25 + t27 + t48 + ((t206 * t246 + t208 * t244) * t236 - (t205 * t246 + t207 * t244) * t239) * t237 ^ 2 + ((-t203 * t239 + t204 * t236 + t222 * t246 + t223 * t244) * t237 + t240 * t221) * t240) * t240; m(4) * t105 + m(5) * t72 + m(6) * t44 + m(7) * t28; (t46 / 0.2e1 + t26 / 0.2e1 + t266) * t240 + (t27 / 0.2e1 + t48 / 0.2e1 + t265) * t225 + (t41 / 0.2e1 + t22 / 0.2e1 + t268) * t216 + (t40 / 0.2e1 + t21 / 0.2e1 + t269) * t215 + m(7) * (t28 * t29 + t30 * t33 + t31 * t32) + m(6) * (t44 * t47 + t53 * t56 + t54 * t55) + m(5) * (t71 * t72 + t79 * t82 + t80 * t81) + m(4) * (t102 * t105 + t107 * t112 + t108 * t111) + ((-t19 / 0.2e1 - t38 / 0.2e1 - t271) * t239 + (t20 / 0.2e1 + t39 / 0.2e1 + t270) * t236) * t237; (t11 + t24 + t46 + t26) * t225 + (t6 + t16 + t39 + t20) * t216 + (t5 + t15 + t19 + t38) * t215 + m(7) * (t28 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t44 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(5) * (t72 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t105 ^ 2 + t111 ^ 2 + t112 ^ 2); t213 * t307; m(7) * (t196 * t32 + t198 * t33 + t213 * t29) + m(6) * (t196 * t55 + t198 * t56 + t213 * t47) + m(5) * (t196 * t79 + t198 * t80 + t213 * t71); m(7) * (t196 * t31 + t198 * t30 + t213 * t28) + m(6) * (t196 * t54 + t198 * t53 + t213 * t44) + m(5) * (t196 * t82 + t198 * t81 + t213 * t72); (t196 ^ 2 + t198 ^ 2 + t213 ^ 2) * t307; m(6) * t83 + m(7) * t45; t267 * t240 + t265 * t213 + t268 * t198 + t269 * t196 + m(7) * (t29 * t45 + t32 * t52 + t33 * t51) + m(6) * (t47 * t83 + t55 * t93 + t56 * t92) + (t236 * t272 - t239 * t273) * t237; t267 * t225 + t272 * t216 + t273 * t215 + t266 * t213 + t270 * t198 + t271 * t196 + m(7) * (t28 * t45 + t30 * t51 + t31 * t52) + m(6) * (t44 * t83 + t53 * t92 + t54 * t93); m(6) * (t196 * t93 + t198 * t92 + t213 * t83) + m(7) * (t196 * t52 + t198 * t51 + t213 * t45); (t10 + t23) * t213 + (t4 + t14) * t198 + (t3 + t13) * t196 + m(7) * (t45 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(6) * (t83 ^ 2 + t92 ^ 2 + t93 ^ 2); m(7) * t57; m(7) * (t29 * t57 + t32 * t74 + t33 * t73) + t7 * t303 + t8 * t302 + t12 * t301 + t240 * t304 + (t236 * t305 - t239 * t1 / 0.2e1) * t237; m(7) * (t28 * t57 + t30 * t73 + t31 * t74) + t215 * t306 + t11 * t301 + t216 * t305 + t6 * t302 + t225 * t304 + t5 * t303; m(7) * (t196 * t74 + t198 * t73 + t213 * t57); t198 * t305 + t196 * t306 + t3 * t303 + t4 * t302 + t213 * t304 + m(7) * (t45 * t57 + t51 * t73 + t52 * t74) + t10 * t301; t176 * t2 + t174 * t1 + t190 * t9 + m(7) * (t57 ^ 2 + t73 ^ 2 + t74 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t109(1) t109(2) t109(4) t109(7) t109(11) t109(16); t109(2) t109(3) t109(5) t109(8) t109(12) t109(17); t109(4) t109(5) t109(6) t109(9) t109(13) t109(18); t109(7) t109(8) t109(9) t109(10) t109(14) t109(19); t109(11) t109(12) t109(13) t109(14) t109(15) t109(20); t109(16) t109(17) t109(18) t109(19) t109(20) t109(21);];
Mq  = res;
