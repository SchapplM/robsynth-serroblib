% Calculate joint inertia matrix for
% S6PRRPRR4
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:21
% EndTime: 2019-03-08 22:10:34
% DurationCPUTime: 6.20s
% Computational Cost: add. (23827->559), mult. (61849->811), div. (0->0), fcn. (80791->12), ass. (0->254)
t222 = sin(pkin(11));
t224 = cos(pkin(11));
t229 = sin(qJ(2));
t225 = cos(pkin(6));
t231 = cos(qJ(2));
t267 = t225 * t231;
t212 = t222 * t229 - t224 * t267;
t214 = t222 * t267 + t224 * t229;
t223 = sin(pkin(6));
t269 = t223 * t231;
t268 = t225 * t229;
t213 = t222 * t231 + t224 * t268;
t228 = sin(qJ(3));
t277 = cos(qJ(3));
t243 = t223 * t277;
t202 = t213 * t228 + t224 * t243;
t270 = t223 * t228;
t203 = t213 * t277 - t224 * t270;
t227 = sin(qJ(5));
t276 = cos(qJ(5));
t160 = -t202 * t276 + t203 * t227;
t161 = t202 * t227 + t203 * t276;
t109 = Icges(6,5) * t161 - Icges(6,6) * t160 - Icges(6,3) * t212;
t111 = Icges(6,4) * t161 - Icges(6,2) * t160 - Icges(6,6) * t212;
t113 = Icges(6,1) * t161 - Icges(6,4) * t160 - Icges(6,5) * t212;
t54 = -t109 * t212 - t111 * t160 + t113 * t161;
t215 = -t222 * t268 + t224 * t231;
t204 = t215 * t228 - t222 * t243;
t205 = t215 * t277 + t222 * t270;
t162 = -t204 * t276 + t205 * t227;
t163 = t204 * t227 + t205 * t276;
t110 = Icges(6,5) * t163 - Icges(6,6) * t162 - Icges(6,3) * t214;
t112 = Icges(6,4) * t163 - Icges(6,2) * t162 - Icges(6,6) * t214;
t114 = Icges(6,1) * t163 - Icges(6,4) * t162 - Icges(6,5) * t214;
t55 = -t110 * t212 - t112 * t160 + t114 * t161;
t216 = -t225 * t277 + t229 * t270;
t217 = t225 * t228 + t229 * t243;
t195 = -t216 * t276 + t217 * t227;
t196 = t216 * t227 + t217 * t276;
t126 = Icges(6,5) * t196 - Icges(6,6) * t195 + Icges(6,3) * t269;
t127 = Icges(6,4) * t196 - Icges(6,2) * t195 + Icges(6,6) * t269;
t128 = Icges(6,1) * t196 - Icges(6,4) * t195 + Icges(6,5) * t269;
t67 = -t126 * t212 - t127 * t160 + t128 * t161;
t14 = -t54 * t212 - t55 * t214 + t269 * t67;
t226 = sin(qJ(6));
t230 = cos(qJ(6));
t129 = -t161 * t226 - t212 * t230;
t130 = t161 * t230 - t212 * t226;
t92 = Icges(7,5) * t130 + Icges(7,6) * t129 + Icges(7,3) * t160;
t94 = Icges(7,4) * t130 + Icges(7,2) * t129 + Icges(7,6) * t160;
t96 = Icges(7,1) * t130 + Icges(7,4) * t129 + Icges(7,5) * t160;
t24 = t129 * t94 + t130 * t96 + t160 * t92;
t131 = -t163 * t226 - t214 * t230;
t132 = t163 * t230 - t214 * t226;
t93 = Icges(7,5) * t132 + Icges(7,6) * t131 + Icges(7,3) * t162;
t95 = Icges(7,4) * t132 + Icges(7,2) * t131 + Icges(7,6) * t162;
t97 = Icges(7,1) * t132 + Icges(7,4) * t131 + Icges(7,5) * t162;
t25 = t129 * t95 + t130 * t97 + t160 * t93;
t171 = -t196 * t226 + t230 * t269;
t172 = t196 * t230 + t226 * t269;
t117 = Icges(7,5) * t172 + Icges(7,6) * t171 + Icges(7,3) * t195;
t118 = Icges(7,4) * t172 + Icges(7,2) * t171 + Icges(7,6) * t195;
t119 = Icges(7,1) * t172 + Icges(7,4) * t171 + Icges(7,5) * t195;
t47 = t117 * t160 + t118 * t129 + t119 * t130;
t4 = -t24 * t212 - t25 * t214 + t269 * t47;
t288 = -t4 - t14;
t56 = -t109 * t214 - t111 * t162 + t113 * t163;
t57 = -t110 * t214 - t112 * t162 + t114 * t163;
t68 = -t126 * t214 - t127 * t162 + t128 * t163;
t16 = -t56 * t212 - t57 * t214 + t269 * t68;
t26 = t131 * t94 + t132 * t96 + t162 * t92;
t27 = t131 * t95 + t132 * t97 + t162 * t93;
t48 = t117 * t162 + t118 * t131 + t119 * t132;
t6 = -t26 * t212 - t27 * t214 + t269 * t48;
t287 = -t6 - t16;
t28 = t171 * t94 + t172 * t96 + t195 * t92;
t29 = t171 * t95 + t172 * t97 + t195 * t93;
t58 = t117 * t195 + t118 * t171 + t119 * t172;
t11 = -t28 * t212 - t29 * t214 + t269 * t58;
t59 = t109 * t269 - t195 * t111 + t196 * t113;
t60 = t110 * t269 - t195 * t112 + t196 * t114;
t70 = t126 * t269 - t195 * t127 + t196 * t128;
t20 = -t59 * t212 - t60 * t214 + t269 * t70;
t286 = t11 + t20;
t283 = m(5) + m(6) + m(7);
t254 = -t4 / 0.2e1 - t14 / 0.2e1;
t285 = t6 / 0.2e1 + t16 / 0.2e1;
t249 = t11 / 0.2e1 + t20 / 0.2e1;
t1 = t160 * t24 + t162 * t25 + t195 * t47;
t278 = t195 / 0.2e1;
t279 = t162 / 0.2e1;
t280 = t160 / 0.2e1;
t2 = t160 * t26 + t162 * t27 + t195 * t48;
t281 = t2 / 0.2e1;
t9 = t160 * t28 + t162 * t29 + t195 * t58;
t284 = -t11 * t278 - t6 * t279 - t4 * t280 + t212 * t1 / 0.2e1 + t214 * t281 - t9 * t269 / 0.2e1;
t98 = rSges(7,1) * t130 + rSges(7,2) * t129 + rSges(7,3) * t160;
t275 = pkin(5) * t161 + pkin(10) * t160 + t98;
t282 = t214 * t275;
t120 = rSges(7,1) * t172 + rSges(7,2) * t171 + rSges(7,3) * t195;
t266 = pkin(5) * t196 + pkin(10) * t195 + t120;
t49 = -t212 * t266 - t269 * t275;
t99 = rSges(7,1) * t132 + rSges(7,2) * t131 + rSges(7,3) * t162;
t274 = pkin(5) * t163 + pkin(10) * t162 + t99;
t115 = rSges(6,1) * t161 - rSges(6,2) * t160 - rSges(6,3) * t212;
t273 = t214 * t115;
t272 = t222 * t223;
t271 = t223 * t224;
t149 = rSges(5,1) * t205 + rSges(5,2) * t214 + rSges(5,3) * t204;
t165 = pkin(3) * t205 + qJ(4) * t204;
t265 = -t149 - t165;
t164 = pkin(3) * t203 + qJ(4) * t202;
t151 = t214 * t164;
t173 = pkin(4) * t203 - pkin(9) * t212;
t264 = t214 * t173 + t151;
t199 = pkin(3) * t217 + qJ(4) * t216;
t263 = t164 * t269 + t212 * t199;
t198 = pkin(2) * t215 + pkin(8) * t214;
t194 = t225 * t198;
t262 = t225 * t165 + t194;
t197 = pkin(2) * t213 + pkin(8) * t212;
t261 = -t164 - t197;
t174 = pkin(4) * t205 - pkin(9) * t214;
t260 = -t165 - t174;
t189 = t217 * rSges(5,1) - rSges(5,2) * t269 + t216 * rSges(5,3);
t259 = -t189 - t199;
t258 = t197 * t272 + t198 * t271;
t206 = t217 * pkin(4) + pkin(9) * t269;
t257 = -t199 - t206;
t17 = t67 * t225 + (t222 * t55 - t224 * t54) * t223;
t7 = t47 * t225 + (t222 * t25 - t224 * t24) * t223;
t252 = -t7 / 0.2e1 - t17 / 0.2e1;
t18 = t68 * t225 + (t222 * t57 - t224 * t56) * t223;
t8 = t48 * t225 + (t222 * t27 - t224 * t26) * t223;
t251 = -t8 / 0.2e1 - t18 / 0.2e1;
t12 = t58 * t225 + (t222 * t29 - t224 * t28) * t223;
t21 = t70 * t225 + (t222 * t60 - t224 * t59) * t223;
t248 = t12 / 0.2e1 + t21 / 0.2e1;
t116 = rSges(6,1) * t163 - rSges(6,2) * t162 - rSges(6,3) * t214;
t247 = -t116 + t260;
t133 = t196 * rSges(6,1) - t195 * rSges(6,2) + rSges(6,3) * t269;
t246 = -t133 + t257;
t245 = t225 * t174 + t262;
t244 = -t173 + t261;
t190 = t217 * rSges(4,1) - t216 * rSges(4,2) - rSges(4,3) * t269;
t218 = (pkin(2) * t229 - pkin(8) * t231) * t223;
t242 = (-t190 - t218) * t223;
t240 = t260 - t274;
t239 = t257 - t266;
t238 = t164 * t272 + t165 * t271 + t258;
t237 = t173 * t269 + t212 * t206 + t263;
t236 = (-t218 + t259) * t223;
t235 = (-t218 + t246) * t223;
t234 = t173 * t272 + t174 * t271 + t238;
t89 = -t115 * t269 - t212 * t133;
t233 = (-t218 + t239) * t223;
t210 = t225 * rSges(3,3) + (rSges(3,1) * t229 + rSges(3,2) * t231) * t223;
t209 = Icges(3,5) * t225 + (Icges(3,1) * t229 + Icges(3,4) * t231) * t223;
t208 = Icges(3,6) * t225 + (Icges(3,4) * t229 + Icges(3,2) * t231) * t223;
t207 = Icges(3,3) * t225 + (Icges(3,5) * t229 + Icges(3,6) * t231) * t223;
t188 = Icges(4,1) * t217 - Icges(4,4) * t216 - Icges(4,5) * t269;
t187 = Icges(5,1) * t217 - Icges(5,4) * t269 + Icges(5,5) * t216;
t186 = Icges(4,4) * t217 - Icges(4,2) * t216 - Icges(4,6) * t269;
t185 = Icges(5,4) * t217 - Icges(5,2) * t269 + Icges(5,6) * t216;
t184 = Icges(4,5) * t217 - Icges(4,6) * t216 - Icges(4,3) * t269;
t183 = Icges(5,5) * t217 - Icges(5,6) * t269 + Icges(5,3) * t216;
t182 = rSges(3,1) * t215 - rSges(3,2) * t214 + rSges(3,3) * t272;
t181 = rSges(3,1) * t213 - rSges(3,2) * t212 - rSges(3,3) * t271;
t180 = Icges(3,1) * t215 - Icges(3,4) * t214 + Icges(3,5) * t272;
t179 = Icges(3,1) * t213 - Icges(3,4) * t212 - Icges(3,5) * t271;
t178 = Icges(3,4) * t215 - Icges(3,2) * t214 + Icges(3,6) * t272;
t177 = Icges(3,4) * t213 - Icges(3,2) * t212 - Icges(3,6) * t271;
t176 = Icges(3,5) * t215 - Icges(3,6) * t214 + Icges(3,3) * t272;
t175 = Icges(3,5) * t213 - Icges(3,6) * t212 - Icges(3,3) * t271;
t154 = -t181 * t225 - t210 * t271;
t153 = t182 * t225 - t210 * t272;
t150 = rSges(4,1) * t205 - rSges(4,2) * t204 + rSges(4,3) * t214;
t148 = rSges(4,1) * t203 - rSges(4,2) * t202 + rSges(4,3) * t212;
t147 = rSges(5,1) * t203 + rSges(5,2) * t212 + rSges(5,3) * t202;
t146 = Icges(4,1) * t205 - Icges(4,4) * t204 + Icges(4,5) * t214;
t145 = Icges(4,1) * t203 - Icges(4,4) * t202 + Icges(4,5) * t212;
t144 = Icges(5,1) * t205 + Icges(5,4) * t214 + Icges(5,5) * t204;
t143 = Icges(5,1) * t203 + Icges(5,4) * t212 + Icges(5,5) * t202;
t142 = Icges(4,4) * t205 - Icges(4,2) * t204 + Icges(4,6) * t214;
t141 = Icges(4,4) * t203 - Icges(4,2) * t202 + Icges(4,6) * t212;
t140 = Icges(5,4) * t205 + Icges(5,2) * t214 + Icges(5,6) * t204;
t139 = Icges(5,4) * t203 + Icges(5,2) * t212 + Icges(5,6) * t202;
t138 = Icges(4,5) * t205 - Icges(4,6) * t204 + Icges(4,3) * t214;
t137 = Icges(4,5) * t203 - Icges(4,6) * t202 + Icges(4,3) * t212;
t136 = Icges(5,5) * t205 + Icges(5,6) * t214 + Icges(5,3) * t204;
t135 = Icges(5,5) * t203 + Icges(5,6) * t212 + Icges(5,3) * t202;
t125 = (t181 * t222 + t182 * t224) * t223;
t122 = -t150 * t269 - t214 * t190;
t121 = t148 * t269 + t212 * t190;
t108 = -t184 * t269 - t216 * t186 + t217 * t188;
t107 = t216 * t183 - t185 * t269 + t217 * t187;
t106 = t148 * t214 - t150 * t212;
t105 = (-t148 - t197) * t225 + t224 * t242;
t104 = t225 * t150 + t222 * t242 + t194;
t103 = t184 * t214 - t186 * t204 + t188 * t205;
t102 = t183 * t204 + t185 * t214 + t187 * t205;
t101 = t184 * t212 - t186 * t202 + t188 * t203;
t100 = t183 * t202 + t185 * t212 + t187 * t203;
t91 = (t148 * t222 + t150 * t224) * t223 + t258;
t90 = t116 * t269 + t214 * t133;
t88 = t214 * t259 + t265 * t269;
t87 = t147 * t269 + t212 * t189 + t263;
t86 = -t138 * t269 - t216 * t142 + t217 * t146;
t85 = -t137 * t269 - t216 * t141 + t217 * t145;
t84 = t216 * t136 - t140 * t269 + t217 * t144;
t83 = t216 * t135 - t139 * t269 + t217 * t143;
t82 = (-t147 + t261) * t225 + t224 * t236;
t81 = t225 * t149 + t222 * t236 + t262;
t80 = t138 * t214 - t142 * t204 + t146 * t205;
t79 = t137 * t214 - t141 * t204 + t145 * t205;
t78 = t136 * t204 + t140 * t214 + t144 * t205;
t77 = t135 * t204 + t139 * t214 + t143 * t205;
t76 = t138 * t212 - t142 * t202 + t146 * t203;
t75 = t137 * t212 - t141 * t202 + t145 * t203;
t74 = t136 * t202 + t140 * t212 + t144 * t203;
t73 = t135 * t202 + t139 * t212 + t143 * t203;
t72 = t116 * t212 - t273;
t71 = t214 * t147 + t212 * t265 + t151;
t69 = (t147 * t222 + t149 * t224) * t223 + t238;
t66 = t214 * t246 + t247 * t269;
t65 = -t89 + t237;
t64 = -t120 * t162 + t195 * t99;
t63 = t120 * t160 - t195 * t98;
t62 = (-t115 + t244) * t225 + t224 * t235;
t61 = t225 * t116 + t222 * t235 + t245;
t53 = t212 * t247 + t264 + t273;
t52 = -t160 * t99 + t162 * t98;
t51 = (t115 * t222 + t116 * t224) * t223 + t234;
t50 = t214 * t266 + t269 * t274;
t46 = t108 * t225 + (t222 * t86 - t224 * t85) * t223;
t45 = t107 * t225 + (t222 * t84 - t224 * t83) * t223;
t44 = -t108 * t269 + t85 * t212 + t86 * t214;
t43 = -t107 * t269 + t83 * t212 + t84 * t214;
t42 = t212 * t274 - t282;
t41 = t103 * t225 + (t222 * t80 - t224 * t79) * t223;
t40 = t102 * t225 + (t222 * t78 - t224 * t77) * t223;
t39 = t101 * t225 + (t222 * t76 - t224 * t75) * t223;
t38 = t100 * t225 + (t222 * t74 - t224 * t73) * t223;
t37 = t214 * t239 + t240 * t269;
t36 = t237 - t49;
t35 = -t103 * t269 + t79 * t212 + t80 * t214;
t34 = -t102 * t269 + t77 * t212 + t78 * t214;
t33 = -t101 * t269 + t75 * t212 + t76 * t214;
t32 = -t100 * t269 + t73 * t212 + t74 * t214;
t31 = (t244 - t275) * t225 + t224 * t233;
t30 = t222 * t233 + t225 * t274 + t245;
t23 = t212 * t240 + t264 + t282;
t22 = (t222 * t275 + t224 * t274) * t223 + t234;
t3 = [m(2) + m(3) + m(4) + t283; m(3) * t125 + m(4) * t91 + m(5) * t69 + m(6) * t51 + m(7) * t22; m(7) * (t22 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t51 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t69 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t104 ^ 2 + t105 ^ 2 + t91 ^ 2) + m(3) * (t125 ^ 2 + t153 ^ 2 + t154 ^ 2) + (t8 + t18 + t41 + t40 + (t176 * t272 - t178 * t214 + t180 * t215) * t272) * t272 + (-t7 - t17 - t39 - t38 + (-t175 * t271 - t177 * t212 + t179 * t213) * t271 + (-t175 * t272 + t176 * t271 + t177 * t214 + t178 * t212 - t179 * t215 - t180 * t213) * t272) * t271 + ((t207 * t272 - t214 * t208 + t215 * t209) * t272 - (-t207 * t271 - t212 * t208 + t213 * t209) * t271 + t12 + t21 + t45 + t46 + ((t178 * t231 + t180 * t229) * t222 - (t177 * t231 + t179 * t229) * t224) * t223 ^ 2 + ((-t175 * t224 + t176 * t222 + t208 * t231 + t209 * t229) * t223 + t225 * t207) * t225) * t225; m(4) * t106 + m(5) * t71 + m(6) * t53 + m(7) * t23; (t44 / 0.2e1 + t43 / 0.2e1 - t249) * t225 + (t40 / 0.2e1 + t41 / 0.2e1 - t251) * t214 + (t38 / 0.2e1 + t39 / 0.2e1 - t252) * t212 + m(7) * (t22 * t23 + t30 * t37 + t31 * t36) + m(6) * (t51 * t53 + t61 * t66 + t62 * t65) + m(5) * (t69 * t71 + t81 * t88 + t82 * t87) + m(4) * (t104 * t122 + t105 * t121 + t106 * t91) + ((-t45 / 0.2e1 - t46 / 0.2e1 - t248) * t231 + (-t32 / 0.2e1 - t33 / 0.2e1 - t254) * t224 + (t34 / 0.2e1 + t35 / 0.2e1 - t285) * t222) * t223; (-t43 - t44 + t286) * t269 + (t35 + t34 + t287) * t214 + (t33 + t32 + t288) * t212 + m(7) * (t23 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t53 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t71 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t106 ^ 2 + t121 ^ 2 + t122 ^ 2); t216 * t283; m(7) * (t202 * t30 + t204 * t31 + t216 * t22) + m(6) * (t202 * t61 + t204 * t62 + t216 * t51) + m(5) * (t202 * t81 + t204 * t82 + t216 * t69); m(7) * (t202 * t37 + t204 * t36 + t216 * t23) + m(6) * (t202 * t66 + t204 * t65 + t216 * t53) + m(5) * (t202 * t88 + t204 * t87 + t216 * t71); (t202 ^ 2 + t204 ^ 2 + t216 ^ 2) * t283; m(6) * t72 + m(7) * t42; t249 * t225 + t251 * t214 + t252 * t212 + m(7) * (t22 * t42 + t30 * t50 + t31 * t49) + m(6) * (t51 * t72 + t61 * t90 + t62 * t89) + (t222 * t285 + t224 * t254 + t231 * t248) * t223; m(7) * (t23 * t42 + t36 * t49 + t37 * t50) + m(6) * (t53 * t72 + t65 * t89 + t66 * t90) - 0.2e1 * t249 * t269 + 0.2e1 * t285 * t214 - 0.2e1 * t254 * t212; m(6) * (t202 * t90 + t204 * t89 + t216 * t72) + m(7) * (t202 * t50 + t204 * t49 + t216 * t42); t286 * t269 + t287 * t214 + t288 * t212 + m(7) * (t42 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t72 ^ 2 + t89 ^ 2 + t90 ^ 2); m(7) * t52; t8 * t279 + t225 * t9 / 0.2e1 + m(7) * (t22 * t52 + t30 * t64 + t31 * t63) + t12 * t278 + t7 * t280 + (t222 * t281 - t224 * t1 / 0.2e1) * t223; m(7) * (t23 * t52 + t36 * t63 + t37 * t64) + t284; m(7) * (t202 * t64 + t204 * t63 + t216 * t52); m(7) * (t42 * t52 + t49 * t63 + t50 * t64) - t284; m(7) * (t52 ^ 2 + t63 ^ 2 + t64 ^ 2) + t162 * t2 + t160 * t1 + t195 * t9;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
