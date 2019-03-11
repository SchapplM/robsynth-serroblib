% Calculate joint inertia matrix for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:48
% EndTime: 2019-03-08 20:27:01
% DurationCPUTime: 6.28s
% Computational Cost: add. (31425->559), mult. (73508->819), div. (0->0), fcn. (97120->14), ass. (0->249)
t230 = sin(pkin(11));
t232 = cos(pkin(11));
t270 = sin(pkin(12));
t271 = cos(pkin(12));
t277 = sin(qJ(2));
t279 = cos(qJ(2));
t219 = -t277 * t270 + t279 * t271;
t233 = cos(pkin(6));
t238 = t233 * t219;
t239 = t270 * t279 + t271 * t277;
t198 = -t230 * t239 + t232 * t238;
t211 = t239 * t233;
t199 = t211 * t232 + t219 * t230;
t231 = sin(pkin(6));
t265 = t231 * t232;
t152 = Icges(4,5) * t199 + Icges(4,6) * t198 - Icges(4,3) * t265;
t252 = t233 * t279;
t214 = -t230 * t277 + t232 * t252;
t251 = t233 * t277;
t215 = t230 * t279 + t232 * t251;
t187 = Icges(3,5) * t215 + Icges(3,6) * t214 - Icges(3,3) * t265;
t289 = -t152 - t187;
t209 = t219 * t231;
t210 = t239 * t231;
t288 = Icges(4,5) * t210 + Icges(4,6) * t209 + (Icges(3,5) * t277 + Icges(3,6) * t279) * t231 + (Icges(4,3) + Icges(3,3)) * t233;
t200 = -t230 * t238 - t232 * t239;
t201 = -t211 * t230 + t219 * t232;
t266 = t230 * t231;
t153 = Icges(4,5) * t201 + Icges(4,6) * t200 + Icges(4,3) * t266;
t216 = -t230 * t252 - t232 * t277;
t217 = -t230 * t251 + t232 * t279;
t188 = Icges(3,5) * t217 + Icges(3,6) * t216 + Icges(3,3) * t266;
t287 = t188 + t153;
t235 = sin(qJ(4));
t278 = cos(qJ(4));
t253 = t231 * t278;
t176 = t199 * t235 + t232 * t253;
t286 = t176 / 0.2e1;
t178 = t201 * t235 - t230 * t253;
t285 = t178 / 0.2e1;
t284 = -t198 / 0.2e1;
t283 = -t200 / 0.2e1;
t202 = t210 * t235 - t233 * t278;
t282 = t202 / 0.2e1;
t281 = -t209 / 0.2e1;
t280 = t233 / 0.2e1;
t276 = t279 * pkin(2);
t236 = cos(qJ(5));
t275 = pkin(5) * t236;
t264 = t231 * t235;
t177 = t199 * t278 - t232 * t264;
t229 = qJ(5) + qJ(6);
t226 = sin(t229);
t227 = cos(t229);
t142 = -t177 * t226 - t198 * t227;
t143 = t177 * t227 - t198 * t226;
t105 = rSges(7,1) * t143 + rSges(7,2) * t142 + rSges(7,3) * t176;
t234 = sin(qJ(5));
t269 = t198 * t234;
t97 = -pkin(5) * t269 + pkin(10) * t176 + t177 * t275;
t273 = t105 + t97;
t179 = t201 * t278 + t230 * t264;
t144 = -t179 * t226 - t200 * t227;
t145 = t179 * t227 - t200 * t226;
t106 = rSges(7,1) * t145 + rSges(7,2) * t144 + rSges(7,3) * t178;
t268 = t200 * t234;
t98 = -pkin(5) * t268 + pkin(10) * t178 + t179 * t275;
t272 = t106 + t98;
t267 = t209 * t234;
t203 = t210 * t278 + t233 * t235;
t172 = -t203 * t226 - t209 * t227;
t173 = t203 * t227 - t209 * t226;
t121 = rSges(7,1) * t173 + rSges(7,2) * t172 + rSges(7,3) * t202;
t122 = -pkin(5) * t267 + pkin(10) * t202 + t203 * t275;
t263 = -t121 - t122;
t164 = pkin(3) * t201 - pkin(8) * t200;
t250 = pkin(2) * t251 - qJ(3) * t231;
t196 = -t230 * t250 + t232 * t276;
t186 = t233 * t196;
t262 = t233 * t164 + t186;
t163 = pkin(3) * t199 - pkin(8) * t198;
t195 = t230 * t276 + t232 * t250;
t261 = -t163 - t195;
t260 = t195 * t266 + t196 * t265;
t220 = pkin(2) * t231 * t277 + t233 * qJ(3);
t259 = -pkin(3) * t210 + pkin(8) * t209 - t220;
t101 = Icges(7,4) * t143 + Icges(7,2) * t142 + Icges(7,6) * t176;
t103 = Icges(7,1) * t143 + Icges(7,4) * t142 + Icges(7,5) * t176;
t99 = Icges(7,5) * t143 + Icges(7,6) * t142 + Icges(7,3) * t176;
t49 = t101 * t144 + t103 * t145 + t178 * t99;
t100 = Icges(7,5) * t145 + Icges(7,6) * t144 + Icges(7,3) * t178;
t102 = Icges(7,4) * t145 + Icges(7,2) * t144 + Icges(7,6) * t178;
t104 = Icges(7,1) * t145 + Icges(7,4) * t144 + Icges(7,5) * t178;
t50 = t100 * t178 + t102 * t144 + t104 * t145;
t118 = Icges(7,5) * t173 + Icges(7,6) * t172 + Icges(7,3) * t202;
t119 = Icges(7,4) * t173 + Icges(7,2) * t172 + Icges(7,6) * t202;
t120 = Icges(7,1) * t173 + Icges(7,4) * t172 + Icges(7,5) * t202;
t66 = t118 * t178 + t119 * t144 + t120 * t145;
t10 = t176 * t49 + t178 * t50 + t202 * t66;
t58 = t101 * t172 + t103 * t173 + t202 * t99;
t59 = t100 * t202 + t102 * t172 + t104 * t173;
t73 = t118 * t202 + t119 * t172 + t120 * t173;
t26 = t176 * t58 + t178 * t59 + t202 * t73;
t47 = t101 * t142 + t103 * t143 + t176 * t99;
t48 = t100 * t176 + t102 * t142 + t104 * t143;
t65 = t118 * t176 + t119 * t142 + t120 * t143;
t9 = t176 * t47 + t178 * t48 + t202 * t65;
t258 = t178 * t10 + t176 * t9 + t202 * t26;
t257 = m(4) + m(5) + m(6) + m(7);
t141 = t179 * pkin(4) + t178 * pkin(9);
t256 = t233 * t141 + t262;
t140 = t177 * pkin(4) + t176 * pkin(9);
t255 = -t140 + t261;
t171 = t203 * pkin(4) + t202 * pkin(9);
t254 = -t171 + t259;
t249 = (-rSges(4,1) * t210 - rSges(4,2) * t209 - rSges(4,3) * t233 - t220) * t231;
t248 = t163 * t266 + t164 * t265 + t260;
t168 = rSges(5,1) * t203 - rSges(5,2) * t202 - rSges(5,3) * t209;
t247 = (-t168 + t259) * t231;
t11 = -t198 * t47 - t200 * t48 - t209 * t65;
t12 = -t198 * t49 - t200 * t50 - t209 * t66;
t28 = -t198 * t58 - t200 * t59 - t209 * t73;
t244 = t10 * t283 + t11 * t286 + t12 * t285 + t26 * t281 + t28 * t282 + t9 * t284;
t17 = t233 * t65 + (t230 * t48 - t232 * t47) * t231;
t18 = t233 * t66 + (t230 * t50 - t232 * t49) * t231;
t31 = t233 * t73 + (t230 * t59 - t232 * t58) * t231;
t243 = t17 * t286 + t18 * t285 + t26 * t280 + t31 * t282 + t10 * t266 / 0.2e1 - t9 * t265 / 0.2e1;
t174 = -t203 * t234 - t209 * t236;
t175 = t203 * t236 - t267;
t135 = rSges(6,1) * t175 + rSges(6,2) * t174 + rSges(6,3) * t202;
t242 = (-t135 + t254) * t231;
t241 = t140 * t266 + t141 * t265 + t248;
t240 = (t254 + t263) * t231;
t208 = t233 * rSges(3,3) + (rSges(3,1) * t277 + rSges(3,2) * t279) * t231;
t207 = Icges(3,5) * t233 + (Icges(3,1) * t277 + Icges(3,4) * t279) * t231;
t206 = Icges(3,6) * t233 + (Icges(3,4) * t277 + Icges(3,2) * t279) * t231;
t194 = rSges(3,1) * t217 + rSges(3,2) * t216 + rSges(3,3) * t266;
t193 = rSges(3,1) * t215 + rSges(3,2) * t214 - rSges(3,3) * t265;
t192 = Icges(3,1) * t217 + Icges(3,4) * t216 + Icges(3,5) * t266;
t191 = Icges(3,1) * t215 + Icges(3,4) * t214 - Icges(3,5) * t265;
t190 = Icges(3,4) * t217 + Icges(3,2) * t216 + Icges(3,6) * t266;
t189 = Icges(3,4) * t215 + Icges(3,2) * t214 - Icges(3,6) * t265;
t184 = Icges(4,1) * t210 + Icges(4,4) * t209 + Icges(4,5) * t233;
t183 = Icges(4,4) * t210 + Icges(4,2) * t209 + Icges(4,6) * t233;
t170 = -t193 * t233 - t208 * t265;
t169 = t194 * t233 - t208 * t266;
t167 = Icges(5,1) * t203 - Icges(5,4) * t202 - Icges(5,5) * t209;
t166 = Icges(5,4) * t203 - Icges(5,2) * t202 - Icges(5,6) * t209;
t165 = Icges(5,5) * t203 - Icges(5,6) * t202 - Icges(5,3) * t209;
t159 = rSges(4,1) * t201 + rSges(4,2) * t200 + rSges(4,3) * t266;
t158 = rSges(4,1) * t199 + rSges(4,2) * t198 - rSges(4,3) * t265;
t157 = Icges(4,1) * t201 + Icges(4,4) * t200 + Icges(4,5) * t266;
t156 = Icges(4,1) * t199 + Icges(4,4) * t198 - Icges(4,5) * t265;
t155 = Icges(4,4) * t201 + Icges(4,2) * t200 + Icges(4,6) * t266;
t154 = Icges(4,4) * t199 + Icges(4,2) * t198 - Icges(4,6) * t265;
t151 = t179 * t236 - t268;
t150 = -t179 * t234 - t200 * t236;
t149 = t177 * t236 - t269;
t148 = -t177 * t234 - t198 * t236;
t147 = (t193 * t230 + t194 * t232) * t231;
t146 = t198 * t171;
t136 = t209 * t141;
t134 = Icges(6,1) * t175 + Icges(6,4) * t174 + Icges(6,5) * t202;
t133 = Icges(6,4) * t175 + Icges(6,2) * t174 + Icges(6,6) * t202;
t132 = Icges(6,5) * t175 + Icges(6,6) * t174 + Icges(6,3) * t202;
t131 = t200 * t140;
t130 = rSges(5,1) * t179 - rSges(5,2) * t178 - rSges(5,3) * t200;
t129 = rSges(5,1) * t177 - rSges(5,2) * t176 - rSges(5,3) * t198;
t128 = Icges(5,1) * t179 - Icges(5,4) * t178 - Icges(5,5) * t200;
t127 = Icges(5,1) * t177 - Icges(5,4) * t176 - Icges(5,5) * t198;
t126 = Icges(5,4) * t179 - Icges(5,2) * t178 - Icges(5,6) * t200;
t125 = Icges(5,4) * t177 - Icges(5,2) * t176 - Icges(5,6) * t198;
t124 = Icges(5,5) * t179 - Icges(5,6) * t178 - Icges(5,3) * t200;
t123 = Icges(5,5) * t177 - Icges(5,6) * t176 - Icges(5,3) * t198;
t117 = t176 * t121;
t116 = (-t158 - t195) * t233 + t232 * t249;
t115 = t159 * t233 + t230 * t249 + t186;
t114 = rSges(6,1) * t151 + rSges(6,2) * t150 + rSges(6,3) * t178;
t113 = rSges(6,1) * t149 + rSges(6,2) * t148 + rSges(6,3) * t176;
t112 = Icges(6,1) * t151 + Icges(6,4) * t150 + Icges(6,5) * t178;
t111 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t176;
t110 = Icges(6,4) * t151 + Icges(6,2) * t150 + Icges(6,6) * t178;
t109 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t176;
t108 = Icges(6,5) * t151 + Icges(6,6) * t150 + Icges(6,3) * t178;
t107 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t176;
t96 = (t158 * t230 + t159 * t232) * t231 + t260;
t95 = -t130 * t209 + t168 * t200;
t94 = t129 * t209 - t168 * t198;
t93 = t202 * t106;
t92 = t178 * t105;
t91 = -t165 * t209 - t166 * t202 + t167 * t203;
t90 = -t129 * t200 + t130 * t198;
t89 = -t165 * t200 - t166 * t178 + t167 * t179;
t88 = -t165 * t198 - t166 * t176 + t167 * t177;
t87 = (-t129 + t261) * t233 + t232 * t247;
t86 = t130 * t233 + t230 * t247 + t262;
t85 = t114 * t202 - t135 * t178;
t84 = -t113 * t202 + t135 * t176;
t83 = -t121 * t178 + t93;
t82 = -t105 * t202 + t117;
t81 = -t124 * t209 - t126 * t202 + t128 * t203;
t80 = -t123 * t209 - t125 * t202 + t127 * t203;
t79 = t132 * t202 + t133 * t174 + t134 * t175;
t78 = (t129 * t230 + t130 * t232) * t231 + t248;
t77 = -t124 * t200 - t126 * t178 + t128 * t179;
t76 = -t123 * t200 - t125 * t178 + t127 * t179;
t75 = -t124 * t198 - t126 * t176 + t128 * t177;
t74 = -t123 * t198 - t125 * t176 + t127 * t177;
t72 = t113 * t178 - t114 * t176;
t71 = -t106 * t176 + t92;
t70 = t132 * t178 + t133 * t150 + t134 * t151;
t69 = t132 * t176 + t133 * t148 + t134 * t149;
t68 = -t114 * t209 - t136 + (t135 + t171) * t200;
t67 = -t135 * t198 - t146 - (-t113 - t140) * t209;
t64 = (-t113 + t255) * t233 + t232 * t242;
t63 = t114 * t233 + t230 * t242 + t256;
t62 = -t113 * t200 - t131 + (t114 + t141) * t198;
t61 = t108 * t202 + t110 * t174 + t112 * t175;
t60 = t107 * t202 + t109 * t174 + t111 * t175;
t57 = t178 * t263 + t202 * t98 + t93;
t56 = t122 * t176 - t202 * t273 + t117;
t55 = t108 * t178 + t110 * t150 + t112 * t151;
t54 = t107 * t178 + t109 * t150 + t111 * t151;
t53 = t108 * t176 + t110 * t148 + t112 * t149;
t52 = t107 * t176 + t109 * t148 + t111 * t149;
t51 = (t113 * t230 + t114 * t232) * t231 + t241;
t46 = -t136 - t272 * t209 + (t171 - t263) * t200;
t45 = -t146 + t263 * t198 - (-t140 - t273) * t209;
t44 = (t255 - t273) * t233 + t232 * t240;
t43 = t230 * t240 + t233 * t272 + t256;
t42 = -t176 * t272 + t178 * t97 + t92;
t41 = t233 * t91 + (t230 * t81 - t232 * t80) * t231;
t40 = -t198 * t80 - t200 * t81 - t209 * t91;
t39 = -t131 - t273 * t200 + (t141 + t272) * t198;
t38 = (t230 * t273 + t232 * t272) * t231 + t241;
t37 = t233 * t89 + (t230 * t77 - t232 * t76) * t231;
t36 = t233 * t88 + (t230 * t75 - t232 * t74) * t231;
t35 = -t198 * t76 - t200 * t77 - t209 * t89;
t34 = -t198 * t74 - t200 * t75 - t209 * t88;
t33 = t233 * t79 + (t230 * t61 - t232 * t60) * t231;
t32 = -t198 * t60 - t200 * t61 - t209 * t79;
t30 = t176 * t60 + t178 * t61 + t202 * t79;
t22 = t233 * t70 + (t230 * t55 - t232 * t54) * t231;
t21 = t233 * t69 + (t230 * t53 - t232 * t52) * t231;
t20 = -t198 * t54 - t200 * t55 - t209 * t70;
t19 = -t198 * t52 - t200 * t53 - t209 * t69;
t16 = t176 * t54 + t178 * t55 + t202 * t70;
t15 = t176 * t52 + t178 * t53 + t202 * t69;
t1 = [m(3) + m(2) + t257; m(3) * t147 + m(4) * t96 + m(5) * t78 + m(6) * t51 + m(7) * t38; m(7) * (t38 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t51 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(5) * (t78 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(4) * (t115 ^ 2 + t116 ^ 2 + t96 ^ 2) + m(3) * (t147 ^ 2 + t169 ^ 2 + t170 ^ 2) + (t31 + t33 + t41 + (t209 * t183 + t210 * t184 + t288 * t233) * t233 + ((t209 * t155 + t210 * t157) * t230 - (t209 * t154 + t210 * t156) * t232 + (-t152 * t232 + t153 * t230 + t206 * t279 + t207 * t277) * t233) * t231) * t233 + (t18 + t22 + t37 + (t200 * t155 + t201 * t157 + t216 * t190 + t217 * t192 + t266 * t287) * t266 + ((t190 * t279 + t192 * t277) * t231 + t216 * t206 + t217 * t207 + t200 * t183 + t201 * t184 + t288 * t266 + t233 * t188) * t233) * t266 + (-t17 - t21 - t36 + (t198 * t154 + t199 * t156 + t214 * t189 + t215 * t191 + t265 * t289) * t265 + (-(t189 * t279 + t191 * t277) * t231 - t198 * t183 - t199 * t184 - t214 * t206 - t215 * t207 + t288 * t265 - t233 * t187) * t233 + (-t200 * t154 - t198 * t155 - t201 * t156 - t199 * t157 - t216 * t189 - t214 * t190 - t217 * t191 - t215 * t192 + t287 * t265 + t266 * t289) * t266) * t265; t257 * t233; m(7) * (t233 * t38 + (t230 * t44 - t232 * t43) * t231) + m(6) * (t233 * t51 + (t230 * t64 - t232 * t63) * t231) + m(5) * (t233 * t78 + (t230 * t87 - t232 * t86) * t231) + m(4) * (t233 * t96 + (-t115 * t232 + t116 * t230) * t231); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t233 ^ 2 + (t230 ^ 2 + t232 ^ 2) * t231 ^ 2); m(5) * t90 + m(6) * t62 + m(7) * t39; (t28 / 0.2e1 + t32 / 0.2e1 + t40 / 0.2e1) * t233 - (t31 / 0.2e1 + t33 / 0.2e1 + t41 / 0.2e1) * t209 + (-t18 / 0.2e1 - t22 / 0.2e1 - t37 / 0.2e1) * t200 + (-t17 / 0.2e1 - t21 / 0.2e1 - t36 / 0.2e1) * t198 + m(7) * (t38 * t39 + t43 * t46 + t44 * t45) + m(6) * (t51 * t62 + t63 * t68 + t64 * t67) + m(5) * (t78 * t90 + t86 * t95 + t87 * t94) + ((-t11 / 0.2e1 - t19 / 0.2e1 - t34 / 0.2e1) * t232 + (t12 / 0.2e1 + t20 / 0.2e1 + t35 / 0.2e1) * t230) * t231; m(5) * (t233 * t90 + (t230 * t94 - t232 * t95) * t231) + m(6) * (t233 * t62 + (t230 * t67 - t232 * t68) * t231) + m(7) * (t233 * t39 + (t230 * t45 - t232 * t46) * t231); -(t28 + t32 + t40) * t209 + (-t12 - t20 - t35) * t200 + (-t11 - t19 - t34) * t198 + m(7) * (t39 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t62 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(5) * (t90 ^ 2 + t94 ^ 2 + t95 ^ 2); m(6) * t72 + m(7) * t42; t22 * t285 + t33 * t282 + t30 * t280 + t21 * t286 + (t230 * t16 / 0.2e1 - t232 * t15 / 0.2e1) * t231 + m(7) * (t38 * t42 + t43 * t57 + t44 * t56) + m(6) * (t51 * t72 + t63 * t85 + t64 * t84) + t243; m(6) * (t233 * t72 + (t230 * t84 - t232 * t85) * t231) + m(7) * (t233 * t42 + (t230 * t56 - t232 * t57) * t231); t19 * t286 + t32 * t282 + t20 * t285 + t16 * t283 + t15 * t284 + t30 * t281 + m(7) * (t39 * t42 + t45 * t56 + t46 * t57) + m(6) * (t62 * t72 + t67 * t84 + t68 * t85) + t244; t176 * t15 + t178 * t16 + t202 * t30 + m(7) * (t42 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(6) * (t72 ^ 2 + t84 ^ 2 + t85 ^ 2) + t258; m(7) * t71; m(7) * (t38 * t71 + t43 * t83 + t44 * t82) + t243; m(7) * (t233 * t71 + (t230 * t82 - t232 * t83) * t231); m(7) * (t39 * t71 + t45 * t82 + t46 * t83) + t244; m(7) * (t42 * t71 + t56 * t82 + t57 * t83) + t258; m(7) * (t71 ^ 2 + t82 ^ 2 + t83 ^ 2) + t258;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
