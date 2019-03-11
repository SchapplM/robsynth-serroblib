% Calculate joint inertia matrix for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:16
% EndTime: 2019-03-09 05:00:23
% DurationCPUTime: 3.27s
% Computational Cost: add. (11236->434), mult. (9274->613), div. (0->0), fcn. (9935->12), ass. (0->225)
t208 = sin(qJ(3));
t291 = Icges(4,5) * t208;
t290 = t291 / 0.2e1;
t203 = qJ(4) + pkin(11);
t195 = sin(t203);
t197 = cos(t203);
t211 = cos(qJ(3));
t147 = -Icges(6,3) * t211 + (Icges(6,5) * t197 - Icges(6,6) * t195) * t208;
t207 = sin(qJ(4));
t210 = cos(qJ(4));
t159 = -Icges(5,3) * t211 + (Icges(5,5) * t210 - Icges(5,6) * t207) * t208;
t289 = -t147 - t159;
t148 = -Icges(6,6) * t211 + (Icges(6,4) * t197 - Icges(6,2) * t195) * t208;
t160 = -Icges(5,6) * t211 + (Icges(5,4) * t210 - Icges(5,2) * t207) * t208;
t288 = -t148 * t195 - t160 * t207;
t204 = qJ(1) + pkin(10);
t196 = sin(t204);
t198 = cos(t204);
t256 = t196 * t211;
t143 = -t195 * t256 - t197 * t198;
t144 = -t195 * t198 + t197 * t256;
t257 = t196 * t208;
t90 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t257;
t92 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t257;
t94 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t257;
t31 = t143 * t92 + t144 * t94 + t257 * t90;
t253 = t198 * t211;
t145 = -t195 * t253 + t196 * t197;
t146 = t195 * t196 + t197 * t253;
t254 = t198 * t208;
t91 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t254;
t93 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t254;
t95 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t254;
t32 = t143 * t93 + t144 * t95 + t257 * t91;
t250 = t207 * t211;
t155 = -t196 * t250 - t198 * t210;
t249 = t210 * t211;
t255 = t198 * t207;
t156 = t196 * t249 - t255;
t102 = Icges(5,5) * t156 + Icges(5,6) * t155 + Icges(5,3) * t257;
t104 = Icges(5,4) * t156 + Icges(5,2) * t155 + Icges(5,6) * t257;
t106 = Icges(5,1) * t156 + Icges(5,4) * t155 + Icges(5,5) * t257;
t39 = t102 * t257 + t104 * t155 + t106 * t156;
t157 = t196 * t210 - t198 * t250;
t258 = t196 * t207;
t158 = t198 * t249 + t258;
t103 = Icges(5,5) * t158 + Icges(5,6) * t157 + Icges(5,3) * t254;
t105 = Icges(5,4) * t158 + Icges(5,2) * t157 + Icges(5,6) * t254;
t107 = Icges(5,1) * t158 + Icges(5,4) * t157 + Icges(5,5) * t254;
t40 = t103 * t257 + t105 * t155 + t107 * t156;
t149 = -Icges(6,5) * t211 + (Icges(6,1) * t197 - Icges(6,4) * t195) * t208;
t57 = t143 * t148 + t144 * t149 + t147 * t257;
t161 = -Icges(5,5) * t211 + (Icges(5,1) * t210 - Icges(5,4) * t207) * t208;
t64 = t155 * t160 + t156 * t161 + t159 * t257;
t287 = (-t64 - t57) * t211 + ((t32 + t40) * t198 + (t31 + t39) * t196) * t208;
t33 = t145 * t92 + t146 * t94 + t254 * t90;
t34 = t145 * t93 + t146 * t95 + t254 * t91;
t41 = t102 * t254 + t104 * t157 + t106 * t158;
t42 = t103 * t254 + t105 * t157 + t107 * t158;
t58 = t145 * t148 + t146 * t149 + t147 * t254;
t65 = t157 * t160 + t158 * t161 + t159 * t254;
t286 = (-t65 - t58) * t211 + ((t34 + t42) * t198 + (t33 + t41) * t196) * t208;
t43 = -t211 * t90 + (-t195 * t92 + t197 * t94) * t208;
t47 = -t102 * t211 + (-t104 * t207 + t106 * t210) * t208;
t285 = -t43 - t47;
t44 = -t211 * t91 + (-t195 * t93 + t197 * t95) * t208;
t48 = -t103 * t211 + (-t105 * t207 + t107 * t210) * t208;
t284 = t44 + t48;
t192 = t210 * pkin(4) + pkin(3);
t167 = pkin(5) * t197 + t192;
t168 = pkin(4) * t207 + pkin(5) * t195;
t199 = qJ(6) + t203;
t190 = sin(t199);
t191 = cos(t199);
t125 = -t190 * t253 + t191 * t196;
t126 = t190 * t196 + t191 * t253;
t88 = t126 * rSges(7,1) + t125 * rSges(7,2) + rSges(7,3) * t254;
t283 = t167 * t253 + t196 * t168 + t88;
t282 = (t149 * t197 + t161 * t210) * t208;
t193 = t196 ^ 2;
t194 = t198 ^ 2;
t281 = m(6) / 0.2e1;
t280 = m(7) / 0.2e1;
t279 = -m(6) - m(7);
t123 = -t190 * t256 - t191 * t198;
t124 = -t190 * t198 + t191 * t256;
t81 = Icges(7,5) * t124 + Icges(7,6) * t123 + Icges(7,3) * t257;
t83 = Icges(7,4) * t124 + Icges(7,2) * t123 + Icges(7,6) * t257;
t85 = Icges(7,1) * t124 + Icges(7,4) * t123 + Icges(7,5) * t257;
t26 = t123 * t83 + t124 * t85 + t257 * t81;
t82 = Icges(7,5) * t126 + Icges(7,6) * t125 + Icges(7,3) * t254;
t84 = Icges(7,4) * t126 + Icges(7,2) * t125 + Icges(7,6) * t254;
t86 = Icges(7,1) * t126 + Icges(7,4) * t125 + Icges(7,5) * t254;
t27 = t123 * t84 + t124 * t86 + t257 * t82;
t127 = -Icges(7,3) * t211 + (Icges(7,5) * t191 - Icges(7,6) * t190) * t208;
t128 = -Icges(7,6) * t211 + (Icges(7,4) * t191 - Icges(7,2) * t190) * t208;
t129 = -Icges(7,5) * t211 + (Icges(7,1) * t191 - Icges(7,4) * t190) * t208;
t53 = t123 * t128 + t124 * t129 + t127 * t257;
t5 = -t211 * t53 + (t196 * t26 + t198 * t27) * t208;
t28 = t125 * t83 + t126 * t85 + t254 * t81;
t29 = t125 * t84 + t126 * t86 + t254 * t82;
t54 = t125 * t128 + t126 * t129 + t127 * t254;
t6 = -t211 * t54 + (t196 * t28 + t198 * t29) * t208;
t278 = t6 * t254 + t5 * t257;
t277 = t196 / 0.2e1;
t276 = -t198 / 0.2e1;
t275 = -t211 / 0.2e1;
t174 = rSges(4,1) * t208 + rSges(4,2) * t211;
t274 = m(4) * t174;
t209 = sin(qJ(1));
t273 = pkin(1) * t209;
t272 = pkin(3) * t211;
t271 = pkin(8) * t208;
t270 = -pkin(3) + t192;
t206 = -qJ(5) - pkin(8);
t202 = -pkin(9) + t206;
t238 = t202 - t206;
t243 = -pkin(4) * t258 - t192 * t253;
t269 = -t238 * t254 + t243 + t283;
t268 = t288 * t208 + t211 * t289 + t282;
t267 = t198 * rSges(4,3);
t117 = t208 * t191 * t129;
t262 = t128 * t190;
t68 = -t211 * t127 - t208 * t262 + t117;
t266 = t68 * t211;
t251 = t206 * t208;
t215 = -t198 * t251 - t243;
t240 = pkin(3) * t253 + pkin(8) * t254;
t109 = t215 - t240;
t97 = t146 * rSges(6,1) + t145 * rSges(6,2) + rSges(6,3) * t254;
t265 = -t109 - t97;
t136 = -rSges(7,3) * t211 + (rSges(7,1) * t191 - rSges(7,2) * t190) * t208;
t222 = -rSges(7,1) * t124 - rSges(7,2) * t123;
t87 = rSges(7,3) * t257 - t222;
t66 = t136 * t257 + t211 * t87;
t263 = Icges(4,4) * t211;
t259 = t168 * t198;
t252 = t202 * t208;
t241 = pkin(4) * t255 + t196 * t251;
t108 = (t211 * t270 - t271) * t196 - t241;
t142 = (pkin(8) + t206) * t211 + t270 * t208;
t248 = t211 * t108 + t142 * t257;
t247 = t193 * (t271 + t272) + t198 * t240;
t150 = -rSges(6,3) * t211 + (rSges(6,1) * t197 - rSges(6,2) * t195) * t208;
t246 = -t142 - t150;
t162 = -rSges(5,3) * t211 + (rSges(5,1) * t210 - rSges(5,2) * t207) * t208;
t181 = pkin(3) * t208 - pkin(8) * t211;
t244 = -t162 - t181;
t242 = t167 - t192;
t239 = t193 + t194;
t237 = -t109 - t269;
t116 = t208 * t242 + t211 * t238;
t236 = -t116 - t136 - t142;
t235 = -t181 + t246;
t111 = t158 * rSges(5,1) + t157 * rSges(5,2) + rSges(5,3) * t254;
t212 = cos(qJ(1));
t201 = t212 * pkin(1);
t234 = t198 * pkin(2) + t196 * pkin(7) + t201;
t233 = t257 / 0.2e1;
t232 = t254 / 0.2e1;
t231 = t198 * pkin(7) - t273;
t37 = -t211 * t81 + (-t190 * t83 + t191 * t85) * t208;
t38 = -t211 * t82 + (-t190 * t84 + t191 * t86) * t208;
t230 = (t37 + t53) * t233 + (t38 + t54) * t232;
t9 = -t266 + (t196 * t37 + t198 * t38) * t208;
t229 = -t211 * t9 + t278;
t228 = t196 * t108 + t198 * t109 + t247;
t227 = -t181 + t236;
t14 = t196 * t27 - t198 * t26;
t15 = t196 * t29 - t198 * t28;
t226 = t14 * t233 + t15 * t232 + t5 * t276 + t6 * t277 + (t38 * t196 - t37 * t198) * t275;
t225 = rSges(4,1) * t211 - rSges(4,2) * t208;
t224 = -rSges(5,1) * t156 - rSges(5,2) * t155;
t223 = -rSges(6,1) * t144 - rSges(6,2) * t143;
t220 = -Icges(4,2) * t208 + t263;
t219 = Icges(4,5) * t211 - Icges(4,6) * t208;
t216 = rSges(4,1) * t253 - rSges(4,2) * t254 + t196 * rSges(4,3);
t214 = t47 / 0.2e1 + t43 / 0.2e1 + t64 / 0.2e1 + t57 / 0.2e1;
t213 = t48 / 0.2e1 + t44 / 0.2e1 + t65 / 0.2e1 + t58 / 0.2e1;
t176 = rSges(2,1) * t212 - t209 * rSges(2,2);
t175 = -t209 * rSges(2,1) - rSges(2,2) * t212;
t171 = Icges(4,6) * t211 + t291;
t164 = rSges(3,1) * t198 - rSges(3,2) * t196 + t201;
t163 = -rSges(3,1) * t196 - rSges(3,2) * t198 - t273;
t131 = Icges(4,3) * t196 + t198 * t219;
t130 = -Icges(4,3) * t198 + t196 * t219;
t115 = t244 * t198;
t114 = t244 * t196;
t113 = t216 + t234;
t112 = t267 + (-pkin(2) - t225) * t196 + t231;
t110 = rSges(5,3) * t257 - t224;
t98 = t108 * t254;
t96 = rSges(6,3) * t257 - t223;
t89 = t198 * t216 + (t196 * t225 - t267) * t196;
t79 = t235 * t198;
t78 = t235 * t196;
t76 = t87 * t254;
t74 = -t259 + (t211 * t242 - t252) * t196 + t241;
t73 = -t111 * t211 - t162 * t254;
t72 = t110 * t211 + t162 * t257;
t71 = t234 + t111 + t240;
t70 = (-t272 - pkin(2) + (-rSges(5,3) - pkin(8)) * t208) * t196 + t224 + t231;
t67 = -t136 * t254 - t211 * t88;
t63 = t215 + t234 + t97;
t62 = (-rSges(6,3) * t208 - t192 * t211 - pkin(2)) * t196 + t223 + t231 + t241;
t61 = t227 * t198;
t60 = t227 * t196;
t59 = (t110 * t198 - t111 * t196) * t208;
t56 = -t198 * t252 + t234 + t283;
t55 = t259 + (-t167 * t211 - pkin(2) + (-rSges(7,3) + t202) * t208) * t196 + t222 + t231;
t52 = -t257 * t88 + t76;
t49 = t110 * t196 + t111 * t198 + t247;
t46 = t211 * t265 + t246 * t254;
t45 = t150 * t257 + t211 * t96 + t248;
t30 = t98 + (t196 * t265 + t198 * t96) * t208;
t25 = t196 * t96 + t198 * t97 + t228;
t24 = t211 * t237 + t236 * t254;
t23 = t116 * t257 + t211 * t74 + t248 + t66;
t22 = t196 * t42 - t198 * t41;
t21 = t196 * t40 - t198 * t39;
t20 = t76 + t98 + (t196 * t237 + t198 * t74) * t208;
t18 = t269 * t198 + (t74 + t87) * t196 + t228;
t17 = t196 * t34 - t198 * t33;
t16 = t196 * t32 - t198 * t31;
t1 = [Icges(2,3) + Icges(3,3) + t117 + (Icges(4,4) * t208 + Icges(4,2) * t211 - t127 + t289) * t211 + (Icges(4,1) * t208 - t262 + t263 + t288) * t208 + m(7) * (t55 ^ 2 + t56 ^ 2) + m(6) * (t62 ^ 2 + t63 ^ 2) + m(5) * (t70 ^ 2 + t71 ^ 2) + m(4) * (t112 ^ 2 + t113 ^ 2) + m(2) * (t175 ^ 2 + t176 ^ 2) + m(3) * (t163 ^ 2 + t164 ^ 2) + t282; 0; m(3) + m(4) + m(5) - t279; m(7) * (t55 * t61 + t56 * t60) + m(6) * (t62 * t79 + t63 * t78) + m(5) * (t114 * t71 + t115 * t70) + (-t37 / 0.2e1 - t53 / 0.2e1 - t112 * t274 + t196 * t220 * t275 - t214 + (-Icges(4,6) * t275 + t290 + t171 / 0.2e1) * t198) * t198 + (t54 / 0.2e1 + t38 / 0.2e1 + t211 * (Icges(4,6) * t196 + t198 * t220) / 0.2e1 + t196 * t290 - t113 * t274 + t171 * t277 + t213) * t196; m(4) * t89 + m(5) * t49 + m(6) * t25 + m(7) * t18; m(7) * (t18 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(6) * (t25 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(5) * (t114 ^ 2 + t115 ^ 2 + t49 ^ 2) + m(4) * (t174 ^ 2 * t239 + t89 ^ 2) + (-t194 * t130 - t14 - t16 - t21) * t198 + (t193 * t131 + t15 + t17 + t22 + (-t196 * t130 + t198 * t131) * t198) * t196; (-t68 - t268) * t211 + m(7) * (t23 * t55 + t24 * t56) + m(6) * (t45 * t62 + t46 * t63) + m(5) * (t70 * t72 + t71 * t73) + (t196 * t214 + t198 * t213) * t208 + t230; m(5) * t59 + m(6) * t30 + m(7) * t20; m(7) * (t18 * t20 + t23 * t61 + t24 * t60) + m(6) * (t30 * t25 + t45 * t79 + t46 * t78) + m(5) * (t114 * t73 + t115 * t72 + t59 * t49) + ((t22 / 0.2e1 + t17 / 0.2e1) * t198 + (t21 / 0.2e1 + t16 / 0.2e1) * t196) * t208 + t226 + t286 * t277 + t287 * t276 + (t284 * t196 + t285 * t198) * t275; m(7) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t30 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t59 ^ 2 + t72 ^ 2 + t73 ^ 2) + (t268 * t211 - t9) * t211 + ((-t284 * t211 + t286) * t198 + (t285 * t211 + t287) * t196) * t208 + t278; 0.2e1 * ((t196 * t56 + t198 * t55) * t280 + (t196 * t63 + t198 * t62) * t281) * t208; t279 * t211; m(7) * (-t18 * t211 + (t196 * t60 + t198 * t61) * t208) + m(6) * (-t211 * t25 + (t196 * t78 + t198 * t79) * t208); m(7) * (-t20 * t211 + (t196 * t24 + t198 * t23) * t208) + m(6) * (-t211 * t30 + (t196 * t46 + t198 * t45) * t208); 0.2e1 * (t281 + t280) * (t208 ^ 2 * t239 + t211 ^ 2); m(7) * (t55 * t66 + t56 * t67) - t266 + t230; m(7) * t52; m(7) * (t18 * t52 + t60 * t67 + t61 * t66) + t226; m(7) * (t20 * t52 + t23 * t66 + t24 * t67) + t229; m(7) * (-t211 * t52 + (t196 * t67 + t198 * t66) * t208); m(7) * (t52 ^ 2 + t66 ^ 2 + t67 ^ 2) + t229;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
