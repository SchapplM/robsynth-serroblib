% Calculate joint inertia matrix for
% S6RPRRPR6
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:26
% EndTime: 2019-03-09 05:15:35
% DurationCPUTime: 3.78s
% Computational Cost: add. (11205->451), mult. (9478->651), div. (0->0), fcn. (10121->12), ass. (0->227)
t201 = pkin(10) + qJ(3);
t193 = sin(t201);
t292 = Icges(4,5) * t193;
t291 = t292 / 0.2e1;
t202 = qJ(4) + pkin(11);
t194 = sin(t202);
t195 = cos(t201);
t196 = cos(t202);
t127 = -Icges(6,3) * t195 + (Icges(6,5) * t196 - Icges(6,6) * t194) * t193;
t209 = sin(qJ(4));
t211 = cos(qJ(4));
t133 = -Icges(5,3) * t195 + (Icges(5,5) * t211 - Icges(5,6) * t209) * t193;
t290 = -t127 - t133;
t128 = -Icges(6,6) * t195 + (Icges(6,4) * t196 - Icges(6,2) * t194) * t193;
t134 = -Icges(5,6) * t195 + (Icges(5,4) * t211 - Icges(5,2) * t209) * t193;
t289 = -t194 * t128 - t209 * t134;
t210 = sin(qJ(1));
t212 = cos(qJ(1));
t253 = t210 * t194;
t259 = t196 * t212;
t154 = -t195 * t253 - t259;
t252 = t210 * t196;
t155 = -t194 * t212 + t195 * t252;
t263 = t193 * t210;
t90 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t263;
t92 = Icges(6,4) * t155 + Icges(6,2) * t154 + Icges(6,6) * t263;
t94 = Icges(6,1) * t155 + Icges(6,4) * t154 + Icges(6,5) * t263;
t31 = t154 * t92 + t155 * t94 + t90 * t263;
t260 = t195 * t212;
t156 = -t194 * t260 + t252;
t157 = t195 * t259 + t253;
t262 = t193 * t212;
t91 = Icges(6,5) * t157 + Icges(6,6) * t156 + Icges(6,3) * t262;
t93 = Icges(6,4) * t157 + Icges(6,2) * t156 + Icges(6,6) * t262;
t95 = Icges(6,1) * t157 + Icges(6,4) * t156 + Icges(6,5) * t262;
t32 = t154 * t93 + t155 * t95 + t91 * t263;
t249 = t211 * t212;
t251 = t210 * t209;
t161 = -t195 * t251 - t249;
t250 = t210 * t211;
t256 = t209 * t212;
t162 = t195 * t250 - t256;
t104 = Icges(5,5) * t162 + Icges(5,6) * t161 + Icges(5,3) * t263;
t106 = Icges(5,4) * t162 + Icges(5,2) * t161 + Icges(5,6) * t263;
t108 = Icges(5,1) * t162 + Icges(5,4) * t161 + Icges(5,5) * t263;
t43 = t104 * t263 + t106 * t161 + t108 * t162;
t163 = -t195 * t256 + t250;
t164 = t195 * t249 + t251;
t105 = Icges(5,5) * t164 + Icges(5,6) * t163 + Icges(5,3) * t262;
t107 = Icges(5,4) * t164 + Icges(5,2) * t163 + Icges(5,6) * t262;
t109 = Icges(5,1) * t164 + Icges(5,4) * t163 + Icges(5,5) * t262;
t44 = t105 * t263 + t107 * t161 + t109 * t162;
t129 = -Icges(6,5) * t195 + (Icges(6,1) * t196 - Icges(6,4) * t194) * t193;
t54 = t127 * t263 + t128 * t154 + t129 * t155;
t135 = -Icges(5,5) * t195 + (Icges(5,1) * t211 - Icges(5,4) * t209) * t193;
t62 = t133 * t263 + t134 * t161 + t135 * t162;
t288 = (-t62 - t54) * t195 + ((t32 + t44) * t212 + (t31 + t43) * t210) * t193;
t33 = t156 * t92 + t157 * t94 + t90 * t262;
t34 = t156 * t93 + t157 * t95 + t91 * t262;
t45 = t104 * t262 + t163 * t106 + t164 * t108;
t46 = t105 * t262 + t163 * t107 + t164 * t109;
t55 = t127 * t262 + t156 * t128 + t157 * t129;
t63 = t133 * t262 + t163 * t134 + t164 * t135;
t287 = (-t63 - t55) * t195 + ((t34 + t46) * t212 + (t33 + t45) * t210) * t193;
t39 = -t195 * t90 + (-t194 * t92 + t196 * t94) * t193;
t47 = -t104 * t195 + (-t106 * t209 + t108 * t211) * t193;
t286 = -t39 - t47;
t40 = -t195 * t91 + (-t194 * t93 + t196 * t95) * t193;
t48 = -t105 * t195 + (-t107 * t209 + t109 * t211) * t193;
t285 = t40 + t48;
t191 = t211 * pkin(4) + pkin(3);
t173 = pkin(5) * t196 + t191;
t174 = pkin(4) * t209 + pkin(5) * t194;
t197 = qJ(6) + t202;
t188 = sin(t197);
t189 = cos(t197);
t254 = t210 * t189;
t141 = -t188 * t260 + t254;
t255 = t210 * t188;
t142 = t189 * t260 + t255;
t88 = t142 * rSges(7,1) + t141 * rSges(7,2) + rSges(7,3) * t262;
t284 = t173 * t260 + t210 * t174 + t88;
t283 = (t129 * t196 + t135 * t211) * t193;
t203 = t210 ^ 2;
t204 = t212 ^ 2;
t282 = m(6) / 0.2e1;
t281 = m(7) / 0.2e1;
t139 = -t189 * t212 - t195 * t255;
t140 = -t188 * t212 + t195 * t254;
t81 = Icges(7,5) * t140 + Icges(7,6) * t139 + Icges(7,3) * t263;
t83 = Icges(7,4) * t140 + Icges(7,2) * t139 + Icges(7,6) * t263;
t85 = Icges(7,1) * t140 + Icges(7,4) * t139 + Icges(7,5) * t263;
t27 = t139 * t83 + t140 * t85 + t81 * t263;
t82 = Icges(7,5) * t142 + Icges(7,6) * t141 + Icges(7,3) * t262;
t84 = Icges(7,4) * t142 + Icges(7,2) * t141 + Icges(7,6) * t262;
t86 = Icges(7,1) * t142 + Icges(7,4) * t141 + Icges(7,5) * t262;
t28 = t139 * t84 + t140 * t86 + t82 * t263;
t123 = -Icges(7,3) * t195 + (Icges(7,5) * t189 - Icges(7,6) * t188) * t193;
t124 = -Icges(7,6) * t195 + (Icges(7,4) * t189 - Icges(7,2) * t188) * t193;
t125 = -Icges(7,5) * t195 + (Icges(7,1) * t189 - Icges(7,4) * t188) * t193;
t51 = t123 * t263 + t124 * t139 + t125 * t140;
t5 = -t51 * t195 + (t210 * t27 + t212 * t28) * t193;
t29 = t141 * t83 + t142 * t85 + t81 * t262;
t30 = t141 * t84 + t142 * t86 + t82 * t262;
t52 = t123 * t262 + t141 * t124 + t142 * t125;
t6 = -t52 * t195 + (t210 * t29 + t212 * t30) * t193;
t280 = t6 * t262 + t5 * t263;
t279 = -t195 / 0.2e1;
t278 = t210 / 0.2e1;
t277 = -t212 / 0.2e1;
t170 = rSges(4,1) * t193 + rSges(4,2) * t195;
t276 = m(4) * t170;
t275 = pkin(3) * t195;
t274 = pkin(8) * t193;
t273 = -pkin(3) + t191;
t207 = -qJ(5) - pkin(8);
t272 = t289 * t193 + t290 * t195 + t283;
t200 = -pkin(9) + t207;
t240 = t200 - t207;
t244 = -pkin(4) * t251 - t191 * t260;
t271 = -t240 * t262 + t244 + t284;
t115 = t193 * t189 * t125;
t264 = t188 * t124;
t59 = -t195 * t123 - t193 * t264 + t115;
t270 = t59 * t195;
t269 = rSges(3,3) + qJ(2);
t215 = -t207 * t262 - t244;
t241 = pkin(3) * t260 + pkin(8) * t262;
t103 = t215 - t241;
t97 = t157 * rSges(6,1) + t156 * rSges(6,2) + rSges(6,3) * t262;
t268 = -t103 - t97;
t242 = pkin(4) * t256 + t207 * t263;
t102 = (t273 * t195 - t274) * t210 - t242;
t122 = (pkin(8) + t207) * t195 + t273 * t193;
t267 = t195 * t102 + t122 * t263;
t126 = -rSges(7,3) * t195 + (rSges(7,1) * t189 - rSges(7,2) * t188) * t193;
t223 = -t140 * rSges(7,1) - t139 * rSges(7,2);
t87 = rSges(7,3) * t263 - t223;
t66 = t126 * t263 + t195 * t87;
t265 = Icges(4,4) * t195;
t208 = -pkin(7) - qJ(2);
t258 = t208 * t212;
t130 = -rSges(6,3) * t195 + (rSges(6,1) * t196 - rSges(6,2) * t194) * t193;
t248 = -t122 - t130;
t138 = -rSges(5,3) * t195 + (rSges(5,1) * t211 - rSges(5,2) * t209) * t193;
t172 = pkin(3) * t193 - pkin(8) * t195;
t247 = -t138 - t172;
t246 = t203 * (t274 + t275) + t212 * t241;
t243 = t173 - t191;
t239 = t203 + t204;
t238 = t282 + t281;
t237 = -t103 - t271;
t112 = t243 * t193 + t240 * t195;
t236 = -t112 - t122 - t126;
t235 = -t172 + t248;
t111 = t164 * rSges(5,1) + t163 * rSges(5,2) + rSges(5,3) * t262;
t234 = t263 / 0.2e1;
t233 = t262 / 0.2e1;
t37 = -t195 * t81 + (-t188 * t83 + t189 * t85) * t193;
t38 = -t195 * t82 + (-t188 * t84 + t189 * t86) * t193;
t232 = (t37 + t51) * t234 + (t38 + t52) * t233;
t206 = cos(pkin(10));
t190 = pkin(2) * t206 + pkin(1);
t231 = t212 * t190 - t210 * t208;
t9 = -t270 + (t210 * t37 + t212 * t38) * t193;
t230 = -t195 * t9 + t280;
t229 = t210 * t102 + t212 * t103 + t246;
t228 = -t172 + t236;
t14 = t28 * t210 - t212 * t27;
t15 = t30 * t210 - t212 * t29;
t227 = t14 * t234 + t15 * t233 + t5 * t277 + t6 * t278 + (t38 * t210 - t37 * t212) * t279;
t226 = rSges(4,1) * t195 - rSges(4,2) * t193;
t225 = -t162 * rSges(5,1) - t161 * rSges(5,2);
t224 = -t155 * rSges(6,1) - t154 * rSges(6,2);
t221 = -Icges(4,2) * t193 + t265;
t220 = Icges(4,5) * t195 - Icges(4,6) * t193;
t217 = rSges(4,1) * t260 - rSges(4,2) * t262 + t210 * rSges(4,3);
t205 = sin(pkin(10));
t216 = rSges(3,1) * t206 - rSges(3,2) * t205 + pkin(1);
t214 = t39 / 0.2e1 + t62 / 0.2e1 + t54 / 0.2e1 + t47 / 0.2e1;
t213 = t48 / 0.2e1 + t40 / 0.2e1 + t63 / 0.2e1 + t55 / 0.2e1;
t177 = rSges(2,1) * t212 - t210 * rSges(2,2);
t176 = -t210 * rSges(2,1) - rSges(2,2) * t212;
t167 = Icges(4,6) * t195 + t292;
t144 = Icges(4,3) * t210 + t220 * t212;
t143 = -Icges(4,3) * t212 + t220 * t210;
t132 = t269 * t210 + t216 * t212;
t131 = -t216 * t210 + t269 * t212;
t120 = t217 + t231;
t119 = (rSges(4,3) - t208) * t212 + (-t190 - t226) * t210;
t114 = t247 * t212;
t113 = t247 * t210;
t110 = rSges(5,3) * t263 - t225;
t99 = t212 * t217 + (-t212 * rSges(4,3) + t226 * t210) * t210;
t96 = rSges(6,3) * t263 - t224;
t89 = t102 * t262;
t79 = t87 * t262;
t77 = -t174 * t212 + (-t193 * t200 + t243 * t195) * t210 + t242;
t76 = t231 + t111 + t241;
t75 = -t258 + (-t275 - t190 + (-rSges(5,3) - pkin(8)) * t193) * t210 + t225;
t74 = t235 * t212;
t73 = t235 * t210;
t72 = -t195 * t111 - t138 * t262;
t71 = t110 * t195 + t138 * t263;
t69 = t215 + t231 + t97;
t68 = -t258 + (-rSges(6,3) * t193 - t191 * t195 - t190) * t210 + t224 + t242;
t67 = -t126 * t262 - t195 * t88;
t65 = (t110 * t212 - t111 * t210) * t193;
t61 = -t200 * t262 + t231 + t284;
t60 = (t174 - t208) * t212 + (-t173 * t195 - t190 + (-rSges(7,3) + t200) * t193) * t210 + t223;
t58 = t228 * t212;
t57 = t228 * t210;
t56 = -t88 * t263 + t79;
t53 = t210 * t110 + t111 * t212 + t246;
t42 = t268 * t195 + t248 * t262;
t41 = t130 * t263 + t195 * t96 + t267;
t26 = t89 + (t268 * t210 + t212 * t96) * t193;
t25 = t210 * t96 + t212 * t97 + t229;
t24 = t237 * t195 + t236 * t262;
t23 = t112 * t263 + t195 * t77 + t267 + t66;
t22 = t46 * t210 - t212 * t45;
t21 = t44 * t210 - t212 * t43;
t20 = t79 + t89 + (t237 * t210 + t212 * t77) * t193;
t19 = t271 * t212 + (t77 + t87) * t210 + t229;
t17 = t34 * t210 - t212 * t33;
t16 = t32 * t210 - t212 * t31;
t1 = [Icges(3,2) * t206 ^ 2 + Icges(2,3) + t115 + (Icges(3,1) * t205 + 0.2e1 * Icges(3,4) * t206) * t205 + (Icges(4,4) * t193 + Icges(4,2) * t195 - t123 + t290) * t195 + (Icges(4,1) * t193 - t264 + t265 + t289) * t193 + m(7) * (t60 ^ 2 + t61 ^ 2) + m(5) * (t75 ^ 2 + t76 ^ 2) + m(6) * (t68 ^ 2 + t69 ^ 2) + m(4) * (t119 ^ 2 + t120 ^ 2) + m(3) * (t131 ^ 2 + t132 ^ 2) + m(2) * (t176 ^ 2 + t177 ^ 2) + t283; m(7) * (t210 * t60 - t212 * t61) + m(5) * (t210 * t75 - t212 * t76) + m(6) * (t210 * t68 - t212 * t69) + m(4) * (t210 * t119 - t120 * t212) + m(3) * (t210 * t131 - t132 * t212); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t238) * t239; m(7) * (t57 * t61 + t58 * t60) + m(5) * (t113 * t76 + t114 * t75) + m(6) * (t68 * t74 + t69 * t73) + (-t37 / 0.2e1 - t51 / 0.2e1 + t221 * t210 * t279 - t119 * t276 - t214 + (-Icges(4,6) * t279 + t291 + t167 / 0.2e1) * t212) * t212 + (t52 / 0.2e1 + t38 / 0.2e1 + (Icges(4,6) * t210 + t221 * t212) * t195 / 0.2e1 + t210 * t291 - t120 * t276 + t167 * t278 + t213) * t210; m(5) * (-t113 * t212 + t114 * t210) + m(6) * (t74 * t210 - t212 * t73) + m(7) * (t58 * t210 - t212 * t57); m(7) * (t19 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(6) * (t25 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(5) * (t113 ^ 2 + t114 ^ 2 + t53 ^ 2) + m(4) * (t239 * t170 ^ 2 + t99 ^ 2) + (-t204 * t143 - t14 - t16 - t21) * t212 + (t203 * t144 + t15 + t17 + t22 + (-t210 * t143 + t212 * t144) * t212) * t210; (-t59 - t272) * t195 + m(7) * (t23 * t60 + t24 * t61) + m(5) * (t71 * t75 + t72 * t76) + m(6) * (t41 * t68 + t42 * t69) + (t210 * t214 + t212 * t213) * t193 + t232; m(5) * (t71 * t210 - t212 * t72) + m(6) * (t41 * t210 - t212 * t42) + m(7) * (t23 * t210 - t212 * t24); m(7) * (t19 * t20 + t23 * t58 + t24 * t57) + m(6) * (t26 * t25 + t41 * t74 + t42 * t73) + m(5) * (t113 * t72 + t114 * t71 + t65 * t53) + ((t22 / 0.2e1 + t17 / 0.2e1) * t212 + (t21 / 0.2e1 + t16 / 0.2e1) * t210) * t193 + t227 + (t285 * t210 + t286 * t212) * t279 + t287 * t278 + t288 * t277; m(7) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t26 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(5) * (t65 ^ 2 + t71 ^ 2 + t72 ^ 2) + (t272 * t195 - t9) * t195 + (t287 * t212 + t288 * t210 + (t286 * t210 - t285 * t212) * t195) * t193 + t280; 0.2e1 * ((t210 * t61 + t212 * t60) * t281 + (t210 * t69 + t212 * t68) * t282) * t193; 0; m(7) * (-t195 * t19 + (t210 * t57 + t212 * t58) * t193) + m(6) * (-t195 * t25 + (t210 * t73 + t212 * t74) * t193); m(7) * (-t195 * t20 + (t210 * t24 + t212 * t23) * t193) + m(6) * (-t195 * t26 + (t210 * t42 + t212 * t41) * t193); 0.2e1 * t238 * (t193 ^ 2 * t239 + t195 ^ 2); -t270 + m(7) * (t60 * t66 + t61 * t67) + t232; m(7) * (t66 * t210 - t212 * t67); m(7) * (t19 * t56 + t57 * t67 + t58 * t66) + t227; m(7) * (t20 * t56 + t23 * t66 + t24 * t67) + t230; m(7) * (-t56 * t195 + (t210 * t67 + t212 * t66) * t193); m(7) * (t56 ^ 2 + t66 ^ 2 + t67 ^ 2) + t230;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
