% Calculate joint inertia matrix for
% S6RPRRRR5
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:08:13
% EndTime: 2019-03-09 07:08:21
% DurationCPUTime: 3.20s
% Computational Cost: add. (13569->422), mult. (10277->612), div. (0->0), fcn. (11063->12), ass. (0->214)
t210 = cos(qJ(5));
t191 = pkin(5) * t210 + pkin(4);
t212 = -pkin(10) - pkin(9);
t208 = sin(qJ(5));
t209 = sin(qJ(1));
t257 = t208 * t209;
t201 = pkin(11) + qJ(3);
t194 = qJ(4) + t201;
t189 = cos(t194);
t211 = cos(qJ(1));
t264 = t189 * t211;
t188 = sin(t194);
t265 = t188 * t211;
t204 = qJ(5) + qJ(6);
t196 = cos(t204);
t260 = t196 * t209;
t195 = sin(t204);
t261 = t195 * t211;
t156 = -t189 * t261 + t260;
t259 = t196 * t211;
t262 = t195 * t209;
t157 = t189 * t259 + t262;
t95 = t157 * rSges(7,1) + t156 * rSges(7,2) + rSges(7,3) * t265;
t290 = pkin(5) * t257 + t191 * t264 - t212 * t265 + t95;
t202 = t209 ^ 2;
t249 = pkin(4) * t264 + pkin(9) * t265;
t289 = -t249 + t290;
t225 = Icges(5,5) * t189 - Icges(5,6) * t188;
t138 = -Icges(5,3) * t211 + t209 * t225;
t139 = Icges(5,3) * t209 + t211 * t225;
t154 = -t189 * t262 - t259;
t155 = t189 * t260 - t261;
t266 = t188 * t209;
t88 = Icges(7,5) * t155 + Icges(7,6) * t154 + Icges(7,3) * t266;
t90 = Icges(7,4) * t155 + Icges(7,2) * t154 + Icges(7,6) * t266;
t92 = Icges(7,1) * t155 + Icges(7,4) * t154 + Icges(7,5) * t266;
t28 = t154 * t90 + t155 * t92 + t266 * t88;
t89 = Icges(7,5) * t157 + Icges(7,6) * t156 + Icges(7,3) * t265;
t91 = Icges(7,4) * t157 + Icges(7,2) * t156 + Icges(7,6) * t265;
t93 = Icges(7,1) * t157 + Icges(7,4) * t156 + Icges(7,5) * t265;
t29 = t154 * t91 + t155 * t93 + t266 * t89;
t15 = t209 * t29 - t211 * t28;
t203 = t211 ^ 2;
t254 = t210 * t211;
t160 = -t189 * t257 - t254;
t255 = t209 * t210;
t256 = t208 * t211;
t161 = t189 * t255 - t256;
t106 = Icges(6,5) * t161 + Icges(6,6) * t160 + Icges(6,3) * t266;
t108 = Icges(6,4) * t161 + Icges(6,2) * t160 + Icges(6,6) * t266;
t110 = Icges(6,1) * t161 + Icges(6,4) * t160 + Icges(6,5) * t266;
t42 = t106 * t266 + t108 * t160 + t110 * t161;
t162 = -t189 * t256 + t255;
t163 = t189 * t254 + t257;
t107 = Icges(6,5) * t163 + Icges(6,6) * t162 + Icges(6,3) * t265;
t109 = Icges(6,4) * t163 + Icges(6,2) * t162 + Icges(6,6) * t265;
t111 = Icges(6,1) * t163 + Icges(6,4) * t162 + Icges(6,5) * t265;
t43 = t107 * t266 + t109 * t160 + t111 * t161;
t21 = t209 * t43 - t211 * t42;
t267 = Icges(5,4) * t189;
t227 = -Icges(5,2) * t188 + t267;
t141 = Icges(5,6) * t209 + t211 * t227;
t268 = Icges(5,4) * t188;
t229 = Icges(5,1) * t189 - t268;
t143 = Icges(5,5) * t209 + t211 * t229;
t223 = -t141 * t188 + t143 * t189;
t140 = -Icges(5,6) * t211 + t209 * t227;
t142 = -Icges(5,5) * t211 + t209 * t229;
t224 = t140 * t188 - t142 * t189;
t288 = -t15 - t21 - t203 * t138 - (t223 * t209 + (-t139 + t224) * t211) * t209;
t124 = -Icges(7,3) * t189 + (Icges(7,5) * t196 - Icges(7,6) * t195) * t188;
t125 = -Icges(7,6) * t189 + (Icges(7,4) * t196 - Icges(7,2) * t195) * t188;
t126 = -Icges(7,5) * t189 + (Icges(7,1) * t196 - Icges(7,4) * t195) * t188;
t55 = t124 * t266 + t125 * t154 + t126 * t155;
t5 = -t189 * t55 + (t209 * t28 + t211 * t29) * t188;
t30 = t156 * t90 + t157 * t92 + t265 * t88;
t31 = t156 * t91 + t157 * t93 + t265 * t89;
t56 = t124 * t265 + t125 * t156 + t126 * t157;
t6 = -t189 * t56 + (t209 * t30 + t211 * t31) * t188;
t287 = t6 * t265 + t5 * t266;
t286 = -t189 / 0.2e1;
t285 = t209 / 0.2e1;
t284 = -t211 / 0.2e1;
t192 = sin(t201);
t283 = pkin(3) * t192;
t282 = pkin(4) * t189;
t281 = -pkin(4) + t191;
t207 = -pkin(7) - qJ(2);
t280 = pkin(9) + t212;
t206 = cos(pkin(11));
t190 = t206 * pkin(2) + pkin(1);
t193 = cos(t201);
t279 = rSges(4,1) * t193;
t278 = rSges(4,2) * t192;
t40 = -t189 * t88 + (-t195 * t90 + t196 * t92) * t188;
t277 = t40 * t211;
t41 = -t189 * t89 + (-t195 * t91 + t196 * t93) * t188;
t276 = t41 * t209;
t48 = -t106 * t189 + (-t108 * t208 + t110 * t210) * t188;
t275 = t48 * t211;
t49 = -t107 * t189 + (-t109 * t208 + t111 * t210) * t188;
t274 = t49 * t209;
t116 = t188 * t196 * t126;
t263 = t195 * t125;
t63 = -t189 * t124 - t188 * t263 + t116;
t273 = t63 * t189;
t272 = rSges(3,3) + qJ(2);
t127 = -rSges(7,3) * t189 + (rSges(7,1) * t196 - rSges(7,2) * t195) * t188;
t231 = -t155 * rSges(7,1) - t154 * rSges(7,2);
t94 = rSges(7,3) * t266 - t231;
t68 = t127 * t266 + t189 * t94;
t270 = Icges(4,4) * t192;
t269 = Icges(4,4) * t193;
t129 = -Icges(6,6) * t189 + (Icges(6,4) * t210 - Icges(6,2) * t208) * t188;
t258 = t208 * t129;
t177 = pkin(3) * t193 + t190;
t170 = t211 * t177;
t253 = t211 * (-t190 * t211 + t170) + (t177 - t190) * t202;
t123 = t188 * t281 + t189 * t280;
t252 = -t123 - t127;
t131 = -rSges(6,3) * t189 + (rSges(6,1) * t210 - rSges(6,2) * t208) * t188;
t169 = t188 * pkin(4) - t189 * pkin(9);
t251 = -t131 - t169;
t218 = rSges(5,1) * t264 - rSges(5,2) * t265 + t209 * rSges(5,3);
t233 = rSges(5,1) * t189 - rSges(5,2) * t188;
t96 = t209 * (-rSges(5,3) * t211 + t209 * t233) + t211 * t218;
t250 = t202 * (pkin(9) * t188 + t282) + t211 * t249;
t248 = t209 * rSges(4,3) + t211 * t279;
t246 = t202 + t203;
t16 = t209 * t31 - t211 * t30;
t44 = t106 * t265 + t108 * t162 + t110 * t163;
t45 = t107 * t265 + t109 * t162 + t111 * t163;
t22 = t209 * t45 - t211 * t44;
t245 = (t202 * t139 + t16 + t22 + (t224 * t211 + (-t138 + t223) * t209) * t211) * t209;
t244 = -t169 + t252;
t113 = t163 * rSges(6,1) + t162 * rSges(6,2) + rSges(6,3) * t265;
t243 = t266 / 0.2e1;
t242 = t265 / 0.2e1;
t168 = rSges(5,1) * t188 + rSges(5,2) * t189;
t241 = -t168 - t283;
t240 = -t169 - t283;
t239 = (t40 + t55) * t243 + (t41 + t56) * t242;
t200 = -pkin(8) + t207;
t238 = -t209 * t200 + t170;
t9 = -t273 + (t209 * t40 + t211 * t41) * t188;
t237 = -t189 * t9 + t287;
t232 = -rSges(6,1) * t161 - rSges(6,2) * t160;
t112 = rSges(6,3) * t266 - t232;
t54 = t209 * t112 + t211 * t113 + t250;
t236 = t15 * t243 + t16 * t242 + t5 * t284 + t6 * t285 + (t276 - t277) * t286;
t235 = -t131 + t240;
t234 = -t278 + t279;
t230 = Icges(4,1) * t193 - t270;
t228 = -Icges(4,2) * t192 + t269;
t226 = Icges(4,5) * t193 - Icges(4,6) * t192;
t166 = Icges(5,2) * t189 + t268;
t167 = Icges(5,1) * t188 + t267;
t220 = -t166 * t188 + t167 * t189;
t219 = t240 + t252;
t101 = -pkin(5) * t256 + (-t188 * t280 + t189 * t281) * t209;
t25 = t250 + t289 * t211 + (t101 + t94) * t209;
t205 = sin(pkin(11));
t217 = rSges(3,1) * t206 - rSges(3,2) * t205 + pkin(1);
t215 = t288 * t211 + t245;
t128 = -Icges(6,3) * t189 + (Icges(6,5) * t210 - Icges(6,6) * t208) * t188;
t130 = -Icges(6,5) * t189 + (Icges(6,1) * t210 - Icges(6,4) * t208) * t188;
t60 = t128 * t266 + t129 * t160 + t130 * t161;
t10 = -t189 * t60 + (t209 * t42 + t211 * t43) * t188;
t61 = t128 * t265 + t129 * t162 + t130 * t163;
t11 = -t189 * t61 + (t209 * t44 + t211 * t45) * t188;
t214 = t10 * t284 + t11 * t285 + t21 * t243 + t22 * t242 + (t274 - t275) * t286 + t236;
t165 = Icges(5,5) * t188 + Icges(5,6) * t189;
t213 = t276 / 0.2e1 - t277 / 0.2e1 + t274 / 0.2e1 - t275 / 0.2e1 + (t141 * t189 + t143 * t188 + t165 * t209 + t211 * t220 + t56 + t61) * t285 + (t140 * t189 + t142 * t188 - t165 * t211 + t209 * t220 + t55 + t60) * t284;
t184 = rSges(2,1) * t211 - rSges(2,2) * t209;
t183 = -rSges(2,1) * t209 - rSges(2,2) * t211;
t176 = rSges(4,1) * t192 + rSges(4,2) * t193;
t149 = Icges(4,3) * t209 + t211 * t226;
t148 = -Icges(4,3) * t211 + t209 * t226;
t137 = t209 * t272 + t211 * t217;
t136 = -t209 * t217 + t211 * t272;
t133 = t241 * t211;
t132 = t241 * t209;
t122 = -t207 * t209 + (t190 - t278) * t211 + t248;
t121 = (rSges(4,3) - t207) * t211 + (-t190 - t234) * t209;
t120 = t188 * t210 * t130;
t115 = t218 + t238;
t114 = (rSges(5,3) - t200) * t211 + (-t177 - t233) * t209;
t105 = t251 * t211;
t104 = t251 * t209;
t103 = t211 * (-t211 * t278 + t248) + (-t211 * rSges(4,3) + t209 * t234) * t209;
t87 = t235 * t211;
t86 = t235 * t209;
t82 = t94 * t265;
t77 = t238 + t113 + t249;
t76 = -t211 * t200 + (-t282 - t177 + (-rSges(6,3) - pkin(9)) * t188) * t209 + t232;
t75 = t244 * t211;
t74 = t244 * t209;
t73 = -t113 * t189 - t131 * t265;
t72 = t112 * t189 + t131 * t266;
t71 = t219 * t211;
t70 = t219 * t209;
t69 = -t127 * t265 - t189 * t95;
t67 = -t189 * t128 - t188 * t258 + t120;
t66 = t238 + t290;
t65 = (pkin(5) * t208 - t200) * t211 + (-t189 * t191 - t177 + (-rSges(7,3) + t212) * t188) * t209 + t231;
t64 = (t112 * t211 - t113 * t209) * t188;
t62 = -t266 * t95 + t82;
t57 = t96 + t253;
t35 = -t189 * t289 + t252 * t265;
t34 = t101 * t189 + t123 * t266 + t68;
t27 = t54 + t253;
t26 = t82 + (t101 * t211 - t209 * t289) * t188;
t24 = t25 + t253;
t1 = [Icges(3,2) * t206 ^ 2 + t193 * (Icges(4,2) * t193 + t270) + t192 * (Icges(4,1) * t192 + t269) + Icges(2,3) + t116 + t120 + (Icges(3,1) * t205 + 0.2e1 * Icges(3,4) * t206) * t205 + (-t124 - t128 + t166) * t189 + (t167 - t258 - t263) * t188 + m(7) * (t65 ^ 2 + t66 ^ 2) + m(6) * (t76 ^ 2 + t77 ^ 2) + m(5) * (t114 ^ 2 + t115 ^ 2) + m(4) * (t121 ^ 2 + t122 ^ 2) + m(3) * (t136 ^ 2 + t137 ^ 2) + m(2) * (t183 ^ 2 + t184 ^ 2); m(7) * (t209 * t65 - t211 * t66) + m(6) * (t209 * t76 - t211 * t77) + m(5) * (t114 * t209 - t115 * t211) + m(4) * (t121 * t209 - t122 * t211) + m(3) * (t136 * t209 - t137 * t211); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t246; ((Icges(4,6) * t209 + t211 * t228) * t193 + (Icges(4,5) * t209 + t211 * t230) * t192) * t285 + t213 + m(7) * (t65 * t71 + t66 * t70) + m(6) * (t76 * t87 + t77 * t86) + m(5) * (t114 * t133 + t115 * t132) + (t202 / 0.2e1 + t203 / 0.2e1) * (Icges(4,5) * t192 + Icges(4,6) * t193) + ((-Icges(4,6) * t211 + t209 * t228) * t193 + (-Icges(4,5) * t211 + t209 * t230) * t192) * t284 + m(4) * (-t121 * t211 - t122 * t209) * t176; m(5) * (-t132 * t211 + t133 * t209) + m(6) * (t209 * t87 - t211 * t86) + m(7) * (t209 * t71 - t211 * t70); m(7) * (t24 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(6) * (t27 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(5) * (t132 ^ 2 + t133 ^ 2 + t57 ^ 2) + m(4) * (t176 ^ 2 * t246 + t103 ^ 2) + t209 * t202 * t149 + t245 + (-t203 * t148 + (-t209 * t148 + t211 * t149) * t209 + t288) * t211; t213 + m(7) * (t65 * t75 + t66 * t74) + m(6) * (t104 * t77 + t105 * t76) + m(5) * (-t114 * t211 - t115 * t209) * t168; m(6) * (-t104 * t211 + t105 * t209) + m(7) * (t209 * t75 - t211 * t74); m(7) * (t25 * t24 + t70 * t74 + t71 * t75) + m(6) * (t104 * t86 + t105 * t87 + t54 * t27) + m(5) * (t57 * t96 + (-t132 * t209 - t133 * t211) * t168) + t215; m(7) * (t25 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(6) * (t104 ^ 2 + t105 ^ 2 + t54 ^ 2) + m(5) * (t168 ^ 2 * t246 + t96 ^ 2) + t215; (-t63 - t67) * t189 + m(7) * (t34 * t65 + t35 * t66) + m(6) * (t72 * t76 + t73 * t77) + ((t61 / 0.2e1 + t49 / 0.2e1) * t211 + (t60 / 0.2e1 + t48 / 0.2e1) * t209) * t188 + t239; m(6) * (t209 * t72 - t211 * t73) + m(7) * (t209 * t34 - t211 * t35); m(7) * (t24 * t26 + t34 * t71 + t35 * t70) + m(6) * (t64 * t27 + t72 * t87 + t73 * t86) + t214; m(7) * (t26 * t25 + t34 * t75 + t35 * t74) + m(6) * (t104 * t73 + t105 * t72 + t64 * t54) + t214; (t67 * t189 - t9) * t189 + m(7) * (t26 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t64 ^ 2 + t72 ^ 2 + t73 ^ 2) + (t211 * t11 + t209 * t10 - t189 * (t48 * t209 + t49 * t211)) * t188 + t287; m(7) * (t65 * t68 + t66 * t69) - t273 + t239; m(7) * (t209 * t68 - t211 * t69); m(7) * (t24 * t62 + t68 * t71 + t69 * t70) + t236; m(7) * (t62 * t25 + t68 * t75 + t69 * t74) + t236; m(7) * (t26 * t62 + t34 * t68 + t35 * t69) + t237; m(7) * (t62 ^ 2 + t68 ^ 2 + t69 ^ 2) + t237;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
