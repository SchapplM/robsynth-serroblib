% Calculate joint inertia matrix for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:07:55
% EndTime: 2019-03-09 10:08:02
% DurationCPUTime: 3.15s
% Computational Cost: add. (8420->413), mult. (6623->600), div. (0->0), fcn. (6884->12), ass. (0->199)
t189 = cos(pkin(11));
t171 = pkin(5) * t189 + pkin(4);
t191 = -pkin(9) - qJ(5);
t188 = sin(pkin(11));
t193 = sin(qJ(1));
t244 = t193 * t188;
t185 = qJ(2) + pkin(10);
t177 = qJ(4) + t185;
t170 = cos(t177);
t195 = cos(qJ(1));
t251 = t170 * t195;
t169 = sin(t177);
t252 = t169 * t195;
t184 = pkin(11) + qJ(6);
t175 = cos(t184);
t245 = t193 * t175;
t173 = sin(t184);
t250 = t173 * t195;
t124 = -t170 * t250 + t245;
t246 = t193 * t173;
t249 = t175 * t195;
t125 = t170 * t249 + t246;
t64 = t125 * rSges(7,1) + t124 * rSges(7,2) + rSges(7,3) * t252;
t287 = pkin(5) * t244 + t171 * t251 - t191 * t252 + t64;
t286 = Icges(3,3) + Icges(4,3);
t174 = sin(t185);
t176 = cos(t185);
t192 = sin(qJ(2));
t194 = cos(qJ(2));
t285 = Icges(3,5) * t194 + Icges(4,5) * t176 - Icges(3,6) * t192 - Icges(4,6) * t174;
t186 = t193 ^ 2;
t284 = t193 * pkin(7);
t225 = rSges(4,1) * t176 - rSges(4,2) * t174;
t247 = t189 * t195;
t134 = -t170 * t244 - t247;
t243 = t193 * t189;
t248 = t188 * t195;
t135 = t170 * t243 - t248;
t223 = -t135 * rSges(6,1) - t134 * rSges(6,2);
t136 = -t170 * t248 + t243;
t137 = t170 * t247 + t244;
t233 = t137 * rSges(6,1) + t136 * rSges(6,2) + rSges(6,3) * t252;
t253 = t169 * t193;
t283 = t193 * (rSges(6,3) * t253 - t223) + t195 * t233;
t242 = qJ(5) + t191;
t270 = -pkin(4) + t171;
t92 = -t170 * rSges(7,3) + (rSges(7,1) * t175 - rSges(7,2) * t173) * t169;
t282 = -t270 * t169 - t242 * t170 - t92;
t281 = -t285 * t193 + t286 * t195;
t280 = t286 * t193 + t285 * t195;
t213 = Icges(5,5) * t170 - Icges(5,6) * t169;
t103 = -Icges(5,3) * t195 + t213 * t193;
t104 = Icges(5,3) * t193 + t213 * t195;
t187 = t195 ^ 2;
t254 = Icges(5,4) * t170;
t216 = -Icges(5,2) * t169 + t254;
t106 = Icges(5,6) * t193 + t216 * t195;
t255 = Icges(5,4) * t169;
t219 = Icges(5,1) * t170 - t255;
t108 = Icges(5,5) * t193 + t219 * t195;
t211 = -t106 * t169 + t108 * t170;
t105 = -Icges(5,6) * t195 + t216 * t193;
t107 = -Icges(5,5) * t195 + t219 * t193;
t212 = t105 * t169 - t107 * t170;
t70 = Icges(6,5) * t135 + Icges(6,6) * t134 + Icges(6,3) * t253;
t71 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t252;
t72 = Icges(6,4) * t135 + Icges(6,2) * t134 + Icges(6,6) * t253;
t73 = Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t252;
t74 = Icges(6,1) * t135 + Icges(6,4) * t134 + Icges(6,5) * t253;
t75 = Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t252;
t122 = -t170 * t246 - t249;
t123 = t170 * t245 - t250;
t57 = Icges(7,5) * t123 + Icges(7,6) * t122 + Icges(7,3) * t253;
t59 = Icges(7,4) * t123 + Icges(7,2) * t122 + Icges(7,6) * t253;
t61 = Icges(7,1) * t123 + Icges(7,4) * t122 + Icges(7,5) * t253;
t16 = t122 * t59 + t123 * t61 + t57 * t253;
t58 = Icges(7,5) * t125 + Icges(7,6) * t124 + Icges(7,3) * t252;
t60 = Icges(7,4) * t125 + Icges(7,2) * t124 + Icges(7,6) * t252;
t62 = Icges(7,1) * t125 + Icges(7,4) * t124 + Icges(7,5) * t252;
t17 = t122 * t60 + t123 * t62 + t58 * t253;
t8 = -t16 * t195 + t17 * t193;
t279 = -t8 - t187 * t103 + (t134 * t72 + t135 * t74 + t70 * t253) * t195 + (-(-t104 + t212) * t195 - t134 * t73 - t135 * t75 - t71 * t253 - t211 * t193) * t193;
t240 = pkin(4) * t251 + qJ(5) * t252;
t222 = -t123 * rSges(7,1) - t122 * rSges(7,2);
t63 = rSges(7,3) * t253 - t222;
t278 = (-t240 + t287) * t195 + (-pkin(5) * t248 + t63 + (-t242 * t169 + t270 * t170) * t193) * t193;
t277 = t170 ^ 2;
t276 = m(6) / 0.2e1;
t275 = m(7) / 0.2e1;
t274 = t193 / 0.2e1;
t273 = -t195 / 0.2e1;
t272 = pkin(2) * t192;
t271 = pkin(4) * t170;
t190 = -qJ(3) - pkin(7);
t172 = t194 * pkin(2) + pkin(1);
t204 = rSges(5,1) * t251 - rSges(5,2) * t252 + t193 * rSges(5,3);
t224 = rSges(5,1) * t170 - rSges(5,2) * t169;
t67 = t193 * (-t195 * rSges(5,3) + t224 * t193) + t195 * t204;
t269 = rSges(3,1) * t194;
t267 = rSges(3,2) * t192;
t90 = -Icges(7,6) * t170 + (Icges(7,4) * t175 - Icges(7,2) * t173) * t169;
t265 = t173 * t90;
t264 = t195 * rSges(3,3);
t22 = -t170 * t57 + (-t173 * t59 + t175 * t61) * t169;
t263 = t22 * t195;
t23 = -t170 * t58 + (-t173 * t60 + t175 * t62) * t169;
t262 = t23 * t193;
t164 = t195 * t172;
t182 = t195 * pkin(7);
t261 = t195 * (-t195 * pkin(1) + t164 - t284) + t193 * (t182 + (-pkin(1) + t172) * t193);
t142 = t169 * pkin(4) - t170 * qJ(5);
t96 = -t170 * rSges(6,3) + (rSges(6,1) * t189 - rSges(6,2) * t188) * t169;
t260 = -t142 - t96;
t259 = Icges(3,4) * t192;
t258 = Icges(3,4) * t194;
t257 = Icges(4,4) * t174;
t256 = Icges(4,4) * t176;
t241 = t186 * (qJ(5) * t169 + t271) + t195 * t240;
t239 = t193 * rSges(3,3) + t195 * t269;
t237 = t186 + t187;
t18 = t124 * t59 + t125 * t61 + t57 * t252;
t19 = t124 * t60 + t125 * t62 + t58 * t252;
t9 = -t18 * t195 + t19 * t193;
t236 = (t9 + t186 * t104 + (t136 * t73 + t137 * t75 + t71 * t252) * t193 + (-t136 * t72 - t137 * t74 - t70 * t252 + (-t103 + t211) * t193 + t212 * t195) * t195) * t193;
t235 = t276 + t275;
t234 = -t142 + t282;
t232 = -rSges(4,1) * t174 - rSges(4,2) * t176 - t272;
t150 = pkin(3) * t176 + t172;
t145 = t195 * t150;
t231 = t195 * (t145 - t164) + t261 + (t150 - t172) * t186;
t183 = -pkin(8) + t190;
t230 = -t193 * t183 + t145;
t89 = -Icges(7,3) * t170 + (Icges(7,5) * t175 - Icges(7,6) * t173) * t169;
t91 = -Icges(7,5) * t170 + (Icges(7,1) * t175 - Icges(7,4) * t173) * t169;
t32 = t122 * t90 + t123 * t91 + t89 * t253;
t3 = -t32 * t170 + (t16 * t193 + t17 * t195) * t169;
t33 = t124 * t90 + t125 * t91 + t89 * t252;
t4 = -t33 * t170 + (t18 * t193 + t19 * t195) * t169;
t229 = t3 * t273 + t4 * t274 - t170 * (t262 - t263) / 0.2e1 + t8 * t253 / 0.2e1 + t9 * t252 / 0.2e1;
t227 = -pkin(3) * t174 - t272;
t226 = -t267 + t269;
t221 = Icges(3,1) * t194 - t259;
t220 = Icges(4,1) * t176 - t257;
t218 = -Icges(3,2) * t192 + t258;
t217 = -Icges(4,2) * t174 + t256;
t140 = Icges(5,2) * t170 + t255;
t141 = Icges(5,1) * t169 + t254;
t206 = -t140 * t169 + t141 * t170;
t205 = t193 * rSges(4,3) + t225 * t195;
t203 = t231 + t241;
t202 = -t142 + t227;
t143 = rSges(5,1) * t169 + rSges(5,2) * t170;
t201 = -t143 + t227;
t199 = t202 - t96;
t198 = t195 * t279 + t236;
t197 = t202 + t282;
t139 = Icges(5,5) * t169 + Icges(5,6) * t170;
t93 = -Icges(6,3) * t170 + (Icges(6,5) * t189 - Icges(6,6) * t188) * t169;
t94 = -Icges(6,6) * t170 + (Icges(6,4) * t189 - Icges(6,2) * t188) * t169;
t95 = -Icges(6,5) * t170 + (Icges(6,1) * t189 - Icges(6,4) * t188) * t169;
t196 = -t263 / 0.2e1 + t262 / 0.2e1 + (t136 * t94 + t137 * t95 + t193 * t139 + t206 * t195 + t93 * t252 + t33 + (-t71 + t106) * t170 + (-t188 * t73 + t189 * t75 + t108) * t169) * t274 + (t134 * t94 + t135 * t95 - t195 * t139 + t206 * t193 + t93 * t253 + t32 + (-t70 + t105) * t170 + (-t188 * t72 + t189 * t74 + t107) * t169) * t273;
t162 = rSges(2,1) * t195 - t193 * rSges(2,2);
t161 = -t193 * rSges(2,1) - rSges(2,2) * t195;
t160 = rSges(3,1) * t192 + rSges(3,2) * t194;
t115 = t232 * t195;
t114 = t232 * t193;
t102 = t284 + (pkin(1) - t267) * t195 + t239;
t101 = t264 + t182 + (-pkin(1) - t226) * t193;
t88 = t201 * t195;
t87 = t201 * t193;
t86 = -t193 * t190 + t164 + t205;
t85 = (rSges(4,3) - t190) * t195 + (-t172 - t225) * t193;
t81 = t195 * (-t195 * t267 + t239) + (t226 * t193 - t264) * t193;
t80 = t169 * t175 * t91;
t79 = t204 + t230;
t78 = (rSges(5,3) - t183) * t195 + (-t150 - t224) * t193;
t77 = t260 * t195;
t76 = t260 * t193;
t54 = t199 * t195;
t53 = t199 * t193;
t48 = t230 + t233 + t240;
t47 = -t195 * t183 + (-t271 - t150 + (-rSges(6,3) - qJ(5)) * t169) * t193 + t223;
t46 = t234 * t195;
t45 = t234 * t193;
t44 = t195 * t205 + (-t195 * rSges(4,3) + t225 * t193) * t193 + t261;
t43 = -t170 * t64 - t92 * t252;
t42 = t170 * t63 + t92 * t253;
t41 = t197 * t195;
t40 = t197 * t193;
t39 = t230 + t287;
t38 = (pkin(5) * t188 - t183) * t195 + (-t170 * t171 - t150 + (-rSges(7,3) + t191) * t169) * t193 + t222;
t37 = -t169 * t265 - t170 * t89 + t80;
t36 = (-t193 * t64 + t195 * t63) * t169;
t31 = t241 + t283;
t26 = t231 + t67;
t15 = t241 + t278;
t14 = t203 + t283;
t11 = t203 + t278;
t1 = [t176 * (Icges(4,2) * t176 + t257) + t174 * (Icges(4,1) * t174 + t256) + t194 * (Icges(3,2) * t194 + t259) + t192 * (Icges(3,1) * t192 + t258) + Icges(2,3) + t80 + (-t89 + t140 - t93) * t170 + (-t188 * t94 + t189 * t95 + t141 - t265) * t169 + m(7) * (t38 ^ 2 + t39 ^ 2) + m(6) * (t47 ^ 2 + t48 ^ 2) + m(5) * (t78 ^ 2 + t79 ^ 2) + m(4) * (t85 ^ 2 + t86 ^ 2) + m(3) * (t101 ^ 2 + t102 ^ 2) + m(2) * (t161 ^ 2 + t162 ^ 2); m(3) * (-t101 * t195 - t102 * t193) * t160 + m(7) * (t38 * t41 + t39 * t40) + m(6) * (t47 * t54 + t48 * t53) + m(5) * (t78 * t88 + t79 * t87) + m(4) * (t114 * t86 + t115 * t85) + t196 + (t176 * (Icges(4,6) * t193 + t217 * t195) + t174 * (Icges(4,5) * t193 + t220 * t195) + t194 * (Icges(3,6) * t193 + t218 * t195) + t192 * (Icges(3,5) * t193 + t221 * t195)) * t274 + (t176 * (-Icges(4,6) * t195 + t217 * t193) + t174 * (-Icges(4,5) * t195 + t220 * t193) + t194 * (-Icges(3,6) * t195 + t218 * t193) + t192 * (-Icges(3,5) * t195 + t221 * t193)) * t273 + (Icges(3,5) * t192 + Icges(4,5) * t174 + Icges(3,6) * t194 + Icges(4,6) * t176) * (t187 / 0.2e1 + t186 / 0.2e1); m(7) * (t11 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t14 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(5) * (t26 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t114 ^ 2 + t115 ^ 2 + t44 ^ 2) + m(3) * (t237 * t160 ^ 2 + t81 ^ 2) + t236 + t280 * t193 * t186 + (t281 * t187 + (t193 * t281 + t195 * t280) * t193 + t279) * t195; m(7) * (t193 * t38 - t195 * t39) + m(6) * (t193 * t47 - t195 * t48) + m(5) * (t193 * t78 - t195 * t79) + m(4) * (t193 * t85 - t195 * t86); m(7) * (t193 * t41 - t195 * t40) + m(6) * (t193 * t54 - t195 * t53) + m(5) * (t193 * t88 - t195 * t87) + m(4) * (-t114 * t195 + t193 * t115); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t235) * t237; m(7) * (t38 * t46 + t39 * t45) + m(6) * (t47 * t77 + t48 * t76) + m(5) * (-t193 * t79 - t195 * t78) * t143 + t196; m(7) * (t11 * t15 + t40 * t45 + t41 * t46) + m(6) * (t14 * t31 + t53 * t76 + t54 * t77) + m(5) * (t67 * t26 + (-t193 * t87 - t195 * t88) * t143) + t198; m(6) * (t77 * t193 - t195 * t76) + m(7) * (t46 * t193 - t195 * t45); m(7) * (t15 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t31 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t237 * t143 ^ 2 + t67 ^ 2) + t198; 0.2e1 * ((t193 * t39 + t195 * t38) * t275 + (t193 * t48 + t195 * t47) * t276) * t169; m(7) * (-t170 * t11 + (t193 * t40 + t195 * t41) * t169) + m(6) * (-t170 * t14 + (t193 * t53 + t195 * t54) * t169); 0; m(7) * (-t170 * t15 + (t193 * t45 + t195 * t46) * t169) + m(6) * (-t170 * t31 + (t193 * t76 + t195 * t77) * t169); 0.2e1 * t235 * (t169 ^ 2 * t237 + t277); -t37 * t170 + m(7) * (t38 * t42 + t39 * t43) + ((t33 / 0.2e1 + t23 / 0.2e1) * t195 + (t22 / 0.2e1 + t32 / 0.2e1) * t193) * t169; m(7) * (t11 * t36 + t40 * t43 + t41 * t42) + t229; m(7) * (t42 * t193 - t195 * t43); m(7) * (t15 * t36 + t42 * t46 + t43 * t45) + t229; m(7) * (-t36 * t170 + (t193 * t43 + t195 * t42) * t169); t277 * t37 + m(7) * (t36 ^ 2 + t42 ^ 2 + t43 ^ 2) + (t195 * t4 + t193 * t3 - t170 * (t193 * t22 + t195 * t23)) * t169;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
