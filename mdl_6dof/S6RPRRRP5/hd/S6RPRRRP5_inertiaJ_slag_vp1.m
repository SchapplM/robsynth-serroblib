% Calculate joint inertia matrix for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:46
% EndTime: 2019-03-09 06:10:55
% DurationCPUTime: 3.72s
% Computational Cost: add. (9000->373), mult. (8193->552), div. (0->0), fcn. (8804->10), ass. (0->184)
t182 = pkin(10) + qJ(3);
t177 = qJ(4) + t182;
t172 = cos(t177);
t190 = cos(qJ(5));
t191 = cos(qJ(1));
t232 = t190 * t191;
t188 = sin(qJ(5));
t189 = sin(qJ(1));
t234 = t189 * t188;
t145 = t172 * t234 + t232;
t233 = t189 * t190;
t235 = t188 * t191;
t146 = t172 * t233 - t235;
t171 = sin(t177);
t239 = t171 * t189;
t82 = Icges(7,5) * t146 + Icges(7,6) * t239 + Icges(7,3) * t145;
t88 = Icges(6,4) * t146 - Icges(6,2) * t145 + Icges(6,6) * t239;
t282 = t82 - t88;
t147 = t172 * t235 - t233;
t148 = t172 * t232 + t234;
t237 = t171 * t191;
t83 = Icges(7,5) * t148 + Icges(7,6) * t237 + Icges(7,3) * t147;
t89 = Icges(6,4) * t148 - Icges(6,2) * t147 + Icges(6,6) * t237;
t281 = t83 - t89;
t84 = Icges(6,5) * t146 - Icges(6,6) * t145 + Icges(6,3) * t239;
t86 = Icges(7,4) * t146 + Icges(7,2) * t239 + Icges(7,6) * t145;
t280 = t84 + t86;
t85 = Icges(6,5) * t148 - Icges(6,6) * t147 + Icges(6,3) * t237;
t87 = Icges(7,4) * t148 + Icges(7,2) * t237 + Icges(7,6) * t147;
t279 = t85 + t87;
t90 = Icges(7,1) * t146 + Icges(7,4) * t239 + Icges(7,5) * t145;
t92 = Icges(6,1) * t146 - Icges(6,4) * t145 + Icges(6,5) * t239;
t278 = t90 + t92;
t91 = Icges(7,1) * t148 + Icges(7,4) * t237 + Icges(7,5) * t147;
t93 = Icges(6,1) * t148 - Icges(6,4) * t147 + Icges(6,5) * t237;
t277 = t91 + t93;
t264 = rSges(7,3) + qJ(6);
t267 = rSges(7,1) + pkin(5);
t276 = -t264 * t145 - t267 * t146;
t275 = t282 * t145 + t278 * t146 + t280 * t239;
t274 = t281 * t145 + t277 * t146 + t279 * t239;
t273 = t282 * t147 + t278 * t148 + t280 * t237;
t272 = t281 * t147 + t277 * t148 + t279 * t237;
t111 = -Icges(7,6) * t172 + (Icges(7,5) * t190 + Icges(7,3) * t188) * t171;
t113 = -Icges(7,2) * t172 + (Icges(7,4) * t190 + Icges(7,6) * t188) * t171;
t115 = -Icges(7,4) * t172 + (Icges(7,1) * t190 + Icges(7,5) * t188) * t171;
t51 = t111 * t145 + t113 * t239 + t115 * t146;
t112 = -Icges(6,3) * t172 + (Icges(6,5) * t190 - Icges(6,6) * t188) * t171;
t114 = -Icges(6,6) * t172 + (Icges(6,4) * t190 - Icges(6,2) * t188) * t171;
t116 = -Icges(6,5) * t172 + (Icges(6,1) * t190 - Icges(6,4) * t188) * t171;
t52 = t112 * t239 - t114 * t145 + t116 * t146;
t271 = -t51 - t52;
t53 = t147 * t111 + t113 * t237 + t148 * t115;
t54 = t112 * t237 - t147 * t114 + t148 * t116;
t270 = -t53 - t54;
t183 = t189 ^ 2;
t269 = t271 * t172 + (t275 * t189 + t274 * t191) * t171;
t268 = t270 * t172 + (t273 * t189 + t272 * t191) * t171;
t266 = t274 * t189 - t275 * t191;
t265 = t272 * t189 - t273 * t191;
t246 = rSges(7,2) * t239 - t276;
t263 = rSges(7,2) * t237 + t264 * t147 + t267 * t148;
t262 = -t112 - t113;
t240 = t171 * t188;
t261 = t111 * t240 + (t115 + t116) * t171 * t190;
t204 = Icges(5,5) * t172 - Icges(5,6) * t171;
t125 = -Icges(5,3) * t191 + t204 * t189;
t126 = Icges(5,3) * t189 + t204 * t191;
t184 = t191 ^ 2;
t241 = Icges(5,4) * t172;
t206 = -Icges(5,2) * t171 + t241;
t128 = Icges(5,6) * t189 + t206 * t191;
t242 = Icges(5,4) * t171;
t208 = Icges(5,1) * t172 - t242;
t130 = Icges(5,5) * t189 + t208 * t191;
t202 = -t128 * t171 + t130 * t172;
t127 = -Icges(5,6) * t191 + t206 * t189;
t129 = -Icges(5,5) * t191 + t208 * t189;
t203 = t127 * t171 - t129 * t172;
t260 = -t184 * t125 - (t202 * t189 + (-t126 + t203) * t191) * t189 - t266;
t258 = t189 / 0.2e1;
t257 = -t191 / 0.2e1;
t175 = sin(t182);
t256 = pkin(3) * t175;
t255 = pkin(4) * t172;
t187 = -pkin(7) - qJ(2);
t186 = cos(pkin(10));
t173 = t186 * pkin(2) + pkin(1);
t254 = -t114 * t240 + t262 * t172 + t261;
t176 = cos(t182);
t253 = rSges(4,1) * t176;
t252 = rSges(4,2) * t175;
t39 = -t172 * t86 + (t188 * t82 + t190 * t90) * t171;
t251 = t39 * t191;
t40 = -t172 * t87 + (t188 * t83 + t190 * t91) * t171;
t250 = t40 * t189;
t41 = -t172 * t84 + (-t188 * t88 + t190 * t92) * t171;
t249 = t41 * t191;
t42 = -t172 * t85 + (-t188 * t89 + t190 * t93) * t171;
t248 = t42 * t189;
t247 = rSges(3,3) + qJ(2);
t244 = Icges(4,4) * t175;
t243 = Icges(4,4) * t176;
t236 = t172 * t191;
t181 = -pkin(8) + t187;
t231 = t191 * t181;
t161 = pkin(3) * t176 + t173;
t155 = t191 * t161;
t230 = t191 * (-t191 * t173 + t155) + (t161 - t173) * t183;
t228 = -t172 * rSges(7,2) + (t264 * t188 + t267 * t190) * t171;
t118 = -t172 * rSges(6,3) + (rSges(6,1) * t190 - rSges(6,2) * t188) * t171;
t154 = pkin(4) * t171 - pkin(9) * t172;
t227 = -t118 - t154;
t197 = rSges(5,1) * t236 - rSges(5,2) * t237 + t189 * rSges(5,3);
t211 = rSges(5,1) * t172 - rSges(5,2) * t171;
t74 = t189 * (-t191 * rSges(5,3) + t211 * t189) + t191 * t197;
t225 = pkin(4) * t236 + pkin(9) * t237;
t226 = t183 * (pkin(9) * t171 + t255) + t191 * t225;
t224 = t189 * rSges(4,3) + t191 * t253;
t222 = t183 + t184;
t221 = (t183 * t126 + (t203 * t191 + (-t125 + t202) * t189) * t191 + t265) * t189;
t220 = -t154 - t228;
t97 = t148 * rSges(6,1) - t147 * rSges(6,2) + rSges(6,3) * t237;
t153 = rSges(5,1) * t171 + rSges(5,2) * t172;
t217 = -t153 - t256;
t216 = -t154 - t256;
t215 = -t161 - t255;
t214 = -t189 * t181 + t155;
t210 = -t146 * rSges(6,1) + t145 * rSges(6,2);
t95 = rSges(6,3) * t239 - t210;
t45 = t189 * t95 + t191 * t97 + t226;
t213 = -t118 + t216;
t212 = -t252 + t253;
t209 = Icges(4,1) * t176 - t244;
t207 = -Icges(4,2) * t175 + t243;
t205 = Icges(4,5) * t176 - Icges(4,6) * t175;
t151 = Icges(5,2) * t172 + t242;
t152 = Icges(5,1) * t171 + t241;
t199 = -t151 * t171 + t152 * t172;
t198 = t216 - t228;
t185 = sin(pkin(10));
t196 = rSges(3,1) * t186 - rSges(3,2) * t185 + pkin(1);
t195 = t214 + t225;
t22 = t246 * t189 + t263 * t191 + t226;
t194 = t260 * t191 + t221;
t193 = -(t250 - t251 + t248 - t249) * t172 / 0.2e1 + t268 * t258 + t269 * t257 + t266 * t239 / 0.2e1 + t265 * t237 / 0.2e1;
t150 = Icges(5,5) * t171 + Icges(5,6) * t172;
t192 = -t251 / 0.2e1 + t250 / 0.2e1 - t249 / 0.2e1 + t248 / 0.2e1 + (t128 * t172 + t130 * t171 + t189 * t150 + t199 * t191 - t270) * t258 + (t127 * t172 + t129 * t171 - t191 * t150 + t199 * t189 - t271) * t257;
t168 = rSges(2,1) * t191 - t189 * rSges(2,2);
t167 = -t189 * rSges(2,1) - rSges(2,2) * t191;
t160 = rSges(4,1) * t175 + rSges(4,2) * t176;
t134 = Icges(4,3) * t189 + t205 * t191;
t133 = -Icges(4,3) * t191 + t205 * t189;
t124 = t247 * t189 + t196 * t191;
t123 = -t196 * t189 + t247 * t191;
t120 = t217 * t191;
t119 = t217 * t189;
t110 = -t189 * t187 + (t173 - t252) * t191 + t224;
t109 = (rSges(4,3) - t187) * t191 + (-t173 - t212) * t189;
t99 = t197 + t214;
t98 = (rSges(5,3) - t181) * t191 + (-t161 - t211) * t189;
t81 = t227 * t191;
t80 = t227 * t189;
t79 = t191 * (-t191 * t252 + t224) + (-t191 * rSges(4,3) + t212 * t189) * t189;
t73 = t213 * t191;
t72 = t213 * t189;
t67 = t220 * t191;
t66 = t220 * t189;
t65 = t198 * t191;
t64 = t198 * t189;
t63 = t195 + t97;
t62 = -t231 + ((-rSges(6,3) - pkin(9)) * t171 + t215) * t189 + t210;
t61 = -t118 * t237 - t172 * t97;
t60 = t118 * t239 + t172 * t95;
t57 = (-t189 * t97 + t191 * t95) * t171;
t56 = t195 + t263;
t55 = -t231 + ((-rSges(7,2) - pkin(9)) * t171 + t215) * t189 + t276;
t46 = t74 + t230;
t44 = -t172 * t263 - t228 * t237;
t43 = t246 * t172 + t228 * t239;
t24 = (-t189 * t263 + t246 * t191) * t171;
t23 = t45 + t230;
t21 = t22 + t230;
t1 = [Icges(3,2) * t186 ^ 2 + t176 * (Icges(4,2) * t176 + t244) + t175 * (Icges(4,1) * t175 + t243) + Icges(2,3) + (Icges(3,1) * t185 + 0.2e1 * Icges(3,4) * t186) * t185 + (-t114 * t188 + t152) * t171 + (t151 + t262) * t172 + m(7) * (t55 ^ 2 + t56 ^ 2) + m(6) * (t62 ^ 2 + t63 ^ 2) + m(5) * (t98 ^ 2 + t99 ^ 2) + m(4) * (t109 ^ 2 + t110 ^ 2) + m(3) * (t123 ^ 2 + t124 ^ 2) + m(2) * (t167 ^ 2 + t168 ^ 2) + t261; m(7) * (t189 * t55 - t191 * t56) + m(6) * (t189 * t62 - t191 * t63) + m(5) * (t189 * t98 - t191 * t99) + m(4) * (t189 * t109 - t110 * t191) + m(3) * (t189 * t123 - t124 * t191); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t222; m(4) * (-t109 * t191 - t110 * t189) * t160 + m(7) * (t55 * t65 + t56 * t64) + m(6) * (t62 * t73 + t63 * t72) + m(5) * (t119 * t99 + t120 * t98) + (t184 / 0.2e1 + t183 / 0.2e1) * (Icges(4,5) * t175 + Icges(4,6) * t176) + t192 + (t176 * (Icges(4,6) * t189 + t207 * t191) + t175 * (Icges(4,5) * t189 + t209 * t191)) * t258 + (t176 * (-Icges(4,6) * t191 + t207 * t189) + t175 * (-Icges(4,5) * t191 + t209 * t189)) * t257; m(5) * (-t119 * t191 + t120 * t189) + m(6) * (t73 * t189 - t191 * t72) + m(7) * (t65 * t189 - t191 * t64); m(7) * (t21 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(6) * (t23 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(5) * (t119 ^ 2 + t120 ^ 2 + t46 ^ 2) + t189 * t183 * t134 + m(4) * (t222 * t160 ^ 2 + t79 ^ 2) + t221 + (-t184 * t133 + (-t189 * t133 + t191 * t134) * t189 + t260) * t191; m(7) * (t55 * t67 + t56 * t66) + m(6) * (t62 * t81 + t63 * t80) + m(5) * (-t189 * t99 - t191 * t98) * t153 + t192; m(6) * (t81 * t189 - t191 * t80) + m(7) * (t67 * t189 - t191 * t66); m(7) * (t21 * t22 + t64 * t66 + t65 * t67) + m(6) * (t23 * t45 + t72 * t80 + t73 * t81) + m(5) * (t74 * t46 + (-t119 * t189 - t120 * t191) * t153) + t194; m(7) * (t22 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(6) * (t45 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t222 * t153 ^ 2 + t74 ^ 2) + t194; -t254 * t172 + m(7) * (t43 * t55 + t44 * t56) + m(6) * (t60 * t62 + t61 * t63) + ((t53 / 0.2e1 + t42 / 0.2e1 + t40 / 0.2e1 + t54 / 0.2e1) * t191 + (t52 / 0.2e1 + t51 / 0.2e1 + t41 / 0.2e1 + t39 / 0.2e1) * t189) * t171; m(6) * (t60 * t189 - t191 * t61) + m(7) * (t43 * t189 - t191 * t44); m(7) * (t21 * t24 + t43 * t65 + t44 * t64) + m(6) * (t23 * t57 + t60 * t73 + t61 * t72) + t193; m(7) * (t22 * t24 + t43 * t67 + t44 * t66) + m(6) * (t45 * t57 + t60 * t81 + t61 * t80) + t193; m(7) * (t24 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t57 ^ 2 + t60 ^ 2 + t61 ^ 2) + t254 * t172 ^ 2 + (t268 * t191 + t269 * t189 + ((-t40 - t42) * t191 + (-t39 - t41) * t189) * t172) * t171; m(7) * (t145 * t56 + t147 * t55); m(7) * (-t145 * t191 + t147 * t189); m(7) * (t145 * t64 + t147 * t65 + t21 * t240); m(7) * (t145 * t66 + t147 * t67 + t22 * t240); m(7) * (t145 * t44 + t147 * t43 + t24 * t240); m(7) * (t171 ^ 2 * t188 ^ 2 + t145 ^ 2 + t147 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
