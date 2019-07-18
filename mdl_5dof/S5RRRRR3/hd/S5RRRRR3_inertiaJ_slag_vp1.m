% Calculate joint inertia matrix for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:11
% EndTime: 2019-07-18 17:17:19
% DurationCPUTime: 2.49s
% Computational Cost: add. (8850->363), mult. (9411->535), div. (0->0), fcn. (10237->10), ass. (0->200)
t199 = cos(qJ(4));
t185 = pkin(3) * t199 + pkin(2);
t196 = sin(qJ(4));
t198 = sin(qJ(1));
t241 = t198 * t196;
t195 = qJ(2) + qJ(3);
t189 = cos(t195);
t201 = cos(qJ(1));
t245 = t189 * t201;
t194 = qJ(4) + qJ(5);
t186 = sin(t194);
t188 = cos(t194);
t242 = t198 * t188;
t143 = -t186 * t245 + t242;
t243 = t198 * t186;
t144 = t188 * t245 + t243;
t187 = sin(t195);
t246 = t187 * t201;
t93 = t144 * rSges(6,1) + t143 * rSges(6,2) + rSges(6,3) * t246;
t274 = pkin(3) * t241 + t185 * t245 + t93;
t179 = pkin(2) * t245;
t273 = -t179 + t274;
t192 = t198 ^ 2;
t193 = t201 ^ 2;
t232 = t192 + t193;
t197 = sin(qJ(2));
t200 = cos(qJ(2));
t272 = rSges(3,1) * t200 - rSges(3,2) * t197;
t212 = Icges(4,5) * t189 - Icges(4,6) * t187;
t131 = -Icges(4,3) * t201 + t198 * t212;
t132 = Icges(4,3) * t198 + t201 * t212;
t141 = -t188 * t201 - t189 * t243;
t142 = -t186 * t201 + t189 * t242;
t247 = t187 * t198;
t86 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t247;
t88 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t247;
t90 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t247;
t26 = t141 * t88 + t142 * t90 + t247 * t86;
t87 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t246;
t89 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t246;
t91 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t246;
t27 = t141 * t89 + t142 * t91 + t247 * t87;
t15 = t27 * t198 - t201 * t26;
t239 = t199 * t201;
t157 = -t189 * t241 - t239;
t240 = t198 * t199;
t244 = t196 * t201;
t158 = t189 * t240 - t244;
t101 = Icges(5,4) * t158 + Icges(5,2) * t157 + Icges(5,6) * t247;
t103 = Icges(5,1) * t158 + Icges(5,4) * t157 + Icges(5,5) * t247;
t99 = Icges(5,5) * t158 + Icges(5,6) * t157 + Icges(5,3) * t247;
t38 = t101 * t157 + t103 * t158 + t247 * t99;
t159 = -t189 * t244 + t240;
t160 = t189 * t239 + t241;
t100 = Icges(5,5) * t160 + Icges(5,6) * t159 + Icges(5,3) * t246;
t102 = Icges(5,4) * t160 + Icges(5,2) * t159 + Icges(5,6) * t246;
t104 = Icges(5,1) * t160 + Icges(5,4) * t159 + Icges(5,5) * t246;
t39 = t100 * t247 + t102 * t157 + t104 * t158;
t21 = t39 * t198 - t201 * t38;
t250 = Icges(4,4) * t189;
t214 = -Icges(4,2) * t187 + t250;
t134 = Icges(4,6) * t198 + t201 * t214;
t251 = Icges(4,4) * t187;
t216 = Icges(4,1) * t189 - t251;
t136 = Icges(4,5) * t198 + t201 * t216;
t210 = -t134 * t187 + t136 * t189;
t133 = -Icges(4,6) * t201 + t198 * t214;
t135 = -Icges(4,5) * t201 + t198 * t216;
t211 = t133 * t187 - t135 * t189;
t271 = -t15 - t21 - t193 * t131 - (t210 * t198 + (-t132 + t211) * t201) * t198;
t119 = -Icges(6,3) * t189 + (Icges(6,5) * t188 - Icges(6,6) * t186) * t187;
t120 = -Icges(6,6) * t189 + (Icges(6,4) * t188 - Icges(6,2) * t186) * t187;
t121 = -Icges(6,5) * t189 + (Icges(6,1) * t188 - Icges(6,4) * t186) * t187;
t55 = t119 * t247 + t120 * t141 + t121 * t142;
t5 = -t55 * t189 + (t198 * t26 + t201 * t27) * t187;
t28 = t143 * t88 + t144 * t90 + t246 * t86;
t29 = t143 * t89 + t144 * t91 + t246 * t87;
t56 = t119 * t246 + t143 * t120 + t144 * t121;
t6 = -t56 * t189 + (t198 * t28 + t201 * t29) * t187;
t270 = t6 * t246 + t5 * t247;
t269 = -t189 / 0.2e1;
t268 = t198 / 0.2e1;
t267 = -t201 / 0.2e1;
t266 = pkin(1) * t197;
t265 = pkin(1) * t200;
t264 = pkin(2) * t189;
t263 = -pkin(2) + t185;
t260 = t201 * rSges(4,3);
t35 = -t189 * t86 + (-t186 * t88 + t188 * t90) * t187;
t259 = t35 * t201;
t36 = -t189 * t87 + (-t186 * t89 + t188 * t91) * t187;
t258 = t36 * t198;
t45 = -t189 * t99 + (-t101 * t196 + t103 * t199) * t187;
t257 = t45 * t201;
t46 = -t189 * t100 + (-t102 * t196 + t104 * t199) * t187;
t256 = t46 * t198;
t110 = t187 * t188 * t121;
t249 = t120 * t186;
t62 = -t189 * t119 - t187 * t249 + t110;
t255 = t62 * t189;
t122 = -t189 * rSges(6,3) + (rSges(6,1) * t188 - rSges(6,2) * t186) * t187;
t218 = -t142 * rSges(6,1) - t141 * rSges(6,2);
t92 = rSges(6,3) * t247 - t218;
t65 = t122 * t247 + t189 * t92;
t253 = Icges(3,4) * t197;
t252 = Icges(3,4) * t200;
t126 = -Icges(5,6) * t189 + (Icges(5,4) * t199 - Icges(5,2) * t196) * t187;
t248 = t126 * t196;
t156 = t263 * t187;
t238 = -t122 - t156;
t205 = rSges(4,1) * t245 - rSges(4,2) * t246 + t198 * rSges(4,3);
t220 = rSges(4,1) * t189 - rSges(4,2) * t187;
t94 = t198 * (t198 * t220 - t260) + t201 * t205;
t128 = -t189 * rSges(5,3) + (rSges(5,1) * t199 - rSges(5,2) * t196) * t187;
t167 = pkin(2) * t187 - pkin(5) * t189;
t237 = -t128 - t167;
t178 = pkin(5) * t246;
t234 = t179 + t178;
t236 = t192 * (pkin(5) * t187 + t264) + t201 * t234;
t233 = t232 * t265;
t16 = t29 * t198 - t201 * t28;
t40 = t159 * t101 + t160 * t103 + t246 * t99;
t41 = t100 * t246 + t159 * t102 + t160 * t104;
t22 = t41 * t198 - t201 * t40;
t231 = (t192 * t132 + t16 + t22 + (t211 * t201 + (-t131 + t210) * t198) * t201) * t198;
t230 = pkin(3) * t244;
t229 = -t167 + t238;
t106 = t160 * rSges(5,1) + t159 * rSges(5,2) + rSges(5,3) * t246;
t228 = t247 / 0.2e1;
t227 = t246 / 0.2e1;
t166 = rSges(4,1) * t187 + rSges(4,2) * t189;
t226 = -t166 - t266;
t225 = -t167 - t266;
t224 = (t35 + t55) * t228 + (t36 + t56) * t227;
t219 = -t158 * rSges(5,1) - t157 * rSges(5,2);
t105 = rSges(5,3) * t247 - t219;
t52 = t198 * t105 + t201 * t106 + t236;
t9 = -t255 + (t198 * t35 + t201 * t36) * t187;
t223 = -t189 * t9 + t270;
t222 = t15 * t228 + t16 * t227 + t5 * t267 + t6 * t268 + (t258 - t259) * t269;
t221 = -t128 + t225;
t217 = Icges(3,1) * t200 - t253;
t215 = -Icges(3,2) * t197 + t252;
t213 = Icges(3,5) * t200 - Icges(3,6) * t197;
t163 = Icges(4,2) * t189 + t251;
t164 = Icges(4,1) * t187 + t250;
t207 = -t163 * t187 + t164 * t189;
t206 = t225 + t238;
t115 = t189 * t198 * t263 - t230;
t25 = t236 + t273 * t201 + (t115 + t92) * t198;
t204 = t201 * t271 + t231;
t125 = -Icges(5,3) * t189 + (Icges(5,5) * t199 - Icges(5,6) * t196) * t187;
t127 = -Icges(5,5) * t189 + (Icges(5,1) * t199 - Icges(5,4) * t196) * t187;
t60 = t125 * t247 + t126 * t157 + t127 * t158;
t10 = -t60 * t189 + (t198 * t38 + t201 * t39) * t187;
t61 = t125 * t246 + t159 * t126 + t160 * t127;
t11 = -t61 * t189 + (t198 * t40 + t201 * t41) * t187;
t203 = t10 * t267 + t11 * t268 + t21 * t228 + t22 * t227 + t222 + (t256 - t257) * t269;
t162 = Icges(4,5) * t187 + Icges(4,6) * t189;
t202 = -t259 / 0.2e1 + t258 / 0.2e1 - t257 / 0.2e1 + t256 / 0.2e1 + (t134 * t189 + t136 * t187 + t198 * t162 + t201 * t207 + t56 + t61) * t268 + (t133 * t189 + t135 * t187 - t201 * t162 + t198 * t207 + t55 + t60) * t267;
t184 = t201 * t265;
t174 = rSges(2,1) * t201 - t198 * rSges(2,2);
t173 = -t198 * rSges(2,1) - rSges(2,2) * t201;
t172 = rSges(3,1) * t197 + rSges(3,2) * t200;
t153 = t198 * rSges(3,3) + t201 * t272;
t151 = -t201 * rSges(3,3) + t198 * t272;
t146 = Icges(3,3) * t198 + t201 * t213;
t145 = -Icges(3,3) * t201 + t198 * t213;
t130 = t226 * t201;
t129 = t226 * t198;
t118 = t184 + t205;
t117 = t260 + (-t220 - t265) * t198;
t112 = t187 * t199 * t127;
t109 = t237 * t201;
t108 = t237 * t198;
t107 = t198 * t151 + t153 * t201;
t96 = t221 * t201;
t95 = t221 * t198;
t80 = t92 * t246;
t77 = t233 + t94;
t76 = t229 * t201;
t75 = t229 * t198;
t74 = t184 + t106 + t234;
t73 = (-t265 - t264 + (-rSges(5,3) - pkin(5)) * t187) * t198 + t219;
t72 = t206 * t201;
t71 = t206 * t198;
t70 = -t189 * t106 - t128 * t246;
t69 = t105 * t189 + t128 * t247;
t68 = t184 + t178 + t274;
t67 = t230 + (-t265 - t185 * t189 + (-rSges(6,3) - pkin(5)) * t187) * t198 + t218;
t66 = -t122 * t246 - t189 * t93;
t64 = -t189 * t125 - t187 * t248 + t112;
t63 = (t105 * t201 - t106 * t198) * t187;
t57 = -t247 * t93 + t80;
t49 = -t189 * t273 + t238 * t246;
t48 = t115 * t189 + t156 * t247 + t65;
t47 = t52 + t233;
t30 = t80 + (t115 * t201 - t198 * t273) * t187;
t24 = t25 + t233;
t1 = [t200 * (Icges(3,2) * t200 + t253) + t197 * (Icges(3,1) * t197 + t252) + Icges(2,3) + t110 + t112 + (-t119 - t125 + t163) * t189 + (t164 - t248 - t249) * t187 + m(6) * (t67 ^ 2 + t68 ^ 2) + m(5) * (t73 ^ 2 + t74 ^ 2) + m(4) * (t117 ^ 2 + t118 ^ 2) + m(3) * (t151 ^ 2 + t153 ^ 2) + m(2) * (t173 ^ 2 + t174 ^ 2); (t193 / 0.2e1 + t192 / 0.2e1) * (Icges(3,5) * t197 + Icges(3,6) * t200) + m(3) * (t151 * t201 - t153 * t198) * t172 + t202 + m(6) * (t67 * t72 + t68 * t71) + m(5) * (t73 * t96 + t74 * t95) + m(4) * (t117 * t130 + t118 * t129) + (t200 * (Icges(3,6) * t198 + t201 * t215) + t197 * (Icges(3,5) * t198 + t201 * t217)) * t268 + (t200 * (-Icges(3,6) * t201 + t198 * t215) + t197 * (-Icges(3,5) * t201 + t198 * t217)) * t267; m(6) * (t24 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(5) * (t47 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(4) * (t129 ^ 2 + t130 ^ 2 + t77 ^ 2) + m(3) * (t172 ^ 2 * t232 + t107 ^ 2) + t198 * t192 * t146 + t231 + (-t193 * t145 + (-t198 * t145 + t201 * t146) * t198 + t271) * t201; m(6) * (t67 * t76 + t68 * t75) + m(5) * (t108 * t74 + t109 * t73) + m(4) * (-t117 * t201 - t118 * t198) * t166 + t202; m(6) * (t25 * t24 + t71 * t75 + t72 * t76) + m(5) * (t108 * t95 + t109 * t96 + t52 * t47) + m(4) * (t94 * t77 + (-t129 * t198 - t130 * t201) * t166) + t204; m(6) * (t25 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(5) * (t108 ^ 2 + t109 ^ 2 + t52 ^ 2) + m(4) * (t166 ^ 2 * t232 + t94 ^ 2) + t204; (-t62 - t64) * t189 + m(6) * (t48 * t67 + t49 * t68) + m(5) * (t69 * t73 + t70 * t74) + ((t61 / 0.2e1 + t46 / 0.2e1) * t201 + (t60 / 0.2e1 + t45 / 0.2e1) * t198) * t187 + t224; m(6) * (t30 * t24 + t48 * t72 + t49 * t71) + m(5) * (t63 * t47 + t69 * t96 + t70 * t95) + t203; m(6) * (t30 * t25 + t48 * t76 + t49 * t75) + m(5) * (t108 * t70 + t109 * t69 + t63 * t52) + t203; (t64 * t189 - t9) * t189 + m(6) * (t30 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t63 ^ 2 + t69 ^ 2 + t70 ^ 2) + (t201 * t11 + t198 * t10 - t189 * (t198 * t45 + t201 * t46)) * t187 + t270; m(6) * (t65 * t67 + t66 * t68) - t255 + t224; m(6) * (t57 * t24 + t65 * t72 + t66 * t71) + t222; m(6) * (t57 * t25 + t65 * t76 + t66 * t75) + t222; m(6) * (t57 * t30 + t48 * t65 + t49 * t66) + t223; m(6) * (t57 ^ 2 + t65 ^ 2 + t66 ^ 2) + t223;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq  = res;
