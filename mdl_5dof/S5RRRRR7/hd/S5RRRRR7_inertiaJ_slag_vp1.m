% Calculate joint inertia matrix for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:23
% EndTime: 2019-12-31 22:21:28
% DurationCPUTime: 1.84s
% Computational Cost: add. (7193->328), mult. (6737->485), div. (0->0), fcn. (7074->10), ass. (0->177)
t161 = sin(qJ(1));
t156 = t161 ^ 2;
t239 = t161 * pkin(6);
t164 = cos(qJ(1));
t158 = qJ(2) + qJ(3);
t147 = sin(t158);
t148 = cos(t158);
t181 = Icges(4,5) * t148 - Icges(4,6) * t147;
t100 = Icges(4,3) * t161 + t164 * t181;
t157 = t164 ^ 2;
t215 = Icges(4,4) * t148;
t184 = -Icges(4,2) * t147 + t215;
t102 = Icges(4,6) * t161 + t164 * t184;
t216 = Icges(4,4) * t147;
t187 = Icges(4,1) * t148 - t216;
t104 = Icges(4,5) * t161 + t164 * t187;
t178 = -t102 * t147 + t104 * t148;
t101 = -Icges(4,6) * t164 + t161 * t184;
t103 = -Icges(4,5) * t164 + t161 * t187;
t179 = t101 * t147 - t103 * t148;
t149 = qJ(4) + t158;
t144 = sin(t149);
t145 = cos(t149);
t213 = Icges(5,4) * t145;
t183 = -Icges(5,2) * t144 + t213;
t92 = Icges(5,6) * t161 + t164 * t183;
t214 = Icges(5,4) * t144;
t186 = Icges(5,1) * t145 - t214;
t94 = Icges(5,5) * t161 + t164 * t186;
t189 = -t144 * t92 + t145 * t94;
t91 = -Icges(5,6) * t164 + t161 * t183;
t93 = -Icges(5,5) * t164 + t161 * t186;
t190 = t144 * t91 - t145 * t93;
t162 = cos(qJ(5));
t206 = t164 * t162;
t159 = sin(qJ(5));
t209 = t161 * t159;
t113 = -t145 * t209 - t206;
t207 = t164 * t159;
t208 = t161 * t162;
t114 = t145 * t208 - t207;
t212 = t144 * t161;
t55 = Icges(6,5) * t114 + Icges(6,6) * t113 + Icges(6,3) * t212;
t57 = Icges(6,4) * t114 + Icges(6,2) * t113 + Icges(6,6) * t212;
t59 = Icges(6,1) * t114 + Icges(6,4) * t113 + Icges(6,5) * t212;
t15 = t113 * t57 + t114 * t59 + t212 * t55;
t115 = -t145 * t207 + t208;
t116 = t145 * t206 + t209;
t211 = t144 * t164;
t56 = Icges(6,5) * t116 + Icges(6,6) * t115 + Icges(6,3) * t211;
t58 = Icges(6,4) * t116 + Icges(6,2) * t115 + Icges(6,6) * t211;
t60 = Icges(6,1) * t116 + Icges(6,4) * t115 + Icges(6,5) * t211;
t16 = t113 * t58 + t114 * t60 + t212 * t56;
t8 = -t15 * t164 + t16 * t161;
t180 = Icges(5,5) * t145 - Icges(5,6) * t144;
t89 = -Icges(5,3) * t164 + t161 * t180;
t90 = Icges(5,3) * t161 + t164 * t180;
t233 = -t157 * t89 - (t189 * t161 + (t190 - t90) * t164) * t161 - t8;
t99 = -Icges(4,3) * t164 + t161 * t181;
t238 = -t157 * t99 - (t178 * t161 + (-t100 + t179) * t164) * t161 + t233;
t210 = t145 * t164;
t64 = rSges(6,1) * t116 + rSges(6,2) * t115 + rSges(6,3) * t211;
t237 = pkin(4) * t210 + pkin(9) * t211 + t64;
t193 = rSges(4,1) * t148 - rSges(4,2) * t147;
t165 = -pkin(7) - pkin(6);
t236 = t161 / 0.2e1;
t235 = -t164 / 0.2e1;
t17 = t115 * t57 + t116 * t59 + t211 * t55;
t18 = t115 * t58 + t116 * t60 + t211 * t56;
t9 = t161 * t18 - t164 * t17;
t234 = (t156 * t90 + t9 + (t190 * t164 + (t189 - t89) * t161) * t164) * t161;
t160 = sin(qJ(2));
t232 = pkin(2) * t160;
t231 = pkin(3) * t147;
t230 = pkin(4) * t145;
t163 = cos(qJ(2));
t146 = t163 * pkin(2) + pkin(1);
t129 = pkin(3) * t148 + t146;
t123 = t164 * t129;
t141 = t164 * t146;
t229 = t164 * (t123 - t141) + (t129 - t146) * t156;
t171 = rSges(5,1) * t210 - rSges(5,2) * t211 + t161 * rSges(5,3);
t192 = rSges(5,1) * t145 - rSges(5,2) * t144;
t51 = t161 * (-t164 * rSges(5,3) + t161 * t192) + t164 * t171;
t154 = t164 * pkin(6);
t228 = t161 * (t154 + (-pkin(1) + t146) * t161) + t164 * (-t164 * pkin(1) + t141 - t239);
t172 = t161 * rSges(4,3) + t164 * t193;
t54 = t161 * (-t164 * rSges(4,3) + t161 * t193) + t164 * t172;
t227 = rSges(3,1) * t163;
t225 = rSges(3,2) * t160;
t76 = -Icges(6,6) * t145 + (Icges(6,4) * t162 - Icges(6,2) * t159) * t144;
t223 = t159 * t76;
t222 = t164 * rSges(3,3);
t24 = -t145 * t55 + (-t159 * t57 + t162 * t59) * t144;
t221 = t24 * t164;
t25 = -t145 * t56 + (-t159 * t58 + t162 * t60) * t144;
t220 = t25 * t161;
t80 = -t145 * rSges(6,3) + (rSges(6,1) * t162 - rSges(6,2) * t159) * t144;
t219 = -pkin(4) * t144 + pkin(9) * t145 - t80;
t218 = Icges(3,4) * t160;
t217 = Icges(3,4) * t163;
t204 = t161 * rSges(3,3) + t164 * t227;
t202 = t156 + t157;
t201 = t161 * (t156 * t100 + (t179 * t164 + (t178 - t99) * t161) * t164) + t234;
t128 = rSges(4,1) * t147 + rSges(4,2) * t148;
t200 = -t128 - t232;
t121 = rSges(5,1) * t144 + rSges(5,2) * t145;
t199 = -t121 - t231;
t29 = t51 + t229;
t191 = -t114 * rSges(6,1) - t113 * rSges(6,2);
t63 = rSges(6,3) * t212 - t191;
t26 = t161 * t63 + t156 * (pkin(9) * t144 + t230) + t237 * t164;
t155 = -pkin(8) + t165;
t198 = -t155 * t161 + t123;
t75 = -Icges(6,3) * t145 + (Icges(6,5) * t162 - Icges(6,6) * t159) * t144;
t77 = -Icges(6,5) * t145 + (Icges(6,1) * t162 - Icges(6,4) * t159) * t144;
t30 = t113 * t76 + t114 * t77 + t212 * t75;
t3 = -t30 * t145 + (t15 * t161 + t16 * t164) * t144;
t31 = t115 * t76 + t116 * t77 + t211 * t75;
t4 = -t31 * t145 + (t161 * t17 + t164 * t18) * t144;
t197 = t8 * t212 / 0.2e1 + t3 * t235 + t4 * t236 - t145 * (t220 - t221) / 0.2e1 + t9 * t211 / 0.2e1;
t196 = t219 - t231;
t195 = -t231 - t232;
t194 = -t225 + t227;
t12 = t26 + t229;
t188 = Icges(3,1) * t163 - t218;
t185 = -Icges(3,2) * t160 + t217;
t182 = Icges(3,5) * t163 - Icges(3,6) * t160;
t119 = Icges(5,2) * t145 + t214;
t120 = Icges(5,1) * t144 + t213;
t175 = -t119 * t144 + t120 * t145;
t126 = Icges(4,2) * t148 + t216;
t127 = Icges(4,1) * t147 + t215;
t174 = -t126 * t147 + t127 * t148;
t173 = t164 * t233 + t234;
t170 = -t121 + t195;
t169 = t195 + t219;
t118 = Icges(5,5) * t144 + Icges(5,6) * t145;
t168 = t220 / 0.2e1 - t221 / 0.2e1 + (t161 * t118 + t144 * t94 + t145 * t92 + t164 * t175 + t31) * t236 + (-t164 * t118 + t144 * t93 + t145 * t91 + t161 * t175 + t30) * t235;
t167 = t164 * t238 + t201;
t125 = Icges(4,5) * t147 + Icges(4,6) * t148;
t166 = t168 + (t102 * t148 + t104 * t147 + t161 * t125 + t164 * t174) * t236 + (t101 * t148 + t103 * t147 - t164 * t125 + t161 * t174) * t235;
t140 = rSges(2,1) * t164 - rSges(2,2) * t161;
t139 = -rSges(2,1) * t161 - rSges(2,2) * t164;
t138 = rSges(3,1) * t160 + rSges(3,2) * t163;
t108 = Icges(3,3) * t161 + t164 * t182;
t107 = -Icges(3,3) * t164 + t161 * t182;
t96 = t200 * t164;
t95 = t200 * t161;
t84 = t239 + (pkin(1) - t225) * t164 + t204;
t83 = t222 + t154 + (-pkin(1) - t194) * t161;
t82 = t199 * t164;
t81 = t199 * t161;
t74 = t170 * t164;
t73 = t170 * t161;
t72 = -t161 * t165 + t141 + t172;
t71 = (rSges(4,3) - t165) * t164 + (-t146 - t193) * t161;
t70 = t144 * t162 * t77;
t67 = t164 * (-t164 * t225 + t204) + (t161 * t194 - t222) * t161;
t66 = t171 + t198;
t65 = (rSges(5,3) - t155) * t164 + (-t129 - t192) * t161;
t62 = t219 * t164;
t61 = t219 * t161;
t50 = t196 * t164;
t49 = t196 * t161;
t44 = t169 * t164;
t43 = t169 * t161;
t38 = t198 + t237;
t37 = -t164 * t155 + (-t230 - t129 + (-rSges(6,3) - pkin(9)) * t144) * t161 + t191;
t36 = -t145 * t64 - t211 * t80;
t35 = t145 * t63 + t212 * t80;
t34 = t54 + t228;
t33 = -t144 * t223 - t145 * t75 + t70;
t32 = (-t161 * t64 + t164 * t63) * t144;
t19 = t29 + t228;
t11 = t12 + t228;
t1 = [t148 * t126 + t147 * t127 + t163 * (Icges(3,2) * t163 + t218) + t160 * (Icges(3,1) * t160 + t217) + Icges(2,3) + t70 + (-t75 + t119) * t145 + (t120 - t223) * t144 + m(6) * (t37 ^ 2 + t38 ^ 2) + m(5) * (t65 ^ 2 + t66 ^ 2) + m(4) * (t71 ^ 2 + t72 ^ 2) + m(3) * (t83 ^ 2 + t84 ^ 2) + m(2) * (t139 ^ 2 + t140 ^ 2); t166 + (t156 / 0.2e1 + t157 / 0.2e1) * (Icges(3,5) * t160 + Icges(3,6) * t163) + m(3) * (-t161 * t84 - t164 * t83) * t138 + m(6) * (t37 * t44 + t38 * t43) + m(5) * (t65 * t74 + t66 * t73) + m(4) * (t71 * t96 + t72 * t95) + (t163 * (Icges(3,6) * t161 + t164 * t185) + t160 * (Icges(3,5) * t161 + t164 * t188)) * t236 + (t163 * (-Icges(3,6) * t164 + t161 * t185) + t160 * (-Icges(3,5) * t164 + t161 * t188)) * t235; m(6) * (t11 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t19 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(4) * (t34 ^ 2 + t95 ^ 2 + t96 ^ 2) + t161 * t156 * t108 + m(3) * (t138 ^ 2 * t202 + t67 ^ 2) + t201 + (-t157 * t107 + (-t161 * t107 + t164 * t108) * t161 + t238) * t164; t166 + m(6) * (t37 * t50 + t38 * t49) + m(5) * (t65 * t82 + t66 * t81) + m(4) * (-t161 * t72 - t164 * t71) * t128; m(6) * (t11 * t12 + t43 * t49 + t44 * t50) + m(5) * (t19 * t29 + t73 * t81 + t74 * t82) + m(4) * (t54 * t34 + (-t161 * t95 - t164 * t96) * t128) + t167; m(6) * (t12 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(5) * (t29 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t128 ^ 2 * t202 + t54 ^ 2) + t167; m(6) * (t37 * t62 + t38 * t61) + m(5) * (-t161 * t66 - t164 * t65) * t121 + t168; m(6) * (t11 * t26 + t43 * t61 + t44 * t62) + m(5) * (t51 * t19 + (-t161 * t73 - t164 * t74) * t121) + t173; m(6) * (t12 * t26 + t49 * t61 + t50 * t62) + m(5) * (t51 * t29 + (-t161 * t81 - t164 * t82) * t121) + t173; m(5) * (t121 ^ 2 * t202 + t51 ^ 2) + m(6) * (t26 ^ 2 + t61 ^ 2 + t62 ^ 2) + t173; m(6) * (t35 * t37 + t36 * t38) - t33 * t145 + ((t31 / 0.2e1 + t25 / 0.2e1) * t164 + (t30 / 0.2e1 + t24 / 0.2e1) * t161) * t144; m(6) * (t11 * t32 + t35 * t44 + t36 * t43) + t197; m(6) * (t12 * t32 + t35 * t50 + t36 * t49) + t197; m(6) * (t26 * t32 + t35 * t62 + t36 * t61) + t197; m(6) * (t32 ^ 2 + t35 ^ 2 + t36 ^ 2) + t145 ^ 2 * t33 + (t164 * t4 + t161 * t3 - t145 * (t161 * t24 + t164 * t25)) * t144;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
