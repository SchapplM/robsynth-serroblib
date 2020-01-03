% Calculate joint inertia matrix for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR15_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR15_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR15_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:12
% EndTime: 2019-12-31 20:41:19
% DurationCPUTime: 2.54s
% Computational Cost: add. (4020->355), mult. (6579->514), div. (0->0), fcn. (7051->8), ass. (0->177)
t238 = Icges(4,1) + Icges(3,3);
t163 = sin(qJ(2));
t166 = cos(qJ(2));
t237 = (-Icges(4,4) + Icges(3,5)) * t166 + (Icges(4,5) - Icges(3,6)) * t163;
t164 = sin(qJ(1));
t236 = -t164 / 0.2e1;
t228 = t164 / 0.2e1;
t167 = cos(qJ(1));
t227 = -t167 / 0.2e1;
t235 = t167 / 0.2e1;
t229 = t163 / 0.2e1;
t234 = t238 * t164 + t237 * t167;
t233 = -t237 * t164 + t238 * t167;
t159 = t164 ^ 2;
t160 = t167 ^ 2;
t232 = m(4) / 0.2e1;
t231 = m(5) / 0.2e1;
t230 = m(6) / 0.2e1;
t136 = rSges(3,1) * t163 + rSges(3,2) * t166;
t226 = m(3) * t136;
t162 = sin(qJ(4));
t225 = pkin(4) * t162;
t168 = -pkin(8) - pkin(7);
t224 = -pkin(7) - t168;
t161 = qJ(4) + qJ(5);
t149 = sin(t161);
t150 = cos(t161);
t212 = t163 * t167;
t111 = -t149 * t164 + t150 * t212;
t112 = t149 * t212 + t150 * t164;
t208 = t166 * t167;
t65 = t112 * rSges(6,1) + t111 * rSges(6,2) + rSges(6,3) * t208;
t165 = cos(qJ(4));
t148 = pkin(4) * t165 + pkin(3);
t197 = t162 * t212;
t169 = pkin(4) * t197 + t164 * t148 - t168 * t208;
t202 = t164 * pkin(3) + pkin(7) * t208;
t80 = t169 - t202;
t223 = t65 + t80;
t213 = t163 * t164;
t113 = t149 * t167 + t150 * t213;
t114 = t149 * t213 - t150 * t167;
t185 = -t114 * rSges(6,1) - t113 * rSges(6,2);
t210 = t164 * t166;
t66 = rSges(6,3) * t210 - t185;
t156 = t167 * pkin(3);
t207 = t167 * t148;
t81 = -t207 + t156 + (t163 * t225 + t224 * t166) * t164;
t222 = t66 + t81;
t221 = rSges(4,1) * t167;
t220 = t167 * rSges(3,3);
t117 = t224 * t163 - t166 * t225;
t91 = rSges(6,3) * t163 + (-rSges(6,1) * t149 - rSges(6,2) * t150) * t166;
t219 = -t117 - t91;
t218 = Icges(3,4) * t163;
t217 = Icges(3,4) * t166;
t216 = Icges(4,6) * t163;
t215 = Icges(4,6) * t166;
t214 = qJ(3) * t163;
t211 = t164 * t165;
t209 = t165 * t167;
t204 = pkin(2) * t208 + qJ(3) * t212;
t206 = t159 * (pkin(2) * t166 + t214) + t167 * t204;
t134 = pkin(2) * t163 - qJ(3) * t166;
t205 = rSges(4,2) * t163 + rSges(4,3) * t166 - t134;
t203 = t167 * pkin(1) + t164 * pkin(6);
t201 = t159 + t160;
t59 = Icges(6,5) * t112 + Icges(6,6) * t111 + Icges(6,3) * t208;
t61 = Icges(6,4) * t112 + Icges(6,2) * t111 + Icges(6,6) * t208;
t63 = Icges(6,1) * t112 + Icges(6,4) * t111 + Icges(6,5) * t208;
t25 = t163 * t59 + (-t149 * t63 - t150 * t61) * t166;
t60 = Icges(6,5) * t114 + Icges(6,6) * t113 + Icges(6,3) * t210;
t62 = Icges(6,4) * t114 + Icges(6,2) * t113 + Icges(6,6) * t210;
t64 = Icges(6,1) * t114 + Icges(6,4) * t113 + Icges(6,5) * t210;
t26 = t163 * t60 + (-t149 * t64 - t150 * t62) * t166;
t89 = Icges(6,6) * t163 + (-Icges(6,4) * t149 - Icges(6,2) * t150) * t166;
t90 = Icges(6,5) * t163 + (-Icges(6,1) * t149 - Icges(6,4) * t150) * t166;
t184 = -t149 * t90 - t150 * t89;
t88 = Icges(6,3) * t163 + (-Icges(6,5) * t149 - Icges(6,6) * t150) * t166;
t85 = t163 * t88;
t46 = (t184 * t166 + t85) * t163;
t19 = t111 * t61 + t112 * t63 + t59 * t208;
t20 = t111 * t62 + t112 * t64 + t60 * t208;
t38 = t111 * t89 + t112 * t90 + t88 * t208;
t5 = t163 * t38 + (t164 * t20 + t167 * t19) * t166;
t21 = t113 * t61 + t114 * t63 + t59 * t210;
t22 = t113 * t62 + t114 * t64 + t60 * t210;
t39 = t113 * t89 + t114 * t90 + t88 * t210;
t6 = t163 * t39 + (t164 * t22 + t167 * t21) * t166;
t200 = t5 * t208 + t6 * t210 + t163 * (t46 + (t164 * t26 + t167 * t25) * t166);
t123 = t162 * t167 + t163 * t211;
t124 = t162 * t213 - t209;
t71 = Icges(5,5) * t124 + Icges(5,6) * t123 + Icges(5,3) * t210;
t73 = Icges(5,4) * t124 + Icges(5,2) * t123 + Icges(5,6) * t210;
t75 = Icges(5,1) * t124 + Icges(5,4) * t123 + Icges(5,5) * t210;
t35 = t163 * t71 + (-t162 * t75 - t165 * t73) * t166;
t101 = Icges(5,5) * t163 + (-Icges(5,1) * t162 - Icges(5,4) * t165) * t166;
t95 = Icges(5,3) * t163 + (-Icges(5,5) * t162 - Icges(5,6) * t165) * t166;
t98 = Icges(5,6) * t163 + (-Icges(5,4) * t162 - Icges(5,2) * t165) * t166;
t44 = t101 * t124 + t123 * t98 + t95 * t210;
t199 = t35 / 0.2e1 + t44 / 0.2e1;
t121 = -t162 * t164 + t163 * t209;
t122 = t197 + t211;
t70 = Icges(5,5) * t122 + Icges(5,6) * t121 + Icges(5,3) * t208;
t72 = Icges(5,4) * t122 + Icges(5,2) * t121 + Icges(5,6) * t208;
t74 = Icges(5,1) * t122 + Icges(5,4) * t121 + Icges(5,5) * t208;
t34 = t163 * t70 + (-t162 * t74 - t165 * t72) * t166;
t43 = t101 * t122 + t121 * t98 + t95 * t208;
t198 = t43 / 0.2e1 + t34 / 0.2e1;
t76 = t122 * rSges(5,1) + t121 * rSges(5,2) + rSges(5,3) * t208;
t196 = t210 / 0.2e1;
t195 = t208 / 0.2e1;
t194 = Icges(3,5) * t229 - Icges(4,4) * t163 / 0.2e1 + (Icges(3,6) / 0.2e1 - Icges(4,5) / 0.2e1) * t166;
t193 = -t163 * pkin(7) - t134;
t192 = t164 * (pkin(7) * t210 - t156) + t167 * t202 + t206;
t191 = t203 + t204;
t12 = t164 * t19 - t167 * t20;
t13 = t164 * t21 - t167 * t22;
t190 = t12 * t195 + t13 * t196 + t6 * t227 + t5 * t228 + (t25 * t164 - t26 * t167) * t229;
t189 = t46 + (t26 + t39) * t196 + (t25 + t38) * t195;
t110 = rSges(5,3) * t163 + (-rSges(5,1) * t162 - rSges(5,2) * t165) * t166;
t188 = -t110 + t193;
t187 = rSges(3,1) * t166 - rSges(3,2) * t163;
t186 = -rSges(5,1) * t124 - rSges(5,2) * t123;
t183 = -t162 * t101 - t165 * t98;
t181 = Icges(3,1) * t166 - t218;
t180 = -Icges(3,2) * t163 + t217;
t177 = -Icges(4,2) * t166 + t216;
t176 = Icges(4,3) * t163 - t215;
t172 = t193 + t219;
t171 = rSges(3,1) * t208 - rSges(3,2) * t212 + t164 * rSges(3,3);
t170 = t164 * rSges(4,1) - rSges(4,2) * t208 + rSges(4,3) * t212;
t155 = t167 * pkin(6);
t138 = rSges(2,1) * t167 - rSges(2,2) * t164;
t137 = -rSges(2,1) * t164 - rSges(2,2) * t167;
t92 = t163 * t95;
t87 = t205 * t167;
t86 = t205 * t164;
t84 = t171 + t203;
t83 = t220 + t155 + (-pkin(1) - t187) * t164;
t82 = t91 * t210;
t79 = t170 + t191;
t78 = t221 + t155 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t166 + (-rSges(4,3) - qJ(3)) * t163) * t164;
t77 = rSges(5,3) * t210 - t186;
t69 = t188 * t167;
t68 = t188 * t164;
t67 = t167 * t171 + (t187 * t164 - t220) * t164;
t58 = t163 * t65;
t57 = t66 * t208;
t56 = t172 * t167;
t55 = t172 * t164;
t54 = -t110 * t208 + t163 * t76;
t53 = t110 * t210 - t163 * t77;
t52 = t191 + t76 + t202;
t51 = t155 + t156 + (-t214 - pkin(1) + (-rSges(5,3) - pkin(2) - pkin(7)) * t166) * t164 + t186;
t50 = (t183 * t166 + t92) * t163;
t49 = t167 * t170 + (-t221 + (-rSges(4,2) * t166 + rSges(4,3) * t163) * t164) * t164 + t206;
t48 = -t91 * t208 + t58;
t47 = -t163 * t66 + t82;
t45 = (-t164 * t76 + t167 * t77) * t166;
t42 = t169 + t191 + t65;
t41 = t207 + t155 + (-pkin(1) + (-qJ(3) - t225) * t163 + (-rSges(6,3) - pkin(2) + t168) * t166) * t164 + t185;
t40 = -t65 * t210 + t57;
t33 = t163 * t80 + t219 * t208 + t58;
t32 = t117 * t210 - t222 * t163 + t82;
t31 = t164 * t77 + t167 * t76 + t192;
t30 = t123 * t73 + t124 * t75 + t71 * t210;
t29 = t123 * t72 + t124 * t74 + t70 * t210;
t28 = t121 * t73 + t122 * t75 + t71 * t208;
t27 = t121 * t72 + t122 * t74 + t70 * t208;
t18 = t57 + (-t223 * t164 + t167 * t81) * t166;
t17 = t222 * t164 + t223 * t167 + t192;
t16 = t164 * t29 - t167 * t30;
t15 = t164 * t27 - t167 * t28;
t9 = t163 * t44 + (t164 * t30 + t167 * t29) * t166;
t8 = t163 * t43 + (t164 * t28 + t167 * t27) * t166;
t1 = [Icges(2,3) + t85 + t92 + m(6) * (t41 ^ 2 + t42 ^ 2) + m(5) * (t51 ^ 2 + t52 ^ 2) + m(4) * (t78 ^ 2 + t79 ^ 2) + m(3) * (t83 ^ 2 + t84 ^ 2) + m(2) * (t137 ^ 2 + t138 ^ 2) + (t183 + t184 + t216 + t218 + (Icges(3,2) + Icges(4,3)) * t166) * t166 + (t215 + t217 + (Icges(3,1) + Icges(4,2)) * t163) * t163; m(6) * (t41 * t56 + t42 * t55) + m(5) * (t51 * t69 + t52 * t68) + m(4) * (t78 * t87 + t79 * t86) + (-t26 / 0.2e1 - t39 / 0.2e1 - t83 * t226 + t194 * t167 + (Icges(4,5) * t227 + Icges(3,6) * t235 + t176 * t228 + t180 * t236) * t166 + (Icges(4,4) * t227 + Icges(3,5) * t235 + t177 * t228 + t181 * t236) * t163 - t199) * t167 + (t38 / 0.2e1 + t25 / 0.2e1 - t84 * t226 + t194 * t164 + (Icges(4,5) * t236 + Icges(3,6) * t228 + t176 * t227 + t180 * t235) * t166 + (Icges(4,4) * t236 + Icges(3,5) * t228 + t177 * t227 + t181 * t235) * t163 + t198) * t164; m(6) * (t17 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t31 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t49 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(3) * (t201 * t136 ^ 2 + t67 ^ 2) + (t233 * t160 - t13 - t16) * t167 + (t12 + t15 + t234 * t159 + (t233 * t164 + t234 * t167) * t167) * t164; 0.2e1 * ((t164 * t42 + t167 * t41) * t230 + (t164 * t52 + t167 * t51) * t231 + (t164 * t79 + t167 * t78) * t232) * t163; m(6) * (-t166 * t17 + (t164 * t55 + t167 * t56) * t163) + m(5) * (-t166 * t31 + (t164 * t68 + t167 * t69) * t163) + m(4) * (-t166 * t49 + (t164 * t86 + t167 * t87) * t163); 0.2e1 * (t232 + t231 + t230) * (t201 * t163 ^ 2 + t166 ^ 2); t50 + m(6) * (t32 * t41 + t33 * t42) + m(5) * (t51 * t53 + t52 * t54) + (t199 * t164 + t198 * t167) * t166 + t189; t8 * t228 + t9 * t227 + (t34 * t164 - t35 * t167) * t229 + (t15 * t235 + t16 * t228) * t166 + m(6) * (t17 * t18 + t32 * t56 + t33 * t55) + m(5) * (t31 * t45 + t53 * t69 + t54 * t68) + t190; m(5) * (-t166 * t45 + (t164 * t54 + t167 * t53) * t163) + m(6) * (-t166 * t18 + (t164 * t33 + t167 * t32) * t163); t163 * t50 + m(6) * (t18 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(5) * (t45 ^ 2 + t53 ^ 2 + t54 ^ 2) + (t167 * t8 + t164 * t9 + t163 * (t164 * t35 + t167 * t34)) * t166 + t200; m(6) * (t41 * t47 + t42 * t48) + t189; m(6) * (t17 * t40 + t47 * t56 + t48 * t55) + t190; m(6) * (-t166 * t40 + (t164 * t48 + t167 * t47) * t163); m(6) * (t18 * t40 + t32 * t47 + t33 * t48) + t200; m(6) * (t40 ^ 2 + t47 ^ 2 + t48 ^ 2) + t200;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
