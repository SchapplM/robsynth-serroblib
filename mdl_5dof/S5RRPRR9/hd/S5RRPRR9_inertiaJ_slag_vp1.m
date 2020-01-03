% Calculate joint inertia matrix for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:03
% EndTime: 2019-12-31 20:20:10
% DurationCPUTime: 2.42s
% Computational Cost: add. (6127->352), mult. (6305->514), div. (0->0), fcn. (6772->10), ass. (0->184)
t244 = Icges(3,3) + Icges(4,3);
t157 = qJ(2) + pkin(9);
t150 = sin(t157);
t151 = cos(t157);
t163 = sin(qJ(2));
t166 = cos(qJ(2));
t243 = Icges(3,5) * t166 + Icges(4,5) * t151 - Icges(3,6) * t163 - Icges(4,6) * t150;
t242 = t150 / 0.2e1;
t241 = t151 / 0.2e1;
t240 = t163 / 0.2e1;
t239 = t166 / 0.2e1;
t165 = cos(qJ(4));
t148 = t165 * pkin(4) + pkin(3);
t168 = -pkin(8) - pkin(7);
t162 = sin(qJ(4));
t164 = sin(qJ(1));
t208 = t164 * t162;
t167 = cos(qJ(1));
t211 = t151 * t167;
t212 = t150 * t167;
t160 = qJ(4) + qJ(5);
t152 = sin(t160);
t206 = t167 * t152;
t153 = cos(t160);
t209 = t164 * t153;
t112 = -t151 * t206 + t209;
t205 = t167 * t153;
t210 = t164 * t152;
t113 = t151 * t205 + t210;
t66 = t113 * rSges(6,1) + t112 * rSges(6,2) + rSges(6,3) * t212;
t238 = pkin(4) * t208 + t148 * t211 - t168 * t212 + t66;
t237 = -t243 * t164 + t244 * t167;
t236 = t244 * t164 + t243 * t167;
t158 = t164 ^ 2;
t159 = t167 ^ 2;
t213 = t150 * t164;
t110 = -t151 * t210 - t205;
t111 = t151 * t209 - t206;
t59 = Icges(6,5) * t111 + Icges(6,6) * t110 + Icges(6,3) * t213;
t61 = Icges(6,4) * t111 + Icges(6,2) * t110 + Icges(6,6) * t213;
t63 = Icges(6,1) * t111 + Icges(6,4) * t110 + Icges(6,5) * t213;
t19 = t110 * t61 + t111 * t63 + t59 * t213;
t60 = Icges(6,5) * t113 + Icges(6,6) * t112 + Icges(6,3) * t212;
t62 = Icges(6,4) * t113 + Icges(6,2) * t112 + Icges(6,6) * t212;
t64 = Icges(6,1) * t113 + Icges(6,4) * t112 + Icges(6,5) * t212;
t20 = t110 * t62 + t111 * t64 + t60 * t213;
t86 = -Icges(6,3) * t151 + (Icges(6,5) * t153 - Icges(6,6) * t152) * t150;
t87 = -Icges(6,6) * t151 + (Icges(6,4) * t153 - Icges(6,2) * t152) * t150;
t88 = -Icges(6,5) * t151 + (Icges(6,1) * t153 - Icges(6,4) * t152) * t150;
t38 = t110 * t87 + t111 * t88 + t86 * t213;
t5 = -t38 * t151 + (t164 * t19 + t167 * t20) * t150;
t21 = t112 * t61 + t113 * t63 + t59 * t212;
t22 = t112 * t62 + t113 * t64 + t60 * t212;
t39 = t112 * t87 + t113 * t88 + t86 * t212;
t6 = -t39 * t151 + (t164 * t21 + t167 * t22) * t150;
t235 = t6 * t212 + t5 * t213;
t234 = -t151 / 0.2e1;
t233 = t164 / 0.2e1;
t232 = -t167 / 0.2e1;
t136 = t163 * rSges(3,1) + t166 * rSges(3,2);
t231 = m(3) * t136;
t230 = pkin(2) * t163;
t229 = pkin(3) * t151;
t228 = -pkin(3) + t148;
t227 = pkin(7) + t168;
t201 = pkin(3) * t211 + pkin(7) * t212;
t226 = -t201 + t238;
t182 = -t111 * rSges(6,1) - t110 * rSges(6,2);
t65 = rSges(6,3) * t213 - t182;
t89 = -t151 * rSges(6,3) + (rSges(6,1) * t153 - rSges(6,2) * t152) * t150;
t49 = t151 * t65 + t89 * t213;
t85 = t228 * t150 + t227 * t151;
t225 = -t85 - t89;
t149 = t166 * pkin(2) + pkin(1);
t144 = t167 * t149;
t156 = t167 * pkin(6);
t161 = -qJ(3) - pkin(6);
t204 = t167 * t161;
t224 = t164 * (t204 + t156 + (-pkin(1) + t149) * t164) + t167 * (-t167 * pkin(1) + t144 + (-pkin(6) - t161) * t164);
t223 = rSges(3,1) * t166;
t222 = rSges(3,2) * t163;
t221 = t152 * t87;
t95 = -Icges(5,6) * t151 + (Icges(5,4) * t165 - Icges(5,2) * t162) * t150;
t220 = t162 * t95;
t219 = t167 * rSges(3,3);
t80 = t150 * t153 * t88;
t43 = -t150 * t221 - t151 * t86 + t80;
t218 = t43 * t151;
t217 = Icges(3,4) * t163;
t216 = Icges(3,4) * t166;
t215 = Icges(4,4) * t150;
t214 = Icges(4,4) * t151;
t207 = t164 * t165;
t203 = t167 * t162;
t202 = t167 * t165;
t200 = t164 * rSges(3,3) + t167 * t223;
t199 = t158 + t159;
t124 = -t151 * t203 + t207;
t125 = t151 * t202 + t208;
t72 = Icges(5,5) * t125 + Icges(5,6) * t124 + Icges(5,3) * t212;
t74 = Icges(5,4) * t125 + Icges(5,2) * t124 + Icges(5,6) * t212;
t76 = Icges(5,1) * t125 + Icges(5,4) * t124 + Icges(5,5) * t212;
t35 = -t151 * t72 + (-t162 * t74 + t165 * t76) * t150;
t94 = -Icges(5,3) * t151 + (Icges(5,5) * t165 - Icges(5,6) * t162) * t150;
t96 = -Icges(5,5) * t151 + (Icges(5,1) * t165 - Icges(5,4) * t162) * t150;
t42 = t124 * t95 + t125 * t96 + t94 * t212;
t198 = t35 / 0.2e1 + t42 / 0.2e1;
t122 = -t151 * t208 - t202;
t123 = t151 * t207 - t203;
t71 = Icges(5,5) * t123 + Icges(5,6) * t122 + Icges(5,3) * t213;
t73 = Icges(5,4) * t123 + Icges(5,2) * t122 + Icges(5,6) * t213;
t75 = Icges(5,1) * t123 + Icges(5,4) * t122 + Icges(5,5) * t213;
t34 = -t151 * t71 + (-t162 * t73 + t165 * t75) * t150;
t41 = t122 * t95 + t123 * t96 + t94 * t213;
t197 = t41 / 0.2e1 + t34 / 0.2e1;
t78 = t125 * rSges(5,1) + t124 * rSges(5,2) + rSges(5,3) * t212;
t196 = t213 / 0.2e1;
t195 = t212 / 0.2e1;
t194 = Icges(3,5) * t240 + Icges(4,5) * t242 + Icges(3,6) * t239 + Icges(4,6) * t241;
t193 = -t150 * rSges(4,1) - t151 * rSges(4,2) - t230;
t192 = -t150 * pkin(3) + t151 * pkin(7) - t230;
t28 = -t151 * t59 + (-t152 * t61 + t153 * t63) * t150;
t29 = -t151 * t60 + (-t152 * t62 + t153 * t64) * t150;
t191 = (t28 + t38) * t196 + (t29 + t39) * t195;
t190 = -t164 * t161 + t144;
t189 = t158 * (pkin(7) * t150 + t229) + t167 * t201 + t224;
t7 = -t218 + (t164 * t28 + t167 * t29) * t150;
t188 = -t151 * t7 + t235;
t12 = t20 * t164 - t19 * t167;
t13 = t22 * t164 - t21 * t167;
t187 = t12 * t196 + t13 * t195 + t5 * t232 + t6 * t233 + (t29 * t164 - t28 * t167) * t234;
t97 = -t151 * rSges(5,3) + (rSges(5,1) * t165 - rSges(5,2) * t162) * t150;
t186 = t192 - t97;
t185 = -t222 + t223;
t184 = rSges(4,1) * t151 - rSges(4,2) * t150;
t183 = -t123 * rSges(5,1) - t122 * rSges(5,2);
t181 = t192 + t225;
t180 = Icges(3,1) * t166 - t217;
t179 = Icges(4,1) * t151 - t215;
t178 = -Icges(3,2) * t163 + t216;
t177 = -Icges(4,2) * t150 + t214;
t170 = rSges(4,1) * t211 - rSges(4,2) * t212 + t164 * rSges(4,3);
t138 = t167 * rSges(2,1) - t164 * rSges(2,2);
t137 = -t164 * rSges(2,1) - t167 * rSges(2,2);
t99 = t193 * t167;
t98 = t193 * t164;
t93 = t164 * pkin(6) + (pkin(1) - t222) * t167 + t200;
t92 = t219 + t156 + (-pkin(1) - t185) * t164;
t84 = t150 * t165 * t96;
t82 = t170 + t190;
t81 = (rSges(4,3) - t161) * t167 + (-t149 - t184) * t164;
t79 = t167 * (-t167 * t222 + t200) + (t185 * t164 - t219) * t164;
t77 = rSges(5,3) * t213 - t183;
t69 = -pkin(4) * t203 + (-t227 * t150 + t228 * t151) * t164;
t68 = t186 * t167;
t67 = t186 * t164;
t57 = t65 * t212;
t56 = t190 + t78 + t201;
t55 = -t204 + (-t229 - t149 + (-rSges(5,3) - pkin(7)) * t150) * t164 + t183;
t54 = -t151 * t78 - t97 * t212;
t53 = t151 * t77 + t97 * t213;
t52 = t181 * t167;
t51 = t181 * t164;
t50 = -t151 * t66 - t89 * t212;
t48 = -t150 * t220 - t151 * t94 + t84;
t47 = t190 + t238;
t46 = (pkin(4) * t162 - t161) * t167 + (-t148 * t151 - t149 + (-rSges(6,3) + t168) * t150) * t164 + t182;
t45 = t167 * t170 + (-t167 * rSges(4,3) + t184 * t164) * t164 + t224;
t44 = (-t164 * t78 + t167 * t77) * t150;
t40 = -t66 * t213 + t57;
t33 = t124 * t74 + t125 * t76 + t72 * t212;
t32 = t124 * t73 + t125 * t75 + t71 * t212;
t31 = t122 * t74 + t123 * t76 + t72 * t213;
t30 = t122 * t73 + t123 * t75 + t71 * t213;
t27 = -t226 * t151 + t225 * t212;
t26 = t151 * t69 + t85 * t213 + t49;
t25 = t164 * t77 + t167 * t78 + t189;
t18 = t57 + (-t226 * t164 + t167 * t69) * t150;
t17 = t226 * t167 + (t65 + t69) * t164 + t189;
t16 = t33 * t164 - t32 * t167;
t15 = t31 * t164 - t30 * t167;
t9 = -t42 * t151 + (t164 * t32 + t167 * t33) * t150;
t8 = -t41 * t151 + (t164 * t30 + t167 * t31) * t150;
t1 = [t166 * (Icges(3,2) * t166 + t217) + t163 * (Icges(3,1) * t163 + t216) + Icges(2,3) + t80 + t84 + (Icges(4,2) * t151 + t215 - t86 - t94) * t151 + (Icges(4,1) * t150 + t214 - t220 - t221) * t150 + m(6) * (t46 ^ 2 + t47 ^ 2) + m(5) * (t55 ^ 2 + t56 ^ 2) + m(4) * (t81 ^ 2 + t82 ^ 2) + m(3) * (t92 ^ 2 + t93 ^ 2) + m(2) * (t137 ^ 2 + t138 ^ 2); m(6) * (t52 * t46 + t51 * t47) + m(5) * (t68 * t55 + t67 * t56) + m(4) * (t99 * t81 + t98 * t82) + (-t38 / 0.2e1 - t28 / 0.2e1 - t92 * t231 - t166 * (-Icges(3,6) * t167 + t178 * t164) / 0.2e1 - t163 * (-Icges(3,5) * t167 + t180 * t164) / 0.2e1 + (-Icges(4,6) * t167 + t177 * t164) * t234 - t150 * (-Icges(4,5) * t167 + t179 * t164) / 0.2e1 + t194 * t167 - t197) * t167 + (t39 / 0.2e1 + t29 / 0.2e1 - t93 * t231 + (Icges(3,6) * t164 + t178 * t167) * t239 + (Icges(3,5) * t164 + t180 * t167) * t240 + (Icges(4,6) * t164 + t177 * t167) * t241 + (Icges(4,5) * t164 + t179 * t167) * t242 + t194 * t164 + t198) * t164; m(6) * (t17 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(5) * (t25 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(4) * (t45 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(3) * (t199 * t136 ^ 2 + t79 ^ 2) + (t237 * t159 - t12 - t15) * t167 + (t13 + t16 + t236 * t158 + (t237 * t164 + t236 * t167) * t167) * t164; m(6) * (t164 * t46 - t167 * t47) + m(5) * (t164 * t55 - t167 * t56) + m(4) * (t164 * t81 - t167 * t82); m(6) * (t164 * t52 - t167 * t51) + m(5) * (t164 * t68 - t167 * t67) + m(4) * (t164 * t99 - t167 * t98); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t199; (-t43 - t48) * t151 + m(6) * (t26 * t46 + t27 * t47) + m(5) * (t53 * t55 + t54 * t56) + (t197 * t164 + t198 * t167) * t150 + t191; t9 * t233 + (t35 * t164 - t34 * t167) * t234 + t8 * t232 + (t15 * t233 + t167 * t16 / 0.2e1) * t150 + m(6) * (t18 * t17 + t26 * t52 + t27 * t51) + m(5) * (t44 * t25 + t53 * t68 + t54 * t67) + t187; m(5) * (t53 * t164 - t54 * t167) + m(6) * (t26 * t164 - t27 * t167); (t48 * t151 - t7) * t151 + m(6) * (t18 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(5) * (t44 ^ 2 + t53 ^ 2 + t54 ^ 2) + (t167 * t9 + t164 * t8 - t151 * (t164 * t34 + t167 * t35)) * t150 + t235; m(6) * (t49 * t46 + t50 * t47) - t218 + t191; m(6) * (t40 * t17 + t49 * t52 + t50 * t51) + t187; m(6) * (t49 * t164 - t50 * t167); m(6) * (t40 * t18 + t49 * t26 + t50 * t27) + t188; m(6) * (t40 ^ 2 + t49 ^ 2 + t50 ^ 2) + t188;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
