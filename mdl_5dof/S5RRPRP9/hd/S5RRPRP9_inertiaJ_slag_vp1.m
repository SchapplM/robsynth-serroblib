% Calculate joint inertia matrix for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:43
% EndTime: 2019-12-31 20:05:49
% DurationCPUTime: 2.17s
% Computational Cost: add. (3838->350), mult. (5485->523), div. (0->0), fcn. (5901->8), ass. (0->162)
t213 = rSges(6,3) + qJ(5);
t215 = rSges(6,1) + pkin(4);
t144 = pkin(8) + qJ(4);
t138 = sin(t144);
t139 = cos(t144);
t154 = cos(qJ(1));
t152 = sin(qJ(1));
t153 = cos(qJ(2));
t180 = t152 * t153;
t98 = t138 * t180 + t139 * t154;
t99 = -t138 * t154 + t139 * t180;
t217 = -t213 * t98 - t215 * t99;
t203 = t152 / 0.2e1;
t216 = t154 / 0.2e1;
t151 = sin(qJ(2));
t82 = -Icges(5,3) * t153 + (Icges(5,5) * t139 - Icges(5,6) * t138) * t151;
t83 = -Icges(6,2) * t153 + (Icges(6,4) * t139 + Icges(6,6) * t138) * t151;
t214 = t82 + t83;
t183 = t151 * t152;
t46 = Icges(6,5) * t99 + Icges(6,6) * t183 + Icges(6,3) * t98;
t50 = Icges(6,4) * t99 + Icges(6,2) * t183 + Icges(6,6) * t98;
t54 = Icges(6,1) * t99 + Icges(6,4) * t183 + Icges(6,5) * t98;
t12 = t50 * t183 + t46 * t98 + t54 * t99;
t179 = t153 * t154;
t100 = t138 * t179 - t152 * t139;
t101 = t152 * t138 + t139 * t179;
t182 = t151 * t154;
t47 = Icges(6,5) * t101 + Icges(6,6) * t182 + Icges(6,3) * t100;
t51 = Icges(6,4) * t101 + Icges(6,2) * t182 + Icges(6,6) * t100;
t55 = Icges(6,1) * t101 + Icges(6,4) * t182 + Icges(6,5) * t100;
t13 = t51 * t183 + t47 * t98 + t55 * t99;
t48 = Icges(5,5) * t99 - Icges(5,6) * t98 + Icges(5,3) * t183;
t52 = Icges(5,4) * t99 - Icges(5,2) * t98 + Icges(5,6) * t183;
t56 = Icges(5,1) * t99 - Icges(5,4) * t98 + Icges(5,5) * t183;
t14 = t48 * t183 - t52 * t98 + t56 * t99;
t49 = Icges(5,5) * t101 - Icges(5,6) * t100 + Icges(5,3) * t182;
t53 = Icges(5,4) * t101 - Icges(5,2) * t100 + Icges(5,6) * t182;
t57 = Icges(5,1) * t101 - Icges(5,4) * t100 + Icges(5,5) * t182;
t15 = t49 * t183 - t53 * t98 + t57 * t99;
t81 = -Icges(6,6) * t153 + (Icges(6,5) * t139 + Icges(6,3) * t138) * t151;
t85 = -Icges(6,4) * t153 + (Icges(6,1) * t139 + Icges(6,5) * t138) * t151;
t29 = t83 * t183 + t81 * t98 + t85 * t99;
t84 = -Icges(5,6) * t153 + (Icges(5,4) * t139 - Icges(5,2) * t138) * t151;
t86 = -Icges(5,5) * t153 + (Icges(5,1) * t139 - Icges(5,4) * t138) * t151;
t30 = t82 * t183 - t84 * t98 + t86 * t99;
t212 = (-t29 - t30) * t153 + ((t13 + t15) * t154 + (t12 + t14) * t152) * t151;
t16 = t100 * t46 + t101 * t54 + t50 * t182;
t17 = t100 * t47 + t101 * t55 + t51 * t182;
t18 = -t100 * t52 + t101 * t56 + t48 * t182;
t19 = -t100 * t53 + t101 * t57 + t49 * t182;
t31 = t100 * t81 + t101 * t85 + t83 * t182;
t32 = -t100 * t84 + t101 * t86 + t82 * t182;
t211 = (-t31 - t32) * t153 + ((t17 + t19) * t154 + (t16 + t18) * t152) * t151;
t20 = -t153 * t50 + (t138 * t46 + t139 * t54) * t151;
t22 = -t153 * t48 + (-t138 * t52 + t139 * t56) * t151;
t210 = -t20 - t22;
t21 = -t153 * t51 + (t138 * t47 + t139 * t55) * t151;
t23 = -t153 * t49 + (-t138 * t53 + t139 * t57) * t151;
t209 = t21 + t23;
t186 = t138 * t151;
t208 = t81 * t186 + (t85 + t86) * t139 * t151;
t146 = t152 ^ 2;
t207 = t153 ^ 2;
t147 = t154 ^ 2;
t206 = m(4) / 0.2e1;
t205 = m(5) / 0.2e1;
t204 = m(6) / 0.2e1;
t201 = -t154 / 0.2e1;
t124 = rSges(3,1) * t151 + rSges(3,2) * t153;
t200 = m(3) * t124;
t199 = pkin(2) * t153;
t149 = cos(pkin(8));
t137 = pkin(3) * t149 + pkin(2);
t198 = -pkin(2) + t137;
t197 = t214 * t153 + t84 * t186 - t208;
t196 = rSges(6,2) * t183 - t217;
t195 = rSges(6,2) * t182 + t213 * t100 + t215 * t101;
t193 = t154 * rSges(3,3);
t192 = -rSges(6,2) * t153 + (t213 * t138 + t215 * t139) * t151;
t123 = pkin(2) * t151 - qJ(3) * t153;
t150 = -pkin(7) - qJ(3);
t191 = -t123 - (qJ(3) + t150) * t153 - t198 * t151;
t148 = sin(pkin(8));
t190 = -t123 + rSges(4,3) * t153 - (rSges(4,1) * t149 - rSges(4,2) * t148) * t151;
t189 = Icges(3,4) * t151;
t188 = Icges(3,4) * t153;
t187 = qJ(3) * t151;
t184 = t148 * t154;
t181 = t152 * t148;
t176 = pkin(2) * t179 + qJ(3) * t182;
t178 = t146 * (t187 + t199) + t154 * t176;
t177 = -pkin(3) * t184 - t150 * t183;
t175 = t154 * pkin(1) + t152 * pkin(6);
t174 = t146 + t147;
t88 = -rSges(5,3) * t153 + (rSges(5,1) * t139 - rSges(5,2) * t138) * t151;
t173 = -t88 + t191;
t61 = t101 * rSges(5,1) - t100 * rSges(5,2) + rSges(5,3) * t182;
t116 = -t148 * t179 + t152 * t149;
t117 = t149 * t179 + t181;
t172 = t117 * rSges(4,1) + t116 * rSges(4,2) + rSges(4,3) * t182;
t142 = t154 * pkin(6);
t171 = t142 - t177;
t170 = -t137 * t153 - pkin(1);
t169 = t191 - t192;
t158 = pkin(3) * t181 + t137 * t179 - t150 * t182;
t168 = t152 * ((t198 * t153 - t187) * t152 + t177) + t154 * (t158 - t176) + t178;
t167 = -rSges(5,1) * t99 + rSges(5,2) * t98;
t166 = rSges(3,1) * t153 - rSges(3,2) * t151;
t114 = -t148 * t180 - t149 * t154;
t115 = t149 * t180 - t184;
t165 = -rSges(4,1) * t115 - rSges(4,2) * t114;
t164 = Icges(3,1) * t153 - t189;
t163 = -Icges(3,2) * t151 + t188;
t162 = Icges(3,5) * t153 - Icges(3,6) * t151;
t159 = rSges(3,1) * t179 - rSges(3,2) * t182 + t152 * rSges(3,3);
t157 = t22 / 0.2e1 + t20 / 0.2e1 + t30 / 0.2e1 + t29 / 0.2e1;
t156 = t23 / 0.2e1 + t21 / 0.2e1 + t32 / 0.2e1 + t31 / 0.2e1;
t155 = t158 + t175;
t145 = t151 ^ 2;
t126 = rSges(2,1) * t154 - t152 * rSges(2,2);
t125 = -t152 * rSges(2,1) - rSges(2,2) * t154;
t120 = Icges(3,5) * t151 + Icges(3,6) * t153;
t103 = Icges(3,3) * t152 + t162 * t154;
t102 = -Icges(3,3) * t154 + t162 * t152;
t96 = -Icges(4,5) * t153 + (Icges(4,1) * t149 - Icges(4,4) * t148) * t151;
t95 = -Icges(4,6) * t153 + (Icges(4,4) * t149 - Icges(4,2) * t148) * t151;
t79 = t159 + t175;
t78 = t193 + t142 + (-pkin(1) - t166) * t152;
t74 = t190 * t154;
t73 = t190 * t152;
t72 = Icges(4,1) * t117 + Icges(4,4) * t116 + Icges(4,5) * t182;
t71 = Icges(4,1) * t115 + Icges(4,4) * t114 + Icges(4,5) * t183;
t70 = Icges(4,4) * t117 + Icges(4,2) * t116 + Icges(4,6) * t182;
t69 = Icges(4,4) * t115 + Icges(4,2) * t114 + Icges(4,6) * t183;
t68 = Icges(4,5) * t117 + Icges(4,6) * t116 + Icges(4,3) * t182;
t67 = Icges(4,5) * t115 + Icges(4,6) * t114 + Icges(4,3) * t183;
t64 = t154 * t159 + (t166 * t152 - t193) * t152;
t59 = rSges(5,3) * t183 - t167;
t45 = t172 + t175 + t176;
t44 = t142 + (-t199 - pkin(1) + (-rSges(4,3) - qJ(3)) * t151) * t152 + t165;
t43 = t173 * t154;
t42 = t173 * t152;
t41 = -t153 * t61 - t88 * t182;
t40 = t153 * t59 + t88 * t183;
t39 = t155 + t61;
t38 = (-rSges(5,3) * t151 + t170) * t152 + t167 + t171;
t37 = t169 * t154;
t36 = t169 * t152;
t33 = (-t152 * t61 + t154 * t59) * t151;
t28 = t152 * (rSges(4,3) * t183 - t165) + t154 * t172 + t178;
t27 = t155 + t195;
t26 = (-rSges(6,2) * t151 + t170) * t152 + t171 + t217;
t25 = -t195 * t153 - t192 * t182;
t24 = t196 * t153 + t192 * t183;
t11 = (-t195 * t152 + t196 * t154) * t151;
t10 = t152 * t59 + t154 * t61 + t168;
t9 = t196 * t152 + t195 * t154 + t168;
t8 = t19 * t152 - t154 * t18;
t7 = t17 * t152 - t154 * t16;
t6 = -t14 * t154 + t15 * t152;
t5 = -t12 * t154 + t13 * t152;
t1 = [Icges(2,3) + (t189 - (Icges(4,5) * t149 - Icges(4,6) * t148) * t151 + (Icges(3,2) + Icges(4,3)) * t153 - t214) * t153 + (Icges(3,1) * t151 - t138 * t84 - t148 * t95 + t149 * t96 + t188) * t151 + m(6) * (t26 ^ 2 + t27 ^ 2) + m(5) * (t38 ^ 2 + t39 ^ 2) + m(3) * (t78 ^ 2 + t79 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + m(2) * (t125 ^ 2 + t126 ^ 2) + t208; m(6) * (t26 * t37 + t27 * t36) + m(5) * (t38 * t43 + t39 * t42) + m(4) * (t44 * t74 + t45 * t73) + (-t78 * t200 - t114 * t95 / 0.2e1 - t115 * t96 / 0.2e1 + t120 * t216 + (Icges(3,6) * t216 - t163 * t152 / 0.2e1 + t67 / 0.2e1) * t153 - t157) * t154 + (-t79 * t200 + t116 * t95 / 0.2e1 + t117 * t96 / 0.2e1 + t120 * t203 + (Icges(3,6) * t203 + t163 * t216 - t68 / 0.2e1) * t153 + t156) * t152 + ((Icges(3,5) * t152 - t148 * t70 + t149 * t72 + t164 * t154) * t203 + (-Icges(3,5) * t154 - t148 * t69 + t149 * t71 + t164 * t152) * t201) * t151; m(6) * (t36 ^ 2 + t37 ^ 2 + t9 ^ 2) + m(5) * (t10 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(4) * (t28 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(3) * (t174 * t124 ^ 2 + t64 ^ 2) + (-t147 * t102 - t5 - t6 + (t114 * t69 + t115 * t71 + t67 * t183) * t154) * t154 + (t8 + t7 + t146 * t103 + (t116 * t70 + t117 * t72 + t68 * t182) * t152 + (-t152 * t102 + t154 * t103 - t114 * t70 - t115 * t72 - t116 * t69 - t117 * t71 - t67 * t182 - t68 * t183) * t154) * t152; 0.2e1 * ((t152 * t27 + t154 * t26) * t204 + (t152 * t39 + t154 * t38) * t205 + (t152 * t45 + t154 * t44) * t206) * t151; m(6) * (-t153 * t9 + (t152 * t36 + t154 * t37) * t151) + m(5) * (-t153 * t10 + (t152 * t42 + t154 * t43) * t151) + m(4) * (-t153 * t28 + (t152 * t73 + t154 * t74) * t151); 0.2e1 * (t206 + t205 + t204) * (t174 * t145 + t207); t197 * t153 + m(6) * (t24 * t26 + t25 * t27) + m(5) * (t38 * t40 + t39 * t41) + (t157 * t152 + t156 * t154) * t151; m(6) * (t11 * t9 + t24 * t37 + t25 * t36) + m(5) * (t10 * t33 + t40 * t43 + t41 * t42) + ((t7 / 0.2e1 + t8 / 0.2e1) * t154 + (t5 / 0.2e1 + t6 / 0.2e1) * t152) * t151 + t211 * t203 - (t209 * t152 + t210 * t154) * t153 / 0.2e1 + t212 * t201; m(5) * (-t33 * t153 + (t152 * t41 + t154 * t40) * t151) + m(6) * (-t11 * t153 + (t152 * t25 + t154 * t24) * t151); m(6) * (t11 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t33 ^ 2 + t40 ^ 2 + t41 ^ 2) - t197 * t207 + ((-t209 * t153 + t211) * t154 + (t210 * t153 + t212) * t152) * t151; m(6) * (t100 * t26 + t27 * t98); m(6) * (t100 * t37 + t9 * t186 + t36 * t98); m(6) * (t100 * t154 - t138 * t153 + t152 * t98) * t151; m(6) * (t100 * t24 + t11 * t186 + t25 * t98); m(6) * (t138 ^ 2 * t145 + t100 ^ 2 + t98 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
