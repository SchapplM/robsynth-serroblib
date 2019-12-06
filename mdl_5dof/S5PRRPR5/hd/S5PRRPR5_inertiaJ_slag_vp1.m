% Calculate joint inertia matrix for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:21
% EndTime: 2019-12-05 16:25:32
% DurationCPUTime: 2.81s
% Computational Cost: add. (11885->415), mult. (21121->623), div. (0->0), fcn. (26534->12), ass. (0->197)
t168 = sin(pkin(5));
t213 = t168 ^ 2;
t167 = sin(pkin(9));
t169 = cos(pkin(9));
t177 = cos(qJ(2));
t170 = cos(pkin(5));
t174 = sin(qJ(2));
t192 = t174 * t170;
t156 = t167 * t177 + t169 * t192;
t188 = qJ(3) + pkin(10);
t166 = sin(t188);
t183 = cos(t188);
t179 = t168 * t183;
t139 = t156 * t166 + t169 * t179;
t158 = -t167 * t192 + t169 * t177;
t141 = t158 * t166 - t167 * t179;
t197 = t168 * t174;
t152 = t166 * t197 - t170 * t183;
t198 = t167 * t168;
t142 = t158 * t183 + t166 * t198;
t191 = t177 * t170;
t157 = t167 * t191 + t169 * t174;
t172 = sin(qJ(5));
t175 = cos(qJ(5));
t113 = -t142 * t172 + t157 * t175;
t114 = t142 * t175 + t157 * t172;
t194 = t169 * t168;
t140 = t156 * t183 - t166 * t194;
t155 = t167 * t174 - t169 * t191;
t111 = -t140 * t172 + t155 * t175;
t112 = t140 * t175 + t155 * t172;
t68 = Icges(6,5) * t112 + Icges(6,6) * t111 + Icges(6,3) * t139;
t70 = Icges(6,4) * t112 + Icges(6,2) * t111 + Icges(6,6) * t139;
t72 = Icges(6,1) * t112 + Icges(6,4) * t111 + Icges(6,5) * t139;
t26 = t113 * t70 + t114 * t72 + t141 * t68;
t69 = Icges(6,5) * t114 + Icges(6,6) * t113 + Icges(6,3) * t141;
t71 = Icges(6,4) * t114 + Icges(6,2) * t113 + Icges(6,6) * t141;
t73 = Icges(6,1) * t114 + Icges(6,4) * t113 + Icges(6,5) * t141;
t27 = t113 * t71 + t114 * t73 + t141 * t69;
t153 = t170 * t166 + t174 * t179;
t195 = t168 * t177;
t143 = -t153 * t172 - t175 * t195;
t144 = t153 * t175 - t172 * t195;
t83 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t152;
t84 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t152;
t85 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t152;
t36 = t113 * t84 + t114 * t85 + t141 * t83;
t2 = t139 * t26 + t141 * t27 + t152 * t36;
t212 = t2 / 0.2e1;
t211 = -m(5) - m(6);
t210 = t139 / 0.2e1;
t209 = t141 / 0.2e1;
t208 = t152 / 0.2e1;
t176 = cos(qJ(3));
t207 = t176 * pkin(3);
t95 = rSges(5,1) * t142 - rSges(5,2) * t141 + rSges(5,3) * t157;
t173 = sin(qJ(3));
t185 = t173 * t198;
t97 = pkin(3) * t185 + qJ(4) * t157 + t158 * t207;
t205 = -t95 - t97;
t74 = rSges(6,1) * t112 + rSges(6,2) * t111 + rSges(6,3) * t139;
t204 = pkin(4) * t140 + pkin(8) * t139 + t74;
t75 = rSges(6,1) * t114 + rSges(6,2) * t113 + rSges(6,3) * t141;
t203 = pkin(4) * t142 + pkin(8) * t141 + t75;
t193 = t170 * t173;
t131 = pkin(3) * t193 + (-qJ(4) * t177 + t174 * t207) * t168;
t184 = t173 * t194;
t96 = -pkin(3) * t184 + qJ(4) * t155 + t156 * t207;
t202 = t155 * t131 + t96 * t195;
t86 = rSges(6,1) * t144 + rSges(6,2) * t143 + rSges(6,3) * t152;
t201 = pkin(4) * t153 + pkin(8) * t152 + t86;
t138 = t158 * pkin(2) + t157 * pkin(7);
t136 = t170 * t138;
t200 = t170 * t97 + t136;
t137 = t156 * pkin(2) + t155 * pkin(7);
t199 = -t137 - t96;
t196 = t168 * t176;
t118 = rSges(5,1) * t153 - rSges(5,2) * t152 - rSges(5,3) * t195;
t190 = -t118 - t131;
t189 = t137 * t198 + t138 * t194;
t187 = -t97 - t203;
t186 = -t131 - t201;
t159 = t170 * t176 - t173 * t197;
t160 = t174 * t196 + t193;
t132 = rSges(4,1) * t160 + rSges(4,2) * t159 - rSges(4,3) * t195;
t161 = (pkin(2) * t174 - pkin(7) * t177) * t168;
t182 = (-t132 - t161) * t168;
t181 = t97 * t194 + t96 * t198 + t189;
t180 = (-t161 + t190) * t168;
t178 = (-t161 + t186) * t168;
t154 = t170 * rSges(3,3) + (rSges(3,1) * t174 + rSges(3,2) * t177) * t168;
t151 = Icges(3,5) * t170 + (Icges(3,1) * t174 + Icges(3,4) * t177) * t168;
t150 = Icges(3,6) * t170 + (Icges(3,4) * t174 + Icges(3,2) * t177) * t168;
t149 = Icges(3,3) * t170 + (Icges(3,5) * t174 + Icges(3,6) * t177) * t168;
t148 = t158 * t176 + t185;
t147 = -t158 * t173 + t167 * t196;
t146 = t156 * t176 - t184;
t145 = -t156 * t173 - t176 * t194;
t130 = Icges(4,1) * t160 + Icges(4,4) * t159 - Icges(4,5) * t195;
t129 = Icges(4,4) * t160 + Icges(4,2) * t159 - Icges(4,6) * t195;
t128 = Icges(4,5) * t160 + Icges(4,6) * t159 - Icges(4,3) * t195;
t127 = rSges(3,1) * t158 - rSges(3,2) * t157 + rSges(3,3) * t198;
t126 = rSges(3,1) * t156 - rSges(3,2) * t155 - rSges(3,3) * t194;
t124 = Icges(3,1) * t158 - Icges(3,4) * t157 + Icges(3,5) * t198;
t123 = Icges(3,1) * t156 - Icges(3,4) * t155 - Icges(3,5) * t194;
t122 = Icges(3,4) * t158 - Icges(3,2) * t157 + Icges(3,6) * t198;
t121 = Icges(3,4) * t156 - Icges(3,2) * t155 - Icges(3,6) * t194;
t120 = Icges(3,5) * t158 - Icges(3,6) * t157 + Icges(3,3) * t198;
t119 = Icges(3,5) * t156 - Icges(3,6) * t155 - Icges(3,3) * t194;
t117 = Icges(5,1) * t153 - Icges(5,4) * t152 - Icges(5,5) * t195;
t116 = Icges(5,4) * t153 - Icges(5,2) * t152 - Icges(5,6) * t195;
t115 = Icges(5,5) * t153 - Icges(5,6) * t152 - Icges(5,3) * t195;
t107 = -t126 * t170 - t154 * t194;
t106 = t127 * t170 - t154 * t198;
t105 = rSges(4,1) * t148 + rSges(4,2) * t147 + rSges(4,3) * t157;
t104 = rSges(4,1) * t146 + rSges(4,2) * t145 + rSges(4,3) * t155;
t103 = Icges(4,1) * t148 + Icges(4,4) * t147 + Icges(4,5) * t157;
t102 = Icges(4,1) * t146 + Icges(4,4) * t145 + Icges(4,5) * t155;
t101 = Icges(4,4) * t148 + Icges(4,2) * t147 + Icges(4,6) * t157;
t100 = Icges(4,4) * t146 + Icges(4,2) * t145 + Icges(4,6) * t155;
t99 = Icges(4,5) * t148 + Icges(4,6) * t147 + Icges(4,3) * t157;
t98 = Icges(4,5) * t146 + Icges(4,6) * t145 + Icges(4,3) * t155;
t94 = rSges(5,1) * t140 - rSges(5,2) * t139 + rSges(5,3) * t155;
t93 = Icges(5,1) * t142 - Icges(5,4) * t141 + Icges(5,5) * t157;
t92 = Icges(5,1) * t140 - Icges(5,4) * t139 + Icges(5,5) * t155;
t91 = Icges(5,4) * t142 - Icges(5,2) * t141 + Icges(5,6) * t157;
t90 = Icges(5,4) * t140 - Icges(5,2) * t139 + Icges(5,6) * t155;
t89 = Icges(5,5) * t142 - Icges(5,6) * t141 + Icges(5,3) * t157;
t88 = Icges(5,5) * t140 - Icges(5,6) * t139 + Icges(5,3) * t155;
t79 = (t126 * t167 + t127 * t169) * t168;
t78 = t157 * t96;
t77 = -t105 * t195 - t132 * t157;
t76 = t104 * t195 + t132 * t155;
t67 = -t128 * t195 + t129 * t159 + t130 * t160;
t66 = t104 * t157 - t105 * t155;
t65 = (-t104 - t137) * t170 + t169 * t182;
t64 = t170 * t105 + t167 * t182 + t136;
t63 = -t115 * t195 - t116 * t152 + t117 * t153;
t62 = t128 * t157 + t129 * t147 + t130 * t148;
t61 = t128 * t155 + t129 * t145 + t130 * t146;
t60 = t115 * t157 - t116 * t141 + t117 * t142;
t59 = t115 * t155 - t116 * t139 + t117 * t140;
t58 = (t104 * t167 + t105 * t169) * t168 + t189;
t57 = t101 * t159 + t103 * t160 - t195 * t99;
t56 = t100 * t159 + t102 * t160 - t195 * t98;
t55 = -t141 * t86 + t152 * t75;
t54 = t139 * t86 - t152 * t74;
t53 = -t152 * t91 + t153 * t93 - t195 * t89;
t52 = -t152 * t90 + t153 * t92 - t195 * t88;
t51 = t157 * t190 + t195 * t205;
t50 = t118 * t155 + t195 * t94 + t202;
t49 = t101 * t147 + t103 * t148 + t157 * t99;
t48 = t100 * t147 + t102 * t148 + t157 * t98;
t47 = t101 * t145 + t103 * t146 + t155 * t99;
t46 = t100 * t145 + t102 * t146 + t155 * t98;
t45 = (-t94 + t199) * t170 + t169 * t180;
t44 = t167 * t180 + t170 * t95 + t200;
t43 = -t141 * t91 + t142 * t93 + t157 * t89;
t42 = -t141 * t90 + t142 * t92 + t157 * t88;
t41 = -t139 * t91 + t140 * t93 + t155 * t89;
t40 = -t139 * t90 + t140 * t92 + t155 * t88;
t39 = t143 * t84 + t144 * t85 + t152 * t83;
t38 = -t139 * t75 + t141 * t74;
t37 = t155 * t205 + t157 * t94 + t78;
t35 = t111 * t84 + t112 * t85 + t139 * t83;
t34 = (t167 * t94 + t169 * t95) * t168 + t181;
t33 = t143 * t71 + t144 * t73 + t152 * t69;
t32 = t143 * t70 + t144 * t72 + t152 * t68;
t31 = t157 * t186 + t187 * t195;
t30 = t155 * t201 + t195 * t204 + t202;
t29 = (t199 - t204) * t170 + t169 * t178;
t28 = t167 * t178 + t170 * t203 + t200;
t25 = t111 * t71 + t112 * t73 + t139 * t69;
t24 = t111 * t70 + t112 * t72 + t139 * t68;
t23 = t155 * t187 + t157 * t204 + t78;
t22 = (t167 * t204 + t169 * t203) * t168 + t181;
t21 = t67 * t170 + (t167 * t57 - t169 * t56) * t168;
t20 = t155 * t56 + t157 * t57 - t195 * t67;
t19 = t63 * t170 + (t167 * t53 - t169 * t52) * t168;
t18 = t62 * t170 + (t167 * t49 - t169 * t48) * t168;
t17 = t61 * t170 + (t167 * t47 - t169 * t46) * t168;
t16 = t155 * t52 + t157 * t53 - t195 * t63;
t15 = t155 * t48 + t157 * t49 - t195 * t62;
t14 = t155 * t46 + t157 * t47 - t195 * t61;
t13 = t60 * t170 + (t167 * t43 - t169 * t42) * t168;
t12 = t59 * t170 + (t167 * t41 - t169 * t40) * t168;
t11 = t155 * t42 + t157 * t43 - t195 * t60;
t10 = t155 * t40 + t157 * t41 - t195 * t59;
t9 = t39 * t170 + (t167 * t33 - t169 * t32) * t168;
t8 = t155 * t32 + t157 * t33 - t195 * t39;
t7 = t139 * t32 + t141 * t33 + t152 * t39;
t6 = t36 * t170 + (t167 * t27 - t169 * t26) * t168;
t5 = t35 * t170 + (t167 * t25 - t169 * t24) * t168;
t4 = t155 * t26 + t157 * t27 - t195 * t36;
t3 = t155 * t24 + t157 * t25 - t195 * t35;
t1 = t139 * t24 + t141 * t25 + t152 * t35;
t80 = [m(2) + m(3) + m(4) - t211; m(3) * t79 + m(4) * t58 + m(5) * t34 + m(6) * t22; m(6) * (t22 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(5) * (t34 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(4) * (t58 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(3) * (t106 ^ 2 + t107 ^ 2 + t79 ^ 2) + (t6 + t18 + t13 + (t120 * t198 - t122 * t157 + t124 * t158) * t198) * t198 + (-t5 - t17 - t12 + (-t119 * t194 - t121 * t155 + t123 * t156) * t194 + (-t119 * t198 + t120 * t194 + t121 * t157 + t122 * t155 - t123 * t158 - t124 * t156) * t198) * t194 + (-(-t149 * t194 - t155 * t150 + t156 * t151) * t194 + (t149 * t198 - t157 * t150 + t158 * t151) * t198 + t9 + t21 + t19 + ((t122 * t177 + t124 * t174) * t167 - (t121 * t177 + t123 * t174) * t169) * t213 + ((-t119 * t169 + t120 * t167 + t150 * t177 + t151 * t174) * t168 + t170 * t149) * t170) * t170; m(4) * t66 + m(5) * t37 + m(6) * t23; (t8 / 0.2e1 + t20 / 0.2e1 + t16 / 0.2e1) * t170 + (t6 / 0.2e1 + t18 / 0.2e1 + t13 / 0.2e1) * t157 + (t5 / 0.2e1 + t17 / 0.2e1 + t12 / 0.2e1) * t155 + m(6) * (t22 * t23 + t28 * t31 + t29 * t30) + m(5) * (t34 * t37 + t44 * t51 + t45 * t50) + m(4) * (t58 * t66 + t64 * t77 + t65 * t76) + ((-t9 / 0.2e1 - t21 / 0.2e1 - t19 / 0.2e1) * t177 + (-t3 / 0.2e1 - t14 / 0.2e1 - t10 / 0.2e1) * t169 + (t4 / 0.2e1 + t15 / 0.2e1 + t11 / 0.2e1) * t167) * t168; (-t16 - t20 - t8) * t195 + (t4 + t15 + t11) * t157 + (t3 + t14 + t10) * t155 + m(6) * (t23 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(5) * (t37 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(4) * (t66 ^ 2 + t76 ^ 2 + t77 ^ 2); t211 * t195; m(6) * (t155 * t28 + t157 * t29 - t195 * t22) + m(5) * (t155 * t44 + t157 * t45 - t195 * t34); m(6) * (t155 * t31 + t157 * t30 - t195 * t23) + m(5) * (t155 * t51 + t157 * t50 - t195 * t37); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t177 ^ 2 * t213 + t155 ^ 2 + t157 ^ 2); m(6) * t38; m(6) * (t22 * t38 + t28 * t55 + t29 * t54) + t170 * t7 / 0.2e1 + t5 * t210 + t6 * t209 + t9 * t208 + (t167 * t212 - t169 * t1 / 0.2e1) * t168; m(6) * (t23 * t38 + t30 * t54 + t31 * t55) + t3 * t210 + t157 * t212 + t4 * t209 + t155 * t1 / 0.2e1 + t8 * t208 - t7 * t195 / 0.2e1; m(6) * (t155 * t55 + t157 * t54 - t195 * t38); m(6) * (t38 ^ 2 + t54 ^ 2 + t55 ^ 2) + t141 * t2 + t139 * t1 + t152 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t80(1), t80(2), t80(4), t80(7), t80(11); t80(2), t80(3), t80(5), t80(8), t80(12); t80(4), t80(5), t80(6), t80(9), t80(13); t80(7), t80(8), t80(9), t80(10), t80(14); t80(11), t80(12), t80(13), t80(14), t80(15);];
Mq = res;
