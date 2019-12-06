% Calculate joint inertia matrix for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:29
% EndTime: 2019-12-05 15:56:36
% DurationCPUTime: 2.48s
% Computational Cost: add. (10613->384), mult. (18150->586), div. (0->0), fcn. (22910->12), ass. (0->172)
t148 = sin(pkin(5));
t188 = t148 ^ 2;
t151 = cos(pkin(5));
t154 = sin(qJ(2));
t156 = cos(qJ(2));
t129 = Icges(3,6) * t151 + (Icges(3,4) * t154 + Icges(3,2) * t156) * t148;
t146 = sin(pkin(10));
t149 = cos(pkin(10));
t173 = t148 * t154;
t134 = -t146 * t173 + t149 * t151;
t176 = t146 * t151;
t135 = t149 * t173 + t176;
t172 = t148 * t156;
t98 = Icges(4,5) * t135 + Icges(4,6) * t134 - Icges(4,3) * t172;
t187 = t129 - t98;
t147 = sin(pkin(9));
t150 = cos(pkin(9));
t171 = t154 * t151;
t137 = t156 * t147 + t150 * t171;
t165 = pkin(10) + qJ(4);
t145 = sin(t165);
t162 = cos(t165);
t158 = t148 * t162;
t118 = t137 * t145 + t150 * t158;
t139 = -t147 * t171 + t156 * t150;
t120 = t139 * t145 - t147 * t158;
t131 = t145 * t173 - t151 * t162;
t174 = t148 * t150;
t119 = t137 * t162 - t145 * t174;
t170 = t156 * t151;
t136 = t154 * t147 - t150 * t170;
t153 = sin(qJ(5));
t155 = cos(qJ(5));
t90 = -t119 * t153 + t136 * t155;
t91 = t119 * t155 + t136 * t153;
t52 = Icges(6,5) * t91 + Icges(6,6) * t90 + Icges(6,3) * t118;
t54 = Icges(6,4) * t91 + Icges(6,2) * t90 + Icges(6,6) * t118;
t56 = Icges(6,1) * t91 + Icges(6,4) * t90 + Icges(6,5) * t118;
t175 = t147 * t148;
t121 = t139 * t162 + t145 * t175;
t138 = t147 * t170 + t150 * t154;
t92 = -t121 * t153 + t138 * t155;
t93 = t121 * t155 + t138 * t153;
t19 = t120 * t52 + t54 * t92 + t56 * t93;
t53 = Icges(6,5) * t93 + Icges(6,6) * t92 + Icges(6,3) * t120;
t55 = Icges(6,4) * t93 + Icges(6,2) * t92 + Icges(6,6) * t120;
t57 = Icges(6,1) * t93 + Icges(6,4) * t92 + Icges(6,5) * t120;
t20 = t120 * t53 + t55 * t92 + t57 * t93;
t132 = t151 * t145 + t154 * t158;
t122 = -t132 * t153 - t155 * t172;
t123 = t132 * t155 - t153 * t172;
t63 = Icges(6,5) * t123 + Icges(6,6) * t122 + Icges(6,3) * t131;
t64 = Icges(6,4) * t123 + Icges(6,2) * t122 + Icges(6,6) * t131;
t65 = Icges(6,1) * t123 + Icges(6,4) * t122 + Icges(6,5) * t131;
t28 = t120 * t63 + t64 * t92 + t65 * t93;
t2 = t118 * t19 + t120 * t20 + t131 * t28;
t186 = t2 / 0.2e1;
t185 = t118 / 0.2e1;
t184 = t120 / 0.2e1;
t183 = t131 / 0.2e1;
t182 = pkin(3) * t149;
t58 = rSges(6,1) * t91 + rSges(6,2) * t90 + rSges(6,3) * t118;
t181 = pkin(4) * t119 + pkin(8) * t118 + t58;
t59 = rSges(6,1) * t93 + rSges(6,2) * t92 + rSges(6,3) * t120;
t180 = pkin(4) * t121 + pkin(8) * t120 + t59;
t67 = rSges(6,1) * t123 + rSges(6,2) * t122 + rSges(6,3) * t131;
t179 = pkin(4) * t132 + pkin(8) * t131 + t67;
t117 = pkin(2) * t139 + qJ(3) * t138;
t115 = t151 * t117;
t164 = t146 * t175;
t77 = pkin(3) * t164 + pkin(7) * t138 + t139 * t182;
t178 = t151 * t77 + t115;
t116 = pkin(2) * t137 + qJ(3) * t136;
t163 = t146 * t174;
t76 = -pkin(3) * t163 + pkin(7) * t136 + t137 * t182;
t177 = -t116 - t76;
t140 = (pkin(2) * t154 - qJ(3) * t156) * t148;
t168 = -pkin(3) * t176 - (-pkin(7) * t156 + t154 * t182) * t148 - t140;
t167 = t116 * t175 + t117 * t174;
t166 = -m(4) - m(5) - m(6);
t161 = (-t135 * rSges(4,1) - t134 * rSges(4,2) + rSges(4,3) * t172 - t140) * t148;
t160 = t77 * t174 + t76 * t175 + t167;
t97 = t132 * rSges(5,1) - t131 * rSges(5,2) - rSges(5,3) * t172;
t159 = (-t97 + t168) * t148;
t157 = (t168 - t179) * t148;
t133 = t151 * rSges(3,3) + (rSges(3,1) * t154 + rSges(3,2) * t156) * t148;
t130 = Icges(3,5) * t151 + (Icges(3,1) * t154 + Icges(3,4) * t156) * t148;
t128 = Icges(3,3) * t151 + (Icges(3,5) * t154 + Icges(3,6) * t156) * t148;
t127 = t139 * t149 + t164;
t126 = -t139 * t146 + t149 * t175;
t125 = t137 * t149 - t163;
t124 = -t137 * t146 - t149 * t174;
t111 = rSges(3,1) * t139 - rSges(3,2) * t138 + rSges(3,3) * t175;
t110 = rSges(3,1) * t137 - rSges(3,2) * t136 - rSges(3,3) * t174;
t106 = Icges(3,1) * t139 - Icges(3,4) * t138 + Icges(3,5) * t175;
t105 = Icges(3,1) * t137 - Icges(3,4) * t136 - Icges(3,5) * t174;
t104 = Icges(3,4) * t139 - Icges(3,2) * t138 + Icges(3,6) * t175;
t103 = Icges(3,4) * t137 - Icges(3,2) * t136 - Icges(3,6) * t174;
t102 = Icges(3,5) * t139 - Icges(3,6) * t138 + Icges(3,3) * t175;
t101 = Icges(3,5) * t137 - Icges(3,6) * t136 - Icges(3,3) * t174;
t100 = Icges(4,1) * t135 + Icges(4,4) * t134 - Icges(4,5) * t172;
t99 = Icges(4,4) * t135 + Icges(4,2) * t134 - Icges(4,6) * t172;
t96 = Icges(5,1) * t132 - Icges(5,4) * t131 - Icges(5,5) * t172;
t95 = Icges(5,4) * t132 - Icges(5,2) * t131 - Icges(5,6) * t172;
t94 = Icges(5,5) * t132 - Icges(5,6) * t131 - Icges(5,3) * t172;
t87 = -t110 * t151 - t133 * t174;
t86 = t111 * t151 - t133 * t175;
t85 = rSges(4,1) * t127 + rSges(4,2) * t126 + rSges(4,3) * t138;
t84 = rSges(4,1) * t125 + rSges(4,2) * t124 + rSges(4,3) * t136;
t83 = Icges(4,1) * t127 + Icges(4,4) * t126 + Icges(4,5) * t138;
t82 = Icges(4,1) * t125 + Icges(4,4) * t124 + Icges(4,5) * t136;
t81 = Icges(4,4) * t127 + Icges(4,2) * t126 + Icges(4,6) * t138;
t80 = Icges(4,4) * t125 + Icges(4,2) * t124 + Icges(4,6) * t136;
t79 = Icges(4,5) * t127 + Icges(4,6) * t126 + Icges(4,3) * t138;
t78 = Icges(4,5) * t125 + Icges(4,6) * t124 + Icges(4,3) * t136;
t75 = rSges(5,1) * t121 - rSges(5,2) * t120 + rSges(5,3) * t138;
t74 = rSges(5,1) * t119 - rSges(5,2) * t118 + rSges(5,3) * t136;
t73 = Icges(5,1) * t121 - Icges(5,4) * t120 + Icges(5,5) * t138;
t72 = Icges(5,1) * t119 - Icges(5,4) * t118 + Icges(5,5) * t136;
t71 = Icges(5,4) * t121 - Icges(5,2) * t120 + Icges(5,6) * t138;
t70 = Icges(5,4) * t119 - Icges(5,2) * t118 + Icges(5,6) * t136;
t69 = Icges(5,5) * t121 - Icges(5,6) * t120 + Icges(5,3) * t138;
t68 = Icges(5,5) * t119 - Icges(5,6) * t118 + Icges(5,3) * t136;
t60 = (t110 * t147 + t111 * t150) * t148;
t51 = -t138 * t97 - t172 * t75;
t50 = t136 * t97 + t172 * t74;
t49 = (-t116 - t84) * t151 + t150 * t161;
t48 = t147 * t161 + t151 * t85 + t115;
t47 = -t131 * t95 + t132 * t96 - t172 * t94;
t46 = -t136 * t75 + t138 * t74;
t45 = -t120 * t95 + t121 * t96 + t138 * t94;
t44 = -t118 * t95 + t119 * t96 + t136 * t94;
t43 = (t147 * t84 + t150 * t85) * t148 + t167;
t42 = -t120 * t67 + t131 * t59;
t41 = t118 * t67 - t131 * t58;
t40 = -t131 * t71 + t132 * t73 - t172 * t69;
t39 = -t131 * t70 + t132 * t72 - t172 * t68;
t38 = -t120 * t71 + t121 * t73 + t138 * t69;
t37 = -t120 * t70 + t121 * t72 + t138 * t68;
t36 = -t118 * t71 + t119 * t73 + t136 * t69;
t35 = -t118 * t70 + t119 * t72 + t136 * t68;
t34 = (-t74 + t177) * t151 + t150 * t159;
t33 = t147 * t159 + t151 * t75 + t178;
t32 = t122 * t64 + t123 * t65 + t131 * t63;
t31 = -t118 * t59 + t120 * t58;
t30 = -t138 * t179 - t172 * t180;
t29 = t136 * t179 + t172 * t181;
t27 = t118 * t63 + t64 * t90 + t65 * t91;
t26 = (t147 * t74 + t150 * t75) * t148 + t160;
t25 = -t136 * t180 + t138 * t181;
t24 = t122 * t55 + t123 * t57 + t131 * t53;
t23 = t122 * t54 + t123 * t56 + t131 * t52;
t22 = (t177 - t181) * t151 + t150 * t157;
t21 = t147 * t157 + t151 * t180 + t178;
t18 = t118 * t53 + t55 * t90 + t57 * t91;
t17 = t118 * t52 + t54 * t90 + t56 * t91;
t16 = (t147 * t181 + t150 * t180) * t148 + t160;
t15 = t151 * t47 + (t147 * t40 - t150 * t39) * t148;
t14 = t39 * t136 + t40 * t138 - t172 * t47;
t13 = t151 * t45 + (t147 * t38 - t150 * t37) * t148;
t12 = t151 * t44 + (t147 * t36 - t150 * t35) * t148;
t11 = t37 * t136 + t38 * t138 - t172 * t45;
t10 = t35 * t136 + t36 * t138 - t172 * t44;
t9 = t151 * t32 + (t147 * t24 - t150 * t23) * t148;
t8 = t23 * t136 + t24 * t138 - t172 * t32;
t7 = t118 * t23 + t120 * t24 + t131 * t32;
t6 = t151 * t28 + (t147 * t20 - t150 * t19) * t148;
t5 = t151 * t27 + (t147 * t18 - t150 * t17) * t148;
t4 = t19 * t136 + t20 * t138 - t172 * t28;
t3 = t17 * t136 + t18 * t138 - t172 * t27;
t1 = t118 * t17 + t120 * t18 + t131 * t27;
t61 = [m(2) + m(3) - t166; m(3) * t60 + m(4) * t43 + m(5) * t26 + m(6) * t16; m(6) * (t16 ^ 2 + t21 ^ 2 + t22 ^ 2) + m(5) * (t26 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(4) * (t43 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(3) * (t60 ^ 2 + t86 ^ 2 + t87 ^ 2) + (t6 + t13 + ((t126 * t81 + t127 * t83 + t138 * t79) * t147 - (t126 * t80 + t127 * t82 + t138 * t78) * t150) * t148 + (t102 * t175 - t104 * t138 + t106 * t139) * t175) * t175 + (-t5 - t12 - ((t124 * t81 + t125 * t83 + t136 * t79) * t147 - (t124 * t80 + t125 * t82 + t136 * t78) * t150) * t148 + (-t101 * t174 - t103 * t136 + t105 * t137) * t174 + (-t101 * t175 + t102 * t174 + t103 * t138 + t104 * t136 - t105 * t139 - t106 * t137) * t175) * t174 + (t9 + t15 + ((t104 * t156 + t106 * t154) * t147 - (t103 * t156 + t105 * t154) * t150) * t188 + ((-t101 * t150 + t102 * t147 + t129 * t156 + t130 * t154) * t148 + t135 * t100 + t134 * t99 - t98 * t172 + t151 * t128) * t151 + (t100 * t127 + t126 * t99 + t128 * t175 + t130 * t139 + t134 * t81 + t135 * t83 - t138 * t187 - t172 * t79) * t175 + (-t100 * t125 - t124 * t99 + t128 * t174 - t130 * t137 - t134 * t80 - t135 * t82 + t136 * t187 + t172 * t78) * t174) * t151; t166 * t172; m(6) * (t136 * t21 + t138 * t22 - t16 * t172) + m(5) * (t136 * t33 + t138 * t34 - t172 * t26) + m(4) * (t136 * t48 + t138 * t49 - t172 * t43); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t188 * t156 ^ 2 + t136 ^ 2 + t138 ^ 2); m(5) * t46 + m(6) * t25; (t8 / 0.2e1 + t14 / 0.2e1) * t151 + (t6 / 0.2e1 + t13 / 0.2e1) * t138 + (t5 / 0.2e1 + t12 / 0.2e1) * t136 + m(6) * (t16 * t25 + t21 * t30 + t22 * t29) + m(5) * (t26 * t46 + t33 * t51 + t34 * t50) + ((-t9 / 0.2e1 - t15 / 0.2e1) * t156 + (-t3 / 0.2e1 - t10 / 0.2e1) * t150 + (t4 / 0.2e1 + t11 / 0.2e1) * t147) * t148; m(5) * (t51 * t136 + t50 * t138 - t172 * t46) + m(6) * (t30 * t136 + t29 * t138 - t172 * t25); (-t14 - t8) * t172 + (t4 + t11) * t138 + (t3 + t10) * t136 + m(6) * (t25 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(5) * (t46 ^ 2 + t50 ^ 2 + t51 ^ 2); m(6) * t31; t151 * t7 / 0.2e1 + m(6) * (t16 * t31 + t21 * t42 + t22 * t41) + t6 * t184 + t9 * t183 + t5 * t185 + (t147 * t186 - t150 * t1 / 0.2e1) * t148; m(6) * (t42 * t136 + t41 * t138 - t172 * t31); m(6) * (t25 * t31 + t29 * t41 + t30 * t42) + t4 * t184 + t8 * t183 + t3 * t185 + t138 * t186 + t136 * t1 / 0.2e1 - t7 * t172 / 0.2e1; m(6) * (t31 ^ 2 + t41 ^ 2 + t42 ^ 2) + t120 * t2 + t118 * t1 + t131 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t61(1), t61(2), t61(4), t61(7), t61(11); t61(2), t61(3), t61(5), t61(8), t61(12); t61(4), t61(5), t61(6), t61(9), t61(13); t61(7), t61(8), t61(9), t61(10), t61(14); t61(11), t61(12), t61(13), t61(14), t61(15);];
Mq = res;
