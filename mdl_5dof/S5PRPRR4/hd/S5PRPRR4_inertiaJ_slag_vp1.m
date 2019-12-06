% Calculate joint inertia matrix for
% S5PRPRR4
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:41
% EndTime: 2019-12-05 15:49:57
% DurationCPUTime: 3.55s
% Computational Cost: add. (11923->378), mult. (30609->582), div. (0->0), fcn. (40176->12), ass. (0->172)
t151 = sin(pkin(9));
t153 = cos(pkin(9));
t154 = cos(pkin(5));
t186 = cos(qJ(2));
t168 = t154 * t186;
t184 = sin(qJ(2));
t139 = -t151 * t184 + t153 * t168;
t167 = t154 * t184;
t140 = t151 * t186 + t153 * t167;
t152 = sin(pkin(5));
t174 = t152 * t153;
t112 = Icges(3,5) * t140 + Icges(3,6) * t139 - Icges(3,3) * t174;
t176 = sin(pkin(10));
t177 = cos(pkin(10));
t144 = -t184 * t176 + t186 * t177;
t158 = t154 * t144;
t159 = t186 * t176 + t184 * t177;
t123 = -t151 * t159 + t153 * t158;
t136 = t159 * t154;
t124 = t153 * t136 + t151 * t144;
t79 = Icges(4,5) * t124 + Icges(4,6) * t123 - Icges(4,3) * t174;
t194 = -t112 - t79;
t141 = -t151 * t168 - t153 * t184;
t142 = -t151 * t167 + t153 * t186;
t175 = t151 * t152;
t113 = Icges(3,5) * t142 + Icges(3,6) * t141 + Icges(3,3) * t175;
t125 = -t151 * t158 - t153 * t159;
t126 = -t151 * t136 + t153 * t144;
t80 = Icges(4,5) * t126 + Icges(4,6) * t125 + Icges(4,3) * t175;
t193 = t113 + t80;
t134 = t144 * t152;
t135 = t159 * t152;
t192 = Icges(4,5) * t135 + Icges(4,6) * t134 + (t184 * Icges(3,5) + t186 * Icges(3,6)) * t152 + (Icges(4,3) + Icges(3,3)) * t154;
t156 = sin(qJ(4));
t185 = cos(qJ(4));
t169 = t152 * t185;
t101 = t124 * t156 + t153 * t169;
t103 = t126 * t156 - t151 * t169;
t127 = t135 * t156 - t154 * t185;
t173 = t152 * t156;
t102 = t124 * t185 - t153 * t173;
t155 = sin(qJ(5));
t157 = cos(qJ(5));
t75 = -t102 * t155 - t123 * t157;
t76 = t102 * t157 - t123 * t155;
t50 = Icges(6,5) * t76 + Icges(6,6) * t75 + Icges(6,3) * t101;
t52 = Icges(6,4) * t76 + Icges(6,2) * t75 + Icges(6,6) * t101;
t54 = Icges(6,1) * t76 + Icges(6,4) * t75 + Icges(6,5) * t101;
t17 = t101 * t50 + t75 * t52 + t76 * t54;
t104 = t126 * t185 + t151 * t173;
t77 = -t104 * t155 - t125 * t157;
t78 = t104 * t157 - t125 * t155;
t51 = Icges(6,5) * t78 + Icges(6,6) * t77 + Icges(6,3) * t103;
t53 = Icges(6,4) * t78 + Icges(6,2) * t77 + Icges(6,6) * t103;
t55 = Icges(6,1) * t78 + Icges(6,4) * t77 + Icges(6,5) * t103;
t18 = t101 * t51 + t75 * t53 + t76 * t55;
t128 = t135 * t185 + t154 * t156;
t100 = t128 * t157 - t134 * t155;
t99 = -t128 * t155 - t134 * t157;
t68 = Icges(6,5) * t100 + Icges(6,6) * t99 + Icges(6,3) * t127;
t69 = Icges(6,4) * t100 + Icges(6,2) * t99 + Icges(6,6) * t127;
t70 = Icges(6,1) * t100 + Icges(6,4) * t99 + Icges(6,5) * t127;
t28 = t101 * t68 + t75 * t69 + t76 * t70;
t1 = t17 * t101 + t18 * t103 + t28 * t127;
t191 = -t1 / 0.2e1;
t21 = t100 * t54 + t127 * t50 + t99 * t52;
t22 = t100 * t55 + t127 * t51 + t99 * t53;
t36 = t100 * t70 + t127 * t68 + t99 * t69;
t7 = t21 * t101 + t22 * t103 + t36 * t127;
t190 = t7 / 0.2e1;
t189 = t101 / 0.2e1;
t188 = t103 / 0.2e1;
t187 = t127 / 0.2e1;
t183 = t186 * pkin(2);
t56 = t76 * rSges(6,1) + t75 * rSges(6,2) + t101 * rSges(6,3);
t182 = t102 * pkin(4) + t101 * pkin(8) + t56;
t57 = t78 * rSges(6,1) + t77 * rSges(6,2) + t103 * rSges(6,3);
t181 = t104 * pkin(4) + t103 * pkin(8) + t57;
t71 = t100 * rSges(6,1) + t99 * rSges(6,2) + t127 * rSges(6,3);
t180 = t128 * pkin(4) + t127 * pkin(8) + t71;
t166 = pkin(2) * t167 - qJ(3) * t152;
t121 = -t166 * t151 + t183 * t153;
t111 = t154 * t121;
t91 = t126 * pkin(3) - t125 * pkin(7);
t179 = t154 * t91 + t111;
t120 = t183 * t151 + t166 * t153;
t90 = t124 * pkin(3) - t123 * pkin(7);
t178 = -t120 - t90;
t172 = t120 * t175 + t121 * t174;
t145 = t152 * t184 * pkin(2) + t154 * qJ(3);
t171 = -t135 * pkin(3) + t134 * pkin(7) - t145;
t170 = m(4) + m(5) + m(6);
t165 = (-t135 * rSges(4,1) - t134 * rSges(4,2) - t154 * rSges(4,3) - t145) * t152;
t164 = t91 * t174 + t90 * t175 + t172;
t95 = t128 * rSges(5,1) - t127 * rSges(5,2) - t134 * rSges(5,3);
t163 = (-t95 + t171) * t152;
t160 = (t171 - t180) * t152;
t133 = t154 * rSges(3,3) + (t184 * rSges(3,1) + t186 * rSges(3,2)) * t152;
t132 = Icges(3,5) * t154 + (t184 * Icges(3,1) + t186 * Icges(3,4)) * t152;
t131 = Icges(3,6) * t154 + (t184 * Icges(3,4) + t186 * Icges(3,2)) * t152;
t119 = t142 * rSges(3,1) + t141 * rSges(3,2) + rSges(3,3) * t175;
t118 = t140 * rSges(3,1) + t139 * rSges(3,2) - rSges(3,3) * t174;
t117 = Icges(3,1) * t142 + Icges(3,4) * t141 + Icges(3,5) * t175;
t116 = Icges(3,1) * t140 + Icges(3,4) * t139 - Icges(3,5) * t174;
t115 = Icges(3,4) * t142 + Icges(3,2) * t141 + Icges(3,6) * t175;
t114 = Icges(3,4) * t140 + Icges(3,2) * t139 - Icges(3,6) * t174;
t109 = Icges(4,1) * t135 + Icges(4,4) * t134 + Icges(4,5) * t154;
t108 = Icges(4,4) * t135 + Icges(4,2) * t134 + Icges(4,6) * t154;
t97 = -t154 * t118 - t133 * t174;
t96 = t154 * t119 - t133 * t175;
t94 = Icges(5,1) * t128 - Icges(5,4) * t127 - Icges(5,5) * t134;
t93 = Icges(5,4) * t128 - Icges(5,2) * t127 - Icges(5,6) * t134;
t92 = Icges(5,5) * t128 - Icges(5,6) * t127 - Icges(5,3) * t134;
t86 = t126 * rSges(4,1) + t125 * rSges(4,2) + rSges(4,3) * t175;
t85 = t124 * rSges(4,1) + t123 * rSges(4,2) - rSges(4,3) * t174;
t84 = Icges(4,1) * t126 + Icges(4,4) * t125 + Icges(4,5) * t175;
t83 = Icges(4,1) * t124 + Icges(4,4) * t123 - Icges(4,5) * t174;
t82 = Icges(4,4) * t126 + Icges(4,2) * t125 + Icges(4,6) * t175;
t81 = Icges(4,4) * t124 + Icges(4,2) * t123 - Icges(4,6) * t174;
t74 = (t118 * t151 + t119 * t153) * t152;
t67 = t104 * rSges(5,1) - t103 * rSges(5,2) - t125 * rSges(5,3);
t66 = t102 * rSges(5,1) - t101 * rSges(5,2) - t123 * rSges(5,3);
t65 = Icges(5,1) * t104 - Icges(5,4) * t103 - Icges(5,5) * t125;
t64 = Icges(5,1) * t102 - Icges(5,4) * t101 - Icges(5,5) * t123;
t63 = Icges(5,4) * t104 - Icges(5,2) * t103 - Icges(5,6) * t125;
t62 = Icges(5,4) * t102 - Icges(5,2) * t101 - Icges(5,6) * t123;
t61 = Icges(5,5) * t104 - Icges(5,6) * t103 - Icges(5,3) * t125;
t60 = Icges(5,5) * t102 - Icges(5,6) * t101 - Icges(5,3) * t123;
t59 = (-t120 - t85) * t154 + t153 * t165;
t58 = t151 * t165 + t154 * t86 + t111;
t49 = (t151 * t85 + t153 * t86) * t152 + t172;
t48 = t125 * t95 - t134 * t67;
t47 = -t123 * t95 + t134 * t66;
t46 = -t127 * t93 + t128 * t94 - t134 * t92;
t45 = t123 * t67 - t125 * t66;
t44 = -t103 * t93 + t104 * t94 - t125 * t92;
t43 = -t101 * t93 + t102 * t94 - t123 * t92;
t42 = (-t66 + t178) * t154 + t153 * t163;
t41 = t151 * t163 + t154 * t67 + t179;
t40 = -t103 * t71 + t127 * t57;
t39 = t101 * t71 - t127 * t56;
t38 = -t127 * t63 + t128 * t65 - t134 * t61;
t37 = -t127 * t62 + t128 * t64 - t134 * t60;
t35 = (t151 * t66 + t153 * t67) * t152 + t164;
t34 = -t103 * t63 + t104 * t65 - t125 * t61;
t33 = -t103 * t62 + t104 * t64 - t125 * t60;
t32 = -t101 * t63 + t102 * t65 - t123 * t61;
t31 = -t101 * t62 + t102 * t64 - t123 * t60;
t30 = -t101 * t57 + t103 * t56;
t29 = t103 * t68 + t77 * t69 + t78 * t70;
t27 = t180 * t125 - t181 * t134;
t26 = -t180 * t123 + t182 * t134;
t25 = (t178 - t182) * t154 + t153 * t160;
t24 = t151 * t160 + t181 * t154 + t179;
t23 = t181 * t123 - t182 * t125;
t20 = t103 * t51 + t77 * t53 + t78 * t55;
t19 = t103 * t50 + t77 * t52 + t78 * t54;
t16 = (t182 * t151 + t181 * t153) * t152 + t164;
t15 = t46 * t154 + (t151 * t38 - t153 * t37) * t152;
t14 = -t37 * t123 - t38 * t125 - t46 * t134;
t13 = t44 * t154 + (t151 * t34 - t153 * t33) * t152;
t12 = t43 * t154 + (t151 * t32 - t153 * t31) * t152;
t11 = -t33 * t123 - t34 * t125 - t44 * t134;
t10 = -t31 * t123 - t32 * t125 - t43 * t134;
t9 = t36 * t154 + (t151 * t22 - t153 * t21) * t152;
t8 = -t21 * t123 - t22 * t125 - t36 * t134;
t6 = t29 * t154 + (t151 * t20 - t153 * t19) * t152;
t5 = t28 * t154 + (t151 * t18 - t153 * t17) * t152;
t4 = -t19 * t123 - t20 * t125 - t29 * t134;
t3 = -t17 * t123 - t18 * t125 - t28 * t134;
t2 = t19 * t101 + t20 * t103 + t29 * t127;
t72 = [m(2) + m(3) + t170; m(3) * t74 + m(4) * t49 + m(5) * t35 + m(6) * t16; m(6) * (t16 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t35 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(4) * (t49 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(3) * (t74 ^ 2 + t96 ^ 2 + t97 ^ 2) + (t9 + t15 + (t134 * t108 + t135 * t109 + t192 * t154) * t154 + ((t134 * t82 + t135 * t84) * t151 - (t134 * t81 + t135 * t83) * t153 + (t186 * t131 + t184 * t132 + t80 * t151 - t79 * t153) * t154) * t152) * t154 + (t6 + t13 + (t141 * t115 + t142 * t117 + t125 * t82 + t126 * t84 + t193 * t175) * t175 + (t141 * t131 + t142 * t132 + t125 * t108 + t126 * t109 + (t186 * t115 + t184 * t117) * t152 + t192 * t175 + t154 * t113) * t154) * t175 + (-t5 - t12 + (t139 * t114 + t140 * t116 + t123 * t81 + t124 * t83 + t194 * t174) * t174 + (-t139 * t131 - t140 * t132 - t123 * t108 - t124 * t109 - (t186 * t114 + t184 * t116) * t152 + t192 * t174 - t154 * t112) * t154 + (-t141 * t114 - t139 * t115 - t142 * t116 - t140 * t117 - t123 * t82 - t124 * t84 - t125 * t81 - t126 * t83 + t193 * t174 + t194 * t175) * t175) * t174; t170 * t154; m(6) * (t154 * t16 + (t151 * t25 - t153 * t24) * t152) + m(5) * (t154 * t35 + (t151 * t42 - t153 * t41) * t152) + m(4) * (t154 * t49 + (t151 * t59 - t153 * t58) * t152); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t154 ^ 2 + (t151 ^ 2 + t153 ^ 2) * t152 ^ 2); m(5) * t45 + m(6) * t23; (t8 / 0.2e1 + t14 / 0.2e1) * t154 - (t9 / 0.2e1 + t15 / 0.2e1) * t134 + (-t6 / 0.2e1 - t13 / 0.2e1) * t125 + (-t5 / 0.2e1 - t12 / 0.2e1) * t123 + m(6) * (t16 * t23 + t24 * t27 + t25 * t26) + m(5) * (t45 * t35 + t48 * t41 + t47 * t42) + ((-t3 / 0.2e1 - t10 / 0.2e1) * t153 + (t4 / 0.2e1 + t11 / 0.2e1) * t151) * t152; m(5) * (t45 * t154 + (t151 * t47 - t153 * t48) * t152) + m(6) * (t23 * t154 + (t151 * t26 - t153 * t27) * t152); -(t8 + t14) * t134 + (-t4 - t11) * t125 + (-t3 - t10) * t123 + m(6) * (t23 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(5) * (t45 ^ 2 + t47 ^ 2 + t48 ^ 2); m(6) * t30; t6 * t188 + t5 * t189 + t9 * t187 + t154 * t190 + m(6) * (t30 * t16 + t40 * t24 + t39 * t25) + (t151 * t2 / 0.2e1 + t153 * t191) * t152; m(6) * (t30 * t154 + (t151 * t39 - t153 * t40) * t152); m(6) * (t30 * t23 + t39 * t26 + t40 * t27) - t134 * t190 + t123 * t191 + t3 * t189 + t8 * t187 + t4 * t188 - t125 * t2 / 0.2e1; m(6) * (t30 ^ 2 + t39 ^ 2 + t40 ^ 2) + t103 * t2 + t101 * t1 + t127 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t72(1), t72(2), t72(4), t72(7), t72(11); t72(2), t72(3), t72(5), t72(8), t72(12); t72(4), t72(5), t72(6), t72(9), t72(13); t72(7), t72(8), t72(9), t72(10), t72(14); t72(11), t72(12), t72(13), t72(14), t72(15);];
Mq = res;
