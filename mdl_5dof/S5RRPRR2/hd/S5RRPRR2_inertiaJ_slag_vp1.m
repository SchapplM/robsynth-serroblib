% Calculate joint inertia matrix for
% S5RRPRR2
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:05
% EndTime: 2019-12-05 18:27:10
% DurationCPUTime: 1.60s
% Computational Cost: add. (3257->252), mult. (2558->363), div. (0->0), fcn. (2324->10), ass. (0->131)
t196 = Icges(3,3) + Icges(4,3);
t119 = qJ(2) + pkin(9);
t108 = sin(t119);
t109 = cos(t119);
t123 = sin(qJ(2));
t125 = cos(qJ(2));
t195 = Icges(3,5) * t125 + Icges(4,5) * t109 - Icges(3,6) * t123 - Icges(4,6) * t108;
t124 = sin(qJ(1));
t120 = t124 ^ 2;
t194 = t124 * pkin(6);
t110 = qJ(4) + t119;
t104 = sin(t110);
t105 = cos(t110);
t153 = rSges(5,1) * t105 - rSges(5,2) * t104;
t154 = rSges(4,1) * t109 - rSges(4,2) * t108;
t106 = qJ(5) + t110;
t100 = cos(t106);
t99 = sin(t106);
t159 = rSges(6,1) * t100 - rSges(6,2) * t99;
t126 = cos(qJ(1));
t193 = -t195 * t124 + t196 * t126;
t192 = t196 * t124 + t195 * t126;
t121 = t126 ^ 2;
t191 = t124 / 0.2e1;
t190 = -t126 / 0.2e1;
t188 = pkin(2) * t123;
t122 = -qJ(3) - pkin(6);
t107 = t125 * pkin(2) + pkin(1);
t132 = t124 * rSges(6,3) + t126 * t159;
t17 = t124 * (-t126 * rSges(6,3) + t159 * t124) + t126 * t132;
t130 = t124 * rSges(5,3) + t126 * t153;
t20 = t124 * (-t126 * rSges(5,3) + t153 * t124) + t126 * t130;
t117 = t126 * pkin(6);
t98 = t126 * t107;
t187 = t124 * (t117 + (-pkin(1) + t107) * t124) + t126 * (-t126 * pkin(1) - t194 + t98);
t186 = rSges(3,1) * t125;
t182 = rSges(3,2) * t123;
t179 = Icges(6,4) * t99;
t178 = t126 * rSges(3,3);
t177 = Icges(3,4) * t123;
t176 = Icges(3,4) * t125;
t175 = Icges(4,4) * t108;
t174 = Icges(4,4) * t109;
t173 = Icges(5,4) * t104;
t172 = Icges(5,4) * t105;
t171 = Icges(6,4) * t100;
t170 = t124 * rSges(3,3) + t126 * t186;
t167 = t120 + t121;
t118 = -pkin(7) + t122;
t87 = pkin(3) * t109 + t107;
t143 = -Icges(6,2) * t99 + t171;
t43 = Icges(6,6) * t124 + t143 * t126;
t144 = Icges(6,1) * t100 - t179;
t45 = Icges(6,5) * t124 + t144 * t126;
t157 = t100 * t45 - t43 * t99;
t42 = -Icges(6,6) * t126 + t143 * t124;
t44 = -Icges(6,5) * t126 + t144 * t124;
t158 = -t100 * t44 + t42 * t99;
t142 = Icges(6,5) * t100 - Icges(6,6) * t99;
t40 = -Icges(6,3) * t126 + t142 * t124;
t41 = Icges(6,3) * t124 + t142 * t126;
t2 = t124 * (t120 * t41 + (t158 * t126 + (t157 - t40) * t124) * t126);
t3 = t121 * t40 + (t157 * t124 + (t158 - t41) * t126) * t124;
t166 = -t126 * t3 + t2;
t165 = -t108 * rSges(4,1) - t109 * rSges(4,2) - t188;
t76 = t99 * rSges(6,1) + t100 * rSges(6,2);
t164 = -pkin(4) * t104 - t76;
t74 = Icges(6,2) * t100 + t179;
t75 = Icges(6,1) * t99 + t171;
t156 = t100 * t75 - t74 * t99;
t73 = Icges(6,5) * t99 + Icges(6,6) * t100;
t163 = (t100 * t43 + t124 * t73 + t156 * t126 + t45 * t99) * t191 + (t100 * t42 + t156 * t124 - t126 * t73 + t44 * t99) * t190;
t65 = pkin(4) * t105 + t87;
t64 = t126 * t65;
t82 = t126 * t87;
t7 = t126 * (t64 - t82) + t17 + (t65 - t87) * t120;
t162 = t126 * (t82 - t98) + t187 + (-t107 + t87) * t120;
t160 = -pkin(3) * t108 - t188;
t155 = -t182 + t186;
t136 = -Icges(5,2) * t104 + t172;
t52 = -Icges(5,6) * t126 + t136 * t124;
t139 = Icges(5,1) * t105 - t173;
t54 = -Icges(5,5) * t126 + t139 * t124;
t152 = t104 * t52 - t105 * t54;
t53 = Icges(5,6) * t124 + t136 * t126;
t55 = Icges(5,5) * t124 + t139 * t126;
t151 = -t104 * t53 + t105 * t55;
t79 = Icges(5,2) * t105 + t173;
t80 = Icges(5,1) * t104 + t172;
t150 = -t104 * t79 + t105 * t80;
t133 = Icges(5,5) * t105 - Icges(5,6) * t104;
t50 = -Icges(5,3) * t126 + t133 * t124;
t51 = Icges(5,3) * t124 + t133 * t126;
t4 = t124 * (t120 * t51 + (t152 * t126 + (t151 - t50) * t124) * t126);
t5 = t121 * t50 + (t151 * t124 + (t152 - t51) * t126) * t124;
t145 = t2 + t4 + (-t5 - t3) * t126;
t141 = Icges(3,1) * t125 - t177;
t140 = Icges(4,1) * t109 - t175;
t138 = -Icges(3,2) * t123 + t176;
t137 = -Icges(4,2) * t108 + t174;
t131 = t124 * rSges(4,3) + t126 * t154;
t81 = t104 * rSges(5,1) + t105 * rSges(5,2);
t129 = t160 - t81;
t78 = Icges(5,5) * t104 + Icges(5,6) * t105;
t128 = t163 + (t104 * t55 + t105 * t53 + t124 * t78 + t150 * t126) * t191 + (t104 * t54 + t105 * t52 + t150 * t124 - t126 * t78) * t190;
t127 = t160 + t164;
t111 = -pkin(8) + t118;
t96 = t126 * rSges(2,1) - t124 * rSges(2,2);
t95 = -t124 * rSges(2,1) - t126 * rSges(2,2);
t94 = t123 * rSges(3,1) + t125 * rSges(3,2);
t57 = t165 * t126;
t56 = t165 * t124;
t49 = t194 + (pkin(1) - t182) * t126 + t170;
t48 = t178 + t117 + (-pkin(1) - t155) * t124;
t35 = t164 * t126;
t34 = t164 * t124;
t33 = t129 * t126;
t32 = t129 * t124;
t31 = -t124 * t122 + t131 + t98;
t30 = (rSges(4,3) - t122) * t126 + (-t107 - t154) * t124;
t27 = t126 * (-t126 * t182 + t170) + (t155 * t124 - t178) * t124;
t26 = -t124 * t118 + t130 + t82;
t25 = (rSges(5,3) - t118) * t126 + (-t153 - t87) * t124;
t24 = t127 * t126;
t23 = t127 * t124;
t19 = -t124 * t111 + t132 + t64;
t18 = (rSges(6,3) - t111) * t126 + (-t159 - t65) * t124;
t8 = t126 * t131 + (-t126 * rSges(4,3) + t154 * t124) * t124 + t187;
t6 = t162 + t20;
t1 = t162 + t7;
t9 = [t100 * t74 + t104 * t80 + t105 * t79 + t108 * (Icges(4,1) * t108 + t174) + t109 * (Icges(4,2) * t109 + t175) + t123 * (Icges(3,1) * t123 + t176) + t125 * (Icges(3,2) * t125 + t177) + t99 * t75 + Icges(2,3) + m(6) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(4) * (t30 ^ 2 + t31 ^ 2) + m(3) * (t48 ^ 2 + t49 ^ 2) + m(2) * (t95 ^ 2 + t96 ^ 2); m(3) * (-t124 * t49 - t126 * t48) * t94 + m(6) * (t18 * t24 + t19 * t23) + m(5) * (t25 * t33 + t26 * t32) + m(4) * (t30 * t57 + t31 * t56) + t128 + (t108 * (Icges(4,5) * t124 + t140 * t126) + t109 * (Icges(4,6) * t124 + t137 * t126) + t123 * (Icges(3,5) * t124 + t141 * t126) + t125 * (Icges(3,6) * t124 + t138 * t126)) * t191 + (t108 * (-Icges(4,5) * t126 + t140 * t124) + t109 * (-Icges(4,6) * t126 + t137 * t124) + t123 * (-Icges(3,5) * t126 + t141 * t124) + t125 * (-Icges(3,6) * t126 + t138 * t124)) * t190 + (Icges(3,5) * t123 + Icges(4,5) * t108 + Icges(3,6) * t125 + Icges(4,6) * t109) * (t121 / 0.2e1 + t120 / 0.2e1); m(6) * (t1 ^ 2 + t23 ^ 2 + t24 ^ 2) + t4 + m(5) * (t32 ^ 2 + t33 ^ 2 + t6 ^ 2) + m(4) * (t56 ^ 2 + t57 ^ 2 + t8 ^ 2) + m(3) * (t167 * t94 ^ 2 + t27 ^ 2) + t166 + t192 * t124 * t120 + (-t5 + t193 * t121 + (t193 * t124 + t192 * t126) * t124) * t126; m(6) * (t124 * t18 - t126 * t19) + m(5) * (t124 * t25 - t126 * t26) + m(4) * (t124 * t30 - t126 * t31); m(6) * (t124 * t24 - t126 * t23) + m(5) * (t124 * t33 - t126 * t32) + m(4) * (t124 * t57 - t126 * t56); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t167; m(6) * (t18 * t35 + t19 * t34) + m(5) * (-t124 * t26 - t126 * t25) * t81 + t128; m(6) * (t1 * t7 + t23 * t34 + t24 * t35) + m(5) * (t20 * t6 + (-t124 * t32 - t126 * t33) * t81) + t145; m(6) * (t35 * t124 - t126 * t34); m(5) * (t167 * t81 ^ 2 + t20 ^ 2) + m(6) * (t34 ^ 2 + t35 ^ 2 + t7 ^ 2) + t145; m(6) * (-t124 * t19 - t126 * t18) * t76 + t163; m(6) * (t17 * t1 + (-t124 * t23 - t126 * t24) * t76) + t166; 0; m(6) * (t17 * t7 + (-t124 * t34 - t126 * t35) * t76) + t166; m(6) * (t167 * t76 ^ 2 + t17 ^ 2) + t166;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
