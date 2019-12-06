% Calculate joint inertia matrix for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:04
% EndTime: 2019-12-05 15:35:09
% DurationCPUTime: 1.31s
% Computational Cost: add. (2601->202), mult. (3584->326), div. (0->0), fcn. (3880->8), ass. (0->106)
t101 = cos(qJ(4));
t95 = qJ(2) + pkin(8);
t91 = sin(t95);
t92 = cos(t95);
t99 = sin(qJ(4));
t55 = -Icges(6,6) * t92 + (Icges(6,5) * t101 + Icges(6,3) * t99) * t91;
t58 = -Icges(5,6) * t92 + (Icges(5,4) * t101 - Icges(5,2) * t99) * t91;
t158 = -t55 + t58;
t59 = -Icges(6,4) * t92 + (Icges(6,1) * t101 + Icges(6,5) * t99) * t91;
t60 = -Icges(5,5) * t92 + (Icges(5,1) * t101 - Icges(5,4) * t99) * t91;
t157 = -t59 - t60;
t96 = sin(pkin(7));
t93 = t96 ^ 2;
t97 = cos(pkin(7));
t94 = t97 ^ 2;
t129 = t93 + t94;
t156 = Icges(3,3) + Icges(4,3);
t100 = sin(qJ(2));
t102 = cos(qJ(2));
t155 = Icges(3,5) * t102 + Icges(4,5) * t92 - Icges(3,6) * t100 - Icges(4,6) * t91;
t56 = -Icges(5,3) * t92 + (Icges(5,5) * t101 - Icges(5,6) * t99) * t91;
t57 = -Icges(6,2) * t92 + (Icges(6,4) * t101 + Icges(6,6) * t99) * t91;
t154 = (-t57 - t56) * t92;
t153 = rSges(6,1) + pkin(4);
t152 = rSges(6,3) + qJ(5);
t140 = t91 * t96;
t124 = Icges(6,6) * t91;
t121 = t97 * t101;
t135 = t96 * t99;
t80 = t92 * t135 + t121;
t122 = t96 * t101;
t134 = t97 * t99;
t81 = t92 * t122 - t134;
t34 = Icges(6,5) * t81 + Icges(6,3) * t80 + t96 * t124;
t126 = Icges(6,2) * t91;
t38 = Icges(6,4) * t81 + Icges(6,6) * t80 + t96 * t126;
t128 = Icges(6,4) * t91;
t42 = Icges(6,1) * t81 + Icges(6,5) * t80 + t96 * t128;
t12 = t38 * t140 + t80 * t34 + t81 * t42;
t82 = t92 * t134 - t122;
t83 = t92 * t121 + t135;
t35 = Icges(6,5) * t83 + Icges(6,3) * t82 + t97 * t124;
t39 = Icges(6,4) * t83 + Icges(6,6) * t82 + t97 * t126;
t43 = Icges(6,1) * t83 + Icges(6,5) * t82 + t97 * t128;
t13 = t39 * t140 + t80 * t35 + t81 * t43;
t123 = Icges(5,3) * t91;
t36 = Icges(5,5) * t81 - Icges(5,6) * t80 + t96 * t123;
t125 = Icges(5,6) * t91;
t40 = Icges(5,4) * t81 - Icges(5,2) * t80 + t96 * t125;
t127 = Icges(5,5) * t91;
t44 = Icges(5,1) * t81 - Icges(5,4) * t80 + t96 * t127;
t14 = t36 * t140 - t80 * t40 + t81 * t44;
t37 = Icges(5,5) * t83 - Icges(5,6) * t82 + t97 * t123;
t41 = Icges(5,4) * t83 - Icges(5,2) * t82 + t97 * t125;
t45 = Icges(5,1) * t83 - Icges(5,4) * t82 + t97 * t127;
t15 = t37 * t140 - t80 * t41 + t81 * t45;
t151 = (t157 * t81 + t158 * t80) * t92 + ((t13 + t15) * t97 + (t12 + t14 + t154) * t96) * t91;
t139 = t91 * t97;
t16 = t38 * t139 + t82 * t34 + t83 * t42;
t17 = t39 * t139 + t82 * t35 + t83 * t43;
t18 = t36 * t139 - t82 * t40 + t83 * t44;
t19 = t37 * t139 - t82 * t41 + t83 * t45;
t150 = (t157 * t83 + t158 * t82) * t92 + ((t17 + t19 + t154) * t97 + (t16 + t18) * t96) * t91;
t149 = -t155 * t96 + t156 * t97;
t148 = t155 * t97 + t156 * t96;
t147 = t92 ^ 2;
t141 = pkin(2) * t100;
t138 = t91 * t99;
t133 = rSges(6,2) * t140 + t152 * t80 + t153 * t81;
t132 = rSges(6,2) * t139 + t152 * t82 + t153 * t83;
t131 = t129 * t102 * pkin(2);
t130 = -t92 * rSges(6,2) + (t153 * t101 + t152 * t99) * t91;
t120 = -t91 * rSges(4,1) - t92 * rSges(4,2) - t141;
t119 = -t91 * pkin(3) + t92 * pkin(6) - t141;
t118 = t131 + t129 * (pkin(3) * t92 + pkin(6) * t91);
t62 = -t92 * rSges(5,3) + (rSges(5,1) * t101 - rSges(5,2) * t99) * t91;
t117 = t119 - t62;
t106 = t119 - t130;
t87 = t100 * rSges(3,1) + t102 * rSges(3,2);
t64 = t120 * t97;
t63 = t120 * t96;
t50 = t129 * (rSges(3,1) * t102 - rSges(3,2) * t100);
t49 = t83 * rSges(5,1) - t82 * rSges(5,2) + rSges(5,3) * t139;
t47 = t81 * rSges(5,1) - t80 * rSges(5,2) + rSges(5,3) * t140;
t33 = t117 * t97;
t32 = t117 * t96;
t31 = t106 * t97;
t30 = t106 * t96;
t29 = -t62 * t139 - t92 * t49;
t28 = t62 * t140 + t92 * t47;
t27 = t131 + t129 * (rSges(4,1) * t92 - rSges(4,2) * t91);
t26 = (t47 * t97 - t49 * t96) * t91;
t25 = -t130 * t139 - t132 * t92;
t24 = t130 * t140 + t133 * t92;
t23 = -t92 * t37 + (t101 * t45 - t41 * t99) * t91;
t22 = -t92 * t36 + (t101 * t44 - t40 * t99) * t91;
t21 = -t92 * t39 + (t101 * t43 + t35 * t99) * t91;
t20 = -t92 * t38 + (t101 * t42 + t34 * t99) * t91;
t11 = t96 * t47 + t97 * t49 + t118;
t10 = (-t132 * t96 + t133 * t97) * t91;
t9 = t132 * t97 + t133 * t96 + t118;
t8 = -t18 * t97 + t19 * t96;
t7 = -t16 * t97 + t17 * t96;
t6 = -t14 * t97 + t15 * t96;
t5 = -t12 * t97 + t13 * t96;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t50 + m(4) * t27 + m(5) * t11 + m(6) * t9; m(6) * (t30 ^ 2 + t31 ^ 2 + t9 ^ 2) + m(5) * (t11 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(4) * (t27 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(3) * (t129 * t87 ^ 2 + t50 ^ 2) + (t148 * t93 + t7 + t8) * t96 + (-t5 - t6 + t149 * t94 + (t148 * t97 + t149 * t96) * t96) * t97; 0; m(6) * (-t97 * t30 + t96 * t31) + m(5) * (-t97 * t32 + t96 * t33) + m(4) * (-t97 * t63 + t96 * t64); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t129; m(5) * t26 + m(6) * t10; m(6) * (t10 * t9 + t24 * t31 + t25 * t30) + m(5) * (t26 * t11 + t28 * t33 + t29 * t32) + ((t8 / 0.2e1 + t7 / 0.2e1) * t97 + (t6 / 0.2e1 + t5 / 0.2e1) * t96) * t91 - ((-t20 - t22) * t97 + (t21 + t23) * t96) * t92 / 0.2e1 + t150 * t96 / 0.2e1 - t151 * t97 / 0.2e1; m(5) * (t28 * t96 - t29 * t97) + m(6) * (t24 * t96 - t25 * t97); m(6) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) - t92 * (t147 * t57 + (t21 * t97 + t20 * t96 - (t101 * t59 + t55 * t99) * t92) * t91) + m(5) * (t26 ^ 2 + t28 ^ 2 + t29 ^ 2) - t92 * (t147 * t56 + (t23 * t97 + t22 * t96 - (t101 * t60 - t58 * t99) * t92) * t91) + t151 * t140 + t150 * t139; m(6) * t138; m(6) * (t9 * t138 + t80 * t30 + t82 * t31); m(6) * (-t80 * t97 + t82 * t96); m(6) * (t10 * t138 + t82 * t24 + t80 * t25); m(6) * (t91 ^ 2 * t99 ^ 2 + t80 ^ 2 + t82 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
