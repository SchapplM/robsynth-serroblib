% Calculate joint inertia matrix for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:21
% EndTime: 2019-12-05 15:40:26
% DurationCPUTime: 1.29s
% Computational Cost: add. (1507->208), mult. (3790->342), div. (0->0), fcn. (4095->6), ass. (0->103)
t100 = cos(qJ(2));
t97 = sin(qJ(4));
t98 = sin(qJ(2));
t99 = cos(qJ(4));
t67 = Icges(6,6) * t98 + (-Icges(6,5) * t97 + Icges(6,3) * t99) * t100;
t70 = Icges(5,6) * t98 + (-Icges(5,4) * t97 - Icges(5,2) * t99) * t100;
t152 = t67 - t70;
t71 = Icges(6,4) * t98 + (-Icges(6,1) * t97 + Icges(6,5) * t99) * t100;
t72 = Icges(5,5) * t98 + (-Icges(5,1) * t97 - Icges(5,4) * t99) * t100;
t151 = t71 + t72;
t150 = Icges(4,1) + Icges(3,3);
t149 = (Icges(4,5) - Icges(3,6)) * t98 + (-Icges(4,4) + Icges(3,5)) * t100;
t68 = Icges(5,3) * t98 + (-Icges(5,5) * t97 - Icges(5,6) * t99) * t100;
t69 = Icges(6,2) * t98 + (-Icges(6,4) * t97 + Icges(6,6) * t99) * t100;
t148 = (t68 + t69) * t98;
t147 = rSges(6,1) + pkin(4);
t95 = sin(pkin(7));
t91 = t95 ^ 2;
t96 = cos(pkin(7));
t92 = t96 ^ 2;
t128 = t91 + t92;
t146 = rSges(6,3) + qJ(5);
t126 = t100 * t96;
t120 = Icges(6,6) * t100;
t134 = t98 * t99;
t80 = -t96 * t134 + t95 * t97;
t137 = t97 * t98;
t81 = t96 * t137 + t95 * t99;
t33 = Icges(6,5) * t81 + Icges(6,3) * t80 + t96 * t120;
t122 = Icges(6,2) * t100;
t37 = Icges(6,4) * t81 + Icges(6,6) * t80 + t96 * t122;
t124 = Icges(6,4) * t100;
t41 = Icges(6,1) * t81 + Icges(6,5) * t80 + t96 * t124;
t11 = t37 * t126 + t80 * t33 + t81 * t41;
t82 = t95 * t134 + t96 * t97;
t84 = t95 * t137 - t96 * t99;
t34 = Icges(6,5) * t84 - Icges(6,3) * t82 + t95 * t120;
t38 = Icges(6,4) * t84 - Icges(6,6) * t82 + t95 * t122;
t42 = Icges(6,1) * t84 - Icges(6,5) * t82 + t95 * t124;
t12 = t38 * t126 + t80 * t34 + t81 * t42;
t119 = Icges(5,3) * t100;
t35 = Icges(5,5) * t81 - Icges(5,6) * t80 + t96 * t119;
t121 = Icges(5,6) * t100;
t39 = Icges(5,4) * t81 - Icges(5,2) * t80 + t96 * t121;
t123 = Icges(5,5) * t100;
t43 = Icges(5,1) * t81 - Icges(5,4) * t80 + t96 * t123;
t13 = t35 * t126 - t80 * t39 + t81 * t43;
t36 = Icges(5,5) * t84 + Icges(5,6) * t82 + t95 * t119;
t40 = Icges(5,4) * t84 + Icges(5,2) * t82 + t95 * t121;
t44 = Icges(5,1) * t84 + Icges(5,4) * t82 + t95 * t123;
t14 = t36 * t126 - t80 * t40 + t81 * t44;
t145 = (t151 * t81 + t152 * t80) * t98 + ((t11 + t13 + t148) * t96 + (t12 + t14) * t95) * t100;
t127 = t100 * t95;
t15 = t37 * t127 - t82 * t33 + t84 * t41;
t16 = t38 * t127 - t82 * t34 + t84 * t42;
t17 = t35 * t127 + t82 * t39 + t84 * t43;
t18 = t36 * t127 + t82 * t40 + t84 * t44;
t144 = (t151 * t84 - t152 * t82) * t98 + ((t15 + t17) * t96 + (t16 + t18 + t148) * t95) * t100;
t143 = -t149 * t95 + t150 * t96;
t142 = t149 * t96 + t150 * t95;
t141 = t98 ^ 2;
t133 = rSges(6,2) * t126 + t146 * t80 + t147 * t81;
t132 = rSges(6,2) * t127 - t146 * t82 + t147 * t84;
t131 = t98 * rSges(6,2) + (t146 * t99 - t147 * t97) * t100;
t130 = t128 * (pkin(2) * t100 + qJ(3) * t98);
t87 = t98 * pkin(2) - t100 * qJ(3);
t129 = t98 * rSges(4,2) + t100 * rSges(4,3) - t87;
t125 = t100 * t99;
t118 = -m(4) - m(5) - m(6);
t117 = -pkin(6) * t98 - t87;
t116 = t130 + (t96 * t126 + t95 * t127) * pkin(6);
t74 = t98 * rSges(5,3) + (-rSges(5,1) * t97 - rSges(5,2) * t99) * t100;
t115 = t117 - t74;
t107 = t117 - t131;
t94 = t100 ^ 2;
t89 = t98 * rSges(3,1) + t100 * rSges(3,2);
t54 = t129 * t96;
t53 = t129 * t95;
t50 = t84 * rSges(5,1) + t82 * rSges(5,2) + rSges(5,3) * t127;
t48 = t81 * rSges(5,1) - t80 * rSges(5,2) + rSges(5,3) * t126;
t46 = t115 * t96;
t45 = t115 * t95;
t32 = t128 * (rSges(3,1) * t100 - rSges(3,2) * t98);
t31 = t107 * t96;
t30 = t107 * t95;
t29 = -t74 * t126 + t98 * t48;
t28 = t74 * t127 - t98 * t50;
t27 = t130 + t128 * (-rSges(4,2) * t100 + rSges(4,3) * t98);
t26 = (-t48 * t95 + t50 * t96) * t100;
t25 = -t131 * t126 + t133 * t98;
t24 = t131 * t127 - t132 * t98;
t23 = t98 * t36 + (-t40 * t99 - t44 * t97) * t100;
t22 = t98 * t35 + (-t39 * t99 - t43 * t97) * t100;
t21 = t98 * t38 + (t34 * t99 - t42 * t97) * t100;
t20 = t98 * t37 + (t33 * t99 - t41 * t97) * t100;
t19 = t96 * t48 + t95 * t50 + t116;
t10 = (t132 * t96 - t133 * t95) * t100;
t9 = t132 * t95 + t133 * t96 + t116;
t8 = t17 * t95 - t18 * t96;
t7 = t15 * t95 - t16 * t96;
t6 = t13 * t95 - t14 * t96;
t5 = t11 * t95 - t12 * t96;
t1 = [m(2) + m(3) - t118; m(3) * t32 + m(4) * t27 + m(5) * t19 + m(6) * t9; m(6) * (t30 ^ 2 + t31 ^ 2 + t9 ^ 2) + m(5) * (t19 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(3) * (t128 * t89 ^ 2 + t32 ^ 2) + m(4) * (t27 ^ 2 + t53 ^ 2 + t54 ^ 2) + (t143 * t92 - t7 - t8) * t96 + (t5 + t6 + t142 * t91 + (t142 * t96 + t143 * t95) * t96) * t95; t118 * t100; m(6) * (-t100 * t9 + (t30 * t95 + t31 * t96) * t98) + m(5) * (-t100 * t19 + (t45 * t95 + t46 * t96) * t98) + m(4) * (-t100 * t27 + (t53 * t95 + t54 * t96) * t98); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t128 * t141 + t94); m(5) * t26 + m(6) * t10; m(6) * (t10 * t9 + t24 * t31 + t25 * t30) + m(5) * (t26 * t19 + t28 * t46 + t29 * t45) + ((t6 / 0.2e1 + t5 / 0.2e1) * t96 + (t7 / 0.2e1 + t8 / 0.2e1) * t95) * t100 + t145 * t95 / 0.2e1 - t144 * t96 / 0.2e1 + ((-t21 - t23) * t96 + (t20 + t22) * t95) * t98 / 0.2e1; m(5) * (-t26 * t100 + (t28 * t96 + t29 * t95) * t98) + m(6) * (-t10 * t100 + (t24 * t96 + t25 * t95) * t98); m(6) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t26 ^ 2 + t28 ^ 2 + t29 ^ 2) + t98 * (t141 * t68 + (t22 * t96 + t23 * t95 + (-t70 * t99 - t72 * t97) * t98) * t100) + t98 * (t141 * t69 + (t20 * t96 + t21 * t95 + (t67 * t99 - t71 * t97) * t98) * t100) + t144 * t127 + t145 * t126; m(6) * t125; m(6) * (t9 * t125 - t82 * t30 + t80 * t31); m(6) * (-t94 * t99 + (t80 * t96 - t82 * t95) * t98); m(6) * (t10 * t125 + t80 * t24 - t82 * t25); m(6) * (t94 * t99 ^ 2 + t80 ^ 2 + t82 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
