% Calculate joint inertia matrix for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:58
% EndTime: 2019-12-05 15:24:03
% DurationCPUTime: 1.03s
% Computational Cost: add. (2288->178), mult. (2416->290), div. (0->0), fcn. (2524->10), ass. (0->102)
t86 = sin(pkin(7));
t81 = t86 ^ 2;
t88 = cos(pkin(7));
t82 = t88 ^ 2;
t119 = t81 + t82;
t84 = qJ(2) + pkin(8);
t78 = sin(t84);
t80 = cos(t84);
t91 = sin(qJ(2));
t92 = cos(qJ(2));
t141 = Icges(3,5) * t92 + Icges(4,5) * t80 - Icges(3,6) * t91 - Icges(4,6) * t78;
t140 = Icges(3,3) + Icges(4,3);
t139 = t140 * t88 - t141 * t86;
t138 = t140 * t86 + t141 * t88;
t137 = t80 ^ 2;
t136 = t86 / 0.2e1;
t135 = -m(5) - m(6);
t134 = pkin(2) * t91;
t131 = t78 * t86;
t130 = t78 * t88;
t83 = pkin(9) + qJ(5);
t77 = sin(t83);
t79 = cos(t83);
t37 = -Icges(6,3) * t80 + (Icges(6,5) * t79 - Icges(6,6) * t77) * t78;
t129 = t80 * t37;
t128 = t86 * t77;
t127 = t86 * t79;
t85 = sin(pkin(9));
t126 = t86 * t85;
t87 = cos(pkin(9));
t125 = t86 * t87;
t124 = t88 * t77;
t123 = t88 * t79;
t122 = t88 * t85;
t121 = t88 * t87;
t120 = t119 * t92 * pkin(2);
t118 = Icges(5,5) * t78;
t117 = Icges(6,5) * t78;
t116 = Icges(5,6) * t78;
t115 = Icges(6,6) * t78;
t114 = Icges(5,3) * t78;
t113 = Icges(6,3) * t78;
t112 = m(5) / 0.2e1 + m(6) / 0.2e1;
t111 = -t78 * pkin(3) + t80 * qJ(4) - t134;
t110 = -t78 * rSges(4,1) - t80 * rSges(4,2) - t134;
t105 = pkin(3) * t80 + qJ(4) * t78;
t109 = t119 * t105 + t120;
t108 = t111 + t80 * rSges(5,3) - (rSges(5,1) * t87 - rSges(5,2) * t85) * t78;
t40 = -t80 * rSges(6,3) + (rSges(6,1) * t79 - rSges(6,2) * t77) * t78;
t74 = t87 * pkin(4) + pkin(3);
t90 = -pkin(6) - qJ(4);
t100 = t111 - (qJ(4) + t90) * t80 - (-pkin(3) + t74) * t78 - t40;
t93 = t74 * t80 - t78 * t90 - t105;
t72 = t91 * rSges(3,1) + t92 * rSges(3,2);
t68 = t80 * t121 + t126;
t67 = -t80 * t122 + t125;
t66 = t80 * t125 - t122;
t65 = -t80 * t126 - t121;
t58 = t80 * t123 + t128;
t57 = -t80 * t124 + t127;
t56 = t80 * t127 - t124;
t55 = -t80 * t128 - t123;
t45 = t110 * t88;
t44 = t110 * t86;
t39 = -Icges(6,5) * t80 + (Icges(6,1) * t79 - Icges(6,4) * t77) * t78;
t38 = -Icges(6,6) * t80 + (Icges(6,4) * t79 - Icges(6,2) * t77) * t78;
t35 = t119 * (rSges(3,1) * t92 - rSges(3,2) * t91);
t34 = Icges(5,1) * t68 + Icges(5,4) * t67 + t88 * t118;
t33 = Icges(5,1) * t66 + Icges(5,4) * t65 + t86 * t118;
t32 = Icges(5,4) * t68 + Icges(5,2) * t67 + t88 * t116;
t31 = Icges(5,4) * t66 + Icges(5,2) * t65 + t86 * t116;
t30 = Icges(5,5) * t68 + Icges(5,6) * t67 + t88 * t114;
t29 = Icges(5,5) * t66 + Icges(5,6) * t65 + t86 * t114;
t28 = t108 * t88;
t27 = t108 * t86;
t26 = t58 * rSges(6,1) + t57 * rSges(6,2) + rSges(6,3) * t130;
t25 = t56 * rSges(6,1) + t55 * rSges(6,2) + rSges(6,3) * t131;
t24 = Icges(6,1) * t58 + Icges(6,4) * t57 + t88 * t117;
t23 = Icges(6,1) * t56 + Icges(6,4) * t55 + t86 * t117;
t22 = Icges(6,4) * t58 + Icges(6,2) * t57 + t88 * t115;
t21 = Icges(6,4) * t56 + Icges(6,2) * t55 + t86 * t115;
t20 = Icges(6,5) * t58 + Icges(6,6) * t57 + t88 * t113;
t19 = Icges(6,5) * t56 + Icges(6,6) * t55 + t86 * t113;
t18 = t100 * t88;
t17 = t100 * t86;
t16 = -t40 * t130 - t80 * t26;
t15 = t40 * t131 + t80 * t25;
t14 = t120 + t119 * (rSges(4,1) * t80 - rSges(4,2) * t78);
t13 = (t25 * t88 - t26 * t86) * t78;
t12 = -t80 * t20 + (-t22 * t77 + t24 * t79) * t78;
t11 = -t80 * t19 + (-t21 * t77 + t23 * t79) * t78;
t10 = t86 * (t66 * rSges(5,1) + t65 * rSges(5,2) + rSges(5,3) * t131) + t88 * (t68 * rSges(5,1) + t67 * rSges(5,2) + rSges(5,3) * t130) + t109;
t9 = t20 * t130 + t57 * t22 + t58 * t24;
t8 = t19 * t130 + t57 * t21 + t58 * t23;
t7 = t20 * t131 + t55 * t22 + t56 * t24;
t6 = t19 * t131 + t55 * t21 + t56 * t23;
t5 = (t93 * t88 + t26) * t88 + (t93 * t86 + t25) * t86 + t109;
t4 = -t8 * t88 + t9 * t86;
t3 = -t6 * t88 + t7 * t86;
t2 = -(t57 * t38 + t58 * t39) * t80 + (t8 * t86 + (t9 - t129) * t88) * t78;
t1 = -(t55 * t38 + t56 * t39) * t80 + (t7 * t88 + (t6 - t129) * t86) * t78;
t36 = [m(2) + m(3) + m(4) - t135; m(3) * t35 + m(4) * t14 + m(5) * t10 + m(6) * t5; m(6) * (t17 ^ 2 + t18 ^ 2 + t5 ^ 2) + m(5) * (t10 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(4) * (t14 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(3) * (t119 * t72 ^ 2 + t35 ^ 2) + (-t3 + (t29 * t131 + t65 * t31 + t66 * t33) * t88 + t139 * t82) * t88 + (t4 + (t30 * t130 + t67 * t32 + t68 * t34) * t86 + t138 * t81 + (-t29 * t130 - t30 * t131 + t138 * t88 + t139 * t86 - t67 * t31 - t65 * t32 - t68 * t33 - t66 * t34) * t88) * t86; 0; m(6) * (-t88 * t17 + t86 * t18) + m(5) * (-t88 * t27 + t86 * t28) + m(4) * (-t88 * t44 + t86 * t45); 0.2e1 * (m(4) / 0.2e1 + t112) * t119; t135 * t80; m(6) * (-t80 * t5 + (t17 * t86 + t18 * t88) * t78) + m(5) * (-t80 * t10 + (t27 * t86 + t28 * t88) * t78); 0; 0.2e1 * t112 * (t119 * t78 ^ 2 + t137); m(6) * t13; t2 * t136 - t88 * t1 / 0.2e1 - t80 * (-t11 * t88 + t12 * t86) / 0.2e1 + m(6) * (t13 * t5 + t15 * t18 + t16 * t17) + (t88 * t4 / 0.2e1 + t3 * t136) * t78; m(6) * (t15 * t86 - t16 * t88); m(6) * (-t13 * t80 + (t15 * t88 + t16 * t86) * t78); m(6) * (t13 ^ 2 + t15 ^ 2 + t16 ^ 2) + t2 * t130 + t1 * t131 - t80 * (t137 * t37 + (t12 * t88 + t11 * t86 - (-t38 * t77 + t39 * t79) * t80) * t78);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t36(1), t36(2), t36(4), t36(7), t36(11); t36(2), t36(3), t36(5), t36(8), t36(12); t36(4), t36(5), t36(6), t36(9), t36(13); t36(7), t36(8), t36(9), t36(10), t36(14); t36(11), t36(12), t36(13), t36(14), t36(15);];
Mq = res;
