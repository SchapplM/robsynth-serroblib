% Calculate joint inertia matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:15
% EndTime: 2022-01-23 09:16:17
% DurationCPUTime: 1.17s
% Computational Cost: add. (3443->255), mult. (3575->379), div. (0->0), fcn. (3767->10), ass. (0->114)
t112 = sin(qJ(1));
t110 = cos(pkin(8));
t113 = cos(qJ(1));
t131 = t110 * t113;
t98 = t112 * qJ(2);
t133 = t113 * pkin(1) + t98;
t108 = sin(pkin(8));
t129 = t113 * t108;
t104 = pkin(9) + qJ(4);
t97 = qJ(5) + t104;
t92 = sin(t97);
t93 = cos(t97);
t71 = t112 * t93 - t92 * t131;
t72 = t112 * t92 + t93 * t131;
t44 = t72 * rSges(6,1) + t71 * rSges(6,2) + rSges(6,3) * t129;
t109 = cos(pkin(9));
t94 = t109 * pkin(3) + pkin(2);
t96 = cos(t104);
t83 = pkin(4) * t96 + t94;
t107 = sin(pkin(9));
t139 = t107 * pkin(3);
t95 = sin(t104);
t84 = pkin(4) * t95 + t139;
t144 = t112 * t84 + t83 * t131 + t133 + t44;
t143 = pkin(1) - (-rSges(4,3) - qJ(3)) * t108 + (rSges(4,1) * t109 - rSges(4,2) * t107 + pkin(2)) * t110;
t142 = m(4) / 0.2e1;
t141 = m(5) / 0.2e1;
t140 = m(6) / 0.2e1;
t60 = -Icges(6,6) * t110 + (Icges(6,4) * t93 - Icges(6,2) * t92) * t108;
t138 = t92 * t60;
t64 = -Icges(5,6) * t110 + (Icges(5,4) * t96 - Icges(5,2) * t95) * t108;
t137 = t95 * t64;
t111 = qJ(3) + pkin(6);
t102 = -pkin(7) - t111;
t80 = t111 * t108 + t94 * t110 + pkin(1);
t122 = -t102 * t108 - t80;
t91 = qJ(2) + t139;
t134 = t91 * t112;
t136 = -t122 * t113 + t134 - t144;
t130 = t112 * t108;
t132 = t110 * t112;
t69 = -t113 * t93 - t92 * t132;
t70 = -t113 * t92 + t93 * t132;
t117 = -t70 * rSges(6,1) - t69 * rSges(6,2);
t43 = rSges(6,3) * t130 - t117;
t62 = -t110 * rSges(6,3) + (rSges(6,1) * t93 - rSges(6,2) * t92) * t108;
t25 = t110 * t43 + t62 * t130;
t61 = -Icges(6,5) * t110 + (Icges(6,1) * t93 - Icges(6,4) * t92) * t108;
t54 = t108 * t93 * t61;
t59 = -Icges(6,3) * t110 + (Icges(6,5) * t93 - Icges(6,6) * t92) * t108;
t23 = -t108 * t138 - t110 * t59 + t54;
t135 = t23 * t110;
t128 = t112 ^ 2 + t113 ^ 2;
t77 = t112 * t96 - t95 * t131;
t78 = t112 * t95 + t96 * t131;
t52 = t78 * rSges(5,1) + t77 * rSges(5,2) + rSges(5,3) * t129;
t38 = Icges(6,5) * t72 + Icges(6,6) * t71 + Icges(6,3) * t129;
t40 = Icges(6,4) * t72 + Icges(6,2) * t71 + Icges(6,6) * t129;
t42 = Icges(6,1) * t72 + Icges(6,4) * t71 + Icges(6,5) * t129;
t10 = -t110 * t38 + (-t40 * t92 + t42 * t93) * t108;
t15 = t59 * t130 + t69 * t60 + t70 * t61;
t16 = t59 * t129 + t71 * t60 + t72 * t61;
t37 = Icges(6,5) * t70 + Icges(6,6) * t69 + Icges(6,3) * t130;
t39 = Icges(6,4) * t70 + Icges(6,2) * t69 + Icges(6,6) * t130;
t41 = Icges(6,1) * t70 + Icges(6,4) * t69 + Icges(6,5) * t130;
t9 = -t110 * t37 + (-t39 * t92 + t41 * t93) * t108;
t125 = (t15 + t9) * t130 / 0.2e1 + (t10 + t16) * t129 / 0.2e1;
t124 = -t110 * t83 - pkin(1);
t120 = -t110 * (-t135 + (t10 * t113 + t112 * t9) * t108) + ((t38 * t130 + t69 * t40 + t70 * t42) * t129 + (t37 * t130 + t69 * t39 + t70 * t41) * t130 - t15 * t110) * t130 + ((t38 * t129 + t71 * t40 + t72 * t42) * t129 + (t37 * t129 + t71 * t39 + t72 * t41) * t130 - t16 * t110) * t129;
t119 = t142 + t141 + t140;
t75 = -t113 * t96 - t95 * t132;
t76 = -t113 * t95 + t96 * t132;
t118 = -t76 * rSges(5,1) - t75 * rSges(5,2);
t116 = rSges(3,1) * t110 - rSges(3,2) * t108;
t114 = t107 * rSges(4,1) + t109 * rSges(4,2);
t100 = t113 * qJ(2);
t87 = t113 * rSges(2,1) - t112 * rSges(2,2);
t86 = -t112 * rSges(2,1) - t113 * rSges(2,2);
t66 = -t110 * rSges(5,3) + (rSges(5,1) * t96 - rSges(5,2) * t95) * t108;
t65 = -Icges(5,5) * t110 + (Icges(5,1) * t96 - Icges(5,4) * t95) * t108;
t63 = -Icges(5,3) * t110 + (Icges(5,5) * t96 - Icges(5,6) * t95) * t108;
t58 = t112 * rSges(3,3) + t116 * t113 + t133;
t57 = t113 * rSges(3,3) + t100 + (-pkin(1) - t116) * t112;
t55 = t108 * t96 * t65;
t53 = (t102 + t111) * t110 + (t83 - t94) * t108;
t51 = rSges(5,3) * t130 - t118;
t50 = Icges(5,1) * t78 + Icges(5,4) * t77 + Icges(5,5) * t129;
t49 = Icges(5,1) * t76 + Icges(5,4) * t75 + Icges(5,5) * t130;
t48 = Icges(5,4) * t78 + Icges(5,2) * t77 + Icges(5,6) * t129;
t47 = Icges(5,4) * t76 + Icges(5,2) * t75 + Icges(5,6) * t130;
t46 = Icges(5,5) * t78 + Icges(5,6) * t77 + Icges(5,3) * t129;
t45 = Icges(5,5) * t76 + Icges(5,6) * t75 + Icges(5,3) * t130;
t35 = t43 * t129;
t34 = -t143 * t112 + t114 * t113 + t100;
t33 = t114 * t112 + t143 * t113 + t98;
t32 = t80 * t113 + t134 + t52;
t31 = t91 * t113 + (-rSges(5,3) * t108 - t80) * t112 + t118;
t29 = -t100 + (-t84 + t91) * t113 + (t122 - t124) * t112;
t28 = -t110 * t52 - t66 * t129;
t27 = t110 * t51 + t66 * t130;
t26 = -t110 * t44 - t62 * t129;
t24 = -t108 * t137 - t110 * t63 + t55;
t22 = -t102 * t129 + t144;
t21 = t113 * t84 + t100 + ((-rSges(6,3) + t102) * t108 + t124) * t112 + t117;
t20 = (-t112 * t52 + t113 * t51) * t108;
t19 = t63 * t129 + t77 * t64 + t78 * t65;
t18 = t63 * t130 + t75 * t64 + t76 * t65;
t17 = -t44 * t130 + t35;
t12 = -t110 * t46 + (-t48 * t95 + t50 * t96) * t108;
t11 = -t110 * t45 + (-t47 * t95 + t49 * t96) * t108;
t6 = t136 * t110 + (-t53 - t62) * t129;
t5 = t110 * t29 + t53 * t130 + t25;
t4 = t35 + (t136 * t112 + t113 * t29) * t108;
t1 = [Icges(2,3) + t54 + t55 + (-t59 - t63 + (Icges(3,2) + Icges(4,3)) * t110) * t110 + m(6) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t31 ^ 2 + t32 ^ 2) + m(3) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(2) * (t86 ^ 2 + t87 ^ 2) + (-t138 - t137 + (Icges(4,1) * t109 ^ 2 + Icges(3,1) + (-0.2e1 * Icges(4,4) * t109 + Icges(4,2) * t107) * t107) * t108 + 0.2e1 * (-Icges(4,5) * t109 + Icges(4,6) * t107 + Icges(3,4)) * t110) * t108; m(6) * (t112 * t21 - t113 * t22) + m(5) * (t112 * t31 - t113 * t32) + m(3) * (t112 * t57 - t113 * t58) + m(4) * (t112 * t34 - t113 * t33); 0.2e1 * (m(3) / 0.2e1 + t119) * t128; 0.2e1 * ((t112 * t22 + t113 * t21) * t140 + (t112 * t32 + t113 * t31) * t141 + (t112 * t33 + t113 * t34) * t142) * t108; 0; 0.2e1 * t119 * (t128 * t108 ^ 2 + t110 ^ 2); (-t23 - t24) * t110 + m(6) * (t5 * t21 + t6 * t22) + m(5) * (t27 * t31 + t28 * t32) + ((t12 / 0.2e1 + t19 / 0.2e1) * t113 + (t11 / 0.2e1 + t18 / 0.2e1) * t112) * t108 + t125; m(5) * (t27 * t112 - t28 * t113) + m(6) * (t5 * t112 - t6 * t113); m(5) * (-t20 * t110 + (t112 * t28 + t113 * t27) * t108) + m(6) * (-t4 * t110 + (t112 * t6 + t113 * t5) * t108); ((t46 * t129 + t77 * t48 + t78 * t50) * t129 + (t45 * t129 + t77 * t47 + t78 * t49) * t130 - t19 * t110) * t129 + ((t46 * t130 + t75 * t48 + t76 * t50) * t129 + (t45 * t130 + t75 * t47 + t76 * t49) * t130 - t18 * t110) * t130 - t110 * (-t24 * t110 + (t11 * t112 + t113 * t12) * t108) + m(5) * (t20 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) + t120; m(6) * (t25 * t21 + t26 * t22) - t135 + t125; m(6) * (t25 * t112 - t26 * t113); m(6) * (-t17 * t110 + (t112 * t26 + t113 * t25) * t108); m(6) * (t17 * t4 + t25 * t5 + t26 * t6) + t120; m(6) * (t17 ^ 2 + t25 ^ 2 + t26 ^ 2) + t120;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
