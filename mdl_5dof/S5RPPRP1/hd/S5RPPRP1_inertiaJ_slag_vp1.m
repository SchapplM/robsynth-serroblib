% Calculate joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:12
% EndTime: 2022-01-23 09:12:13
% DurationCPUTime: 0.79s
% Computational Cost: add. (1661->204), mult. (2045->288), div. (0->0), fcn. (2125->8), ass. (0->101)
t80 = sin(pkin(8));
t81 = cos(pkin(8));
t83 = sin(qJ(4));
t85 = cos(qJ(4));
t48 = -Icges(6,3) * t81 + (Icges(6,5) * t85 - Icges(6,6) * t83) * t80;
t49 = -Icges(5,3) * t81 + (Icges(5,5) * t85 - Icges(5,6) * t83) * t80;
t130 = -t48 - t49;
t50 = -Icges(6,6) * t81 + (Icges(6,4) * t85 - Icges(6,2) * t83) * t80;
t51 = -Icges(5,6) * t81 + (Icges(5,4) * t85 - Icges(5,2) * t83) * t80;
t129 = (-t50 - t51) * t83;
t52 = -Icges(6,5) * t81 + (Icges(6,1) * t85 - Icges(6,4) * t83) * t80;
t53 = -Icges(5,5) * t81 + (Icges(5,1) * t85 - Icges(5,4) * t83) * t80;
t128 = (t52 + t53) * t80 * t85;
t79 = qJ(1) + pkin(7);
t76 = cos(t79);
t116 = t76 * t83;
t113 = t81 * t83;
t75 = sin(t79);
t56 = -t75 * t113 - t76 * t85;
t112 = t81 * t85;
t57 = t75 * t112 - t116;
t127 = t57 * rSges(6,1) + t56 * rSges(6,2) - pkin(4) * t116;
t117 = t76 * t80;
t118 = t75 * t83;
t72 = t85 * pkin(4) + pkin(3);
t120 = t72 * t81;
t58 = -t76 * t113 + t75 * t85;
t59 = t76 * t112 + t118;
t126 = -t59 * rSges(6,1) - t58 * rSges(6,2) - rSges(6,3) * t117 - pkin(4) * t118 - t76 * t120;
t125 = t81 ^ 2;
t124 = pkin(3) * t81;
t84 = sin(qJ(1));
t123 = t84 * pkin(1);
t122 = -pkin(3) + t72;
t82 = -qJ(5) - pkin(6);
t121 = pkin(6) + t82;
t119 = t75 * t80;
t111 = t80 * t129 + t130 * t81 + t128;
t94 = Icges(6,3) * t80;
t21 = Icges(6,5) * t57 + Icges(6,6) * t56 + t75 * t94;
t95 = Icges(5,3) * t80;
t23 = Icges(5,5) * t57 + Icges(5,6) * t56 + t75 * t95;
t110 = t23 + t21;
t22 = Icges(6,5) * t59 + Icges(6,6) * t58 + t76 * t94;
t24 = Icges(5,5) * t59 + Icges(5,6) * t58 + t76 * t95;
t109 = t24 + t22;
t96 = Icges(6,6) * t80;
t26 = Icges(6,4) * t59 + Icges(6,2) * t58 + t76 * t96;
t97 = Icges(5,6) * t80;
t28 = Icges(5,4) * t59 + Icges(5,2) * t58 + t76 * t97;
t108 = t26 + t28;
t25 = Icges(6,4) * t57 + Icges(6,2) * t56 + t75 * t96;
t27 = Icges(5,4) * t57 + Icges(5,2) * t56 + t75 * t97;
t107 = t27 + t25;
t98 = Icges(6,5) * t80;
t30 = Icges(6,1) * t59 + Icges(6,4) * t58 + t76 * t98;
t99 = Icges(5,5) * t80;
t32 = Icges(5,1) * t59 + Icges(5,4) * t58 + t76 * t99;
t106 = t30 + t32;
t29 = Icges(6,1) * t57 + Icges(6,4) * t56 + t75 * t98;
t31 = Icges(5,1) * t57 + Icges(5,4) * t56 + t75 * t99;
t105 = t31 + t29;
t91 = t121 * t80;
t104 = (t122 * t81 - t91) * t75 + rSges(6,3) * t119 + t127;
t103 = -(-t91 - t124) * t76 + t126;
t102 = (t121 - rSges(6,3)) * t81 + (rSges(6,1) * t85 - rSges(6,2) * t83 + t122) * t80;
t100 = t75 ^ 2 + t76 ^ 2;
t38 = t59 * rSges(5,1) + t58 * rSges(5,2) + rSges(5,3) * t117;
t86 = cos(qJ(1));
t77 = t86 * pkin(1);
t92 = t76 * pkin(2) + t75 * qJ(3) + t77;
t90 = t76 * qJ(3) - t123;
t89 = rSges(4,1) * t81 - rSges(4,2) * t80;
t88 = -t57 * rSges(5,1) - t56 * rSges(5,2);
t67 = t86 * rSges(2,1) - t84 * rSges(2,2);
t66 = -t84 * rSges(2,1) - t86 * rSges(2,2);
t61 = t76 * rSges(3,1) - t75 * rSges(3,2) + t77;
t60 = -t75 * rSges(3,1) - t76 * rSges(3,2) - t123;
t55 = -t81 * rSges(5,3) + (rSges(5,1) * t85 - rSges(5,2) * t83) * t80;
t40 = t75 * rSges(4,3) + t89 * t76 + t92;
t39 = t76 * rSges(4,3) + (-pkin(2) - t89) * t75 + t90;
t36 = rSges(5,3) * t119 - t88;
t18 = -t55 * t117 - t81 * t38;
t17 = t55 * t119 + t81 * t36;
t16 = (pkin(6) * t80 + t124) * t76 + t92 + t38;
t15 = (-t124 - pkin(2) + (-rSges(5,3) - pkin(6)) * t80) * t75 + t88 + t90;
t14 = -t82 * t117 - t126 + t92;
t13 = (-t120 - pkin(2) + (-rSges(6,3) + t82) * t80) * t75 + t90 - t127;
t12 = t49 * t117 + t58 * t51 + t59 * t53;
t11 = t48 * t117 + t58 * t50 + t59 * t52;
t10 = t49 * t119 + t56 * t51 + t57 * t53;
t9 = t48 * t119 + t56 * t50 + t57 * t52;
t8 = (t36 * t76 - t38 * t75) * t80;
t7 = -t102 * t117 + t103 * t81;
t6 = t102 * t119 + t104 * t81;
t5 = -t81 * t24 + (-t28 * t83 + t32 * t85) * t80;
t4 = -t81 * t23 + (-t27 * t83 + t31 * t85) * t80;
t3 = -t81 * t22 + (-t26 * t83 + t30 * t85) * t80;
t2 = -t81 * t21 + (-t25 * t83 + t29 * t85) * t80;
t1 = (t103 * t75 + t104 * t76) * t80;
t19 = [Icges(2,3) + Icges(3,3) + (Icges(4,2) * t81 + t130) * t81 + (Icges(4,1) * t80 + 0.2e1 * Icges(4,4) * t81 + t129) * t80 + m(6) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(4) * (t39 ^ 2 + t40 ^ 2) + m(3) * (t60 ^ 2 + t61 ^ 2) + m(2) * (t66 ^ 2 + t67 ^ 2) + t128; 0; m(3) + m(4) + m(5) + m(6); m(6) * (t75 * t13 - t76 * t14) + m(5) * (t75 * t15 - t76 * t16) + m(4) * (t75 * t39 - t76 * t40); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t100; -t111 * t81 + m(6) * (t6 * t13 + t7 * t14) + m(5) * (t17 * t15 + t18 * t16) + ((t12 / 0.2e1 + t11 / 0.2e1 + t5 / 0.2e1 + t3 / 0.2e1) * t76 + (t4 / 0.2e1 + t2 / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1) * t75) * t80; m(5) * t8 + m(6) * t1; m(5) * (t17 * t75 - t18 * t76) + m(6) * (t6 * t75 - t7 * t76); m(5) * (t17 ^ 2 + t18 ^ 2 + t8 ^ 2) + m(6) * (t1 ^ 2 + t6 ^ 2 + t7 ^ 2) + t111 * t125 + (((t106 * t59 + t108 * t58 + t109 * t117) * t117 + (-t11 - t12 - t3 - t5) * t81) * t76 + ((t105 * t57 + t107 * t56 + t110 * t119) * t119 + (-t10 - t4 - t9 - t2) * t81 + ((t109 * t75 + t110 * t76) * t80 + t105 * t59 + t107 * t58 + t106 * t57 + t108 * t56) * t117) * t75) * t80; m(6) * (t13 * t76 + t14 * t75) * t80; -m(6) * t81; 0; m(6) * (-t81 * t1 + (t6 * t76 + t7 * t75) * t80); m(6) * (t100 * t80 ^ 2 + t125);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
