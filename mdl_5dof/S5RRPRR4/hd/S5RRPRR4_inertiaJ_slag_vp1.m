% Calculate joint inertia matrix for
% S5RRPRR4
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:49
% EndTime: 2020-01-03 12:01:50
% DurationCPUTime: 0.61s
% Computational Cost: add. (2310->155), mult. (1388->219), div. (0->0), fcn. (1200->10), ass. (0->91)
t95 = qJ(1) + qJ(2);
t87 = pkin(9) + t95;
t82 = sin(t87);
t83 = cos(t87);
t137 = t82 * t83;
t94 = qJ(4) + qJ(5);
t88 = sin(t94);
t90 = cos(t94);
t139 = -rSges(6,1) * t90 + rSges(6,2) * t88;
t96 = sin(qJ(4));
t98 = cos(qJ(4));
t138 = -rSges(5,1) * t98 + rSges(5,2) * t96;
t29 = -t82 * rSges(6,3) + t139 * t83;
t86 = t98 * pkin(4) + pkin(3);
t136 = t83 * t86 - t29;
t77 = t82 ^ 2;
t78 = t83 ^ 2;
t135 = -t82 / 0.2e1;
t134 = -t83 / 0.2e1;
t65 = t96 * rSges(5,1) + t98 * rSges(5,2);
t133 = m(5) * t65;
t50 = t88 * rSges(6,1) + t90 * rSges(6,2);
t132 = m(6) * t50;
t100 = -pkin(8) - pkin(7);
t127 = t83 * t100 + t82 * t86;
t126 = -t83 * pkin(3) - t82 * pkin(7);
t125 = t77 + t78;
t89 = sin(t95);
t91 = cos(t95);
t51 = t89 * rSges(3,1) + t91 * rSges(3,2);
t124 = Icges(5,4) * t96;
t123 = Icges(5,4) * t98;
t122 = Icges(6,4) * t88;
t121 = Icges(6,4) * t90;
t84 = pkin(2) * t89;
t40 = t82 * rSges(4,1) + t83 * rSges(4,2) + t84;
t107 = -Icges(6,2) * t88 + t121;
t109 = Icges(6,1) * t90 - t122;
t48 = Icges(6,2) * t90 + t122;
t49 = Icges(6,1) * t88 + t121;
t112 = t48 * t88 - t49 * t90;
t47 = Icges(6,5) * t88 + Icges(6,6) * t90;
t120 = (t112 * t83 + t90 * (-Icges(6,6) * t82 - t107 * t83) + t88 * (-Icges(6,5) * t82 - t109 * t83) - t82 * t47) * t135 + (-t112 * t82 + t90 * (-Icges(6,6) * t83 + t107 * t82) + t88 * (-Icges(6,5) * t83 + t109 * t82) - t83 * t47) * t134;
t119 = pkin(4) * t96 + t50;
t52 = t91 * rSges(3,1) - t89 * rSges(3,2);
t118 = t138 * t82;
t85 = pkin(2) * t91;
t41 = t83 * rSges(4,1) - t82 * rSges(4,2) + t85;
t105 = Icges(6,5) * t90 - Icges(6,6) * t88;
t23 = -Icges(6,3) * t83 + t105 * t82;
t24 = -Icges(6,3) * t82 - t105 * t83;
t117 = -t83 * (t24 * t137 + t78 * t23) - t82 * (t23 * t137 + t77 * t24);
t63 = Icges(5,2) * t98 + t124;
t64 = Icges(5,1) * t96 + t123;
t111 = t63 * t96 - t64 * t98;
t110 = Icges(5,1) * t98 - t124;
t108 = -Icges(5,2) * t96 + t123;
t106 = Icges(5,5) * t98 - Icges(5,6) * t96;
t104 = -t82 * rSges(5,3) + t138 * t83;
t103 = t90 * t48 + t88 * t49 + t98 * t63 + t96 * t64 + Icges(3,3) + Icges(4,3);
t62 = Icges(5,5) * t96 + Icges(5,6) * t98;
t102 = t120 + (t111 * t83 + t98 * (-Icges(5,6) * t82 - t108 * t83) + t96 * (-Icges(5,5) * t82 - t110 * t83) - t82 * t62) * t135 + (-t111 * t82 + t98 * (-Icges(5,6) * t83 + t108 * t82) + t96 * (-Icges(5,5) * t83 + t110 * t82) - t83 * t62) * t134;
t101 = -t83 * rSges(6,3) - t139 * t82;
t21 = -t104 + t85 - t126;
t16 = t101 + t84 + t127;
t17 = -t82 * t100 + t136 + t85;
t75 = t82 * pkin(3);
t20 = t75 + t84 + (-rSges(5,3) - pkin(7)) * t83 - t118;
t99 = cos(qJ(1));
t97 = sin(qJ(1));
t93 = t99 * pkin(1);
t92 = t97 * pkin(1);
t67 = t99 * rSges(2,1) - t97 * rSges(2,2);
t66 = t97 * rSges(2,1) + t99 * rSges(2,2);
t45 = t52 + t93;
t44 = t92 + t51;
t39 = t41 + t93;
t38 = t92 + t40;
t33 = -Icges(5,3) * t82 - t106 * t83;
t32 = -Icges(5,3) * t83 + t106 * t82;
t31 = t119 * t83;
t30 = t119 * t82;
t22 = t82 * t101;
t19 = t21 + t93;
t18 = t92 + t20;
t13 = t17 + t93;
t12 = t16 + t92;
t11 = -t83 * t104 + t82 * (-t83 * rSges(5,3) - t118);
t6 = -t83 * t29 + t22;
t3 = t82 * (-t75 + t127) + t22 + ((pkin(7) - t100) * t82 + t126 + t136) * t83;
t1 = [Icges(2,3) + m(2) * (t66 ^ 2 + t67 ^ 2) + m(3) * (t44 ^ 2 + t45 ^ 2) + m(4) * (t38 ^ 2 + t39 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(6) * (t12 ^ 2 + t13 ^ 2) + t103; m(3) * (t51 * t44 + t52 * t45) + m(4) * (t40 * t38 + t41 * t39) + m(5) * (t20 * t18 + t21 * t19) + m(6) * (t16 * t12 + t17 * t13) + t103; m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t40 ^ 2 + t41 ^ 2) + m(3) * (t51 ^ 2 + t52 ^ 2) + t103; 0; 0; m(4) + m(5) + m(6); m(6) * (t31 * t12 - t30 * t13) + (t18 * t83 - t19 * t82) * t133 + t102; m(6) * (t31 * t16 - t30 * t17) + (t20 * t83 - t21 * t82) * t133 + t102; m(5) * t11 + m(6) * t3; m(5) * (t125 * t65 ^ 2 + t11 ^ 2) - t83 * (t33 * t137 + t78 * t32) - t82 * (t32 * t137 + t77 * t33) + m(6) * (t3 ^ 2 + t30 ^ 2 + t31 ^ 2) + t117; (t12 * t83 - t13 * t82) * t132 + t120; (t16 * t83 - t17 * t82) * t132 + t120; m(6) * t6; m(6) * (t6 * t3 + (t30 * t82 + t31 * t83) * t50) + t117; m(6) * (t125 * t50 ^ 2 + t6 ^ 2) + t117;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
