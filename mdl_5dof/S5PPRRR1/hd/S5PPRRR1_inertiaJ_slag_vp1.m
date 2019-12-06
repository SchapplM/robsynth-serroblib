% Calculate joint inertia matrix for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:37
% EndTime: 2019-12-05 15:12:40
% DurationCPUTime: 0.78s
% Computational Cost: add. (3366->150), mult. (3174->255), div. (0->0), fcn. (3384->8), ass. (0->84)
t87 = sin(pkin(8));
t84 = t87 ^ 2;
t88 = cos(pkin(8));
t85 = t88 ^ 2;
t129 = t84 + t85;
t86 = pkin(9) + qJ(3);
t81 = qJ(4) + t86;
t76 = sin(t81);
t77 = cos(t81);
t95 = Icges(5,4) * t77 - Icges(5,2) * t76;
t97 = Icges(5,1) * t77 - Icges(5,4) * t76;
t101 = -(Icges(5,6) * t87 + t95 * t88) * t76 + (Icges(5,5) * t87 + t97 * t88) * t77;
t102 = (-Icges(5,6) * t88 + t95 * t87) * t76 - (-Icges(5,5) * t88 + t97 * t87) * t77;
t93 = Icges(5,5) * t77 - Icges(5,6) * t76;
t51 = -Icges(5,3) * t88 + t93 * t87;
t52 = Icges(5,3) * t87 + t93 * t88;
t123 = t76 * t87;
t109 = Icges(6,3) * t76;
t91 = cos(qJ(5));
t117 = t88 * t91;
t90 = sin(qJ(5));
t120 = t87 * t90;
t65 = -t77 * t120 - t117;
t118 = t88 * t90;
t119 = t87 * t91;
t66 = t77 * t119 - t118;
t31 = Icges(6,5) * t66 + Icges(6,6) * t65 + t87 * t109;
t110 = Icges(6,6) * t76;
t33 = Icges(6,4) * t66 + Icges(6,2) * t65 + t87 * t110;
t111 = Icges(6,5) * t76;
t35 = Icges(6,1) * t66 + Icges(6,4) * t65 + t87 * t111;
t14 = t31 * t123 + t65 * t33 + t66 * t35;
t67 = -t77 * t118 + t119;
t68 = t77 * t117 + t120;
t32 = Icges(6,5) * t68 + Icges(6,6) * t67 + t88 * t109;
t34 = Icges(6,4) * t68 + Icges(6,2) * t67 + t88 * t110;
t36 = Icges(6,1) * t68 + Icges(6,4) * t67 + t88 * t111;
t15 = t32 * t123 + t65 * t34 + t66 * t36;
t8 = -t14 * t88 + t15 * t87;
t127 = -t85 * t51 - (t101 * t87 + (t102 - t52) * t88) * t87 - t8;
t79 = sin(t86);
t80 = cos(t86);
t94 = Icges(4,5) * t80 - Icges(4,6) * t79;
t126 = t87 * (Icges(4,3) * t87 + t94 * t88);
t125 = pkin(3) * t79;
t122 = t76 * t88;
t16 = t31 * t122 + t67 * t33 + t68 * t35;
t17 = t32 * t122 + t67 * t34 + t68 * t36;
t9 = -t16 * t88 + t17 * t87;
t124 = (t84 * t52 + t9 + (t102 * t88 + (t101 - t51) * t87) * t88) * t87;
t47 = -Icges(6,3) * t77 + (Icges(6,5) * t91 - Icges(6,6) * t90) * t76;
t121 = t77 * t47;
t116 = t129 * pkin(3) * t80;
t27 = t129 * (rSges(5,1) * t77 - rSges(5,2) * t76);
t50 = -t77 * rSges(6,3) + (rSges(6,1) * t91 - rSges(6,2) * t90) * t76;
t115 = -t76 * pkin(4) + t77 * pkin(7) - t50;
t70 = t76 * rSges(5,1) + t77 * rSges(5,2);
t108 = -t70 - t125;
t37 = t66 * rSges(6,1) + t65 * rSges(6,2) + rSges(6,3) * t123;
t38 = t68 * rSges(6,1) + t67 * rSges(6,2) + rSges(6,3) * t122;
t20 = t87 * t37 + t88 * t38 + t129 * (pkin(4) * t77 + pkin(7) * t76);
t18 = -t77 * t31 + (-t33 * t90 + t35 * t91) * t76;
t19 = -t77 * t32 + (-t34 * t90 + t36 * t91) * t76;
t48 = -Icges(6,6) * t77 + (Icges(6,4) * t91 - Icges(6,2) * t90) * t76;
t49 = -Icges(6,5) * t77 + (Icges(6,1) * t91 - Icges(6,4) * t90) * t76;
t3 = -(t65 * t48 + t66 * t49) * t77 + (t15 * t88 + (t14 - t121) * t87) * t76;
t4 = -(t67 * t48 + t68 * t49) * t77 + (t16 * t87 + (t17 - t121) * t88) * t76;
t107 = t87 * t4 / 0.2e1 - t77 * (-t18 * t88 + t19 * t87) / 0.2e1 - t88 * t3 / 0.2e1 + t8 * t123 / 0.2e1 + t9 * t122 / 0.2e1;
t106 = t115 - t125;
t92 = t127 * t88 + t124;
t73 = t79 * rSges(4,1) + t80 * rSges(4,2);
t46 = t108 * t88;
t45 = t108 * t87;
t40 = t115 * t88;
t39 = t115 * t87;
t30 = t129 * (rSges(4,1) * t80 - rSges(4,2) * t79);
t26 = t106 * t88;
t25 = t106 * t87;
t24 = -t50 * t122 - t77 * t38;
t23 = t50 * t123 + t77 * t37;
t22 = (t37 * t88 - t38 * t87) * t76;
t21 = t27 + t116;
t11 = t20 + t116;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t129; m(4) * t30 + m(5) * t21 + m(6) * t11; m(5) * (-t45 * t88 + t46 * t87) + m(6) * (-t25 * t88 + t26 * t87); m(6) * (t11 ^ 2 + t25 ^ 2 + t26 ^ 2) + t84 * t126 + m(4) * (t129 * t73 ^ 2 + t30 ^ 2) + m(5) * (t21 ^ 2 + t45 ^ 2 + t46 ^ 2) + t124 + (t88 * t126 - t129 * (-Icges(4,3) * t88 + t94 * t87) + t127) * t88; m(5) * t27 + m(6) * t20; m(6) * (-t39 * t88 + t40 * t87); m(6) * (t20 * t11 + t39 * t25 + t40 * t26) + m(5) * (t27 * t21 + (-t45 * t87 - t46 * t88) * t70) + t92; m(5) * (t129 * t70 ^ 2 + t27 ^ 2) + m(6) * (t20 ^ 2 + t39 ^ 2 + t40 ^ 2) + t92; m(6) * t22; m(6) * (t23 * t87 - t24 * t88); m(6) * (t22 * t11 + t23 * t26 + t24 * t25) + t107; m(6) * (t22 * t20 + t23 * t40 + t24 * t39) + t107; m(6) * (t22 ^ 2 + t23 ^ 2 + t24 ^ 2) + t4 * t122 + t3 * t123 - t77 * (t77 ^ 2 * t47 + (t19 * t88 + t18 * t87 - (-t48 * t90 + t49 * t91) * t77) * t76);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
