% Calculate joint inertia matrix for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:39
% EndTime: 2019-12-05 15:06:42
% DurationCPUTime: 1.09s
% Computational Cost: add. (2623->192), mult. (3361->310), div. (0->0), fcn. (3620->6), ass. (0->96)
t83 = pkin(8) + qJ(3);
t79 = sin(t83);
t80 = cos(t83);
t87 = sin(qJ(4));
t88 = cos(qJ(4));
t55 = -Icges(6,6) * t80 + (Icges(6,4) * t88 - Icges(6,2) * t87) * t79;
t56 = -Icges(5,6) * t80 + (Icges(5,4) * t88 - Icges(5,2) * t87) * t79;
t128 = -t55 - t56;
t57 = -Icges(6,5) * t80 + (Icges(6,1) * t88 - Icges(6,4) * t87) * t79;
t58 = -Icges(5,5) * t80 + (Icges(5,1) * t88 - Icges(5,4) * t87) * t79;
t127 = -t57 - t58;
t53 = -Icges(6,3) * t80 + (Icges(6,5) * t88 - Icges(6,6) * t87) * t79;
t54 = -Icges(5,3) * t80 + (Icges(5,5) * t88 - Icges(5,6) * t87) * t79;
t126 = (-t53 - t54) * t80;
t84 = sin(pkin(7));
t81 = t84 ^ 2;
t85 = cos(pkin(7));
t82 = t85 ^ 2;
t104 = t81 + t82;
t117 = t79 * t84;
t110 = t85 * t88;
t113 = t84 * t87;
t69 = -t80 * t113 - t110;
t111 = t85 * t87;
t112 = t84 * t88;
t70 = t80 * t112 - t111;
t98 = Icges(6,3) * t79;
t34 = Icges(6,5) * t70 + Icges(6,6) * t69 + t84 * t98;
t100 = Icges(6,6) * t79;
t38 = Icges(6,4) * t70 + Icges(6,2) * t69 + t84 * t100;
t102 = Icges(6,5) * t79;
t42 = Icges(6,1) * t70 + Icges(6,4) * t69 + t84 * t102;
t11 = t34 * t117 + t69 * t38 + t70 * t42;
t71 = -t80 * t111 + t112;
t72 = t80 * t110 + t113;
t35 = Icges(6,5) * t72 + Icges(6,6) * t71 + t85 * t98;
t39 = Icges(6,4) * t72 + Icges(6,2) * t71 + t85 * t100;
t43 = Icges(6,1) * t72 + Icges(6,4) * t71 + t85 * t102;
t12 = t35 * t117 + t69 * t39 + t70 * t43;
t99 = Icges(5,3) * t79;
t36 = Icges(5,5) * t70 + Icges(5,6) * t69 + t84 * t99;
t101 = Icges(5,6) * t79;
t40 = Icges(5,4) * t70 + Icges(5,2) * t69 + t84 * t101;
t103 = Icges(5,5) * t79;
t44 = Icges(5,1) * t70 + Icges(5,4) * t69 + t84 * t103;
t13 = t36 * t117 + t69 * t40 + t70 * t44;
t37 = Icges(5,5) * t72 + Icges(5,6) * t71 + t85 * t99;
t41 = Icges(5,4) * t72 + Icges(5,2) * t71 + t85 * t101;
t45 = Icges(5,1) * t72 + Icges(5,4) * t71 + t85 * t103;
t14 = t37 * t117 + t69 * t41 + t70 * t45;
t125 = (t127 * t70 + t128 * t69) * t80 + ((t12 + t14) * t85 + (t11 + t13 + t126) * t84) * t79;
t116 = t79 * t85;
t15 = t34 * t116 + t71 * t38 + t72 * t42;
t16 = t35 * t116 + t71 * t39 + t72 * t43;
t17 = t36 * t116 + t71 * t40 + t72 * t44;
t18 = t37 * t116 + t71 * t41 + t72 * t45;
t124 = (t127 * t72 + t128 * t71) * t80 + ((t16 + t18 + t126) * t85 + (t15 + t17) * t84) * t79;
t123 = t80 ^ 2;
t119 = t88 * pkin(4);
t89 = qJ(5) * t79 + t119 * t80;
t109 = t70 * rSges(6,1) + t69 * rSges(6,2) + rSges(6,3) * t117 - pkin(4) * t111 + t89 * t84;
t108 = t72 * rSges(6,1) + t71 * rSges(6,2) + rSges(6,3) * t116 + pkin(4) * t113 + t89 * t85;
t107 = (-qJ(5) - rSges(6,3)) * t80 + (rSges(6,1) * t88 - rSges(6,2) * t87 + t119) * t79;
t60 = -t80 * rSges(5,3) + (rSges(5,1) * t88 - rSges(5,2) * t87) * t79;
t75 = t79 * pkin(3) - t80 * pkin(6);
t106 = -t60 - t75;
t105 = t104 * (pkin(3) * t80 + pkin(6) * t79);
t97 = -t75 - t107;
t90 = Icges(4,5) * t80 - Icges(4,6) * t79;
t74 = t79 * rSges(4,1) + t80 * rSges(4,2);
t62 = Icges(4,3) * t84 + t90 * t85;
t61 = -Icges(4,3) * t85 + t90 * t84;
t51 = t106 * t85;
t50 = t106 * t84;
t49 = t72 * rSges(5,1) + t71 * rSges(5,2) + rSges(5,3) * t116;
t47 = t70 * rSges(5,1) + t69 * rSges(5,2) + rSges(5,3) * t117;
t31 = t104 * (rSges(4,1) * t80 - rSges(4,2) * t79);
t30 = t97 * t85;
t29 = t97 * t84;
t28 = -t60 * t116 - t80 * t49;
t27 = t60 * t117 + t80 * t47;
t26 = (t47 * t85 - t49 * t84) * t79;
t25 = t84 * t47 + t85 * t49 + t105;
t24 = -t80 * t37 + (-t41 * t87 + t45 * t88) * t79;
t23 = -t80 * t36 + (-t40 * t87 + t44 * t88) * t79;
t22 = -t80 * t35 + (-t39 * t87 + t43 * t88) * t79;
t21 = -t80 * t34 + (-t38 * t87 + t42 * t88) * t79;
t20 = -t107 * t116 - t108 * t80;
t19 = t107 * t117 + t109 * t80;
t10 = (-t108 * t84 + t109 * t85) * t79;
t9 = t108 * t85 + t109 * t84 + t105;
t8 = -t17 * t85 + t18 * t84;
t7 = -t15 * t85 + t16 * t84;
t6 = -t13 * t85 + t14 * t84;
t5 = -t11 * t85 + t12 * t84;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t104; m(4) * t31 + m(5) * t25 + m(6) * t9; m(5) * (-t50 * t85 + t51 * t84) + m(6) * (-t29 * t85 + t30 * t84); m(6) * (t29 ^ 2 + t30 ^ 2 + t9 ^ 2) + m(5) * (t25 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(4) * (t104 * t74 ^ 2 + t31 ^ 2) + (t81 * t62 + t7 + t8) * t84 + (-t82 * t61 - t5 - t6 + (-t84 * t61 + t85 * t62) * t84) * t85; m(5) * t26 + m(6) * t10; m(5) * (t27 * t84 - t28 * t85) + m(6) * (t19 * t84 - t20 * t85); m(6) * (t10 * t9 + t19 * t30 + t20 * t29) + m(5) * (t26 * t25 + t27 * t51 + t28 * t50) + ((t8 / 0.2e1 + t7 / 0.2e1) * t85 + (t6 / 0.2e1 + t5 / 0.2e1) * t84) * t79 - ((-t21 - t23) * t85 + (t22 + t24) * t84) * t80 / 0.2e1 + t124 * t84 / 0.2e1 - t125 * t85 / 0.2e1; m(6) * (t10 ^ 2 + t19 ^ 2 + t20 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) - t80 * (t123 * t53 + (t22 * t85 + t21 * t84 - (-t55 * t87 + t57 * t88) * t80) * t79) - t80 * (t123 * t54 + (t24 * t85 + t23 * t84 - (-t56 * t87 + t58 * t88) * t80) * t79) + t125 * t117 + t124 * t116; -m(6) * t80; 0; m(6) * (-t80 * t9 + (t29 * t84 + t30 * t85) * t79); m(6) * (-t80 * t10 + (t19 * t85 + t20 * t84) * t79); m(6) * (t104 * t79 ^ 2 + t123);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
