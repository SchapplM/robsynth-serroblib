% Calculate joint inertia matrix for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP6_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP6_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:13
% EndTime: 2019-12-31 16:30:15
% DurationCPUTime: 1.01s
% Computational Cost: add. (1230->176), mult. (3092->289), div. (0->0), fcn. (3426->6), ass. (0->94)
t79 = sin(qJ(3));
t80 = sin(qJ(2));
t81 = cos(qJ(3));
t82 = cos(qJ(2));
t58 = -Icges(5,6) * t82 + (Icges(5,5) * t81 + Icges(5,3) * t79) * t80;
t61 = -Icges(4,6) * t82 + (Icges(4,4) * t81 - Icges(4,2) * t79) * t80;
t122 = -t58 + t61;
t62 = -Icges(5,4) * t82 + (Icges(5,1) * t81 + Icges(5,5) * t79) * t80;
t63 = -Icges(4,5) * t82 + (Icges(4,1) * t81 - Icges(4,4) * t79) * t80;
t121 = -t62 - t63;
t59 = -Icges(4,3) * t82 + (Icges(4,5) * t81 - Icges(4,6) * t79) * t80;
t60 = -Icges(5,2) * t82 + (Icges(5,4) * t81 + Icges(5,6) * t79) * t80;
t120 = (-t59 - t60) * t82;
t119 = rSges(5,1) + pkin(3);
t118 = rSges(5,3) + qJ(4);
t78 = cos(pkin(6));
t113 = t78 ^ 2;
t77 = sin(pkin(6));
t114 = t77 ^ 2;
t115 = t113 + t114;
t108 = t77 * t80;
t105 = t79 * t82;
t68 = t77 * t105 + t78 * t81;
t104 = t81 * t82;
t69 = t77 * t104 - t78 * t79;
t92 = Icges(5,6) * t80;
t32 = Icges(5,5) * t69 + Icges(5,3) * t68 + t77 * t92;
t94 = Icges(5,2) * t80;
t36 = Icges(5,4) * t69 + Icges(5,6) * t68 + t77 * t94;
t96 = Icges(5,4) * t80;
t40 = Icges(5,1) * t69 + Icges(5,5) * t68 + t77 * t96;
t11 = t36 * t108 + t68 * t32 + t69 * t40;
t70 = t78 * t105 - t77 * t81;
t71 = t78 * t104 + t77 * t79;
t33 = Icges(5,5) * t71 + Icges(5,3) * t70 + t78 * t92;
t37 = Icges(5,4) * t71 + Icges(5,6) * t70 + t78 * t94;
t41 = Icges(5,1) * t71 + Icges(5,5) * t70 + t78 * t96;
t12 = t37 * t108 + t68 * t33 + t69 * t41;
t91 = Icges(4,3) * t80;
t34 = Icges(4,5) * t69 - Icges(4,6) * t68 + t77 * t91;
t93 = Icges(4,6) * t80;
t38 = Icges(4,4) * t69 - Icges(4,2) * t68 + t77 * t93;
t95 = Icges(4,5) * t80;
t42 = Icges(4,1) * t69 - Icges(4,4) * t68 + t77 * t95;
t13 = t34 * t108 - t68 * t38 + t69 * t42;
t35 = Icges(4,5) * t71 - Icges(4,6) * t70 + t78 * t91;
t39 = Icges(4,4) * t71 - Icges(4,2) * t70 + t78 * t93;
t43 = Icges(4,1) * t71 - Icges(4,4) * t70 + t78 * t95;
t14 = t35 * t108 - t68 * t39 + t69 * t43;
t117 = (t121 * t69 + t122 * t68) * t82 + ((t12 + t14) * t78 + (t11 + t13 + t120) * t77) * t80;
t107 = t78 * t80;
t15 = t36 * t107 + t70 * t32 + t71 * t40;
t16 = t37 * t107 + t70 * t33 + t71 * t41;
t17 = t34 * t107 - t70 * t38 + t71 * t42;
t18 = t35 * t107 - t70 * t39 + t71 * t43;
t116 = (t121 * t71 + t122 * t70) * t82 + ((t16 + t18 + t120) * t78 + (t15 + t17) * t77) * t80;
t112 = t82 ^ 2;
t106 = t79 * t80;
t101 = rSges(5,2) * t108 + t118 * t68 + t119 * t69;
t100 = rSges(5,2) * t107 + t118 * t70 + t119 * t71;
t99 = -t82 * rSges(5,2) + (t118 * t79 + t119 * t81) * t80;
t65 = -t82 * rSges(4,3) + (rSges(4,1) * t81 - rSges(4,2) * t79) * t80;
t75 = t80 * pkin(2) - t82 * pkin(5);
t98 = -t65 - t75;
t97 = t115 * (pkin(2) * t82 + pkin(5) * t80);
t90 = -t75 - t99;
t83 = Icges(3,5) * t82 - Icges(3,6) * t80;
t74 = t80 * rSges(3,1) + t82 * rSges(3,2);
t53 = Icges(3,3) * t77 + t78 * t83;
t52 = -Icges(3,3) * t78 + t77 * t83;
t49 = t98 * t78;
t48 = t98 * t77;
t47 = t71 * rSges(4,1) - t70 * rSges(4,2) + rSges(4,3) * t107;
t45 = t69 * rSges(4,1) - t68 * rSges(4,2) + rSges(4,3) * t108;
t31 = t115 * (rSges(3,1) * t82 - rSges(3,2) * t80);
t30 = t90 * t78;
t29 = t90 * t77;
t28 = -t65 * t107 - t82 * t47;
t27 = t65 * t108 + t82 * t45;
t26 = (t45 * t78 - t47 * t77) * t80;
t25 = t77 * t45 + t78 * t47 + t97;
t24 = -t100 * t82 - t99 * t107;
t23 = t101 * t82 + t99 * t108;
t22 = -t82 * t35 + (-t39 * t79 + t43 * t81) * t80;
t21 = -t82 * t34 + (-t38 * t79 + t42 * t81) * t80;
t20 = -t82 * t37 + (t33 * t79 + t41 * t81) * t80;
t19 = -t82 * t36 + (t32 * t79 + t40 * t81) * t80;
t10 = (-t100 * t77 + t101 * t78) * t80;
t9 = t100 * t78 + t101 * t77 + t97;
t8 = -t17 * t78 + t18 * t77;
t7 = -t15 * t78 + t16 * t77;
t6 = -t13 * t78 + t14 * t77;
t5 = -t11 * t78 + t12 * t77;
t1 = [m(2) + m(3) + m(4) + m(5); m(3) * t31 + m(4) * t25 + m(5) * t9; m(5) * (t29 ^ 2 + t30 ^ 2 + t9 ^ 2) + m(4) * (t25 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(3) * (t115 * t74 ^ 2 + t31 ^ 2) + (t114 * t53 + t7 + t8) * t77 + (-t113 * t52 - t5 - t6 + (-t77 * t52 + t78 * t53) * t77) * t78; m(4) * t26 + m(5) * t10; m(5) * (t10 * t9 + t23 * t30 + t24 * t29) + m(4) * (t26 * t25 + t27 * t49 + t28 * t48) + ((t8 / 0.2e1 + t7 / 0.2e1) * t78 + (t6 / 0.2e1 + t5 / 0.2e1) * t77) * t80 + t116 * t77 / 0.2e1 - t117 * t78 / 0.2e1 - ((-t19 - t21) * t78 + (t20 + t22) * t77) * t82 / 0.2e1; m(4) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(5) * (t10 ^ 2 + t23 ^ 2 + t24 ^ 2) - t82 * (t112 * t59 + (t22 * t78 + t21 * t77 - (-t61 * t79 + t63 * t81) * t82) * t80) - t82 * (t112 * t60 + (t20 * t78 + t19 * t77 - (t58 * t79 + t62 * t81) * t82) * t80) + t117 * t108 + t116 * t107; m(5) * t106; m(5) * (t9 * t106 + t68 * t29 + t70 * t30); m(5) * (t10 * t106 + t70 * t23 + t68 * t24); m(5) * (t80 ^ 2 * t79 ^ 2 + t68 ^ 2 + t70 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
