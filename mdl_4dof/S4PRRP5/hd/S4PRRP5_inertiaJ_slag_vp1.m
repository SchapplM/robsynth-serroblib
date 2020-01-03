% Calculate joint inertia matrix for
% S4PRRP5
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP5_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP5_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:51
% EndTime: 2019-12-31 16:28:54
% DurationCPUTime: 1.03s
% Computational Cost: add. (1339->181), mult. (3155->294), div. (0->0), fcn. (3438->6), ass. (0->95)
t83 = sin(qJ(3));
t84 = sin(qJ(2));
t85 = cos(qJ(3));
t86 = cos(qJ(2));
t61 = -Icges(5,6) * t86 + (Icges(5,4) * t85 - Icges(5,2) * t83) * t84;
t62 = -Icges(4,6) * t86 + (Icges(4,4) * t85 - Icges(4,2) * t83) * t84;
t126 = -t61 - t62;
t63 = -Icges(5,5) * t86 + (Icges(5,1) * t85 - Icges(5,4) * t83) * t84;
t64 = -Icges(4,5) * t86 + (Icges(4,1) * t85 - Icges(4,4) * t83) * t84;
t125 = -t63 - t64;
t59 = -Icges(5,3) * t86 + (Icges(5,5) * t85 - Icges(5,6) * t83) * t84;
t60 = -Icges(4,3) * t86 + (Icges(4,5) * t85 - Icges(4,6) * t83) * t84;
t124 = (-t60 - t59) * t86;
t80 = sin(pkin(6));
t77 = t80 ^ 2;
t81 = cos(pkin(6));
t78 = t81 ^ 2;
t102 = t77 + t78;
t114 = t80 * t84;
t111 = t83 * t86;
t69 = -t111 * t80 - t81 * t85;
t110 = t85 * t86;
t113 = t81 * t83;
t70 = t110 * t80 - t113;
t96 = Icges(5,3) * t84;
t32 = Icges(5,5) * t70 + Icges(5,6) * t69 + t80 * t96;
t98 = Icges(5,6) * t84;
t36 = Icges(5,4) * t70 + Icges(5,2) * t69 + t80 * t98;
t100 = Icges(5,5) * t84;
t40 = Icges(5,1) * t70 + Icges(5,4) * t69 + t100 * t80;
t11 = t114 * t32 + t36 * t69 + t40 * t70;
t71 = -t111 * t81 + t80 * t85;
t115 = t80 * t83;
t72 = t110 * t81 + t115;
t33 = Icges(5,5) * t72 + Icges(5,6) * t71 + t81 * t96;
t37 = Icges(5,4) * t72 + Icges(5,2) * t71 + t81 * t98;
t41 = Icges(5,1) * t72 + Icges(5,4) * t71 + t100 * t81;
t12 = t114 * t33 + t37 * t69 + t41 * t70;
t97 = Icges(4,3) * t84;
t34 = Icges(4,5) * t70 + Icges(4,6) * t69 + t80 * t97;
t99 = Icges(4,6) * t84;
t38 = Icges(4,4) * t70 + Icges(4,2) * t69 + t80 * t99;
t101 = Icges(4,5) * t84;
t42 = Icges(4,1) * t70 + Icges(4,4) * t69 + t101 * t80;
t13 = t114 * t34 + t38 * t69 + t42 * t70;
t35 = Icges(4,5) * t72 + Icges(4,6) * t71 + t81 * t97;
t39 = Icges(4,4) * t72 + Icges(4,2) * t71 + t81 * t99;
t43 = Icges(4,1) * t72 + Icges(4,4) * t71 + t101 * t81;
t14 = t114 * t35 + t39 * t69 + t43 * t70;
t123 = (t125 * t70 + t126 * t69) * t86 + ((t12 + t14) * t81 + (t11 + t13 + t124) * t80) * t84;
t112 = t81 * t84;
t15 = t112 * t32 + t36 * t71 + t40 * t72;
t16 = t112 * t33 + t37 * t71 + t41 * t72;
t17 = t112 * t34 + t38 * t71 + t42 * t72;
t18 = t112 * t35 + t39 * t71 + t43 * t72;
t122 = (t125 * t72 + t126 * t71) * t86 + ((t16 + t18 + t124) * t81 + (t15 + t17) * t80) * t84;
t121 = t86 ^ 2;
t117 = t85 * pkin(3);
t87 = qJ(4) * t84 + t117 * t86;
t107 = rSges(5,1) * t70 + rSges(5,2) * t69 + rSges(5,3) * t114 - pkin(3) * t113 + t80 * t87;
t106 = rSges(5,1) * t72 + rSges(5,2) * t71 + rSges(5,3) * t112 + pkin(3) * t115 + t81 * t87;
t105 = (-qJ(4) - rSges(5,3)) * t86 + (rSges(5,1) * t85 - rSges(5,2) * t83 + t117) * t84;
t66 = -t86 * rSges(4,3) + (rSges(4,1) * t85 - rSges(4,2) * t83) * t84;
t75 = t84 * pkin(2) - t86 * pkin(5);
t104 = -t66 - t75;
t103 = t102 * (pkin(2) * t86 + pkin(5) * t84);
t95 = -t75 - t105;
t88 = Icges(3,5) * t86 - Icges(3,6) * t84;
t74 = rSges(3,1) * t84 + rSges(3,2) * t86;
t54 = Icges(3,3) * t80 + t81 * t88;
t53 = -Icges(3,3) * t81 + t80 * t88;
t51 = t104 * t81;
t50 = t104 * t80;
t49 = rSges(4,1) * t72 + rSges(4,2) * t71 + rSges(4,3) * t112;
t47 = rSges(4,1) * t70 + rSges(4,2) * t69 + rSges(4,3) * t114;
t31 = t102 * (rSges(3,1) * t86 - rSges(3,2) * t84);
t30 = t95 * t81;
t29 = t95 * t80;
t28 = -t112 * t66 - t49 * t86;
t27 = t114 * t66 + t47 * t86;
t26 = (t47 * t81 - t49 * t80) * t84;
t25 = t47 * t80 + t49 * t81 + t103;
t24 = -t86 * t35 + (-t39 * t83 + t43 * t85) * t84;
t23 = -t86 * t34 + (-t38 * t83 + t42 * t85) * t84;
t22 = -t86 * t33 + (-t37 * t83 + t41 * t85) * t84;
t21 = -t86 * t32 + (-t36 * t83 + t40 * t85) * t84;
t20 = -t105 * t112 - t106 * t86;
t19 = t105 * t114 + t107 * t86;
t10 = (-t106 * t80 + t107 * t81) * t84;
t9 = t106 * t81 + t107 * t80 + t103;
t8 = -t17 * t81 + t18 * t80;
t7 = -t15 * t81 + t16 * t80;
t6 = -t13 * t81 + t14 * t80;
t5 = -t11 * t81 + t12 * t80;
t1 = [m(2) + m(3) + m(4) + m(5); m(3) * t31 + m(4) * t25 + m(5) * t9; m(5) * (t29 ^ 2 + t30 ^ 2 + t9 ^ 2) + m(4) * (t25 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(3) * (t102 * t74 ^ 2 + t31 ^ 2) + (t77 * t54 + t7 + t8) * t80 + (-t78 * t53 - t5 - t6 + (-t80 * t53 + t81 * t54) * t80) * t81; m(4) * t26 + m(5) * t10; m(5) * (t10 * t9 + t19 * t30 + t20 * t29) + m(4) * (t25 * t26 + t27 * t51 + t28 * t50) + ((t8 / 0.2e1 + t7 / 0.2e1) * t81 + (t6 / 0.2e1 + t5 / 0.2e1) * t80) * t84 + t122 * t80 / 0.2e1 - t123 * t81 / 0.2e1 - ((-t21 - t23) * t81 + (t22 + t24) * t80) * t86 / 0.2e1; m(5) * (t10 ^ 2 + t19 ^ 2 + t20 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) - t86 * (t121 * t60 + (t24 * t81 + t23 * t80 - (-t62 * t83 + t64 * t85) * t86) * t84) - t86 * (t121 * t59 + (t22 * t81 + t21 * t80 - (-t61 * t83 + t63 * t85) * t86) * t84) + t123 * t114 + t122 * t112; -m(5) * t86; m(5) * (-t86 * t9 + (t29 * t80 + t30 * t81) * t84); m(5) * (-t86 * t10 + (t19 * t81 + t20 * t80) * t84); m(5) * (t102 * t84 ^ 2 + t121);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
