% Calculate joint inertia matrix for
% S5PPRRP2
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:33
% EndTime: 2019-12-05 15:08:37
% DurationCPUTime: 1.08s
% Computational Cost: add. (2417->188), mult. (3300->308), div. (0->0), fcn. (3614->6), ass. (0->97)
t82 = pkin(8) + qJ(3);
t78 = sin(t82);
t79 = cos(t82);
t85 = sin(qJ(4));
t86 = cos(qJ(4));
t52 = -Icges(6,6) * t79 + (Icges(6,5) * t86 + Icges(6,3) * t85) * t78;
t55 = -Icges(5,6) * t79 + (Icges(5,4) * t86 - Icges(5,2) * t85) * t78;
t126 = -t52 + t55;
t56 = -Icges(6,4) * t79 + (Icges(6,1) * t86 + Icges(6,5) * t85) * t78;
t57 = -Icges(5,5) * t79 + (Icges(5,1) * t86 - Icges(5,4) * t85) * t78;
t125 = -t56 - t57;
t53 = -Icges(5,3) * t79 + (Icges(5,5) * t86 - Icges(5,6) * t85) * t78;
t54 = -Icges(6,2) * t79 + (Icges(6,4) * t86 + Icges(6,6) * t85) * t78;
t124 = (-t53 - t54) * t79;
t123 = rSges(6,1) + pkin(4);
t83 = sin(pkin(7));
t80 = t83 ^ 2;
t84 = cos(pkin(7));
t81 = t84 ^ 2;
t101 = t80 + t81;
t122 = rSges(6,3) + qJ(5);
t115 = t78 * t83;
t107 = t84 * t86;
t110 = t83 * t85;
t69 = t79 * t110 + t107;
t108 = t84 * t85;
t109 = t83 * t86;
t70 = t79 * t109 - t108;
t96 = Icges(6,6) * t78;
t32 = Icges(6,5) * t70 + Icges(6,3) * t69 + t83 * t96;
t98 = Icges(6,2) * t78;
t36 = Icges(6,4) * t70 + Icges(6,6) * t69 + t83 * t98;
t100 = Icges(6,4) * t78;
t40 = Icges(6,1) * t70 + Icges(6,5) * t69 + t83 * t100;
t11 = t36 * t115 + t69 * t32 + t70 * t40;
t71 = t79 * t108 - t109;
t72 = t79 * t107 + t110;
t33 = Icges(6,5) * t72 + Icges(6,3) * t71 + t84 * t96;
t37 = Icges(6,4) * t72 + Icges(6,6) * t71 + t84 * t98;
t41 = Icges(6,1) * t72 + Icges(6,5) * t71 + t84 * t100;
t12 = t37 * t115 + t69 * t33 + t70 * t41;
t95 = Icges(5,3) * t78;
t34 = Icges(5,5) * t70 - Icges(5,6) * t69 + t83 * t95;
t97 = Icges(5,6) * t78;
t38 = Icges(5,4) * t70 - Icges(5,2) * t69 + t83 * t97;
t99 = Icges(5,5) * t78;
t42 = Icges(5,1) * t70 - Icges(5,4) * t69 + t83 * t99;
t13 = t34 * t115 - t69 * t38 + t70 * t42;
t35 = Icges(5,5) * t72 - Icges(5,6) * t71 + t84 * t95;
t39 = Icges(5,4) * t72 - Icges(5,2) * t71 + t84 * t97;
t43 = Icges(5,1) * t72 - Icges(5,4) * t71 + t84 * t99;
t14 = t35 * t115 - t69 * t39 + t70 * t43;
t121 = (t125 * t70 + t126 * t69) * t79 + ((t12 + t14) * t84 + (t11 + t13 + t124) * t83) * t78;
t114 = t78 * t84;
t15 = t36 * t114 + t71 * t32 + t72 * t40;
t16 = t37 * t114 + t71 * t33 + t72 * t41;
t17 = t34 * t114 - t71 * t38 + t72 * t42;
t18 = t35 * t114 - t71 * t39 + t72 * t43;
t120 = (t125 * t72 + t126 * t71) * t79 + ((t16 + t18 + t124) * t84 + (t15 + t17) * t83) * t78;
t119 = t79 ^ 2;
t113 = t78 * t85;
t106 = rSges(6,2) * t115 + t122 * t69 + t123 * t70;
t105 = rSges(6,2) * t114 + t122 * t71 + t123 * t72;
t104 = -t79 * rSges(6,2) + (t122 * t85 + t123 * t86) * t78;
t59 = -t79 * rSges(5,3) + (rSges(5,1) * t86 - rSges(5,2) * t85) * t78;
t75 = t78 * pkin(3) - t79 * pkin(6);
t103 = -t59 - t75;
t102 = t101 * (pkin(3) * t79 + pkin(6) * t78);
t94 = -t75 - t104;
t87 = Icges(4,5) * t79 - Icges(4,6) * t78;
t74 = t78 * rSges(4,1) + t79 * rSges(4,2);
t61 = Icges(4,3) * t83 + t87 * t84;
t60 = -Icges(4,3) * t84 + t87 * t83;
t49 = t103 * t84;
t48 = t103 * t83;
t47 = t72 * rSges(5,1) - t71 * rSges(5,2) + rSges(5,3) * t114;
t45 = t70 * rSges(5,1) - t69 * rSges(5,2) + rSges(5,3) * t115;
t31 = t101 * (rSges(4,1) * t79 - rSges(4,2) * t78);
t30 = t94 * t84;
t29 = t94 * t83;
t28 = -t59 * t114 - t79 * t47;
t27 = t59 * t115 + t79 * t45;
t26 = (t45 * t84 - t47 * t83) * t78;
t25 = t83 * t45 + t84 * t47 + t102;
t24 = -t104 * t114 - t105 * t79;
t23 = t104 * t115 + t106 * t79;
t22 = -t79 * t35 + (-t39 * t85 + t43 * t86) * t78;
t21 = -t79 * t34 + (-t38 * t85 + t42 * t86) * t78;
t20 = -t79 * t37 + (t33 * t85 + t41 * t86) * t78;
t19 = -t79 * t36 + (t32 * t85 + t40 * t86) * t78;
t10 = (-t105 * t83 + t106 * t84) * t78;
t9 = t105 * t84 + t106 * t83 + t102;
t8 = -t17 * t84 + t18 * t83;
t7 = -t15 * t84 + t16 * t83;
t6 = -t13 * t84 + t14 * t83;
t5 = -t11 * t84 + t12 * t83;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t101; m(4) * t31 + m(5) * t25 + m(6) * t9; m(5) * (-t48 * t84 + t49 * t83) + m(6) * (-t29 * t84 + t30 * t83); m(6) * (t29 ^ 2 + t30 ^ 2 + t9 ^ 2) + m(5) * (t25 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t101 * t74 ^ 2 + t31 ^ 2) + (t80 * t61 + t7 + t8) * t83 + (-t81 * t60 - t5 - t6 + (-t83 * t60 + t84 * t61) * t83) * t84; m(5) * t26 + m(6) * t10; m(5) * (t27 * t83 - t28 * t84) + m(6) * (t23 * t83 - t24 * t84); m(6) * (t10 * t9 + t23 * t30 + t24 * t29) + m(5) * (t26 * t25 + t27 * t49 + t28 * t48) + ((t8 / 0.2e1 + t7 / 0.2e1) * t84 + (t6 / 0.2e1 + t5 / 0.2e1) * t83) * t78 - ((-t19 - t21) * t84 + (t20 + t22) * t83) * t79 / 0.2e1 + t120 * t83 / 0.2e1 - t121 * t84 / 0.2e1; m(6) * (t10 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) - t79 * (t119 * t53 + (t22 * t84 + t21 * t83 - (-t55 * t85 + t57 * t86) * t79) * t78) - t79 * (t119 * t54 + (t20 * t84 + t19 * t83 - (t52 * t85 + t56 * t86) * t79) * t78) + t121 * t115 + t120 * t114; m(6) * t113; m(6) * (-t69 * t84 + t71 * t83); m(6) * (t9 * t113 + t69 * t29 + t71 * t30); m(6) * (t10 * t113 + t71 * t23 + t69 * t24); m(6) * (t78 ^ 2 * t85 ^ 2 + t69 ^ 2 + t71 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
