% Calculate joint inertia matrix for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:51
% EndTime: 2019-03-09 01:29:53
% DurationCPUTime: 0.86s
% Computational Cost: add. (2080->197), mult. (2433->301), div. (0->0), fcn. (2492->8), ass. (0->101)
t86 = cos(qJ(5));
t126 = Icges(6,5) * t86;
t125 = t126 / 0.2e1;
t124 = m(5) + m(7);
t81 = qJ(1) + pkin(9);
t78 = sin(t81);
t76 = t78 ^ 2;
t79 = cos(t81);
t77 = t79 ^ 2;
t60 = t76 + t77;
t123 = -t78 / 0.2e1;
t83 = sin(qJ(5));
t122 = t83 / 0.2e1;
t58 = m(6) * t60;
t121 = pkin(5) * t83;
t84 = sin(qJ(1));
t120 = t84 * pkin(1);
t119 = -rSges(6,3) - pkin(7);
t118 = rSges(7,3) + pkin(8);
t117 = t78 * t86;
t116 = t79 * t83;
t115 = t79 * t86;
t82 = sin(qJ(6));
t85 = cos(qJ(6));
t52 = Icges(7,6) * t83 + (Icges(7,4) * t85 - Icges(7,2) * t82) * t86;
t114 = t82 * t52;
t113 = t82 * t83;
t112 = t83 * t85;
t111 = -pkin(2) - qJ(4);
t51 = Icges(7,3) * t83 + (Icges(7,5) * t85 - Icges(7,6) * t82) * t86;
t53 = Icges(7,5) * t83 + (Icges(7,1) * t85 - Icges(7,4) * t82) * t86;
t110 = t86 * t85 * t53 + t83 * t51;
t49 = -t79 * t113 - t78 * t85;
t50 = t79 * t112 - t78 * t82;
t109 = t50 * rSges(7,1) + t49 * rSges(7,2);
t54 = t83 * rSges(7,3) + (rSges(7,1) * t85 - rSges(7,2) * t82) * t86;
t108 = t86 * pkin(5) + t83 * pkin(8) + t54;
t107 = rSges(6,1) * t116 + rSges(6,2) * t115;
t106 = Icges(6,4) * t83;
t104 = Icges(7,5) * t86;
t103 = Icges(7,6) * t86;
t102 = Icges(7,3) * t86;
t101 = t124 * t60 + t58;
t87 = cos(qJ(1));
t80 = t87 * pkin(1);
t100 = t79 * pkin(2) + t78 * qJ(3) + t80;
t47 = -t78 * t113 + t79 * t85;
t48 = t78 * t112 + t79 * t82;
t13 = -t51 * t117 + t47 * t52 + t48 * t53;
t23 = Icges(7,5) * t48 + Icges(7,6) * t47 - t78 * t102;
t25 = Icges(7,4) * t48 + Icges(7,2) * t47 - t78 * t103;
t27 = Icges(7,1) * t48 + Icges(7,4) * t47 - t78 * t104;
t9 = t83 * t23 + (-t25 * t82 + t27 * t85) * t86;
t99 = -t9 / 0.2e1 - t13 / 0.2e1;
t24 = Icges(7,5) * t50 + Icges(7,6) * t49 - t79 * t102;
t26 = Icges(7,4) * t50 + Icges(7,2) * t49 - t79 * t103;
t28 = Icges(7,1) * t50 + Icges(7,4) * t49 - t79 * t104;
t10 = t83 * t24 + (-t26 * t82 + t28 * t85) * t86;
t14 = -t51 * t115 + t49 * t52 + t50 * t53;
t98 = -t10 / 0.2e1 - t14 / 0.2e1;
t97 = t79 * qJ(3) - t120;
t96 = t79 * qJ(4) + t100;
t95 = -rSges(6,1) * t83 - rSges(6,2) * t86;
t94 = -t48 * rSges(7,1) - t47 * rSges(7,2);
t90 = Icges(6,2) * t86 + t106;
t89 = Icges(6,5) * t83 + Icges(6,6) * t86;
t21 = t119 * t79 + (t95 + t111) * t78 + t97;
t22 = t119 * t78 + t107 + t96;
t88 = m(6) * (t79 * t21 + t78 * t22);
t71 = pkin(5) * t116;
t67 = t87 * rSges(2,1) - t84 * rSges(2,2);
t66 = t86 * rSges(6,1) - t83 * rSges(6,2);
t65 = -t84 * rSges(2,1) - t87 * rSges(2,2);
t56 = t79 * rSges(3,1) - t78 * rSges(3,2) + t80;
t55 = -t78 * rSges(3,1) - t79 * rSges(3,2) - t120;
t37 = Icges(6,3) * t79 + t89 * t78;
t36 = -t79 * rSges(4,2) + t78 * rSges(4,3) + t100;
t35 = t79 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t78 + t97;
t34 = t78 * rSges(5,2) + t79 * rSges(5,3) + t96;
t33 = t79 * rSges(5,2) + (-rSges(5,3) + t111) * t78 + t97;
t32 = t108 * t79;
t31 = t108 * t78;
t30 = -rSges(7,3) * t115 + t109;
t29 = -rSges(7,3) * t117 - t94;
t20 = -t79 * t107 + t95 * t76;
t19 = (-t86 * t114 + t110) * t83;
t18 = t54 * t115 + t83 * t30;
t17 = -t54 * t117 - t83 * t29;
t16 = -t78 * pkin(7) - t118 * t115 + t109 + t71 + t96;
t15 = -t79 * pkin(7) + (t118 * t86 + t111 - t121) * t78 + t94 + t97;
t12 = (-t29 * t79 + t30 * t78) * t86;
t11 = (pkin(8) * t115 - t30 - t71) * t79 + (-t29 + (pkin(8) * t86 - t121) * t78) * t78;
t8 = -t24 * t115 + t49 * t26 + t50 * t28;
t7 = -t23 * t115 + t49 * t25 + t50 * t27;
t6 = -t24 * t117 + t47 * t26 + t48 * t28;
t5 = -t23 * t117 + t47 * t25 + t48 * t27;
t4 = t7 * t79 - t8 * t78;
t3 = t5 * t79 - t6 * t78;
t2 = t14 * t83 + (-t7 * t78 - t79 * t8) * t86;
t1 = t13 * t83 + (-t5 * t78 - t6 * t79) * t86;
t38 = [-t83 * (Icges(6,4) * t86 - Icges(6,2) * t83) + Icges(4,1) + Icges(5,1) + Icges(2,3) + Icges(3,3) + (Icges(6,1) * t86 - t106 - t114) * t86 + m(7) * (t15 ^ 2 + t16 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2) + m(5) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t55 ^ 2 + t56 ^ 2) + m(2) * (t65 ^ 2 + t67 ^ 2) + t110; 0; m(3) + m(4) + m(6) + t124; m(7) * (t78 * t15 - t79 * t16) + m(6) * (t78 * t21 - t79 * t22) + m(4) * (t78 * t35 - t79 * t36) + m(5) * (t78 * t33 - t79 * t34); 0; m(4) * t60 + t101; m(7) * (t79 * t15 + t78 * t16) + t88 + m(5) * (t79 * t33 + t78 * t34); 0; 0; t101; (-t83 * (Icges(6,6) * t79 + t90 * t78) / 0.2e1 + t79 * t125 - t99) * t79 + ((-Icges(6,6) * t78 + t90 * t79) * t122 + t78 * t125 + t98) * t78 + m(7) * (t32 * t15 + t31 * t16) + t66 * t88 + (t77 / 0.2e1 + t76 / 0.2e1) * (-Icges(6,6) * t83 + t126); m(6) * t20 + m(7) * t11; m(7) * (-t31 * t79 + t32 * t78); m(7) * (t31 * t78 + t32 * t79) + t66 * t58; m(6) * (t60 * t66 ^ 2 + t20 ^ 2) + m(7) * (t11 ^ 2 + t31 ^ 2 + t32 ^ 2) + (t77 * t37 + t3) * t79 + (t78 * t37 * t79 - t4 - t60 * (-Icges(6,3) * t78 + t89 * t79)) * t78; m(7) * (t17 * t15 + t18 * t16) + t19 + (t99 * t78 + t98 * t79) * t86; m(7) * t12; m(7) * (t17 * t78 - t18 * t79); m(7) * (t17 * t79 + t18 * t78); t2 * t123 + t79 * t1 / 0.2e1 + (-t10 * t78 + t9 * t79) * t122 + m(7) * (t12 * t11 + t17 * t32 + t18 * t31) + (-t79 * t4 / 0.2e1 + t3 * t123) * t86; t83 * t19 + m(7) * (t12 ^ 2 + t17 ^ 2 + t18 ^ 2) + (-t79 * t2 - t78 * t1 + t83 * (-t10 * t79 - t78 * t9)) * t86;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t38(1) t38(2) t38(4) t38(7) t38(11) t38(16); t38(2) t38(3) t38(5) t38(8) t38(12) t38(17); t38(4) t38(5) t38(6) t38(9) t38(13) t38(18); t38(7) t38(8) t38(9) t38(10) t38(14) t38(19); t38(11) t38(12) t38(13) t38(14) t38(15) t38(20); t38(16) t38(17) t38(18) t38(19) t38(20) t38(21);];
Mq  = res;
