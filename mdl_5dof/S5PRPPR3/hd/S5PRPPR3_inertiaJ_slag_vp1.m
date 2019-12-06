% Calculate joint inertia matrix for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:19
% EndTime: 2019-12-05 15:26:22
% DurationCPUTime: 0.91s
% Computational Cost: add. (1594->138), mult. (2192->236), div. (0->0), fcn. (2256->8), ass. (0->78)
t76 = sin(pkin(7));
t73 = t76 ^ 2;
t77 = cos(pkin(7));
t74 = t77 ^ 2;
t111 = t73 + t74;
t75 = qJ(2) + pkin(8);
t71 = sin(t75);
t72 = cos(t75);
t80 = sin(qJ(2));
t82 = cos(qJ(2));
t129 = Icges(3,5) * t82 - Icges(3,6) * t80 + (-Icges(5,4) + Icges(4,5)) * t72 + (Icges(5,5) - Icges(4,6)) * t71;
t128 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t127 = t128 * t77 - t129 * t76;
t126 = t128 * t76 + t129 * t77;
t125 = t71 ^ 2;
t124 = t76 / 0.2e1;
t123 = -m(5) - m(6);
t122 = pkin(2) * t80;
t79 = sin(qJ(5));
t81 = cos(qJ(5));
t32 = Icges(6,3) * t71 + (-Icges(6,5) * t79 - Icges(6,6) * t81) * t72;
t119 = t71 * t32;
t118 = t72 * t76;
t117 = t72 * t77;
t116 = t76 * t79;
t115 = t76 * t81;
t114 = t77 * t79;
t113 = t77 * t81;
t112 = t111 * t82 * pkin(2);
t110 = Icges(6,5) * t72;
t109 = Icges(6,6) * t72;
t108 = Icges(6,3) * t72;
t107 = m(5) / 0.2e1 + m(6) / 0.2e1;
t106 = -t71 * pkin(3) + t72 * qJ(4) - t122;
t105 = -t71 * rSges(4,1) - t72 * rSges(4,2) - t122;
t104 = t112 + t111 * (pkin(3) * t72 + qJ(4) * t71);
t103 = t71 * rSges(5,2) + t72 * rSges(5,3) + t106;
t35 = t71 * rSges(6,3) + (-rSges(6,1) * t79 - rSges(6,2) * t81) * t72;
t83 = -pkin(6) * t71 + t106 - t35;
t67 = t80 * rSges(3,1) + t82 * rSges(3,2);
t62 = t71 * t116 - t113;
t61 = t71 * t115 + t114;
t60 = t71 * t114 + t115;
t59 = t71 * t113 - t116;
t37 = t105 * t77;
t36 = t105 * t76;
t34 = Icges(6,5) * t71 + (-Icges(6,1) * t79 - Icges(6,4) * t81) * t72;
t33 = Icges(6,6) * t71 + (-Icges(6,4) * t79 - Icges(6,2) * t81) * t72;
t29 = t103 * t77;
t28 = t103 * t76;
t27 = t111 * (rSges(3,1) * t82 - rSges(3,2) * t80);
t26 = t62 * rSges(6,1) + t61 * rSges(6,2) + rSges(6,3) * t118;
t25 = t60 * rSges(6,1) + t59 * rSges(6,2) + rSges(6,3) * t117;
t24 = Icges(6,1) * t62 + Icges(6,4) * t61 + t76 * t110;
t23 = Icges(6,1) * t60 + Icges(6,4) * t59 + t77 * t110;
t22 = Icges(6,4) * t62 + Icges(6,2) * t61 + t76 * t109;
t21 = Icges(6,4) * t60 + Icges(6,2) * t59 + t77 * t109;
t20 = Icges(6,5) * t62 + Icges(6,6) * t61 + t76 * t108;
t19 = Icges(6,5) * t60 + Icges(6,6) * t59 + t77 * t108;
t18 = t83 * t77;
t17 = t83 * t76;
t16 = -t35 * t117 + t71 * t25;
t15 = t35 * t118 - t71 * t26;
t14 = t112 + t111 * (rSges(4,1) * t72 - rSges(4,2) * t71);
t13 = (-t25 * t76 + t26 * t77) * t72;
t12 = t104 + t111 * (-rSges(5,2) * t72 + rSges(5,3) * t71);
t11 = t71 * t20 + (-t22 * t81 - t24 * t79) * t72;
t10 = t71 * t19 + (-t21 * t81 - t23 * t79) * t72;
t9 = t20 * t118 + t61 * t22 + t62 * t24;
t8 = t19 * t118 + t61 * t21 + t62 * t23;
t7 = t20 * t117 + t59 * t22 + t60 * t24;
t6 = t19 * t117 + t59 * t21 + t60 * t23;
t5 = t111 * t72 * pkin(6) + t77 * t25 + t76 * t26 + t104;
t4 = t8 * t76 - t9 * t77;
t3 = t6 * t76 - t7 * t77;
t2 = (t61 * t33 + t62 * t34) * t71 + (t8 * t77 + (t9 + t119) * t76) * t72;
t1 = (t59 * t33 + t60 * t34) * t71 + (t7 * t76 + (t6 + t119) * t77) * t72;
t30 = [m(2) + m(3) + m(4) - t123; m(3) * t27 + m(4) * t14 + m(5) * t12 + m(6) * t5; m(6) * (t17 ^ 2 + t18 ^ 2 + t5 ^ 2) + m(5) * (t12 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(4) * (t14 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(3) * (t111 * t67 ^ 2 + t27 ^ 2) + (t127 * t74 - t4) * t77 + (t3 + t126 * t73 + (t126 * t77 + t127 * t76) * t77) * t76; 0; m(6) * (-t77 * t17 + t76 * t18) + m(5) * (-t77 * t28 + t76 * t29) + m(4) * (-t77 * t36 + t76 * t37); 0.2e1 * (m(4) / 0.2e1 + t107) * t111; t123 * t72; m(6) * (-t72 * t5 + (t17 * t76 + t18 * t77) * t71) + m(5) * (-t72 * t12 + (t28 * t76 + t29 * t77) * t71); 0; 0.2e1 * t107 * (t111 * t125 + t72 ^ 2); m(6) * t13; -t77 * t2 / 0.2e1 + m(6) * (t13 * t5 + t15 * t18 + t16 * t17) + t71 * (t10 * t76 - t11 * t77) / 0.2e1 + t1 * t124 + (t77 * t3 / 0.2e1 + t4 * t124) * t72; m(6) * (t15 * t76 - t16 * t77); m(6) * (-t13 * t72 + (t15 * t77 + t16 * t76) * t71); m(6) * (t13 ^ 2 + t15 ^ 2 + t16 ^ 2) + t1 * t117 + t2 * t118 + t71 * (t125 * t32 + (t10 * t77 + t11 * t76 + (-t33 * t81 - t34 * t79) * t71) * t72);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t30(1), t30(2), t30(4), t30(7), t30(11); t30(2), t30(3), t30(5), t30(8), t30(12); t30(4), t30(5), t30(6), t30(9), t30(13); t30(7), t30(8), t30(9), t30(10), t30(14); t30(11), t30(12), t30(13), t30(14), t30(15);];
Mq = res;
