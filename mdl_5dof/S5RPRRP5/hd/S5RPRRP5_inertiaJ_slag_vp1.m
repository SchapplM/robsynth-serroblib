% Calculate joint inertia matrix for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:42
% EndTime: 2019-12-31 18:40:43
% DurationCPUTime: 0.61s
% Computational Cost: add. (1465->111), mult. (1059->152), div. (0->0), fcn. (905->8), ass. (0->67)
t81 = cos(qJ(4));
t133 = t81 ^ 2;
t132 = Icges(6,4) + Icges(5,5);
t131 = Icges(5,6) - Icges(6,6);
t79 = sin(qJ(4));
t129 = t132 * t79;
t130 = t131 * t81;
t124 = t129 + t130;
t77 = qJ(1) + pkin(8);
t75 = qJ(3) + t77;
t71 = sin(t75);
t68 = t71 ^ 2;
t72 = cos(t75);
t69 = t72 ^ 2;
t128 = -t131 * t79 + t132 * t81;
t126 = Icges(6,2) + Icges(5,3);
t120 = rSges(6,1) + pkin(4);
t122 = t120 * t81;
t121 = rSges(6,3) + qJ(5);
t119 = t126 * t72 - t128 * t71;
t118 = t126 * t71 + t128 * t72;
t57 = t79 * rSges(5,1) + t81 * rSges(5,2);
t115 = m(5) * t57;
t114 = m(6) * t79;
t80 = sin(qJ(1));
t113 = t80 * pkin(1);
t112 = rSges(5,1) * t81;
t111 = t72 * t79;
t110 = t72 * t81;
t109 = t71 * t79 * rSges(5,2) + t72 * rSges(5,3);
t108 = -t120 * t79 + t121 * t81;
t107 = t72 * pkin(3) + t71 * pkin(7);
t106 = t68 + t69;
t74 = cos(t77);
t82 = cos(qJ(1));
t76 = t82 * pkin(1);
t105 = pkin(2) * t74 + t76;
t100 = qJ(5) * t79;
t38 = t72 * rSges(4,1) - t71 * rSges(4,2);
t99 = t71 * rSges(6,2) + rSges(6,3) * t111 + t72 * t100 + t120 * t110;
t73 = sin(t77);
t98 = -pkin(2) * t73 - t113;
t37 = -t71 * rSges(4,1) - t72 * rSges(4,2);
t85 = rSges(5,1) * t110 - rSges(5,2) * t111 + t71 * rSges(5,3);
t84 = Icges(4,3) + (Icges(5,2) + Icges(6,3)) * t133 + ((Icges(5,1) + Icges(6,1)) * t79 + (2 * Icges(5,4) - 2 * Icges(6,5)) * t81) * t79;
t10 = t99 + t107;
t83 = t124 * t69 + (t130 / 0.2e1 + t129 / 0.2e1 + t124 / 0.2e1) * t68;
t18 = t85 + t107;
t66 = t72 * pkin(7);
t17 = t66 + (-pkin(3) - t112) * t71 + t109;
t63 = t72 * rSges(6,2);
t9 = t63 + t66 + (-t121 * t79 - pkin(3) - t122) * t71;
t59 = t82 * rSges(2,1) - t80 * rSges(2,2);
t58 = -t80 * rSges(2,1) - t82 * rSges(2,2);
t36 = t74 * rSges(3,1) - t73 * rSges(3,2) + t76;
t35 = -t73 * rSges(3,1) - t74 * rSges(3,2) - t113;
t34 = t38 + t105;
t33 = t37 + t98;
t20 = t108 * t72;
t19 = t108 * t71;
t16 = t18 + t105;
t15 = t17 + t98;
t8 = t71 * (t71 * t112 - t109) + t72 * t85;
t3 = t10 + t105;
t2 = t9 + t98;
t1 = t99 * t72 + (-t63 + (rSges(6,3) * t79 + t100 + t122) * t71) * t71;
t4 = [Icges(2,3) + Icges(3,3) + m(6) * (t2 ^ 2 + t3 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t35 ^ 2 + t36 ^ 2) + m(2) * (t58 ^ 2 + t59 ^ 2) + t84; 0; m(3) + m(4) + m(5) + m(6); m(6) * (t10 * t3 + t9 * t2) + m(5) * (t17 * t15 + t18 * t16) + m(4) * (t37 * t33 + t38 * t34) + t84; 0; m(6) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + t84; m(6) * (t19 * t3 + t20 * t2) + (-t15 * t72 - t16 * t71) * t115 + t83; m(5) * t8 + m(6) * t1; m(6) * (t19 * t10 + t20 * t9) + (-t17 * t72 - t18 * t71) * t115 + t83; m(5) * (t106 * t57 ^ 2 + t8 ^ 2) + m(6) * (t1 ^ 2 + t19 ^ 2 + t20 ^ 2) + t118 * t71 * t68 + (t119 * t69 + (t118 * t72 + t119 * t71) * t71) * t72; (t2 * t72 + t3 * t71) * t114; -m(6) * t81; (t10 * t71 + t72 * t9) * t114; m(6) * (-t81 * t1 + (t19 * t71 + t20 * t72) * t79); m(6) * (t106 * t79 ^ 2 + t133);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
