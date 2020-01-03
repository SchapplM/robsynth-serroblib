% Calculate joint inertia matrix for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:07
% EndTime: 2019-12-31 21:49:09
% DurationCPUTime: 0.72s
% Computational Cost: add. (2198->137), mult. (1570->191), div. (0->0), fcn. (1336->8), ass. (0->76)
t90 = cos(qJ(4));
t142 = t90 ^ 2;
t141 = Icges(6,4) + Icges(5,5);
t140 = Icges(5,6) - Icges(6,6);
t88 = sin(qJ(4));
t138 = t141 * t88;
t139 = t140 * t90;
t133 = t138 + t139;
t87 = qJ(1) + qJ(2);
t84 = qJ(3) + t87;
t80 = sin(t84);
t76 = t80 ^ 2;
t81 = cos(t84);
t77 = t81 ^ 2;
t137 = -t140 * t88 + t141 * t90;
t135 = Icges(6,2) + Icges(5,3);
t129 = rSges(6,1) + pkin(4);
t132 = t129 * t90;
t131 = rSges(6,3) + qJ(5);
t128 = t135 * t81 - t137 * t80;
t127 = t135 * t80 + t137 * t81;
t65 = t88 * rSges(5,1) + t90 * rSges(5,2);
t124 = m(5) * t65;
t123 = m(6) * t88;
t82 = sin(t87);
t122 = pkin(2) * t82;
t89 = sin(qJ(1));
t121 = t89 * pkin(1);
t120 = rSges(5,1) * t90;
t119 = t81 * t88;
t118 = t81 * t90;
t117 = t80 * t88 * rSges(5,2) + t81 * rSges(5,3);
t116 = -t129 * t88 + t131 * t90;
t115 = t81 * pkin(3) + t80 * pkin(8);
t114 = t76 + t77;
t109 = qJ(5) * t88;
t83 = cos(t87);
t46 = t83 * rSges(3,1) - t82 * rSges(3,2);
t44 = t81 * rSges(4,1) - t80 * rSges(4,2);
t79 = pkin(2) * t83;
t40 = t44 + t79;
t108 = t80 * rSges(6,2) + rSges(6,3) * t119 + t81 * t109 + t129 * t118;
t45 = -t82 * rSges(3,1) - t83 * rSges(3,2);
t43 = -t80 * rSges(4,1) - t81 * rSges(4,2);
t95 = rSges(5,1) * t118 - rSges(5,2) * t119 + t80 * rSges(5,3);
t94 = Icges(4,3) + (Icges(5,2) + Icges(6,3)) * t142 + ((Icges(5,1) + Icges(6,1)) * t88 + (2 * Icges(5,4) - 2 * Icges(6,5)) * t90) * t88;
t12 = t108 + t115;
t93 = t133 * t77 + (t139 / 0.2e1 + t138 / 0.2e1 + t133 / 0.2e1) * t76;
t22 = t95 + t115;
t39 = t43 - t122;
t10 = t79 + t12;
t92 = Icges(3,3) + t94;
t20 = t22 + t79;
t74 = t81 * pkin(8);
t21 = t74 + (-pkin(3) - t120) * t80 + t117;
t19 = t21 - t122;
t71 = t81 * rSges(6,2);
t11 = t71 + t74 + (-t131 * t88 - pkin(3) - t132) * t80;
t9 = t11 - t122;
t91 = cos(qJ(1));
t85 = t91 * pkin(1);
t67 = t91 * rSges(2,1) - t89 * rSges(2,2);
t66 = -t89 * rSges(2,1) - t91 * rSges(2,2);
t42 = t46 + t85;
t41 = t45 - t121;
t38 = t40 + t85;
t37 = t39 - t121;
t24 = t116 * t81;
t23 = t116 * t80;
t18 = t20 + t85;
t17 = t19 - t121;
t8 = t80 * (t120 * t80 - t117) + t81 * t95;
t3 = t85 + t10;
t2 = t9 - t121;
t1 = t108 * t81 + (-t71 + (rSges(6,3) * t88 + t109 + t132) * t80) * t80;
t4 = [Icges(2,3) + m(6) * (t2 ^ 2 + t3 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + m(3) * (t41 ^ 2 + t42 ^ 2) + m(2) * (t66 ^ 2 + t67 ^ 2) + t92; m(6) * (t10 * t3 + t9 * t2) + m(5) * (t19 * t17 + t20 * t18) + m(4) * (t39 * t37 + t40 * t38) + m(3) * (t45 * t41 + t46 * t42) + t92; m(6) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2) + m(4) * (t39 ^ 2 + t40 ^ 2) + m(3) * (t45 ^ 2 + t46 ^ 2) + t92; m(6) * (t11 * t2 + t12 * t3) + m(5) * (t21 * t17 + t22 * t18) + m(4) * (t43 * t37 + t44 * t38) + t94; m(6) * (t12 * t10 + t11 * t9) + m(5) * (t21 * t19 + t22 * t20) + m(4) * (t43 * t39 + t44 * t40) + t94; m(6) * (t11 ^ 2 + t12 ^ 2) + m(5) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t43 ^ 2 + t44 ^ 2) + t94; m(6) * (t24 * t2 + t23 * t3) + (-t17 * t81 - t18 * t80) * t124 + t93; m(6) * (t23 * t10 + t24 * t9) + (-t19 * t81 - t20 * t80) * t124 + t93; m(6) * (t24 * t11 + t23 * t12) + (-t21 * t81 - t22 * t80) * t124 + t93; m(5) * (t114 * t65 ^ 2 + t8 ^ 2) + m(6) * (t1 ^ 2 + t23 ^ 2 + t24 ^ 2) + t128 * t81 * t77 + (t127 * t76 + (t127 * t81 + t128 * t80) * t81) * t80; (t2 * t81 + t3 * t80) * t123; (t10 * t80 + t81 * t9) * t123; (t11 * t81 + t12 * t80) * t123; m(6) * (-t90 * t1 + (t23 * t80 + t24 * t81) * t88); m(6) * (t114 * t88 ^ 2 + t142);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
