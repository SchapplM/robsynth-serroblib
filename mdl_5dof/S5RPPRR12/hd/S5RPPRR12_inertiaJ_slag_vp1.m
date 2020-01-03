% Calculate joint inertia matrix for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR12_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR12_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:57
% EndTime: 2019-12-31 18:06:59
% DurationCPUTime: 0.89s
% Computational Cost: add. (1764->184), mult. (2353->292), div. (0->0), fcn. (2424->8), ass. (0->99)
t79 = pkin(8) + qJ(4);
t74 = cos(t79);
t131 = Icges(5,5) * t74;
t73 = sin(t79);
t130 = Icges(5,6) * t73;
t129 = t131 / 0.2e1 - t130 / 0.2e1;
t86 = sin(qJ(1));
t80 = t86 ^ 2;
t88 = cos(qJ(1));
t81 = t88 ^ 2;
t69 = t80 + t81;
t128 = (rSges(5,1) * t73 + rSges(5,2) * t74) * t88;
t125 = t88 / 0.2e1;
t62 = m(5) * t69;
t82 = sin(pkin(8));
t124 = pkin(3) * t82;
t123 = pkin(4) * t73;
t122 = t73 * t86;
t121 = t74 * t86;
t120 = t74 * t88;
t85 = sin(qJ(5));
t87 = cos(qJ(5));
t38 = Icges(6,6) * t73 + (Icges(6,4) * t87 - Icges(6,2) * t85) * t74;
t119 = t85 * t38;
t118 = t86 * t85;
t117 = t86 * t87;
t116 = t88 * t85;
t115 = t88 * t87;
t37 = Icges(6,3) * t73 + (Icges(6,5) * t87 - Icges(6,6) * t85) * t74;
t39 = Icges(6,5) * t73 + (Icges(6,1) * t87 - Icges(6,4) * t85) * t74;
t114 = t74 * t87 * t39 + t73 * t37;
t40 = t73 * rSges(6,3) + (rSges(6,1) * t87 - rSges(6,2) * t85) * t74;
t113 = t74 * pkin(4) + t73 * pkin(7) + t40;
t51 = -t73 * t118 + t115;
t52 = t73 * t117 + t116;
t112 = t52 * rSges(6,1) + t51 * rSges(6,2);
t111 = t88 * pkin(1) + t86 * qJ(2);
t108 = Icges(6,5) * t74;
t107 = Icges(6,6) * t74;
t106 = Icges(6,3) * t74;
t105 = rSges(4,3) + qJ(3);
t104 = t62 + (m(4) + m(6)) * t69;
t103 = rSges(5,1) * t122 + rSges(5,2) * t121 + t88 * rSges(5,3);
t76 = t88 * qJ(2);
t84 = -pkin(6) - qJ(3);
t102 = t88 * t124 + t86 * t84 + t76;
t12 = -t37 * t121 + t51 * t38 + t52 * t39;
t21 = Icges(6,5) * t52 + Icges(6,6) * t51 - t86 * t106;
t23 = Icges(6,4) * t52 + Icges(6,2) * t51 - t86 * t107;
t25 = Icges(6,1) * t52 + Icges(6,4) * t51 - t86 * t108;
t9 = t73 * t21 + (-t23 * t85 + t25 * t87) * t74;
t101 = -t9 / 0.2e1 - t12 / 0.2e1;
t53 = t73 * t116 + t117;
t54 = -t73 * t115 + t118;
t22 = Icges(6,5) * t54 + Icges(6,6) * t53 + t88 * t106;
t24 = Icges(6,4) * t54 + Icges(6,2) * t53 + t88 * t107;
t26 = Icges(6,1) * t54 + Icges(6,4) * t53 + t88 * t108;
t10 = t73 * t22 + (-t24 * t85 + t26 * t87) * t74;
t13 = t37 * t120 + t53 * t38 + t54 * t39;
t100 = t10 / 0.2e1 + t13 / 0.2e1;
t99 = (-rSges(6,3) - pkin(7)) * t74;
t83 = cos(pkin(8));
t98 = rSges(4,1) * t82 + rSges(4,2) * t83;
t96 = -t54 * rSges(6,1) - t53 * rSges(6,2);
t91 = Icges(5,5) * t73 + Icges(5,6) * t74;
t90 = t86 * t124 - t88 * t84 + t111;
t31 = t128 + (-rSges(5,3) - pkin(1)) * t86 + t102;
t32 = t90 + t103;
t89 = m(5) * (t86 * t31 - t88 * t32);
t68 = pkin(4) * t122;
t65 = t88 * rSges(2,1) - t86 * rSges(2,2);
t64 = -t86 * rSges(2,1) - t88 * rSges(2,2);
t59 = t74 * rSges(5,1) - t73 * rSges(5,2);
t50 = -t88 * rSges(3,2) + t86 * rSges(3,3) + t111;
t49 = t88 * rSges(3,3) + t76 + (rSges(3,2) - pkin(1)) * t86;
t41 = Icges(5,3) * t88 + t91 * t86;
t35 = t105 * t88 + t98 * t86 + t111;
t34 = t76 + t98 * t88 + (-pkin(1) - t105) * t86;
t30 = t113 * t88;
t29 = t113 * t86;
t28 = rSges(6,3) * t120 - t96;
t27 = -rSges(6,3) * t121 + t112;
t20 = -t86 * t103 + (t86 * rSges(5,3) - t128) * t88;
t19 = t86 * t99 + t112 + t68 + t90;
t18 = -t86 * pkin(1) + (t99 + t123) * t88 + t96 + t102;
t17 = t40 * t120 - t73 * t28;
t16 = t40 * t121 + t73 * t27;
t15 = (-t74 * t119 + t114) * t73;
t14 = (-t27 * t88 - t28 * t86) * t74;
t11 = (pkin(7) * t121 - t27 - t68) * t86 + (t28 + (pkin(7) * t74 - t123) * t88) * t88;
t8 = t22 * t120 + t53 * t24 + t54 * t26;
t7 = t21 * t120 + t53 * t23 + t54 * t25;
t6 = -t22 * t121 + t51 * t24 + t52 * t26;
t5 = -t21 * t121 + t51 * t23 + t52 * t25;
t4 = t7 * t88 + t8 * t86;
t3 = t5 * t88 + t6 * t86;
t2 = t13 * t73 + (-t7 * t86 + t8 * t88) * t74;
t1 = t12 * t73 + (-t5 * t86 + t6 * t88) * t74;
t33 = [Icges(4,1) * t83 ^ 2 + Icges(3,1) + Icges(2,3) + (-0.2e1 * Icges(4,4) * t83 + Icges(4,2) * t82) * t82 + (Icges(5,1) * t74 - t119) * t74 + m(6) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t31 ^ 2 + t32 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2) + m(2) * (t64 ^ 2 + t65 ^ 2) + m(3) * (t49 ^ 2 + t50 ^ 2) + t114 + (-0.2e1 * Icges(5,4) * t74 + Icges(5,2) * t73) * t73; m(6) * (t86 * t18 - t88 * t19) + t89 + m(4) * (t86 * t34 - t88 * t35) + m(3) * (t86 * t49 - t88 * t50); m(3) * t69 + t104; m(6) * (t88 * t18 + t86 * t19) + m(5) * (t88 * t31 + t86 * t32) + m(4) * (t88 * t34 + t86 * t35); 0; t104; m(6) * (t29 * t18 - t30 * t19) + t59 * t89 + (t81 / 0.2e1 + t80 / 0.2e1) * (-t130 + t131) + (t129 * t88 - t101) * t88 + (t129 * t86 + t100) * t86; m(6) * (t29 * t86 + t30 * t88) + t59 * t62; m(6) * (t29 * t88 - t30 * t86); m(5) * (t69 * t59 ^ 2 + t20 ^ 2) + m(6) * (t11 ^ 2 + t29 ^ 2 + t30 ^ 2) + (t81 * t41 + t3) * t88 + (t86 * t41 * t88 + t4 + t69 * (Icges(5,3) * t86 - t91 * t88)) * t86; t15 + m(6) * (t16 * t19 + t17 * t18) + (t100 * t88 + t101 * t86) * t74; m(6) * (-t16 * t88 + t17 * t86); m(6) * (t16 * t86 + t17 * t88); m(6) * (t14 * t11 - t16 * t30 + t17 * t29) + t1 * t125 + t86 * t2 / 0.2e1 + t73 * (t10 * t86 + t9 * t88) / 0.2e1 + (-t86 * t3 / 0.2e1 + t4 * t125) * t74; m(6) * (t14 ^ 2 + t16 ^ 2 + t17 ^ 2) + t73 * t15 + (-t86 * t1 + t88 * t2 + t73 * (t10 * t88 - t86 * t9)) * t74;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t33(1), t33(2), t33(4), t33(7), t33(11); t33(2), t33(3), t33(5), t33(8), t33(12); t33(4), t33(5), t33(6), t33(9), t33(13); t33(7), t33(8), t33(9), t33(10), t33(14); t33(11), t33(12), t33(13), t33(14), t33(15);];
Mq = res;
