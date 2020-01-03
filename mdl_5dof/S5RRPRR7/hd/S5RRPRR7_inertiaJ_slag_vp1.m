% Calculate joint inertia matrix for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:28
% EndTime: 2019-12-31 20:15:30
% DurationCPUTime: 0.64s
% Computational Cost: add. (1835->163), mult. (1534->236), div. (0->0), fcn. (1308->8), ass. (0->91)
t90 = qJ(1) + qJ(2);
t85 = sin(t90);
t87 = cos(t90);
t137 = t85 * t87;
t89 = qJ(4) + qJ(5);
t84 = sin(t89);
t86 = cos(t89);
t138 = rSges(6,1) * t84 + rSges(6,2) * t86;
t91 = sin(qJ(4));
t128 = t85 * t91;
t31 = t87 * rSges(6,3) + t138 * t85;
t136 = -pkin(4) * t128 - t31;
t82 = t85 ^ 2;
t83 = t87 ^ 2;
t135 = t85 / 0.2e1;
t134 = t87 / 0.2e1;
t102 = Icges(6,5) * t84 + Icges(6,6) * t86;
t25 = Icges(6,3) * t87 + t102 * t85;
t26 = Icges(6,3) * t85 - t102 * t87;
t133 = t85 * (t25 * t137 + t82 * t26) + t87 * (t26 * t137 + t83 * t25);
t92 = sin(qJ(1));
t132 = t92 * pkin(1);
t93 = cos(qJ(4));
t130 = rSges(5,2) * t93;
t127 = t87 * t91;
t126 = t138 * t87;
t125 = rSges(5,1) * t127 + t87 * t130;
t95 = -pkin(8) - pkin(7);
t124 = -pkin(4) * t127 - t85 * t95;
t123 = t87 * pkin(2) + t85 * qJ(3);
t122 = t82 + t83;
t121 = Icges(5,4) * t91;
t120 = Icges(5,4) * t93;
t119 = Icges(6,4) * t84;
t118 = Icges(6,4) * t86;
t117 = rSges(5,1) * t128 + t87 * rSges(5,3) + t85 * t130;
t104 = Icges(6,2) * t86 + t119;
t106 = Icges(6,1) * t84 + t118;
t47 = -Icges(6,2) * t84 + t118;
t48 = Icges(6,1) * t86 - t119;
t109 = t47 * t86 + t48 * t84;
t46 = Icges(6,5) * t86 - Icges(6,6) * t84;
t116 = (-t109 * t87 - t84 * (Icges(6,6) * t85 - t104 * t87) + t86 * (Icges(6,5) * t85 - t106 * t87) + t85 * t46) * t135 + (t109 * t85 - t84 * (Icges(6,6) * t87 + t104 * t85) + t86 * (Icges(6,5) * t87 + t106 * t85) + t87 * t46) * t134;
t50 = t86 * rSges(6,1) - t84 * rSges(6,2);
t115 = pkin(4) * t93 + t50;
t51 = t87 * rSges(3,1) - t85 * rSges(3,2);
t49 = -t85 * rSges(3,1) - t87 * rSges(3,2);
t32 = t115 * t85;
t33 = t115 * t87;
t112 = t32 * t85 + t33 * t87;
t60 = -Icges(5,2) * t91 + t120;
t61 = Icges(5,1) * t93 - t121;
t108 = t60 * t93 + t61 * t91;
t73 = t87 * qJ(3);
t34 = t87 * rSges(4,3) + t73 + (rSges(4,2) - pkin(2)) * t85;
t107 = Icges(5,1) * t91 + t120;
t105 = Icges(5,2) * t93 + t121;
t103 = Icges(5,5) * t91 + Icges(5,6) * t93;
t35 = -t87 * rSges(4,2) + t85 * rSges(4,3) + t123;
t80 = t87 * pkin(7);
t21 = t80 + t117 + t123;
t20 = t73 + (-rSges(5,3) - pkin(2) - pkin(7)) * t85 + t125;
t18 = t20 - t132;
t94 = cos(qJ(1));
t88 = t94 * pkin(1);
t19 = t88 + t21;
t101 = m(5) * (t85 * t18 - t87 * t19);
t100 = m(5) * (t85 * t20 - t87 * t21);
t14 = t73 + (-rSges(6,3) - pkin(2)) * t85 - t124 + t126;
t11 = t14 - t132;
t15 = -t87 * t95 + t123 - t136;
t12 = t15 + t88;
t99 = m(6) * (t85 * t11 - t87 * t12);
t98 = m(6) * (t85 * t14 - t87 * t15);
t59 = Icges(5,5) * t93 - Icges(5,6) * t91;
t97 = t116 + (-t108 * t87 - t91 * (Icges(5,6) * t85 - t105 * t87) + t93 * (Icges(5,5) * t85 - t107 * t87) + t85 * t59) * t135 + (t108 * t85 - t91 * (Icges(5,6) * t87 + t105 * t85) + t93 * (Icges(5,5) * t87 + t107 * t85) + t87 * t59) * t134;
t96 = -t84 * t47 + t86 * t48 - t91 * t60 + t93 * t61 + Icges(4,1) + Icges(3,3);
t64 = t94 * rSges(2,1) - t92 * rSges(2,2);
t63 = t93 * rSges(5,1) - t91 * rSges(5,2);
t62 = -t92 * rSges(2,1) - t94 * rSges(2,2);
t44 = t51 + t88;
t43 = t49 - t132;
t37 = Icges(5,3) * t85 - t103 * t87;
t36 = Icges(5,3) * t87 + t103 * t85;
t24 = t35 + t88;
t23 = t34 - t132;
t22 = t87 * (t85 * rSges(6,3) - t126);
t13 = t87 * (t85 * rSges(5,3) - t125) - t85 * t117;
t8 = -t85 * t31 + t22;
t3 = t87 * t124 + t22 + (t80 + (-pkin(7) + t95) * t87 + t136) * t85;
t1 = [Icges(2,3) + m(6) * (t11 ^ 2 + t12 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(4) * (t23 ^ 2 + t24 ^ 2) + m(3) * (t43 ^ 2 + t44 ^ 2) + m(2) * (t62 ^ 2 + t64 ^ 2) + t96; m(6) * (t14 * t11 + t15 * t12) + m(5) * (t20 * t18 + t21 * t19) + m(4) * (t34 * t23 + t35 * t24) + m(3) * (t49 * t43 + t51 * t44) + t96; m(6) * (t14 ^ 2 + t15 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2) + m(3) * (t49 ^ 2 + t51 ^ 2) + t96; t99 + t101 + m(4) * (t85 * t23 - t87 * t24); t98 + t100 + m(4) * (t85 * t34 - t87 * t35); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t122; m(6) * (t32 * t11 - t33 * t12) + t63 * t101 + t97; m(6) * (t32 * t14 - t33 * t15) + t63 * t100 + t97; m(5) * t122 * t63 + m(6) * t112; m(5) * (t122 * t63 ^ 2 + t13 ^ 2) + t87 * (t37 * t137 + t83 * t36) + t85 * (t36 * t137 + t82 * t37) + m(6) * (t3 ^ 2 + t32 ^ 2 + t33 ^ 2) + t133; t50 * t99 + t116; t50 * t98 + t116; m(6) * t122 * t50; m(6) * (t112 * t50 + t8 * t3) + t133; m(6) * (t122 * t50 ^ 2 + t8 ^ 2) + t133;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
