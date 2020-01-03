% Calculate joint inertia matrix for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:54
% EndTime: 2020-01-03 11:58:56
% DurationCPUTime: 0.64s
% Computational Cost: add. (1446->117), mult. (1014->153), div. (0->0), fcn. (852->8), ass. (0->66)
t133 = Icges(5,5) + Icges(6,5);
t132 = -Icges(5,6) - Icges(6,6);
t90 = cos(qJ(4));
t128 = t132 * t90;
t88 = sin(qJ(4));
t129 = t133 * t88;
t127 = t128 - t129;
t86 = qJ(1) + qJ(2);
t81 = pkin(8) + t86;
t76 = sin(t81);
t71 = t76 ^ 2;
t77 = cos(t81);
t72 = t77 ^ 2;
t131 = t132 * t88 + t133 * t90;
t130 = Icges(5,3) + Icges(6,3);
t125 = t130 * t77 - t131 * t76;
t124 = t130 * t76 + t131 * t77;
t60 = rSges(5,1) * t88 + rSges(5,2) * t90;
t121 = m(5) * t60;
t120 = t76 * t88;
t119 = t76 * t90;
t118 = t77 * t88;
t117 = t77 * t90;
t116 = pkin(3) * t77 + pkin(7) * t76;
t115 = t72 + t71;
t82 = sin(t86);
t83 = cos(t86);
t39 = rSges(3,1) * t82 + rSges(3,2) * t83;
t78 = pkin(2) * t82;
t35 = rSges(4,1) * t76 + rSges(4,2) * t77 + t78;
t110 = t90 * rSges(6,2) + (rSges(6,1) + pkin(4)) * t88;
t40 = rSges(3,1) * t83 - rSges(3,2) * t82;
t109 = rSges(5,1) * t119 - rSges(5,2) * t120;
t79 = pkin(2) * t83;
t36 = rSges(4,1) * t77 - rSges(4,2) * t76 + t79;
t96 = -rSges(5,1) * t117 + rSges(5,2) * t118 - rSges(5,3) * t76;
t95 = Icges(3,3) + Icges(4,3) + (Icges(5,2) + Icges(6,2)) * t90 ^ 2 + ((Icges(5,1) + Icges(6,1)) * t88 + (2 * Icges(5,4) + 2 * Icges(6,4)) * t90) * t88;
t80 = pkin(4) * t90 + pkin(3);
t87 = -qJ(5) - pkin(7);
t94 = rSges(6,1) * t119 - rSges(6,2) * t120 + t76 * t80 + t77 * t87;
t93 = rSges(6,1) * t117 - rSges(6,2) * t118 + rSges(6,3) * t76 + t77 * t80;
t92 = -t127 * t72 + (-t128 / 0.2e1 + t129 / 0.2e1 - t127 / 0.2e1) * t71;
t18 = t79 - t96 + t116;
t13 = -rSges(6,3) * t77 + t78 + t94;
t14 = -t76 * t87 + t79 + t93;
t69 = t76 * pkin(3);
t17 = t69 + t78 + (-rSges(5,3) - pkin(7)) * t77 + t109;
t91 = cos(qJ(1));
t89 = sin(qJ(1));
t85 = t91 * pkin(1);
t84 = t89 * pkin(1);
t62 = rSges(2,1) * t91 - rSges(2,2) * t89;
t61 = rSges(2,1) * t89 + rSges(2,2) * t91;
t38 = t40 + t85;
t37 = t84 + t39;
t34 = t36 + t85;
t33 = t84 + t35;
t32 = t110 * t77;
t31 = t110 * t76;
t16 = t18 + t85;
t15 = t84 + t17;
t12 = t14 + t85;
t11 = t13 + t84;
t6 = -t77 * t96 + t76 * (-rSges(5,3) * t77 + t109);
t1 = (t93 - t116) * t77 + (-t69 + (-rSges(6,3) + pkin(7) - t87) * t77 + t94) * t76;
t2 = [Icges(2,3) + m(2) * (t61 ^ 2 + t62 ^ 2) + m(3) * (t37 ^ 2 + t38 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(6) * (t11 ^ 2 + t12 ^ 2) + t95; m(3) * (t37 * t39 + t38 * t40) + m(4) * (t33 * t35 + t34 * t36) + m(5) * (t15 * t17 + t16 * t18) + m(6) * (t11 * t13 + t12 * t14) + t95; m(6) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2) + m(3) * (t39 ^ 2 + t40 ^ 2) + t95; 0; 0; m(4) + m(5) + m(6); m(6) * (t11 * t32 - t12 * t31) + (t15 * t77 - t16 * t76) * t121 + t92; m(6) * (t13 * t32 - t14 * t31) + (t17 * t77 - t18 * t76) * t121 + t92; m(5) * t6 + m(6) * t1; m(5) * (t115 * t60 ^ 2 + t6 ^ 2) + m(6) * (t1 ^ 2 + t31 ^ 2 + t32 ^ 2) + t125 * t77 * t72 + (t124 * t71 + (t124 * t77 + t125 * t76) * t77) * t76; m(6) * (-t11 * t76 - t12 * t77); m(6) * (-t13 * t76 - t14 * t77); 0; m(6) * (t31 * t77 - t32 * t76); m(6) * t115;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
