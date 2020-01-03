% Calculate joint inertia matrix for
% S5RRPRP2
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:26
% EndTime: 2019-12-31 19:49:28
% DurationCPUTime: 0.61s
% Computational Cost: add. (1526->116), mult. (1103->158), div. (0->0), fcn. (937->8), ass. (0->69)
t84 = cos(qJ(4));
t136 = t84 ^ 2;
t135 = Icges(6,4) + Icges(5,5);
t134 = Icges(5,6) - Icges(6,6);
t82 = sin(qJ(4));
t132 = t135 * t82;
t133 = t134 * t84;
t127 = t132 + t133;
t81 = qJ(1) + qJ(2);
t76 = pkin(8) + t81;
t73 = sin(t76);
t70 = t73 ^ 2;
t74 = cos(t76);
t71 = t74 ^ 2;
t131 = -t134 * t82 + t135 * t84;
t129 = Icges(6,2) + Icges(5,3);
t123 = rSges(6,1) + pkin(4);
t125 = t123 * t84;
t124 = rSges(6,3) + qJ(5);
t122 = t129 * t74 - t131 * t73;
t121 = t129 * t73 + t131 * t74;
t59 = t82 * rSges(5,1) + t84 * rSges(5,2);
t118 = m(5) * t59;
t117 = m(6) * t82;
t77 = sin(t81);
t116 = pkin(2) * t77;
t83 = sin(qJ(1));
t115 = t83 * pkin(1);
t114 = rSges(5,1) * t84;
t113 = t74 * t82;
t112 = t74 * t84;
t111 = t73 * t82 * rSges(5,2) + t74 * rSges(5,3);
t110 = -t123 * t82 + t124 * t84;
t109 = t70 + t71;
t104 = qJ(5) * t82;
t78 = cos(t81);
t75 = pkin(2) * t78;
t103 = t74 * pkin(3) + t73 * pkin(7) + t75;
t102 = t74 * pkin(7) - t116;
t40 = t78 * rSges(3,1) - t77 * rSges(3,2);
t36 = t74 * rSges(4,1) - t73 * rSges(4,2) + t75;
t101 = t73 * rSges(6,2) + rSges(6,3) * t113 + t74 * t104 + t123 * t112;
t39 = -t77 * rSges(3,1) - t78 * rSges(3,2);
t88 = rSges(5,1) * t112 - rSges(5,2) * t113 + t73 * rSges(5,3);
t87 = t127 * t71 + (t133 / 0.2e1 + t132 / 0.2e1 + t127 / 0.2e1) * t70;
t35 = -t73 * rSges(4,1) - t74 * rSges(4,2) - t116;
t10 = t101 + t103;
t86 = Icges(3,3) + Icges(4,3) + (Icges(5,2) + Icges(6,3)) * t136 + ((Icges(5,1) + Icges(6,1)) * t82 + (2 * Icges(5,4) - 2 * Icges(6,5)) * t84) * t82;
t18 = t88 + t103;
t17 = (-pkin(3) - t114) * t73 + t102 + t111;
t65 = t74 * rSges(6,2);
t9 = t65 + (-t124 * t82 - pkin(3) - t125) * t73 + t102;
t85 = cos(qJ(1));
t79 = t85 * pkin(1);
t61 = t85 * rSges(2,1) - t83 * rSges(2,2);
t60 = -t83 * rSges(2,1) - t85 * rSges(2,2);
t38 = t40 + t79;
t37 = t39 - t115;
t34 = t36 + t79;
t33 = t35 - t115;
t20 = t110 * t74;
t19 = t110 * t73;
t16 = t18 + t79;
t15 = t17 - t115;
t8 = t73 * (t73 * t114 - t111) + t74 * t88;
t3 = t79 + t10;
t2 = t9 - t115;
t1 = t101 * t74 + (-t65 + (rSges(6,3) * t82 + t104 + t125) * t73) * t73;
t4 = [Icges(2,3) + m(6) * (t2 ^ 2 + t3 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t37 ^ 2 + t38 ^ 2) + m(2) * (t60 ^ 2 + t61 ^ 2) + t86; m(6) * (t10 * t3 + t9 * t2) + m(5) * (t17 * t15 + t18 * t16) + m(4) * (t35 * t33 + t36 * t34) + m(3) * (t39 * t37 + t40 * t38) + t86; m(6) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2) + m(3) * (t39 ^ 2 + t40 ^ 2) + t86; 0; 0; m(4) + m(5) + m(6); m(6) * (t19 * t3 + t20 * t2) + (-t15 * t74 - t16 * t73) * t118 + t87; m(6) * (t19 * t10 + t20 * t9) + (-t17 * t74 - t18 * t73) * t118 + t87; m(5) * t8 + m(6) * t1; m(5) * (t109 * t59 ^ 2 + t8 ^ 2) + m(6) * (t1 ^ 2 + t19 ^ 2 + t20 ^ 2) + t122 * t74 * t71 + (t121 * t70 + (t121 * t74 + t122 * t73) * t74) * t73; (t2 * t74 + t3 * t73) * t117; (t10 * t73 + t74 * t9) * t117; -m(6) * t84; m(6) * (-t84 * t1 + (t19 * t73 + t20 * t74) * t82); m(6) * (t109 * t82 ^ 2 + t136);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
