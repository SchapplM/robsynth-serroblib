% Calculate joint inertia matrix for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:39
% EndTime: 2019-12-31 18:08:41
% DurationCPUTime: 0.95s
% Computational Cost: add. (1312->128), mult. (1077->171), div. (0->0), fcn. (933->8), ass. (0->65)
t73 = qJ(3) + pkin(8);
t70 = cos(t73);
t139 = t70 ^ 2;
t74 = qJ(1) + pkin(7);
t69 = sin(t74);
t66 = t69 ^ 2;
t71 = cos(t74);
t67 = t71 ^ 2;
t110 = t66 + t67;
t137 = Icges(6,4) + Icges(5,5);
t136 = Icges(5,6) - Icges(6,6);
t68 = sin(t73);
t76 = sin(qJ(3));
t78 = cos(qJ(3));
t132 = Icges(4,5) * t78 - Icges(4,6) * t76 - t136 * t68 + t137 * t70;
t131 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t127 = t69 * pkin(6);
t124 = rSges(6,1) + pkin(4);
t126 = t124 * t70;
t125 = rSges(6,3) + qJ(5);
t123 = t131 * t71 - t132 * t69;
t122 = t131 * t69 + t132 * t71;
t119 = pkin(3) * t76;
t77 = sin(qJ(1));
t118 = t77 * pkin(1);
t117 = rSges(4,1) * t78;
t116 = rSges(4,2) * t76;
t115 = t68 * t71;
t114 = t70 * t71;
t113 = t71 * rSges(4,3);
t64 = t78 * pkin(3) + pkin(2);
t51 = t71 * t64;
t63 = t71 * pkin(6);
t112 = t69 * (t63 + (-pkin(2) + t64) * t69) + t71 * (-t71 * pkin(2) - t127 + t51);
t111 = t69 * rSges(4,3) + t71 * t117;
t103 = qJ(5) * t68;
t102 = -t68 * rSges(5,1) - t70 * rSges(5,2) - t119;
t101 = -t124 * t68 + t125 * t70 - t119;
t79 = cos(qJ(1));
t72 = t79 * pkin(1);
t75 = -qJ(4) - pkin(6);
t99 = -t69 * t75 + t51 + t72;
t98 = t69 * rSges(6,2) + rSges(6,3) * t115 + t71 * t103 + t124 * t114;
t97 = -t116 + t117;
t96 = rSges(5,1) * t70 - rSges(5,2) * t68;
t80 = rSges(5,1) * t114 - rSges(5,2) * t115 + t69 * rSges(5,3);
t58 = t79 * rSges(2,1) - t77 * rSges(2,2);
t57 = -t77 * rSges(2,1) - t79 * rSges(2,2);
t56 = t76 * rSges(4,1) + t78 * rSges(4,2);
t35 = t71 * rSges(3,1) - t69 * rSges(3,2) + t72;
t34 = -t69 * rSges(3,1) - t71 * rSges(3,2) - t118;
t27 = t102 * t71;
t26 = t102 * t69;
t11 = t127 + t72 + (pkin(2) - t116) * t71 + t111;
t10 = t113 - t118 + t63 + (-pkin(2) - t97) * t69;
t9 = t101 * t71;
t8 = t101 * t69;
t7 = t80 + t99;
t6 = -t118 + (rSges(5,3) - t75) * t71 + (-t64 - t96) * t69;
t5 = t71 * (-t71 * t116 + t111) + (t97 * t69 - t113) * t69;
t4 = t98 + t99;
t3 = -t118 + (rSges(6,2) - t75) * t71 + (-t125 * t68 - t126 - t64) * t69;
t2 = t71 * t80 + (-t71 * rSges(5,3) + t96 * t69) * t69 + t112;
t1 = t98 * t71 + (-t71 * rSges(6,2) + (rSges(6,3) * t68 + t103 + t126) * t69) * t69 + t112;
t12 = [Icges(4,2) * t78 ^ 2 + Icges(2,3) + Icges(3,3) + m(2) * (t57 ^ 2 + t58 ^ 2) + m(3) * (t34 ^ 2 + t35 ^ 2) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2) + (Icges(5,2) + Icges(6,3)) * t139 + (Icges(4,1) * t76 + 0.2e1 * Icges(4,4) * t78) * t76 + (0.2e1 * (-Icges(6,5) + Icges(5,4)) * t70 + (Icges(5,1) + Icges(6,1)) * t68) * t68; 0; m(3) + m(4) + m(5) + m(6); m(5) * (t26 * t7 + t27 * t6) + m(6) * (t9 * t3 + t8 * t4) + m(4) * (-t10 * t71 - t11 * t69) * t56 + t110 * (Icges(4,5) * t76 + Icges(4,6) * t78 + t136 * t70 + t137 * t68); m(4) * t5 + m(5) * t2 + m(6) * t1; m(6) * (t1 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t2 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(4) * (t110 * t56 ^ 2 + t5 ^ 2) + t122 * t69 * t66 + (t123 * t67 + (t122 * t71 + t123 * t69) * t69) * t71; m(5) * (t69 * t6 - t71 * t7) + m(6) * (t69 * t3 - t71 * t4); 0; m(6) * (t69 * t9 - t71 * t8) + m(5) * (-t71 * t26 + t69 * t27); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t110; m(6) * (t3 * t71 + t4 * t69) * t68; -m(6) * t70; m(6) * (-t70 * t1 + (t69 * t8 + t71 * t9) * t68); 0; m(6) * (t110 * t68 ^ 2 + t139);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t12(1), t12(2), t12(4), t12(7), t12(11); t12(2), t12(3), t12(5), t12(8), t12(12); t12(4), t12(5), t12(6), t12(9), t12(13); t12(7), t12(8), t12(9), t12(10), t12(14); t12(11), t12(12), t12(13), t12(14), t12(15);];
Mq = res;
