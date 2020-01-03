% Calculate joint inertia matrix for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:58:55
% DurationCPUTime: 0.83s
% Computational Cost: add. (450->122), mult. (1022->177), div. (0->0), fcn. (898->4), ass. (0->63)
t137 = -Icges(5,4) - Icges(4,5);
t136 = Icges(3,5) + Icges(4,4);
t135 = Icges(3,6) + Icges(5,6);
t134 = Icges(3,1) + Icges(4,1) + Icges(5,1);
t67 = sin(qJ(2));
t133 = (Icges(3,4) + t137) * t67;
t69 = cos(qJ(2));
t131 = (Icges(5,5) - t136) * t69 + (-Icges(4,6) + t135) * t67;
t130 = t134 * t69 - t133;
t129 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t118 = -Icges(4,6) / 0.2e1;
t128 = t118 + Icges(3,6) / 0.2e1 + Icges(5,6) / 0.2e1;
t120 = -Icges(5,5) / 0.2e1;
t127 = t120 + Icges(3,5) / 0.2e1 + Icges(4,4) / 0.2e1;
t68 = sin(qJ(1));
t126 = -t68 / 0.2e1;
t70 = cos(qJ(1));
t123 = t70 / 0.2e1;
t110 = rSges(5,1) + pkin(3);
t65 = t68 ^ 2;
t66 = t70 ^ 2;
t101 = t65 + t66;
t114 = t129 * t70 + t131 * t68;
t113 = t129 * t68 - t131 * t70;
t112 = m(4) / 0.2e1;
t111 = m(5) / 0.2e1;
t109 = t67 * t70;
t108 = t69 * t70;
t107 = t70 * rSges(4,2);
t106 = t70 * rSges(3,3);
t94 = qJ(3) * t67;
t103 = pkin(2) * t108 + t70 * t94;
t105 = t65 * (pkin(2) * t69 + t94) + t70 * t103;
t45 = t67 * pkin(2) - t69 * qJ(3);
t104 = -t67 * rSges(4,1) + t69 * rSges(4,3) - t45;
t102 = t70 * pkin(1) + t68 * pkin(5);
t99 = Icges(3,4) * t69;
t93 = -rSges(5,3) - qJ(4);
t92 = rSges(5,2) * t109 + t110 * t108;
t91 = rSges(4,1) * t108 + t68 * rSges(4,2) + rSges(4,3) * t109;
t90 = t102 + t103;
t89 = t69 * rSges(5,2) - t110 * t67 - t45;
t87 = rSges(3,1) * t69 - rSges(3,2) * t67;
t77 = -Icges(3,2) * t67 + t99;
t71 = rSges(3,1) * t108 - rSges(3,2) * t109 + t68 * rSges(3,3);
t62 = t70 * pkin(5);
t50 = t70 * rSges(2,1) - t68 * rSges(2,2);
t49 = -t68 * rSges(2,1) - t70 * rSges(2,2);
t48 = t67 * rSges(3,1) + t69 * rSges(3,2);
t13 = t104 * t70;
t12 = t104 * t68;
t11 = t71 + t102;
t10 = t106 + t62 + (-pkin(1) - t87) * t68;
t9 = t89 * t70;
t8 = t89 * t68;
t7 = t90 + t91;
t6 = t107 + t62 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t69 + (-rSges(4,3) - qJ(3)) * t67) * t68;
t5 = t70 * t71 + (t68 * t87 - t106) * t68;
t4 = t68 * t93 + t90 + t92;
t3 = t62 + t93 * t70 + (-pkin(1) + (-rSges(5,2) - qJ(3)) * t67 + (-pkin(2) - t110) * t69) * t68;
t2 = t70 * t91 + (-t107 + (rSges(4,1) * t69 + rSges(4,3) * t67) * t68) * t68 + t105;
t1 = t92 * t70 + (rSges(5,2) * t67 + t110 * t69) * t65 + t105;
t14 = [Icges(2,3) + m(2) * (t49 ^ 2 + t50 ^ 2) + m(3) * (t10 ^ 2 + t11 ^ 2) + m(4) * (t6 ^ 2 + t7 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + (t134 * t67 + t99) * t67 + (t133 + t137 * t67 + (Icges(3,2) + Icges(5,2) + Icges(4,3)) * t69) * t69; m(4) * (t12 * t7 + t13 * t6) + m(5) * (t9 * t3 + t8 * t4) + m(3) * (-t10 * t70 - t11 * t68) * t48 + ((t77 * t126 + t128 * t70) * t70 + (t77 * t123 + t128 * t68) * t68) * t69 + ((t130 * t126 + t127 * t70) * t70 + (t130 * t123 + t127 * t68) * t68) * t67 + t101 * ((t118 + t135 / 0.2e1) * t69 + (t120 + t136 / 0.2e1) * t67); m(4) * (t12 ^ 2 + t13 ^ 2 + t2 ^ 2) + m(5) * (t1 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(3) * (t101 * t48 ^ 2 + t5 ^ 2) + t113 * t68 * t65 + (t114 * t66 + (t113 * t70 + t114 * t68) * t68) * t70; 0.2e1 * ((t6 * t70 + t68 * t7) * t112 + (t3 * t70 + t4 * t68) * t111) * t67; m(4) * (-t69 * t2 + (t12 * t68 + t13 * t70) * t67) + m(5) * (-t69 * t1 + (t68 * t8 + t70 * t9) * t67); 0.2e1 * (t112 + t111) * (t101 * t67 ^ 2 + t69 ^ 2); m(5) * (-t68 * t3 + t70 * t4); m(5) * (-t68 * t9 + t70 * t8); 0; m(5) * t101;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t14(1), t14(2), t14(4), t14(7); t14(2), t14(3), t14(5), t14(8); t14(4), t14(5), t14(6), t14(9); t14(7), t14(8), t14(9), t14(10);];
Mq = res;
