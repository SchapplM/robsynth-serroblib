% Calculate joint inertia matrix for
% S4RRPP5
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP5_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP5_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:10
% EndTime: 2019-12-31 17:00:12
% DurationCPUTime: 0.80s
% Computational Cost: add. (467->132), mult. (1076->190), div. (0->0), fcn. (944->4), ass. (0->73)
t143 = Icges(3,1) + Icges(5,3);
t142 = Icges(3,5) + Icges(5,5);
t141 = Icges(4,5) + Icges(5,4);
t140 = Icges(5,2) + Icges(4,3);
t71 = cos(qJ(2));
t139 = (Icges(4,6) - Icges(5,6)) * t71;
t69 = sin(qJ(2));
t138 = (Icges(3,4) - Icges(5,6)) * t69;
t137 = t143 * t71 - t138;
t136 = t140 * t69 - t139;
t135 = (-Icges(4,4) + t142) * t71 + (-Icges(3,6) + t141) * t69;
t134 = Icges(4,1) + Icges(5,1) + Icges(3,3);
t122 = Icges(3,6) / 0.2e1;
t133 = t122 - Icges(4,5) / 0.2e1 - Icges(5,4) / 0.2e1;
t127 = -Icges(4,4) / 0.2e1;
t132 = t127 + Icges(3,5) / 0.2e1 + Icges(5,5) / 0.2e1;
t70 = sin(qJ(1));
t131 = -t70 / 0.2e1;
t130 = t70 / 0.2e1;
t72 = cos(qJ(1));
t129 = -t72 / 0.2e1;
t128 = t72 / 0.2e1;
t121 = rSges(5,1) + pkin(3);
t97 = rSges(5,3) + qJ(4);
t66 = t70 ^ 2;
t68 = t72 ^ 2;
t105 = t66 + t68;
t118 = t121 * t72;
t117 = t134 * t70 + t135 * t72;
t116 = t134 * t72 - t135 * t70;
t115 = m(4) / 0.2e1;
t113 = t69 * t72;
t112 = t71 * t72;
t111 = t72 * rSges(4,1);
t110 = t72 * rSges(3,3);
t98 = qJ(3) * t69;
t107 = pkin(2) * t112 + t72 * t98;
t109 = t66 * (pkin(2) * t71 + t98) + t72 * t107;
t45 = t69 * pkin(2) - t71 * qJ(3);
t108 = t69 * rSges(4,2) + t71 * rSges(4,3) - t45;
t106 = t72 * pkin(1) + t70 * pkin(5);
t103 = Icges(3,4) * t71;
t102 = Icges(4,6) * t69;
t96 = t106 + t107;
t95 = t71 * rSges(5,2) - t97 * t69 - t45;
t93 = rSges(5,2) * t113 + t97 * t112 + t121 * t70;
t8 = t95 * t70;
t9 = t95 * t72;
t92 = t70 * t8 + t72 * t9;
t91 = rSges(3,1) * t71 - rSges(3,2) * t69;
t83 = -Icges(3,2) * t69 + t103;
t79 = -Icges(4,2) * t71 + t102;
t75 = rSges(3,1) * t112 - rSges(3,2) * t113 + t70 * rSges(3,3);
t74 = t70 * rSges(4,1) - rSges(4,2) * t112 + rSges(4,3) * t113;
t63 = t72 * pkin(5);
t3 = t63 + t118 + (-pkin(1) + (-rSges(5,2) - qJ(3)) * t69 + (-pkin(2) - t97) * t71) * t70;
t4 = t93 + t96;
t73 = m(5) * (t3 * t72 + t4 * t70);
t67 = t71 ^ 2;
t65 = t69 ^ 2;
t50 = t72 * rSges(2,1) - t70 * rSges(2,2);
t48 = -t70 * rSges(2,1) - t72 * rSges(2,2);
t47 = t69 * rSges(3,1) + t71 * rSges(3,2);
t13 = t108 * t72;
t12 = t108 * t70;
t11 = t75 + t106;
t10 = t110 + t63 + (-pkin(1) - t91) * t70;
t7 = t74 + t96;
t6 = t111 + t63 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t71 + (-rSges(4,3) - qJ(3)) * t69) * t70;
t5 = t72 * t75 + (t70 * t91 - t110) * t70;
t2 = t72 * t74 + (-t111 + (-rSges(4,2) * t71 + rSges(4,3) * t69) * t70) * t70 + t109;
t1 = t93 * t72 + ((rSges(5,2) * t69 + t71 * t97) * t70 - t118) * t70 + t109;
t14 = [Icges(2,3) + m(2) * (t48 ^ 2 + t50 ^ 2) + m(3) * (t10 ^ 2 + t11 ^ 2) + m(4) * (t6 ^ 2 + t7 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + (t102 + (Icges(3,2) + t140) * t71 + t138) * t71 + (t103 + (Icges(4,2) + t143) * t69 + t139) * t69; m(4) * (t12 * t7 + t13 * t6) + m(5) * (t9 * t3 + t8 * t4) + m(3) * (-t10 * t72 - t11 * t70) * t47 + ((t136 * t130 + t83 * t131 + t133 * t72) * t72 + (t83 * t128 + t136 * t129 + t133 * t70) * t70) * t71 + ((t79 * t130 + t137 * t131 + t132 * t72) * t72 + (t137 * t128 + t79 * t129 + t132 * t70) * t70) * t69 + t105 * ((t122 - t141 / 0.2e1) * t71 + (t127 + t142 / 0.2e1) * t69); m(4) * (t12 ^ 2 + t13 ^ 2 + t2 ^ 2) + m(5) * (t1 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(3) * (t105 * t47 ^ 2 + t5 ^ 2) + t116 * t72 * t68 + (t117 * t66 + (t116 * t70 + t117 * t72) * t72) * t70; 0.2e1 * ((t6 * t72 + t7 * t70) * t115 + t73 / 0.2e1) * t69; m(4) * (-t71 * t2 + (t12 * t70 + t13 * t72) * t69) + m(5) * (-t71 * t1 + t92 * t69); 0.2e1 * (t115 + m(5) / 0.2e1) * (t105 * t65 + t67); t71 * t73; m(5) * (t69 * t1 + t71 * t92); m(5) * (-0.1e1 + t105) * t71 * t69; m(5) * (t105 * t67 + t65);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t14(1), t14(2), t14(4), t14(7); t14(2), t14(3), t14(5), t14(8); t14(4), t14(5), t14(6), t14(9); t14(7), t14(8), t14(9), t14(10);];
Mq = res;
