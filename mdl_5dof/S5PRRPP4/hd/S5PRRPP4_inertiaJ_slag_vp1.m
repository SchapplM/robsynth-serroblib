% Calculate joint inertia matrix for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:47
% EndTime: 2019-12-31 17:40:49
% DurationCPUTime: 0.81s
% Computational Cost: add. (1026->129), mult. (1106->181), div. (0->0), fcn. (964->4), ass. (0->65)
t139 = -Icges(6,4) - Icges(5,5);
t138 = Icges(4,5) + Icges(5,4);
t137 = Icges(4,6) + Icges(6,6);
t136 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t70 = sin(qJ(3));
t135 = (Icges(4,4) + t139) * t70;
t71 = cos(qJ(3));
t133 = (Icges(6,5) - t138) * t71 + (-Icges(5,6) + t137) * t70;
t132 = t136 * t71 - t135;
t131 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t120 = -Icges(5,6) / 0.2e1;
t130 = t120 + Icges(4,6) / 0.2e1 + Icges(6,6) / 0.2e1;
t122 = -Icges(6,5) / 0.2e1;
t129 = t122 + Icges(4,5) / 0.2e1 + Icges(5,4) / 0.2e1;
t68 = pkin(7) + qJ(2);
t66 = sin(t68);
t128 = -t66 / 0.2e1;
t67 = cos(t68);
t125 = t67 / 0.2e1;
t111 = rSges(6,1) + pkin(4);
t64 = t66 ^ 2;
t65 = t67 ^ 2;
t102 = t64 + t65;
t116 = t131 * t67 + t133 * t66;
t115 = t131 * t66 - t133 * t67;
t114 = m(5) / 0.2e1;
t113 = m(6) / 0.2e1;
t112 = -m(5) - m(6);
t110 = t67 * rSges(5,2);
t109 = t67 * rSges(4,3);
t108 = t67 * t70;
t107 = t67 * t71;
t95 = qJ(4) * t70;
t104 = pkin(3) * t107 + t67 * t95;
t106 = t64 * (pkin(3) * t71 + t95) + t67 * t104;
t47 = t70 * pkin(3) - t71 * qJ(4);
t105 = -t70 * rSges(5,1) + t71 * rSges(5,3) - t47;
t103 = t67 * pkin(2) + t66 * pkin(6);
t100 = Icges(4,4) * t71;
t94 = -rSges(6,3) - qJ(5);
t93 = rSges(6,2) * t108 + t111 * t107;
t92 = rSges(5,1) * t107 + t66 * rSges(5,2) + rSges(5,3) * t108;
t91 = t103 + t104;
t90 = t71 * rSges(6,2) - t111 * t70 - t47;
t88 = rSges(4,1) * t71 - rSges(4,2) * t70;
t78 = -Icges(4,2) * t70 + t100;
t72 = rSges(4,1) * t107 - rSges(4,2) * t108 + t66 * rSges(4,3);
t62 = t67 * pkin(6);
t50 = t70 * rSges(4,1) + t71 * rSges(4,2);
t36 = t67 * rSges(3,1) - t66 * rSges(3,2);
t35 = -t66 * rSges(3,1) - t67 * rSges(3,2);
t13 = t105 * t67;
t12 = t105 * t66;
t11 = t72 + t103;
t10 = t109 + t62 + (-pkin(2) - t88) * t66;
t9 = t90 * t67;
t8 = t90 * t66;
t7 = t91 + t92;
t6 = t110 + t62 + (-pkin(2) + (-rSges(5,1) - pkin(3)) * t71 + (-rSges(5,3) - qJ(4)) * t70) * t66;
t5 = t67 * t72 + (t88 * t66 - t109) * t66;
t4 = t94 * t66 + t91 + t93;
t3 = t62 + t94 * t67 + (-pkin(2) + (-rSges(6,2) - qJ(4)) * t70 + (-pkin(3) - t111) * t71) * t66;
t2 = t67 * t92 + (-t110 + (rSges(5,1) * t71 + rSges(5,3) * t70) * t66) * t66 + t106;
t1 = t93 * t67 + (rSges(6,2) * t70 + t111 * t71) * t64 + t106;
t14 = [m(2) + m(3) + m(4) - t112; 0; Icges(3,3) + m(5) * (t6 ^ 2 + t7 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(3) * (t35 ^ 2 + t36 ^ 2) + (t136 * t70 + t100) * t70 + (t135 + t139 * t70 + (Icges(4,2) + Icges(6,2) + Icges(5,3)) * t71) * t71; m(4) * t5 + m(5) * t2 + m(6) * t1; m(5) * (t12 * t7 + t13 * t6) + m(6) * (t9 * t3 + t8 * t4) + m(4) * (-t10 * t67 - t11 * t66) * t50 + ((t78 * t128 + t130 * t67) * t67 + (t78 * t125 + t130 * t66) * t66) * t71 + ((t132 * t128 + t129 * t67) * t67 + (t132 * t125 + t129 * t66) * t66) * t70 + t102 * ((t120 + t137 / 0.2e1) * t71 + (t122 + t138 / 0.2e1) * t70); m(6) * (t1 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2 + t2 ^ 2) + m(4) * (t102 * t50 ^ 2 + t5 ^ 2) + t115 * t66 * t64 + (t116 * t65 + (t115 * t67 + t116 * t66) * t66) * t67; t112 * t71; 0.2e1 * ((t6 * t67 + t66 * t7) * t114 + (t3 * t67 + t4 * t66) * t113) * t70; m(6) * (-t71 * t1 + (t66 * t8 + t67 * t9) * t70) + m(5) * (-t71 * t2 + (t12 * t66 + t13 * t67) * t70); 0.2e1 * (t114 + t113) * (t102 * t70 ^ 2 + t71 ^ 2); 0; m(6) * (-t66 * t3 + t67 * t4); m(6) * (-t66 * t9 + t67 * t8); 0; m(6) * t102;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t14(1), t14(2), t14(4), t14(7), t14(11); t14(2), t14(3), t14(5), t14(8), t14(12); t14(4), t14(5), t14(6), t14(9), t14(13); t14(7), t14(8), t14(9), t14(10), t14(14); t14(11), t14(12), t14(13), t14(14), t14(15);];
Mq = res;
