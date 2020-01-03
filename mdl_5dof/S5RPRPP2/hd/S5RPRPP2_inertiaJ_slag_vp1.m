% Calculate joint inertia matrix for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:41
% EndTime: 2019-12-31 18:10:43
% DurationCPUTime: 0.84s
% Computational Cost: add. (1050->137), mult. (1136->190), div. (0->0), fcn. (988->6), ass. (0->71)
t146 = -Icges(6,4) - Icges(5,5);
t145 = Icges(4,5) + Icges(5,4);
t144 = Icges(4,6) + Icges(6,6);
t143 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t73 = sin(qJ(3));
t142 = (Icges(4,4) + t146) * t73;
t75 = cos(qJ(3));
t140 = (Icges(6,5) - t145) * t75 + (-Icges(5,6) + t144) * t73;
t139 = t143 * t75 - t142;
t138 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t127 = -Icges(5,6) / 0.2e1;
t137 = t127 + Icges(4,6) / 0.2e1 + Icges(6,6) / 0.2e1;
t129 = -Icges(6,5) / 0.2e1;
t136 = t129 + Icges(4,5) / 0.2e1 + Icges(5,4) / 0.2e1;
t71 = qJ(1) + pkin(7);
t68 = sin(t71);
t135 = -t68 / 0.2e1;
t69 = cos(t71);
t132 = t69 / 0.2e1;
t117 = rSges(6,1) + pkin(4);
t66 = t68 ^ 2;
t67 = t69 ^ 2;
t109 = t66 + t67;
t123 = t138 * t69 + t140 * t68;
t122 = t138 * t68 - t140 * t69;
t121 = m(5) / 0.2e1;
t120 = m(6) / 0.2e1;
t119 = -m(5) - m(6);
t74 = sin(qJ(1));
t118 = t74 * pkin(1);
t116 = t69 * rSges(5,2);
t115 = t69 * rSges(4,3);
t114 = t69 * t73;
t113 = t69 * t75;
t102 = qJ(4) * t73;
t110 = pkin(3) * t113 + t69 * t102;
t112 = t66 * (pkin(3) * t75 + t102) + t69 * t110;
t47 = t73 * pkin(3) - t75 * qJ(4);
t111 = -t73 * rSges(5,1) + t75 * rSges(5,3) - t47;
t107 = Icges(4,4) * t75;
t101 = -rSges(6,3) - qJ(5);
t100 = rSges(6,2) * t114 + t117 * t113;
t99 = rSges(5,1) * t113 + t68 * rSges(5,2) + rSges(5,3) * t114;
t76 = cos(qJ(1));
t70 = t76 * pkin(1);
t98 = t69 * pkin(2) + t68 * pkin(6) + t70;
t97 = t69 * pkin(6) - t118;
t96 = t75 * rSges(6,2) - t117 * t73 - t47;
t94 = t98 + t110;
t93 = rSges(4,1) * t75 - rSges(4,2) * t73;
t83 = -Icges(4,2) * t73 + t107;
t77 = rSges(4,1) * t113 - rSges(4,2) * t114 + t68 * rSges(4,3);
t52 = t76 * rSges(2,1) - t74 * rSges(2,2);
t51 = -t74 * rSges(2,1) - t76 * rSges(2,2);
t50 = t73 * rSges(4,1) + t75 * rSges(4,2);
t36 = t69 * rSges(3,1) - t68 * rSges(3,2) + t70;
t35 = -t68 * rSges(3,1) - t69 * rSges(3,2) - t118;
t13 = t111 * t69;
t12 = t111 * t68;
t11 = t96 * t69;
t10 = t96 * t68;
t9 = t77 + t98;
t8 = t115 + (-pkin(2) - t93) * t68 + t97;
t7 = t94 + t99;
t6 = t116 + (-pkin(2) + (-rSges(5,1) - pkin(3)) * t75 + (-rSges(5,3) - qJ(4)) * t73) * t68 + t97;
t5 = t69 * t77 + (t68 * t93 - t115) * t68;
t4 = t101 * t68 + t100 + t94;
t3 = t101 * t69 + (-pkin(2) + (-rSges(6,2) - qJ(4)) * t73 + (-pkin(3) - t117) * t75) * t68 + t97;
t2 = t69 * t99 + (-t116 + (rSges(5,1) * t75 + rSges(5,3) * t73) * t68) * t68 + t112;
t1 = t100 * t69 + (rSges(6,2) * t73 + t117 * t75) * t66 + t112;
t14 = [Icges(2,3) + Icges(3,3) + m(6) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) + m(2) * (t51 ^ 2 + t52 ^ 2) + m(3) * (t35 ^ 2 + t36 ^ 2) + m(4) * (t8 ^ 2 + t9 ^ 2) + (t143 * t73 + t107) * t73 + (t142 + t146 * t73 + (Icges(4,2) + Icges(6,2) + Icges(5,3)) * t75) * t75; 0; m(3) + m(4) - t119; m(6) * (t10 * t4 + t11 * t3) + m(5) * (t12 * t7 + t13 * t6) + m(4) * (-t68 * t9 - t69 * t8) * t50 + ((t83 * t135 + t137 * t69) * t69 + (t83 * t132 + t137 * t68) * t68) * t75 + ((t139 * t135 + t136 * t69) * t69 + (t139 * t132 + t136 * t68) * t68) * t73 + t109 * ((t127 + t144 / 0.2e1) * t75 + (t129 + t145 / 0.2e1) * t73); m(4) * t5 + m(5) * t2 + m(6) * t1; m(5) * (t12 ^ 2 + t13 ^ 2 + t2 ^ 2) + m(6) * (t1 ^ 2 + t10 ^ 2 + t11 ^ 2) + m(4) * (t109 * t50 ^ 2 + t5 ^ 2) + t122 * t68 * t66 + (t123 * t67 + (t122 * t69 + t123 * t68) * t68) * t69; 0.2e1 * ((t3 * t69 + t4 * t68) * t120 + (t6 * t69 + t68 * t7) * t121) * t73; t119 * t75; m(5) * (-t75 * t2 + (t12 * t68 + t13 * t69) * t73) + m(6) * (-t75 * t1 + (t10 * t68 + t11 * t69) * t73); 0.2e1 * (t121 + t120) * (t109 * t73 ^ 2 + t75 ^ 2); m(6) * (-t68 * t3 + t69 * t4); 0; m(6) * (t69 * t10 - t68 * t11); 0; m(6) * t109;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t14(1), t14(2), t14(4), t14(7), t14(11); t14(2), t14(3), t14(5), t14(8), t14(12); t14(4), t14(5), t14(6), t14(9), t14(13); t14(7), t14(8), t14(9), t14(10), t14(14); t14(11), t14(12), t14(13), t14(14), t14(15);];
Mq = res;
