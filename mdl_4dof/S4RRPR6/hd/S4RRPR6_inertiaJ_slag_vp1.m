% Calculate joint inertia matrix for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR6_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR6_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:23
% EndTime: 2019-12-31 17:04:25
% DurationCPUTime: 0.78s
% Computational Cost: add. (1143->159), mult. (1200->235), div. (0->0), fcn. (1068->8), ass. (0->81)
t78 = qJ(2) + pkin(7);
t69 = sin(t78);
t70 = cos(t78);
t82 = sin(qJ(2));
t84 = cos(qJ(2));
t136 = Icges(3,5) * t84 + Icges(4,5) * t70 - Icges(3,6) * t82 - Icges(4,6) * t69;
t135 = Icges(3,3) + Icges(4,3);
t83 = sin(qJ(1));
t134 = t83 * pkin(5);
t85 = cos(qJ(1));
t133 = t85 * t83;
t71 = qJ(4) + t78;
t66 = sin(t71);
t67 = cos(t71);
t105 = rSges(5,1) * t67 - rSges(5,2) * t66;
t106 = rSges(4,1) * t70 - rSges(4,2) * t69;
t132 = t135 * t85 - t136 * t83;
t131 = t135 * t83 + t136 * t85;
t79 = t83 ^ 2;
t80 = t85 ^ 2;
t130 = t83 / 0.2e1;
t129 = -t85 / 0.2e1;
t128 = pkin(2) * t82;
t68 = t84 * pkin(2) + pkin(1);
t127 = rSges(3,1) * t84;
t124 = rSges(3,2) * t82;
t121 = t85 * rSges(3,3);
t81 = -qJ(3) - pkin(5);
t87 = t83 * rSges(5,3) + t105 * t85;
t9 = t83 * (-t85 * rSges(5,3) + t105 * t83) + t85 * t87;
t63 = t85 * t68;
t76 = t85 * pkin(5);
t120 = t83 * (t76 + (-pkin(1) + t68) * t83) + t85 * (-t85 * pkin(1) - t134 + t63);
t119 = t83 * rSges(3,3) + t85 * t127;
t118 = t79 + t80;
t117 = Icges(3,4) * t82;
t116 = Icges(3,4) * t84;
t115 = Icges(4,4) * t69;
t114 = Icges(4,4) * t70;
t113 = Icges(5,4) * t66;
t112 = Icges(5,4) * t67;
t44 = Icges(5,5) * t66 + Icges(5,6) * t67;
t92 = -Icges(5,2) * t66 + t112;
t95 = Icges(5,1) * t67 - t113;
t45 = Icges(5,2) * t67 + t113;
t46 = Icges(5,1) * t66 + t112;
t98 = -t45 * t66 + t46 * t67;
t111 = (t67 * (Icges(5,6) * t83 + t92 * t85) + t66 * (Icges(5,5) * t83 + t95 * t85) + t83 * t44 + t98 * t85) * t130 + (t67 * (-Icges(5,6) * t85 + t92 * t83) + t66 * (-Icges(5,5) * t85 + t95 * t83) - t85 * t44 + t98 * t83) * t129;
t89 = Icges(5,5) * t67 - Icges(5,6) * t66;
t23 = -Icges(5,3) * t85 + t89 * t83;
t24 = Icges(5,3) * t83 + t89 * t85;
t110 = -t85 * (-t24 * t133 + t80 * t23) + t83 * (-t23 * t133 + t79 * t24);
t109 = -t69 * rSges(4,1) - t70 * rSges(4,2) - t128;
t107 = -t124 + t127;
t97 = Icges(3,1) * t84 - t117;
t96 = Icges(4,1) * t70 - t115;
t94 = -Icges(3,2) * t82 + t116;
t93 = -Icges(4,2) * t69 + t114;
t88 = t83 * rSges(4,3) + t106 * t85;
t47 = t66 * rSges(5,1) + t67 * rSges(5,2);
t86 = -pkin(3) * t69 - t128 - t47;
t77 = -pkin(6) + t81;
t61 = t85 * rSges(2,1) - t83 * rSges(2,2);
t60 = -t83 * rSges(2,1) - t85 * rSges(2,2);
t59 = t82 * rSges(3,1) + t84 * rSges(3,2);
t53 = pkin(3) * t70 + t68;
t48 = t85 * t53;
t30 = t109 * t85;
t29 = t109 * t83;
t22 = t134 + (pkin(1) - t124) * t85 + t119;
t21 = t121 + t76 + (-pkin(1) - t107) * t83;
t16 = t86 * t85;
t15 = t86 * t83;
t14 = -t83 * t81 + t63 + t88;
t13 = (rSges(4,3) - t81) * t85 + (-t106 - t68) * t83;
t12 = t85 * (-t85 * t124 + t119) + (t107 * t83 - t121) * t83;
t11 = -t83 * t77 + t48 + t87;
t10 = (rSges(5,3) - t77) * t85 + (-t105 - t53) * t83;
t4 = t85 * t88 + (-t85 * rSges(4,3) + t106 * t83) * t83 + t120;
t3 = t85 * (t48 - t63) + (t53 - t68) * t79 + t120 + t9;
t1 = [t67 * t45 + t66 * t46 + t70 * (Icges(4,2) * t70 + t115) + t69 * (Icges(4,1) * t69 + t114) + t84 * (Icges(3,2) * t84 + t117) + t82 * (Icges(3,1) * t82 + t116) + Icges(2,3) + m(2) * (t60 ^ 2 + t61 ^ 2) + m(3) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2); m(4) * (t30 * t13 + t29 * t14) + m(5) * (t16 * t10 + t15 * t11) + m(3) * (-t21 * t85 - t22 * t83) * t59 + t111 + (t70 * (Icges(4,6) * t83 + t93 * t85) + t69 * (Icges(4,5) * t83 + t96 * t85) + t84 * (Icges(3,6) * t83 + t94 * t85) + t82 * (Icges(3,5) * t83 + t97 * t85)) * t130 + (t70 * (-Icges(4,6) * t85 + t93 * t83) + t69 * (-Icges(4,5) * t85 + t96 * t83) + t84 * (-Icges(3,6) * t85 + t94 * t83) + t82 * (-Icges(3,5) * t85 + t97 * t83)) * t129 + (Icges(3,5) * t82 + Icges(4,5) * t69 + Icges(3,6) * t84 + Icges(4,6) * t70) * (t79 / 0.2e1 + t80 / 0.2e1); m(5) * (t15 ^ 2 + t16 ^ 2 + t3 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2 + t4 ^ 2) + m(3) * (t118 * t59 ^ 2 + t12 ^ 2) + t110 + t132 * t85 * t80 + (t131 * t79 + (t131 * t85 + t132 * t83) * t85) * t83; m(4) * (t83 * t13 - t85 * t14) + m(5) * (t83 * t10 - t85 * t11); m(5) * (-t85 * t15 + t83 * t16) + m(4) * (-t85 * t29 + t83 * t30); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t118; m(5) * (-t10 * t85 - t11 * t83) * t47 + t111; m(5) * (t9 * t3 + (-t15 * t83 - t16 * t85) * t47) + t110; 0; m(5) * (t118 * t47 ^ 2 + t9 ^ 2) + t110;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
