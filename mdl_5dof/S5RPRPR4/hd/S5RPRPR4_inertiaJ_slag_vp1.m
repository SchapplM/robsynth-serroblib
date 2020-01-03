% Calculate joint inertia matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:15
% EndTime: 2020-01-03 11:38:21
% DurationCPUTime: 1.00s
% Computational Cost: add. (1929->181), mult. (1340->249), div. (0->0), fcn. (1180->10), ass. (0->92)
t88 = qJ(3) + pkin(9);
t79 = sin(t88);
t81 = cos(t88);
t91 = sin(qJ(3));
t93 = cos(qJ(3));
t144 = Icges(4,5) * t93 + Icges(5,5) * t81 - Icges(4,6) * t91 - Icges(5,6) * t79;
t143 = Icges(4,3) + Icges(5,3);
t89 = qJ(1) + pkin(8);
t80 = sin(t89);
t82 = cos(t89);
t142 = t80 * t82;
t83 = qJ(5) + t88;
t74 = sin(t83);
t75 = cos(t83);
t113 = -rSges(6,1) * t75 + rSges(6,2) * t74;
t141 = rSges(5,1) * t81 - rSges(5,2) * t79;
t140 = t143 * t82 - t144 * t80;
t139 = t143 * t80 + t144 * t82;
t77 = t80 ^ 2;
t78 = t82 ^ 2;
t138 = -t80 / 0.2e1;
t137 = -t82 / 0.2e1;
t136 = pkin(3) * t91;
t76 = t93 * pkin(3) + pkin(2);
t135 = rSges(4,1) * t93;
t132 = rSges(4,2) * t91;
t90 = -qJ(4) - pkin(6);
t56 = pkin(4) * t81 + t76;
t87 = -pkin(7) + t90;
t129 = t80 * t56 + t82 * t87;
t128 = t80 * t76 + t82 * t90;
t127 = t82 * pkin(2) + t80 * pkin(6);
t126 = t78 + t77;
t125 = Icges(4,4) * t91;
t124 = Icges(4,4) * t93;
t123 = Icges(5,4) * t79;
t122 = Icges(5,4) * t81;
t121 = Icges(6,4) * t74;
t120 = Icges(6,4) * t75;
t100 = -Icges(6,2) * t74 + t120;
t103 = Icges(6,1) * t75 - t121;
t48 = Icges(6,2) * t75 + t121;
t49 = Icges(6,1) * t74 + t120;
t106 = t48 * t74 - t49 * t75;
t47 = Icges(6,5) * t74 + Icges(6,6) * t75;
t119 = (t106 * t82 + t75 * (-Icges(6,6) * t80 - t100 * t82) + t74 * (-Icges(6,5) * t80 - t103 * t82) - t80 * t47) * t138 + (-t106 * t80 + t75 * (-Icges(6,6) * t82 + t100 * t80) + t74 * (-Icges(6,5) * t82 + t103 * t80) - t82 * t47) * t137;
t50 = t74 * rSges(6,1) + t75 * rSges(6,2);
t118 = pkin(4) * t79 + t50;
t116 = t141 * t80;
t97 = Icges(6,5) * t75 - Icges(6,6) * t74;
t22 = -Icges(6,3) * t82 + t80 * t97;
t23 = -Icges(6,3) * t80 - t82 * t97;
t115 = -t82 * (t23 * t142 + t78 * t22) - t80 * (t22 * t142 + t77 * t23);
t114 = -t132 + t135;
t105 = Icges(4,1) * t93 - t125;
t104 = Icges(5,1) * t81 - t123;
t102 = -Icges(4,2) * t91 + t124;
t101 = -Icges(5,2) * t79 + t122;
t96 = t141 * t82;
t95 = -t82 * rSges(6,3) - t113 * t80;
t94 = cos(qJ(1));
t92 = sin(qJ(1));
t86 = t94 * pkin(1);
t84 = t92 * pkin(1);
t70 = t82 * t136;
t68 = t80 * t135;
t67 = t94 * rSges(2,1) - t92 * rSges(2,2);
t66 = t92 * rSges(2,1) + t94 * rSges(2,2);
t65 = t91 * rSges(4,1) + t93 * rSges(4,2);
t60 = t82 * t76;
t54 = t79 * rSges(5,1) + t81 * rSges(5,2);
t44 = t82 * rSges(3,1) - t80 * rSges(3,2) + t86;
t43 = t80 * rSges(3,1) + t82 * rSges(3,2) + t84;
t36 = t82 * t54 + t70;
t35 = (-t54 - t136) * t80;
t28 = -t80 * rSges(6,3) + t113 * t82;
t21 = t80 * t90 + t127 - t60;
t20 = t80 * t95;
t19 = t80 * (-t80 * pkin(2) + t82 * pkin(6) + t128);
t18 = t80 * rSges(4,3) + t114 * t82 + t127 + t86;
t17 = t68 + t84 + (-rSges(4,3) - pkin(6)) * t82 + (pkin(2) - t132) * t80;
t16 = t118 * t82 + t70;
t15 = (-t118 - t136) * t80;
t14 = t60 + t86 + t96 + (rSges(5,3) - t90) * t80;
t13 = -t82 * rSges(5,3) + t116 + t128 + t84;
t12 = t80 * (-t80 * t132 + t68) + t114 * t78;
t11 = t86 + (rSges(6,3) - t87) * t80 + (-t113 + t56) * t82;
t10 = t84 + t95 + t129;
t9 = -t82 * t28 + t20;
t4 = t19 + t80 * t116 + (-t21 + t96) * t82;
t3 = t19 + t80 * (-t128 + t129) + t20 + (t82 * t56 - t21 - t28 - t60 + (-t87 + t90) * t80) * t82;
t1 = [t75 * t48 + t74 * t49 + t81 * (Icges(5,2) * t81 + t123) + t79 * (Icges(5,1) * t79 + t122) + t93 * (Icges(4,2) * t93 + t125) + t91 * (Icges(4,1) * t91 + t124) + Icges(2,3) + Icges(3,3) + m(2) * (t66 ^ 2 + t67 ^ 2) + m(3) * (t43 ^ 2 + t44 ^ 2) + m(4) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2); 0; m(3) + m(4) + m(5) + m(6); m(5) * (t36 * t13 + t35 * t14) + m(6) * (t16 * t10 + t15 * t11) + m(4) * (t17 * t82 - t18 * t80) * t65 + t119 + (t81 * (-Icges(5,6) * t80 - t101 * t82) + t79 * (-Icges(5,5) * t80 - t104 * t82) + t93 * (-Icges(4,6) * t80 - t102 * t82) + t91 * (-Icges(4,5) * t80 - t105 * t82)) * t138 + (t81 * (-Icges(5,6) * t82 + t101 * t80) + t79 * (-Icges(5,5) * t82 + t104 * t80) + t93 * (-Icges(4,6) * t82 + t102 * t80) + t91 * (-Icges(4,5) * t82 + t105 * t80)) * t137 + (Icges(4,5) * t91 + Icges(5,5) * t79 + Icges(4,6) * t93 + Icges(5,6) * t81) * (t78 / 0.2e1 + t77 / 0.2e1); m(4) * t12 + m(5) * t4 + m(6) * t3; m(6) * (t15 ^ 2 + t16 ^ 2 + t3 ^ 2) + m(5) * (t35 ^ 2 + t36 ^ 2 + t4 ^ 2) + m(4) * (t126 * t65 ^ 2 + t12 ^ 2) + t115 + t140 * t82 * t78 + (t139 * t77 + (t139 * t82 + t140 * t80) * t82) * t80; m(5) * (-t80 * t13 - t82 * t14) + m(6) * (-t80 * t10 - t82 * t11); 0; m(6) * (-t82 * t15 - t80 * t16) + m(5) * (-t82 * t35 - t80 * t36); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t126; m(6) * (t10 * t82 - t11 * t80) * t50 + t119; m(6) * t9; m(6) * (t9 * t3 + (-t15 * t80 + t16 * t82) * t50) + t115; 0; m(6) * (t126 * t50 ^ 2 + t9 ^ 2) + t115;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
