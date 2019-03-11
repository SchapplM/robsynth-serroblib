% Calculate joint inertia matrix for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:39:08
% EndTime: 2019-03-09 01:39:11
% DurationCPUTime: 1.42s
% Computational Cost: add. (4561->290), mult. (3287->431), div. (0->0), fcn. (3401->12), ass. (0->135)
t111 = qJ(1) + pkin(9);
t104 = sin(t111);
t160 = t104 / 0.2e1;
t107 = cos(t111);
t158 = t107 / 0.2e1;
t110 = pkin(10) + qJ(4);
t106 = cos(t110);
t141 = t106 * t107;
t112 = sin(pkin(11));
t143 = t104 * t112;
t103 = sin(t110);
t116 = -pkin(8) - qJ(5);
t145 = t103 * t116;
t146 = t103 * t107;
t109 = pkin(11) + qJ(6);
t102 = sin(t109);
t105 = cos(t109);
t69 = -t102 * t141 + t104 * t105;
t70 = t102 * t104 + t105 * t141;
t33 = t70 * rSges(7,1) + t69 * rSges(7,2) + rSges(7,3) * t146;
t114 = cos(pkin(11));
t97 = pkin(5) * t114 + pkin(4);
t165 = pkin(5) * t143 - t107 * t145 + t97 * t141 + t33;
t100 = t104 ^ 2;
t164 = t106 ^ 2;
t101 = t107 ^ 2;
t163 = m(6) / 0.2e1;
t162 = m(7) / 0.2e1;
t161 = -m(6) - m(7);
t159 = -t107 / 0.2e1;
t118 = sin(qJ(1));
t157 = pkin(1) * t118;
t156 = pkin(4) * t106;
t128 = qJ(5) * t103 + t156;
t153 = pkin(4) * t141 + qJ(5) * t146;
t155 = t100 * t128 + t107 * t153;
t84 = pkin(4) * t103 - qJ(5) * t106;
t154 = rSges(6,3) * t106 - (rSges(6,1) * t114 - rSges(6,2) * t112) * t103 - t84;
t51 = -Icges(7,6) * t106 + (Icges(7,4) * t105 - Icges(7,2) * t102) * t103;
t152 = t102 * t51;
t151 = t106 * t97;
t150 = rSges(4,3) + qJ(3);
t149 = Icges(5,4) * t103;
t148 = Icges(5,4) * t106;
t147 = t103 * t104;
t144 = t104 * t106;
t142 = t104 * t114;
t140 = t107 * t112;
t139 = t107 * t114;
t138 = t100 + t101;
t137 = t163 + t162;
t57 = -rSges(7,3) * t106 + (rSges(7,1) * t105 - rSges(7,2) * t102) * t103;
t136 = -(qJ(5) + t116) * t106 - (-pkin(4) + t97) * t103 - t57 - t84;
t75 = -t106 * t140 + t142;
t76 = t106 * t139 + t143;
t135 = t76 * rSges(6,1) + t75 * rSges(6,2) + rSges(6,3) * t146;
t27 = Icges(7,5) * t70 + Icges(7,6) * t69 + Icges(7,3) * t146;
t29 = Icges(7,4) * t70 + Icges(7,2) * t69 + Icges(7,6) * t146;
t31 = Icges(7,1) * t70 + Icges(7,4) * t69 + Icges(7,5) * t146;
t11 = -t106 * t27 + (-t102 * t29 + t105 * t31) * t103;
t48 = -Icges(7,3) * t106 + (Icges(7,5) * t105 - Icges(7,6) * t102) * t103;
t54 = -Icges(7,5) * t106 + (Icges(7,1) * t105 - Icges(7,4) * t102) * t103;
t15 = t48 * t146 + t51 * t69 + t54 * t70;
t134 = t11 / 0.2e1 + t15 / 0.2e1;
t67 = -t102 * t144 - t105 * t107;
t68 = -t102 * t107 + t105 * t144;
t26 = Icges(7,5) * t68 + Icges(7,6) * t67 + Icges(7,3) * t147;
t28 = Icges(7,4) * t68 + Icges(7,2) * t67 + Icges(7,6) * t147;
t30 = Icges(7,1) * t68 + Icges(7,4) * t67 + Icges(7,5) * t147;
t10 = -t106 * t26 + (-t102 * t28 + t105 * t30) * t103;
t14 = t48 * t147 + t51 * t67 + t54 * t68;
t133 = t14 / 0.2e1 + t10 / 0.2e1;
t119 = cos(qJ(1));
t108 = t119 * pkin(1);
t117 = -pkin(7) - qJ(3);
t115 = cos(pkin(10));
t98 = pkin(3) * t115 + pkin(2);
t132 = -t104 * t117 + t107 * t98 + t108;
t73 = -t106 * t143 - t139;
t74 = t106 * t142 - t140;
t131 = -rSges(6,1) * t74 - rSges(6,2) * t73;
t130 = -rSges(7,1) * t68 - rSges(7,2) * t67;
t129 = rSges(5,1) * t106 - rSges(5,2) * t103;
t125 = Icges(5,1) * t106 - t149;
t124 = -Icges(5,2) * t103 + t148;
t123 = Icges(5,5) * t106 - Icges(5,6) * t103;
t122 = rSges(5,1) * t141 - rSges(5,2) * t146 + t104 * rSges(5,3);
t113 = sin(pkin(10));
t120 = rSges(4,1) * t115 - rSges(4,2) * t113 + pkin(2);
t94 = rSges(2,1) * t119 - t118 * rSges(2,2);
t93 = -t118 * rSges(2,1) - rSges(2,2) * t119;
t85 = rSges(5,1) * t103 + rSges(5,2) * t106;
t81 = Icges(5,5) * t103 + Icges(5,6) * t106;
t78 = rSges(3,1) * t107 - rSges(3,2) * t104 + t108;
t77 = -rSges(3,1) * t104 - rSges(3,2) * t107 - t157;
t65 = -Icges(6,5) * t106 + (Icges(6,1) * t114 - Icges(6,4) * t112) * t103;
t64 = -Icges(6,6) * t106 + (Icges(6,4) * t114 - Icges(6,2) * t112) * t103;
t50 = Icges(5,3) * t104 + t123 * t107;
t49 = -Icges(5,3) * t107 + t123 * t104;
t46 = t103 * t105 * t54;
t45 = t150 * t104 + t120 * t107 + t108;
t44 = -t120 * t104 + t150 * t107 - t157;
t43 = t154 * t107;
t42 = t154 * t104;
t41 = t122 + t132;
t40 = -t157 + (rSges(5,3) - t117) * t107 + (-t129 - t98) * t104;
t39 = Icges(6,1) * t76 + Icges(6,4) * t75 + Icges(6,5) * t146;
t38 = Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t147;
t37 = Icges(6,4) * t76 + Icges(6,2) * t75 + Icges(6,6) * t146;
t36 = Icges(6,4) * t74 + Icges(6,2) * t73 + Icges(6,6) * t147;
t35 = Icges(6,5) * t76 + Icges(6,6) * t75 + Icges(6,3) * t146;
t34 = Icges(6,5) * t74 + Icges(6,6) * t73 + Icges(6,3) * t147;
t32 = rSges(7,3) * t147 - t130;
t25 = t107 * t122 + (-t107 * rSges(5,3) + t129 * t104) * t104;
t24 = t136 * t107;
t23 = t136 * t104;
t22 = t132 + t135 + t153;
t21 = -t157 - t107 * t117 + (-t156 - t98 + (-rSges(6,3) - qJ(5)) * t103) * t104 + t131;
t20 = -t106 * t33 - t57 * t146;
t19 = t106 * t32 + t57 * t147;
t18 = -t103 * t152 - t106 * t48 + t46;
t17 = t132 + t165;
t16 = -t157 + (pkin(5) * t112 - t117) * t107 + (-t151 - t98 + (-rSges(7,3) + t116) * t103) * t104 + t130;
t13 = (-t104 * t33 + t107 * t32) * t103;
t12 = t104 * (rSges(6,3) * t147 - t131) + t107 * t135 + t155;
t9 = t27 * t146 + t29 * t69 + t31 * t70;
t8 = t26 * t146 + t28 * t69 + t30 * t70;
t7 = t27 * t147 + t29 * t67 + t31 * t68;
t6 = t26 * t147 + t28 * t67 + t30 * t68;
t5 = (-t153 + t165) * t107 + (-pkin(5) * t140 + t32 + (-t128 - t145 + t151) * t104) * t104 + t155;
t4 = t104 * t9 - t107 * t8;
t3 = t104 * t7 - t107 * t6;
t2 = -t106 * t15 + (t104 * t8 + t107 * t9) * t103;
t1 = -t106 * t14 + (t104 * t6 + t107 * t7) * t103;
t47 = [Icges(4,2) * t115 ^ 2 + Icges(2,3) + Icges(3,3) + t46 + (Icges(4,1) * t113 + 0.2e1 * Icges(4,4) * t115) * t113 + (-t48 + t149 - (Icges(6,5) * t114 - Icges(6,6) * t112) * t103 + (Icges(5,2) + Icges(6,3)) * t106) * t106 + (Icges(5,1) * t103 - t112 * t64 + t114 * t65 + t148 - t152) * t103 + m(7) * (t16 ^ 2 + t17 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t40 ^ 2 + t41 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + m(3) * (t77 ^ 2 + t78 ^ 2) + m(2) * (t93 ^ 2 + t94 ^ 2); 0; m(3) + m(4) + m(5) - t161; m(7) * (t104 * t16 - t107 * t17) + m(6) * (t104 * t21 - t107 * t22) + m(5) * (t104 * t40 - t107 * t41) + m(4) * (t104 * t44 - t107 * t45); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t137) * t138; (-t73 * t64 / 0.2e1 - t74 * t65 / 0.2e1 + t81 * t158 - t133) * t107 + (t75 * t64 / 0.2e1 + t76 * t65 / 0.2e1 + t81 * t160 + t134) * t104 + m(7) * (t16 * t24 + t17 * t23) + m(6) * (t21 * t43 + t22 * t42) + m(5) * (-t104 * t41 - t107 * t40) * t85 + ((Icges(5,6) * t158 - t124 * t104 / 0.2e1 + t34 / 0.2e1) * t107 + (Icges(5,6) * t160 + t124 * t158 - t35 / 0.2e1) * t104) * t106 + ((Icges(5,5) * t104 + t125 * t107 - t112 * t37 + t114 * t39) * t160 + (-Icges(5,5) * t107 + t125 * t104 - t112 * t36 + t114 * t38) * t159) * t103; m(5) * t25 + m(6) * t12 + m(7) * t5; m(6) * (t104 * t43 - t107 * t42) + m(7) * (t104 * t24 - t107 * t23); m(7) * (t23 ^ 2 + t24 ^ 2 + t5 ^ 2) + m(6) * (t12 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t138 * t85 ^ 2 + t25 ^ 2) + (-t101 * t49 - t3 + (t34 * t147 + t73 * t36 + t74 * t38) * t107) * t107 + (t4 + t100 * t50 + (t35 * t146 + t75 * t37 + t76 * t39) * t104 + (-t104 * t49 + t107 * t50 - t34 * t146 - t35 * t147 - t36 * t75 - t37 * t73 - t38 * t76 - t39 * t74) * t107) * t104; 0.2e1 * ((t104 * t17 + t107 * t16) * t162 + (t104 * t22 + t107 * t21) * t163) * t103; t161 * t106; 0; m(7) * (-t106 * t5 + (t104 * t23 + t107 * t24) * t103) + m(6) * (-t106 * t12 + (t104 * t42 + t107 * t43) * t103); 0.2e1 * t137 * (t138 * t103 ^ 2 + t164); m(7) * (t16 * t19 + t17 * t20) - t18 * t106 + (t133 * t104 + t134 * t107) * t103; m(7) * t13; m(7) * (t104 * t19 - t107 * t20); t2 * t160 + t1 * t159 - t106 * (-t10 * t107 + t11 * t104) / 0.2e1 + m(7) * (t13 * t5 + t19 * t24 + t20 * t23) + (t4 * t158 + t3 * t160) * t103; m(7) * (-t106 * t13 + (t104 * t20 + t107 * t19) * t103); t164 * t18 + m(7) * (t13 ^ 2 + t19 ^ 2 + t20 ^ 2) + (t107 * t2 + t104 * t1 - t106 * (t10 * t104 + t107 * t11)) * t103;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t47(1) t47(2) t47(4) t47(7) t47(11) t47(16); t47(2) t47(3) t47(5) t47(8) t47(12) t47(17); t47(4) t47(5) t47(6) t47(9) t47(13) t47(18); t47(7) t47(8) t47(9) t47(10) t47(14) t47(19); t47(11) t47(12) t47(13) t47(14) t47(15) t47(20); t47(16) t47(17) t47(18) t47(19) t47(20) t47(21);];
Mq  = res;
