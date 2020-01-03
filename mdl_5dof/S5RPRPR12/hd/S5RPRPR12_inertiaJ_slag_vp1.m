% Calculate joint inertia matrix for
% S5RPRPR12
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR12_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR12_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:27
% EndTime: 2019-12-31 18:29:31
% DurationCPUTime: 1.42s
% Computational Cost: add. (3075->270), mult. (3099->417), div. (0->0), fcn. (3229->10), ass. (0->126)
t112 = sin(qJ(1));
t152 = t112 / 0.2e1;
t113 = cos(qJ(1));
t150 = t113 / 0.2e1;
t110 = -pkin(7) - qJ(4);
t106 = sin(pkin(9));
t134 = t112 * t106;
t103 = pkin(8) + qJ(3);
t100 = cos(t103);
t137 = t100 * t113;
t98 = sin(t103);
t141 = t113 * t98;
t102 = pkin(9) + qJ(5);
t97 = sin(t102);
t99 = cos(t102);
t69 = t112 * t99 - t137 * t97;
t70 = t112 * t97 + t137 * t99;
t32 = t70 * rSges(6,1) + t69 * rSges(6,2) + rSges(6,3) * t141;
t108 = cos(pkin(9));
t94 = pkin(4) * t108 + pkin(3);
t156 = pkin(4) * t134 - t110 * t141 + t94 * t137 + t32;
t155 = t100 ^ 2;
t104 = t112 ^ 2;
t105 = t113 ^ 2;
t154 = m(5) / 0.2e1;
t153 = m(6) / 0.2e1;
t151 = -t113 / 0.2e1;
t149 = pkin(3) * t100;
t47 = -Icges(6,6) * t100 + (Icges(6,4) * t99 - Icges(6,2) * t97) * t98;
t148 = t47 * t97;
t82 = t98 * pkin(3) - t100 * qJ(4);
t147 = t100 * rSges(5,3) - (rSges(5,1) * t108 - rSges(5,2) * t106) * t98 - t82;
t122 = qJ(4) * t98 + t149;
t145 = pkin(3) * t137 + qJ(4) * t141;
t146 = t104 * t122 + t113 * t145;
t144 = Icges(4,4) * t98;
t143 = t100 * t94;
t142 = t112 * t98;
t140 = rSges(3,3) + qJ(2);
t139 = Icges(4,4) * t100;
t138 = t100 * t112;
t136 = t106 * t113;
t135 = t108 * t113;
t133 = t112 * t108;
t132 = t104 + t105;
t131 = t154 + t153;
t49 = -t100 * rSges(6,3) + (rSges(6,1) * t99 - rSges(6,2) * t97) * t98;
t130 = -(-pkin(3) + t94) * t98 - (qJ(4) + t110) * t100 - t49 - t82;
t75 = -t100 * t136 + t133;
t76 = t100 * t135 + t134;
t129 = t76 * rSges(5,1) + t75 * rSges(5,2) + rSges(5,3) * t141;
t67 = -t113 * t99 - t138 * t97;
t68 = -t113 * t97 + t138 * t99;
t25 = Icges(6,5) * t68 + Icges(6,6) * t67 + Icges(6,3) * t142;
t27 = Icges(6,4) * t68 + Icges(6,2) * t67 + Icges(6,6) * t142;
t29 = Icges(6,1) * t68 + Icges(6,4) * t67 + Icges(6,5) * t142;
t10 = -t100 * t25 + (-t27 * t97 + t29 * t99) * t98;
t46 = -Icges(6,3) * t100 + (Icges(6,5) * t99 - Icges(6,6) * t97) * t98;
t48 = -Icges(6,5) * t100 + (Icges(6,1) * t99 - Icges(6,4) * t97) * t98;
t13 = t142 * t46 + t47 * t67 + t48 * t68;
t128 = t13 / 0.2e1 + t10 / 0.2e1;
t26 = Icges(6,5) * t70 + Icges(6,6) * t69 + Icges(6,3) * t141;
t28 = Icges(6,4) * t70 + Icges(6,2) * t69 + Icges(6,6) * t141;
t30 = Icges(6,1) * t70 + Icges(6,4) * t69 + Icges(6,5) * t141;
t11 = -t100 * t26 + (-t28 * t97 + t30 * t99) * t98;
t14 = t141 * t46 + t47 * t69 + t48 * t70;
t127 = t14 / 0.2e1 + t11 / 0.2e1;
t111 = -pkin(6) - qJ(2);
t109 = cos(pkin(8));
t95 = pkin(2) * t109 + pkin(1);
t126 = -t111 * t112 + t113 * t95;
t73 = -t100 * t134 - t135;
t74 = t100 * t133 - t136;
t125 = -t74 * rSges(5,1) - t73 * rSges(5,2);
t124 = -t68 * rSges(6,1) - t67 * rSges(6,2);
t123 = rSges(4,1) * t100 - rSges(4,2) * t98;
t119 = Icges(4,1) * t100 - t144;
t118 = -Icges(4,2) * t98 + t139;
t117 = Icges(4,5) * t100 - Icges(4,6) * t98;
t116 = rSges(4,1) * t137 - rSges(4,2) * t141 + t112 * rSges(4,3);
t107 = sin(pkin(8));
t114 = rSges(3,1) * t109 - rSges(3,2) * t107 + pkin(1);
t85 = rSges(2,1) * t113 - rSges(2,2) * t112;
t84 = -rSges(2,1) * t112 - rSges(2,2) * t113;
t83 = rSges(4,1) * t98 + rSges(4,2) * t100;
t78 = Icges(4,5) * t98 + Icges(4,6) * t100;
t57 = Icges(4,3) * t112 + t113 * t117;
t56 = -Icges(4,3) * t113 + t112 * t117;
t54 = -Icges(5,5) * t100 + (Icges(5,1) * t108 - Icges(5,4) * t106) * t98;
t53 = -Icges(5,6) * t100 + (Icges(5,4) * t108 - Icges(5,2) * t106) * t98;
t51 = t112 * t140 + t113 * t114;
t50 = -t112 * t114 + t113 * t140;
t44 = t116 + t126;
t43 = (rSges(4,3) - t111) * t113 + (-t123 - t95) * t112;
t42 = t98 * t99 * t48;
t41 = t147 * t113;
t40 = t147 * t112;
t39 = Icges(5,1) * t76 + Icges(5,4) * t75 + Icges(5,5) * t141;
t38 = Icges(5,1) * t74 + Icges(5,4) * t73 + Icges(5,5) * t142;
t37 = Icges(5,4) * t76 + Icges(5,2) * t75 + Icges(5,6) * t141;
t36 = Icges(5,4) * t74 + Icges(5,2) * t73 + Icges(5,6) * t142;
t35 = Icges(5,5) * t76 + Icges(5,6) * t75 + Icges(5,3) * t141;
t34 = Icges(5,5) * t74 + Icges(5,6) * t73 + Icges(5,3) * t142;
t33 = t113 * t116 + (-t113 * rSges(4,3) + t112 * t123) * t112;
t31 = rSges(6,3) * t142 - t124;
t24 = t126 + t129 + t145;
t23 = -t113 * t111 + (-t149 - t95 + (-rSges(5,3) - qJ(4)) * t98) * t112 + t125;
t22 = t130 * t113;
t21 = t130 * t112;
t20 = -t100 * t32 - t141 * t49;
t19 = t100 * t31 + t142 * t49;
t18 = t126 + t156;
t17 = (pkin(4) * t106 - t111) * t113 + (-t143 - t95 + (-rSges(6,3) + t110) * t98) * t112 + t124;
t16 = -t100 * t46 - t148 * t98 + t42;
t15 = (-t112 * t32 + t113 * t31) * t98;
t12 = t112 * (rSges(5,3) * t142 - t125) + t113 * t129 + t146;
t9 = t141 * t26 + t28 * t69 + t30 * t70;
t8 = t141 * t25 + t27 * t69 + t29 * t70;
t7 = t142 * t26 + t28 * t67 + t30 * t68;
t6 = t142 * t25 + t27 * t67 + t29 * t68;
t5 = (-t145 + t156) * t113 + (-pkin(4) * t136 + t31 + (-t110 * t98 - t122 + t143) * t112) * t112 + t146;
t4 = t112 * t9 - t113 * t8;
t3 = t112 * t7 - t113 * t6;
t2 = -t14 * t100 + (t112 * t8 + t113 * t9) * t98;
t1 = -t13 * t100 + (t112 * t6 + t113 * t7) * t98;
t45 = [Icges(3,2) * t109 ^ 2 + Icges(2,3) + t42 + (Icges(3,1) * t107 + 0.2e1 * Icges(3,4) * t109) * t107 + (-t46 + t144 - (Icges(5,5) * t108 - Icges(5,6) * t106) * t98 + (Icges(4,2) + Icges(5,3)) * t100) * t100 + (Icges(4,1) * t98 - t106 * t53 + t108 * t54 + t139 - t148) * t98 + m(6) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2) + m(4) * (t43 ^ 2 + t44 ^ 2) + m(3) * (t50 ^ 2 + t51 ^ 2) + m(2) * (t84 ^ 2 + t85 ^ 2); m(6) * (t112 * t17 - t113 * t18) + m(5) * (t112 * t23 - t113 * t24) + m(4) * (t112 * t43 - t113 * t44) + m(3) * (t112 * t50 - t113 * t51); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t131) * t132; (-t73 * t53 / 0.2e1 - t74 * t54 / 0.2e1 + t78 * t150 - t128) * t113 + (t75 * t53 / 0.2e1 + t76 * t54 / 0.2e1 + t78 * t152 + t127) * t112 + m(6) * (t17 * t22 + t18 * t21) + m(5) * (t23 * t41 + t24 * t40) + m(4) * (-t112 * t44 - t113 * t43) * t83 + ((Icges(4,6) * t150 - t112 * t118 / 0.2e1 + t34 / 0.2e1) * t113 + (Icges(4,6) * t152 + t118 * t150 - t35 / 0.2e1) * t112) * t100 + ((Icges(4,5) * t112 - t106 * t37 + t108 * t39 + t113 * t119) * t152 + (-Icges(4,5) * t113 - t106 * t36 + t108 * t38 + t112 * t119) * t151) * t98; m(5) * (t112 * t41 - t113 * t40) + m(6) * (t112 * t22 - t113 * t21); m(6) * (t21 ^ 2 + t22 ^ 2 + t5 ^ 2) + m(5) * (t12 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(4) * (t132 * t83 ^ 2 + t33 ^ 2) + (-t105 * t56 - t3 + (t34 * t142 + t73 * t36 + t74 * t38) * t113) * t113 + (t4 + t104 * t57 + (t35 * t141 + t75 * t37 + t76 * t39) * t112 + (-t112 * t56 + t113 * t57 - t141 * t34 - t142 * t35 - t36 * t75 - t37 * t73 - t38 * t76 - t39 * t74) * t113) * t112; 0.2e1 * ((t112 * t18 + t113 * t17) * t153 + (t112 * t24 + t113 * t23) * t154) * t98; 0; m(6) * (-t100 * t5 + (t112 * t21 + t113 * t22) * t98) + m(5) * (-t100 * t12 + (t112 * t40 + t113 * t41) * t98); 0.2e1 * t131 * (t132 * t98 ^ 2 + t155); -t16 * t100 + m(6) * (t17 * t19 + t18 * t20) + (t112 * t128 + t113 * t127) * t98; m(6) * (t112 * t19 - t113 * t20); -t100 * (-t10 * t113 + t11 * t112) / 0.2e1 + m(6) * (t15 * t5 + t19 * t22 + t20 * t21) + t2 * t152 + t1 * t151 + (t150 * t4 + t152 * t3) * t98; m(6) * (-t15 * t100 + (t112 * t20 + t113 * t19) * t98); m(6) * (t15 ^ 2 + t19 ^ 2 + t20 ^ 2) + t155 * t16 + (t113 * t2 + t112 * t1 - t100 * (t10 * t112 + t11 * t113)) * t98;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t45(1), t45(2), t45(4), t45(7), t45(11); t45(2), t45(3), t45(5), t45(8), t45(12); t45(4), t45(5), t45(6), t45(9), t45(13); t45(7), t45(8), t45(9), t45(10), t45(14); t45(11), t45(12), t45(13), t45(14), t45(15);];
Mq = res;
