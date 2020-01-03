% Calculate joint inertia matrix for
% S5RPRPR8
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:18
% EndTime: 2019-12-31 18:21:22
% DurationCPUTime: 1.27s
% Computational Cost: add. (3135->259), mult. (2979->389), div. (0->0), fcn. (3131->10), ass. (0->128)
t105 = qJ(1) + pkin(8);
t100 = sin(t105);
t154 = t100 / 0.2e1;
t102 = cos(t105);
t152 = t102 / 0.2e1;
t109 = -pkin(7) - qJ(4);
t110 = sin(qJ(3));
t132 = t109 * t110;
t112 = cos(qJ(3));
t135 = t102 * t112;
t107 = sin(pkin(9));
t140 = t100 * t107;
t136 = t102 * t110;
t104 = pkin(9) + qJ(5);
t101 = cos(t104);
t99 = sin(t104);
t58 = t100 * t101 - t99 * t135;
t59 = t100 * t99 + t101 * t135;
t33 = t59 * rSges(6,1) + t58 * rSges(6,2) + rSges(6,3) * t136;
t108 = cos(pkin(9));
t96 = t108 * pkin(4) + pkin(3);
t159 = pkin(4) * t140 - t102 * t132 + t96 * t135 + t33;
t97 = t100 ^ 2;
t98 = t102 ^ 2;
t158 = t112 ^ 2;
t157 = m(5) / 0.2e1;
t156 = m(6) / 0.2e1;
t155 = -m(5) - m(6);
t153 = -t102 / 0.2e1;
t151 = pkin(3) * t112;
t111 = sin(qJ(1));
t150 = t111 * pkin(1);
t61 = -Icges(6,6) * t112 + (Icges(6,4) * t101 - Icges(6,2) * t99) * t110;
t149 = t99 * t61;
t121 = qJ(4) * t110 + t151;
t146 = pkin(3) * t135 + qJ(4) * t136;
t148 = t102 * t146 + t97 * t121;
t82 = t110 * pkin(3) - t112 * qJ(4);
t147 = t112 * rSges(5,3) - (rSges(5,1) * t108 - rSges(5,2) * t107) * t110 - t82;
t145 = t97 + t98;
t144 = t102 * rSges(4,3);
t143 = t112 * t96;
t142 = Icges(4,4) * t110;
t141 = Icges(4,4) * t112;
t139 = t100 * t110;
t138 = t100 * t112;
t137 = t102 * t107;
t134 = t107 * t112;
t133 = t108 * t112;
t63 = -t112 * rSges(6,3) + (rSges(6,1) * t101 - rSges(6,2) * t99) * t110;
t131 = -(qJ(4) + t109) * t112 - (-pkin(3) + t96) * t110 - t63 - t82;
t69 = t100 * t108 - t102 * t134;
t70 = t102 * t133 + t140;
t130 = t70 * rSges(5,1) + t69 * rSges(5,2) + rSges(5,3) * t136;
t129 = pkin(4) * t137;
t113 = cos(qJ(1));
t103 = t113 * pkin(1);
t128 = t102 * pkin(2) + t100 * pkin(6) + t103;
t27 = Icges(6,5) * t59 + Icges(6,6) * t58 + Icges(6,3) * t136;
t29 = Icges(6,4) * t59 + Icges(6,2) * t58 + Icges(6,6) * t136;
t31 = Icges(6,1) * t59 + Icges(6,4) * t58 + Icges(6,5) * t136;
t11 = -t112 * t27 + (t101 * t31 - t29 * t99) * t110;
t60 = -Icges(6,3) * t112 + (Icges(6,5) * t101 - Icges(6,6) * t99) * t110;
t62 = -Icges(6,5) * t112 + (Icges(6,1) * t101 - Icges(6,4) * t99) * t110;
t15 = t60 * t136 + t58 * t61 + t59 * t62;
t127 = t11 / 0.2e1 + t15 / 0.2e1;
t56 = -t102 * t101 - t99 * t138;
t57 = t101 * t138 - t102 * t99;
t26 = Icges(6,5) * t57 + Icges(6,6) * t56 + Icges(6,3) * t139;
t28 = Icges(6,4) * t57 + Icges(6,2) * t56 + Icges(6,6) * t139;
t30 = Icges(6,1) * t57 + Icges(6,4) * t56 + Icges(6,5) * t139;
t10 = -t112 * t26 + (t101 * t30 - t28 * t99) * t110;
t14 = t60 * t139 + t56 * t61 + t57 * t62;
t126 = t14 / 0.2e1 + t10 / 0.2e1;
t125 = t102 * pkin(6) - t150;
t67 = -t100 * t134 - t102 * t108;
t68 = t100 * t133 - t137;
t124 = -t68 * rSges(5,1) - t67 * rSges(5,2);
t123 = -t57 * rSges(6,1) - t56 * rSges(6,2);
t122 = rSges(4,1) * t112 - rSges(4,2) * t110;
t118 = Icges(4,1) * t112 - t142;
t117 = -Icges(4,2) * t110 + t141;
t116 = Icges(4,5) * t112 - Icges(4,6) * t110;
t115 = rSges(4,1) * t135 - rSges(4,2) * t136 + t100 * rSges(4,3);
t85 = t113 * rSges(2,1) - t111 * rSges(2,2);
t84 = -t111 * rSges(2,1) - t113 * rSges(2,2);
t83 = t110 * rSges(4,1) + t112 * rSges(4,2);
t79 = Icges(4,5) * t110 + Icges(4,6) * t112;
t76 = t102 * rSges(3,1) - t100 * rSges(3,2) + t103;
t75 = -t100 * rSges(3,1) - t102 * rSges(3,2) - t150;
t73 = -Icges(5,5) * t112 + (Icges(5,1) * t108 - Icges(5,4) * t107) * t110;
t72 = -Icges(5,6) * t112 + (Icges(5,4) * t108 - Icges(5,2) * t107) * t110;
t46 = Icges(4,3) * t100 + t116 * t102;
t45 = -Icges(4,3) * t102 + t116 * t100;
t44 = t110 * t101 * t62;
t43 = t147 * t102;
t42 = t147 * t100;
t41 = t115 + t128;
t40 = t144 + (-pkin(2) - t122) * t100 + t125;
t39 = Icges(5,1) * t70 + Icges(5,4) * t69 + Icges(5,5) * t136;
t38 = Icges(5,1) * t68 + Icges(5,4) * t67 + Icges(5,5) * t139;
t37 = Icges(5,4) * t70 + Icges(5,2) * t69 + Icges(5,6) * t136;
t36 = Icges(5,4) * t68 + Icges(5,2) * t67 + Icges(5,6) * t139;
t35 = Icges(5,5) * t70 + Icges(5,6) * t69 + Icges(5,3) * t136;
t34 = Icges(5,5) * t68 + Icges(5,6) * t67 + Icges(5,3) * t139;
t32 = rSges(6,3) * t139 - t123;
t25 = t102 * t115 + (t122 * t100 - t144) * t100;
t24 = t131 * t102;
t23 = t131 * t100;
t22 = t128 + t130 + t146;
t21 = (-t151 - pkin(2) + (-rSges(5,3) - qJ(4)) * t110) * t100 + t124 + t125;
t20 = -t110 * t149 - t112 * t60 + t44;
t19 = -t112 * t33 - t63 * t136;
t18 = t112 * t32 + t63 * t139;
t17 = t128 + t159;
t16 = t129 + (-t143 - pkin(2) + (-rSges(6,3) + t109) * t110) * t100 + t123 + t125;
t13 = (-t100 * t33 + t102 * t32) * t110;
t12 = t100 * (rSges(5,3) * t139 - t124) + t102 * t130 + t148;
t9 = t27 * t136 + t58 * t29 + t59 * t31;
t8 = t26 * t136 + t58 * t28 + t59 * t30;
t7 = t27 * t139 + t56 * t29 + t57 * t31;
t6 = t26 * t139 + t56 * t28 + t57 * t30;
t5 = (-t146 + t159) * t102 + (-t129 + t32 + (-t121 - t132 + t143) * t100) * t100 + t148;
t4 = t9 * t100 - t8 * t102;
t3 = t7 * t100 - t6 * t102;
t2 = -t15 * t112 + (t100 * t8 + t102 * t9) * t110;
t1 = -t14 * t112 + (t100 * t6 + t102 * t7) * t110;
t47 = [Icges(2,3) + Icges(3,3) + t44 + (-t60 + t142 - (Icges(5,5) * t108 - Icges(5,6) * t107) * t110 + (Icges(4,2) + Icges(5,3)) * t112) * t112 + (Icges(4,1) * t110 - t107 * t72 + t108 * t73 + t141 - t149) * t110 + m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t40 ^ 2 + t41 ^ 2) + m(3) * (t75 ^ 2 + t76 ^ 2) + m(2) * (t84 ^ 2 + t85 ^ 2); 0; m(3) + m(4) - t155; (-t67 * t72 / 0.2e1 - t68 * t73 / 0.2e1 + t79 * t152 - t126) * t102 + (t69 * t72 / 0.2e1 + t70 * t73 / 0.2e1 + t79 * t154 + t127) * t100 + m(6) * (t24 * t16 + t23 * t17) + m(5) * (t43 * t21 + t42 * t22) + m(4) * (-t100 * t41 - t102 * t40) * t83 + ((Icges(4,6) * t152 - t117 * t100 / 0.2e1 + t34 / 0.2e1) * t102 + (Icges(4,6) * t154 + t117 * t152 - t35 / 0.2e1) * t100) * t112 + ((Icges(4,5) * t100 + t118 * t102 - t107 * t37 + t108 * t39) * t154 + (-Icges(4,5) * t102 + t118 * t100 - t107 * t36 + t108 * t38) * t153) * t110; m(4) * t25 + m(5) * t12 + m(6) * t5; m(6) * (t23 ^ 2 + t24 ^ 2 + t5 ^ 2) + m(5) * (t12 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(4) * (t145 * t83 ^ 2 + t25 ^ 2) + (-t98 * t45 - t3 + (t34 * t139 + t67 * t36 + t68 * t38) * t102) * t102 + (t4 + t97 * t46 + (t35 * t136 + t69 * t37 + t70 * t39) * t100 + (-t100 * t45 + t102 * t46 - t34 * t136 - t35 * t139 - t69 * t36 - t67 * t37 - t70 * t38 - t68 * t39) * t102) * t100; 0.2e1 * ((t100 * t17 + t102 * t16) * t156 + (t100 * t22 + t102 * t21) * t157) * t110; t155 * t112; m(6) * (-t112 * t5 + (t100 * t23 + t102 * t24) * t110) + m(5) * (-t112 * t12 + (t100 * t42 + t102 * t43) * t110); 0.2e1 * (t157 + t156) * (t145 * t110 ^ 2 + t158); m(6) * (t18 * t16 + t19 * t17) - t20 * t112 + (t126 * t100 + t127 * t102) * t110; m(6) * t13; m(6) * (t13 * t5 + t18 * t24 + t19 * t23) + t1 * t153 - t112 * (-t10 * t102 + t11 * t100) / 0.2e1 + t2 * t154 + (t4 * t152 + t3 * t154) * t110; m(6) * (-t13 * t112 + (t100 * t19 + t102 * t18) * t110); m(6) * (t13 ^ 2 + t18 ^ 2 + t19 ^ 2) + t158 * t20 + (t102 * t2 + t100 * t1 - t112 * (t10 * t100 + t102 * t11)) * t110;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t47(1), t47(2), t47(4), t47(7), t47(11); t47(2), t47(3), t47(5), t47(8), t47(12); t47(4), t47(5), t47(6), t47(9), t47(13); t47(7), t47(8), t47(9), t47(10), t47(14); t47(11), t47(12), t47(13), t47(14), t47(15);];
Mq = res;
