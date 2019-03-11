% Calculate joint inertia matrix for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:32
% EndTime: 2019-03-09 01:48:35
% DurationCPUTime: 1.67s
% Computational Cost: add. (2187->288), mult. (3427->429), div. (0->0), fcn. (3497->8), ass. (0->131)
t118 = cos(qJ(4));
t119 = cos(qJ(1));
t141 = t118 * t119;
t116 = sin(qJ(4));
t147 = t116 * t119;
t113 = sin(pkin(9));
t114 = cos(pkin(9));
t117 = sin(qJ(1));
t143 = t117 * t114;
t76 = -t113 * t147 - t143;
t144 = t117 * t113;
t77 = t114 * t147 - t144;
t168 = -t77 * rSges(6,1) - t76 * rSges(6,2) - pkin(4) * t147 + (rSges(6,3) + qJ(5)) * t141;
t161 = -t117 / 0.2e1;
t167 = -t119 / 0.2e1;
t160 = t119 / 0.2e1;
t107 = t119 * qJ(2);
t74 = t114 * t119 - t116 * t144;
t75 = t113 * t119 + t116 * t143;
t133 = -t75 * rSges(6,1) - t74 * rSges(6,2);
t157 = -pkin(1) - qJ(3);
t158 = pkin(4) * t116;
t142 = t117 * t118;
t98 = qJ(5) * t142;
t21 = -pkin(7) * t119 + t107 + t98 + (rSges(6,3) * t118 + t157 - t158) * t117 + t133;
t140 = t119 * pkin(1) + t117 * qJ(2);
t136 = t119 * qJ(3) + t140;
t22 = -t117 * pkin(7) + t136 - t168;
t165 = m(6) * (t117 * t22 + t119 * t21);
t115 = -pkin(8) - qJ(5);
t109 = pkin(9) + qJ(6);
t104 = cos(t109);
t103 = sin(t109);
t146 = t117 * t103;
t59 = t104 * t119 - t116 * t146;
t145 = t117 * t104;
t60 = t103 * t119 + t116 * t145;
t132 = -t60 * rSges(7,1) - t59 * rSges(7,2);
t135 = -pkin(5) * t113 - pkin(7);
t102 = pkin(5) * t114 + pkin(4);
t148 = t102 * t116;
t16 = t107 + t135 * t119 + (-t148 + (rSges(7,3) - t115) * t118 + t157) * t117 + t132;
t61 = -t103 * t147 - t145;
t62 = t104 * t147 - t146;
t32 = t62 * rSges(7,1) + t61 * rSges(7,2) - rSges(7,3) * t141;
t163 = -t102 * t147 - t115 * t141 - t32;
t17 = t135 * t117 + t136 - t163;
t164 = m(7) * (t117 * t17 + t119 * t16);
t110 = t117 ^ 2;
t112 = t119 ^ 2;
t97 = t110 + t112;
t85 = m(5) * t97;
t159 = -rSges(5,3) - pkin(7);
t47 = Icges(7,3) * t116 + (Icges(7,5) * t104 - Icges(7,6) * t103) * t118;
t49 = Icges(7,5) * t116 + (Icges(7,1) * t104 - Icges(7,4) * t103) * t118;
t156 = t118 * t104 * t49 + t116 * t47;
t50 = t116 * rSges(7,3) + (rSges(7,1) * t104 - rSges(7,2) * t103) * t118;
t155 = (-pkin(4) + t102) * t118 + (-qJ(5) - t115) * t116 + t50;
t48 = Icges(7,6) * t116 + (Icges(7,4) * t104 - Icges(7,2) * t103) * t118;
t152 = t103 * t48;
t151 = rSges(5,1) * t147 + rSges(5,2) * t141;
t150 = Icges(5,4) * t116;
t149 = Icges(5,4) * t118;
t139 = m(6) / 0.2e1 + m(7) / 0.2e1;
t25 = Icges(7,5) * t60 + Icges(7,6) * t59 - Icges(7,3) * t142;
t27 = Icges(7,4) * t60 + Icges(7,2) * t59 - Icges(7,6) * t142;
t29 = Icges(7,1) * t60 + Icges(7,4) * t59 - Icges(7,5) * t142;
t10 = t116 * t25 + (-t103 * t27 + t104 * t29) * t118;
t13 = -t47 * t142 + t48 * t59 + t49 * t60;
t138 = -t10 / 0.2e1 - t13 / 0.2e1;
t26 = Icges(7,5) * t62 + Icges(7,6) * t61 - Icges(7,3) * t141;
t28 = Icges(7,4) * t62 + Icges(7,2) * t61 - Icges(7,6) * t141;
t30 = Icges(7,1) * t62 + Icges(7,4) * t61 - Icges(7,5) * t141;
t11 = t116 * t26 + (-t103 * t28 + t104 * t30) * t118;
t14 = -t47 * t141 + t61 * t48 + t62 * t49;
t137 = -t11 / 0.2e1 - t14 / 0.2e1;
t134 = t85 + (m(4) + m(6) + m(7)) * t97;
t131 = -rSges(5,1) * t116 - rSges(5,2) * t118;
t31 = -rSges(7,3) * t142 - t132;
t19 = -t116 * t31 - t50 * t142;
t20 = t116 * t32 + t50 * t141;
t127 = t20 * t117 + t119 * t19;
t93 = pkin(4) * t118 + qJ(5) * t116;
t80 = t117 * t93;
t23 = t155 * t117 + t80;
t81 = t119 * t93;
t24 = t155 * t119 + t81;
t125 = t23 * t117 + t24 * t119;
t58 = t116 * rSges(6,3) + (rSges(6,1) * t114 - rSges(6,2) * t113) * t118;
t42 = t117 * t58 + t80;
t43 = t119 * t58 + t81;
t124 = t42 * t117 + t119 * t43;
t123 = Icges(5,1) * t116 + t149;
t122 = Icges(5,2) * t118 + t150;
t121 = Icges(5,5) * t116 + Icges(5,6) * t118;
t40 = t107 + t159 * t119 + (t131 + t157) * t117;
t41 = t159 * t117 + t136 + t151;
t120 = m(5) * (t117 * t41 + t119 * t40);
t95 = rSges(2,1) * t119 - t117 * rSges(2,2);
t94 = rSges(5,1) * t118 - rSges(5,2) * t116;
t92 = -t117 * rSges(2,1) - rSges(2,2) * t119;
t89 = Icges(5,5) * t118 - Icges(5,6) * t116;
t78 = t117 * t158 - t98;
t70 = -rSges(3,2) * t119 + t117 * rSges(3,3) + t140;
t69 = rSges(3,3) * t119 + t107 + (rSges(3,2) - pkin(1)) * t117;
t64 = -Icges(5,3) * t117 + t121 * t119;
t63 = Icges(5,3) * t119 + t121 * t117;
t57 = Icges(6,5) * t116 + (Icges(6,1) * t114 - Icges(6,4) * t113) * t118;
t56 = Icges(6,6) * t116 + (Icges(6,4) * t114 - Icges(6,2) * t113) * t118;
t52 = t117 * rSges(4,2) + rSges(4,3) * t119 + t136;
t51 = rSges(4,2) * t119 + t107 + (-rSges(4,3) + t157) * t117;
t39 = Icges(6,1) * t77 + Icges(6,4) * t76 - Icges(6,5) * t141;
t38 = Icges(6,1) * t75 + Icges(6,4) * t74 - Icges(6,5) * t142;
t37 = Icges(6,4) * t77 + Icges(6,2) * t76 - Icges(6,6) * t141;
t36 = Icges(6,4) * t75 + Icges(6,2) * t74 - Icges(6,6) * t142;
t35 = Icges(6,5) * t77 + Icges(6,6) * t76 - Icges(6,3) * t141;
t34 = Icges(6,5) * t75 + Icges(6,6) * t74 - Icges(6,3) * t142;
t33 = t131 * t110 - t119 * t151;
t18 = (-t118 * t152 + t156) * t116;
t15 = (t117 * t32 - t119 * t31) * t118;
t12 = t168 * t119 + (rSges(6,3) * t142 + t133 - t78) * t117;
t9 = -t26 * t141 + t61 * t28 + t62 * t30;
t8 = -t25 * t141 + t61 * t27 + t62 * t29;
t7 = -t26 * t142 + t28 * t59 + t30 * t60;
t6 = -t25 * t142 + t27 * t59 + t29 * t60;
t5 = t163 * t119 + (-t31 - t78 - t98 + (-t115 * t118 - t148 + t158) * t117) * t117;
t4 = -t9 * t117 + t119 * t8;
t3 = -t7 * t117 + t119 * t6;
t2 = t14 * t116 + (-t117 * t8 - t119 * t9) * t118;
t1 = t13 * t116 + (-t117 * t6 - t119 * t7) * t118;
t44 = [Icges(3,1) + Icges(4,1) + Icges(2,3) + (-t149 + (Icges(6,5) * t114 - Icges(6,6) * t113) * t118 + (Icges(5,2) + Icges(6,3)) * t116) * t116 + (Icges(5,1) * t118 - t113 * t56 + t114 * t57 - t150 - t152) * t118 + m(7) * (t16 ^ 2 + t17 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t40 ^ 2 + t41 ^ 2) + m(4) * (t51 ^ 2 + t52 ^ 2) + m(3) * (t69 ^ 2 + t70 ^ 2) + m(2) * (t92 ^ 2 + t95 ^ 2) + t156; m(7) * (t117 * t16 - t119 * t17) + m(6) * (t117 * t21 - t119 * t22) + m(5) * (t117 * t40 - t119 * t41) + m(4) * (t117 * t51 - t119 * t52) + m(3) * (t117 * t69 - t119 * t70); m(3) * t97 + t134; t164 + t165 + t120 + m(4) * (t117 * t52 + t119 * t51); 0; t134; (t74 * t56 / 0.2e1 + t75 * t57 / 0.2e1 + t89 * t160 - t138) * t119 + (-t76 * t56 / 0.2e1 - t77 * t57 / 0.2e1 + t117 * t89 / 0.2e1 + t137) * t117 + m(7) * (t16 * t24 + t17 * t23) + m(6) * (t21 * t43 + t22 * t42) + t94 * t120 + ((Icges(5,6) * t167 + t122 * t161 + t34 / 0.2e1) * t119 + (Icges(5,6) * t161 + t122 * t160 - t35 / 0.2e1) * t117) * t116 + ((-Icges(5,5) * t117 - t113 * t37 + t114 * t39 + t123 * t119) * t161 + (Icges(5,5) * t119 - t113 * t36 + t114 * t38 + t123 * t117) * t160) * t118; m(6) * (t43 * t117 - t119 * t42) + m(7) * (t24 * t117 - t119 * t23); m(6) * t124 + m(7) * t125 + t94 * t85; m(7) * (t23 ^ 2 + t24 ^ 2 + t5 ^ 2) + m(6) * (t12 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t97 * t94 ^ 2 + t33 ^ 2) + (t112 * t63 + t3 + (-t34 * t142 + t74 * t36 + t75 * t38) * t119) * t119 + (-t4 - t110 * t64 + (-t35 * t141 + t76 * t37 + t77 * t39) * t117 + (t117 * t63 - t119 * t64 + t34 * t141 + t35 * t142 - t76 * t36 - t37 * t74 - t77 * t38 - t39 * t75) * t119) * t117; 0.2e1 * (-t164 / 0.2e1 - t165 / 0.2e1) * t118; 0; -0.2e1 * t139 * t97 * t118; m(7) * (t116 * t5 - t125 * t118) + m(6) * (t116 * t12 - t124 * t118); 0.2e1 * t139 * (t97 * t118 ^ 2 + t116 ^ 2); t18 + m(7) * (t16 * t19 + t17 * t20) + (t138 * t117 + t137 * t119) * t118; m(7) * (t19 * t117 - t119 * t20); m(7) * t127; t2 * t161 + t1 * t160 + m(7) * (t15 * t5 + t19 * t24 + t20 * t23) + t116 * (t10 * t119 - t11 * t117) / 0.2e1 + (t3 * t161 + t4 * t167) * t118; m(7) * (t15 * t116 - t127 * t118); t116 * t18 + m(7) * (t15 ^ 2 + t19 ^ 2 + t20 ^ 2) + (-t119 * t2 - t117 * t1 + t116 * (-t10 * t117 - t11 * t119)) * t118;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t44(1) t44(2) t44(4) t44(7) t44(11) t44(16); t44(2) t44(3) t44(5) t44(8) t44(12) t44(17); t44(4) t44(5) t44(6) t44(9) t44(13) t44(18); t44(7) t44(8) t44(9) t44(10) t44(14) t44(19); t44(11) t44(12) t44(13) t44(14) t44(15) t44(20); t44(16) t44(17) t44(18) t44(19) t44(20) t44(21);];
Mq  = res;
