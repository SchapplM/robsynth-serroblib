% Calculate joint inertia matrix for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:52:54
% EndTime: 2019-03-09 01:52:57
% DurationCPUTime: 1.39s
% Computational Cost: add. (3363->293), mult. (3481->441), div. (0->0), fcn. (3551->10), ass. (0->143)
t125 = sin(qJ(1));
t183 = -t125 / 0.2e1;
t126 = cos(qJ(1));
t172 = t126 / 0.2e1;
t182 = rSges(6,3) + qJ(5);
t116 = pkin(9) + qJ(4);
t108 = sin(t116);
t121 = cos(pkin(10));
t152 = t125 * t121;
t119 = sin(pkin(10));
t157 = t119 * t126;
t78 = t108 * t157 + t152;
t153 = t125 * t119;
t156 = t121 * t126;
t79 = -t108 * t156 + t153;
t142 = -t79 * rSges(6,1) - t78 * rSges(6,2);
t110 = cos(t116);
t144 = t110 * t182;
t112 = t126 * qJ(2);
t124 = -pkin(7) - qJ(3);
t120 = sin(pkin(9));
t171 = pkin(3) * t120;
t145 = t125 * t124 + t126 * t171 + t112;
t161 = t108 * t126;
t99 = pkin(4) * t161;
t21 = -t125 * pkin(1) - t126 * t144 + t142 + t145 + t99;
t151 = t126 * pkin(1) + t125 * qJ(2);
t128 = -t124 * t126 + t125 * t171 + t151;
t162 = t108 * t125;
t76 = -t108 * t153 + t156;
t77 = t108 * t152 + t157;
t179 = -t77 * rSges(6,1) - t76 * rSges(6,2) - pkin(4) * t162;
t22 = -t125 * t144 + t128 - t179;
t181 = m(6) * (t125 * t21 - t126 * t22);
t123 = -pkin(8) - qJ(5);
t115 = pkin(10) + qJ(6);
t107 = sin(t115);
t109 = cos(t115);
t154 = t125 * t109;
t69 = t107 * t161 + t154;
t155 = t125 * t107;
t160 = t109 * t126;
t70 = -t108 * t160 + t155;
t141 = -t70 * rSges(7,1) - t69 * rSges(7,2);
t104 = pkin(5) * t121 + pkin(4);
t163 = t104 * t108;
t17 = (-pkin(5) * t119 - pkin(1)) * t125 + (t163 + (-rSges(7,3) + t123) * t110) * t126 + t141 + t145;
t159 = t110 * t125;
t67 = -t108 * t155 + t160;
t68 = t107 * t126 + t108 * t154;
t31 = t68 * rSges(7,1) + t67 * rSges(7,2) - rSges(7,3) * t159;
t178 = -pkin(5) * t157 - t104 * t162 - t123 * t159 - t31;
t18 = t128 - t178;
t180 = m(7) * (t125 * t17 - t126 * t18);
t177 = (rSges(5,1) * t108 + rSges(5,2) * t110) * t126;
t176 = -qJ(5) - t123;
t117 = t125 ^ 2;
t118 = t126 ^ 2;
t54 = Icges(6,6) * t108 + (Icges(6,4) * t121 - Icges(6,2) * t119) * t110;
t175 = t54 / 0.2e1;
t55 = Icges(6,5) * t108 + (Icges(6,1) * t121 - Icges(6,4) * t119) * t110;
t174 = t55 / 0.2e1;
t173 = t125 / 0.2e1;
t100 = t117 + t118;
t91 = m(5) * t100;
t49 = Icges(7,3) * t108 + (Icges(7,5) * t109 - Icges(7,6) * t107) * t110;
t51 = Icges(7,5) * t108 + (Icges(7,1) * t109 - Icges(7,4) * t107) * t110;
t170 = t110 * t109 * t51 + t108 * t49;
t52 = rSges(7,3) * t108 + (rSges(7,1) * t109 - rSges(7,2) * t107) * t110;
t169 = (-pkin(4) + t104) * t110 + t176 * t108 + t52;
t50 = Icges(7,6) * t108 + (Icges(7,4) * t109 - Icges(7,2) * t107) * t110;
t167 = t107 * t50;
t166 = rSges(4,3) + qJ(3);
t165 = Icges(5,4) * t108;
t164 = Icges(5,4) * t110;
t158 = t110 * t126;
t150 = m(6) / 0.2e1 + m(7) / 0.2e1;
t148 = rSges(5,1) * t162 + rSges(5,2) * t159 + t126 * rSges(5,3);
t25 = Icges(7,5) * t68 + Icges(7,6) * t67 - Icges(7,3) * t159;
t27 = Icges(7,4) * t68 + Icges(7,2) * t67 - Icges(7,6) * t159;
t29 = Icges(7,1) * t68 + Icges(7,4) * t67 - Icges(7,5) * t159;
t10 = t108 * t25 + (-t107 * t27 + t109 * t29) * t110;
t13 = -t49 * t159 + t50 * t67 + t51 * t68;
t147 = -t10 / 0.2e1 - t13 / 0.2e1;
t26 = Icges(7,5) * t70 + Icges(7,6) * t69 + Icges(7,3) * t158;
t28 = Icges(7,4) * t70 + Icges(7,2) * t69 + Icges(7,6) * t158;
t30 = Icges(7,1) * t70 + Icges(7,4) * t69 + Icges(7,5) * t158;
t11 = t108 * t26 + (-t107 * t28 + t109 * t30) * t110;
t14 = t49 * t158 + t69 * t50 + t70 * t51;
t146 = t14 / 0.2e1 + t11 / 0.2e1;
t143 = t91 + (m(4) + m(6) + m(7)) * t100;
t122 = cos(pkin(9));
t140 = rSges(4,1) * t120 + rSges(4,2) * t122;
t19 = t108 * t31 + t52 * t159;
t32 = rSges(7,3) * t158 - t141;
t20 = -t108 * t32 + t52 * t158;
t135 = t20 * t125 - t126 * t19;
t87 = pkin(4) * t110 + qJ(5) * t108;
t80 = t125 * t87;
t23 = t169 * t125 + t80;
t24 = (-t87 - t169) * t126;
t133 = t23 * t125 - t126 * t24;
t56 = rSges(6,3) * t108 + (rSges(6,1) * t121 - rSges(6,2) * t119) * t110;
t40 = t125 * t56 + t80;
t41 = (-t56 - t87) * t126;
t132 = t40 * t125 - t126 * t41;
t131 = Icges(5,1) * t108 + t164;
t130 = Icges(5,2) * t110 + t165;
t129 = Icges(5,5) * t108 + Icges(5,6) * t110;
t42 = t177 + (-rSges(5,3) - pkin(1)) * t125 + t145;
t43 = t128 + t148;
t127 = m(5) * (t125 * t42 - t126 * t43);
t95 = rSges(2,1) * t126 - t125 * rSges(2,2);
t94 = -t125 * rSges(2,1) - rSges(2,2) * t126;
t88 = rSges(5,1) * t110 - rSges(5,2) * t108;
t83 = Icges(5,5) * t110 - Icges(5,6) * t108;
t75 = -rSges(3,2) * t126 + t125 * rSges(3,3) + t151;
t74 = rSges(3,3) * t126 + t112 + (rSges(3,2) - pkin(1)) * t125;
t66 = t126 * (qJ(5) * t158 - t99);
t58 = Icges(5,3) * t125 - t129 * t126;
t57 = Icges(5,3) * t126 + t129 * t125;
t48 = t140 * t125 + t166 * t126 + t151;
t47 = t112 + t140 * t126 + (-pkin(1) - t166) * t125;
t39 = Icges(6,1) * t79 + Icges(6,4) * t78 + Icges(6,5) * t158;
t38 = Icges(6,1) * t77 + Icges(6,4) * t76 - Icges(6,5) * t159;
t37 = Icges(6,4) * t79 + Icges(6,2) * t78 + Icges(6,6) * t158;
t36 = Icges(6,4) * t77 + Icges(6,2) * t76 - Icges(6,6) * t159;
t35 = Icges(6,5) * t79 + Icges(6,6) * t78 + Icges(6,3) * t158;
t34 = Icges(6,5) * t77 + Icges(6,6) * t76 - Icges(6,3) * t159;
t33 = -t125 * t148 + (t125 * rSges(5,3) - t177) * t126;
t16 = (-t110 * t167 + t170) * t108;
t15 = (-t125 * t32 - t126 * t31) * t110;
t12 = t66 + t126 * (rSges(6,3) * t158 - t142) + (t159 * t182 + t179) * t125;
t9 = t26 * t158 + t69 * t28 + t70 * t30;
t8 = t25 * t158 + t69 * t27 + t70 * t29;
t7 = -t26 * t159 + t28 * t67 + t30 * t68;
t6 = -t25 * t159 + t27 * t67 + t29 * t68;
t5 = t66 + t178 * t125 + (pkin(5) * t153 + t32 + t99 + (t176 * t110 - t163) * t126) * t126;
t4 = t9 * t125 + t126 * t8;
t3 = t7 * t125 + t126 * t6;
t2 = t14 * t108 + (-t125 * t8 + t126 * t9) * t110;
t1 = t13 * t108 + (-t125 * t6 + t126 * t7) * t110;
t44 = [Icges(4,1) * t122 ^ 2 + Icges(3,1) + Icges(2,3) + (-0.2e1 * Icges(4,4) * t122 + Icges(4,2) * t120) * t120 + (-t164 + (Icges(6,5) * t121 - Icges(6,6) * t119) * t110 + (Icges(5,2) + Icges(6,3)) * t108) * t108 + (Icges(5,1) * t110 - t119 * t54 + t121 * t55 - t165 - t167) * t110 + m(7) * (t17 ^ 2 + t18 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t47 ^ 2 + t48 ^ 2) + m(3) * (t74 ^ 2 + t75 ^ 2) + m(2) * (t94 ^ 2 + t95 ^ 2) + t170; t180 + t181 + t127 + m(4) * (t125 * t47 - t126 * t48) + m(3) * (t125 * t74 - t126 * t75); m(3) * t100 + t143; m(7) * (t125 * t18 + t126 * t17) + m(6) * (t125 * t22 + t126 * t21) + m(5) * (t125 * t43 + t126 * t42) + m(4) * (t125 * t48 + t126 * t47); 0; t143; (t83 * t172 + t77 * t174 + t76 * t175 - t147) * t126 + (t83 * t173 + t79 * t174 + t78 * t175 + t146) * t125 + m(7) * (t17 * t23 + t18 * t24) + m(6) * (t21 * t40 + t22 * t41) + t88 * t127 + ((-Icges(5,6) * t126 / 0.2e1 + t130 * t183 + t34 / 0.2e1) * t126 + (Icges(5,6) * t183 + t130 * t172 + t35 / 0.2e1) * t125) * t108 + ((Icges(5,5) * t125 - t119 * t37 + t121 * t39 - t131 * t126) * t173 + (Icges(5,5) * t126 - t119 * t36 + t121 * t38 + t131 * t125) * t172) * t110; m(6) * t132 + m(7) * t133 + t88 * t91; m(6) * (t41 * t125 + t126 * t40) + m(7) * (t24 * t125 + t126 * t23); m(7) * (t23 ^ 2 + t24 ^ 2 + t5 ^ 2) + m(6) * (t12 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(5) * (t100 * t88 ^ 2 + t33 ^ 2) + (t118 * t57 + t3 + (-t34 * t159 + t76 * t36 + t77 * t38) * t126) * t126 + (t4 + t117 * t58 + (t35 * t158 + t78 * t37 + t79 * t39) * t125 + (t125 * t57 + t126 * t58 + t34 * t158 - t35 * t159 + t78 * t36 + t37 * t76 + t79 * t38 + t39 * t77) * t126) * t125; 0.2e1 * (-t180 / 0.2e1 - t181 / 0.2e1) * t110; -0.2e1 * t150 * t100 * t110; 0; m(7) * (t108 * t5 - t110 * t133) + m(6) * (t108 * t12 - t110 * t132); 0.2e1 * t150 * (t100 * t110 ^ 2 + t108 ^ 2); t16 + m(7) * (t17 * t20 + t18 * t19) + (t125 * t147 + t126 * t146) * t110; m(7) * t135; m(7) * (t125 * t19 + t126 * t20); t108 * (t10 * t126 + t11 * t125) / 0.2e1 + m(7) * (t15 * t5 + t19 * t24 + t20 * t23) + t1 * t172 + t2 * t173 + (t4 * t172 + t183 * t3) * t110; m(7) * (t108 * t15 - t110 * t135); t108 * t16 + m(7) * (t15 ^ 2 + t19 ^ 2 + t20 ^ 2) + (-t125 * t1 + t126 * t2 + t108 * (-t10 * t125 + t11 * t126)) * t110;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t44(1) t44(2) t44(4) t44(7) t44(11) t44(16); t44(2) t44(3) t44(5) t44(8) t44(12) t44(17); t44(4) t44(5) t44(6) t44(9) t44(13) t44(18); t44(7) t44(8) t44(9) t44(10) t44(14) t44(19); t44(11) t44(12) t44(13) t44(14) t44(15) t44(20); t44(16) t44(17) t44(18) t44(19) t44(20) t44(21);];
Mq  = res;
