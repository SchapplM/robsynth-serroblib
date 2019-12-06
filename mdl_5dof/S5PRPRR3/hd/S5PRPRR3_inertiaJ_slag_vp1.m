% Calculate joint inertia matrix for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:38
% EndTime: 2019-12-05 15:46:43
% DurationCPUTime: 1.45s
% Computational Cost: add. (4808->225), mult. (4994->361), div. (0->0), fcn. (5423->10), ass. (0->124)
t174 = Icges(3,3) + Icges(4,3);
t117 = sin(pkin(8));
t113 = t117 ^ 2;
t118 = cos(pkin(8));
t114 = t118 ^ 2;
t146 = t113 + t114;
t115 = qJ(2) + pkin(9);
t109 = sin(t115);
t110 = cos(t115);
t121 = sin(qJ(2));
t123 = cos(qJ(2));
t173 = Icges(3,5) * t123 + Icges(4,5) * t110 - Icges(3,6) * t121 - Icges(4,6) * t109;
t172 = -t173 * t117 + t174 * t118;
t171 = t174 * t117 + t173 * t118;
t170 = t110 ^ 2;
t169 = -t110 / 0.2e1;
t168 = t117 / 0.2e1;
t167 = -t118 / 0.2e1;
t166 = pkin(2) * t121;
t122 = cos(qJ(4));
t164 = pkin(4) * t122;
t125 = pkin(7) * t109 + t110 * t164;
t120 = sin(qJ(4));
t150 = t117 * t120;
t155 = t109 * t118;
t116 = qJ(4) + qJ(5);
t112 = cos(t116);
t152 = t112 * t117;
t111 = sin(t116);
t153 = t111 * t118;
t90 = -t110 * t153 + t152;
t151 = t112 * t118;
t154 = t111 * t117;
t91 = t110 * t151 + t154;
t52 = rSges(6,1) * t91 + rSges(6,2) * t90 + rSges(6,3) * t155;
t161 = pkin(4) * t150 + t118 * t125 + t52;
t156 = t109 * t117;
t88 = -t110 * t154 - t151;
t89 = t110 * t152 - t153;
t51 = rSges(6,1) * t89 + rSges(6,2) * t88 + rSges(6,3) * t156;
t71 = -rSges(6,3) * t110 + (rSges(6,1) * t112 - rSges(6,2) * t111) * t109;
t37 = t110 * t51 + t71 * t156;
t67 = -pkin(7) * t110 + t109 * t164;
t160 = -t67 - t71;
t159 = t146 * pkin(2) * t123;
t68 = -Icges(6,3) * t110 + (Icges(6,5) * t112 - Icges(6,6) * t111) * t109;
t158 = t110 * t68;
t74 = -Icges(5,3) * t110 + (Icges(5,5) * t122 - Icges(5,6) * t120) * t109;
t157 = t110 * t74;
t149 = t117 * t122;
t148 = t118 * t120;
t147 = t118 * t122;
t145 = -rSges(4,1) * t109 - rSges(4,2) * t110 - t166;
t144 = -t109 * pkin(3) + t110 * pkin(6) - t166;
t143 = t159 + t146 * (pkin(3) * t110 + pkin(6) * t109);
t45 = Icges(6,5) * t89 + Icges(6,6) * t88 + Icges(6,3) * t156;
t47 = Icges(6,4) * t89 + Icges(6,2) * t88 + Icges(6,6) * t156;
t49 = Icges(6,1) * t89 + Icges(6,4) * t88 + Icges(6,5) * t156;
t26 = -t110 * t45 + (-t111 * t47 + t112 * t49) * t109;
t46 = Icges(6,5) * t91 + Icges(6,6) * t90 + Icges(6,3) * t155;
t48 = Icges(6,4) * t91 + Icges(6,2) * t90 + Icges(6,6) * t155;
t50 = Icges(6,1) * t91 + Icges(6,4) * t90 + Icges(6,5) * t155;
t27 = -t110 * t46 + (-t111 * t48 + t112 * t50) * t109;
t19 = t156 * t45 + t47 * t88 + t49 * t89;
t20 = t156 * t46 + t48 * t88 + t50 * t89;
t69 = -Icges(6,6) * t110 + (Icges(6,4) * t112 - Icges(6,2) * t111) * t109;
t70 = -Icges(6,5) * t110 + (Icges(6,1) * t112 - Icges(6,4) * t111) * t109;
t5 = -(t69 * t88 + t70 * t89) * t110 + (t20 * t118 + (t19 - t158) * t117) * t109;
t21 = t155 * t45 + t47 * t90 + t49 * t91;
t22 = t155 * t46 + t48 * t90 + t50 * t91;
t6 = -(t69 * t90 + t70 * t91) * t110 + (t21 * t117 + (t22 - t158) * t118) * t109;
t142 = -t110 * (t170 * t68 + (t27 * t118 + t26 * t117 - (-t111 * t69 + t112 * t70) * t110) * t109) + t6 * t155 + t5 * t156;
t12 = t117 * t20 - t118 * t19;
t13 = t117 * t22 - t118 * t21;
t141 = t12 * t156 / 0.2e1 + t5 * t167 + t6 * t168 + t13 * t155 / 0.2e1 + (t117 * t27 - t118 * t26) * t169;
t77 = -rSges(5,3) * t110 + (rSges(5,1) * t122 - rSges(5,2) * t120) * t109;
t140 = t144 - t77;
t132 = t144 + t160;
t105 = rSges(3,1) * t121 + rSges(3,2) * t123;
t101 = t110 * t147 + t150;
t100 = -t110 * t148 + t149;
t99 = t110 * t149 - t148;
t98 = -t110 * t150 - t147;
t79 = t145 * t118;
t78 = t145 * t117;
t76 = -Icges(5,5) * t110 + (Icges(5,1) * t122 - Icges(5,4) * t120) * t109;
t75 = -Icges(5,6) * t110 + (Icges(5,4) * t122 - Icges(5,2) * t120) * t109;
t65 = t146 * (rSges(3,1) * t123 - rSges(3,2) * t121);
t64 = rSges(5,1) * t101 + rSges(5,2) * t100 + rSges(5,3) * t155;
t63 = rSges(5,1) * t99 + rSges(5,2) * t98 + rSges(5,3) * t156;
t62 = Icges(5,1) * t101 + Icges(5,4) * t100 + Icges(5,5) * t155;
t61 = Icges(5,1) * t99 + Icges(5,4) * t98 + Icges(5,5) * t156;
t60 = Icges(5,4) * t101 + Icges(5,2) * t100 + Icges(5,6) * t155;
t59 = Icges(5,4) * t99 + Icges(5,2) * t98 + Icges(5,6) * t156;
t58 = Icges(5,5) * t101 + Icges(5,6) * t100 + Icges(5,3) * t155;
t57 = Icges(5,5) * t99 + Icges(5,6) * t98 + Icges(5,3) * t156;
t55 = -pkin(4) * t148 + t117 * t125;
t54 = t140 * t118;
t53 = t140 * t117;
t43 = t51 * t155;
t42 = -t110 * t64 - t155 * t77;
t41 = t110 * t63 + t156 * t77;
t40 = t132 * t118;
t39 = t132 * t117;
t38 = -t110 * t52 - t155 * t71;
t36 = t159 + t146 * (rSges(4,1) * t110 - rSges(4,2) * t109);
t35 = (-t117 * t64 + t118 * t63) * t109;
t34 = -t156 * t52 + t43;
t33 = -t110 * t58 + (-t120 * t60 + t122 * t62) * t109;
t32 = -t110 * t57 + (-t120 * t59 + t122 * t61) * t109;
t31 = t100 * t60 + t101 * t62 + t155 * t58;
t30 = t100 * t59 + t101 * t61 + t155 * t57;
t29 = t156 * t58 + t60 * t98 + t62 * t99;
t28 = t156 * t57 + t59 * t98 + t61 * t99;
t25 = -t110 * t161 + t155 * t160;
t24 = t110 * t55 + t156 * t67 + t37;
t23 = t117 * t63 + t118 * t64 + t143;
t18 = t43 + (-t117 * t161 + t118 * t55) * t109;
t17 = t161 * t118 + (t51 + t55) * t117 + t143;
t16 = t117 * t31 - t118 * t30;
t15 = t117 * t29 - t118 * t28;
t8 = -(t100 * t75 + t101 * t76) * t110 + (t30 * t117 + (t31 - t157) * t118) * t109;
t7 = -(t75 * t98 + t76 * t99) * t110 + (t29 * t118 + (t28 - t157) * t117) * t109;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t65 + m(4) * t36 + m(5) * t23 + m(6) * t17; m(6) * (t17 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t23 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(4) * (t36 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(3) * (t105 ^ 2 * t146 + t65 ^ 2) + (t172 * t114 - t12 - t15) * t118 + (t13 + t16 + t171 * t113 + (t172 * t117 + t171 * t118) * t118) * t117; 0; m(6) * (t117 * t40 - t118 * t39) + m(5) * (t117 * t54 - t118 * t53) + m(4) * (t117 * t79 - t118 * t78); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t146; m(5) * t35 + m(6) * t18; t8 * t168 + (t117 * t33 - t118 * t32) * t169 + t7 * t167 + (t118 * t16 / 0.2e1 + t15 * t168) * t109 + m(6) * (t17 * t18 + t24 * t40 + t25 * t39) + m(5) * (t23 * t35 + t41 * t54 + t42 * t53) + t141; m(5) * (t117 * t41 - t118 * t42) + m(6) * (t117 * t24 - t118 * t25); t7 * t156 + t8 * t155 + m(6) * (t18 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t35 ^ 2 + t41 ^ 2 + t42 ^ 2) - t110 * (t170 * t74 + (t33 * t118 + t32 * t117 - (-t120 * t75 + t122 * t76) * t110) * t109) + t142; m(6) * t34; m(6) * (t17 * t34 + t37 * t40 + t38 * t39) + t141; m(6) * (t117 * t37 - t118 * t38); m(6) * (t18 * t34 + t24 * t37 + t25 * t38) + t142; m(6) * (t34 ^ 2 + t37 ^ 2 + t38 ^ 2) + t142;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
