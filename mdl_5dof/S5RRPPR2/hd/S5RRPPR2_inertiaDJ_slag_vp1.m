% Calculate time derivative of joint inertia matrix for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:04
% EndTime: 2019-12-05 18:20:10
% DurationCPUTime: 1.76s
% Computational Cost: add. (5693->270), mult. (4836->393), div. (0->0), fcn. (4330->10), ass. (0->141)
t126 = sin(qJ(5));
t128 = cos(qJ(5));
t124 = sin(pkin(9));
t125 = cos(pkin(9));
t161 = Icges(6,4) * t128;
t83 = -Icges(6,6) * t125 + (-Icges(6,2) * t126 + t161) * t124;
t162 = Icges(6,4) * t126;
t84 = -Icges(6,5) * t125 + (Icges(6,1) * t128 - t162) * t124;
t179 = -t126 * t84 - t128 * t83;
t178 = 2 * m(3);
t177 = 2 * m(4);
t176 = 2 * m(5);
t175 = 2 * m(6);
t127 = sin(qJ(1));
t174 = pkin(1) * t127;
t129 = cos(qJ(1));
t173 = pkin(1) * t129;
t123 = qJ(1) + qJ(2);
t120 = sin(t123);
t172 = pkin(2) * t120;
t121 = cos(t123);
t171 = pkin(2) * t121;
t122 = qJD(1) + qJD(2);
t119 = pkin(8) + t123;
t116 = sin(t119);
t117 = cos(t119);
t152 = t125 * t128;
t134 = t116 * t152 - t117 * t126;
t153 = t125 * t126;
t135 = -t116 * t128 + t117 * t153;
t61 = t134 * qJD(5) + t135 * t122;
t78 = t116 * t153 + t117 * t128;
t81 = t116 * t126 + t117 * t152;
t62 = t78 * qJD(5) - t81 * t122;
t170 = t62 * rSges(6,1) + t61 * rSges(6,2);
t169 = -rSges(6,1) * t134 + t78 * rSges(6,2);
t168 = pkin(1) * qJD(1);
t151 = qJD(5) * t124;
t89 = (-Icges(6,2) * t128 - t162) * t151;
t166 = t126 * t89;
t88 = (-Icges(6,5) * t126 - Icges(6,6) * t128) * t151;
t90 = (-Icges(6,1) * t126 - t161) * t151;
t141 = t124 * t128 * t90 - t125 * t88;
t26 = (t179 * qJD(5) - t166) * t124 + t141;
t164 = t26 * t125;
t163 = -rSges(5,3) - qJ(4);
t160 = t116 * t124;
t159 = t117 * t124;
t158 = t120 * t122;
t157 = t121 * t122;
t156 = t122 * t116;
t155 = t122 * t117;
t154 = t122 * t124;
t85 = rSges(3,1) * t158 + rSges(3,2) * t157;
t150 = t129 * t168;
t105 = rSges(5,2) * t160;
t149 = t116 * t154;
t148 = t125 * t156;
t147 = t117 * t154;
t110 = pkin(2) * t158;
t67 = rSges(4,1) * t156 + rSges(4,2) * t155 + t110;
t143 = -rSges(5,1) * t125 - pkin(3);
t142 = t117 * qJ(4) - t172;
t100 = -t121 * rSges(3,1) + t120 * rSges(3,2);
t86 = -rSges(3,1) * t157 + rSges(3,2) * t158;
t138 = -t81 * rSges(6,1) + rSges(6,2) * t135;
t137 = -t117 * rSges(4,1) - t171;
t136 = pkin(3) * t156 - qJD(4) * t116 + t110;
t99 = -rSges(3,1) * t120 - rSges(3,2) * t121;
t59 = -t81 * qJD(5) + t78 * t122;
t60 = -t135 * qJD(5) - t134 * t122;
t33 = rSges(6,1) * t60 + t59 * rSges(6,2) - rSges(6,3) * t149;
t77 = t116 * rSges(4,2) + t137;
t133 = -pkin(4) * t125 - pkin(3) + (-rSges(6,3) - pkin(7)) * t124;
t76 = -rSges(4,1) * t116 - rSges(4,2) * t117 - t172;
t68 = rSges(4,2) * t156 + t137 * t122;
t47 = -Icges(6,5) * t134 + Icges(6,6) * t78 - Icges(6,3) * t160;
t49 = -Icges(6,4) * t134 + Icges(6,2) * t78 - Icges(6,6) * t160;
t51 = -Icges(6,1) * t134 + Icges(6,4) * t78 - Icges(6,5) * t160;
t18 = -t125 * t47 + (-t126 * t49 + t128 * t51) * t124;
t48 = Icges(6,5) * t81 - Icges(6,6) * t135 + Icges(6,3) * t159;
t50 = Icges(6,4) * t81 - Icges(6,2) * t135 + Icges(6,6) * t159;
t52 = Icges(6,1) * t81 - Icges(6,4) * t135 + Icges(6,5) * t159;
t19 = -t125 * t48 + (-t126 * t50 + t128 * t52) * t124;
t28 = Icges(6,5) * t62 + Icges(6,6) * t61 - Icges(6,3) * t147;
t30 = Icges(6,4) * t62 + Icges(6,2) * t61 - Icges(6,6) * t147;
t32 = Icges(6,1) * t62 + Icges(6,4) * t61 - Icges(6,5) * t147;
t3 = -t125 * t28 + (-t126 * t30 + t128 * t32 + (-t126 * t51 - t128 * t49) * qJD(5)) * t124;
t82 = -Icges(6,3) * t125 + (Icges(6,5) * t128 - Icges(6,6) * t126) * t124;
t35 = -t134 * t84 - t82 * t160 + t78 * t83;
t36 = -t135 * t83 + t82 * t159 + t81 * t84;
t27 = Icges(6,5) * t60 + Icges(6,6) * t59 - Icges(6,3) * t149;
t29 = Icges(6,4) * t60 + Icges(6,2) * t59 - Icges(6,6) * t149;
t31 = Icges(6,1) * t60 + Icges(6,4) * t59 - Icges(6,5) * t149;
t4 = -t125 * t27 + (-t126 * t29 + t128 * t31 + (-t126 * t52 - t128 * t50) * qJD(5)) * t124;
t8 = t59 * t83 + t60 * t84 - t135 * t89 + t81 * t90 + (t117 * t88 - t156 * t82) * t124;
t9 = t61 * t83 + t62 * t84 + t78 * t89 - t134 * t90 + (-t116 * t88 - t155 * t82) * t124;
t132 = -t164 - (t3 + t9) * t160 / 0.2e1 + (t4 + t8) * t159 / 0.2e1 - ((t19 + t36) * t116 + (t18 + t35) * t117) * t154 / 0.2e1;
t63 = t117 * rSges(5,3) + t143 * t116 + t105 + t142;
t131 = t163 * t116 + t143 * t117 - t171;
t64 = rSges(5,2) * t159 + t131;
t130 = -t116 * qJ(4) + t133 * t117 - t171;
t39 = t133 * t116 + t142 + t169;
t45 = rSges(5,1) * t148 + (t163 * t117 - t105) * t122 + t136;
t20 = pkin(4) * t148 + pkin(7) * t149 - qJ(4) * t155 + t136 - t33;
t111 = qJD(4) * t117;
t46 = rSges(5,2) * t147 + t131 * t122 + t111;
t40 = t130 + t138;
t21 = t130 * t122 + t111 + t170;
t118 = t127 * t168;
t93 = t100 - t173;
t92 = t99 - t174;
t91 = (-rSges(6,1) * t126 - rSges(6,2) * t128) * t151;
t87 = -rSges(6,3) * t125 + (rSges(6,1) * t128 - rSges(6,2) * t126) * t124;
t72 = t86 - t150;
t71 = t118 + t85;
t70 = t77 - t173;
t69 = t76 - t174;
t66 = t68 - t150;
t65 = t118 + t67;
t56 = t64 - t173;
t55 = t63 - t174;
t54 = rSges(6,3) * t159 - t138;
t53 = -rSges(6,3) * t160 + t169;
t44 = t46 - t150;
t43 = t118 + t45;
t42 = t125 * t54 + t87 * t159;
t41 = -t125 * t53 + t87 * t160;
t38 = t40 - t173;
t37 = t39 - t174;
t34 = -rSges(6,3) * t147 + t170;
t23 = -t125 * t34 + (t116 * t91 + t155 * t87) * t124;
t22 = t125 * t33 + (t117 * t91 - t156 * t87) * t124;
t17 = t21 - t150;
t16 = t118 + t20;
t13 = -t135 * t50 + t159 * t48 + t52 * t81;
t12 = -t135 * t49 + t159 * t47 + t51 * t81;
t11 = -t134 * t52 - t160 * t48 + t50 * t78;
t10 = -t134 * t51 - t160 * t47 + t49 * t78;
t5 = ((-t122 * t54 - t34) * t117 + (t122 * t53 - t33) * t116) * t124;
t1 = [(t71 * t93 + t72 * t92) * t178 + (t65 * t70 + t66 * t69) * t177 + (t43 * t56 + t44 * t55) * t176 - t124 * t166 + (t16 * t38 + t17 * t37) * t175 + t141 + t179 * t151; m(3) * (t100 * t71 + t72 * t99 + t85 * t93 + t86 * t92) + m(4) * (t65 * t77 + t66 * t76 + t67 * t70 + t68 * t69) + m(5) * (t43 * t64 + t44 * t63 + t45 * t56 + t46 * t55) + m(6) * (t16 * t40 + t17 * t39 + t20 * t38 + t21 * t37) + t26; (t20 * t40 + t21 * t39) * t175 + (t67 * t77 + t68 * t76) * t177 + (t45 * t64 + t46 * t63) * t176 + (t100 * t85 + t86 * t99) * t178 + t26; 0; 0; 0; m(5) * ((t122 * t55 + t43) * t117 + (-t122 * t56 + t44) * t116) + m(6) * ((t122 * t37 + t16) * t117 + (-t122 * t38 + t17) * t116); m(6) * ((t122 * t39 + t20) * t117 + (-t122 * t40 + t21) * t116) + m(5) * ((t122 * t63 + t45) * t117 + (-t122 * t64 + t46) * t116); 0; 0; m(6) * (t16 * t42 + t17 * t41 + t22 * t38 + t23 * t37) + t132; m(6) * (t20 * t42 + t21 * t41 + t22 * t40 + t23 * t39) + t132; m(6) * t5; m(6) * ((t122 * t41 + t22) * t117 + (-t122 * t42 + t23) * t116); (t42 * t22 + t41 * t23 + (-t116 * t54 - t117 * t53) * t5 * t124) * t175 - t125 * (-t164 + ((-t122 * t18 + t4) * t117 + (-t122 * t19 - t3) * t116) * t124) - (-t125 * t35 + (-t10 * t116 + t11 * t117) * t124) * t147 - (-t9 * t125 + (-(-t134 * t32 + t78 * t30 + t61 * t49 + t62 * t51) * t116 - t10 * t155 + (-t134 * t31 + t78 * t29 + t61 * t50 + t62 * t52) * t117 - t11 * t156 + (-(-t116 * t28 - t155 * t47) * t116 + (-t116 * t27 - t155 * t48) * t117) * t124) * t124) * t160 - (-t125 * t36 + (-t116 * t12 + t117 * t13) * t124) * t149 + (-t8 * t125 + (-(-t135 * t30 + t81 * t32 + t59 * t49 + t60 * t51) * t116 - t12 * t155 + (-t135 * t29 + t81 * t31 + t59 * t50 + t60 * t52) * t117 - t13 * t156 + (-(t117 * t28 - t156 * t47) * t116 + (t117 * t27 - t156 * t48) * t117) * t124) * t124) * t159;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
