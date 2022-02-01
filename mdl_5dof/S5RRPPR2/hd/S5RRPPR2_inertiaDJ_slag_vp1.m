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
% m [6x1]
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:05:27
% EndTime: 2022-01-20 10:05:33
% DurationCPUTime: 1.68s
% Computational Cost: add. (5693->268), mult. (4836->389), div. (0->0), fcn. (4330->10), ass. (0->138)
t124 = qJ(1) + qJ(2);
t120 = sin(t124);
t121 = cos(t124);
t97 = t121 * rSges(3,1) - rSges(3,2) * t120;
t119 = pkin(8) + t124;
t116 = sin(t119);
t118 = pkin(2) * t121;
t176 = t116 * qJ(4) + t118;
t117 = cos(t119);
t175 = -t117 * rSges(4,1) - t118;
t127 = sin(qJ(5));
t129 = cos(qJ(5));
t125 = sin(pkin(9));
t126 = cos(pkin(9));
t158 = Icges(6,4) * t129;
t83 = -Icges(6,6) * t126 + (-Icges(6,2) * t127 + t158) * t125;
t159 = Icges(6,4) * t127;
t84 = -Icges(6,5) * t126 + (Icges(6,1) * t129 - t159) * t125;
t174 = -t127 * t84 - t129 * t83;
t173 = 2 * m(3);
t172 = 2 * m(4);
t171 = 2 * m(5);
t170 = 2 * m(6);
t128 = sin(qJ(1));
t169 = pkin(1) * t128;
t168 = pkin(2) * t120;
t123 = qJD(1) + qJD(2);
t151 = t126 * t127;
t135 = t116 * t151 + t117 * t129;
t150 = t126 * t129;
t81 = t116 * t127 + t117 * t150;
t59 = -t81 * qJD(5) + t135 * t123;
t79 = t116 * t150 - t117 * t127;
t80 = t116 * t129 - t117 * t151;
t60 = t80 * qJD(5) - t79 * t123;
t167 = t60 * rSges(6,1) + t59 * rSges(6,2);
t165 = pkin(1) * qJD(1);
t149 = qJD(5) * t125;
t89 = (-Icges(6,2) * t129 - t159) * t149;
t163 = t127 * t89;
t88 = (-Icges(6,5) * t127 - Icges(6,6) * t129) * t149;
t90 = (-Icges(6,1) * t127 - t158) * t149;
t140 = t125 * t129 * t90 - t126 * t88;
t26 = (t174 * qJD(5) - t163) * t125 + t140;
t161 = t26 * t126;
t153 = t123 * t117;
t160 = qJ(4) * t153 + qJD(4) * t116;
t157 = t116 * t125;
t156 = t117 * t125;
t155 = t117 * t126;
t154 = t123 * t116;
t152 = t123 * t125;
t54 = t81 * rSges(6,1) + t80 * rSges(6,2) + rSges(6,3) * t156;
t148 = t128 * t165;
t130 = cos(qJ(1));
t147 = t130 * t165;
t146 = t116 * t152;
t145 = t117 * t152;
t144 = t117 * pkin(3) + t176;
t141 = -rSges(5,1) * t126 - pkin(3);
t86 = t97 * t123;
t77 = -rSges(4,2) * t116 - t175;
t61 = -t79 * qJD(5) + t80 * t123;
t62 = -t135 * qJD(5) + t81 * t123;
t137 = -t62 * rSges(6,1) - t61 * rSges(6,2);
t136 = -rSges(6,1) * t79 + rSges(6,2) * t135;
t96 = -rSges(3,1) * t120 - rSges(3,2) * t121;
t85 = t96 * t123;
t134 = -pkin(4) * t126 - pkin(3) + (-rSges(6,3) - pkin(7)) * t125;
t76 = -rSges(4,1) * t116 - rSges(4,2) * t117 - t168;
t133 = t141 * t116 - t168;
t68 = rSges(4,2) * t154 + t175 * t123;
t40 = pkin(4) * t155 + pkin(7) * t156 + t144 + t54;
t64 = rSges(5,1) * t155 - rSges(5,2) * t156 + t116 * rSges(5,3) + t144;
t67 = t76 * t123;
t47 = Icges(6,5) * t79 - Icges(6,6) * t135 + Icges(6,3) * t157;
t49 = Icges(6,4) * t79 - Icges(6,2) * t135 + Icges(6,6) * t157;
t51 = Icges(6,1) * t79 - Icges(6,4) * t135 + Icges(6,5) * t157;
t18 = -t126 * t47 + (-t127 * t49 + t129 * t51) * t125;
t48 = Icges(6,5) * t81 + Icges(6,6) * t80 + Icges(6,3) * t156;
t50 = Icges(6,4) * t81 + Icges(6,2) * t80 + Icges(6,6) * t156;
t52 = Icges(6,1) * t81 + Icges(6,4) * t80 + Icges(6,5) * t156;
t19 = -t126 * t48 + (-t127 * t50 + t129 * t52) * t125;
t28 = Icges(6,5) * t62 + Icges(6,6) * t61 + Icges(6,3) * t145;
t30 = Icges(6,4) * t62 + Icges(6,2) * t61 + Icges(6,6) * t145;
t32 = Icges(6,1) * t62 + Icges(6,4) * t61 + Icges(6,5) * t145;
t3 = -t126 * t28 + (-t127 * t30 + t129 * t32 + (-t127 * t51 - t129 * t49) * qJD(5)) * t125;
t82 = -Icges(6,3) * t126 + (Icges(6,5) * t129 - Icges(6,6) * t127) * t125;
t35 = -t135 * t83 + t82 * t157 + t79 * t84;
t36 = t82 * t156 + t80 * t83 + t81 * t84;
t27 = Icges(6,5) * t60 + Icges(6,6) * t59 - Icges(6,3) * t146;
t29 = Icges(6,4) * t60 + Icges(6,2) * t59 - Icges(6,6) * t146;
t31 = Icges(6,1) * t60 + Icges(6,4) * t59 - Icges(6,5) * t146;
t4 = -t126 * t27 + (-t127 * t29 + t129 * t31 + (-t127 * t52 - t129 * t50) * qJD(5)) * t125;
t8 = t59 * t83 + t60 * t84 + t80 * t89 + t81 * t90 + (t117 * t88 - t82 * t154) * t125;
t9 = t61 * t83 + t62 * t84 - t135 * t89 + t79 * t90 + (t116 * t88 + t82 * t153) * t125;
t132 = -t161 + (t3 + t9) * t157 / 0.2e1 - (t36 + t19) * t146 / 0.2e1 + (t4 + t8 + (t18 + t35) * t123) * t156 / 0.2e1;
t110 = t117 * qJ(4);
t63 = rSges(5,2) * t157 + t117 * rSges(5,3) + t110 + t133;
t131 = t134 * t116 - t168;
t45 = rSges(5,2) * t146 + rSges(5,3) * t153 + t133 * t123 + t160;
t108 = qJD(4) * t117;
t46 = t108 + rSges(5,2) * t145 + (-t118 + t141 * t117 + (-rSges(5,3) - qJ(4)) * t116) * t123;
t39 = t110 + t131 + t136;
t20 = t131 * t123 + t160 + t167;
t21 = t108 + (t134 * t117 - t176) * t123 + t137;
t122 = t130 * pkin(1);
t93 = t122 + t97;
t92 = t96 - t169;
t91 = (-rSges(6,1) * t127 - rSges(6,2) * t129) * t149;
t87 = -rSges(6,3) * t126 + (rSges(6,1) * t129 - rSges(6,2) * t127) * t125;
t72 = -t86 - t147;
t71 = t85 - t148;
t70 = t122 + t77;
t69 = t76 - t169;
t66 = t68 - t147;
t65 = t67 - t148;
t56 = t122 + t64;
t55 = t63 - t169;
t53 = rSges(6,3) * t157 - t136;
t44 = t46 - t147;
t43 = t45 - t148;
t42 = -t126 * t54 - t87 * t156;
t41 = t126 * t53 + t87 * t157;
t38 = t122 + t40;
t37 = t39 - t169;
t34 = rSges(6,3) * t145 - t137;
t33 = -rSges(6,3) * t146 + t167;
t23 = t126 * t34 + (t116 * t91 + t87 * t153) * t125;
t22 = -t126 * t33 + (-t117 * t91 + t87 * t154) * t125;
t17 = t21 - t147;
t16 = t20 - t148;
t13 = t48 * t156 + t50 * t80 + t52 * t81;
t12 = t47 * t156 + t49 * t80 + t51 * t81;
t11 = -t135 * t50 + t48 * t157 + t52 * t79;
t10 = -t135 * t49 + t47 * t157 + t51 * t79;
t5 = ((-t123 * t54 + t34) * t117 + (-t123 * t53 - t33) * t116) * t125;
t1 = [(t16 * t38 + t17 * t37) * t170 - t125 * t163 + (t65 * t70 + t66 * t69) * t172 + (t43 * t56 + t44 * t55) * t171 + (t71 * t93 + t72 * t92) * t173 + t140 + t174 * t149; m(6) * (t16 * t40 + t17 * t39 + t20 * t38 + t21 * t37) + m(4) * (t65 * t77 + t66 * t76 + t67 * t70 + t68 * t69) + m(5) * (t43 * t64 + t44 * t63 + t45 * t56 + t46 * t55) + m(3) * (t71 * t97 + t72 * t96 + t85 * t93 - t86 * t92) + t26; (t20 * t40 + t21 * t39) * t170 + (t67 * t77 + t68 * t76) * t172 + (t45 * t64 + t46 * t63) * t171 + (t85 * t97 - t86 * t96) * t173 + t26; 0; 0; 0; m(6) * ((t123 * t37 - t16) * t117 + (t123 * t38 + t17) * t116) + m(5) * ((t123 * t55 - t43) * t117 + (t123 * t56 + t44) * t116); m(6) * ((t123 * t39 - t20) * t117 + (t123 * t40 + t21) * t116) + m(5) * ((t123 * t63 - t45) * t117 + (t123 * t64 + t46) * t116); 0; 0; m(6) * (t16 * t42 + t17 * t41 + t22 * t38 + t23 * t37) + t132; m(6) * (t20 * t42 + t21 * t41 + t22 * t40 + t23 * t39) + t132; m(6) * t5; m(6) * ((t123 * t41 - t22) * t117 + (t123 * t42 + t23) * t116); (t42 * t22 + t41 * t23 + (-t116 * t54 + t117 * t53) * t5 * t125) * t170 - (-t36 * t126 + (t116 * t12 + t117 * t13) * t125) * t146 + (-t8 * t126 + ((t80 * t29 + t81 * t31 + t59 * t50 + t60 * t52) * t117 - t13 * t154 + (t80 * t30 + t81 * t32 + t59 * t49 + t60 * t51) * t116 + t12 * t153 + ((t117 * t27 - t48 * t154) * t117 + (t117 * t28 - t47 * t154) * t116) * t125) * t125) * t156 + (-t126 * t35 + (t10 * t116 + t11 * t117) * t125) * t145 + (-t9 * t126 + ((-t135 * t29 + t79 * t31 + t61 * t50 + t62 * t52) * t117 - t11 * t154 + (-t135 * t30 + t79 * t32 + t61 * t49 + t62 * t51) * t116 + t10 * t153 + ((t116 * t27 + t48 * t153) * t117 + (t116 * t28 + t47 * t153) * t116) * t125) * t125) * t157 - t126 * (-t161 + ((t123 * t18 + t4) * t117 + (-t123 * t19 + t3) * t116) * t125);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
