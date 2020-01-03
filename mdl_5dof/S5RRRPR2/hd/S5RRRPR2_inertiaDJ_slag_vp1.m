% Calculate time derivative of joint inertia matrix for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:22
% EndTime: 2020-01-03 12:07:26
% DurationCPUTime: 1.72s
% Computational Cost: add. (5296->221), mult. (3006->300), div. (0->0), fcn. (2068->10), ass. (0->130)
t128 = sin(qJ(5));
t160 = qJD(5) * t128;
t130 = cos(qJ(5));
t159 = qJD(5) * t130;
t126 = qJD(1) + qJD(2);
t120 = qJD(3) + t126;
t142 = Icges(6,5) * t130 - Icges(6,6) * t128;
t170 = Icges(6,4) * t130;
t143 = -Icges(6,2) * t128 + t170;
t171 = Icges(6,4) * t128;
t144 = Icges(6,1) * t130 - t171;
t92 = Icges(6,2) * t130 + t171;
t93 = Icges(6,1) * t128 + t170;
t145 = t128 * t92 - t130 * t93;
t190 = -t142 * qJD(5) + (t128 * t144 + t130 * t143 - t145) * t120;
t127 = qJ(1) + qJ(2);
t123 = qJ(3) + t127;
t116 = pkin(9) + t123;
t107 = sin(t116);
t108 = cos(t116);
t47 = -Icges(6,6) * t107 - t143 * t108;
t49 = -Icges(6,5) * t107 - t144 * t108;
t147 = t128 * t47 - t130 * t49;
t187 = t147 * t107;
t117 = sin(t123);
t118 = cos(t123);
t72 = t117 * rSges(4,1) + t118 * rSges(4,2);
t121 = sin(t127);
t122 = cos(t127);
t81 = t121 * rSges(3,1) + t122 * rSges(3,2);
t109 = pkin(3) * t117;
t54 = t107 * rSges(5,1) + t108 * rSges(5,2) + t109;
t91 = Icges(6,5) * t128 + Icges(6,6) * t130;
t186 = -Icges(6,3) * t120 + t91 * qJD(5);
t182 = 2 * m(3);
t181 = 2 * m(4);
t180 = 2 * m(5);
t179 = 2 * m(6);
t166 = t107 * t120;
t175 = rSges(6,1) * t130;
t90 = t108 * t175;
t176 = rSges(6,3) * t166 + t120 * t90;
t174 = rSges(4,2) * t117;
t173 = rSges(6,2) * t128;
t172 = pkin(1) * qJD(1);
t165 = t108 * t120;
t164 = t108 * t128;
t163 = t118 * t120;
t162 = t121 * t126;
t161 = t122 * t126;
t158 = pkin(2) * t162;
t157 = t107 * t173;
t129 = sin(qJ(1));
t156 = t129 * t172;
t114 = pkin(2) * t121;
t62 = t114 + t72;
t82 = t122 * rSges(3,1) - rSges(3,2) * t121;
t73 = t118 * rSges(4,1) - t174;
t52 = t114 + t54;
t65 = rSges(3,1) * t161 - rSges(3,2) * t162;
t61 = rSges(4,1) * t163 - t120 * t174;
t153 = t107 * t175 - t157;
t115 = pkin(2) * t122;
t63 = t115 + t73;
t110 = pkin(3) * t118;
t55 = t108 * rSges(5,1) - rSges(5,2) * t107 + t110;
t96 = rSges(6,1) * t128 + rSges(6,2) * t130;
t46 = -Icges(6,6) * t108 + t143 * t107;
t48 = -Icges(6,5) * t108 + t144 * t107;
t150 = t128 * t48 + t130 * t46;
t149 = t128 * t46 - t130 * t48;
t148 = t128 * t49 + t130 * t47;
t95 = pkin(2) * t161;
t43 = t61 + t95;
t88 = pkin(3) * t163;
t37 = rSges(5,1) * t165 - rSges(5,2) * t166 + t88;
t51 = rSges(6,2) * t164 - t107 * rSges(6,3) - t90;
t53 = t115 + t55;
t35 = t37 + t95;
t141 = (t144 - t92) * t160 + (t143 + t93) * t159;
t140 = -(-t147 * qJD(5) + t190 * t107) * t107 / 0.2e1 - (-t149 * qJD(5) + t190 * t108) * t108 / 0.2e1 + (-t145 * t107 - t108 * t91 + t150) * t166 / 0.2e1 - (-t107 * t91 + t145 * t108 + t148) * t165 / 0.2e1;
t64 = t81 * t126;
t60 = t72 * t120;
t139 = t149 * t108;
t138 = qJD(5) * t96;
t135 = t142 * t120;
t31 = t108 * pkin(4) + t107 * pkin(8) + t110 - t51;
t36 = t54 * t120;
t29 = t115 + t31;
t132 = rSges(6,3) * t165 - t108 * t138 + t120 * t157;
t30 = t107 * pkin(4) + t109 + (-rSges(6,3) - pkin(8)) * t108 + t153;
t42 = -t60 - t158;
t28 = t114 + t30;
t34 = -t36 - t158;
t15 = -t107 * rSges(6,1) * t160 + pkin(8) * t166 + pkin(4) * t165 + t88 + (-t107 * t159 - t120 * t164) * rSges(6,2) + t176;
t13 = t95 + t15;
t14 = pkin(8) * t165 + (-t109 + (-pkin(4) - t175) * t107) * t120 + t132;
t12 = t14 - t158;
t131 = cos(qJ(1));
t125 = t131 * pkin(1);
t124 = t129 * pkin(1);
t119 = t131 * t172;
t87 = (-t173 + t175) * qJD(5);
t67 = t125 + t82;
t66 = t124 + t81;
t59 = t119 + t65;
t58 = -t64 - t156;
t57 = t125 + t63;
t56 = t124 + t62;
t50 = -rSges(6,3) * t108 + t153;
t45 = -Icges(6,3) * t107 - t142 * t108;
t44 = -Icges(6,3) * t108 + t142 * t107;
t41 = t125 + t53;
t40 = t124 + t52;
t39 = t119 + t43;
t38 = t42 - t156;
t33 = t119 + t35;
t32 = t34 - t156;
t27 = t125 + t29;
t26 = t124 + t28;
t19 = -t186 * t107 + t108 * t135;
t18 = t107 * t135 + t186 * t108;
t11 = t119 + t13;
t10 = t12 - t156;
t9 = -t107 * t45 + t147 * t108;
t8 = -t107 * t44 + t139;
t7 = -t108 * t45 - t187;
t6 = -t149 * t107 - t108 * t44;
t1 = (t120 * t50 + t132) * t108 + (-t107 * t138 + (t51 + (-t173 - t175) * t108) * t120 + t176) * t107;
t2 = [(t58 * t67 + t59 * t66) * t182 + (t38 * t57 + t39 * t56) * t181 + (t32 * t41 + t33 * t40) * t180 + (t10 * t27 + t11 * t26) * t179 + t141; m(3) * (t58 * t82 + t59 * t81 - t64 * t67 + t65 * t66) + m(4) * (t38 * t63 + t39 * t62 + t42 * t57 + t43 * t56) + m(5) * (t32 * t53 + t33 * t52 + t34 * t41 + t35 * t40) + m(6) * (t10 * t29 + t11 * t28 + t12 * t27 + t13 * t26) + t141; (t12 * t29 + t13 * t28) * t179 + (t34 * t53 + t35 * t52) * t180 + (t42 * t63 + t43 * t62) * t181 + (-t64 * t82 + t65 * t81) * t182 + t141; m(4) * (t38 * t73 + t39 * t72 + t56 * t61 - t57 * t60) + m(5) * (t32 * t55 + t33 * t54 - t36 * t41 + t37 * t40) + m(6) * (t10 * t31 + t11 * t30 + t14 * t27 + t15 * t26) + t141; m(6) * (t12 * t31 + t13 * t30 + t14 * t29 + t15 * t28) + m(5) * (t34 * t55 + t35 * t54 - t36 * t53 + t37 * t52) + m(4) * (t42 * t73 + t43 * t72 - t60 * t63 + t61 * t62) + t141; (-t60 * t73 + t61 * t72) * t181 + (-t36 * t55 + t37 * t54) * t180 + (t14 * t31 + t15 * t30) * t179 + t141; 0; 0; 0; 0; m(6) * ((-t107 * t27 + t108 * t26) * t87 + ((-t120 * t27 + t11) * t108 + (-t120 * t26 - t10) * t107) * t96) + t140; m(6) * ((-t107 * t29 + t108 * t28) * t87 + ((-t120 * t29 + t13) * t108 + (-t120 * t28 - t12) * t107) * t96) + t140; m(6) * ((-t107 * t31 + t108 * t30) * t87 + ((-t120 * t31 + t15) * t108 + (-t120 * t30 - t14) * t107) * t96) + t140; m(6) * t1; ((t107 * t50 - t108 * t51) * t1 + (t107 ^ 2 + t108 ^ 2) * t96 * t87) * t179 + (-t107 * t7 - t108 * t6) * t166 - t108 * ((t108 * t19 + (-t7 + t139) * t120) * t108 + (t6 * t120 + (t47 * t159 + t49 * t160) * t107 + (t150 * qJD(5) + t147 * t120 + t18) * t108) * t107) - (-t107 * t9 - t108 * t8) * t165 - t107 * ((t107 * t18 + (t8 + t187) * t120) * t107 + (-t9 * t120 + (-t46 * t159 - t48 * t160) * t108 + (-t148 * qJD(5) + t149 * t120 + t19) * t107) * t108);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
