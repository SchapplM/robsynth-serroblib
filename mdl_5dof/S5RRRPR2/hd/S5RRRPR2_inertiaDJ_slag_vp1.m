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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:30:33
% EndTime: 2022-01-20 11:30:37
% DurationCPUTime: 1.53s
% Computational Cost: add. (5296->216), mult. (3006->294), div. (0->0), fcn. (2068->10), ass. (0->131)
t116 = sin(qJ(5));
t151 = qJD(5) * t116;
t118 = cos(qJ(5));
t150 = qJD(5) * t118;
t115 = qJ(1) + qJ(2);
t112 = qJ(3) + t115;
t107 = sin(t112);
t108 = cos(t112);
t73 = t108 * rSges(4,1) - rSges(4,2) * t107;
t106 = pkin(9) + t112;
t102 = cos(t106);
t103 = pkin(3) * t108;
t175 = -t102 * rSges(5,1) - t103;
t101 = sin(t106);
t156 = Icges(6,4) * t118;
t132 = -Icges(6,2) * t116 + t156;
t127 = t132 * t102;
t47 = Icges(6,6) * t101 + t127;
t157 = Icges(6,4) * t116;
t133 = Icges(6,1) * t118 - t157;
t128 = t133 * t102;
t49 = Icges(6,5) * t101 + t128;
t135 = t116 * t47 - t118 * t49;
t174 = t135 * t101;
t46 = -Icges(6,6) * t102 + t132 * t101;
t48 = -Icges(6,5) * t102 + t133 * t101;
t137 = t116 * t46 - t118 * t48;
t173 = t137 * t102;
t114 = qJD(1) + qJD(2);
t109 = qJD(3) + t114;
t131 = Icges(6,5) * t118 - Icges(6,6) * t116;
t90 = Icges(6,2) * t118 + t157;
t91 = Icges(6,1) * t116 + t156;
t134 = t116 * t90 - t118 * t91;
t172 = t131 * qJD(5) + t134 * t109;
t171 = 2 * m(3);
t170 = 2 * m(4);
t169 = 2 * m(5);
t168 = 2 * m(6);
t117 = sin(qJ(1));
t165 = pkin(1) * t117;
t110 = sin(t115);
t164 = pkin(2) * t110;
t163 = pkin(3) * t107;
t159 = rSges(6,2) * t116;
t87 = t101 * t159;
t162 = t102 * rSges(6,3) + t87;
t161 = rSges(6,1) * t118;
t158 = pkin(1) * qJD(1);
t94 = t101 * rSges(6,3);
t155 = t101 * t109;
t154 = t102 * t109;
t153 = t110 * t114;
t111 = cos(t115);
t152 = t111 * t114;
t146 = t102 * t159;
t149 = -t109 * t146 + (-rSges(6,1) * t151 - rSges(6,2) * t150) * t101;
t148 = pkin(2) * t153;
t147 = pkin(2) * t152;
t145 = t117 * t158;
t119 = cos(qJ(1));
t144 = t119 * t158;
t141 = -pkin(4) - t161;
t79 = t111 * rSges(3,1) - rSges(3,2) * t110;
t65 = -rSges(3,1) * t152 + rSges(3,2) * t153;
t61 = t73 * t109;
t55 = -rSges(5,2) * t101 - t175;
t105 = pkin(2) * t111;
t63 = t105 + t73;
t78 = -rSges(3,1) * t110 - rSges(3,2) * t111;
t72 = -rSges(4,1) * t107 - rSges(4,2) * t108;
t93 = rSges(6,1) * t116 + rSges(6,2) * t118;
t138 = t116 * t48 + t118 * t46;
t136 = t116 * t49 + t118 * t47;
t89 = Icges(6,5) * t116 + Icges(6,6) * t118;
t51 = t102 * t161 - t146 + t94;
t53 = t105 + t55;
t130 = (t133 - t90) * t151 + (t132 + t91) * t150;
t129 = (-t135 * qJD(5) + t172 * t101 + (-t116 * t133 - t118 * t132) * t155) * t101 / 0.2e1 - (-t137 * qJD(5) - t172 * t102 + (t116 * t128 + t118 * t127) * t109) * t102 / 0.2e1 + (-t134 * t101 - t102 * t89 + t138) * t155 / 0.2e1 + (t101 * t89 - t134 * t102 + t136) * t154 / 0.2e1;
t64 = t78 * t114;
t60 = t72 * t109;
t126 = t131 * t102;
t62 = t72 - t164;
t54 = -rSges(5,1) * t101 - rSges(5,2) * t102 - t163;
t125 = t141 * t101 - t163;
t43 = -t61 - t147;
t37 = rSges(5,2) * t155 + t175 * t109;
t31 = t102 * pkin(4) + t101 * pkin(8) + t103 + t51;
t36 = t54 * t109;
t29 = t105 + t31;
t52 = t54 - t164;
t124 = -t93 * t102 * qJD(5) + rSges(6,3) * t154 + t109 * t87;
t121 = Icges(6,3) * t109 - t89 * qJD(5);
t30 = t102 * pkin(8) + t125 + t162;
t42 = t60 - t148;
t35 = t37 - t147;
t28 = t30 - t164;
t34 = t36 - t148;
t15 = (-t103 + t141 * t102 + (-rSges(6,3) - pkin(8)) * t101) * t109 - t149;
t13 = t15 - t147;
t14 = pkin(8) * t154 + t125 * t109 + t124;
t12 = t14 - t148;
t113 = t119 * pkin(1);
t86 = (-t159 + t161) * qJD(5);
t67 = t113 + t79;
t66 = t78 - t165;
t59 = t65 - t144;
t58 = t64 - t145;
t57 = t113 + t63;
t56 = t62 - t165;
t50 = t101 * t161 - t162;
t45 = Icges(6,3) * t101 + t126;
t44 = -Icges(6,3) * t102 + t131 * t101;
t41 = t113 + t53;
t40 = t52 - t165;
t39 = t43 - t144;
t38 = t42 - t145;
t33 = t35 - t144;
t32 = t34 - t145;
t27 = t113 + t29;
t26 = t28 - t165;
t19 = t121 * t101 + t109 * t126;
t18 = t121 * t102 - t131 * t155;
t11 = t13 - t144;
t10 = t12 - t145;
t9 = t101 * t45 - t135 * t102;
t8 = t101 * t44 - t173;
t7 = -t102 * t45 - t174;
t6 = -t137 * t101 - t102 * t44;
t1 = ((-t51 + t94) * t109 + t149) * t101 + (t109 * t50 + t124) * t102;
t2 = [(t32 * t41 + t33 * t40) * t169 + (t10 * t27 + t11 * t26) * t168 + (t38 * t57 + t39 * t56) * t170 + (t58 * t67 + t59 * t66) * t171 + t130; m(5) * (t32 * t53 + t33 * t52 + t34 * t41 + t35 * t40) + m(6) * (t10 * t29 + t11 * t28 + t12 * t27 + t13 * t26) + m(4) * (t38 * t63 + t39 * t62 + t42 * t57 + t43 * t56) + m(3) * (t58 * t79 + t59 * t78 + t64 * t67 + t65 * t66) + t130; (t12 * t29 + t13 * t28) * t168 + (t34 * t53 + t35 * t52) * t169 + (t42 * t63 + t43 * t62) * t170 + (t64 * t79 + t65 * t78) * t171 + t130; m(5) * (t32 * t55 + t33 * t54 + t36 * t41 + t37 * t40) + m(6) * (t10 * t31 + t11 * t30 + t14 * t27 + t15 * t26) + m(4) * (t38 * t73 + t39 * t72 - t56 * t61 + t57 * t60) + t130; m(6) * (t12 * t31 + t13 * t30 + t14 * t29 + t15 * t28) + m(5) * (t34 * t55 + t35 * t54 + t36 * t53 + t37 * t52) + m(4) * (t42 * t73 + t43 * t72 + t60 * t63 - t61 * t62) + t130; (t60 * t73 - t61 * t72) * t170 + (t36 * t55 + t37 * t54) * t169 + (t14 * t31 + t15 * t30) * t168 + t130; 0; 0; 0; 0; m(6) * ((-t101 * t27 - t102 * t26) * t86 + ((-t109 * t27 - t11) * t102 + (t109 * t26 - t10) * t101) * t93) + t129; m(6) * ((-t101 * t29 - t102 * t28) * t86 + ((-t109 * t29 - t13) * t102 + (t109 * t28 - t12) * t101) * t93) + t129; m(6) * ((-t101 * t31 - t102 * t30) * t86 + ((-t109 * t31 - t15) * t102 + (t109 * t30 - t14) * t101) * t93) + t129; m(6) * t1; ((t101 * t50 + t102 * t51) * t1 + (t101 ^ 2 + t102 ^ 2) * t93 * t86) * t168 + (t101 * t9 - t102 * t8) * t154 + t101 * ((t101 * t18 + (t8 + t174) * t109) * t101 + (t9 * t109 + (t46 * t150 + t48 * t151) * t102 + (-t136 * qJD(5) - t137 * t109 - t19) * t101) * t102) + (t101 * t7 - t102 * t6) * t155 - t102 * ((t102 * t19 + (t7 + t173) * t109) * t102 + (t6 * t109 + (-t47 * t150 - t151 * t49) * t101 + (t138 * qJD(5) - t135 * t109 - t18) * t102) * t101);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
