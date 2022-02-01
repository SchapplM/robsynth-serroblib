% Calculate time derivative of joint inertia matrix for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:04
% EndTime: 2022-01-23 08:59:08
% DurationCPUTime: 1.66s
% Computational Cost: add. (4249->313), mult. (11317->477), div. (0->0), fcn. (12809->10), ass. (0->152)
t134 = sin(qJ(1));
t136 = cos(qJ(1));
t127 = sin(pkin(9));
t130 = cos(pkin(9));
t115 = pkin(4) * t130 + pkin(6) * t127 + pkin(3);
t129 = sin(pkin(7));
t131 = cos(pkin(8));
t132 = cos(pkin(7));
t128 = sin(pkin(8));
t168 = t128 * qJ(4) + pkin(2);
t75 = (t115 * t131 + t168) * t132 + pkin(1) + (pkin(4) * t127 - pkin(6) * t130 + qJ(3)) * t129;
t143 = qJ(4) * t131 - qJ(2);
t95 = -t115 * t128 + t143;
t177 = t134 * t95 - t75 * t136;
t176 = -t75 * t134 - t136 * t95;
t138 = rSges(3,1) * t132 - rSges(3,2) * t129 + pkin(1);
t175 = t138 * t136;
t125 = t134 * qJ(2);
t169 = t129 * qJ(3) + pkin(1);
t144 = rSges(4,3) * t129 + pkin(2) * t132 + t169;
t174 = -t136 * t144 - t125;
t173 = 2 * m(6);
t172 = 0.2e1 * t129;
t171 = m(5) / 0.2e1;
t170 = m(6) / 0.2e1;
t167 = rSges(5,1) * t130;
t113 = -pkin(3) * t128 + t143;
t164 = t113 * t134;
t163 = t113 * t136;
t133 = sin(qJ(5));
t162 = t128 * t133;
t135 = cos(qJ(5));
t161 = t128 * t135;
t160 = t129 * t127;
t159 = t129 * t131;
t158 = t129 * t134;
t157 = t129 * t136;
t156 = t132 * t130;
t155 = t134 * t128;
t154 = t134 * t131;
t153 = t136 * t128;
t152 = t136 * t131;
t147 = t129 * qJD(3);
t112 = qJD(4) * t128 * t132 + t147;
t117 = qJD(4) * t131 - qJD(2);
t151 = t112 * t136 - t117 * t134;
t148 = qJD(1) * t136;
t150 = qJ(2) * t148 + qJD(2) * t134;
t149 = qJD(1) * t134;
t146 = t171 + t170;
t107 = t130 * t161 - t131 * t133;
t106 = t130 * t162 + t131 * t135;
t105 = t131 * t156 + t160;
t89 = -t105 * t133 + t132 * t161;
t74 = -t106 * t134 + t136 * t89;
t88 = t105 * t135 + t132 * t162;
t47 = (t107 * t136 - t134 * t88) * qJD(1) + t74 * qJD(5);
t73 = t106 * t136 + t134 * t89;
t83 = t88 * qJD(5);
t97 = t107 * qJD(5);
t50 = -qJD(1) * t73 - t134 * t97 - t83 * t136;
t145 = qJD(1) * t129 * t130;
t109 = t132 * t154 - t153;
t99 = t109 * qJD(1);
t79 = -t127 * t99 + t134 * t145;
t27 = rSges(6,1) * t47 + rSges(6,2) * t50 + t79 * rSges(6,3);
t72 = t107 * t134 + t136 * t88;
t111 = t132 * t152 + t155;
t91 = t111 * t127 - t130 * t157;
t40 = rSges(6,1) * t72 + rSges(6,2) * t74 + t91 * rSges(6,3);
t103 = t127 * t159 + t156;
t104 = -t127 * t132 + t130 * t159;
t85 = t104 * t135 + t129 * t162;
t86 = -t104 * t133 + t129 * t161;
t59 = Icges(6,4) * t85 + Icges(6,2) * t86 + Icges(6,6) * t103;
t60 = Icges(6,1) * t85 + Icges(6,4) * t86 + Icges(6,5) * t103;
t81 = t86 * qJD(5);
t82 = t85 * qJD(5);
t63 = Icges(6,5) * t81 - Icges(6,6) * t82;
t64 = Icges(6,4) * t81 - Icges(6,2) * t82;
t65 = Icges(6,1) * t81 - Icges(6,4) * t82;
t142 = t103 * t63 - t59 * t82 + t60 * t81 + t64 * t86 + t65 * t85;
t94 = (pkin(3) * t131 + t168) * t132 + t169;
t141 = -rSges(5,1) * t160 - t94;
t140 = t144 * t134;
t139 = -t112 * t134 - t117 * t136;
t110 = t132 * t153 - t154;
t108 = t132 * t155 + t152;
t100 = t110 * qJD(1);
t87 = t105 * t134 - t130 * t153;
t48 = t100 * t133 + (t108 * t135 - t133 * t87) * qJD(5) + (t105 * t136 + t130 * t155) * t135 * qJD(1);
t49 = qJD(1) * t74 - t83 * t134 + t136 * t97;
t101 = t111 * qJD(1);
t80 = t101 * t127 - t136 * t145;
t28 = t48 * rSges(6,1) + t49 * rSges(6,2) + t80 * rSges(6,3);
t71 = t108 * t133 + t135 * t87;
t90 = t109 * t127 - t130 * t158;
t39 = rSges(6,1) * t71 + rSges(6,2) * t73 + rSges(6,3) * t90;
t137 = t136 * rSges(3,3) - t134 * t138;
t126 = t136 * qJ(2);
t123 = qJD(2) * t136;
t98 = t108 * qJD(1);
t93 = t134 * rSges(3,3) + t125 + t175;
t92 = t126 + t137;
t77 = t123 + ((-rSges(3,3) - qJ(2)) * t134 - t175) * qJD(1);
t76 = qJD(1) * t137 + t150;
t68 = -t109 * rSges(4,1) + t108 * rSges(4,2) + t126 - t140;
t67 = t111 * rSges(4,1) - t110 * rSges(4,2) - t174;
t66 = rSges(6,1) * t81 - rSges(6,2) * t82;
t61 = rSges(6,1) * t85 + rSges(6,2) * t86 + rSges(6,3) * t103;
t58 = Icges(6,5) * t85 + Icges(6,6) * t86 + Icges(6,3) * t103;
t57 = -t99 * rSges(4,1) + t98 * rSges(4,2) - qJD(1) * t140 + t136 * t147 + t150;
t56 = -t101 * rSges(4,1) + t100 * rSges(4,2) + qJD(1) * t174 - t134 * t147 + t123;
t42 = t94 * t136 - t164 + (t111 * t130 + t127 * t157) * rSges(5,1) - t91 * rSges(5,2) + t110 * rSges(5,3);
t41 = -t94 * t134 - t163 - (t109 * t130 + t127 * t158) * rSges(5,1) + t90 * rSges(5,2) - t108 * rSges(5,3);
t38 = Icges(6,1) * t72 + Icges(6,4) * t74 + Icges(6,5) * t91;
t37 = Icges(6,1) * t71 + Icges(6,4) * t73 + Icges(6,5) * t90;
t36 = Icges(6,4) * t72 + Icges(6,2) * t74 + Icges(6,6) * t91;
t35 = Icges(6,4) * t71 + Icges(6,2) * t73 + Icges(6,6) * t90;
t34 = Icges(6,5) * t72 + Icges(6,6) * t74 + Icges(6,3) * t91;
t33 = Icges(6,5) * t71 + Icges(6,6) * t73 + Icges(6,3) * t90;
t32 = -t99 * t167 - t79 * rSges(5,2) - t98 * rSges(5,3) + (t134 * t141 - t163) * qJD(1) + t151;
t31 = -t101 * t167 + t80 * rSges(5,2) - t100 * rSges(5,3) + (t136 * t141 + t164) * qJD(1) + t139;
t30 = -t177 + t40;
t29 = t176 - t39;
t26 = Icges(6,1) * t48 + Icges(6,4) * t49 + Icges(6,5) * t80;
t25 = Icges(6,1) * t47 + Icges(6,4) * t50 + Icges(6,5) * t79;
t24 = Icges(6,4) * t48 + Icges(6,2) * t49 + Icges(6,6) * t80;
t23 = Icges(6,4) * t47 + Icges(6,2) * t50 + Icges(6,6) * t79;
t22 = Icges(6,5) * t48 + Icges(6,6) * t49 + Icges(6,3) * t80;
t21 = Icges(6,5) * t47 + Icges(6,6) * t50 + Icges(6,3) * t79;
t20 = t103 * t40 - t61 * t91;
t19 = -t103 * t39 + t61 * t90;
t18 = qJD(1) * t176 + t151 + t27;
t17 = qJD(1) * t177 + t139 - t28;
t16 = t58 * t91 + t59 * t74 + t60 * t72;
t15 = t58 * t90 + t59 * t73 + t60 * t71;
t14 = -t103 * t28 + t61 * t80 + t66 * t90;
t13 = t103 * t27 - t61 * t79 - t66 * t91;
t12 = t103 * t34 + t36 * t86 + t38 * t85;
t11 = t103 * t33 + t35 * t86 + t37 * t85;
t10 = t142 * t103;
t9 = t34 * t91 + t36 * t74 + t38 * t72;
t8 = t33 * t91 + t35 * t74 + t37 * t72;
t7 = t34 * t90 + t36 * t73 + t38 * t71;
t6 = t33 * t90 + t35 * t73 + t37 * t71;
t5 = t48 * t60 + t49 * t59 + t58 * t80 + t63 * t90 + t64 * t73 + t65 * t71;
t4 = t47 * t60 + t50 * t59 + t58 * t79 + t63 * t91 + t64 * t74 + t65 * t72;
t3 = -t27 * t90 + t28 * t91 + t39 * t79 - t40 * t80;
t2 = t103 * t21 + t23 * t86 + t25 * t85 - t36 * t82 + t38 * t81;
t1 = t103 * t22 + t24 * t86 + t26 * t85 - t35 * t82 + t37 * t81;
t43 = [(t17 * t29 + t18 * t30) * t173 + 0.2e1 * m(5) * (t31 * t41 + t32 * t42) + 0.2e1 * m(4) * (t56 * t68 + t57 * t67) + 0.2e1 * m(3) * (t76 * t93 + t77 * t92) + t142; m(6) * (t134 * t17 - t136 * t18 + (t134 * t30 + t136 * t29) * qJD(1)) + m(5) * (t134 * t31 - t136 * t32 + (t134 * t42 + t136 * t41) * qJD(1)) + m(4) * (t134 * t56 - t136 * t57 + (t134 * t67 + t136 * t68) * qJD(1)) + m(3) * (t134 * t77 - t136 * t76 + (t134 * t93 + t136 * t92) * qJD(1)); 0; ((t134 * t18 + t136 * t17 + t148 * t30 - t149 * t29) * t170 + (t134 * t32 + t136 * t31 + t148 * t42 - t149 * t41) * t171 + m(4) * (t134 * t57 + t136 * t56 + t148 * t67 - t149 * t68) / 0.2e1) * t172; 0; 0; m(6) * (t100 * t30 + t108 * t18 + t110 * t17 - t29 * t98) + m(5) * (t100 * t42 + t108 * t32 + t110 * t31 - t41 * t98); 0.2e1 * t146 * (-t100 * t136 - t98 * t134 + (t108 * t134 + t110 * t136) * qJD(1)); t146 * (t100 * t134 - t136 * t98 + (t108 * t136 - t110 * t134) * qJD(1)) * t172; 0.4e1 * t146 * (t100 * t108 - t110 * t98); t10 + m(6) * (t13 * t30 + t14 * t29 + t17 * t19 + t18 * t20) + (t2 / 0.2e1 + t4 / 0.2e1) * t91 + (t5 / 0.2e1 + t1 / 0.2e1) * t90 + (t15 / 0.2e1 + t11 / 0.2e1) * t80 + (t12 / 0.2e1 + t16 / 0.2e1) * t79; m(6) * (-t13 * t136 + t14 * t134 + (t134 * t20 + t136 * t19) * qJD(1)); m(6) * (-t3 * t132 + (t13 * t134 + t136 * t14 + (-t134 * t19 + t136 * t20) * qJD(1)) * t129); m(6) * (t128 * t129 * t3 + t100 * t20 + t108 * t13 + t110 * t14 - t19 * t98); (t19 * t14 + t20 * t13 + (t39 * t91 - t40 * t90) * t3) * t173 + t79 * (t16 * t103 + t8 * t90 + t9 * t91) + t91 * ((t21 * t91 + t23 * t74 + t25 * t72 + t34 * t79 + t36 * t50 + t38 * t47) * t91 + t9 * t79 + (t22 * t91 + t24 * t74 + t26 * t72 + t33 * t79 + t35 * t50 + t37 * t47) * t90 + t8 * t80 + t4 * t103) + t80 * (t15 * t103 + t6 * t90 + t7 * t91) + t90 * ((t21 * t90 + t23 * t73 + t25 * t71 + t34 * t80 + t36 * t49 + t38 * t48) * t91 + t7 * t79 + (t22 * t90 + t24 * t73 + t26 * t71 + t33 * t80 + t35 * t49 + t37 * t48) * t90 + t6 * t80 + t5 * t103) + t103 * (t1 * t90 + t11 * t80 + t12 * t79 + t2 * t91 + t10);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t43(1), t43(2), t43(4), t43(7), t43(11); t43(2), t43(3), t43(5), t43(8), t43(12); t43(4), t43(5), t43(6), t43(9), t43(13); t43(7), t43(8), t43(9), t43(10), t43(14); t43(11), t43(12), t43(13), t43(14), t43(15);];
Mq = res;
