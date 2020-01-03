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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:22:21
% EndTime: 2020-01-03 11:22:26
% DurationCPUTime: 1.47s
% Computational Cost: add. (4249->314), mult. (12011->458), div. (0->0), fcn. (13377->10), ass. (0->140)
t179 = 2 * m(6);
t130 = sin(pkin(7));
t178 = 0.2e1 * t130;
t177 = m(5) / 0.2e1;
t176 = m(6) / 0.2e1;
t175 = -rSges(6,3) - pkin(6);
t132 = cos(pkin(7));
t174 = pkin(2) * t132;
t173 = rSges(5,3) + qJ(4);
t172 = cos(pkin(8));
t129 = sin(pkin(8));
t171 = t129 * t130;
t134 = sin(qJ(1));
t170 = t130 * t134;
t136 = cos(qJ(1));
t169 = t130 * t136;
t168 = t134 * t129;
t167 = t136 * t129;
t162 = qJD(1) * t136;
t163 = qJD(1) * t134;
t166 = pkin(1) * t162 + qJ(2) * t163;
t165 = qJ(2) * t162 + qJD(2) * t134;
t164 = t136 * pkin(1) + t134 * qJ(2);
t161 = qJD(2) * t136;
t160 = qJD(3) * t130;
t159 = t177 + t176;
t133 = sin(qJ(5));
t135 = cos(qJ(5));
t151 = t136 * t172;
t104 = t132 * t168 + t151;
t152 = t134 * t172;
t105 = t132 * t152 - t167;
t128 = sin(pkin(9));
t131 = cos(pkin(9));
t88 = t105 * t131 + t128 * t170;
t72 = t104 * t133 + t88 * t135;
t155 = t130 * t162;
t140 = t132 * t151 + t168;
t98 = t140 * qJD(1);
t81 = t128 * t155 + t98 * t131;
t154 = t132 * t162;
t97 = -qJD(1) * t152 + t129 * t154;
t49 = -t72 * qJD(5) - t81 * t133 + t97 * t135;
t71 = t104 * t135 - t88 * t133;
t50 = t71 * qJD(5) + t97 * t133 + t81 * t135;
t80 = t98 * t128 - t131 * t155;
t28 = t50 * rSges(6,1) + t49 * rSges(6,2) + t80 * rSges(6,3);
t87 = t105 * t128 - t131 * t170;
t39 = t72 * rSges(6,1) + t71 * rSges(6,2) + t87 * rSges(6,3);
t158 = t136 * t160 + t165;
t157 = -pkin(1) - t174;
t156 = t130 * t163;
t153 = t130 * t172;
t150 = qJ(3) * t169 + t136 * t174 + t164;
t102 = t128 * t153 + t132 * t131;
t103 = -t132 * t128 + t131 * t153;
t85 = -t103 * t133 + t135 * t171;
t86 = t103 * t135 + t133 * t171;
t59 = Icges(6,4) * t86 + Icges(6,2) * t85 + Icges(6,6) * t102;
t60 = Icges(6,1) * t86 + Icges(6,4) * t85 + Icges(6,5) * t102;
t82 = t85 * qJD(5);
t83 = t86 * qJD(5);
t63 = Icges(6,5) * t82 - Icges(6,6) * t83;
t64 = Icges(6,4) * t82 - Icges(6,2) * t83;
t65 = Icges(6,1) * t82 - Icges(6,4) * t83;
t149 = t102 * t63 - t83 * t59 + t82 * t60 + t85 * t64 + t86 * t65;
t106 = t132 * t167 - t152;
t90 = -t128 * t169 - t131 * t140;
t144 = t106 * t133 - t90 * t135;
t96 = t105 * qJD(1);
t79 = t128 * t156 + t96 * t131;
t95 = t104 * qJD(1);
t47 = t144 * qJD(5) - t79 * t133 + t95 * t135;
t73 = -t106 * t135 - t90 * t133;
t48 = t73 * qJD(5) + t95 * t133 + t79 * t135;
t148 = -t48 * rSges(6,1) - t47 * rSges(6,2);
t147 = rSges(6,1) * t144 - t73 * rSges(6,2);
t146 = rSges(3,1) * t132 - rSges(3,2) * t130;
t145 = pkin(2) * t154 + qJ(3) * t155 + t134 * t160 + t166;
t126 = t134 * pkin(1);
t143 = -t136 * qJ(2) + qJ(3) * t170 + t134 * t174 + t126;
t142 = t105 * pkin(3) + t143;
t141 = pkin(3) * t140 + t106 * qJ(4) + t150;
t139 = t134 * rSges(3,3) + t146 * t136;
t138 = t98 * pkin(3) + t104 * qJD(4) + t145 - t161;
t137 = t106 * qJD(4) + (-qJ(3) * t130 + t157) * t163 - t95 * qJ(4) - t96 * pkin(3) + t158;
t92 = t139 + t164;
t91 = t126 + (-rSges(3,3) - qJ(2)) * t136 + t146 * t134;
t89 = -t128 * t140 + t131 * t169;
t78 = t96 * t128 - t131 * t156;
t76 = t139 * qJD(1) - t161 + t166;
t75 = (t136 * rSges(3,3) + (-pkin(1) - t146) * t134) * qJD(1) + t165;
t68 = rSges(4,1) * t140 - t106 * rSges(4,2) + rSges(4,3) * t169 + t150;
t67 = t105 * rSges(4,1) - t104 * rSges(4,2) + rSges(4,3) * t170 + t143;
t66 = t82 * rSges(6,1) - t83 * rSges(6,2);
t61 = t86 * rSges(6,1) + t85 * rSges(6,2) + t102 * rSges(6,3);
t58 = Icges(6,5) * t86 + Icges(6,6) * t85 + Icges(6,3) * t102;
t55 = t98 * rSges(4,1) - t97 * rSges(4,2) + (rSges(4,3) * qJD(1) * t130 - qJD(2)) * t136 + t145;
t54 = -t96 * rSges(4,1) + t95 * rSges(4,2) + ((-rSges(4,3) - qJ(3)) * t130 + t157) * t163 + t158;
t42 = -t90 * rSges(5,1) + t89 * rSges(5,2) + t106 * rSges(5,3) + t141;
t41 = t88 * rSges(5,1) - t87 * rSges(5,2) + t173 * t104 + t142;
t40 = t89 * rSges(6,3) - t147;
t38 = -Icges(6,1) * t144 + Icges(6,4) * t73 + Icges(6,5) * t89;
t37 = Icges(6,1) * t72 + Icges(6,4) * t71 + Icges(6,5) * t87;
t36 = -Icges(6,4) * t144 + Icges(6,2) * t73 + Icges(6,6) * t89;
t35 = Icges(6,4) * t72 + Icges(6,2) * t71 + Icges(6,6) * t87;
t34 = -Icges(6,5) * t144 + Icges(6,6) * t73 + Icges(6,3) * t89;
t33 = Icges(6,5) * t72 + Icges(6,6) * t71 + Icges(6,3) * t87;
t32 = t81 * rSges(5,1) - t80 * rSges(5,2) + t173 * t97 + t138;
t31 = -t79 * rSges(5,1) + t78 * rSges(5,2) - t95 * rSges(5,3) + t137;
t30 = -t90 * pkin(4) + t175 * t89 + t141 + t147;
t29 = t88 * pkin(4) + t87 * pkin(6) + t104 * qJ(4) + t142 + t39;
t27 = t78 * rSges(6,3) - t148;
t26 = Icges(6,1) * t50 + Icges(6,4) * t49 + Icges(6,5) * t80;
t25 = Icges(6,1) * t48 + Icges(6,4) * t47 + Icges(6,5) * t78;
t24 = Icges(6,4) * t50 + Icges(6,2) * t49 + Icges(6,6) * t80;
t23 = Icges(6,4) * t48 + Icges(6,2) * t47 + Icges(6,6) * t78;
t22 = Icges(6,5) * t50 + Icges(6,6) * t49 + Icges(6,3) * t80;
t21 = Icges(6,5) * t48 + Icges(6,6) * t47 + Icges(6,3) * t78;
t20 = -t102 * t40 + t89 * t61;
t19 = t102 * t39 - t87 * t61;
t18 = -t144 * t60 + t89 * t58 + t73 * t59;
t17 = t87 * t58 + t71 * t59 + t72 * t60;
t16 = t81 * pkin(4) + t80 * pkin(6) + t97 * qJ(4) + t138 + t28;
t15 = -t79 * pkin(4) + t175 * t78 + t137 + t148;
t14 = t102 * t28 - t80 * t61 - t87 * t66;
t13 = -t102 * t27 + t78 * t61 + t89 * t66;
t12 = t102 * t34 + t85 * t36 + t86 * t38;
t11 = t102 * t33 + t85 * t35 + t86 * t37;
t10 = t149 * t102;
t9 = -t144 * t38 + t89 * t34 + t73 * t36;
t8 = -t144 * t37 + t89 * t33 + t73 * t35;
t7 = t87 * t34 + t71 * t36 + t72 * t38;
t6 = t87 * t33 + t71 * t35 + t72 * t37;
t5 = t49 * t59 + t50 * t60 + t80 * t58 + t87 * t63 + t71 * t64 + t72 * t65;
t4 = -t144 * t65 + t47 * t59 + t48 * t60 + t78 * t58 + t89 * t63 + t73 * t64;
t3 = t87 * t27 - t89 * t28 - t78 * t39 + t80 * t40;
t2 = t102 * t21 + t85 * t23 + t86 * t25 - t83 * t36 + t82 * t38;
t1 = t102 * t22 + t85 * t24 + t86 * t26 - t83 * t35 + t82 * t37;
t43 = [0.2e1 * m(3) * (t92 * t75 + t91 * t76) + 0.2e1 * m(4) * (t68 * t54 + t67 * t55) + 0.2e1 * m(5) * (t42 * t31 + t41 * t32) + (t30 * t15 + t29 * t16) * t179 + t149; m(3) * (-t134 * t76 - t136 * t75 + (t134 * t92 - t136 * t91) * qJD(1)) + m(4) * (-t134 * t55 - t136 * t54 + (t134 * t68 - t136 * t67) * qJD(1)) + m(5) * (-t134 * t32 - t136 * t31 + (t134 * t42 - t136 * t41) * qJD(1)) + m(6) * (-t134 * t16 - t136 * t15 + (t134 * t30 - t136 * t29) * qJD(1)); 0; (m(4) * (t134 * t54 - t136 * t55 + t68 * t162 + t67 * t163) / 0.2e1 + (t134 * t31 - t136 * t32 + t42 * t162 + t41 * t163) * t177 + (t134 * t15 - t136 * t16 + t30 * t162 + t29 * t163) * t176) * t178; 0; 0; m(5) * (t104 * t31 - t106 * t32 + t95 * t41 + t97 * t42) + m(6) * (t104 * t15 - t106 * t16 + t95 * t29 + t97 * t30); 0.2e1 * t159 * (-t95 * t134 - t97 * t136 + (t104 * t134 + t106 * t136) * qJD(1)); t159 * (t134 * t97 - t136 * t95 + (t104 * t136 - t106 * t134) * qJD(1)) * t178; 0.4e1 * t159 * (t104 * t97 - t106 * t95); t10 + m(6) * (t13 * t30 + t14 * t29 + t20 * t15 + t19 * t16) + (t2 / 0.2e1 + t4 / 0.2e1) * t89 + (t1 / 0.2e1 + t5 / 0.2e1) * t87 + (t11 / 0.2e1 + t17 / 0.2e1) * t80 + (t12 / 0.2e1 + t18 / 0.2e1) * t78; m(6) * (-t13 * t136 - t14 * t134 + (t134 * t20 - t136 * t19) * qJD(1)); m(6) * (-t3 * t132 + (t13 * t134 - t136 * t14 + (t134 * t19 + t136 * t20) * qJD(1)) * t130); m(6) * (t13 * t104 - t14 * t106 + t3 * t171 + t19 * t95 + t20 * t97); ((-t89 * t39 + t87 * t40) * t3 + t20 * t13 + t19 * t14) * t179 + t102 * (t1 * t87 + t11 * t80 + t12 * t78 + t2 * t89 + t10) + t80 * (t17 * t102 + t6 * t87 + t7 * t89) + t87 * (t5 * t102 + (t87 * t22 + t71 * t24 + t72 * t26 + t80 * t33 + t49 * t35 + t50 * t37) * t87 + t6 * t80 + (t87 * t21 + t71 * t23 + t72 * t25 + t80 * t34 + t49 * t36 + t50 * t38) * t89 + t7 * t78) + t78 * (t18 * t102 + t8 * t87 + t9 * t89) + t89 * (t4 * t102 + (-t144 * t26 + t89 * t22 + t73 * t24 + t78 * t33 + t47 * t35 + t48 * t37) * t87 + t8 * t80 + (-t144 * t25 + t89 * t21 + t73 * t23 + t78 * t34 + t47 * t36 + t48 * t38) * t89 + t9 * t78);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t43(1), t43(2), t43(4), t43(7), t43(11); t43(2), t43(3), t43(5), t43(8), t43(12); t43(4), t43(5), t43(6), t43(9), t43(13); t43(7), t43(8), t43(9), t43(10), t43(14); t43(11), t43(12), t43(13), t43(14), t43(15);];
Mq = res;
