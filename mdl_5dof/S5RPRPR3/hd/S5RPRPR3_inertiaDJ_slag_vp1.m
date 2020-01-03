% Calculate time derivative of joint inertia matrix for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:58
% EndTime: 2020-01-03 11:36:07
% DurationCPUTime: 2.16s
% Computational Cost: add. (5523->256), mult. (4724->377), div. (0->0), fcn. (4266->10), ass. (0->116)
t127 = qJ(1) + pkin(8);
t150 = pkin(2) * cos(t127) + cos(qJ(1)) * pkin(1);
t123 = qJ(3) + t127;
t118 = sin(t123);
t119 = cos(t123);
t86 = t118 * rSges(4,1) + t119 * rSges(4,2);
t151 = pkin(2) * sin(t127) + sin(qJ(1)) * pkin(1);
t130 = sin(qJ(5));
t132 = cos(qJ(5));
t128 = sin(pkin(9));
t129 = cos(pkin(9));
t163 = Icges(6,4) * t132;
t79 = -Icges(6,6) * t129 + (-Icges(6,2) * t130 + t163) * t128;
t164 = Icges(6,4) * t130;
t80 = -Icges(6,5) * t129 + (Icges(6,1) * t132 - t164) * t128;
t174 = -t130 * t80 - t132 * t79;
t173 = 2 * m(4);
t172 = 2 * m(5);
t171 = 2 * m(6);
t126 = qJD(1) + qJD(3);
t157 = t126 * t119;
t158 = t126 * t118;
t170 = pkin(3) * t157 + qJ(4) * t158;
t149 = qJD(5) * t128;
t83 = (-Icges(6,2) * t132 - t164) * t149;
t168 = t130 * t83;
t82 = (-Icges(6,5) * t130 - Icges(6,6) * t132) * t149;
t84 = (-Icges(6,1) * t130 - t163) * t149;
t142 = t128 * t132 * t84 - t129 * t82;
t26 = (t174 * qJD(5) - t168) * t128 + t142;
t166 = t26 * t129;
t165 = qJ(4) * t157 + qJD(4) * t118;
t162 = t118 * t128;
t161 = t118 * t129;
t160 = t119 * t128;
t159 = t119 * t129;
t156 = t126 * t128;
t155 = t129 * t130;
t154 = t129 * t132;
t153 = t119 * pkin(3) + t118 * qJ(4);
t152 = t150 * qJD(1);
t147 = t119 * t156;
t75 = t118 * t154 - t119 * t130;
t76 = -t118 * t132 + t119 * t155;
t61 = -t75 * qJD(5) - t76 * t126;
t136 = t118 * t130 + t119 * t154;
t74 = -t118 * t155 - t119 * t132;
t62 = t74 * qJD(5) + t136 * t126;
t34 = t62 * rSges(6,1) + t61 * rSges(6,2) + rSges(6,3) * t147;
t53 = t75 * rSges(6,1) + t74 * rSges(6,2) + rSges(6,3) * t162;
t148 = t118 * t156;
t146 = t129 * t157;
t87 = t119 * rSges(4,1) - rSges(4,2) * t118;
t73 = rSges(4,1) * t157 - rSges(4,2) * t158;
t59 = t136 * qJD(5) + t74 * t126;
t60 = t76 * qJD(5) + t75 * t126;
t139 = -rSges(6,1) * t60 - t59 * rSges(6,2);
t54 = -rSges(6,1) * t136 + rSges(6,2) * t76 - rSges(6,3) * t160;
t135 = t151 * qJD(1);
t72 = t86 * t126;
t64 = rSges(5,1) * t159 - rSges(5,2) * t160 + t118 * rSges(5,3) + t153;
t114 = t118 * pkin(3);
t39 = pkin(4) * t161 + pkin(7) * t162 - qJ(4) * t119 + t114 + t53;
t23 = pkin(4) * t146 + pkin(7) * t147 - qJD(4) * t119 + t170 + t34;
t47 = Icges(6,5) * t75 + Icges(6,6) * t74 + Icges(6,3) * t162;
t49 = Icges(6,4) * t75 + Icges(6,2) * t74 + Icges(6,6) * t162;
t51 = Icges(6,1) * t75 + Icges(6,4) * t74 + Icges(6,5) * t162;
t18 = -t129 * t47 + (-t130 * t49 + t132 * t51) * t128;
t48 = -Icges(6,5) * t136 + Icges(6,6) * t76 - Icges(6,3) * t160;
t50 = -Icges(6,4) * t136 + Icges(6,2) * t76 - Icges(6,6) * t160;
t52 = -Icges(6,1) * t136 + Icges(6,4) * t76 - Icges(6,5) * t160;
t19 = -t129 * t48 + (-t130 * t50 + t132 * t52) * t128;
t28 = Icges(6,5) * t62 + Icges(6,6) * t61 + Icges(6,3) * t147;
t30 = Icges(6,4) * t62 + Icges(6,2) * t61 + Icges(6,6) * t147;
t32 = Icges(6,1) * t62 + Icges(6,4) * t61 + Icges(6,5) * t147;
t3 = -t129 * t28 + (-t130 * t30 + t132 * t32 + (-t130 * t51 - t132 * t49) * qJD(5)) * t128;
t78 = -Icges(6,3) * t129 + (Icges(6,5) * t132 - Icges(6,6) * t130) * t128;
t35 = t78 * t162 + t74 * t79 + t75 * t80;
t36 = -t136 * t80 - t78 * t160 + t76 * t79;
t27 = Icges(6,5) * t60 + Icges(6,6) * t59 + Icges(6,3) * t148;
t29 = Icges(6,4) * t60 + Icges(6,2) * t59 + Icges(6,6) * t148;
t31 = Icges(6,1) * t60 + Icges(6,4) * t59 + Icges(6,5) * t148;
t4 = -t129 * t27 + (-t130 * t29 + t132 * t31 + (-t130 * t52 - t132 * t50) * qJD(5)) * t128;
t8 = t59 * t79 + t60 * t80 + t76 * t83 - t136 * t84 + (-t119 * t82 + t78 * t158) * t128;
t9 = t61 * t79 + t62 * t80 + t74 * t83 + t75 * t84 + (t118 * t82 + t78 * t157) * t128;
t134 = -t166 + (t3 + t9) * t162 / 0.2e1 - (t4 + t8) * t160 / 0.2e1 + ((t19 + t36) * t118 + (t18 + t35) * t119) * t156 / 0.2e1;
t63 = -rSges(5,2) * t162 + rSges(5,1) * t161 + t114 + (-rSges(5,3) - qJ(4)) * t119;
t40 = pkin(4) * t159 + pkin(7) * t160 + t153 - t54;
t45 = rSges(5,2) * t148 + rSges(5,3) * t157 + (-rSges(5,1) * t129 - pkin(3)) * t158 + t165;
t46 = rSges(5,1) * t146 + rSges(5,3) * t158 + (-rSges(5,2) * t156 - qJD(4)) * t119 + t170;
t22 = (-pkin(4) * t129 - pkin(3) + (-rSges(6,3) - pkin(7)) * t128) * t158 + t139 + t165;
t85 = (-rSges(6,1) * t130 - rSges(6,2) * t132) * t149;
t81 = -rSges(6,3) * t129 + (rSges(6,1) * t132 - rSges(6,2) * t130) * t128;
t68 = t87 + t150;
t67 = t151 + t86;
t66 = t73 + t152;
t65 = -t72 - t135;
t56 = t64 + t150;
t55 = t63 + t151;
t44 = t46 + t152;
t43 = -t135 + t45;
t42 = t129 * t54 - t81 * t160;
t41 = -t129 * t53 - t81 * t162;
t38 = t40 + t150;
t37 = t39 + t151;
t33 = rSges(6,3) * t148 - t139;
t21 = -t129 * t34 + (-t118 * t85 - t81 * t157) * t128;
t20 = t129 * t33 + (-t119 * t85 + t81 * t158) * t128;
t17 = t23 + t152;
t16 = -t135 + t22;
t13 = -t136 * t52 - t48 * t160 + t50 * t76;
t12 = -t136 * t51 - t47 * t160 + t49 * t76;
t11 = t48 * t162 + t50 * t74 + t52 * t75;
t10 = t47 * t162 + t49 * t74 + t51 * t75;
t5 = ((t126 * t54 + t34) * t119 + (-t126 * t53 + t33) * t118) * t128;
t1 = [(t65 * t68 + t66 * t67) * t173 + (t43 * t56 + t44 * t55) * t172 - t128 * t168 + (t16 * t38 + t17 * t37) * t171 + t142 + t174 * t149; 0; 0; m(4) * (t65 * t87 + t66 * t86 + t67 * t73 - t68 * t72) + m(5) * (t43 * t64 + t44 * t63 + t45 * t56 + t46 * t55) + m(6) * (t16 * t40 + t17 * t39 + t22 * t38 + t23 * t37) + t26; 0; (t22 * t40 + t23 * t39) * t171 + (t45 * t64 + t46 * t63) * t172 + (-t72 * t87 + t73 * t86) * t173 + t26; m(5) * ((-t126 * t55 - t43) * t119 + (t126 * t56 - t44) * t118) + m(6) * ((-t126 * t37 - t16) * t119 + (t126 * t38 - t17) * t118); 0; m(6) * ((-t126 * t39 - t22) * t119 + (t126 * t40 - t23) * t118) + m(5) * ((-t126 * t63 - t45) * t119 + (t126 * t64 - t46) * t118); 0; m(6) * (t16 * t42 + t17 * t41 + t20 * t38 + t21 * t37) + t134; m(6) * t5; m(6) * (t20 * t40 + t21 * t39 + t22 * t42 + t23 * t41) + t134; m(6) * ((-t126 * t41 - t20) * t119 + (t126 * t42 - t21) * t118); (t42 * t20 + t41 * t21 + (t118 * t54 + t119 * t53) * t5 * t128) * t171 - t129 * (-t166 + ((t126 * t18 - t4) * t119 + (t126 * t19 + t3) * t118) * t128) + (-t129 * t35 + (t10 * t118 - t11 * t119) * t128) * t147 + (-t9 * t129 + ((t74 * t30 + t75 * t32 + t61 * t49 + t62 * t51) * t118 + t10 * t157 - (t74 * t29 + t75 * t31 + t61 * t50 + t62 * t52) * t119 + t11 * t158 + ((t118 * t28 + t47 * t157) * t118 - (t118 * t27 + t48 * t157) * t119) * t128) * t128) * t162 + (-t129 * t36 + (t118 * t12 - t119 * t13) * t128) * t148 - (-t8 * t129 + ((-t136 * t32 + t76 * t30 + t59 * t49 + t60 * t51) * t118 + t12 * t157 - (-t136 * t31 + t76 * t29 + t59 * t50 + t60 * t52) * t119 + t13 * t158 + ((-t119 * t28 + t47 * t158) * t118 - (-t119 * t27 + t48 * t158) * t119) * t128) * t128) * t160;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
