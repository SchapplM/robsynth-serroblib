% Calculate time derivative of joint inertia matrix for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR7_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR7_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:41
% DurationCPUTime: 2.33s
% Computational Cost: add. (2163->241), mult. (6305->411), div. (0->0), fcn. (5919->6), ass. (0->131)
t107 = sin(qJ(2));
t109 = cos(qJ(2));
t193 = ((Icges(4,5) - Icges(3,6)) * t109 + (Icges(4,4) - Icges(3,5)) * t107) * qJD(2);
t104 = sin(pkin(6));
t102 = t104 ^ 2;
t105 = cos(pkin(6));
t103 = t105 ^ 2;
t183 = t102 + t103;
t182 = t183 * qJD(2);
t192 = m(3) * (t107 * rSges(3,1) + rSges(3,2) * t109) * t182;
t191 = -t183 + 0.1e1;
t190 = pkin(5) * t183;
t189 = t193 * t104;
t188 = t193 * t105;
t179 = 2 * m(5);
t178 = m(4) / 0.2e1;
t177 = m(5) / 0.2e1;
t176 = t104 / 0.2e1;
t175 = -t105 / 0.2e1;
t99 = t107 * pkin(2) - qJ(3) * t109;
t174 = t183 * (-t99 * qJD(2) + qJD(3) * t107);
t140 = pkin(2) * t109 + qJ(3) * t107;
t173 = t183 * t140;
t142 = -rSges(4,2) * t109 + rSges(4,3) * t107;
t92 = t140 * qJD(2) - qJD(3) * t109;
t172 = -t142 * qJD(2) - t92;
t156 = t105 * t109;
t158 = t104 * t109;
t106 = sin(qJ(4));
t108 = cos(qJ(4));
t154 = t107 * t108;
t95 = t104 * t154 + t105 * t106;
t155 = t106 * t107;
t96 = t104 * t155 - t105 * t108;
t42 = Icges(5,5) * t96 + Icges(5,6) * t95 + Icges(5,3) * t158;
t44 = Icges(5,4) * t96 + Icges(5,2) * t95 + Icges(5,6) * t158;
t46 = Icges(5,1) * t96 + Icges(5,4) * t95 + Icges(5,5) * t158;
t93 = -t104 * t106 + t105 * t154;
t94 = t104 * t108 + t105 * t155;
t16 = t42 * t156 + t93 * t44 + t94 * t46;
t170 = t104 * t16;
t41 = Icges(5,5) * t94 + Icges(5,6) * t93 + Icges(5,3) * t156;
t43 = Icges(5,4) * t94 + Icges(5,2) * t93 + Icges(5,6) * t156;
t45 = Icges(5,1) * t94 + Icges(5,4) * t93 + Icges(5,5) * t156;
t17 = t41 * t158 + t95 * t43 + t96 * t45;
t169 = t105 * t17;
t124 = Icges(5,5) * t106 + Icges(5,6) * t108;
t150 = qJD(4) * t109;
t168 = t107 * ((-Icges(5,5) * t108 + Icges(5,6) * t106) * t150 + (Icges(5,3) * t109 + t124 * t107) * qJD(2));
t74 = Icges(5,3) * t107 - t124 * t109;
t167 = t107 * t74;
t141 = t107 * rSges(4,2) + rSges(4,3) * t109;
t166 = t141 - t99;
t163 = Icges(5,4) * t106;
t162 = Icges(5,4) * t108;
t159 = t104 * t107;
t157 = t105 * t107;
t152 = qJD(2) * t107;
t151 = qJD(2) * t109;
t149 = t104 * t152;
t148 = t105 * t152;
t147 = t106 * t151;
t146 = t108 * t151;
t143 = rSges(5,1) * t106 + rSges(5,2) * t108;
t77 = t107 * rSges(5,3) - t143 * t109;
t145 = -pkin(5) * t107 - t77 - t99;
t137 = t106 * t45 + t108 * t43;
t20 = t107 * t41 - t137 * t109;
t136 = t106 * t46 + t108 * t44;
t21 = t107 * t42 - t136 * t109;
t139 = t21 * t104 + t20 * t105;
t49 = t94 * rSges(5,1) + t93 * rSges(5,2) + rSges(5,3) * t156;
t50 = t96 * rSges(5,1) + t95 * rSges(5,2) + rSges(5,3) * t158;
t138 = t104 * t49 - t105 * t50;
t125 = Icges(5,2) * t108 + t163;
t75 = Icges(5,6) * t107 - t125 * t109;
t127 = Icges(5,1) * t106 + t162;
t76 = Icges(5,5) * t107 - t127 * t109;
t135 = t106 * t76 + t108 * t75;
t54 = (-rSges(5,1) * t108 + rSges(5,2) * t106) * t150 + (rSges(5,3) * t109 + t143 * t107) * qJD(2);
t121 = -pkin(5) * t151 - t54 - t92;
t59 = t93 * qJD(4) + t105 * t147;
t60 = -t94 * qJD(4) + t105 * t146;
t32 = Icges(5,5) * t59 + Icges(5,6) * t60 - Icges(5,3) * t148;
t120 = t109 * t32 - t41 * t152;
t61 = t95 * qJD(4) + t104 * t147;
t62 = -t96 * qJD(4) + t104 * t146;
t33 = Icges(5,5) * t61 + Icges(5,6) * t62 - Icges(5,3) * t149;
t119 = t109 * t33 - t42 * t152;
t65 = t166 * t105;
t64 = t166 * t104;
t56 = t172 * t105;
t55 = t172 * t104;
t53 = (-Icges(5,1) * t108 + t163) * t150 + (Icges(5,5) * t109 + t127 * t107) * qJD(2);
t52 = (Icges(5,2) * t106 - t162) * t150 + (Icges(5,6) * t109 + t125 * t107) * qJD(2);
t48 = t145 * t105;
t47 = t145 * t104;
t39 = rSges(5,1) * t61 + rSges(5,2) * t62 - rSges(5,3) * t149;
t38 = rSges(5,1) * t59 + rSges(5,2) * t60 - rSges(5,3) * t148;
t37 = Icges(5,1) * t61 + Icges(5,4) * t62 - Icges(5,5) * t149;
t36 = Icges(5,1) * t59 + Icges(5,4) * t60 - Icges(5,5) * t148;
t35 = Icges(5,4) * t61 + Icges(5,2) * t62 - Icges(5,6) * t149;
t34 = Icges(5,4) * t59 + Icges(5,2) * t60 - Icges(5,6) * t148;
t31 = t107 * t49 - t77 * t156;
t30 = -t107 * t50 + t77 * t158;
t29 = t121 * t105;
t28 = t121 * t104;
t27 = -t135 * t109 + t167;
t26 = t183 * t142 + t173;
t25 = t74 * t158 + t95 * t75 + t96 * t76;
t24 = t74 * t156 + t93 * t75 + t94 * t76;
t23 = t138 * t109;
t22 = t141 * t182 + t174;
t19 = t104 * t50 + t105 * t49 + t109 * t190 + t173;
t18 = t42 * t158 + t95 * t44 + t96 * t46;
t15 = t41 * t156 + t93 * t43 + t94 * t45;
t14 = -t54 * t156 + t107 * t38 + (t109 * t49 + t77 * t157) * qJD(2);
t13 = t54 * t158 - t107 * t39 + (-t109 * t50 - t77 * t159) * qJD(2);
t12 = t104 * t39 + t105 * t38 - t152 * t190 + t174;
t11 = (-t104 * t38 + t105 * t39) * t109 + t138 * t152;
t10 = t119 * t104 + t95 * t35 + t96 * t37 + t62 * t44 + t61 * t46;
t9 = t120 * t104 + t95 * t34 + t96 * t36 + t62 * t43 + t61 * t45;
t8 = t119 * t105 + t93 * t35 + t94 * t37 + t60 * t44 + t59 * t46;
t7 = t120 * t105 + t93 * t34 + t94 * t36 + t60 * t43 + t59 * t45;
t6 = (t136 * qJD(2) + t33) * t107 + (qJD(2) * t42 - t106 * t37 - t108 * t35 + (t106 * t44 - t108 * t46) * qJD(4)) * t109;
t5 = (t137 * qJD(2) + t32) * t107 + (qJD(2) * t41 - t106 * t36 - t108 * t34 + (t106 * t43 - t108 * t45) * qJD(4)) * t109;
t4 = -t10 * t105 + t104 * t9;
t3 = t104 * t7 - t105 * t8;
t2 = (t95 * t52 + t96 * t53 + t61 * t76 + t62 * t75) * t107 + (t9 * t105 + (t10 + t168) * t104) * t109 + (t25 * t109 + (-t169 + (-t18 - t167) * t104) * t107) * qJD(2);
t1 = (t93 * t52 + t94 * t53 + t59 * t76 + t60 * t75) * t107 + (t8 * t104 + (t7 + t168) * t105) * t109 + (t24 * t109 + (-t170 + (-t15 - t167) * t105) * t107) * qJD(2);
t40 = [0; m(4) * t22 + m(5) * t12 - t192; 0.2e1 * m(4) * (t26 * t22 + t55 * t64 + t56 * t65) + (t12 * t19 + t28 * t47 + t29 * t48) * t179 + (-t189 * t103 - t4) * t105 + (t3 + t188 * t102 + (-t189 * t104 + t188 * t105) * t105) * t104 + 0.2e1 * t191 * (rSges(3,1) * t109 - rSges(3,2) * t107) * t192; (m(4) + m(5)) * t152; m(4) * (-t109 * t22 + t56 * t157 + t55 * t159) + m(5) * (-t109 * t12 + t29 * t157 + t28 * t159) + 0.2e1 * ((t107 * t26 + t65 * t156 + t64 * t158) * t178 + (t107 * t19 + t48 * t156 + t47 * t158) * t177) * qJD(2); -0.4e1 * (t178 + t177) * t191 * t107 * t151; m(5) * t11; t107 * (t104 * t5 - t105 * t6) / 0.2e1 + t1 * t176 + t2 * t175 + m(5) * (t11 * t19 - t12 * t23 + t13 * t48 + t14 * t47 + t28 * t31 + t29 * t30) + (t105 * t3 / 0.2e1 + t4 * t176) * t109 + (t109 * (t104 * t20 - t105 * t21) / 0.2e1 + ((t104 * t15 - t105 * t16) * t175 - t104 * (t104 * t17 - t105 * t18) / 0.2e1) * t107) * qJD(2); m(5) * (-t109 * t11 + (t104 * t14 + t105 * t13) * t107 + (-t107 * t23 + (t104 * t31 + t105 * t30) * t109) * qJD(2)); (-t11 * t23 + t13 * t30 + t14 * t31) * t179 - (t24 * t107 + (t105 * t15 + t170) * t109) * t148 + t1 * t156 - (t25 * t107 + (t104 * t18 + t169) * t109) * t149 + t2 * t158 + (t27 * t107 + t139 * t109) * t151 + t107 * ((t168 + (t135 * t107 - t139) * qJD(2)) * t107 + (t5 * t105 + t6 * t104 + (qJD(2) * t74 - t106 * t53 - t108 * t52 + (t106 * t75 - t108 * t76) * qJD(4)) * t107 + t27 * qJD(2)) * t109);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t40(1), t40(2), t40(4), t40(7); t40(2), t40(3), t40(5), t40(8); t40(4), t40(5), t40(6), t40(9); t40(7), t40(8), t40(9), t40(10);];
Mq = res;
