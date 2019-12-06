% Calculate time derivative of joint inertia matrix for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:00
% EndTime: 2019-12-05 15:03:06
% DurationCPUTime: 2.60s
% Computational Cost: add. (4254->246), mult. (6509->418), div. (0->0), fcn. (6071->6), ass. (0->134)
t106 = pkin(8) + qJ(3);
t102 = sin(t106);
t103 = cos(t106);
t196 = ((Icges(5,5) - Icges(4,6)) * t103 + (Icges(5,4) - Icges(4,5)) * t102) * qJD(3);
t107 = sin(pkin(7));
t104 = t107 ^ 2;
t108 = cos(pkin(7));
t105 = t108 ^ 2;
t186 = t104 + t105;
t185 = t186 * qJD(3);
t195 = m(4) * (rSges(4,1) * t102 + rSges(4,2) * t103) * t185;
t194 = -t186 + 0.1e1;
t193 = pkin(6) * t186;
t192 = t196 * t107;
t191 = t196 * t108;
t182 = 2 * m(6);
t181 = m(5) / 0.2e1;
t180 = m(6) / 0.2e1;
t179 = t107 / 0.2e1;
t178 = -t108 / 0.2e1;
t99 = pkin(3) * t102 - qJ(4) * t103;
t177 = t186 * (-t99 * qJD(3) + qJD(4) * t102);
t143 = -rSges(5,2) * t103 + rSges(5,3) * t102;
t141 = pkin(3) * t103 + qJ(4) * t102;
t78 = t141 * qJD(3) - qJD(4) * t103;
t176 = -t143 * qJD(3) - t78;
t175 = t186 * t141;
t109 = sin(qJ(5));
t110 = cos(qJ(5));
t125 = Icges(6,5) * t109 + Icges(6,6) * t110;
t151 = qJD(5) * t103;
t173 = t102 * ((-Icges(6,5) * t110 + Icges(6,6) * t109) * t151 + (Icges(6,3) * t103 + t125 * t102) * qJD(3));
t66 = Icges(6,3) * t102 - t125 * t103;
t172 = t102 * t66;
t159 = t103 * t108;
t160 = t103 * t107;
t156 = t108 * t109;
t157 = t107 * t110;
t95 = t102 * t157 + t156;
t155 = t108 * t110;
t158 = t107 * t109;
t96 = t102 * t158 - t155;
t48 = Icges(6,5) * t96 + Icges(6,6) * t95 + Icges(6,3) * t160;
t50 = Icges(6,4) * t96 + Icges(6,2) * t95 + Icges(6,6) * t160;
t52 = Icges(6,1) * t96 + Icges(6,4) * t95 + Icges(6,5) * t160;
t93 = t102 * t155 - t158;
t94 = t102 * t156 + t157;
t16 = t48 * t159 + t50 * t93 + t52 * t94;
t171 = t107 * t16;
t47 = Icges(6,5) * t94 + Icges(6,6) * t93 + Icges(6,3) * t159;
t49 = Icges(6,4) * t94 + Icges(6,2) * t93 + Icges(6,6) * t159;
t51 = Icges(6,1) * t94 + Icges(6,4) * t93 + Icges(6,5) * t159;
t17 = t47 * t160 + t49 * t95 + t51 * t96;
t170 = t108 * t17;
t142 = rSges(5,2) * t102 + rSges(5,3) * t103;
t169 = t142 - t99;
t166 = Icges(6,4) * t109;
t165 = Icges(6,4) * t110;
t162 = t102 * t107;
t161 = t102 * t108;
t153 = qJD(3) * t102;
t152 = qJD(3) * t103;
t150 = t107 * t153;
t149 = t108 * t153;
t148 = t109 * t152;
t147 = t110 * t152;
t144 = rSges(6,1) * t109 + rSges(6,2) * t110;
t69 = t102 * rSges(6,3) - t144 * t103;
t146 = -pkin(6) * t102 - t69 - t99;
t132 = t109 * t51 + t110 * t49;
t20 = t102 * t47 - t132 * t103;
t131 = t109 * t52 + t110 * t50;
t21 = t102 * t48 - t131 * t103;
t134 = t21 * t107 + t20 * t108;
t53 = rSges(6,1) * t94 + rSges(6,2) * t93 + rSges(6,3) * t159;
t54 = rSges(6,1) * t96 + rSges(6,2) * t95 + rSges(6,3) * t160;
t133 = t107 * t53 - t108 * t54;
t126 = Icges(6,2) * t110 + t166;
t67 = Icges(6,6) * t102 - t126 * t103;
t128 = Icges(6,1) * t109 + t165;
t68 = Icges(6,5) * t102 - t128 * t103;
t130 = t109 * t68 + t110 * t67;
t46 = (-rSges(6,1) * t110 + rSges(6,2) * t109) * t151 + (rSges(6,3) * t103 + t144 * t102) * qJD(3);
t122 = -pkin(6) * t152 - t46 - t78;
t62 = t93 * qJD(5) + t108 * t148;
t63 = -t94 * qJD(5) + t108 * t147;
t32 = Icges(6,5) * t62 + Icges(6,6) * t63 - Icges(6,3) * t149;
t121 = t103 * t32 - t47 * t153;
t64 = t95 * qJD(5) + t107 * t148;
t65 = -t96 * qJD(5) + t107 * t147;
t33 = Icges(6,5) * t64 + Icges(6,6) * t65 - Icges(6,3) * t150;
t120 = t103 * t33 - t48 * t153;
t61 = t169 * t108;
t60 = t169 * t107;
t56 = t176 * t108;
t55 = t176 * t107;
t45 = (-Icges(6,1) * t110 + t166) * t151 + (Icges(6,5) * t103 + t128 * t102) * qJD(3);
t44 = (Icges(6,2) * t109 - t165) * t151 + (Icges(6,6) * t103 + t126 * t102) * qJD(3);
t41 = t146 * t108;
t40 = t146 * t107;
t39 = rSges(6,1) * t64 + rSges(6,2) * t65 - rSges(6,3) * t150;
t38 = rSges(6,1) * t62 + rSges(6,2) * t63 - rSges(6,3) * t149;
t37 = Icges(6,1) * t64 + Icges(6,4) * t65 - Icges(6,5) * t150;
t36 = Icges(6,1) * t62 + Icges(6,4) * t63 - Icges(6,5) * t149;
t35 = Icges(6,4) * t64 + Icges(6,2) * t65 - Icges(6,6) * t150;
t34 = Icges(6,4) * t62 + Icges(6,2) * t63 - Icges(6,6) * t149;
t31 = t102 * t53 - t69 * t159;
t30 = -t102 * t54 + t69 * t160;
t29 = t122 * t108;
t28 = t122 * t107;
t27 = -t130 * t103 + t172;
t26 = t186 * t143 + t175;
t25 = t133 * t103;
t24 = t66 * t160 + t67 * t95 + t68 * t96;
t23 = t66 * t159 + t67 * t93 + t68 * t94;
t22 = t142 * t185 + t177;
t19 = t103 * t193 + t107 * t54 + t108 * t53 + t175;
t18 = t48 * t160 + t50 * t95 + t52 * t96;
t15 = t47 * t159 + t49 * t93 + t51 * t94;
t14 = -t46 * t159 + t102 * t38 + (t103 * t53 + t69 * t161) * qJD(3);
t13 = t46 * t160 - t102 * t39 + (-t103 * t54 - t69 * t162) * qJD(3);
t12 = t107 * t39 + t108 * t38 - t153 * t193 + t177;
t11 = (-t107 * t38 + t108 * t39) * t103 + t133 * t153;
t10 = t120 * t107 + t35 * t95 + t37 * t96 + t50 * t65 + t52 * t64;
t9 = t121 * t107 + t34 * t95 + t36 * t96 + t49 * t65 + t51 * t64;
t8 = t120 * t108 + t35 * t93 + t37 * t94 + t50 * t63 + t52 * t62;
t7 = t121 * t108 + t34 * t93 + t36 * t94 + t49 * t63 + t51 * t62;
t6 = (t131 * qJD(3) + t33) * t102 + (qJD(3) * t48 - t109 * t37 - t110 * t35 + (t109 * t50 - t110 * t52) * qJD(5)) * t103;
t5 = (t132 * qJD(3) + t32) * t102 + (qJD(3) * t47 - t109 * t36 - t110 * t34 + (t109 * t49 - t110 * t51) * qJD(5)) * t103;
t4 = -t10 * t108 + t107 * t9;
t3 = t107 * t7 - t108 * t8;
t2 = (t44 * t95 + t45 * t96 + t64 * t68 + t65 * t67) * t102 + (t9 * t108 + (t10 + t173) * t107) * t103 + (t24 * t103 + (-t170 + (-t18 - t172) * t107) * t102) * qJD(3);
t1 = (t44 * t93 + t45 * t94 + t62 * t68 + t63 * t67) * t102 + (t8 * t107 + (t7 + t173) * t108) * t103 + (t23 * t103 + (-t171 + (-t15 - t172) * t108) * t102) * qJD(3);
t42 = [0; 0; 0; m(5) * t22 + m(6) * t12 - t195; m(5) * (t107 * t56 - t108 * t55) + m(6) * (t107 * t29 - t108 * t28); 0.2e1 * m(5) * (t22 * t26 + t55 * t60 + t56 * t61) + (t12 * t19 + t28 * t40 + t29 * t41) * t182 + (-t192 * t105 - t4) * t108 + (t3 + t191 * t104 + (-t192 * t107 + t191 * t108) * t108) * t107 + 0.2e1 * t194 * (rSges(4,1) * t103 - rSges(4,2) * t102) * t195; (m(5) + m(6)) * t153; 0; m(5) * (-t103 * t22 + t56 * t161 + t55 * t162) + m(6) * (-t103 * t12 + t29 * t161 + t28 * t162) + 0.2e1 * ((t102 * t26 + t61 * t159 + t60 * t160) * t181 + (t102 * t19 + t41 * t159 + t40 * t160) * t180) * qJD(3); -0.4e1 * (t181 + t180) * t194 * t102 * t152; m(6) * t11; m(6) * (t107 * t13 - t108 * t14); t102 * (t107 * t5 - t108 * t6) / 0.2e1 + t1 * t179 + t2 * t178 + m(6) * (t11 * t19 - t12 * t25 + t13 * t41 + t14 * t40 + t28 * t31 + t29 * t30) + (t108 * t3 / 0.2e1 + t4 * t179) * t103 + (t103 * (t107 * t20 - t108 * t21) / 0.2e1 + ((t107 * t15 - t108 * t16) * t178 - t107 * (t107 * t17 - t108 * t18) / 0.2e1) * t102) * qJD(3); m(6) * (-t103 * t11 + (t107 * t14 + t108 * t13) * t102 + (-t102 * t25 + (t107 * t31 + t108 * t30) * t103) * qJD(3)); (-t11 * t25 + t13 * t30 + t14 * t31) * t182 - (t102 * t23 + (t108 * t15 + t171) * t103) * t149 + t1 * t159 - (t102 * t24 + (t107 * t18 + t170) * t103) * t150 + t2 * t160 + (t102 * t27 + t134 * t103) * t152 + t102 * ((t173 + (t130 * t102 - t134) * qJD(3)) * t102 + (t5 * t108 + t6 * t107 + (qJD(3) * t66 - t109 * t45 - t110 * t44 + (t109 * t67 - t110 * t68) * qJD(5)) * t102 + t27 * qJD(3)) * t103);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t42(1), t42(2), t42(4), t42(7), t42(11); t42(2), t42(3), t42(5), t42(8), t42(12); t42(4), t42(5), t42(6), t42(9), t42(13); t42(7), t42(8), t42(9), t42(10), t42(14); t42(11), t42(12), t42(13), t42(14), t42(15);];
Mq = res;
