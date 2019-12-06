% Calculate time derivative of joint inertia matrix for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:19
% EndTime: 2019-12-05 15:26:28
% DurationCPUTime: 3.48s
% Computational Cost: add. (4467->266), mult. (7001->451), div. (0->0), fcn. (6457->8), ass. (0->151)
t131 = sin(qJ(2));
t133 = cos(qJ(2));
t127 = sin(pkin(7));
t124 = t127 ^ 2;
t128 = cos(pkin(7));
t125 = t128 ^ 2;
t248 = t125 + t124;
t247 = qJD(2) * t248;
t251 = m(3) * (rSges(3,1) * t131 + rSges(3,2) * t133) * t247;
t250 = 0.1e1 - t248;
t126 = qJ(2) + pkin(8);
t122 = sin(t126);
t123 = cos(t126);
t249 = (-Icges(3,5) * t131 - Icges(3,6) * t133 + (Icges(5,5) - Icges(4,6)) * t123 + (Icges(5,4) - Icges(4,5)) * t122) * qJD(2);
t246 = t249 * t127;
t245 = t249 * t128;
t231 = 2 * m(6);
t230 = m(5) / 0.2e1;
t229 = m(6) / 0.2e1;
t228 = t127 / 0.2e1;
t227 = -t128 / 0.2e1;
t226 = pkin(2) * t131;
t225 = pkin(2) * t133;
t224 = pkin(6) * t123;
t221 = t248 * t225;
t220 = pkin(2) * qJD(2);
t130 = sin(qJ(5));
t132 = cos(qJ(5));
t157 = Icges(6,5) * t130 + Icges(6,6) * t132;
t194 = qJD(5) * t123;
t219 = t122 * ((-Icges(6,5) * t132 + Icges(6,6) * t130) * t194 + (Icges(6,3) * t123 + t122 * t157) * qJD(2));
t71 = Icges(6,3) * t122 - t123 * t157;
t218 = t122 * t71;
t199 = t128 * t132;
t202 = t127 * t130;
t108 = t122 * t199 - t202;
t200 = t128 * t130;
t201 = t127 * t132;
t109 = t122 * t200 + t201;
t203 = t123 * t128;
t110 = t122 * t201 + t200;
t111 = t122 * t202 - t199;
t204 = t123 * t127;
t50 = Icges(6,5) * t111 + Icges(6,6) * t110 + Icges(6,3) * t204;
t52 = Icges(6,4) * t111 + Icges(6,2) * t110 + Icges(6,6) * t204;
t54 = Icges(6,1) * t111 + Icges(6,4) * t110 + Icges(6,5) * t204;
t17 = t108 * t52 + t109 * t54 + t203 * t50;
t217 = t127 * t17;
t49 = Icges(6,5) * t109 + Icges(6,6) * t108 + Icges(6,3) * t203;
t51 = Icges(6,4) * t109 + Icges(6,2) * t108 + Icges(6,6) * t203;
t53 = Icges(6,1) * t109 + Icges(6,4) * t108 + Icges(6,5) * t203;
t18 = t110 * t51 + t111 * t53 + t204 * t49;
t216 = t128 * t18;
t211 = Icges(6,4) * t130;
t210 = Icges(6,4) * t132;
t206 = t122 * t127;
t205 = t122 * t128;
t198 = t248 * t131 * t220;
t196 = qJD(2) * t122;
t195 = qJD(2) * t123;
t192 = t133 * t220;
t191 = t127 * t196;
t190 = t128 * t196;
t189 = t130 * t195;
t188 = t132 * t195;
t116 = pkin(3) * t122 - qJ(4) * t123;
t187 = -t116 - t226;
t118 = rSges(4,1) * t122 + rSges(4,2) * t123;
t186 = -t118 - t226;
t176 = pkin(3) * t123 + qJ(4) * t122;
t185 = t176 * t248 + t221;
t184 = -t198 + t248 * (-qJD(2) * t116 + qJD(4) * t122);
t177 = rSges(5,2) * t122 + rSges(5,3) * t123;
t183 = t177 + t187;
t180 = rSges(4,1) * t123 - rSges(4,2) * t122;
t182 = -qJD(2) * t180 - t192;
t179 = rSges(6,1) * t130 + rSges(6,2) * t132;
t178 = -rSges(5,2) * t123 + rSges(5,3) * t122;
t167 = t130 * t53 + t132 * t51;
t20 = t122 * t49 - t123 * t167;
t166 = t130 * t54 + t132 * t52;
t21 = t122 * t50 - t123 * t166;
t169 = t21 * t127 + t20 * t128;
t55 = rSges(6,1) * t109 + rSges(6,2) * t108 + rSges(6,3) * t203;
t56 = rSges(6,1) * t111 + rSges(6,2) * t110 + rSges(6,3) * t204;
t168 = t127 * t55 - t128 * t56;
t158 = Icges(6,2) * t132 + t211;
t72 = Icges(6,6) * t122 - t123 * t158;
t161 = Icges(6,1) * t130 + t210;
t73 = Icges(6,5) * t122 - t123 * t161;
t165 = t130 * t73 + t132 * t72;
t83 = qJD(2) * t176 - qJD(4) * t123;
t152 = -qJD(2) * t178 - t192 - t83;
t63 = qJD(5) * t108 + t128 * t189;
t64 = -qJD(5) * t109 + t128 * t188;
t33 = Icges(6,5) * t63 + Icges(6,6) * t64 - Icges(6,3) * t190;
t151 = t123 * t33 - t196 * t49;
t65 = qJD(5) * t110 + t127 * t189;
t66 = -qJD(5) * t111 + t127 * t188;
t34 = Icges(6,5) * t65 + Icges(6,6) * t66 - Icges(6,3) * t191;
t150 = t123 * t34 - t196 * t50;
t32 = -t118 * t247 - t198;
t149 = t180 * t32;
t74 = rSges(6,3) * t122 - t123 * t179;
t146 = -pkin(6) * t122 + t187 - t74;
t48 = (-rSges(6,1) * t132 + rSges(6,2) * t130) * t194 + (rSges(6,3) * t123 + t122 * t179) * qJD(2);
t135 = -t48 - t83 + (-t224 - t225) * qJD(2);
t68 = t182 * t128;
t67 = t182 * t127;
t59 = t183 * t128;
t58 = t183 * t127;
t47 = (-Icges(6,1) * t132 + t211) * t194 + (Icges(6,5) * t123 + t122 * t161) * qJD(2);
t46 = (Icges(6,2) * t130 - t210) * t194 + (Icges(6,6) * t123 + t122 * t158) * qJD(2);
t44 = t152 * t128;
t43 = t152 * t127;
t42 = t146 * t128;
t41 = t146 * t127;
t40 = rSges(6,1) * t65 + rSges(6,2) * t66 - rSges(6,3) * t191;
t39 = rSges(6,1) * t63 + rSges(6,2) * t64 - rSges(6,3) * t190;
t38 = Icges(6,1) * t65 + Icges(6,4) * t66 - Icges(6,5) * t191;
t37 = Icges(6,1) * t63 + Icges(6,4) * t64 - Icges(6,5) * t190;
t36 = Icges(6,4) * t65 + Icges(6,2) * t66 - Icges(6,6) * t191;
t35 = Icges(6,4) * t63 + Icges(6,2) * t64 - Icges(6,6) * t190;
t31 = t122 * t55 - t203 * t74;
t30 = -t122 * t56 + t204 * t74;
t29 = -t123 * t165 + t218;
t28 = t135 * t128;
t27 = t135 * t127;
t26 = t168 * t123;
t25 = t110 * t72 + t111 * t73 + t204 * t71;
t24 = t108 * t72 + t109 * t73 + t203 * t71;
t23 = t177 * t247 + t184;
t22 = t178 * t248 + t185;
t19 = t110 * t52 + t111 * t54 + t204 * t50;
t16 = t108 * t51 + t109 * t53 + t203 * t49;
t15 = -t48 * t203 + t122 * t39 + (t123 * t55 + t205 * t74) * qJD(2);
t14 = t48 * t204 - t122 * t40 + (-t123 * t56 - t206 * t74) * qJD(2);
t13 = t127 * t56 + t128 * t55 + t224 * t248 + t185;
t12 = (-t127 * t39 + t128 * t40) * t123 + t168 * t196;
t11 = -pkin(6) * t196 * t248 + t127 * t40 + t128 * t39 + t184;
t10 = t110 * t36 + t111 * t38 + t127 * t150 + t52 * t66 + t54 * t65;
t9 = t110 * t35 + t111 * t37 + t127 * t151 + t51 * t66 + t53 * t65;
t8 = t108 * t36 + t109 * t38 + t128 * t150 + t52 * t64 + t54 * t63;
t7 = t108 * t35 + t109 * t37 + t128 * t151 + t51 * t64 + t53 * t63;
t6 = (qJD(2) * t166 + t34) * t122 + (qJD(2) * t50 - t130 * t38 - t132 * t36 + (t130 * t52 - t132 * t54) * qJD(5)) * t123;
t5 = (qJD(2) * t167 + t33) * t122 + (qJD(2) * t49 - t130 * t37 - t132 * t35 + (t130 * t51 - t132 * t53) * qJD(5)) * t123;
t4 = -t10 * t128 + t127 * t9;
t3 = t127 * t7 - t128 * t8;
t2 = (t110 * t46 + t111 * t47 + t65 * t73 + t66 * t72) * t122 + (t9 * t128 + (t10 + t219) * t127) * t123 + (t25 * t123 + (-t216 + (-t19 - t218) * t127) * t122) * qJD(2);
t1 = (t108 * t46 + t109 * t47 + t63 * t73 + t64 * t72) * t122 + (t8 * t127 + (t7 + t219) * t128) * t123 + (t24 * t123 + (-t217 + (-t16 - t218) * t128) * t122) * qJD(2);
t45 = [0; m(4) * t32 + m(5) * t23 + m(6) * t11 - t251; -t128 * t4 + (t11 * t13 + t27 * t41 + t28 * t42) * t231 + 0.2e1 * m(5) * (t22 * t23 + t43 * t58 + t44 * t59) + 0.2e1 * m(4) * (t221 * t32 + (t128 * t149 + t186 * t68) * t128) - t246 * t128 * t125 + (t3 + 0.2e1 * (t127 * t149 + t186 * t67) * m(4) + t245 * t124 + (-t127 * t246 + t128 * t245) * t128) * t127 + 0.2e1 * t250 * (rSges(3,1) * t133 - rSges(3,2) * t131) * t251; 0; m(6) * (t127 * t28 - t128 * t27) + m(5) * (t127 * t44 - t128 * t43) + m(4) * (t127 * t68 - t128 * t67); 0; (m(5) + m(6)) * t196; m(6) * (-t11 * t123 + t205 * t28 + t206 * t27) + m(5) * (-t123 * t23 + t205 * t44 + t206 * t43) + 0.2e1 * ((t122 * t13 + t203 * t42 + t204 * t41) * t229 + (t122 * t22 + t203 * t59 + t204 * t58) * t230) * qJD(2); 0; -0.4e1 * (t230 + t229) * t250 * t122 * t195; m(6) * t12; t2 * t227 + m(6) * (-t11 * t26 + t12 * t13 + t14 * t42 + t15 * t41 + t27 * t31 + t28 * t30) + t122 * (t127 * t5 - t128 * t6) / 0.2e1 + t1 * t228 + (t128 * t3 / 0.2e1 + t4 * t228) * t123 + (t123 * (t127 * t20 - t128 * t21) / 0.2e1 + ((t127 * t16 - t128 * t17) * t227 - t127 * (t18 * t127 - t19 * t128) / 0.2e1) * t122) * qJD(2); m(6) * (t127 * t14 - t128 * t15); m(6) * (-t12 * t123 + (t127 * t15 + t128 * t14) * t122 + (-t122 * t26 + (t127 * t31 + t128 * t30) * t123) * qJD(2)); (-t12 * t26 + t14 * t30 + t15 * t31) * t231 - (t122 * t24 + (t128 * t16 + t217) * t123) * t190 + t1 * t203 - (t122 * t25 + (t127 * t19 + t216) * t123) * t191 + t2 * t204 + (t29 * t122 + t123 * t169) * t195 + t122 * ((t219 + (t122 * t165 - t169) * qJD(2)) * t122 + (t5 * t128 + t6 * t127 + (qJD(2) * t71 - t130 * t47 - t132 * t46 + (t130 * t72 - t132 * t73) * qJD(5)) * t122 + t29 * qJD(2)) * t123);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t45(1), t45(2), t45(4), t45(7), t45(11); t45(2), t45(3), t45(5), t45(8), t45(12); t45(4), t45(5), t45(6), t45(9), t45(13); t45(7), t45(8), t45(9), t45(10), t45(14); t45(11), t45(12), t45(13), t45(14), t45(15);];
Mq = res;
