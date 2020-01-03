% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP11_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:32
% EndTime: 2019-12-31 18:54:40
% DurationCPUTime: 2.65s
% Computational Cost: add. (3658->304), mult. (9672->382), div. (0->0), fcn. (7040->6), ass. (0->146)
t106 = cos(pkin(8));
t100 = -t106 * pkin(2) - pkin(1);
t93 = t100 * qJD(1) + qJD(2);
t205 = qJD(2) + t93;
t105 = sin(pkin(8));
t108 = sin(qJ(3));
t184 = cos(qJ(3));
t122 = -t108 * t105 + t184 * t106;
t204 = t122 * qJD(1);
t113 = qJD(3) * t204;
t78 = qJD(4) - t204;
t203 = qJD(4) * qJD(3) + t113;
t107 = sin(qJ(4));
t109 = cos(qJ(4));
t155 = qJD(4) * t109;
t156 = qJD(4) * t107;
t115 = t122 * qJD(2);
t169 = pkin(6) + qJ(2);
t95 = t169 * t105;
t91 = qJD(1) * t95;
t96 = t169 * t106;
t92 = qJD(1) * t96;
t197 = -t108 * t92 - t184 * t91;
t34 = qJD(1) * t115 + qJD(3) * t197;
t90 = t184 * t105 + t108 * t106;
t84 = t90 * qJD(1);
t44 = -pkin(3) * t204 - t84 * pkin(7) + t93;
t87 = t90 * qJD(3);
t77 = qJD(1) * t87;
t49 = t77 * pkin(3) - pkin(7) * t113;
t63 = -t108 * t91 + t184 * t92;
t58 = qJD(3) * pkin(7) + t63;
t146 = t107 * t34 - t109 * t49 + t58 * t155 + t44 * t156;
t21 = t107 * t44 + t109 * t58;
t181 = t21 * t78;
t199 = -t146 + t181;
t42 = t203 * t107 + t84 * t155;
t161 = t42 * t109;
t41 = -t203 * t109 + t84 * t156;
t162 = t41 * t107;
t69 = -t109 * qJD(3) + t107 * t84;
t71 = t107 * qJD(3) + t109 * t84;
t86 = t122 * qJD(3);
t202 = ((t107 * t69 - t109 * t71) * qJD(4) - t161 + t162) * t90 - (t107 * t71 + t109 * t69) * t86;
t16 = t78 * qJ(5) + t21;
t183 = t16 * t78;
t185 = t77 * pkin(4);
t2 = t146 - t185;
t201 = -t2 + t183;
t20 = -t107 * t58 + t109 * t44;
t5 = t107 * t49 + t109 * t34 + t44 * t155 - t58 * t156;
t200 = -t20 * t78 + t5;
t196 = -t108 * t96 - t184 * t95;
t195 = t84 * qJD(3);
t186 = pkin(7) * t77;
t57 = -qJD(3) * pkin(3) - t197;
t22 = t69 * pkin(4) - t71 * qJ(5) + t57;
t193 = t22 * t78 - t186;
t170 = t90 * t77;
t138 = t78 * t86 + t170;
t149 = t90 * t155;
t192 = t138 * t107 - t122 * t42 + t78 * t149 + t87 * t69;
t61 = -pkin(3) * t122 - t90 * pkin(7) + t100;
t68 = -t108 * t95 + t184 * t96;
t166 = t107 * t61 + t109 * t68;
t45 = t196 * qJD(3) + t115;
t60 = t87 * pkin(3) - t86 * pkin(7);
t10 = -qJD(4) * t166 - t107 * t45 + t109 * t60;
t190 = t71 ^ 2;
t189 = t78 ^ 2;
t188 = t84 ^ 2;
t187 = pkin(7) * t71;
t118 = t90 * qJD(2);
t35 = qJD(1) * t118 + t63 * qJD(3);
t180 = t35 * t196;
t179 = t35 * t90;
t178 = t69 * t78;
t177 = t69 * t204;
t176 = t71 * t69;
t175 = t71 * t84;
t174 = t77 * t122;
t173 = t78 * t84;
t172 = t84 * t69;
t171 = t84 * t204;
t136 = pkin(4) * t107 - qJ(5) * t109;
t168 = t107 * qJD(5) - t78 * t136 + t63;
t167 = -t107 * t42 - t69 * t155;
t59 = t84 * pkin(3) - pkin(7) * t204;
t26 = t107 * t59 + t109 * t197;
t74 = t107 * t77;
t165 = t107 * t204;
t76 = t109 * t77;
t163 = t109 * t78;
t160 = t77 * qJ(5);
t158 = qJD(5) - t20;
t157 = t105 ^ 2 + t106 ^ 2;
t153 = pkin(7) * qJD(4) * t78;
t152 = t69 ^ 2 - t190;
t150 = t90 * t156;
t65 = t71 * t156;
t145 = t107 * t78;
t144 = t157 * qJD(1) ^ 2;
t7 = t42 * pkin(4) + t41 * qJ(5) - t71 * qJD(5) + t35;
t143 = -t7 - t153;
t142 = t35 + t153;
t141 = t71 * t165 - t65;
t140 = -t22 * t86 - t7 * t90;
t139 = t57 * t86 + t179;
t137 = (qJD(4) * t69 - t41) * pkin(7);
t15 = -t78 * pkin(4) + t158;
t135 = t107 * t16 - t109 * t15;
t134 = t107 * t21 + t109 * t20;
t25 = -t107 * t197 + t109 * t59;
t28 = -t107 * t68 + t109 * t61;
t128 = t78 * t155 - t163 * t204 + t74;
t127 = t76 + (-t156 + t165) * t78;
t126 = 0.2e1 * t157 * qJD(2) * qJD(1);
t125 = t22 * t71 + t146;
t123 = t78 * t57 - t186;
t9 = t107 * t60 + t109 * t45 + t61 * t155 - t68 * t156;
t121 = (t41 + t177) * t109 + t167;
t120 = t69 * t145 - t161;
t119 = -t127 - t172;
t112 = t69 * t149 + (t42 * t90 + t69 * t86) * t107;
t46 = t68 * qJD(3) + t118;
t111 = t71 * t78 - t42;
t94 = -t109 * pkin(4) - t107 * qJ(5) - pkin(3);
t80 = t204 ^ 2;
t37 = t71 * pkin(4) + t69 * qJ(5);
t36 = pkin(7) * t161;
t31 = t136 * t90 - t196;
t30 = t78 * t87 - t174;
t24 = pkin(4) * t122 - t28;
t23 = -qJ(5) * t122 + t166;
t19 = -t41 + t178;
t18 = -t84 * pkin(4) - t25;
t17 = t84 * qJ(5) + t26;
t14 = t128 - t175;
t13 = t71 * t163 - t162;
t12 = (qJ(5) * qJD(4) * t90 + pkin(4) * t86) * t107 + (-qJ(5) * t86 + (pkin(4) * qJD(4) - qJD(5)) * t90) * t109 + t46;
t11 = -t90 * t65 + (-t41 * t90 + t71 * t86) * t109;
t8 = -t87 * pkin(4) - t10;
t4 = t87 * qJ(5) - qJD(5) * t122 + t9;
t3 = t138 * t109 + t122 * t41 - t78 * t150 + t71 * t87;
t1 = t78 * qJD(5) + t160 + t5;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, qJ(2) * t126, t90 * t113 + t84 * t86, t113 * t122 + t204 * t86 - t84 * t87 - t170, t86 * qJD(3), -t204 * t87 - t174, -t87 * qJD(3), 0, -t46 * qJD(3) + t100 * t77 + t93 * t87, -t45 * qJD(3) + t100 * t113 + t93 * t86, -t113 * t196 + t122 * t34 - t197 * t86 + t204 * t45 + t46 * t84 - t63 * t87 - t68 * t77 + t179, -t197 * t46 + t34 * t68 + t63 * t45 - t180, t11, t202, t3, t112, -t192, t30, t10 * t78 + t139 * t107 + t122 * t146 + t57 * t149 - t196 * t42 + t20 * t87 + t28 * t77 + t46 * t69, t139 * t109 + t122 * t5 - t57 * t150 - t166 * t77 + t196 * t41 - t21 * t87 + t46 * t71 - t9 * t78, -t10 * t71 + t28 * t41 - t166 * t42 - t9 * t69 - t134 * t86 + (-t5 * t107 + t146 * t109 + (t107 * t20 - t109 * t21) * qJD(4)) * t90, t20 * t10 - t146 * t28 + t166 * t5 + t21 * t9 + t57 * t46 - t180, t11, t3, -t202, t30, t192, t112, -t140 * t107 + t12 * t69 + t122 * t2 + t149 * t22 - t15 * t87 - t24 * t77 + t31 * t42 - t8 * t78, -t23 * t42 - t24 * t41 - t4 * t69 + t8 * t71 - t135 * t86 + (-t1 * t107 + t2 * t109 + (-t107 * t15 - t109 * t16) * qJD(4)) * t90, -t1 * t122 + t109 * t140 - t12 * t71 + t150 * t22 + t16 * t87 + t23 * t77 + t31 * t41 + t4 * t78, t1 * t23 + t22 * t12 + t15 * t8 + t16 * t4 + t2 * t24 + t7 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, -qJ(2) * t144, 0, 0, 0, 0, 0, 0, 0.2e1 * t195, 0.2e1 * t113, -t80 - t188, t197 * t84 - t204 * t63, 0, 0, 0, 0, 0, 0, t127 - t172, -t109 * t189 - t175 - t74, t71 * t145 + t121, t200 * t107 + t199 * t109 - t57 * t84, 0, 0, 0, 0, 0, 0, -t107 * t189 - t172 + t76, t121 - t141, t128 + t175, -t22 * t84 + t201 * t109 + (t15 * t78 + t1) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, -t80 + t188, 0, t171, 0, 0, -t205 * t84, -t205 * t204, 0, 0, t13, (-t41 + t177) * t109 + t141 + t167, t14, t120, -t119, -t173, -pkin(3) * t42 + t123 * t107 - t142 * t109 - t20 * t84 - t25 * t78 - t63 * t69, pkin(3) * t41 + t142 * t107 + t123 * t109 + t21 * t84 + t26 * t78 - t63 * t71, t25 * t71 + t26 * t69 - t36 + (t20 * t204 + t5 + (-t20 + t187) * qJD(4)) * t109 + (t137 - t199) * t107, -t35 * pkin(3) - t20 * t25 - t21 * t26 - t57 * t63 + (-t134 * qJD(4) + t107 * t146 + t5 * t109) * pkin(7), t13, t14, t65 + (-t204 * t71 + t42) * t107 + (t41 + t178) * t109, -t173, t119, t120, t193 * t107 + t143 * t109 + t15 * t84 - t168 * t69 + t18 * t78 + t94 * t42, t17 * t69 - t18 * t71 - t36 + (-t15 * t204 + t1 + (t15 + t187) * qJD(4)) * t109 + (t137 - t201) * t107, t143 * t107 - t193 * t109 - t16 * t84 + t168 * t71 - t17 * t78 + t94 * t41, -t15 * t18 - t16 * t17 + t7 * t94 - t168 * t22 + (-qJD(4) * t135 + t1 * t109 + t2 * t107) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, -t152, t19, -t176, t111, t77, -t57 * t71 + t199, t57 * t69 - t200, 0, 0, t176, t19, t152, t77, -t111, -t176, -t37 * t69 - t125 + t181 + 0.2e1 * t185, pkin(4) * t41 - t42 * qJ(5) + (t16 - t21) * t71 + (t15 - t158) * t69, 0.2e1 * t160 - t22 * t69 + t37 * t71 + (0.2e1 * qJD(5) - t20) * t78 + t5, -t2 * pkin(4) + t1 * qJ(5) - t15 * t21 + t158 * t16 - t22 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176 - t195, t19, -t189 - t190, t125 - t183 - t185;];
tauc_reg = t6;
