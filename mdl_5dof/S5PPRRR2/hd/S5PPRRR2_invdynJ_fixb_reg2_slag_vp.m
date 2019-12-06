% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:56
% EndTime: 2019-12-05 15:14:59
% DurationCPUTime: 1.52s
% Computational Cost: add. (1869->237), mult. (4254->316), div. (0->0), fcn. (3342->14), ass. (0->149)
t114 = qJD(4) ^ 2;
t104 = sin(pkin(8));
t106 = cos(pkin(8));
t140 = g(1) * t106 + g(2) * t104;
t98 = pkin(9) + qJ(3);
t89 = sin(t98);
t90 = cos(t98);
t126 = -g(3) * t90 + t140 * t89;
t103 = sin(pkin(9));
t105 = cos(pkin(9));
t109 = sin(qJ(3));
t112 = cos(qJ(3));
t57 = t103 * t112 + t105 * t109;
t50 = t57 * qJD(1);
t121 = qJD(3) * t50 + t126;
t178 = qJDD(3) * pkin(3);
t165 = qJD(1) * qJD(3);
t150 = t109 * t165;
t161 = qJDD(1) * t112;
t202 = qJDD(1) * t109 + t112 * t165;
t145 = (t150 - t161) * t105 + t202 * t103;
t32 = t145 - t178;
t205 = -pkin(6) * t114 + t121 + t178 - t32;
t107 = sin(qJ(5));
t110 = cos(qJ(5));
t111 = cos(qJ(4));
t159 = t111 * qJDD(3);
t108 = sin(qJ(4));
t160 = t108 * qJDD(3);
t59 = t107 * t111 + t108 * t110;
t99 = qJD(4) + qJD(5);
t204 = t99 * t59;
t21 = qJD(3) * t204 + t107 * t160 - t110 * t159;
t157 = -t103 * t161 - t202 * t105;
t130 = -t103 * t150 - t157;
t31 = qJDD(3) * pkin(6) + t130;
t203 = qJD(4) * qJD(2) + t31;
t175 = t110 * t111;
t176 = t107 * t108;
t58 = -t175 + t176;
t29 = t58 * t57;
t152 = -g(1) * t104 + g(2) * t106;
t199 = g(3) * t89;
t127 = t140 * t90 + t199;
t168 = qJD(5) * t107;
t48 = qJD(3) * pkin(6) + t50;
t147 = pkin(7) * qJD(3) + t48;
t93 = t111 * qJD(2);
t34 = -t147 * t108 + t93;
t33 = qJD(4) * pkin(4) + t34;
t166 = t108 * qJD(2);
t35 = t147 * t111 + t166;
t92 = t111 * qJDD(2);
t6 = qJDD(4) * pkin(4) + t92 + (-pkin(7) * qJDD(3) - t31) * t108 - t35 * qJD(4);
t158 = t108 * qJDD(2) + t203 * t111;
t169 = qJD(4) * t108;
t10 = -t48 * t169 + t158;
t164 = qJD(3) * qJD(4);
t131 = t108 * t164 - t159;
t7 = -t131 * pkin(7) + t10;
t1 = (qJD(5) * t33 + t7) * t110 + t107 * t6 - t35 * t168;
t200 = pkin(7) + pkin(6);
t153 = qJD(3) * t175;
t170 = qJD(3) * t108;
t154 = t107 * t170;
t53 = -t153 + t154;
t55 = t59 * qJD(3);
t193 = t55 * t53;
t67 = t200 * t108;
t68 = t200 * t111;
t40 = -t107 * t68 - t110 * t67;
t174 = t112 * t105;
t177 = t103 * t109;
t78 = qJD(1) * t177;
t49 = qJD(1) * t174 - t78;
t155 = qJD(4) * t200;
t60 = t108 * t155;
t61 = t111 * t155;
t192 = t40 * qJD(5) - t107 * t61 - t110 * t60 + t58 * t49;
t41 = -t107 * t67 + t110 * t68;
t191 = -t41 * qJD(5) + t107 * t60 - t110 * t61 + t59 * t49;
t134 = t99 * t176;
t167 = qJD(5) * t110;
t36 = -qJD(4) * t175 - t111 * t167 + t134;
t190 = -t59 * t21 + t36 * t53;
t189 = qJD(3) * pkin(3);
t102 = qJ(4) + qJ(5);
t94 = sin(t102);
t188 = t104 * t94;
t95 = cos(t102);
t187 = t104 * t95;
t186 = t106 * t94;
t185 = t106 * t95;
t184 = t107 * t35;
t183 = t108 * t48;
t182 = t110 * t35;
t56 = -t174 + t177;
t179 = qJD(3) * t56;
t180 = qJD(3) * t179;
t39 = t111 * t48 + t166;
t173 = t39 * qJD(4);
t100 = t108 ^ 2;
t101 = t111 ^ 2;
t172 = t100 - t101;
t171 = t100 + t101;
t115 = qJD(3) ^ 2;
t156 = t108 * t115 * t111;
t88 = pkin(4) * t111 + pkin(3);
t148 = t111 * t164;
t144 = -qJD(5) * t153 - t107 * t159 + (-t148 - t160) * t110;
t143 = t171 * qJDD(3);
t142 = t108 * t148;
t141 = pkin(4) * t169 - t50;
t20 = qJD(3) * t134 + t144;
t138 = -t20 * t58 + t204 * t55;
t97 = qJDD(4) + qJDD(5);
t137 = t36 * t99 - t59 * t97;
t9 = t107 * t33 + t182;
t38 = t93 - t183;
t133 = t108 * t38 - t111 * t39;
t52 = t57 * qJD(3);
t132 = -qJD(3) * t52 - qJDD(3) * t56;
t129 = t114 * t57 - t132;
t128 = 0.2e1 * t179 * qJD(4) - qJDD(4) * t57;
t47 = -t49 - t189;
t123 = -pkin(6) * qJDD(4) + (t47 + t49 - t189) * qJD(4);
t2 = -t9 * qJD(5) - t107 * t7 + t110 * t6;
t122 = -t47 * qJD(3) + t127;
t11 = -t108 * t31 - t173 + t92;
t120 = t10 * t111 - t11 * t108 + (-t108 * t39 - t111 * t38) * qJD(4);
t42 = -t88 * qJD(3) - t49;
t118 = -g(1) * (-t90 * t185 - t188) - g(2) * (-t90 * t187 + t186) + t42 * t53 + t95 * t199 - t1;
t117 = -g(1) * (-t90 * t186 + t187) - g(2) * (-t90 * t188 - t185) - t42 * t55 + t2 + t94 * t199;
t116 = t120 - t127;
t66 = qJDD(4) * t111 - t108 * t114;
t65 = qJDD(4) * t108 + t111 * t114;
t62 = qJDD(2) + t152;
t28 = t59 * t57;
t23 = -t53 ^ 2 + t55 ^ 2;
t22 = t131 * pkin(4) + t32;
t19 = -t204 * t99 - t58 * t97;
t15 = t55 * t99 - t21;
t14 = -t144 + (-t154 + t53) * t99;
t13 = t110 * t34 - t184;
t12 = -t107 * t34 - t182;
t8 = t110 * t33 - t184;
t4 = t179 * t59 + t99 * t29;
t3 = t179 * t58 - t204 * t57;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) + (t103 ^ 2 + t105 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t132, -qJDD(3) * t57 + t180, 0, t130 * t57 + t145 * t56 - t179 * t50 - t49 * t52 - g(3), 0, 0, 0, 0, 0, 0, t128 * t108 - t129 * t111, t129 * t108 + t128 * t111, t57 * t143 - t171 * t180, t120 * t57 + t133 * t179 + t32 * t56 + t47 * t52 - g(3), 0, 0, 0, 0, 0, 0, t21 * t56 - t28 * t97 + t4 * t99 + t52 * t53, -t20 * t56 + t29 * t97 - t3 * t99 + t52 * t55, -t20 * t28 + t21 * t29 - t3 * t53 - t4 * t55, -t1 * t29 - t2 * t28 + t22 * t56 + t3 * t9 + t4 * t8 + t42 * t52 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, t66, -t65, 0, -t133 * qJD(4) + t10 * t108 + t11 * t111 + t152, 0, 0, 0, 0, 0, 0, t19, t137, t138 + t190, t1 * t59 - t2 * t58 - t204 * t8 - t36 * t9 + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t121 - t145, (t49 + t78) * qJD(3) + t127 + t157, 0, 0, qJDD(3) * t100 + 0.2e1 * t142, 0.2e1 * t108 * t159 - 0.2e1 * t172 * t164, t65, qJDD(3) * t101 - 0.2e1 * t142, t66, 0, t123 * t108 + t111 * t205, -t205 * t108 + t123 * t111, -t171 * t49 * qJD(3) + pkin(6) * t143 + t116, -t47 * t50 + t133 * t49 + (-t32 + t126) * pkin(3) + t116 * pkin(6), -t20 * t59 - t36 * t55, -t138 + t190, -t137, t204 * t53 + t21 * t58, t19, 0, t126 * t95 + t141 * t53 + t191 * t99 + t204 * t42 - t21 * t88 + t22 * t58 + t40 * t97, -t126 * t94 + t141 * t55 - t192 * t99 + t20 * t88 + t22 * t59 - t36 * t42 - t41 * t97, -t1 * t58 - t191 * t55 - t192 * t53 - t2 * t59 + t20 * t40 - t204 * t9 - t21 * t41 + t36 * t8 - t127, t1 * t41 + t2 * t40 - t22 * t88 - g(3) * (t200 * t89 + t88 * t90) + t192 * t9 + t191 * t8 + t141 * t42 + t140 * (-t200 * t90 + t88 * t89); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, t172 * t115, t160, t156, t159, qJDD(4), t173 + t92 + (-qJD(4) * t48 + t152) * t111 + (t122 - t203) * t108, -t152 * t108 + t122 * t111 + (t38 + t183) * qJD(4) - t158, 0, 0, t193, t23, t14, -t193, t15, t97, -t12 * t99 + (t110 * t97 - t99 * t168 - t53 * t170) * pkin(4) + t117, t13 * t99 + (-t107 * t97 - t99 * t167 - t55 * t170) * pkin(4) + t118, (t12 + t9) * t55 + (t13 - t8) * t53 + (-t107 * t21 + t110 * t20 + (t107 * t55 - t110 * t53) * qJD(5)) * pkin(4), -t8 * t12 - t9 * t13 + (t1 * t107 + t2 * t110 + t152 * t111 + (-t42 * qJD(3) + t127) * t108 + (-t107 * t8 + t110 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t23, t14, -t193, t15, t97, t9 * t99 + t117, t8 * t99 + t118, 0, 0;];
tau_reg = t5;
