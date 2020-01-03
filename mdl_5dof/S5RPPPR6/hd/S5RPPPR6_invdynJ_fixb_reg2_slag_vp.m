% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:55
% EndTime: 2019-12-31 17:47:59
% DurationCPUTime: 1.56s
% Computational Cost: add. (1779->265), mult. (3984->385), div. (0->0), fcn. (2752->8), ass. (0->158)
t197 = pkin(3) + qJ(2);
t112 = sin(qJ(5));
t109 = sin(pkin(7));
t166 = qJD(1) * t109;
t156 = t112 * t166;
t108 = sin(pkin(8));
t111 = cos(pkin(7));
t114 = cos(qJ(5));
t170 = t111 * t114;
t157 = t108 * t170;
t45 = qJD(1) * t157 - t156;
t172 = t111 * t112;
t179 = t109 * t114;
t48 = t108 * t172 + t179;
t42 = t48 * qJD(5);
t93 = t109 * qJDD(1);
t19 = qJD(1) * t42 - qJDD(1) * t157 + t112 * t93;
t41 = t48 * qJD(1);
t110 = cos(pkin(8));
t165 = qJD(1) * t111;
t64 = t110 * t165 + qJD(5);
t196 = -t41 * t64 + t19;
t20 = t45 * qJD(5) + t48 * qJDD(1);
t195 = -t45 * t64 + t20;
t107 = t111 ^ 2;
t161 = qJD(1) * qJD(2);
t154 = t107 * t161;
t159 = t111 * qJDD(1);
t38 = t111 * t161 + t197 * t159 + qJDD(4);
t66 = t197 * t111;
t194 = (qJDD(1) * t66 + t38) * t111 + t154;
t149 = -t109 * qJ(3) - pkin(1);
t54 = (-pkin(2) - qJ(4)) * t111 + t149;
t34 = t54 * qJD(1) + qJD(2);
t71 = qJ(2) * t166 + qJD(3);
t56 = pkin(3) * t166 + t71;
t16 = t108 * t56 + t110 * t34;
t14 = pkin(6) * t166 + t16;
t130 = (pkin(4) * t110 + pkin(6) * t108) * t111;
t57 = t197 * t165 + qJD(4);
t27 = qJD(1) * t130 + t57;
t138 = t112 * t14 - t114 * t27;
t21 = qJDD(1) * t130 + t38;
t61 = -qJD(3) * t109 - qJD(4) * t111;
t24 = t61 * qJD(1) + t54 * qJDD(1) + qJDD(2);
t55 = qJ(2) * t93 + t109 * t161 + qJDD(3);
t37 = pkin(3) * t93 + t55;
t12 = t108 * t37 + t110 * t24;
t8 = pkin(6) * t93 + t12;
t1 = -t138 * qJD(5) + t112 * t21 + t114 * t8;
t113 = sin(qJ(1));
t193 = g(1) * t113;
t115 = cos(qJ(1));
t103 = g(2) * t115;
t192 = g(3) * t109;
t191 = t41 * t45;
t65 = t197 * t109;
t26 = t108 * t65 + t110 * t54;
t117 = qJ(2) ^ 2;
t78 = 0.2e1 * t154;
t94 = t107 * qJDD(1);
t188 = qJ(2) * t78 + t117 * t94;
t187 = t114 * t64;
t98 = t115 * qJ(2);
t186 = t115 * pkin(3) + t98;
t185 = t115 * pkin(1) + t113 * qJ(2);
t184 = qJD(5) * t64;
t183 = qJDD(1) * pkin(1);
t182 = t108 * t113;
t116 = qJD(1) ^ 2;
t181 = t108 * t116;
t180 = t109 * t112;
t178 = t109 * t115;
t177 = t109 * t116;
t176 = t110 * t111;
t175 = t110 * t113;
t174 = t110 * t115;
t173 = t110 * t116;
t171 = t111 * t113;
t169 = t111 * t115;
t104 = t108 ^ 2;
t106 = t110 ^ 2;
t168 = -t104 - t106;
t105 = t109 ^ 2;
t167 = t105 + t107;
t164 = qJD(2) * t109;
t163 = qJD(2) * t111;
t59 = -pkin(2) * t111 + t149;
t162 = qJDD(1) * t59;
t160 = qJD(1) * qJD(3);
t158 = 0.2e1 * t94;
t155 = -t113 * pkin(1) + t98;
t153 = t108 * t93;
t152 = t109 * t160;
t151 = t109 * t159;
t150 = t110 * t93;
t148 = -t103 + t193;
t147 = t168 * t116;
t60 = t167 * t116;
t146 = pkin(2) * t169 + qJ(3) * t178 + t185;
t144 = -0.2e1 * t151;
t91 = qJDD(2) - t183;
t50 = -t109 * t174 + t182;
t52 = t108 * t115 + t109 * t175;
t142 = g(1) * t52 + g(2) * t50;
t141 = g(1) * t115 + g(2) * t113;
t6 = t112 * t27 + t114 * t14;
t140 = t112 * t138 + t114 * t6;
t11 = -t108 * t24 + t110 * t37;
t15 = -t108 * t34 + t110 * t56;
t25 = -t108 * t54 + t110 * t65;
t23 = pkin(6) * t109 + t26;
t31 = t130 + t66;
t10 = t112 * t31 + t114 * t23;
t9 = -t112 * t23 + t114 * t31;
t137 = t113 * pkin(3) + qJ(4) * t169 + t146;
t136 = t183 - t91 - t103;
t124 = qJDD(2) + t162;
t33 = t124 - t152;
t135 = -t162 - t33 - t103;
t63 = t110 * t159 + qJDD(5);
t134 = -t112 * t63 - t114 * t184;
t133 = -t112 * t184 + t114 * t63;
t53 = -t109 * t182 + t174;
t132 = t112 * t53 + t113 * t170;
t131 = t112 * t171 - t114 * t53;
t129 = t140 * t110;
t49 = t157 - t180;
t35 = t108 * t61 - t110 * t164;
t128 = -qJD(1) * t35 + qJDD(1) * t25 + t11;
t36 = t108 * t164 + t110 * t61;
t127 = -qJD(1) * t36 - qJDD(1) * t26 - t12;
t126 = qJ(2) * t158 - t141 + t78;
t125 = (qJ(2) * qJDD(1) + t161) * t105;
t123 = t54 * t193;
t122 = -t41 * t166 + t134;
t121 = -t45 * t166 - t133;
t2 = -t6 * qJD(5) - t112 * t8 + t114 * t21;
t120 = (-t108 * t15 + t110 * t16) * qJD(1) - t141;
t119 = t112 * t19 + t114 * t20 + (-t112 * t41 - t114 * t45) * qJD(5);
t13 = -pkin(4) * t166 - t15;
t118 = t13 * t166 + t1 * t114 - t112 * t2 + (-t112 * t6 + t114 * t138) * qJD(5);
t99 = g(3) * t111;
t92 = t105 * qJDD(1);
t80 = g(1) * t171;
t70 = 0.2e1 * t151;
t51 = t108 * t178 + t175;
t47 = t59 * qJD(1) + qJD(2);
t44 = (t108 * t179 + t172) * qJD(1);
t43 = t49 * qJD(5);
t40 = (-t108 * t180 + t170) * qJD(1);
t29 = t112 * t169 + t114 * t51;
t28 = -t112 * t51 + t114 * t169;
t22 = -pkin(4) * t109 - t25;
t7 = -pkin(4) * t93 - t11;
t4 = -t10 * qJD(5) - t112 * t36 + t114 * t163;
t3 = t9 * qJD(5) + t112 * t163 + t114 * t36;
t5 = [0, 0, 0, 0, 0, qJDD(1), t148, t141, 0, 0, t92, t70, 0, t94, 0, 0, t136 * t111 + t80, (-t136 - t193) * t109, 0.2e1 * t125 + t126, -t91 * pkin(1) - g(1) * t155 - g(2) * t185 + (0.2e1 * qJ(2) * t161 + qJDD(1) * t117) * t105 + t188, 0, 0, 0, t92, t70, t94, t55 * t109 + t125 + t126, -t80 + (-t135 - t152) * t111, t105 * t160 + (t135 + t193) * t109, t33 * t59 - g(1) * (-pkin(2) * t171 + t155) - g(2) * t146 + (qJ(2) * t55 + qJ(3) * t193 + qJD(2) * t71 - qJD(3) * t47) * t109 + t188, t104 * t94, t108 * t110 * t158, t108 * t144, t106 * t94, t110 * t144, t92, -g(1) * t53 - g(2) * t51 + t128 * t109 + t194 * t110, -t194 * t108 + t127 * t109 + t142, t80 + (t108 * t128 + t110 * t127 - t103) * t111, -g(1) * t186 - g(2) * t137 + t11 * t25 + t12 * t26 - t15 * t35 + t16 * t36 + t57 * t163 + t38 * t66 - t123, -t19 * t49 - t42 * t45, t19 * t48 - t20 * t49 + t41 * t42 - t43 * t45, t176 * t19 + t42 * t64 - t49 * t63, t20 * t48 + t41 * t43, t176 * t20 + t43 * t64 + t48 * t63, t63 * t176, g(1) * t131 - g(2) * t29 - t13 * t43 + t176 * t2 - t22 * t20 - t35 * t41 + t4 * t64 - t7 * t48 + t9 * t63, g(1) * t132 - g(2) * t28 - t1 * t176 - t10 * t63 + t13 * t42 + t22 * t19 - t3 * t64 - t35 * t45 - t7 * t49, t1 * t48 + t10 * t20 + t138 * t42 - t19 * t9 + t2 * t49 + t3 * t41 + t4 * t45 + t43 * t6 - t142, t1 * t10 + t6 * t3 + t2 * t9 - t138 * t4 + t7 * t22 + t13 * t35 - g(1) * (t53 * pkin(4) + t52 * pkin(6) + t186) - g(2) * (pkin(4) * t51 + pkin(6) * t50 + t137) - t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159, t93, -t60, -qJ(2) * t60 - t148 + t91, 0, 0, 0, 0, 0, 0, -t60, t159, -t93, -t107 * t116 * qJ(2) + (-qJD(3) - t71) * t166 + t124 - t148, 0, 0, 0, 0, 0, 0, -t110 * t60 - t153, t167 * t181 - t150, t168 * t159, -t11 * t108 + t12 * t110 + (-t111 * t57 + (-t108 * t16 - t110 * t15) * t109) * qJD(1) - t148, 0, 0, 0, 0, 0, 0, -t108 * t20 + t110 * t122 - t40 * t64, t108 * t19 + t110 * t121 + t44 * t64, t110 * t119 - t40 * t45 - t44 * t41, t7 * t108 + t110 * t118 + t138 * t40 - t6 * t44 - t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t111 * t177, -t105 * t116, t99 + (qJD(1) * t47 - t141) * t109 + t55, 0, 0, 0, 0, 0, 0, -t105 * t181 + t150, -t105 * t173 - t153, t111 * t109 * t147, t12 * t108 + t109 * t120 + t11 * t110 + t99, 0, 0, 0, 0, 0, 0, (-t156 * t64 + t20) * t110 + t122 * t108, (-t166 * t187 - t19) * t110 + t121 * t108, (-t112 * t45 + t114 * t41) * t110 * t166 + t119 * t108, -t7 * t110 + t99 + (qJD(1) * t129 - t141) * t109 + t118 * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (qJDD(1) * t110 - t108 * t177) * t111, (-qJDD(1) * t108 - t109 * t173) * t111, t107 * t147, t111 * t120 - t192 + t38, 0, 0, 0, 0, 0, 0, (-t110 * t112 * t64 - t108 * t41) * t165 + t133, (-t108 * t45 - t110 * t187) * t165 + t134, t195 * t112 - t196 * t114, -t192 + t1 * t112 + t2 * t114 + t140 * qJD(5) + ((t108 * t13 + t129) * qJD(1) - t141) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, -t41 ^ 2 + t45 ^ 2, t196, -t191, t195, t63, -g(1) * t28 - g(2) * t132 - g(3) * t48 + t13 * t45 + t6 * t64 + t2, g(1) * t29 + g(2) * t131 - g(3) * t49 - t13 * t41 - t138 * t64 - t1, 0, 0;];
tau_reg = t5;
