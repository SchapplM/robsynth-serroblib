% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:41
% EndTime: 2020-01-03 11:20:46
% DurationCPUTime: 1.62s
% Computational Cost: add. (1859->225), mult. (3932->326), div. (0->0), fcn. (2778->14), ass. (0->144)
t106 = sin(pkin(8));
t112 = sin(qJ(5));
t108 = cos(pkin(9));
t114 = cos(qJ(5));
t162 = t114 * t108;
t146 = t106 * t162;
t105 = sin(pkin(9));
t152 = qJDD(1) * t105;
t71 = t114 * t105 + t112 * t108;
t189 = t71 * qJD(5);
t19 = -qJDD(1) * t146 + (qJD(1) * t189 + t112 * t152) * t106;
t104 = qJ(1) + pkin(7);
t93 = sin(t104);
t95 = cos(t104);
t190 = g(2) * t95 + g(3) * t93;
t110 = cos(pkin(7));
t90 = -t110 * pkin(1) - pkin(2);
t155 = qJDD(1) * t90;
t74 = qJDD(3) + t155;
t192 = t190 + t74;
t109 = cos(pkin(8));
t82 = t109 * qJD(1) - qJD(5);
t160 = qJD(1) * t106;
t145 = t105 * t160;
t138 = t112 * t145;
t144 = qJD(5) * t162;
t20 = -qJD(5) * t138 + (qJD(1) * t144 + t71 * qJDD(1)) * t106;
t166 = pkin(1) * qJDD(1);
t107 = sin(pkin(7));
t84 = t107 * pkin(1) + qJ(3);
t68 = qJD(1) * qJD(3) + qJDD(1) * t84;
t183 = pkin(6) * t106;
t147 = t108 * t183;
t126 = -t109 * pkin(4) - t147;
t157 = qJD(4) * t106;
t67 = -t109 * pkin(3) - t106 * qJ(4) + t90;
t29 = -qJD(1) * t157 + t67 * qJDD(1) + qJDD(3);
t42 = t106 * qJDD(2) + t109 * t68;
t14 = -t105 * t42 + t108 * t29;
t10 = t126 * qJDD(1) + t14;
t150 = t106 * qJDD(1);
t143 = t105 * t150;
t15 = t105 * t29 + t108 * t42;
t11 = -pkin(6) * t143 + t15;
t43 = t67 * qJD(1) + qJD(3);
t78 = t84 * qJD(1);
t56 = t106 * qJD(2) + t109 * t78;
t17 = -t105 * t56 + t108 * t43;
t13 = t126 * qJD(1) + t17;
t18 = t105 * t43 + t108 * t56;
t16 = -pkin(6) * t145 + t18;
t127 = t112 * t16 - t114 * t13;
t1 = -t127 * qJD(5) + t112 * t10 + t114 * t11;
t49 = qJD(1) * t146 - t138;
t187 = t49 ^ 2;
t186 = g(2) * t93;
t185 = g(3) * t95;
t184 = pkin(4) * t105;
t182 = g(1) * t106;
t122 = qJD(1) * t71;
t46 = t106 * t122;
t181 = t46 * t82;
t180 = t49 * t46;
t179 = t49 * t82;
t52 = t106 * t189;
t163 = t112 * t105;
t70 = t162 - t163;
t58 = t70 * t106;
t178 = -t58 * t20 + t52 * t46;
t53 = (-qJD(5) * t163 + t144) * t106;
t57 = t71 * t106;
t149 = t109 * qJDD(1);
t81 = -qJDD(5) + t149;
t177 = t53 * t82 + t57 * t81;
t176 = -t109 * t122 + t189;
t175 = t82 * t70;
t171 = t109 * t84;
t28 = t105 * t67 + t108 * t171;
t113 = sin(qJ(1));
t174 = t113 * pkin(1) + t93 * pkin(2);
t172 = t109 * t20;
t170 = t109 * t93;
t169 = t109 * t95;
t168 = t19 * t109;
t101 = t108 ^ 2;
t99 = t105 ^ 2;
t167 = -t101 - t99;
t165 = t106 * (-pkin(6) - qJ(4));
t116 = qJD(1) ^ 2;
t164 = t109 * t116;
t100 = t106 ^ 2;
t102 = t109 ^ 2;
t161 = t100 + t102;
t159 = qJD(3) * t106;
t158 = qJD(3) * t109;
t151 = t100 * qJDD(1);
t115 = cos(qJ(1));
t148 = t115 * pkin(1) + t95 * pkin(2) + t93 * qJ(3);
t141 = t108 * t149;
t140 = t161 * t116;
t55 = t109 * qJD(2) - t106 * t78;
t41 = t109 * qJDD(2) - t106 * t68;
t139 = 0.2e1 * t106 * t149;
t137 = -t95 * qJ(3) + t174;
t134 = t185 - t186;
t133 = -g(2) * t115 - g(3) * t113;
t132 = -t57 * t19 + t49 * t53;
t131 = -t52 * t82 + t58 * t81;
t54 = qJD(4) - t55;
t129 = -t41 * t106 + t42 * t109;
t128 = t55 * t106 - t56 * t109;
t6 = t112 * t13 + t114 * t16;
t60 = t108 * t67;
t21 = -t147 + t60 + (-t105 * t84 - pkin(4)) * t109;
t23 = -t105 * t183 + t28;
t7 = -t112 * t23 + t114 * t21;
t8 = t112 * t21 + t114 * t23;
t36 = qJDD(4) - t41;
t27 = -t105 * t171 + t60;
t63 = -t105 * t158 - t108 * t157;
t125 = -t63 * qJD(1) - t27 * qJDD(1) - t14;
t64 = -t105 * t157 + t108 * t158;
t124 = t64 * qJD(1) + t28 * qJDD(1) + t15;
t123 = g(1) * t109 + t106 * t185 + t36;
t120 = t155 + t192;
t2 = -t6 * qJD(5) + t114 * t10 - t112 * t11;
t119 = t68 * t100 + t36 * t106 + t134;
t103 = pkin(9) + qJ(5);
t94 = cos(t103);
t92 = sin(t103);
t91 = t102 * qJDD(1);
t89 = t108 * pkin(4) + pkin(3);
t80 = pkin(4) * t143;
t61 = (t84 + t184) * t106;
t45 = t46 ^ 2;
t40 = t94 * t169 + t93 * t92;
t39 = t92 * t169 - t93 * t94;
t38 = t94 * t170 - t95 * t92;
t37 = -t92 * t170 - t95 * t94;
t31 = pkin(4) * t145 + t54;
t26 = t36 + t80;
t4 = -t8 * qJD(5) - t112 * t64 + t114 * t63;
t3 = t7 * qJD(5) + t112 * t63 + t114 * t64;
t5 = [0, 0, 0, 0, 0, qJDD(1), t133, g(2) * t113 - g(3) * t115, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t110 * t166 - t190, -0.2e1 * t107 * t166 - t134, 0, (t133 + (t107 ^ 2 + t110 ^ 2) * t166) * pkin(1), t151, t139, 0, t91, 0, 0, -t120 * t109, t120 * t106, t68 * t161 + t129 + t134, -g(2) * t148 - g(3) * t137 - t128 * qJD(3) + t129 * t84 + t74 * t90, t101 * t151, -0.2e1 * t108 * t105 * t151, -0.2e1 * t106 * t141, t99 * t151, t105 * t139, t91, (-t108 * t190 + t125) * t109 + t119 * t105, (t105 * t190 + t124) * t109 + t119 * t108, (-t124 * t105 + t125 * t108 - t190) * t106, t15 * t28 + t18 * t64 + t14 * t27 + t17 * t63 - g(2) * (pkin(3) * t169 + t148) - g(3) * (pkin(3) * t170 + t137) + (-qJ(4) * t190 + t54 * qJD(3) + t36 * t84) * t106, -t19 * t58 - t49 * t52, -t132 + t178, -t131 + t168, t20 * t57 + t46 * t53, t172 + t177, t81 * t109, -g(2) * t40 - g(3) * t38 - t2 * t109 + t46 * t159 + t61 * t20 + t26 * t57 + t31 * t53 - t4 * t82 - t7 * t81, g(2) * t39 - g(3) * t37 + t1 * t109 + t49 * t159 - t61 * t19 + t26 * t58 + t3 * t82 - t31 * t52 + t8 * t81, -t1 * t57 - t106 * t190 - t127 * t52 + t7 * t19 - t2 * t58 - t8 * t20 - t3 * t46 - t4 * t49 - t6 * t53, t1 * t8 + t6 * t3 + t2 * t7 - t127 * t4 + t26 * t61 + t31 * t159 - g(2) * (t93 * t184 + t148) - g(3) * (-t93 * t165 + t89 * t170 + t174) + (-g(2) * (t109 * t89 - t165) - g(3) * (-qJ(3) - t184)) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t106 + t41 * t109 - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36 * t109 - g(1) + (-t105 * t14 + t108 * t15) * t106, 0, 0, 0, 0, 0, 0, -t172 + t177, t131 + t168, t132 + t178, t1 * t58 - t26 * t109 + t127 * t53 - t2 * t57 - t6 * t52 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t150, -t140, t128 * qJD(1) + t192, 0, 0, 0, 0, 0, 0, -t105 * t140 - t141, t105 * t149 - t108 * t140, t167 * t150, t15 * t105 + t14 * t108 + (-t106 * t54 + (t105 * t17 - t108 * t18) * t109) * qJD(1) + t190, 0, 0, 0, 0, 0, 0, -t46 * t160 + t176 * t82 - t70 * t81, -t49 * t160 - t175 * t82 + t71 * t81, t175 * t46 + t176 * t49 + t70 * t19 - t71 * t20, t1 * t71 + t127 * t176 - t31 * t160 - t175 * t6 + t2 * t70 + t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t108 * t164 + t152) * t106, (qJDD(1) * t108 + t105 * t164) * t106, t167 * t116 * t100, (-t186 + (t105 * t18 + t108 * t17) * qJD(1)) * t106 + t123, 0, 0, 0, 0, 0, 0, t20 - t179, -t19 + t181, -t45 - t187, -t106 * t186 - t127 * t49 + t6 * t46 + t123 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, -t45 + t187, -t19 - t181, -t180, -t179 - t20, -t81, -g(2) * t37 - g(3) * t39 + t92 * t182 - t31 * t49 - t6 * t82 + t2, g(2) * t38 - g(3) * t40 + t127 * t82 + t94 * t182 + t31 * t46 - t1, 0, 0;];
tau_reg = t5;
