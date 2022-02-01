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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
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
% StartTime: 2022-01-20 09:12:44
% EndTime: 2022-01-20 09:12:51
% DurationCPUTime: 1.74s
% Computational Cost: add. (1859->227), mult. (3932->326), div. (0->0), fcn. (2778->14), ass. (0->145)
t104 = sin(pkin(8));
t110 = sin(qJ(5));
t106 = cos(pkin(9));
t112 = cos(qJ(5));
t162 = t112 * t106;
t146 = t104 * t162;
t103 = sin(pkin(9));
t152 = qJDD(1) * t103;
t71 = t112 * t103 + t110 * t106;
t187 = t71 * qJD(5);
t19 = -qJDD(1) * t146 + (qJD(1) * t187 + t110 * t152) * t104;
t107 = cos(pkin(8));
t82 = t107 * qJD(1) - qJD(5);
t161 = qJD(1) * t104;
t144 = t103 * t161;
t136 = t110 * t144;
t143 = qJD(5) * t162;
t20 = -qJD(5) * t136 + t104 * (qJD(1) * t143 + t71 * qJDD(1));
t102 = qJ(1) + pkin(7);
t92 = sin(t102);
t184 = g(1) * t92;
t94 = cos(t102);
t87 = g(2) * t94;
t145 = t87 - t184;
t166 = pkin(1) * qJDD(1);
t105 = sin(pkin(7));
t85 = t105 * pkin(1) + qJ(3);
t68 = qJD(1) * qJD(3) + qJDD(1) * t85;
t182 = pkin(6) * t104;
t147 = t106 * t182;
t124 = -t107 * pkin(4) - t147;
t158 = qJD(4) * t104;
t108 = cos(pkin(7));
t89 = -t108 * pkin(1) - pkin(2);
t67 = -t107 * pkin(3) - t104 * qJ(4) + t89;
t29 = -qJD(1) * t158 + t67 * qJDD(1) + qJDD(3);
t42 = t104 * qJDD(2) + t107 * t68;
t14 = -t103 * t42 + t106 * t29;
t10 = t124 * qJDD(1) + t14;
t150 = t104 * qJDD(1);
t141 = t103 * t150;
t15 = t103 * t29 + t106 * t42;
t11 = -pkin(6) * t141 + t15;
t43 = t67 * qJD(1) + qJD(3);
t78 = t85 * qJD(1);
t56 = t104 * qJD(2) + t107 * t78;
t17 = -t103 * t56 + t106 * t43;
t13 = t124 * qJD(1) + t17;
t18 = t103 * t43 + t106 * t56;
t16 = -pkin(6) * t144 + t18;
t126 = t110 * t16 - t112 * t13;
t1 = -t126 * qJD(5) + t110 * t10 + t112 * t11;
t49 = qJD(1) * t146 - t136;
t185 = t49 ^ 2;
t183 = pkin(4) * t103;
t181 = g(3) * t104;
t119 = qJD(1) * t71;
t46 = t104 * t119;
t180 = t46 * t82;
t179 = t49 * t46;
t178 = t49 * t82;
t52 = t104 * t187;
t163 = t110 * t103;
t70 = t162 - t163;
t58 = t70 * t104;
t177 = -t58 * t20 + t52 * t46;
t53 = (-qJD(5) * t163 + t143) * t104;
t57 = t71 * t104;
t149 = t107 * qJDD(1);
t81 = -qJDD(5) + t149;
t176 = t53 * t82 + t57 * t81;
t175 = -t107 * t119 + t187;
t174 = t82 * t70;
t171 = t107 * t85;
t28 = t103 * t67 + t106 * t171;
t97 = t103 ^ 2;
t99 = t106 ^ 2;
t173 = -t97 - t99;
t172 = t107 * t20;
t170 = t107 * t92;
t169 = t107 * t94;
t168 = t19 * t107;
t100 = t107 ^ 2;
t98 = t104 ^ 2;
t167 = t100 + t98;
t165 = t104 * (-pkin(6) - qJ(4));
t114 = qJD(1) ^ 2;
t164 = t107 * t114;
t160 = qJD(3) * t104;
t159 = qJD(3) * t107;
t156 = qJDD(1) * t89;
t154 = t98 * qJDD(1);
t151 = qJDD(1) * t106;
t113 = cos(qJ(1));
t148 = t113 * pkin(1) + t94 * pkin(2) + t92 * qJ(3);
t111 = sin(qJ(1));
t142 = -t111 * pkin(1) + t94 * qJ(3);
t139 = t106 * t149;
t138 = t167 * t114;
t55 = t107 * qJD(2) - t104 * t78;
t41 = t107 * qJDD(2) - t104 * t68;
t137 = 0.2e1 * t104 * t149;
t134 = -g(1) * t94 - g(2) * t92;
t132 = g(1) * t111 - g(2) * t113;
t131 = -t57 * t19 + t49 * t53;
t130 = -t52 * t82 + t58 * t81;
t54 = qJD(4) - t55;
t128 = -t41 * t104 + t42 * t107;
t127 = t55 * t104 - t56 * t107;
t6 = t110 * t13 + t112 * t16;
t60 = t106 * t67;
t21 = -t147 + t60 + (-t103 * t85 - pkin(4)) * t107;
t23 = -t103 * t182 + t28;
t7 = -t110 * t23 + t112 * t21;
t8 = t110 * t21 + t112 * t23;
t36 = qJDD(4) - t41;
t125 = -t92 * pkin(2) + t142;
t74 = qJDD(3) + t156;
t123 = t156 + t74 + t87;
t122 = g(3) * t107 + t36;
t27 = -t103 * t171 + t60;
t63 = -t103 * t159 - t106 * t158;
t121 = -t63 * qJD(1) - t27 * qJDD(1) - t14;
t64 = -t103 * t158 + t106 * t159;
t120 = t64 * qJD(1) + t28 * qJDD(1) + t15;
t2 = -t6 * qJD(5) + t112 * t10 - t110 * t11;
t117 = t36 * t104 + t68 * t98 + t134;
t101 = pkin(9) + qJ(5);
t93 = cos(t101);
t91 = sin(t101);
t90 = t100 * qJDD(1);
t88 = t106 * pkin(4) + pkin(3);
t80 = pkin(4) * t141;
t79 = t104 * t184;
t61 = (t85 + t183) * t104;
t45 = t46 ^ 2;
t40 = t93 * t169 + t92 * t91;
t39 = -t91 * t169 + t92 * t93;
t38 = -t93 * t170 + t94 * t91;
t37 = t91 * t170 + t94 * t93;
t31 = pkin(4) * t144 + t54;
t26 = t36 + t80;
t4 = -t8 * qJD(5) - t110 * t64 + t112 * t63;
t3 = t7 * qJD(5) + t110 * t63 + t112 * t64;
t5 = [0, 0, 0, 0, 0, qJDD(1), t132, g(1) * t113 + g(2) * t111, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t108 * t166 - t145, -0.2e1 * t105 * t166 - t134, 0, (t132 + (t105 ^ 2 + t108 ^ 2) * t166) * pkin(1), t154, t137, 0, t90, 0, 0, (-t123 + t184) * t107, t123 * t104 - t79, t167 * t68 + t128 + t134, -g(1) * t125 - g(2) * t148 - t127 * qJD(3) + t128 * t85 + t74 * t89, t99 * t154, -0.2e1 * t98 * t103 * t151, -0.2e1 * t104 * t139, t97 * t154, t103 * t137, t90, (-t106 * t145 + t121) * t107 + t117 * t103, (t103 * t145 + t120) * t107 + t117 * t106, t79 + (-t120 * t103 + t121 * t106 - t87) * t104, t15 * t28 + t18 * t64 + t14 * t27 + t17 * t63 - g(1) * (-pkin(3) * t170 + t125) - g(2) * (pkin(3) * t169 + t148) + (-qJ(4) * t145 + t54 * qJD(3) + t36 * t85) * t104, -t19 * t58 - t49 * t52, -t131 + t177, -t130 + t168, t20 * t57 + t46 * t53, t172 + t176, t81 * t107, -g(1) * t38 - g(2) * t40 - t2 * t107 + t46 * t160 + t61 * t20 + t26 * t57 + t31 * t53 - t4 * t82 - t7 * t81, -g(1) * t37 - g(2) * t39 + t1 * t107 + t49 * t160 - t61 * t19 + t26 * t58 + t3 * t82 - t31 * t52 + t8 * t81, -t1 * t57 - t104 * t87 - t126 * t52 + t7 * t19 - t2 * t58 - t8 * t20 - t3 * t46 - t4 * t49 - t6 * t53 + t79, t1 * t8 + t6 * t3 + t2 * t7 - t126 * t4 + t26 * t61 + t31 * t160 - g(1) * (t94 * t183 + t142) - g(2) * (-t94 * t165 + t88 * t169 + t148) + (-g(1) * (-t107 * t88 - pkin(2) + t165) - g(2) * t183) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t104 + t41 * t107 - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36 * t107 - g(3) + (-t103 * t14 + t106 * t15) * t104, 0, 0, 0, 0, 0, 0, -t172 + t176, t130 + t168, t131 + t177, t1 * t58 - t26 * t107 + t126 * t53 - t2 * t57 - t6 * t52 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t150, -t138, t127 * qJD(1) + t145 + t74, 0, 0, 0, 0, 0, 0, -t103 * t138 - t139, t103 * t149 - t106 * t138, t173 * t150, t15 * t103 + t14 * t106 + (-t104 * t54 + (t103 * t17 - t106 * t18) * t107) * qJD(1) + t145, 0, 0, 0, 0, 0, 0, -t46 * t161 + t175 * t82 - t70 * t81, -t49 * t161 - t174 * t82 + t71 * t81, t174 * t46 + t175 * t49 + t70 * t19 - t71 * t20, t1 * t71 + t126 * t175 - t31 * t161 - t174 * t6 + t2 * t70 + t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t106 * t164 + t152) * t104, (t103 * t164 + t151) * t104, t173 * t98 * t114, ((t103 * t18 + t106 * t17) * qJD(1) + t134) * t104 + t122, 0, 0, 0, 0, 0, 0, t20 - t178, -t19 + t180, -t45 - t185, t104 * t134 - t126 * t49 + t6 * t46 + t122 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t45 + t185, -t19 - t180, -t179, -t178 - t20, -t81, -g(1) * t39 + g(2) * t37 + t91 * t181 - t31 * t49 - t6 * t82 + t2, g(1) * t40 - g(2) * t38 + t126 * t82 + t93 * t181 + t31 * t46 - t1, 0, 0;];
tau_reg = t5;
