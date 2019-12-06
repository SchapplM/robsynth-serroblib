% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPPR1
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:16
% EndTime: 2019-12-05 15:22:19
% DurationCPUTime: 1.44s
% Computational Cost: add. (1526->218), mult. (3385->309), div. (0->0), fcn. (2491->10), ass. (0->139)
t101 = sin(pkin(8));
t105 = sin(qJ(5));
t102 = cos(pkin(9));
t106 = cos(qJ(5));
t151 = t106 * t102;
t137 = t101 * t151;
t100 = sin(pkin(9));
t142 = qJDD(2) * t100;
t64 = t106 * t100 + t105 * t102;
t177 = t64 * qJD(5);
t16 = -qJDD(2) * t137 + (qJD(2) * t177 + t105 * t142) * t101;
t103 = cos(pkin(8));
t146 = t103 * qJD(2);
t77 = -qJD(5) + t146;
t150 = qJD(2) * t101;
t134 = t100 * t150;
t127 = t105 * t134;
t133 = qJD(5) * t151;
t17 = -qJD(5) * t127 + t101 * (qJD(2) * t133 + t64 * qJDD(2));
t98 = pkin(7) + qJ(2);
t89 = sin(t98);
t174 = g(1) * t89;
t91 = cos(t98);
t83 = g(2) * t91;
t178 = -t83 + t174;
t143 = qJD(2) * qJD(3);
t115 = qJ(3) * qJDD(2) + t143;
t172 = pkin(6) * t101;
t138 = t102 * t172;
t117 = -t103 * pkin(4) - t138;
t147 = qJD(4) * t101;
t69 = -t103 * pkin(3) - t101 * qJ(4) - pkin(2);
t33 = -qJD(2) * t147 + t69 * qJDD(2) + qJDD(3);
t53 = t101 * qJDD(1) + t115 * t103;
t14 = -t100 * t53 + t102 * t33;
t10 = t117 * qJDD(2) + t14;
t54 = t69 * qJD(2) + qJD(3);
t67 = qJ(3) * t146 + t101 * qJD(1);
t20 = -t100 * t67 + t102 * t54;
t13 = t117 * qJD(2) + t20;
t21 = t100 * t54 + t102 * t67;
t18 = -pkin(6) * t134 + t21;
t118 = t105 * t18 - t106 * t13;
t140 = t101 * qJDD(2);
t131 = t100 * t140;
t15 = t100 * t33 + t102 * t53;
t12 = -pkin(6) * t131 + t15;
t1 = -t118 * qJD(5) + t105 * t10 + t106 * t12;
t42 = qJD(2) * t137 - t127;
t175 = t42 ^ 2;
t173 = pkin(4) * t100;
t171 = g(3) * t101;
t111 = qJD(2) * t64;
t39 = t101 * t111;
t170 = t39 * t77;
t169 = t42 * t39;
t168 = t42 * t77;
t45 = t101 * t177;
t152 = t105 * t100;
t63 = t151 - t152;
t51 = t63 * t101;
t167 = -t51 * t17 + t45 * t39;
t46 = (-qJD(5) * t152 + t133) * t101;
t50 = t64 * t101;
t139 = t103 * qJDD(2);
t76 = -qJDD(5) + t139;
t166 = t46 * t77 + t50 * t76;
t165 = -t103 * t111 + t177;
t164 = t77 * t63;
t156 = qJ(3) * t103;
t35 = t100 * t69 + t102 * t156;
t163 = t91 * pkin(2) + t89 * qJ(3);
t93 = t100 ^ 2;
t95 = t102 ^ 2;
t162 = -t93 - t95;
t94 = t101 ^ 2;
t96 = t103 ^ 2;
t161 = t94 + t96;
t160 = t103 * t17;
t159 = t103 * t89;
t158 = t103 * t91;
t157 = t16 * t103;
t155 = qJDD(2) * pkin(2);
t154 = t101 * (-pkin(6) - qJ(4));
t107 = qJD(2) ^ 2;
t153 = t103 * t107;
t149 = qJD(3) * t101;
t148 = qJD(3) * t103;
t145 = t94 * qJDD(2);
t141 = qJDD(2) * t102;
t79 = t91 * qJ(3);
t136 = -t89 * pkin(2) + t79;
t132 = t161 * t107;
t129 = t102 * t139;
t66 = -qJ(3) * t150 + t103 * qJD(1);
t128 = 0.2e1 * t101 * t139;
t86 = qJDD(3) - t155;
t125 = g(1) * t91 + g(2) * t89;
t123 = -t50 * t16 + t42 * t46;
t122 = -t45 * t77 + t51 * t76;
t52 = -qJ(3) * t140 + t103 * qJDD(1) - t101 * t143;
t62 = qJD(4) - t66;
t120 = -t52 * t101 + t53 * t103;
t119 = t66 * t101 - t67 * t103;
t6 = t105 * t13 + t106 * t18;
t61 = t102 * t69;
t22 = -t138 + t61 + (-qJ(3) * t100 - pkin(4)) * t103;
t23 = -t100 * t172 + t35;
t7 = -t105 * t23 + t106 * t22;
t8 = t105 * t22 + t106 * t23;
t116 = -t155 + t86 + t83;
t49 = qJDD(4) - t52;
t34 = -t100 * t156 + t61;
t55 = -t100 * t148 - t102 * t147;
t114 = -t55 * qJD(2) - t34 * qJDD(2) - t14;
t56 = -t100 * t147 + t102 * t148;
t113 = t56 * qJD(2) + t35 * qJDD(2) + t15;
t112 = g(3) * t103 + t49;
t2 = -t6 * qJD(5) + t106 * t10 - t105 * t12;
t109 = t49 * t101 + t115 * t94 - t125;
t99 = qJDD(1) - g(3);
t97 = pkin(9) + qJ(5);
t90 = cos(t97);
t88 = sin(t97);
t87 = t96 * qJDD(2);
t85 = t102 * pkin(4) + pkin(3);
t74 = pkin(4) * t131;
t73 = t101 * t174;
t65 = (qJ(3) + t173) * t101;
t38 = pkin(4) * t134 + t62;
t36 = t39 ^ 2;
t32 = t90 * t158 + t89 * t88;
t31 = -t88 * t158 + t89 * t90;
t30 = -t90 * t159 + t91 * t88;
t29 = t88 * t159 + t91 * t90;
t25 = t49 + t74;
t4 = -t8 * qJD(5) - t105 * t56 + t106 * t55;
t3 = t7 * qJD(5) + t105 * t55 + t106 * t56;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t101 + t52 * t103 - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49 * t103 - g(3) + (-t100 * t14 + t102 * t15) * t101, 0, 0, 0, 0, 0, 0, -t160 + t166, t122 + t157, t123 + t167, t1 * t51 - t25 * t103 + t118 * t46 - t2 * t50 - t6 * t45 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t178, t125, 0, 0, t145, t128, 0, t87, 0, 0, (-t116 + t174) * t103, t116 * t101 - t73, t115 * t161 + t120 - t125, -t86 * pkin(2) - g(1) * t136 - g(2) * t163 + t120 * qJ(3) - t119 * qJD(3), t95 * t145, -0.2e1 * t94 * t100 * t141, -0.2e1 * t101 * t129, t93 * t145, t100 * t128, t87, (t102 * t178 + t114) * t103 + t109 * t100, (-t100 * t178 + t113) * t103 + t109 * t102, t73 + (-t113 * t100 + t114 * t102 - t83) * t101, t15 * t35 + t21 * t56 + t14 * t34 + t20 * t55 - g(1) * (-pkin(3) * t159 + t136) - g(2) * (pkin(3) * t158 + t163) + (t49 * qJ(3) + qJ(4) * t178 + t62 * qJD(3)) * t101, -t16 * t51 - t42 * t45, -t123 + t167, -t122 + t157, t17 * t50 + t39 * t46, t160 + t166, t76 * t103, -g(1) * t30 - g(2) * t32 - t2 * t103 + t149 * t39 + t65 * t17 + t25 * t50 + t38 * t46 - t4 * t77 - t7 * t76, -g(1) * t29 - g(2) * t31 + t1 * t103 + t149 * t42 - t65 * t16 + t25 * t51 + t3 * t77 - t38 * t45 + t8 * t76, -t1 * t50 - t101 * t83 - t118 * t45 + t7 * t16 - t8 * t17 - t2 * t51 - t3 * t39 - t4 * t42 - t6 * t46 + t73, t1 * t8 + t6 * t3 + t2 * t7 - t118 * t4 + t25 * t65 + t38 * t149 - g(1) * (t91 * t173 + t79) - g(2) * (-t91 * t154 + t85 * t158 + t163) + (-g(1) * (-t103 * t85 - pkin(2) + t154) - g(2) * t173) * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, t140, -t132, t119 * qJD(2) - t178 + t86, 0, 0, 0, 0, 0, 0, -t100 * t132 - t129, t100 * t139 - t102 * t132, t162 * t140, t15 * t100 + t14 * t102 + (-t101 * t62 + (t100 * t20 - t102 * t21) * t103) * qJD(2) - t178, 0, 0, 0, 0, 0, 0, -t39 * t150 + t165 * t77 - t63 * t76, -t42 * t150 - t164 * t77 + t64 * t76, t63 * t16 + t164 * t39 + t165 * t42 - t64 * t17, t1 * t64 + t118 * t165 - t38 * t150 - t164 * t6 + t2 * t63 - t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t102 * t153 + t142) * t101, (t100 * t153 + t141) * t101, t162 * t94 * t107, ((t100 * t21 + t102 * t20) * qJD(2) - t125) * t101 + t112, 0, 0, 0, 0, 0, 0, t17 - t168, -t16 + t170, -t36 - t175, -t101 * t125 - t118 * t42 + t6 * t39 + t112 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, -t36 + t175, -t16 - t170, -t169, -t168 - t17, -t76, -g(1) * t31 + g(2) * t29 + t88 * t171 - t38 * t42 - t6 * t77 + t2, g(1) * t32 - g(2) * t30 + t118 * t77 + t90 * t171 + t38 * t39 - t1, 0, 0;];
tau_reg = t5;
