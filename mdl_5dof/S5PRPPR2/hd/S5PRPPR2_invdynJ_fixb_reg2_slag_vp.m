% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPPR2
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:49
% EndTime: 2019-12-05 15:24:51
% DurationCPUTime: 1.36s
% Computational Cost: add. (1354->207), mult. (2829->274), div. (0->0), fcn. (2228->14), ass. (0->125)
t102 = cos(pkin(7));
t99 = sin(pkin(7));
t126 = g(1) * t102 + g(2) * t99;
t96 = qJ(2) + pkin(8);
t87 = sin(t96);
t89 = cos(t96);
t171 = -g(3) * t89 + t126 * t87;
t101 = cos(pkin(8));
t105 = sin(qJ(2));
t107 = cos(qJ(2));
t98 = sin(pkin(8));
t60 = t101 * t105 + t98 * t107;
t46 = t60 * qJD(1);
t168 = t46 * qJD(2) + t171;
t138 = qJD(1) * qJD(2);
t170 = t105 * qJDD(1) + t107 * t138;
t150 = qJDD(2) * pkin(2);
t90 = t107 * qJDD(1);
t57 = -t105 * t138 + t150 + t90;
t135 = -t101 * t57 + t170 * t98;
t127 = qJDD(4) + t135;
t149 = qJDD(2) * pkin(3);
t20 = t127 - t149;
t160 = t101 * pkin(2);
t81 = -pkin(3) - t160;
t172 = -qJDD(2) * t81 + t168 - t20;
t140 = t107 * qJD(1);
t143 = qJD(1) * t105;
t73 = t98 * t143;
t50 = t101 * t140 - t73;
t145 = qJD(4) - t50;
t100 = cos(pkin(9));
t106 = cos(qJ(5));
t104 = sin(qJ(5));
t97 = sin(pkin(9));
t151 = t104 * t97;
t120 = t106 * t100 - t151;
t22 = t120 * t60;
t61 = t104 * t100 + t106 * t97;
t51 = t61 * qJD(2);
t164 = g(3) * t87;
t112 = -t126 * t89 - t164;
t166 = t51 ^ 2;
t76 = t98 * pkin(2) + qJ(4);
t163 = pkin(6) + t76;
t162 = pkin(2) * t105;
t142 = qJD(2) * t100;
t130 = t106 * t142;
t133 = qJD(2) * t151;
t47 = -t130 + t133;
t159 = t51 * t47;
t95 = pkin(9) + qJ(5);
t86 = sin(t95);
t158 = t99 * t86;
t88 = cos(t95);
t157 = t99 * t88;
t55 = t163 * t97;
t56 = t163 * t100;
t25 = -t104 * t56 - t106 * t55;
t156 = t25 * qJD(5) + t145 * t120;
t26 = -t104 * t55 + t106 * t56;
t155 = -t26 * qJD(5) - t145 * t61;
t28 = t170 * t101 + t98 * t57;
t19 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t28;
t13 = t97 * qJDD(3) + t100 * t19;
t137 = t100 * qJDD(2);
t139 = t97 * qJDD(2);
t124 = t104 * t139 - t106 * t137;
t54 = t61 * qJD(5);
t30 = qJD(2) * t54 + t124;
t53 = t120 * qJD(5);
t154 = -t61 * t30 - t53 * t47;
t67 = qJD(2) * pkin(2) + t140;
t41 = t101 * t143 + t98 * t67;
t39 = qJD(2) * qJ(4) + t41;
t34 = t97 * qJD(3) + t100 * t39;
t93 = t97 ^ 2;
t94 = t100 ^ 2;
t153 = t93 + t94;
t152 = t102 * t89;
t119 = t101 * t107 - t98 * t105;
t49 = t119 * qJD(2);
t146 = t49 * qJD(2);
t144 = qJDD(1) - g(3);
t134 = qJD(5) * t130 + t104 * t137 + t106 * t139;
t80 = t100 * pkin(4) + pkin(3);
t132 = -g(1) * t99 + g(2) * t102;
t83 = t100 * qJDD(3);
t8 = t83 + (-pkin(6) * qJDD(2) - t19) * t97;
t9 = pkin(6) * t137 + t13;
t131 = -t104 * t9 + t106 * t8;
t40 = t101 * t67 - t73;
t128 = qJDD(2) * t153;
t125 = qJD(4) - t40;
t123 = t104 * t8 + t106 * t9;
t29 = qJD(5) * t133 - t134;
t122 = t120 * t29 + t51 * t54;
t121 = t80 * qJDD(2);
t85 = t100 * qJD(3);
t23 = t85 + (-pkin(6) * qJD(2) - t39) * t97;
t24 = pkin(6) * t142 + t34;
t5 = -t104 * t24 + t106 * t23;
t6 = t104 * t23 + t106 * t24;
t45 = t60 * qJD(2);
t118 = t45 * qJD(2) - qJDD(2) * t119;
t21 = t61 * t60;
t111 = -g(3) * t107 + t126 * t105;
t109 = t127 - t171;
t108 = qJD(2) ^ 2;
t103 = -pkin(6) - qJ(4);
t92 = t107 * pkin(2);
t63 = -t80 - t160;
t44 = t47 ^ 2;
t38 = -qJD(2) * pkin(3) + t125;
t36 = -t80 * qJD(2) + t125;
t33 = -t97 * t39 + t85;
t32 = -t54 * qJD(5) + qJDD(5) * t120;
t31 = t53 * qJD(5) + t61 * qJDD(5);
t15 = -t121 + t127;
t12 = -t97 * t19 + t83;
t4 = -qJD(5) * t22 - t61 * t49;
t3 = -qJD(5) * t21 + t120 * t49;
t2 = -t6 * qJD(5) + t131;
t1 = t5 * qJD(5) + t123;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t144, 0, 0, 0, 0, 0, 0, t107 * qJDD(2) - t108 * t105, -qJDD(2) * t105 - t108 * t107, 0, -g(3) + (t105 ^ 2 + t107 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t118, -t60 * qJDD(2) - t146, 0, -t119 * t135 + t28 * t60 - t40 * t45 + t41 * t49 - g(3), 0, 0, 0, 0, 0, 0, -t118 * t100, t118 * t97, t60 * t128 + t153 * t146, -t20 * t119 + t38 * t45 - g(3) + (-t12 * t60 - t33 * t49) * t97 + (t13 * t60 + t34 * t49) * t100, 0, 0, 0, 0, 0, 0, t4 * qJD(5) - t21 * qJDD(5) - t119 * t30 + t45 * t47, -t3 * qJD(5) - t22 * qJDD(5) + t119 * t29 + t45 * t51, -t21 * t29 - t22 * t30 - t3 * t47 - t4 * t51, t1 * t22 - t119 * t15 - t2 * t21 + t6 * t3 + t36 * t45 + t5 * t4 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t90 + t111, -t144 * t105 + t126 * t107, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t101 * t150 - t135 + t168, t50 * qJD(2) - t98 * t150 - t112 - t28, 0, t40 * t46 - t41 * t50 + (-t101 * t135 + t28 * t98 + t111) * pkin(2), t93 * qJDD(2), 0.2e1 * t97 * t137, 0, t94 * qJDD(2), 0, 0, t172 * t100, -t172 * t97, t145 * qJD(2) * t153 + t13 * t100 - t12 * t97 + t76 * t128 + t112, t20 * t81 - t38 * t46 - g(3) * (t89 * pkin(3) + t87 * qJ(4) + t92) + (-t12 * t76 - t145 * t33) * t97 + (t13 * t76 + t145 * t34) * t100 + t126 * (pkin(3) * t87 - qJ(4) * t89 + t162), -t29 * t61 + t51 * t53, -t122 + t154, t31, -t120 * t30 + t47 * t54, t32, 0, t155 * qJD(5) + t25 * qJDD(5) - t15 * t120 + t171 * t88 + t63 * t30 + t36 * t54 - t46 * t47, -t156 * qJD(5) - t26 * qJDD(5) + t15 * t61 - t171 * t86 - t63 * t29 + t36 * t53 - t46 * t51, t1 * t120 - t155 * t51 - t156 * t47 - t2 * t61 + t25 * t29 - t26 * t30 - t5 * t53 - t6 * t54 + t112, t1 * t26 + t2 * t25 + t15 * t63 - t36 * t46 - g(3) * (-t87 * t103 + t89 * t80 + t92) + t156 * t6 + t155 * t5 + t126 * (t103 * t89 + t80 * t87 + t162); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + t132, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t100 + t13 * t97 + t132, 0, 0, 0, 0, 0, 0, t32, -t31, t122 + t154, t1 * t61 + t120 * t2 - t5 * t54 + t6 * t53 + t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, t139, -t153 * t108, -t149 + (-t34 * t100 + t33 * t97) * qJD(2) + t109, 0, 0, 0, 0, 0, 0, 0.2e1 * t51 * qJD(5) + t124, (-t47 - t133) * qJD(5) + t134, -t44 - t166, t6 * t47 + t5 * t51 + t109 - t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, -t44 + t166, (t47 - t133) * qJD(5) + t134, -t159, -t124, qJDD(5), -t36 * t51 - g(1) * (-t86 * t152 + t157) - g(2) * (-t102 * t88 - t89 * t158) + t86 * t164 + t131, t36 * t47 - g(1) * (-t88 * t152 - t158) - g(2) * (t102 * t86 - t89 * t157) + t88 * t164 - t123, 0, 0;];
tau_reg = t7;
