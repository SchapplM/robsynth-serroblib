% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:38
% EndTime: 2018-11-14 13:45:39
% DurationCPUTime: 0.60s
% Computational Cost: add. (533->175), mult. (1352->241), div. (0->0), fcn. (981->10), ass. (0->141)
t82 = sin(pkin(4));
t78 = t82 ^ 2;
t83 = cos(pkin(6));
t160 = t78 * t83;
t84 = cos(pkin(4));
t159 = pkin(1) * t84;
t85 = sin(qJ(1));
t158 = g(1) * t85;
t81 = sin(pkin(6));
t77 = t81 ^ 2;
t157 = t77 * t78;
t156 = t81 * t82;
t155 = t82 * t83;
t86 = cos(qJ(1));
t154 = t82 * t86;
t87 = qJD(1) ^ 2;
t153 = t84 * t87;
t140 = qJDD(1) * t82;
t122 = qJ(2) * t140;
t138 = qJD(1) * qJD(2);
t120 = t82 * t138;
t49 = t81 * t120;
t152 = t81 * t122 + t49;
t146 = qJD(1) * t82;
t127 = qJ(2) * t146;
t145 = qJD(1) * t84;
t129 = t81 * t145;
t31 = pkin(1) * t129 + t83 * t127;
t150 = qJ(2) * t82;
t36 = t83 * t150 + t81 * t159;
t151 = t86 * pkin(1) + t85 * t150;
t149 = t84 * qJ(3);
t100 = pkin(3) * t155 + t149;
t8 = t100 * qJD(1) + qJD(4) + t31;
t148 = -qJD(4) - t8;
t144 = qJD(2) * t82;
t41 = t84 * qJD(3) + t83 * t144;
t147 = qJD(1) * t41;
t52 = t81 * t127;
t143 = qJD(3) + t52;
t142 = qJDD(1) * t78;
t141 = qJDD(1) * t81;
t139 = t84 * qJDD(1);
t137 = qJD(1) * qJD(3);
t76 = pkin(4) - pkin(6);
t68 = cos(t76) / 0.2e1;
t75 = pkin(4) + pkin(6);
t72 = cos(t75);
t136 = t68 + t72 / 0.2e1;
t67 = -sin(t76) / 0.2e1;
t71 = sin(t75);
t135 = t71 / 0.2e1 + t67;
t134 = t83 * t159;
t133 = t81 * t160;
t132 = t83 * t153;
t131 = pkin(1) * t139;
t51 = t83 * t120;
t12 = t83 * t122 + t81 * t131 + t51;
t130 = t81 * t146;
t128 = -pkin(1) * t83 - pkin(2);
t126 = t81 * t140;
t125 = t83 * t140;
t124 = qJDD(3) + t152;
t123 = -t85 * pkin(1) + t86 * t150;
t121 = t78 * t138;
t119 = t81 * t137;
t118 = -qJ(3) * t81 - pkin(1);
t117 = t87 * t133;
t110 = t128 * t84;
t61 = t81 * t150;
t25 = t61 + t110;
t104 = t128 * qJDD(1);
t6 = t84 * t104 + t124;
t116 = qJDD(1) * t25 + t6;
t101 = -pkin(2) * t83 + t118;
t26 = t101 * t82;
t92 = t101 * qJDD(1);
t7 = qJDD(2) + (t92 - t119) * t82;
t115 = qJDD(1) * t26 + t7;
t114 = qJDD(1) * t133;
t113 = t84 * t126;
t112 = t84 * t125;
t111 = -qJ(4) + t128;
t4 = -qJ(3) * t139 - t84 * t137 - t12;
t20 = -t86 * t136 + t85 * t81;
t22 = t85 * t136 + t86 * t81;
t109 = g(1) * t20 - g(2) * t22;
t21 = t86 * t135 + t85 * t83;
t23 = -t85 * t135 + t86 * t83;
t108 = g(1) * t21 - g(2) * t23;
t107 = g(1) * t86 + g(2) * t85;
t106 = g(2) * t154 - g(3) * t84 + qJDD(2);
t105 = (qJD(1) * t134 - t52) * t81 - t31 * t83;
t103 = t23 * pkin(2) + t22 * qJ(3) + t151;
t102 = -qJD(3) * t81 - qJD(4) * t83;
t99 = t111 * qJDD(1);
t58 = -pkin(1) * t140 + qJDD(2);
t98 = pkin(1) * t142 - t58 * t82;
t54 = pkin(3) * t126;
t1 = t54 + (-qJD(1) * qJD(4) + t99) * t84 + t124;
t89 = pkin(3) * t156 + t111 * t84;
t10 = t61 + t89;
t39 = -t84 * qJD(4) + t81 * t144;
t97 = qJD(1) * t39 + qJDD(1) * t10 + t1;
t13 = t100 + t36;
t3 = pkin(3) * t125 + qJDD(4) - t4;
t96 = qJDD(1) * t13 + t147 + t3;
t24 = -t149 - t36;
t95 = -qJDD(1) * t24 + t147 - t4;
t94 = -t21 * pkin(2) - t20 * qJ(3) + t123;
t93 = (-pkin(2) - qJ(4)) * t83 + t118;
t14 = t93 * t82;
t90 = t93 * qJDD(1);
t2 = qJDD(2) + (t102 * qJD(1) + t90) * t82;
t34 = t102 * t82;
t91 = (-qJD(1) * t34 - qJDD(1) * t14 - t2) * t82;
t88 = -g(1) * t22 - g(2) * t20 - g(3) * (-t71 / 0.2e1 + t67) + t124;
t80 = t84 ^ 2;
t79 = t83 ^ 2;
t73 = t80 * qJDD(1);
t60 = t79 * t142;
t59 = t77 * t142;
t47 = t77 * t121;
t46 = t153 * t156;
t45 = -0.2e1 * t112;
t44 = 0.2e1 * t113;
t43 = 0.2e1 * t114;
t40 = (-t80 - t157) * t87;
t35 = -t61 + t134;
t33 = (-t77 - t79) * t87 * t78;
t32 = t117 + t139;
t29 = -t46 + t125;
t28 = (-t132 + t141) * t82;
t27 = (t132 + t141) * t82;
t17 = -qJ(3) * t145 - t31;
t16 = qJD(1) * t26 + qJD(2);
t15 = qJD(1) * t110 + t143;
t11 = t83 * t131 - t152;
t9 = qJD(1) * t14 + qJD(2);
t5 = t89 * qJD(1) + t143;
t18 = [0, 0, 0, 0, 0, qJDD(1), -g(2) * t86 + t158, t107, 0, 0, t59, t43, t44, t60, 0.2e1 * t112, t73, t98 * t83 + (qJDD(1) * t35 + t11 - t49) * t84 + t108, -t98 * t81 + (-qJDD(1) * t36 - t12 - t51) * t84 - t109, t79 * t121 + t47 + (-t11 * t81 + t12 * t83 + (-t35 * t81 + t36 * t83) * qJDD(1) - t107) * t82, t12 * t36 + t11 * t35 - g(1) * t123 - g(2) * t151 + (-t58 * pkin(1) - t105 * qJD(2)) * t82, t73, -0.2e1 * t113, t45, t59, t43, t60, t47 + (t116 * t81 + t95 * t83 - t107) * t82, -t119 * t160 + t116 * t84 + (qJD(2) * t129 + t115 * t83) * t82 - t108, -t115 * t156 + t137 * t157 + t95 * t84 + t109, t7 * t26 + t4 * t24 - t17 * t41 + t6 * t25 - g(1) * t94 - g(2) * t103 + (qJD(2) * t15 - qJD(3) * t16) * t156, t73, t45, t44, t60, -0.2e1 * t114, t59 (t97 * t81 + t96 * t83 - t107) * t82, t81 * t91 + t96 * t84 + t109, t83 * t91 - t97 * t84 + t108, t2 * t14 + t9 * t34 + t1 * t10 + t5 * t39 + t3 * t13 + t8 * t41 - g(1) * (pkin(3) * t154 - t21 * qJ(4) + t94) - g(2) * (t85 * t82 * pkin(3) + t23 * qJ(4) + t103); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t27, t33 (-pkin(1) * qJDD(1) + t105 * qJD(1) - t158) * t82 + t106, 0, 0, 0, 0, 0, 0, t33, t29, -t27 (-t158 + t92 + (t17 * t83 + (-qJD(3) - t15) * t81) * qJD(1)) * t82 + t106, 0, 0, 0, 0, 0, 0, t33, -t27, -t29 (-t158 + t90 + (t148 * t83 + (-qJD(3) - t5) * t81) * qJD(1)) * t82 + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t32, t40, t16 * t130 + (qJD(1) * t17 + t104) * t84 + t88, 0, 0, 0, 0, 0, 0, t28, t40, -t32, t9 * t130 + t54 + (t148 * qJD(1) + t99) * t84 + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 + t125, -t117 + t139 (-t78 * t79 - t80) * t87, -g(1) * t23 - g(2) * t21 - g(3) * (t68 - t72 / 0.2e1) + (t9 * t155 + t5 * t84) * qJD(1) + t3;];
tau_reg  = t18;
