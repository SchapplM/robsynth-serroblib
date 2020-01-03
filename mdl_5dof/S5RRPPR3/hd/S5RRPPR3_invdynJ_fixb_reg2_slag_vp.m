% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:40
% EndTime: 2019-12-31 19:26:41
% DurationCPUTime: 0.75s
% Computational Cost: add. (1358->190), mult. (1974->212), div. (0->0), fcn. (1086->12), ass. (0->132)
t134 = pkin(1) * qJD(2);
t118 = qJD(1) * t134;
t79 = sin(qJ(2));
t129 = qJDD(1) * t79;
t82 = cos(qJ(2));
t158 = pkin(1) * t129 + t82 * t118;
t75 = qJ(1) + qJ(2);
t65 = pkin(8) + t75;
t56 = sin(t65);
t57 = cos(t65);
t140 = g(1) * t56 - g(2) * t57;
t135 = pkin(1) * qJD(1);
t121 = t79 * t135;
t120 = t82 * t135;
t71 = qJD(1) + qJD(2);
t36 = t71 * pkin(2) + t120;
t76 = sin(pkin(8));
t77 = cos(pkin(8));
t20 = t77 * t121 + t76 * t36;
t17 = t71 * qJ(4) + t20;
t145 = t17 * t71;
t157 = -t140 - t145;
t103 = -g(1) * t57 - g(2) * t56;
t66 = sin(t75);
t67 = cos(t75);
t156 = g(1) * t66 - g(2) * t67;
t43 = t76 * t121;
t19 = t77 * t36 - t43;
t101 = qJD(4) - t19;
t155 = -pkin(3) - pkin(7);
t14 = t155 * t71 + t101;
t78 = sin(qJ(5));
t81 = cos(qJ(5));
t11 = t81 * qJD(3) + t78 * t14;
t131 = t11 * qJD(5);
t10 = -t78 * qJD(3) + t81 * t14;
t132 = t10 * qJD(5);
t147 = t82 * pkin(1);
t64 = qJDD(1) * t147;
t70 = qJDD(1) + qJDD(2);
t25 = t70 * pkin(2) - t79 * t118 + t64;
t123 = t158 * t76 - t77 * t25;
t107 = qJDD(4) + t123;
t6 = t155 * t70 + t107;
t1 = t81 * qJDD(3) + t78 * t6 + t132;
t151 = t1 * t78;
t3 = t81 * t6;
t2 = -t78 * qJDD(3) - t131 + t3;
t89 = -(t2 + t131) * t81 + t78 * t132 + t140 - t151;
t69 = t71 ^ 2;
t154 = pkin(1) * t79;
t153 = pkin(2) * t66;
t150 = t76 * pkin(2);
t149 = t77 * pkin(2);
t80 = sin(qJ(1));
t148 = t80 * pkin(1);
t133 = qJD(5) * t81;
t13 = t158 * t77 + t76 * t25;
t7 = t70 * qJ(4) + t71 * qJD(4) + t13;
t146 = t17 * t133 + t7 * t78;
t142 = t77 * t79;
t95 = pkin(1) * (t76 * t82 + t142);
t30 = qJD(1) * t95;
t144 = t30 * t71;
t31 = qJD(2) * t95;
t143 = t31 * t71;
t141 = t81 * t70;
t139 = g(1) * t67 + g(2) * t66;
t84 = qJD(5) ^ 2;
t138 = -t69 - t84;
t72 = t78 ^ 2;
t73 = t81 ^ 2;
t137 = t72 - t73;
t136 = t72 + t73;
t32 = t77 * t120 - t43;
t130 = qJD(4) - t32;
t74 = qJDD(3) - g(3);
t53 = t76 * t154;
t63 = pkin(2) + t147;
t34 = t77 * t63 - t53;
t29 = -pkin(3) - t34;
t24 = -pkin(7) + t29;
t128 = qJDD(5) * t24;
t58 = -pkin(3) - t149;
t54 = -pkin(7) + t58;
t127 = qJDD(5) * t54;
t126 = qJDD(5) * t78;
t125 = qJDD(5) * t81;
t124 = t81 * t69 * t78;
t62 = pkin(2) * t67;
t122 = t57 * pkin(3) + t56 * qJ(4) + t62;
t117 = t57 * qJ(4) - t153;
t116 = t136 * t70;
t35 = pkin(1) * t142 + t76 * t63;
t28 = qJ(4) + t35;
t115 = t28 * t71 + t31;
t55 = qJ(4) + t150;
t114 = t55 * t71 - t30;
t83 = cos(qJ(1));
t68 = t83 * pkin(1);
t112 = t68 + t122;
t111 = qJD(1) * (-qJD(2) + t71);
t110 = qJD(2) * (-qJD(1) - t71);
t109 = t71 * t78 * t133;
t106 = t64 + t156;
t105 = t123 - t140;
t104 = -t13 - t103;
t102 = g(1) * t80 - g(2) * t83;
t33 = t77 * t82 * t134 - qJD(2) * t53;
t23 = qJD(4) + t33;
t100 = t17 * t23 + t7 * t28;
t99 = t10 * t81 + t11 * t78;
t98 = -qJD(5) * t14 - t74;
t96 = -t56 * pkin(3) + t117;
t9 = -t70 * pkin(3) + t107;
t94 = t130 * t17 + t7 * t55;
t93 = t155 * t56 + t117;
t92 = t105 - t144;
t91 = t105 + t143;
t90 = -qJD(3) * qJD(5) + t157;
t88 = t151 + t2 * t81 + (-t10 * t78 + t11 * t81) * qJD(5);
t87 = t23 * t71 - t24 * t84 + t28 * t70 + t103;
t86 = t130 * t71 - t54 * t84 + t55 * t70 + t103;
t47 = t57 * pkin(7);
t38 = -t84 * t78 + t125;
t37 = -t84 * t81 - t126;
t27 = t73 * t70 - 0.2e1 * t109;
t26 = t72 * t70 + 0.2e1 * t109;
t18 = 0.2e1 * t137 * t71 * qJD(5) - 0.2e1 * t78 * t141;
t16 = -t71 * pkin(3) + t101;
t5 = t7 * t81;
t4 = [0, 0, 0, 0, 0, qJDD(1), t102, g(1) * t83 + g(2) * t80, 0, 0, 0, 0, 0, 0, 0, t70, (t79 * t110 + t70 * t82) * pkin(1) + t106, ((-qJDD(1) - t70) * t79 + t82 * t110) * pkin(1) + t139, 0, (t102 + (t79 ^ 2 + t82 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t70, t34 * t70 - t91, -t33 * t71 - t35 * t70 + t104, 0, t13 * t35 + t20 * t33 - t123 * t34 - t19 * t31 - g(1) * (-t148 - t153) - g(2) * (t62 + t68), t70, 0, 0, 0, 0, 0, 0, qJDD(4) + (-pkin(3) + t29) * t70 + t91, (qJD(4) + t23) * t71 + (qJ(4) + t28) * t70 - t104, t9 * t29 + t16 * t31 - g(1) * (t96 - t148) - g(2) * t112 + t100, t27, t18, t38, t26, t37, 0, (t115 * qJD(5) + t128) * t81 + t87 * t78 + t146, t5 + (-t128 + (-t115 - t17) * qJD(5)) * t78 + t87 * t81, -t24 * t116 - t136 * t143 + t89, -g(1) * (t93 - t148) - g(2) * (t47 + t112) + t99 * t31 + t88 * t24 + t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t111 * t154 + t106, (t82 * t111 - t129) * pkin(1) + t139, 0, 0, 0, 0, 0, 0, 0, t70, t70 * t149 - t92, -t70 * t150 + t32 * t71 + t104, 0, t19 * t30 - t20 * t32 + (-t123 * t77 + t13 * t76 + t156) * pkin(2), t70, 0, 0, 0, 0, 0, 0, qJDD(4) + (-pkin(3) + t58) * t70 + t92, (0.2e1 * qJD(4) - t32) * t71 + (qJ(4) + t55) * t70 - t104, -g(1) * t96 - g(2) * t122 - t16 * t30 + t9 * t58 + t94, t27, t18, t38, t26, t37, 0, (t114 * qJD(5) + t127) * t81 + t86 * t78 + t146, t5 + (-t127 + (-t114 - t17) * qJD(5)) * t78 + t86 * t81, -t54 * t116 + t136 * t144 + t89, -g(1) * t93 - g(2) * (t47 + t122) - t99 * t30 + t88 * t54 + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, t37, -t38, 0, -t99 * qJD(5) + t1 * t81 - t2 * t78 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t69, t9 + t157, 0, 0, 0, 0, 0, 0, t138 * t78 + t125, t138 * t81 - t126, -t116, -t145 - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, -t137 * t69, t141, -t124, -t78 * t70, qJDD(5), t98 * t78 + t90 * t81 + t131 + t3, t132 + t98 * t81 + (-t6 - t90) * t78, 0, 0;];
tau_reg = t4;
