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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
% StartTime: 2019-03-08 18:26:28
% EndTime: 2019-03-08 18:26:28
% DurationCPUTime: 0.60s
% Computational Cost: add. (440->169), mult. (1321->235), div. (0->0), fcn. (981->6), ass. (0->137)
t76 = sin(pkin(4));
t72 = t76 ^ 2;
t77 = cos(pkin(6));
t156 = t72 * t77;
t78 = cos(pkin(4));
t155 = pkin(1) * t78;
t79 = sin(qJ(1));
t154 = g(1) * t79;
t75 = sin(pkin(6));
t71 = t75 ^ 2;
t153 = t71 * t72;
t152 = t75 * t76;
t151 = t76 * t77;
t80 = cos(qJ(1));
t150 = t76 * t80;
t81 = qJD(1) ^ 2;
t149 = t78 * t81;
t148 = t79 * t75;
t147 = t79 * t77;
t146 = t80 * t75;
t145 = t80 * t77;
t132 = qJDD(1) * t76;
t116 = qJ(2) * t132;
t130 = qJD(1) * qJD(2);
t114 = t76 * t130;
t46 = t75 * t114;
t144 = t75 * t116 + t46;
t138 = qJD(1) * t76;
t121 = qJ(2) * t138;
t137 = qJD(1) * t78;
t123 = t75 * t137;
t25 = pkin(1) * t123 + t77 * t121;
t142 = qJ(2) * t76;
t36 = t77 * t142 + t75 * t155;
t143 = t80 * pkin(1) + t79 * t142;
t141 = t78 * qJ(3);
t94 = pkin(3) * t151 + t141;
t8 = t94 * qJD(1) + qJD(4) + t25;
t140 = -qJD(4) - t8;
t136 = qJD(2) * t76;
t39 = t78 * qJD(3) + t77 * t136;
t139 = qJD(1) * t39;
t49 = t75 * t121;
t135 = qJD(3) + t49;
t134 = qJDD(1) * t72;
t133 = qJDD(1) * t75;
t131 = t78 * qJDD(1);
t129 = qJD(1) * qJD(3);
t128 = t77 * t155;
t127 = t75 * t156;
t126 = t77 * t149;
t125 = pkin(1) * t131;
t48 = t77 * t114;
t12 = t77 * t116 + t75 * t125 + t48;
t124 = t75 * t138;
t122 = -pkin(1) * t77 - pkin(2);
t120 = t75 * t132;
t119 = t77 * t132;
t118 = qJDD(3) + t144;
t117 = -t79 * pkin(1) + t80 * t142;
t115 = t72 * t130;
t113 = t75 * t129;
t112 = -qJ(3) * t75 - pkin(1);
t111 = t81 * t127;
t104 = t122 * t78;
t60 = t75 * t142;
t19 = t60 + t104;
t98 = t122 * qJDD(1);
t6 = t78 * t98 + t118;
t110 = qJDD(1) * t19 + t6;
t95 = -pkin(2) * t77 + t112;
t20 = t95 * t76;
t86 = t95 * qJDD(1);
t7 = qJDD(2) + (t86 - t113) * t76;
t109 = qJDD(1) * t20 + t7;
t108 = qJDD(1) * t127;
t107 = t78 * t120;
t106 = t78 * t119;
t105 = -qJ(4) + t122;
t4 = -qJ(3) * t131 - t78 * t129 - t12;
t31 = -t78 * t145 + t148;
t33 = t78 * t147 + t146;
t103 = g(1) * t31 - g(2) * t33;
t32 = t78 * t146 + t147;
t34 = -t78 * t148 + t145;
t102 = g(1) * t32 - g(2) * t34;
t101 = g(1) * t80 + g(2) * t79;
t100 = g(2) * t150 - g(3) * t78 + qJDD(2);
t99 = (qJD(1) * t128 - t49) * t75 - t25 * t77;
t97 = t34 * pkin(2) + t33 * qJ(3) + t143;
t96 = -qJD(3) * t75 - qJD(4) * t77;
t93 = t105 * qJDD(1);
t57 = -pkin(1) * t132 + qJDD(2);
t92 = pkin(1) * t134 - t57 * t76;
t53 = pkin(3) * t120;
t1 = t53 + (-qJD(1) * qJD(4) + t93) * t78 + t118;
t82 = pkin(3) * t152 + t105 * t78;
t10 = t60 + t82;
t37 = -t78 * qJD(4) + t75 * t136;
t91 = qJD(1) * t37 + qJDD(1) * t10 + t1;
t13 = t94 + t36;
t3 = pkin(3) * t119 + qJDD(4) - t4;
t90 = qJDD(1) * t13 + t139 + t3;
t18 = -t141 - t36;
t89 = -qJDD(1) * t18 + t139 - t4;
t88 = -t32 * pkin(2) - t31 * qJ(3) + t117;
t87 = (-pkin(2) - qJ(4)) * t77 + t112;
t14 = t87 * t76;
t83 = t87 * qJDD(1);
t2 = qJDD(2) + (t96 * qJD(1) + t83) * t76;
t30 = t96 * t76;
t85 = (-qJD(1) * t30 - qJDD(1) * t14 - t2) * t76;
t84 = -g(1) * t33 - g(2) * t31 + g(3) * t151 + t118;
t74 = t78 ^ 2;
t73 = t77 ^ 2;
t69 = t74 * qJDD(1);
t59 = t73 * t134;
t58 = t71 * t134;
t44 = t71 * t115;
t43 = t149 * t152;
t42 = -0.2e1 * t106;
t41 = 0.2e1 * t107;
t40 = 0.2e1 * t108;
t38 = (-t74 - t153) * t81;
t35 = -t60 + t128;
t29 = (-t71 - t73) * t81 * t72;
t28 = t111 + t131;
t23 = -t43 + t119;
t22 = (-t126 + t133) * t76;
t21 = (t126 + t133) * t76;
t17 = -qJ(3) * t137 - t25;
t16 = qJD(1) * t20 + qJD(2);
t15 = qJD(1) * t104 + t135;
t11 = t77 * t125 - t144;
t9 = qJD(1) * t14 + qJD(2);
t5 = t82 * qJD(1) + t135;
t24 = [0, 0, 0, 0, 0, qJDD(1), -g(2) * t80 + t154, t101, 0, 0, t58, t40, t41, t59, 0.2e1 * t106, t69, t92 * t77 + (qJDD(1) * t35 + t11 - t46) * t78 + t102, -t92 * t75 + (-qJDD(1) * t36 - t12 - t48) * t78 - t103, t73 * t115 + t44 + (-t11 * t75 + t12 * t77 + (-t35 * t75 + t36 * t77) * qJDD(1) - t101) * t76, t12 * t36 + t11 * t35 - g(1) * t117 - g(2) * t143 + (-t57 * pkin(1) - t99 * qJD(2)) * t76, t69, -0.2e1 * t107, t42, t58, t40, t59, t44 + (t110 * t75 + t89 * t77 - t101) * t76, -t113 * t156 + t110 * t78 + (qJD(2) * t123 + t109 * t77) * t76 - t102, -t109 * t152 + t129 * t153 + t89 * t78 + t103, t7 * t20 + t4 * t18 - t17 * t39 + t6 * t19 - g(1) * t88 - g(2) * t97 + (qJD(2) * t15 - qJD(3) * t16) * t152, t69, t42, t41, t59, -0.2e1 * t108, t58 (t91 * t75 + t90 * t77 - t101) * t76, t75 * t85 + t90 * t78 + t103, t77 * t85 - t91 * t78 + t102, t2 * t14 + t9 * t30 + t1 * t10 + t5 * t37 + t3 * t13 + t8 * t39 - g(1) * (pkin(3) * t150 - t32 * qJ(4) + t88) - g(2) * (t79 * t76 * pkin(3) + t34 * qJ(4) + t97); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t21, t29 (-pkin(1) * qJDD(1) + t99 * qJD(1) - t154) * t76 + t100, 0, 0, 0, 0, 0, 0, t29, t23, -t21 (-t154 + t86 + (t17 * t77 + (-qJD(3) - t15) * t75) * qJD(1)) * t76 + t100, 0, 0, 0, 0, 0, 0, t29, -t21, -t23 (-t154 + t83 + (t140 * t77 + (-qJD(3) - t5) * t75) * qJD(1)) * t76 + t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t28, t38, t16 * t124 + (qJD(1) * t17 + t98) * t78 + t84, 0, 0, 0, 0, 0, 0, t22, t38, -t28, t9 * t124 + t53 + (t140 * qJD(1) + t93) * t78 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 + t119, -t111 + t131 (-t72 * t73 - t74) * t81, -g(3) * t152 - g(1) * t34 - g(2) * t32 + (t9 * t151 + t5 * t78) * qJD(1) + t3;];
tau_reg  = t24;
