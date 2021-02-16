% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:56:38
% EndTime: 2021-01-15 10:56:43
% DurationCPUTime: 1.05s
% Computational Cost: add. (1200->220), mult. (2908->311), div. (0->0), fcn. (2051->10), ass. (0->118)
t84 = sin(qJ(2));
t124 = qJD(1) * t84;
t125 = cos(pkin(7));
t87 = cos(qJ(2));
t112 = t125 * t87;
t66 = qJD(1) * t112;
t81 = sin(pkin(7));
t48 = t124 * t81 - t66;
t42 = qJD(4) + t48;
t82 = -qJ(3) - pkin(5);
t116 = t82 * t84;
t109 = qJD(2) * t82;
t96 = -t84 * qJD(3) + t109 * t87;
t24 = qJDD(2) * pkin(2) + qJD(1) * t96 + qJDD(1) * t116;
t47 = t87 * qJD(3) + t109 * t84;
t63 = t82 * t87;
t31 = qJD(1) * t47 - qJDD(1) * t63;
t7 = t125 * t24 - t81 * t31;
t5 = -qJDD(2) * pkin(3) - t7;
t113 = t125 * t84;
t58 = t81 * t87 + t113;
t51 = t58 * qJD(1);
t69 = pkin(2) * t81 + pkin(6);
t85 = sin(qJ(1));
t88 = cos(qJ(1));
t105 = g(1) * t88 + g(2) * t85;
t78 = qJ(2) + pkin(7);
t73 = sin(t78);
t74 = cos(t78);
t95 = -g(3) * t74 + t105 * t73;
t147 = -(pkin(2) * t124 + pkin(3) * t51 + pkin(6) * t48 + qJD(4) * t69) * t42 - t5 + t95;
t118 = qJD(1) * qJD(2);
t115 = t84 * t118;
t91 = qJDD(1) * t58 - t115 * t81;
t30 = qJD(2) * t66 + t91;
t83 = sin(qJ(4));
t86 = cos(qJ(4));
t38 = qJD(2) * t83 + t51 * t86;
t10 = qJD(4) * t38 - qJDD(2) * t86 + t83 * t30;
t72 = pkin(2) * t87 + pkin(1);
t62 = -qJD(1) * t72 + qJD(3);
t11 = t48 * pkin(3) - t51 * pkin(6) + t62;
t8 = t125 * t31 + t24 * t81;
t111 = qJDD(2) * pkin(6) + qJD(4) * t11 + t8;
t13 = t125 * t47 + t81 * t96;
t61 = qJD(1) * t63;
t133 = t81 * t61;
t126 = qJD(2) * pkin(2);
t60 = qJD(1) * t116;
t56 = t60 + t126;
t27 = t125 * t56 + t133;
t19 = -qJD(2) * pkin(3) - t27;
t120 = t84 * qJDD(1);
t103 = -qJDD(1) * t112 + t120 * t81;
t50 = t58 * qJD(2);
t29 = qJD(1) * t50 + t103;
t25 = qJDD(4) + t29;
t99 = -t81 * t84 + t112;
t26 = -pkin(3) * t99 - pkin(6) * t58 - t72;
t35 = t116 * t81 - t125 * t63;
t53 = t99 * qJD(2);
t145 = -(qJD(4) * t26 + t13) * t42 + t111 * t99 + t19 * t53 - t35 * t25 + t5 * t58;
t144 = pkin(2) * t84;
t143 = g(3) * t73;
t141 = g(3) * t87;
t121 = t86 * qJD(2);
t123 = qJD(4) * t83;
t9 = qJD(4) * t121 + qJDD(2) * t83 - t123 * t51 + t30 * t86;
t140 = t9 * t83;
t139 = t19 * t58;
t138 = t26 * t25;
t36 = t51 * t83 - t121;
t137 = t36 * t42;
t136 = t38 * t42;
t135 = t38 * t51;
t134 = t51 * t36;
t132 = t83 * t25;
t131 = t85 * t83;
t130 = t85 * t86;
t18 = t86 * t25;
t129 = t88 * t83;
t128 = t88 * t86;
t54 = t125 * t61;
t28 = t56 * t81 - t54;
t79 = t84 ^ 2;
t127 = -t87 ^ 2 + t79;
t122 = qJD(4) * t86;
t119 = t87 * qJDD(1);
t117 = t84 * t126;
t20 = qJD(2) * pkin(6) + t28;
t41 = pkin(2) * t115 - qJDD(1) * t72 + qJDD(3);
t4 = t29 * pkin(3) - t30 * pkin(6) + t41;
t110 = qJD(4) * t20 - t4;
t108 = t42 * t86;
t104 = g(1) * t85 - g(2) * t88;
t102 = t18 + (-t48 * t83 - t123) * t42;
t101 = -t111 + t143;
t100 = -t123 * t58 + t53 * t86;
t98 = -0.2e1 * pkin(1) * t118 - pkin(5) * qJDD(2);
t33 = t125 * t60 + t133;
t94 = -t69 * t25 + (t19 + t33) * t42;
t89 = qJD(2) ^ 2;
t93 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t89 + t104;
t90 = qJD(1) ^ 2;
t92 = pkin(1) * t90 - pkin(5) * qJDD(1) + t105;
t70 = -pkin(2) * t125 - pkin(3);
t46 = t128 * t74 + t131;
t45 = -t129 * t74 + t130;
t44 = -t130 * t74 + t129;
t43 = t131 * t74 + t128;
t34 = -t113 * t82 - t63 * t81;
t32 = t60 * t81 - t54;
t15 = pkin(3) * t50 - pkin(6) * t53 + t117;
t12 = -t125 * t96 + t47 * t81;
t3 = t11 * t83 + t20 * t86;
t2 = t11 * t86 - t20 * t83;
t1 = t86 * t4;
t6 = [qJDD(1), t104, t105, qJDD(1) * t79 + 0.2e1 * t115 * t87, -0.2e1 * t118 * t127 + 0.2e1 * t119 * t84, qJDD(2) * t84 + t87 * t89, qJDD(2) * t87 - t84 * t89, 0, t84 * t98 + t87 * t93, -t84 * t93 + t87 * t98, -t34 * qJDD(2) - t72 * t29 - t41 * t99 + t62 * t50 + t104 * t74 + (t144 * t48 - t12) * qJD(2), -t35 * qJDD(2) - t72 * t30 + t41 * t58 + t62 * t53 - t104 * t73 + (t144 * t51 - t13) * qJD(2), t12 * t51 - t13 * t48 - t27 * t53 - t28 * t50 - t29 * t35 + t30 * t34 - t58 * t7 + t8 * t99 - t105, t8 * t35 + t28 * t13 - t7 * t34 - t27 * t12 - t41 * t72 + t62 * t117 - g(1) * (-t72 * t85 - t82 * t88) - g(2) * (t72 * t88 - t82 * t85), t9 * t86 * t58 + t100 * t38, (-t36 * t86 - t38 * t83) * t53 + (-t10 * t86 - t140 + (t36 * t83 - t38 * t86) * qJD(4)) * t58, t100 * t42 + t18 * t58 + t38 * t50 - t9 * t99, -t58 * t132 + t10 * t99 - t36 * t50 + (-t122 * t58 - t53 * t83) * t42, -t25 * t99 + t42 * t50, -g(1) * t44 - g(2) * t46 - t1 * t99 + t34 * t10 + t12 * t36 + t2 * t50 + (t15 * t42 + t138 + (t20 * t99 - t35 * t42 + t139) * qJD(4)) * t86 + t145 * t83, -g(1) * t43 - g(2) * t45 + t12 * t38 - t3 * t50 + t34 * t9 + (-(-qJD(4) * t35 + t15) * t42 - t138 - t110 * t99 - qJD(4) * t139) * t83 + t145 * t86; 0, 0, 0, -t84 * t90 * t87, t127 * t90, t120, t119, qJDD(2), t84 * t92 - t141, g(3) * t84 + t87 * t92, t32 * qJD(2) - t62 * t51 + (qJDD(2) * t125 - t124 * t48) * pkin(2) + t95 + t7, t143 + t33 * qJD(2) + t62 * t48 + t105 * t74 + (-qJDD(2) * t81 - t124 * t51) * pkin(2) - t8, (t28 - t32) * t51 + (-t27 + t33) * t48 + (-t125 * t30 - t29 * t81) * pkin(2), t27 * t32 - t28 * t33 + (t125 * t7 - t141 + t8 * t81 + (-qJD(1) * t62 + t105) * t84) * pkin(2), t108 * t38 + t140, (t9 - t137) * t86 + (-t10 - t136) * t83, t108 * t42 + t132 - t135, t102 + t134, -t42 * t51, t70 * t10 + t147 * t86 - t2 * t51 - t32 * t36 + t94 * t83, -t147 * t83 + t3 * t51 - t32 * t38 + t70 * t9 + t94 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t51 * qJD(2) + t103, (t66 - t48) * qJD(2) + t91, -t48 ^ 2 - t51 ^ 2, t27 * t51 + t28 * t48 - t104 + t41, 0, 0, 0, 0, 0, t102 - t134, -t42 ^ 2 * t86 - t132 - t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t36, -t36 ^ 2 + t38 ^ 2, t9 + t137, -t10 + t136, t25, -g(1) * t45 + g(2) * t43 + t101 * t83 - t122 * t20 - t19 * t38 + t3 * t42 + t1, g(1) * t46 - g(2) * t44 + t101 * t86 + t110 * t83 + t19 * t36 + t2 * t42;];
tau_reg = t6;
