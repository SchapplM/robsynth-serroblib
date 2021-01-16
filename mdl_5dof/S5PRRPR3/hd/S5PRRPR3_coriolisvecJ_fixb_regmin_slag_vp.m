% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:09
% EndTime: 2021-01-15 15:42:14
% DurationCPUTime: 0.78s
% Computational Cost: add. (949->158), mult. (2501->236), div. (0->0), fcn. (1789->6), ass. (0->105)
t94 = sin(qJ(5));
t115 = qJD(5) * t94;
t92 = sin(pkin(9));
t93 = cos(pkin(9));
t95 = sin(qJ(3));
t97 = cos(qJ(3));
t75 = t92 * t97 + t93 * t95;
t118 = qJD(2) * t75;
t125 = t93 * t97;
t110 = qJD(2) * t125;
t116 = qJD(2) * t95;
t68 = t92 * t116 - t110;
t96 = cos(qJ(5));
t59 = t96 * t68;
t70 = t75 * qJD(3);
t62 = qJD(2) * t70;
t114 = qJD(2) * qJD(3);
t107 = t97 * t114;
t108 = t95 * t114;
t82 = t92 * t108;
t63 = t93 * t107 - t82;
t1 = -qJD(5) * t59 - t115 * t118 - t94 * t62 + t96 * t63;
t29 = -t118 * t94 - t59;
t89 = qJD(3) + qJD(5);
t127 = t29 * t89;
t141 = t1 - t127;
t101 = -t118 * t96 + t94 * t68;
t100 = qJD(5) * t101 - t96 * t62 - t94 * t63;
t128 = t101 * t89;
t140 = t100 - t128;
t139 = t101 * t29;
t138 = t101 ^ 2 - t29 ^ 2;
t131 = t68 * pkin(7);
t121 = -qJ(4) - pkin(6);
t81 = t121 * t97;
t66 = t95 * qJD(1) - qJD(2) * t81;
t126 = t93 * t66;
t119 = qJD(3) * pkin(3);
t80 = t121 * t95;
t64 = t97 * qJD(1) + qJD(2) * t80;
t58 = t64 + t119;
t19 = t92 * t58 + t126;
t11 = t19 - t131;
t87 = -t97 * pkin(3) - pkin(2);
t117 = qJD(2) * t87;
t79 = qJD(4) + t117;
t36 = t68 * pkin(4) + t79;
t137 = t11 * t115 - t36 * t29;
t113 = qJD(3) * qJD(1);
t105 = qJD(3) * t121;
t65 = t97 * qJD(4) + t95 * t105;
t38 = qJD(2) * t65 + t97 * t113;
t67 = -t95 * qJD(4) + t105 * t97;
t39 = qJD(2) * t67 - t95 * t113;
t12 = -t92 * t38 + t93 * t39;
t4 = -t63 * pkin(7) + t12;
t13 = t93 * t38 + t92 * t39;
t5 = -t62 * pkin(7) + t13;
t136 = t36 * t101 + t96 * t4 - t94 * t5;
t135 = -0.2e1 * t114;
t134 = qJD(5) - t89;
t133 = pkin(3) * t92;
t132 = pkin(3) * t95;
t74 = t92 * t95 - t125;
t31 = t96 * t74 + t94 * t75;
t73 = t74 * qJD(3);
t7 = -qJD(5) * t31 - t94 * t70 - t96 * t73;
t130 = t7 * t89;
t129 = t118 * pkin(7);
t51 = t92 * t66;
t99 = qJD(2) ^ 2;
t124 = t97 * t99;
t98 = qJD(3) ^ 2;
t123 = t98 * t95;
t122 = t98 * t97;
t22 = t93 * t64 - t51;
t23 = t93 * t65 + t92 * t67;
t41 = t92 * t80 - t93 * t81;
t120 = t95 ^ 2 - t97 ^ 2;
t112 = t95 * t119;
t111 = pkin(3) * t116;
t18 = t93 * t58 - t51;
t20 = -t92 * t64 - t126;
t21 = -t92 * t65 + t93 * t67;
t40 = t93 * t80 + t92 * t81;
t104 = 0.2e1 * t118;
t103 = pkin(2) * t135;
t10 = qJD(3) * pkin(4) - t129 + t18;
t102 = -t94 * t10 - t96 * t11;
t32 = -t94 * t74 + t96 * t75;
t86 = t93 * pkin(3) + pkin(4);
t84 = pkin(3) * t108;
t50 = t74 * pkin(4) + t87;
t43 = t70 * pkin(4) + t112;
t42 = pkin(4) * t118 + t111;
t37 = t62 * pkin(4) + t84;
t25 = -t74 * pkin(7) + t41;
t24 = -t75 * pkin(7) + t40;
t17 = -t70 * pkin(7) + t23;
t16 = t22 - t129;
t15 = t73 * pkin(7) + t21;
t14 = t20 + t131;
t8 = qJD(5) * t32 + t96 * t70 - t94 * t73;
t6 = t8 * t89;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, -t122, -t70 * qJD(3), t73 * qJD(3), t118 * t70 - t75 * t62 + t74 * t63 + t73 * t68, -t12 * t74 + t13 * t75 - t18 * t70 - t19 * t73, 0, 0, 0, 0, 0, -t6, -t130; 0, 0, 0, 0, 0.2e1 * t95 * t107, t120 * t135, t122, -t123, 0, -pkin(6) * t122 + t103 * t95, pkin(6) * t123 + t103 * t97, t87 * t62 + t79 * t70 + (t21 + (qJD(2) * t74 + t68) * t132) * qJD(3), t87 * t63 - t79 * t73 + (t104 * t132 - t23) * qJD(3), -t118 * t21 - t12 * t75 - t13 * t74 + t18 * t73 - t19 * t70 - t23 * t68 - t40 * t63 - t41 * t62, t12 * t40 + t13 * t41 + t18 * t21 + t19 * t23 + (t79 + t117) * t112, t1 * t32 - t101 * t7, -t1 * t31 + t100 * t32 + t101 * t8 + t29 * t7, t130, -t6, 0, -t43 * t29 - t50 * t100 + t37 * t31 + t36 * t8 + (t96 * t15 - t94 * t17 + (-t24 * t94 - t25 * t96) * qJD(5)) * t89, -t43 * t101 + t50 * t1 + t37 * t32 + t36 * t7 - (t94 * t15 + t96 * t17 + (t24 * t96 - t25 * t94) * qJD(5)) * t89; 0, 0, 0, 0, -t95 * t124, t120 * t99, 0, 0, 0, t99 * pkin(2) * t95, pkin(2) * t124, -t20 * qJD(3) - t68 * t111 - t118 * t79 + t12, t22 * qJD(3) - t111 * t118 + t79 * t68 - t13, (t19 + t20) * t118 + (-t18 + t22) * t68 + (-t62 * t92 - t63 * t93) * pkin(3), -t18 * t20 - t19 * t22 + (-t79 * t116 + t12 * t93 + t13 * t92) * pkin(3), t139, t138, t141, t140, 0, t42 * t29 - (t96 * t14 - t94 * t16) * t89 + ((-t96 * t133 - t86 * t94) * t89 + t102) * qJD(5) + t136, -t96 * t5 - t94 * t4 + t42 * t101 + (t94 * t14 + t96 * t16) * t89 + (-(-t94 * t133 + t86 * t96) * t89 - t96 * t10) * qJD(5) + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104 * qJD(3), -t82 + (-t68 + t110) * qJD(3), -t118 ^ 2 - t68 ^ 2, t118 * t18 + t19 * t68 + t84, 0, 0, 0, 0, 0, -t100 - t128, t1 + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t138, t141, t140, 0, t134 * t102 + t136, (-t11 * t89 - t4) * t94 + (-t134 * t10 - t5) * t96 + t137;];
tauc_reg = t2;
