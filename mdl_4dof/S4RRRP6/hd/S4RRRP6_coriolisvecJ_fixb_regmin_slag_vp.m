% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauc_reg [4x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:15
% EndTime: 2021-01-15 14:39:20
% DurationCPUTime: 0.96s
% Computational Cost: add. (965->208), mult. (2504->311), div. (0->0), fcn. (1471->4), ass. (0->110)
t92 = (qJD(1) * qJD(2));
t136 = -2 * t92;
t65 = sin(qJ(2));
t101 = qJD(1) * t65;
t64 = sin(qJ(3));
t88 = t64 * t101;
t66 = cos(qJ(3));
t95 = t66 * qJD(2);
t37 = t88 - t95;
t67 = cos(qJ(2));
t93 = t67 * qJD(1);
t55 = -qJD(3) + t93;
t124 = t37 * t55;
t87 = t67 * t92;
t91 = qJD(2) * qJD(3);
t16 = qJD(3) * t88 + (-t87 - t91) * t66;
t135 = -t16 + t124;
t96 = t64 * qJD(2);
t39 = t66 * t101 + t96;
t123 = t39 * t55;
t98 = qJD(3) * t66;
t89 = t65 * t98;
t72 = t67 * t96 + t89;
t17 = t72 * qJD(1) + t64 * t91;
t134 = t17 - t123;
t133 = t39 ^ 2;
t44 = -t67 * pkin(2) - t65 * pkin(6) - pkin(1);
t29 = t44 * qJD(1);
t60 = pkin(5) * t93;
t50 = qJD(2) * pkin(6) + t60;
t12 = t66 * t29 - t64 * t50;
t6 = -t39 * qJ(4) + t12;
t5 = -t55 * pkin(3) + t6;
t132 = t5 - t6;
t131 = pkin(3) * t64;
t130 = t37 * pkin(3);
t118 = t64 * t29;
t13 = t66 * t50 + t118;
t7 = -t37 * qJ(4) + t13;
t129 = t7 * t55;
t9 = t17 * pkin(3) + pkin(5) * t87;
t128 = t9 * t64;
t127 = t9 * t66;
t76 = pkin(2) * t65 - pkin(6) * t67;
t41 = t76 * qJD(1);
t106 = pkin(5) * t88 + t66 * t41;
t102 = t66 * qJ(4);
t74 = pkin(3) * t65 - t67 * t102;
t111 = qJ(4) + pkin(6);
t83 = qJD(3) * t111;
t126 = -t74 * qJD(1) - t64 * qJD(4) - t66 * t83 - t106;
t125 = t16 * t64;
t122 = t39 * t64;
t49 = -qJD(2) * pkin(2) + pkin(5) * t101;
t121 = t49 * t64;
t120 = t49 * t66;
t119 = t55 * t66;
t117 = t64 * t67;
t116 = t65 * t66;
t115 = t66 * t67;
t69 = qJD(1) ^ 2;
t114 = t67 * t69;
t68 = qJD(2) ^ 2;
t113 = t68 * t65;
t112 = t68 * t67;
t103 = t64 * qJ(4);
t25 = t64 * t41;
t94 = t66 * qJD(4);
t110 = t25 + (-pkin(5) * t116 - t67 * t103) * qJD(1) + t64 * t83 - t94;
t42 = t76 * qJD(2);
t30 = qJD(1) * t42;
t85 = t65 * t92;
t78 = pkin(5) * t85;
t109 = -t66 * t30 - t64 * t78;
t108 = t64 * t42 + t44 * t98;
t107 = t65 * pkin(5) * t96 + t66 * t42;
t56 = pkin(5) * t115;
t105 = t64 * t44 + t56;
t62 = t65 ^ 2;
t104 = -t67 ^ 2 + t62;
t100 = qJD(2) * t67;
t99 = qJD(3) * t64;
t97 = t49 * qJD(3);
t90 = t65 * t99;
t84 = -qJD(4) - t130;
t82 = -t29 * t98 - t64 * t30 + t50 * t99 + t66 * t78;
t81 = t37 + t95;
t80 = -t39 + t96;
t79 = pkin(1) * t136;
t77 = pkin(3) * t85;
t75 = qJD(1) * t62 - t55 * t67;
t73 = t17 * qJ(4) + t82;
t71 = -t13 * qJD(3) - t109;
t70 = t16 * qJ(4) + t71;
t58 = -t66 * pkin(3) - pkin(2);
t48 = t111 * t66;
t47 = t111 * t64;
t43 = (pkin(5) + t131) * t65;
t36 = t66 * t44;
t34 = t37 ^ 2;
t32 = t93 * t131 + t60;
t18 = t72 * pkin(3) + pkin(5) * t100;
t15 = t49 - t84;
t14 = -t65 * t103 + t105;
t11 = -t65 * t102 + t36 + (-pkin(5) * t64 - pkin(3)) * t67;
t4 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t116 + (-qJD(4) * t65 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t67) * t64 + t108;
t3 = -t65 * t94 + t74 * qJD(2) + (-t56 + (qJ(4) * t65 - t44) * t64) * qJD(3) + t107;
t2 = -t37 * qJD(4) - t73;
t1 = -t39 * qJD(4) + t70 + t77;
t8 = [0, 0, 0, 0.2e1 * t67 * t85, t104 * t136, t112, -t113, 0, -pkin(5) * t112 + t65 * t79, pkin(5) * t113 + t67 * t79, -t16 * t116 + (t67 * t95 - t90) * t39, (-t37 * t66 - t122) * t100 + (t125 - t17 * t66 + (t37 * t64 - t39 * t66) * qJD(3)) * t65, t55 * t90 + t16 * t67 + (t39 * t65 + t75 * t66) * qJD(2), t55 * t89 + t17 * t67 + (-t37 * t65 - t75 * t64) * qJD(2), (-t55 - t93) * t65 * qJD(2), -(-t44 * t99 + t107) * t55 + (t66 * t97 + pkin(5) * t17 + (qJD(1) * t36 + t12) * qJD(2)) * t65 + ((pkin(5) * t37 + t121) * qJD(2) + (t118 + (pkin(5) * t55 + t50) * t66) * qJD(3) + t109) * t67, (-t67 * pkin(5) * t99 + t108) * t55 - t82 * t67 + (-pkin(5) * t16 - t64 * t97) * t65 + ((pkin(5) * t39 + t120) * t67 + (-t105 * qJD(1) - t13 + (-t55 + t93) * t66 * pkin(5)) * t65) * qJD(2), t43 * t17 + t18 * t37 - t3 * t55 + (t15 * t96 - t1) * t67 + (t15 * t98 + t128 + (qJD(1) * t11 + t5) * qJD(2)) * t65, -t43 * t16 + t18 * t39 + t4 * t55 + (t15 * t95 + t2) * t67 + (-t15 * t99 + t127 + (-qJD(1) * t14 - t7) * qJD(2)) * t65, t11 * t16 - t14 * t17 - t3 * t39 - t4 * t37 + (-t5 * t66 - t64 * t7) * t100 + (-t1 * t66 - t2 * t64 + (t5 * t64 - t66 * t7) * qJD(3)) * t65, t1 * t11 + t2 * t14 + t15 * t18 + t5 * t3 + t7 * t4 + t9 * t43; 0, 0, 0, -t65 * t114, t104 * t69, 0, 0, 0, t69 * pkin(1) * t65, pkin(1) * t114, -t39 * t119 - t125, -t134 * t64 + t135 * t66, -t55 * t98 + (t55 * t115 + t80 * t65) * qJD(1), t55 * t99 + (-t55 * t117 + t81 * t65) * qJD(1), t55 * t101, -pkin(2) * t17 + t106 * t55 + (pkin(6) * t119 + t121) * qJD(3) + ((-pkin(6) * t96 - t12) * t65 + (-t81 * pkin(5) - t121) * t67) * qJD(1), pkin(2) * t16 - t25 * t55 + (-t64 * pkin(6) * t55 + t120) * qJD(3) + (-t49 * t115 + (-pkin(6) * t95 + t13) * t65 + (t55 * t116 + t80 * t67) * pkin(5)) * qJD(1), t58 * t17 - t32 * t37 - t127 - t126 * t55 + (t15 + t130) * t99 + (-t15 * t117 + (-qJD(2) * t47 - t5) * t65) * qJD(1), -t58 * t16 - t32 * t39 + t128 - t110 * t55 + (pkin(3) * t122 + t15 * t66) * qJD(3) + (-t15 * t115 + (-qJD(2) * t48 + t7) * t65) * qJD(1), -t47 * t16 - t48 * t17 - t126 * t39 + t110 * t37 + (t55 * t5 + t2) * t66 + (-t1 + t129) * t64, -t1 * t47 + t2 * t48 + t9 * t58 - t110 * t7 + t126 * t5 + (pkin(3) * t99 - t32) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t37, -t34 + t133, -t16 - t124, -t123 - t17, t85, -t13 * t55 - t49 * t39 + t71, -t12 * t55 + t49 * t37 + t82, 0.2e1 * t77 - t129 + (-t15 + t84) * t39 + t70, -t133 * pkin(3) - t6 * t55 + (qJD(4) + t15) * t37 + t73, t16 * pkin(3) - t132 * t37, t132 * t7 + (-t15 * t39 + t1) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t135, -t34 - t133, t7 * t37 + t5 * t39 + t9;];
tauc_reg = t8;
