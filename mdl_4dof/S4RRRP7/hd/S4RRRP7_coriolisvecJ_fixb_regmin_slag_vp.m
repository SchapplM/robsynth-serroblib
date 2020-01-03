% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRP7
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:04
% EndTime: 2019-12-31 17:21:07
% DurationCPUTime: 1.01s
% Computational Cost: add. (956->208), mult. (2444->309), div. (0->0), fcn. (1399->4), ass. (0->107)
t58 = cos(qJ(2));
t94 = t58 * qJD(1);
t46 = -qJD(3) + t94;
t92 = qJD(1) * qJD(2);
t129 = -0.2e1 * t92;
t56 = sin(qJ(2));
t79 = t56 * t92;
t73 = pkin(3) * t79;
t41 = -t58 * pkin(2) - t56 * pkin(6) - pkin(1);
t29 = t41 * qJD(1);
t72 = pkin(2) * t56 - pkin(6) * t58;
t39 = t72 * qJD(2);
t30 = qJD(1) * t39;
t51 = pkin(5) * t94;
t44 = qJD(2) * pkin(6) + t51;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t74 = pkin(5) * t79;
t98 = qJD(3) * t57;
t99 = qJD(3) * t55;
t78 = t29 * t99 - t57 * t30 + t44 * t98 - t55 * t74;
t3 = -t73 + t78;
t12 = t55 * t29 + t57 * t44;
t9 = -t46 * qJ(4) + t12;
t128 = t9 * t46 + t3;
t84 = t56 * t98;
t91 = qJD(2) * qJD(3);
t96 = t55 * qJD(2);
t19 = (t58 * t96 + t84) * qJD(1) + t55 * t91;
t102 = qJD(1) * t57;
t36 = t56 * t102 + t96;
t127 = t36 ^ 2;
t126 = pkin(5) * t57;
t81 = t58 * t92;
t103 = qJD(1) * t56;
t83 = t55 * t103;
t18 = qJD(3) * t83 + (-t81 - t91) * t57;
t2 = t19 * pkin(3) + pkin(5) * t81 + t18 * qJ(4) - t36 * qJD(4);
t125 = t2 * t55;
t124 = t2 * t57;
t95 = t57 * qJD(2);
t34 = t83 - t95;
t43 = -qJD(2) * pkin(2) + pkin(5) * t103;
t10 = t34 * pkin(3) - t36 * qJ(4) + t43;
t122 = t10 * t36;
t121 = t18 * t55;
t120 = t34 * t46;
t119 = t36 * t34;
t118 = t36 * t46;
t117 = t43 * t57;
t116 = t46 * t55;
t115 = t46 * t57;
t114 = t55 * t58;
t38 = t72 * qJD(1);
t113 = t57 * t38;
t112 = t57 * t56;
t111 = t57 * t58;
t60 = qJD(1) ^ 2;
t110 = t58 * t60;
t59 = qJD(2) ^ 2;
t109 = t59 * t56;
t108 = t59 * t58;
t69 = pkin(3) * t55 - qJ(4) * t57;
t107 = t55 * qJD(4) + t46 * t69 + t51;
t106 = t55 * t39 + t41 * t98;
t105 = pkin(5) * t111 + t55 * t41;
t53 = t56 ^ 2;
t104 = -t58 ^ 2 + t53;
t101 = qJD(2) * t56;
t100 = qJD(2) * t58;
t97 = qJD(3) * t58;
t11 = t57 * t29 - t55 * t44;
t93 = qJD(4) - t11;
t90 = pkin(6) * t116;
t89 = pkin(6) * t115;
t88 = pkin(6) * t101;
t87 = pkin(6) * t95;
t86 = t56 * t99;
t85 = t55 * t97;
t82 = pkin(5) * t55 + pkin(3);
t77 = t34 + t95;
t76 = -t36 + t96;
t75 = pkin(1) * t129;
t8 = t46 * pkin(3) + t93;
t71 = -t55 * t9 + t57 * t8;
t70 = t57 * pkin(3) + t55 * qJ(4);
t68 = qJD(1) * t53 - t46 * t58;
t67 = t29 * t98 + t55 * t30 - t44 * t99;
t66 = pkin(5) + t69;
t65 = (qJ(4) - t126) * t103;
t64 = t57 * t39 - t41 * t99;
t63 = -t12 * t46 - t78;
t61 = t56 * t96 - t57 * t97;
t40 = -pkin(2) - t70;
t27 = t55 * t38;
t23 = t66 * t56;
t17 = -t57 * t41 + t82 * t58;
t16 = -t58 * qJ(4) + t105;
t15 = t36 * pkin(3) + t34 * qJ(4);
t14 = -t82 * t103 - t113;
t13 = t27 + t65;
t7 = -t18 - t120;
t6 = (t70 * qJD(3) - qJD(4) * t57) * t56 + t66 * t100;
t5 = -pkin(3) * t101 - t61 * pkin(5) - t64;
t4 = qJ(4) * t101 - t58 * qJD(4) + (-t56 * t95 - t85) * pkin(5) + t106;
t1 = qJD(2) * t65 - t46 * qJD(4) + t67;
t20 = [0, 0, 0, 0.2e1 * t58 * t79, t104 * t129, t108, -t109, 0, -pkin(5) * t108 + t56 * t75, pkin(5) * t109 + t58 * t75, -t18 * t112 + (t58 * t95 - t86) * t36, (-t34 * t57 - t36 * t55) * t100 + (t121 - t19 * t57 + (t34 * t55 - t36 * t57) * qJD(3)) * t56, t46 * t86 + t18 * t58 + (t36 * t56 + t68 * t57) * qJD(2), t46 * t84 + t19 * t58 + (-t34 * t56 - t68 * t55) * qJD(2), (-t46 - t94) * t101, -t64 * t46 + t78 * t58 + t43 * t84 + (t43 * t114 + (t41 * t102 + t11) * t56) * qJD(2) + (t34 * t100 + t56 * t19 - t61 * t46) * pkin(5), (-pkin(5) * t85 + t106) * t46 + t67 * t58 + (-pkin(5) * t18 - t43 * t99) * t56 + ((pkin(5) * t36 + t117) * t58 + (-pkin(5) * t115 - t105 * qJD(1) - t12) * t56) * qJD(2), t23 * t19 + t6 * t34 + t5 * t46 + (t10 * t96 + t3) * t58 + (t10 * t98 + t125 + (-qJD(1) * t17 - t8) * qJD(2)) * t56, -t16 * t19 - t17 * t18 - t4 * t34 + t5 * t36 + t71 * t100 + (-t1 * t55 + t3 * t57 + (-t55 * t8 - t57 * t9) * qJD(3)) * t56, t23 * t18 - t6 * t36 - t4 * t46 + (-t10 * t95 - t1) * t58 + (t10 * t99 - t124 + (qJD(1) * t16 + t9) * qJD(2)) * t56, t1 * t16 + t10 * t6 + t3 * t17 + t2 * t23 + t9 * t4 + t8 * t5; 0, 0, 0, -t56 * t110, t104 * t60, 0, 0, 0, t60 * pkin(1) * t56, pkin(1) * t110, -t36 * t115 - t121, (-t18 + t120) * t57 + (-t19 + t118) * t55, -t46 * t98 + (t46 * t111 + t76 * t56) * qJD(1), t46 * t99 + (-t46 * t114 + t77 * t56) * qJD(1), t46 * t103, t46 * t113 - pkin(2) * t19 + (t43 * t55 + t89) * qJD(3) + (-t11 * t56 + (-t43 * t58 - t88) * t55 + (t56 * t116 - t77 * t58) * pkin(5)) * qJD(1), pkin(2) * t18 - t27 * t46 + (-t90 + t117) * qJD(3) + (-t43 * t111 + (t12 - t87) * t56 + (t46 * t112 + t76 * t58) * pkin(5)) * qJD(1), -t14 * t46 + t40 * t19 - t124 - t107 * t34 + (t10 * t55 + t89) * qJD(3) + (t56 * t8 + (-t10 * t58 - t88) * t55) * qJD(1), t13 * t34 - t14 * t36 + (t1 - t46 * t8 + (qJD(3) * t36 - t19) * pkin(6)) * t57 + ((qJD(3) * t34 - t18) * pkin(6) + t128) * t55, t13 * t46 + t40 * t18 - t125 + t107 * t36 + (-t10 * t57 + t90) * qJD(3) + (t10 * t111 + (-t9 + t87) * t56) * qJD(1), -t9 * t13 - t8 * t14 + t2 * t40 - t107 * t10 + (t71 * qJD(3) + t1 * t57 + t3 * t55) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -t34 ^ 2 + t127, t7, -t118 - t19, t79, -t43 * t36 + t63, -t11 * t46 + t43 * t34 + t57 * t74 - t67, -t15 * t34 - t122 + t63 + 0.2e1 * t73, pkin(3) * t18 - t19 * qJ(4) + (-t12 + t9) * t36 + (t8 - t93) * t34, -t10 * t34 + t15 * t36 + (-0.2e1 * qJD(4) + t11) * t46 + (0.2e1 * qJ(4) - t126) * t79 + t67, -t3 * pkin(3) + t1 * qJ(4) - t10 * t15 - t8 * t12 + t93 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 + t119, t7, -t46 ^ 2 - t127, t122 + t128;];
tauc_reg = t20;
