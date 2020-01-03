% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:20
% EndTime: 2019-12-31 18:57:23
% DurationCPUTime: 0.98s
% Computational Cost: add. (1066->193), mult. (2286->297), div. (0->0), fcn. (1276->4), ass. (0->110)
t89 = 2 * qJD(1);
t59 = sin(qJ(3));
t62 = -pkin(1) - pkin(6);
t119 = t59 * t62;
t61 = cos(qJ(3));
t46 = pkin(3) * t59 - pkin(7) * t61 + qJ(2);
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t109 = t60 * t119 + t58 * t46;
t96 = t60 * qJD(3);
t84 = t59 * t96;
t99 = qJD(4) * t61;
t86 = t58 * t99;
t66 = t84 + t86;
t16 = qJD(1) * t66 - qJD(4) * t96;
t104 = qJD(1) * t61;
t83 = t60 * t104;
t98 = t58 * qJD(3);
t41 = t83 + t98;
t97 = t59 * qJD(1);
t88 = t58 * t97;
t17 = -qJD(3) * t88 + qJD(4) * t41;
t132 = t41 ^ 2;
t53 = qJD(4) + t97;
t28 = t46 * qJD(1);
t52 = qJD(1) * t62 + qJD(2);
t45 = t59 * t52;
t31 = qJD(3) * pkin(7) + t45;
t11 = t60 * t28 - t31 * t58;
t6 = -qJ(5) * t41 + t11;
t5 = pkin(4) * t53 + t6;
t131 = t5 - t6;
t106 = t60 * qJ(5);
t73 = pkin(3) * t61 + pkin(7) * t59;
t44 = t73 * qJD(1);
t30 = t60 * t44;
t112 = -qJ(5) - pkin(7);
t78 = qJD(4) * t112;
t115 = t61 * t52;
t91 = t58 * t115;
t130 = -t58 * qJD(5) + t60 * t78 + t91 - t30 - (pkin(4) * t61 + t106 * t59) * qJD(1);
t129 = t16 * t58;
t128 = t16 * t59;
t127 = t17 * t59;
t32 = -qJD(3) * pkin(3) - t115;
t126 = t32 * t58;
t39 = t104 * t58 - t96;
t125 = t39 * t53;
t124 = t41 * t53;
t123 = t53 * t60;
t122 = t58 * t53;
t121 = t58 * t61;
t120 = t58 * t62;
t118 = t60 * t61;
t117 = t61 * t16;
t116 = t61 * t17;
t63 = qJD(3) ^ 2;
t114 = t63 * t59;
t113 = t63 * t61;
t110 = t60 * t115 + t58 * t44;
t95 = t60 * qJD(5);
t111 = -qJ(5) * t88 + t58 * t78 - t110 + t95;
t57 = t61 ^ 2;
t108 = t59 ^ 2 - t57;
t64 = qJD(1) ^ 2;
t107 = -t63 - t64;
t105 = t64 * qJ(2);
t103 = qJD(3) * t59;
t102 = qJD(3) * t61;
t101 = qJD(4) * t58;
t100 = qJD(4) * t60;
t94 = qJ(2) * qJD(3);
t93 = qJD(1) * qJD(3);
t92 = t58 * t119;
t38 = qJD(3) * t73 + qJD(2);
t87 = t61 * t96;
t90 = t46 * t100 + t58 * t38 + t62 * t87;
t85 = t60 * t99;
t82 = qJD(2) * t89;
t9 = pkin(4) * t17 + t52 * t103;
t81 = pkin(4) - t120;
t80 = pkin(4) * t58 - t62;
t79 = t61 * t93;
t77 = t53 + t97;
t76 = t39 + t96;
t75 = -t41 + t98;
t74 = qJD(4) * t59 + qJD(1);
t12 = t28 * t58 + t31 * t60;
t7 = -qJ(5) * t39 + t12;
t72 = t5 * t60 + t58 * t7;
t71 = t5 * t58 - t60 * t7;
t70 = qJD(1) * t57 - t53 * t59;
t68 = -pkin(7) * t102 + t32 * t59;
t22 = t38 * qJD(1);
t67 = t28 * t100 - t101 * t31 + t22 * t58 + t52 * t87;
t19 = t60 * t22;
t65 = -t12 * qJD(4) + t19;
t49 = t112 * t60;
t48 = t112 * t58;
t37 = t39 ^ 2;
t36 = t60 * t46;
t26 = t60 * t38;
t15 = -qJ(5) * t121 + t109;
t14 = pkin(4) * t39 + qJD(5) + t32;
t13 = -t106 * t61 + t59 * t81 + t36;
t4 = -qJ(5) * t85 + (-qJD(5) * t61 + (qJ(5) * qJD(3) - qJD(4) * t62) * t59) * t58 + t90;
t3 = qJ(5) * t84 + t26 - t109 * qJD(4) + (qJ(5) * t101 + qJD(3) * t81 - t95) * t61;
t2 = -qJ(5) * t17 - qJD(5) * t39 + t67;
t1 = t16 * qJ(5) - t41 * qJD(5) + (pkin(4) * qJD(1) - t52 * t58) * t102 + t65;
t8 = [0, 0, 0, 0, t82, qJ(2) * t82, -0.2e1 * t59 * t79, 0.2e1 * t108 * t93, -t114, -t113, 0, -t62 * t114 + (qJD(2) * t59 + t61 * t94) * t89, -t62 * t113 + (qJD(2) * t61 - t59 * t94) * t89, -t117 * t60 - t41 * t66, (t39 * t60 + t41 * t58) * t103 + (t129 - t17 * t60 + (t39 * t58 - t41 * t60) * qJD(4)) * t61, -t53 * t86 - t128 + (t41 * t61 + t60 * t70) * qJD(3), -t53 * t85 - t127 + (-t39 * t61 - t58 * t70) * qJD(3), t77 * t102, -t62 * t116 + t19 * t59 + t26 * t53 + (-t109 * t53 + t118 * t32 - t12 * t59) * qJD(4) + ((t39 * t62 - t126) * t59 + (-t53 * t120 + (t36 - t92) * qJD(1) + t11) * t61) * qJD(3), -(-qJD(4) * t92 + t90) * t53 - t67 * t59 + (-t101 * t32 + t62 * t16) * t61 + ((-qJD(1) * t109 - t12) * t61 + (t62 * t41 + (-t32 + t115) * t60) * t59) * qJD(3), t13 * t16 - t15 * t17 - t3 * t41 - t4 * t39 + t72 * t103 + (qJD(4) * t71 - t1 * t60 - t2 * t58) * t61, t1 * t13 + t2 * t15 + t5 * t3 + t7 * t4 - t14 * t80 * t103 + (pkin(4) * t100 * t14 + t80 * t9) * t61; 0, 0, 0, 0, -t64, -t105, 0, 0, 0, 0, 0, t107 * t59, t107 * t61, 0, 0, 0, 0, 0, -t116 - t74 * t123 + (-t121 * t77 + t59 * t39) * qJD(3), t117 + t74 * t122 + (-t53 * t118 + (t41 - t83) * t59) * qJD(3), (-t102 * t39 + t41 * t74 - t127) * t60 + (t102 * t41 + t39 * t74 - t128) * t58, -t72 * qJD(1) + (-qJD(3) * t71 - t9) * t61 + (qJD(3) * t14 - qJD(4) * t72 - t1 * t58 + t2 * t60) * t59; 0, 0, 0, 0, 0, 0, t61 * t64 * t59, -t108 * t64, 0, 0, 0, -t61 * t105, t59 * t105, t123 * t41 - t129, (-t16 - t125) * t60 + (-t17 - t124) * t58, t53 * t100 + (t123 * t59 + t61 * t75) * qJD(1), -t53 * t101 + (-t122 * t59 + t61 * t76) * qJD(1), -t53 * t104, -pkin(3) * t17 - t30 * t53 + (t121 * t53 - t59 * t76) * t52 + (-pkin(7) * t123 + t126) * qJD(4) + (-t11 * t61 + t58 * t68) * qJD(1), pkin(3) * t16 + t110 * t53 + t75 * t45 + (pkin(7) * t122 + t32 * t60) * qJD(4) + (t12 * t61 + t60 * t68) * qJD(1), t48 * t16 + t49 * t17 - t130 * t41 - t111 * t39 + (-t5 * t53 + t2) * t60 + (-t53 * t7 - t1) * t58, -t2 * t49 + t1 * t48 + t9 * (-pkin(4) * t60 - pkin(3)) + t111 * t7 + t130 * t5 + (pkin(4) * t122 - t45) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t39, -t37 + t132, t125 - t16, t124 - t17, t79, -qJD(3) * t91 + t12 * t53 - t32 * t41 + t65, t11 * t53 + t32 * t39 - t67, pkin(4) * t16 - t131 * t39, t131 * t7 + (-t14 * t41 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 - t132, t39 * t7 + t41 * t5 + t9;];
tauc_reg = t8;
