% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRR5
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:49
% EndTime: 2019-12-31 16:33:50
% DurationCPUTime: 0.62s
% Computational Cost: add. (722->134), mult. (1372->192), div. (0->0), fcn. (961->10), ass. (0->97)
t59 = sin(qJ(4));
t54 = t59 ^ 2;
t62 = cos(qJ(4));
t55 = t62 ^ 2;
t107 = t54 + t55;
t52 = qJDD(2) + qJDD(3);
t63 = cos(qJ(3));
t105 = qJD(3) * t63;
t101 = qJD(1) * qJD(2);
t64 = cos(qJ(2));
t48 = t64 * qJDD(1);
t61 = sin(qJ(2));
t23 = qJDD(2) * pkin(2) - t61 * t101 + t48;
t40 = qJD(2) * pkin(2) + t64 * qJD(1);
t60 = sin(qJ(3));
t77 = -t61 * qJDD(1) - t64 * t101;
t106 = qJD(1) * t61;
t89 = qJD(3) * t106;
t85 = -t40 * t105 + t77 * t63 + (-t23 + t89) * t60;
t4 = t52 * pkin(6) - t85;
t94 = t107 * t4;
t56 = qJ(2) + qJ(3);
t49 = sin(t56);
t58 = cos(pkin(7));
t113 = t49 * t58;
t57 = sin(pkin(7));
t114 = t49 * t57;
t125 = g(1) * t113 + g(2) * t114;
t50 = cos(t56);
t119 = g(3) * t50;
t118 = t52 * pkin(3);
t20 = t63 * t23;
t7 = -(qJD(3) * t40 - t77) * t60 - t63 * t89 + t20;
t5 = -t118 - t7;
t124 = t5 + t119;
t25 = t60 * t61 - t63 * t64;
t53 = qJD(2) + qJD(3);
t123 = t25 * t53;
t22 = t25 * qJD(1);
t122 = pkin(2) * t105 + t22;
t121 = -t119 + t125;
t116 = t60 * pkin(2);
t46 = pkin(6) + t116;
t47 = -t63 * pkin(2) - pkin(3);
t65 = qJD(4) ^ 2;
t26 = t60 * t64 + t63 * t61;
t21 = t26 * qJD(1);
t86 = qJD(3) * t116 - t21;
t120 = t46 * t65 + t47 * t52 + t86 * t53;
t44 = g(3) * t49;
t117 = t53 * pkin(3);
t16 = t63 * t106 + t60 * t40;
t115 = t16 * t53;
t112 = t50 * t57;
t111 = t50 * t58;
t110 = t62 * t52;
t109 = t50 * pkin(3) + t49 * pkin(6);
t108 = t54 - t55;
t104 = qJD(4) * t62;
t103 = qJDD(1) - g(3);
t15 = -t60 * t106 + t63 * t40;
t13 = -t15 - t117;
t100 = t13 * t104 + t124 * t59;
t51 = t53 ^ 2;
t99 = t59 * t51 * t62;
t98 = t13 * qJD(4) * t59 + t125 * t62;
t97 = -g(1) * t111 - g(2) * t112 - t44;
t93 = t107 * t123;
t92 = t107 * t15;
t91 = t107 * t52;
t88 = t59 * t53 * t104;
t87 = t97 + t94;
t84 = -pkin(2) * t61 - pkin(3) * t49;
t83 = g(1) * t58 + g(2) * t57;
t82 = g(1) * t57 - g(2) * t58;
t9 = t53 * t26;
t81 = -t25 * t52 - t9 * t53;
t78 = t83 * t49;
t76 = pkin(6) * t65 - t115 - t118;
t75 = t26 * t65 - t81;
t74 = t85 - t97;
t73 = -pkin(6) * qJDD(4) + (t15 - t117) * qJD(4);
t72 = 0.2e1 * t123 * qJD(4) - qJDD(4) * t26;
t71 = -g(3) * t64 + t83 * t61;
t69 = t122 * t107;
t68 = -qJDD(4) * t46 + (t47 * t53 - t122) * qJD(4);
t67 = -t13 * t53 + t83 * t50 - t4 + t44;
t66 = qJD(2) ^ 2;
t32 = pkin(6) * t111;
t31 = pkin(6) * t112;
t30 = qJDD(4) * t62 - t65 * t59;
t29 = qJDD(4) * t59 + t65 * t62;
t18 = t55 * t52 - 0.2e1 * t88;
t17 = t54 * t52 + 0.2e1 * t88;
t14 = t53 * pkin(6) + t16;
t12 = -0.2e1 * t108 * t53 * qJD(4) + 0.2e1 * t59 * t110;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, 0, 0, 0, 0, 0, t64 * qJDD(2) - t66 * t61, -qJDD(2) * t61 - t66 * t64, 0, -g(3) + (t61 ^ 2 + t64 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t81, t123 * t53 - t26 * t52, 0, -t123 * t16 - t15 * t9 - t7 * t25 - t26 * t85 - g(3), 0, 0, 0, 0, 0, 0, t72 * t59 - t75 * t62, t75 * t59 + t72 * t62, t26 * t91 - t53 * t93, t13 * t9 - t14 * t93 + t5 * t25 + t26 * t94 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t48 + t71, -t103 * t61 + t83 * t64, 0, 0, 0, 0, 0, 0, 0, t52, t21 * t53 + t20 + (pkin(2) * t52 - t89) * t63 + ((-pkin(2) * t53 - t40) * qJD(3) + t77) * t60 + t121, -t22 * t53 + (-t53 * t105 - t52 * t60) * pkin(2) + t74, 0, t15 * t21 + t16 * t22 + (-t85 * t60 + t63 * t7 + (-t15 * t60 + t16 * t63) * qJD(3) + t71) * pkin(2), t17, t12, t29, t18, t30, 0, t68 * t59 + (-t120 - t124) * t62 + t98, t68 * t62 + (t120 - t78) * t59 + t100, t46 * t91 + t69 * t53 + t87, t5 * t47 - g(1) * (t84 * t58 + t32) - g(2) * (t84 * t57 + t31) - g(3) * (t64 * pkin(2) + t109) + t46 * t94 + t86 * t13 + t69 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t115 + t7 + t121, t15 * t53 + t74, 0, 0, t17, t12, t29, t18, t30, 0, t73 * t59 + (-t76 - t124) * t62 + t98, t73 * t62 + (-t78 + t76) * t59 + t100, pkin(6) * t91 - t53 * t92 + t87, -t5 * pkin(3) - t13 * t16 - g(1) * (-pkin(3) * t113 + t32) - g(2) * (-pkin(3) * t114 + t31) - g(3) * t109 - t14 * t92 + pkin(6) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t108 * t51, t59 * t52, t99, t110, qJDD(4), t67 * t59 - t82 * t62, t82 * t59 + t67 * t62, 0, 0;];
tau_reg = t1;
