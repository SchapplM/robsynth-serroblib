% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPR2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:37
% EndTime: 2019-12-05 16:17:41
% DurationCPUTime: 0.85s
% Computational Cost: add. (852->175), mult. (1173->245), div. (0->0), fcn. (709->10), ass. (0->112)
t56 = pkin(8) + qJ(2);
t51 = qJ(3) + t56;
t43 = cos(t51);
t129 = g(2) * t43;
t53 = qJDD(2) + qJDD(3);
t128 = t53 * pkin(3);
t65 = cos(qJ(3));
t127 = t65 * pkin(2);
t103 = pkin(2) * qJD(2);
t63 = sin(qJ(3));
t89 = t63 * t103;
t106 = -qJD(3) * t89 + qJDD(2) * t127;
t83 = qJDD(4) - t106;
t24 = t83 - t128;
t132 = t24 + t129;
t101 = t53 * qJ(4);
t57 = qJD(2) + qJD(3);
t93 = qJDD(2) * t63;
t98 = qJD(3) * t65;
t15 = t101 + t57 * qJD(4) + (qJD(2) * t98 + t93) * pkin(2);
t60 = sin(pkin(9));
t61 = cos(pkin(9));
t12 = t60 * qJDD(1) + t61 * t15;
t11 = -t61 * qJDD(1) + t60 * t15;
t126 = t11 * t60;
t74 = t12 * t61 + t126;
t42 = sin(t51);
t39 = g(1) * t42;
t131 = t106 + t39;
t100 = qJD(2) * t65;
t75 = -pkin(2) * t100 + qJD(4);
t28 = t57 * qJ(4) + t89;
t22 = -t61 * qJD(1) + t60 * t28;
t125 = t22 * t60;
t23 = t60 * qJD(1) + t61 * t28;
t117 = t61 * t57;
t33 = -qJD(5) + t117;
t95 = qJD(5) + t33;
t130 = -t57 * t125 - t23 * t95;
t118 = t61 * t53;
t32 = -qJDD(5) + t118;
t124 = t32 * t61;
t36 = pkin(2) * t98 + qJD(4);
t123 = t36 * t57;
t62 = sin(qJ(5));
t122 = t36 * t62;
t121 = t53 * t62;
t52 = t57 ^ 2;
t54 = t60 ^ 2;
t120 = t54 * t52;
t119 = t60 * t62;
t116 = t61 * t62;
t64 = cos(qJ(5));
t115 = t61 * t64;
t114 = t61 * t65;
t113 = t62 * t32;
t112 = t62 * t64;
t111 = t64 * t32;
t110 = t64 * t33;
t109 = t132 * t60;
t108 = t43 * pkin(3) + t42 * qJ(4);
t107 = g(1) * t43 + g(2) * t42;
t105 = t61 ^ 2 + t54;
t59 = t64 ^ 2;
t104 = t62 ^ 2 - t59;
t102 = qJ(4) * t62;
t99 = qJD(3) * t63;
t97 = qJD(5) * t62;
t96 = qJD(5) * t64;
t94 = qJ(4) * qJD(5);
t91 = t33 * t103;
t90 = pkin(2) * t99;
t87 = t57 * t99;
t86 = t57 * t96;
t85 = t33 * t97;
t84 = t64 * t94;
t82 = -t42 * pkin(3) + t43 * qJ(4);
t80 = t105 * t53;
t79 = t32 - t118;
t78 = t32 + t118;
t77 = t57 * t95;
t76 = t57 * t89;
t73 = t23 * t61 + t125;
t72 = -t107 + t74;
t29 = -t61 * pkin(4) - t60 * pkin(7) - pkin(3);
t10 = t29 * t53 + t83;
t13 = t29 * t57 + t75;
t18 = t42 * t116 + t43 * t64;
t20 = -t43 * t116 + t42 * t64;
t71 = -g(1) * t18 - g(2) * t20 + (t62 * t10 + t64 * t12 + (t64 * t13 - t62 * t23) * qJD(5)) * t61 + t64 * t126;
t19 = -t42 * t115 + t43 * t62;
t21 = t43 * t115 + t42 * t62;
t70 = -g(1) * t19 - g(2) * t21 + t11 * t119 + t96 * t125;
t69 = t76 - t129;
t45 = -pkin(3) - t127;
t68 = pkin(2) * t87 + t45 * t53;
t67 = g(3) * t60 - t13 * t95 - t12;
t66 = -t33 ^ 2 - t120;
t50 = cos(t56);
t49 = sin(t56);
t44 = t63 * pkin(2) + qJ(4);
t31 = t61 * t39;
t27 = -t57 * pkin(3) + t75;
t26 = t60 * t97 * t117;
t25 = t29 - t127;
t16 = (t53 * t59 - 0.2e1 * t62 * t86) * t54;
t8 = 0.2e1 * (qJD(5) * t104 * t57 - t53 * t112) * t54;
t7 = t64 * t10;
t4 = (t78 * t62 + (t33 + t117) * t96) * t60;
t3 = t26 + (-t64 * t78 + t85) * t60;
t2 = -t62 * t12 + t7 + (-t62 * t13 - t64 * t23) * qJD(5);
t1 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t61 + t12 * t60 - g(3), 0, 0, 0, 0, 0, (t79 * t62 + (t33 - t117) * t96) * t60, t26 + (t64 * t79 - t85) * t60; 0, qJDD(2), g(1) * t49 - g(2) * t50, g(1) * t50 + g(2) * t49, t53, -t129 + (t53 * t65 - t87) * pkin(2) + t131, ((-qJDD(2) - t53) * t63 + (-qJD(2) - t57) * t98) * pkin(2) + t107, t31 + (-t68 - t132) * t61, (t68 - t39) * t60 + t109, t105 * t123 + t44 * t80 + t72, t24 * t45 + t27 * t90 - g(1) * (-pkin(2) * t49 + t82) - g(2) * (pkin(2) * t50 + t108) + t74 * t44 + t73 * t36, t16, t8, t3, t4, t124, -(-t25 * t97 + t64 * t90) * t33 - t25 * t111 + (-(-t44 * t96 - t122) * t33 + t44 * t113 - t2) * t61 + (t57 * t122 + (t86 + t121) * t44) * t54 + t70, (t36 * t115 + t62 * t90) * t33 + (t44 * t115 + t62 * t25) * t32 + (t44 * t53 + t123) * t64 * t54 + (t25 * t110 + (-t125 + (-t33 * t61 - t54 * t57) * t44) * t62) * qJD(5) + t71; 0, 0, 0, 0, t53, t69 + t131, (-t93 + (-qJD(3) + t57) * t100) * pkin(2) + t107, t31 + (-t24 + t69 + t128) * t61, (-t128 - t76 - t39) * t60 + t109, t75 * t57 * t105 + qJ(4) * t80 + t72, -t24 * pkin(3) - g(1) * t82 - g(2) * t108 + t73 * qJD(4) + t74 * qJ(4) + (-t27 * t63 - t65 * t73) * t103, t16, t8, t3, t4, t124, (t85 - t111) * t29 + (-(-qJD(4) * t62 - t84) * t33 + t32 * t102 - t2) * t61 + (-t62 * t114 + t63 * t64) * t91 + (t62 * t101 + (t62 * t75 + t84) * t57) * t54 + t70, t61 * qJD(4) * t110 + (qJ(4) * t115 + t62 * t29) * t32 + ((-t102 * t61 + t64 * t29) * t33 - t22 * t119) * qJD(5) - (t64 * t114 + t62 * t63) * t91 + (t64 * t101 + (-t62 * t94 + t64 * t75) * t57) * t54 + t71; 0, 0, 0, 0, 0, 0, 0, -t118, t60 * t53, -t105 * t52, -t57 * t73 + t132 - t39, 0, 0, 0, 0, 0, t62 * t66 - t111, t64 * t66 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112 * t120, -t104 * t120, (t53 * t64 - t62 * t77) * t60, (-t64 * t77 - t121) * t60, -t32, -g(1) * t20 + g(2) * t18 + t130 * t64 + t67 * t62 + t7, g(1) * t21 - g(2) * t19 + t67 * t64 + (-t10 - t130) * t62;];
tau_reg = t1;
