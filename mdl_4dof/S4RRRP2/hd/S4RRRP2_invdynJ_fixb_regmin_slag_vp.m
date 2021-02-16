% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRP2
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:29
% EndTime: 2021-01-15 11:04:31
% DurationCPUTime: 0.74s
% Computational Cost: add. (785->176), mult. (1169->219), div. (0->0), fcn. (597->8), ass. (0->118)
t135 = qJ(4) + pkin(6);
t65 = sin(qJ(3));
t60 = qJD(1) + qJD(2);
t116 = pkin(1) * qJD(1);
t66 = sin(qJ(2));
t99 = t66 * t116;
t90 = t135 * t60 + t99;
t12 = t90 * t65;
t68 = cos(qJ(3));
t13 = t90 * t68;
t63 = qJ(1) + qJ(2);
t54 = sin(t63);
t46 = g(2) * t54;
t55 = cos(t63);
t48 = g(1) * t55;
t120 = t46 + t48;
t50 = t68 * pkin(3) + pkin(2);
t59 = qJDD(1) + qJDD(2);
t127 = t50 * t59;
t109 = qJD(3) * t65;
t96 = t60 * t109;
t134 = -pkin(3) * t96 + t127;
t61 = t65 ^ 2;
t133 = pkin(3) * t61;
t47 = g(1) * t54;
t132 = g(2) * t55;
t131 = g(3) * t68;
t130 = t59 * pkin(2);
t69 = cos(qJ(2));
t129 = t69 * pkin(1);
t115 = qJD(3) * pkin(3);
t9 = -t12 + t115;
t128 = t12 + t9;
t126 = t54 * t68;
t125 = t55 * t65;
t58 = t60 ^ 2;
t124 = t58 * t68;
t123 = t60 * t65;
t122 = t60 * t68;
t44 = t65 * t59;
t121 = t68 * t59;
t119 = -qJD(2) * t99 + qJDD(1) * t129;
t62 = t68 ^ 2;
t118 = -t61 - t62;
t117 = t61 - t62;
t49 = t66 * pkin(1) + pkin(6);
t114 = -qJ(4) - t49;
t113 = qJD(1) * t69;
t112 = qJD(2) * t66;
t111 = qJD(2) * t69;
t110 = qJD(3) * t60;
t108 = qJD(3) * t68;
t107 = qJDD(3) * pkin(3);
t106 = qJDD(1) * t66;
t17 = t59 * pkin(6) + (qJD(1) * t111 + t106) * pkin(1);
t91 = -qJ(4) * t59 - t17;
t76 = qJD(4) * t60 - t91;
t78 = qJD(3) * t90;
t3 = -t65 * t78 + t68 * t76;
t105 = t3 * t68 - t120;
t98 = pkin(1) * t113;
t15 = -t50 * t60 + qJD(4) - t98;
t40 = g(2) * t125;
t5 = qJDD(4) - t119 - t134;
t104 = t15 * t108 + t5 * t65 + t40;
t16 = -t119 - t130;
t28 = -t60 * pkin(2) - t98;
t103 = t28 * t108 + t16 * t65 + t40;
t41 = g(1) * t126;
t85 = qJD(3) * t98;
t86 = t60 * t99;
t102 = t65 * t85 + t68 * t86 + t41;
t101 = g(2) * t126 + g(3) * t65 + t68 * t48;
t100 = pkin(1) * t111;
t97 = t60 * t112;
t10 = t15 * t109;
t95 = -t5 - t132;
t94 = -t16 - t132;
t93 = -t28 * t60 - t17;
t92 = t135 * t54 + t55 * t50;
t89 = qJD(3) * t135;
t88 = 0.2e1 * t60 * t108;
t87 = qJD(3) * t114;
t84 = g(1) * t125 + t65 * t46 - t131;
t71 = qJD(3) ^ 2;
t83 = -pkin(6) * t71 + t130;
t82 = -t13 * t68 + t65 * t9;
t25 = pkin(1) * t112 + pkin(3) * t109;
t32 = -t50 - t129;
t81 = t25 * t60 + t32 * t59;
t80 = t135 * t55 - t54 * t50;
t79 = -t119 - t47 + t132;
t77 = -pkin(2) * t110 - pkin(6) * qJDD(3);
t75 = -t86 - t47;
t74 = (-qJD(4) - t15) * t60 + t91;
t51 = -pkin(2) - t129;
t73 = pkin(1) * t97 + t49 * t71 + t51 * t59;
t72 = -qJDD(3) * t49 + (t51 * t60 - t100) * qJD(3);
t70 = cos(qJ(1));
t67 = sin(qJ(1));
t56 = t68 * qJ(4);
t53 = t68 * qJD(4);
t36 = t68 * pkin(6) + t56;
t35 = t135 * t65;
t34 = t68 * t85;
t31 = qJDD(3) * t68 - t71 * t65;
t30 = qJDD(3) * t65 + t71 * t68;
t24 = t68 * t49 + t56;
t23 = t114 * t65;
t21 = t28 * t109;
t20 = -t65 * qJD(4) - t68 * t89;
t19 = -t65 * t89 + t53;
t18 = t61 * t59 + t65 * t88;
t8 = -0.2e1 * t110 * t117 + 0.2e1 * t65 * t121;
t7 = (-qJD(4) - t100) * t65 + t68 * t87;
t6 = t100 * t68 + t65 * t87 + t53;
t2 = -t65 * t76 - t68 * t78 + t107;
t1 = [qJDD(1), g(1) * t67 - g(2) * t70, g(1) * t70 + g(2) * t67, t59, (t59 * t69 - t97) * pkin(1) - t79, ((-qJDD(1) - t59) * t66 + (-qJD(1) - t60) * t111) * pkin(1) + t120, t18, t8, t30, t31, 0, t21 + t41 + t72 * t65 + (-t73 + t94) * t68, t72 * t68 + (t73 - t47) * t65 + t103, t23 * qJDD(3) + t10 + t41 + (t32 * t123 + t7) * qJD(3) + (-t81 + t95) * t68, -t24 * qJDD(3) + (t32 * t122 - t6) * qJD(3) + (t81 - t47) * t65 + t104, (t24 * t59 + t6 * t60 + (-t23 * t60 - t9) * qJD(3)) * t68 + (-t23 * t59 - t60 * t7 - t2 + (-t24 * t60 - t13) * qJD(3)) * t65 + t105, t3 * t24 + t13 * t6 + t2 * t23 + t9 * t7 + t5 * t32 + t15 * t25 - g(1) * (-t67 * pkin(1) + t80) - g(2) * (t70 * pkin(1) + t92); 0, 0, 0, t59, -t79 + t86, (-t106 + (-qJD(2) + t60) * t113) * pkin(1) + t120, t18, t8, t30, t31, 0, t21 + t77 * t65 + (t83 + t94) * t68 + t102, t34 + t77 * t68 + (t75 - t83) * t65 + t103, -t35 * qJDD(3) + t10 + (-t50 * t123 + t20) * qJD(3) + (t95 + t134) * t68 + t102, -t36 * qJDD(3) + t34 + (t75 - t127) * t65 + (-t19 + (-t50 * t68 + t133) * t60) * qJD(3) + t104, (-qJD(3) * t9 + t36 * t59) * t68 + (-t13 * qJD(3) + t35 * t59 - t2) * t65 + (t19 * t68 - t20 * t65 + (t35 * t68 - t36 * t65) * qJD(3) + t118 * t98) * t60 + t105, t3 * t36 + t13 * t19 - t2 * t35 + t9 * t20 - t5 * t50 + pkin(3) * t10 - g(1) * t80 - g(2) * t92 + (-t15 * t66 + t69 * t82) * t116; 0, 0, 0, 0, 0, 0, -t65 * t124, t117 * t58, t44, t121, qJDD(3), t65 * t93 + t84, t68 * t93 + t101, 0.2e1 * t107 + (pkin(3) * t124 + t74) * t65 + t84, -t58 * t133 + t74 * t68 + t101, -pkin(3) * t44 + (-t115 + t128) * t122, t128 * t13 + (-t131 + t2 + (-t15 * t60 + t120) * t65) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t96 - t121, t44 + t88, t118 * t58, t60 * t82 - t47 - t95;];
tau_reg = t1;
