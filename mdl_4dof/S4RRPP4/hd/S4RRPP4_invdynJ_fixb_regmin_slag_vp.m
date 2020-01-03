% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPP4
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tau_reg [4x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:19
% EndTime: 2019-12-31 16:59:21
% DurationCPUTime: 0.73s
% Computational Cost: add. (491->186), mult. (1056->218), div. (0->0), fcn. (526->4), ass. (0->119)
t68 = cos(qJ(2));
t120 = qJ(3) * t68;
t132 = pkin(2) + pkin(3);
t66 = sin(qJ(2));
t135 = t132 * t66;
t80 = t120 - t135;
t48 = t66 * qJ(3);
t95 = pkin(1) + t48;
t46 = t66 * qJDD(1);
t105 = qJD(1) * qJD(2);
t96 = t68 * t105;
t138 = t96 + t46;
t108 = qJ(4) * qJD(1);
t117 = qJD(1) * t66;
t44 = pkin(5) * t117;
t17 = t66 * t108 - t44;
t109 = qJD(3) - t17;
t98 = t132 * qJD(2);
t8 = -t98 + t109;
t94 = t132 * qJDD(2);
t121 = pkin(5) * qJD(1);
t45 = t68 * t121;
t19 = -t68 * t108 + t45;
t62 = qJD(2) * qJ(3);
t14 = t19 + t62;
t67 = sin(qJ(1));
t69 = cos(qJ(1));
t137 = g(1) * t69 + g(2) * t67;
t6 = qJD(4) + (t132 * t68 + t95) * qJD(1);
t119 = qJD(4) + t6;
t136 = t119 * t66;
t129 = t66 * t69;
t130 = t66 * t67;
t134 = -g(1) * t129 - g(2) * t130 + g(3) * t68;
t118 = pkin(5) * qJDD(2);
t53 = t68 * pkin(2);
t81 = -t95 - t53;
t15 = t81 * qJD(1);
t125 = t53 + t48;
t21 = -pkin(1) - t125;
t133 = (qJD(1) * t21 + t15) * qJD(2) - t118;
t58 = g(1) * t67;
t131 = g(2) * t69;
t52 = t68 * pkin(3);
t72 = qJD(1) ^ 2;
t128 = t66 * t72;
t127 = t68 * t69;
t126 = pkin(5) - qJ(4);
t63 = t66 ^ 2;
t64 = t68 ^ 2;
t123 = t63 - t64;
t122 = t63 + t64;
t25 = t126 * t68;
t116 = qJD(2) * t25;
t115 = qJD(2) * t66;
t65 = qJDD(1) * pkin(1);
t114 = qJDD(2) * pkin(2);
t113 = t14 * qJD(2);
t112 = t66 * qJD(3);
t111 = t66 * qJD(4);
t110 = t68 * qJD(4);
t107 = qJ(4) * qJD(2);
t47 = t68 * qJDD(1);
t106 = qJ(4) * qJDD(1);
t104 = t68 * t128;
t41 = pkin(5) * t47;
t60 = qJDD(2) * qJ(3);
t61 = qJD(2) * qJD(3);
t103 = t41 + 0.2e1 * t60 + 0.2e1 * t61;
t102 = t41 + t60 + t61;
t101 = t52 + t125;
t32 = pkin(5) * t96;
t40 = pkin(5) * t46;
t100 = qJDD(3) + t32 + t40;
t99 = t58 - t131;
t97 = t66 * t105;
t16 = pkin(1) + t101;
t93 = qJD(1) * t16 + t6;
t92 = pkin(2) * t127 + t67 * pkin(5) + t95 * t69;
t91 = t40 + t134;
t89 = -qJD(2) * pkin(2) + qJD(3);
t88 = 0.2e1 * t96;
t87 = pkin(2) * t47 + t138 * qJ(3) + qJD(1) * t112 + t65;
t71 = qJD(2) ^ 2;
t86 = pkin(5) * t71 + t131;
t84 = pkin(2) * t66 - t120;
t20 = t44 + t89;
t22 = t45 + t62;
t83 = t20 * t68 - t22 * t66;
t82 = -qJDD(3) - t91;
t9 = t100 - t114;
t79 = -0.2e1 * pkin(1) * t105 - t118;
t78 = pkin(3) * t47 + qJDD(4) + t87;
t77 = -t86 + 0.2e1 * t65;
t1 = -t132 * t97 + t78;
t5 = t80 * qJD(2) + t112;
t76 = qJD(1) * t5 + qJDD(1) * t16 + t1 - t131;
t75 = -t66 * t106 + t32 - t82;
t13 = t84 * qJD(2) - t112;
t4 = pkin(2) * t97 - t87;
t74 = -qJD(1) * t13 - qJDD(1) * t21 - t4 - t86;
t7 = -pkin(5) * t97 + t102;
t73 = t83 * qJD(2) + t9 * t66 + t7 * t68;
t54 = t69 * pkin(5);
t36 = t68 * t58;
t35 = g(1) * t130;
t31 = t69 * t120;
t29 = t67 * t120;
t27 = qJ(4) * t97;
t26 = -t63 * t72 - t71;
t24 = t126 * t66;
t23 = -qJDD(2) - t104;
t18 = t84 * qJD(1);
t12 = -t111 + t116;
t11 = -t126 * t115 - t110;
t10 = t80 * qJD(1);
t3 = -t68 * t106 + t27 + (-pkin(5) * t115 - t110) * qJD(1) + t102;
t2 = -t138 * qJ(4) - qJD(1) * t111 + t100 - t94;
t28 = [qJDD(1), t99, t137, t63 * qJDD(1) + t66 * t88, -0.2e1 * t123 * t105 + 0.2e1 * t66 * t47, qJDD(2) * t66 + t71 * t68, qJDD(2) * t68 - t71 * t66, 0, t79 * t66 + t77 * t68 + t36, -t77 * t66 + t79 * t68 - t35, t133 * t66 + t74 * t68 + t36, t122 * qJDD(1) * pkin(5) - t137 + t73, -t133 * t68 + t74 * t66 + t35, t73 * pkin(5) - g(1) * t54 - g(2) * t92 + t15 * t13 + t4 * t21 - t81 * t58, -t24 * qJDD(2) + t36 + (-t93 * t66 - t12) * qJD(2) + t76 * t68, t25 * qJDD(2) + t35 + (t93 * t68 + t11) * qJD(2) + t76 * t66, (-qJD(2) * t8 - qJDD(1) * t25 - t3 + (-qJD(2) * t24 - t11) * qJD(1)) * t68 + (t113 - qJDD(1) * t24 - t2 + (-t12 + t116) * qJD(1)) * t66 + t137, t3 * t25 + t14 * t11 + t2 * t24 + t8 * t12 + t1 * t16 + t6 * t5 - g(1) * (-t69 * qJ(4) + t54) - g(2) * (pkin(3) * t127 + t92) + (-g(1) * (t81 - t52) + g(2) * qJ(4)) * t67; 0, 0, 0, -t104, t123 * t72, t46, t47, qJDD(2), pkin(1) * t128 - t91, g(3) * t66 - t41 + (pkin(1) * t72 + t137) * t68, 0.2e1 * t114 + (-t15 * t66 + t18 * t68) * qJD(1) + t82, -t84 * qJDD(1) + ((t22 - t62) * t66 + (-t20 + t89) * t68) * qJD(1), (qJD(1) * t18 - g(3)) * t66 + (qJD(1) * t15 - t137) * t68 + t103, t7 * qJ(3) + t22 * qJD(3) - t9 * pkin(2) - t15 * t18 - g(1) * (-pkin(2) * t129 + t31) - g(2) * (-pkin(2) * t130 + t29) - g(3) * t125 - t83 * t121, t19 * qJD(2) + 0.2e1 * t94 + ((-t10 + t107) * t68 + t136) * qJD(1) - t75, -t17 * qJD(2) + t27 + (-g(3) + (-pkin(5) * qJD(2) - t10) * qJD(1)) * t66 + (-t119 * qJD(1) - t106 - t137) * t68 + t103, -t80 * qJDD(1), -g(1) * t31 - g(2) * t29 - g(3) * t101 + t3 * qJ(3) - t6 * t10 + t109 * t14 - t132 * t2 + t137 * t135 - t8 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t46, t26, -t22 * qJD(2) + t15 * t117 + t134 + t9, t23, t26, -t46, -t113 - t94 + (-t68 * t107 - t136) * qJD(1) + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 - 0.2e1 * t97, t46 + t88, -t122 * t72, (t14 * t68 + (t8 - t98) * t66) * qJD(1) + t78 + t99;];
tau_reg = t28;
