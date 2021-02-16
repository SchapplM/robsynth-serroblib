% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:14:40
% EndTime: 2021-01-15 16:14:44
% DurationCPUTime: 0.75s
% Computational Cost: add. (950->199), mult. (1330->233), div. (0->0), fcn. (695->8), ass. (0->126)
t142 = qJDD(1) - g(3);
t74 = cos(qJ(4));
t68 = qJD(2) + qJD(3);
t121 = qJ(5) * t68;
t122 = pkin(2) * qJD(2);
t73 = sin(qJ(3));
t105 = t73 * t122;
t29 = t68 * pkin(7) + t105;
t97 = t29 + t121;
t87 = t97 * t74;
t67 = pkin(8) + qJ(2);
t60 = qJ(3) + t67;
t51 = sin(t60);
t46 = g(2) * t51;
t52 = cos(t60);
t48 = g(1) * t52;
t126 = t46 + t48;
t72 = sin(qJ(4));
t113 = t72 * qJD(4);
t102 = t68 * t113;
t54 = t74 * pkin(4) + pkin(3);
t66 = qJDD(2) + qJDD(3);
t131 = t54 * t66;
t141 = -pkin(4) * t102 + t131;
t69 = t72 ^ 2;
t140 = pkin(4) * t69;
t47 = g(1) * t51;
t139 = g(2) * t52;
t138 = g(3) * t74;
t137 = t66 * pkin(3);
t75 = cos(qJ(3));
t136 = t75 * pkin(2);
t62 = t74 * qJD(1);
t10 = -t97 * t72 + t62;
t120 = qJD(4) * pkin(4);
t6 = t10 + t120;
t135 = -t10 + t6;
t118 = qJD(2) * t75;
t104 = pkin(2) * t118;
t30 = -t68 * pkin(3) - t104;
t134 = t30 * t68;
t133 = t51 * t74;
t132 = t52 * t72;
t65 = t68 ^ 2;
t130 = t65 * t74;
t129 = t68 * t72;
t128 = t68 * t74;
t49 = t72 * t66;
t127 = t74 * t66;
t71 = qJ(5) + pkin(7);
t125 = -qJD(3) * t105 + qJDD(2) * t136;
t70 = t74 ^ 2;
t124 = -t69 - t70;
t123 = t69 - t70;
t53 = t73 * pkin(2) + pkin(7);
t119 = -qJ(5) - t53;
t117 = qJD(3) * t73;
t116 = qJD(3) * t75;
t115 = qJD(4) * t68;
t114 = qJDD(4) * pkin(4);
t112 = t74 * qJD(4);
t111 = qJDD(2) * t73;
t23 = t29 * t113;
t17 = t66 * pkin(7) + (qJD(2) * t116 + t111) * pkin(2);
t92 = -qJD(4) * qJD(1) - t17;
t83 = -qJ(5) * t66 + t92;
t79 = qJD(5) * t68 - t83;
t3 = -t23 + (-qJ(5) * t115 + qJDD(1)) * t72 + t79 * t74;
t110 = t3 * t74 - t126;
t15 = -t54 * t68 + qJD(5) - t104;
t35 = g(2) * t132;
t5 = qJDD(5) - t125 - t141;
t109 = t15 * t112 + t5 * t72 + t35;
t16 = -t125 - t137;
t108 = t30 * t112 + t16 * t72 + t35;
t36 = g(1) * t133;
t91 = qJD(4) * t104;
t93 = t68 * t105;
t107 = t72 * t91 + t74 * t93 + t36;
t106 = pkin(2) * t116;
t103 = t68 * t117;
t12 = t15 * t113;
t101 = -t5 - t139;
t100 = -t16 - t139;
t99 = t51 * t71 + t52 * t54;
t98 = -t51 * t54 + t71 * t52;
t96 = qJD(4) * t71;
t95 = 0.2e1 * t68 * t112;
t94 = qJD(4) * t119;
t76 = qJD(4) ^ 2;
t90 = -pkin(7) * t76 + t137;
t11 = t72 * qJD(1) + t87;
t89 = t11 * t74 - t6 * t72;
t27 = pkin(2) * t117 + pkin(4) * t113;
t39 = -t54 - t136;
t88 = t27 * t68 + t39 * t66;
t86 = -t125 - t47 + t139;
t59 = t74 * qJDD(1);
t85 = g(1) * t132 + t72 * t46 - t138 + t59;
t84 = -pkin(3) * t115 - pkin(7) * qJDD(4);
t82 = g(2) * t133 - t142 * t72 + t74 * t48 + t23;
t81 = -t93 - t47;
t55 = -pkin(3) - t136;
t80 = pkin(2) * t103 + t53 * t76 + t55 * t66;
t78 = (-qJD(5) - t15) * t68 + t83;
t77 = -qJDD(4) * t53 + (t55 * t68 - t106) * qJD(4);
t63 = t74 * qJ(5);
t61 = t74 * qJD(5);
t58 = cos(t67);
t57 = sin(t67);
t44 = t74 * pkin(7) + t63;
t43 = t71 * t72;
t41 = t74 * t91;
t32 = qJDD(4) * t74 - t76 * t72;
t31 = qJDD(4) * t72 + t76 * t74;
t25 = t74 * t53 + t63;
t24 = t119 * t72;
t21 = t30 * t113;
t20 = -t72 * qJD(5) - t74 * t96;
t19 = -t72 * t96 + t61;
t18 = t69 * t66 + t72 * t95;
t9 = -0.2e1 * t123 * t115 + 0.2e1 * t72 * t127;
t8 = (-qJD(5) - t106) * t72 + t74 * t94;
t7 = t74 * t106 + t72 * t94 + t61;
t2 = -qJD(4) * t87 - t79 * t72 + t114 + t59;
t1 = [t142, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, t32, -t31, 0, t89 * qJD(4) + t2 * t74 + t3 * t72 - g(3); 0, qJDD(2), g(1) * t57 - g(2) * t58, g(1) * t58 + g(2) * t57, t66, (t66 * t75 - t103) * pkin(2) - t86, ((-qJDD(2) - t66) * t73 + (-qJD(2) - t68) * t116) * pkin(2) + t126, t18, t9, t31, t32, 0, t21 + t36 + t77 * t72 + (t100 - t80) * t74, t77 * t74 + (t80 - t47) * t72 + t108, t24 * qJDD(4) + t12 + t36 + (t39 * t129 + t8) * qJD(4) + (t101 - t88) * t74, -t25 * qJDD(4) + (t39 * t128 - t7) * qJD(4) + (t88 - t47) * t72 + t109, (t25 * t66 + t68 * t7 + (-t24 * t68 - t6) * qJD(4)) * t74 + (-t24 * t66 - t68 * t8 - t2 + (-t25 * t68 - t11) * qJD(4)) * t72 + t110, t3 * t25 + t11 * t7 + t2 * t24 + t6 * t8 + t5 * t39 + t15 * t27 - g(1) * (-pkin(2) * t57 + t98) - g(2) * (pkin(2) * t58 + t99); 0, 0, 0, 0, t66, -t86 + t93, (-t111 + (-qJD(3) + t68) * t118) * pkin(2) + t126, t18, t9, t31, t32, 0, t21 + t84 * t72 + (t100 + t90) * t74 + t107, t41 + t84 * t74 + (t81 - t90) * t72 + t108, -t43 * qJDD(4) + t12 + (-t54 * t129 + t20) * qJD(4) + (t101 + t141) * t74 + t107, -t44 * qJDD(4) + t41 + (t81 - t131) * t72 + (-t19 + (-t54 * t74 + t140) * t68) * qJD(4) + t109, (-qJD(4) * t6 + t44 * t66) * t74 + (-t11 * qJD(4) + t43 * t66 - t2) * t72 + (t19 * t74 - t20 * t72 + (t43 * t74 - t44 * t72) * qJD(4) + t124 * t104) * t68 + t110, t3 * t44 + t11 * t19 - t2 * t43 + t6 * t20 - t5 * t54 + pkin(4) * t12 - g(1) * t98 - g(2) * t99 + (-t15 * t73 - t89 * t75) * t122; 0, 0, 0, 0, 0, 0, 0, -t72 * t130, t123 * t65, t49, t127, qJDD(4), (-t17 - t134) * t72 + t85, (-t72 * t29 + t62) * qJD(4) + (t92 - t134) * t74 + t82, 0.2e1 * t114 + (t11 - t87) * qJD(4) + (pkin(4) * t130 + t78) * t72 + t85, -t65 * t140 + (t72 * t121 + t10) * qJD(4) + t78 * t74 + t82, -pkin(4) * t49 + (-t120 + t135) * t128, t135 * t11 + (-t138 + t2 + (-t15 * t68 + t126) * t72) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t102 - t127, t49 + t95, t124 * t65, -t89 * t68 - t101 - t47;];
tau_reg = t1;
