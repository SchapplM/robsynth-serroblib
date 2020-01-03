% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:53
% EndTime: 2019-12-31 18:17:55
% DurationCPUTime: 0.69s
% Computational Cost: add. (1258->160), mult. (1964->181), div. (0->0), fcn. (1062->12), ass. (0->117)
t64 = qJ(1) + pkin(8);
t59 = qJ(3) + t64;
t52 = sin(t59);
t53 = cos(t59);
t141 = g(1) * t53 + g(2) * t52;
t68 = sin(pkin(8));
t140 = pkin(1) * t68;
t112 = qJD(3) * t140;
t74 = cos(qJ(3));
t122 = qJD(3) * t74;
t69 = cos(pkin(8));
t54 = t69 * pkin(1) + pkin(2);
t71 = sin(qJ(3));
t20 = -t71 * t112 + t54 * t122;
t19 = -qJD(4) - t20;
t28 = t74 * t140 + t71 * t54;
t25 = qJ(4) + t28;
t62 = qJDD(1) + qJDD(3);
t63 = qJD(1) + qJD(3);
t144 = -t19 * t63 + t25 * t62 - t141;
t143 = qJD(1) * t112 - t54 * qJDD(1);
t129 = g(1) * t52 - g(2) * t53;
t124 = t63 * qJ(4);
t113 = qJD(1) * t140;
t39 = t54 * qJD(1);
t132 = t71 * t39;
t18 = t74 * t113 + t132;
t16 = t18 + t124;
t135 = t16 * t63;
t142 = -t129 - t135;
t123 = pkin(1) * qJDD(1);
t76 = -pkin(3) - pkin(7);
t17 = t71 * t113 - t74 * t39;
t90 = qJD(4) + t17;
t12 = t76 * t63 + t90;
t70 = sin(qJ(5));
t73 = cos(qJ(5));
t11 = t73 * qJD(2) + t70 * t12;
t119 = t11 * qJD(5);
t10 = -t70 * qJD(2) + t73 * t12;
t120 = t10 * qJD(5);
t110 = t68 * t123;
t105 = qJD(3) * t132 + t71 * t110 + t143 * t74;
t93 = qJDD(4) + t105;
t6 = t76 * t62 + t93;
t1 = t73 * qJDD(2) + t70 * t6 + t120;
t139 = t1 * t70;
t3 = t73 * t6;
t2 = -t70 * qJDD(2) - t119 + t3;
t83 = -(t2 + t119) * t73 + t70 * t120 + t129 - t139;
t61 = t63 ^ 2;
t138 = t52 * pkin(3);
t137 = t62 * pkin(3);
t121 = qJD(5) * t73;
t55 = t62 * qJ(4);
t56 = t63 * qJD(4);
t86 = -t74 * t110 - t39 * t122 + t143 * t71;
t7 = -t55 - t56 + t86;
t136 = t16 * t121 - t7 * t70;
t134 = t18 * t63;
t21 = t28 * qJD(3);
t133 = t21 * t63;
t131 = t73 * t62;
t130 = t53 * pkin(3) + t52 * qJ(4);
t58 = cos(t64);
t75 = cos(qJ(1));
t128 = t75 * pkin(1) + pkin(2) * t58;
t77 = qJD(5) ^ 2;
t127 = -t61 - t77;
t65 = t70 ^ 2;
t66 = t73 ^ 2;
t126 = t65 - t66;
t125 = t65 + t66;
t67 = qJDD(2) - g(3);
t27 = -t71 * t140 + t74 * t54;
t26 = -pkin(3) - t27;
t22 = -pkin(7) + t26;
t118 = qJDD(5) * t22;
t117 = qJDD(5) * t70;
t116 = qJDD(5) * t73;
t115 = qJDD(5) * t76;
t114 = t73 * t61 * t70;
t111 = t76 * t52;
t109 = t125 * t62;
t108 = t25 * t63 + t21;
t107 = -t18 + t124;
t104 = t128 + t130;
t103 = t63 * t70 * t121;
t57 = sin(t64);
t72 = sin(qJ(1));
t101 = -t72 * pkin(1) - pkin(2) * t57;
t99 = g(1) * t72 - g(2) * t75;
t98 = -t16 * t19 - t7 * t25;
t97 = t10 * t73 + t11 * t70;
t95 = -qJD(5) * t12 - t67;
t92 = -t105 + t129;
t43 = t53 * qJ(4);
t91 = t101 + t43;
t89 = -t7 * qJ(4) + t90 * t16;
t8 = t93 - t137;
t88 = t92 + t134;
t87 = -t92 + t133;
t85 = -qJD(2) * qJD(5) + t142;
t84 = t86 + t141;
t82 = t139 + t2 * t73 + (-t10 * t70 + t11 * t73) * qJD(5);
t81 = -t22 * t77 + t144;
t80 = -t17 * t63 + t84;
t79 = t90 * t63 - t76 * t77 - t141 + t55;
t44 = t53 * pkin(7);
t36 = -t77 * t70 + t116;
t35 = -t77 * t73 - t117;
t24 = t66 * t62 - 0.2e1 * t103;
t23 = t65 * t62 + 0.2e1 * t103;
t15 = -t63 * pkin(3) + t90;
t14 = 0.2e1 * t126 * t63 * qJD(5) - 0.2e1 * t70 * t131;
t5 = t7 * t73;
t4 = [0, 0, 0, 0, 0, qJDD(1), t99, g(1) * t75 + g(2) * t72, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(1) * t57 - g(2) * t58 + 0.2e1 * t69 * t123, g(1) * t58 + g(2) * t57 - 0.2e1 * t110, 0, (t99 + (t68 ^ 2 + t69 ^ 2) * t123) * pkin(1), 0, 0, 0, 0, 0, t62, t27 * t62 - t87, -t20 * t63 - t28 * t62 + t84, 0, -g(1) * t101 - g(2) * t128 - t105 * t27 + t17 * t21 + t18 * t20 - t86 * t28, t62, 0, 0, 0, 0, 0, 0, qJDD(4) + (-pkin(3) + t26) * t62 + t87, -t7 + t144, t8 * t26 + t15 * t21 - g(1) * (t91 - t138) - g(2) * t104 + t98, t24, t14, t36, t23, t35, 0, (t108 * qJD(5) + t118) * t73 + t81 * t70 + t136, -t5 + (-t118 + (-t108 - t16) * qJD(5)) * t70 + t81 * t73, -t22 * t109 - t125 * t133 + t83, -g(1) * (t111 + t91) - g(2) * (t44 + t104) + t97 * t21 + t82 * t22 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, t35, -t36, 0, -t97 * qJD(5) + t1 * t73 - t2 * t70 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t88, t80, 0, 0, t62, 0, 0, 0, 0, 0, 0, qJDD(4) - t88 - 0.2e1 * t137, 0.2e1 * t55 + 0.2e1 * t56 - t80, -t8 * pkin(3) - t15 * t18 - g(1) * (t43 - t138) - g(2) * t130 + t89, t24, t14, t36, t23, t35, 0, (t107 * qJD(5) + t115) * t73 + t79 * t70 + t136, -t5 + (-t115 + (-t107 - t16) * qJD(5)) * t70 + t79 * t73, -t76 * t109 + t125 * t134 + t83, -g(1) * (t43 + t111) - g(2) * (t44 + t130) - t97 * t18 + t82 * t76 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t61, t8 + t142, 0, 0, 0, 0, 0, 0, t127 * t70 + t116, t127 * t73 - t117, -t109, -t135 - t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t126 * t61, t131, -t114, -t70 * t62, qJDD(5), t95 * t70 + t85 * t73 + t119 + t3, t120 + t95 * t73 + (-t6 - t85) * t70, 0, 0;];
tau_reg = t4;
