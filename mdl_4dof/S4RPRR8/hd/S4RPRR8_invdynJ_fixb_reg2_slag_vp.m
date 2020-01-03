% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRR8
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:18
% EndTime: 2019-12-31 16:55:19
% DurationCPUTime: 0.77s
% Computational Cost: add. (1080->175), mult. (2040->225), div. (0->0), fcn. (1213->8), ass. (0->106)
t72 = sin(qJ(1));
t75 = cos(qJ(1));
t125 = g(1) * t72 - g(2) * t75;
t79 = qJD(1) ^ 2;
t136 = -t79 * qJ(2) - t125;
t109 = qJD(1) * qJD(3);
t74 = cos(qJ(3));
t111 = t74 * qJDD(1);
t71 = sin(qJ(3));
t135 = t71 * t109 - t111;
t108 = qJDD(1) * qJ(2);
t67 = t71 ^ 2;
t68 = t74 ^ 2;
t123 = t67 + t68;
t77 = -pkin(1) - pkin(5);
t44 = t77 * qJDD(1) + qJDD(2);
t104 = t123 * t44;
t70 = sin(qJ(4));
t73 = cos(qJ(4));
t33 = t70 * t74 + t73 * t71;
t28 = t33 * qJD(1);
t65 = qJD(3) + qJD(4);
t115 = qJD(4) * t70;
t117 = qJD(3) * t71;
t35 = t74 * t44;
t45 = t77 * qJD(1) + qJD(2);
t15 = qJDD(3) * pkin(3) + t135 * pkin(6) - t45 * t117 + t35;
t116 = qJD(3) * t74;
t102 = t74 * t109;
t112 = t71 * qJDD(1);
t87 = t102 + t112;
t16 = -t87 * pkin(6) + t45 * t116 + t71 * t44;
t118 = qJD(1) * t74;
t26 = -pkin(6) * t118 + t74 * t45;
t23 = qJD(3) * pkin(3) + t26;
t119 = qJD(1) * t71;
t25 = -pkin(6) * t119 + t71 * t45;
t1 = (qJD(4) * t23 + t16) * t73 - t25 * t115 + t70 * t15;
t134 = g(3) * t71;
t133 = t71 * pkin(3);
t132 = pkin(6) - t77;
t106 = t70 * t119;
t30 = t73 * t118 - t106;
t131 = t30 * t28;
t130 = t70 * t25;
t129 = t73 * t25;
t18 = t65 * t33;
t34 = -t70 * t71 + t73 * t74;
t64 = qJDD(3) + qJDD(4);
t128 = -t18 * t65 + t34 * t64;
t110 = qJD(1) * qJD(2);
t105 = 0.2e1 * t110;
t127 = (t105 + t108) * qJ(2);
t126 = t75 * pkin(1) + t72 * qJ(2);
t124 = t67 - t68;
t78 = qJD(3) ^ 2;
t122 = -t78 - t79;
t120 = pkin(1) * qJDD(1);
t51 = qJ(2) + t133;
t39 = t51 * qJD(1);
t114 = t39 * qJD(1);
t113 = qJDD(3) * t71;
t107 = t74 * t79 * t71;
t38 = t132 * t74;
t99 = t123 * qJDD(1);
t98 = t65 * t74;
t97 = qJDD(2) - t120;
t96 = t71 * t102;
t94 = g(1) * t75 + g(2) * t72;
t92 = -t73 * t111 + t70 * t112;
t7 = t18 * qJD(1) + t92;
t91 = -t30 * t18 - t34 * t7;
t19 = -t71 * t115 - t70 * t117 + t73 * t98;
t88 = -qJD(4) * t106 - t135 * t70;
t8 = (qJD(1) * t98 + t112) * t73 + t88;
t90 = t19 * t28 + t33 * t8;
t89 = -t19 * t65 - t33 * t64;
t11 = t70 * t23 + t129;
t37 = t132 * t71;
t21 = -t73 * t37 - t70 * t38;
t20 = t70 * t37 - t73 * t38;
t86 = 0.2e1 * qJ(2) * t109 + qJDD(3) * t77;
t2 = -qJD(4) * t11 + t73 * t15 - t70 * t16;
t84 = t105 - t94 + 0.2e1 * t108;
t10 = t73 * t23 - t130;
t83 = t1 * t33 - t10 * t18 + t11 * t19 + t2 * t34 - t125;
t82 = -t77 * t78 + t84;
t69 = qJ(3) + qJ(4);
t57 = sin(t69);
t58 = cos(t69);
t81 = g(3) * t58 + t125 * t57 + t39 * t28 - t1;
t80 = g(3) * t57 - t125 * t58 - t39 * t30 + t2;
t76 = -pkin(6) - pkin(5);
t60 = t75 * qJ(2);
t56 = qJDD(3) * t74;
t46 = pkin(3) * t116 + qJD(2);
t32 = qJD(3) * t38;
t31 = t132 * t117;
t24 = t87 * pkin(3) + t108 + t110;
t14 = t73 * t26 - t130;
t13 = -t70 * t26 - t129;
t9 = -t28 ^ 2 + t30 ^ 2;
t6 = -t21 * qJD(4) + t73 * t31 + t70 * t32;
t5 = t20 * qJD(4) + t70 * t31 - t73 * t32;
t4 = t30 * t65 + (-t65 * t118 - t112) * t73 - t88;
t3 = [0, 0, 0, 0, 0, qJDD(1), t125, t94, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - 0.2e1 * t120 - t125, t84, -t97 * pkin(1) - g(1) * (-t72 * pkin(1) + t60) - g(2) * t126 + t127, t68 * qJDD(1) - 0.2e1 * t96, 0.2e1 * t124 * t109 - 0.2e1 * t71 * t111, -t78 * t71 + t56, t67 * qJDD(1) + 0.2e1 * t96, -t78 * t74 - t113, 0, t82 * t71 + t86 * t74, -t86 * t71 + t82 * t74, -t77 * t99 - t104 + t125, -g(1) * (t77 * t72 + t60) - g(2) * (t75 * pkin(5) + t126) + t77 * t104 + t127, t91, t18 * t28 - t30 * t19 + t7 * t33 - t34 * t8, t128, t90, t89, 0, t39 * t19 + t20 * t64 + t24 * t33 + t46 * t28 + t51 * t8 - t94 * t57 + t6 * t65, -t39 * t18 - t21 * t64 + t24 * t34 + t46 * t30 - t5 * t65 - t51 * t7 - t94 * t58, t20 * t7 - t21 * t8 - t5 * t28 - t6 * t30 - t83, t1 * t21 + t11 * t5 + t2 * t20 + t10 * t6 + t24 * t51 + t39 * t46 - g(1) * (t75 * t133 + t60 + (-pkin(1) + t76) * t72) - g(2) * (t72 * t133 - t75 * t76 + t126); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t79, t136 + t97, 0, 0, 0, 0, 0, 0, t122 * t71 + t56, t122 * t74 - t113, -t99, t104 + t136, 0, 0, 0, 0, 0, 0, -qJD(1) * t28 + t128, -qJD(1) * t30 + t89, -t90 - t91, t83 - t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, -t124 * t79, t111, -t107, -t112, qJDD(3), t136 * t74 + t134 + t35, g(3) * t74 + (-t44 - t136) * t71, 0, 0, t131, t9, -t92, -t131, t4, t64, -t13 * t65 + (-t65 * t115 - t28 * t118 + t64 * t73) * pkin(3) + t80, t14 * t65 + (-qJD(4) * t65 * t73 - t30 * t118 - t64 * t70) * pkin(3) + t81, (t11 + t13) * t30 + (-t10 + t14) * t28 + (t7 * t73 - t70 * t8 + (-t28 * t73 + t30 * t70) * qJD(4)) * pkin(3), -t10 * t13 - t11 * t14 + (t134 + t1 * t70 + t2 * t73 + (-t10 * t70 + t11 * t73) * qJD(4) + (-t125 - t114) * t74) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t9, -t92, -t131, t4, t64, t11 * t65 + t80, t10 * t65 + t81, 0, 0;];
tau_reg = t3;
