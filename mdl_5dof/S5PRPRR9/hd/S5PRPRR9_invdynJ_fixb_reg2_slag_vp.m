% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRR9
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:49
% EndTime: 2019-12-31 17:39:50
% DurationCPUTime: 0.66s
% Computational Cost: add. (1162->148), mult. (1609->170), div. (0->0), fcn. (835->6), ass. (0->94)
t54 = qJDD(2) - qJDD(4);
t62 = cos(qJ(4));
t111 = t62 * t54;
t98 = qJD(2) - qJD(4);
t125 = t98 ^ 2;
t60 = sin(qJ(4));
t64 = qJD(5) ^ 2;
t126 = (t64 + t125) * t60 + t111;
t102 = qJ(3) * qJD(2);
t63 = -pkin(2) - pkin(3);
t35 = t63 * qJD(2) + qJD(3);
t17 = t62 * t102 + t60 * t35;
t114 = t17 * t98;
t118 = t98 * pkin(4);
t16 = -t60 * t102 + t62 * t35;
t10 = -t16 + t118;
t124 = t98 * t10;
t115 = t98 * t16;
t123 = -qJD(4) * t102 + t63 * qJDD(2) + qJDD(3);
t103 = qJD(5) * t98;
t101 = qJ(3) * qJDD(2);
t99 = qJD(2) * qJD(3);
t121 = qJD(4) * t35 + t101 + t99;
t120 = t62 * t125;
t28 = t62 * qJ(3) + t60 * t63;
t11 = -pkin(7) * t98 + t17;
t59 = sin(qJ(5));
t61 = cos(qJ(5));
t7 = -t61 * qJD(1) - t59 * t11;
t106 = t7 * qJD(5);
t5 = t121 * t62 + t123 * t60;
t3 = -t54 * pkin(7) + t5;
t1 = -t59 * qJDD(1) + t61 * t3 + t106;
t100 = qJD(1) * qJD(5);
t47 = t59 * t100;
t88 = -qJD(5) * t11 - qJDD(1);
t2 = -t59 * t3 + t88 * t61 + t47;
t8 = -t59 * qJD(1) + t61 * t11;
t71 = -(t59 * t8 + t61 * t7) * qJD(5) + t1 * t61 - t2 * t59;
t97 = pkin(8) + qJ(2);
t49 = cos(t97);
t91 = sin(t97);
t18 = -t49 * t62 - t91 * t60;
t19 = t49 * t60 - t91 * t62;
t84 = g(1) * t18 + g(2) * t19;
t66 = t71 + t84;
t119 = t54 * pkin(4);
t27 = -t60 * qJ(3) + t62 * t63;
t14 = t62 * qJD(3) + t27 * qJD(4);
t117 = t14 * t98;
t15 = t60 * qJD(3) + t28 * qJD(4);
t116 = t15 * t98;
t113 = t59 * t61;
t112 = t61 * t54;
t110 = t49 * pkin(2) + t91 * qJ(3);
t109 = g(1) * t91 - g(2) * t49;
t56 = t59 ^ 2;
t57 = t61 ^ 2;
t108 = t56 - t57;
t107 = t56 + t57;
t105 = pkin(2) * qJDD(2);
t96 = t125 * t113;
t95 = t49 * pkin(3) + t110;
t94 = 0.2e1 * t99;
t93 = t107 * t54;
t89 = qJDD(3) - t105;
t87 = -0.2e1 * t103 * t113;
t86 = t121 * t60 - t123 * t62;
t85 = g(1) * t19 - g(2) * t18;
t81 = t7 * t59 - t8 * t61;
t80 = -g(3) - t88;
t4 = t86 + t119;
t79 = -t4 + t85;
t78 = -t91 * pkin(2) + t49 * qJ(3);
t77 = -t3 - t84 + t124;
t76 = g(1) * t49 + g(2) * t91;
t75 = t85 - t86;
t74 = -pkin(7) * qJDD(5) + (t10 + t16 + t118) * qJD(5);
t22 = pkin(4) - t27;
t23 = -pkin(7) + t28;
t73 = -qJDD(5) * t23 + (-t22 * t98 - t10 - t14) * qJD(5);
t72 = -qJDD(5) * t60 + 0.2e1 * t62 * t103;
t70 = -t91 * pkin(3) + t78;
t69 = pkin(7) * t64 + t114 + t119 - t79;
t68 = -t22 * t54 + t23 * t64 - t116 + t79;
t67 = -t5 - t84;
t65 = qJD(2) ^ 2;
t58 = qJDD(1) - g(3);
t30 = qJDD(5) * t61 - t64 * t59;
t29 = qJDD(5) * t59 + t64 * t61;
t13 = t57 * t54 + t87;
t12 = -t56 * t54 + t87;
t9 = t108 * t103 - t59 * t112;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, -t30, t29, 0, t81 * qJD(5) - t1 * t59 - t2 * t61 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t109, t76, 0, 0, 0, 0, 0, qJDD(2), 0, 0, -qJDD(3) + 0.2e1 * t105 + t109, 0, -t76 + t94 + 0.2e1 * t101, -t89 * pkin(2) - g(1) * t78 - g(2) * t110 + (t94 + t101) * qJ(3), 0, 0, 0, 0, 0, t54, -t27 * t54 + t116 - t75, t28 * t54 + t117 - t67, 0, -g(1) * t70 - g(2) * t95 + t17 * t14 - t16 * t15 - t27 * t86 + t5 * t28, -t12, -0.2e1 * t9, -t29, t13, -t30, 0, t73 * t59 - t68 * t61, t68 * t59 + t73 * t61, -t107 * t117 - t23 * t93 - t66, t4 * t22 + t10 * t15 - g(1) * (t19 * pkin(4) + t18 * pkin(7) + t70) - g(2) * (-t18 * pkin(4) + t19 * pkin(7) + t95) - t81 * t14 + t71 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t65, -t65 * qJ(3) - t109 + t89, 0, 0, 0, 0, 0, 0, -t125 * t60 - t111, t60 * t54 - t120, 0, (-t86 - t114) * t62 + (t5 + t115) * t60 - t109, 0, 0, 0, 0, 0, 0, -t126 * t61 + t72 * t59, t126 * t59 + t72 * t61, t107 * t120 - t60 * t93, (t98 * t81 - t4) * t62 + (t71 - t124) * t60 - t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t75 - t114, t67 - t115, 0, 0, t12, 0.2e1 * t9, t29, -t13, t30, 0, t74 * t59 - t69 * t61, t69 * t59 + t74 * t61, -pkin(7) * t93 + t107 * t115 + t66, t79 * pkin(4) + t66 * pkin(7) - t10 * t17 + t81 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, t108 * t125, -t59 * t54, t96, -t112, qJDD(5), t8 * qJD(5) + t77 * t59 - t80 * t61 + t47, t106 + t80 * t59 + (t77 + t100) * t61, 0, 0;];
tau_reg = t6;
