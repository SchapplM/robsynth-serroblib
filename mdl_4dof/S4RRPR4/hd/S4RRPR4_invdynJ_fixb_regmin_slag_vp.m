% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR4
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:34
% EndTime: 2019-12-31 17:02:36
% DurationCPUTime: 0.59s
% Computational Cost: add. (690->135), mult. (1028->174), div. (0->0), fcn. (669->12), ass. (0->92)
t70 = sin(pkin(7));
t71 = cos(pkin(7));
t111 = t70 ^ 2 + t71 ^ 2;
t69 = qJ(1) + qJ(2);
t61 = cos(t69);
t122 = g(2) * t61;
t60 = sin(t69);
t54 = g(1) * t60;
t126 = t122 - t54;
t64 = qJDD(1) + qJDD(2);
t121 = t64 * pkin(2);
t110 = pkin(1) * qJD(1);
t73 = sin(qJ(2));
t101 = t73 * t110;
t76 = cos(qJ(2));
t120 = t76 * pkin(1);
t112 = -qJD(2) * t101 + qJDD(1) * t120;
t98 = qJDD(3) - t112;
t21 = t98 - t121;
t125 = t21 + t122;
t75 = cos(qJ(4));
t116 = t75 * t71;
t72 = sin(qJ(4));
t119 = t72 * t70;
t29 = -t116 + t119;
t30 = t75 * t70 + t72 * t71;
t68 = qJD(1) + qJD(2);
t24 = t30 * t68;
t124 = g(1) * t61 + g(2) * t60;
t123 = t112 + t54;
t26 = t30 * qJD(4);
t109 = qJD(1) * t76;
t87 = -pkin(1) * t109 + qJD(3);
t115 = t125 * t70;
t114 = t61 * pkin(2) + t60 * qJ(3);
t108 = qJD(2) * t73;
t107 = qJD(2) * t76;
t106 = qJDD(1) * t73;
t105 = t68 * t119;
t104 = t68 * t116;
t103 = qJD(4) * t104 + t30 * t64;
t102 = pkin(1) * t108;
t99 = t68 * t108;
t51 = -t71 * pkin(3) - pkin(2);
t97 = -t60 * pkin(2) + t61 * qJ(3);
t16 = t64 * qJ(3) + t68 * qJD(3) + (qJD(1) * t107 + t106) * pkin(1);
t96 = pkin(6) * t64 + t16;
t45 = pkin(1) * t107 + qJD(3);
t94 = t111 * t45;
t50 = t73 * pkin(1) + qJ(3);
t93 = t111 * t50;
t92 = qJ(3) * t111;
t91 = t111 * t16 - t124;
t90 = t68 * t101;
t89 = t29 * t64;
t27 = (-pkin(6) - t50) * t70;
t62 = t71 * pkin(6);
t28 = t71 * t50 + t62;
t86 = t75 * t27 - t72 * t28;
t85 = t72 * t27 + t75 * t28;
t37 = (-pkin(6) - qJ(3)) * t70;
t38 = t71 * qJ(3) + t62;
t84 = t75 * t37 - t72 * t38;
t83 = t72 * t37 + t75 * t38;
t13 = t51 * t64 + t98;
t20 = t51 * t68 + t87;
t25 = t29 * qJD(4);
t67 = pkin(7) + qJ(4);
t58 = sin(t67);
t82 = t126 * t58 + t13 * t30 - t20 * t25;
t59 = cos(t67);
t81 = -t126 * t59 + t13 * t29 + t20 * t26;
t80 = t90 - t122;
t56 = -pkin(2) - t120;
t79 = pkin(1) * t99 + t56 * t64;
t78 = t87 * t111;
t77 = cos(qJ(1));
t74 = sin(qJ(1));
t44 = t71 * t54;
t36 = t51 - t120;
t32 = t68 * qJ(3) + t101;
t31 = -t68 * pkin(2) + t87;
t22 = -t104 + t105;
t12 = -t26 * qJD(4) - t29 * qJDD(4);
t11 = -t25 * qJD(4) + t30 * qJDD(4);
t8 = t96 * t71;
t7 = t96 * t70;
t4 = t26 * t68 + t89;
t3 = -qJD(4) * t105 + t103;
t2 = -t24 * t25 + t3 * t30;
t1 = t25 * t22 - t24 * t26 - t3 * t29 - t30 * t4;
t5 = [qJDD(1), g(1) * t74 - g(2) * t77, g(1) * t77 + g(2) * t74, t64, -t122 + (t64 * t76 - t99) * pkin(1) + t123, ((-qJDD(1) - t64) * t73 + (-qJD(1) - t68) * t107) * pkin(1) + t124, t44 + (-t79 - t125) * t71, (t79 - t54) * t70 + t115, t64 * t93 + t68 * t94 + t91, t21 * t56 + t31 * t102 - g(1) * (-t74 * pkin(1) + t97) - g(2) * (t77 * pkin(1) + t114) + t32 * t94 + t16 * t93, t2, t1, t11, t12, 0, t22 * t102 + t36 * t4 + t86 * qJDD(4) + (-qJD(4) * t85 - t30 * t45) * qJD(4) + t81, t24 * t102 + t36 * t3 - t85 * qJDD(4) + (-qJD(4) * t86 + t29 * t45) * qJD(4) + t82; 0, 0, 0, t64, t80 + t123, (-t106 + (-qJD(2) + t68) * t109) * pkin(1) + t124, t44 + (-t21 + t80 + t121) * t71, (-t121 - t90 - t54) * t70 + t115, t64 * t92 + t68 * t78 + t91, -t21 * pkin(2) - g(1) * t97 - g(2) * t114 - t31 * t101 + t16 * t92 + t78 * t32, t2, t1, t11, t12, 0, t51 * t4 + t84 * qJDD(4) + (-qJD(3) * t30 - qJD(4) * t83) * qJD(4) + (-t73 * t22 + t76 * t26) * t110 + t81, t51 * t3 - t83 * qJDD(4) + (qJD(3) * t29 - qJD(4) * t84) * qJD(4) + (-t73 * t24 - t25 * t76) * t110 + t82; 0, 0, 0, 0, 0, 0, -t71 * t64, t70 * t64, -t111 * t68 ^ 2, -t111 * t68 * t32 + t125 - t54, 0, 0, 0, 0, 0, 0.2e1 * t24 * qJD(4) + t89, (-t22 - t105) * qJD(4) + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t22, -t22 ^ 2 + t24 ^ 2, (t22 - t105) * qJD(4) + t103, -t89, qJDD(4), -g(3) * t59 + t124 * t58 - t20 * t24 - t75 * t7 - t72 * t8, g(3) * t58 + t124 * t59 + t20 * t22 + t72 * t7 - t75 * t8;];
tau_reg = t5;
