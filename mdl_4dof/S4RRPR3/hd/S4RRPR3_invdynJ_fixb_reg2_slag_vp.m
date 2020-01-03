% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPR3
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
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:36
% EndTime: 2019-12-31 17:01:37
% DurationCPUTime: 0.54s
% Computational Cost: add. (848->142), mult. (1487->188), div. (0->0), fcn. (853->12), ass. (0->108)
t67 = qJ(1) + qJ(2);
t58 = pkin(7) + t67;
t48 = cos(t58);
t128 = g(2) * t48;
t63 = qJDD(1) + qJDD(2);
t71 = sin(qJ(2));
t105 = qJDD(1) * t71;
t130 = pkin(1) * t71;
t101 = qJD(1) * t130;
t74 = cos(qJ(2));
t123 = t74 * pkin(1);
t56 = qJDD(1) * t123;
t20 = t63 * pkin(2) - qJD(2) * t101 + t56;
t68 = sin(pkin(7));
t69 = cos(pkin(7));
t108 = qJD(1) * t74;
t98 = qJD(2) * t108;
t7 = -(t98 + t105) * pkin(1) * t68 + t69 * t20;
t5 = -t63 * pkin(3) - t7;
t100 = -t5 - t128;
t64 = qJD(1) + qJD(2);
t32 = pkin(1) * t108 + t64 * pkin(2);
t17 = t69 * t101 + t68 * t32;
t15 = t64 * pkin(6) + t17;
t73 = cos(qJ(4));
t116 = t73 * t15;
t70 = sin(qJ(4));
t10 = t70 * qJD(3) + t116;
t9 = t73 * qJD(3) - t70 * t15;
t109 = t9 * qJD(4);
t117 = t69 * t71;
t46 = pkin(1) * t117;
t8 = t69 * pkin(1) * t98 + qJDD(1) * t46 + t68 * t20;
t6 = t63 * pkin(6) + t8;
t2 = t70 * qJDD(3) + t73 * t6 + t109;
t57 = t73 * qJDD(3);
t3 = -t10 * qJD(4) - t70 * t6 + t57;
t80 = t2 * t73 - t3 * t70 + (-t10 * t70 - t73 * t9) * qJD(4);
t47 = sin(t58);
t113 = -g(1) * t48 - g(2) * t47;
t59 = sin(t67);
t60 = cos(t67);
t132 = g(1) * t59 - g(2) * t60;
t129 = pkin(2) * t59;
t44 = g(1) * t47;
t126 = t68 * pkin(2);
t125 = t69 * pkin(2);
t72 = sin(qJ(1));
t124 = t72 * pkin(1);
t89 = pkin(1) * (t68 * t74 + t117);
t25 = qJD(1) * t89;
t122 = t25 * t64;
t26 = qJD(2) * t89;
t121 = t26 * t64;
t118 = t68 * t71;
t88 = pkin(1) * (t69 * t74 - t118);
t27 = qJD(1) * t88;
t120 = t27 * t64;
t28 = qJD(2) * t88;
t119 = t28 * t64;
t115 = t73 * t63;
t16 = -t68 * t101 + t69 * t32;
t14 = -t64 * pkin(3) - t16;
t114 = t14 * qJD(4) * t70 + t73 * t44;
t55 = pkin(2) + t123;
t30 = t68 * t55 + t46;
t112 = g(1) * t60 + g(2) * t59;
t65 = t70 ^ 2;
t66 = t73 ^ 2;
t111 = t65 - t66;
t110 = t65 + t66;
t107 = qJD(4) * t73;
t106 = qJDD(3) - g(3);
t104 = -t100 * t70 + t14 * t107;
t62 = t64 ^ 2;
t103 = t70 * t62 * t73;
t54 = pkin(2) * t60;
t102 = t48 * pkin(3) + t47 * pkin(6) + t54;
t99 = t110 * t63;
t97 = qJD(1) * (-qJD(2) + t64);
t96 = qJD(2) * (-qJD(1) - t64);
t95 = t70 * t64 * t107;
t94 = t56 + t132;
t93 = -t8 - t113;
t75 = cos(qJ(1));
t92 = g(1) * t72 - g(2) * t75;
t91 = t10 * t73 - t9 * t70;
t90 = -t47 * pkin(3) + t48 * pkin(6) - t129;
t29 = -pkin(1) * t118 + t69 * t55;
t23 = -pkin(3) - t29;
t24 = pkin(6) + t30;
t76 = qJD(4) ^ 2;
t86 = t23 * t63 + t24 * t76 + t121;
t49 = pkin(6) + t126;
t50 = -pkin(3) - t125;
t85 = t49 * t76 + t50 * t63 - t122;
t84 = -qJDD(4) * t24 + (t23 * t64 - t28) * qJD(4);
t83 = -qJDD(4) * t49 + (t50 * t64 + t27) * qJD(4);
t81 = -qJD(3) * qJD(4) - t14 * t64 - t113 - t6;
t79 = t113 + t80;
t78 = t44 + t7 - t128;
t61 = t75 * pkin(1);
t34 = qJDD(4) * t73 - t76 * t70;
t33 = qJDD(4) * t70 + t76 * t73;
t22 = t66 * t63 - 0.2e1 * t95;
t21 = t65 * t63 + 0.2e1 * t95;
t13 = -0.2e1 * t111 * t64 * qJD(4) + 0.2e1 * t70 * t115;
t1 = [0, 0, 0, 0, 0, qJDD(1), t92, g(1) * t75 + g(2) * t72, 0, 0, 0, 0, 0, 0, 0, t63, (t63 * t74 + t71 * t96) * pkin(1) + t94, ((-qJDD(1) - t63) * t71 + t74 * t96) * pkin(1) + t112, 0, (t92 + (t71 ^ 2 + t74 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t63, t29 * t63 - t121 + t78, -t30 * t63 - t119 + t93, 0, t8 * t30 + t17 * t28 + t7 * t29 - t16 * t26 - g(1) * (-t124 - t129) - g(2) * (t54 + t61), t21, t13, t33, t22, t34, 0, t84 * t70 + (t100 - t86) * t73 + t114, t84 * t73 + (t86 - t44) * t70 + t104, t110 * t119 + t24 * t99 + t79, t5 * t23 + t14 * t26 - g(1) * (t90 - t124) - g(2) * (t61 + t102) + t91 * t28 + t80 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t97 * t130 + t94, (t74 * t97 - t105) * pkin(1) + t112, 0, 0, 0, 0, 0, 0, 0, t63, t63 * t125 + t122 + t78, -t63 * t126 + t120 + t93, 0, t16 * t25 - t17 * t27 + (t68 * t8 + t69 * t7 + t132) * pkin(2), t21, t13, t33, t22, t34, 0, t83 * t70 + (t100 - t85) * t73 + t114, t83 * t73 + (t85 - t44) * t70 + t104, -t110 * t120 + t49 * t99 + t79, -g(1) * t90 - g(2) * t102 - t14 * t25 - t91 * t27 + t80 * t49 + t5 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, 0, 0, 0, 0, 0, t34, -t33, 0, t91 * qJD(4) + t2 * t70 + t3 * t73 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, t111 * t62, t70 * t63, t103, t115, qJDD(4), -g(3) * t73 + t57 + (t10 - t116) * qJD(4) + t81 * t70, t109 + (qJD(4) * t15 - t106) * t70 + t81 * t73, 0, 0;];
tau_reg = t1;
