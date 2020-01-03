% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPPR6
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:46
% EndTime: 2019-12-31 16:40:47
% DurationCPUTime: 0.60s
% Computational Cost: add. (562->157), mult. (1339->209), div. (0->0), fcn. (925->6), ass. (0->95)
t97 = (qJD(1) * qJD(2));
t98 = (qJ(2) * qJDD(1));
t119 = t97 + t98;
t67 = cos(pkin(6));
t66 = sin(pkin(6));
t104 = t66 * qJ(3);
t88 = pkin(1) + t104;
t27 = (pkin(2) + pkin(3)) * t67 + t88;
t96 = qJD(1) * qJD(3);
t89 = t66 * t96;
t9 = qJDD(1) * t27 - qJDD(2) + t89;
t69 = sin(qJ(1));
t71 = cos(qJ(1));
t105 = g(1) * t71 + g(2) * t69;
t68 = sin(qJ(4));
t70 = cos(qJ(4));
t30 = t66 * t68 + t67 * t70;
t24 = t30 * qJD(1);
t111 = t67 * t68;
t31 = t66 * t70 - t111;
t95 = qJD(1) * t111;
t76 = qJD(4) * t95 - qJDD(1) * t30;
t117 = t24 ^ 2;
t103 = qJD(1) * t66;
t94 = t70 * t103;
t26 = t94 - t95;
t116 = t26 ^ 2;
t114 = g(1) * t69;
t62 = g(2) * t71;
t113 = t26 * t24;
t110 = t67 * t69;
t109 = t67 * t71;
t108 = -pkin(5) + qJ(2);
t65 = t67 ^ 2;
t56 = t65 * qJDD(1);
t74 = qJ(2) ^ 2;
t93 = 2 * t97;
t84 = qJ(2) * t93;
t107 = t74 * t56 + t65 * t84;
t106 = t71 * pkin(1) + t69 * qJ(2);
t102 = qJD(3) * t66;
t101 = qJDD(1) * pkin(1);
t38 = qJ(2) * t103 + qJD(3);
t33 = -t67 * pkin(2) - t88;
t100 = qJDD(1) * t33;
t55 = t66 * qJDD(1);
t99 = t67 * qJDD(1);
t92 = t66 * t99;
t28 = t119 * t66 + qJDD(3);
t58 = t71 * qJ(2);
t91 = -t69 * pkin(1) + t58;
t90 = -t62 + t114;
t36 = t108 * t67;
t64 = t66 ^ 2;
t73 = qJD(1) ^ 2;
t34 = (-t64 - t65) * t73;
t15 = -pkin(5) * t55 + t28;
t16 = (qJDD(1) * t108 + t97) * t67;
t87 = t70 * t15 - t68 * t16;
t86 = -t105 + (t93 + 2 * t98) * t65;
t85 = pkin(2) * t109 + t71 * t104 + t106;
t53 = qJDD(2) - t101;
t82 = -t70 * t55 + t68 * t99;
t81 = t68 * t15 + t70 * t16;
t29 = -pkin(5) * t103 + t38;
t32 = qJD(1) * t36;
t7 = t70 * t29 - t68 * t32;
t8 = t68 * t29 + t70 * t32;
t35 = t108 * t66;
t10 = t70 * t35 - t68 * t36;
t11 = t68 * t35 + t70 * t36;
t79 = t101 - t53 - t62;
t75 = qJDD(2) + t100;
t13 = t75 - t89;
t78 = -t100 - t13 - t62;
t22 = t30 * qJD(4);
t77 = t119 * t64;
t72 = qJD(4) ^ 2;
t59 = g(3) * t67;
t54 = t64 * qJDD(1);
t47 = g(1) * t110;
t23 = t31 * qJD(4);
t21 = qJD(1) * t33 + qJD(2);
t20 = t30 * t71;
t19 = t31 * t71;
t18 = t30 * t69;
t17 = t31 * t69;
t14 = qJD(1) * t27 - qJD(2);
t6 = qJD(4) * t94 - t76;
t5 = qJD(1) * t22 + t82;
t4 = qJD(2) * t31 - qJD(4) * t11;
t3 = qJD(2) * t30 + qJD(4) * t10;
t2 = -t8 * qJD(4) + t87;
t1 = t7 * qJD(4) + t81;
t12 = [0, 0, 0, 0, 0, qJDD(1), t90, t105, 0, 0, t54, 0.2e1 * t92, 0, t56, 0, 0, t67 * t79 + t47, (-t79 - t114) * t66, 0.2e1 * t77 + t86, -t53 * pkin(1) - g(1) * t91 - g(2) * t106 + (qJDD(1) * t74 + t84) * t64 + t107, t54, 0, -0.2e1 * t92, 0, 0, t56, t47 + (t78 + t89) * t67, t28 * t66 + t77 + t86, t64 * t96 + (t78 + t114) * t66, t13 * t33 - g(1) * (-pkin(2) * t110 + t91) - g(2) * t85 + (t28 * qJ(2) + qJ(3) * t114 + t38 * qJD(2) - t21 * qJD(3)) * t66 + t107, -t26 * t22 - t5 * t31, t22 * t24 - t26 * t23 + t5 * t30 - t31 * t6, -t22 * qJD(4) + t31 * qJDD(4), t24 * t23 + t6 * t30, -t23 * qJD(4) - t30 * qJDD(4), 0, g(1) * t18 - g(2) * t20 + t4 * qJD(4) + t10 * qJDD(4) + t102 * t24 + t14 * t23 + t27 * t6 + t9 * t30, g(1) * t17 - g(2) * t19 - t3 * qJD(4) - t11 * qJDD(4) + t102 * t26 - t14 * t22 - t27 * t5 + t9 * t31, -t1 * t30 + t10 * t5 - t11 * t6 - t2 * t31 + t7 * t22 - t8 * t23 - t3 * t24 - t4 * t26 + t105, t1 * t11 + t8 * t3 + t2 * t10 + t7 * t4 + t9 * t27 + t14 * t102 - g(1) * (-t71 * pkin(5) + t58) - g(2) * (pkin(3) * t109 + t85) + (-g(1) * (-t67 * pkin(3) + t33) + g(2) * pkin(5)) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t55, t34, qJ(2) * t34 + t53 - t90, 0, 0, 0, 0, 0, 0, -t99, t34, -t55, -t65 * t73 * qJ(2) + (-qJD(3) - t38) * t103 + t75 - t90, 0, 0, 0, 0, 0, 0, (-t26 - t94) * qJD(4) + t76, 0.2e1 * t24 * qJD(4) + t82, t116 + t117, -t8 * t24 - t7 * t26 - t9 - t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66 * t73 * t67, t55, -t64 * t73, t59 + (qJD(1) * t21 - t105) * t66 + t28, 0, 0, 0, 0, 0, 0, qJDD(4) * t70 - t103 * t24 - t72 * t68, -qJDD(4) * t68 - t103 * t26 - t72 * t70, t70 * t5 - t68 * t6 + (-t24 * t70 + t26 * t68) * qJD(4), t1 * t68 + t2 * t70 + t59 + (-t68 * t7 + t70 * t8) * qJD(4) + (-qJD(1) * t14 - t105) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t116 - t117, -t82, -t113, (t26 - t94) * qJD(4) + t76, qJDD(4), -g(1) * t19 - g(2) * t17 + g(3) * t30 - t14 * t26 + t87, g(1) * t20 + g(2) * t18 + g(3) * t31 + t14 * t24 - t81, 0, 0;];
tau_reg = t12;
