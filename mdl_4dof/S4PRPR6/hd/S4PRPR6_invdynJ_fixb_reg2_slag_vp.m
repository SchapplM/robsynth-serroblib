% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRPR6
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:44
% EndTime: 2019-12-31 16:24:45
% DurationCPUTime: 0.66s
% Computational Cost: add. (667->146), mult. (1557->204), div. (0->0), fcn. (1143->10), ass. (0->92)
t62 = cos(qJ(2));
t50 = g(3) * t62;
t60 = sin(qJ(2));
t55 = sin(pkin(6));
t57 = cos(pkin(6));
t80 = g(1) * t57 + g(2) * t55;
t67 = t80 * t60 - t50;
t96 = t62 * qJD(1);
t81 = qJD(3) - t96;
t56 = cos(pkin(7));
t61 = cos(qJ(4));
t102 = t61 * t56;
t54 = sin(pkin(7));
t59 = sin(qJ(4));
t103 = t59 * t54;
t33 = -t102 + t103;
t23 = t33 * t60;
t115 = t80 * t62;
t34 = t61 * t54 + t59 * t56;
t28 = t34 * qJD(2);
t91 = qJD(1) * qJD(2);
t92 = t62 * qJDD(1);
t75 = t60 * t91 + qJDD(3) - t92;
t98 = qJDD(2) * pkin(2);
t25 = t75 - t98;
t113 = (t80 + t91) * t60 - t25 + t98 - t50;
t112 = t28 ^ 2;
t109 = g(3) * t60;
t101 = pkin(5) + qJ(3);
t35 = t101 * t54;
t36 = t101 * t56;
t14 = -t61 * t35 - t59 * t36;
t70 = t33 * t62;
t108 = qJD(1) * t70 - t33 * qJD(3) + t14 * qJD(4);
t15 = -t59 * t35 + t61 * t36;
t107 = -t15 * qJD(4) - t81 * t34;
t88 = qJD(2) * t102;
t89 = qJD(2) * t103;
t26 = -t88 + t89;
t106 = t28 * t26;
t105 = t55 * t62;
t104 = t57 * t62;
t51 = t54 ^ 2;
t52 = t56 ^ 2;
t100 = t51 + t52;
t99 = qJD(2) * t60;
t97 = t60 * qJD(1);
t95 = qJDD(1) - g(3);
t94 = t54 * qJDD(2);
t93 = t56 * qJDD(2);
t90 = qJD(4) * t88 + t59 * t93 + t61 * t94;
t46 = t56 * pkin(3) + pkin(2);
t21 = qJDD(2) * qJ(3) + t60 * qJDD(1) + (qJD(3) + t96) * qJD(2);
t87 = t100 * t21;
t86 = t100 * t62;
t83 = pkin(5) * qJDD(2) + t21;
t12 = t83 * t54;
t13 = t83 * t56;
t85 = -t61 * t12 - t59 * t13;
t39 = qJD(2) * qJ(3) + t97;
t84 = pkin(5) * qJD(2) + t39;
t82 = t100 * qJDD(2);
t79 = t59 * t94 - t61 * t93;
t78 = -t59 * t12 + t61 * t13;
t19 = t84 * t54;
t20 = t84 * t56;
t5 = -t61 * t19 - t59 * t20;
t6 = -t59 * t19 + t61 * t20;
t76 = t46 * qJDD(2);
t63 = qJD(2) ^ 2;
t74 = t62 * qJDD(2) - t63 * t60;
t22 = t34 * t60;
t31 = t34 * qJD(4);
t68 = -t109 - t115;
t66 = t81 * t100;
t65 = t75 - t67;
t64 = t87 + t68;
t53 = pkin(7) + qJ(4);
t49 = cos(t53);
t48 = sin(t53);
t37 = -qJD(2) * pkin(2) + t81;
t32 = -t46 * qJD(2) + t81;
t30 = t33 * qJD(4);
t24 = t26 ^ 2;
t16 = -t76 + t75;
t10 = qJD(2) * t31 + t79;
t9 = qJD(4) * t89 - t90;
t8 = qJD(4) * t23 - t62 * t28;
t7 = -qJD(2) * t70 - qJD(4) * t22;
t2 = -t6 * qJD(4) + t85;
t1 = t5 * qJD(4) + t78;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, t74, -qJDD(2) * t60 - t63 * t62, 0, -g(3) + (t60 ^ 2 + t62 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t74 * t56, -t74 * t54, t60 * t82 + t63 * t86, -t25 * t62 - g(3) + t60 * t87 + (t37 * t60 + t39 * t86) * qJD(2), 0, 0, 0, 0, 0, 0, t8 * qJD(4) - t22 * qJDD(4) - t62 * t10 + t26 * t99, -t7 * qJD(4) + t23 * qJDD(4) + t28 * t99 + t62 * t9, t23 * t10 - t22 * t9 - t7 * t26 - t8 * t28, -t1 * t23 - t16 * t62 - t2 * t22 + t32 * t99 + t5 * t8 + t6 * t7 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t92 + t67, -t95 * t60 + t115, 0, 0, t51 * qJDD(2), 0.2e1 * t54 * t93, 0, t52 * qJDD(2), 0, 0, t113 * t56, -t113 * t54, qJ(3) * t82 + t66 * qJD(2) + t64, -t37 * t97 + (-t25 + t67) * pkin(2) + t64 * qJ(3) + t66 * t39, -t28 * t30 - t9 * t34, -t34 * t10 + t30 * t26 - t28 * t31 + t9 * t33, -t30 * qJD(4) + t34 * qJDD(4), t10 * t33 + t26 * t31, -t31 * qJD(4) - t33 * qJDD(4), 0, t107 * qJD(4) + t14 * qJDD(4) - t46 * t10 + t16 * t33 - t26 * t97 + t32 * t31 + t67 * t49, -t108 * qJD(4) - t15 * qJDD(4) + t16 * t34 - t28 * t97 - t32 * t30 + t46 * t9 - t48 * t67, -t1 * t33 - t15 * t10 - t107 * t28 - t108 * t26 + t14 * t9 - t2 * t34 + t5 * t30 - t6 * t31 + t68, t1 * t15 + t2 * t14 - t16 * t46 - t32 * t97 - g(3) * (t101 * t60 + t62 * t46) + t108 * t6 + t107 * t5 + t80 * (-t101 * t62 + t46 * t60); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t94, -t100 * t63, -t100 * t39 * qJD(2) + t65 - t98, 0, 0, 0, 0, 0, 0, 0.2e1 * t28 * qJD(4) + t79, (-t26 - t89) * qJD(4) + t90, -t24 - t112, t6 * t26 + t5 * t28 + t65 - t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, -t24 + t112, (t26 - t89) * qJD(4) + t90, -t106, -t79, qJDD(4), -t32 * t28 - g(1) * (-t48 * t104 + t55 * t49) - g(2) * (-t48 * t105 - t57 * t49) + t48 * t109 + t85, t32 * t26 - g(1) * (-t49 * t104 - t55 * t48) - g(2) * (-t49 * t105 + t57 * t48) + t49 * t109 - t78, 0, 0;];
tau_reg = t3;
