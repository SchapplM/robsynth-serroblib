% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:04:41
% EndTime: 2021-01-15 17:04:43
% DurationCPUTime: 0.59s
% Computational Cost: add. (606->153), mult. (969->165), div. (0->0), fcn. (491->8), ass. (0->96)
t53 = cos(pkin(7));
t34 = -pkin(1) * t53 - pkin(2);
t27 = -pkin(6) + t34;
t17 = qJD(1) * t27 + qJD(3);
t91 = qJ(5) * qJD(1);
t115 = t17 - t91;
t57 = cos(qJ(4));
t116 = t115 * t57;
t52 = sin(pkin(7));
t31 = t52 * pkin(1) + qJ(3);
t101 = t31 * qJDD(1);
t46 = qJ(1) + pkin(7);
t39 = sin(t46);
t107 = g(1) * t39;
t40 = cos(t46);
t33 = g(2) * t40;
t81 = -t33 + t107;
t55 = sin(qJ(4));
t97 = qJ(5) - t27;
t12 = t97 * t55;
t114 = qJD(4) * t12;
t51 = qJDD(2) - g(3);
t48 = qJD(3) * qJD(1);
t18 = t48 + t101;
t86 = qJD(1) * qJD(4);
t79 = t57 * t86;
t88 = t55 * qJDD(1);
t6 = qJDD(5) + t18 + (t79 + t88) * pkin(4);
t73 = -g(1) * t40 - g(2) * t39;
t113 = t6 + t73;
t105 = t55 * pkin(4);
t19 = t31 + t105;
t96 = qJD(1) * t19;
t11 = qJD(5) + t96;
t112 = -t11 * qJD(1) - t81;
t111 = (qJD(5) + t11) * qJD(1);
t110 = t34 * qJDD(1);
t22 = qJD(1) * t31;
t109 = 0.2e1 * qJD(4) * t22 + qJDD(4) * t27;
t4 = -t55 * qJD(2) + t116;
t98 = qJD(4) * pkin(4);
t3 = t4 + t98;
t108 = t3 - t4;
t44 = g(3) * t55;
t104 = t39 * t55;
t103 = t55 * t17;
t60 = qJD(1) ^ 2;
t102 = t60 * t55;
t49 = t55 ^ 2;
t50 = t57 ^ 2;
t100 = -t49 - t50;
t99 = t49 - t50;
t13 = t97 * t57;
t95 = qJD(4) * t13;
t94 = qJDD(4) * pkin(4);
t93 = t22 * qJD(1);
t89 = qJDD(4) * t55;
t41 = t57 * qJDD(1);
t87 = qJ(5) * qJDD(1);
t85 = qJD(1) * qJD(5);
t84 = qJD(2) * qJD(4);
t58 = cos(qJ(1));
t83 = pkin(1) * t58 + pkin(2) * t40 + qJ(3) * t39;
t56 = sin(qJ(1));
t82 = -t56 * pkin(1) + qJ(3) * t40;
t80 = t55 * t86;
t78 = t11 + t96;
t59 = qJD(4) ^ 2;
t21 = qJDD(4) * t57 - t55 * t59;
t75 = -0.2e1 * t80;
t16 = qJDD(1) * t27 + qJDD(3);
t74 = -t16 + t87;
t72 = g(1) * t56 - g(2) * t58;
t5 = qJD(2) * t57 - t55 * t91 + t103;
t71 = t3 * t57 + t5 * t55;
t69 = -t93 - t107;
t9 = t57 * t16;
t68 = -t55 * qJDD(2) + t33 * t57 + t44 + t9;
t36 = t55 * t84;
t67 = g(1) * t104 - t51 * t57 + t36;
t66 = -t84 - t87;
t65 = qJDD(3) + t110;
t64 = t73 + t101;
t25 = t57 * t98 + qJD(3);
t63 = qJD(1) * t25 + qJDD(1) * t19 + t113;
t62 = -t27 * t59 + t18 + t48 + t64;
t54 = -qJ(5) - pkin(6);
t26 = qJ(5) * t80;
t20 = -t57 * t59 - t89;
t15 = t21 - t102;
t14 = -t89 + (-t59 - t60) * t57;
t8 = -t55 * qJD(5) - t95;
t7 = -t57 * qJD(5) + t114;
t2 = -t36 + (qJD(4) * t115 + qJDD(2)) * t57 + (-t74 - t85) * t55;
t1 = t94 + t26 + t9 + (-qJD(4) * t17 - qJDD(2)) * t55 + (t66 - t85) * t57;
t10 = [qJDD(1), t72, g(1) * t58 + g(2) * t56, (t72 + (t52 ^ 2 + t53 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(3) - t81 + 0.2e1 * t110, 0.2e1 * t48 + t64 + t101, t18 * t31 + t22 * qJD(3) + t65 * t34 - g(1) * (-pkin(2) * t39 + t82) - g(2) * t83, qJDD(1) * t50 + t57 * t75, -0.2e1 * t41 * t55 + 0.2e1 * t86 * t99, t21, t20, 0, t109 * t57 + t62 * t55, -t109 * t55 + t62 * t57, -t13 * qJDD(4) + (t57 * t78 + t7) * qJD(4) + t63 * t55, t12 * qJDD(4) + (-t55 * t78 - t8) * qJD(4) + t63 * t57, (-t5 * qJD(4) + qJDD(1) * t13 - t1 + (-t7 + t114) * qJD(1)) * t57 + (qJD(4) * t3 + qJDD(1) * t12 - t2 + (-t8 - t95) * qJD(1)) * t55 + t81, -t2 * t12 + t5 * t8 - t1 * t13 + t3 * t7 + t6 * t19 + t11 * t25 - g(1) * (t40 * t105 + (-pkin(2) + t54) * t39 + t82) - g(2) * (pkin(4) * t104 - t40 * t54 + t83); 0, 0, 0, t51, 0, 0, t51, 0, 0, 0, 0, 0, t20, -t21, t20, -t21, 0, -qJD(4) * t71 - t1 * t55 + t2 * t57 - g(3); 0, 0, 0, 0, qJDD(1), -t60, t33 + t65 + t69, 0, 0, 0, 0, 0, t15, t14, t15, t14, t100 * qJDD(1), t1 * t57 + t2 * t55 + (-t3 * t55 + t5 * t57) * qJD(4) + t112; 0, 0, 0, 0, 0, 0, 0, t57 * t102, -t99 * t60, t41, -t88, qJDD(4), t57 * t69 + t68, (-t16 - t84 + t93 - t33) * t55 + t67, 0.2e1 * t94 + t26 + (t5 - t103) * qJD(4) + (-pkin(4) * t102 - t107 - t111 + t66) * t57 + t68, -t50 * t60 * pkin(4) + (t4 - t116) * qJD(4) + (-t33 + t74 + t111) * t55 + t67, -pkin(4) * t41 + (t98 - t108) * t55 * qJD(1), t108 * t5 + (t112 * t57 + t1 + t44) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t79 + t88, t41 + t75, t100 * t60, qJD(1) * t71 + t113;];
tau_reg = t10;
