% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:57
% EndTime: 2019-12-31 17:36:59
% DurationCPUTime: 0.49s
% Computational Cost: add. (402->128), mult. (776->163), div. (0->0), fcn. (569->6), ass. (0->81)
t67 = sin(pkin(8));
t64 = t67 ^ 2;
t68 = cos(pkin(8));
t110 = t68 ^ 2 + t64;
t109 = qJ(3) * qJDD(2) + qJD(2) * qJD(3);
t66 = pkin(7) + qJ(2);
t62 = sin(t66);
t108 = g(1) * t62;
t63 = cos(t66);
t55 = g(2) * t63;
t84 = -t55 + t108;
t69 = sin(qJ(5));
t70 = cos(qJ(5));
t31 = t67 * t69 + t68 * t70;
t26 = t31 * qJD(2);
t101 = t68 * t69;
t32 = t67 * t70 - t101;
t87 = qJD(2) * t101;
t75 = qJD(5) * t87 - t31 * qJDD(2);
t107 = t68 * pkin(3);
t20 = t67 * qJDD(1) + t109 * t68;
t106 = t20 * t67 - g(3);
t93 = qJ(3) * qJD(2);
t104 = (t67 * qJD(1) + t68 * t93) * t68;
t16 = t20 * t68;
t105 = qJ(3) * t16 + qJD(3) * t104;
t103 = t62 * t68;
t100 = -pkin(6) + qJ(3);
t99 = t109 * t67;
t98 = t63 * pkin(2) + t62 * qJ(3);
t97 = qJD(2) * t67;
t96 = qJD(4) * t67;
t95 = qJDD(2) * pkin(2);
t94 = qJDD(1) - g(3);
t82 = t67 * qJ(4) + pkin(2);
t35 = -t82 - t107;
t92 = qJDD(2) * t35;
t59 = t67 * qJDD(2);
t91 = t68 * qJDD(2);
t88 = qJD(2) * qJD(4);
t86 = t70 * t97;
t85 = -t62 * pkin(2) + t63 * qJ(3);
t83 = t67 * t88;
t33 = t68 * qJD(1) - t67 * t93;
t58 = qJDD(3) - t95;
t81 = g(1) * t63 + g(2) * t62;
t80 = t70 * t59 - t69 * t91;
t19 = t68 * qJDD(1) - t99;
t37 = t100 * t67;
t38 = t100 * t68;
t79 = t70 * t37 - t69 * t38;
t78 = t69 * t37 + t70 * t38;
t77 = -t58 + t95 - t55;
t73 = qJDD(3) + t92;
t8 = t73 - t83;
t76 = -t8 - t92 - t55;
t17 = qJDD(4) - t19;
t24 = t31 * qJD(5);
t29 = (pkin(3) + pkin(4)) * t68 + t82;
t74 = t110 * t109 + t16 - t81;
t72 = qJD(2) ^ 2;
t71 = qJD(5) ^ 2;
t39 = g(1) * t103;
t36 = t110 * t72;
t30 = qJD(4) - t33;
t28 = t86 - t87;
t25 = t32 * qJD(5);
t22 = t35 * qJD(2) + qJD(3);
t14 = t31 * t63;
t13 = t32 * t63;
t12 = t31 * t62;
t11 = t32 * t62;
t10 = t29 * qJD(2) - qJD(3);
t7 = -pkin(6) * t91 + t20;
t6 = -pkin(6) * t59 + t17;
t5 = t29 * qJDD(2) - qJDD(3) + t83;
t4 = -t25 * qJD(5) - t31 * qJDD(5);
t3 = -t24 * qJD(5) + t32 * qJDD(5);
t2 = qJD(5) * t86 - t75;
t1 = -qJD(2) * t24 + t80;
t9 = [t94, 0, 0, 0, 0, 0, 0, t19 * t68 + t106, 0, 0, 0, -t17 * t68 + t106, 0, 0, 0, 0, 0, t4, -t3; 0, qJDD(2), t84, t81, t77 * t68 + t39, (-t77 - t108) * t67, -t19 * t67 + t74, -t58 * pkin(2) - g(1) * t85 - g(2) * t98 + (-t19 * qJ(3) - t33 * qJD(3)) * t67 + t105, t39 + (t76 + t83) * t68, t17 * t67 + t74, t64 * t88 + (t76 + t108) * t67, t8 * t35 - g(1) * (-pkin(3) * t103 + t85) - g(2) * (t63 * t107 + t98) + (t17 * qJ(3) + t84 * qJ(4) + t30 * qJD(3) - t22 * qJD(4)) * t67 + t105, t1 * t32 - t28 * t24, -t1 * t31 - t32 * t2 + t24 * t26 - t28 * t25, t3, t4, 0, t26 * t96 + t29 * t2 + t5 * t31 + t10 * t25 + t79 * qJDD(5) + g(1) * t12 - g(2) * t14 + (t32 * qJD(3) - t78 * qJD(5)) * qJD(5), t28 * t96 + t29 * t1 + t5 * t32 - t10 * t24 - t78 * qJDD(5) + g(1) * t11 - g(2) * t13 + (-t31 * qJD(3) - t79 * qJD(5)) * qJD(5); 0, 0, 0, 0, -t91, t59, -t36, (t33 * t67 - t104) * qJD(2) + t58 - t84, -t91, -t36, -t59, (-t104 + (-qJD(4) - t30) * t67) * qJD(2) + t73 - t84, 0, 0, 0, 0, 0, (-t28 - t86) * qJD(5) + t75, 0.2e1 * t26 * qJD(5) - t80; 0, 0, 0, 0, 0, 0, 0, 0, -t67 * t72 * t68, t59, -t64 * t72, qJDD(4) - t94 * t68 + (qJD(2) * t22 - t81) * t67 + t99, 0, 0, 0, 0, 0, qJDD(5) * t70 - t26 * t97 - t71 * t69, -qJDD(5) * t69 - t28 * t97 - t71 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t26, -t26 ^ 2 + t28 ^ 2, t80, (t28 - t86) * qJD(5) + t75, qJDD(5), -g(1) * t13 - g(2) * t11 + g(3) * t31 - t10 * t28 + t70 * t6 - t69 * t7, g(1) * t14 + g(2) * t12 + g(3) * t32 + t10 * t26 - t69 * t6 - t70 * t7;];
tau_reg = t9;
