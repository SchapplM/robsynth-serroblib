% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:50
% EndTime: 2019-12-31 17:47:53
% DurationCPUTime: 0.61s
% Computational Cost: add. (437->123), mult. (1209->223), div. (0->0), fcn. (821->6), ass. (0->90)
t108 = pkin(3) + qJ(2);
t56 = cos(pkin(8));
t57 = cos(pkin(7));
t93 = qJD(1) * t57;
t35 = t56 * t93 + qJD(5);
t107 = t35 * t56;
t55 = sin(pkin(7));
t51 = t55 ^ 2;
t60 = qJD(1) ^ 2;
t106 = t51 * t60;
t54 = sin(pkin(8));
t105 = t54 * t60;
t104 = t55 * t56;
t103 = t55 * t57;
t58 = sin(qJ(5));
t102 = t55 * t58;
t59 = cos(qJ(5));
t101 = t55 * t59;
t100 = t56 * t57;
t99 = t57 * t58;
t98 = t57 * t59;
t75 = -t55 * qJ(3) - pkin(1);
t27 = (-pkin(2) - qJ(4)) * t57 + t75;
t12 = t27 * qJD(1) + qJD(2);
t87 = qJ(2) * qJD(1);
t41 = t55 * t87 + qJD(3);
t94 = qJD(1) * t55;
t28 = pkin(3) * t94 + t41;
t4 = t56 * t12 + t54 * t28;
t36 = t108 * t55;
t97 = t56 * t27 + t54 * t36;
t82 = t54 * t93;
t73 = qJD(5) * t82;
t78 = qJD(5) * t94;
t15 = t58 * t73 + t59 * t78;
t96 = t108 * t57;
t53 = t57 ^ 2;
t95 = t51 + t53;
t92 = qJD(2) * t53;
t91 = qJD(2) * t55;
t90 = qJD(2) * t57;
t89 = qJD(5) * t35;
t88 = t55 * qJD(3);
t86 = qJD(1) * qJD(2);
t85 = qJD(1) * qJD(3);
t84 = t60 * t103;
t83 = 0.2e1 * t92;
t81 = t58 * t89;
t80 = t59 * t89;
t29 = pkin(3) * t93 + t57 * t87 + qJD(4);
t79 = (-t54 ^ 2 - t56 ^ 2) * t60;
t33 = t95 * t60;
t77 = t55 * t86;
t76 = t57 * t86;
t74 = qJ(2) * t86;
t2 = pkin(6) * t94 + t4;
t65 = (pkin(4) * t56 + pkin(6) * t54) * t57;
t7 = qJD(1) * t65 + t29;
t72 = -t59 * t2 - t58 * t7;
t71 = t58 * t2 - t59 * t7;
t3 = -t54 * t12 + t56 * t28;
t70 = -t3 * t54 + t4 * t56;
t34 = -t57 * qJD(4) - t88;
t31 = t34 * qJD(1);
t8 = t54 * t31 - t56 * t77;
t9 = t56 * t31 + t54 * t77;
t69 = -t8 * t54 - t9 * t56;
t68 = -t54 * t27 + t56 * t36;
t67 = -t57 * pkin(2) + t75;
t66 = t54 * t98 - t102;
t24 = t54 * t99 + t101;
t16 = -t58 * t78 + t59 * t73;
t64 = t58 * t76 + t59 * t9;
t18 = t24 * qJD(1);
t63 = qJD(1) * (-t58 * t107 - t18 * t54);
t21 = -t58 * t94 + t59 * t82;
t62 = qJD(1) * (-t59 * t107 - t21 * t54);
t61 = t72 * qJD(5) - t58 * t9 + t59 * t76;
t43 = t51 * t74;
t26 = 0.2e1 * t95 * t86;
t23 = t67 * qJD(1) + qJD(2);
t20 = t66 * qJD(5);
t19 = t24 * qJD(5);
t14 = t56 * t34 + t54 * t91;
t13 = t54 * t34 - t56 * t91;
t11 = t65 + t96;
t6 = t55 * pkin(6) + t97;
t5 = -t55 * pkin(4) - t68;
t1 = -pkin(4) * t94 - t3;
t10 = [0, 0, 0, 0, 0, t26, 0.2e1 * t53 * t74 + 0.2e1 * t43, t26, -0.2e1 * t85 * t103, 0.2e1 * t51 * t85, t43 + (t41 * qJD(2) - t23 * qJD(3)) * t55 + (qJ(2) * t83 - t67 * t88) * qJD(1), -t8 * t55 + (-t13 * t55 + t56 * t83) * qJD(1), -t9 * t55 + (-t14 * t55 - 0.2e1 * t54 * t92) * qJD(1), ((-t13 * t54 - t14 * t56) * qJD(1) + t69) * t57, t9 * t97 + t4 * t14 - t8 * t68 - t3 * t13 + (qJD(1) * t96 + t29) * t90, -t15 * t66 - t21 * t19, t15 * t24 - t16 * t66 + t19 * t18 - t21 * t20, t15 * t100 + t19 * t35, t16 * t100 + t20 * t35, 0, (t59 * t90 - t58 * t14 + (-t11 * t58 - t59 * t6) * qJD(5)) * t35 + t61 * t100 - t13 * t18 - t5 * t16 - t8 * t24 - t1 * t20, -(t58 * t90 + t59 * t14 + (t11 * t59 - t58 * t6) * qJD(5)) * t35 - (-qJD(5) * t71 + t64) * t100 - t13 * t21 + t5 * t15 - t8 * t66 + t1 * t19; 0, 0, 0, 0, 0, -t33, -qJ(2) * t33, -t33, 0, 0, -t53 * t60 * qJ(2) + (-qJD(3) - t41) * t94, -t56 * t33, t95 * t105, 0, (-t29 * t57 + (-t3 * t56 - t4 * t54) * t55) * qJD(1) - t69, 0, 0, 0, 0, 0, -t56 * t80 - t54 * t16 + (-(-t54 * t102 + t98) * t35 - t18 * t104) * qJD(1), t56 * t81 + t54 * t15 + ((t54 * t101 + t99) * t35 - t21 * t104) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, t84, -t106, (qJD(2) + t23) * t94, -t51 * t105, -t56 * t106, t79 * t103, t9 * t54 - t8 * t56 + t70 * t94, 0, 0, 0, 0, 0, t56 * t16 - t54 * t80 + t55 * t63, -t56 * t15 + t54 * t81 + t55 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54 * t84, -t56 * t84, t53 * t79, (qJD(2) + t70) * t93, 0, 0, 0, 0, 0, t57 * t63 - t81, t57 * t62 - t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t18, -t18 ^ 2 + t21 ^ 2, -t18 * t35 + t15, -t21 * t35 + t16, 0, t1 * t21 - t35 * t72 + t61, -t1 * t18 - t64 + (qJD(5) - t35) * t71;];
tauc_reg = t10;
