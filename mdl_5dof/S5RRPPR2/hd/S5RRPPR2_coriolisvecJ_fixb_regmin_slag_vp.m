% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:36
% EndTime: 2020-01-03 11:57:37
% DurationCPUTime: 0.43s
% Computational Cost: add. (598->101), mult. (1320->166), div. (0->0), fcn. (762->8), ass. (0->88)
t52 = sin(pkin(8));
t56 = sin(qJ(2));
t86 = pkin(1) * qJD(1);
t79 = t56 * t86;
t37 = t52 * t79;
t54 = cos(pkin(8));
t58 = cos(qJ(2));
t78 = t58 * t86;
t28 = t54 * t78 - t37;
t23 = t28 * qJD(2);
t48 = qJD(1) + qJD(2);
t16 = t48 * qJD(4) + t23;
t51 = sin(pkin(9));
t46 = t51 ^ 2;
t53 = cos(pkin(9));
t88 = t53 ^ 2 + t46;
t106 = t16 * t88;
t105 = qJD(4) - t28;
t85 = pkin(1) * qJD(2);
t92 = t54 * t56;
t27 = (t52 * t58 + t92) * t85;
t22 = qJD(1) * t27;
t61 = -t53 * pkin(4) - t51 * pkin(7) - pkin(3);
t33 = t48 * pkin(2) + t78;
t18 = t54 * t33 - t37;
t66 = qJD(4) - t18;
t3 = t61 * t48 + t66;
t55 = sin(qJ(5));
t57 = cos(qJ(5));
t38 = t54 * t79;
t19 = t52 * t33 + t38;
t14 = t48 * qJ(4) + t19;
t6 = t51 * qJD(3) + t53 * t14;
t93 = t53 * t57;
t97 = t46 * t57;
t104 = (t16 * t93 + t55 * t22 + (t57 * t3 - t55 * t6) * qJD(5)) * t53 + t16 * t97;
t5 = -t53 * qJD(3) + t51 * t14;
t103 = t5 * t51;
t102 = t54 * pkin(2);
t40 = t52 * t56 * pkin(1);
t62 = t54 * t58 * t85 - qJD(2) * t40;
t24 = qJD(4) + t62;
t101 = t24 * t48;
t45 = t48 ^ 2;
t100 = t46 * t45;
t99 = t46 * t48;
t98 = t46 * t55;
t96 = t48 * t51;
t95 = t53 * t48;
t94 = t53 * t55;
t36 = -qJD(5) + t95;
t91 = t57 * t36;
t43 = t58 * pkin(1) + pkin(2);
t89 = pkin(1) * t92 + t52 * t43;
t87 = t55 ^ 2 - t57 ^ 2;
t84 = qJD(5) * t55;
t83 = qJD(5) * t57;
t82 = qJD(5) + t36;
t81 = t5 * t96;
t80 = t48 * t97;
t77 = qJD(5) * t99;
t76 = t51 * t84;
t75 = t51 * t83;
t73 = t54 * t43 - t40;
t68 = -t55 * t3 - t57 * t6;
t2 = t68 * qJD(5) - t16 * t94 + t57 * t22;
t71 = t16 * t98 - t2 * t53 + t5 * t75;
t70 = t36 * t76;
t69 = t82 * t96;
t67 = t6 * t53 + t103;
t65 = (-qJD(2) + t48) * t86;
t64 = (-qJD(1) - t48) * t85;
t63 = t36 * t53 + t99;
t42 = t52 * pkin(2) + qJ(4);
t60 = t105 * t55 + t42 * t83;
t59 = -t36 ^ 2 - t100;
t32 = t76 * t95;
t31 = -0.2e1 * t57 * t55 * t77;
t29 = t61 - t102;
t26 = t52 * t78 + t38;
t25 = qJ(4) + t89;
t20 = t22 * t51;
t17 = 0.2e1 * t87 * t77;
t15 = t61 - t73;
t11 = -t48 * pkin(3) + t66;
t10 = (t36 + t95) * t75;
t9 = t32 + t70;
t1 = [0, 0, 0, 0, t56 * t64, t58 * t64, -t18 * t27 + t19 * t62 - t22 * t73 + t23 * t89, (-t27 * t48 - t22) * t53, t27 * t96 + t20, t88 * t101 + t106, t22 * (-pkin(3) - t73) + t11 * t27 + t67 * t24 + t25 * t106, t31, t17, t9, t10, 0, -(-t24 * t94 + t57 * t27) * t36 + t98 * t101 + (-(-t15 * t55 - t25 * t93) * t36 + t25 * t80) * qJD(5) + t71, (t24 * t93 + t55 * t27) * t36 + t24 * t80 + (t15 * t91 + (-t63 * t25 - t103) * t55) * qJD(5) + t104; 0, 0, 0, 0, t56 * t65, t58 * t65, t18 * t26 - t19 * t28 + (-t22 * t54 + t23 * t52) * pkin(2), (t26 * t48 - t22) * t53, -t26 * t96 + t20, t105 * t48 * t88 + t106, t22 * (-pkin(3) - t102) - t11 * t26 + t42 * t106 + t105 * t67, t31, t17, t9, t10, 0, (t57 * t26 + t29 * t84 + t60 * t53) * t36 + t60 * t99 + t71, -t55 * t26 * t36 + t105 * t63 * t57 + (t29 * t91 + (-t63 * t42 - t103) * t55) * qJD(5) + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t36 - t95) * t75, t32 - t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88 * t45, -t67 * t48 + t22, 0, 0, 0, 0, 0, t59 * t55, t59 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t55 * t97, -t87 * t100, -t55 * t69, -t57 * t69, 0, t68 * t36 - t57 * t81 + t2, (-t53 * t16 - t82 * t3) * t57 + (t82 * t6 - t22 + t81) * t55;];
tauc_reg = t1;
