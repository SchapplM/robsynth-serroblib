% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:42
% EndTime: 2019-12-31 17:42:43
% DurationCPUTime: 0.35s
% Computational Cost: add. (517->84), mult. (1183->130), div. (0->0), fcn. (872->8), ass. (0->73)
t48 = qJD(2) + qJD(3);
t54 = sin(qJ(3));
t58 = cos(qJ(2));
t55 = sin(qJ(2));
t57 = cos(qJ(3));
t88 = t57 * t55;
t35 = t54 * t58 + t88;
t64 = t35 * qJD(2);
t62 = (-qJD(3) * t88 - t64) * qJD(1);
t77 = t58 * qJD(1);
t41 = qJD(2) * pkin(2) + t77;
t78 = qJD(3) * t41;
t61 = -t54 * t78 + t62;
t59 = qJD(5) ^ 2;
t45 = t57 * pkin(2) + pkin(3);
t51 = sin(pkin(9));
t52 = cos(pkin(9));
t90 = t52 * t54;
t83 = pkin(2) * t90 + t51 * t45;
t32 = t35 * qJD(1);
t34 = -t54 * t55 + t57 * t58;
t33 = t34 * qJD(1);
t81 = pkin(2) * qJD(3);
t86 = t52 * t32 + t51 * t33 - (t51 * t57 + t90) * t81;
t96 = -t86 * t48 + (pkin(7) + t83) * t59;
t72 = qJD(2) * t77;
t79 = qJD(1) * t55;
t74 = t54 * t79;
t84 = t48 * t74;
t15 = (t72 + t78) * t57 - t84;
t2 = t51 * t15 - t52 * t61;
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t25 = t57 * t41 - t74;
t23 = t48 * pkin(3) + t25;
t26 = t54 * t41 + t57 * t79;
t93 = t51 * t26;
t10 = t52 * t23 - t93;
t8 = -t48 * pkin(4) - t10;
t80 = qJD(5) * t8;
t95 = t2 * t53 + t56 * t80;
t92 = t51 * t54;
t91 = t52 * t26;
t89 = t53 * t56;
t87 = t59 * t53;
t85 = t51 * t32 - t52 * t33 + (t52 * t57 - t92) * t81;
t82 = t53 ^ 2 - t56 ^ 2;
t76 = 0.2e1 * qJD(5) * t48;
t3 = t52 * t15 + t61 * t51;
t73 = -t8 * t48 - t3;
t19 = t51 * t34 + t52 * t35;
t20 = t48 * t34;
t21 = -t35 * qJD(3) - t64;
t4 = t51 * t20 - t52 * t21;
t71 = t19 * t59 + t4 * t48;
t12 = t51 * t25 + t91;
t70 = -t12 * t48 + (t51 * pkin(3) + pkin(7)) * t59;
t18 = -t52 * t34 + t51 * t35;
t5 = t52 * t20 + t51 * t21;
t69 = qJD(5) * (t18 * t48 - t5);
t68 = (-pkin(2) * t48 - t41) * qJD(3);
t13 = t52 * t25 - t93;
t67 = qJD(5) * ((-t52 * pkin(3) - pkin(4)) * t48 + t13);
t66 = -pkin(2) * t92 + t52 * t45;
t65 = qJD(5) * ((-pkin(4) - t66) * t48 - t85);
t60 = qJD(2) ^ 2;
t47 = t48 ^ 2;
t46 = t59 * t56;
t37 = t76 * t89;
t27 = t82 * t76;
t11 = t51 * t23 + t91;
t6 = t53 * t80;
t1 = [0, 0, -t60 * t55, -t60 * t58, 0, t21 * t48, -t20 * t48, -t10 * t4 + t11 * t5 + t2 * t18 + t3 * t19, 0, 0, 0, 0, 0, t53 * t69 - t71 * t56, t71 * t53 + t56 * t69; 0, 0, 0, 0, 0, t32 * t48 + t54 * t68 + t62, t33 * t48 + (t68 - t72) * t57 + t84, t86 * t10 + t85 * t11 - t2 * t66 + t3 * t83, t37, -t27, t46, -t87, 0, t6 + t53 * t65 + (-t2 - t96) * t56, t96 * t53 + t56 * t65 + t95; 0, 0, 0, 0, 0, t26 * t48 + t61, t25 * t48 - t15, t10 * t12 - t11 * t13 + (-t2 * t52 + t3 * t51) * pkin(3), t37, -t27, t46, -t87, 0, t6 + t53 * t67 + (-t2 - t70) * t56, t70 * t53 + t56 * t67 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t46; 0, 0, 0, 0, 0, 0, 0, 0, -t47 * t89, t82 * t47, 0, 0, 0, t73 * t53, t73 * t56;];
tauc_reg = t1;
