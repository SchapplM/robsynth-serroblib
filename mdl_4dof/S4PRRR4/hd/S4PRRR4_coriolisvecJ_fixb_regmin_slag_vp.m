% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:39
% EndTime: 2019-12-31 16:32:40
% DurationCPUTime: 0.29s
% Computational Cost: add. (240->65), mult. (655->111), div. (0->0), fcn. (426->4), ass. (0->58)
t33 = qJD(3) + qJD(4);
t68 = qJD(4) - t33;
t55 = (qJD(2) * qJD(3));
t70 = -2 * t55;
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t67 = pkin(5) + pkin(6);
t49 = qJD(2) * t67;
t14 = t39 * qJD(1) - t37 * t49;
t36 = sin(qJ(4));
t38 = cos(qJ(4));
t20 = t36 * t39 + t38 * t37;
t69 = qJD(2) * t20;
t19 = t36 * t37 - t38 * t39;
t7 = t33 * t19;
t66 = t7 * t33;
t57 = qJD(2) * t39;
t51 = t38 * t57;
t58 = qJD(2) * t37;
t52 = t36 * t58;
t16 = -t51 + t52;
t18 = -t36 * t57 - t38 * t58;
t65 = t18 * t16;
t25 = t67 * t39;
t56 = t37 * qJD(1);
t15 = qJD(2) * t25 + t56;
t64 = t38 * t15;
t41 = qJD(2) ^ 2;
t63 = t39 * t41;
t40 = qJD(3) ^ 2;
t62 = t40 * t37;
t61 = t40 * t39;
t60 = t37 ^ 2 - t39 ^ 2;
t59 = qJD(3) * pkin(3);
t54 = pkin(3) * t58;
t53 = t37 * t59;
t31 = -t39 * pkin(3) - pkin(2);
t13 = t14 + t59;
t50 = -pkin(3) * t33 - t13;
t48 = qJD(3) * t67;
t47 = t39 * t55;
t46 = pkin(2) * t70;
t11 = t14 * qJD(3);
t12 = (-t39 * t49 - t56) * qJD(3);
t23 = t31 * qJD(2);
t43 = -t36 * t11 + t38 * t12 + t23 * t18;
t4 = qJD(4) * t51 - t33 * t52 + t38 * t47;
t42 = t23 * t16 + (t68 * t15 - t12) * t36;
t8 = t33 * t20;
t24 = t67 * t37;
t22 = t39 * t48;
t21 = t37 * t48;
t6 = t8 * t33;
t5 = t8 * qJD(2);
t3 = -t16 ^ 2 + t18 ^ 2;
t2 = (-t18 - t69) * t33;
t1 = t16 * t33 + t4;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t61, 0, 0, 0, 0, 0, -t6, t66; 0, 0, 0, 0, 0.2e1 * t37 * t47, t60 * t70, t61, -t62, 0, -pkin(5) * t61 + t37 * t46, pkin(5) * t62 + t39 * t46, t18 * t7 + t4 * t20, t7 * t16 + t18 * t8 - t4 * t19 - t20 * t5, -t66, -t6, 0, t31 * t5 + t23 * t8 + (t36 * t21 - t38 * t22 + (t24 * t36 - t25 * t38) * qJD(4)) * t33 + (qJD(2) * t19 + t16) * t53, t31 * t4 - t23 * t7 - (-t38 * t21 - t36 * t22 + (-t24 * t38 - t25 * t36) * qJD(4)) * t33 + (-t18 + t69) * t53; 0, 0, 0, 0, -t37 * t63, t60 * t41, 0, 0, 0, t41 * pkin(2) * t37, pkin(2) * t63, -t65, t3, t1, t2, 0, -t16 * t54 - (-t36 * t14 - t64) * t33 + (t50 * t36 - t64) * qJD(4) + t43, t18 * t54 + (t50 * qJD(4) + t14 * t33 - t11) * t38 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t3, t1, t2, 0, t43 + t68 * (-t36 * t13 - t64), (-t68 * t13 - t11) * t38 + t42;];
tauc_reg = t9;
