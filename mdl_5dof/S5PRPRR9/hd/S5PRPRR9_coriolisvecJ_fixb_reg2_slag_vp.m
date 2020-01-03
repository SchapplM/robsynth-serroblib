% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:45
% EndTime: 2019-12-31 17:39:47
% DurationCPUTime: 0.41s
% Computational Cost: add. (710->73), mult. (1188->112), div. (0->0), fcn. (533->4), ass. (0->67)
t32 = sin(qJ(4));
t36 = qJD(5) ^ 2;
t57 = qJD(2) - qJD(4);
t80 = t57 ^ 2;
t82 = t32 * (t36 + t80);
t35 = -pkin(2) - pkin(3);
t19 = t35 * qJD(2) + qJD(3);
t34 = cos(qJ(4));
t59 = qJD(2) * qJ(3);
t52 = t32 * t59;
t58 = qJD(2) * qJD(3);
t63 = qJD(4) * t34;
t39 = qJD(4) * t52 - t19 * t63 - t34 * t58;
t11 = t34 * t19 - t52;
t69 = t57 * t11;
t81 = t39 - t69;
t73 = t57 * pkin(4);
t7 = -t11 + t73;
t79 = t57 * t7;
t67 = t32 * t19;
t12 = t34 * t59 + t67;
t60 = t32 * qJD(3);
t6 = (qJ(3) * t63 + t60) * qJD(2) + qJD(4) * t67;
t78 = -t12 * t57 - t6;
t76 = qJD(5) * t57;
t66 = t34 * qJ(3) + t32 * t35;
t10 = t66 * qJD(4) + t60;
t75 = t10 * t57 + t6;
t74 = t34 * t80;
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t8 = -pkin(7) * t57 + t12;
t3 = -t33 * qJD(1) - t31 * t8;
t1 = t3 * qJD(5) - t33 * t39;
t61 = t31 * qJD(1);
t21 = qJD(5) * t61;
t62 = qJD(5) * t33;
t2 = t31 * t39 - t8 * t62 + t21;
t72 = t33 * t8;
t4 = -t61 + t72;
t38 = -(t3 * t33 + t31 * t4) * qJD(5) + t1 * t33 - t2 * t31;
t46 = -t32 * qJ(3) + t34 * t35;
t9 = t34 * qJD(3) + t46 * qJD(4);
t71 = t9 * t57;
t29 = t31 ^ 2;
t30 = t33 ^ 2;
t65 = t29 - t30;
t64 = t29 + t30;
t56 = t31 * t80 * t33;
t54 = 0.2e1 * t58;
t53 = t39 + t79;
t50 = t31 * t57 * t62;
t47 = t3 * t31 - t33 * t4;
t45 = pkin(7) * t36 - t78;
t16 = -pkin(7) + t66;
t44 = t16 * t36 - t75;
t43 = qJD(5) * (t11 + t7 + t73);
t15 = pkin(4) - t46;
t42 = qJD(5) * (-t15 * t57 - t7 - t9);
t41 = 0.2e1 * t34 * t76;
t37 = qJD(2) ^ 2;
t25 = t36 * t33;
t24 = t36 * t31;
t18 = 0.2e1 * t50;
t17 = -0.2e1 * t50;
t13 = t65 * t76;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t25, 0, t47 * qJD(5) - t1 * t31 - t2 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, qJ(3) * t54, 0, 0, 0, 0, 0, 0, t75, -t39 + t71, 0, -t11 * t10 + t12 * t9 - t39 * t66 - t6 * t46, t18, -0.2e1 * t13, -t25, t17, t24, 0, t31 * t42 - t44 * t33, t44 * t31 + t33 * t42, -t64 * t71 - t38, t7 * t10 + t6 * t15 + t38 * t16 - t47 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t37 * qJ(3), 0, 0, 0, 0, 0, 0, -t32 * t80, -t74, 0, -t81 * t32 + t78 * t34, 0, 0, 0, 0, 0, 0, t31 * t41 - t33 * t82, t31 * t82 + t33 * t41, t64 * t74, (t57 * t47 - t6) * t34 + (t38 - t79) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t81, 0, 0, t17, 0.2e1 * t13, t25, t18, -t24, 0, t31 * t43 - t45 * t33, t45 * t31 + t33 * t43, t64 * t69 + t38, -t6 * pkin(4) + t38 * pkin(7) + t47 * t11 - t7 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t65 * t80, 0, t56, 0, 0, t21 + t53 * t31 + (t4 - t72) * qJD(5), t53 * t33, 0, 0;];
tauc_reg = t5;
