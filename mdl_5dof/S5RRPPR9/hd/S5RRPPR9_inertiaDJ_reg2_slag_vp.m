% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:48
% EndTime: 2019-12-31 19:41:51
% DurationCPUTime: 0.85s
% Computational Cost: add. (572->98), mult. (1181->189), div. (0->0), fcn. (811->4), ass. (0->73)
t36 = qJ(3) + pkin(4);
t38 = sin(qJ(2));
t40 = cos(qJ(2));
t78 = pkin(2) + pkin(3);
t63 = pkin(7) + t78;
t53 = t63 * t40;
t85 = t36 * t38 + t53;
t84 = t78 * t40;
t83 = t36 * t40 - t63 * t38;
t82 = pkin(6) - qJ(4);
t66 = t38 * qJD(3);
t81 = t83 * qJD(2) + t66;
t71 = t38 * qJ(3);
t55 = pkin(1) + t71;
t80 = t38 * pkin(4) + t53 + t55;
t37 = sin(qJ(5));
t32 = t37 ^ 2;
t39 = cos(qJ(5));
t34 = t39 ^ 2;
t74 = t32 - t34;
t54 = qJD(5) * t74;
t26 = t38 * qJD(2);
t61 = pkin(6) * t26;
t64 = qJ(4) * qJD(2);
t11 = t40 * qJD(4) - t38 * t64 + t61;
t79 = t83 * qJD(5) + t11;
t42 = 0.2e1 * qJD(3);
t21 = t82 * t40;
t77 = t21 * t11;
t35 = t40 ^ 2;
t73 = t38 ^ 2 - t35;
t72 = qJ(3) * t40;
t70 = qJD(2) * t39;
t28 = qJD(5) * t37;
t69 = qJD(5) * t39;
t68 = qJD(5) * t40;
t67 = t21 * qJD(3);
t27 = t40 * qJD(2);
t65 = t40 * qJD(3);
t62 = -0.2e1 * pkin(1) * qJD(2);
t60 = t37 * t68;
t59 = t39 * t68;
t58 = t37 * t69;
t57 = t38 * t27;
t56 = t39 * t26;
t19 = t73 * qJD(2);
t52 = t37 * t56;
t51 = t35 * t58;
t20 = t82 * t38;
t4 = -t37 * t20 + t80 * t39;
t5 = t39 * t20 + t80 * t37;
t50 = -t37 * t5 - t39 * t4;
t49 = -t40 * pkin(2) - t71;
t46 = t49 * qJD(2) + t65;
t25 = pkin(6) * t27;
t12 = -t38 * qJD(4) - t40 * t64 + t25;
t2 = t20 * t28 - t39 * t12 - t37 * t81 - (pkin(1) + t85) * t69;
t3 = -qJD(5) * t5 - t37 * t12 + t81 * t39;
t44 = -t2 * t37 + t3 * t39 + (-t37 * t4 + t39 * t5) * qJD(5);
t1 = t50 * qJD(5) - t2 * t39 - t3 * t37;
t43 = t85 * qJD(2) - qJD(5) * t21 - t65;
t30 = qJ(3) * t42;
t23 = -0.2e1 * t57;
t22 = 0.2e1 * t57;
t18 = -pkin(1) + t49;
t17 = -0.2e1 * t19;
t16 = t55 + t84;
t15 = -t37 * t27 - t38 * t69;
t14 = -t39 * t27 + t38 * t28;
t13 = -t66 + (pkin(2) * t38 - t72) * qJD(2);
t9 = t66 + (-t78 * t38 + t72) * qJD(2);
t8 = t40 * t54 + t52;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t17, 0, t23, 0, 0, t38 * t62, t40 * t62, 0, 0, t22, 0, 0.2e1 * t19, 0, 0, t23, -0.2e1 * t13 * t40 + 0.2e1 * t18 * t26, 0, -0.2e1 * t13 * t38 - 0.2e1 * t18 * t27, 0.2e1 * t18 * t13, t23, t17, 0, t22, 0, 0, 0.2e1 * t16 * t27 + 0.2e1 * t9 * t38, 0.2e1 * t16 * t26 - 0.2e1 * t9 * t40, 0.2e1 * t11 * t40 - 0.2e1 * t12 * t38 + 0.2e1 * (-t20 * t40 + t21 * t38) * qJD(2), 0.2e1 * t20 * t12 + 0.2e1 * t16 * t9 - 0.2e1 * t77, -0.2e1 * t34 * t57 - 0.2e1 * t51, 0.2e1 * t35 * t54 + 0.4e1 * t40 * t52, 0.2e1 * t38 * t60 + 0.2e1 * t73 * t70, -0.2e1 * t32 * t57 + 0.2e1 * t51, -0.2e1 * t37 * t19 + 0.2e1 * t38 * t59, t22, 0.2e1 * (qJD(2) * t21 * t37 + t3) * t38 + 0.2e1 * (qJD(2) * t4 + t11 * t37 - t21 * t69) * t40, 0.2e1 * (t21 * t70 + t2) * t38 + 0.2e1 * (-qJD(2) * t5 + t11 * t39 + t21 * t28) * t40, 0.2e1 * t50 * t26 + 0.2e1 * t44 * t40, -0.2e1 * t5 * t2 + 0.2e1 * t4 * t3 - 0.2e1 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t26, 0, -t25, t61, 0, 0, 0, t27, 0, 0, t26, 0, -t25, t46, -t61, t46 * pkin(6), 0, 0, -t26, 0, t27, 0, -t11, t12, -t65 + (t71 + t84) * qJD(2), -t11 * qJ(3) - t12 * t78 + t67, -t8, t74 * t26 - 0.4e1 * t40 * t58, t15, t8, t14, 0, t43 * t37 - t79 * t39, t79 * t37 + t43 * t39, -t1, -t1 * t63 - t11 * t36 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t30, 0, 0, 0, 0, 0, 0, t42, 0, 0, t30, 0.2e1 * t58, -0.2e1 * t54, 0, -0.2e1 * t58, 0, 0, 0.2e1 * qJD(3) * t39 - 0.2e1 * t36 * t28, -0.2e1 * qJD(3) * t37 - 0.2e1 * t36 * t69, 0, t36 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t12, 0, 0, 0, 0, 0, 0, t15, t14, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t26, 0, t9, 0, 0, 0, 0, 0, 0, -t14, t15, (-t32 - t34) * t26, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 + t60, 0, -t37 * t26 + t59, t27, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, t28, 0, t63 * t69, -t63 * t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t69, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
