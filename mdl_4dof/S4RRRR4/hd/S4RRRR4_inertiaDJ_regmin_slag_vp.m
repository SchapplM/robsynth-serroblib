% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:15
% EndTime: 2019-12-31 17:26:16
% DurationCPUTime: 0.39s
% Computational Cost: add. (479->85), mult. (1247->170), div. (0->0), fcn. (1069->6), ass. (0->69)
t42 = cos(qJ(4));
t38 = t42 ^ 2;
t39 = sin(qJ(4));
t70 = t39 ^ 2 - t38;
t55 = t70 * qJD(4);
t79 = qJD(2) + qJD(3);
t78 = pkin(5) + pkin(6);
t41 = sin(qJ(2));
t28 = t78 * t41;
t44 = cos(qJ(2));
t29 = t78 * t44;
t40 = sin(qJ(3));
t43 = cos(qJ(3));
t17 = -t40 * t28 + t43 * t29;
t57 = qJD(2) * t78;
t53 = t44 * t57;
t54 = t41 * t57;
t8 = t17 * qJD(3) - t40 * t54 + t43 * t53;
t6 = t8 * t39;
t16 = t43 * t28 + t40 * t29;
t36 = qJD(4) * t42;
t77 = t16 * t36 + t6;
t24 = t40 * t41 - t43 * t44;
t14 = t79 * t24;
t25 = t40 * t44 + t43 * t41;
t76 = t25 * t14;
t75 = t25 * t42;
t15 = t79 * t25;
t74 = t39 * t15;
t73 = t42 * t14;
t72 = t42 * t15;
t34 = -t43 * pkin(2) - pkin(3);
t69 = qJD(3) * t40;
t60 = pkin(2) * t69;
t71 = t34 * t36 + t39 * t60;
t68 = qJD(3) * t43;
t67 = qJD(4) * t39;
t66 = t41 * qJD(2);
t65 = t44 * qJD(2);
t64 = -0.2e1 * pkin(1) * qJD(2);
t63 = pkin(3) * t67;
t62 = pkin(3) * t36;
t61 = pkin(2) * t66;
t59 = pkin(2) * t68;
t58 = t39 * t36;
t35 = -t44 * pkin(2) - pkin(1);
t56 = -0.4e1 * t39 * t75;
t11 = t24 * pkin(3) - t25 * pkin(7) + t35;
t52 = t42 * t11 - t39 * t17;
t51 = t39 * t11 + t42 * t17;
t33 = t40 * pkin(2) + pkin(7);
t50 = t24 * t33 - t25 * t34;
t49 = t34 * t67 - t42 * t60;
t48 = -t39 * t14 + t25 * t36;
t47 = t25 * t67 + t73;
t46 = t24 * t67 - t72;
t45 = -t14 * t34 - t15 * t33 + (-t24 * t43 + t25 * t40) * qJD(3) * pkin(2);
t31 = 0.2e1 * t58;
t23 = -0.2e1 * t55;
t22 = t25 ^ 2;
t12 = t16 * t67;
t10 = t24 * t36 + t74;
t7 = t28 * t68 + t29 * t69 + t40 * t53 + t43 * t54;
t5 = t15 * pkin(3) + t14 * pkin(7) + t61;
t4 = -t25 * t55 - t39 * t73;
t3 = qJD(4) * t56 + t70 * t14;
t2 = -t51 * qJD(4) + t39 * t7 + t42 * t5;
t1 = -t52 * qJD(4) - t39 * t5 + t42 * t7;
t9 = [0, 0, 0, 0.2e1 * t41 * t65, 0.2e1 * (-t41 ^ 2 + t44 ^ 2) * qJD(2), 0, 0, 0, t41 * t64, t44 * t64, -0.2e1 * t76, 0.2e1 * t14 * t24 - 0.2e1 * t25 * t15, 0, 0, 0, 0.2e1 * t35 * t15 + 0.2e1 * t24 * t61, -0.2e1 * t35 * t14 + 0.2e1 * t25 * t61, -0.2e1 * t22 * t58 - 0.2e1 * t38 * t76, -t14 * t56 + 0.2e1 * t22 * t55, -0.2e1 * t47 * t24 + 0.2e1 * t25 * t72, -0.2e1 * t48 * t24 - 0.2e1 * t25 * t74, 0.2e1 * t24 * t15, 0.2e1 * t52 * t15 + 0.2e1 * t48 * t16 + 0.2e1 * t2 * t24 + 0.2e1 * t25 * t6, 0.2e1 * t1 * t24 - 0.2e1 * t51 * t15 - 0.2e1 * t47 * t16 + 0.2e1 * t8 * t75; 0, 0, 0, 0, 0, t65, -t66, 0, -pkin(5) * t65, pkin(5) * t66, 0, 0, -t14, -t15, 0, -t8, t7, t4, t3, t10, -t46, 0, t12 + (-t50 * qJD(4) - t8) * t42 + t45 * t39, t45 * t42 + t50 * t67 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t60, -0.2e1 * t59, t31, t23, 0, 0, 0, 0.2e1 * t49, 0.2e1 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, -t8, t7, t4, t3, t10, -t46, 0, t12 + (pkin(3) * t14 - pkin(7) * t15) * t39 + (-t8 + (-pkin(3) * t25 - pkin(7) * t24) * qJD(4)) * t42, t47 * pkin(3) + t46 * pkin(7) + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t59, t31, t23, 0, 0, 0, t49 - t63, -t62 + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t23, 0, 0, 0, -0.2e1 * t63, -0.2e1 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t48, t15, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t67, 0, -t33 * t36 - t39 * t59, t33 * t67 - t42 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t67, 0, -pkin(7) * t36, pkin(7) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
