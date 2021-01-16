% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:59:01
% EndTime: 2021-01-15 19:59:03
% DurationCPUTime: 0.45s
% Computational Cost: add. (543->93), mult. (1301->183), div. (0->0), fcn. (1122->6), ass. (0->65)
t38 = sin(qJ(5));
t34 = t38 ^ 2;
t40 = cos(qJ(5));
t65 = -t40 ^ 2 + t34;
t56 = t65 * qJD(5);
t71 = 2 * qJD(4);
t70 = pkin(3) + pkin(7);
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t39 = sin(qJ(2));
t41 = cos(qJ(2));
t22 = t36 * t39 - t37 * t41;
t23 = t36 * t41 + t37 * t39;
t19 = t23 * qJD(2);
t66 = -qJ(3) - pkin(6);
t57 = qJD(2) * t66;
t18 = t41 * qJD(3) + t39 * t57;
t45 = -t39 * qJD(3) + t41 * t57;
t9 = t37 * t18 + t36 * t45;
t5 = -t19 * pkin(4) + t9;
t69 = t5 * t22;
t68 = t38 * t19;
t67 = t40 * t19;
t64 = qJD(5) * t38;
t63 = qJD(5) * t40;
t62 = t39 * qJD(2);
t61 = t41 * qJD(2);
t60 = -0.2e1 * pkin(1) * qJD(2);
t59 = t38 * t67;
t33 = pkin(2) * t62;
t58 = t38 * t63;
t32 = -t41 * pkin(2) - pkin(1);
t31 = -t37 * pkin(2) - pkin(3);
t8 = t36 * t18 - t37 * t45;
t25 = t66 * t39;
t26 = t66 * t41;
t15 = -t37 * t25 - t36 * t26;
t16 = t36 * t25 - t37 * t26;
t55 = t15 * t8 + t16 * t9;
t10 = t23 * pkin(4) + t15;
t51 = -t23 * qJ(4) + t32;
t7 = t70 * t22 + t51;
t54 = t40 * t10 - t38 * t7;
t53 = t38 * t10 + t40 * t7;
t29 = t36 * pkin(2) + qJ(4);
t52 = -qJD(4) * t22 - t29 * t19;
t50 = t22 * t63 + t68;
t49 = t22 * t64 - t67;
t20 = -t36 * t62 + t37 * t61;
t48 = t38 * t20 + t23 * t63;
t47 = -t40 * t20 + t23 * t64;
t46 = -t20 * qJ(4) - t23 * qJD(4) + t33;
t28 = -pkin(7) + t31;
t44 = t5 + (t22 * t29 - t23 * t28) * qJD(5);
t11 = -t22 * pkin(4) + t16;
t43 = -qJD(5) * t11 - t20 * t28 - t52;
t42 = 0.2e1 * t15 * t20 - 0.2e1 * t16 * t19 - 0.2e1 * t9 * t22 + 0.2e1 * t8 * t23;
t21 = t22 ^ 2;
t14 = t22 * pkin(3) + t51;
t6 = t19 * pkin(3) + t46;
t4 = t20 * pkin(4) + t8;
t3 = t70 * t19 + t46;
t2 = -t53 * qJD(5) - t38 * t3 + t40 * t4;
t1 = -t54 * qJD(5) - t40 * t3 - t38 * t4;
t12 = [0, 0, 0, 0.2e1 * t39 * t61, 0.2e1 * (-t39 ^ 2 + t41 ^ 2) * qJD(2), 0, 0, 0, t39 * t60, t41 * t60, 0.2e1 * t32 * t19 + 0.2e1 * t22 * t33, 0.2e1 * t32 * t20 + 0.2e1 * t23 * t33, t42, 0.2e1 * t32 * t33 + 0.2e1 * t55, t42, -0.2e1 * t14 * t19 - 0.2e1 * t6 * t22, -0.2e1 * t14 * t20 - 0.2e1 * t6 * t23, 0.2e1 * t14 * t6 + 0.2e1 * t55, 0.2e1 * t34 * t22 * t19 + 0.2e1 * t21 * t58, -0.2e1 * t21 * t56 + 0.4e1 * t22 * t59, 0.2e1 * t48 * t22 + 0.2e1 * t23 * t68, -0.2e1 * t47 * t22 + 0.2e1 * t23 * t67, 0.2e1 * t23 * t20, 0.2e1 * t49 * t11 + 0.2e1 * t2 * t23 + 0.2e1 * t54 * t20 - 0.2e1 * t40 * t69, 0.2e1 * t1 * t23 + 0.2e1 * t50 * t11 - 0.2e1 * t53 * t20 + 0.2e1 * t38 * t69; 0, 0, 0, 0, 0, t61, -t62, 0, -pkin(6) * t61, pkin(6) * t62, -t8, -t9, (-t19 * t36 - t20 * t37) * pkin(2), (t36 * t9 - t37 * t8) * pkin(2), t31 * t20 + t52, t8, t9, t16 * qJD(4) + t9 * t29 + t8 * t31, -t22 * t56 + t59, -t65 * t19 - 0.4e1 * t22 * t58, -t47, -t48, 0, t44 * t38 - t43 * t40, t43 * t38 + t44 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t29 * t71, -0.2e1 * t58, 0.2e1 * t56, 0, 0, 0, 0.2e1 * qJD(4) * t38 + 0.2e1 * t29 * t63, 0.2e1 * qJD(4) * t40 - 0.2e1 * t29 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, t33, 0, -t19, -t20, t6, 0, 0, 0, 0, 0, -t48, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, t8, 0, 0, 0, 0, 0, -t47, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t49, t20, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63, 0, -t28 * t64, -t28 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
