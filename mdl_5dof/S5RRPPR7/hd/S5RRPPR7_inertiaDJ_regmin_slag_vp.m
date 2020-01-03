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
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 19:36:22
% EndTime: 2019-12-31 19:36:23
% DurationCPUTime: 0.41s
% Computational Cost: add. (523->89), mult. (1243->175), div. (0->0), fcn. (1078->6), ass. (0->65)
t38 = sin(qJ(5));
t34 = t38 ^ 2;
t40 = cos(qJ(5));
t65 = -t40 ^ 2 + t34;
t55 = t65 * qJD(5);
t71 = 2 * qJD(4);
t70 = pkin(3) + pkin(7);
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t39 = sin(qJ(2));
t41 = cos(qJ(2));
t23 = t36 * t39 - t37 * t41;
t24 = t36 * t41 + t37 * t39;
t20 = t24 * qJD(2);
t66 = -qJ(3) - pkin(6);
t56 = qJD(2) * t66;
t18 = t41 * qJD(3) + t39 * t56;
t19 = -t39 * qJD(3) + t41 * t56;
t9 = t37 * t18 + t36 * t19;
t5 = -t20 * pkin(4) + t9;
t69 = t5 * t23;
t68 = t38 * t20;
t67 = t40 * t20;
t64 = qJD(5) * t38;
t63 = qJD(5) * t40;
t62 = t39 * qJD(2);
t61 = t41 * qJD(2);
t60 = -0.2e1 * pkin(1) * qJD(2);
t59 = t38 * t67;
t33 = pkin(2) * t62;
t58 = t38 * t63;
t57 = -t41 * pkin(2) - pkin(1);
t32 = -t37 * pkin(2) - pkin(3);
t8 = t36 * t18 - t37 * t19;
t26 = t66 * t39;
t27 = t66 * t41;
t15 = -t37 * t26 - t36 * t27;
t16 = t36 * t26 - t37 * t27;
t54 = t15 * t8 + t16 * t9;
t10 = t24 * pkin(4) + t15;
t50 = -t24 * qJ(4) + t57;
t7 = t70 * t23 + t50;
t53 = t40 * t10 - t38 * t7;
t52 = t38 * t10 + t40 * t7;
t30 = t36 * pkin(2) + qJ(4);
t51 = -qJD(4) * t23 - t30 * t20;
t49 = t23 * t63 + t68;
t48 = t23 * t64 - t67;
t21 = -t36 * t62 + t37 * t61;
t47 = t38 * t21 + t24 * t63;
t46 = -t40 * t21 + t24 * t64;
t45 = -t21 * qJ(4) - t24 * qJD(4) + t33;
t29 = -pkin(7) + t32;
t44 = t5 + (t23 * t30 - t24 * t29) * qJD(5);
t11 = -t23 * pkin(4) + t16;
t43 = -qJD(5) * t11 - t21 * t29 - t51;
t42 = 0.2e1 * t15 * t21 - 0.2e1 * t16 * t20 - 0.2e1 * t9 * t23 + 0.2e1 * t8 * t24;
t22 = t23 ^ 2;
t14 = t23 * pkin(3) + t50;
t6 = t20 * pkin(3) + t45;
t4 = t21 * pkin(4) + t8;
t3 = t70 * t20 + t45;
t2 = -t52 * qJD(5) - t38 * t3 + t40 * t4;
t1 = -t53 * qJD(5) - t40 * t3 - t38 * t4;
t12 = [0, 0, 0, 0.2e1 * t39 * t61, 0.2e1 * (-t39 ^ 2 + t41 ^ 2) * qJD(2), 0, 0, 0, t39 * t60, t41 * t60, t42, 0.2e1 * t57 * t33 + 0.2e1 * t54, t42, -0.2e1 * t14 * t20 - 0.2e1 * t6 * t23, -0.2e1 * t14 * t21 - 0.2e1 * t6 * t24, 0.2e1 * t14 * t6 + 0.2e1 * t54, 0.2e1 * t34 * t23 * t20 + 0.2e1 * t22 * t58, -0.2e1 * t22 * t55 + 0.4e1 * t23 * t59, 0.2e1 * t47 * t23 + 0.2e1 * t24 * t68, -0.2e1 * t46 * t23 + 0.2e1 * t24 * t67, 0.2e1 * t24 * t21, 0.2e1 * t48 * t11 + 0.2e1 * t2 * t24 + 0.2e1 * t53 * t21 - 0.2e1 * t40 * t69, 0.2e1 * t1 * t24 + 0.2e1 * t49 * t11 - 0.2e1 * t52 * t21 + 0.2e1 * t38 * t69; 0, 0, 0, 0, 0, t61, -t62, 0, -pkin(6) * t61, pkin(6) * t62, (-t20 * t36 - t21 * t37) * pkin(2), (t36 * t9 - t37 * t8) * pkin(2), t32 * t21 + t51, t8, t9, t16 * qJD(4) + t9 * t30 + t8 * t32, -t23 * t55 + t59, -t65 * t20 - 0.4e1 * t23 * t58, -t46, -t47, 0, t44 * t38 - t43 * t40, t43 * t38 + t44 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t30 * t71, -0.2e1 * t58, 0.2e1 * t55, 0, 0, 0, 0.2e1 * qJD(4) * t38 + 0.2e1 * t30 * t63, 0.2e1 * qJD(4) * t40 - 0.2e1 * t30 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t20, -t21, t6, 0, 0, 0, 0, 0, -t47, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, t8, 0, 0, 0, 0, 0, -t46, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t48, t21, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63, 0, -t29 * t64, -t29 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
