% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR13_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:52
% EndTime: 2019-12-31 18:32:54
% DurationCPUTime: 0.39s
% Computational Cost: add. (466->85), mult. (1138->164), div. (0->0), fcn. (1048->6), ass. (0->61)
t37 = sin(qJ(5));
t33 = t37 ^ 2;
t39 = cos(qJ(5));
t64 = -t39 ^ 2 + t33;
t55 = t64 * qJD(5);
t35 = sin(pkin(8));
t66 = pkin(6) + qJ(2);
t25 = t66 * t35;
t36 = cos(pkin(8));
t26 = t66 * t36;
t38 = sin(qJ(3));
t70 = cos(qJ(3));
t56 = qJD(3) * t70;
t57 = t70 * t36;
t8 = (qJD(2) * t35 + qJD(3) * t26) * t38 - qJD(2) * t57 + t25 * t56;
t68 = t38 * t35;
t18 = qJD(3) * t68 - t36 * t56;
t74 = -0.2e1 * t18;
t73 = 2 * qJD(4);
t72 = pkin(3) + pkin(7);
t23 = -t57 + t68;
t24 = t35 * t70 + t36 * t38;
t19 = t24 * qJD(3);
t4 = -t19 * pkin(4) - t8;
t71 = t4 * t23;
t69 = t37 * t19;
t67 = t39 * t19;
t63 = qJD(5) * t37;
t62 = qJD(5) * t39;
t61 = qJD(5) * t72;
t60 = qJ(4) * qJD(5);
t59 = t37 * t67;
t58 = t37 * t62;
t30 = -t36 * pkin(2) - pkin(1);
t54 = 0.2e1 * (t35 ^ 2 + t36 ^ 2) * qJD(2);
t16 = t25 * t70 + t26 * t38;
t10 = pkin(4) * t24 + t16;
t48 = -t24 * qJ(4) + t30;
t7 = t23 * t72 + t48;
t53 = t10 * t39 - t37 * t7;
t52 = t10 * t37 + t39 * t7;
t51 = t18 * qJ(4) - t24 * qJD(4);
t50 = -qJ(4) * t19 - qJD(4) * t23;
t47 = -t18 * t37 + t24 * t62;
t46 = t18 * t39 + t24 * t63;
t45 = t23 * t62 + t69;
t44 = t23 * t63 - t67;
t43 = -t25 * t38 + t26 * t70;
t42 = t4 + (qJ(4) * t23 + t24 * t72) * qJD(5);
t11 = -pkin(4) * t23 + t43;
t41 = -qJD(5) * t11 - t18 * t72 - t50;
t9 = qJD(2) * t24 + qJD(3) * t43;
t21 = t23 ^ 2;
t15 = t24 * t74;
t14 = pkin(3) * t23 + t48;
t6 = pkin(3) * t19 + t51;
t5 = -t18 * pkin(4) + t9;
t3 = t19 * t72 + t51;
t2 = -qJD(5) * t52 - t37 * t3 + t39 * t5;
t1 = -qJD(5) * t53 - t39 * t3 - t37 * t5;
t12 = [0, 0, 0, 0, 0, t54, qJ(2) * t54, t15, 0.2e1 * t18 * t23 - 0.2e1 * t19 * t24, 0, 0, 0, 0.2e1 * t30 * t19, t30 * t74, -0.2e1 * t16 * t18 - 0.2e1 * t19 * t43 + 0.2e1 * t23 * t8 + 0.2e1 * t24 * t9, -0.2e1 * t14 * t19 - 0.2e1 * t23 * t6, 0.2e1 * t14 * t18 - 0.2e1 * t24 * t6, 0.2e1 * t14 * t6 + 0.2e1 * t16 * t9 - 0.2e1 * t43 * t8, 0.2e1 * t19 * t23 * t33 + 0.2e1 * t21 * t58, -0.2e1 * t21 * t55 + 0.4e1 * t23 * t59, 0.2e1 * t23 * t47 + 0.2e1 * t24 * t69, -0.2e1 * t23 * t46 + 0.2e1 * t24 * t67, t15, 0.2e1 * t11 * t44 - 0.2e1 * t18 * t53 + 0.2e1 * t2 * t24 - 0.2e1 * t39 * t71, 0.2e1 * t1 * t24 + 0.2e1 * t11 * t45 + 0.2e1 * t18 * t52 + 0.2e1 * t37 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, -t19, t18, t6, 0, 0, 0, 0, 0, -t47, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, -t9, t8, pkin(3) * t18 + t50, t9, -t8, -pkin(3) * t9 - qJ(4) * t8 + qJD(4) * t43, -t23 * t55 + t59, -t19 * t64 - 0.4e1 * t23 * t58, -t46, -t47, 0, t37 * t42 - t39 * t41, t37 * t41 + t39 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, qJ(4) * t73, -0.2e1 * t58, 0.2e1 * t55, 0, 0, 0, 0.2e1 * qJD(4) * t37 + 0.2e1 * t39 * t60, 0.2e1 * qJD(4) * t39 - 0.2e1 * t37 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, t9, 0, 0, 0, 0, 0, -t46, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t44, -t18, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t62, 0, t37 * t61, t39 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
