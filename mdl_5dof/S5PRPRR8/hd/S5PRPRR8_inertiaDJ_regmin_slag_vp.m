% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:26
% EndTime: 2019-12-05 16:04:28
% DurationCPUTime: 0.37s
% Computational Cost: add. (165->75), mult. (519->155), div. (0->0), fcn. (425->8), ass. (0->68)
t69 = 2 * qJD(3);
t27 = cos(qJ(4));
t68 = t27 * pkin(8);
t21 = sin(pkin(5));
t25 = sin(qJ(2));
t67 = t21 * t25;
t28 = cos(qJ(2));
t66 = t21 * t28;
t24 = sin(qJ(4));
t29 = -pkin(2) - pkin(7);
t65 = t24 * t29;
t64 = t27 * t29;
t26 = cos(qJ(5));
t19 = t26 ^ 2;
t23 = sin(qJ(5));
t63 = t23 ^ 2 - t19;
t18 = t24 ^ 2;
t20 = t27 ^ 2;
t62 = t18 - t20;
t61 = t18 + t20;
t22 = cos(pkin(5));
t7 = t22 * t24 + t27 * t66;
t60 = qJD(4) * t7;
t59 = qJD(2) * t28;
t58 = qJD(4) * t26;
t57 = qJD(5) * t20;
t56 = qJD(5) * t23;
t55 = qJD(5) * t26;
t54 = qJD(5) * t27;
t53 = t24 * qJD(4);
t52 = t27 * qJD(4);
t51 = qJ(3) * qJD(4);
t50 = -0.2e1 * pkin(4) * qJD(5);
t49 = t26 * t65;
t48 = t23 * t54;
t47 = t26 * t54;
t16 = qJD(2) * t67;
t46 = t21 * t59;
t45 = t23 * t55;
t44 = t26 * t53;
t43 = t24 * t52;
t42 = t29 * t52;
t41 = t63 * qJD(5);
t40 = t62 * qJD(4);
t39 = t23 * t42;
t38 = t23 * t44;
t37 = t26 * t42;
t36 = pkin(4) * t27 + pkin(8) * t24;
t35 = t24 * pkin(4) - t68;
t8 = t22 * t27 - t24 * t66;
t34 = -t8 * t23 + t26 * t67;
t33 = t23 * t67 + t8 * t26;
t6 = t8 * qJD(4) - t27 * t16;
t32 = t6 * t23 + t7 * t55;
t31 = -t6 * t26 + t7 * t56;
t14 = qJ(3) + t35;
t30 = t14 * t52 - t29 * t57;
t13 = t36 * qJD(4) + qJD(3);
t12 = t23 * t53 - t47;
t11 = t23 * t52 + t24 * t55;
t10 = -t44 - t48;
t9 = t24 * t56 - t26 * t52;
t5 = -t24 * t16 + t60;
t4 = -t39 + t26 * t13 + (-t23 * t14 - t49) * qJD(5);
t3 = -t37 - t23 * t13 + (-t26 * t14 + t23 * t65) * qJD(5);
t2 = t34 * qJD(5) + t23 * t46 - t5 * t26;
t1 = -t33 * qJD(5) + t5 * t23 + t26 * t46;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t16, -t46, t16, t46, (qJD(3) * t25 + (-pkin(2) * t25 + qJ(3) * t28) * qJD(2)) * t21, 0, 0, 0, 0, 0, (t24 * t59 + t25 * t52) * t21, (-t25 * t53 + t27 * t59) * t21, 0, 0, 0, 0, 0, (-t23 * t60 + t1) * t24 + (t34 * qJD(4) + t32) * t27, (-t7 * t58 - t2) * t24 + (-t33 * qJD(4) - t31) * t27; 0, 0, 0, 0, 0, t69, qJ(3) * t69, -0.2e1 * t43, 0.2e1 * t40, 0, 0, 0, 0.2e1 * qJD(3) * t24 + 0.2e1 * t27 * t51, 0.2e1 * qJD(3) * t27 - 0.2e1 * t24 * t51, -0.2e1 * t19 * t43 - 0.2e1 * t20 * t45, 0.4e1 * t27 * t38 + 0.2e1 * t63 * t57, -0.2e1 * t24 * t48 - 0.2e1 * t62 * t58, 0.2e1 * t23 * t40 - 0.2e1 * t24 * t47, 0.2e1 * t43, 0.2e1 * t30 * t26 + 0.2e1 * (t4 + t39) * t24, 0.2e1 * (t3 + t37) * t24 - 0.2e1 * t30 * t23; 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 * t55, t61 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, t31, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t52, 0, -t29 * t53, -t42, -t27 * t41 - t38, -0.4e1 * t27 * t45 + t63 * t53, t11, -t9, 0, (-t23 * t64 - t36 * t26) * qJD(5) + (t35 * t23 - t49) * qJD(4), (t36 * t23 - t26 * t64) * qJD(5) + (-t26 * t68 + (pkin(4) * t26 + t23 * t29) * t24) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t52, 0, 0, 0, 0, 0, t10, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t45, -0.2e1 * t41, 0, 0, 0, t23 * t50, t26 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12, t52, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t56, 0, -pkin(8) * t55, pkin(8) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
