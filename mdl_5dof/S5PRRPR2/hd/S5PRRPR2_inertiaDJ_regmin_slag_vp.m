% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPR2
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
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:35
% EndTime: 2019-12-05 16:17:36
% DurationCPUTime: 0.24s
% Computational Cost: add. (136->48), mult. (431->94), div. (0->0), fcn. (267->6), ass. (0->53)
t34 = cos(qJ(3));
t62 = t34 * pkin(2);
t52 = pkin(2) * qJD(3);
t45 = t34 * t52;
t19 = qJD(4) + t45;
t29 = sin(pkin(9));
t27 = t29 ^ 2;
t14 = t27 * t19;
t32 = sin(qJ(3));
t24 = t32 * pkin(2) + qJ(4);
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t46 = t32 * t52;
t30 = cos(pkin(9));
t57 = t30 * t33;
t58 = t30 * t31;
t17 = -t30 * pkin(4) - t29 * pkin(7) - pkin(3);
t8 = t17 - t62;
t2 = -t31 * t46 - t19 * t57 + (t24 * t58 - t33 * t8) * qJD(5);
t61 = t33 * t14 - t2 * t30;
t47 = qJD(5) * t33;
t43 = t27 * t47;
t60 = t31 * t14 + t24 * t43;
t25 = t27 * qJD(4);
t50 = qJD(4) * t30;
t51 = qJ(4) * t30;
t5 = -t33 * t50 + (-t17 * t33 + t31 * t51) * qJD(5);
t59 = t33 * t25 - t5 * t30;
t28 = t30 ^ 2;
t56 = t28 * t19 + t14;
t49 = qJD(5) * t27;
t40 = qJ(4) * t49;
t55 = t31 * t25 + t33 * t40;
t54 = t28 * qJD(4) + t25;
t53 = t27 + t28;
t48 = qJD(5) * t31;
t44 = t27 * t48;
t42 = t29 * t48;
t41 = t29 * t47;
t39 = t53 * t19;
t38 = t53 * qJD(4);
t37 = 0.2e1 * qJD(5) * t29 * t30;
t36 = t29 * t46;
t35 = t30 * t46;
t23 = t30 * t47;
t22 = t30 * t48;
t16 = -0.2e1 * t31 * t43;
t13 = t33 * t37;
t12 = t31 * t37;
t7 = 0.2e1 * (t31 ^ 2 - t33 ^ 2) * t49;
t6 = -t31 * t50 + (-t17 * t31 - t33 * t51) * qJD(5);
t3 = t33 * t46 - t19 * t58 + (-t24 * t57 - t31 * t8) * qJD(5);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -0.2e1 * t46, -0.2e1 * t45, -0.2e1 * t35, 0.2e1 * t36, 0.2e1 * t56, 0.2e1 * (-pkin(3) - t62) * t46 + 0.2e1 * t24 * t39, t16, t7, t12, t13, 0, -0.2e1 * t3 * t30 + 0.2e1 * t60, -0.2e1 * t24 * t44 + 0.2e1 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t46, -t45, -t35, t36, t54 + t56, -pkin(3) * t46 + qJ(4) * t39 + t24 * t38, t16, t7, t12, t13, 0, (-t3 - t6) * t30 + t55 + t60, (-qJ(4) - t24) * t44 + t59 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t54, 0.2e1 * qJ(4) * t38, t16, t7, t12, t13, 0, -0.2e1 * t6 * t30 + 0.2e1 * t55, -0.2e1 * t31 * t40 + 0.2e1 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, t22, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t41, 0, t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t41, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
