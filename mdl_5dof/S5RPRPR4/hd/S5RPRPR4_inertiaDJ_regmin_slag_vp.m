% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:57
% EndTime: 2019-12-05 17:53:59
% DurationCPUTime: 0.32s
% Computational Cost: add. (419->61), mult. (937->120), div. (0->0), fcn. (862->8), ass. (0->44)
t42 = sin(pkin(9));
t43 = cos(pkin(9));
t46 = sin(qJ(3));
t47 = cos(qJ(3));
t35 = -t42 * t46 + t43 * t47;
t36 = t42 * t47 + t43 * t46;
t45 = sin(qJ(5));
t56 = cos(qJ(5));
t13 = -t56 * t35 + t45 * t36;
t57 = pkin(3) * t42;
t38 = sin(pkin(8)) * pkin(1) + pkin(6);
t54 = qJ(4) + t38;
t49 = qJD(3) * t54;
t20 = t47 * qJD(4) - t46 * t49;
t21 = -t46 * qJD(4) - t47 * t49;
t8 = t43 * t20 + t42 * t21;
t33 = t54 * t46;
t34 = t54 * t47;
t12 = -t42 * t33 + t43 * t34;
t53 = t46 * qJD(3);
t52 = t47 * qJD(3);
t51 = 0.2e1 * t52;
t41 = pkin(3) * t53;
t40 = -cos(pkin(8)) * pkin(1) - pkin(2);
t7 = -t42 * t20 + t43 * t21;
t11 = -t43 * t33 - t42 * t34;
t48 = -t47 * pkin(3) + t40;
t14 = t45 * t35 + t56 * t36;
t39 = t43 * pkin(3) + pkin(4);
t32 = -t42 * t53 + t43 * t52;
t31 = t36 * qJD(3);
t23 = (-t39 * t45 - t56 * t57) * qJD(5);
t22 = (-t56 * t39 + t45 * t57) * qJD(5);
t19 = t31 * pkin(4) + t41;
t18 = -t35 * pkin(4) + t48;
t10 = t35 * pkin(7) + t12;
t9 = -t36 * pkin(7) + t11;
t6 = -t31 * pkin(7) + t8;
t5 = -t32 * pkin(7) + t7;
t4 = t14 * qJD(5) + t56 * t31 + t45 * t32;
t3 = t13 * qJD(5) + t45 * t31 - t56 * t32;
t2 = t56 * t5 - t45 * t6 + (-t56 * t10 - t45 * t9) * qJD(5);
t1 = -t56 * t6 - t45 * t5 + (t10 * t45 - t56 * t9) * qJD(5);
t15 = [0, 0, 0, 0, t46 * t51, 0.2e1 * (-t46 ^ 2 + t47 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t40 * t53, t40 * t51, -0.2e1 * t11 * t32 - 0.2e1 * t12 * t31 + 0.2e1 * t8 * t35 - 0.2e1 * t7 * t36, 0.2e1 * t11 * t7 + 0.2e1 * t12 * t8 + 0.2e1 * t48 * t41, -0.2e1 * t14 * t3, 0.2e1 * t3 * t13 - 0.2e1 * t14 * t4, 0, 0, 0, 0.2e1 * t19 * t13 + 0.2e1 * t18 * t4, 0.2e1 * t19 * t14 - 0.2e1 * t18 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t31 + t12 * t32 + t7 * t35 + t8 * t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t35 * t31 + 0.2e1 * t36 * t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t52, -t53, 0, -t38 * t52, t38 * t53, (-t31 * t42 - t32 * t43) * pkin(3), (t42 * t8 + t43 * t7) * pkin(3), 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t52, 0, (-t31 * t43 + t32 * t42) * pkin(3), 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t23, 0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
