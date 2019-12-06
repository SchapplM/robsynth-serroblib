% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPR4
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
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:22
% EndTime: 2019-12-05 16:23:24
% DurationCPUTime: 0.36s
% Computational Cost: add. (414->74), mult. (1106->156), div. (0->0), fcn. (1038->8), ass. (0->53)
t48 = sin(qJ(2));
t44 = sin(pkin(9));
t45 = cos(pkin(9));
t47 = sin(qJ(3));
t50 = cos(qJ(3));
t52 = t44 * t47 - t45 * t50;
t28 = t52 * t48;
t62 = pkin(3) * t44;
t61 = -qJ(4) - pkin(6);
t54 = qJD(3) * t61;
t29 = t50 * qJD(4) + t47 * t54;
t30 = -t47 * qJD(4) + t50 * t54;
t12 = t45 * t29 + t44 * t30;
t38 = t61 * t47;
t39 = t61 * t50;
t18 = t44 * t38 - t45 * t39;
t60 = t47 * qJD(3);
t59 = t48 * qJD(2);
t58 = t50 * qJD(3);
t51 = cos(qJ(2));
t57 = t51 * qJD(2);
t56 = -0.2e1 * pkin(2) * qJD(3);
t43 = pkin(3) * t60;
t55 = t51 * t60;
t42 = -t50 * pkin(3) - pkin(2);
t11 = -t44 * t29 + t45 * t30;
t17 = t45 * t38 + t44 * t39;
t34 = t44 * t50 + t45 * t47;
t46 = sin(qJ(5));
t49 = cos(qJ(5));
t53 = -t46 * t34 - t49 * t52;
t16 = t49 * t34 - t46 * t52;
t31 = t34 * qJD(3);
t41 = t45 * pkin(3) + pkin(4);
t32 = -t44 * t60 + t45 * t58;
t27 = t34 * t48;
t25 = (-t41 * t46 - t49 * t62) * qJD(5);
t24 = (-t41 * t49 + t46 * t62) * qJD(5);
t20 = pkin(4) * t52 + t42;
t19 = t31 * pkin(4) + t43;
t14 = -pkin(7) * t52 + t18;
t13 = -t34 * pkin(7) + t17;
t10 = -t48 * t31 - t52 * t57;
t9 = qJD(3) * t28 - t34 * t57;
t8 = -t31 * pkin(7) + t12;
t7 = -t32 * pkin(7) + t11;
t6 = t16 * qJD(5) + t49 * t31 + t46 * t32;
t5 = t53 * qJD(5) - t46 * t31 + t49 * t32;
t4 = -t46 * t10 + t49 * t9 + (t27 * t46 + t28 * t49) * qJD(5);
t3 = -t49 * t10 - t46 * t9 + (t27 * t49 - t28 * t46) * qJD(5);
t2 = -t46 * t8 + t49 * t7 + (-t13 * t46 - t14 * t49) * qJD(5);
t1 = -t46 * t7 - t49 * t8 + (-t13 * t49 + t14 * t46) * qJD(5);
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t28 * t10 - 0.2e1 * t27 * t9 - 0.2e1 * t48 * t57, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t59, -t57, 0, 0, 0, 0, 0, -t50 * t59 - t55, t47 * t59 - t51 * t58, -t10 * t52 + t27 * t32 + t28 * t31 - t9 * t34, -pkin(3) * t55 + t10 * t18 - t27 * t11 - t28 * t12 + t9 * t17 + t42 * t59, 0, 0, 0, 0, 0, -t51 * t6 - t53 * t59, t16 * t59 - t51 * t5; 0, 0, 0, 0, 0.2e1 * t47 * t58, 0.2e1 * (-t47 ^ 2 + t50 ^ 2) * qJD(3), 0, 0, 0, t47 * t56, t50 * t56, -0.2e1 * t11 * t34 - 0.2e1 * t12 * t52 - 0.2e1 * t17 * t32 - 0.2e1 * t18 * t31, 0.2e1 * t17 * t11 + 0.2e1 * t18 * t12 + 0.2e1 * t42 * t43, 0.2e1 * t16 * t5, -0.2e1 * t16 * t6 + 0.2e1 * t5 * t53, 0, 0, 0, -0.2e1 * t19 * t53 + 0.2e1 * t20 * t6, 0.2e1 * t19 * t16 + 0.2e1 * t20 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47 * t57 - t48 * t58, t48 * t60 - t50 * t57, 0, (t10 * t44 + t45 * t9) * pkin(3), 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, t58, -t60, 0, -pkin(6) * t58, pkin(6) * t60, (-t31 * t44 - t32 * t45) * pkin(3), (t11 * t45 + t12 * t44) * pkin(3), 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t25, 0.2e1 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
