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
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 15:52:41
% EndTime: 2021-01-15 15:52:44
% DurationCPUTime: 0.38s
% Computational Cost: add. (446->80), mult. (1206->169), div. (0->0), fcn. (1126->8), ass. (0->55)
t44 = sin(pkin(9));
t45 = cos(pkin(9));
t47 = sin(qJ(3));
t50 = cos(qJ(3));
t31 = t44 * t47 - t45 * t50;
t48 = sin(qJ(2));
t26 = t31 * t48;
t62 = pkin(3) * t44;
t61 = qJ(4) + pkin(6);
t37 = t61 * t47;
t38 = t61 * t50;
t18 = -t44 * t37 + t45 * t38;
t60 = t47 * qJD(3);
t59 = t48 * qJD(2);
t58 = t50 * qJD(3);
t51 = cos(qJ(2));
t57 = t51 * qJD(2);
t56 = -0.2e1 * pkin(2) * qJD(3);
t43 = pkin(3) * t60;
t55 = t51 * t60;
t54 = t47 * t57;
t53 = t50 * t57;
t42 = -pkin(3) * t50 - pkin(2);
t52 = qJD(3) * t61;
t27 = t50 * qJD(4) - t47 * t52;
t28 = -t47 * qJD(4) - t50 * t52;
t11 = -t27 * t44 + t45 * t28;
t17 = -t45 * t37 - t38 * t44;
t12 = t27 * t45 + t28 * t44;
t32 = t44 * t50 + t45 * t47;
t46 = sin(qJ(5));
t49 = cos(qJ(5));
t15 = t31 * t49 + t32 * t46;
t16 = -t31 * t46 + t32 * t49;
t29 = t32 * qJD(3);
t41 = pkin(3) * t45 + pkin(4);
t30 = -t44 * t60 + t45 * t58;
t25 = t32 * t48;
t23 = (-t41 * t46 - t49 * t62) * qJD(5);
t22 = (-t41 * t49 + t46 * t62) * qJD(5);
t20 = pkin(4) * t31 + t42;
t19 = pkin(4) * t29 + t43;
t14 = -pkin(7) * t31 + t18;
t13 = -pkin(7) * t32 + t17;
t10 = t29 * t48 + t44 * t54 - t45 * t53;
t9 = qJD(3) * t26 - t32 * t57;
t8 = -pkin(7) * t29 + t12;
t7 = -pkin(7) * t30 + t11;
t6 = qJD(5) * t16 + t49 * t29 + t46 * t30;
t5 = -qJD(5) * t15 - t46 * t29 + t49 * t30;
t4 = t46 * t10 + t49 * t9 + (t25 * t46 + t26 * t49) * qJD(5);
t3 = t49 * t10 - t46 * t9 + (t25 * t49 - t26 * t46) * qJD(5);
t2 = -t46 * t8 + t49 * t7 + (-t13 * t46 - t14 * t49) * qJD(5);
t1 = -t46 * t7 - t49 * t8 + (-t13 * t49 + t14 * t46) * qJD(5);
t21 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t26 - 0.2e1 * t25 * t9 - 0.2e1 * t48 * t57, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t59, -t57, 0, 0, 0, 0, 0, -t50 * t59 - t55, t47 * t59 - t51 * t58, -t29 * t51 + t31 * t59, -t30 * t51 + t32 * t59, t10 * t31 + t25 * t30 + t26 * t29 - t32 * t9, -pkin(3) * t55 - t10 * t18 - t11 * t25 - t12 * t26 + t17 * t9 + t42 * t59, 0, 0, 0, 0, 0, t15 * t59 - t51 * t6, t16 * t59 - t5 * t51; 0, 0, 0, 0, 0.2e1 * t47 * t58, 0.2e1 * (-t47 ^ 2 + t50 ^ 2) * qJD(3), 0, 0, 0, t47 * t56, t50 * t56, 0.2e1 * t29 * t42 + 0.2e1 * t31 * t43, 0.2e1 * t30 * t42 + 0.2e1 * t32 * t43, -0.2e1 * t11 * t32 - 0.2e1 * t12 * t31 - 0.2e1 * t17 * t30 - 0.2e1 * t18 * t29, 0.2e1 * t11 * t17 + 0.2e1 * t12 * t18 + 0.2e1 * t42 * t43, 0.2e1 * t16 * t5, -0.2e1 * t15 * t5 - 0.2e1 * t16 * t6, 0, 0, 0, 0.2e1 * t15 * t19 + 0.2e1 * t20 * t6, 0.2e1 * t16 * t19 + 0.2e1 * t20 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t58 - t54, t48 * t60 - t53, t9, t10, 0, (-t10 * t44 + t45 * t9) * pkin(3), 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, t58, -t60, 0, -pkin(6) * t58, pkin(6) * t60, t11, -t12, (-t29 * t44 - t30 * t45) * pkin(3), (t11 * t45 + t12 * t44) * pkin(3), 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t23, 0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t30, 0, t43, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t21;
