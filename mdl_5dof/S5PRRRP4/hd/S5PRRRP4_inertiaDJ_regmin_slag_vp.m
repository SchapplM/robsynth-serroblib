% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:20
% EndTime: 2019-12-05 16:46:22
% DurationCPUTime: 0.36s
% Computational Cost: add. (284->66), mult. (761->105), div. (0->0), fcn. (601->6), ass. (0->51)
t34 = sin(qJ(3));
t35 = sin(qJ(2));
t52 = cos(qJ(3));
t53 = cos(qJ(2));
t18 = t34 * t35 - t52 * t53;
t49 = pkin(2) * qJD(3);
t60 = t49 * t52;
t33 = sin(qJ(4));
t31 = t33 ^ 2;
t36 = cos(qJ(4));
t32 = t36 ^ 2;
t59 = t31 + t32;
t58 = -pkin(4) * t33 + qJ(5) * t36;
t57 = qJD(2) + qJD(3);
t56 = 2 * qJD(5);
t30 = t36 * qJD(4);
t47 = t33 * qJD(4);
t14 = pkin(4) * t47 - qJ(5) * t30 - t33 * qJD(5);
t44 = t34 * t49;
t8 = t14 + t44;
t54 = -t14 - t8;
t41 = t52 * pkin(2);
t28 = -t41 - pkin(3);
t50 = t28 * t30 + t33 * t44;
t46 = pkin(3) * t47;
t45 = pkin(3) * t30;
t43 = pkin(7) * t47;
t42 = pkin(7) * t30;
t6 = t18 * t57;
t1 = t59 * t6;
t40 = qJD(3) * t41;
t38 = -t36 * pkin(4) - t33 * qJ(5);
t37 = t28 * t47 - t36 * t44;
t22 = -pkin(3) + t38;
t19 = t34 * t53 + t52 * t35;
t13 = t38 * qJD(4) + t36 * qJD(5);
t12 = t59 * t60;
t27 = pkin(2) * t34 + pkin(7);
t24 = 0.2e1 * t33 * t30;
t17 = 0.2e1 * (-t31 + t32) * qJD(4);
t16 = -t41 + t22;
t15 = t22 * t47;
t11 = t16 * t47;
t10 = t27 * t30 + t33 * t40;
t9 = t27 * t47 - t36 * t40;
t7 = t57 * t19;
t5 = t18 * t47 - t36 * t7;
t4 = t18 * t30 + t7 * t33;
t3 = t19 * t30 - t33 * t6;
t2 = t19 * t47 + t36 * t6;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t19 * t1 + 0.2e1 * t18 * t7; 0, 0, -t35 * qJD(2), -t53 * qJD(2), 0, -t7, t6, 0, 0, 0, 0, 0, t5, t4, t5, -t1, -t4, -t27 * t1 + t19 * t12 + t7 * t16 + t18 * t8; 0, 0, 0, 0, 0, -0.2e1 * t44, -0.2e1 * t40, t24, t17, 0, 0, 0, 0.2e1 * t37, 0.2e1 * t50, -0.2e1 * t36 * t8 + 0.2e1 * t11, 0.2e1 * t12, -0.2e1 * t16 * t30 - 0.2e1 * t8 * t33, 0.2e1 * t27 * t12 + 0.2e1 * t16 * t8; 0, 0, 0, 0, 0, -t7, t6, 0, 0, 0, 0, 0, t5, t4, t5, -t1, -t4, -pkin(7) * t1 + t14 * t18 + t22 * t7; 0, 0, 0, 0, 0, -t44, -t40, t24, t17, 0, 0, 0, t37 - t46, -t45 + t50, t54 * t36 + t11 + t15, t12, t54 * t33 + (-t16 - t22) * t30, pkin(7) * t12 + t16 * t14 + t8 * t22; 0, 0, 0, 0, 0, 0, 0, t24, t17, 0, 0, 0, -0.2e1 * t46, -0.2e1 * t45, -0.2e1 * t14 * t36 + 0.2e1 * t15, 0, -0.2e1 * t14 * t33 - 0.2e1 * t22 * t30, 0.2e1 * t22 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2, -t3, 0, -t2, t13 * t19 - t58 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t47, 0, -t10, t9, -t10, t13, -t9, t13 * t27 + t58 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t47, 0, -t42, t43, -t42, t13, -t43, t13 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, qJ(5) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t20;
