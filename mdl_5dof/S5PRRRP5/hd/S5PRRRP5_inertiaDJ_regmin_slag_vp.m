% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRP5
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
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:33:42
% EndTime: 2021-01-15 16:33:45
% DurationCPUTime: 0.38s
% Computational Cost: add. (505->77), mult. (1347->151), div. (0->0), fcn. (1189->6), ass. (0->51)
t38 = cos(qJ(4));
t39 = cos(qJ(3));
t35 = sin(qJ(4));
t36 = sin(qJ(3));
t57 = t35 * t36;
t23 = -t38 * t39 + t57;
t37 = sin(qJ(2));
t17 = t23 * t37;
t60 = qJD(3) + qJD(4);
t59 = pkin(6) + pkin(7);
t58 = t38 * pkin(3);
t56 = qJD(4) * t35;
t55 = qJD(4) * t38;
t54 = t36 * qJD(3);
t53 = t37 * qJD(2);
t52 = t39 * qJD(3);
t40 = cos(qJ(2));
t51 = t40 * qJD(2);
t50 = -0.2e1 * pkin(2) * qJD(3);
t49 = pkin(3) * t54;
t48 = pkin(3) * t56;
t47 = pkin(3) * t55;
t46 = t36 * t51;
t45 = t39 * t51;
t34 = -t39 * pkin(3) - pkin(2);
t44 = qJD(3) * t59;
t43 = t36 * t44;
t42 = t39 * t44;
t26 = t59 * t36;
t27 = t59 * t39;
t41 = t35 * t26 - t38 * t27;
t24 = t35 * t39 + t38 * t36;
t5 = t26 * t55 + t27 * t56 + t35 * t42 + t38 * t43;
t6 = t41 * qJD(4) + t35 * t43 - t38 * t42;
t14 = t60 * t24;
t33 = pkin(4) + t58;
t31 = -0.2e1 * t47;
t30 = -0.2e1 * t48;
t16 = t24 * t37;
t15 = t23 * pkin(4) + t34;
t13 = -t38 * t52 - t39 * t55 + t60 * t57;
t11 = t14 * pkin(4) + t49;
t10 = -t23 * qJ(5) - t41;
t9 = -t24 * qJ(5) - t38 * t26 - t35 * t27;
t8 = t40 * t13 + t24 * t53;
t7 = -t40 * t14 + t23 * t53;
t4 = t60 * t17 - t24 * t51;
t3 = t14 * t37 + t35 * t46 - t38 * t45;
t2 = t13 * qJ(5) - t24 * qJD(5) + t6;
t1 = t14 * qJ(5) + t23 * qJD(5) + t5;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t16 * t4 + 0.2e1 * t17 * t3 - 0.2e1 * t37 * t51; 0, 0, -t53, -t51, 0, 0, 0, 0, 0, -t39 * t53 - t40 * t54, t36 * t53 - t40 * t52, 0, 0, 0, 0, 0, t7, t8, t7, t8, -t16 * t13 + t17 * t14 + t3 * t23 - t4 * t24, t17 * t1 - t3 * t10 - t40 * t11 + t15 * t53 - t16 * t2 + t4 * t9; 0, 0, 0, 0, 0.2e1 * t36 * t52, 0.2e1 * (-t36 ^ 2 + t39 ^ 2) * qJD(3), 0, 0, 0, t36 * t50, t39 * t50, -0.2e1 * t24 * t13, 0.2e1 * t13 * t23 - 0.2e1 * t24 * t14, 0, 0, 0, 0.2e1 * t34 * t14 + 0.2e1 * t23 * t49, -0.2e1 * t34 * t13 + 0.2e1 * t24 * t49, 0.2e1 * t11 * t23 + 0.2e1 * t15 * t14, 0.2e1 * t11 * t24 - 0.2e1 * t15 * t13, 0.2e1 * t1 * t23 - 0.2e1 * t10 * t14 + 0.2e1 * t9 * t13 - 0.2e1 * t2 * t24, -0.2e1 * t10 * t1 + 0.2e1 * t15 * t11 + 0.2e1 * t9 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t52 - t46, t37 * t54 - t45, 0, 0, 0, 0, 0, t4, t3, t4, t3, 0, t4 * t33 + (-t3 * t35 + (t16 * t35 - t17 * t38) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, t52, -t54, 0, -pkin(6) * t52, pkin(6) * t54, 0, 0, -t13, -t14, 0, t6, t5, t2, t1, t33 * t13 + (-t14 * t35 + (-t23 * t38 + t24 * t35) * qJD(4)) * pkin(3), t2 * t33 + (-t1 * t35 + (t10 * t38 - t35 * t9) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t31, t30, t31, 0, 0.2e1 * (-t33 + t58) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, t4, t3, 0, t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, t6, t5, t2, t1, t13 * pkin(4), t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, -t48, -t47, 0, -pkin(4) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
