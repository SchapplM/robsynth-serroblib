% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:46:06
% EndTime: 2021-01-15 12:46:08
% DurationCPUTime: 0.33s
% Computational Cost: add. (530->69), mult. (1164->121), div. (0->0), fcn. (1005->6), ass. (0->43)
t54 = qJD(3) + qJD(4);
t33 = sin(qJ(4));
t34 = sin(qJ(3));
t35 = cos(qJ(3));
t51 = cos(qJ(4));
t22 = t33 * t35 + t51 * t34;
t12 = t54 * t22;
t53 = t12 * pkin(4);
t28 = sin(pkin(8)) * pkin(1) + pkin(6);
t52 = pkin(7) + t28;
t40 = t51 * qJD(4);
t41 = t51 * qJD(3);
t49 = t33 * t34;
t11 = t54 * t49 + (-t41 - t40) * t35;
t50 = t22 * t11;
t48 = qJD(4) * t33;
t47 = t34 * qJD(3);
t46 = t35 * qJD(3);
t45 = 0.2e1 * t46;
t44 = pkin(3) * t47;
t43 = pkin(3) * t48;
t42 = t51 * pkin(3);
t29 = -cos(pkin(8)) * pkin(1) - pkin(2);
t39 = pkin(3) * t40;
t38 = qJD(3) * t33 * t52;
t23 = -t35 * pkin(3) + t29;
t37 = t52 * t41;
t19 = t52 * t34;
t20 = t52 * t35;
t36 = t33 * t19 - t51 * t20;
t3 = t19 * t40 + t20 * t48 + t34 * t37 + t35 * t38;
t4 = t36 * qJD(4) + t34 * t38 - t35 * t37;
t31 = t42 + pkin(4);
t27 = -0.2e1 * t39;
t26 = -0.2e1 * t43;
t21 = -t51 * t35 + t49;
t13 = t21 * pkin(4) + t23;
t9 = t44 + t53;
t6 = -t21 * qJ(5) - t36;
t5 = -t22 * qJ(5) - t51 * t19 - t33 * t20;
t2 = t11 * qJ(5) - t22 * qJD(5) + t4;
t1 = t12 * qJ(5) + t21 * qJD(5) + t3;
t7 = [0, 0, 0, 0, t34 * t45, 0.2e1 * (-t34 ^ 2 + t35 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t29 * t47, t29 * t45, -0.2e1 * t50, 0.2e1 * t21 * t11 - 0.2e1 * t22 * t12, 0, 0, 0, 0.2e1 * t23 * t12 + 0.2e1 * t21 * t44, -0.2e1 * t23 * t11 + 0.2e1 * t22 * t44, 0.2e1 * t13 * t12 + 0.2e1 * t9 * t21, -0.2e1 * t13 * t11 + 0.2e1 * t9 * t22, 0.2e1 * t1 * t21 + 0.2e1 * t5 * t11 - 0.2e1 * t6 * t12 - 0.2e1 * t2 * t22, -0.2e1 * t6 * t1 + 0.2e1 * t13 * t9 + 0.2e1 * t5 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t22 - t6 * t11 - t5 * t12 - t2 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t21 * t12 - 0.2e1 * t50; 0, 0, 0, 0, 0, 0, t46, -t47, 0, -t28 * t46, t28 * t47, 0, 0, -t11, -t12, 0, t4, t3, t2, t1, t31 * t11 + (-t12 * t33 + (-t51 * t21 + t22 * t33) * qJD(4)) * pkin(3), t2 * t31 + (-t1 * t33 + (-t33 * t5 + t51 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t46, 0, 0, 0, 0, 0, -t12, t11, -t12, t11, 0, -t12 * t31 + (-t11 * t33 + (t21 * t33 + t51 * t22) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t27, t26, t27, 0, 0.2e1 * (t42 - t31) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, t4, t3, t2, t1, t11 * pkin(4), t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, -t12, t11, 0, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t39, -t43, -t39, 0, -pkin(4) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
