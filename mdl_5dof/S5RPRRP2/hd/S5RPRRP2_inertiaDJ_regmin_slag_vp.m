% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP2
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
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:36:27
% EndTime: 2021-01-15 12:36:29
% DurationCPUTime: 0.29s
% Computational Cost: add. (328->72), mult. (724->115), div. (0->0), fcn. (500->6), ass. (0->48)
t29 = cos(pkin(8)) * pkin(1) + pkin(2);
t40 = sin(qJ(3));
t42 = cos(qJ(3));
t55 = pkin(1) * sin(pkin(8));
t57 = -t42 * t29 + t40 * t55;
t18 = -pkin(3) + t57;
t41 = cos(qJ(4));
t54 = t41 * pkin(4);
t14 = t18 - t54;
t35 = t41 * qJD(4);
t39 = sin(qJ(4));
t43 = t40 * t29 + t42 * t55;
t17 = t43 * qJD(3);
t33 = t39 * qJD(4);
t32 = pkin(4) * t33;
t9 = t32 + t17;
t56 = t14 * t35 + t9 * t39;
t52 = -qJ(5) - pkin(7);
t51 = t17 * t39 + t18 * t35;
t30 = -pkin(3) - t54;
t37 = t39 ^ 2;
t50 = qJD(4) * t37 * pkin(4) + t30 * t35;
t19 = pkin(7) + t43;
t49 = qJ(5) + t19;
t48 = pkin(3) * t33;
t47 = pkin(3) * t35;
t46 = pkin(4) * t35;
t5 = t14 * t33;
t23 = t30 * t33;
t45 = t39 * t35;
t44 = -t17 * t41 + t18 * t33;
t36 = t41 * qJ(5);
t34 = t41 * qJD(5);
t27 = 0.2e1 * t45;
t26 = t41 * pkin(7) + t36;
t25 = t52 * t39;
t22 = 0.2e1 * (t41 ^ 2 - t37) * qJD(4);
t21 = -t39 * qJD(5) + t52 * t35;
t20 = -t52 * t33 - t34;
t16 = t57 * qJD(3);
t15 = t20 * t41;
t10 = t41 * t16;
t8 = t41 * t19 + t36;
t7 = t49 * t39;
t3 = (-qJD(5) + t16) * t39 - t49 * t35;
t2 = t49 * t33 + t10 - t34;
t1 = t2 * t41;
t4 = [0, 0, 0, 0, 0, -0.2e1 * t17, 0.2e1 * t16, t27, t22, 0, 0, 0, 0.2e1 * t44, 0.2e1 * t51, -0.2e1 * t9 * t41 + 0.2e1 * t5, 0.2e1 * t56, -0.2e1 * t3 * t39 - 0.2e1 * t1 + 0.2e1 * (-t39 * t8 + t41 * t7) * qJD(4), 0.2e1 * t14 * t9 - 0.2e1 * t8 * t2 - 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t39 + t3 * t41 + (t39 * t7 + t41 * t8) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t17, t16, t27, t22, 0, 0, 0, t44 - t48, -t47 + t51, t23 + t5 + (-t9 - t32) * t41, t50 + t56, -t1 - t15 + (-t21 - t3) * t39 + ((-t25 + t7) * t41 + (-t26 - t8) * t39) * qJD(4), pkin(4) * t5 - t2 * t26 - t8 * t20 - t7 * t21 + t3 * t25 + t9 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39 * t20 + t41 * t21 + (-t25 * t39 + t26 * t41) * qJD(4); 0, 0, 0, 0, 0, 0, 0, t27, t22, 0, 0, 0, -0.2e1 * t48, -0.2e1 * t47, -0.2e1 * pkin(4) * t45 + 0.2e1 * t23, 0.2e1 * t50, -0.2e1 * t21 * t39 - 0.2e1 * t15 + 0.2e1 * (-t25 * t41 - t26 * t39) * qJD(4), 0.2e1 * pkin(4) * t23 - 0.2e1 * t26 * t20 + 0.2e1 * t25 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t33, 0, t39 * t16 - t19 * t35, t19 * t33 + t10, t3, t2, -t46, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t35, -t33, -t35, 0, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t33, 0, -pkin(7) * t35, pkin(7) * t33, t21, t20, -t46, t21 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t35, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t35, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
