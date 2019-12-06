% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:51
% EndTime: 2019-12-05 15:12:53
% DurationCPUTime: 0.43s
% Computational Cost: add. (301->46), mult. (842->88), div. (0->0), fcn. (837->8), ass. (0->46)
t26 = sin(qJ(5));
t24 = t26 ^ 2;
t28 = cos(qJ(5));
t25 = t28 ^ 2;
t53 = t24 + t25;
t52 = qJD(3) + qJD(4);
t43 = sin(pkin(9));
t44 = cos(pkin(9));
t47 = sin(qJ(3));
t49 = cos(qJ(3));
t14 = t49 * t43 + t47 * t44;
t30 = t14 * qJD(3);
t27 = sin(qJ(4));
t32 = t47 * t43 - t49 * t44;
t31 = t27 * t32;
t48 = cos(qJ(4));
t36 = qJD(4) * t48;
t5 = t14 * t36 + t48 * t30 - t52 * t31;
t29 = t48 * t32;
t6 = t27 * t14 + t29;
t51 = t6 * t5;
t50 = t27 * t6;
t22 = -t48 * pkin(3) - pkin(4);
t23 = t28 * qJD(5);
t42 = qJD(4) * t27;
t38 = pkin(3) * t42;
t46 = t22 * t23 + t26 * t38;
t45 = pkin(3) * qJD(4);
t41 = t26 * qJD(5);
t40 = pkin(4) * t41;
t39 = pkin(4) * t23;
t37 = t26 * t23;
t4 = t14 * t42 + t27 * t30 + t52 * t29;
t1 = t53 * t4;
t35 = pkin(3) * t36;
t34 = t22 * t41 - t28 * t38;
t33 = t53 * t48;
t21 = t27 * pkin(3) + pkin(7);
t20 = -0.2e1 * t37;
t19 = 0.2e1 * t37;
t15 = 0.2e1 * (-t24 + t25) * qJD(5);
t12 = t33 * t45;
t7 = t48 * t14 - t31;
t3 = -t5 * t28 + t6 * t41;
t2 = t6 * t23 + t5 * t26;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t7 * t4 + 0.2e1 * t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t7 * t1 + 0.2e1 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t32 * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4, 0, (-t48 * t5 - t27 * t4 + (t48 * t7 + t50) * qJD(4)) * pkin(3), 0, 0, 0, 0, 0, 0, t3, t2, -t1, t5 * t22 - t21 * t1 + (t33 * t7 + t50) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t38, -0.2e1 * t35, 0, 0, t19, t15, 0, t20, 0, 0, 0.2e1 * t34, 0.2e1 * t46, 0.2e1 * t12, 0.2e1 * (t33 * t21 + t22 * t27) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3, t2, -t1, -t5 * pkin(4) - pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t35, 0, 0, t19, t15, 0, t20, 0, 0, t34 - t40, -t39 + t46, t12, (-pkin(4) * t27 + t33 * pkin(7)) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t15, 0, t20, 0, 0, -0.2e1 * t40, -0.2e1 * t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t23 + t26 * t4, t28 * t4 + t7 * t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, -t41, 0, -t21 * t23 - t26 * t35, t21 * t41 - t28 * t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, -t41, 0, -pkin(7) * t23, pkin(7) * t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
