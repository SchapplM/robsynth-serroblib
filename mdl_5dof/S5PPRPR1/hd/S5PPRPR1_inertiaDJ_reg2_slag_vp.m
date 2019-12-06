% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRPR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:20
% EndTime: 2019-12-05 15:01:22
% DurationCPUTime: 0.37s
% Computational Cost: add. (286->52), mult. (813->108), div. (0->0), fcn. (813->8), ass. (0->41)
t31 = sin(pkin(8));
t34 = sin(qJ(3));
t46 = cos(pkin(8));
t55 = cos(qJ(3));
t21 = t55 * t31 + t34 * t46;
t30 = sin(pkin(9));
t32 = cos(pkin(9));
t33 = sin(qJ(5));
t54 = cos(qJ(5));
t35 = -t33 * t30 + t54 * t32;
t6 = t35 * t21;
t19 = t34 * t31 - t55 * t46;
t50 = t33 * t32;
t20 = t54 * t30 + t50;
t14 = t20 * qJD(5);
t53 = t35 * t14;
t16 = t21 * qJD(3);
t8 = t19 * t16;
t41 = qJD(5) * t54;
t45 = qJD(5) * t33;
t13 = t30 * t45 - t32 * t41;
t52 = t20 * t13;
t48 = pkin(6) + qJ(4);
t47 = t30 ^ 2 + t32 ^ 2;
t43 = t48 * t30;
t15 = t19 * qJD(3);
t42 = t47 * t15;
t40 = t54 * qJD(4);
t39 = t47 * qJD(4);
t38 = 0.2e1 * t39;
t36 = t54 * t43;
t27 = -t32 * pkin(4) - pkin(3);
t22 = t48 * t32;
t11 = t54 * t22 - t33 * t43;
t10 = -t33 * t22 - t36;
t5 = t20 * t21;
t4 = -t22 * t41 - qJD(4) * t50 + (t48 * t45 - t40) * t30;
t3 = qJD(5) * t36 - t32 * t40 + (qJD(4) * t30 + qJD(5) * t22) * t33;
t2 = -qJD(5) * t6 + t20 * t15;
t1 = t21 * t14 + t35 * t15;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t21 * t15 + 0.2e1 * t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t21 * t42 + 0.2e1 * t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t6 * t1 - 0.2e1 * t5 * t2 + 0.2e1 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t20 - t6 * t13 + t5 * t14 + t2 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t52 - 0.2e1 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t32, t16 * t30, -t42, -t16 * pkin(3) - qJ(4) * t42 + t21 * t39, 0, 0, 0, 0, 0, 0, t19 * t14 - t16 * t35, -t19 * t13 + t16 * t20, -t1 * t35 - t5 * t13 - t6 * t14 - t2 * t20, -t1 * t11 + t2 * t10 + t16 * t27 - t6 * t3 - t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t10 - t13 * t11 - t20 * t3 + t35 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, qJ(4) * t38, -0.2e1 * t52, -0.2e1 * t13 * t35 - 0.2e1 * t20 * t14, 0, -0.2e1 * t53, 0, 0, 0.2e1 * t27 * t14, -0.2e1 * t27 * t13, 0.2e1 * t10 * t13 - 0.2e1 * t11 * t14 - 0.2e1 * t4 * t20 - 0.2e1 * t3 * t35, 0.2e1 * t10 * t4 - 0.2e1 * t11 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, -t14, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
