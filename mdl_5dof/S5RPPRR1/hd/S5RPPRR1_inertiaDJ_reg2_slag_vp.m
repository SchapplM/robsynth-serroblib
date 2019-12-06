% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:15
% EndTime: 2019-12-05 17:38:17
% DurationCPUTime: 0.44s
% Computational Cost: add. (367->63), mult. (724->104), div. (0->0), fcn. (585->4), ass. (0->46)
t23 = sin(qJ(4));
t24 = cos(qJ(5));
t25 = cos(qJ(4));
t48 = sin(qJ(5));
t36 = t48 * t25;
t10 = t24 * t23 + t36;
t37 = t48 * t23;
t47 = t24 * t25;
t11 = -t37 + t47;
t53 = qJD(4) + qJD(5);
t3 = t53 * t10;
t35 = qJD(5) * t48;
t4 = -qJD(4) * t37 - t23 * t35 + t53 * t47;
t54 = (qJD(5) * (-t10 * t24 + t48 * t11) + t24 * t3 - t48 * t4) * pkin(4);
t51 = t10 * t4;
t50 = t11 * t3;
t21 = -pkin(6) + qJ(2);
t49 = pkin(7) - t21;
t22 = pkin(1) + qJ(3);
t46 = (t22 * qJD(3));
t45 = t23 * qJD(2);
t44 = t23 * qJD(4);
t43 = t25 * qJD(2);
t42 = t25 * qJD(4);
t41 = qJ(2) * qJD(2);
t40 = qJD(5) * t24 * pkin(4);
t39 = t23 * t42;
t38 = t49 * t25;
t19 = t23 ^ 2;
t20 = t25 ^ 2;
t13 = (t19 + t20) * qJD(2);
t34 = t24 * t38;
t33 = pkin(4) * t35;
t32 = -t50 + t51;
t29 = -qJD(4) * t38 + t45;
t28 = t49 * t44 + t43;
t12 = t49 * t23;
t6 = -t24 * t12 - t49 * t36;
t1 = qJD(5) * t34 - t12 * t35 - t24 * t29 - t48 * t28;
t2 = -qJD(5) * t6 + t24 * t28 - t48 * t29;
t5 = t48 * t12 - t34;
t27 = t1 * t10 - t2 * t11 + t5 * t3 - t6 * t4;
t26 = 0.2e1 * qJD(2);
t17 = t23 * pkin(4) + t22;
t14 = pkin(4) * t42 + qJD(3);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0.2e1 * t41, 0, 0, 0, 0, 0, 0, 0, t26, 0.2e1 * qJD(3), 0.2e1 * t41 + (2 * t46), -0.2e1 * t39, 0.2e1 * (t19 - t20) * qJD(4), 0, 0.2e1 * t39, 0, 0, 0.2e1 * qJD(3) * t23 + 0.2e1 * t22 * t42, 0.2e1 * qJD(3) * t25 - 0.2e1 * t22 * t44, -0.2e1 * t13, 0.2e1 * t13 * t21 + (2 * t46), -0.2e1 * t50, 0.2e1 * t3 * t10 - 0.2e1 * t11 * t4, 0, 0.2e1 * t51, 0, 0, 0.2e1 * t14 * t10 + 0.2e1 * t17 * t4, 0.2e1 * t14 * t11 - 0.2e1 * t17 * t3, 0.2e1 * t27, -0.2e1 * t6 * t1 + 0.2e1 * t17 * t14 + 0.2e1 * t5 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, 0, 0, 0, 0, 0, -t42, t44, 0, -qJD(3), 0, 0, 0, 0, 0, 0, -t4, t3, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t32, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, 0, -t42, 0, -t21 * t44 + t43, -t21 * t42 - t45, 0, 0, 0, 0, -t3, 0, -t4, 0, t2, t1, t54, (-t48 * t1 + t2 * t24 + (t24 * t6 - t48 * t5) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t42, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t33, -0.2e1 * t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, -t4, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
