% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPR3
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:41
% EndTime: 2019-12-05 16:19:42
% DurationCPUTime: 0.28s
% Computational Cost: add. (352->59), mult. (870->118), div. (0->0), fcn. (795->6), ass. (0->42)
t40 = sin(pkin(9));
t41 = cos(pkin(9));
t43 = sin(qJ(3));
t44 = cos(qJ(3));
t30 = -t40 * t43 + t41 * t44;
t31 = t40 * t44 + t41 * t43;
t42 = sin(qJ(5));
t53 = cos(qJ(5));
t11 = -t53 * t30 + t42 * t31;
t54 = pkin(3) * t40;
t51 = -qJ(4) - pkin(6);
t45 = qJD(3) * t51;
t26 = t44 * qJD(4) + t43 * t45;
t27 = -t43 * qJD(4) + t44 * t45;
t8 = t41 * t26 + t40 * t27;
t35 = t51 * t43;
t36 = t51 * t44;
t14 = t40 * t35 - t41 * t36;
t50 = t43 * qJD(3);
t49 = t44 * qJD(3);
t48 = -0.2e1 * pkin(2) * qJD(3);
t39 = pkin(3) * t50;
t47 = -t44 * pkin(3) - pkin(2);
t7 = -t40 * t26 + t41 * t27;
t13 = t41 * t35 + t40 * t36;
t12 = t42 * t30 + t53 * t31;
t38 = t41 * pkin(3) + pkin(4);
t29 = -t40 * t50 + t41 * t49;
t28 = t31 * qJD(3);
t21 = (-t38 * t42 - t53 * t54) * qJD(5);
t20 = (-t53 * t38 + t42 * t54) * qJD(5);
t16 = -t30 * pkin(4) + t47;
t15 = t28 * pkin(4) + t39;
t10 = t30 * pkin(7) + t14;
t9 = -t31 * pkin(7) + t13;
t6 = -t28 * pkin(7) + t8;
t5 = -t29 * pkin(7) + t7;
t4 = t12 * qJD(5) + t53 * t28 + t42 * t29;
t3 = t11 * qJD(5) + t42 * t28 - t53 * t29;
t2 = t53 * t5 - t42 * t6 + (-t53 * t10 - t42 * t9) * qJD(5);
t1 = -t53 * t6 - t42 * t5 + (t10 * t42 - t53 * t9) * qJD(5);
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t30 * t28 + 0.2e1 * t31 * t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t13 + t29 * t14 + t30 * t7 + t31 * t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0.2e1 * t43 * t49, 0.2e1 * (-t43 ^ 2 + t44 ^ 2) * qJD(3), 0, 0, 0, t43 * t48, t44 * t48, -0.2e1 * t13 * t29 - 0.2e1 * t14 * t28 + 0.2e1 * t8 * t30 - 0.2e1 * t7 * t31, 0.2e1 * t13 * t7 + 0.2e1 * t14 * t8 + 0.2e1 * t47 * t39, -0.2e1 * t12 * t3, 0.2e1 * t3 * t11 - 0.2e1 * t12 * t4, 0, 0, 0, 0.2e1 * t15 * t11 + 0.2e1 * t16 * t4, 0.2e1 * t15 * t12 - 0.2e1 * t16 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t49, 0, (-t28 * t41 + t29 * t40) * pkin(3), 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, t49, -t50, 0, -pkin(6) * t49, pkin(6) * t50, (-t28 * t40 - t29 * t41) * pkin(3), (t40 * t8 + t41 * t7) * pkin(3), 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t21, 0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
