% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:54
% EndTime: 2019-12-31 18:12:55
% DurationCPUTime: 0.44s
% Computational Cost: add. (402->61), mult. (990->112), div. (0->0), fcn. (876->4), ass. (0->43)
t29 = sin(pkin(7));
t47 = pkin(6) + qJ(2);
t21 = t47 * t29;
t30 = cos(pkin(7));
t22 = t47 * t30;
t32 = sin(qJ(3));
t50 = cos(qJ(3));
t43 = qJD(3) * t50;
t44 = t50 * t30;
t6 = (qJD(2) * t29 + qJD(3) * t22) * t32 - qJD(2) * t44 + t21 * t43;
t49 = t32 * t29;
t15 = qJD(3) * t49 - t30 * t43;
t52 = -0.2e1 * t15;
t20 = t50 * t29 + t32 * t30;
t16 = t20 * qJD(3);
t51 = 0.2e1 * t16;
t48 = pkin(3) + qJ(5);
t45 = (qJ(4) * qJD(4));
t19 = -t44 + t49;
t11 = t19 * t51;
t12 = t20 * t52;
t26 = -t30 * pkin(2) - pkin(1);
t13 = t50 * t21 + t32 * t22;
t42 = 0.2e1 * (t29 ^ 2 + t30 ^ 2) * qJD(2);
t14 = -t32 * t21 + t50 * t22;
t7 = t20 * qJD(2) + t14 * qJD(3);
t41 = t13 * t7 - t14 * t6;
t40 = t15 * t19 - t16 * t20;
t39 = t15 * qJ(4) - t20 * qJD(4);
t38 = -qJ(4) * t16 - qJD(4) * t19;
t36 = -t20 * qJ(4) + t26;
t35 = 0.2e1 * t40;
t34 = -0.2e1 * t13 * t15 - 0.2e1 * t14 * t16 + 0.2e1 * t6 * t19 + 0.2e1 * t7 * t20;
t33 = 0.2e1 * qJD(4);
t10 = t19 * pkin(3) + t36;
t9 = -t19 * pkin(4) + t14;
t8 = t20 * pkin(4) + t13;
t5 = t48 * t19 + t36;
t4 = t16 * pkin(3) + t39;
t3 = -t15 * pkin(4) + t7;
t2 = -t16 * pkin(4) - t6;
t1 = t19 * qJD(5) + t48 * t16 + t39;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, qJ(2) * t42, t12, t35, 0, t11, 0, 0, t26 * t51, t26 * t52, t34, 0.2e1 * t41, 0, 0, 0, t12, t35, t11, t34, -0.2e1 * t10 * t16 - 0.2e1 * t4 * t19, 0.2e1 * t10 * t15 - 0.2e1 * t4 * t20, 0.2e1 * t10 * t4 + 0.2e1 * t41, 0, 0, 0, t11, -0.2e1 * t40, t12, -0.2e1 * t8 * t15 - 0.2e1 * t9 * t16 - 0.2e1 * t2 * t19 + 0.2e1 * t3 * t20, -0.2e1 * t1 * t20 + 0.2e1 * t5 * t15, 0.2e1 * t1 * t19 + 0.2e1 * t5 * t16, 0.2e1 * t5 * t1 + 0.2e1 * t9 * t2 + 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, t4, 0, 0, 0, 0, 0, 0, 0, t15, t16, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, -t16, 0, -t7, t6, 0, 0, 0, t15, t16, 0, 0, 0, pkin(3) * t15 + t38, t7, -t6, -t7 * pkin(3) - t6 * qJ(4) + t14 * qJD(4), 0, t16, -t15, 0, 0, 0, -qJD(5) * t20 + t15 * t48 + t38, t2, -t3, t2 * qJ(4) + t9 * qJD(4) - t8 * qJD(5) - t3 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 2 * t45, 0, 0, 0, 0, 0, 0, 0, t33, 0.2e1 * qJD(5), 0.2e1 * qJD(5) * t48 + (2 * t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, t7, 0, 0, 0, 0, 0, 0, -t15, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
