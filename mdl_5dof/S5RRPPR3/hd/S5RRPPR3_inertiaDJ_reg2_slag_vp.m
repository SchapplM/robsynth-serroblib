% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:36
% EndTime: 2019-12-31 19:26:37
% DurationCPUTime: 0.35s
% Computational Cost: add. (148->46), mult. (441->74), div. (0->0), fcn. (277->6), ass. (0->40)
t29 = 2 * qJD(4);
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t28 = cos(qJ(2));
t41 = pkin(1) * qJD(2);
t36 = t28 * t41;
t26 = sin(qJ(2));
t37 = t26 * t41;
t9 = -t23 * t37 + t24 * t36;
t5 = qJD(4) + t9;
t18 = t28 * pkin(1) + pkin(2);
t43 = t24 * t26;
t30 = pkin(1) * t43 + t23 * t18;
t7 = qJ(4) + t30;
t45 = t7 * t5;
t25 = sin(qJ(5));
t27 = cos(qJ(5));
t39 = t27 * qJD(5);
t44 = t5 * t25 + t7 * t39;
t17 = t23 * pkin(2) + qJ(4);
t42 = qJD(4) * t25 + t17 * t39;
t40 = t25 * qJD(5);
t38 = t17 * t29;
t35 = t25 * t39;
t34 = -t24 * pkin(2) - pkin(3);
t21 = t25 ^ 2;
t22 = t27 ^ 2;
t8 = (t23 * t28 + t43) * t41;
t1 = (t21 + t22) * t8;
t33 = -t23 * t26 * pkin(1) + t24 * t18;
t32 = -pkin(3) - t33;
t31 = t7 * qJD(4) + t5 * t17;
t20 = qJD(4) * t27;
t16 = -pkin(7) + t34;
t14 = -0.2e1 * t35;
t13 = 0.2e1 * t35;
t10 = 0.2e1 * (t21 - t22) * qJD(5);
t6 = -pkin(7) + t32;
t4 = t5 * t27;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t37, -0.2e1 * t36, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t8, -0.2e1 * t9, 0, 0.2e1 * t30 * t9 - 0.2e1 * t33 * t8, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t8, 0.2e1 * t5, 0.2e1 * t32 * t8 + 0.2e1 * t45, t14, t10, 0, t13, 0, 0, 0.2e1 * t44, -0.2e1 * t7 * t40 + 0.2e1 * t4, -0.2e1 * t1, 0.2e1 * t6 * t1 + 0.2e1 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t36, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, 0, (t23 * t9 - t24 * t8) * pkin(2), 0, 0, 0, 0, 0, 0, 0, t8, t29 + t9, t8 * t34 + t31, t14, t10, 0, t13, 0, 0, t42 + t44, t20 + t4 + (-t17 - t7) * t40, -t1, t16 * t1 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t38, t14, t10, 0, t13, 0, 0, 0.2e1 * t42, -0.2e1 * t17 * t40 + 0.2e1 * t20, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, 0, -t39, 0, t27 * t8 - t6 * t40, -t25 * t8 - t6 * t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, 0, -t39, 0, -t16 * t40, -t16 * t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
