% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:01
% EndTime: 2020-01-03 11:34:04
% DurationCPUTime: 0.27s
% Computational Cost: add. (195->41), mult. (501->71), div. (0->0), fcn. (388->8), ass. (0->40)
t31 = sin(pkin(9));
t33 = cos(pkin(9));
t43 = t31 ^ 2 + t33 ^ 2;
t58 = t43 * qJD(4);
t59 = 0.2e1 * t58;
t24 = cos(pkin(8)) * pkin(1) + pkin(2);
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t54 = pkin(1) * sin(pkin(8));
t38 = -t37 * t24 + t35 * t54;
t13 = t38 * qJD(3);
t12 = qJD(4) - t13;
t57 = t43 * t12;
t34 = sin(qJ(5));
t36 = cos(qJ(5));
t18 = t34 * t31 - t36 * t33;
t40 = -pkin(3) + t38;
t53 = t33 * pkin(4);
t11 = t40 - t53;
t39 = t35 * t24 + t37 * t54;
t14 = t39 * qJD(3);
t19 = t36 * t31 + t34 * t33;
t17 = t19 * qJD(5);
t56 = t11 * t17 + t14 * t18;
t16 = t18 * qJD(5);
t55 = -t11 * t16 + t14 * t19;
t51 = t14 * t31;
t50 = t14 * t33;
t25 = -pkin(3) - t53;
t49 = t25 * t16;
t48 = t25 * t17;
t28 = t33 * pkin(7);
t21 = t33 * qJ(4) + t28;
t20 = (-pkin(7) - qJ(4)) * t31;
t15 = qJ(4) + t39;
t8 = t33 * t15 + t28;
t7 = (-pkin(7) - t15) * t31;
t6 = -0.2e1 * t19 * t16;
t1 = 0.2e1 * t16 * t18 - 0.2e1 * t19 * t17;
t2 = [0, 0, 0, 0, 0, -0.2e1 * t14, 0.2e1 * t13, -0.2e1 * t50, 0.2e1 * t51, 0.2e1 * t57, 0.2e1 * t40 * t14 + 0.2e1 * t15 * t57, t6, t1, 0, 0, 0, 0.2e1 * t56, 0.2e1 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t14, t13, -t50, t51, t58 + t57, -t14 * pkin(3) + qJ(4) * t57 + t15 * t58, t6, t1, 0, 0, 0, t48 + t56, -t49 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, qJ(4) * t59, t6, t1, 0, 0, 0, 0.2e1 * t48, -0.2e1 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, t17, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, -t19 * t12 + (-t34 * t7 - t36 * t8) * qJD(5), t18 * t12 + (t34 * t8 - t36 * t7) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, (-t20 * t34 - t21 * t36) * qJD(5) - t19 * qJD(4), (-t20 * t36 + t21 * t34) * qJD(5) + t18 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
