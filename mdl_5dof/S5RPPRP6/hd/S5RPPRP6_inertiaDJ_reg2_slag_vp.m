% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRP6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:17
% EndTime: 2019-12-31 17:55:18
% DurationCPUTime: 0.38s
% Computational Cost: add. (366->52), mult. (773->86), div. (0->0), fcn. (690->4), ass. (0->38)
t30 = sin(pkin(7));
t31 = cos(pkin(7));
t52 = sin(qJ(4));
t42 = qJD(4) * t52;
t32 = cos(qJ(4));
t47 = qJD(4) * t32;
t17 = -t30 * t42 + t31 * t47;
t19 = t32 * t30 + t52 * t31;
t14 = t19 * t17;
t16 = t19 * qJD(4);
t20 = -t52 * t30 + t32 * t31;
t50 = t20 * t16;
t59 = 0.2e1 * t50 - 0.2e1 * t14;
t48 = pkin(1) + qJ(3);
t45 = -pkin(6) - t48;
t21 = t45 * t30;
t41 = t45 * t31;
t38 = t32 * t41;
t11 = t52 * t21 - t38;
t35 = t52 * t41;
t12 = t32 * t21 + t35;
t6 = t19 * qJD(3) - qJD(4) * t38 + t21 * t42;
t7 = t20 * qJD(3) + qJD(4) * t35 + t21 * t47;
t57 = -t11 * t16 - t12 * t17 + t6 * t19 + t7 * t20;
t22 = (t30 ^ 2 + t31 ^ 2) * qJD(3);
t56 = 0.2e1 * t14;
t55 = 2 * qJD(2);
t54 = 2 * qJD(5);
t25 = t30 * pkin(3) + qJ(2);
t46 = qJ(2) * qJD(2);
t44 = t11 * t7 - t12 * t6;
t39 = t16 * t19 - t20 * t17;
t34 = t16 * pkin(4) - t17 * qJ(5) - t19 * qJD(5);
t33 = 0.2e1 * t57;
t13 = -0.2e1 * t50;
t10 = t19 * pkin(4) - t20 * qJ(5) + t25;
t5 = t17 * pkin(4) + t16 * qJ(5) - t20 * qJD(5) + qJD(2);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0.2e1 * t46, 0, 0, 0, 0, 0, 0, t30 * t55, t31 * t55, 0.2e1 * t22, 0.2e1 * t48 * t22 + 0.2e1 * t46, t13, 0.2e1 * t39, 0, t56, 0, 0, 0.2e1 * qJD(2) * t19 + 0.2e1 * t25 * t17, 0.2e1 * qJD(2) * t20 - 0.2e1 * t25 * t16, t33, 0.2e1 * t25 * qJD(2) + 0.2e1 * t44, t13, 0, -0.2e1 * t39, 0, 0, t56, 0.2e1 * t10 * t17 + 0.2e1 * t5 * t19, t33, 0.2e1 * t10 * t16 - 0.2e1 * t5 * t20, 0.2e1 * t10 * t5 + 0.2e1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t57, 0, 0, 0, 0, 0, 0, 0, t59, 0, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, t17, -t16, 0, qJD(2), 0, 0, 0, 0, 0, 0, t17, 0, t16, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, -t17, 0, -t7, t6, 0, 0, 0, -t16, 0, 0, t17, 0, -t7, t34, -t6, -t7 * pkin(4) - t6 * qJ(5) + t12 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, t17, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, qJ(5) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
