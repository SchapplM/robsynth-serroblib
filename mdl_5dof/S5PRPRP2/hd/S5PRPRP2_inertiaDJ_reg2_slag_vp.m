% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:05
% EndTime: 2019-12-05 15:31:07
% DurationCPUTime: 0.48s
% Computational Cost: add. (243->65), mult. (647->120), div. (0->0), fcn. (488->4), ass. (0->51)
t29 = sin(qJ(4));
t30 = cos(qJ(4));
t28 = cos(pkin(8));
t44 = qJ(3) * qJD(4);
t38 = t29 * t44;
t27 = sin(pkin(8));
t39 = t28 * pkin(3) + pkin(2);
t31 = -t27 * pkin(6) - t39;
t13 = t30 * t31;
t47 = t28 * qJD(3);
t53 = -qJD(4) * t13 - t30 * t47;
t4 = t28 * t38 + t53;
t46 = t29 * qJD(3);
t40 = t28 * t46;
t23 = t30 * t28 * qJ(3);
t8 = t29 * t31 + t23;
t5 = -qJD(4) * t8 - t40;
t52 = qJ(3) * t29;
t43 = t28 * t52;
t7 = t13 - t43;
t57 = t4 * t29 - t5 * t30 + (t29 * t7 - t30 * t8) * qJD(4);
t51 = t27 * qJ(5);
t42 = t30 * t51;
t48 = t27 * qJD(5);
t1 = t29 * t48 + (t42 + t43) * qJD(4) + t53;
t2 = -t40 - t30 * t48 + (-t23 + ((pkin(6) + qJ(5)) * t27 + t39) * t29) * qJD(4);
t3 = -t42 + t13 + (-pkin(4) - t52) * t28;
t6 = -t29 * t51 + t8;
t56 = t1 * t29 - t2 * t30 + (t29 * t3 - t30 * t6) * qJD(4);
t55 = 0.2e1 * t27;
t54 = 0.2e1 * qJD(4);
t50 = qJD(4) * t29;
t49 = qJD(4) * t30;
t45 = qJ(3) * qJD(3);
t20 = t27 * t50;
t41 = t27 * t49;
t37 = t27 * t28 * t54;
t25 = t27 ^ 2;
t36 = t25 * t29 * t49;
t26 = t28 ^ 2;
t24 = t25 * t45;
t22 = t28 * t49;
t21 = t28 * t50;
t18 = -0.2e1 * t36;
t17 = 0.2e1 * t36;
t16 = t30 * t37;
t15 = t29 * t37;
t14 = (pkin(4) * t29 + qJ(3)) * t27;
t11 = (pkin(4) * t49 + qJD(3)) * t27;
t9 = (t29 ^ 2 - t30 ^ 2) * t25 * t54;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t47 - t29 * t5 - t30 * t4 + (-t29 * t8 - t30 * t7) * qJD(4)) * t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t11 + (-t1 * t30 - t2 * t29 + (-t29 * t6 - t3 * t30) * qJD(4)) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t25 + t26) * qJD(3), 0.2e1 * t26 * t45 + 0.2e1 * t24, t18, t9, t15, t17, t16, 0, -0.2e1 * t5 * t28 + 0.2e1 * (t30 * t44 + t46) * t25, -0.2e1 * t4 * t28 + 0.2e1 * (qJD(3) * t30 - t38) * t25, t57 * t55, -0.2e1 * t8 * t4 + 0.2e1 * t7 * t5 + 0.2e1 * t24, t18, t9, t15, t17, t16, 0, -0.2e1 * t2 * t28 + 0.2e1 * (t11 * t29 + t14 * t49) * t27, -0.2e1 * t1 * t28 + 0.2e1 * (t11 * t30 - t14 * t50) * t27, t56 * t55, -0.2e1 * t6 * t1 + 0.2e1 * t14 * t11 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22, 0, -t57, 0, 0, 0, 0, 0, 0, t21, t22, 0, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t20, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t20, 0, -pkin(4) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, -t41, 0, t5, t4, 0, 0, 0, 0, -t20, 0, -t41, 0, t2, t1, pkin(4) * t20, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t49, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t49, 0, -pkin(4) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t20, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
