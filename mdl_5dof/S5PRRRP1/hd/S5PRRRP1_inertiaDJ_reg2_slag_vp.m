% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:13
% EndTime: 2019-12-05 16:40:17
% DurationCPUTime: 0.49s
% Computational Cost: add. (214->75), mult. (596->123), div. (0->0), fcn. (354->4), ass. (0->49)
t40 = cos(qJ(3));
t55 = t40 * pkin(2);
t37 = sin(qJ(4));
t31 = t37 * qJD(4);
t30 = pkin(4) * t31;
t38 = sin(qJ(3));
t50 = pkin(2) * qJD(3);
t46 = t38 * t50;
t14 = t30 + t46;
t39 = cos(qJ(4));
t28 = -t39 * pkin(4) - pkin(3);
t19 = t28 - t55;
t33 = t39 * qJD(4);
t54 = t14 * t37 + t19 * t33;
t53 = -qJ(5) - pkin(7);
t27 = -pkin(3) - t55;
t52 = t27 * t33 + t37 * t46;
t35 = t37 ^ 2;
t51 = qJD(4) * t35 * pkin(4) + t28 * t33;
t26 = t38 * pkin(2) + pkin(7);
t49 = qJ(5) + t26;
t48 = pkin(3) * t31;
t47 = pkin(3) * t33;
t45 = t40 * t50;
t44 = pkin(4) * t33;
t9 = t19 * t31;
t17 = t28 * t31;
t43 = t37 * t33;
t36 = t39 ^ 2;
t42 = (t35 + t36) * t40;
t41 = t27 * t31 - t39 * t46;
t34 = t39 * qJ(5);
t32 = t39 * qJD(5);
t25 = -0.2e1 * t43;
t24 = 0.2e1 * t43;
t22 = t39 * t45;
t21 = t39 * pkin(7) + t34;
t20 = t53 * t37;
t13 = 0.2e1 * (-t35 + t36) * qJD(4);
t12 = t39 * t26 + t34;
t11 = t49 * t37;
t7 = -t37 * qJD(5) + t53 * t33;
t6 = -t53 * t31 - t32;
t5 = t42 * t50;
t4 = t6 * t39;
t3 = (-qJD(5) - t45) * t37 - t49 * t33;
t2 = t49 * t31 - t22 - t32;
t1 = t2 * t39;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t2 + t39 * t3 + (t11 * t37 + t12 * t39) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t46, -0.2e1 * t45, 0, 0, t24, t13, 0, t25, 0, 0, 0.2e1 * t41, 0.2e1 * t52, 0.2e1 * t5, 0.2e1 * (t26 * t42 + t27 * t38) * t50, t24, t13, 0, t25, 0, 0, -0.2e1 * t14 * t39 + 0.2e1 * t9, 0.2e1 * t54, -0.2e1 * t3 * t37 - 0.2e1 * t1 + 0.2e1 * (t11 * t39 - t12 * t37) * qJD(4), -0.2e1 * t11 * t3 - 0.2e1 * t12 * t2 + 0.2e1 * t19 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t6 + t39 * t7 + (-t20 * t37 + t21 * t39) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t45, 0, 0, t24, t13, 0, t25, 0, 0, t41 - t48, -t47 + t52, t5, (-pkin(3) * t38 + pkin(7) * t42) * t50, t24, t13, 0, t25, 0, 0, t17 + t9 + (-t14 - t30) * t39, t51 + t54, -t1 - t4 + (-t3 - t7) * t37 + ((t11 - t20) * t39 + (-t12 - t21) * t37) * qJD(4), pkin(4) * t9 - t11 * t7 - t12 * t6 + t14 * t28 - t2 * t21 + t3 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t13, 0, t25, 0, 0, -0.2e1 * t48, -0.2e1 * t47, 0, 0, t24, t13, 0, t25, 0, 0, -0.2e1 * pkin(4) * t43 + 0.2e1 * t17, 0.2e1 * t51, -0.2e1 * t7 * t37 - 0.2e1 * t4 + 0.2e1 * (-t20 * t39 - t21 * t37) * qJD(4), 0.2e1 * pkin(4) * t17 + 0.2e1 * t20 * t7 - 0.2e1 * t21 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t33, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t33, 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t31, 0, -t26 * t33 - t37 * t45, t26 * t31 - t22, 0, 0, 0, 0, t33, 0, -t31, 0, t3, t2, -t44, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t31, 0, -pkin(7) * t33, pkin(7) * t31, 0, 0, 0, 0, t33, 0, -t31, 0, t7, t6, -t44, t7 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t33, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t33, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
