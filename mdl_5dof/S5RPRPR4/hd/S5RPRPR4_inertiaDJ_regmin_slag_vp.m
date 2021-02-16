% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR4
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
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:45
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:44:14
% EndTime: 2021-01-15 11:44:17
% DurationCPUTime: 0.31s
% Computational Cost: add. (447->66), mult. (1009->128), div. (0->0), fcn. (920->8), ass. (0->44)
t42 = sin(pkin(9));
t43 = cos(pkin(9));
t46 = sin(qJ(3));
t47 = cos(qJ(3));
t33 = t42 * t46 - t43 * t47;
t34 = t42 * t47 + t43 * t46;
t45 = sin(qJ(5));
t55 = cos(qJ(5));
t13 = t55 * t33 + t45 * t34;
t56 = pkin(3) * t42;
t38 = sin(pkin(8)) * pkin(1) + pkin(6);
t53 = qJ(4) + t38;
t31 = t53 * t46;
t32 = t53 * t47;
t12 = -t42 * t31 + t43 * t32;
t52 = t46 * qJD(3);
t51 = t47 * qJD(3);
t50 = 0.2e1 * t51;
t41 = pkin(3) * t52;
t40 = -cos(pkin(8)) * pkin(1) - pkin(2);
t48 = qJD(3) * t53;
t18 = t47 * qJD(4) - t46 * t48;
t19 = -t46 * qJD(4) - t47 * t48;
t7 = -t42 * t18 + t43 * t19;
t11 = -t43 * t31 - t42 * t32;
t8 = t43 * t18 + t42 * t19;
t35 = -t47 * pkin(3) + t40;
t14 = -t45 * t33 + t55 * t34;
t39 = t43 * pkin(3) + pkin(4);
t30 = -t42 * t52 + t43 * t51;
t29 = t34 * qJD(3);
t21 = (-t39 * t45 - t55 * t56) * qJD(5);
t20 = (-t55 * t39 + t45 * t56) * qJD(5);
t17 = t29 * pkin(4) + t41;
t16 = t33 * pkin(4) + t35;
t10 = -t33 * pkin(7) + t12;
t9 = -t34 * pkin(7) + t11;
t6 = -t29 * pkin(7) + t8;
t5 = -t30 * pkin(7) + t7;
t4 = t14 * qJD(5) + t55 * t29 + t45 * t30;
t3 = t13 * qJD(5) + t45 * t29 - t55 * t30;
t2 = t55 * t5 - t45 * t6 + (-t55 * t10 - t45 * t9) * qJD(5);
t1 = -t55 * t6 - t45 * t5 + (t10 * t45 - t55 * t9) * qJD(5);
t15 = [0, 0, 0, 0, t46 * t50, 0.2e1 * (-t46 ^ 2 + t47 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t40 * t52, t40 * t50, 0.2e1 * t35 * t29 + 0.2e1 * t33 * t41, 0.2e1 * t35 * t30 + 0.2e1 * t34 * t41, -0.2e1 * t11 * t30 - 0.2e1 * t12 * t29 - 0.2e1 * t8 * t33 - 0.2e1 * t7 * t34, 0.2e1 * t11 * t7 + 0.2e1 * t12 * t8 + 0.2e1 * t35 * t41, -0.2e1 * t14 * t3, 0.2e1 * t3 * t13 - 0.2e1 * t14 * t4, 0, 0, 0, 0.2e1 * t17 * t13 + 0.2e1 * t16 * t4, 0.2e1 * t17 * t14 - 0.2e1 * t16 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t29 + t12 * t30 - t7 * t33 + t8 * t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t33 * t29 + 0.2e1 * t34 * t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t51, -t52, 0, -t38 * t51, t38 * t52, t7, -t8, (-t29 * t42 - t30 * t43) * pkin(3), (t42 * t8 + t43 * t7) * pkin(3), 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t51, -t29, -t30, 0, (-t29 * t43 + t30 * t42) * pkin(3), 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t21, 0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t30, 0, t41, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
