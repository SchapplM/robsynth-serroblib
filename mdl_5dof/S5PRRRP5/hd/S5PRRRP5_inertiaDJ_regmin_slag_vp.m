% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRP5
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
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:18
% EndTime: 2019-12-05 16:49:20
% DurationCPUTime: 0.36s
% Computational Cost: add. (391->71), mult. (1045->141), div. (0->0), fcn. (921->6), ass. (0->47)
t31 = cos(qJ(4));
t32 = cos(qJ(3));
t28 = sin(qJ(4));
t29 = sin(qJ(3));
t49 = t28 * t29;
t16 = -t31 * t32 + t49;
t30 = sin(qJ(2));
t14 = t16 * t30;
t52 = qJD(3) + qJD(4);
t51 = pkin(6) + pkin(7);
t37 = qJD(3) * t51;
t18 = t29 * t37;
t19 = t32 * t37;
t21 = t51 * t29;
t22 = t51 * t32;
t35 = -t31 * t21 - t28 * t22;
t5 = -t35 * qJD(4) + t31 * t18 + t28 * t19;
t50 = t31 * pkin(3);
t48 = qJD(4) * t31;
t47 = t29 * qJD(3);
t46 = t30 * qJD(2);
t45 = t32 * qJD(3);
t33 = cos(qJ(2));
t44 = t33 * qJD(2);
t43 = -0.2e1 * pkin(2) * qJD(3);
t42 = pkin(3) * t47;
t41 = qJD(4) * t28 * pkin(3);
t40 = pkin(3) * t48;
t39 = t29 * t44;
t38 = t32 * t44;
t27 = -t32 * pkin(3) - pkin(2);
t34 = t28 * t21 - t31 * t22;
t17 = t28 * t32 + t31 * t29;
t6 = t34 * qJD(4) + t28 * t18 - t31 * t19;
t11 = t52 * t17;
t26 = pkin(4) + t50;
t13 = t17 * t30;
t12 = t16 * pkin(4) + t27;
t10 = -t31 * t45 - t32 * t48 + t52 * t49;
t9 = t11 * pkin(4) + t42;
t8 = -t16 * qJ(5) - t34;
t7 = -t17 * qJ(5) + t35;
t4 = t52 * t14 - t17 * t44;
t3 = t11 * t30 + t28 * t39 - t31 * t38;
t2 = t10 * qJ(5) - t17 * qJD(5) + t6;
t1 = -t11 * qJ(5) - t16 * qJD(5) - t5;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t13 * t4 + 0.2e1 * t14 * t3 - 0.2e1 * t30 * t44; 0, 0, -t46, -t44, 0, 0, 0, 0, 0, -t32 * t46 - t33 * t47, t29 * t46 - t33 * t45, 0, 0, 0, 0, 0, -t33 * t11 + t16 * t46, t33 * t10 + t17 * t46, -t13 * t10 + t14 * t11 + t3 * t16 - t4 * t17, -t14 * t1 + t12 * t46 - t13 * t2 - t3 * t8 - t33 * t9 + t4 * t7; 0, 0, 0, 0, 0.2e1 * t29 * t45, 0.2e1 * (-t29 ^ 2 + t32 ^ 2) * qJD(3), 0, 0, 0, t29 * t43, t32 * t43, -0.2e1 * t17 * t10, 0.2e1 * t10 * t16 - 0.2e1 * t17 * t11, 0, 0, 0, 0.2e1 * t27 * t11 + 0.2e1 * t16 * t42, -0.2e1 * t27 * t10 + 0.2e1 * t17 * t42, -0.2e1 * t1 * t16 + 0.2e1 * t7 * t10 - 0.2e1 * t8 * t11 - 0.2e1 * t2 * t17, 0.2e1 * t8 * t1 + 0.2e1 * t12 * t9 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t45 - t39, t30 * t47 - t38, 0, 0, 0, 0, 0, t4, t3, 0, t4 * t26 + (-t28 * t3 + (t13 * t28 - t14 * t31) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, t45, -t47, 0, -pkin(6) * t45, pkin(6) * t47, 0, 0, -t10, -t11, 0, t6, t5, t26 * t10 + (-t11 * t28 + (-t16 * t31 + t17 * t28) * qJD(4)) * pkin(3), t2 * t26 + (t1 * t28 + (-t28 * t7 + t31 * t8) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t41, -0.2e1 * t40, 0, 0.2e1 * (-t26 + t50) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, t6, t5, pkin(4) * t10, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t40, 0, -pkin(4) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
