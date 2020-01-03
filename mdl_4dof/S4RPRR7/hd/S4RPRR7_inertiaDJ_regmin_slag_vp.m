% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:04
% EndTime: 2019-12-31 16:54:06
% DurationCPUTime: 0.26s
% Computational Cost: add. (257->58), mult. (692->124), div. (0->0), fcn. (627->6), ass. (0->51)
t29 = cos(qJ(4));
t24 = t29 ^ 2;
t27 = sin(qJ(4));
t45 = t27 ^ 2 - t24;
t39 = t45 * qJD(4);
t26 = cos(pkin(7));
t20 = -t26 * pkin(2) - pkin(1);
t56 = 0.2e1 * t20;
t25 = sin(pkin(7));
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t15 = t30 * t25 + t28 * t26;
t46 = pkin(5) + qJ(2);
t16 = t46 * t25;
t17 = t46 * t26;
t9 = -t28 * t16 + t30 * t17;
t4 = t15 * qJD(2) + t9 * qJD(3);
t55 = t4 * t27;
t54 = t4 * t29;
t47 = t30 * t26;
t14 = t28 * t25 - t47;
t10 = t14 * qJD(3);
t53 = t15 * t10;
t52 = t27 * t10;
t11 = t15 * qJD(3);
t51 = t27 * t11;
t50 = t29 * t10;
t49 = t29 * t11;
t48 = t30 * t16;
t44 = qJD(4) * t27;
t43 = qJD(4) * t29;
t42 = -0.2e1 * pkin(3) * qJD(4);
t41 = t27 * t43;
t40 = -0.4e1 * t15 * t27 * t29;
t38 = 0.2e1 * (t25 ^ 2 + t26 ^ 2) * qJD(2);
t37 = pkin(3) * t10 - pkin(6) * t11;
t36 = pkin(3) * t15 + pkin(6) * t14;
t7 = t14 * pkin(3) - t15 * pkin(6) + t20;
t35 = t27 * t9 - t29 * t7;
t34 = t27 * t7 + t29 * t9;
t33 = -t15 * t43 + t52;
t32 = -t15 * t44 - t50;
t31 = t14 * t43 + t51;
t13 = t15 ^ 2;
t8 = t28 * t17 + t48;
t6 = t11 * pkin(3) + t10 * pkin(6);
t5 = -t14 * t44 + t49;
t3 = qJD(3) * t48 - qJD(2) * t47 + (qJD(2) * t25 + qJD(3) * t17) * t28;
t2 = -t34 * qJD(4) + t27 * t3 + t29 * t6;
t1 = t35 * qJD(4) - t27 * t6 + t29 * t3;
t12 = [0, 0, 0, 0, 0, t38, qJ(2) * t38, -0.2e1 * t53, 0.2e1 * t10 * t14 - 0.2e1 * t15 * t11, 0, 0, 0, t11 * t56, -t10 * t56, -0.2e1 * t13 * t41 - 0.2e1 * t24 * t53, -t10 * t40 + 0.2e1 * t13 * t39, 0.2e1 * t32 * t14 + 0.2e1 * t15 * t49, 0.2e1 * t33 * t14 - 0.2e1 * t15 * t51, 0.2e1 * t14 * t11, 0.2e1 * t2 * t14 - 0.2e1 * t35 * t11 - 0.2e1 * t8 * t52 + 0.2e1 * (t8 * t43 + t55) * t15, 0.2e1 * t1 * t14 - 0.2e1 * t34 * t11 - 0.2e1 * t8 * t50 + 0.2e1 * (-t8 * t44 + t54) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, 0, 0, 0, 0, 0, t5, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, -t4, t3, -t15 * t39 - t27 * t50, qJD(4) * t40 + t45 * t10, t31, t5, 0, -t54 + t37 * t27 + (t27 * t8 - t36 * t29) * qJD(4), t55 + t37 * t29 + (t36 * t27 + t29 * t8) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t41, -0.2e1 * t39, 0, 0, 0, t27 * t42, t29 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t33, t11, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t44, 0, -pkin(6) * t43, pkin(6) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
