% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR3
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
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:55
% EndTime: 2022-01-23 09:20:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (242->50), mult. (621->95), div. (0->0), fcn. (425->8), ass. (0->54)
t28 = cos(pkin(8)) * pkin(1) + pkin(2);
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t65 = pkin(1) * sin(pkin(8));
t40 = -t39 * t28 + t37 * t65;
t14 = t40 * qJD(3);
t13 = qJD(4) - t14;
t33 = sin(pkin(9));
t31 = t33 ^ 2;
t11 = t31 * t13;
t41 = t37 * t28 + t39 * t65;
t15 = t41 * qJD(3);
t17 = qJ(4) + t41;
t36 = sin(qJ(5));
t38 = cos(qJ(5));
t35 = cos(pkin(9));
t60 = t35 * t38;
t61 = t35 * t36;
t21 = -t35 * pkin(4) - t33 * pkin(7) - pkin(3);
t9 = t21 + t40;
t2 = -t13 * t60 - t36 * t15 + (t17 * t61 - t38 * t9) * qJD(5);
t66 = t38 * t11 - t2 * t35;
t50 = qJD(5) * t38;
t48 = t31 * t50;
t64 = t36 * t11 + t17 * t48;
t29 = t31 * qJD(4);
t53 = qJD(4) * t35;
t54 = qJ(4) * t35;
t5 = -t38 * t53 + (-t21 * t38 + t36 * t54) * qJD(5);
t63 = t38 * t29 - t5 * t35;
t62 = t15 * t35;
t32 = t35 ^ 2;
t58 = t32 * t13 + t11;
t52 = qJD(5) * t31;
t45 = qJ(4) * t52;
t57 = t36 * t29 + t38 * t45;
t56 = t32 * qJD(4) + t29;
t55 = t31 + t32;
t51 = qJD(5) * t36;
t49 = t31 * t51;
t47 = t33 * t51;
t46 = t33 * t50;
t44 = t55 * t13;
t43 = t55 * qJD(4);
t42 = 0.2e1 * qJD(5) * t33 * t35;
t26 = t35 * t50;
t25 = t35 * t51;
t20 = -0.2e1 * t36 * t48;
t19 = t38 * t42;
t18 = t36 * t42;
t16 = 0.2e1 * (t36 ^ 2 - t38 ^ 2) * t52;
t6 = -t36 * t53 + (-t21 * t36 - t38 * t54) * qJD(5);
t3 = -t13 * t61 + t38 * t15 + (-t17 * t60 - t36 * t9) * qJD(5);
t1 = [0, 0, 0, 0, 0, -0.2e1 * t15, 0.2e1 * t14, -0.2e1 * t62, 0.2e1 * t58, 0.2e1 * (-pkin(3) + t40) * t15 + 0.2e1 * t17 * t44, t20, t16, t18, t19, 0, -0.2e1 * t3 * t35 + 0.2e1 * t64, -0.2e1 * t17 * t49 + 0.2e1 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t15, t14, -t62, t56 + t58, -t15 * pkin(3) + qJ(4) * t44 + t17 * t43, t20, t16, t18, t19, 0, (-t3 - t6) * t35 + t57 + t64, (-qJ(4) - t17) * t49 + t63 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t56, 0.2e1 * qJ(4) * t43, t20, t16, t18, t19, 0, -0.2e1 * t6 * t35 + 0.2e1 * t57, -0.2e1 * t36 * t45 + 0.2e1 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, t25, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t46, 0, t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t46, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
