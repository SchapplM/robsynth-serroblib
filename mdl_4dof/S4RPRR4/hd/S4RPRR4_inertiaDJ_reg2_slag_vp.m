% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RPRR4
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
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:35
% EndTime: 2019-12-31 16:50:37
% DurationCPUTime: 0.39s
% Computational Cost: add. (294->64), mult. (720->139), div. (0->0), fcn. (518->6), ass. (0->58)
t23 = sin(qJ(3));
t61 = -0.4e1 * t23;
t22 = sin(qJ(4));
t18 = t22 ^ 2;
t24 = cos(qJ(4));
t20 = t24 ^ 2;
t37 = qJD(4) * (t18 - t20);
t19 = t23 ^ 2;
t25 = cos(qJ(3));
t55 = -t25 ^ 2 + t19;
t36 = t55 * qJD(3);
t60 = pkin(3) * t23;
t59 = pkin(6) * t25;
t16 = sin(pkin(7)) * pkin(1) + pkin(5);
t58 = t22 * t16;
t57 = t24 * t25;
t54 = qJD(4) * t22;
t53 = qJD(4) * t24;
t52 = qJD(4) * t25;
t51 = t23 * qJD(3);
t50 = t25 * qJD(3);
t49 = -0.2e1 * pkin(3) * qJD(4);
t48 = t25 * t58;
t47 = t16 * t57;
t17 = -cos(pkin(7)) * pkin(1) - pkin(2);
t46 = 0.2e1 * qJD(3) * t17;
t45 = qJD(4) * t16 * t19;
t44 = t22 * t52;
t43 = t24 * t52;
t42 = t18 * t50;
t41 = t22 * t53;
t40 = t23 * t50;
t39 = t24 * t50;
t38 = t16 * t50;
t35 = t22 * t39;
t34 = t19 * t41;
t33 = -t25 * pkin(3) - t23 * pkin(6);
t32 = -t59 + t60;
t28 = -t17 - t33;
t27 = t24 * t28;
t3 = -t27 - t48;
t4 = -t22 * t28 + t47;
t31 = -t22 * t4 - t24 * t3;
t30 = t22 * t3 - t24 * t4;
t29 = t32 * t22;
t8 = t24 * t51 + t44;
t1 = -qJD(3) * t29 + qJD(4) * t27 + t8 * t16;
t2 = -qJD(4) * t4 + (t23 * t58 + t24 * t32) * qJD(3);
t26 = qJD(4) * t31 - t1 * t24 - t2 * t22;
t15 = t20 * t50;
t14 = -0.2e1 * t40;
t13 = t20 * t40;
t12 = t18 * t40;
t10 = t22 * t51 - t43;
t9 = -t22 * t50 - t23 * t53;
t7 = t23 * t54 - t39;
t5 = t23 * t37 - t35;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t40, -0.2e1 * t36, 0, t14, 0, 0, t23 * t46, t25 * t46, 0, 0, 0.2e1 * t13 - 0.2e1 * t34, 0.2e1 * t19 * t37 + t35 * t61, 0.2e1 * t23 * t44 + 0.2e1 * t24 * t36, 0.2e1 * t12 + 0.2e1 * t34, -0.2e1 * t22 * t36 + 0.2e1 * t23 * t43, t14, 0.2e1 * t24 * t45 - 0.2e1 * t2 * t25 + 0.2e1 * (t3 + 0.2e1 * t48) * t51, -0.2e1 * t22 * t45 - 0.2e1 * t1 * t25 + 0.2e1 * (-t4 + 0.2e1 * t47) * t51, 0.2e1 * t31 * t50 + 0.2e1 * (qJD(4) * t30 + t1 * t22 - t2 * t24) * t23, 0.2e1 * t16 ^ 2 * t40 - 0.2e1 * t1 * t4 + 0.2e1 * t2 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t23 + (t16 * t55 - t25 * t30) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t12 + 0.2e1 * t13 - 0.2e1 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t51, 0, -t38, t16 * t51, 0, 0, -t5, t41 * t61 + t15 - t42, t10, t5, t8, 0, (pkin(6) * t57 + (-pkin(3) * t24 + t58) * t23) * qJD(4) + (t22 * t33 - t47) * qJD(3), (t16 * t23 * t24 + t29) * qJD(4) + (t24 * t33 + t48) * qJD(3), t26, -pkin(3) * t38 + pkin(6) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t10, t15 + t42, (-t60 + (t18 + t20) * t59) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t41, -0.2e1 * t37, 0, -0.2e1 * t41, 0, 0, t22 * t49, t24 * t49, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, t9, t51, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, -t54, 0, -pkin(6) * t53, pkin(6) * t54, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
