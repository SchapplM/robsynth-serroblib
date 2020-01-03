% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:28
% EndTime: 2019-12-31 16:56:30
% DurationCPUTime: 0.40s
% Computational Cost: add. (280->69), mult. (649->150), div. (0->0), fcn. (438->4), ass. (0->61)
t21 = cos(qJ(4));
t20 = sin(qJ(3));
t22 = cos(qJ(3));
t30 = pkin(3) * t22 + pkin(6) * t20;
t64 = t21 * t30;
t19 = sin(qJ(4));
t15 = t19 ^ 2;
t17 = t21 ^ 2;
t58 = t15 - t17;
t34 = qJD(4) * t58;
t16 = t20 ^ 2;
t18 = t22 ^ 2;
t33 = (t16 - t18) * qJD(3);
t63 = 2 * qJD(2);
t62 = t20 * pkin(3);
t61 = t22 * pkin(6);
t23 = -pkin(1) - pkin(5);
t60 = t20 * t23;
t59 = t22 * t23;
t57 = t15 + t17;
t55 = t16 + t18;
t54 = qJD(4) * t19;
t53 = qJD(4) * t21;
t52 = qJD(4) * t22;
t51 = qJD(4) * t23;
t50 = t20 * qJD(3);
t49 = t22 * qJD(3);
t48 = qJ(2) * qJD(3);
t47 = -0.2e1 * pkin(3) * qJD(4);
t46 = t19 * t60;
t45 = t19 * t59;
t44 = t21 * t60;
t43 = t19 * t52;
t42 = t19 * t51;
t41 = t21 * t52;
t40 = t19 * t53;
t39 = t23 * t50;
t38 = t21 * t50;
t37 = t20 * t49;
t36 = t23 * t49;
t35 = t57 * t22;
t13 = 0.2e1 * t37;
t32 = t19 * t38;
t31 = t18 * t40;
t29 = -t61 + t62;
t26 = qJ(2) + t29;
t25 = t21 * t26;
t4 = t25 - t46;
t5 = t19 * t26 + t44;
t28 = t19 * t5 + t21 * t4;
t27 = t19 * t4 - t21 * t5;
t1 = t20 * t42 - qJD(4) * t25 - t21 * t36 - t19 * (t30 * qJD(3) + qJD(2));
t2 = t21 * qJD(2) - qJD(4) * t5 + (-t45 + t64) * qJD(3);
t24 = -t28 * qJD(4) - t1 * t21 - t2 * t19;
t14 = qJ(2) * t63;
t10 = t19 * t50 - t41;
t9 = t19 * t49 + t20 * t53;
t8 = -t38 - t43;
t7 = t20 * t54 - t21 * t49;
t3 = t22 * t34 + t32;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t14, -0.2e1 * t37, 0.2e1 * t33, 0, t13, 0, 0, 0.2e1 * qJD(2) * t20 + 0.2e1 * t22 * t48, 0.2e1 * qJD(2) * t22 - 0.2e1 * t20 * t48, 0, t14, -0.2e1 * t17 * t37 - 0.2e1 * t31, 0.2e1 * t18 * t34 + 0.4e1 * t22 * t32, -0.2e1 * t20 * t43 - 0.2e1 * t21 * t33, -0.2e1 * t15 * t37 + 0.2e1 * t31, 0.2e1 * t19 * t33 - 0.2e1 * t20 * t41, t13, -0.2e1 * t18 * t21 * t51 + 0.2e1 * t2 * t20 + 0.2e1 * (t4 + 0.2e1 * t46) * t49, 0.2e1 * t18 * t42 + 0.2e1 * t1 * t20 + 0.2e1 * (-t5 + 0.2e1 * t44) * t49, 0.2e1 * t28 * t50 + 0.2e1 * (t27 * qJD(4) + t1 * t19 - t2 * t21) * t22, -0.2e1 * t23 ^ 2 * t37 - 0.2e1 * t5 * t1 + 0.2e1 * t4 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 * t53, t55 * t54, 0, -t27 * t49 + (t24 - 0.2e1 * t36) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t57) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, 0, -t49, 0, -t39, -t36, 0, 0, -t3, -0.4e1 * t22 * t40 + t58 * t50, t9, t3, -t7, 0, (-t45 - t64) * qJD(4) + (t29 * t19 - t44) * qJD(3), (t30 * t19 - t21 * t59) * qJD(4) + (-t21 * t61 + (t21 * pkin(3) + t19 * t23) * t20) * qJD(3), t24, -pkin(3) * t39 + t24 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t49, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10, qJD(3) * t35, (pkin(6) * t35 - t62) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t40, -0.2e1 * t34, 0, -0.2e1 * t40, 0, 0, t19 * t47, t21 * t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, t10, t49, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, -t54, 0, -pkin(6) * t53, pkin(6) * t54, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
