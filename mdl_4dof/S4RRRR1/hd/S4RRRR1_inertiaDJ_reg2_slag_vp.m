% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:14
% EndTime: 2019-12-31 17:22:15
% DurationCPUTime: 0.36s
% Computational Cost: add. (258->54), mult. (802->95), div. (0->0), fcn. (493->6), ass. (0->49)
t35 = cos(qJ(3));
t58 = cos(qJ(2));
t46 = t58 * pkin(1);
t39 = t46 + pkin(2);
t37 = t35 * t39;
t40 = qJD(2) * t46;
t32 = sin(qJ(3));
t33 = sin(qJ(2));
t59 = pkin(1) * t33;
t52 = t32 * t59;
t61 = qJD(2) + qJD(3);
t5 = -qJD(3) * t37 - t35 * t40 + t61 * t52;
t31 = sin(qJ(4));
t29 = t31 ^ 2;
t34 = cos(qJ(4));
t30 = t34 ^ 2;
t55 = t29 + t30;
t43 = t55 * t5;
t57 = t33 * t35;
t63 = t61 * pkin(1) * (t58 * t32 + t57);
t28 = t34 * qJD(4);
t54 = qJD(3) * pkin(2);
t48 = t32 * t54;
t6 = t48 + t63;
t12 = t37 - t52;
t9 = -pkin(3) - t12;
t60 = t9 * t28 + t6 * t31;
t47 = t35 * t54;
t11 = t55 * t47;
t27 = -t35 * pkin(2) - pkin(3);
t56 = t27 * t28 + t31 * t48;
t13 = pkin(1) * t57 + t32 * t39;
t53 = t31 * qJD(4);
t51 = pkin(3) * t53;
t50 = pkin(3) * t28;
t49 = qJD(2) * t59;
t45 = t31 * t28;
t7 = t9 * t53;
t44 = -t6 * t34 + t7;
t42 = t55 * t35;
t41 = -0.2e1 * t48;
t15 = t27 * t53;
t38 = -t34 * t48 + t15;
t26 = t32 * pkin(2) + pkin(7);
t24 = -0.2e1 * t45;
t23 = 0.2e1 * t45;
t14 = 0.2e1 * (-t29 + t30) * qJD(4);
t10 = pkin(7) + t13;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t49, -0.2e1 * t40, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t6, 0.2e1 * t5, 0, -0.2e1 * t12 * t6 - 0.2e1 * t13 * t5, t23, t14, 0, t24, 0, 0, 0.2e1 * t44, 0.2e1 * t60, -0.2e1 * t43, -0.2e1 * t10 * t43 + 0.2e1 * t9 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t40, 0, 0, 0, 0, 0, 0, 0, 0, t41 - t63, t5 - t47, 0, (-t32 * t5 - t35 * t6 + (-t12 * t32 + t13 * t35) * qJD(3)) * pkin(2), t23, t14, 0, t24, 0, 0, t15 + t7 + (-t6 - t48) * t34, t56 + t60, t11 - t43, t6 * t27 - t26 * t43 + (t10 * t42 + t32 * t9) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -0.2e1 * t47, 0, 0, t23, t14, 0, t24, 0, 0, 0.2e1 * t38, 0.2e1 * t56, 0.2e1 * t11, 0.2e1 * (t26 * t42 + t27 * t32) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, t23, t14, 0, t24, 0, 0, t44 - t51, -t50 + t60, -t43, -t6 * pkin(3) - pkin(7) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, 0, 0, t23, t14, 0, t24, 0, 0, t38 - t51, -t50 + t56, t11, (-pkin(3) * t32 + pkin(7) * t42) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t14, 0, t24, 0, 0, -0.2e1 * t51, -0.2e1 * t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t53, 0, -t10 * t28 + t31 * t5, t10 * t53 + t34 * t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t53, 0, -t26 * t28 - t31 * t47, t26 * t53 - t34 * t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t53, 0, -pkin(7) * t28, pkin(7) * t53, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
