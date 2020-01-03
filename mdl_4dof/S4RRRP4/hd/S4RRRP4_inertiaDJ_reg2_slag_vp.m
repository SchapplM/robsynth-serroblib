% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:41
% EndTime: 2019-12-31 17:15:43
% DurationCPUTime: 0.36s
% Computational Cost: add. (442->67), mult. (1112->129), div. (0->0), fcn. (906->4), ass. (0->48)
t37 = sin(qJ(3));
t38 = sin(qJ(2));
t56 = t37 * t38;
t59 = qJD(2) + qJD(3);
t60 = t59 * t56;
t58 = -pkin(6) - pkin(5);
t57 = cos(qJ(3));
t39 = cos(qJ(2));
t25 = t37 * t39 + t57 * t38;
t14 = t59 * t25;
t44 = t57 * t39;
t24 = -t44 + t56;
t43 = t57 * qJD(3);
t42 = pkin(2) * t43;
t55 = -t37 * pkin(2) * t14 - t24 * t42;
t29 = t58 * t39;
t48 = t37 * t58;
t16 = -t57 * t29 + t38 * t48;
t54 = qJD(3) * t37;
t53 = t38 * qJD(2);
t52 = t39 * qJD(2);
t51 = -0.2e1 * pkin(1) * qJD(2);
t50 = pkin(2) * t53;
t49 = pkin(2) * t54;
t47 = t57 * pkin(2);
t46 = t25 * t54;
t45 = t38 * t52;
t36 = -t39 * pkin(2) - pkin(1);
t41 = t58 * t57;
t27 = t38 * t41;
t15 = t37 * t29 + t27;
t40 = qJD(2) * t41;
t4 = -qJD(3) * t27 - t29 * t54 - t38 * t40 - t48 * t52;
t5 = t29 * t43 + t39 * t40 - t58 * t60;
t35 = t47 + pkin(3);
t33 = -0.2e1 * t42;
t32 = -0.2e1 * t49;
t17 = t24 * pkin(3) + t36;
t13 = -qJD(2) * t44 - t39 * t43 + t60;
t10 = t14 * pkin(3) + t50;
t9 = -t24 * qJ(4) + t16;
t8 = -t25 * qJ(4) + t15;
t7 = -0.2e1 * t25 * t13;
t6 = 0.2e1 * t24 * t14;
t3 = 0.2e1 * t13 * t24 - 0.2e1 * t25 * t14;
t2 = t13 * qJ(4) - t25 * qJD(4) + t5;
t1 = t14 * qJ(4) + t24 * qJD(4) + t4;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t45, 0.2e1 * (-t38 ^ 2 + t39 ^ 2) * qJD(2), 0, -0.2e1 * t45, 0, 0, t38 * t51, t39 * t51, 0, 0, t7, t3, 0, t6, 0, 0, 0.2e1 * t36 * t14 + 0.2e1 * t24 * t50, -0.2e1 * t36 * t13 + 0.2e1 * t25 * t50, 0.2e1 * t15 * t13 - 0.2e1 * t16 * t14 + 0.2e1 * t4 * t24 - 0.2e1 * t5 * t25, 0.2e1 * t15 * t5 - 0.2e1 * t16 * t4 + 0.2e1 * t36 * t50, t7, t3, 0, t6, 0, 0, 0.2e1 * t10 * t24 + 0.2e1 * t17 * t14, 0.2e1 * t10 * t25 - 0.2e1 * t17 * t13, 0.2e1 * t1 * t24 + 0.2e1 * t8 * t13 - 0.2e1 * t9 * t14 - 0.2e1 * t2 * t25, -0.2e1 * t9 * t1 + 0.2e1 * t17 * t10 + 0.2e1 * t8 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, -t53, 0, -pkin(5) * t52, pkin(5) * t53, 0, 0, 0, 0, -t13, 0, -t14, 0, t5, t4, (t57 * t13 + t46) * pkin(2) + t55, (t57 * t5 - t37 * t4 + (-t15 * t37 + t57 * t16) * qJD(3)) * pkin(2), 0, 0, -t13, 0, -t14, 0, t2, t1, pkin(2) * t46 + t35 * t13 + t55, t2 * t35 + (-t1 * t37 + (-t37 * t8 + t57 * t9) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t33, 0, 0, 0, 0, 0, 0, 0, 0, t32, t33, 0, 0.2e1 * (t47 - t35) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, -t14, 0, t5, t4, 0, 0, 0, 0, -t13, 0, -t14, 0, t2, t1, t13 * pkin(3), t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t42, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t42, 0, -pkin(3) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
