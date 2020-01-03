% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RPRR6
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:41
% EndTime: 2019-12-31 16:52:42
% DurationCPUTime: 0.32s
% Computational Cost: add. (577->60), mult. (1336->125), div. (0->0), fcn. (1277->6), ass. (0->44)
t38 = sin(pkin(7));
t39 = cos(pkin(7));
t41 = sin(qJ(3));
t42 = cos(qJ(3));
t30 = t42 * t38 + t41 * t39;
t25 = t30 * qJD(3);
t61 = 0.2e1 * t25;
t55 = pkin(5) + qJ(2);
t31 = t55 * t38;
t32 = t55 * t39;
t17 = -t41 * t31 + t42 * t32;
t53 = qJD(3) * t42;
t56 = t42 * t39;
t10 = (qJD(2) * t38 + qJD(3) * t32) * t41 - qJD(2) * t56 + t31 * t53;
t57 = t41 * t38;
t24 = qJD(3) * t57 - t39 * t53;
t60 = -0.2e1 * t24;
t59 = t25 * pkin(3);
t58 = cos(qJ(4));
t40 = sin(qJ(4));
t52 = qJD(4) * t40;
t51 = pkin(3) * t52;
t35 = -t39 * pkin(2) - pkin(1);
t16 = -t42 * t31 - t41 * t32;
t50 = qJD(4) * t58;
t49 = pkin(3) * t50;
t48 = 0.2e1 * (t38 ^ 2 + t39 ^ 2) * qJD(2);
t47 = -t56 + t57;
t12 = -t30 * pkin(6) + t16;
t13 = -t47 * pkin(6) + t17;
t4 = t40 * t12 + t58 * t13;
t45 = t58 * t47;
t15 = t58 * t30 - t40 * t47;
t44 = t25 * pkin(6) + t10;
t11 = -t30 * qJD(2) - t17 * qJD(3);
t43 = t24 * pkin(6) + t11;
t18 = t47 * pkin(3) + t35;
t14 = t40 * t30 + t45;
t6 = qJD(4) * t15 - t40 * t24 + t58 * t25;
t5 = qJD(4) * t45 + t58 * t24 + t40 * t25 + t30 * t52;
t3 = t58 * t12 - t40 * t13;
t2 = -t4 * qJD(4) + t40 * t44 + t58 * t43;
t1 = -t12 * t50 + t13 * t52 - t40 * t43 + t58 * t44;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, qJ(2) * t48, t30 * t60, 0.2e1 * t24 * t47 - 0.2e1 * t30 * t25, 0, t47 * t61, 0, 0, t35 * t61, t35 * t60, 0.2e1 * t10 * t47 - 0.2e1 * t11 * t30 + 0.2e1 * t16 * t24 - 0.2e1 * t17 * t25, -0.2e1 * t17 * t10 + 0.2e1 * t16 * t11, -0.2e1 * t15 * t5, 0.2e1 * t5 * t14 - 0.2e1 * t15 * t6, 0, 0.2e1 * t14 * t6, 0, 0, 0.2e1 * t14 * t59 + 0.2e1 * t18 * t6, 0.2e1 * t15 * t59 - 0.2e1 * t18 * t5, 0.2e1 * t1 * t14 - 0.2e1 * t2 * t15 + 0.2e1 * t3 * t5 - 0.2e1 * t4 * t6, -0.2e1 * t4 * t1 + 0.2e1 * t18 * t59 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, -t25, 0, t11, t10, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, (t58 * t5 - t40 * t6 + (-t58 * t14 + t15 * t40) * qJD(4)) * pkin(3), (t58 * t2 - t1 * t40 + (-t3 * t40 + t58 * t4) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t51, -0.2e1 * t49, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t49, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
