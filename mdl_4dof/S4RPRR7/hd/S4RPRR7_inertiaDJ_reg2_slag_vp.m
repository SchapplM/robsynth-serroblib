% Calculate inertial parameters regressor of joint inertia matrix time derivative for
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
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:05
% EndTime: 2019-12-31 16:54:07
% DurationCPUTime: 0.49s
% Computational Cost: add. (645->75), mult. (1529->158), div. (0->0), fcn. (1400->6), ass. (0->58)
t28 = sin(pkin(7));
t29 = cos(pkin(7));
t31 = sin(qJ(3));
t65 = cos(qJ(3));
t18 = t65 * t28 + t31 * t29;
t72 = -0.4e1 * t18;
t37 = -t31 * t28 + t65 * t29;
t59 = pkin(5) + qJ(2);
t19 = t59 * t28;
t20 = t59 * t29;
t38 = -t65 * t19 - t31 * t20;
t33 = t37 * qJD(2) + t38 * qJD(3);
t23 = -t29 * pkin(2) - pkin(1);
t36 = -pkin(3) * t37 - t18 * pkin(6) + t23;
t71 = -qJD(4) * t36 - t33;
t30 = sin(qJ(4));
t26 = t30 ^ 2;
t32 = cos(qJ(4));
t27 = t32 ^ 2;
t58 = t26 - t27;
t50 = qJD(4) * t58;
t12 = -t31 * t19 + t65 * t20;
t13 = t37 * qJD(3);
t14 = t18 * qJD(3);
t47 = t14 * pkin(3) - t13 * pkin(6);
t57 = qJD(4) * t30;
t1 = t12 * t57 - t30 * t47 + t32 * t71;
t56 = qJD(4) * t32;
t2 = -t12 * t56 + t71 * t30 + t32 * t47;
t3 = -t30 * t12 + t32 * t36;
t4 = t32 * t12 + t30 * t36;
t70 = t1 * t30 - t2 * t32 + (t3 * t30 - t32 * t4) * qJD(4);
t69 = 0.2e1 * t14;
t7 = t18 * qJD(2) + t12 * qJD(3);
t68 = t38 * t7;
t67 = t7 * t30;
t66 = t7 * t32;
t64 = t18 * t13;
t63 = t30 * t14;
t61 = t32 * t13;
t60 = t32 * t14;
t55 = t37 * t69;
t54 = -0.2e1 * pkin(3) * qJD(4);
t53 = t30 * t61;
t52 = t30 * t56;
t15 = t18 ^ 2;
t49 = t15 * t52;
t48 = 0.2e1 * (t28 ^ 2 + t29 ^ 2) * qJD(2);
t46 = -pkin(3) * t13 - pkin(6) * t14;
t45 = pkin(3) * t18 - pkin(6) * t37;
t43 = t3 * t32 + t30 * t4;
t41 = t30 * t13 + t18 * t56;
t40 = -t18 * t57 + t61;
t39 = -t37 * t56 + t63;
t34 = -t43 * qJD(4) - t1 * t32 - t2 * t30;
t10 = t37 * t57 + t60;
t5 = t18 * t50 - t53;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, qJ(2) * t48, 0.2e1 * t64, 0.2e1 * t13 * t37 - 0.2e1 * t18 * t14, 0, -t55, 0, 0, t23 * t69, 0.2e1 * t23 * t13, -0.2e1 * t12 * t14 - 0.2e1 * t13 * t38 + 0.2e1 * t7 * t18 + 0.2e1 * t33 * t37, 0.2e1 * t12 * t33 - 0.2e1 * t68, 0.2e1 * t27 * t64 - 0.2e1 * t49, 0.2e1 * t15 * t50 + t53 * t72, 0.2e1 * t18 * t60 - 0.2e1 * t37 * t40, 0.2e1 * t26 * t64 + 0.2e1 * t49, -0.2e1 * t18 * t63 + 0.2e1 * t37 * t41, -t55, 0.2e1 * t3 * t14 + 0.2e1 * t18 * t67 - 0.2e1 * t2 * t37 - 0.2e1 * t38 * t41, -0.2e1 * t1 * t37 - 0.2e1 * t4 * t14 + 0.2e1 * t18 * t66 - 0.2e1 * t38 * t40, -0.2e1 * t43 * t13 + 0.2e1 * t18 * t70, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 - 0.2e1 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t13, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t39, -(t26 + t27) * t13, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, -t14, 0, -t7, -t33, 0, 0, -t5, -t58 * t13 + t52 * t72, t39, t5, t10, 0, -t66 + t46 * t30 + (-t30 * t38 - t45 * t32) * qJD(4), t67 + t46 * t32 + (t45 * t30 - t32 * t38) * qJD(4), t34, -t7 * pkin(3) + t34 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t52, -0.2e1 * t50, 0, -0.2e1 * t52, 0, 0, t30 * t54, t32 * t54, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t41, t14, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t56, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, -t57, 0, -pkin(6) * t56, pkin(6) * t57, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
