% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:45
% EndTime: 2019-12-05 16:00:49
% DurationCPUTime: 0.64s
% Computational Cost: add. (492->86), mult. (1157->147), div. (0->0), fcn. (1014->6), ass. (0->56)
t55 = qJD(4) + qJD(5);
t32 = sin(qJ(5));
t33 = sin(qJ(4));
t35 = cos(qJ(5));
t36 = cos(qJ(4));
t15 = t32 * t36 + t35 * t33;
t65 = t35 * t36;
t43 = t32 * t33 - t65;
t5 = t55 * t15;
t60 = t33 * qJD(4);
t62 = qJD(5) * t32;
t6 = -t32 * t60 - t33 * t62 + t55 * t65;
t73 = ((t15 * t35 + t32 * t43) * qJD(5) + t32 * t6 - t35 * t5) * pkin(4);
t34 = sin(qJ(2));
t28 = t34 * qJD(2);
t30 = t33 ^ 2;
t31 = t36 ^ 2;
t63 = t30 + t31;
t13 = t63 * t28;
t68 = -pkin(2) - pkin(6);
t56 = pkin(7) - t68;
t47 = t56 * t36;
t72 = t55 * t47;
t37 = cos(qJ(2));
t59 = t36 * qJD(4);
t71 = t33 * t28 - t37 * t59;
t69 = 2 * qJD(3);
t67 = t15 * t6;
t66 = t43 * t5;
t58 = t37 * qJD(2);
t64 = qJ(3) * t58 + t34 * qJD(3);
t61 = qJD(5) * t35;
t57 = qJ(3) * qJD(4);
t54 = pkin(4) * t62;
t53 = pkin(4) * t61;
t51 = t36 * t28;
t22 = t34 * t58;
t50 = t33 * t59;
t48 = t68 * qJD(4);
t46 = t66 + t67;
t40 = t56 * t60;
t17 = t56 * t33;
t1 = -t17 * t62 - t32 * t40 + t72 * t35;
t2 = t17 * t61 + t72 * t32 + t35 * t40;
t7 = t32 * t17 - t35 * t47;
t8 = -t35 * t17 - t32 * t47;
t39 = t1 * t15 + t2 * t43 + t7 * t5 - t8 * t6;
t10 = t15 * t37;
t3 = -t36 * t37 * t61 + (t55 * t37 * t33 + t51) * t32 + t71 * t35;
t4 = -t43 * t28 + t5 * t37;
t9 = t43 * t37;
t38 = -t10 * t6 + t3 * t15 - t4 * t43 - t9 * t5;
t29 = qJ(3) * t69;
t24 = t33 * pkin(4) + qJ(3);
t20 = pkin(4) * t59 + qJD(3);
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-t63 + 0.1e1) * t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t10 * t3 + 0.2e1 * t9 * t4 + 0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t58, -pkin(2) * t28 + t64, 0, 0, 0, 0, 0, 0, t33 * t58 + t34 * t59, -t34 * t60 + t36 * t58, -t13, t68 * t13 + t64, 0, 0, 0, 0, 0, 0, t15 * t58 + t34 * t6, -t34 * t5 - t43 * t58, -t38, t10 * t1 + t9 * t2 + t34 * t20 + t24 * t58 + t3 * t8 + t4 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t29, -0.2e1 * t50, 0.2e1 * (t30 - t31) * qJD(4), 0, 0.2e1 * t50, 0, 0, 0.2e1 * qJD(3) * t33 + 0.2e1 * t36 * t57, 0.2e1 * qJD(3) * t36 - 0.2e1 * t33 * t57, 0, t29, 0.2e1 * t66, 0.2e1 * t5 * t15 + 0.2e1 * t43 * t6, 0, 0.2e1 * t67, 0, 0, 0.2e1 * t20 * t15 + 0.2e1 * t24 * t6, -0.2e1 * t20 * t43 - 0.2e1 * t24 * t5, 0.2e1 * t39, -0.2e1 * t8 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t24 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t46, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t60 + t51, -t71, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, (t3 * t32 + t35 * t4 + (-t10 * t35 - t32 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, 0, -t59, 0, -t33 * t48, -t36 * t48, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, -t73, (-t1 * t32 + t2 * t35 + (-t32 * t7 + t35 * t8) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t59, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t54, -0.2e1 * t53, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t53, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
