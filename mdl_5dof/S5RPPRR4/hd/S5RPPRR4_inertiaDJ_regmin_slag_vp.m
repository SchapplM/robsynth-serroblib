% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:55
% EndTime: 2019-12-05 17:44:57
% DurationCPUTime: 0.44s
% Computational Cost: add. (634->80), mult. (1641->165), div. (0->0), fcn. (1596->8), ass. (0->66)
t48 = sin(pkin(9));
t50 = cos(pkin(9));
t55 = cos(qJ(4));
t67 = qJD(4) * t55;
t53 = sin(qJ(4));
t68 = qJD(4) * t53;
t80 = -t48 * t68 + t50 * t67;
t51 = cos(pkin(8));
t79 = -0.2e1 * t51;
t78 = 0.2e1 * t51;
t77 = t51 * pkin(4);
t36 = t55 * t48 + t53 * t50;
t49 = sin(pkin(8));
t27 = t36 * t49;
t38 = -t51 * pkin(2) - t49 * qJ(3) - pkin(1);
t34 = t50 * t38;
t17 = -t50 * t49 * pkin(6) + t34 + (-qJ(2) * t48 - pkin(3)) * t51;
t71 = qJ(2) * t51;
t73 = t48 * t38 + t50 * t71;
t74 = t48 * t49;
t19 = -pkin(6) * t74 + t73;
t57 = -t53 * t17 - t55 * t19;
t11 = -t27 * pkin(7) - t57;
t52 = sin(qJ(5));
t76 = t11 * t52;
t54 = cos(qJ(5));
t75 = t11 * t54;
t37 = pkin(3) * t74 + t49 * qJ(2);
t72 = pkin(4) * qJD(5);
t70 = qJD(2) * t51;
t69 = qJD(3) * t49;
t46 = t49 ^ 2;
t66 = t46 * qJD(2);
t44 = t49 * qJD(2);
t65 = qJ(2) * qJD(2);
t64 = t52 * t72;
t63 = t54 * t72;
t23 = t80 * t49;
t29 = -t48 * t70 - t50 * t69;
t30 = -t48 * t69 + t50 * t70;
t6 = -t17 * t67 + t19 * t68 - t53 * t29 - t55 * t30;
t4 = -t23 * pkin(7) - t6;
t32 = t36 * qJD(4);
t22 = t49 * t32;
t7 = t57 * qJD(4) + t55 * t29 - t53 * t30;
t5 = t22 * pkin(7) + t7;
t60 = -t52 * t4 + t54 * t5;
t35 = -t53 * t48 + t55 * t50;
t28 = t35 * t49;
t10 = -t28 * pkin(7) + t55 * t17 - t53 * t19 - t77;
t59 = -t10 + t77;
t58 = -t54 * t4 - t52 * t5;
t14 = t54 * t27 + t52 * t28;
t15 = -t52 * t27 + t54 * t28;
t56 = t29 * t50 + t30 * t48;
t47 = t51 ^ 2;
t43 = t46 * t65;
t20 = t23 * pkin(4) + t44;
t18 = t27 * pkin(4) + t37;
t13 = -t52 * t80 - t54 * t32 + (-t35 * t52 - t36 * t54) * qJD(5);
t12 = -t54 * t80 + t52 * t32 + (-t35 * t54 + t36 * t52) * qJD(5);
t9 = t15 * qJD(5) - t52 * t22 + t54 * t23;
t8 = -t14 * qJD(5) - t54 * t22 - t52 * t23;
t2 = (-t10 * t52 - t75) * qJD(5) + t60;
t1 = (-t10 * t54 + t76) * qJD(5) + t58;
t3 = [0, 0, 0, 0, 0, 0.2e1 * (t46 + t47) * qJD(2), 0.2e1 * t47 * t65 + 0.2e1 * t43, -0.2e1 * t29 * t51 + 0.2e1 * t48 * t66, 0.2e1 * t30 * t51 + 0.2e1 * t50 * t66, -0.2e1 * t56 * t49, 0.2e1 * t73 * t30 + 0.2e1 * (-t48 * t71 + t34) * t29 + 0.2e1 * t43, -0.2e1 * t28 * t22, 0.2e1 * t22 * t27 - 0.2e1 * t28 * t23, -t22 * t79, t23 * t78, 0, 0.2e1 * t37 * t23 + 0.2e1 * t27 * t44 - 0.2e1 * t7 * t51, -0.2e1 * t37 * t22 + 0.2e1 * t28 * t44 - 0.2e1 * t6 * t51, 0.2e1 * t15 * t8, -0.2e1 * t8 * t14 - 0.2e1 * t15 * t9, t8 * t79, t9 * t78, 0, 0.2e1 * t20 * t14 + 0.2e1 * t18 * t9 - 0.2e1 * t2 * t51, -0.2e1 * t1 * t51 + 0.2e1 * t20 * t15 + 0.2e1 * t18 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, t32 * t51, t80 * t51, 0, 0, 0, 0, 0, -t13 * t51, -t12 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, t23, -t22, 0, 0, 0, 0, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, 0, t7, t6, 0, 0, t8, -t9, 0, (t59 * t52 - t75) * qJD(5) + t60, (t59 * t54 + t76) * qJD(5) + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t80, 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t64, -0.2e1 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
