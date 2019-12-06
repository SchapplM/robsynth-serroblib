% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:59
% EndTime: 2019-12-05 17:13:01
% DurationCPUTime: 0.45s
% Computational Cost: add. (617->79), mult. (1613->153), div. (0->0), fcn. (1492->8), ass. (0->64)
t40 = sin(qJ(4));
t41 = sin(qJ(3));
t44 = cos(qJ(4));
t45 = cos(qJ(3));
t27 = t40 * t41 - t44 * t45;
t42 = sin(qJ(2));
t25 = t27 * t42;
t70 = qJD(3) + qJD(4);
t69 = pkin(6) + pkin(7);
t52 = qJD(3) * t69;
t29 = t41 * t52;
t30 = t45 * t52;
t32 = t69 * t41;
t33 = t69 * t45;
t49 = -t44 * t32 - t40 * t33;
t11 = -t49 * qJD(4) + t44 * t29 + t40 * t30;
t39 = sin(qJ(5));
t67 = pkin(3) * qJD(4);
t58 = t40 * t67;
t66 = qJD(5) * t39;
t68 = t40 * pkin(3) * t66 + t39 * t58;
t43 = cos(qJ(5));
t65 = qJD(5) * t43;
t64 = t41 * qJD(3);
t63 = t42 * qJD(2);
t62 = t45 * qJD(3);
t46 = cos(qJ(2));
t61 = t46 * qJD(2);
t60 = -0.2e1 * pkin(2) * qJD(3);
t59 = pkin(3) * t64;
t57 = t44 * t67;
t56 = pkin(4) * t66;
t55 = pkin(4) * t65;
t54 = t41 * t61;
t53 = t45 * t61;
t38 = -t45 * pkin(3) - pkin(2);
t37 = t44 * pkin(3) + pkin(4);
t51 = (-pkin(4) - t37) * qJD(5);
t28 = t40 * t45 + t44 * t41;
t19 = t43 * t27 + t39 * t28;
t20 = -t39 * t27 + t43 * t28;
t48 = t40 * t32 - t44 * t33;
t12 = t48 * qJD(4) + t40 * t29 - t44 * t30;
t47 = (-t40 * t65 + (-t39 * t44 - t40 * t43) * qJD(4)) * pkin(3);
t22 = t70 * t28;
t24 = t28 * t42;
t23 = t27 * pkin(4) + t38;
t21 = t70 * t27;
t17 = -t37 * t66 + t47;
t16 = (-qJD(5) * t37 - t57) * t43 + t68;
t15 = t22 * pkin(4) + t59;
t14 = -t27 * pkin(8) - t48;
t13 = -t28 * pkin(8) + t49;
t10 = t70 * t25 - t28 * t61;
t9 = t22 * t42 + t40 * t54 - t44 * t53;
t8 = t21 * pkin(8) + t12;
t7 = -t22 * pkin(8) - t11;
t6 = t20 * qJD(5) - t39 * t21 + t43 * t22;
t5 = -t19 * qJD(5) - t43 * t21 - t39 * t22;
t4 = t43 * t10 + t39 * t9 + (t24 * t39 + t25 * t43) * qJD(5);
t3 = -t39 * t10 + t43 * t9 + (t24 * t43 - t25 * t39) * qJD(5);
t2 = -t39 * t7 + t43 * t8 + (-t13 * t39 - t14 * t43) * qJD(5);
t1 = -t39 * t8 - t43 * t7 + (-t13 * t43 + t14 * t39) * qJD(5);
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t63, -t61, 0, 0, 0, 0, 0, -t45 * t63 - t46 * t64, t41 * t63 - t46 * t62, 0, 0, 0, 0, 0, -t46 * t22 + t27 * t63, t46 * t21 + t28 * t63, 0, 0, 0, 0, 0, t19 * t63 - t46 * t6, t20 * t63 - t46 * t5; 0, 0, 0, 0, 0.2e1 * t41 * t62, 0.2e1 * (-t41 ^ 2 + t45 ^ 2) * qJD(3), 0, 0, 0, t41 * t60, t45 * t60, -0.2e1 * t28 * t21, 0.2e1 * t21 * t27 - 0.2e1 * t28 * t22, 0, 0, 0, 0.2e1 * t38 * t22 + 0.2e1 * t27 * t59, -0.2e1 * t38 * t21 + 0.2e1 * t28 * t59, 0.2e1 * t20 * t5, -0.2e1 * t5 * t19 - 0.2e1 * t20 * t6, 0, 0, 0, 0.2e1 * t15 * t19 + 0.2e1 * t23 * t6, 0.2e1 * t15 * t20 + 0.2e1 * t23 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42 * t62 - t54, t42 * t64 - t53, 0, 0, 0, 0, 0, t10, t9, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, t62, -t64, 0, -pkin(6) * t62, pkin(6) * t64, 0, 0, -t21, -t22, 0, t12, t11, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t58, -0.2e1 * t57, 0, 0, 0, 0, 0, 0.2e1 * t17, 0.2e1 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22, 0, t12, t11, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t57, 0, 0, 0, 0, 0, t39 * t51 + t47, (t51 - t57) * t43 + t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t56, -0.2e1 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t18;
