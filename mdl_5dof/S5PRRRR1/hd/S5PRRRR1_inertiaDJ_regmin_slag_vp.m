% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_inertiaDJ_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:49
% EndTime: 2019-07-18 13:28:50
% DurationCPUTime: 0.45s
% Computational Cost: add. (347->72), mult. (1115->163), div. (0->0), fcn. (1073->8), ass. (0->65)
t32 = sin(qJ(2));
t35 = cos(qJ(3));
t36 = cos(qJ(2));
t53 = t36 * qJD(2);
t49 = t35 * t53;
t31 = sin(qJ(3));
t56 = t31 * qJD(3);
t70 = -t32 * t56 + t49;
t33 = cos(qJ(5));
t28 = t33 ^ 2;
t29 = sin(qJ(5));
t62 = t29 ^ 2 - t28;
t45 = t62 * qJD(5);
t69 = qJD(3) + qJD(4);
t68 = 2 * pkin(2);
t30 = sin(qJ(4));
t34 = cos(qJ(4));
t20 = t30 * t35 + t31 * t34;
t12 = t69 * t20;
t67 = t12 * t35;
t19 = t30 * t31 - t34 * t35;
t11 = t69 * t19;
t66 = t20 * t11;
t65 = t20 * t29;
t64 = t20 * t33;
t63 = t29 * t33;
t61 = qJD(4) * t30;
t60 = qJD(4) * t34;
t59 = qJD(5) * t29;
t26 = qJD(5) * t33;
t58 = qJD(5) * t34;
t57 = qJD(5) * t35;
t55 = t32 * qJD(2);
t54 = t35 * qJD(3);
t52 = pkin(2) * t61;
t51 = pkin(2) * t60;
t50 = t31 * t53;
t48 = t29 * t26;
t46 = -0.4e1 * t20 * t63;
t14 = t19 * t32;
t44 = -t33 * t14 - t29 * t36;
t43 = -t29 * t14 + t33 * t36;
t42 = t19 * t30 + t20 * t34;
t41 = -t11 * t29 + t20 * t26;
t40 = -t11 * t33 - t20 * t59;
t39 = t29 * t57 + t33 * t56;
t38 = -t29 * t56 + t33 * t57;
t37 = t11 * t34 - t12 * t30 + (-t19 * t34 + t20 * t30) * qJD(4);
t24 = 0.2e1 * t48;
t18 = -0.2e1 * t45;
t17 = t20 ^ 2;
t16 = (-t29 * t58 - t33 * t61) * pkin(2);
t15 = (t29 * t61 - t33 * t58) * pkin(2);
t13 = t20 * t32;
t10 = t12 * t29 + t19 * t26;
t9 = t12 * t33 - t19 * t59;
t8 = -t31 * t32 * t61 + (t69 * t35 * t32 + t50) * t34 + t70 * t30;
t7 = t12 * t32 + t30 * t50 - t34 * t49;
t6 = -t11 * t63 - t20 * t45;
t5 = t13 * t59 - t33 * t8;
t4 = t13 * t26 + t29 * t8;
t3 = qJD(5) * t46 + t62 * t11;
t2 = -t44 * qJD(5) + t29 * t7 + t33 * t55;
t1 = t43 * qJD(5) - t29 * t55 + t33 * t7;
t21 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t55, -t53, 0, 0, 0, 0, 0, -t35 * t55 - t36 * t56, t31 * t55 - t36 * t54, 0, 0, 0, 0, 0, -t12 * t36 + t19 * t55, t11 * t36 + t20 * t55, 0, 0, 0, 0, 0, -t43 * t12 + t41 * t13 + t2 * t19 + t8 * t65, t1 * t19 - t44 * t12 + t40 * t13 + t8 * t64; 0, 0, 0, 0, 0.2e1 * t31 * t54, 0.2e1 * (-t31 ^ 2 + t35 ^ 2) * qJD(3), 0, 0, 0, 0, 0, -0.2e1 * t66, 0.2e1 * t11 * t19 - 0.2e1 * t12 * t20, 0, 0, 0, (t19 * t56 - t67) * t68, (t11 * t35 + t20 * t56) * t68, -0.2e1 * t17 * t48 - 0.2e1 * t28 * t66, -t11 * t46 + 0.2e1 * t17 * t45, 0.2e1 * t12 * t64 + 0.2e1 * t40 * t19, -0.2e1 * t12 * t65 - 0.2e1 * t41 * t19, 0.2e1 * t19 * t12, (t39 * t19 - t33 * t67) * t68, (t38 * t19 + t29 * t67) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32 * t54 - t50, -t70, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t5, t4; 0, 0, 0, 0, 0, 0, t54, -t56, 0, 0, 0, 0, 0, -t11, -t12, 0, 0, 0, t6, t3, t10, t9, 0, (-t42 * t26 + t37 * t29) * pkin(2), (t37 * t33 + t42 * t59) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t52, -0.2e1 * t51, t24, t18, 0, 0, 0, 0.2e1 * t16, 0.2e1 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, 0, 0, t6, t3, t10, t9, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t51, t24, t18, 0, 0, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t18, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t41, t12, t39 * pkin(2), t38 * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t59, 0, (-t30 * t26 - t29 * t60) * pkin(2), (t30 * t59 - t33 * t60) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t59, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t21;
