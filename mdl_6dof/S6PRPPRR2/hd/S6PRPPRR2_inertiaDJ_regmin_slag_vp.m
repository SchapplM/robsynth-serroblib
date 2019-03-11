% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x22]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPPRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:01
% EndTime: 2019-03-08 19:20:02
% DurationCPUTime: 0.46s
% Computational Cost: add. (266->87), mult. (809->174), div. (0->0), fcn. (744->10), ass. (0->72)
t35 = cos(qJ(6));
t26 = t35 ^ 2;
t32 = sin(qJ(6));
t71 = t32 ^ 2 - t26;
t47 = t71 * qJD(6);
t75 = 2 * qJD(4);
t30 = cos(pkin(11));
t48 = -t30 * pkin(2) - pkin(3);
t22 = -pkin(8) + t48;
t74 = t22 * t32;
t33 = sin(qJ(5));
t73 = t22 * t33;
t36 = cos(qJ(5));
t72 = t35 * t36;
t25 = t33 ^ 2;
t27 = t36 ^ 2;
t70 = t25 - t27;
t69 = t25 + t27;
t29 = sin(pkin(6));
t68 = qJD(2) * t29;
t67 = qJD(5) * t33;
t66 = qJD(5) * t35;
t65 = qJD(5) * t36;
t64 = qJD(6) * t32;
t63 = qJD(6) * t35;
t62 = qJD(6) * t36;
t28 = sin(pkin(11));
t34 = sin(qJ(2));
t37 = cos(qJ(2));
t14 = (t28 * t37 + t30 * t34) * t29;
t61 = t14 * qJD(5);
t60 = -0.2e1 * pkin(5) * qJD(6);
t59 = t35 * t73;
t58 = t32 * t62;
t57 = t35 * t62;
t56 = t34 * t68;
t55 = t37 * t68;
t54 = t32 * t65;
t53 = t32 * t63;
t52 = t33 * t66;
t51 = t35 * t65;
t50 = t33 * t65;
t49 = t22 * t65;
t23 = t28 * pkin(2) + qJ(4);
t46 = t70 * qJD(5);
t45 = t35 * t50;
t44 = pkin(5) * t36 + pkin(9) * t33;
t43 = t33 * pkin(5) - t36 * pkin(9);
t13 = (t28 * t34 - t30 * t37) * t29;
t31 = cos(pkin(6));
t8 = t13 * t33 + t31 * t36;
t42 = t14 * t35 - t8 * t32;
t41 = t14 * t32 + t8 * t35;
t11 = qJD(2) * t14;
t6 = qJD(5) * t8 - t11 * t36;
t7 = -t13 * t36 + t31 * t33;
t40 = t6 * t32 + t7 * t63;
t39 = -t6 * t35 + t7 * t64;
t12 = -t28 * t56 + t30 * t55;
t38 = 0.2e1 * t13 * t11 + 0.2e1 * t14 * t12;
t20 = qJD(5) * t44 + qJD(4);
t19 = t23 + t43;
t18 = t32 * t67 - t57;
t17 = t33 * t63 + t54;
t16 = t52 + t58;
t15 = t33 * t64 - t51;
t5 = -t11 * t33 - t13 * t65 + t31 * t67;
t4 = -t32 * t49 + t35 * t20 + (-t32 * t19 - t59) * qJD(6);
t3 = -t35 * t49 - t32 * t20 + (-t35 * t19 + t32 * t73) * qJD(6);
t2 = t42 * qJD(6) + t12 * t32 - t5 * t35;
t1 = -t41 * qJD(6) + t12 * t35 + t5 * t32;
t9 = [0, 0, 0, 0, t38, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t56, -t55 (-t11 * t30 + t12 * t28) * pkin(2), t11, t12, t14 * qJD(4) + t11 * t48 + t12 * t23, 0, 0, 0, 0, 0, t12 * t33 + t36 * t61, t12 * t36 - t33 * t61, 0, 0, 0, 0, 0 (-t7 * t32 * qJD(5) + t1) * t33 + (qJD(5) * t42 + t40) * t36 (-t7 * t66 - t2) * t33 + (-qJD(5) * t41 - t39) * t36; 0, 0, 0, 0, 0, 0, t75, t23 * t75, -0.2e1 * t50, 0.2e1 * t46, 0, 0, 0, 0.2e1 * qJD(4) * t33 + 0.2e1 * t23 * t65, 0.2e1 * qJD(4) * t36 - 0.2e1 * t23 * t67, -0.2e1 * t26 * t50 - 0.2e1 * t27 * t53, 0.2e1 * t27 * t47 + 0.4e1 * t32 * t45, -0.2e1 * t33 * t58 - 0.2e1 * t70 * t66, 0.2e1 * t32 * t46 - 0.2e1 * t33 * t57, 0.2e1 * t50, 0.2e1 * t19 * t51 + 0.2e1 * t4 * t33 + 0.2e1 * (-t27 * t63 + t32 * t50) * t22, -0.2e1 * t19 * t54 + 0.2e1 * t3 * t33 + 0.2e1 * (t27 * t64 + t45) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 * t63, t69 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, t39, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t65, 0, -t22 * t67, -t49, -t32 * t52 - t36 * t47, -0.4e1 * t36 * t53 + t71 * t67, t17, -t15, 0 (-t35 * t44 - t36 * t74) * qJD(6) + (t32 * t43 - t59) * qJD(5) (-t22 * t72 + t32 * t44) * qJD(6) + (-pkin(9) * t72 + (pkin(5) * t35 + t74) * t33) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t67, 0, 0, 0, 0, 0, t15, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t65, 0, 0, 0, 0, 0, -t16, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t53, -0.2e1 * t47, 0, 0, 0, t32 * t60, t35 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t18, t65, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t64, 0, -pkin(9) * t63, pkin(9) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
