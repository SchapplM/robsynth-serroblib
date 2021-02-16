% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR14_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:03
% EndTime: 2021-01-15 12:17:06
% DurationCPUTime: 0.44s
% Computational Cost: add. (499->89), mult. (1086->172), div. (0->0), fcn. (959->6), ass. (0->65)
t34 = sin(qJ(3));
t36 = cos(qJ(3));
t63 = sin(pkin(8));
t49 = qJD(3) * t63;
t64 = cos(pkin(8));
t50 = qJD(3) * t64;
t14 = t34 * t49 - t36 * t50;
t15 = -t34 * t50 - t36 * t49;
t75 = (-t63 * t14 + t64 * t15) * pkin(3);
t53 = t64 * t36;
t18 = -t63 * t34 + t53;
t17 = t18 ^ 2;
t74 = 2 * qJD(2);
t73 = 2 * qJD(5);
t59 = t36 * qJD(3);
t37 = -pkin(1) - pkin(6);
t65 = qJ(4) - t37;
t13 = -t34 * qJD(4) - t65 * t59;
t60 = t34 * qJD(3);
t40 = -t36 * qJD(4) + t65 * t60;
t3 = t63 * t13 - t64 * t40;
t33 = sin(qJ(5));
t72 = t3 * t33;
t35 = cos(qJ(5));
t71 = t3 * t35;
t70 = t18 * t15;
t52 = t63 * t36;
t19 = t64 * t34 + t52;
t69 = t19 * t14;
t68 = t33 * t15;
t67 = t35 * t15;
t32 = t35 ^ 2;
t66 = t33 ^ 2 - t32;
t27 = t34 * pkin(3) + qJ(2);
t62 = qJD(5) * t33;
t61 = qJD(5) * t35;
t21 = pkin(3) * t59 + qJD(2);
t58 = qJ(2) * qJD(3);
t28 = -t64 * pkin(3) - pkin(4);
t57 = t28 * t73;
t56 = t33 * t61;
t55 = t19 ^ 2 + t17;
t54 = -0.4e1 * t18 * t33 * t35;
t51 = t66 * qJD(5);
t20 = t65 * t34;
t10 = -t64 * t20 - t65 * t52;
t8 = t19 * pkin(4) - t18 * pkin(7) + t27;
t48 = t35 * t10 + t33 * t8;
t47 = t33 * t10 - t35 * t8;
t46 = t69 - t70;
t26 = t63 * pkin(3) + pkin(7);
t45 = t14 * t26 + t15 * t28;
t44 = t18 * t28 - t19 * t26;
t7 = -t33 * t14 + t19 * t61;
t6 = t35 * t14 + t19 * t62;
t43 = -t18 * t61 - t68;
t42 = -t18 * t62 + t67;
t41 = 0.2e1 * t46;
t4 = t64 * t13 + t63 * t40;
t9 = -t63 * t20 + t65 * t53;
t38 = t10 * t14 + t9 * t15 + t3 * t18 - t4 * t19;
t5 = -t14 * pkin(4) - t15 * pkin(7) + t21;
t2 = -t48 * qJD(5) - t33 * t4 + t35 * t5;
t1 = t47 * qJD(5) - t33 * t5 - t35 * t4;
t11 = [0, 0, 0, 0, t74, qJ(2) * t74, -0.2e1 * t34 * t59, 0.2e1 * (t34 ^ 2 - t36 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t34 + 0.2e1 * t36 * t58, 0.2e1 * qJD(2) * t36 - 0.2e1 * t34 * t58, -0.2e1 * t27 * t14 + 0.2e1 * t21 * t19, 0.2e1 * t27 * t15 + 0.2e1 * t21 * t18, 0.2e1 * t38, 0.2e1 * t10 * t4 + 0.2e1 * t27 * t21 + 0.2e1 * t9 * t3, -0.2e1 * t17 * t56 + 0.2e1 * t32 * t70, t66 * t17 * t73 + t15 * t54, -0.2e1 * t18 * t6 + 0.2e1 * t19 * t67, -0.2e1 * t18 * t7 - 0.2e1 * t19 * t68, -0.2e1 * t69, 0.2e1 * t2 * t19 + 0.2e1 * t47 * t14 + 0.2e1 * t9 * t68 + 0.2e1 * (t9 * t61 + t72) * t18, 0.2e1 * t1 * t19 + 0.2e1 * t48 * t14 + 0.2e1 * t9 * t67 + 0.2e1 * (-t9 * t62 + t71) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t38, 0, 0, 0, 0, 0, t33 * t41 - t55 * t61, t35 * t41 + t55 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t46, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t59, 0, -t37 * t60, -t37 * t59, -t3, -t4, -t75, (-t64 * t3 + t63 * t4) * pkin(3), -t18 * t51 + t33 * t67, qJD(5) * t54 - t66 * t15, t7, -t6, 0, -t71 + t45 * t33 + (t33 * t9 + t44 * t35) * qJD(5), t72 + t45 * t35 + (-t44 * t33 + t35 * t9) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t59, t15, t14, 0, t75, 0, 0, 0, 0, 0, t42, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t56, -0.2e1 * t51, 0, 0, 0, t33 * t57, t35 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t15, 0, t21, 0, 0, 0, 0, 0, -t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t43, -t14, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t62, 0, -t26 * t61, t26 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
