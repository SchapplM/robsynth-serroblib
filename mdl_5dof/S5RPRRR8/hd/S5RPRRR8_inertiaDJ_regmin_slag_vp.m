% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:06
% EndTime: 2019-12-31 19:06:07
% DurationCPUTime: 0.40s
% Computational Cost: add. (404->75), mult. (882->134), div. (0->0), fcn. (736->6), ass. (0->69)
t37 = sin(qJ(5));
t38 = sin(qJ(4));
t40 = cos(qJ(5));
t41 = cos(qJ(4));
t21 = t37 * t41 + t40 * t38;
t75 = qJD(4) + qJD(5);
t10 = t75 * t21;
t67 = t37 * t38;
t20 = -t40 * t41 + t67;
t59 = t41 * qJD(4);
t61 = qJD(5) * t40;
t9 = -t40 * t59 - t41 * t61 + t75 * t67;
t44 = 0.2e1 * t21 * t10 - 0.2e1 * t9 * t20;
t74 = 2 * qJD(2);
t73 = pkin(7) + pkin(8);
t72 = t41 * pkin(4);
t71 = t9 * t21;
t39 = sin(qJ(3));
t42 = cos(qJ(3));
t43 = -pkin(1) - pkin(2);
t47 = t42 * qJ(2) + t39 * t43;
t23 = -pkin(7) + t47;
t70 = pkin(8) - t23;
t15 = t39 * qJD(2) + t47 * qJD(3);
t69 = t15 * t38;
t68 = t15 * t41;
t62 = qJD(3) * t42;
t63 = qJD(3) * t39;
t14 = qJ(2) * t63 - t42 * qJD(2) - t43 * t62;
t66 = t38 * t14;
t65 = t41 * t14;
t22 = t39 * qJ(2) - t42 * t43 + pkin(3);
t16 = t22 + t72;
t34 = -pkin(3) - t72;
t64 = t16 - t34;
t60 = t38 * qJD(4);
t58 = -0.2e1 * t71;
t57 = -0.2e1 * pkin(3) * qJD(4);
t56 = pkin(4) * t60;
t55 = qJD(5) * t37 * pkin(4);
t54 = pkin(4) * t61;
t53 = t38 * t59;
t52 = qJD(4) * t73;
t51 = qJD(4) * (pkin(3) + t22);
t50 = qJD(4) * t70;
t11 = t15 - t56;
t49 = t11 - t56;
t46 = t21 * t63 + t42 * t9;
t45 = -t42 * t10 + t20 * t63;
t29 = 0.2e1 * t53;
t28 = t73 * t41;
t27 = t73 * t38;
t26 = (-t38 ^ 2 + t41 ^ 2) * qJD(4);
t25 = t41 * t52;
t24 = t38 * t52;
t19 = 0.2e1 * t26;
t18 = t41 * t63 + t42 * t60;
t17 = t38 * t63 - t42 * t59;
t13 = t70 * t41;
t12 = t70 * t38;
t8 = t41 * t50 + t66;
t7 = t38 * t50 - t65;
t6 = t37 * t24 - t40 * t25 + (t27 * t37 - t28 * t40) * qJD(5);
t5 = t40 * t24 + t37 * t25 + (t27 * t40 + t28 * t37) * qJD(5);
t4 = t75 * t39 * t20 - t21 * t62;
t3 = t10 * t39 + t20 * t62;
t2 = -t37 * t7 + t40 * t8 + (-t12 * t37 + t13 * t40) * qJD(5);
t1 = -t37 * t8 - t40 * t7 + (-t12 * t40 - t13 * t37) * qJD(5);
t30 = [0, 0, 0, 0, t74, qJ(2) * t74, 0, 0.2e1 * t15, -0.2e1 * t14, t29, t19, 0, 0, 0, -0.2e1 * t22 * t60 + 0.2e1 * t68, -0.2e1 * t22 * t59 - 0.2e1 * t69, t58, -t44, 0, 0, 0, -0.2e1 * t16 * t10 - 0.2e1 * t11 * t20, -0.2e1 * t11 * t21 + 0.2e1 * t16 * t9; 0, 0, 0, 0, 0, 0, 0, t63, t62, 0, 0, 0, 0, 0, t18, -t17, 0, 0, 0, 0, 0, -t45, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t15, t14, -0.2e1 * t53, -0.2e1 * t26, 0, 0, 0, t38 * t51 - t68, t41 * t51 + t69, 0.2e1 * t71, t44, 0, 0, 0, t64 * t10 + t49 * t20, t49 * t21 - t64 * t9; 0, 0, 0, 0, 0, 0, 0, -t63, -t62, 0, 0, 0, 0, 0, -t18, t17, 0, 0, 0, 0, 0, t45, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t19, 0, 0, 0, t38 * t57, t41 * t57, t58, -t44, 0, 0, 0, 0.2e1 * t34 * t10 + 0.2e1 * t20 * t56, 0.2e1 * t21 * t56 - 0.2e1 * t34 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t60, 0, -t23 * t59 + t66, t23 * t60 + t65, 0, 0, t9, t10, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38 * t62 - t39 * t59, t39 * t60 - t41 * t62, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t60, 0, -pkin(7) * t59, pkin(7) * t60, 0, 0, -t9, -t10, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t55, -0.2e1 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t30;
