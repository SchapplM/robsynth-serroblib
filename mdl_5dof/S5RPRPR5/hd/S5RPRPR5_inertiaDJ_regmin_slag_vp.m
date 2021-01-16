% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:56
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:55:12
% EndTime: 2021-01-15 11:55:16
% DurationCPUTime: 0.57s
% Computational Cost: add. (832->98), mult. (2138->203), div. (0->0), fcn. (1942->8), ass. (0->69)
t60 = sin(qJ(3));
t56 = sin(pkin(8));
t58 = cos(pkin(8));
t44 = -pkin(2) * t58 - t56 * pkin(6) - pkin(1);
t77 = qJ(4) * t56;
t65 = -t44 + t77;
t62 = cos(qJ(3));
t70 = t62 * t58 * qJ(2);
t84 = t65 * t60 - t70;
t83 = 2 * qJD(3);
t55 = sin(pkin(9));
t82 = pkin(3) * t55;
t81 = t56 * t60;
t57 = cos(pkin(9));
t80 = t57 * t62;
t78 = qJ(2) * t60;
t24 = -t65 * t62 + (-pkin(3) - t78) * t58;
t15 = t55 * t24 - t57 * t84;
t73 = qJD(3) * t62;
t75 = qJD(2) * t62;
t79 = t44 * t73 + t58 * t75;
t68 = t56 * t73;
t40 = pkin(3) * t68 + t56 * qJD(2);
t43 = pkin(3) * t81 + t56 * qJ(2);
t76 = qJD(2) * t60;
t74 = qJD(3) * t60;
t72 = qJD(4) * t56;
t71 = qJ(2) * qJD(3);
t69 = t56 * t74;
t67 = t58 * t76;
t66 = t60 * t71;
t17 = -t60 * t72 + (-t58 * t78 - t62 * t77) * qJD(3) + t79;
t18 = t84 * qJD(3) - t62 * t72 - t67;
t6 = -t55 * t17 + t57 * t18;
t14 = t57 * t24 + t55 * t84;
t64 = t56 * t58 * t83;
t53 = t56 ^ 2;
t63 = 0.2e1 * (t58 ^ 2 + t53) * qJD(2);
t7 = t57 * t17 + t55 * t18;
t42 = t55 * t62 + t57 * t60;
t35 = t42 * t56;
t36 = -t55 * t81 + t56 * t80;
t59 = sin(qJ(5));
t61 = cos(qJ(5));
t19 = t61 * t35 + t59 * t36;
t20 = -t59 * t35 + t61 * t36;
t41 = -t55 * t60 + t80;
t10 = -t58 * pkin(4) - t36 * pkin(7) + t14;
t11 = -t35 * pkin(7) + t15;
t32 = qJD(3) * t35;
t4 = t32 * pkin(7) + t6;
t31 = t55 * t69 - t57 * t68;
t5 = t31 * pkin(7) + t7;
t2 = -t59 * t5 + t61 * t4 + (-t10 * t59 - t11 * t61) * qJD(5);
t1 = -t59 * t4 - t61 * t5 + (-t10 * t61 + t11 * t59) * qJD(5);
t50 = t57 * pkin(3) + pkin(4);
t38 = t41 * qJD(3);
t37 = t42 * qJD(3);
t34 = (-t50 * t59 - t61 * t82) * qJD(5);
t33 = (-t50 * t61 + t59 * t82) * qJD(5);
t28 = -t67 + (-t60 * t44 - t70) * qJD(3);
t27 = t58 * t66 - t79;
t25 = t35 * pkin(4) + t43;
t21 = -t31 * pkin(4) + t40;
t13 = -t61 * t37 - t59 * t38 + (-t41 * t59 - t42 * t61) * qJD(5);
t12 = t59 * t37 - t61 * t38 + (-t41 * t61 + t42 * t59) * qJD(5);
t9 = t20 * qJD(5) - t61 * t31 - t59 * t32;
t8 = -t19 * qJD(5) + t59 * t31 - t61 * t32;
t3 = [0, 0, 0, 0, t63, qJ(2) * t63, -0.2e1 * t53 * t60 * t73, (t60 ^ 2 - t62 ^ 2) * t53 * t83, t60 * t64, t62 * t64, 0, -0.2e1 * t28 * t58 + 0.2e1 * (t62 * t71 + t76) * t53, -0.2e1 * t27 * t58 + 0.2e1 * (-t66 + t75) * t53, -0.2e1 * t43 * t31 + 0.2e1 * t40 * t35 - 0.2e1 * t6 * t58, -0.2e1 * t43 * t32 + 0.2e1 * t40 * t36 + 0.2e1 * t7 * t58, 0.2e1 * t14 * t32 + 0.2e1 * t15 * t31 - 0.2e1 * t7 * t35 - 0.2e1 * t6 * t36, 0.2e1 * t14 * t6 + 0.2e1 * t15 * t7 + 0.2e1 * t43 * t40, 0.2e1 * t20 * t8, -0.2e1 * t8 * t19 - 0.2e1 * t20 * t9, -0.2e1 * t8 * t58, 0.2e1 * t9 * t58, 0, 0.2e1 * t21 * t19 - 0.2e1 * t2 * t58 + 0.2e1 * t25 * t9, -0.2e1 * t1 * t58 + 0.2e1 * t21 * t20 + 0.2e1 * t25 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t74, t58 * t73, t37 * t58, t38 * t58, t42 * t31 + t41 * t32 - t38 * t35 + t37 * t36, -t14 * t37 + t15 * t38 + t6 * t41 + t7 * t42, 0, 0, 0, 0, 0, -t13 * t58, -t12 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t41 * t37 + 0.2e1 * t42 * t38, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t68, 0, t28, t27, t6, -t7, (t31 * t55 + t32 * t57) * pkin(3), (t55 * t7 + t57 * t6) * pkin(3), 0, 0, t8, -t9, 0, -t34 * t58 + t2, -t33 * t58 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t73, -t37, -t38, 0, (-t37 * t57 + t38 * t55) * pkin(3), 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t34, 0.2e1 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t32, 0, t40, 0, 0, 0, 0, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
