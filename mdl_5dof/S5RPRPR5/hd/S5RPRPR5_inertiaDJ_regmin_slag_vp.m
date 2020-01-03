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
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:43:02
% EndTime: 2020-01-03 11:43:06
% DurationCPUTime: 0.52s
% Computational Cost: add. (761->88), mult. (1942->189), div. (0->0), fcn. (1778->8), ass. (0->69)
t60 = sin(qJ(3));
t56 = sin(pkin(8));
t58 = cos(pkin(8));
t45 = -t58 * pkin(2) - t56 * pkin(6) - pkin(1);
t79 = qJ(4) * t56;
t67 = -t45 + t79;
t62 = cos(qJ(3));
t71 = t62 * t58 * qJ(2);
t87 = t67 * t60 - t71;
t86 = 0.2e1 * t58;
t55 = sin(pkin(9));
t85 = pkin(3) * t55;
t84 = t56 * t60;
t57 = cos(pkin(9));
t83 = t57 * t62;
t73 = qJD(4) * t56;
t80 = qJ(2) * t60;
t74 = qJD(3) * t62;
t77 = qJD(2) * t62;
t82 = t45 * t74 + t58 * t77;
t19 = -t60 * t73 + (-t58 * t80 - t62 * t79) * qJD(3) + t82;
t78 = qJD(2) * t60;
t69 = t58 * t78;
t20 = t87 * qJD(3) - t62 * t73 - t69;
t7 = t57 * t19 + t55 * t20;
t26 = -t67 * t62 + (-pkin(3) - t80) * t58;
t15 = t55 * t26 - t57 * t87;
t70 = t56 * t74;
t42 = pkin(3) * t70 + t56 * qJD(2);
t81 = pkin(3) * t84 + t56 * qJ(2);
t76 = qJD(3) * t56;
t75 = qJD(3) * t60;
t72 = qJ(2) * qJD(3);
t68 = t60 * t72;
t6 = -t55 * t19 + t57 * t20;
t14 = t57 * t26 + t55 * t87;
t66 = t76 * t86;
t53 = t56 ^ 2;
t65 = 0.2e1 * (t58 ^ 2 + t53) * qJD(2);
t44 = t55 * t62 + t57 * t60;
t37 = t44 * t56;
t38 = -t55 * t84 + t56 * t83;
t59 = sin(qJ(5));
t61 = cos(qJ(5));
t64 = -t61 * t37 - t59 * t38;
t22 = -t59 * t37 + t61 * t38;
t63 = t55 * t60 - t83;
t10 = -t58 * pkin(4) - t38 * pkin(7) + t14;
t11 = -t37 * pkin(7) + t15;
t34 = t44 * t76;
t4 = t34 * pkin(7) + t6;
t33 = t63 * t76;
t5 = t33 * pkin(7) + t7;
t2 = -t59 * t5 + t61 * t4 + (-t10 * t59 - t11 * t61) * qJD(5);
t1 = -t59 * t4 - t61 * t5 + (-t10 * t61 + t11 * t59) * qJD(5);
t50 = t57 * pkin(3) + pkin(4);
t40 = t63 * qJD(3);
t39 = t44 * qJD(3);
t36 = (-t50 * t59 - t61 * t85) * qJD(5);
t35 = (-t50 * t61 + t59 * t85) * qJD(5);
t30 = -t69 + (-t60 * t45 - t71) * qJD(3);
t29 = t58 * t68 - t82;
t27 = t37 * pkin(4) + t81;
t23 = -t33 * pkin(4) + t42;
t13 = -t61 * t39 + t59 * t40 + (-t44 * t61 + t59 * t63) * qJD(5);
t12 = t59 * t39 + t61 * t40 + (t44 * t59 + t61 * t63) * qJD(5);
t9 = t22 * qJD(5) - t61 * t33 - t59 * t34;
t8 = t64 * qJD(5) + t59 * t33 - t61 * t34;
t3 = [0, 0, 0, 0, 0, t65, qJ(2) * t65, -0.2e1 * t53 * t60 * t74, 0.2e1 * (t60 ^ 2 - t62 ^ 2) * t53 * qJD(3), t60 * t66, t62 * t66, 0, -0.2e1 * t30 * t58 + 0.2e1 * (t62 * t72 + t78) * t53, -0.2e1 * t29 * t58 + 0.2e1 * (-t68 + t77) * t53, 0.2e1 * t14 * t34 + 0.2e1 * t15 * t33 - 0.2e1 * t7 * t37 - 0.2e1 * t6 * t38, 0.2e1 * t14 * t6 + 0.2e1 * t15 * t7 + 0.2e1 * t81 * t42, 0.2e1 * t22 * t8, -0.2e1 * t22 * t9 + 0.2e1 * t64 * t8, -0.2e1 * t8 * t58, t9 * t86, 0, -0.2e1 * t2 * t58 - 0.2e1 * t23 * t64 + 0.2e1 * t27 * t9, -0.2e1 * t1 * t58 + 0.2e1 * t23 * t22 + 0.2e1 * t27 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t75, t58 * t74, t44 * t33 - t34 * t63 + t40 * t37 + t39 * t38, -t14 * t39 - t15 * t40 + t7 * t44 - t6 * t63, 0, 0, 0, 0, 0, -t13 * t58, -t12 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t39 * t63 - 0.2e1 * t44 * t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56 * t75, -t70, 0, t30, t29, (t33 * t55 + t34 * t57) * pkin(3), (t55 * t7 + t57 * t6) * pkin(3), 0, 0, t8, -t9, 0, -t36 * t58 + t2, -t35 * t58 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t74, 0, (-t39 * t57 - t40 * t55) * pkin(3), 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t36, 0.2e1 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
