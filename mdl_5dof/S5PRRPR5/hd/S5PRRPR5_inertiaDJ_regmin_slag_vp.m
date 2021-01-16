% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:42
% EndTime: 2021-01-15 16:04:46
% DurationCPUTime: 0.48s
% Computational Cost: add. (549->110), mult. (1537->236), div. (0->0), fcn. (1453->10), ass. (0->80)
t86 = 2 * qJD(5);
t37 = sin(pkin(10));
t43 = cos(qJ(3));
t40 = sin(qJ(3));
t74 = cos(pkin(10));
t60 = t74 * t40;
t26 = t37 * t43 + t60;
t39 = sin(qJ(5));
t85 = t26 * t39;
t42 = cos(qJ(5));
t84 = t26 * t42;
t83 = t37 * t40;
t38 = sin(pkin(5));
t41 = sin(qJ(2));
t82 = t38 * t41;
t44 = cos(qJ(2));
t81 = t38 * t44;
t21 = t26 * qJD(3);
t80 = t39 * t21;
t79 = t42 * t21;
t59 = t74 * t43;
t69 = t40 * qJD(3);
t22 = qJD(3) * t59 - t37 * t69;
t78 = t42 * t22;
t77 = qJ(4) + pkin(7);
t36 = t42 ^ 2;
t76 = t39 ^ 2 - t36;
t75 = cos(pkin(5));
t73 = qJD(2) * t41;
t72 = qJD(2) * t44;
t71 = qJD(5) * t39;
t70 = qJD(5) * t42;
t68 = t43 * qJD(3);
t67 = -0.2e1 * pkin(2) * qJD(3);
t32 = -t74 * pkin(3) - pkin(4);
t66 = t32 * t86;
t34 = pkin(3) * t69;
t65 = t44 * t69;
t64 = t38 * t73;
t63 = t38 * t72;
t62 = t39 * t70;
t33 = -t43 * pkin(3) - pkin(2);
t61 = -0.4e1 * t39 * t84;
t58 = qJD(3) * t77;
t57 = t76 * qJD(5);
t25 = -t59 + t83;
t13 = t25 * pkin(4) - t26 * pkin(8) + t33;
t28 = t77 * t43;
t17 = t74 * t28 - t77 * t83;
t56 = t42 * t13 - t39 * t17;
t55 = t39 * t13 + t42 * t17;
t31 = t37 * pkin(3) + pkin(8);
t54 = -t21 * t31 + t22 * t32;
t53 = t25 * t31 - t26 * t32;
t23 = t75 * t40 + t43 * t82;
t47 = -t40 * t82 + t75 * t43;
t11 = t74 * t23 + t37 * t47;
t52 = t39 * t11 + t42 * t81;
t51 = -t42 * t11 + t39 * t81;
t50 = t25 * t70 + t80;
t49 = t39 * t22 + t26 * t70;
t48 = -t26 * t71 + t78;
t46 = -t40 * qJD(4) - t43 * t58;
t45 = -t23 * qJD(3) - t40 * t63;
t24 = t26 ^ 2;
t20 = t43 * qJD(4) - t40 * t58;
t16 = t37 * t28 + t77 * t60;
t15 = t47 * qJD(3) + t43 * t63;
t12 = -t25 * t71 + t79;
t10 = t37 * t23 - t74 * t47;
t9 = t21 * pkin(4) - t22 * pkin(8) + t34;
t8 = t74 * t20 + t37 * t46;
t7 = t37 * t20 - t74 * t46;
t6 = t74 * t15 + t37 * t45;
t5 = t37 * t15 - t74 * t45;
t4 = t51 * qJD(5) - t39 * t6 + t42 * t64;
t3 = t52 * qJD(5) - t39 * t64 - t42 * t6;
t2 = -t55 * qJD(5) - t39 * t8 + t42 * t9;
t1 = -t56 * qJD(5) - t39 * t9 - t42 * t8;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t38 ^ 2 * t41 * t72 + 0.2e1 * t10 * t5 + 0.2e1 * t11 * t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t64, -t63, 0, 0, 0, 0, 0, (-t43 * t73 - t65) * t38, (t40 * t73 - t44 * t68) * t38, (-t21 * t44 + t25 * t73) * t38, (-t22 * t44 + t26 * t73) * t38, t10 * t22 - t11 * t21 - t6 * t25 + t5 * t26, t10 * t7 + t11 * t8 + t5 * t16 + t6 * t17 + (-pkin(3) * t65 + t33 * t73) * t38, 0, 0, 0, 0, 0, t10 * t49 - t21 * t52 + t4 * t25 + t5 * t85, t10 * t48 + t21 * t51 + t3 * t25 + t5 * t84; 0, 0, 0, 0, 0.2e1 * t40 * t68, 0.2e1 * (-t40 ^ 2 + t43 ^ 2) * qJD(3), 0, 0, 0, t40 * t67, t43 * t67, 0.2e1 * t33 * t21 + 0.2e1 * t25 * t34, 0.2e1 * t33 * t22 + 0.2e1 * t26 * t34, 0.2e1 * t16 * t22 - 0.2e1 * t17 * t21 - 0.2e1 * t8 * t25 + 0.2e1 * t7 * t26, 0.2e1 * t16 * t7 + 0.2e1 * t17 * t8 + 0.2e1 * t33 * t34, 0.2e1 * t36 * t26 * t22 - 0.2e1 * t24 * t62, t24 * t76 * t86 + t22 * t61, 0.2e1 * t25 * t48 + 0.2e1 * t26 * t79, -0.2e1 * t25 * t49 - 0.2e1 * t26 * t80, 0.2e1 * t25 * t21, 0.2e1 * t16 * t49 + 0.2e1 * t2 * t25 + 0.2e1 * t21 * t56 + 0.2e1 * t7 * t85, 0.2e1 * t1 * t25 + 0.2e1 * t16 * t48 - 0.2e1 * t21 * t55 + 0.2e1 * t7 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t15, -t5, -t6, 0, (t37 * t6 - t5 * t74) * pkin(3), 0, 0, 0, 0, 0, t10 * t71 - t5 * t42, t10 * t70 + t5 * t39; 0, 0, 0, 0, 0, 0, t68, -t69, 0, -pkin(7) * t68, pkin(7) * t69, -t7, -t8, (-t21 * t37 - t22 * t74) * pkin(3), (t37 * t8 - t7 * t74) * pkin(3), -t26 * t57 + t39 * t78, qJD(5) * t61 - t22 * t76, t50, t12, 0, -t7 * t42 + t54 * t39 + (t16 * t39 - t42 * t53) * qJD(5), t7 * t39 + t54 * t42 + (t16 * t42 + t39 * t53) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t62, -0.2e1 * t57, 0, 0, 0, t39 * t66, t42 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22, 0, t34, 0, 0, 0, 0, 0, t12, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t49, t21, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t71, 0, -t31 * t70, t31 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
