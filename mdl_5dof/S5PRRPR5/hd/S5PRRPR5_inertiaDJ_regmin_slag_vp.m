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
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 16:27:54
% EndTime: 2019-12-05 16:27:55
% DurationCPUTime: 0.47s
% Computational Cost: add. (513->102), mult. (1421->222), div. (0->0), fcn. (1349->10), ass. (0->80)
t86 = 2 * qJD(5);
t37 = sin(pkin(10));
t44 = cos(qJ(3));
t41 = sin(qJ(3));
t75 = cos(pkin(10));
t61 = t75 * t41;
t26 = t37 * t44 + t61;
t40 = sin(qJ(5));
t85 = t26 * t40;
t43 = cos(qJ(5));
t84 = t26 * t43;
t83 = t37 * t41;
t38 = sin(pkin(5));
t42 = sin(qJ(2));
t82 = t38 * t42;
t45 = cos(qJ(2));
t81 = t38 * t45;
t21 = t26 * qJD(3);
t80 = t40 * t21;
t79 = t43 * t21;
t60 = t75 * t44;
t70 = t41 * qJD(3);
t22 = qJD(3) * t60 - t37 * t70;
t78 = t43 * t22;
t77 = -qJ(4) - pkin(7);
t36 = t43 ^ 2;
t76 = t40 ^ 2 - t36;
t74 = qJD(2) * t42;
t73 = qJD(2) * t45;
t72 = qJD(5) * t40;
t71 = qJD(5) * t43;
t69 = t44 * qJD(3);
t68 = -0.2e1 * pkin(2) * qJD(3);
t32 = -t75 * pkin(3) - pkin(4);
t67 = t32 * t86;
t34 = pkin(3) * t70;
t66 = t45 * t70;
t65 = t38 * t74;
t64 = t38 * t73;
t63 = t40 * t71;
t33 = -t44 * pkin(3) - pkin(2);
t62 = -0.4e1 * t40 * t84;
t59 = qJD(3) * t77;
t58 = t76 * qJD(5);
t25 = -t60 + t83;
t13 = t25 * pkin(4) - t26 * pkin(8) + t33;
t28 = t77 * t44;
t17 = -t75 * t28 + t77 * t83;
t57 = t43 * t13 - t40 * t17;
t56 = t40 * t13 + t43 * t17;
t31 = t37 * pkin(3) + pkin(8);
t55 = -t21 * t31 + t22 * t32;
t54 = t25 * t31 - t26 * t32;
t39 = cos(pkin(5));
t23 = t39 * t41 + t44 * t82;
t51 = t39 * t44 - t41 * t82;
t11 = t75 * t23 + t37 * t51;
t53 = t40 * t11 + t43 * t81;
t52 = -t43 * t11 + t40 * t81;
t50 = t25 * t71 + t80;
t49 = t40 * t22 + t26 * t71;
t48 = -t26 * t72 + t78;
t47 = -t41 * qJD(4) + t44 * t59;
t46 = -t23 * qJD(3) - t41 * t64;
t24 = t26 ^ 2;
t20 = t44 * qJD(4) + t41 * t59;
t16 = -t37 * t28 - t77 * t61;
t15 = t51 * qJD(3) + t44 * t64;
t12 = -t25 * t72 + t79;
t10 = t37 * t23 - t75 * t51;
t9 = t21 * pkin(4) - t22 * pkin(8) + t34;
t8 = t75 * t20 + t37 * t47;
t7 = t37 * t20 - t75 * t47;
t6 = t75 * t15 + t37 * t46;
t5 = t37 * t15 - t75 * t46;
t4 = t52 * qJD(5) - t40 * t6 + t43 * t65;
t3 = t53 * qJD(5) - t40 * t65 - t43 * t6;
t2 = -t56 * qJD(5) - t40 * t8 + t43 * t9;
t1 = -t57 * qJD(5) - t40 * t9 - t43 * t8;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t38 ^ 2 * t42 * t73 + 0.2e1 * t10 * t5 + 0.2e1 * t11 * t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t65, -t64, 0, 0, 0, 0, 0, (-t44 * t74 - t66) * t38, (t41 * t74 - t45 * t69) * t38, t10 * t22 - t11 * t21 - t6 * t25 + t5 * t26, t10 * t7 + t11 * t8 + t5 * t16 + t6 * t17 + (-pkin(3) * t66 + t33 * t74) * t38, 0, 0, 0, 0, 0, t49 * t10 - t53 * t21 + t4 * t25 + t5 * t85, t48 * t10 + t52 * t21 + t3 * t25 + t5 * t84; 0, 0, 0, 0, 0.2e1 * t41 * t69, 0.2e1 * (-t41 ^ 2 + t44 ^ 2) * qJD(3), 0, 0, 0, t41 * t68, t44 * t68, 0.2e1 * t16 * t22 - 0.2e1 * t17 * t21 - 0.2e1 * t8 * t25 + 0.2e1 * t7 * t26, 0.2e1 * t16 * t7 + 0.2e1 * t17 * t8 + 0.2e1 * t33 * t34, 0.2e1 * t36 * t26 * t22 - 0.2e1 * t24 * t63, t76 * t24 * t86 + t22 * t62, 0.2e1 * t48 * t25 + 0.2e1 * t26 * t79, -0.2e1 * t49 * t25 - 0.2e1 * t26 * t80, 0.2e1 * t25 * t21, 0.2e1 * t49 * t16 + 0.2e1 * t2 * t25 + 0.2e1 * t57 * t21 + 0.2e1 * t7 * t85, 0.2e1 * t1 * t25 + 0.2e1 * t48 * t16 - 0.2e1 * t56 * t21 + 0.2e1 * t7 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t15, 0, (t37 * t6 - t75 * t5) * pkin(3), 0, 0, 0, 0, 0, t10 * t72 - t5 * t43, t10 * t71 + t5 * t40; 0, 0, 0, 0, 0, 0, t69, -t70, 0, -pkin(7) * t69, pkin(7) * t70, (-t21 * t37 - t75 * t22) * pkin(3), (t37 * t8 - t75 * t7) * pkin(3), -t26 * t58 + t40 * t78, qJD(5) * t62 - t76 * t22, t50, t12, 0, -t7 * t43 + t55 * t40 + (t16 * t40 - t54 * t43) * qJD(5), t7 * t40 + t55 * t43 + (t16 * t43 + t54 * t40) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t63, -0.2e1 * t58, 0, 0, 0, t40 * t67, t43 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, t12, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t49, t21, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t72, 0, -t31 * t71, t31 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
