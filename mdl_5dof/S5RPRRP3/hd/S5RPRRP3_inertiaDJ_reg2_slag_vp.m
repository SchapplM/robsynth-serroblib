% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:54
% EndTime: 2019-12-05 18:03:57
% DurationCPUTime: 0.60s
% Computational Cost: add. (770->89), mult. (1699->153), div. (0->0), fcn. (1474->6), ass. (0->60)
t78 = qJD(3) + qJD(4);
t77 = 2 * qJD(3);
t45 = sin(qJ(4));
t76 = pkin(3) * t45;
t46 = sin(qJ(3));
t47 = cos(qJ(3));
t74 = cos(qJ(4));
t35 = t45 * t47 + t74 * t46;
t21 = t78 * t35;
t75 = t21 * pkin(4);
t57 = t74 * t47;
t71 = t45 * t46;
t34 = -t57 + t71;
t73 = t34 * t21;
t56 = t74 * qJD(4);
t20 = -qJD(3) * t57 - t47 * t56 + t78 * t71;
t72 = t35 * t20;
t55 = pkin(3) * t56;
t69 = -t21 * t76 - t34 * t55;
t68 = -t20 * t76 + t35 * t55;
t54 = sin(pkin(8)) * pkin(1) + pkin(6);
t53 = pkin(7) + t54;
t51 = t45 * t53;
t29 = t46 * t51;
t33 = t53 * t47;
t16 = t74 * t33 - t29;
t67 = qJD(4) * t45;
t66 = t46 * qJD(3);
t65 = t47 * qJD(3);
t41 = -cos(pkin(8)) * pkin(1) - pkin(2);
t64 = t41 * t77;
t63 = pkin(3) * t66;
t62 = pkin(3) * t67;
t61 = t74 * pkin(3);
t60 = t34 * t67;
t59 = t35 * t67;
t58 = t46 * t65;
t49 = t53 * t74;
t30 = t46 * t49;
t15 = -t45 * t33 - t30;
t36 = -t47 * pkin(3) + t41;
t48 = qJD(3) * t49;
t50 = qJD(3) * t51;
t6 = qJD(4) * t30 + t33 * t67 + t46 * t48 + t47 * t50;
t52 = t54 * qJD(3);
t7 = qJD(4) * t29 - t33 * t56 + t46 * t50 - t47 * t48;
t43 = t61 + pkin(4);
t40 = -0.2e1 * t55;
t39 = -0.2e1 * t62;
t22 = t34 * pkin(4) + t36;
t14 = t63 + t75;
t11 = -0.2e1 * t72;
t10 = 0.2e1 * t73;
t9 = -t34 * qJ(5) + t16;
t8 = -t35 * qJ(5) + t15;
t5 = -0.2e1 * t72 + 0.2e1 * t73;
t4 = 0.2e1 * t34 * t20 - 0.2e1 * t35 * t21;
t3 = t20 * qJ(5) - t35 * qJD(5) + t7;
t2 = t21 * qJ(5) + t34 * qJD(5) + t6;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t58, (-t46 ^ 2 + t47 ^ 2) * t77, 0, -0.2e1 * t58, 0, 0, t46 * t64, t47 * t64, 0, 0, t11, t4, 0, t10, 0, 0, 0.2e1 * t36 * t21 + 0.2e1 * t34 * t63, -0.2e1 * t36 * t20 + 0.2e1 * t35 * t63, 0.2e1 * t15 * t20 - 0.2e1 * t16 * t21 + 0.2e1 * t6 * t34 - 0.2e1 * t7 * t35, 0.2e1 * t15 * t7 - 0.2e1 * t16 * t6 + 0.2e1 * t36 * t63, t11, t4, 0, t10, 0, 0, 0.2e1 * t14 * t34 + 0.2e1 * t22 * t21, 0.2e1 * t14 * t35 - 0.2e1 * t22 * t20, 0.2e1 * t2 * t34 + 0.2e1 * t8 * t20 - 0.2e1 * t9 * t21 - 0.2e1 * t3 * t35, 0.2e1 * t22 * t14 - 0.2e1 * t9 * t2 + 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t21 - t16 * t20 - t7 * t34 - t6 * t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t35 - t9 * t20 - t8 * t21 - t3 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, -t66, 0, -t47 * t52, t46 * t52, 0, 0, 0, 0, -t20, 0, -t21, 0, t7, t6, (t74 * t20 + t59) * pkin(3) + t69, (t74 * t7 - t45 * t6 + (-t15 * t45 + t74 * t16) * qJD(4)) * pkin(3), 0, 0, -t20, 0, -t21, 0, t3, t2, pkin(3) * t59 + t43 * t20 + t69, t3 * t43 + (-t2 * t45 + (-t45 * t8 + t74 * t9) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t20, 0, (-t74 * t21 + t60) * pkin(3) + t68, 0, 0, 0, 0, 0, 0, -t21, t20, 0, pkin(3) * t60 - t21 * t43 + t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t40, 0, 0, 0, 0, 0, 0, 0, 0, t39, t40, 0, 0.2e1 * (t61 - t43) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, -t21, 0, t7, t6, 0, 0, 0, 0, -t20, 0, -t21, 0, t3, t2, t20 * pkin(4), t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t20, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t20, 0, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t55, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t55, 0, -pkin(4) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
