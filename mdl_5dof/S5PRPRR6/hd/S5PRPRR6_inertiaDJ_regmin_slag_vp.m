% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:54
% EndTime: 2019-12-05 15:57:56
% DurationCPUTime: 0.48s
% Computational Cost: add. (397->89), mult. (1148->184), div. (0->0), fcn. (1142->10), ass. (0->76)
t42 = cos(qJ(5));
t34 = t42 ^ 2;
t39 = sin(qJ(5));
t70 = t39 ^ 2 - t34;
t58 = t70 * qJD(5);
t37 = cos(pkin(10));
t30 = -t37 * pkin(3) - pkin(2);
t83 = 0.2e1 * t30;
t35 = sin(pkin(10));
t72 = pkin(7) + qJ(3);
t26 = t72 * t35;
t27 = t72 * t37;
t40 = sin(qJ(4));
t43 = cos(qJ(4));
t15 = -t40 * t26 + t43 * t27;
t24 = t43 * t35 + t40 * t37;
t8 = t24 * qJD(3) + t15 * qJD(4);
t82 = t8 * t39;
t73 = t43 * t37;
t23 = t40 * t35 - t73;
t19 = t23 * qJD(4);
t81 = t24 * t19;
t80 = t24 * t42;
t36 = sin(pkin(5));
t41 = sin(qJ(2));
t79 = t36 * t41;
t44 = cos(qJ(2));
t78 = t36 * t44;
t77 = t39 * t19;
t20 = t24 * qJD(4);
t76 = t39 * t20;
t75 = t42 * t19;
t74 = t42 * t20;
t71 = t35 ^ 2 + t37 ^ 2;
t69 = qJD(2) * t36;
t68 = qJD(2) * t41;
t67 = qJD(4) * t43;
t66 = qJD(5) * t39;
t65 = qJD(5) * t42;
t64 = -0.2e1 * pkin(4) * qJD(5);
t63 = t36 * t68;
t62 = t44 * t69;
t61 = t39 * t65;
t60 = t71 * t44;
t59 = -0.4e1 * t39 * t80;
t57 = 0.2e1 * t71 * qJD(3);
t56 = pkin(4) * t19 - pkin(8) * t20;
t55 = pkin(4) * t24 + pkin(8) * t23;
t13 = t23 * pkin(4) - t24 * pkin(8) + t30;
t54 = t42 * t13 - t39 * t15;
t53 = t39 * t13 + t42 * t15;
t38 = cos(pkin(5));
t17 = -t35 * t79 + t38 * t37;
t18 = t38 * t35 + t37 * t79;
t52 = -t17 * t35 + t18 * t37;
t10 = t40 * t17 + t43 * t18;
t6 = t10 * qJD(4) + t24 * t62;
t9 = -t43 * t17 + t40 * t18;
t51 = t6 * t39 + t9 * t65;
t50 = -t6 * t42 + t9 * t66;
t49 = t39 * t10 + t42 * t78;
t48 = -t42 * t10 + t39 * t78;
t47 = t24 * t65 - t77;
t46 = -t24 * t66 - t75;
t45 = t23 * t65 + t76;
t22 = t24 ^ 2;
t14 = t43 * t26 + t40 * t27;
t12 = t20 * pkin(4) + t19 * pkin(8);
t11 = -t23 * t66 + t74;
t7 = t26 * t67 - qJD(3) * t73 + (qJD(3) * t35 + qJD(4) * t27) * t40;
t5 = -t17 * t67 - t62 * t73 + (qJD(4) * t18 + t35 * t62) * t40;
t4 = t48 * qJD(5) + t39 * t5 + t42 * t63;
t3 = t49 * qJD(5) - t39 * t63 + t42 * t5;
t2 = -t53 * qJD(5) + t42 * t12 + t39 * t7;
t1 = -t54 * qJD(5) - t39 * t12 + t42 * t7;
t16 = [0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t52 - t79) * t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t63, -t62, -t37 * t63, t35 * t63, t60 * t69, t52 * qJD(3) + (-pkin(2) * t41 + qJ(3) * t60) * t69, 0, 0, 0, 0, 0, (-t20 * t44 + t23 * t68) * t36, (t19 * t44 + t24 * t68) * t36, 0, 0, 0, 0, 0, -t49 * t20 + t4 * t23 + t51 * t24 - t9 * t77, t48 * t20 + t3 * t23 - t50 * t24 - t9 * t75; 0, 0, 0, 0, 0, 0, t57, qJ(3) * t57, -0.2e1 * t81, 0.2e1 * t19 * t23 - 0.2e1 * t24 * t20, 0, 0, 0, t20 * t83, -t19 * t83, -0.2e1 * t22 * t61 - 0.2e1 * t34 * t81, -t19 * t59 + 0.2e1 * t22 * t58, 0.2e1 * t46 * t23 + 0.2e1 * t24 * t74, -0.2e1 * t47 * t23 - 0.2e1 * t24 * t76, 0.2e1 * t23 * t20, 0.2e1 * t47 * t14 + 0.2e1 * t2 * t23 + 0.2e1 * t54 * t20 + 0.2e1 * t24 * t82, 0.2e1 * t1 * t23 + 0.2e1 * t46 * t14 - 0.2e1 * t53 * t20 + 0.2e1 * t8 * t80; 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, 0, 0, 0, 0, t11, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, t50, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, -t8, t7, -t24 * t58 - t39 * t75, qJD(5) * t59 + t70 * t19, t45, t11, 0, -t8 * t42 + t56 * t39 + (t14 * t39 - t55 * t42) * qJD(5), t82 + t56 * t42 + (t14 * t42 + t55 * t39) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t61, -0.2e1 * t58, 0, 0, 0, t39 * t64, t42 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t47, t20, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t66, 0, -pkin(8) * t65, pkin(8) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t16;
