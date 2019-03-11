% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPPRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:02
% EndTime: 2019-03-08 19:23:04
% DurationCPUTime: 0.56s
% Computational Cost: add. (327->100), mult. (931->212), div. (0->0), fcn. (859->10), ass. (0->82)
t40 = cos(qJ(6));
t31 = t40 ^ 2;
t37 = sin(qJ(6));
t86 = t37 ^ 2 - t31;
t58 = qJD(6) * t86;
t88 = 2 * qJD(3);
t33 = sin(pkin(11));
t41 = cos(qJ(5));
t87 = t33 * t41;
t35 = cos(pkin(11));
t43 = -pkin(2) - pkin(3);
t24 = t35 * qJ(3) + t33 * t43;
t38 = sin(qJ(5));
t30 = t38 ^ 2;
t85 = -t41 ^ 2 + t30;
t34 = sin(pkin(6));
t84 = qJD(2) * t34;
t83 = qJD(5) * t37;
t82 = qJD(5) * t38;
t81 = qJD(5) * t40;
t80 = qJD(5) * t41;
t79 = qJD(6) * t37;
t78 = qJD(6) * t38;
t77 = qJD(6) * t40;
t76 = qJD(6) * t41;
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t14 = (t33 * t39 + t35 * t42) * t34;
t75 = t14 * qJD(5);
t74 = t33 * qJD(3);
t73 = t35 * qJD(3);
t72 = t35 * qJD(5);
t71 = -0.2e1 * pkin(5) * qJD(6);
t70 = pkin(9) * t77;
t69 = t37 * t76;
t68 = t40 * t76;
t67 = t39 * t84;
t66 = t42 * t84;
t65 = t40 * t73;
t64 = t37 * t82;
t63 = t37 * t77;
t62 = t38 * t80;
t61 = t33 * t82;
t60 = t38 * t81;
t59 = t40 * t80;
t23 = -t33 * qJ(3) + t35 * t43;
t57 = t85 * qJD(5);
t56 = t38 * t59;
t21 = pkin(4) - t23;
t55 = t41 * pkin(5) + t38 * pkin(9);
t54 = -pkin(5) * t38 + pkin(9) * t41;
t22 = -pkin(8) + t24;
t53 = -t54 * qJD(5) + t22 * t76 - t74;
t15 = (-t33 * t42 + t35 * t39) * t34;
t36 = cos(pkin(6));
t8 = t15 * t41 - t36 * t38;
t52 = t14 * t40 - t8 * t37;
t51 = t14 * t37 + t8 * t40;
t7 = t15 * t38 + t36 * t41;
t12 = qJD(2) * t14;
t5 = t8 * qJD(5) + t12 * t38;
t50 = t5 * t37 + t7 * t77;
t49 = -t5 * t40 + t7 * t79;
t48 = t22 * t82 - t41 * t73;
t47 = -t22 * t80 - t38 * t73;
t17 = t37 * t78 - t59;
t19 = t37 * t80 + t38 * t77;
t46 = t30 * t79 - t56;
t45 = -t30 * t77 - t37 * t62;
t13 = t21 + t55;
t44 = -qJD(6) * t13 + t48;
t20 = -t64 + t68;
t18 = -t60 - t69;
t11 = t33 * t66 - t35 * t67;
t10 = t37 * t61 + (t37 * t35 - t40 * t87) * qJD(6);
t9 = t33 * t60 + (t40 * t35 + t37 * t87) * qJD(6);
t6 = -qJD(5) * t7 + t12 * t41;
t4 = t44 * t37 - t53 * t40;
t3 = t53 * t37 + t44 * t40;
t2 = t52 * qJD(6) + t11 * t37 + t6 * t40;
t1 = -t51 * qJD(6) + t11 * t40 - t6 * t37;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t14 * t11 + 0.2e1 * t15 * t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t67, -t66, -t67, t66 (qJD(3) * t39 + (-pkin(2) * t39 + qJ(3) * t42) * qJD(2)) * t34, t11, t12, -t11 * t23 + t12 * t24 + (t14 * t33 + t15 * t35) * qJD(3), 0, 0, 0, 0, 0, t11 * t41 - t38 * t75, -t11 * t38 - t41 * t75, 0, 0, 0, 0, 0 (-t7 * t83 + t1) * t41 + (-t52 * qJD(5) - t50) * t38 (-t7 * t81 - t2) * t41 + (t51 * qJD(5) + t49) * t38; 0, 0, 0, 0, 0, t88, qJ(3) * t88, 0.2e1 * t74, 0.2e1 * t73 (-t23 * t33 + t24 * t35) * t88, 0.2e1 * t62, -0.2e1 * t57, 0, 0, 0, -0.2e1 * t21 * t82 + 0.2e1 * t41 * t74, -0.2e1 * t21 * t80 - 0.2e1 * t38 * t74, -0.2e1 * t30 * t63 + 0.2e1 * t31 * t62, 0.2e1 * t30 * t58 - 0.4e1 * t37 * t56, 0.2e1 * t38 * t69 + 0.2e1 * t85 * t81, -0.2e1 * t37 * t57 + 0.2e1 * t38 * t68, -0.2e1 * t62, -0.2e1 * t30 * t37 * t73 - 0.2e1 * t13 * t60 + 0.2e1 * t45 * t22 + 0.2e1 * t4 * t41, 0.2e1 * t13 * t64 + 0.2e1 * t46 * t22 + 0.2e1 * t3 * t41 - 0.2e1 * t30 * t65; 0, 0, 0, 0, 0, 0, t67, 0, 0, -t11 * t35 + t12 * t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t72, t41 * t72, 0, 0, 0, 0, 0, t10 * t41 + t45 * t33 + t35 * t60, t46 * t33 - t35 * t64 + t9 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, t49, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t82, 0, t47, t48, -t37 * t59 + t38 * t58, 0.4e1 * t38 * t63 + t86 * t80, t20, t18, 0 (-t70 + (pkin(5) * t37 - t22 * t40) * qJD(5)) * t41 + (pkin(9) * t83 - t65 + (pkin(5) * t40 + t22 * t37) * qJD(6)) * t38 (t55 * qJD(5) + t22 * t78) * t40 + (t54 * qJD(6) - t47) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33 * t80, t61, 0, 0, 0, 0, 0, t17 * t33, t19 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t80, 0, 0, 0, 0, 0, t18, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t63, -0.2e1 * t58, 0, 0, 0, t37 * t71, t40 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t19, -t82, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t79, 0, -t70, pkin(9) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t16;
