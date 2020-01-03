% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:58
% EndTime: 2020-01-03 12:03:59
% DurationCPUTime: 0.50s
% Computational Cost: add. (688->80), mult. (1632->126), div. (0->0), fcn. (1558->8), ass. (0->79)
t58 = sin(pkin(9));
t59 = cos(pkin(9));
t83 = t58 ^ 2 + t59 ^ 2;
t99 = t83 * qJD(3);
t100 = 0.2e1 * t99;
t65 = cos(qJ(2));
t82 = pkin(1) * qJD(2);
t80 = t65 * t82;
t48 = qJD(3) + t80;
t98 = t83 * t48;
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t41 = t61 * t58 - t64 * t59;
t62 = sin(qJ(2));
t50 = t62 * pkin(1) + qJ(3);
t38 = (-pkin(7) - t50) * t58;
t55 = t59 * pkin(7);
t39 = t59 * t50 + t55;
t71 = t64 * t38 - t61 * t39;
t97 = -t71 * qJD(4) + t41 * t48;
t46 = (-pkin(7) - qJ(3)) * t58;
t47 = t59 * qJ(3) + t55;
t69 = t64 * t46 - t61 * t47;
t96 = t41 * qJD(3) - t69 * qJD(4);
t36 = t41 * qJD(4);
t95 = t36 * pkin(8);
t42 = t64 * t58 + t61 * t59;
t37 = t42 * qJD(4);
t94 = t37 * pkin(4);
t93 = t42 * pkin(8);
t92 = t65 * pkin(1);
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t25 = -t60 * t41 + t63 * t42;
t14 = t25 * qJD(5) - t60 * t36 + t63 * t37;
t24 = t63 * t41 + t60 * t42;
t51 = -t59 * pkin(3) - pkin(2);
t31 = t41 * pkin(4) + t51;
t29 = t31 - t92;
t52 = t62 * t82;
t30 = t52 + t94;
t91 = t29 * t14 + t30 * t24;
t13 = -t24 * qJD(5) - t63 * t36 - t60 * t37;
t90 = t29 * t13 + t30 * t25;
t89 = t51 * t36;
t88 = t51 * t37;
t45 = t51 - t92;
t87 = t45 * t37 + t41 * t52;
t86 = -t45 * t36 + t42 * t52;
t81 = pkin(4) * qJD(5);
t79 = t60 * t81;
t78 = t63 * t81;
t75 = t58 * t52;
t74 = t59 * t52;
t73 = t31 * t13 + t25 * t94;
t72 = t31 * t14 + t24 * t94;
t70 = -t61 * t38 - t64 * t39;
t68 = -t61 * t46 - t64 * t47;
t67 = t70 * qJD(4) - t42 * t48;
t66 = -t42 * qJD(3) + t68 * qJD(4);
t40 = t41 * pkin(8);
t35 = t37 * pkin(8);
t26 = -0.2e1 * t42 * t36;
t23 = -t40 - t68;
t22 = t69 - t93;
t21 = -t40 - t70;
t20 = t71 - t93;
t17 = 0.2e1 * t36 * t41 - 0.2e1 * t42 * t37;
t16 = t66 + t95;
t15 = -t35 - t96;
t12 = t67 + t95;
t11 = -t35 - t97;
t6 = 0.2e1 * t25 * t13;
t5 = -0.2e1 * t13 * t24 - 0.2e1 * t25 * t14;
t4 = -t60 * t15 + t63 * t16 + (-t22 * t60 - t23 * t63) * qJD(5);
t3 = -t63 * t15 - t60 * t16 + (-t22 * t63 + t23 * t60) * qJD(5);
t2 = -t60 * t11 + t63 * t12 + (-t20 * t60 - t21 * t63) * qJD(5);
t1 = -t63 * t11 - t60 * t12 + (-t20 * t63 + t21 * t60) * qJD(5);
t7 = [0, 0, 0, 0, -0.2e1 * t52, -0.2e1 * t80, -0.2e1 * t74, 0.2e1 * t75, 0.2e1 * t98, 0.2e1 * (-pkin(2) - t92) * t52 + 0.2e1 * t50 * t98, t26, t17, 0, 0, 0, 0.2e1 * t87, 0.2e1 * t86, t6, t5, 0, 0, 0, 0.2e1 * t91, 0.2e1 * t90; 0, 0, 0, 0, -t52, -t80, -t74, t75, t99 + t98, -pkin(2) * t52 + qJ(3) * t98 + t50 * t99, t26, t17, 0, 0, 0, t87 + t88, t86 - t89, t6, t5, 0, 0, 0, t72 + t91, t73 + t90; 0, 0, 0, 0, 0, 0, 0, 0, t100, qJ(3) * t100, t26, t17, 0, 0, 0, 0.2e1 * t88, -0.2e1 * t89, t6, t5, 0, 0, 0, 0.2e1 * t72, 0.2e1 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, t37, -t36, 0, 0, 0, 0, 0, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, 0, 0, 0, 0, 0, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, t67, t97, 0, 0, t13, -t14, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, t66, t96, 0, 0, t13, -t14, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t79, -0.2e1 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
