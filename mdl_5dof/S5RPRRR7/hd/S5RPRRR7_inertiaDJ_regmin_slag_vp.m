% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:09
% EndTime: 2019-12-31 19:04:12
% DurationCPUTime: 0.76s
% Computational Cost: add. (666->128), mult. (1725->237), div. (0->0), fcn. (1457->8), ass. (0->93)
t57 = sin(qJ(3));
t110 = -0.4e1 * t57;
t46 = -cos(pkin(9)) * pkin(1) - pkin(2);
t60 = cos(qJ(3));
t68 = -t60 * pkin(3) - t57 * pkin(7);
t31 = t46 + t68;
t56 = sin(qJ(4));
t22 = t56 * t31;
t45 = sin(pkin(9)) * pkin(1) + pkin(6);
t59 = cos(qJ(4));
t99 = t59 * t60;
t35 = t45 * t99;
t109 = -t22 - t35;
t90 = t60 * qJD(3);
t76 = t59 * t90;
t94 = qJD(4) * t56;
t108 = -t57 * t94 + t76;
t51 = t57 ^ 2;
t70 = (-t60 ^ 2 + t51) * qJD(3);
t52 = t59 ^ 2;
t96 = t56 ^ 2 - t52;
t71 = t96 * qJD(4);
t107 = qJD(4) + qJD(5);
t49 = t57 * qJD(3);
t92 = qJD(4) * t60;
t82 = t56 * t92;
t26 = t59 * t49 + t82;
t67 = pkin(3) * t57 - pkin(7) * t60;
t38 = t67 * qJD(3);
t93 = qJD(4) * t59;
t6 = t26 * t45 - t31 * t93 - t56 * t38;
t106 = pkin(7) + pkin(8);
t80 = t56 * t90;
t61 = t57 * t93 + t80;
t5 = -t61 * pkin(8) - t6;
t58 = cos(qJ(5));
t105 = t58 * t5;
t104 = t45 * t56;
t102 = t56 * t57;
t14 = -pkin(8) * t102 - t109;
t55 = sin(qJ(5));
t103 = t55 * t14;
t101 = t57 * t59;
t100 = t58 * t14;
t77 = t45 * t49;
t97 = t59 * t38 + t56 * t77;
t91 = qJD(5) * t55;
t89 = -0.2e1 * pkin(3) * qJD(4);
t88 = 0.2e1 * qJD(3) * t46;
t87 = pkin(4) * t94;
t86 = pkin(4) * t49;
t85 = pkin(4) * t91;
t84 = qJD(5) * t58 * pkin(4);
t81 = t59 * t92;
t79 = t56 * t93;
t78 = t57 * t90;
t75 = t45 * t90;
t4 = (pkin(4) * t57 - pkin(8) * t99) * qJD(3) + (-t35 + (pkin(8) * t57 - t31) * t56) * qJD(4) + t97;
t74 = t58 * t4 - t55 * t5;
t23 = t59 * t31;
t13 = -pkin(8) * t101 + t23 + (-pkin(4) - t104) * t60;
t73 = t60 * pkin(4) - t13;
t72 = qJD(4) * t106;
t69 = t56 * t76;
t66 = t58 * t13 - t103;
t65 = t55 * t13 + t100;
t42 = t106 * t56;
t43 = t106 * t59;
t64 = -t58 * t42 - t55 * t43;
t63 = -t55 * t42 + t58 * t43;
t34 = t55 * t59 + t58 * t56;
t33 = t55 * t56 - t58 * t59;
t16 = t107 * t34;
t62 = -t60 * t16 + t33 * t49;
t48 = -t59 * pkin(4) - pkin(3);
t44 = -0.2e1 * t78;
t37 = t59 * t72;
t36 = t56 * t72;
t28 = t56 * t49 - t81;
t24 = (pkin(4) * t56 + t45) * t57;
t21 = t33 * t57;
t20 = t34 * t57;
t17 = t61 * pkin(4) + t75;
t15 = t107 * t33;
t12 = t15 * t60 + t34 * t49;
t11 = -t63 * qJD(5) + t55 * t36 - t58 * t37;
t10 = -t64 * qJD(5) + t58 * t36 + t55 * t37;
t9 = -t91 * t102 + (t107 * t101 + t80) * t58 + t108 * t55;
t8 = t16 * t57 + t55 * t80 - t58 * t76;
t7 = t109 * qJD(4) + t97;
t2 = -t65 * qJD(5) + t74;
t1 = -t66 * qJD(5) - t55 * t4 - t105;
t3 = [0, 0, 0, 0, 0.2e1 * t78, -0.2e1 * t70, 0, 0, 0, t57 * t88, t60 * t88, -0.2e1 * t51 * t79 + 0.2e1 * t52 * t78, t69 * t110 + 0.2e1 * t51 * t71, 0.2e1 * t57 * t82 + 0.2e1 * t59 * t70, -0.2e1 * t56 * t70 + 0.2e1 * t57 * t81, t44, 0.2e1 * t23 * t49 - 0.2e1 * t7 * t60 + 0.2e1 * (t51 * t93 + t56 * t78) * t45, -0.2e1 * t51 * t45 * t94 - 0.2e1 * t6 * t60 + 0.2e1 * (-t22 + t35) * t49, 0.2e1 * t21 * t8, 0.2e1 * t8 * t20 + 0.2e1 * t21 * t9, -0.2e1 * t21 * t49 + 0.2e1 * t8 * t60, -0.2e1 * t20 * t49 + 0.2e1 * t60 * t9, t44, 0.2e1 * t17 * t20 - 0.2e1 * t2 * t60 + 0.2e1 * t24 * t9 + 0.2e1 * t66 * t49, -0.2e1 * t1 * t60 - 0.2e1 * t17 * t21 - 0.2e1 * t24 * t8 - 0.2e1 * t65 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t90, -t49, 0, -t75, t77, -t57 * t71 + t69, t79 * t110 - t96 * t90, t28, t26, 0, (pkin(7) * t99 + (-pkin(3) * t59 + t104) * t57) * qJD(4) + (t68 * t56 - t35) * qJD(3), (t45 * t101 + t67 * t56) * qJD(4) + (t60 * t104 + t68 * t59) * qJD(3), t21 * t15 - t8 * t34, t15 * t20 + t21 * t16 + t8 * t33 - t34 * t9, t12, -t62, 0, -t11 * t60 + t24 * t16 + t17 * t33 + t20 * t87 + t48 * t9 + t64 * t49, -t10 * t60 - t24 * t15 + t17 * t34 - t21 * t87 - t48 * t8 - t63 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t90, 0, 0, 0, 0, 0, -t26, t28, 0, 0, 0, 0, 0, t62, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t79, -0.2e1 * t71, 0, 0, 0, t56 * t89, t59 * t89, -0.2e1 * t34 * t15, 0.2e1 * t15 * t33 - 0.2e1 * t34 * t16, 0, 0, 0, 0.2e1 * t48 * t16 + 0.2e1 * t33 * t87, -0.2e1 * t48 * t15 + 0.2e1 * t34 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, -t61, t49, t7, t6, 0, 0, -t8, -t9, t49, t58 * t86 + (t73 * t55 - t100) * qJD(5) + t74, -t105 + (-t4 - t86) * t55 + (t73 * t58 + t103) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t108, 0, 0, 0, 0, 0, -t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, -t94, 0, -pkin(7) * t93, pkin(7) * t94, 0, 0, -t15, -t16, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t85, -0.2e1 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, t49, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
