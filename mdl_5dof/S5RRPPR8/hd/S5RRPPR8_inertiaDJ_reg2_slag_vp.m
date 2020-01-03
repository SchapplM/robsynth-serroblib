% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:15
% EndTime: 2019-12-31 19:39:18
% DurationCPUTime: 0.85s
% Computational Cost: add. (966->117), mult. (2079->227), div. (0->0), fcn. (1851->6), ass. (0->74)
t59 = sin(pkin(8));
t60 = cos(pkin(8));
t62 = sin(qJ(2));
t63 = cos(qJ(2));
t36 = t62 * t59 + t63 * t60;
t72 = t36 * qJD(2);
t90 = pkin(6) - qJ(4);
t46 = t90 * t63;
t98 = t90 * t62;
t23 = t60 * t46 + t59 * t98;
t77 = t59 * t46;
t22 = t60 * t98 - t77;
t97 = -t63 * pkin(2) - t62 * qJ(3);
t92 = t62 * t60;
t75 = t59 * t63 - t92;
t96 = t75 * qJD(4);
t95 = 2 * qJD(3);
t64 = -pkin(2) - pkin(3);
t94 = cos(qJ(5));
t55 = t63 * qJD(2);
t89 = qJ(3) * t55 + qJD(3) * t62;
t61 = sin(qJ(5));
t88 = qJD(5) * t61;
t87 = t59 * qJD(3);
t86 = t62 * qJD(2);
t85 = -0.2e1 * pkin(1) * qJD(2);
t43 = -pkin(1) + t97;
t84 = -pkin(7) + t90;
t83 = pkin(6) * t86;
t82 = pkin(6) * t55;
t81 = t62 * t55;
t80 = t94 * t60;
t78 = qJD(5) * t94;
t34 = t63 * pkin(3) - t43;
t51 = t59 * t55;
t76 = t60 * t86 - t51;
t41 = -qJ(3) * t59 + t60 * t64;
t74 = -pkin(4) + t41;
t18 = -t36 * t61 - t75 * t94;
t38 = t59 * t94 + t60 * t61;
t73 = t36 * qJD(4);
t71 = pkin(7) * t75 + t22;
t70 = t61 * t71;
t69 = t94 * t74;
t68 = t94 * t71;
t67 = qJD(2) * t97 + t63 * qJD(3);
t42 = qJ(3) * t60 + t59 * t64;
t20 = t42 * t94 + t61 * t74;
t66 = t72 * t84 + t96;
t65 = -t51 * pkin(7) - t73 + (-t84 * t92 + t77) * qJD(2);
t48 = -0.2e1 * t81;
t47 = 0.2e1 * t81;
t44 = (-t62 ^ 2 + t63 ^ 2) * qJD(2);
t35 = -t59 * t61 + t80;
t33 = t38 * qJD(5);
t32 = t59 * t88 - t60 * t78;
t30 = pkin(2) * t86 - t89;
t24 = t64 * t86 + t89;
t21 = pkin(4) * t36 + t34;
t19 = -t42 * t61 + t69;
t17 = t36 * t94 - t61 * t75;
t16 = t51 * pkin(4) + (-pkin(4) * t60 + t64) * t86 + t89;
t15 = -pkin(7) * t36 + t23;
t14 = -qJD(2) * t22 - t73;
t13 = -qJD(2) * t23 - t96;
t11 = qJD(3) * t38 + qJD(5) * t20;
t10 = -qJD(5) * t69 - qJD(3) * t80 + (qJD(5) * t42 + t87) * t61;
t6 = qJD(5) * t18 + t61 * t72 - t76 * t94;
t5 = t36 * t78 - t61 * t76 - t72 * t94 - t75 * t88;
t4 = t15 * t94 + t70;
t3 = -t15 * t61 + t68;
t2 = qJD(5) * t70 + t15 * t78 + t61 * t65 - t66 * t94;
t1 = -qJD(5) * t68 + t15 * t88 - t61 * t66 - t65 * t94;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0.2e1 * t44, 0, t48, 0, 0, t62 * t85, t63 * t85, 0, 0, t47, 0, -0.2e1 * t44, 0, 0, t48, -0.2e1 * t30 * t63 + 0.2e1 * t43 * t86, 0, -0.2e1 * t30 * t62 - 0.2e1 * t43 * t55, 0.2e1 * t43 * t30, -0.2e1 * t75 * t72, 0.2e1 * t75 * t51 + 0.2e1 * (-t36 ^ 2 - t75 * t92) * qJD(2), 0, -0.2e1 * t36 * t76, 0, 0, 0.2e1 * t24 * t36 - 0.2e1 * t34 * t76, -0.2e1 * t24 * t75 + 0.2e1 * t34 * t72, -0.2e1 * t13 * t75 - 0.2e1 * t14 * t36 - 0.2e1 * t23 * t51 + 0.2e1 * (-t22 * t36 + t23 * t92) * qJD(2), -0.2e1 * t13 * t22 + 0.2e1 * t14 * t23 + 0.2e1 * t24 * t34, -0.2e1 * t18 * t5, 0.2e1 * t17 * t5 - 0.2e1 * t18 * t6, 0, 0.2e1 * t17 * t6, 0, 0, 0.2e1 * t16 * t17 + 0.2e1 * t21 * t6, 0.2e1 * t16 * t18 - 0.2e1 * t21 * t5, 0.2e1 * t1 * t17 + 0.2e1 * t18 * t2 + 0.2e1 * t3 * t5 - 0.2e1 * t4 * t6, -0.2e1 * t1 * t4 + 0.2e1 * t16 * t21 - 0.2e1 * t2 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, -t86, 0, -t82, t83, 0, 0, 0, t55, 0, 0, t86, 0, -t82, t67, -t83, t67 * pkin(6), 0, 0, -t72, 0, -t76, 0, t13, t14, -t42 * t51 + (-t36 * t60 - t59 * t75) * qJD(3) + (-t36 * t41 + t42 * t92) * qJD(2), -t13 * t41 + t14 * t42 + (-t22 * t59 + t23 * t60) * qJD(3), 0, 0, t5, 0, t6, 0, t2, -t1, t10 * t17 + t11 * t18 + t19 * t5 - t20 * t6, -t1 * t20 - t10 * t4 - t11 * t3 - t19 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, qJ(3) * t95, 0, 0, 0, 0, 0, 0, 0.2e1 * t87, t60 * t95, 0, (-t41 * t59 + t42 * t60) * t95, 0, 0, 0, 0, 0, 0, 0.2e1 * t11, -0.2e1 * t10, 0, -0.2e1 * t10 * t20 - 0.2e1 * t11 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, -t55 * t60 ^ 2 - t59 * t51, -t13 * t60 + t14 * t59, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t32 + t18 * t33 + t35 * t5 - t38 * t6, -t1 * t38 - t2 * t35 - t3 * t33 - t32 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, 0, -t10 * t38 - t11 * t35 - t19 * t33 - t20 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t32 * t38 - 0.2e1 * t33 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t72, 0, t24, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t6, 0, -t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
