% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR13_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:27
% EndTime: 2019-12-31 19:15:31
% DurationCPUTime: 0.92s
% Computational Cost: add. (667->133), mult. (1690->258), div. (0->0), fcn. (1415->6), ass. (0->102)
t111 = qJD(4) + qJD(5);
t53 = sin(qJ(5));
t54 = sin(qJ(4));
t56 = cos(qJ(5));
t57 = cos(qJ(4));
t32 = t53 * t54 - t56 * t57;
t15 = t111 * t32;
t55 = sin(qJ(3));
t114 = t15 * t55;
t33 = t53 * t57 + t56 * t54;
t58 = cos(qJ(3));
t113 = t33 * t58;
t66 = t55 * pkin(3) - t58 * pkin(7);
t36 = qJ(2) + t66;
t28 = t54 * t36;
t103 = t55 * t57;
t59 = -pkin(1) - pkin(6);
t44 = t59 * t103;
t112 = -t28 - t44;
t91 = t55 * qJD(3);
t78 = t57 * t91;
t94 = qJD(4) * t58;
t84 = t54 * t94;
t24 = -t78 - t84;
t50 = t55 ^ 2;
t52 = t58 ^ 2;
t70 = (t50 - t52) * qJD(3);
t51 = t57 ^ 2;
t99 = t54 ^ 2 - t51;
t71 = t99 * qJD(4);
t110 = 2 * qJD(2);
t109 = pkin(7) + pkin(8);
t67 = pkin(3) * t58 + pkin(7) * t55;
t30 = t67 * qJD(3) + qJD(2);
t48 = t58 * qJD(3);
t76 = t59 * t48;
t93 = qJD(4) * t59;
t83 = t54 * t93;
t95 = qJD(4) * t57;
t12 = -t54 * t30 - t36 * t95 + t55 * t83 - t57 * t76;
t81 = t54 * t91;
t82 = t57 * t94;
t26 = t81 - t82;
t5 = t26 * pkin(8) - t12;
t108 = t56 * t5;
t16 = t111 * t33;
t107 = t16 * t55;
t105 = t54 * t58;
t17 = -pkin(8) * t105 - t112;
t106 = t53 * t17;
t104 = t54 * t59;
t102 = t56 * t17;
t101 = t57 * t58;
t100 = t58 * t59;
t97 = t50 + t52;
t96 = qJD(4) * t54;
t92 = qJD(5) * t53;
t90 = qJ(2) * qJD(3);
t89 = -0.2e1 * pkin(3) * qJD(4);
t88 = pkin(4) * t96;
t87 = pkin(4) * t48;
t86 = pkin(4) * t92;
t85 = qJD(5) * t56 * pkin(4);
t80 = t54 * t95;
t79 = t59 * t91;
t77 = t55 * t48;
t22 = t57 * t30;
t73 = pkin(4) - t104;
t4 = t22 + (-t44 + (pkin(8) * t58 - t36) * t54) * qJD(4) + (pkin(8) * t103 + t73 * t58) * qJD(3);
t75 = t56 * t4 - t53 * t5;
t29 = t57 * t36;
t14 = -pkin(8) * t101 + t73 * t55 + t29;
t74 = -t55 * pkin(4) - t14;
t72 = qJD(4) * t109;
t69 = t54 * t76;
t68 = t54 * t78;
t65 = t56 * t14 - t106;
t64 = t53 * t14 + t102;
t42 = t109 * t54;
t43 = t109 * t57;
t63 = -t56 * t42 - t53 * t43;
t62 = -t53 * t42 + t56 * t43;
t20 = t32 * t58;
t60 = qJD(3) * t32;
t47 = -t57 * pkin(4) - pkin(3);
t45 = 0.2e1 * t77;
t35 = t57 * t72;
t34 = t54 * t72;
t31 = (pkin(4) * t54 - t59) * t58;
t25 = t54 * t48 + t55 * t95;
t23 = -t57 * t48 + t55 * t96;
t18 = -t26 * pkin(4) + t79;
t13 = t112 * qJD(4) + t22 - t69;
t11 = -t62 * qJD(5) + t53 * t34 - t56 * t35;
t10 = -t63 * qJD(5) + t56 * t34 + t53 * t35;
t9 = -t92 * t105 + (t111 * t101 - t81) * t56 + t24 * t53;
t8 = -qJD(3) * t113 + t114;
t7 = -t111 * t113 + t55 * t60;
t6 = t58 * t60 + t107;
t2 = -t64 * qJD(5) + t75;
t1 = -t65 * qJD(5) - t53 * t4 - t108;
t3 = [0, 0, 0, 0, t110, qJ(2) * t110, -0.2e1 * t77, 0.2e1 * t70, 0, 0, 0, 0.2e1 * qJD(2) * t55 + 0.2e1 * t58 * t90, 0.2e1 * qJD(2) * t58 - 0.2e1 * t55 * t90, -0.2e1 * t51 * t77 - 0.2e1 * t52 * t80, 0.2e1 * t52 * t71 + 0.4e1 * t58 * t68, -0.2e1 * t55 * t84 - 0.2e1 * t57 * t70, 0.2e1 * t54 * t70 - 0.2e1 * t55 * t82, t45, -0.2e1 * t52 * t57 * t93 + 0.2e1 * t29 * t48 + 0.2e1 * (t13 + t69) * t55, 0.2e1 * t52 * t83 + 0.2e1 * t12 * t55 + 0.2e1 * (-t28 + t44) * t48, -0.2e1 * t20 * t7, -0.2e1 * t113 * t7 + 0.2e1 * t20 * t9, -0.2e1 * t20 * t48 + 0.2e1 * t7 * t55, -0.2e1 * t113 * t48 - 0.2e1 * t9 * t55, t45, 0.2e1 * t113 * t18 + 0.2e1 * t2 * t55 + 0.2e1 * t31 * t9 + 0.2e1 * t65 * t48, 0.2e1 * t1 * t55 - 0.2e1 * t18 * t20 + 0.2e1 * t31 * t7 - 0.2e1 * t64 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97 * t95, t97 * t96, 0, 0, 0, 0, 0, t8 * t55 - t58 * t9, t6 * t55 - t58 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t48, 0, -t79, -t76, -t58 * t71 - t68, -0.4e1 * t58 * t80 + t99 * t91, t25, -t23, 0, (-t54 * t100 - t67 * t57) * qJD(4) + (t66 * t54 - t44) * qJD(3), (-t57 * t100 + t67 * t54) * qJD(4) + (-pkin(7) * t101 + (pkin(3) * t57 + t104) * t55) * qJD(3), t20 * t15 + t7 * t33, t113 * t15 + t20 * t16 - t7 * t32 - t33 * t9, t33 * t48 - t114, -t32 * t48 - t107, 0, t11 * t55 + t113 * t88 + t31 * t16 + t18 * t32 + t47 * t9 + t63 * t48, t10 * t55 - t31 * t15 + t18 * t33 - t20 * t88 + t47 * t7 - t62 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t48, 0, 0, 0, 0, 0, t24, t26, 0, 0, 0, 0, 0, -t58 * t16 + t32 * t91, t58 * t15 + t33 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t80, -0.2e1 * t71, 0, 0, 0, t54 * t89, t57 * t89, -0.2e1 * t33 * t15, 0.2e1 * t15 * t32 - 0.2e1 * t33 * t16, 0, 0, 0, 0.2e1 * t47 * t16 + 0.2e1 * t32 * t88, -0.2e1 * t47 * t15 + 0.2e1 * t33 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t26, t48, t13, t12, 0, 0, t7, -t9, t48, t56 * t87 + (t74 * t53 - t102) * qJD(5) + t75, -t108 + (-t4 - t87) * t53 + (t74 * t56 + t106) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t23, 0, 0, 0, 0, 0, t8, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t96, 0, -pkin(7) * t95, pkin(7) * t96, 0, 0, -t15, -t16, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t86, -0.2e1 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t9, t48, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
