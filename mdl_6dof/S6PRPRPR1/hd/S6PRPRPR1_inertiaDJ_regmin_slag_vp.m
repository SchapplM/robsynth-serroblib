% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x21]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:27:57
% EndTime: 2019-03-08 19:27:59
% DurationCPUTime: 0.64s
% Computational Cost: add. (840->125), mult. (2261->248), div. (0->0), fcn. (2290->12), ass. (0->92)
t102 = 2 * qJD(6);
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t48 = sin(pkin(11));
t78 = t48 * pkin(2) + pkin(8);
t74 = qJ(5) + t78;
t62 = qJD(4) * t74;
t26 = t55 * qJD(5) - t52 * t62;
t47 = sin(pkin(12));
t59 = -t52 * qJD(5) - t55 * t62;
t89 = cos(pkin(12));
t11 = t47 * t26 - t89 * t59;
t51 = sin(qJ(6));
t101 = t11 * t51;
t49 = sin(pkin(6));
t50 = cos(pkin(11));
t53 = sin(qJ(2));
t56 = cos(qJ(2));
t60 = (t48 * t56 + t50 * t53) * t49;
t27 = qJD(2) * t60;
t66 = t48 * t53 - t50 * t56;
t30 = t66 * t49;
t100 = t30 * t27;
t37 = t47 * t55 + t89 * t52;
t32 = t37 * qJD(4);
t76 = t89 * t55;
t36 = t47 * t52 - t76;
t99 = t36 * t32;
t85 = t52 * qJD(4);
t33 = qJD(4) * t76 - t47 * t85;
t98 = t37 * t33;
t54 = cos(qJ(6));
t97 = t37 * t54;
t96 = t51 * t32;
t95 = t51 * t33;
t94 = t54 * t32;
t93 = t54 * t33;
t92 = t36 * t93 + t37 * t94;
t46 = t54 ^ 2;
t91 = t51 ^ 2 - t46;
t90 = cos(pkin(6));
t88 = qJD(2) * t49;
t87 = qJD(6) * t51;
t86 = qJD(6) * t54;
t84 = t55 * qJD(4);
t42 = -t89 * pkin(4) - pkin(5);
t83 = t42 * t102;
t82 = 0.2e1 * t84;
t44 = pkin(4) * t85;
t81 = t37 * t87;
t80 = t30 * t85;
t79 = t51 * t86;
t43 = -t50 * pkin(2) - pkin(3);
t77 = -0.4e1 * t51 * t97;
t75 = t91 * qJD(6);
t24 = t90 * t52 + t55 * t60;
t58 = -t52 * t60 + t90 * t55;
t8 = t89 * t24 + t47 * t58;
t73 = t30 * t54 - t51 * t8;
t72 = t30 * t51 + t54 * t8;
t71 = t78 * qJD(4);
t38 = -t55 * pkin(4) + t43;
t14 = t36 * pkin(5) - t37 * pkin(9) + t38;
t34 = t74 * t55;
t65 = t74 * t52;
t20 = t89 * t34 - t47 * t65;
t70 = t54 * t14 - t51 * t20;
t69 = t51 * t14 + t54 * t20;
t41 = t47 * pkin(4) + pkin(9);
t68 = -t32 * t41 + t33 * t42;
t67 = t36 * t41 - t37 * t42;
t28 = t66 * t88;
t10 = t58 * qJD(4) - t28 * t55;
t57 = -t24 * qJD(4) + t28 * t52;
t5 = t47 * t10 - t89 * t57;
t7 = t47 * t24 - t89 * t58;
t64 = t5 * t51 + t7 * t86;
t63 = -t5 * t54 + t7 * t87;
t17 = t36 * t86 + t96;
t61 = t37 * t86 + t95;
t16 = t81 - t93;
t35 = t37 ^ 2;
t19 = t47 * t34 + t89 * t65;
t15 = t36 * t87 - t94;
t13 = t32 * pkin(5) - t33 * pkin(9) + t44;
t12 = t89 * t26 + t47 * t59;
t6 = t89 * t10 + t47 * t57;
t4 = -t69 * qJD(6) - t51 * t12 + t54 * t13;
t3 = -t70 * qJD(6) - t54 * t12 - t51 * t13;
t2 = -t72 * qJD(6) + t27 * t54 - t51 * t6;
t1 = -t73 * qJD(6) - t27 * t51 - t54 * t6;
t9 = [0, 0, 0, 0, -0.2e1 * t28 * t60 + 0.2e1 * t100, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t7 * t5 + 0.2e1 * t8 * t6 + 0.2e1 * t100, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t53 * t88, -t56 * t88 (-t27 * t50 - t28 * t48) * pkin(2), 0, 0, 0, 0, 0, -t27 * t55 + t80, t27 * t52 + t30 * t84, -t8 * t32 + t7 * t33 - t6 * t36 + t5 * t37, pkin(4) * t80 + t7 * t11 + t8 * t12 + t5 * t19 + t6 * t20 + t27 * t38, 0, 0, 0, 0, 0, t2 * t36 + t73 * t32 + t64 * t37 + t7 * t95, t1 * t36 - t72 * t32 - t63 * t37 + t7 * t93; 0, 0, 0, 0, 0, t52 * t82, 0.2e1 * (-t52 ^ 2 + t55 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * t43 * t85, t43 * t82, 0.2e1 * t11 * t37 - 0.2e1 * t12 * t36 + 0.2e1 * t19 * t33 - 0.2e1 * t20 * t32, 0.2e1 * t19 * t11 + 0.2e1 * t20 * t12 + 0.2e1 * t38 * t44, -0.2e1 * t35 * t79 + 0.2e1 * t46 * t98, t91 * t35 * t102 + t33 * t77, -0.2e1 * t36 * t81 + 0.2e1 * t92, -0.2e1 * t36 * t61 - 0.2e1 * t37 * t96, 0.2e1 * t99, 0.2e1 * t37 * t101 + 0.2e1 * t61 * t19 + 0.2e1 * t70 * t32 + 0.2e1 * t4 * t36, 0.2e1 * t11 * t97 - 0.2e1 * t16 * t19 + 0.2e1 * t3 * t36 - 0.2e1 * t69 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t32 + t8 * t33 + t5 * t36 + t6 * t37, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t36 + t12 * t37 + t19 * t32 + t20 * t33, 0, 0, 0, 0, 0, 0 (-t37 * t32 - t33 * t36) * t54 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t98 + 0.2e1 * t99, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t10, 0 (t47 * t6 - t89 * t5) * pkin(4), 0, 0, 0, 0, 0, t63, t64; 0, 0, 0, 0, 0, 0, 0, t84, -t85, 0, -t55 * t71, t52 * t71 (-t32 * t47 - t89 * t33) * pkin(4) (-t89 * t11 + t12 * t47) * pkin(4), -t37 * t75 + t51 * t93, qJD(6) * t77 - t91 * t33, t17, -t15, 0, -t11 * t54 + t68 * t51 + (t19 * t51 - t67 * t54) * qJD(6), t101 + t68 * t54 + (t19 * t54 + t67 * t51) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t84, 0 (-t89 * t32 + t33 * t47) * pkin(4), 0, 0, 0, 0, 0, t15, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t79, -0.2e1 * t75, 0, 0, 0, t51 * t83, t54 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, -t15, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t61, t32, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t87, 0, -t41 * t86, t41 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
