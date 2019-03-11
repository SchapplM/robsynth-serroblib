% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:06:36
% EndTime: 2019-03-08 21:06:38
% DurationCPUTime: 0.69s
% Computational Cost: add. (858->134), mult. (2255->260), div. (0->0), fcn. (2144->10), ass. (0->92)
t52 = sin(qJ(6));
t46 = t52 ^ 2;
t55 = cos(qJ(6));
t92 = -t55 ^ 2 + t46;
t79 = t92 * qJD(6);
t101 = 2 * qJD(5);
t100 = pkin(4) + pkin(9);
t48 = sin(pkin(11));
t50 = cos(pkin(11));
t53 = sin(qJ(3));
t56 = cos(qJ(3));
t32 = t48 * t53 - t50 * t56;
t99 = t32 * t52;
t98 = t32 * t55;
t49 = sin(pkin(6));
t54 = sin(qJ(2));
t97 = t49 * t54;
t57 = cos(qJ(2));
t96 = t49 * t57;
t33 = t48 * t56 + t50 * t53;
t28 = t33 * qJD(3);
t95 = t52 * t28;
t94 = t55 * t28;
t93 = -qJ(4) - pkin(8);
t91 = qJD(2) * t54;
t90 = qJD(2) * t57;
t89 = qJD(6) * t52;
t88 = qJD(6) * t55;
t87 = t53 * qJD(3);
t86 = t56 * qJD(3);
t85 = -0.2e1 * pkin(2) * qJD(3);
t84 = t52 * t94;
t44 = pkin(3) * t87;
t83 = t57 * t87;
t38 = t49 * t91;
t82 = t49 * t90;
t81 = t52 * t88;
t43 = -t56 * pkin(3) - pkin(2);
t42 = -pkin(3) * t50 - pkin(4);
t80 = qJD(3) * t93;
t27 = t56 * qJD(4) + t53 * t80;
t65 = -t53 * qJD(4) + t56 * t80;
t12 = t27 * t48 - t50 * t65;
t35 = t93 * t53;
t36 = t93 * t56;
t23 = -t35 * t50 - t36 * t48;
t74 = -t33 * qJ(5) + t43;
t11 = t100 * t32 + t74;
t14 = pkin(5) * t33 + t23;
t78 = t11 * t55 + t14 * t52;
t77 = t11 * t52 - t14 * t55;
t13 = t50 * t27 + t48 * t65;
t24 = t35 * t48 - t36 * t50;
t76 = t12 * t23 + t13 * t24;
t40 = pkin(3) * t48 + qJ(5);
t75 = -qJD(5) * t32 - t28 * t40;
t51 = cos(pkin(6));
t30 = t51 * t53 + t56 * t97;
t71 = t51 * t56 - t53 * t97;
t16 = t30 * t48 - t50 * t71;
t73 = -t16 * t52 + t55 * t96;
t72 = t16 * t55 + t52 * t96;
t70 = t32 * t88 + t95;
t69 = t32 * t89 - t94;
t29 = -t48 * t87 + t50 * t86;
t68 = t29 * t52 + t33 * t88;
t67 = -t29 * t55 + t33 * t89;
t66 = -t29 * qJ(5) - t33 * qJD(5) + t44;
t39 = -pkin(9) + t42;
t9 = -t28 * pkin(5) + t13;
t64 = t9 + (t32 * t40 - t33 * t39) * qJD(6);
t17 = t30 * t50 + t48 * t71;
t22 = qJD(3) * t71 + t56 * t82;
t58 = -qJD(3) * t30 - t53 * t82;
t6 = t22 * t48 - t50 * t58;
t7 = t50 * t22 + t48 * t58;
t63 = t16 * t12 + t17 * t13 + t6 * t23 + t7 * t24;
t62 = t16 * t29 - t17 * t28 - t7 * t32 + t6 * t33;
t15 = -pkin(5) * t32 + t24;
t61 = -qJD(6) * t15 - t29 * t39 - t75;
t60 = -0.2e1 * t49 ^ 2 * t54 * t90 + 0.2e1 * t16 * t6 + 0.2e1 * t17 * t7;
t59 = 0.2e1 * t12 * t33 - 0.2e1 * t13 * t32 + 0.2e1 * t23 * t29 - 0.2e1 * t24 * t28;
t31 = t32 ^ 2;
t20 = pkin(4) * t32 + t74;
t10 = pkin(4) * t28 + t66;
t8 = pkin(5) * t29 + t12;
t5 = t100 * t28 + t66;
t4 = qJD(6) * t72 + t38 * t55 + t52 * t6;
t3 = qJD(6) * t73 - t38 * t52 + t55 * t6;
t2 = -qJD(6) * t78 - t52 * t5 + t55 * t8;
t1 = qJD(6) * t77 - t55 * t5 - t52 * t8;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t38, -t82, 0, 0, 0, 0, 0 (-t56 * t91 - t83) * t49 (t53 * t91 - t57 * t86) * t49, t62 (-pkin(3) * t83 + t43 * t91) * t49 + t63, t62 (t28 * t57 - t32 * t91) * t49 (t29 * t57 - t33 * t91) * t49 (-t10 * t57 + t20 * t91) * t49 + t63, 0, 0, 0, 0, 0, t17 * t69 + t29 * t72 + t3 * t33 - t7 * t98, t17 * t70 + t29 * t73 - t4 * t33 + t7 * t99; 0, 0, 0, 0, 0.2e1 * t53 * t86, 0.2e1 * (-t53 ^ 2 + t56 ^ 2) * qJD(3), 0, 0, 0, t53 * t85, t56 * t85, t59, 0.2e1 * t43 * t44 + 0.2e1 * t76, t59, -0.2e1 * t10 * t32 - 0.2e1 * t20 * t28, -0.2e1 * t10 * t33 - 0.2e1 * t20 * t29, 0.2e1 * t10 * t20 + 0.2e1 * t76, 0.2e1 * t28 * t32 * t46 + 0.2e1 * t31 * t81, -0.2e1 * t31 * t79 + 0.4e1 * t32 * t84, 0.2e1 * t32 * t68 + 0.2e1 * t33 * t95, -0.2e1 * t32 * t67 + 0.2e1 * t33 * t94, 0.2e1 * t33 * t29, 0.2e1 * t15 * t69 + 0.2e1 * t2 * t33 - 0.2e1 * t29 * t77 - 0.2e1 * t9 * t98, 0.2e1 * t1 * t33 + 0.2e1 * t15 * t70 - 0.2e1 * t29 * t78 + 0.2e1 * t9 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t22, 0 (t48 * t7 - t50 * t6) * pkin(3), 0, t6, t7, qJD(5) * t17 + t40 * t7 + t42 * t6, 0, 0, 0, 0, 0, t17 * t88 + t52 * t7, -t17 * t89 + t55 * t7; 0, 0, 0, 0, 0, 0, t86, -t87, 0, -pkin(8) * t86, pkin(8) * t87 (-t28 * t48 - t29 * t50) * pkin(3) (-t12 * t50 + t13 * t48) * pkin(3), t29 * t42 + t75, t12, t13, qJD(5) * t24 + t12 * t42 + t13 * t40, -t32 * t79 + t84, -t28 * t92 - 0.4e1 * t32 * t81, -t67, -t68, 0, t52 * t64 - t55 * t61, t52 * t61 + t55 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t40 * t101, -0.2e1 * t81, 0.2e1 * t79, 0, 0, 0, 0.2e1 * qJD(5) * t52 + 0.2e1 * t40 * t88, 0.2e1 * qJD(5) * t55 - 0.2e1 * t40 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t28, -t29, t10, 0, 0, 0, 0, 0, -t68, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, t12, 0, 0, 0, 0, 0, -t67, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t69, t29, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t88, 0, -t39 * t89, -t39 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t18;
