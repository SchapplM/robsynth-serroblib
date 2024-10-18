% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR12_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:20
% EndTime: 2024-09-28 18:08:21
% DurationCPUTime: 0.62s
% Computational Cost: add. (745->96), mult. (2484->169), div. (0->0), fcn. (2278->12), ass. (0->80)
t51 = sin(pkin(5));
t56 = sin(qJ(3));
t57 = sin(qJ(2));
t60 = cos(qJ(3));
t61 = cos(qJ(2));
t30 = (t56 * t57 - t60 * t61) * t51;
t31 = (t56 * t61 + t57 * t60) * t51;
t94 = qJD(2) + qJD(3);
t50 = sin(pkin(6));
t48 = t50 * pkin(10);
t47 = pkin(2) * t60 + pkin(3);
t55 = sin(qJ(4));
t59 = cos(qJ(4));
t82 = qJD(4) * t59;
t87 = t56 * t59;
t62 = (-t56 * t82 + (-t55 * t60 - t87) * qJD(3)) * pkin(2);
t83 = qJD(4) * t55;
t24 = -t47 * t83 + t62;
t52 = cos(pkin(6));
t85 = pkin(2) * qJD(3);
t77 = t60 * t85;
t78 = t56 * t85;
t79 = t55 * t56 * pkin(2);
t86 = qJD(4) * t79 + t55 * t78;
t23 = (-qJD(4) * t47 - t77) * t59 + t86;
t28 = pkin(2) * t87 + t47 * t55 + t48;
t34 = t47 * t59 + pkin(4) - t79;
t54 = sin(qJ(5));
t58 = cos(qJ(5));
t88 = t52 * t58;
t89 = t52 * t54;
t9 = t24 * t88 + t54 * t23 + (-t28 * t58 - t34 * t89) * qJD(5);
t49 = t50 ^ 2;
t90 = t49 * t58;
t93 = t24 * t90 + t9 * t52;
t92 = t24 * t54;
t91 = t49 * t54;
t84 = qJD(2) * t51;
t81 = qJD(5) * t54;
t80 = qJD(5) * t58;
t76 = pkin(3) * t83;
t75 = pkin(3) * t82;
t74 = t49 * t81;
t73 = t49 * t80;
t72 = t52 * t80;
t71 = t50 * t81;
t45 = t50 * t80;
t70 = qJD(5) * (-pkin(4) - t34);
t46 = pkin(3) * t59 + pkin(4);
t69 = qJD(5) * (-pkin(4) - t46);
t68 = (-pkin(3) - t47) * qJD(4);
t67 = qJD(5) * (-t34 - t46);
t66 = t58 * t76;
t16 = -t30 * t59 - t31 * t55;
t53 = cos(pkin(5));
t65 = -t16 * t52 - t50 * t53;
t17 = -t30 * t55 + t31 * t59;
t40 = pkin(3) * t55 + t48;
t39 = 0.2e1 * t54 * t73;
t38 = t76 * t91;
t37 = 0.2e1 * t52 * t45;
t36 = -0.2e1 * t52 * t71;
t33 = (-pkin(4) * t89 - t48 * t58) * qJD(5);
t32 = -pkin(4) * t72 + pkin(10) * t71;
t29 = 0.2e1 * (-t54 ^ 2 + t58 ^ 2) * t49 * qJD(5);
t27 = t33 * t52;
t21 = t94 * t31;
t20 = t94 * t30;
t13 = (-t40 * t58 - t46 * t89) * qJD(5) + (-t54 * t59 - t55 * t88) * qJD(4) * pkin(3);
t12 = -t46 * t72 - t58 * t75 + (qJD(5) * t40 + t52 * t76) * t54;
t11 = t13 * t52;
t10 = -t16 * t50 + t52 * t53;
t8 = t23 * t58 - t24 * t89 + t28 * t81 - t34 * t72;
t6 = -qJD(4) * t17 + t55 * t20 - t59 * t21;
t5 = t20 * t59 + t21 * t55 + t30 * t82 + t31 * t83;
t4 = t6 * t88 + t54 * t5 + (-t17 * t58 + t54 * t65) * qJD(5);
t3 = -t6 * t89 + t58 * t5 + (t17 * t54 + t58 * t65) * qJD(5);
t2 = t10 * t71 + t4 * t52 + t6 * t90;
t1 = t10 * t45 + t3 * t52 - t6 * t91;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t57 * t84, -t61 * t84, 0, -t21, t20, 0, t6, t5, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, -0.2e1 * t78, -0.2e1 * t77, 0, 0.2e1 * t24, 0.2e1 * t23, t39, t29, t37, t36, 0, -0.2e1 * t34 * t74 + 0.2e1 * t93, 0.2e1 * t8 * t52 + 0.2e1 * (-t34 * t80 - t92) * t49; 0, 0, 0, 0, 0, -t21, t20, 0, t6, t5, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, -t78, -t77, 0, t55 * t68 + t62, (t68 - t77) * t59 + t86, t39, t29, t37, t36, 0, t11 + (t54 * t67 - t66) * t49 + t93, t38 + (t12 + t8) * t52 + (t58 * t67 - t92) * t49; 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t76, -0.2e1 * t75, t39, t29, t37, t36, 0, 0.2e1 * t11 + 0.2e1 * (-t46 * t81 - t66) * t49, 0.2e1 * t12 * t52 - 0.2e1 * t46 * t73 + 0.2e1 * t38; 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, t24, t23, t39, t29, t37, t36, 0, t70 * t91 + t27 + t93, (t32 + t8) * t52 + (t58 * t70 - t92) * t49; 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t75, t39, t29, t37, t36, 0, t11 + t27 + (t54 * t69 - t66) * t49, t38 + (t12 + t32) * t52 + t69 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t29, t37, t36, 0, -0.2e1 * pkin(4) * t74 + 0.2e1 * t27, -0.2e1 * pkin(4) * t73 + 0.2e1 * t32 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t71, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t71, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t71, 0, t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
