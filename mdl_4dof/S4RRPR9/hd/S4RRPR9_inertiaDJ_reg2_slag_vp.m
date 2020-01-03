% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:04
% EndTime: 2019-12-31 17:10:07
% DurationCPUTime: 0.72s
% Computational Cost: add. (707->117), mult. (1860->252), div. (0->0), fcn. (1552->6), ass. (0->79)
t45 = sin(pkin(7));
t46 = cos(pkin(7));
t89 = cos(qJ(4));
t66 = qJD(4) * t89;
t47 = sin(qJ(4));
t83 = qJD(4) * t47;
t97 = -t45 * t83 + t46 * t66;
t68 = t89 * t46;
t96 = -t47 * t45 + t68;
t48 = sin(qJ(2));
t49 = cos(qJ(2));
t95 = (t48 ^ 2 - t49 ^ 2) * qJD(2);
t94 = 0.2e1 * t97;
t27 = t45 * t89 + t46 * t47;
t22 = t27 * qJD(4);
t93 = 0.2e1 * t22;
t92 = pkin(2) * t48;
t91 = pkin(5) * t45;
t90 = t46 * pkin(5);
t88 = t45 * t48;
t86 = pkin(6) + qJ(3);
t59 = -pkin(2) * t49 - qJ(3) * t48;
t29 = -pkin(1) + t59;
t37 = t49 * t90;
t18 = t45 * t29 + t37;
t84 = qJD(3) * t45;
t82 = t46 * qJD(3);
t81 = t48 * qJD(2);
t80 = t48 * qJD(3);
t79 = t49 * qJD(2);
t78 = t49 * t91;
t77 = t48 * t90;
t76 = -0.2e1 * pkin(1) * qJD(2);
t40 = pkin(5) * t79;
t74 = t45 * t79;
t73 = t45 * t80;
t72 = t46 * t79;
t71 = t46 * t80;
t70 = t48 * t79;
t69 = pkin(3) + t91;
t67 = t86 * t45;
t65 = t89 * qJD(3);
t64 = 0.2e1 * t70;
t63 = t45 * t72;
t41 = t45 ^ 2;
t42 = t46 ^ 2;
t62 = 0.2e1 * (t41 + t42) * qJD(3);
t61 = -0.2e1 * t95;
t58 = -qJ(3) * t49 + t92;
t12 = -t71 + (pkin(5) * t88 + t46 * t58) * qJD(2);
t13 = -t73 + (t45 * t58 - t77) * qJD(2);
t57 = -t12 * t45 + t13 * t46;
t56 = t89 * t67;
t55 = -t49 * t86 + t92;
t54 = -t69 * t49 + (-pkin(6) * t48 + t29) * t46;
t53 = t47 * t54;
t52 = t89 * t54;
t51 = t73 - (t45 * t55 - t77) * qJD(2);
t50 = -t71 + (t46 * t55 + t48 * t69) * qJD(2);
t39 = -pkin(3) * t46 - pkin(2);
t34 = -0.2e1 * t70;
t32 = t86 * t46;
t28 = (pkin(3) * t45 + pkin(5)) * t48;
t23 = pkin(3) * t74 + t40;
t20 = t96 * t48;
t19 = t27 * t48;
t17 = t46 * t29 - t78;
t16 = t32 * t89 - t47 * t67;
t15 = -t47 * t32 - t56;
t14 = -pkin(6) * t88 + t18;
t9 = t27 * t79 + t48 * t97;
t8 = t22 * t48 + t47 * t74 - t68 * t79;
t7 = -t32 * t66 - t47 * t82 + (t83 * t86 - t65) * t45;
t6 = qJD(4) * t56 - t46 * t65 + (qJD(4) * t32 + t84) * t47;
t4 = t14 * t89 + t53;
t3 = -t47 * t14 + t52;
t2 = -qJD(4) * t53 - t14 * t66 + t47 * t51 + t50 * t89;
t1 = -qJD(4) * t52 + t14 * t83 - t47 * t50 + t51 * t89;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t61, 0, t34, 0, 0, t48 * t76, t49 * t76, 0, 0, t42 * t64, -0.4e1 * t48 * t63, 0.2e1 * t46 * t95, t41 * t64, t45 * t61, t34, -0.2e1 * t12 * t49 + 0.2e1 * (t17 + 0.2e1 * t78) * t81, 0.2e1 * t13 * t49 + 0.2e1 * (-t18 + 0.2e1 * t37) * t81, 0.2e1 * (-t12 * t46 - t13 * t45) * t48 + 0.2e1 * (-t17 * t46 - t18 * t45) * t79, 0.2e1 * pkin(5) ^ 2 * t70 + 0.2e1 * t12 * t17 + 0.2e1 * t13 * t18, -0.2e1 * t20 * t8, 0.2e1 * t19 * t8 - 0.2e1 * t20 * t9, 0.2e1 * t20 * t81 + 0.2e1 * t49 * t8, 0.2e1 * t19 * t9, -0.2e1 * t19 * t81 + 0.2e1 * t49 * t9, t34, 0.2e1 * t19 * t23 - 0.2e1 * t2 * t49 + 0.2e1 * t28 * t9 + 0.2e1 * t3 * t81, -0.2e1 * t1 * t49 + 0.2e1 * t20 * t23 - 0.2e1 * t28 * t8 - 0.2e1 * t4 * t81, 0.2e1 * t1 * t19 - 0.2e1 * t2 * t20 + 0.2e1 * t3 * t8 - 0.2e1 * t4 * t9, -0.2e1 * t1 * t4 + 0.2e1 * t2 * t3 + 0.2e1 * t23 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, -t81, 0, -t40, pkin(5) * t81, 0, 0, t63, (-t41 + t42) * t79, t45 * t81, -t63, t46 * t81, 0, t49 * t84 + (t45 * t59 - t37) * qJD(2), t49 * t82 + (t46 * t59 + t78) * qJD(2), t57, -pkin(2) * t40 + (-t17 * t45 + t18 * t46) * qJD(3) + t57 * qJ(3), t20 * t97 - t27 * t8, -t19 * t97 - t20 * t22 - t27 * t9 - t8 * t96, t27 * t81 - t49 * t97, t19 * t22 - t9 * t96, t22 * t49 + t81 * t96, 0, t15 * t81 + t22 * t28 - t23 * t96 + t39 * t9 - t49 * t7, -t16 * t81 + t23 * t27 + t28 * t97 - t39 * t8 - t49 * t6, -t1 * t96 + t15 * t8 - t16 * t9 + t19 * t6 - t2 * t27 - t20 * t7 - t22 * t4 - t3 * t97, -t1 * t16 + t15 * t2 + t23 * t39 + t3 * t7 - t4 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, qJ(3) * t62, t27 * t94, -0.2e1 * t22 * t27 + 0.2e1 * t96 * t97, 0, -t96 * t93, 0, 0, t39 * t93, t39 * t94, -0.2e1 * t15 * t97 - 0.2e1 * t16 * t22 - 0.2e1 * t27 * t7 - 0.2e1 * t6 * t96, 0.2e1 * t15 * t7 - 0.2e1 * t16 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t72, 0, t40, 0, 0, 0, 0, 0, 0, t9, -t8, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t97, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, -t9, t81, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, -t22, 0, t7, t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
