% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:39
% EndTime: 2019-03-08 19:37:40
% DurationCPUTime: 0.61s
% Computational Cost: add. (403->109), mult. (1197->205), div. (0->0), fcn. (1120->10), ass. (0->90)
t46 = cos(qJ(4));
t48 = -pkin(4) - pkin(9);
t96 = t46 * t48;
t37 = t46 ^ 2;
t43 = sin(qJ(4));
t67 = qJD(4) * (t43 ^ 2 - t37);
t42 = sin(qJ(6));
t34 = t42 ^ 2;
t45 = cos(qJ(6));
t90 = -t45 ^ 2 + t34;
t66 = t90 * qJD(6);
t38 = sin(pkin(11));
t29 = t38 * pkin(2) + pkin(8);
t92 = pkin(5) + t29;
t26 = t92 * t46;
t76 = t26 * qJD(6);
t83 = qJD(5) * t46;
t87 = t43 * qJ(5);
t95 = (t87 - t96) * qJD(4) - t76 - t83;
t94 = 0.2e1 * qJD(5);
t93 = t46 * pkin(4);
t39 = sin(pkin(6));
t40 = cos(pkin(11));
t44 = sin(qJ(2));
t47 = cos(qJ(2));
t53 = (t38 * t47 + t40 * t44) * t39;
t13 = qJD(2) * t53;
t59 = t38 * t44 - t40 * t47;
t15 = t59 * t39;
t91 = t15 * t13;
t88 = qJ(5) * t46;
t86 = qJD(2) * t39;
t85 = qJD(4) * t42;
t32 = qJD(4) * t43;
t84 = qJD(4) * t45;
t33 = qJD(4) * t46;
t82 = qJD(6) * t42;
t81 = qJD(6) * t45;
t80 = qJD(6) * t46;
t79 = qJD(6) * t48;
t41 = cos(pkin(6));
t10 = t41 * t43 + t46 * t53;
t78 = t10 * qJD(6);
t77 = t15 * qJD(4);
t75 = qJ(5) * qJD(6);
t30 = -t40 * pkin(2) - pkin(3);
t74 = 0.2e1 * qJD(4) * t30;
t73 = t42 * t80;
t72 = t45 * t80;
t71 = t43 * t33;
t70 = t29 * t32;
t69 = t43 * t84;
t68 = t42 * t81;
t65 = pkin(4) * t32 - t43 * qJD(5);
t64 = t42 * t69;
t52 = t43 * t53;
t9 = -t41 * t46 + t52;
t63 = t15 * t45 + t9 * t42;
t62 = -t15 * t42 + t9 * t45;
t56 = t30 - t87;
t16 = t56 + t96;
t25 = t92 * t43;
t61 = t45 * t16 + t42 * t25;
t60 = t42 * t16 - t45 * t25;
t14 = t59 * t86;
t6 = -qJD(4) * t52 + (t41 * qJD(4) - t14) * t46;
t55 = t6 * t42 + t45 * t78;
t54 = -t42 * t78 + t6 * t45;
t17 = t92 * t32;
t51 = -t17 + (-t43 * t48 - t88) * qJD(6);
t50 = t83 + (-t87 - t93) * qJD(4);
t5 = qJD(4) * t10 - t14 * t43;
t49 = t5 * t43 + t6 * t46 + (-t10 * t43 + t46 * t9) * qJD(4);
t28 = 0.2e1 * t71;
t27 = t29 * t33;
t24 = t56 - t93;
t23 = -t42 * t32 + t72;
t22 = t42 * t33 + t43 * t81;
t21 = t69 + t73;
t20 = t45 * t33 - t43 * t82;
t19 = qJ(5) * t33 - t65;
t18 = pkin(5) * t33 + t27;
t12 = (pkin(9) * t43 - t88) * qJD(4) + t65;
t8 = t13 * t46 - t43 * t77;
t7 = t13 * t43 + t46 * t77;
t4 = -t61 * qJD(6) - t42 * t12 + t45 * t18;
t3 = t60 * qJD(6) - t45 * t12 - t42 * t18;
t2 = t62 * qJD(6) + t13 * t45 + t5 * t42;
t1 = -t63 * qJD(6) - t13 * t42 + t5 * t45;
t11 = [0, 0, 0, 0, -0.2e1 * t14 * t53 + 0.2e1 * t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t6 + 0.2e1 * t9 * t5 + 0.2e1 * t91, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t44 * t86, -t47 * t86 (-t13 * t40 - t14 * t38) * pkin(2), 0, 0, 0, 0, 0, -t8, t7, t49, t8, -t7, t13 * t24 - t15 * t19 + t29 * t49, 0, 0, 0, 0, 0 (-t10 * t84 + t1) * t43 + (t62 * qJD(4) + t54) * t46 (t10 * t85 - t2) * t43 + (-t63 * qJD(4) - t55) * t46; 0, 0, 0, 0, 0, t28, -0.2e1 * t67, 0, 0, 0, t43 * t74, t46 * t74, 0, -0.2e1 * t19 * t46 - 0.2e1 * t24 * t32, 0.2e1 * t19 * t43 - 0.2e1 * t24 * t33, -0.2e1 * t24 * t19, -0.2e1 * t34 * t71 + 0.2e1 * t37 * t68, -0.2e1 * t37 * t66 - 0.4e1 * t46 * t64, 0.2e1 * t42 * t67 - 0.2e1 * t43 * t72, 0.2e1 * t43 * t73 + 0.2e1 * t45 * t67, t28, 0.2e1 * (-t26 * t84 + t4) * t43 + 0.2e1 * (-t60 * qJD(4) - t17 * t45 - t42 * t76) * t46, 0.2e1 * (t26 * t85 + t3) * t43 + 0.2e1 * (-t61 * qJD(4) + t17 * t42 - t45 * t76) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t43 - t5 * t46 + (t10 * t46 + t43 * t9) * qJD(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, t5, t6, -t5 * pkin(4) + t6 * qJ(5) + t10 * qJD(5), 0, 0, 0, 0, 0, t55, t54; 0, 0, 0, 0, 0, 0, 0, t33, -t32, 0, -t27, t70, t50, t27, -t70, t50 * t29, t46 * t66 + t64, -t90 * t32 + 0.4e1 * t46 * t68, t20, -t22, 0, t51 * t42 - t95 * t45, t95 * t42 + t51 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, t32, t33, t19, 0, 0, 0, 0, 0, t22, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, qJ(5) * t94, -0.2e1 * t68, 0.2e1 * t66, 0, 0, 0, 0.2e1 * qJD(5) * t42 + 0.2e1 * t45 * t75, 0.2e1 * qJD(5) * t45 - 0.2e1 * t42 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, t27, 0, 0, 0, 0, 0, t20, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t21, t33, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t81, 0, -t42 * t79, -t45 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t11;
