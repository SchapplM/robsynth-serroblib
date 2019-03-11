% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:49:34
% EndTime: 2019-03-08 19:49:37
% DurationCPUTime: 0.94s
% Computational Cost: add. (715->159), mult. (1930->306), div. (0->0), fcn. (1792->10), ass. (0->99)
t55 = sin(pkin(11));
t57 = cos(pkin(11));
t59 = sin(qJ(6));
t62 = cos(qJ(6));
t38 = t59 * t55 - t62 * t57;
t31 = t38 * qJD(6);
t60 = sin(qJ(4));
t113 = t31 * t60;
t63 = cos(qJ(4));
t28 = t38 * t63;
t39 = t62 * t55 + t59 * t57;
t112 = t39 * t63;
t52 = -t57 * pkin(5) - pkin(4);
t111 = 0.2e1 * t52;
t110 = 2 * qJD(3);
t109 = t60 * pkin(4);
t56 = sin(pkin(6));
t64 = cos(qJ(2));
t103 = t56 * t64;
t58 = cos(pkin(6));
t34 = -t60 * t103 + t58 * t63;
t61 = sin(qJ(2));
t104 = t56 * t61;
t49 = qJD(2) * t104;
t22 = t34 * qJD(4) - t63 * t49;
t108 = t22 * t63;
t32 = t39 * qJD(6);
t107 = t32 * t60;
t33 = t63 * t103 + t58 * t60;
t106 = t33 * t60;
t105 = t55 * t63;
t102 = t57 * t60;
t101 = t57 * t63;
t65 = -pkin(2) - pkin(8);
t99 = t60 * t65;
t97 = pkin(9) + qJ(5);
t29 = -t63 * qJD(5) + qJD(3) + (pkin(4) * t63 + qJ(5) * t60) * qJD(4);
t91 = t63 * qJD(4);
t84 = t65 * t91;
t17 = t55 * t29 + t57 * t84;
t77 = -t63 * qJ(5) + t109;
t40 = qJ(3) + t77;
t48 = t57 * t99;
t24 = t55 * t40 + t48;
t96 = t55 ^ 2 + t57 ^ 2;
t95 = qJD(2) * t64;
t94 = qJD(5) * t60;
t92 = t60 * qJD(4);
t90 = qJ(3) * qJD(4);
t89 = t55 * t99;
t88 = t55 * t92;
t87 = t56 * t95;
t86 = t57 * t92;
t50 = t65 * t92;
t85 = t60 * t91;
t83 = -t55 * t65 + pkin(5);
t82 = t96 * t63;
t81 = t96 * qJD(5);
t80 = 0.2e1 * t85;
t79 = 0.2e1 * t81;
t21 = t33 * qJD(4) - t60 * t49;
t8 = t21 * t55 + t57 * t87;
t9 = -t21 * t57 + t55 * t87;
t78 = -t8 * t55 + t9 * t57;
t36 = t57 * t40;
t15 = -pkin(9) * t101 + t83 * t60 + t36;
t18 = -pkin(9) * t105 + t24;
t76 = t62 * t15 - t59 * t18;
t75 = t59 * t15 + t62 * t18;
t26 = t57 * t29;
t16 = -t55 * t84 + t26;
t74 = -t16 * t55 + t17 * t57;
t19 = t57 * t104 - t34 * t55;
t20 = t55 * t104 + t34 * t57;
t73 = -t19 * t55 + t20 * t57;
t72 = t62 * t19 - t59 * t20;
t71 = t59 * t19 + t62 * t20;
t23 = t36 - t89;
t70 = -t23 * t55 + t24 * t57;
t45 = t97 * t55;
t46 = t97 * t57;
t69 = -t62 * t45 - t59 * t46;
t68 = -t59 * t45 + t62 * t46;
t66 = qJD(4) * t38;
t37 = (pkin(5) * t55 - t65) * t63;
t30 = -pkin(5) * t88 + t50;
t14 = pkin(9) * t88 + t17;
t13 = -qJD(6) * t28 - t59 * t86 - t62 * t88;
t12 = -qJD(4) * t112 + t113;
t11 = -qJD(6) * t112 + t60 * t66;
t10 = t63 * t66 + t107;
t7 = -t39 * qJD(5) - t68 * qJD(6);
t6 = t38 * qJD(5) - t69 * qJD(6);
t5 = t26 + (pkin(9) * t102 + t83 * t63) * qJD(4);
t4 = -t71 * qJD(6) - t59 * t9 + t62 * t8;
t3 = -t72 * qJD(6) - t59 * t8 - t62 * t9;
t2 = -t75 * qJD(6) - t59 * t14 + t62 * t5;
t1 = -t76 * qJD(6) - t62 * t14 - t59 * t5;
t25 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t19 * t8 + 0.2e1 * t20 * t9 + 0.2e1 * t33 * t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t49, -t87, t49, t87 (qJD(3) * t61 + (-pkin(2) * t61 + qJ(3) * t64) * qJD(2)) * t56, 0, 0, 0, 0, 0 (t60 * t95 + t61 * t91) * t56 (-t61 * t92 + t63 * t95) * t56, t22 * t105 + t8 * t60 + (-t55 * t106 + t19 * t63) * qJD(4), t22 * t101 - t9 * t60 + (-t33 * t102 - t20 * t63) * qJD(4) (-t55 * t9 - t57 * t8) * t63 + (t19 * t57 + t20 * t55) * t92, t19 * t16 + t20 * t17 + t8 * t23 + t9 * t24 + (t33 * t92 - t108) * t65, 0, 0, 0, 0, 0, t112 * t22 + t33 * t13 + t4 * t60 + t72 * t91, t33 * t11 - t22 * t28 + t3 * t60 - t71 * t91; 0, 0, 0, 0, 0, t110, qJ(3) * t110, -0.2e1 * t85, 0.2e1 * (t60 ^ 2 - t63 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t60 + 0.2e1 * t63 * t90, 0.2e1 * qJD(3) * t63 - 0.2e1 * t60 * t90, 0.2e1 * t16 * t60 + 0.2e1 * (t23 + 0.2e1 * t89) * t91, -0.2e1 * t17 * t60 + 0.2e1 * (-t24 + 0.2e1 * t48) * t91, 0.2e1 * (-t16 * t57 - t17 * t55) * t63 + 0.2e1 * (t23 * t57 + t24 * t55) * t92, -0.2e1 * t65 ^ 2 * t85 + 0.2e1 * t23 * t16 + 0.2e1 * t24 * t17, -0.2e1 * t28 * t11, -0.2e1 * t11 * t112 + 0.2e1 * t28 * t13, 0.2e1 * t11 * t60 - 0.2e1 * t28 * t91, -0.2e1 * t112 * t91 - 0.2e1 * t13 * t60, t80, 0.2e1 * t112 * t30 + 0.2e1 * t37 * t13 + 0.2e1 * t2 * t60 + 0.2e1 * t76 * t91, 0.2e1 * t1 * t60 + 0.2e1 * t37 * t11 - 0.2e1 * t30 * t28 - 0.2e1 * t75 * t91; 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108 + t78 * t60 + (t73 * t63 + t106) * qJD(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t60 + (t70 - 0.2e1 * t99) * t91, 0, 0, 0, 0, 0, t12 * t60 - t63 * t13, t10 * t60 - t63 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t96) * t80, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t21, -t22 * t57, t22 * t55, t78, -t22 * pkin(4) + t78 * qJ(5) + t73 * qJD(5), 0, 0, 0, 0, 0, t22 * t38 + t33 * t32, t22 * t39 - t33 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t91, 0, -t50, -t84, -t55 * t94 + (t77 * t55 - t48) * qJD(4), -t57 * t94 + (t77 * t57 + t89) * qJD(4), t74, -pkin(4) * t50 + t74 * qJ(5) + t70 * qJD(5), t11 * t39 + t28 * t31, -t11 * t38 + t112 * t31 - t39 * t13 + t28 * t32, t39 * t91 - t113, -t38 * t91 - t107, 0, t52 * t13 + t30 * t38 + t37 * t32 + t7 * t60 + t69 * t91, t52 * t11 + t30 * t39 - t37 * t31 + t6 * t60 - t68 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t91, -t86, t88, qJD(4) * t82, t60 * t81 + (qJ(5) * t82 - t109) * qJD(4), 0, 0, 0, 0, 0, -t63 * t32 + t38 * t92, t63 * t31 + t39 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, qJ(5) * t79, -0.2e1 * t39 * t31, 0.2e1 * t31 * t38 - 0.2e1 * t39 * t32, 0, 0, 0, t32 * t111, -t31 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t86, 0, t50, 0, 0, 0, 0, 0, t13, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t13, t91, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t32, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t25;
