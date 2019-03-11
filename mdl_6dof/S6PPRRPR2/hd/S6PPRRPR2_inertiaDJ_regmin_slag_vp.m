% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPRRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:50:54
% EndTime: 2019-03-08 18:50:57
% DurationCPUTime: 0.82s
% Computational Cost: add. (597->142), mult. (1940->270), div. (0->0), fcn. (1979->12), ass. (0->109)
t56 = cos(qJ(4));
t45 = t56 ^ 2;
t53 = sin(qJ(4));
t81 = qJD(4) * (t53 ^ 2 - t45);
t52 = sin(qJ(6));
t42 = t52 ^ 2;
t55 = cos(qJ(6));
t112 = -t55 ^ 2 + t42;
t80 = t112 * qJD(6);
t109 = t53 * qJ(5);
t58 = -pkin(4) - pkin(10);
t72 = -t56 * t58 + t109;
t102 = qJD(5) * t56;
t116 = pkin(5) + pkin(9);
t36 = t116 * t56;
t94 = t36 * qJD(6);
t118 = t72 * qJD(4) - t102 - t94;
t117 = 0.2e1 * qJD(5);
t47 = sin(pkin(7));
t54 = sin(qJ(3));
t115 = t47 * t54;
t57 = cos(qJ(3));
t114 = t47 * t57;
t49 = cos(pkin(12));
t50 = cos(pkin(7));
t113 = t49 * t50;
t110 = qJ(5) * t56;
t108 = qJD(3) * t54;
t107 = qJD(3) * t57;
t106 = qJD(4) * t52;
t105 = qJD(4) * t53;
t104 = qJD(4) * t55;
t41 = qJD(4) * t56;
t103 = qJD(4) * t57;
t101 = qJD(6) * t52;
t100 = qJD(6) * t55;
t99 = qJD(6) * t56;
t98 = qJD(6) * t58;
t48 = sin(pkin(6));
t51 = cos(pkin(6));
t23 = -t48 * t49 * t47 + t51 * t50;
t46 = sin(pkin(12));
t62 = t51 * t115 + (t54 * t113 + t46 * t57) * t48;
t12 = t23 * t53 + t62 * t56;
t97 = t12 * qJD(6);
t91 = t57 * t113;
t16 = -t51 * t114 + (t46 * t54 - t91) * t48;
t96 = t16 * qJD(4);
t25 = t56 * t115 + t53 * t50;
t95 = t25 * qJD(6);
t93 = qJ(5) * qJD(6);
t92 = -0.2e1 * pkin(3) * qJD(4);
t90 = t53 * t115;
t89 = pkin(9) * t105;
t88 = t52 * t99;
t87 = t55 * t99;
t86 = t47 * t108;
t85 = t47 * t107;
t84 = t53 * t41;
t83 = t53 * t104;
t82 = t52 * t100;
t79 = pkin(4) * t105 - t53 * qJD(5);
t78 = t52 * t83;
t77 = -t56 * pkin(4) - t109;
t59 = t62 * t53;
t11 = -t23 * t56 + t59;
t76 = t11 * t55 - t16 * t52;
t75 = t11 * t52 + t16 * t55;
t28 = -pkin(3) - t72;
t35 = t116 * t53;
t74 = t55 * t28 + t52 * t35;
t73 = t52 * t28 - t55 * t35;
t24 = -t56 * t50 + t90;
t70 = t55 * t114 - t52 * t24;
t69 = t52 * t114 + t55 * t24;
t14 = -t51 * t85 + (-qJD(3) * t91 + t108 * t46) * t48;
t4 = -qJD(4) * t59 + (t23 * qJD(4) - t14) * t56;
t68 = t4 * t52 + t55 * t97;
t67 = t4 * t55 - t52 * t97;
t17 = qJD(4) * t90 - t50 * t41 - t56 * t85;
t66 = -t17 * t52 + t55 * t95;
t65 = -t17 * t55 - t52 * t95;
t30 = t116 * t105;
t64 = -t30 + (-t53 * t58 - t110) * qJD(6);
t63 = t77 * qJD(4) + t102;
t3 = qJD(4) * t12 - t14 * t53;
t61 = t3 * t53 + t4 * t56 + (t11 * t56 - t12 * t53) * qJD(4);
t18 = qJD(4) * t25 + t53 * t85;
t60 = -t17 * t56 + t18 * t53 + (t24 * t56 - t25 * t53) * qJD(4);
t39 = pkin(9) * t41;
t37 = 0.2e1 * t84;
t33 = -pkin(3) + t77;
t31 = pkin(5) * t41 + t39;
t27 = -t53 * t100 - t52 * t41;
t26 = -t53 * t101 + t55 * t41;
t22 = -qJ(5) * t41 + t79;
t21 = (t53 * t103 + t56 * t108) * t47;
t20 = (-t56 * t103 + t53 * t108) * t47;
t19 = (pkin(10) * t53 - t110) * qJD(4) + t79;
t15 = t62 * qJD(3);
t10 = t70 * qJD(6) + t55 * t18 - t52 * t86;
t9 = t69 * qJD(6) + t52 * t18 + t55 * t86;
t8 = -t74 * qJD(6) - t52 * t19 + t55 * t31;
t7 = t73 * qJD(6) - t55 * t19 - t52 * t31;
t6 = t15 * t56 - t53 * t96;
t5 = t15 * t53 + t56 * t96;
t2 = t76 * qJD(6) + t15 * t55 + t3 * t52;
t1 = -t75 * qJD(6) - t15 * t52 + t3 * t55;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t11 * t3 + 0.2e1 * t12 * t4 + 0.2e1 * t16 * t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t18 - t12 * t17 + t3 * t24 + t4 * t25 + (t16 * t108 - t15 * t57) * t47, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t47 ^ 2 * t54 * t107 - 0.2e1 * t25 * t17 + 0.2e1 * t24 * t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t15, t14, 0, 0, 0, 0, 0, -t6, t5, t61, t6, -t5, t61 * pkin(9) + t15 * t33 + t16 * t22, 0, 0, 0, 0, 0 (-t12 * t104 + t1) * t53 + (t76 * qJD(4) + t67) * t56 (t12 * t106 - t2) * t53 + (-t75 * qJD(4) - t68) * t56; 0, 0, 0, -t86, -t85, 0, 0, 0, 0, 0, -t21, t20, t60, t21, -t20 (t33 * t108 - t22 * t57) * t47 + t60 * pkin(9), 0, 0, 0, 0, 0 (-t25 * t104 + t10) * t53 + (t69 * qJD(4) + t65) * t56 (t25 * t106 - t9) * t53 + (t70 * qJD(4) - t66) * t56; 0, 0, 0, 0, 0, t37, -0.2e1 * t81, 0, 0, 0, t53 * t92, t56 * t92, 0, -0.2e1 * t33 * t105 + 0.2e1 * t22 * t56, -0.2e1 * t22 * t53 - 0.2e1 * t33 * t41, 0.2e1 * t33 * t22, -0.2e1 * t42 * t84 + 0.2e1 * t45 * t82, -0.2e1 * t45 * t80 - 0.4e1 * t56 * t78, 0.2e1 * t52 * t81 - 0.2e1 * t53 * t87, 0.2e1 * t53 * t88 + 0.2e1 * t55 * t81, t37, 0.2e1 * (-t36 * t104 + t8) * t53 + 0.2e1 * (-t73 * qJD(4) - t30 * t55 - t52 * t94) * t56, 0.2e1 * (t36 * t106 + t7) * t53 + 0.2e1 * (-t74 * qJD(4) + t30 * t52 - t55 * t94) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t3, t4, -t3 * pkin(4) + t4 * qJ(5) + t12 * qJD(5), 0, 0, 0, 0, 0, t68, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, 0, t18, -t17, -t18 * pkin(4) - t17 * qJ(5) + t25 * qJD(5), 0, 0, 0, 0, 0, t66, t65; 0, 0, 0, 0, 0, 0, 0, t41, -t105, 0, -t39, t89, t63, t39, -t89, t63 * pkin(9), t56 * t80 + t78, -t112 * t105 + 0.4e1 * t56 * t82, t26, t27, 0, -t118 * t55 + t64 * t52, t118 * t52 + t64 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, qJ(5) * t117, -0.2e1 * t82, 0.2e1 * t80, 0, 0, 0, 0.2e1 * qJD(5) * t52 + 0.2e1 * t55 * t93, 0.2e1 * qJD(5) * t55 - 0.2e1 * t52 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, t39, 0, 0, 0, 0, 0, t26, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t105 - t87, t83 + t88, t41, t8, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t100, 0, -t52 * t98, -t55 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t13;
