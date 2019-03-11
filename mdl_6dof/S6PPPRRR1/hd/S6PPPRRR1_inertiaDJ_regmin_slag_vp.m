% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x20]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPPRRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_inertiaDJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:35
% EndTime: 2019-03-08 18:40:37
% DurationCPUTime: 0.78s
% Computational Cost: add. (785->123), mult. (2570->255), div. (0->0), fcn. (2985->16), ass. (0->107)
t59 = cos(qJ(6));
t44 = t59 ^ 2;
t56 = sin(qJ(6));
t106 = t56 ^ 2 - t44;
t81 = t106 * qJD(6);
t114 = pkin(10) * t56;
t52 = cos(pkin(13));
t54 = cos(pkin(7));
t108 = t52 * t54;
t49 = sin(pkin(7));
t55 = cos(pkin(6));
t110 = t49 * t55;
t46 = sin(pkin(14));
t47 = sin(pkin(13));
t50 = sin(pkin(6));
t51 = cos(pkin(14));
t26 = t51 * t110 + (t51 * t108 - t46 * t47) * t50;
t53 = cos(pkin(8));
t113 = t26 * t53;
t48 = sin(pkin(8));
t58 = sin(qJ(4));
t112 = t48 * t58;
t61 = cos(qJ(4));
t111 = t48 * t61;
t109 = t51 * t53;
t60 = cos(qJ(5));
t107 = t59 * t60;
t57 = sin(qJ(5));
t43 = t57 ^ 2;
t105 = -t60 ^ 2 + t43;
t104 = qJD(4) * t58;
t103 = qJD(4) * t61;
t102 = qJD(5) * t56;
t101 = qJD(5) * t59;
t100 = qJD(6) * t56;
t99 = qJD(6) * t59;
t98 = qJD(6) * t60;
t29 = t54 * t112 + (t58 * t109 + t46 * t61) * t49;
t34 = -t48 * t49 * t51 + t53 * t54;
t18 = t57 * t29 - t60 * t34;
t97 = t18 * qJD(6);
t91 = t57 * t112;
t36 = -t60 * t53 + t91;
t96 = t36 * qJD(6);
t95 = t57 * qJD(5);
t94 = t60 * qJD(5);
t93 = -0.2e1 * pkin(4) * qJD(5);
t92 = -0.2e1 * pkin(5) * qJD(6);
t90 = t56 * t95;
t89 = t56 * t98;
t88 = t48 * t104;
t87 = t48 * t103;
t86 = t56 * t99;
t85 = t57 * t94;
t84 = t59 * t95;
t83 = t59 * t94;
t82 = t59 * t98;
t80 = t105 * qJD(5);
t79 = t57 * t83;
t78 = -t60 * pkin(5) - t57 * pkin(11);
t77 = pkin(5) * t57 - pkin(11) * t60;
t27 = t50 * t47 * t51 + (t50 * t108 + t110) * t46;
t35 = -t50 * t52 * t49 + t55 * t54;
t72 = t35 * t48 + t113;
t11 = t27 * t58 - t72 * t61;
t12 = t27 * t61 + t72 * t58;
t17 = -t26 * t48 + t35 * t53;
t8 = t12 * t60 + t17 * t57;
t76 = t11 * t59 - t8 * t56;
t75 = t11 * t56 + t8 * t59;
t7 = t12 * t57 - t17 * t60;
t19 = t60 * t29 + t57 * t34;
t28 = -t54 * t111 + (-t61 * t109 + t46 * t58) * t49;
t74 = t59 * t19 + t56 * t28;
t73 = t56 * t19 - t59 * t28;
t9 = -t103 * t113 + t27 * t104 - t35 * t87;
t3 = t8 * qJD(5) - t9 * t57;
t71 = t3 * t56 + t7 * t99;
t70 = t7 * t100 - t3 * t59;
t37 = t60 * t112 + t57 * t53;
t69 = t59 * t111 + t56 * t37;
t68 = t56 * t111 - t59 * t37;
t24 = t28 * qJD(4);
t14 = t19 * qJD(5) - t57 * t24;
t67 = t14 * t56 + t59 * t97;
t66 = -t14 * t59 + t56 * t97;
t31 = t37 * qJD(5) + t57 * t87;
t65 = t31 * t56 + t59 * t96;
t64 = -t31 * t59 + t56 * t96;
t63 = t84 + t89;
t62 = -t82 + t90;
t40 = -pkin(4) + t78;
t38 = t77 * qJD(5);
t30 = qJD(5) * t91 - t53 * t94 - t60 * t87;
t25 = t29 * qJD(4);
t21 = t62 * pkin(10) - t40 * t100 + t59 * t38;
t20 = t63 * pkin(10) - t56 * t38 - t40 * t99;
t16 = t68 * qJD(6) + t56 * t30 + t59 * t88;
t15 = t69 * qJD(6) + t59 * t30 - t56 * t88;
t13 = t60 * t24 + t29 * t95 - t34 * t94;
t10 = t12 * qJD(4);
t6 = -t74 * qJD(6) + t56 * t13 + t59 * t25;
t5 = t73 * qJD(6) + t59 * t13 - t56 * t25;
t4 = -qJD(5) * t7 - t9 * t60;
t2 = t76 * qJD(6) + t10 * t56 + t4 * t59;
t1 = -t75 * qJD(6) + t10 * t59 - t4 * t56;
t22 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t10, t9, 0, 0, 0, 0, 0, -t10 * t60 + t11 * t95, t10 * t57 + t11 * t94, 0, 0, 0, 0, 0 (t102 * t7 - t1) * t60 + (qJD(5) * t76 + t71) * t57 (t101 * t7 + t2) * t60 + (-qJD(5) * t75 - t70) * t57; 0, 0, 0, 0, -t25, t24, 0, 0, 0, 0, 0, -t25 * t60 + t28 * t95, t25 * t57 + t28 * t94, 0, 0, 0, 0, 0 (t102 * t18 - t6) * t60 + (-qJD(5) * t73 + t67) * t57 (t101 * t18 - t5) * t60 + (-qJD(5) * t74 - t66) * t57; 0, 0, 0, 0, -t88, -t87, 0, 0, 0, 0, 0 (-t104 * t60 - t61 * t95) * t48 (t104 * t57 - t61 * t94) * t48, 0, 0, 0, 0, 0 (t102 * t36 - t16) * t60 + (-qJD(5) * t69 + t65) * t57 (t101 * t36 - t15) * t60 + (qJD(5) * t68 - t64) * t57; 0, 0, 0, 0, 0, 0, 0.2e1 * t85, -0.2e1 * t80, 0, 0, 0, t57 * t93, t60 * t93, -0.2e1 * t43 * t86 + 0.2e1 * t44 * t85, 0.2e1 * t43 * t81 - 0.4e1 * t56 * t79, 0.2e1 * t101 * t105 + 0.2e1 * t57 * t89, -0.2e1 * t56 * t80 + 0.2e1 * t57 * t82, -0.2e1 * t85, 0.2e1 * t40 * t84 - 0.2e1 * t21 * t60 + 0.2e1 * (t43 * t99 + t56 * t85) * pkin(10), -0.2e1 * t40 * t90 - 0.2e1 * t20 * t60 + 0.2e1 * (-t100 * t43 + t79) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, t70, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t13, 0, 0, 0, 0, 0, t66, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t30, 0, 0, 0, 0, 0, t64, t65; 0, 0, 0, 0, 0, 0, 0, 0, t94, -t95, 0, -pkin(10) * t94, pkin(10) * t95, t56 * t83 - t57 * t81, -t106 * t94 - 0.4e1 * t57 * t86, t62, t63, 0 (pkin(11) * t107 + (-pkin(5) * t59 + t114) * t57) * qJD(6) + (-pkin(10) * t107 + t56 * t78) * qJD(5) (pkin(10) * t57 * t59 + t56 * t77) * qJD(6) + (t60 * t114 + t59 * t78) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t86, -0.2e1 * t81, 0, 0, 0, t56 * t92, t59 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100 * t57 + t83, -t56 * t94 - t57 * t99, t95, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, -t100, 0, -pkin(11) * t99, pkin(11) * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t22;
