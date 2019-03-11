% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:43:58
% EndTime: 2019-03-08 21:44:01
% DurationCPUTime: 1.00s
% Computational Cost: add. (791->165), mult. (1980->297), div. (0->0), fcn. (1662->8), ass. (0->106)
t117 = pkin(4) + pkin(8);
t58 = sin(qJ(3));
t107 = t58 * qJ(4);
t61 = cos(qJ(3));
t63 = -pkin(3) - pkin(9);
t73 = -t61 * t63 + t107;
t31 = -pkin(2) - t73;
t41 = t117 * t58;
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t112 = t60 * t31 + t57 * t41;
t54 = t61 ^ 2;
t80 = qJD(3) * (t58 ^ 2 - t54);
t51 = t57 ^ 2;
t111 = -t60 ^ 2 + t51;
t79 = t111 * qJD(5);
t101 = qJD(4) * t61;
t42 = t117 * t61;
t96 = t42 * qJD(5);
t119 = t73 * qJD(3) - t101 - t96;
t118 = 0.2e1 * qJD(4);
t116 = pkin(5) * t60;
t56 = cos(pkin(6));
t62 = cos(qJ(2));
t104 = qJD(2) * t62;
t55 = sin(pkin(6));
t85 = t55 * t104;
t59 = sin(qJ(2));
t114 = t55 * t59;
t91 = t58 * t114;
t11 = -qJD(3) * t91 + (qJD(3) * t56 + t85) * t61;
t26 = t61 * t114 + t56 * t58;
t115 = t26 * t11;
t113 = t55 * t62;
t109 = qJ(4) * t61;
t108 = qJ(6) * t61;
t106 = qJ(6) - t63;
t105 = qJD(2) * t59;
t103 = qJD(3) * t57;
t102 = qJD(3) * t60;
t100 = qJD(5) * t57;
t99 = qJD(5) * t60;
t98 = qJD(5) * t61;
t97 = qJD(5) * t63;
t95 = t58 * qJD(3);
t94 = t60 * qJD(6);
t48 = t61 * qJD(3);
t93 = qJ(4) * qJD(5);
t92 = -0.2e1 * pkin(2) * qJD(3);
t90 = pkin(5) * t100;
t89 = pkin(8) * t95;
t88 = t57 * t98;
t87 = t60 * t98;
t86 = t55 * t105;
t84 = t58 * t48;
t83 = t60 * t95;
t82 = t57 * t99;
t38 = t106 * t60;
t81 = -t31 + t108;
t78 = pkin(3) * t95 - t58 * qJD(4);
t77 = t57 * t83;
t10 = -t60 * t108 + t112;
t34 = t60 * t41;
t9 = t58 * pkin(5) + t81 * t57 + t34;
t76 = t10 * t60 - t57 * t9;
t75 = -t61 * pkin(3) - t107;
t25 = -t56 * t61 + t91;
t13 = t57 * t113 + t25 * t60;
t71 = t60 * t113 - t25 * t57;
t74 = -t13 * t57 - t60 * t71;
t70 = t11 * t57 + t26 * t99;
t69 = -t26 * t100 + t11 * t60;
t17 = (pkin(9) * t58 - t109) * qJD(3) + t78;
t46 = pkin(8) * t48;
t36 = pkin(4) * t48 + t46;
t5 = t31 * t100 - t60 * t17 - t57 * t36 - t41 * t99;
t68 = t83 + t88;
t67 = t57 * t95 - t87;
t35 = t117 * t95;
t66 = -t35 + (-t58 * t63 - t109) * qJD(5);
t65 = t75 * qJD(3) + t101;
t12 = qJD(3) * t26 + t58 * t85;
t64 = t11 * t61 + t12 * t58 + (t25 * t61 - t26 * t58) * qJD(3);
t45 = t57 * pkin(5) + qJ(4);
t44 = pkin(5) * t99 + qJD(4);
t43 = 0.2e1 * t84;
t39 = -pkin(2) + t75;
t37 = t106 * t57;
t30 = t60 * t36;
t28 = -t57 * t48 - t58 * t99;
t27 = -t58 * t100 + t60 * t48;
t24 = t61 * t116 + t42;
t23 = -qJ(4) * t48 + t78;
t21 = -qJD(5) * t38 - t57 * qJD(6);
t20 = t106 * t100 - t94;
t19 = (t61 * t105 + t62 * t95) * t55;
t18 = (t58 * t105 - t62 * t48) * t55;
t15 = -pkin(5) * t88 + (-t116 - t117) * t95;
t8 = t13 * qJD(5) + t12 * t57 + t60 * t86;
t7 = t71 * qJD(5) + t12 * t60 - t57 * t86;
t6 = -t112 * qJD(5) - t57 * t17 + t30;
t4 = t20 * t60 + t21 * t57 + (-t37 * t60 + t38 * t57) * qJD(5);
t3 = t68 * qJ(6) - t61 * t94 - t5;
t2 = pkin(5) * t48 + t30 + t81 * t99 + (-qJ(6) * t95 - qJD(5) * t41 + qJD(6) * t61 - t17) * t57;
t1 = t74 * qJD(5) + t8 * t57 + t7 * t60;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t55 ^ 2 * t59 * t104 + 0.2e1 * t25 * t12 + 0.2e1 * t115, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t13 * t7 - 0.2e1 * t71 * t8 + 0.2e1 * t115; 0, 0, -t86, -t85, 0, 0, 0, 0, 0, -t19, t18, t64, t19, -t18 (t39 * t105 - t23 * t62) * t55 + t64 * pkin(8), 0, 0, 0, 0, 0 (-t26 * t102 + t7) * t58 + (qJD(3) * t13 + t69) * t61 (t26 * t103 - t8) * t58 + (qJD(3) * t71 - t70) * t61, t74 * t95 + (t57 * t7 - t60 * t8 + (t13 * t60 - t57 * t71) * qJD(5)) * t61, t8 * t10 + t11 * t24 + t13 * t2 + t26 * t15 - t3 * t71 + t7 * t9; 0, 0, 0, 0, t43, -0.2e1 * t80, 0, 0, 0, t58 * t92, t61 * t92, 0, 0.2e1 * t23 * t61 - 0.2e1 * t39 * t95, -0.2e1 * t23 * t58 - 0.2e1 * t39 * t48, 0.2e1 * t39 * t23, -0.2e1 * t51 * t84 + 0.2e1 * t54 * t82, -0.2e1 * t54 * t79 - 0.4e1 * t61 * t77, 0.2e1 * t57 * t80 - 0.2e1 * t58 * t87, 0.2e1 * t58 * t88 + 0.2e1 * t60 * t80, t43, 0.2e1 * (-t42 * t102 + t6) * t58 + 0.2e1 * ((-t57 * t31 + t34) * qJD(3) - t35 * t60 - t57 * t96) * t61, 0.2e1 * (t42 * t103 + t5) * t58 + 0.2e1 * (-t112 * qJD(3) + t35 * t57 - t60 * t96) * t61, 0.2e1 * t76 * t95 + 0.2e1 * (t2 * t57 - t3 * t60 + (t10 * t57 + t60 * t9) * qJD(5)) * t61, 0.2e1 * t10 * t3 + 0.2e1 * t24 * t15 + 0.2e1 * t9 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t11, 0, t12, t11, -t12 * pkin(3) + t11 * qJ(4) + t26 * qJD(4), 0, 0, 0, 0, 0, t70, t69, -t1, t11 * t45 + t13 * t20 - t21 * t71 + t26 * t44 - t8 * t37 - t7 * t38; 0, 0, 0, 0, 0, 0, t48, -t95, 0, -t46, t89, t65, t46, -t89, t65 * pkin(8), t61 * t79 + t77, -t111 * t95 + 0.4e1 * t61 * t82, t27, t28, 0, -t119 * t60 + t66 * t57, t119 * t57 + t66 * t60 (-t37 * t95 - t21 * t61 - t2 + (-t38 * t61 - t10) * qJD(5)) * t60 + (t38 * t95 + t20 * t61 - t3 + (-t37 * t61 + t9) * qJD(5)) * t57, t10 * t21 + t15 * t45 - t2 * t38 + t9 * t20 + t24 * t44 - t3 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, qJ(4) * t118, -0.2e1 * t82, 0.2e1 * t79, 0, 0, 0, 0.2e1 * qJD(4) * t57 + 0.2e1 * t60 * t93, 0.2e1 * qJD(4) * t60 - 0.2e1 * t57 * t93, -0.2e1 * t4, -0.2e1 * t38 * t20 - 0.2e1 * t37 * t21 + 0.2e1 * t45 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, t46, 0, 0, 0, 0, 0, t27, t28, 0, t76 * qJD(5) + t2 * t60 + t3 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t68, t48, t6, t5, -t67 * pkin(5), t2 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t99, 0, -t57 * t97, -t60 * t97, t90, t20 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t99, 0, -t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t14;
