% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:06
% EndTime: 2019-03-08 20:25:08
% DurationCPUTime: 0.77s
% Computational Cost: add. (1016->135), mult. (2679->243), div. (0->0), fcn. (2689->12), ass. (0->101)
t75 = cos(qJ(6));
t66 = t75 ^ 2;
t71 = sin(qJ(6));
t108 = t71 ^ 2 - t66;
t92 = t108 * qJD(6);
t121 = qJD(4) + qJD(5);
t67 = sin(pkin(12));
t59 = pkin(2) * t67 + pkin(8);
t120 = pkin(9) + t59;
t119 = cos(qJ(5));
t72 = sin(qJ(5));
t73 = sin(qJ(4));
t114 = t72 * t73;
t76 = cos(qJ(4));
t93 = t119 * qJD(5);
t94 = qJD(4) * t119;
t34 = t121 * t114 + (-t94 - t93) * t76;
t51 = t119 * t73 + t72 * t76;
t118 = t51 * t34;
t117 = t51 * t71;
t116 = t51 * t75;
t35 = t121 * t51;
t115 = t71 * t35;
t113 = t75 * t34;
t112 = t75 * t35;
t46 = t120 * t73;
t47 = t120 * t76;
t31 = t119 * t47 - t46 * t72;
t81 = t120 * t94;
t90 = qJD(4) * t72 * t120;
t12 = qJD(5) * t31 - t73 * t90 + t76 * t81;
t30 = t119 * t46 + t47 * t72;
t64 = qJD(6) * t75;
t111 = t12 * t71 + t30 * t64;
t50 = -t119 * t76 + t114;
t110 = t112 * t51 - t113 * t50;
t63 = -pkin(4) * t119 - pkin(5);
t106 = qJD(5) * t72;
t97 = pkin(4) * t106;
t109 = t63 * t64 + t71 * t97;
t68 = sin(pkin(6));
t107 = qJD(2) * t68;
t105 = qJD(6) * t71;
t104 = t73 * qJD(4);
t103 = t76 * qJD(4);
t102 = t71 * t113;
t101 = 0.2e1 * t103;
t100 = pkin(5) * t105;
t99 = pkin(5) * t64;
t98 = pkin(4) * t104;
t96 = t51 * t105;
t95 = t71 * t64;
t69 = cos(pkin(12));
t60 = -pkin(2) * t69 - pkin(3);
t91 = pkin(4) * t93;
t74 = sin(qJ(2));
t77 = cos(qJ(2));
t44 = (t67 * t77 + t69 * t74) * t68;
t70 = cos(pkin(6));
t36 = -t44 * t73 + t70 * t76;
t37 = t44 * t76 + t70 * t73;
t16 = t119 * t37 + t36 * t72;
t83 = t67 * t74 - t69 * t77;
t43 = t83 * t68;
t89 = t16 * t75 + t43 * t71;
t88 = t16 * t71 - t43 * t75;
t54 = -pkin(4) * t76 + t60;
t29 = pkin(5) * t50 - pkin(10) * t51 + t54;
t87 = t29 * t75 - t31 * t71;
t86 = t29 * t71 + t31 * t75;
t85 = t34 * t50 - t35 * t51;
t62 = pkin(4) * t72 + pkin(10);
t84 = t50 * t62 - t51 * t63;
t82 = t105 * t63 - t75 * t97;
t80 = -t34 * t71 + t51 * t64;
t22 = t96 + t113;
t21 = t105 * t50 - t112;
t39 = t83 * t107;
t79 = qJD(4) * t37 - t39 * t73;
t78 = -t34 * t63 - t35 * t62 + (-t119 * t50 + t51 * t72) * qJD(5) * pkin(4);
t56 = 0.2e1 * t95;
t49 = -0.2e1 * t92;
t48 = t51 ^ 2;
t38 = qJD(2) * t44;
t27 = t30 * t105;
t23 = t50 * t64 + t115;
t20 = qJD(4) * t36 - t39 * t76;
t15 = -t119 * t36 + t37 * t72;
t14 = pkin(5) * t35 + pkin(10) * t34 + t98;
t13 = -t51 * t92 - t102;
t11 = t106 * t47 + t46 * t93 + t73 * t81 + t76 * t90;
t9 = t108 * t34 - 0.4e1 * t51 * t95;
t8 = qJD(5) * t16 + t119 * t79 + t72 * t20;
t7 = t106 * t37 - t119 * t20 - t36 * t93 + t72 * t79;
t6 = -qJD(6) * t86 + t71 * t11 + t75 * t14;
t5 = -qJD(6) * t87 + t75 * t11 - t71 * t14;
t4 = t105 * t15 - t75 * t8;
t3 = t15 * t64 + t71 * t8;
t2 = -qJD(6) * t89 + t38 * t75 + t71 * t7;
t1 = qJD(6) * t88 - t38 * t71 + t75 * t7;
t10 = [0, 0, 0, 0, 0.2e1 * t38 * t43 - 0.2e1 * t39 * t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t74 * t107, -t77 * t107 (-t38 * t69 - t39 * t67) * pkin(2), 0, 0, 0, 0, 0, t104 * t43 - t38 * t76, t103 * t43 + t38 * t73, 0, 0, 0, 0, 0, t35 * t43 + t38 * t50, -t34 * t43 + t38 * t51, 0, 0, 0, 0, 0, t117 * t8 + t15 * t80 + t2 * t50 - t35 * t88, t1 * t50 + t116 * t8 - t15 * t22 - t35 * t89; 0, 0, 0, 0, 0, t73 * t101, 0.2e1 * (-t73 ^ 2 + t76 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * t60 * t104, t60 * t101, -0.2e1 * t118, 0.2e1 * t85, 0, 0, 0, 0.2e1 * t35 * t54 + 0.2e1 * t50 * t98, -0.2e1 * t34 * t54 + 0.2e1 * t51 * t98, -0.2e1 * t118 * t66 - 0.2e1 * t48 * t95, 0.4e1 * t102 * t51 + 0.2e1 * t48 * t92, -0.2e1 * t50 * t96 + 0.2e1 * t110, -0.2e1 * t115 * t51 - 0.2e1 * t50 * t80, 0.2e1 * t50 * t35, 0.2e1 * t117 * t12 + 0.2e1 * t30 * t80 + 0.2e1 * t35 * t87 + 0.2e1 * t6 * t50, 0.2e1 * t116 * t12 - 0.2e1 * t22 * t30 - 0.2e1 * t35 * t86 + 0.2e1 * t5 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75 * t85 + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t20, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, t103, -t104, 0, -t59 * t103, t59 * t104, 0, 0, -t34, -t35, 0, -t12, t11, t13, t9, t23, -t21, 0, t27 + (-qJD(6) * t84 - t12) * t75 + t78 * t71, t105 * t84 + t75 * t78 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t103, 0, 0, 0, 0, 0, -t35, t34, 0, 0, 0, 0, 0, t21, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t97, -0.2e1 * t91, t56, t49, 0, 0, 0, 0.2e1 * t82, 0.2e1 * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t35, 0, -t12, t11, t13, t9, t23, -t21, 0, t27 + (pkin(5) * t34 - pkin(10) * t35) * t71 + (-t12 + (-pkin(5) * t51 - pkin(10) * t50) * qJD(6)) * t75, pkin(5) * t22 + pkin(10) * t21 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t34, 0, 0, 0, 0, 0, t21, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t91, t56, t49, 0, 0, 0, t82 - t100, -t99 + t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t49, 0, 0, 0, -0.2e1 * t100, -0.2e1 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t80, t35, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t105, 0, -t62 * t64 - t71 * t91, t105 * t62 - t75 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t105, 0, -pkin(10) * t64, pkin(10) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
