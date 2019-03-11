% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:19
% EndTime: 2019-03-09 02:06:23
% DurationCPUTime: 1.14s
% Computational Cost: add. (1062->171), mult. (2063->295), div. (0->0), fcn. (1674->6), ass. (0->110)
t50 = sin(pkin(9));
t51 = cos(pkin(9));
t56 = -pkin(1) - pkin(2);
t114 = t51 * qJ(2) + t50 * t56;
t31 = -pkin(7) + t114;
t53 = sin(qJ(4));
t44 = t53 * qJD(4);
t105 = t51 * qJD(2);
t55 = cos(qJ(4));
t83 = t55 * t105;
t127 = t31 * t44 - t83;
t80 = -t50 * qJ(2) + t51 * t56;
t30 = pkin(3) - t80;
t121 = t53 * pkin(8);
t75 = t55 * pkin(4) + t121;
t16 = t30 + t75;
t52 = sin(qJ(5));
t13 = t52 * t16;
t54 = cos(qJ(5));
t116 = t54 * t55;
t20 = t31 * t116;
t66 = t20 + t13;
t46 = t52 ^ 2;
t48 = t54 ^ 2;
t79 = qJD(5) * (t46 - t48);
t103 = t55 * qJD(4);
t62 = -t31 * t103 - t53 * t105;
t122 = pkin(8) * t55;
t74 = pkin(4) * t53 - t122;
t126 = t74 * qJD(5) + t62;
t106 = t50 * qJD(2);
t22 = -t74 * qJD(4) + t106;
t4 = -t66 * qJD(5) + t127 * t52 + t54 * t22;
t125 = 0.2e1 * qJD(2);
t124 = 2 * qJD(6);
t123 = pkin(5) * t52;
t104 = t54 * qJD(6);
t82 = qJ(6) * t103;
t115 = t53 * t104 + t54 * t82;
t70 = -t54 * pkin(5) - t52 * qJ(6);
t63 = t70 * qJD(5);
t6 = (t31 - t123) * t103 + (t63 + t105) * t53 + t115;
t120 = t6 * t54;
t69 = -qJ(6) * t54 + t123;
t21 = t69 * qJD(5) - t52 * qJD(6);
t119 = t21 * t53;
t118 = t31 * t54;
t32 = -pkin(4) + t70;
t117 = t32 * t53;
t47 = t53 ^ 2;
t112 = -t55 ^ 2 + t47;
t9 = (t31 - t69) * t53;
t111 = qJD(4) * t9;
t110 = qJD(4) * t54;
t109 = qJD(5) * t52;
t108 = qJD(5) * t53;
t45 = qJD(5) * t54;
t107 = qJD(5) * t55;
t102 = -0.2e1 * pkin(4) * qJD(5);
t101 = t50 * t116;
t100 = t16 * t45 + t52 * t22 + t54 * t83;
t99 = pkin(5) * t44;
t98 = pkin(8) * t109;
t97 = pkin(8) * t45;
t96 = t31 * t108;
t95 = t52 * t107;
t94 = t54 * t107;
t93 = t31 * t109;
t92 = t47 * t45;
t91 = t48 * t103;
t90 = t52 * t44;
t89 = t52 * t45;
t88 = t53 * t103;
t86 = t50 * t44;
t85 = t54 * t44;
t84 = t54 * t103;
t81 = qJ(6) * t108;
t78 = t112 * qJD(4);
t77 = t52 * t88;
t76 = t52 * t84;
t7 = t55 * qJ(6) + t66;
t8 = -t54 * t16 + (t31 * t52 - pkin(5)) * t55;
t73 = t52 * t8 + t54 * t7;
t72 = t52 * t7 - t54 * t8;
t71 = -t117 - t122;
t23 = t52 * t55 * t50 + t54 * t51;
t24 = -t52 * t51 + t101;
t68 = t23 * t54 - t24 * t52;
t67 = t23 * t52 + t24 * t54;
t65 = t32 * t103 + t119;
t25 = t52 * t108 - t84;
t61 = t85 + t95;
t27 = t52 * t103 + t53 * t45;
t60 = t63 + t104;
t1 = (qJD(6) - t93) * t55 + (-qJ(6) - t118) * t44 + t100;
t2 = -t4 + t99;
t59 = -t72 * qJD(5) + t1 * t54 + t2 * t52;
t11 = -t51 * t109 + t50 * t94 - t52 * t86;
t58 = -t11 * t55 + t23 * t44 + (-0.2e1 * t77 - t92) * t50;
t10 = t23 * qJD(5) + t50 * t85;
t57 = t68 * qJD(5) - t10 * t54 + t11 * t52;
t39 = t46 * t103;
t37 = pkin(8) * t90;
t35 = t48 * t88;
t28 = -t90 + t94;
t15 = t25 * t50;
t14 = t27 * t50;
t5 = -t47 * t50 * t109 - t10 * t55 + (-t24 + 0.2e1 * t101) * t44;
t3 = t61 * t31 - t100;
t12 = [0, 0, 0, 0, t125, qJ(2) * t125, 0.2e1 * t106, 0.2e1 * t105 (t114 * t51 - t80 * t50) * t125, 0.2e1 * t88, -0.2e1 * t78, 0, 0, 0, 0.2e1 * t55 * t106 - 0.2e1 * t30 * t44, -0.2e1 * t30 * t103 - 0.2e1 * t53 * t106, -0.2e1 * t47 * t89 + 0.2e1 * t35, 0.2e1 * t47 * t79 - 0.4e1 * t53 * t76, 0.2e1 * t112 * t110 + 0.2e1 * t53 * t95, -0.2e1 * t52 * t78 + 0.2e1 * t53 * t94, -0.2e1 * t88, -0.2e1 * t47 * t52 * t105 - 0.2e1 * t16 * t85 + 0.2e1 * t4 * t55 + 0.2e1 * (-t77 - t92) * t31, 0.2e1 * t3 * t55 + 0.2e1 * (-t54 * t105 + t93) * t47 + 0.2e1 * (-t20 + t13) * t44, 0.2e1 * (-t52 * t111 - t2) * t55 + 0.2e1 * (qJD(4) * t8 - t9 * t45 - t6 * t52) * t53, 0.2e1 * t72 * t103 + 0.2e1 * (t73 * qJD(5) + t1 * t52 - t2 * t54) * t53, 0.2e1 * (t9 * t110 + t1) * t55 + 0.2e1 * (-qJD(4) * t7 - t9 * t109 + t120) * t53, 0.2e1 * t7 * t1 + 0.2e1 * t8 * t2 + 0.2e1 * t9 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51 * t44, t51 * t103, 0, 0, 0, 0, 0, t58, -t5, t58, -t68 * t103 + (t67 * qJD(5) - t10 * t52 - t11 * t54) * t53, t5, t1 * t24 - t7 * t10 + t8 * t11 + t2 * t23 + (t103 * t9 + t53 * t6) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t50 ^ 2 * t88 - 0.2e1 * t24 * t10 + 0.2e1 * t23 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t73 * qJD(4) - t6) * t55 + (t59 + t111) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t53 + (t112 * t50 + t67 * t55) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t35 + 0.2e1 * (t46 - 0.1e1) * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, t44, 0, t62, t127, t53 * t79 - t76, 0.4e1 * t53 * t89 + t39 - t91, t28, -t61, 0, t37 + (pkin(4) * t103 + t96) * t52 + t126 * t54 (t75 * qJD(4) + t96) * t54 - t126 * t52, -t120 + t37 - t65 * t52 + (t52 * t9 + t71 * t54) * qJD(5), t59 (t71 * qJD(5) - t6) * t52 + (-qJD(5) * t9 + t119 + (t32 * t55 - t121) * qJD(4)) * t54, pkin(8) * t59 + t9 * t21 + t6 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50 * t103, t86, 0, 0, 0, 0, 0, t15, t14, t15, t57, -t14, pkin(8) * t57 + t50 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t103, 0, 0, 0, 0, 0, -t61, -t28, -t61, t39 + t91, t28, -t55 * t21 + (t117 + (t46 + t48) * t122) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t89, -0.2e1 * t79, 0, 0, 0, t52 * t102, t54 * t102, 0.2e1 * t32 * t109 - 0.2e1 * t21 * t54, 0, -0.2e1 * t21 * t52 - 0.2e1 * t32 * t45, 0.2e1 * t32 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t27, -t44, t4, t3, t4 - 0.2e1 * t99 (pkin(5) * t103 + t81) * t54 + (t82 + (-pkin(5) * qJD(5) + qJD(6)) * t53) * t52 (t124 - t93) * t55 + (-0.2e1 * qJ(6) - t118) * t44 + t100, -t2 * pkin(5) + t1 * qJ(6) + t7 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, -t11, 0, -t10, -t11 * pkin(5) - t10 * qJ(6) + t24 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t25, -t27, 0, -t25, -pkin(5) * t27 - t52 * t81 + t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t109, 0, -t97, t98, -t97, t60, -t98, t60 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, qJ(6) * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t25, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
