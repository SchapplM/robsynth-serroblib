% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:30:07
% EndTime: 2019-03-09 15:30:10
% DurationCPUTime: 1.01s
% Computational Cost: add. (1321->161), mult. (2984->270), div. (0->0), fcn. (2584->6), ass. (0->98)
t69 = cos(qJ(6));
t64 = t69 ^ 2;
t66 = sin(qJ(6));
t105 = t66 ^ 2 - t64;
t89 = qJD(6) * t105;
t115 = qJD(2) + qJD(3);
t72 = 2 * qJD(4);
t114 = -pkin(4) - pkin(9);
t113 = -pkin(8) - pkin(7);
t67 = sin(qJ(3));
t68 = sin(qJ(2));
t108 = t67 * t68;
t111 = cos(qJ(3));
t70 = cos(qJ(2));
t94 = t111 * t70;
t38 = -t94 + t108;
t104 = qJD(3) * t67;
t107 = t67 * t70;
t43 = t113 * t70;
t87 = t113 * t111;
t79 = qJD(2) * t87;
t82 = t68 * t87;
t93 = t113 * qJD(2);
t12 = -qJD(3) * t82 - t43 * t104 - t93 * t107 - t68 * t79;
t39 = t111 * t68 + t107;
t25 = t115 * t39;
t6 = -t25 * qJ(5) - t38 * qJD(5) + t12;
t112 = t6 * t38;
t5 = t6 * t69;
t56 = -t70 * pkin(2) - pkin(1);
t92 = t111 * qJD(3);
t88 = pkin(2) * t92;
t45 = t88 + qJD(4);
t53 = t67 * pkin(2) + qJ(4);
t110 = t53 * t45;
t109 = t66 * t25;
t106 = t69 * t25;
t27 = t113 * t108 - t111 * t43;
t19 = t38 * qJ(5) + t27;
t103 = qJD(6) * t19;
t58 = qJD(6) * t66;
t102 = qJD(6) * t69;
t101 = t38 * qJD(4);
t100 = t68 * qJD(2);
t99 = t70 * qJD(2);
t98 = -0.2e1 * pkin(1) * qJD(2);
t97 = t66 * t106;
t96 = pkin(2) * t100;
t57 = pkin(2) * t104;
t95 = t66 * t102;
t91 = -t38 * pkin(3) + t39 * qJ(4) - t56;
t50 = pkin(5) + t53;
t65 = qJ(4) + pkin(5);
t90 = qJD(6) * (-t50 - t65);
t55 = -t111 * pkin(2) - pkin(3);
t10 = t39 * pkin(5) + t114 * t38 + t91;
t26 = -t67 * t43 - t82;
t18 = -t39 * qJ(5) + t26;
t86 = t69 * t10 - t66 * t18;
t85 = t66 * t10 + t69 * t18;
t51 = -pkin(4) + t55;
t47 = -pkin(9) + t51;
t84 = t38 * t50 - t39 * t47;
t62 = -pkin(3) + t114;
t83 = t38 * t65 - t39 * t62;
t81 = -t25 * qJ(4) - t101;
t80 = t45 * qJ(4) + t53 * qJD(4);
t24 = -qJD(2) * t94 + t108 * t115 - t70 * t92;
t15 = -t39 * t102 + t66 * t24;
t14 = t69 * t24 + t39 * t58;
t78 = t38 * t102 + t109;
t77 = -t38 * t58 + t106;
t9 = t25 * pkin(3) + t24 * qJ(4) - t39 * qJD(4) + t96;
t76 = -t45 * t38 + t39 * t57;
t75 = t53 * t25 - t76;
t74 = t24 * t62 + t25 * t65 + t101 - t103;
t73 = t24 * t47 + t25 * t50 - t103 - t76;
t13 = -t43 * t92 - t70 * t79 + (qJD(3) * t113 + t93) * t108;
t71 = -pkin(3) - pkin(4);
t61 = qJ(4) * t72;
t59 = qJD(4) * t69;
t52 = -0.2e1 * t57;
t46 = t72 + t88;
t44 = 0.2e1 * t95;
t42 = 0.2e1 * t45;
t41 = t45 * t69;
t37 = -0.2e1 * t89;
t36 = t38 ^ 2;
t17 = -0.2e1 * t39 * t24;
t16 = -t38 * pkin(4) + t91;
t11 = t38 * t89 - t97;
t8 = t105 * t25 + 0.4e1 * t38 * t95;
t7 = t24 * qJ(5) - t39 * qJD(5) + t13;
t4 = -t25 * pkin(4) - t9;
t3 = -t24 * pkin(5) + t114 * t25 - t9;
t2 = -t85 * qJD(6) + t69 * t3 - t66 * t7;
t1 = -t86 * qJD(6) - t66 * t3 - t69 * t7;
t20 = [0, 0, 0, 0.2e1 * t68 * t99, 0.2e1 * (-t68 ^ 2 + t70 ^ 2) * qJD(2), 0, 0, 0, t68 * t98, t70 * t98, t17, 0.2e1 * t24 * t38 - 0.2e1 * t39 * t25, 0, 0, 0, 0.2e1 * t56 * t25 + 0.2e1 * t38 * t96, -0.2e1 * t56 * t24 + 0.2e1 * t39 * t96, -0.2e1 * t25 * t91 + 0.2e1 * t9 * t38, 0.2e1 * t12 * t38 + 0.2e1 * t13 * t39 - 0.2e1 * t26 * t24 - 0.2e1 * t27 * t25, -0.2e1 * t24 * t91 - 0.2e1 * t9 * t39, -0.2e1 * t27 * t12 + 0.2e1 * t26 * t13 - 0.2e1 * t9 * t91, -0.2e1 * t16 * t24 + 0.2e1 * t4 * t39, 0.2e1 * t16 * t25 + 0.2e1 * t4 * t38, 0.2e1 * t18 * t24 + 0.2e1 * t19 * t25 - 0.2e1 * t7 * t39 - 0.2e1 * t112, 0.2e1 * t16 * t4 + 0.2e1 * t18 * t7 - 0.2e1 * t19 * t6, 0.2e1 * t64 * t38 * t25 - 0.2e1 * t36 * t95, 0.2e1 * t36 * t89 - 0.4e1 * t38 * t97, 0.2e1 * t39 * t106 - 0.2e1 * t14 * t38, -0.2e1 * t39 * t109 + 0.2e1 * t15 * t38, t17, -0.2e1 * t66 * t112 + 0.2e1 * t78 * t19 + 0.2e1 * t2 * t39 - 0.2e1 * t86 * t24, 0.2e1 * t1 * t39 + 0.2e1 * t77 * t19 + 0.2e1 * t85 * t24 - 0.2e1 * t38 * t5; 0, 0, 0, 0, 0, t99, -t100, 0, -pkin(7) * t99, pkin(7) * t100, 0, 0, -t24, -t25, 0, -t13, t12, -t13, -t55 * t24 - t75, -t12, -t12 * t53 + t13 * t55 + t26 * t57 + t27 * t45, -t6, t7, t51 * t24 + t75, t18 * t57 + t19 * t45 + t7 * t51 - t6 * t53, t11, t8, t15, t14, 0, t84 * t102 + t73 * t66 - t5 (-t84 * qJD(6) + t6) * t66 + t73 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -0.2e1 * t88, t52, 0, t42, 0.2e1 * t55 * t57 + 0.2e1 * t110, t42, 0.2e1 * t57, 0, 0.2e1 * t51 * t57 + 0.2e1 * t110, t44, t37, 0, 0, 0, -0.2e1 * t50 * t58 + 0.2e1 * t41, -0.2e1 * t50 * t102 - 0.2e1 * t45 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, 0, -t13, t12, -t13, pkin(3) * t24 + t81, -t12, -t13 * pkin(3) - t12 * qJ(4) + t27 * qJD(4), -t6, t7, t71 * t24 - t81, -t6 * qJ(4) + t19 * qJD(4) + t7 * t71, t11, t8, t15, t14, 0, t83 * t102 + t74 * t66 - t5 (-t83 * qJD(6) + t6) * t66 + t74 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t88, -t57, 0, t46, -pkin(3) * t57 + t80, t46, t57, 0, t71 * t57 + t80, t44, t37, 0, 0, 0, t66 * t90 + t41 + t59 (-qJD(4) - t45) * t66 + t69 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t61, t72, 0, 0, t61, t44, t37, 0, 0, 0, -0.2e1 * t65 * t58 + 0.2e1 * t59, -0.2e1 * qJD(4) * t66 - 0.2e1 * t65 * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, t13, 0, 0, t24, t7, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t25, 0, t4, 0, 0, 0, 0, 0, -t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t78, -t24, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, t58, 0, -t47 * t102 - t66 * t57, t47 * t58 - t69 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, t58, 0, -t62 * t102, t62 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t20;
