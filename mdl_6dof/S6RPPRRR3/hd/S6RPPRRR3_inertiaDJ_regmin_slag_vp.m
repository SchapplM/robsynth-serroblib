% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:48
% EndTime: 2019-03-09 02:23:51
% DurationCPUTime: 0.98s
% Computational Cost: add. (822->143), mult. (1961->257), div. (0->0), fcn. (1690->8), ass. (0->100)
t116 = qJD(5) + qJD(6);
t60 = sin(qJ(6));
t61 = sin(qJ(5));
t63 = cos(qJ(6));
t64 = cos(qJ(5));
t38 = t60 * t61 - t63 * t64;
t18 = t116 * t38;
t62 = sin(qJ(4));
t120 = t18 * t62;
t39 = t60 * t64 + t63 * t61;
t65 = cos(qJ(4));
t119 = t39 * t65;
t51 = sin(pkin(10)) * pkin(1) + qJ(3);
t72 = t62 * pkin(4) - t65 * pkin(8);
t34 = t51 + t72;
t23 = t61 * t34;
t108 = t62 * t64;
t50 = -cos(pkin(10)) * pkin(1) - pkin(2) - pkin(7);
t37 = t50 * t108;
t118 = -t23 - t37;
t96 = t62 * qJD(4);
t83 = t64 * t96;
t98 = qJD(5) * t65;
t88 = t61 * t98;
t31 = t83 + t88;
t19 = t116 * t39;
t54 = t65 * qJD(4);
t117 = t62 * t19 + t38 * t54;
t56 = t62 ^ 2;
t58 = t65 ^ 2;
t75 = (t56 - t58) * qJD(4);
t57 = t64 ^ 2;
t103 = t61 ^ 2 - t57;
t76 = t103 * qJD(5);
t115 = 2 * qJD(3);
t114 = pkin(8) + pkin(9);
t86 = t61 * t96;
t87 = t64 * t98;
t33 = t86 - t87;
t73 = pkin(4) * t65 + pkin(8) * t62;
t36 = t73 * qJD(4) + qJD(3);
t81 = t50 * t54;
t100 = qJD(5) * t61;
t89 = t62 * t100;
t99 = qJD(5) * t64;
t6 = -t34 * t99 - t61 * t36 + t50 * t89 - t64 * t81;
t5 = t33 * pkin(9) - t6;
t113 = t63 * t5;
t112 = t50 * t61;
t110 = t61 * t65;
t17 = -pkin(9) * t110 - t118;
t111 = t60 * t17;
t107 = t63 * t17;
t106 = t64 * t65;
t105 = t65 * t19;
t104 = t65 * t50;
t101 = t56 + t58;
t97 = qJD(6) * t60;
t95 = -0.2e1 * pkin(4) * qJD(5);
t94 = pkin(5) * t100;
t93 = pkin(5) * t54;
t92 = pkin(5) * t97;
t91 = qJD(6) * t63 * pkin(5);
t85 = t61 * t99;
t84 = t50 * t96;
t82 = t62 * t54;
t29 = t64 * t36;
t78 = pkin(5) - t112;
t4 = t29 + (-t37 + (pkin(9) * t65 - t34) * t61) * qJD(5) + (pkin(9) * t108 + t78 * t65) * qJD(4);
t80 = t63 * t4 - t60 * t5;
t24 = t64 * t34;
t16 = -pkin(9) * t106 + t78 * t62 + t24;
t79 = -t62 * pkin(5) - t16;
t77 = qJD(5) * t114;
t74 = t61 * t83;
t71 = t63 * t16 - t111;
t70 = t60 * t16 + t107;
t47 = t114 * t61;
t48 = t114 * t64;
t69 = -t63 * t47 - t60 * t48;
t68 = -t60 * t47 + t63 * t48;
t26 = t38 * t65;
t53 = -t64 * pkin(5) - pkin(4);
t49 = 0.2e1 * t82;
t41 = t64 * t77;
t40 = t61 * t77;
t32 = t61 * t54 + t62 * t99;
t30 = -t64 * t54 + t89;
t27 = (pkin(5) * t61 - t50) * t65;
t20 = -t33 * pkin(5) + t84;
t15 = t39 * t54 - t120;
t14 = -t68 * qJD(6) + t60 * t40 - t63 * t41;
t13 = -t69 * qJD(6) + t63 * t40 + t60 * t41;
t12 = -t97 * t110 + (t116 * t106 - t86) * t63 - t31 * t60;
t11 = -qJD(4) * t119 + t120;
t10 = -t60 * t86 + t63 * t83 + t105;
t7 = t118 * qJD(5) - t61 * t81 + t29;
t2 = -t70 * qJD(6) + t80;
t1 = -t71 * qJD(6) - t60 * t4 - t113;
t3 = [0, 0, 0, 0, 0, t115, t51 * t115, -0.2e1 * t82, 0.2e1 * t75, 0, 0, 0, 0.2e1 * qJD(3) * t62 + 0.2e1 * t51 * t54, 0.2e1 * qJD(3) * t65 - 0.2e1 * t51 * t96, -0.2e1 * t57 * t82 - 0.2e1 * t58 * t85, 0.2e1 * t58 * t76 + 0.4e1 * t65 * t74, -0.2e1 * t62 * t88 - 0.2e1 * t64 * t75, 0.2e1 * t61 * t75 - 0.2e1 * t62 * t87, t49, 0.2e1 * t24 * t54 + 0.2e1 * t7 * t62 + 0.2e1 * (-t58 * t99 + t61 * t82) * t50, 0.2e1 * t58 * t50 * t100 + 0.2e1 * t6 * t62 + 0.2e1 * (-t23 + t37) * t54, 0.2e1 * t26 * t10, 0.2e1 * t10 * t119 + 0.2e1 * t26 * t12, -0.2e1 * t10 * t62 - 0.2e1 * t26 * t54, -0.2e1 * t119 * t54 - 0.2e1 * t62 * t12, t49, 0.2e1 * t119 * t20 + 0.2e1 * t27 * t12 + 0.2e1 * t2 * t62 + 0.2e1 * t71 * t54, 0.2e1 * t1 * t62 - 0.2e1 * t27 * t10 - 0.2e1 * t20 * t26 - 0.2e1 * t70 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 * t99, t101 * t100, 0, 0, 0, 0, 0, t11 * t62 - t65 * t12, t65 * t10 + t117 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t54, 0, -t84, -t81, -t65 * t76 - t74, t103 * t96 - 0.4e1 * t65 * t85, t32, -t30, 0 (-t61 * t104 - t73 * t64) * qJD(5) + (t72 * t61 - t37) * qJD(4) (-t64 * t104 + t73 * t61) * qJD(5) + (-pkin(8) * t106 + (pkin(4) * t64 + t112) * t62) * qJD(4), -t10 * t39 + t26 * t18, t10 * t38 + t119 * t18 - t39 * t12 + t26 * t19, t15, -t117, 0, t119 * t94 + t53 * t12 + t14 * t62 + t27 * t19 + t20 * t38 + t69 * t54, -t53 * t10 + t13 * t62 - t27 * t18 + t20 * t39 - t26 * t94 - t68 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t96, 0, 0, 0, 0, 0, t30, t32, 0, 0, 0, 0, 0, t117, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t54, 0, 0, 0, 0, 0, -t31, t33, 0, 0, 0, 0, 0, t38 * t96 - t105, t65 * t18 + t39 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t85, -0.2e1 * t76, 0, 0, 0, t61 * t95, t64 * t95, -0.2e1 * t39 * t18, 0.2e1 * t18 * t38 - 0.2e1 * t39 * t19, 0, 0, 0, 0.2e1 * t53 * t19 + 0.2e1 * t38 * t94, -0.2e1 * t53 * t18 + 0.2e1 * t39 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t33, t54, t7, t6, 0, 0, -t10, -t12, t54, t63 * t93 + (t79 * t60 - t107) * qJD(6) + t80, -t113 + (-t4 - t93) * t60 + (t79 * t63 + t111) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t31, 0, 0, 0, 0, 0, -t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t30, 0, 0, 0, 0, 0, t11, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, -t100, 0, -pkin(8) * t99, pkin(8) * t100, 0, 0, -t18, -t19, 0, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t92, -0.2e1 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t12, t54, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
