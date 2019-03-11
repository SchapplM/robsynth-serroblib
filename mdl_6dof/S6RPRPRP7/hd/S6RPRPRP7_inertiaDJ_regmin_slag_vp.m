% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:46
% EndTime: 2019-03-09 03:22:49
% DurationCPUTime: 0.80s
% Computational Cost: add. (1351->149), mult. (2671->260), div. (0->0), fcn. (2421->6), ass. (0->92)
t54 = sin(qJ(3));
t56 = cos(qJ(3));
t97 = sin(pkin(9));
t77 = qJD(3) * t97;
t98 = cos(pkin(9));
t78 = qJD(3) * t98;
t30 = -t54 * t78 - t56 * t77;
t81 = t98 * t56;
t35 = -t97 * t54 + t81;
t110 = t35 * t30;
t29 = t54 * t77 - t56 * t78;
t80 = t97 * t56;
t36 = -t98 * t54 - t80;
t111 = t29 * t36;
t61 = 0.2e1 * t110 + 0.2e1 * t111;
t116 = (-t97 * t29 + t98 * t30) * pkin(3);
t101 = t54 * pkin(3) + qJ(2);
t20 = -t36 * pkin(4) - t35 * pkin(8) + t101;
t57 = -pkin(1) - pkin(7);
t100 = qJ(4) - t57;
t38 = t100 * t54;
t23 = -t100 * t80 - t98 * t38;
t55 = cos(qJ(5));
t21 = t55 * t23;
t53 = sin(qJ(5));
t105 = t53 * t20 + t21;
t47 = t97 * pkin(3) + pkin(8);
t99 = qJ(6) + t47;
t76 = qJD(5) * t99;
t24 = t55 * qJD(6) - t53 * t76;
t25 = -t53 * qJD(6) - t55 * t76;
t31 = t99 * t53;
t32 = t99 * t55;
t115 = (-t31 * t55 + t32 * t53) * qJD(5) - t24 * t55 + t25 * t53;
t102 = qJ(6) * t35;
t66 = -qJ(6) * t30 - qJD(6) * t35;
t93 = t56 * qJD(3);
t28 = -t54 * qJD(4) - t100 * t93;
t94 = t54 * qJD(3);
t60 = -t56 * qJD(4) + t100 * t94;
t13 = t98 * t28 + t97 * t60;
t40 = pkin(3) * t93 + qJD(2);
t14 = -t29 * pkin(4) - t30 * pkin(8) + t40;
t83 = -t53 * t13 + t55 * t14;
t1 = -t29 * pkin(5) + t66 * t55 + (-t21 + (-t20 + t102) * t53) * qJD(5) + t83;
t95 = qJD(5) * t55;
t87 = t35 * t95;
t91 = t55 * t13 + t53 * t14 + t20 * t95;
t2 = -qJ(6) * t87 + (-qJD(5) * t23 + t66) * t53 + t91;
t82 = t55 * t20 - t53 * t23;
t5 = -t36 * pkin(5) - t55 * t102 + t82;
t6 = -t53 * t102 + t105;
t73 = t5 * t53 - t55 * t6;
t114 = t73 * qJD(5) - t1 * t55 - t2 * t53;
t34 = t35 ^ 2;
t113 = 0.2e1 * qJD(2);
t112 = 0.2e1 * qJD(5);
t109 = t35 * t53;
t108 = t35 * t55;
t107 = t53 * t30;
t106 = t55 * t30;
t51 = t53 ^ 2;
t52 = t55 ^ 2;
t104 = t51 - t52;
t103 = t51 + t52;
t96 = qJD(5) * t53;
t92 = qJ(2) * qJD(3);
t48 = -t98 * pkin(3) - pkin(4);
t90 = t48 * t112;
t89 = pkin(5) * t96;
t88 = t35 * t96;
t86 = t53 * t95;
t85 = t36 ^ 2 + t34;
t84 = -0.4e1 * t53 * t108;
t79 = t104 * qJD(5);
t74 = t5 * t55 + t53 * t6;
t70 = t29 * t47 + t30 * t48;
t68 = -t31 * t53 - t32 * t55;
t67 = t35 * t48 + t36 * t47;
t12 = t97 * t28 - t98 * t60;
t22 = t100 * t81 - t97 * t38;
t65 = t53 * t29 + t36 * t95;
t64 = -t55 * t29 + t36 * t96;
t63 = t87 + t107;
t62 = t88 - t106;
t58 = t12 * t35 + t13 * t36 + t22 * t30 + t23 * t29;
t39 = -t55 * pkin(5) + t48;
t9 = pkin(5) * t109 + t22;
t7 = t63 * pkin(5) + t12;
t4 = -t105 * qJD(5) + t83;
t3 = t23 * t96 - t91;
t8 = [0, 0, 0, 0, t113, qJ(2) * t113, -0.2e1 * t54 * t93, 0.2e1 * (t54 ^ 2 - t56 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t54 + 0.2e1 * t56 * t92, 0.2e1 * qJD(2) * t56 - 0.2e1 * t54 * t92, 0.2e1 * t58, 0.2e1 * t101 * t40 + 0.2e1 * t22 * t12 + 0.2e1 * t23 * t13, 0.2e1 * t52 * t110 - 0.2e1 * t34 * t86, t104 * t34 * t112 + t30 * t84, -0.2e1 * t36 * t106 + 0.2e1 * t64 * t35, 0.2e1 * t36 * t107 + 0.2e1 * t65 * t35, 0.2e1 * t111, 0.2e1 * t12 * t109 + 0.2e1 * t63 * t22 - 0.2e1 * t82 * t29 - 0.2e1 * t4 * t36, 0.2e1 * t105 * t29 + 0.2e1 * t12 * t108 - 0.2e1 * t62 * t22 - 0.2e1 * t3 * t36, 0.2e1 * t114 * t35 - 0.2e1 * t74 * t30, 0.2e1 * t5 * t1 + 0.2e1 * t6 * t2 + 0.2e1 * t9 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t58, 0, 0, 0, 0, 0, -t53 * t61 - t85 * t95, -t55 * t61 + t85 * t96, 0, -t9 * t30 - t7 * t35 + t73 * t29 + (t74 * qJD(5) + t1 * t53 - t2 * t55) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t103 * t111 + 0.2e1 * t110; 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t93, 0, -t57 * t94, -t57 * t93, -t116 (-t98 * t12 + t97 * t13) * pkin(3), t53 * t106 - t35 * t79, qJD(5) * t84 - t104 * t30, -t65, t64, 0, -t12 * t55 + t70 * t53 + (t22 * t53 + t67 * t55) * qJD(5), t12 * t53 + t70 * t55 + (t22 * t55 - t67 * t53) * qJD(5) (-t25 * t35 + t30 * t31 + t2 + (-t32 * t35 - t5) * qJD(5)) * t55 + (-t24 * t35 - t30 * t32 - t1 + (-t31 * t35 - t6) * qJD(5)) * t53, -t1 * t31 + t2 * t32 + t6 * t24 + t5 * t25 + t7 * t39 + t9 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t93, 0, t116, 0, 0, 0, 0, 0, -t62, -t63, -t103 * t29, -pkin(5) * t88 + t115 * t36 + t68 * t29 - t30 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t86, -0.2e1 * t79, 0, 0, 0, t53 * t90, t55 * t90, -0.2e1 * t115, 0.2e1 * t32 * t24 - 0.2e1 * t31 * t25 + 0.2e1 * t39 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, t64, t65, -t103 * t30, -t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68 * qJD(5) + t24 * t53 + t25 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t63, -t29, t4, t3, t62 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t64, 0, t65 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t96, 0, -t47 * t95, t47 * t96, -pkin(5) * t95, t25 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t95, 0, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
