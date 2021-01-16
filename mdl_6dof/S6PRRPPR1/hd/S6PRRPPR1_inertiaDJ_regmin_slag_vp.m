% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:06:04
% EndTime: 2021-01-16 02:06:12
% DurationCPUTime: 1.59s
% Computational Cost: add. (1549->173), mult. (4089->347), div. (0->0), fcn. (4022->12), ass. (0->110)
t72 = sin(pkin(12));
t75 = cos(pkin(12));
t76 = sin(qJ(6));
t79 = cos(qJ(6));
t132 = -t76 * t72 + t79 * t75;
t58 = t79 * t72 + t76 * t75;
t51 = t58 * qJD(6);
t131 = t132 * qJD(6);
t130 = 0.2e1 * t131;
t73 = sin(pkin(11));
t64 = t73 * pkin(3) + qJ(5);
t129 = pkin(9) + t64;
t115 = cos(pkin(11));
t81 = cos(qJ(2));
t113 = qJD(2) * t81;
t74 = sin(pkin(6));
t104 = t74 * t113;
t80 = cos(qJ(3));
t116 = cos(pkin(6));
t78 = sin(qJ(2));
t123 = t74 * t78;
t77 = sin(qJ(3));
t85 = t116 * t80 - t77 * t123;
t41 = t85 * qJD(3) + t80 * t104;
t52 = t116 * t77 + t80 * t123;
t82 = -t52 * qJD(3) - t77 * t104;
t17 = -t115 * t82 + t73 * t41;
t34 = -t115 * t85 + t73 * t52;
t128 = t34 * t17;
t118 = qJ(4) + pkin(8);
t101 = qJD(3) * t118;
t47 = t80 * qJD(4) - t77 * t101;
t83 = -t77 * qJD(4) - t80 * t101;
t30 = -t115 * t83 + t73 * t47;
t103 = t115 * t77;
t61 = t118 * t80;
t42 = t118 * t103 + t73 * t61;
t127 = t42 * t30;
t57 = t73 * t80 + t103;
t126 = t57 * t72;
t102 = t115 * t80;
t112 = qJD(3) * t77;
t49 = qJD(3) * t102 - t73 * t112;
t125 = t72 * t49;
t124 = t73 * t77;
t122 = t74 * t81;
t121 = t75 * t49;
t48 = t57 * qJD(3);
t69 = pkin(3) * t112;
t23 = t48 * pkin(4) - t49 * qJ(5) - t57 * qJD(5) + t69;
t31 = t115 * t47 + t73 * t83;
t8 = t72 * t23 + t75 * t31;
t55 = -t102 + t124;
t68 = -t80 * pkin(3) - pkin(2);
t38 = t55 * pkin(4) - t57 * qJ(5) + t68;
t43 = t115 * t61 - t118 * t124;
t16 = t72 * t38 + t75 * t43;
t117 = t72 ^ 2 + t75 ^ 2;
t114 = qJD(2) * t78;
t111 = qJD(3) * t80;
t110 = qJD(3) * t81;
t109 = -0.2e1 * pkin(2) * qJD(3);
t108 = t77 * t110;
t105 = t74 * t114;
t7 = t75 * t23 - t72 * t31;
t15 = t75 * t38 - t72 * t43;
t100 = 0.2e1 * t117 * qJD(5);
t99 = t7 * t75 + t8 * t72;
t98 = -t7 * t72 + t8 * t75;
t67 = -t115 * pkin(3) - pkin(4);
t10 = -pkin(9) * t126 + t16;
t9 = -t75 * t57 * pkin(9) + t55 * pkin(5) + t15;
t97 = t79 * t10 + t76 * t9;
t96 = t76 * t10 - t79 * t9;
t18 = t115 * t41 + t73 * t82;
t13 = t75 * t105 - t72 * t18;
t14 = t72 * t105 + t75 * t18;
t95 = t13 * t75 + t14 * t72;
t94 = -t13 * t72 + t14 * t75;
t93 = t17 * t42 + t34 * t30;
t92 = t17 * t57 + t34 * t49;
t35 = t115 * t52 + t73 * t85;
t26 = -t75 * t122 - t72 * t35;
t27 = -t72 * t122 + t75 * t35;
t91 = t79 * t26 - t76 * t27;
t90 = t76 * t26 + t79 * t27;
t89 = t30 * t57 + t42 * t49;
t88 = t131 * t55 + t58 * t48;
t53 = t129 * t72;
t54 = t129 * t75;
t87 = -t79 * t53 - t76 * t54;
t86 = -t76 * t53 + t79 * t54;
t84 = -qJD(5) * t55 - t48 * t64 + t49 * t67;
t60 = -t75 * pkin(5) + t67;
t33 = t132 * t57;
t32 = t58 * t57;
t29 = pkin(5) * t126 + t42;
t25 = -t58 * qJD(5) - t86 * qJD(6);
t24 = -qJD(5) * t132 - t87 * qJD(6);
t22 = t132 * t48 - t51 * t55;
t19 = pkin(5) * t125 + t30;
t12 = t131 * t57 + t58 * t49;
t11 = t132 * t49 - t57 * t51;
t6 = -pkin(9) * t125 + t8;
t5 = t48 * pkin(5) - pkin(9) * t121 + t7;
t4 = -t90 * qJD(6) + t79 * t13 - t76 * t14;
t3 = -t91 * qJD(6) - t76 * t13 - t79 * t14;
t2 = -t97 * qJD(6) + t79 * t5 - t76 * t6;
t1 = t96 * qJD(6) - t76 * t5 - t79 * t6;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t74 ^ 2 * t78 * t113 + 0.2e1 * t35 * t18 + 0.2e1 * t128, 0, 0, 0, 0.2e1 * t26 * t13 + 0.2e1 * t27 * t14 + 0.2e1 * t128, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t105, -t104, 0, 0, 0, 0, 0, (-t80 * t114 - t108) * t74, (-t80 * t110 + t77 * t114) * t74, (t55 * t114 - t48 * t81) * t74, (t57 * t114 - t49 * t81) * t74, -t18 * t55 - t35 * t48 + t92, t18 * t43 + t35 * t31 + (-pkin(3) * t108 + t68 * t114) * t74 + t93, t13 * t55 + t26 * t48 + t92 * t72, -t14 * t55 - t27 * t48 + t92 * t75, -t95 * t57 + (-t26 * t75 - t27 * t72) * t49, t13 * t15 + t14 * t16 + t26 * t7 + t27 * t8 + t93, 0, 0, 0, 0, 0, t34 * t12 + t17 * t32 + t4 * t55 + t48 * t91, t34 * t11 + t17 * t33 + t3 * t55 - t48 * t90; 0, 0, 0, 0, 0.2e1 * t77 * t111, 0.2e1 * (-t77 ^ 2 + t80 ^ 2) * qJD(3), 0, 0, 0, t77 * t109, t80 * t109, 0.2e1 * t68 * t48 + 0.2e1 * t55 * t69, 0.2e1 * t68 * t49 + 0.2e1 * t57 * t69, -0.2e1 * t31 * t55 - 0.2e1 * t43 * t48 + 0.2e1 * t89, 0.2e1 * t43 * t31 + 0.2e1 * t68 * t69 + 0.2e1 * t127, 0.2e1 * t15 * t48 + 0.2e1 * t7 * t55 + 0.2e1 * t89 * t72, -0.2e1 * t16 * t48 - 0.2e1 * t8 * t55 + 0.2e1 * t89 * t75, -0.2e1 * t99 * t57 + 0.2e1 * (-t15 * t75 - t16 * t72) * t49, 0.2e1 * t15 * t7 + 0.2e1 * t16 * t8 + 0.2e1 * t127, 0.2e1 * t33 * t11, -0.2e1 * t11 * t32 - 0.2e1 * t33 * t12, 0.2e1 * t11 * t55 + 0.2e1 * t33 * t48, -0.2e1 * t12 * t55 - 0.2e1 * t32 * t48, 0.2e1 * t55 * t48, 0.2e1 * t29 * t12 + 0.2e1 * t19 * t32 + 0.2e1 * t2 * t55 - 0.2e1 * t48 * t96, 0.2e1 * t1 * t55 + 0.2e1 * t29 * t11 + 0.2e1 * t19 * t33 - 0.2e1 * t48 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t41, -t17, -t18, 0, (-t115 * t17 + t18 * t73) * pkin(3), -t17 * t75, t17 * t72, t94, t17 * t67 + t94 * t64 + (-t26 * t72 + t27 * t75) * qJD(5), 0, 0, 0, 0, 0, -t132 * t17 + t34 * t51, t131 * t34 + t17 * t58; 0, 0, 0, 0, 0, 0, t111, -t112, 0, -pkin(8) * t111, pkin(8) * t112, -t30, -t31, (-t115 * t49 - t48 * t73) * pkin(3), (-t115 * t30 + t31 * t73) * pkin(3), -t30 * t75 + t84 * t72, t30 * t72 + t84 * t75, t98, t30 * t67 + t98 * t64 + (-t15 * t72 + t16 * t75) * qJD(5), t11 * t58 + t131 * t33, t11 * t132 - t58 * t12 - t131 * t32 - t33 * t51, t88, t22, 0, t60 * t12 - t132 * t19 + t25 * t55 + t29 * t51 + t48 * t87, t60 * t11 + t131 * t29 + t19 * t58 + t24 * t55 - t48 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t64 * t100, t58 * t130, 0.2e1 * t131 * t132 - 0.2e1 * t58 * t51, 0, 0, 0, 0.2e1 * t60 * t51, t60 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t49, 0, t69, t75 * t48, -t72 * t48, -t117 * t49, t99, 0, 0, 0, 0, 0, t22, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t121, 0, t30, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t12, t48, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, -t51, 0, t25, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t20;
