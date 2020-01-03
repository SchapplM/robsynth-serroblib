% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:16
% EndTime: 2019-12-31 21:17:20
% DurationCPUTime: 1.12s
% Computational Cost: add. (1608->143), mult. (3869->273), div. (0->0), fcn. (3544->8), ass. (0->105)
t84 = sin(pkin(9));
t85 = cos(pkin(9));
t118 = t84 ^ 2 + t85 ^ 2;
t137 = t118 * qJD(4);
t139 = 0.2e1 * t137;
t90 = cos(qJ(3));
t116 = qJD(3) * t90;
t109 = pkin(2) * t116;
t72 = qJD(4) + t109;
t138 = t118 * t72;
t86 = sin(qJ(5));
t89 = cos(qJ(5));
t136 = -t86 * t84 + t89 * t85;
t135 = qJD(2) + qJD(3);
t134 = -pkin(7) - pkin(6);
t81 = t85 * pkin(8);
t133 = t90 * pkin(2);
t108 = qJD(2) * t134;
t91 = cos(qJ(2));
t102 = t91 * t108;
t88 = sin(qJ(2));
t103 = t88 * t108;
t70 = t134 * t88;
t71 = t134 * t91;
t87 = sin(qJ(3));
t48 = t87 * t70 - t90 * t71;
t32 = t48 * qJD(3) - t90 * t102 + t87 * t103;
t132 = t32 * t85;
t63 = t87 * t91 + t90 * t88;
t131 = t63 * t84;
t52 = t136 * qJD(5);
t75 = -t85 * pkin(4) - pkin(3);
t130 = t75 * t52;
t61 = t89 * t84 + t86 * t85;
t53 = t61 * qJD(5);
t129 = t75 * t53;
t62 = t87 * t88 - t90 * t91;
t42 = t135 * t62;
t128 = t84 * t42;
t127 = t85 * t42;
t15 = -pkin(4) * t128 + t32;
t47 = -t90 * t70 - t87 * t71;
t35 = pkin(4) * t131 + t47;
t124 = -t136 * t15 + t35 * t53;
t123 = t15 * t61 + t35 * t52;
t114 = t88 * qJD(2);
t111 = pkin(2) * t114;
t43 = t135 * t63;
t19 = t43 * pkin(3) + t42 * qJ(4) - t63 * qJD(4) + t111;
t117 = qJD(3) * t87;
t31 = -t87 * t102 - t90 * t103 - t70 * t116 - t71 * t117;
t9 = t84 * t19 - t85 * t31;
t78 = -t91 * pkin(2) - pkin(1);
t40 = t62 * pkin(3) - t63 * qJ(4) + t78;
t25 = t84 * t40 + t85 * t48;
t110 = pkin(2) * t117;
t67 = t75 - t133;
t122 = -t110 * t136 + t67 * t53;
t121 = t61 * t110 + t67 * t52;
t115 = qJD(5) * t63;
t113 = t91 * qJD(2);
t112 = -0.2e1 * pkin(1) * qJD(2);
t8 = t85 * t19 + t84 * t31;
t24 = t85 * t40 - t84 * t48;
t105 = t84 * t110;
t104 = t85 * t110;
t4 = -t8 * t84 + t9 * t85;
t14 = t62 * pkin(4) - t63 * t81 + t24;
t18 = -pkin(8) * t131 + t25;
t101 = t89 * t14 - t86 * t18;
t100 = t86 * t14 + t89 * t18;
t99 = -t24 * t84 + t25 * t85;
t98 = t32 * t63 - t42 * t47;
t74 = t87 * pkin(2) + qJ(4);
t57 = (-pkin(8) - t74) * t84;
t58 = t85 * t74 + t81;
t97 = t89 * t57 - t86 * t58;
t96 = t86 * t57 + t89 * t58;
t68 = (-pkin(8) - qJ(4)) * t84;
t69 = t85 * qJ(4) + t81;
t95 = t89 * t68 - t86 * t69;
t94 = t86 * t68 + t89 * t69;
t93 = pkin(3) * t42 - qJ(4) * t43 - qJD(4) * t62;
t77 = -pkin(3) - t133;
t92 = t63 * t110 - t42 * t77 - t43 * t74 - t62 * t72;
t41 = 0.2e1 * t61 * t52;
t37 = t136 * t63;
t36 = t61 * t63;
t34 = -t61 * qJD(4) - t94 * qJD(5);
t33 = -qJD(4) * t136 - t95 * qJD(5);
t28 = t32 * t84;
t26 = 0.2e1 * t136 * t52 - 0.2e1 * t61 * t53;
t23 = -t96 * qJD(5) - t61 * t72;
t22 = -t97 * qJD(5) - t136 * t72;
t21 = t136 * t43 - t53 * t62;
t20 = t61 * t43 + t52 * t62;
t11 = t136 * t115 - t61 * t42;
t10 = -t61 * t115 - t136 * t42;
t7 = t10 * t61 + t37 * t52;
t6 = pkin(8) * t128 + t9;
t5 = t43 * pkin(4) + pkin(8) * t127 + t8;
t3 = t10 * t136 - t61 * t11 - t52 * t36 - t37 * t53;
t2 = -t100 * qJD(5) + t89 * t5 - t86 * t6;
t1 = -t101 * qJD(5) - t86 * t5 - t89 * t6;
t12 = [0, 0, 0, 0.2e1 * t88 * t113, 0.2e1 * (-t88 ^ 2 + t91 ^ 2) * qJD(2), 0, 0, 0, t88 * t112, t91 * t112, -0.2e1 * t63 * t42, 0.2e1 * t42 * t62 - 0.2e1 * t63 * t43, 0, 0, 0, 0.2e1 * t62 * t111 + 0.2e1 * t78 * t43, 0.2e1 * t63 * t111 - 0.2e1 * t78 * t42, 0.2e1 * t24 * t43 + 0.2e1 * t8 * t62 + 0.2e1 * t98 * t84, -0.2e1 * t25 * t43 - 0.2e1 * t9 * t62 + 0.2e1 * t98 * t85, 0.2e1 * (-t8 * t85 - t84 * t9) * t63 - 0.2e1 * (-t24 * t85 - t25 * t84) * t42, 0.2e1 * t24 * t8 + 0.2e1 * t25 * t9 + 0.2e1 * t47 * t32, 0.2e1 * t37 * t10, -0.2e1 * t10 * t36 - 0.2e1 * t37 * t11, 0.2e1 * t10 * t62 + 0.2e1 * t37 * t43, -0.2e1 * t11 * t62 - 0.2e1 * t36 * t43, 0.2e1 * t62 * t43, 0.2e1 * t101 * t43 + 0.2e1 * t35 * t11 + 0.2e1 * t15 * t36 + 0.2e1 * t2 * t62, 0.2e1 * t1 * t62 + 0.2e1 * t35 * t10 - 0.2e1 * t100 * t43 + 0.2e1 * t15 * t37; 0, 0, 0, 0, 0, t113, -t114, 0, -pkin(6) * t113, pkin(6) * t114, 0, 0, -t42, -t43, 0, -t32, t31, t92 * t84 - t132, t92 * t85 + t28, t4, t47 * t110 + t32 * t77 + t4 * t74 + t99 * t72, t7, t3, t20, t21, 0, t67 * t11 + t36 * t110 + t23 * t62 + t97 * t43 + t124, t67 * t10 + t37 * t110 + t22 * t62 - t96 * t43 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t110, -0.2e1 * t109, -0.2e1 * t104, 0.2e1 * t105, 0.2e1 * t138, 0.2e1 * t77 * t110 + 0.2e1 * t138 * t74, t41, t26, 0, 0, 0, 0.2e1 * t122, 0.2e1 * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t43, 0, -t32, t31, t93 * t84 - t132, t93 * t85 + t28, t4, -t32 * pkin(3) + t4 * qJ(4) + t99 * qJD(4), t7, t3, t20, t21, 0, t75 * t11 + t34 * t62 + t95 * t43 + t124, t75 * t10 + t33 * t62 - t94 * t43 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t109, -t104, t105, t137 + t138, -pkin(3) * t110 + qJ(4) * t138 + t137 * t74, t41, t26, 0, 0, 0, t122 + t129, t121 + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, qJ(4) * t139, t41, t26, 0, 0, 0, 0.2e1 * t129, 0.2e1 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -t127, 0, t32, 0, 0, 0, 0, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, 0, 0, 0, 0, 0, t53, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t43, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t53, 0, t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t53, 0, t34, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
