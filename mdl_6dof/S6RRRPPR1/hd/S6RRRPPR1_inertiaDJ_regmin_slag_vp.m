% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:23:06
% EndTime: 2019-03-09 15:23:12
% DurationCPUTime: 1.87s
% Computational Cost: add. (4336->188), mult. (9881->353), div. (0->0), fcn. (9591->10), ass. (0->135)
t104 = sin(pkin(11));
t106 = cos(pkin(11));
t143 = t104 ^ 2 + t106 ^ 2;
t169 = qJD(5) * t143;
t170 = 0.2e1 * t169;
t105 = sin(pkin(10));
t111 = cos(qJ(3));
t155 = pkin(2) * qJD(3);
t137 = t111 * t155;
t108 = sin(qJ(3));
t138 = t108 * t155;
t146 = cos(pkin(10));
t72 = -t105 * t138 + t146 * t137;
t68 = qJD(5) + t72;
t168 = t68 * t143;
t107 = sin(qJ(6));
t110 = cos(qJ(6));
t144 = t110 * t104;
t81 = t107 * t106 + t144;
t77 = t81 * qJD(6);
t109 = sin(qJ(2));
t167 = pkin(7) + pkin(8);
t87 = t167 * t109;
t112 = cos(qJ(2));
t88 = t167 * t112;
t123 = t108 * t88 + t111 * t87;
t82 = t108 * t112 + t111 * t109;
t117 = -t82 * qJ(4) - t123;
t166 = t106 * pkin(5);
t101 = t106 * pkin(9);
t122 = t108 * t87 - t111 * t88;
t136 = qJD(2) * t167;
t83 = t109 * t136;
t84 = t112 * t136;
t47 = t122 * qJD(3) + t108 * t83 - t111 * t84;
t121 = t108 * t109 - t111 * t112;
t62 = (-qJD(2) - qJD(3)) * t121;
t113 = -t62 * qJ(4) - t82 * qJD(4) + t47;
t119 = t82 * qJD(2);
t124 = t108 * t84 + t111 * t83;
t34 = -qJ(4) * t119 + t117 * qJD(3) - t121 * qJD(4) - t124;
t19 = t105 * t34 - t146 * t113;
t52 = -t121 * qJ(4) - t122;
t39 = t105 * t52 - t146 * t117;
t165 = t39 * t19;
t133 = t146 * t108;
t71 = (t105 * t111 + t133) * t155;
t164 = t39 * t71;
t58 = -t105 * t121 + t146 * t82;
t163 = t71 * t58;
t145 = t107 * t104;
t135 = qJD(6) * t145;
t140 = qJD(6) * t110;
t76 = -t106 * t140 + t135;
t95 = -t146 * pkin(3) - pkin(4);
t85 = t95 - t166;
t162 = t85 * t76;
t161 = t85 * t77;
t114 = -t82 * qJD(3) - t119;
t44 = t105 * t114 + t146 * t62;
t154 = t104 * t44;
t12 = pkin(5) * t154 + t19;
t153 = t104 * t58;
t30 = pkin(5) * t153 + t39;
t80 = -t110 * t106 + t145;
t160 = t12 * t80 + t30 * t77;
t159 = t12 * t81 - t30 * t76;
t20 = t105 * t113 + t146 * t34;
t43 = t105 * t62 - t146 * t114;
t142 = qJD(2) * t109;
t98 = pkin(2) * t142;
t55 = -t114 * pkin(3) + t98;
t23 = t43 * pkin(4) - t44 * qJ(5) - t58 * qJD(5) + t55;
t8 = t104 * t23 + t106 * t20;
t97 = -t112 * pkin(2) - pkin(1);
t115 = t121 * pkin(3) + t97;
t57 = t105 * t82 + t146 * t121;
t38 = t57 * pkin(4) - t58 * qJ(5) + t115;
t40 = t105 * t117 + t146 * t52;
t25 = t104 * t38 + t106 * t40;
t96 = t111 * pkin(2) + pkin(3);
t73 = -t105 * t108 * pkin(2) + t146 * t96;
t70 = -pkin(4) - t73;
t67 = t70 - t166;
t158 = t67 * t77 + t71 * t80;
t157 = -t67 * t76 + t71 * t81;
t74 = pkin(2) * t133 + t105 * t96;
t152 = t106 * t44;
t151 = t19 * t106;
t150 = t71 * t104;
t149 = t71 * t106;
t141 = qJD(2) * t112;
t139 = -0.2e1 * pkin(1) * qJD(2);
t7 = -t104 * t20 + t106 * t23;
t24 = -t104 * t40 + t106 * t38;
t3 = -t7 * t104 + t8 * t106;
t131 = t19 * t58 + t39 * t44;
t28 = t81 * t43 - t76 * t57;
t13 = t57 * pkin(5) - t58 * t101 + t24;
t16 = -pkin(9) * t153 + t25;
t130 = t107 * t16 - t110 * t13;
t129 = t107 * t13 + t110 * t16;
t69 = qJ(5) + t74;
t63 = (-pkin(9) - t69) * t104;
t64 = t106 * t69 + t101;
t128 = t107 * t64 - t110 * t63;
t127 = t107 * t63 + t110 * t64;
t93 = t105 * pkin(3) + qJ(5);
t78 = (-pkin(9) - t93) * t104;
t79 = t106 * t93 + t101;
t126 = t107 * t79 - t110 * t78;
t125 = t107 * t78 + t110 * t79;
t120 = t97 * t82;
t118 = -qJD(5) * t57 - t43 * t93 + t44 * t95;
t116 = -t43 * t69 + t44 * t70 - t57 * t68 + t163;
t61 = -0.2e1 * t81 * t76;
t49 = -t81 * qJD(5) - t125 * qJD(6);
t48 = t80 * qJD(5) + t126 * qJD(6);
t46 = t123 * qJD(3) + t124;
t45 = 0.2e1 * t76 * t80 - 0.2e1 * t81 * t77;
t42 = t80 * t58;
t41 = t81 * t58;
t32 = -t127 * qJD(6) - t81 * t68;
t31 = t128 * qJD(6) + t80 * t68;
t29 = -t80 * t43 - t77 * t57;
t18 = t19 * t104;
t15 = t44 * t144 - t58 * t135 + (t107 * t44 + t58 * t140) * t106;
t14 = -t80 * t44 - t58 * t77;
t9 = t14 * t81 + t42 * t76;
t6 = -pkin(9) * t154 + t8;
t5 = t43 * pkin(5) - pkin(9) * t152 + t7;
t4 = -t14 * t80 - t81 * t15 + t76 * t41 + t42 * t77;
t2 = -t129 * qJD(6) - t107 * t6 + t110 * t5;
t1 = t130 * qJD(6) - t107 * t5 - t110 * t6;
t10 = [0, 0, 0, 0.2e1 * t109 * t141, 0.2e1 * (-t109 ^ 2 + t112 ^ 2) * qJD(2), 0, 0, 0, t109 * t139, t112 * t139, 0.2e1 * t82 * t62, 0.2e1 * t82 * t114 - 0.2e1 * t62 * t121, 0, 0, 0, 0.2e1 * qJD(3) * t120 + 0.2e1 * (t109 * pkin(2) * t121 + t120) * qJD(2), 0.2e1 * t97 * t62 + 0.2e1 * t82 * t98, -0.2e1 * t20 * t57 - 0.2e1 * t40 * t43 + 0.2e1 * t131, 0.2e1 * t115 * t55 + 0.2e1 * t40 * t20 + 0.2e1 * t165, 0.2e1 * t131 * t104 + 0.2e1 * t24 * t43 + 0.2e1 * t7 * t57, 0.2e1 * t131 * t106 - 0.2e1 * t25 * t43 - 0.2e1 * t8 * t57, 0.2e1 * (-t24 * t44 - t58 * t7) * t106 + 0.2e1 * (-t25 * t44 - t58 * t8) * t104, 0.2e1 * t24 * t7 + 0.2e1 * t25 * t8 + 0.2e1 * t165, -0.2e1 * t42 * t14, -0.2e1 * t14 * t41 + 0.2e1 * t42 * t15, 0.2e1 * t14 * t57 - 0.2e1 * t42 * t43, -0.2e1 * t15 * t57 - 0.2e1 * t41 * t43, 0.2e1 * t57 * t43, 0.2e1 * t12 * t41 - 0.2e1 * t130 * t43 + 0.2e1 * t30 * t15 + 0.2e1 * t2 * t57, 0.2e1 * t1 * t57 - 0.2e1 * t12 * t42 - 0.2e1 * t129 * t43 + 0.2e1 * t30 * t14; 0, 0, 0, 0, 0, t141, -t142, 0, -pkin(7) * t141, pkin(7) * t142, 0, 0, t62, t114, 0, t47, t46, -t74 * t43 - t73 * t44 - t72 * t57 + t163, -t19 * t73 + t20 * t74 + t40 * t72 + t164, t116 * t104 - t151, t116 * t106 + t18, t3, t19 * t70 + t164 + (t25 * t68 + t69 * t8) * t106 + (-t24 * t68 - t69 * t7) * t104, t9, t4, t28, t29, 0, -t128 * t43 + t67 * t15 + t32 * t57 + t71 * t41 + t160, -t127 * t43 + t67 * t14 + t31 * t57 - t71 * t42 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t138, -0.2e1 * t137, 0, -0.2e1 * t73 * t71 + 0.2e1 * t74 * t72, -0.2e1 * t149, 0.2e1 * t150, 0.2e1 * t168, 0.2e1 * t168 * t69 + 0.2e1 * t70 * t71, t61, t45, 0, 0, 0, 0.2e1 * t158, 0.2e1 * t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t114, 0, t47, t46 (-t105 * t43 - t146 * t44) * pkin(3) (t105 * t20 - t146 * t19) * pkin(3), t118 * t104 - t151, t118 * t106 + t18, t3, t19 * t95 + t3 * t93 + (-t104 * t24 + t106 * t25) * qJD(5), t9, t4, t28, t29, 0, -t126 * t43 + t85 * t15 + t49 * t57 + t160, -t125 * t43 + t85 * t14 + t48 * t57 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, -t137, 0 (t105 * t72 - t146 * t71) * pkin(3), -t149, t150, t169 + t168, t168 * t93 + t169 * t69 + t71 * t95, t61, t45, 0, 0, 0, t158 + t161, t157 - t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t93 * t170, t61, t45, 0, 0, 0, 0.2e1 * t161, -0.2e1 * t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t106 * t43, -t104 * t43, -t143 * t44, t8 * t104 + t7 * t106, 0, 0, 0, 0, 0, t29, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t152, 0, t19, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, 0, 0, 0, t77, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, t43, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t77, 0, t32, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t77, 0, t49, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
