% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:16:26
% EndTime: 2019-03-08 21:16:30
% DurationCPUTime: 1.37s
% Computational Cost: add. (883->204), mult. (2552->369), div. (0->0), fcn. (2310->10), ass. (0->119)
t76 = sin(pkin(11));
t74 = t76 ^ 2;
t78 = cos(pkin(11));
t144 = (t78 ^ 2 + t74) * qJD(4);
t80 = sin(qJ(3));
t122 = qJD(3) * t80;
t83 = cos(qJ(3));
t143 = qJ(5) * t122 - t83 * qJD(5);
t125 = t80 * qJ(4);
t102 = t76 * qJ(5) + pkin(3);
t52 = -t78 * pkin(4) - t102;
t142 = (-t52 * t83 + t125) * qJD(3);
t82 = cos(qJ(6));
t117 = qJD(6) * t82;
t79 = sin(qJ(6));
t118 = qJD(6) * t79;
t44 = t76 * t117 - t78 * t118;
t51 = t82 * t76 - t79 * t78;
t141 = pkin(4) + pkin(5);
t81 = sin(qJ(2));
t123 = qJD(2) * t81;
t77 = sin(pkin(6));
t106 = t77 * t123;
t84 = cos(qJ(2));
t134 = t77 * t84;
t105 = qJD(2) * t134;
t124 = cos(pkin(6));
t135 = t77 * t81;
t45 = -t124 * t83 + t80 * t135;
t30 = -t45 * qJD(3) + t83 * t105;
t15 = t76 * t106 + t30 * t78;
t10 = t15 * t78;
t46 = t124 * t80 + t83 * t135;
t29 = -t76 * t134 + t46 * t78;
t140 = t29 * t78 * qJD(4) + qJ(4) * t10;
t31 = t46 * qJD(3) + t80 * t105;
t139 = t31 * t76;
t138 = t31 * t78;
t137 = t76 * t80;
t136 = t76 * t83;
t42 = -t80 * qJD(4) + (pkin(3) * t80 - qJ(4) * t83) * qJD(3);
t133 = t78 * t42;
t132 = t78 * t80;
t131 = t78 * t83;
t128 = -pkin(9) + qJ(4);
t99 = -t83 * pkin(3) - t125;
t54 = -pkin(2) + t99;
t67 = pkin(8) * t131;
t37 = t76 * t54 + t67;
t127 = qJ(4) * t144;
t126 = qJ(5) * t78;
t121 = qJD(3) * t83;
t120 = qJD(3) * t84;
t119 = qJD(4) * t83;
t116 = qJD(6) * t83;
t115 = t74 * qJD(5);
t114 = t76 * qJD(4);
t113 = t76 * qJD(5);
t66 = pkin(8) * t136;
t110 = -0.2e1 * pkin(2) * qJD(3);
t109 = pkin(8) * t122;
t72 = pkin(8) * t121;
t63 = t76 * t121;
t64 = t78 * t121;
t104 = t80 * t121;
t103 = -pkin(8) * t76 - pkin(4);
t14 = -t78 * t106 + t30 * t76;
t101 = t14 * t76 + t10;
t36 = t78 * t54 - t66;
t33 = -t83 * qJ(5) + t37;
t100 = -qJD(5) * t132 + t72;
t98 = pkin(4) * t76 - t126;
t73 = t83 * pkin(4);
t20 = t83 * pkin(5) + t66 + t73 + (-pkin(9) * t80 - t54) * t78;
t24 = pkin(9) * t137 + t33;
t97 = t82 * t20 - t79 * t24;
t96 = t79 * t20 + t82 * t24;
t26 = t76 * t109 + t133;
t38 = t76 * t42;
t27 = -t78 * t109 + t38;
t95 = -t26 * t76 + t27 * t78;
t28 = t78 * t134 + t46 * t76;
t94 = t28 * t82 - t29 * t79;
t93 = t28 * t79 + t29 * t82;
t59 = t128 * t76;
t60 = t128 * t78;
t92 = t82 * t59 - t79 * t60;
t91 = t79 * t59 + t82 * t60;
t50 = t79 * t76 + t82 * t78;
t89 = qJ(4) * t14 + qJD(4) * t28;
t88 = -t141 * t76 + t126;
t87 = 0.2e1 * t28 * t14 + 0.2e1 * t29 * t15 + 0.2e1 * t45 * t31;
t43 = t50 * qJD(6);
t86 = -t28 * t122 + t31 * t137 + t14 * t83 + t45 * t63;
t85 = (t14 * t78 - t15 * t76) * t80 + (t28 * t78 - t29 * t76) * t121;
t62 = t83 * t114;
t49 = 0.2e1 * t144;
t47 = t141 * t78 + t102;
t41 = t50 * t80;
t40 = t51 * t80;
t39 = (pkin(8) + t98) * t80;
t34 = -t36 + t73;
t32 = (-pkin(8) + t88) * t80;
t25 = t98 * t121 + t100;
t21 = t103 * t122 - t133;
t19 = t88 * t121 - t100;
t18 = t27 + t143;
t17 = t80 * t43 - t82 * t63 + t79 * t64;
t16 = -t44 * t80 - t79 * t63 - t82 * t64;
t12 = t51 * qJD(4) - t91 * qJD(6);
t11 = -t50 * qJD(4) - t92 * qJD(6);
t8 = t38 + (-pkin(8) * t132 + pkin(9) * t136) * qJD(3) + t143;
t7 = -t133 + (-pkin(9) * t131 + (-pkin(5) + t103) * t80) * qJD(3);
t5 = t31 * t132 + t15 * t83 + (t45 * t131 - t29 * t80) * qJD(3);
t4 = t94 * qJD(6) + t14 * t79 + t15 * t82;
t3 = -t93 * qJD(6) + t14 * t82 - t15 * t79;
t2 = -t96 * qJD(6) + t82 * t7 - t79 * t8;
t1 = -t97 * qJD(6) - t79 * t7 - t82 * t8;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, 0, 0, t87, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t106, -t105, 0, 0, 0, 0, 0 (-t80 * t120 - t83 * t123) * t77 (-t83 * t120 + t80 * t123) * t77, t86, t5, t85, -t14 * t36 + t15 * t37 - t28 * t26 + t29 * t27 + (t45 * t121 + t31 * t80) * pkin(8), t86, t85, -t5, t14 * t34 + t15 * t33 + t29 * t18 + t28 * t21 + t45 * t25 + t31 * t39, 0, 0, 0, 0, 0, -t94 * t122 - t45 * t17 + t3 * t83 + t31 * t40, t93 * t122 + t45 * t16 - t31 * t41 - t4 * t83; 0, 0, 0, 0, 0.2e1 * t104, 0.2e1 * (-t80 ^ 2 + t83 ^ 2) * qJD(3), 0, 0, 0, t80 * t110, t83 * t110, -0.2e1 * t26 * t83 + 0.2e1 * (t36 + 0.2e1 * t66) * t122, 0.2e1 * t27 * t83 + 0.2e1 * (-t37 + 0.2e1 * t67) * t122, 0.2e1 * (-t26 * t78 - t27 * t76) * t80 + 0.2e1 * (-t36 * t78 - t37 * t76) * t121, 0.2e1 * pkin(8) ^ 2 * t104 + 0.2e1 * t36 * t26 + 0.2e1 * t37 * t27, 0.2e1 * t25 * t137 + 0.2e1 * t21 * t83 + 0.2e1 * (t39 * t136 - t34 * t80) * qJD(3), 0.2e1 * (-t18 * t76 + t21 * t78) * t80 + 0.2e1 * (-t33 * t76 + t34 * t78) * t121, -0.2e1 * t25 * t132 - 0.2e1 * t18 * t83 + 0.2e1 * (-t39 * t131 + t33 * t80) * qJD(3), 0.2e1 * t33 * t18 + 0.2e1 * t34 * t21 + 0.2e1 * t39 * t25, -0.2e1 * t41 * t16, -0.2e1 * t16 * t40 - 0.2e1 * t41 * t17, -0.2e1 * t41 * t122 - 0.2e1 * t16 * t83, -0.2e1 * t40 * t122 - 0.2e1 * t17 * t83, -0.2e1 * t104, -0.2e1 * t97 * t122 + 0.2e1 * t32 * t17 - 0.2e1 * t19 * t40 + 0.2e1 * t2 * t83, 0.2e1 * t1 * t83 + 0.2e1 * t96 * t122 - 0.2e1 * t32 * t16 + 0.2e1 * t19 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t30, -t138, t139, t101, -t31 * pkin(3) + t89 * t76 + t140, -t138, t101, -t139, t31 * t52 + (-qJD(5) * t45 + t89) * t76 + t140, 0, 0, 0, 0, 0, -t31 * t50 - t45 * t44, -t31 * t51 + t45 * t43; 0, 0, 0, 0, 0, 0, t121, -t122, 0, -t72, t109, t62 + (t99 * t76 - t67) * qJD(3), t78 * t119 + (t99 * t78 + t66) * qJD(3), t95, -pkin(3) * t72 + (-t36 * t76 + t37 * t78) * qJD(4) + t95 * qJ(4), -t80 * t115 - t76 * t142 - t25 * t78 + t62, t18 * t78 + t21 * t76, -t25 * t76 + (t80 * t113 - t119 + t142) * t78, t25 * t52 + (qJ(4) * t18 + qJD(4) * t33) * t78 + (qJ(4) * t21 + qJD(4) * t34 - qJD(5) * t39) * t76, -t16 * t51 - t41 * t43, t16 * t50 - t51 * t17 - t43 * t40 - t41 * t44, -t51 * t122 - t43 * t83, t50 * t122 - t44 * t83, 0, -t40 * t113 + t12 * t83 - t92 * t122 + t47 * t17 + t19 * t50 + t32 * t44, t11 * t83 + t41 * t113 + t91 * t122 - t47 * t16 + t19 * t51 - t32 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0.2e1 * t127, 0.2e1 * t78 * t113, t49, 0.2e1 * t115, -0.2e1 * t52 * t113 + 0.2e1 * t127, -0.2e1 * t51 * t43, 0.2e1 * t43 * t50 - 0.2e1 * t51 * t44, 0, 0, 0, 0.2e1 * t50 * t113 + 0.2e1 * t47 * t44, 0.2e1 * t51 * t113 - 0.2e1 * t47 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t64, 0, t72, t63, 0, -t64, t25, 0, 0, 0, 0, 0, -t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, 0, 0, 0, 0, 0, -t44, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t64, 0, t21, 0, 0, 0, 0, 0, -t79 * t116 - t82 * t122, -t82 * t116 + t79 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, -t122, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, -t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
