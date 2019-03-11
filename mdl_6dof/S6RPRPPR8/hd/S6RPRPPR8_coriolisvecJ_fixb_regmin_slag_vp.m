% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:59:59
% EndTime: 2019-03-09 03:00:03
% DurationCPUTime: 1.31s
% Computational Cost: add. (1007->229), mult. (1979->314), div. (0->0), fcn. (964->4), ass. (0->141)
t68 = cos(qJ(3));
t121 = t68 * qJD(1);
t51 = qJD(6) + t121;
t160 = qJD(6) - t51;
t116 = qJ(5) * qJD(1);
t70 = -pkin(1) - pkin(7);
t47 = t70 * qJD(1) + qJD(2);
t26 = (t47 + t116) * t68;
t159 = qJD(4) - t26;
t69 = -pkin(3) - pkin(4);
t13 = qJD(3) * t69 + t159;
t140 = t68 * t47;
t66 = sin(qJ(3));
t38 = t66 * t47;
t61 = qJD(3) * qJ(4);
t31 = t38 + t61;
t150 = t31 * t68;
t28 = (qJD(4) + t140) * qJD(3);
t133 = qJD(3) * pkin(3);
t98 = -qJD(4) + t133;
t29 = -t98 - t140;
t158 = ((-t29 + t140) * t66 - t150) * qJD(3) - t28 * t66;
t111 = 0.2e1 * qJD(1);
t114 = qJD(1) * qJD(5);
t115 = qJD(1) * qJD(3);
t107 = t66 * t115;
t129 = qJD(3) * t66;
t36 = t47 * t129;
t137 = qJ(5) * t107 + t36;
t12 = -t68 * t114 + t137;
t54 = t66 * t116;
t21 = -t54 - t31;
t15 = qJD(3) * pkin(5) - t21;
t58 = t68 * qJ(4);
t60 = -pkin(8) + t69;
t79 = t68 * pkin(5) + t60 * t66 - qJ(2);
t19 = t58 + t79;
t131 = qJ(5) + t70;
t99 = qJD(3) * t131;
t22 = -t68 * qJD(5) + t66 * t99;
t55 = qJ(4) * t121;
t118 = qJD(5) + t55;
t6 = t79 * qJD(1) + t118;
t156 = -(qJD(6) * t19 + t22) * t51 + (t15 * qJD(3) - qJD(6) * t6 - t12) * t68;
t67 = cos(qJ(6));
t127 = qJD(6) * t67;
t109 = t66 * t127;
t65 = sin(qJ(6));
t123 = t65 * qJD(3);
t17 = (t68 * t123 + t109) * qJD(1) - qJD(6) * t123;
t8 = t66 * t114 + (qJD(4) + t26) * qJD(3);
t154 = t8 * t65;
t153 = t8 * t67;
t106 = t68 * t115;
t122 = t67 * qJD(3);
t130 = qJD(1) * t66;
t33 = t65 * t130 + t122;
t16 = -t33 * qJD(6) + t67 * t106;
t152 = t16 * t65;
t149 = t33 * t51;
t148 = t33 * t66;
t35 = t67 * t130 - t123;
t147 = t35 * t51;
t146 = t35 * t66;
t145 = t51 * t65;
t144 = t51 * t67;
t143 = t51 * t68;
t142 = t66 * t16;
t141 = t67 * t68;
t71 = qJD(3) ^ 2;
t139 = t71 * t66;
t138 = t71 * t68;
t64 = qJ(4) + pkin(5);
t62 = t66 ^ 2;
t63 = t68 ^ 2;
t136 = t62 - t63;
t72 = qJD(1) ^ 2;
t135 = t71 + t72;
t134 = qJ(4) * t66;
t132 = t72 * qJ(2);
t128 = qJD(6) * t65;
t126 = qJD(6) * t68;
t125 = t15 * qJD(6);
t124 = t21 * qJD(3);
t57 = t68 * qJD(4);
t92 = t69 * t66 - qJ(2);
t20 = t92 * qJD(1) + t118;
t119 = qJD(5) + t20;
t117 = qJ(2) * qJD(3);
t113 = t65 * t143;
t112 = t51 * t141;
t110 = t66 * t128;
t108 = qJD(2) * t111;
t105 = t66 * pkin(3) + qJ(2);
t83 = t69 * t68 - t134;
t75 = t83 * qJD(3) - qJD(2);
t18 = t57 + t75;
t52 = qJD(1) * t57;
t7 = t75 * qJD(1) + t52;
t103 = qJD(1) * t18 + t7;
t102 = t119 * t68;
t32 = t58 + t92;
t101 = qJD(1) * t32 + t20;
t96 = -0.2e1 * t107;
t95 = qJD(1) + t126;
t44 = t65 * t107;
t94 = -t51 * t127 + t44;
t45 = t67 * t107;
t93 = -t51 * t128 - t45;
t11 = qJD(3) * t60 + t159;
t2 = t67 * t11 + t65 * t6;
t91 = t65 * t11 - t67 * t6;
t90 = pkin(3) * t68 + t134;
t40 = t131 * t68;
t78 = t60 * t68 - t64 * t66;
t74 = t78 * qJD(3) - qJD(2);
t89 = (qJD(6) * t40 + t57 + t74) * t51;
t88 = (-t51 - t121) * t66;
t87 = qJD(1) * t62 - t143;
t30 = t105 * qJD(1) - t55;
t41 = t105 - t58;
t86 = qJD(3) * (qJD(1) * t41 + t30);
t77 = t90 * qJD(3) + qJD(2);
t14 = t77 * qJD(1) - t52;
t24 = -t57 + t77;
t84 = qJD(1) * t24 - t70 * t71 + t14;
t76 = t8 * t66 + t13 * t129 + (-t12 - t124) * t68;
t59 = 0.2e1 * qJD(3) * qJD(4);
t50 = t68 * t72 * t66;
t48 = -t63 * t72 - t71;
t43 = t135 * t68;
t42 = t135 * t66;
t39 = t131 * t66;
t37 = t90 * qJD(1);
t27 = t83 * qJD(1);
t25 = t38 + t54;
t23 = t66 * qJD(5) + t68 * t99;
t10 = t78 * qJD(1);
t4 = t74 * qJD(1) + t52;
t3 = t67 * t4;
t1 = [0, 0, 0, 0, t108, qJ(2) * t108, t68 * t96, 0.2e1 * t136 * t115, -t139, -t138, 0, -t70 * t139 + (qJD(2) * t66 + t68 * t117) * t111, -t70 * t138 + (qJD(2) * t68 - t66 * t117) * t111, t84 * t66 + t68 * t86, t158, t66 * t86 - t68 * t84, t14 * t41 - t158 * t70 + t30 * t24, t103 * t68 + (-t101 * t66 + t23) * qJD(3), t103 * t66 + (t101 * t68 + t22) * qJD(3) (-t22 * t68 + t23 * t66 + (t39 * t68 - t40 * t66) * qJD(3)) * qJD(1) + t76, -t12 * t40 + t13 * t22 + t20 * t18 - t21 * t23 + t7 * t32 + t8 * t39, t67 * t142 + (t68 * t122 - t110) * t35 (-t33 * t67 - t35 * t65) * t68 * qJD(3) + (-t152 - t17 * t67 + (t33 * t65 - t35 * t67) * qJD(6)) * t66, -t51 * t110 + t16 * t68 + (-t67 * t87 - t146) * qJD(3), -t51 * t109 - t17 * t68 + (t65 * t87 + t148) * qJD(3), qJD(3) * t88, t39 * t17 + t23 * t33 + t3 * t68 + (-t11 * t126 + t89) * t67 + t156 * t65 + (t67 * t125 + t154 + (-(t67 * t19 + t65 * t40) * qJD(1) + t91) * qJD(3)) * t66, t39 * t16 + t23 * t35 + (-t89 - (-qJD(6) * t11 + t4) * t68) * t65 + t156 * t67 + (-t65 * t125 + t153 + ((t65 * t19 - t67 * t40) * qJD(1) + t2) * qJD(3)) * t66; 0, 0, 0, 0, -t72, -t132, 0, 0, 0, 0, 0, -t42, -t43, -t42, 0, t43, -t30 * qJD(1) - t158, t43, t42, 0, t20 * qJD(1) + t76, 0, 0, 0, 0, 0, t66 * t17 + t95 * t144 + (t33 * t68 + t65 * t88) * qJD(3), t142 - t95 * t145 + (t35 * t68 + t67 * t88) * qJD(3); 0, 0, 0, 0, 0, 0, t50, -t136 * t72, 0, 0, 0, -t68 * t132, t66 * t132 (-t30 * t68 - t37 * t66) * qJD(1) ((t31 - t61) * t68 + (t29 + t98) * t66) * qJD(1), t59 + (-t30 * t66 + t37 * t68) * qJD(1), t28 * qJ(4) + t31 * qJD(4) - t30 * t37 + (-t150 + (-t29 - t133) * t66) * t47, t59 + (-t26 + t140) * qJD(3) + ((qJ(5) * qJD(3) - t27) * t68 + t119 * t66) * qJD(1), -t25 * qJD(3) + (-t27 * t66 - t102) * qJD(1) + t137 (t21 + t25 + t61) * t121, t8 * qJ(4) + t12 * t69 - t13 * t25 - t159 * t21 - t20 * t27, -t144 * t35 - t152 (-t16 + t149) * t67 + (t17 + t147) * t65 (-t112 + t146) * qJD(1) + t94 (t113 - t148) * qJD(1) - t93, t51 * t130, t64 * t17 + t153 - (t67 * t10 - t65 * t25) * t51 + t159 * t33 + (-t60 * t144 - t15 * t65) * qJD(6) + (-t91 * t66 + (t60 * t129 - t15 * t68) * t65) * qJD(1), t64 * t16 - t154 + (t65 * t10 + t67 * t25) * t51 + t159 * t35 + (t60 * t145 - t15 * t67) * qJD(6) + (-t15 * t141 + (t60 * t122 - t2) * t66) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t48, -t31 * qJD(3) + t30 * t121 + t36, t48, -t50, 0, -qJD(1) * t102 + t124 + t137, 0, 0, 0, 0, 0, -qJD(3) * t33 - t144 * t51 + t44, t51 ^ 2 * t65 - qJD(3) * t35 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0.2e1 * t106 (-t62 - t63) * t72, t52 + (t13 * t68 + t21 * t66 + t75) * qJD(1), 0, 0, 0, 0, 0 (-t113 - t148) * qJD(1) + t93 (-t112 - t146) * qJD(1) + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t33, -t33 ^ 2 + t35 ^ 2, t16 + t149, t147 - t17, -t107, -t65 * t12 - t15 * t35 - t160 * t2 + t3, -t67 * t12 + t15 * t33 + t160 * t91 - t65 * t4;];
tauc_reg  = t1;
