% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:45:56
% EndTime: 2019-03-09 16:46:02
% DurationCPUTime: 1.75s
% Computational Cost: add. (2799->218), mult. (6005->350), div. (0->0), fcn. (5412->6), ass. (0->138)
t96 = sin(qJ(3));
t99 = cos(qJ(2));
t160 = t96 * t99;
t167 = cos(qJ(3));
t97 = sin(qJ(2));
t63 = t167 * t97 + t160;
t85 = -pkin(2) * t99 - pkin(1);
t113 = -qJ(4) * t63 + t85;
t172 = pkin(3) + pkin(9);
t132 = t167 * t99;
t161 = t96 * t97;
t62 = -t132 + t161;
t31 = t172 * t62 + t113;
t171 = -pkin(8) - pkin(7);
t72 = t171 * t97;
t73 = t171 * t99;
t42 = -t167 * t72 - t96 * t73;
t36 = t63 * pkin(4) + t42;
t95 = sin(qJ(5));
t98 = cos(qJ(5));
t177 = t98 * t31 + t95 * t36;
t43 = t167 * t73 - t96 * t72;
t121 = t95 * pkin(5) - qJ(6) * t98;
t93 = t95 ^ 2;
t94 = t98 ^ 2;
t152 = t93 - t94;
t126 = t152 * qJD(5);
t176 = qJD(2) + qJD(3);
t137 = qJ(6) * qJD(5);
t140 = t98 * qJD(6);
t144 = qJD(5) * t98;
t175 = pkin(5) * t144 + t95 * t137 - t140;
t123 = qJD(2) * t132;
t127 = t167 * qJD(3);
t39 = -t99 * t127 + t176 * t161 - t123;
t141 = t97 * qJD(2);
t88 = pkin(2) * t141;
t110 = qJ(4) * t39 - qJD(4) * t63 + t88;
t40 = t176 * t63;
t10 = t172 * t40 + t110;
t129 = qJD(2) * t171;
t67 = t97 * t129;
t26 = -t43 * qJD(3) - t171 * t123 + t96 * t67;
t18 = -t39 * pkin(4) + t26;
t6 = -qJD(5) * t177 - t10 * t95 + t98 * t18;
t142 = qJD(6) * t95;
t174 = t121 * qJD(5) - t142;
t101 = 0.2e1 * qJD(4);
t173 = 0.2e1 * qJD(6);
t170 = pkin(2) * t96;
t169 = pkin(5) * t39;
t122 = pkin(5) * t98 + qJ(6) * t95;
t112 = -pkin(4) - t122;
t23 = t112 * t62 - t43;
t147 = qJD(3) * t96;
t25 = -t72 * t127 - t129 * t160 - t73 * t147 - t167 * t67;
t8 = t112 * t40 + t174 * t62 - t25;
t168 = t23 * t144 + t8 * t95;
t84 = -t167 * pkin(2) - pkin(3);
t79 = -pkin(9) + t84;
t166 = t39 * t79;
t165 = t62 * t95;
t164 = t62 * t98;
t163 = t63 * t79;
t162 = t95 * t40;
t159 = t98 * t40;
t17 = -pkin(4) * t40 - t25;
t37 = -pkin(4) * t62 - t43;
t158 = t37 * t144 + t17 * t95;
t124 = pkin(2) * t127;
t76 = t124 + qJD(4);
t44 = t76 + t175;
t71 = qJ(4) + t121;
t58 = t71 + t170;
t156 = t58 * t144 + t44 * t95;
t50 = qJD(4) + t175;
t155 = t71 * t144 + t50 * t95;
t80 = qJ(4) + t170;
t154 = t80 * t144 + t76 * t95;
t138 = qJ(4) * qJD(5);
t153 = qJD(4) * t95 + t98 * t138;
t151 = qJ(6) * t39;
t149 = t172 * t39;
t148 = t172 * t63;
t146 = qJD(5) * t37;
t145 = qJD(5) * t95;
t143 = qJD(6) * t63;
t139 = t99 * qJD(2);
t136 = qJD(5) * t172;
t135 = -0.2e1 * pkin(1) * qJD(2);
t134 = t95 * t159;
t87 = pkin(2) * t147;
t133 = t95 * t144;
t131 = t95 * t136;
t130 = t98 * t136;
t125 = t63 * t87;
t11 = qJ(6) * t63 + t177;
t117 = -t31 * t95 + t36 * t98;
t12 = -pkin(5) * t63 - t117;
t120 = t11 * t98 + t12 * t95;
t119 = t11 * t95 - t12 * t98;
t116 = t58 * t62 - t163;
t115 = -t62 * t71 - t148;
t51 = (t93 + t94) * t87;
t114 = -qJ(4) * t40 - qJD(4) * t62;
t30 = t63 * t144 - t39 * t95;
t28 = -t63 * t145 - t39 * t98;
t29 = t62 * t144 + t162;
t111 = t62 * t145 - t159;
t5 = -t98 * t10 - t36 * t144 + t31 * t145 - t95 * t18;
t109 = qJD(5) * (t62 * t80 - t163);
t108 = qJD(5) * (qJ(4) * t62 + t148);
t107 = -t40 * t71 - t50 * t62 + t149;
t106 = t114 + t149;
t105 = -t40 * t80 - t62 * t76 + t125;
t103 = -t40 * t58 - t44 * t62 + t125 - t166;
t102 = t105 - t166;
t2 = t143 - t5 - t151;
t4 = t169 - t6;
t1 = t120 * qJD(5) + t2 * t95 - t4 * t98;
t91 = qJD(4) * t98;
t75 = -0.2e1 * t133;
t69 = t76 * t98;
t61 = 0.2e1 * t126;
t60 = t62 ^ 2;
t53 = t71 * t145;
t52 = -pkin(5) * t145 + t98 * t137 + t142;
t48 = t58 * t145;
t47 = t79 * t144 + t95 * t87;
t46 = t79 * t145 - t98 * t87;
t38 = pkin(3) * t62 + t113;
t34 = -0.2e1 * t63 * t39;
t24 = -t62 * t126 + t134;
t21 = t23 * t145;
t20 = pkin(3) * t40 + t110;
t19 = -0.4e1 * t62 * t133 - t152 * t40;
t16 = t17 * t98;
t3 = [0, 0, 0, 0.2e1 * t97 * t139, 0.2e1 * (-t97 ^ 2 + t99 ^ 2) * qJD(2), 0, 0, 0, t97 * t135, t99 * t135, t34, 0.2e1 * t39 * t62 - 0.2e1 * t40 * t63, 0, 0, 0, 0.2e1 * t40 * t85 + 0.2e1 * t62 * t88, -0.2e1 * t39 * t85 + 0.2e1 * t63 * t88, 0.2e1 * t25 * t62 + 0.2e1 * t26 * t63 - 0.2e1 * t39 * t42 + 0.2e1 * t40 * t43, -0.2e1 * t20 * t62 - 0.2e1 * t38 * t40, -0.2e1 * t20 * t63 + 0.2e1 * t38 * t39, 0.2e1 * t20 * t38 + 0.2e1 * t25 * t43 + 0.2e1 * t26 * t42, 0.2e1 * t40 * t62 * t93 + 0.2e1 * t60 * t133, -0.2e1 * t60 * t126 + 0.4e1 * t62 * t134, 0.2e1 * t63 * t162 + 0.2e1 * t30 * t62, 0.2e1 * t63 * t159 + 0.2e1 * t28 * t62, t34, 0.2e1 * t111 * t37 - 0.2e1 * t117 * t39 - 0.2e1 * t17 * t164 + 0.2e1 * t6 * t63, 0.2e1 * t17 * t165 + 0.2e1 * t177 * t39 + 0.2e1 * t29 * t37 + 0.2e1 * t5 * t63, 0.2e1 * t111 * t23 + 0.2e1 * t12 * t39 - 0.2e1 * t8 * t164 - 0.2e1 * t4 * t63, 0.2e1 * t120 * t40 + 0.2e1 * (-t119 * qJD(5) + t2 * t98 + t4 * t95) * t62, -0.2e1 * t11 * t39 - 0.2e1 * t8 * t165 + 0.2e1 * t2 * t63 - 0.2e1 * t23 * t29, 0.2e1 * t11 * t2 + 0.2e1 * t12 * t4 + 0.2e1 * t23 * t8; 0, 0, 0, 0, 0, t139, -t141, 0, -pkin(7) * t139, pkin(7) * t141, 0, 0, -t39, -t40, 0, -t26, t25, -t39 * t84 + t105, t26, -t25, -t25 * t80 + t26 * t84 + t42 * t87 - t43 * t76, t24, t19, t28, -t30, 0, t102 * t98 + t109 * t95 + t158, t16 + t98 * t109 + (-t102 - t146) * t95, t103 * t98 + t116 * t145 + t168, -t1, t21 + (-qJD(5) * t116 - t8) * t98 + t103 * t95, t1 * t79 + t119 * t87 + t23 * t44 + t58 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t87, -0.2e1 * t124, 0, 0.2e1 * t87, 0.2e1 * t76, 0.2e1 * t76 * t80 + 0.2e1 * t84 * t87, t75, t61, 0, 0, 0, 0.2e1 * t154, -0.2e1 * t80 * t145 + 0.2e1 * t69, 0.2e1 * t156, -0.2e1 * t51, -0.2e1 * t44 * t98 + 0.2e1 * t48, 0.2e1 * t44 * t58 + 0.2e1 * t51 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t40, 0, -t26, t25, pkin(3) * t39 + t114, t26, -t25, -pkin(3) * t26 - qJ(4) * t25 - qJD(4) * t43, t24, t19, t28, -t30, 0, t106 * t98 + t108 * t95 + t158, t16 + t98 * t108 + (-t106 - t146) * t95, t107 * t98 - t115 * t145 + t168, -t1, t21 + (qJD(5) * t115 - t8) * t98 + t107 * t95, -t1 * t172 + t23 * t50 + t71 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t124, 0, t87, t101 + t124, -pkin(3) * t87 + qJ(4) * t76 + qJD(4) * t80, t75, t61, 0, 0, 0, t153 + t154, t69 + t91 + (-qJ(4) - t80) * t145, t155 + t156, -t51, t48 + t53 + (-t44 - t50) * t98, -t172 * t51 + t44 * t71 + t50 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, qJ(4) * t101, t75, t61, 0, 0, 0, 0.2e1 * t153, -0.2e1 * t95 * t138 + 0.2e1 * t91, 0.2e1 * t155, 0, -0.2e1 * t50 * t98 + 0.2e1 * t53, 0.2e1 * t71 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, 0, 0, t26, 0, 0, 0, 0, 0, t28, -t30, t28, 0, t30, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t111, -t39, t6, t5, t6 - 0.2e1 * t169, -t121 * t40 + (-qJD(5) * t122 + t140) * t62, 0.2e1 * t143 - t5 - 0.2e1 * t151, -pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t144, 0, -t46, -t47, -t46, -t52, t47, t122 * t87 - t174 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t144, 0, t131, t130, t131, -t52, -t130, t174 * t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t144, -t145, 0, t144, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, qJ(6) * t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t29, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, 0, -t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
