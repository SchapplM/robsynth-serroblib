% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:13:22
% EndTime: 2019-03-09 18:13:28
% DurationCPUTime: 1.92s
% Computational Cost: add. (2663->217), mult. (5813->353), div. (0->0), fcn. (5509->8), ass. (0->133)
t162 = sin(qJ(3));
t164 = cos(qJ(3));
t90 = sin(qJ(2));
t92 = cos(qJ(2));
t113 = -t162 * t90 + t164 * t92;
t57 = t162 * t92 + t164 * t90;
t170 = -pkin(3) * t113 - t57 * qJ(4);
t91 = cos(qJ(6));
t87 = t91 ^ 2;
t88 = sin(qJ(6));
t156 = t88 ^ 2 - t87;
t141 = qJD(6) * t156;
t169 = qJD(2) + qJD(3);
t166 = -pkin(8) - pkin(7);
t137 = t162 * t166;
t123 = qJD(2) * t137;
t139 = t166 * t164;
t125 = qJD(2) * t139;
t144 = t164 * qJD(3);
t66 = t166 * t92;
t21 = (qJD(3) * t137 + t123) * t90 - t92 * t125 - t66 * t144;
t168 = -0.2e1 * t141;
t93 = 2 * qJD(4);
t167 = -pkin(3) - pkin(4);
t165 = t90 * pkin(2);
t77 = -t92 * pkin(2) - pkin(1);
t163 = cos(qJ(5));
t114 = t163 * t113;
t89 = sin(qJ(5));
t83 = qJD(5) * t89;
t108 = t113 * qJD(3);
t96 = t113 * qJD(2) + t108;
t109 = t57 * qJD(3);
t97 = t57 * qJD(2) + t109;
t13 = qJD(5) * t114 - t163 * t96 + t57 * t83 - t89 * t97;
t117 = t89 * t113;
t34 = t163 * t57 - t117;
t161 = t34 * t13;
t160 = t34 * t88;
t159 = t34 * t91;
t158 = t91 * t13;
t38 = t90 * t137 - t164 * t66;
t138 = t163 * t167;
t157 = t163 * qJD(4) + qJD(5) * t138;
t126 = t90 * t139;
t37 = -t162 * t66 - t126;
t110 = -t57 * pkin(9) + t37;
t105 = t163 * t110;
t27 = -pkin(9) * t113 + t38;
t17 = t89 * t27 - t105;
t155 = qJD(6) * t17;
t82 = qJD(6) * t88;
t154 = qJD(6) * t91;
t153 = t90 * qJD(2);
t152 = t92 * qJD(2);
t151 = -0.2e1 * pkin(1) * qJD(2);
t76 = -t164 * pkin(2) - pkin(3);
t127 = -pkin(4) + t76;
t111 = t163 * t127;
t143 = t162 * qJD(3);
t140 = pkin(2) * t143;
t78 = pkin(2) * t144;
t70 = t78 + qJD(4);
t150 = qJD(5) * t111 + t89 * t140 + t163 * t70;
t149 = t162 * pkin(2);
t148 = t88 * t154;
t146 = 0.4e1 * t88 * t159;
t84 = qJD(5) * t163;
t145 = qJD(6) * t163;
t28 = t170 + t77;
t74 = t149 + qJ(4);
t40 = t89 * t74 + pkin(5) - t111;
t61 = t89 * qJ(4) + pkin(5) - t138;
t142 = qJD(6) * (-t40 - t61);
t14 = -qJD(5) * t117 - t163 * t97 + t57 * t84 + t89 * t96;
t136 = pkin(5) * t13 - pkin(10) * t14;
t33 = t57 * t89 + t114;
t135 = pkin(5) * t34 + pkin(10) * t33;
t106 = t89 * t127 + t163 * t74;
t115 = t163 * t140;
t26 = t106 * qJD(5) + t89 * t70 - t115;
t23 = t26 * t91;
t134 = -t40 * t82 + t23;
t112 = t163 * qJ(4) + t89 * t167;
t43 = t89 * qJD(4) + t112 * qJD(5);
t39 = t43 * t91;
t133 = -t61 * t82 + t39;
t22 = pkin(4) * t113 - t28;
t11 = t33 * pkin(5) - t34 * pkin(10) + t22;
t107 = t89 * t110;
t18 = t163 * t27 + t107;
t131 = t11 * t91 - t18 * t88;
t130 = t11 * t88 + t18 * t91;
t41 = -pkin(10) + t106;
t129 = t33 * t41 - t34 * t40;
t62 = -pkin(10) + t112;
t128 = t33 * t62 - t34 * t61;
t121 = -t88 * t13 + t34 * t154;
t120 = -t34 * t82 - t158;
t119 = -t40 * t154 - t26 * t88;
t118 = -t61 * t154 - t43 * t88;
t116 = t163 * t34 + t33 * t89;
t98 = t169 * t57;
t19 = pkin(2) * t153 + t98 * pkin(3) - t96 * qJ(4) - t57 * qJD(4);
t20 = -qJD(3) * t126 - t92 * t123 - t90 * t125 - t66 * t143;
t25 = t74 * t83 - t150;
t104 = -t13 * t40 - t14 * t41 + t25 * t33 + t26 * t34 - t155;
t42 = qJ(4) * t83 - t157;
t103 = -t13 * t61 - t14 * t62 + t33 * t42 + t34 * t43 - t155;
t102 = t163 * t13 - t14 * t89 + (-t163 * t33 + t34 * t89) * qJD(5);
t100 = t76 * t113 - t74 * t57;
t12 = -pkin(4) * t97 - t19;
t95 = t97 * pkin(9) - t20;
t94 = -t96 * pkin(9) + t21;
t80 = pkin(5) * t154;
t79 = pkin(5) * t82;
t73 = -0.2e1 * t140;
t69 = -0.2e1 * t148;
t68 = 0.2e1 * t148;
t47 = t88 * t145 + t91 * t83;
t46 = -t91 * t145 + t88 * t83;
t32 = t34 ^ 2;
t10 = t14 * t88 + t33 * t154;
t9 = -t14 * t91 + t33 * t82;
t8 = t34 * t141 + t88 * t158;
t7 = qJD(6) * t146 - t156 * t13;
t6 = qJD(5) * t107 - t163 * t94 + t27 * t84 + t89 * t95;
t5 = -qJD(5) * t105 - t163 * t95 + t27 * t83 - t89 * t94;
t4 = t6 * t91;
t3 = t14 * pkin(5) + t13 * pkin(10) + t12;
t2 = -t130 * qJD(6) + t3 * t91 + t5 * t88;
t1 = -t131 * qJD(6) - t3 * t88 + t5 * t91;
t15 = [0, 0, 0, 0.2e1 * t90 * t152, 0.2e1 * (-t90 ^ 2 + t92 ^ 2) * qJD(2), 0, 0, 0, t90 * t151, t92 * t151, 0.2e1 * t57 * t96, 0.2e1 * t113 * t96 - 0.2e1 * t57 * t98, 0, 0, 0, 0.2e1 * t77 * t109 + 0.2e1 * (-t113 * t165 + t77 * t57) * qJD(2), 0.2e1 * t77 * t108 + 0.2e1 * (t77 * t113 + t57 * t165) * qJD(2), -0.2e1 * t113 * t19 + 0.2e1 * t28 * t97, -0.2e1 * t20 * t113 + 0.2e1 * t21 * t57 + 0.2e1 * t169 * (t37 * t113 - t38 * t57) -0.2e1 * t28 * t169 * t113 - 0.2e1 * t19 * t57, 0.2e1 * t19 * t28 - 0.2e1 * t20 * t38 + 0.2e1 * t21 * t37, -0.2e1 * t161, 0.2e1 * t13 * t33 - 0.2e1 * t14 * t34, 0, 0, 0, 0.2e1 * t12 * t33 + 0.2e1 * t14 * t22, 0.2e1 * t12 * t34 - 0.2e1 * t13 * t22, -0.2e1 * t32 * t148 - 0.2e1 * t87 * t161, t13 * t146 + 0.2e1 * t32 * t141, 0.2e1 * t120 * t33 + 0.2e1 * t14 * t159, -0.2e1 * t121 * t33 - 0.2e1 * t14 * t160, 0.2e1 * t33 * t14, 0.2e1 * t121 * t17 + 0.2e1 * t131 * t14 + 0.2e1 * t6 * t160 + 0.2e1 * t2 * t33, 0.2e1 * t1 * t33 + 0.2e1 * t120 * t17 - 0.2e1 * t130 * t14 + 0.2e1 * t6 * t159; 0, 0, 0, 0, 0, t152, -t153, 0, -pkin(7) * t152, pkin(7) * t153, 0, 0, t96, -t98, 0, -t21, t20, -t21, t70 * t113 + t100 * qJD(2) + (t57 * t149 + t100) * qJD(3), -t20, t37 * t140 - t20 * t74 + t21 * t76 + t38 * t70, 0, 0, t13, t14, 0, t6, -t5, t8, t7, -t10, t9, 0, t104 * t88 - t129 * t154 + t4 (qJD(6) * t129 - t6) * t88 + t104 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -0.2e1 * t78, t73, 0, 0.2e1 * t70, 0.2e1 * t76 * t140 + 0.2e1 * t74 * t70, 0, 0, 0, 0, 0, 0.2e1 * t26, -0.2e1 * t25, t68, t168, 0, 0, 0, 0.2e1 * t134, 0.2e1 * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, -t98, 0, -t21, t20, -t21, t113 * qJD(4) + t169 * t170, -t20, -pkin(3) * t21 - qJ(4) * t20 + qJD(4) * t38, 0, 0, t13, t14, 0, t6, -t5, t8, t7, -t10, t9, 0, t103 * t88 - t128 * t154 + t4 (qJD(6) * t128 - t6) * t88 + t103 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, -t78, -t140, 0, t93 + t78, -pkin(3) * t140 + t70 * qJ(4) + t74 * qJD(4), 0, 0, 0, 0, 0, -t115 + (qJD(4) + t70) * t89 + (t106 + t112) * qJD(5) (-qJ(4) - t74) * t83 + t150 + t157, t68, t168, 0, 0, 0, t88 * t142 + t23 + t39 (-t26 - t43) * t88 + t91 * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, qJ(4) * t93, 0, 0, 0, 0, 0, 0.2e1 * t43, -0.2e1 * t42, t68, t168, 0, 0, 0, 0.2e1 * t133, 0.2e1 * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * t88 - t116 * t154, t102 * t91 + t116 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, 0, 0, 0, 0, 0, t83, t84, 0, 0, 0, 0, 0, t47, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t84, 0, 0, 0, 0, 0, t47, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, -t6, t5, -t8, -t7, t10, -t9, 0, -t4 + t136 * t88 + (-t135 * t91 + t17 * t88) * qJD(6), t6 * t88 + t136 * t91 + (t135 * t88 + t17 * t91) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t25, t69, -t168, 0, 0, 0, -t134 + t79, -t119 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t42, t69, -t168, 0, 0, 0, -t133 + t79, -t118 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t84, 0, 0, 0, 0, 0, -t47, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t168, 0, 0, 0, -0.2e1 * t79, -0.2e1 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t121, t14, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t82, 0, -t41 * t154 + t25 * t88, t25 * t91 + t41 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t82, 0, -t62 * t154 + t42 * t88, t42 * t91 + t62 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89 * t154 - t88 * t84, t89 * t82 - t84 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, -t82, 0, -pkin(10) * t154, pkin(10) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t15;
