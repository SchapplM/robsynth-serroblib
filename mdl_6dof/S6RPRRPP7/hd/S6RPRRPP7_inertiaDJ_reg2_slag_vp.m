% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPP7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:12
% EndTime: 2019-03-09 04:52:20
% DurationCPUTime: 2.83s
% Computational Cost: add. (1411->262), mult. (3097->387), div. (0->0), fcn. (2177->4), ass. (0->146)
t81 = sin(qJ(4));
t152 = t81 * qJ(5);
t167 = pkin(4) + pkin(5);
t83 = cos(qJ(4));
t100 = t167 * t83 + t152;
t142 = t83 * qJD(5);
t179 = qJD(4) * t100 - t142;
t77 = t81 ^ 2;
t79 = t83 ^ 2;
t157 = t77 + t79;
t84 = cos(qJ(3));
t72 = t84 * qJD(3);
t34 = t157 * t72;
t139 = qJ(5) * qJD(3);
t69 = t84 * t139;
t82 = sin(qJ(3));
t175 = qJD(5) * t82 + t69;
t178 = t167 * qJD(3);
t177 = 0.4e1 * t84;
t176 = 0.2e1 * t175;
t158 = t77 - t79;
t86 = -pkin(1) - pkin(7);
t154 = qJ(5) * t83;
t97 = -t167 * t81 + t154;
t174 = t86 + t97;
t173 = t158 * qJD(4);
t107 = pkin(4) * t81 - t154;
t145 = t81 * qJD(5);
t30 = qJD(4) * t107 - t145;
t163 = t84 * t30;
t165 = pkin(8) * t84;
t98 = t107 - t86;
t28 = t98 * t84;
t108 = pkin(4) * t83 + t152;
t47 = -pkin(3) - t108;
t171 = qJD(3) * (t47 * t82 + t165) - qJD(4) * t28 - t163;
t141 = t83 * qJD(6);
t149 = qJD(4) * t81;
t161 = pkin(8) - qJ(6);
t29 = -t149 * t161 - t141;
t144 = t81 * qJD(6);
t53 = t161 * t83;
t31 = qJD(4) * t53 - t144;
t52 = t161 * t81;
t170 = -(t52 * t83 - t53 * t81) * qJD(4) - t29 * t83 - t31 * t81;
t169 = qJD(4) * t108 - t142;
t143 = t82 * qJD(3);
t125 = t81 * t143;
t116 = t83 * t125;
t80 = t84 ^ 2;
t12 = t116 * t177 + 0.2e1 * t173 * t80;
t78 = t82 ^ 2;
t156 = t78 - t80;
t118 = t156 * qJD(3);
t147 = qJD(4) * t84;
t129 = t81 * t147;
t14 = 0.2e1 * t118 * t83 + 0.2e1 * t129 * t82;
t127 = t83 * t147;
t151 = qJD(3) * t81;
t15 = 0.2e1 * t127 * t82 - 0.2e1 * t151 * t156;
t45 = -0.2e1 * t173;
t168 = 2 * qJD(2);
t87 = 0.2e1 * qJD(5);
t166 = pkin(8) * t82;
t164 = t82 * t86;
t162 = t84 * t86;
t160 = t142 * t82 + t69 * t83;
t112 = pkin(3) * t82 - t165;
t46 = qJ(2) + t112;
t57 = t83 * t164;
t26 = t46 * t81 + t57;
t159 = pkin(8) * t34;
t155 = t78 + t80;
t153 = qJ(6) * t84;
t150 = qJD(3) * t83;
t148 = qJD(4) * t82;
t73 = qJD(4) * t83;
t146 = qJD(4) * t86;
t140 = qJ(2) * qJD(3);
t138 = -0.2e1 * pkin(3) * qJD(4);
t56 = t81 * t164;
t121 = t86 * t72;
t113 = pkin(3) * t84 + t166;
t94 = qJD(3) * t113 + qJD(2);
t137 = t121 * t83 + t46 * t73 + t81 * t94;
t17 = qJ(5) * t82 + t26;
t136 = pkin(4) * t72;
t134 = pkin(8) * t149;
t133 = pkin(8) * t73;
t132 = t81 * t72;
t131 = t83 * t143;
t130 = t81 * t148;
t128 = t81 * t146;
t126 = t83 * t146;
t124 = t81 * t73;
t123 = t86 * t143;
t122 = t82 * t72;
t120 = t82 * t139;
t119 = qJ(5) * t147;
t25 = t46 * t83 - t56;
t62 = 0.2e1 * t122;
t117 = t121 * t81 + t126 * t82 + t149 * t46 - t83 * t94;
t115 = t80 * t124;
t3 = -t143 * t174 - t179 * t84;
t43 = pkin(3) + t100;
t114 = -t147 * t43 + t3;
t10 = t153 * t81 + t17;
t9 = t56 + (-t46 - t153) * t83 - t167 * t82;
t110 = t10 * t83 + t81 * t9;
t109 = -t10 * t81 + t83 * t9;
t18 = -t82 * pkin(4) - t25;
t106 = t17 * t83 + t18 * t81;
t105 = t17 * t81 - t18 * t83;
t104 = t25 * t83 + t26 * t81;
t103 = t25 * t81 - t26 * t83;
t96 = -qJ(6) * t149 + t141;
t39 = t73 * t82 + t132;
t6 = t128 * t82 - t137;
t8 = -t143 * t98 + t169 * t84;
t95 = -t8 + (t47 * t84 - t166) * qJD(4);
t93 = -qJ(6) * t131 - t117;
t16 = t174 * t84;
t27 = qJD(4) * t97 + t145;
t92 = -qJD(4) * t16 + t143 * t43 - t27 * t84;
t4 = -t6 + t175;
t5 = t117 - t136;
t90 = -qJD(4) * t105 + t4 * t83 + t5 * t81;
t89 = -qJD(4) * t104 + t117 * t81 - t6 * t83;
t88 = qJ(6) * t127 + t84 * t144 + (-qJ(6) * qJD(3) - t146) * t82 * t81 + t137;
t76 = qJ(2) * t168;
t75 = qJ(5) * t87;
t61 = -0.2e1 * t124;
t60 = 0.2e1 * t124;
t40 = -t125 + t127;
t38 = t155 * t73;
t37 = t129 + t131;
t36 = -t72 * t83 + t130;
t35 = t155 * t149;
t24 = -0.2e1 * t122 * t79 - 0.2e1 * t115;
t23 = -0.2e1 * t122 * t77 + 0.2e1 * t115;
t22 = t147 * t158 + t116;
t19 = t124 * t177 - t143 * t158;
t11 = (-0.1e1 + t157) * t62;
t2 = t88 + t175;
t1 = (-t96 - t178) * t84 - t93;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t76, -0.2e1 * t122, 0.2e1 * t118, 0, t62, 0, 0, 0.2e1 * qJD(2) * t82 + 0.2e1 * t140 * t84, 0.2e1 * qJD(2) * t84 - 0.2e1 * t140 * t82, 0, t76, t24, t12, -t14, t23, -t15, t62, -0.2e1 * t80 * t126 - 0.2e1 * t117 * t82 + 0.2e1 * (t25 + 0.2e1 * t56) * t72, 0.2e1 * t80 * t128 + 0.2e1 * t6 * t82 + 0.2e1 * (-t26 + 0.2e1 * t57) * t72, 0.2e1 * t104 * t143 + 0.2e1 * (qJD(4) * t103 + t117 * t83 + t6 * t81) * t84, -0.2e1 * t122 * t86 ^ 2 - 0.2e1 * t117 * t25 - 0.2e1 * t26 * t6, t24, -t14, -t12, t62, t15, t23, 0.2e1 * (-t151 * t28 - t5) * t82 + 0.2e1 * (-qJD(3) * t18 + t28 * t73 + t8 * t81) * t84, 0.2e1 * t105 * t143 + 0.2e1 * (-qJD(4) * t106 - t4 * t81 + t5 * t83) * t84, 0.2e1 * (t150 * t28 + t4) * t82 + 0.2e1 * (qJD(3) * t17 + t149 * t28 - t8 * t83) * t84, 0.2e1 * t17 * t4 + 0.2e1 * t18 * t5 + 0.2e1 * t28 * t8, t24, -t12, t14, t23, -t15, t62, 0.2e1 * (t151 * t16 - t1) * t82 + 0.2e1 * (-qJD(3) * t9 - t16 * t73 - t3 * t81) * t84, 0.2e1 * (-t150 * t16 + t2) * t82 + 0.2e1 * (qJD(3) * t10 - t149 * t16 + t3 * t83) * t84, 0.2e1 * t109 * t143 + 0.2e1 * (qJD(4) * t110 - t1 * t83 + t2 * t81) * t84, 0.2e1 * t1 * t9 + 0.2e1 * t10 * t2 + 0.2e1 * t16 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t35, 0, -t103 * t72 + (t89 - 0.2e1 * t121) * t82, 0, 0, 0, 0, 0, 0, -t38, 0, -t35 (qJD(3) * t106 - t8) * t84 + (qJD(3) * t28 + t90) * t82, 0, 0, 0, 0, 0, 0, -t38, -t35, 0 (qJD(3) * t110 + t3) * t84 + (-qJD(3) * t16 + qJD(4) * t109 + t1 * t81 + t2 * t83) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, 0, -t72, 0, -t123, -t121, 0, 0, -t22, -t19, t39, t22, -t36, 0 (-t113 * t83 - t162 * t81) * qJD(4) + (t112 * t81 - t57) * qJD(3) (t113 * t81 - t162 * t83) * qJD(4) + (-t83 * t165 + (pkin(3) * t83 + t81 * t86) * t82) * qJD(3), t89, -pkin(3) * t123 + pkin(8) * t89, -t22, t39, t19, 0, t36, t22, -t171 * t81 + t83 * t95, t90, t171 * t83 + t81 * t95, pkin(8) * t90 + t28 * t30 + t8 * t47, -t22, t19, -t39, t22, -t36, 0, t114 * t83 - t31 * t82 - t52 * t72 + t81 * t92, t114 * t81 + t29 * t82 + t53 * t72 - t83 * t92 (t52 * t143 - t31 * t84 - t2 + (t53 * t84 - t9) * qJD(4)) * t83 + (-t53 * t143 + t29 * t84 - t1 + (t52 * t84 + t10) * qJD(4)) * t81, t1 * t52 + t10 * t29 + t16 * t27 + t2 * t53 + t3 * t43 + t31 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, -t72, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t40, t34, -pkin(3) * t143 + t159, 0, 0, 0, 0, 0, 0, -t37, t34, t40, t143 * t47 + t159 - t163, 0, 0, 0, 0, 0, 0, -t37, t40, -t34 (t27 + (t52 * t81 + t53 * t83) * qJD(3)) * t84 + (-qJD(3) * t43 - t170) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t45, 0, t61, 0, 0, t81 * t138, t83 * t138, 0, 0, t60, 0, -t45, 0, 0, t61, 0.2e1 * t149 * t47 - 0.2e1 * t30 * t83, 0, -0.2e1 * t30 * t81 - 0.2e1 * t47 * t73, 0.2e1 * t47 * t30, t60, -t45, 0, t61, 0, 0, -0.2e1 * t149 * t43 + 0.2e1 * t27 * t83, 0.2e1 * t27 * t81 + 0.2e1 * t43 * t73, 0.2e1 * t170, 0.2e1 * t27 * t43 + 0.2e1 * t29 * t53 + 0.2e1 * t31 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, -t40, t72, -t117, t6, 0, 0, 0, -t37, 0, t72, t40, 0, -t117 + 0.2e1 * t136 (pkin(4) * t143 - t119) * t83 + (t120 + (pkin(4) * qJD(4) - qJD(5)) * t84) * t81, -t6 + t176, -pkin(4) * t5 + qJ(5) * t4 + qJD(5) * t17, 0, 0, t37, 0, -t40, t72 (t96 + 0.2e1 * t178) * t84 + t93, t88 + t176 (-t143 * t167 + t119) * t83 + (-t120 + (-qJD(4) * t167 + qJD(5)) * t84) * t81, qJ(5) * t2 + qJD(5) * t10 - t1 * t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t36, 0, 0, 0, 0, 0, 0, 0, 0, -t39, 0, -t36, -pkin(4) * t39 - qJ(5) * t130 + t160, 0, 0, 0, 0, 0, 0, -t39, -t36, 0, -t100 * t148 - t132 * t167 + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, -t149, 0, -t133, t134, 0, 0, 0, t73, 0, 0, t149, 0, -t133, -t169, -t134, -t169 * pkin(8), 0, 0, -t73, 0, -t149, 0, -t31, t29, t179, qJ(5) * t29 + qJD(5) * t53 - t167 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t75, 0, 0, 0, 0, 0, 0, 0, t87, 0, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t37, 0, t5, 0, 0, 0, 0, 0, 0, -t72, 0, t37, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, t133, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t37, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t73, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
