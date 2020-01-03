% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR15_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:24
% EndTime: 2019-12-31 18:37:30
% DurationCPUTime: 1.96s
% Computational Cost: add. (2587->278), mult. (5700->420), div. (0->0), fcn. (3649->6), ass. (0->157)
t118 = sin(qJ(5));
t120 = cos(qJ(5));
t116 = sin(pkin(8));
t121 = cos(qJ(3));
t168 = qJD(1) * t121;
t148 = t116 * t168;
t117 = cos(pkin(8));
t166 = qJD(3) * t117;
t83 = -t148 + t166;
t154 = t117 * t168;
t167 = qJD(3) * t116;
t84 = t154 + t167;
t34 = t118 * t84 - t120 * t83;
t202 = t34 ^ 2;
t37 = t118 * t83 + t120 * t84;
t201 = t37 ^ 2;
t119 = sin(qJ(3));
t169 = qJD(1) * t119;
t108 = qJD(5) + t169;
t200 = t108 * t34;
t159 = 0.2e1 * qJD(1);
t176 = t120 * t117;
t86 = t116 * t118 - t176;
t63 = t86 * t119;
t161 = qJD(5) * t120;
t162 = qJD(5) * t118;
t199 = -t116 * t162 + t117 * t161;
t122 = -pkin(1) - pkin(6);
t198 = qJD(1) * t122;
t114 = t119 ^ 2;
t115 = t121 ^ 2;
t197 = qJD(1) * (0.2e1 * t114 - t115);
t158 = pkin(7) * t117 * t119;
t126 = (pkin(4) * t121 + t158) * qJD(1);
t136 = pkin(3) * t121 + qJ(4) * t119;
t65 = t136 * qJD(3) - t121 * qJD(4) + qJD(2);
t48 = t65 * qJD(1);
t105 = qJD(2) + t198;
t182 = t105 * t121;
t67 = (qJD(4) + t182) * qJD(3);
t18 = -t116 * t67 + t117 * t48;
t12 = qJD(3) * t126 + t18;
t94 = pkin(3) * t119 - qJ(4) * t121 + qJ(2);
t77 = t94 * qJD(1);
t93 = t119 * t105;
t78 = qJD(3) * qJ(4) + t93;
t30 = -t116 * t78 + t117 * t77;
t15 = pkin(4) * t169 - pkin(7) * t84 + t30;
t31 = t116 * t77 + t117 * t78;
t17 = pkin(7) * t83 + t31;
t131 = t118 * t17 - t120 * t15;
t160 = qJD(1) * qJD(3);
t147 = t119 * t160;
t140 = t116 * t147;
t19 = t116 * t48 + t117 * t67;
t16 = pkin(7) * t140 + t19;
t1 = -t131 * qJD(5) + t118 * t12 + t120 * t16;
t181 = t116 * t121;
t89 = t136 * qJD(1);
t42 = -t105 * t181 + t117 * t89;
t23 = t126 + t42;
t155 = t116 * t169;
t179 = t117 * t121;
t43 = t105 * t179 + t116 * t89;
t29 = pkin(7) * t155 + t43;
t193 = pkin(7) + qJ(4);
t100 = t193 * t116;
t101 = t193 * t117;
t46 = -t100 * t120 - t101 * t118;
t196 = -t86 * qJD(4) + t46 * qJD(5) - t118 * t23 - t120 * t29;
t47 = -t100 * t118 + t101 * t120;
t87 = t116 * t120 + t117 * t118;
t195 = -t87 * qJD(4) - t47 * qJD(5) + t118 * t29 - t120 * t23;
t194 = t37 * t34;
t61 = t87 * t119;
t64 = t86 * t121;
t72 = t87 * qJD(1);
t192 = -qJD(3) * t64 - qJD(5) * t61 - t72;
t164 = qJD(3) * t121;
t191 = t86 * qJD(1) + qJD(5) * t63 - t87 * t164;
t74 = t87 * qJD(5);
t190 = t119 * t72 + t74;
t189 = -t118 * t155 + t169 * t176 + t199;
t138 = t120 * t147;
t139 = t118 * t147;
t188 = -t116 * t138 - t117 * t139;
t187 = qJD(3) * pkin(3);
t186 = t116 * t83;
t185 = t116 * t84;
t184 = t117 * t83;
t183 = t117 * t84;
t163 = qJD(3) * t122;
t151 = t121 * t163;
t39 = t116 * t65 + t117 * t151;
t177 = t119 * t122;
t50 = t116 * t94 + t117 * t177;
t124 = qJD(1) ^ 2;
t180 = t116 * t124;
t178 = t117 * t124;
t123 = qJD(3) ^ 2;
t175 = t123 * t119;
t174 = t123 * t121;
t173 = t124 * qJ(2);
t172 = t114 - t115;
t170 = -t123 - t124;
t165 = qJD(3) * t119;
t157 = qJD(2) * t159;
t156 = t121 * t124 * t119;
t153 = t118 * t165;
t152 = t120 * t165;
t146 = t121 * t160;
t145 = -t116 * t122 + pkin(4);
t144 = pkin(4) * t116 - t122;
t143 = -qJD(4) + t187;
t142 = -t84 + t167;
t141 = -t83 + t166;
t137 = t119 * t146;
t135 = -t116 * t18 + t117 * t19;
t134 = -t116 * t31 - t117 * t30;
t133 = -t116 * t30 + t117 * t31;
t6 = t118 * t15 + t120 * t17;
t81 = t117 * t94;
t32 = -pkin(7) * t179 + t145 * t119 + t81;
t41 = -pkin(7) * t181 + t50;
t9 = -t118 * t41 + t120 * t32;
t10 = t118 * t32 + t120 * t41;
t130 = t141 * t119;
t129 = qJD(1) * t141;
t69 = -t143 - t182;
t128 = t142 * t169;
t88 = t105 * t165;
t51 = -pkin(4) * t140 + t88;
t127 = -t69 + (t105 + t198) * t121;
t14 = t37 * qJD(5) + t188;
t125 = -qJ(4) * t164 + (t143 + t69) * t119;
t2 = -t6 * qJD(5) - t118 * t16 + t120 * t12;
t13 = -t116 * t139 + t117 * t138 - t83 * t161 + t84 * t162;
t113 = t117 ^ 2;
t112 = t116 ^ 2;
t111 = qJ(2) * t157;
t110 = -pkin(4) * t117 - pkin(3);
t103 = 0.2e1 * t137;
t82 = t144 * t121;
t70 = t144 * t165;
t62 = t87 * t121;
t56 = -pkin(4) * t155 + t93;
t53 = t117 * t65;
t49 = -t116 * t177 + t81;
t40 = -pkin(4) * t83 + t69;
t38 = -t116 * t151 + t53;
t28 = pkin(7) * t116 * t165 + t39;
t27 = -t116 * t152 - t117 * t153 + t199 * t121;
t25 = -t116 * t153 + t117 * t152 + t121 * t74;
t20 = t53 + (t145 * t121 + t158) * qJD(3);
t4 = -t10 * qJD(5) - t118 * t28 + t120 * t20;
t3 = t9 * qJD(5) + t118 * t20 + t120 * t28;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, t111, -0.2e1 * t137, 0.2e1 * t172 * t160, -t175, t103, -t174, 0, -t122 * t175 + (qJ(2) * t164 + qJD(2) * t119) * t159, -t122 * t174 + (-qJ(2) * t165 + qJD(2) * t121) * t159, 0, t111, (-t113 * t168 - t183) * t165, (-t184 + (t84 + 0.2e1 * t154) * t116) * t165, (-t117 * t197 + t121 * t84) * qJD(3), (-t112 * t168 + t186) * t165, (t116 * t197 + t121 * t83) * qJD(3), t103, (qJD(1) * t38 + t18) * t119 + ((qJD(1) * t49 + t30) * t121 + (t116 * t127 - t122 * t83) * t119) * qJD(3), (-qJD(1) * t39 - t19) * t119 + ((-qJD(1) * t50 - t31) * t121 + (t117 * t127 + t122 * t84) * t119) * qJD(3), -t38 * t84 + t39 * t83 + (-t116 * t19 - t117 * t18) * t121 + ((t116 * t50 + t117 * t49) * qJD(1) - t134) * t165, t18 * t49 + t19 * t50 + t30 * t38 + t31 * t39 + (t69 - t182) * t119 * t163, t13 * t64 - t25 * t37, t13 * t62 + t14 * t64 + t25 * t34 - t27 * t37, -t25 * t108 - t13 * t119 + (-qJD(1) * t64 + t37) * t164, t14 * t62 + t27 * t34, -t27 * t108 - t14 * t119 + (-qJD(1) * t62 - t34) * t164, (t108 + t169) * t164, t108 * t4 + t119 * t2 + t14 * t82 + t27 * t40 - t34 * t70 + t51 * t62 + (qJD(1) * t9 - t131) * t164, -t1 * t119 - t108 * t3 - t13 * t82 - t25 * t40 - t37 * t70 - t51 * t64 + (-qJD(1) * t10 - t6) * t164, -t1 * t62 - t10 * t14 + t13 * t9 - t131 * t25 + t2 * t64 - t27 * t6 - t3 * t34 - t37 * t4, t1 * t10 - t131 * t4 + t2 * t9 + t3 * t6 - t40 * t70 + t51 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t173, 0, 0, 0, 0, 0, 0, t170 * t119, t170 * t121, 0, -t173, 0, 0, 0, 0, 0, 0, (-t178 + (-t83 - t148) * qJD(3)) * t119, (t180 + (t84 - t154) * qJD(3)) * t119, (t184 + t185) * t164 + (t183 - t186) * qJD(1), t135 * t119 + t134 * qJD(1) + (t119 * t69 + (t133 - t93) * t121) * qJD(3), 0, 0, 0, 0, 0, 0, -t121 * t14 + t191 * t108 + (t119 * t34 - t61 * t168) * qJD(3), t121 * t13 - t192 * t108 + (t119 * t37 + t63 * t168) * qJD(3), -t61 * t13 + t63 * t14 - t191 * t37 - t192 * t34, -t1 * t63 - t121 * t51 - t131 * t191 + t40 * t165 + t192 * t6 - t2 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, -t172 * t124, 0, -t156, 0, 0, -t121 * t173, t119 * t173, 0, 0, -t117 * t128, (-t185 + t184 + (t112 - t113) * qJD(3)) * t169, t114 * t178 + t142 * t168, t119 * t116 * t129, -t114 * t180 + t121 * t129, -t156, -t105 * t130 + (t116 * t125 - t119 * t42 - t121 * t30) * qJD(1), t142 * t93 + (t117 * t125 + t119 * t43 + t121 * t31) * qJD(1), t42 * t84 - t43 * t83 + (qJD(4) * t83 - t30 * t169 + t19) * t117 + (qJD(4) * t84 - t31 * t169 - t18) * t116, -t30 * t42 - t31 * t43 + (-t69 - t187) * t93 + t133 * qJD(4) + t135 * qJ(4), -t13 * t87 + t189 * t37, t13 * t86 - t14 * t87 - t189 * t34 - t190 * t37, t189 * t108 + (qJD(3) * t87 - t37) * t168, t14 * t86 + t190 * t34, -t190 * t108 + (-qJD(3) * t86 + t34) * t168, -t108 * t168, t110 * t14 - t34 * t56 + t51 * t86 + t190 * t40 + t195 * t108 + (qJD(3) * t46 + t131) * t168, -t110 * t13 - t37 * t56 + t51 * t87 + t189 * t40 - t196 * t108 + (-qJD(3) * t47 + t6) * t168, -t1 * t86 + t13 * t46 + t131 * t189 - t14 * t47 - t190 * t6 - t195 * t37 - t196 * t34 - t2 * t87, t1 * t47 + t110 * t51 - t131 * t195 + t196 * t6 + t2 * t46 - t40 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -qJD(1) * t130, -t83 ^ 2 - t84 ^ 2, t30 * t84 - t31 * t83 + t88, 0, 0, 0, 0, 0, 0, t37 * t108 + t14, -t13 - t200, -t201 - t202, -t131 * t37 + t34 * t6 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, t201 - t202, -t13 + t200, -t194, -t188 + (-qJD(5) + t108) * t37, t146, t108 * t6 - t37 * t40 + t2, -t108 * t131 + t34 * t40 - t1, 0, 0;];
tauc_reg = t5;
