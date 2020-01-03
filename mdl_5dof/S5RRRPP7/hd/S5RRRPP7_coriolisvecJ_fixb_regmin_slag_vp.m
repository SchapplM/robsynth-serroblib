% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:54
% EndTime: 2019-12-31 21:06:01
% DurationCPUTime: 2.34s
% Computational Cost: add. (2035->334), mult. (4946->437), div. (0->0), fcn. (2877->4), ass. (0->165)
t103 = cos(qJ(3));
t102 = sin(qJ(2));
t164 = qJD(1) * t102;
t146 = t103 * t164;
t101 = sin(qJ(3));
t162 = qJD(2) * t101;
t67 = t146 + t162;
t104 = cos(qJ(2));
t163 = qJD(1) * t104;
t88 = -qJD(3) + t163;
t196 = t67 * t88;
t158 = qJD(3) * t103;
t143 = t102 * t158;
t152 = qJD(2) * qJD(3);
t160 = qJD(2) * t104;
t40 = qJD(1) * (t101 * t160 + t143) + t101 * t152;
t215 = -t40 + t196;
t153 = qJD(1) * qJD(2);
t140 = t102 * t153;
t214 = qJ(4) * t140 - t88 * qJD(4);
t203 = pkin(3) + pkin(4);
t213 = t203 * t140;
t155 = t103 * qJD(2);
t65 = t101 * t164 - t155;
t212 = t40 * qJ(5) + t65 * qJD(5);
t211 = 0.2e1 * t214;
t74 = -pkin(2) * t104 - pkin(7) * t102 - pkin(1);
t55 = t74 * qJD(1);
t96 = pkin(6) * t163;
t80 = qJD(2) * pkin(7) + t96;
t28 = t101 * t55 + t103 * t80;
t83 = t88 * qJ(4);
t22 = -t83 + t28;
t132 = pkin(3) * t140;
t133 = pkin(6) * t140;
t159 = qJD(3) * t101;
t129 = pkin(2) * t102 - pkin(7) * t104;
t72 = t129 * qJD(2);
t58 = qJD(1) * t72;
t138 = -t101 * t133 - t103 * t58 + t80 * t158 + t55 * t159;
t8 = -t132 + t138;
t210 = t22 * t88 + t8;
t63 = t67 ^ 2;
t209 = -t88 ^ 2 - t63;
t79 = -qJD(2) * pkin(2) + pkin(6) * t164;
t121 = t67 * qJ(4) - t79;
t13 = -t203 * t65 + qJD(5) + t121;
t208 = (qJD(5) + t13) * t67;
t207 = -t101 * qJD(4) - t96;
t198 = t65 * t88;
t142 = t104 * t153;
t145 = t102 * t159;
t39 = qJD(1) * t145 + (-t142 - t152) * t103;
t206 = -t39 + t198;
t27 = -t101 * t80 + t103 * t55;
t166 = qJD(4) - t27;
t175 = t101 * qJ(4);
t205 = -t203 * t103 - t175;
t204 = t65 ^ 2;
t202 = pkin(6) * t101;
t17 = qJ(5) * t65 + t28;
t12 = t17 - t83;
t201 = t12 * t88;
t24 = pkin(3) * t65 - t121;
t199 = t24 * t67;
t197 = t67 * t65;
t7 = t40 * pkin(3) + pkin(6) * t142 + t39 * qJ(4) - t67 * qJD(4);
t195 = t7 * t101;
t194 = t7 * t103;
t193 = pkin(7) - qJ(5);
t147 = -pkin(3) - t202;
t134 = -pkin(4) + t147;
t171 = t103 * t104;
t71 = t129 * qJD(1);
t182 = t103 * t71;
t78 = t193 * t103;
t192 = -t182 + (-qJ(5) * t171 + t102 * t134) * qJD(1) - qJD(3) * t78 + t101 * qJD(5);
t154 = t103 * qJD(5);
t172 = t102 * t103;
t173 = t101 * t104;
t53 = t101 * t71;
t188 = qJ(4) * t164 + t53;
t191 = (-pkin(6) * t172 + qJ(5) * t173) * qJD(1) + t188 + t159 * t193 + t154;
t176 = qJ(4) * t103;
t116 = -t203 * t101 + t176;
t190 = t88 * t116 + t207;
t127 = pkin(3) * t101 - t176;
t189 = t88 * t127 - t207;
t187 = t101 * t72 + t74 * t158;
t91 = pkin(6) * t171;
t186 = qJD(3) * t91 + t74 * t159;
t185 = t101 * t74 + t91;
t184 = qJ(4) * t40;
t183 = t101 * t79;
t181 = t103 * t79;
t180 = t103 * t88;
t179 = t39 * t101;
t178 = t65 * qJ(4);
t99 = t102 ^ 2;
t177 = t104 ^ 2 - t99;
t174 = t101 * t102;
t107 = qJD(1) ^ 2;
t170 = t104 * t107;
t106 = qJD(2) ^ 2;
t169 = t106 * t102;
t168 = t106 * t104;
t16 = t67 * qJ(5) + t27;
t167 = qJD(4) - t16;
t161 = qJD(2) * t102;
t157 = qJD(4) * t103;
t151 = pkin(7) * t101 * t88;
t150 = pkin(7) * t180;
t149 = pkin(7) * t161;
t148 = pkin(7) * t155;
t144 = t104 * t159;
t90 = pkin(6) * t173;
t139 = t103 * t74 - t90;
t137 = -0.2e1 * pkin(1) * t153;
t136 = -t67 + t162;
t135 = t65 + t155;
t131 = t103 * t72 - t186;
t37 = -qJ(4) * t104 + t185;
t130 = t147 * t102;
t128 = pkin(3) * t103 + t175;
t21 = pkin(3) * t88 + t166;
t126 = -t101 * t22 + t103 * t21;
t125 = qJD(1) * t99 - t104 * t88;
t124 = t101 * t58 + t55 * t158 - t159 * t80;
t123 = qJ(4) * t161 - t104 * qJD(4) + t187;
t122 = pkin(6) + t127;
t120 = -t28 * t88 - t138;
t3 = -pkin(4) * t40 - t7;
t119 = -t3 * t101 - t13 * t158;
t118 = t3 * t103 - t13 * t159;
t117 = t39 * qJ(5) + t138;
t115 = -t140 + t197;
t114 = -pkin(6) + t116;
t112 = t39 + t198;
t111 = t103 * t133 - t124;
t110 = t117 - t213;
t6 = -t111 + t214;
t109 = -t27 * t88 + t111;
t98 = t104 * pkin(3);
t77 = t193 * t101;
t73 = -pkin(2) - t128;
t61 = pkin(2) - t205;
t45 = t122 * t102;
t38 = -t139 + t98;
t36 = t114 * t102;
t32 = pkin(3) * t67 + t178;
t31 = qJD(1) * t130 - t182;
t30 = -pkin(6) * t146 + t188;
t26 = qJ(5) * t174 + t37;
t25 = t104 * pkin(4) + t90 + t98 + (-qJ(5) * t102 - t74) * t103;
t19 = -t203 * t67 - t178;
t15 = (qJD(3) * t128 - t157) * t102 + t122 * t160;
t14 = qJD(2) * t130 - t131;
t11 = (-t102 * t155 - t144) * pkin(6) + t123;
t10 = (t205 * qJD(3) + t157) * t102 + t114 * t160;
t9 = t203 * t88 + t167;
t5 = (-pkin(6) * qJD(2) + qJ(5) * qJD(3)) * t172 + (qJD(5) * t102 + (-pkin(6) * qJD(3) + qJ(5) * qJD(2)) * t104) * t101 + t123;
t4 = (-qJ(5) * t160 - t72) * t103 + (qJ(5) * t159 + t134 * qJD(2) - t154) * t102 + t186;
t2 = -t67 * qJD(5) + t110;
t1 = t6 + t212;
t18 = [0, 0, 0, 0.2e1 * t104 * t140, 0.2e1 * t177 * t153, t168, -t169, 0, -pkin(6) * t168 + t102 * t137, pkin(6) * t169 + t104 * t137, t67 * t104 * t155 + (-t39 * t103 - t159 * t67) * t102, (-t101 * t67 - t103 * t65) * t160 + (t179 - t103 * t40 + (t101 * t65 - t103 * t67) * qJD(3)) * t102, t88 * t145 + t39 * t104 + (t102 * t67 + t103 * t125) * qJD(2), t88 * t143 + t40 * t104 + (-t101 * t125 - t102 * t65) * qJD(2), (-t88 - t163) * t161, -t131 * t88 + t138 * t104 + (pkin(6) * t40 + t158 * t79) * t102 + ((pkin(6) * t65 + t183) * t104 + (t139 * qJD(1) + t27 + (-t88 + t163) * t202) * t102) * qJD(2), (-pkin(6) * t144 + t187) * t88 + t124 * t104 + (-pkin(6) * t39 - t159 * t79) * t102 + ((pkin(6) * t67 + t181) * t104 + (-pkin(6) * t180 - qJD(1) * t185 - t28) * t102) * qJD(2), t14 * t88 + t15 * t65 + t45 * t40 + (t162 * t24 + t8) * t104 + (t24 * t158 + t195 + (-qJD(1) * t38 - t21) * qJD(2)) * t102, -t11 * t65 + t14 * t67 - t37 * t40 - t38 * t39 + t126 * t160 + (-t101 * t6 + t103 * t8 + (-t101 * t21 - t103 * t22) * qJD(3)) * t102, -t11 * t88 - t15 * t67 + t45 * t39 + (-t155 * t24 - t6) * t104 + (t24 * t159 - t194 + (qJD(1) * t37 + t22) * qJD(2)) * t102, t11 * t22 + t14 * t21 + t15 * t24 + t37 * t6 + t38 * t8 + t45 * t7, -t10 * t65 - t36 * t40 + t4 * t88 + (-t13 * t162 + t2) * t104 + ((-qJD(1) * t25 - t9) * qJD(2) + t119) * t102, t10 * t67 - t36 * t39 - t5 * t88 + (t13 * t155 - t1) * t104 + ((qJD(1) * t26 + t12) * qJD(2) + t118) * t102, t25 * t39 + t26 * t40 - t4 * t67 + t5 * t65 + (t101 * t12 - t103 * t9) * t160 + (t1 * t101 - t103 * t2 + (t101 * t9 + t103 * t12) * qJD(3)) * t102, t1 * t26 + t10 * t13 + t12 * t5 + t2 * t25 + t3 * t36 + t4 * t9; 0, 0, 0, -t102 * t170, -t177 * t107, 0, 0, 0, t107 * pkin(1) * t102, pkin(1) * t170, -t180 * t67 - t179, t215 * t101 + t206 * t103, -t88 * t158 + (t102 * t136 + t171 * t88) * qJD(1), t88 * t159 + (t102 * t135 - t173 * t88) * qJD(1), t88 * t164, t71 * t180 - pkin(2) * t40 + (t150 + t183) * qJD(3) + (-t27 * t102 + (-t104 * t79 - t149) * t101 + (-t104 * t135 + t174 * t88) * pkin(6)) * qJD(1), pkin(2) * t39 - t53 * t88 + (-t151 + t181) * qJD(3) + (-t79 * t171 + (t28 - t148) * t102 + (t104 * t136 + t172 * t88) * pkin(6)) * qJD(1), -t194 - t31 * t88 + t73 * t40 - t189 * t65 + (t101 * t24 + t150) * qJD(3) + (t102 * t21 + (-t104 * t24 - t149) * t101) * qJD(1), t30 * t65 - t31 * t67 + (t6 - t88 * t21 + (qJD(3) * t67 - t40) * pkin(7)) * t103 + ((qJD(3) * t65 - t39) * pkin(7) + t210) * t101, -t195 + t30 * t88 + t39 * t73 + t189 * t67 + (-t103 * t24 + t151) * qJD(3) + (t24 * t171 + (-t22 + t148) * t102) * qJD(1), -t21 * t31 - t22 * t30 + t7 * t73 - t189 * t24 + (qJD(3) * t126 + t8 * t101 + t6 * t103) * pkin(7), -t40 * t61 - t192 * t88 + t190 * t65 + (t13 * t173 + (-qJD(2) * t77 + t9) * t102) * qJD(1) + t118, -t39 * t61 + t191 * t88 - t190 * t67 + (-t13 * t171 + (qJD(2) * t78 - t12) * t102) * qJD(1) - t119, t39 * t77 + t40 * t78 + t192 * t67 - t191 * t65 + (t88 * t9 - t1) * t103 + (-t2 - t201) * t101, t1 * t78 - t12 * t191 - t13 * t190 - t192 * t9 + t2 * t77 + t3 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t63 - t204, -t112, -t40 - t196, t140, -t67 * t79 + t120, t65 * t79 + t109, -t32 * t65 + t120 + 0.2e1 * t132 - t199, pkin(3) * t39 - t184 + (t22 - t28) * t67 + (t21 - t166) * t65, -t24 * t65 + t32 * t67 - t109 + t211, -t8 * pkin(3) + t6 * qJ(4) + t166 * t22 - t21 * t28 - t24 * t32, -t17 * t88 + t19 * t65 - t117 + t208 + 0.2e1 * t213, t13 * t65 + t16 * t88 - t19 * t67 - t111 + t211 + t212, t184 - t203 * t39 + (-t12 + t17) * t67 + (-t9 + t167) * t65, t1 * qJ(4) + t12 * t167 - t13 * t19 - t9 * t17 - t2 * t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, -t112, t209, t199 + t210, t115, t209, t112, t110 + t201 - t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, t206, -t63 - t204, -t12 * t65 + t67 * t9 + t3;];
tauc_reg = t18;
