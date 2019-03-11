% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:03
% EndTime: 2019-03-08 19:45:11
% DurationCPUTime: 2.81s
% Computational Cost: add. (4277->331), mult. (10985->442), div. (0->0), fcn. (8631->10), ass. (0->182)
t118 = sin(pkin(11));
t120 = cos(pkin(11));
t123 = sin(qJ(4));
t221 = cos(qJ(4));
t99 = t221 * t118 + t123 * t120;
t136 = qJD(2) * t99;
t241 = qJD(6) + t136;
t122 = sin(qJ(6));
t125 = cos(qJ(6));
t175 = t221 * t120;
t161 = qJD(2) * t175;
t184 = qJD(2) * t118;
t170 = t123 * t184;
t88 = -t161 + t170;
t65 = qJD(4) * t122 - t125 * t88;
t168 = t241 * t65;
t180 = qJD(6) * t125;
t181 = qJD(6) * t122;
t93 = t99 * qJD(4);
t80 = qJD(2) * t93;
t42 = qJD(4) * t181 - t122 * t80 - t88 * t180;
t243 = t42 - t168;
t119 = sin(pkin(6));
t126 = cos(qJ(2));
t197 = t119 * t126;
t134 = t99 * t197;
t212 = pkin(8) + qJ(3);
t104 = t212 * t118;
t105 = t212 * t120;
t62 = -t123 * t104 + t221 * t105;
t208 = -qJD(1) * t134 + qJD(3) * t99 + qJD(4) * t62;
t124 = sin(qJ(2));
t188 = qJD(1) * t119;
t173 = t124 * t188;
t169 = qJD(4) * t221;
t182 = qJD(4) * t123;
t92 = t118 * t182 - t120 * t169;
t247 = -t92 * qJ(5) + t99 * qJD(5) + t173;
t186 = qJD(1) * t126;
t171 = t119 * t186;
t97 = (qJD(3) + t171) * qJD(2);
t137 = t99 * t97;
t103 = qJD(2) * qJ(3) + t173;
t121 = cos(pkin(6));
t187 = qJD(1) * t121;
t109 = t120 * t187;
t207 = pkin(8) * qJD(2);
t59 = t109 + (-t103 - t207) * t118;
t69 = t120 * t103 + t118 * t187;
t60 = t120 * t207 + t69;
t32 = t123 * t59 + t221 * t60;
t16 = t32 * qJD(4) + t137;
t107 = qJD(4) * t161;
t79 = qJD(4) * t170 - t107;
t13 = -t79 * pkin(5) + t16;
t183 = qJD(2) * t124;
t172 = t119 * t183;
t106 = qJD(1) * t172;
t165 = t79 * qJ(5) + t106;
t141 = -qJD(5) * t136 + t165;
t226 = pkin(4) + pkin(9);
t20 = t226 * t80 + t141;
t31 = t123 * t60 - t221 * t59;
t192 = -qJD(5) - t31;
t191 = pkin(5) * t136 - t192;
t19 = -t226 * qJD(4) + t191;
t153 = qJD(3) - t171;
t114 = -t120 * pkin(3) - pkin(2);
t185 = qJD(2) * t114;
t81 = t153 + t185;
t130 = -qJ(5) * t136 + t81;
t24 = t226 * t88 + t130;
t6 = t122 * t19 + t125 * t24;
t2 = -qJD(6) * t6 - t122 * t20 + t125 * t13;
t231 = t241 * t6 + t2;
t150 = t122 * t24 - t125 * t19;
t1 = -qJD(6) * t150 + t122 * t13 + t125 * t20;
t159 = t150 * t241 + t1;
t246 = -t92 * pkin(5) + t208;
t235 = t125 * t241;
t67 = qJD(4) * t125 + t122 * t88;
t245 = t67 * t235;
t244 = -t226 * t93 + t247;
t164 = t122 * t241;
t145 = -t125 * t79 - t164 * t241;
t195 = t123 * t118;
t98 = -t175 + t195;
t242 = t136 * qJD(4);
t227 = t88 ^ 2;
t86 = t136 ^ 2;
t239 = -t227 - t86;
t238 = -t227 + t86;
t58 = t79 * t99;
t237 = -t136 * t92 - t58;
t152 = t118 * (-t103 * t118 + t109) - t120 * t69;
t234 = t152 * t126;
t233 = t208 * qJD(4);
t148 = t175 * t197;
t209 = qJD(1) * t148 + (qJD(3) * t118 + qJD(4) * t105) * t123 - t171 * t195 - qJD(3) * t175 + t104 * t169;
t232 = t209 * qJD(4);
t224 = t88 * pkin(5);
t26 = -qJD(4) * qJ(5) - t32;
t21 = -t26 - t224;
t230 = t21 * t241 + t226 * t79;
t198 = t119 * t124;
t140 = t118 * t198 - t120 * t121;
t87 = t118 * t121 + t120 * t198;
t51 = -t123 * t140 + t221 * t87;
t28 = qJD(2) * t134 + qJD(4) * t51;
t229 = t119 * (-t126 * t80 + t88 * t183) - t28 * qJD(4);
t135 = t221 * t140;
t27 = -qJD(2) * t148 + qJD(4) * t135 + (qJD(4) * t87 + t184 * t197) * t123;
t228 = t119 * (-t126 * t79 - t136 * t183) - t27 * qJD(4);
t225 = t80 * pkin(4);
t146 = -t99 * qJ(5) + t114;
t41 = t226 * t98 + t146;
t61 = t221 * t104 + t105 * t123;
t46 = pkin(5) * t99 + t61;
t17 = -t122 * t41 + t125 * t46;
t223 = qJD(6) * t17 + t246 * t122 - t244 * t125;
t18 = t122 * t46 + t125 * t41;
t222 = -qJD(6) * t18 + t244 * t122 + t246 * t125;
t50 = t123 * t87 + t135;
t220 = t16 * t50;
t219 = t16 * t61;
t218 = t67 * t65;
t217 = t67 * t88;
t216 = t79 * t98;
t215 = t88 * t65;
t214 = t88 * t136;
t211 = -t26 - t32;
t210 = -t93 * pkin(5) - t209;
t206 = qJD(2) * pkin(2);
t205 = t125 * t42;
t72 = t125 * t80;
t43 = qJD(6) * t67 - t72;
t203 = t43 * t122;
t202 = t88 * qJ(5);
t200 = qJD(4) * t92;
t199 = qJD(4) * t93;
t128 = qJD(2) ^ 2;
t196 = t119 * t128;
t189 = t118 ^ 2 + t120 ^ 2;
t179 = t98 * t181;
t178 = t98 * t180;
t177 = t124 * t196;
t176 = t126 * t196;
t174 = t119 ^ 2 * t186;
t167 = t189 * t97;
t166 = -t59 * t169 + t60 * t182 + t97 * t98;
t14 = -qJD(4) * qJD(5) + t166;
t11 = -pkin(5) * t80 - t14;
t160 = qJD(6) * t226 * t241 + t11;
t158 = -pkin(4) * t93 + t247;
t157 = t11 * t98 + t21 * t93;
t156 = t122 * t150 + t125 * t6;
t155 = t241 * t93 - t216;
t154 = t80 * t98 + t88 * t93;
t142 = -t122 * t50 + t125 * t197;
t37 = t122 * t197 + t125 * t50;
t138 = t122 * t79 - t235 * t241;
t133 = t136 * t28 + t27 * t88 - t50 * t79 - t51 * t80;
t132 = -t136 * t93 - t80 * t99 + t88 * t92 + t216;
t35 = t88 * pkin(4) + t130;
t131 = t136 * t35 + t137;
t129 = t136 * t208 + t16 * t99 + t209 * t88 - t61 * t79 - t62 * t80;
t102 = t153 - t206;
t82 = qJD(4) * t88;
t54 = pkin(4) * t98 + t146;
t53 = pkin(4) * t136 + t202;
t52 = t79 - t82;
t47 = -t98 * pkin(5) + t62;
t39 = t125 * t43;
t36 = t136 * t226 + t202;
t30 = t141 + t225;
t25 = -qJD(4) * pkin(4) - t192;
t23 = t32 - t224;
t10 = qJD(6) * t37 + t122 * t28 + t125 * t172;
t9 = qJD(6) * t142 - t122 * t172 + t125 * t28;
t8 = t122 * t23 + t125 * t36;
t7 = -t122 * t36 + t125 * t23;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, -t176, 0, 0, 0, 0, 0, 0, 0, 0, -t120 * t177, t118 * t177, t189 * t176 (t118 * t140 + t120 * t87) * t97 + (-t124 * t174 + (t102 * t124 - t234) * t119) * qJD(2), 0, 0, 0, 0, 0, 0, t229, -t228, t133, -t166 * t51 + t220 - t32 * t27 + t31 * t28 + (t119 * t81 - t174) * t183, 0, 0, 0, 0, 0, 0, t133, -t229, t228, -t14 * t51 + t220 + t25 * t28 + t26 * t27 + (-t126 * t30 + t35 * t183) * t119, 0, 0, 0, 0, 0, 0, t241 * t9 - t27 * t65 - t37 * t79 + t43 * t51, -t10 * t241 - t142 * t79 - t27 * t67 - t42 * t51, -t10 * t65 + t142 * t43 + t37 * t42 - t67 * t9, -t1 * t142 + t10 * t6 + t11 * t51 - t150 * t9 + t2 * t37 - t21 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t153 * t189 + t167, -t152 * qJD(3) + qJ(3) * t167 + (t234 + (-t102 - t206) * t124) * t188, t237, t132, -t200, t154, -t199, 0, t114 * t80 + t81 * t93 - t233 + (qJD(2) * t98 - t88) * t173, -t114 * t79 - t81 * t92 + t232, t166 * t98 - t31 * t92 - t32 * t93 + t129, -t166 * t62 + t219 - t209 * t32 + t208 * t31 + (-t81 + t185) * t173, 0, t200, t199, t237, t132, t154, t14 * t98 - t25 * t92 + t26 * t93 + t129, t158 * t88 - t30 * t98 - t35 * t93 - t54 * t80 + t233, t136 * t158 - t30 * t99 + t35 * t92 + t54 * t79 - t232, -t14 * t62 - t158 * t35 + t208 * t25 + t209 * t26 + t30 * t54 + t219, t67 * t178 + (-t42 * t98 + t67 * t93) * t122 (-t122 * t65 + t125 * t67) * t93 + (-t203 - t205 + (-t122 * t67 - t125 * t65) * qJD(6)) * t98, t122 * t155 + t178 * t241 - t42 * t99 - t67 * t92, t65 * t179 + (-t43 * t98 - t65 * t93) * t125, t125 * t155 - t179 * t241 - t43 * t99 + t65 * t92, -t241 * t92 - t58, -t157 * t125 + t150 * t92 - t17 * t79 + t21 * t179 + t2 * t99 + t210 * t65 + t222 * t241 + t47 * t43, -t1 * t99 + t157 * t122 + t21 * t178 + t18 * t79 + t210 * t67 - t223 * t241 - t47 * t42 + t6 * t92, t17 * t42 - t18 * t43 + t156 * t93 - t222 * t67 - t223 * t65 + (t1 * t125 - t2 * t122 + (-t122 * t6 + t125 * t150) * qJD(6)) * t98, t1 * t18 + t11 * t47 - t150 * t222 + t2 * t17 + t210 * t21 + t223 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189 * t128, qJD(2) * t152 + t106, 0, 0, 0, 0, 0, 0, 0.2e1 * t242, t107 + (-t88 - t170) * qJD(4), t239, -t136 * t31 + t32 * t88 + t106, 0, 0, 0, 0, 0, 0, t239, -0.2e1 * t242, t79 + t82, t225 - t26 * t88 + (-qJD(5) - t25) * t136 + t165, 0, 0, 0, 0, 0, 0, t138 + t215, t217 - t145, -t243 * t122 + t245 - t39, -t122 * t231 + t159 * t125 + t21 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, t238, t107 + (t88 - t170) * qJD(4), -t214, 0, 0, -t136 * t81 - t137, -qJD(4) * t31 + t81 * t88 + t166, 0, 0, 0, t52, 0, t214, t238, -t214, pkin(4) * t79 - qJ(5) * t80 + t211 * t136 + (t25 + t192) * t88, t53 * t88 + t131, -t35 * t88 + t53 * t136 + (0.2e1 * qJD(5) + t31) * qJD(4) - t166, -t16 * pkin(4) - t14 * qJ(5) + t192 * t26 - t25 * t32 - t35 * t53, -t67 * t164 - t205, -t39 - t245 + (t42 + t168) * t122, t145 + t217, t235 * t65 + t203, t138 - t215, t241 * t88, qJ(5) * t43 + t160 * t122 + t125 * t230 - t150 * t88 + t191 * t65 - t241 * t7, -qJ(5) * t42 - t122 * t230 + t160 * t125 + t191 * t67 + t241 * t8 - t6 * t88, t8 * t65 + t7 * t67 + (-t226 * t42 - t6 * t136 - t2 + (t226 * t65 - t6) * qJD(6)) * t125 + (t226 * t43 - t150 * t136 - t1 + (-t226 * t67 - t150) * qJD(6)) * t122, t11 * qJ(5) + t150 * t7 - t6 * t8 + t191 * t21 - (qJD(6) * t156 + t1 * t122 + t2 * t125) * t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t214, -qJD(4) ^ 2 - t86, -t211 * qJD(4) + t131, 0, 0, 0, 0, 0, 0, -qJD(4) * t65 + t145, -qJD(4) * t67 + t138, t243 * t125 + (t241 * t67 - t43) * t122, -t21 * qJD(4) + t159 * t122 + t125 * t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, -t65 ^ 2 + t67 ^ 2, -t243, -t218, t72 + (-qJD(6) + t241) * t67, -t79, -t21 * t67 + t231, t21 * t65 - t159, 0, 0;];
tauc_reg  = t3;
