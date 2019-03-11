% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:12
% EndTime: 2019-03-08 21:02:20
% DurationCPUTime: 2.92s
% Computational Cost: add. (3588->325), mult. (9609->487), div. (0->0), fcn. (7589->12), ass. (0->181)
t159 = sin(qJ(2));
t154 = sin(pkin(6));
t210 = qJD(1) * t154;
t196 = t159 * t210;
t158 = sin(qJ(3));
t204 = t158 * qJD(3);
t250 = pkin(3) * t204 - t196;
t161 = cos(qJ(3));
t226 = cos(pkin(11));
t190 = t226 * t161;
t142 = qJD(2) * t190;
t153 = sin(pkin(11));
t208 = qJD(2) * t158;
t119 = t153 * t208 - t142;
t114 = qJD(6) + t119;
t132 = t153 * t161 + t226 * t158;
t121 = t132 * qJD(3);
t169 = -t153 * t158 + t190;
t124 = t169 * qJD(3);
t249 = t121 * pkin(4) - t124 * qJ(5) - t132 * qJD(5) + t250;
t241 = -qJ(4) - pkin(8);
t191 = qJD(3) * t241;
t115 = t161 * qJD(4) + t158 * t191;
t116 = -t158 * qJD(4) + t161 * t191;
t162 = cos(qJ(2));
t195 = t162 * t210;
t227 = t226 * t115 + t153 * t116 - t169 * t195;
t248 = qJD(6) - t114;
t152 = sin(pkin(12));
t155 = cos(pkin(12));
t157 = sin(qJ(6));
t160 = cos(qJ(6));
t133 = t160 * t152 + t157 * t155;
t228 = t114 * t133;
t122 = t132 * qJD(2);
t105 = t152 * qJD(3) + t155 * t122;
t106 = t155 * qJD(3) - t152 * t122;
t175 = -t160 * t105 - t157 * t106;
t247 = t114 * t175;
t246 = t160 * t106;
t238 = -t227 * t152 + t155 * t249;
t237 = t152 * t249 + t227 * t155;
t235 = t153 * t115 - t226 * t116 - t132 * t195;
t111 = qJD(2) * t121;
t214 = t160 * t155;
t215 = t157 * t152;
t131 = -t214 + t215;
t229 = t114 * t131;
t245 = -t133 * t111 + t229 * t114;
t203 = qJD(2) * qJD(3);
t192 = t158 * t203;
t112 = qJD(3) * t142 - t153 * t192;
t22 = -t175 * qJD(6) + t133 * t112;
t118 = t119 ^ 2;
t244 = t155 * pkin(9);
t135 = qJD(2) * pkin(8) + t196;
t182 = qJD(4) + t195;
t156 = cos(pkin(6));
t209 = qJD(1) * t156;
t194 = t158 * t209;
t165 = (-t161 * t135 - t194) * qJD(3) + (-qJD(3) * t161 * qJ(4) - t182 * t158) * qJD(2);
t143 = t161 * t209;
t65 = (-t158 * t135 + t143) * qJD(3) + (-qJ(4) * t204 + t182 * t161) * qJD(2);
t23 = t153 * t65 - t226 * t165;
t219 = t154 * t159;
t127 = t156 * t158 + t161 * t219;
t173 = t156 * t161 - t158 * t219;
t79 = t153 * t127 - t226 * t173;
t243 = t23 * t79;
t145 = t153 * pkin(3) + qJ(5);
t242 = pkin(9) + t145;
t240 = t121 * pkin(5) - t124 * t244 + t238;
t224 = t124 * t152;
t239 = pkin(9) * t224 - t237;
t24 = t153 * t165 + t226 * t65;
t20 = qJD(3) * qJD(5) + t24;
t117 = pkin(3) * t192 + qJD(2) * t196;
t37 = t111 * pkin(4) - t112 * qJ(5) - t122 * qJD(5) + t117;
t7 = t152 * t37 + t155 * t20;
t189 = qJ(4) * qJD(2) + t135;
t97 = t189 * t161 + t194;
t193 = t226 * t97;
t96 = -t189 * t158 + t143;
t91 = qJD(3) * pkin(3) + t96;
t45 = t153 * t91 + t193;
t40 = qJD(3) * qJ(5) + t45;
t200 = -t161 * pkin(3) - pkin(2);
t113 = t200 * qJD(2) + qJD(4) - t195;
t62 = t119 * pkin(4) - t122 * qJ(5) + t113;
t14 = t152 * t62 + t155 * t40;
t87 = t153 * t97;
t49 = t226 * t96 - t87;
t74 = pkin(3) * t208 + t122 * pkin(4) + t119 * qJ(5);
t19 = t152 * t74 + t155 * t49;
t236 = pkin(5) * t224 + t235;
t138 = t241 * t158;
t139 = t241 * t161;
t102 = t153 * t138 - t226 * t139;
t83 = -pkin(4) * t169 - t132 * qJ(5) + t200;
t42 = t155 * t102 + t152 * t83;
t234 = qJD(2) * pkin(2);
t53 = t157 * t105 - t246;
t233 = t122 * t53;
t101 = -t226 * t138 - t153 * t139;
t232 = t23 * t101;
t231 = t175 * t122;
t205 = qJD(6) * t160;
t230 = t106 * t205 + t112 * t214;
t225 = t119 * t152;
t223 = t132 * t152;
t222 = t132 * t155;
t220 = t152 * t112;
t218 = t154 * t162;
t164 = qJD(2) ^ 2;
t217 = t154 * t164;
t216 = t155 * t112;
t163 = qJD(3) ^ 2;
t213 = t163 * t158;
t212 = t163 * t161;
t211 = t158 ^ 2 - t161 ^ 2;
t207 = qJD(2) * t159;
t206 = qJD(6) * t132;
t201 = t159 * t217;
t6 = -t152 * t20 + t155 * t37;
t4 = t111 * pkin(5) - pkin(9) * t216 + t6;
t5 = -pkin(9) * t220 + t7;
t199 = -t157 * t5 + t160 * t4;
t198 = t154 * t207;
t197 = qJD(2) * t218;
t13 = -t152 * t40 + t155 * t62;
t18 = -t152 * t49 + t155 * t74;
t47 = t153 * t96 + t193;
t41 = -t152 * t102 + t155 * t83;
t188 = t158 * t197;
t187 = t161 * t197;
t186 = -t131 * t111 - t228 * t114;
t148 = -t226 * pkin(3) - pkin(4);
t185 = t157 * t4 + t160 * t5;
t8 = t119 * pkin(5) - t105 * pkin(9) + t13;
t9 = pkin(9) * t106 + t14;
t184 = t157 * t9 - t160 * t8;
t2 = t157 * t8 + t160 * t9;
t44 = t226 * t91 - t87;
t181 = -t13 * t152 + t14 * t155;
t27 = -pkin(5) * t169 - pkin(9) * t222 + t41;
t29 = -pkin(9) * t223 + t42;
t180 = -t157 * t29 + t160 * t27;
t179 = t157 * t27 + t160 * t29;
t80 = t226 * t127 + t153 * t173;
t63 = -t152 * t80 - t155 * t218;
t64 = -t152 * t218 + t155 * t80;
t178 = -t157 * t64 + t160 * t63;
t177 = t157 * t63 + t160 * t64;
t176 = t101 * t112 + t23 * t132;
t174 = -qJD(6) * t105 - t220;
t129 = t242 * t155;
t172 = t122 * pkin(5) + qJD(5) * t152 + qJD(6) * t129 + t119 * t244 + t18;
t128 = t242 * t152;
t171 = pkin(9) * t225 - qJD(5) * t155 + qJD(6) * t128 + t19;
t170 = t234 * qJD(2);
t39 = -qJD(3) * pkin(4) + qJD(5) - t44;
t168 = t124 * t39 + t176;
t167 = -0.2e1 * qJD(3) * t234;
t21 = t174 * t157 + t230;
t166 = -t111 * t145 + t112 * t148 + (-qJD(5) + t39) * t119;
t137 = -t155 * pkin(5) + t148;
t95 = -t127 * qJD(3) - t188;
t94 = t173 * qJD(3) + t187;
t78 = t131 * t132;
t77 = t133 * t132;
t71 = pkin(5) * t223 + t101;
t48 = t153 * t95 + t226 * t94;
t46 = t153 * t94 - t226 * t95;
t34 = t152 * t198 + t155 * t48;
t33 = -t152 * t48 + t155 * t198;
t32 = -pkin(5) * t225 + t47;
t31 = t133 * t124 + t205 * t222 - t206 * t215;
t30 = -t131 * t124 - t133 * t206;
t28 = -pkin(5) * t106 + t39;
t15 = pkin(5) * t220 + t23;
t1 = [0, 0, -t201, -t162 * t217, 0, 0, 0, 0, 0, -t161 * t201 + (t95 - t188) * qJD(3), t158 * t201 + (-t94 - t187) * qJD(3), -t80 * t111 + t79 * t112 - t48 * t119 + t46 * t122, t243 + t24 * t80 - t44 * t46 + t45 * t48 + (t113 * t207 - t117 * t162) * t154, -t106 * t46 + t63 * t111 + t33 * t119 + t79 * t220, t46 * t105 - t64 * t111 - t34 * t119 + t79 * t216, t34 * t106 - t33 * t105 + (-t152 * t64 - t155 * t63) * t112, t13 * t33 + t14 * t34 + t39 * t46 + t6 * t63 + t7 * t64 + t243, 0, 0, 0, 0, 0 (-qJD(6) * t177 - t157 * t34 + t160 * t33) * t114 + t178 * t111 + t46 * t53 + t79 * t22 -(qJD(6) * t178 + t157 * t33 + t160 * t34) * t114 - t177 * t111 - t46 * t175 + t79 * t21; 0, 0, 0, 0, 0.2e1 * t161 * t192, -0.2e1 * t211 * t203, t212, -t213, 0, -pkin(8) * t212 + t158 * t167, pkin(8) * t213 + t161 * t167, -t102 * t111 - t227 * t119 - t45 * t121 + t235 * t122 - t44 * t124 + t169 * t24 + t176, t24 * t102 + t113 * t250 + t117 * t200 + t227 * t45 - t235 * t44 + t232, -t106 * t235 + t41 * t111 + t238 * t119 + t13 * t121 + t168 * t152 - t169 * t6, t235 * t105 - t42 * t111 - t237 * t119 - t14 * t121 + t168 * t155 + t169 * t7, -t238 * t105 + t237 * t106 + (-t112 * t41 - t124 * t13 - t132 * t6) * t155 + (-t112 * t42 - t124 * t14 - t132 * t7) * t152, t238 * t13 + t237 * t14 + t235 * t39 + t6 * t41 + t7 * t42 + t232, -t175 * t30 - t21 * t78, t175 * t31 - t21 * t77 + t78 * t22 - t30 * t53, -t78 * t111 + t30 * t114 - t121 * t175 - t169 * t21, -t77 * t111 - t31 * t114 - t53 * t121 + t169 * t22, -t111 * t169 + t114 * t121, t180 * t111 - t199 * t169 - t184 * t121 + t71 * t22 + t15 * t77 + t28 * t31 + t236 * t53 + (t239 * t157 + t240 * t160) * t114 + (-t114 * t179 + t169 * t2) * qJD(6), -t179 * t111 + t185 * t169 - t2 * t121 + t71 * t21 - t15 * t78 + t28 * t30 - t236 * t175 + (-t240 * t157 + t239 * t160) * t114 + (-t114 * t180 - t169 * t184) * qJD(6); 0, 0, 0, 0, -t158 * t164 * t161, t211 * t164, 0, 0, 0, t158 * t170, t161 * t170 (t45 - t47) * t122 + (-t44 + t49) * t119 + (-t111 * t153 - t226 * t112) * pkin(3), t44 * t47 - t45 * t49 + (-t113 * t208 + t153 * t24 - t226 * t23) * pkin(3), t106 * t47 - t18 * t119 - t13 * t122 + t166 * t152 - t23 * t155, -t47 * t105 + t19 * t119 + t14 * t122 + t23 * t152 + t166 * t155, -t19 * t106 + t18 * t105 + (qJD(5) * t106 - t119 * t13 + t7) * t155 + (qJD(5) * t105 - t119 * t14 - t6) * t152, -t13 * t18 - t14 * t19 + t23 * t148 - t39 * t47 + (-t6 * t152 + t7 * t155) * t145 + t181 * qJD(5), t21 * t133 + t175 * t229, -t21 * t131 - t133 * t22 + t175 * t228 + t229 * t53, t231 - t245, t186 + t233, -t114 * t122 (-t160 * t128 - t157 * t129) * t111 + t137 * t22 + t15 * t131 + t184 * t122 - t32 * t53 + t228 * t28 + (t157 * t171 - t160 * t172) * t114 -(-t157 * t128 + t160 * t129) * t111 + t137 * t21 + t15 * t133 + t2 * t122 + t32 * t175 - t229 * t28 + (t157 * t172 + t160 * t171) * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122 ^ 2 - t118, t45 * t119 + t44 * t122 + t117, t106 * t122 + t155 * t111 - t152 * t118, -t122 * t105 - t152 * t111 - t155 * t118 (t105 * t152 + t106 * t155) * t119 + (-t152 ^ 2 - t155 ^ 2) * t112, t119 * t181 - t39 * t122 + t7 * t152 + t6 * t155, 0, 0, 0, 0, 0, t186 - t233, t231 + t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t119 + t220, t106 * t119 + t216, -t105 ^ 2 - t106 ^ 2, t105 * t13 - t14 * t106 + t23, 0, 0, 0, 0, 0, t22 - t247, t114 * t246 + (-t105 * t114 + t174) * t157 + t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175 * t53, t175 ^ 2 - t53 ^ 2, t53 * t114 + t21, -t22 - t247, t111, t28 * t175 - t2 * t248 + t199, t184 * t248 + t28 * t53 - t185;];
tauc_reg  = t1;
