% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:06
% EndTime: 2019-03-08 20:03:14
% DurationCPUTime: 3.30s
% Computational Cost: add. (3863->373), mult. (9764->494), div. (0->0), fcn. (7348->10), ass. (0->203)
t125 = sin(qJ(4));
t128 = cos(qJ(4));
t167 = pkin(4) * t125 - pkin(9) * t128;
t103 = t167 * qJD(4);
t121 = sin(pkin(6));
t120 = sin(pkin(11));
t122 = cos(pkin(11));
t126 = sin(qJ(2));
t129 = cos(qJ(2));
t156 = t120 * t129 + t122 * t126;
t86 = t156 * t121;
t77 = qJD(1) * t86;
t264 = t103 - t77;
t124 = sin(qJ(5));
t127 = cos(qJ(5));
t200 = t128 * qJD(2);
t110 = -qJD(5) + t200;
t202 = t124 * qJD(4);
t207 = qJD(2) * t125;
t99 = t127 * t207 + t202;
t219 = t99 * t110;
t180 = t124 * t207;
t201 = t127 * qJD(4);
t97 = t180 - t201;
t220 = t97 * t110;
t199 = qJD(2) * qJD(4);
t179 = t127 * t199;
t198 = qJD(4) * qJD(5);
t73 = qJD(5) * t180 - t127 * t198 - t128 * t179;
t203 = qJD(5) * t127;
t184 = t125 * t203;
t188 = t128 * t202;
t74 = qJD(2) * (t184 + t188) + t124 * t198;
t263 = (t74 - t219) * t124 + t127 * (t73 - t220);
t123 = cos(pkin(6));
t109 = t123 * qJD(1) + qJD(3);
t216 = t125 * t109;
t208 = qJD(1) * t121;
t181 = t129 * t208;
t104 = qJD(2) * pkin(2) + t181;
t189 = t126 * t208;
t68 = t120 * t104 + t122 * t189;
t65 = qJD(2) * pkin(8) + t68;
t49 = t128 * t65 + t216;
t44 = qJD(4) * pkin(9) + t49;
t151 = -t128 * pkin(4) - t125 * pkin(9) - pkin(3);
t105 = t120 * t189;
t67 = t122 * t104 - t105;
t52 = qJD(2) * t151 - t67;
t12 = t124 * t52 + t127 * t44;
t256 = t128 * t109 - t125 * t65;
t85 = (t120 * t126 - t122 * t129) * t121;
t79 = qJD(2) * t85;
t72 = qJD(1) * t79;
t18 = qJD(4) * t256 - t128 * t72;
t204 = qJD(5) * t124;
t53 = (t103 + t77) * qJD(2);
t176 = t124 * t18 - t127 * t53 + t44 * t203 + t52 * t204;
t149 = -t12 * t110 - t176;
t186 = t110 * t204;
t217 = t124 * t128;
t260 = (t125 * (t97 + t201) - t110 * t217) * qJD(2) + t186;
t10 = -t110 * qJ(6) + t12;
t115 = t125 * t199;
t171 = pkin(5) * t115;
t2 = -t171 + t176;
t259 = t10 * t110 + t2;
t113 = t120 * pkin(2) + pkin(8);
t214 = t127 * t128;
t249 = t122 * pkin(2);
t95 = t151 - t249;
t258 = t113 * t214 + t124 * t95;
t80 = t122 * t181 - t105;
t257 = t124 * t264 + t95 * t203 - t80 * t214;
t255 = t74 + t219;
t253 = qJD(5) * t258 - t264 * t127 - t80 * t217;
t252 = t99 ^ 2;
t251 = pkin(9) * t99;
t250 = pkin(9) * t110;
t43 = -qJD(4) * pkin(4) - t256;
t17 = t97 * pkin(5) - t99 * qJ(6) + t43;
t248 = t17 * t99;
t157 = t123 * t128 - t86 * t125;
t19 = t49 * qJD(4) - t125 * t72;
t247 = t19 * t157;
t78 = qJD(2) * t86;
t71 = qJD(1) * t78;
t246 = t71 * t85;
t8 = t74 * pkin(5) + t73 * qJ(6) - t99 * qJD(6) + t19;
t245 = t8 * t124;
t244 = t8 * t127;
t243 = t80 * t97;
t242 = t80 * t99;
t241 = t99 * t97;
t206 = qJD(4) * t125;
t240 = (-t113 * t204 - qJD(6)) * t128 + (-t113 * t127 + qJ(6)) * t206 + t257;
t177 = t113 * t124 + pkin(5);
t239 = -t177 * t206 + t253;
t238 = (-t125 * t201 - t128 * t204) * t113 + t257;
t237 = t125 * t113 * t202 - t253;
t162 = pkin(5) * t124 - qJ(6) * t127;
t236 = t216 + (t162 * qJD(2) + t65) * t128 - t162 * qJD(5) + t124 * qJD(6);
t102 = t167 * qJD(2);
t28 = t124 * t102 + t127 * t256;
t221 = t74 * t127;
t187 = t128 * t201;
t82 = t97 * t187;
t235 = -t125 * t221 - t82;
t231 = t124 * t43;
t230 = t125 * t80;
t229 = t127 * t43;
t228 = t127 * t95;
t227 = t128 * t73;
t226 = t128 * t74;
t225 = t19 * t124;
t224 = t19 * t125;
t223 = t19 * t127;
t222 = t73 * t124;
t131 = qJD(2) ^ 2;
t218 = t121 * t131;
t130 = qJD(4) ^ 2;
t212 = t130 * t125;
t211 = t130 * t128;
t11 = -t124 * t44 + t127 * t52;
t210 = qJD(6) - t11;
t118 = t125 ^ 2;
t119 = t128 ^ 2;
t209 = t118 - t119;
t205 = qJD(4) * t128;
t84 = t99 * t184;
t197 = -t125 * t222 + t99 * t188 + t84;
t196 = t124 * t250;
t195 = t127 * t250;
t194 = pkin(9) * t206;
t193 = pkin(9) * t201;
t192 = t97 ^ 2 - t252;
t190 = t125 * t131 * t128;
t185 = t125 * t204;
t183 = t110 * t203;
t182 = t110 * t207;
t64 = -qJD(2) * pkin(3) - t67;
t175 = -qJD(2) * t64 + t72;
t174 = qJD(5) * t97 - t73;
t107 = t118 * t179;
t173 = t107 - t227;
t170 = t125 * t183;
t169 = qJ(6) * t115;
t168 = t128 * t115;
t166 = t174 * pkin(9);
t9 = t110 * pkin(5) + t210;
t165 = t10 * t127 + t124 * t9;
t164 = -t10 * t124 + t127 * t9;
t163 = t127 * pkin(5) + t124 * qJ(6);
t161 = -t11 * t127 - t12 * t124;
t160 = t11 * t124 - t12 * t127;
t63 = t123 * t125 + t86 * t128;
t35 = t85 * t124 + t63 * t127;
t34 = t63 * t124 - t85 * t127;
t159 = t125 * t256 - t128 * t49;
t27 = t127 * t102 - t124 * t256;
t152 = qJD(2) * t118 - t110 * t128;
t150 = t113 + t162;
t147 = t77 * qJD(2) - t113 * t130 - t71;
t33 = qJD(4) * t157 - t79 * t128;
t6 = t35 * qJD(5) + t33 * t124 - t78 * t127;
t7 = -t34 * qJD(5) + t78 * t124 + t33 * t127;
t146 = -t34 * t73 - t35 * t74 + t6 * t99 - t7 * t97;
t114 = -pkin(3) - t249;
t145 = qJD(4) * (qJD(2) * t114 + t64 + t80);
t144 = -t124 * t53 - t127 * t18 - t52 * t203 + t44 * t204;
t143 = t97 * t185 + t235;
t142 = t152 * t124;
t141 = t170 + t226;
t32 = qJD(4) * t63 - t79 * t125;
t139 = t6 * t110 - t34 * t115 - t157 * t74 + t32 * t97;
t138 = -t124 * t220 - t221;
t1 = -t110 * qJD(6) - t144 + t169;
t137 = t164 * qJD(5) + t1 * t127 + t2 * t124;
t136 = -t7 * t110 + t35 * t115 - t157 * t73 - t32 * t99;
t135 = t161 * qJD(5) + t124 * t176 - t127 * t144;
t134 = t224 + t18 * t128 + (-t125 * t49 - t128 * t256) * qJD(4);
t133 = t97 * t184 + (t125 * t74 + t97 * t205) * t124;
t90 = t97 * t206;
t132 = -qJD(4) * t142 + t170 - t226 + t90;
t106 = -pkin(4) - t163;
t96 = t110 * t187;
t89 = t99 * t206;
t81 = (-t110 - t200) * t206;
t75 = t150 * t125;
t66 = pkin(9) * t221;
t61 = t99 * pkin(5) + t97 * qJ(6);
t57 = -t113 * t217 + t228;
t55 = t177 * t128 - t228;
t54 = -t128 * qJ(6) + t258;
t47 = -t73 - t220;
t40 = (t163 * qJD(5) - qJD(6) * t127) * t125 + t150 * t205;
t39 = -t183 + (t110 * t214 + (-t99 + t202) * t125) * qJD(2);
t26 = -t127 * t219 - t222;
t24 = -pkin(5) * t207 - t27;
t23 = qJ(6) * t207 + t28;
t21 = t99 * t187 + (-t73 * t127 - t99 * t204) * t125;
t13 = t110 * t185 + t107 + t227 + t89 - t96;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 * t218, -t129 * t218, 0, 0, 0, 0, 0, 0, 0, 0, -t78 * qJD(2), t79 * qJD(2), 0, -t67 * t78 - t68 * t79 - t72 * t86 + t246, 0, 0, 0, 0, 0, 0, -t32 * qJD(4) + (-t128 * t78 + t85 * t206) * qJD(2), -t33 * qJD(4) + (t125 * t78 + t205 * t85) * qJD(2) (t125 * t32 + t128 * t33 + (-t125 * t63 - t128 * t157) * qJD(4)) * qJD(2), t18 * t63 - t256 * t32 + t49 * t33 + t64 * t78 + t246 - t247, 0, 0, 0, 0, 0, 0, t139, -t136, t146, -t11 * t6 + t12 * t7 - t144 * t35 + t176 * t34 + t43 * t32 - t247, 0, 0, 0, 0, 0, 0, t139, t146, t136, t1 * t35 + t10 * t7 - t157 * t8 + t17 * t32 + t2 * t34 + t9 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t156 * t208 + t77) * qJD(2) (qJD(1) * t85 + t80) * qJD(2), 0, t67 * t77 - t68 * t80 + (-t120 * t72 - t122 * t71) * pkin(2), 0.2e1 * t168, -0.2e1 * t209 * t199, t211, -0.2e1 * t168, -t212, 0, t125 * t145 + t128 * t147, -t125 * t147 + t128 * t145 (-t118 - t119) * t80 * qJD(2) + t134, t113 * t134 + t71 * t114 + t159 * t80 - t64 * t77, t21, t143 - t197, t13, t133 (-t125 * t97 - t142) * qJD(4) + t141, t81, -t237 * t110 + (t176 + (t113 * t97 + t231) * qJD(4)) * t128 + (t43 * t203 + t113 * t74 + t225 - t243 + (qJD(2) * t57 + t11) * qJD(4)) * t125, t238 * t110 + (-t144 + (t113 * t99 + t229) * qJD(4)) * t128 + (-t43 * t204 - t113 * t73 + t223 - t242 + (-qJD(2) * t258 - t12) * qJD(4)) * t125, t57 * t73 - t258 * t74 - t237 * t99 - t238 * t97 + t161 * t205 + (qJD(5) * t160 + t124 * t144 + t127 * t176) * t125, -t43 * t230 - t144 * t258 - t176 * t57 + t238 * t12 + (t205 * t43 + t224) * t113 + t237 * t11, t21, t13, t82 + (-t204 * t97 + t221) * t125 + t197, t81, t152 * t202 - t141 + t90, t133, t40 * t97 + t75 * t74 + (t17 * t202 + t2) * t128 + t239 * t110 + (t17 * t203 + t245 - t243 + (-qJD(2) * t55 - t9) * qJD(4)) * t125, -t54 * t74 - t55 * t73 + t239 * t99 - t240 * t97 + t164 * t205 + (-qJD(5) * t165 - t1 * t124 + t127 * t2) * t125, -t40 * t99 + t75 * t73 + (-t17 * t201 - t1) * t128 - t240 * t110 + (t17 * t204 - t244 + t242 + (qJD(2) * t54 + t10) * qJD(4)) * t125, t1 * t54 + t2 * t55 + t8 * t75 + t239 * t9 + (t40 - t230) * t17 + t240 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t212, -t211, 0, -qJD(4) * t159 + t18 * t125 - t19 * t128, 0, 0, 0, 0, 0, 0, t132, t89 + (-t185 + t187) * t110 - t173, t84 + (t125 * t174 + t205 * t99) * t124 + t235 (-qJD(4) * t160 - t19) * t128 + (qJD(4) * t43 + t135) * t125, 0, 0, 0, 0, 0, 0, t132, t143 + t197, -t96 + (-qJD(4) * t99 + t186) * t125 + t173 (qJD(4) * t165 - t8) * t128 + (qJD(4) * t17 + t137) * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, t209 * t131, 0, t190, 0, 0, t175 * t125, t175 * t128, 0, 0, t26, -t263, t39, t138, t260, t182, -pkin(4) * t74 + t27 * t110 - t223 - t49 * t97 + (t195 + t231) * qJD(5) + (-t11 * t125 + (-t128 * t43 - t194) * t124) * qJD(2), pkin(4) * t73 - t28 * t110 + t225 - t49 * t99 + (-t196 + t229) * qJD(5) + (-t43 * t214 + (t12 - t193) * t125) * qJD(2), t27 * t99 + t28 * t97 - t66 + (t11 * t200 - t144 + (-t11 + t251) * qJD(5)) * t127 + (t166 - t149) * t124, -t19 * pkin(4) + pkin(9) * t135 - t11 * t27 - t12 * t28 - t43 * t49, t26, t39, t263, t182, -t260, t138, t106 * t74 - t24 * t110 - t244 - t236 * t97 + (t124 * t17 + t195) * qJD(5) + (t125 * t9 + (-t128 * t17 - t194) * t124) * qJD(2), t23 * t97 - t24 * t99 - t66 + (-t9 * t200 + t1 + (t9 + t251) * qJD(5)) * t127 + (t166 + t259) * t124, t106 * t73 + t23 * t110 - t245 + t236 * t99 + (-t127 * t17 + t196) * qJD(5) + (t17 * t214 + (-t10 + t193) * t125) * qJD(2), t137 * pkin(9) - t10 * t23 + t8 * t106 - t236 * t17 - t9 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, -t192, t47, -t241, -t255, t115, -t43 * t99 + t149, -t11 * t110 + t43 * t97 + t144, 0, 0, t241, t47, t192, t115, t255, -t241, -t61 * t97 + t149 + 0.2e1 * t171 - t248, pkin(5) * t73 - t74 * qJ(6) + (t10 - t12) * t99 + (t9 - t210) * t97, 0.2e1 * t169 - t17 * t97 + t61 * t99 + (-0.2e1 * qJD(6) + t11) * t110 - t144, -t2 * pkin(5) + t1 * qJ(6) + t10 * t210 - t9 * t12 - t17 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115 + t241, t47, -t110 ^ 2 - t252, t248 + t259;];
tauc_reg  = t3;
