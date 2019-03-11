% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:32:54
% EndTime: 2019-03-08 19:33:05
% DurationCPUTime: 4.26s
% Computational Cost: add. (4728->364), mult. (12142->547), div. (0->0), fcn. (9528->12), ass. (0->202)
t159 = sin(qJ(4));
t162 = cos(qJ(4));
t189 = pkin(4) * t159 - qJ(5) * t162;
t112 = t189 * qJD(4) - t159 * qJD(5);
t154 = sin(pkin(6));
t153 = sin(pkin(11));
t156 = cos(pkin(11));
t160 = sin(qJ(2));
t163 = cos(qJ(2));
t183 = t153 * t163 + t156 * t160;
t104 = t183 * t154;
t93 = qJD(1) * t104;
t262 = t112 - t93;
t152 = sin(pkin(12));
t145 = pkin(2) * t153 + pkin(8);
t214 = qJD(4) * t159;
t206 = t145 * t214;
t120 = t152 * t206;
t155 = cos(pkin(12));
t230 = t152 * t162;
t220 = qJD(1) * t154;
t210 = t160 * t220;
t131 = t153 * t210;
t209 = t163 * t220;
t96 = t156 * t209 - t131;
t245 = -t155 * t262 - t96 * t230 - t120;
t227 = t155 * t162;
t261 = t152 * t262 - t96 * t227;
t215 = qJD(4) * t155;
t218 = qJD(2) * t159;
t121 = -t152 * t218 + t215;
t207 = t155 * t218;
t216 = qJD(4) * t152;
t122 = t207 + t216;
t158 = sin(qJ(6));
t161 = cos(qJ(6));
t74 = -t161 * t121 + t122 * t158;
t260 = t74 ^ 2;
t77 = t121 * t158 + t122 * t161;
t259 = t77 ^ 2;
t179 = pkin(5) * t159 - pkin(9) * t227;
t172 = t179 * qJD(4);
t258 = -t172 + t245;
t228 = t155 * t159;
t257 = -(-pkin(9) * t230 - t145 * t228) * qJD(4) - t261;
t217 = qJD(2) * t162;
t142 = -qJD(6) + t217;
t256 = t74 * t142;
t103 = (t153 * t160 - t156 * t163) * t154;
t212 = qJD(6) * t161;
t231 = t152 * t158;
t255 = -qJD(6) * t231 + t155 * t212;
t211 = qJD(2) * qJD(4);
t200 = t162 * t211;
t192 = t152 * t200;
t213 = qJD(4) * t162;
t204 = t155 * t213;
t193 = t161 * t204;
t40 = (qJD(6) * t122 + t192) * t158 - qJD(2) * t193 - t121 * t212;
t157 = cos(pkin(6));
t141 = qJD(1) * t157 + qJD(3);
t95 = qJD(2) * t103;
t91 = qJD(1) * t95;
t234 = t141 * t213 - t162 * t91;
t130 = qJD(2) * pkin(2) + t209;
t89 = t153 * t130 + t156 * t210;
t84 = qJD(2) * pkin(8) + t89;
t78 = t159 * t84;
t31 = (qJD(5) - t78) * qJD(4) + t234;
t59 = (t93 + t112) * qJD(2);
t14 = -t152 * t31 + t155 * t59;
t12 = qJD(2) * t172 + t14;
t15 = t152 * t59 + t155 * t31;
t13 = -pkin(9) * t192 + t15;
t226 = t159 * t141;
t63 = t162 * t84 + t226;
t56 = qJD(4) * qJ(5) + t63;
t178 = -pkin(4) * t162 - qJ(5) * t159 - pkin(3);
t88 = t156 * t130 - t131;
t68 = t178 * qJD(2) - t88;
t22 = -t152 * t56 + t155 * t68;
t18 = -pkin(5) * t217 - pkin(9) * t122 + t22;
t23 = t152 * t68 + t155 * t56;
t21 = pkin(9) * t121 + t23;
t185 = t158 * t21 - t161 * t18;
t1 = -t185 * qJD(6) + t158 * t12 + t161 * t13;
t251 = pkin(2) * t156;
t119 = t178 - t251;
t108 = t155 * t119;
t65 = -pkin(9) * t228 + t108 + (-t145 * t152 - pkin(5)) * t162;
t80 = t152 * t119 + t145 * t227;
t70 = -pkin(9) * t152 * t159 + t80;
t25 = t158 * t65 + t161 * t70;
t254 = t25 * qJD(6) - t158 * t257 + t161 * t258;
t24 = -t158 * t70 + t161 * t65;
t253 = -t24 * qJD(6) + t158 * t258 + t161 * t257;
t125 = t152 * t161 + t155 * t158;
t127 = t189 * qJD(2);
t62 = t141 * t162 - t78;
t37 = t155 * t127 - t152 * t62;
t32 = t179 * qJD(2) + t37;
t208 = t152 * t217;
t38 = t152 * t127 + t155 * t62;
t36 = -pkin(9) * t208 + t38;
t248 = pkin(9) + qJ(5);
t134 = t248 * t152;
t135 = t248 * t155;
t87 = -t134 * t158 + t135 * t161;
t252 = -t125 * qJD(5) - t87 * qJD(6) + t158 * t36 - t161 * t32;
t35 = t63 * qJD(4) - t159 * t91;
t81 = t104 * t159 - t157 * t162;
t250 = t35 * t81;
t249 = t77 * t74;
t124 = -t161 * t155 + t231;
t86 = -t134 * t161 - t135 * t158;
t247 = t124 * qJD(5) - t86 * qJD(6) + t158 * t32 + t161 * t36;
t106 = t124 * t159;
t173 = t125 * t162;
t169 = qJD(4) * t173;
t41 = qJD(2) * t169 + t77 * qJD(6);
t115 = t125 * qJD(6);
t205 = t152 * t213;
t66 = t159 * t115 + t158 * t205 - t193;
t246 = t106 * t41 + t66 * t74;
t194 = t155 * t206;
t244 = -t194 + t261;
t105 = t125 * t159;
t199 = t159 * t211;
t67 = t255 * t159 + t169;
t243 = -t105 * t199 + t67 * t142;
t94 = qJD(2) * t104;
t90 = qJD(1) * t94;
t242 = t103 * t90;
t241 = t159 * t96;
t240 = t35 * t152;
t239 = t35 * t155;
t238 = t35 * t159;
t237 = t35 * t162;
t236 = -t124 * t217 - t255;
t235 = -qJD(2) * t173 + t115;
t233 = t121 * t159;
t151 = t162 ^ 2;
t165 = qJD(2) ^ 2;
t232 = t151 * t165;
t229 = t154 * t165;
t164 = qJD(4) ^ 2;
t225 = t164 * t159;
t224 = t164 * t162;
t150 = t159 ^ 2;
t222 = t150 - 0.2e1 * t151;
t221 = t150 - t151;
t219 = qJD(2) * t152;
t201 = t155 * t211;
t198 = pkin(5) * t152 + t145;
t197 = t162 * t40 + t77 * t214;
t196 = -qJD(4) * pkin(4) + qJD(5);
t195 = -t121 + t215;
t191 = t162 * t199;
t190 = -t105 * t40 + t67 * t77;
t187 = -t14 * t152 + t15 * t155;
t186 = -t152 * t22 + t155 * t23;
t6 = t158 * t18 + t161 * t21;
t82 = t104 * t162 + t157 * t159;
t44 = t103 * t155 - t152 * t82;
t45 = t103 * t152 + t155 * t82;
t16 = -t158 * t45 + t161 * t44;
t17 = t158 * t44 + t161 * t45;
t184 = t159 * t62 - t162 * t63;
t181 = qJD(2) * t195;
t180 = qJD(2) * (-t122 + t216);
t177 = t162 * t41 - t74 * t214;
t176 = t93 * qJD(2) - t145 * t164 - t90;
t146 = -pkin(3) - t251;
t83 = -qJD(2) * pkin(3) - t88;
t175 = qJD(4) * (qJD(2) * t146 + t83 + t96);
t174 = t162 * t180;
t54 = t196 - t62;
t171 = t106 * t199 - t142 * t66;
t168 = -qJ(5) * t214 + (t196 - t54) * t162;
t2 = -t6 * qJD(6) + t161 * t12 - t158 * t13;
t28 = pkin(5) * t192 + t35;
t34 = -t84 * t214 + t234;
t166 = t238 + t34 * t162 + (-t159 * t63 - t162 * t62) * qJD(4);
t149 = t155 ^ 2;
t148 = t152 ^ 2;
t147 = -pkin(5) * t155 - pkin(4);
t140 = t159 * t165 * t162;
t137 = -0.2e1 * t191;
t111 = t198 * t159;
t110 = t122 * t214;
t102 = t198 * t213;
t99 = t121 * t204;
t79 = -t145 * t230 + t108;
t48 = t226 + (pkin(5) * t219 + t84) * t162;
t43 = -t81 * qJD(4) - t95 * t162;
t42 = t82 * qJD(4) - t95 * t159;
t39 = -pkin(5) * t121 + t54;
t27 = t152 * t94 + t155 * t43;
t26 = -t152 * t43 + t155 * t94;
t4 = -t17 * qJD(6) - t158 * t27 + t161 * t26;
t3 = t16 * qJD(6) + t158 * t26 + t161 * t27;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160 * t229, -t163 * t229, 0, 0, 0, 0, 0, 0, 0, 0, -t94 * qJD(2), t95 * qJD(2), 0, -t104 * t91 - t88 * t94 - t89 * t95 + t242, 0, 0, 0, 0, 0, 0, -t42 * qJD(4) + (t103 * t214 - t162 * t94) * qJD(2), -t43 * qJD(4) + (t103 * t213 + t159 * t94) * qJD(2) (t159 * t42 + t162 * t43 + (-t159 * t82 + t162 * t81) * qJD(4)) * qJD(2), t34 * t82 - t42 * t62 + t43 * t63 + t83 * t94 + t242 + t250, 0, 0, 0, 0, 0, 0, -t42 * t121 + (-t162 * t26 + (t159 * t44 + t230 * t81) * qJD(4)) * qJD(2), t42 * t122 + (t162 * t27 + (-t159 * t45 + t227 * t81) * qJD(4)) * qJD(2), t27 * t121 - t26 * t122 + (-t152 * t45 - t155 * t44) * t200, t14 * t44 + t15 * t45 + t22 * t26 + t23 * t27 + t42 * t54 + t250, 0, 0, 0, 0, 0, 0, -t142 * t4 + t16 * t199 + t41 * t81 + t42 * t74, t142 * t3 - t17 * t199 - t40 * t81 + t42 * t77, t16 * t40 - t17 * t41 - t3 * t74 - t4 * t77, t1 * t17 + t16 * t2 - t185 * t4 + t28 * t81 + t3 * t6 + t39 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t183 * t220 + t93) * qJD(2) (qJD(1) * t103 + t96) * qJD(2), 0, t88 * t93 - t89 * t96 + (-t153 * t91 - t156 * t90) * pkin(2), 0.2e1 * t191, -0.2e1 * t221 * t211, t224, t137, -t225, 0, t159 * t175 + t176 * t162, -t176 * t159 + t162 * t175 (-t150 - t151) * t96 * qJD(2) + t166, t166 * t145 + t90 * t146 + t184 * t96 - t83 * t93 (t122 * t155 + t149 * t218) * t213, t99 + (-t122 - 0.2e1 * t207) * t205, t201 * t222 + t110 (-t121 * t152 + t148 * t218) * t213 (-t219 * t222 + t233) * qJD(4), t137 (t121 * t96 + t240 + (qJD(2) * t79 + t22) * qJD(4)) * t159 + (-t14 + (-t121 * t145 + t152 * t54) * qJD(4) + (t120 + t245) * qJD(2)) * t162 (-t122 * t96 + t239 + (-qJD(2) * t80 - t23) * qJD(4)) * t159 + (t15 + (t122 * t145 + t155 * t54) * qJD(4) + (t194 + t244) * qJD(2)) * t162 (-t14 * t155 - t15 * t152) * t159 + t245 * t122 + t244 * t121 + (-t152 * t23 - t155 * t22 + (-t152 * t80 - t155 * t79) * qJD(2)) * t213, -t54 * t241 + t14 * t79 + t15 * t80 + t244 * t23 - t245 * t22 + (t213 * t54 + t238) * t145, t106 * t40 - t66 * t77, -t190 + t246, -t171 + t197, t105 * t41 + t67 * t74, t177 + t243 (-t142 - t217) * t214, t102 * t74 + t105 * t28 + t111 * t41 - t162 * t2 + t39 * t67 + t254 * t142 + (-t74 * t96 + (qJD(2) * t24 - t185) * qJD(4)) * t159, t1 * t162 + t102 * t77 - t28 * t106 - t111 * t40 - t39 * t66 - t253 * t142 + (-t77 * t96 + (-qJD(2) * t25 - t6) * qJD(4)) * t159, -t1 * t105 + t2 * t106 - t185 * t66 + t24 * t40 - t25 * t41 + t253 * t74 + t254 * t77 - t6 * t67, t1 * t25 + t28 * t111 + t2 * t24 - t253 * t6 + t254 * t185 + (t102 - t241) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225, -t224, 0, -t184 * qJD(4) + t34 * t159 - t237, 0, 0, 0, 0, 0, 0 (-t150 * t219 - t233) * qJD(4), -t150 * t201 + t110, t122 * t205 + t99, -t237 + t187 * t159 + (t159 * t54 + t162 * t186) * qJD(4), 0, 0, 0, 0, 0, 0, -t177 + t243, t171 + t197, t190 + t246, -t1 * t106 - t105 * t2 - t162 * t28 + t185 * t67 + t214 * t39 - t6 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, t221 * t165, 0, t140, 0, 0 (-qJD(2) * t83 + t91) * t159, -t83 * t217 + (t62 + t78) * qJD(4) - t234, 0, 0, t155 * t174 (-t121 * t155 + t122 * t152 + (-t148 + t149) * qJD(4)) * t217, t155 * t232 + t159 * t180, -t195 * t208, -t152 * t232 + t159 * t181, t140, t63 * t121 - t239 + (t152 * t168 - t159 * t22 + t162 * t37) * qJD(2), -t63 * t122 + t240 + (t155 * t168 + t159 * t23 - t162 * t38) * qJD(2), -t38 * t121 + t37 * t122 + (qJD(5) * t121 + t217 * t22 + t15) * t155 + (qJD(5) * t122 + t217 * t23 - t14) * t152, -t35 * pkin(4) + qJ(5) * t187 + qJD(5) * t186 - t22 * t37 - t23 * t38 - t54 * t63, -t40 * t125 - t236 * t77, t40 * t124 - t125 * t41 - t235 * t77 + t236 * t74, t236 * t142 + (qJD(4) * t125 - t77) * t218, t41 * t124 + t235 * t74, t235 * t142 + (-qJD(4) * t124 + t74) * t218, t142 * t218, t28 * t124 + t147 * t41 - t48 * t74 + t235 * t39 - t252 * t142 + (qJD(4) * t86 + t185) * t218, t28 * t125 - t147 * t40 - t48 * t77 - t236 * t39 - t247 * t142 + (-qJD(4) * t87 + t6) * t218, -t1 * t124 - t125 * t2 - t185 * t236 - t235 * t6 + t247 * t74 - t252 * t77 + t40 * t86 - t41 * t87, t1 * t87 + t28 * t147 - t185 * t252 + t2 * t86 - t247 * t6 - t39 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t162 * t181, -t121 ^ 2 - t122 ^ 2, -t23 * t121 + t22 * t122 + t35, 0, 0, 0, 0, 0, 0, -t77 * t142 + t41, -t40 + t256, -t259 - t260, -t185 * t77 + t6 * t74 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, t259 - t260, -t40 - t256, -t249, -t125 * t200 + (-qJD(6) - t142) * t77, t199, -t6 * t142 - t39 * t77 + t2, t142 * t185 + t39 * t74 - t1, 0, 0;];
tauc_reg  = t5;
