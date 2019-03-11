% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:11:46
% EndTime: 2019-03-08 19:11:55
% DurationCPUTime: 3.99s
% Computational Cost: add. (4811->389), mult. (13448->609), div. (0->0), fcn. (12131->16), ass. (0->203)
t135 = sin(pkin(8));
t144 = sin(qJ(4));
t241 = t135 * t144;
t129 = pkin(10) * t241;
t148 = cos(qJ(4));
t139 = cos(pkin(8));
t233 = t139 * t148;
t234 = t139 * t144;
t141 = cos(pkin(6));
t128 = qJD(1) * t141 + qJD(2);
t136 = sin(pkin(7));
t149 = cos(qJ(3));
t237 = t136 * t149;
t115 = t128 * t237;
t140 = cos(pkin(7));
t138 = cos(pkin(14));
t137 = sin(pkin(6));
t226 = qJD(1) * t137;
t205 = t138 * t226;
t134 = sin(pkin(14));
t145 = sin(qJ(3));
t243 = t134 * t145;
t73 = t140 * t149 * t205 - t226 * t243 + t115;
t238 = t136 * t145;
t235 = t138 * t140;
t263 = (t134 * t149 + t145 * t235) * t137;
t74 = -qJD(1) * t263 - t128 * t238;
t269 = -t148 * t73 - t74 * t234 + (pkin(3) * t233 - t129) * qJD(4);
t143 = sin(qJ(5));
t147 = cos(qJ(5));
t224 = qJD(3) * t139;
t192 = qJD(4) + t224;
t225 = qJD(3) * t135;
t204 = t144 * t225;
t268 = -t143 * t204 + t147 * t192;
t92 = qJD(6) - t268;
t80 = t141 * t238 + t263;
t223 = qJD(3) * t148;
t203 = t135 * t223;
t267 = qJD(5) - t203;
t239 = t135 * t148;
t212 = pkin(10) * t239;
t248 = -t144 * t73 + t74 * t233 + (pkin(3) * t234 + t212) * qJD(4);
t103 = t212 + (pkin(3) * t144 + pkin(11)) * t139;
t183 = -pkin(4) * t148 - pkin(11) * t144;
t104 = (-pkin(3) + t183) * t135;
t182 = pkin(4) * t144 - pkin(11) * t148;
t160 = t182 * qJD(4);
t106 = t135 * t160;
t216 = qJD(5) * t147;
t218 = qJD(5) * t143;
t242 = t135 * t143;
t266 = -t103 * t218 + t104 * t216 + t143 * t106 + t269 * t147 + t74 * t242;
t240 = t135 * t147;
t256 = t147 * t103 + t143 * t104;
t265 = t256 * qJD(5) - t147 * t106 + t269 * t143 - t74 * t240;
t68 = qJD(3) * pkin(3) + t73;
t95 = t128 * t140 - t136 * t205;
t178 = t135 * t95 + t139 * t68;
t65 = pkin(10) * t225 - t74;
t33 = -t144 * t65 + t178 * t148;
t264 = t144 * t148;
t63 = t148 * t65;
t34 = t68 * t234 + t95 * t241 + t63;
t29 = pkin(11) * t192 + t34;
t88 = t139 * t95;
t42 = t88 + (qJD(3) * t183 - t68) * t135;
t14 = t143 * t42 + t147 * t29;
t165 = t149 * t235 - t243;
t157 = t165 * t137;
t66 = (qJD(1) * t157 + t115) * qJD(3);
t67 = qJD(3) * t74;
t16 = t33 * qJD(4) + t148 * t66 + t67 * t234;
t51 = (qJD(3) * t160 - t67) * t135;
t153 = -qJD(5) * t14 - t143 * t16 + t147 * t51;
t213 = qJD(3) * qJD(4);
t199 = t135 * t213;
t187 = t144 * t199;
t4 = -pkin(5) * t187 - t153;
t98 = t143 * t192 + t147 * t204;
t262 = (pkin(5) * t98 + t92 * pkin(12)) * t92 + t4;
t186 = t148 * t199;
t76 = qJD(5) * t98 + t143 * t186;
t109 = -t136 * t137 * t138 + t140 * t141;
t79 = t141 * t237 + t157;
t172 = t109 * t135 + t139 * t79;
t261 = -t80 * t144 + t148 * t172;
t142 = sin(qJ(6));
t146 = cos(qJ(6));
t170 = -t142 * t267 - t146 * t98;
t75 = t268 * qJD(5) + t147 * t186;
t41 = -qJD(6) * t170 + t142 * t75 - t146 * t187;
t60 = t67 * t233;
t17 = t144 * t66 - t60 + (t144 * t178 + t63) * qJD(4);
t10 = pkin(5) * t76 - pkin(12) * t75 + t17;
t12 = pkin(12) * t267 + t14;
t28 = -t192 * pkin(4) - t33;
t20 = -pkin(5) * t268 - t98 * pkin(12) + t28;
t179 = t12 * t142 - t146 * t20;
t161 = t143 * t51 + t147 * t16 + t42 * t216 - t29 * t218;
t3 = pkin(12) * t187 + t161;
t1 = -qJD(6) * t179 + t10 * t142 + t146 * t3;
t150 = qJD(3) ^ 2;
t70 = t142 * t98 - t146 * t267;
t260 = t70 * t92;
t259 = t170 * t92;
t105 = t182 * t225;
t258 = t143 * t105 + t147 * t33;
t222 = qJD(4) * t144;
t202 = t135 * t222;
t257 = -pkin(5) * t202 + t265;
t255 = t267 * t268;
t131 = t135 ^ 2;
t78 = t80 * qJD(3);
t254 = t131 * t78;
t253 = t142 * t76;
t252 = t146 * t76;
t214 = qJD(6) * t146;
t215 = qJD(6) * t142;
t40 = t142 * t187 + t146 * t75 + t214 * t267 - t98 * t215;
t251 = t40 * t142;
t249 = t98 * t267;
t247 = -t34 + t267 * (pkin(5) * t143 - pkin(12) * t147);
t246 = t267 * t143;
t245 = t267 * t147;
t244 = t131 * t150;
t236 = t136 * t150;
t232 = t144 * t145;
t231 = t144 * t149;
t230 = t145 * t148;
t229 = t147 * t148;
t228 = t148 * t149;
t227 = t144 ^ 2 - t148 ^ 2;
t221 = qJD(4) * t147;
t220 = qJD(4) * t148;
t219 = qJD(5) * t142;
t217 = qJD(5) * t146;
t211 = t92 * t219;
t210 = t92 * t217;
t208 = t142 * t239;
t206 = t145 * t236;
t201 = t135 * t220;
t200 = t135 * t139 * t150;
t195 = t146 * t92;
t126 = -pkin(5) * t147 - pkin(12) * t143 - pkin(4);
t194 = pkin(12) * t204 - qJD(6) * t126 + t258;
t193 = 0.2e1 * t131 * t213;
t191 = qJD(4) + 0.2e1 * t224;
t190 = t131 * t206;
t189 = t225 * t238;
t102 = t129 + (-pkin(3) * t148 - pkin(4)) * t139;
t111 = -t147 * t139 + t143 * t241;
t112 = t139 * t143 + t144 * t240;
t57 = pkin(5) * t111 - pkin(12) * t112 + t102;
t185 = -pkin(12) * t202 - qJD(6) * t57 - t266;
t59 = -pkin(12) * t239 + t256;
t83 = -qJD(5) * t111 + t147 * t201;
t84 = qJD(5) * t112 + t143 * t201;
t184 = -pkin(5) * t84 + pkin(12) * t83 + qJD(6) * t59 - t248;
t6 = t12 * t146 + t142 * t20;
t38 = t144 * t172 + t148 * t80;
t54 = t109 * t139 - t135 * t79;
t22 = t143 * t54 + t147 * t38;
t177 = -t142 * t261 + t146 * t22;
t176 = -t142 * t22 - t146 * t261;
t110 = -t135 * t237 + t139 * t140;
t162 = t139 * t231 + t230;
t82 = t136 * t162 + t140 * t241;
t56 = t110 * t143 + t147 * t82;
t163 = t139 * t228 - t232;
t81 = -t136 * t163 - t140 * t239;
t175 = t142 * t81 + t146 * t56;
t174 = -t142 * t56 + t146 * t81;
t13 = -t143 * t29 + t147 * t42;
t21 = t143 * t38 - t147 * t54;
t173 = t105 * t147 - t143 * t33;
t171 = t110 * t147 - t143 * t82;
t169 = -t103 * t143 + t104 * t147;
t167 = -t92 * t214 - t253;
t166 = -t92 * t215 + t252;
t85 = t112 * t142 + t146 * t239;
t11 = -pkin(5) * t267 - t13;
t154 = -pkin(12) * t76 + (t11 + t13) * t92;
t2 = -qJD(6) * t6 + t146 * t10 - t142 * t3;
t48 = -t135 * t68 + t88;
t152 = -qJD(4) * t178 - t48 * t225 - t66;
t91 = (t142 * t144 + t146 * t229) * t225;
t90 = t142 * t147 * t203 - t146 * t204;
t86 = t112 * t146 - t208;
t77 = t79 * qJD(3);
t58 = pkin(5) * t239 - t169;
t53 = t140 * t202 + (t162 * qJD(4) + (t139 * t230 + t231) * qJD(3)) * t136;
t52 = t140 * t201 + (t163 * qJD(4) + (-t139 * t232 + t228) * qJD(3)) * t136;
t47 = -qJD(6) * t208 + t112 * t214 + t142 * t83 - t146 * t202;
t46 = -qJD(6) * t85 + t142 * t202 + t146 * t83;
t31 = qJD(5) * t56 + t143 * t52 - t147 * t189;
t30 = qJD(5) * t171 + t143 * t189 + t147 * t52;
t23 = -pkin(5) * t204 - t173;
t19 = qJD(4) * t261 + t148 * t77 - t78 * t234;
t18 = qJD(4) * t38 + t144 * t77 + t78 * t233;
t8 = -qJD(5) * t21 + t147 * t19 + t78 * t242;
t7 = qJD(5) * t22 + t143 * t19 - t78 * t240;
t5 = [0, 0, 0, -t78 * qJD(3), -t77 * qJD(3), 0, 0, 0, 0, 0, -t18 * qJD(4) + (-t139 * t18 - t148 * t254 + t202 * t54) * qJD(3), -t19 * qJD(4) + (-t139 * t19 + t144 * t254 + t201 * t54) * qJD(3), 0, 0, 0, 0, 0, -t18 * t268 - t187 * t21 - t261 * t76 - t267 * t7, t18 * t98 - t187 * t22 - t261 * t75 - t267 * t8, 0, 0, 0, 0, 0 (-qJD(6) * t177 - t142 * t8 + t146 * t18) * t92 + t176 * t76 + t7 * t70 + t21 * t41 -(qJD(6) * t176 + t142 * t18 + t146 * t8) * t92 - t177 * t76 - t7 * t170 + t21 * t40; 0, 0, 0, -t206, -t149 * t236, 0, 0, 0, 0, 0, t110 * t187 - t148 * t190 - t192 * t53, t110 * t186 + t144 * t190 - t192 * t52, 0, 0, 0, 0, 0, t171 * t187 - t267 * t31 - t268 * t53 + t81 * t76, -t187 * t56 - t267 * t30 + t53 * t98 + t75 * t81, 0, 0, 0, 0, 0 (-qJD(6) * t175 - t142 * t30 + t146 * t53) * t92 + t174 * t76 + t31 * t70 - t171 * t41 -(qJD(6) * t174 + t142 * t53 + t146 * t30) * t92 - t175 * t76 - t31 * t170 - t171 * t40; 0, 0, 0, 0 (-t165 * t226 - t115 + t73) * qJD(3), t193 * t264, -t227 * t193, t191 * t201, -t191 * t202, 0, t48 * t202 - t17 * t139 + (t67 * t148 + (-pkin(3) * t222 - t148 * t74) * qJD(3)) * t131 - t248 * t192, t48 * t201 - t16 * t139 + (-t144 * t67 + (-pkin(3) * t220 + t144 * t74) * qJD(3)) * t131 - t269 * t192, t112 * t75 + t83 * t98, -t111 * t75 - t112 * t76 + t268 * t83 - t84 * t98, t267 * t83 + (-t148 * t75 + (qJD(3) * t112 + t98) * t222) * t135, -t267 * t84 + (t148 * t76 + (-qJD(3) * t111 + t268) * t222) * t135 (-t131 * t223 + t135 * t267) * t222, t102 * t76 + t17 * t111 + t28 * t84 - t248 * t268 - t265 * t267 + (-t153 * t148 + (qJD(3) * t169 + t13) * t222) * t135, t102 * t75 + t17 * t112 + t28 * t83 + t248 * t98 - t266 * t267 + (t161 * t148 + (-qJD(3) * t256 - t14) * t222) * t135, -t170 * t46 + t40 * t86, t170 * t47 - t40 * t85 - t41 * t86 - t46 * t70, t111 * t40 - t170 * t84 + t46 * t92 + t76 * t86, -t111 * t41 - t47 * t92 - t70 * t84 - t76 * t85, t111 * t76 + t84 * t92 (-t142 * t59 + t146 * t57) * t76 + t2 * t111 - t179 * t84 + t58 * t41 + t4 * t85 + t11 * t47 + (t142 * t185 - t146 * t184) * t92 + t257 * t70 -(t142 * t57 + t146 * t59) * t76 - t1 * t111 - t6 * t84 + t58 * t40 + t4 * t86 + t11 * t46 + (t142 * t184 + t146 * t185) * t92 - t257 * t170; 0, 0, 0, 0, 0, -t244 * t264, t227 * t244, -t148 * t200, t144 * t200, 0, t144 * t152 + t192 * t34 - t220 * t65 + t60, t33 * t192 + (qJD(4) * t65 - t139 * t67) * t144 + t152 * t148, t143 * t75 + t245 * t98 (t75 + t255) * t147 + (-t249 - t76) * t143, t267 * t216 + (-t267 * t229 + (t143 * qJD(4) - t98) * t144) * t225, -t267 * t218 + (t148 * t246 + (-t268 + t221) * t144) * t225, -t267 * t204, -pkin(4) * t76 - t17 * t147 - t173 * t267 + t34 * t268 + (-pkin(11) * t245 + t143 * t28) * qJD(5) + (-t13 * t144 + (-pkin(11) * t222 - t148 * t28) * t143) * t225, -pkin(4) * t75 + t17 * t143 + t258 * t267 - t34 * t98 + (pkin(11) * t246 + t147 * t28) * qJD(5) + (-t28 * t229 + (-pkin(11) * t221 + t14) * t144) * t225, t143 * t146 * t40 - (-t143 * t215 + t146 * t216 - t91) * t170, t70 * t91 - t170 * t90 + (t142 * t170 - t146 * t70) * t216 + (-t251 - t146 * t41 + (t142 * t70 + t146 * t170) * qJD(6)) * t143, -t91 * t92 + (-t40 + t210) * t147 + (-t170 * t267 + t166) * t143, t90 * t92 + (t41 - t211) * t147 + (-t267 * t70 + t167) * t143, -t147 * t76 + t246 * t92, t126 * t252 - t11 * t90 - t23 * t70 + (t142 * t194 + t146 * t247) * t92 + (t11 * t219 - t2 + (qJD(5) * t70 + t167) * pkin(11)) * t147 + (t11 * t214 + t4 * t142 - t267 * t179 + (t41 + t211) * pkin(11)) * t143, -t126 * t253 - t11 * t91 + t23 * t170 + (-t142 * t247 + t146 * t194) * t92 + (t11 * t217 + t1 + (-qJD(5) * t170 - t166) * pkin(11)) * t147 + (-t11 * t215 + t4 * t146 - t267 * t6 + (t40 + t210) * pkin(11)) * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98 * t268, -t268 ^ 2 + t98 ^ 2, t75 - t255, -t76 + t249, t187, t14 * t267 - t28 * t98 + t153, t13 * t267 - t268 * t28 - t161, -t170 * t195 + t251 (t40 - t260) * t146 + (-t41 + t259) * t142, t170 * t98 + t195 * t92 + t253, -t142 * t92 ^ 2 + t70 * t98 + t252, -t92 * t98, -pkin(5) * t41 - t14 * t70 + t154 * t142 - t146 * t262 + t179 * t98, -pkin(5) * t40 + t14 * t170 + t142 * t262 + t154 * t146 + t6 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170 * t70, t170 ^ 2 - t70 ^ 2, t40 + t260, -t41 - t259, t76, t11 * t170 + t6 * t92 + t2, t11 * t70 - t179 * t92 - t1;];
tauc_reg  = t5;
