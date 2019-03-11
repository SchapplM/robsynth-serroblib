% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% tauc_reg [6x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:19:57
% EndTime: 2019-03-09 04:20:04
% DurationCPUTime: 2.66s
% Computational Cost: add. (2139->327), mult. (4479->456), div. (0->0), fcn. (2780->6), ass. (0->185)
t132 = sin(qJ(6));
t133 = sin(qJ(5));
t135 = cos(qJ(6));
t136 = cos(qJ(5));
t215 = t135 * t136;
t159 = t132 * t133 - t215;
t134 = sin(qJ(3));
t206 = qJD(1) * t134;
t181 = t136 * t206;
t84 = qJD(3) * t133 - t181;
t203 = qJD(3) * t136;
t86 = t133 * t206 + t203;
t27 = t132 * t86 + t135 * t84;
t253 = (t159 * t206 - t27) * qJD(3);
t139 = -pkin(1) - pkin(7);
t238 = pkin(4) - t139;
t137 = cos(qJ(3));
t205 = qJD(1) * t137;
t241 = qJD(5) + qJD(6);
t89 = t132 * t136 + t133 * t135;
t251 = t241 * t89;
t236 = -t89 * t205 - t251;
t161 = t132 * t84 - t135 * t86;
t252 = t161 * t27;
t62 = t159 * t134;
t250 = t161 ^ 2 - t27 ^ 2;
t117 = qJD(5) + t205;
t108 = qJD(6) + t117;
t196 = qJD(6) * t135;
t197 = qJD(6) * t132;
t195 = qJD(1) * qJD(3);
t177 = t137 * t195;
t199 = qJD(5) * t133;
t44 = -qJD(3) * t199 + qJD(5) * t181 + t133 * t177;
t45 = t86 * qJD(5) - t136 * t177;
t8 = -t132 * t45 + t135 * t44 - t84 * t196 - t86 * t197;
t249 = t108 * t27 + t8;
t138 = -pkin(3) - pkin(8);
t201 = qJD(3) * t138;
t109 = t139 * qJD(1) + qJD(2);
t93 = t137 * t109;
t242 = qJD(4) - t93;
t210 = pkin(4) * t205 + t242;
t42 = t210 + t201;
t220 = qJ(4) * t137;
t163 = pkin(8) * t134 - t220;
t209 = pkin(3) * t206 + qJD(1) * qJ(2);
t50 = t163 * qJD(1) + t209;
t18 = t133 * t42 + t136 * t50;
t13 = -pkin(9) * t84 + t18;
t11 = t13 * t197;
t128 = qJD(3) * qJ(4);
t92 = t134 * t109;
t64 = -pkin(4) * t206 + t92;
t52 = t128 + t64;
t24 = pkin(5) * t84 + t52;
t248 = t24 * t27 + t11;
t160 = (qJD(3) * pkin(8) - qJD(4)) * t137;
t127 = qJD(1) * qJD(2);
t176 = t134 * t195;
t189 = pkin(3) * t177 + qJ(4) * t176 + t127;
t31 = qJD(1) * t160 + t189;
t204 = qJD(3) * t134;
t87 = t109 * t204;
t57 = -pkin(4) * t176 + t87;
t172 = -t133 * t31 + t136 * t57;
t145 = -t18 * qJD(5) + t172;
t4 = -pkin(5) * t176 - pkin(9) * t44 + t145;
t198 = qJD(5) * t136;
t194 = -t133 * t57 - t136 * t31 - t42 * t198;
t152 = -t50 * t199 - t194;
t5 = -pkin(9) * t45 + t152;
t187 = -t132 * t5 + t135 * t4;
t17 = -t133 * t50 + t136 * t42;
t12 = -pkin(9) * t86 + t17;
t10 = pkin(5) * t117 + t12;
t229 = t13 * t135;
t2 = t10 * t132 + t229;
t247 = -t2 * qJD(6) + t24 * t161 + t187;
t144 = t161 * qJD(6) - t132 * t44 - t135 * t45;
t245 = -t108 * t161 + t144;
t75 = -t92 - t128;
t224 = t137 * t75;
t126 = qJD(3) * qJD(4);
t202 = qJD(3) * t137;
t68 = -t109 * t202 - t126;
t232 = qJD(3) * pkin(3);
t167 = -qJD(4) + t232;
t71 = -t167 - t93;
t244 = ((-t71 + t93) * t134 + t224) * qJD(3) + t134 * t68;
t150 = t159 * t108;
t243 = t134 * pkin(3) + qJ(2);
t239 = 0.2e1 * t127;
t237 = pkin(9) - t138;
t186 = t136 * t205;
t217 = t133 * t137;
t191 = t132 * t217;
t235 = -qJD(1) * t191 - t132 * t199 - t133 * t197 + t135 * t186 + t241 * t215;
t91 = pkin(3) * t205 + qJ(4) * t206;
t61 = pkin(8) * t205 + t91;
t234 = t133 * t64 + t136 * t61;
t77 = t163 + t243;
t98 = t238 * t137;
t78 = t133 * t98;
t233 = t136 * t77 + t78;
t231 = t117 * t84;
t230 = t117 * t86;
t228 = t134 * t44;
t226 = t134 * t84;
t225 = t136 * t44;
t47 = -pkin(4) * t177 - t68;
t223 = t47 * t133;
t222 = t47 * t136;
t188 = -pkin(5) * t136 - pkin(4);
t221 = pkin(5) * t198 - t188 * t205 + t242;
t149 = qJD(3) * t89;
t219 = t117 * t138;
t216 = t134 * t136;
t214 = t136 * t137;
t140 = qJD(3) ^ 2;
t213 = t140 * t134;
t212 = t140 * t137;
t141 = qJD(1) ^ 2;
t211 = t141 * qJ(2);
t130 = t134 ^ 2;
t131 = t137 ^ 2;
t208 = t130 - t131;
t207 = t140 + t141;
t200 = qJD(4) * t137;
t193 = 0.2e1 * qJD(1);
t192 = t117 * t217;
t190 = t137 * t141 * t134;
t185 = t136 * t202;
t184 = t117 * t199;
t183 = t134 * t199;
t182 = t134 * t198;
t180 = pkin(3) * t202 + qJ(4) * t204 + qJD(2);
t179 = -pkin(9) * t134 - t77;
t97 = t237 * t136;
t178 = t235 * t108;
t175 = qJD(6) * t10 + t5;
t174 = t236 * t108;
t46 = t160 + t180;
t80 = t238 * t204;
t171 = -t133 * t46 - t136 * t80;
t170 = -t133 * t61 + t136 * t64;
t74 = -qJ(4) * t205 + t209;
t94 = -t220 + t243;
t169 = qJD(1) * t94 + t74;
t166 = t117 + t205;
t154 = -pkin(5) * t134 - pkin(9) * t217;
t95 = t237 * t133;
t165 = t154 * qJD(1) - qJD(6) * t95 - t237 * t199 + t170;
t164 = pkin(9) * t186 + t241 * t97 + t234;
t81 = t238 * t202;
t79 = t136 * t98;
t21 = pkin(5) * t137 + t179 * t133 + t79;
t22 = pkin(9) * t216 + t233;
t162 = t132 * t21 + t135 * t22;
t157 = t166 * t134;
t156 = -qJD(1) * t130 + t117 * t137;
t155 = t117 * (qJD(5) * t137 + qJD(1));
t43 = -qJD(1) * t200 + t189;
t59 = t180 - t200;
t153 = -qJD(1) * t59 + t139 * t140 - t43;
t151 = -t133 * t80 + t136 * t46 + t98 * t198 - t77 * t199;
t148 = t89 * qJD(1);
t147 = -t183 + t185;
t121 = t134 * t139;
t118 = pkin(5) * t133 + qJ(4);
t103 = t133 * t176;
t100 = t207 * t137;
t99 = t207 * t134;
t96 = -pkin(4) * t134 + t121;
t66 = t188 * t134 + t121;
t63 = t89 * t134;
t58 = t74 * t205;
t34 = -t147 * pkin(5) - t81;
t20 = pkin(5) * t45 + t47;
t15 = qJD(3) * t191 + t134 * t251 - t135 * t185;
t14 = t137 * t149 - t241 * t62;
t7 = t147 * pkin(9) + t151;
t6 = t154 * qJD(3) + (t179 * t136 - t78) * qJD(5) + t171;
t1 = t10 * t135 - t13 * t132;
t3 = [0, 0, 0, 0, t239, qJ(2) * t239, -0.2e1 * t137 * t176, 0.2e1 * t208 * t195, -t213, -t212, 0, -t139 * t213 + (qJ(2) * t202 + qJD(2) * t134) * t193, -t139 * t212 + (-qJ(2) * t204 + qJD(2) * t137) * t193, t244, t134 * t153 - t169 * t202, t137 * t153 + t169 * t204, -t139 * t244 + t43 * t94 + t59 * t74, t86 * t182 + (t202 * t86 + t228) * t133 (-t133 * t84 + t136 * t86) * t202 + (-t133 * t45 + t225 + (-t133 * t86 - t136 * t84) * qJD(5)) * t134, t117 * t182 + t137 * t44 + (t133 * t156 - t134 * t86) * qJD(3), -t117 * t183 - t137 * t45 + (t136 * t156 + t226) * qJD(3), -qJD(3) * t157, t171 * t117 - t81 * t84 + t96 * t45 + (-t203 * t52 + t172) * t137 + (-t233 * t117 - t18 * t137) * qJD(5) + (t52 * t199 - t222 + (-(-t133 * t77 + t79) * qJD(1) - t17) * qJD(3)) * t134, -t151 * t117 - t81 * t86 + t96 * t44 + ((qJD(3) * t52 + qJD(5) * t50) * t133 + t194) * t137 + (t52 * t198 + t223 + (t233 * qJD(1) + t18) * qJD(3)) * t134, -t14 * t161 + t63 * t8, -t14 * t27 + t144 * t63 + t15 * t161 - t62 * t8, t108 * t14 + t137 * t8 + (-qJD(1) * t63 + t161) * t204, -t108 * t15 + t137 * t144 + (qJD(1) * t62 + t27) * t204 (-t108 - t205) * t204 (-t132 * t7 + t135 * t6) * t108 + t187 * t137 + t34 * t27 - t66 * t144 + t20 * t62 + t24 * t15 + (-t108 * t162 - t137 * t2) * qJD(6) + (-(-t132 * t22 + t135 * t21) * qJD(1) - t1) * t204, t11 * t137 + t24 * t14 + t20 * t63 - t34 * t161 + t66 * t8 + (-(-qJD(6) * t22 + t6) * t108 - t4 * t137) * t132 + (-(qJD(6) * t21 + t7) * t108 - t175 * t137) * t135 + (qJD(1) * t162 + t2) * t204; 0, 0, 0, 0, -t141, -t211, 0, 0, 0, 0, 0, -t99, -t100, 0, t99, t100, -qJD(1) * t74 - t244, 0, 0, 0, 0, 0, t134 * t45 + t133 * t155 + (t137 * t84 + t166 * t216) * qJD(3), t228 + t136 * t155 + (-t133 * t157 + t137 * t86) * qJD(3), 0, 0, 0, 0, 0, t108 * t148 + (-qJD(3) * t150 - t144) * t134 + (t108 * t251 - t253) * t137, -qJD(1) * t150 + (-t108 * t149 + t8) * t134 + ((-t134 * t148 - t161) * qJD(3) - t241 * t150) * t137; 0, 0, 0, 0, 0, 0, t190, -t208 * t141, 0, 0, 0, -t137 * t211, t134 * t211 ((-t75 - t128) * t137 + (t167 + t71) * t134) * qJD(1), t206 * t91 + t58, 0.2e1 * t126 + (-t134 * t74 + t137 * t91) * qJD(1), -qJ(4) * t68 - qJD(4) * t75 - t74 * t91 + (t224 + (-t71 - t232) * t134) * t109, -t133 * t230 + t225 (-t45 - t230) * t136 + (-t44 + t231) * t133, -t184 + (-t192 + (t86 - t203) * t134) * qJD(1), -t117 * t198 + t103 + (-t117 * t214 - t226) * qJD(1), t117 * t206, qJ(4) * t45 + t223 - t170 * t117 + t210 * t84 + (-t133 * t219 + t136 * t52) * qJD(5) + (t52 * t214 + (-t136 * t201 + t17) * t134) * qJD(1), qJ(4) * t44 + t222 + t234 * t117 + t210 * t86 + (-t133 * t52 - t136 * t219) * qJD(5) + (-t18 * t134 + (t134 * t201 - t137 * t52) * t133) * qJD(1), -t159 * t8 - t161 * t236, -t144 * t159 + t161 * t235 - t236 * t27 - t8 * t89 (qJD(3) * t159 - t161) * t206 + t174, -t178 + (-t27 + t149) * t206, t108 * t206, -t118 * t144 + t20 * t89 + t221 * t27 + t235 * t24 + (t132 * t164 - t135 * t165) * t108 + (-(t132 * t95 - t135 * t97) * qJD(3) + t1) * t206, t118 * t8 - t20 * t159 - t221 * t161 + t236 * t24 + (t132 * t165 + t135 * t164) * t108 + ((-t132 * t97 - t135 * t95) * qJD(3) - t2) * t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, -t131 * t141 - t140, qJD(3) * t75 + t58 + t87, 0, 0, 0, 0, 0, -t184 - qJD(3) * t84 + (-t134 * t203 - t192) * qJD(1), -t117 ^ 2 * t136 - qJD(3) * t86 + t103, 0, 0, 0, 0, 0, t174 + t253, -t178 + (t206 * t89 + t161) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86 * t84, -t84 ^ 2 + t86 ^ 2, t44 + t231, t230 - t45, -t176, t117 * t18 - t52 * t86 + t145, t117 * t17 + t52 * t84 - t152, -t252, t250, t249, t245, -t176 -(-t12 * t132 - t229) * t108 + (-t108 * t197 - t135 * t176 - t27 * t86) * pkin(5) + t247 (-t108 * t13 - t4) * t132 + (t108 * t12 - t175) * t135 + (-t108 * t196 + t132 * t176 + t161 * t86) * pkin(5) + t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t252, t250, t249, t245, -t176, t108 * t2 + t247, t1 * t108 - t132 * t4 - t135 * t175 + t248;];
tauc_reg  = t3;
