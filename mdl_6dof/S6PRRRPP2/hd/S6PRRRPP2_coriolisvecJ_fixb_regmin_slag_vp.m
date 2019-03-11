% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRPP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:52:40
% EndTime: 2019-03-08 22:52:49
% DurationCPUTime: 3.17s
% Computational Cost: add. (3164->398), mult. (7787->507), div. (0->0), fcn. (5340->8), ass. (0->197)
t137 = cos(qJ(3));
t199 = t137 * qJD(2);
t119 = -qJD(4) + t199;
t136 = cos(qJ(4));
t133 = sin(qJ(4));
t203 = t133 * qJD(3);
t134 = sin(qJ(3));
t209 = qJD(2) * t134;
t97 = t136 * t209 + t203;
t233 = t97 * t119;
t205 = qJD(4) * t136;
t184 = t134 * t205;
t197 = qJD(3) * qJD(4);
t69 = qJD(2) * (t137 * t203 + t184) + t133 * t197;
t277 = -t69 + t233;
t105 = -t137 * pkin(3) - t134 * pkin(9) - pkin(2);
t220 = t136 * t137;
t122 = pkin(8) * t220;
t135 = sin(qJ(2));
t206 = qJD(4) * t133;
t131 = sin(pkin(6));
t212 = qJD(1) * t131;
t138 = cos(qJ(2));
t219 = t137 * t138;
t276 = (t133 * t219 - t135 * t136) * t212 - qJD(4) * t122 - t105 * t206;
t168 = pkin(3) * t134 - pkin(9) * t137;
t100 = t168 * qJD(3);
t275 = -(t133 * t135 + t136 * t219) * t212 + t133 * t100 + t105 * t205;
t198 = qJD(2) * qJD(3);
t181 = t134 * t198;
t259 = pkin(4) + pkin(5);
t274 = t259 * t181;
t273 = qJ(5) * t181 - t119 * qJD(5);
t112 = t119 * qJ(5);
t190 = t135 * t212;
t102 = qJD(2) * pkin(8) + t190;
t132 = cos(pkin(6));
t223 = t132 * t134;
t115 = qJD(1) * t223;
t74 = t137 * t102 + t115;
t61 = qJD(3) * pkin(9) + t74;
t183 = t138 * t212;
t75 = qJD(2) * t105 - t183;
t25 = t133 * t75 + t136 * t61;
t19 = -t112 + t25;
t175 = pkin(4) * t181;
t210 = qJD(2) * t131;
t188 = t138 * t210;
t208 = qJD(3) * t134;
t211 = qJD(1) * t137;
t39 = -t102 * t208 + (qJD(3) * t132 + t188) * t211;
t72 = (t100 + t190) * qJD(2);
t179 = -t133 * t39 + t136 * t72 - t61 * t205 - t75 * t206;
t4 = -t175 - t179;
t272 = t119 * t19 + t4;
t201 = t136 * qJD(3);
t95 = t133 * t209 - t201;
t271 = t69 * qJ(6) + t95 * qJD(6);
t73 = -t134 * t102 + t132 * t211;
t170 = qJD(3) * pkin(3) + t73;
t150 = qJ(5) * t97 + t170;
t15 = -t259 * t95 + qJD(6) + t150;
t270 = (qJD(6) + t15) * t97;
t189 = t135 * t210;
t172 = t137 * t188;
t226 = t131 * t135;
t82 = -t132 * t137 + t134 * t226;
t48 = -qJD(3) * t82 + t172;
t225 = t131 * t138;
t83 = t137 * t226 + t223;
t51 = -t133 * t225 + t83 * t136;
t13 = qJD(4) * t51 + t48 * t133 - t136 * t189;
t173 = t134 * t188;
t49 = qJD(3) * t83 + t173;
t50 = t83 * t133 + t136 * t225;
t269 = t13 * t119 - t50 * t181 + t49 * t95 + t82 * t69;
t93 = t97 ^ 2;
t268 = -t119 ^ 2 - t93;
t267 = qJ(5) * t208 - t137 * qJD(5) + t275;
t14 = -qJD(4) * t50 + t133 * t189 + t48 * t136;
t185 = t134 * t206;
t187 = t137 * t201;
t68 = -t136 * t197 + (t185 - t187) * qJD(2);
t266 = t14 * t119 - t51 * t181 + t49 * t97 - t68 * t82;
t265 = -t136 * t100 - t276;
t264 = -t133 * qJD(5) - t74;
t238 = t119 * t95;
t263 = -t68 + t238;
t24 = -t133 * t61 + t136 * t75;
t215 = qJD(5) - t24;
t262 = 0.2e1 * t273;
t222 = t133 * qJ(5);
t261 = -t259 * t136 - t222;
t260 = t95 ^ 2;
t191 = -pkin(8) * t133 - pkin(4);
t200 = t136 * qJD(6);
t207 = qJD(3) * t137;
t258 = -(-qJ(6) * t207 - t100) * t136 - (qJ(6) * t206 - t200 + (-pkin(5) + t191) * qJD(3)) * t134 + t276;
t257 = -(-pkin(8) * qJD(3) + qJ(6) * qJD(4)) * t136 * t134 - (qJD(6) * t134 + (-pkin(8) * qJD(4) + qJ(6) * qJD(3)) * t137) * t133 - t267;
t26 = pkin(4) * t95 - t150;
t256 = t26 * t97;
t171 = t134 * t183;
t40 = qJD(2) * t171 + qJD(3) * t115 + t102 * t207;
t6 = t69 * pkin(4) + t68 * qJ(5) - t97 * qJD(5) + t40;
t255 = t6 * t133;
t254 = t6 * t136;
t253 = t97 * t95;
t252 = pkin(9) - qJ(6);
t108 = t252 * t136;
t99 = t168 * qJD(2);
t180 = -t133 * t73 + t136 * t99;
t251 = (-qJ(6) * t220 - t259 * t134) * qJD(2) - t180 - qJD(4) * t108 + t133 * qJD(6);
t250 = (-t134 * t201 - t137 * t206) * pkin(8) + t267;
t245 = t133 * t99 + t136 * t73;
t28 = qJ(5) * t209 + t245;
t249 = qJ(6) * t133 * t199 + t206 * t252 + t200 + t28;
t248 = t191 * t208 + t265;
t230 = qJ(5) * t136;
t155 = -t259 * t133 + t230;
t247 = t119 * t155 + t264;
t165 = pkin(4) * t133 - t230;
t246 = t119 * t165 - t264;
t243 = qJ(5) * t69;
t242 = qJ(5) * t95;
t241 = qJD(2) * pkin(2);
t17 = qJ(6) * t95 + t25;
t12 = -t112 + t17;
t240 = t119 * t12;
t237 = t40 * t133;
t236 = t40 * t136;
t235 = t170 * t133;
t234 = t68 * t133;
t231 = t133 * t105 + t122;
t229 = qJ(6) * t134;
t228 = t119 * t133;
t227 = t119 * t136;
t141 = qJD(2) ^ 2;
t224 = t131 * t141;
t221 = t133 * t137;
t140 = qJD(3) ^ 2;
t218 = t140 * t134;
t217 = t140 * t137;
t16 = qJ(6) * t97 + t24;
t216 = qJD(5) - t16;
t129 = t134 ^ 2;
t213 = -t137 ^ 2 + t129;
t204 = qJD(5) * t136;
t196 = pkin(9) * t228;
t195 = pkin(9) * t227;
t194 = pkin(9) * t208;
t193 = pkin(9) * t201;
t192 = t135 * t224;
t186 = t119 * t206;
t121 = pkin(8) * t221;
t178 = t136 * t105 - t121;
t177 = t95 * t183;
t176 = t97 * t183;
t66 = -t137 * qJ(5) + t231;
t103 = -t183 - t241;
t167 = -t103 - t183;
t166 = t136 * pkin(4) + t222;
t18 = pkin(4) * t119 + t215;
t164 = -t133 * t19 + t136 * t18;
t163 = qJD(2) * t129 - t119 * t137;
t160 = pkin(8) + t165;
t5 = -pkin(5) * t69 - t6;
t159 = -t5 * t133 - t15 * t205;
t158 = t5 * t136 - t15 * t206;
t157 = t68 * qJ(6) - t179;
t156 = -t119 * t25 + t179;
t154 = -t181 + t253;
t153 = -t133 * t72 - t136 * t39 - t75 * t205 + t206 * t61;
t149 = -pkin(8) + t155;
t3 = -t153 + t273;
t147 = -t119 * t24 + t153;
t146 = qJD(3) * (-t167 - t241);
t145 = t13 * t97 - t14 * t95 - t50 * t68 - t51 * t69;
t144 = t68 + t238;
t143 = t157 - t274;
t128 = t137 * pkin(4);
t107 = t252 * t133;
t104 = -pkin(3) - t166;
t92 = pkin(3) - t261;
t77 = t160 * t134;
t67 = t128 - t178;
t65 = t149 * t134;
t45 = pkin(4) * t97 + t242;
t43 = t133 * t229 + t66;
t41 = t137 * pkin(5) + t121 + t128 + (-t105 - t229) * t136;
t33 = -t259 * t97 - t242;
t31 = (qJD(4) * t166 - t204) * t134 + t160 * t207;
t30 = -pkin(4) * t209 - t180;
t20 = (t261 * qJD(4) + t204) * t134 + t149 * t207;
t7 = t259 * t119 + t216;
t2 = -t97 * qJD(6) + t143;
t1 = t3 + t271;
t8 = [0, 0, -t192, -t138 * t224, 0, 0, 0, 0, 0, -t137 * t192 + (-t49 - t173) * qJD(3), t134 * t192 + (-t48 - t172) * qJD(3), 0, 0, 0, 0, 0, t269, t266, t269, t145, -t266, t13 * t18 + t14 * t19 + t26 * t49 + t3 * t51 + t4 * t50 + t6 * t82, t269, -t266, -t145, t1 * t51 + t12 * t14 + t13 * t7 - t15 * t49 + t2 * t50 - t5 * t82; 0, 0, 0, 0, 0.2e1 * t137 * t181, -0.2e1 * t213 * t198, t217, -t218, 0, -pkin(8) * t217 + t134 * t146, pkin(8) * t218 + t137 * t146, t97 * t187 + (-t68 * t136 - t206 * t97) * t134 (-t133 * t97 - t136 * t95) * t207 + (t234 - t136 * t69 + (t133 * t95 - t136 * t97) * qJD(4)) * t134, t119 * t185 + t68 * t137 + (t134 * t97 + t136 * t163) * qJD(3), t119 * t184 + t69 * t137 + (-t133 * t163 - t134 * t95) * qJD(3) (-t119 - t199) * t208, t265 * t119 + ((pkin(8) * t95 - t235) * qJD(3) - t179) * t137 + (-t177 - t170 * t205 + pkin(8) * t69 + t237 + (-pkin(8) * t228 + qJD(2) * t178 + t24) * qJD(3)) * t134, t275 * t119 + (-t170 * t201 + (qJD(3) * t97 - t186) * pkin(8) - t153) * t137 + (-t176 + t170 * t206 - pkin(8) * t68 + t236 + (-pkin(8) * t227 - qJD(2) * t231 - t25) * qJD(3)) * t134, t31 * t95 + t77 * t69 + (t203 * t26 + t4) * t137 + t248 * t119 + (-t177 + t26 * t205 + t255 + (-qJD(2) * t67 - t18) * qJD(3)) * t134, -t66 * t69 - t67 * t68 + t248 * t97 - t250 * t95 + t164 * t207 + (-t133 * t3 + t136 * t4 + (-t133 * t18 - t136 * t19) * qJD(4)) * t134, -t31 * t97 + t77 * t68 + (-t201 * t26 - t3) * t137 - t250 * t119 + (t176 + t26 * t206 - t254 + (qJD(2) * t66 + t19) * qJD(3)) * t134, t3 * t66 + t4 * t67 + t6 * t77 + (t31 - t171) * t26 + t250 * t19 + t248 * t18, -t20 * t95 - t65 * t69 + (-t15 * t203 + t2) * t137 - t258 * t119 + (-t177 + (-qJD(2) * t41 - t7) * qJD(3) + t159) * t134, t20 * t97 - t65 * t68 + (t15 * t201 - t1) * t137 + t257 * t119 + (t176 + (qJD(2) * t43 + t12) * qJD(3) + t158) * t134, t41 * t68 + t43 * t69 + t258 * t97 - t257 * t95 + (t12 * t133 - t136 * t7) * t207 + (t1 * t133 - t136 * t2 + (t12 * t136 + t133 * t7) * qJD(4)) * t134, t1 * t43 + t2 * t41 + t5 * t65 - t258 * t7 + (t20 + t171) * t15 - t257 * t12; 0, 0, 0, 0, -t134 * t141 * t137, t213 * t141, 0, 0, 0, qJD(3) * t74 - t103 * t209 - t40, t167 * t199, -t227 * t97 - t234, t133 * t277 + t263 * t136, -t119 * t205 + (t119 * t220 + (-t97 + t203) * t134) * qJD(2), t186 + (-t119 * t221 + (t95 + t201) * t134) * qJD(2), t119 * t209, -pkin(3) * t69 - t236 + t180 * t119 - t74 * t95 + (t195 - t235) * qJD(4) + (-t24 * t134 + (t137 * t170 - t194) * t133) * qJD(2), pkin(3) * t68 + t237 - t245 * t119 - t74 * t97 + (-t136 * t170 - t196) * qJD(4) + (t170 * t220 + (t25 - t193) * t134) * qJD(2), t104 * t69 - t30 * t119 - t254 - t246 * t95 + (t133 * t26 + t195) * qJD(4) + (t134 * t18 + (-t137 * t26 - t194) * t133) * qJD(2), t28 * t95 - t30 * t97 + (t3 - t119 * t18 + (qJD(4) * t97 - t69) * pkin(9)) * t136 + ((qJD(4) * t95 - t68) * pkin(9) + t272) * t133, t104 * t68 + t28 * t119 - t255 + t246 * t97 + (-t136 * t26 + t196) * qJD(4) + (t26 * t220 + (-t19 + t193) * t134) * qJD(2), t6 * t104 - t18 * t30 - t19 * t28 - t246 * t26 + (qJD(4) * t164 + t4 * t133 + t3 * t136) * pkin(9), -t92 * t69 + t247 * t95 - t251 * t119 + (t15 * t221 + (-qJD(3) * t107 + t7) * t134) * qJD(2) + t158, -t92 * t68 - t247 * t97 + t249 * t119 + (-t15 * t220 + (qJD(3) * t108 - t12) * t134) * qJD(2) - t159, t107 * t68 + t108 * t69 + t251 * t97 - t249 * t95 + (t119 * t7 - t1) * t136 + (-t2 - t240) * t133, t1 * t108 + t107 * t2 - t12 * t249 - t15 * t247 - t251 * t7 + t5 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, t93 - t260, -t144, -t69 - t233, t181, t170 * t97 + t156, -t170 * t95 + t147, -t45 * t95 + t156 + 0.2e1 * t175 - t256, pkin(4) * t68 - t243 + (t19 - t25) * t97 + (t18 - t215) * t95, -t26 * t95 + t45 * t97 - t147 + t262, -pkin(4) * t4 + qJ(5) * t3 - t18 * t25 + t19 * t215 - t26 * t45, -t17 * t119 + t33 * t95 - t157 + t270 + 0.2e1 * t274, t119 * t16 + t15 * t95 - t33 * t97 - t153 + t262 + t271, t243 - t259 * t68 + (-t12 + t17) * t97 + (-t7 + t216) * t95, t1 * qJ(5) + t12 * t216 - t15 * t33 - t7 * t17 - t2 * t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, -t144, t268, t256 + t272, t154, t268, t144, t143 + t240 - t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t277, t263, -t93 - t260, -t12 * t95 + t7 * t97 + t5;];
tauc_reg  = t8;
