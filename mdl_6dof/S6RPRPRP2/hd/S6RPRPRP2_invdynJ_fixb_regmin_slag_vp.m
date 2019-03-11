% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:16
% EndTime: 2019-03-09 03:06:22
% DurationCPUTime: 2.79s
% Computational Cost: add. (5110->407), mult. (10858->504), div. (0->0), fcn. (7590->14), ass. (0->204)
t149 = sin(pkin(9));
t128 = pkin(1) * t149 + pkin(7);
t228 = qJ(4) + t128;
t156 = cos(qJ(3));
t236 = cos(pkin(10));
t204 = t236 * t156;
t123 = qJD(1) * t204;
t148 = sin(pkin(10));
t153 = sin(qJ(3));
t222 = qJD(1) * t153;
t99 = -t148 * t222 + t123;
t96 = qJD(5) - t99;
t150 = cos(pkin(9));
t130 = -t150 * pkin(1) - pkin(2);
t277 = qJDD(1) * t130;
t110 = t148 * t156 + t236 * t153;
t100 = t110 * qJD(3);
t178 = -t148 * t153 + t204;
t152 = sin(qJ(5));
t155 = cos(qJ(5));
t218 = qJD(1) * qJD(3);
t208 = t153 * t218;
t269 = qJD(3) * t123 + qJDD(1) * t110 - t148 * t208;
t161 = -t155 * qJDD(3) + t152 * t269;
t101 = t110 * qJD(1);
t84 = qJD(3) * t152 + t101 * t155;
t33 = qJD(5) * t84 + t161;
t82 = -t155 * qJD(3) + t101 * t152;
t276 = -t100 * t82 + t178 * t33;
t219 = qJD(5) * t155;
t216 = t153 * qJDD(1);
t189 = qJDD(1) * t204 - t148 * t216;
t73 = qJD(1) * t100 + qJDD(5) - t189;
t66 = t152 * t73;
t275 = -t96 * t219 - t66;
t220 = qJD(5) * t152;
t211 = t96 * t220;
t67 = t155 * t73;
t274 = t211 - t67;
t144 = qJ(3) + pkin(10);
t133 = sin(t144);
t145 = qJ(1) + pkin(9);
t134 = sin(t145);
t136 = cos(t145);
t197 = g(1) * t136 + g(2) * t134;
t273 = t133 * t197;
t206 = -g(1) * t134 + g(2) * t136;
t142 = t156 * pkin(3);
t272 = t130 - t142;
t135 = cos(t144);
t262 = pkin(8) * t133;
t271 = pkin(4) * t135 + t262;
t201 = t228 * qJD(1);
t87 = qJD(2) * t153 + t156 * t201;
t77 = t148 * t87;
t244 = qJD(3) * pkin(3);
t86 = t156 * qJD(2) - t153 * t201;
t80 = t86 + t244;
t42 = t236 * t80 - t77;
t36 = -qJD(3) * pkin(4) - t42;
t19 = t82 * pkin(5) - t84 * qJ(6) + t36;
t264 = pkin(3) * t148;
t127 = pkin(8) + t264;
t241 = t127 * t73;
t270 = t96 * t19 - t241;
t268 = t84 ^ 2;
t267 = t96 ^ 2;
t266 = pkin(5) * t73;
t209 = t236 * t87;
t43 = t148 * t80 + t209;
t37 = qJD(3) * pkin(8) + t43;
t98 = qJD(1) * t272 + qJD(4);
t54 = -pkin(4) * t99 - pkin(8) * t101 + t98;
t18 = t152 * t54 + t155 * t37;
t8 = qJ(6) * t96 + t18;
t265 = t8 * t96;
t258 = g(3) * t133;
t257 = g(3) * t135;
t256 = g(3) * t156;
t254 = t18 * t96;
t253 = t82 * t99;
t252 = t84 * t82;
t205 = t84 * t96;
t103 = t178 * qJD(3);
t229 = t155 * t103;
t235 = t110 * t155;
t251 = -t82 * t229 - t33 * t235;
t250 = -t152 * t33 - t82 * t219;
t139 = t156 * qJDD(2);
t114 = t128 * qJDD(1);
t172 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(3) * qJD(2) + t114;
t185 = t201 * qJD(3);
t40 = qJDD(3) * pkin(3) - t153 * t172 - t156 * t185 + t139;
t44 = (qJDD(2) - t185) * t153 + t172 * t156;
t14 = t148 * t40 + t236 * t44;
t46 = t236 * t86 - t77;
t60 = pkin(3) * t222 + pkin(4) * t101 - pkin(8) * t99;
t249 = t152 * t60 + t155 * t46;
t217 = qJD(3) * qJD(5);
t32 = -t152 * qJDD(3) + t101 * t220 + (-t269 - t217) * t155;
t248 = t84 * t100 + t178 * t32;
t65 = -pkin(4) * t178 - pkin(8) * t110 + t272;
t105 = t228 * t153;
t106 = t228 * t156;
t70 = -t148 * t105 + t236 * t106;
t247 = t152 * t65 + t155 * t70;
t190 = pkin(5) * t152 - qJ(6) * t155;
t45 = t148 * t86 + t209;
t246 = -qJD(6) * t152 + t190 * t96 - t45;
t245 = qJ(6) * t73;
t243 = t101 * t82;
t242 = t101 * t84;
t240 = t152 * t32;
t239 = t152 * t82;
t238 = t152 * t96;
t237 = t155 * t84;
t234 = t134 * t152;
t233 = t134 * t155;
t232 = t136 * t152;
t231 = t136 * t155;
t230 = t152 * t103;
t17 = -t152 * t37 + t155 * t54;
t227 = qJD(6) - t17;
t226 = qJDD(2) - g(3);
t225 = (g(1) * t231 + g(2) * t233) * t133;
t132 = t142 + pkin(2);
t157 = cos(qJ(1));
t224 = t157 * pkin(1) + t136 * t132;
t146 = t153 ^ 2;
t223 = -t156 ^ 2 + t146;
t117 = qJD(1) * t130;
t221 = qJD(5) * t127;
t215 = t156 * qJDD(1);
t214 = pkin(3) * t208 + qJDD(4);
t213 = t153 * t244;
t212 = t84 * t230;
t72 = t96 * t229;
t210 = t236 * pkin(3);
t12 = qJDD(3) * pkin(8) + t14;
t74 = -qJD(3) * t101 + t189;
t26 = -pkin(3) * t215 - t74 * pkin(4) - pkin(8) * t269 + t214 + t277;
t207 = t152 * t12 - t155 * t26 + t37 * t219 + t54 * t220;
t202 = qJD(3) * t228;
t88 = qJD(4) * t156 - t153 * t202;
t89 = -qJD(4) * t153 - t156 * t202;
t52 = t148 * t88 - t236 * t89;
t69 = t236 * t105 + t106 * t148;
t91 = t135 * t234 + t231;
t93 = t135 * t232 - t233;
t199 = -g(1) * t91 + g(2) * t93;
t92 = t135 * t233 - t232;
t94 = t135 * t231 + t234;
t198 = g(1) * t92 - g(2) * t94;
t154 = sin(qJ(1));
t195 = g(1) * t154 - g(2) * t157;
t7 = -pkin(5) * t96 + t227;
t194 = -t152 * t8 + t155 * t7;
t193 = t152 * t7 + t155 * t8;
t13 = -t148 * t44 + t236 * t40;
t151 = -qJ(4) - pkin(7);
t192 = -pkin(1) * t154 - t136 * t151;
t191 = t155 * pkin(5) + t152 * qJ(6);
t188 = -t155 * t96 * t99 - t275;
t186 = t99 * t238 - t274;
t184 = pkin(4) + t191;
t183 = t221 * t96 + t257;
t182 = t155 * t12 + t152 * t26 + t54 * t219 - t220 * t37;
t53 = t148 * t89 + t236 * t88;
t61 = pkin(4) * t100 - pkin(8) * t103 + t213;
t181 = t152 * t61 + t155 * t53 + t65 * t219 - t220 * t70;
t180 = -t110 * t211 + t73 * t235 + t72;
t179 = t96 * t36 - t241;
t11 = -qJDD(3) * pkin(4) - t13;
t3 = t33 * pkin(5) + t32 * qJ(6) - t84 * qJD(6) + t11;
t177 = -t183 - t3;
t175 = -t117 * qJD(1) - t114 + t197;
t174 = g(1) * t93 + g(2) * t91 + t152 * t258 - t207;
t173 = 0.2e1 * qJD(3) * t117 - qJDD(3) * t128;
t171 = qJDD(1) * t272 + t214;
t158 = qJD(3) ^ 2;
t170 = -t128 * t158 - t206 - 0.2e1 * t277;
t1 = qJD(6) * t96 + t182 + t245;
t2 = qJDD(6) + t207 - t266;
t169 = qJD(5) * t194 + t1 * t155 + t2 * t152;
t168 = t110 * t275 - t96 * t230;
t167 = t19 * t84 + qJDD(6) - t174;
t165 = t168 - t276;
t164 = -g(1) * t94 - g(2) * t92 - t155 * t258 + t182;
t159 = qJD(1) ^ 2;
t129 = -t210 - pkin(4);
t113 = qJDD(3) * t156 - t153 * t158;
t112 = qJDD(3) * t153 + t156 * t158;
t104 = -t210 - t184;
t50 = pkin(5) * t84 + qJ(6) * t82;
t31 = t110 * t190 + t69;
t24 = pkin(5) * t178 + t152 * t70 - t155 * t65;
t21 = -qJ(6) * t178 + t247;
t20 = t82 * t96 - t32;
t16 = -pkin(5) * t101 + t152 * t46 - t155 * t60;
t15 = qJ(6) * t101 + t249;
t6 = t190 * t103 + (qJD(5) * t191 - qJD(6) * t155) * t110 + t52;
t5 = -pkin(5) * t100 + t247 * qJD(5) + t152 * t53 - t155 * t61;
t4 = qJ(6) * t100 - qJD(6) * t178 + t181;
t9 = [qJDD(1), t195, g(1) * t157 + g(2) * t154 (t195 + (t149 ^ 2 + t150 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t146 + 0.2e1 * t156 * t208, 0.2e1 * t153 * t215 - 0.2e1 * t223 * t218, t112, t113, 0, t153 * t173 + t156 * t170, -t153 * t170 + t156 * t173, -t43 * t100 + t52 * t101 - t42 * t103 - t13 * t110 + t14 * t178 + t269 * t69 + t53 * t99 + t70 * t74 - t197, t14 * t70 + t43 * t53 - t13 * t69 - t42 * t52 + t171 * t272 + t98 * t213 - g(1) * (-t132 * t134 + t192) - g(2) * (-t134 * t151 + t224) t84 * t229 + (-t155 * t32 - t220 * t84) * t110, -t212 + (t240 + (-t237 + t239) * qJD(5)) * t110 + t251, t180 + t248, t168 + t276, t100 * t96 - t178 * t73, t207 * t178 + t17 * t100 + t52 * t82 + t69 * t33 + ((-qJD(5) * t70 + t61) * t96 + t65 * t73 + t36 * qJD(5) * t110) * t155 + ((-qJD(5) * t65 - t53) * t96 - t70 * t73 + t11 * t110 + t36 * t103) * t152 + t198, -t181 * t96 - t247 * t73 + t182 * t178 - t18 * t100 + t52 * t84 - t69 * t32 + t36 * t229 + (t11 * t155 - t220 * t36) * t110 + t199, t19 * t230 - t100 * t7 + t178 * t2 - t24 * t73 + t31 * t33 - t5 * t96 + t6 * t82 + (t152 * t3 + t19 * t219) * t110 + t198, -t21 * t33 - t24 * t32 - t4 * t82 + t5 * t84 - t206 * t133 + t194 * t103 + (-qJD(5) * t193 - t1 * t152 + t155 * t2) * t110, -t19 * t229 - t1 * t178 + t100 * t8 + t21 * t73 + t31 * t32 + t4 * t96 - t6 * t84 + (-t155 * t3 + t19 * t220) * t110 - t199, t1 * t21 + t8 * t4 + t3 * t31 + t19 * t6 + t2 * t24 + t7 * t5 - g(1) * (-pkin(5) * t92 - qJ(6) * t91 + t192) - g(2) * (pkin(5) * t94 + qJ(6) * t93 + t136 * t271 + t224) + (-g(1) * (-t132 - t271) + g(2) * t151) * t134; 0, 0, 0, t226, 0, 0, 0, 0, 0, t113, -t112, t100 * t101 + t103 * t99 + t110 * t74 - t178 * t269, -t100 * t42 + t103 * t43 + t110 * t14 + t13 * t178 - g(3), 0, 0, 0, 0, 0, t165, t110 * t274 + t248 - t72, t165, t212 + (-t240 + (t237 + t239) * qJD(5)) * t110 + t251, t180 - t248, t100 * t19 + t103 * t193 + t110 * t169 - t178 * t3 - g(3); 0, 0, 0, 0, -t153 * t159 * t156, t223 * t159, t216, t215, qJDD(3), t153 * t175 + t139 - t256, -t226 * t153 + t175 * t156, -t269 * t210 + t74 * t264 + (-t46 + t42) * t99 - (-t43 + t45) * t101, t42 * t45 - t43 * t46 + (t236 * t13 - t256 + t14 * t148 + (-qJD(1) * t98 + t197) * t153) * pkin(3), t155 * t205 - t240 (-t32 + t253) * t155 - t84 * t238 + t250, t188 - t242, t186 + t243, -t96 * t101, -t17 * t101 + t129 * t33 - t45 * t82 + (-t257 - t11 + (-t60 - t221) * t96) * t155 + (t46 * t96 + t179) * t152 + t225, -t129 * t32 + t249 * t96 + t18 * t101 - t45 * t84 + t179 * t155 + (t11 + t183 - t273) * t152, t101 * t7 + t104 * t33 + t152 * t270 + t177 * t155 + t16 * t96 + t246 * t82 + t225, -t258 + t15 * t82 - t16 * t84 - t197 * t135 + (-t127 * t33 - t7 * t99 + t1 + (t127 * t84 + t7) * qJD(5)) * t155 + (-t127 * t32 + t8 * t99 + t2 + (t127 * t82 - t8) * qJD(5)) * t152, -t101 * t8 + t104 * t32 - t15 * t96 - t246 * t84 - t270 * t155 + (t177 + t273) * t152, t3 * t104 - t8 * t15 - t7 * t16 - g(3) * (t142 + t262) + t246 * t19 - t184 * t257 + t169 * t127 + t197 * (pkin(3) * t153 - pkin(8) * t135 + t133 * t184); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 ^ 2 - t99 ^ 2, t101 * t42 - t43 * t99 + t171 + t206, 0, 0, 0, 0, 0, t186 - t243, -t267 * t155 - t242 - t66, -t238 * t96 - t243 + t67 (t32 + t253) * t155 + t152 * t205 + t250, t188 + t242, -t101 * t19 + (-t2 + t265) * t155 + (t96 * t7 + t1) * t152 + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, -t82 ^ 2 + t268, t20, -t101 * t219 - t152 * t217 - t161 + t205, t73, -t36 * t84 + t174 + t254, t17 * t96 + t36 * t82 - t164, -t50 * t82 - t167 + t254 + 0.2e1 * t266, pkin(5) * t32 - qJ(6) * t33 + (-t18 + t8) * t84 + (t7 - t227) * t82, 0.2e1 * t245 - t19 * t82 + t50 * t84 + (0.2e1 * qJD(6) - t17) * t96 + t164, t1 * qJ(6) - t2 * pkin(5) - t19 * t50 - t7 * t18 - g(1) * (-pkin(5) * t93 + qJ(6) * t94) - g(2) * (-pkin(5) * t91 + qJ(6) * t92) + t227 * t8 + t190 * t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(5) + t74 + t252, t20, -t267 - t268, t167 - t265 - t266;];
tau_reg  = t9;
