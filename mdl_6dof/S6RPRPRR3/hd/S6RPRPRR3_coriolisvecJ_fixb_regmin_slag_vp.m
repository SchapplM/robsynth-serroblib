% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:42:26
% EndTime: 2019-03-09 03:42:36
% DurationCPUTime: 3.41s
% Computational Cost: add. (4119->353), mult. (10021->511), div. (0->0), fcn. (7441->10), ass. (0->188)
t192 = cos(qJ(3));
t250 = qJD(1) * t192;
t172 = -qJD(5) + t250;
t169 = -qJD(6) + t172;
t187 = sin(qJ(6));
t190 = cos(qJ(6));
t183 = sin(pkin(11));
t185 = cos(pkin(11));
t239 = t185 * qJD(3);
t189 = sin(qJ(3));
t251 = qJD(1) * t189;
t145 = t183 * t251 - t239;
t248 = qJD(3) * t183;
t147 = t185 * t251 + t248;
t188 = sin(qJ(5));
t191 = cos(qJ(5));
t92 = t145 * t188 - t147 * t191;
t93 = t191 * t145 + t147 * t188;
t276 = t187 * t92 - t190 * t93;
t287 = t169 * t276;
t235 = t183 * t250;
t219 = pkin(3) * t189 - qJ(4) * t192;
t157 = t219 * qJD(1);
t174 = sin(pkin(10)) * pkin(1) + pkin(7);
t164 = t174 * qJD(1);
t274 = qJD(2) * t192 - t189 * t164;
t85 = t183 * t157 + t185 * t274;
t67 = -pkin(8) * t235 + t85;
t286 = qJD(4) * t185 - t67;
t215 = t187 * t93 + t190 * t92;
t285 = t215 * t276;
t284 = t92 * t172;
t283 = t93 * t172;
t282 = t169 * t215;
t259 = t191 * t185;
t264 = t183 * t188;
t152 = -t259 + t264;
t204 = t152 * t192;
t256 = qJD(1) * t204 - t152 * qJD(5);
t153 = t183 * t191 + t185 * t188;
t205 = t153 * t192;
t255 = -qJD(1) * t205 + t153 * qJD(5);
t281 = t215 ^ 2 - t276 ^ 2;
t180 = t189 * qJD(2);
t131 = t192 * t164 + t180;
t115 = qJD(3) * qJ(4) + t131;
t175 = -cos(pkin(10)) * pkin(1) - pkin(2);
t141 = -pkin(3) * t192 - qJ(4) * t189 + t175;
t118 = t141 * qJD(1);
t57 = -t115 * t183 + t185 * t118;
t37 = -pkin(4) * t250 - pkin(8) * t147 + t57;
t58 = t185 * t115 + t183 * t118;
t40 = -pkin(8) * t145 + t58;
t16 = t188 * t37 + t191 * t40;
t13 = -pkin(9) * t93 + t16;
t241 = qJD(6) * t187;
t11 = t13 * t241;
t224 = -qJD(3) * pkin(3) + qJD(4);
t112 = -t274 + t224;
t86 = pkin(4) * t145 + t112;
t36 = pkin(5) * t93 + t86;
t280 = -t276 * t36 + t11;
t240 = qJD(6) * t190;
t238 = qJD(1) * qJD(3);
t231 = t192 * t238;
t222 = t183 * t231;
t242 = qJD(5) * t191;
t48 = -t145 * t242 + t231 * t259 + (-qJD(5) * t147 - t222) * t188;
t200 = qJD(3) * t205;
t275 = qJD(5) * t92;
t49 = qJD(1) * t200 - t275;
t8 = -t187 * t49 + t190 * t48 - t93 * t240 + t241 * t92;
t279 = t8 + t287;
t15 = -t188 * t40 + t191 * t37;
t12 = pkin(9) * t92 + t15;
t10 = -pkin(5) * t172 + t12;
t266 = t13 * t190;
t2 = t10 * t187 + t266;
t178 = t189 * t238;
t260 = t185 * t192;
t210 = pkin(4) * t189 - pkin(8) * t260;
t202 = t210 * qJD(3);
t109 = (qJD(4) + t274) * qJD(3);
t134 = qJD(3) * t219 - t189 * qJD(4);
t119 = t134 * qJD(1);
t52 = -t183 * t109 + t185 * t119;
t38 = qJD(1) * t202 + t52;
t53 = t185 * t109 + t183 * t119;
t41 = -pkin(8) * t222 + t53;
t228 = -t188 * t41 + t191 * t38;
t196 = -qJD(5) * t16 + t228;
t4 = pkin(5) * t178 - t48 * pkin(9) + t196;
t244 = qJD(5) * t188;
t207 = t188 * t38 + t191 * t41 + t37 * t242 - t244 * t40;
t5 = -pkin(9) * t49 + t207;
t236 = -t187 * t5 + t190 * t4;
t278 = -qJD(6) * t2 + t36 * t215 + t236;
t195 = qJD(6) * t215 - t187 * t48 - t190 * t49;
t277 = t195 + t282;
t272 = pkin(8) + qJ(4);
t161 = t272 * t183;
t162 = t272 * t185;
t254 = -t188 * t161 + t191 * t162;
t211 = qJD(4) * t183 + qJD(5) * t162;
t84 = t185 * t157 - t183 * t274;
t54 = qJD(1) * t210 + t84;
t273 = -t161 * t242 + t286 * t191 + (-t211 - t54) * t188;
t126 = t153 * t189;
t127 = t152 * t189;
t73 = -t126 * t187 - t127 * t190;
t243 = qJD(5) * t189;
t79 = -qJD(3) * t204 - t153 * t243;
t261 = t185 * t189;
t80 = t242 * t261 - t243 * t264 + t200;
t18 = qJD(6) * t73 + t187 * t79 + t190 * t80;
t72 = t190 * t126 - t127 * t187;
t271 = t18 * t169 - t72 * t178;
t97 = t190 * t152 + t153 * t187;
t270 = -qJD(6) * t97 - t255 * t187 + t256 * t190;
t98 = -t152 * t187 + t153 * t190;
t269 = qJD(6) * t98 + t256 * t187 + t255 * t190;
t129 = t185 * t141;
t75 = -pkin(8) * t261 + t129 + (-t174 * t183 - pkin(4)) * t192;
t156 = t174 * t260;
t100 = t183 * t141 + t156;
t263 = t183 * t189;
t83 = -pkin(8) * t263 + t100;
t267 = t188 * t75 + t191 * t83;
t265 = -t126 * t178 + t80 * t172;
t262 = t183 * t192;
t166 = t189 * t174;
t193 = qJD(3) ^ 2;
t258 = t193 * t189;
t257 = t193 * t192;
t246 = qJD(3) * t192;
t120 = qJD(3) * t180 + t164 * t246;
t247 = qJD(3) * t189;
t234 = t174 * t247;
t90 = t185 * t134 + t183 * t234;
t124 = (pkin(4) * t183 + t174) * t246;
t133 = pkin(4) * t263 + t166;
t181 = t189 ^ 2;
t253 = -t192 ^ 2 + t181;
t165 = qJD(1) * t175;
t252 = qJD(1) * t181;
t237 = t174 * t262;
t101 = pkin(4) * t222 + t120;
t103 = pkin(4) * t235 + t131;
t176 = -pkin(4) * t185 - pkin(3);
t233 = t255 * pkin(5) - t103;
t232 = -t192 * t8 - t215 * t247;
t230 = qJD(6) * t10 + t5;
t65 = t202 + t90;
t122 = t183 * t134;
t74 = t122 + (-pkin(8) * t262 - t185 * t166) * qJD(3);
t227 = -t188 * t74 + t191 * t65;
t226 = -t188 * t83 + t191 * t75;
t225 = -t192 * t48 - t247 * t92;
t223 = -t191 * t161 - t162 * t188;
t51 = t191 * t54;
t77 = -pkin(9) * t152 + t254;
t221 = pkin(5) * t251 + t256 * pkin(9) + t153 * qJD(4) + t254 * qJD(5) + qJD(6) * t77 - t188 * t67 + t51;
t76 = -pkin(9) * t153 + t223;
t220 = -t255 * pkin(9) + qJD(6) * t76 + t273;
t218 = -t52 * t183 + t53 * t185;
t217 = -t183 * t57 + t185 * t58;
t21 = -pkin(5) * t192 + pkin(9) * t127 + t226;
t22 = -pkin(9) * t126 + t267;
t216 = t187 * t21 + t190 * t22;
t212 = 0.2e1 * qJD(3) * t165;
t209 = -t192 * t195 + t247 * t276;
t208 = t192 * t49 - t247 * t93;
t206 = t188 * t65 + t191 * t74 + t75 * t242 - t244 * t83;
t17 = -qJD(6) * t72 - t187 * t80 + t190 * t79;
t203 = t169 * t17 - t178 * t73;
t201 = t127 * t178 + t172 * t79;
t197 = -qJ(4) * t247 + (-t112 + t224) * t192;
t194 = qJD(1) ^ 2;
t121 = pkin(5) * t152 + t176;
t99 = t129 - t237;
t91 = -t185 * t234 + t122;
t89 = pkin(5) * t126 + t133;
t45 = pkin(5) * t80 + t124;
t23 = pkin(5) * t49 + t101;
t7 = -pkin(9) * t80 + t206;
t6 = pkin(5) * t247 - t79 * pkin(9) - qJD(5) * t267 + t227;
t1 = t10 * t190 - t13 * t187;
t3 = [0, 0, 0, 0, 0.2e1 * t192 * t178, -0.2e1 * t253 * t238, t257, -t258, 0, -t174 * t257 + t189 * t212, t174 * t258 + t192 * t212, t120 * t263 + (-qJD(1) * t90 - t52) * t192 + ((t112 * t183 + t145 * t174) * t192 + (t57 + (t99 + t237) * qJD(1)) * t189) * qJD(3), t120 * t261 + (qJD(1) * t91 + t53) * t192 + ((t112 * t185 + t147 * t174) * t192 + (-t58 + (-t100 + t156) * qJD(1)) * t189) * qJD(3), -t91 * t145 - t90 * t147 + (-t183 * t53 - t185 * t52) * t189 + (-t183 * t58 - t185 * t57 + (-t100 * t183 - t185 * t99) * qJD(1)) * t246, t53 * t100 + t52 * t99 + t57 * t90 + t58 * t91 + (t112 * t246 + t120 * t189) * t174, -t127 * t48 - t79 * t92, -t126 * t48 + t127 * t49 - t79 * t93 + t80 * t92, -t201 + t225, t208 + t265 (-t172 - t250) * t247, -t227 * t172 - t228 * t192 + t124 * t93 + t133 * t49 + t101 * t126 + t86 * t80 + (t16 * t192 + t172 * t267) * qJD(5) + (qJD(1) * t226 + t15) * t247, t206 * t172 + t207 * t192 - t124 * t92 + t133 * t48 - t101 * t127 + t86 * t79 + (-t267 * qJD(1) - t16) * t247, -t17 * t215 + t73 * t8, t17 * t276 + t18 * t215 + t195 * t73 - t72 * t8, -t203 + t232, t209 + t271 (-t169 - t250) * t247 -(-t187 * t7 + t190 * t6) * t169 - t236 * t192 - t45 * t276 - t89 * t195 + t23 * t72 + t36 * t18 + (t169 * t216 + t192 * t2) * qJD(6) + ((-t187 * t22 + t190 * t21) * qJD(1) + t1) * t247, -t11 * t192 + t36 * t17 + t23 * t73 - t45 * t215 + t89 * t8 + ((-qJD(6) * t22 + t6) * t169 + t4 * t192) * t187 + ((qJD(6) * t21 + t7) * t169 + t230 * t192) * t190 + (-qJD(1) * t216 - t2) * t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t258, -t257 (t145 * t189 - t183 * t252) * qJD(3) (t147 * t189 - t185 * t252) * qJD(3) (-t145 * t185 + t147 * t183) * t246, -t120 * t192 + t218 * t189 + (t112 * t189 + t192 * t217) * qJD(3), 0, 0, 0, 0, 0, -t208 + t265, t201 + t225, 0, 0, 0, 0, 0, -t209 + t271, t203 + t232; 0, 0, 0, 0, -t189 * t194 * t192, t253 * t194, 0, 0, 0, qJD(3) * t131 - t165 * t251 - t120, -t165 * t250, -t120 * t185 - t131 * t145 + (t183 * t197 - t189 * t57 + t192 * t84) * qJD(1), t120 * t183 - t131 * t147 + (t185 * t197 + t189 * t58 - t192 * t85) * qJD(1), t85 * t145 + t84 * t147 + (-qJD(4) * t145 + t250 * t57 + t53) * t185 + (qJD(4) * t147 + t250 * t58 - t52) * t183, -t120 * pkin(3) + qJ(4) * t218 + qJD(4) * t217 - t112 * t131 - t57 * t84 - t58 * t85, t48 * t153 - t256 * t92, -t48 * t152 - t153 * t49 + t255 * t92 - t256 * t93, -t256 * t172 + (qJD(3) * t153 + t92) * t251, t255 * t172 + (-qJD(3) * t152 + t93) * t251, t172 * t251, t101 * t152 - t103 * t93 + t176 * t49 + t255 * t86 + (t51 + t211 * t191 + (-qJD(5) * t161 + t286) * t188) * t172 + (qJD(3) * t223 - t15) * t251, t101 * t153 + t103 * t92 + t176 * t48 + t256 * t86 + t273 * t172 + (-t254 * qJD(3) + t16) * t251, -t215 * t270 + t8 * t98, t195 * t98 + t215 * t269 + t270 * t276 - t8 * t97, -t270 * t169 + (qJD(3) * t98 + t215) * t251, t269 * t169 + (-qJD(3) * t97 - t276) * t251, t169 * t251, -t121 * t195 + t23 * t97 + t269 * t36 - t233 * t276 + (t187 * t220 + t190 * t221) * t169 + ((-t187 * t77 + t190 * t76) * qJD(3) - t1) * t251, t121 * t8 + t23 * t98 + t270 * t36 - t233 * t215 + (-t187 * t221 + t190 * t220) * t169 + (-(t187 * t76 + t190 * t77) * qJD(3) + t2) * t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t147 + t248) * t250 (t145 + t239) * t250, -t145 ^ 2 - t147 ^ 2, t145 * t58 + t147 * t57 + t120, 0, 0, 0, 0, 0, t49 + t284, t48 + t283, 0, 0, 0, 0, 0, -t195 + t282, t8 - t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92 * t93, t92 ^ 2 - t93 ^ 2, t48 - t283, -t153 * t231 + t275 + t284, t178, -t16 * t172 + t86 * t92 + t196, -t15 * t172 + t86 * t93 - t207, t285, t281, t279, t277, t178 (-t12 * t187 - t266) * t169 + (t169 * t241 + t178 * t190 - t276 * t92) * pkin(5) + t278 (t13 * t169 - t4) * t187 + (-t12 * t169 - t230) * t190 + (t169 * t240 - t178 * t187 - t215 * t92) * pkin(5) + t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, t281, t279, t277, t178, -t2 * t169 + t278, -t1 * t169 - t187 * t4 - t190 * t230 + t280;];
tauc_reg  = t3;
