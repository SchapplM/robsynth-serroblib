% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:19
% EndTime: 2020-01-03 12:15:27
% DurationCPUTime: 2.71s
% Computational Cost: add. (6965->324), mult. (11818->433), div. (0->0), fcn. (7995->8), ass. (0->212)
t198 = qJD(1) + qJD(2);
t202 = sin(qJ(4));
t203 = sin(qJ(3));
t266 = t202 * t203;
t242 = t198 * t266;
t206 = cos(qJ(4));
t207 = cos(qJ(3));
t263 = t206 * t207;
t139 = -t198 * t263 + t242;
t162 = t202 * t207 + t203 * t206;
t141 = t162 * t198;
t201 = sin(qJ(5));
t205 = cos(qJ(5));
t81 = t205 * t139 + t141 * t201;
t190 = -pkin(3) * t207 - pkin(2);
t208 = cos(qJ(2));
t278 = pkin(1) * qJD(1);
t244 = t208 * t278;
t142 = t190 * t198 - t244;
t88 = t139 * pkin(4) + t142;
t286 = t88 * t81;
t204 = sin(qJ(2));
t246 = t204 * t278;
t171 = pkin(7) * t198 + t246;
t232 = pkin(8) * t198 + t171;
t130 = t232 * t203;
t123 = qJD(3) * pkin(3) - t130;
t131 = t232 * t207;
t252 = qJD(4) * t206;
t253 = qJD(4) * t202;
t217 = qJD(3) * t232;
t277 = pkin(1) * qJD(2);
t243 = qJD(1) * t277;
t222 = t208 * t243;
t97 = -t203 * t217 + t207 * t222;
t98 = -t203 * t222 - t207 * t217;
t223 = -t123 * t252 + t131 * t253 - t202 * t98 - t206 * t97;
t197 = qJD(3) + qJD(4);
t111 = t197 * t162;
t86 = t111 * t198;
t19 = -pkin(9) * t86 - t223;
t228 = -t202 * t97 + t206 * t98;
t122 = t206 * t131;
t68 = t123 * t202 + t122;
t34 = -qJD(4) * t68 + t228;
t218 = t197 * t266;
t255 = qJD(3) * t207;
t297 = -t206 * t255 - t207 * t252;
t260 = t297 * t198;
t85 = t198 * t218 + t260;
t20 = t85 * pkin(9) + t34;
t251 = qJD(5) * t201;
t293 = pkin(9) * t139;
t58 = t68 - t293;
t230 = t201 * t20 - t58 * t251;
t136 = t141 * pkin(9);
t120 = t202 * t131;
t67 = t206 * t123 - t120;
t57 = -t136 + t67;
t55 = pkin(4) * t197 + t57;
t4 = (qJD(5) * t55 + t19) * t205 + t230;
t302 = t286 - t4;
t216 = t139 * t201 - t205 * t141;
t287 = t88 * t216;
t275 = t205 * t58;
t22 = t201 * t55 + t275;
t231 = -t201 * t19 + t205 * t20;
t5 = -qJD(5) * t22 + t231;
t301 = t287 + t5;
t296 = -pkin(8) - pkin(7);
t180 = t296 * t203;
t195 = t207 * pkin(8);
t181 = pkin(7) * t207 + t195;
t125 = t202 * t180 + t206 * t181;
t240 = qJD(3) * t296;
t166 = t203 * t240;
t167 = t207 * t240;
t272 = qJD(4) * t125 - t162 * t244 + t202 * t166 - t206 * t167;
t161 = -t263 + t266;
t271 = -t161 * t244 - t206 * t166 - t202 * t167 - t180 * t252 + t181 * t253;
t300 = t81 * t216;
t110 = t218 + t297;
t291 = t110 * pkin(9);
t299 = t291 - t272;
t108 = t111 * pkin(9);
t298 = -t108 - t271;
t199 = t203 ^ 2;
t200 = t207 ^ 2;
t257 = t199 + t200;
t17 = t216 ^ 2 - t81 ^ 2;
t193 = qJD(5) + t197;
t250 = qJD(5) * t205;
t29 = t139 * t250 + t141 * t251 + t201 * t86 + t205 * t85;
t15 = t193 * t81 - t29;
t210 = qJD(5) * t216 + t201 * t85 - t205 * t86;
t16 = -t193 * t216 + t210;
t187 = pkin(1) * t204 + pkin(7);
t285 = -pkin(8) - t187;
t158 = t285 * t203;
t159 = t187 * t207 + t195;
t100 = t202 * t158 + t206 * t159;
t295 = pkin(1) * t208;
t294 = pkin(4) * t141;
t292 = pkin(9) * t162;
t124 = t206 * t180 - t181 * t202;
t89 = t124 - t292;
t157 = t161 * pkin(9);
t90 = -t157 + t125;
t53 = -t201 * t90 + t205 * t89;
t284 = qJD(5) * t53 + t299 * t201 + t298 * t205;
t54 = t201 * t89 + t205 * t90;
t283 = -qJD(5) * t54 - t298 * t201 + t299 * t205;
t188 = pkin(3) * t206 + pkin(4);
t265 = t202 * t205;
t69 = t130 * t202 - t122;
t59 = t69 + t293;
t70 = -t206 * t130 - t120;
t60 = -t136 + t70;
t282 = -t201 * t60 + t205 * t59 + t188 * t251 - (-t202 * t250 + (-t201 * t206 - t265) * qJD(4)) * pkin(3);
t267 = t201 * t202;
t281 = t201 * t59 + t205 * t60 - t188 * t250 - (-t202 * t251 + (t205 * t206 - t267) * qJD(4)) * pkin(3);
t106 = t205 * t161 + t162 * t201;
t107 = -t161 * t201 + t162 * t205;
t44 = qJD(5) * t107 - t201 * t110 + t205 * t111;
t185 = t204 * t243;
t256 = qJD(3) * t203;
t239 = t198 * t256;
t148 = pkin(3) * t239 + t185;
t63 = pkin(4) * t86 + t148;
t280 = t63 * t106 + t88 * t44;
t43 = t205 * t110 + t201 * t111 + t161 * t250 + t162 * t251;
t279 = t63 * t107 - t88 * t43;
t276 = t201 * t58;
t274 = t142 * t111 + t148 * t161;
t273 = -t142 * t110 + t148 * t162;
t270 = t141 * t139;
t269 = t142 * t141;
t268 = t198 * t203;
t264 = t204 * t207;
t209 = qJD(3) ^ 2;
t262 = t209 * t203;
t194 = t209 * t207;
t172 = -pkin(2) * t198 - t244;
t261 = t172 * t255 + t203 * t185;
t259 = t257 * t222;
t258 = t199 - t200;
t254 = qJD(3) * t208;
t249 = -qJD(1) - t198;
t248 = -qJD(2) + t198;
t247 = pkin(3) * t268;
t245 = t208 * t277;
t192 = t204 * t277;
t191 = pkin(3) * t256;
t196 = t198 ^ 2;
t241 = t203 * t196 * t207;
t238 = t198 * t255;
t237 = t203 * t254;
t21 = t205 * t55 - t276;
t234 = -t4 * t106 - t5 * t107 + t21 * t43 - t22 * t44;
t233 = -pkin(4) * t193 - t55;
t94 = pkin(4) * t111 + t191;
t227 = qJD(3) * t285;
t226 = t67 * t110 - t68 * t111 + t161 * t223 - t34 * t162;
t99 = t206 * t158 - t159 * t202;
t224 = t257 * qJD(2);
t221 = t203 * t238;
t220 = t94 - t246;
t219 = -t21 * t81 - t216 * t22;
t75 = t99 - t292;
t76 = -t157 + t100;
t41 = -t201 * t76 + t205 * t75;
t42 = t201 * t75 + t205 * t76;
t137 = pkin(4) * t161 + t190;
t215 = t142 * t139 + t223;
t126 = t203 * t227 + t207 * t245;
t127 = -t203 * t245 + t207 * t227;
t47 = t206 * t126 + t202 * t127 + t158 * t252 - t159 * t253;
t214 = t191 - t246;
t213 = -t172 * t198 - t222;
t212 = -t204 * t268 + t207 * t254;
t48 = -qJD(4) * t100 - t202 * t126 + t206 * t127;
t189 = -pkin(2) - t295;
t175 = t190 - t295;
t170 = -0.2e1 * t221;
t169 = 0.2e1 * t221;
t168 = t192 + t191;
t151 = t172 * t256;
t150 = pkin(3) * t265 + t188 * t201;
t149 = -pkin(3) * t267 + t188 * t205;
t138 = -0.2e1 * t258 * t198 * qJD(3);
t129 = t137 - t295;
t114 = t247 + t294;
t102 = t111 * t197;
t101 = t110 * t197;
t87 = t192 + t94;
t64 = -t139 ^ 2 + t141 ^ 2;
t61 = -t260 + (t139 - t242) * t197;
t40 = t44 * t193;
t39 = t43 * t193;
t38 = t111 * t139 + t161 * t86;
t37 = -t110 * t141 - t162 * t85;
t36 = t48 + t291;
t35 = -t108 + t47;
t24 = t205 * t57 - t276;
t23 = -t201 * t57 - t275;
t14 = t110 * t139 - t111 * t141 + t161 * t85 - t162 * t86;
t9 = -t106 * t210 + t44 * t81;
t8 = -t107 * t29 + t216 * t43;
t7 = -qJD(5) * t42 - t201 * t35 + t205 * t36;
t6 = qJD(5) * t41 + t201 * t36 + t205 * t35;
t1 = t106 * t29 + t107 * t210 + t216 * t44 + t43 * t81;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192 * t198 - t185, t249 * t245, 0, 0, t169, t138, t194, t170, -t262, 0, t189 * t239 - t187 * t194 + t151 + (t249 * t264 - t237) * t277, t187 * t262 + t189 * t238 - t212 * t277 + t261, t198 * t224 * t295 + t259, ((qJD(1) * t189 + t172) * t204 + (qJD(1) * t187 + t171) * t208 * t257) * t277, t37, t14, -t101, t38, -t102, 0, t139 * t168 + t175 * t86 + t197 * t48 + t274, t141 * t168 - t175 * t85 - t197 * t47 + t273, -t100 * t86 - t139 * t47 - t141 * t48 + t85 * t99 + t226, -t100 * t223 + t142 * t168 + t148 * t175 + t34 * t99 + t47 * t68 + t48 * t67, t8, t1, -t39, t9, -t40, 0, -t129 * t210 + t193 * t7 + t81 * t87 + t280, -t129 * t29 - t193 * t6 - t216 * t87 + t279, t210 * t42 + t216 * t7 + t29 * t41 - t6 * t81 + t234, t129 * t63 + t21 * t7 + t22 * t6 + t4 * t42 + t41 * t5 + t87 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198 * t246 - t185, t248 * t244, 0, 0, t169, t138, t194, t170, -t262, 0, -pkin(2) * t239 - pkin(7) * t194 + t151 + (t248 * t264 + t237) * t278, -pkin(2) * t238 + pkin(7) * t262 + t212 * t278 + t261, -t198 * t244 * t257 + t259, ((-pkin(2) * qJD(2) - t172) * t204 + (pkin(7) * t224 - t171 * t257) * t208) * t278, t37, t14, -t101, t38, -t102, 0, t214 * t139 + t190 * t86 - t272 * t197 + t274, t214 * t141 - t190 * t85 + t271 * t197 + t273, t124 * t85 - t125 * t86 + t271 * t139 + t272 * t141 + t226, t34 * t124 - t125 * t223 + t214 * t142 + t148 * t190 - t271 * t68 - t272 * t67, t8, t1, -t39, t9, -t40, 0, -t137 * t210 + t193 * t283 + t220 * t81 + t280, -t137 * t29 - t193 * t284 - t216 * t220 + t279, t210 * t54 + t216 * t283 - t284 * t81 + t53 * t29 + t234, t63 * t137 + t21 * t283 + t22 * t284 + t220 * t88 + t4 * t54 + t5 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, t258 * t196, 0, t241, 0, 0, t213 * t203, t213 * t207, 0, 0, t270, t64, t61, -t270, 0, 0, -t139 * t247 - t269 - t69 * t197 + (-t122 + (-pkin(3) * t197 - t123) * t202) * qJD(4) + t228, t70 * t197 + (-t141 * t268 - t197 * t252) * pkin(3) + t215, (t68 + t69) * t141 + (-t67 + t70) * t139 + (-t202 * t86 + t206 * t85 + (-t139 * t206 + t141 * t202) * qJD(4)) * pkin(3), -t67 * t69 - t68 * t70 + (-t142 * t268 - t202 * t223 + t206 * t34 + (-t202 * t67 + t206 * t68) * qJD(4)) * pkin(3), -t300, t17, t15, t300, t16, 0, -t114 * t81 - t193 * t282 + t301, t114 * t216 + t193 * t281 + t302, t149 * t29 + t150 * t210 - t216 * t282 + t281 * t81 + t219, -t88 * t114 + t5 * t149 + t4 * t150 - t21 * t282 - t22 * t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t270, t64, t61, -t270, 0, 0, t68 * t197 - t269 + t34, t197 * t67 + t215, 0, 0, -t300, t17, t15, t300, t16, 0, -t81 * t294 - t23 * t193 + t287 + (t201 * t233 - t275) * qJD(5) + t231, t216 * t294 + t24 * t193 + t286 + (qJD(5) * t233 - t19) * t205 - t230, -t23 * t216 + t24 * t81 + (t201 * t210 + t205 * t29 + (-t201 * t216 - t205 * t81) * qJD(5)) * pkin(4) + t219, -t21 * t23 - t22 * t24 + (-t141 * t88 + t201 * t4 + t205 * t5 + (-t201 * t21 + t205 * t22) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t300, t17, t15, t300, t16, 0, t22 * t193 + t301, t21 * t193 + t302, 0, 0;];
tauc_reg = t2;
