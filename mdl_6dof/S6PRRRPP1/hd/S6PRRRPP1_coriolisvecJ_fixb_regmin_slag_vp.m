% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:46:54
% EndTime: 2019-03-08 22:47:03
% DurationCPUTime: 3.65s
% Computational Cost: add. (4779->384), mult. (12032->531), div. (0->0), fcn. (8707->10), ass. (0->200)
t175 = sin(qJ(3));
t178 = cos(qJ(3));
t194 = pkin(3) * t175 - pkin(9) * t178;
t139 = t194 * qJD(3);
t145 = -pkin(3) * t178 - pkin(9) * t175 - pkin(2);
t174 = sin(qJ(4));
t176 = sin(qJ(2));
t177 = cos(qJ(4));
t228 = qJD(4) * t177;
t171 = sin(pkin(6));
t238 = qJD(1) * t171;
t179 = cos(qJ(2));
t249 = t178 * t179;
t290 = -(t174 * t176 + t177 * t249) * t238 + t174 * t139 + t145 * t228;
t234 = qJD(2) * t178;
t214 = t174 * t234;
t230 = qJD(4) * t174;
t289 = -t214 + t230;
t232 = qJD(3) * t175;
t280 = pkin(8) * t174;
t288 = t177 * t139 + t232 * t280 - (-t174 * t249 + t176 * t177) * t238;
t250 = t177 * t178;
t160 = pkin(8) * t250;
t189 = pkin(4) * t175 - qJ(5) * t250;
t226 = t177 * qJD(5);
t287 = -t175 * t226 + t189 * qJD(3) + (-t160 + (qJ(5) * t175 - t145) * t174) * qJD(4) + t288;
t251 = t175 * t177;
t286 = -(-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t251 - (-qJD(5) * t175 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t178) * t174 - t290;
t157 = -qJD(4) + t234;
t235 = qJD(2) * t175;
t205 = t174 * t235;
t227 = t177 * qJD(3);
t135 = t205 - t227;
t233 = qJD(3) * t174;
t137 = t177 * t235 + t233;
t170 = sin(pkin(11));
t172 = cos(pkin(11));
t79 = t172 * t135 + t137 * t170;
t285 = t157 * t79;
t131 = t170 * t177 + t172 * t174;
t118 = t131 * qJD(4);
t245 = t131 * t234 - t118;
t256 = t172 * t177;
t244 = t289 * t170 - t172 * t228 + t234 * t256;
t231 = qJD(3) * t178;
t213 = t174 * t231;
t284 = t175 * t228 + t213;
t191 = -t135 * t170 + t172 * t137;
t283 = t191 ^ 2;
t275 = t286 * t170 + t287 * t172;
t274 = t287 * t170 - t286 * t172;
t278 = -qJ(5) - pkin(9);
t200 = qJD(4) * t278;
t115 = t174 * t200 + t226;
t184 = -t174 * qJD(5) + t177 * t200;
t207 = t176 * t238;
t141 = qJD(2) * pkin(8) + t207;
t173 = cos(pkin(6));
t237 = qJD(1) * t178;
t107 = -t175 * t141 + t173 * t237;
t138 = t194 * qJD(2);
t199 = -t174 * t107 + t177 * t138;
t60 = qJD(2) * t189 + t199;
t264 = t177 * t107 + t174 * t138;
t67 = -qJ(5) * t214 + t264;
t270 = (t184 - t60) * t172 + (-t115 + t67) * t170;
t99 = -qJD(3) * pkin(3) - t107;
t72 = pkin(4) * t135 + qJD(5) + t99;
t31 = pkin(5) * t79 - qJ(6) * t191 + t72;
t282 = t31 * t191;
t255 = t173 * t175;
t155 = qJD(1) * t255;
t108 = t178 * t141 + t155;
t281 = t289 * pkin(4) - t108;
t224 = qJD(3) * qJD(4);
t103 = qJD(2) * t284 + t174 * t224;
t100 = qJD(3) * pkin(9) + t108;
t106 = (t139 + t207) * qJD(2);
t206 = t179 * t238;
t109 = qJD(2) * t145 - t206;
t236 = qJD(2) * t171;
t215 = t179 * t236;
t74 = -t141 * t232 + (qJD(3) * t173 + t215) * t237;
t186 = -t100 * t230 + t174 * t106 + t109 * t228 + t177 * t74;
t12 = -qJ(5) * t103 - qJD(5) * t135 + t186;
t212 = t178 * t227;
t102 = qJD(2) * t212 - qJD(4) * t205 + t177 * t224;
t201 = -t177 * t106 + t174 * t74;
t254 = t174 * t109;
t64 = t100 * t177 + t254;
t182 = -qJD(4) * t64 - t201;
t225 = qJD(2) * qJD(3);
t203 = t175 * t225;
t9 = pkin(4) * t203 - t102 * qJ(5) - t137 * qJD(5) + t182;
t3 = -t170 * t12 + t172 * t9;
t4 = t172 * t12 + t170 * t9;
t50 = -qJ(5) * t135 + t64;
t44 = t172 * t50;
t62 = -t100 * t174 + t177 * t109;
t49 = -qJ(5) * t137 + t62;
t23 = t170 * t49 + t44;
t279 = t23 * t191;
t277 = qJ(6) * t232 - qJD(6) * t178 + t274;
t276 = -pkin(5) * t232 - t275;
t30 = t170 * t60 + t172 * t67;
t27 = qJ(6) * t235 + t30;
t71 = t172 * t115 + t170 * t184;
t273 = t27 - t71;
t272 = pkin(5) * t235 - t270;
t271 = pkin(5) * t245 - qJ(6) * t244 + qJD(6) * t131 - t281;
t36 = -pkin(4) * t157 + t49;
t18 = t170 * t36 + t44;
t133 = t177 * t145;
t84 = -qJ(5) * t251 + t133 + (-pkin(4) - t280) * t178;
t241 = t174 * t145 + t160;
t253 = t174 * t175;
t90 = -qJ(5) * t253 + t241;
t52 = t170 * t84 + t172 * t90;
t269 = qJD(2) * pkin(2);
t268 = t170 * t50;
t267 = t174 * t99;
t196 = t175 * t215;
t75 = qJD(1) * t196 + qJD(3) * t155 + t141 * t231;
t266 = t75 * t174;
t265 = t75 * t177;
t263 = t102 * t174;
t262 = t135 * t157;
t261 = t137 * t157;
t260 = t157 * t177;
t259 = t171 * t176;
t258 = t171 * t179;
t181 = qJD(2) ^ 2;
t257 = t171 * t181;
t252 = t174 * t178;
t180 = qJD(3) ^ 2;
t248 = t180 * t175;
t247 = t180 * t178;
t24 = t172 * t49 - t268;
t246 = qJD(6) - t24;
t240 = pkin(4) * t253 + t175 * pkin(8);
t168 = t175 ^ 2;
t239 = -t178 ^ 2 + t168;
t229 = qJD(4) * t175;
t223 = qJ(6) * t203 + t4;
t221 = t176 * t257;
t219 = pkin(4) * t284 + pkin(8) * t231;
t218 = -pkin(4) * t177 - pkin(3);
t216 = t176 * t236;
t211 = t157 * t230;
t210 = t174 * t229;
t208 = t157 * t228;
t204 = t278 * t174;
t58 = t102 * t170 + t172 * t103;
t59 = t102 * t172 - t103 * t170;
t146 = t278 * t177;
t96 = -t170 * t146 - t172 * t204;
t97 = -t172 * t146 + t170 * t204;
t198 = -t97 * t58 + t96 * t59 - t71 * t79;
t197 = t175 * t206;
t195 = t178 * t215;
t193 = -t79 ^ 2 - t283;
t142 = -t206 - t269;
t192 = -t142 - t206;
t17 = t172 * t36 - t268;
t51 = -t170 * t90 + t172 * t84;
t130 = t170 * t174 - t256;
t190 = qJD(2) * t168 - t157 * t178;
t55 = pkin(4) * t103 + t75;
t122 = t178 * t259 + t255;
t88 = -t122 * t174 - t177 * t258;
t188 = -t122 * t177 + t174 * t258;
t121 = -t173 * t178 + t175 * t259;
t2 = -pkin(5) * t203 - t3;
t86 = -qJD(3) * t121 + t195;
t42 = qJD(4) * t188 - t86 * t174 + t177 * t216;
t43 = qJD(4) * t88 + t174 * t216 + t86 * t177;
t20 = t170 * t43 - t172 * t42;
t22 = t170 * t42 + t172 * t43;
t53 = -t170 * t188 - t172 * t88;
t54 = t170 * t88 - t172 * t188;
t187 = t191 * t20 - t22 * t79 + t53 * t59 - t54 * t58;
t183 = qJD(3) * (-t192 - t269);
t5 = pkin(5) * t58 - qJ(6) * t59 - qJD(6) * t191 + t55;
t163 = -pkin(4) * t172 - pkin(5);
t161 = pkin(4) * t170 + qJ(6);
t114 = -t170 * t253 + t172 * t251;
t113 = t131 * t175;
t87 = qJD(3) * t122 + t196;
t76 = pkin(5) * t130 - qJ(6) * t131 + t218;
t69 = t118 * t175 + t170 * t213 - t172 * t212;
t68 = t130 * t229 - t131 * t231;
t66 = pkin(5) * t113 - qJ(6) * t114 + t240;
t47 = pkin(5) * t178 - t51;
t46 = -qJ(6) * t178 + t52;
t32 = pkin(4) * t137 + pkin(5) * t191 + qJ(6) * t79;
t25 = -pkin(5) * t68 + qJ(6) * t69 - qJD(6) * t114 + t219;
t15 = -qJ(6) * t157 + t18;
t13 = pkin(5) * t157 + qJD(6) - t17;
t1 = -qJD(6) * t157 + t223;
t6 = [0, 0, -t221, -t179 * t257, 0, 0, 0, 0, 0, -t178 * t221 + (-t87 - t196) * qJD(3), t175 * t221 + (-t86 - t195) * qJD(3), 0, 0, 0, 0, 0, t103 * t121 + t135 * t87 - t157 * t42 + t203 * t88, t102 * t121 + t137 * t87 + t157 * t43 + t188 * t203, t187, t121 * t55 - t17 * t20 + t18 * t22 - t3 * t53 + t4 * t54 + t72 * t87, t121 * t58 + t157 * t20 - t203 * t53 + t79 * t87, t187, -t121 * t59 - t157 * t22 - t191 * t87 + t203 * t54, t1 * t54 + t121 * t5 + t13 * t20 + t15 * t22 + t2 * t53 + t31 * t87; 0, 0, 0, 0, 0.2e1 * t178 * t203, -0.2e1 * t239 * t225, t247, -t248, 0, -pkin(8) * t247 + t175 * t183, pkin(8) * t248 + t178 * t183, t102 * t251 + (-t210 + t212) * t137 (-t135 * t177 - t137 * t174) * t231 + (-t263 - t103 * t177 + (t135 * t174 - t137 * t177) * qJD(4)) * t175, t157 * t210 - t102 * t178 + (t137 * t175 + t177 * t190) * qJD(3), t175 * t208 + t103 * t178 + (-t135 * t175 - t174 * t190) * qJD(3) (-t157 - t234) * t232 (t145 * t230 - t288) * t157 + ((pkin(8) * t135 + t267) * qJD(3) + (t254 + (pkin(8) * t157 + t100) * t177) * qJD(4) + t201) * t178 + (-t135 * t206 + t99 * t228 + pkin(8) * t103 + t266 + ((-pkin(8) * t252 + t133) * qJD(2) + t62) * qJD(3)) * t175, t290 * t157 + (t99 * t227 + (qJD(3) * t137 - t211) * pkin(8) + t186) * t178 + (-t137 * t206 - t99 * t230 + pkin(8) * t102 + t265 + (-pkin(8) * t260 - t241 * qJD(2) - t64) * qJD(3)) * t175, -t4 * t113 - t3 * t114 + t17 * t69 + t18 * t68 - t191 * t275 - t274 * t79 - t51 * t59 - t52 * t58, t4 * t52 + t3 * t51 + t55 * t240 + (-t197 + t219) * t72 + t274 * t18 + t275 * t17, t5 * t113 + t2 * t178 + t25 * t79 - t31 * t68 + t66 * t58 + t276 * t157 + (-t79 * t206 + (-qJD(2) * t47 - t13) * qJD(3)) * t175, -t1 * t113 + t2 * t114 - t13 * t69 + t15 * t68 + t191 * t276 - t277 * t79 - t46 * t58 + t47 * t59, -t1 * t178 - t5 * t114 - t25 * t191 + t31 * t69 - t66 * t59 - t277 * t157 + (t191 * t206 + (qJD(2) * t46 + t15) * qJD(3)) * t175, t1 * t46 + t2 * t47 + t5 * t66 + (t25 - t197) * t31 + t277 * t15 + t276 * t13; 0, 0, 0, 0, -t175 * t181 * t178, t239 * t181, 0, 0, 0, qJD(3) * t108 - t142 * t235 - t75, t192 * t234, -t137 * t260 + t263 (t102 + t262) * t177 + (-t103 + t261) * t174, -t208 + (t157 * t250 + (-t137 + t233) * t175) * qJD(2), t211 + (-t157 * t252 + (t135 + t227) * t175) * qJD(2), t157 * t235, -pkin(3) * t103 - t265 + t199 * t157 - t108 * t135 + (pkin(9) * t260 + t267) * qJD(4) + (-t62 * t175 + (-pkin(9) * t232 - t178 * t99) * t174) * qJD(2), -pkin(3) * t102 + t266 - t264 * t157 - t108 * t137 + (-pkin(9) * t157 * t174 + t177 * t99) * qJD(4) + (-t99 * t250 + (-pkin(9) * t227 + t64) * t175) * qJD(2), -t4 * t130 - t3 * t131 + t244 * t17 + t245 * t18 - t191 * t270 + t30 * t79 + t198, t4 * t97 - t3 * t96 + t55 * t218 + t281 * t72 + (t71 - t30) * t18 + t270 * t17, t5 * t130 + t76 * t58 - t271 * t79 - t245 * t31 + t272 * t157 + (-qJD(3) * t96 + t13) * t235, -t1 * t130 - t13 * t244 + t2 * t131 + t15 * t245 + t191 * t272 + t27 * t79 + t198, -t5 * t131 - t76 * t59 + t271 * t191 + t244 * t31 + t273 * t157 + (qJD(3) * t97 - t15) * t235, t1 * t97 + t13 * t272 - t15 * t273 + t2 * t96 - t271 * t31 + t5 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137 * t135, -t135 ^ 2 + t137 ^ 2, t102 - t262, -t103 - t261, t203, -t99 * t137 - t64 * t157 + t182, t135 * t99 - t157 * t62 - t186, t18 * t191 - t279 + (-t170 * t58 - t172 * t59) * pkin(4) + (-t17 + t24) * t79, t17 * t23 - t18 * t24 + (-t137 * t72 + t170 * t4 + t172 * t3) * pkin(4), -t23 * t157 - t282 - t32 * t79 + (pkin(5) - t163) * t203 + t3, t15 * t191 - t161 * t58 + t163 * t59 - t279 + (t13 - t246) * t79, t161 * t203 - t31 * t79 + t32 * t191 + (-0.2e1 * qJD(6) + t24) * t157 + t223, t1 * t161 - t13 * t23 + t15 * t246 + t2 * t163 - t31 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t17 * t191 + t18 * t79 + t55, -t157 * t191 + t58, t193, -t59 - t285, -t13 * t191 + t15 * t79 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191 * t79 - t203, t59 - t285, -t157 ^ 2 - t283, t15 * t157 + t2 + t282;];
tauc_reg  = t6;
