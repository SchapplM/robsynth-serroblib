% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPPR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:33
% EndTime: 2019-03-09 02:57:40
% DurationCPUTime: 4.02s
% Computational Cost: add. (5437->469), mult. (10669->539), div. (0->0), fcn. (7183->10), ass. (0->244)
t157 = sin(qJ(3));
t273 = sin(pkin(9));
t222 = t273 * t157;
t116 = qJD(1) * t222;
t154 = cos(pkin(9));
t160 = cos(qJ(3));
t98 = t154 * t157 + t160 * t273;
t169 = -qJD(3) * t116 + qJDD(1) * t98;
t246 = qJD(1) * t160;
t226 = t154 * t246;
t92 = -t116 + t226;
t328 = t169 + qJD(3) * (-t92 + t226);
t327 = t169 + qJD(3) * (t92 + t226);
t94 = t98 * qJD(3);
t99 = t154 * t160 - t222;
t281 = -qJD(3) * t94 + qJDD(3) * t99;
t89 = t98 * qJD(1);
t326 = qJD(1) * t89 - t281;
t159 = cos(qJ(6));
t81 = qJD(6) + t92;
t218 = t159 * t81;
t156 = sin(qJ(6));
t64 = qJD(3) * t159 + t156 * t89;
t325 = t64 * t218;
t146 = qJ(3) + pkin(9);
t136 = sin(t146);
t158 = sin(qJ(1));
t161 = cos(qJ(1));
t205 = g(1) * t161 + g(2) * t158;
t182 = t205 * t136;
t245 = qJD(3) * t157;
t162 = -pkin(1) - pkin(7);
t258 = qJ(4) - t162;
t74 = -qJD(4) * t160 + t245 * t258;
t104 = t258 * t160;
t75 = -qJD(3) * t104 - t157 * qJD(4);
t45 = -t154 * t74 + t273 * t75;
t103 = t258 * t157;
t60 = -t103 * t273 + t104 * t154;
t324 = -t45 * qJD(3) - t60 * qJDD(3) - t182;
t247 = qJD(1) * t157;
t102 = pkin(3) * t247 + qJD(1) * qJ(2) + qJD(4);
t144 = g(2) * t161;
t315 = g(1) * t158 - t144;
t323 = -qJD(1) * t102 - t315;
t164 = qJD(1) ^ 2;
t322 = -t164 * qJ(2) - t315;
t149 = qJDD(1) * qJ(2);
t113 = qJD(1) * t162 + qJD(2);
t79 = -qJ(4) * t247 + t113 * t157;
t68 = t273 * t79;
t80 = -qJ(4) * t246 + t113 * t160;
t71 = qJD(3) * pkin(3) + t80;
t42 = t154 * t71 - t68;
t207 = qJD(5) - t42;
t301 = t92 * pkin(5);
t305 = pkin(4) + pkin(8);
t22 = -qJD(3) * t305 + t207 + t301;
t185 = -t92 * qJ(5) + t102;
t26 = t305 * t89 + t185;
t10 = t156 * t22 + t159 * t26;
t236 = t160 * qJDD(1);
t181 = -qJDD(1) * t222 + t154 * t236;
t59 = qJD(1) * t94 - t181;
t274 = t59 * qJ(5);
t150 = qJD(1) * qJD(2);
t240 = qJD(1) * qJD(3);
t225 = t160 * t240;
t237 = t157 * qJDD(1);
t67 = qJDD(4) + t149 + t150 + (t225 + t237) * pkin(3);
t168 = -t92 * qJD(5) + t274 + t67;
t58 = -t154 * t225 - t169;
t5 = -t305 * t58 + t168;
t112 = qJDD(1) * t162 + qJDD(2);
t100 = t160 * t112;
t239 = qJD(1) * qJD(4);
t41 = -t160 * t239 - t113 * t245 + qJDD(3) * pkin(3) + t100 + (t157 * t240 - t236) * qJ(4);
t244 = qJD(3) * t160;
t47 = (-qJ(4) * qJD(1) + t113) * t244 + (-qJ(4) * qJDD(1) + t112 - t239) * t157;
t282 = -t154 * t41 + t273 * t47;
t230 = qJDD(5) + t282;
t7 = -t59 * pkin(5) - qJDD(3) * t305 + t230;
t2 = -qJD(6) * t10 - t156 * t5 + t159 * t7;
t321 = t10 * t81 + t2;
t62 = qJD(3) * t156 - t159 * t89;
t223 = t62 * t81;
t242 = qJD(6) * t159;
t243 = qJD(6) * t156;
t23 = qJD(3) * t243 - qJDD(3) * t159 + t156 * t58 - t242 * t89;
t320 = t23 - t223;
t24 = qJD(6) * t64 + t156 * qJDD(3) + t159 * t58;
t319 = t64 * t81 - t24;
t307 = t89 ^ 2;
t88 = t92 ^ 2;
t318 = -t307 - t88;
t317 = -t307 + t88;
t152 = t157 ^ 2;
t153 = t160 ^ 2;
t250 = t152 + t153;
t217 = t250 * t112;
t51 = t154 * t80 - t68;
t257 = -qJD(5) + t51;
t142 = t157 * pkin(3);
t155 = -qJ(4) - pkin(7);
t316 = t142 * t161 + t155 * t158;
t234 = 0.2e1 * t150;
t312 = 0.2e1 * t149 + t234 - t205;
t127 = -pkin(3) * t154 - pkin(4);
t118 = -pkin(8) + t127;
t303 = t89 * pkin(5);
t278 = t154 * t79;
t43 = t273 * t71 + t278;
t38 = -qJD(3) * qJ(5) - t43;
t25 = -t38 - t303;
t56 = -qJDD(6) + t59;
t311 = -t118 * t56 + t25 * t81;
t137 = cos(t146);
t174 = -g(3) * t137 - t136 * t315;
t46 = t154 * t75 + t273 * t74;
t61 = -t103 * t154 - t104 * t273;
t310 = -t46 * qJD(3) - t61 * qJDD(3) - t137 * t205;
t19 = t154 * t47 + t273 * t41;
t91 = qJD(3) * t222 - t154 * t244;
t309 = -t19 * t98 + t282 * t99 + t42 * t94 + t43 * t91;
t235 = qJDD(3) * qJ(5) + t19;
t14 = -qJD(3) * qJD(5) - t235;
t268 = qJDD(3) * pkin(4);
t15 = t230 - t268;
t36 = -qJD(3) * pkin(4) + t207;
t308 = t14 * t98 + t15 * t99 - t36 * t94 - t38 * t91;
t304 = t58 * pkin(4);
t9 = -t156 * t26 + t159 * t22;
t302 = t9 * t81;
t300 = pkin(3) * t160;
t299 = pkin(4) * t136;
t125 = g(3) * t136;
t297 = g(3) * t157;
t287 = t64 * t62;
t285 = t64 * t89;
t284 = t89 * t62;
t283 = t89 * t92;
t277 = t156 * t56;
t276 = t159 * t23;
t53 = t159 * t56;
t275 = t24 * t156;
t272 = pkin(1) * qJDD(1);
t270 = qJD(3) * t89;
t267 = t136 * t158;
t123 = t137 * qJ(5);
t266 = t137 * t158;
t265 = t137 * t161;
t264 = t156 * t158;
t263 = t156 * t161;
t262 = t158 * t159;
t261 = t159 * t161;
t50 = t273 * t80 + t278;
t259 = t50 * qJD(3);
t128 = qJ(2) + t142;
t256 = t301 - t257;
t255 = t123 - t142;
t254 = (t234 + t149) * qJ(2);
t253 = pkin(1) * t161 + qJ(2) * t158;
t251 = t152 - t153;
t163 = qJD(3) ^ 2;
t249 = -t163 - t164;
t115 = pkin(3) * t244 + qJD(2);
t238 = qJDD(3) * t157;
t233 = t98 * t243;
t232 = t98 * t242;
t231 = t160 * t164 * t157;
t229 = pkin(4) * t266 + qJ(5) * t267 + t158 * t300;
t228 = g(1) * t266 - g(2) * t265 - t125;
t227 = t142 * t158 + t253;
t141 = t161 * qJ(2);
t224 = -t158 * pkin(1) + t141;
t221 = pkin(3) * t246 + qJ(5) * t89;
t220 = t156 * t81;
t216 = qJD(6) * t99 + qJD(1);
t214 = t250 * qJDD(1);
t213 = qJDD(2) - t272;
t212 = t157 * t225;
t211 = -t2 * t99 + t9 * t94;
t1 = qJD(6) * t9 + t156 * t7 + t159 * t5;
t210 = -t1 * t99 + t10 * t94;
t8 = pkin(5) * t58 - t14;
t209 = t25 * t91 - t8 * t98;
t208 = -t99 * qJ(5) + t128;
t206 = t1 - t302;
t201 = -t98 * t23 - t91 * t64;
t200 = -t23 * t99 - t64 * t94;
t199 = -t98 * t24 + t91 * t62;
t198 = t24 * t99 - t62 * t94;
t197 = -t56 * t98 - t81 * t91;
t196 = -t58 * t98 - t89 * t91;
t195 = -t59 * t99 - t92 * t94;
t193 = -t228 - t282;
t192 = t10 * t159 - t156 * t9;
t191 = -qJ(5) * t136 - t300;
t190 = pkin(8) * t136 - t123;
t37 = t305 * t98 + t208;
t48 = pkin(5) * t99 + t60;
t17 = t156 * t48 + t159 * t37;
t16 = -t156 * t37 + t159 * t48;
t187 = qJD(3) * t91 - qJDD(3) * t98;
t186 = t224 + t316;
t184 = -t161 * t155 + t227;
t183 = -t220 * t81 - t53;
t180 = 0.2e1 * qJ(2) * t240 + qJDD(3) * t162;
t178 = t94 * qJ(5) - t99 * qJD(5) + t115;
t177 = -t218 * t81 + t277;
t176 = qJD(1) * t92 - t187;
t44 = pkin(4) * t89 + t185;
t175 = t44 * t92 + qJDD(5) - t193;
t173 = t58 * t99 + t59 * t98 + t89 * t94 + t91 * t92;
t172 = -t195 - t196;
t171 = -t205 + t67;
t167 = t45 * t92 - t46 * t89 + t58 * t61 - t59 * t60 + t315;
t166 = -t162 * t163 + t312;
t165 = -qJD(6) * t118 * t81 + t174 + t8;
t139 = qJDD(3) * t160;
t124 = pkin(3) * t273 + qJ(5);
t109 = t161 * t299;
t107 = pkin(4) * t267;
t85 = -t137 * t264 + t261;
t84 = -t137 * t262 - t263;
t83 = -t137 * t263 - t262;
t82 = -t137 * t261 + t264;
t57 = pkin(4) * t98 + t208;
t52 = pkin(4) * t92 + t221;
t49 = -pkin(5) * t98 + t61;
t33 = -t270 + t59;
t32 = -pkin(4) * t91 + t178;
t31 = t305 * t92 + t221;
t29 = t50 - t303;
t28 = pkin(5) * t91 + t46;
t27 = -pkin(5) * t94 + t45;
t21 = t159 * t24;
t20 = -t305 * t91 + t178;
t13 = t168 - t304;
t12 = t156 * t29 + t159 * t31;
t11 = -t156 * t31 + t159 * t29;
t4 = -qJD(6) * t17 - t156 * t20 + t159 * t27;
t3 = qJD(6) * t16 + t156 * t27 + t159 * t20;
t6 = [0, 0, 0, 0, 0, qJDD(1), t315, t205, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t315 - 0.2e1 * t272, t312, -pkin(1) * t213 - g(1) * t224 - g(2) * t253 + t254, qJDD(1) * t153 - 0.2e1 * t212, -0.2e1 * t157 * t236 + 0.2e1 * t240 * t251, -t157 * t163 + t139, qJDD(1) * t152 + 0.2e1 * t212, -t160 * t163 - t238, 0, t157 * t166 + t160 * t180, -t157 * t180 + t160 * t166, -t162 * t214 - t217 + t315, -g(1) * (t158 * t162 + t141) - g(2) * (pkin(7) * t161 + t253) + t162 * t217 + t254, t195, t173, t281, t196, t187, 0, -t102 * t91 + t115 * t89 - t128 * t58 + t67 * t98 + t324, -t102 * t94 + t115 * t92 - t128 * t59 + t67 * t99 + t310, t167 + t309, -g(1) * t186 - g(2) * t184 + t102 * t115 + t128 * t67 + t19 * t61 + t282 * t60 - t42 * t45 + t43 * t46, 0, -t281, -t187, t195, t173, t196, t167 + t308, -t13 * t98 - t32 * t89 + t44 * t91 + t57 * t58 - t324, -t13 * t99 - t32 * t92 + t44 * t94 + t57 * t59 - t310, t13 * t57 + t44 * t32 - t14 * t61 - t38 * t46 + t15 * t60 + t36 * t45 - g(1) * (-qJ(5) * t265 + t109 + t186) - g(2) * (-qJ(5) * t266 + t107 + t184) t156 * t201 + t232 * t64 (t156 * t62 - t159 * t64) * t91 + (-t275 - t276 + (-t156 * t64 - t159 * t62) * qJD(6)) * t98, t156 * t197 + t232 * t81 + t200, t159 * t199 + t233 * t62, t159 * t197 - t233 * t81 - t198, -t56 * t99 - t81 * t94, -g(1) * t83 - g(2) * t85 + t159 * t209 - t16 * t56 + t233 * t25 + t49 * t24 + t28 * t62 + t4 * t81 - t211, -g(1) * t82 - g(2) * t84 - t156 * t209 + t17 * t56 - t49 * t23 + t232 * t25 + t28 * t64 - t3 * t81 + t210, t16 * t23 - t17 * t24 - t3 * t62 - t4 * t64 - t192 * t91 - t182 + (t1 * t159 - t2 * t156 + (-t10 * t156 - t159 * t9) * qJD(6)) * t98, t1 * t17 + t10 * t3 + t2 * t16 + t9 * t4 + t8 * t49 + t25 * t28 - g(1) * (t109 + t141 + t316) - g(2) * (t107 + t227) + (-g(1) * t190 - g(2) * (pkin(5) - t155)) * t161 + (-g(1) * (-pkin(1) - pkin(5)) - g(2) * t190) * t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t164, t322 + t213, 0, 0, 0, 0, 0, 0, t157 * t249 + t139, t160 * t249 - t238, -t214, t217 + t322, 0, 0, 0, 0, 0, 0, -t326, -t176, t172, -t309 + t323, 0, 0, 0, 0, 0, 0, t172, t326, t176, -qJD(1) * t44 - t308 - t315, 0, 0, 0, 0, 0, 0, t99 * t53 + (t156 * t216 + t159 * t94) * t81 - t199, -t99 * t277 + (-t156 * t94 + t159 * t216) * t81 + t201 (t216 * t62 + t200) * t159 + (-t216 * t64 + t198) * t156 (-t10 * t216 + t211) * t159 + (t216 * t9 + t210) * t156 - t209 - t315; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, -t251 * t164, t236, -t231, -t237, qJDD(3), t160 * t322 + t100 + t297, g(3) * t160 + (-t112 - t322) * t157, 0, 0, t283, t317, t181, -t283, -t328, qJDD(3), t259 - t102 * t92 + (qJDD(3) * t154 - t246 * t89) * pkin(3) + t193, t51 * qJD(3) + t102 * t89 + (-qJDD(3) * t273 - t246 * t92) * pkin(3) - t19 - t174 (t43 - t50) * t92 + (-t42 + t51) * t89 + (t154 * t59 + t273 * t58) * pkin(3), t42 * t50 - t43 * t51 + (-t154 * t282 + t160 * t323 + t19 * t273 + t297) * pkin(3), qJDD(3), t33, t328, t283, t317, -t283, t124 * t58 - t127 * t59 + (-t38 - t50) * t92 + (t36 + t257) * t89, -t259 + t52 * t89 + (-pkin(4) + t127) * qJDD(3) + t175, t124 * qJDD(3) - t44 * t89 + t52 * t92 + (0.2e1 * qJD(5) - t51) * qJD(3) + t174 + t235, -t14 * t124 + t15 * t127 - t44 * t52 - t36 * t50 - g(1) * t229 - g(3) * (t255 - t299) + t257 * t38 - (-pkin(4) * t137 + t191) * t144, -t220 * t64 - t276, -t21 - t325 + (t23 + t223) * t156, t183 + t285, t218 * t62 + t275, t177 - t284, t81 * t89, -t11 * t81 + t124 * t24 + t165 * t156 + t159 * t311 + t256 * t62 + t9 * t89, -t10 * t89 + t12 * t81 - t124 * t23 - t156 * t311 + t165 * t159 + t256 * t64, t11 * t64 + t12 * t62 + (-t10 * t92 + t118 * t23 - t2 + (-t118 * t62 - t10) * qJD(6)) * t159 + (-t118 * t24 + t9 * t92 - t1 + (t118 * t64 + t9) * qJD(6)) * t156 - t228, t8 * t124 - t10 * t12 - t9 * t11 - g(1) * (pkin(8) * t266 + t229) - g(3) * (-t136 * t305 + t255) + t256 * t25 - (-t137 * t305 + t191) * t144 + (qJD(6) * t192 + t1 * t156 + t2 * t159) * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t327, t181 - 0.2e1 * t270, t318, t42 * t92 + t43 * t89 + t171, 0, 0, 0, 0, 0, 0, t318, -t327, t270 + t59, -t304 + t274 - t38 * t89 + (-qJD(5) - t36) * t92 + t171, 0, 0, 0, 0, 0, 0, t177 + t284, t156 * t81 ^ 2 + t285 + t53, -t156 * t320 - t21 + t325, -t156 * t321 + t159 * t206 + t25 * t89 - t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, qJDD(3) - t283, -t88 - t163, qJD(3) * t38 + t175 - t268, 0, 0, 0, 0, 0, 0, -qJD(3) * t62 + t183, -qJD(3) * t64 + t177, t156 * t319 + t159 * t320, -t25 * qJD(3) + t156 * t206 + t159 * t321 + t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t287, -t62 ^ 2 + t64 ^ 2, -t320, -t287, t319, -t56, -g(1) * t84 + g(2) * t82 - t125 * t159 - t25 * t64 + t321, g(1) * t85 - g(2) * t83 + t25 * t62 + t302 + (-qJD(6) * t22 - t5) * t159 + (qJD(6) * t26 + t125 - t7) * t156, 0, 0;];
tau_reg  = t6;
