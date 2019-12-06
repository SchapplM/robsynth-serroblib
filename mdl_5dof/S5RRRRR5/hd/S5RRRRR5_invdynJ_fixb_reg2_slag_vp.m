% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:35
% EndTime: 2019-12-05 18:58:41
% DurationCPUTime: 2.19s
% Computational Cost: add. (5359->345), mult. (7735->438), div. (0->0), fcn. (4670->16), ass. (0->215)
t166 = sin(qJ(4));
t161 = t166 ^ 2;
t171 = cos(qJ(4));
t162 = t171 ^ 2;
t247 = t161 + t162;
t158 = qJDD(1) + qJDD(2);
t149 = qJDD(3) + t158;
t160 = qJD(1) + qJD(2);
t173 = cos(qJ(2));
t266 = pkin(1) * qJD(1);
t236 = t173 * t266;
t105 = pkin(2) * t160 + t236;
t167 = sin(qJ(3));
t172 = cos(qJ(3));
t245 = qJD(3) * t172;
t234 = qJD(2) * t266;
t168 = sin(qJ(2));
t240 = qJDD(1) * t168;
t289 = -pkin(1) * t240 - t173 * t234;
t275 = t173 * pkin(1);
t146 = qJDD(1) * t275;
t239 = t168 * t266;
t290 = -pkin(2) * t158 + qJD(3) * t239 + t168 * t234 - t146;
t200 = -t105 * t245 + t167 * t290 + t172 * t289;
t29 = pkin(8) * t149 - t200;
t219 = t247 * t29;
t164 = qJ(1) + qJ(2);
t155 = qJ(3) + t164;
t140 = sin(t155);
t141 = cos(t155);
t205 = g(2) * t141 + g(3) * t140;
t131 = g(2) * t140;
t249 = -g(3) * t141 + t131;
t170 = cos(qJ(5));
t241 = qJD(5) * t170;
t243 = qJD(4) * t171;
t288 = -t170 * t243 - t171 * t241;
t150 = qJD(3) + t160;
t165 = sin(qJ(5));
t252 = t170 * t166;
t95 = t165 * t171 + t252;
t77 = t95 * t150;
t238 = pkin(2) * t245;
t125 = t167 * t239;
t87 = t172 * t236 - t125;
t287 = t238 - t87;
t246 = qJD(3) * t167;
t253 = t168 * t172;
t195 = t167 * t173 + t253;
t86 = t195 * t266;
t207 = pkin(2) * t246 - t86;
t152 = sin(t164);
t154 = cos(t164);
t286 = g(2) * t154 + g(3) * t152;
t159 = qJD(4) + qJD(5);
t142 = pkin(2) * t167 + pkin(8);
t276 = t172 * pkin(2);
t144 = -pkin(3) - t276;
t176 = qJD(4) ^ 2;
t285 = t142 * t176 + t144 * t149 + t150 * t207;
t228 = t150 * t243;
t255 = t166 * t149;
t74 = t105 * t167 + t172 * t239;
t66 = pkin(8) * t150 + t74;
t11 = -t66 * t243 + qJDD(4) * pkin(4) - t166 * t29 + (-t228 - t255) * pkin(9);
t244 = qJD(4) * t166;
t250 = t171 * t149;
t192 = t150 * t244 - t250;
t12 = -pkin(9) * t192 + t171 * t29 - t244 * t66;
t242 = qJD(5) * t165;
t225 = pkin(9) * t150 + t66;
t52 = t225 * t166;
t51 = qJD(4) * pkin(4) - t52;
t53 = t225 * t171;
t3 = (qJD(5) * t51 + t12) * t170 + t165 * t11 - t53 * t242;
t175 = -pkin(9) - pkin(8);
t145 = pkin(2) + t275;
t89 = pkin(1) * t253 + t145 * t167;
t85 = pkin(8) + t89;
t284 = -pkin(9) - t85;
t283 = pkin(2) * t152;
t282 = pkin(2) * t154;
t281 = g(1) * t171;
t280 = t149 * pkin(3);
t279 = t150 * pkin(3);
t169 = sin(qJ(1));
t278 = t169 * pkin(1);
t277 = t171 * pkin(4);
t174 = cos(qJ(1));
t274 = t174 * pkin(1);
t256 = t165 * t166;
t233 = t150 * t256;
t251 = t170 * t171;
t75 = -t150 * t251 + t233;
t273 = t77 * t75;
t272 = -pkin(9) - t142;
t92 = t272 * t166;
t156 = t171 * pkin(9);
t93 = t142 * t171 + t156;
t58 = -t165 * t93 + t170 * t92;
t221 = qJD(4) * t272;
t71 = t166 * t221 + t171 * t238;
t72 = -t166 * t238 + t171 * t221;
t94 = -t251 + t256;
t271 = qJD(5) * t58 + t165 * t72 + t170 * t71 + t87 * t94;
t59 = t165 * t92 + t170 * t93;
t270 = -qJD(5) * t59 - t165 * t71 + t170 * t72 + t87 * t95;
t39 = -t105 * t246 + t167 * t289 - t172 * t290;
t30 = -t39 - t280;
t73 = t172 * t105 - t125;
t65 = -t73 - t279;
t269 = t166 * t30 + t243 * t65;
t119 = t175 * t166;
t120 = pkin(8) * t171 + t156;
t69 = t119 * t170 - t120 * t165;
t229 = qJD(4) * t175;
t98 = t166 * t229;
t99 = t171 * t229;
t268 = qJD(5) * t69 + t165 * t99 + t170 * t98 + t73 * t94;
t70 = t119 * t165 + t120 * t170;
t267 = -qJD(5) * t70 - t165 * t98 + t170 * t99 + t73 * t95;
t265 = t165 * t53;
t264 = t170 * t53;
t254 = t167 * t168;
t56 = t145 * t245 + (-t168 * t246 + (t172 * t173 - t254) * qJD(2)) * pkin(1);
t263 = t56 * t150;
t57 = t145 * t246 + (qJD(2) * t195 + t168 * t245) * pkin(1);
t262 = t57 * t150;
t261 = t74 * t150;
t147 = pkin(4) * t244;
t260 = t147 + t207;
t163 = qJ(4) + qJ(5);
t153 = cos(t163);
t259 = t140 * t153;
t258 = t141 * t153;
t257 = t150 * t166;
t248 = t161 - t162;
t235 = t171 * t205 + t65 * t244;
t148 = t150 ^ 2;
t232 = t166 * t148 * t171;
t230 = t146 + t286;
t143 = pkin(3) + t277;
t224 = qJD(4) * t284;
t223 = -g(2) * t152 + g(3) * t154;
t220 = t73 * t247;
t217 = t247 * t149;
t216 = t140 * t175 - t141 * t143;
t88 = -pkin(1) * t254 + t145 * t172;
t23 = pkin(4) * t192 + t30;
t55 = -t143 * t150 - t73;
t64 = t159 * t95;
t215 = g(2) * t258 + g(3) * t259 + t23 * t94 + t55 * t64;
t214 = -t149 * t252 + t150 * t288 - t165 * t250;
t213 = t249 + t219;
t212 = qJD(1) * (-qJD(2) + t160);
t211 = qJD(2) * (-qJD(1) - t160);
t208 = t166 * t228;
t206 = -t74 + t147;
t84 = -pkin(3) - t88;
t203 = g(2) * t174 + g(3) * t169;
t202 = t165 * t255 - t170 * t250;
t199 = t159 * t256;
t20 = t165 * t51 + t264;
t67 = t284 * t166;
t68 = t171 * t85 + t156;
t44 = -t165 * t68 + t170 * t67;
t45 = t165 * t67 + t170 * t68;
t19 = t170 * t51 - t265;
t4 = -qJD(5) * t20 + t11 * t170 - t165 * t12;
t63 = t199 + t288;
t198 = t19 * t63 - t20 * t64 - t3 * t94 - t4 * t95 + t249;
t197 = -t140 * t143 - t141 * t175;
t128 = t141 * pkin(8);
t194 = -pkin(3) * t140 + t128 - t283;
t193 = t216 - t282;
t191 = t200 - t249;
t190 = t39 + t205;
t189 = -pkin(3) * t141 - pkin(8) * t140 - t282;
t188 = -pkin(8) * t176 + t261 + t280;
t187 = -t149 * t84 - t176 * t85 - t262;
t186 = t197 - t283;
t185 = -t150 * t65 - t249 - t29;
t184 = -pkin(8) * qJDD(4) + (t73 - t279) * qJD(4);
t183 = -qJDD(4) * t85 + (t150 * t84 - t56) * qJD(4);
t151 = sin(t163);
t182 = -t151 * t205 + t23 * t95 - t55 * t63;
t181 = t287 * t247;
t180 = -qJDD(4) * t142 + (t144 * t150 - t287) * qJD(4);
t179 = g(1) * t151 - g(2) * t259 + g(3) * t258 + t55 * t75 - t3;
t178 = -g(1) * t153 - t151 * t249 - t55 * t77 + t4;
t157 = qJDD(4) + qJDD(5);
t112 = -t143 - t276;
t111 = qJDD(4) * t171 - t166 * t176;
t110 = qJDD(4) * t166 + t171 * t176;
t82 = t149 * t162 - 0.2e1 * t208;
t81 = t149 * t161 + 0.2e1 * t208;
t80 = t84 - t277;
t62 = -0.2e1 * qJD(4) * t150 * t248 + 0.2e1 * t166 * t250;
t54 = t147 + t57;
t47 = -t157 * t94 - t159 * t64;
t46 = t157 * t95 - t159 * t63;
t37 = -t166 * t56 + t171 * t224;
t36 = t166 * t224 + t171 * t56;
t35 = -t75 ^ 2 + t77 ^ 2;
t32 = t150 * t64 + t202;
t31 = t150 * t199 + t214;
t25 = -t170 * t52 - t265;
t24 = t165 * t52 - t264;
t17 = -t214 + (-t233 + t75) * t159;
t9 = t32 * t94 + t64 * t75;
t8 = -t31 * t95 - t63 * t77;
t7 = -qJD(5) * t45 - t165 * t36 + t170 * t37;
t6 = qJD(5) * t44 + t165 * t37 + t170 * t36;
t5 = t31 * t94 - t32 * t95 + t63 * t75 - t64 * t77;
t1 = [0, 0, 0, 0, 0, qJDD(1), t203, -g(2) * t169 + g(3) * t174, 0, 0, 0, 0, 0, 0, 0, t158, (t158 * t173 + t168 * t211) * pkin(1) + t230, ((-qJDD(1) - t158) * t168 + t173 * t211) * pkin(1) + t223, 0, (t203 + (t168 ^ 2 + t173 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t149, t149 * t88 + t190 - t262, -t149 * t89 + t191 - t263, 0, -t200 * t89 + t74 * t56 + t39 * t88 - t73 * t57 - g(2) * (-t274 - t282) - g(3) * (-t278 - t283), t81, t62, t110, t82, t111, 0, t183 * t166 + (t187 - t30) * t171 + t235, t183 * t171 + (-t187 - t205) * t166 + t269, t217 * t85 + t247 * t263 + t213, t30 * t84 + t65 * t57 - g(2) * (t189 - t274) - g(3) * (t194 - t278) + t247 * (t29 * t85 + t56 * t66), t8, t5, t46, t9, t47, 0, t157 * t44 + t159 * t7 + t32 * t80 + t54 * t75 + t215, -t45 * t157 - t6 * t159 - t80 * t31 + t54 * t77 + t182, t31 * t44 - t32 * t45 - t6 * t75 - t7 * t77 + t198, t3 * t45 + t20 * t6 + t4 * t44 + t19 * t7 + t23 * t80 + t55 * t54 - g(2) * (t193 - t274) - g(3) * (t186 - t278); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, pkin(1) * t168 * t212 + t230, (t173 * t212 - t240) * pkin(1) + t223, 0, 0, 0, 0, 0, 0, 0, t149, t86 * t150 + (t149 * t172 - t150 * t246) * pkin(2) + t190, t87 * t150 + (-t149 * t167 - t150 * t245) * pkin(2) + t191, 0, t73 * t86 - t74 * t87 + (-t167 * t200 + t172 * t39 + (-t167 * t73 + t172 * t74) * qJD(3) + t286) * pkin(2), t81, t62, t110, t82, t111, 0, t180 * t166 + (-t285 - t30) * t171 + t235, t180 * t171 + (-t205 + t285) * t166 + t269, t142 * t217 + t150 * t181 + t213, -g(2) * t189 - g(3) * t194 + t142 * t219 + t30 * t144 + t181 * t66 + t207 * t65, t8, t5, t46, t9, t47, 0, t112 * t32 + t58 * t157 + t159 * t270 + t260 * t75 + t215, -t112 * t31 - t59 * t157 - t159 * t271 + t260 * t77 + t182, -t270 * t77 - t271 * t75 + t58 * t31 - t59 * t32 + t198, -g(2) * t193 - g(3) * t186 + t23 * t112 + t19 * t270 + t20 * t271 + t260 * t55 + t3 * t59 + t4 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, t190 + t261, t150 * t73 + t191, 0, 0, t81, t62, t110, t82, t111, 0, t184 * t166 + (t188 - t30) * t171 + t235, t184 * t171 + (-t188 - t205) * t166 + t269, pkin(8) * t217 - t150 * t220 + t213, -g(3) * t128 - t65 * t74 - t66 * t220 + (t205 - t30) * pkin(3) + (t219 + t131) * pkin(8), t8, t5, t46, t9, t47, 0, -t143 * t32 + t69 * t157 + t159 * t267 + t206 * t75 + t215, t143 * t31 - t70 * t157 - t159 * t268 + t206 * t77 + t182, -t267 * t77 - t268 * t75 + t69 * t31 - t70 * t32 + t198, -g(2) * t216 - g(3) * t197 - t23 * t143 + t19 * t267 + t20 * t268 + t206 * t55 + t3 * t70 + t4 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t232, t248 * t148, t255, t232, t250, qJDD(4), t166 * t185 - t281, g(1) * t166 + t171 * t185, 0, 0, t273, t35, t17, -t273, -t202, t157, -t24 * t159 + (t157 * t170 - t159 * t242 - t257 * t75) * pkin(4) + t178, t25 * t159 + (-t157 * t165 - t159 * t241 - t257 * t77) * pkin(4) + t179, (t20 + t24) * t77 + (-t19 + t25) * t75 + (-t165 * t32 + t170 * t31 + (t165 * t77 - t170 * t75) * qJD(5)) * pkin(4), -t19 * t24 - t20 * t25 + (-t281 + t165 * t3 + t170 * t4 + (-t165 * t19 + t170 * t20) * qJD(5) + (-t150 * t55 - t249) * t166) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, t35, t17, -t273, -t202, t157, t20 * t159 + t178, t19 * t159 + t179, 0, 0;];
tau_reg = t1;
