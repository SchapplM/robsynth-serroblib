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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:13:28
% EndTime: 2020-01-03 12:13:32
% DurationCPUTime: 2.20s
% Computational Cost: add. (5359->342), mult. (7735->438), div. (0->0), fcn. (4670->16), ass. (0->217)
t170 = sin(qJ(4));
t165 = t170 ^ 2;
t175 = cos(qJ(4));
t166 = t175 ^ 2;
t253 = t165 + t166;
t162 = qJDD(1) + qJDD(2);
t151 = qJDD(3) + t162;
t164 = qJD(1) + qJD(2);
t177 = cos(qJ(2));
t274 = pkin(1) * qJD(1);
t242 = t177 * t274;
t107 = t164 * pkin(2) + t242;
t171 = sin(qJ(3));
t176 = cos(qJ(3));
t251 = qJD(3) * t176;
t240 = qJD(2) * t274;
t172 = sin(qJ(2));
t246 = qJDD(1) * t172;
t294 = pkin(1) * t246 + t177 * t240;
t281 = t177 * pkin(1);
t148 = qJDD(1) * t281;
t245 = t172 * t274;
t295 = -t162 * pkin(2) + qJD(3) * t245 + t172 * t240 - t148;
t202 = -t107 * t251 + t295 * t171 - t294 * t176;
t29 = t151 * pkin(8) - t202;
t225 = t253 * t29;
t168 = qJ(1) + qJ(2);
t157 = qJ(3) + t168;
t142 = sin(t157);
t143 = cos(t157);
t208 = -g(2) * t143 - g(3) * t142;
t252 = qJD(3) * t171;
t201 = t107 * t252 + t294 * t171 + t295 * t176;
t285 = t151 * pkin(3);
t30 = t201 - t285;
t196 = t208 - t30;
t133 = g(3) * t143;
t291 = g(2) * t142 - t133;
t174 = cos(qJ(5));
t247 = qJD(5) * t174;
t249 = qJD(4) * t175;
t293 = -t174 * t249 - t175 * t247;
t152 = qJD(3) + t164;
t169 = sin(qJ(5));
t259 = t174 * t170;
t95 = t169 * t175 + t259;
t77 = t95 * t152;
t244 = pkin(2) * t251;
t127 = t171 * t245;
t87 = t176 * t242 - t127;
t292 = t244 - t87;
t260 = t172 * t176;
t197 = t171 * t177 + t260;
t86 = t197 * t274;
t212 = pkin(2) * t252 - t86;
t163 = qJD(4) + qJD(5);
t144 = t171 * pkin(2) + pkin(8);
t282 = t176 * pkin(2);
t146 = -pkin(3) - t282;
t180 = qJD(4) ^ 2;
t290 = t144 * t180 + t146 * t151 + t212 * t152;
t234 = t152 * t249;
t262 = t170 * t151;
t74 = t171 * t107 + t176 * t245;
t66 = t152 * pkin(8) + t74;
t11 = -t66 * t249 + qJDD(4) * pkin(4) - t170 * t29 + (-t234 - t262) * pkin(9);
t250 = qJD(4) * t170;
t257 = t175 * t151;
t194 = t152 * t250 - t257;
t12 = -t194 * pkin(9) + t175 * t29 - t66 * t250;
t248 = qJD(5) * t169;
t231 = pkin(9) * t152 + t66;
t52 = t231 * t170;
t51 = qJD(4) * pkin(4) - t52;
t53 = t231 * t175;
t3 = (qJD(5) * t51 + t12) * t174 + t169 * t11 - t53 * t248;
t179 = -pkin(9) - pkin(8);
t147 = pkin(2) + t281;
t89 = pkin(1) * t260 + t171 * t147;
t85 = pkin(8) + t89;
t289 = -pkin(9) - t85;
t288 = g(1) * t175;
t284 = t152 * pkin(3);
t283 = t175 * pkin(4);
t263 = t169 * t170;
t239 = t152 * t263;
t258 = t174 * t175;
t75 = -t152 * t258 + t239;
t280 = t77 * t75;
t279 = -pkin(9) - t144;
t92 = t279 * t170;
t159 = t175 * pkin(9);
t93 = t175 * t144 + t159;
t58 = -t169 * t93 + t174 * t92;
t227 = qJD(4) * t279;
t71 = t170 * t227 + t175 * t244;
t72 = -t170 * t244 + t175 * t227;
t94 = -t258 + t263;
t278 = t58 * qJD(5) + t169 * t72 + t174 * t71 + t94 * t87;
t59 = t169 * t92 + t174 * t93;
t277 = -t59 * qJD(5) - t169 * t71 + t174 * t72 + t95 * t87;
t235 = qJD(4) * t179;
t100 = t170 * t235;
t101 = t175 * t235;
t121 = t179 * t170;
t122 = t175 * pkin(8) + t159;
t69 = t174 * t121 - t169 * t122;
t73 = t176 * t107 - t127;
t276 = t69 * qJD(5) + t174 * t100 + t169 * t101 + t94 * t73;
t70 = t169 * t121 + t174 * t122;
t275 = -t70 * qJD(5) - t169 * t100 + t174 * t101 + t95 * t73;
t273 = t169 * t53;
t272 = t174 * t53;
t261 = t171 * t172;
t56 = t147 * t251 + (-t172 * t252 + (t176 * t177 - t261) * qJD(2)) * pkin(1);
t271 = t56 * t152;
t57 = t147 * t252 + (t197 * qJD(2) + t172 * t251) * pkin(1);
t270 = t57 * t152;
t269 = t74 * t152;
t149 = pkin(4) * t250;
t268 = t149 + t212;
t145 = pkin(3) + t283;
t267 = t142 * t145 + t143 * t179;
t167 = qJ(4) + qJ(5);
t153 = sin(t167);
t266 = t142 * t153;
t265 = t143 * t153;
t264 = t152 * t170;
t256 = t143 * pkin(3) + t142 * pkin(8);
t254 = t165 - t166;
t154 = sin(t168);
t139 = pkin(2) * t154;
t241 = t139 + t267;
t150 = t152 ^ 2;
t238 = t170 * t150 * t175;
t156 = cos(t168);
t140 = pkin(2) * t156;
t236 = t140 + t256;
t230 = qJD(4) * t289;
t229 = g(2) * t154 - g(3) * t156;
t226 = t73 * t253;
t223 = -t142 * t179 + t143 * t145;
t222 = t253 * t151;
t88 = -pkin(1) * t261 + t176 * t147;
t23 = t194 * pkin(4) + t30;
t55 = -t145 * t152 - t73;
t200 = t163 * t263;
t63 = t200 + t293;
t221 = g(2) * t265 + g(3) * t266 + t23 * t95 - t55 * t63;
t220 = -t151 * t259 + t293 * t152 - t169 * t257;
t65 = -t73 - t284;
t219 = -t196 * t170 + t65 * t249;
t218 = -t291 + t225;
t217 = qJD(1) * (-qJD(2) + t164);
t216 = qJD(2) * (-qJD(1) - t164);
t213 = t170 * t234;
t211 = -t74 + t149;
t84 = -pkin(3) - t88;
t131 = t142 * pkin(3);
t210 = -t143 * pkin(8) + t131 + t139;
t209 = t140 + t223;
t206 = -g(2) * t156 - g(3) * t154;
t173 = sin(qJ(1));
t178 = cos(qJ(1));
t205 = -g(2) * t178 - g(3) * t173;
t204 = t169 * t262 - t174 * t257;
t20 = t169 * t51 + t272;
t67 = t289 * t170;
t68 = t175 * t85 + t159;
t44 = -t169 * t68 + t174 * t67;
t45 = t169 * t67 + t174 * t68;
t19 = t174 * t51 - t273;
t4 = -t20 * qJD(5) + t174 * t11 - t169 * t12;
t64 = t163 * t95;
t199 = t19 * t63 - t20 * t64 - t3 * t94 - t4 * t95 - t291;
t195 = t148 + t206;
t193 = t202 + t291;
t192 = pkin(8) * t180 - t269 - t285;
t191 = t151 * t84 + t180 * t85 + t270;
t190 = -t152 * t65 - t29 + t291;
t189 = -pkin(8) * qJDD(4) + (t73 - t284) * qJD(4);
t188 = -qJDD(4) * t85 + (t152 * t84 - t56) * qJD(4);
t155 = cos(t167);
t187 = t208 * t155 + t23 * t94 + t55 * t64;
t186 = -t201 + t208;
t185 = t292 * t253;
t184 = -qJDD(4) * t144 + (t146 * t152 - t292) * qJD(4);
t183 = g(1) * t153 + t291 * t155 + t55 * t75 - t3;
t182 = -g(1) * t155 + g(2) * t266 - g(3) * t265 - t55 * t77 + t4;
t161 = qJDD(4) + qJDD(5);
t160 = t178 * pkin(1);
t158 = t173 * pkin(1);
t114 = -t145 - t282;
t113 = qJDD(4) * t175 - t180 * t170;
t112 = qJDD(4) * t170 + t180 * t175;
t82 = t166 * t151 - 0.2e1 * t213;
t81 = t165 * t151 + 0.2e1 * t213;
t80 = t84 - t283;
t62 = -0.2e1 * t254 * t152 * qJD(4) + 0.2e1 * t170 * t257;
t60 = t65 * t250;
t54 = t149 + t57;
t47 = -t94 * t161 - t64 * t163;
t46 = t95 * t161 - t63 * t163;
t37 = -t170 * t56 + t175 * t230;
t36 = t170 * t230 + t175 * t56;
t35 = -t75 ^ 2 + t77 ^ 2;
t32 = t64 * t152 + t204;
t31 = t152 * t200 + t220;
t25 = -t174 * t52 - t273;
t24 = t169 * t52 - t272;
t17 = -t220 + (-t239 + t75) * t163;
t9 = t32 * t94 + t75 * t64;
t8 = -t31 * t95 - t77 * t63;
t7 = -t45 * qJD(5) - t169 * t36 + t174 * t37;
t6 = t44 * qJD(5) + t169 * t37 + t174 * t36;
t5 = t31 * t94 - t95 * t32 + t63 * t75 - t77 * t64;
t1 = [0, 0, 0, 0, 0, qJDD(1), t205, g(2) * t173 - g(3) * t178, 0, 0, 0, 0, 0, 0, 0, t162, (t162 * t177 + t172 * t216) * pkin(1) + t195, ((-qJDD(1) - t162) * t172 + t177 * t216) * pkin(1) + t229, 0, (t205 + (t172 ^ 2 + t177 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t151, t88 * t151 + t186 - t270, -t89 * t151 + t193 - t271, 0, -t202 * t89 + t74 * t56 - t201 * t88 - t73 * t57 - g(2) * (t140 + t160) - g(3) * (t139 + t158), t81, t62, t112, t82, t113, 0, t60 + t188 * t170 + (-t191 + t196) * t175, t191 * t170 + t188 * t175 + t219, t85 * t222 + t253 * t271 + t218, t30 * t84 + t65 * t57 - g(2) * (t160 + t236) - g(3) * (t158 + t210) + t253 * (t29 * t85 + t66 * t56), t8, t5, t46, t9, t47, 0, t44 * t161 + t7 * t163 + t80 * t32 + t54 * t75 + t187, -t45 * t161 - t6 * t163 - t80 * t31 + t54 * t77 + t221, t44 * t31 - t45 * t32 - t6 * t75 - t7 * t77 + t199, t3 * t45 + t20 * t6 + t4 * t44 + t19 * t7 + t23 * t80 + t55 * t54 - g(2) * (t160 + t209) - g(3) * (t158 + t241); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, t172 * pkin(1) * t217 + t195, (t177 * t217 - t246) * pkin(1) + t229, 0, 0, 0, 0, 0, 0, 0, t151, t86 * t152 + (t151 * t176 - t152 * t252) * pkin(2) + t186, t87 * t152 + (-t151 * t171 - t152 * t251) * pkin(2) + t193, 0, t73 * t86 - t74 * t87 + (-t171 * t202 - t176 * t201 + (-t171 * t73 + t176 * t74) * qJD(3) + t206) * pkin(2), t81, t62, t112, t82, t113, 0, t60 + t184 * t170 + (t196 - t290) * t175, t290 * t170 + t184 * t175 + t219, t144 * t222 + t185 * t152 + t218, -g(2) * t236 - g(3) * t210 + t144 * t225 + t30 * t146 + t185 * t66 + t212 * t65, t8, t5, t46, t9, t47, 0, t114 * t32 + t58 * t161 + t277 * t163 + t268 * t75 + t187, -t114 * t31 - t59 * t161 - t278 * t163 + t268 * t77 + t221, -t277 * t77 - t278 * t75 + t58 * t31 - t59 * t32 + t199, -g(2) * t209 - g(3) * t241 + t23 * t114 + t277 * t19 + t278 * t20 + t268 * t55 + t3 * t59 + t4 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t186 + t269, t73 * t152 + t193, 0, 0, t81, t62, t112, t82, t113, 0, t60 + t189 * t170 + (-t192 + t196) * t175, t192 * t170 + t189 * t175 + t219, pkin(8) * t222 - t152 * t226 + t218, -t30 * pkin(3) - t65 * t74 - g(2) * t256 - g(3) * t131 - t66 * t226 + (t225 + t133) * pkin(8), t8, t5, t46, t9, t47, 0, -t145 * t32 + t69 * t161 + t275 * t163 + t211 * t75 + t187, t145 * t31 - t70 * t161 - t276 * t163 + t211 * t77 + t221, -t275 * t77 - t276 * t75 + t69 * t31 - t70 * t32 + t199, -g(2) * t223 - g(3) * t267 - t23 * t145 + t275 * t19 + t276 * t20 + t211 * t55 + t3 * t70 + t4 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, t254 * t150, t262, t238, t257, qJDD(4), t190 * t170 - t288, g(1) * t170 + t190 * t175, 0, 0, t280, t35, t17, -t280, -t204, t161, -t24 * t163 + (t161 * t174 - t163 * t248 - t264 * t75) * pkin(4) + t182, t25 * t163 + (-t161 * t169 - t163 * t247 - t264 * t77) * pkin(4) + t183, (t20 + t24) * t77 + (-t19 + t25) * t75 + (-t169 * t32 + t174 * t31 + (t169 * t77 - t174 * t75) * qJD(5)) * pkin(4), -t19 * t24 - t20 * t25 + (-t288 + t169 * t3 + t174 * t4 + (-t169 * t19 + t174 * t20) * qJD(5) + (-t152 * t55 + t291) * t170) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t280, t35, t17, -t280, -t204, t161, t20 * t163 + t182, t19 * t163 + t183, 0, 0;];
tau_reg = t1;
