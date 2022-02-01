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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:02:07
% EndTime: 2022-01-20 12:02:15
% DurationCPUTime: 2.43s
% Computational Cost: add. (5359->346), mult. (7735->443), div. (0->0), fcn. (4670->16), ass. (0->219)
t171 = sin(qJ(4));
t166 = t171 ^ 2;
t176 = cos(qJ(4));
t167 = t176 ^ 2;
t253 = t166 + t167;
t163 = qJDD(1) + qJDD(2);
t153 = qJDD(3) + t163;
t165 = qJD(1) + qJD(2);
t178 = cos(qJ(2));
t277 = pkin(1) * qJD(1);
t241 = t178 * t277;
t108 = t165 * pkin(2) + t241;
t172 = sin(qJ(3));
t177 = cos(qJ(3));
t251 = qJD(3) * t177;
t240 = qJD(2) * t277;
t173 = sin(qJ(2));
t246 = qJDD(1) * t173;
t300 = pkin(1) * t246 + t178 * t240;
t284 = t178 * pkin(1);
t150 = qJDD(1) * t284;
t244 = t173 * t277;
t301 = -t163 * pkin(2) + qJD(3) * t244 + t173 * t240 - t150;
t205 = -t108 * t251 + t301 * t172 - t300 * t177;
t29 = t153 * pkin(8) - t205;
t224 = t253 * t29;
t169 = qJ(1) + qJ(2);
t159 = qJ(3) + t169;
t145 = cos(t159);
t292 = g(2) * t145;
t252 = qJD(3) * t172;
t204 = t108 * t252 + t300 * t172 + t301 * t177;
t289 = t153 * pkin(3);
t30 = t204 - t289;
t230 = -t30 - t292;
t175 = cos(qJ(5));
t247 = qJD(5) * t175;
t249 = qJD(4) * t176;
t299 = -t175 * t249 - t176 * t247;
t154 = qJD(3) + t165;
t170 = sin(qJ(5));
t260 = t175 * t171;
t95 = t170 * t176 + t260;
t77 = t95 * t154;
t243 = pkin(2) * t251;
t127 = t172 * t244;
t87 = t177 * t241 - t127;
t298 = t243 - t87;
t261 = t173 * t177;
t199 = t172 * t178 + t261;
t86 = t199 * t277;
t213 = pkin(2) * t252 - t86;
t144 = sin(t159);
t297 = g(1) * t145 + g(2) * t144;
t156 = sin(t169);
t158 = cos(t169);
t296 = g(1) * t156 - g(2) * t158;
t164 = qJD(4) + qJD(5);
t146 = t172 * pkin(2) + pkin(8);
t285 = t177 * pkin(2);
t148 = -pkin(3) - t285;
t181 = qJD(4) ^ 2;
t295 = t146 * t181 + t148 * t153 + t154 * t213;
t234 = t154 * t249;
t263 = t171 * t153;
t74 = t172 * t108 + t177 * t244;
t66 = t154 * pkin(8) + t74;
t11 = -t66 * t249 + qJDD(4) * pkin(4) - t171 * t29 + (-t234 - t263) * pkin(9);
t250 = qJD(4) * t171;
t258 = t176 * t153;
t195 = t154 * t250 - t258;
t12 = -pkin(9) * t195 + t176 * t29 - t250 * t66;
t248 = qJD(5) * t170;
t231 = pkin(9) * t154 + t66;
t52 = t231 * t171;
t51 = qJD(4) * pkin(4) - t52;
t53 = t231 * t176;
t3 = (qJD(5) * t51 + t12) * t175 + t170 * t11 - t53 * t248;
t180 = -pkin(9) - pkin(8);
t149 = pkin(2) + t284;
t89 = pkin(1) * t261 + t172 * t149;
t85 = pkin(8) + t89;
t294 = -pkin(9) - t85;
t293 = pkin(2) * t156;
t134 = g(1) * t144;
t290 = g(3) * t176;
t288 = t154 * pkin(3);
t174 = sin(qJ(1));
t287 = t174 * pkin(1);
t286 = t176 * pkin(4);
t264 = t170 * t171;
t239 = t154 * t264;
t259 = t175 * t176;
t75 = -t154 * t259 + t239;
t283 = t77 * t75;
t282 = -pkin(9) - t146;
t92 = t282 * t171;
t160 = t176 * pkin(9);
t93 = t176 * t146 + t160;
t58 = -t170 * t93 + t175 * t92;
t226 = qJD(4) * t282;
t71 = t171 * t226 + t176 * t243;
t72 = -t171 * t243 + t176 * t226;
t94 = -t259 + t264;
t281 = qJD(5) * t58 + t170 * t72 + t175 * t71 + t94 * t87;
t59 = t170 * t92 + t175 * t93;
t280 = -qJD(5) * t59 - t170 * t71 + t175 * t72 + t95 * t87;
t235 = qJD(4) * t180;
t100 = t176 * t235;
t122 = t180 * t171;
t123 = t176 * pkin(8) + t160;
t69 = t175 * t122 - t170 * t123;
t73 = t177 * t108 - t127;
t99 = t171 * t235;
t279 = qJD(5) * t69 + t170 * t100 + t175 * t99 + t94 * t73;
t70 = t170 * t122 + t175 * t123;
t278 = -qJD(5) * t70 + t175 * t100 - t170 * t99 + t95 * t73;
t276 = t170 * t53;
t275 = t175 * t53;
t262 = t172 * t173;
t56 = t149 * t251 + (-t173 * t252 + (t177 * t178 - t262) * qJD(2)) * pkin(1);
t274 = t56 * t154;
t57 = t149 * t252 + (qJD(2) * t199 + t173 * t251) * pkin(1);
t273 = t57 * t154;
t272 = t74 * t154;
t151 = pkin(4) * t250;
t271 = t151 + t213;
t65 = -t73 - t288;
t270 = t176 * t134 + t65 * t250;
t168 = qJ(4) + qJ(5);
t155 = sin(t168);
t269 = t144 * t155;
t157 = cos(t168);
t268 = t144 * t157;
t267 = t145 * t155;
t266 = t145 * t157;
t265 = t154 * t171;
t257 = t145 * pkin(3) + t144 * pkin(8);
t255 = g(1) * t158 + g(2) * t156;
t254 = t166 - t167;
t245 = -t230 * t171 + t65 * t249;
t152 = t154 ^ 2;
t238 = t171 * t152 * t176;
t142 = pkin(2) * t158;
t236 = t142 + t257;
t147 = pkin(3) + t286;
t229 = qJD(4) * t294;
t228 = -t144 * pkin(3) + t145 * pkin(8);
t225 = t253 * t73;
t222 = -t144 * t180 + t145 * t147;
t221 = t253 * t153;
t88 = -pkin(1) * t262 + t177 * t149;
t220 = -t153 * t260 + t299 * t154 - t170 * t258;
t219 = -t297 + t224;
t218 = qJD(1) * (-qJD(2) + t165);
t217 = qJD(2) * (-qJD(1) - t165);
t214 = t171 * t234;
t212 = -t74 + t151;
t84 = -pkin(3) - t88;
t211 = t150 + t296;
t210 = t142 + t222;
t179 = cos(qJ(1));
t208 = g(1) * t174 - g(2) * t179;
t207 = t170 * t263 - t175 * t258;
t203 = t164 * t264;
t20 = t170 * t51 + t275;
t67 = t294 * t171;
t68 = t176 * t85 + t160;
t44 = -t170 * t68 + t175 * t67;
t45 = t170 * t67 + t175 * t68;
t19 = t175 * t51 - t276;
t4 = -qJD(5) * t20 + t175 * t11 - t170 * t12;
t63 = t203 + t299;
t64 = t164 * t95;
t202 = t19 * t63 - t20 * t64 - t3 * t94 - t4 * t95 - t297;
t201 = -t144 * t147 - t145 * t180;
t198 = t228 - t293;
t23 = pkin(4) * t195 + t30;
t55 = -t147 * t154 - t73;
t197 = -g(1) * t269 + g(2) * t267 + t23 * t95 - t55 * t63;
t196 = g(1) * t268 - g(2) * t266 + t23 * t94 + t55 * t64;
t194 = t205 + t297;
t193 = pkin(8) * t181 - t272 - t289;
t192 = t153 * t84 + t181 * t85 + t273;
t191 = t201 - t293;
t190 = -t154 * t65 - t29 + t297;
t189 = -pkin(8) * qJDD(4) + (t73 - t288) * qJD(4);
t188 = t134 - t204 - t292;
t187 = -qJDD(4) * t85 + (t154 * t84 - t56) * qJD(4);
t186 = t298 * t253;
t185 = -qJDD(4) * t146 + (t148 * t154 - t298) * qJD(4);
t184 = g(1) * t266 + g(2) * t268 + g(3) * t155 + t55 * t75 - t3;
t183 = g(1) * t267 + g(2) * t269 - g(3) * t157 - t55 * t77 + t4;
t162 = qJDD(4) + qJDD(5);
t161 = t179 * pkin(1);
t115 = -t147 - t285;
t114 = qJDD(4) * t176 - t181 * t171;
t113 = qJDD(4) * t171 + t181 * t176;
t82 = t167 * t153 - 0.2e1 * t214;
t81 = t166 * t153 + 0.2e1 * t214;
t80 = t84 - t286;
t62 = -0.2e1 * qJD(4) * t154 * t254 + 0.2e1 * t171 * t258;
t54 = t151 + t57;
t47 = -t94 * t162 - t64 * t164;
t46 = t95 * t162 - t63 * t164;
t37 = -t171 * t56 + t176 * t229;
t36 = t171 * t229 + t176 * t56;
t35 = -t75 ^ 2 + t77 ^ 2;
t32 = t154 * t64 + t207;
t31 = t154 * t203 + t220;
t25 = -t175 * t52 - t276;
t24 = t170 * t52 - t275;
t17 = -t220 + (-t239 + t75) * t164;
t9 = t32 * t94 + t75 * t64;
t8 = -t31 * t95 - t77 * t63;
t7 = -qJD(5) * t45 - t170 * t36 + t175 * t37;
t6 = qJD(5) * t44 + t170 * t37 + t175 * t36;
t5 = t31 * t94 - t95 * t32 + t63 * t75 - t77 * t64;
t1 = [0, 0, 0, 0, 0, qJDD(1), t208, g(1) * t179 + g(2) * t174, 0, 0, 0, 0, 0, 0, 0, t163, (t163 * t178 + t173 * t217) * pkin(1) + t211, ((-qJDD(1) - t163) * t173 + t178 * t217) * pkin(1) + t255, 0, (t208 + (t173 ^ 2 + t178 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t153, t88 * t153 + t188 - t273, -t89 * t153 + t194 - t274, 0, -t205 * t89 + t74 * t56 - t204 * t88 - t73 * t57 - g(1) * (-t287 - t293) - g(2) * (t142 + t161), t81, t62, t113, t82, t114, 0, t187 * t171 + (-t192 + t230) * t176 + t270, t187 * t176 + (t192 - t134) * t171 + t245, t221 * t85 + t253 * t274 + t219, t30 * t84 + t65 * t57 - g(1) * (t198 - t287) - g(2) * (t161 + t236) + t253 * (t29 * t85 + t66 * t56), t8, t5, t46, t9, t47, 0, t44 * t162 + t7 * t164 + t80 * t32 + t54 * t75 + t196, -t45 * t162 - t6 * t164 - t80 * t31 + t54 * t77 + t197, t44 * t31 - t45 * t32 - t6 * t75 - t7 * t77 + t202, t3 * t45 + t20 * t6 + t4 * t44 + t19 * t7 + t23 * t80 + t55 * t54 - g(1) * (t191 - t287) - g(2) * (t161 + t210); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, pkin(1) * t173 * t218 + t211, (t178 * t218 - t246) * pkin(1) + t255, 0, 0, 0, 0, 0, 0, 0, t153, t86 * t154 + (t153 * t177 - t154 * t252) * pkin(2) + t188, t87 * t154 + (-t153 * t172 - t154 * t251) * pkin(2) + t194, 0, t73 * t86 - t74 * t87 + (-t172 * t205 - t177 * t204 + (-t172 * t73 + t177 * t74) * qJD(3) + t296) * pkin(2), t81, t62, t113, t82, t114, 0, t185 * t171 + (t230 - t295) * t176 + t270, t185 * t176 + (-t134 + t295) * t171 + t245, t146 * t221 + t154 * t186 + t219, -g(1) * t198 - g(2) * t236 + t146 * t224 + t30 * t148 + t186 * t66 + t213 * t65, t8, t5, t46, t9, t47, 0, t115 * t32 + t58 * t162 + t280 * t164 + t271 * t75 + t196, -t115 * t31 - t59 * t162 - t281 * t164 + t271 * t77 + t197, -t280 * t77 - t281 * t75 + t58 * t31 - t59 * t32 + t202, -g(1) * t191 - g(2) * t210 + t23 * t115 + t280 * t19 + t281 * t20 + t271 * t55 + t3 * t59 + t4 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t188 + t272, t73 * t154 + t194, 0, 0, t81, t62, t113, t82, t114, 0, t189 * t171 + (-t193 + t230) * t176 + t270, t189 * t176 + (t193 - t134) * t171 + t245, pkin(8) * t221 - t154 * t225 + t219, -t30 * pkin(3) + pkin(8) * t224 - g(1) * t228 - g(2) * t257 - t225 * t66 - t65 * t74, t8, t5, t46, t9, t47, 0, -t147 * t32 + t69 * t162 + t164 * t278 + t212 * t75 + t196, t147 * t31 - t70 * t162 - t164 * t279 + t212 * t77 + t197, -t278 * t77 - t279 * t75 + t69 * t31 - t70 * t32 + t202, -g(1) * t201 - g(2) * t222 - t23 * t147 + t278 * t19 + t20 * t279 + t212 * t55 + t3 * t70 + t4 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, t254 * t152, t263, t238, t258, qJDD(4), t171 * t190 - t290, g(3) * t171 + t176 * t190, 0, 0, t283, t35, t17, -t283, -t207, t162, -t24 * t164 + (t162 * t175 - t164 * t248 - t265 * t75) * pkin(4) + t183, t25 * t164 + (-t162 * t170 - t164 * t247 - t265 * t77) * pkin(4) + t184, (t20 + t24) * t77 + (-t19 + t25) * t75 + (-t170 * t32 + t175 * t31 + (t170 * t77 - t175 * t75) * qJD(5)) * pkin(4), -t19 * t24 - t20 * t25 + (-t290 + t170 * t3 + t175 * t4 + (-t170 * t19 + t175 * t20) * qJD(5) + (-t154 * t55 + t297) * t171) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, t35, t17, -t283, -t207, t162, t20 * t164 + t183, t19 * t164 + t184, 0, 0;];
tau_reg = t1;
