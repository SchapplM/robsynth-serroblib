% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRPPRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPPRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:09
% EndTime: 2019-03-08 19:23:15
% DurationCPUTime: 3.81s
% Computational Cost: add. (4214->430), mult. (8278->584), div. (0->0), fcn. (6318->12), ass. (0->209)
t144 = sin(pkin(11));
t147 = cos(pkin(11));
t156 = cos(qJ(2));
t146 = sin(pkin(6));
t241 = qJD(1) * t146;
t211 = t156 * t241;
t153 = sin(qJ(2));
t212 = t153 * t241;
t65 = t144 * t211 - t147 * t212;
t287 = t144 * qJD(3) - t65;
t157 = -pkin(2) - pkin(3);
t117 = qJD(2) * t212;
t229 = qJDD(1) * t146;
t120 = t156 * t229;
t208 = qJDD(3) + t117 - t120;
t62 = t157 * qJDD(2) + t208;
t119 = t153 * t229;
t228 = qJDD(2) * qJ(3);
t63 = t228 + t119 + (qJD(3) + t211) * qJD(2);
t20 = -t144 * t63 + t147 * t62;
t18 = qJDD(2) * pkin(4) - t20;
t155 = cos(qJ(5));
t137 = t155 * qJDD(2);
t152 = sin(qJ(5));
t230 = qJD(2) * qJD(5);
t209 = t152 * t230;
t279 = -t209 + t137;
t210 = t155 * t230;
t225 = t152 * qJDD(2);
t283 = t210 + t225;
t14 = t279 * pkin(5) + pkin(9) * t283 + t18;
t151 = sin(qJ(6));
t154 = cos(qJ(6));
t149 = cos(pkin(6));
t126 = -t149 * qJD(1) + qJD(4);
t104 = qJD(2) * qJ(3) + t212;
t189 = qJD(3) - t211;
t93 = t157 * qJD(2) + t189;
t54 = t147 * t104 + t144 * t93;
t50 = -qJD(2) * pkin(8) + t54;
t35 = t152 * t126 + t155 * t50;
t27 = qJD(5) * pkin(9) + t35;
t277 = t155 * pkin(5);
t196 = t152 * pkin(9) + t277;
t53 = -t144 * t104 + t147 * t93;
t49 = qJD(2) * pkin(4) - t53;
t36 = t196 * qJD(2) + t49;
t187 = t151 * t27 - t154 * t36;
t122 = -t149 * qJDD(1) + qJDD(4);
t21 = t144 * t62 + t147 * t63;
t19 = -qJDD(2) * pkin(8) + t21;
t238 = qJD(5) * t155;
t224 = -t152 * t122 - t126 * t238 - t155 * t19;
t239 = qJD(5) * t152;
t7 = -t50 * t239 - t224;
t5 = qJDD(5) * pkin(9) + t7;
t1 = -t187 * qJD(6) + t151 * t14 + t154 * t5;
t231 = t155 * qJD(2);
t130 = qJD(6) + t231;
t286 = t187 * t130 + t1;
t106 = t147 * qJ(3) + t144 * t157;
t101 = -pkin(8) + t106;
t195 = -pkin(5) * t152 + pkin(9) * t155;
t285 = -qJD(6) * t101 * t155 + t195 * qJD(5) + t287;
t145 = sin(pkin(10));
t148 = cos(pkin(10));
t248 = t149 * t156;
t89 = t145 * t248 + t148 * t153;
t249 = t149 * t153;
t90 = -t145 * t249 + t148 * t156;
t204 = t90 * t144 - t89 * t147;
t87 = t145 * t153 - t148 * t248;
t88 = t145 * t156 + t148 * t249;
t205 = t88 * t144 - t87 * t147;
t76 = (t144 * t153 + t147 * t156) * t146;
t172 = g(1) * t204 + g(2) * t205 + g(3) * t76;
t284 = t287 * qJD(2) - t172;
t234 = t147 * qJD(3);
t68 = qJD(1) * t76;
t201 = -t68 + t234;
t282 = qJD(2) * t201;
t10 = t151 * t36 + t154 * t27;
t2 = -qJD(6) * t10 + t154 * t14 - t151 * t5;
t281 = -t10 * t130 - t2;
t105 = -t144 * qJ(3) + t147 * t157;
t100 = pkin(4) - t105;
t158 = qJD(5) ^ 2;
t280 = qJDD(2) * t100 - t101 * t158 + t18 + t284;
t251 = t146 * t156;
t253 = t146 * t153;
t278 = t144 * t251 - t147 * t253;
t108 = t155 * t122;
t8 = -t35 * qJD(5) - t152 * t19 + t108;
t266 = t152 * t50;
t34 = t155 * t126 - t266;
t162 = -(t152 * t35 + t155 * t34) * qJD(5) - t8 * t152 + t7 * t155;
t232 = t154 * qJD(5);
t240 = qJD(2) * t152;
t98 = t151 * t240 + t232;
t233 = t151 * qJD(5);
t99 = t154 * t240 - t233;
t275 = t98 * t99;
t274 = t149 ^ 2 * qJDD(1) - g(3);
t69 = t100 + t196;
t165 = qJD(6) * t69 - t101 * t239 + t155 * t234;
t246 = t154 * t155;
t273 = t151 * t285 + t165 * t154 - t68 * t246;
t247 = t151 * t155;
t272 = -t165 * t151 + t154 * t285 + t68 * t247;
t214 = t144 * t239;
t91 = -t144 * t247 - t154 * t147;
t271 = t91 * qJD(6) - t154 * t214 - (t144 * t151 + t147 * t246) * qJD(2);
t92 = t144 * t246 - t151 * t147;
t270 = -t92 * qJD(6) + t151 * t214 - (t144 * t154 - t147 * t247) * qJD(2);
t268 = t151 * t98;
t267 = t151 * t99;
t265 = t154 * t98;
t264 = t154 * t99;
t256 = qJD(6) * t98;
t51 = t151 * qJDD(5) - t154 * t283 + t256;
t263 = t51 * t151;
t236 = qJD(6) * t154;
t213 = t152 * t236;
t52 = t151 * t225 + t154 * qJDD(5) - qJD(6) * t233 + (t155 * t233 + t213) * qJD(2);
t262 = t52 * t154;
t261 = t98 * t130;
t260 = t99 * t130;
t259 = qJD(2) * t49;
t258 = qJD(5) * t98;
t257 = qJD(5) * t99;
t255 = qJDD(2) * pkin(2);
t254 = t146 * t152;
t252 = t146 * t155;
t159 = qJD(2) ^ 2;
t250 = t147 * t159;
t245 = pkin(2) * t251 + qJ(3) * t253;
t141 = t152 ^ 2;
t142 = t155 ^ 2;
t244 = t141 - t142;
t243 = t141 + t142;
t242 = t158 + t159;
t237 = qJD(6) * t151;
t227 = qJDD(5) * t152;
t226 = qJDD(5) * t155;
t223 = t99 * t238;
t220 = t152 * t159 * t155;
t26 = -qJD(5) * pkin(5) - t34;
t219 = t26 * t238;
t218 = pkin(3) * t251 + t245;
t217 = t147 * t240;
t216 = t130 * t233;
t215 = t130 * t232;
t207 = -t87 * pkin(2) + t88 * qJ(3);
t206 = -t89 * pkin(2) + t90 * qJ(3);
t44 = t87 * t144 + t88 * t147;
t48 = t89 * t144 + t90 * t147;
t200 = qJDD(2) * t243;
t199 = -t87 * pkin(3) + t207;
t198 = -t89 * pkin(3) + t206;
t197 = t155 * t209;
t192 = t10 * t154 + t151 * t187;
t191 = t10 * t151 - t154 * t187;
t188 = t53 * t144 - t54 * t147;
t56 = -t149 * t152 - t155 * t278;
t25 = t76 * t151 + t56 * t154;
t24 = -t56 * t151 + t76 * t154;
t185 = t34 * t152 - t35 * t155;
t184 = t104 * t156 + t153 * (-qJD(2) * pkin(2) + t189);
t55 = t149 * t155 - t152 * t278;
t182 = -g(1) * t89 - g(2) * t87 + g(3) * t251;
t181 = g(3) * t149 + (g(1) * t145 - g(2) * t148) * t146;
t180 = t152 * t52 + t98 * t238;
t179 = t76 * pkin(4) + pkin(8) * t278 + t218;
t96 = -qJDD(6) - t279;
t178 = t130 * t236 - t151 * t96;
t177 = t130 * t237 + t154 * t96;
t176 = t120 - t182;
t175 = -g(1) * (-t145 * t252 - t48 * t152) - g(2) * (t148 * t252 - t44 * t152) + g(3) * t55;
t31 = t148 * t254 + t44 * t155;
t33 = -t145 * t254 + t48 * t155;
t174 = -g(1) * t33 - g(2) * t31 - g(3) * t56;
t173 = -g(1) * t48 - g(2) * t44 + g(3) * t278;
t6 = -qJDD(5) * pkin(5) - t8;
t171 = t175 - t6;
t170 = pkin(4) * t205 - pkin(8) * t44 + t199;
t169 = pkin(4) * t204 - pkin(8) * t48 + t198;
t168 = -qJDD(3) + t176;
t167 = pkin(9) * t96 + t130 * t26;
t166 = g(1) * t90 + g(2) * t88 + g(3) * t253 - t119;
t164 = pkin(9) * qJD(6) * t130 - t171;
t163 = -t191 * qJD(6) + t1 * t154 - t2 * t151;
t161 = -qJDD(5) * t101 + (-qJD(2) * t100 - t201 - t49) * qJD(5);
t112 = -t158 * t152 + t226;
t111 = -t158 * t155 - t227;
t102 = t195 * qJD(2);
t86 = (qJDD(2) * t153 + t156 * t159) * t146;
t85 = (qJDD(2) * t156 - t153 * t159) * t146;
t67 = qJD(2) * t76;
t66 = t278 * qJD(2);
t64 = t208 - t255;
t38 = t101 * t246 + t151 * t69;
t37 = -t101 * t247 + t154 * t69;
t23 = -t55 * qJD(5) + t67 * t155;
t22 = -t149 * t239 + t67 * t152 - t238 * t278;
t16 = t151 * t102 + t154 * t34;
t15 = t154 * t102 - t151 * t34;
t4 = t24 * qJD(6) + t66 * t151 + t23 * t154;
t3 = -t25 * qJD(6) - t23 * t151 + t66 * t154;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, t85, -t86, 0 (t153 ^ 2 + t156 ^ 2) * t146 ^ 2 * qJDD(1) + t274, 0, 0, 0, 0, 0, 0, t85, 0, t86 (qJD(2) * t184 + t153 * t63 - t156 * t64) * t146 + t274, 0, 0, 0, 0, 0, 0, t66 * qJD(2) + t76 * qJDD(2), t67 * qJD(2) - qJDD(2) * t278, 0, -t122 * t149 - t20 * t76 - t21 * t278 - t53 * t66 + t54 * t67 - g(3), 0, 0, 0, 0, 0, 0, t76 * t137 - t22 * qJD(5) - t55 * qJDD(5) + (t155 * t66 - t239 * t76) * qJD(2), -t76 * t225 - t23 * qJD(5) - t56 * qJDD(5) + (-t152 * t66 - t238 * t76) * qJD(2) (-t152 * t55 - t155 * t56) * qJDD(2) + (-t152 * t22 - t155 * t23 + (t152 * t56 - t155 * t55) * qJD(5)) * qJD(2), t18 * t76 - t34 * t22 + t35 * t23 + t49 * t66 - t8 * t55 + t7 * t56 - g(3), 0, 0, 0, 0, 0, 0, t3 * t130 - t22 * t98 - t24 * t96 - t55 * t52, -t4 * t130 - t22 * t99 + t25 * t96 + t55 * t51, -t24 * t51 + t25 * t52 + t3 * t99 + t4 * t98, t1 * t25 + t10 * t4 - t187 * t3 + t2 * t24 + t26 * t22 + t6 * t55 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t176, t166, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t168 + 0.2e1 * t255, 0, 0.2e1 * qJD(2) * qJD(3) - t166 + 0.2e1 * t228, -t64 * pkin(2) - g(1) * t206 - g(2) * t207 - g(3) * t245 + t63 * qJ(3) + t104 * qJD(3) - t184 * t241, 0, 0, 0, 0, 0, qJDD(2), -t105 * qJDD(2) - t20 + t284, t106 * qJDD(2) + t173 + t21 + t282, 0, -g(1) * t198 - g(2) * t199 - g(3) * t218 - qJD(3) * t188 + t20 * t105 + t21 * t106 + t53 * t65 - t54 * t68, t141 * qJDD(2) + 0.2e1 * t197, 0.2e1 * t152 * t137 - 0.2e1 * t244 * t230, t111, t142 * qJDD(2) - 0.2e1 * t197, -t112, 0, t161 * t152 + t280 * t155, -t280 * t152 + t161 * t155, -t101 * t200 - t243 * t282 - t162 - t173, t18 * t100 - t49 * t65 - g(1) * t169 - g(2) * t170 - g(3) * t179 + t185 * t68 + (t49 * t144 - t147 * t185) * qJD(3) + t162 * t101, t154 * t223 + (-t51 * t154 - t237 * t99) * t152 (-t265 - t267) * t238 + (t263 - t262 + (-t264 + t268) * qJD(6)) * t152 (t51 - t215) * t155 + (t177 + t257) * t152, t151 * t180 + t213 * t98 (t52 + t216) * t155 + (t178 - t258) * t152, -t130 * t239 - t96 * t155, -t37 * t96 + t272 * t130 - t173 * t151 + (t2 + (-t101 * t98 - t151 * t26) * qJD(5) - t172 * t154) * t155 + (qJD(5) * t187 - t101 * t52 - t6 * t151 - t201 * t98 - t236 * t26) * t152, t38 * t96 - t273 * t130 - t173 * t154 + (-t1 + (-t101 * t99 - t154 * t26) * qJD(5) + t172 * t151) * t155 + (t10 * qJD(5) + t101 * t51 - t6 * t154 - t201 * t99 + t237 * t26) * t152, -t37 * t51 + t38 * t52 + t272 * t99 + t273 * t98 + t191 * t238 + (qJD(6) * t192 + t1 * t151 + t154 * t2 - t172) * t152, t1 * t38 + t2 * t37 + t101 * t219 - g(1) * (t204 * t277 + t169) - g(2) * (t205 * t277 + t170) - g(3) * (t76 * t277 + t179) - t272 * t187 + t273 * t10 + (-pkin(9) * t172 + t6 * t101 + t201 * t26) * t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t159, -t104 * qJD(2) + t117 - t168 - t255, 0, 0, 0, 0, 0, 0, -t147 * qJDD(2) - t144 * t159, t144 * qJDD(2) - t250, 0, qJD(2) * t188 + t21 * t144 + t20 * t147 + t182, 0, 0, 0, 0, 0, 0 (0.2e1 * t209 - t137) * t147 + (-t242 * t155 - t227) * t144 (0.2e1 * t210 + t225) * t147 + (t242 * t152 - t226) * t144, -t144 * t200 + t243 * t250 (qJD(2) * t185 - t18) * t147 + (t162 - t259) * t144 + t182, 0, 0, 0, 0, 0, 0, t270 * t130 - t180 * t144 + t98 * t217 - t91 * t96, t99 * t217 + t92 * t96 + (t152 * t51 - t223) * t144 - t271 * t130, t270 * t99 + t271 * t98 - t91 * t51 + t92 * t52, -t26 * t217 + t1 * t92 + t2 * t91 - t270 * t187 + (t6 * t152 + t219) * t144 + t271 * t10 + t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122 + t181, 0, 0, 0, 0, 0, 0, t112, t111, 0, -t185 * qJD(5) + t7 * t152 + t8 * t155 + t181, 0, 0, 0, 0, 0, 0 (t52 - t216) * t155 + (-t178 - t258) * t152 (-t51 - t215) * t155 + (t177 - t257) * t152 (t265 - t267) * t238 + (t263 + t262 + (-t264 - t268) * qJD(6)) * t152 (qJD(5) * t192 - t6) * t155 + (qJD(5) * t26 + t163) * t152 + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t244 * t159, -t225, t220, -t137, qJDD(5), t108 + (-t19 + t259) * t152 + t175, t49 * t231 + (t34 + t266) * qJD(5) - t174 + t224, 0, 0, -t154 * t260 + t263 (t51 + t261) * t154 + (t52 + t260) * t151 (t130 * t246 - t152 * t99) * qJD(2) + t178, -t151 * t261 + t262 (-t130 * t247 + t152 * t98) * qJD(2) - t177, t130 * t240, pkin(5) * t52 - t15 * t130 + t151 * t167 - t154 * t164 - t187 * t240 + t35 * t98, -pkin(5) * t51 - t10 * t240 + t16 * t130 + t151 * t164 + t154 * t167 + t35 * t99, -t15 * t99 - t16 * t98 + ((-qJD(6) * t99 + t52) * pkin(9) + t286) * t154 + ((t51 - t256) * pkin(9) + t281) * t151 + t174, -t10 * t16 + t187 * t15 - t26 * t35 + t171 * pkin(5) + (t163 + t174) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, -t98 ^ 2 + t99 ^ 2, t51 - t261, -t275, t52 - t260, -t96, t26 * t99 - g(1) * (-t33 * t151 + t154 * t204) - g(2) * (-t31 * t151 + t154 * t205) - g(3) * t24 - t281, -t26 * t98 - g(1) * (-t151 * t204 - t33 * t154) - g(2) * (-t151 * t205 - t31 * t154) + g(3) * t25 - t286, 0, 0;];
tau_reg  = t9;
