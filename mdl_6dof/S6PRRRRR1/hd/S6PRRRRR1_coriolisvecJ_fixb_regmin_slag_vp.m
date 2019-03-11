% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:40:01
% EndTime: 2019-03-09 00:40:09
% DurationCPUTime: 3.68s
% Computational Cost: add. (5998->337), mult. (15282->484), div. (0->0), fcn. (12258->12), ass. (0->212)
t170 = cos(qJ(6));
t238 = qJD(6) * t170;
t167 = sin(qJ(4));
t168 = sin(qJ(3));
t245 = qJD(2) * t168;
t226 = t167 * t245;
t172 = cos(qJ(4));
t173 = cos(qJ(3));
t243 = qJD(2) * t173;
t227 = t172 * t243;
t128 = -t226 + t227;
t129 = -t167 * t243 - t172 * t245;
t166 = sin(qJ(5));
t171 = cos(qJ(5));
t90 = t128 * t171 + t129 * t166;
t299 = t170 * t90;
t304 = t238 - t299;
t156 = -pkin(3) * t173 - pkin(2);
t174 = cos(qJ(2));
t163 = sin(pkin(6));
t248 = qJD(1) * t163;
t231 = t174 * t248;
t122 = qJD(2) * t156 - t231;
t98 = -pkin(4) * t128 + t122;
t280 = t90 * t98;
t169 = sin(qJ(2));
t232 = t169 * t248;
t288 = pkin(8) + pkin(9);
t215 = qJD(2) * t288 + t232;
t164 = cos(pkin(6));
t247 = qJD(1) * t164;
t112 = t168 * t247 + t173 * t215;
t242 = qJD(4) * t167;
t246 = qJD(2) * t163;
t224 = qJD(1) * t246;
t203 = t174 * t224;
t77 = -qJD(3) * t112 - t168 * t203;
t216 = -t112 * t242 + t167 * t77;
t111 = -t168 * t215 + t173 * t247;
t271 = qJD(3) * pkin(3);
t106 = t111 + t271;
t76 = qJD(3) * t111 + t173 * t203;
t292 = (qJD(4) * t106 + t76) * t172;
t135 = t167 * t173 + t168 * t172;
t160 = qJD(3) + qJD(4);
t108 = t160 * t135;
t96 = t108 * qJD(2);
t18 = -pkin(10) * t96 + t216 + t292;
t105 = t172 * t112;
t196 = -t106 * t167 - t105;
t217 = -t167 * t76 + t172 * t77;
t181 = qJD(4) * t196 + t217;
t237 = qJD(2) * qJD(3);
t223 = t173 * t237;
t95 = qJD(4) * t227 - t160 * t226 + t172 * t223;
t19 = -pkin(10) * t95 + t181;
t241 = qJD(5) * t166;
t285 = pkin(10) * t128;
t60 = -t196 + t285;
t219 = t166 * t19 - t241 * t60;
t123 = t129 * pkin(10);
t103 = t167 * t112;
t213 = t106 * t172 - t103;
t59 = t123 + t213;
t52 = pkin(4) * t160 + t59;
t3 = (qJD(5) * t52 + t18) * t171 + t219;
t303 = -t280 - t3;
t233 = qJD(3) * t288;
t136 = t168 * t233;
t137 = t173 * t233;
t141 = t288 * t168;
t142 = t288 * t173;
t194 = t141 * t167 - t142 * t172;
t302 = qJD(4) * t194 + t135 * t231 + t167 * t136 - t137 * t172;
t134 = t167 * t168 - t172 * t173;
t259 = t141 * t172;
t301 = qJD(4) * t259 - t134 * t231 + t136 * t172 + t137 * t167 + t142 * t242;
t251 = -qJD(6) + t90;
t300 = qJD(6) + t251;
t165 = sin(qJ(6));
t159 = qJD(5) + t160;
t195 = t128 * t166 - t129 * t171;
t239 = qJD(6) * t165;
t240 = qJD(5) * t171;
t41 = t128 * t240 + t129 * t241 - t166 * t96 + t171 * t95;
t31 = t159 * t238 + t170 * t41 - t195 * t239;
t27 = t31 * t165;
t73 = t159 * t165 + t170 * t195;
t7 = t304 * t73 + t27;
t42 = qJD(5) * t195 + t166 * t95 + t171 * t96;
t274 = t165 * t42 - t238 * t251;
t6 = -t195 * t73 + t251 * t299 + t274;
t296 = t165 * t251;
t40 = t170 * t42;
t71 = -t159 * t170 + t165 * t195;
t5 = t195 * t71 - t251 * t296 + t40;
t28 = t31 * t170;
t32 = qJD(6) * t73 + t165 * t41;
t1 = -t165 * t32 + t296 * t73 - t304 * t71 + t28;
t267 = t166 * t60;
t24 = t171 * t52 - t267;
t22 = -pkin(5) * t159 - t24;
t284 = t22 * t90;
t279 = t195 * t90;
t107 = t160 * t134;
t298 = -pkin(10) * t107 - t302;
t297 = -pkin(10) * t108 - t301;
t281 = t195 * t98;
t220 = t166 * t18 - t171 * t19;
t263 = t171 * t60;
t25 = t166 * t52 + t263;
t4 = qJD(5) * t25 + t220;
t295 = -t4 - t281;
t186 = -t168 * t271 + t232;
t37 = t195 ^ 2 - t90 ^ 2;
t63 = pkin(5) * t195 - pkin(11) * t90;
t35 = -t159 * t90 + t41;
t282 = t251 * t195;
t200 = pkin(4) * t108 - t186;
t23 = pkin(11) * t159 + t25;
t43 = -pkin(5) * t90 - pkin(11) * t195 + t98;
t198 = t165 * t23 - t170 * t43;
t234 = t195 * t198 + t22 * t239;
t9 = t165 * t43 + t170 * t23;
t206 = t165 * t4 + t195 * t9 + t22 * t238;
t101 = t134 * t171 + t135 * t166;
t102 = -t134 * t166 + t135 * t171;
t79 = -pkin(10) * t135 - t142 * t167 - t259;
t80 = -pkin(10) * t134 - t194;
t55 = t166 * t80 - t171 * t79;
t278 = qJD(5) * t55 + t166 * t298 - t171 * t297;
t47 = -qJD(5) * t101 - t107 * t171 - t108 * t166;
t117 = pkin(4) * t134 + t156;
t54 = pkin(5) * t101 - pkin(11) * t102 + t117;
t56 = t166 * t79 + t171 * t80;
t289 = -(qJD(6) * t43 + t3) * t101 + t22 * t47 + t4 * t102 - (-qJD(6) * t54 + t278) * t251 - t56 * t42;
t36 = t159 * t195 - t42;
t286 = pkin(4) * t129;
t283 = t54 * t42;
t277 = qJD(5) * t56 + t166 * t297 + t171 * t298;
t155 = pkin(3) * t172 + pkin(4);
t255 = t166 * t167;
t207 = -t111 * t167 - t105;
t61 = t207 - t285;
t250 = t111 * t172 - t103;
t62 = t123 + t250;
t275 = t166 * t61 + t171 * t62 - t155 * t240 - (-t167 * t241 + (t171 * t172 - t255) * qJD(4)) * pkin(3);
t254 = t167 * t171;
t273 = -t166 * t62 + t171 * t61 + t155 * t241 + (t167 * t240 + (t166 * t172 + t254) * qJD(4)) * pkin(3);
t272 = qJD(2) * pkin(2);
t270 = t102 * t22;
t266 = t170 * t73;
t261 = t122 * t129;
t260 = t129 * t128;
t258 = t163 * t169;
t257 = t163 * t174;
t176 = qJD(2) ^ 2;
t256 = t163 * t176;
t175 = qJD(3) ^ 2;
t253 = t175 * t168;
t252 = t175 * t173;
t157 = pkin(3) * t245;
t124 = qJD(3) * t157 + t169 * t224;
t249 = t168 ^ 2 - t173 ^ 2;
t244 = qJD(2) * t169;
t235 = t169 * t256;
t229 = t163 * t244;
t228 = t174 * t246;
t225 = -pkin(4) * t159 - t52;
t222 = -pkin(3) * t160 - t106;
t121 = pkin(3) * t254 + t155 * t166 + pkin(11);
t49 = -t286 + t63;
t209 = qJD(6) * t121 + t157 + t49;
t153 = pkin(4) * t166 + pkin(11);
t208 = qJD(6) * t153 + t49;
t205 = t168 * t228;
t204 = t173 * t228;
t70 = pkin(4) * t96 + t124;
t29 = t166 * t59 + t263;
t202 = pkin(4) * t241 - t29;
t48 = qJD(5) * t102 - t107 * t166 + t108 * t171;
t201 = pkin(5) * t48 - pkin(11) * t47 + t200;
t125 = t164 * t173 - t168 * t258;
t126 = t164 * t168 + t173 * t258;
t82 = t125 * t172 - t126 * t167;
t83 = t125 * t167 + t126 * t172;
t57 = t166 * t83 - t171 * t82;
t58 = t166 * t82 + t171 * t83;
t193 = -t122 * t128 - t216;
t191 = -t239 * t251 - t40;
t190 = -t165 * t58 - t170 * t257;
t189 = t165 * t257 - t170 * t58;
t188 = t272 * qJD(2);
t185 = -t121 * t42 - t251 * t275 - t284;
t184 = -0.2e1 * qJD(3) * t272;
t30 = t171 * t59 - t267;
t179 = -t153 * t42 - t284 - (-pkin(4) * t240 + t30) * t251;
t154 = -pkin(4) * t171 - pkin(5);
t120 = pkin(3) * t255 - t155 * t171 - pkin(5);
t115 = t157 - t286;
t110 = -qJD(3) * t126 - t205;
t109 = qJD(3) * t125 + t204;
t74 = -t128 ^ 2 + t129 ^ 2;
t67 = (-qJD(2) * t135 - t129) * t160;
t66 = -t128 * t160 + t95;
t45 = -qJD(4) * t83 - t109 * t167 + t110 * t172;
t44 = qJD(4) * t82 + t109 * t172 + t110 * t167;
t13 = pkin(5) * t42 - pkin(11) * t41 + t70;
t12 = t170 * t13;
t11 = qJD(5) * t58 + t166 * t44 - t171 * t45;
t10 = -qJD(5) * t57 + t166 * t45 + t171 * t44;
t2 = [0, 0, -t235, -t174 * t256, 0, 0, 0, 0, 0, -t173 * t235 + (t110 - t205) * qJD(3), t168 * t235 + (-t109 - t204) * qJD(3), 0, 0, 0, 0, 0, t160 * t45 + (-t128 * t244 - t174 * t96) * t163, -t160 * t44 + (-t129 * t244 - t174 * t95) * t163, 0, 0, 0, 0, 0, -t11 * t159 + (-t174 * t42 - t244 * t90) * t163, -t10 * t159 + (-t174 * t41 + t195 * t244) * t163, 0, 0, 0, 0, 0 -(qJD(6) * t189 - t10 * t165 + t170 * t229) * t251 + t190 * t42 + t11 * t71 + t57 * t32 (qJD(6) * t190 + t10 * t170 + t165 * t229) * t251 + t189 * t42 + t11 * t73 + t57 * t31; 0, 0, 0, 0, 0.2e1 * t168 * t223, -0.2e1 * t249 * t237, t252, -t253, 0, -pkin(8) * t252 + t168 * t184, pkin(8) * t253 + t173 * t184, t107 * t129 + t135 * t95, -t107 * t128 + t108 * t129 - t134 * t95 - t135 * t96, -t107 * t160, -t108 * t160, 0, t122 * t108 + t124 * t134 + t128 * t186 + t156 * t96 + t160 * t302, -t122 * t107 + t124 * t135 + t129 * t186 + t156 * t95 + t160 * t301, t102 * t41 + t195 * t47, -t101 * t41 - t102 * t42 - t195 * t48 + t47 * t90, t47 * t159, -t48 * t159, 0, t101 * t70 + t117 * t42 - t159 * t277 - t200 * t90 + t48 * t98, t102 * t70 + t117 * t41 + t159 * t278 + t195 * t200 + t47 * t98, t47 * t266 + (-t239 * t73 + t28) * t102 (-t165 * t73 - t170 * t71) * t47 + (-t27 - t170 * t32 + (t165 * t71 - t266) * qJD(6)) * t102, -t170 * t251 * t47 + t101 * t31 - t102 * t191 + t48 * t73, -t101 * t32 - t102 * t274 + t296 * t47 - t48 * t71, t101 * t42 - t251 * t48, t12 * t101 + t55 * t32 - t198 * t48 + t277 * t71 + (t283 - t201 * t251 + (-t101 * t23 + t251 * t56 + t270) * qJD(6)) * t170 + t289 * t165, t55 * t31 - t9 * t48 + t277 * t73 + (-t283 - (-qJD(6) * t23 + t13) * t101 - qJD(6) * t270 - (qJD(6) * t56 - t201) * t251) * t165 + t289 * t170; 0, 0, 0, 0, -t168 * t176 * t173, t249 * t176, 0, 0, 0, t168 * t188, t173 * t188, t260, t74, t66, t67, 0, -t207 * t160 + t128 * t157 + t261 + (t167 * t222 - t105) * qJD(4) + t217, t250 * t160 + t129 * t157 + (qJD(4) * t222 - t76) * t172 + t193, -t279, t37, t35, t36, 0, t115 * t90 - t159 * t273 + t295, -t115 * t195 + t159 * t275 + t303, t7, t1, t6, t5, t282, t120 * t32 + t273 * t71 + (t209 * t251 - t4) * t170 + t185 * t165 + t234, t120 * t31 + t170 * t185 - t209 * t296 + t273 * t73 + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t74, t66, t67, 0, -t160 * t196 + t181 + t261, t160 * t213 + t193 - t292, -t279, t37, t35, t36, 0, -t90 * t286 + t159 * t29 - t281 + (t166 * t225 - t263) * qJD(5) - t220, t195 * t286 + t159 * t30 - t280 + (qJD(5) * t225 - t18) * t171 - t219, t7, t1, t6, t5, t282, t154 * t32 + t202 * t71 + (t208 * t251 - t4) * t170 + t179 * t165 + t234, t154 * t31 + t170 * t179 + t202 * t73 - t208 * t296 + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t279, t37, t35, t36, 0, t159 * t25 + t295, t159 * t24 + t303, t7, t1, t6, t5, t282, -pkin(5) * t32 - t4 * t170 + (-t165 * t24 + t170 * t63) * t251 - t25 * t71 - t165 * t284 - t274 * pkin(11) + t234, -pkin(5) * t31 - (t165 * t63 + t170 * t24) * t251 - t25 * t73 - t22 * t299 + t191 * pkin(11) + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, -t251 * t71 + t31, -t251 * t73 - t32, t42, -t165 * t3 - t22 * t73 - t300 * t9 + t12, -t13 * t165 - t170 * t3 + t198 * t300 + t22 * t71;];
tauc_reg  = t2;
