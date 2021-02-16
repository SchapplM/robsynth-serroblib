% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:01
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:01:02
% EndTime: 2021-01-16 00:01:09
% DurationCPUTime: 1.92s
% Computational Cost: add. (3288->303), mult. (4938->357), div. (0->0), fcn. (3237->12), ass. (0->202)
t157 = qJ(1) + qJ(2);
t143 = sin(t157);
t145 = cos(t157);
t192 = g(2) * t145 + g(3) * t143;
t160 = sin(qJ(2));
t229 = qJD(2) * t160;
t140 = pkin(1) * t229;
t163 = cos(qJ(2));
t269 = pkin(1) * t163;
t232 = -qJD(1) * t140 + qJDD(1) * t269;
t150 = qJDD(1) + qJDD(2);
t268 = pkin(2) * t150;
t75 = -t232 - t268;
t286 = t75 + t192;
t271 = cos(qJ(4));
t158 = sin(qJ(4));
t162 = cos(qJ(3));
t153 = qJD(1) + qJD(2);
t270 = pkin(1) * t160;
t221 = qJD(1) * t270;
t102 = pkin(7) * t153 + t221;
t208 = pkin(8) * t153 + t102;
t63 = t208 * t162;
t55 = t158 * t63;
t159 = sin(qJ(3));
t62 = t208 * t159;
t58 = qJD(3) * pkin(3) - t62;
t207 = t271 * t58 - t55;
t215 = t271 * t159;
t94 = t158 * t162 + t215;
t72 = t94 * t153;
t246 = t72 * qJ(5);
t19 = t207 - t246;
t210 = qJD(4) * t271;
t214 = t271 * t162;
t285 = -qJD(3) * t214 - t162 * t210;
t147 = t162 * pkin(3);
t284 = -pkin(2) - t147;
t273 = -pkin(7) - pkin(8);
t216 = qJD(3) * t273;
t100 = t162 * t216;
t230 = qJD(1) * t163;
t219 = pkin(1) * t230;
t118 = t273 * t159;
t146 = t162 * pkin(8);
t119 = pkin(7) * t162 + t146;
t253 = t158 * t118 + t271 * t119;
t99 = t159 * t216;
t283 = -t253 * qJD(4) + t271 * t100 - t158 * t99 + t94 * t219;
t239 = t158 * t159;
t184 = t214 - t239;
t225 = qJD(4) * t158;
t282 = -t158 * t100 - t118 * t210 + t119 * t225 + t184 * t219 - t271 * t99;
t226 = qJD(3) * t162;
t211 = t153 * t226;
t238 = t159 * t150;
t224 = qJDD(1) * t160;
t228 = qJD(2) * t163;
t76 = t150 * pkin(7) + (qJD(1) * t228 + t224) * pkin(1);
t29 = -t102 * t226 + qJDD(3) * pkin(3) - t159 * t76 + (-t211 - t238) * pkin(8);
t227 = qJD(3) * t159;
t212 = t153 * t227;
t237 = t162 * t150;
t30 = -t102 * t227 + t162 * t76 + (-t212 + t237) * pkin(8);
t281 = -t158 * t30 + t271 * t29;
t149 = qJDD(3) + qJDD(4);
t152 = qJD(3) + qJD(4);
t266 = pkin(3) * t152;
t280 = -t158 * pkin(3) * t149 - t210 * t266;
t279 = g(2) * t143 - g(3) * t145;
t139 = pkin(3) * t227;
t278 = t139 - t221;
t141 = t149 * pkin(4);
t189 = t152 * t239;
t196 = -t150 * t215 + t285 * t153 - t158 * t237;
t25 = t153 * t189 + t196;
t250 = t25 * qJ(5);
t277 = t141 + t250;
t165 = qJD(3) ^ 2;
t276 = pkin(7) * t165 - t268;
t156 = qJ(3) + qJ(4);
t142 = sin(t156);
t243 = t142 * t145;
t244 = t142 * t143;
t144 = cos(t156);
t264 = g(1) * t144;
t275 = g(2) * t244 - g(3) * t243 - t264;
t274 = t72 ^ 2;
t272 = pkin(4) * t184;
t267 = pkin(2) * t153;
t218 = t153 * t239;
t70 = -t153 * t214 + t218;
t261 = t72 * t70;
t134 = pkin(7) + t270;
t260 = -pkin(8) - t134;
t47 = t152 * t94;
t256 = -t47 * qJ(5) + qJD(5) * t184;
t259 = t256 - t282;
t46 = t189 + t285;
t188 = t46 * qJ(5) - t94 * qJD(5);
t258 = t188 + t283;
t18 = pkin(4) * t152 + t19;
t257 = t18 - t19;
t255 = -t271 * t62 - t55;
t90 = t260 * t159;
t91 = t134 * t162 + t146;
t254 = t158 * t90 + t271 * t91;
t252 = qJ(5) * t94;
t190 = -t150 * t214 + t158 * t238;
t26 = t47 * t153 + t190;
t249 = t26 * qJ(5);
t248 = t70 * qJ(5);
t247 = t70 * t152;
t242 = t143 * t144;
t241 = t144 * t145;
t240 = t153 * t159;
t209 = t70 * pkin(4) + qJD(5);
t74 = t153 * t284 - t219;
t40 = t209 + t74;
t236 = qJD(5) + t40;
t235 = g(2) * t243 + g(3) * t244;
t233 = pkin(4) * t144 + t147;
t154 = t159 ^ 2;
t231 = -t162 ^ 2 + t154;
t223 = pkin(3) * t240;
t220 = pkin(1) * t228;
t57 = t271 * t63;
t136 = -pkin(2) - t269;
t213 = t153 * t229;
t43 = pkin(4) * t47 + t139;
t206 = t158 * t62 - t57;
t205 = -t158 * t91 + t271 * t90;
t204 = qJD(3) * t260;
t203 = t271 * t118 - t119 * t158;
t151 = -qJ(5) + t273;
t98 = -pkin(2) - t233;
t202 = -t143 * t98 + t145 * t151;
t201 = -t143 * t151 - t145 * t98;
t44 = pkin(3) * t212 + t150 * t284 - t232;
t15 = t26 * pkin(4) + qJDD(5) + t44;
t200 = t15 * t94 - t40 * t46 + t235;
t199 = t44 * t94 - t74 * t46 + t235;
t103 = -t219 - t267;
t198 = t103 * t226 + t286 * t159;
t197 = t153 * t221;
t113 = t136 - t147;
t193 = t43 - t221;
t186 = -t158 * t58 - t57;
t20 = -t186 - t248;
t174 = t186 * qJD(4) + t281;
t3 = -t72 * qJD(5) + t174 + t277;
t172 = t158 * t29 + t58 * t210 - t63 * t225 + t271 * t30;
t4 = -t70 * qJD(5) + t172 - t249;
t187 = t18 * t46 + t184 * t4 - t20 * t47 - t3 * t94 - t279;
t59 = t159 * t204 + t162 * t220;
t60 = -t159 * t220 + t162 * t204;
t185 = t158 * t60 + t90 * t210 - t91 * t225 + t271 * t59;
t182 = -t103 * t153 + t279 - t76;
t181 = -t192 * t144 - t15 * t184 + t40 * t47;
t180 = -g(2) * t241 - g(3) * t242 - t184 * t44 + t74 * t47;
t179 = -t192 + t197;
t178 = pkin(1) * t213 + t134 * t165 + t136 * t150;
t177 = -t152 * t218 - t196;
t176 = -pkin(7) * qJDD(3) + (t219 - t267) * qJD(3);
t175 = -qJDD(3) * t134 + (t136 * t153 - t220) * qJD(3);
t173 = -t254 * qJD(4) - t158 * t59 + t271 * t60;
t170 = g(1) * t142 + g(2) * t242 - g(3) * t241 - t172;
t169 = t174 + t275;
t168 = t74 * t70 + t170;
t167 = -t74 * t72 + t169;
t166 = t236 * t70 + t170 + t249;
t164 = cos(qJ(1));
t161 = sin(qJ(1));
t148 = t153 ^ 2;
t135 = t271 * pkin(3) + pkin(4);
t112 = qJDD(3) * t162 - t159 * t165;
t111 = qJDD(3) * t159 + t162 * t165;
t101 = t140 + t139;
t89 = t184 * qJ(5);
t83 = t103 * t227;
t77 = t150 * t154 + 0.2e1 * t159 * t211;
t69 = t284 - t272;
t68 = t70 ^ 2;
t61 = t113 - t272;
t49 = -0.2e1 * t231 * t153 * qJD(3) + 0.2e1 * t159 * t237;
t48 = pkin(4) * t72 + t223;
t42 = t89 + t253;
t41 = t203 - t252;
t39 = t140 + t43;
t34 = t89 + t254;
t33 = t205 - t252;
t32 = t149 * t184 - t152 * t47;
t31 = t149 * t94 - t152 * t46;
t28 = -t68 + t274;
t22 = -t246 + t255;
t21 = t206 + t248;
t16 = t177 + t247;
t8 = -t25 * t94 - t46 * t72;
t7 = t173 + t188;
t6 = t185 + t256;
t5 = -t184 * t25 - t26 * t94 + t46 * t70 - t47 * t72;
t1 = [qJDD(1), -g(2) * t164 - g(3) * t161, g(2) * t161 - g(3) * t164, t150, (t150 * t163 - t213) * pkin(1) - t192 + t232, ((-qJDD(1) - t150) * t160 + (-qJD(1) - t153) * t228) * pkin(1) + t279, t77, t49, t111, t112, 0, t83 + t175 * t159 + (-t178 - t286) * t162, t159 * t178 + t162 * t175 + t198, t8, t5, t31, t32, 0, t101 * t70 + t113 * t26 + t205 * t149 + t173 * t152 + t180, t101 * t72 - t113 * t25 - t254 * t149 - t185 * t152 + t199, t149 * t33 + t152 * t7 + t26 * t61 + t39 * t70 + t181, -t149 * t34 - t152 * t6 - t25 * t61 + t39 * t72 + t200, t25 * t33 - t26 * t34 - t6 * t70 - t7 * t72 + t187, t4 * t34 + t20 * t6 + t3 * t33 + t18 * t7 + t15 * t61 + t40 * t39 - g(2) * (pkin(1) * t164 + t201) - g(3) * (pkin(1) * t161 + t202); 0, 0, 0, t150, t179 + t232, (-t224 + (-qJD(2) + t153) * t230) * pkin(1) + t279, t77, t49, t111, t112, 0, t83 + t176 * t159 + (t179 - t75 - t276) * t162, t176 * t162 + (-t197 + t276) * t159 + t198, t8, t5, t31, t32, 0, t203 * t149 + t283 * t152 + t26 * t284 + t278 * t70 + t180, -t253 * t149 + t282 * t152 - t25 * t284 + t278 * t72 + t199, t149 * t41 + t258 * t152 + t193 * t70 + t26 * t69 + t181, -t149 * t42 - t259 * t152 + t193 * t72 - t25 * t69 + t200, t25 * t41 - t258 * t72 - t259 * t70 - t26 * t42 + t187, -g(2) * t201 - g(3) * t202 + t15 * t69 + t258 * t18 + t193 * t40 + t259 * t20 + t3 * t41 + t4 * t42; 0, 0, 0, 0, 0, 0, -t159 * t148 * t162, t231 * t148, t238, t237, qJDD(3), -g(1) * t162 + t159 * t182, g(1) * t159 + t162 * t182, t261, t28, t16, -t190, t149, -t206 * t152 + (t271 * t149 - t152 * t225 - t70 * t240) * pkin(3) + t167, t255 * t152 - t72 * t223 + t168 + t280, t135 * t149 - t21 * t152 - t48 * t70 - t236 * t72 + (-t57 + (-t58 - t266) * t158) * qJD(4) + t275 + t277 + t281, t22 * t152 - t48 * t72 + t166 + t280, t135 * t25 + (t20 + t21) * t72 + (-t18 + t22) * t70 + (-t158 * t26 + (t158 * t72 - t271 * t70) * qJD(4)) * pkin(3), t3 * t135 - t20 * t22 - t18 * t21 - t40 * t48 - g(1) * t233 + t279 * (pkin(3) * t159 + pkin(4) * t142) + (t4 * t158 + (-t158 * t18 + t271 * t20) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t261, t28, t16, -t190, t149, -t152 * t186 + t167, t207 * t152 + t168, t250 + t20 * t152 + 0.2e1 * t141 + (-t209 - t40) * t72 + t169, -t274 * pkin(4) + t19 * t152 + t166, t25 * pkin(4) - t257 * t70, t257 * t20 + (t142 * t279 - t40 * t72 - t264 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t152 + t26, t177 - t247, -t68 - t274, t18 * t72 + t20 * t70 + t15 + t192;];
tau_reg = t1;
