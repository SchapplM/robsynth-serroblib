% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% tau_reg [6x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:49:37
% EndTime: 2019-03-08 19:49:45
% DurationCPUTime: 3.78s
% Computational Cost: add. (2646->425), mult. (5662->599), div. (0->0), fcn. (4402->14), ass. (0->212)
t143 = sin(pkin(11));
t146 = cos(pkin(11));
t148 = sin(qJ(6));
t151 = cos(qJ(6));
t286 = -t143 * t148 + t151 * t146;
t284 = t286 * qJD(6);
t152 = cos(qJ(4));
t236 = qJD(2) * t152;
t214 = t143 * t236;
t231 = t146 * qJD(4);
t101 = t214 - t231;
t211 = t146 * t236;
t235 = qJD(4) * t143;
t103 = t211 + t235;
t35 = t151 * t101 + t103 * t148;
t149 = sin(qJ(4));
t237 = qJD(2) * t149;
t132 = qJD(6) + t237;
t180 = t101 * t148 - t103 * t151;
t289 = t132 * t180;
t153 = cos(qJ(2));
t154 = -pkin(2) - pkin(8);
t233 = qJD(4) * t152;
t209 = t154 * t233;
t145 = sin(pkin(6));
t242 = qJD(1) * t145;
t150 = sin(qJ(2));
t250 = t149 * t150;
t189 = pkin(4) * t152 + qJ(5) * t149;
t79 = qJD(4) * t189 - qJD(5) * t152 + qJD(3);
t272 = t143 * t79 + t146 * t209 - (t143 * t153 + t146 * t250) * t242;
t288 = t146 * t79 - (-t143 * t250 + t146 * t153) * t242;
t147 = cos(pkin(6));
t251 = t147 * t149;
t217 = t153 * t242;
t187 = qJD(3) - t217;
t97 = t154 * qJD(2) + t187;
t287 = -qJD(1) * t251 + t152 * t97;
t252 = t145 * t153;
t265 = cos(pkin(10));
t196 = t265 * t153;
t144 = sin(pkin(10));
t255 = t144 * t150;
t87 = -t147 * t196 + t255;
t197 = t265 * t150;
t254 = t144 * t153;
t89 = t147 * t254 + t197;
t162 = g(1) * t89 + g(2) * t87 - g(3) * t252;
t88 = t147 * t197 + t254;
t90 = -t147 * t255 + t196;
t285 = -g(1) * t90 - g(2) * t88;
t107 = t143 * t151 + t146 * t148;
t166 = t107 * qJD(6);
t283 = t132 - qJD(6);
t238 = qJD(2) * t145;
t213 = t150 * t238;
t120 = qJD(1) * t213;
t229 = qJDD(1) * t145;
t203 = t153 * t229;
t175 = qJDD(3) + t120 - t203;
t59 = t154 * qJDD(2) + t175;
t282 = -qJDD(4) * pkin(4) - t152 * t59 + qJDD(5);
t246 = qJDD(1) - g(3);
t253 = t145 * t150;
t281 = -t246 * t253 - t285;
t218 = t150 * t242;
t240 = qJD(2) * qJ(3);
t110 = t218 + t240;
t280 = qJD(4) * (-t110 + t218 - t240) - qJDD(4) * t154;
t222 = t152 * qJDD(2);
t194 = -qJDD(4) * t146 + t143 * t222;
t230 = qJD(2) * qJD(4);
t206 = t149 * t230;
t66 = t143 * t206 - t194;
t245 = t143 * qJDD(4) + t146 * t222;
t67 = t146 * t206 - t245;
t10 = -qJD(6) * t180 - t148 * t67 - t151 * t66;
t277 = pkin(9) * t152;
t276 = pkin(9) + qJ(5);
t228 = qJDD(1) * t147;
t202 = t152 * t228;
t14 = t202 + qJDD(4) * qJ(5) + t149 * t59 + (qJD(5) + t287) * qJD(4);
t188 = pkin(4) * t149 - qJ(5) * t152;
t111 = qJ(3) + t188;
t204 = t150 * t229;
t21 = t204 + t111 * qJDD(2) + (t79 + t217) * qJD(2);
t7 = t146 * t14 + t143 * t21;
t201 = -t143 * t154 + pkin(5);
t221 = pkin(9) * t146 * t149;
t275 = (t201 * t152 + t221) * qJD(4) + t288;
t234 = qJD(4) * t149;
t210 = t143 * t234;
t274 = -pkin(9) * t210 - t272;
t241 = qJD(1) * t152;
t128 = t147 * t241;
t53 = t149 * t97 + t128;
t43 = qJD(4) * qJ(5) + t53;
t65 = qJD(2) * t111 + t218;
t17 = t143 * t65 + t146 * t43;
t273 = -t143 * t209 + t288;
t109 = t189 * qJD(2);
t24 = t143 * t109 + t146 * t287;
t169 = t286 * t149;
t271 = qJD(2) * t169 + t284;
t167 = qJD(2) * t107;
t270 = t149 * t167 + t166;
t269 = t132 * t35;
t220 = -qJD(4) * t128 - t149 * t228 - t97 * t234;
t15 = -t220 + t282;
t268 = t15 * t152;
t249 = t149 * t154;
t64 = t143 * t111 + t146 * t249;
t264 = qJDD(2) * pkin(2);
t205 = t152 * t230;
t223 = t149 * qJDD(2);
t168 = t205 + t223;
t105 = qJDD(6) + t168;
t261 = t105 * t286;
t260 = t105 * t107;
t140 = pkin(11) + qJ(6);
t137 = sin(t140);
t259 = t137 * t149;
t138 = cos(t140);
t258 = t138 * t149;
t256 = t144 * t145;
t156 = qJD(2) ^ 2;
t247 = t153 * t156;
t141 = t149 ^ 2;
t142 = t152 ^ 2;
t244 = t141 - t142;
t155 = qJD(4) ^ 2;
t243 = -t155 - t156;
t239 = qJD(2) * t110;
t227 = qJDD(2) * qJ(3);
t226 = qJDD(2) * t141;
t225 = qJDD(4) * t149;
t6 = -t14 * t143 + t146 * t21;
t2 = pkin(5) * t168 + pkin(9) * t67 + t6;
t5 = pkin(9) * t66 + t7;
t219 = -t148 * t5 + t151 * t2;
t216 = t150 * t241;
t215 = t143 * t237;
t212 = t153 * t238;
t208 = g(3) * (pkin(2) * t252 + qJ(3) * t253);
t207 = qJ(5) * t223;
t200 = pkin(5) * t143 - t154;
t16 = -t143 * t43 + t146 * t65;
t23 = t146 * t109 - t143 * t287;
t198 = t145 * t265;
t195 = -t59 + t239;
t191 = -t143 * t6 + t146 * t7;
t190 = t148 * t2 + t151 * t5;
t11 = pkin(5) * t237 - pkin(9) * t103 + t16;
t13 = -pkin(9) * t101 + t17;
t3 = t11 * t151 - t13 * t148;
t4 = t11 * t148 + t13 * t151;
t186 = t143 * t17 + t146 * t16;
t185 = -t143 * t16 + t146 * t17;
t99 = t146 * t111;
t33 = -t146 * t277 + t201 * t149 + t99;
t42 = -t143 * t277 + t64;
t184 = -t148 * t42 + t151 * t33;
t183 = t148 * t33 + t151 * t42;
t94 = t147 * t152 - t149 * t252;
t48 = -t143 * t94 + t146 * t253;
t49 = t143 * t253 + t146 * t94;
t182 = -t148 * t49 + t151 * t48;
t181 = t148 * t48 + t151 * t49;
t179 = (-qJD(2) * pkin(2) + t187) * t150 + t110 * t153;
t177 = qJDD(2) * t150 + t247;
t93 = t152 * t252 + t251;
t115 = t276 * t143;
t174 = pkin(9) * t215 - qJD(5) * t146 + qJD(6) * t115 + t24;
t116 = t276 * t146;
t173 = qJD(5) * t143 + qJD(6) * t116 + (pkin(5) * t152 + t221) * qJD(2) + t23;
t9 = -t35 * qJD(6) + t148 * t66 - t151 * t67;
t172 = g(1) * (t149 * t256 - t89 * t152) - g(2) * (t149 * t198 + t87 * t152) + g(3) * t93;
t45 = t149 * t89 + t152 * t256;
t47 = -t87 * t149 + t152 * t198;
t171 = g(1) * t45 - g(2) * t47 + g(3) * t94;
t170 = t132 * t286;
t40 = -qJD(4) * pkin(4) + qJD(5) - t287;
t165 = -t15 + t172;
t163 = g(3) * t253 - t285;
t161 = -qJ(5) * t233 + (-qJD(5) + t40) * t149;
t160 = t162 + t203;
t159 = t172 + t220;
t158 = qJDD(3) - t160;
t60 = t204 + t227 + (qJD(3) + t217) * qJD(2);
t157 = qJD(2) * t187 - t154 * t155 - t163 + t227 + t60;
t139 = qJDD(4) * t152;
t134 = -pkin(5) * t146 - pkin(4);
t100 = t200 * t152;
t86 = t177 * t145;
t85 = (-qJDD(2) * t153 + t150 * t156) * t145;
t84 = t200 * t234;
t82 = t89 * pkin(2);
t81 = t87 * pkin(2);
t78 = t286 * t152;
t77 = t107 * t152;
t68 = t175 - t264;
t63 = -t143 * t249 + t99;
t51 = qJD(4) * t94 - t152 * t213;
t50 = -qJD(4) * t93 + t149 * t213;
t32 = -pkin(5) * t215 + t53;
t30 = -t148 * t149 * t231 - t151 * t210 + t152 * t284;
t29 = -qJD(4) * t169 - t152 * t166;
t28 = t143 * t212 + t146 * t50;
t27 = -t143 * t50 + t146 * t212;
t26 = pkin(5) * t101 + t40;
t8 = -pkin(5) * t66 + t15;
t1 = [t246, 0, -t85, -t86, t85, t86, qJDD(1) * t147 ^ 2 - g(3) + (qJD(2) * t179 + t150 * t60 - t153 * t68) * t145, 0, 0, 0, 0, 0, -qJD(4) * t51 - qJDD(4) * t93 + (t149 * t177 + t150 * t205) * t145, -qJD(4) * t50 - qJDD(4) * t94 + (t152 * t247 + (-t206 + t222) * t150) * t145, t48 * t223 + t101 * t51 - t66 * t93 + (t149 * t27 + t48 * t233) * qJD(2), -t49 * t223 + t103 * t51 - t67 * t93 + (-t149 * t28 - t49 * t233) * qJD(2), -t101 * t28 - t103 * t27 + t48 * t67 + t49 * t66, t15 * t93 + t16 * t27 + t17 * t28 + t40 * t51 + t48 * t6 + t49 * t7 - g(3), 0, 0, 0, 0, 0 (-qJD(6) * t181 - t148 * t28 + t151 * t27) * t132 + t182 * t105 + t51 * t35 + t93 * t10 -(qJD(6) * t182 + t148 * t27 + t151 * t28) * t132 - t181 * t105 - t51 * t180 + t93 * t9; 0, qJDD(2), t160, t281, t158 - 0.2e1 * t264, 0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t227 - t281, t60 * qJ(3) + t110 * qJD(3) - t68 * pkin(2) - g(1) * (qJ(3) * t90 - t82) - g(2) * (qJ(3) * t88 - t81) - t208 - t179 * t242, qJDD(2) * t142 - 0.2e1 * t149 * t205, -0.2e1 * t149 * t222 + 0.2e1 * t244 * t230, -t149 * t155 + t139, -t152 * t155 - t225, 0, t157 * t149 - t152 * t280, t149 * t280 + t157 * t152, t162 * t143 + (t101 * t218 + t15 * t143 + t154 * t66 + (qJD(2) * t63 + t16) * qJD(4)) * t152 + (t63 * qJDD(2) + t6 + (t101 * t154 - t143 * t40) * qJD(4) + t273 * qJD(2) - t163 * t146) * t149, t162 * t146 + (t103 * t218 + t15 * t146 + t154 * t67 + (-qJD(2) * t64 - t17) * qJD(4)) * t152 + (-t64 * qJDD(2) - t7 + (t103 * t154 - t146 * t40) * qJD(4) - t272 * qJD(2) + t163 * t143) * t149, t63 * t67 + t64 * t66 - t273 * t103 - t272 * t101 + t186 * t234 + (-t143 * t7 - t146 * t6 + t163) * t152, t7 * t64 + t6 * t63 - g(1) * (-pkin(8) * t89 - t82) - g(2) * (-pkin(8) * t87 - t81) - t208 + t272 * t17 + t273 * t16 + (t40 * t234 - t268) * t154 + (-g(3) * pkin(8) * t153 + (-g(3) * t188 + t40 * t241) * t150) * t145 + t285 * t111, -t180 * t29 + t78 * t9, -t10 * t78 + t180 * t30 - t29 * t35 - t77 * t9, t105 * t78 + t132 * t29 + t149 * t9 - t180 * t233, -t10 * t149 - t105 * t77 - t132 * t30 - t233 * t35, t105 * t149 + t132 * t233, t184 * t105 + t219 * t149 + t3 * t233 - t84 * t35 + t100 * t10 + t8 * t77 + t26 * t30 - g(1) * (-t137 * t89 + t258 * t90) - g(2) * (-t137 * t87 + t258 * t88) + (t148 * t274 + t151 * t275) * t132 + (-t132 * t183 - t149 * t4) * qJD(6) + (t35 * t216 - g(3) * (t137 * t153 + t138 * t250)) * t145, -t183 * t105 - t190 * t149 - t4 * t233 + t84 * t180 + t100 * t9 + t8 * t78 + t26 * t29 - g(1) * (-t138 * t89 - t259 * t90) - g(2) * (-t138 * t87 - t259 * t88) + (-t148 * t275 + t151 * t274) * t132 + (-t132 * t184 - t149 * t3) * qJD(6) + (-t180 * t216 - g(3) * (-t137 * t250 + t138 * t153)) * t145; 0, 0, 0, 0, qJDD(2), -t156, t120 + t158 - t239 - t264, 0, 0, 0, 0, 0, t243 * t149 + t139, t243 * t152 - t225, -t143 * t226 + t152 * t66 + (-t146 * t156 + (t101 - 0.2e1 * t214) * qJD(4)) * t149, -t146 * t226 + t152 * t67 + (t143 * t156 + (t103 - 0.2e1 * t211) * qJD(4)) * t149 (qJD(2) * t103 - t101 * t233 + t149 * t66) * t146 + (qJD(2) * t101 + t103 * t233 - t149 * t67) * t143, -t268 + t191 * t149 - t186 * qJD(2) + (t149 * t40 + t152 * t185) * qJD(4) - t162, 0, 0, 0, 0, 0, -qJD(2) * t170 + (-qJD(4) * t107 * t132 - t10) * t152 + (qJD(4) * t35 - t132 * t284 - t260) * t149, t132 * t167 + (-qJD(4) * t170 - t9) * t152 + (-qJD(4) * t180 + t132 * t166 - t261) * t149; 0, 0, 0, 0, 0, 0, 0, t152 * t156 * t149, -t244 * t156, t222, -t223, qJDD(4), qJD(4) * t53 - t152 * t195 + t159, t195 * t149 + t171 - t202, -t143 * t207 + pkin(4) * t66 - t101 * t53 + t165 * t146 + (t143 * t161 - t149 * t23 - t152 * t16) * qJD(2), -t146 * t207 + pkin(4) * t67 - t103 * t53 - t165 * t143 + (t146 * t161 + t149 * t24 + t152 * t17) * qJD(2), t101 * t24 + t103 * t23 + (qJ(5) * t66 - qJD(5) * t101 - t16 * t237 + t7) * t146 + (-qJ(5) * t67 + qJD(5) * t103 - t17 * t237 - t6) * t143 - t171, -t16 * t23 - t17 * t24 - t40 * t53 + t185 * qJD(5) + t165 * pkin(4) + (-t171 + t191) * qJ(5), t107 * t9 - t180 * t271, -t10 * t107 + t180 * t270 - t271 * t35 + t286 * t9, t132 * t271 + t180 * t236 + t260, -t132 * t270 + t236 * t35 + t261, -t132 * t236 (-t115 * t151 - t116 * t148) * t105 + t134 * t10 - t8 * t286 - t3 * t236 - t32 * t35 + t270 * t26 + (t148 * t174 - t151 * t173) * t132 + t172 * t138 -(-t115 * t148 + t116 * t151) * t105 + t134 * t9 + t8 * t107 + t4 * t236 + t32 * t180 + t271 * t26 + (t148 * t173 + t151 * t174) * t132 - t172 * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t103 - t235) * t237 + t194 (-t101 - t231) * t237 + t245, -t101 ^ 2 - t103 ^ 2, t101 * t17 + t103 * t16 - t159 + t282, 0, 0, 0, 0, 0, t10 - t289, t9 - t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180 * t35, t180 ^ 2 - t35 ^ 2, t9 + t269, -t10 - t289, t105, t26 * t180 - g(1) * (-t137 * t45 + t138 * t90) - g(2) * (t137 * t47 + t138 * t88) - g(3) * (-t137 * t94 + t138 * t253) + t219 + t283 * t4, t26 * t35 - g(1) * (-t137 * t90 - t138 * t45) - g(2) * (-t137 * t88 + t138 * t47) - g(3) * (-t137 * t253 - t138 * t94) - t190 + t283 * t3;];
tau_reg  = t1;
