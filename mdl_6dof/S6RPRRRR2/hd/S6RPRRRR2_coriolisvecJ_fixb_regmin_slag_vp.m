% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:58:30
% EndTime: 2019-03-09 06:58:41
% DurationCPUTime: 3.73s
% Computational Cost: add. (5421->345), mult. (12440->479), div. (0->0), fcn. (9208->10), ass. (0->208)
t180 = qJD(3) + qJD(4);
t187 = sin(qJ(4));
t191 = cos(qJ(3));
t291 = cos(qJ(4));
t235 = qJD(1) * t291;
t188 = sin(qJ(3));
t255 = qJD(1) * t188;
t309 = -t187 * t255 + t191 * t235;
t105 = t309 * t180;
t294 = qJD(5) + qJD(6);
t310 = t309 - t294;
t234 = qJD(4) * t291;
t171 = sin(pkin(11)) * pkin(1) + pkin(7);
t289 = pkin(8) + t171;
t228 = t289 * qJD(1);
t133 = t188 * qJD(2) + t228 * t191;
t124 = t187 * t133;
t132 = t191 * qJD(2) - t228 * t188;
t81 = t291 * t132 - t124;
t305 = pkin(3) * t234 - t81;
t186 = sin(qJ(5));
t190 = cos(qJ(5));
t254 = qJD(1) * t191;
t145 = -t187 * t254 - t188 * t235;
t107 = -pkin(4) * t145 - pkin(9) * t309;
t95 = pkin(3) * t255 + t107;
t308 = -t186 * t305 - t190 * t95;
t127 = -t145 * t186 - t190 * t180;
t185 = sin(qJ(6));
t189 = cos(qJ(6));
t209 = t145 * t190 - t180 * t186;
t210 = t127 * t185 + t189 * t209;
t70 = t189 * t127 - t185 * t209;
t307 = t210 * t70;
t151 = t185 * t186 - t189 * t190;
t277 = t310 * t151;
t262 = t185 * t190;
t153 = t186 * t189 + t262;
t306 = t310 * t153;
t304 = t210 ^ 2 - t70 ^ 2;
t250 = qJD(6) * t185;
t125 = t291 * t133;
t281 = qJD(3) * pkin(3);
t126 = t132 + t281;
t76 = t187 * t126 + t125;
t67 = pkin(9) * t180 + t76;
t172 = -cos(pkin(11)) * pkin(1) - pkin(2);
t158 = -pkin(3) * t191 + t172;
t146 = t158 * qJD(1);
t85 = -pkin(4) * t309 + t145 * pkin(9) + t146;
t41 = t186 * t85 + t190 * t67;
t26 = -pkin(10) * t127 + t41;
t23 = t26 * t250;
t75 = t291 * t126 - t124;
t66 = -t180 * pkin(4) - t75;
t44 = t127 * pkin(5) + t66;
t303 = t44 * t70 + t23;
t142 = qJD(5) - t309;
t137 = qJD(6) + t142;
t249 = qJD(6) * t189;
t251 = qJD(5) * t190;
t252 = qJD(5) * t186;
t62 = t190 * t105 + t145 * t252 + t180 * t251;
t63 = -t209 * qJD(5) + t186 * t105;
t16 = -t127 * t249 - t185 * t63 + t189 * t62 + t209 * t250;
t302 = t137 * t70 + t16;
t154 = t187 * t191 + t291 * t188;
t117 = t180 * t154;
t106 = t117 * qJD(1);
t122 = t132 * qJD(3);
t123 = t133 * qJD(3);
t253 = qJD(4) * t187;
t34 = t291 * t122 - t187 * t123 + t126 * t234 - t133 * t253;
t247 = qJD(1) * qJD(3);
t233 = t188 * t247;
t49 = pkin(3) * t233 + pkin(4) * t106 - pkin(9) * t105;
t47 = t190 * t49;
t197 = -t41 * qJD(5) - t186 * t34 + t47;
t3 = pkin(5) * t106 - pkin(10) * t62 + t197;
t206 = t186 * t49 + t190 * t34 + t85 * t251 - t67 * t252;
t5 = -pkin(10) * t63 + t206;
t241 = -t185 * t5 + t189 * t3;
t40 = -t186 * t67 + t190 * t85;
t25 = pkin(10) * t209 + t40;
t20 = pkin(5) * t142 + t25;
t280 = t189 * t26;
t8 = t185 * t20 + t280;
t301 = -t8 * qJD(6) + t44 * t210 + t241;
t196 = t210 * qJD(6) - t185 * t62 - t189 * t63;
t300 = -t137 * t210 + t196;
t35 = t187 * t122 + t291 * t123 + t126 * t253 + t133 * t234;
t299 = t35 * t190 - t66 * t252;
t80 = t187 * t132 + t125;
t223 = pkin(3) * t253 - t80;
t101 = t153 * t154;
t268 = t309 * t186;
t298 = (t252 - t268) * pkin(5);
t147 = t289 * t188;
t148 = t289 * t191;
t297 = -t291 * t147 - t187 * t148;
t296 = t186 * t95 - t305 * t190;
t295 = qJD(1) * t154;
t237 = t154 * t252;
t205 = -t187 * t188 + t291 * t191;
t116 = t180 * t205;
t273 = t116 * t190;
t203 = t237 - t273;
t259 = t190 * t106;
t293 = -t142 * t203 + t154 * t259;
t292 = -pkin(9) - pkin(10);
t290 = t190 * pkin(5);
t175 = pkin(3) * t187 + pkin(9);
t288 = -pkin(10) - t175;
t287 = -t117 * t210 - t16 * t205;
t264 = t154 * t190;
t265 = t154 * t186;
t274 = t116 * t186;
t22 = t116 * t262 - t185 * t237 - t250 * t265 + (t264 * t294 + t274) * t189;
t286 = -t101 * t106 - t22 * t137;
t285 = -t117 * t209 - t205 * t62;
t284 = t186 * t107 + t190 * t75;
t104 = -t187 * t147 + t291 * t148;
t96 = t190 * t104;
t98 = -pkin(4) * t205 - pkin(9) * t154 + t158;
t282 = t186 * t98 + t96;
t278 = t62 * t186;
t275 = t298 + t223;
t272 = t127 * t142;
t271 = t209 * t142;
t270 = t137 * t145;
t269 = t142 * t145;
t267 = t309 * t190;
t266 = t145 * t309;
t263 = t180 * t116;
t260 = t186 * t106;
t192 = qJD(3) ^ 2;
t258 = t192 * t188;
t257 = t192 * t191;
t256 = t188 ^ 2 - t191 ^ 2;
t161 = qJD(1) * t172;
t246 = pkin(10) * t268;
t245 = t188 * t281;
t240 = qJD(5) * t292;
t236 = t154 * t251;
t232 = qJD(6) * t20 + t5;
t230 = qJD(3) * t289;
t229 = qJD(5) * t288;
t227 = t190 * t107 - t186 * t75;
t226 = t142 * t190;
t176 = -t291 * pkin(3) - pkin(4);
t224 = -t41 * t145 + t35 * t186 + t66 * t251;
t222 = -t76 + t298;
t221 = -t145 * pkin(5) - pkin(10) * t267;
t149 = t288 * t186;
t220 = -qJD(6) * t149 - t186 * t229 - t246 + t296;
t179 = t190 * pkin(10);
t150 = t175 * t190 + t179;
t219 = qJD(6) * t150 - t190 * t229 + t221 - t308;
t162 = t292 * t186;
t218 = -qJD(6) * t162 - t186 * t240 - t246 + t284;
t163 = pkin(9) * t190 + t179;
t217 = qJD(6) * t163 - t190 * t240 + t221 + t227;
t216 = -t117 * t70 - t196 * t205;
t102 = t151 * t154;
t21 = -t101 * t294 - t151 * t116;
t213 = t102 * t106 - t137 * t21;
t212 = -t106 * t175 - t309 * t66;
t211 = -t117 * t127 + t205 * t63;
t208 = 0.2e1 * qJD(3) * t161;
t207 = t40 * t145 - t299;
t204 = t236 + t274;
t139 = t188 * t230;
t140 = t191 * t230;
t51 = qJD(4) * t297 - t291 * t139 - t187 * t140;
t61 = pkin(4) * t117 - pkin(9) * t116 + t245;
t202 = -t104 * t252 + t186 * t61 + t190 * t51 + t98 * t251;
t18 = pkin(5) * t63 + t35;
t7 = -t185 * t26 + t189 * t20;
t201 = t7 * t145 + t18 * t151 - t306 * t44;
t200 = -t8 * t145 + t18 * t153 + t277 * t44;
t199 = t145 * t146 - t35;
t195 = -t204 * t142 - t154 * t260;
t194 = -t146 * t309 - t34;
t52 = t104 * qJD(4) - t187 * t139 + t291 * t140;
t193 = qJD(1) ^ 2;
t177 = -pkin(4) - t290;
t159 = t176 - t290;
t108 = t180 * t117;
t94 = t190 * t98;
t87 = t145 ^ 2 - t309 ^ 2;
t86 = t106 * t205;
t83 = pkin(5) * t265 - t297;
t79 = (-t145 - t295) * t180;
t57 = t190 * t61;
t43 = -pkin(10) * t265 + t282;
t42 = -pkin(5) * t205 - pkin(10) * t264 - t104 * t186 + t94;
t31 = t204 * pkin(5) + t52;
t29 = t142 * t226 - t145 * t209 + t260;
t28 = -t142 ^ 2 * t186 - t127 * t145 + t259;
t24 = -t209 * t226 + t278;
t12 = -t151 * t106 + t137 * t306 - t70 * t145;
t11 = t153 * t106 + t277 * t137 - t145 * t210;
t10 = -t204 * pkin(10) + t202;
t9 = (t62 - t272) * t190 + (-t63 + t271) * t186;
t6 = -pkin(10) * t273 + pkin(5) * t117 - t186 * t51 + t57 + (-t96 + (pkin(10) * t154 - t98) * t186) * qJD(5);
t4 = t16 * t153 - t210 * t277;
t1 = -t151 * t16 + t153 * t196 - t210 * t306 - t277 * t70;
t2 = [0, 0, 0, 0, 0.2e1 * t191 * t233, -0.2e1 * t256 * t247, t257, -t258, 0, -t171 * t257 + t188 * t208, t171 * t258 + t191 * t208, t105 * t154 - t116 * t145, t105 * t205 - t106 * t154 + t116 * t309 + t117 * t145, t263, -t108, 0, t158 * t106 + t146 * t117 - t52 * t180 + (-qJD(1) * t205 - t309) * t245, t158 * t105 + t146 * t116 - t51 * t180 + (-t145 + t295) * t245, t203 * t209 + t62 * t264 (-t127 * t190 + t186 * t209) * t116 + (-t278 - t190 * t63 + (t127 * t186 + t190 * t209) * qJD(5)) * t154, t285 + t293, t195 + t211, t117 * t142 - t86 (-t104 * t251 + t57) * t142 + t94 * t106 - (-t67 * t251 + t47) * t205 + t40 * t117 + t52 * t127 - t297 * t63 + t66 * t236 + ((-qJD(5) * t98 - t51) * t142 - t104 * t106 - (-qJD(5) * t85 - t34) * t205 + t35 * t154 + t66 * t116) * t186, -t282 * t106 - t41 * t117 - t202 * t142 + t154 * t299 + t205 * t206 - t209 * t52 + t66 * t273 - t297 * t62, -t102 * t16 - t21 * t210, -t101 * t16 - t102 * t196 - t21 * t70 + t210 * t22, -t213 + t287, t216 + t286, t117 * t137 - t86 (-t185 * t10 + t189 * t6) * t137 + (-t185 * t43 + t189 * t42) * t106 - t241 * t205 + t7 * t117 + t31 * t70 - t83 * t196 + t18 * t101 + t44 * t22 + ((-t185 * t42 - t189 * t43) * t137 + t8 * t205) * qJD(6), -t18 * t102 - t8 * t117 - t23 * t205 + t83 * t16 + t44 * t21 - t31 * t210 + (-(-qJD(6) * t43 + t6) * t137 - t42 * t106 + t3 * t205) * t185 + (-(qJD(6) * t42 + t10) * t137 - t43 * t106 + t232 * t205) * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t258, -t257, 0, 0, 0, 0, 0, -t108, -t263, 0, 0, 0, 0, 0, t195 - t211, t285 - t293, 0, 0, 0, 0, 0, -t216 + t286, t213 + t287; 0, 0, 0, 0, -t188 * t193 * t191, t256 * t193, 0, 0, 0, -t161 * t255, -t161 * t254, t266, t87, 0, t79, 0, t180 * t80 + (-t180 * t253 + t255 * t309) * pkin(3) + t199, t81 * t180 + (t145 * t255 - t180 * t234) * pkin(3) + t194, t24, t9, t29, t28, t269, t176 * t63 + t212 * t186 + t223 * t127 + (-t175 * t251 + t308) * t142 + t207, t176 * t62 + t212 * t190 - t223 * t209 + (t175 * t252 + t296) * t142 + t224, t4, t1, t11, t12, t270 (t149 * t189 - t150 * t185) * t106 - t159 * t196 + t275 * t70 + (t185 * t220 - t189 * t219) * t137 + t201 -(t149 * t185 + t150 * t189) * t106 + t159 * t16 - t275 * t210 + (t185 * t219 + t189 * t220) * t137 + t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266, t87, 0, t79, 0, t180 * t76 + t199, t75 * t180 + t194, t24, t9, t29, t28, t269, -pkin(4) * t63 - t227 * t142 - t76 * t127 - t66 * t268 + (-t142 * t251 - t260) * pkin(9) + t207, -pkin(4) * t62 + t284 * t142 + t76 * t209 - t66 * t267 + (t142 * t252 - t259) * pkin(9) + t224, t4, t1, t11, t12, t270 (t162 * t189 - t163 * t185) * t106 - t177 * t196 + t222 * t70 + (t185 * t218 - t189 * t217) * t137 + t201 -(t162 * t185 + t163 * t189) * t106 + t177 * t16 - t222 * t210 + (t185 * t217 + t189 * t218) * t137 + t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209 * t127, -t127 ^ 2 + t209 ^ 2, t62 + t272, -t63 - t271, t106, t142 * t41 + t209 * t66 + t197, t127 * t66 + t142 * t40 - t206, -t307, t304, t302, t300, t106 -(-t185 * t25 - t280) * t137 + (t106 * t189 - t137 * t250 + t209 * t70) * pkin(5) + t301 (-t137 * t26 - t3) * t185 + (t137 * t25 - t232) * t189 + (-t106 * t185 - t137 * t249 - t209 * t210) * pkin(5) + t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t307, t304, t302, t300, t106, t137 * t8 + t301, t137 * t7 - t185 * t3 - t189 * t232 + t303;];
tauc_reg  = t2;
