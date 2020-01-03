% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPP7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:57
% EndTime: 2019-12-31 21:06:07
% DurationCPUTime: 3.76s
% Computational Cost: add. (2764->449), mult. (6032->528), div. (0->0), fcn. (3782->6), ass. (0->205)
t150 = cos(qJ(2));
t222 = qJD(1) * t150;
t289 = qJD(3) - t222;
t149 = cos(qJ(3));
t147 = sin(qJ(2));
t223 = qJD(1) * t147;
t203 = t149 * t223;
t146 = sin(qJ(3));
t221 = qJD(2) * t146;
t83 = t203 + t221;
t244 = t83 * t289;
t211 = t147 * qJDD(1);
t288 = qJD(3) * t223 - qJDD(2);
t29 = ((qJD(3) + t222) * qJD(2) + t211) * t146 + t288 * t149;
t291 = -t29 - t244;
t134 = t150 * qJDD(1);
t212 = qJD(1) * qJD(2);
t281 = -t147 * t212 + t134;
t79 = qJDD(3) - t281;
t290 = t79 * qJ(4) + qJD(4) * t289;
t151 = cos(qJ(1));
t234 = t147 * t151;
t148 = sin(qJ(1));
t267 = g(2) * t148;
t283 = g(1) * t234 + t147 * t267;
t238 = t146 * qJ(4);
t273 = pkin(3) + pkin(4);
t172 = -t273 * t149 - t238;
t76 = pkin(2) - t172;
t214 = t149 * qJD(2);
t81 = t146 * t223 - t214;
t287 = t29 * qJ(5) + t81 * qJD(5);
t286 = 0.2e1 * t290;
t78 = t83 ^ 2;
t285 = -t289 ^ 2 - t78;
t249 = t289 * t81;
t198 = t150 * t212;
t28 = -qJD(3) * t214 + (-t198 - t211) * t149 + t288 * t146;
t284 = -t28 - t249;
t136 = t147 * pkin(7);
t139 = t150 * pkin(2);
t206 = -pkin(1) - t139;
t179 = t206 - t136;
t70 = t179 * qJD(1);
t131 = pkin(6) * t222;
t99 = qJD(2) * pkin(7) + t131;
t34 = -t146 * t99 + t149 * t70;
t227 = qJD(4) - t34;
t282 = t146 * qJD(4) + t131;
t280 = g(1) * t151 + t267;
t71 = t79 * pkin(3);
t279 = t71 - qJDD(4);
t98 = -qJD(2) * pkin(2) + pkin(6) * t223;
t174 = t83 * qJ(4) - t98;
t26 = pkin(3) * t81 - t174;
t271 = pkin(7) * t79;
t278 = -t26 * t289 + t271;
t14 = -t273 * t81 + qJD(5) + t174;
t217 = qJD(3) * t149;
t218 = qJD(3) * t146;
t184 = pkin(2) * t147 - pkin(7) * t150;
t88 = t184 * qJD(2);
t40 = qJD(1) * t88 + t179 * qJDD(1);
t60 = t281 * pkin(6) + qJDD(2) * pkin(7);
t194 = t146 * t60 - t149 * t40 + t99 * t217 + t70 * t218;
t237 = t146 * t147;
t232 = t148 * t150;
t63 = t146 * t232 + t149 * t151;
t229 = t151 * t146;
t65 = -t148 * t149 + t150 * t229;
t162 = g(1) * t65 + g(2) * t63 + g(3) * t237 - t194;
t161 = t162 + t279;
t253 = qJ(5) * t28;
t277 = (qJD(5) + t14) * t83 + t161 - t253;
t275 = -0.2e1 * pkin(1);
t274 = t81 ^ 2;
t272 = pkin(4) * t79;
t270 = pkin(3) * t146;
t269 = g(1) * t148;
t266 = g(2) * t151;
t265 = g(3) * t150;
t264 = t83 * t81;
t263 = pkin(7) - qJ(5);
t210 = t273 * t146;
t240 = qJ(4) * t149;
t173 = -t210 + t240;
t262 = t289 * t173 + t282;
t213 = t149 * qJD(5);
t235 = t147 * t149;
t236 = t146 * t150;
t87 = t184 * qJD(1);
t67 = t146 * t87;
t241 = qJ(4) * t223 + t67;
t261 = -t263 * t218 - t213 - (-pkin(6) * t235 + qJ(5) * t236) * qJD(1) - t241;
t181 = -t240 + t270;
t260 = t289 * t181 - t282;
t205 = -pkin(6) * t146 - pkin(3);
t192 = -pkin(4) + t205;
t231 = t149 * t150;
t248 = t149 * t87;
t93 = t263 * t149;
t259 = qJD(3) * t93 - t146 * qJD(5) + t248 - (-qJ(5) * t231 + t192 * t147) * qJD(1);
t225 = -t136 - t139;
t90 = -pkin(1) + t225;
t258 = t146 * t88 + t90 * t217;
t35 = t146 * t70 + t149 * t99;
t257 = t283 * t146;
t256 = t283 * t149;
t255 = pkin(7) * qJD(3);
t254 = qJ(4) * t29;
t102 = t289 * qJ(4);
t18 = qJ(5) * t81 + t35;
t13 = t102 + t18;
t252 = t289 * t13;
t22 = t102 + t35;
t251 = t289 * t22;
t250 = t289 * t35;
t247 = t150 * t83;
t246 = t28 * t146;
t245 = t81 * qJ(4);
t113 = pkin(6) * t231;
t243 = qJD(3) * t113 + t90 * t218;
t242 = t146 * t90 + t113;
t239 = qJD(2) * t83;
t154 = qJD(1) ^ 2;
t233 = t147 * t154;
t230 = t150 * t151;
t17 = t83 * qJ(5) + t34;
t228 = qJD(4) - t17;
t143 = t147 ^ 2;
t224 = -t150 ^ 2 + t143;
t220 = qJD(2) * t147;
t219 = qJD(2) * t150;
t216 = qJD(4) * t149;
t209 = t14 * t218;
t208 = t14 * t217;
t129 = pkin(6) * t211;
t61 = -qJDD(2) * pkin(2) + pkin(6) * t198 + t129;
t207 = g(1) * t230 + g(2) * t232 + g(3) * t147;
t5 = t29 * pkin(3) + t28 * qJ(4) - t83 * qJD(4) + t61;
t3 = -pkin(4) * t29 + qJDD(5) - t5;
t204 = t3 - t265;
t202 = t289 * t221;
t201 = t289 * t214;
t200 = t289 * t218;
t64 = t148 * t231 - t229;
t196 = -t63 * pkin(3) + qJ(4) * t64;
t66 = t146 * t148 + t149 * t230;
t195 = -t65 * pkin(3) + qJ(4) * t66;
t111 = pkin(6) * t236;
t193 = t149 * t90 - t111;
t190 = pkin(3) * t231 + qJ(4) * t236 - t225;
t189 = g(1) * t63 - g(2) * t65;
t188 = g(1) * t64 - g(2) * t66;
t187 = -t64 * pkin(3) + t151 * pkin(6) - t63 * qJ(4);
t186 = t149 * t88 - t243;
t42 = -qJ(4) * t150 + t242;
t185 = t205 * t147;
t183 = qJD(3) * t98 - t271;
t182 = pkin(3) * t149 + t238;
t21 = -pkin(3) * t289 + t227;
t180 = -t146 * t22 + t149 * t21;
t178 = qJ(4) * t220 - t150 * qJD(4) + t258;
t6 = t194 - t279;
t177 = pkin(2) + t182;
t176 = -t255 * t289 - t265;
t171 = -pkin(6) * qJDD(2) + t212 * t275;
t170 = t146 * t79 + t217 * t289;
t169 = t149 * t79 - t200;
t168 = t146 * t40 + t149 * t60 + t70 * t217 - t99 * t218;
t167 = t176 - t5;
t153 = qJD(2) ^ 2;
t165 = pkin(6) * t153 + qJDD(1) * t275 - t269;
t164 = t151 * pkin(1) + pkin(2) * t230 + t66 * pkin(3) + t148 * pkin(6) + pkin(7) * t234 + qJ(4) * t65;
t4 = t168 + t290;
t163 = -t79 + t264;
t160 = t28 - t249;
t159 = t26 * t83 - t161;
t157 = g(1) * t66 + g(2) * t64 + g(3) * t235 - t168;
t156 = t289 * t34 + t157;
t138 = t150 * pkin(3);
t121 = g(2) * t234;
t116 = pkin(7) * t230;
t112 = pkin(7) * t232;
t106 = qJ(4) * t235;
t92 = t263 * t146;
t50 = -t106 + (pkin(6) + t270) * t147;
t43 = t138 - t193;
t41 = t106 + (-pkin(6) - t210) * t147;
t39 = pkin(3) * t83 + t245;
t38 = qJD(1) * t185 - t248;
t37 = -pkin(6) * t203 + t241;
t33 = qJ(5) * t237 + t42;
t30 = t150 * pkin(4) + t111 + t138 + (-qJ(5) * t147 - t90) * t149;
t19 = -t273 * t83 - t245;
t16 = (t182 * qJD(3) - t216) * t147 + (pkin(6) + t181) * t219;
t15 = qJD(2) * t185 - t186;
t11 = (-t147 * t214 - t150 * t218) * pkin(6) + t178;
t10 = (t172 * qJD(3) + t216) * t147 + (-pkin(6) + t173) * t219;
t9 = -t273 * t289 + t228;
t8 = (-pkin(6) * qJD(2) + qJ(5) * qJD(3)) * t235 + (qJD(5) * t147 + (-pkin(6) * qJD(3) + qJ(5) * qJD(2)) * t150) * t146 + t178;
t7 = (-qJ(5) * t219 - t88) * t149 + (qJ(5) * t218 + t192 * qJD(2) - t213) * t147 + t243;
t2 = t4 + t287;
t1 = -qJD(5) * t83 + t253 - t272 + t6;
t12 = [qJDD(1), -t266 + t269, t280, qJDD(1) * t143 + 0.2e1 * t147 * t198, 0.2e1 * t147 * t134 - 0.2e1 * t224 * t212, qJDD(2) * t147 + t150 * t153, qJDD(2) * t150 - t147 * t153, 0, t171 * t147 + (-t165 - t266) * t150, t165 * t147 + t171 * t150 + t121, t214 * t247 + (-t28 * t149 - t83 * t218) * t147, (-t146 * t83 - t149 * t81) * t219 + (t246 - t149 * t29 + (t146 * t81 - t149 * t83) * qJD(3)) * t147, (t28 + t201) * t150 + (t169 + t239) * t147, (t29 - t202) * t150 + (-qJD(2) * t81 - t170) * t147, -t150 * t79 + t220 * t289, t186 * t289 + t193 * t79 + ((pkin(6) * t81 + t146 * t98) * qJD(2) + t194) * t150 + (t98 * t217 + t34 * qJD(2) + t61 * t146 + (t29 + t202) * pkin(6)) * t147 + t188, -t258 * t289 - t242 * t79 + (t98 * t214 + (t200 + t239) * pkin(6) + t168) * t150 + (-t98 * t218 - t35 * qJD(2) + t61 * t149 + (-t28 + t201) * pkin(6)) * t147 - t189, -t15 * t289 + t16 * t81 + t50 * t29 - t43 * t79 + (t26 * t221 + t6) * t150 + (-qJD(2) * t21 + t5 * t146 + t26 * t217) * t147 + t188, -t11 * t81 + t15 * t83 - t43 * t28 - t42 * t29 - t121 + t180 * t219 + (t269 - t146 * t4 + t149 * t6 + (-t146 * t21 - t149 * t22) * qJD(3)) * t147, t11 * t289 - t16 * t83 + t50 * t28 + t42 * t79 + (-t26 * t214 - t4) * t150 + (qJD(2) * t22 - t5 * t149 + t26 * t218) * t147 + t189, -g(1) * t187 - g(2) * t164 + t22 * t11 + t21 * t15 + t26 * t16 - t179 * t269 + t4 * t42 + t6 * t43 + t5 * t50, -t10 * t81 - t7 * t289 - t41 * t29 - t30 * t79 + (-t14 * t221 + t1) * t150 + (-qJD(2) * t9 - t3 * t146 - t208) * t147 + t188, t10 * t83 + t8 * t289 - t41 * t28 + t33 * t79 + (t14 * t214 - t2) * t150 + (qJD(2) * t13 + t149 * t3 - t209) * t147 + t189, t30 * t28 + t33 * t29 - t7 * t83 + t8 * t81 + t121 + (t13 * t146 - t149 * t9) * t219 + (-t269 - t1 * t149 + t146 * t2 + (t13 * t149 + t146 * t9) * qJD(3)) * t147, t2 * t33 + t13 * t8 + t1 * t30 + t9 * t7 + t3 * t41 + t14 * t10 - g(1) * (-t64 * pkin(4) + t187) - g(2) * (pkin(4) * t66 - qJ(5) * t234 + t164) - (-t263 * t147 + t206) * t269; 0, 0, 0, -t150 * t233, t224 * t154, t211, t134, qJDD(2), pkin(1) * t233 - t129 - t265 + t283, (pkin(1) * t154 - pkin(6) * qJDD(1)) * t150 + t207, t149 * t244 - t246, t291 * t146 + t284 * t149, (-t147 * t83 - t231 * t289) * qJD(1) + t170, (t147 * t81 + t236 * t289) * qJD(1) + t169, -t289 * t223, -pkin(2) * t29 + t183 * t146 + (-t265 - t61 - (t87 + t255) * t289) * t149 + (-t98 * t236 - t34 * t147 + (-t150 * t81 - t237 * t289) * pkin(6)) * qJD(1) + t256, pkin(2) * t28 + t67 * t289 + t183 * t149 + (-t176 + t61) * t146 + (-t98 * t231 + t35 * t147 + (-t235 * t289 - t247) * pkin(6)) * qJD(1) - t257, -t278 * t146 + t167 * t149 - t177 * t29 + t21 * t223 + t260 * t81 + t289 * t38 + t256, t37 * t81 - t38 * t83 + (t4 + t289 * t21 + (qJD(3) * t83 - t29) * pkin(7)) * t149 + (t6 - t251 + (qJD(3) * t81 - t28) * pkin(7)) * t146 - t207, t167 * t146 + t278 * t149 - t177 * t28 - t22 * t223 - t260 * t83 - t289 * t37 + t257, -t22 * t37 - t21 * t38 - g(1) * t116 - g(2) * t112 - g(3) * t190 + t260 * t26 + (t180 * qJD(3) + t6 * t146 + t4 * t149) * pkin(7) + (t280 * t147 - t5) * t177, -t209 - t29 * t76 - t79 * t92 - t262 * t81 + t204 * t149 - t259 * t289 + (t14 * t236 + t147 * t9) * qJD(1) + t256, t208 - t76 * t28 + t93 * t79 + t262 * t83 + t204 * t146 + t261 * t289 + (-t13 * t147 - t14 * t231) * qJD(1) + t257, t28 * t92 + t29 * t93 - t259 * t83 + t261 * t81 + (-t289 * t9 - t2) * t149 + (-t1 + t252) * t146 + t207, t2 * t93 + t1 * t92 + t3 * t76 - g(1) * (-qJ(5) * t230 + t116) - g(2) * (-qJ(5) * t232 + t112) - g(3) * (pkin(4) * t231 + t190) + t259 * t9 + t262 * t14 + t261 * t13 + (g(3) * qJ(5) + t280 * t76) * t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, t78 - t274, -t160, -t29 + t244, t79, -t83 * t98 + t162 + t250, t81 * t98 + t156, -t39 * t81 - t159 + t250 + t71, pkin(3) * t28 - t254 + (t22 - t35) * t83 + (t21 - t227) * t81, -t26 * t81 + t39 * t83 - t156 + t286, t4 * qJ(4) - t6 * pkin(3) - t26 * t39 - t21 * t35 - g(1) * t195 - g(2) * t196 - g(3) * (-pkin(3) * t237 + t106) + t227 * t22, t289 * t18 + t19 * t81 + (pkin(4) + t273) * t79 + t277, t14 * t81 - t17 * t289 - t19 * t83 - t157 + t286 + t287, t254 - t273 * t28 + (-t13 + t18) * t83 + (-t9 + t228) * t81, t2 * qJ(4) - t1 * t273 - t9 * t18 - t14 * t19 - g(1) * (-pkin(4) * t65 + t195) - g(2) * (-pkin(4) * t63 + t196) - g(3) * (-t147 * t210 + t106) + t228 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, -t160, t285, t159 - t251, t163, t285, t160, -t252 - t272 - t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291, t284, -t78 - t274, -t13 * t81 + t83 * t9 + t204 + t283;];
tau_reg = t12;
