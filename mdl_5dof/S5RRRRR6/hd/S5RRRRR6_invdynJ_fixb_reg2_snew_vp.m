% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRR6
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
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:21
% EndTime: 2020-01-03 12:15:29
% DurationCPUTime: 2.88s
% Computational Cost: add. (23822->306), mult. (30790->429), div. (0->0), fcn. (20952->10), ass. (0->214)
t223 = sin(qJ(5));
t221 = qJD(1) + qJD(2);
t224 = sin(qJ(4));
t228 = cos(qJ(4));
t229 = cos(qJ(3));
t225 = sin(qJ(3));
t270 = t221 * t225;
t187 = -t228 * t229 * t221 + t224 * t270;
t189 = (t229 * t224 + t225 * t228) * t221;
t227 = cos(qJ(5));
t152 = t227 * t187 + t223 * t189;
t154 = -t223 * t187 + t227 * t189;
t120 = t154 * t152;
t218 = qJDD(3) + qJDD(4);
t213 = qJDD(5) + t218;
t287 = -t120 + t213;
t292 = t223 * t287;
t165 = t189 * t187;
t286 = -t165 + t218;
t291 = t224 * t286;
t290 = t227 * t287;
t289 = t228 * t286;
t220 = qJD(3) + qJD(4);
t216 = qJD(5) + t220;
t146 = t216 * t152;
t256 = qJD(3) * t221;
t250 = t229 * t256;
t219 = qJDD(1) + qJDD(2);
t264 = t225 * t219;
t196 = t250 + t264;
t251 = t225 * t256;
t259 = t229 * t219;
t237 = -t251 + t259;
t242 = t224 * t196 - t228 * t237;
t138 = -t189 * qJD(4) - t242;
t139 = -t187 * qJD(4) + t228 * t196 + t224 * t237;
t90 = -t152 * qJD(5) + t223 * t138 + t227 * t139;
t288 = -t146 + t90;
t281 = sin(qJ(1));
t282 = cos(qJ(1));
t235 = t281 * g(2) - t282 * g(3);
t202 = -qJD(1) ^ 2 * pkin(1) - t235;
t226 = sin(qJ(2));
t230 = cos(qJ(2));
t236 = -t282 * g(2) - t281 * g(3);
t232 = qJDD(1) * pkin(1) + t236;
t168 = t230 * t202 + t226 * t232;
t284 = t221 ^ 2;
t163 = -t284 * pkin(2) + t219 * pkin(7) + t168;
t266 = t225 * t163;
t148 = t229 * g(1) + t266;
t280 = t225 * g(1);
t149 = t229 * t163 - t280;
t111 = t225 * t148 + t229 * t149;
t182 = t220 * t187;
t285 = -t182 + t139;
t132 = t182 + t139;
t150 = t152 ^ 2;
t151 = t154 ^ 2;
t185 = t187 ^ 2;
t186 = t189 ^ 2;
t212 = t216 ^ 2;
t217 = t220 ^ 2;
t283 = t229 ^ 2;
t167 = -t226 * t202 + t230 * t232;
t162 = -t219 * pkin(2) - t284 * pkin(7) - t167;
t255 = pkin(8) * t270;
t205 = qJD(3) * pkin(3) - t255;
t211 = t283 * t284;
t126 = -t237 * pkin(3) - pkin(8) * t211 + t205 * t270 + t162;
t178 = t220 * pkin(4) - t189 * pkin(9);
t77 = -t138 * pkin(4) - t185 * pkin(9) + t189 * t178 + t126;
t279 = t223 * t77;
t121 = qJDD(3) * pkin(3) - t196 * pkin(8) - t266 + (-g(1) + (pkin(3) * t270 + pkin(8) * qJD(3)) * t221) * t229;
t258 = t229 * t284;
t122 = -t280 + (-t205 - t255) * qJD(3) + (-pkin(3) * t258 + t219 * pkin(8) + t163) * t229;
t87 = -t228 * t121 + t224 * t122;
t60 = t286 * pkin(4) - t132 * pkin(9) - t87;
t88 = t224 * t121 + t228 * t122;
t61 = -t185 * pkin(4) + t138 * pkin(9) - t220 * t178 + t88;
t36 = t223 * t61 - t227 * t60;
t37 = t223 * t60 + t227 * t61;
t20 = t223 * t37 - t227 * t36;
t278 = t224 * t20;
t52 = t224 * t88 - t228 * t87;
t277 = t225 * t52;
t276 = t227 * t77;
t275 = t228 * t20;
t274 = t216 * t223;
t273 = t216 * t227;
t272 = t220 * t224;
t271 = t220 * t228;
t114 = t120 + t213;
t269 = t223 * t114;
t268 = t224 * t126;
t159 = t165 + t218;
t267 = t224 * t159;
t207 = t225 * t258;
t265 = t225 * (qJDD(3) + t207);
t263 = t227 * t114;
t262 = t228 * t126;
t261 = t228 * t159;
t260 = t229 * (qJDD(3) - t207);
t257 = -pkin(2) * t162 + pkin(7) * t111;
t222 = t225 ^ 2;
t210 = t222 * t284;
t231 = qJD(3) ^ 2;
t176 = -t260 - t225 * (-t210 - t231);
t195 = 0.2e1 * t250 + t264;
t254 = -pkin(2) * t195 + pkin(7) * t176 + t225 * t162;
t175 = t229 * (-t211 - t231) - t265;
t197 = -0.2e1 * t251 + t259;
t253 = pkin(2) * t197 + pkin(7) * t175 - t229 * t162;
t21 = t223 * t36 + t227 * t37;
t17 = -pkin(4) * t77 + pkin(9) * t21;
t7 = t224 * t21 + t275;
t8 = t228 * t21 - t278;
t4 = -t225 * t7 + t229 * t8;
t252 = t225 * (-pkin(8) * t7 - pkin(9) * t275 - t224 * t17) + t229 * (-pkin(3) * t77 + pkin(8) * t8 - pkin(9) * t278 + t228 * t17) - pkin(2) * t77 + pkin(7) * t4;
t106 = -t150 - t151;
t243 = -t227 * t138 + t223 * t139;
t234 = (-qJD(5) + t216) * t154 - t243;
t74 = t146 + t90;
t44 = t223 * t234 - t227 * t74;
t46 = t223 * t74 + t227 * t234;
t24 = t224 * t46 + t228 * t44;
t25 = -t224 * t44 + t228 * t46;
t11 = -t225 * t24 + t229 * t25;
t12 = -pkin(4) * t106 + pkin(9) * t46 + t21;
t14 = -pkin(9) * t44 - t20;
t249 = t225 * (-pkin(8) * t24 - t224 * t12 + t228 * t14) + t229 * (-pkin(3) * t106 + pkin(8) * t25 + t228 * t12 + t224 * t14) - pkin(2) * t106 + pkin(7) * t11;
t53 = t224 * t87 + t228 * t88;
t112 = -t212 - t150;
t79 = t223 * t112 + t290;
t80 = t227 * t112 - t292;
t49 = t224 * t80 + t228 * t79;
t50 = -t224 * t79 + t228 * t80;
t27 = -t225 * t49 + t229 * t50;
t70 = (qJD(5) + t216) * t154 + t243;
t34 = -pkin(4) * t70 + pkin(9) * t80 - t276;
t48 = -pkin(9) * t79 + t279;
t248 = t225 * (-pkin(8) * t49 - t224 * t34 + t228 * t48) + t229 * (-pkin(3) * t70 + pkin(8) * t50 + t224 * t48 + t228 * t34) - pkin(2) * t70 + pkin(7) * t27;
t141 = -t151 - t212;
t97 = t227 * t141 - t269;
t98 = -t223 * t141 - t263;
t57 = t224 * t98 + t228 * t97;
t58 = -t224 * t97 + t228 * t98;
t31 = -t225 * t57 + t229 * t58;
t39 = -pkin(4) * t288 + pkin(9) * t98 + t279;
t51 = -pkin(9) * t97 + t276;
t247 = t225 * (-pkin(8) * t57 - t224 * t39 + t228 * t51) + t229 * (-pkin(3) * t288 + pkin(8) * t58 + t224 * t51 + t228 * t39) - pkin(2) * t288 + pkin(7) * t31;
t140 = -t185 - t186;
t233 = (-qJD(4) + t220) * t189 - t242;
t91 = -t228 * t132 + t224 * t233;
t92 = t224 * t132 + t228 * t233;
t56 = -t225 * t91 + t229 * t92;
t246 = t225 * (-pkin(8) * t91 - t52) + t229 * (-pkin(3) * t140 + pkin(8) * t92 + t53) - pkin(2) * t140 + pkin(7) * t56;
t157 = -t217 - t185;
t116 = t224 * t157 + t289;
t117 = t228 * t157 - t291;
t127 = (qJD(4) + t220) * t189 + t242;
t82 = -t225 * t116 + t229 * t117;
t245 = t225 * (-pkin(8) * t116 + t268) + t229 * (-pkin(3) * t127 + pkin(8) * t117 - t262) - pkin(2) * t127 + pkin(7) * t82;
t177 = -t186 - t217;
t133 = t228 * t177 - t267;
t134 = -t224 * t177 - t261;
t95 = -t225 * t133 + t229 * t134;
t244 = t225 * (-pkin(8) * t133 + t262) + t229 * (-pkin(3) * t285 + pkin(8) * t134 + t268) - pkin(2) * t285 + pkin(7) * t95;
t200 = (t222 + t283) * t219;
t201 = t210 + t211;
t241 = pkin(2) * t201 + pkin(7) * t200 + t111;
t240 = pkin(4) * t79 - t36;
t239 = pkin(4) * t97 - t37;
t29 = t229 * t53 - t277;
t238 = pkin(7) * t29 - pkin(8) * t277 - pkin(2) * t126 + t229 * (-pkin(3) * t126 + pkin(8) * t53);
t180 = -t186 + t217;
t179 = t185 - t217;
t174 = t265 + t229 * (-t210 + t231);
t173 = t225 * (t211 - t231) + t260;
t170 = (t196 + t250) * t225;
t169 = t197 * t229;
t166 = t229 * t195 + t225 * t197;
t164 = t186 - t185;
t145 = -t151 + t212;
t144 = t150 - t212;
t118 = t151 - t150;
t109 = (-t152 * t227 + t154 * t223) * t216;
t108 = (-t152 * t223 - t154 * t227) * t216;
t107 = (t225 * (-t187 * t228 + t189 * t224) + t229 * (-t187 * t224 - t189 * t228)) * t220;
t104 = t225 * (t228 * t179 - t267) + t229 * (t224 * t179 + t261);
t103 = t225 * (-t224 * t180 + t289) + t229 * (t228 * t180 + t291);
t102 = t227 * t144 - t269;
t101 = -t223 * t145 + t290;
t100 = t223 * t144 + t263;
t99 = t227 * t145 + t292;
t89 = -t154 * qJD(5) - t243;
t84 = t225 * (t228 * t139 - t189 * t272) + t229 * (t224 * t139 + t189 * t271);
t83 = t225 * (-t224 * t138 + t187 * t271) + t229 * (t228 * t138 + t187 * t272);
t67 = -t154 * t274 + t227 * t90;
t66 = t154 * t273 + t223 * t90;
t65 = t152 * t273 - t223 * t89;
t64 = t152 * t274 + t227 * t89;
t55 = t225 * (-t228 * t127 - t224 * t285) + t229 * (-t224 * t127 + t228 * t285);
t47 = t225 * (-t224 * t108 + t228 * t109) + t229 * (t228 * t108 + t224 * t109);
t45 = -t223 * t288 - t227 * t70;
t43 = -t223 * t70 + t227 * t288;
t42 = pkin(4) * t44;
t33 = t225 * (-t224 * t100 + t228 * t102) + t229 * (t228 * t100 + t224 * t102);
t32 = t225 * (t228 * t101 - t224 * t99) + t229 * (t224 * t101 + t228 * t99);
t23 = t225 * (-t224 * t66 + t228 * t67) + t229 * (t224 * t67 + t228 * t66);
t22 = t225 * (-t224 * t64 + t228 * t65) + t229 * (t224 * t65 + t228 * t64);
t19 = pkin(4) * t20;
t10 = t225 * (-t224 * t43 + t228 * t45) + t229 * (t224 * t45 + t228 * t43);
t1 = [0, 0, 0, 0, 0, qJDD(1), t236, t235, 0, 0, 0, 0, 0, 0, 0, t219, pkin(1) * (t230 * t219 - t226 * t284) + t167, pkin(1) * (-t226 * t219 - t230 * t284) - t168, 0, pkin(1) * (t230 * t167 + t226 * t168), t170, t166, t174, t169, t173, 0, pkin(1) * (t226 * t175 + t230 * t197) + t253, pkin(1) * (t226 * t176 - t230 * t195) + t254, pkin(1) * (t226 * t200 + t230 * t201) + t241, pkin(1) * (t226 * t111 - t230 * t162) + t257, t84, t55, t103, t83, t104, t107, pkin(1) * (-t230 * t127 + t226 * t82) + t245, pkin(1) * (t226 * t95 - t230 * t285) + t244, pkin(1) * (-t230 * t140 + t226 * t56) + t246, pkin(1) * (-t230 * t126 + t226 * t29) + t238, t23, t10, t32, t22, t33, t47, pkin(1) * (t226 * t27 - t230 * t70) + t248, pkin(1) * (t226 * t31 - t230 * t288) + t247, pkin(1) * (-t230 * t106 + t226 * t11) + t249, pkin(1) * (t226 * t4 - t230 * t77) + t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219, t167, -t168, 0, 0, t170, t166, t174, t169, t173, 0, t253, t254, t241, t257, t84, t55, t103, t83, t104, t107, t245, t244, t246, t238, t23, t10, t32, t22, t33, t47, t248, t247, t249, t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t207, t210 - t211, t264, t207, t259, qJDD(3), -t148, -t149, 0, 0, t165, t164, t132, -t165, t233, t218, pkin(3) * t116 - t87, pkin(3) * t133 - t88, pkin(3) * t91, pkin(3) * t52, t120, t118, t74, -t120, t234, t213, pkin(3) * t49 + t240, pkin(3) * t57 + t239, pkin(3) * t24 + t42, pkin(3) * t7 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, t164, t132, -t165, t233, t218, -t87, -t88, 0, 0, t120, t118, t74, -t120, t234, t213, t240, t239, t42, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t118, t74, -t120, t234, t213, -t36, -t37, 0, 0;];
tauJ_reg = t1;
